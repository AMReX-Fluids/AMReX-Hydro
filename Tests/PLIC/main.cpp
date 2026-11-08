// Verification of the PLIC library, as a small solver would use it. The run fails (non-zero exit)
// when a check does not hold.
//
//   test.kind = kernels   the plane geometry, the normals, the plane fit, the wall ghost cells and
//                         the face density
//   test.kind = advect    a prescribed velocity on an AMR hierarchy (see inputs.*), advanced with
//                         PLIC::ComputeEdgeStateAndFluxes in the solver's own MFIter loop, or with
//                         HydroUtils::ComputeFluxesOnBoxFromState and advection_type "PLIC" if
//                         test.use_hydroutils = 1

#include <PLICTestUtils.H>

#include <hydro_plic.H>
#include <hydro_plic_normals_K.H>
#include <hydro_utils.H>

#include <AMReX.H>
#include <AMReX_AmrCore.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_TagBox.H>

#include <random>
#include <sstream>
#include <string>

using namespace amrex;
using namespace amrex::literals;

namespace {

int n_failures = 0;

void check (bool ok, std::string const& what)
{
    Print() << "  " << (ok ? "ok    " : "FAILED") << " " << what << "\n";
    if (!ok) { ++n_failures; }
}

#ifdef AMREX_USE_FLOAT
constexpr Real roundoff_tol = 2.0e-3_rt;
#else
constexpr Real roundoff_tol = 1.0e-10_rt;
#endif

enum class Prob { Zalesak, Vortex, Translate };

struct ProbParm
{
    Prob type = Prob::Vortex;
    Real radius = 0.15_rt;
    GpuArray<Real, AMREX_SPACEDIM> center = {AMREX_D_DECL(0.5_rt, 0.75_rt, 0.5_rt)};
    GpuArray<Real, AMREX_SPACEDIM> vel = {AMREX_D_DECL(1.0_rt, 1.0_rt, 1.0_rt)};
    Real slot_width = 0.05_rt;
    Real slot_length = 0.25_rt;
    Real period = 4.0_rt;
};

// The vector potential of the prescribed velocity. Back to the initial shape at t = period.
struct Potential
{
    ProbParm p;

    AMREX_GPU_HOST_DEVICE
    void operator() (const Real* x, Real time, Real* A) const noexcept
    {
        constexpr Real pi = Math::pi<Real>();
        A[0] = A[1] = A[2] = 0.0_rt;
        if (p.type == Prob::Zalesak) {
            const Real omega = 2.0_rt * pi / p.period;
            A[2] = -0.5_rt * omega * ((x[0] - 0.5_rt) * (x[0] - 0.5_rt) + (x[1] - 0.5_rt) * (x[1] - 0.5_rt));
        } else if (p.type == Prob::Vortex) {
            const Real ct = std::cos(pi * time / p.period);
            const Real sx2 = std::sin(pi * x[0]) * std::sin(pi * x[0]);
            const Real sy2 = std::sin(pi * x[1]) * std::sin(pi * x[1]);
#if AMREX_SPACEDIM == 2
            A[2] = sx2 * sy2 / pi * ct;
#else
            const Real sz2 = std::sin(pi * x[2]) * std::sin(pi * x[2]);
            A[1] = -sx2 * sz2 / pi * std::sin(2.0_rt * pi * x[1]) * ct;
            A[2] = sx2 * sy2 / pi * std::sin(2.0_rt * pi * x[2]) * ct;
#endif
        } else {
            auto const& U = p.vel;
#if AMREX_SPACEDIM == 2
            A[2] = U[0] * x[1] - U[1] * x[0];
#else
            A[0] = 0.5_rt * (U[1] * x[2] - U[2] * x[1]);
            A[1] = 0.5_rt * (U[2] * x[0] - U[0] * x[2]);
            A[2] = 0.5_rt * (U[0] * x[1] - U[1] * x[0]);
#endif
        }
    }
};

// phi > 0 in the liquid. The Zalesak disk is extruded along z.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
initial_fraction (int i, int j, int k, GeometryData const& geomdata, ProbParm const& p) noexcept
{
    auto const c = p.center;
    const Real R = p.radius;
    const bool slotted = (p.type == Prob::Zalesak);
    const Real hw = 0.5_rt * p.slot_width;
    const Real slot_top = c[1] - R + p.slot_length;
    auto phi = [=] (Real x, Real y, Real z) -> Real {
        Real r2 = (x - c[0]) * (x - c[0]) + (y - c[1]) * (y - c[1]);
#if AMREX_SPACEDIM == 3
        if (!slotted) { r2 += (z - c[2]) * (z - c[2]); }
#else
        amrex::ignore_unused(z);
#endif
        if (r2 > R * R) { return -1.0_rt; }
        if (slotted && std::abs(x - c[0]) < hw && y < slot_top) { return -1.0_rt; }
        return 1.0_rt;
    };
    return Clamp(PLIC::cell_fraction(i, j, k, geomdata, phi, 16), 0.0_rt, 1.0_rt);
}

// The ghost cells beyond the domain are all extrapolated or periodic: nothing is prescribed
auto make_bndry (Geometry const& geom, Vector<BCRec> const& bcrec)
{
    return PhysBCFunct<GpuBndryFuncFab<FabFillNoOp>>(geom, bcrec, GpuBndryFuncFab<FabFillNoOp>{FabFillNoOp{}});
}

// The part of a solver that the PLIC library leaves to it: the mesh and the ghost cells
class AdvectTest final : public AmrCore
{
public:
    AdvectTest ()
    {
        ParmParse ppp("plic");
        ppp.query("reconstruction", m_reconstruction);

        ParmParse pp("prob");
        std::string type = "vortex";
        pp.query("type", type);
        if (type == "zalesak") {
            m_prob.type = Prob::Zalesak;
        } else if (type == "vortex") {
            m_prob.type = Prob::Vortex;
        } else if (type == "translate") {
            m_prob.type = Prob::Translate;
        } else {
            Abort("prob.type must be zalesak, vortex or translate");
        }
        pp.query("radius", m_prob.radius);
        pp.query("period", m_prob.period);
        pp.query("slot_width", m_prob.slot_width);
        pp.query("slot_length", m_prob.slot_length);
        Vector<Real> v;
        if (pp.queryarr("center", v)) { for (int d = 0; d < AMREX_SPACEDIM; ++d) { m_prob.center[d] = v[d]; } }
        v.clear();
        if (pp.queryarr("vel", v)) { for (int d = 0; d < AMREX_SPACEDIM; ++d) { m_prob.vel[d] = v[d]; } }

        ParmParse ppa("amr");
        ppa.query("regrid_int", m_regrid_int);
        ppa.query("fixed_dt", m_fixed_dt);
        ppa.query("cfl", m_cfl);
        ppa.query("stop_time", m_stop_time);
        ppa.query("max_step", m_max_step);

        ParmParse ppt("test");
        ppt.query("refine_band", m_refine_band);
        ppt.query("plot_file", m_plot_file);
        ppt.query("max_l1", m_max_l1);
        ppt.query("max_volume_change", m_max_volume_change);
        ppt.query("max_f_excursion", m_max_f_excursion);
        ppt.query("use_hydroutils", m_use_hydroutils);
        if (m_use_hydroutils) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_reconstruction == "myc_height",
                                             "test.use_hydroutils needs plic.reconstruction = myc_height");
        }

        // All the non-periodic domain faces are walls: extrapolated, which FillWallGhosts then extends
        m_bcrec.resize(1);
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            const bool periodic = Geom(0).isPeriodic(d);
            m_bcrec[0].setLo(d, periodic ? BCType::int_dir : BCType::foextrap);
            m_bcrec[0].setHi(d, periodic ? BCType::int_dir : BCType::foextrap);
        }

        m_bcrec_d.resize(1);
        m_iconserv_d.resize(1);
        Gpu::copyAsync(Gpu::hostToDevice, m_bcrec.begin(), m_bcrec.end(), m_bcrec_d.begin());
        const int iconserv = 1;
        Gpu::copyAsync(Gpu::hostToDevice, &iconserv, &iconserv + 1, m_iconserv_d.begin());
        Gpu::streamSynchronize();

        m_F.resize(max_level + 1);
    }

    void Run ()
    {
        InitFromScratch(0.0_rt);
        for (int lev = finest_level; lev > 0; --lev) {
            average_down(m_F[lev], m_F[lev - 1], 0, 1, refRatio(lev - 1));
        }
        for (int lev = 0; lev <= finest_level; ++lev) {
            FillLevel(lev, m_F[lev]);
        }

        Print() << "  reconstruction " << m_reconstruction
                << (m_use_hydroutils ? ", through HydroUtils" : "") << "\n";

        const Real vol0 = Volume();
        Real max_change = 0.0_rt;
        Real time = 0.0_rt;
        Real prev_dt = 0.0_rt;
        int nstep = 0;
        while (nstep < m_max_step && time < m_stop_time * (1.0_rt - 1.0e-6_rt)) {
            if (m_regrid_int > 0 && max_level > 0 && nstep > 0 && nstep % m_regrid_int == 0) {
                regrid(0, time);
            }
            const Real dt = ComputeDt(time, prev_dt, nstep);

            const Real vol_old = Volume();
            Advance(time, dt);
            max_change = max(max_change, std::abs(Volume() - vol_old) / max(vol_old, std::numeric_limits<Real>::min()));

            time += dt;
            prev_dt = dt;
            ++nstep;
        }

        Real f_min = std::numeric_limits<Real>::max();
        Real f_max = std::numeric_limits<Real>::lowest();
        for (int lev = 0; lev <= finest_level; ++lev) {
            f_min = min(f_min, m_F[lev].min(0));
            f_max = max(f_max, m_F[lev].max(0));
        }
        const Real vol1 = Volume();
        const Real f_out = max(0.0_rt - f_min, f_max - 1.0_rt, 0.0_rt);
        const Real l1 = ShapeError();
        Print() << std::setprecision(6)
                << "  steps " << nstep << ", time " << time << ", finest level " << finest_level << "\n"
                << "  L1 shape error " << l1 << "\n"
                << "  volume " << vol1 << ", total relative change " << (vol1 - vol0) / vol0
                << ", largest per step " << max_change << "\n"
                << "  F left [0,1] by " << f_out << "\n";

        if (m_max_l1 > 0.0_rt) { check(l1 <= m_max_l1, "L1 shape error within " + std::to_string(m_max_l1)); }
        if (m_max_volume_change > 0.0_rt) {
            std::ostringstream os; os << "volume change per step within " << m_max_volume_change;
            check(max_change <= m_max_volume_change, os.str());
        }
        if (m_max_f_excursion > 0.0_rt) {
            std::ostringstream os; os << "F within " << m_max_f_excursion << " of [0,1]";
            check(f_out <= m_max_f_excursion, os.str());
        }

        if (!m_plot_file.empty()) {
            Vector<MultiFab const*> mfs(finest_level + 1);
            Vector<int> steps(finest_level + 1, nstep);
            for (int lev = 0; lev <= finest_level; ++lev) { mfs[lev] = &m_F[lev]; }
            WriteMultiLevelPlotfile(m_plot_file, finest_level + 1, mfs, {"vof"}, Geom(), time, steps, refRatio());
        }
    }

private:
    [[nodiscard]] Real Volume () const
    {
        Vector<MultiFab const*> F(finest_level + 1);
        for (int lev = 0; lev <= finest_level; ++lev) { F[lev] = &m_F[lev]; }
        return testutils::LiquidVolume(F, Geom(), refRatio());
    }

    [[nodiscard]] Real ShapeError () const
    {
        // |F - F(0)| on the same cells, not counting covered cells
        Vector<MultiFab> diff(finest_level + 1);
        Vector<MultiFab const*> d(finest_level + 1);
        for (int lev = 0; lev <= finest_level; ++lev) {
            diff[lev].define(grids[lev], dmap[lev], 1, 0);
            auto const& dm = diff[lev].arrays();
            auto const& Fm = m_F[lev].const_arrays();
            const auto geomdata = Geom(lev).data();
            const auto prob = m_prob;
            ParallelFor(diff[lev], [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept {
                dm[b](i, j, k) = std::abs(Fm[b](i, j, k) - initial_fraction(i, j, k, geomdata, prob));
            });
            Gpu::streamSynchronize();
            d[lev] = &diff[lev];
        }
        return testutils::LiquidVolume(d, Geom(), refRatio());
    }

    [[nodiscard]] Real ComputeDt (Real time, Real prev_dt, int nstep) const
    {
        Real dt = m_fixed_dt;
        if (dt <= 0.0_rt) {
            // The cell-centered velocity at this time
            Vector<MultiFab> vel(finest_level + 1);
            Vector<MultiFab const*> vptr(finest_level + 1);
            for (int lev = 0; lev <= finest_level; ++lev) {
                Array<MultiFab, AMREX_SPACEDIM> um;
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    um[d].define(convert(grids[lev], IntVect::TheDimensionVector(d)), dmap[lev], 1, 0);
                }
                testutils::FaceVelocityFromPotential(GetArrOfPtrs(um), Geom(lev), time, Potential{m_prob});
                vel[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 0);
                average_face_to_cellcenter(vel[lev], 0, GetArrOfConstPtrs(um));
                vptr[lev] = &vel[lev];
            }
            dt = testutils::EstimateDt(vptr, 0, Geom(), m_cfl);
        }
        if (nstep > 0) {
            dt = min(dt, prev_dt * 1.1_rt);
            // Ensure ~O(dt) last step by checking a little in advance
            const Real left = m_stop_time - time;
            if (2.0_rt * dt > left && left > dt) {
                dt = 0.5_rt * left;
            } else {
                dt = min(dt, left);
            }
        }
        return dt;
    }

    // One step on all levels with the same time step, as a solver without subcycling would: the fluxes
    // of each box, the fine fluxes onto the coarse faces under them, the conservative update, and the
    // ghost cells. A subcycling solver would advance each level on its own and reflux with the same
    // fluxes instead.
    void Advance (Real time, Real dt)
    {
        const int nlevs = finest_level + 1;
        Vector<Array<MultiFab, AMREX_SPACEDIM>> umac(nlevs);
        Vector<Array<MultiFab, AMREX_SPACEDIM>> flux(nlevs);
        Vector<Array<MultiFab, AMREX_SPACEDIM>> face(nlevs);
        for (int lev = 0; lev < nlevs; ++lev) {
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                const BoxArray fba = convert(grids[lev], IntVect::TheDimensionVector(d));
                umac[lev][d].define(fba, dmap[lev], 1, 1);
                flux[lev][d].define(fba, dmap[lev], 1, 0);
                face[lev][d].define(fba, dmap[lev], 1, 0);
            }
            // umac at t^{n+1/2}, with a ghost layer
            testutils::FaceVelocityFromPotential(GetArrOfPtrs(umac[lev]), Geom(lev), time + 0.5_rt * dt,
                                                 Potential{m_prob});
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(m_F[lev]); mfi.isValid(); ++mfi) {
                if (m_use_hydroutils) {
                    // As for any other conservative scalar, with the boundary conditions in the
                    // ghost cells of F. PLIC reads neither divu nor a forcing.
                    HydroUtils::ComputeFluxesOnBoxFromState(
                        mfi.validbox(), 1, mfi, m_F[lev].const_array(mfi),
                        AMREX_D_DECL(flux[lev][0].array(mfi), flux[lev][1].array(mfi), flux[lev][2].array(mfi)),
                        AMREX_D_DECL(face[lev][0].array(mfi), face[lev][1].array(mfi), face[lev][2].array(mfi)),
                        false,
                        AMREX_D_DECL(umac[lev][0].const_array(mfi), umac[lev][1].const_array(mfi),
                                     umac[lev][2].const_array(mfi)),
                        Array4<Real const>{}, Array4<Real const>{}, Geom(lev), dt,
                        m_bcrec, m_bcrec_d.data(), m_iconserv_d.data(),
                        false, false, false, false, "PLIC");
                    continue;
                }
                PLIC::ComputeEdgeStateAndFluxes(
                    mfi.validbox(),
                    AMREX_D_DECL(flux[lev][0].array(mfi), flux[lev][1].array(mfi), flux[lev][2].array(mfi)),
                    AMREX_D_DECL(face[lev][0].array(mfi), face[lev][1].array(mfi), face[lev][2].array(mfi)),
                    m_F[lev].const_array(mfi),
                    AMREX_D_DECL(umac[lev][0].const_array(mfi), umac[lev][1].const_array(mfi),
                                 umac[lev][2].const_array(mfi)),
                    Geom(lev), dt, m_reconstruction);
            }
        }
        for (int lev = nlevs - 1; lev > 0; --lev) {
            average_down_faces(GetArrOfConstPtrs(flux[lev]), GetArrOfPtrs(flux[lev - 1]), refRatio(lev - 1),
                               Geom(lev - 1));
        }
        for (int lev = 0; lev < nlevs; ++lev) {
            const auto dxinv = Geom(lev).InvCellSizeArray();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(m_F[lev]); mfi.isValid(); ++mfi) {
                auto const& Fa = m_F[lev].array(mfi);
                AMREX_D_TERM(auto const& fx = flux[lev][0].const_array(mfi);,
                             auto const& fy = flux[lev][1].const_array(mfi);,
                             auto const& fz = flux[lev][2].const_array(mfi););
                ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Fa(i, j, k) -= dt * (AMREX_D_TERM((fx(i + 1, j, k) - fx(i, j, k)) * dxinv[0],
                                                      + (fy(i, j + 1, k) - fy(i, j, k)) * dxinv[1],
                                                      + (fz(i, j, k + 1) - fz(i, j, k)) * dxinv[2]));
                });
            }
        }
        Gpu::streamSynchronize();
        for (int lev = nlevs - 1; lev > 0; --lev) {
            average_down(m_F[lev], m_F[lev - 1], 0, 1, refRatio(lev - 1));
        }
        for (int lev = 0; lev < nlevs; ++lev) {
            FillLevel(lev, m_F[lev]);
        }
    }

    // The ghost cells: periodic, coarse/fine, physical, and then the wall ghost cells
    void FillLevel (int lev, MultiFab& dst, Real time = 0.0_rt)
    {
        const IntVect ng(dst.nGrow());
        auto bc = make_bndry(Geom(lev), m_bcrec);
        if (lev == 0) {
            FillPatchSingleLevel(dst, ng, time, {&m_F[0]}, {time}, 0, 0, 1, Geom(0), bc, 0);
        } else {
            auto cbc = make_bndry(Geom(lev - 1), m_bcrec);
            FillPatchTwoLevels(dst, ng, time, {&m_F[lev - 1]}, {time}, {&m_F[lev]}, {time}, 0, 0, 1,
                               Geom(lev - 1), Geom(lev), cbc, 0, bc, 0, refRatio(lev - 1),
                               &cell_cons_interp, m_bcrec, 0);
        }
        dst.EnforcePeriodicity(Geom(lev).periodicity());
        PLIC::FillWallGhosts(dst, Geom(lev), m_bcrec[0], m_reconstruction);
    }

protected:
    void MakeNewLevelFromScratch (int lev, Real /*time*/, const BoxArray& ba, const DistributionMapping& dm) override
    {
        m_F[lev].define(ba, dm, 1, PLIC::nGrowF(m_reconstruction));
        const auto geomdata = Geom(lev).data();
        const auto prob = m_prob;
        auto const& F = m_F[lev].arrays();
        ParallelFor(m_F[lev], [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept {
            F[b](i, j, k) = initial_fraction(i, j, k, geomdata, prob);
        });
        Gpu::streamSynchronize();
        FillLevel(lev, m_F[lev]);
    }

    void MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba, const DistributionMapping& dm) override
    {
        m_F[lev].define(ba, dm, 1, PLIC::nGrowF(m_reconstruction));
        auto cbc = make_bndry(Geom(lev - 1), m_bcrec);
        auto bc = make_bndry(Geom(lev), m_bcrec);
        InterpFromCoarseLevel(m_F[lev], IntVect(m_F[lev].nGrow()), time, m_F[lev - 1], 0, 0, 1, Geom(lev - 1),
                              Geom(lev), cbc, 0, bc, 0, refRatio(lev - 1), &cell_cons_interp, m_bcrec, 0);
        PLIC::FillWallGhosts(m_F[lev], Geom(lev), m_bcrec[0], m_reconstruction);
    }

    void RemakeLevel (int lev, Real /*time*/, const BoxArray& ba, const DistributionMapping& dm) override
    {
        MultiFab tmp(ba, dm, 1, PLIC::nGrowF(m_reconstruction));
        FillLevel(lev, tmp);
        std::swap(tmp, m_F[lev]);
    }

    void ClearLevel (int lev) override { m_F[lev].clear(); }

    void ErrorEst (int lev, TagBoxArray& tags, Real /*time*/, int /*ngrow*/) override
    {
        if (lev < max_level) { testutils::TagInterface(tags, m_F[lev], m_refine_band); }
    }

private:
    std::string m_reconstruction = "myc_height";
    ProbParm m_prob;
    Vector<MultiFab> m_F;
    Vector<BCRec> m_bcrec;
    Gpu::DeviceVector<BCRec> m_bcrec_d;
    Gpu::DeviceVector<int> m_iconserv_d;
    bool m_use_hydroutils = false;
    int m_regrid_int = 4;
    int m_refine_band = 2;
    int m_max_step = 1000000;
    Real m_fixed_dt = -1.0_rt;
    Real m_cfl = 0.45_rt;
    Real m_stop_time = 4.0_rt;
    Real m_max_l1 = -1.0_rt;
    Real m_max_volume_change = -1.0_rt;
    Real m_max_f_excursion = -1.0_rt;
    std::string m_plot_file;
};

// ---------------------------------------------------------------------------------------------

template <typename Normal>
void check_planar_fit (std::string const& name, bool exact_normal)
{
    // A planar interface on a periodic-free box, with the plane well inside it
    const int n = 16;
    Box domain(IntVect(0), IntVect(n - 1));
    RealBox rb({AMREX_D_DECL(0.0_rt, 0.0_rt, 0.0_rt)}, {AMREX_D_DECL(1.0_rt, 1.0_rt, 1.0_rt)});
    Geometry geom(domain, rb, 0, {AMREX_D_DECL(0, 0, 0)});
    const auto geomdata = geom.data();
    const Box bx = grow(domain, 4);
    FArrayBox fab(bx, 1);
    auto const& F = fab.array();
    const Real nrm[AMREX_SPACEDIM] = {AMREX_D_DECL(0.3_rt, 1.0_rt, 0.2_rt)};
    const Real d = AMREX_D_TERM(0.5_rt * nrm[0], + 0.5_rt * nrm[1], + 0.5_rt * nrm[2]);
    Real exact[3] = {0.0_rt, 0.0_rt, 0.0_rt};
    Real norm = 0.0_rt;
    for (int a = 0; a < AMREX_SPACEDIM; ++a) { exact[a] = nrm[a] * geom.CellSize(a); norm += std::abs(exact[a]); }
    for (Real& e : exact) { e /= norm; }
    amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
        F(i, j, k) = PLIC::plane_fraction(i, j, k, geomdata, nrm, d);
    });
    Array4<Real const> const Fc = fab.const_array();

    Real err_volume = 0.0_rt;
    Real err_normal = 0.0_rt;
    amrex::LoopOnCpu(grow(domain, -1), [&] (int i, int j, int k) {
        const Real f = Fc(i, j, k);
        if (PLIC::is_mixed(f)) {
            Real m[3], alpha;
            PLIC::fit_plane<Normal>(i, j, k, Fc, m[0], m[1], m[2], alpha);
            err_volume = max(err_volume, std::abs(PLIC::cut_volume(m[0], m[1], m[2], alpha, 0.0_rt, 1.0_rt) - f));
            for (int a = 0; a < 3; ++a) { err_normal = max(err_normal, std::abs(m[a] - exact[a])); }
        }
    });
    Print() << "  " << name << ": largest normal error " << err_normal << "\n";
    check(err_volume < roundoff_tol, name + " plane holds the volume of its cell");
    if (exact_normal) { check(err_normal < roundoff_tol, name + " normal is exact for a plane"); }
}

// FitPlanes on a level of several boxes: the planes of the mixed cells, valid and in the first ghost
// layer, hold the volume of their cell, and the other cells are zeroed
void check_fit_planes ()
{
    const int n = 32;
    Box domain(IntVect(0), IntVect(n - 1));
    RealBox rb({AMREX_D_DECL(0.0_rt, 0.0_rt, 0.0_rt)}, {AMREX_D_DECL(1.0_rt, 1.0_rt, 1.0_rt)});
    Geometry geom(domain, rb, 0, {AMREX_D_DECL(1, 1, 1)});
    PLIC::CheckMesh(MultiFab(), geom);
    const auto geomdata = geom.data();
    BoxArray ba(domain);
    ba.maxSize(8);
    DistributionMapping dm(ba);
    for (const char* name : {"youngs", "myc", "myc_height"}) {
        const std::string reconstruction(name);
        const int ng = PLIC::nGrowF(reconstruction);
        MultiFab F(ba, dm, 1, ng);
        const GpuArray<Real, AMREX_SPACEDIM> nrm = {AMREX_D_DECL(0.7_rt, -0.4_rt, 0.2_rt)};
        const Real d = 0.1_rt;
        auto const& fm = F.arrays();
        ParallelFor(F, IntVect(ng), [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept {
            fm[b](i, j, k) = PLIC::plane_fraction(i, j, k, geomdata, nrm.data(), d);
        });
        Gpu::streamSynchronize();
        MultiFab plane(ba, dm, 4, 1);
        plane.setVal(-1.0_rt);
        PLIC::FitPlanes(plane, F, 1, reconstruction);
        Real err = 0.0_rt;
        for (MFIter mfi(F); mfi.isValid(); ++mfi) {
            auto const& f = F.const_array(mfi);
            auto const& p = plane.const_array(mfi);
            amrex::LoopOnCpu(grow(mfi.validbox(), 1), [&] (int i, int j, int k) {
                const Real v = f(i, j, k);
                const bool mixed = PLIC::is_mixed(v);
                const Real got = mixed ? PLIC::cut_volume(p(i, j, k, 0), p(i, j, k, 1), p(i, j, k, 2),
                                                          p(i, j, k, 3), 0.0_rt, 1.0_rt)
                                       : std::abs(p(i, j, k, 0)) + std::abs(p(i, j, k, 1))
                                         + std::abs(p(i, j, k, 2)) + std::abs(p(i, j, k, 3));
                err = max(err, std::abs(got - (mixed ? v : 0.0_rt)));
            });
        }
        ParallelDescriptor::ReduceRealMax(err);
        check(err < roundoff_tol, "FitPlanes with " + reconstruction + ": planes hold the volume, pure cells zeroed");
    }
}

void check_kernels ()
{
    Print() << "kernels\n";

    // volume_intercept and cut_volume invert each other, including grid-aligned planes
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> uni(0.0, 1.0);
    Real err = 0.0_rt;
    for (int t = 0; t < 200000; ++t) {
        Real m[3] = {static_cast<Real>(uni(gen)), static_cast<Real>(uni(gen)), static_cast<Real>(uni(gen))};
        if (t % 5 == 0) { m[t % 3] = 0.0_rt; }
        if (t % 25 == 0) { m[(t + 1) % 3] = 0.0_rt; }
        const Real sum = m[0] + m[1] + m[2];
        for (Real& a : m) { a /= sum; }
        const Real V = static_cast<Real>(uni(gen));
        const Real alpha = PLIC::volume_intercept(m[0], m[1], m[2], V);
        err = max(err, std::abs(PLIC::cut_volume(m[0], m[1], m[2], alpha, 0.0_rt, 1.0_rt) - V));
    }
    Print() << "  largest volume -> plane -> volume error " << err << "\n";
    check(err < roundoff_tol, "volume_intercept and cut_volume invert each other");

    check_planar_fit<PLIC::Youngs>("youngs", false);
    check_planar_fit<PLIC::MixedYoungsCentered>("myc", false);
    check_planar_fit<PLIC::HeightFunctionMYC>("myc_height", true);

    check_fit_planes();

    // Boxes, ghost cells and the wall fill: a planar interface meeting the low and high y walls,
    // in a domain periodic in the other directions, must continue as the same plane
    {
        const int n = 32;
        Box domain(IntVect(0), IntVect(n - 1));
        RealBox rb({AMREX_D_DECL(0.0_rt, 0.0_rt, 0.0_rt)}, {AMREX_D_DECL(1.0_rt, 1.0_rt, 1.0_rt)});
        Geometry geom(domain, rb, 0, {AMREX_D_DECL(0, 0, 0)});
        const auto geomdata = geom.data();
        BoxArray ba(domain);
        ba.maxSize(8);
        DistributionMapping dm(ba);
        const std::string reconstruction = "myc_height";
        const int ng = PLIC::nGrowF(reconstruction);
        MultiFab F(ba, dm, 1, ng);
        const Real nrm[AMREX_SPACEDIM] = {AMREX_D_DECL(1.0_rt, 0.4_rt, 0.1_rt)};
        const Real d = AMREX_D_TERM(0.5_rt * nrm[0], + 0.5_rt * nrm[1], + 0.5_rt * nrm[2]);
        auto const& fm = F.arrays();
        const GpuArray<Real, AMREX_SPACEDIM> nrm_g = {AMREX_D_DECL(nrm[0], nrm[1], nrm[2])};
        ParallelFor(F, IntVect(ng), [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept {
            fm[b](i, j, k) = PLIC::plane_fraction(i, j, k, geomdata, nrm_g.data(), d);
        });
        Gpu::streamSynchronize();
        // What a solver starts from in the wall ghost cells: extrapolated values
        for (MFIter mfi(F); mfi.isValid(); ++mfi) {
            auto const& a = F.array(mfi);
            amrex::LoopOnCpu(F[mfi].box(), [&] (int i, int j, int k) {
                if (!domain.contains(IntVect(AMREX_D_DECL(i, j, k)))) {
                    a(i, j, k) = PLIC::plane_fraction(Clamp(i, 0, n - 1), Clamp(j, 0, n - 1),
                                                      AMREX_D_PICK(k, k, Clamp(k, 0, n - 1)),
                                                      geomdata, nrm_g.data(), d);
                }
            });
        }
        F.FillBoundary(geom.periodicity());

        // Inflow and symmetry faces keep the ghost values of the solver: ext_dir on the low x side,
        // reflect_even on the high x side, and walls elsewhere
        {
            MultiFab Fkeep(ba, dm, 1, ng);
            MultiFab::Copy(Fkeep, F, 0, 0, 1, ng);
            BCRec bc_keep;
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                bc_keep.setLo(dir, BCType::foextrap);
                bc_keep.setHi(dir, BCType::foextrap);
            }
            bc_keep.setLo(0, BCType::ext_dir);
            bc_keep.setHi(0, BCType::reflect_even);
            PLIC::FillWallGhosts(Fkeep, geom, bc_keep, reconstruction);
            Real keep_err = 0.0_rt;
            for (MFIter mfi(F); mfi.isValid(); ++mfi) {
                auto const& a = F.const_array(mfi);
                auto const& b = Fkeep.const_array(mfi);
                const Box gbx = F[mfi].box();
                Box lo = gbx; lo.setBig(0, -1);
                Box hi = gbx; hi.setSmall(0, n);
                for (Box const& x : {lo, hi}) {
                    if (x.ok()) {
                        amrex::LoopOnCpu(x, [&] (int i, int j, int k) {
                            keep_err = max(keep_err, std::abs(a(i, j, k) - b(i, j, k)));
                        });
                    }
                }
            }
            ParallelDescriptor::ReduceRealMax(keep_err);
            Print() << "  largest change beyond ext_dir and reflect_even faces " << keep_err << "\n";
            check(keep_err == 0.0_rt, "ext_dir and reflect_even ghost cells are kept");
        }

        BCRec bc_wall;
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            bc_wall.setLo(dir, BCType::foextrap);
            bc_wall.setHi(dir, BCType::foextrap);
        }
        PLIC::FillWallGhosts(F, geom, bc_wall, reconstruction);

        // Beyond a wall face, away from the other walls: the stencils near an edge read the ghost
        // cells of the other wall, and the search along the wall stops at it
        Real ghost_err = 0.0_rt;
        for (MFIter mfi(F); mfi.isValid(); ++mfi) {
            auto const& a = F.const_array(mfi);
            const Box gbx = F[mfi].box();
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                Box inner = grow(domain, -(ng + PLIC::wall_nsearch));
                inner.setRange(dir, -ng, n + 2 * ng);
                Box lo = gbx & inner; lo.setBig(dir, -1);
                Box hi = gbx & inner; hi.setSmall(dir, n);
                for (Box const& b : {lo, hi}) {
                    if (b.ok()) {
                        amrex::LoopOnCpu(b, [&] (int i, int j, int k) {
                            ghost_err = max(ghost_err, std::abs(a(i, j, k) - PLIC::plane_fraction(i, j, k, geomdata, nrm_g.data(), d)));
                        });
                    }
                }
            }
        }
        ParallelDescriptor::ReduceRealMax(ghost_err);
        Print() << "  largest wall ghost error for a plane " << ghost_err << "\n";
        check(ghost_err < 100.0_rt * roundoff_tol, "wall ghost cells continue a planar interface");

        // The face density is a discrete gradient for gravity normal to the plane: its curl vanishes
        Array<MultiFab, AMREX_SPACEDIM> rho;
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            rho[dir].define(convert(ba, IntVect::TheDimensionVector(dir)), dm, 1, 1);
        }
        const Real rho_l = 1000.0_rt;
        const Real rho_g = 1.0_rt;
        PLIC::SharpFaceDensity(GetArrOfPtrs(rho), {&F}, geom, rho_l, rho_g, reconstruction);
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) { rho[dir].FillBoundary(geom.periodicity()); }
        Real curl = 0.0_rt;
        const Real h = geom.CellSize(0);
        for (MFIter mfi(F); mfi.isValid(); ++mfi) {
            auto const& rx = rho[0].const_array(mfi);
            auto const& ry = rho[1].const_array(mfi);
            Box nb = surroundingNodes(mfi.validbox(), 0);
            nb.surroundingNodes(1);
            nb &= surroundingNodes(surroundingNodes(grow(domain, -6), 0), 1);
            if (nb.ok()) {
                amrex::LoopOnCpu(nb, [&] (int i, int j, int k) {
                    // d(rho_x g_x)/dy - d(rho_y g_y)/dx around a node of the xy plane
                    const Real c = (nrm[0] * (rx(i, j, k) - rx(i, j - 1, k)) - nrm[1] * (ry(i, j, k) - ry(i - 1, j, k))) / h;
                    curl = max(curl, std::abs(c));
                });
            }
        }
        ParallelDescriptor::ReduceRealMax(curl);
        Print() << "  largest curl of rho_face g, relative to rho_l g / dx " << curl * h / rho_l << "\n";
        check(curl * h / rho_l < 100.0_rt * roundoff_tol, "face density balances gravity normal to a plane");
    }
}

} // namespace

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        ParmParse pp("test");
        std::string kind = "advect";
        pp.query("kind", kind);
        if (kind == "kernels") {
            check_kernels();
        } else if (kind == "advect") {
            AdvectTest test;
            test.Run();
        } else {
            Abort("test.kind must be kernels or advect");
        }
        Print() << ((n_failures == 0) ? "PASSED" : "FAILED") << "\n";
    }
    const int failures = n_failures;
    amrex::Finalize();
    return (failures == 0) ? 0 : 1;
}

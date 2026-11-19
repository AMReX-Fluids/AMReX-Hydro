// Advection: a prescribed velocity on an AMR hierarchy (see inputs.*), advanced with
// PLIC::ComputeEdgeStateAndFluxes in the solver's own MFIter loop

#include "test_util.H"

#include <PLICTestUtils.H>

#include <AMReX_AmrCore.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>
#include <AMReX_TagBox.H>

#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <utility>

namespace plic_test {

namespace {

enum class Prob { Zalesak, Vortex };

struct ProbParm
{
    Prob type = Prob::Vortex;
    Real radius = 0.15_rt;
    GpuArray<Real, AMREX_SPACEDIM> center = {AMREX_D_DECL(0.5_rt, 0.75_rt, 0.5_rt)};
    Real slot_width = 0.05_rt;
    Real slot_length = 0.25_rt;
    Real period = 4.0_rt;
};

// The vector potential of the prescribed velocity. Back to the initial shape at t = period. It is
// its value at t = 0 times time_factor(t), so that the velocity is too.
struct Potential
{
    ProbParm p;

    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    Real time_factor (Real time) const noexcept
    {
        return (p.type == Prob::Vortex) ? std::cos(Math::pi<Real>() * time / p.period) : 1.0_rt;
    }

    AMREX_GPU_HOST_DEVICE
    void operator() (const Real* x, Real time, Real* A) const noexcept
    {
        constexpr Real pi = Math::pi<Real>();
        A[0] = A[1] = A[2] = 0.0_rt;
        if (p.type == Prob::Zalesak) {
            const Real omega = 2.0_rt * pi / p.period;
            A[2] = -0.5_rt * omega * ((x[0] - 0.5_rt) * (x[0] - 0.5_rt) + (x[1] - 0.5_rt) * (x[1] - 0.5_rt));
        } else {
            const Real ct = time_factor(time);
            const Real sx2 = std::sin(pi * x[0]) * std::sin(pi * x[0]);
            const Real sy2 = std::sin(pi * x[1]) * std::sin(pi * x[1]);
#if AMREX_SPACEDIM == 2
            A[2] = sx2 * sy2 / pi * ct;
#else
            const Real sz2 = std::sin(pi * x[2]) * std::sin(pi * x[2]);
            A[1] = -sx2 * sz2 / pi * std::sin(2.0_rt * pi * x[1]) * ct;
            A[2] = sx2 * sy2 / pi * std::sin(2.0_rt * pi * x[2]) * ct;
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
        } else {
            Abort("prob.type must be zalesak or vortex");
        }
        pp.query("radius", m_prob.radius);
        pp.query("period", m_prob.period);
        pp.query("slot_width", m_prob.slot_width);
        pp.query("slot_length", m_prob.slot_length);
        Vector<Real> v;
        if (pp.queryarr("center", v)) { for (int d = 0; d < AMREX_SPACEDIM; ++d) { m_prob.center[d] = v[d]; } }

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

        // All the non-periodic domain faces are walls: extrapolated, which FillWallGhosts then extends
        m_bcrec.resize(1);
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            const bool periodic = Geom(0).isPeriodic(d);
            m_bcrec[0].setLo(d, periodic ? BCType::int_dir : BCType::foextrap);
            m_bcrec[0].setHi(d, periodic ? BCType::int_dir : BCType::foextrap);
        }

        m_F.resize(max_level + 1);
        m_umac0.resize(max_level + 1);
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

        Print() << "  reconstruction " << m_reconstruction << "\n";

        const Real vol0 = Volume();
        Real vol_now = vol0;
        Real max_change = 0.0_rt;
        Real time = 0.0_rt;
        Real prev_dt = 0.0_rt;
        int nstep = 0;
        while (nstep < m_max_step && time < m_stop_time * (1.0_rt - 1.0e-6_rt)) {
            if (m_regrid_int > 0 && max_level > 0 && nstep > 0 && nstep % m_regrid_int == 0) {
                regrid(0, time);
                vol_now = Volume();
            }
            const Real dt = ComputeDt(time, prev_dt, nstep);

            Advance(time, dt);
            const Real vol_new = Volume();
            max_change = max(max_change, std::abs(vol_new - vol_now) / max(vol_now, std::numeric_limits<Real>::min()));
            vol_now = vol_new;

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
        Real edge_out = m_edge_out;
        Real edge_mismatch = m_edge_mismatch;
        ParallelDescriptor::ReduceRealMax(edge_out);
        ParallelDescriptor::ReduceRealMax(edge_mismatch);
        const Real vol1 = vol_now;
        const Real f_out = max(0.0_rt - f_min, f_max - 1.0_rt, 0.0_rt);
        const Real l1 = ShapeError();
        Print() << std::setprecision(6)
                << "  steps " << nstep << ", time " << time << ", finest level " << finest_level << "\n"
                << "  L1 shape error " << l1 << "\n"
                << "  volume " << vol1 << ", total relative change " << (vol1 - vol0) / vol0
                << ", largest per step " << max_change << "\n"
                << "  F left [0,1] by " << f_out << "\n"
                << "  edge state left [0,1] by " << edge_out << ", largest |edge state - Clamp(flux / u, 0, 1)| "
                << edge_mismatch << "\n";
        // Exactly: an edge state the clamp misses leaves [0,1] by round-off only
        check(edge_out == 0.0_rt, "edge state within [0,1]");
        check(edge_mismatch <= roundoff_tol, "edge state = flux / u, held to [0,1], on every face with a flow");

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

    // The member functions below are public only because nvcc rejects an extended __device__
    // lambda in a private or protected member function
    [[nodiscard]] Real Volume () const
    {
        Vector<MultiFab const*> F(finest_level + 1);
        for (int lev = 0; lev <= finest_level; ++lev) { F[lev] = &m_F[lev]; }
        return volumeWeightedSum(F, 0, Geom(), refRatio());
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
        return volumeWeightedSum(d, 0, Geom(), refRatio());
    }

    [[nodiscard]] Real ComputeDt (Real time, Real prev_dt, int nstep) const
    {
        Real dt = m_fixed_dt;
        if (dt <= 0.0_rt) {
            // The cell-centered velocity at this time
            Vector<MultiFab> vel(finest_level + 1);
            for (int lev = 0; lev <= finest_level; ++lev) {
                const auto um = Velocity(lev, time);
                vel[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 0);
                average_face_to_cellcenter(vel[lev], 0, GetArrOfConstPtrs(um));
            }
            dt = testutils::EstimateDt(vel, Geom(), m_cfl);
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

    // The edge states of one step on one level, on this rank: on every face with a flow, the liquid
    // fraction of the flux volume (the flux over the face velocity), held to [0,1]
    void CheckEdgeStates (Array<MultiFab, AMREX_SPACEDIM> const& flux, Array<MultiFab, AMREX_SPACEDIM> const& face,
                          Array<MultiFab, AMREX_SPACEDIM> const& umac)
    {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            auto const& fl = flux[d].const_arrays();
            auto const& fc = face[d].const_arrays();
            auto const& um = umac[d].const_arrays();
            auto r = ParReduce(TypeList<ReduceOpMax, ReduceOpMax>{}, TypeList<Real, Real>{}, flux[d], IntVect(0),
                               [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept -> GpuTuple<Real, Real> {
                const Real f = fc[b](i, j, k);
                const Real u = um[b](i, j, k);
                const Real out = max(-f, f - 1.0_rt, 0.0_rt);
                const Real mismatch = (u != 0.0_rt) ? std::abs(f - Clamp(fl[b](i, j, k) / u, 0.0_rt, 1.0_rt)) : 0.0_rt;
                return {nan_as_max(out), nan_as_max(mismatch)};
            });
            m_edge_out = max(m_edge_out, amrex::get<0>(r));
            m_edge_mismatch = max(m_edge_mismatch, amrex::get<1>(r));
        }
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
                flux[lev][d].define(fba, dmap[lev], 1, 0);
                face[lev][d].define(fba, dmap[lev], 1, 0);
            }
            // umac at t^{n+1/2}, with a ghost layer
            umac[lev] = Velocity(lev, time + 0.5_rt * dt);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(m_F[lev]); mfi.isValid(); ++mfi) {
                PLIC::ComputeEdgeStateAndFluxes(
                    mfi.validbox(),
                    AMREX_D_DECL(flux[lev][0].array(mfi), flux[lev][1].array(mfi), flux[lev][2].array(mfi)),
                    AMREX_D_DECL(face[lev][0].array(mfi), face[lev][1].array(mfi), face[lev][2].array(mfi)),
                    m_F[lev].const_array(mfi),
                    AMREX_D_DECL(umac[lev][0].const_array(mfi), umac[lev][1].const_array(mfi),
                                 umac[lev][2].const_array(mfi)),
                    Geom(lev), dt, m_reconstruction);
            }
            CheckEdgeStates(flux[lev], face[lev], umac[lev]);
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
        for (int lev = nlevs - 1; lev > 0; --lev) {
            average_down(m_F[lev], m_F[lev - 1], 0, 1, refRatio(lev - 1));
        }
        for (int lev = 0; lev < nlevs; ++lev) {
            FillLevel(lev, m_F[lev]);
        }
    }

    // The face velocity of a level at a time, with a ghost layer: the one at t = 0 scaled by the time
    // factor of the potential, as integrating the potential over the edges of every face is costly
    [[nodiscard]] Array<MultiFab, AMREX_SPACEDIM> Velocity (int lev, Real time) const
    {
        const Real factor = Potential{m_prob}.time_factor(time);
        Array<MultiFab, AMREX_SPACEDIM> um;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            um[d].define(m_umac0[lev][d].boxArray(), m_umac0[lev][d].DistributionMap(), 1, 1);
            MultiFab::Copy(um[d], m_umac0[lev][d], 0, 0, 1, 1);
            um[d].mult(factor, 1);
        }
        return um;
    }

    void MakeVelocity (int lev, BoxArray const& ba, DistributionMapping const& dm)
    {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            m_umac0[lev][d].define(convert(ba, IntVect::TheDimensionVector(d)), dm, 1, 1);
        }
        testutils::FaceVelocityFromPotential(GetArrOfPtrs(m_umac0[lev]), Geom(lev), 0.0_rt, Potential{m_prob});
    }

    // The ghost cells: periodic, coarse/fine, physical, and then the wall ghost cells. F is held at
    // one time only.
    void FillLevel (int lev, MultiFab& dst)
    {
        constexpr Real time = 0.0_rt;
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

    void MakeNewLevelFromScratch (int lev, Real /*time*/, const BoxArray& ba, const DistributionMapping& dm) override
    {
        MakeVelocity(lev, ba, dm);
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
        MakeVelocity(lev, ba, dm);
        m_F[lev].define(ba, dm, 1, PLIC::nGrowF(m_reconstruction));
        auto cbc = make_bndry(Geom(lev - 1), m_bcrec);
        auto bc = make_bndry(Geom(lev), m_bcrec);
        InterpFromCoarseLevel(m_F[lev], IntVect(m_F[lev].nGrow()), time, m_F[lev - 1], 0, 0, 1, Geom(lev - 1),
                              Geom(lev), cbc, 0, bc, 0, refRatio(lev - 1), &cell_cons_interp, m_bcrec, 0);
        PLIC::FillWallGhosts(m_F[lev], Geom(lev), m_bcrec[0], m_reconstruction);
    }

    void RemakeLevel (int lev, Real /*time*/, const BoxArray& ba, const DistributionMapping& dm) override
    {
        MakeVelocity(lev, ba, dm);
        MultiFab tmp(ba, dm, 1, PLIC::nGrowF(m_reconstruction));
        FillLevel(lev, tmp);
        std::swap(tmp, m_F[lev]);
    }

    void ClearLevel (int lev) override
    {
        m_F[lev].clear();
        for (auto& u : m_umac0[lev]) { u.clear(); }
    }

    void ErrorEst (int lev, TagBoxArray& tags, Real /*time*/, int /*ngrow*/) override
    {
        if (lev < max_level) { testutils::TagInterface(tags, m_F[lev], m_refine_band); }
    }

private:
    std::string m_reconstruction = "myc_height";
    ProbParm m_prob;
    Vector<MultiFab> m_F;
    // The face velocity at t = 0 on each level, with a ghost layer (see Velocity)
    Vector<Array<MultiFab, AMREX_SPACEDIM>> m_umac0;
    Vector<BCRec> m_bcrec;
    int m_regrid_int = 4;
    int m_refine_band = 2;
    int m_max_step = 1000000;
    Real m_fixed_dt = -1.0_rt;
    Real m_cfl = 0.45_rt;
    Real m_stop_time = 4.0_rt;
    Real m_max_l1 = -1.0_rt;
    Real m_max_volume_change = -1.0_rt;
    Real m_max_f_excursion = -1.0_rt;
    // On this rank
    Real m_edge_out = 0.0_rt;
    Real m_edge_mismatch = 0.0_rt;
    std::string m_plot_file;
};

} // namespace

void check_advect ()
{
    AdvectTest test;
    test.Run();
}

} // namespace plic_test

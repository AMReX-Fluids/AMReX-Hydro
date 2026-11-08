// Verification of the PLIC library. The run fails (non-zero exit) when a check does not hold.
//
//   test.kind = kernels   the plane geometry, the normals, the plane fit, the wall ghost cells and
//                         the face density

#include <hydro_plic.H>
#include <hydro_plic_normals_K.H>

#include <AMReX.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <random>
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
        std::string kind = "kernels";
        pp.query("kind", kind);
        if (kind == "kernels") {
            check_kernels();
        } else {
            Abort("test.kind must be kernels");
        }
        Print() << ((n_failures == 0) ? "PASSED" : "FAILED") << "\n";
    }
    const int failures = n_failures;
    amrex::Finalize();
    return (failures == 0) ? 0 : 1;
}

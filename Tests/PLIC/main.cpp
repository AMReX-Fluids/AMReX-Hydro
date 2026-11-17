// Verification of the PLIC library. The run fails (non-zero exit) when a check does not hold.
//
//   test.kind = kernels   the plane geometry, the normals, the plane fit, the wall ghost cells and
//                         the face density

#include <hydro_plic.H>
#include <hydro_plic_normals_K.H>

#include <AMReX.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Math.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <cmath>
#include <limits>
#include <random>
#include <string>
#include <utility>

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

// A NaN counts as the largest error
AMREX_GPU_HOST_DEVICE Real nan_as_max (Real e) noexcept
{
    return amrex::isnan(e) ? std::numeric_limits<Real>::max() : e;
}

// Largest volume and normal errors of the planes fitted to the interface n.x = n.(1/2, 1/2, 1/2)
template <typename Normal>
std::pair<Real, Real> planar_fit_errors (GpuArray<Real, AMREX_SPACEDIM> const& nrm)
{
    const int n = 16;
    Box domain(IntVect(0), IntVect(n - 1));
    RealBox rb({AMREX_D_DECL(0.0_rt, 0.0_rt, 0.0_rt)}, {AMREX_D_DECL(1.0_rt, 1.0_rt, 1.0_rt)});
    Geometry geom(domain, rb, 0, {AMREX_D_DECL(0, 0, 0)});
    const auto geomdata = geom.data();
    const Box bx = grow(domain, 4);
    FArrayBox fab(bx, 1);
    auto const& F = fab.array();
    const Real d = AMREX_D_TERM(0.5_rt * nrm[0], + 0.5_rt * nrm[1], + 0.5_rt * nrm[2]);
    Real exact[3] = {0.0_rt, 0.0_rt, 0.0_rt};
    Real norm = 0.0_rt;
    for (int a = 0; a < AMREX_SPACEDIM; ++a) { exact[a] = nrm[a] * geom.CellSize(a); norm += std::abs(exact[a]); }
    for (Real& e : exact) { e /= norm; }
    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        F(i, j, k) = PLIC::plane_fraction(i, j, k, geomdata, nrm.data(), d);
    });
    Array4<Real const> const Fc = fab.const_array();

    ReduceOps<ReduceOpMax, ReduceOpMax> reduce_op;
    ReduceData<Real, Real> reduce_data(reduce_op);
    reduce_op.eval(grow(domain, -1), reduce_data,
                   [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> GpuTuple<Real, Real> {
        const Real f = Fc(i, j, k);
        Real e_volume = 0.0_rt;
        Real e_normal = 0.0_rt;
        if (PLIC::is_mixed(f)) {
            Real m[3], alpha;
            PLIC::fit_plane<Normal>(i, j, k, Fc, m[0], m[1], m[2], alpha);
            e_volume = nan_as_max(std::abs(PLIC::cut_volume(m[0], m[1], m[2], alpha, 0.0_rt, 1.0_rt) - f));
            for (int a = 0; a < 3; ++a) { e_normal = max(e_normal, nan_as_max(std::abs(m[a] - exact[a]))); }
        }
        return {e_volume, e_normal};
    });
    const auto hv = reduce_data.value(reduce_op);
    return {amrex::get<0>(hv), amrex::get<1>(hv)};
}

template <typename Normal>
void check_planar_fit (std::string const& name, bool exact_normal)
{
    const auto [err_volume, err_normal] =
        planar_fit_errors<Normal>({AMREX_D_DECL(0.3_rt, 1.0_rt, 0.2_rt)});
    Print() << "  " << name << ": largest normal error " << err_normal << "\n";
    check(err_volume < roundoff_tol, name + " plane holds the volume of its cell");
    if (exact_normal) { check(err_normal < roundoff_tol, name + " normal is exact for a plane"); }
}

// The height-function normal in each orientation: every dominant direction, either side liquid
void check_height_orientations ()
{
    const Real base[3] = {0.3_rt, 1.0_rt, 0.2_rt};
    Real err_volume = 0.0_rt;
    Real err_normal = 0.0_rt;
    for (int perm = 0; perm < AMREX_SPACEDIM; ++perm) {
        for (int signs = 0; signs < (1 << AMREX_SPACEDIM); ++signs) {
            GpuArray<Real, AMREX_SPACEDIM> nrm;
            for (int a = 0; a < AMREX_SPACEDIM; ++a) {
                nrm[a] = ((signs >> a) & 1) ? -base[(a + perm) % 3] : base[(a + perm) % 3];
            }
            const auto [ev, en] = planar_fit_errors<PLIC::HeightFunctionMYC>(nrm);
            err_volume = max(err_volume, ev);
            err_normal = max(err_normal, en);
        }
    }
    check(err_volume < roundoff_tol && err_normal < roundoff_tol,
          "myc_height normal is exact for a plane in every orientation");
}

#if AMREX_SPACEDIM == 3
// The height-function normal over random orientations in 3D, among them those near 45 degrees with a
// tilt, where a neighboring column crosses the interface up to 2.5 cells from the row of the center
// column, which a column centered on the cell does not hold
void check_height_random_orientations ()
{
    std::mt19937 gen(5);
    std::normal_distribution<Real> normal_dist;
    Real err_volume = 0.0_rt;
    Real err_normal = 0.0_rt;
    for (int p = 0; p < 200; ++p) {
        GpuArray<Real, AMREX_SPACEDIM> nrm{normal_dist(gen), normal_dist(gen), normal_dist(gen)};
        if (p == 0) {
            const Real th = 44.0_rt * Math::pi<Real>() / 180.0_rt;
            nrm = {std::sin(th), std::cos(th), 0.3_rt * std::sin(th)};
        }
        const auto [ev, en] = planar_fit_errors<PLIC::HeightFunctionMYC>(nrm);
        err_volume = max(err_volume, ev);
        err_normal = max(err_normal, en);
    }
    Print() << "  myc_height over random orientations: largest normal error " << err_normal << "\n";
    check(err_volume < roundoff_tol && err_normal < roundoff_tol,
          "myc_height normal is exact for a plane in random orientations, 44 degrees with a tilt among them");
}
#endif

// A liquid sliver above a thin gas gap in a height column: myc_height falls back to myc
void check_height_thin_gap ()
{
    const Box bx(IntVect(-4), IntVect(4));
    for (const int col : {0, 1}) {
        FArrayBox fab(bx, 1);
        auto const& F = fab.array();
        // Interface in row y = 0, sliver in row y = 2
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            F(i, j, k) = (j < 0) ? 1.0_rt : ((j == 0) ? 0.4_rt : 0.0_rt);
            if (i == col && j == 2) { F(i, j, k) = 0.3_rt; }
        });
        Array4<Real const> const Fc = fab.const_array();
        ReduceOps<ReduceOpMax> reduce_op;
        ReduceData<Real> reduce_data(reduce_op);
        reduce_op.eval(Box(IntVect(0), IntVect(0)), reduce_data,
                       [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> GpuTuple<Real> {
            Real mh[3], mm[3];
            PLIC::HeightFunctionMYC::normal(i, j, k, Fc, mh[0], mh[1], mh[2]);
            PLIC::MixedYoungsCentered::normal(i, j, k, Fc, mm[0], mm[1], mm[2]);
            Real e = 0.0_rt;
            for (int a = 0; a < 3; ++a) { e = max(e, nan_as_max(std::abs(mh[a] - mm[a]))); }
            return {e};
        });
        const Real err = amrex::get<0>(reduce_data.value(reduce_op));
        check(err == 0.0_rt, std::string("myc_height falls back to myc for a liquid sliver in the ")
                             + (col == 0 ? "center" : "neighboring") + " height column");
    }
}

// A symmetry of a domain of n cells: swapping axes 0 and b, or with b = 0 mirroring in x. Maps a
// cell and a normal to their images.
struct Symmetry {
    const char* name;
    int b;
    AMREX_GPU_HOST_DEVICE IntVect cell (IntVect iv, int n) const noexcept {
        if (b == 0) {
            iv[0] = n - 1 - iv[0];
        } else {
            amrex::Swap(iv[0], iv[b]);
        }
        return iv;
    }
    AMREX_GPU_HOST_DEVICE void normal (Real* m) const noexcept {
        if (b == 0) {
            m[0] = -m[0];
        } else {
            amrex::Swap(m[0], m[b]);
        }
    }
};

const Symmetry symmetries[] = {
    {"swapping x and y", 1},
#if AMREX_SPACEDIM == 3
    {"swapping x and z", 2},
#endif
    {"mirroring in x", 0},
};

// The largest difference of the normal from its image under each symmetry, over spheres centered in
// the box: on the planes of symmetry, components of the normal and the candidates of MYC tie
template <typename Normal>
void check_normal_symmetry (std::string const& name)
{
    const int n = 24;
    const Box bx(IntVect(-4), IntVect(n + 3));
    FArrayBox fab(bx, 1);
    auto const& F = fab.array();
    Array4<Real const> const Fc = fab.const_array();
    // Radius and ramp width in cells; the first two broke ties by direction
    const Real spheres[][2] = {{3.37_rt, 1.0_rt}, {9.29_rt, 1.0_rt}, {9.6_rt, 2.0_rt}, {6.0_rt, 3.0_rt}};
    for (Symmetry const& sym : symmetries) {
        Real err = 0.0_rt;
        for (auto const& sp : spheres) {
            const Real radius = sp[0];
            const Real ramp = sp[1];
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                amrex::ignore_unused(k);
                const IntVect iv(AMREX_D_DECL(i, j, k));
                Real r2 = 0.0_rt;
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    const Real x = static_cast<Real>(iv[d]) + 0.5_rt - 0.5_rt * static_cast<Real>(n);
                    r2 += x * x;
                }
                F(i, j, k) = Clamp(0.5_rt + (radius - std::sqrt(r2)) / ramp, 0.0_rt, 1.0_rt);
            });
            ReduceOps<ReduceOpMax> reduce_op;
            ReduceData<Real> reduce_data(reduce_op);
            reduce_op.eval(Box(IntVect(0), IntVect(n - 1)), reduce_data,
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> GpuTuple<Real> {
                const IntVect iv(AMREX_D_DECL(i, j, k));
                if (!PLIC::is_mixed(Fc(iv))) { return {0.0_rt}; }
                Real m[3];
                Real mi[3];
                Normal::normal(i, j, k, Fc, m[0], m[1], m[2]);
                sym.normal(m);
                const auto c = sym.cell(iv, n).dim3();
                Normal::normal(c.x, c.y, c.z, Fc, mi[0], mi[1], mi[2]);
                const Real nm = std::abs(m[0]) + std::abs(m[1]) + std::abs(m[2]);
                const Real ni = std::abs(mi[0]) + std::abs(mi[1]) + std::abs(mi[2]);
                Real e = 0.0_rt;
                for (int a = 0; a < 3; ++a) { e = max(e, nan_as_max(std::abs(m[a] / nm - mi[a] / ni))); }
                return {e};
            });
            err = max(err, amrex::get<0>(reduce_data.value(reduce_op)));
        }
        Print() << "  largest difference of the " << name << " normal from its image, " << sym.name << " " << err << "\n";
        check(err < roundoff_tol, name + " normal keeps the symmetries of spheres, " + sym.name);
    }
}

// A structure one cell thick and symmetric across it, where the normals vanish: the plane lies
// across it and holds the volume of the cell. Films and gaps normal to each direction, filaments
// along each, and a drop (every direction).
template <typename Normal>
void check_degenerate_normals (std::string const& name)
{
    const Box bx(IntVect(-4), IntVect(4));
    FArrayBox fab(bx, 1);
    auto const& F = fab.array();
    Array4<Real const> const Fc = fab.const_array();
    Real err = 0.0_rt;
    // Each structure: F = 0.3 in the cells whose components along the directions in mask are 0,
    // and gas around it, or liquid for a gap
    for (int mask = 1; mask < (1 << AMREX_SPACEDIM); ++mask) {
        for (const Real outside : {0.0_rt, 1.0_rt}) {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                amrex::ignore_unused(k);
                const IntVect iv(AMREX_D_DECL(i, j, k));
                bool in = true;
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    if ((mask & (1 << d)) != 0 && iv[d] != 0) { in = false; }
                }
                F(i, j, k) = in ? 0.3_rt : outside;
            });
            ReduceOps<ReduceOpMax> reduce_op;
            ReduceData<Real> reduce_data(reduce_op);
            reduce_op.eval(Box(IntVect(0), IntVect(0)), reduce_data,
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> GpuTuple<Real> {
                Real m[3];
                Real alpha;
                PLIC::fit_plane<Normal>(i, j, k, Fc, m[0], m[1], m[2], alpha);
                // Across the structure: |m| equal along the masked directions, 0 along the others
                int nacross = 0;
                for (int d = 0; d < AMREX_SPACEDIM; ++d) { nacross += (mask >> d) & 1; }
                Real e = nan_as_max(std::abs(PLIC::cut_volume(m[0], m[1], m[2], alpha, 0.0_rt, 1.0_rt)
                                             - Fc(i, j, k)));
                for (int d = 0; d < 3; ++d) {
                    const Real expected = (d < AMREX_SPACEDIM && ((mask >> d) & 1) != 0)
                                          ? 1.0_rt / static_cast<Real>(nacross) : 0.0_rt;
                    e = max(e, nan_as_max(std::abs(std::abs(m[d]) - expected)));
                }
                return {e};
            });
            err = max(err, amrex::get<0>(reduce_data.value(reduce_op)));
        }
    }
    Print() << "  " << name << " on structures one cell thick: largest error " << err << "\n";
    check(err < roundoff_tol, name + " plane lies across films, gaps, filaments and drops one cell thick");
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
        auto const& fa = F.const_arrays();
        auto const& pa = plane.const_arrays();
        Real err = ParReduce(TypeList<ReduceOpMax>{}, TypeList<Real>{}, F, IntVect(1),
                             [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept -> GpuTuple<Real> {
            auto const& p = pa[b];
            const Real v = fa[b](i, j, k);
            const bool mixed = PLIC::is_mixed(v);
            const Real got = mixed ? PLIC::cut_volume(p(i, j, k, 0), p(i, j, k, 1), p(i, j, k, 2),
                                                      p(i, j, k, 3), 0.0_rt, 1.0_rt)
                                   : std::abs(p(i, j, k, 0)) + std::abs(p(i, j, k, 1))
                                     + std::abs(p(i, j, k, 2)) + std::abs(p(i, j, k, 3));
            return {nan_as_max(std::abs(got - (mixed ? v : 0.0_rt)))};
        });
        ParallelDescriptor::ReduceRealMax(err);
        check(err < roundoff_tol, "FitPlanes with " + reconstruction + ": planes hold the volume, pure cells zeroed");
    }
}

// A fine level covering part of the domain along the y walls: past its coarse/fine boundaries, F
// holds only nGrowF ghost layers of data (here the plane, as coarse interpolation would give), and
// the wall fill must not take the layers beyond them for data. The covered range in x is shifted so
// that the contact line on the low wall passes the coarse/fine boundary.
void check_fine_level_walls ()
{
    const int n = 32;
    Box domain(IntVect(0), IntVect(n - 1));
    RealBox rb({AMREX_D_DECL(0.0_rt, 0.0_rt, 0.0_rt)}, {AMREX_D_DECL(1.0_rt, 1.0_rt, 1.0_rt)});
    // Not periodic: the plane is not
    Geometry geom(domain, rb, 0, {AMREX_D_DECL(0, 0, 0)});
    const auto geomdata = geom.data();
    const std::string reconstruction = "myc_height";
    const int ng = PLIC::nGrowF(reconstruction);
    // Liquid towards high x, so that layers read as zeros would be a gas region
    const GpuArray<Real, AMREX_SPACEDIM> nrm_g = {AMREX_D_DECL(-1.0_rt, -0.4_rt, -0.5_rt)};
    const Real d = AMREX_D_TERM(0.5_rt * nrm_g[0], + 0.5_rt * nrm_g[1], + 0.5_rt * nrm_g[2]);
    BCRec bc;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        bc.setLo(dir, BCType::ext_dir);
        bc.setHi(dir, BCType::ext_dir);
    }
    bc.setLo(1, BCType::foextrap);
    bc.setHi(1, BCType::foextrap);

    Real err = 0.0_rt;
    for (int xhi = 15; xhi <= 25; ++xhi) {
        Box covered = domain;
        covered.setSmall(0, xhi - 15);
        covered.setBig(0, xhi);
        BoxArray ba(covered);
        ba.maxSize(8);
        DistributionMapping dm(ba);
        MultiFab F(ba, dm, 1, ng);
        auto const& fm = F.arrays();
        // The plane in the domain, extrapolated beyond the walls
        ParallelFor(F, IntVect(ng), [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept {
            fm[b](i, j, k) = PLIC::plane_fraction(i, Clamp(j, 0, n - 1), k, geomdata, nrm_g.data(), d);
        });
        Gpu::streamSynchronize();
        PLIC::FillWallGhosts(F, geom, bc, reconstruction);

        ReduceOps<ReduceOpMax> op;
        ReduceData<Real> data(op);
        for (MFIter mfi(F); mfi.isValid(); ++mfi) {
            auto const& a = F.const_array(mfi);
            Box beyond = mfi.validbox();
            beyond.grow(1, ng);
            for (Box const& b : {Box(beyond).setBig(1, -1), Box(beyond).setSmall(1, n)}) {
                op.eval(b, data, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> GpuTuple<Real> {
                    return {nan_as_max(std::abs(a(i, j, k) - PLIC::plane_fraction(i, j, k, geomdata, nrm_g.data(), d)))};
                });
            }
        }
        err = max(err, amrex::get<0>(data.value(op)));
    }
    ParallelDescriptor::ReduceRealMax(err);
    Print() << "  largest wall ghost error for a plane, fine level " << err << "\n";
    check(err < 100.0_rt * roundoff_tol, "wall ghost cells of a fine level do not read past its ghost cells");
}

// Periodic in x, with the boundary condition of x left unset (bogus), as a solver may do: the fill
// must not read it, and must still continue a plane through the y walls
void check_periodic_unset_bc ()
{
    const int n = 32;
    Box domain(IntVect(0), IntVect(n - 1));
    RealBox rb({AMREX_D_DECL(0.0_rt, 0.0_rt, 0.0_rt)}, {AMREX_D_DECL(1.0_rt, 1.0_rt, 1.0_rt)});
    Geometry geom(domain, rb, 0, {AMREX_D_DECL(1, 0, 0)});
    const auto geomdata = geom.data();
    const std::string reconstruction = "myc_height";
    const int ng = PLIC::nGrowF(reconstruction);
    // Independent of x, so periodic in x
    const GpuArray<Real, AMREX_SPACEDIM> nrm_g = {AMREX_D_DECL(0.0_rt, 0.4_rt, 1.0_rt)};
    const Real d = AMREX_D_TERM(0.5_rt * nrm_g[0], + 0.5_rt * nrm_g[1], + 0.5_rt * nrm_g[2]);
    BCRec bc;
    bc.setLo(1, BCType::foextrap);
    bc.setHi(1, BCType::foextrap);
#if AMREX_SPACEDIM == 3
    bc.setLo(2, BCType::ext_dir);
    bc.setHi(2, BCType::ext_dir);
#endif

    BoxArray ba(domain);
    ba.maxSize(8);
    DistributionMapping dm(ba);
    MultiFab F(ba, dm, 1, ng);
    auto const& fm = F.arrays();
    ParallelFor(F, IntVect(ng), [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept {
        fm[b](i, j, k) = PLIC::plane_fraction(i, Clamp(j, 0, n - 1), k, geomdata, nrm_g.data(), d);
    });
    Gpu::streamSynchronize();
    PLIC::FillWallGhosts(F, geom, bc, reconstruction);

    ReduceOps<ReduceOpMax> op;
    ReduceData<Real> data(op);
    for (MFIter mfi(F); mfi.isValid(); ++mfi) {
        auto const& a = F.const_array(mfi);
        Box beyond = mfi.validbox();
        beyond.grow(1, ng);
        for (Box const& b : {Box(beyond).setBig(1, -1), Box(beyond).setSmall(1, n)}) {
            op.eval(b, data, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> GpuTuple<Real> {
                return {nan_as_max(std::abs(a(i, j, k) - PLIC::plane_fraction(i, j, k, geomdata, nrm_g.data(), d)))};
            });
        }
    }
    Real err = amrex::get<0>(data.value(op));
    ParallelDescriptor::ReduceRealMax(err);
    Print() << "  largest wall ghost error for a plane, x periodic with an unset condition " << err << "\n";
    check(err < 100.0_rt * roundoff_tol, "the condition of a periodic direction is not read");
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
        err = max(err, nan_as_max(std::abs(PLIC::cut_volume(m[0], m[1], m[2], alpha, 0.0_rt, 1.0_rt) - V)));
    }
    Print() << "  largest volume -> plane -> volume error " << err << "\n";
    check(err < roundoff_tol, "volume_intercept and cut_volume invert each other");

    // Planes near a face diagonal: one tiny coefficient, the other two nearly equal
    err = 0.0_rt;
    for (int t = 0; t < 200000; ++t) {
        const Real eps = static_cast<Real>(std::pow(10.0, -4.0 - 10.0 * uni(gen)));
        Real m[3] = {eps * static_cast<Real>(uni(gen)), 0.5_rt,
                     0.5_rt + eps * static_cast<Real>(uni(gen) - 0.5)};
        const Real sum = m[0] + m[1] + m[2];
        for (Real& a : m) { a /= sum; }
        std::swap(m[0], m[t % 3]);
        // alpha in the thin range between the two larger coefficients
        const Real lo = min(m[(t + 1) % 3], m[(t + 2) % 3]);
        const Real width = max(m[t % 3], std::abs(m[(t + 1) % 3] - m[(t + 2) % 3]));
        const Real alpha0 = lo + 1.2_rt * width * static_cast<Real>(uni(gen));
        const Real V = PLIC::cut_volume(m[0], m[1], m[2], alpha0, 0.0_rt, 1.0_rt);
        const Real alpha = PLIC::volume_intercept(m[0], m[1], m[2], V);
        err = max(err, nan_as_max(std::abs(alpha - alpha0)));
        err = max(err, nan_as_max(std::abs(PLIC::cut_volume(m[0], m[1], m[2], alpha, 0.0_rt, 1.0_rt) - V)));
    }
    Print() << "  largest error, near-diagonal planes " << err << "\n";
    check(err < roundoff_tol, "volume_intercept and cut_volume invert each other, near-diagonal planes");

    // Two components near the underflow threshold
    {
        const Real s = 1.0e4_rt * std::numeric_limits<Real>::min();
        const Real alpha = 0.1_rt * s;
        const Real v = PLIC::cut_volume(s, s, 1.0_rt - 2.0_rt * s, alpha, 0.0_rt, 1.0_rt);
        check(std::isfinite(v) && std::abs(v - alpha) < roundoff_tol,
              "cut_volume with components near the underflow threshold");
    }

    check_planar_fit<PLIC::Youngs>("youngs", false);
    check_planar_fit<PLIC::MixedYoungsCentered>("myc", false);
    check_planar_fit<PLIC::HeightFunctionMYC>("myc_height", true);
    check_height_orientations();

#if AMREX_SPACEDIM == 3
    check_height_random_orientations();
#endif
    check_height_thin_gap();
    check_degenerate_normals<PLIC::Youngs>("youngs");
    check_degenerate_normals<PLIC::MixedYoungsCentered>("myc");
    check_degenerate_normals<PLIC::HeightFunctionMYC>("myc_height");
    check_normal_symmetry<PLIC::MixedYoungsCentered>("myc");
    check_normal_symmetry<PLIC::HeightFunctionMYC>("myc_height");
    check_fit_planes();
    check_fine_level_walls();
    check_periodic_unset_bc();

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
            ParallelFor(F[mfi].box(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
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
            ReduceOps<ReduceOpMax> keep_op;
            ReduceData<Real> keep_data(keep_op);
            for (MFIter mfi(F); mfi.isValid(); ++mfi) {
                auto const& a = F.const_array(mfi);
                auto const& b = Fkeep.const_array(mfi);
                const Box gbx = F[mfi].box();
                Box lo = gbx; lo.setBig(0, -1);
                Box hi = gbx; hi.setSmall(0, n);
                for (Box const& x : {lo, hi}) {
                    if (x.ok()) {
                        keep_op.eval(x, keep_data, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> GpuTuple<Real> {
                            return {nan_as_max(std::abs(a(i, j, k) - b(i, j, k)))};
                        });
                    }
                }
            }
            Real keep_err = amrex::get<0>(keep_data.value(keep_op));
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
        ReduceOps<ReduceOpMax> ghost_op;
        ReduceData<Real> ghost_data(ghost_op);
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
                        ghost_op.eval(b, ghost_data, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> GpuTuple<Real> {
                            return {nan_as_max(std::abs(a(i, j, k) - PLIC::plane_fraction(i, j, k, geomdata, nrm_g.data(), d)))};
                        });
                    }
                }
            }
        }
        Real ghost_err = amrex::get<0>(ghost_data.value(ghost_op));
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
        // The ghost faces must be left as they are
        const Real sentinel = -7.0_rt;
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) { rho[dir].setVal(sentinel); }
        PLIC::SharpFaceDensity(GetArrOfPtrs(rho), {&F}, geom, rho_l, rho_g, reconstruction);
        {
            ReduceOps<ReduceOpMax> gh_op;
            ReduceData<Real> gh_data(gh_op);
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                for (MFIter mfi(rho[dir]); mfi.isValid(); ++mfi) {
                    auto const& r = rho[dir].const_array(mfi);
                    const Box vbx = mfi.validbox();
                    gh_op.eval(rho[dir][mfi].box(), gh_data, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> GpuTuple<Real> {
                        const bool ghost = !vbx.contains(IntVect(AMREX_D_DECL(i, j, k)));
                        return {ghost ? nan_as_max(std::abs(r(i, j, k) - sentinel)) : 0.0_rt};
                    });
                }
            }
            Real gh_err = amrex::get<0>(gh_data.value(gh_op));
            ParallelDescriptor::ReduceRealMax(gh_err);
            check(gh_err == 0.0_rt, "SharpFaceDensity leaves the ghost faces of rho_face as they are");
        }
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) { rho[dir].FillBoundary(geom.periodicity()); }
        ReduceOps<ReduceOpMax> curl_op;
        ReduceData<Real> curl_data(curl_op);
        const Real h = geom.CellSize(0);
        for (MFIter mfi(F); mfi.isValid(); ++mfi) {
            auto const& rx = rho[0].const_array(mfi);
            auto const& ry = rho[1].const_array(mfi);
            Box nb = surroundingNodes(mfi.validbox(), 0);
            nb.surroundingNodes(1);
            nb &= surroundingNodes(surroundingNodes(grow(domain, -6), 0), 1);
            if (nb.ok()) {
                curl_op.eval(nb, curl_data, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> GpuTuple<Real> {
                    // d(rho_x g_x)/dy - d(rho_y g_y)/dx around a node of the xy plane
                    const Real c = (nrm_g[0] * (rx(i, j, k) - rx(i, j - 1, k)) - nrm_g[1] * (ry(i, j, k) - ry(i - 1, j, k))) / h;
                    return {nan_as_max(std::abs(c))};
                });
            }
        }
        Real curl = amrex::get<0>(curl_data.value(curl_op));
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

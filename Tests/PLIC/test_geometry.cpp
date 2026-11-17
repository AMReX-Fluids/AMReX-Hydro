// The plane geometry, the normals and the plane fit

#include "test_util.H"

#include <AMReX_FArrayBox.H>
#include <AMReX_Math.H>
#include <AMReX_Print.H>

#include <cmath>
#include <random>
#include <utility>

namespace plic_test {

namespace {

// volume_intercept and cut_volume invert each other, including grid-aligned planes, planes near a
// face diagonal, and components near the underflow threshold
void check_volume_relations ()
{
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
}

// Largest volume and normal errors of the planes fitted to the interface through the center of a
// box of 16^dim cells and the given lengths. The pure cells of every other cell (checkerboard)
// overshoot [0, 1] by overshoot, as advection may leave them.
template <typename Normal>
std::pair<Real, Real> planar_fit_errors (GpuArray<Real, AMREX_SPACEDIM> const& nrm,
                                         GpuArray<Real, AMREX_SPACEDIM> const& lengths
                                             = {AMREX_D_DECL(1.0_rt, 1.0_rt, 1.0_rt)},
                                         Real overshoot = 0.0_rt)
{
    const int n = 16;
    Box domain(IntVect(0), IntVect(n - 1));
    RealBox rb({AMREX_D_DECL(0.0_rt, 0.0_rt, 0.0_rt)}, {AMREX_D_DECL(lengths[0], lengths[1], lengths[2])});
    Geometry geom(domain, rb, 0, {AMREX_D_DECL(0, 0, 0)});
    const auto geomdata = geom.data();
    const Box bx = grow(domain, 4);
    FArrayBox fab(bx, 1);
    auto const& F = fab.array();
    const Real d = AMREX_D_TERM(0.5_rt * nrm[0] * lengths[0], + 0.5_rt * nrm[1] * lengths[1],
                                + 0.5_rt * nrm[2] * lengths[2]);
    Real exact[3] = {0.0_rt, 0.0_rt, 0.0_rt};
    Real norm = 0.0_rt;
    for (int a = 0; a < AMREX_SPACEDIM; ++a) { exact[a] = nrm[a] * geom.CellSize(a); norm += std::abs(exact[a]); }
    for (Real& e : exact) { e /= norm; }
    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real f = PLIC::plane_fraction(i, j, k, geomdata, nrm.data(), d);
        if (((i + j + k) & 1) != 0 && !PLIC::is_mixed(f)) {
            f += (f > 0.5_rt) ? overshoot : -overshoot;
        }
        F(i, j, k) = f;
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

// The myc_height normal over planes in every orientation and random ones, on cells of different
// lengths in each direction (aspect ratios 1:4:2): the normals work in index space, where the
// interface is still a plane
void check_height_anisotropic ()
{
    const GpuArray<Real, AMREX_SPACEDIM> lengths = {AMREX_D_DECL(1.0_rt, 4.0_rt, 2.0_rt)};
    const Real base[3] = {0.3_rt, 1.0_rt, 0.2_rt};
    std::mt19937 gen(7);
    std::normal_distribution<Real> normal_dist;
    Real err_volume = 0.0_rt;
    Real err_normal = 0.0_rt;
    for (int p = 0; p < 50 + AMREX_SPACEDIM * (1 << AMREX_SPACEDIM); ++p) {
        GpuArray<Real, AMREX_SPACEDIM> nrm;
        if (p < AMREX_SPACEDIM * (1 << AMREX_SPACEDIM)) {
            const int perm = p >> AMREX_SPACEDIM;
            const int signs = p & ((1 << AMREX_SPACEDIM) - 1);
            for (int a = 0; a < AMREX_SPACEDIM; ++a) {
                nrm[a] = ((signs >> a) & 1) ? -base[(a + perm) % 3] : base[(a + perm) % 3];
            }
        } else {
            for (int a = 0; a < AMREX_SPACEDIM; ++a) { nrm[a] = normal_dist(gen); }
        }
        const auto [ev, en] = planar_fit_errors<PLIC::HeightFunctionMYC>(nrm, lengths);
        err_volume = max(err_volume, ev);
        err_normal = max(err_normal, en);
    }
    Print() << "  myc_height on cells of aspect ratio 1:4:2: largest normal error " << err_normal << "\n";
    check(err_volume < roundoff_tol && err_normal < roundoff_tol,
          "myc_height normal is exact for a plane on cells of different lengths in each direction");
}

// Pure cells that overshoot [0, 1]: by round-off, myc_height stays exact for a plane; by 1e-3, as
// a poor limiter may leave them, every reconstruction still gives planes that hold the volume
void check_overshoot ()
{
    const Real base[3] = {0.3_rt, 1.0_rt, 0.2_rt};
    const GpuArray<Real, AMREX_SPACEDIM> nrm = {AMREX_D_DECL(base[0], base[1], base[2])};
    const GpuArray<Real, AMREX_SPACEDIM> unit = {AMREX_D_DECL(1.0_rt, 1.0_rt, 1.0_rt)};
    const Real roundoff = 100.0_rt * std::numeric_limits<Real>::epsilon();
    const auto [ev, en] = planar_fit_errors<PLIC::HeightFunctionMYC>(nrm, unit, roundoff);
    Print() << "  myc_height with pure cells overshooting by round-off: largest normal error " << en << "\n";
    check(ev < roundoff_tol && en < roundoff_tol,
          "myc_height normal is exact for a plane with pure cells overshooting by round-off");
    Real err_volume = 0.0_rt;
    Real err_normal = 0.0_rt;
    for (int p = 0; p < 3; ++p) {
        std::pair<Real, Real> e;
        if (p == 0) { e = planar_fit_errors<PLIC::Youngs>(nrm, unit, 1.0e-3_rt); }
        if (p == 1) { e = planar_fit_errors<PLIC::MixedYoungsCentered>(nrm, unit, 1.0e-3_rt); }
        if (p == 2) { e = planar_fit_errors<PLIC::HeightFunctionMYC>(nrm, unit, 1.0e-3_rt); }
        err_volume = max(err_volume, e.first);
        err_normal = max(err_normal, e.second);
    }
    Print() << "  pure cells overshooting by 1e-3: largest volume error " << err_volume
            << ", largest normal error " << err_normal << "\n";
    check(err_volume < roundoff_tol, "planes hold the volume with pure cells overshooting by 1e-3");
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

// A structure one cell thick and symmetric across it: the cells c with c.e = 0 for each of its
// constraints e, and the normal expected (up to sign) across it
struct Thin {
    const char* name;
    int nc;
    int e[3][3];
    Real normal[3];
};

const Thin thin_structures[] = {
    {"film normal to x", 1, {{1, 0, 0}}, {1, 0, 0}},
    {"film normal to y", 1, {{0, 1, 0}}, {0, 1, 0}},
    {"diagonal film", 1, {{1, -1, 0}}, {1, -1, 0}},
    {"anti-diagonal film", 1, {{1, 1, 0}}, {1, 1, 0}},
#if AMREX_SPACEDIM == 3
    {"film normal to z", 1, {{0, 0, 1}}, {0, 0, 1}},
    {"diagonal film in x and z", 1, {{1, 0, -1}}, {1, 0, -1}},
    {"film normal to (1, 1, 1)", 1, {{1, 1, 1}}, {1, 1, 1}},
    {"filament along z", 2, {{1, 0, 0}, {0, 1, 0}}, {1, 1, 0}},
    {"drop", 3, {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, {1, 1, 1}},
#else
    {"drop", 2, {{1, 0, 0}, {0, 1, 0}}, {1, 1, 0}},
#endif
};

// Where the normals vanish, on structures one cell thick (liquid in gas, and gaps in liquid), the
// plane lies across the structure and holds the volume of the cell
template <typename Normal>
void check_degenerate_normals (std::string const& name)
{
    const Box bx(IntVect(-4), IntVect(4));
    FArrayBox fab(bx, 1);
    auto const& F = fab.array();
    Array4<Real const> const Fc = fab.const_array();
    Real err = 0.0_rt;
    for (Thin const& th : thin_structures) {
        for (const Real outside : {0.0_rt, 1.0_rt}) {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                const int c[3] = {i, j, k};
                bool in = true;
                for (int q = 0; q < th.nc; ++q) {
                    in = in && (th.e[q][0] * c[0] + th.e[q][1] * c[1] + th.e[q][2] * c[2] == 0);
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
                const Real nn = std::abs(th.normal[0]) + std::abs(th.normal[1]) + std::abs(th.normal[2]);
                Real e_plus = 0.0_rt;
                Real e_minus = 0.0_rt;
                for (int d = 0; d < 3; ++d) {
                    e_plus = max(e_plus, std::abs(m[d] - th.normal[d] / nn));
                    e_minus = max(e_minus, std::abs(m[d] + th.normal[d] / nn));
                }
                const Real e_vol = std::abs(PLIC::cut_volume(m[0], m[1], m[2], alpha, 0.0_rt, 1.0_rt) - Fc(i, j, k));
                return {nan_as_max(max(min(e_plus, e_minus), e_vol))};
            });
            const Real e = amrex::get<0>(reduce_data.value(reduce_op));
            if (e >= roundoff_tol) {
                Print() << "    " << name << ", " << th.name << (outside > 0.5_rt ? " (gap)" : "") << ": " << e << "\n";
            }
            err = max(err, e);
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

} // namespace

void check_geometry ()
{
    Print() << "plane geometry, normals and plane fit\n";
    check_volume_relations();
    check_planar_fit<PLIC::Youngs>("youngs", false);
    check_planar_fit<PLIC::MixedYoungsCentered>("myc", false);
    check_planar_fit<PLIC::HeightFunctionMYC>("myc_height", true);
    check_height_orientations();
#if AMREX_SPACEDIM == 3
    check_height_random_orientations();
#endif
    check_height_anisotropic();
    check_overshoot();
    check_height_thin_gap();
    check_degenerate_normals<PLIC::Youngs>("youngs");
    check_degenerate_normals<PLIC::MixedYoungsCentered>("myc");
    check_degenerate_normals<PLIC::HeightFunctionMYC>("myc_height");
    check_normal_symmetry<PLIC::MixedYoungsCentered>("myc");
    check_normal_symmetry<PLIC::HeightFunctionMYC>("myc_height");
    check_fit_planes();
}

} // namespace plic_test

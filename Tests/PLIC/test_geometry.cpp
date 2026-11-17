// The plane geometry, the normals and the plane fit

#include "test_util.H"

#include <AMReX_FArrayBox.H>
#include <AMReX_Math.H>
#include <AMReX_Print.H>

#include <cmath>
#include <random>
#include <utility>
#include <vector>

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

using Normals = std::vector<GpuArray<Real, AMREX_SPACEDIM>>;

// Plane normals in every orientation: each dominant direction, either side liquid
Normals all_orientations ()
{
    const Real base[3] = {0.3_rt, 1.0_rt, 0.2_rt};
    Normals normals;
    for (int perm = 0; perm < AMREX_SPACEDIM; ++perm) {
        for (int signs = 0; signs < (1 << AMREX_SPACEDIM); ++signs) {
            GpuArray<Real, AMREX_SPACEDIM> nrm;
            for (int a = 0; a < AMREX_SPACEDIM; ++a) {
                nrm[a] = ((signs >> a) & 1) ? -base[(a + perm) % 3] : base[(a + perm) % 3];
            }
            normals.push_back(nrm);
        }
    }
    return normals;
}

// n random plane normals
Normals random_normals (int n, unsigned seed)
{
    std::mt19937 gen(seed);
    std::normal_distribution<Real> normal_dist;
    Normals normals(n);
    for (auto& nrm : normals) {
        for (int a = 0; a < AMREX_SPACEDIM; ++a) { nrm[a] = normal_dist(gen); }
    }
    return normals;
}

// Whether the myc_height normal is exact, and its planes hold the volume, for planes of the given
// normals on cells of the given lengths; prints the largest normal error under label
bool height_exact (Normals const& normals, GpuArray<Real, AMREX_SPACEDIM> const& lengths,
                   std::string const& label)
{
    Real err_volume = 0.0_rt;
    Real err_normal = 0.0_rt;
    for (auto const& nrm : normals) {
        const auto [ev, en] = planar_fit_errors<PLIC::HeightFunctionMYC>(nrm, lengths);
        err_volume = max(err_volume, ev);
        err_normal = max(err_normal, en);
    }
    if (!label.empty()) {
        Print() << "  myc_height " << label << ": largest normal error " << err_normal << "\n";
    }
    return err_volume < roundoff_tol && err_normal < roundoff_tol;
}

void check_height_orientations ()
{
    const GpuArray<Real, AMREX_SPACEDIM> unit = {AMREX_D_DECL(1.0_rt, 1.0_rt, 1.0_rt)};
    check(height_exact(all_orientations(), unit, ""),
          "myc_height normal is exact for a plane in every orientation");

#if AMREX_SPACEDIM == 3
    // Among random orientations, one near 45 degrees with a tilt, where a neighboring column
    // crosses the interface up to 2.5 cells from the row of the center column, which a column
    // centered on the cell does not hold
    Normals normals = random_normals(200, 5);
    const Real th = 44.0_rt * Math::pi<Real>() / 180.0_rt;
    normals[0] = {std::sin(th), std::cos(th), 0.3_rt * std::sin(th)};
    check(height_exact(normals, unit, "over random orientations"),
          "myc_height normal is exact for a plane in random orientations, 44 degrees with a tilt among them");
#endif

    // On cells of aspect ratio 1:4:2: the normals work in index space, where the interface is
    // still a plane
    Normals aniso = all_orientations();
    const Normals random = random_normals(50, 7);
    aniso.insert(aniso.end(), random.begin(), random.end());
    check(height_exact(aniso, {AMREX_D_DECL(1.0_rt, 4.0_rt, 2.0_rt)}, "on cells of aspect ratio 1:4:2"),
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
    for (auto const& e : {planar_fit_errors<PLIC::Youngs>(nrm, unit, 1.0e-3_rt),
                          planar_fit_errors<PLIC::MixedYoungsCentered>(nrm, unit, 1.0e-3_rt),
                          planar_fit_errors<PLIC::HeightFunctionMYC>(nrm, unit, 1.0e-3_rt)}) {
        err_volume = max(err_volume, e.first);
        err_normal = max(err_normal, e.second);
    }
    Print() << "  pure cells overshooting by 1e-3: largest volume error " << err_volume
            << ", largest normal error " << err_normal << "\n";
    check(err_volume < roundoff_tol, "planes hold the volume with pure cells overshooting by 1e-3");
}


// A liquid sliver above a thin gas gap in a height column: myc_height falls back to myc. The
// interface crosses a neighboring column beyond the 3x3 block of myc, where the heights differ
// from it, so that a column check that let the sliver through would show.
void check_height_thin_gap ()
{
    const Box bx(IntVect(-4), IntVect(4));
    for (const int col : {0, 1}) {
        Real err[2] = {0.0_rt, 0.0_rt};
        for (const bool sliver : {false, true}) {
            FArrayBox fab(bx, 1);
            auto const& F = fab.array();
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                // Column x = -1 crosses in row 2, beyond the 3x3 block of myc
                const Real fl = (j < 2) ? 1.0_rt : ((j == 2) ? 0.1_rt : 0.0_rt);
                const Real fr = (j < 0) ? 1.0_rt : ((j == 0) ? 0.4_rt : 0.0_rt);
                F(i, j, k) = (i < 0) ? fl : fr;
                if (sliver && i == col && j == 2) { F(i, j, k) = 0.3_rt; }
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
            err[sliver ? 1 : 0] = amrex::get<0>(reduce_data.value(reduce_op));
        }
        Print() << "  myc_height - myc without the sliver " << err[0] << ", with it " << err[1] << "\n";
        check(err[0] > 1.0e-3_rt && err[1] == 0.0_rt,
              std::string("myc_height falls back to myc for a liquid sliver in the ")
              + (col == 0 ? "center" : "neighboring") + " height column");
    }
}

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

// On structures one cell thick (liquid in gas, and gaps in liquid), in the cells of the structure
// next to the origin, the plane holds the volume of the cell and lies across the structure:
// - exactly symmetric across it, along the expected normal (to round-off);
// - with noise of 1e-3 along it (not odd about any cell, so that it reaches the Hessian) and of 1e-9
//   in the pure cells, within
//   2e-2 of the expected normal for a film, and normal to a filament; the cells of a film agree,
//   orientation included.
template <typename Normal>
void check_degenerate_normals (std::string const& name)
{
    const Box bx(IntVect(-5), IntVect(5));
    FArrayBox fab(bx, 1);
    auto const& F = fab.array();
    Array4<Real const> const Fc = fab.const_array();
    const Real noise_tol = 2.0e-2_rt;
    Real err_exact = 0.0_rt;
    Real err_noise = 0.0_rt;
    for (Thin const& th : thin_structures) {
        for (int variant = 0; variant < 4; ++variant) {
            const Real outside = (variant & 1) ? 1.0_rt : 0.0_rt;
            const Real noise = (variant & 2) ? 1.0e-3_rt : 0.0_rt;
            const auto in = [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                const int c[3] = {i, j, k};
                bool on = true;
                for (int q = 0; q < th.nc; ++q) {
                    on = on && (th.e[q][0] * c[0] + th.e[q][1] * c[1] + th.e[q][2] * c[2] == 0);
                }
                return on;
            };
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                const Real r = std::sin(12.9898_rt * Real(i) + 78.233_rt * Real(j) + 37.719_rt * Real(k) + 0.5_rt);
                // and pure cells off by round-off, so that the diagonal of the Hessian is not exact
                const Real pure = (outside > 0.5_rt) ? outside - 1.0e-9_rt * (1.0_rt + r) : 1.0e-9_rt * (1.0_rt + r);
                F(i, j, k) = in(i, j, k) ? 0.3_rt + noise * r : ((noise > 0.0_rt) ? pure : outside);
            });
            // The error, and the largest and smallest of each component of m on a film (its spread)
            ReduceOps<ReduceOpMax, ReduceOpMax, ReduceOpMax, ReduceOpMax, ReduceOpMin, ReduceOpMin, ReduceOpMin> reduce_op;
            ReduceData<Real, Real, Real, Real, Real, Real, Real> reduce_data(reduce_op);
            const Real lowest = std::numeric_limits<Real>::lowest();
            const Real highest = std::numeric_limits<Real>::max();
            reduce_op.eval(Box(IntVect(-1), IntVect(1)), reduce_data,
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                           -> GpuTuple<Real, Real, Real, Real, Real, Real, Real> {
                if (!in(i, j, k)) {
                    return {0.0_rt, lowest, lowest, lowest, highest, highest, highest};
                }
                Real m[3];
                Real alpha;
                PLIC::fit_plane<Normal>(i, j, k, Fc, m[0], m[1], m[2], alpha);
                Real e = std::abs(PLIC::cut_volume(m[0], m[1], m[2], alpha, 0.0_rt, 1.0_rt) - Fc(i, j, k));
                if (noise == 0.0_rt || th.nc == 1) {
                    const Real nn = std::abs(th.normal[0]) + std::abs(th.normal[1]) + std::abs(th.normal[2]);
                    Real e_plus = 0.0_rt;
                    Real e_minus = 0.0_rt;
                    for (int d = 0; d < 3; ++d) {
                        e_plus = max(e_plus, std::abs(m[d] - th.normal[d] / nn));
                        e_minus = max(e_minus, std::abs(m[d] + th.normal[d] / nn));
                    }
                    e = max(e, min(e_plus, e_minus));
                } else if (th.nc == 2 && AMREX_SPACEDIM == 3) {
                    // Normal to the filament, along e0 x e1
                    const int (&a)[3] = th.e[0];
                    const int (&b)[3] = th.e[1];
                    const Real ax[3] = {Real(a[1] * b[2] - a[2] * b[1]), Real(a[2] * b[0] - a[0] * b[2]),
                                        Real(a[0] * b[1] - a[1] * b[0])};
                    e = max(e, std::abs(m[0] * ax[0] + m[1] * ax[1] + m[2] * ax[2]));
                }
                const bool film = th.nc == 1;
                return {nan_as_max(e), film ? m[0] : lowest, film ? m[1] : lowest, film ? m[2] : lowest,
                        film ? m[0] : highest, film ? m[1] : highest, film ? m[2] : highest};
            });
            const auto r = reduce_data.value(reduce_op);
            Real e = amrex::get<0>(r);
            const Real hi[3] = {amrex::get<1>(r), amrex::get<2>(r), amrex::get<3>(r)};
            const Real lo[3] = {amrex::get<4>(r), amrex::get<5>(r), amrex::get<6>(r)};
            for (int d = 0; d < 3; ++d) {
                if (hi[d] > lowest) { e = max(e, hi[d] - lo[d]); }
            }
            const Real tol = (noise > 0.0_rt) ? noise_tol : roundoff_tol;
            if (e >= tol) {
                Print() << "    " << name << ", " << th.name << (outside > 0.5_rt ? " (gap)" : "")
                        << (noise > 0.0_rt ? " (noise)" : "") << ": " << e << "\n";
            }
            Real& err = (noise > 0.0_rt) ? err_noise : err_exact;
            err = max(err, e);
        }
    }
    Print() << "  " << name << " on structures one cell thick: largest error " << err_exact
            << ", with noise " << err_noise << "\n";
    check(err_exact < roundoff_tol && err_noise < noise_tol,
          name + " plane lies across films, gaps, filaments and drops one cell thick");
}

// On spheres of radius 0.6 to 6 cells (volume fractions from PLIC::cell_fraction), the
// tips of the interface are not taken for structures one cell thick: fit_plane keeps the normal
template <typename Normal>
void check_spheres_not_thin (std::string const& name)
{
    Real err = 0.0_rt;
    for (const Real R : {0.6_rt, 0.8_rt, 1.0_rt, 1.3_rt, 1.6_rt, 2.0_rt, 3.0_rt, 6.0_rt}) {
        const int n = static_cast<int>(R) + 4;
        const Box bx(IntVect(-n - 2), IntVect(n + 2));
        FArrayBox fab(bx, 1);
        auto const& F = fab.array();
        // Unit cells, the center of the sphere at c from the low corner of cell 0
        const Real c[3] = {0.31_rt, 0.17_rt, 0.07_rt};
        const Geometry geom(Box(IntVect(0), IntVect(1)),
                            RealBox({AMREX_D_DECL(-c[0], -c[1], -c[2])},
                                    {AMREX_D_DECL(2.0_rt - c[0], 2.0_rt - c[1], 2.0_rt - c[2])}),
                            0, {AMREX_D_DECL(0, 0, 0)});
        const auto geomdata = geom.data();
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            F(i, j, k) = PLIC::cell_fraction(i, j, k, geomdata, [=] AMREX_GPU_DEVICE (Real x, Real y, Real z) {
                return R * R - (x * x + y * y + (AMREX_SPACEDIM == 3 ? z * z : 0.0_rt));
            });
        });
        Array4<Real const> const Fc = fab.const_array();
        ReduceOps<ReduceOpMax> reduce_op;
        ReduceData<Real> reduce_data(reduce_op);
        reduce_op.eval(grow(bx, -Normal::stencil), reduce_data,
                       [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> GpuTuple<Real> {
            Real e = 0.0_rt;
            if (PLIC::is_mixed(Fc(i, j, k))) {
                Real m[3];
                Real mn[3];
                Real alpha;
                PLIC::fit_plane<Normal>(i, j, k, Fc, m[0], m[1], m[2], alpha);
                Normal::normal(i, j, k, Fc, mn[0], mn[1], mn[2]);
                const Real t = std::abs(mn[0]) + std::abs(mn[1]) + std::abs(mn[2]);
                for (int d = 0; d < 3; ++d) { e = max(e, nan_as_max(std::abs(m[d] - mn[d] / t))); }
            }
            return {e};
        });
        err = max(err, amrex::get<0>(reduce_data.value(reduce_op)));
    }
    check(err < roundoff_tol, name + " keeps its normal on spheres of radius 0.6 to 6 cells");
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
    check_overshoot();
    check_height_thin_gap();
    check_degenerate_normals<PLIC::Youngs>("youngs");
    check_degenerate_normals<PLIC::MixedYoungsCentered>("myc");
    check_degenerate_normals<PLIC::HeightFunctionMYC>("myc_height");
    check_spheres_not_thin<PLIC::Youngs>("youngs");
    check_spheres_not_thin<PLIC::MixedYoungsCentered>("myc");
    check_spheres_not_thin<PLIC::HeightFunctionMYC>("myc_height");
    check_normal_symmetry<PLIC::MixedYoungsCentered>("myc");
    check_normal_symmetry<PLIC::HeightFunctionMYC>("myc_height");
    check_fit_planes();
}

} // namespace plic_test

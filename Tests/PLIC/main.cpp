// Verification of the PLIC library. The run fails (non-zero exit) when a check does not hold.
//
//   test.kind = kernels   the plane geometry, the normals and the plane fit

#include <hydro_plic.H>
#include <hydro_plic_normals_K.H>

#include <AMReX.H>
#include <AMReX_FArrayBox.H>
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

    check_height_thin_gap();
    check_fit_planes();
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

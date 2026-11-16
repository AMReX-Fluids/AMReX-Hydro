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

// An error to take the max of. A NaN would drop out of a max, so it counts as the largest error.
AMREX_GPU_HOST_DEVICE Real nan_as_max (Real e) noexcept
{
    return amrex::isnan(e) ? std::numeric_limits<Real>::max() : e;
}

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
    const GpuArray<Real, AMREX_SPACEDIM> nrm = {AMREX_D_DECL(0.3_rt, 1.0_rt, 0.2_rt)};
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
    const Real err_volume = amrex::get<0>(hv);
    const Real err_normal = amrex::get<1>(hv);
    Print() << "  " << name << ": largest normal error " << err_normal << "\n";
    check(err_volume < roundoff_tol, name + " plane holds the volume of its cell");
    if (exact_normal) { check(err_normal < roundoff_tol, name + " normal is exact for a plane"); }
}

// A liquid sliver above a thin gas gap, in the center height column or in a neighboring one: the
// column crosses the interface three times, so it holds no height, and myc_height must give the myc
// normal instead
void check_height_thin_gap ()
{
    const Box bx(IntVect(-4), IntVect(4));
    for (const int col : {0, 1}) {
        FArrayBox fab(bx, 1);
        auto const& F = fab.array();
        // Liquid below a flat interface in cell row y = 0, F = 0.4 there. Above it gas at y = 1, a
        // liquid sliver at y = 2 and gas at y = 3, so the column still ends pure.
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

    // The same for a plane close to a diagonal of a cell face: one tiny coefficient, the other two
    // within it of each other, and the plane through the thin range of alpha between the two larger
    // coefficients. A diagonal interface with a little noise in F gives this.
    err = 0.0_rt;
    for (int t = 0; t < 200000; ++t) {
        const Real eps = static_cast<Real>(std::pow(10.0, -4.0 - 10.0 * uni(gen)));
        Real m[3] = {eps * static_cast<Real>(uni(gen)), 0.5_rt,
                     0.5_rt + eps * static_cast<Real>(uni(gen) - 0.5)};
        const Real sum = m[0] + m[1] + m[2];
        for (Real& a : m) { a /= sum; }
        std::swap(m[0], m[t % 3]);
        // alpha in [lo, lo + 1.2 max(m1, |m2 - m3|)], lo = min(m2, m3). The volume in that range
        // grows at a rate near 1 with alpha, so an error in cut_volume shows up in alpha.
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

    check_planar_fit<PLIC::Youngs>("youngs", false);
    check_planar_fit<PLIC::MixedYoungsCentered>("myc", false);
    check_planar_fit<PLIC::HeightFunctionMYC>("myc_height", true);

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

#include "test_util.H"

#include <AMReX_Print.H>

namespace plic_test {

namespace {
int n_failures = 0;
} // namespace

void check (bool ok, std::string const& what)
{
    Print() << "  " << (ok ? "ok    " : "FAILED") << " " << what << "\n";
    if (!ok) { ++n_failures; }
}

int failures ()
{
    return n_failures;
}

Plane centered_plane (GpuArray<Real, AMREX_SPACEDIM> const& n)
{
    return {n, AMREX_D_TERM(0.5_rt * n[0], + 0.5_rt * n[1], + 0.5_rt * n[2])};
}

Geometry unit_geometry (int n, IntVect const& periodic)
{
    const Box domain(IntVect(0), IntVect(n - 1));
    const RealBox rb({AMREX_D_DECL(0.0_rt, 0.0_rt, 0.0_rt)}, {AMREX_D_DECL(1.0_rt, 1.0_rt, 1.0_rt)});
    return Geometry(domain, rb, 0, {AMREX_D_DECL(periodic[0], periodic[1], periodic[2])});
}

BCRec uniform_bc (int bc)
{
    BCRec rec;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        rec.setLo(dir, bc);
        rec.setHi(dir, bc);
    }
    return rec;
}

void fill_plane (MultiFab& F, Geometry const& geom, Plane const& p, IntVect const& extrapolate)
{
    const auto geomdata = geom.data();
    const IntVect lo = geom.Domain().smallEnd();
    const IntVect hi = geom.Domain().bigEnd();
    auto const& fm = F.arrays();
    ParallelFor(F, F.nGrowVect(), [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept {
        IntVect iv(AMREX_D_DECL(i, j, k));
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            if (extrapolate[dir] != 0) {
                iv[dir] = Clamp(iv[dir], lo[dir], hi[dir]);
            }
        }
        const auto c = iv.dim3();
        fm[b](i, j, k) = PLIC::plane_fraction(c.x, c.y, c.z, geomdata, p.n.data(), p.d);
    });
    Gpu::streamSynchronize();
}

Vector<Box> beyond_faces (Box const& box, Box const& domain, int dir)
{
    Box lo = box;
    Box hi = box;
    lo.setBig(dir, domain.smallEnd(dir) - 1);
    hi.setSmall(dir, domain.bigEnd(dir) + 1);
    return {lo, hi};
}

} // namespace plic_test

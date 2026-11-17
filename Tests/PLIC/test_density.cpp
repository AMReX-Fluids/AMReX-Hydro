// The face density: SharpFaceDensity

#include "test_util.H"

#include <AMReX_Print.H>

#include <random>

namespace plic_test {

namespace {

// Largest curl of rho_face g over the edges away from the domain faces, relative to rho_l g / dx,
// for gravity g: every component in 3D
Real max_curl (Array<MultiFab, AMREX_SPACEDIM>& rho, Geometry const& geom,
               GpuArray<Real, AMREX_SPACEDIM> const& g, Real rho_l)
{
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) { rho[dir].FillBoundary(geom.periodicity()); }
    const Box inner = grow(geom.Domain(), -6);
    const Real h = geom.CellSize(0);
    ReduceOps<ReduceOpMax> op;
    ReduceData<Real> data(op);
    for (int a = 0; a < AMREX_SPACEDIM; ++a) {
        for (int b = a + 1; b < AMREX_SPACEDIM; ++b) {
            const IntVect ea = IntVect::TheDimensionVector(a);
            const IntVect eb = IntVect::TheDimensionVector(b);
            const Real ga = g[a];
            const Real gb = g[b];
            for (MFIter mfi(rho[0]); mfi.isValid(); ++mfi) {
                auto const& ra = rho[a].const_array(mfi);
                auto const& rb = rho[b].const_array(mfi);
                const Box cells = enclosedCells(mfi.validbox());
                const Box nb = surroundingNodes(surroundingNodes(cells & inner, a), b);
                // A reduction over an empty box launches no blocks, which fails on the GPU
                if (!nb.ok()) {
                    continue;
                }
                op.eval(nb, data, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> GpuTuple<Real> {
                    amrex::ignore_unused(k);
                    const IntVect iv(AMREX_D_DECL(i, j, k));
                    // d(rho_a g_a)/db - d(rho_b g_b)/da around the edge
                    const Real c = (ga * (ra(iv) - ra(iv - eb)) - gb * (rb(iv) - rb(iv - ea))) / h;
                    return {nan_as_max(std::abs(c))};
                });
            }
        }
    }
    Real curl = amrex::get<0>(data.value(op));
    ParallelDescriptor::ReduceRealMax(curl);
    return curl * h / rho_l;
}

// The face density balances gravity normal to random planes: every interface cell has the plane,
// so a normal error anywhere shows as a curl
void check_random_planes ()
{
    const Geometry geom = unit_geometry(32, IntVect(0));
    BoxArray ba(geom.Domain());
    ba.maxSize(8);
    DistributionMapping dm(ba);
    const std::string reconstruction = "myc_height";
    const Real rho_l = 1000.0_rt;
    const Real rho_g = 1.0_rt;
    std::mt19937 gen(13);
    std::normal_distribution<Real> normal_dist;
    std::uniform_real_distribution<Real> offset(-0.1_rt, 0.1_rt);
    Real curl = 0.0_rt;
    for (int p = 0; p < 40; ++p) {
        const GpuArray<Real, AMREX_SPACEDIM> nrm{AMREX_D_DECL(normal_dist(gen), normal_dist(gen), normal_dist(gen))};
        Plane pl = centered_plane(nrm);
        pl.d += offset(gen);
        MultiFab F(ba, dm, 1, PLIC::nGrowF(reconstruction));
        fill_plane(F, geom, pl, IntVect(0));
        Array<MultiFab, AMREX_SPACEDIM> rho;
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            rho[dir].define(convert(ba, IntVect::TheDimensionVector(dir)), dm, 1, 1);
            rho[dir].setVal(0.0_rt);
        }
        PLIC::SharpFaceDensity(GetArrOfPtrs(rho), {&F}, geom, rho_l, rho_g, reconstruction);
        curl = max(curl, max_curl(rho, geom, nrm, rho_l));
    }
    Print() << "  largest curl of rho_face g over random planes, relative to rho_l g / dx " << curl << "\n";
    check(curl < 100.0_rt * roundoff_tol, "face density balances gravity normal to random planes");
}

} // namespace

void check_face_density ()
{
    Print() << "face density\n";

    const Geometry geom = unit_geometry(32, IntVect(0));
    const Box domain = geom.Domain();
    BoxArray ba(domain);
    ba.maxSize(8);
    DistributionMapping dm(ba);
    const std::string reconstruction = "myc_height";
    const GpuArray<Real, AMREX_SPACEDIM> nrm = {AMREX_D_DECL(1.0_rt, 0.4_rt, 0.1_rt)};
    const Plane p = centered_plane(nrm);
    MultiFab F(ba, dm, 1, PLIC::nGrowF(reconstruction));
    fill_plane(F, geom, p, IntVect(0));

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

    Real ghost_err = 0.0_rt;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        ghost_err = max(ghost_err, max_over(rho[dir],
            [] (Box const& valid, Box const& fab) { return boxDiff(fab, valid).data(); },
            [=] AMREX_GPU_DEVICE (Array4<Real const> const& a, int i, int j, int k) noexcept -> Real {
                return std::abs(a(i, j, k) - sentinel);
            }));
    }
    check(ghost_err == 0.0_rt, "SharpFaceDensity leaves the ghost faces of rho_face as they are");

    // The face density is a discrete gradient for gravity normal to the plane: its curl vanishes
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
        // A reduction over an empty box launches no blocks, which fails on the GPU
        if (nb.ok()) {
            curl_op.eval(nb, curl_data, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> GpuTuple<Real> {
                // d(rho_x g_x)/dy - d(rho_y g_y)/dx around a node of the xy plane
                const Real c = (nrm[0] * (rx(i, j, k) - rx(i, j - 1, k)) - nrm[1] * (ry(i, j, k) - ry(i - 1, j, k))) / h;
                return {nan_as_max(std::abs(c))};
            });
        }
    }
    Real curl = amrex::get<0>(curl_data.value(curl_op));
    ParallelDescriptor::ReduceRealMax(curl);
    Print() << "  largest curl of rho_face g, relative to rho_l g / dx " << curl * h / rho_l << "\n";
    check(curl * h / rho_l < 100.0_rt * roundoff_tol, "face density balances gravity normal to a plane");

    // Over several F, the segment fractions are averaged before the mixture, which is linear in
    // them: the density of {F, G} is the mean of the densities of F and of G
    MultiFab G(ba, dm, 1, PLIC::nGrowF(reconstruction));
    fill_plane(G, geom, Plane{{AMREX_D_DECL(-0.3_rt, 1.0_rt, 0.5_rt)}, 0.2_rt}, IntVect(0));
    Array<MultiFab, AMREX_SPACEDIM> rho_G;
    Array<MultiFab, AMREX_SPACEDIM> rho_both;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        rho_G[dir].define(rho[dir].boxArray(), dm, 1, 0);
        rho_both[dir].define(rho[dir].boxArray(), dm, 1, 0);
    }
    PLIC::SharpFaceDensity(GetArrOfPtrs(rho_G), {&G}, geom, rho_l, rho_g, reconstruction);
    PLIC::SharpFaceDensity(GetArrOfPtrs(rho_both), {&F, &G}, geom, rho_l, rho_g, reconstruction);
    Real mean_err = 0.0_rt;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        // rho_both - (rho_F + rho_G) / 2, in rho_both
        MultiFab::Saxpy(rho_both[dir], -0.5_rt, rho[dir], 0, 0, 1, 0);
        MultiFab::Saxpy(rho_both[dir], -0.5_rt, rho_G[dir], 0, 0, 1, 0);
        mean_err = max(mean_err, max_abs(rho_both[dir], valid_only));
    }
    Print() << "  largest error of the density over two F, relative to rho_l " << mean_err / rho_l << "\n";
    check(mean_err / rho_l < roundoff_tol, "face density over several F is the mean of their densities");

    check_random_planes();
}

} // namespace plic_test

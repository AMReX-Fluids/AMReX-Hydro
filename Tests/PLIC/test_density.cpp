// The face density: SharpFaceDensity

#include "test_util.H"

#include <AMReX_Print.H>

namespace plic_test {

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
}

} // namespace plic_test

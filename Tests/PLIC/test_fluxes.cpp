// The fluxes of Owkes & Desjardins: the flux volume of a face, and one step of ComputeEdgeStateAndFluxes

#include "test_util.H"

#include <hydro_plic_owkes_desjardins_K.H>

#include <AMReX_Print.H>

#include <sstream>

namespace plic_test {

namespace {

// Node displacements around the low x face of cell 0: x_disp along x for every node, and y_disp0,
// y_disp1 along y for the nodes at y = 0 and y = 1
FArrayBox face_displacements (Real x_disp, Real y_disp0, Real y_disp1)
{
    FArrayBox disp(Box(IntVect(-1), IntVect(2)), 3, The_Pinned_Arena());
    auto const& d = disp.array();
    amrex::LoopOnCpu(disp.box(), [&] (int i, int j, int k) noexcept {
        d(i, j, k, 0) = x_disp;
        d(i, j, k, 1) = (j <= 0) ? y_disp0 : y_disp1;
        d(i, j, k, 2) = 0.0_rt;
    });
    return disp;
}

// Largest distance of a vertex of the flux volume of the low x face of cell 0 beyond the cells
// around the face (x in [-1, 1], the tangential directions in [-1, 2]), and its signed volume
std::pair<Real, Real> flux_volume_extent (FArrayBox const& disp, Real courant)
{
    PLIC::Tet tets[PLIC::max_flux_tets];
    const int ntet = PLIC::flux_volume_tets(0, 0, 0, 0, courant, disp.const_array(), tets);
    Real beyond = 0.0_rt;
    Real vol = 0.0_rt;
    for (int n = 0; n < ntet; ++n) {
        for (auto const& v : tets[n].v) {
            beyond = std::max({beyond, -1.0_rt - v.x[0], v.x[0] - 1.0_rt});
            for (int d = 1; d < AMREX_SPACEDIM; ++d) {
                beyond = std::max({beyond, -1.0_rt - v.x[d], v.x[d] - 2.0_rt});
            }
        }
        vol += PLIC::signed_volume(tets[n].v[0], tets[n].v[1], tets[n].v[2], tets[n].v[3]);
    }
    return {beyond, vol};
}

// The flux volume matches the Courant number under shear, and stays in the cells around the face
// when its traced face has almost no extent normal to x, where the correction apex lies far off
void check_flux_volume ()
{
    {
        const Real courant = 0.28_rt;
        auto [beyond, vol] = flux_volume_extent(face_displacements(-0.3_rt, 0.2_rt, -0.1_rt), courant);
        check(beyond == 0.0_rt && std::abs(vol - courant) <= roundoff_tol,
              "flux volume under shear holds the Courant number, in the cells around the face");
    }
    {
        // The traced nodes at y = 0 and y = 1 meet at y = 1/2, to 1e-9
        const Real h = 0.5_rt - 0.5e-9_rt;
        auto [beyond, vol] = flux_volume_extent(face_displacements(-0.3_rt, h, -h), 0.3_rt);
        Print() << "  degenerate traced face: flux volume reaches " << beyond
                << " cells beyond the cells around the face\n";
        check(beyond == 0.0_rt, "flux volume of a degenerate traced face stays in the cells around the face");
    }
}

// Largest excursion of F from [0,1], away from the domain faces, after one step of a uniform flow of
// Courant numbers c across a plane, with the pure cells next to the interface off by delta
Real excursion_after_step (GpuArray<Real, AMREX_SPACEDIM> const& c, Real delta)
{
    constexpr int n = 16;
    const std::string reconstruction = "myc_height";
    const Geometry geom = unit_geometry(n, IntVect(1));
    const BoxArray ba(geom.Domain());
    const DistributionMapping dm(ba);
    MultiFab F(ba, dm, 1, PLIC::nGrowF(reconstruction), MFInfo().SetArena(The_Pinned_Arena()));
    fill_plane(F, geom, centered_plane({AMREX_D_DECL(1.0_rt, 0.7_rt, 0.4_rt)}), IntVect(0));

    // Off by delta: the pure cells with a mixed neighbor
    for (MFIter mfi(F); mfi.isValid(); ++mfi) {
        auto const& a = F.array(mfi);
        FArrayBox orig(F[mfi], amrex::make_alias, 0, 1);
        FArrayBox copy(orig.box(), 1, The_Pinned_Arena());
        copy.copy<RunOn::Host>(orig);
        auto const& o = copy.const_array();
        amrex::LoopOnCpu(grow(orig.box(), -1), [&] (int i, int j, int k) noexcept {
            if (PLIC::is_mixed(o(i, j, k))) { return; }
            bool near = false;
            amrex::LoopOnCpu(Box(IntVect(-1), IntVect(1)), [&] (int ii, int jj, int kk) noexcept {
                near = near || PLIC::is_mixed(o(i + ii, j + jj, k + kk));
            });
            if (near) { a(i, j, k) = (o(i, j, k) > 0.5_rt) ? 1.0_rt - delta : delta; }
        });
    }

    const Real dt = 1.0_rt;
    Array<MultiFab, AMREX_SPACEDIM> umac, flux, face;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        const BoxArray fba = convert(ba, IntVect::TheDimensionVector(d));
        umac[d].define(fba, dm, 1, 1, MFInfo().SetArena(The_Pinned_Arena()));
        flux[d].define(fba, dm, 1, 0, MFInfo().SetArena(The_Pinned_Arena()));
        face[d].define(fba, dm, 1, 0, MFInfo().SetArena(The_Pinned_Arena()));
        umac[d].setVal(c[d] * geom.CellSize(d) / dt);
    }
    MultiFab Fnew(ba, dm, 1, 0, MFInfo().SetArena(The_Pinned_Arena()));
    for (MFIter mfi(F); mfi.isValid(); ++mfi) {
        PLIC::ComputeEdgeStateAndFluxes(
            mfi.validbox(),
            AMREX_D_DECL(flux[0].array(mfi), flux[1].array(mfi), flux[2].array(mfi)),
            AMREX_D_DECL(face[0].array(mfi), face[1].array(mfi), face[2].array(mfi)),
            F.const_array(mfi),
            AMREX_D_DECL(umac[0].const_array(mfi), umac[1].const_array(mfi), umac[2].const_array(mfi)),
            geom, dt, reconstruction);
    }
    Gpu::streamSynchronize();
    const auto dxinv = geom.InvCellSizeArray();
    Real excursion = 0.0_rt;
    for (MFIter mfi(F); mfi.isValid(); ++mfi) {
        auto const& a = F.const_array(mfi);
        AMREX_D_TERM(auto const& fx = flux[0].const_array(mfi);,
                     auto const& fy = flux[1].const_array(mfi);,
                     auto const& fz = flux[2].const_array(mfi););
        amrex::LoopOnCpu(grow(geom.Domain(), -4), [&] (int i, int j, int k) noexcept {
            const Real f = a(i, j, k) - dt * (AMREX_D_TERM((fx(i + 1, j, k) - fx(i, j, k)) * dxinv[0],
                                                           + (fy(i, j + 1, k) - fy(i, j, k)) * dxinv[1],
                                                           + (fz(i, j, k + 1) - fz(i, j, k)) * dxinv[2]));
            excursion = std::max({excursion, -f, f - 1.0_rt});
        });
    }
    return excursion;
}

// Cells within pure_tol of 0 or 1 count as pure: F then leaves [0,1] by no more than the amount
// they are off by
void check_near_pure_cells ()
{
    const Real delta = 0.5_rt * PLIC::pure_tol;
    for (Real s : {1.0_rt, -1.0_rt}) {
        const GpuArray<Real, AMREX_SPACEDIM> c = {AMREX_D_DECL(0.45_rt * s, 0.35_rt * s, 0.25_rt * s)};
        const Real e0 = excursion_after_step(c, 0.0_rt);
        const Real e1 = excursion_after_step(c, delta);
        Print() << "  one step " << (s > 0 ? "into the gas" : "into the liquid") << ": F leaves [0,1] by "
                << e0 << ", and by " << e1 << " with the pure cells near the interface off by " << delta << "\n";
        std::ostringstream os;
        os << "near-pure cells move F out of [0,1] by no more than they are off by ("
           << (s > 0 ? "into the gas" : "into the liquid") << ")";
        check(e0 <= roundoff_tol && e1 <= delta + roundoff_tol, os.str());
    }
}

} // namespace

void check_fluxes ()
{
    check_flux_volume();
    check_near_pure_cells();
}

} // namespace plic_test

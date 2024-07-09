/**
 * \file hydro_redistribution.cpp
 * \addtogroup Redistribution
 * @{
 *
 */

#include <hydro_redistribution.H>
#include <AMReX_EB_utils.H>

using namespace amrex;

void Redistribution::Apply ( Box const& bx, int ncomp,
                             Array4<Real      > const& dUdt_out,
                             Array4<Real      > const& dUdt_in,
                             Array4<Real const> const& U_in,
                             Array4<Real> const& scratch,
                             Array4<EBCellFlag const> const& flag,
                             AMREX_D_DECL(Array4<Real const> const& apx,
                                          Array4<Real const> const& apy,
                                          Array4<Real const> const& apz),
                             Array4<amrex::Real const> const& vfrac,
                             AMREX_D_DECL(Array4<Real const> const& fcx,
                                          Array4<Real const> const& fcy,
                                          Array4<Real const> const& fcz),
                             Array4<Real const> const& ccc,
                             amrex::BCRec  const* d_bcrec_ptr,
                             Geometry const& lev_geom, Real dt,
                             std::string redistribution_type,
                             amrex::Real target_volfrac)
{
    // redistribution_type = "NoRedist";       // no redistribution
    // redistribution_type = "FluxRedist"      // flux_redistribute
    // redistribution_type = "StateRedist";    // (weighted) state redistribute

    amrex::ParallelFor(bx,ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            dUdt_out(i,j,k,n) = 0.;
        });

    if (redistribution_type == "FluxRedist")
    {
        int icomp = 0;
        apply_flux_redistribution (bx, dUdt_out, dUdt_in, scratch, icomp, ncomp, flag, vfrac, lev_geom);

    } else if (redistribution_type == "StateRedist") {

        Box const& bxg1 = grow(bx,1);
        Box const& bxg2 = grow(bx,2);
        Box const& bxg3 = grow(bx,3);
        Box const& bxg4 = grow(bx,4);

#if (AMREX_SPACEDIM == 2)
        // We assume that in 2D a cell will only need at most 3 neighbors to merge with, and we
        //    use the first component of this for the number of neighbors
        IArrayBox itracker(bxg4,4);
        // How many nbhds is a cell in
#else
        // We assume that in 3D a cell will only need at most 7 neighbors to merge with, and we
        //    use the first component of this for the number of neighbors
        IArrayBox itracker(bxg4,8);
#endif
        FArrayBox nrs_fab(bxg3,1);
        FArrayBox alpha_fab(bxg3,2);

        // Total volume of all cells in my nbhd
        FArrayBox nbhd_vol_fab(bxg2,1);

        // Centroid of my nbhd
        FArrayBox cent_hat_fab     (bxg3,AMREX_SPACEDIM);

        Elixir eli_itr = itracker.elixir();
        Array4<int> itr = itracker.array();
        Array4<int const> itr_const = itracker.const_array();

        Elixir eli_nrs = nrs_fab.elixir();
        Array4<Real      > nrs       = nrs_fab.array();
        Array4<Real const> nrs_const = nrs_fab.const_array();

        Elixir eli_alpha = alpha_fab.elixir();
        Array4<Real      > alpha       = alpha_fab.array();
        Array4<Real const> alpha_const = alpha_fab.const_array();

        Elixir eli_nbf = nbhd_vol_fab.elixir();
        Array4<Real      > nbhd_vol       = nbhd_vol_fab.array();
        Array4<Real const> nbhd_vol_const = nbhd_vol_fab.const_array();

        Elixir eli_chf = cent_hat_fab.elixir();
        Array4<Real      > cent_hat       = cent_hat_fab.array();
        Array4<Real const> cent_hat_const = cent_hat_fab.const_array();

        Box domain_per_grown = lev_geom.Domain();
        AMREX_D_TERM(if (lev_geom.isPeriodic(0)) domain_per_grown.grow(0,1);,
                     if (lev_geom.isPeriodic(1)) domain_per_grown.grow(1,1);,
                     if (lev_geom.isPeriodic(2)) domain_per_grown.grow(2,1););

        // At any external Dirichlet domain boundaries we need to set dUdt_in to 0
        //    in the cells just outside the domain because those values will be used
        //    in the slope computation in state redistribution.  We assume here that
        //    the ext_dir values of U_in itself have already been set.
        if (!domain_per_grown.contains(bxg1))
            amrex::ParallelFor(bxg1,ncomp,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (!domain_per_grown.contains(IntVect(AMREX_D_DECL(i,j,k))))
                        dUdt_in(i,j,k,n) = 0.;
                });

        amrex::ParallelFor(Box(scratch), ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                scratch(i,j,k,n) = U_in(i,j,k,n) + dt * dUdt_in(i,j,k,n);
            }
        );

        MakeITracker(bx, AMREX_D_DECL(apx, apy, apz), vfrac, itr, lev_geom, target_volfrac);

        MakeStateRedistUtils(bx, flag, vfrac, ccc, itr, nrs, alpha, nbhd_vol, cent_hat,
                             lev_geom, target_volfrac);

        StateRedistribute(bx, ncomp, dUdt_out, scratch, flag, vfrac,
                          AMREX_D_DECL(fcx, fcy, fcz), ccc,  d_bcrec_ptr,
                          itr_const, nrs_const, alpha_const, nbhd_vol_const, cent_hat_const, lev_geom);

        amrex::ParallelFor(bx, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                // Only update the values which actually changed -- this makes
                // the results insensitive to tiling -- otherwise cells that aren't
                // changed but are in a tile on which StateRedistribute gets called
                // will have precision-level changes due to adding/subtracting U_in
                // and multiplying/dividing by dt.   Here we test on whether (i,j,k)
                // has at least one neighbor and/or whether (i,j,k) is in the
                // neighborhood of another cell -- if either of those is true the
                // value may have changed

                if (itr(i,j,k,0) > 0 || nrs(i,j,k) > 1.)
                {
                   dUdt_out(i,j,k,n) = (dUdt_out(i,j,k,n) - U_in(i,j,k,n)) / dt;
                }
                else
                {
                   dUdt_out(i,j,k,n) = dUdt_in(i,j,k,n);
                }
            }
        );

    } else if (redistribution_type == "NoRedist") {
        amrex::ParallelFor(bx, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                dUdt_out(i,j,k,n) = dUdt_in(i,j,k,n);
            }
        );

    } else {
       amrex::Error("Not a legit redist_type");
    }
}

void
Redistribution::ApplyToInitialData ( Box const& bx, int ncomp,
                                     Array4<Real      > const& U_out,
                                     Array4<Real      > const& U_in,
                                     Array4<EBCellFlag const> const& flag,
                                     AMREX_D_DECL(amrex::Array4<amrex::Real const> const& apx,
                                                  amrex::Array4<amrex::Real const> const& apy,
                                                  amrex::Array4<amrex::Real const> const& apz),
                                     amrex::Array4<amrex::Real const> const& vfrac,
                                     AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                                                  amrex::Array4<amrex::Real const> const& fcy,
                                                  amrex::Array4<amrex::Real const> const& fcz),
                                     amrex::Array4<amrex::Real const> const& ccc,
                                     amrex::BCRec  const* d_bcrec_ptr,
                                     Geometry& lev_geom, std::string redistribution_type,
                                     amrex::Real target_volfrac)
{
    if (redistribution_type == "StateRedist") {
        amrex::Error("Redistribution::ApplyToInitialData: Shouldn't be here with this redist type");
    }

    Box const& bxg2 = grow(bx,2);
    Box const& bxg3 = grow(bx,3);
    Box const& bxg4 = grow(bx,4);

#if (AMREX_SPACEDIM == 2)
    // We assume that in 2D a cell will only need at most 3 neighbors to merge with, and we
    //    use the first component of this for the number of neighbors
    IArrayBox itracker(bxg4,4);
#else
    // We assume that in 3D a cell will only need at most 7 neighbors to merge with, and we
    //    use the first component of this for the number of neighbors
    IArrayBox itracker(bxg4,8);
#endif
    FArrayBox nrs_fab(bxg3,1);
    FArrayBox alpha_fab(bxg3,2);

    // Total volume of all cells in my nbhd
    FArrayBox nbhd_vol_fab(bxg2,1);

    // Centroid of my nbhd
    FArrayBox cent_hat_fab  (bxg3,AMREX_SPACEDIM);

    Elixir eli_itr = itracker.elixir();
    Array4<int> itr = itracker.array();
    Array4<int const> itr_const = itracker.const_array();

    Elixir eli_nrs = nrs_fab.elixir();
    Array4<Real      > nrs       = nrs_fab.array();
    Array4<Real const> nrs_const = nrs_fab.const_array();

    Elixir eli_alpha = alpha_fab.elixir();
    Array4<Real      > alpha       = alpha_fab.array();
    Array4<Real const> alpha_const = alpha_fab.const_array();

    Elixir eli_nbf = nbhd_vol_fab.elixir();
    Array4<Real      > nbhd_vol       = nbhd_vol_fab.array();
    Array4<Real const> nbhd_vol_const = nbhd_vol_fab.const_array();

    Elixir eli_chf = cent_hat_fab.elixir();
    Array4<Real      > cent_hat       = cent_hat_fab.array();
    Array4<Real const> cent_hat_const = cent_hat_fab.const_array();

    amrex::ParallelFor(bx,ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        U_out(i,j,k,n) = 0.;
    });

    MakeITracker(bx, AMREX_D_DECL(apx, apy, apz), vfrac, itr, lev_geom, target_volfrac);

    MakeStateRedistUtils(bx, flag, vfrac, ccc, itr, nrs, alpha, nbhd_vol, cent_hat,
                            lev_geom, target_volfrac);

    StateRedistribute(bx, ncomp, U_out, U_in, flag, vfrac,
                         AMREX_D_DECL(fcx, fcy, fcz), ccc,  d_bcrec_ptr,
                         itr_const, nrs_const, alpha_const, nbhd_vol_const, cent_hat_const, lev_geom);
}

void
Redistribution::Make1DProfile ( Box const& bx, int ncomp,
                                Array4<Real> const& U_in,
                                Array4<EBCellFlag const> const& flag,
                                Array4<Real const> const& vfrac,
                                AMREX_D_DECL(Array4<Real const> const& fcx,
                                             Array4<Real const> const& fcy,
                                             Array4<Real const> const& fcz),
                                Array4<Real const> const& bcent,
                                Array4<Real const> const& ccent,
                                amrex::BCRec  const* d_bcrec_ptr,
                                Geometry const& lev_geom)
{
    amrex::Vector<amrex::Real> profile;
    amrex::Vector<amrex::Real> xloc;
    profile.resize(64,0.0);
    xloc.resize(64,0.0);

    amrex::Real* prof_ptr = profile.data();
    amrex::Real* xloc_ptr = xloc.data();

    const Box domain = lev_geom.Domain();
    const int domain_ilo = domain.smallEnd(0);
    const int domain_ihi = domain.bigEnd(0);
    const int domain_jlo = domain.smallEnd(1);
    const int domain_jhi = domain.bigEnd(1);
#if (AMREX_SPACEDIM == 3)
    const int domain_klo = domain.smallEnd(2);
    const int domain_khi = domain.bigEnd(2);
#endif

    AMREX_D_TERM(const auto& is_periodic_x = lev_geom.isPeriodic(0);,
                 const auto& is_periodic_y = lev_geom.isPeriodic(1);,
                 const auto& is_periodic_z = lev_geom.isPeriodic(2););

    Box const& bxg1 = amrex::grow(bx,1);
    Box const& bxg2 = amrex::grow(bx,2);
    Box const& bxg3 = amrex::grow(bx,3);

    Box domain_per_grown = domain;
    if (is_periodic_x) domain_per_grown.grow(0,2);
    if (is_periodic_y) domain_per_grown.grow(1,2);
#if (AMREX_SPACEDIM == 3)
    if (is_periodic_z) domain_per_grown.grow(2,2);
#endif

    amrex::ParallelFor(bxg1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (vfrac(i,j,k) > 0.0 && vfrac(i,j,k) < 1.0)
        {
            int max_order = 2;

                for (int n = 0; n < ncomp; n++)
                {
                    bool extdir_ilo = (d_bcrec_ptr[n].lo(0) == amrex::BCType::ext_dir ||
                                       d_bcrec_ptr[n].lo(0) == amrex::BCType::hoextrap);
                    bool extdir_ihi = (d_bcrec_ptr[n].hi(0) == amrex::BCType::ext_dir ||
                                       d_bcrec_ptr[n].hi(0) == amrex::BCType::hoextrap);
                    bool extdir_jlo = (d_bcrec_ptr[n].lo(1) == amrex::BCType::ext_dir ||
                                       d_bcrec_ptr[n].lo(1) == amrex::BCType::hoextrap);
                    bool extdir_jhi = (d_bcrec_ptr[n].hi(1) == amrex::BCType::ext_dir ||
                                       d_bcrec_ptr[n].hi(1) == amrex::BCType::hoextrap);
#if (AMREX_SPACEDIM == 3)
                    bool extdir_klo = (d_bcrec_ptr[n].lo(2) == amrex::BCType::ext_dir ||
                                       d_bcrec_ptr[n].lo(2) == amrex::BCType::hoextrap);
                    bool extdir_khi = (d_bcrec_ptr[n].hi(2) == amrex::BCType::ext_dir ||
                                       d_bcrec_ptr[n].hi(2) == amrex::BCType::hoextrap);
#endif
                    // Initialize so that the slope stencil goes from -1:1 in each direction
                    int nx = 1; int ny = 1; int nz = 1;

                    // Do we have enough extent in each coordinate direction to use the 3x3x3 stencil
                    //    or do we need to enlarge it?
                    AMREX_D_TERM(Real x_max = -1.e30; Real x_min = 1.e30;,
                                 Real y_max = -1.e30; Real y_min = 1.e30;,
                                 Real z_max = -1.e30; Real z_min = 1.e30;);

                    Real slope_stencil_min_width = 0.5;
#if (AMREX_SPACEDIM == 2)
                    int kk = 0;
#elif (AMREX_SPACEDIM == 3)
                    for(int kk(-1); kk<=1; kk++)
#endif
                    {
                     for(int jj(-1); jj<=1; jj++)
                      for(int ii(-1); ii<=1; ii++)
                        if (flag(i,j,k).isConnected(ii,jj,kk))
                        {
                            int r = i+ii; int s = j+jj; int t = k+kk;

                            x_max = std::max(x_max, ccent(r,s,t,0)+static_cast<Real>(ii));
                            x_min = std::min(x_min, ccent(r,s,t,0)+static_cast<Real>(ii));
                            y_max = std::max(y_max, ccent(r,s,t,1)+static_cast<Real>(jj));
                            y_min = std::min(y_min, ccent(r,s,t,1)+static_cast<Real>(jj));
#if (AMREX_SPACEDIM == 3)
                            z_max = std::max(z_max, ccent(r,s,t,2)+static_cast<Real>(kk));
                            z_min = std::min(z_min, ccent(r,s,t,2)+static_cast<Real>(kk));
#endif
                        }
                    }
                    // If we need to grow the stencil, we let it be -nx:nx in the x-direction,
                    //    for example.   Note that nx,ny,nz are either 1 or 2
                    if ( (x_max-x_min) < slope_stencil_min_width ) nx = 2;
                    if ( (y_max-y_min) < slope_stencil_min_width ) ny = 2;
#if (AMREX_SPACEDIM == 3)
                    if ( (z_max-z_min) < slope_stencil_min_width ) nz = 2;
#endif

                    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> slopes_eb;
                    if (nx*ny*nz == 1)
                        // Compute slope using 3x3x3 stencil
                        slopes_eb = amrex_calc_slopes_extdir_eb(
                                                    i,j,k,n,U_in,ccent,vfrac,
                                                    AMREX_D_DECL(fcx,fcy,fcz),flag,
                                                    AMREX_D_DECL(extdir_ilo, extdir_jlo, extdir_klo),
                                                    AMREX_D_DECL(extdir_ihi, extdir_jhi, extdir_khi),
                                                    AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
                                                    AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi),
                                                    max_order);
                    else
                    {
                        // Compute slope using grown stencil (no larger than 5x5x5)
                        slopes_eb = amrex_calc_slopes_extdir_eb_grown(
                                                    i,j,k,n,AMREX_D_DECL(nx,ny,nz),
                                                    U_in,ccent,vfrac,
                                                    AMREX_D_DECL(fcx,fcy,fcz),flag,
                                                    AMREX_D_DECL(extdir_ilo, extdir_jlo, extdir_klo),
                                                    AMREX_D_DECL(extdir_ihi, extdir_jhi, extdir_khi),
                                                    AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
                                                    AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi),
                                                    max_order);
                    }

                    // We do the limiting separately because this limiter limits the slope based on the values
                    //    extrapolated to the cell centroid locations - unlike the limiter in amrex
                    //    which bases the limiting on values extrapolated to the face centroids.
                    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> lim_slope =
                        amrex_calc_centroid_limiter(i,j,k,n,U_in,flag,slopes_eb,ccent);

                    AMREX_D_TERM(lim_slope[0] *= slopes_eb[0];,
                                 lim_slope[1] *= slopes_eb[1];,
                                 lim_slope[2] *= slopes_eb[2];);

                    // Add to the cell itself
                    if (bx.contains(IntVect(AMREX_D_DECL(i,j,k))))
                    {
                        if (n == 0 and j < 80)
                        {
                        prof_ptr[i] =  U_in(i,j,k,n) +
                            lim_slope[0] * (bcent(i,j,k,0)-ccent(i,j,k,0)) +
                            lim_slope[1] * (bcent(i,j,k,1)-ccent(i,j,k,1));
                        xloc_ptr[i] = (i+0.5+bcent(i,j,k,0));
                        }

                    } // if bx contains
                } // n
        } // vfrac
    });

    for (int i = 0; i < 64; i++)
    {
        amrex::Print() << xloc[i] << " " << prof_ptr[i] << std::endl;
    }

}
/** @} */

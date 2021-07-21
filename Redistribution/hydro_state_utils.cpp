/**
 * \file hydro_state_utils.cpp
 * \addtogroup Redistribution
 * @{
 *
 */

#include <hydro_redistribution.H>

using namespace amrex;

void
Redistribution::MakeStateRedistUtils ( Box const& bx,
                                       Array4<EBCellFlag const> const& flag,
                                       Array4<Real const> const& vfrac,
                                       Array4<Real const> const& ccent,
                                       Array4<int  const> const& itracker,
                                       Array4<Real      > const& nrs,
                                       Array4<Real      > const& nbhd_vol,
                                       Array4<Real      > const& cent_hat,
                                       Geometry const& lev_geom)
{
    // Note that itracker has {4 in 2D, 8 in 3D} components and all are initialized to zero
    // We will add to the first component every time this cell is included in a merged neighborhood,
    //    either by merging or being merged
    //
    // In 2D, we identify the cells in the remaining three components with the following ordering
    //
    // ^  6 7 8
    // |  4   5
    // j  1 2 3
    //   i --->
    //
    // In 3D, We identify the cells in the remaining three components with the following ordering
    //
    //    at k-1   |   at k  |   at k+1
    //
    // ^  15 16 17 |  6 7 8  |  24 25 26
    // |  12 13 14 |  4   5  |  21 22 23
    // j  9  10 11 |  1 2 3  |  18 19 20
    //   i --->
    //
    // Note the first component of each of these arrays should never be used
    //
#if (AMREX_SPACEDIM == 2)
    Array<int,9> imap{0,-1, 0, 1,-1, 1,-1, 0, 1};
    Array<int,9> jmap{0,-1,-1,-1, 0, 0, 1, 1, 1};
    Array<int,9> kmap{0, 0, 0, 0, 0, 0, 0, 0, 0};
#else
    Array<int,27>    imap{0,-1, 0, 1,-1, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1};
    Array<int,27>    jmap{0,-1,-1,-1, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1};
    Array<int,27>    kmap{0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
#endif

    AMREX_D_TERM(const auto& is_periodic_x = lev_geom.isPeriodic(0);,
                 const auto& is_periodic_y = lev_geom.isPeriodic(1);,
                 const auto& is_periodic_z = lev_geom.isPeriodic(2););

    Box const& bxg2 = amrex::grow(bx,2);
    Box const& bxg3 = amrex::grow(bx,3);
    Box const& bxg4 = amrex::grow(bx,4);

    const Box domain = lev_geom.Domain();

    Box domain_per_grown = domain;
    if (is_periodic_x) domain_per_grown.grow(0,1);
    if (is_periodic_y) domain_per_grown.grow(1,1);
#if (AMREX_SPACEDIM == 3)
    if (is_periodic_z) domain_per_grown.grow(2,1);
#endif

    amrex::ParallelFor(bxg3,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        // Everyone is in their own neighborhood at least
        nrs(i,j,k) = 1.;
    });

    // nrs captures how many neighborhoods (r,s) is in
    amrex::ParallelFor(bxg4,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        // This loops over the neighbors of (i,j,k), and doesn't include (i,j,k) itself
        for (int i_nbor = 1; i_nbor <= itracker(i,j,k,0); i_nbor++)
        {
            int r = i+imap[itracker(i,j,k,i_nbor)];
            int s = j+jmap[itracker(i,j,k,i_nbor)];
            int t = k+kmap[itracker(i,j,k,i_nbor)];
            if ( domain_per_grown.contains(IntVect(AMREX_D_DECL(r,s,t))) &&
                 bxg3.contains(IntVect(AMREX_D_DECL(r,s,t))) )
            {
                amrex::Gpu::Atomic::Add(&nrs(r,s,t),1.);
            }
        }
    });

    amrex::ParallelFor(bxg2,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (!flag(i,j,k).isCovered())
        {
            Real unwted_vol = 0.; // This is used as a diagnostic to make sure we don't miss any small cells

            // Start with the vfrac of (i,j,k)
            nbhd_vol(i,j,k) = vfrac(i,j,k) / nrs(i,j,k);
            unwted_vol      = vfrac(i,j,k);

            // This loops over the neighbors of (i,j,k), and doesn't include (i,j,k) itself
            for (int i_nbor = 1; i_nbor <= itracker(i,j,k,0); i_nbor++)
            {
                int r = i+imap[itracker(i,j,k,i_nbor)];
                int s = j+jmap[itracker(i,j,k,i_nbor)];
                int t = k+kmap[itracker(i,j,k,i_nbor)];
                nbhd_vol(i,j,k) += vfrac(r,s,t) / nrs(r,s,t);
                unwted_vol      += vfrac(r,s,t);
            }

            // We know stability is guaranteed with unwted_vol > 0.5, but we don't know for sure that it will
            // be unstable for unwted_vol < 0.5.  Here we arbitrarily issue an error if < 0.3 and a warning 
            // if > 0.3 but < 0.5
            if (domain_per_grown.contains(IntVect(AMREX_D_DECL(i,j,k))))
            {
                if (unwted_vol < 0.3) {
                    amrex::Abort("NBHD VOL < 0.3 -- this may be too small");
                } else if (unwted_vol < 0.5) {
                    amrex::Warning("NBHD VOL > 0.3 but < 0.5 -- this may be too small");
                }
            }
        } else {
            nbhd_vol(i,j,k) = 0.;
        }
    });

    // Define xhat,yhat,zhat (from Berger and Guliani)
    amrex::ParallelFor(bxg2,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (vfrac(i,j,k) > 0.5)
        {
            AMREX_D_TERM(cent_hat(i,j,k,0) = ccent(i,j,k,0);,
                         cent_hat(i,j,k,1) = ccent(i,j,k,1);,
                         cent_hat(i,j,k,2) = ccent(i,j,k,2););

        } else if (vfrac(i,j,k) > 0.0 && domain_per_grown.contains(IntVect(AMREX_D_DECL(i,j,k)))) {

            AMREX_D_TERM(cent_hat(i,j,k,0) = ccent(i,j,k,0) * vfrac(i,j,k) / nrs(i,j,k);,
                         cent_hat(i,j,k,1) = ccent(i,j,k,1) * vfrac(i,j,k) / nrs(i,j,k);,
                         cent_hat(i,j,k,2) = ccent(i,j,k,2) * vfrac(i,j,k) / nrs(i,j,k););

            // This loops over the neighbors of (i,j,k), and doesn't include (i,j,k) itself
            for (int i_nbor = 1; i_nbor <= itracker(i,j,k,0); i_nbor++)
            {
                int ii = imap[itracker(i,j,k,i_nbor)]; int r = i+ii;
                int jj = jmap[itracker(i,j,k,i_nbor)]; int s = j+jj;
                int kk = kmap[itracker(i,j,k,i_nbor)]; int t = k+kk;

                AMREX_D_TERM(cent_hat(i,j,k,0) += (ccent(r,s,t,0) + ii) * vfrac(r,s,t) / nrs(r,s,t);,
                             cent_hat(i,j,k,1) += (ccent(r,s,t,1) + jj) * vfrac(r,s,t) / nrs(r,s,t);,
                             cent_hat(i,j,k,2) += (ccent(r,s,t,2) + kk) * vfrac(r,s,t) / nrs(r,s,t););
            }

            AMREX_D_TERM(cent_hat(i,j,k,0) /= nbhd_vol(i,j,k);,
                         cent_hat(i,j,k,1) /= nbhd_vol(i,j,k);,
                         cent_hat(i,j,k,2) /= nbhd_vol(i,j,k););
        } else {

            AMREX_D_TERM(cent_hat(i,j,k,0) = 0.;,
                         cent_hat(i,j,k,1) = 0.;,
                         cent_hat(i,j,k,2) = 0.;);
        }
    });
}
/** @} */

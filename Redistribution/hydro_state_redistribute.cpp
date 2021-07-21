/**
 * \file hydro_state_redistribute.cpp
 * \addtogroup Redistribution
 * @{
 *
 */

#include <hydro_redistribution.H>
#include <AMReX_EB_slopes_K.H>

using namespace amrex;

void
Redistribution::StateRedistribute ( Box const& bx, int ncomp,
                                    Array4<Real> const& U_out,
                                    Array4<Real> const& U_in,
                                    Array4<EBCellFlag const> const& flag,
                                    Array4<Real const> const& vfrac,
                                    AMREX_D_DECL(Array4<Real const> const& fcx,
                                                 Array4<Real const> const& fcy,
                                                 Array4<Real const> const& fcz),
                                    Array4<Real const> const& ccent,
                                    amrex::BCRec  const* d_bcrec_ptr,
                                    Array4< int const> const& itracker,
                                    Array4<Real const> const& nrs,
                                    Array4<Real const> const& nbhd_vol,
                                    Array4<Real const> const& cent_hat,
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

    // amrex::Print() << " IN STATE_REDISTRIBUTE DOING BOX " << bx << " with ncomp " << ncomp << std::endl;
    // amrex::Print() << " Box(U_in) " << Box(U_in) << std::endl;
    // amrex::Print() << " Box(U_out) " << Box(U_out) << std::endl;

    Box const& bxg1 = amrex::grow(bx,1);
    Box const& bxg2 = amrex::grow(bx,2);

    Box domain_per_grown = domain;
    if (is_periodic_x) domain_per_grown.grow(0,1);
    if (is_periodic_y) domain_per_grown.grow(1,1);
#if (AMREX_SPACEDIM == 3)
    if (is_periodic_z) domain_per_grown.grow(2,1);
#endif

    // Solution at the centroid of my nbhd
    FArrayBox    soln_hat_fab (bxg2,ncomp);
    Array4<Real> soln_hat = soln_hat_fab.array();
    Elixir   eli_soln_hat = soln_hat_fab.elixir();

    // Define Qhat (from Berger and Guliani)
    amrex::ParallelFor(bxg2, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        soln_hat(i,j,k,n) = 1.e40; // This is just here for diagnostic purposes

        // This is needed to retain Dirichlet values on domain faces for use in the slope routine
        if (vfrac(i,j,k) > 0.0 && !domain_per_grown.contains(IntVect(AMREX_D_DECL(i,j,k))))
            soln_hat(i,j,k,n) = U_in(i,j,k,n);

        if (vfrac(i,j,k) > 0.5)
        {
            soln_hat(i,j,k,n) = U_in(i,j,k,n);

        } else if (vfrac(i,j,k) > 0.0 && domain_per_grown.contains(IntVect(AMREX_D_DECL(i,j,k)))) {

            // Start with U_in(i,j,k) itself
            soln_hat(i,j,k,n) = U_in(i,j,k,n) * vfrac(i,j,k) / nrs(i,j,k);

            // This loops over the neighbors of (i,j,k), and doesn't include (i,j,k) itself
            for (int i_nbor = 1; i_nbor <= itracker(i,j,k,0); i_nbor++)
            {
                int r = i+imap[itracker(i,j,k,i_nbor)];
                int s = j+jmap[itracker(i,j,k,i_nbor)];
                int t = k+kmap[itracker(i,j,k,i_nbor)];

                if (domain_per_grown.contains(IntVect(AMREX_D_DECL(r,s,t))))
                {
                    soln_hat(i,j,k,n) += U_in(r,s,t,n) * vfrac(r,s,t) / nrs(r,s,t);
                }
            }
            soln_hat(i,j,k,n) /= nbhd_vol(i,j,k);
        }
    });

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (vfrac(i,j,k) > 0.0)
        {
            int max_order = 2;
            for (int n = 0; n < ncomp; n++)
            {
                U_out(i,j,k,n) +=soln_hat(i,j,k,n);

                // Recall that above if (vfrac > 0.5) we set cent_hat equal to ccent ... which means
                // even if we compute the slopes we are multiplying them by zero ... so no point computing them
                if (vfrac(i,j,k) <= 0.5)
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

                    const auto& slopes_eb = amrex_lim_slopes_extdir_eb(i,j,k,n,soln_hat,cent_hat,vfrac,
                                                                AMREX_D_DECL(fcx,fcy,fcz),flag,
                                                                AMREX_D_DECL(extdir_ilo, extdir_jlo, extdir_klo), 
                                                                AMREX_D_DECL(extdir_ihi, extdir_jhi, extdir_khi), 
                                                                AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo), 
                                                                AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi), 
                                                                max_order);

                    AMREX_D_TERM(U_out(i,j,k,n) += slopes_eb[0] * (ccent(i,j,k,0)-cent_hat(i,j,k,0));,
                                 U_out(i,j,k,n) += slopes_eb[1] * (ccent(i,j,k,1)-cent_hat(i,j,k,1));,
                                 U_out(i,j,k,n) += slopes_eb[2] * (ccent(i,j,k,2)-cent_hat(i,j,k,2)););
                }
            } // n
        } // vfrac
    });

    // Note -- we can't merge this loop with the above loop because that was over bx and this is 
    //         over bxg1
    amrex::ParallelFor(bxg1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (vfrac(i,j,k) > 0.0)
        {
            int max_order = 2;
            int num_nbors = itracker(i,j,k,0);

            if (num_nbors > 0)
            {
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

                    const auto& slopes_eb = amrex_lim_slopes_extdir_eb(i,j,k,n,soln_hat,cent_hat,vfrac,
                                                                       AMREX_D_DECL(fcx,fcy,fcz),flag,
                                                                       AMREX_D_DECL(extdir_ilo, extdir_jlo, extdir_klo), 
                                                                       AMREX_D_DECL(extdir_ihi, extdir_jhi, extdir_khi), 
                                                                       AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo), 
                                                                       AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi), 
                                                                       max_order);

                    // This loops over the neighbors of (i,j,k), and doesn't include (i,j,k) itself
                    for (int i_nbor = 1; i_nbor <= num_nbors; i_nbor++)
                    {
                        int r = i+imap[itracker(i,j,k,i_nbor)];
                        int s = j+jmap[itracker(i,j,k,i_nbor)];
                        int t = k+kmap[itracker(i,j,k,i_nbor)];

                        if (bx.contains(IntVect(AMREX_D_DECL(r,s,t))))
                        {
                            Real update = soln_hat(i,j,k,n);
                            AMREX_D_TERM(update += slopes_eb[0] * (ccent(r,s,t,0)-cent_hat(i,j,k,0));,
                                         update += slopes_eb[1] * (ccent(r,s,t,1)-cent_hat(i,j,k,1));,
                                         update += slopes_eb[2] * (ccent(r,s,t,2)-cent_hat(i,j,k,2)););
		    	    amrex::Gpu::Atomic::Add(&U_out(r,s,t,n),update);
                        } // if bx contains
                    } // i_nbor
                } // n
            } // num_nbors
        } // vfrac
    });

    amrex::ParallelFor(bx,ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (!flag(i,j,k).isCovered())
        {
            // This seems to help with a compiler issue ...
            Real denom = 1. / (nrs(i,j,k) + 1.e-40);
            U_out(i,j,k,n) *= denom;
        }
        else
        {
            U_out(i,j,k,n) = 1.e40;
        }
    });

#if 0
    //
    // This tests whether the redistribution procedure was conservative -- this is only relevant
    //      if bx is the whole domain
    //
    { // STRT:SUM OF FINAL DUDT
      for (int n = 0; n < ncomp; n++)
      {
        Real sum1(0);
        Real sum2(0);
#if (AMREX_SPACEDIM == 2)
        int k = 0;
#else
        for (int k = bx.smallEnd(2); k <= domain.bigEnd(2); k++)
#endif
        for (int j = bx.smallEnd(1); j <= domain.bigEnd(1); j++)
        for (int i = bx.smallEnd(0); i <= domain.bigEnd(0); i++)
        {
            sum1 += vfrac(i,j,k)*U_in(i,j,k,n);
            sum2 += vfrac(i,j,k)*U_out(i,j,k,n);
        }
        if (std::abs(sum1-sum2) > 1.e-8 * sum1 && std::abs(sum1-sum2) > 1.e-8)
        {
           printf("SUMS DO NOT MATCH IN STATE REDIST : %f %f ",sum1,sum2);
           amrex::Abort();
        }
      }
    } //  END:SUM OF FINAL DUDT
#endif
}
/** @} */

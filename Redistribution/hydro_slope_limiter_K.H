/**
 * \file hydro_slope_limiter_K.H
 * \addtogroup Redistribution
 * @{
 *
 */

using namespace amrex;

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real
amrex_calc_alpha_stencil(Real q_hat, Real q_max, Real q_min, Real state) noexcept
{
    Real small = 1.0e-12*std::abs(q_max);
    Real alpha;

    if ((q_hat-state) > small) {
        alpha = std::min(1.0,(q_max-state)/(q_hat-state));
    } else if ((q_hat-state) < -small) {
        alpha = std::min(1.0,(q_min-state)/(q_hat-state));
    } else {
        alpha = 1.0;
    }
    return alpha;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>
amrex_calc_centroid_limiter(int i, int j, int k, int n,
                         amrex::Array4<amrex::Real const> const& state,
                         amrex::Array4<amrex::EBCellFlag const> const& flag,
                         const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>& slopes,
                         amrex::Array4<amrex::Real const> const& ccent) noexcept
{
    AMREX_D_TERM(int cuts_x = 0;,
                 int cuts_y = 0;,
                 int cuts_z = 0;);

#if (AMREX_SPACEDIM == 2)
    int kk = 0;
#endif

    Real q_min = state(i,j,k,n);
    Real q_max = state(i,j,k,n);

    // 1) Compute how many cut or regular faces we have in 3x3(x3) block
    // 2) Set the max and min of the values at the centroids in this block
#if (AMREX_SPACEDIM == 3)
    for(int kk(-1); kk<=1; kk++)
#endif
    {
        for(int jj(-1); jj<=1; jj++){
          for(int ii(-1); ii<=1; ii++){
            if( flag(i,j,k).isConnected(ii,jj,kk) and
                not (ii==0 && jj==0 && kk==0)) {
                if ((ii==-1 || ii==1) && jj==0 && kk==0) cuts_x++;
                if ((jj==-1 || jj==1) && ii==0 && kk==0) cuts_y++;
#if (AMREX_SPACEDIM == 3)
                if ((kk==-1 || kk==1) && ii==0 && jj==0) cuts_z++;
#endif

                q_max = std::max(state(i+ii,j+jj,k+kk,n),q_max);
                q_min = std::min(state(i+ii,j+jj,k+kk,n),q_min);
            }
          }
        }
    }

    amrex::Real alpha = 1.0;

    // Compute the limiter needed to keep the predicted q_hat between the max and min
#if (AMREX_SPACEDIM == 3)
    for(int kk(-1); kk<=1; kk++)
#endif
    {
     for(int jj(-1); jj<=1; jj++){
      for(int ii(-1); ii<=1; ii++){
        if (flag(i,j,k).isConnected(ii,jj,kk))
        {
            Real delta_x = ccent(i+ii,j+jj,k+kk,0) -  ccent(i,j,k,0) + static_cast<Real>(ii);
            Real delta_y = ccent(i+ii,j+jj,k+kk,1) -  ccent(i,j,k,1) + static_cast<Real>(jj);
            Real delta_z = ccent(i+ii,j+jj,k+kk,2) -  ccent(i,j,k,2) + static_cast<Real>(kk);

            Real q_hat = state(i,j,k,n) + delta_x * slopes[0] + delta_y * slopes[1] + delta_z * slopes[2];

            if (q_hat > (q_max+1e-12) || q_hat < (q_min-1.e-12))
            {
                alpha = std::min(alpha,amrex_calc_alpha_stencil(q_hat, q_max, q_min, state(i,j,k,n)));
            }
        }
      }
     }
    }

    amrex::Real xalpha = alpha;
    amrex::Real yalpha = alpha;
    amrex::Real zalpha = alpha;

    return {AMREX_D_DECL(xalpha,yalpha,zalpha)};
}
/** @} */

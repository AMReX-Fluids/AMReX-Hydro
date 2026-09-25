/**
 * \file hydro_ebgodunov_transverse_2D_K.H
 * \addtogroup EBGodunov
 * @{
 *
 */

#ifndef HYDRO_EBGODUNOV_TRANSVERSE_2D_K_H_
#define HYDRO_EBGODUNOV_TRANSVERSE_2D_K_H_

#include <AMReX_Array4.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

namespace EBGodunovTransverse {

AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real
linear_interp (amrex::Real xout,
               amrex::Real x1, amrex::Real val1,
               amrex::Real x2, amrex::Real val2)
{
   return (val1 + (xout-x1)/(x2-x1) * (val2-val1));
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void create_transverse_terms_for_xface(int ic, int j, int k,
                                       amrex::Array4<amrex::Real const> const& v_ad,
                                       amrex::Array4<amrex::Real const> const& yhat,
                                       amrex::Array4<amrex::Real const> const& apy,
                                       amrex::Array4<amrex::Real const> const& fcy,
                                       amrex::Real& trans_y,
                                       const amrex::Real dy)
{
    using namespace amrex;

    Real v_tmp_j, v_tmp_jp1;
    Real y_hat_j, y_hat_jp1;

    // Tangential extrapolation in the x-direction at j-1/2
    if (apy(ic,j,k) == amrex::Real(1.0))
    {
        v_tmp_j = v_ad(ic,j,k);
        y_hat_j = yhat(ic,j,k);
    }
    else
    {
        // Tangential extrapolation in x-direction at j-1/2
        int jc    = j;
        int ii    = (fcy(ic,jc,k,0) <= amrex::Real(0.0)) ? ic-1 : ic+1;
        Real xoff = (fcy(ic,jc,k,0) <= amrex::Real(0.0)) ?  -amrex::Real(1.) : amrex::Real(1.);

        v_tmp_j = linear_interp(amrex::Real(0.0), fcy(ic,jc,k), v_ad(ic,jc,k), (xoff+fcy(ii,jc,k)), v_ad(ii,jc,k));
        y_hat_j = linear_interp(amrex::Real(0.0), fcy(ic,jc,k), yhat(ic,jc,k), (xoff+fcy(ii,jc,k)), yhat(ii,jc,k));
    }

    // Tangential extrapolation in the x-direction at j+1/2
    if (apy(ic,j+1,k) == amrex::Real(1.0))
    {
        v_tmp_jp1 = v_ad(ic,j+1,k);
        y_hat_jp1 = yhat(ic,j+1,k);
    }
    else
    {
        // Tangential extrapolation in x-direction at j+1/2
        int jc    = j+1;
        int ii    = (fcy(ic,jc,k,0) <= amrex::Real(0.0)) ? ic-1 : ic+1;
        Real xoff = (fcy(ic,jc,k,0) <= amrex::Real(0.0)) ?  -amrex::Real(1.) : amrex::Real(1.);

        v_tmp_jp1 = linear_interp(amrex::Real(0.0), fcy(ic,jc,k), v_ad(ic,jc,k), (xoff+fcy(ii,jc,k)), v_ad(ii,jc,k));
        y_hat_jp1 = linear_interp(amrex::Real(0.0), fcy(ic,jc,k), yhat(ic,jc,k), (xoff+fcy(ii,jc,k)), yhat(ii,jc,k));
    }

    trans_y = amrex::Real(0.5) * (v_tmp_jp1 + v_tmp_j) * (y_hat_jp1 - y_hat_j) / dy;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void create_transverse_terms_for_yface(int i, int jc, int k,
                                       amrex::Array4<amrex::Real const> const& u_ad,
                                       amrex::Array4<amrex::Real const> const& xhat,
                                       amrex::Array4<amrex::Real const> const& apx,
                                       amrex::Array4<amrex::Real const> const& fcx,
                                       amrex::Real& trans_x,
                                       const amrex::Real dx)
{
    using namespace amrex;

    Real u_tmp_i, u_tmp_ip1;
    Real x_hat_i, x_hat_ip1;

    // Tangential extrapolation in the y-direction at i-1/2
    if (apx(i,jc,k) == amrex::Real(1.0))
    {
        u_tmp_i = u_ad(i,jc,k);
        x_hat_i = xhat(i,jc,k);
    }
    else
    {
        // Tangential extrapolation in y-direction at i-1/2
        int ic    = i;
        int jj    = (fcx(ic,jc,k,0) <= amrex::Real(0.0)) ? jc-1 : jc+1;
        Real yoff = (fcx(ic,jc,k,0) <= amrex::Real(0.0)) ?  -amrex::Real(1.) : amrex::Real(1.);

        u_tmp_i = linear_interp(amrex::Real(0.0), fcx(ic,jc,k), u_ad(ic,jc,k), (yoff+fcx(ic,jj,k)), u_ad(ic,jj,k));
        x_hat_i = linear_interp(amrex::Real(0.0), fcx(ic,jc,k), xhat(ic,jc,k), (yoff+fcx(ic,jj,k)), xhat(ic,jj,k));
    }

    // Tangential extrapolation in the x-direction at i+1/2
    if (apx(i+1,jc,k) == amrex::Real(1.0))
    {
        u_tmp_ip1 = u_ad(i+1,jc,k);
        x_hat_ip1 = xhat(i+1,jc,k);
    }
    else
    {
        // Tangential extrapolation in y-direction at i+1/2
        int ic = i+1;
        int jj = (fcx(ic,jc,k,0) <= amrex::Real(0.0)) ? jc-1 : jc+1;
        Real yoff = (fcx(ic,jc,k,0) <= amrex::Real(0.0)) ?  -amrex::Real(1.) : amrex::Real(1.);

        u_tmp_ip1 = linear_interp(amrex::Real(0.0), fcx(ic,jc,k), u_ad(ic,jc,k), (yoff+fcx(ic,jj,k)), u_ad(ic,jj,k));
        x_hat_ip1 = linear_interp(amrex::Real(0.0), fcx(ic,jc,k), xhat(ic,jc,k), (yoff+fcx(ic,jj,k)), xhat(ic,jj,k));
    }

    trans_x = amrex::Real(0.5) * (u_tmp_ip1 + u_tmp_i) * (x_hat_ip1 - x_hat_i) / dx;
}

} // namespace EBGodunovTransverse

#endif
/** @} */

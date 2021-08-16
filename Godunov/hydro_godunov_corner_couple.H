/**
 * \file hydro_godunov_corner_couple.H
 *
 * \addtogroup Godunov
 *  @{
 */

#ifndef HYDRO_GODUNOV_CORNER_COUPLE_H
#define HYDRO_GODUNOV_CORNER_COUPLE_H

#include <AMReX_Gpu.H>
#include <AMReX_BCRec.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Array.H>
#include <iomanip>
#include <hydro_constants.H>

namespace GodunovCornerCouple {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void AddCornerCoupleTermYX ( amrex::Real& lo1, amrex::Real& hi1,
                             int i, int j, int k, int n, amrex::Real dt, amrex::Real dx,
                             bool iconserv, amrex::Real lo, amrex::Real hi,
                             amrex::Array4<amrex::Real const> const& s,
                             amrex::Array4<amrex::Real const> const& divu_cc,
                             amrex::Array4<amrex::Real const> const& mac,
                             amrex::Array4<amrex::Real const> const& state )
{
    // Modify state on y-faces with x-derivatives to be used for computing state on z-faces

    lo1 = lo;
    hi1 = hi;

    // Here we add  - dt/3 ( (q u)_x - q u_x ) = -dt/3 (u q_x)
    lo1 += - dt/(3.*dx)*( state(i+1,j-1,k,n)*mac(i+1,j-1,k)
                         -state(i  ,j-1,k,n)*mac(i  ,j-1,k) )
           + dt/(3.*dx)*s(i,j-1,k,n)*(mac(i+1,j-1,k)-mac(i,j-1,k));
        hi1 += - dt/(3.*dx)*( state(i+1,j  ,k,n)*mac(i+1,j  ,k)
                             -state(i  ,j  ,k,n)*mac(i  ,j  ,k) )
           + dt/(3.*dx)*s(i,j  ,k,n)*(mac(i+1,j  ,k) - mac(i,j  ,k));

    // Here we add  - dt/3 ( q divu ) if conservative
    lo1 += (iconserv) ? - dt/(3.) * s(i,j-1,k,n)*divu_cc(i,j-1,k) : 0.;
    hi1 += (iconserv) ? - dt/(3.) * s(i,j  ,k,n)*divu_cc(i,j  ,k) : 0.;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void AddCornerCoupleTermZX ( amrex::Real& lo1, amrex::Real& hi1,
                             int i, int j, int k, int n, amrex::Real dt, amrex::Real dx,
                             bool iconserv, amrex::Real lo, amrex::Real hi,
                             amrex::Array4<amrex::Real const> const& s,
                             amrex::Array4<amrex::Real const> const& divu_cc,
                             amrex::Array4<amrex::Real const> const& mac,
                             amrex::Array4<amrex::Real const> const& state )
{
    // Modify state on z-faces with x-derivatives to be used for computing state on y-faces

    lo1 = lo;
    hi1 = hi;

    // Here we add  - dt/3 ( (q u)_x - q u_x ) = -dt/3 (u q_x)
    lo1 += - dt/(3.*dx)*(state(i+1,j,k-1,n)*mac(i+1,j,k-1)
                       - state(i  ,j,k-1,n)*mac(i  ,j,k-1))
           + dt/(3.*dx)*s(i,j,k-1,n)*(mac(i+1,j,k-1)-mac(i,j,k-1));

    hi1 += - dt/(3.*dx)*( state(i+1,j,k  ,n)*mac(i+1,j,k  )
                         -state(i  ,j,k  ,n)*mac(i  ,j,k  ) )
           + dt/(3.*dx)*s(i,j,k  ,n)*(mac(i+1,j,k  ) - mac(i,j,k  ));

    // Here we add  - dt/3 ( q divu ) if conservative
    lo1 += (iconserv) ? - dt/(3.) * s(i,j,k-1,n)*divu_cc(i,j,k-1) : 0.;
    hi1 += (iconserv) ? - dt/(3.) * s(i,j,k  ,n)*divu_cc(i,j,k  ) : 0.;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void AddCornerCoupleTermXY ( amrex::Real& lo1, amrex::Real& hi1,
                             int i, int j, int k, int n, amrex::Real dt, amrex::Real dy,
                             bool iconserv, amrex::Real lo, amrex::Real hi,
                             amrex::Array4<amrex::Real const> const& s,
                             amrex::Array4<amrex::Real const> const& divu_cc,
                             amrex::Array4<amrex::Real const> const& mac,
                             amrex::Array4<amrex::Real const> const& state )
{
    // Modify state on x-faces with y-derivatives to be used for computing state on z-faces

    lo1 = lo;
    hi1 = hi;

    // Here we add  - dt/3 ( (q v)_y - q v_y ) = -dt/3 (v q_y)
    lo1 += - dt/(3.*dy)*( state(i-1,j+1,k,n)*mac(i-1,j+1,k)
                           -state(i-1,j  ,k,n)*mac(i-1,j  ,k) )
           + dt/(3.*dy)*s(i-1,j,k,n)*(mac(i-1,j+1,k)-mac(i-1,j,k));
    hi1 += - dt/(3.*dy)*( state(i  ,j+1,k,n)*mac(i  ,j+1,k)
                           -state(i  ,j  ,k,n)*mac(i  ,j  ,k) )
           + dt/(3.*dy)*s(i,j,k,n)*(mac(i,j+1,k)-mac(i,j,k));

    // Here we add  - dt/3 ( q divu ) if conservative
    lo1 += (iconserv) ? - dt/(3.) * s(i-1,j,k,n)*divu_cc(i-1,j,k) : 0.;
    hi1 += (iconserv) ? - dt/(3.) * s(i  ,j,k,n)*divu_cc(i  ,j,k) : 0.;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void AddCornerCoupleTermZY ( amrex::Real& lo1, amrex::Real& hi1,
                             int i, int j, int k, int n, amrex::Real dt, amrex::Real dy,
                             bool iconserv, amrex::Real lo, amrex::Real hi,
                             amrex::Array4<amrex::Real const> const& s,
                             amrex::Array4<amrex::Real const> const& divu_cc,
                             amrex::Array4<amrex::Real const> const& mac,
                             amrex::Array4<amrex::Real const> const& state )
{
    // Modify state on z-faces with y-derivatives to be used for computing state on x-faces

    lo1 = lo;
    hi1 = hi;

    // Here we add  - dt/3 ( (q v)_y - q v_y ) = -dt/3 (v q_y)
    lo1 += - dt/(3.*dy)*( state(i,j+1,k-1,n)*mac(i,j+1,k-1)
                         -state(i,j  ,k-1,n)*mac(i,j  ,k-1) )
           + dt/(3.*dy)*s(i,j,k-1,n)*(mac(i,j+1,k-1)-mac(i,j,k-1));
    hi1 += - dt/(3.*dy)*( state(i,j+1,k  ,n)*mac(i,j+1,k  )
                         -state(i,j  ,k  ,n)*mac(i,j  ,k  ) )
           + dt/(3.*dy)*s(i,j,k,n)*(mac(i,j+1,k) - mac(i,j,k));

    // Here we add  - dt/3 ( q divu ) if conservative
    lo1 += (iconserv) ? - dt/(3.) * s(i,j,k-1,n)*divu_cc(i,j,k-1) : 0.;
    hi1 += (iconserv) ? - dt/(3.) * s(i,j,k  ,n)*divu_cc(i,j,k  ) : 0.;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void AddCornerCoupleTermXZ ( amrex::Real& lo1, amrex::Real& hi1,
                             int i, int j, int k, int n, amrex::Real dt, amrex::Real dz,
                             bool iconserv, amrex::Real lo, amrex::Real hi,
                             amrex::Array4<amrex::Real const> const& s,
                             amrex::Array4<amrex::Real const> const& divu_cc,
                             amrex::Array4<amrex::Real const> const& mac,
                             amrex::Array4<amrex::Real const> const& state)
{
    // Modify state on x-faces with z-derivatives to be used for computing state on y-faces

    lo1 = lo;
    hi1 = hi;

    // Here we add  - dt/3 ( (q w)_z - q w_z ) = -dt/3 (w q_z)
    lo1 += - dt/(3.*dz)*( state(i-1,j,k+1,n)*mac(i-1,j,k+1)
                         -state(i-1,j,k  ,n)*mac(i-1,j,k) )
           + dt/(3.*dz)*s(i-1,j,k,n)*(mac(i-1,j,k+1) - mac(i-1,j,k));
    hi1 += - dt/(3.*dz)*( state(i  ,j,k+1,n)*mac(i  ,j,k+1)
                         -state(i  ,j,k  ,n)*mac(i  ,j,k) )
           + dt/(3.*dz)*s(i,j,k,n)*(mac(i,j,k+1) - mac(i,j,k));

    // Here we add  - dt/3 ( q divu ) if conservative
    lo1 += (iconserv) ? - dt/(3.) * s(i-1,j,k,n)*divu_cc(i-1,j,k) : 0.;
    hi1 += (iconserv) ? - dt/(3.) * s(i  ,j,k,n)*divu_cc(i  ,j,k) : 0.;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void AddCornerCoupleTermYZ ( amrex::Real& lo1, amrex::Real& hi1,
                             int i, int j, int k, int n, amrex::Real dt, amrex::Real dz,
                             bool iconserv, amrex::Real lo, amrex::Real hi,
                             amrex::Array4<amrex::Real const> const& s,
                             amrex::Array4<amrex::Real const> const& divu_cc,
                             amrex::Array4<amrex::Real const> const& mac,
                             amrex::Array4<amrex::Real const> const& state)
{
    // Modify state on y-faces with z-derivatives to be used for computing state on x-faces

    lo1 = lo;
    hi1 = hi;

    // Here we add  - dt/3 ( (q w)_z - q w_z ) = -dt/3 (w q_z)
    lo1 += - dt/(3.*dz)*( state(i,j-1,k+1,n)*mac(i,j-1,k+1)
                         -state(i,j-1,k  ,n)*mac(i,j-1,k  ) )
           + dt/(3.*dz)*s(i,j-1,k,n)*(mac(i,j-1,k+1)-mac(i,j-1,k));
    hi1 += - dt/(3.*dz)*( state(i,j  ,k+1,n)*mac(i,j  ,k+1)
                         -state(i,j  ,k  ,n)*mac(i,j  ,k  ) )
           + dt/(3.*dz)*s(i,j,k,n)*(mac(i,j,k+1) - mac(i,j,k));

    // Here we add  - dt/3 ( q divu ) if conservative
    lo1 += (iconserv) ? - dt/(3.) * s(i,j-1,k,n)*divu_cc(i,j-1,k) : 0.;
    hi1 += (iconserv) ? - dt/(3.) * s(i,j  ,k,n)*divu_cc(i,j  ,k) : 0.;
}
}
#endif
/** @} */

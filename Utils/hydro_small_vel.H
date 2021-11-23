/**
 * \file hydro_small_vel.H
 */


/** \addtogroup Utilities
 * @{
 */

#ifndef HYDRO_SMALL_VEL_H_
#define HYDRO_SMALL_VEL_H_

#include <AMReX_REAL.H>
#include <AMReX_Math.H>
#include <algorithm>


namespace HydroSmallVel {
    
/**
 * \var small_vel_eps
 *
 * Used in determining what's considered an ~zero velocity when upwinding, either
 * as small_vel_eps*vel_scale for computing edge states or in calc_small_vel for
 * computing the advective velocity.
 *
 */
static constexpr amrex::Real small_vel_eps = 1.e-8;
/**
 * \var vel_scale
 *
 * Order of magnitude of the expected velocity field. 
 * Used in determining what's considered an ~zero velocity when upwinding, either
 * as small_vel_eps*vel_scale for computing edge states or in calc_small_vel for
 * computing the advective velocity.
 * 
 */
static constexpr amrex::Real vel_scale = 1.;
static constexpr amrex::Real small_vel = small_vel_eps*vel_scale;

/**
 * \brief Determines what's considered an ~zero velocity when upwinding during
 * the compuation of the advective velocity (Umac).
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real calc_small_vel (amrex::Real lo, amrex::Real hi)
{
    // This routine is for comparing two predicted edge states (predicted
    // from the low side and high side), so enlarge epsilon by an order of
    // magnitude
    amrex::Real eps = 10*small_vel_eps;
    amrex::Real small_v;
    amrex::Real vmax = std::max(amrex::Math::abs(lo),amrex::Math::abs(hi));

    small_v = vmax < eps*vel_scale ?
		     eps*vel_scale : vmax*eps;
    return small_v;
}

}

#endif
/** @}*/

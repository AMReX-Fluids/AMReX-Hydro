/**
 * \file hydro_constants.H
 */


/** \addtogroup Utilities
 * @{
 */

#ifndef HYDRO_CONSTANTS_H_
#define HYDRO_CONSTANTS_H_

#include <AMReX_REAL.H>
#include <AMReX_Math.H>
#include <algorithm>

/**
 * \def XVEL
 * Index for the X velocity component.
 */
/**
 * \def YVEL
 * Index for the Y velocity component.
 */
/**
 * \def ZVEL
 * Index for the Z velocity component.
 */
#define XVEL 0
#define YVEL 1
#define ZVEL 2

/**
 * \var small_vel
 *
 * What we consider a ~zero velocity when upwinding.
 */
static constexpr amrex::Real small_vel = 1.e-8;
static constexpr amrex::Real small_vel_eps = 1.e-8;
//fixme - need to make this a runtime parameter
static constexpr amrex::Real vel_scale = 1.;

//fixme - need to edit to take vel array and i,j,k. Use mag vel to determine vmax
// this way limiting is a bit more consistent... maybe...
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real calc_small_vel (amrex::Real lo, amrex::Real hi)
{
    amrex::Real small_v;
    amrex::Real vmax = std::max(amrex::Math::abs(lo),amrex::Math::abs(hi));

    small_v = vmax < small_vel_eps*vel_scale ?
                     small_vel_eps*vel_scale : vmax*small_vel_eps;
    return small_v;
}


/**
 * \var covered_val
 *
 * Value to use in covered cells.
 */
static constexpr amrex::Real covered_val = 1.0e40;

#endif
/** @}*/

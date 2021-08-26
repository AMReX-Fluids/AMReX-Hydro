/**
 * \file hydro_constants.H
 */


/** \addtogroup Utilities
 * @{
 */

#ifndef HYDRO_CONSTANTS_H_
#define HYDRO_CONSTANTS_H_

#include <AMReX_REAL.H>
#include <cmath>

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

/**
 * \var covered_val
 *
 * Value to use in covered cells.
 */
static constexpr amrex::Real covered_val = 1.0e40;

#endif
/** @}*/

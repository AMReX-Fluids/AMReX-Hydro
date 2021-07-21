/** \addtogroup Utilities
 * @{
 */

#ifndef HYDRO_CONSTANTS_H_
#define HYDRO_CONSTANTS_H_

#include <AMReX_REAL.H>
#include <cmath>

// Indexing for velocity components
#define XVEL 0
#define YVEL 1
#define ZVEL 2

// What we consider a ~zero velocity when upwinding
static constexpr amrex::Real small_vel = 1.e-8;

// Value to use in covered cells
static constexpr amrex::Real covered_val = 1.0e40;

#endif
/** @}*/

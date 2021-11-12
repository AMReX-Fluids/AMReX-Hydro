
#include <hydro_small_vel.H>
#include <AMReX_ParmParse.H>

namespace HydroSmallVel {

    // Default value
    vel_scale = 1.;

/* application code can override with 
    ParmParse pp("hydro");
    pp.query("vel_scale", HydroSmallVel::vel_scale);
*/
}

    

/**
 * \file hydro_plic.cpp
 * \addtogroup PLIC
 * @{
 */

#include <hydro_plic.H>

namespace PLIC {

// Neither the plane geometry nor the flux volumes account for cut cells or for R-Z volumes
void CheckMesh (amrex::MultiFab const& F, amrex::Geometry const& geom)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(geom.IsCartesian(), "PLIC requires Cartesian coordinates");
#ifdef AMREX_USE_EB
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!F.hasEBFabFactory(), "PLIC does not support embedded boundaries");
#else
    amrex::ignore_unused(F);
#endif
}

} // namespace PLIC
/** @}*/

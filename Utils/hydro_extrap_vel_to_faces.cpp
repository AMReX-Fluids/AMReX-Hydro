/** \addtogroup Utilities
 * @{
 */

#include <hydro_godunov.H>
#include <hydro_mol.H>
#include <hydro_utils.H>

#ifdef AMREX_USE_EB
#include <hydro_ebgodunov.H>
#include <hydro_ebmol.H>
#endif

using namespace amrex;

#ifdef AMREX_USE_EB
void
HydroUtils::ExtrapVelToFaces ( amrex::MultiFab const& vel,
                               amrex::MultiFab const& vel_forces,
                               AMREX_D_DECL(amrex::MultiFab& u_mac,
                                            amrex::MultiFab& v_mac,
                                            amrex::MultiFab& w_mac),
                               amrex::Vector<amrex::BCRec> const& h_bcrec,
                               amrex::BCRec  const* d_bcrec,
                               const amrex::Geometry& geom,
                               amrex::Real dt,
                               const EBFArrayBoxFactory& ebfact,
                               bool godunov_ppm, bool godunov_use_forces_in_trans,
                               std::string advection_type)
{
   ExtrapVelToFaces(vel, vel_forces, AMREX_D_DECL(u_mac,v_mac,w_mac),
                    h_bcrec, d_bcrec, geom, dt,
                    ebfact, /*velocity_on_eb_inflow*/ nullptr,
                    godunov_ppm, godunov_use_forces_in_trans,
                    advection_type);
}
#endif

void
HydroUtils::ExtrapVelToFaces ( amrex::MultiFab const& vel,
                               amrex::MultiFab const& vel_forces,
                               AMREX_D_DECL(amrex::MultiFab& u_mac,
                                            amrex::MultiFab& v_mac,
                                            amrex::MultiFab& w_mac),
                               amrex::Vector<amrex::BCRec> const& h_bcrec,
                               amrex::BCRec  const* d_bcrec,
                               const amrex::Geometry& geom,
                               amrex::Real dt,
#ifdef AMREX_USE_EB
                               const EBFArrayBoxFactory& ebfact,
                               amrex::MultiFab const* velocity_on_eb_inflow,
#endif
                               bool godunov_ppm, bool godunov_use_forces_in_trans,
                               std::string advection_type)
{
    if (advection_type == "Godunov") {
#ifdef AMREX_USE_EB
        if (!ebfact.isAllRegular())
            EBGodunov::ExtrapVelToFaces(vel, vel_forces,
                                        AMREX_D_DECL(u_mac, v_mac, w_mac),
                                        h_bcrec, d_bcrec, geom, dt,
                                        velocity_on_eb_inflow);  // Note that PPM not supported for EB
        else
#endif
            Godunov::ExtrapVelToFaces(vel, vel_forces,
                                      AMREX_D_DECL(u_mac, v_mac, w_mac),
                                      h_bcrec, d_bcrec,
                                      geom, dt, godunov_ppm, godunov_use_forces_in_trans);

    } else if (advection_type == "MOL") {

#ifdef AMREX_USE_EB
        if (!ebfact.isAllRegular())
            EBMOL::ExtrapVelToFaces(vel, AMREX_D_DECL(u_mac, v_mac, w_mac), geom, h_bcrec, d_bcrec);
        else
#endif
            MOL::ExtrapVelToFaces(vel, AMREX_D_DECL(u_mac, v_mac, w_mac), geom, h_bcrec, d_bcrec);
    } else {
        amrex::Abort("Dont know this advection_type in HydroUtils::ExtrapVelToFaces");
    }
}
/** @}*/

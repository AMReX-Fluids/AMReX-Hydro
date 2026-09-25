/** \addtogroup Utilities
 * @{
 */

#include <hydro_godunov.H>
#include <hydro_mol.H>
#include <hydro_utils.H>

#if defined(AMREX_USE_EB) && !defined(HYDRO_NO_EB)
#include <hydro_ebgodunov.H>
#include <hydro_ebmol.H>
#endif

#include <mutex>

using namespace amrex;

#if defined(AMREX_USE_EB) && !defined(HYDRO_NO_EB)
namespace {
    // EBGodunov has no PPM and no forces-in-transverse option, so on an EB level those
    // two flags are dropped. Say so once, rather than silently changing the scheme.
    void WarnOnceAboutDroppedGodunovOptions ()
    {
        static std::once_flag once;
        std::call_once(once, [] {
            amrex::Warning("HydroUtils::ExtrapVelToFaces: on a level with cut cells the "
                           "velocity is extrapolated with EBGodunov, which uses PLM and adds "
                           "the forces after the transverse terms, so godunov_ppm and "
                           "godunov_use_forces_in_trans are ignored.");
        });
    }
}
#endif

#if defined(AMREX_USE_EB) && !defined(HYDRO_NO_EB)
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
                               std::string const& advection_type,
                               int limiter_type,
                               bool allow_inflow_on_outflow)
{
   ExtrapVelToFaces(vel, vel_forces, AMREX_D_DECL(u_mac,v_mac,w_mac),
                    h_bcrec, d_bcrec, geom, dt,
                    ebfact, /*velocity_on_eb_inflow*/ nullptr,
                    godunov_ppm, godunov_use_forces_in_trans,
                    advection_type, limiter_type, allow_inflow_on_outflow);
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
#if defined(AMREX_USE_EB) && !defined(HYDRO_NO_EB)
                               const EBFArrayBoxFactory& ebfact,
                               amrex::MultiFab const* velocity_on_eb_inflow,
#endif
                               bool godunov_ppm, bool godunov_use_forces_in_trans,
                               std::string const& advection_type,
                               int limiter_type,
                               bool allow_inflow_on_outflow,
                               iMultiFab* BC_MF)
{
    // Only (EB)Godunov reads the position-dependent boundary conditions. MOL, EBMOL
    // and BDS see only h_bcrec/d_bcrec, so silently accepting BC_MF for them would
    // replace the mixed boundary condition by the blanket BCRec without any warning.
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(BC_MF == nullptr || advection_type == "Godunov",
                                     "HydroUtils::ExtrapVelToFaces: BC_MF is only supported with (EB)Godunov");

    // BDS extrapolates the velocity to faces with Godunov PLM (see Docs/source/BDS.rst),
    // so that the same advection_type string works here and in ComputeFluxesOnBoxFromState.
    const bool use_ppm = godunov_ppm && (advection_type != "BDS");

    if (advection_type == "Godunov" || advection_type == "BDS") {
#if defined(AMREX_USE_EB) && !defined(HYDRO_NO_EB)
        if (!ebfact.isAllRegular()) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(advection_type != "BDS", "BDS is not available with EB");
            if (use_ppm || godunov_use_forces_in_trans) {
                WarnOnceAboutDroppedGodunovOptions();
            }
            EBGodunov::ExtrapVelToFaces(vel, vel_forces,
                                        AMREX_D_DECL(u_mac, v_mac, w_mac),
                                        h_bcrec, d_bcrec, geom, dt,
                                        velocity_on_eb_inflow,
                                        // Note that PPM is not supported for EB
                                        allow_inflow_on_outflow, BC_MF);
        }
        else
#endif
            Godunov::ExtrapVelToFaces(vel, vel_forces,
                                      AMREX_D_DECL(u_mac, v_mac, w_mac),
                                      h_bcrec, d_bcrec,
                                      geom, dt, use_ppm, godunov_use_forces_in_trans,
                                      limiter_type, allow_inflow_on_outflow, BC_MF);

    } else if (advection_type == "MOL") {

#if defined(AMREX_USE_EB) && !defined(HYDRO_NO_EB)
        if (!ebfact.isAllRegular()) {
            EBMOL::ExtrapVelToFaces(vel, AMREX_D_DECL(u_mac, v_mac, w_mac), geom, h_bcrec, d_bcrec,
                                    allow_inflow_on_outflow);
        }
        else
#endif
            MOL::ExtrapVelToFaces(vel, AMREX_D_DECL(u_mac, v_mac, w_mac), geom, h_bcrec, d_bcrec, allow_inflow_on_outflow);
    } else {
        amrex::Abort("HydroUtils::ExtrapVelToFaces: unknown advection_type " + advection_type
                     + " (expected Godunov, MOL or BDS)");
    }
}
/** @}*/

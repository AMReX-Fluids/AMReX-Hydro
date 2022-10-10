/** \addtogroup Utilities
 * @{
 */

#include <hydro_godunov.H>
#include <hydro_bds.H>
#include <hydro_mol.H>
#include <hydro_utils.H>

#ifdef AMREX_USE_EB
#include <hydro_ebgodunov.H>
#include <hydro_ebmol.H>
#endif

using namespace amrex;

namespace {
    // Limit this function to this file
    void
    ComputeEdgeState (Box const& bx, int ncomp, MFIter& mfi,
                      Array4<Real const> const& q,
                      AMREX_D_DECL(Array4<Real> const& face_x,
                                   Array4<Real> const& face_y,
                                   Array4<Real> const& face_z),
                      AMREX_D_DECL(Array4<Real const> const& u_mac,
                                   Array4<Real const> const& v_mac,
                                   Array4<Real const> const& w_mac),
                      Array4<Real const> const& divu,
                      Array4<Real const> const& fq,
                      Geometry geom, Real l_dt,
                      Vector<BCRec> const& h_bcrec,
                      const BCRec* d_bcrec,
                      int const* iconserv,
#ifdef AMREX_USE_EB
                      const EBFArrayBoxFactory& ebfact,
                      Array4<Real const> const& values_on_eb_inflow,
                      bool regular,
#endif
                      bool godunov_use_ppm, bool godunov_use_forces_in_trans,
                      bool is_velocity,
                      std::string& advection_type)

    {
#ifdef AMREX_USE_EB
        if (!regular)
        {
            EBCellFlagFab const& flagfab = ebfact.getMultiEBCellFlagFab()[mfi];
            Array4<EBCellFlag const> const& flag = flagfab.const_array();

            const auto& vfrac = ebfact.getVolFrac().const_array(mfi);
            const auto& ccc   = ebfact.getCentroid().const_array(mfi);

            AMREX_D_TERM( const auto& apx = ebfact.getAreaFrac()[0]->const_array(mfi);,
                          const auto& apy = ebfact.getAreaFrac()[1]->const_array(mfi);,
                          const auto& apz = ebfact.getAreaFrac()[2]->const_array(mfi););

            AMREX_D_TERM( const auto& fcx = ebfact.getFaceCent()[0]->const_array(mfi);,
                          const auto& fcy = ebfact.getFaceCent()[1]->const_array(mfi);,
                          const auto& fcz = ebfact.getFaceCent()[2]->const_array(mfi););

            if (advection_type == "MOL")
            {
                EBMOL::ComputeEdgeState( bx,
                                         AMREX_D_DECL(face_x,face_y,face_z),
                                         q, ncomp,
                                         AMREX_D_DECL(u_mac,v_mac,w_mac),
                                         geom.Domain(), h_bcrec, d_bcrec,
                                         AMREX_D_DECL(fcx,fcy,fcz),
                                         ccc, vfrac, flag,
                                         is_velocity);
            }
            else if (advection_type == "Godunov")
            {
                int ngrow = 4; // NOT SURE ABOUT THIS
                FArrayBox tmpfab_v(amrex::grow(bx,ngrow),  (4*AMREX_SPACEDIM + 2)*ncomp,
                                   The_Async_Arena());
                EBGodunov::ComputeEdgeState(bx, ncomp, q,
                                            AMREX_D_DECL(face_x,face_y,face_z),
                                            AMREX_D_DECL(u_mac,v_mac,w_mac),
                                            divu, fq,
                                            geom, l_dt,
                                            h_bcrec, d_bcrec, iconserv,
                                            tmpfab_v.dataPtr(), flag,
                                            AMREX_D_DECL(apx,apy,apz), vfrac,
                                            AMREX_D_DECL(fcx,fcy,fcz), ccc,
                                            is_velocity,
                                            values_on_eb_inflow);
            }
            else if (advection_type == "BDS")
            {
                Abort("BDS is not available with EB");
            }
            else
            {
                Abort("Unknown advection_type: "+advection_type);
            }
        }
        else
#endif
        {
            if (advection_type == "MOL")
            {
                MOL::ComputeEdgeState( bx,
                                       AMREX_D_DECL(face_x,face_y,face_z),
                                       q, ncomp,
                                       AMREX_D_DECL(u_mac,v_mac,w_mac),
                                       geom.Domain(), h_bcrec, d_bcrec,
                                       is_velocity);
            }
            else if (advection_type == "Godunov")
            {
                Godunov::ComputeEdgeState(bx, ncomp, q,
                                          AMREX_D_DECL(face_x,face_y,face_z),
                                          AMREX_D_DECL(u_mac,v_mac,w_mac),
                                          divu, fq,
                                          geom,
                                          l_dt, d_bcrec, iconserv,
                                          godunov_use_ppm, godunov_use_forces_in_trans,
                                          is_velocity);
            }
            else if (advection_type == "BDS")
            {
                BDS::ComputeEdgeState( bx, ncomp, q,
                                       AMREX_D_DECL(face_x,face_y,face_z),
                                       AMREX_D_DECL(u_mac,v_mac,w_mac),
                                       divu, fq, geom,
                                       l_dt, d_bcrec, iconserv,
                                       is_velocity);
            }
            else
            {
                Abort("Unknown advection_type: "+advection_type);
            }
        }
    }
}

#ifdef AMREX_USE_EB
void
HydroUtils::ComputeFluxesOnBoxFromState (Box const& bx, int ncomp, MFIter& mfi,
                                         Array4<Real const> const& q,
                                         AMREX_D_DECL(Array4<Real> const& flux_x,
                                                      Array4<Real> const& flux_y,
                                                      Array4<Real> const& flux_z),
                                         AMREX_D_DECL(Array4<Real> const& face_x,
                                                      Array4<Real> const& face_y,
                                                      Array4<Real> const& face_z),
                                         bool knownFaceState,
                                         AMREX_D_DECL(Array4<Real const> const& u_mac,
                                                      Array4<Real const> const& v_mac,
                                                      Array4<Real const> const& w_mac),
                                         Array4<Real const> const& divu,
                                         Array4<Real const> const& fq,
                                         Geometry geom, Real l_dt,
                                         Vector<BCRec> const& h_bcrec,
                                         const BCRec* d_bcrec,
                                         int const* iconserv,
                                         const EBFArrayBoxFactory& ebfact,
                                         bool godunov_use_ppm, bool godunov_use_forces_in_trans,
                                         bool is_velocity, bool fluxes_are_area_weighted,
                                         std::string& advection_type)

{
    ComputeFluxesOnBoxFromState(bx, ncomp, mfi, q,
                                AMREX_D_DECL(flux_x, flux_y, flux_z),
                                AMREX_D_DECL(face_x, face_y, face_z),
                                knownFaceState,
                                AMREX_D_DECL(u_mac, v_mac, w_mac),
                                divu, fq, geom, l_dt, h_bcrec, d_bcrec, iconserv,
                                ebfact, /*values_on_eb_inflow*/ Array4<Real const>{},
                                godunov_use_ppm, godunov_use_forces_in_trans,
                                is_velocity, fluxes_are_area_weighted, advection_type);

}
#endif


void
HydroUtils::ComputeFluxesOnBoxFromState (Box const& bx, int ncomp, MFIter& mfi,
                                         Array4<Real const> const& q,
                                         AMREX_D_DECL(Array4<Real> const& flux_x,
                                                      Array4<Real> const& flux_y,
                                                      Array4<Real> const& flux_z),
                                         AMREX_D_DECL(Array4<Real> const& face_x,
                                                      Array4<Real> const& face_y,
                                                      Array4<Real> const& face_z),
                                         bool knownFaceState,
                                         AMREX_D_DECL(Array4<Real const> const& u_mac,
                                                      Array4<Real const> const& v_mac,
                                                      Array4<Real const> const& w_mac),
                                         Array4<Real const> const& divu,
                                         Array4<Real const> const& fq,
                                         Geometry geom, Real l_dt,
                                         Vector<BCRec> const& h_bcrec,
                                         const BCRec* d_bcrec,
                                         int const* iconserv,
#ifdef AMREX_USE_EB
                                         const EBFArrayBoxFactory& ebfact,
                                         Array4<Real const> const& values_on_eb_inflow,
#endif
                                         bool godunov_use_ppm, bool godunov_use_forces_in_trans,
                                         bool is_velocity, bool fluxes_are_area_weighted,
                                         std::string& advection_type)

{
    ComputeFluxesOnBoxFromState(bx, ncomp, mfi, q,
                                AMREX_D_DECL(flux_x, flux_y, flux_z),
                                AMREX_D_DECL(face_x, face_y, face_z),
                                knownFaceState,
                                AMREX_D_DECL(u_mac, v_mac, w_mac),
                                AMREX_D_DECL(u_mac, v_mac, w_mac),
                                divu, fq, geom, l_dt, h_bcrec, d_bcrec, iconserv,
#ifdef AMREX_USE_EB
                                ebfact, values_on_eb_inflow,
#endif
                                godunov_use_ppm, godunov_use_forces_in_trans,
                                is_velocity, fluxes_are_area_weighted, advection_type);

}

void
HydroUtils::ComputeFluxesOnBoxFromState (Box const& bx, int ncomp, MFIter& mfi,
                                         Array4<Real const> const& q,
                                         AMREX_D_DECL(Array4<Real> const& flux_x,
                                                      Array4<Real> const& flux_y,
                                                      Array4<Real> const& flux_z),
                                         AMREX_D_DECL(Array4<Real> const& face_x,
                                                      Array4<Real> const& face_y,
                                                      Array4<Real> const& face_z),
                                         bool knownFaceState,
                                         AMREX_D_DECL(Array4<Real const> const& u_mac,
                                                      Array4<Real const> const& v_mac,
                                                      Array4<Real const> const& w_mac),
                                         AMREX_D_DECL(Array4<Real const> const& u_flux,
                                                      Array4<Real const> const& v_flux,
                                                      Array4<Real const> const& w_flux),
                                         Array4<Real const> const& divu,
                                         Array4<Real const> const& fq,
                                         Geometry geom, Real l_dt,
                                         Vector<BCRec> const& h_bcrec,
                                         const BCRec* d_bcrec,
                                         int const* iconserv,
#ifdef AMREX_USE_EB
                                         const EBFArrayBoxFactory& ebfact,
                                         Array4<Real const> const& values_on_eb_inflow,
#endif
                                         bool godunov_use_ppm, bool godunov_use_forces_in_trans,
                                         bool is_velocity, bool fluxes_are_area_weighted,
                                         std::string& advection_type)

{
#ifdef AMREX_USE_EB
    EBCellFlagFab const& flagfab = ebfact.getMultiEBCellFlagFab()[mfi];
    Array4<EBCellFlag const> const& flag = flagfab.const_array();

    // If entire box is covered, don't do anything and return
    if (flagfab.getType(bx) == FabType::covered)
        return;

    //FIXME? -- Godunov needs to check on grow 3, but MOL only needs 2
    bool regular = (flagfab.getType(amrex::grow(bx,3)) == FabType::regular);
#endif

    // Compute edge state if needed
    if (!knownFaceState) {
        ComputeEdgeState(bx, ncomp, mfi, q,
                         AMREX_D_DECL(face_x,face_y,face_z),
                         AMREX_D_DECL(u_mac,v_mac,w_mac),
                         divu, fq,
                         geom, l_dt,
                         h_bcrec, d_bcrec, iconserv,
#ifdef AMREX_USE_EB
                         ebfact, values_on_eb_inflow, regular,
#endif
                         godunov_use_ppm, godunov_use_forces_in_trans,
                         is_velocity, advection_type);
    }

    // Compute fluxes.
    // For a typical advection step, the velocity here is the u_mac above.
    // For a multilevel synchronization, the velocity here is the "corrective" velocity.
#ifdef AMREX_USE_EB
    if (!regular) {
        AMREX_D_TERM(Array4<Real const> const& apx = ebfact.getAreaFrac()[0]->const_array(mfi);,
                     Array4<Real const> const& apy = ebfact.getAreaFrac()[1]->const_array(mfi);,
                     Array4<Real const> const& apz = ebfact.getAreaFrac()[2]->const_array(mfi););

        HydroUtils::EB_ComputeFluxes( bx,
                                      AMREX_D_DECL(flux_x,flux_y,flux_z),
                                      AMREX_D_DECL(u_flux,v_flux,w_flux),
                                      AMREX_D_DECL(face_x,face_y,face_z),
                                      AMREX_D_DECL(apx,apy,apz),
                                      geom, ncomp,
                                      flag, fluxes_are_area_weighted);
    } else
#endif
    {
        HydroUtils::ComputeFluxes( bx,
                                   AMREX_D_DECL(flux_x,flux_y,flux_z),
                                   AMREX_D_DECL(u_flux,v_flux,w_flux),
                                   AMREX_D_DECL(face_x,face_y,face_z),
                                   geom, ncomp, fluxes_are_area_weighted );
    }
}

/** @}*/

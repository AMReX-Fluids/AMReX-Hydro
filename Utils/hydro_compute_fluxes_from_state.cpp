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

void
HydroUtils::ComputeFluxesOnBoxFromState (
                            Box const& bx, int ncomp, MFIter& mfi,
                            Array4<Real const> const& q,
                            AMREX_D_DECL(Array4<Real> const& flux_x,
                                         Array4<Real> const& flux_y,
                                         Array4<Real> const& flux_z),
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
#endif
                            bool godunov_use_ppm, bool godunov_use_forces_in_trans,
                            bool is_velocity, bool fluxes_are_area_weighted,
                            std::string& advection_type)

{
#ifdef AMREX_USE_EB
            EBCellFlagFab const& flagfab = ebfact.getMultiEBCellFlagFab()[mfi];
            Array4<EBCellFlag const> const& flag = flagfab.const_array();

            Array4<Real const> AMREX_D_DECL(fcx, fcy, fcz), ccc, vfrac, AMREX_D_DECL(apx, apy, apz);

            // If entire box is covered, don't do anything and return
            if (flagfab.getType(bx) == FabType::covered)
               return;

	    //FIXME? -- Godunov needs to check on grow 3, but MOL only needs 2
            bool regular = (flagfab.getType(amrex::grow(bx,3)) == FabType::regular);

            if (!regular)
            {
                vfrac = ebfact.getVolFrac().const_array(mfi);
                ccc   = ebfact.getCentroid().const_array(mfi);

                AMREX_D_TERM( apx = ebfact.getAreaFrac()[0]->const_array(mfi);,
                              apy = ebfact.getAreaFrac()[1]->const_array(mfi);,
                              apz = ebfact.getAreaFrac()[2]->const_array(mfi););

                AMREX_D_TERM( fcx = ebfact.getFaceCent()[0]->const_array(mfi);,
                              fcy = ebfact.getFaceCent()[1]->const_array(mfi);,
                              fcz = ebfact.getFaceCent()[2]->const_array(mfi););
            }
#endif
            if (advection_type == "MOL")
            {
#ifdef AMREX_USE_EB
              if (!regular) 
                EBMOL::ComputeEdgeState( bx, 
                                       AMREX_D_DECL(face_x,face_y,face_z),
                                       q, ncomp, 
                                       AMREX_D_DECL(u_mac,v_mac,w_mac),
                                       geom.Domain(), h_bcrec, d_bcrec,
                                       AMREX_D_DECL(fcx,fcy,fcz),
                                       ccc, vfrac, flag,
                                       is_velocity);
            else
#endif
                MOL::ComputeEdgeState( bx, 
                                       AMREX_D_DECL(face_x,face_y,face_z),
                                       q, ncomp, 
                                       AMREX_D_DECL(u_mac,v_mac,w_mac),
                                       geom.Domain(), h_bcrec, d_bcrec,
                                       is_velocity);

            } else if (advection_type == "Godunov") {

              int ngrow = 4; // NOT SURE ABOUT THIS
              FArrayBox tmpfab_v(amrex::grow(bx,ngrow),  (4*AMREX_SPACEDIM + 2)*ncomp);
              Elixir    eli = tmpfab_v.elixir();
#ifdef AMREX_USE_EB
              if (!regular) 
                EBGodunov::ComputeEdgeState(
                                       bx, ncomp, q,
                                       AMREX_D_DECL(face_x,face_y,face_z),
                                       AMREX_D_DECL(u_mac,v_mac,w_mac),
                                       divu, fq,
                                       geom, l_dt, 
                                       h_bcrec, d_bcrec, iconserv,
                                       tmpfab_v.dataPtr(), flag,
                                       AMREX_D_DECL(apx,apy,apz), vfrac,
                                       AMREX_D_DECL(fcx,fcy,fcz), ccc,
                                       is_velocity);
              else
#endif

                Godunov::ComputeEdgeState( 
                                       bx, ncomp, q,
                                       AMREX_D_DECL(face_x,face_y,face_z),
                                       AMREX_D_DECL(u_mac,v_mac,w_mac),
                                       divu, fq, 
                                       geom,
                                       l_dt, d_bcrec, iconserv,
                                       godunov_use_ppm, godunov_use_forces_in_trans,
                                       is_velocity);
            } // Godunov

            // Compute fluxes
#ifdef AMREX_USE_EB
            if (!regular) 
                HydroUtils::EB_ComputeFluxes( bx,
                                             AMREX_D_DECL(flux_x,flux_y,flux_z),
                                             AMREX_D_DECL(u_mac,v_mac,w_mac),
                                             AMREX_D_DECL(face_x,face_y,face_z),
                                             AMREX_D_DECL(apx,apy,apz),
                                             geom, ncomp,
                                             flag, fluxes_are_area_weighted);
            else
#endif

                HydroUtils::ComputeFluxes( bx,
                                           AMREX_D_DECL(flux_x,flux_y,flux_z),
                                           AMREX_D_DECL(u_mac,v_mac,w_mac),
                                           AMREX_D_DECL(face_x,face_y,face_z),
                                           geom, ncomp, fluxes_are_area_weighted );
}
/** @}*/

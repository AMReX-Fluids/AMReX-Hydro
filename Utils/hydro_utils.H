/** \addtogroup Utilities
 * @{
 */

#ifndef HYDRO_UTILS_H
#define HYDRO_UTILS_H

#include <AMReX_MultiFabUtil.H>
#include <AMReX_BCRec.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBFabFactory.H>
#include <AMReX_EBMultiFabUtil.H>
#endif

/**
 * Collection of hydro utility functions
 *
 */

namespace HydroUtils {

/**
 * \brief Compute edge state and flux. Most general version for use with multilevel synchonization.
 *
 */
void
ComputeFluxesOnBoxFromState ( amrex::Box const& bx, int ncomp, amrex::MFIter& mfi,
                              amrex::Array4<amrex::Real const> const& q,
                              AMREX_D_DECL(amrex::Array4<amrex::Real> const& flux_x,
                                           amrex::Array4<amrex::Real> const& flux_y,
                                           amrex::Array4<amrex::Real> const& flux_z),
                              AMREX_D_DECL(amrex::Array4<amrex::Real> const& xface,
                                           amrex::Array4<amrex::Real> const& yface,
                                           amrex::Array4<amrex::Real> const& zface),
                              bool knownFaceState,
                              AMREX_D_DECL(amrex::Array4<amrex::Real const> const& umac,
                                           amrex::Array4<amrex::Real const> const& vmac,
                                           amrex::Array4<amrex::Real const> const& wmac),
                              AMREX_D_DECL(amrex::Array4<amrex::Real const> const& uflux,
                                           amrex::Array4<amrex::Real const> const& vflux,
                                           amrex::Array4<amrex::Real const> const& wflux),
                              amrex::Array4<amrex::Real const> const& divu,
                              amrex::Array4<amrex::Real const> const& fq,
                              amrex::Geometry geom,
                              amrex::Real l_dt,
                              amrex::Vector<amrex::BCRec> const& h_bcrec,
                              const amrex::BCRec* d_bcrec,
                              int const* iconserv,
#ifdef AMREX_USE_EB
                              const amrex::EBFArrayBoxFactory& ebfact,
                              amrex::Array4<amrex::Real const> const& values_on_eb_inflow,
#endif
                              bool godunov_use_ppm, bool godunov_use_forces_in_trans,
                              bool is_velocity, bool fluxes_are_area_weighted,
                              std::string& advection_type);
/**
 * \brief Compute edge state and flux. For typical advection, and also allows for inflow on EB.
 *
 */
void
ComputeFluxesOnBoxFromState ( amrex::Box const& bx, int ncomp, amrex::MFIter& mfi,
                              amrex::Array4<amrex::Real const> const& q,
                              AMREX_D_DECL(amrex::Array4<amrex::Real> const& flux_x,
                                           amrex::Array4<amrex::Real> const& flux_y,
                                           amrex::Array4<amrex::Real> const& flux_z),
                              AMREX_D_DECL(amrex::Array4<amrex::Real> const& xface,
                                           amrex::Array4<amrex::Real> const& yface,
                                           amrex::Array4<amrex::Real> const& zface),
                              bool knownFaceState,
                              AMREX_D_DECL(amrex::Array4<amrex::Real const> const& umac,
                                           amrex::Array4<amrex::Real const> const& vmac,
                                           amrex::Array4<amrex::Real const> const& wmac),
                              amrex::Array4<amrex::Real const> const& divu,
                              amrex::Array4<amrex::Real const> const& fq,
                              amrex::Geometry geom,
                              amrex::Real l_dt,
                              amrex::Vector<amrex::BCRec> const& h_bcrec,
                              const amrex::BCRec* d_bcrec,
                              int const* iconserv,
#ifdef AMREX_USE_EB
                              const amrex::EBFArrayBoxFactory& ebfact,
                              amrex::Array4<amrex::Real const> const& values_on_eb_inflow,
#endif
                              bool godunov_use_ppm, bool godunov_use_forces_in_trans,
                              bool is_velocity, bool fluxes_are_area_weighted,
                              std::string& advection_type);

/**
 * \brief Compute edge state and flux. For typical advection, but no inflow through EB.
 *
 */
#ifdef AMREX_USE_EB
void
ComputeFluxesOnBoxFromState ( amrex::Box const& bx, int ncomp, amrex::MFIter& mfi,
                             amrex::Array4<amrex::Real const> const& q,
                             AMREX_D_DECL(amrex::Array4<amrex::Real> const& flux_x,
                                          amrex::Array4<amrex::Real> const& flux_y,
                                          amrex::Array4<amrex::Real> const& flux_z),
                             AMREX_D_DECL(amrex::Array4<amrex::Real> const& xface,
                                          amrex::Array4<amrex::Real> const& yface,
                                          amrex::Array4<amrex::Real> const& zface),
                             bool knownFaceState,
                             AMREX_D_DECL(amrex::Array4<amrex::Real const> const& umac,
                                          amrex::Array4<amrex::Real const> const& vmac,
                                          amrex::Array4<amrex::Real const> const& wmac),
                             amrex::Array4<amrex::Real const> const& divu,
                             amrex::Array4<amrex::Real const> const& fq,
                             amrex::Geometry geom,
                             amrex::Real l_dt,
                             amrex::Vector<amrex::BCRec> const& h_bcrec,
                             const amrex::BCRec* d_bcrec,
                             int const* iconserv,
                             const amrex::EBFArrayBoxFactory& ebfact,
                             bool godunov_use_ppm, bool godunov_use_forces_in_trans,
                             bool is_velocity, bool fluxes_are_area_weighted,
                             std::string& advection_type);
#endif

#ifdef AMREX_USE_EB
void
ExtrapVelToFaces ( amrex::MultiFab const& vel,
                   amrex::MultiFab const& vel_forces,
                   AMREX_D_DECL(amrex::MultiFab& u_mac,
                                amrex::MultiFab& v_mac,
                                amrex::MultiFab& w_mac),
                   amrex::Vector<amrex::BCRec> const& h_bcrec,
                   amrex::BCRec  const* d_bcrec,
                   const amrex::Geometry& geom,
                   amrex::Real l_dt,
                   const amrex::EBFArrayBoxFactory& ebfact,
                   bool godunov_use_ppm, bool godunov_use_forces_in_trans,
                   std::string advection_type);
#endif

void
ExtrapVelToFaces ( amrex::MultiFab const& vel,
                   amrex::MultiFab const& vel_forces,
                   AMREX_D_DECL(amrex::MultiFab& u_mac,
                                amrex::MultiFab& v_mac,
                                amrex::MultiFab& w_mac),
                   amrex::Vector<amrex::BCRec> const& h_bcrec,
                   amrex::BCRec  const* d_bcrec,
                   const amrex::Geometry& geom,
                   amrex::Real l_dt,
#ifdef AMREX_USE_EB
                   const amrex::EBFArrayBoxFactory& ebfact,
                   amrex::MultiFab const* velocity_on_eb_inflow,
#endif
                   bool godunov_use_ppm, bool godunov_use_forces_in_trans,
                   std::string advection_type);

/**
 * \brief If convective, compute convTerm = u dot grad q = div (u q) - q div(u).
 *
 */
void
ComputeConvectiveTerm ( amrex::Box const& bx, int ncomp, amrex::MFIter& mfi,
                        amrex::Array4<amrex::Real const> const& q,
                        AMREX_D_DECL(amrex::Array4<amrex::Real const> const& xface,
                                     amrex::Array4<amrex::Real const> const& yface,
                                     amrex::Array4<amrex::Real const> const& zface),
                        amrex::Array4<amrex::Real const> const& divu,
                        amrex::Array4<amrex::Real> const& convTerm,
                        int const* iconserv,
#ifdef AMREX_USE_EB
                        const amrex::EBFArrayBoxFactory& ebfact,
#endif
                        std::string& advection_type);

/**
 * \brief Compute Fluxes.
 *
 */
void ComputeFluxes ( amrex::Box const& bx,
                     AMREX_D_DECL( amrex::Array4<amrex::Real> const& flux_x,
                                   amrex::Array4<amrex::Real> const& flux_y,
                                   amrex::Array4<amrex::Real> const& flux_z),
                     AMREX_D_DECL( amrex::Array4<amrex::Real const> const& umac,
                                   amrex::Array4<amrex::Real const> const& vmac,
                                   amrex::Array4<amrex::Real const> const& wmac),
                     AMREX_D_DECL( amrex::Array4<amrex::Real const> const& xface,
                                   amrex::Array4<amrex::Real const> const& yface,
                                   amrex::Array4<amrex::Real const> const& zface),
                     amrex::Geometry const& geom, const int ncomp,
                     bool fluxes_are_area_weighted);

/**
 * \brief Compute divergence.
 *
 */

void ComputeDivergence ( amrex::Box const& bx,
                         amrex::Array4<amrex::Real> const& div,
                         AMREX_D_DECL( amrex::Array4<amrex::Real const> const& flux_x,
                                       amrex::Array4<amrex::Real const> const& flux_y,
                                       amrex::Array4<amrex::Real const> const& flux_z),
                         const int ncomp, amrex::Geometry const& geom,
                         const amrex::Real mult,
                         bool fluxes_are_area_weighted);

#ifdef AMREX_USE_EB

void EB_ComputeFluxes ( amrex::Box const& bx,
                        AMREX_D_DECL( amrex::Array4<amrex::Real> const& fx,
                                      amrex::Array4<amrex::Real> const& fy,
                                      amrex::Array4<amrex::Real> const& fz),
                        AMREX_D_DECL( amrex::Array4<amrex::Real const> const& umac,
                                      amrex::Array4<amrex::Real const> const& vmac,
                                      amrex::Array4<amrex::Real const> const& wmac),
                        AMREX_D_DECL( amrex::Array4<amrex::Real const> const& xedge,
                                      amrex::Array4<amrex::Real const> const& yedge,
                                      amrex::Array4<amrex::Real const> const& zedge),
                        AMREX_D_DECL( amrex::Array4<amrex::Real const> const& apx,
                                      amrex::Array4<amrex::Real const> const& apy,
                                      amrex::Array4<amrex::Real const> const& apz),
                        amrex::Geometry const& geom, const int ncomp,
                        amrex::Array4<amrex::EBCellFlag const> const& flag,
                        bool fluxes_are_area_weighted);


void EB_ComputeDivergence ( amrex::Box const& bx,
                            amrex::Array4<amrex::Real> const& div,
                            AMREX_D_DECL( amrex::Array4<amrex::Real const> const& fx,
                                          amrex::Array4<amrex::Real const> const& fy,
                                          amrex::Array4<amrex::Real const> const& fz),
                            amrex::Array4<amrex::Real const> const& vfrac,
                            const int ncomp, amrex::Geometry const& geom,
                            const amrex::Real mult,
                            bool fluxes_are_area_weighted);


void EB_ComputeDivergence ( amrex::Box const& bx,
                            amrex::Array4<amrex::Real> const& div,
                            AMREX_D_DECL( amrex::Array4<amrex::Real const> const& fx,
                                          amrex::Array4<amrex::Real const> const& fy,
                                          amrex::Array4<amrex::Real const> const& fz),
                            amrex::Array4<amrex::Real const> const& vfrac,
                            const int ncomp, amrex::Geometry const& geom,
                            const amrex::Real mult,
                            bool fluxes_are_area_weighted,
                            amrex::Array4<amrex::Real const> const& eb_velocity,
                            amrex::Array4<amrex::Real const> const& values_on_eb_inflow,
                            amrex::Array4<amrex::EBCellFlag const> const& flag_arr,
                            amrex::Array4<amrex::Real const> const& barea,
                            amrex::Array4<amrex::Real const> const& bnorm);
#endif
}

#endif
/** @}*/

/** \addtogroup Utilities
 * @{
 */

#ifndef HYDRO_UTILS_H
#define HYDRO_UTILS_H

#include <AMReX_MultiFabUtil.H>
#include <AMReX_BCRec.H>
#include <hydro_godunov_ppm.H>

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
 * \brief Enforces solvablity by scaling outflow to match with inflow.
 *
 */
void enforceInOutSolvability (
    const amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>>& vels_vec,
    amrex::BCRec const* bc_type,
    const amrex::Vector<amrex::Geometry>& geom,
    bool include_bndry_corners = false);

/**
 * \brief Compute edge state and flux. Most general version for use with multilevel synchronization.
 *        All other versions ultimately call this one.
 */
void
ComputeFluxesOnBoxFromState ( amrex::Box const& bx, int ncomp, amrex::MFIter& mfi,
                              amrex::Array4<amrex::Real const> const& q,
                              AMREX_D_DECL(amrex::Array4<amrex::Real> const& flux_x,
                                           amrex::Array4<amrex::Real> const& flux_y,
                                           amrex::Array4<amrex::Real> const& flux_z),
                              AMREX_D_DECL(amrex::Array4<amrex::Real> const& face_x,
                                           amrex::Array4<amrex::Real> const& face_y,
                                           amrex::Array4<amrex::Real> const& face_z),
                              bool knownFaceState,
                              AMREX_D_DECL(amrex::Array4<amrex::Real const> const& u_mac,
                                           amrex::Array4<amrex::Real const> const& v_mac,
                                           amrex::Array4<amrex::Real const> const& w_mac),
                              AMREX_D_DECL(amrex::Array4<amrex::Real const> const& u_flux,
                                           amrex::Array4<amrex::Real const> const& v_flux,
                                           amrex::Array4<amrex::Real const> const& w_flux),
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
                              std::string const& advection_type,
                              int limiter_type = PPM::default_limiter,
                              amrex::Array4<int const> const& bc_arr = {});
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
                              AMREX_D_DECL(amrex::Array4<amrex::Real> const& face_x,
                                           amrex::Array4<amrex::Real> const& face_y,
                                           amrex::Array4<amrex::Real> const& face_z),
                              bool knownFaceState,
                              AMREX_D_DECL(amrex::Array4<amrex::Real const> const& u_mac,
                                           amrex::Array4<amrex::Real const> const& v_mac,
                                           amrex::Array4<amrex::Real const> const& w_mac),
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
                              std::string const& advection_type,
                              int limiter_type = PPM::default_limiter,
                              amrex::Array4<int const> const& bc_arr = {});

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
                             AMREX_D_DECL(amrex::Array4<amrex::Real> const& face_x,
                                          amrex::Array4<amrex::Real> const& face_y,
                                          amrex::Array4<amrex::Real> const& face_z),
                             bool knownFaceState,
                             AMREX_D_DECL(amrex::Array4<amrex::Real const> const& u_mac,
                                          amrex::Array4<amrex::Real const> const& v_mac,
                                          amrex::Array4<amrex::Real const> const& w_mac),
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
                             std::string const& advection_type,
                              int limiter_type = PPM::default_limiter,
                              amrex::Array4<int const> const& bc_arr = {});
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
                   amrex::Real dt,
                   const amrex::EBFArrayBoxFactory& ebfact,
                   bool godunov_ppm, bool godunov_use_forces_in_trans,
                   std::string const& advection_type,
                   int limiter_type = PPM::default_limiter);
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
                   amrex::Real dt,
#ifdef AMREX_USE_EB
                   const amrex::EBFArrayBoxFactory& ebfact,
                   amrex::MultiFab const* velocity_on_eb_inflow,
#endif
                   bool godunov_ppm, bool godunov_use_forces_in_trans,
                   std::string const& advection_type,
                   int limiter_type = PPM::default_limiter,
                   amrex::iMultiFab* BC_MF = nullptr);

/**
 * \brief If convective, compute convTerm = u dot grad q = div (u q) - q div(u).
 *
 */
void
ComputeConvectiveTerm ( amrex::Box const& bx, int num_comp, amrex::MFIter& mfi,
                        amrex::Array4<amrex::Real const> const& q,
                        AMREX_D_DECL(amrex::Array4<amrex::Real const> const& q_on_face_x,
                                     amrex::Array4<amrex::Real const> const& q_on_face_y,
                                     amrex::Array4<amrex::Real const> const& q_on_face_z),
                        amrex::Array4<amrex::Real const> const& divu,
                        amrex::Array4<amrex::Real> const& convTerm,
                        int const* iconserv,
#ifdef AMREX_USE_EB
                        const amrex::EBFArrayBoxFactory& ebfact,
#endif
                        std::string const& advection_type);

/**
 * \brief Compute Fluxes.
 *
 */
void ComputeFluxes ( amrex::Box const& bx,
                     AMREX_D_DECL( amrex::Array4<amrex::Real> const& fx,
                                   amrex::Array4<amrex::Real> const& fy,
                                   amrex::Array4<amrex::Real> const& fz),
                     AMREX_D_DECL( amrex::Array4<amrex::Real const> const& umac,
                                   amrex::Array4<amrex::Real const> const& vmac,
                                   amrex::Array4<amrex::Real const> const& wmac),
                     AMREX_D_DECL( amrex::Array4<amrex::Real const> const& xed,
                                   amrex::Array4<amrex::Real const> const& yed,
                                   amrex::Array4<amrex::Real const> const& zed),
                     amrex::Geometry const& geom, int ncomp,
                     bool fluxes_are_area_weighted);

/**
 * \brief Compute divergence.
 *
 */

void ComputeDivergence ( amrex::Box const& bx,
                         amrex::Array4<amrex::Real> const& div,
                         AMREX_D_DECL( amrex::Array4<amrex::Real const> const& fx,
                                       amrex::Array4<amrex::Real const> const& fy,
                                       amrex::Array4<amrex::Real const> const& fz),
                         int ncomp, amrex::Geometry const& geom,
                         amrex::Real mult,
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
                        amrex::Geometry const& geom, int ncomp,
                        amrex::Array4<amrex::EBCellFlag const> const& flag,
                        bool fluxes_are_area_weighted);


void EB_ComputeDivergence ( amrex::Box const& bx,
                            amrex::Array4<amrex::Real> const& div,
                            AMREX_D_DECL( amrex::Array4<amrex::Real const> const& fx,
                                          amrex::Array4<amrex::Real const> const& fy,
                                          amrex::Array4<amrex::Real const> const& fz),
                            amrex::Array4<amrex::Real const> const& vfrac,
                            int ncomp, amrex::Geometry const& geom,
                            amrex::Real mult,
                            bool fluxes_are_area_weighted);


void EB_ComputeDivergence ( amrex::Box const& bx,
                            amrex::Array4<amrex::Real> const& div,
                            AMREX_D_DECL( amrex::Array4<amrex::Real const> const& fx,
                                          amrex::Array4<amrex::Real const> const& fy,
                                          amrex::Array4<amrex::Real const> const& fz),
                            amrex::Array4<amrex::Real const> const& vfrac,
                            int ncomp, amrex::Geometry const& geom,
                            amrex::Real mult,
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

/** \addtogroup Utilities
 * @{
 */

#ifndef HYDRO_UTILS_H
#define HYDRO_UTILS_H

#include <AMReX_MultiFabUtil.H>
#include <AMReX_BCRec.H>
#include <hydro_godunov_ppm.H>

#if defined(AMREX_USE_EB) && !defined(HYDRO_NO_EB)
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
 * Only boundaries with BCType::direction_dependent take part. The inflow and
 * outflow are integrated over the tagged boundary faces, and the outflow is then
 * scaled by influx/outflux. The test for "no flow at all" compares the mean
 * normal speed over the tagged faces against small_vel, so it does not depend on
 * the length scale of the problem.
 *
 * \param[in,out] vels_vec  Per level, the normal velocity in each direction. Each
 *                          MultiFab is read and written in its boundary ghost
 *                          cells, so it needs (filled) ghost cells: one in the
 *                          normal direction if it is cell-centered, and one in each
 *                          tangential direction if include_bndry_corners is true.
 *                          This is asserted.
 * \param[in]     bc_type   Boundary conditions for the normal velocity.
 * \param[in]     geom      Per level geometry.
 * \param[in]     include_bndry_corners  Also tag and correct the cells in the
 *                          corners of the domain boundary.
 */
void enforceInOutSolvability (
    const amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>>& vels_vec,
    amrex::BCRec const* bc_type,
    const amrex::Vector<amrex::Geometry>& geom,
    bool include_bndry_corners = false);

/**
 * \brief Compute edge state and flux. Most general version for use with multilevel synchronization.
 *        All other versions ultimately call this one.
 *        Allows time-depdendent inflow and inflow through EB.
 *
 * \note See ComputeFluxes for the meaning of iconserv, in particular the
 *       iconserv[n] == -1 convention, which means the face array already holds
 *       the final flux and bypasses every geometric factor. It requires
 *       knownFaceState = true.
 *
 * \note Ghost-cell requirements, when knownFaceState is false and the edge state
 *       is predicted here: q needs 3 filled ghost cells (4 for 3D EBGodunov, and
 *       BDS needs the corners as well), the MAC velocities need 1, and divu needs
 *       1. divu must be a valid Array4 in 3D Godunov for every component, and for
 *       BDS whenever any component is conservative.
 */
void
ComputeFluxesOnBoxFromState ( amrex::Box const& bx, int ncomp, amrex::MFIter& mfi,
                              amrex::Array4<amrex::Real const> const& q,
                              amrex::Array4<amrex::Real const> const& qnph,
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
#if defined(AMREX_USE_EB) && !defined(HYDRO_NO_EB)
                              const amrex::EBFArrayBoxFactory& ebfact,
                              amrex::Array4<amrex::Real const> const& values_on_eb_inflow,
#endif
                              bool godunov_use_ppm, bool godunov_use_forces_in_trans,
                              bool is_velocity, bool fluxes_are_area_weighted,
                              std::string const& advection_type,
                              int limiter_type = PPM::default_limiter,
                              bool allow_inflow_on_outflow = false,
                              amrex::Array4<int const> const& bc_arr = {});

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
#if defined(AMREX_USE_EB) && !defined(HYDRO_NO_EB)
                              const amrex::EBFArrayBoxFactory& ebfact,
                              amrex::Array4<amrex::Real const> const& values_on_eb_inflow,
#endif
                              bool godunov_use_ppm, bool godunov_use_forces_in_trans,
                              bool is_velocity, bool fluxes_are_area_weighted,
                              std::string const& advection_type,
                              int limiter_type = PPM::default_limiter,
                              bool allow_inflow_on_outflow = false,
                              amrex::Array4<int const> const& bc_arr = {});
/**
 * \brief Compute edge state and flux. Allows inflow on EB but NOT time-dependent inflow.
 *
 *
 * \note See ComputeFluxes for the meaning of iconserv, in particular the
 *       iconserv[n] == -1 convention, which means the face array already holds
 *       the final flux and bypasses every geometric factor. It requires
 *       knownFaceState = true.
 *
 * \note Ghost-cell requirements, when knownFaceState is false and the edge state
 *       is predicted here: q needs 3 filled ghost cells (4 for 3D EBGodunov, and
 *       BDS needs the corners as well), the MAC velocities need 1, and divu needs
 *       1. divu must be a valid Array4 in 3D Godunov for every component, and for
 *       BDS whenever any component is conservative.
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
#if defined(AMREX_USE_EB) && !defined(HYDRO_NO_EB)
                              const amrex::EBFArrayBoxFactory& ebfact,
                              amrex::Array4<amrex::Real const> const& values_on_eb_inflow,
#endif
                              bool godunov_use_ppm, bool godunov_use_forces_in_trans,
                              bool is_velocity, bool fluxes_are_area_weighted,
                              std::string const& advection_type,
                              int limiter_type = PPM::default_limiter,
                              bool allow_inflow_on_outflow = false,
                              amrex::Array4<int const> const& bc_arr = {});

#if defined(AMREX_USE_EB) && !defined(HYDRO_NO_EB)
/**
 * \brief Compute edge state and flux. Allows time-depdendent inflow but no inflow through EB.
 *
 *
 * \note See ComputeFluxes for the meaning of iconserv, in particular the
 *       iconserv[n] == -1 convention, which means the face array already holds
 *       the final flux and bypasses every geometric factor. It requires
 *       knownFaceState = true.
 *
 * \note Ghost-cell requirements, when knownFaceState is false and the edge state
 *       is predicted here: q needs 3 filled ghost cells (4 for 3D EBGodunov, and
 *       BDS needs the corners as well), the MAC velocities need 1, and divu needs
 *       1. divu must be a valid Array4 in 3D Godunov for every component, and for
 *       BDS whenever any component is conservative.
 */
void
ComputeFluxesOnBoxFromState (amrex::Box const& bx, int ncomp, amrex::MFIter& mfi,
                             amrex::Array4<amrex::Real const> const& q,
                             amrex::Array4<amrex::Real const> const& qnph,
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
                             /*amrex::Array4<amrex::Real const> const& values_on_eb_inflow,*/
                             bool godunov_use_ppm, bool godunov_use_forces_in_trans,
                             bool is_velocity, bool fluxes_are_area_weighted,
                             std::string const& advection_type,
                             int limiter_type = PPM::default_limiter,
                             bool allow_inflow_on_outflow = false,
                             amrex::Array4<int const> const& bc_arr = {});
#endif

/**
 * \brief Compute edge state and flux. Does NOT allow time-depdendent inflow or inflow through EB.
 *
 *
 * \note See ComputeFluxes for the meaning of iconserv, in particular the
 *       iconserv[n] == -1 convention, which means the face array already holds
 *       the final flux and bypasses every geometric factor. It requires
 *       knownFaceState = true.
 *
 * \note Ghost-cell requirements, when knownFaceState is false and the edge state
 *       is predicted here: q needs 3 filled ghost cells (4 for 3D EBGodunov, and
 *       BDS needs the corners as well), the MAC velocities need 1, and divu needs
 *       1. divu must be a valid Array4 in 3D Godunov for every component, and for
 *       BDS whenever any component is conservative.
 */
#if defined(AMREX_USE_EB) && !defined(HYDRO_NO_EB)
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
                             /*amrex::Array4<amrex::Real const> const& values_on_eb_inflow,*/
                             bool godunov_use_ppm, bool godunov_use_forces_in_trans,
                             bool is_velocity, bool fluxes_are_area_weighted,
                             std::string const& advection_type,
                             int limiter_type = PPM::default_limiter,
                             bool allow_inflow_on_outflow = false,
                             amrex::Array4<int const> const& bc_arr = {});
#endif

/**
 * \brief Compute edge state and flux. Allows time-depdendent inflow and inflow through EB.
 *
 *
 * \note See ComputeFluxes for the meaning of iconserv, in particular the
 *       iconserv[n] == -1 convention, which means the face array already holds
 *       the final flux and bypasses every geometric factor. It requires
 *       knownFaceState = true.
 *
 * \note Ghost-cell requirements, when knownFaceState is false and the edge state
 *       is predicted here: q needs 3 filled ghost cells (4 for 3D EBGodunov, and
 *       BDS needs the corners as well), the MAC velocities need 1, and divu needs
 *       1. divu must be a valid Array4 in 3D Godunov for every component, and for
 *       BDS whenever any component is conservative.
 */
void
ComputeFluxesOnBoxFromState (amrex::Box const& bx, int ncomp, amrex::MFIter& mfi,
                             amrex::Array4<amrex::Real const> const& q,
                             amrex::Array4<amrex::Real const> const& qnph,
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
#if defined(AMREX_USE_EB) && !defined(HYDRO_NO_EB)
                             const amrex::EBFArrayBoxFactory& ebfact,
                             amrex::Array4<amrex::Real const> const& values_on_eb_inflow,
#endif
                             bool godunov_use_ppm, bool godunov_use_forces_in_trans,
                             bool is_velocity, bool fluxes_are_area_weighted,
                             std::string const& advection_type,
                             int limiter_type = PPM::default_limiter,
                             bool allow_inflow_on_outflow = false,
                             amrex::Array4<int const> const& bc_arr = {});

#if defined(AMREX_USE_EB) && !defined(HYDRO_NO_EB)
/**
 * \brief Extrapolate the velocity to faces, for use as the pre-MAC velocity.
 *
 * \note Ghost-cell requirements: vel needs 3 filled ghost cells, or 4 for
 *       3D EBGodunov; vel_forces needs 1; and velocity_on_eb_inflow, when given,
 *       needs 2.
 */
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
                   int limiter_type = PPM::default_limiter,
                   bool allow_inflow_on_outflow = false);
#endif

/**
 * \brief Extrapolate the velocity to faces, for use as the pre-MAC velocity.
 *
 * \note Ghost-cell requirements: vel needs 3 filled ghost cells, or 4 for
 *       3D EBGodunov; vel_forces needs 1; and velocity_on_eb_inflow, when given,
 *       needs 2.
 */
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
#if defined(AMREX_USE_EB) && !defined(HYDRO_NO_EB)
                   const amrex::EBFArrayBoxFactory& ebfact,
                   amrex::MultiFab const* velocity_on_eb_inflow,
#endif
                   bool godunov_ppm, bool godunov_use_forces_in_trans,
                   std::string const& advection_type,
                   int limiter_type = PPM::default_limiter,
                   bool allow_inflow_on_outflow = false,
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
#if defined(AMREX_USE_EB) && !defined(HYDRO_NO_EB)
                        const amrex::EBFArrayBoxFactory& ebfact,
#endif
                        std::string const& advection_type);

/**
 * \brief Compute Fluxes.
 *
 * \param[in] iconserv  Array of length ncomp giving, per component:
 *                     1  conservative,
 *                     0  convective,
 *                    -1  the face array already holds the final flux.
 *                    For -1 the face value is copied to the flux unchanged: no
 *                    velocity, no face area, no RZ \f$2\pi r\f$ and, in the EB
 *                    version, no aperture factor is applied. A caller that uses -1
 *                    together with fluxes_are_area_weighted must therefore supply a
 *                    face array that is already weighted by all of those factors.
 *                    -1 is only meaningful with knownFaceState = true; with
 *                    knownFaceState = false the edge-state kernels treat -1 as
 *                    conservative and write a predicted state into the face array,
 *                    which is then copied out as the flux with no velocity factor.
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
                     bool fluxes_are_area_weighted,
                     int const* iconserv);

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

#if defined(AMREX_USE_EB) && !defined(HYDRO_NO_EB)

/**
 * \brief Compute Fluxes in the presence of EB.
 *
 * \param[in] iconserv  Array of length ncomp giving, per component:
 *                     1  conservative,
 *                     0  convective,
 *                    -1  the face array already holds the final flux.
 *                    For -1 the face value is copied to the flux unchanged: no
 *                    velocity, no face area, no RZ \f$2\pi r\f$ and, in the EB
 *                    version, no aperture factor is applied. A caller that uses -1
 *                    together with fluxes_are_area_weighted must therefore supply a
 *                    face array that is already weighted by all of those factors.
 *                    -1 is only meaningful with knownFaceState = true; with
 *                    knownFaceState = false the edge-state kernels treat -1 as
 *                    conservative and write a predicted state into the face array,
 *                    which is then copied out as the flux with no velocity factor.
 *
 * \note EB advection requires Cartesian geometry; the routine asserts on this.
 */
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
                        bool fluxes_are_area_weighted,
                        int const* iconserv);


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

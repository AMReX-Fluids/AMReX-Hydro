/**
 * \file hydro_mol.H
 * \brief A documented file.
 *
 * Here is where detailed documentation for the file can go.
 *
 */

/** \addtogroup MOL
 *  @{
 */





#ifndef HYDRO_MOL_H_
#define HYDRO_MOL_H_

#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>


/**
 * \namespace MOL
 * The procedures for computing MAC velocities and edge states with MOL does
 * not involve time derivatives. All slope computations use second-order
 * limited slopes as described in
 * <a href="https://amrex-codes.github.io/amrex/hydro_html/MOL.html#sec:slopes">[sec:slopes]</a>.
 *
 */


namespace MOL {

/**
 *
 * \brief Compute upwinded face-centered states by extrapolating cell-centered values in SPACE ONLY.
 * \param [in]     a_bx        Current grid patch
 * \param [in,out] xedge       Array4 containing x-edges, starting at component of interest
 * \param [in,out] yedge       Array4 containing y-edges, starting at component of interest
 * \param [in,out] zedge       Array4 containing z-edges, starting at component of interest
 * \param [in]     q           Array4 of state, starting at component of interest
 * \param [in]     ncomp       Number of components to work on
 * \param [in]     umac        x-Face velocities.
 * \param [in]     vmac        y-Face velocities.
 * \param [in]     wmac        z-Face velocities.
 * \param [in]     domain      Problem domain.
 * \param [in]     bcs         Vector of boundary conditions
 * \param [in]     d_bcrec_ptr BCRec
 * \param [in]     is_velocity Indicates a component is velocity so boundary conditions can
 *                             be properly addressed. The header hydro_constants.H
 *                             defines the component positon by [XYZ]VEL macro.
 *
 * Compute upwinded FC states by extrapolating CC values in SPACE ONLY.
 * This is NOT a Godunov type extrapolation: there is NO dependence on time!
 * The resulting FC states are computed at the CENTROID of the face.
 *
 */

void ComputeEdgeState (amrex::Box const& a_bx,
                       AMREX_D_DECL( amrex::Array4<amrex::Real> const& xedge,
                                     amrex::Array4<amrex::Real> const& yedge,
                                     amrex::Array4<amrex::Real> const& zedge),
                       amrex::Array4<amrex::Real const> const& q,
                       const int ncomp,
                       AMREX_D_DECL( amrex::Array4<amrex::Real const> const& umac,
                                     amrex::Array4<amrex::Real const> const& vmac,
                                     amrex::Array4<amrex::Real const> const& wmac),
                       amrex::Box const&       domain,
                       amrex::Vector<amrex::BCRec> const& bcs,
                       amrex::BCRec  const* d_bcrec_ptr,
                       bool is_velocity);


/**
 *
 * \brief Compute upwinded FC velocities by extrapolating CC values in SPACE ONLY.
 * \param vel   The cell-centered velocity field
 * \param umac  x-component of velocity on centroid of x-face
 * \param vmac  y-component of velocity on centroid of y-face
 * \param wmac  z-component of velocity on centroid of z-face
 * \param a_geom  Geometry object at this level
 * \param h_bcrec Host version of BCRec
 * \param d_bcrec Device version of BCRec
 *
 * Compute upwinded FC velocities by extrapolating CC values in SPACE ONLY.
 * This is NOT a Godunov type extrapolation: there is NO dependence on time!
 * The resulting FC velocities are computed at the CENTROID of the face.
 *
 */


void ExtrapVelToFaces ( const amrex::MultiFab&  vel,
                        AMREX_D_DECL(amrex::MultiFab& umac,
                                     amrex::MultiFab& vmac,
                                     amrex::MultiFab& wmac ),
                        const amrex::Geometry&  a_geom,
                        amrex::Vector<amrex::BCRec> const& h_bcrec,
                        const amrex::BCRec* d_bcrec );

/**
 * \brief For Computing the pre-MAC edge states to be MAC-projected.
 *
 * For detailed info on the algorithm, see
 * <a href="https://amrex-codes.github.io/amrex/hydro_html/MOL.html#pre-mac-mol-extrapveltofacesbox"> MOL</a>.
 *
 */
void ExtrapVelToFacesBox ( AMREX_D_DECL( amrex::Box const& ubx,
                                         amrex::Box const& vbx,
                                         amrex::Box const& wbx ),
                           AMREX_D_DECL( amrex::Array4<amrex::Real> const& u,
                                         amrex::Array4<amrex::Real> const& v,
                                         amrex::Array4<amrex::Real> const& w ),
                           amrex::Array4<amrex::Real const > const& vcc,
                           const amrex::Geometry&  geom,
                           amrex::Vector<amrex::BCRec> const& h_bcrec,
                           const amrex::BCRec* d_bcrec );


}
#endif
/** @}*/

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
 *  \brief Example documentation
 *
 *  \param aofs My favorite variable.
 *  \param aofs_comp My second variable.
 *
 *  Doxygen demo docs.
 *
 */


void ComputeAofs ( amrex::MultiFab& aofs, int aofs_comp, int ncomp,
                   amrex::MultiFab const& state, int state_comp,
                   AMREX_D_DECL( amrex::MultiFab const& umac,
                                 amrex::MultiFab const& vmac,
                                 amrex::MultiFab const& wmac),
                   AMREX_D_DECL( amrex::MultiFab& xedge,
                                 amrex::MultiFab& yedge,
                                 amrex::MultiFab& zedge),
                   int  edge_comp,
                   bool known_edgestate,
                   AMREX_D_DECL( amrex::MultiFab& xfluxes,
                                 amrex::MultiFab& yfluxes,
                                 amrex::MultiFab& zfluxes),
                   int fluxes_comp,
                   amrex::MultiFab const& divu,
                   amrex::Vector<amrex::BCRec> const& bcs,
                   amrex::BCRec  const* d_bcrec_ptr,
                   amrex::Gpu::DeviceVector<int>& iconserv,
                   amrex::Geometry const& geom,
                   bool is_velocity);

/**
 *  <A ID="ComputeSyncAofs"></A> 
 *  
 *  \brief Example documentation
 *
 *  \param[in] aofs This is an `in` variable.
 *  \param[in,out] aofs_comp This is an `in,out` variable.
 *
 *  This is an example of the documentation capabilities.  
 *
 *  The html above 
 *  creates an anchor that can be used to link directly to the documentation
 *  for this function. Navigating to
 *
 *      https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceMOL.html#ComputeSyncAofs
 *
 *  will bring you to the documentation for this function specifically.
 *
 *  The tag below allows for cross referencing of other functions in Doxygen.
 *
 *  \ref ExtrapVelToFaces
 *
 *
 */


void ComputeSyncAofs ( amrex::MultiFab& aofs, int aofs_comp, int ncomp,
                       amrex::MultiFab const& state, int state_comp,
                       AMREX_D_DECL( amrex::MultiFab const& umac,
                                     amrex::MultiFab const& vmac,
                                     amrex::MultiFab const& wmac),
                       AMREX_D_DECL( amrex::MultiFab const& ucorr,
                                     amrex::MultiFab const& vcorr,
                                     amrex::MultiFab const& wcorr),
                       AMREX_D_DECL( amrex::MultiFab& xedge,
                                     amrex::MultiFab& yedge,
                                     amrex::MultiFab& zedge),
                       int  edge_comp,
                       bool known_edgestate,
                       AMREX_D_DECL( amrex::MultiFab& xfluxes,
                                     amrex::MultiFab& yfluxes,
                                     amrex::MultiFab& zfluxes),
                       int fluxes_comp,
                       amrex::Vector<amrex::BCRec> const& bcs,
                       amrex::BCRec  const* d_bcrec_ptr,
                       amrex::Geometry const& geom,
                       bool is_velocity);

/**
 *  \brief ComputeEdgeState description.
 *  
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

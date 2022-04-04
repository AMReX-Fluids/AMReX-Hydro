/**
 * \file hydro_bds.H
 *
 * \addtogroup BDS
 *  @{
 */

#ifndef HYDRO_BDS_H
#define HYDRO_BDS_H

#include <AMReX_MultiFabUtil.H>
#include <AMReX_BCRec.H>

/**
 * Collection of routines for the BDS (Bell-Dawson-Shubin) algorithm.
 */
namespace BDS {

/**
 * Compute advection of a scalar (s).
 *
 * \param [in]     aofs             Multifab containing the scalar field to advect.
 * \param [in]     aofs_comp        Integer representing location of first component to advect.
 * \param [in]     ncomp            Number of components to advect.
 * \param [in]     state            Multifab of state.
 * \param [in]     state_comp       The location of the component of interest in the Multifab of state.
 * \param [in]     umac             x-Face velocities.
 * \param [in]     vmac             y-Face velocities.
 * \param [in]     wmac             z-Face velocities.
 * \param [in,out] xedge            x-Edge velocities.
 * \param [in,out] yedge            y-Edge velocities.
 * \param [in,out] zedge            z-Edge velocities.
 * \param [in]     edge_comp        Integer location of the edge component of the edge MultiFabs.
 * \param [in]     known_edgestate  Indicates whether the edge state is known.
 * \param [in,out] xfluxes          x-fluxes.
 * \param [in,out] yfluxes          y-fluxes.
 * \param [in,out] zfluxes          z-fluxes.
 * \param [in]     fluxes_comp      Integer location of the flux component of the flux MultiFabs.
 * \param [in]     fq               MultiFab for forces.
 * \param [in]     fq_comp          Integer location of the force component of the forces MultiFab.
 * \param [in,out] divu             MultiFab for the divergence of the cell-centered velocity field.
 * \param [in]     d_bc             Boundary conditions.
 * \param [in]     geom             Level geometry.
 * \param [in]     iconserv         Indicates conservative dimensions.
 * \param [in]     dt               Time step.
 * \param [in]     is_velocity      Indicates a component is velocity so boundary conditions can
 *                                  be properly addressed. The header hydro_constants.H
 *                                  defines the component positon by [XYZ]VEL macro.
 */

void ComputeAofs ( amrex::MultiFab& aofs, const int aofs_comp, const int ncomp,
                   amrex::MultiFab const& state, const int state_comp,
                   AMREX_D_DECL( amrex::MultiFab const& umac,
                                 amrex::MultiFab const& vmac,
                                 amrex::MultiFab const& wmac),
                   AMREX_D_DECL( amrex::MultiFab& xedge,
                                 amrex::MultiFab& yedge,
                                 amrex::MultiFab& zedge),
                   const int  edge_comp,
                   const bool known_edgestate,
                   AMREX_D_DECL( amrex::MultiFab& xfluxes,
                                 amrex::MultiFab& yfluxes,
                                 amrex::MultiFab& zfluxes),
                   int fluxes_comp,
                   amrex::MultiFab const& fq,
                   const int fq_comp,
                   amrex::MultiFab const& divu,
                   amrex::BCRec const* d_bc,
                   amrex::Geometry const& geom,
                   amrex::Vector<int>& iconserv,
                   const amrex::Real dt,
                   const bool is_velocity);
/**
 * Synchronize the advection of a scalar (s) across levels.
 *
 * \param [in]     aofs             Multifab containing the scalar field to advect.
 * \param [in]     aofs_comp        Integer representing location of first component to advect.
 * \param [in]     ncomp            Number of components to advect.
 * \param [in]     state            Multifab of state.
 * \param [in]     state_comp       The location of the component of interest in the Multifab of state.
 * \param [in]     umac             x-Face velocities.
 * \param [in]     vmac             y-Face velocities.
 * \param [in]     wmac             z-Face velocities.
 * \param [in]     ucorr            x-Face velocity corrections.
 * \param [in]     vcorr            y-Face velocity corrections.
 * \param [in]     wcorr            z-Face velocity corrections.
 * \param [in,out] xedge            x-Edge velocities.
 * \param [in,out] yedge            y-Edge velocities.
 * \param [in,out] zedge            z-Edge velocities.
 * \param [in]     edge_comp        Integer location of the edge component of the edge MultiFabs.
 * \param [in]     known_edgestate  Indicates whether the edge state is known.
 * \param [in,out] xfluxes          x-fluxes.
 * \param [in,out] yfluxes          y-fluxes.
 * \param [in,out] zfluxes          z-fluxes.
 * \param [in]     fluxes_comp      Integer location of the flux component of the flux MultiFabs.
 * \param [in]     fq               MultiFab for forces.
 * \param [in]     fq_comp          Integer location of the force component of the forces MultiFab.
 * \param [in,out] divu             MultiFab for the divergence of the cell-centered velocity field.
 * \param [in]     d_bc             Boundary conditions.
 * \param [in]     geom             Level geometry.
 * \param [in]     iconserv         Indicates conservative dimensions.
 * \param [in]     dt               Time step.
 * \param [in]     is_velocity      Indicates a component is velocity so boundary conditions can
 *                                  be properly addressed. The header hydro_constants.H
 *                                  defines the component positon by [XYZ]VEL macro.
 */

void ComputeSyncAofs ( amrex::MultiFab& aofs, const int aofs_comp, const int ncomp,
                       amrex::MultiFab const& state, const int state_comp,
                       AMREX_D_DECL( amrex::MultiFab const& umac,
                                     amrex::MultiFab const& vmac,
                                     amrex::MultiFab const& wmac),
                       AMREX_D_DECL( amrex::MultiFab const& ucorr,
                                     amrex::MultiFab const& vcorr,
                                     amrex::MultiFab const& wcorr),
                       AMREX_D_DECL( amrex::MultiFab& xedge,
                                     amrex::MultiFab& yedge,
                                     amrex::MultiFab& zedge),
                       const int  edge_comp,
                       const bool known_edgestate,
                       AMREX_D_DECL( amrex::MultiFab& xfluxes,
                                     amrex::MultiFab& yfluxes,
                                     amrex::MultiFab& zfluxes),
                       int fluxes_comp,
                       amrex::MultiFab const& fq,
                       const int fq_comp,
                       amrex::MultiFab const& divu,
                       amrex::BCRec const* d_bc,
                       amrex::Geometry const& geom,
                       amrex::Gpu::DeviceVector<int>& iconserv,
                       const amrex::Real dt,
                       const bool is_velocity);

/**
 * Uses the Bell-Dawson-Shubin (BDS) algorithm, a higher order Godunov
 * method for scalar conservation laws in three dimensions, to compute
 * edge states.
 *
 * \param [in]     bx          Current grid patch
 * \param [in]     ncomp       Number of components to work on
 * \param [in]     q           Array4 of state, starting at component of interest
 * \param [in,out] xedge       Array4 containing x-edges, starting at component of interest
 * \param [in,out] yedge       Array4 containing y-edges, starting at component of interest
 * \param [in,out] zedge       Array4 containing z-edges, starting at component of interest
 * \param [in]     umac        x-Face velocities.
 * \param [in]     vmac        y-Face velocities.
 * \param [in]     wmac        z-Face velocities.
 * \param [in]     fq          Array4 for forces, starting at component of interest
 * \param [in]     geom        Level geometry.
 * \param [in]     l_dt        Time step.
 * \param [in]     iconserv    Indicates conservative dimensions.
 * \param [in]     is_velocity Indicates a component is velocity so boundary conditions can
 *                             be properly addressed. The header hydro_constants.H
 *                             defines the component positon by [XYZ]VEL macro.
 */

void ComputeEdgeState ( amrex::Box const& bx, int ncomp,
                        amrex::Array4<amrex::Real const> const& q,
                        AMREX_D_DECL(amrex::Array4<amrex::Real> const& xedge,
                                     amrex::Array4<amrex::Real> const& yedge,
                                     amrex::Array4<amrex::Real> const& zedge),
                        AMREX_D_DECL(amrex::Array4<amrex::Real const> const& umac,
                                     amrex::Array4<amrex::Real const> const& vmac,
                                     amrex::Array4<amrex::Real const> const& wmac),
                        amrex::Array4<amrex::Real const> const& divu,
                        amrex::Array4<amrex::Real const> const& fq,
                        amrex::Geometry geom,
                        amrex::Real l_dt,
                        amrex::BCRec const* pbc,
                        int const* iconserv,
                        const bool is_velocity);

/**
 * Compute bilinear slopes for BDS algorithm.
 *
 * \param [in]  bx      Current grid patch
 * \param [in]  geom    Level geometry.
 * \param [in]  icomp   Component of the state Array4.
 * \param [in]  s       Array4<const> of state vector.
 * \param [out] slopes  Array4 to store slope information.
 *
 */

void ComputeSlopes ( amrex::Box const& bx,
                     const amrex::Geometry& geom,
                     int icomp,
                     amrex::Array4<amrex::Real const> const& s,
                     amrex::Array4<amrex::Real      > const& slopes,
                     amrex::BCRec const* pbc);

/**
 * Compute Conc for BDS algorithm.
 *
 * \param [in]     bx          Current grid patch
 * \param [in]     geom        Level geometry.
 * \param [in]     icomp       Component of the Array4s.
 * \param [in]     s           Array4 of state.
 * \param [in,out] sedgex      Array4 containing x-edges.
 * \param [in,out] sedgey      Array4 containing y-edges.
 * \param [in,out] sedgez      Array4 containing z-edges.
 * \param [in]     slopes      Array4 containing slope information.
 * \param [in]     umac        Array4 for u-face velocity.
 * \param [in]     vmac        Array4 for v-face velocity.
 * \param [in]     wmac        Array4 for z-face velocity.
 * \param [in]     force       Array4 for forces.
 * \param [in]     iconserv    Indicates conservative dimensions.
 * \param [in]     dt          Time step.
 * \param [in]     is_velocity Indicates a component is velocity so boundary conditions can
 *                             be properly addressed. The header hydro_constants.H
 *                             defines the component positon by [XYZ]VEL macro.
 *
 *
 */

void ComputeConc ( amrex::Box const& bx,
                   const amrex::Geometry& geom,
                   int icomp,
                   amrex::Array4<amrex::Real const> const& s,
                   AMREX_D_DECL(amrex::Array4<amrex::Real      > const& sedgex,
                                amrex::Array4<amrex::Real      > const& sedgey,
                                amrex::Array4<amrex::Real      > const& sedgez),
                   amrex::Array4<amrex::Real const> const& slopes,
                   AMREX_D_DECL(amrex::Array4<amrex::Real const> const& umac,
                                amrex::Array4<amrex::Real const> const& vmac,
                                amrex::Array4<amrex::Real const> const& wmac),
                   amrex::Array4<amrex::Real const> const& divu,
                   amrex::Array4<amrex::Real const> const& force,
                   int const* iconserv,
                   const amrex::Real dt,
                   amrex::BCRec const* pbc,
                   const bool is_velocity);

}
#endif
/** @} */

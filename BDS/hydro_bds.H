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
 * \param [in]     iconserv    If true, use conservative form, otherwise use convective.
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
 * \param [in]     iconserv    If true, use conservative form, otherwise use convective.
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

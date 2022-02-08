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

namespace BDS {


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
                   amrex::Geometry const& geom,
                   amrex::Vector<int>& iconserv,
                   const amrex::Real dt);

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
                       amrex::Geometry const& geom,
                       amrex::Gpu::DeviceVector<int>& iconserv,
                       const amrex::Real dt);

void ComputeEdgeState ( amrex::Box const& bx, int ncomp,
                        amrex::Array4<amrex::Real const> const& q,
                        AMREX_D_DECL(amrex::Array4<amrex::Real> const& xedge,
                                     amrex::Array4<amrex::Real> const& yedge,
                                     amrex::Array4<amrex::Real> const& zedge),
                        AMREX_D_DECL(amrex::Array4<amrex::Real const> const& umac,
                                     amrex::Array4<amrex::Real const> const& vmac,
                                     amrex::Array4<amrex::Real const> const& wmac),
                        amrex::Array4<amrex::Real const> const& fq,
                        amrex::Geometry geom,
                        amrex::Real l_dt,
                        /*BCRec const* pbc,*/
                        int const* iconserv);

void ComputeSlopes ( amrex::Box const& bx,
                     const amrex::Geometry& geom,
                     int icomp,
                     amrex::Array4<amrex::Real const> const& s,
                     amrex::Array4<amrex::Real      > const& slopes);

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
                   amrex::Array4<amrex::Real const> const& force,
                   int const* iconserv,
                   const amrex::Real dt);

}
#endif
/** @} */

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

void ComputeEdgeState (amrex::MultiFab const& state,
                          const int state_comp,
                          amrex::Geometry const& geom,
                          AMREX_D_DECL(amrex::MultiFab& xedge,
                                       amrex::MultiFab& yedge,
                                       amrex::MultiFab& zedge),
                          const int edge_comp,
                          AMREX_D_DECL(amrex::MultiFab const& umac,
                                       amrex::MultiFab const& vmac,
                                       amrex::MultiFab const& wmac),
                          amrex::MultiFab const& fq,
                          const int fq_comp,
                          const int is_conservative,
                          const amrex::Real dt);

void ComputeSlopes (amrex::MultiFab const& s_mf,
                    const amrex::Geometry& geom,
                    amrex::MultiFab& slope_mf,
                    const int state_comp);

void ComputeConc (const amrex::MultiFab& s_mf,
                  const int state_comp,
                  const amrex::Geometry& geom,
                  AMREX_D_DECL(amrex::MultiFab& xedge,
                               amrex::MultiFab& yedge,
                               amrex::MultiFab& zedge),
                  const int edge_comp,
                  const amrex::MultiFab& slope_mf,
                  AMREX_D_DECL(amrex::MultiFab const& umac,
                               amrex::MultiFab const& vmac,
                               amrex::MultiFab const& wmac),
                  amrex::MultiFab const& fq,
                  const int fq_comp,
                  const int is_conservative,
                  const amrex::Real dt);

}
#endif
/** @} */

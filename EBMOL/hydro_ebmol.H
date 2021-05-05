#ifndef HYDRO_EBMOL_H_
#define HYDRO_EBMOL_H_

#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>

namespace EBMOL {

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
                   amrex::Real dt,
                   std::string redistribution_type );

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
                       amrex::Real dt,
                       std::string redistribution_type );

void ComputeEdgeState ( amrex::Box const& bx,
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
                        AMREX_D_DECL( amrex::Array4<amrex::Real const> const& fcx,
                                      amrex::Array4<amrex::Real const> const& fcy,
                                      amrex::Array4<amrex::Real const> const& fcz),
                        amrex::Array4<amrex::Real const> const& ccc,
                        amrex::Array4<amrex::Real const> const& vfrac,
                        amrex::Array4<amrex::EBCellFlag const> const& flag);

void ExtrapVelToFaces ( const amrex::MultiFab&  vel,
                        AMREX_D_DECL(amrex::MultiFab& umac,
                                     amrex::MultiFab& vmac,
                                     amrex::MultiFab& wmac ),
                        const amrex::Geometry&  a_geom,
                        amrex::Vector<amrex::BCRec> const& h_bcrec,
                        const amrex::BCRec* d_bcrec );

void ExtrapVelToFacesBox( AMREX_D_DECL( amrex::Box const& ubx,
                                        amrex::Box const& vbx,
                                        amrex::Box const& wbx ),
                          AMREX_D_DECL( amrex::Array4<amrex::Real> const& u,
                                        amrex::Array4<amrex::Real> const& v,
                                        amrex::Array4<amrex::Real> const& w ),
                          amrex::Array4<amrex::Real const > const& vcc,
                          amrex::Array4<amrex::EBCellFlag const> const& flag,
                          AMREX_D_DECL( amrex::Array4<amrex::Real const> const& fcx,
                                        amrex::Array4<amrex::Real const> const& fcy,
                                        amrex::Array4<amrex::Real const> const& fcz),
                          amrex::Array4<amrex::Real const> const& ccc,
                          amrex::Array4<amrex::Real const> const& vfrac,
                          const amrex::Geometry&  geom,
                          amrex::Vector<amrex::BCRec> const& h_bcrec,
                          const amrex::BCRec* d_bcrec );

}

#endif

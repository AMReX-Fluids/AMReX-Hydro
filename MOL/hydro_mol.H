#ifndef HYDRO_MOL_H_
#define HYDRO_MOL_H_

#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>

namespace MOL {

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
                   amrex::Vector<amrex::BCRec> const& bcs,
                   amrex::BCRec  const* d_bcrec_ptr,
                   amrex::Gpu::DeviceVector<int>& iconserv,
                   amrex::Geometry const& geom,
                   const amrex::Real dt,
                   bool is_velocity);

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
                       const amrex::Real dt,
                       bool is_velocity);

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


void ExtrapVelToFaces ( const amrex::MultiFab&  vel,
                        AMREX_D_DECL(amrex::MultiFab& umac,
                                     amrex::MultiFab& vmac,
                                     amrex::MultiFab& wmac ),
                        const amrex::Geometry&  a_geom,
                        amrex::Vector<amrex::BCRec> const& h_bcrec,
                        const amrex::BCRec* d_bcrec );

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

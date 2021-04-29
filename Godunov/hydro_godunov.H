#ifndef HYDRO_GODUNOV_H
#define HYDRO_GODUNOV_H

#include <AMReX_MultiFabUtil.H>
#include <AMReX_BCRec.H>

namespace Godunov {


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
                   amrex::Gpu::DeviceVector<int>& iconserv,
                   const amrex::Real dt,
                   const bool use_ppm,
                   const bool use_forces_in_trans,
                   const bool is_velocity  );

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
                       const bool use_ppm,
                       const bool use_forces_in_trans,
                       const bool is_velocity  );


    // void ComputeFluxes ( amrex::Box const& bx,
    //                      AMREX_D_DECL( amrex::Array4<amrex::Real> const& fx,
    //                                    amrex::Array4<amrex::Real> const& fy,
    //                                    amrex::Array4<amrex::Real> const& fz),
    //                      AMREX_D_DECL( amrex::Array4<amrex::Real const> const& umac,
    //                                    amrex::Array4<amrex::Real const> const& vmac,
    //                                    amrex::Array4<amrex::Real const> const& wmac),
    //                      AMREX_D_DECL( amrex::Array4<amrex::Real const> const& xedge,
    //                                    amrex::Array4<amrex::Real const> const& yedge,
    //                                    amrex::Array4<amrex::Real const> const& zedge),
    //                      amrex::Geometry const& geom, const int ncomp );


    // void ComputeDivergence ( amrex::Box const& bx,
    //                          amrex::Array4<amrex::Real> const& div,
    //                          AMREX_D_DECL( amrex::Array4<amrex::Real const> const& fx,
    //                                        amrex::Array4<amrex::Real const> const& fy,
    //                                        amrex::Array4<amrex::Real const> const& fz),
    //                          AMREX_D_DECL( amrex::Array4<amrex::Real const> const& xed,
    //                                        amrex::Array4<amrex::Real const> const& yed,
    //                                        amrex::Array4<amrex::Real const> const& zed),
    //                          AMREX_D_DECL( amrex::Array4<amrex::Real const> const& umac,
    //                                        amrex::Array4<amrex::Real const> const& vmac,
    //                                        amrex::Array4<amrex::Real const> const& wmac),
    //                          const int ncomp, amrex::Geometry const& geom,
    //                          int const* iconserv );

    // void ComputeSyncDivergence ( amrex::Box const& bx,
    //                              amrex::Array4<amrex::Real> const& div,
    //                              AMREX_D_DECL( amrex::Array4<amrex::Real const> const& fx,
    //                                            amrex::Array4<amrex::Real const> const& fy,
    //                                            amrex::Array4<amrex::Real const> const& fz),
    //                              const int ncomp, amrex::Geometry const& geom );



    // void ComputeSyncDivergence_rz ( amrex::Box const& bx,
    //     			    amrex::Array4<amrex::Real> const& div,
    //     			    amrex::Array4<amrex::Real const> const& fx,
    //     			    amrex::Array4<amrex::Real const> const& fy,
    //     			    amrex::Array4<amrex::Real const> const& vol,
    //     			    const int ncomp );

void ExtrapVelToFaces ( amrex::MultiFab const& a_vel,
                        amrex::MultiFab const& a_forces,
                        AMREX_D_DECL( amrex::MultiFab& a_umac,
                                      amrex::MultiFab& a_vmac,
                                      amrex::MultiFab& a_wmac ),
                        const amrex::Vector<amrex::BCRec> & h_bcrec,
                        const               amrex::BCRec  * d_bcrec,
                        const amrex::Geometry& geom, amrex::Real l_dt,
                        bool use_ppm, bool use_forces_in_trans);

void ComputeAdvectiveVel (AMREX_D_DECL(amrex::Box const& xbx,
                                       amrex::Box const& ybx,
                                       amrex::Box const& zbx),
                          AMREX_D_DECL(amrex::Array4<amrex::Real> const& u_ad,
                                       amrex::Array4<amrex::Real> const& v_ad,
                                       amrex::Array4<amrex::Real> const& w_ad),
                          AMREX_D_DECL(amrex::Array4<amrex::Real const> const& Imx,
                                       amrex::Array4<amrex::Real const> const& Imy,
                                       amrex::Array4<amrex::Real const> const& Imz),
                          AMREX_D_DECL(amrex::Array4<amrex::Real const> const& Ipx,
                                       amrex::Array4<amrex::Real const> const& Ipy,
                                       amrex::Array4<amrex::Real const> const& Ipz),
                          amrex::Array4<amrex::Real const> const& vel,
                          amrex::Array4<amrex::Real const> const& f,
                          const amrex::Box& domain,
                          amrex::Real dt,
                          amrex::BCRec const* d_bcrec,
                          bool use_forces_in_trans);

void ExtrapVelToFacesOnBox (amrex::Box const& bx, int ncomp,
                            AMREX_D_DECL(amrex::Box const& xbx,
                                         amrex::Box const& ybx,
                                         amrex::Box const& zbx),
                            AMREX_D_DECL(amrex::Array4<amrex::Real> const& qx,
                                         amrex::Array4<amrex::Real> const& qy,
                                         amrex::Array4<amrex::Real> const& qz),
                            amrex::Array4<amrex::Real const> const& q,
                            AMREX_D_DECL(amrex::Array4<amrex::Real const> const& u_ad,
                                         amrex::Array4<amrex::Real const> const& v_ad,
                                         amrex::Array4<amrex::Real const> const& w_ad),
                            AMREX_D_DECL(amrex::Array4<amrex::Real> const& Imx,
                                         amrex::Array4<amrex::Real> const& Imy,
                                         amrex::Array4<amrex::Real> const& Imz),
                            AMREX_D_DECL(amrex::Array4<amrex::Real> const& Ipx,
                                         amrex::Array4<amrex::Real> const& Ipy,
                                         amrex::Array4<amrex::Real> const& Ipz),
                            amrex::Array4<amrex::Real const> const& f,
                            const amrex::Box& domain,
                            const amrex::Real* dx,
                            amrex::Real dt,
                            amrex::BCRec const* d_bcrec,
                            bool use_forces_in_trans,
                            amrex::Real* p);

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
                        amrex::Real dt,
                        amrex::BCRec const* d_bcrec,
                        int const* iconserv,
                        const bool use_ppm, bool is_velocity,
                        const bool use_forces_in_trans);

}

#endif
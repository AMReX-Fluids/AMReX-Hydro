/**
 * \file hydro_ebgodunov.H
 * \addtogroup EBGodunov
 * @{
 *
 */

#ifndef HYDRO_EBGODUNOV_H
#define HYDRO_EBGODUNOV_H

#include <AMReX_MultiFabUtil.H>
#include <AMReX_MultiCutFab.H>


namespace EBGodunov {


    void ExtrapVelToFaces ( amrex::MultiFab const& vel,
                            amrex::MultiFab const& vel_forces,
                            AMREX_D_DECL(amrex::MultiFab& u_mac,
                                         amrex::MultiFab& v_mac,
                                         amrex::MultiFab& w_mac),
                            amrex::Vector<amrex::BCRec> const& h_bcrec,
                            amrex::BCRec  const* d_bcrec,
                            const amrex::Geometry& geom,
                            amrex::Real dt, amrex::MultiFab const* velocity_on_eb_inflow = nullptr);


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
                              amrex::Array4<amrex::EBCellFlag const> const& flag,
                              const amrex::Box& domain,
                              amrex::BCRec const* d_bcrec);

    void ExtrapVelToFacesOnBox ( amrex::Box const& bx, int ncomp,
                                 AMREX_D_DECL(amrex::Box const& xbx,
                                              amrex::Box const& ybx,
                                              amrex::Box const& zbx),
                                 AMREX_D_DECL(amrex::Box const& xebx,
                                              amrex::Box const& yebx,
                                              amrex::Box const& zebx),
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
                                 amrex::Array4<amrex::EBCellFlag const> const& flag,
                                 AMREX_D_DECL(amrex::Array4<amrex::Real const> const& apx,
                                              amrex::Array4<amrex::Real const> const& apy,
                                              amrex::Array4<amrex::Real const> const& apz),
#if (AMREX_SPACEDIM == 3)
                                 amrex::Array4<amrex::Real const> const& vfrac,
#endif
                                 AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                                              amrex::Array4<amrex::Real const> const& fcy,
                                              amrex::Array4<amrex::Real const> const& fcz),
                                 amrex::Real* p, amrex::Array4<amrex::Real const> const& velocity_on_eb_inflow);


    void ComputeEdgeState ( amrex::Box const& bx, int ncomp,
                            amrex::Array4<amrex::Real const> const& q,
                            AMREX_D_DECL( amrex::Array4<amrex::Real> const& xedge,
                                          amrex::Array4<amrex::Real> const& yedge,
                                          amrex::Array4<amrex::Real> const& zedge),
                            AMREX_D_DECL( amrex::Array4<amrex::Real const> const& umac,
                                          amrex::Array4<amrex::Real const> const& vmac,
                                          amrex::Array4<amrex::Real const> const& wmac),
                            amrex::Array4<amrex::Real const> const& divu,
                            amrex::Array4<amrex::Real const> const& fq,
                            amrex::Geometry const& geom,
                            amrex::Real l_dt,
                            amrex::Vector<amrex::BCRec> const& h_bcrec,
                            amrex::BCRec const*  pbc,
                            int const* iconserv,
                            amrex::Real* p,
                            amrex::Array4<amrex::EBCellFlag const> const& flag_arr,
                            AMREX_D_DECL(amrex::Array4<amrex::Real const> const& apx,
                                         amrex::Array4<amrex::Real const> const& apy,
                                         amrex::Array4<amrex::Real const> const& apz),
                            amrex::Array4<amrex::Real const> const& vfrac_arr,
                            AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                                         amrex::Array4<amrex::Real const> const& fcy,
                                         amrex::Array4<amrex::Real const> const& fcz),
                            amrex::Array4<amrex::Real const> const& ccent_arr,
                            bool is_velocity,
                            amrex::Array4<amrex::Real const> const& values_on_eb_inflow);

} // namespace ebgodunov

#endif /* EBGodunov_H */
/** @} */

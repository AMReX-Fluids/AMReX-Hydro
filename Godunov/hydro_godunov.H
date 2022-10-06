/**
 * \file hydro_godunov.H
 *
 * \addtogroup Godunov
 *  @{
 */

#ifndef HYDRO_GODUNOV_H
#define HYDRO_GODUNOV_H

#include <AMReX_MultiFabUtil.H>
#include <AMReX_BCRec.H>

namespace Godunov {


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
/** @} */

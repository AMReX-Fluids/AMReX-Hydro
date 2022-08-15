/**
 * \file hydro_ebmol.H
 * \addtogroup EBMOL
 * @{
 *
 */

#ifndef HYDRO_EBMOL_H_
#define HYDRO_EBMOL_H_

#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>

/**
 * \namespace EBMOL
 *
 * The procedure for computing MAC velocities and edge states with EB-aware MOL does not
 * involve any time derivatives. All slope computations use second-order limited slopes
 * as described in
 * <a href="https://amrex-codes.github.io/amrex/hydro_html/EBMOL.html#sec:slopes">[sec:slopes]</a>.
 *
 */
namespace EBMOL {

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
                        amrex::Array4<amrex::EBCellFlag const> const& flag,
                        const bool is_velocity );

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
/** @} */

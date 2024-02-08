/**
 * \file hydro_ebgodunov_plm.H
 * \addtogroup EBGodunov
 * @{
 *
 */

#ifndef HYDRO_PLM_EBGODUNOV_H
#define HYDRO_PLM_EBGODUNOV_H

#include <AMReX_MultiFabUtil.H>
#include <AMReX_MultiCutFab.H>

namespace EBPLM {

void PredictVelOnXFace ( amrex::Box const& xebox,
                         amrex::Array4<amrex::Real> const& Imx,
                         amrex::Array4<amrex::Real> const& Ipx,
                         amrex::Array4<amrex::Real const> const& q,
                         amrex::Array4<amrex::Real  const> const& ccvel,
                         amrex::Array4<amrex::EBCellFlag const> const& flag,
                         amrex::Array4<amrex::Real const> const& vfrac,
                         AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                                      amrex::Array4<amrex::Real const> const& fcy,
                                      amrex::Array4<amrex::Real const> const& fcz),
                         amrex::Array4<amrex::Real const> const& ccc,
                         const amrex::Geometry& geom,
                         amrex::Real dt,
                         amrex::Vector<amrex::BCRec> const& h_bcrec,
                         amrex::BCRec const* pbc,
                         Array4<int const> const& bc_arr = {});

void PredictVelOnYFace ( amrex::Box const& yebox,
                         amrex::Array4<amrex::Real> const& Imy,
                         amrex::Array4<amrex::Real> const& Ipy,
                         amrex::Array4<amrex::Real const> const& q,
                         amrex::Array4<amrex::Real  const> const& ccvel,
                         amrex::Array4<amrex::EBCellFlag const> const& flag,
                         amrex::Array4<amrex::Real const> const& vfrac,
                         AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                                      amrex::Array4<amrex::Real const> const& fcy,
                                      amrex::Array4<amrex::Real const> const& fcz),
                         amrex::Array4<amrex::Real const> const& ccc,
                         const amrex::Geometry& geom,
                         amrex::Real dt,
                         amrex::Vector<amrex::BCRec> const& h_bcrec,
                         amrex::BCRec const* pbc,
                         Array4<int const> const& bc_arr = {});

#if (AMREX_SPACEDIM==3)
void PredictVelOnZFace ( amrex::Box const& zebox,
                         amrex::Array4<amrex::Real> const& Imz,
                         amrex::Array4<amrex::Real> const& Ipz,
                         amrex::Array4<amrex::Real const> const& q,
                         amrex::Array4<amrex::Real  const> const& ccvel,
                         amrex::Array4<amrex::EBCellFlag const> const& flag,
                         amrex::Array4<amrex::Real const> const& vfrac,
                         AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                                      amrex::Array4<amrex::Real const> const& fcy,
                                      amrex::Array4<amrex::Real const> const& fcz),
                         amrex::Array4<amrex::Real const> const& ccc,
                         const amrex::Geometry& geom,
                         amrex::Real dt,
                         amrex::Vector<amrex::BCRec> const& h_bcrec,
                         amrex::BCRec const* pbc,
                         Array4<int const> const& bc_arr = {});
#endif



void PredictStateOnXFace ( amrex::Box const& xebox, int ncomp,
                           amrex::Array4<amrex::Real> const& Imx, amrex::Array4<amrex::Real> const& Ipx,
                           amrex::Array4<amrex::Real const> const& q,
                           amrex::Array4<amrex::Real  const> const& umac,
                           amrex::Array4<amrex::EBCellFlag const> const& flag,
                           amrex::Array4<amrex::Real const> const& vfrac,
                           AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                                        amrex::Array4<amrex::Real const> const& fcy,
                                        amrex::Array4<amrex::Real const> const& fcz),
                           amrex::Array4<amrex::Real const> const& ccc,
                           amrex::Geometry const& geom,
                           amrex::Real dt,
                           amrex::Vector<amrex::BCRec> const& h_bcrec,
                           amrex::BCRec const* pbc, bool is_velocity);

void PredictStateOnYFace ( amrex::Box const& yebox, int ncomp,
                           amrex::Array4<amrex::Real> const& Imy, amrex::Array4<amrex::Real> const& Ipy,
                           amrex::Array4<amrex::Real const> const& q,
                           amrex::Array4<amrex::Real  const> const& vmac,
                           amrex::Array4<amrex::EBCellFlag const> const& flag,
                           amrex::Array4<amrex::Real const> const& vfrac,
                           AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                                        amrex::Array4<amrex::Real const> const& fcy,
                                        amrex::Array4<amrex::Real const> const& fcz),
                           amrex::Array4<amrex::Real const> const& ccc,
                           amrex::Geometry const& geom,
                           amrex::Real dt,
                           amrex::Vector<amrex::BCRec> const& h_bcrec,
                           amrex::BCRec const* pbc, bool is_velocity);

#if (AMREX_SPACEDIM == 3)
void PredictStateOnZFace ( amrex::Box const& zebox, int ncomp,
                           amrex::Array4<amrex::Real> const& Imz, amrex::Array4<amrex::Real> const& Ipz,
                           amrex::Array4<amrex::Real const> const& q,
                           amrex::Array4<amrex::Real  const> const& wmac,
                           amrex::Array4<amrex::EBCellFlag const> const& flag,
                           amrex::Array4<amrex::Real const> const& vfrac,
                           AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                                        amrex::Array4<amrex::Real const> const& fcy,
                                        amrex::Array4<amrex::Real const> const& fcz),
                           amrex::Array4<amrex::Real const> const& ccc,
                           amrex::Geometry const& geom,
                           amrex::Real dt,
                           amrex::Vector<amrex::BCRec> const& h_bcrec,
                           amrex::BCRec const* pbc, bool is_velocity);
#endif

}
#endif
/** @} */

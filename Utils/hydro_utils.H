#ifndef HYDRO_UTILS_H
#define HYDRO_UTILS_H

#include <AMReX_MultiFabUtil.H>
#include <AMReX_BCRec.H>

namespace HydroUtils {


void ComputeFluxes ( amrex::Box const& bx,
                     AMREX_D_DECL( amrex::Array4<amrex::Real> const& fx,
                                   amrex::Array4<amrex::Real> const& fy,
                                   amrex::Array4<amrex::Real> const& fz),
                     AMREX_D_DECL( amrex::Array4<amrex::Real const> const& umac,
                                   amrex::Array4<amrex::Real const> const& vmac,
                                   amrex::Array4<amrex::Real const> const& wmac),
                     AMREX_D_DECL( amrex::Array4<amrex::Real const> const& xedge,
                                   amrex::Array4<amrex::Real const> const& yedge,
                                   amrex::Array4<amrex::Real const> const& zedge),
                     amrex::Geometry const& geom, const int ncomp );

void ComputeDivergence ( amrex::Box const& bx,
                         amrex::Array4<amrex::Real> const& div,
                         AMREX_D_DECL( amrex::Array4<amrex::Real const> const& fx,
                                       amrex::Array4<amrex::Real const> const& fy,
                                       amrex::Array4<amrex::Real const> const& fz),
                         AMREX_D_DECL( amrex::Array4<amrex::Real const> const& xed,
                                       amrex::Array4<amrex::Real const> const& yed,
                                       amrex::Array4<amrex::Real const> const& zed),
                         AMREX_D_DECL( amrex::Array4<amrex::Real const> const& umac,
                                       amrex::Array4<amrex::Real const> const& vmac,
                                       amrex::Array4<amrex::Real const> const& wmac),
                         const int ncomp, amrex::Geometry const& geom,
                         int const* iconserv );


void ComputeSyncDivergence ( amrex::Box const& bx,
                             amrex::Array4<amrex::Real> const& div,
                             AMREX_D_DECL( amrex::Array4<amrex::Real const> const& fx,
                                           amrex::Array4<amrex::Real const> const& fy,
                                           amrex::Array4<amrex::Real const> const& fz),
                             const int ncomp, amrex::Geometry const& geom );


#ifdef AMREX_USE_EB

void EB_ComputeFluxes ( amrex::Box const& bx,
                        AMREX_D_DECL( amrex::Array4<amrex::Real> const& fx,
                                      amrex::Array4<amrex::Real> const& fy,
                                      amrex::Array4<amrex::Real> const& fz),
                        AMREX_D_DECL( amrex::Array4<amrex::Real const> const& umac,
                                      amrex::Array4<amrex::Real const> const& vmac,
                                      amrex::Array4<amrex::Real const> const& wmac),
                        AMREX_D_DECL( amrex::Array4<amrex::Real const> const& xedge,
                                      amrex::Array4<amrex::Real const> const& yedge,
                                      amrex::Array4<amrex::Real const> const& zedge),
                        AMREX_D_DECL( amrex::Array4<amrex::Real const> const& apx,
                                      amrex::Array4<amrex::Real const> const& apy,
                                      amrex::Array4<amrex::Real const> const& apz),
                        amrex::Geometry const& geom, const int ncomp,
                        amrex::Array4<amrex::EBCellFlag const> const& flag);


void EB_ComputeDivergence ( amrex::Box const& bx,
                            amrex::Array4<amrex::Real> const& div,
                            AMREX_D_DECL( amrex::Array4<amrex::Real const> const& fx,
                                          amrex::Array4<amrex::Real const> const& fy,
                                          amrex::Array4<amrex::Real const> const& fz),
                            amrex::Array4<amrex::Real const> const& vfrac,
                            const int ncomp, amrex::Geometry const& geom );

#endif
}
#endif

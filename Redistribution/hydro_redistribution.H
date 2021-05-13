#ifndef IAMR_REDISTRIBUTION_H_
#define IAMR_REDISTRIBUTION_H_

#include <AMReX_MultiFabUtil.H>
#include <AMReX_MultiCutFab.H>

namespace Redistribution {

    void Apply ( amrex::Box const& bx, int ncomp,
                 amrex::Array4<amrex::Real>       const& dUdt_out,
                 amrex::Array4<amrex::Real>       const& dUdt_in,
                 amrex::Array4<amrex::Real const> const& U_in,
                 amrex::Array4<amrex::Real> const& scratch,
                 amrex::Array4<amrex::EBCellFlag const> const& flag,
                 AMREX_D_DECL(amrex::Array4<amrex::Real const> const& apx,
                              amrex::Array4<amrex::Real const> const& apy,
                              amrex::Array4<amrex::Real const> const& apz),
                 amrex::Array4<amrex::Real const> const& vfrac,
                 AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                              amrex::Array4<amrex::Real const> const& fcy,
                              amrex::Array4<amrex::Real const> const& fcz),
                 amrex::Array4<amrex::Real const> const& ccent,
                 AMREX_D_DECL(bool extdir_ilo, bool extdir_jlo, bool extdir_klo),
                 AMREX_D_DECL(bool extdir_ihi, bool extdir_jhi, bool extdir_khi),
                 amrex::Geometry const& geom, 
                 amrex::Real dt, std::string redistribution_type);

    void ApplyToInitialData ( amrex::Box const& bx, int ncomp,
                              amrex::Array4<amrex::Real                  > const& U_out,
                              amrex::Array4<amrex::Real                  > const& U_in,
                              amrex::Array4<amrex::EBCellFlag const> const& flag,
                              AMREX_D_DECL(amrex::Array4<amrex::Real const> const& apx,
                                           amrex::Array4<amrex::Real const> const& apy,
                                           amrex::Array4<amrex::Real const> const& apz),
                              amrex::Array4<amrex::Real const> const& vfrac,
                              AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                                           amrex::Array4<amrex::Real const> const& fcy,
                                           amrex::Array4<amrex::Real const> const& fcz),
                              amrex::Array4<amrex::Real const> const& ccent,
                              AMREX_D_DECL(bool extdir_ilo, bool extdir_jlo, bool extdir_klo),
                              AMREX_D_DECL(bool extdir_ihi, bool extdir_jhi, bool extdir_khi),
                              amrex::Geometry& geom, std::string redistribution_type);

    void FluxRedistribute ( amrex::Box const& bx, int ncomp,
                            amrex::Array4<amrex::Real>       const& dUdt_out,
                            amrex::Array4<amrex::Real const> const& dUdt_in,
                            amrex::Array4<amrex::Real> const& scratch,
                            amrex::Array4<amrex::EBCellFlag const> const& flag,
                            amrex::Array4<amrex::Real const> const& vfrac,
                            amrex::Geometry const& geom);

    void StateRedistribute ( amrex::Box const& bx, int ncomp,
                             amrex::Array4<amrex::Real> const& dUdt_out,
                             amrex::Array4<amrex::Real> const& dUdt_in,
                             amrex::Array4<amrex::EBCellFlag const> const& flag,
                             amrex::Array4<amrex::Real const> const& vfrac,
                             AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                                          amrex::Array4<amrex::Real const> const& fcy,
                                          amrex::Array4<amrex::Real const> const& fcz),
                             amrex::Array4<amrex::Real const> const& ccent,
                             AMREX_D_DECL(bool extdir_ilo, bool extdir_jlo, bool extdir_klo),
                             AMREX_D_DECL(bool extdir_ihi, bool extdir_jhi, bool extdir_khi),
                             amrex::Array4<int> const& itracker,
                             amrex::Geometry const& geom);

    void MakeITracker ( amrex::Box const& bx,
                        AMREX_D_DECL(amrex::Array4<amrex::Real const> const& apx,
                                     amrex::Array4<amrex::Real const> const& apy,
                                     amrex::Array4<amrex::Real const> const& apz),
                        amrex::Array4<amrex::Real const> const& vfrac,
                        amrex::Array4<int> const& itracker,
                        amrex::Geometry const& geom,
                        std::string redist_type);

} // namespace redistribution

#endif

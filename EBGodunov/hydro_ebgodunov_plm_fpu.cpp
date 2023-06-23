/**
 * \file hydro_ebgodunov_plm_fpu.cpp
 * \addtogroup EBGodunov
 * @{
 *
 */

#include <hydro_ebgodunov_plm.H>
#include <AMReX_Slopes_K.H>
#if (AMREX_SPACEDIM == 2)
#include <AMReX_EB_Slopes_2D_K.H>
#elif (AMREX_SPACEDIM == 3)
#include <AMReX_EB_Slopes_3D_K.H>
#endif

using namespace amrex;

namespace {
    std::pair<bool,bool> has_extdir_or_ho (BCRec const* bcrec, int ncomp, int dir)
    {
        std::pair<bool,bool> r{false,false};
        for (int n = 0; n < ncomp; ++n) {
            r.first = r.first
                 || (bcrec[n].lo(dir) == BCType::ext_dir)
                 || (bcrec[n].lo(dir) == BCType::hoextrap);
            r.second = r.second
                 || (bcrec[n].hi(dir) == BCType::ext_dir)
                 || (bcrec[n].hi(dir) == BCType::hoextrap);
        }
        return r;
    }
}

// This version is called after the MAC projection
void
EBPLM::PredictStateOnXFace (Box const& xebox, int ncomp,
                            Array4<Real> const& Imx, Array4<Real> const& Ipx,
                            Array4<Real const> const& q,
                            Array4<Real const> const& umac,
                            Array4<EBCellFlag const> const& flag,
                            Array4<Real const> const& vfrac,
                            AMREX_D_DECL(Array4<Real const> const& fcx,
                                         Array4<Real const> const& fcy,
                                         Array4<Real const> const& fcz),
                            Array4<Real const> const& ccc,
                            Geometry const& geom,
                            Real dt,
                            Vector<BCRec> const& h_bcrec,
                            BCRec const* pbc, bool is_velocity)
{
    const Real dx = geom.CellSize(0);
    const Real dtdx = dt/dx;

    const Box& domain_box = geom.Domain();
    const int domain_ilo = domain_box.smallEnd(0);
    const int domain_ihi = domain_box.bigEnd(0);
    const int domain_jlo = domain_box.smallEnd(1);
    const int domain_jhi = domain_box.bigEnd(1);

    // At an ext_dir boundary, the boundary value is on the face, not cell center.
    auto extdir_lohi_x = has_extdir_or_ho(h_bcrec.data(), ncomp, static_cast<int>(Direction::x));
    auto extdir_lohi_y = has_extdir_or_ho(h_bcrec.data(), ncomp, static_cast<int>(Direction::y));

    bool has_extdir_or_ho_lo_x = extdir_lohi_x.first;
    bool has_extdir_or_ho_hi_x = extdir_lohi_x.second;
    bool has_extdir_or_ho_lo_y = extdir_lohi_y.first;
    bool has_extdir_or_ho_hi_y = extdir_lohi_y.second;

#if (AMREX_SPACEDIM == 3)
    const int domain_klo = domain_box.smallEnd(2);
    const int domain_khi = domain_box.bigEnd(2);
    auto extdir_lohi_z = has_extdir_or_ho(h_bcrec.data(), ncomp, static_cast<int>(Direction::z));
    bool has_extdir_or_ho_lo_z = extdir_lohi_z.first;
    bool has_extdir_or_ho_hi_z = extdir_lohi_z.second;
#endif

    if ( (has_extdir_or_ho_lo_x && domain_ilo >= xebox.smallEnd(0)-1) ||
         (has_extdir_or_ho_hi_x && domain_ihi <= xebox.bigEnd(0)    ) ||
#if (AMREX_SPACEDIM == 3)
         (has_extdir_or_ho_lo_z && domain_klo >= xebox.smallEnd(2)-1) ||
         (has_extdir_or_ho_hi_z && domain_khi <= xebox.bigEnd(2)    ) ||
#endif
         (has_extdir_or_ho_lo_y && domain_jlo >= xebox.smallEnd(1)-1) ||
         (has_extdir_or_ho_hi_y && domain_jhi <= xebox.bigEnd(1)    )  )
    {
        amrex::ParallelFor(xebox, ncomp, [q,umac,AMREX_D_DECL(domain_ilo,domain_jlo,domain_klo),
                                                 AMREX_D_DECL(domain_ihi,domain_jhi,domain_khi),
                                          Imx,Ipx,dtdx,pbc,flag,vfrac,ccc,AMREX_D_DECL(fcx,fcy,fcz),
                                          is_velocity]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real qpls(0.);
            Real qmns(0.);

            // This means apx(i,j,k) > 0 and we have un-covered cells on both sides
            if (flag(i,j,k).isConnected(-1,0,0))
            {
                const auto& bc = pbc[n];
                bool extdir_or_ho_ilo = (bc.lo(0) == BCType::ext_dir) ||
                                        (bc.lo(0) == BCType::hoextrap);
                bool extdir_or_ho_ihi = (bc.hi(0) == BCType::ext_dir) ||
                                        (bc.hi(0) == BCType::hoextrap);
                bool extdir_or_ho_jlo = (bc.lo(1) == BCType::ext_dir) ||
                                        (bc.lo(1) == BCType::hoextrap);
                bool extdir_or_ho_jhi = (bc.hi(1) == BCType::ext_dir) ||
                                        (bc.hi(1) == BCType::hoextrap);
#if (AMREX_SPACEDIM == 3)
                bool extdir_or_ho_klo = (bc.lo(2) == BCType::ext_dir) ||
                                        (bc.lo(2) == BCType::hoextrap);
                bool extdir_or_ho_khi = (bc.hi(2) == BCType::ext_dir) ||
                                        (bc.hi(2) == BCType::hoextrap);
#endif

                // *************************************************
                // Making qpls
                // *************************************************

                // We have enough cells to do 4th order slopes centered on (i,j,k) with all values at cell centers
                if (vfrac(i,j,k) == 1. && vfrac(i-1,j,k) == 1. && vfrac(i-2,j,k) == 1. &&
                                           vfrac(i+1,j,k) == 1. && vfrac(i+2,j,k) == 1.)
                {
                    int order = 4;
                    qpls = q(i,j,k,n) + Real(0.5) * (-Real(1.0) - umac(i,j,k,0) * dtdx) *
                        amrex_calc_xslope_extdir(i  ,j,k,n,order,q,extdir_or_ho_ilo,extdir_or_ho_ihi,domain_ilo,domain_ihi);

                // We have enough cells to do 2nd order slopes with all values at cell centers
                } else if (vfrac(i,j,k) == 1. && vfrac(i-1,j,k) == 1. && vfrac(i+1,j,k) == 1.) {

                    int order = 2;
                    qpls = q(i,j,k,n) + Real(0.5) * (-Real(1.0) - umac(i,j,k,0) * dtdx) *
                        amrex_calc_xslope_extdir(i  ,j,k,n,order,q,extdir_or_ho_ilo,extdir_or_ho_ihi,domain_ilo,domain_ihi);

                // We need to use LS slopes
                } else {

                   Real yf = fcx(i,j,k,0); // local (y,z) of centroid of x-face we are extrapolating to
#if (AMREX_SPACEDIM == 3)
                   Real zf = fcx(i,j,k,1);
#endif
                   AMREX_D_TERM(Real delta_x = Real(0.5) + ccc(i,j,k,0);,
                                Real delta_y = yf  - ccc(i,j,k,1);,
                                Real delta_z = zf  - ccc(i,j,k,2););

                   Real qcc_max = amrex::max(q(i,j,k,n), q(i-1,j,k,n));
                   Real qcc_min = amrex::min(q(i,j,k,n), q(i-1,j,k,n));

                   // This will be used in the EB slope routine only if the slope can be computed without LS
                   int max_order = 2;

                   const auto& slopes_eb_hi = amrex_lim_slopes_extdir_eb(i,j,k,n,q,ccc,vfrac,
                                              AMREX_D_DECL(fcx,fcy,fcz), flag,
                                              AMREX_D_DECL(extdir_or_ho_ilo, extdir_or_ho_jlo, extdir_or_ho_klo),
                                              AMREX_D_DECL(extdir_or_ho_ihi, extdir_or_ho_jhi, extdir_or_ho_khi),
                                              AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
                                              AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi),
                                              max_order);

#if (AMREX_SPACEDIM == 3)
                   qpls = q(i,j,k,n) - delta_x * slopes_eb_hi[0]
                                     + delta_y * slopes_eb_hi[1]
                                     + delta_z * slopes_eb_hi[2];
#else
                   qpls = q(i,j,k,n) - delta_x * slopes_eb_hi[0]
                                     + delta_y * slopes_eb_hi[1];
#endif
                   qpls = amrex::max(amrex::min(qpls, qcc_max), qcc_min);
                   qpls -= Real(0.5) * dtdx * umac(i,j,k) * slopes_eb_hi[0];

                }  // end of making qpls

                // Only over-write normal velocity with Dirichlet bc at lo face
                if (i == domain_ilo && (bc.lo(0) == BCType::ext_dir))
                    if (is_velocity && n == 0) qpls = q(i-1,j,k,n);

                // Over-write all with Dirichlet bc at hi face
                if (i == domain_ihi+1 && (bc.hi(0) == BCType::ext_dir))
                    qpls = q(i,j,k,n);

                // *************************************************
                // Making qmns
                // *************************************************

                // We have enough cells to do 4th order slopes centered on (i-1,j,k) with all values at cell centers
                if (vfrac(i-1,j,k) == 1. && vfrac(i-2,j,k) == 1. && vfrac(i-3,j,k) == 1. &&
                                            vfrac(i  ,j,k) == 1. && vfrac(i+1,j,k) == 1.)
                {
                    int order = 4;
                    qmns = q(i-1,j,k,n) + Real(0.5) * ( Real(1.0) - umac(i,j,k) * dtdx) *
                        amrex_calc_xslope_extdir(i-1,j,k,n,order,q,extdir_or_ho_ilo,extdir_or_ho_ihi,domain_ilo,domain_ihi);

                // We have enough cells to do 2nd order slopes with all values at cell centers
                } else if (vfrac(i-1,j,k) == 1. && vfrac(i-2,j,k) == 1. && vfrac(i  ,j,k) == 1.)
                {
                    int order = 2;
                    qmns = q(i-1,j,k,n) + Real(0.5) * ( Real(1.0) - umac(i,j,k) * dtdx) *
                        amrex_calc_xslope_extdir(i-1,j,k,n,order,q,extdir_or_ho_ilo,extdir_or_ho_ihi,domain_ilo,domain_ihi);

                // We need to use LS slopes
                } else {

                   Real yf = fcx(i,j,k,0); // local (y,z) of centroid of x-face we are extrapolating to
#if (AMREX_SPACEDIM == 3)
                   Real zf = fcx(i,j,k,1);
#endif
                   AMREX_D_TERM(Real delta_x = Real(0.5) - ccc(i-1,j,k,0);,
                                Real delta_y = yf  - ccc(i-1,j,k,1);,
                                Real delta_z = zf  - ccc(i-1,j,k,2););

                   Real qcc_max = amrex::max(q(i,j,k,n), q(i-1,j,k,n));
                   Real qcc_min = amrex::min(q(i,j,k,n), q(i-1,j,k,n));

                   // This will be used in the EB slope routine only if the slope can be computed without LS
                   int max_order = 2;

                   const auto& slopes_eb_lo = amrex_lim_slopes_extdir_eb(i-1,j,k,n,q,ccc,vfrac,
                                              AMREX_D_DECL(fcx,fcy,fcz), flag,
                                              AMREX_D_DECL(extdir_or_ho_ilo, extdir_or_ho_jlo, extdir_or_ho_klo),
                                              AMREX_D_DECL(extdir_or_ho_ihi, extdir_or_ho_jhi, extdir_or_ho_khi),
                                              AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
                                              AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi),
                                              max_order);

#if (AMREX_SPACEDIM == 3)
                   qmns = q(i-1,j,k,n) + delta_x * slopes_eb_lo[0]
                                       + delta_y * slopes_eb_lo[1]
                                       + delta_z * slopes_eb_lo[2];
#else
                   qmns = q(i-1,j,k,n) + delta_x * slopes_eb_lo[0]
                                       + delta_y * slopes_eb_lo[1];
#endif
                   qmns = amrex::max(amrex::min(qmns, qcc_max), qcc_min);
                   qmns -= Real(0.5) * dtdx * umac(i,j,k) * slopes_eb_lo[0];
                }  // end of making qmns

                // Over-write all with Dirichlet bc at lo face
                if (i == domain_ilo && (bc.lo(0) == BCType::ext_dir))
                    qmns = q(i-1,j,k,n);

                // Only over-write normal velocity with Dirichlet bc at hi face
                if (i == domain_ihi+1 && (bc.hi(0) == BCType::ext_dir))
                    if (is_velocity && n == 0) qmns = q(i,j,k,n);
            }

            Ipx(i-1,j,k,n) = qmns;
            Imx(i  ,j,k,n) = qpls;
        });
    }
    else // The cases below are not near any domain boundary
    {
        amrex::ParallelFor(xebox, ncomp, [q,umac,AMREX_D_DECL(domain_ilo,domain_jlo,domain_klo),
                                                AMREX_D_DECL(domain_ihi,domain_jhi,domain_khi),
                                          Imx,Ipx,dtdx,pbc,flag,vfrac,ccc,AMREX_D_DECL(fcx,fcy,fcz),
                                          is_velocity]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real qpls(0.);
            Real qmns(0.);

            // This means apx(i,j,k) > 0 and we have un-covered cells on both sides
            if (flag(i,j,k).isConnected(-1,0,0))
            {
                // *************************************************
                // Making qpls
                // *************************************************

                // We have enough cells to do 4th order slopes centered on (i,j,k) with all values at cell centers
                if (vfrac(i,j,k) == 1. && vfrac(i-1,j,k) == 1. && vfrac(i-2,j,k) == 1. &&
                                          vfrac(i+1,j,k) == 1. && vfrac(i+2,j,k) == 1.)
                {
                    int order = 4;
                    qpls = q(i  ,j,k,n) + Real(0.5) * (-Real(1.0) - umac(i,j,k,0) * dtdx) *
                        amrex_calc_xslope(i  ,j,k,n,order,q);

                // We have enough cells to do 2nd order slopes with all values at cell centers
                } else if (vfrac(i,j,k) == 1. && vfrac(i-1,j,k) == 1. && vfrac(i+1,j,k) == 1.) {

                    int order = 2;
                    qpls = q(i  ,j,k,n) + Real(0.5) * (-Real(1.0) - umac(i,j,k,0) * dtdx) *
                        amrex_calc_xslope(i  ,j,k,n,order,q);

                // We need to use LS slopes
                } else {

                   Real yf = fcx(i,j,k,0); // local (y,z) of centroid of x-face we are extrapolating to
#if (AMREX_SPACEDIM == 3)
                   Real zf = fcx(i,j,k,1);
#endif
                   AMREX_D_TERM(Real delta_x = Real(0.5) + ccc(i,j,k,0);,
                                Real delta_y = yf  - ccc(i,j,k,1);,
                                Real delta_z = zf  - ccc(i,j,k,2););

                   Real qcc_max = amrex::max(q(i,j,k,n), q(i-1,j,k,n));
                   Real qcc_min = amrex::min(q(i,j,k,n), q(i-1,j,k,n));

                   // This will be used in the EB slope routine only if the slope can be computed without LS
                   int max_order = 2;

                   const auto& slopes_eb_hi = amrex_lim_slopes_eb(i,j,k,n,q,ccc,vfrac,
                                                                  AMREX_D_DECL(fcx,fcy,fcz),flag,max_order);

#if (AMREX_SPACEDIM == 3)
                   qpls = q(i,j,k,n) - delta_x * slopes_eb_hi[0]
                                     + delta_y * slopes_eb_hi[1]
                                     + delta_z * slopes_eb_hi[2];
#else
                   qpls = q(i,j,k,n) - delta_x * slopes_eb_hi[0]
                                     + delta_y * slopes_eb_hi[1];
#endif
                   qpls = amrex::max(amrex::min(qpls, qcc_max), qcc_min);
                   qpls -= Real(0.5) * dtdx * umac(i,j,k) * slopes_eb_hi[0];
                }  // end of making qpls

                // *************************************************
                // Making qmns
                // *************************************************

                // We have enough cells to do 4th order slopes centered on (i-1,j,k) with all values at cell centers
                if (vfrac(i-1,j,k) == 1. && vfrac(i-2,j,k) == 1. && vfrac(i-3,j,k) == 1. &&
                                            vfrac(i  ,j,k) == 1. && vfrac(i+1,j,k) == 1.)
                {
                    int order = 4;
                    qmns = q(i-1,j,k,n) + Real(0.5) * ( Real(1.0) - umac(i,j,k) * dtdx) *
                        amrex_calc_xslope(i-1,j,k,n,order,q);

                // We have enough cells to do 2nd order slopes with all values at cell centers
                } else if (vfrac(i-1,j,k) == 1. && vfrac(i-2,j,k) == 1. && vfrac(i  ,j,k) == 1.)
                {
                    int order = 2;
                    qmns = q(i-1,j,k,n) + Real(0.5) * ( Real(1.0) - umac(i,j,k) * dtdx) *
                        amrex_calc_xslope(i-1,j,k,n,order,q);

                // We need to use LS slopes
                } else {

                   Real yf = fcx(i,j,k,0); // local (y,z) of centroid of x-face we are extrapolating to
#if (AMREX_SPACEDIM == 3)
                   Real zf = fcx(i,j,k,1);
#endif
                   AMREX_D_TERM(Real delta_x = Real(0.5) - ccc(i-1,j,k,0);,
                                Real delta_y = yf  - ccc(i-1,j,k,1);,
                                Real delta_z = zf  - ccc(i-1,j,k,2););

                   Real qcc_max = amrex::max(q(i,j,k,n), q(i-1,j,k,n));
                   Real qcc_min = amrex::min(q(i,j,k,n), q(i-1,j,k,n));

                   // This will be used in the EB slope routine only if the slope can be computed without LS
                   int max_order = 2;

                   const auto& slopes_eb_lo = amrex_lim_slopes_eb(i-1,j,k,n,q,ccc,vfrac,
                                                                  AMREX_D_DECL(fcx,fcy,fcz),flag,max_order);

#if (AMREX_SPACEDIM == 3)
                   qmns = q(i-1,j,k,n) + delta_x * slopes_eb_lo[0]
                                       + delta_y * slopes_eb_lo[1]
                                       + delta_z * slopes_eb_lo[2];
#else
                   qmns = q(i-1,j,k,n) + delta_x * slopes_eb_lo[0]
                                       + delta_y * slopes_eb_lo[1];
#endif
                   qmns = amrex::max(amrex::min(qmns, qcc_max), qcc_min);
                   qmns -= Real(0.5) * dtdx * umac(i,j,k) * slopes_eb_lo[0];

                }  // end of making qmns
            }

            Ipx(i-1,j,k,n) = qmns;
            Imx(i  ,j,k,n) = qpls;
        });
    }
}

// This version is called after the MAC projection
void
EBPLM::PredictStateOnYFace ( Box const& yebox, int ncomp,
                             Array4<Real> const& Imy, Array4<Real> const& Ipy,
                             Array4<Real const> const& q,
                             Array4<Real const> const& vmac,
                             Array4<EBCellFlag const> const& flag,
                             Array4<Real const> const& vfrac,
                             AMREX_D_DECL(Array4<Real const> const& fcx,
                                          Array4<Real const> const& fcy,
                                          Array4<Real const> const& fcz),
                             Array4<Real const> const& ccc,
                             Geometry const& geom,
                             Real dt,
                             Vector<BCRec> const& h_bcrec,
                             BCRec const* pbc, bool is_velocity)
{
    const Real dy = geom.CellSize(1);
    const Real dtdy = dt/dy;

    const Box& domain_box = geom.Domain();
    const int domain_ilo = domain_box.smallEnd(0);
    const int domain_ihi = domain_box.bigEnd(0);
    const int domain_jlo = domain_box.smallEnd(1);
    const int domain_jhi = domain_box.bigEnd(1);

    // At an ext_dir boundary, the boundary value is on the face, not cell center.
    auto extdir_lohi_x = has_extdir_or_ho(h_bcrec.data(), ncomp, static_cast<int>(Direction::x));
    auto extdir_lohi_y = has_extdir_or_ho(h_bcrec.data(), ncomp, static_cast<int>(Direction::y));

    bool has_extdir_or_ho_lo_x = extdir_lohi_x.first;
    bool has_extdir_or_ho_hi_x = extdir_lohi_x.second;
    bool has_extdir_or_ho_lo_y = extdir_lohi_y.first;
    bool has_extdir_or_ho_hi_y = extdir_lohi_y.second;

#if (AMREX_SPACEDIM == 3)
    const int domain_klo = domain_box.smallEnd(2);
    const int domain_khi = domain_box.bigEnd(2);
    auto extdir_lohi_z = has_extdir_or_ho(h_bcrec.data(), ncomp, static_cast<int>(Direction::z));
    bool has_extdir_or_ho_lo_z = extdir_lohi_z.first;
    bool has_extdir_or_ho_hi_z = extdir_lohi_z.second;
#endif

    if ( (has_extdir_or_ho_lo_x && domain_ilo >= yebox.smallEnd(0)-1) ||
         (has_extdir_or_ho_hi_x && domain_ihi <= yebox.bigEnd(0)    ) ||
#if (AMREX_SPACEDIM == 3)
         (has_extdir_or_ho_lo_z && domain_klo >= yebox.smallEnd(2)-1) ||
         (has_extdir_or_ho_hi_z && domain_khi <= yebox.bigEnd(2)    ) ||
#endif
         (has_extdir_or_ho_lo_y && domain_jlo >= yebox.smallEnd(1)-1) ||
         (has_extdir_or_ho_hi_y && domain_jhi <= yebox.bigEnd(1)    )  )
    {
        amrex::ParallelFor(yebox, ncomp, [q,vmac,AMREX_D_DECL(domain_ilo,domain_jlo,domain_klo),
                                                AMREX_D_DECL(domain_ihi,domain_jhi,domain_khi),
                                          Imy,Ipy,dtdy,pbc,flag,vfrac,ccc,AMREX_D_DECL(fcx,fcy,fcz),
                                          is_velocity]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real qpls(0.);
            Real qmns(0.);

            // This means apy(i,j,k) > 0 and we have un-covered cells on both sides
            if (flag(i,j,k).isConnected(0,-1,0))
            {
                const auto& bc = pbc[n];
                bool extdir_or_ho_ilo = (bc.lo(0) == BCType::ext_dir) ||
                                        (bc.lo(0) == BCType::hoextrap);
                bool extdir_or_ho_ihi = (bc.hi(0) == BCType::ext_dir) ||
                                        (bc.hi(0) == BCType::hoextrap);
                bool extdir_or_ho_jlo = (bc.lo(1) == BCType::ext_dir) ||
                                        (bc.lo(1) == BCType::hoextrap);
                bool extdir_or_ho_jhi = (bc.hi(1) == BCType::ext_dir) ||
                                        (bc.hi(1) == BCType::hoextrap);
#if (AMREX_SPACEDIM == 3)
                bool extdir_or_ho_klo = (bc.lo(2) == BCType::ext_dir) ||
                                        (bc.lo(2) == BCType::hoextrap);
                bool extdir_or_ho_khi = (bc.hi(2) == BCType::ext_dir) ||
                                        (bc.hi(2) == BCType::hoextrap);
#endif

                // *************************************************
                // Making qpls
                // *************************************************

                // We have enough cells to do 4th order slopes centered on (i,j,k) with all values at cell centers
                if (vfrac(i,j,k) == 1. && vfrac(i,j-1,k) == 1. && vfrac(i,j-2,k) == 1. &&
                                          vfrac(i,j+1,k) == 1. && vfrac(i,j+2,k) == 1.)
                {
                    int order = 4;
                    qpls = q(i,j  ,k,n) + Real(0.5) * (-Real(1.0) - vmac(i,j,k) * dtdy) *
                        amrex_calc_yslope_extdir(i,j,k,n,order,q,extdir_or_ho_jlo,extdir_or_ho_jhi,domain_jlo,domain_jhi);

                // We have enough cells to do 2nd order slopes with all values at cell centers
                } else if (vfrac(i,j,k) == 1. && vfrac(i,j-1,k) == 1. && vfrac(i,j+1,k) == 1.) {

                    int order = 2;
                    qpls = q(i,j  ,k,n) + Real(0.5) * (-Real(1.0) - vmac(i,j,k) * dtdy) *
                        amrex_calc_yslope_extdir(i,j,k,n,order,q,extdir_or_ho_jlo,extdir_or_ho_jhi,domain_jlo,domain_jhi);

                // We need to use LS slopes
                } else {

                   Real xf = fcy(i,j,k,0); // local (x,z) of centroid of y-face we are extrapolating to
#if (AMREX_SPACEDIM == 3)
                   Real zf = fcy(i,j,k,1);
#endif
                   AMREX_D_TERM(Real delta_y = Real(0.5) + ccc(i,j,k,1);,
                                Real delta_x = xf  - ccc(i,j,k,0);,
                                Real delta_z = zf  - ccc(i,j,k,2););

                   Real qcc_max = amrex::max(q(i,j,k,n), q(i,j-1,k,n));
                   Real qcc_min = amrex::min(q(i,j,k,n), q(i,j-1,k,n));

                   // This will be used in the EB slope routine only if the slope can be computed without LS
                   int max_order = 2;

                   const auto& slopes_eb_hi = amrex_lim_slopes_extdir_eb(i,j,k,n,q,ccc,vfrac,
                                              AMREX_D_DECL(fcx,fcy,fcz), flag,
                                              AMREX_D_DECL(extdir_or_ho_ilo, extdir_or_ho_jlo, extdir_or_ho_klo),
                                              AMREX_D_DECL(extdir_or_ho_ihi, extdir_or_ho_jhi, extdir_or_ho_khi),
                                              AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
                                              AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi),
                                              max_order);

#if (AMREX_SPACEDIM == 3)
                   qpls = q(i,j,k,n) - delta_y * slopes_eb_hi[1]
                                     + delta_x * slopes_eb_hi[0]
                                     + delta_z * slopes_eb_hi[2];
#else
                   qpls = q(i,j,k,n) - delta_y * slopes_eb_hi[1]
                                     + delta_x * slopes_eb_hi[0];
#endif
                   qpls = amrex::max(amrex::min(qpls, qcc_max), qcc_min);
                   qpls -= Real(0.5) * dtdy * vmac(i,j,k) * slopes_eb_hi[1];

                }  // end of making qpls

                // Only over-write normal velocity with Dirichlet bc at lo face
                if (j == domain_jlo && (bc.lo(1) == BCType::ext_dir))
                    if (is_velocity && n == 1) qpls = q(i,j-1,k,n);

                // Over-write all with Dirichlet bc at hi face
                if (j == domain_jhi+1 && (bc.hi(1) == BCType::ext_dir))
                    qpls = q(i,j,k,n);

                // *************************************************
                // Making qmns
                // *************************************************

                // We have enough cells to do 4th order slopes centered on (i,j-1,k) with all values at cell centers
                if (vfrac(i,j-1,k) == 1. && vfrac(i,j-2,k) == 1. && vfrac(i,j-3,k) == 1. &&
                                            vfrac(i,j  ,k) == 1. && vfrac(i,j+1,k) == 1.)
                {
                    int order = 4;
                    qmns = q(i,j-1,k,n) + Real(0.5) * ( Real(1.0) - vmac(i,j,k) * dtdy) *
                        amrex_calc_yslope_extdir(i,j-1,k,n,order,q,extdir_or_ho_jlo,extdir_or_ho_jhi,domain_jlo,domain_jhi);

                // We have enough cells to do 2nd order slopes with all values at cell centers
                } else if (vfrac(i,j-1,k) == 1. && vfrac(i,j-2,k) == 1. && vfrac(i,j  ,k) == 1.)
                {
                    int order = 2;
                    qmns = q(i,j-1,k,n) + Real(0.5) * ( Real(1.0) - vmac(i,j,k) * dtdy) *
                        amrex_calc_yslope_extdir(i,j-1,k,n,order,q,extdir_or_ho_jlo,extdir_or_ho_jhi,domain_jlo,domain_jhi);

                // We need to use LS slopes
                } else {

                   Real xf = fcy(i,j,k,0); // local (x,z) of centroid of y-face we are extrapolating to
#if (AMREX_SPACEDIM == 3)
                   Real zf = fcy(i,j,k,1);
#endif
                   AMREX_D_TERM(Real delta_y = Real(0.5) - ccc(i,j-1,k,1);,
                                Real delta_x = xf  - ccc(i,j-1,k,0);,
                                Real delta_z = zf  - ccc(i,j-1,k,2););

                   Real qcc_max = amrex::max(q(i,j,k,n), q(i,j-1,k,n));
                   Real qcc_min = amrex::min(q(i,j,k,n), q(i,j-1,k,n));

                   // This will be used in the EB slope routine only if the slope can be computed without LS
                   int max_order = 2;

                   const auto& slopes_eb_lo = amrex_lim_slopes_extdir_eb(i,j-1,k,n,q,ccc,vfrac,
                                              AMREX_D_DECL(fcx,fcy,fcz), flag,
                                              AMREX_D_DECL(extdir_or_ho_ilo, extdir_or_ho_jlo, extdir_or_ho_klo),
                                              AMREX_D_DECL(extdir_or_ho_ihi, extdir_or_ho_jhi, extdir_or_ho_khi),
                                              AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
                                              AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi),
                                              max_order);


#if (AMREX_SPACEDIM == 3)
                   qmns = q(i,j-1,k,n) + delta_x * slopes_eb_lo[0]
                                       + delta_y * slopes_eb_lo[1]
                                       + delta_z * slopes_eb_lo[2];
#else
                   qmns = q(i,j-1,k,n) + delta_x * slopes_eb_lo[0]
                                       + delta_y * slopes_eb_lo[1];
#endif
                   qmns = amrex::max(amrex::min(qmns, qcc_max), qcc_min);
                   qmns -= Real(0.5) * dtdy * vmac(i,j,k) * slopes_eb_lo[1];

                }  // end of making qmns

                // Over-write all with Dirichlet bc at lo face
                if (j == domain_jlo && (bc.lo(1) == BCType::ext_dir))
                    qmns = q(i,j-1,k,n);

                // Only over-write normal velocity with Dirichlet bc at hi face
                if (j == domain_jhi+1 && (bc.hi(1) == BCType::ext_dir))
                    if (is_velocity && n == 1) qmns = q(i,j,k,n);
            }

            Ipy(i,j-1,k,n) = qmns;
            Imy(i,j  ,k,n) = qpls;
        });
    }
    else // The cases below are not near any domain boundary
    {
        amrex::ParallelFor(yebox, ncomp, [q,vmac,AMREX_D_DECL(domain_ilo,domain_jlo,domain_klo),
                                                AMREX_D_DECL(domain_ihi,domain_jhi,domain_khi),
                                          Imy,Ipy,dt,dtdy,pbc,flag,vfrac,ccc,AMREX_D_DECL(fcx,fcy,fcz),
                                          is_velocity]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real qpls(0.);
            Real qmns(0.);

            // This means apy(i,j,k) > 0 and we have un-covered cells on both sides
            if (flag(i,j,k).isConnected(0,-1,0))
            {
                // *************************************************
                // Making qpls
                // *************************************************

                // We have enough cells to do 4th order slopes centered on (i,j,k) with all values at cell centers
                if (vfrac(i,j,k) == 1. && vfrac(i,j-1,k) == 1. && vfrac(i,j-2,k) == 1. &&
                                          vfrac(i,j+1,k) == 1. && vfrac(i,j+2,k) == 1.)
                {
                    int order = 4;
                    qpls = q(i,j,k,n) + Real(0.5) * (-Real(1.0) - vmac(i,j,k) * dtdy) *
                        amrex_calc_yslope(i,j,k,n,order,q);

                // We have enough cells to do 2nd order slopes with all values at cell centers
                } else if (vfrac(i,j,k) == 1. && vfrac(i,j-1,k) == 1. && vfrac(i,j+1,k) == 1.) {

                    int order = 2;
                    qpls = q(i,j,k,n) + Real(0.5) * (-Real(1.0) - vmac(i,j,k) * dtdy) *
                        amrex_calc_yslope(i,j,k,n,order,q);

                // We need to use LS slopes
                } else {

                   Real xf = fcy(i,j,k,0); // local (x,z) of centroid of y-face we are extrapolating to
#if (AMREX_SPACEDIM == 3)
                   Real zf = fcy(i,j,k,1);
#endif
                   AMREX_D_TERM(Real delta_y = Real(0.5) + ccc(i,j,k,1);,
                                Real delta_x = xf  - ccc(i,j,k,0);,
                                Real delta_z = zf  - ccc(i,j,k,2););

                   Real qcc_max = amrex::max(q(i,j,k,n), q(i,j-1,k,n));
                   Real qcc_min = amrex::min(q(i,j,k,n), q(i,j-1,k,n));

                   // This will be used in the EB slope routine only if the slope can be computed without LS
                   int max_order = 2;

                   const auto& slopes_eb_hi = amrex_lim_slopes_eb(i,j,k,n,q,ccc,vfrac,
                                                                  AMREX_D_DECL(fcx,fcy,fcz),flag,max_order);

#if (AMREX_SPACEDIM == 3)
                   qpls = q(i,j,k,n) - delta_y * slopes_eb_hi[1]
                                     + delta_x * slopes_eb_hi[0]
                                     + delta_z * slopes_eb_hi[2];
#else
                   qpls = q(i,j,k,n) - delta_y * slopes_eb_hi[1]
                                     + delta_x * slopes_eb_hi[0];
#endif
                   qpls = amrex::max(amrex::min(qpls, qcc_max), qcc_min);
                   qpls -= Real(0.5) * dtdy * vmac(i,j,k) * slopes_eb_hi[1];

                }  // end of making qpls

                // *************************************************
                // Making qmns
                // *************************************************

                // We have enough cells to do 4th order slopes centered on (i,j-1,k) with all values at cell centers
                if (vfrac(i,j-1,k) == 1. && vfrac(i,j-2,k) == 1. && vfrac(i,j-3,k) == 1. &&
                                            vfrac(i,j  ,k) == 1. && vfrac(i,j+1,k) == 1.)
                {
                    int order = 4;
                    qmns = q(i,j-1,k,n) + Real(0.5) * ( Real(1.0) - vmac(i,j,k) * dtdy) *
                        amrex_calc_yslope(i,j-1,k,n,order,q);

                // We have enough cells to do 2nd order slopes with all values at cell centers
                } else if (vfrac(i,j-1,k) == 1. && vfrac(i,j-2,k) == 1. && vfrac(i,j  ,k) == 1.)
                {
                    int order = 2;
                    qmns = q(i,j-1,k,n) + Real(0.5) * ( Real(1.0) - vmac(i,j,k) * dtdy) *
                        amrex_calc_yslope(i,j-1,k,n,order,q);

                // We need to use LS slopes
                } else {

                   Real xf = fcy(i,j,k,0); // local (x,z) of centroid of y-face we are extrapolating to
#if (AMREX_SPACEDIM == 3)
                   Real zf = fcy(i,j,k,1);
#endif
                   AMREX_D_TERM(Real delta_y = Real(0.5) - ccc(i,j-1,k,1);,
                                Real delta_x = xf  - ccc(i,j-1,k,0);,
                                Real delta_z = zf  - ccc(i,j-1,k,2););

                   Real qcc_max = amrex::max(q(i,j,k,n), q(i,j-1,k,n));
                   Real qcc_min = amrex::min(q(i,j,k,n), q(i,j-1,k,n));

                   // This will be used in the EB slope routine only if the slope can be computed without LS
                   int max_order = 2;

                   const auto& slopes_eb_lo = amrex_lim_slopes_eb(i,j-1,k,n,q,ccc,vfrac,
                                                                  AMREX_D_DECL(fcx,fcy,fcz),flag,max_order);

#if (AMREX_SPACEDIM == 3)
                   qmns = q(i,j-1,k,n) + delta_x * slopes_eb_lo[0]
                                       + delta_y * slopes_eb_lo[1]
                                       + delta_z * slopes_eb_lo[2];
#else
                   qmns = q(i,j-1,k,n) + delta_x * slopes_eb_lo[0]
                                       + delta_y * slopes_eb_lo[1];
#endif
                   qmns = amrex::max(amrex::min(qmns, qcc_max), qcc_min);
                   qmns -= Real(0.5) * dtdy * vmac(i,j,k) * slopes_eb_lo[1];

                }  // end of making qmns
            }

            Ipy(i,j-1,k,n) = qmns;
            Imy(i,j  ,k,n) = qpls;
        });
    }
}

#if (AMREX_SPACEDIM == 3)
// This version is called after the MAC projection
void
EBPLM::PredictStateOnZFace ( Box const& zebox, int ncomp,
                             Array4<Real> const& Imz, Array4<Real> const& Ipz,
                             Array4<Real const> const& q,
                             Array4<Real const> const& wmac,
                             Array4<EBCellFlag const> const& flag,
                             Array4<Real const> const& vfrac,
                             AMREX_D_DECL(Array4<Real const> const& fcx,
                                          Array4<Real const> const& fcy,
                                          Array4<Real const> const& fcz),
                             Array4<Real const> const& ccc,
                             Geometry const& geom,
                             Real dt,
                             Vector<BCRec> const& h_bcrec,
                             BCRec const* pbc, bool is_velocity)
{
    const Real dz = geom.CellSize(1);
    const Real dtdz = dt/dz;

    const Box& domain_box = geom.Domain();
    const int domain_ilo = domain_box.smallEnd(0);
    const int domain_ihi = domain_box.bigEnd(0);
    const int domain_jlo = domain_box.smallEnd(1);
    const int domain_jhi = domain_box.bigEnd(1);
    const int domain_klo = domain_box.smallEnd(2);
    const int domain_khi = domain_box.bigEnd(2);

    // At an ext_dir boundary, the boundary value is on the face, not cell center.
    auto extdir_lohi_x = has_extdir_or_ho(h_bcrec.data(), ncomp, static_cast<int>(Direction::x));
    auto extdir_lohi_y = has_extdir_or_ho(h_bcrec.data(), ncomp, static_cast<int>(Direction::y));
    auto extdir_lohi_z = has_extdir_or_ho(h_bcrec.data(), ncomp, static_cast<int>(Direction::z));

    bool has_extdir_or_ho_lo_x = extdir_lohi_x.first;
    bool has_extdir_or_ho_hi_x = extdir_lohi_x.second;
    bool has_extdir_or_ho_lo_y = extdir_lohi_y.first;
    bool has_extdir_or_ho_hi_y = extdir_lohi_y.second;
    bool has_extdir_or_ho_lo_z = extdir_lohi_z.first;
    bool has_extdir_or_ho_hi_z = extdir_lohi_z.second;

    if ( (has_extdir_or_ho_lo_x && domain_ilo >= zebox.smallEnd(0)-1) ||
         (has_extdir_or_ho_hi_x && domain_ihi <= zebox.bigEnd(0)    ) ||
         (has_extdir_or_ho_lo_z && domain_klo >= zebox.smallEnd(2)-1) ||
         (has_extdir_or_ho_hi_z && domain_khi <= zebox.bigEnd(2)    ) ||
         (has_extdir_or_ho_lo_y && domain_jlo >= zebox.smallEnd(1)-1) ||
         (has_extdir_or_ho_hi_y && domain_jhi <= zebox.bigEnd(1)    )  )
    {
        amrex::ParallelFor(zebox, ncomp, [q,wmac,AMREX_D_DECL(domain_ilo,domain_jlo,domain_klo),
                                                 AMREX_D_DECL(domain_ihi,domain_jhi,domain_khi),
                                          Imz,Ipz,dtdz,pbc,flag,vfrac,ccc,AMREX_D_DECL(fcx,fcy,fcz),
                                          is_velocity]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real qpls(0.);
            Real qmns(0.);

            // This means apz(i,j,k) > 0 and we have un-covered cells on both sides
            if (flag(i,j,k).isConnected(0,0,-1))
            {
                const auto& bc = pbc[n];
                bool extdir_or_ho_ilo = (bc.lo(0) == BCType::ext_dir) ||
                                        (bc.lo(0) == BCType::hoextrap);
                bool extdir_or_ho_ihi = (bc.hi(0) == BCType::ext_dir) ||
                                        (bc.hi(0) == BCType::hoextrap);
                bool extdir_or_ho_jlo = (bc.lo(1) == BCType::ext_dir) ||
                                        (bc.lo(1) == BCType::hoextrap);
                bool extdir_or_ho_jhi = (bc.hi(1) == BCType::ext_dir) ||
                                        (bc.hi(1) == BCType::hoextrap);
                bool extdir_or_ho_klo = (bc.lo(2) == BCType::ext_dir) ||
                                        (bc.lo(2) == BCType::hoextrap);
                bool extdir_or_ho_khi = (bc.hi(2) == BCType::ext_dir) ||
                                        (bc.hi(2) == BCType::hoextrap);

                // *************************************************
                // Making qpls
                // *************************************************

                // We have enough cells to do 4th order slopes centered on (i,j,k) with all values at cell centers
                if (vfrac(i,j,k) == 1. && vfrac(i,j,k-1) == 1. && vfrac(i,j,k-2) == 1. &&
                                          vfrac(i,j,k+1) == 1. && vfrac(i,j,k+2) == 1.)
                {
                    int order = 4;
                    qpls = q(i,j,k,n) + Real(0.5) * (-Real(1.0) - wmac(i,j,k) * dtdz) *
                        amrex_calc_zslope_extdir(i,j,k,n,order,q,extdir_or_ho_klo,extdir_or_ho_khi,domain_klo,domain_khi);

                // We have enough cells to do 2nd order slopes with all values at cell centers
                } else if (vfrac(i,j,k) == 1. && vfrac(i,j,k-1) == 1. && vfrac(i,j,k+1) == 1.) {

                    int order = 2;
                    qpls = q(i,j,k,n) + Real(0.5) * (-Real(1.0) - wmac(i,j,k) * dtdz) *
                        amrex_calc_zslope_extdir(i,j,k,n,order,q,extdir_or_ho_klo,extdir_or_ho_khi,domain_klo,domain_khi);

                // We need to use LS slopes
                } else {

                   Real xf = fcz(i,j,k,0); // local (x,y) of centroid of z-face we are extrapolating to
                   Real yf = fcz(i,j,k,1);

                   AMREX_D_TERM(Real delta_z = Real(0.5) + ccc(i,j,k,2);,
                                Real delta_x = xf  - ccc(i,j,k,0);,
                                Real delta_y = yf  - ccc(i,j,k,1););

                   Real qcc_max = amrex::max(q(i,j,k,n), q(i,j,k-1,n));
                   Real qcc_min = amrex::min(q(i,j,k,n), q(i,j,k-1,n));

                   // This will be used in the EB slope routine only if the slope can be computed without LS
                   int max_order = 2;

                   const auto& slopes_eb_hi = amrex_lim_slopes_extdir_eb(i,j,k,n,q,ccc,vfrac,
                                              AMREX_D_DECL(fcx,fcy,fcz), flag,
                                              AMREX_D_DECL(extdir_or_ho_ilo, extdir_or_ho_jlo, extdir_or_ho_klo),
                                              AMREX_D_DECL(extdir_or_ho_ihi, extdir_or_ho_jhi, extdir_or_ho_khi),
                                              AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
                                              AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi),
                                              max_order);

                   qpls = q(i,j,k,n) - delta_z * slopes_eb_hi[2]
                                     + delta_x * slopes_eb_hi[0]
                                     + delta_y * slopes_eb_hi[1];

                   qpls = amrex::max(amrex::min(qpls, qcc_max), qcc_min);
                   qpls -= Real(0.5) * dtdz * wmac(i,j,k) * slopes_eb_hi[2];

                }  // end of making qpls

                // Only over-write normal velocity with Dirichlet bc at lo face
                if (k == domain_klo && (bc.lo(2) == BCType::ext_dir))
                    if (is_velocity && n == 2) qpls = q(i,j,k-1,n);

                // Over-write all with Dirichlet bc at hi face
                if (k == domain_khi+1 && (bc.hi(2) == BCType::ext_dir))
                    qpls = q(i,j,k,n);

                // *************************************************
                // Making qmns
                // *************************************************

                // We have enough cells to do 4th order slopes centered on (i,j,k-1) with all values at cell centers
                if (vfrac(i,j,k-1) == 1. && vfrac(i,j,k-2) == 1. && vfrac(i,j,k-3) == 1. &&
                                            vfrac(i,j,k  ) == 1. && vfrac(i,j,k+1) == 1.)
                {
                    int order = 4;
                    qmns = q(i,j,k-1,n) + Real(0.5) * ( Real(1.0) - wmac(i,j,k) * dtdz) *
                        amrex_calc_zslope_extdir(i,j,k-1,n,order,q,extdir_or_ho_klo,extdir_or_ho_khi,domain_klo,domain_khi);

                // We have enough cells to do 2nd order slopes with all values at cell centers
                } else if (vfrac(i,j,k-1) == 1. && vfrac(i,j,k-2) == 1. && vfrac(i,j,k  ) == 1.)
                {
                    int order = 2;
                    qmns = q(i,j,k-1,n) + Real(0.5) * ( Real(1.0) - wmac(i,j,k) * dtdz) *
                        amrex_calc_zslope_extdir(i,j,k-1,n,order,q,extdir_or_ho_klo,extdir_or_ho_khi,domain_klo,domain_khi);

                // We need to use LS slopes
                } else {

                   Real xf = fcz(i,j,k,0); // local (x,y) of centroid of z-face we are extrapolating to
                   Real yf = fcz(i,j,k,1);

                   AMREX_D_TERM(Real delta_z = Real(0.5) - ccc(i,j,k-1,2);,
                                Real delta_x = xf  - ccc(i,j,k-1,0);,
                                Real delta_y = yf  - ccc(i,j,k-1,1););

                   Real qcc_max = amrex::max(q(i,j,k,n), q(i,j,k-1,n));
                   Real qcc_min = amrex::min(q(i,j,k,n), q(i,j,k-1,n));

                   // This will be used in the EB slope routine only if the slope can be computed without LS
                   int max_order = 2;

                   const auto& slopes_eb_lo = amrex_lim_slopes_extdir_eb(i,j,k-1,n,q,ccc,vfrac,
                                              AMREX_D_DECL(fcx,fcy,fcz), flag,
                                              AMREX_D_DECL(extdir_or_ho_ilo, extdir_or_ho_jlo, extdir_or_ho_klo),
                                              AMREX_D_DECL(extdir_or_ho_ihi, extdir_or_ho_jhi, extdir_or_ho_khi),
                                              AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
                                              AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi),
                                              max_order);


                   qmns = q(i,j,k-1,n) + delta_x * slopes_eb_lo[0]
                                       + delta_y * slopes_eb_lo[1]
                                       + delta_z * slopes_eb_lo[2];

                   qmns = amrex::max(amrex::min(qmns, qcc_max), qcc_min);
                   qmns -= Real(0.5) * dtdz * wmac(i,j,k) * slopes_eb_lo[2];

                }  // end of making qmns

                // Over-write all with Dirichlet bc at lo face
                if (k == domain_klo && (bc.lo(2) == BCType::ext_dir))
                    qmns = q(i,j,k-1,n);

                // Only over-write normal velocity with Dirichlet bc at hi face
                if (k == domain_khi+1 && (bc.hi(2) == BCType::ext_dir))
                    if (is_velocity && n == 2) qmns = q(i,j,k,n);
            }

            Ipz(i,j,k-1,n) = qmns;
            Imz(i,j,k  ,n) = qpls;
        });
    }
    else // The cases below are not near any domain boundary
    {
        amrex::ParallelFor(zebox, ncomp, [q,wmac,AMREX_D_DECL(domain_ilo,domain_jlo,domain_klo),
                                                 AMREX_D_DECL(domain_ihi,domain_jhi,domain_khi),
                                          Imz,Ipz,dt,dtdz,pbc,flag,vfrac,ccc,AMREX_D_DECL(fcx,fcy,fcz),
                                          is_velocity]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real qpls(0.);
            Real qmns(0.);

            // This means apz(i,j,k) > 0 and we have un-covered cells on both sides
            if (flag(i,j,k).isConnected(0,0,-1))
            {
                // *************************************************
                // Making qpls
                // *************************************************

                // We have enough cells to do 4th order slopes centered on (i,j,k) with all values at cell centers
                if (vfrac(i,j,k) == 1. && vfrac(i,j,k-1) == 1. && vfrac(i,j,k-2) == 1. &&
                                          vfrac(i,j,k+1) == 1. && vfrac(i,j,k+2) == 1.)
                {
                    int order = 4;
                    qpls = q(i,j,k,n) + Real(0.5) * (-Real(1.0) - wmac(i,j,k) * dtdz) *
                        amrex_calc_zslope(i,j,k,n,order,q);

                // We have enough cells to do 2nd order slopes with all values at cell centers
                } else if (vfrac(i,j,k) == 1. && vfrac(i,j,k-1) == 1. && vfrac(i,j,k+1) == 1.) {

                    int order = 2;
                    qpls = q(i,j,k,n) + Real(0.5) * (-Real(1.0) - wmac(i,j,k) * dtdz) *
                        amrex_calc_zslope(i,j,k,n,order,q);

                // We need to use LS slopes
                } else {

                   Real xf = fcz(i,j,k,0); // local (x,y) of centroid of z-face we are extrapolating to
                   Real yf = fcz(i,j,k,1);

                   AMREX_D_TERM(Real delta_z = Real(0.5) + ccc(i,j,k,2);,
                                Real delta_x = xf  - ccc(i,j,k,0);,
                                Real delta_y = yf  - ccc(i,j,k,1););

                   Real qcc_max = amrex::max(q(i,j,k,n), q(i,j,k-1,n));
                   Real qcc_min = amrex::min(q(i,j,k,n), q(i,j,k-1,n));

                   // This will be used in the EB slope routine only if the slope can be computed without LS
                   int max_order = 2;

                   const auto& slopes_eb_hi = amrex_lim_slopes_eb(i,j,k,n,q,ccc,vfrac,
                                                                  AMREX_D_DECL(fcx,fcy,fcz),flag,max_order);

                   qpls = q(i,j,k,n) - delta_z * slopes_eb_hi[2]
                                     + delta_x * slopes_eb_hi[0]
                                     + delta_y * slopes_eb_hi[1];

                   qpls = amrex::max(amrex::min(qpls, qcc_max), qcc_min);
                   qpls -= Real(0.5) * dtdz * wmac(i,j,k) * slopes_eb_hi[2];

                }  // end of making qpls

                // *************************************************
                // Making qmns
                // *************************************************

                // We have enough cells to do 4th order slopes centered on (i,j,k-1) with all values at cell centers
                if (vfrac(i,j,k-1) == 1. && vfrac(i,j,k-2) == 1. && vfrac(i,j,k-3) == 1. &&
                                            vfrac(i,j,k  ) == 1. && vfrac(i,j,k+1) == 1.)
                {
                    int order = 4;
                    qmns = q(i,j,k-1,n) + Real(0.5) * ( Real(1.0) - wmac(i,j,k) * dtdz) *
                        amrex_calc_zslope(i,j,k-1,n,order,q);

                // We have enough cells to do 2nd order slopes with all values at cell centers
                } else if (vfrac(i,j,k-1) == 1. && vfrac(i,j,k-2) == 1. && vfrac(i,j,k  ) == 1.)
                {
                    int order = 2;
                    qmns = q(i,j,k-1,n) + Real(0.5) * ( Real(1.0) - wmac(i,j,k) * dtdz) *
                        amrex_calc_zslope(i,j,k-1,n,order,q);

                // We need to use LS slopes
                } else {

                   Real xf = fcz(i,j,k,0); // local (x,y) of centroid of z-face we are extrapolating to
                   Real yf = fcz(i,j,k,1);

                   AMREX_D_TERM(Real delta_z = Real(0.5) - ccc(i,j,k-1,2);,
                                Real delta_x = xf  - ccc(i,j,k-1,0);,
                                Real delta_y = yf  - ccc(i,j,k-1,1););

                   Real qcc_max = amrex::max(q(i,j,k,n), q(i,j,k-1,n));
                   Real qcc_min = amrex::min(q(i,j,k,n), q(i,j,k-1,n));

                   // This will be used in the EB slope routine only if the slope can be computed without LS
                   int max_order = 2;

                   const auto& slopes_eb_lo = amrex_lim_slopes_eb(i,j,k-1,n,q,ccc,vfrac,
                                                                  AMREX_D_DECL(fcx,fcy,fcz),flag,max_order);

                   qmns = q(i,j,k-1,n) + delta_x * slopes_eb_lo[0]
                                       + delta_y * slopes_eb_lo[1]
                                       + delta_z * slopes_eb_lo[2];

                   qmns = amrex::max(amrex::min(qmns, qcc_max), qcc_min);
                   qmns -= Real(0.5) * dtdz * wmac(i,j,k) * slopes_eb_lo[2];

                }  // end of making qmns
            }

            Ipz(i,j,k-1,n) = qmns;
            Imz(i,j,k  ,n) = qpls;
        });
    }
}

#endif
/** @} */

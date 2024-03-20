/**
 * \file hydro_ebgodunov_edge_state_3D.cpp
 * \addtogroup EBGodunov
 * @{
 *
 */

#include <hydro_bcs_K.H>
#include <hydro_ebgodunov.H>
#include <hydro_ebgodunov_plm.H>
#include <hydro_ebgodunov_corner_couple.H>

using namespace amrex;

void
EBGodunov::ComputeEdgeState ( Box const& bx, int ncomp,
                              Array4<Real const> const& q,
                              Array4<Real> const& xedge,
                              Array4<Real> const& yedge,
                              Array4<Real> const& zedge,
                              Array4<Real const> const& u_mac,
                              Array4<Real const> const& v_mac,
                              Array4<Real const> const& w_mac,
                              Array4<Real const> const& divu,
                              Array4<Real const> const& fq,
                              Geometry const& geom,
                              Real l_dt,
                              Vector<BCRec> const& h_bcrec,
                              BCRec const*  pbc,
                              int const* iconserv,
                              Real* p,
                              Array4<EBCellFlag const> const& flag_arr,
                              Array4<Real const> const& apx,
                              Array4<Real const> const& apy,
                              Array4<Real const> const& apz,
                              Array4<Real const> const& vfrac_arr,
                              Array4<Real const> const& fcx,
                              Array4<Real const> const& fcy,
                              Array4<Real const> const& fcz,
                              Array4<Real const> const& ccent_arr,
                              bool is_velocity,
                              Array4<Real const> const& values_on_eb_inflow,
                              Array4<int const> const& bc_arr)
{

    // bx is the cell-centered box on which we want to compute the advective update
    Box const& xbx = amrex::surroundingNodes(bx,0);
    Box const& ybx = amrex::surroundingNodes(bx,1);
    Box const& zbx = amrex::surroundingNodes(bx,2);

    // Start with above and grow 1 tangentially
    Box xebx = Box(xbx).grow(1,1).grow(2,1);
    Box yebx = Box(ybx).grow(0,1).grow(2,1);
    Box zebx = Box(zbx).grow(0,1).grow(1,1);

    Box const& bxg2 = amrex::grow(bx,2);

    const Real dx = geom.CellSize(0);
    const Real dy = geom.CellSize(1);
    const Real dz = geom.CellSize(2);
    Real dtdx = l_dt/dx;
    Real dtdy = l_dt/dy;
    Real dtdz = l_dt/dz;

    Box const& domain = geom.Domain();
    const auto dlo = amrex::lbound(domain);
    const auto dhi = amrex::ubound(domain);

    Array4<Real> Imx = makeArray4(p, bxg2, ncomp);
    p +=         Imx.size();
    Array4<Real> Ipx = makeArray4(p, bxg2, ncomp);
    p +=         Ipx.size();
    Array4<Real> Imy = makeArray4(p, bxg2, ncomp);
    p +=         Imy.size();
    Array4<Real> Ipy = makeArray4(p, bxg2, ncomp);
    p +=         Ipy.size();
    Array4<Real> Imz = makeArray4(p, bxg2, ncomp);
    p +=         Imz.size();
    Array4<Real> Ipz = makeArray4(p, bxg2, ncomp);
    p +=         Ipz.size();

    Array4<Real> xlo = makeArray4(p, xebx, ncomp);
    p +=         xlo.size();
    Array4<Real> xhi = makeArray4(p, xebx, ncomp);
    p +=         xhi.size();
    Array4<Real> ylo = makeArray4(p, yebx, ncomp);
    p +=         ylo.size();
    Array4<Real> yhi = makeArray4(p, yebx, ncomp);
    p +=         yhi.size();
    Array4<Real> zlo = makeArray4(p, zebx, ncomp);
    p +=         zlo.size();
    Array4<Real> zhi = makeArray4(p, zebx, ncomp);
    p +=         zhi.size();

    Array4<Real> xyzlo = makeArray4(p, bxg2, ncomp);
    p +=         xyzlo.size();
    Array4<Real> xyzhi = makeArray4(p, bxg2, ncomp);

    // Initialize this way out of an abundance of paranoia
    amrex::ParallelFor(
        Box(Imx), ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
           Imx(i,j,k,n) = i*Real(1.e10) + j*Real(1.e20) + k*Real(1.30) + n*Real(1.e4);
        });
    amrex::ParallelFor(
        Box(Imy), ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
           Imy(i,j,k,n) = i*Real(1.e10) + j*Real(1.e20) + k*Real(1.30) + n*Real(1.e4);
        });
    amrex::ParallelFor(
        Box(Imz), ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
           Imz(i,j,k,n) = i*Real(1.e10) + j*Real(1.e20) + k*Real(1.30) + n*Real(1.e4);
        });

    EBPLM::PredictStateOnXFace( xebx, ncomp, Imx, Ipx, q, u_mac,
                                flag_arr, vfrac_arr,
                                AMREX_D_DECL(fcx,fcy,fcz),ccent_arr,
                                geom, l_dt, h_bcrec, pbc, is_velocity, bc_arr);

    EBPLM::PredictStateOnYFace( yebx, ncomp, Imy, Ipy, q, v_mac,
                                flag_arr, vfrac_arr,
                                AMREX_D_DECL(fcx,fcy,fcz),ccent_arr,
                                geom, l_dt, h_bcrec, pbc, is_velocity, bc_arr);

    EBPLM::PredictStateOnZFace( zebx, ncomp, Imz, Ipz, q, w_mac,
                                flag_arr, vfrac_arr,
                                AMREX_D_DECL(fcx,fcy,fcz),ccent_arr,
                                geom, l_dt, h_bcrec, pbc, is_velocity, bc_arr);

    amrex::ParallelFor(
        xebx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real lo = Ipx(i-1,j,k,n);
            Real hi = Imx(i  ,j,k,n);

            Real uad = u_mac(i,j,k);

            const auto bc = HydroBC::getBC(i, j, k, n, domain, pbc, bc_arr);

            HydroBC::SetXEdgeBCs(i, j, k, n, q, lo, hi, bc.lo(0), dlo.x, bc.hi(0), dhi.x, is_velocity);

            xlo(i,j,k,n) = lo;
            xhi(i,j,k,n) = hi;

            Real st = (uad >= 0.) ? lo : hi;
            Real fux = (amrex::Math::abs(uad) < small_vel)? Real(0.0) : Real(1.0);
            Imx(i,j,k,n) = fux*st + (Real(1.0) - fux)*Real(0.5)*(hi + lo);
        },
        yebx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real lo = Ipy(i,j-1,k,n);
            Real hi = Imy(i,j  ,k,n);

            Real vad = v_mac(i,j,k);

            const auto bc = HydroBC::getBC(i, j, k, n, domain, pbc, bc_arr);

            HydroBC::SetYEdgeBCs(i, j, k, n, q, lo, hi, bc.lo(1), dlo.y, bc.hi(1), dhi.y, is_velocity);

            ylo(i,j,k,n) = lo;
            yhi(i,j,k,n) = hi;

            Real st = (vad >= 0.) ? lo : hi;
            Real fuy = (amrex::Math::abs(vad) < small_vel)? Real(0.0) : Real(1.0);
            Imy(i,j,k,n) = fuy*st + (Real(1.0) - fuy)*Real(0.5)*(hi + lo);
        },
        zebx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real lo = Ipz(i,j,k-1,n);
            Real hi = Imz(i,j,k  ,n);

            Real wad = w_mac(i,j,k);

            const auto bc = HydroBC::getBC(i, j, k, n, domain, pbc, bc_arr);

            HydroBC::SetZEdgeBCs(i, j, k, n, q, lo, hi, bc.lo(2), dlo.z, bc.hi(2), dhi.z, is_velocity);

            zlo(i,j,k,n) = lo;
            zhi(i,j,k,n) = hi;

            Real st = (wad >= 0.) ? lo : hi;
            Real fuz = (amrex::Math::abs(wad) < small_vel) ? 0. : 1.;
            Imz(i,j,k,n) = fuz*st + (Real(1.0) - fuz)*Real(0.5)*(hi + lo);
        });

    // We can reuse the space in Ipx, Ipy and Ipz.
    Array4<Real> xed = Imx;
    Array4<Real> yed = Imy;
    Array4<Real> zed = Imz;

    //
    // Grow in x-direction
    //
    Box const& xbxtmp = amrex::grow(bx,0,1);
    Array4<Real> yzlo = makeArray4(xyzlo.dataPtr(), amrex::surroundingNodes(xbxtmp,1), ncomp);
    Array4<Real> zylo = makeArray4(xyzhi.dataPtr(), amrex::surroundingNodes(xbxtmp,2), ncomp);

    amrex::ParallelFor(
    Box(zylo), ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        const auto bc = HydroBC::getBC(i, j, k, n, domain, pbc, bc_arr);
        Real l_zylo, l_zyhi;
        EBGodunovCornerCouple::EBGodunov_corner_couple_zy(l_zylo, l_zyhi,
                                   i, j, k, n, l_dt, dy, iconserv[n],
                                   zlo(i,j,k,n), zhi(i,j,k,n),
                                   q, divu, apx, apy, apz, vfrac_arr, v_mac, yed);

        Real wad = w_mac(i,j,k);
        HydroBC::SetZEdgeBCs(i, j, k, n, q, l_zylo, l_zyhi, bc.lo(2), dlo.z, bc.hi(2), dhi.z, is_velocity);

        Real st = (wad >= 0.) ? l_zylo : l_zyhi;
        Real fu = (amrex::Math::abs(wad) < small_vel) ? Real(0.0) : Real(1.0);
        zylo(i,j,k,n) = fu*st + (Real(1.0) - fu) * Real(0.5) * (l_zyhi + l_zylo);
    },
    Box(yzlo), ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        const auto bc = HydroBC::getBC(i, j, k, n, domain, pbc, bc_arr);
        Real l_yzlo, l_yzhi;
        EBGodunovCornerCouple::EBGodunov_corner_couple_yz(l_yzlo, l_yzhi,
                                   i, j, k, n, l_dt, dz, iconserv[n],
                                   ylo(i,j,k,n), yhi(i,j,k,n),
                                   q, divu, apx, apy, apz, vfrac_arr, w_mac, zed);

        Real vad = v_mac(i,j,k);
        HydroBC::SetYEdgeBCs(i, j, k, n, q, l_yzlo, l_yzhi, bc.lo(1), dlo.y, bc.hi(1), dhi.y, is_velocity);

        Real st = (vad >= 0.) ? l_yzlo : l_yzhi;
        Real fu = (amrex::Math::abs(vad) < small_vel) ? Real(0.0) : Real(1.0);
        yzlo(i,j,k,n) = fu*st + (Real(1.0) - fu) * Real(0.5) * (l_yzhi + l_yzlo);
    });
    //

    amrex::ParallelFor(xbx, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (apx(i,j,k) > 0.)
        {
            Real stl = xlo(i,j,k,n);
            Real sth = xhi(i,j,k,n);

#ifdef AMREX_USE_MOVING_EB
            const int no_eb_flow_xlo = !(values_on_eb_inflow) ? 1 : 0;
#else
            const int no_eb_flow_xlo = !(values_on_eb_inflow) ? 1 :
                ((Math::abs(values_on_eb_inflow(i  ,j,k,n)) > 0. ||
                  Math::abs(values_on_eb_inflow(i-1,j,k,n)) > 0.) ? 0 : 1);
#endif

            // If we can't compute good transverse terms, don't use any d/dt terms at all
            if (apy(i-1,j+1,k) > 0. && apy(i-1,j  ,k) > 0. && apz(i-1,j,k+1) > 0. && apz(i-1,j,k) > 0. && no_eb_flow_xlo)
            {
                // Here we add  dt/2 (-q u_x - (v q)_y - (w q)_z) to the term that is already
                //     q + dx/2 q_x + dt/2 (-u q_x) to get
                //     q + dx/2 q_x - dt/2 (u q_x  + q u_x + (v q)_y + (w q)_z) which is equivalent to
                // --> q + dx/2 q_x - dt/2 ( div (uvec q) )
                Real quxl = (apx(i,j,k)*u_mac(i,j,k) - apx(i-1,j,k)*u_mac(i-1,j,k)) * q(i-1,j,k,n);

                stl += ( - (Real(0.5)*dtdx) * quxl
                         - (Real(0.5)*dtdy)*(apy(i-1,j+1,k  )*yzlo(i-1,j+1,k  ,n)*v_mac(i-1,j+1,k  )
                                     - apy(i-1,j  ,k  )*yzlo(i-1,j  ,k  ,n)*v_mac(i-1,j  ,k  ))
                         - (Real(0.5)*dtdz)*(apz(i-1,j  ,k+1)*zylo(i-1,j  ,k+1,n)*w_mac(i-1,j  ,k+1)
                                     - apz(i-1,j  ,k  )*zylo(i-1,j  ,k  ,n)*w_mac(i-1,j  ,k  )) ) / vfrac_arr(i-1,j,k);

                // Here we adjust for non-conservative by removing the q divu contribution to get
                //     q + dx/2 q_x - dt/2 ( div (uvec q) - q divu ) which is equivalent to
                // --> q + dx/2 q_x - dt/2 ( uvec dot grad q)
                stl += (!iconserv[n]) ? Real(0.5)*l_dt* q(i-1,j,k,n)*divu(i-1,j,k) : Real(0.0);

                stl += (fq)           ? Real(0.5)*l_dt*fq(i-1,j,k,n) : Real(0.0);
            }

#ifdef AMREX_USE_MOVING_EB
            const int no_eb_flow_xhi = !(values_on_eb_inflow) ? 1 : 0;
#else
            const int no_eb_flow_xhi = !(values_on_eb_inflow) ? 1 :
                ((Math::abs(values_on_eb_inflow(i+1,j,k,n)) > 0. ||
                  Math::abs(values_on_eb_inflow(i  ,j,k,n)) > 0.) ? 0 : 1);
#endif

            // If we can't compute good transverse terms, don't use any d/dt terms at all
            if (apy(i,j+1,k) > 0. && apy(i,j  ,k) > 0. && apz(i,j,k+1) > 0. && apz(i,j,k) > 0. && no_eb_flow_xhi)
            {
                // Here we add  dt/2 (-q u_x - (v q)_y - (w q)_z) to the term that is already
                //     q + dx/2 q_x + dt/2 (-u q_x) to get
                //     q + dx/2 q_x - dt/2 (u q_x  + q u_x + (v q)_y + (w q)_z) which is equivalent to
                // --> q + dx/2 q_x - dt/2 ( div (uvec q) )
                Real quxh = (apx(i+1,j,k)*u_mac(i+1,j,k) - apx(i,j,k)*u_mac(i,j,k)) * q(i,j,k,n);

                sth += ( - (Real(0.5)*dtdx) * quxh
                         - (Real(0.5)*dtdy)*(apy(i,j+1,k  )*yzlo(i,j+1,k  ,n)*v_mac(i,j+1,k  )
                                     - apy(i,j  ,k  )*yzlo(i,j  ,k  ,n)*v_mac(i,j  ,k  ))
                         - (Real(0.5)*dtdz)*(apz(i,j  ,k+1)*zylo(i,j  ,k+1,n)*w_mac(i,j  ,k+1)
                                     - apz(i,j  ,k  )*zylo(i,j  ,k  ,n)*w_mac(i,j  ,k  )) ) / vfrac_arr(i,j,k);

                // Here we adjust for non-conservative by removing the q divu contribution to get
                //     q + dx/2 q_x - dt/2 ( div (uvec q) - q divu ) which is equivalent to
                // --> q + dx/2 q_x - dt/2 ( uvec dot grad q)
                sth += (!iconserv[n]) ? Real(0.5)*l_dt* q(i  ,j,k,n)*divu(i,j,k) : Real(0.0);

                sth += (fq)           ? Real(0.5)*l_dt*fq(i  ,j,k,n) : Real(0.0);
            }

            const auto bc = HydroBC::getBC(i, j, k, n, domain, pbc, bc_arr);
            HydroBC::SetXEdgeBCs(i, j, k, n, q, stl, sth, bc.lo(0), dlo.x, bc.hi(0), dhi.x, is_velocity);

            if ( (i==dlo.x) && (bc.lo(0) == BCType::foextrap || bc.lo(0) == BCType::hoextrap) )
            {
                if ( u_mac(i,j,k) >= 0. && n==XVEL && is_velocity )  sth = amrex::min(sth,0.0_rt);
                stl = sth;
            }
            if ( (i==dhi.x+1) && (bc.hi(0) == BCType::foextrap || bc.hi(0) == BCType::hoextrap) )
            {
                if ( u_mac(i,j,k) <= 0. && n==XVEL && is_velocity ) stl = amrex::max(stl,0.0_rt);
                sth = stl;
            }

            Real temp = (u_mac(i,j,k) >= 0.) ? stl : sth;
            temp = (amrex::Math::abs(u_mac(i,j,k)) < small_vel) ? Real(0.5)*(stl + sth) : temp;
            xedge(i,j,k,n) = temp;

        } else {
           xedge(i,j,k,n) = Real(0.0);
        }
    });

    //
    // Grow in y-direction
    //
    Box const& ybxtmp = amrex::grow(bx,1,1);
    Array4<Real> xzlo = makeArray4(xyzlo.dataPtr(), amrex::surroundingNodes(ybxtmp,0), ncomp);
    Array4<Real> zxlo = makeArray4(xyzhi.dataPtr(), amrex::surroundingNodes(ybxtmp,2), ncomp);

    amrex::ParallelFor(
    Box(xzlo), ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        const auto bc = HydroBC::getBC(i, j, k, n, domain, pbc, bc_arr);
        Real l_xzlo, l_xzhi;
        EBGodunovCornerCouple::EBGodunov_corner_couple_xz(l_xzlo, l_xzhi,
                                   i, j, k, n, l_dt, dz, iconserv[n],
                                   xlo(i,j,k,n),  xhi(i,j,k,n),
                                   q, divu, apx, apy, apz, vfrac_arr, w_mac, zed);

        Real uad = u_mac(i,j,k);
        HydroBC::SetXEdgeBCs(i, j, k, n, q, l_xzlo, l_xzhi, bc.lo(0), dlo.x, bc.hi(0), dhi.x, is_velocity);

        Real st = (uad >= 0.) ? l_xzlo : l_xzhi;
        Real fu = (amrex::Math::abs(uad) < small_vel) ? 0.0 : 1.0;
        xzlo(i,j,k,n) = fu*st + (Real(1.0) - fu) * Real(0.5) * (l_xzhi + l_xzlo);
    },
    Box(zxlo), ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        const auto bc = HydroBC::getBC(i, j, k, n, domain, pbc, bc_arr);
        Real l_zxlo, l_zxhi;
        EBGodunovCornerCouple::EBGodunov_corner_couple_zx(l_zxlo, l_zxhi,
                                   i, j, k, n, l_dt, dx, iconserv[n],
                                   zlo(i,j,k,n), zhi(i,j,k,n),
                                   q, divu, apx, apy, apz, vfrac_arr, u_mac, xed);

        Real wad = w_mac(i,j,k);
        HydroBC::SetZEdgeBCs(i, j, k, n, q, l_zxlo, l_zxhi, bc.lo(2), dlo.z, bc.hi(2), dhi.z, is_velocity);

        Real st = (wad >= 0.) ? l_zxlo : l_zxhi;
        Real fu = (amrex::Math::abs(wad) < small_vel) ? 0.0 : 1.0;
        zxlo(i,j,k,n) = fu*st + (Real(1.0) - fu) * Real(0.5) * (l_zxhi + l_zxlo);
    });
    //

    amrex::ParallelFor(ybx, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (apy(i,j,k) > 0.)
        {
            Real stl = ylo(i,j,k,n);
            Real sth = yhi(i,j,k,n);

#ifdef AMREX_USE_MOVING_EB
            const int no_eb_flow_ylo = !(values_on_eb_inflow) ? 1 : 0;
#else
            const int no_eb_flow_ylo = !(values_on_eb_inflow) ? 1 :
                ((Math::abs(values_on_eb_inflow(i,j  ,k,n)) > 0. ||
                  Math::abs(values_on_eb_inflow(i,j-1,k,n)) > 0.) ? 0 : 1);
#endif

            // If we can't compute good transverse terms, don't use any d/dt terms at all
            if (apx(i+1,j-1,k) > 0. && apx(i,j-1,k) > 0. && apz(i,j-1,k+1) > 0. && apz(i,j-1,k) > 0. && no_eb_flow_ylo)
            {
                // Here we add  dt/2 (-q v_y - (u q)_x - (w q)_z) to the term that is already
                //     q + dy/2 q_y + dt/2 (-v q_y) to get
                //     q + dy/2 q_y - dt/2 (v q_y  + q v_y + (u q)_x + (w q)_z) which is equivalent to
                // --> q + dy/2 q_y - dt/2 ( div (uvec q) )
                Real qvyl = (apy(i,j,k)*v_mac(i,j,k) - apy(i,j-1,k)*v_mac(i,j-1,k)) * q(i,j-1,k,n);

                stl += ( - (Real(0.5)*dtdy) * qvyl
                         - (Real(0.5)*dtdx)*(apx(i+1,j-1,k  )*xzlo(i+1,j-1,k  ,n)*u_mac(i+1,j-1,k  )
                                     - apx(i  ,j-1,k  )*xzlo(i  ,j-1,k  ,n)*u_mac(i  ,j-1,k  ))
                         - (Real(0.5)*dtdz)*(apz(i  ,j-1,k+1)*zxlo(i  ,j-1,k+1,n)*w_mac(i  ,j-1,k+1)
                                     - apz(i  ,j-1,k  )*zxlo(i  ,j-1,k  ,n)*w_mac(i  ,j-1,k  )) ) / vfrac_arr(i,j-1,k);

                // Here we adjust for non-conservative by removing the q divu contribution to get
                //     q + dy/2 q_y - dt/2 ( div (uvec q) - q divu ) which is equivalent to
                // --> q + dy/2 q_y - dt/2 ( uvec dot grad q)
                stl += (!iconserv[n]) ? Real(0.5)*l_dt* q(i,j-1,k,n)*divu(i,j-1,k) : Real(0.0);

                stl += (fq)           ? Real(0.5)*l_dt*fq(i,j-1,k,n) : Real(0.0);
            }

#ifdef AMREX_USE_MOVING_EB
            const int no_eb_flow_yhi = !(values_on_eb_inflow) ? 1 : 0;
#else
            const int no_eb_flow_yhi = !(values_on_eb_inflow) ? 1 :
                ((Math::abs(values_on_eb_inflow(i,j+1,k,n)) > 0. ||
                  Math::abs(values_on_eb_inflow(i,j  ,k,n)) > 0.) ? 0 : 1);
#endif

            // If we can't compute good transverse terms, don't use any d/dt terms at all
            if (apx(i+1,j,k) > 0. && apx(i,j,k) > 0. && apz(i,j,k+1) > 0. && apz(i,j,k) > 0. && no_eb_flow_yhi)
            {
                // Here we add  dt/2 (-q v_y - (u q)_x - (w q)_z) to the term that is already
                //     q + dy/2 q_y + dt/2 (-v q_y) to get
                //     q + dy/2 q_y - dt/2 (v q_y  + q v_y + (u q)_x + (w q)_z) which is equivalent to
                // --> q + dy/2 q_y - dt/2 ( div (uvec q) )
                Real qvyh = (apy(i,j+1,k)*v_mac(i,j+1,k) - apy(i,j,k)*v_mac(i,j,k)) * q(i,j,k,n);

                sth += ( - (Real(0.5)*dtdy) * qvyh
                         - (Real(0.5)*dtdx)*(apx(i+1,j,k  )*xzlo(i+1,j,k  ,n)*u_mac(i+1,j,k  )
                                     - apx(i  ,j,k  )*xzlo(i  ,j,k  ,n)*u_mac(i  ,j,k  ))
                         - (Real(0.5)*dtdz)*(apz(i  ,j,k+1)*zxlo(i  ,j,k+1,n)*w_mac(i  ,j,k+1)
                                     - apz(i  ,j,k  )*zxlo(i  ,j,k  ,n)*w_mac(i  ,j,k  )) ) / vfrac_arr(i,j,k);

                // Here we adjust for non-conservative by removing the q divu contribution to get
                //     q + dy/2 q_y - dt/2 ( div (uvec q) - q divu ) which is equivalent to
                // --> q + dy/2 q_y - dt/2 ( uvec dot grad q)
                sth += (!iconserv[n]) ? Real(0.5)*l_dt* q(i,j,k,n)*divu(i,j,k) : Real(0.0);

                sth += (fq)           ? Real(0.5)*l_dt*fq(i,j,k,n) : Real(0.0);
            }

            const auto bc = HydroBC::getBC(i, j, k, n, domain, pbc, bc_arr);
            HydroBC::SetYEdgeBCs(i, j, k, n, q, stl, sth, bc.lo(1), dlo.y, bc.hi(1), dhi.y, is_velocity);

            if ( (j==dlo.y) && (bc.lo(1) == BCType::foextrap || bc.lo(1) == BCType::hoextrap) )
            {
                if ( v_mac(i,j,k) >= 0. && n==YVEL && is_velocity ) sth = amrex::min(sth,0.0_rt);
                stl = sth;
            }
            if ( (j==dhi.y+1) && (bc.hi(1) == BCType::foextrap || bc.hi(1) == BCType::hoextrap) )
            {
                if ( v_mac(i,j,k) <= 0. && n==YVEL && is_velocity ) stl = amrex::max(stl,0.0_rt);
                sth = stl;
            }
            Real temp = (v_mac(i,j,k) >= 0.) ? stl : sth;
            temp = (amrex::Math::abs(v_mac(i,j,k)) < small_vel) ? Real(0.5)*(stl + sth) : temp;
            yedge(i,j,k,n) = temp;

        } else {
            yedge(i,j,k,n) = Real(0.0);
        }
    });

    //
    // z-direction
    //
    Box const& zbxtmp = amrex::grow(bx,2,1);
    Array4<Real> xylo = makeArray4(xyzlo.dataPtr(), amrex::surroundingNodes(zbxtmp,0), ncomp);
    Array4<Real> yxlo = makeArray4(xyzhi.dataPtr(), amrex::surroundingNodes(zbxtmp,1), ncomp);
    amrex::ParallelFor(
    Box(xylo), ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        const auto bc = HydroBC::getBC(i, j, k, n, domain, pbc, bc_arr);
        Real l_xylo, l_xyhi;
        EBGodunovCornerCouple::EBGodunov_corner_couple_xy(l_xylo, l_xyhi,
                                   i, j, k, n, l_dt, dy, iconserv[n],
                                   xlo(i,j,k,n), xhi(i,j,k,n),
                                   q, divu, apx, apy, apz, vfrac_arr, v_mac, yed);

        Real uad = u_mac(i,j,k);
        HydroBC::SetXEdgeBCs(i, j, k, n, q, l_xylo, l_xyhi, bc.lo(0), dlo.x, bc.hi(0), dhi.x, is_velocity);

        Real st = (uad >= 0.) ? l_xylo : l_xyhi;
        Real fu = (amrex::Math::abs(uad) < small_vel) ? 0.0 : 1.0;
        xylo(i,j,k,n) = fu*st + (Real(1.0) - fu) * Real(0.5) * (l_xyhi + l_xylo);
    },
    Box(yxlo), ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        const auto bc = HydroBC::getBC(i, j, k, n, domain, pbc, bc_arr);
        Real l_yxlo, l_yxhi;
        EBGodunovCornerCouple::EBGodunov_corner_couple_yx(l_yxlo, l_yxhi,
                                   i, j, k, n, l_dt, dx, iconserv[n],
                                   ylo(i,j,k,n), yhi(i,j,k,n),
                                   q, divu, apx, apy, apz, vfrac_arr, u_mac, xed);

        Real vad = v_mac(i,j,k);
        HydroBC::SetYEdgeBCs(i, j, k, n, q, l_yxlo, l_yxhi, bc.lo(1), dlo.y, bc.hi(1), dhi.y, is_velocity);

        Real st = (vad >= 0.) ? l_yxlo : l_yxhi;
        Real fu = (amrex::Math::abs(vad) < small_vel) ? 0.0 : 1.0;
        yxlo(i,j,k,n) = fu*st + (Real(1.0) - fu) * Real(0.5) * (l_yxhi + l_yxlo);
    });
    //
    amrex::ParallelFor(zbx, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (apz(i,j,k) > 0.)
        {
            Real stl = zlo(i,j,k,n);
            Real sth = zhi(i,j,k,n);

#ifdef AMREX_USE_MOVING_EB
            const int no_eb_flow_zlo = !(values_on_eb_inflow) ? 1 : 0;
#else
            const int no_eb_flow_zlo = !(values_on_eb_inflow) ? 1 :
                ((Math::abs(values_on_eb_inflow(i,j,k  ,n)) > 0. ||
                  Math::abs(values_on_eb_inflow(i,j,k-1,n)) > 0.) ? 0 : 1);
#endif

            // If we can't compute good transverse terms, don't use any d/dt terms at all
            if (apx(i+1,j,k-1) > 0. && apx(i,j,k-1) > 0. && apy(i,j+1,k-1) > 0. && apy(i,j,k-1) > 0. && no_eb_flow_zlo)
            {
                // Here we add  dt/2 (-q w_z - (u q)_x - (v q)_y) to the term that is already
                //     q + dz/2 q_z + dt/2 (-w q_z) to get
                //     q + dz/2 q_z - dt/2 (w q_z  + q w_z + (u q)_x + (v q)_y) which is equivalent to
                // --> q + dz/2 q_z - dt/2 ( div (uvec q) )
                Real qwzl = (apz(i,j,k)*w_mac(i,j,k) - apz(i,j,k-1)*w_mac(i,j,k-1)) * q(i,j,k-1,n);

                stl += ( - (Real(0.5)*dtdz) * qwzl
                         - (Real(0.5)*dtdx)*(apx(i+1,j  ,k-1)*xylo(i+1,j  ,k-1,n)*u_mac(i+1,j  ,k-1)
                                      -apx(i  ,j  ,k-1)*xylo(i  ,j  ,k-1,n)*u_mac(i  ,j  ,k-1))
                         - (Real(0.5)*dtdy)*(apy(i  ,j+1,k-1)*yxlo(i  ,j+1,k-1,n)*v_mac(i  ,j+1,k-1)
                                      -apy(i  ,j  ,k-1)*yxlo(i  ,j  ,k-1,n)*v_mac(i  ,j  ,k-1)) ) / vfrac_arr(i,j,k-1);

                // Here we adjust for non-conservative by removing the q divu contribution to get
                //     q + dz/2 q_z - dt/2 ( div (uvec q) - q divu ) which is equivalent to
                // --> q + dz/2 q_z - dt/2 ( uvec dot grad q)
                stl += (!iconserv[n]) ? Real(0.5)*l_dt* q(i,j,k-1,n)*divu(i,j,k-1) : Real(0.0);

                stl += (fq)           ? Real(0.5)*l_dt*fq(i,j,k-1,n) : Real(0.0);
            }

#ifdef AMREX_USE_MOVING_EB
            const int no_eb_flow_zhi = !(values_on_eb_inflow) ? 1 : 0;
#else
            const int no_eb_flow_zhi = !(values_on_eb_inflow) ? 1 :
                ((Math::abs(values_on_eb_inflow(i,j,k+1,n)) > 0. ||
                  Math::abs(values_on_eb_inflow(i,j,k  ,n)) > 0.) ? 0 : 1);
#endif

            // If we can't compute good transverse terms, don't use any d/dt terms at all
            if (apx(i+1,j,k) > 0. && apx(i,j,k) > 0. && apy(i,j+1,k) > 0. && apy(i,j,k) > 0. && no_eb_flow_zhi)
            {
                // Here we add  dt/2 (-q w_z - (u q)_x - (v q)_y) to the term that is already
                //     q + dz/2 q_z + dt/2 (-w q_z) to get
                //     q + dz/2 q_z - dt/2 (w q_z  + q w_z + (u q)_x + (v q)_y) which is equivalent to
                // --> q + dz/2 q_z - dt/2 ( div (uvec q) )
                Real qwzh = (apz(i,j,k+1)*w_mac(i,j,k+1) - apz(i,j,k)*w_mac(i,j,k)) * q(i,j,k,n);

                sth += ( - (Real(0.5)*dtdz) * qwzh
                         - (Real(0.5)*dtdx)*(apx(i+1,j  ,k)*xylo(i+1,j  ,k,n)*u_mac(i+1,j  ,k)
                                      -apx(i  ,j  ,k)*xylo(i  ,j  ,k,n)*u_mac(i  ,j  ,k))
                         - (Real(0.5)*dtdy)*(apy(i  ,j+1,k)*yxlo(i  ,j+1,k,n)*v_mac(i  ,j+1,k)
                                      -apy(i  ,j  ,k)*yxlo(i  ,j  ,k,n)*v_mac(i  ,j  ,k)) ) / vfrac_arr(i,j,k);

                // Here we adjust for non-conservative by removing the q divu contribution to get
                //     q + dz/2 q_z - dt/2 ( div (uvec q) - q divu ) which is equivalent to
                // --> q + dz/2 q_z - dt/2 ( uvec dot grad q)
                sth += (!iconserv[n]) ? Real(0.5)*l_dt* q(i,j,k,n)*divu(i,j,k) : Real(0.0);

                sth += (fq)           ? Real(0.5)*l_dt*fq(i,j,k,n) : Real(0.0);
            }

            const auto bc = HydroBC::getBC(i, j, k, n, domain, pbc, bc_arr);
            HydroBC::SetZEdgeBCs(i, j, k, n, q, stl, sth, bc.lo(2), dlo.z, bc.hi(2), dhi.z, is_velocity);

            if ( (k==dlo.z) && (bc.lo(2) == BCType::foextrap || bc.lo(2) == BCType::hoextrap) )
            {
                if ( w_mac(i,j,k) >= 0. && n==ZVEL && is_velocity ) sth = amrex::min(sth,0.0_rt);
                stl = sth;
            }
            if ( (k==dhi.z+1) && (bc.hi(2) == BCType::foextrap || bc.hi(2) == BCType::hoextrap) )
            {
                if ( w_mac(i,j,k) <= 0. && n==ZVEL && is_velocity ) stl = amrex::max(stl,0.0_rt);
                sth = stl;
            }
            Real temp = (w_mac(i,j,k) >= 0.) ? stl : sth;
            temp = (amrex::Math::abs(w_mac(i,j,k)) < small_vel) ? Real(0.5)*(stl + sth) : temp;
            zedge(i,j,k,n) = temp;

        } else {
            zedge(i,j,k,n) = Real(0.0);
        }
    });

}
/** @} */

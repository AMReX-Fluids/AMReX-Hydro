/**
 * \file hydro_ebgodunov_edge_state_2D.cpp
 * \addtogroup EBGodunov
 * @{
 *
 */

#include <hydro_bcs_K.H>
#include <hydro_ebgodunov.H>
#include <hydro_ebgodunov_plm.H>
#include <AMReX_MultiCutFab.H>
#include <AMReX_EBMultiFabUtil_2D_C.H>

using namespace amrex;

void
EBGodunov::ComputeEdgeState ( Box const& bx, int ncomp,
                              Array4<Real const> const& q,
                              Array4<Real> const& xedge,
                              Array4<Real> const& yedge,
                              Array4<Real const> const& u_mac,
                              Array4<Real const> const& v_mac,
                              Array4<Real const> const& divu,
                              Array4<Real const> const& fq,
                              Geometry const& geom,
                              Real l_dt,
                              Vector<amrex::BCRec> const& h_bcrec,
                              BCRec const*  pbc,
                              int const* iconserv,
                              Real* p,
                              Array4<EBCellFlag const> const& flag_arr,
                              Array4<Real const> const& apx,
                              Array4<Real const> const& apy,
                              Array4<Real const> const& vfrac_arr,
                              Array4<Real const> const& fcx,
                              Array4<Real const> const& fcy,
                              Array4<Real const> const& ccent_arr,
                              bool is_velocity,
                              Array4<Real const> const& values_on_eb_inflow)
{
    Box const& xbx = amrex::surroundingNodes(bx,0);
    Box const& ybx = amrex::surroundingNodes(bx,1);
    Box const& bxg1 = amrex::grow(bx,1);

    // Start with above and grow 1 tangentially
    Box xebx = Box(xbx).grow(1,1);
    Box yebx = Box(ybx).grow(0,1);

    const Real dx = geom.CellSize(0);
    const Real dy = geom.CellSize(1);
    Real dtdx = l_dt/dx;
    Real dtdy = l_dt/dy;

    Box const& domain = geom.Domain();
    const auto dlo = amrex::lbound(domain);
    const auto dhi = amrex::ubound(domain);

    Array4<Real> Imx = makeArray4(p, bxg1, ncomp);
    p +=         Imx.size();
    Array4<Real> Ipx = makeArray4(p, bxg1, ncomp);
    p +=         Ipx.size();
    Array4<Real> Imy = makeArray4(p, bxg1, ncomp);
    p +=         Imy.size();
    Array4<Real> Ipy = makeArray4(p, bxg1, ncomp);
    p +=         Ipy.size();

    Array4<Real> xlo = makeArray4(p, xebx, ncomp);
    p +=         xlo.size();
    Array4<Real> xhi = makeArray4(p, xebx, ncomp);
    p +=         xhi.size();
    Array4<Real> ylo = makeArray4(p, yebx, ncomp);
    p +=         ylo.size();
    Array4<Real> yhi = makeArray4(p, yebx, ncomp);
    p +=         yhi.size();

    Array4<Real> xyzlo = makeArray4(p, bxg1, ncomp);
    p +=         xyzlo.size();
    Array4<Real> xyzhi = makeArray4(p, bxg1, ncomp);
    p +=         xyzhi.size();


    // Initialize this way out of an abundance of paranoia
    amrex::ParallelFor(
        Box(Imx), ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
           Imx(i,j,k,n) = i*1.e10 + j*1.e20 + k*1.30 + n*1.e4;
        });
    amrex::ParallelFor(
        Box(Imy), ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
           Imy(i,j,k,n) = i*1.e10 + j*1.e20 + k*1.30 + n*1.e4;
        });


    EBPLM::PredictStateOnXFace( xebx, ncomp, Imx, Ipx, q, u_mac,
                                flag_arr, vfrac_arr,
                                AMREX_D_DECL(fcx,fcy,fcz),ccent_arr,
                                geom, l_dt, h_bcrec, pbc, is_velocity);

    EBPLM::PredictStateOnYFace( yebx, ncomp, Imy, Ipy, q, v_mac,
                                flag_arr, vfrac_arr,
                                AMREX_D_DECL(fcx,fcy,fcz),ccent_arr,
                                geom, l_dt, h_bcrec, pbc, is_velocity);

    // No need for the intermediate upwinded edge state like in 3D. Could combine
    // with x/yzlo loop below ...
    amrex::ParallelFor(
        xebx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real lo = Ipx(i-1,j,k,n);
            Real hi = Imx(i  ,j,k,n);

            auto bc = pbc[n];

            HydroBC::SetXEdgeBCs(i, j, k, n, q, lo, hi, bc.lo(0), dlo.x, bc.hi(0), dhi.x, is_velocity);

            xlo(i,j,k,n) = lo;
            xhi(i,j,k,n) = hi;
        },
        yebx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real lo = Ipy(i,j-1,k,n);
            Real hi = Imy(i,j  ,k,n);

            auto bc = pbc[n];

            HydroBC::SetYEdgeBCs(i, j, k, n, q, lo, hi, bc.lo(1), dlo.y, bc.hi(1), dhi.y, is_velocity);

            ylo(i,j,k,n) = lo;
            yhi(i,j,k,n) = hi;
        });

    //
    // Upwinding on y-faces to use as transverse terms for x-faces
    //
    Array4<Real> yzlo = makeArray4(xyzlo.dataPtr(), yebx, ncomp);
    amrex::ParallelFor(
    Box(yzlo), ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (apy(i,j,k) > 0.)
        {
            const auto bc = pbc[n];
            Real l_yzlo, l_yzhi;

            l_yzlo = ylo(i,j,k,n);
            l_yzhi = yhi(i,j,k,n);
            Real vad = v_mac(i,j,k);
            HydroBC::SetYEdgeBCs(i, j, k, n, q, l_yzlo, l_yzhi, bc.lo(1), dlo.y, bc.hi(1), dhi.y, is_velocity);

            Real st = (vad >= 0.) ? l_yzlo : l_yzhi;
            Real fu = (amrex::Math::abs(vad) < small_vel) ? 0.0 : 1.0;
            yzlo(i,j,k,n) = fu*st + (1.0 - fu) * 0.5 * (l_yzhi + l_yzlo);
        } else {
            yzlo(i,j,k,n) = 0.;
        }
    });
    //
    amrex::ParallelFor(xbx, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        Real stl, sth;

        if (apx(i,j,k) > 0.)
        {
            stl = xlo(i,j,k,n);
            sth = xhi(i,j,k,n);

#ifdef AMREX_USE_MOVING_EB
            const int no_eb_flow_xlo = !(values_on_eb_inflow) ? 1 : 0;
#else
            const int no_eb_flow_xlo = !(values_on_eb_inflow) ? 1 :
                ((Math::abs(values_on_eb_inflow(i  ,j,k,n)) > 0. ||
                  Math::abs(values_on_eb_inflow(i-1,j,k,n)) > 0.) ? 0 : 1);
#endif

            // If we can't compute good transverse terms, don't use any d/dt terms at all
            if (apy(i-1,j,k) > 0. && apy(i-1,j+1,k) > 0. && no_eb_flow_xlo)
            {
                // Here we add  dt/2 (-q u_x - (v q)_y) to the term that is already
                //     q + dx/2 q_x + dt/2 (-u q_x) to get
                //     q + dx/2 q_x - dt/2 (u q_x  + q u_x + (v q)_y) which is equivalent to
                // --> q + dx/2 q_x - dt/2 ( div (uvec q) )
                Real quxl = (apx(i,j,k)*u_mac(i,j,k) - apx(i-1,j,k)*u_mac(i-1,j,k)) * q(i-1,j,k,n);
                stl += ( - (0.5*dtdx) * quxl
                         - (0.5*dtdy) * (apy(i-1,j+1,k)*yzlo(i-1,j+1,k  ,n)*v_mac(i-1,j+1,k  )
                                        -apy(i-1,j  ,k)*yzlo(i-1,j  ,k  ,n)*v_mac(i-1,j  ,k  )) ) / vfrac_arr(i-1,j,k);

                // Here we adjust for non-conservative by removing the q divu contribution to get
                //     q + dx/2 q_x - dt/2 ( div (uvec q) - q divu ) which is equivalent to
                // --> q + dx/2 q_x - dt/2 ( uvec dot grad q)
                stl += (!iconserv[n]) ?  0.5*l_dt* q(i-1,j,k,n)*divu(i-1,j,k) : 0.;

                stl += (fq)           ?  0.5*l_dt*fq(i-1,j,k,n) : 0.;
            }

#ifdef AMREX_USE_MOVING_EB
            const int no_eb_flow_xhi = !(values_on_eb_inflow) ? 1 : 0;
#else
            const int no_eb_flow_xhi = !(values_on_eb_inflow) ? 1 :
                ((Math::abs(values_on_eb_inflow(i+1,j,k,n)) > 0. ||
                  Math::abs(values_on_eb_inflow(i  ,j,k,n)) > 0.) ? 0 : 1);
#endif

            // If we can't compute good transverse terms, don't use any d/dt terms at all
            if (apy(i,j,k) > 0. && apy(i,j+1,k) > 0. && no_eb_flow_xhi)
            {
                // Here we add  dt/2 (-q u_x - (v q)_y) to the term that is already
                //     q + dx/2 q_x + dt/2 (-u q_x) to get
                //     q + dx/2 q_x - dt/2 (u q_x  + q u_x + (v q)_y)  which is equivalent to
                // --> q + dx/2 q_x - dt/2 ( div (uvec q) )
                Real quxh = (apx(i+1,j,k)*u_mac(i+1,j,k) - apx(i,j,k)*u_mac(i,j,k)) * q(i,j,k,n);
                sth += ( - (0.5*dtdx) * quxh
                         - (0.5*dtdy)*(apy(i,j+1,k)*yzlo(i,j+1,k,n)*v_mac(i,j+1,k)
                                      -apy(i,j  ,k)*yzlo(i,j  ,k,n)*v_mac(i,j  ,k)) ) / vfrac_arr(i,j,k);

                // Here we adjust for non-conservative by removing the q divu contribution to get
                //     q + dx/2 q_x - dt/2 ( div (uvec q) - q divu ) which is equivalent to
                // --> q + dx/2 q_x - dt/2 ( uvec dot grad q)
                sth += (!iconserv[n]) ? 0.5*l_dt* q(i  ,j,k,n)*divu(i,j,k) : 0.;

                sth += (fq)           ? 0.5*l_dt*fq(i  ,j,k,n) : 0.;
            }

            auto bc = pbc[n];
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
            temp = (amrex::Math::abs(u_mac(i,j,k)) < small_vel) ? 0.5*(stl + sth) : temp;
            xedge(i,j,k,n) = temp;

        } else {
           xedge(i,j,k,n) = 0.;
        }
    });

    //
    // Upwinding on x-faces to use as transverse terms for y-faces
    //
    Array4<Real> xzlo = makeArray4(xyzlo.dataPtr(), xebx, ncomp);
    amrex::ParallelFor(
    Box(xzlo), ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (apx(i,j,k) > 0.)
        {
            const auto bc = pbc[n];
            Real l_xzlo, l_xzhi;

            l_xzlo = xlo(i,j,k,n);
            l_xzhi = xhi(i,j,k,n);

            Real uad = u_mac(i,j,k);
            HydroBC::SetXEdgeBCs(i, j, k, n, q, l_xzlo, l_xzhi, bc.lo(0), dlo.x, bc.hi(0), dhi.x, is_velocity);

            Real st = (uad >= 0.) ? l_xzlo : l_xzhi;
            Real fu = (amrex::Math::abs(uad) < small_vel) ? 0.0 : 1.0;
            xzlo(i,j,k,n) = fu*st + (1.0 - fu) * 0.5 * (l_xzhi + l_xzlo);
        } else {
            xzlo(i,j,k,n) = 0.;
        }
    });
    //
    amrex::ParallelFor(ybx, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        Real stl, sth;

        if (apy(i,j,k) > 0.)
        {
            stl = ylo(i,j,k,n);
            sth = yhi(i,j,k,n);

#ifdef AMREX_USE_MOVING_EB
            const int no_eb_flow_ylo = !(values_on_eb_inflow) ? 1 : 0;
#else
            const int no_eb_flow_ylo = !(values_on_eb_inflow) ? 1 :
                ((Math::abs(values_on_eb_inflow(i,j  ,k,n)) > 0. ||
                  Math::abs(values_on_eb_inflow(i,j-1,k,n)) > 0.) ? 0 : 1);
#endif

            // If we can't compute good transverse terms, don't use any d/dt terms at all
            if (apx(i,j-1,k) > 0. && apx(i+1,j-1,k) > 0. && no_eb_flow_ylo)
            {
                // Here we add  dt/2 (-q v_y - (u q)_x) to the term that is already
                //     q + dy/2 q_y + dt/2 (-v q_y) to get
                //     q + dy/2 q_y - dt/2 (v q_y  + q v_y + (u q)_x) which is equivalent to
                // --> q + dy/2 q_y - dt/2 ( div (uvec q) )
                Real qvyl = (apy(i,j,k)*v_mac(i,j,k) - apy(i,j-1,k)*v_mac(i,j-1,k)) * q(i,j-1,k,n);
                stl += ( - (0.5*dtdy)*qvyl
                         - (0.5*dtdx)*(apx(i+1,j-1,k)*xzlo(i+1,j-1,k  ,n)*u_mac(i+1,j-1,k  )
                                      -apx(i  ,j-1,k)*xzlo(i  ,j-1,k  ,n)*u_mac(i  ,j-1,k  )) ) / vfrac_arr(i,j-1,k);

                // Here we adjust for non-conservative by removing the q divu contribution to get
                //     q + dy/2 q_y - dt/2 ( div (uvec q) - q divu ) which is equivalent to
                // --> q + dy/2 q_y - dt/2 ( uvec dot grad q)
                stl += (!iconserv[n]) ? 0.5*l_dt* q(i,j-1,k,n)*divu(i,j-1,k) : 0.;

                stl += (fq)           ? 0.5*l_dt*fq(i,j-1,k,n) : 0.;
            }

#ifdef AMREX_USE_MOVING_EB
            const int no_eb_flow_yhi = !(values_on_eb_inflow) ? 1 : 0;
#else
            const int no_eb_flow_yhi = !(values_on_eb_inflow) ? 1 :
                ((Math::abs(values_on_eb_inflow(i,j+1,k,n)) > 0. ||
                  Math::abs(values_on_eb_inflow(i,j  ,k,n)) > 0.) ? 0 : 1);
#endif

            // If we can't compute good transverse terms, don't use any d/dt terms at all
             if (apx(i,j,k) > 0. && apx(i+1,j,k) > 0. && no_eb_flow_yhi)
            {
                // Here we add  dt/2 (-q v_y - (u q)_x) to the term that is already
                //     q + dy/2 q_y + dt/2 (-v q_y) to get
                //     q + dy/2 q_y - dt/2 (v q_y  + q v_y + (u q)_x) which is equivalent to
                // --> q + dy/2 q_y - dt/2 ( div (uvec q) )
                Real qvyh = (apy(i,j+1,k)*v_mac(i,j+1,k) - apy(i,j,k)*v_mac(i,j,k)) * q(i,j,k,n);
                sth += ( - (0.5*dtdy)*qvyh
                         - (0.5*dtdx)*(apx(i+1,j,k)*xzlo(i+1,j,k  ,n)*u_mac(i+1,j,k  )
                                      -apx(i  ,j,k)*xzlo(i  ,j,k  ,n)*u_mac(i  ,j,k  )) ) / vfrac_arr(i,j  ,k);

                // Here we adjust for non-conservative by removing the q divu contribution to get
                //     q + dy/2 q_y - dt/2 ( div (uvec q) - q divu ) which is equivalent to
                // --> q + dy/2 q_y - dt/2 ( uvec dot grad q)
                sth += (!iconserv[n]) ? 0.5*l_dt* q(i,j,k,n)*divu(i,j,k) : 0.;

                sth += (fq)           ? 0.5*l_dt*fq(i,j,k,n) : 0.;
            }

            auto bc = pbc[n];
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
            temp = (amrex::Math::abs(v_mac(i,j,k)) < small_vel) ? 0.5*(stl + sth) : temp;
            yedge(i,j,k,n) = temp;

        } else {
            yedge(i,j,k,n) = 0.;
        }
    });

}
/** @} */

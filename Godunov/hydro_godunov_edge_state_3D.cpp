/**
 * \file hydro_godunov_edge_state_3D.cpp
 *
 * \addtogroup Godunov
 *  @{
 */

#include <hydro_godunov_plm.H>
#include <hydro_godunov_ppm.H>
#include <hydro_godunov.H>
#include <hydro_godunov_corner_couple.H>
#include <hydro_bcs_K.H>

using namespace amrex;

void
Godunov::ComputeEdgeState (Box const& bx, int ncomp,
                           Array4<Real const> const& q,
                           Array4<Real> const& xedge,
                           Array4<Real> const& yedge,
                           Array4<Real> const& zedge,
                           Array4<Real const> const& umac,
                           Array4<Real const> const& vmac,
                           Array4<Real const> const& wmac,
                           Array4<Real const> const& divu,
                           Array4<Real const> const& fq,
                           Geometry geom,
                           Real l_dt,
                           BCRec const* pbc, int const* iconserv,
                           const bool use_ppm,
                           const bool use_forces_in_trans,
                           const bool is_velocity,
                           const int limiter_type)
{
    Box const& xbx = amrex::surroundingNodes(bx,0);
    Box const& ybx = amrex::surroundingNodes(bx,1);
    Box const& zbx = amrex::surroundingNodes(bx,2);

    Box const& bxg1 = amrex::grow(bx,1);

    FArrayBox tmpfab(amrex::grow(bx,1),  (4*AMREX_SPACEDIM + 2)*ncomp,The_Async_Arena());
    Real* p   = tmpfab.dataPtr();

    Box xebox = Box(xbx).grow(1,1).grow(2,1);
    Box yebox = Box(ybx).grow(0,1).grow(2,1);
    Box zebox = Box(zbx).grow(0,1).grow(1,1);

    const Real dx = geom.CellSize(0);
    const Real dy = geom.CellSize(1);
    const Real dz = geom.CellSize(2);

    Real dtdx = l_dt/dx;
    Real dtdy = l_dt/dy;
    Real dtdz = l_dt/dz;

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
    Array4<Real> Imz = makeArray4(p, bxg1, ncomp);
    p +=         Imz.size();
    Array4<Real> Ipz = makeArray4(p, bxg1, ncomp);
    p +=         Ipz.size();
    Array4<Real> xlo = makeArray4(p, xebox, ncomp);
    p +=         xlo.size();
    Array4<Real> xhi = makeArray4(p, xebox, ncomp);
    p +=         xhi.size();
    Array4<Real> ylo = makeArray4(p, yebox, ncomp);
    p +=         ylo.size();
    Array4<Real> yhi = makeArray4(p, yebox, ncomp);
    p +=         yhi.size();
    Array4<Real> zlo = makeArray4(p, zebox, ncomp);
    p +=         zlo.size();
    Array4<Real> zhi = makeArray4(p, zebox, ncomp);
    p +=         zhi.size();
    Array4<Real> xyzlo = makeArray4(p, bxg1, ncomp);
    p +=         xyzlo.size();
    Array4<Real> xyzhi = makeArray4(p, bxg1, ncomp);

    // Use PPM to generate Im and Ip */
    if (use_ppm)
    {
        if(limiter_type == PPM::VanLeer) {
            auto limiter = PPM::vanleer();
            PPM::PredictStateOnFaces(bxg1,AMREX_D_DECL(Imx,Imy,Imz),
                                     AMREX_D_DECL(Ipx,Ipy,Ipz),
                                     AMREX_D_DECL(umac,vmac,wmac),
                                     q,geom,l_dt,pbc,ncomp,limiter);
        } else if ( limiter_type == PPM::WENOZ) {
            auto limiter = PPM::wenoz();
            PPM::PredictStateOnFaces(bxg1,AMREX_D_DECL(Imx,Imy,Imz),
                                     AMREX_D_DECL(Ipx,Ipy,Ipz),
                                     AMREX_D_DECL(umac,vmac,wmac),
                                     q,geom,l_dt,pbc,ncomp,limiter);
        } else {
            auto limiter = PPM::nolimiter();
            PPM::PredictStateOnFaces(bxg1,AMREX_D_DECL(Imx,Imy,Imz),
                                     AMREX_D_DECL(Ipx,Ipy,Ipz),
                                     AMREX_D_DECL(umac,vmac,wmac),
                                     q,geom,l_dt,pbc,ncomp,limiter);
        }
    // Use PLM to generate Im and Ip */
    }
    else
    {

        amrex::ParallelFor(xebox, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            PLM::PredictStateOnXFace(i, j, k, n, l_dt, dx, Imx(i,j,k,n), Ipx(i-1,j,k,n),
                                     q, umac(i,j,k), pbc[n], dlo.x, dhi.x, is_velocity);
        });

        amrex::ParallelFor(yebox, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            PLM::PredictStateOnYFace(i, j, k, n, l_dt, dy, Imy(i,j,k,n), Ipy(i,j-1,k,n),
                                     q, vmac(i,j,k), pbc[n], dlo.y, dhi.y, is_velocity);
        });
        amrex::ParallelFor(zebox, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            PLM::PredictStateOnZFace(i, j, k, n, l_dt, dz, Imz(i,j,k,n), Ipz(i,j,k-1,n),
                                     q, wmac(i,j,k), pbc[n], dlo.z, dhi.z, is_velocity);
        });
    }


    amrex::ParallelFor(
    xebox, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        Real uad = umac(i,j,k);
        Real fux = (amrex::Math::abs(uad) < small_vel)? 0. : 1.;
        bool uval = uad >= 0.;
        Real lo = Ipx(i-1,j,k,n);
        Real hi = Imx(i  ,j,k,n);

        if (use_forces_in_trans && fq)
        {
            lo += 0.5*l_dt*fq(i-1,j,k,n);
            hi += 0.5*l_dt*fq(i  ,j,k,n);
        }

        auto bc = pbc[n];

        HydroBC::SetXEdgeBCs(i, j, k, n, q, lo, hi, bc.lo(0), dlo.x, bc.hi(0), dhi.x, is_velocity);
        xlo(i,j,k,n) = lo;
        xhi(i,j,k,n) = hi;
        Real st = (uval) ? lo : hi;
        Imx(i,j,k,n) = fux*st + (1. - fux)*0.5*(hi + lo);

    },
    yebox, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        Real vad = vmac(i,j,k);
        Real fuy = (amrex::Math::abs(vad) < small_vel)? 0. : 1.;
        bool vval = vad >= 0.;
        Real lo = Ipy(i,j-1,k,n);
        Real hi = Imy(i,j  ,k,n);

        if (use_forces_in_trans && fq)
        {
            lo += 0.5*l_dt*fq(i,j-1,k,n);
            hi += 0.5*l_dt*fq(i,j  ,k,n);
        }

        auto bc = pbc[n];

        HydroBC::SetYEdgeBCs(i, j, k, n, q, lo, hi, bc.lo(1), dlo.y, bc.hi(1), dhi.y, is_velocity);

        ylo(i,j,k,n) = lo;
        yhi(i,j,k,n) = hi;
        Real st = (vval) ? lo : hi;
        Imy(i,j,k,n) = fuy*st + (1. - fuy)*0.5*(hi + lo);
    },
    zebox, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        Real wad = wmac(i,j,k);
        Real fuz = (amrex::Math::abs(wad) < small_vel) ? 0. : 1.;
        bool wval = wad >= 0.;
        Real lo = Ipz(i,j,k-1,n);
        Real hi = Imz(i,j,k  ,n);

        if (use_forces_in_trans && fq)
        {
            lo += 0.5*l_dt*fq(i,j,k-1,n);
            hi += 0.5*l_dt*fq(i,j,k  ,n);
        }

        auto bc = pbc[n];

        HydroBC::SetZEdgeBCs(i, j, k, n, q, lo, hi, bc.lo(2), dlo.z, bc.hi(2), dhi.z, is_velocity);

        zlo(i,j,k,n) = lo;
        zhi(i,j,k,n) = hi;
        Real st = (wval) ? lo : hi;
        Imz(i,j,k,n) = fuz*st + (1. - fuz)*0.5*(hi + lo);
    }
    );



    //
    // x-direction
    //
    Box const& xbxtmp = amrex::grow(bx,0,1);
    Array4<Real> yzlo = makeArray4(xyzlo.dataPtr(), amrex::surroundingNodes(xbxtmp,1), ncomp);
    Array4<Real> zylo = makeArray4(xyzhi.dataPtr(), amrex::surroundingNodes(xbxtmp,2), ncomp);
    amrex::ParallelFor(
    Box(zylo), ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        const auto bc = pbc[n];
        Real l_zylo, l_zyhi;
        GodunovCornerCouple::AddCornerCoupleTermZY(l_zylo, l_zyhi,
                              i, j, k, n, l_dt, dy, iconserv[n],
                              zlo(i,j,k,n), zhi(i,j,k,n),
                              q, divu, vmac, Imy);

        Real wad = wmac(i,j,k);
        HydroBC::SetZEdgeBCs(i, j, k, n, q, l_zylo, l_zyhi, bc.lo(2), dlo.z, bc.hi(2), dhi.z, is_velocity);

        Real st = (wad >= 0.) ? l_zylo : l_zyhi;
        Real fu = (amrex::Math::abs(wad) < small_vel) ? 0.0 : 1.0;
        zylo(i,j,k,n) = fu*st + (1.0 - fu) * 0.5 * (l_zyhi + l_zylo);
    },
    Box(yzlo), ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        const auto bc = pbc[n];
        Real l_yzlo, l_yzhi;
        GodunovCornerCouple::AddCornerCoupleTermYZ(l_yzlo, l_yzhi,
                              i, j, k, n, l_dt, dz, iconserv[n],
                              ylo(i,j,k,n), yhi(i,j,k,n),
                              q, divu, wmac, Imz);

        Real vad = vmac(i,j,k);
        HydroBC::SetYEdgeBCs(i, j, k, n, q, l_yzlo, l_yzhi, bc.lo(1), dlo.y, bc.hi(1), dhi.y, is_velocity);

        Real st = (vad >= 0.) ? l_yzlo : l_yzhi;
        Real fu = (amrex::Math::abs(vad) < small_vel) ? 0.0 : 1.0;
        yzlo(i,j,k,n) = fu*st + (1.0 - fu) * 0.5 * (l_yzhi + l_yzlo);
    });


    //
    amrex::ParallelFor(xbx, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
    Real stl = xlo(i,j,k,n);
    Real sth = xhi(i,j,k,n);

        // To match EBGodunov
        // Here we add  dt/2 (-q u_x - (v q)_y - (w q)_z) to the term that is already
    //     q + dx/2 q_x + dt/2 (-u q_x) to get
    //     q + dx/2 q_x - dt/2 (u q_x  + q u_x + (v q)_y + (w q)_z) which is equivalent to
    // --> q + dx/2 q_x - dt/2 ( div (uvec q) )
    Real quxl = (umac(i,j,k) - umac(i-1,j,k)) * q(i-1,j,k,n);
    stl += ( - (0.5*dtdx) * quxl
         - (0.5*dtdy)*(yzlo(i-1,j+1,k  ,n)*vmac(i-1,j+1,k  )
                     - yzlo(i-1,j  ,k  ,n)*vmac(i-1,j  ,k  ))
         - (0.5*dtdz)*(zylo(i-1,j  ,k+1,n)*wmac(i-1,j  ,k+1)
                   - zylo(i-1,j  ,k  ,n)*wmac(i-1,j  ,k  )) );

    // Here we adjust for non-conservative by removing the q divu contribution to get
    //     q + dx/2 q_x - dt/2 ( div (uvec q) - q divu ) which is equivalent to
    // --> q + dx/2 q_x - dt/2 ( uvec dot grad q)
    stl += (!iconserv[n])               ? 0.5*l_dt* q(i-1,j,k,n)*divu(i-1,j,k) : 0.;

    stl += (!use_forces_in_trans && fq) ? 0.5*l_dt*fq(i-1,j,k,n) : 0.;

    // High side
    Real quxh = (umac(i+1,j,k) - umac(i,j,k)) * q(i,j,k,n);
    sth += ( - (0.5*dtdx) * quxh
         - (0.5*dtdy)*(yzlo(i,j+1,k  ,n)*vmac(i,j+1,k  )
                     - yzlo(i,j  ,k  ,n)*vmac(i,j  ,k  ))
         - (0.5*dtdz)*(zylo(i,j  ,k+1,n)*wmac(i,j  ,k+1)
                         - zylo(i,j  ,k  ,n)*wmac(i,j  ,k  )) );

    sth += (!iconserv[n])               ? 0.5*l_dt* q(i  ,j,k,n)*divu(i,j,k) : 0.;

    sth += (!use_forces_in_trans && fq) ? 0.5*l_dt*fq(i  ,j,k,n) : 0.;


    auto bc = pbc[n];
    HydroBC::SetXEdgeBCs(i, j, k, n, q, stl, sth, bc.lo(0), dlo.x, bc.hi(0), dhi.x, is_velocity);

        if ( (i==dlo.x) && (bc.lo(0) == BCType::foextrap || bc.lo(0) == BCType::hoextrap) )
        {
            if ( umac(i,j,k) >= 0. && n==XVEL && is_velocity )  sth = amrex::min(sth,0.0_rt);
            stl = sth;
        }
        if ( (i==dhi.x+1) && (bc.hi(0) == BCType::foextrap || bc.hi(0) == BCType::hoextrap) )
        {
            if ( umac(i,j,k) <= 0. && n==XVEL && is_velocity ) stl = amrex::max(stl,0.0_rt);
             sth = stl;
        }

        Real temp = (umac(i,j,k) >= 0.) ? stl : sth;
        temp = (amrex::Math::abs(umac(i,j,k)) < small_vel) ? 0.5*(stl + sth) : temp;
        xedge(i,j,k,n) = temp;
    });

    //
    // y-direction
    //
    Box const& ybxtmp = amrex::grow(bx,1,1);
    Array4<Real> xzlo = makeArray4(xyzlo.dataPtr(), amrex::surroundingNodes(ybxtmp,0), ncomp);
    Array4<Real> zxlo = makeArray4(xyzhi.dataPtr(), amrex::surroundingNodes(ybxtmp,2), ncomp);
    amrex::ParallelFor(
    Box(xzlo), ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        const auto bc = pbc[n];
        Real l_xzlo, l_xzhi;
        GodunovCornerCouple::AddCornerCoupleTermXZ(l_xzlo, l_xzhi,
                              i, j, k, n, l_dt, dz, iconserv[n],
                              xlo(i,j,k,n),  xhi(i,j,k,n),
                              q, divu, wmac, Imz);

        Real uad = umac(i,j,k);
        HydroBC::SetXEdgeBCs(i, j, k, n, q, l_xzlo, l_xzhi, bc.lo(0), dlo.x, bc.hi(0), dhi.x, is_velocity);

        Real st = (uad >= 0.) ? l_xzlo : l_xzhi;
        Real fu = (amrex::Math::abs(uad) < small_vel) ? 0.0 : 1.0;
        xzlo(i,j,k,n) = fu*st + (1.0 - fu) * 0.5 * (l_xzhi + l_xzlo);
    },
    Box(zxlo), ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        const auto bc = pbc[n];
        Real l_zxlo, l_zxhi;
        GodunovCornerCouple::AddCornerCoupleTermZX(l_zxlo, l_zxhi,
                              i, j, k, n, l_dt, dx, iconserv[n],
                              zlo(i,j,k,n), zhi(i,j,k,n),
                              q, divu, umac, Imx);

        Real wad = wmac(i,j,k);
        HydroBC::SetZEdgeBCs(i, j, k, n, q, l_zxlo, l_zxhi, bc.lo(2), dlo.z, bc.hi(2), dhi.z, is_velocity);

        Real st = (wad >= 0.) ? l_zxlo : l_zxhi;
        Real fu = (amrex::Math::abs(wad) < small_vel) ? 0.0 : 1.0;
        zxlo(i,j,k,n) = fu*st + (1.0 - fu) * 0.5 * (l_zxhi + l_zxlo);
    });

    //
    amrex::ParallelFor(ybx, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
    Real stl = ylo(i,j,k,n);
    Real sth = yhi(i,j,k,n);

    // Here we add  dt/2 (-q v_y - (u q)_x - (w q)_z) to the term that is already
    //     q + dy/2 q_y + dt/2 (-v q_y) to get
    //     q + dy/2 q_y - dt/2 (v q_y  + q v_y + (u q)_x + (w q)_z) which is equivalent to
    // --> q + dy/2 q_y - dt/2 ( div (uvec q) )
    Real qvyl = (vmac(i,j,k) - vmac(i,j-1,k)) * q(i,j-1,k,n);
    stl += ( - (0.5*dtdy) * qvyl
         - (0.5*dtdx)*(xzlo(i+1,j-1,k  ,n)*umac(i+1,j-1,k  )
                 - xzlo(i  ,j-1,k  ,n)*umac(i  ,j-1,k  ))
         - (0.5*dtdz)*(zxlo(i  ,j-1,k+1,n)*wmac(i  ,j-1,k+1)
                 - zxlo(i  ,j-1,k  ,n)*wmac(i  ,j-1,k  )) );

    // Here we adjust for non-conservative by removing the q divu contribution to get
    //     q + dy/2 q_y - dt/2 ( div (uvec q) - q divu ) which is equivalent to
    // --> q + dy/2 q_y - dt/2 ( uvec dot grad q)
    stl += (!iconserv[n]) ? 0.5*l_dt* q(i,j-1,k,n)*divu(i,j-1,k) : 0.;

    stl += (!use_forces_in_trans && fq)           ? 0.5*l_dt*fq(i,j-1,k,n) : 0.;

    // High side
    Real qvyh = (vmac(i,j+1,k) - vmac(i,j,k)) * q(i,j,k,n);
    sth += ( - (0.5*dtdy) * qvyh
         - (0.5*dtdx)*(xzlo(i+1,j,k  ,n)*umac(i+1,j,k  )
                 - xzlo(i  ,j,k  ,n)*umac(i  ,j,k  ))
         - (0.5*dtdz)*(zxlo(i  ,j,k+1,n)*wmac(i  ,j,k+1)
                     - zxlo(i  ,j,k  ,n)*wmac(i  ,j,k  )) );

    sth += (!iconserv[n])               ? 0.5*l_dt* q(i,j,k,n)*divu(i,j,k) : 0.;

    sth += (!use_forces_in_trans && fq) ? 0.5*l_dt*fq(i,j,k,n) : 0.;


        auto bc = pbc[n];
        HydroBC::SetYEdgeBCs(i, j, k, n, q, stl, sth, bc.lo(1), dlo.y, bc.hi(1), dhi.y, is_velocity);

        if ( (j==dlo.y) && (bc.lo(1) == BCType::foextrap || bc.lo(1) == BCType::hoextrap) )
        {
            if ( vmac(i,j,k) >= 0. && n==YVEL && is_velocity ) sth = amrex::min(sth,0.0_rt);
            stl = sth;
        }
        if ( (j==dhi.y+1) && (bc.hi(1) == BCType::foextrap || bc.hi(1) == BCType::hoextrap) )
        {
            if ( vmac(i,j,k) <= 0. && n==YVEL && is_velocity ) stl = amrex::max(stl,0.0_rt);
            sth = stl;
        }

        Real temp = (vmac(i,j,k) >= 0.) ? stl : sth;
        temp = (amrex::Math::abs(vmac(i,j,k)) < small_vel) ? 0.5*(stl + sth) : temp;
        yedge(i,j,k,n) = temp;
    });

    //
    // z-direcion
    //
    Box const& zbxtmp = amrex::grow(bx,2,1);
    Array4<Real> xylo = makeArray4(xyzlo.dataPtr(), amrex::surroundingNodes(zbxtmp,0), ncomp);
    Array4<Real> yxlo = makeArray4(xyzhi.dataPtr(), amrex::surroundingNodes(zbxtmp,1), ncomp);
    amrex::ParallelFor(
    Box(xylo), ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        const auto bc = pbc[n];
        Real l_xylo, l_xyhi;
        GodunovCornerCouple::AddCornerCoupleTermXY(l_xylo, l_xyhi,
                              i, j, k, n, l_dt, dy, iconserv[n],
                              xlo(i,j,k,n), xhi(i,j,k,n),
                              q, divu, vmac, Imy);

        Real uad = umac(i,j,k);
        HydroBC::SetXEdgeBCs(i, j, k, n, q, l_xylo, l_xyhi, bc.lo(0), dlo.x, bc.hi(0), dhi.x, is_velocity);

        Real st = (uad >= 0.) ? l_xylo : l_xyhi;
        Real fu = (amrex::Math::abs(uad) < small_vel) ? 0.0 : 1.0;
        xylo(i,j,k,n) = fu*st + (1.0 - fu) * 0.5 * (l_xyhi + l_xylo);
    },
    Box(yxlo), ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        const auto bc = pbc[n];
        Real l_yxlo, l_yxhi;
        GodunovCornerCouple::AddCornerCoupleTermYX(l_yxlo, l_yxhi,
                              i, j, k, n, l_dt, dx, iconserv[n],
                              ylo(i,j,k,n), yhi(i,j,k,n),
                              q, divu, umac, Imx);

        Real vad = vmac(i,j,k);
        HydroBC::SetYEdgeBCs(i, j, k, n, q, l_yxlo, l_yxhi, bc.lo(1), dlo.y, bc.hi(1), dhi.y, is_velocity);

        Real st = (vad >= 0.) ? l_yxlo : l_yxhi;
        Real fu = (amrex::Math::abs(vad) < small_vel) ? 0.0 : 1.0;
        yxlo(i,j,k,n) = fu*st + (1.0 - fu) * 0.5 * (l_yxhi + l_yxlo);
    });
    //

    amrex::ParallelFor(zbx, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        Real stl = zlo(i,j,k,n);
    Real sth = zhi(i,j,k,n);

    // Here we add  dt/2 (-q w_z - (u q)_x - (v q)_y) to the term that is already
    //     q + dz/2 q_z + dt/2 (-w q_z) to get
    //     q + dz/2 q_z - dt/2 (w q_z  + q w_z + (u q)_x + (v q)_y) which is equivalent to
    // --> q + dz/2 q_z - dt/2 ( div (uvec q) )
    Real qwzl = (wmac(i,j,k) - wmac(i,j,k-1)) * q(i,j,k-1,n);
    stl += ( - (0.5*dtdz) * qwzl
         - (0.5*dtdx)*(xylo(i+1,j  ,k-1,n)*umac(i+1,j  ,k-1)
                  -xylo(i  ,j  ,k-1,n)*umac(i  ,j  ,k-1))
         - (0.5*dtdy)*(yxlo(i  ,j+1,k-1,n)*vmac(i  ,j+1,k-1)
                  -yxlo(i  ,j  ,k-1,n)*vmac(i  ,j  ,k-1)) );

    // Here we adjust for non-conservative by removing the q divu contribution to get
    //     q + dz/2 q_z - dt/2 ( div (uvec q) - q divu ) which is equivalent to
    // --> q + dz/2 q_z - dt/2 ( uvec dot grad q)
    stl += (!iconserv[n])               ? 0.5*l_dt* q(i,j,k-1,n)*divu(i,j,k-1) : 0.;

    stl += (!use_forces_in_trans && fq) ? 0.5*l_dt*fq(i,j,k-1,n) : 0.;

    // High side
    Real qwzh = (wmac(i,j,k+1) - wmac(i,j,k)) * q(i,j,k,n);
    sth += ( - (0.5*dtdz) * qwzh
         - (0.5*dtdx)*(xylo(i+1,j  ,k,n)*umac(i+1,j  ,k)
                  -xylo(i  ,j  ,k,n)*umac(i  ,j  ,k))
         - (0.5*dtdy)*(yxlo(i  ,j+1,k,n)*vmac(i  ,j+1,k)
                  -yxlo(i  ,j  ,k,n)*vmac(i  ,j  ,k)) );

    sth += (!iconserv[n])               ? 0.5*l_dt* q(i,j,k,n)*divu(i,j,k) : 0.;

    sth += (!use_forces_in_trans && fq) ? 0.5*l_dt*fq(i,j,k,n) : 0.;



        auto bc = pbc[n];
        HydroBC::SetZEdgeBCs(i, j, k, n, q, stl, sth, bc.lo(2), dlo.z, bc.hi(2), dhi.z, is_velocity);

        if ( (k==dlo.z) && (bc.lo(2) == BCType::foextrap || bc.lo(2) == BCType::hoextrap) )
        {
            if ( wmac(i,j,k) >= 0. && n==ZVEL && is_velocity ) sth = amrex::min(sth,0.0_rt);
            stl = sth;
        }
        if ( (k==dhi.z+1) && (bc.hi(2) == BCType::foextrap || bc.hi(2) == BCType::hoextrap) )
        {
            if ( wmac(i,j,k) <= 0. && n==ZVEL && is_velocity ) stl = amrex::max(stl,0.0_rt);
            sth = stl;
        }

        Real temp = (wmac(i,j,k) >= 0.) ? stl : sth;
        temp = (amrex::Math::abs(wmac(i,j,k)) < small_vel) ? 0.5*(stl + sth) : temp;
        zedge(i,j,k,n) = temp;
    });

}
/** @} */

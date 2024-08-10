/**
 * \file hydro_godunov_extrap_vel_to_faces_2D.cpp
 *
 * \addtogroup Godunov
 *  @{
 */

#include <hydro_godunov_ppm.H>
#include <hydro_godunov_plm.H>
#include <hydro_godunov.H>
#include <hydro_bcs_K.H>

using namespace amrex;

void
Godunov::ExtrapVelToFaces ( MultiFab const& a_vel,
                            MultiFab const& a_forces,
                            MultiFab& a_umac,
                            MultiFab& a_vmac,
                            const Vector<BCRec> & h_bcrec,
                            const        BCRec  * d_bcrec,
                            const Geometry& geom, Real l_dt,
                            const bool use_ppm, const bool use_forces_in_trans,
                            const int limiter_type,
                            const bool allow_inflow_on_outflow,
                            iMultiFab const* BC_MF)
{
    Box const& domain = geom.Domain();
    const Real* dx    = geom.CellSize();
    const int ncomp = AMREX_SPACEDIM;
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        FArrayBox scratch;
        for (MFIter mfi(a_vel,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();
            Box const& bxg1 = amrex::grow(bx,1);

            Box const& xbx = mfi.nodaltilebox(0);
            Box const& ybx = mfi.nodaltilebox(1);

            Array4<Real> const& umac = a_umac.array(mfi);
            Array4<Real> const& vmac = a_vmac.array(mfi);

            Array4<Real const> const& vel = a_vel.const_array(mfi);
            Array4<Real const> const& f   = a_forces.const_array(mfi);
            Array4<int  const> const& bc_arr = BC_MF ? BC_MF->const_array(mfi)
                                                     : Array4<int const> {};

            scratch.resize(bxg1, (ncomp*4 + 1)*AMREX_SPACEDIM);
            Real* p = scratch.dataPtr();

            Array4<Real> Imx = makeArray4(p,bxg1,ncomp);
            p +=         Imx.size();
            Array4<Real> Ipx = makeArray4(p,bxg1,ncomp);
            p +=         Ipx.size();
            Array4<Real> Imy = makeArray4(p,bxg1,ncomp);
            p +=         Imy.size();
            Array4<Real> Ipy = makeArray4(p,bxg1,ncomp);
            p +=         Ipy.size();
            Array4<Real> u_ad = makeArray4(p,Box(bx).grow(1,1).surroundingNodes(0),1);
            p +=         u_ad.size();
            Array4<Real> v_ad = makeArray4(p,Box(bx).grow(0,1).surroundingNodes(1),1);
            p +=         v_ad.size();

            if (use_ppm)
            {
                if(limiter_type == PPM::VanLeer) {
                    auto limiter = PPM::vanleer();
                    PPM::PredictVelOnFaces( bxg1,
                                            Imx, Imy, Ipx, Ipy,
                                            vel, vel,
                                            geom, l_dt, d_bcrec, limiter);
                } else if (limiter_type == PPM::WENOZ) {
                    auto limiter = PPM::wenoz();
                    PPM::PredictVelOnFaces( bxg1,
                                            Imx, Imy, Ipx, Ipy,
                                            vel, vel,
                                            geom, l_dt, d_bcrec, limiter);
                } else if (limiter_type == PPM::WENO_JS) {
                    auto limiter = PPM::weno_js();
                    PPM::PredictVelOnFaces( bxg1,
                                            Imx, Imy, Ipx, Ipy,
                                            vel, vel,
                                            geom, l_dt, d_bcrec, limiter);
                } else {
                    auto limiter = PPM::nolimiter();
                    PPM::PredictVelOnFaces( bxg1,
                                            Imx, Imy, Ipx, Ipy,
                                            vel, vel,
                                            geom, l_dt, d_bcrec, limiter);
                }
            }
            else
            {
                PLM::PredictVelOnXFace( Box(u_ad), AMREX_SPACEDIM, Imx, Ipx, vel, vel,
                                        geom, l_dt, h_bcrec, d_bcrec, bc_arr);
                PLM::PredictVelOnYFace( Box(v_ad), AMREX_SPACEDIM, Imy, Ipy, vel, vel,
                                        geom, l_dt, h_bcrec, d_bcrec, bc_arr);
            }

            ComputeAdvectiveVel( Box(u_ad), Box(v_ad),
                                 u_ad, v_ad,
                                 Imx, Imy, Ipx, Ipy,
                                 vel, f, domain, l_dt, d_bcrec, use_forces_in_trans,
                                 bc_arr);

            ExtrapVelToFacesOnBox( bx, ncomp, xbx, ybx,
                                   umac, vmac, vel,
                                   u_ad, v_ad,
                                   Imx, Imy, Ipx, Ipy,
                                   f, domain, dx, l_dt, d_bcrec, use_forces_in_trans, p,
                                   allow_inflow_on_outflow, bc_arr);

            Gpu::streamSynchronize();  // otherwise we might be using too much memory
        }
    }
}

void
Godunov::ComputeAdvectiveVel ( Box const& xbx,
                               Box const& ybx,
                               Array4<Real> const& u_ad,
                               Array4<Real> const& v_ad,
                               Array4<Real const> const& Imx,
                               Array4<Real const> const& Imy,
                               Array4<Real const> const& Ipx,
                               Array4<Real const> const& Ipy,
                               Array4<Real const> const& vel,
                               Array4<Real const> const& f,
                               const Box& domain,
                               Real l_dt,
                               BCRec  const* pbc,
                               bool l_use_forces_in_trans,
                               Array4<int const> const& bc_arr)
{
    const Dim3 dlo = amrex::lbound(domain);
    const Dim3 dhi = amrex::ubound(domain);

    amrex::ParallelFor(xbx, ybx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        // We only care about x-velocity on x-faces here
        constexpr int n = 0;

        Real lo = Ipx(i-1,j,k,n);
        Real hi = Imx(i  ,j,k,n);

        if (l_use_forces_in_trans)
        {
            lo += 0.5*l_dt*f(i-1,j,k,n);
            hi += 0.5*l_dt*f(i  ,j,k,n);
        }

        const auto bc = HydroBC::getBC(i, j, k, n, domain, pbc, bc_arr);
        HydroBC::SetExtrapVelBCsLo(0, i, j, k, n, vel, lo, hi, bc.lo(0), dlo.x);
        HydroBC::SetExtrapVelBCsHi(0, i, j, k, n, vel, lo, hi, bc.hi(0), dhi.x);

        Real st = ( (lo+hi) >= 0.) ? lo : hi;
        bool ltm = ( (lo <= 0. && hi >= 0.) || (amrex::Math::abs(lo+hi) < small_vel) );
        u_ad(i,j,k) = ltm ? 0. : st;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        // We only care about y-velocity on y-faces here
        constexpr int n = 1;

        Real lo = Ipy(i,j-1,k,n);
        Real hi = Imy(i,j  ,k,n);

        if (l_use_forces_in_trans)
        {
            lo += 0.5*l_dt*f(i,j-1,k,n);
            hi += 0.5*l_dt*f(i,j  ,k,n);
        }

        const auto bc = HydroBC::getBC(i, j, k, n, domain, pbc, bc_arr);
        HydroBC::SetExtrapVelBCsLo(1, i, j, k, n, vel, lo, hi, bc.lo(1), dlo.y);
        HydroBC::SetExtrapVelBCsHi(1, i, j, k, n, vel, lo, hi, bc.hi(1), dhi.y);

        Real st = ( (lo+hi) >= 0.) ? lo : hi;
        bool ltm = ( (lo <= 0. && hi >= 0.) || (amrex::Math::abs(lo+hi) < small_vel) );
        v_ad(i,j,k) = ltm ? 0. : st;
    }
    );
}

void
Godunov::ExtrapVelToFacesOnBox (Box const& bx, int ncomp,
                                Box const& xbx,
                                Box const& ybx,
                                Array4<Real> const& qx,
                                Array4<Real> const& qy,
                                Array4<Real const> const& q,
                                Array4<Real const> const& u_ad,
                                Array4<Real const> const& v_ad,
                                Array4<Real> const& Imx,
                                Array4<Real> const& Imy,
                                Array4<Real> const& Ipx,
                                Array4<Real> const& Ipy,
                                Array4<Real const> const& f,
                                const Box& domain,
                                const Real* dx_arr,
                                Real l_dt,
                                BCRec  const* pbc,
                                bool l_use_forces_in_trans,
                                Real* p,
                                bool allow_inflow_on_outflow,
                                Array4<int const> const& bc_arr)
{

    const Dim3 dlo = amrex::lbound(domain);
    const Dim3 dhi = amrex::ubound(domain);

    Real dx = dx_arr[0];
    Real dy = dx_arr[1];

    Box xebox = Box(bx).grow(1,1).surroundingNodes(0);
    Box yebox = Box(bx).grow(0,1).surroundingNodes(1);

    Array4<Real> xlo = makeArray4(p, xebox, ncomp);
    p += xlo.size();
    Array4<Real> xhi = makeArray4(p, xebox, ncomp);
    p += xhi.size();
    Array4<Real> ylo = makeArray4(p, yebox, ncomp);
    p += ylo.size();
    Array4<Real> yhi = makeArray4(p, yebox, ncomp);

    // Don't need to save an intermediate upwinded edgestate in 2D like we do
    // in 3D. Could potentially combine these loops with making yz/xzlo, but
    // then would need to rework scratch, because we can't reuse the space in
    // Ipy in that case...
    amrex::ParallelFor(
    xebox, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        Real lo = Ipx(i-1,j,k,n);
        Real hi = Imx(i  ,j,k,n);

        if (l_use_forces_in_trans)
        {
            lo += 0.5*l_dt*f(i-1,j,k,n);
            hi += 0.5*l_dt*f(i  ,j,k,n);
        }

        const auto bc = HydroBC::getBC(i, j, k, n, domain, pbc, bc_arr);
        HydroBC::SetExtrapVelBCsLo(0, i, j, k, n, q, lo, hi, bc.lo(0), dlo.x);
        HydroBC::SetExtrapVelBCsHi(0, i, j, k, n, q, lo, hi, bc.hi(0), dhi.x);

        xlo(i,j,k,n) = lo;
        xhi(i,j,k,n) = hi;
    },
    yebox, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        Real lo = Ipy(i,j-1,k,n);
        Real hi = Imy(i,j  ,k,n);

        if (l_use_forces_in_trans)
        {
            lo += 0.5*l_dt*f(i,j-1,k,n);
            hi += 0.5*l_dt*f(i,j  ,k,n);
        }

        const auto bc = HydroBC::getBC(i, j, k, n, domain, pbc, bc_arr);
        HydroBC::SetExtrapVelBCsLo(1, i, j, k, n, q, lo, hi, bc.lo(1), dlo.y);
        HydroBC::SetExtrapVelBCsHi(1, i, j, k, n, q, lo, hi, bc.hi(1), dhi.y);

        ylo(i,j,k,n) = lo;
        yhi(i,j,k,n) = hi;
    }
    );

    //
    // X-Flux
    //
    Box const  xbxtmp = Box(xbx).enclosedCells().grow(0,1);
    Array4<Real> yzlo = makeArray4(Ipy.dataPtr(), amrex::surroundingNodes(xbxtmp,1), 1);

    amrex::ParallelFor(Box(yzlo),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        constexpr int n = 0;
        const auto bc = HydroBC::getBC(i, j, k, n, domain, pbc, bc_arr);
        Real l_yzlo, l_yzhi;

        l_yzlo = ylo(i,j,k,n);
        l_yzhi = yhi(i,j,k,n);

        HydroBC::SetExtrapVelBCsLo(1, i, j, k, n, q, l_yzlo, l_yzhi, bc.lo(1), dlo.y);
        HydroBC::SetExtrapVelBCsHi(1, i, j, k, n, q, l_yzlo, l_yzhi, bc.hi(1), dhi.y);

        Real vad = v_ad(i,j,k);
        Real st = (vad >= 0.) ? l_yzlo : l_yzhi;
        Real fu = (amrex::Math::abs(vad) < small_vel) ? 0.0 : 1.0;
        yzlo(i,j,k) = fu*st + (1.0 - fu) * 0.5 * (l_yzhi + l_yzlo);
    });

    amrex::ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        constexpr int n = 0;
        const auto bc = HydroBC::getBC(i, j, k, n, domain, pbc, bc_arr);
        Real stl = xlo(i,j,k,n) - (0.25*l_dt/dy)*(v_ad(i-1,j+1,k  )+v_ad(i-1,j,k))*
                                                 (yzlo(i-1,j+1,k  )-yzlo(i-1,j,k));
        Real sth = xhi(i,j,k,n) - (0.25*l_dt/dy)*(v_ad(i  ,j+1,k  )+v_ad(i  ,j,k))*
                                                 (yzlo(i  ,j+1,k  )-yzlo(i  ,j,k));

        if (!l_use_forces_in_trans)
        {
            stl += 0.5 * l_dt * f(i-1,j,k,n);
            sth += 0.5 * l_dt * f(i  ,j,k,n);
        }

        HydroBC::SetExtrapVelBCsLo(0, i, j, k, n, q, stl, sth, bc.lo(0), dlo.x);
        HydroBC::SetExtrapVelBCsHi(0, i, j, k, n, q, stl, sth, bc.hi(0), dhi.x);

        if (!allow_inflow_on_outflow) {
            if ( (i==dlo.x) && (bc.lo(0) == BCType::foextrap || bc.lo(0) == BCType::hoextrap) )
            {
                sth = amrex::min(sth,0.0_rt);
                stl = sth;
            }
            if ( (i==dhi.x+1) && (bc.hi(0) == BCType::foextrap || bc.hi(0) == BCType::hoextrap) )
            {
                 stl = amrex::max(stl,0.0_rt);
                 sth = stl;
            }
        }

        Real st = ( (stl+sth) >= 0.) ? stl : sth;
        bool ltm = ( (stl <= 0. && sth >= 0.) || (amrex::Math::abs(stl+sth) < small_vel) );
        qx(i,j,k) = ltm ? 0. : st;
    });

    //
    // Y-Flux
    //
    Box const ybxtmp  = Box(ybx).enclosedCells().grow(1,1);
    Array4<Real> xzlo = makeArray4(Ipy.dataPtr(), amrex::surroundingNodes(ybxtmp,0), 1);
    amrex::ParallelFor(Box(xzlo),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        constexpr int n = 1;
        const auto bc = HydroBC::getBC(i, j, k, n, domain, pbc, bc_arr);
        Real l_xzlo, l_xzhi;

        l_xzlo = xlo(i,j,k,n);
        l_xzhi = xhi(i,j,k,n);

        HydroBC::SetExtrapVelBCsLo(0,i, j, k, n, q, l_xzlo, l_xzhi, bc.lo(0), dlo.x);
        HydroBC::SetExtrapVelBCsHi(0,i, j, k, n, q, l_xzlo, l_xzhi, bc.hi(0), dhi.x);

        Real uad = u_ad(i,j,k);
        Real st = (uad >= 0.) ? l_xzlo : l_xzhi;
        Real fu = (amrex::Math::abs(uad) < small_vel) ? 0.0 : 1.0;
        xzlo(i,j,k) = fu*st + (1.0 - fu) * 0.5 * (l_xzhi + l_xzlo);
    });

    amrex::ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        constexpr int n = 1;
        const auto bc = HydroBC::getBC(i, j, k, n, domain, pbc, bc_arr);

        Real stl = ylo(i,j,k,n) - (0.25*l_dt/dx)*(u_ad(i+1,j-1,k  )+u_ad(i,j-1,k))*
                                                 (xzlo(i+1,j-1,k  )-xzlo(i,j-1,k));
        Real sth = yhi(i,j,k,n) - (0.25*l_dt/dx)*(u_ad(i+1,j  ,k  )+u_ad(i,j  ,k))*
                                                 (xzlo(i+1,j  ,k  )-xzlo(i,j  ,k));

        if (!l_use_forces_in_trans)
        {
           stl += 0.5 * l_dt * f(i,j-1,k,n);
           sth += 0.5 * l_dt * f(i,j  ,k,n);
        }

        HydroBC::SetExtrapVelBCsLo(1, i, j, k, n, q, stl, sth, bc.lo(1), dlo.y);
        HydroBC::SetExtrapVelBCsHi(1, i, j, k, n, q, stl, sth, bc.hi(1), dhi.y);

        if (!allow_inflow_on_outflow) {
            if ( (j==dlo.y) && (bc.lo(1) == BCType::foextrap || bc.lo(1) == BCType::hoextrap) )
            {
                sth = amrex::min(sth,0.0_rt);
                stl = sth;
            }
            if ( (j==dhi.y+1) && (bc.hi(1) == BCType::foextrap || bc.hi(1) == BCType::hoextrap) )
            {
                stl = amrex::max(stl,0.0_rt);
                sth = stl;
            }
        }

        Real st = ( (stl+sth) >= 0.) ? stl : sth;
        bool ltm = ( (stl <= 0. && sth >= 0.) || (amrex::Math::abs(stl+sth) < small_vel) );
        qy(i,j,k) = ltm ? 0. : st;
    });
}
/** @} */

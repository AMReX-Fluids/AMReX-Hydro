/** \addtogroup Utilities
 * @{
 */

#include <hydro_utils.H>

using namespace amrex;


void
HydroUtils::ComputeFluxes ( Box const& bx,
                            AMREX_D_DECL( Array4<Real> const& fx,
                                          Array4<Real> const& fy,
                                          Array4<Real> const& fz),
                            AMREX_D_DECL( Array4<Real const> const& umac,
                                          Array4<Real const> const& vmac,
                                          Array4<Real const> const& wmac),
                            AMREX_D_DECL( Array4<Real const> const& xed,
                                          Array4<Real const> const& yed,
                                          Array4<Real const> const& zed),
                            Geometry const& geom, const int ncomp,
                            const bool fluxes_are_area_weighted )
{

    const auto dx = geom.CellSizeArray();

    GpuArray<Real,AMREX_SPACEDIM> area;
#if ( AMREX_SPACEDIM == 3 )
    if (fluxes_are_area_weighted)
    {
        area[0] = dx[1]*dx[2];
        area[1] = dx[0]*dx[2];
        area[2] = dx[0]*dx[1];
    } else {
        area[0] = 1.; area[1] = 1.; area[2] = 1.;
    }
#else
    if (fluxes_are_area_weighted)
    {
        area[0] = dx[1];
        area[1] = dx[0];
    } else {
        area[0] = 1.; area[1] = 1.;
    }
#endif


    //
    //  X flux
    //
    const Box& xbx = amrex::surroundingNodes(bx,0);

    amrex::ParallelFor(xbx, ncomp, [fx, umac, xed, area]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        fx(i,j,k,n) = xed(i,j,k,n) * umac(i,j,k) * area[0];
    });

    //
    //  y flux
    //
    const Box& ybx = amrex::surroundingNodes(bx,1);

    amrex::ParallelFor(ybx, ncomp, [fy, vmac, yed, area]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        fy(i,j,k,n) = yed(i,j,k,n) * vmac(i,j,k) * area[1];
    });

#if (AMREX_SPACEDIM==3)
    //
    //  z flux
    //
    const Box& zbx = amrex::surroundingNodes(bx,2);

    amrex::ParallelFor(zbx, ncomp, [fz, wmac, zed, area]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        fz(i,j,k,n) = zed(i,j,k,n) * wmac(i,j,k) * area[2];
    });

#endif

}



void
HydroUtils::ComputeDivergence ( Box const& bx,
                                Array4<Real> const& div,
                                AMREX_D_DECL( Array4<Real const> const& fx,
                                              Array4<Real const> const& fy,
                                              Array4<Real const> const& fz),
                                const int ncomp, Geometry const& geom,
                                const Real mult,
                                const bool fluxes_are_area_weighted )
{

    const auto dxinv = geom.InvCellSizeArray();

    AMREX_D_TERM(Real fact_x = mult;,
                 Real fact_y = mult;,
                 Real fact_z = mult;);

    if (fluxes_are_area_weighted)
    {
        Real qvol;

#if (AMREX_SPACEDIM==3)
        qvol = dxinv[0] * dxinv[1] * dxinv[2];
#else
        qvol = dxinv[0] * dxinv[1];
#endif

    AMREX_D_TERM(fact_x *= qvol;,
                 fact_y *= qvol;,
                 fact_z *= qvol;);
    }
    else
    {
        AMREX_D_TERM(fact_x *= dxinv[0];,
                     fact_y *= dxinv[1];,
                     fact_z *= dxinv[2];);
    }

    amrex::ParallelFor(bx, ncomp,[=]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        div(i,j,k,n) =
              fact_x * ( fx(i+1,j,k,n) - fx(i,j,k,n) )
            + fact_y * ( fy(i,j+1,k,n) - fy(i,j,k,n) )
#if (AMREX_SPACEDIM==3)
            + fact_z * ( fz(i,j,k+1,n) - fz(i,j,k,n) )
#endif
            ;
    });
}

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//   RZ routines                                                         //
//                                                                       //
///////////////////////////////////////////////////////////////////////////
#if (AMREX_SPACEDIM==2)


void
HydroUtils::ComputeFluxesRZ ( Box const& bx,
                              Array4<Real> const& fx,
                              Array4<Real> const& fy,
                              Array4<Real const> const& umac,
                              Array4<Real const> const& vmac,
                              Array4<Real const> const& xed,
                              Array4<Real const> const& yed,
                              Array4<Real const> const& areax,
                              Array4<Real const> const& areay,
                              const int ncomp,
                              const bool fluxes_are_area_weighted )
{
    //
    //  X flux
    //
    const Box& xbx = amrex::surroundingNodes(bx,0);

    amrex::ParallelFor(xbx, ncomp, [fx, umac, xed, areax, fluxes_are_area_weighted]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (fluxes_are_area_weighted)
            fx(i,j,k,n) = xed(i,j,k,n) * umac(i,j,k) * areax(i,j,k);
        else
            fx(i,j,k,n) = xed(i,j,k,n) * umac(i,j,k);
    });

    //
    //  y flux
    //
    const Box& ybx = amrex::surroundingNodes(bx,1);

    amrex::ParallelFor(ybx, ncomp, [fy, vmac, yed, areay, fluxes_are_area_weighted]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (fluxes_are_area_weighted)
            fy(i,j,k,n) = yed(i,j,k,n) * vmac(i,j,k) * areay(i,j,k);
        else
            fy(i,j,k,n) = yed(i,j,k,n) * vmac(i,j,k);
    });
}

void
HydroUtils::ComputeDivergenceRZ ( Box const& bx,
                                  Array4<Real> const& div,
                                  Array4<Real const> const& fx,
                                  Array4<Real const> const& fy,
                                  Array4<Real const> const& vol,
                                  const int ncomp,
                                  const Real mult,
                                  const bool fluxes_are_area_weighted )
{
    AMREX_ALWAYS_ASSERT(fluxes_are_area_weighted);

    amrex::ParallelFor(bx, ncomp,[=]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (fluxes_are_area_weighted)
            div(i,j,k,n) = ( fx(i+1,j,k,n) -  fx(i,j,k,n) +
                             fy(i,j+1,k,n) -  fy(i,j,k,n) ) * mult / vol(i,j,k);
        else
            amrex::Abort("RZ Divergence not implemented for fluxes not area-weighted");

    });
}

#endif


///////////////////////////////////////////////////////////////////////////
//                                                                       //
//   EB routines                                                         //
//                                                                       //
///////////////////////////////////////////////////////////////////////////
#ifdef AMREX_USE_EB

void
HydroUtils::EB_ComputeDivergence ( Box const& bx,
                                   Array4<Real> const& div,
                                   AMREX_D_DECL( Array4<Real const> const& fx,
                                                 Array4<Real const> const& fy,
                                                 Array4<Real const> const& fz),
                                   Array4<Real const> const& vfrac,
                                   const int ncomp, Geometry const& geom,
                                   const Real mult,
                                   const bool fluxes_are_area_weighted )
{
    const auto dxinv = geom.InvCellSizeArray();

#if (AMREX_SPACEDIM==3)
    Real qvol = dxinv[0] * dxinv[1] * dxinv[2];
#else
    Real qvol = dxinv[0] * dxinv[1];
#endif

    amrex::ParallelFor(bx, ncomp, [=]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if ( vfrac(i,j,k) > 0.)
        {
            if (fluxes_are_area_weighted)
                div(i,j,k,n) =  mult * qvol / vfrac(i,j,k) *
                    (
                             fx(i+1,j,k,n) -  fx(i,j,k,n)
                           + fy(i,j+1,k,n) -  fy(i,j,k,n)
#if (AMREX_SPACEDIM==3)
                           + fz(i,j,k+1,n) -  fz(i,j,k,n)
#endif
                    );
            else
                div(i,j,k,n) =  mult / vfrac(i,j,k) *
                    (
                             (fx(i+1,j,k,n) -  fx(i,j,k,n)) * dxinv[0]
                           + (fy(i,j+1,k,n) -  fy(i,j,k,n)) * dxinv[1]
#if (AMREX_SPACEDIM==3)
                           + (fz(i,j,k+1,n) -  fz(i,j,k,n)) * dxinv[1]
#endif
                    );
        }
        else
        {
            div(i,j,k,n) = 0.0;
        }

    });
}


void
HydroUtils::EB_ComputeFluxes ( Box const& bx,
                               AMREX_D_DECL( Array4<Real> const& fx,
                                             Array4<Real> const& fy,
                                             Array4<Real> const& fz),
                               AMREX_D_DECL( Array4<Real const> const& umac,
                                             Array4<Real const> const& vmac,
                                             Array4<Real const> const& wmac),
                               AMREX_D_DECL( Array4<Real const> const& xed,
                                             Array4<Real const> const& yed,
                                             Array4<Real const> const& zed),
                               AMREX_D_DECL( Array4<Real const> const& apx,
                                             Array4<Real const> const& apy,
                                             Array4<Real const> const& apz),
                               Geometry const& geom, const int ncomp,
                               Array4<EBCellFlag const> const& flag,
                               const bool fluxes_are_area_weighted )
{

    const auto dx = geom.CellSizeArray();

    GpuArray<Real,AMREX_SPACEDIM> area;

#if ( AMREX_SPACEDIM == 3 )
    if (fluxes_are_area_weighted)
    {
        area[0] = dx[1]*dx[2];
        area[1] = dx[0]*dx[2];
        area[2] = dx[0]*dx[1];
    } else {
        area[0] = 1.; area[1] = 1.; area[2] = 1.;
    }
#else
    if (fluxes_are_area_weighted)
    {
        area[0] = dx[1];
        area[1] = dx[0];
    } else {
        area[0] = 1.; area[1] = 1.;
    }
#endif

    //
    //  X flux
    //
    const Box& xbx = amrex::surroundingNodes(bx,0);

    amrex::ParallelFor(xbx, ncomp, [fx, umac, xed, area, apx, flag]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (flag(i,j,k).isConnected(-1,0,0))
            fx(i,j,k,n) = xed(i,j,k,n) * umac(i,j,k) * apx(i,j,k) * area[0];
        else
            fx(i,j,k,n) = 0.;
    });

    //
    //  y flux
    //
    const Box& ybx = amrex::surroundingNodes(bx,1);

    amrex::ParallelFor(ybx, ncomp, [fy, vmac, yed, area, apy, flag]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (flag(i,j,k).isConnected(0,-1,0))
            fy(i,j,k,n) = yed(i,j,k,n) * vmac(i,j,k) * apy(i,j,k) * area[1];
        else
            fy(i,j,k,n) = 0.;
    });

#if (AMREX_SPACEDIM==3)
    //
    //  z flux
    //
    const Box& zbx = amrex::surroundingNodes(bx,2);

    amrex::ParallelFor(zbx, ncomp, [fz, wmac, zed, area, apz, flag]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (flag(i,j,k).isConnected(0,0,-1))
            fz(i,j,k,n) = zed(i,j,k,n) * wmac(i,j,k) * apz(i,j,k) * area[2];
        else
            fz(i,j,k,n) = 0.;
    });
#endif

}

#endif
/** @}*/

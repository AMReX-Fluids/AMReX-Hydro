/** \addtogroup Utilities
 * @{
 */

#include <hydro_utils.H>

#ifdef AMREX_USE_EB
#include <AMReX_MultiCutFab.H>
#endif

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
#if (AMREX_SPACEDIM == 2)
    if (geom.IsRZ()) {
        // Need metrics when using RZ
        Array<FArrayBox,AMREX_SPACEDIM> area;
        if (fluxes_are_area_weighted) {
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
            {
                const Box ebx = surroundingNodes(bx,dir);
                area[dir].resize(ebx,1,The_Async_Arena());
                geom.SetFaceArea(area[dir],ebx,dir);
            }
        }
        const auto& ax = (fluxes_are_area_weighted) ? area[0].const_array()
                                                    : Array4<const Real>{};
        const auto& ay = (fluxes_are_area_weighted) ? area[1].const_array()
                                                    : Array4<const Real>{};
        //
        //  X flux
        //
        const Box& xbx = amrex::surroundingNodes(bx,0);
        amrex::ParallelFor(xbx, ncomp, [fx, umac, xed, ax, fluxes_are_area_weighted]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            if (fluxes_are_area_weighted) {
                fx(i,j,k,n) = xed(i,j,k,n) * umac(i,j,k) * ax(i,j,k);
            } else {
                fx(i,j,k,n) = xed(i,j,k,n) * umac(i,j,k);
            }
        });

        //
        //  Y flux
        //
        const Box& ybx = amrex::surroundingNodes(bx,1);
        amrex::ParallelFor(ybx, ncomp, [fy, vmac, yed, ay, fluxes_are_area_weighted]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            if (fluxes_are_area_weighted) {
                fy(i,j,k,n) = yed(i,j,k,n) * vmac(i,j,k) * ay(i,j,k);
            } else {
                fy(i,j,k,n) = yed(i,j,k,n) * vmac(i,j,k);
            }
        });
    } else
#endif
    {
        const auto dx = geom.CellSizeArray();

        GpuArray<Real,AMREX_SPACEDIM> area;
#if ( AMREX_SPACEDIM == 3 )
        if (fluxes_are_area_weighted) {
            area[0] = dx[1]*dx[2];
            area[1] = dx[0]*dx[2];
            area[2] = dx[0]*dx[1];
        } else {
            area[0] = 1.; area[1] = 1.; area[2] = 1.;
        }
#else
        if (fluxes_are_area_weighted) {
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
        //  Y flux
        //
        const Box& ybx = amrex::surroundingNodes(bx,1);
        amrex::ParallelFor(ybx, ncomp, [fy, vmac, yed, area]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            fy(i,j,k,n) = yed(i,j,k,n) * vmac(i,j,k) * area[1];
        });

#if ( AMREX_SPACEDIM == 3 )
        //
        //  Z flux
        //
        const Box& zbx = amrex::surroundingNodes(bx,2);
        amrex::ParallelFor(zbx, ncomp, [fz, wmac, zed, area]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            fz(i,j,k,n) = zed(i,j,k,n) * wmac(i,j,k) * area[2];
        });
#endif
    }

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
#if (AMREX_SPACEDIM == 2)
    if (geom.IsRZ())
    {
        FArrayBox vol_fab(bx,1,The_Async_Arena());
        geom.SetVolume(vol_fab,bx);
        Array<FArrayBox,AMREX_SPACEDIM> area;
        if (!fluxes_are_area_weighted) {
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
            {
                const Box ebx = surroundingNodes(bx,dir);
                area[dir].resize(ebx,1,The_Async_Arena());
                geom.SetFaceArea(area[dir],ebx,dir);
            }
        }
        const auto& vol = vol_fab.const_array();
        const auto& ax  = (fluxes_are_area_weighted) ? Array4<const Real>{}
                                                     : area[0].const_array();
        const auto& ay  = (fluxes_are_area_weighted) ? Array4<const Real>{}
                                                      : area[1].const_array();
        amrex::ParallelFor(bx, ncomp,[=]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            if (fluxes_are_area_weighted) {
                div(i,j,k,n) = ( fx(i+1,j,k,n) -  fx(i,j,k,n) +
                                 fy(i,j+1,k,n) -  fy(i,j,k,n) ) * mult / vol(i,j,k);
            } else {
                div(i,j,k,n) = ( ax(i+1,j,k)*fx(i+1,j,k,n) -  ax(i,j,k)*fx(i,j,k,n) +
                                 ay(i,j+1,k)*fy(i,j+1,k,n) -  ay(i,j,k)*fy(i,j,k,n) ) * mult / vol(i,j,k);
            }
        });
    } else
#endif
    {
        const auto dxinv = geom.InvCellSizeArray();
        AMREX_D_TERM(Real fact_x = mult;,
                     Real fact_y = mult;,
                     Real fact_z = mult;);

        if (fluxes_are_area_weighted) {
            Real qvol = AMREX_D_TERM(dxinv[0],*dxinv[1],*dxinv[2]);

            AMREX_D_TERM(fact_x *= qvol;,
                         fact_y *= qvol;,
                         fact_z *= qvol;);
        } else {
            AMREX_D_TERM(fact_x *= dxinv[0];,
                         fact_y *= dxinv[1];,
                         fact_z *= dxinv[2];);
        }

        amrex::ParallelFor(bx, ncomp,[=]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            div(i,j,k,n) = AMREX_D_TERM(  fact_x * ( fx(i+1,j,k,n) - fx(i,j,k,n) ),
                                        + fact_y * ( fy(i,j+1,k,n) - fy(i,j,k,n) ),
                                        + fact_z * ( fz(i,j,k+1,n) - fz(i,j,k,n) ));
        });
    }
}

void
HydroUtils::ComputeConvectiveTerm(Box const& bx, int num_comp, MFIter& mfi,
                                  Array4<Real const> const& q,
                                  AMREX_D_DECL(Array4<Real const> const& q_on_face_x,
                                               Array4<Real const> const& q_on_face_y,
                                               Array4<Real const> const& q_on_face_z),
                                  Array4<Real const> const& divu,
                                  Array4<Real> const& convTerm,
                                  int const* iconserv,
#ifdef AMREX_USE_EB
                                  const EBFArrayBoxFactory& ebfact,
#endif
                                  std::string& advection_type)
{
    //
    // If convective, we define convTerm = u dot grad q = div (u q) - q div(u)
    //
    if (advection_type == "MOL")
    {
        // Here we use q at the same time as the velocity
        amrex::ParallelFor(bx, num_comp, [=]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            if (!iconserv[n])
                convTerm(i,j,k,n) += q(i,j,k,n)*divu(i,j,k);
        });
    }
    else if (advection_type == "Godunov" || advection_type == "BDS")
    {
        bool regular = true;
#ifdef AMREX_USE_EB
        EBCellFlagFab const& flagfab = ebfact.getMultiEBCellFlagFab()[mfi];
        regular = (flagfab.getType(bx) == FabType::regular);
#endif
        // Here we want to use q predicted to t^{n+1/2}
        if (regular)
        {
            amrex::ParallelFor(bx, num_comp, [=]
            AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                if (!iconserv[n])
                {
                    Real qavg  = q_on_face_x(i,j,k,n) + q_on_face_x(i+1,j,k,n)
                               + q_on_face_y(i,j,k,n) + q_on_face_y(i,j+1,k,n);
#if (AMREX_SPACEDIM == 2)
                    qavg *= Real(0.25);
#else
                    qavg += q_on_face_z(i,j,k,n) + q_on_face_z(i,j,k+1,n);
                    qavg /= 6.0;
#endif
                    // Note that because we define convTerm as MINUS div(u u), here we add u div (u)
                    convTerm(i,j,k,n) += qavg*divu(i,j,k);
                }
            });
        }
#ifdef AMREX_USE_EB
        else {
            if (flagfab.getType(bx) != FabType::covered) {
                auto const& vfrac_arr            = ebfact.getVolFrac().const_array(mfi);
                AMREX_D_TERM(auto const& apx_arr = ebfact.getAreaFrac()[0]->const_array(mfi);,
                             auto const& apy_arr = ebfact.getAreaFrac()[1]->const_array(mfi);,
                             auto const& apz_arr = ebfact.getAreaFrac()[2]->const_array(mfi););

                amrex::ParallelFor(bx, num_comp, [=]
                AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (!iconserv[n] && vfrac_arr(i,j,k) > 0.)
                    {
                        Real qavg  = apx_arr(i,j,k)*q_on_face_x(i,j,k,n) + apx_arr(i+1,j,k)*q_on_face_x(i+1,j,k,n);
                        qavg += apy_arr(i,j,k)*q_on_face_y(i,j,k,n) + apy_arr(i,j+1,k)*q_on_face_y(i,j+1,k,n);
#if (AMREX_SPACEDIM == 2)
                        qavg *= 1.0 / (apx_arr(i,j,k) + apx_arr(i+1,j,k) + apy_arr(i,j,k) + apy_arr(i,j+1,k));
#else
                        qavg += apz_arr(i,j,k)*q_on_face_z(i,j,k,n) + apz_arr(i,j,k+1)*q_on_face_z(i,j,k+1,n);
                        qavg *= 1.0 / ( apx_arr(i,j,k) + apx_arr(i+1,j,k) + apy_arr(i,j,k) + apy_arr(i,j+1,k)
                                        +apz_arr(i,j,k) + apz_arr(i,j,k+1) );
#endif

                        // Note that because we define convTerm as MINUS div(u u), here we add u div (u)
                        convTerm(i,j,k,n) += qavg*divu(i,j,k);
                    }
                });
            }
        }
#endif
    }
    else {
        amrex::Abort("HydroUtils::ComputeConvectiveTerm: unknown advection_type ");
    }
}

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
                           + (fz(i,j,k+1,n) -  fz(i,j,k,n)) * dxinv[2]
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
HydroUtils::EB_ComputeDivergence ( Box const& bx,
                                   Array4<Real> const& div,
                                   AMREX_D_DECL( Array4<Real const> const& fx,
                                                 Array4<Real const> const& fy,
                                                 Array4<Real const> const& fz),
                                   Array4<Real const> const& vfrac,
                                   const int ncomp, Geometry const& geom,
                                   const Real mult,
                                   const bool fluxes_are_area_weighted,
                                   Array4<Real const> const& eb_velocity,
                                   Array4<Real const> const& values_on_eb_inflow,
                                   Array4<EBCellFlag const> const& flag_arr,
                                   Array4<Real const> const& barea,
                                   Array4<Real const> const& bnorm)
{

    // Compute the standard EB divergence term
    EB_ComputeDivergence(bx, div, AMREX_D_DECL(fx, fy, fz),
                         vfrac, ncomp, geom, mult, fluxes_are_area_weighted);

    if( eb_velocity ) {

        AMREX_ASSERT( values_on_eb_inflow );

        const auto &dxinv = geom.InvCellSizeArray();

        amrex::ParallelFor(bx, ncomp, [div,eb_velocity,values_on_eb_inflow,
          flag_arr,barea,vfrac,bnorm,dxinv,mult]
          AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
         if (flag_arr(i,j,k).isSingleValued()) {

           Real Ueb_dot_n = (AMREX_D_TERM(
                                eb_velocity(i,j,k,0)*bnorm(i,j,k,0),
                              + eb_velocity(i,j,k,1)*bnorm(i,j,k,1),
                              + eb_velocity(i,j,k,2)*bnorm(i,j,k,2)));

           div(i,j,k,n) += (mult / vfrac(i,j,k)) *
              values_on_eb_inflow(i,j,k,n) * Ueb_dot_n * barea(i,j,k) * dxinv[0];

         }
        });
     }
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

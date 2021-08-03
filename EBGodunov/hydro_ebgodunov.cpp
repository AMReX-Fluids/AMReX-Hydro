/**
 * \file hydro_ebgodunov.cpp
 * \addtogroup EBGodunov
 * @{
 *
 */

#include <hydro_ebgodunov.H>
#include <hydro_godunov.H>
#include <hydro_redistribution.H>
#include <hydro_utils.H>
#include <hydro_constants.H>

using namespace amrex;


void
EBGodunov::ComputeAofs ( MultiFab& aofs, const int aofs_comp, const int ncomp,
                         MultiFab const& state, const int state_comp,
                         AMREX_D_DECL( MultiFab const& umac,
                                       MultiFab const& vmac,
                                       MultiFab const& wmac),
                         AMREX_D_DECL( MultiFab& xedge,
                                       MultiFab& yedge,
                                       MultiFab& zedge),
                         const int  edge_comp,
                         const bool known_edgestate,
                         AMREX_D_DECL( MultiFab& xfluxes,
                                       MultiFab& yfluxes,
                                       MultiFab& zfluxes),
                         int fluxes_comp,
                         MultiFab const& fq,
                         const int fq_comp,
                         MultiFab const& divu,
                         Vector<BCRec> const& h_bc,
                         BCRec const* d_bc,
                         Geometry const& geom,
                         Gpu::DeviceVector<int>& iconserv,
                         const Real dt,
                         const bool is_velocity,
                         std::string redistribution_type)
{
    BL_PROFILE("EBGodunov::ComputeAofs()");

    bool fluxes_are_area_weighted = true;
    int const* iconserv_ptr = iconserv.data();

    AMREX_ALWAYS_ASSERT(state.hasEBFabFactory());

    auto const& ebfact= dynamic_cast<EBFArrayBoxFactory const&>(state.Factory());
    auto const& flags = ebfact.getMultiEBCellFlagFab();
    auto const& fcent = ebfact.getFaceCent();
    auto const& ccent = ebfact.getCentroid();
    auto const& vfrac = ebfact.getVolFrac();
    auto const& areafrac = ebfact.getAreaFrac();

    // Create temporary holder for advection term. Needed so we can call FillBoundary.
    MultiFab advc(state.boxArray(),state.DistributionMap(),ncomp,3,MFInfo(),ebfact);
    advc.setVal(0.);

    // if we need convective form, we must also compute
    // div(u_mac)
    MultiFab divu_mac(state.boxArray(),state.DistributionMap(),1,4, MFInfo(), ebfact);
    for (long unsigned i = 0; i < iconserv.size(); ++i)
    {
        if (!iconserv[i])
        {
            Array<MultiFab const*,AMREX_SPACEDIM> u;
            AMREX_D_TERM(u[0] = &umac;,
                         u[1] = &vmac;,
                         u[2] = &wmac;);

            if (!ebfact.isAllRegular())
                amrex::EB_computeDivergence(divu_mac,u,geom,true);
            else
                amrex::computeDivergence(divu_mac,u,geom);

            divu_mac.FillBoundary(geom.periodicity());

            break;
        }
    }

    // Compute -div instead of computing div -- this is just for consistency
    // with the way we HAVE to do it for EB (because redistribution operates on
    // -div rather than div
    Real mult = -1.0;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(aofs, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {

        const Box& bx   = mfi.tilebox();

        auto const& flagfab = ebfact.getMultiEBCellFlagFab()[mfi];
	// A regular box uses 3 ghost cells:
	// We predict the state on the edge based box; that calls slopes on
	// i & i-1; slopes then looks at (i-1)-2 for 4th order slopes
	// => test on bx grow 3
        bool regular = (flagfab.getType(amrex::grow(bx,3)) == FabType::regular);

        // Get handlers to Array4
        //
        AMREX_D_TERM( const auto& fx = xfluxes.array(mfi,fluxes_comp);,
                      const auto& fy = yfluxes.array(mfi,fluxes_comp);,
                      const auto& fz = zfluxes.array(mfi,fluxes_comp););

        AMREX_D_TERM( const auto& xed = xedge.array(mfi,edge_comp);,
                      const auto& yed = yedge.array(mfi,edge_comp);,
                      const auto& zed = zedge.array(mfi,edge_comp););

        AMREX_D_TERM( const auto& u = umac.const_array(mfi);,
                      const auto& v = vmac.const_array(mfi);,
                      const auto& w = wmac.const_array(mfi););

	Array4<Real> advc_arr = advc.array(mfi);

        if (flagfab.getType(bx) == FabType::covered)
        {
	    AMREX_D_TERM( const Box& xbx = mfi.nodaltilebox(0);,
			  const Box& ybx = mfi.nodaltilebox(1);,
			  const Box& zbx = mfi.nodaltilebox(2); );

            auto const& aofs_arr = aofs.array(mfi, aofs_comp);

            amrex::ParallelFor(
                bx, ncomp, [aofs_arr] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                { aofs_arr( i, j, k, n ) = covered_val;},

                xbx, ncomp, [fx,xed] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                { fx( i, j, k, n ) = 0.0; xed( i, j, k, n ) = covered_val;},

                ybx, ncomp, [fy,yed] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                { fy( i, j, k, n ) = 0.0; yed( i, j, k, n ) = covered_val;});

#if (AMREX_SPACEDIM==3)
            amrex::ParallelFor(
                zbx, ncomp, [fz,zed]AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                { fz( i, j, k, n ) = 0.0; zed( i, j, k, n ) = covered_val;});
#endif
        }
        else
        {
        if (regular)   // Plain Godunov
        {
            if (not known_edgestate)
            {
                Godunov::ComputeEdgeState( bx, ncomp,
                                           state.array(mfi,state_comp),
                                           AMREX_D_DECL( xed, yed, zed ),
                                           AMREX_D_DECL( u, v, w ),
                                           divu.array(mfi),
                                           fq.array(mfi,fq_comp),
                                           geom, dt, d_bc,
                                           iconserv.data(),
                                           false,
                                           false,
                                           is_velocity );
            }

            HydroUtils::ComputeFluxes( bx,
                                       AMREX_D_DECL( fx, fy, fz ),
                                       AMREX_D_DECL( u, v, w ),
                                       AMREX_D_DECL( xed, yed, zed ),
                                       geom, ncomp, fluxes_are_area_weighted );

            HydroUtils::ComputeDivergence( bx, advc_arr,
                                           AMREX_D_DECL( fx, fy, fz ),
                                           ncomp, geom,
                                           mult, fluxes_are_area_weighted);

            // Compute the convective form if needed by accounting for extra term
            auto const& divu_arr  = divu_mac.array(mfi);
            amrex::ParallelFor(bx, ncomp, [=]
            AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                if (!iconserv_ptr[n])
                {
                    Real q = xed(i,j,k,n) + xed(i+1,j,k,n)
                           + yed(i,j,k,n) + yed(i,j+1,k,n);
#if (AMREX_SPACEDIM == 2)
                    q *= 0.25;
#else
                    q += zed(i,j,k,n) + zed(i,j,k+1,n);
                    q /= 6.0;
#endif
                    advc_arr(i,j,k,n) += q*divu_arr(i,j,k);
                }
            });

        }
        else     // EB Godunov
        {
            AMREX_D_TERM(Array4<Real const> const& fcx = fcent[0]->const_array(mfi);,
                         Array4<Real const> const& fcy = fcent[1]->const_array(mfi);,
                         Array4<Real const> const& fcz = fcent[2]->const_array(mfi););

            AMREX_D_TERM(Array4<Real const> const& apx = areafrac[0]->const_array(mfi);,
                         Array4<Real const> const& apy = areafrac[1]->const_array(mfi);,
                         Array4<Real const> const& apz = areafrac[2]->const_array(mfi););

            Array4<Real const> const& ccent_arr = ccent.const_array(mfi);
            Array4<Real const> const& vfrac_arr = vfrac.const_array(mfi);
            auto const& flags_arr  = flags.const_array(mfi);

	    //FIXME - compare to HydroUtils which hard codes 4 ghost cells for all
            int ngrow = 4;

            if (redistribution_type=="StateRedist")
                ++ngrow;

            FArrayBox tmpfab(amrex::grow(bx,ngrow),  (4*AMREX_SPACEDIM + 2)*ncomp);
            Elixir    eli = tmpfab.elixir();


            if (not known_edgestate)
            {
                EBGodunov::ComputeEdgeState( bx, ncomp,
                                             state.array(mfi,state_comp),
                                             AMREX_D_DECL( xed, yed, zed ),
                                             AMREX_D_DECL( u, v, w ),
                                             divu.array(mfi),
                                             fq.array(mfi,fq_comp),
                                             geom, dt, h_bc, d_bc,
                                             iconserv.data(),
                                             tmpfab.dataPtr(),
                                             flags_arr,
                                             AMREX_D_DECL( apx, apy, apz ),
                                             vfrac_arr,
                                             AMREX_D_DECL( fcx, fcy, fcz ),
                                             ccent_arr,
                                             is_velocity );
            }

            HydroUtils::EB_ComputeFluxes( bx,
                                          AMREX_D_DECL( fx, fy, fz ),
                                          AMREX_D_DECL( u, v, w ),
                                          AMREX_D_DECL( xed, yed, zed ),
                                          AMREX_D_DECL( apx, apy, apz ),
                                          geom, ncomp, flags_arr, fluxes_are_area_weighted );

            // div at ncomp*3 to make space for the 3 redistribute temporaries
            HydroUtils::EB_ComputeDivergence( bx,
                                              advc_arr,
                                              AMREX_D_DECL( fx, fy, fz ),
                                              vfrac_arr, ncomp, geom,
                                              mult, fluxes_are_area_weighted);

            // Compute the convective form if needed by accounting for extra term
            auto const& divu_arr  = divu_mac.array(mfi);
            amrex::ParallelFor(bx, ncomp, [=]
            AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                if (!iconserv_ptr[n])
                {
		  if ( vfrac_arr(i,j,k) != 0 )
                  {
                    Real q = xed(i,j,k,n)*apx(i,j,k) + xed(i+1,j,k,n)*apx(i+1,j,k)
                           + yed(i,j,k,n)*apy(i,j,k) + yed(i,j+1,k,n)*apy(i,j+1,k);
#if (AMREX_SPACEDIM == 2)
                    q /= (apx(i,j,k)+apx(i+1,j,k)+apy(i,j,k)+apy(i,j+1,k));
#else
                    q += zed(i,j,k,n)*apz(i,j,k) + zed(i,j,k+1,n)*apz(i,j,k+1);
                    q /= (apx(i,j,k)+apx(i+1,j,k)+apy(i,j,k)+apy(i,j+1,k)+apz(i,j,k)+apz(i,j,k+1));
#endif
                    advc_arr(i,j,k,n) += q*divu_arr(i,j,k);
		  }
                 }
            });
	  }
	}
    }

    advc.FillBoundary(geom.periodicity());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(aofs, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        auto const& bx = mfi.tilebox();

        auto const& flagfab = ebfact.getMultiEBCellFlagFab()[mfi];
        auto const& flag    = flagfab.const_array();
	auto const& aofs_arr = aofs.array(mfi, aofs_comp);
	auto const& advc_arr = advc.array(mfi);

        if (flagfab.getType(bx) != FabType::covered )
	{
	  // FIXME? not sure if 4 is really needed or if 3 could do
	  // But this is a safe choice
	  if (flagfab.getType(grow(bx,4)) != FabType::regular)
	  {
	    //
	    // Redistribute
	    //
	    AMREX_D_TERM( auto apx = ebfact.getAreaFrac()[0]->const_array(mfi);,
			  auto apy = ebfact.getAreaFrac()[1]->const_array(mfi);,
			  auto apz = ebfact.getAreaFrac()[2]->const_array(mfi); );

	    AMREX_D_TERM( Array4<Real const> fcx = ebfact.getFaceCent()[0]->const_array(mfi);,
			  Array4<Real const> fcy = ebfact.getFaceCent()[1]->const_array(mfi);,
			  Array4<Real const> fcz = ebfact.getFaceCent()[2]->const_array(mfi););

	    Array4<Real const> ccc = ebfact.getCentroid().const_array(mfi);
            Array4<Real const> const& vfrac_arr = vfrac.const_array(mfi);

	    // This is scratch space if calling StateRedistribute,
            //  but is used as the weights (here set to 1) if calling
            //  FluxRedistribute
	    Box gbx = bx;

	    if (redistribution_type == "StateRedist")
	      gbx.grow(3);
	    else if (redistribution_type == "FluxRedist")
	      gbx.grow(2);

	    FArrayBox tmpfab(gbx, ncomp);
	    Elixir eli = tmpfab.elixir();
            Array4<Real> scratch = tmpfab.array(0);
            if (redistribution_type == "FluxRedist")
            {
                amrex::ParallelFor(Box(scratch),
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                { scratch(i,j,k) = 1.;});
            }

	    Redistribution::Apply( bx, ncomp, aofs_arr, advc.array(mfi),
				   state.const_array(mfi, state_comp), scratch, flag,
                                   AMREX_D_DECL(apx,apy,apz), vfrac_arr,
				   AMREX_D_DECL(fcx,fcy,fcz), ccc, d_bc,
                                   geom, dt, redistribution_type );

	    // Change sign because we computed -div for all cases
	    amrex::ParallelFor(bx, ncomp, [aofs_arr]
	    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	    { aofs_arr( i, j, k, n ) *=  - 1.0; });
	  }
	  else
	  {
	    // Change sign because for EB we computed -div
            amrex::ParallelFor(bx, ncomp, [aofs_arr, advc_arr]
	    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	    { aofs_arr( i, j, k, n ) =  -advc_arr(i,j,k,n); });
	  }
	}
    }
}


void
EBGodunov::ComputeSyncAofs ( MultiFab& aofs, const int aofs_comp, const int ncomp,
                             MultiFab const& state, const int state_comp,
                             AMREX_D_DECL( MultiFab const& umac,
                                           MultiFab const& vmac,
                                           MultiFab const& wmac),
                             AMREX_D_DECL( MultiFab const& ucorr,
                                           MultiFab const& vcorr,
                                           MultiFab const& wcorr),
                             AMREX_D_DECL( MultiFab& xedge,
                                           MultiFab& yedge,
                                           MultiFab& zedge),
                             const int  edge_comp,
                             const bool known_edgestate,
                             AMREX_D_DECL( MultiFab& xfluxes,
                                           MultiFab& yfluxes,
                                           MultiFab& zfluxes),
                             int fluxes_comp,
                             MultiFab const& fq,
                             const int fq_comp,
                             MultiFab const& divu,
                             Vector<BCRec> const& h_bc,
                             BCRec const* d_bc,
                             Geometry const& geom,
                             Gpu::DeviceVector<int>& iconserv,
                             const Real dt,
                             const bool is_velocity,
                             std::string redistribution_type )
{
    BL_PROFILE("EBGodunov::ComputeSyncAofs()");

    bool fluxes_are_area_weighted = true;

    AMREX_ALWAYS_ASSERT(state.hasEBFabFactory());

    auto const& ebfact= dynamic_cast<EBFArrayBoxFactory const&>(state.Factory());
    auto const& flags = ebfact.getMultiEBCellFlagFab();
    auto const& fcent = ebfact.getFaceCent();
    auto const& ccent = ebfact.getCentroid();
    auto const& vfrac = ebfact.getVolFrac();
    auto const& areafrac = ebfact.getAreaFrac();

    // Create temporary holder for advection term. Needed so we can call FillBoundary.
    MultiFab advc(state.boxArray(),state.DistributionMap(),ncomp,3,MFInfo(),ebfact);
    advc.setVal(0.);

    // Compute -div instead of computing div -- this is just for consistency
    // with the way we HAVE to do it for EB (because redistribution operates on
    // -div rather than div
    Real mult = -1.0;


#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(aofs,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {

        const Box& bx   = mfi.tilebox();

        auto const& flagfab = ebfact.getMultiEBCellFlagFab()[mfi];
        bool regular = (flagfab.getType(amrex::grow(bx,3)) == FabType::regular);

        //
        // Get handlers to Array4
        //
        AMREX_D_TERM( const auto& fx = xfluxes.array(mfi,fluxes_comp);,
                      const auto& fy = yfluxes.array(mfi,fluxes_comp);,
                      const auto& fz = zfluxes.array(mfi,fluxes_comp););

        AMREX_D_TERM( const auto& xed = xedge.array(mfi,edge_comp);,
                      const auto& yed = yedge.array(mfi,edge_comp);,
                      const auto& zed = zedge.array(mfi,edge_comp););

        AMREX_D_TERM( const auto& uc = ucorr.const_array(mfi);,
                      const auto& vc = vcorr.const_array(mfi);,
                      const auto& wc = wcorr.const_array(mfi););

	Array4<Real> advc_arr = advc.array(mfi);

        if (flagfab.getType(bx) == FabType::covered)
        {
	    AMREX_D_TERM( const Box& xbx = mfi.nodaltilebox(0);,
			  const Box& ybx = mfi.nodaltilebox(1);,
			  const Box& zbx = mfi.nodaltilebox(2); );

            auto const& aofs_arr = aofs.array(mfi, aofs_comp);

            amrex::ParallelFor(
                bx, ncomp, [aofs_arr] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                { aofs_arr( i, j, k, n ) = covered_val;},

                xbx, ncomp, [fx,xed] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                { fx( i, j, k, n ) = 0.0; xed( i, j, k, n ) = covered_val;},

                ybx, ncomp, [fy,yed] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                { fy( i, j, k, n ) = 0.0; yed( i, j, k, n ) = covered_val;});

#if (AMREX_SPACEDIM==3)
            amrex::ParallelFor(
                zbx, ncomp, [fz,zed]AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                { fz( i, j, k, n ) = 0.0; zed( i, j, k, n ) = covered_val;});
#endif
        }
        else
        {
        if (regular) // Plain Godunov
        {
            if (not known_edgestate)
            {
                AMREX_D_TERM( const auto& u = umac.const_array(mfi);,
                              const auto& v = vmac.const_array(mfi);,
                              const auto& w = wmac.const_array(mfi););

                Godunov::ComputeEdgeState( bx, ncomp,
                                           state.array(mfi,state_comp),
                                           AMREX_D_DECL( xed, yed, zed ),
                                           AMREX_D_DECL( u, v, w ),
                                           divu.array(mfi),
                                           fq.array(mfi,fq_comp),
                                           geom, dt, d_bc,
                                           iconserv.data(),
                                           false,
                                           false,
                                           is_velocity );
            }



            HydroUtils::ComputeFluxes( bx,
                                       AMREX_D_DECL( fx, fy, fz ),
                                       AMREX_D_DECL( uc, vc, wc ),
                                       AMREX_D_DECL( xed, yed, zed ),
                                       geom, ncomp, fluxes_are_area_weighted );


            HydroUtils::ComputeDivergence( bx, advc_arr,
                                           AMREX_D_DECL( fx, fy, fz ),
                                           ncomp, geom,
                                           mult, fluxes_are_area_weighted);
        }
        else  // EB Godunov
        {
            AMREX_D_TERM(Array4<Real const> const& fcx = fcent[0]->const_array(mfi);,
                         Array4<Real const> const& fcy = fcent[1]->const_array(mfi);,
                         Array4<Real const> const& fcz = fcent[2]->const_array(mfi););

            AMREX_D_TERM(Array4<Real const> const& apx = areafrac[0]->const_array(mfi);,
                         Array4<Real const> const& apy = areafrac[1]->const_array(mfi);,
                         Array4<Real const> const& apz = areafrac[2]->const_array(mfi););

            Array4<Real const> const& ccent_arr = ccent.const_array(mfi);
            Array4<Real const> const& vfrac_arr = vfrac.const_array(mfi);
            auto const& flags_arr  = flags.const_array(mfi);

            int ngrow = 4;
            FArrayBox tmpfab(amrex::grow(bx,ngrow),  (4*AMREX_SPACEDIM + 2)*ncomp);
            Elixir    eli = tmpfab.elixir();


            if (not known_edgestate)
            {
                AMREX_D_TERM( const auto& u = umac.const_array(mfi);,
                              const auto& v = vmac.const_array(mfi);,
                              const auto& w = wmac.const_array(mfi););

                EBGodunov::ComputeEdgeState( bx, ncomp,
                                             state.array(mfi,state_comp),
                                             AMREX_D_DECL( xed, yed, zed ),
                                             AMREX_D_DECL( u, v, w ),
                                             divu.array(mfi),
                                             fq.array(mfi,fq_comp),
                                             geom, dt, h_bc, d_bc,
                                             iconserv.data(),
                                             tmpfab.dataPtr(),
                                             flags_arr,
                                             AMREX_D_DECL( apx, apy, apz ),
                                             vfrac_arr,
                                             AMREX_D_DECL( fcx, fcy, fcz ),
                                             ccent_arr,
                                             is_velocity );
            }

            HydroUtils::EB_ComputeFluxes( bx,
                                          AMREX_D_DECL( fx, fy, fz ),
                                          AMREX_D_DECL( uc, vc, wc ),
                                          AMREX_D_DECL( xed, yed, zed ),
                                          AMREX_D_DECL( apx, apy, apz ),
                                          geom, ncomp, flags_arr, fluxes_are_area_weighted );

            HydroUtils::EB_ComputeDivergence( bx,
                                              advc_arr,
                                              AMREX_D_DECL( fx, fy, fz ),
                                              vfrac_arr, ncomp, geom, mult, fluxes_are_area_weighted );

	  }
	}
    }

    advc.FillBoundary(geom.periodicity());

    MultiFab* sstate;
    if (redistribution_type == "StateRedist")
    {
      // Create temporary holder for sync "state" passed in via aofs
      // Do this so we're not overwriting the "state" as we go through the redistribution
      // process.
      sstate = new MultiFab(state.boxArray(),state.DistributionMap(),ncomp,state.nGrow(),
			    MFInfo(),ebfact);
      MultiFab::Copy(*sstate,aofs,aofs_comp,0,ncomp,state.nGrow());
    }
    else
    {
      // Doesn't matter what we put here, sstate only gets used for StateRedist
      sstate = &aofs;
    }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(aofs, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        auto const& bx = mfi.tilebox();

        auto const& flagfab   = ebfact.getMultiEBCellFlagFab()[mfi];
        auto const& flags_arr = flagfab.const_array();

        if (flagfab.getType(bx) != FabType::covered )
	{
	  auto const& aofs_arr = aofs.array(mfi, aofs_comp);
	  auto const& advc_arr = advc.array(mfi);

	  // FIXME? not sure if 4 is really needed or if 3 could do
	  // But this is a safe choice
	  if (flagfab.getType(grow(bx,4)) != FabType::regular)
	  {
	    //
	    // Redistribute
	    //
	    AMREX_D_TERM( auto apx = ebfact.getAreaFrac()[0]->const_array(mfi);,
			  auto apy = ebfact.getAreaFrac()[1]->const_array(mfi);,
			  auto apz = ebfact.getAreaFrac()[2]->const_array(mfi); );

	    AMREX_D_TERM( Array4<Real const> fcx = ebfact.getFaceCent()[0]->const_array(mfi);,
			  Array4<Real const> fcy = ebfact.getFaceCent()[1]->const_array(mfi);,
			  Array4<Real const> fcz = ebfact.getFaceCent()[2]->const_array(mfi););

	    Array4<Real const> ccent_arr = ebfact.getCentroid().const_array(mfi);
            Array4<Real const> const& vfrac_arr = vfrac.const_array(mfi);

	    // This is scratch space if calling StateRedistribute,
            //  but is used as the weights (here set to 1) if calling
            //  FluxRedistribute
	    Box gbx = bx;

	    if (redistribution_type == "StateRedist")
	      gbx.grow(3);
	    else if (redistribution_type == "FluxRedist")
	      gbx.grow(2);

	    FArrayBox tmpfab(gbx, ncomp*2);
	    Elixir eli = tmpfab.elixir();
            Array4<Real> scratch = tmpfab.array(0);
            if (redistribution_type == "FluxRedist")
            {
                amrex::ParallelFor(Box(scratch),
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                { scratch(i,j,k) = 1.;});
            }
	    Array4<Real> divtmp_redist_arr = tmpfab.array(ncomp);

	    // Redistribute
	    //
	    // For StateRedistribution, we use the Sync as the "state".
	    // This may lead to oversmoothing.
	    //
            Redistribution::Apply( bx, ncomp, divtmp_redist_arr, advc_arr,
                                   sstate->const_array(mfi, 0), scratch, flags_arr,
                                   AMREX_D_DECL(apx,apy,apz), vfrac_arr,
                                   AMREX_D_DECL(fcx,fcy,fcz), ccent_arr, d_bc,
                                   geom, dt, redistribution_type );

            // Subtract contribution to sync aofs -- sign of divergence is aofs is opposite
            // of sign to div computed by EB_ComputeDivergence, thus it must be subtracted.
            amrex::ParallelFor(bx, ncomp, [aofs_arr, divtmp_redist_arr]
            AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            { aofs_arr( i, j, k, n ) -= divtmp_redist_arr( i, j, k, n ); });
	  }
	  else
	  {
	    // Subtract contribution to sync aofs -- sign of divergence is aofs is opposite
            // of sign to div computed by EB_ComputeDivergence, thus it must be subtracted.
            amrex::ParallelFor(bx, ncomp, [aofs_arr, advc_arr]
            AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            { aofs_arr( i, j, k, n ) -= advc_arr( i, j, k, n ); });
	  }
	}
    }
}
/** @} */

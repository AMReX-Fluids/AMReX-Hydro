/**
 * \file hydro_ebmol.cpp
 * \addtogroup EBMOL
 * @{
 *
 */

#include <hydro_mol.H>
#include <hydro_ebmol.H>
#include <hydro_constants.H>
#include <hydro_utils.H>
// #include <AMReX_MultiCutFab.H>
#include <hydro_redistribution.H>

using namespace amrex;


void
EBMOL::ComputeAofs ( MultiFab& aofs, int aofs_comp, int ncomp,
                     MultiFab const& state, int state_comp,
                     AMREX_D_DECL( MultiFab const& umac,
                                   MultiFab const& vmac,
                                   MultiFab const& wmac),
                     AMREX_D_DECL( MultiFab& xedge,
                                   MultiFab& yedge,
                                   MultiFab& zedge),
                     int  edge_comp,
                     bool known_edgestate,
                     AMREX_D_DECL( MultiFab& xfluxes,
                                   MultiFab& yfluxes,
                                   MultiFab& zfluxes),
                     int fluxes_comp,
                     MultiFab const& divu,
                     Vector<BCRec> const& bcs,
                     BCRec  const* d_bcrec_ptr,
                     Gpu::DeviceVector<int>& iconserv,
                     Geometry const&  geom,
                     const Real dt,
                     const bool is_velocity,
                     std::string redistribution_type )
{
    BL_PROFILE("EBMOL::ComputeAofs()");

    bool fluxes_are_area_weighted = true;

    int const* iconserv_ptr = iconserv.data();

    AMREX_ALWAYS_ASSERT(aofs.nComp()  >= aofs_comp  + ncomp);
    AMREX_ALWAYS_ASSERT(state.nComp() >= state_comp + ncomp);
    AMREX_D_TERM( AMREX_ALWAYS_ASSERT(xedge.nComp() >= edge_comp  + ncomp);,
                  AMREX_ALWAYS_ASSERT(yedge.nComp() >= edge_comp  + ncomp);,
                  AMREX_ALWAYS_ASSERT(zedge.nComp() >= edge_comp  + ncomp););
    AMREX_D_TERM( AMREX_ALWAYS_ASSERT(xfluxes.nComp() >= fluxes_comp  + ncomp);,
                  AMREX_ALWAYS_ASSERT(yfluxes.nComp() >= fluxes_comp  + ncomp);,
                  AMREX_ALWAYS_ASSERT(zfluxes.nComp() >= fluxes_comp  + ncomp););
    AMREX_ALWAYS_ASSERT(aofs.nGrow() == 0);
    AMREX_D_TERM( AMREX_ALWAYS_ASSERT(xfluxes.nGrow() == xedge.nGrow());,
                  AMREX_ALWAYS_ASSERT(yfluxes.nGrow() == yedge.nGrow());,
                  AMREX_ALWAYS_ASSERT(zfluxes.nGrow() == zedge.nGrow()););

    // To compute edge states, need at least 2 ghost cells in state
    if ( !known_edgestate )
        AMREX_ALWAYS_ASSERT(state.nGrow() >= 2);

    // If !known_edgestate, need 2 additional cells in state to compute
    //  the slopes needed to compute the edge state, since MOL uses slope
    //  order==2.
    int halo = known_edgestate ? 0 : 2;

    AMREX_ALWAYS_ASSERT(state.hasEBFabFactory());
    auto const& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(state.Factory());

    // Create temporary holder for advection term. Needed so we can call FillBoundary.
    MultiFab advc(state.boxArray(),state.DistributionMap(),ncomp,3,MFInfo(),ebfactory);
    advc.setVal(0.);

    Box  const& domain = geom.Domain();
    MFItInfo mfi_info;

    if (Gpu::notInLaunchRegion())  mfi_info.EnableTiling().SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(aofs,mfi_info); mfi.isValid(); ++mfi)
    {
        auto const& bx = mfi.tilebox();

	AMREX_D_TERM( const Box& xbx = mfi.nodaltilebox(0);,
                      const Box& ybx = mfi.nodaltilebox(1);,
                      const Box& zbx = mfi.nodaltilebox(2); );

        AMREX_D_TERM( Array4<Real> fx = xfluxes.array(mfi,fluxes_comp);,
                      Array4<Real> fy = yfluxes.array(mfi,fluxes_comp);,
                      Array4<Real> fz = zfluxes.array(mfi,fluxes_comp););

	AMREX_D_TERM( Array4<Real> xed = xedge.array(mfi,edge_comp);,
                      Array4<Real> yed = yedge.array(mfi,edge_comp);,
                      Array4<Real> zed = zedge.array(mfi,edge_comp););

        // Initialize covered cells
        auto const& flagfab = ebfactory.getMultiEBCellFlagFab()[mfi];
        auto const& flag    = flagfab.const_array();

        if (flagfab.getType(bx) == FabType::covered)
        {
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

            AMREX_D_TERM( Array4<Real const> u = umac.const_array(mfi);,
                          Array4<Real const> v = vmac.const_array(mfi);,
                          Array4<Real const> w = wmac.const_array(mfi););

	    Array4<Real> advc_arr = advc.array(mfi);

            bool regular = flagfab.getType(amrex::grow(bx,halo)) == FabType::regular;

            if (!regular)
            {
                AMREX_D_TERM( Array4<Real const> fcx = ebfactory.getFaceCent()[0]->const_array(mfi);,
                              Array4<Real const> fcy = ebfactory.getFaceCent()[1]->const_array(mfi);,
                              Array4<Real const> fcz = ebfactory.getFaceCent()[2]->const_array(mfi););

                AMREX_D_TERM( auto apx = ebfactory.getAreaFrac()[0]->const_array(mfi);,
                              auto apy = ebfactory.getAreaFrac()[1]->const_array(mfi);,
                              auto apz = ebfactory.getAreaFrac()[2]->const_array(mfi); );

                Array4<Real const> ccc = ebfactory.getCentroid().const_array(mfi);
                auto vfrac = ebfactory.getVolFrac().const_array(mfi);

                // Compute edge state if needed
                if (!known_edgestate)
                {
		    Array4<Real const> const q = state.const_array(mfi,state_comp);

                    EBMOL::ComputeEdgeState( bx, AMREX_D_DECL(xed,yed,zed),
                                             q, ncomp,
                                             AMREX_D_DECL(u,v,w),
                                             domain, bcs, d_bcrec_ptr,
                                             AMREX_D_DECL(fcx,fcy,fcz),
                                             ccc, vfrac, flag, is_velocity );
                }

                // Compute fluxes
                HydroUtils::EB_ComputeFluxes(bx, AMREX_D_DECL(fx,fy,fz),
                                             AMREX_D_DECL(u,v,w),
                                             AMREX_D_DECL(xed,yed,zed),
                                             AMREX_D_DECL(apx,apy,apz),
                                             geom, ncomp, flag, fluxes_are_area_weighted );

                //
                // Compute divergence
                //

		// Compute -div because that's what redistribution needs
                Real mult = -1.0;
                HydroUtils::EB_ComputeDivergence(bx, advc_arr,
                                                 AMREX_D_DECL(fx,fy,fz),
                                                 vfrac, ncomp, geom,
                                                 mult, fluxes_are_area_weighted );
                // Account for extra term needed for convective differencing
		// Don't forget we're mutliplying by -1.0 here...
                auto const& q = state.array(mfi, state_comp);
                auto const& divu_arr  = divu.array(mfi);
		amrex::ParallelFor(bx, ncomp, [=]
                AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (!iconserv_ptr[n])
                        advc_arr( i, j, k, n ) += q(i,j,k,n)*divu_arr(i,j,k);

                });
            }
            else
            {
                // Compute edge state if needed
                if (!known_edgestate)
                {
                    Array4<Real const> const q = state.const_array(mfi,state_comp);
                    MOL::ComputeEdgeState( bx,
                                           AMREX_D_DECL( xed, yed, zed ),
                                           q, ncomp,
                                           AMREX_D_DECL( u, v, w ),
                                           domain, bcs, d_bcrec_ptr,
                                           is_velocity);
                }

                // Compute fluxes
                HydroUtils::ComputeFluxes( bx,
                                           AMREX_D_DECL(fx,fy,fz),
                                           AMREX_D_DECL(u,v,w),
                                           AMREX_D_DECL(xed,yed,zed),
                                           geom, ncomp, fluxes_are_area_weighted  );

                // Compute divergence
                // We use minus sign, i.e. -div, for consistency with EB above
                Real mult = -1.0;
                HydroUtils::ComputeDivergence(bx, advc_arr,
                                              AMREX_D_DECL(fx,fy,fz),
                                              ncomp, geom,
                                              mult, fluxes_are_area_weighted );

                // Account for extra term needed for convective differencing
                auto const& q = state.array(mfi, state_comp);
                auto const& divu_arr  = divu.array(mfi);
                amrex::ParallelFor(bx, ncomp, [=]
                AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (!iconserv_ptr[n])
		      advc_arr( i, j, k, n ) += q(i,j,k,n)*divu_arr(i,j,k);
                });
	    }
        }
    }

    advc.FillBoundary(geom.periodicity());

    if (Gpu::notInLaunchRegion())  mfi_info.EnableTiling().SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(aofs,mfi_info); mfi.isValid(); ++mfi)
    {
        auto const& bx = mfi.tilebox();

        auto const& flagfab = ebfactory.getMultiEBCellFlagFab()[mfi];
        auto const& flag    = flagfab.const_array();
	auto const& aofs_arr = aofs.array(mfi, aofs_comp);

        if (flagfab.getType(bx) != FabType::covered )
	{
	  // FIXME? not sure if 4 is really needed or if 3 could do
	  // But this is a safe choice
	  if (flagfab.getType(grow(bx,4)) != FabType::regular)
	  {
	    //
	    // Redistribute
	    //
	    AMREX_D_TERM( auto apx = ebfactory.getAreaFrac()[0]->const_array(mfi);,
			  auto apy = ebfactory.getAreaFrac()[1]->const_array(mfi);,
			  auto apz = ebfactory.getAreaFrac()[2]->const_array(mfi); );

	    AMREX_D_TERM( Array4<Real const> fcx = ebfactory.getFaceCent()[0]->const_array(mfi);,
			  Array4<Real const> fcy = ebfactory.getFaceCent()[1]->const_array(mfi);,
			  Array4<Real const> fcz = ebfactory.getFaceCent()[2]->const_array(mfi););

	    Array4<Real const> ccc = ebfactory.getCentroid().const_array(mfi);
	    auto vfrac = ebfactory.getVolFrac().const_array(mfi);

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
				   AMREX_D_DECL(apx,apy,apz), vfrac,
				   AMREX_D_DECL(fcx,fcy,fcz), ccc, d_bcrec_ptr,
				   geom, dt, redistribution_type );

	    // Change sign because we computed -div for all cases
	    amrex::ParallelFor(bx, ncomp, [aofs_arr]
	    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	    { aofs_arr( i, j, k, n ) *=  - 1.0; });
	  }
	  else
	  {
	    // Change sign because we computed -div for all cases
	    auto const& advc_arr = advc.array(mfi);
	    amrex::ParallelFor(bx, ncomp, [aofs_arr, advc_arr]
	    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	    { aofs_arr( i, j, k, n ) =  - advc_arr(i,j,k,n); });
	  }
	}
    }
}

void
EBMOL::ComputeSyncAofs ( MultiFab& aofs, int aofs_comp, int ncomp,
                         MultiFab const& state, int state_comp,
                         AMREX_D_DECL( MultiFab const& umac,
                                       MultiFab const& vmac,
                                       MultiFab const& wmac),
                         AMREX_D_DECL( MultiFab const& ucorr,
                                       MultiFab const& vcorr,
                                       MultiFab const& wcorr),
                         AMREX_D_DECL( MultiFab& xedge,
                                       MultiFab& yedge,
                                       MultiFab& zedge),
                         int  edge_comp,
                         bool known_edgestate,
                         AMREX_D_DECL( MultiFab& xfluxes,
                                       MultiFab& yfluxes,
                                       MultiFab& zfluxes),
                         int fluxes_comp,
                         Vector<BCRec> const& bcs,
                         BCRec  const* d_bcrec_ptr,
                         Geometry const&  geom,
                         const Real dt,
                         const bool is_velocity,
                         std::string redistribution_type )
{
    BL_PROFILE("EBMOL::ComputeSyncAofs()");

    bool fluxes_are_area_weighted = true;

    AMREX_ALWAYS_ASSERT(state.nComp() >= state_comp + ncomp);
    AMREX_ALWAYS_ASSERT(aofs.nComp()  >= aofs_comp  + ncomp);
    AMREX_D_TERM( AMREX_ALWAYS_ASSERT(xedge.nComp() >= edge_comp  + ncomp);,
                  AMREX_ALWAYS_ASSERT(yedge.nComp() >= edge_comp  + ncomp);,
                  AMREX_ALWAYS_ASSERT(zedge.nComp() >= edge_comp  + ncomp););
    AMREX_D_TERM( AMREX_ALWAYS_ASSERT(xfluxes.nComp() >= fluxes_comp  + ncomp);,
                  AMREX_ALWAYS_ASSERT(yfluxes.nComp() >= fluxes_comp  + ncomp);,
                  AMREX_ALWAYS_ASSERT(zfluxes.nComp() >= fluxes_comp  + ncomp););
    AMREX_D_TERM( AMREX_ALWAYS_ASSERT(xfluxes.nGrow() == xedge.nGrow());,
                  AMREX_ALWAYS_ASSERT(yfluxes.nGrow() == yedge.nGrow());,
                  AMREX_ALWAYS_ASSERT(zfluxes.nGrow() == zedge.nGrow()););


    // To compute edge states, need at least 2 ghost cells in state
    if ( !known_edgestate )
        AMREX_ALWAYS_ASSERT(state.nGrow() >= 2);

    AMREX_ALWAYS_ASSERT(state.hasEBFabFactory());
    auto const& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(state.Factory());

    // Need 2 grow cells in state to compute the slopes needed to compute the edge state.
    int halo = known_edgestate ? 0 : 2;

    // Create temporary holder for advection term. Needed to fill ghost cells.
    MultiFab advc(state.boxArray(),state.DistributionMap(),ncomp,3,MFInfo(),ebfactory);
    advc.setVal(0.);

    Box  const& domain = geom.Domain();
    MFItInfo mfi_info;

    if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(aofs,mfi_info); mfi.isValid(); ++mfi)
    {
        auto const& bx = mfi.tilebox();
	AMREX_D_TERM( const Box& xbx = mfi.nodaltilebox(0);,
                      const Box& ybx = mfi.nodaltilebox(1);,
                      const Box& zbx = mfi.nodaltilebox(2); );

        AMREX_D_TERM( Array4<Real> fx = xfluxes.array(mfi,fluxes_comp);,
                      Array4<Real> fy = yfluxes.array(mfi,fluxes_comp);,
                      Array4<Real> fz = zfluxes.array(mfi,fluxes_comp););

	AMREX_D_TERM( Array4<Real> xed = xedge.array(mfi,edge_comp);,
                      Array4<Real> yed = yedge.array(mfi,edge_comp);,
                      Array4<Real> zed = zedge.array(mfi,edge_comp););

        // Initialize covered cells
        auto const& flagfab = ebfactory.getMultiEBCellFlagFab()[mfi];
        auto const& flag    = flagfab.const_array();

        if (flagfab.getType(bx) == FabType::covered)
        {
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
            AMREX_D_TERM( Array4<Real const> uc = ucorr.const_array(mfi);,
                          Array4<Real const> vc = vcorr.const_array(mfi);,
                          Array4<Real const> wc = wcorr.const_array(mfi););

	    Array4<Real> advc_arr = advc.array(mfi);

            bool regular = flagfab.getType(amrex::grow(bx,halo)) == FabType::regular;

            if (!regular)
            {
                AMREX_D_TERM( Array4<Real const> fcx = ebfactory.getFaceCent()[0]->const_array(mfi);,
                              Array4<Real const> fcy = ebfactory.getFaceCent()[1]->const_array(mfi);,
                              Array4<Real const> fcz = ebfactory.getFaceCent()[2]->const_array(mfi););

                Array4<Real const> ccc = ebfactory.getCentroid().const_array(mfi);

                auto vfrac = ebfactory.getVolFrac().const_array(mfi);

                AMREX_D_TERM( auto apx = ebfactory.getAreaFrac()[0]->const_array(mfi);,
                              auto apy = ebfactory.getAreaFrac()[1]->const_array(mfi);,
                              auto apz = ebfactory.getAreaFrac()[2]->const_array(mfi); );

                // Compute edge state if needed
                if (!known_edgestate)
                {
		    Array4<Real const> const q = state.const_array(mfi,state_comp);

		    AMREX_D_TERM( Array4<Real const> u = umac.const_array(mfi);,
                                  Array4<Real const> v = vmac.const_array(mfi);,
                                  Array4<Real const> w = wmac.const_array(mfi););

                    EBMOL::ComputeEdgeState( bx,
                                             AMREX_D_DECL(xed,yed,zed),
                                             q, ncomp,
                                             AMREX_D_DECL(u,v,w),
                                             domain, bcs, d_bcrec_ptr,
                                             AMREX_D_DECL(fcx,fcy,fcz),
                                             ccc, vfrac, flag,
                                             is_velocity );
                }

                // Compute fluxes
                HydroUtils::EB_ComputeFluxes(bx,
                                             AMREX_D_DECL(fx,fy,fz),
                                             AMREX_D_DECL(uc,vc,wc),
                                             AMREX_D_DECL(xed,yed,zed),
                                             AMREX_D_DECL(apx,apy,apz),
                                             geom, ncomp, flag, fluxes_are_area_weighted  );

                // Compute divergence
		// Compute -div because that's what redistribution needs
                Real mult = -1.0;
                HydroUtils::EB_ComputeDivergence(bx, advc_arr,
                                                 AMREX_D_DECL(fx,fy,fz), vfrac,
                                                 ncomp, geom, mult, fluxes_are_area_weighted  );
            }
            else
            {
                // Compute edge state if needed
                if (!known_edgestate)
                {
		    Array4<Real const> const q = state.const_array(mfi,state_comp);

		    AMREX_D_TERM( Array4<Real const> u = umac.const_array(mfi);,
                                  Array4<Real const> v = vmac.const_array(mfi);,
                                  Array4<Real const> w = wmac.const_array(mfi););

                    MOL::ComputeEdgeState( bx,
                                           AMREX_D_DECL( xed, yed, zed ),
                                           q, ncomp,
                                           AMREX_D_DECL( u, v, w ),
                                           domain, bcs, d_bcrec_ptr,
                                           is_velocity);

                }

                // Compute fluxes
                HydroUtils::ComputeFluxes(bx,
                                          AMREX_D_DECL(fx,fy,fz),
                                          AMREX_D_DECL(uc,vc,wc),
                                          AMREX_D_DECL(xed,yed,zed),
                                          geom, ncomp, fluxes_are_area_weighted  );

                // Compute divergence
                // We use minus sign, i.e. -div, for consistency with EB above
                Real mult = -1.0;
                HydroUtils::ComputeDivergence(bx, advc_arr,
                                              AMREX_D_DECL(fx,fy,fz),
                                              ncomp, geom,
                                              mult, fluxes_are_area_weighted );
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
			    MFInfo(),ebfactory);
      MultiFab::Copy(*sstate,aofs,aofs_comp,0,ncomp,state.nGrow());
    }
    else
    {
      // Doesn't matter what we put here, sstate only gets used for StateRedist
      sstate = &aofs;
    }

    if (Gpu::notInLaunchRegion())  mfi_info.EnableTiling().SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(aofs,mfi_info); mfi.isValid(); ++mfi)
    {
        auto const& bx = mfi.tilebox();

        auto const& flagfab = ebfactory.getMultiEBCellFlagFab()[mfi];
        auto const& flag    = flagfab.const_array();

        if (flagfab.getType(bx) != FabType::covered)
	{
	  // FIXME? not sure if 4 is really needed or if 3 could do
	  // But this is a safe choice
	  if (flagfab.getType(grow(bx,4)) != FabType::regular)
	  {
	    //
	    // Redistribute
	    //
	    AMREX_D_TERM( auto apx = ebfactory.getAreaFrac()[0]->const_array(mfi);,
			  auto apy = ebfactory.getAreaFrac()[1]->const_array(mfi);,
			  auto apz = ebfactory.getAreaFrac()[2]->const_array(mfi); );

	    AMREX_D_TERM( Array4<Real const> fcx = ebfactory.getFaceCent()[0]->const_array(mfi);,
			  Array4<Real const> fcy = ebfactory.getFaceCent()[1]->const_array(mfi);,
			  Array4<Real const> fcz = ebfactory.getFaceCent()[2]->const_array(mfi););

	    Array4<Real const> ccc = ebfactory.getCentroid().const_array(mfi);
	    auto vfrac = ebfactory.getVolFrac().const_array(mfi);

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
	    Redistribution::Apply( bx, ncomp,  divtmp_redist_arr, advc.array(mfi),
				   sstate->const_array(mfi, 0), scratch, flag,
				   AMREX_D_DECL(apx,apy,apz), vfrac,
				   AMREX_D_DECL(fcx,fcy,fcz), ccc, d_bcrec_ptr,
				   geom, dt, redistribution_type );

	    // Subtract contribution to sync aofs -- sign of divergence in aofs is opposite
	    // of sign of div as computed by EB_ComputeDivergence, thus it must be subtracted.
	    auto const& aofs_arr = aofs.array(mfi, aofs_comp);

	    amrex::ParallelFor(bx, ncomp, [aofs_arr, divtmp_redist_arr]
	    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            { aofs_arr( i, j, k, n ) -=  divtmp_redist_arr( i, j, k, n ); });
	  }
	  else
	  {
	    // Subtract contribution to sync aofs -- sign of divergence in aofs is opposite
	    // of sign of div as computed here, and thus it must be subtracted.
	    auto const& aofs_arr = aofs.array(mfi, aofs_comp);
	    Array4<Real> advc_arr = advc.array(mfi);

	    amrex::ParallelFor(bx, ncomp, [aofs_arr, advc_arr]
            AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            { aofs_arr( i, j, k, n ) -= advc_arr( i, j, k, n ); });
	  }
	}
    }
}
/** @} */

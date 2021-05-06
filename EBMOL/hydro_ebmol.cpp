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

    // To compute edge states, need at least 2 more ghost cells in state than in
    //  xedge
    if ( !known_edgestate )
        AMREX_ALWAYS_ASSERT(state.nGrow() >= xedge.nGrow()+2);

    AMREX_ALWAYS_ASSERT(state.hasEBFabFactory());
    auto const& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(state.Factory());

    Box  const& domain = geom.Domain();
    MFItInfo mfi_info;

    if (Gpu::notInLaunchRegion())  mfi_info.EnableTiling().SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(aofs,mfi_info); mfi.isValid(); ++mfi)
    {
        auto const& bx = mfi.tilebox();
	int ng_f = xfluxes.nGrow();
	AMREX_D_TERM( const Box& xbx = mfi.grownnodaltilebox(0,ng_f);,
                      const Box& ybx = mfi.grownnodaltilebox(1,ng_f);,
                      const Box& zbx = mfi.grownnodaltilebox(2,ng_f); );

        AMREX_D_TERM( Array4<Real> fx = xfluxes.array(mfi,fluxes_comp);,
                      Array4<Real> fy = yfluxes.array(mfi,fluxes_comp);,
                      Array4<Real> fz = zfluxes.array(mfi,fluxes_comp););

	AMREX_D_TERM( Array4<Real> xed = xedge.array(mfi,edge_comp);,
                      Array4<Real> yed = yedge.array(mfi,edge_comp);,
                      Array4<Real> zed = zedge.array(mfi,edge_comp););

        // Initialize covered cells
        auto const& flagfab = ebfactory.getMultiEBCellFlagFab()[mfi];
        auto const& flag    = flagfab.const_array();
	auto const& gtbx    = mfi.growntilebox(ng_f);

        if (flagfab.getType(gtbx) == FabType::covered)
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

	    // Grown box on which to compute the edge states and fluxes for regular boxes
	    Box gbx = mfi.growntilebox(ng_f);

	    Box tmpbox = amrex::surroundingNodes(gbx);
	    // Space for fluxes
	    int tmpcomp = ncomp*AMREX_SPACEDIM;

	    // PeleLM needs valid flux on all ng_f cells. If !known_edgestate, need 2 additional
	    //  cells in state to compute the slopes needed to compute the edge state.
	    int halo = known_edgestate ? 0 : ng_f+2;
            bool regular = flagfab.getType(amrex::grow(bx,halo)) == FabType::regular;
	    if (!regular)
            {
                // Grown box on which to compute the edge states and fluxes for EB containing boxes
                // need at least 2 filled ghost cells all around for redistribution
                gbx = amrex::grow(bx,ng_f);
                tmpbox = amrex::surroundingNodes(gbx);
                // Not sure if we really need 3(incflo) here or 2
                int ng_diff = 3-ng_f;
                if ( ng_diff>0 )
                    tmpbox.grow(ng_diff);

                // Add space for the temporaries needed by Redistribute
#if (AMREX_SPACEDIM == 3)
                tmpcomp += ncomp;
#else
                tmpcomp += 2*ncomp;
#endif
	    }

	    FArrayBox tmpfab(tmpbox, tmpcomp);
	    Elixir eli = tmpfab.elixir();

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

                    EBMOL::ComputeEdgeState( gbx, AMREX_D_DECL(xed,yed,zed),
                                             q, ncomp,
                                             AMREX_D_DECL(u,v,w),
                                             domain, bcs, d_bcrec_ptr,
                                             AMREX_D_DECL(fcx,fcy,fcz),
                                             ccc, vfrac, flag, is_velocity );
                }

                // Compute fluxes
                HydroUtils::EB_ComputeFluxes(gbx, AMREX_D_DECL(fx,fy,fz),
                                             AMREX_D_DECL(u,v,w),
                                             AMREX_D_DECL(xed,yed,zed),
                                             AMREX_D_DECL(apx,apy,apz),
                                             geom, ncomp, flag, fluxes_are_area_weighted );

                //
                // Compute divergence and redistribute
                //
		// div at ncomp*3 to make space for the 3 redistribute temporaries
                Array4<Real> divtmp_arr = tmpfab.array(ncomp*3);

                // Compute conservative divergence
                // Redistribute needs 2 ghost cells in div
	        Box g2bx = amrex::grow(bx,2);
                Real mult = -1.0;
                HydroUtils::EB_ComputeDivergence(g2bx, divtmp_arr,
                                                 AMREX_D_DECL(fx,fy,fz),
                                                 vfrac, ncomp, geom,
                                                 mult, fluxes_are_area_weighted );

                // Account for extra term needed for convective differencing
                auto const& q = state.array(mfi, state_comp);
                auto const& divu_arr  = divu.array(mfi);
                amrex::ParallelFor(g2bx, ncomp, [=]
                AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (!iconserv_ptr[n])
                        divtmp_arr( i, j, k, n ) -= q(i,j,k,n)*divu_arr(i,j,k);
                });

                // Redistribute
		Array4<Real> scratch = tmpfab.array(0);
                auto const& aofs_arr = aofs.array(mfi, aofs_comp);
                Redistribution::Apply( bx, ncomp, aofs_arr, divtmp_arr,
                                       state.const_array(mfi, state_comp), scratch, flag,
                                       AMREX_D_DECL(apx,apy,apz), vfrac,
                                       AMREX_D_DECL(fcx,fcy,fcz), ccc, geom, dt,
                                       redistribution_type );

                // Change sign because for EB redistribution we compute -div
                amrex::ParallelFor(bx, ncomp, [aofs_arr]
                AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                { aofs_arr( i, j, k, n ) *=  - 1.0; });
            }
            else
            {
                // Compute edge state if needed
                if (!known_edgestate)
                {
                    Array4<Real const> const q = state.const_array(mfi,state_comp);
                    MOL::ComputeEdgeState( gbx,
                                           AMREX_D_DECL( xed, yed, zed ),
                                           q, ncomp,
                                           AMREX_D_DECL( u, v, w ),
                                           domain, bcs, d_bcrec_ptr,
                                           is_velocity);
                }

                // Compute fluxes
                HydroUtils::ComputeFluxes( gbx,
                                           AMREX_D_DECL(fx,fy,fz),
                                           AMREX_D_DECL(u,v,w),
                                           AMREX_D_DECL(xed,yed,zed),
                                           geom, ncomp, fluxes_are_area_weighted  );

                // Compute divergence -- always use conservative form
                // If convetive form is required, the next parallel for
                // will take care of it.
                Real mult = 1.0;
                std::vector<int>  div_iconserv(ncomp,1);
                HydroUtils::ComputeDivergence(bx, aofs.array(mfi, aofs_comp),
                                              AMREX_D_DECL(fx,fy,fz),
                                              AMREX_D_DECL( xed, yed, zed ),
                                              AMREX_D_DECL( u, v, w ),
                                              ncomp, geom, div_iconserv.data(),
                                              mult, fluxes_are_area_weighted );

                // Account for extra term needed for convective differencing
                auto const& aofs_arr  = aofs.array(mfi, aofs_comp);
                auto const& q = state.array(mfi, state_comp);
                auto const& divu_arr  = divu.array(mfi);
                amrex::ParallelFor(bx, ncomp, [=]
                AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (!iconserv_ptr[n])
                        aofs_arr( i, j, k, n ) -= q(i,j,k,n)*divu_arr(i,j,k);
                });

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


    // To compute edge states, need at least 2 more ghost cells in state than in
    //  xedge
    if ( !known_edgestate )
        AMREX_ALWAYS_ASSERT(state.nGrow() >= xedge.nGrow()+2);

    // Only conservative scheme for EB
    std::vector<int> iconserv(ncomp,1);

    AMREX_ALWAYS_ASSERT(state.hasEBFabFactory());
    auto const& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(state.Factory());

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

	    Box tmpbox = amrex::surroundingNodes(bx);
	    int tmpcomp = ncomp*(AMREX_SPACEDIM+1);

	    Box gbx = bx;
	    // Need 2 grow cells in state to compute the slopes needed to compute the edge state.
	    int halo = known_edgestate ? 0 : 2;
            bool regular = flagfab.getType(amrex::grow(bx,halo)) == FabType::regular;
	    if (!regular)
            {
                // Grown box on which to compute the fluxes and divergence.
                gbx.grow(2);
                tmpbox.grow(3);

                // Add space for the temporaries needed by Redistribute
#if (AMREX_SPACEDIM == 3)
                tmpcomp += ncomp;
#else
                tmpcomp += 2*ncomp;
#endif
	    }

	    FArrayBox tmpfab(tmpbox, tmpcomp);
	    Elixir eli = tmpfab.elixir();

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

                    EBMOL::ComputeEdgeState( gbx,
                                             AMREX_D_DECL(xed,yed,zed),
                                             q, ncomp,
                                             AMREX_D_DECL(u,v,w),
                                             domain, bcs, d_bcrec_ptr,
                                             AMREX_D_DECL(fcx,fcy,fcz),
                                             ccc, vfrac, flag,
                                             is_velocity );
                }

                // Compute fluxes
                HydroUtils::EB_ComputeFluxes(gbx,
                                             AMREX_D_DECL(fx,fy,fz),
                                             AMREX_D_DECL(uc,vc,wc),
                                             AMREX_D_DECL(xed,yed,zed),
                                             AMREX_D_DECL(apx,apy,apz),
                                             geom, ncomp, flag, fluxes_are_area_weighted  );

                //
                // Compute divergence and redistribute
                //
		// div at ncomp*3 to make space for the 3 redistribute temporaries
                Array4<Real> divtmp_arr = tmpfab.array(ncomp*3);
                Array4<Real> divtmp_redist_arr = tmpfab.array(ncomp*4);

                // Compute conservative divergence
                Box g2bx = amrex::grow(bx,2);
                Real mult = -1.0;
                HydroUtils::EB_ComputeDivergence(g2bx, divtmp_arr,
                                                 AMREX_D_DECL(fx,fy,fz), vfrac,
                                                 ncomp, geom, mult, fluxes_are_area_weighted  );

                // Redistribute
		Array4<Real> scratch = tmpfab.array(0);
                Redistribution::Apply( bx, ncomp,  divtmp_redist_arr, divtmp_arr,
                                       state.const_array(mfi, state_comp), scratch, flag,
                                       AMREX_D_DECL(apx,apy,apz), vfrac,
                                       AMREX_D_DECL(fcx,fcy,fcz), ccc, geom, dt,
                                       redistribution_type );

                // Subtract contribution to sync aofs -- sign of divergence in aofs is opposite
                // of sign of div as computed by EB_ComputeDivergence, thus it must be subtracted.
                auto const& aofs_arr = aofs.array(mfi, aofs_comp);

                amrex::ParallelFor(bx, ncomp, [aofs_arr, divtmp_redist_arr]
                AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                { aofs_arr( i, j, k, n ) +=  -divtmp_redist_arr( i, j, k, n ); });

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
                Array4<Real> divtmp_arr = tmpfab.array(ncomp*AMREX_SPACEDIM);
                Real mult = 1.0;
                HydroUtils::ComputeDivergence(bx, divtmp_arr,
                                              AMREX_D_DECL(fx,fy,fz),
                                              AMREX_D_DECL( xed, yed, zed ),
                                              AMREX_D_DECL( uc, vc, wc ),
                                              ncomp, geom, iconserv.data(),
                                              mult, fluxes_are_area_weighted );

                // Sum contribution to sync aofs
                auto const& aofs_arr = aofs.array(mfi, aofs_comp);

                amrex::ParallelFor(bx, ncomp, [aofs_arr, divtmp_arr]
                AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                { aofs_arr( i, j, k, n ) += divtmp_arr( i, j, k, n ); });
            }
        }
    }
}

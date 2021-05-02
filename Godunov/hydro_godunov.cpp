#include <hydro_godunov.H>
#include <hydro_utils.H>

using namespace amrex;


void
Godunov::ComputeAofs ( MultiFab& aofs, const int aofs_comp, const int ncomp,
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
                       BCRec const* d_bc,
                       Geometry const& geom,
                       Gpu::DeviceVector<int>& iconserv,
                       const Real dt,
                       const bool use_ppm,
                       const bool use_forces_in_trans,
                       const bool is_velocity  )
{
    BL_PROFILE("Godunov::ComputeAofs()");

#if (AMREX_SPACEDIM==2)
    MultiFab* volume;
    MultiFab* area[AMREX_SPACEDIM];

    if ( geom.IsRZ() )
    {
        const DistributionMapping& dmap = aofs.DistributionMap();
        const BoxArray& grids = aofs.boxArray();
        const int ngrow_vol = aofs.nGrow();

        volume = new MultiFab(grids,dmap,1,ngrow_vol);
        geom.GetVolume(*volume);

        const int ngrow_area = xfluxes.nGrow();

        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
        {
            BoxArray edge_ba(grids);
            area[dir] = new MultiFab(edge_ba.surroundingNodes(dir),dmap,1,ngrow_area);
            geom.GetFaceArea(*area[dir],dir);
        }
    }
#endif

    //FIXME - check on adding tiling here
    for (MFIter mfi(aofs); mfi.isValid(); ++mfi)
    {

        const Box& bx   = mfi.tilebox();

        //
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

        if (not known_edgestate)
        {
            ComputeEdgeState( bx, ncomp,
                              state.array(mfi,state_comp),
                              AMREX_D_DECL( xed, yed, zed ),
                              AMREX_D_DECL( u, v, w ),
                              divu.array(mfi),
                              fq.array(mfi,fq_comp),
                              geom, dt, d_bc,
                              iconserv.data(),
                              use_ppm,
                              use_forces_in_trans,
                              is_velocity );
        }

        Real mult = 1.;
#if (AMREX_SPACEDIM == 2)
	if ( geom.IsRZ() )
	{
            const auto& areax = area[0]->array(mfi);
            const auto& areay = area[1]->array(mfi);
            const auto& vol   = volume->array(mfi);

            HydroUtils::ComputeFluxesRZ( bx,
                                         AMREX_D_DECL( fx, fy, fz ),
                                         AMREX_D_DECL( u, v, w ),
                                         AMREX_D_DECL( xed, yed, zed ),
                                         areax, areay,
                                         ncomp );


            HydroUtils::ComputeDivergenceRZ( bx,
                                             aofs.array(mfi,aofs_comp),
                                             AMREX_D_DECL( fx, fy, fz ),
                                             AMREX_D_DECL( xed, yed, zed ),
                                             AMREX_D_DECL( u, v, w ),
                                             areax, areay, vol,
                                             ncomp, iconserv.data(),
                                             mult);

	}
	else
#endif
	{
            HydroUtils::ComputeFluxes( bx,
                                       AMREX_D_DECL( fx, fy, fz ),
                                       AMREX_D_DECL( u, v, w ),
                                       AMREX_D_DECL( xed, yed, zed ),
                                       geom, ncomp );

            HydroUtils::ComputeDivergence( bx,
                                           aofs.array(mfi,aofs_comp),
                                           AMREX_D_DECL( fx, fy, fz ),
                                           AMREX_D_DECL( xed, yed, zed ),
                                           AMREX_D_DECL( u, v, w ),
                                           ncomp, geom, iconserv.data(),
                                           mult);
	}

	//
	// NOTE this sync cannot protect temporaries in ComputeEdgeState, ComputeFluxes
	// or ComputeDivergence, since functions have their own scope. As soon as the
	// CPU hits the end of the function, it will call the destructor for all
	// temporaries created in that function.
	//
        Gpu::streamSynchronize();  // otherwise we might be using too much memory
    }

}



void
Godunov::ComputeSyncAofs ( MultiFab& aofs, const int aofs_comp, const int ncomp,
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
                           BCRec const* d_bc,
                           Geometry const& geom,
                           Gpu::DeviceVector<int>& iconserv,
                           const Real dt,
                           const bool use_ppm,
                           const bool use_forces_in_trans,
                           const bool is_velocity  )
{
    BL_PROFILE("Godunov::ComputeSyncAofs()");

#if (AMREX_SPACEDIM==2)
    MultiFab* volume;
    MultiFab* area[AMREX_SPACEDIM];

    if ( geom.IsRZ() )
    {
        const DistributionMapping& dmap = aofs.DistributionMap();
        const BoxArray& grids = aofs.boxArray();
        const int ngrow_vol = aofs.nGrow();

        volume = new MultiFab(grids,dmap,1,ngrow_vol);
        geom.GetVolume(*volume);

        const int ngrow_area = xfluxes.nGrow();

        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
        {
            BoxArray edge_ba(grids);
            area[dir] = new MultiFab(edge_ba.surroundingNodes(dir),dmap,1,ngrow_area);
            geom.GetFaceArea(*area[dir],dir);
        }
    }
#endif

    // Sync divergence computation is always conservative
    std::vector<int> div_iconserv(ncomp,1);

    //FIXME - check on adding tiling here
    for (MFIter mfi(aofs); mfi.isValid(); ++mfi)
    {

        const Box& bx   = mfi.tilebox();

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

        if (not known_edgestate)
        {

            AMREX_D_TERM( const auto& u = umac.const_array(mfi);,
                          const auto& v = vmac.const_array(mfi);,
                          const auto& w = wmac.const_array(mfi););

            ComputeEdgeState( bx, ncomp,
                              state.array(mfi,state_comp),
                              AMREX_D_DECL( xed, yed, zed ),
                              AMREX_D_DECL( u, v, w ),
                              divu.array(mfi),
                              fq.array(mfi,fq_comp),
                              geom, dt, d_bc,
                              iconserv.data(),
                              use_ppm,
                              use_forces_in_trans,
                              is_velocity );
        }

        // Temporary divergence
        Box tmpbox = amrex::surroundingNodes(bx);
        int tmpcomp = ncomp*AMREX_SPACEDIM;
        FArrayBox tmpfab(tmpbox, tmpcomp);
        Elixir eli = tmpfab.elixir();
        Array4<Real> divtmp_arr = tmpfab.array();

        Real mult = 1.0; 

#if (AMREX_SPACEDIM == 2)
	if ( geom.IsRZ() )
	{
            const auto& areax = area[0]->array(mfi);
            const auto& areay = area[1]->array(mfi);
            const auto& vol   = volume->array(mfi);

            HydroUtils::ComputeFluxesRZ( bx,
                                         AMREX_D_DECL( fx, fy, fz ),
                                         AMREX_D_DECL( uc, vc, wc ),
                                         AMREX_D_DECL( xed, yed, zed ),
                                         areax, areay,
                                         ncomp );

            HydroUtils::ComputeDivergenceRZ( bx, divtmp_arr,
                                             AMREX_D_DECL( fx, fy, fz ),
                                             AMREX_D_DECL( xed, yed, zed ),
                                             AMREX_D_DECL( uc, vc, wc ),
                                             areax, areay, vol,
                                             ncomp, div_iconserv.data(),
                                             mult);

	}
	else
#endif
	{
            HydroUtils::ComputeFluxes( bx,
                                       AMREX_D_DECL( fx, fy, fz ),
                                       AMREX_D_DECL( uc, vc, wc ),
                                       AMREX_D_DECL( xed, yed, zed ),
                                       geom, ncomp );

            HydroUtils::ComputeDivergence( bx, divtmp_arr,
                                           AMREX_D_DECL( fx, fy, fz ),
                                           AMREX_D_DECL( xed, yed, zed ),
                                           AMREX_D_DECL( uc, vc, wc ),
                                           ncomp, geom, div_iconserv.data(),
                                           mult);
	}

        // Sum contribution to sync aofs
        auto const& aofs_arr = aofs.array(mfi, aofs_comp);

        amrex::ParallelFor(bx, ncomp, [aofs_arr, divtmp_arr]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        { aofs_arr( i, j, k, n ) += divtmp_arr( i, j, k, n ); });

        Gpu::streamSynchronize();  // otherwise we might be using too much memory
    }

}

#include <hydro_mol.H>
#include <hydro_constants.H>
#include <hydro_utils.H>
#include <AMReX_MultiFab.H>

using namespace amrex;


void
MOL::ComputeAofs ( MultiFab& aofs, int aofs_comp, int ncomp,
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
                   Vector<BCRec> const& bcs,
                   BCRec  const* d_bcrec_ptr,
                   int const* iconserv,
                   Geometry const&  geom,
                   Real dt )
{
    BL_PROFILE("MOL::ComputeAofs()");

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

        AMREX_D_TERM( Array4<Real> fx = xfluxes.array(mfi,fluxes_comp);,
                      Array4<Real> fy = yfluxes.array(mfi,fluxes_comp);,
                      Array4<Real> fz = zfluxes.array(mfi,fluxes_comp););

	AMREX_D_TERM( Array4<Real> xed = xedge.array(mfi,edge_comp);,
                      Array4<Real> yed = yedge.array(mfi,edge_comp);,
                      Array4<Real> zed = zedge.array(mfi,edge_comp););

        AMREX_D_TERM( Array4<Real const> u = umac.const_array(mfi);,
                      Array4<Real const> v = vmac.const_array(mfi);,
                      Array4<Real const> w = wmac.const_array(mfi););

        // Grown box on which to compute the edge states and fluxes
        Box gbx = mfi.growntilebox(ng_f);

        // Compute edge state if needed
        if (!known_edgestate)
        {
            Array4<Real const> const q = state.const_array(mfi,state_comp);
            ComputeEdgeState( gbx, AMREX_D_DECL( xed, yed, zed ), q, ncomp,
                              AMREX_D_DECL( u, v, w ), domain, bcs, d_bcrec_ptr);

        }

        // Compute fluxes
        HydroUtils::ComputeFluxes(gbx, AMREX_D_DECL(fx,fy,fz),
                                  AMREX_D_DECL(u,v,w), AMREX_D_DECL(xed,yed,zed),
                                  geom, ncomp );

        // Compute divergence
        HydroUtils::ComputeDivergence(bx, aofs.array(mfi, aofs_comp), AMREX_D_DECL(fx,fy,fz),
                                      AMREX_D_DECL( xed, yed, zed ), AMREX_D_DECL( u, v, w ),
                                      ncomp, geom, iconserv);

    }
}


void
MOL::ComputeSyncAofs ( MultiFab& aofs, int aofs_comp, int ncomp,
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
                       Real dt )

{
    BL_PROFILE("MOL::ComputeSyncAofs()");

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

    // Sync divergence computation is always conservative
    std::vector<int>  iconserv(ncomp,1);

    Box  const& domain = geom.Domain();

    MFItInfo mfi_info;

    if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(aofs,mfi_info); mfi.isValid(); ++mfi)
    {
        auto const& bx = mfi.tilebox();

        AMREX_D_TERM( Array4<Real> fx = xfluxes.array(mfi,fluxes_comp);,
                      Array4<Real> fy = yfluxes.array(mfi,fluxes_comp);,
                      Array4<Real> fz = zfluxes.array(mfi,fluxes_comp););

	AMREX_D_TERM( Array4<Real> xed = xedge.array(mfi,edge_comp);,
                      Array4<Real> yed = yedge.array(mfi,edge_comp);,
                      Array4<Real> zed = zedge.array(mfi,edge_comp););

        AMREX_D_TERM( Array4<Real const> uc = ucorr.const_array(mfi);,
                      Array4<Real const> vc = vcorr.const_array(mfi);,
                      Array4<Real const> wc = wcorr.const_array(mfi););

        // Compute edge state if needed
        if (!known_edgestate)
        {
            Array4<Real const> const q = state.const_array(mfi,state_comp);

            AMREX_D_TERM( Array4<Real const> u = umac.const_array(mfi);,
                          Array4<Real const> v = vmac.const_array(mfi);,
                          Array4<Real const> w = wmac.const_array(mfi););

            ComputeEdgeState( bx, AMREX_D_DECL( xed, yed, zed ), q, ncomp,
                              AMREX_D_DECL( u, v, w ), domain, bcs,
                              d_bcrec_ptr);

        }

        // Compute fluxes
        HydroUtils::ComputeFluxes(bx, AMREX_D_DECL(fx,fy,fz),
                                  AMREX_D_DECL(uc,vc,wc),
                                  AMREX_D_DECL(xed,yed,zed), geom, ncomp );

        // Compute divergence
        Box tmpbox = amrex::surroundingNodes(bx);
        int tmpcomp = ncomp*AMREX_SPACEDIM;
        FArrayBox tmpfab(tmpbox, tmpcomp);
        Elixir eli = tmpfab.elixir();
        Array4<Real> divtmp_arr = tmpfab.array(ncomp*AMREX_SPACEDIM);
        HydroUtils::ComputeDivergence(bx, divtmp_arr,
                                      AMREX_D_DECL(fx,fy,fz),
                                      AMREX_D_DECL( xed, yed, zed ),
                                      AMREX_D_DECL( uc, vc, wc ),
                                      ncomp, geom, iconserv.data() );

        // Sum contribution to sync aofs
        auto const& aofs_arr = aofs.array(mfi, aofs_comp);

        amrex::ParallelFor(bx, ncomp, [aofs_arr, divtmp_arr]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            aofs_arr( i, j, k, n ) += divtmp_arr( i, j, k, n );
        });
    }
}

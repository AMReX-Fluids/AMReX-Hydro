/** \addtogroup Utilities
 * @{
 */

#include <hydro_utils.H>
#include <hydro_utils_K.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_VisMF.H>
#ifdef AMREX_USE_EB
#include <AMReX_EBMultiFabUtil.H>
#endif

using namespace amrex;

void
HydroUtils::create_umac_grown (int lev, int nGrow, BoxArray& fine_grids,
                               const Geometry& fine_geom,
                               const Array<MultiFab*,AMREX_SPACEDIM> u_mac_crse,
                               const Array<MultiFab*,AMREX_SPACEDIM> u_mac_fine,
                               const IntVect& crse_ratio)
{
    BL_PROFILE("HydroUtils::create_umac_grown()");

    if (lev > 0)
    {
        BoxList bl = amrex::GetBndryCells(fine_grids,nGrow);

        BoxArray f_bnd_ba(std::move(bl));

        BoxArray c_bnd_ba = f_bnd_ba; c_bnd_ba.coarsen(crse_ratio);

        c_bnd_ba.maxSize(32);

        f_bnd_ba = c_bnd_ba; f_bnd_ba.refine(crse_ratio);

        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            //
            // crse_src & fine_src must have same parallel distribution.
            // We'll use the KnapSack distribution for the fine_src_ba.
            // Since fine_src_ba should contain more points, this'll lead
            // to a better distribution.
            //
            BoxArray crse_src_ba(c_bnd_ba), fine_src_ba(f_bnd_ba);

            crse_src_ba.surroundingNodes(idim);
            fine_src_ba.surroundingNodes(idim);

            const int N = fine_src_ba.size();

            std::vector<long> wgts(N);

#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (int i = 0; i < N; i++)
                wgts[i] = fine_src_ba[i].numPts();

            DistributionMapping dm;
            // This DM won't be put into the cache.
            dm.KnapSackProcessorMap(wgts,ParallelDescriptor::NProcs());

            // FIXME
            // Declaring in this way doesn't work. I think it's because the box arrays
            // have been changed and each src box is not completely contained within a
            // single box in the Factory's BA
            // For now, coarse-fine boundary doesn't intersect EB, so should be okay...
            // MultiFab crse_src(crse_src_ba, dm, 1, 0, MFInfo(), getLevel(lev-1).Factory());
            // MultiFab fine_src(fine_src_ba, dm, 1, 0, MFInfo(), Factory());
            MultiFab crse_src(crse_src_ba, dm, 1, 0);
            MultiFab fine_src(fine_src_ba, dm, 1, 0);

            crse_src.setVal(1.e200);
            fine_src.setVal(1.e200);
            //
            // We want to fill crse_src from lower level u_mac including u_mac's grow cells.
            //
            const MultiFab& u_macLL = *u_mac_crse[idim];
            crse_src.ParallelCopy(u_macLL,0,0,1,u_macLL.nGrow(),0);

	    const amrex::GpuArray<int,AMREX_SPACEDIM> c_ratio = {AMREX_D_DECL(crse_ratio[0],crse_ratio[1],crse_ratio[2])};

	    //
	    // Fill fine values with piecewise-constant interp of coarse data.
	    // Operate only on faces that overlap--ie, only fill the fine faces that make up each
	    // coarse face, leave the in-between faces alone.
	    //
#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(crse_src); mfi.isValid(); ++mfi)
            {
                const Box& box       = crse_src[mfi].box();
                auto const& crs_arr  = crse_src.array(mfi);
                auto const& fine_arr = fine_src.array(mfi);

                ParallelFor(box,[crs_arr,fine_arr,idim,c_ratio]
                AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
		   int idx[3] = {AMREX_D_DECL(i*c_ratio[0],j*c_ratio[1],k*c_ratio[2])};
#if ( AMREX_SPACEDIM == 2 )
                   // dim1 are the complement of idim
                   int dim1 = ( idim == 0 ) ? 1 : 0;
                   for (int n1 = 0; n1 < c_ratio[dim1]; n1++) {
                      int id[3] = {idx[0],idx[1]};
                      id[dim1] += n1;
                      fine_arr(id[0],id[1],0) = crs_arr(i,j,k);
                   }
#elif ( AMREX_SPACEDIM == 3 )
                   // dim1 and dim2 are the complements of idim
                   int dim1 = ( idim != 0 ) ? 0 : 1 ;
                   int dim2 = ( idim != 0 ) ? ( ( idim == 2 ) ? 1 : 2 ) : 2 ;
                   for (int n1 = 0; n1 < c_ratio[dim1]; n1++) {
                      for (int n2 = 0; n2 < c_ratio[dim2]; n2++) {
                         int id[3] = {idx[0],idx[1],idx[2]};
                         id[dim1] += n1;
                         id[dim2] += n2;
                         fine_arr(id[0],id[1],id[2]) = crs_arr(i,j,k);
                      }
                   }
#endif
                });
            }
            crse_src.clear();
            //
            // Replace pc-interpd fine data with preferred u_mac data at
            // this level u_mac valid only on surrounding faces of valid
            // region - this op will not fill grow region.
            //
            fine_src.ParallelCopy(*u_mac_fine[idim]);
            //
            // Interpolate unfilled grow cells using best data from
            // surrounding faces of valid region, and pc-interpd data
            // on fine faces overlaying coarse edges.
            //
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(fine_src); mfi.isValid(); ++mfi)
            {
                const int  nComp = 1;
                const Box& fbox  = fine_src[mfi].box();
                auto const& fine_arr = fine_src.array(mfi);

		if (fbox.type(0) == IndexType::NODE)
	        {
		  AMREX_HOST_DEVICE_PARALLEL_FOR_4D(fbox,nComp,i,j,k,n,
		  {
		    face_interp_x(i,j,k,n,fine_arr,c_ratio);
		  });
		}
		else if (fbox.type(1) == IndexType::NODE)
		{
		  AMREX_HOST_DEVICE_PARALLEL_FOR_4D(fbox,nComp,i,j,k,n,
		  {
		    face_interp_y(i,j,k,n,fine_arr,c_ratio);
		  });
		}
#if (AMREX_SPACEDIM == 3)
		else
		{
		  AMREX_HOST_DEVICE_PARALLEL_FOR_4D(fbox,nComp,i,j,k,n,
                  {
		    face_interp_z(i,j,k,n,fine_arr,c_ratio);
		  });
		}
#endif
            }

            MultiFab u_mac_save(u_mac_fine[idim]->boxArray(),u_mac_fine[idim]->DistributionMap(),1,0,MFInfo(),
                                u_mac_fine[idim]->Factory());
            u_mac_save.ParallelCopy(*u_mac_fine[idim]);
            u_mac_fine[idim]->ParallelCopy(fine_src,0,0,1,0,nGrow);
            u_mac_fine[idim]->ParallelCopy(u_mac_save);
        }
    }

    for (int n = 0; n < BL_SPACEDIM; ++n)
    {
	u_mac_fine[n]->FillBoundary(fine_geom.periodicity());
    }
}

void
HydroUtils::create_constrained_umac_grown (int lev, int nGrow, BoxArray& fine_grids,
                                           const Geometry* crse_geom, const Geometry* fine_geom,
                                           const Array<MultiFab*,AMREX_SPACEDIM> u_mac_crse,
                                           const Array<MultiFab*,AMREX_SPACEDIM> u_mac_fine,
                                           const MultiFab* divu,
                                           const IntVect& crse_ratio)
{
    BL_PROFILE("HydroUtils::create_constrained_umac_grown()");

    if ( lev > 0) {

       // Check for divU
       int has_divu = (divu != nullptr);

       // Set interpolator
       Interpolater* mapper = &face_divfree_interp;

       // Set BCRec for Umac
       Vector<BCRec> bcrec(1);
       for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
            bcrec[0].setLo(idim,BCType::foextrap);
            bcrec[0].setHi(idim,BCType::foextrap);
       }
       Array<Vector<BCRec>,AMREX_SPACEDIM> bcrecArr = {AMREX_D_DECL(bcrec,bcrec,bcrec)};

       // Set PhysBCFunct. FCFillExtDirDummy is defined in Hydro_utils.H and is not expected to be called.
       PhysBCFunct<GpuBndryFuncFab<FCFillExtDirDummy>> crse_bndry_func(*crse_geom, bcrec, FCFillExtDirDummy{});
       Array<PhysBCFunct<GpuBndryFuncFab<FCFillExtDirDummy>>,AMREX_SPACEDIM> cbndyFuncArr = {AMREX_D_DECL(crse_bndry_func,crse_bndry_func,crse_bndry_func)};

       PhysBCFunct<GpuBndryFuncFab<FCFillExtDirDummy>> fine_bndry_func(*fine_geom, bcrec, FCFillExtDirDummy{});
       Array<PhysBCFunct<GpuBndryFuncFab<FCFillExtDirDummy>>,AMREX_SPACEDIM> fbndyFuncArr = {AMREX_D_DECL(fine_bndry_func,fine_bndry_func,fine_bndry_func)};

       // Squirel away the interior umac
       Array<MultiFab,AMREX_SPACEDIM> u_mac_tmp;
       for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
          u_mac_tmp[idim].define(u_mac_fine[idim]->boxArray(), u_mac_fine[idim]->DistributionMap(),
                                 u_mac_fine[idim]->nComp(), 0);
          MultiFab::Copy(u_mac_tmp[idim], *u_mac_fine[idim], 0, 0, 1, 0);
       }

       // Interpolate umac_fine from coarse data, including ghost cells
       InterpFromCoarseLevel(u_mac_fine, IntVect(nGrow), 0.0,
                             u_mac_crse, 0, 0, 1,
                             *crse_geom, *fine_geom,
                             cbndyFuncArr, 0, fbndyFuncArr, 0,
                             crse_ratio, mapper, bcrecArr, 0);

       // Overwritte valid region with data squireled above
       for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
          MultiFab::Copy(*u_mac_fine[idim], u_mac_tmp[idim], 0, 0, 1, 0);
       }

       // Fill boundary before going further
       for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
       {
           u_mac_fine[idim]->FillBoundary(fine_geom->periodicity());
       }

       // When copying back valid data, the umac divergence of the first ghost cell was
       // modified. Correct this first ghost cell outer (w.r. to the valid region) face velocity to
       // recover the divu at CF boundary.

       // Use mask to flgg C-F boundary ghost cells.
       iMultiFab mask(fine_grids, u_mac_fine[0]->DistributionMap(), 1, 1, MFInfo(),
                      DefaultFabFactory<IArrayBox>());
       // Flags
       int finebnd = 0;
       int crsebnd = 1;
       int physbnd = 0;
       int interior = 0;
       mask.BuildMask(fine_geom->Domain(), fine_geom->periodicity(),
                      finebnd, crsebnd, physbnd, interior);

       const GpuArray<Real,AMREX_SPACEDIM> dxinv = fine_geom->InvCellSizeArray();
       const GpuArray<Real,AMREX_SPACEDIM> dx = fine_geom->CellSizeArray();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(mask,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            auto const& divuarr = (has_divu) ? divu->const_array(mfi) : u_mac_fine[0]->array(mfi);
            auto const& maskarr = mask.const_array(mfi);
            Array<Array4<Real>, AMREX_SPACEDIM> const &umac_arr = {AMREX_D_DECL(u_mac_fine[0]->array(mfi),
                                                                                u_mac_fine[1]->array(mfi),
                                                                                u_mac_fine[2]->array(mfi))};
            for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
                // Grow the box in the direction of the faces we will correct
                const Box& gbx = mfi.growntilebox(IntVect::TheDimensionVector(idim));
                amrex::ParallelFor(gbx, [idim, bx, divuarr, maskarr, crsebnd, umac_arr, dx, dxinv, has_divu]
                AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // Only works on ghost cells flagged as C-F boundaries
                    if ( !bx.contains(i,j,k) && (maskarr(i,j,k) == crsebnd) ) {
                        // Get the divU components from the transverse velocities (!= idim)
                        GpuArray<int,3> idx = {AMREX_D_DECL(i,j,k)};
                        Real transverseTerm = 0.0;
                        for (int trdim = 0; trdim < AMREX_SPACEDIM; trdim++) {
                            if (trdim != idim) {
                                GpuArray<int,3> idxp1 = {AMREX_D_DECL(i,j,k)};
                                idxp1[trdim]++;
                                transverseTerm +=  (  umac_arr[trdim](idxp1[0], idxp1[1], idxp1[2])
                                                    - umac_arr[trdim](idx[0], idx[1], idx[2]) ) * dxinv[trdim];
                            }
                        }
                        // Correct the outer umac face
                        GpuArray<int,3> idxp1 = {AMREX_D_DECL(i,j,k)};
                        idxp1[idim]++;
                        if ( idx[idim] < bx.smallEnd(idim) ) {
                            umac_arr[idim](i,j,k) =  umac_arr[idim](idxp1[0], idxp1[1], idxp1[2])
                                                   + dx[idim] * transverseTerm;
                            if (has_divu) umac_arr[idim](i,j,k) -= dx[idim] * divuarr(i,j,k);
                        } else if (idx[idim] > bx.bigEnd(idim)) {
                            umac_arr[idim](idxp1[0], idxp1[1], idxp1[2]) =  umac_arr[idim](i,j,k)
                                                                          - dx[idim] * transverseTerm;
                            if (has_divu) umac_arr[idim](idxp1[0], idxp1[1], idxp1[2]) += dx[idim] * divuarr(i,j,k);
                        }
                    }
                });
            }
        }
    }

    // Fill boundary for all the levels
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        u_mac_fine[idim]->FillBoundary(fine_geom->periodicity());
    }

}
/** @}*/

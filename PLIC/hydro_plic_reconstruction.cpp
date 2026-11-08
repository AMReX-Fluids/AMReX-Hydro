/**
 * \file hydro_plic_reconstruction.cpp
 * \addtogroup PLIC
 * @{
 */

#include <hydro_plic.H>
#include <hydro_plic_normals_K.H>

#include <string>

using namespace amrex::literals;

namespace PLIC {

template <typename Normal>
void FitPlanesImpl (amrex::MultiFab& plane, amrex::MultiFab const& F, int ngrow)
{
    BL_PROFILE("PLIC::FitPlanes()");
    AMREX_ALWAYS_ASSERT(plane.nComp() == 4 && plane.nGrow() >= ngrow);
    AMREX_ALWAYS_ASSERT(F.nGrow() >= ngrow + Normal::stencil);

    plane.setVal(0.0_rt);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(plane); mfi.isValid(); ++mfi) {
        const amrex::Box gbx = amrex::grow(mfi.validbox(), ngrow);
        auto const& Fa = F.const_array(mfi);
        auto const& pl = plane.array(mfi);
        amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (is_mixed(Fa(i, j, k))) {
                amrex::Real mx, my, mz, alpha;
                fit_plane<Normal>(i, j, k, Fa, mx, my, mz, alpha);
                pl(i, j, k, 0) = mx;
                pl(i, j, k, 1) = my;
                pl(i, j, k, 2) = mz;
                pl(i, j, k, 3) = alpha;
            }
        });
    }
    amrex::Gpu::streamSynchronize();
}

// One more than the stencil, for the first ghost layer of the planes. The advection schemes read no
// further than that.
int nGrowF (std::string const& reconstruction)
{
    return dispatch_normal(reconstruction, [] (auto normal) { return decltype(normal)::stencil + 1; });
}

void FitPlanes (amrex::MultiFab& plane, amrex::MultiFab const& F, int ngrow,
                std::string const& reconstruction)
{
    dispatch_normal(reconstruction, [&] (auto normal) {
        FitPlanesImpl<decltype(normal)>(plane, F, ngrow);
    });
}

} // namespace PLIC
/** @}*/

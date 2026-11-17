/**
 * \file hydro_plic_reconstruction.cpp
 * \addtogroup PLIC
 * @{
 */

#include <hydro_plic.H>
#include <hydro_plic_normals_K.H>

namespace PLIC {

template <typename Normal>
void FitPlanesImpl (amrex::MultiFab& plane, amrex::MultiFab const& F, int ngrow)
{
    BL_PROFILE("PLIC::FitPlanes()");
    AMREX_ALWAYS_ASSERT(plane.nComp() == 4 && plane.nGrow() >= ngrow);
    AMREX_ALWAYS_ASSERT(F.nGrow() >= ngrow + Normal::stencil);
    AMREX_ALWAYS_ASSERT(plane.boxArray() == F.boxArray() &&
                        plane.DistributionMap() == F.DistributionMap());

    // The ghost layers beyond ngrow; the kernel writes every cell within it
    if (plane.nGrow() > ngrow) {
        plane.setVal(0.0_rt);
    }
    auto const& Fa = F.const_arrays();
    auto const& pl = plane.arrays();
    amrex::ParallelFor(plane, amrex::IntVect(ngrow),
                       [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept
    {
        amrex::Real mx = 0.0_rt;
        amrex::Real my = 0.0_rt;
        amrex::Real mz = 0.0_rt;
        amrex::Real alpha = 0.0_rt;
        if (is_mixed(Fa[b](i, j, k))) {
            fit_plane<Normal>(i, j, k, Fa[b], mx, my, mz, alpha);
        }
        pl[b](i, j, k, 0) = mx;
        pl[b](i, j, k, 1) = my;
        pl[b](i, j, k, 2) = mz;
        pl[b](i, j, k, 3) = alpha;
    });
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

/**
 * \file hydro_plic.cpp
 * \addtogroup PLIC
 * @{
 */

#include <hydro_plic.H>

using namespace amrex::literals;

namespace PLIC {

// Neither the plane geometry nor the flux volumes account for cut cells or for R-Z volumes
void CheckMesh (amrex::MultiFab const& F, amrex::Geometry const& geom)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(geom.IsCartesian(), "PLIC requires Cartesian coordinates");
#ifdef AMREX_USE_EB
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!F.hasEBFabFactory(), "PLIC does not support embedded boundaries");
#else
    amrex::ignore_unused(F);
#endif
}

void
SharpFaceDensity (FacePtrs const& rho_face,
                  amrex::Vector<amrex::MultiFab const*> const& F,
                  amrex::Geometry const& geom,
                  amrex::Real rho_l, amrex::Real rho_g,
                  std::string const& reconstruction)
{
    BL_PROFILE("PLIC::SharpFaceDensity()");
    AMREX_ALWAYS_ASSERT(!F.empty());

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        rho_face[dir]->setVal(0.0_rt);
    }

    // Segment liquid fraction, averaged over the given F
    const amrex::Real weight = 1.0_rt / static_cast<amrex::Real>(F.size());
    amrex::MultiFab plane(F[0]->boxArray(), F[0]->DistributionMap(), 4, 1);
    for (auto const* Fi : F) {
        CheckMesh(*Fi, geom);
        // The planes of the cells on either side of each face
        FitPlanes(plane, *Fi, 1, reconstruction);
        auto const& F_ma = Fi->const_arrays();
        auto const& pl_ma = plane.const_arrays();
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            auto const& fc_ma = rho_face[dir]->arrays();
            amrex::ParallelFor(*rho_face[dir], [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
            {
                fc_ma[box_no](i, j, k) += weight * face_segment_fraction(i, j, k, dir, F_ma[box_no], pl_ma[box_no]);
            });
        }
        amrex::Gpu::streamSynchronize();
    }

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        auto const& fc_ma = rho_face[dir]->arrays();
        amrex::ParallelFor(*rho_face[dir], [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
        {
            fc_ma[box_no](i, j, k) = mixture(fc_ma[box_no](i, j, k), rho_l, rho_g);
        });
        amrex::Gpu::streamSynchronize();
        // A face must have one coefficient, or the projection operator is neither conservative nor
        // symmetric. Keep one copy against round-off.
        rho_face[dir]->OverrideSync(geom.periodicity());
    }
}

} // namespace PLIC
/** @}*/

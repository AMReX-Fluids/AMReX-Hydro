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

    // The kernels index F and the planes by the box numbers of rho_face
    amrex::BoxArray const& ba = F[0]->boxArray();
    amrex::DistributionMapping const& dm = F[0]->DistributionMap();
    for (auto const* Fi : F) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(Fi->boxArray() == ba && Fi->DistributionMap() == dm,
                                         "PLIC::SharpFaceDensity: the F must share one layout");
    }
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            rho_face[dir]->boxArray() == amrex::convert(ba, amrex::IntVect::TheDimensionVector(dir)) &&
            rho_face[dir]->DistributionMap() == dm,
            "PLIC::SharpFaceDensity: rho_face must be on the faces of the boxes of F");
    }

    // Segment liquid fraction, averaged over the given F, then the mixture, on the valid faces only:
    // the ghost faces are left as they are
    const auto nF = static_cast<int>(F.size());
    const amrex::Real weight = 1.0_rt / static_cast<amrex::Real>(nF);
    amrex::MultiFab plane(ba, dm, 4, 1);
    for (int n = 0; n < nF; ++n) {
        CheckMesh(*F[n], geom);
        // The planes of the cells on either side of each face
        FitPlanes(plane, *F[n], 1, reconstruction);
        auto const& F_ma = F[n]->const_arrays();
        auto const& pl_ma = plane.const_arrays();
        const bool first = n == 0;
        const bool last = n + 1 == nF;
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            auto const& fc_ma = rho_face[dir]->arrays();
            amrex::ParallelFor(*rho_face[dir], [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
            {
                const amrex::Real sum = (first ? 0.0_rt : fc_ma[box_no](i, j, k))
                    + weight * face_segment_fraction(i, j, k, dir, F_ma[box_no], pl_ma[box_no]);
                fc_ma[box_no](i, j, k) = last ? mixture(sum, rho_l, rho_g) : sum;
            });
        }
    }

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        // A face must have one coefficient, or the projection operator is neither conservative nor
        // symmetric. Keep one copy against round-off.
        rho_face[dir]->OverrideSync(geom.periodicity());
    }
}

} // namespace PLIC
/** @}*/

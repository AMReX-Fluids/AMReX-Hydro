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
                  PhaseDensity const& rho_l, PhaseDensity const& rho_g,
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
    for (auto const* rho : {&rho_l, &rho_g}) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            !rho->field || (rho->field->boxArray() == ba && rho->field->DistributionMap() == dm &&
                            rho->field->nGrowVect().allGE(amrex::IntVect(1))),
            "PLIC::SharpFaceDensity: a density field must have the layout of F and a ghost layer");
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
        // A phase density on a face: its value, or the mean of the two cells beside the face
        const amrex::Real rl_value = rho_l.value;
        const amrex::Real rg_value = rho_g.value;
        const auto rl_ma = rho_l.field ? rho_l.field->const_arrays() : amrex::MultiArray4<amrex::Real const>{};
        const auto rg_ma = rho_g.field ? rho_g.field->const_arrays() : amrex::MultiArray4<amrex::Real const>{};
        const bool first = n == 0;
        const bool last = n + 1 == nF;
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            auto const& fc_ma = rho_face[dir]->arrays();
            const amrex::IntVect e = amrex::IntVect::TheDimensionVector(dir);
            amrex::ParallelFor(*rho_face[dir], [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
            {
                const amrex::Real sum = (first ? 0.0_rt : fc_ma[box_no](i, j, k))
                    + weight * face_segment_fraction(i, j, k, dir, F_ma[box_no], pl_ma[box_no]);
                if (!last) {
                    fc_ma[box_no](i, j, k) = sum;
                    return;
                }
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                const amrex::Real rl = rl_ma ? 0.5_rt * (rl_ma[box_no](iv - e) + rl_ma[box_no](iv)) : rl_value;
                const amrex::Real rg = rg_ma ? 0.5_rt * (rg_ma[box_no](iv - e) + rg_ma[box_no](iv)) : rg_value;
                fc_ma[box_no](i, j, k) = mixture(sum, rl, rg);
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

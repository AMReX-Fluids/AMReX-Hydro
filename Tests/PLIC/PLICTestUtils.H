#ifndef PLIC_TEST_UTILS_H_
#define PLIC_TEST_UTILS_H_

// What the test solver needs around the PLIC library and a solver would do its own way: prescribed
// velocities, the time step, refinement and the liquid volume.

#include <hydro_plic.H>

#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParReduce.H>
#include <AMReX_TagBox.H>

#include <limits>

namespace testutils {

using namespace amrex::literals;

/**
 * \brief Velocity normal to the low face of cell iv in direction dir, from a vector potential
 *        u = curl(A): the circulation of A around the face divided by its area, so that the face
 *        velocities are discretely divergence-free to round-off. In 3D the edge integrals use
 *        3-point Gauss quadrature, which makes coarse face fluxes equal to the sum of the
 *        underlying fine face fluxes for potentials up to degree 5.
 *
 * \param potential potential(x, time, A) with Real x[AMREX_SPACEDIM] and Real A[3]. In 2D, A[2] is
 *                  the stream function and the other components are not read.
 */
template <typename Potential>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
face_velocity_from_potential(const int dir,
                             const amrex::IntVect& iv,
                             const amrex::Real time,
                             amrex::GeometryData const& geomdata,
                             Potential const& potential) {
    const amrex::Real* plo = geomdata.ProbLo();
    const amrex::Real* dx = geomdata.CellSize();
    amrex::Real A[3] = {0.0_rt, 0.0_rt, 0.0_rt};
#if AMREX_SPACEDIM == 2
    // Stream function at the two nodes bounding the face
    const int t = 1 - dir;
    amrex::Real x[2];
    x[dir] = plo[dir] + static_cast<amrex::Real>(iv[dir]) * dx[dir];
    x[t] = plo[t] + static_cast<amrex::Real>(iv[t]) * dx[t];
    potential(x, time, A);
    const amrex::Real psi_lo = A[2];
    x[t] += dx[t];
    potential(x, time, A);
    const amrex::Real psi_hi = A[2];
    // u = dpsi/dy, v = -dpsi/dx
    return (dir == 0 ? 1.0_rt : -1.0_rt) * (psi_hi - psi_lo) / dx[t];
#else
    // Right-handed tangential directions
    const int t1 = (dir + 1) % 3;
    const int t2 = (dir + 2) % 3;
    constexpr amrex::Real gx[3] = {-0.7745966692414834_rt, 0.0_rt, 0.7745966692414834_rt};
    constexpr amrex::Real gw[3] = {5.0_rt / 18.0_rt, 8.0_rt / 18.0_rt, 5.0_rt / 18.0_rt};
    amrex::Real lo[3];
    for (int d = 0; d < 3; ++d) {
        lo[d] = plo[d] + static_cast<amrex::Real>(iv[d]) * dx[d];
    }
    // Integral of A_c along the edge parallel to c, offset by `off` (0 or 1 cell) in the other
    // tangential direction o
    auto edge = [&](const int c, const int o, const int off) {
        amrex::Real x[3] = {lo[0], lo[1], lo[2]};
        x[o] += static_cast<amrex::Real>(off) * dx[o];
        amrex::Real sum = 0.0_rt;
        for (int q = 0; q < 3; ++q) {
            x[c] = lo[c] + 0.5_rt * (1.0_rt + gx[q]) * dx[c];
            potential(x, time, A);
            sum += gw[q] * A[c];
        }
        return sum * dx[c];
    };
    const amrex::Real circ = edge(t1, t2, 0) - edge(t1, t2, 1) + edge(t2, t1, 1) - edge(t2, t1, 0);
    return circ / (dx[t1] * dx[t2]);
#endif
}

/**
 * \brief Prescribed face velocities u = curl(A), on all the faces of umac including ghost faces.
 *        They are discretely divergence-free to round-off and consistent across coarse/fine faces,
 *        which makes them suited to verifying a scheme.
 *
 * \param potential potential(x, time, A), callable on the device
 */
template <typename Potential>
void FaceVelocityFromPotential (PLIC::FacePtrs const& umac, amrex::Geometry const& geom, amrex::Real time,
                                Potential const& potential)
{
    const auto geomdata = geom.data();
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        auto const& um = umac[dir]->arrays();
        amrex::ParallelFor(*umac[dir], umac[dir]->nGrowVect(),
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
        {
            um[box_no](i, j, k) = face_velocity_from_potential(
                dir, amrex::IntVect(AMREX_D_DECL(i, j, k)), time, geomdata, potential);
        });
    }
    amrex::Gpu::streamSynchronize();
}

/**
 * \brief Tag the cells within nband cells of the interface, to keep it on the finest level: the flux
 *        volumes of a geometric unsplit scheme do not tile across coarse/fine boundaries, which must
 *        stay in pure cells.
 *
 * \param F  volume fraction with nband filled ghost cells
 */
inline void
TagInterface (amrex::TagBoxArray& tags, amrex::MultiFab const& F, int nband = 2)
{
    AMREX_ALWAYS_ASSERT(F.nGrow() >= nband);
    auto const& tag_ma = tags.arrays();
    auto const& F_ma = F.const_arrays();
    amrex::ParallelFor(tags, [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
    {
        if (PLIC::interface_band(i, j, k, F_ma[box_no], nband)) {
            tag_ma[box_no](i, j, k) = amrex::TagBox::SET;
        }
    });
    amrex::Gpu::streamSynchronize();
}

//! The liquid volume on the composite grid, not counting covered cells
inline amrex::Real
LiquidVolume (amrex::Vector<amrex::MultiFab const*> const& F,
              amrex::Vector<amrex::Geometry> const& geom,
              amrex::Vector<amrex::IntVect> const& ref_ratio)
{
    amrex::Real volume = 0.0_rt;
    const int nlevs = static_cast<int>(F.size());
    for (int lev = 0; lev < nlevs; ++lev) {
        amrex::iMultiFab mask;
        if (lev < nlevs - 1) {
            mask = amrex::makeFineMask(*F[lev], F[lev + 1]->boxArray(), ref_ratio[lev], 1, 0);
        } else {
            mask.define(F[lev]->boxArray(), F[lev]->DistributionMap(), 1, 0);
            mask.setVal(1);
        }
        auto const& F_ma = F[lev]->const_arrays();
        auto const& mask_ma = mask.const_arrays();
        amrex::Real sum = amrex::ParReduce(
            amrex::TypeList<amrex::ReduceOpSum>{}, amrex::TypeList<amrex::Real>{}, *F[lev],
            amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) -> amrex::GpuTuple<amrex::Real>
            {
                return {(mask_ma[box_no](i, j, k) != 0) ? F_ma[box_no](i, j, k) : 0.0_rt};
            });
        amrex::ParallelDescriptor::ReduceRealSum(sum);
        const amrex::Real* dx = geom[lev].CellSize();
        volume += sum * AMREX_D_TERM(dx[0], *dx[1], *dx[2]);
    }
    return volume;
}

/**
 * \brief The largest time step for a Courant number cfl, from the cell-centered velocity. The
 *        advection uses the face velocity at t^{n+1/2}, which can be larger: keep cfl well below 1.
 */
inline amrex::Real
EstimateDt (amrex::Vector<amrex::MultiFab const*> const& vel, int vel_comp,
            amrex::Vector<amrex::Geometry> const& geom, amrex::Real cfl)
{
    amrex::Real estdt = std::numeric_limits<amrex::Real>::max();
    for (int lev = 0; lev < static_cast<int>(vel.size()); ++lev) {
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            const amrex::Real umax = vel[lev]->norm0(vel_comp + dir, 0, true);
            if (umax > 0.0_rt) {
                estdt = amrex::min(estdt, cfl * geom[lev].CellSize(dir) / umax);
            }
        }
    }
    amrex::ParallelDescriptor::ReduceRealMin(estdt);
    return estdt;
}

} // namespace testutils

#endif

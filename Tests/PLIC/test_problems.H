#ifndef PLIC_TEST_PROBLEMS_H_
#define PLIC_TEST_PROBLEMS_H_

// The prescribed flows and initial shapes of the advection and momentum checks

#include "test_util.H"

#include <AMReX_GpuQualifiers.H>
#include <AMReX_Math.H>

#include <cmath>

namespace plic_test {

enum class Prob { Zalesak, Vortex, Translate };

struct ProbParm
{
    Prob type = Prob::Vortex;
    Real radius = 0.15_rt;
    GpuArray<Real, AMREX_SPACEDIM> center = {AMREX_D_DECL(0.5_rt, 0.75_rt, 0.5_rt)};
    GpuArray<Real, AMREX_SPACEDIM> vel = {AMREX_D_DECL(1.0_rt, 1.0_rt, 1.0_rt)};
    Real slot_width = 0.05_rt;
    Real slot_length = 0.25_rt;
    Real period = 4.0_rt;
};

// The vector potential of the prescribed velocity. Back to the initial shape at t = period. It is
// its value at t = 0 times time_factor(t), so that the velocity is too.
struct Potential
{
    ProbParm p;

    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    Real time_factor (Real time) const noexcept
    {
        return (p.type == Prob::Vortex) ? std::cos(Math::pi<Real>() * time / p.period) : 1.0_rt;
    }

    AMREX_GPU_HOST_DEVICE
    void operator() (const Real* x, Real time, Real* A) const noexcept
    {
        constexpr Real pi = Math::pi<Real>();
        A[0] = A[1] = A[2] = 0.0_rt;
        if (p.type == Prob::Zalesak) {
            const Real omega = 2.0_rt * pi / p.period;
            A[2] = -0.5_rt * omega * ((x[0] - 0.5_rt) * (x[0] - 0.5_rt) + (x[1] - 0.5_rt) * (x[1] - 0.5_rt));
        } else if (p.type == Prob::Vortex) {
            const Real ct = time_factor(time);
            const Real sx2 = std::sin(pi * x[0]) * std::sin(pi * x[0]);
            const Real sy2 = std::sin(pi * x[1]) * std::sin(pi * x[1]);
#if AMREX_SPACEDIM == 2
            A[2] = sx2 * sy2 / pi * ct;
#else
            const Real sz2 = std::sin(pi * x[2]) * std::sin(pi * x[2]);
            A[1] = -sx2 * sz2 / pi * std::sin(2.0_rt * pi * x[1]) * ct;
            A[2] = sx2 * sy2 / pi * std::sin(2.0_rt * pi * x[2]) * ct;
#endif
        } else {
            auto const& U = p.vel;
#if AMREX_SPACEDIM == 2
            A[2] = U[0] * x[1] - U[1] * x[0];
#else
            A[0] = 0.5_rt * (U[1] * x[2] - U[2] * x[1]);
            A[1] = 0.5_rt * (U[2] * x[0] - U[0] * x[2]);
            A[2] = 0.5_rt * (U[0] * x[1] - U[1] * x[0]);
#endif
        }
    }
};

// phi > 0 in the liquid. The Zalesak disk is extruded along z.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
initial_fraction (int i, int j, int k, GeometryData const& geomdata, ProbParm const& p) noexcept
{
    auto const c = p.center;
    const Real R = p.radius;
    const bool slotted = (p.type == Prob::Zalesak);
    const Real hw = 0.5_rt * p.slot_width;
    const Real slot_top = c[1] - R + p.slot_length;
    auto phi = [=] (Real x, Real y, Real z) -> Real {
        Real r2 = (x - c[0]) * (x - c[0]) + (y - c[1]) * (y - c[1]);
#if AMREX_SPACEDIM == 3
        if (!slotted) { r2 += (z - c[2]) * (z - c[2]); }
#else
        amrex::ignore_unused(z);
#endif
        if (r2 > R * R) { return -1.0_rt; }
        if (slotted && std::abs(x - c[0]) < hw && y < slot_top) { return -1.0_rt; }
        return 1.0_rt;
    };
    return Clamp(PLIC::cell_fraction(i, j, k, geomdata, phi, 16), 0.0_rt, 1.0_rt);
}

} // namespace plic_test

#endif

// Mass and momentum fluxes: a drop in a uniform and in a vortex flow, with a density ratio of 1000

#include "test_problems.H"
#include "test_util.H"

#include <PLICTestUtils.H>

#include <AMReX_MultiFabUtil.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <string>

namespace plic_test {

namespace {

// The face arrays of a box, in each direction, as the flux functions take them
#define PLIC_FACES(set, acc) AMREX_D_DECL(set[0].acc(mfi), set[1].acc(mfi), set[2].acc(mfi))

using FaceMFs = Array<MultiFab, AMREX_SPACEDIM>;

// The fluxes of one way of calling the library
struct FaceSet
{
    FaceMFs flux, face, mflux, pflux;

    FaceSet (BoxArray const& ba, DistributionMapping const& dm)
    {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            const BoxArray fba = convert(ba, IntVect::TheDimensionVector(d));
            flux[d].define(fba, dm, 1, 0);
            face[d].define(fba, dm, 1, 0);
            mflux[d].define(fba, dm, 1, 0);
            pflux[d].define(fba, dm, AMREX_SPACEDIM, 0);
        }
    }
};

// Largest |a - b| over the valid cells and the components of two MultiFabs, with a NaN as the largest
Real max_diff (MultiFab const& a, MultiFab const& b)
{
    const int nc = a.nComp();
    MultiFab d(a.boxArray(), a.DistributionMap(), nc, 0);
    MultiFab::Copy(d, a, 0, 0, nc, 0);
    MultiFab::Subtract(d, b, 0, 0, nc, 0);
    return max_over(d, valid_only, [=] AMREX_GPU_DEVICE (Array4<Real const> const& x, int i, int j, int k) noexcept {
        Real e = 0.0_rt;
        for (int c = 0; c < nc; ++c) {
            const Real ec = std::abs(x(i, j, k, c));
            if (!(ec <= e)) { e = ec; }
        }
        return e;
    });
}

Real max_diff (FaceMFs const& a, FaceMFs const& b)
{
    Real diff = 0.0_rt;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        diff = std::max(diff, max_diff(a[d], b[d]));
    }
    return diff;
}

// Largest difference between the fluxes of two sets, the liquid ones only where both have them
Real max_diff (FaceSet const& a, FaceSet const& b, bool liquid)
{
    Real diff = std::max(max_diff(a.mflux, b.mflux), max_diff(a.pflux, b.pflux));
    if (liquid) {
        diff = std::max({diff, max_diff(a.flux, b.flux), max_diff(a.face, b.face)});
    }
    return diff;
}

} // namespace

// A drop in a uniform and in a vortex flow, with a density ratio of 1000, on one periodic level of
// several boxes: F, the density and the momentum advanced together, as a solver with momentum in
// conservative form would. The phase densities are values, then fields that vary across the domain.
void check_momentum ()
{
    Print() << "momentum\n";
    const int n = 64;
    const Geometry geom = unit_geometry(n, IntVect(1));
    BoxArray ba(geom.Domain());
    ba.maxSize(16);
    const DistributionMapping dm(ba);
    const std::string reconstruction = "myc_height";
    const Real rho_l = 1000.0_rt;
    const Real rho_g = 1.0_rt;
    const Real dt = 0.25_rt / static_cast<Real>(n);
    const int nsteps = 8;
    const auto dxinv = geom.InvCellSizeArray();

    for (int test_case = 0; test_case < 4; ++test_case) {
        const bool uniform = (test_case % 2 == 0);
        const bool fields = (test_case >= 2);
        ProbParm prob;
        prob.type = uniform ? Prob::Translate : Prob::Vortex;
        prob.center = {AMREX_D_DECL(0.5_rt, 0.5_rt, 0.5_rt)};

        MultiFab F(ba, dm, 1, PLIC::nGrowF(reconstruction));
        const auto geomdata = geom.data();
        auto const& Fm = F.arrays();
        ParallelFor(F, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept {
            Fm[b](i, j, k) = initial_fraction(i, j, k, geomdata, prob);
        });
        Gpu::streamSynchronize();
        F.FillBoundary(geom.periodicity());

        // The phase density fields: varying, or uniform at the values, which must give the fluxes of
        // the values to the bit
        MultiFab rl_field(ba, dm, 1, 1);
        MultiFab rg_field(ba, dm, 1, 1);
        {
            auto const& rl = rl_field.arrays();
            auto const& rg = rg_field.arrays();
            ParallelFor(rl_field, IntVect(1), [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept {
                amrex::ignore_unused(k);
                const Real x = geomdata.ProbLo(0) + (static_cast<Real>(i) + 0.5_rt) * geomdata.CellSize(0);
                const Real y = geomdata.ProbLo(1) + (static_cast<Real>(j) + 0.5_rt) * geomdata.CellSize(1);
                rl[b](i, j, k) = fields ? rho_l * (1.0_rt + 0.2_rt * x) : rho_l;
                rg[b](i, j, k) = fields ? rho_g * (1.0_rt + 0.5_rt * y) : rho_g;
            });
            Gpu::streamSynchronize();
        }
        const PLIC::PhaseDensity pd_l = fields ? PLIC::PhaseDensity(rl_field) : PLIC::PhaseDensity(rho_l);
        const PLIC::PhaseDensity pd_g = fields ? PLIC::PhaseDensity(rg_field) : PLIC::PhaseDensity(rho_g);

        // The fluxes of the separate calls (sep), of one set of cuts completed by
        // FarFieldMomentumFluxes after the update of F (cuts), of the values from uniform fields
        // (uniform, mass and momentum only) and of one call with vel_face (one_call)
        FaceSet sep(ba, dm), cuts(ba, dm), uniform_fields(ba, dm), one_call(ba, dm);
        FaceMFs um, vface;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            const BoxArray fba = convert(ba, IntVect::TheDimensionVector(d));
            um[d].define(fba, dm, 1, 1);
            vface[d].define(fba, dm, AMREX_SPACEDIM, 0);
        }
        testutils::FaceVelocityFromPotential(GetArrOfPtrs(um), geom, 0.0_rt, Potential{prob});

        // The velocity, density and momentum at the start
        MultiFab vel(ba, dm, AMREX_SPACEDIM, 1);
        average_face_to_cellcenter(vel, 0, GetArrOfConstPtrs(um));
        vel.FillBoundary(geom.periodicity());
        MultiFab vel0(ba, dm, AMREX_SPACEDIM, 0);
        MultiFab::Copy(vel0, vel, 0, 0, AMREX_SPACEDIM, 0);
        MultiFab rho(ba, dm, 1, 0);
        MultiFab mom(ba, dm, AMREX_SPACEDIM, 0);
        {
            auto const& r = rho.arrays();
            auto const& m = mom.arrays();
            auto const& v = vel.const_arrays();
            auto const& Fc = F.const_arrays();
            auto const& rl = rl_field.const_arrays();
            auto const& rg = rg_field.const_arrays();
            ParallelFor(rho, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept {
                r[b](i, j, k) = PLIC::mixture(Fc[b](i, j, k), rl[b](i, j, k), rg[b](i, j, k));
                for (int c = 0; c < AMREX_SPACEDIM; ++c) { m[b](i, j, k, c) = r[b](i, j, k) * v[b](i, j, k, c); }
            });
            Gpu::streamSynchronize();
        }
        Array<Real, AMREX_SPACEDIM> mom0{};
        for (int c = 0; c < AMREX_SPACEDIM; ++c) { mom0[c] = mom.sum(c); }

        Real rho_err = 0.0_rt;
        Real two_stage_diff = 0.0_rt;
        Real field_diff = 0.0_rt;
        Real one_call_diff = 0.0_rt;
        for (int step = 0; step < nsteps; ++step) {
            for (MFIter mfi(F); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.validbox();
                // The upwind velocity, for the faces away from the interface
                auto const& v = vel.const_array(mfi);
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    auto const& vf = vface[d].array(mfi);
                    auto const& u = um[d].const_array(mfi);
                    const IntVect e = IntVect::TheDimensionVector(d);
                    ParallelFor(surroundingNodes(bx, d), AMREX_SPACEDIM,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int c) noexcept {
                        amrex::ignore_unused(k);
                        const IntVect iv(AMREX_D_DECL(i, j, k));
                        vf(iv, c) = (u(iv) >= 0.0_rt) ? v(iv - e, c) : v(iv, c);
                    });
                }
                auto const& Fa = F.const_array(mfi);
                auto const& va = vel.const_array(mfi);
                auto const rl = pd_l.const_array(mfi);
                auto const rg = pd_g.const_array(mfi);
                PLIC::ComputeEdgeStateAndFluxes(bx, PLIC_FACES(sep.flux, array), PLIC_FACES(sep.face, array), Fa,
                                                PLIC_FACES(um, const_array), geom, dt, reconstruction);
                PLIC::ComputeMassMomentumFluxes(bx, PLIC_FACES(sep.mflux, array), PLIC_FACES(sep.pflux, array),
                                                Fa, va, PLIC_FACES(vface, const_array),
                                                PLIC_FACES(um, const_array), geom, dt, rl, rg, reconstruction);

                // One set of cuts, before the face velocity is known; the faces away from the
                // interface follow the update of F
                PLIC::ComputeEdgeStateAndFluxes(bx, PLIC_FACES(cuts.flux, array), PLIC_FACES(cuts.face, array),
                                                PLIC_FACES(cuts.mflux, array), PLIC_FACES(cuts.pflux, array),
                                                Fa, va, PLIC_FACES(um, const_array), geom, dt, rl, rg,
                                                reconstruction);

                PLIC::ComputeEdgeStateAndFluxes(bx, PLIC_FACES(one_call.flux, array),
                                                PLIC_FACES(one_call.face, array),
                                                PLIC_FACES(one_call.mflux, array),
                                                PLIC_FACES(one_call.pflux, array), Fa, va,
                                                PLIC_FACES(vface, const_array), PLIC_FACES(um, const_array),
                                                geom, dt, rl, rg, reconstruction);

                if (!fields) {
                    PLIC::ComputeMassMomentumFluxes(bx, PLIC_FACES(uniform_fields.mflux, array),
                                                    PLIC_FACES(uniform_fields.pflux, array), Fa, va,
                                                    PLIC_FACES(vface, const_array), PLIC_FACES(um, const_array),
                                                    geom, dt, rl_field.const_array(mfi),
                                                    rg_field.const_array(mfi), reconstruction);
                }
            }
            Gpu::streamSynchronize();
            one_call_diff = std::max(one_call_diff, max_diff(sep, one_call, true));
            if (!fields) {
                field_diff = std::max(field_diff, max_diff(sep, uniform_fields, false));
            }
            for (MFIter mfi(F); mfi.isValid(); ++mfi) {
                auto const& Fa = F.array(mfi);
                auto const& r = rho.array(mfi);
                auto const& m = mom.array(mfi);
                auto const& v = vel.array(mfi);
                AMREX_D_TERM(auto const& fx = sep.flux[0].const_array(mfi);,
                             auto const& fy = sep.flux[1].const_array(mfi);,
                             auto const& fz = sep.flux[2].const_array(mfi););
                AMREX_D_TERM(auto const& mx = sep.mflux[0].const_array(mfi);,
                             auto const& my = sep.mflux[1].const_array(mfi);,
                             auto const& mz = sep.mflux[2].const_array(mfi););
                AMREX_D_TERM(auto const& px = sep.pflux[0].const_array(mfi);,
                             auto const& py = sep.pflux[1].const_array(mfi);,
                             auto const& pz = sep.pflux[2].const_array(mfi););
                ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Fa(i, j, k) -= dt * (AMREX_D_TERM((fx(i + 1, j, k) - fx(i, j, k)) * dxinv[0],
                                                      + (fy(i, j + 1, k) - fy(i, j, k)) * dxinv[1],
                                                      + (fz(i, j, k + 1) - fz(i, j, k)) * dxinv[2]));
                    r(i, j, k) -= dt * (AMREX_D_TERM((mx(i + 1, j, k) - mx(i, j, k)) * dxinv[0],
                                                     + (my(i, j + 1, k) - my(i, j, k)) * dxinv[1],
                                                     + (mz(i, j, k + 1) - mz(i, j, k)) * dxinv[2]));
                    for (int c = 0; c < AMREX_SPACEDIM; ++c) {
                        m(i, j, k, c) -= dt * (AMREX_D_TERM((px(i + 1, j, k, c) - px(i, j, k, c)) * dxinv[0],
                                                            + (py(i, j + 1, k, c) - py(i, j, k, c)) * dxinv[1],
                                                            + (pz(i, j, k + 1, c) - pz(i, j, k, c)) * dxinv[2]));
                        v(i, j, k, c) = m(i, j, k, c) / r(i, j, k);
                    }
                });
            }
            Gpu::streamSynchronize();

            // The faces away from the interface of the one set of cuts, once F has been updated in
            // place, as a solver that advances F first would
            for (MFIter mfi(F); mfi.isValid(); ++mfi) {
                PLIC::FarFieldMomentumFluxes(mfi.validbox(), PLIC_FACES(cuts.pflux, array),
                                             PLIC_FACES(cuts.mflux, const_array), PLIC_FACES(vface, const_array));
            }
            Gpu::streamSynchronize();
            two_stage_diff = std::max(two_stage_diff, max_diff(sep, cuts, true));
            F.FillBoundary(geom.periodicity());
            vel.FillBoundary(geom.periodicity());

            // With the values, the density from the mass fluxes is the mixture density of the new F
            if (fields) { continue; }
            auto const& r = rho.const_arrays();
            auto const& Fc = F.const_arrays();
            const Real err = ParReduce(TypeList<ReduceOpMax>{}, TypeList<Real>{}, rho, IntVect(0),
                [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) -> GpuTuple<Real> {
                    return {std::abs(r[b](i, j, k) - PLIC::mixture(Fc[b](i, j, k), rho_l, rho_g))};
                });
            rho_err = max(rho_err, err / rho_l);
        }
        ParallelDescriptor::ReduceRealMax(rho_err);

        Print() << (uniform ? "  uniform flow" : "  vortex flow") << ", density ratio 1000, "
                << (fields ? "phase densities varying fields, " : "") << nsteps << " steps\n";
        if (!fields) {
            Print() << "    largest |rho - mixture(F)| / rho_l " << rho_err << "\n";
            check(rho_err < 100.0_rt * roundoff_tol, "the mass fluxes give the mixture density of the new F");
            Print() << "    largest difference of the fluxes from uniform phase density fields " << field_diff << "\n";
            check(field_diff == 0.0_rt, "uniform phase density fields give the fluxes of the values, to the bit");
        }
        Print() << "    largest difference of the fluxes from one set of cuts " << two_stage_diff << "\n";
        check(two_stage_diff == 0.0_rt, "one set of cuts gives the fluxes of the separate calls, to the bit");
        Print() << "    largest difference of the fluxes from one call " << one_call_diff << "\n";
        check(one_call_diff == 0.0_rt, "one call with vel_face gives the fluxes of the separate calls, to the bit");
        if (uniform) {
            const Real du = max_diff(vel, vel0);
            Print() << "    largest change of the velocity " << du << "\n";
            check(du < 100.0_rt * roundoff_tol, "a uniform velocity stays uniform across the interface");
        } else {
            Real dm_rel = 0.0_rt;
            for (int c = 0; c < AMREX_SPACEDIM; ++c) {
                dm_rel = max(dm_rel, std::abs(mom.sum(c) - mom0[c]) / (rho_l * static_cast<Real>(ba.numPts())));
            }
            Print() << "    total momentum change, per cell and relative to rho_l " << dm_rel << "\n";
            check(dm_rel < roundoff_tol, "the momentum is conserved");
        }
    }
}

#undef PLIC_FACES

} // namespace plic_test

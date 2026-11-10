/**
 * \file hydro_plic_owkes_desjardins.cpp
 * \addtogroup PLIC
 * @{
 */

#include <hydro_plic.H>
#include <hydro_plic_normals_K.H>
#include <hydro_plic_owkes_desjardins_K.H>

#include <AMReX_IArrayBox.H>

using namespace amrex::literals;

namespace PLIC {

namespace detail {

namespace {

bool covers (amrex::Array4<amrex::Real const> const& a, amrex::Box const& b)
{
    const auto alo = amrex::lbound(a);
    const auto ahi = amrex::ubound(a);
    const auto blo = amrex::lbound(b);
    const auto bhi = amrex::ubound(b);
    return alo.x <= blo.x && alo.y <= blo.y && alo.z <= blo.z &&
           ahi.x >= bhi.x && ahi.y >= bhi.y && ahi.z >= bhi.z;
}

// Face area, or 1, times dx / dt: from a volume per unit face area over the step to a flux
amrex::Real ToFlux (amrex::Geometry const& geom, int dir, amrex::Real dt, bool fluxes_are_area_weighted)
{
    amrex::Real area = 1.0_rt;
    if (fluxes_are_area_weighted) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            if (d != dir) { area *= geom.CellSize(d); }
        }
    }
    return area * geom.CellSize(dir) / dt;
}

} // namespace

// The cells of bx and their first ghost layer hold the planes and the interface band read by the
// flux volumes of the faces of bx, which stay in the cells around the face for a Courant number
// below 1. The mesh nodes of bx are traced back from the face velocities around them.
template <typename Normal>
void FluxVolumeData (amrex::Box const& bx,
                     amrex::Array4<amrex::Real const> const& F,
                     amrex::Array<amrex::Array4<amrex::Real const>, AMREX_SPACEDIM> const& umac,
                     amrex::Geometry const& geom, amrex::Real dt,
                     amrex::IArrayBox& band_fab, amrex::FArrayBox& plane_fab, amrex::FArrayBox& disp_fab)
{
    const amrex::Box cells = amrex::grow(bx, 1);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(covers(F, amrex::grow(cells, Normal::stencil)),
                                     "PLIC: F needs PLIC::nGrowF() filled ghost cells");

    // Cells near the interface, and the plane of each mixed cell
    band_fab.resize(cells, 1, amrex::The_Async_Arena());
    plane_fab.resize(cells, 4, amrex::The_Async_Arena());
    auto const& band = band_fab.array();
    auto const& plane = plane_fab.array();
    amrex::ParallelFor(cells, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        band(i, j, k) = static_cast<int>(interface_band(i, j, k, F));
        amrex::Real mx = 0.0_rt, my = 0.0_rt, mz = 0.0_rt, alpha = 0.0_rt;
        if (is_mixed(F(i, j, k))) {
            fit_plane<Normal>(i, j, k, F, mx, my, mz, alpha);
        }
        plane(i, j, k, 0) = mx;
        plane(i, j, k, 1) = my;
        plane(i, j, k, 2) = mz;
        plane(i, j, k, 3) = alpha;
    });

    // Displacement of each node of bx over the step. It depends only on the face velocities around
    // the node, so that the boxes sharing a node trace it the same way and the flux volumes tile
    // across boxes. At coarse/fine boundaries the solver's ghost faces are used.
    const amrex::Box& domain = geom.Domain();
    const amrex::Box nodes = amrex::surroundingNodes(bx);
    amrex::GpuArray<amrex::Array4<amrex::Real const>, AMREX_SPACEDIM> um_all{AMREX_D_DECL(umac[0], umac[1], umac[2])};
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dtdx_all{};
    amrex::GpuArray<int, AMREX_SPACEDIM> dom_lo{}, dom_hi{}, periodic{};
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        dtdx_all[d] = dt * geom.InvCellSize(d);
        dom_lo[d] = domain.smallEnd(d);
        dom_hi[d] = domain.bigEnd(d);
        periodic[d] = static_cast<int>(geom.isPeriodic(d));
    }
    disp_fab.resize(nodes, 3, amrex::The_Async_Arena());
    auto const& disp = disp_fab.array();
    amrex::ParallelFor(nodes, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        amrex::Real dsp[3];
        node_displacement(i, j, k, um_all, dtdx_all, dom_lo, dom_hi, periodic, dsp);
        for (int d = 0; d < 3; ++d) {
            disp(i, j, k, d) = dsp[d];
        }
    });
}

// The fluxes of the faces of bx from one set of cuts of the flux volumes. Each output is optional
// (an empty Array4): the liquid flux and edge state of the volume fraction, and with MassMomentum
// the mass and momentum fluxes. vel_face, the face velocity away from the interface, is optional
// too: without it the momentum flux of those faces is left for FarFieldMomentumFluxes.
template <typename Normal, bool MassMomentum>
void FluxesImpl (amrex::Box const& bx,
                 amrex::Array<amrex::Array4<amrex::Real>, AMREX_SPACEDIM> const& flux,
                 amrex::Array<amrex::Array4<amrex::Real>, AMREX_SPACEDIM> const& face,
                 amrex::Array<amrex::Array4<amrex::Real>, AMREX_SPACEDIM> const& mflux,
                 amrex::Array<amrex::Array4<amrex::Real>, AMREX_SPACEDIM> const& pflux,
                 amrex::Array4<amrex::Real const> const& F,
                 amrex::Array4<amrex::Real const> const& vel,
                 amrex::Array<amrex::Array4<amrex::Real const>, AMREX_SPACEDIM> const& vel_face,
                 amrex::Array<amrex::Array4<amrex::Real const>, AMREX_SPACEDIM> const& umac,
                 amrex::Geometry const& geom, amrex::Real dt,
                 amrex::Real rho_l, amrex::Real rho_g,
                 bool fluxes_are_area_weighted)
{
    if constexpr (MassMomentum) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(covers(vel, amrex::grow(bx, 1)) && vel.nComp() >= AMREX_SPACEDIM,
                                         "PLIC: vel needs AMREX_SPACEDIM components and one filled ghost layer");
    }
    amrex::IArrayBox band_fab;
    amrex::FArrayBox plane_fab, disp_fab;
    FluxVolumeData<Normal>(bx, F, umac, geom, dt, band_fab, plane_fab, disp_fab);
    const amrex::Box cells = band_fab.box();
    auto const& plane = plane_fab.const_array();
    auto const& disp = disp_fab.const_array();
    auto const& band = band_fab.const_array();

    // Near the interface the fluxes come from the cuts of the flux volume. Away from it the flux
    // volume, which stays in the cells around the face for a Courant number below 1, is all in one
    // phase, and carries the given face velocity.
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        const amrex::Box fbx = amrex::surroundingNodes(bx, dir);
        const amrex::IntVect e = amrex::IntVect::TheDimensionVector(dir);
        const amrex::Real dtdx = dt * geom.InvCellSize(dir);
        const amrex::Real to_flux = ToFlux(geom, dir, dt, fluxes_are_area_weighted);
        auto const& um = umac[dir];
        auto const& fl = flux[dir];
        auto const& fc = face[dir];
        auto const& mf = mflux[dir];
        auto const& pf = pflux[dir];
        auto const& vf = vel_face[dir];
        const bool has_flux = static_cast<bool>(fl);
        const bool has_face = static_cast<bool>(fc);
        const bool has_vel_face = static_cast<bool>(vf);
        amrex::ParallelFor(fbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // nvcc forbids an extended lambda from first capturing a variable inside if constexpr
            amrex::ignore_unused(dir, F, plane, disp, cells, vel, rho_l, rho_g, to_flux, mf, pf, vf,
                                 has_vel_face);
            const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
            const amrex::Real courant = um(iv) * dtdx;
            AMREX_ASSERT(std::abs(courant) < 1.0_rt);
            amrex::Real liquid;
            if (band(iv) != 0 || band(iv - e) != 0) {
                if constexpr (MassMomentum) {
                    amrex::Real m = 0.0_rt;
                    amrex::Real p[AMREX_SPACEDIM];
                    liquid = face_mass_momentum_flux(i, j, k, dir, courant, F, plane, disp, vel, cells,
                                                     rho_l, rho_g, m, p);
                    mf(iv) = to_flux * m;
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        pf(iv, n) = to_flux * p[n];
                    }
                } else {
                    liquid = face_liquid_flux(i, j, k, dir, courant, F, plane, disp, cells);
                }
            } else {
                const bool in_liquid = F(iv) > 0.5_rt;
                liquid = in_liquid ? courant : 0.0_rt;
                if constexpr (MassMomentum) {
                    const amrex::Real m = (in_liquid ? rho_l : rho_g) * courant;
                    mf(iv) = to_flux * m;
                    if (has_vel_face) {
                        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                            pf(iv, n) = mf(iv) * vf(iv, n);
                        }
                    }
                }
            }
            if (has_flux) { fl(iv) = to_flux * liquid; }
            // The liquid fraction of the flux volume, or the upwind value without a flow
            if (has_face) {
                if (courant != 0.0_rt) {
                    fc(iv) = liquid / courant;
                } else {
                    fc(iv) = (um(iv) >= 0.0_rt) ? F(iv - e) : F(iv);
                }
            }
        });
    }
}

} // namespace detail

namespace {

using FaceArrays = amrex::Array<amrex::Array4<amrex::Real>, AMREX_SPACEDIM>;
using ConstFaceArrays = amrex::Array<amrex::Array4<amrex::Real const>, AMREX_SPACEDIM>;

} // namespace

void
ComputeEdgeStateAndFluxes (amrex::Box const& bx,
                           AMREX_D_DECL(amrex::Array4<amrex::Real> const& flux_x,
                                        amrex::Array4<amrex::Real> const& flux_y,
                                        amrex::Array4<amrex::Real> const& flux_z),
                           AMREX_D_DECL(amrex::Array4<amrex::Real> const& face_x,
                                        amrex::Array4<amrex::Real> const& face_y,
                                        amrex::Array4<amrex::Real> const& face_z),
                           amrex::Array4<amrex::Real const> const& F,
                           AMREX_D_DECL(amrex::Array4<amrex::Real const> const& umac,
                                        amrex::Array4<amrex::Real const> const& vmac,
                                        amrex::Array4<amrex::Real const> const& wmac),
                           amrex::Geometry const& geom, amrex::Real dt,
                           std::string const& reconstruction,
                           bool fluxes_are_area_weighted)
{
    BL_PROFILE("PLIC::ComputeEdgeStateAndFluxes()");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(geom.IsCartesian(), "PLIC requires Cartesian coordinates");
    const FaceArrays flux{AMREX_D_DECL(flux_x, flux_y, flux_z)};
    const FaceArrays face{AMREX_D_DECL(face_x, face_y, face_z)};
    const ConstFaceArrays u{AMREX_D_DECL(umac, vmac, wmac)};
    dispatch_normal(reconstruction, [&] (auto normal) {
        detail::FluxesImpl<decltype(normal), false>(bx, flux, face, FaceArrays{}, FaceArrays{}, F, {},
                                                    ConstFaceArrays{}, u, geom, dt, 0.0_rt, 0.0_rt,
                                                    fluxes_are_area_weighted);
    });
}

void
ComputeEdgeStateAndFluxes (amrex::Box const& bx,
                           AMREX_D_DECL(amrex::Array4<amrex::Real> const& flux_x,
                                        amrex::Array4<amrex::Real> const& flux_y,
                                        amrex::Array4<amrex::Real> const& flux_z),
                           AMREX_D_DECL(amrex::Array4<amrex::Real> const& face_x,
                                        amrex::Array4<amrex::Real> const& face_y,
                                        amrex::Array4<amrex::Real> const& face_z),
                           AMREX_D_DECL(amrex::Array4<amrex::Real> const& mflux_x,
                                        amrex::Array4<amrex::Real> const& mflux_y,
                                        amrex::Array4<amrex::Real> const& mflux_z),
                           AMREX_D_DECL(amrex::Array4<amrex::Real> const& pflux_x,
                                        amrex::Array4<amrex::Real> const& pflux_y,
                                        amrex::Array4<amrex::Real> const& pflux_z),
                           amrex::Array4<amrex::Real const> const& F,
                           amrex::Array4<amrex::Real const> const& vel,
                           AMREX_D_DECL(amrex::Array4<amrex::Real const> const& umac,
                                        amrex::Array4<amrex::Real const> const& vmac,
                                        amrex::Array4<amrex::Real const> const& wmac),
                           amrex::Geometry const& geom, amrex::Real dt,
                           amrex::Real rho_l, amrex::Real rho_g,
                           std::string const& reconstruction,
                           bool fluxes_are_area_weighted)
{
    BL_PROFILE("PLIC::ComputeEdgeStateAndFluxes(mass,momentum)");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(geom.IsCartesian(), "PLIC requires Cartesian coordinates");
    const FaceArrays flux{AMREX_D_DECL(flux_x, flux_y, flux_z)};
    const FaceArrays face{AMREX_D_DECL(face_x, face_y, face_z)};
    const FaceArrays mflux{AMREX_D_DECL(mflux_x, mflux_y, mflux_z)};
    const FaceArrays pflux{AMREX_D_DECL(pflux_x, pflux_y, pflux_z)};
    const ConstFaceArrays u{AMREX_D_DECL(umac, vmac, wmac)};
    dispatch_normal(reconstruction, [&] (auto normal) {
        detail::FluxesImpl<decltype(normal), true>(bx, flux, face, mflux, pflux, F, vel, ConstFaceArrays{},
                                                   u, geom, dt, rho_l, rho_g, fluxes_are_area_weighted);
    });
}

void
FarFieldMomentumFluxes (amrex::Box const& bx,
                        AMREX_D_DECL(amrex::Array4<amrex::Real> const& pflux_x,
                                     amrex::Array4<amrex::Real> const& pflux_y,
                                     amrex::Array4<amrex::Real> const& pflux_z),
                        AMREX_D_DECL(amrex::Array4<amrex::Real const> const& mflux_x,
                                     amrex::Array4<amrex::Real const> const& mflux_y,
                                     amrex::Array4<amrex::Real const> const& mflux_z),
                        amrex::Array4<amrex::Real const> const& F,
                        AMREX_D_DECL(amrex::Array4<amrex::Real const> const& vel_face_x,
                                     amrex::Array4<amrex::Real const> const& vel_face_y,
                                     amrex::Array4<amrex::Real const> const& vel_face_z))
{
    BL_PROFILE("PLIC::FarFieldMomentumFluxes()");
    // The faces away from the interface are those that the cuts skipped: neither cell next to the
    // face is in the interface band of F
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(detail::covers(F, amrex::grow(bx, 2)),
                                     "PLIC::FarFieldMomentumFluxes: F needs two filled ghost layers");
    const FaceArrays pflux{AMREX_D_DECL(pflux_x, pflux_y, pflux_z)};
    const ConstFaceArrays mflux{AMREX_D_DECL(mflux_x, mflux_y, mflux_z)};
    const ConstFaceArrays vel_face{AMREX_D_DECL(vel_face_x, vel_face_y, vel_face_z)};
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        const amrex::IntVect e = amrex::IntVect::TheDimensionVector(dir);
        auto const& pf = pflux[dir];
        auto const& mf = mflux[dir];
        auto const& vf = vel_face[dir];
        amrex::ParallelFor(amrex::surroundingNodes(bx, dir), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
            const amrex::IntVect lo = iv - e;
            if (!interface_band(i, j, k, F) && !interface_band(lo[0], lo[1], AMREX_D_PICK(0, 0, lo[2]), F)) {
                for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                    pf(iv, n) = mf(iv) * vf(iv, n);
                }
            }
        });
    }
}

void
ComputeMassMomentumFluxes (amrex::Box const& bx,
                           AMREX_D_DECL(amrex::Array4<amrex::Real> const& mflux_x,
                                        amrex::Array4<amrex::Real> const& mflux_y,
                                        amrex::Array4<amrex::Real> const& mflux_z),
                           AMREX_D_DECL(amrex::Array4<amrex::Real> const& pflux_x,
                                        amrex::Array4<amrex::Real> const& pflux_y,
                                        amrex::Array4<amrex::Real> const& pflux_z),
                           amrex::Array4<amrex::Real const> const& F,
                           amrex::Array4<amrex::Real const> const& vel,
                           AMREX_D_DECL(amrex::Array4<amrex::Real const> const& vel_face_x,
                                        amrex::Array4<amrex::Real const> const& vel_face_y,
                                        amrex::Array4<amrex::Real const> const& vel_face_z),
                           AMREX_D_DECL(amrex::Array4<amrex::Real const> const& umac,
                                        amrex::Array4<amrex::Real const> const& vmac,
                                        amrex::Array4<amrex::Real const> const& wmac),
                           amrex::Geometry const& geom, amrex::Real dt,
                           amrex::Real rho_l, amrex::Real rho_g,
                           std::string const& reconstruction,
                           bool fluxes_are_area_weighted)
{
    BL_PROFILE("PLIC::ComputeMassMomentumFluxes()");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(geom.IsCartesian(), "PLIC requires Cartesian coordinates");
    const FaceArrays mflux{AMREX_D_DECL(mflux_x, mflux_y, mflux_z)};
    const FaceArrays pflux{AMREX_D_DECL(pflux_x, pflux_y, pflux_z)};
    const ConstFaceArrays vface{AMREX_D_DECL(vel_face_x, vel_face_y, vel_face_z)};
    const ConstFaceArrays u{AMREX_D_DECL(umac, vmac, wmac)};
    dispatch_normal(reconstruction, [&] (auto normal) {
        detail::FluxesImpl<decltype(normal), true>(bx, FaceArrays{}, FaceArrays{}, mflux, pflux, F, vel, vface,
                                                   u, geom, dt, rho_l, rho_g, fluxes_are_area_weighted);
    });
}

} // namespace PLIC
/** @}*/

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

// Face area, or 1
amrex::Real FaceArea (amrex::Geometry const& geom, int dir, bool fluxes_are_area_weighted)
{
    amrex::Real area = 1.0_rt;
    if (fluxes_are_area_weighted) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            if (d != dir) { area *= geom.CellSize(d); }
        }
    }
    return area;
}

// The faces of the velocity normal to each direction that node_displacement may read: the domain,
// grown by one face across periodic boundaries. Across a wall the nearest face on the domain boundary
// stands in for the faces beyond it.
amrex::GpuArray<amrex::Box, AMREX_SPACEDIM> VelocityFaces (amrex::Geometry const& geom)
{
    const amrex::Box domain = geom.growPeriodicDomain(1);
    amrex::GpuArray<amrex::Box, AMREX_SPACEDIM> faces{};
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        faces[d] = amrex::surroundingNodes(domain, d);
    }
    return faces;
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
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(amrex::Box(F).contains(amrex::grow(cells, Normal::stencil)),
                                     "PLIC: F needs PLIC::nGrowF() filled ghost cells");
    // The nodes of bx read the faces up to one beyond them (see node_displacement)
    const auto vel_faces = VelocityFaces(geom);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            amrex::Box(umac[d]).contains(amrex::grow(amrex::surroundingNodes(bx, d), 1) & vel_faces[d]),
            "PLIC: the face velocities need 1 filled ghost face");
    }

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
    // across boxes. At coarse/fine boundaries the solver's ghost faces are used. Only the faces with a
    // band cell on either side trace their corners, so only the nodes with a band cell around them are
    // traced; the others are left at 0, unread.
    const amrex::Box nodes = amrex::surroundingNodes(bx);
    amrex::GpuArray<amrex::Array4<amrex::Real const>, AMREX_SPACEDIM> um_all{AMREX_D_DECL(umac[0], umac[1], umac[2])};
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dtdx_all{};
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        dtdx_all[d] = dt * geom.InvCellSize(d);
    }
    disp_fab.resize(nodes, 3, amrex::The_Async_Arena());
    auto const& disp = disp_fab.array();
    auto const& band_c = band_fab.const_array();
    amrex::ParallelFor(nodes, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        bool near = false;
#if AMREX_SPACEDIM == 3
        for (int kk = k - 1; kk <= k; ++kk) {
#else
        {
            const int kk = k;
#endif
            for (int jj = j - 1; jj <= j; ++jj) {
                for (int ii = i - 1; ii <= i; ++ii) {
                    near = near || (band_c(ii, jj, kk) != 0);
                }
            }
        }
        if (!near) {
            for (int d = 0; d < 3; ++d) {
                disp(i, j, k, d) = 0.0_rt;
            }
            return;
        }
        amrex::Real dsp[3];
        node_displacement(i, j, k, um_all, dtdx_all, vel_faces, dsp);
        for (int d = 0; d < 3; ++d) {
            // What keeps every flux volume in the cells around its face, and the cuts within the
            // stack: checked in every build, as a larger step reads beyond the data
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::abs(dsp[d]) < 1.0_rt,
                                             "PLIC: a mesh node moves a cell or more in one step");
            disp(i, j, k, d) = dsp[d];
        }
    });
}

template <typename Normal>
void ComputeEdgeStateAndFluxesImpl (amrex::Box const& bx,
                                    amrex::Array<amrex::Array4<amrex::Real>, AMREX_SPACEDIM> const& flux,
                                    amrex::Array<amrex::Array4<amrex::Real>, AMREX_SPACEDIM> const& face,
                                    amrex::Array4<amrex::Real const> const& F,
                                    amrex::Array<amrex::Array4<amrex::Real const>, AMREX_SPACEDIM> const& umac,
                                    amrex::Geometry const& geom, amrex::Real dt,
                                    bool fluxes_are_area_weighted)
{
    amrex::IArrayBox band_fab;
    amrex::FArrayBox plane_fab, disp_fab;
    FluxVolumeData<Normal>(bx, F, umac, geom, dt, band_fab, plane_fab, disp_fab);
    const amrex::Box cells = band_fab.box();
    auto const& plane = plane_fab.const_array();
    auto const& disp = disp_fab.const_array();
    auto const& band = band_fab.const_array();

    // Liquid crossing each face over the step. Near the interface it comes from the cuts of the flux
    // volume. Away from it the flux volume, which stays in the cells around the face for a Courant
    // number below 1, is all in one phase.
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        const amrex::Box fbx = amrex::surroundingNodes(bx, dir);
        const amrex::IntVect e = amrex::IntVect::TheDimensionVector(dir);
        const amrex::Real dtdx = dt * geom.InvCellSize(dir);
        // From a volume per unit face area over the step to a flux
        const amrex::Real to_flux = FaceArea(geom, dir, fluxes_are_area_weighted) * geom.CellSize(dir) / dt;
        auto const& um = umac[dir];
        auto const& fl = flux[dir];
        auto const& fc = face[dir];
        amrex::ParallelFor(fbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
            const amrex::Real courant = um(iv) * dtdx;
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::abs(courant) < 1.0_rt,
                                             "PLIC: the face Courant number must be below 1");
            amrex::Real liquid;
            if (band(iv) != 0 || band(iv - e) != 0) {
                liquid = face_liquid_flux(i, j, k, dir, courant, F, plane, disp, cells);
            } else {
                liquid = (F(iv) > 0.5_rt) ? courant : 0.0_rt;
            }
            // The liquid fraction of the flux volume, or the upwind value without a flow.
            // Where shear folds the flux volume onto itself its net volume, the Courant number, can
            // be far smaller than the liquid in it: the fraction is then held to [0,1].
            if (courant != 0.0_rt) {
                fc(iv) = amrex::Clamp(liquid / courant, 0.0_rt, 1.0_rt);
            } else {
                fc(iv) = (um(iv) >= 0.0_rt) ? F(iv - e) : F(iv);
            }
            fl(iv) = to_flux * liquid;
        });
    }
}

// Without a step there is no flux volume: the edge state is the upwind value and the flux the face
// velocity times it, as from the other advection schemes
void UpwindEdgeStateAndFluxes (amrex::Box const& bx,
                               amrex::Array<amrex::Array4<amrex::Real>, AMREX_SPACEDIM> const& flux,
                               amrex::Array<amrex::Array4<amrex::Real>, AMREX_SPACEDIM> const& face,
                               amrex::Array4<amrex::Real const> const& F,
                               amrex::Array<amrex::Array4<amrex::Real const>, AMREX_SPACEDIM> const& umac,
                               amrex::Geometry const& geom, bool fluxes_are_area_weighted)
{
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        const amrex::IntVect e = amrex::IntVect::TheDimensionVector(dir);
        const amrex::Real area = FaceArea(geom, dir, fluxes_are_area_weighted);
        auto const& um = umac[dir];
        auto const& fl = flux[dir];
        auto const& fc = face[dir];
        amrex::ParallelFor(amrex::surroundingNodes(bx, dir), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            amrex::ignore_unused(k);
            const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
            fc(iv) = (um(iv) >= 0.0_rt) ? F(iv - e) : F(iv);
            fl(iv) = area * um(iv) * fc(iv);
        });
    }
}

} // namespace detail

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
    const amrex::Array<amrex::Array4<amrex::Real>, AMREX_SPACEDIM> flux{AMREX_D_DECL(flux_x, flux_y, flux_z)};
    const amrex::Array<amrex::Array4<amrex::Real>, AMREX_SPACEDIM> face{AMREX_D_DECL(face_x, face_y, face_z)};
    const amrex::Array<amrex::Array4<amrex::Real const>, AMREX_SPACEDIM> vel{AMREX_D_DECL(umac, vmac, wmac)};
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(dt >= 0.0_rt, "PLIC: the time step must not be negative");
    if (dt == 0.0_rt) {
        detail::UpwindEdgeStateAndFluxes(bx, flux, face, F, vel, geom, fluxes_are_area_weighted);
        return;
    }
    dispatch_normal(reconstruction, [&] (auto normal) {
        detail::ComputeEdgeStateAndFluxesImpl<decltype(normal)>(bx, flux, face, F, vel, geom, dt,
                                                                fluxes_are_area_weighted);
    });
}

} // namespace PLIC
/** @}*/

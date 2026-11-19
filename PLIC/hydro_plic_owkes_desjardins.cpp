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
// grown by one face across periodic boundaries. Beyond any other domain face, whatever its boundary
// condition, the nearest face on the domain boundary stands in for the faces beyond it.
amrex::GpuArray<amrex::Box, AMREX_SPACEDIM> VelocityFaces (amrex::Geometry const& geom)
{
    const amrex::Box domain = geom.growPeriodicDomain(1);
    amrex::GpuArray<amrex::Box, AMREX_SPACEDIM> faces{};
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        faces[d] = amrex::surroundingNodes(domain, d);
    }
    return faces;
}

// The faces normal to dir that an Array4 of them holds. Box(Array4) gives their indices as cells, and
// convert would add a face at the top: keep the indices and set the type only.
template <typename T>
amrex::Box FaceBox (amrex::Array4<T> const& a, int dir)
{
    const amrex::Box b(a);
    return amrex::Box(b.smallEnd(), b.bigEnd(), amrex::IntVect::TheDimensionVector(dir));
}

// Liquid crossing the faces normal to dir over the step, and their edge states. Near the interface
// the liquid comes from the cuts of the flux volume. Away from it the flux volume, which stays in the
// cells around the face for a Courant number below 1, is all in one phase.
struct FaceFluxes
{
    amrex::Array4<amrex::Real> fl;
    amrex::Array4<amrex::Real> fc;
    amrex::Array4<amrex::Real const> F;
    amrex::Array4<amrex::Real const> plane;
    amrex::Array4<amrex::Real const> disp;
    amrex::Array4<amrex::Real const> um;
    amrex::Array4<int const> band;
    amrex::Box cells;
    int dir = 0;
    amrex::Real dtdx = 0.0_rt;
    amrex::Real to_flux = 0.0_rt;

    AMREX_GPU_DEVICE void operator() (int i, int j, int k) const noexcept
    {
        const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
        const amrex::Real courant = um(iv) * dtdx;
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::abs(courant) < 1.0_rt,
                                         "PLIC: the face Courant number must be below 1");
        amrex::Real liquid;
        if (band(iv) != 0 || band(iv - amrex::IntVect::TheDimensionVector(dir)) != 0) {
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
            fc(iv) = upwind_fraction(F, iv, dir, um(iv));
        }
        fl(iv) = to_flux * liquid;
    }
};

// The upwind edge states and fluxes of the faces normal to dir
struct UpwindFaceFluxes
{
    amrex::Array4<amrex::Real> fl;
    amrex::Array4<amrex::Real> fc;
    amrex::Array4<amrex::Real const> F;
    amrex::Array4<amrex::Real const> um;
    int dir = 0;
    amrex::Real area = 0.0_rt;

    AMREX_GPU_DEVICE void operator() (int i, int j, int k) const noexcept
    {
        amrex::ignore_unused(k);
        const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
        const amrex::Real f = upwind_fraction(F, iv, dir, um(iv));
        fc(iv) = f;
        fl(iv) = area * um(iv) * f;
    }
};

// One launch over the faces of bx in every direction, kernel[dir] on those normal to dir
template <typename Kernel>
void ForEachFace (amrex::Box const& bx, amrex::Array<Kernel, AMREX_SPACEDIM> const& kernel)
{
    amrex::ParallelFor(AMREX_D_DECL(amrex::surroundingNodes(bx, 0), amrex::surroundingNodes(bx, 1),
                                    amrex::surroundingNodes(bx, 2)),
                       AMREX_D_DECL(kernel[0], kernel[1], kernel[2]));
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
    // The nodes of bx read the faces up to one beyond them (see node_displacement)
    const auto vel_faces = VelocityFaces(geom);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            FaceBox(umac[d], d).contains(amrex::grow(amrex::surroundingNodes(bx, d), 1) & vel_faces[d]),
            "PLIC: the face velocities need 1 filled ghost face");
    }

    amrex::GpuArray<amrex::Array4<amrex::Real const>, AMREX_SPACEDIM> um_all{AMREX_D_DECL(umac[0], umac[1], umac[2])};
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dtdx_all{};
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        dtdx_all[d] = dt * geom.InvCellSize(d);
    }

    // Cells near the interface, and the plane of each mixed cell. Also the Courant number of the
    // faces the nodes of bx read beyond the faces of bx, which the face kernel checks: each from the
    // ghost cell above it (and the last layer from the cell below). A node moves by a mean of the
    // Courant numbers of the faces it reads, so this keeps every node within a cell of where it was,
    // and every flux volume in the cells around its face, near the interface or not. Checked in every
    // build, as a larger step reads beyond the data.
    const amrex::Box valid = bx;
    band_fab.resize(cells, 1, amrex::The_Async_Arena());
    plane_fab.resize(cells, 4, amrex::The_Async_Arena());
    auto const& band = band_fab.array();
    auto const& plane = plane_fab.array();
    amrex::ParallelFor(cells, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        band(i, j, k) = static_cast<int>(interface_band(i, j, k, F));
        store_plane<Normal>(i, j, k, F, plane);
        const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
        if (valid.contains(iv)) {
            return;
        }
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            const int top = (iv[d] == cells.bigEnd(d)) ? 1 : 0;
            for (int s = 0; s <= top; ++s) {
                const amrex::IntVect f = iv + s * amrex::IntVect::TheDimensionVector(d);
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!vel_faces[d].contains(f) || std::abs(um_all[d](f) * dtdx_all[d]) < 1.0_rt,
                                                 "PLIC: the face Courant number must be below 1");
            }
        }
    });

    // Displacement of each node of bx over the step. It depends only on the face velocities around
    // the node, so that the boxes sharing a node trace it the same way and the flux volumes tile
    // across boxes. At coarse/fine boundaries the solver's ghost faces are used. Only the faces with a
    // band cell on either side trace their corners, so only the nodes with a band cell around them are
    // traced; the others are left at 0, unread.
    const amrex::Box nodes = amrex::surroundingNodes(bx);
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
            // stack. It follows from the Courant numbers checked by the face kernel and above.
            AMREX_ASSERT_WITH_MESSAGE(std::abs(dsp[d]) < 1.0_rt,
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

    amrex::Array<FaceFluxes, AMREX_SPACEDIM> kernels;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        // From a volume per unit face area over the step to a flux
        const amrex::Real to_flux = FaceArea(geom, dir, fluxes_are_area_weighted) * geom.CellSize(dir) / dt;
        kernels[dir] = FaceFluxes{flux[dir], face[dir], F, plane, disp, umac[dir], band, cells, dir,
                                  dt * geom.InvCellSize(dir), to_flux};
    }
    ForEachFace(bx, kernels);
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
    amrex::Array<UpwindFaceFluxes, AMREX_SPACEDIM> kernels;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        kernels[dir] = UpwindFaceFluxes{flux[dir], face[dir], F, umac[dir], dir,
                                        FaceArea(geom, dir, fluxes_are_area_weighted)};
    }
    ForEachFace(bx, kernels);
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
    // Checked with dt = 0 too, so that a call at initialization finds what the first step would:
    // nGrowF rejects an unknown reconstruction
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(amrex::Box(F).contains(amrex::grow(bx, nGrowF(reconstruction))),
                                     "PLIC: F needs PLIC::nGrowF(reconstruction) filled ghost cells");
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        const amrex::Box faces = amrex::surroundingNodes(bx, d);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(detail::FaceBox(flux[d], d).contains(faces) &&
                                         detail::FaceBox(face[d], d).contains(faces) &&
                                         detail::FaceBox(vel[d], d).contains(faces),
                                         "PLIC: the fluxes, edge states and face velocities need the faces of bx");
    }
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

/**
 * \file hydro_plic_owkes_desjardins.cpp
 * \addtogroup PLIC
 * @{
 */

#include <hydro_plic.H>
#include <hydro_plic_normals_K.H>
#include <hydro_plic_owkes_desjardins_K.H>

#include <AMReX_IArrayBox.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Reduce.H>

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

// The faces of bx and one face beyond them, normal to dir, that its flux volumes and traced nodes
// read, within the faces node_displacement may read
amrex::Box ReadFaces (amrex::Box const& bx, int dir, amrex::GpuArray<amrex::Box, AMREX_SPACEDIM> const& vel_faces)
{
    return amrex::grow(amrex::surroundingNodes(bx, dir), 1) & vel_faces[dir];
}

// Add to reduce_data the largest face Courant number, |u| dt/dx, of the faces read by bx normal to dir
void ReduceCourantNumber (amrex::Box const& read, amrex::Array4<amrex::Real const> const& um, amrex::Real dtdx,
                          amrex::ReduceOps<amrex::ReduceOpMax>& reduce_op,
                          amrex::ReduceData<amrex::Real>& reduce_data)
{
    reduce_op.eval(read, reduce_data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> amrex::GpuTuple<amrex::Real>
    {
        return {std::abs(um(i, j, k) * dtdx)};
    });
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
        if (fl) { fl(iv) = to_flux * liquid; }
    }
};

} // namespace

// The cells of bx and their first ghost layer hold the planes and the interface band read by the
// flux volumes of the faces of bx, which stay in the cells around the face for a Courant number
// below 1. The mesh nodes of bx are traced back from the face velocities around them.
template <typename Normal>
void FluxVolumeData (amrex::Box const& bx,
                     amrex::Array4<amrex::Real const> const& F,
                     amrex::Array<amrex::Array4<amrex::Real const>, AMREX_SPACEDIM> const& umac,
                     amrex::GpuArray<amrex::Box, AMREX_SPACEDIM> const& vel_faces,
                     amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dtdx,
                     amrex::IArrayBox& band_fab, amrex::FArrayBox& plane_fab, amrex::FArrayBox& disp_fab)
{
    const amrex::Box cells = amrex::grow(bx, 1);
    amrex::GpuArray<amrex::Array4<amrex::Real const>, AMREX_SPACEDIM> um_all{AMREX_D_DECL(umac[0], umac[1], umac[2])};

    // Cells near the interface, and the plane of each mixed cell
    band_fab.resize(cells, 1, amrex::The_Async_Arena());
    plane_fab.resize(cells, 4, amrex::The_Async_Arena());
    auto const& band = band_fab.array();
    auto const& plane = plane_fab.array();
    amrex::ParallelFor(cells, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        band(i, j, k) = static_cast<int>(interface_band(i, j, k, F));
        store_plane<Normal>(i, j, k, F, plane);
    });

    // Displacement of each node of bx over the step. It depends only on the face velocities around
    // the node, so that the boxes sharing a node trace it the same way and the flux volumes tile
    // across boxes. At coarse/fine boundaries the solver's ghost faces are used. Only the faces with a
    // band cell on either side trace their corners, so only the nodes with a band cell around them are
    // traced; the others, which no face reads, are set to 0.
    const amrex::Box nodes = amrex::surroundingNodes(bx);
    disp_fab.resize(nodes, 3, amrex::The_Async_Arena());
    auto const& disp = disp_fab.array();
    auto const& band_c = band_fab.const_array();
    amrex::ParallelFor(nodes, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        amrex::ignore_unused(k);
        const amrex::IntVect node(AMREX_D_DECL(i, j, k));
        bool near = false;
        amrex::Loop(amrex::Box(node - 1, node), [&] (amrex::IntVect const& c) {
            near = near || (band_c(c) != 0);
        });
        if (!near) {
            for (int d = 0; d < 3; ++d) {
                disp(i, j, k, d) = 0.0_rt;
            }
            return;
        }
        amrex::Real dsp[3];
        // Below a cell in each direction, as a mean of Courant numbers below 1
        node_displacement(i, j, k, um_all, dtdx, vel_faces, dsp);
        for (int d = 0; d < 3; ++d) {
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
                                    amrex::GpuArray<amrex::Box, AMREX_SPACEDIM> const& vel_faces,
                                    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dtdx,
                                    bool fluxes_are_area_weighted)
{
    amrex::IArrayBox band_fab;
    amrex::FArrayBox plane_fab, disp_fab;
    FluxVolumeData<Normal>(bx, F, umac, vel_faces, dtdx, band_fab, plane_fab, disp_fab);
    const amrex::Box cells = band_fab.box();
    auto const& plane = plane_fab.const_array();
    auto const& disp = disp_fab.const_array();
    auto const& band = band_fab.const_array();

    amrex::Array<FaceFluxes, AMREX_SPACEDIM> kernels;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        // From a volume per unit face area over the step to a flux
        const amrex::Real to_flux = FaceArea(geom, dir, fluxes_are_area_weighted) * geom.CellSize(dir) / dt;
        kernels[dir] = FaceFluxes{flux[dir], face[dir], F, plane, disp, umac[dir], band, cells, dir,
                                  dtdx[dir], to_flux};
    }
    // One launch over the faces of every direction
    amrex::ParallelFor(AMREX_D_DECL(amrex::surroundingNodes(bx, 0), amrex::surroundingNodes(bx, 1),
                                    amrex::surroundingNodes(bx, 2)),
                       AMREX_D_DECL(kernels[0], kernels[1], kernels[2]));
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
                           bool fluxes_are_area_weighted,
                           bool check_courant)
{
    BL_PROFILE("PLIC::ComputeEdgeStateAndFluxes()");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(geom.IsCartesian(), "PLIC requires Cartesian coordinates");
    const amrex::Array<amrex::Array4<amrex::Real>, AMREX_SPACEDIM> flux{AMREX_D_DECL(flux_x, flux_y, flux_z)};
    const amrex::Array<amrex::Array4<amrex::Real>, AMREX_SPACEDIM> face{AMREX_D_DECL(face_x, face_y, face_z)};
    const amrex::Array<amrex::Array4<amrex::Real const>, AMREX_SPACEDIM> vel{AMREX_D_DECL(umac, vmac, wmac)};
    // Without a step there is no flux volume, and the flux would divide by zero
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(dt > 0.0_rt, "PLIC: ComputeEdgeStateAndFluxes needs a positive time step");
    // A face or nodal box would lose its top faces, and the geometry of another level would give
    // wrong cell sizes and velocity faces
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(bx.cellCentered() && geom.Domain().contains(bx),
                                     "PLIC: bx must be cell-centered and in the domain of geom");
    // nGrowF rejects an unknown reconstruction
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(amrex::Box(F).contains(amrex::grow(bx, nGrowF(reconstruction))),
                                     "PLIC: F needs PLIC::nGrowF(reconstruction) filled ghost cells");
    // The Courant number of every face read must be below 1. A node then moves by less than a cell,
    // as a mean of the Courant numbers of the faces it reads, and every flux volume stays in the
    // cells around its face, near the interface or not: a larger step would read beyond the data.
    // It is reduced to the host, as an assertion in a kernel is compiled out of a GPU release build,
    // unless the caller has checked the level with CheckCourantNumber: the reduction waits for the
    // GPU, so that boxes checked one by one do not overlap.
    const auto vel_faces = detail::VelocityFaces(geom);
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dtdx{};
    amrex::ReduceOps<amrex::ReduceOpMax> reduce_op;
    amrex::ReduceData<amrex::Real> reduce_data(reduce_op);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        const amrex::Box faces = amrex::surroundingNodes(bx, d);
        // An empty flux array asks for the edge states only
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE((!flux[d] || detail::FaceBox(flux[d], d).contains(faces)) &&
                                         detail::FaceBox(face[d], d).contains(faces),
                                         "PLIC: the fluxes and edge states need the faces of bx");
        // The nodes of bx read the faces up to one beyond them (see node_displacement)
        const amrex::Box read = detail::ReadFaces(bx, d, vel_faces);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(detail::FaceBox(vel[d], d).contains(read),
                                         "PLIC: the face velocities need the faces of bx and 1 filled ghost face");
        dtdx[d] = dt * geom.InvCellSize(d);
        if (check_courant) {
            detail::ReduceCourantNumber(read, vel[d], dtdx[d], reduce_op, reduce_data);
        }
    }
    if (check_courant) {
        const amrex::Real courant = amrex::get<0>(reduce_data.value(reduce_op));
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(courant < 1.0_rt, "PLIC: the face Courant number must be below 1");
    }

    dispatch_normal(reconstruction, [&] (auto normal) {
        detail::ComputeEdgeStateAndFluxesImpl<decltype(normal)>(bx, flux, face, F, vel, geom, dt, vel_faces, dtdx,
                                                                fluxes_are_area_weighted);
    });
}

void
CheckCourantNumber (AMREX_D_DECL(amrex::MultiFab const& umac,
                                 amrex::MultiFab const& vmac,
                                 amrex::MultiFab const& wmac),
                    amrex::Geometry const& geom, amrex::Real dt)
{
    BL_PROFILE("PLIC::CheckCourantNumber()");
    const amrex::Array<amrex::MultiFab const*, AMREX_SPACEDIM> vel{AMREX_D_DECL(&umac, &vmac, &wmac)};
    const auto vel_faces = detail::VelocityFaces(geom);
    amrex::ReduceOps<amrex::ReduceOpMax> reduce_op;
    amrex::ReduceData<amrex::Real> reduce_data(reduce_op);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(vel[d]->ixType() == amrex::IndexType(amrex::IntVect::TheDimensionVector(d)) &&
                                         vel[d]->nGrowVect().allGE(1),
                                         "PLIC: the face velocities need 1 filled ghost face");
        const amrex::Real dtdx = dt * geom.InvCellSize(d);
        // The same faces as ComputeEdgeStateAndFluxes reads for each box, all in one wait for the GPU
        for (amrex::MFIter mfi(*vel[d]); mfi.isValid(); ++mfi) {
            const amrex::Box read = detail::ReadFaces(amrex::enclosedCells(mfi.validbox()), d, vel_faces);
            detail::ReduceCourantNumber(read, vel[d]->const_array(mfi), dtdx, reduce_op, reduce_data);
        }
    }
    amrex::Real courant = amrex::get<0>(reduce_data.value(reduce_op));
    amrex::ParallelDescriptor::ReduceRealMax(courant);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(courant < 1.0_rt, "PLIC: the face Courant number must be below 1");
}

} // namespace PLIC
/** @}*/

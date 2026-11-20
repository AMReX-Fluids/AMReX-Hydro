/**
 * \file hydro_plic_owkes_desjardins.cpp
 * \addtogroup PLIC
 * @{
 */

#include <hydro_plic.H>
#include <hydro_plic_normals_K.H>
#include <hydro_plic_owkes_desjardins_K.H>

#include <AMReX_IArrayBox.H>
#include <AMReX_Reduce.H>

#include <string>
#include <type_traits>

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

// Whether a holds ncomp components on the faces of bx normal to dir
template <typename T>
bool HoldsFaces (amrex::Array4<T> const& a, amrex::Box const& bx, int dir, int ncomp = 1)
{
    return a && a.nComp() >= ncomp && FaceBox(a, dir).contains(amrex::surroundingNodes(bx, dir));
}

// The faces of bx and one face beyond them, normal to dir, that its flux volumes and traced nodes
// read, within the faces node_displacement may read
amrex::Box ReadFaces (amrex::Box const& bx, int dir, amrex::GpuArray<amrex::Box, AMREX_SPACEDIM> const& vel_faces)
{
    return amrex::grow(amrex::surroundingNodes(bx, dir), 1) & vel_faces[dir];
}

// What the flux kernel reads and writes for the mass and momentum fluxes of one direction
struct MassMomentumArgs
{
    amrex::Array4<amrex::Real const> vel;
    PhaseDensityArray rho_l{0.0_rt};
    PhaseDensityArray rho_g{0.0_rt};
    amrex::Array4<amrex::Real> mf;
    amrex::Array4<amrex::Real> pf;
    amrex::Array4<amrex::Real const> vf;
};
struct NoMassMomentumArgs {};

// Liquid crossing the faces normal to dir over the step, and their edge states, and with
// MassMomentum the mass and momentum fluxes. Near the interface they come from the cuts of the flux
// volume. Away from it the flux volume, which stays in the cells around the face for a Courant number
// below 1, is all in the phase of the upwind cell, and carries its density and the given face
// velocity, or far_field_mark without one. Each output may be an empty Array4.
template <bool MassMomentum>
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
    std::conditional_t<MassMomentum, MassMomentumArgs, NoMassMomentumArgs> mm;

    AMREX_GPU_DEVICE void operator() (int i, int j, int k) const noexcept
    {
        const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
        const amrex::Real courant = um(iv) * dtdx;
        amrex::Real liquid;
        if (band(iv) != 0 || band(iv - amrex::IntVect::TheDimensionVector(dir)) != 0) {
            if constexpr (MassMomentum) {
                amrex::Real m = 0.0_rt;
                amrex::Real p[AMREX_SPACEDIM];
                liquid = face_mass_momentum_flux(i, j, k, dir, courant, F, plane, disp, mm.vel, cells,
                                                 mm.rho_l, mm.rho_g, m, p);
                mm.mf(iv) = to_flux * m;
                for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                    mm.pf(iv, n) = to_flux * p[n];
                }
            } else {
                liquid = face_liquid_flux(i, j, k, dir, courant, F, plane, disp, cells);
            }
        } else {
            // Both cells, and the cells around them, are in one phase
            const amrex::IntVect up = upwind_cell(iv, dir, courant);
            const bool in_liquid = F(up) > 0.5_rt;
            liquid = in_liquid ? courant : 0.0_rt;
            if constexpr (MassMomentum) {
                const auto u3 = up.dim3();
                const amrex::Real m =
                    to_flux * (in_liquid ? mm.rho_l(u3.x, u3.y, u3.z) : mm.rho_g(u3.x, u3.y, u3.z)) * courant;
                mm.mf(iv) = m;
                for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                    mm.pf(iv, n) = mm.vf ? m * mm.vf(iv, n) : far_field_mark;
                }
            }
        }
        // The liquid fraction of the flux volume, or the upwind value without a flow.
        // Where shear folds the flux volume onto itself its net volume, the Courant number, can
        // be far smaller than the liquid in it: the fraction is then held to [0,1].
        if (fc) {
            if (courant != 0.0_rt) {
                fc(iv) = amrex::Clamp(liquid / courant, 0.0_rt, 1.0_rt);
            } else {
                fc(iv) = upwind_fraction(F, iv, dir, um(iv));
            }
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

// The fluxes of the faces of bx from one set of cuts of the flux volumes (see FaceFluxes)
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
                 amrex::GpuArray<amrex::Box, AMREX_SPACEDIM> const& vel_faces,
                 amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dtdx,
                 PhaseDensityArray const& rho_l, PhaseDensityArray const& rho_g,
                 bool fluxes_are_area_weighted)
{
    amrex::IArrayBox band_fab;
    amrex::FArrayBox plane_fab, disp_fab;
    FluxVolumeData<Normal>(bx, F, umac, vel_faces, dtdx, band_fab, plane_fab, disp_fab);
    const amrex::Box cells = band_fab.box();
    auto const& plane = plane_fab.const_array();
    auto const& disp = disp_fab.const_array();
    auto const& band = band_fab.const_array();

    amrex::Array<FaceFluxes<MassMomentum>, AMREX_SPACEDIM> kernels;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        // From a volume per unit face area over the step to a flux
        const amrex::Real to_flux = FaceArea(geom, dir, fluxes_are_area_weighted) * geom.CellSize(dir) / dt;
        kernels[dir] = FaceFluxes<MassMomentum>{flux[dir], face[dir], F, plane, disp, umac[dir], band, cells, dir,
                                                dtdx[dir], to_flux, {}};
        if constexpr (MassMomentum) {
            kernels[dir].mm = MassMomentumArgs{vel, rho_l, rho_g, mflux[dir], pflux[dir], vel_face[dir]};
        } else {
            amrex::ignore_unused(mflux, pflux, vel, vel_face, rho_l, rho_g);
        }
    }
    // One launch over the faces of every direction
    amrex::ParallelFor(AMREX_D_DECL(amrex::surroundingNodes(bx, 0), amrex::surroundingNodes(bx, 1),
                                    amrex::surroundingNodes(bx, 2)),
                       AMREX_D_DECL(kernels[0], kernels[1], kernels[2]));
}

namespace {

using FaceArrays = amrex::Array<amrex::Array4<amrex::Real>, AMREX_SPACEDIM>;
using ConstFaceArrays = amrex::Array<amrex::Array4<amrex::Real const>, AMREX_SPACEDIM>;

// What every flux needs: a step whose flux volumes stay in the data, through faces of bx read from
// the face velocities with one ghost face. Returns the velocity faces node_displacement may read and
// dt/dx. nGrowF rejects an unknown reconstruction.
//
// The Courant number of every face read must be below 1. A node then moves by less than a cell, as a
// mean of the Courant numbers of the faces it reads, and every flux volume stays in the cells around
// its face, near the interface or not: a larger step would read beyond the data. It is reduced to
// the host, as an assertion in a kernel is compiled out of a GPU release build.
void CheckStep (amrex::Box const& bx, amrex::Array4<amrex::Real const> const& F, ConstFaceArrays const& umac,
                amrex::Geometry const& geom, amrex::Real dt, std::string const& reconstruction,
                amrex::GpuArray<amrex::Box, AMREX_SPACEDIM>& vel_faces,
                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dtdx)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(geom.IsCartesian(), "PLIC requires Cartesian coordinates");
    // Without a step there is no flux volume, and the flux would divide by zero
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(dt > 0.0_rt, "PLIC: the fluxes need a positive time step");
    // A face or nodal box would lose its top faces, and the geometry of another level would give
    // wrong cell sizes and velocity faces
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(bx.cellCentered() && geom.Domain().contains(bx),
                                     "PLIC: bx must be cell-centered and in the domain of geom");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(amrex::Box(F).contains(amrex::grow(bx, nGrowF(reconstruction))),
                                     "PLIC: F needs PLIC::nGrowF(reconstruction) filled ghost cells");
    vel_faces = VelocityFaces(geom);
    amrex::ReduceOps<amrex::ReduceOpMax> reduce_op;
    amrex::ReduceData<amrex::Real> reduce_data(reduce_op);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        // The nodes of bx read the faces up to one beyond them (see node_displacement)
        const amrex::Box read = ReadFaces(bx, d, vel_faces);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(umac[d] && FaceBox(umac[d], d).contains(read),
                                         "PLIC: the face velocities need the faces of bx and 1 filled ghost face");
        dtdx[d] = dt * geom.InvCellSize(d);
        auto const& um = umac[d];
        const amrex::Real c = dtdx[d];
        reduce_op.eval(read, reduce_data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> amrex::GpuTuple<amrex::Real>
        {
            return {std::abs(um(i, j, k) * c)};
        });
    }
    const amrex::Real courant = amrex::get<0>(reduce_data.value(reduce_op));
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(courant < 1.0_rt, "PLIC: the face Courant number must be below 1");
}

// The fluxes and edge states of the volume fraction
void CheckFluxInputs (amrex::Box const& bx, FaceArrays const& flux, FaceArrays const& face)
{
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(HoldsFaces(flux[d], bx, d) && HoldsFaces(face[d], bx, d),
                                         "PLIC: the fluxes and edge states need the faces of bx");
    }
}

// What the mass and momentum fluxes need besides
void CheckMassMomentumInputs (amrex::Box const& bx, FaceArrays const& mflux, FaceArrays const& pflux,
                              amrex::Array4<amrex::Real const> const& vel, ConstFaceArrays const& vel_face,
                              PhaseDensityArray const& rho_l, PhaseDensityArray const& rho_g)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(vel && amrex::Box(vel).contains(amrex::grow(bx, 1)) &&
                                     vel.nComp() >= AMREX_SPACEDIM,
                                     "PLIC: vel needs AMREX_SPACEDIM components and one filled ghost layer");
    for (auto const* rho : {&rho_l, &rho_g}) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!rho->field || amrex::Box(rho->field).contains(amrex::grow(bx, 1)),
                                         "PLIC: a phase density field needs one filled ghost layer");
    }
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(HoldsFaces(mflux[d], bx, d) && HoldsFaces(pflux[d], bx, d, AMREX_SPACEDIM),
                                         "PLIC: mflux needs 1 component, and pflux AMREX_SPACEDIM, on the faces of bx");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!vel_face[d] || HoldsFaces(vel_face[d], bx, d, AMREX_SPACEDIM),
                                         "PLIC: vel_face needs AMREX_SPACEDIM components on the faces of bx");
    }
}

// The face velocity of the functions that need it on every face
void RequireVelFace (ConstFaceArrays const& vel_face, char const* name)
{
    for (auto const& vf : vel_face) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(vf, (std::string(name) + ": vel_face is needed").c_str());
    }
}

} // namespace

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
    const detail::FaceArrays flux{AMREX_D_DECL(flux_x, flux_y, flux_z)};
    const detail::FaceArrays face{AMREX_D_DECL(face_x, face_y, face_z)};
    const detail::ConstFaceArrays u{AMREX_D_DECL(umac, vmac, wmac)};
    amrex::GpuArray<amrex::Box, AMREX_SPACEDIM> vel_faces;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dtdx;
    detail::CheckStep(bx, F, u, geom, dt, reconstruction, vel_faces, dtdx);
    detail::CheckFluxInputs(bx, flux, face);
    dispatch_normal(reconstruction, [&] (auto normal) {
        detail::FluxesImpl<decltype(normal), false>(bx, flux, face, detail::FaceArrays{}, detail::FaceArrays{}, F,
                                                    {}, detail::ConstFaceArrays{}, u, geom, dt, vel_faces, dtdx,
                                                    0.0_rt, 0.0_rt, fluxes_are_area_weighted);
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
                           PhaseDensityArray const& rho_l, PhaseDensityArray const& rho_g,
                           std::string const& reconstruction,
                           bool fluxes_are_area_weighted)
{
    BL_PROFILE("PLIC::ComputeEdgeStateAndFluxes(mass,momentum)");
    const detail::FaceArrays flux{AMREX_D_DECL(flux_x, flux_y, flux_z)};
    const detail::FaceArrays face{AMREX_D_DECL(face_x, face_y, face_z)};
    const detail::FaceArrays mflux{AMREX_D_DECL(mflux_x, mflux_y, mflux_z)};
    const detail::FaceArrays pflux{AMREX_D_DECL(pflux_x, pflux_y, pflux_z)};
    const detail::ConstFaceArrays u{AMREX_D_DECL(umac, vmac, wmac)};
    amrex::GpuArray<amrex::Box, AMREX_SPACEDIM> vel_faces;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dtdx;
    detail::CheckStep(bx, F, u, geom, dt, reconstruction, vel_faces, dtdx);
    detail::CheckFluxInputs(bx, flux, face);
    detail::CheckMassMomentumInputs(bx, mflux, pflux, vel, detail::ConstFaceArrays{}, rho_l, rho_g);
    dispatch_normal(reconstruction, [&] (auto normal) {
        detail::FluxesImpl<decltype(normal), true>(bx, flux, face, mflux, pflux, F, vel, detail::ConstFaceArrays{},
                                                   u, geom, dt, vel_faces, dtdx, rho_l, rho_g, fluxes_are_area_weighted);
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
                           AMREX_D_DECL(amrex::Array4<amrex::Real const> const& vel_face_x,
                                        amrex::Array4<amrex::Real const> const& vel_face_y,
                                        amrex::Array4<amrex::Real const> const& vel_face_z),
                           AMREX_D_DECL(amrex::Array4<amrex::Real const> const& umac,
                                        amrex::Array4<amrex::Real const> const& vmac,
                                        amrex::Array4<amrex::Real const> const& wmac),
                           amrex::Geometry const& geom, amrex::Real dt,
                           PhaseDensityArray const& rho_l, PhaseDensityArray const& rho_g,
                           std::string const& reconstruction,
                           bool fluxes_are_area_weighted)
{
    BL_PROFILE("PLIC::ComputeEdgeStateAndFluxes(mass,momentum,vel_face)");
    const detail::FaceArrays flux{AMREX_D_DECL(flux_x, flux_y, flux_z)};
    const detail::FaceArrays face{AMREX_D_DECL(face_x, face_y, face_z)};
    const detail::FaceArrays mflux{AMREX_D_DECL(mflux_x, mflux_y, mflux_z)};
    const detail::FaceArrays pflux{AMREX_D_DECL(pflux_x, pflux_y, pflux_z)};
    const detail::ConstFaceArrays vface{AMREX_D_DECL(vel_face_x, vel_face_y, vel_face_z)};
    const detail::ConstFaceArrays u{AMREX_D_DECL(umac, vmac, wmac)};
    detail::RequireVelFace(vface, "PLIC::ComputeEdgeStateAndFluxes");
    amrex::GpuArray<amrex::Box, AMREX_SPACEDIM> vel_faces;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dtdx;
    detail::CheckStep(bx, F, u, geom, dt, reconstruction, vel_faces, dtdx);
    detail::CheckFluxInputs(bx, flux, face);
    detail::CheckMassMomentumInputs(bx, mflux, pflux, vel, vface, rho_l, rho_g);
    dispatch_normal(reconstruction, [&] (auto normal) {
        detail::FluxesImpl<decltype(normal), true>(bx, flux, face, mflux, pflux, F, vel, vface,
                                                   u, geom, dt, vel_faces, dtdx, rho_l, rho_g, fluxes_are_area_weighted);
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
                        AMREX_D_DECL(amrex::Array4<amrex::Real const> const& vel_face_x,
                                     amrex::Array4<amrex::Real const> const& vel_face_y,
                                     amrex::Array4<amrex::Real const> const& vel_face_z))
{
    BL_PROFILE("PLIC::FarFieldMomentumFluxes()");
    const detail::FaceArrays pflux{AMREX_D_DECL(pflux_x, pflux_y, pflux_z)};
    const detail::ConstFaceArrays mflux{AMREX_D_DECL(mflux_x, mflux_y, mflux_z)};
    const detail::ConstFaceArrays vel_face{AMREX_D_DECL(vel_face_x, vel_face_y, vel_face_z)};
    // The faces away from the interface are those that the cuts marked
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(detail::HoldsFaces(pflux[dir], bx, dir, AMREX_SPACEDIM) &&
                                         detail::HoldsFaces(mflux[dir], bx, dir) &&
                                         detail::HoldsFaces(vel_face[dir], bx, dir, AMREX_SPACEDIM),
                                         "PLIC::FarFieldMomentumFluxes: pflux and vel_face need AMREX_SPACEDIM "
                                         "components, and mflux one, on the faces of bx");
        auto const& pf = pflux[dir];
        auto const& mf = mflux[dir];
        auto const& vf = vel_face[dir];
        amrex::ParallelFor(amrex::surroundingNodes(bx, dir), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            amrex::ignore_unused(k); // 2D
            const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
            if (pf(iv, 0) == far_field_mark) {
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
                           PhaseDensityArray const& rho_l, PhaseDensityArray const& rho_g,
                           std::string const& reconstruction,
                           bool fluxes_are_area_weighted)
{
    BL_PROFILE("PLIC::ComputeMassMomentumFluxes()");
    const detail::FaceArrays mflux{AMREX_D_DECL(mflux_x, mflux_y, mflux_z)};
    const detail::FaceArrays pflux{AMREX_D_DECL(pflux_x, pflux_y, pflux_z)};
    const detail::ConstFaceArrays vface{AMREX_D_DECL(vel_face_x, vel_face_y, vel_face_z)};
    const detail::ConstFaceArrays u{AMREX_D_DECL(umac, vmac, wmac)};
    detail::RequireVelFace(vface, "PLIC::ComputeMassMomentumFluxes");
    amrex::GpuArray<amrex::Box, AMREX_SPACEDIM> vel_faces;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dtdx;
    detail::CheckStep(bx, F, u, geom, dt, reconstruction, vel_faces, dtdx);
    detail::CheckMassMomentumInputs(bx, mflux, pflux, vel, vface, rho_l, rho_g);
    dispatch_normal(reconstruction, [&] (auto normal) {
        detail::FluxesImpl<decltype(normal), true>(bx, detail::FaceArrays{}, detail::FaceArrays{}, mflux, pflux, F,
                                                   vel, vface, u, geom, dt, vel_faces, dtdx, rho_l, rho_g,
                                                   fluxes_are_area_weighted);
    });
}

} // namespace PLIC
/** @}*/

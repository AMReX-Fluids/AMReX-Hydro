/**
 * \file hydro_plic.H
 *
 * Geometric volume-of-fluid transport with a piecewise-linear interface calculation (PLIC), for any
 * AMReX-based solver.
 *
 * The interface in each mixed cell is a plane, fitted to the volume fractions F (the liquid fraction
 * of each cell) around it. The reconstruction, how the normal of the plane is computed, is selected
 * by name: "youngs", "myc" or "myc_height". Like Godunov, MOL and BDS, the library provides spatial
 * operators and leaves the update, the ghost cells and the synchronization of the levels to the
 * solver.
 *
 * Requires Cartesian coordinates, and has no embedded boundary support.
 *
 * \addtogroup PLIC
 *  @{
 */

#ifndef HYDRO_PLIC_H
#define HYDRO_PLIC_H

#include <hydro_plic_K.H>

#include <AMReX_Array.H>
#include <AMReX_Array4.H>
#include <AMReX_BCRec.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

#include <string>

namespace PLIC {

/** \brief Pointers to a face-centered quantity, one MultiFab per direction. */
using FacePtrs = amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>;

/**
 * \brief The ghost cells of F that the operators of this library read, with the given
 *        reconstruction. Allocate and fill at least this many.
 */
[[nodiscard]] int nGrowF (std::string const& reconstruction);

/**
 * \brief Abort unless the mesh is supported: Cartesian coordinates and no embedded boundaries.
 */
void CheckMesh (amrex::MultiFab const& F, amrex::Geometry const& geom);

/**
 * \brief Fit the plane of each mixed cell, and zero the others.
 *
 * \param plane           4 components (m_x, m_y, m_z, alpha), in unit-cube cell coordinates, with the
 *                        liquid on the side m.x <= alpha and sum|m| = 1
 * \param F               volume fraction, with its ghost cells filled
 * \param ngrow           ghost layers of plane to fit
 * \param reconstruction  "youngs", "myc" or "myc_height"
 */
void FitPlanes (amrex::MultiFab& plane, amrex::MultiFab const& F, int ngrow,
                std::string const& reconstruction);

/**
 * \brief Fluxes and edge states of the volume fraction on the faces of bx, from the unsplit geometric
 *        fluxes of Owkes & Desjardins (JCP 270, 2014), for the solver's MFIter loop as
 *        Godunov::ComputeEdgeState is.
 *
 * The flux through a face is the liquid inside its flux volume, the region bounded by the face and
 * its corners traced back along the flow over dt, divided by dt: F^{n+1} = F^n - dt div(flux)
 * conserves the liquid volume to round-off and keeps F in [0,1] up to round-off, without clipping.
 * The edge state is the liquid fraction of the flux volume, so that flux = edge state * u (times the
 * face area if fluxes_are_area_weighted), and the upwind F where u = 0.
 *
 * Faces shared by boxes, and periodic images, get the same flux to the bit, given the same face
 * velocities and ghost cells. On several levels, the coarse faces under fine faces take the fine
 * fluxes (average_down_faces with one time step, or refluxing with subcycling), and the interface
 * must stay away from coarse/fine boundaries, as the flux volumes do not tile across them: tag the
 * cells of PLIC::interface_band for refinement.
 *
 * \param bx              cells whose faces get fluxes
 * \param flux_x          output flux on the x faces of bx (likewise y, z), or an empty Array4 for the
 *                        edge states only, as HydroUtils::ComputeFluxesOnBoxFromState asks with
 *                        advection_type "PLIC"
 * \param face_x          output edge state on the x faces of bx (likewise y, z)
 * \param F               volume fraction at t^n, with nGrowF(reconstruction) filled ghost cells
 * \param umac            face velocity at t^{n+1/2}, discretely divergence-free, with one filled ghost
 *                        layer (likewise vmac, wmac). The face Courant number must be below 1.
 * \param geom            geometry of the level
 * \param dt              the time step
 * \param reconstruction  "youngs", "myc" or "myc_height"
 * \param fluxes_are_area_weighted  multiply the fluxes by the face area
 */
void ComputeEdgeStateAndFluxes (amrex::Box const& bx,
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
                                bool fluxes_are_area_weighted = false);

/**
 * \brief Mass and momentum fluxes on the faces of bx, from the same flux volumes as
 *        ComputeEdgeStateAndFluxes, for momentum advected in conservative form with a large density
 *        ratio (Rudman, JCP 1998).
 *
 * Each piece of a flux volume near the interface carries the density of its phase and the velocity of
 * the cell it lies in. Away from the interface the flux volume is in one phase, and carries vel_face.
 * With rho^n = PLIC::mixture(F^n, rho_l, rho_g), the solver updates
 *
 *     rho^{n+1} = rho^n - dt div(mflux),   (rho u)^{n+1} = (rho u)^n - dt div(pflux)
 *
 * and rho^{n+1} is the mixture density of the F^{n+1} of ComputeEdgeStateAndFluxes, to round-off, so
 * that a uniform velocity stays uniform across the interface whatever the density ratio. The fluxes
 * follow the conventions of ComputeEdgeStateAndFluxes, and are the same to the bit on faces shared by
 * boxes.
 *
 * \param bx              cells whose faces get fluxes
 * \param mflux_x         output mass flux on the x faces of bx, 1 component (likewise y, z)
 * \param pflux_x         output momentum flux on the x faces of bx, AMREX_SPACEDIM components
 * \param F               volume fraction at t^n, as given to ComputeEdgeStateAndFluxes
 * \param vel             cell-centered velocity at t^n, AMREX_SPACEDIM components from component 0
 *                        (offset an Array4 for another component), with one filled ghost layer
 * \param vel_face_x      velocity on the x faces of bx, AMREX_SPACEDIM components, used away from the
 *                        interface: an upwind value, or the edge state of the solver's velocity scheme
 * \param umac            face velocity, as given to ComputeEdgeStateAndFluxes (likewise vmac, wmac)
 * \param geom            geometry of the level
 * \param dt              the time step
 * \param rho_l           liquid density
 * \param rho_g           gas density
 * \param reconstruction  as given to ComputeEdgeStateAndFluxes
 * \param fluxes_are_area_weighted  multiply the fluxes by the face area
 */
void ComputeMassMomentumFluxes (amrex::Box const& bx,
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
                                std::string const& reconstruction, bool fluxes_are_area_weighted = false);

/**
 * \brief Overwrite the ghost cells of F beyond the domain faces where F is extrapolated, by extending
 *        the planes of the interior interface. At a wall the interface then meets the wall at whatever
 *        angle it has: a free contact line. Extrapolated values impose 90 degrees instead, and the
 *        inconsistent normals near the wall unbalance gravity at a tilted interface. At an outflow the
 *        interface continues smoothly. The result does not depend on the box layout.
 *
 * The boundary condition of F on each domain face selects what happens there:
 *  - foextrap, hoextrap: the ghost cells are overwritten (a wall or an outflow);
 *  - ext_dir: the ghost cells are kept, as the volume fraction that flows in;
 *  - reflect_even: the ghost cells are kept, as the mirror image of a symmetry plane (90 degrees);
 *  - int_dir: periodic, nothing to do;
 *  - reflect_odd is not a valid condition for a volume fraction, and aborts.
 *
 * \param F               volume fraction with all its ghost cells filled, including those beyond the
 *                        domain (first-order extrapolation is enough at foextrap and hoextrap faces)
 *                        and at coarse/fine boundaries
 * \param geom            geometry of this level
 * \param bcrec           boundary condition of F, as for its ghost cells
 * \param reconstruction  "youngs", "myc" or "myc_height"
 */
void FillWallGhosts (amrex::MultiFab& F, amrex::Geometry const& geom, amrex::BCRec const& bcrec,
                     std::string const& reconstruction);

/**
 * \brief Face density rho_g + (rho_l - rho_g) * <face_segment_fraction>, with the average over the
 *        given volume fractions. For a planar interface normal to gravity it makes rho_face * g a
 *        discrete gradient, which is the hydrostatic balance of a balanced-force projection (Francois
 *        et al., JCP 213, 2006). Each face gets one value.
 *
 * \param rho_face        output, face-centered with any ghost cells, which are not set
 * \param F               the volume fractions to average, e.g. at t^n and t^{n+1} for t^{n+1/2}, with
 *                        nGrowF() filled ghost cells
 * \param geom            geometry of this level
 * \param rho_l           liquid density
 * \param rho_g           gas density
 * \param reconstruction  "youngs", "myc" or "myc_height"
 */
void SharpFaceDensity (FacePtrs const& rho_face,
                       amrex::Vector<amrex::MultiFab const*> const& F,
                       amrex::Geometry const& geom,
                       amrex::Real rho_l, amrex::Real rho_g,
                       std::string const& reconstruction);

} // namespace PLIC

#endif
/** @}*/

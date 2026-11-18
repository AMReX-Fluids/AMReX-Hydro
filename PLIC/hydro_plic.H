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
 * conserves the liquid volume to round-off and, as the flux volumes tile, keeps F in [0,1] up to
 * round-off, without clipping. The exception is a face whose traced face has almost no extent normal
 * to it: its flux volume is held near the face and falls short of the Courant number, and its liquid
 * is the Courant number times the liquid fraction of what remains. A cell full of one phase stays
 * full there, but next to the interface F can leave [0,1] by more than round-off.
 * The edge state is the liquid fraction of the flux volume, in [0,1], so that flux = edge state * u
 * (times the face area if fluxes_are_area_weighted), and the upwind F where u = 0. Where shear folds
 * the flux volume onto itself, so that its net volume nearly cancels, the fraction is held to [0,1]
 * and the flux, which is what conserves F, can differ from edge state * u.
 *
 * Faces shared by boxes, and periodic images, get the same flux to the bit, given the same face
 * velocities and ghost cells. On several levels, the coarse faces under fine faces take the fine
 * fluxes (average_down_faces with one time step, or refluxing with subcycling), and the interface
 * must stay away from coarse/fine boundaries, as the flux volumes do not tile across them: tag the
 * cells of PLIC::interface_band for refinement.
 *
 * \param bx              cells whose faces get fluxes, cell-centered and in the domain of geom. Every
 *                        face of bx is written, so in a tiled MFIter loop two tiles write the faces
 *                        between them (the same values, as with Godunov::ComputeEdgeState): with
 *                        OpenMP, give each tile its own face arrays, or loop without tiling.
 * \param flux_x          output flux on the x faces of bx (likewise y, z), or an empty Array4 for the
 *                        edge states only
 * \param face_x          output edge state on the x faces of bx (likewise y, z)
 * \param F               volume fraction at t^n, with nGrowF(reconstruction) filled ghost cells
 * \param umac            face velocity at t^{n+1/2}, discretely divergence-free, with one filled ghost
 *                        layer (likewise vmac, wmac). No ghost face beyond a domain face that is not
 *                        periodic is read: there the faces on the domain boundary stand in for them,
 *                        whatever the boundary condition. The face Courant number must be below 1, on
 *                        the ghost faces too: checked in every build.
 * \param geom            geometry of the level
 * \param dt              the time step, positive
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
 * \brief Overwrite the ghost cells of F beyond the domain faces where F is extrapolated, by extending
 *        the planes of the interior interface. At a wall the interface then meets the wall at whatever
 *        angle it has: a free contact line. Extrapolated values impose 90 degrees instead, and the
 *        inconsistent normals near the wall unbalance gravity at a tilted interface. At an outflow the
 *        interface continues smoothly. The result does not depend on the box layout, or on the order
 *        of the faces: a ghost cell beyond several (an edge or a corner) takes the mean of their
 *        values, over those that extend a plane to it, which for a planar interface is the plane.
 *
 * The planes are extended beyond the wall cells within wall_nsearch cells, along the face, of the
 * interface; beyond the others a ghost cell takes the value of its wall cell. At a shallow angle the
 * extended plane would cross the ghost layers further from the contact line than that, where no
 * reconstruction stencil of an interface cell reaches. The search reads no cell beyond a face whose
 * ghost cells are kept, and a ghost cell beyond a wall and such a face keeps its value.
 *
 * The boundary condition of F on each non-periodic domain face selects what happens there (that of
 * a periodic direction is not read):
 *  - foextrap, hoextrap, hoextrapcc: the ghost cells are overwritten (a wall or an outflow);
 *  - ext_dir, ext_dir_cc: the ghost cells are kept, as the volume fraction that flows in;
 *  - reflect_even: the ghost cells are kept, as the mirror image of a symmetry plane (90 degrees);
 *  - direction_dependent: the ghost cells are kept (inflow or outflow, face by face);
 *  - int_dir: nothing to do;
 *  - reflect_odd is not a valid condition for a volume fraction, and aborts.
 *
 * \param F               volume fraction with all its ghost cells filled, including those beyond the
 *                        domain (first-order extrapolation is enough at foextrap, hoextrap and
 *                        hoextrapcc faces) and at coarse/fine boundaries
 * \param geom            geometry of this level
 * \param bcrec           boundary condition of F, as for its ghost cells
 * \param reconstruction  "youngs", "myc" or "myc_height"
 */
void FillWallGhosts (amrex::MultiFab& F, amrex::Geometry const& geom, amrex::BCRec const& bcrec,
                     std::string const& reconstruction);

/**
 * \brief The density of one phase for SharpFaceDensity: a value, or a cell-centered field.
 *
 * The field is the density of that phase, not the mixture density: it is defined in every cell,
 * extended across the interface by the solver (a low Mach solver may give the gas density from its
 * temperature and composition everywhere). It has the layout of the F given to SharpFaceDensity and
 * one filled ghost layer, and a face takes the arithmetic mean of the two cells beside it. A phase
 * density that varies along gravity is then balanced exactly away from the interface; on a face cut
 * by the interface rho_face is off by at most |grad rho| h / 8, as the phase's part of the segment
 * between the cell centers has its mean density at its own middle rather than at the face.
 */
struct PhaseDensity
{
    PhaseDensity (amrex::Real a_value) noexcept : value(a_value) {}
    PhaseDensity (amrex::MultiFab const& a_field) noexcept : field(&a_field) {}
    // It keeps a pointer to the field: not to a temporary
    PhaseDensity (amrex::MultiFab&& a_field) = delete;

    amrex::Real value = 0.0;
    amrex::MultiFab const* field = nullptr;
};

/**
 * \brief Face density rho_g + (rho_l - rho_g) * <face_segment_fraction>, with the average over the
 *        given volume fractions, and each phase density a value or a field (PhaseDensity). For a
 *        planar interface normal to gravity it makes rho_face * g a discrete gradient, which is the
 *        hydrostatic balance of a balanced-force projection (Francois et al., JCP 213, 2006). Each
 *        face gets one value. On several AMR levels, call it on each: the MAC projection averages the
 *        coefficients of a fine level down onto the coarse faces it covers.
 *
 * \param rho_face        output, on the faces of the boxes of F with its DistributionMapping; its
 *                        ghost faces, if any, are left as they are
 * \param F               the volume fractions to average, e.g. at t^n and t^{n+1} for t^{n+1/2}, with
 *                        nGrowF(reconstruction) filled ghost cells
 * \param geom            geometry of this level
 * \param rho_l           liquid density, a value or a field
 * \param rho_g           gas density, a value or a field
 * \param reconstruction  "youngs", "myc" or "myc_height"
 */
void SharpFaceDensity (FacePtrs const& rho_face,
                       amrex::Vector<amrex::MultiFab const*> const& F,
                       amrex::Geometry const& geom,
                       PhaseDensity const& rho_l, PhaseDensity const& rho_g,
                       std::string const& reconstruction);

} // namespace PLIC

#endif
/** @}*/

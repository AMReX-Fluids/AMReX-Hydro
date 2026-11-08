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

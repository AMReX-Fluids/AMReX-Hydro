/**
 * \file hydro_plic_normals_K.H
 *
 * Interface normals of the PLIC reconstruction, and the plane fit that takes one of them as a
 * template parameter. The comment above the schemes below says how to add one.
 *
 * The Mixed Youngs-Centered normal is that of Aulisa et al. (JCP 225, 2007). Derived from
 * Kynema-SGF (formerly AMR-Wind), src/equation_systems/vof, BSD 3-Clause License, Copyright (c)
 * 2020, Alliance for Energy Innovation and Regents of the University of California.
 *
 * \addtogroup PLIC
 *  @{
 */

#ifndef HYDRO_PLIC_NORMALS_K_H
#define HYDRO_PLIC_NORMALS_K_H

#include <hydro_plic_K.H>

#include <AMReX.H>

#include <cmath>
#include <exception>
#include <limits>
#include <string>
#include <utility>

namespace PLIC {

// Half-length of the height-function columns (7 cells)
static constexpr int hf_half = 3;
// A height column ends in cells within this of 0 or 1
#ifdef AMREX_USE_FLOAT
static constexpr amrex::Real hf_pure_tol = 1.0e-3_rt;
#else
static constexpr amrex::Real hf_pure_tol = 1.0e-6_rt;
#endif
// Relative difference below which two values tie, to round-off (see tied)
static constexpr amrex::Real tie_tol = 1000.0_rt * std::numeric_limits<amrex::Real>::epsilon();

//! Whether a ties with the largest of its kind amax (to round-off)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool tied(const amrex::Real a, const amrex::Real amax) noexcept {
    return a >= (1.0_rt - tie_tol) * amax;
}

/**
 * \brief The MYC normal from its centered candidates m[c] (height direction c, c < n, height
 *        component +-1) and the Youngs normal m[n] (at any scale): the centered candidate whose
 *        height component dominates after L1 normalization, or the mean of those that tie, so that
 *        no direction is preferred, unless Youngs' largest component is smaller (an interface not
 *        aligned with a grid plane).
 */
template <int n>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void myc_choose(amrex::Real (&m)[n + 1][3],
                                                         amrex::Real& mx,
                                                         amrex::Real& my,
                                                         amrex::Real& mz) {
    // The L1 norm of each candidate, and the components compared after normalization: only the
    // candidates kept are normalized in full
    amrex::Real t[n + 1];
    amrex::Real hc[n];
    amrex::Real t0 = 0.0_rt;
    for (int c = 0; c <= n; ++c) {
        t[c] = std::abs(m[c][0]) + std::abs(m[c][1]) + std::abs(m[c][2]);
        if (c < n) {
            hc[c] = std::abs(m[c][c]) / t[c];
            t0 = amrex::max(t0, hc[c]);
        }
    }
    const amrex::Real ty =
        (t[n] > 0.0_rt)
            ? amrex::max(std::abs(m[n][0]), std::abs(m[n][1]), std::abs(m[n][2])) / t[n]
            : 0.0_rt;
    amrex::Real mc[3] = {0.0_rt, 0.0_rt, 0.0_rt};
    if (ty > 0.0_rt && t0 > ty) {
        for (int d = 0; d < 3; ++d) {
            mc[d] = m[n][d] / t[n];
        }
    } else {
        for (int c = 0; c < n; ++c) {
            if (tied(hc[c], t0)) {
                for (int d = 0; d < 3; ++d) {
                    mc[d] += m[c][d] / t[c];
                }
            }
        }
    }
    // Tied candidates may cancel: a null normal, which fit_plane replaces
    const amrex::Real tc = std::abs(mc[0]) + std::abs(mc[1]) + std::abs(mc[2]);
    const amrex::Real scale = (tc > 0.0_rt) ? 1.0_rt / tc : 0.0_rt;
    mx = mc[0] * scale;
    my = mc[1] * scale;
    mz = mc[2] * scale;
}

#if AMREX_SPACEDIM == 3
/**
 * \brief Youngs' finite-difference normal on the 27-point stencil, scaled
 *        by -32 (i.e. m = -32 grad F).
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
youngs_normal(const int i,
              const int j,
              const int k,
              amrex::Array4<amrex::Real const> const& F,
              amrex::Real& mx,
              amrex::Real& my,
              amrex::Real& mz) {
    amrex::Real mm1 = F(i - 1, j - 1, k - 1) + F(i - 1, j - 1, k + 1) + F(i - 1, j + 1, k - 1) +
                      F(i - 1, j + 1, k + 1) +
                      2.0_rt * (F(i - 1, j - 1, k) + F(i - 1, j + 1, k) + F(i - 1, j, k - 1) +
                                F(i - 1, j, k + 1)) +
                      4.0_rt * F(i - 1, j, k);
    amrex::Real mm2 = F(i + 1, j - 1, k - 1) + F(i + 1, j - 1, k + 1) + F(i + 1, j + 1, k - 1) +
                      F(i + 1, j + 1, k + 1) +
                      2.0_rt * (F(i + 1, j - 1, k) + F(i + 1, j + 1, k) + F(i + 1, j, k - 1) +
                                F(i + 1, j, k + 1)) +
                      4.0_rt * F(i + 1, j, k);
    mx = mm1 - mm2;

    mm1 = F(i - 1, j - 1, k - 1) + F(i - 1, j - 1, k + 1) + F(i + 1, j - 1, k - 1) +
          F(i + 1, j - 1, k + 1) +
          2.0_rt *
              (F(i - 1, j - 1, k) + F(i + 1, j - 1, k) + F(i, j - 1, k - 1) + F(i, j - 1, k + 1)) +
          4.0_rt * F(i, j - 1, k);
    mm2 = F(i - 1, j + 1, k - 1) + F(i - 1, j + 1, k + 1) + F(i + 1, j + 1, k - 1) +
          F(i + 1, j + 1, k + 1) +
          2.0_rt *
              (F(i - 1, j + 1, k) + F(i + 1, j + 1, k) + F(i, j + 1, k - 1) + F(i, j + 1, k + 1)) +
          4.0_rt * F(i, j + 1, k);
    my = mm1 - mm2;

    mm1 = F(i - 1, j - 1, k - 1) + F(i - 1, j + 1, k - 1) + F(i + 1, j - 1, k - 1) +
          F(i + 1, j + 1, k - 1) +
          2.0_rt *
              (F(i - 1, j, k - 1) + F(i + 1, j, k - 1) + F(i, j - 1, k - 1) + F(i, j + 1, k - 1)) +
          4.0_rt * F(i, j, k - 1);
    mm2 = F(i - 1, j - 1, k + 1) + F(i - 1, j + 1, k + 1) + F(i + 1, j - 1, k + 1) +
          F(i + 1, j + 1, k + 1) +
          2.0_rt *
              (F(i - 1, j, k + 1) + F(i + 1, j, k + 1) + F(i, j - 1, k + 1) + F(i, j + 1, k + 1)) +
          4.0_rt * F(i, j, k + 1);
    mz = mm1 - mm2;
}

/**
 * \brief Mixed Youngs-Centered normal (3D).  Returns an L1-normalized
 *        normal pointing out of the liquid.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void myc_normal(const int i,
                                                         const int j,
                                                         const int k,
                                                         amrex::Array4<amrex::Real const> const& F,
                                                         amrex::Real& mx,
                                                         amrex::Real& my,
                                                         amrex::Real& mz) {
    // m[c][d]: candidate c (0-2 centered height-function schemes along x, y,
    // z, 3 Youngs), component d
    amrex::Real m[4][3];

    // Centered scheme with x as the height direction:
    //   sgn(mx) X = my Y + mz Z + alpha
    amrex::Real m1 = F(i - 1, j, k - 1) + F(i - 1, j, k + 1) + F(i - 1, j - 1, k) +
                     F(i - 1, j + 1, k) + F(i - 1, j, k);
    amrex::Real m2 = F(i + 1, j, k - 1) + F(i + 1, j, k + 1) + F(i + 1, j - 1, k) +
                     F(i + 1, j + 1, k) + F(i + 1, j, k);
    m[0][0] = (m1 > m2) ? 1.0_rt : -1.0_rt;
    m1 = F(i - 1, j - 1, k) + F(i + 1, j - 1, k) + F(i, j - 1, k);
    m2 = F(i - 1, j + 1, k) + F(i + 1, j + 1, k) + F(i, j + 1, k);
    m[0][1] = 0.5_rt * (m1 - m2);
    m1 = F(i - 1, j, k - 1) + F(i + 1, j, k - 1) + F(i, j, k - 1);
    m2 = F(i - 1, j, k + 1) + F(i + 1, j, k + 1) + F(i, j, k + 1);
    m[0][2] = 0.5_rt * (m1 - m2);

    // y as the height direction
    m1 = F(i - 1, j - 1, k) + F(i - 1, j + 1, k) + F(i - 1, j, k);
    m2 = F(i + 1, j - 1, k) + F(i + 1, j + 1, k) + F(i + 1, j, k);
    m[1][0] = 0.5_rt * (m1 - m2);
    m1 = F(i, j - 1, k - 1) + F(i, j - 1, k + 1) + F(i + 1, j - 1, k) + F(i - 1, j - 1, k) +
         F(i, j - 1, k);
    m2 = F(i, j + 1, k - 1) + F(i, j + 1, k + 1) + F(i + 1, j + 1, k) + F(i - 1, j + 1, k) +
         F(i, j + 1, k);
    m[1][1] = (m1 > m2) ? 1.0_rt : -1.0_rt;
    m1 = F(i, j - 1, k - 1) + F(i, j, k - 1) + F(i, j + 1, k - 1);
    m2 = F(i, j - 1, k + 1) + F(i, j, k + 1) + F(i, j + 1, k + 1);
    m[1][2] = 0.5_rt * (m1 - m2);

    // z as the height direction
    m1 = F(i - 1, j, k - 1) + F(i - 1, j, k + 1) + F(i - 1, j, k);
    m2 = F(i + 1, j, k - 1) + F(i + 1, j, k + 1) + F(i + 1, j, k);
    m[2][0] = 0.5_rt * (m1 - m2);
    m1 = F(i, j - 1, k - 1) + F(i, j - 1, k + 1) + F(i, j - 1, k);
    m2 = F(i, j + 1, k - 1) + F(i, j + 1, k + 1) + F(i, j + 1, k);
    m[2][1] = 0.5_rt * (m1 - m2);
    m1 = F(i - 1, j, k - 1) + F(i + 1, j, k - 1) + F(i, j - 1, k - 1) + F(i, j + 1, k - 1) +
         F(i, j, k - 1);
    m2 = F(i - 1, j, k + 1) + F(i + 1, j, k + 1) + F(i, j - 1, k + 1) + F(i, j + 1, k + 1) +
         F(i, j, k + 1);
    m[2][2] = (m1 > m2) ? 1.0_rt : -1.0_rt;

    youngs_normal(i, j, k, F, m[3][0], m[3][1], m[3][2]);
    myc_choose<3>(m, mx, my, mz);
}
#else
/**
 * \brief Youngs' finite-difference normal on the 9-point stencil, scaled
 *        by -8 (i.e. m = -8 grad F).
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
youngs_normal(const int i,
              const int j,
              const int k,
              amrex::Array4<amrex::Real const> const& F,
              amrex::Real& mx,
              amrex::Real& my,
              amrex::Real& mz) {
    mx = (F(i - 1, j - 1, k) + 2.0_rt * F(i - 1, j, k) + F(i - 1, j + 1, k)) -
         (F(i + 1, j - 1, k) + 2.0_rt * F(i + 1, j, k) + F(i + 1, j + 1, k));
    my = (F(i - 1, j - 1, k) + 2.0_rt * F(i, j - 1, k) + F(i + 1, j - 1, k)) -
         (F(i - 1, j + 1, k) + 2.0_rt * F(i, j + 1, k) + F(i + 1, j + 1, k));
    mz = 0.0_rt;
}

/**
 * \brief Mixed Youngs-Centered normal (2D).  Returns an L1-normalized
 *        normal pointing out of the liquid.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void myc_normal(const int i,
                                                         const int j,
                                                         const int k,
                                                         amrex::Array4<amrex::Real const> const& F,
                                                         amrex::Real& mx,
                                                         amrex::Real& my,
                                                         amrex::Real& mz) {
    // m[c][d]: candidate c (0-1 centered height-function schemes along x, y,
    // 2 Youngs), component d
    amrex::Real m[3][3] = {{0.0_rt, 0.0_rt, 0.0_rt}, {0.0_rt, 0.0_rt, 0.0_rt}, {0.0_rt, 0.0_rt, 0.0_rt}};

    // Column / row sums
    const amrex::Real cL = F(i - 1, j - 1, k) + F(i - 1, j, k) + F(i - 1, j + 1, k);
    const amrex::Real cR = F(i + 1, j - 1, k) + F(i + 1, j, k) + F(i + 1, j + 1, k);
    const amrex::Real rB = F(i - 1, j - 1, k) + F(i, j - 1, k) + F(i + 1, j - 1, k);
    const amrex::Real rT = F(i - 1, j + 1, k) + F(i, j + 1, k) + F(i + 1, j + 1, k);

    // x as the height direction: sgn(mx) X = my Y + alpha
    m[0][0] = (cL > cR) ? 1.0_rt : -1.0_rt;
    m[0][1] = 0.5_rt * (rB - rT);
    // y as the height direction: sgn(my) Y = mx X + alpha
    m[1][0] = 0.5_rt * (cL - cR);
    m[1][1] = (rB > rT) ? 1.0_rt : -1.0_rt;

    youngs_normal(i, j, k, F, m[2][0], m[2][1], m[2][2]);
    myc_choose<2>(m, mx, my, mz);
}
#endif

//! +1 if the liquid is on the low side of the column through iv along ed, else -1
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
height_sigma(const amrex::IntVect& iv,
             const amrex::IntVect& ed,
             amrex::Array4<amrex::Real const> const& F) {
    return (F(iv - hf_half * ed) > F(iv + hf_half * ed)) ? 1.0_rt : -1.0_rt;
}

//! The column through c along ed over the rows r - hf_half to r + hf_half: false unless F is
//! monotone along it, from pure liquid to pure gas in the orientation sigma. H is the sum of F.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
height_column(const amrex::IntVect& c,
              const amrex::IntVect& ed,
              const amrex::Real sigma,
              const int r,
              amrex::Array4<amrex::Real const> const& F,
              amrex::Real& H) {
    const auto liquid = [&](const int s) {
        const amrex::Real f = F(c + s * ed);
        return (sigma > 0.0_rt) ? f : 1.0_rt - f;
    };
    if (liquid(r - hf_half) < 1.0_rt - hf_pure_tol || liquid(r + hf_half) > hf_pure_tol) {
        return false;
    }
    H = 0.0_rt;
    amrex::Real prev = 1.0_rt;
    for (int s = r - hf_half; s <= r + hf_half; ++s) {
        const amrex::Real l = liquid(s);
        if (l > prev + hf_pure_tol) {
            return false;
        }
        prev = l;
        H += F(c + s * ed);
    }
    return true;
}

/**
 * \brief The row r, along d from (i,j,k), on which height_normal centers the neighboring columns:
 *        where the column through (i,j,k) crosses the interface. The crossing lies Lc - hf_half -
 *        1/2 from the center of (i,j,k), with Lc the sum of the phase on the low side of the
 *        column, and within 1.5 as the cell is mixed: r is the nearest row, in [-1, 1].
 * \return false if the column through (i,j,k) is not a height; sigma as height_sigma
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
height_row(const amrex::IntVect& iv,
           const int d,
           amrex::Array4<amrex::Real const> const& F,
           int& r,
           amrex::Real& sigma) {
    const amrex::IntVect ed = amrex::IntVect::TheDimensionVector(d);
    sigma = height_sigma(iv, ed, F);
    amrex::Real Hc = 0.0_rt;
    if (!height_column(iv, ed, sigma, 0, F, Hc)) {
        return false;
    }
    const amrex::Real Lc = (sigma > 0.0_rt) ? Hc : static_cast<amrex::Real>(2 * hf_half + 1) - Hc;
    r = amrex::max(-1, amrex::min(1, static_cast<int>(
                                 std::floor(Lc - static_cast<amrex::Real>(hf_half)))));
    return true;
}

/**
 * \brief Height-function normal with height direction d: the liquid heights of the 7-cell columns
 *        through the tangential neighbors of (i,j,k), differenced centrally. The columns are
 *        centered on the row where the column through (i,j,k) crosses the interface, at most one
 *        row from it, so they read up to hf_half + 1 cells along d. Exact for a planar interface
 *        whose slopes relative to d are at most 1, which holds for the dominant normal direction:
 *        the crossing of a neighboring column then lies within 2.5 cells of that row, and its
 *        end cells are pure. Returns false unless each column, the center one
 *        included, crosses the interface once: from pure liquid to pure gas in the orientation of
 *        the center column, with F monotone in between. A thin gas gap with a liquid sliver inside
 *        a column is not a height, and the caller falls back to another normal.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
height_normal(const int i,
              const int j,
              const int k,
              const int d,
              amrex::Array4<amrex::Real const> const& F,
              amrex::Real m[3]) {
    amrex::ignore_unused(k);
    const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
    const amrex::IntVect ed = amrex::IntVect::TheDimensionVector(d);
    // sigma: +1 for liquid on the low side
    int r = 0;
    amrex::Real sigma = 0.0_rt;
    if (!height_row(iv, d, F, r, sigma)) {
        return false;
    }
    m[0] = m[1] = m[2] = 0.0_rt;
    for (int t = 0; t < AMREX_SPACEDIM; ++t) {
        if (t == d) {
            continue;
        }
        const amrex::IntVect et = amrex::IntVect::TheDimensionVector(t);
        amrex::Real Hm = 0.0_rt;
        amrex::Real Hp = 0.0_rt;
        if (!height_column(iv - et, ed, sigma, r, F, Hm) ||
            !height_column(iv + et, ed, sigma, r, F, Hp)) {
            return false;
        }
        // m = grad(x_d - h) (sigma = 1) or grad(h - x_d): m_t = (H- - H+) / 2
        m[t] = 0.5_rt * (Hm - Hp);
    }
    m[d] = sigma;
    return true;
}

//-----------------------------------------------------------------------------
// Normal schemes.
//
// To add one, write a struct like these and give it a name in dispatch_normal below. Nothing else
// in the PLIC directory depends on which scheme is in use. A scheme has
//
//   static constexpr int stencil;
//     The largest distance, in cells along any direction, from (i,j,k) that normal() reads.  It
//     sets the ghost cells needed in F (PLIC::nGrowF).
//
//   static void normal(i, j, k, F, mx, my, mz);
//     The normal in cell (i,j,k), pointing out of the liquid, in unit-cube cell coordinates, at any
//     scale, but not far below 1, as a sum|m| below machine epsilon counts as null.  It is called
//     in mixed cells only, and not on structures one cell thick, where
//     fit_plane takes the principal direction of the Hessian of F instead (thin_structure,
//     degenerate_normal), as it does for a null normal.
//
//   static int reach(i, j, k, dir, sgn, F);   (only with stencil > 1)
//     The distance along direction dir, toward side sgn, that normal() reads in this cell, at most
//     stencil. It reads only cells within that distance, so that it holds whatever lies beyond.
//     The wall ghost cells take the plane of a cell near the wall that does not read ghost cells,
//     found with this (see wall_contact_planes).
//-----------------------------------------------------------------------------

/** \brief Youngs' finite-difference normal. */
struct Youngs {
    static constexpr int stencil = 1;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static void
    normal(const int i,
           const int j,
           const int k,
           amrex::Array4<amrex::Real const> const& F,
           amrex::Real& mx,
           amrex::Real& my,
           amrex::Real& mz) noexcept {
        youngs_normal(i, j, k, F, mx, my, mz);
    }
};

/** \brief Mixed Youngs-Centered normal (Aulisa et al., JCP 225, 2007). */
struct MixedYoungsCentered {
    static constexpr int stencil = 1;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static void
    normal(const int i,
           const int j,
           const int k,
           amrex::Array4<amrex::Real const> const& F,
           amrex::Real& mx,
           amrex::Real& my,
           amrex::Real& mz) noexcept {
        myc_normal(i, j, k, F, mx, my, mz);
    }
};

/**
 * \brief The height-function normal along the dominant direction of the Mixed Youngs-Centered
 *        normal, or the latter where the heights are not defined (thin or strongly curved
 *        structures). Exact for a planar interface. Where components tie, the mean over the tied
 *        directions, so that no direction is preferred (the slopes are then 1 relative to each,
 *        where the heights are exact), and the MYC normal unless the heights hold along each.
 */
struct HeightFunctionMYC {
    static constexpr int stencil = hf_half + 1;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static void
    normal(const int i,
           const int j,
           const int k,
           amrex::Array4<amrex::Real const> const& F,
           amrex::Real& mx,
           amrex::Real& my,
           amrex::Real& mz) noexcept {
        myc_normal(i, j, k, F, mx, my, mz);
        const amrex::Real mm[3] = {mx, my, mz};
        if (amrex::max(std::abs(mx), std::abs(my), std::abs(mz)) <= 0.0_rt) {
            return;
        }
        amrex::Real m[3] = {0.0_rt, 0.0_rt, 0.0_rt};
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            if (!height_direction(mm, d)) {
                continue;
            }
            amrex::Real mh[3];
            if (!height_normal(i, j, k, d, F, mh)) {
                return;
            }
            const amrex::Real norm = std::abs(mh[0]) + std::abs(mh[1]) + std::abs(mh[2]);
            for (int a = 0; a < 3; ++a) {
                m[a] += mh[a] / norm;
            }
        }
        mx = m[0];
        my = m[1];
        mz = m[2];
    }

    //! Whether normal() takes heights along d, from the MYC normal mm: its component along d ties
    //! with the largest
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static bool
    height_direction(const amrex::Real mm[3], const int d) noexcept {
        return tied(std::abs(mm[d]), amrex::max(std::abs(mm[0]), std::abs(mm[1]), std::abs(mm[2])));
    }

    // Height columns along dir reach hf_half cells, one more toward side sgn when they are centered
    // a row that way, and the other stencils one
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static int
    reach(const int i,
          const int j,
          const int k,
          const int dir,
          const int sgn,
          amrex::Array4<amrex::Real const> const& F) noexcept {
        amrex::Real mm[3];
        myc_normal(i, j, k, F, mm[0], mm[1], mm[2]);
        if (!height_direction(mm, dir)) {
            return 1;
        }
        int r = 0;
        amrex::Real sigma = 0.0_rt;
        if (!height_row(amrex::IntVect(AMREX_D_DECL(i, j, k)), dir, F, r, sigma)) {
            return hf_half;
        }
        return hf_half + amrex::max(0, r * sgn);
    }
};

// See wall_nsearch
static_assert(wall_nsearch >= 2 * HeightFunctionMYC::stencil + 1,
              "wall_nsearch must cover the stencil of HeightFunctionMYC: see its comment");

/**
 * \brief Whether cell (i,j,k) lies on a structure one cell thick (a film, a gap, a filament or a
 *        drop), where the gradient is noise and the normal comes from the Hessian instead
 *        (degenerate_normal): along some line through it (an axis or a diagonal) its neighbors on
 *        both sides are pure in the same phase, and it differs from them by more than 8 times as
 *        much as F changes across it along any line. The margin keeps the gradient at the tip of
 *        an interface, as on spheres down to a radius of 0.6 cells.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool thin_structure(const int i,
                                                             const int j,
                                                             const int k,
                                                             amrex::Array4<amrex::Real const> const& F) {
    amrex::ignore_unused(k);
    const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
    const amrex::Real fc = F(iv);
    amrex::Real depth = 0.0_rt;
    amrex::Real slope = 0.0_rt;
    // The lines through the cell: the offsets e in {-1, 0, 1}^dim after the center in lexicographic
    // order, one of each pair +-e
    constexpr int noff = AMREX_D_TERM(3, *3, *3);
    for (int o = noff / 2 + 1; o < noff; ++o) {
        const amrex::IntVect e(AMREX_D_DECL(o % 3 - 1, (o / 3) % 3 - 1, o / 9 - 1));
        const amrex::Real fp = F(iv + e);
        const amrex::Real fm = F(iv - e);
        // Pure as at the ends of a height column
        const bool gas = fp <= hf_pure_tol && fm <= hf_pure_tol;
        const bool liquid = fp >= 1.0_rt - hf_pure_tol && fm >= 1.0_rt - hf_pure_tol;
        if (gas || liquid) {
            depth = amrex::max(depth, amrex::min(std::abs(fp - fc), std::abs(fm - fc)));
        }
        slope = amrex::max(slope, 0.5_rt * std::abs(fp - fm));
    }
    return depth > 8.0_rt * slope;
}

/**
 * \brief The orientation of a structure one cell thick (a film, a gap, a filament or a drop), where
 *        the normals vanish or are noise: across it, along the eigenvector of the
 *        largest eigenvalue (in magnitude) of the Hessian of F, from central differences. Where
 *        eigenvalues tie, the sum of their eigenvectors; for a structure aligned with the grid the
 *        Hessian is diagonal, and that is the sum of the tied axes. Each eigenvector is
 *        oriented by a fixed reference direction, so that the cells along a structure agree.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void degenerate_normal(const int i,
                                                                const int j,
                                                                const int k,
                                                                amrex::Array4<amrex::Real const> const& F,
                                                                amrex::Real& mx,
                                                                amrex::Real& my,
                                                                amrex::Real& mz) {
    amrex::ignore_unused(k);
    const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
    constexpr int n = AMREX_SPACEDIM;
    // The Hessian h, and the eigenvectors v (columns) of its Jacobi diagonalization
    amrex::Real h[3][3] = {{0.0_rt, 0.0_rt, 0.0_rt}, {0.0_rt, 0.0_rt, 0.0_rt}, {0.0_rt, 0.0_rt, 0.0_rt}};
    amrex::Real v[3][3] = {{1.0_rt, 0.0_rt, 0.0_rt}, {0.0_rt, 1.0_rt, 0.0_rt}, {0.0_rt, 0.0_rt, 1.0_rt}};
    for (int a = 0; a < n; ++a) {
        const amrex::IntVect ea = amrex::IntVect::TheDimensionVector(a);
        h[a][a] = F(iv + ea) + F(iv - ea) - 2.0_rt * F(iv);
        for (int b = a + 1; b < n; ++b) {
            const amrex::IntVect eb = amrex::IntVect::TheDimensionVector(b);
            h[a][b] = 0.25_rt * (F(iv + ea + eb) - F(iv + ea - eb) - F(iv - ea + eb) + F(iv - ea - eb));
            h[b][a] = h[a][b];
        }
    }
    constexpr amrex::Real eps = std::numeric_limits<amrex::Real>::epsilon();
    bool rotated = true;
    for (int sweep = 0; sweep < 10 && rotated; ++sweep) {
        rotated = false;
        for (int p = 0; p < n; ++p) {
            for (int q = p + 1; q < n; ++q) {
                if (std::abs(h[p][q]) <= eps * (std::abs(h[p][p]) + std::abs(h[q][q])) ||
                    std::abs(h[p][q]) <= std::numeric_limits<amrex::Real>::min()) {
                    continue;
                }
                rotated = true;
                // The rotation in the (p, q) plane that zeroes h[p][q]
                const amrex::Real theta = (h[q][q] - h[p][p]) / (2.0_rt * h[p][q]);
                const amrex::Real t = ((theta >= 0.0_rt) ? 1.0_rt : -1.0_rt) /
                                      (std::abs(theta) + std::sqrt(theta * theta + 1.0_rt));
                const amrex::Real c = 1.0_rt / std::sqrt(t * t + 1.0_rt);
                const amrex::Real sn = t * c;
                // Over all three rows: in 2D the third is zero (of v, the identity row) and stays so
                for (auto& row : h) {
                    const amrex::Real hp = row[p];
                    const amrex::Real hq = row[q];
                    row[p] = c * hp - sn * hq;
                    row[q] = sn * hp + c * hq;
                }
                for (int r = 0; r < n; ++r) {
                    const amrex::Real hp = h[p][r];
                    const amrex::Real hq = h[q][r];
                    h[p][r] = c * hp - sn * hq;
                    h[q][r] = sn * hp + c * hq;
                }
                for (auto& row : v) {
                    const amrex::Real vp = row[p];
                    const amrex::Real vq = row[q];
                    row[p] = c * vp - sn * vq;
                    row[q] = sn * vp + c * vq;
                }
                h[p][q] = 0.0_rt;
                h[q][p] = 0.0_rt;
            }
        }
    }
    amrex::Real lmax = 0.0_rt;
    for (int a = 0; a < n; ++a) {
        lmax = amrex::max(lmax, std::abs(h[a][a]));
    }
    // Each tied eigenvector on the positive side of a reference direction, so that the sum does not
    // cancel and the cells along a structure agree. Its components (1, sqrt(2) - 1, sqrt(3) - 3/2)
    // are rationally independent, so no direction of integer components is normal to it.
    constexpr amrex::Real ref[3] = {1.0_rt, 0.41421356237309515_rt, 0.2320508075688772_rt};
    amrex::Real m[3] = {0.0_rt, 0.0_rt, 0.0_rt};
    for (int a = 0; a < n; ++a) {
        if (!tied(std::abs(h[a][a]), lmax)) {
            continue;
        }
        amrex::Real proj = 0.0_rt;
        for (int r = 0; r < n; ++r) {
            proj += ref[r] * v[r][a];
        }
        const amrex::Real sgn = (proj >= 0.0_rt) ? 1.0_rt : -1.0_rt;
        for (int r = 0; r < n; ++r) {
            m[r] += sgn * v[r][a];
        }
    }
    mx = m[0];
    my = m[1];
    mz = m[2];
}

/**
 * \brief Normalize the normal (mx, my, mz) to sum|m| = 1 and give the constant alpha of the plane
 *        m.x = alpha, in unit-cube cell coordinates, that holds the volume fraction f of the cell.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void plane_constant(const amrex::Real f,
                                                             amrex::Real& mx,
                                                             amrex::Real& my,
                                                             amrex::Real& mz,
                                                             amrex::Real& alpha) {
    // A null normal (fit_plane replaces one by degenerate_normal): pick an arbitrary one
    amrex::Real norm = std::abs(mx) + std::abs(my) + std::abs(mz);
    if (norm <= std::numeric_limits<amrex::Real>::epsilon()) {
        mx = 1.0_rt;
        my = 0.0_rt;
        mz = 0.0_rt;
        norm = 1.0_rt;
    }
    mx /= norm;
    my /= norm;
    mz /= norm;

    alpha = volume_intercept(std::abs(mx), std::abs(my), std::abs(mz), f);

    // Back from the reflected (all positive) frame to the cell frame
    alpha += amrex::min(0.0_rt, mx) + amrex::min(0.0_rt, my) + amrex::min(0.0_rt, mz);
}

/**
 * \brief Reconstruct the PLIC plane m.x = alpha in cell (i,j,k), in unit-cube cell coordinates with
 *        sum|m| = 1. Reads up to Normal::stencil cells away. Only meaningful in mixed cells
 *        (is_mixed).
 *
 * \tparam Normal the normal scheme
 */
template <typename Normal>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void fit_plane(const int i,
                                                        const int j,
                                                        const int k,
                                                        amrex::Array4<amrex::Real const> const& F,
                                                        amrex::Real& mx,
                                                        amrex::Real& my,
                                                        amrex::Real& mz,
                                                        amrex::Real& alpha) {
    // On a structure one cell thick the normal would be null or noise
    bool degenerate = thin_structure(i, j, k, F);
    if (!degenerate) {
        Normal::normal(i, j, k, F, mx, my, mz);
        degenerate = std::abs(mx) + std::abs(my) + std::abs(mz) <= std::numeric_limits<amrex::Real>::epsilon();
    }
    if (degenerate) {
        degenerate_normal(i, j, k, F, mx, my, mz);
    }
    plane_constant(F(i, j, k), mx, my, mz, alpha);
}

/**
 * \brief Normal of an interface that runs inside the wall layer, from a one-sided height function:
 *        with the cell inside the wall layer pure, and the ghost cell beyond it taken as pure of the
 *        other phase, the volume fraction of that phase in a wall cell is the height of the
 *        interface at the cell center, measured from the wall face, and centered differences of
 *        those heights along the wall give its slopes. Exact for a planar interface that stays
 *        within the wall layer across the cell and its tangential neighbors, which is checked: the
 *        layer inside must be pure, of one phase, and each height plus or minus half the change of
 *        height across a cell must lie in [0, 1]. Reads the wall layer and the layer inside it,
 *        one cell away along the wall, and no ghost cells.
 *
 * \param c   the wall cell
 * \param sgn -1 for a low wall, +1 for a high wall
 * \return false where the height function does not hold
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
wall_layer_normal(const amrex::IntVect& c,
                  const int dir,
                  const int sgn,
                  amrex::Array4<amrex::Real const> const& F,
                  amrex::Real* m) {
    const amrex::IntVect e_in = (-sgn) * amrex::IntVect::TheDimensionVector(dir);
    const amrex::Real f_in = F(c + e_in);
    if (is_mixed(f_in)) {
        return false;
    }
    const bool liquid_in = f_in >= 0.5_rt;
    // The tangential directions in increasing order, with a dummy one in 2D
    const int t[2] = {(dir == 0) ? 1 : 0, (AMREX_SPACEDIM == 3) ? ((dir == 2) ? 1 : 2) : 0};
    const int b_max = (AMREX_SPACEDIM == 3) ? 1 : 0;
    // Heights of the wall-side phase, h[b+1][a+1] at offsets a along t[0] and b along t[1]
    amrex::Real h[3][3] = {{0.0_rt, 0.0_rt, 0.0_rt}, {0.0_rt, 0.0_rt, 0.0_rt}, {0.0_rt, 0.0_rt, 0.0_rt}};
    for (int b = -b_max; b <= b_max; ++b) {
        for (int a = -1; a <= 1; ++a) {
            amrex::IntVect cc = c;
            cc[t[0]] += a;
            if (AMREX_SPACEDIM == 3) {
                cc[t[1]] += b;
            }
            const amrex::Real f = F(cc);
            const amrex::Real p = F(cc + e_in);
            if (!is_mixed(f) || is_mixed(p) || ((p >= 0.5_rt) != liquid_in)) {
                return false;
            }
            h[b + 1][a + 1] = liquid_in ? 1.0_rt - f : f;
        }
    }
    const amrex::Real s0 = 0.5_rt * (h[1][2] - h[1][0]);
    const amrex::Real s1 = (AMREX_SPACEDIM == 3) ? 0.5_rt * (h[2][1] - h[0][1]) : 0.0_rt;
    const amrex::Real spread = 0.5_rt * (std::abs(s0) + std::abs(s1));
    for (int b = -b_max; b <= b_max; ++b) {
        for (int a = -1; a <= 1; ++a) {
            const amrex::Real hc = h[b + 1][a + 1];
            if (hc - spread < 0.0_rt || hc + spread > 1.0_rt) {
                return false;
            }
        }
    }
    // m.x - alpha < 0 in the liquid: the liquid is inside the interface if liquid_in
    const amrex::Real sigma = liquid_in ? 1.0_rt : -1.0_rt;
    m[0] = m[1] = m[2] = 0.0_rt;
    m[t[0]] = sigma * s0;
    if (AMREX_SPACEDIM == 3) {
        m[t[1]] = sigma * s1;
    }
    m[dir] = sigma * static_cast<amrex::Real>(sgn);
    return true;
}

[[noreturn]] inline void unknown_reconstruction(std::string const& reconstruction) {
    amrex::Abort("Unknown reconstruction: " + reconstruction +
                 " (expected youngs, myc or myc_height)");
    std::terminate();
}

/**
 * \brief Call fn with the normal scheme of a reconstruction name, as fn(Normal{}), and return what it
 *        returns. The one place that names the schemes; it aborts on an unknown name. Each path
 *        calls fn once, so it is forwarded there.
 */
template <typename Fn>
decltype(auto) dispatch_normal(std::string const& reconstruction, Fn&& fn) {
    if (reconstruction == "youngs") {
        return std::forward<Fn>(fn)(Youngs{});
    }
    if (reconstruction == "myc") {
        return std::forward<Fn>(fn)(MixedYoungsCentered{});
    }
    if (reconstruction != "myc_height") {
        unknown_reconstruction(reconstruction);
    }
    return std::forward<Fn>(fn)(HeightFunctionMYC{});
}

//! Call f(offset, distance) for the offsets along the wall (zero along dir) within wall_nsearch cells
//! in each tangential direction, with distance the sum of their absolute components, by increasing
//! distance. Stops when f returns false.
template <typename Fn>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void for_each_wall_offset(const int dir, Fn const& f) {
    constexpr int ns = wall_nsearch;
#if AMREX_SPACEDIM == 3
    // The tangential directions
    const int to = (dir == 2) ? 1 : 2;
    const int ti = (dir == 0) ? 1 : 0;
    for (int t = 0; t <= 2 * ns; ++t) {
        for (int b = -ns; b <= ns; ++b) {
            const int rem = t - std::abs(b);
            if (rem < 0 || rem > ns) {
                continue;
            }
            // -rem, then +rem unless it is the same
            for (int a = -rem; a <= rem; a += amrex::max(1, 2 * rem)) {
                amrex::IntVect off(0);
                off[ti] = a;
                off[to] = b;
                if (!f(off, t)) {
                    return;
                }
            }
        }
    }
#else
    const int ti = 1 - dir;
    for (int t = 0; t <= ns; ++t) {
        for (int a = -t; a <= t; a += amrex::max(1, 2 * t)) {
            amrex::IntVect off(0);
            off[ti] = a;
            if (!f(off, t)) {
                return;
            }
        }
    }
#endif
}

//! What a cell offers the search along a wall (see wall_candidate)
constexpr int wc_none = 0;    //!< nothing
constexpr int wc_myc = 1;     //!< below the wall layer, mixed: its MYC plane reads only data and no ghost cells
constexpr int wc_plane = 2;   //!< as wc_myc, and the plane of the normal scheme does too
constexpr int wc_height = 3;  //!< in the wall layer, where the height function of wall_layer_normal holds

/**
 * \brief What cell c, up to Normal::stencil layers inside the wall, offers the search of
 *        wall_contact_planes, and the plane it offers in m and alpha. Computed once per cell, as the
 *        searches of neighboring wall cells overlap.
 *
 * \param fit_near nonzero on the cells of F whose 3x3(x3) block holds only data
 * \param fit_ok   nonzero on the cells of F whose reconstruction stencil reads only data
 */
template <typename Normal>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int
wall_candidate(const amrex::IntVect& c,
               const int dir,
               const int iwall,
               const int sgn,
               amrex::Array4<int const> const& fit_near,
               amrex::Array4<int const> const& fit_ok,
               amrex::Array4<amrex::Real const> const& F,
               amrex::Real* m,
               amrex::Real& alpha) {
    if (!fit_near.contains(c) || fit_near(c) == 0) {
        return wc_none;
    }
    const int layer = (iwall - c[dir]) * sgn;
    if (layer == 0) {
        if (!wall_layer_normal(c, dir, sgn, F, m)) {
            return wc_none;
        }
        plane_constant(F(c), m[0], m[1], m[2], alpha);
        return wc_height;
    }
    if (!is_mixed(F(c))) {
        return wc_none;
    }
    // The stencil of the scheme would read cells without data, or wall-normal ghost cells
    const auto cell = c.dim3();
    if constexpr (Normal::stencil > 1) {
        if (fit_ok(c) == 0 ||
            (layer < Normal::stencil && Normal::reach(cell.x, cell.y, cell.z, dir, sgn, F) > layer)) {
            fit_plane<MixedYoungsCentered>(cell.x, cell.y, cell.z, F, m[0], m[1], m[2], alpha);
            return wc_myc;
        }
    }
    fit_plane<Normal>(cell.x, cell.y, cell.z, F, m[0], m[1], m[2], alpha);
    return wc_plane;
}

/**
 * \brief The PLIC planes that the ghost cells beyond the wall cell iv_w extend, when the interface
 *        meets the wall within wall_nsearch cells (tangentially) of iv_w: a free contact line,
 *        exact for a planar interface and for one uniform normal to the wall. The planes are those
 *        of the nearest interface cells, up to Normal::stencil layers inside the wall, whose
 *        reconstruction does not read ghost cells (Normal::reach). Below about 26.6 degrees from
 *        the wall, or near another wall, where myc_height finds none, the MYC planes of the nearest
 *        interface cells below the wall layer, which are exact there (see wall_nsearch), down to about
 *        atan(1/wall_nsearch). Below that, the interface runs inside the wall layer near the
 *        contact line, and the planes are those of the nearest wall cells where the one-sided
 *        height function of wall_layer_normal holds. All the equally near cells give their plane,
 *        so that no direction is preferred. No plane reads ghost cells, so one pass of the fill is
 *        enough.
 *
 * \tparam Normal the normal scheme
 * \param iwall index along dir of the wall-adjacent interior cell
 * \param sgn   -1 for a low wall, +1 for a high wall
 * \param cand  wall_candidate of the cells searched
 * \param pl    the plane each candidate offers (m_x, m_y, m_z, alpha)
 * \param m     the sum of the planes taken, in the coordinates of iv_w, with alpha
 * \return the number of planes
 */
template <typename Normal>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int
wall_contact_planes(const amrex::IntVect& iv_w,
                    const int dir,
                    const int iwall,
                    const int sgn,
                    amrex::Array4<int const> const& cand,
                    amrex::Array4<amrex::Real const> const& pl,
                    amrex::Real* m,
                    amrex::Real& alpha) {
    constexpr int nfar = 1000;

    // The sum of the planes of the cells at the nearest distance found so far
    struct Sum {
        amrex::Real m[3] = {0.0_rt, 0.0_rt, 0.0_rt};
        amrex::Real alpha = 0.0_rt;
        int n = 0;
        int dist = nfar;
        AMREX_GPU_HOST_DEVICE void add (amrex::IntVect const& iv, amrex::IntVect const& ref,
                                        amrex::Array4<amrex::Real const> const& p, int d_ref) {
            if (d_ref < dist) {
                *this = Sum{};
                dist = d_ref;
            }
            alpha += p(ref, 3);
            for (int d = 0; d < 3; ++d) {
                m[d] += p(ref, d);
            }
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                alpha -= p(ref, d) * static_cast<amrex::Real>(iv[d] - ref[d]);
            }
            ++n;
        }
    };
    const auto take = [&] (Sum const& sum) {
        for (int d = 0; d < 3; ++d) {
            m[d] = sum.m[d];
        }
        alpha = sum.alpha;
        return sum.n;
    };

    // The nearest cells below the wall layer: with a plane of the normal scheme that reads no ghost
    // cells, or else, as a fallback, with any plane, which is then a MYC plane (only for a normal with
    // a wider stencil, as with stencil 1 any such cell already gives a plane of the scheme). The
    // rings are visited in order, so each sum adds its cells in the same order.
    Sum plane;
    Sum any;
    for_each_wall_offset(dir, [&] (amrex::IntVect const& off, int tdist) {
        // The cells searched lie a layer or more inside the wall: none further can be nearer
        if (tdist + 1 > plane.dist) {
            return false;
        }
        for (int layer = 1; layer <= Normal::stencil; ++layer) {
            const int dist = tdist + layer;
            amrex::IntVect c = iv_w + off;
            c[dir] = iwall - layer * sgn;
            const int kind = cand(c);
            // With stencil 1 every such cell has a plane of the scheme
            if constexpr (Normal::stencil > 1) {
                if (kind != wc_none && dist <= any.dist) {
                    any.add(iv_w, c, pl, dist);
                }
            }
            if (kind == wc_plane && dist <= plane.dist) {
                plane.add(iv_w, c, pl, dist);
            }
        }
        return true;
    });
    if (plane.n > 0) {
        return take(plane);
    }
    if (any.n > 0) {
        return take(any);
    }

    // The interface inside the wall layer: the nearest wall cells where its height function holds
    Sum layer0;
    for_each_wall_offset(dir, [&] (amrex::IntVect const& off, int tdist) {
        if (tdist > layer0.dist) {
            return false;
        }
        const amrex::IntVect c = iv_w + off;
        if (cand(c) == wc_height) {
            layer0.add(iv_w, c, pl, tdist);
        }
        return true;
    });
    return take(layer0);
}

} // namespace PLIC

#endif
/** @}*/

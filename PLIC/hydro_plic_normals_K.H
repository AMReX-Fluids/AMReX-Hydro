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
// Relative difference below which two values tie (to round-off): the sums of F differenced by the
// normals, the components of the MYC normal for the height direction, and the height components of
// the centered candidates of MYC
static constexpr amrex::Real hf_tie_tol = 1000.0_rt * std::numeric_limits<amrex::Real>::epsilon();

//! lo - hi for two sums of F, or 0 where they tie: the same cells summed in another order differ by
//! round-off, which would give a structure symmetric across the cell a spurious gradient
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real sum_difference(const amrex::Real lo,
                                                                    const amrex::Real hi) noexcept {
    const amrex::Real d = lo - hi;
    return (std::abs(d) <= hf_tie_tol * (std::abs(lo) + std::abs(hi))) ? 0.0_rt : d;
}

//! The sign of the height component of a centered MYC candidate from the liquid on its low and high
//! sides: 0 where they tie, as then the interface does not cross that direction
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real crossing_sign(const amrex::Real lo,
                                                                   const amrex::Real hi) noexcept {
    const amrex::Real d = sum_difference(lo, hi);
    return (d > 0.0_rt) ? 1.0_rt : ((d < 0.0_rt) ? -1.0_rt : 0.0_rt);
}

//! Whether a ties with the largest of its kind amax (to round-off)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool tied(const amrex::Real a, const amrex::Real amax) noexcept {
    return a >= (1.0_rt - hf_tie_tol) * amax;
}

/**
 * \brief The MYC normal from its centered candidates m[c] (height direction c, c < n) and the Youngs
 *        normal m[n], each at any scale: the centered candidate whose height component dominates
 *        after L1 normalization, or the mean of those that tie, so that no direction is preferred,
 *        unless Youngs' largest component is smaller (an interface not aligned with a grid plane).
 *        A candidate without a crossing (height component 0) is never kept, and without any the
 *        normal is Youngs, or null (see degenerate_normal).
 */
template <int n>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void myc_choose(amrex::Real (&m)[n + 1][3],
                                                         amrex::Real& mx,
                                                         amrex::Real& my,
                                                         amrex::Real& mz) {
    amrex::Real t0 = 0.0_rt;
    for (int c = 0; c <= n; ++c) {
        const amrex::Real t = std::abs(m[c][0]) + std::abs(m[c][1]) + std::abs(m[c][2]);
        if (t > 0.0_rt) {
            for (int d = 0; d < 3; ++d) {
                m[c][d] /= t;
            }
        }
        if (c < n) {
            t0 = amrex::max(t0, std::abs(m[c][c]));
        }
    }
    const amrex::Real ty = amrex::max(std::abs(m[n][0]), std::abs(m[n][1]), std::abs(m[n][2]));
    amrex::Real mc[3] = {0.0_rt, 0.0_rt, 0.0_rt};
    if (ty > 0.0_rt && (t0 > ty || t0 == 0.0_rt)) {
        for (int d = 0; d < 3; ++d) {
            mc[d] = m[n][d];
        }
    } else if (t0 > 0.0_rt) {
        for (int c = 0; c < n; ++c) {
            if (tied(std::abs(m[c][c]), t0)) {
                for (int d = 0; d < 3; ++d) {
                    mc[d] += m[c][d];
                }
            }
        }
    }
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
    mx = sum_difference(mm1, mm2);

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
    my = sum_difference(mm1, mm2);

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
    mz = sum_difference(mm1, mm2);
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
    m[0][0] = crossing_sign(m1, m2);
    m1 = F(i - 1, j - 1, k) + F(i + 1, j - 1, k) + F(i, j - 1, k);
    m2 = F(i - 1, j + 1, k) + F(i + 1, j + 1, k) + F(i, j + 1, k);
    m[0][1] = 0.5_rt * sum_difference(m1, m2);
    m1 = F(i - 1, j, k - 1) + F(i + 1, j, k - 1) + F(i, j, k - 1);
    m2 = F(i - 1, j, k + 1) + F(i + 1, j, k + 1) + F(i, j, k + 1);
    m[0][2] = 0.5_rt * sum_difference(m1, m2);

    // y as the height direction
    m1 = F(i - 1, j - 1, k) + F(i - 1, j + 1, k) + F(i - 1, j, k);
    m2 = F(i + 1, j - 1, k) + F(i + 1, j + 1, k) + F(i + 1, j, k);
    m[1][0] = 0.5_rt * sum_difference(m1, m2);
    m1 = F(i, j - 1, k - 1) + F(i, j - 1, k + 1) + F(i + 1, j - 1, k) + F(i - 1, j - 1, k) +
         F(i, j - 1, k);
    m2 = F(i, j + 1, k - 1) + F(i, j + 1, k + 1) + F(i + 1, j + 1, k) + F(i - 1, j + 1, k) +
         F(i, j + 1, k);
    m[1][1] = crossing_sign(m1, m2);
    m1 = F(i, j - 1, k - 1) + F(i, j, k - 1) + F(i, j + 1, k - 1);
    m2 = F(i, j - 1, k + 1) + F(i, j, k + 1) + F(i, j + 1, k + 1);
    m[1][2] = 0.5_rt * sum_difference(m1, m2);

    // z as the height direction
    m1 = F(i - 1, j, k - 1) + F(i - 1, j, k + 1) + F(i - 1, j, k);
    m2 = F(i + 1, j, k - 1) + F(i + 1, j, k + 1) + F(i + 1, j, k);
    m[2][0] = 0.5_rt * sum_difference(m1, m2);
    m1 = F(i, j - 1, k - 1) + F(i, j - 1, k + 1) + F(i, j - 1, k);
    m2 = F(i, j + 1, k - 1) + F(i, j + 1, k + 1) + F(i, j + 1, k);
    m[2][1] = 0.5_rt * sum_difference(m1, m2);
    m1 = F(i - 1, j, k - 1) + F(i + 1, j, k - 1) + F(i, j - 1, k - 1) + F(i, j + 1, k - 1) +
         F(i, j, k - 1);
    m2 = F(i - 1, j, k + 1) + F(i + 1, j, k + 1) + F(i, j - 1, k + 1) + F(i, j + 1, k + 1) +
         F(i, j, k + 1);
    m[2][2] = crossing_sign(m1, m2);

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
    mx = sum_difference(F(i - 1, j - 1, k) + 2.0_rt * F(i - 1, j, k) + F(i - 1, j + 1, k),
                        F(i + 1, j - 1, k) + 2.0_rt * F(i + 1, j, k) + F(i + 1, j + 1, k));
    my = sum_difference(F(i - 1, j - 1, k) + 2.0_rt * F(i, j - 1, k) + F(i + 1, j - 1, k),
                        F(i - 1, j + 1, k) + 2.0_rt * F(i, j + 1, k) + F(i + 1, j + 1, k));
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
    m[0][0] = crossing_sign(cL, cR);
    m[0][1] = 0.5_rt * sum_difference(rB, rT);
    // y as the height direction: sgn(my) Y = mx X + alpha
    m[1][0] = 0.5_rt * sum_difference(cL, cR);
    m[1][1] = crossing_sign(rB, rT);

    youngs_normal(i, j, k, F, m[2][0], m[2][1], m[2][2]);
    myc_choose<2>(m, mx, my, mz);
}
#endif

/**
 * \brief Height-function normal with height direction d: the liquid heights of the 7-cell columns
 *        through the tangential neighbors of (i,j,k), differenced centrally. The columns are
 *        centered on the row where the column through (i,j,k) crosses the interface, at most one
 *        row from it, so they read up to hf_half + 1 cells along d. Exact for a planar interface
 *        whose slopes relative to d are at most 1, which holds for the dominant normal direction:
 *        the crossing of a neighboring column then lies within 2.5 cells of that row, and its
 *        end cells are pure. Returns false (m untouched) unless each column, the center one
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
    // +1 for liquid on the low side
    const amrex::Real sigma = (F(iv - hf_half * ed) > F(iv + hf_half * ed)) ? 1.0_rt : -1.0_rt;
    // The column through c over the rows r - hf_half to r + hf_half: false unless F is monotone
    // along it, from pure liquid to pure gas. H is the sum of F.
    const auto column = [&](const amrex::IntVect& c, const int r, amrex::Real& H) {
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
    };
    amrex::Real Hc = 0.0_rt;
    if (!column(iv, 0, Hc)) {
        return false;
    }
    // The column crosses Lc - hf_half - 1/2 along d from the center of (i,j,k), with Lc the sum of
    // the phase on its low side, and within 1.5 as the cell is mixed: the nearest row
    const amrex::Real Lc = (sigma > 0.0_rt) ? Hc : static_cast<amrex::Real>(2 * hf_half + 1) - Hc;
    const int r = amrex::max(-1, amrex::min(1, static_cast<int>(
                                     std::floor(Lc - static_cast<amrex::Real>(hf_half)))));
    amrex::Real mt[3] = {0.0_rt, 0.0_rt, 0.0_rt};
    for (int t = 0; t < AMREX_SPACEDIM; ++t) {
        if (t == d) {
            continue;
        }
        const amrex::IntVect et = amrex::IntVect::TheDimensionVector(t);
        amrex::Real Hm = 0.0_rt;
        amrex::Real Hp = 0.0_rt;
        if (!column(iv - et, r, Hm) || !column(iv + et, r, Hp)) {
            return false;
        }
        // m = grad(x_d - h) (sigma = 1) or grad(h - x_d): m_t = (H- - H+) / 2
        mt[t] = 0.5_rt * (Hm - Hp);
    }
    mt[d] = sigma;
    for (int a = 0; a < 3; ++a) {
        m[a] = mt[a];
    }
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
//     scale.  It is called in mixed cells only.  A null normal is replaced by an arbitrary one.
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
        const amrex::Real am[3] = {std::abs(mx), std::abs(my), std::abs(mz)};
        const amrex::Real amax = amrex::max(am[0], am[1], am[2]);
        if (amax <= 0.0_rt) {
            return;
        }
        amrex::Real m[3] = {0.0_rt, 0.0_rt, 0.0_rt};
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            if (!tied(am[d], amax)) {
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
};

/**
 * \brief The orientation of a structure one cell thick and symmetric across it (a film, a gap, a
 *        filament or a drop), where the normals vanish: across it, along the eigenvector of the
 *        largest eigenvalue (in magnitude) of the Hessian of F, from central differences. Where
 *        eigenvalues tie, the sum of their eigenvectors; for a structure aligned with the grid the
 *        Hessian is diagonal, and that is the sum of the tied axes. The sign is arbitrary.
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
    for (int sweep = 0; sweep < 10; ++sweep) {
        for (int p = 0; p < n; ++p) {
            for (int q = p + 1; q < n; ++q) {
                if (std::abs(h[p][q]) <= std::numeric_limits<amrex::Real>::min()) {
                    continue;
                }
                // The rotation in the (p, q) plane that zeroes h[p][q]
                const amrex::Real theta = (h[q][q] - h[p][p]) / (2.0_rt * h[p][q]);
                const amrex::Real t = ((theta >= 0.0_rt) ? 1.0_rt : -1.0_rt) /
                                      (std::abs(theta) + std::sqrt(theta * theta + 1.0_rt));
                const amrex::Real c = 1.0_rt / std::sqrt(t * t + 1.0_rt);
                const amrex::Real sn = t * c;
                for (int r = 0; r < n; ++r) {
                    const amrex::Real hp = h[r][p];
                    const amrex::Real hq = h[r][q];
                    h[r][p] = c * hp - sn * hq;
                    h[r][q] = sn * hp + c * hq;
                }
                for (int r = 0; r < n; ++r) {
                    const amrex::Real hp = h[p][r];
                    const amrex::Real hq = h[q][r];
                    h[p][r] = c * hp - sn * hq;
                    h[q][r] = sn * hp + c * hq;
                }
                for (int r = 0; r < n; ++r) {
                    const amrex::Real vp = v[r][p];
                    const amrex::Real vq = v[r][q];
                    v[r][p] = c * vp - sn * vq;
                    v[r][q] = sn * vp + c * vq;
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
    // Each tied eigenvector with its largest component positive, so that the sum does not cancel
    amrex::Real m[3] = {0.0_rt, 0.0_rt, 0.0_rt};
    for (int a = 0; a < n; ++a) {
        if (!tied(std::abs(h[a][a]), lmax)) {
            continue;
        }
        int big = 0;
        for (int r = 1; r < n; ++r) {
            if (std::abs(v[r][a]) > std::abs(v[big][a])) {
                big = r;
            }
        }
        const amrex::Real sgn = (v[big][a] >= 0.0_rt) ? 1.0_rt : -1.0_rt;
        for (int r = 0; r < n; ++r) {
            m[r] += sgn * v[r][a];
        }
    }
    mx = m[0];
    my = m[1];
    mz = m[2];
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
    Normal::normal(i, j, k, F, mx, my, mz);

    // A structure one cell thick and symmetric across it gives a null normal
    amrex::Real norm = std::abs(mx) + std::abs(my) + std::abs(mz);
    if (norm <= std::numeric_limits<amrex::Real>::epsilon()) {
        degenerate_normal(i, j, k, F, mx, my, mz);
        norm = std::abs(mx) + std::abs(my) + std::abs(mz);
    }
    mx /= norm;
    my /= norm;
    mz /= norm;

    alpha = volume_intercept(std::abs(mx), std::abs(my), std::abs(mz), F(i, j, k));

    // Back from the reflected (all positive) frame to the cell frame
    alpha += amrex::min(0.0_rt, mx) + amrex::min(0.0_rt, my) + amrex::min(0.0_rt, mz);
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

} // namespace PLIC

#endif
/** @}*/

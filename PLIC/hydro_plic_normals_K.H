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

#include <exception>
#include <string>

namespace PLIC {

// Half-length of the height-function columns (7 cells)
static constexpr int hf_half = 3;
// A height column ends in cells within this of 0 or 1
#ifdef AMREX_USE_FLOAT
static constexpr amrex::Real hf_pure_tol = 1.0e-3_rt;
#else
static constexpr amrex::Real hf_pure_tol = 1.0e-6_rt;
#endif

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

    // L1-normalize the centered candidates and keep the one whose height
    // direction dominates
    for (int c = 0; c < 3; ++c) {
        const amrex::Real t = std::abs(m[c][0]) + std::abs(m[c][1]) + std::abs(m[c][2]);
        for (int d = 0; d < 3; ++d) {
            m[c][d] /= t;
        }
    }
    int cn = 0;
    amrex::Real t0 = std::abs(m[0][0]);
    if (std::abs(m[1][1]) > t0) {
        t0 = std::abs(m[1][1]);
        cn = 1;
    }
    if (std::abs(m[2][2]) > t0) {
        cn = 2;
    }

    // Youngs candidate
    constexpr amrex::Real tiny = std::numeric_limits<amrex::Real>::epsilon();
    youngs_normal(i, j, k, F, m[3][0], m[3][1], m[3][2]);
    amrex::Real t = std::abs(m[3][0]) + std::abs(m[3][1]) + std::abs(m[3][2]) + tiny;
    for (int d = 0; d < 3; ++d) {
        m[3][d] /= t;
    }

    // Youngs wins when its largest component is smaller than the centered
    // scheme's height component (interface not aligned with a grid plane)
    t = amrex::max(std::abs(m[3][0]), std::abs(m[3][1]), std::abs(m[3][2]));
    if (std::abs(m[cn][cn]) > t && t > 0.0_rt) {
        cn = 3;
    }

    mx = m[cn][0];
    my = m[cn][1];
    mz = m[cn][2];
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
    amrex::Real m[3][2];

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

    for (int c = 0; c < 2; ++c) {
        const amrex::Real t = std::abs(m[c][0]) + std::abs(m[c][1]);
        m[c][0] /= t;
        m[c][1] /= t;
    }
    int cn = (std::abs(m[1][1]) > std::abs(m[0][0])) ? 1 : 0;

    constexpr amrex::Real tiny = std::numeric_limits<amrex::Real>::epsilon();
    amrex::Real mzy = 0.0_rt;
    youngs_normal(i, j, k, F, m[2][0], m[2][1], mzy);
    amrex::Real t = std::abs(m[2][0]) + std::abs(m[2][1]) + tiny;
    m[2][0] /= t;
    m[2][1] /= t;

    t = amrex::max(std::abs(m[2][0]), std::abs(m[2][1]));
    if (std::abs(m[cn][cn]) > t && t > 0.0_rt) {
        cn = 2;
    }

    mx = m[cn][0];
    my = m[cn][1];
    mz = 0.0_rt;
}
#endif

/**
 * \brief Height-function normal with height direction d: the liquid heights
 *        of the 7-cell columns through (i,j,k) and its neighbors in the
 *        tangential directions, differenced centrally.  Exact for a planar
 *        interface whose slope relative to d is at most 1, which holds for
 *        the dominant normal direction.  Returns false (m untouched) when a
 *        column does not span from pure liquid to pure gas in the same
 *        orientation as the others.
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
    // Orientation from the center column: +1 liquid on the low side
    const amrex::Real f_lo = F(iv - hf_half * ed);
    const amrex::Real f_hi = F(iv + hf_half * ed);
    amrex::Real sigma = 0.0_rt;
    if (f_lo >= 1.0_rt - hf_pure_tol && f_hi <= hf_pure_tol) {
        sigma = 1.0_rt;
    } else if (f_lo <= hf_pure_tol && f_hi >= 1.0_rt - hf_pure_tol) {
        sigma = -1.0_rt;
    } else {
        return false;
    }
    const auto column = [&](const amrex::IntVect& c, amrex::Real& H) {
        const amrex::Real a = F(c - hf_half * ed);
        const amrex::Real b = F(c + hf_half * ed);
        const amrex::Real a_pure = (sigma > 0.0_rt) ? 1.0_rt - a : a;
        const amrex::Real b_pure = (sigma > 0.0_rt) ? b : 1.0_rt - b;
        if (a_pure > hf_pure_tol || b_pure > hf_pure_tol) {
            return false;
        }
        H = 0.0_rt;
        for (int s = -hf_half; s <= hf_half; ++s) {
            H += F(c + s * ed);
        }
        return true;
    };
    amrex::Real mt[3] = {0.0_rt, 0.0_rt, 0.0_rt};
    for (int t = 0; t < AMREX_SPACEDIM; ++t) {
        if (t == d) {
            continue;
        }
        const amrex::IntVect et = amrex::IntVect::TheDimensionVector(t);
        amrex::Real Hm = 0.0_rt;
        amrex::Real Hp = 0.0_rt;
        if (!column(iv - et, Hm) || !column(iv + et, Hp)) {
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
//
//   static int reach(i, j, k, dir, F);
//     The distance along direction dir that normal() reads in this cell, at most stencil.  Near a
//     domain boundary it tells whether the normal of a cell reads ghost cells.
//
// All must be callable on the device.
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

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static int
    reach(const int /*i*/,
          const int /*j*/,
          const int /*k*/,
          const int /*dir*/,
          amrex::Array4<amrex::Real const> const& /*F*/) noexcept {
        return 1;
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

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static int
    reach(const int /*i*/,
          const int /*j*/,
          const int /*k*/,
          const int /*dir*/,
          amrex::Array4<amrex::Real const> const& /*F*/) noexcept {
        return 1;
    }
};

/**
 * \brief The height-function normal along the dominant direction of the Mixed Youngs-Centered
 *        normal, or the latter where the heights are not defined (thin or strongly curved
 *        structures). Exact for a planar interface wherever the heights are defined.
 */
struct HeightFunctionMYC {
    static constexpr int stencil = hf_half;

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
        int d = (am[1] > am[0]) ? 1 : 0;
        if (AMREX_SPACEDIM == 3 && am[2] > am[d]) {
            d = 2;
        }
        amrex::Real mh[3];
        if (height_normal(i, j, k, d, F, mh)) {
            mx = mh[0];
            my = mh[1];
            mz = mh[2];
        }
    }

    // Height columns along dir reach hf_half cells, and the other stencils one
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static int
    reach(const int i,
          const int j,
          const int k,
          const int dir,
          amrex::Array4<amrex::Real const> const& F) noexcept {
        amrex::Real mm[3];
        myc_normal(i, j, k, F, mm[0], mm[1], mm[2]);
        for (int t = 0; t < AMREX_SPACEDIM; ++t) {
            if (t != dir && std::abs(mm[t]) > std::abs(mm[dir])) {
                return 1;
            }
        }
        return hf_half;
    }
};

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

    // An isolated blob can produce a null normal: pick an arbitrary one
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
 *        returns. The one place that names the schemes; it aborts on an unknown name.
 */
template <typename Fn>
decltype(auto) dispatch_normal(std::string const& reconstruction, Fn&& fn) {
    if (reconstruction == "youngs") {
        return fn(Youngs{});
    }
    if (reconstruction == "myc") {
        return fn(MixedYoungsCentered{});
    }
    if (reconstruction != "myc_height") {
        unknown_reconstruction(reconstruction);
    }
    return fn(HeightFunctionMYC{});
}

} // namespace PLIC

#endif
/** @}*/

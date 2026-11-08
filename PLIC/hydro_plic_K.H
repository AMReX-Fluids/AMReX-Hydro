/**
 * \file hydro_plic_K.H
 *
 * Geometric kernels of the piecewise-linear interface calculation (PLIC).
 *
 * The interface in a cell is the plane m.x = alpha in unit-cube cell coordinates, with the liquid
 * on the side m.x <= alpha (m points out of the liquid). The volume <-> plane-constant relations
 * are the analytic ones of Scardovelli & Zaleski (JCP 164, 2000).
 *
 * Derived from Kynema-SGF (formerly AMR-Wind),
 * src/equation_systems/vof, BSD 3-Clause License, Copyright (c) 2020, Alliance for Energy
 * Innovation and Regents of the University of California. Changes: 2D variants, per-face flux
 * evaluation, AMReX-Hydro conventions, normal schemes as template parameters.
 *
 * \addtogroup PLIC
 *  @{
 */

#ifndef HYDRO_PLIC_K_H
#define HYDRO_PLIC_K_H

#include <AMReX_Algorithm.H>
#include <AMReX_Array4.H>
#include <AMReX_Geometry.H>
#include <AMReX_REAL.H>
#include <cmath>
#include <limits>

namespace PLIC {

using namespace amrex::literals;

// Cells within this of 0 or 1 count as pure. It must stay well above the noise in F: otherwise a
// cell flips between mixed and pure from one neighbour to the next, which breaks the symmetry of
// the flow. The noise is about 1e-11 in double precision.
#ifdef AMREX_USE_FLOAT
static constexpr amrex::Real pure_tol = 1.0e-4_rt;
#else
static constexpr amrex::Real pure_tol = 1.0e-8_rt;
#endif

/** \brief Whether a cell with volume fraction f holds the interface. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool is_mixed(const amrex::Real f) noexcept {
    return f > pure_tol && f < 1.0_rt - pure_tol;
}

/** \brief a, b and c in ascending order, as lo <= mid <= hi. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void sort3(const amrex::Real a,
                                                    const amrex::Real b,
                                                    const amrex::Real c,
                                                    amrex::Real& lo,
                                                    amrex::Real& mid,
                                                    amrex::Real& hi) noexcept {
    lo = amrex::min(a, b);
    hi = amrex::max(a, b);
    mid = c;
    if (mid < lo) {
        amrex::Swap(lo, mid);
    } else if (mid > hi) {
        amrex::Swap(mid, hi);
    }
}

/**
 * \brief The root of the cubic of volume_intercept in its two upper volume ranges, in the
 *        trigonometric form of Scardovelli & Zaleski, without the constant offset of the range.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real trig_root(const amrex::Real p,
                                                               const amrex::Real q) noexcept {
    const amrex::Real p12 = std::sqrt(amrex::max(0.0_rt, p));
    const amrex::Real pp = p * p12;
    const amrex::Real teta =
        std::acos((pp > 0.0_rt) ? amrex::Clamp(q / pp, -1.0_rt, 1.0_rt) : 1.0_rt) / 3.0_rt;
    const amrex::Real cs = std::cos(teta);
    return p12 * (std::sqrt(3.0_rt * (1.0_rt - cs * cs)) - cs);
}

// Cells searched along a wall for the interface meeting it, on either side of a wall ghost cell
// (see wall_ghost_fraction)
static constexpr int wall_nsearch = 8;

/**
 * \brief Analytic plane constant alpha such that the region
 *        m1*x1 + m2*x2 + m3*x3 <= alpha of the unit cube has volume volF.
 *        Requires m1,m2,m3 >= 0 and m1+m2+m3 = 1.  Any one (2D) or two
 *        (grid-aligned) coefficients may be zero.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real volume_intercept(const amrex::Real b1,
                                                                      const amrex::Real b2,
                                                                      const amrex::Real b3,
                                                                      const amrex::Real volF) {
    constexpr amrex::Real tiny = std::numeric_limits<amrex::Real>::epsilon();

    amrex::Real m1, m2, m3;
    sort3(b1, b2, b3, m1, m2, m3);

    // Volume at the breakpoints V1 < V2 < V3
    const amrex::Real m12 = m1 + m2;
    const amrex::Real pr = 6.0_rt * m1 * m2 * m3;
    // Without a corner region (m1 = 0, a 2D or grid-aligned plane), V1 = 0
    const amrex::Real V1 = (pr > 0.0_rt) ? m1 * m1 * m1 / pr : 0.0_rt;
    const amrex::Real V2 = V1 + 0.5_rt * (m2 - m1) / (m3 + tiny);
    amrex::Real mm = 0.0_rt;
    amrex::Real V3 = 0.0_rt;
    if (m3 < m12) {
        mm = m3;
        V3 = (m3 * m3 * (3.0_rt * m12 - m3) + m1 * m1 * (m1 - 3.0_rt * m3) +
              m2 * m2 * (m2 - 3.0_rt * m3)) /
             pr;
    } else {
        mm = m12;
        V3 = 0.5_rt * mm / (m3 + tiny);
    }

    // Use symmetry to work with ch <= 1/2
    const amrex::Real ch = amrex::min(volF, 1.0_rt - volF);

    amrex::Real alpha = 0.0_rt;
    if (ch < V1) {
        alpha = std::cbrt(pr * ch);
    } else if (ch < V2) {
        alpha = 0.5_rt * (m1 + std::sqrt(m1 * m1 + 8.0_rt * m2 * m3 * (ch - V1)));
    } else if (ch < V3) {
        const amrex::Real p = 2.0_rt * m1 * m2;
        const amrex::Real q = 1.5_rt * m1 * m2 * (m12 - 2.0_rt * m3 * ch);
        alpha = trig_root(p, q) + m12;
    } else if (m12 < m3) {
        alpha = m3 * ch + 0.5_rt * mm;
    } else {
        const amrex::Real p = m1 * (m2 + m3) + m2 * m3 - 0.25_rt;
        const amrex::Real q = 1.5_rt * m1 * m2 * m3 * (0.5_rt - ch);
        alpha = trig_root(p, q) + 0.5_rt;
    }

    if (volF > 0.5_rt) {
        alpha = 1.0_rt - alpha;
    }
    return alpha;
}

/**
 * \brief Liquid volume fraction of the sub-box [r0, r0+dr0] x [0,1] x [0,1]
 *        of the unit cube cut by m1*x1 + m2*x2 + m3*x3 <= alpha.  The
 *        result is normalized by the sub-box volume dr0.  Coefficients may
 *        have any sign; m3 = 0 gives the 2D result.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real cut_volume(const amrex::Real m1,
                                                                const amrex::Real m2,
                                                                const amrex::Real m3,
                                                                const amrex::Real alpha,
                                                                const amrex::Real r0,
                                                                const amrex::Real dr0) {
    constexpr amrex::Real tiny = std::numeric_limits<amrex::Real>::epsilon();

    // Move origin to r0, then reflect so that all coefficients are positive
    amrex::Real al = alpha - m1 * r0;
    al += amrex::max(0.0_rt, -m1 * dr0) + amrex::max(0.0_rt, -m2) + amrex::max(0.0_rt, -m3);

    // Normalize to the unit cube: n1*y1 + n2*y2 + n3*y3 = al, sum(n) = 1
    const amrex::Real tmp = amrex::max(tiny, std::abs(m1) * dr0 + std::abs(m2) + std::abs(m3));
    const amrex::Real n1 = std::abs(m1) / tmp;
    const amrex::Real n2 = std::abs(m2) / tmp;
    const amrex::Real n3 = std::abs(m3) / tmp;
    al = amrex::max(0.0_rt, amrex::min(1.0_rt, al / tmp));
    const amrex::Real al0 = amrex::min(al, 1.0_rt - al);

    amrex::Real b1, b2, b3;
    sort3(n1 * dr0, n2, n3, b1, b2, b3);
    const amrex::Real b12 = b1 + b2;

    // Scardovelli & Zaleski (JCP 2000), as arranged by Kawano (C&F 2016)
    const amrex::Real a = al0;
    amrex::Real v = 0.0_rt;
    if (a > 0.0_rt) {
        if (a < b1) {
            v = a * a * a / (6.0_rt * b1 * b2 * b3);
        } else if (a < b2) {
            v = a * (a - b1) / (2.0_rt * b2 * b3) + b1 * b1 / (6.0_rt * b2 * b3 + tiny);
        } else if (a < amrex::min(b12, b3)) {
            v = (a * a * (3.0_rt * b12 - a) + b1 * b1 * (b1 - 3.0_rt * a) +
                 b2 * b2 * (b2 - 3.0_rt * a)) /
                (6.0_rt * b1 * b2 * b3);
        } else if (b3 < b12) {
            v = (a * a * (3.0_rt - 2.0_rt * a) + b1 * b1 * (b1 - 3.0_rt * a) +
                 b2 * b2 * (b2 - 3.0_rt * a) + b3 * b3 * (b3 - 3.0_rt * a)) /
                (6.0_rt * b1 * b2 * b3);
        } else {
            v = (a - 0.5_rt * b12) / b3;
        }
    }

    return (al <= 0.5_rt) ? v : 1.0_rt - v;
}

/**
 * \brief Liquid volume fraction of cell (i,j,k) for an implicit surface
 *        phi(x) > 0 inside the liquid, by midpoint sampling on nsub^dim
 *        sub-cells.  Intended for problem initialization.
 */
template <typename Phi>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
cell_fraction(const int i,
              const int j,
              const int k,
              amrex::GeometryData const& geomdata,
              Phi const& phi,
              const int nsub = 8) {
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* dx = geomdata.CellSize();
    const amrex::Real h = 1.0_rt / static_cast<amrex::Real>(nsub);
    int inside = 0;
#if AMREX_SPACEDIM == 3
    for (int kk = 0; kk < nsub; ++kk) {
        const amrex::Real z =
            prob_lo[2] +
            (static_cast<amrex::Real>(k) + (static_cast<amrex::Real>(kk) + 0.5_rt) * h) * dx[2];
#else
    amrex::ignore_unused(k);
    {
        const amrex::Real z = 0.0_rt;
#endif
        for (int jj = 0; jj < nsub; ++jj) {
            const amrex::Real y =
                prob_lo[1] +
                (static_cast<amrex::Real>(j) + (static_cast<amrex::Real>(jj) + 0.5_rt) * h) * dx[1];
            for (int ii = 0; ii < nsub; ++ii) {
                const amrex::Real x = prob_lo[0] + (static_cast<amrex::Real>(i) +
                                                    (static_cast<amrex::Real>(ii) + 0.5_rt) * h) *
                                                       dx[0];
                if (phi(x, y, z) > 0.0_rt) {
                    ++inside;
                }
            }
        }
    }
    return static_cast<amrex::Real>(inside) /
           static_cast<amrex::Real>(AMREX_D_TERM(nsub, *nsub, *nsub));
}

/**
 * \brief True if cell (i,j,k) is within nband cells of the interface.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
interface_band(const int i,
               const int j,
               const int k,
               amrex::Array4<amrex::Real const> const& F,
               const int nband = 1) {
    amrex::Real Fmax = 0.0_rt;
    amrex::Real Fmin = 1.0_rt;
#if AMREX_SPACEDIM == 3
    for (int kk = -nband; kk <= nband; ++kk) {
#else
    {
        constexpr int kk = 0;
#endif
        for (int jj = -nband; jj <= nband; ++jj) {
            for (int ii = -nband; ii <= nband; ++ii) {
                const amrex::Real f = F(i + ii, j + jj, k + kk);
                if (is_mixed(f)) {
                    return true;
                }
                Fmax = amrex::max(Fmax, f);
                Fmin = amrex::min(Fmin, f);
            }
        }
    }
    // Interface lying exactly on a cell face
    return (Fmax - Fmin) > 0.5_rt;
}

/**
 * \brief Exact liquid volume fraction of cell (i,j,k) for the half space
 *        n.x <= d (physical coordinates, n need not be normalized).
 *        Intended for problem initialization.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
plane_fraction(const int i,
               const int j,
               const int k,
               amrex::GeometryData const& geomdata,
               const amrex::Real n[AMREX_SPACEDIM],
               const amrex::Real d) {
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* dx = geomdata.CellSize();
    const int idx[3] = {i, j, k};
    amrex::Real m[3] = {0.0_rt, 0.0_rt, 0.0_rt};
    amrex::Real alpha = d;
    amrex::Real norm = 0.0_rt;
    for (int a = 0; a < AMREX_SPACEDIM; ++a) {
        m[a] = n[a] * dx[a];
        alpha -= n[a] * (prob_lo[a] + static_cast<amrex::Real>(idx[a]) * dx[a]);
        norm += std::abs(m[a]);
    }
    if (norm <= std::numeric_limits<amrex::Real>::min()) {
        return (alpha >= 0.0_rt) ? 1.0_rt : 0.0_rt;
    }
    return amrex::Clamp(
        cut_volume(m[0] / norm, m[1] / norm, m[2] / norm, alpha / norm, 0.0_rt, 1.0_rt), 0.0_rt,
        1.0_rt);
}

/**
 * \brief Volume fraction of the ghost cell g beyond the wall cell iv_w, from the plane (m, alpha) of
 *        cell ref (see wall_contact_plane). The ghost takes the change of the plane's volume from the
 *        wall cell: the same as the plane's volume at g for a planar interface, and the value of the
 *        wall cell where the plane does not vary across the wall, whichever cell the plane comes
 *        from, so that a flow uniform normal to the wall stays uniform.
 *
 * \param F_w volume fraction of the wall cell
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real wall_ghost_fraction(const amrex::IntVect& g,
                                                                         const amrex::IntVect& iv_w,
                                                                         const amrex::IntVect& ref,
                                                                         const amrex::Real* m,
                                                                         const amrex::Real alpha,
                                                                         const amrex::Real F_w) {
    // Plane constant in the unit-cube coordinates of g and of the wall cell
    amrex::Real alpha_g = alpha;
    amrex::Real alpha_w = alpha;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        alpha_g -= m[d] * static_cast<amrex::Real>(g[d] - ref[d]);
        alpha_w -= m[d] * static_cast<amrex::Real>(iv_w[d] - ref[d]);
    }
    const amrex::Real V_g = cut_volume(m[0], m[1], m[2], alpha_g, 0.0_rt, 1.0_rt);
    const amrex::Real V_w = cut_volume(m[0], m[1], m[2], alpha_w, 0.0_rt, 1.0_rt);
    return amrex::Clamp(F_w + V_g - V_w, 0.0_rt, 1.0_rt);
}

/**
 * \brief Liquid fraction of the half segment from the center of cell (i,j,k) to the center of its
 *        face in direction dir on side sgn (+1 high, -1 low), from the plane of the cell.
 *
 * \param plane the PLIC plane of each mixed cell (m_x, m_y, m_z, alpha)
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
half_segment_fraction(const int i,
                      const int j,
                      const int k,
                      const int dir,
                      const int sgn,
                      amrex::Array4<amrex::Real const> const& F,
                      amrex::Array4<amrex::Real const> const& plane) {
    const amrex::Real Fc = F(i, j, k);
    if (Fc <= pure_tol) {
        return 0.0_rt;
    }
    if (Fc >= 1.0_rt - pure_tol) {
        return 1.0_rt;
    }
    const amrex::Real m[3] = {plane(i, j, k, 0), plane(i, j, k, 1), plane(i, j, k, 2)};
    const amrex::Real alpha = plane(i, j, k, 3);
    // psi = m.x - alpha < 0 in the liquid, linear along the half segment
    const amrex::Real psi_c = 0.5_rt * (m[0] + m[1] + m[2]) - alpha;
    const amrex::Real psi_f = psi_c + 0.5_rt * static_cast<amrex::Real>(sgn) * m[dir];
    if (psi_c <= 0.0_rt && psi_f <= 0.0_rt) {
        return 1.0_rt;
    }
    if (psi_c > 0.0_rt && psi_f > 0.0_rt) {
        return 0.0_rt;
    }
    const amrex::Real t = psi_c / (psi_c - psi_f);
    return (psi_c <= 0.0_rt) ? t : 1.0_rt - t;
}

/**
 * \brief Liquid fraction of the segment joining the centers of the cells on either side of the low
 *        face of (i,j,k) in direction dir. Exact for a planar interface, which makes
 *        rho_face * g a discrete gradient (hydrostatic balance) when
 *        rho_face = rho_g + (rho_l - rho_g) * this.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
face_segment_fraction(const int i,
                      const int j,
                      const int k,
                      const int dir,
                      amrex::Array4<amrex::Real const> const& F,
                      amrex::Array4<amrex::Real const> const& plane) {
    const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
    const amrex::IntVect ivm = iv - amrex::IntVect::TheDimensionVector(dir);
#if AMREX_SPACEDIM == 2
    const int km = 0;
#else
    const int km = ivm[2];
#endif
    return 0.5_rt * (half_segment_fraction(ivm[0], ivm[1], km, dir, 1, F, plane) +
                     half_segment_fraction(i, j, k, dir, -1, F, plane));
}

/**
 * \brief Mixture property (density, viscosity) of a cell with liquid fraction F: the arithmetic, or
 *        if harmonic the harmonic, volume average of the phase values
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real mixture(const amrex::Real F,
                                                             const amrex::Real liq,
                                                             const amrex::Real gas,
                                                             const bool harmonic = false) noexcept {
    const amrex::Real f = amrex::Clamp(F, 0.0_rt, 1.0_rt);
    if (harmonic) {
        return liq * gas / (f * gas + (1.0_rt - f) * liq);
    }
    return gas + (liq - gas) * f;
}

} // namespace PLIC

#endif
/** @}*/

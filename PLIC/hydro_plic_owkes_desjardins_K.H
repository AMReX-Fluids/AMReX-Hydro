/**
 * \file hydro_plic_owkes_desjardins_K.H
 *
 * Kernels of the unsplit geometric advection of Owkes & Desjardins (JCP 270, 2014).
 *
 * \addtogroup PLIC
 *  @{
 */

#ifndef HYDRO_PLIC_OWKES_DESJARDINS_K_H
#define HYDRO_PLIC_OWKES_DESJARDINS_K_H

#include <hydro_plic_K.H>

namespace PLIC {

//-----------------------------------------------------------------------------
// Unsplit geometric flux kernels for the PLIC volume of fluid, after Owkes &
// Desjardins (JCP 270, 2014).
//
// The liquid crossing a face in a step is the liquid inside the flux volume,
// the region bounded by the face and its corners traced back along the flow.
// Everything here works in cell units: cell (i,j,k) is the unit cube
// [i,i+1] x [j,j+1] x [k,k+1], and face velocities are Courant numbers.  In 2D
// the mesh is one cell thick in z, so that one kernel serves both dimensions.
//
// The flux of a face is computed in coordinates relative to the face (see
// face_liquid_flux), from the displacement of each mesh node over the step.
// Two boxes sharing a face, or a face and its periodic image, then do the same
// arithmetic and get the same flux to the bit, so the volume is conserved
// exactly without synchronizing the faces.
//
// Differences from the paper, none of which change the method:
// - The flux volume is the signed sum of the tetrahedra joining a corner of the
//   face to the triangles of its boundary, rather than the tabulated partition
//   of its Fig. 9.  The two are equivalent: any closed oriented triangulated
//   surface gives the right signed integral, including crossed flux volumes.
// - The side faces are shared by the flux volumes of neighboring faces, so they
//   must be triangulated identically by each of them.  The diagonal joins the
//   traced lower node of the mesh edge to its upper node, which depends only
//   on the mesh edge.  This is what makes the fluxes tile without overlaps.
// - Simplices are cut using an explicit stack, for GPU use, and the normals
//   come from the reconstruction in use rather than ELVIRA, as the paper
//   allows.
//-----------------------------------------------------------------------------

// Tolerances in cell units: positions on a mesh plane (above the round-off in positions of order
// one), and a traced face with no extent normal to the face.
#ifdef AMREX_USE_FLOAT
static constexpr amrex::Real cut_eps = 1.0e-3_rt;
static constexpr amrex::Real apex_tol = 1.0e-5_rt;
#else
static constexpr amrex::Real cut_eps = 1.0e-12_rt;
static constexpr amrex::Real apex_tol = 1.0e-10_rt;
#endif

struct Vec3 {
    amrex::Real x[3];
};

struct Tet {
    Vec3 v[4];
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Vec3 lerp(const Vec3& a,
                                                   const Vec3& b,
                                                   const amrex::Real t) noexcept {
    return Vec3{{a.x[0] + t * (b.x[0] - a.x[0]), a.x[1] + t * (b.x[1] - a.x[1]),
                 a.x[2] + t * (b.x[2] - a.x[2])}};
}

/**
 * \brief Signed volume of the tetrahedron with apex o over the triangle
 *        (a,b,c): positive when the triangle's right-handed normal points away
 *        from o.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
signed_volume(const Vec3& o, const Vec3& a, const Vec3& b, const Vec3& c) noexcept {
    const amrex::Real ax = a.x[0] - o.x[0], ay = a.x[1] - o.x[1], az = a.x[2] - o.x[2];
    const amrex::Real bx = b.x[0] - o.x[0], by = b.x[1] - o.x[1], bz = b.x[2] - o.x[2];
    const amrex::Real cx = c.x[0] - o.x[0], cy = c.x[1] - o.x[1], cz = c.x[2] - o.x[2];
    return (ax * (by * cz - bz * cy) + ay * (bz * cx - bx * cz) + az * (bx * cy - by * cx)) /
           6.0_rt;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real abs_volume(const Tet& t) noexcept {
    return std::abs(signed_volume(t.v[0], t.v[1], t.v[2], t.v[3]));
}

/** \brief Tetrahedra filling the prism with triangles p, q and edges p-q */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void prism_to_tets(const Vec3& p0,
                                                            const Vec3& p1,
                                                            const Vec3& p2,
                                                            const Vec3& q0,
                                                            const Vec3& q1,
                                                            const Vec3& q2,
                                                            Tet* out) noexcept {
    out[0] = Tet{{p0, p1, p2, q2}};
    out[1] = Tet{{p0, p1, q2, q1}};
    out[2] = Tet{{p0, q0, q1, q2}};
}

/**
 * \brief Split a tetrahedron by the plane where the vertex values d vanish.
 *
 * The pieces with d <= 0 are written to out first, then those with d > 0.
 *
 * \param out         room for 6 tetrahedra
 * \param n_neg,n_pos number of pieces with d <= 0 and with d > 0, at most 3 each
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
split_tet(const Tet& t, const amrex::Real* d, Tet* out, int& n_neg, int& n_pos) noexcept {
    int vin[4];
    int vout[4];
    int ni = 0;
    int no = 0;
    for (int v = 0; v < 4; ++v) {
        if (d[v] <= 0.0_rt) {
            vin[ni++] = v;
        } else {
            vout[no++] = v;
        }
    }
    n_neg = 0;
    n_pos = 0;
    if (no == 0) {
        out[0] = t;
        n_neg = 1;
        return;
    }
    if (ni == 0) {
        out[0] = t;
        n_pos = 1;
        return;
    }

    // Point where the edge joining an inside and an outside vertex crosses
    auto cross = [&](const int a, const int b) {
        return lerp(t.v[a], t.v[b], d[a] / (d[a] - d[b]));
    };
    if (ni == 1) {
        const int a = vin[0];
        const Vec3 p0 = cross(a, vout[0]);
        const Vec3 p1 = cross(a, vout[1]);
        const Vec3 p2 = cross(a, vout[2]);
        out[0] = Tet{{t.v[a], p0, p1, p2}};
        n_neg = 1;
        prism_to_tets(t.v[vout[0]], t.v[vout[1]], t.v[vout[2]], p0, p1, p2, out + 1);
        n_pos = 3;
    } else if (no == 1) {
        const int a = vout[0];
        const Vec3 p0 = cross(vin[0], a);
        const Vec3 p1 = cross(vin[1], a);
        const Vec3 p2 = cross(vin[2], a);
        prism_to_tets(t.v[vin[0]], t.v[vin[1]], t.v[vin[2]], p0, p1, p2, out);
        n_neg = 3;
        out[3] = Tet{{t.v[a], p0, p1, p2}};
        n_pos = 1;
    } else {
        const int a = vin[0];
        const int b = vin[1];
        const int c = vout[0];
        const int e = vout[1];
        const Vec3 ac = cross(a, c);
        const Vec3 ae = cross(a, e);
        const Vec3 bc = cross(b, c);
        const Vec3 be = cross(b, e);
        prism_to_tets(t.v[a], ac, ae, t.v[b], bc, be, out);
        n_neg = 3;
        prism_to_tets(t.v[c], ac, bc, t.v[e], ae, be, out + 3);
        n_pos = 3;
    }
}

/**
 * \brief Volume of the part of a tetrahedron of volume vol on the liquid side of the plane
 *        m.(x - x_lo) <= alpha of the cell with low corner x_lo.
 *
 * With one vertex on the liquid side, that part is the corner of the tetrahedron at the vertex,
 * cut at a fraction r of each edge from it: vol times the product of the r.  With one vertex on the
 * gas side, it is the rest of the tetrahedron.  With two on each side, it is the prism of split_tet
 * on the liquid side.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real tet_liquid_volume(const Tet& t,
                                                                       const amrex::Real vol,
                                                                       const int i,
                                                                       const int j,
                                                                       const int k,
                                                                       const amrex::Real* plane) noexcept {
    amrex::Real d[4];
    int vin[4];
    int vout[4];
    int ni = 0;
    int no = 0;
    for (int v = 0; v < 4; ++v) {
        d[v] = plane[0] * (t.v[v].x[0] - static_cast<amrex::Real>(i)) +
               plane[1] * (t.v[v].x[1] - static_cast<amrex::Real>(j)) +
               plane[2] * (t.v[v].x[2] - static_cast<amrex::Real>(k)) - plane[3];
        if (d[v] <= 0.0_rt) {
            vin[ni++] = v;
        } else {
            vout[no++] = v;
        }
    }
    if (no == 0) {
        return vol;
    }
    if (ni == 0) {
        return 0.0_rt;
    }
    // Fraction of the edge from vertex a, on one side, to vertex b, on the other, before the plane.
    // d[a] and d[b] have opposite signs, so the fraction is in [0,1].
    auto frac = [&](const int a, const int b) { return d[a] / (d[a] - d[b]); };
    if (ni == 1) {
        const int a = vin[0];
        return vol * (frac(a, vout[0]) * frac(a, vout[1]) * frac(a, vout[2]));
    }
    if (no == 1) {
        const int b = vout[0];
        return vol * (1.0_rt - frac(b, vin[0]) * frac(b, vin[1]) * frac(b, vin[2]));
    }
    const int a = vin[0];
    const int b = vin[1];
    const int c = vout[0];
    const int e = vout[1];
    Tet prism[3];
    prism_to_tets(t.v[a], lerp(t.v[a], t.v[c], frac(a, c)), lerp(t.v[a], t.v[e], frac(a, e)), t.v[b],
                  lerp(t.v[b], t.v[c], frac(b, c)), lerp(t.v[b], t.v[e], frac(b, e)), prism);
    return abs_volume(prism[0]) + abs_volume(prism[1]) + abs_volume(prism[2]);
}

/** \brief Liquid volume of a piece of tetrahedron lying within cell (i,j,k). */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
piece_liquid(const Tet& t,
             const amrex::Real vol,
             const int i,
             const int j,
             const int k,
             amrex::Array4<amrex::Real const> const& F,
             amrex::Array4<amrex::Real const> const& plane) noexcept {
    const amrex::Real f = F(i, j, k);
    if (f <= pure_tol) {
        return 0.0_rt;
    }
    if (f >= 1.0_rt - pure_tol) {
        return vol;
    }
    const amrex::Real pl[4] = {plane(i, j, k, 0), plane(i, j, k, 1), plane(i, j, k, 2),
                               plane(i, j, k, 3)};
    return tet_liquid_volume(t, vol, i, j, k, pl);
}

// Tetrahedra held while cutting by the mesh planes.  Each cut replaces one
// tetrahedron by at most 6 and a depth-first path crosses each mesh plane at
// most once.  For a Courant number below 1 a flux volume lies within the cells
// around its face (no node moves a cell, and the correction apex is held there,
// see flux_volume): it crosses at most 1 mesh plane along the face normal
// and 2 along each other direction, 5 in all, and the stack holds at most
// 1 + 5 * 5 = 26 tetrahedra.  Kept small, as each GPU thread holds its own.
static constexpr int cut_stack_size = 32;

/**
 * \brief Liquid volume of a tetrahedron: cut along the mesh planes into one
 *        piece per cell, each with the liquid of its cell (piece_liquid).
 *
 * \param bx cells where F and plane are valid; pieces beyond it are charged to
 *           the nearest cell inside it
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
cut_tet_liquid(const Tet& tet,
               const amrex::Box& bx,
               amrex::Array4<amrex::Real const> const& F,
               amrex::Array4<amrex::Real const> const& plane) noexcept {
    amrex::Real liquid = 0.0_rt;
    Tet stack[cut_stack_size];
    int top = 0;
    stack[top++] = tet;
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    while (top > 0) {
        const Tet t = stack[--top];

        // The first mesh plane crossing this tetrahedron, if any.  The z planes
        // play no part in 2D, whose mesh is one cell thick.
        int cut_dir = -1;
        amrex::Real cut_at = 0.0_rt;
        for (int d = 0; d < AMREX_SPACEDIM && cut_dir < 0; ++d) {
            amrex::Real tmin = t.v[0].x[d];
            amrex::Real tmax = t.v[0].x[d];
            for (int v = 1; v < 4; ++v) {
                tmin = amrex::min(tmin, t.v[v].x[d]);
                tmax = amrex::max(tmax, t.v[v].x[d]);
            }
            const amrex::Real p = std::floor(tmin + cut_eps) + 1.0_rt;
            if (tmax > p + cut_eps) {
                cut_dir = d;
                cut_at = p;
            }
        }

        // A full stack means the flux volume reaches beyond the cells around its face, which the
        // check on the Courant numbers rules out in every build. Should it happen anyway, a release
        // build charges the uncut piece to the cell of its centroid.
        AMREX_ASSERT_WITH_MESSAGE(cut_dir < 0 || top + 6 <= cut_stack_size,
                                  "PLIC: flux volume needs more cuts than the stack holds");
        if (cut_dir >= 0 && top + 6 <= cut_stack_size) {
            const amrex::Real d[4] = {t.v[0].x[cut_dir] - cut_at, t.v[1].x[cut_dir] - cut_at,
                                      t.v[2].x[cut_dir] - cut_at, t.v[3].x[cut_dir] - cut_at};
            int n_neg = 0;
            int n_pos = 0;
            split_tet(t, d, stack + top, n_neg, n_pos);
            top += n_neg + n_pos;
            continue;
        }

        // Within a single cell, found from the centroid
        const amrex::Real vol = abs_volume(t);
        if (vol <= 0.0_rt) {
            continue;
        }
        int c[3] = {0, 0, 0};
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            c[d] = static_cast<int>(
                std::floor(0.25_rt * (t.v[0].x[d] + t.v[1].x[d] + t.v[2].x[d] + t.v[3].x[d])));
        }
        const int ci = amrex::Clamp(c[0], lo.x, hi.x);
        const int cj = amrex::Clamp(c[1], lo.y, hi.y);
        const int ck = amrex::Clamp(c[2], lo.z, hi.z);
        liquid += piece_liquid(t, vol, ci, cj, ck, F, plane);
    }
    return liquid;
}

/**
 * \brief Apex o of the two tetrahedra (o,a,c,e) and (o,a,b,c) added on a
 *        traced face (a,b,c,e) to give them a total signed volume vcor.  o lies
 *        over the barycenter of the face; the volume is linear in o along dir.
 *
 * \param tri the triangles (a,c,e) and (a,b,c) of the traced face, with +dir normals
 * \return false if the face has no extent normal to dir (no correction)
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool correction_apex(const Vec3* const (&tri)[2][3],
                                                              const int dir,
                                                              const amrex::Real vcor,
                                                              Vec3& o) noexcept {
    const Vec3& a = *tri[0][0];
    const Vec3& b = *tri[1][1];
    const Vec3& c = *tri[0][1];
    const Vec3& e = *tri[0][2];
    for (int d = 0; d < 3; ++d) {
        o.x[d] = 0.25_rt * (a.x[d] + b.x[d] + c.x[d] + e.x[d]);
    }
    amrex::Real s0 = 0.0_rt; // 6 V at o_dir = 0
    amrex::Real sd = 0.0_rt; // d(6 V)/d(o_dir)
    for (const auto& n : tri) {
        const Vec3 &p = *n[0], &q = *n[1], &r = *n[2];
        const amrex::Real u[3] = {q.x[0] - p.x[0], q.x[1] - p.x[1], q.x[2] - p.x[2]};
        const amrex::Real w[3] = {r.x[0] - p.x[0], r.x[1] - p.x[1], r.x[2] - p.x[2]};
        const amrex::Real nrm[3] = {u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2],
                                    u[0] * w[1] - u[1] * w[0]};
        // V = (p - o).N / 6
        for (int d = 0; d < 3; ++d) {
            const amrex::Real od = (d == dir) ? 0.0_rt : o.x[d];
            s0 += (p.x[d] - od) * nrm[d];
        }
        sd -= nrm[dir];
    }
    if (std::abs(sd) <= apex_tol) {
        return false;
    }
    o.x[dir] = (6.0_rt * vcor - s0) / sd;
    return true;
}

// Tetrahedra of a flux volume: the fan and the two of the correction
static constexpr int max_flux_tets = 10;

/**
 * \brief A flux volume as signed tetrahedra over 9 points: the face corners a, b, c, e
 *        counterclockwise about +dir (points 0-3), their traced positions a', b', c', e' (4-7),
 *        and the apex o of the correction (8).
 */
struct FluxVolume {
    Vec3 p[9];
    //! Signed volume of each tetrahedron
    amrex::Real vol[max_flux_tets];
    //! Number of tetrahedra: the 8 of the fan, and the 2 of the correction if there is one
    int ntet;
    //! Whether the volume of the tetrahedra is the Courant number: it falls short where the apex
    //! is held near the face, or there is no correction
    bool matches_courant;

    //! Number of points in use: o only with a correction
    [[nodiscard]] AMREX_GPU_HOST_DEVICE int npoint() const noexcept { return (ntet > 8) ? 9 : 8; }

    /**
     * \brief Tetrahedron n: the fan of the boundary from a, with outward triangles, then the
     *        correction from o over the traced face.  Triangles touching a vanish: the face
     *        itself and the side faces through a.  The diagonals of the side faces depend only on
     *        the mesh edge (see above).
     */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE Tet tet(const int n) const noexcept {
        constexpr int v[max_flux_tets][4] = {{0, 4, 5, 1},  // side a-b, diagonal a' - b
                                             {0, 3, 7, 4},  // side e-a, diagonal e - a'
                                             {0, 1, 5, 2},  // side b-c, diagonal b' - c
                                             {0, 5, 6, 2},  //
                                             {0, 2, 6, 7},  // side c-e, diagonal e' - c
                                             {0, 2, 7, 3},  //
                                             {0, 4, 7, 6},  // traced face, diagonal a' - c'
                                             {0, 4, 6, 5},  //
                                             {8, 4, 6, 7},  // correction, +dir normals
                                             {8, 4, 5, 6}};
        return Tet{{p[v[n][0]], p[v[n][1]], p[v[n][2]], p[v[n][3]]}};
    }
};

/**
 * \brief The flux volume of the low face of cell (i,j,k) in direction dir.
 *
 * \param disp    displacement of each mesh node over the step, from its traced
 *                position at the start of the step, in cell units (3 components)
 * \param courant the face velocity times dt/dx, which the flux volume is
 *                corrected to match
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
flux_volume(const int i,
            const int j,
            const int k,
            const int dir,
            const amrex::Real courant,
            amrex::Array4<amrex::Real const> const& disp,
            FluxVolume& fv) noexcept {
    // Right-handed tangential directions.  The face spans 0..1 along each, and
    // in 2D the node arrays have no extent in z: the k+1 nodes are the k nodes
    // shifted by one cell.
    const int t1 = (dir + 1) % 3;
    const int t2 = (dir + 2) % 3;

    int off1[3] = {0, 0, 0};
    int off2[3] = {0, 0, 0};
    off1[t1] = 1;
    off2[t2] = 1;

    auto node = [&](const int s1, const int s2) {
        return Vec3{{static_cast<amrex::Real>(i + s1 * off1[0] + s2 * off2[0]),
                     static_cast<amrex::Real>(j + s1 * off1[1] + s2 * off2[1]),
                     static_cast<amrex::Real>(k + s1 * off1[2] + s2 * off2[2])}};
    };
    auto traced = [&](const int s1, const int s2) {
        const Vec3 p = node(s1, s2);
        const int ni = i + s1 * off1[0] + s2 * off2[0];
        const int nj = j + s1 * off1[1] + s2 * off2[1];
        // Nodes of a 2D mesh: z stays as it was, and only the in-plane offset selects the node
        constexpr bool is3d = (AMREX_SPACEDIM == 3);
        const int nk = is3d ? k + s1 * off1[2] + s2 * off2[2] : 0;
        return Vec3{{p.x[0] + disp(ni, nj, nk, 0), p.x[1] + disp(ni, nj, nk, 1),
                     is3d ? p.x[2] + disp(ni, nj, nk, 2) : p.x[2]}};
    };

    // Face corners counterclockwise about +dir, and their traced positions
    Vec3* const p = fv.p;
    p[0] = node(0, 0);
    p[1] = node(1, 0);
    p[2] = node(1, 1);
    p[3] = node(0, 1);
    p[4] = traced(0, 0);
    p[5] = traced(1, 0);
    p[6] = traced(1, 1);
    p[7] = traced(0, 1);

    // Fan of the boundary from a.  With outward triangles the signed volumes
    // add up to the flux volume, positive for a flow along +dir.
    amrex::Real vol = 0.0_rt;
    for (int n = 0; n < 8; ++n) {
        const Tet t = fv.tet(n);
        fv.vol[n] = signed_volume(t.v[0], t.v[1], t.v[2], t.v[3]);
        vol += fv.vol[n];
    }

    // Extra volume so that the flux volume matches the face velocity exactly:
    // two tetrahedra from the traced face to a point o above its barycenter, in
    // the direction of dir (Eq. 27-28).  Their volume is linear in o along dir.
    const Vec3* const ttri[2][3] = {{&p[4], &p[6], &p[7]}, {&p[4], &p[5], &p[6]}}; // +dir normals
    Vec3& o = p[8];
    const bool has_cor = correction_apex(ttri, dir, courant - vol, o);
    // The rest of the flux volume lies within a cell of the face, as no node moves a cell. Where the
    // traced face has almost no extent normal to dir, the apex would lie far beyond, so it is held
    // there too: the flux volume then falls short of the Courant number (see flux_volume_liquid), but
    // every face still reads only the cells around it, and the boxes sharing a face clamp it the same
    // way.
    const amrex::Real o_dir = o.x[dir];
    o.x[dir] = amrex::Clamp(o.x[dir], p[0].x[dir] - 1.0_rt, p[0].x[dir] + 1.0_rt);
    fv.matches_courant = has_cor && o.x[dir] == o_dir;

    fv.ntet = has_cor ? max_flux_tets : 8;
    for (int n = 8; n < fv.ntet; ++n) {
        const Tet t = fv.tet(n);
        fv.vol[n] = signed_volume(t.v[0], t.v[1], t.v[2], t.v[3]);
    }
}

/**
 * \brief Whether the cells [lo, hi], limited to bx, are all gas or all liquid.
 *
 * \return true if they are pure, with a_liquid set to 1 or 0 accordingly
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool pure_cells(const int* lo,
                                                         const int* hi,
                                                         amrex::Array4<amrex::Real const> const& F,
                                                         const amrex::Box& bx,
                                                         bool& a_liquid) noexcept {
    const auto blo = amrex::lbound(bx);
    const auto bhi = amrex::ubound(bx);
    bool all_gas = true;
    bool all_liq = true;
    // In 2D bx is one cell thick, at k = 0
    const int klo = amrex::Clamp(lo[2], blo.z, bhi.z);
    const int khi = amrex::Clamp(hi[2], blo.z, bhi.z);
    const int ilo = amrex::Clamp(lo[0], blo.x, bhi.x);
    const int ihi = amrex::Clamp(hi[0], blo.x, bhi.x);
    const int jlo = amrex::Clamp(lo[1], blo.y, bhi.y);
    const int jhi = amrex::Clamp(hi[1], blo.y, bhi.y);
    for (int kk = klo; kk <= khi; ++kk) {
        for (int jj = jlo; jj <= jhi; ++jj) {
            for (int ii = ilo; ii <= ihi; ++ii) {
                const amrex::Real f = F(ii, jj, kk);
                all_gas = all_gas && (f <= pure_tol);
                all_liq = all_liq && (f >= 1.0_rt - pure_tol);
            }
            if (!all_gas && !all_liq) {
                return false;
            }
        }
    }
    a_liquid = all_liq;
    return true;
}

/**
 * \brief Whether a flux volume lies in pure cells only, from its bounding box:
 *        it then needs no cutting.
 *
 * \param liquid its signed liquid volume, if so
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
pure_flux_volume(const FluxVolume& fv,
                 const amrex::Real courant,
                 amrex::Array4<amrex::Real const> const& F,
                 const amrex::Box& bx,
                 amrex::Real& liquid) noexcept {
    amrex::Real pmin[3] = {fv.p[0].x[0], fv.p[0].x[1], fv.p[0].x[2]};
    amrex::Real pmax[3] = {pmin[0], pmin[1], pmin[2]};
    for (int n = 1; n < fv.npoint(); ++n) {
        for (int d = 0; d < 3; ++d) {
            pmin[d] = amrex::min(pmin[d], fv.p[n].x[d]);
            pmax[d] = amrex::max(pmax[d], fv.p[n].x[d]);
        }
    }
    // The cells the bounding box overlaps: a vertex on a mesh plane does not bring in the cell
    // beyond it, as the downwind cell of a flux volume that ends on its face
    int lo[3];
    int hi[3];
    for (int d = 0; d < 3; ++d) {
        lo[d] = static_cast<int>(std::floor(pmin[d]));
        hi[d] = amrex::max(lo[d], static_cast<int>(std::ceil(pmax[d])) - 1);
    }
    bool is_liquid = false;
    if (!pure_cells(lo, hi, F, bx, is_liquid)) {
        return false;
    }
    liquid = is_liquid ? courant : 0.0_rt;
    return true;
}

/**
 * \brief F of the upwind cell of the low face of cell iv in direction dir, for a flow u along dir:
 *        the cell below the face for u >= 0.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
upwind_fraction(amrex::Array4<amrex::Real const> const& F,
                amrex::IntVect iv,
                const int dir,
                const amrex::Real u) noexcept {
    if (u >= 0.0_rt) {
        iv[dir] -= 1;
    }
    return F(iv);
}

/**
 * \brief Liquid volume in the flux volume of the low face of cell (i,j,k) in direction dir, per
 *        cell volume.
 *
 * Where the flux volume falls short of the Courant number (see flux_volume), the liquid is its
 * liquid fraction times the Courant number, so that a face in liquid only carries the Courant number,
 * as on the faces that need no cuts: a cell full of liquid stays full. Without a volume to take a
 * fraction of, the fraction is F of the upwind cell.
 *
 * \param bx   cells where F and plane are valid
 * \return signed liquid volume crossing the face
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
flux_volume_liquid(const int i,
                    const int j,
                    const int k,
                    const int dir,
                    const amrex::Real courant,
                    amrex::Array4<amrex::Real const> const& F,
                    amrex::Array4<amrex::Real const> const& plane,
                    amrex::Array4<amrex::Real const> const& disp,
                    const amrex::Box& bx) noexcept {
    FluxVolume fv;
    flux_volume(i, j, k, dir, courant, disp, fv);
    amrex::Real liquid = 0.0_rt;
    // All in one phase: no cuts needed
    if (pure_flux_volume(fv, courant, F, bx, liquid)) {
        return liquid;
    }
    amrex::Real volume = 0.0_rt;
    for (int n = 0; n < fv.ntet; ++n) {
        const amrex::Real v = fv.vol[n];
        volume += v;
        if (v == 0.0_rt) {
            continue;
        }
        const amrex::Real s = (v > 0.0_rt) ? 1.0_rt : -1.0_rt;
        liquid += s * cut_tet_liquid(fv.tet(n), bx, F, plane);
    }
    if (!fv.matches_courant) {
        amrex::Real fraction = 0.0_rt;
        if (std::abs(volume) > apex_tol) {
            fraction = liquid / volume;
        } else {
            fraction = upwind_fraction(F, amrex::IntVect(AMREX_D_DECL(i, j, k)), dir, courant);
        }
        liquid = courant * amrex::Clamp(fraction, 0.0_rt, 1.0_rt);
    }
    return liquid;
}

/**
 * \brief The accessor a with its origin moved to cell (di,dj,dk): the result at (0,0,0) is a at
 *        (di,dj,dk).
 */
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Array4<T>
shift_origin(amrex::Array4<T> const& a, const int di, const int dj, const int dk) noexcept {
    const auto lo = amrex::lbound(a);
    const auto hi = amrex::ubound(a);
    return amrex::Array4<T>(a.p, amrex::Dim3{lo.x - di, lo.y - dj, lo.z - dk},
                            amrex::Dim3{hi.x + 1 - di, hi.y + 1 - dj, hi.z + 1 - dk}, a.nComp());
}

/**
 * \brief Signed liquid volume crossing the low face of cell (i,j,k) in direction dir, per cell
 *        volume, in coordinates relative to the face.
 *
 * \param disp displacement of each mesh node over the step (see flux_volume)
 * \param bx   cells where F and plane are valid
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
face_liquid_flux(const int i,
                 const int j,
                 const int k,
                 const int dir,
                 const amrex::Real courant,
                 amrex::Array4<amrex::Real const> const& F,
                 amrex::Array4<amrex::Real const> const& plane,
                 amrex::Array4<amrex::Real const> const& disp,
                 const amrex::Box& bx) noexcept {
    const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
    return flux_volume_liquid(0, 0, 0, dir, courant, shift_origin(F, i, j, k),
                              shift_origin(plane, i, j, k), shift_origin(disp, i, j, k),
                              amrex::shift(bx, -iv));
}

/**
 * \brief Displacement of mesh node (i,j,k) over the step, in cell units, by the midpoint rule: the
 *        velocity at the node, then the velocity half a step back along it. Velocities between
 *        faces are multilinear interpolations of the face velocities themselves, so that a node of
 *        a box reads no faces beyond the first ghost layer, and the arithmetic is relative to the
 *        node, so that every box holding it, and its periodic images, get the same displacement.
 *
 * \param umac     face velocities in each direction
 * \param dtdx     dt / dx in each direction, to Courant numbers
 * \param faces    faces normal to each direction that may be read: a face beyond them, across a
 *                 wall, is replaced by the nearest face on the domain boundary. The node reads the
 *                 faces up to one beyond it in every direction.
 * \param disp     output, 3 components (the third is 0 in 2D)
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
node_displacement(const int i,
                  const int j,
                  const int k,
                  amrex::GpuArray<amrex::Array4<amrex::Real const>, AMREX_SPACEDIM> const& umac,
                  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dtdx,
                  amrex::GpuArray<amrex::Box, AMREX_SPACEDIM> const& faces,
                  amrex::Real* disp) noexcept {
    const int node[3] = {i, j, k};
    // The Courant number of the face normal to d at node + off (normal) and cell node + off
    // (tangential)
    auto face = [&](const int d, const int* off) {
        int f[3] = {node[0] + off[0], node[1] + off[1], node[2] + off[2]};
        for (int t = 0; t < AMREX_SPACEDIM; ++t) {
            f[t] = amrex::Clamp(f[t], faces[d].smallEnd(t), faces[d].bigEnd(t));
        }
        // In 2D, k and off[2] are 0
        return umac[d](f[0], f[1], f[2]) * dtdx[d];
    };
    // Velocity at x, relative to the node, |x| < 1/2 in each direction. Face centers normal to d lie
    // at integer positions along d and half-integer positions along the others.
    auto velocity = [&](const amrex::Real* x, amrex::Real* u) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            const int b = (x[d] < 0.0_rt) ? -1 : 0;
            const amrex::Real w[3] = {
                (d == 0) ? x[0] - static_cast<amrex::Real>(b) : x[0] + 0.5_rt,
                (d == 1) ? x[1] - static_cast<amrex::Real>(b) : x[1] + 0.5_rt,
                (d == 2) ? x[2] - static_cast<amrex::Real>(b) : x[2] + 0.5_rt};
            const int base[3] = {(d == 0) ? b : -1, (d == 1) ? b : -1, (d == 2) ? b : -1};
            amrex::Real sum = 0.0_rt;
#if AMREX_SPACEDIM == 2
            for (int s1 = 0; s1 <= 1; ++s1) {
                for (int s0 = 0; s0 <= 1; ++s0) {
                    const int off[3] = {base[0] + s0, base[1] + s1, 0};
                    sum += (s0 != 0 ? w[0] : 1.0_rt - w[0]) * (s1 != 0 ? w[1] : 1.0_rt - w[1]) *
                           face(d, off);
                }
            }
#else
            for (int s2 = 0; s2 <= 1; ++s2) {
                for (int s1 = 0; s1 <= 1; ++s1) {
                    for (int s0 = 0; s0 <= 1; ++s0) {
                        const int off[3] = {base[0] + s0, base[1] + s1, base[2] + s2};
                        sum += (s0 != 0 ? w[0] : 1.0_rt - w[0]) * (s1 != 0 ? w[1] : 1.0_rt - w[1]) *
                               (s2 != 0 ? w[2] : 1.0_rt - w[2]) * face(d, off);
                    }
                }
            }
#endif
            u[d] = sum;
        }
    };
    const amrex::Real x0[3] = {0.0_rt, 0.0_rt, 0.0_rt};
    amrex::Real u0[3] = {0.0_rt, 0.0_rt, 0.0_rt};
    velocity(x0, u0);
    const amrex::Real xm[3] = {-0.5_rt * u0[0], -0.5_rt * u0[1], -0.5_rt * u0[2]};
    amrex::Real um[3] = {0.0_rt, 0.0_rt, 0.0_rt};
    velocity(xm, um);
    for (int d = 0; d < 3; ++d) {
        disp[d] = -um[d];
    }
}

} // namespace PLIC

#endif
/** @}*/

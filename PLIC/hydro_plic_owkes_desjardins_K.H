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
 * Pieces with d <= 0 go to neg and those with d > 0 to pos.
 *
 * \param n_neg,n_pos number of pieces written, at most 3 each
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
split_tet(const Tet& t, const amrex::Real* d, Tet* neg, int& n_neg, Tet* pos, int& n_pos) noexcept {
    int in[4];
    int out[4];
    int ni = 0;
    int no = 0;
    for (int v = 0; v < 4; ++v) {
        if (d[v] <= 0.0_rt) {
            in[ni++] = v;
        } else {
            out[no++] = v;
        }
    }
    n_neg = 0;
    n_pos = 0;
    if (no == 0) {
        neg[0] = t;
        n_neg = 1;
        return;
    }
    if (ni == 0) {
        pos[0] = t;
        n_pos = 1;
        return;
    }

    // Point where the edge joining an inside and an outside vertex crosses
    auto cross = [&](const int a, const int b) {
        return lerp(t.v[a], t.v[b], d[a] / (d[a] - d[b]));
    };
    if (ni == 1) {
        const int a = in[0];
        const Vec3 p0 = cross(a, out[0]);
        const Vec3 p1 = cross(a, out[1]);
        const Vec3 p2 = cross(a, out[2]);
        neg[0] = Tet{{t.v[a], p0, p1, p2}};
        n_neg = 1;
        prism_to_tets(t.v[out[0]], t.v[out[1]], t.v[out[2]], p0, p1, p2, pos);
        n_pos = 3;
    } else if (no == 1) {
        const int a = out[0];
        const Vec3 p0 = cross(in[0], a);
        const Vec3 p1 = cross(in[1], a);
        const Vec3 p2 = cross(in[2], a);
        pos[0] = Tet{{t.v[a], p0, p1, p2}};
        n_pos = 1;
        prism_to_tets(t.v[in[0]], t.v[in[1]], t.v[in[2]], p0, p1, p2, neg);
        n_neg = 3;
    } else {
        const int a = in[0];
        const int b = in[1];
        const int c = out[0];
        const int e = out[1];
        const Vec3 ac = cross(a, c);
        const Vec3 ae = cross(a, e);
        const Vec3 bc = cross(b, c);
        const Vec3 be = cross(b, e);
        prism_to_tets(t.v[a], ac, ae, t.v[b], bc, be, neg);
        n_neg = 3;
        prism_to_tets(t.v[c], ac, bc, t.v[e], ae, be, pos);
        n_pos = 3;
    }
}

/**
 * \brief Volume of the part of a tetrahedron on the liquid side of the plane
 *        m.(x - x_lo) <= alpha of the cell with low corner x_lo.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real tet_liquid_volume(
    const Tet& t, const int i, const int j, const int k, const amrex::Real* plane) noexcept {
    amrex::Real d[4];
    for (int v = 0; v < 4; ++v) {
        d[v] = plane[0] * (t.v[v].x[0] - static_cast<amrex::Real>(i)) +
               plane[1] * (t.v[v].x[1] - static_cast<amrex::Real>(j)) +
               plane[2] * (t.v[v].x[2] - static_cast<amrex::Real>(k)) - plane[3];
    }
    Tet neg[3];
    Tet pos[3];
    int n_neg = 0;
    int n_pos = 0;
    split_tet(t, d, neg, n_neg, pos, n_pos);
    amrex::Real vol = 0.0_rt;
    for (int n = 0; n < n_neg; ++n) {
        vol += abs_volume(neg[n]);
    }
    return vol;
}

// Tetrahedra held while cutting by the mesh planes.  Each cut replaces one
// tetrahedron by at most 6 and a depth-first path crosses few planes for a
// Courant number below 1: 64 is well above what can be reached.
static constexpr int cut_stack_size = 64;

/**
 * \brief Cut a tetrahedron along the mesh planes into one piece per cell and
 *        call visit(piece, volume, i, j, k) on each.
 *
 * \param bx cells where the data are valid; pieces beyond it are charged to
 *           the nearest cell inside it
 */
template <typename Visit>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
for_each_cell_piece(const Tet& tet, const amrex::Box& bx, Visit const& visit) noexcept {
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

        if (cut_dir >= 0 && top + 6 <= cut_stack_size) {
            const amrex::Real d[4] = {t.v[0].x[cut_dir] - cut_at, t.v[1].x[cut_dir] - cut_at,
                                      t.v[2].x[cut_dir] - cut_at, t.v[3].x[cut_dir] - cut_at};
            Tet neg[3];
            Tet pos[3];
            int n_neg = 0;
            int n_pos = 0;
            split_tet(t, d, neg, n_neg, pos, n_pos);
            for (int n = 0; n < n_neg; ++n) {
                stack[top++] = neg[n];
            }
            for (int n = 0; n < n_pos; ++n) {
                stack[top++] = pos[n];
            }
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
        visit(t, vol, ci, cj, ck);
    }
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
    return tet_liquid_volume(t, i, j, k, pl);
}

/**
 * \brief Apex o of the two tetrahedra (o,a,c,e) and (o,a,b,c) added on a
 *        traced face (a,b,c,e) to give them a total signed volume vcor.  o lies
 *        over the barycenter of the face; the volume is linear in o along dir.
 *
 * \return false if the face has no extent normal to dir (no correction)
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool correction_apex(const Vec3& a,
                                                              const Vec3& b,
                                                              const Vec3& c,
                                                              const Vec3& e,
                                                              const int dir,
                                                              const amrex::Real vcor,
                                                              Vec3& o) noexcept {
    const Vec3* tri[2][3] = {{&a, &c, &e}, {&a, &b, &c}}; // +dir normals
    for (int d = 0; d < 3; ++d) {
        o.x[d] = 0.25_rt * (a.x[d] + b.x[d] + c.x[d] + e.x[d]);
    }
    amrex::Real s0 = 0.0_rt; // 6 V at o_dir = 0
    amrex::Real sd = 0.0_rt; // d(6 V)/d(o_dir)
    for (auto& n : tri) {
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
 * \brief Signed tetrahedra making up the flux volume of the low face of cell
 *        (i,j,k) in direction dir.
 *
 * \param disp    displacement of each mesh node over the step, from its traced
 *                position at the start of the step, in cell units (3 components)
 * \param courant the face velocity times dt/dx, which the flux volume is
 *                corrected to match
 * \return the number of tetrahedra in tets
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int
flux_volume_tets(const int i,
                 const int j,
                 const int k,
                 const int dir,
                 const amrex::Real courant,
                 amrex::Array4<amrex::Real const> const& disp,
                 Tet* tets) noexcept {
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
#if AMREX_SPACEDIM == 2
        // Nodes of a 2D mesh: z stays as it was, and only the in-plane offset
        // selects the node
        const int ni = i + s1 * off1[0] + s2 * off2[0];
        const int nj = j + s1 * off1[1] + s2 * off2[1];
        const auto z = static_cast<amrex::Real>(k + s1 * off1[2] + s2 * off2[2]);
        return Vec3{{static_cast<amrex::Real>(ni) + disp(ni, nj, 0, 0),
                     static_cast<amrex::Real>(nj) + disp(ni, nj, 0, 1), z}};
#else
        const int ni = i + s1 * off1[0] + s2 * off2[0];
        const int nj = j + s1 * off1[1] + s2 * off2[1];
        const int nk = k + s1 * off1[2] + s2 * off2[2];
        return Vec3{{static_cast<amrex::Real>(ni) + disp(ni, nj, nk, 0),
                     static_cast<amrex::Real>(nj) + disp(ni, nj, nk, 1),
                     static_cast<amrex::Real>(nk) + disp(ni, nj, nk, 2)}};
#endif
    };

    // Face corners counterclockwise about +dir, and their traced positions
    const Vec3 a = node(0, 0);
    const Vec3 b = node(1, 0);
    const Vec3 c = node(1, 1);
    const Vec3 e = node(0, 1);
    const Vec3 ap = traced(0, 0);
    const Vec3 bp = traced(1, 0);
    const Vec3 cp = traced(1, 1);
    const Vec3 ep = traced(0, 1);

    // Fan of the boundary from a.  With outward triangles the signed volumes
    // add up to the flux volume, positive for a flow along +dir.  Triangles
    // touching a vanish: the face itself and the side faces through a.  The
    // diagonals of the side faces depend only on the mesh edge (see above).
    constexpr int ntri = 8;
    const Vec3* tri[ntri][3] = {{&ap, &bp, &b},  // side a-b, diagonal a' - b
                                {&e, &ep, &ap},  // side e-a, diagonal e - a'
                                {&b, &bp, &c},   // side b-c, diagonal b' - c
                                {&bp, &cp, &c},  //
                                {&c, &cp, &ep},  // side c-e, diagonal e' - c
                                {&c, &ep, &e},   //
                                {&ap, &ep, &cp}, // traced face, diagonal a' - c'
                                {&ap, &cp, &bp}};

    amrex::Real vol = 0.0_rt;
    for (auto const& n : tri) {
        vol += signed_volume(a, *n[0], *n[1], *n[2]);
    }

    // Extra volume so that the flux volume matches the face velocity exactly:
    // two tetrahedra from the traced face to a point o above its barycenter, in
    // the direction of dir (Eq. 27-28).  Their volume is linear in o along dir.
    const Vec3* ttri[2][3] = {{&ap, &cp, &ep}, {&ap, &bp, &cp}}; // +dir normals
    Vec3 o;
    const bool has_cor = correction_apex(ap, bp, cp, ep, dir, courant - vol, o);

    int ntet = 0;
    for (auto const& n : tri) {
        tets[ntet++] = Tet{{a, *n[0], *n[1], *n[2]}};
    }
    if (has_cor) {
        for (auto& n : ttri) {
            tets[ntet++] = Tet{{o, *n[0], *n[1], *n[2]}};
        }
    }
    return ntet;
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
#if AMREX_SPACEDIM == 2
    const int klo = 0;
    const int khi = 0;
#else
    const int klo = amrex::Clamp(lo[2], blo.z, bhi.z);
    const int khi = amrex::Clamp(hi[2], blo.z, bhi.z);
#endif
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
pure_flux_volume(const Tet* tets,
                 const int ntet,
                 const amrex::Real courant,
                 amrex::Array4<amrex::Real const> const& F,
                 const amrex::Box& bx,
                 amrex::Real& liquid) noexcept {
    amrex::Real pmin[3] = {tets[0].v[0].x[0], tets[0].v[0].x[1], tets[0].v[0].x[2]};
    amrex::Real pmax[3] = {pmin[0], pmin[1], pmin[2]};
    for (int n = 0; n < ntet; ++n) {
        for (const auto& vert : tets[n].v) {
            for (int d = 0; d < 3; ++d) {
                pmin[d] = amrex::min(pmin[d], vert.x[d]);
                pmax[d] = amrex::max(pmax[d], vert.x[d]);
            }
        }
    }
    int lo[3];
    int hi[3];
    for (int d = 0; d < 3; ++d) {
        lo[d] = static_cast<int>(std::floor(pmin[d]));
        hi[d] = static_cast<int>(std::floor(pmax[d]));
    }
    bool is_liquid = false;
    if (!pure_cells(lo, hi, F, bx, is_liquid)) {
        return false;
    }
    liquid = is_liquid ? courant : 0.0_rt;
    return true;
}

/**
 * \brief Liquid volume in the flux volume of the low face of cell (i,j,k) in direction dir, per
 *        cell volume, calling visit(liquid, gas, ci, cj, ck) with the signed liquid and gas volumes
 *        of each piece of it lying in cell (ci,cj,ck).
 *
 * \param bx   cells where F and plane are valid
 * \return signed liquid volume crossing the face
 */
template <typename Visit>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
for_each_flux_piece(const int i,
                    const int j,
                    const int k,
                    const int dir,
                    const amrex::Real courant,
                    amrex::Array4<amrex::Real const> const& F,
                    amrex::Array4<amrex::Real const> const& plane,
                    amrex::Array4<amrex::Real const> const& disp,
                    const amrex::Box& bx,
                    Visit const& visit) noexcept {
    Tet tets[max_flux_tets];
    const int ntet = flux_volume_tets(i, j, k, dir, courant, disp, tets);
    amrex::Real pure_liquid = 0.0_rt;
    const bool pure = pure_flux_volume(tets, ntet, courant, F, bx, pure_liquid);
    amrex::Real liquid = 0.0_rt;
    for (int n = 0; n < ntet; ++n) {
        const Tet& tet = tets[n];
        const amrex::Real v = signed_volume(tet.v[0], tet.v[1], tet.v[2], tet.v[3]);
        if (v == 0.0_rt) {
            continue;
        }
        const amrex::Real s = (v > 0.0_rt) ? 1.0_rt : -1.0_rt;
        amrex::Real tet_liquid = 0.0_rt;
        for_each_cell_piece(
            tet, bx, [&](const Tet& t, const amrex::Real vol, int ci, int cj, int ck) noexcept {
                const amrex::Real l = piece_liquid(t, vol, ci, cj, ck, F, plane);
                tet_liquid += l;
                visit(s * l, s * (vol - l), ci, cj, ck);
            });
        liquid += s * tet_liquid;
    }
    return pure ? pure_liquid : liquid;
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
 * \param disp displacement of each mesh node over the step (see flux_volume_tets)
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
    return for_each_flux_piece(0, 0, 0, dir, courant, shift_origin(F, i, j, k),
                               shift_origin(plane, i, j, k), shift_origin(disp, i, j, k),
                               amrex::shift(bx, -iv),
                               [](const amrex::Real, const amrex::Real, int, int, int) noexcept {});
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
 * \param dom_lo   low cell of the domain
 * \param dom_hi   high cell of the domain. Faces beyond a non-periodic domain boundary are
 *                 replaced by the nearest face inside it.
 * \param periodic whether each direction is periodic
 * \param disp     output, 3 components (the third is 0 in 2D)
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
node_displacement(const int i,
                  const int j,
                  const int k,
                  amrex::GpuArray<amrex::Array4<amrex::Real const>, AMREX_SPACEDIM> const& umac,
                  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dtdx,
                  amrex::GpuArray<int, AMREX_SPACEDIM> const& dom_lo,
                  amrex::GpuArray<int, AMREX_SPACEDIM> const& dom_hi,
                  amrex::GpuArray<int, AMREX_SPACEDIM> const& periodic,
                  amrex::Real* disp) noexcept {
    const int node[3] = {i, j, k};
    // The Courant number of the face normal to d at node + off (normal) and cell node + off
    // (tangential)
    auto face = [&](const int d, const int* off) {
        int f[3] = {node[0] + off[0], node[1] + off[1], node[2] + off[2]};
        for (int t = 0; t < AMREX_SPACEDIM; ++t) {
            if (periodic[t] == 0) {
                f[t] = amrex::Clamp(f[t], dom_lo[t], (t == d) ? dom_hi[t] + 1 : dom_hi[t]);
            }
        }
#if AMREX_SPACEDIM == 2
        return umac[d](f[0], f[1], 0) * dtdx[d];
#else
        return umac[d](f[0], f[1], f[2]) * dtdx[d];
#endif
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

/**
 * \file hydro_plic_reconstruction.cpp
 * \addtogroup PLIC
 * @{
 */

#include <hydro_plic.H>
#include <hydro_plic_normals_K.H>

#include <AMReX_BC_TYPES.H>
#include <AMReX_IArrayBox.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_GpuBuffer.H>
#include <AMReX_TagParallelFor.H>

#include <limits>
#include <string>
#include <utility>
#include <vector>

using namespace amrex::literals;

namespace PLIC {

namespace {

//! The domain faces whose ghost cells FillWallGhosts overwrites, nonzero on the low (lo) and high
//! (hi) side of each direction. Never those of a periodic direction.
struct Walls {
    amrex::IntVect lo{0};
    amrex::IntVect hi{0};
};

//! The part of box b beyond the wall on side sgn of direction dir
amrex::Box beyond_wall (amrex::Box b, amrex::Box const& domain, int dir, int sgn)
{
    if (sgn < 0) {
        b.setBig(dir, domain.smallEnd(dir) - 1);
    } else {
        b.setSmall(dir, domain.bigEnd(dir) + 1);
    }
    return b;
}

} // namespace

namespace detail {

// The boxes and arrays of the fused kernels of FillWallGhosts: each kernel runs over the cells of the
// box of every tag at once (amrex::ParallelFor over tags), rather than one launch per box

//! a mask a filtered over its box bx, with p as scratch on the same box
struct FilterTag {
    amrex::Array4<int> a;
    amrex::Array4<int> p;
    amrex::Box bx;
    [[nodiscard]] AMREX_GPU_HOST_DEVICE amrex::Box const& box () const noexcept { return bx; }
};

// The steps of the fill of one wall, each over bx in one box of the copy w

//! The wall cells with a mixed cell near: near
struct NearTag {
    amrex::Array4<amrex::Real const> w;
    amrex::Array4<int> near;
    amrex::Box bx;
    [[nodiscard]] AMREX_GPU_HOST_DEVICE amrex::Box const& box () const noexcept { return bx; }
};

//! The cells that hold data (covered by a box, or in the box of F with its ghost cells own) and lie
//! in the box inside: fit_ok
struct FitTag {
    amrex::Array4<int> fit_ok;
    amrex::Array4<int const> covered;
    amrex::Box bx;
    amrex::Box own;
    amrex::Box inside;
    [[nodiscard]] AMREX_GPU_HOST_DEVICE amrex::Box const& box () const noexcept { return bx; }
};

//! What each cell offers the search, cand, with its plane pl
struct CandTag {
    amrex::Array4<amrex::Real const> w;
    amrex::Array4<int const> fit_ok;
    amrex::Array4<int> cand;
    amrex::Array4<amrex::Real> pl;
    amrex::Box bx;
    [[nodiscard]] AMREX_GPU_HOST_DEVICE amrex::Box const& box () const noexcept { return bx; }
};

//! The ghost values gbx beyond the wall cells bx, added to the sums acc
struct PlaneTag {
    amrex::Array4<amrex::Real const> w;
    amrex::Array4<int const> near;
    amrex::Array4<int const> cand;
    amrex::Array4<amrex::Real const> pl;
    amrex::Array4<amrex::Real> acc;
    amrex::Box bx;
    amrex::Box gbx;
    [[nodiscard]] AMREX_GPU_HOST_DEVICE amrex::Box const& box () const noexcept { return bx; }
};

//! One box of the copy w and of the sums of the ghost values acc, over its valid cells bx
struct MeanTag {
    amrex::Array4<amrex::Real const> w;
    amrex::Array4<amrex::Real> acc;
    amrex::Box bx;
    [[nodiscard]] AMREX_GPU_HOST_DEVICE amrex::Box const& box () const noexcept { return bx; }
};

//! The walls, other than those of direction dir (none for dir = -1): whether a cell lies beyond one,
//! and the nearest cell that does not
struct OtherWalls {
    amrex::IntVect wall_lo, wall_hi, dlo, dhi;
    int dir;
    [[nodiscard]] AMREX_GPU_HOST_DEVICE bool beyond (amrex::IntVect const& c) const noexcept {
        return clamp(c) != c;
    }
    [[nodiscard]] AMREX_GPU_HOST_DEVICE amrex::IntVect clamp (amrex::IntVect c) const noexcept {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            if (d != dir) {
                if (wall_lo[d] != 0) { c[d] = amrex::max(c[d], dlo[d]); }
                if (wall_hi[d] != 0) { c[d] = amrex::min(c[d], dhi[d]); }
            }
        }
        return c;
    }
};

// Replace each cell of the mask of each tag, one direction t with dirs[t] != 0 at a time, by whether
// all (all = true) or any (all = false) of the cells within radius along t are nonzero. Cells outside
// its box are zero.
void SeparableFilter (amrex::Vector<FilterTag> const& tags, int radius, amrex::IntVect const& dirs,
                      bool all)
{
    if (tags.empty()) {
        return;
    }
    const amrex::TagVector<FilterTag> tv(tags);
    for (int t = 0; t < AMREX_SPACEDIM; ++t) {
        if (dirs[t] == 0) {
            continue;
        }
        const amrex::IntVect et = amrex::IntVect::TheDimensionVector(t);
        amrex::ParallelFor(tv, [=] AMREX_GPU_DEVICE (int i, int j, int k, FilterTag const& g) noexcept
        {
            g.p(i, j, k) = g.a(i, j, k);
        });
        amrex::ParallelFor(tv, [=] AMREX_GPU_DEVICE (int i, int j, int k, FilterTag const& g) noexcept
        {
            amrex::ignore_unused(k);
            const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
            // Stops at the first cell that decides: a zero for all, a nonzero for any
            int r = all ? 1 : 0;
            for (int a = -radius; a <= radius && (r != 0) == all; ++a) {
                const amrex::IntVect c = iv + a * et;
                r = (g.bx.contains(c) && g.p(c) != 0) ? 1 : 0;
            }
            g.a(i, j, k) = r;
        });
    }
}

} // namespace detail

template <typename Normal>
void FitPlanesImpl (amrex::MultiFab& plane, amrex::MultiFab const& F, int ngrow)
{
    BL_PROFILE("PLIC::FitPlanes()");
    AMREX_ALWAYS_ASSERT(plane.nComp() == 4 && plane.nGrow() >= ngrow);
    AMREX_ALWAYS_ASSERT(F.nGrow() >= ngrow + Normal::stencil);
    AMREX_ALWAYS_ASSERT(plane.boxArray() == F.boxArray() &&
                        plane.DistributionMap() == F.DistributionMap());

    // The ghost layers beyond ngrow; the kernel writes every cell within it
    if (plane.nGrow() > ngrow) {
        plane.setVal(0.0_rt);
    }
    auto const& Fa = F.const_arrays();
    auto const& pl = plane.arrays();
    amrex::ParallelFor(plane, amrex::IntVect(ngrow),
                       [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept
    {
        amrex::Real mx = 0.0_rt;
        amrex::Real my = 0.0_rt;
        amrex::Real mz = 0.0_rt;
        amrex::Real alpha = 0.0_rt;
        if (is_mixed(Fa[b](i, j, k))) {
            fit_plane<Normal>(i, j, k, Fa[b], mx, my, mz, alpha);
        }
        pl[b](i, j, k, 0) = mx;
        pl[b](i, j, k, 1) = my;
        pl[b](i, j, k, 2) = mz;
        pl[b](i, j, k, 3) = alpha;
    });
}

// A ghost value depends on cells further along the wall than F has ghost cells, so that F alone
// would give different values in neighbouring boxes. The fill therefore works on a copy of the boxes
// that touch a wall, extended through it: each wall ghost cell belongs to one box, and its ghost
// layers cover the whole search. Each wall reads only cells in the domain, beyond itself or beyond
// a periodic face, so the walls do not depend on one another or on their order, and a ghost cell beyond several walls
// (an edge or a corner) takes the mean over those that extend a plane to it. The result does not
// depend on the box layout, except near a coarse/fine boundary: past it the copy holds no more ghost layers of data
// than F, and the search only takes planes whose stencil lies on data.
template <typename Normal>
void FillWallGhostsImpl (amrex::MultiFab& a_F, amrex::Geometry const& geom, Walls const& walls)
{
    BL_PROFILE("PLIC::FillWallGhosts()");
    const int nGhost = a_F.nGrow();
    const amrex::IntVect wall_lo = walls.lo;
    const amrex::IntVect wall_hi = walls.hi;
    const auto is_wall = [&] (const int dir, const int sgn) {
        return ((sgn < 0) ? wall_lo : wall_hi)[dir] != 0;
    };
    if (wall_lo == amrex::IntVect(0) && wall_hi == amrex::IntVect(0)) {
        return;
    }
    AMREX_ALWAYS_ASSERT(nGhost >= Normal::stencil + 1);

    // The ghost layers read: the search along the wall and the reconstruction stencil of the plane
    // it finds, or the ghost layers of F themselves if there are more
    const int nread = amrex::max(wall_nsearch + Normal::stencil, nGhost);
    const amrex::Box& domain = geom.Domain();
    const auto periodicity = geom.periodicity();

    // The boxes that touch a wall, extended through it. They stay disjoint, and keep the owner of
    // the box of F they extend, f_index. With the other boxes of F, all_boxes.
    amrex::BoxList all_boxes;
    amrex::BoxList boxes;
    amrex::Vector<int> owners;
    amrex::Vector<int> f_index;
    for (int ib = 0; ib < a_F.boxArray().size(); ++ib) {
        const amrex::Box b = a_F.boxArray()[ib];
        amrex::Box ext = b;
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            if (is_wall(dir, -1) && b.smallEnd(dir) == domain.smallEnd(dir)) {
                ext.growLo(dir, nGhost);
            }
            if (is_wall(dir, 1) && b.bigEnd(dir) == domain.bigEnd(dir)) {
                ext.growHi(dir, nGhost);
            }
        }
        all_boxes.push_back(ext);
        if (ext != b) {
            boxes.push_back(ext);
            owners.push_back(a_F.DistributionMap()[ib]);
            f_index.push_back(ib);
        }
    }
    if (boxes.isEmpty()) {
        return;
    }
    const amrex::BoxArray ext_ba(std::move(boxes));
    const amrex::BoxArray all_ba(std::move(all_boxes));
    const amrex::DistributionMapping ext_dm(std::move(owners));

    // The copy holds, in each cell, the value of the box whose valid region covers it (marked in
    // covered), or else the ghost cell of F of its own box (beyond the wall or over a coarse/fine
    // boundary). Other cells hold no data.
    amrex::MultiFab work(ext_ba, ext_dm, 1, nread);
    amrex::iMultiFab covered(ext_ba, ext_dm, 1, nread);
    work.setVal(0.0_rt);
    covered.setVal(0);
    {
        // The cells of the boxes of F and of their extensions, with their periodic images, are known
        // from the box lists, so no exchange is needed
        const std::vector<amrex::IntVect> shifts = periodicity.shiftIntVect(amrex::IntVect(nread));
        amrex::Vector<amrex::Array4PairTag<amrex::Real>> own_tags;
        amrex::Vector<amrex::Array4BoxTag<int>> cover_tags;
        std::vector<std::pair<int, amrex::Box>> isects;
        for (amrex::MFIter mfi(work); mfi.isValid(); ++mfi) {
            const int fb = f_index[mfi.index()];
            own_tags.push_back({work.array(mfi), a_F.const_array(fb), a_F[fb].box()});
            for (amrex::IntVect const& shift : shifts) {
                all_ba.intersections(work[mfi].box() + shift, isects);
                for (auto const& is : isects) {
                    cover_tags.push_back({covered.array(mfi), is.second - shift});
                }
            }
        }
        amrex::ParallelFor(own_tags, [=] AMREX_GPU_DEVICE (int i, int j, int k,
                                                           amrex::Array4PairTag<amrex::Real> const& g) noexcept
        {
            g.dfab(i, j, k) = g.sfab(i, j, k);
        });
        amrex::ParallelFor(cover_tags, [=] AMREX_GPU_DEVICE (int i, int j, int k,
                                                             amrex::Array4BoxTag<int> const& g) noexcept
        {
            g.dfab(i, j, k) = 1;
        });
    }
    work.ParallelCopy(a_F, 0, 0, 1, 0, nread, periodicity);
    work.FillBoundary(periodicity);

    // The domain and its ghost layers, without bounds in periodic directions
    amrex::Box extent = domain;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        extent.grow(dir, geom.isPeriodic(dir) ? 1000000 : nGhost);
    }

    // On the valid cells of the copy, which hold every wall ghost cell: the sum of the ghost values
    // of the walls that extend a plane to the cell, and their number
    const amrex::IntVect dlo = domain.smallEnd();
    const amrex::IntVect dhi = domain.bigEnd();
    amrex::MultiFab acc(ext_ba, ext_dm, 2, 0);
    acc.setVal(0.0_rt);
    auto const& acc_ma = acc.arrays();

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        for (int sgn = -1; sgn <= 1; sgn += 2) {
            if (!is_wall(dir, sgn)) {
                continue;
            }
            const int iwall = (sgn < 0) ? domain.smallEnd(dir) : domain.bigEnd(dir);
            amrex::IntVect along(1);
            along[dir] = 0;
            const detail::OtherWalls ow{wall_lo, wall_hi, dlo, dhi, dir};
            // The cells the fill of this wall may read: in the domain, or beyond this wall or a
            // periodic face. Beyond another face F holds data only in the ghost cells of some
            // boxes, which would make the fill depend on the box layout.
            amrex::Box inside = domain;
            for (int t = 0; t < AMREX_SPACEDIM; ++t) {
                if (geom.isPeriodic(t)) {
                    inside.grow(t, 1000000);
                }
            }
            if (sgn < 0) {
                inside.growLo(dir, 1000000);
            } else {
                inside.growHi(dir, 1000000);
            }

            // The boxes of the copy that reach beyond this wall, with the regions of each step:
            //   scan  the wall cells within the search of the valid wall cells
            //   wbx   the valid wall cells, each with one plane shared by the ghost cells beyond it
            //   slab  the cells searched from them, up to Normal::stencil layers inside the wall
            //   fitbx the slab grown by the stencil, where the data mask is eroded
            //   gbx   the ghost cells beyond the valid wall cells
            struct Regions {
                int li;
                amrex::Box scan, wbx, slab, fitbx, gbx, own;
            };
            amrex::Vector<Regions> regions;
            for (amrex::MFIter mfi(work); mfi.isValid(); ++mfi) {
                const amrex::Box vb = mfi.validbox();
                Regions r;
                r.li = mfi.LocalIndex();
                r.gbx = beyond_wall(vb, domain, dir, sgn);
                if (!r.gbx.ok()) {
                    continue;
                }
                r.scan = amrex::grow(vb, wall_nsearch);
                r.scan.setRange(dir, iwall);
                r.scan &= (work[mfi].box() & extent);
                r.wbx = vb;
                r.wbx.setRange(dir, iwall);
                r.slab = r.wbx;
                for (int t = 0; t < AMREX_SPACEDIM; ++t) {
                    if (t != dir) {
                        r.slab.grow(t, wall_nsearch);
                    }
                }
                r.slab.setRange(dir, (sgn < 0) ? iwall : iwall - Normal::stencil, Normal::stencil + 1);
                r.fitbx = amrex::grow(r.slab, Normal::stencil) & work[mfi].box();
                r.own = a_F[f_index[mfi.index()]].box();
                regions.push_back(r);
            }
            const auto nr = static_cast<int>(regions.size());
            auto const& w_ma = work.const_arrays();
            auto const& c_ma = covered.const_arrays();
            amrex::Vector<amrex::IArrayBox> near(nr);
            amrex::Vector<amrex::IArrayBox> near_scratch(nr);
            for (int n = 0; n < nr; ++n) {
                near[n].resize(regions[n].scan, 1, amrex::The_Async_Arena());
                near_scratch[n].resize(regions[n].scan, 1, amrex::The_Async_Arena());
            }

            // The wall cells near the interface: those with a mixed wall cell inside, within the
            // search. Elsewhere the ghost cells take the wall cell value. Without
            // a mixed wall cell in its boxes, which the filter only spreads, a rank extends no plane
            // on this wall.
            {
                amrex::Gpu::Buffer<int> flag({0});
                int* p_flag = flag.data();
                amrex::Vector<detail::NearTag> tags;
                amrex::Vector<detail::FilterTag> ftags;
                for (int n = 0; n < nr; ++n) {
                    tags.push_back({w_ma[regions[n].li], near[n].array(), regions[n].scan});
                    ftags.push_back({near[n].array(), near_scratch[n].array(), regions[n].scan});
                }
                amrex::ParallelFor(tags, [=] AMREX_GPU_DEVICE (int i, int j, int k, detail::NearTag const& g) noexcept
                {
                    amrex::ignore_unused(k);
                    const amrex::IntVect c(AMREX_D_DECL(i, j, k));
                    const bool mixed = is_mixed(g.w(c)) && inside.contains(c);
                    g.near(i, j, k) = static_cast<int>(mixed);
                    if (mixed) {
                        amrex::Gpu::Atomic::Max(p_flag, 1);
                    }
                });
                if (*flag.copyToHost() == 0) {
                    continue;
                }
                detail::SeparableFilter(ftags, wall_nsearch, along, false);
            }

            // The cells whose reconstruction stencil reads only data inside
            // (on a fine level, the layers past the ghost cells of F over a coarse/fine boundary hold
            // no data), and what each cell of the slab offers the search with its plane, computed once
            // for the overlapping searches
            amrex::Vector<amrex::IArrayBox> fit_ok(nr);
            amrex::Vector<amrex::IArrayBox> fit_scratch(nr);
            amrex::Vector<amrex::IArrayBox> cand(nr);
            amrex::Vector<amrex::FArrayBox> cand_plane(nr);
            for (int n = 0; n < nr; ++n) {
                fit_ok[n].resize(regions[n].fitbx, 1, amrex::The_Async_Arena());
                fit_scratch[n].resize(regions[n].fitbx, 1, amrex::The_Async_Arena());
                cand[n].resize(regions[n].slab, 1, amrex::The_Async_Arena());
                cand_plane[n].resize(regions[n].slab, 4, amrex::The_Async_Arena());
            }
            amrex::Vector<detail::FitTag> fit_tags;
            amrex::Vector<detail::FilterTag> fit_ftags;
            amrex::Vector<detail::CandTag> cand_tags;
            amrex::Vector<detail::PlaneTag> plane_tags;
            for (int n = 0; n < nr; ++n) {
                Regions const& r = regions[n];
                auto const& Fw = w_ma[r.li];
                fit_tags.push_back({fit_ok[n].array(), c_ma[r.li], r.fitbx, r.own, inside});
                fit_ftags.push_back({fit_ok[n].array(), fit_scratch[n].array(), r.fitbx});
                cand_tags.push_back({Fw, fit_ok[n].const_array(), cand[n].array(), cand_plane[n].array(), r.slab});
                plane_tags.push_back({Fw, near[n].const_array(), cand[n].const_array(), cand_plane[n].const_array(),
                                      acc_ma[r.li], r.wbx, r.gbx});
            }
            amrex::ParallelFor(fit_tags, [=] AMREX_GPU_DEVICE (int i, int j, int k, detail::FitTag const& g) noexcept
            {
                amrex::ignore_unused(k);
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                g.fit_ok(i, j, k) = ((g.covered(i, j, k) != 0 || g.own.contains(iv)) && g.inside.contains(iv)) ? 1 : 0;
            });
            detail::SeparableFilter(fit_ftags, Normal::stencil, amrex::IntVect(1), true);
            amrex::ParallelFor(cand_tags, [=] AMREX_GPU_DEVICE (int i, int j, int k, detail::CandTag const& g) noexcept
            {
                const amrex::IntVect c(AMREX_D_DECL(i, j, k));
                amrex::Real m[3] = {0.0_rt, 0.0_rt, 0.0_rt};
                amrex::Real alpha = 0.0_rt;
                g.cand(i, j, k) = wall_candidate<Normal>(c, dir, iwall, sgn, g.fit_ok, g.w, m, alpha);
                for (int d = 0; d < 3; ++d) {
                    g.pl(i, j, k, d) = m[d];
                }
                g.pl(i, j, k, 3) = alpha;
            });

            // The plane of each wall cell, the mean of the equally near planes found, extended to the
            // ghost cells beyond it. Each ghost value keeps the offset of the wall cell from the plane
            // (exact for a uniform normal to the wall), from the nearest wall cell not beyond another
            // wall.
            amrex::ParallelFor(plane_tags, [=] AMREX_GPU_DEVICE (int i, int j, int k, detail::PlaneTag const& g) noexcept
            {
                amrex::ignore_unused(k);
                const amrex::IntVect iv_w(AMREX_D_DECL(i, j, k));
                if (g.near(iv_w) == 0) {
                    return;
                }
                // The sum of the planes, in the coordinates of iv_w
                amrex::Real m[3] = {0.0_rt, 0.0_rt, 0.0_rt};
                amrex::Real alpha = 0.0_rt;
                const int n = wall_contact_planes<Normal>(iv_w, dir, iwall, sgn, g.cand,
                    [&] (amrex::IntVect const& ref) {
                        alpha += g.pl(ref, 3);
                        for (int d = 0; d < 3; ++d) {
                            m[d] += g.pl(ref, d);
                        }
                        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                            alpha -= g.pl(ref, d) * static_cast<amrex::Real>(iv_w[d] - ref[d]);
                        }
                    });
                // Without a plane, or with planes that cancel out
                const amrex::Real norm = std::abs(m[0]) + std::abs(m[1]) + std::abs(m[2]);
                if (norm <= std::numeric_limits<amrex::Real>::epsilon() * static_cast<amrex::Real>(n)) {
                    return;
                }
                for (amrex::Real& md : m) {
                    md /= norm;
                }
                alpha /= norm;
                const amrex::IntVect anchor = ow.clamp(iv_w);
                const amrex::Real Fa = g.w(anchor);
                amrex::IntVect gc = iv_w;
                for (int l = g.gbx.smallEnd(dir); l <= g.gbx.bigEnd(dir); ++l) {
                    gc[dir] = l;
                    g.acc(gc, 0) += wall_ghost_fraction(gc, anchor, iv_w, m, alpha, Fa);
                    g.acc(gc, 1) += 1.0_rt;
                }
            });
        }
    }

    // The mean of the values of the planes extended to each ghost cell, or without one the value of
    // its nearest cell not beyond a wall, which is the wall cell of each wall it lies beyond
    const detail::OtherWalls all_walls{wall_lo, wall_hi, dlo, dhi, -1};
    {
        amrex::Vector<detail::MeanTag> tags;
        for (amrex::MFIter mfi(acc); mfi.isValid(); ++mfi) {
            tags.push_back({work.const_array(mfi), acc_ma[mfi.LocalIndex()], mfi.validbox()});
        }
        amrex::ParallelFor(tags, [=] AMREX_GPU_DEVICE (int i, int j, int k, detail::MeanTag const& g) noexcept
        {
            amrex::ignore_unused(k);
            const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
            const amrex::Real nplane = g.acc(iv, 1);
            if (nplane > 0.5_rt) {
                g.acc(iv, 0) /= nplane;
            } else {
                g.acc(iv, 0) = g.w(all_walls.clamp(iv));
            }
        });
    }

    // Back to the ghost cells of F beyond the walls, in the boxes of F whose ghost cells reach
    // beyond one. Each takes the copy as seen from it: the values of the boxes that cover its cells,
    // else its own.
    amrex::BoxList back_boxes;
    amrex::Vector<int> back_owners;
    amrex::Vector<int> back_index;
    for (int ib = 0; ib < a_F.boxArray().size(); ++ib) {
        const amrex::Box b = a_F.boxArray()[ib];
        if (all_walls.beyond(b.smallEnd() - nGhost) || all_walls.beyond(b.bigEnd() + nGhost)) {
            back_boxes.push_back(b);
            back_owners.push_back(a_F.DistributionMap()[ib]);
            back_index.push_back(ib);
        }
    }
    amrex::MultiFab back(amrex::BoxArray(std::move(back_boxes)),
                         amrex::DistributionMapping(std::move(back_owners)), 1, nGhost);
    amrex::Vector<amrex::Array4PairTag<amrex::Real>> in_tags;
    amrex::Vector<amrex::Array4PairTag<amrex::Real>> out_tags;
    for (amrex::MFIter mfi(back); mfi.isValid(); ++mfi) {
        const int fb = back_index[mfi.index()];
        in_tags.push_back({back.array(mfi), a_F.const_array(fb), back[mfi].box()});
        out_tags.push_back({a_F.array(fb), back.const_array(mfi), back[mfi].box()});
    }
    amrex::ParallelFor(in_tags, [=] AMREX_GPU_DEVICE (int i, int j, int k,
                                                      amrex::Array4PairTag<amrex::Real> const& g) noexcept
    {
        g.dfab(i, j, k) = g.sfab(i, j, k);
    });
    back.ParallelCopy(acc, 0, 0, 1, 0, nGhost, periodicity);
    amrex::ParallelFor(out_tags, [=] AMREX_GPU_DEVICE (int i, int j, int k,
                                                       amrex::Array4PairTag<amrex::Real> const& g) noexcept
    {
        amrex::ignore_unused(k);
        if (all_walls.beyond(amrex::IntVect(AMREX_D_DECL(i, j, k)))) {
            g.dfab(i, j, k) = g.sfab(i, j, k);
        }
    });
}

namespace {

// What the boundary condition of F on a domain face asks of FillWallGhosts: true to overwrite the
// ghost cells with the extended interface, false to keep the solver's values
bool extends_interface (int bc, int dir, int sgn)
{
    switch (bc) {
    case amrex::BCType::foextrap:
    case amrex::BCType::hoextrap:
    case amrex::BCType::hoextrapcc:
        return true;
    case amrex::BCType::int_dir:
    case amrex::BCType::ext_dir:
    case amrex::BCType::ext_dir_cc:
    case amrex::BCType::reflect_even:
    // Prescribed where the flow enters and extrapolated where it leaves, face by face: keep both, as
    // overwriting would also replace the volume fraction that flows in
    case amrex::BCType::direction_dependent:
        return false;
    default:
        amrex::Abort("PLIC::FillWallGhosts: boundary condition " + std::to_string(bc) + " on the "
                     + (sgn < 0 ? "low" : "high") + " side of direction " + std::to_string(dir)
                     + " is not valid for a volume fraction (expected int_dir, foextrap, hoextrap,"
                     + " hoextrapcc, ext_dir, ext_dir_cc, reflect_even or direction_dependent)");
        return false;
    }
}

// The boundary condition of a periodic direction is not used, and not checked
Walls wall_faces (amrex::BCRec const& bcrec, amrex::Geometry const& geom)
{
    Walls walls;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (geom.isPeriodic(dir)) {
            continue;
        }
        walls.lo[dir] = static_cast<int>(extends_interface(bcrec.lo(dir), dir, -1));
        walls.hi[dir] = static_cast<int>(extends_interface(bcrec.hi(dir), dir, 1));
    }
    return walls;
}

} // namespace

// One more than the stencil, for the first ghost layer of the planes. The advection schemes read no
// further than that.
int nGrowF (std::string const& reconstruction)
{
    return dispatch_normal(reconstruction, [] (auto normal) { return decltype(normal)::stencil + 1; });
}

void FitPlanes (amrex::MultiFab& plane, amrex::MultiFab const& F, int ngrow,
                std::string const& reconstruction)
{
    dispatch_normal(reconstruction, [&] (auto normal) {
        FitPlanesImpl<decltype(normal)>(plane, F, ngrow);
    });
}

void FillWallGhosts (amrex::MultiFab& F, amrex::Geometry const& geom, amrex::BCRec const& bcrec,
                     std::string const& reconstruction)
{
    const Walls walls = wall_faces(bcrec, geom);
    dispatch_normal(reconstruction, [&] (auto normal) {
        FillWallGhostsImpl<decltype(normal)>(F, geom, walls);
    });
}

} // namespace PLIC
/** @}*/

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

//! What each cell offers the search, cand, with its plane pl
struct CandTag {
    amrex::Array4<amrex::Real const> w;
    amrex::Array4<int const> fit_near;
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

//! The nearest cell of b to c
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::IntVect nearest_in (amrex::Box const& b, amrex::IntVect const& c) noexcept
{
    return amrex::max(amrex::min(c, b.bigEnd()), b.smallEnd());
}

//! Set the boxes of the tags to their values
template <typename T>
void SetBoxes (amrex::TagVector<amrex::Array4BoxValTag<T>> const& tv)
{
    amrex::ParallelFor(tv, [=] AMREX_GPU_DEVICE (int i, int j, int k, amrex::Array4BoxValTag<T> const& g) noexcept
    {
        g.dfab(i, j, k) = g.val;
    });
}

//! Copy the boxes of the tags from their sources
void CopyBoxes (amrex::TagVector<amrex::Array4PairTag<amrex::Real>> const& tv)
{
    amrex::ParallelFor(tv, [=] AMREX_GPU_DEVICE (int i, int j, int k, amrex::Array4PairTag<amrex::Real> const& g) noexcept
    {
        g.dfab(i, j, k) = g.sfab(i, j, k);
    });
}

//! Where SeparableFilter leaves its result: in the scratch p after an odd number of directions
constexpr bool filtered_in_scratch = (AMREX_SPACEDIM - 1) % 2 == 1;

// Replace each cell of the mask a of each tag, one direction t other than skip_dir at a time, by
// whether any of the cells within radius along t is nonzero. Cells outside its box are zero. The
// directions alternate between a and the scratch p (see filtered_in_scratch).
void SeparableFilter (amrex::TagVector<FilterTag> const& tv, int radius, int skip_dir)
{
    bool in_p = false;
    for (int t = 0; t < AMREX_SPACEDIM; ++t) {
        if (t == skip_dir) {
            continue;
        }
        const amrex::IntVect et = amrex::IntVect::TheDimensionVector(t);
        amrex::ParallelFor(tv, [=] AMREX_GPU_DEVICE (int i, int j, int k, FilterTag const& g) noexcept
        {
            amrex::ignore_unused(k);
            auto const& src = in_p ? g.p : g.a;
            auto const& dst = in_p ? g.a : g.p;
            const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
            int r = 0;
            for (int a = -radius; a <= radius && r == 0; ++a) {
                const amrex::IntVect c = iv + a * et;
                r = (g.bx.contains(c) && src(c) != 0) ? 1 : 0;
            }
            dst(i, j, k) = r;
        });
        in_p = !in_p;
    }
}

} // namespace detail

template <typename Normal>
void FitPlanesImpl (amrex::MultiFab& plane, amrex::MultiFab const& F, int ngrow)
{
    BL_PROFILE("PLIC::FitPlanes()");
    // In every direction: nGrow() is that of the first
    AMREX_ALWAYS_ASSERT(plane.nComp() == 4 && plane.nGrowVect().allGE(amrex::IntVect(ngrow)));
    AMREX_ALWAYS_ASSERT(F.nGrowVect().allGE(amrex::IntVect(ngrow + Normal::stencil)));
    AMREX_ALWAYS_ASSERT(plane.boxArray() == F.boxArray() &&
                        plane.DistributionMap() == F.DistributionMap());

    // The ghost layers beyond ngrow; the kernel writes every cell within it
    if (plane.nGrowVect() != amrex::IntVect(ngrow)) {
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
// a periodic face, so the walls do not depend on one another or on their order, and a ghost cell
// beyond several walls (an edge or a corner) takes the mean over those that extend a plane to it.
// The result does not depend on the box layout, except near a coarse/fine boundary: past it the
// copy holds no more ghost layers of data than F, and the search only takes planes whose stencil
// lies on data.
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
    // The boxes of the copy grow by nGhost in every direction
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(a_F.nGrowVect() == amrex::IntVect(nGhost),
                                     "PLIC::FillWallGhosts: F needs as many ghost cells in every direction");
    AMREX_ALWAYS_ASSERT(nGhost >= Normal::stencil + 1);

    // The ghost layers read: the search along the wall and the reconstruction stencil of the plane
    // it finds, or the ghost layers of F themselves if there are more
    const int nread = amrex::max(wall_nsearch + Normal::stencil, nGhost);
    const amrex::Box& domain = geom.Domain();
    const auto periodicity = geom.periodicity();

    // The boxes that touch a wall, extended through it, with the other boxes of F: all_boxes. The
    // copy holds them, cut to the layers the fill reads when they touch walls in one direction only
    // (the search reads Normal::stencil layers inside the wall). They stay disjoint, and keep the
    // owner of the box of F they extend, f_index.
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
            // The smallest box holding, for each wall the box touches, the layers beyond it and
            // Normal::stencil layers inside
            amrex::Box keep;
            const auto add = [&] (amrex::Box const& s) { keep = keep.ok() ? amrex::minBox(keep, s) : s; };
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                if (ext.smallEnd(dir) != b.smallEnd(dir)) {
                    amrex::Box s = ext;
                    s.setBig(dir, amrex::min(b.bigEnd(dir), domain.smallEnd(dir) + Normal::stencil));
                    add(s);
                }
                if (ext.bigEnd(dir) != b.bigEnd(dir)) {
                    amrex::Box s = ext;
                    s.setSmall(dir, amrex::max(b.smallEnd(dir), domain.bigEnd(dir) - Normal::stencil));
                    add(s);
                }
            }
            boxes.push_back(keep);
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

    // The copy holds, in each cell, the value of the box whose valid region covers it, or else the
    // ghost cell of F of its own box (beyond the wall or over a coarse/fine boundary): the boxes of
    // data, by local index. Other cells hold no data.
    amrex::MultiFab work(ext_ba, ext_dm, 1, nread);
    amrex::Vector<amrex::BoxList> data_boxes(work.local_size());
    {
        // The cells of the boxes of F and of their extensions, with their periodic images, are known
        // from the box lists, so no exchange is needed
        const std::vector<amrex::IntVect> shifts = periodicity.shiftIntVect(amrex::IntVect(nread));
        amrex::Vector<amrex::Array4PairTag<amrex::Real>> own_tags;
        amrex::Vector<amrex::Array4BoxValTag<amrex::Real>> zero_tags;
        std::vector<std::pair<int, amrex::Box>> isects;
        for (amrex::MFIter mfi(work); mfi.isValid(); ++mfi) {
            const int fb = f_index[mfi.index()];
            own_tags.push_back({work.array(mfi), a_F.const_array(fb), a_F[fb].box() & work[mfi].box()});
            amrex::BoxList& data = data_boxes[mfi.LocalIndex()];
            data.push_back(a_F[fb].box());
            for (amrex::IntVect const& shift : shifts) {
                all_ba.intersections(work[mfi].box() + shift, isects);
                for (auto const& is : isects) {
                    data.push_back(is.second - shift);
                }
            }
            // The cells without data, read but never used
            for (amrex::Box const& b : amrex::complementIn(work[mfi].box(), data)) {
                zero_tags.push_back({work.array(mfi), b, 0.0_rt});
            }
        }
        const amrex::TagVector<amrex::Array4BoxValTag<amrex::Real>> zero_tv(zero_tags);
        const amrex::TagVector<amrex::Array4PairTag<amrex::Real>> own_tv(own_tags);
        detail::SetBoxes(zero_tv);
        detail::CopyBoxes(own_tv);
    }
    work.ParallelCopy(a_F, 0, 0, 1, 0, nread, periodicity);
    work.FillBoundary(periodicity);

    // Further than any box reaches: a box grown by it has no bound on that side
    constexpr int unbounded = 1000000;

    // The domain without bounds in periodic directions, and with its ghost layers
    const amrex::Box periodic_domain = geom.growPeriodicDomain(unbounded);
    const amrex::Box extent = amrex::grow(periodic_domain, nGhost);

    // The cells beyond no wall, and the parts of a box beyond some wall
    amrex::Box nowall = domain;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (!is_wall(dir, -1)) {
            nowall.growLo(dir, unbounded);
        }
        if (!is_wall(dir, 1)) {
            nowall.growHi(dir, unbounded);
        }
    }
    const auto beyond_walls = [&] (amrex::Box const& b) { return amrex::boxDiff(b, nowall); };

    // On the valid cells of the copy beyond the walls, which hold every wall ghost cell: the sum of
    // the ghost values of the walls that extend a plane to the cell, and their number
    amrex::MultiFab acc(ext_ba, ext_dm, 2, 0);
    auto const& acc_ma = acc.arrays();
    amrex::Vector<detail::MeanTag> mean_tags;
    for (amrex::MFIter mfi(acc); mfi.isValid(); ++mfi) {
        for (amrex::Box const& b : beyond_walls(mfi.validbox())) {
            mean_tags.push_back({work.const_array(mfi), acc.array(mfi), b});
        }
    }
    const amrex::TagVector<detail::MeanTag> mean_tv(mean_tags);
    amrex::ParallelFor(mean_tv, [=] AMREX_GPU_DEVICE (int i, int j, int k, detail::MeanTag const& g) noexcept
    {
        g.acc(i, j, k, 0) = 0.0_rt;
        g.acc(i, j, k, 1) = 0.0_rt;
    });

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        for (int sgn = -1; sgn <= 1; sgn += 2) {
            if (!is_wall(dir, sgn)) {
                continue;
            }
            const int iwall = (sgn < 0) ? domain.smallEnd(dir) : domain.bigEnd(dir);
            // The cells beyond none of the other walls
            amrex::Box nowall_other = nowall;
            nowall_other.grow(dir, unbounded);
            // The cells the fill of this wall may read: in the domain, or beyond this wall or a
            // periodic face. Beyond another face F holds data only in the ghost cells of some
            // boxes, which would make the fill depend on the box layout.
            amrex::Box inside = periodic_domain;
            if (sgn < 0) {
                inside.growLo(dir, unbounded);
            } else {
                inside.growHi(dir, unbounded);
            }

            // The boxes of the copy that reach beyond this wall, with the regions of each step:
            //   scan  the wall cells within the search of the valid wall cells
            //   wbx   the valid wall cells, each with one plane shared by the ghost cells beyond it
            //   slab  the cells searched from them, up to Normal::stencil layers inside the wall
            //   gbx   the ghost cells beyond the valid wall cells
            struct Regions {
                int li;
                amrex::Box scan, wbx, slab, gbx;
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
                // The ghost layers of the copy hold the stencils of the slab
                AMREX_ASSERT(work[mfi].box().contains(amrex::grow(r.slab, Normal::stencil)));
                regions.push_back(r);
            }
            const auto nr = static_cast<int>(regions.size());
            auto const& w_ma = work.const_arrays();
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
                const amrex::TagVector<detail::NearTag> tv(tags);
                amrex::ParallelFor(tv, [=] AMREX_GPU_DEVICE (int i, int j, int k, detail::NearTag const& g) noexcept
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
                detail::SeparableFilter(amrex::TagVector<detail::FilterTag>(ftags), wall_nsearch, dir);
            }

            // The cells whose 3x3(x3) block, fit_near, and whose reconstruction stencil, fit_ok, read
            // only data inside (on a fine level, the layers past the ghost cells of F over a
            // coarse/fine boundary hold no data), and what each cell of the slab offers the search
            // with its plane, computed once for the overlapping searches
            amrex::Vector<amrex::IArrayBox> fit_near(nr);
            amrex::Vector<amrex::IArrayBox> fit_ok(nr);
            amrex::Vector<amrex::IArrayBox> cand(nr);
            amrex::Vector<amrex::FArrayBox> cand_plane(nr);
            // Each mask is one on the slab, and zero within its radius of a cell without data: the
            // cells the stencils of the slab read that lie in no box of data inside
            amrex::Vector<amrex::Array4BoxValTag<int>> ones;
            amrex::Vector<amrex::Array4BoxValTag<int>> zeros;
            amrex::Vector<detail::CandTag> cand_tags;
            amrex::Vector<detail::PlaneTag> plane_tags;
            for (int n = 0; n < nr; ++n) {
                Regions const& r = regions[n];
                fit_ok[n].resize(r.slab, 1, amrex::The_Async_Arena());
                cand[n].resize(r.slab, 1, amrex::The_Async_Arena());
                cand_plane[n].resize(r.slab, 4, amrex::The_Async_Arena());
                if (Normal::stencil > 1) {
                    fit_near[n].resize(r.slab, 1, amrex::The_Async_Arena());
                }
                const amrex::Box fitbx = amrex::grow(r.slab, Normal::stencil);
                amrex::BoxList have;
                for (amrex::Box const& b : data_boxes[r.li]) {
                    const amrex::Box x = b & fitbx & inside;
                    if (x.ok()) {
                        have.push_back(x);
                    }
                }
                const amrex::BoxList lack = amrex::complementIn(fitbx, have);
                const auto set_mask = [&] (amrex::IArrayBox& mask, int radius) {
                    ones.push_back({mask.array(), r.slab, 1});
                    for (amrex::Box const& b : lack) {
                        const amrex::Box z = amrex::grow(b, radius) & r.slab;
                        if (z.ok()) {
                            zeros.push_back({mask.array(), z, 0});
                        }
                    }
                };
                set_mask(fit_ok[n], Normal::stencil);
                if (Normal::stencil > 1) {
                    set_mask(fit_near[n], 1);
                }
                amrex::IArrayBox const& near_mask = (Normal::stencil > 1) ? fit_near[n] : fit_ok[n];
                auto const& Fw = w_ma[r.li];
                cand_tags.push_back({Fw, near_mask.const_array(), fit_ok[n].const_array(), cand[n].array(),
                                     cand_plane[n].array(), r.slab});
                amrex::IArrayBox const& near_n = detail::filtered_in_scratch ? near_scratch[n] : near[n];
                plane_tags.push_back({Fw, near_n.const_array(), cand[n].const_array(), cand_plane[n].const_array(),
                                      acc_ma[r.li], r.wbx, r.gbx});
            }
            // Kept until the end of the wall, so that the kernels queue without a wait between them
            const amrex::TagVector<amrex::Array4BoxValTag<int>> ones_tv(ones);
            const amrex::TagVector<amrex::Array4BoxValTag<int>> zeros_tv(zeros);
            const amrex::TagVector<detail::CandTag> cand_tv(cand_tags);
            const amrex::TagVector<detail::PlaneTag> plane_tv(plane_tags);
            detail::SetBoxes(ones_tv);
            detail::SetBoxes(zeros_tv);
            amrex::ParallelFor(cand_tv, [=] AMREX_GPU_DEVICE (int i, int j, int k, detail::CandTag const& g) noexcept
            {
                const amrex::IntVect c(AMREX_D_DECL(i, j, k));
                amrex::Real m[3] = {0.0_rt, 0.0_rt, 0.0_rt};
                amrex::Real alpha = 0.0_rt;
                g.cand(i, j, k) = wall_candidate<Normal>(c, dir, iwall, sgn, g.fit_near, g.fit_ok, g.w, m, alpha);
                for (int d = 0; d < 3; ++d) {
                    g.pl(i, j, k, d) = m[d];
                }
                g.pl(i, j, k, 3) = alpha;
            });

            // The plane of each wall cell, the mean of the equally near planes found, extended to the
            // ghost cells beyond it. Each ghost value keeps the offset of the wall cell from the plane
            // (exact for a uniform normal to the wall), from the nearest wall cell not beyond another
            // wall.
            amrex::ParallelFor(plane_tv, [=] AMREX_GPU_DEVICE (int i, int j, int k, detail::PlaneTag const& g) noexcept
            {
                amrex::ignore_unused(k);
                const amrex::IntVect iv_w(AMREX_D_DECL(i, j, k));
                if (g.near(iv_w) == 0) {
                    return;
                }
                // The sum of the planes, in the coordinates of iv_w
                amrex::Real m[3] = {0.0_rt, 0.0_rt, 0.0_rt};
                amrex::Real alpha = 0.0_rt;
                const int n = wall_contact_planes<Normal>(iv_w, dir, iwall, sgn, g.cand, g.pl, m, alpha);
                // Without a plane, or with planes that cancel out
                const amrex::Real norm = std::abs(m[0]) + std::abs(m[1]) + std::abs(m[2]);
                if (norm <= std::numeric_limits<amrex::Real>::epsilon() * static_cast<amrex::Real>(n)) {
                    return;
                }
                for (amrex::Real& md : m) {
                    md /= norm;
                }
                alpha /= norm;
                const amrex::IntVect anchor = detail::nearest_in(nowall_other, iv_w);
                const amrex::Real Fa = g.w(anchor);
                const amrex::Real Va = plane_volume_in(anchor, iv_w, m, alpha);
                amrex::IntVect gc = iv_w;
                for (int l = g.gbx.smallEnd(dir); l <= g.gbx.bigEnd(dir); ++l) {
                    gc[dir] = l;
                    g.acc(gc, 0) += wall_ghost_fraction(gc, iv_w, m, alpha, Fa, Va);
                    g.acc(gc, 1) += 1.0_rt;
                }
            });
        }
    }

    // The mean of the values of the planes extended to each ghost cell, or without one the value of
    // its nearest cell not beyond a wall, which is the wall cell of each wall it lies beyond
    amrex::ParallelFor(mean_tv, [=] AMREX_GPU_DEVICE (int i, int j, int k, detail::MeanTag const& g) noexcept
    {
        amrex::ignore_unused(k);
        const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
        const amrex::Real nplane = g.acc(iv, 1);
        if (nplane > 0.5_rt) {
            g.acc(iv, 0) /= nplane;
        } else {
            g.acc(iv, 0) = g.w(detail::nearest_in(nowall, iv));
        }
    });

    // Back to the ghost cells of F beyond the walls: the parts of the boxes of F with their ghost
    // cells beyond a wall, which may overlap. Each takes the copy as seen from its box: the values of
    // the boxes that cover its cells, else its own.
    amrex::BoxList back_boxes;
    amrex::Vector<int> back_owners;
    amrex::Vector<int> back_index;
    for (int ib = 0; ib < a_F.boxArray().size(); ++ib) {
        for (amrex::Box const& b : beyond_walls(amrex::grow(a_F.boxArray()[ib], nGhost))) {
            back_boxes.push_back(b);
            back_owners.push_back(a_F.DistributionMap()[ib]);
            back_index.push_back(ib);
        }
    }
    amrex::MultiFab back(amrex::BoxArray(std::move(back_boxes)),
                         amrex::DistributionMapping(std::move(back_owners)), 1, 0);
    amrex::Vector<amrex::Array4PairTag<amrex::Real>> in_tags;
    amrex::Vector<amrex::Array4PairTag<amrex::Real>> out_tags;
    for (amrex::MFIter mfi(back); mfi.isValid(); ++mfi) {
        const int fb = back_index[mfi.index()];
        in_tags.push_back({back.array(mfi), a_F.const_array(fb), mfi.validbox()});
        out_tags.push_back({a_F.array(fb), back.const_array(mfi), mfi.validbox()});
    }
    detail::CopyBoxes(amrex::TagVector<amrex::Array4PairTag<amrex::Real>>(in_tags));
    back.ParallelCopy(acc, 0, 0, 1, 0, 0, periodicity);
    detail::CopyBoxes(amrex::TagVector<amrex::Array4PairTag<amrex::Real>>(out_tags));
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

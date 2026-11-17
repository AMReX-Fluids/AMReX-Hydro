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

//! The domain faces whose ghost cells FillWallGhosts overwrites: extend[dir][0] low, extend[dir][1]
//! high. Never those of a periodic direction.
using ExtendFlags = amrex::Array<amrex::Array<bool, 2>, AMREX_SPACEDIM>;

// Components of the plane shared by the ghost cells beyond a wall cell: whether there is one, its
// m (3) and alpha, and the cell it belongs to (AMREX_SPACEDIM)

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

//! dst = src over bx
struct CopyTag {
    amrex::Array4<amrex::Real> dst;
    amrex::Array4<amrex::Real const> src;
    amrex::Box bx;
    [[nodiscard]] AMREX_GPU_HOST_DEVICE amrex::Box const& box () const noexcept { return bx; }
};

//! a mask m set to 1 over bx
struct MaskTag {
    amrex::Array4<int> m;
    amrex::Box bx;
    [[nodiscard]] AMREX_GPU_HOST_DEVICE amrex::Box const& box () const noexcept { return bx; }
};

//! a mask a filtered over its box bx, with p as scratch on the same box
struct FilterTag {
    amrex::Array4<int> a;
    amrex::Array4<int> p;
    amrex::Box bx;
    [[nodiscard]] AMREX_GPU_HOST_DEVICE amrex::Box const& box () const noexcept { return bx; }
};

//! the arrays of one box of the copy for a step of the fill of one wall, over bx: the copy w, an
//! int array written (out) and up to two read (in, in2), the planes of the wall cells pl, the box
//! of F with its ghost cells own, the sum of the ghost values of the walls acc, and the ghost cells
//! of the box beyond the wall gbx
struct WallTag {
    amrex::Array4<amrex::Real const> w;
    amrex::Array4<int> out;
    amrex::Array4<int const> in;
    amrex::Array4<int const> in2;
    amrex::Array4<amrex::Real> pl;
    amrex::Box bx;
    amrex::Box own;
    amrex::Array4<amrex::Real> acc;
    amrex::Box gbx;
    [[nodiscard]] AMREX_GPU_HOST_DEVICE amrex::Box const& box () const noexcept { return bx; }
};

//! The walls other than (dir, sgn): whether a cell lies beyond one, and the nearest cell that does not
struct OtherWalls {
    amrex::IntVect wall_lo, wall_hi, dlo, dhi;
    int dir, sgn;
    [[nodiscard]] AMREX_GPU_HOST_DEVICE bool beyond (amrex::IntVect const& c) const noexcept {
        bool b = false;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            b = b || (wall_lo[d] != 0 && c[d] < dlo[d] && !(d == dir && sgn < 0)) ||
                (wall_hi[d] != 0 && c[d] > dhi[d] && !(d == dir && sgn > 0));
        }
        return b;
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

// Directions are processed in turn, so edge and corner ghosts see the ghosts of earlier directions.
//
// A ghost value depends on cells further along the wall than F has ghost cells, so that F alone
// would give different values in neighbouring boxes. The fill therefore works on a copy of the boxes
// that touch a wall, extended through it: each wall ghost cell belongs to one box, and its ghost
// layers cover the whole search. Each wall reads only cells that are not beyond another wall, so
// the walls do not depend on one another or on their order, and a ghost cell beyond several walls
// (an edge or a corner) takes the mean over those that extend a plane to it. The result does not
// depend on the box
// layout, except near a coarse/fine boundary: past it the copy holds no more ghost layers of data
// than F, and the search only takes planes whose stencil lies on data.
template <typename Normal>
void FillWallGhostsImpl (amrex::MultiFab& a_F, amrex::Geometry const& geom, ExtendFlags const& walls)
{
    BL_PROFILE("PLIC::FillWallGhosts()");
    const int nGhost = a_F.nGrow();
    AMREX_ALWAYS_ASSERT(nGhost >= Normal::stencil + 1);

    amrex::IntVect wall_lo(0);
    amrex::IntVect wall_hi(0);
    bool has_wall = false;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        wall_lo[dir] = static_cast<int>(walls[dir][0]);
        wall_hi[dir] = static_cast<int>(walls[dir][1]);
        has_wall = has_wall || walls[dir][0] || walls[dir][1];
    }
    const auto is_wall = [&] (const int dir, const int sgn) {
        return ((sgn < 0) ? wall_lo : wall_hi)[dir] != 0;
    };
    if (!has_wall) {
        return;
    }

    // The ghost layers read: the search along the wall and the reconstruction stencil of the plane
    // it finds, or the ghost layers of F themselves if there are more
    const int nwall = nGhost;
    const int nread = amrex::max(wall_nsearch + Normal::stencil, nGhost);
    const amrex::Box& domain = geom.Domain();
    const auto periodicity = geom.periodicity();

    // The boxes that touch a wall, extended through it. They stay disjoint, and keep the owner of
    // the box of F they extend, f_index.
    amrex::BoxList boxes;
    amrex::Vector<int> owners;
    amrex::Vector<int> f_index;
    for (int ib = 0; ib < a_F.boxArray().size(); ++ib) {
        const amrex::Box b = a_F.boxArray()[ib];
        amrex::Box ext = b;
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            if (is_wall(dir, -1) && b.smallEnd(dir) == domain.smallEnd(dir)) {
                ext.growLo(dir, nwall);
            }
            if (is_wall(dir, 1) && b.bigEnd(dir) == domain.bigEnd(dir)) {
                ext.growHi(dir, nwall);
            }
        }
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
        amrex::Vector<detail::CopyTag> own_tags;
        amrex::Vector<detail::MaskTag> cover_tags;
        std::vector<std::pair<int, amrex::Box>> isects;
        for (amrex::MFIter mfi(work); mfi.isValid(); ++mfi) {
            const int fb = f_index[mfi.index()];
            own_tags.push_back({work.array(mfi), a_F.const_array(fb), a_F[fb].box()});
            for (amrex::IntVect const& shift : shifts) {
                for (amrex::BoxArray const* ba : {&a_F.boxArray(), &ext_ba}) {
                    ba->intersections(work[mfi].box() + shift, isects);
                    for (auto const& is : isects) {
                        cover_tags.push_back({covered.array(mfi), is.second - shift});
                    }
                }
            }
        }
        amrex::ParallelFor(own_tags, [=] AMREX_GPU_DEVICE (int i, int j, int k, detail::CopyTag const& g) noexcept
        {
            g.dst(i, j, k) = g.src(i, j, k);
        });
        amrex::ParallelFor(cover_tags, [=] AMREX_GPU_DEVICE (int i, int j, int k, detail::MaskTag const& g) noexcept
        {
            g.m(i, j, k) = 1;
        });
        amrex::Gpu::streamSynchronize();
    }
    work.ParallelCopy(a_F, 0, 0, 1, 0, nread, periodicity);
    work.FillBoundary(periodicity);

    // The domain and its ghost layers, without bounds in periodic directions
    amrex::Box extent = domain;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        extent.grow(dir, geom.isPeriodic(dir) ? 1000000 : nwall);
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
            const detail::OtherWalls ow{wall_lo, wall_hi, dlo, dhi, dir, sgn};

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

            // The wall cells near the interface: those with a mixed wall cell, not beyond another
            // wall, within the search. Elsewhere the ghost cells take the wall cell value.
            {
                amrex::Vector<detail::WallTag> tags;
                amrex::Vector<detail::FilterTag> ftags;
                for (int n = 0; n < nr; ++n) {
                    tags.push_back({w_ma[regions[n].li], near[n].array(), {}, {}, {}, regions[n].scan, {}, {}, {}});
                    ftags.push_back({near[n].array(), near_scratch[n].array(), regions[n].scan});
                }
                amrex::ParallelFor(tags, [=] AMREX_GPU_DEVICE (int i, int j, int k, detail::WallTag const& g) noexcept
                {
                    amrex::ignore_unused(k);
                    const amrex::IntVect c(AMREX_D_DECL(i, j, k));
                    g.out(i, j, k) = static_cast<int>(is_mixed(g.w(c)) && !ow.beyond(c));
                });
                detail::SeparableFilter(ftags, wall_nsearch, along, false);
            }
            int wall_contact = 0;
            {
                amrex::Gpu::Buffer<int> flag({0});
                int* p_flag = flag.data();
                amrex::Vector<detail::MaskTag> tags;
                for (int n = 0; n < nr; ++n) {
                    tags.push_back({near[n].array(), regions[n].scan});
                }
                amrex::ParallelFor(tags, [=] AMREX_GPU_DEVICE (int i, int j, int k, detail::MaskTag const& g) noexcept
                {
                    if (g.m(i, j, k) != 0) {
                        amrex::Gpu::Atomic::Max(p_flag, 1);
                    }
                });
                wall_contact = *flag.copyToHost();
                amrex::ParallelDescriptor::ReduceIntMax(wall_contact);
            }

            // Without a contact line the wall extends no plane
            if (wall_contact == 0) {
                continue;
            }

            // The cells whose reconstruction stencil reads only data and no cell beyond another wall
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
            amrex::Vector<detail::WallTag> fit_tags;
            amrex::Vector<detail::FilterTag> fit_ftags;
            amrex::Vector<detail::WallTag> cand_tags;
            amrex::Vector<detail::WallTag> plane_tags;
            for (int n = 0; n < nr; ++n) {
                Regions const& r = regions[n];
                auto const& Fw = w_ma[r.li];
                fit_tags.push_back({Fw, fit_ok[n].array(), c_ma[r.li], {}, {}, r.fitbx, r.own, {}, {}});
                fit_ftags.push_back({fit_ok[n].array(), fit_scratch[n].array(), r.fitbx});
                cand_tags.push_back({Fw, cand[n].array(), fit_ok[n].const_array(), {}, cand_plane[n].array(), r.slab,
                                     {}, {}, {}});
                plane_tags.push_back({Fw, {}, cand[n].const_array(), near[n].const_array(), cand_plane[n].array(),
                                      r.wbx, {}, acc_ma[r.li], r.gbx});
            }
            amrex::ParallelFor(fit_tags, [=] AMREX_GPU_DEVICE (int i, int j, int k, detail::WallTag const& g) noexcept
            {
                amrex::ignore_unused(k);
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                g.out(i, j, k) = ((g.in(i, j, k) != 0 || g.own.contains(iv)) && !ow.beyond(iv)) ? 1 : 0;
            });
            detail::SeparableFilter(fit_ftags, Normal::stencil, amrex::IntVect(1), true);
            amrex::ParallelFor(cand_tags, [=] AMREX_GPU_DEVICE (int i, int j, int k, detail::WallTag const& g) noexcept
            {
                const amrex::IntVect c(AMREX_D_DECL(i, j, k));
                amrex::Real m[3] = {0.0_rt, 0.0_rt, 0.0_rt};
                amrex::Real alpha = 0.0_rt;
                g.out(i, j, k) = wall_candidate<Normal>(c, dir, iwall, sgn, g.in, g.w, m, alpha);
                for (int d = 0; d < 3; ++d) {
                    g.pl(i, j, k, d) = m[d];
                }
                g.pl(i, j, k, 3) = alpha;
            });

            // The plane of each wall cell, the mean of the equally near planes found, extended to the
            // ghost cells beyond it. Each ghost value keeps the offset of the wall cell from the plane
            // (exact for a uniform normal to the wall), from the nearest wall cell not beyond another
            // wall.
            amrex::ParallelFor(plane_tags, [=] AMREX_GPU_DEVICE (int i, int j, int k, detail::WallTag const& g) noexcept
            {
                amrex::ignore_unused(k);
                const amrex::IntVect iv_w(AMREX_D_DECL(i, j, k));
                if (g.in2(iv_w) == 0) {
                    return;
                }
                // The sum of the planes, in the coordinates of iv_w
                amrex::Real m[3] = {0.0_rt, 0.0_rt, 0.0_rt};
                amrex::Real alpha = 0.0_rt;
                const int n = wall_contact_planes<Normal>(iv_w, dir, iwall, sgn, g.in, g.pl,
                    [&] (amrex::IntVect const& ref, amrex::Real const* mr, amrex::Real alpha_r) {
                        alpha += alpha_r;
                        for (int d = 0; d < 3; ++d) {
                            m[d] += mr[d];
                        }
                        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                            alpha -= mr[d] * static_cast<amrex::Real>(iv_w[d] - ref[d]);
                        }
                    });
                // Without a plane, or with planes that cancel out
                const amrex::Real norm = std::abs(m[0]) + std::abs(m[1]) + std::abs(m[2]);
                if (n == 0 || norm <= std::numeric_limits<amrex::Real>::epsilon() * static_cast<amrex::Real>(n)) {
                    return;
                }
                for (int d = 0; d < 3; ++d) {
                    m[d] /= norm;
                }
                alpha /= norm;
                const amrex::IntVect iv_a = ow.clamp(iv_w);
                const amrex::Real Fa = g.w(iv_a);
                const int glo = (sgn < 0) ? g.gbx.smallEnd(dir) : iwall + 1;
                const int ghi = (sgn < 0) ? iwall - 1 : g.gbx.bigEnd(dir);
                amrex::IntVect gc = iv_w;
                for (int l = glo; l <= ghi; ++l) {
                    gc[dir] = l;
                    g.acc(gc, 0) += wall_ghost_fraction(gc, iv_a, iv_w, m, alpha, Fa);
                    g.acc(gc, 1) += 1.0_rt;
                }
            });
            amrex::Gpu::streamSynchronize();
        }
    }

    // The mean of the planes extended to each ghost cell, or without one the value of its nearest
    // cell not beyond a wall, which is the wall cell of each wall it lies beyond
    {
        amrex::Vector<detail::WallTag> tags;
        for (amrex::MFIter mfi(acc); mfi.isValid(); ++mfi) {
            tags.push_back({work.const_array(mfi), {}, {}, {}, {}, mfi.validbox(), {}, acc_ma[mfi.LocalIndex()], {}});
        }
        amrex::ParallelFor(tags, [=] AMREX_GPU_DEVICE (int i, int j, int k, detail::WallTag const& g) noexcept
        {
            amrex::ignore_unused(k);
            const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
            const amrex::Real nplane = g.acc(iv, 1);
            if (nplane > 0.5_rt) {
                g.acc(iv, 0) /= nplane;
            } else {
                amrex::IntVect c = iv;
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    if (wall_lo[d] != 0) { c[d] = amrex::max(c[d], dlo[d]); }
                    if (wall_hi[d] != 0) { c[d] = amrex::min(c[d], dhi[d]); }
                }
                g.acc(iv, 0) = g.w(c);
            }
        });
        amrex::Gpu::streamSynchronize();
    }

    // Back to the ghost cells of F beyond the walls, in the boxes of F whose ghost cells reach
    // beyond one. Each takes the copy as seen from it: the values of the boxes that cover its cells,
    // else its own.
    amrex::BoxList back_boxes;
    amrex::Vector<int> back_owners;
    amrex::Vector<int> back_index;
    for (int ib = 0; ib < a_F.boxArray().size(); ++ib) {
        const amrex::Box b = a_F.boxArray()[ib];
        bool reaches = false;
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            reaches = reaches || (wall_lo[dir] != 0 && b.smallEnd(dir) - nGhost < dlo[dir]) ||
                      (wall_hi[dir] != 0 && b.bigEnd(dir) + nGhost > dhi[dir]);
        }
        if (reaches) {
            back_boxes.push_back(b);
            back_owners.push_back(a_F.DistributionMap()[ib]);
            back_index.push_back(ib);
        }
    }
    amrex::MultiFab back(amrex::BoxArray(std::move(back_boxes)),
                         amrex::DistributionMapping(std::move(back_owners)), 1, nGhost);
    amrex::Vector<detail::CopyTag> in_tags;
    amrex::Vector<detail::CopyTag> out_tags;
    for (amrex::MFIter mfi(back); mfi.isValid(); ++mfi) {
        const int fb = back_index[mfi.index()];
        in_tags.push_back({back.array(mfi), a_F.const_array(fb), back[mfi].box()});
        out_tags.push_back({a_F.array(fb), back.const_array(mfi), back[mfi].box()});
    }
    amrex::ParallelFor(in_tags, [=] AMREX_GPU_DEVICE (int i, int j, int k, detail::CopyTag const& g) noexcept
    {
        g.dst(i, j, k) = g.src(i, j, k);
    });
    amrex::Gpu::streamSynchronize();
    back.ParallelCopy(acc, 0, 0, 1, 0, nGhost, periodicity);
    amrex::ParallelFor(out_tags, [=] AMREX_GPU_DEVICE (int i, int j, int k, detail::CopyTag const& g) noexcept
    {
        amrex::ignore_unused(k);
        const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
        bool beyond = false;
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            beyond = beyond || (wall_lo[dir] != 0 && iv[dir] < dlo[dir]) ||
                     (wall_hi[dir] != 0 && iv[dir] > dhi[dir]);
        }
        if (beyond) {
            g.dst(i, j, k) = g.src(i, j, k);
        }
    });
    amrex::Gpu::streamSynchronize();
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
ExtendFlags extend_flags (amrex::BCRec const& bcrec, amrex::Geometry const& geom)
{
    ExtendFlags extend{};
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (geom.isPeriodic(dir)) {
            continue;
        }
        extend[dir][0] = extends_interface(bcrec.lo(dir), dir, -1);
        extend[dir][1] = extends_interface(bcrec.hi(dir), dir, 1);
    }
    return extend;
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
    const ExtendFlags walls = extend_flags(bcrec, geom);
    dispatch_normal(reconstruction, [&] (auto normal) {
        FillWallGhostsImpl<decltype(normal)>(F, geom, walls);
    });
}

} // namespace PLIC
/** @}*/

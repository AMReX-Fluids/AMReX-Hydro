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
#include <AMReX_LayoutData.H>

#include <string>
#include <utility>

using namespace amrex::literals;

namespace PLIC {

namespace {

//! The domain faces whose ghost cells FillWallGhosts overwrites: extend[dir][0] low, extend[dir][1]
//! high. Never those of a periodic direction.
using ExtendFlags = amrex::Array<amrex::Array<bool, 2>, AMREX_SPACEDIM>;

// Components of the plane shared by the ghost cells beyond a wall cell: whether there is one, its
// m (3) and alpha, and the cell it belongs to (AMREX_SPACEDIM)
constexpr int wp_found = 0;
constexpr int wp_m = 1;
constexpr int wp_alpha = 4;
constexpr int wp_ref = 5;
constexpr int wp_ncomp = 8;

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

// Replace each cell of fab, one direction t with dirs[t] != 0 at a time, by whether all (all = true)
// or any (all = false) of the cells within radius along t are nonzero. Cells outside the fab are
// zero.
void SeparableFilter (amrex::IArrayBox& fab, int radius, amrex::IntVect const& dirs, bool all)
{
    const amrex::Box bx = fab.box();
    amrex::IArrayBox prev(bx, 1, amrex::The_Async_Arena());
    auto const& A = fab.array();
    auto const& P = prev.array();
    for (int t = 0; t < AMREX_SPACEDIM; ++t) {
        if (dirs[t] == 0) {
            continue;
        }
        const amrex::IntVect et = amrex::IntVect::TheDimensionVector(t);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            P(i, j, k) = A(i, j, k);
        });
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            amrex::ignore_unused(k);
            const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
            // Stops at the first cell that decides: a zero for all, a nonzero for any
            int r = all ? 1 : 0;
            for (int a = -radius; a <= radius && (r != 0) == all; ++a) {
                const amrex::IntVect c = iv + a * et;
                r = (bx.contains(c) && P(c) != 0) ? 1 : 0;
            }
            A(i, j, k) = r;
        });
    }
}

// Mark, in near, the cells of the wall-adjacent layer that have a mixed cell of that layer within
// the search along the wall, in the boxes of work that reach beyond the wall (dir, sgn). A ghost cell
// beyond any other wall cell takes the value of the wall cell, without a search. Returns whether any
// cell is marked.
bool MarkWallContact (amrex::LayoutData<amrex::IArrayBox>& near, amrex::MultiFab const& work,
                      amrex::Box const& domain, amrex::Box const& extent, int dir, int sgn)
{
    const int iwall = (sgn < 0) ? domain.smallEnd(dir) : domain.bigEnd(dir);
    amrex::IntVect along(1);
    along[dir] = 0;
    int wall_contact = 0;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion()) reduction(max : wall_contact)
#endif
    for (amrex::MFIter mfi(work); mfi.isValid(); ++mfi) {
        if (!beyond_wall(mfi.validbox(), domain, dir, sgn).ok()) {
            continue;
        }
        amrex::Box scan = amrex::grow(mfi.validbox(), wall_nsearch);
        scan.setRange(dir, iwall);
        scan &= (work[mfi].box() & extent);
        near[mfi].resize(scan, 1);
        auto const& Fw = work.const_array(mfi);
        auto const& N = near[mfi].array();
        amrex::ParallelFor(scan, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            N(i, j, k) = static_cast<int>(is_mixed(Fw(i, j, k)));
        });
        SeparableFilter(near[mfi], wall_nsearch, along, false);
        wall_contact = amrex::max(wall_contact, near[mfi].template max<amrex::RunOn::Device>());
    }
    amrex::Gpu::streamSynchronize();
    amrex::ParallelDescriptor::ReduceIntMax(wall_contact);
    return wall_contact != 0;
}

// Set fit_ok, over slab grown by nstencil within fab_box, to 1 on the cells whose cube of half-width
// nstencil lies on cells that hold data: those covered by a box, and those of own (the box of F with
// its ghost cells). Outside fab_box there is no data.
void MarkFitOk (amrex::IArrayBox& fit_ok, amrex::Box const& slab, amrex::Box const& fab_box,
                amrex::Array4<int const> const& covered, amrex::Box const& own, int nstencil)
{
    const amrex::Box bx = amrex::grow(slab, nstencil) & fab_box;
    fit_ok.resize(bx, 1);
    auto const& O = fit_ok.array();
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        amrex::ignore_unused(k);
        const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
        O(i, j, k) = (covered(i, j, k) != 0 || own.contains(iv)) ? 1 : 0;
    });
    SeparableFilter(fit_ok, nstencil, amrex::IntVect(1), true);
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
// that touch a wall, extended through it: each wall ghost cell belongs to one box, its ghost layers
// cover the whole search, and they are exchanged after each wall with a contact line (without one,
// each box fills its ghost layers itself). The result does not depend on the box
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
    covered.setVal(1, 0);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(work); mfi.isValid(); ++mfi) {
        const int fb = f_index[mfi.index()];
        auto const& W = work.array(mfi);
        auto const& Fa = a_F.const_array(fb);
        amrex::ParallelFor(a_F[fb].box(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            W(i, j, k) = Fa(i, j, k);
        });
    }
    {
        amrex::iMultiFab f_cells(a_F.boxArray(), a_F.DistributionMap(), 1, 0);
        f_cells.setVal(1);
        work.ParallelCopy(a_F, 0, 0, 1, 0, nread, periodicity);
        covered.ParallelCopy(f_cells, 0, 0, 1, 0, nread, periodicity);
        work.FillBoundary(periodicity);
        covered.FillBoundary(periodicity);
    }

    // The domain and its ghost layers, without bounds in periodic directions
    amrex::Box extent = domain;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        extent.grow(dir, geom.isPeriodic(dir) ? 1000000 : nwall);
    }

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        for (int sgn = -1; sgn <= 1; sgn += 2) {
            if (!is_wall(dir, sgn)) {
                continue;
            }
            const int iwall = (sgn < 0) ? domain.smallEnd(dir) : domain.bigEnd(dir);

            // The wall cells near the interface: elsewhere the ghost cells take the wall cell value
            amrex::LayoutData<amrex::IArrayBox> near(ext_ba, ext_dm);
            const bool wall_contact = detail::MarkWallContact(near, work, domain, extent, dir, sgn);

            if (!wall_contact) {
                // Each ghost cell takes the value of its wall cell, which every box holds: filling
                // the covered cells of each fab gives what the exchange would, without it
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for (amrex::MFIter mfi(work); mfi.isValid(); ++mfi) {
                    const amrex::Box gbx = beyond_wall(work[mfi].box(), domain, dir, sgn);
                    if (!gbx.ok()) {
                        continue;
                    }
                    auto const& Fw = work.array(mfi);
                    auto const& C = covered.const_array(mfi);
                    amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        amrex::ignore_unused(k);
                        const amrex::IntVect g(AMREX_D_DECL(i, j, k));
                        if (C(g) != 0) {
                            amrex::IntVect iv_w = g;
                            iv_w[dir] = iwall;
                            Fw(g) = Fw(iv_w);
                        }
                    });
                }
                amrex::Gpu::streamSynchronize();
                continue;
            }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi(work); mfi.isValid(); ++mfi) {
                const amrex::Box gbx = beyond_wall(mfi.validbox(), domain, dir, sgn);
                if (!gbx.ok()) {
                    continue;
                }
                auto const& Fw = work.array(mfi);
                auto const& N = near[mfi].const_array();

                // One plane for each wall cell, shared by the ghost cells beyond it. The planes read
                // no ghost cells of this wall, so they are all found before any is written.
                amrex::Box wbx = mfi.validbox();
                wbx.setRange(dir, iwall);

                // What each cell searched offers, computed once for the overlapping searches
                amrex::Box slab = wbx;
                for (int t = 0; t < AMREX_SPACEDIM; ++t) {
                    if (t != dir) {
                        slab.grow(t, wall_nsearch);
                    }
                }
                slab.setRange(dir, (sgn < 0) ? iwall : iwall - Normal::stencil, Normal::stencil + 1);
                // The cells whose reconstruction stencil reads only data. On a fine level, the layers
                // past the ghost cells of F over a coarse/fine boundary hold none.
                amrex::IArrayBox fit_ok(amrex::The_Async_Arena());
                detail::MarkFitOk(fit_ok, slab, work[mfi].box(), covered.const_array(mfi),
                                  a_F[f_index[mfi.index()]].box(), Normal::stencil);
                auto const& OK = fit_ok.const_array();
                amrex::IArrayBox cand_fab(slab, 1, amrex::The_Async_Arena());
                auto const& CA = cand_fab.array();
                amrex::ParallelFor(slab, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    const amrex::IntVect c(AMREX_D_DECL(i, j, k));
                    CA(i, j, k) = wall_candidate<Normal>(c, dir, iwall, sgn, OK, Fw);
                });
                auto const& CAc = cand_fab.const_array();
                amrex::FArrayBox wall_plane(wbx, wp_ncomp, amrex::The_Async_Arena());
                auto const& PL = wall_plane.array();
                amrex::ParallelFor(wbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    amrex::ignore_unused(k);
                    const amrex::IntVect iv_w(AMREX_D_DECL(i, j, k));
                    amrex::IntVect ref = iv_w;
                    amrex::Real m[3] = {0.0_rt, 0.0_rt, 0.0_rt};
                    amrex::Real alpha = 0.0_rt;
                    const bool found = (N(iv_w) != 0) &&
                        wall_contact_plane<Normal>(iv_w, dir, iwall, sgn, CAc, Fw, ref, m, alpha);
                    PL(i, j, k, wp_found) = found ? 1.0_rt : 0.0_rt;
                    for (int n = 0; n < 3; ++n) {
                        PL(i, j, k, wp_m + n) = m[n];
                    }
                    PL(i, j, k, wp_alpha) = alpha;
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        PL(i, j, k, wp_ref + n) = static_cast<amrex::Real>(ref[n]);
                    }
                });

                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    amrex::ignore_unused(k);
                    const amrex::IntVect g(AMREX_D_DECL(i, j, k));
                    amrex::IntVect iv_w = g;
                    iv_w[dir] = iwall;
                    amrex::Real val = Fw(iv_w);
                    if (PL(iv_w, wp_found) > 0.5_rt) {
                        const amrex::Real m[3] = {PL(iv_w, wp_m), PL(iv_w, wp_m + 1), PL(iv_w, wp_m + 2)};
                        amrex::IntVect ref;
                        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                            ref[n] = static_cast<int>(PL(iv_w, wp_ref + n));
                        }
                        val = wall_ghost_fraction(g, iv_w, ref, m, PL(iv_w, wp_alpha), Fw(iv_w));
                    }
                    Fw(g) = val;
                });
            }
            amrex::Gpu::streamSynchronize();
            // The walls that follow read these ghost cells at edges and corners
            work.FillBoundary(periodicity);
        }
    }

    // Back to the ghost cells of F beyond the walls
    const amrex::IntVect dlo = domain.smallEnd();
    const amrex::IntVect dhi = domain.bigEnd();
    // The copy as seen from each box of F: the values of the boxes that cover its cells, else its own
    amrex::MultiFab back(a_F.boxArray(), a_F.DistributionMap(), 1, nGhost);
    amrex::MultiFab::Copy(back, a_F, 0, 0, 1, nGhost);
    back.ParallelCopy(work, 0, 0, 1, 0, nGhost, periodicity);
    auto const& f_ma = a_F.arrays();
    auto const& w_ma = back.const_arrays();
    amrex::ParallelFor(a_F, amrex::IntVect(nGhost),
                       [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept
    {
        amrex::ignore_unused(k);
        const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
        bool beyond = false;
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            beyond = beyond || (wall_lo[dir] != 0 && iv[dir] < dlo[dir]) ||
                     (wall_hi[dir] != 0 && iv[dir] > dhi[dir]);
        }
        if (beyond) {
            f_ma[b](i, j, k) = w_ma[b](i, j, k);
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

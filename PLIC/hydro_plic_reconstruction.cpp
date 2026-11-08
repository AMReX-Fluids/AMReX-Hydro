/**
 * \file hydro_plic_reconstruction.cpp
 * \addtogroup PLIC
 * @{
 */

#include <hydro_plic.H>
#include <hydro_plic_normals_K.H>

#include <AMReX_BC_TYPES.H>
#include <AMReX_IArrayBox.H>
#include <AMReX_LayoutData.H>
#include <AMReX_ParReduce.H>

#include <string>
#include <utility>

using namespace amrex::literals;

namespace PLIC {

namespace {

//! The domain faces whose ghost cells FillWallGhosts overwrites: extend[dir][0] low, extend[dir][1] high
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

// Mark, in near, the cells of the wall-adjacent layer that have a mixed cell of that layer within
// the search along the wall, in the boxes of work that reach beyond the wall (dir, sgn). A ghost cell
// beyond any other wall cell takes the value of the wall cell, without a search. The square of the
// search is dilated one direction at a time. Returns whether any cell is marked.
bool MarkWallContact (amrex::LayoutData<amrex::IArrayBox>& near, amrex::MultiFab const& work,
                      amrex::Box const& domain, amrex::Box const& extent, int dir, int sgn)
{
    const int iwall = (sgn < 0) ? domain.smallEnd(dir) : domain.bigEnd(dir);
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
        amrex::IArrayBox prev(scan, 1, amrex::The_Async_Arena());
        auto const& Fw = work.const_array(mfi);
        auto const& N = near[mfi].array();
        auto const& P = prev.array();
        amrex::ParallelFor(scan, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            N(i, j, k) = static_cast<int>(is_mixed(Fw(i, j, k)));
        });
        for (int t = 0; t < AMREX_SPACEDIM; ++t) {
            if (t == dir) {
                continue;
            }
            const amrex::IntVect et = amrex::IntVect::TheDimensionVector(t);
            amrex::ParallelFor(scan, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                P(i, j, k) = N(i, j, k);
            });
            amrex::ParallelFor(scan, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                amrex::ignore_unused(k);
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                int any = 0;
                for (int a = -wall_nsearch; a <= wall_nsearch && any == 0; ++a) {
                    const amrex::IntVect c = iv + a * et;
                    any = (scan.contains(c) && P(c) != 0) ? 1 : 0;
                }
                N(i, j, k) = any;
            });
        }
        wall_contact = amrex::max(wall_contact, near[mfi].template max<amrex::RunOn::Device>());
    }
    amrex::Gpu::streamSynchronize();
    amrex::ParallelDescriptor::ReduceIntMax(wall_contact);
    return wall_contact != 0;
}

} // namespace detail

template <typename Normal>
void FitPlanesImpl (amrex::MultiFab& plane, amrex::MultiFab const& F, int ngrow)
{
    BL_PROFILE("PLIC::FitPlanes()");
    AMREX_ALWAYS_ASSERT(plane.nComp() == 4 && plane.nGrow() >= ngrow);
    AMREX_ALWAYS_ASSERT(F.nGrow() >= ngrow + Normal::stencil);

    plane.setVal(0.0_rt);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(plane); mfi.isValid(); ++mfi) {
        const amrex::Box gbx = amrex::grow(mfi.validbox(), ngrow);
        auto const& Fa = F.const_array(mfi);
        auto const& pl = plane.array(mfi);
        amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (is_mixed(Fa(i, j, k))) {
                amrex::Real mx, my, mz, alpha;
                fit_plane<Normal>(i, j, k, Fa, mx, my, mz, alpha);
                pl(i, j, k, 0) = mx;
                pl(i, j, k, 1) = my;
                pl(i, j, k, 2) = mz;
                pl(i, j, k, 3) = alpha;
            }
        });
    }
    amrex::Gpu::streamSynchronize();
}

// Directions are processed in turn, so edge and corner ghosts see the ghosts of earlier directions.
//
// A ghost value depends on cells further along the wall than F has ghost cells, so that F alone
// would give different values in neighbouring boxes. The fill therefore works on a copy where the
// boxes extend through the walls: each wall ghost cell belongs to one box, its ghost layers cover
// the whole search, and they are exchanged between passes. The result does not depend on the box
// layout.
template <typename Normal>
void FillWallGhostsImpl (amrex::MultiFab& a_F, amrex::Geometry const& geom, ExtendFlags const& walls)
{
    BL_PROFILE("PLIC::FillWallGhosts()");
    const int nGhost = a_F.nGrow();
    AMREX_ALWAYS_ASSERT(nGhost >= Normal::stencil + 1);

    const auto is_wall = [&] (const int dir, const int sgn) {
        return !geom.isPeriodic(dir) && walls[dir][(sgn < 0) ? 0 : 1];
    };
    bool has_wall = false;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        has_wall = has_wall || is_wall(dir, -1) || is_wall(dir, 1);
    }
    if (!has_wall) {
        return;
    }

    // Passes beyond the first refine the ghost values at shallow contact angles, where the only
    // available planes read ghost cells
    constexpr int npass = 4;
    // The ghost layers read: the search along the wall, the reconstruction stencil of the plane it
    // finds, and the ghost layers of F themselves
    const int nwall = nGhost;
    constexpr int nread = wall_nsearch + Normal::stencil + 1;
    const amrex::Box& domain = geom.Domain();
    const auto periodicity = geom.periodicity();

    // Boxes extended through the walls they touch. They stay disjoint.
    amrex::BoxList boxes;
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
        boxes.push_back(ext);
    }
    amrex::BoxArray ext_ba(std::move(boxes));
    amrex::MultiFab work(ext_ba, a_F.DistributionMap(), 1, nread);
    amrex::MultiFab snap(ext_ba, a_F.DistributionMap(), 1, nread);
    work.setVal(0.0_rt);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(a_F); mfi.isValid(); ++mfi) {
        const amrex::Box bx = amrex::grow(mfi.validbox(), nGhost);
        work[mfi].template copy<amrex::RunOn::Device>(a_F[mfi], bx, 0, bx, 0, 1);
    }
    work.FillBoundary(periodicity);

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

            // The wall cells near the interface. The passes of this wall do not change them, and a
            // wall with none needs no more than one pass.
            amrex::LayoutData<amrex::IArrayBox> near(ext_ba, a_F.DistributionMap());
            const bool wall_contact = detail::MarkWallContact(near, work, domain, extent, dir, sgn);

            // Passes beyond the first can read ghost values, the same way each time: once one
            // changes nothing, neither would the next.
            const int wall_npass = wall_contact ? npass : 1;
            for (int pass = 0; pass < wall_npass; ++pass) {
                if (wall_contact) {
                    amrex::MultiFab::Copy(snap, work, 0, 0, 1, nread);
                }
                const bool allow_ghost = (pass > 0);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for (amrex::MFIter mfi(work); mfi.isValid(); ++mfi) {
                    const amrex::Box gbx = beyond_wall(mfi.validbox(), domain, dir, sgn);
                    if (!gbx.ok()) {
                        continue;
                    }
                    const amrex::Box limit = work[mfi].box() & extent;
                    auto const& Fw = work.array(mfi);
                    auto const& S = snap.const_array(mfi);
                    auto const& N = near[mfi].const_array();

                    // One plane for each wall cell, shared by the ghost cells beyond it
                    amrex::Box wbx = mfi.validbox();
                    wbx.setRange(dir, iwall);
                    amrex::FArrayBox wall_plane(wbx, wp_ncomp, amrex::The_Async_Arena());
                    auto const& PL = wall_plane.array();
                    if (wall_contact) {
                        amrex::ParallelFor(wbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            amrex::ignore_unused(k);
                            const amrex::IntVect iv_w(AMREX_D_DECL(i, j, k));
                            amrex::IntVect ref = iv_w;
                            amrex::Real m[3] = {0.0_rt, 0.0_rt, 0.0_rt};
                            amrex::Real alpha = 0.0_rt;
                            const bool found = (N(iv_w) != 0) &&
                                wall_contact_plane<Normal>(iv_w, dir, iwall, sgn, limit, S,
                                                           allow_ghost, ref, m, alpha);
                            PL(i, j, k, wp_found) = found ? 1.0_rt : 0.0_rt;
                            for (int n = 0; n < 3; ++n) {
                                PL(i, j, k, wp_m + n) = m[n];
                            }
                            PL(i, j, k, wp_alpha) = alpha;
                            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                                PL(i, j, k, wp_ref + n) = static_cast<amrex::Real>(ref[n]);
                            }
                        });
                    } else {
                        wall_plane.template setVal<amrex::RunOn::Device>(0.0_rt);
                    }

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
                if (pass > 0) {
                    // The pass only changes this wall's ghost cells
                    auto const& work_ma = work.const_arrays();
                    auto const& snap_ma = snap.const_arrays();
                    int changed = amrex::ParReduce(
                        amrex::TypeList<amrex::ReduceOpMax>{}, amrex::TypeList<int>{}, work,
                        amrex::IntVect(0),
                        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) -> amrex::GpuTuple<int>
                        {
                            return {static_cast<int>(work_ma[box_no](i, j, k) != snap_ma[box_no](i, j, k))};
                        });
                    amrex::ParallelDescriptor::ReduceIntMax(changed);
                    if (changed == 0) {
                        break;
                    }
                }
                work.FillBoundary(periodicity);
            }
        }
    }

    // Back to the ghost cells of F beyond the walls
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(a_F); mfi.isValid(); ++mfi) {
        const amrex::Box bx = amrex::grow(mfi.validbox(), nGhost);
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            for (int sgn = -1; sgn <= 1; sgn += 2) {
                if (!is_wall(dir, sgn)) {
                    continue;
                }
                const amrex::Box gbx = beyond_wall(bx, domain, dir, sgn);
                if (gbx.ok()) {
                    a_F[mfi].template copy<amrex::RunOn::Device>(work[mfi], gbx, 0, gbx, 0, 1);
                }
            }
        }
    }
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
                     + " ext_dir, reflect_even or direction_dependent)");
        return false;
    }
}

ExtendFlags extend_flags (amrex::BCRec const& bcrec)
{
    ExtendFlags extend{};
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
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
    const ExtendFlags walls = extend_flags(bcrec);
    dispatch_normal(reconstruction, [&] (auto normal) {
        FillWallGhostsImpl<decltype(normal)>(F, geom, walls);
    });
}

} // namespace PLIC
/** @}*/

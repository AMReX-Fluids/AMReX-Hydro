#ifndef PLIC_TEST_UTIL_H_
#define PLIC_TEST_UTIL_H_

// Helpers shared by the checks of the PLIC library

#include <hydro_plic.H>
#include <hydro_plic_normals_K.H>

#include <AMReX_BCRec.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Reduce.H>

#include <cmath>
#include <functional>
#include <limits>
#include <string>

namespace plic_test {

using namespace amrex;
using namespace amrex::literals;

// The groups of checks
void check_geometry ();
void check_walls ();
void check_face_density ();

//! Print the outcome of a check, and count it if it fails
void check (bool ok, std::string const& what);
int failures ();

#ifdef AMREX_USE_FLOAT
constexpr Real roundoff_tol = 2.0e-3_rt;
#else
constexpr Real roundoff_tol = 1.0e-10_rt;
#endif

//! A NaN counts as the largest error
AMREX_GPU_HOST_DEVICE inline Real nan_as_max (Real e) noexcept
{
    return amrex::isnan(e) ? std::numeric_limits<Real>::max() : e;
}

//! The plane n.x = d; the liquid is on the side n.x < d
struct Plane {
    GpuArray<Real, AMREX_SPACEDIM> n;
    Real d;
};

//! The plane of normal n through the center of the unit domain
Plane centered_plane (GpuArray<Real, AMREX_SPACEDIM> const& n);

//! The unit domain of n^dim cells, periodic where periodic[dir] != 0
Geometry unit_geometry (int n, IntVect const& periodic);

//! The same boundary condition on every face
BCRec uniform_bc (int bc);

//! F = the volume fraction of the plane in every cell of its fabs, ghost cells included. Beyond the
//! domain faces of the directions where extrapolate[dir] != 0, a ghost cell takes the value of the
//! nearest domain cell, which is what first-order extrapolation by a solver gives.
void fill_plane (MultiFab& F, Geometry const& geom, Plane const& p, IntVect const& extrapolate);

//! The cells of box beyond the low and the high domain face of direction dir (either may be empty)
Vector<Box> beyond_faces (Box const& box, Box const& domain, int dir);

//! The boxes of a fab to check, from its valid box and its box with ghost cells. Not a template
//! parameter: nvcc rejects a host lambda type as the template argument of a function that holds an
//! extended device lambda.
using Region = std::function<Vector<Box>(Box const& valid, Box const& fab)>;

//! Largest err(a, i, j, k) over the boxes region(valid box, fab box) of each fab a of mf, across
//! ranks, with a NaN as the largest; 0 if there are no cells
template <typename Err>
Real max_over (MultiFab const& mf, Region const& region, Err const& err)
{
    ReduceOps<ReduceOpMax> op;
    ReduceData<Real> data(op);
    bool any = false;
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        auto const& a = mf.const_array(mfi);
        for (Box const& b : region(mfi.validbox(), mf[mfi].box())) {
            // A reduction over an empty box launches no blocks, which fails on the GPU
            if (!b.ok()) {
                continue;
            }
            any = true;
            op.eval(b, data, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> GpuTuple<Real> {
                return {nan_as_max(err(a, i, j, k))};
            });
        }
    }
    Real result = any ? amrex::get<0>(data.value(op)) : 0.0_rt;
    ParallelDescriptor::ReduceRealMax(result);
    return result;
}

//! The valid box of a fab, as a Region
inline Vector<Box> valid_only (Box const& valid, Box const& /*fab*/)
{
    return {valid};
}

//! The domain with its ng ghost layers along dir, less the cells within shrink of the other faces
inline Box away_from_other_faces (Box const& domain, int dir, int shrink, int ng)
{
    Box b = grow(domain, -shrink);
    b.setRange(dir, domain.smallEnd(dir) - ng, domain.length(dir) + 2 * ng);
    return b;
}

//! Largest |a| over the boxes region(valid box, fab box) of mf
inline Real max_abs (MultiFab const& mf, Region const& region)
{
    return max_over(mf, region,
                    [=] AMREX_GPU_DEVICE (Array4<Real const> const& a, int i, int j, int k) noexcept -> Real {
        return std::abs(a(i, j, k));
    });
}

//! Largest |F - the volume fraction of the plane| over the boxes region(valid box, fab box)
inline Real max_plane_error (MultiFab const& F, Geometry const& geom, Plane const& p, Region const& region)
{
    const auto geomdata = geom.data();
    return max_over(F, region,
                    [=] AMREX_GPU_DEVICE (Array4<Real const> const& a, int i, int j, int k) noexcept -> Real {
        return std::abs(a(i, j, k) - PLIC::plane_fraction(i, j, k, geomdata, p.n.data(), p.d));
    });
}

} // namespace plic_test

#endif

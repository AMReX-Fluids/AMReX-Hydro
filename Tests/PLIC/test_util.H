#ifndef PLIC_TEST_UTIL_H_
#define PLIC_TEST_UTIL_H_

// Helpers shared by the checks of the PLIC library

#include <hydro_plic.H>
#include <hydro_plic_normals_K.H>

#include <AMReX_MultiFab.H>
#include <AMReX_Reduce.H>

#include <limits>
#include <string>

namespace plic_test {

using namespace amrex;
using namespace amrex::literals;

// The groups of checks
void check_geometry ();

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

} // namespace plic_test

#endif

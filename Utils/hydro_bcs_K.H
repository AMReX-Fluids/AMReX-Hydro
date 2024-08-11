/**
 * \file hydro_bcs_K.H
 *
 * This header file contains the inlined __host__ __device__ functions required for
 * the boundary condition routines for Godunov and MOL. It also contains function declarations
 * for controlling host functions.
 *
 */

/** \addtogroup Utilities
 * @{
 */


#ifndef HYDRO_BCS_K_H
#define HYDRO_BCS_K_H

#include <AMReX_Gpu.H>
#include <AMReX_BCRec.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Array.H>
#include <iomanip>
#include <hydro_constants.H>

/**
 * \namespace HydroBC
 *
 */
namespace HydroBC{
//
// Choose between single BC per domain face or position dependent BC array
//
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::BCRec getBC (const int i, const int j, const int k, const int n,
                    const amrex::Box& m_domain, const amrex::BCRec* bcr,
                    amrex::Array4<int const> const& bca)
{
    if ( !bca ) { return bcr[n]; }

    int lo[AMREX_SPACEDIM];
    int hi[AMREX_SPACEDIM];
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
        int index[] = {i,j,k};

        for ( int dim = 0; dim < AMREX_SPACEDIM; dim++){
            if (index[dim] < m_domain.smallEnd(dim)) { index[dim] = m_domain.smallEnd(dim); }
            if (index[dim] > m_domain.bigEnd(dim))   { index[dim] = m_domain.bigEnd(dim); }
            if (dim == dir) { index[dim] = m_domain.smallEnd(dim)-1; }
        }
        lo[dir] = bca.contains(index[0], index[1], index[2]) ?
            bca(index[0], index[1], index[2], n) : 0;
// FIXME?? if we don't contain (i,j,k) then it doesn't matter what the bc is there, because we don't touch it

        index[dir] = m_domain.bigEnd(dir)+1;
        hi[dir] = bca.contains(index[0], index[1], index[2]) ?
            bca(index[0], index[1], index[2], n) : 0;
    }

    amrex::BCRec bc(lo, hi);
    return bc;
}

/**
 *
 * <A ID="SetExtrapVelBCsLo"></A>
 *
 * \brief Boundary condition effects
 *
 *
 * For a detailed discussion see the AMReX-Hydro Guide
 * <A HREF="https://amrex-fluids.github.io/amrex-hydro/docs_html"></A>.
 *
 */

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void SetExtrapVelBCsLo (int edge_dir, int i, int j, int k, int n,
                        const amrex::Array4<const amrex::Real> &s,
                        amrex::Real &lo, amrex::Real &hi,
                        int bclo, int domlo)
{
    using namespace amrex;

    //
    // This is called in the extrapolation of normal velocity components only
    //

    GpuArray<int,3> iv = {i,j,k};

    if (iv[edge_dir] == domlo)
    {
        iv[edge_dir] = domlo-1;
        Real s_ext = s(iv[0],iv[1],iv[2],n);

        bool is_vel_in_dir = (n == edge_dir);

        if ( (bclo == BCType::ext_dir) ||
             (bclo == BCType::direction_dependent && s_ext > 0.0) )
        {
            if (is_vel_in_dir) {
                lo = s_ext;
                hi = s_ext;
            } else {
                lo = s_ext;
            }
        }
        else if (bclo == BCType::reflect_odd)
        {
            hi = Real(0.);
            lo = Real(0.);
        }
        else if ( bclo == BCType::foextrap || bclo == BCType::hoextrap ||
                  bclo == BCType::reflect_even ||
                 (bclo == BCType::direction_dependent && s_ext <= 0.0))
        {
            lo = hi;
        }
    }
}

/**
 *
 * <A ID="SetExtrapVelBCsHi"></A>
 *
 * \brief Boundary condition effects
 *
 *
 * For a detailed discussion see the AMReX-Hydro Guide
 * <A HREF="https://amrex-fluids.github.io/amrex-hydro/docs_html"></A>.
 *
 */

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void SetExtrapVelBCsHi (int edge_dir, int i, int j, int k, int n,
                   const amrex::Array4<const amrex::Real> &s,
                   amrex::Real &lo, amrex::Real &hi,
                   int bchi, int domhi)
{
    using namespace amrex;

    //
    // This is called in the extrapolation of normal velocity components only
    //

    GpuArray<int,3> iv = {i,j,k};

    if (iv[edge_dir] == domhi+1)
    {
        iv[edge_dir] = domhi+1;
        Real s_ext = s(iv[0],iv[1],iv[2],n);

        bool is_vel_in_dir = (n == edge_dir);

        if ( (bchi == BCType::ext_dir) ||
             (bchi == BCType::direction_dependent && s_ext < 0.0) )
        {
            if (is_vel_in_dir) {
                hi = s_ext;
                lo = s_ext;
            } else {
                hi = s_ext;
            }
        }
        else if (bchi == BCType::reflect_odd)
        {
            lo = Real(0.);
            hi = Real(0.);
        }
        else if ( bchi == BCType::foextrap || bchi == BCType::hoextrap ||
                  bchi == BCType::reflect_even ||
                 (bchi == BCType::direction_dependent && s_ext >= 0.0))
        {
            hi = lo;
        }
    }
}

/**
 *
 * <A ID="SetEdgeBCsLo"></A>
 *
 * \brief Boundary condition effects
 *
 *
 * For a detailed discussion see the AMReX-Hydro Guide
 * <A HREF="https://amrex-fluids.github.io/amrex-hydro/docs_html"></A>.
 *
 */

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void SetEdgeBCsLo (int edge_dir, int i, int j, int k, int n,
                   const amrex::Array4<const amrex::Real> &s,
                   amrex::Real &lo, amrex::Real &hi,
                   amrex::Real macvel,
                   int bclo, int domlo,
                   bool is_velocity )
{
    using namespace amrex;

    //
    // This is called for all quantities in the construction of edge states
    //

    GpuArray<int,3> iv = {i,j,k};

    if (iv[edge_dir] == domlo)
    {
        iv[edge_dir] = domlo-1;
        Real s_ext = s(iv[0],iv[1],iv[2],n);

        bool is_vel_in_dir = is_velocity && (n == edge_dir);

        if (bclo == BCType::direction_dependent && macvel > 0.0)
        {
            lo = s_ext;
            hi = s_ext;
        }
        else if (bclo == BCType::ext_dir)
        {
            lo = s_ext;

            // For turbulent inflow, there are times when the inflow face
            // may have a predicted outflowing velocity. Here, we preserve
            // the normal component of the Dirichlet BC, but allow the
            // tangential components to transport values from the interior.
            if (is_vel_in_dir) {
                hi = s_ext;
            }
        }
        else if (bclo == BCType::reflect_odd)
        {
            hi = Real(0.);
            lo = Real(0.);
        }
        else if ( bclo == BCType::foextrap || bclo == BCType::hoextrap ||
                  bclo == BCType::reflect_even ||
                 (bclo == BCType::direction_dependent && macvel <= 0.0))
        {
            lo = hi;
        }
    }
}

/**
 *
 * <A ID="SetEdgeBCsHi"></A>
 *
 * \brief Boundary condition effects
 *
 *
 * For a detailed discussion see the AMReX-Hydro Guide
 * <A HREF="https://amrex-fluids.github.io/amrex-hydro/docs_html"></A>.
 *
 */

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void SetEdgeBCsHi (int edge_dir, int i, int j, int k, int n,
                   const amrex::Array4<const amrex::Real> &s,
                   amrex::Real &lo, amrex::Real &hi,
                   amrex::Real macvel,
                   int bchi, int domhi,
                   bool is_velocity )
{
    using namespace amrex;

    //
    // This is called for all quantities in the construction of edge states
    //

    GpuArray<int,3> iv = {i,j,k};

    if (iv[edge_dir] == domhi+1)
    {
        iv[edge_dir] = domhi+1;
        Real s_ext = s(iv[0],iv[1],iv[2],n);

        bool is_vel_in_dir = is_velocity && (n == edge_dir);

        if (bchi == BCType::direction_dependent && macvel < 0.0)
        {
            hi = s_ext;
            lo = s_ext;
        }
        else if (bchi == BCType::ext_dir)
        {
            hi = s_ext;

            // For turbulent inflow, there are times when the inflow face
            // may have a predicted outflowing velocity. Here, we preserve
            // the normal component of the Dirichlet BC, but allow the
            // tangential components to transport values from the interior.
            if (is_vel_in_dir) {
                lo = s_ext;
            }
        }
        else if (bchi == BCType::reflect_odd)
        {
            lo = Real(0.);
            hi = Real(0.);
        }
        else if ( bchi == BCType::foextrap || bchi == BCType::hoextrap ||
                  bchi == BCType::reflect_even ||
                 (bchi == BCType::direction_dependent && macvel >= 0.0))
        {
            hi = lo;
        }
    }
}

} // namespace
#endif
/** @}*/

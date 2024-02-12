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
    // [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    //  const int lobc (const int dir, const int i, const int j, const int k, const int n,
    //                 const Box& m_domain, const BCRec* bcr,
    //                 Array4<int const> const& bca) ) const& noexcept
    // {
    //     if ( !bca ) { return bcr[n].lo(dir); }

    //     IntVect index {AMREX_D_DECL(i,j,k)};
    //     for ( int dim = 0; dim < AMREX_SPACEDIM; dim++){
    //         if (index[dim] < m_domain.smallEnd(dim)) { index[dim] = m_domain.smallEnd(dim); }
    //         if (index[dim] > m_domain.bigEnd(dim))   { index[dim] = m_domain.bigEnd(dim); }
    //         if (dim == dir) { index[dim] = m_domain.smallEnd(dim)-1; }
    //     }
    //     return bca(index, n);
    // };
    // [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    // const int hibc (const int dir, const int i, const int j, const int k, const int n,
    //                 const Box& m_domain, const BCRec* bcr,
    //                 Array4<int const> const& bca) ) const& noexcept
    // {
    //     if ( !bca ) { return bcr[n].hi(dir); }

    //     IntVect index {AMREX_D_DECL(i,j,k)};
    //     for ( int dim = 0; dim < AMREX_SPACEDIM; dim++){
    //         if (index[dim] < m_domain.smallEnd(dim)) { index[dim] = m_domain.smallEnd(dim); }
    //         if (index[dim] > m_domain.bigEnd(dim))   { index[dim] = m_domain.bigEnd(dim); }
    //         if (dim == dir) { index[dim] = m_domain.bigEnd(dim)+1; }
    //     }
    //     return bca(index, n);
    // };
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    const amrex::BCRec getBC (const int i, const int j, const int k, const int n,
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
    };

/**
 *
 * <A ID="SetXEdgeBCs"></A>
 *
 * \brief Boundary condition effects
 *
 *
 * For a detailed discussion see the AMReX-Hydro Guide
 * <A HREF="https://amrex-codes.github.io/amrex/hydro_html/EBMOL.html#sec:slopes">
 * section for EBMOL</A>.
 *
 */

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void SetXEdgeBCs( int i, int j, int k, int n,
                  const amrex::Array4<const amrex::Real> &s,
                  amrex::Real &lo, amrex::Real &hi,
                  int bclo, int domlo,
                  int bchi, int domhi,
                  bool is_velocity )
{

    using namespace amrex;

    // This function is for setting BCs ON domain faces only.
    AMREX_ASSERT( domlo <= i && i <= domhi+1 );

    if (i == domlo)
    {
        if (bclo == BCType::ext_dir)
        {
            lo = s(domlo-1, j, k, n);
            // For turbulent inflow, there are times when the inflow face
            // may have a predicted outflowing velocity. Here, we preserve
            // the normal component of the Dirichlet BC, but allow the
            // tangential components to transport values from the interior.
            if ( n==XVEL && is_velocity ) hi=lo;
        }
        else if ( bclo == BCType::foextrap || bclo == BCType::hoextrap ||
                  bclo == BCType::reflect_even )
        {
            lo = hi;
        }
        else if (bclo == BCType::reflect_odd)
        {
            hi = 0.;
            lo = hi;
        }
    }
    else if (i == domhi+1)
    {
        if (bchi == BCType::ext_dir)
        {
            hi = s(domhi+1, j, k, n);
            if ( n==XVEL && is_velocity ) lo = hi;
        }
        else if ( bchi == BCType::foextrap || bchi == BCType::hoextrap ||
                  bchi == BCType::reflect_even )
        {
            hi = lo;
        }
        else if (bchi == BCType::reflect_odd)
        {
            lo = 0.;
            hi = lo;
        }
    }
    else
        return;
}

/*
 * \fn SetYEdgeBCs
 */
/**
 * \brief Boundary condition effects
 *
 * For a detailed discussion see the AMReX-Hydro Guide section EBMOL.
 */

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void SetYEdgeBCs ( int i, int j, int k, int n,
                   const amrex::Array4<const amrex::Real> &s,
                   amrex::Real &lo, amrex::Real &hi,
                   int bclo, int domlo,
                   int bchi, int domhi,
                   bool is_velocity )
{
    using namespace amrex;

    // This function is for setting BCs ON domain faces only.
    AMREX_ASSERT( domlo <= j && j <= domhi+1 );

    if (j == domlo)
    {
        if (bclo == BCType::ext_dir)
        {
            lo = s(i, domlo-1, k, n);
            if ( n==YVEL && is_velocity ) hi=lo;
        }
        else if ( bclo == BCType::foextrap || bclo == BCType::hoextrap ||
                  bclo == BCType::reflect_even )
        {
            lo = hi;
        }
        else if (bclo == BCType::reflect_odd)
        {
            hi = 0.;
            lo = hi;
        }
    }
    else  if (j == domhi+1)
    {
        if (bchi == BCType::ext_dir)
        {
            hi = s(i, domhi+1, k, n);
            if ( n==YVEL && is_velocity ) lo=hi;
        }
        else if ( bchi == BCType::foextrap || bchi == BCType::hoextrap ||
                  bchi == BCType::reflect_even )
        {
            hi = lo;
        }
        else if(bchi == BCType::reflect_odd)
        {
            lo = 0.;
            hi = lo;
        }
    }
    else
        return;
}


//                                                                             //
//   =============================== 3D ONLY ===============================   //
//                                                                             //
#if (AMREX_SPACEDIM==3)

/**
 *
 * \brief Boundary condition effects
 *
 * For a detailed discussion see the AMReX-Hydro Guide.
 * For 3D only.
 */

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void SetZEdgeBCs ( int i, int j, int k, int n,
                   const amrex::Array4<const amrex::Real> &s,
                   amrex::Real &lo, amrex::Real &hi,
                   int bclo, int domlo,
                   int bchi, int domhi,
                   bool is_velocity )
{
    using namespace amrex;

    // This function is for setting BCs ON domain faces only.
    AMREX_ASSERT( domlo <= k && k <= domhi+1 );

    if (k == domlo)
    {
        if (bclo == BCType::ext_dir)
        {
            lo = s(i, j, domlo-1, n);
            if ( n==ZVEL && is_velocity ) hi=lo;
        }
        else if ( bclo == BCType::foextrap || bclo == BCType::hoextrap ||
                  bclo == BCType::reflect_even )
        {
            lo = hi;
        }
        else if(bclo == BCType::reflect_odd)
        {
            hi = 0.;
            lo = hi;
        }
    }
    else if (k == domhi+1)
    {
        if (bchi == BCType::ext_dir)
        {
            hi = s(i,j,domhi+1, n);
            if ( n==ZVEL && is_velocity ) lo=hi;
        }
        else if ( bchi == BCType::foextrap || bchi == BCType::hoextrap ||
                  bchi == BCType::reflect_even )
        {
            hi = lo;
        }
        else if(bchi == BCType::reflect_odd)
        {
            lo = 0.;
            hi = lo;
        }
    }
    else
        return;
}
#endif
} // namespace
#endif
/** @}*/

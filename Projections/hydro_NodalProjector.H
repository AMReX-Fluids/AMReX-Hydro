#ifndef NODAL_PROJECTION_H
#define NODAL_PROJECTION_H
#include <AMReX_Config.H>

#include <AMReX.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_MLMG.H>

//
//
// ***************************  DEFAULT MODE  ***************************
//
// Solves
//
//     div(sigma*grad(phi)) = div(vel) + S_nd + S_cc
//
// and  performs the projection
//
//     vel = vel - sigma * grad(phi)
//
// vel, sigma, and S_cc are cell-centered variables, while
// phi and S_nd are nodal-centered variables.
//
// ***************************  CUSTOM MODE   ***************************
//
// Solves
//
//     div(sigma*grad(phi)) = rhs
//
// and performs the projection
//
//     vel = vel - (sigma/alpha) * grad(phi)
//
// alpha is a cell-centered variable, while rhs is nodal-centered.
//
// In this mode, the user provides rhs and alpha
//
// By default alpha is assumed to be 1. Use setAlpha to change the default.
//
// Use setCustomRHS to provide a custom RHS, else "div(vel) + S_nd + S_cc"
// is used.
//
// Example: rhs = div(alpha*vel)
//
namespace Hydro {

class NodalProjector
{

public:

    NodalProjector ( const amrex::Vector<amrex::MultiFab*>&       a_vel,
                     const amrex::Vector<const amrex::MultiFab*>& a_sigma,
                     const amrex::Vector<amrex::Geometry>&                      a_geom,
                     const amrex::LPInfo&                                a_lpinfo,
                     const amrex::Vector<amrex::MultiFab*>&       a_S_cc = {},
                     const amrex::Vector<const amrex::MultiFab*>& a_S_nd = {} );

    NodalProjector ( const amrex::Vector<amrex::MultiFab*>&       a_vel,
                     const amrex::Vector<const amrex::MultiFab*>& a_sigma,
                     const amrex::Vector<amrex::Geometry>&                      a_geom,
                     const amrex::Vector<amrex::MultiFab*>&       a_S_cc = {},
                     const amrex::Vector<const amrex::MultiFab*>& a_S_nd = {} )
    : NodalProjector(a_vel, a_sigma, a_geom, amrex::LPInfo(), a_S_cc, a_S_nd ) {}


    NodalProjector ( const amrex::Vector<amrex::MultiFab*>&       a_vel,
                     const amrex::Real                            a_const_sigma,
                     const amrex::Vector<amrex::Geometry>&                      a_geom,
                     const amrex::LPInfo&                                a_lpinfo,
                     const amrex::Vector<amrex::MultiFab*>&       a_S_cc = {},
                     const amrex::Vector<const amrex::MultiFab*>& a_S_nd = {} );

    void project ( amrex::Real a_rtol = amrex::Real(1.0e-11), amrex::Real a_atol = amrex::Real(1.0e-14) );
    void project ( const amrex::Vector<amrex::MultiFab*>& a_phi, amrex::Real a_rtol = amrex::Real(1.0e-11),
                   amrex::Real a_atol = amrex::Real(1.0e-14) );

    amrex::Vector<       amrex::MultiFab* > getGradPhi      ()       {return GetVecOfPtrs(m_fluxes);}
    amrex::Vector< const amrex::MultiFab* > getGradPhiConst () const {return GetVecOfConstPtrs(m_fluxes);}
    amrex::Vector<       amrex::MultiFab* > getPhi          ()       {return GetVecOfPtrs(m_phi);}
    amrex::Vector< const amrex::MultiFab* > getPhiConst     () const {return GetVecOfConstPtrs(m_phi);}

    void computeRHS ( const amrex::Vector<amrex::MultiFab*>&       a_rhs,
                      const amrex::Vector<amrex::MultiFab*>&       a_vel,
                      const amrex::Vector<amrex::MultiFab*>&       a_S_cc = {},
                      const amrex::Vector<const amrex::MultiFab*>& a_S_nd = {} );

    void setAlpha     (const amrex::Vector<const amrex::MultiFab*> a_alpha)
        {m_alpha=a_alpha;m_has_alpha=true;}
    void setCustomRHS (const amrex::Vector<const amrex::MultiFab*> a_rhs);


    // Methods to set verbosity
    void setVerbose (int  v) noexcept { m_verbose = v; }


    // Set domain BC
    void setDomainBC ( std::array<amrex::LinOpBCType,AMREX_SPACEDIM> a_bc_lo,
                       std::array<amrex::LinOpBCType,AMREX_SPACEDIM> a_bc_hi );

    // Methods to get underlying objects
    // Use these to modify properties of MLMG and linear operator
    amrex::MLNodeLaplacian& getLinOp () noexcept { return *m_linop; }
    amrex::MLMG&            getMLMG  () noexcept { return *m_mlmg;  }

    // Methods to set MF for sync
    void setSyncResidualFine (amrex::MultiFab* a_sync_resid_fine) {m_sync_resid_fine=a_sync_resid_fine;}
    void setSyncResidualCrse (amrex::MultiFab* a_sync_resid_crse, amrex::IntVect a_ref_ratio, amrex::BoxArray a_fine_grids )
        {m_sync_resid_crse=a_sync_resid_crse; m_ref_ratio=a_ref_ratio; m_fine_grids=a_fine_grids;}

private:

    void setOptions ();
    void setCoarseBoundaryVelocityForSync ();
    void computeSyncResidual ();
    void averageDown (const amrex::Vector<amrex::MultiFab*> a_var);
    void define (amrex::LPInfo const& a_lpinfo);

    bool m_has_rhs   = false;
    bool m_has_alpha = false;
    bool m_need_bcs  = true;

    // Verbosity
    int  m_verbose        = 0;

    // amrex::Geometry
    amrex::Vector<amrex::Geometry>               m_geom;

    // EB factory if any
#ifdef AMREX_USE_EB
    amrex::Vector<amrex::EBFArrayBoxFactory const *>  m_ebfactory;
#endif

    // Cell-centered data
    amrex::Vector<amrex::MultiFab*>        m_vel;
    amrex::Vector<amrex::MultiFab>         m_fluxes;
    amrex::Vector<const amrex::MultiFab*>  m_alpha;
    amrex::Vector<amrex::MultiFab*>        m_S_cc;
    amrex::Vector<const amrex::MultiFab*>  m_sigma;
    amrex::Real                     m_const_sigma = 0.0;

    // Node-centered data
    amrex::Vector<amrex::MultiFab>         m_phi;
    amrex::Vector<amrex::MultiFab>         m_rhs;
    amrex::Vector<const amrex::MultiFab*>  m_S_nd;

    // Linear operator
    std::unique_ptr< amrex::MLNodeLaplacian > m_linop;

    // Solver
    std::unique_ptr< amrex::MLMG > m_mlmg;

     // Boundary conditions
    std::array<amrex::LinOpBCType,AMREX_SPACEDIM>  m_bc_lo;
    std::array<amrex::LinOpBCType,AMREX_SPACEDIM>  m_bc_hi;

    // Members for synchronization
    amrex::IntVect   m_ref_ratio;
    amrex::BoxArray  m_fine_grids;  // Grid at crse level + 1
    amrex::MultiFab* m_sync_resid_crse = nullptr;
    amrex::MultiFab* m_sync_resid_fine = nullptr;

    void printInfo ();
};

}

#endif

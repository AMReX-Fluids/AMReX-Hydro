#ifndef NODAL_PROJECTION_H
#define NODAL_PROJECTION_H
#include <AMReX_Config.H>

#include <AMReX.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_MLMG.H>
#ifdef AMREX_USE_HYPRE
#include <AMReX_HypreMLNodeLap.H>
#endif
#include <utility>

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

    NodalProjector (amrex::Vector<amrex::MultiFab*>       a_vel,
                    amrex::Vector<const amrex::MultiFab*> a_sigma,
                    amrex::Vector<amrex::Geometry>        a_geom,
                    const amrex::LPInfo&                  a_lpinfo,
                    amrex::Vector<amrex::MultiFab*>       a_S_cc = {},
                    amrex::Vector<const amrex::MultiFab*> a_S_nd = {});

    NodalProjector (amrex::Vector<amrex::MultiFab*>       a_vel,
                    amrex::Vector<const amrex::MultiFab*> a_sigma,
                    amrex::Vector<amrex::Geometry>        a_geom,
                    amrex::Vector<amrex::MultiFab*>       a_S_cc = {},
                    amrex::Vector<const amrex::MultiFab*> a_S_nd = {});

    NodalProjector (amrex::Vector<amrex::MultiFab*>       a_vel,
                    amrex::Real                           a_const_sigma,
                    amrex::Vector<amrex::Geometry>        a_geom,
                    const amrex::LPInfo&                  a_lpinfo,
                    amrex::Vector<amrex::MultiFab*>       a_S_cc = {},
                    amrex::Vector<const amrex::MultiFab*> a_S_nd = {} );

    void project ( amrex::Real a_rtol = amrex::Real(1.0e-11), amrex::Real a_atol = amrex::Real(1.0e-14) );
    void project ( const amrex::Vector<amrex::MultiFab*>& a_phi, amrex::Real a_rtol = amrex::Real(1.0e-11),
                   amrex::Real a_atol = amrex::Real(1.0e-14) );

    amrex::Vector<       amrex::MultiFab* > getGradPhi      ()       {return GetVecOfPtrs(m_fluxes);}
    amrex::Vector< const amrex::MultiFab* > getGradPhiConst () const {return GetVecOfConstPtrs(m_fluxes);}
    amrex::Vector<       amrex::MultiFab* > getPhi          ()       {return GetVecOfPtrs(m_phi);}
    amrex::Vector< const amrex::MultiFab* > getPhiConst     () const {return GetVecOfConstPtrs(m_phi);}
    amrex::Vector<       amrex::MultiFab* > getRHS          ()       {return GetVecOfPtrs(m_rhs);}
    amrex::Vector< const amrex::MultiFab* > getRHSConst     () const {return GetVecOfConstPtrs(m_rhs);}

    amrex::Vector<       amrex::MultiFab* > calcGradPhi     (const amrex::Vector<amrex::MultiFab*>& a_phi);

    void computeRHS ( const amrex::Vector<amrex::MultiFab*>&       a_rhs,
                      const amrex::Vector<amrex::MultiFab*>&       a_vel,
                      const amrex::Vector<amrex::MultiFab*>&       a_S_cc = {},
                      const amrex::Vector<const amrex::MultiFab*>& a_S_nd = {} );

    //! Set alpha for the custom mode.  alpha is stored as raw pointers, so it
    //! must outlive the projector, and it persists until it is reset.  Pass an
    //! empty vector to go back to the default alpha = 1.
    void setAlpha     (const amrex::Vector<const amrex::MultiFab*>& a_alpha)
        {
            if (a_alpha.empty()) { m_alpha.clear(); m_has_alpha = false; return; }
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(a_alpha.size() == m_phi.size(),
                "NodalProjector::setAlpha: alpha must have one entry per level");
            for (auto const* a : a_alpha) {
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(a != nullptr,
                    "NodalProjector::setAlpha: alpha entries must not be null");
            }
            m_alpha=a_alpha;m_has_alpha=true;
        }
    //! Set a custom RHS.  It persists until it is reset; pass an empty vector to
    //! go back to computing "div(vel) + S_nd + S_cc" in project().
    void setCustomRHS (amrex::Vector<const amrex::MultiFab*> a_rhs);

    //
    // Methods to point a projector that is kept across time steps at new
    // data.  The data must be on the grids (and EB) of the constructor; the
    // linear operator and, with use_mlhypre, the composite solver and its
    // setup are kept (see mlhypre_reuse_setup).  Like the constructor
    // arguments, the data are stored as raw pointers and must outlive the
    // projector.
    //

    //! Set the velocity to project (one entry per level)
    void setVelocity (amrex::Vector<amrex::MultiFab*> a_vel);
    //! Set sigma (one entry per level).  Only for a projector constructed
    //! with a variable sigma; a constant sigma cannot be changed.
    void setSigma (amrex::Vector<const amrex::MultiFab*> a_sigma);
    //! Set the cell- and node-centered source terms (one entry per level, or
    //! an empty vector for none)
    void setSourceTerms (amrex::Vector<amrex::MultiFab*>       a_S_cc,
                         amrex::Vector<const amrex::MultiFab*> a_S_nd);


    // Methods to set verbosity
    void setVerbose (int  v) noexcept
        { m_verbose = v;
          if (m_mlmg) { m_mlmg->setVerbose(m_verbose); } }


    // Set domain BC
    void setDomainBC ( std::array<amrex::LinOpBCType,AMREX_SPACEDIM> a_bc_lo,
                       std::array<amrex::LinOpBCType,AMREX_SPACEDIM> a_bc_hi );

    //! Throw std::runtime_error instead of aborting when the solve fails
    void setThrowException (bool t)
        { m_throw_exception = t;
          if (m_mlmg) { m_mlmg->setThrowException(t); }
#ifdef AMREX_USE_HYPRE
          if (m_hypremlnodelap) { m_hypremlnodelap->setThrowException(t); }
#endif
        }

    //! Iterations of the last projection (MLMG, or hypre with use_mlhypre)
    [[nodiscard]] int getNumIters () const noexcept { return m_num_iters; }

    // Methods to get underlying objects
    // Use these to modify properties of MLMG and linear operator
    amrex::MLNodeLaplacian& getLinOp () noexcept { return *m_linop; }
    amrex::MLMG&            getMLMG  () noexcept { return *m_mlmg;  }

    // Methods to set MF for sync
    void setSyncResidualFine (amrex::MultiFab* a_sync_resid_fine) {m_sync_resid_fine=a_sync_resid_fine;}
    void setSyncResidualCrse (amrex::MultiFab* a_sync_resid_crse, amrex::IntVect a_ref_ratio, amrex::BoxArray a_fine_grids )
        {m_sync_resid_crse=a_sync_resid_crse; m_ref_ratio=a_ref_ratio; m_fine_grids=std::move(a_fine_grids);}

private:

    void setOptions ();
    void queryMLHypre ();
    void setCoarseBoundaryVelocityForSync ();
    void computeSyncResidual ();
    void averageDown (amrex::Vector<amrex::MultiFab*> const& a_var);
    void define (amrex::LPInfo const& a_lpinfo);

    bool m_has_rhs   = false;
    bool m_has_alpha = false;
    bool m_need_bcs  = true;

    // Verbosity
    int  m_verbose        = 0;

    // amrex::Geometry
    amrex::Vector<amrex::Geometry>               m_geom;

    // EB factory if any
#if defined(AMREX_USE_EB) && !defined(HYDRO_NO_EB)
    amrex::Vector<amrex::EBFArrayBoxFactory const *>  m_ebfactory;
#endif

    // Cell-centered data
    amrex::Vector<amrex::MultiFab*>        m_vel;
    amrex::Vector<amrex::MultiFab>         m_fluxes;
    amrex::Vector<const amrex::MultiFab*>  m_alpha;
    amrex::Vector<amrex::MultiFab*>        m_S_cc;
    amrex::Vector<const amrex::MultiFab*>  m_sigma;
    int                                    m_sigma_ncomp = 0;
    amrex::Real                     m_const_sigma = amrex::Real(0.0);

    // Node-centered data
    amrex::Vector<amrex::MultiFab>         m_phi;
    amrex::Vector<amrex::MultiFab>         m_rhs;
    amrex::Vector<const amrex::MultiFab*>  m_S_nd;

    // Linear operator
    std::unique_ptr< amrex::MLNodeLaplacian > m_linop;

    // Solver
    std::unique_ptr< amrex::MLMG > m_mlmg;
    int m_maxiter   = 100;
    int m_num_iters = 0;

    bool m_throw_exception = false;

    // Composite (all AMR levels) solve with hypre instead of MLMG.  The
    // HypreMLNodeLap is kept across project() calls: with a variable sigma
    // its matrix is updated in place from the current sigma before each
    // solve, with a constant sigma it is assembled once.  setDomainBC
    // discards it (the matrix pattern depends on the BC types).
    bool m_use_mlhypre = false;
    std::string m_mlhypre_prefix = "hypre_mlnodelap";
    // Keep the solver setup (hypre) or the AlgMG hierarchy when sigma changes
    bool m_mlhypre_reuse_setup = false;
#ifdef AMREX_USE_HYPRE
    std::unique_ptr< amrex::HypreMLNodeLap > m_hypremlnodelap;
#endif

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

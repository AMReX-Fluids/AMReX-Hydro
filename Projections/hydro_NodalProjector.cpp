#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_MLMG.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>

#include <hydro_NodalProjector.H>

using namespace amrex;

namespace Hydro {

//
// Add defaults for boundary conditions???
//
NodalProjector::NodalProjector (amrex::Vector<amrex::MultiFab*>       a_vel,
                                amrex::Vector<const amrex::MultiFab*> a_sigma,
                                amrex::Vector<amrex::Geometry>        a_geom,
                                const LPInfo&                         a_lpinfo,
                                amrex::Vector<amrex::MultiFab*>       a_S_cc,
                                amrex::Vector<const amrex::MultiFab*> a_S_nd )
    : m_geom(std::move(a_geom)),
      m_vel(std::move(a_vel)),
      m_S_cc(std::move(a_S_cc)),
      m_sigma(std::move(a_sigma)),
      m_S_nd(std::move(a_S_nd))
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!m_sigma.empty(),
        "NodalProjector: sigma must not be empty; use the constant-sigma constructor instead");
    define(a_lpinfo);
}

NodalProjector::NodalProjector (amrex::Vector<amrex::MultiFab*>       a_vel,
                                amrex::Vector<const amrex::MultiFab*> a_sigma,
                                amrex::Vector<amrex::Geometry>        a_geom,
                                amrex::Vector<amrex::MultiFab*>       a_S_cc,
                                amrex::Vector<const amrex::MultiFab*> a_S_nd )
    : m_geom(std::move(a_geom)),
      m_vel(std::move(a_vel)),
      m_S_cc(std::move(a_S_cc)),
      m_sigma(std::move(a_sigma)),
      m_S_nd(std::move(a_S_nd))
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!m_sigma.empty(),
        "NodalProjector: sigma must not be empty; use the constant-sigma constructor instead");
    define(amrex::LPInfo());
}

NodalProjector::NodalProjector (amrex::Vector<amrex::MultiFab*>       a_vel,
                                amrex::Real                           a_const_sigma,
                                amrex::Vector<amrex::Geometry>        a_geom,
                                const LPInfo&                         a_lpinfo,
                                amrex::Vector<amrex::MultiFab*>       a_S_cc,
                                amrex::Vector<const amrex::MultiFab*> a_S_nd )
    : m_geom(std::move(a_geom)),
      m_vel(std::move(a_vel)),
      m_S_cc(std::move(a_S_cc)),
      m_const_sigma(a_const_sigma),
      m_S_nd(std::move(a_S_nd))
{
    // MLNodeLaplacian treats a constant sigma of 0 as the sentinel for
    // "variable sigma", which would leave the operator identically zero here.
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(a_const_sigma > amrex::Real(0.0),
        "NodalProjector: const sigma must be > 0");
    define(a_lpinfo);
}

void NodalProjector::define (LPInfo const& a_lpinfo)
{
    auto nlevs = int(m_vel.size());

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nlevs > 0,
        "NodalProjector: vel must have at least one level");
    for (auto const* v : m_vel) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(v != nullptr,
            "NodalProjector: vel entries must not be null");
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(int(m_geom.size()) == nlevs,
        "NodalProjector: geom must have the same number of levels as vel");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_sigma.empty() || int(m_sigma.size()) == nlevs,
        "NodalProjector: sigma must have the same number of levels as vel");
    for (auto const* sig : m_sigma) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(sig != nullptr,
            "NodalProjector: sigma entries must not be null");
    }

    Vector<BoxArray> ba(nlevs);
    Vector<DistributionMapping> dm(nlevs);
    for (int lev = 0; lev < nlevs; ++lev)
    {
        ba[lev] = m_vel[lev]->boxArray();
        dm[lev] = m_vel[lev]->DistributionMap();
    }

    // Resize member data
    m_phi.resize(nlevs);
    m_fluxes.resize(nlevs);
    m_rhs.resize(nlevs);


#if defined(AMREX_USE_EB) && !defined(HYDRO_NO_EB)
    bool has_eb = m_vel[0] -> hasEBFabFactory();
    for (int lev = 1; lev < nlevs; ++lev) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_vel[lev]->hasEBFabFactory() == has_eb,
            "NodalProjector: all levels must have the same EB-ness");
    }
    if (has_eb)
    {
        m_ebfactory.resize(nlevs,nullptr);
        for (int lev = 0; lev < nlevs; ++lev )
        {
            m_ebfactory[lev] = dynamic_cast<EBFArrayBoxFactory const*>(&(m_vel[lev]->Factory()));

            // Cell-centered data
            m_fluxes[lev].define(ba[lev], dm[lev], AMREX_SPACEDIM, 0, MFInfo(), m_vel[lev]->Factory());

            // Node-centered data
            auto tmp = ba[lev];
            const auto& ba_nd = tmp.surroundingNodes();
            m_phi[lev].define(ba_nd, dm[lev], 1, 1, MFInfo(), m_vel[lev]->Factory());
            m_rhs[lev].define(ba_nd, dm[lev], 1, 0, MFInfo(), m_vel[lev]->Factory());
        }
    }
    else
#endif
    {
        for (int lev = 0; lev < nlevs; ++lev )
        {
            // Cell-centered data
            m_fluxes[lev].define(ba[lev], dm[lev], AMREX_SPACEDIM, 0);

            // Node-centered data
            BoxArray tmp = ba[lev];
            const auto& ba_nd = tmp.surroundingNodes();
            m_phi[lev].define(ba_nd, dm[lev], 1, 1);
            m_rhs[lev].define(ba_nd, dm[lev], 1, 0);
        }
    }

    // Initialize all variables
    for (int lev(0); lev < m_phi.size(); ++lev)
    {
        m_phi[lev].setVal(Real(0.0));
        m_fluxes[lev].setVal(Real(0.0));
        m_rhs[lev].setVal(Real(0.0));
    }

    //
    // Setup linear operator
    //
    if (m_sigma.empty()) {
#if defined(AMREX_USE_EB) && !defined(HYDRO_NO_EB)
        m_linop = std::make_unique<MLNodeLaplacian>(m_geom, ba, dm, a_lpinfo, m_ebfactory,
                                                    m_const_sigma);
#else
        m_linop = std::make_unique<MLNodeLaplacian>(m_geom, ba, dm, a_lpinfo,
                                                    Vector<FabFactory<FArrayBox> const*>{},
                                                    m_const_sigma);
#endif
    } else {
#if defined(AMREX_USE_EB) && !defined(HYDRO_NO_EB)
        m_linop = std::make_unique<MLNodeLaplacian>(m_geom, ba, dm, a_lpinfo, m_ebfactory);
#else
        m_linop = std::make_unique<MLNodeLaplacian>(m_geom, ba, dm, a_lpinfo);
#endif
    }

    bool use_gauss_seidel = true;
    {
        ParmParse pp("nodal_proj");
        pp.query("use_gauss_seidel", use_gauss_seidel);
    }
    m_linop->setGaussSeidel(use_gauss_seidel);
    m_linop->setHarmonicAverage(false);

    //
    // Setup solver
    //
    m_mlmg = std::make_unique<MLMG>(*m_linop);

    setOptions();
}


//
// Set options by using default values and values read in input file
//
void
NodalProjector::setOptions ()
{
    // Default values
    int          bottom_verbose(0);
    int          maxiter(100);
    int          bottom_maxiter(100);
    Real         bottom_rtol(1.0e-4_rt);
    Real         bottom_atol(-1.0_rt);
    std::string  bottom_solver("bicgcg");

    int          num_pre_smooth (2);
    int          num_post_smooth(2);
    int          num_final_smooth(8);

    Real         normalization_threshold(-Real(1.));

    // Read from input file
    ParmParse pp("nodal_proj");
    pp.query( "verbose"       , m_verbose );
    pp.query( "bottom_verbose", bottom_verbose );
    pp.query( "maxiter"       , maxiter );
    pp.query( "bottom_maxiter", bottom_maxiter );
    pp.query( "bottom_rtol"   , bottom_rtol );
    pp.query( "bottom_atol"   , bottom_atol );
    pp.query( "bottom_solver" , bottom_solver );

    pp.query( "normalization_threshold" , normalization_threshold);

    pp.query( "num_pre_smooth"  , num_pre_smooth );
    pp.query( "num_post_smooth" , num_post_smooth );
    pp.query( "num_final_smooth" , num_final_smooth );

    // This is only used by the Krylov solvers but we pass it through the nodal operator
    //      if it is set here.  Otherwise we use the default set in AMReX_NodeLaplacian.H
    if (normalization_threshold > Real(0.))
        m_linop->setNormalizationThreshold(normalization_threshold);

    // Set default/input values
    m_mlmg->setVerbose(m_verbose);
    m_mlmg->setBottomVerbose(bottom_verbose);
    m_mlmg->setMaxIter(maxiter);
    m_mlmg->setBottomMaxIter(bottom_maxiter);
    m_mlmg->setBottomTolerance(bottom_rtol);
    m_mlmg->setBottomToleranceAbs(bottom_atol);

    m_mlmg->setPreSmooth(num_pre_smooth);
    m_mlmg->setPostSmooth(num_post_smooth);
    m_mlmg->setFinalSmooth(num_final_smooth);

    if (bottom_solver == "smoother")
    {
        m_mlmg->setBottomSolver(MLMG::BottomSolver::smoother);
    }
    else if (bottom_solver == "bicg")
    {
        m_mlmg->setBottomSolver(MLMG::BottomSolver::bicgstab);
    }
    else if (bottom_solver == "cg")
    {
        m_mlmg->setBottomSolver(MLMG::BottomSolver::cg);
    }
    else if (bottom_solver == "bicgcg")
    {
        m_mlmg->setBottomSolver(MLMG::BottomSolver::bicgcg);
    }
    else if (bottom_solver == "cgbicg")
    {
        m_mlmg->setBottomSolver(MLMG::BottomSolver::cgbicg);
    }
#ifdef AMREX_USE_HYPRE
    else if (bottom_solver == "hypre")
    {
        m_mlmg->setBottomSolver(MLMG::BottomSolver::hypre);
    }
#endif
}

void
NodalProjector::setDomainBC ( std::array<LinOpBCType,AMREX_SPACEDIM> a_bc_lo,
                              std::array<LinOpBCType,AMREX_SPACEDIM> a_bc_hi )
{
    m_bc_lo=a_bc_lo;
    m_bc_hi=a_bc_hi;
    m_linop->setDomainBC(m_bc_lo,m_bc_hi);
    m_need_bcs = false;
}

void
NodalProjector::setCustomRHS (amrex::Vector<const amrex::MultiFab*> a_rhs)
{
    // An empty vector resets to the default RHS, "div(vel) + S_nd + S_cc", so
    // that a projector reused across time steps does not keep a stale RHS.
    if (a_rhs.empty()) { m_has_rhs = false; return; }

    AMREX_ALWAYS_ASSERT(m_rhs.size()==a_rhs.size());

    for (int lev=0; lev < m_rhs.size(); ++lev)
    {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(a_rhs[lev] != nullptr,
            "NodalProjector::setCustomRHS: rhs entries must not be null");
        MultiFab::Copy(m_rhs[lev], *a_rhs[lev], 0, 0, 1, 0);
    }

    m_has_rhs = true;
}


void
NodalProjector::project ( Real a_rtol, Real a_atol )
{
    BL_PROFILE("NodalProjector::project");
    AMREX_ALWAYS_ASSERT(!m_need_bcs);

    if (m_verbose > 0)
        amrex::Print() << "Nodal Projection:" << '\n';

    //
    // Average fine grid velocity values down to the coarse grid
    // By doing this operation at this time we:
    //
    //   1) fill regions covered by a finer grid with some "valid" data, i.e. not NaNs
    //      for example. This is required internally by MLMG.
    //   2) make the velocity field ready for projection on the entirety of each level.
    //
    averageDown(m_vel);

    // Set matrix coefficients
    for (int lev = 0; lev < m_sigma.size(); ++lev)
    {
        m_linop -> setSigma(lev, *m_sigma[lev]);
    }

    // Compute RHS if necessary
    if (!m_has_rhs)
    {
        computeRHS( GetVecOfPtrs(m_rhs), m_vel, m_S_cc, m_S_nd );
    }

    // Print diagnostics
    if (m_verbose > 0)
    {
        amrex::Print() << " >> Before projection:" << '\n';
        printInfo();
        amrex::Print() << '\n';
    }

    // Solve
    // phi comes out already averaged-down and ready to be used by caller if needed
    m_mlmg -> solve( GetVecOfPtrs(m_phi), GetVecOfConstPtrs(m_rhs), a_rtol, a_atol );

    // Get fluxes -- fluxes = - sigma * grad(phi)
    m_mlmg -> getFluxes( GetVecOfPtrs(m_fluxes) );

    // Compute sync residual BEFORE performing projection
    computeSyncResidual();

    // Perform projection
    for (int lev(0); lev < m_phi.size(); ++lev)
    {
        if (m_has_alpha)
        {
            // fluxes -> fluxes/alpha = - ( sigma / alpha ) * grad(phi)
            for (int n = 0; n < AMREX_SPACEDIM; ++n)
            {
                MultiFab::Divide( m_fluxes[lev], *m_alpha[lev], 0, n, 1, 0 );
            }
        }

        //
        // vel = vel + fluxes = vel - ( sigma / alpha ) * grad(phi)
        //
        MultiFab::Add( *m_vel[lev], m_fluxes[lev], 0, 0, AMREX_SPACEDIM, 0);

        // set m_fluxes = grad(phi)
        m_linop->compGrad(lev,m_fluxes[lev],m_phi[lev]);
    }

    //
    // At this time, results are "correct" only on regions not covered by finer grids.
    // We average them down so that they are "correct" everywhere in each level.
    //
    averageDown(GetVecOfPtrs(m_fluxes));
    averageDown(m_vel);


    // Print diagnostics
    if ( (m_verbose > 0) && (!m_has_rhs))
    {
        computeRHS( GetVecOfPtrs(m_rhs), m_vel, m_S_cc, m_S_nd );
        amrex::Print() << " >> After projection:" << '\n';
        printInfo();
        amrex::Print() << '\n';
    }

}

void
NodalProjector::project ( const Vector<MultiFab*>& a_phi, Real a_rtol, Real a_atol )
{
    AMREX_ALWAYS_ASSERT(a_phi.size()==m_phi.size());

    // Only the valid nodes are needed: MLMG refills the ghost nodes of the
    // solution, and compGrad reads valid nodes only.  Copying with m_phi's
    // ghost width would run outside a_phi if the caller's phi has no ghosts.
    for (int lev=0; lev < m_phi.size(); ++lev )
    {
        MultiFab::Copy(m_phi[lev],*a_phi[lev],0,0,1,0);
    }

    project(a_rtol, a_atol);

    for (int lev=0; lev < m_phi.size(); ++lev )
    {
        MultiFab::Copy(*a_phi[lev],m_phi[lev],0,0,1,
                       amrex::min(a_phi[lev]->nGrowVect(), m_phi[lev].nGrowVect()));
    }
}

amrex::Vector<amrex::MultiFab*>
NodalProjector::calcGradPhi ( const Vector<MultiFab*>& a_phi )
{
    BL_PROFILE("NodalProjector::calcGradPhi");
    AMREX_ALWAYS_ASSERT(!m_need_bcs);

    // Copy input
    AMREX_ALWAYS_ASSERT(a_phi.size()==m_phi.size());
    for (int lev=0; lev < m_phi.size(); ++lev )
    {
        MultiFab::Copy(m_phi[lev],*a_phi[lev],0,0,1,0);
    }

    // Set coeffs involved with fluxes
    for (int lev = 0; lev < m_sigma.size(); ++lev)
    {
        m_linop -> setSigma(lev, *m_sigma[lev]);
    }

    // Perform projection
    for (int lev(0); lev < m_phi.size(); ++lev)
    {
        // set m_fluxes = grad(phi)
        m_linop->compGrad(lev,m_fluxes[lev],m_phi[lev]);
    }

    // Average down
    averageDown(GetVecOfPtrs(m_fluxes));

    // Output fluxes
    return GetVecOfPtrs(m_fluxes);
}


//
// Compute RHS: div(u) + S_nd + S_cc
//
void
NodalProjector::computeRHS (  const amrex::Vector<amrex::MultiFab*>&       a_rhs,
                              const amrex::Vector<amrex::MultiFab*>&       a_vel,
                              const amrex::Vector<amrex::MultiFab*>&       a_S_cc,
                              const amrex::Vector<const amrex::MultiFab*>& a_S_nd )
{
    AMREX_ALWAYS_ASSERT(!m_need_bcs); // This is needed to use linop
    AMREX_ALWAYS_ASSERT(a_rhs.size()==m_phi.size());
    AMREX_ALWAYS_ASSERT(a_vel.size()==m_phi.size());
    AMREX_ALWAYS_ASSERT(a_S_cc.empty() || (a_S_cc.size()==m_phi.size()) );
    AMREX_ALWAYS_ASSERT(a_S_nd.empty() || (a_S_nd.size()==m_phi.size()) );

    BL_PROFILE("NodalProjector::computeRHS");

    bool has_S_cc(!a_S_cc.empty());
    bool has_S_nd(!a_S_nd.empty());

    // Check the type of grids used
    for (int lev(0); lev < m_phi.size(); ++lev)
    {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(a_rhs[lev]->ixType().nodeCentered(), "a_rhs is not node centered");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(a_vel[lev]->ixType().cellCentered(), "a_vel is not cell centered");

        if (has_S_cc && (a_S_cc[lev]!=nullptr) )
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(a_S_cc[lev]->ixType().cellCentered(),"a_S_cc is not cell centered");

        if (has_S_nd && (a_S_nd[lev]!=nullptr) )
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(a_S_nd[lev]->ixType().nodeCentered(),"a_S_nd is not node centered");
    }

    m_linop -> compRHS( a_rhs, a_vel, a_S_nd, a_S_cc );
}


void
NodalProjector::printInfo ()
{
    for (int lev(0); lev < m_rhs.size(); ++lev)
    {
        amrex::Print() << "  * On lev " << lev
                       << " max(abs(rhs)) = "
                       << m_rhs[lev].norm0(0,0,false,true)
                       << '\n';
    }
}


void
NodalProjector::computeSyncResidual ()
{
    BL_PROFILE("NodalProjector::computeSyncResidual");

    if ( (m_sync_resid_fine != nullptr) || (m_sync_resid_crse != nullptr) )
    {
        int c_lev = 0;

        setCoarseBoundaryVelocityForSync();

        if (m_sync_resid_fine != nullptr)
        {
            MultiFab* rhptr = nullptr;
            if (!m_S_cc.empty())
                rhptr = m_S_cc[c_lev];
            m_linop->compSyncResidualFine(*m_sync_resid_fine, m_phi[c_lev], *m_vel[c_lev], rhptr);
        }

        if (m_sync_resid_crse != nullptr)
        {
            MultiFab* rhptr = nullptr;
            if (!m_S_cc.empty())
                rhptr = m_S_cc[c_lev];

            // this requires sigma to have 2 ghost cells (valid at -2)
            m_linop->compSyncResidualCoarse(*m_sync_resid_crse, m_phi[c_lev], *m_vel[c_lev],
                                           rhptr, m_fine_grids, m_ref_ratio);
        }


    }

}


// Set velocity in ghost cells to zero except for inflow
// 1) At non-inflow faces, the normal component of velocity will be completely zero'd
// 2) The normal velocity at corners -- even periodic corners -- just outside inflow faces will be zero'd
void
NodalProjector::setCoarseBoundaryVelocityForSync ()
{

    const BoxArray& grids = m_vel[0]->boxArray();
    const Box& domainBox  = m_geom[0].Domain();

    for (int idir=0; idir < AMREX_SPACEDIM; idir++)
    {
        if (m_bc_lo[idir] != LinOpBCType::inflow && m_bc_hi[idir] != LinOpBCType::inflow)
        {
            m_vel[0]->setBndry(Real(0.0), idir, 1);
        }
        else
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(*m_vel[0]); mfi.isValid(); ++mfi)
            {
                int i = mfi.index();
                FArrayBox& v_fab = (*m_vel[0])[mfi];

                const Box& reg = grids[i];
                const Box& bxg1 = amrex::grow(reg,1);
                BoxList bxlist(reg);

                if (m_bc_lo[idir] == LinOpBCType::inflow && reg.smallEnd(idir) == domainBox.smallEnd(idir))
                {
                    Box bx;                // bx is the region we *protect* from zero'ing
                    bx = amrex::adjCellLo(reg, idir);
                    bxlist.push_back(bx);
                }

                if (m_bc_hi[idir] == LinOpBCType::inflow && reg.bigEnd(idir) == domainBox.bigEnd(idir)) {
                    Box bx;                // bx is the region we *protect* from zero'ing
                    bx = amrex::adjCellHi(reg, idir);
                    bxlist.push_back(bx);
                }

                BoxList bxlist2 = amrex::complementIn(bxg1, bxlist);

                for (auto const& bit : bxlist2)
                {
                    Box ovlp = bit & v_fab.box();
                    if (ovlp.ok())
                    {
                        v_fab.setVal<RunOn::Device>(Real(0.0), ovlp, idir, 1);
                    }
                }
            }
        }
    }

}


void
NodalProjector::averageDown (const amrex::Vector<amrex::MultiFab*>& a_var)
{

    int f_lev = int(a_var.size())-1;
    int c_lev = 0;

    for (int lev = f_lev-1; lev >= c_lev; --lev)
    {
        IntVect rr   = m_geom[lev+1].Domain().size() / m_geom[lev].Domain().size();

#if defined(AMREX_USE_EB) && !defined(HYDRO_NO_EB)
        const auto* ebf = dynamic_cast<EBFArrayBoxFactory const*>(&(a_var[lev+1]->Factory()));
        if (ebf)
        {
            amrex::MultiFab volume(a_var[lev+1]->boxArray(),a_var[lev+1]->DistributionMap(),1,0);
            m_geom[lev+1].GetVolume(volume);

            EB_average_down(*a_var[lev+1], *a_var[lev], volume, ebf->getVolFrac(),
                            0, a_var[lev]->nComp(), rr);
        }
        else
#endif
        {
            average_down(*a_var[lev+1], *a_var[lev], m_geom[lev+1], m_geom[lev],
                         0, a_var[lev]->nComp(), rr);
        }
    } // lev
}

}

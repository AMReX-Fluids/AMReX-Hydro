#ifndef AMREX_MAC_PROJECTOR_H_
#define AMREX_MAC_PROJECTOR_H_
#include <AMReX_Config.H>

#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_MLABecLaplacian.H>

#ifdef AMREX_USE_EB
#include <AMReX_MLEBABecLap.H>
#endif

namespace Hydro {

class MacProjector
{
public:
    MacProjector(
        const amrex::Vector<amrex::Geometry>& a_geom,
#ifndef AMREX_USE_EB
        amrex::MLMG::Location a_umac_loc = amrex::MLMG::Location::FaceCenter,
        amrex::MLMG::Location a_beta_loc = amrex::MLMG::Location::FaceCenter,
        amrex::MLMG::Location a_phi_loc = amrex::MLMG::Location::CellCenter,
        amrex::MLMG::Location a_divu_loc = amrex::MLMG::Location::CellCenter
#else
        amrex::MLMG::Location a_umac_loc,
        amrex::MLMG::Location a_beta_loc,
        amrex::MLMG::Location a_phi_loc,
        amrex::MLMG::Location a_divu_loc = amrex::MLMG::Location::CellCenter
#endif
    );

    //
    // Constructors
    //
    MacProjector (const amrex::Vector<amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM> >& a_umac,
                  amrex::MLMG::Location a_umac_loc,
                  const amrex::Vector<amrex::Array<amrex::MultiFab const*,AMREX_SPACEDIM> >& a_beta,
                  amrex::MLMG::Location a_beta_loc,
                  amrex::MLMG::Location  a_phi_loc,
                  const amrex::Vector<amrex::Geometry>& a_geom,
                  const amrex::LPInfo& a_lpinfo,
                  const amrex::Vector<amrex::MultiFab const*>& a_divu,
                  amrex::MLMG::Location a_divu_loc,
                  const amrex::Vector<amrex::iMultiFab const*>& a_overset_mask = {});

    MacProjector (const amrex::Vector<amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM> >& a_umac,
                  amrex::MLMG::Location a_umac_loc,
                  const amrex::Vector<amrex::Array<amrex::MultiFab const*,AMREX_SPACEDIM> >& a_beta,
                  amrex::MLMG::Location a_beta_loc,
                  amrex::MLMG::Location  a_phi_loc,
                  const amrex::Vector<amrex::Geometry>& a_geom,
                  const amrex::LPInfo& a_lpinfo,
                  const amrex::Vector<amrex::MultiFab const*>& a_divu = {})
     : MacProjector(a_umac, a_umac_loc, a_beta, a_beta_loc, a_phi_loc,
                    a_geom, a_lpinfo, a_divu, amrex::MLMG::Location::CellCenter) {}

#ifndef AMREX_USE_EB
    MacProjector (const amrex::Vector<amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM> >& a_umac,
                  const amrex::Vector<amrex::Array<amrex::MultiFab const*,AMREX_SPACEDIM> >& a_beta,
                  const amrex::Vector<amrex::Geometry>& a_geom,
                  const amrex::LPInfo& a_lpinfo,
                  const amrex::Vector<amrex::MultiFab const*>& a_divu = {})
     : MacProjector(a_umac, amrex::MLMG::Location::FaceCenter,
                    a_beta, amrex::MLMG::Location::FaceCenter, amrex::MLMG::Location::CellCenter,
                    a_geom, a_lpinfo, a_divu, amrex::MLMG::Location::CellCenter) {}

    MacProjector (const amrex::Vector<amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM> >& a_umac,
                  const amrex::Vector<amrex::Array<amrex::MultiFab const*,AMREX_SPACEDIM> >& a_beta,
                  const amrex::Vector<amrex::Geometry>& a_geom,
                  const amrex::Vector<amrex::MultiFab const*>& a_divu = {})
     : MacProjector(a_umac, amrex::MLMG::Location::FaceCenter,
                    a_beta, amrex::MLMG::Location::FaceCenter, amrex::MLMG::Location::CellCenter,
                    a_geom, amrex::LPInfo(), a_divu, amrex::MLMG::Location::CellCenter) {}

    MacProjector (const amrex::Vector<amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM> >& a_umac,
                  const amrex::Real a_const_beta,
                  const amrex::Vector<amrex::Geometry>& a_geom,
                  const amrex::LPInfo& a_lpinfo,
                  const amrex::Vector<amrex::iMultiFab const*>& a_overset_mask = {},
                  const amrex::Vector<amrex::MultiFab const*>& a_divu = {});
#endif

    /** Initialize the underlying linear operator and MLMG instances
     */
    void initProjector (
        const amrex::LPInfo& a_lpinfo,
        const amrex::Vector<amrex::Array<amrex::MultiFab const*,AMREX_SPACEDIM> >& a_beta,
        const amrex::Vector<amrex::iMultiFab const*>& a_overset_mask = {});

    //! Update Bcoeffs for the linear operator
    void updateBeta (const amrex::Vector<amrex::Array<amrex::MultiFab const*,AMREX_SPACEDIM> >&);

#ifndef AMREX_USE_EB
    void initProjector (amrex::Vector<amrex::BoxArray> const& a_grids,
                        amrex::Vector<amrex::DistributionMapping> const& a_dmap,
                        const amrex::LPInfo& a_lpinfo, amrex::Real const a_const_beta,
                        const amrex::Vector<amrex::iMultiFab const*>& a_overset_mask = {});

    void updateBeta (amrex::Real a_const_beta);
#endif

    //! Set Umac before calling the projection step
    void setUMAC(const amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM> >&);

    //! Set div(U)
    void setDivU(const amrex::Vector<amrex::MultiFab const*>&);

#ifdef AMREX_USE_EB
    void setEBInflowVelocity (int amrlev, const amrex::MultiFab& eb_vel);
#endif

    //
    // Methods to set BCs and coarse/fine values
    //
    // These methods are wrappers of the linop methods of the same name
    // However, use of these is preferred to make sure operations
    // are performed in the correct order
    //
    void setDomainBC (const amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM>& lobc,
                      const amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM>& hibc);

    void setLevelBC  (int amrlev, const amrex::MultiFab* levelbcdata);

    void setCoarseFineBC (const amrex::MultiFab* crse, int crse_ratio)
        { m_linop->setCoarseFineBC(crse, crse_ratio);}

    //
    // Methods to perform projection
    //
    void project (const amrex::Vector<amrex::MultiFab*>& phi_in, amrex::Real reltol, amrex::Real atol);
    void project (amrex::Real reltol, amrex::Real atol);

    //
    // Get Fluxes.  DO NOT USE LinOp to get fluxes!!!
    //
    void getFluxes (const amrex::Vector<amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM> >& a_flux,
                    const amrex::Vector<amrex::MultiFab*>& a_sol, amrex::MLMG::Location a_loc) const;

    //
    // Setters and getters
    //
    void setVerbose            (int  v) noexcept
       { m_verbose = v;
         m_mlmg->setVerbose(m_verbose); }

    // Methods to get underlying objects
    // Use these to modify properties of MLMG and linear operator
    amrex::MLLinOp& getLinOp () noexcept { return *m_linop; }
    amrex::MLMG&    getMLMG  () noexcept { return *m_mlmg;  }

    bool needInitialization()  const noexcept { return m_needs_init; }

private:
    void setOptions ();

    void averageDownVelocity ();

    std::unique_ptr<amrex::MLPoisson> m_poisson;
    std::unique_ptr<amrex::MLABecLaplacian> m_abeclap;
#ifdef AMREX_USE_EB
    std::unique_ptr<amrex::MLEBABecLap> m_eb_abeclap;
    amrex::Vector<amrex::EBFArrayBoxFactory const*> m_eb_factory;
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > m_eb_vel;
#endif
    amrex::MLLinOp* m_linop = nullptr;

    amrex::Real m_const_beta = 0.;

    std::unique_ptr<amrex::MLMG> m_mlmg;

    amrex::Vector<amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM> > m_umac;
    amrex::Vector<amrex::MultiFab> m_rhs;
    amrex::Vector<amrex::MultiFab> m_phi;
    amrex::Vector<amrex::MultiFab> m_divu;
    amrex::Vector<amrex::Array<amrex::MultiFab,AMREX_SPACEDIM> > m_fluxes;

    amrex::Vector<amrex::Geometry> m_geom;

    int m_verbose = 0;

    bool m_needs_domain_bcs = true;
    amrex::Vector<int> m_needs_level_bcs;

    // Location of umac -- face center vs face centroid
    amrex::MLMG::Location m_umac_loc;

    amrex::MLMG::Location m_beta_loc;

    amrex::MLMG::Location m_phi_loc;

    // Location of divu (RHS -- optional) -- cell center vs cell centroid
    amrex::MLMG::Location m_divu_loc;

    bool m_needs_init = true;
};

}

#endif

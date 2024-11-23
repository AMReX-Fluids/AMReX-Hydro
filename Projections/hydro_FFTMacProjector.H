#ifndef AMREX_FFT_MAC_PROJECTOR_H_
#define AMREX_FFT_MAC_PROJECTOR_H_

#include <AMReX_FFT_Poisson.H>
#include <AMReX_LO_BCTYPES.H>

namespace Hydro {

class FFTMacProjector
{
public:

    FFTMacProjector (amrex::Geometry const& geom,
                     amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM> const& lobc,
                     amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM> const& hibc);

    void setUMAC (amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM> const& umac);

    void project ();

private:
    amrex::Geometry m_geom;
    amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM> m_lobc;
    amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM> m_hibc;
    amrex::FFT::Poisson<amrex::MultiFab> m_fft;
    amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM> m_umac;
};

}
#endif

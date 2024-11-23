#include <hydro_FFTMacProjector.H>
#include <AMReX_MultiFabUtil.H>

using namespace amrex;

namespace Hydro {

namespace {
    Array<std::pair<FFT::Boundary,FFT::Boundary>,AMREX_SPACEDIM>
    make_fft_bc (amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM> const& lobc,
                 amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM> const& hibc)
    {
        Array<std::pair<FFT::Boundary,FFT::Boundary>,AMREX_SPACEDIM> fft_bc;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            if ((lobc[idim] == LinOpBCType::Periodic &&
                 hibc[idim] != LinOpBCType::Periodic) ||
                (lobc[idim] != LinOpBCType::Periodic &&
                 hibc[idim] == LinOpBCType::Periodic))
            {
                amrex::Abort("FFTMacProjector: wrong periodic BC type");
            }
            for (int face = 0; face < 2; ++face) {
                auto& lbc = (face == 0) ? lobc[idim] : hibc[idim];
                auto& fbc = (face == 0) ? fft_bc[idim].first : fft_bc[idim].second;
                switch (lbc)
                {
                case LinOpBCType::Periodic:
                    fbc = FFT::Boundary::periodic;
                    break;
                case LinOpBCType::Dirichlet:
                    fbc = FFT::Boundary::odd;
                    break;
                case LinOpBCType::Neumann:
                    fbc = FFT::Boundary::even;
                    break;
                default:
                    amrex::Abort("FFTMacProjector: only Periodic, Dirichlet & Neumann are supported");
                }
            }
        }
        return fft_bc;
    }
}

FFTMacProjector::FFTMacProjector (
    amrex::Geometry const& geom,
    amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM> const& lobc,
    amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM> const& hibc)
    : m_geom(geom),
      m_lobc(lobc),
      m_hibc(hibc),
      m_fft(geom, make_fft_bc(lobc,hibc))
{}

void FFTMacProjector::setUMAC (Array<MultiFab*,AMREX_SPACEDIM> const& umac)
{
    m_umac = umac;
}

void FFTMacProjector::project ()
{
    MultiFab rhs(amrex::convert(m_umac[0]->boxArray(), IntVect(0)),
                 m_umac[0]->DistributionMap(), 1, 0);
    amrex::computeDivergence(rhs, GetArrOfConstPtrs(m_umac), m_geom);

    AMREX_ALWAYS_ASSERT(m_geom.Domain().numPts() == rhs.boxArray().numPts());

    bool has_dirichlet = false;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        has_dirichlet = has_dirichlet ||
            m_lobc[idim] == LinOpBCType::Dirichlet ||
            m_hibc[idim] == LinOpBCType::Dirichlet;
    }
    if (! has_dirichlet) {
        auto rhosum = rhs.sum(0);
        rhs.plus(-rhosum/m_geom.Domain().d_numPts(), 0, 1);
    }

    MultiFab phi(rhs.boxArray(), rhs.DistributionMap(), 1, 1);
    m_fft.solve(phi, rhs);
    phi.FillBoundary(m_geom.periodicity());

    auto const& phima = phi.const_arrays();
    Box const& domain = m_geom.growPeriodicDomain(1);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        auto const& uma = m_umac[idim]->arrays();
        int domlo = domain.smallEnd(idim);
        int domhi = domain.bigEnd(idim)+1;
        auto const bclo = m_lobc[idim];
        auto const bchi = m_hibc[idim];
        Real const dxinv = m_geom.InvCellSize(idim);
        ParallelFor(*m_umac[idim], [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
        {
            IntVect iv(AMREX_D_DECL(i,j,k));
            IntVect miv = iv - IntVect::TheDimensionVector(idim);
            auto const& u = uma[b];
            auto const& p = phima[b];
            if (iv[idim] == domlo) {
                if (bclo == LinOpBCType::Dirichlet) {
                    u(i,j,k) -= Real(2.0)*dxinv*p(iv);
                }
            } else if (iv[idim] == domhi) {
                if (bchi == LinOpBCType::Dirichlet) {
                    u(i,j,k) += Real(2.0)*dxinv*p(miv);
                }
            } else {
                u(i,j,k) -= dxinv*(p(iv)-p(miv));
            }
        });
    }
    Gpu::streamSynchronize();
}

}

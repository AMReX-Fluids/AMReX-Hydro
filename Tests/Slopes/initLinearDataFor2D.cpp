#include "MyTest.H"
#include <AMReX_Config.H>
#include <AMReX_EB2.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_MultiFabUtil.H>

using namespace amrex;

#if (AMREX_SPACEDIM == 2)
void MyTest::initializeLinearDataFor2D(int ilev) {
  GpuArray<Real, AMREX_SPACEDIM> dx = geom[ilev].CellSizeArray();
  for (MFIter mfi(phi[ilev]); mfi.isValid(); ++mfi) {
    const Box &bx = mfi.fabbox();
    Array4<Real> const &fab = phi[ilev].array(mfi);
    Array4<Real> const &fab_gx = grad_x_analytic[ilev].array(mfi);
    Array4<Real> const &fab_gy = grad_y_analytic[ilev].array(mfi);

    const FabArray<EBCellFlagFab> *flags =
        &(factory[ilev]->getMultiEBCellFlagFab());
    Array4<EBCellFlag const> const &flag = flags->const_array(mfi);

    Array4<Real const> const &ccent = (factory[ilev]->getCentroid()).array(mfi);
    Array4<Real const> const &fcx =
        (factory[ilev]->getFaceCent())[0]->const_array(mfi);
    Array4<Real const> const &fcy =
        (factory[ilev]->getFaceCent())[1]->const_array(mfi);

    GpuArray<int, AMREX_SPACEDIM> dlo;
    GpuArray<int, AMREX_SPACEDIM> dhi;
    for (int n=0; n < AMREX_SPACEDIM; ++n){
      dlo[n] = geom[ilev].Domain().loVect()[n];
      dhi[n] = geom[ilev].Domain().hiVect()[n];
    }

    Real H = linear_1d_height;
    Real t = (linear_1d_rotation / 180.) * M_PI;

    Real a = std::tan(t);
    Real b = -1.0;
    Real c = linear_1d_pt_on_top_wall[1] -
               std::tan(t) * linear_1d_pt_on_top_wall[0];

    GpuArray<const int, 3> is_periodic_tmp = {is_periodic[0], is_periodic[1]};

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

      Real rx = (i + 0.5 + ccent(i, j, k, 0)) * dx[0];
      Real ry = (j + 0.5 + ccent(i, j, k, 1)) * dx[1];

      // if not periodic, set the ghost cell values to corr. domain face values
      if (i < dlo[0] and not is_periodic_tmp[0]) {
        rx = dlo[0] * dx[0];
        ry = (j + 0.5 + fcx(i, j, k, 0)) * dx[1];
      }
      if (i > dhi[0] and not is_periodic_tmp[0]) {
        rx = (dhi[0] + 1) * dx[0];
        ry = (j + 0.5 + fcx(i, j, k, 0)) * dx[1];
      }
      if (j < dlo[1] and not is_periodic_tmp[1]) {
        rx = (i + 0.5 + fcy(i, j, k, 0)) * dx[0];
        ry = dlo[1] * dx[1];
      }
      if (j > dhi[1] and not is_periodic_tmp[1]) {
        rx = (i + 0.5 + fcy(i, j, k, 0)) * dx[0];
        ry = (dhi[1] + 1) * dx[1];
      }

      auto d = std::fabs(a * rx + b * ry + c) / std::sqrt(a * a + b * b);
      auto phi_mag = (!flag(i, j, k).isCovered()) ? (H - d) : 0.0;
      fab(i, j, k, 0) = phi_mag * std::cos(t);
      fab(i, j, k, 1) = phi_mag * std::sin(t);

      if (flag(i, j, k).isCovered()) {
        fab_gx(i, j, k, 0) = 0.0;
        fab_gx(i, j, k, 1) = 0.0;
        fab_gy(i, j, k, 0) = 0.0;
        fab_gy(i, j, k, 1) = 0.0;
      } else {

        Real fac =
          -1.0 ;
        fab_gx(i, j, k, 0) =
                (a * std::cos(t) / std::sqrt(a * a + b * b)) * fac * dx[0];
        fab_gx(i, j, k, 1) =
                 (a * std::sin(t) / std::sqrt(a * a + b * b)) * fac * dx[0];

        fac = -1.0 ;
        fab_gy(i, j, k, 0) =
                (b * std::cos(t) / std::sqrt(a * a + b * b)) * fac * dx[1];
        fab_gy(i, j, k, 1) =
                 (b * std::sin(t) / std::sqrt(a * a + b * b)) * fac * dx[1];
      }

    });
  }
}

#else
void MyTest::initializeLinearDataFor2D(int ilev) {
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(0, "Calling 2D function for 3D case. Error!!");
}
#endif

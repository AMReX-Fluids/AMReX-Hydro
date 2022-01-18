#include <AMReX_Array.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

#include "bds.H"


const bool limit_slopes = true;
constexpr amrex::Real eps = 1.0e-10;


using namespace amrex;


void bdsslope ( MultiFab const& s_mf, const Geometry& geom, MultiFab& slope_mf, int comp){

    BoxArray ba = s_mf.boxArray();
    DistributionMapping dmap = s_mf.DistributionMap();
    GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    // local variables
    MultiFab sint_mf(convert(ba,IntVect(AMREX_D_DECL(1,1,1))), dmap, 1, 1);

#if (AMREX_SPACEDIM ==2)
    Real hx = dx[0];
    Real hy = dx[1];


    for ( MFIter mfi(sint_mf); mfi.isValid(); ++mfi){

        const Box& bx = mfi.growntilebox(1);
        Array4<const Real> const& s    = s_mf.array(mfi, comp);
        Array4<      Real> const& sint = sint_mf.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k){

            // bicubic interpolation to corner points
            // (i,j,k) refers to lower corner of cell
            // Added k index -- placeholder for 2d
            sint(i,j,k) = (s(i-2,j-2,k) + s(i-2,j+1,k) + s(i+1,j-2,k) + s(i+1,j+1,k)
                    - 7.0*(s(i-2,j-1,k) + s(i-2,j  ,k) + s(i-1,j-2,k) + s(i  ,j-2,k) +
                           s(i-1,j+1,k) + s(i  ,j+1,k) + s(i+1,j-1,k) + s(i+1,j  ,k))
                   + 49.0*(s(i-1,j-1,k) + s(i  ,j-1,k) + s(i-1,j  ,k) + s(i  ,j  ,k)) ) / 144.0;

        });

    }

    for ( MFIter mfi(s_mf); mfi.isValid(); ++mfi){

        const Box& bx = mfi.growntilebox(1);

        Array4<const Real> const& s     = s_mf.array(mfi, comp);
        Array4<const Real> const& sint  = sint_mf.array(mfi);
        Array4<      Real> const& slope = slope_mf.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k){

            // compute initial estimates of slopes from unlimited corner points

            // local variables
            Real sumloc, redfac, redmax, div, kdp, sumdif, sgndif;

            Array1D<Real, 1, 4> diff;
            Array1D<Real, 1, 4> smin;
            Array1D<Real, 1, 4> smax;
            Array1D<Real, 1, 4> sc;

            // sx
            slope(i,j,k,0) = 0.5*(sint(i+1,j+1,k) + sint(i+1,j,k) - sint(i,j+1,k) - sint(i,j,k)) / hx;
            // sy
            slope(i,j,k,1) = 0.5*(sint(i+1,j+1,k) - sint(i+1,j,k) + sint(i,j+1,k) - sint(i,j,k)) / hy;
            // sxy
            slope(i,j,k,2) =     (sint(i+1,j+1,k) - sint(i+1,j,k) - sint(i,j+1,k) + sint(i,j,k)) / (hx*hy);

            if (limit_slopes) {

                // ++ / sint(i+1,j+1)
                sc(4) = s(i,j,k) + 0.5*(hx*slope(i,j,k,0) + hy*slope(i,j,k,1)) + 0.25*hx*hy*slope(i,j,k,2);

                // +- / sint(i+1,j  )
                sc(3) = s(i,j,k) + 0.5*(hx*slope(i,j,k,0) - hy*slope(i,j,k,1)) - 0.25*hx*hy*slope(i,j,k,2);

                // -+ / sint(i  ,j+1)
                sc(2) = s(i,j,k) - 0.5*(hx*slope(i,j,k,0) - hy*slope(i,j,k,1)) - 0.25*hx*hy*slope(i,j,k,2);

                // -- / sint(i  ,j  )
                sc(1) = s(i,j,k) - 0.5*(hx*slope(i,j,k,0) + hy*slope(i,j,k,1)) + 0.25*hx*hy*slope(i,j,k,2);

                // enforce max/min bounds
                smin(4) = min(s(i,j,k), s(i+1,j,k), s(i,j+1,k), s(i+1,j+1,k));
                smax(4) = max(s(i,j,k), s(i+1,j,k), s(i,j+1,k), s(i+1,j+1,k));

                smin(3) = min(s(i,j,k), s(i+1,j,k), s(i,j-1,k), s(i+1,j-1,k));
                smax(3) = max(s(i,j,k), s(i+1,j,k), s(i,j-1,k), s(i+1,j-1,k));

                smin(2) = min(s(i,j,k), s(i-1,j,k), s(i,j+1,k), s(i-1,j+1,k));
                smax(2) = max(s(i,j,k), s(i-1,j,k), s(i,j+1,k), s(i-1,j+1,k));

                smin(1) = min(s(i,j,k), s(i-1,j,k), s(i,j-1,k), s(i-1,j-1,k));
                smax(1) = max(s(i,j,k), s(i-1,j,k), s(i,j-1,k), s(i-1,j-1,k));

                for(int mm=1; mm<=4; ++mm){
                   sc(mm) = max(min(sc(mm), smax(mm)), smin(mm));
                }

                // iterative loop
                for(int ll=1; ll<=3; ++ll){
                   sumloc = 0.25*(sc(4) + sc(3) + sc(2) + sc(1));
                   sumdif = (sumloc - s(i,j,k))*4.0;
                   sgndif = std::copysign(1.0,sumdif);

                   for(int mm=1; mm<=4; ++mm){
                      diff(mm) = (sc(mm) - s(i,j,k))*sgndif;
                   }

                   kdp = 0;

                   for(int mm=1; mm<=4; ++mm){
                      if (diff(mm) > eps) {
                         kdp = kdp+1;
                      }
                   }

                   for(int mm=1; mm<=4; ++mm){
                      if (kdp<1) {
                         div = 1.0;
                      } else {
                         div = kdp;
                      }

                      if (diff(mm)>eps) {
                         redfac = sumdif*sgndif/div;
                         kdp = kdp-1;
                      } else {
                         redfac = 0.0;
                      }

                      if (sgndif > 0.0) {
                         redmax = sc(mm) - smin(mm);
                      } else {
                         redmax = smax(mm) - sc(mm);
                      }

                      redfac = min(redfac,redmax);
                      sumdif = sumdif - redfac*sgndif;
                      sc(mm) = sc(mm) - redfac*sgndif;
                   }
                }

                // final slopes
                // sx
                slope(i,j,k,0) = 0.5*( sc(4) + sc(3) -sc(1) - sc(2) )/hx;
                // sy
                slope(i,j,k,1) = 0.5*( sc(4) + sc(2) -sc(1) - sc(3) )/hy;
                // sxy
                slope(i,j,k,2) =     ( sc(1) + sc(4) -sc(2) - sc(3) )/(hx*hy);

            }


        });
    }
}

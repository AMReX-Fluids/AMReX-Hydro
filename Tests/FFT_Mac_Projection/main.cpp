#include <AMReX.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>

#include <hydro_FFTMacProjector.H>

using namespace amrex;
using namespace Hydro;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        int n_cell = 128;
        int max_grid_size = 32;
        {
            ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
        }

        Geometry geom;
        BoxArray grids;
        DistributionMapping dmap;
        {
            RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(1.,1.,1.)});

            Array<int,AMREX_SPACEDIM> isp{AMREX_D_DECL(0,1,0)};
            Box domain(IntVect(0), IntVect(n_cell-1));
            geom.define(domain, rb, 0, isp);

            grids.define(domain);
            grids.maxSize(max_grid_size);

            dmap.define(grids);
        }

        Array<MultiFab,AMREX_SPACEDIM> vel;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            vel[idim].define(amrex::convert(grids,IntVect::TheDimensionVector(idim)),
                             dmap, 1, 0);
            vel[idim].setVal(0);
            Box region = geom.Domain();
            region.growLo(idim, -geom.Domain().length(idim)/4);
            region.growHi(idim, -geom.Domain().length(idim)/4);
            vel[idim].setVal(Real(1.0), region, 0, 1);
        }

        MultiFab divu(grids, dmap, 1, 0);
        amrex::computeDivergence(divu, GetArrOfConstPtrs(vel), geom);
        amrex::Print() << "Divergence before projection: "
                       << divu.min(0) << " " << divu.max(0) << "\n";

        FFTMacProjector macproj(geom,
                                {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                              LinOpBCType::Periodic,
                                              LinOpBCType::Neumann)},
                                {AMREX_D_DECL(LinOpBCType::Neumann,
                                              LinOpBCType::Periodic,
                                              LinOpBCType::Dirichlet)});

        macproj.setUMAC(GetArrOfPtrs(vel));

        macproj.project();

        amrex::computeDivergence(divu, GetArrOfConstPtrs(vel), geom);
        amrex::Print() << "Divergence after projection: "
                       << divu.min(0) << " " << divu.max(0) << "\n";
    }
    amrex::Finalize();
}

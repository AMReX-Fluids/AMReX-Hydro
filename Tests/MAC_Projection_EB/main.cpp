#include <AMReX.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_TagBox.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>

#include <hydro_MacProjector.H>

using namespace amrex;
using namespace Hydro;

void write_plotfile(const Geometry& geom, const MultiFab& plotmf, int regtest)
{
    std::stringstream sstream;
    sstream << "plt00000";
    std::string plotfile_name = sstream.str();

    amrex::Print() << "Writing " << plotfile_name << '\n';

#if (AMREX_SPACEDIM == 2)
       EB_WriteSingleLevelPlotfile(plotfile_name, plotmf,
                                   { "proc" ,"xvel", "yvel" },
                                     geom, 0.0, 0);
#elif (AMREX_SPACEDIM == 3)

       if (regtest == 1)
       {
           MultiFab small_mf(plotmf.boxArray(), plotmf.DistributionMap(), AMREX_SPACEDIM, 0, MFInfo(), plotmf.Factory());
           MultiFab::Copy(small_mf, plotmf, 0, 0, small_mf.nComp(),0);

           EB_WriteSingleLevelPlotfile(plotfile_name, small_mf,
                                       { "proc", "xvel", "yvel" },
                                         geom, 0.0, 0);
       } else {
           EB_WriteSingleLevelPlotfile(plotfile_name, plotmf,
                                       { "proc", "xvel", "yvel", "zvel" },
                                         geom, 0.0, 0);
       }
#endif
}

int main (int argc, char* argv[])
{
    // Turn off amrex-related output
    amrex::SetVerbose(0);

    amrex::Initialize(argc, argv);

    auto strt_time = amrex::second();

    {
        int mg_verbose = 0;
        int bottom_verbose = 0;
        int n_cell = 128;
        int max_grid_size = 32;
        int use_hypre  = 0;
        int regtest   = 0;
        // With hypre: project again with mac_proj.use_mlhypre (the composite
        // HypreMLABecLap) and compare with the MLMG projection; then project
        // mlhypre_num_solves times in all with beta scaled by 2 each time
        // (the projected velocity does not depend on a uniform scaling of
        // beta), which exercises updateBeta and mac_proj.mlhypre_reuse_setup
        int mlhypre_check = 0;
        int mlhypre_num_solves = 1;
        int loc_centroid = 0;  // umac and beta at the face centroids instead of the face centers
        Real mlhypre_max_diff = 1.e-6;

        Real obstacle_radius = 0.10;
        Real reltol = 1.e-8; // Define the relative tolerance
        Real abstol = 1.e-15; // Define the absolute tolerance; note that this argument is optional

        // read parameters
        {
            ParmParse pp;
            pp.query("mg_verbose", mg_verbose);
            pp.query("bottom_verbose", bottom_verbose);
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("use_hypre", use_hypre);
            pp.query("regtest", regtest);
            pp.query("abstol", abstol);
            pp.query("reltol", reltol);
            pp.query("mlhypre_check", mlhypre_check);
            pp.query("mlhypre_num_solves", mlhypre_num_solves);
            pp.query("mlhypre_max_diff", mlhypre_max_diff);
            pp.query("loc_centroid", loc_centroid);
        }
        const auto face_loc = loc_centroid ? MLMG::Location::FaceCentroid : MLMG::Location::FaceCenter;

#ifndef AMREX_USE_HYPRE
        if (use_hypre == 1 || mlhypre_check == 1)
           amrex::Abort("Can't use hypre if we dont build with USE_HYPRE=TRUE");
#endif

#if (AMREX_SPACEDIM == 3)
        // Only n_cell_z, which is used in 3D only, needs n_cell to be divisible by 8
        if (n_cell%8 != 0)
           amrex::Abort("n_cell must be a multiple of 8");
#endif

        int n_cell_y =   n_cell;
        int n_cell_x = 2*n_cell;
        int n_cell_z =   n_cell/8;

        Real ylen = 1.0;
        Real xlen = 2.0 * ylen;
        Real zlen = ylen / 8.0;

        Geometry geom;
        BoxArray grids;
        DistributionMapping dmap;
        {
            RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(xlen,ylen,zlen)});

            Array<int,AMREX_SPACEDIM> isp{AMREX_D_DECL(0,1,1)};
            Geometry::Setup(&rb, 0, isp.data());
            Box domain(IntVect{AMREX_D_DECL(0,0,0)},
                       IntVect{AMREX_D_DECL(n_cell_x-1,n_cell_y-1,n_cell_z-1)});
            geom.define(domain);

            grids.define(domain);
            grids.maxSize(max_grid_size);

            dmap.define(grids);
        }

        Array<MultiFab,AMREX_SPACEDIM> vel;
        Array<MultiFab,AMREX_SPACEDIM> beta;
        MultiFab plotfile_mf;

        int required_coarsening_level = 0; // typically the same as the max AMR level index
        int max_coarsening_level = 100;    // typically a huge number so MG coarsens as much as possible

        amrex::Vector<amrex::RealArray> obstacle_center = {
            {AMREX_D_DECL(0.3,0.2,0.5)},
            {AMREX_D_DECL(0.3,0.5,0.5)},
            {AMREX_D_DECL(0.3,0.8,0.5)},
            {AMREX_D_DECL(0.7,0.25,0.5)},
            {AMREX_D_DECL(0.7,0.60,0.5)},
            {AMREX_D_DECL(0.7,0.85,0.5)},
            {AMREX_D_DECL(1.1,0.2,0.5)},
            {AMREX_D_DECL(1.1,0.5,0.5)},
            {AMREX_D_DECL(1.1,0.8,0.5)}};

        // The "false" below is the boolean that determines if the fluid is inside ("true") or
        //     outside ("false") the object(s)

#if (AMREX_SPACEDIM == 3)
        int direction =  2;
        Real height   = -1.0;  // Putting a negative number for height means it extends beyond the domain

        auto make_obstacle = [&] (Real radius, int n) {
            return EB2::CylinderIF(radius, height, direction, obstacle_center[n], false); };
#else
        // In 2D the z-aligned cylinders degenerate into circles in the x-y plane
        auto make_obstacle = [&] (Real radius, int n) {
            return EB2::SphereIF(radius, obstacle_center[n], false); };
#endif

        Array<decltype(make_obstacle(Real(0.0),0)),9> obstacles{
            make_obstacle(    obstacle_radius, 0),
            make_obstacle(    obstacle_radius, 1),
            make_obstacle(    obstacle_radius, 2),
            make_obstacle(0.9*obstacle_radius, 3),
            make_obstacle(0.9*obstacle_radius, 4),
            make_obstacle(0.9*obstacle_radius, 5),
            make_obstacle(    obstacle_radius, 6),
            make_obstacle(    obstacle_radius, 7),
            make_obstacle(    obstacle_radius, 8)};

        auto group_1 = EB2::makeUnion(obstacles[0],obstacles[1],obstacles[2]);
        auto group_2 = EB2::makeUnion(obstacles[3],obstacles[4],obstacles[5]);
        auto group_3 = EB2::makeUnion(obstacles[6],obstacles[7],obstacles[8]);
        auto all     = EB2::makeUnion(group_1,group_2,group_3);
        auto gshop9  = EB2::makeShop(all);
        EB2::Build(gshop9, geom, required_coarsening_level, max_coarsening_level);

        const EB2::IndexSpace& eb_is = EB2::IndexSpace::top();
        const EB2::Level& eb_level = eb_is.getLevel(geom);

        // options are basic, volume, or full
        EBSupport ebs = EBSupport::full;

        // number of ghost cells for each of the 3 EBSupport types
        Vector<int> ng_ebs = {2,2,2};

        // This object provides access to the EB database in the format of basic AMReX objects
        // such as BaseFab, FArrayBox, FabArray, and MultiFab
        EBFArrayBoxFactory factory(eb_level, geom, grids, dmap, ng_ebs, ebs);

        // allocate face-centered velocities and face-centered beta coefficient
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            vel[idim].define (amrex::convert(grids,IntVect::TheDimensionVector(idim)), dmap, 1, 1,
                              MFInfo(), factory);
            beta[idim].define(amrex::convert(grids,IntVect::TheDimensionVector(idim)), dmap, 1, 0,
                              MFInfo(), factory);
            beta[idim].setVal(1.0);  // set beta to 1
        }

        // If we want to supply a non-zero S we must allocate and fill it outside the solver
        // MultiFab S(grids, dmap, 1, 0, MFInfo(), factory);
        // Set S here ...

        // store plotfile variables; velocity and processor id
        plotfile_mf.define(grids, dmap, AMREX_SPACEDIM+1, 0, MFInfo(), factory);

        // set initial velocity to u=(1,0,0)
        AMREX_D_TERM(vel[0].setVal(1.0);,
                     vel[1].setVal(0.0);,
                     vel[2].setVal(0.0););

        LPInfo lp_info;

        // If we want to use hypre to solve the full problem we need to not coarsen inside AMReX
        if (use_hypre)
            lp_info.setMaxCoarseningLevel(0);


        // Note that when we build with USE_EB = TRUE, we must specify whether the quantities are located
        //  at centers (MLMG::Location::CellCenter, MLMG::Location::FaceCenter) or
        //  centroids (MLMG::Location::CellCentroid, MLMG::Location::FaceCentroid).
        MacProjector macproj({amrex::GetArrOfPtrs(vel)},       // mac velocity
                             face_loc,                         // Location of vel
                             {amrex::GetArrOfConstPtrs(beta)}, // beta
                             face_loc,                         // Location of beta
                             MLMG::Location::CellCenter,       // Location of solution variable phi
                             {geom},                           // the geometry object
                             lp_info);                         // structure for passing info to the operator

        // Here we specify the desired divergence S
        // MacProjector macproj({amrex::GetArrOfPtrs(vel)},       // mac velocity
        //                      MLMG::Location::FaceCenter,       // Location of vel
        //                      {amrex::GetArrOfConstPtrs(beta)}, // beta
        //                      MLMG::Location::FaceCenter,       // Location of beta
        //                      MLMG::Location::CellCenter,       // Location of solution variable phi
        //                      {geom},                           // the geometry object
        //                      lp_info,                          // structure for passing info to the operator
        //                      {&S},                             // defines the specified RHS divergence
        //                      MLMG::Location::CellCenter);      // Location of S


        // Set bottom-solver to use hypre instead of native BiCGStab
        if (use_hypre)
            macproj.getMLMG().setBottomSolver(MLMG::BottomSolver::hypre);

        // Hard-wire the boundary conditions to be Neumann on the low x-face, Dirichlet
        // on the high x-face, and periodic in the other two directions
        // (the first argument is for the low end, the second is for the high end)
        macproj.setDomainBC({AMREX_D_DECL(LinOpBCType::Neumann,
                                          LinOpBCType::Periodic,
                                          LinOpBCType::Periodic)},
                            {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                          LinOpBCType::Periodic,
                                          LinOpBCType::Periodic)});

        macproj.setVerbose(mg_verbose);
        macproj.getMLMG().setBottomVerbose(bottom_verbose);


        amrex::Print() << " \n********************************************************************" << '\n';
        amrex::Print() << " Let's project the initial velocity to find " << '\n';
        amrex::Print() << "   the flow field around the obstacles ... " << '\n';
        amrex::Print() << " The domain has " << n_cell_x << " cells in the x-direction "          << '\n';
        amrex::Print() << " The maximum grid size is " << max_grid_size                             << '\n';
        amrex::Print() << "******************************************************************** \n" << '\n';

        // Solve for phi and subtract from the velocity to make it divergence-free
        // Note that the normal velocities are at face centers (not centroids) unless loc_centroid
        macproj.project(reltol,abstol);

        // If we want to use phi elsewhere, we can pass in an array in which to return the solution.
    // This can also be used to pass in an initial guess for phi if one is available; here we
    // simply set it to zero.
    //
        // MultiFab phi_inout(grids, dmap, 1, 1, MFInfo(), factory);
    // phi_inout.setVal(0.);
        // macproj.project_center_vels({&phi_inout},reltol,abstol,MLMG::Location::FaceCenter);

        amrex::Print() << " \n********************************************************************" << '\n';
        amrex::Print() << " Done!" << '\n';
        amrex::Print() << "******************************************************************** \n" << '\n';

        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            vel[idim].FillBoundary(geom.periodicity());
        }

#ifdef AMREX_USE_HYPRE
        if (mlhypre_check) {
            // The same projection with the composite solver
            ParmParse ppm("mac_proj");
            ppm.add("use_mlhypre", 1);

            Array<MultiFab,AMREX_SPACEDIM> vel2;
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                vel2[idim].define(vel[idim].boxArray(), vel[idim].DistributionMap(), 1, 1,
                                  MFInfo(), factory);
                vel2[idim].setVal(idim == 0 ? Real(1.0) : Real(0.0));
            }

            MacProjector macproj2({amrex::GetArrOfPtrs(vel2)}, face_loc,
                                  {amrex::GetArrOfConstPtrs(beta)}, face_loc,
                                  MLMG::Location::CellCenter, {geom}, lp_info);
            macproj2.setDomainBC({AMREX_D_DECL(LinOpBCType::Neumann,
                                               LinOpBCType::Periodic,
                                               LinOpBCType::Periodic)},
                                 {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                               LinOpBCType::Periodic,
                                               LinOpBCType::Periodic)});
            macproj2.setVerbose(mg_verbose);

            auto compare = [&] (std::string const& what) {
                Real max_diff = Real(0.0);
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                    MultiFab diff(vel[idim].boxArray(), vel[idim].DistributionMap(), 1, 0,
                                  MFInfo(), factory);
                    MultiFab::Copy(diff, vel2[idim], 0, 0, 1, 0);
                    MultiFab::Subtract(diff, vel[idim], 0, 0, 1, 0);
                    max_diff = std::max(max_diff, diff.norm0(0, 0, false, true));
                }
                amrex::Print() << what << ": max |u_mlhypre - u_mlmg| = " << max_diff
                               << " (" << macproj2.getNumIters() << " iterations)\n";
                if (max_diff > mlhypre_max_diff) {
                    amrex::Abort(what + ": the composite and MLMG projections differ");
                }
            };

            macproj2.project(reltol, abstol);
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                vel2[idim].FillBoundary(geom.periodicity());
            }
            compare("mlhypre");

            for (int k = 1; k < mlhypre_num_solves; ++k) {
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                    beta[idim].mult(Real(2.0));
                    vel2[idim].setVal(idim == 0 ? Real(1.0) : Real(0.0));
                }
                macproj2.updateBeta({amrex::GetArrOfConstPtrs(beta)});
                macproj2.project(reltol, abstol);
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                    vel2[idim].FillBoundary(geom.periodicity());
                }
                compare("mlhypre solve " + std::to_string(k+1));
            }
        }
#endif

        // Copy processor id into plotfile_mf
        plotfile_mf.setVal(ParallelDescriptor::MyProc(), 0, 1);

        // Copy velocity into plotfile (starting at component 1)
        average_face_to_cellcenter(plotfile_mf,1,amrex::GetArrOfConstPtrs(vel));

        write_plotfile(geom, plotfile_mf, regtest);
    }

    auto stop_time = amrex::second() - strt_time;
    amrex::Print() << "Total run time " << stop_time << '\n';

    amrex::Finalize();
}

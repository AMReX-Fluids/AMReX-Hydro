#include <AMReX.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_MLMG.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>

#include <hydro_NodalProjector.H>

using namespace amrex;

void write_plotfile(const Geometry& geom, const MultiFab& plotmf)
{
    std::string plotfile_name("plt00000");

    amrex::Print() << "Writing " << plotfile_name << '\n';

#if (AMREX_SPACEDIM == 2)
       EB_WriteSingleLevelPlotfile(plotfile_name, plotmf,
                                   { "proc" ,"xvel", "yvel" },
                                     geom, 0.0, 0);
#elif (AMREX_SPACEDIM == 3)
       EB_WriteSingleLevelPlotfile(plotfile_name, plotmf,
                                   { "proc", "xvel", "yvel", "zvel" },
                                     geom, 0.0, 0);
#endif
}

int main (int argc, char* argv[])
{
    // Turn off amrex-related output
    amrex::SetVerbose(0);

    amrex::Initialize(argc, argv);

    auto strt_time = amrex::second();

    {
        int n_cell = 128;
        int max_grid_size = 32;
        int use_hypre  = 0;
        // With hypre: project again with nodal_proj.use_mlhypre (the composite
        // HypreMLNodeLap) and compare with the MLMG projection; then project
        // mlhypre_num_solves times in all with the same projector, sigma
        // scaled by 2 each time (the projected velocity does not depend on a
        // uniform scaling of sigma), which exercises the matrix update of the
        // kept composite solver and nodal_proj.mlhypre_reuse_setup
        int mlhypre_check = 0;
        int mlhypre_num_solves = 1;
        Real mlhypre_max_diff = 1.e-6;

        Real obstacle_radius = 0.10;
        Real reltol = 1.e-8;    // Define the relative tolerance
        Real abstol = 1.e-15;   // Define the absolute tolerance; note that this argument is optional


        // read parameters
        {
            ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("use_hypre", use_hypre);
            pp.query("mlhypre_check", mlhypre_check);
            pp.query("mlhypre_num_solves", mlhypre_num_solves);
            pp.query("mlhypre_max_diff", mlhypre_max_diff);
        }


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
            Box domain(IntVect{AMREX_D_DECL(0,0,0)},
                       IntVect{AMREX_D_DECL(n_cell_x-1,n_cell_y-1,n_cell_z-1)});
            geom.define(domain, rb, CoordSys::cartesian, isp);

            grids.define(domain);
            grids.maxSize(max_grid_size);

            dmap.define(grids);
        }

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

        //
        // Given a cell-centered velocity (vel) field, a cell-centered
        // scalar field (sigma) field, and a source term S (either node-
        // or cell-centered) solve:
        //
        //   div( sigma * grad(phi) ) = div(vel) - S
        //
        // and then perform the projection:
        //
        //     vel = vel - sigma * grad(phi)
        //

        //
        //  Create the cell-centered velocity field we want to project.
        //  Set velocity field to (1,0,0) including ghost cells for this example.
        //
        MultiFab vel(grids, dmap, AMREX_SPACEDIM, 1, MFInfo(), factory);
        vel.setVal(1.0, 0, 1, 1);
        vel.setVal(0.0, 1, AMREX_SPACEDIM-1, 1);

        //
        // Create the cell-centered sigma field and set it to 1 for this example
        //
        MultiFab sigma(grids, dmap, 1, 1, MFInfo(), factory);
        sigma.setVal(1.0);

        //
        // Create cell-centered contributions to RHS and set it to zero for this example
        //
        MultiFab S_cc(grids, dmap, 1, 1, MFInfo(), factory);
        S_cc.setVal(0.0);

        //
        // Create node-centered contributions to RHS and set it to zero for this example
        //
        const BoxArray & nd_grids = amrex::convert(grids, IntVect::TheNodeVector()); // nodal grids
        MultiFab S_nd(nd_grids, dmap, 1, 1, MFInfo(), factory);
        S_nd.setVal(0.0);

        //
        // Setup linear operator, AKA the nodal Laplacian
        //
        LPInfo lp_info;

        // If we want to use hypre to solve the full problem we do not need to coarsen the GMG stencils
        // if (use_hypre_as_full_solver)
        //    lp_info.setMaxCoarseningLevel(0);

        // Setup nodal projector object
        Hydro::NodalProjector nodal_proj({&vel}, {&sigma}, {geom}, lp_info, {&S_cc}, {&S_nd});

        // Set boundary conditions.
        // Here we use Neumann on the low x-face, Dirichlet on the high x-face,
        // and periodic in the other two directions
        // (the first argument is for the low end, the second is for the high end)
        // Note that Dirichlet boundary conditions are assumed to be homogeneous (i.e. phi = 0)
        nodal_proj.setDomainBC({AMREX_D_DECL(LinOpBCType::Neumann,
                                         LinOpBCType::Periodic,
                                         LinOpBCType::Periodic)},
                           {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                          LinOpBCType::Periodic,
                                         LinOpBCType::Periodic)});


        amrex::Print() << " \n********************************************************************" << '\n';
        amrex::Print() << " Let's project the initial velocity to find " << '\n';
        amrex::Print() << "   the flow field around the obstacles ... " << '\n';
        amrex::Print() << " The domain has " << n_cell_x << " cells in the x-direction "          << '\n';
        amrex::Print() << " The maximum grid size is " << max_grid_size                             << '\n';
        amrex::Print() << "******************************************************************** \n" << '\n';

        //
        // Solve div( sigma * grad(phi) ) = RHS
        //
        nodal_proj.project( reltol, abstol);

        // Optionally, the projection can return the resulting phi and/or phi can be used to provide
        // an initial guess if available.
        //
        // MultiFab phi(nd_grids, dmap, 1, 1, MFInfo(), factory);
        // phi.setVal(0.0); // Must initialize phi; we simply set to 0 for this example.
        // nodal_proj.project( {&phi}, reltol, abstol);

        amrex::Print() << " \n********************************************************************" << '\n';
        amrex::Print() << " Projection complete " << '\n';
        amrex::Print() << "******************************************************************** \n" << '\n';

#ifdef AMREX_USE_HYPRE
        if (mlhypre_check) {
            // The same projection with the composite solver
            ParmParse ppn("nodal_proj");
            ppn.add("use_mlhypre", 1);

            MultiFab vel2(grids, dmap, AMREX_SPACEDIM, 1, MFInfo(), factory);
            vel2.setVal(1.0, 0, 1, 1);
            vel2.setVal(0.0, 1, AMREX_SPACEDIM-1, 1);
            MultiFab sigma2(grids, dmap, 1, 1, MFInfo(), factory);
            MultiFab::Copy(sigma2, sigma, 0, 0, 1, 1);

            Hydro::NodalProjector nodal_proj2({&vel2}, {&sigma2}, {geom}, lp_info, {&S_cc}, {&S_nd});
            nodal_proj2.setDomainBC({AMREX_D_DECL(LinOpBCType::Neumann,
                                                  LinOpBCType::Periodic,
                                                  LinOpBCType::Periodic)},
                                    {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                                  LinOpBCType::Periodic,
                                                  LinOpBCType::Periodic)});

            auto compare = [&] (MultiFab const& u, std::string const& what) {
                MultiFab diff(grids, dmap, AMREX_SPACEDIM, 0, MFInfo(), factory);
                MultiFab::Copy(diff, u, 0, 0, AMREX_SPACEDIM, 0);
                MultiFab::Subtract(diff, vel, 0, 0, AMREX_SPACEDIM, 0);
                Real max_diff = Real(0.0);
                for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                    max_diff = std::max(max_diff, diff.norm0(n, 0, false, true));
                }
                amrex::Print() << what << ": max |u_mlhypre - u_mlmg| = " << max_diff
                               << " (" << nodal_proj2.getNumIters() << " iterations)\n";
                if (max_diff > mlhypre_max_diff) {
                    amrex::Abort(what + ": the composite and MLMG projections differ");
                }
            };

            nodal_proj2.project(reltol, abstol);
            compare(vel2, "mlhypre");

            // The kept projector, pointed at new data on the same grids
            MultiFab vel3(grids, dmap, AMREX_SPACEDIM, 1, MFInfo(), factory);
            for (int k = 1; k < mlhypre_num_solves; ++k) {
                MultiFab& u = (k % 2 == 1) ? vel3 : vel2;
                u.setVal(1.0, 0, 1, 1);
                u.setVal(0.0, 1, AMREX_SPACEDIM-1, 1);
                sigma2.mult(Real(2.0), 0, 1, 1);
                nodal_proj2.setVelocity({&u});
                nodal_proj2.setSigma({&sigma2});
                nodal_proj2.project(reltol, abstol);
                compare(u, "mlhypre solve " + std::to_string(k+1));
            }
        }
#endif

        // Store plotfile variables; velocity and processor id
        MultiFab plotfile_mf(grids, dmap, AMREX_SPACEDIM+1, 0, MFInfo(), factory);

        // copy processor id into plotfile_mf
        plotfile_mf.setVal(ParallelDescriptor::MyProc(), 0, 1);
        plotfile_mf.setVal(ParallelDescriptor::MyProc(), 0, 1);

        // copy velocity into plotfile
        MultiFab::Copy(plotfile_mf, vel, 0, 1, AMREX_SPACEDIM, 0);

        write_plotfile(geom, plotfile_mf);
    }

    auto stop_time = amrex::second() - strt_time;
    amrex::Print() << "Total run time " << stop_time << '\n';

    amrex::Finalize();
}

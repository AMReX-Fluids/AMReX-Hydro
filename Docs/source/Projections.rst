.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran


Projections
===========

MAC Projection
--------------

Some codes define a velocity field :math:`U = (u,v,w)` on faces, i.e.
:math:`u` is defined on x-faces, :math:`v` is defined on y-faces,
and :math:`w` is defined on z-faces.   We refer to the exact projection
of this velocity field as a MAC projection, in which we solve

.. math::

   D( \beta \nabla \phi) = D(U^*) - S

for :math:`\phi` and then set

.. math::

   U = U^* - \beta \nabla \phi


where :math:`U^*` is a vector field (typically velocity) that we want to satisfy
:math:`D(U) = S`.  For incompressible flow,  :math:`S = 0`.

The MacProjection class can be defined and used to perform the MAC projection without explicitly
calling the solver directly.  In addition to solving the variable coefficient Poisson equation,
the MacProjector internally computes the divergence of the vector field, :math:`D(U^*)`,
to compute the right-hand-side, and after the solve, subtracts the weighted gradient term to
make the vector field result satisfy the divergence constraint.

In the simplest form of the call, :math:`S` is assumed to be zero and does not need to be specified.
Typically, the user does not allocate the solution array, but it is also possible to create and pass
in the solution array and have :math:`\phi` returned as well as :math:`U`.

Caveat:  Currently the MAC projection only works when the base level covers the full domain; it does
not yet have the interface to pass boundary conditions for a fine level that come from coarser data.

Also note that any Dirichlet or Neumann boundary conditions at domain boundaries
are assumed to be homogeneous.  The call to the :cpp:`MLLinOp` member function
:cpp:`setLevelBC` occurs inside the MacProjection class; one does not need to call that
explicitly when using the MacProjection class.

The code below is taken from the file
``Tests/MAC_Projection_EB/main.cpp`` in AMReX-Hydro, and demonstrates how to set up
the MACProjector object and use it to perform a MAC projection.

.. highlight:: c++

::

    EBFArrayBoxFactory factory(eb_level, geom, grids, dmap, ng_ebs, ebs);

    // allocate face-centered velocities and face-centered beta coefficient
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        vel[idim].define (amrex::convert(grids,IntVect::TheDimensionVector(idim)), dmap, 1, 1,
                          MFInfo(), factory);
        beta[idim].define(amrex::convert(grids,IntVect::TheDimensionVector(idim)), dmap, 1, 0,
                          MFInfo(), factory);
        beta[idim].setVal(1.0);  // set beta to 1
    }

    // If we want to use phi elsewhere, we must create an array in which to return the solution
    // MultiFab phi_inout(grids, dmap, 1, 1, MFInfo(), factory);

    // If we want to supply a non-zero S we must allocate and fill it outside the solver
    // MultiFab S(grids, dmap, 1, 0, MFInfo(), factory);
    // Set S here ...

    // set initial velocity to U=(1,0,0)
    AMREX_D_TERM(vel[0].setVal(1.0);,
                 vel[1].setVal(0.0);,
                 vel[2].setVal(0.0););

    LPInfo lp_info;

    // If we want to use hypre to solve the full problem we do not need to coarsen the GMG stencils
    if (use_hypre_as_full_solver)
        lp_info.setMaxCoarseningLevel(0);

    MacProjector macproj({amrex::GetArrOfPtrs(vel)},       // face-based velocity
                         {amrex::GetArrOfConstPtrs(beta)}, // beta
                         {geom},                           // the geometry object
                         lp_info);                         // structure for passing info to the operator

    // Here we specify the desired divergence S
    // MacProjector macproj({amrex::GetArrOfPtrs(vel)},       // face-based velocity
    //                      {amrex::GetArrOfConstPtrs(beta)}, // beta
    //                      {geom},                           // the geometry object
    //                      lp_info,                          // structure for passing info to the operator
    //                      {&S});                            // defines the specified RHS divergence

    // Set bottom-solver to use hypre instead of native BiCGStab
    if (use_hypre_as_full_solver || use_hypre_as_bottom_solver)
       macproj.setBottomSolver(MLMG::BottomSolver::hypre);

    // Set boundary conditions.
    //  Here we use Neumann on the low x-face, Dirichlet on the high x-face,
    //  and periodic in the other two directions
    //  (the first argument is for the low end, the second is for the high end)
    // Note that Dirichlet and Neumann boundary conditions are assumed to be homogeneous.
    macproj.setDomainBC({AMREX_D_DECL(LinOpBCType::Neumann,
                                      LinOpBCType::Periodic,
                                      LinOpBCType::Periodic)},
                        {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                      LinOpBCType::Periodic,
                                      LinOpBCType::Periodic)});

    macproj.setVerbose(mg_verbose);
    macproj.setBottomVerbose(bottom_verbose);

    // Define the relative tolerance
    Real reltol = 1.e-8;

    // Define the absolute tolerance; note that this argument is optional
    Real abstol = 1.e-15;

    // Solve for phi and subtract from the velocity to make it divergence-free
    // Note that when we build with USE_EB = TRUE, we must specify whether the velocities live
    //  at face centers (MLMG::Location::FaceCenter) or face centroids (MLMG::Location::FaceCentroid)
    macproj.project(reltol,abstol,MLMG::Location::FaceCenter);

    // If we want to use phi elsewhere, we can pass in an array in which to return the solution
    // macproj.project({&phi_inout},reltol,abstol,MLMG::Location::FaceCenter);

Nodal Projection
----------------

Some codes define a velocity field :math:`U = (u,v,w)` with all
components co-located on cell centers.  The nodal solver in AMReX
can be used to compute an approximate projection of the cell-centered
velocity field, with pressure and velocity divergence defined on nodes.
When we use the nodal solver this way, and subtract only the cell average
of the gradient from the velocity, it is effectively an approximate projection.

As with the MAC projection, consider that we want to solve

.. math::

   D( \beta \nabla \phi) = D(U^*) - S

for :math:`\phi` and then set

.. math::

   U = U^* - \beta \nabla \phi

where :math:`U^*` is a vector field defined on cell centers and we want to satisfy
:math:`D(U) = S`.  For incompressible flow,  :math:`S = 0`.

Currently this nodal approximate projection does not exist in a separate
operator like the MAC projection; instead we demonstrate below the steps needed
to compute the approximate projection.  This means we must compute explicitly the
right-hand-side , including the the divergence of the vector field, :math:`D(U^*)`,
solve the variable coefficient Poisson equation, then subtract the weighted
gradient term to make the vector field result satisfy the divergence constraint.

The code below is taken from the file
``Tests/Nodal_Projection_EB/main.cpp`` in AMReX-Hydro, and demonstrates how to set up
the NodalProjector object and use it to perform a nodal projection.

.. highlight:: c++

::

   //
   // Given a cell-centered velocity (vel) field, a cell-centered
   // scalar field (sigma) field, and a source term S (either node-
   // or cell-centered )solve:
   //
   //   div( sigma * grad(phi) ) = div(vel) - S
   //
   // and then perform the projection:
   //
   //     vel = vel - sigma * grad(phi)
   //

   //
   // Create the EB factory
   //
   EBFArrayBoxFactory factory(eb_level, geom, grids, dmap, ng_ebs, ebs);

   //
   //  Create the cell-centered velocity field we want to project
   //
   MultiFab vel(grids, dmap, AMREX_SPACEDIM, 1, MFInfo(), factory);

   // Set velocity field to (1,0,0) including ghost cells for this example
   vel.setVal(1.0, 0, 1, 1);
   vel.setVal(0.0, 1, AMREX_SPACEDIM-1, 1);

   //
   // Setup linear operator, AKA the nodal Laplacian
   //
   LPInfo lp_info;

   // If we want to use hypre to solve the full problem we do not need to coarsen the GMG stencils
   // if (use_hypre_as_full_solver)
   //     lp_info.setMaxCoarseningLevel(0);

   MLNodeLaplacian matrix({geom}, {grids}, {dmap}, lp_info,
                          Vector<EBFArrayBoxFactory const*>{&factory});

   // Set boundary conditions.
   // Here we use Neumann on the low x-face, Dirichlet on the high x-face,
   // and periodic in the other two directions
   // (the first argument is for the low end, the second is for the high end)
   // Note that Dirichlet boundary conditions are assumed to be homogeneous (i.e. phi = 0)
   matrix.setDomainBC({AMREX_D_DECL(LinOpBCType::Neumann,
                                    LinOpBCType::Periodic,
                                    LinOpBCType::Periodic)},
                      {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                    LinOpBCType::Periodic,
                                    LinOpBCType::Periodic)});

   // Set matrix attributes to be used by MLMG solver
   matrix.setGaussSeidel(true);
   matrix.setHarmonicAverage(false);

   //
   // Compute RHS
   //
   // NOTE: it's up to the user to compute the RHS. as opposed
   //       to the MAC projection case !!!
   //
   // NOTE: do this operation AFTER setting up the linear operator so
   //       that compRHS method can be used
   //

   // RHS is nodal
   const BoxArray & nd_grids = amrex::convert(grids, IntVect{1,1,1}); // nodal grids

   // MultiFab to host RHS
   MultiFab rhs(nd_grids, dmap, 1, 1, MFInfo(), factory);

   // Cell-centered contributions to RHS
   MultiFab S_cc(grids, dmap, 1, 1, MFInfo(), factory);
   S_cc.setVal(0.0); // Set it to zero for this example

   // Node-centered contributions to RHS
   MultiFab S_nd(nd_grids, dmap, 1, 1, MFInfo(), factory);
   S_nd.setVal(0.0); // Set it to zero for this example

   // Compute RHS -- vel must be cell-centered
   matrix.compRHS({&rhs}, {&vel}, {&S_nd}, {&S_cc});

   //
   // Create the cell-centered sigma field and set it to 1 for this example
   //
   MultiFab sigma(grids, dmap, 1, 1, MFInfo(), factory);
   sigma.setVal(1.0);

   // Set sigma
   matrix.setSigma(0, sigma);

   //
   // Create node-centered phi
   //
   MultiFab phi(nd_grids, dmap, 1, 1, MFInfo(), factory);
   phi.setVal(0.0);

   //
   // Setup MLMG solver
   //
   MLMG nodal_solver(matrix);

   // We can specify the maximum number of iterations
   nodal_solver.setMaxIter(mg_maxiter);
   nodal_solver.setBottomMaxIter(mg_bottom_maxiter);

   nodal_solver.setVerbose(mg_verbose);
   nodal_solver.setBottomVerbose(mg_bottom_verbose);

   // Set bottom-solver to use hypre instead of native BiCGStab
   //   ( we could also have set this to cg, bicgcg, cgbicg)
   // if (use_hypre_as_full_solver || use_hypre_as_bottom_solver)
   //     nodal_solver.setBottomSolver(MLMG::BottomSolver::hypre);

   // Define the relative tolerance
   Real reltol = 1.e-8;

   // Define the absolute tolerance; note that this argument is optional
   Real abstol = 1.e-15;

   //
   // Solve div( sigma * grad(phi) ) = RHS
   //
   nodal_solver.solve( {&phi}, {&rhs}, reltol, abstol);

   //
   // Create cell-centered MultiFab to hold value of -sigma*grad(phi) at cell-centers
   //
   //
   MultiFab fluxes(grids, dmap, AMREX_SPACEDIM, 1, MFInfo(), factory);
   fluxes.setVal(0.0);

   // Get fluxes from solver
   nodal_solver.getFluxes( {&fluxes} );

   //
   // Apply projection explicitly --  vel = vel - sigma * grad(phi)
   //
   MultiFab::Add( *vel, *fluxes, 0, 0, AMREX_SPACEDIM, 0);


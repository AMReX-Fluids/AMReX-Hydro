.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran

.. _projections:

Projection Methods
==================

Here, we first include a brief discussion of projection methodology for incompressible and low Mach number flow.
Details of AMReX-Hydro's implementation of projections follow.

The compressible Navier-Stokes equations can be written in the form:

.. math:: {{\bf U}}_t + \nabla \cdot F({{\bf U}}) = S

where :math:`{{\bf U}}` is a vector of conserved quantities, :math:`{{\bf U}}= (\rho, \rho u,
\rho E)`, with :math:`\rho` the density, :math:`u` the velocity, :math:`E` the total
energy per unit mass, and :math:`S` are source terms. This system
can be expressed as a coupled set of advection/diffusion equations:

.. math:: {\bf q}_t + A({\bf q}) \nabla {\bf q} + D = {\cal S}

where :math:`{\bf q}` are called the primitive variables, :math:`A` is the advective
flux Jacobian, :math:`A \equiv \partial F / \partial U`, :math:`D` are diffusion terms,
and :math:`{\cal S}` are the transformed sources. The eigenvalues of the
matrix :math:`A` are the characteristic speeds—the real-valued speeds at which
information propagates in the system, :math:`u` and :math:`u
\pm c`, where :math:`c` is the sound speed. Solution methods for the
compressible equations that are strictly conservative make use of this wave-nature to compute advective fluxes
at the interfaces of grid cells. Diffusive fluxes can be computed
either implicit or explicit in time, and are added to the advective fluxes,
and used, along with the source terms to update the state in time. An
excellent introduction to these methods is provided by LeVeque’s book
:cite:`leveque`. The timestep for these methods is limited by all three processes
and their numerical implementation. Typically, advection terms are treated
time-explicitly, and the time step will be constrained by the time
it takes for the maximum characteristic speed to traverse one grid cell.
However, in low speed flow applications, it can be shown the acoustics
transport very little energy in the system. As a result, the time-step
restrictions arising from numerical treatement of the advection terms
can be unnecessarily limited, even if A-stable methods are used to incorporate
the diffusion and source terms.

In contrast, solving incompressible or low Mach number systems
typically involves a stage where one or more
advection-like equations are solved (representing, e.g. conservation of mass and
momentum), and coupling that advance with a divergence constraint on the velocity field.
For example, the equations of invicid constant-density incompressible flow
are:

.. math::

   \begin{aligned}
   {{\bf U}}_t = & -{{\bf U}}\cdot \nabla {{\bf U}}- \frac{1}{\rho}\nabla p \label{eq:incompressible_u} \\
   \nabla \cdot {{\bf U}} = &\  0 \end{aligned}

Here, :math:`{{\bf U}}` represents the velocity vector
and :math:`p` is the dynamical pressure. The time-evolution equation for
the velocity (Eq. [eq:incompressible\_u]) can be solved using
techniques similar to those developed for compressible hydrodynamics,
updating the old velocity, :math:`{{\bf U}}^n`, to the new time-level, :math:`{{\bf U}}^\star`.
Here the “:math:`^\star`” indicates that the updated velocity does not, in
general, satisfy the divergence constraint. A projection method will
take this updated velocity and force it to obey the constraint
equation. The basic idea follows from the fact that any vector
field can be expressed as the sum of a divergence-free quantity and
the gradient of a scalar. For the velocity, we can write:

.. math:: {{\bf U}}^\star = {{\bf U}}^d + \nabla \phi \label{eq:decomposition}

where :math:`{{\bf U}}^d` is the divergence free portion of the velocity vector,
:math:`{{\bf U}}^\star`, and :math:`\phi` is a scalar. Taking the divergence of
Eq. ([eq:decomposition]), we have

.. math:: \nabla^2 \phi = \nabla \cdot {{\bf U}}^\star

(where we used :math:`\nabla \cdot {{\bf U}}^d = 0`).
With appropriate boundary conditions, this Poisson equation can be
solved for :math:`\phi`, and the final, divergence-free velocity can
be computed as

.. math:: {{\bf U}}^{n+1} = {{\bf U}}^\star - \nabla \phi

Because soundwaves are filtered, the timestep constraint now depends only
on :math:`|{{\bf U}}|`.

Extensions to variable-density incompressible
flows :cite:`bellMarcus:1992b` involve a slightly different
decomposition of the velocity field and, as a result, a slightly
different, variable-coefficient Poisson equation.
There are also a variety of different ways
to express what is being projected :cite:`almgren:bell:crutchfield`,
and different discretizations of the divergence and gradient operators
lead to slightly different mathematical properties of the methods
(leading to “approximate
projections” :cite:`almgrenBellSzymczak:1996`).

For second-order methods, two projections are typically done per timestep.
First, the ‘MAC’ projection :cite:`bellColellaHowell:1991`
operates on the half-time, edge-centered advective velocities, making
sure that they satisfy the divergence constraint. These advective
velocities are used to construct the fluxes through the interfaces to
advance the solution to the new time. The second projection
operates on the cell-centered velocities at the new time, again
enforcing the divergence constraint.

AMReX-Hydro provides two projection classes: the ``MacProjector`` class
for face-centered velocity fields and ``NodalProjector`` for cell-centered
velocity fields. The projection classes use AMReX's linear solvers internally.
Both classes provide member functions ``getLinOp`` and ``getMLMG`` to
access the underlying objects and allow for modification of the linear operator
and multigrid properties if needed.
Details of the linear solver implementations are in the :ref:`amrex:Chap:LinearSolvers`
section of AMReX's documentation.

Both Projector classes provide the following parameters, which can be set in an
inputs file or on the command line. For the MacProjector, these must be preceeded by
"mac_proj.", or for the NodalProjector, "nodal_proj."

+-------------------+-----------------------------------------------------------------------+-------------+--------------+
|                   |  Description                                                          |   Type      | Default      |
+-------------------+-----------------------------------------------------------------------+-------------+--------------+
| verbose           |  Verbosity in nodal projection                                        |    Int      |   0          |
+-------------------+-----------------------------------------------------------------------+-------------+--------------+
| bottom_verbose    |  Verbosity of the bottom solver in nodal projection                   |    Int      |   0          |
+-------------------+-----------------------------------------------------------------------+-------------+--------------+
| maxiter           |  Maximum number of iterations                                         |    Int      |  MAC: 200    |
|                   |                                                                       |             |  Nodal: 100  |
+-------------------+-----------------------------------------------------------------------+-------------+--------------+
| bottom_maxiter    |  Maximum number of iterations in the bottom solver                    |    Int      |  MAC: 200    |
|                   |  if using bicg, cg, bicgcg or cgbicg                                  |             |  Nodal: 100  |
+-------------------+-----------------------------------------------------------------------+-------------+--------------+
| bottom_solver     |  Which bottom solver to use.                                          |  String     |   bicgcg     |
|                   |  Options are bicgcg, bicgstab, cg, cgbicg, smoother or hypre          |             |              |
+-------------------+-----------------------------------------------------------------------+-------------+--------------+
| bottom_rtol       |  Relative tolerance                                                   |   Real      |   1.0e-4     |
+-------------------+-----------------------------------------------------------------------+-------------+--------------+
| bottom_atol       |  Absolute tolerance, a negative number means it won't be used         |   Real      |   -1.0       |
+-------------------+-----------------------------------------------------------------------+-------------+--------------+
| num_pre_smooth    |  Number of smoother iterations when going down the V-cycle            |    Int      |   2          |
+-------------------+-----------------------------------------------------------------------+-------------+--------------+
| num_post_smooth   |  Number of smoother iterations when going up the V-cycle              |    Int      |   2          |
+-------------------+-----------------------------------------------------------------------+-------------+--------------+



.. _mac_proj:

MAC Projection
--------------

For a velocity field :math:`U = (u,v,w)` defined on faces, i.e.
:math:`u` is defined on x-faces, :math:`v` is defined on y-faces,
and :math:`w` is defined on z-faces, AMReX-Hydro provides an exact projection
we refer to as a MAC projection. For this we solve

.. math::

   D( \beta \nabla \phi) = D(U^*) - S

for :math:`\phi` and then set

.. math::

   U = U^* - \beta \nabla \phi


where :math:`U^*` is a vector field (typically velocity) that we want to satisfy
:math:`D(U) = S`.  For incompressible flow,  :math:`S = 0`.

The ``MacProjector`` class can be defined and used to perform the MAC projection without explicitly
calling the solver directly.  In addition to solving the Poisson equation (either variable or
constant coefficient),
the MacProjector internally computes the divergence of the vector field, :math:`D(U^*)`,
to compute the right-hand-side, and after the solve, subtracts the weighted gradient term to
make the vector field result satisfy the divergence constraint.

.. Note that passing ``nullptr`` for :math:`D(U^*)` is used for the MAC synchronization step in time-subcycling AMR (and more specifically, IAMR), where we want to solve for the correction velocity field which accounts for the mis-match in the advection velocity at the coarse-fine interface resulting from solving for the advection velocity on single levels rather than on the composite grid. In this case, currently, only the Poisson solve is done. Might make more sense to have MacProjector allocate and pass out -beta grad phi? or +beta grad phi?

In the simplest form of the call, :math:`S` is assumed to be zero and does not need to be specified.
Typically, the user does not allocate the solution array, but it is also possible to create and pass
in the solution array and have :math:`\phi` returned as well as :math:`U`.

The MacProjector class defaults to homogeneous Dirichlet or Neumann boundary conditions at domain
boundaries; for this case nothing further needs to be done.
Non-homogeneous Dirichlet or Neumann boundary conditions at domain boundaries are set with
member function ``void setLevelBC  (int amrlev, const amrex::MultiFab* levelbcdata)``.

If the MAC projection base level doesn't cover the full domain, one must pass boundary conditions
that come from coarser data with member function
``void setCoarseFineBC (const amrex::MultiFab* crse, int crse_ratio)``

The code below is taken from ``AMReX-Hydro/Tests/MAC_Projection_EB/main.cpp``,
and demonstrates how to set up the MACProjector object and use it to perform a MAC projection.

.. collapse:: Code Example - MacProjector object setup and MAC projection.

   .. code-block:: c++

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

      // Note that when we build with USE_EB = TRUE, we must specify whether the quantities are located
      //  at centers (MLMG::Location::CellCenter, MLMG::Location::FaceCenter) or
      //  centroids (MLMG::Location::CellCentroid, MLMG::Location::FaceCentroid).
      MacProjector macproj({amrex::GetArrOfPtrs(vel)},       // mac velocity
                           MLMG::Location::FaceCenter,       // Location of vel
                           {amrex::GetArrOfConstPtrs(beta)}, // beta
                           MLMG::Location::FaceCenter,       // Location of beta
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
      // Here, we specify that velocities are at face centers
      macproj.project(reltol,abstol,MLMG::Location::FaceCenter);

      // If we want to use phi elsewhere, we can pass in an array in which to return the solution
      // macproj.project({&phi_inout},reltol,abstol,MLMG::Location::FaceCenter);

|
|

.. _nodal_proj:

Nodal Projection
----------------

For a velocity field :math:`U = (u,v,w)` defined with all components co-located on cell centers,
AMReX-Hydro provides an approximate projection we refer to as a nodal projection.
Velocity divergence and pressure are defined on nodes, and the pressure gradient is defined
at cell centers as the cell average of face-based values. It is the use of this cell-averaged
pressure gradient that makes this projection approximate rather than exact.

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
to compute the approximate projection.  This means we must

The ``NodalProjector`` class can be used to solve the nodal projection without explicitly
calling the linear solver. In addtion to solving the nodal variable coefficient Poisson
equation, it internally computes the right-hand-side,
including the the divergence of the vector field, :math:`D(U^*)`,
and also subtracts the weighted gradient term to make the vector field result satisfy the
divergence constraint.

The NodalProjector class does not provide defaults for domain boundary conditions, and thus
member function ``void setLevelBC  (int amrlev, const amrex::MultiFab* levelbcdata)``
must always be called.

The code below is taken from ``AMReX-Hydro/Tests/Nodal_Projection_EB/main.cpp``,
and demonstrates how to set up the NodalProjector object and use it to perform a nodal projection.


.. collapse:: Example Code - NodalProjector object setup and nodal projection.

   .. code-block:: c++

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
      // Create the cell-centered velocity field we want to project.
      // Set velocity field to (1,0,0) including ghost cells for this example
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
      Hydro::NodalProjector nodal_proj({vel}, {sigma}, {geom}, lp_info, {rhs_cc}, {rhs_nd});

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

|
|

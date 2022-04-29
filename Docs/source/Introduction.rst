
.. _intro:

Introduction
============

AMReX-Hydro is set of routines that support the construction of convective
terms for incompressible and low Mach number flow modeling
in cartesian coordinates with (or without) embedded boundaries and R-Z coordinate systems.
It is not a stand-alone code, but is used in several application codes, such as
`incflo <https://amrex-codes.github.io/incflo/docs_html/>`_
(a variable density incompressible Navier-Stokes solver with adaptive mesh refinement (AMR)),
`IAMR <https://amrex-codes.github.io/IAMR/docs_html/index.html>`_
(a variable density incompressible Navier-Stokes solver with time subcycling AMR),
and `MFIX-Exa <https://amrex-codes.github.io/MFIX-Exa/docs_html/>`_
(a multiphase computational fluid dynamics modeling tool).

In application codes, the general procedure for constructing convective terms from cell-centered data
is as follows:

1. Construct values of the normal velocity at the centroid on each cell face using chosen advection scheme

2. MAC-project this face-based velocity field

3. Use MAC-projected velocities to help construct values on faces

4. Define fluxes on faces

5. Difference fluxes to create the convective term.

AMReX-Hydro provides routines to support all of these steps.
Here we group the AMReX-Hydro routines into a few general categories and map them to the step(s) they address:

* :ref:`Schemes`: the fundamental algorithm is either a Method-of-Lines (MOL) or Godunov approach
  (used in steps 1 and 3).

* :ref:`projections` frameworks:

  + MAC Projection - enforces a divergence condition on an edge-based velocity field (used in step 2).

  + Nodal Projection - can be used to compute an approximate projection of a cell-centered
    velocity field, with pressure and velocity divergence defined on nodes
    (not generally used as part of computing the convective term, but used in application codes to define a
    cell-centered velocity update that approximately obeys a divergence constraint).

* :ref:`Redistribution:Redistribution` schemes: to address the "small cell problem" associated with explicit cut
  cell algorithms (part of step 4 for problems with embedded boundaries).

* :ref:`utilities`: to do things like compute slopes, create fluxes from face-centered values, and
  create the convective term from fluxes (used in all but step 2).



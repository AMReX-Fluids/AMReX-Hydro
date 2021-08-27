
.. _amrex_hydro_indx:

AMReX-Hydro
===========

AMReX-Hydro is set of routines that support the construction of convective 
terms for incompressible and low Mach number flow modeling.  There are different
routines depending on whether the grid being operated on contains cut cells.

The fundamental algorithm is either a Method-of-Lines (MOL) or Godunov approach.
We use one of these to construct values of the normal velocity the centroid on each cell face,
then MAC-project this velocity field, then use these MAC-projected velocities to help
construct values on faces which will be used to define fluxes on faces that will
be differenced to create the convective term.

AMReX-Hydro also includes the flux redistribution and state redistribution algorithms
which are used to address the "small cell problem" associated with explicit cut
cell algorithms.

We organize the documentation here by chapters with the same headings as the
directories within AMReX-Hydro.


.. toctree::
   :maxdepth: 1
   :caption: Contents:

   EBGodunov
   EBMOL
   Godunov
   MOL
   Redistribution
   Slopes


API documentation can be found in the `Doxygen Technical Reference`_.

.. _`Doxygen Technical Reference`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/

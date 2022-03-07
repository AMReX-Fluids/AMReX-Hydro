.. _schemes:

Advection schemes
=================

In AMReX-Hydro, the fundamental algorithm is either a Method-of-Lines (MOL) or Godunov approach.
Each method provides functions for two separate usages:

1. **Construct values of the normal velocity at the centroid on each cell face, termed "Pre-MAC" or
   extrapolated velocity.**
   Typically, this velocity is later MAC projected before being used as the advective velocity.
   (Information on the MAC projection is in the :ref:`mac_proj` section.)
   
2. **Construct states on faces, termed "edge states" or "Post-MAC."**
   These are typically later used to make fluxes which are
   then differenced to create the advective term.
   (Information on how fluxes and the convective term are constructed from edge states is given in the
   :ref:`fluxes` and :ref:`advective_term` sections.)

How domain boundary conditions affect the computation of these pre- and post-MAC states is
the same for all advection methods, and is described in the :ref:`bcs` section.

Next, we provide some notation, and 
then we detail the available methods in EB-regular, as well as EB-aware form when available. 

.. note::

   If a cell and all of its neighbors have volume fraction of 1 (i.e. they
   are not cut or covered cells), the EB methodology will return exactly the same answer (to machine
   precision) as the non-EB methodology. Here we define neighbor to mean any cell that would be
   involved in the calculation of the face-based state, and the extent of the resulting neighborhood
   varies depending on the order of the slope used.


Notation
---------

We define :math:`\varepsilon = 1.e-8` in `Utils/hydro_constants.H <https://amrex-codes.github.io/amrex-hydro/Doxygen/html/group__Utilities.html#ga57d5ce9bc3bca16e249c611342f3c550>`_. This is a empirically determined constant that works well for flows where velocities are on the order of 1.
     
.. include:: MOL.rst

.. include:: Godunov.rst

.. include:: BDS.rst



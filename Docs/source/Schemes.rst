.. _schemes:

Advection schemes
^^^^^^^^^^^^^^^^^

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


First, we provide details that are common to all the advection schemes.
Then we detail the available methods.

Domain boundary conditions are described in the :ref:`bcs` section.
Note that the boundary conditions are imposed before the upwinding descirbed below.


Throughout all? advection schemes, we define :math:`\varepsilon = 1.e-8` in `hydro_constants.H`_

.. _`hydro_constants.H`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/group__Utilities.html#ga57d5ce9bc3bca16e249c611342f3c550

.. We define :math:`\varepsilon = 1.e-8` in **Utils / hydro_constants.H**


.. note::

   Note: if a cell and all of its immediate neighbors have volume fraction of 1 (i.e. they
   are not cut or covered cells), the EB methodology will return exactly the same answer (to machine
   precision) as the non-EB methodology.


.. include:: MOL.rst

.. include:: Godunov.rst

.. include:: BDS.rst



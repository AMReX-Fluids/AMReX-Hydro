.. include:: CustomCommands.rst

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

Domain boundary conditions affect the computation of these pre- and post-MAC states in
the same way for all advection methods, and this is described in the :ref:`bcs` section.
All schemes also use the same routines to construct fluxes and then the convective term.

Next, we provide notation, and
then detail the available advection schemes in EB-regular, as well as EB-aware form when available.

.. note::

   If a cell and all of its neighbors have volume fraction of 1 (i.e. they
   are not cut or covered cells), the EB methodology will return exactly the same answer (to machine
   precision) as the non-EB methodology. Here we define neighbor to mean any cell that would be
   involved in the calculation of the face-based state, and the extent of the resulting neighborhood
   varies depending on the order of the slope used.


Notation
---------

.. Maybe this would be more clear presented in a table

Here we use :math:`(i,j,k)` to denote cell centers (or centroids for EB),
and thus :math:`(i-\frac{1}{2},j,k)` denotes the lower x-face of the :math:`(i,j,k)`-th cell,
:math:`(i,j+\frac{1}{2},k)` denotes the upper y-face of the :math:`(i,j,k)`-th cell, etc.

Super- or subscript :math:`L` (for left) indicates a state that has been extrapolated from values at lower x indices.
:math:`R` (for right) indicates extrapolation higher x indices.
For example, for the x-face located at :math:`(i+1/2,j,k)`, :math:`L` indicates extrapolation from
the :math:`(i,j,k)`-th cell center/centroid, and
:math:`R` extrapolation from the :math:`(i+1,j,k)`-th cell center/centroid.

Similarly, for the y-dimension, :math:`F` (for forward) indicates a state that has been extrapolated from values at lower y indices and
:math:`B` (for back) indicates extrapolation from higher y indices.
And for the third dimension,
:math:`D` (for down) indicates a state that has been extrapolated from values at lower z indices.
:math:`U` (for up) indicates extrapolation from higher z indices.

:math:`\U^{MAC}` is the MAC-projected velocity at face centers (or centroids for EB).

We define :math:`\varepsilon = 1.e-8` in `Utils/hydro_constants.H <https://amrex-codes.github.io/amrex-hydro/Doxygen/html/group__Utilities.html#ga57d5ce9bc3bca16e249c611342f3c550>`_. This is an empirically determined constant that works well for flows where velocities are on the order of 1.


.. include:: MOL.rst

.. include:: Godunov.rst

.. include:: BDS.rst


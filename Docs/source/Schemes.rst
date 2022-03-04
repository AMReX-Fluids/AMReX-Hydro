.. _schemes:

Advection schemes
=================

In AMReX-Hydro, the fundamental algorithm is either a Method-of-Lines (MOL) or Godunov approach.
Each method provides functions for two separate usages:

1. Construct values of the normal velocity at the centroid on each cell face, termed "Pre-MAC" or
   extrapolated velocity. Typically, this velocity is MAC projected before being used as the advective velocity.

2. Construct states on faces, termed "edge states." These are typically later used to make fluxes which are
   then differenced to create the advective term.

The available methods are

.. toctree::
   :maxdepth: 1

   MOL
   Godunov
   BDS


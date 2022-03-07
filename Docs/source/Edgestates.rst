.. _edgestates:

Computing Values on Cell Faces
==============================

In AMReX-Hydro, the fundamental algorithm is either a Method-of-Lines (MOL) or Godunov approach.
Each method provides functions for two separate usages:

1. Construct values of the normal velocity at the centroid on each cell face, termed "Pre-MAC" or
   extrapolated velocity. Typically, this is MAC projected and then used as the advective velocity.

2. Construct states on faces, termed "edge states." These are typically later differenced to create the advective term.

The available methods are

.. toctree::
   :maxdepth: 1

   MOL
   Godunov
   BDS
   EBGodunov
   EBMOL


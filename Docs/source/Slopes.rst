.. _slopes:

Slopes
------

AMReX-Hydro includes implementations of several different slope routines along with options to apply limiters.
For cells where this calculation would involve all regular cells (i.e. no cut or covered cells),
there are second-order and fourth-order stencils.

For (EB)Godunov, the default is monotonicity-limited fourth-order slopes.
For (EB)MOL, the default is monotonicity-limited second-order slopes.
Default limiting is as described in Colella (1985) :cite:`colglaz`,
where limiting is done on each component of the velocity individually.


.. _EBslopes:

EB Slopes
---------

The procedure for problems with embedded boundaries
is detailed below, and attempts to use standard (non-EB) stencils wherever possible.

* First, AMReX-Hydro attempts to compute the slope of the desired order (4 for Godunov and 2 for MOL)
  using standard stencils and limiters.

* For cases where a fourth-order slope is desired but the stencil would require the use of cut cells
  (as happens in EBGodunov), the code next attempts to use the standard second-order slope and limiter.

* If the standard second-order slope calculation
  would require the use of cut cells, then the slope computation will use a least squares approach,
  involving a linear fit to the at-most 26 (or 8 in 2D) nearest neighbors, with the function
  going through the centroid of cell (i,j,k) to the face centroid. This does not assume that the
  cell centroids, where the data is assumed to live, are the same as cell centers.
  This least-squares slope is then multiplied by a limiter based on the work of Barth-Jespersen
  that enforces no new maxima or minima when the state is predicted to the face centroids.

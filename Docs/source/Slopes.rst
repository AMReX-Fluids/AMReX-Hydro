.. _slopes:

Slopes
------

AMReX-Hydro's implementation of the piecewise linear method provides several options, and
leverages slopes routines from AMReX.
For cells where this calculation would involve all regular cells (i.e. no cut or covered cells),
there are second-order and fourth-order stencils along with options to apply limiters.
Note that the piecewise parabolic and BDS methods have their own routines for formulating slopes
(and these are housed within AMReX-Hydro).

For (EB)Godunov, the default is monotonicity-limited fourth-order slopes.
Default limiting is as described in Colella (1985) :cite:`colglaz`,
where limiting is done on each component of the velocity individually.

For (EB)MOL, the default is monotonicity-limited second-order slopes.
Default is the second order Monotonized Central (MC)
limiter (van Leer, 1977 :cite:`vanleer`), where limiting is applied direction by direction.

.. The scheme is described below for the u-velocity.

   The limiter computes the slope at cell `i` by combining the left, central
   and right u-variation `du`:

   .. code:: shell

             du_l = u(i) - u(i-1)               = left variation
             du_c = 0.5 * ( u(i+1) - u(i-1) )   = central (umlimited) variation
             du_r = u(i+1) - u(i)               = right variation

             Finally, the u-variation at cell `i` is given by :

             .. code:: shell

                       du(i) = sign(du_c) * min(2|du_l|, |du_c|, 2|du_r|)) if du_l*du_r > 0
                       du(i) = 0                                           otherwise



NOTE ON BOUNDARY CONDITIONS:

When periodic or Neumann BCs are imposed, schemes can be applied
without any change since the ghost cells outside the domain are filled
by either periodicity or by extrapolation.

For Dirichlet BCs, the BC value stored in the first ghost cell outside the domain
is considered as located directly on the boundary,
despite the fact that it is stored in what is otherwise considered a cell-centered
array. We then utilize one-sided differencing schemes that
use ONLY values from inside or on the domain boundary.



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

.. _fluxes:


Computing Fluxes (`ComputeFluxes`_)
===================================

.. _`ComputeFluxes`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceHydroUtils.html#ab70f040557a658e70ba076c9d105bab7

AMReX-Hydro has the option to compute intesive or extensive, i.e. area-weighted, fluxes.
Extensive fluxes are always used for problems using R-Z geometry or embedded boundaries.

Intensive fluxes are computed from the advective velocity, :math:`\U^{MAC}`, and edge states, :math:`s`,
by defining

.. math::

   F_{i-\frac{1}{2},j,k}^{x,n+\frac{1}{2}} = u^{MAC}_{i-\frac{1}{2},j,k}\; s_{i-\frac{1}{2},j,k}^{n+\frac{1}{2}}

on all x-faces,

.. math::

   F_{i,j-\frac{1}{2},k}^{y,n+\frac{1}{2}} = v^{MAC}_{i,j-\frac{1}{2},k}\; s_{i,j-\frac{1}{2},k}^{n+\frac{1}{2}}

on all y-faces, and

.. math::

   F_{i,j,k-\frac{1}{2}}^{z,n+\frac{1}{2}} = w^{MAC}_{i,j,k-\frac{1}{2}}\; s_{i,j,k-\frac{1}{2}}^{n+\frac{1}{2}}

on all z-faces.

Extensive fluxes are computed as

.. math::
   :label: fluxebg-eq1

   F_{i-\frac{1}{2},j,k}^{x,n+\frac{1}{2}} = a_{i-\frac{1}{2},j,k} \; u^{MAC}_{i-\frac{1}{2},j,k} \; s_{i-\frac{1}{2},j,k}^{n+\frac{1}{2}}

on all x-faces,

.. math::
   :label: fluxebg-eq2

   F_{i,j-\frac{1}{2},k}^{y,n+\frac{1}{2}} = a_{i,j-\frac{1}{2},k} \; v^{MAC}_{i,j-\frac{1}{2},k} \; s_{i,j-\frac{1}{2},k}^{n+\frac{1}{2}}

on all y-faces,

.. math::
   :label: fluxebg-eq3

   F_{i,j,k-\frac{1}{2}}^{z,n+\frac{1}{2}} = a_{i,j,k-\frac{1}{2}} \; w^{MAC}_{i,j,k-\frac{1}{2}}\; s_{i,j,k-\frac{1}{2}}^{n+\frac{1}{2}}

on all z-faces.

where :math:`a_{i-\frac{1}{2},j,k}` is the area of the lower x-face of cell-centered element `(i,j,k)`, etc.



.. _EBfluxes:


Computing EB Fluxes (`EB_ComputeFluxes`_)
=========================================

.. _`EB_ComputeFluxes`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceHydroUtils.html#ab70f040557a658e70ba076c9d105bab7

The fluxes are computed from the edge states above by defining, e.g.,

.. math::
   :label: fluxebg-eq4

   F_{i-\frac{1}{2},j,k}^{x,n+\frac{1}{2}} = a_{i-\frac{1}{2},j,k} \; u^{MAC}_{i-\frac{1}{2},j,k} \; s_{i-\frac{1}{2},j,k}^{n+\frac{1}{2}} \; \Delta_y \; \Delta_z

on all x-faces with non-zero area fraction,

.. math::
   :label: fluxebg-eq5

   F_{i,j-\frac{1}{2},k}^{y,n+\frac{1}{2}} = a_{i,j-\frac{1}{2},k} \; v^{MAC}_{i,j-\frac{1}{2},k} \; s_{i,j-\frac{1}{2},k}^{n+\frac{1}{2}} \; \Delta_x \; \Delta_z

on all y-faces with non-zero area fraction, and

.. math::
   :label: fluxebg-eq6

   F_{i,j,k-\frac{1}{2}}^{z,n+\frac{1}{2}} = a_{i,j,k-\frac{1}{2}} \; w^{MAC}_{i,j,k-\frac{1}{2}}\; s_{i,j,k-\frac{1}{2}}^{n+\frac{1}{2}} \; \Delta_x \; \Delta_y

on all z-faces with non-zero area fraction.

Here :math:`\Delta_x, \Delta_y,` and :math:`\Delta_z` are the cell sizes in the 3 directions.

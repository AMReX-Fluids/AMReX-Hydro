.. include:: CustomCommands.rst

.. _fluxes:


Computing Fluxes
----------------

Doxygen links
`ComputeFluxes <https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceHydroUtils.html#ab70f040557a658e70ba076c9d105bab7>`_
and
`EB_ComputeFluxes <https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceHydroUtils.html#ab70f040557a658e70ba076c9d105bab7>`_ .

AMReX-Hydro has the option to compute intesive or extensive, i.e. area-weighted, fluxes.
Extensive fluxes are always used for problems using R-Z geometry,
and AMReX has functions for computing the position dependent cell volume and face areas
(see `Geometry::GetFaceArea <https://amrex-codes.github.io/amrex/doxygen/classamrex_1_1Geometry.html#a8e4aa2c2e88a46cb31b1ce0f9590c350>`_ and
`Geometry::GetVolume <https://amrex-codes.github.io/amrex/doxygen/classamrex_1_1Geometry.html#af7996fe47b1e82704565102c15df47c9>`_).
We first give the formulas for the EB-regular case,
and then those used when embedded boundaries are present.


Intensive
~~~~~~~~~

Intensive fluxes are computed from the advective velocity, :math:`\U^{MAC}`, and edge states, :math:`s`,
by defining

.. math::
   :label: flux-eq1

   F_{i-\frac{1}{2},j,k}^{x,n+\frac{1}{2}} = u^{MAC}_{i-\frac{1}{2},j,k}\; s_{i-\frac{1}{2},j,k}^{n+\frac{1}{2}}

on all x-faces,

.. math::
   :label: flux-eq2

   F_{i,j-\frac{1}{2},k}^{y,n+\frac{1}{2}} = v^{MAC}_{i,j-\frac{1}{2},k}\; s_{i,j-\frac{1}{2},k}^{n+\frac{1}{2}}

on all y-faces, and

.. math::
   :label: flux-eq3

   F_{i,j,k-\frac{1}{2}}^{z,n+\frac{1}{2}} = w^{MAC}_{i,j,k-\frac{1}{2}}\; s_{i,j,k-\frac{1}{2}}^{n+\frac{1}{2}}

on all z-faces.

|

When embedded boundaries are present, intensive fluxes are computed as

.. math::

   F_{i-\frac{1}{2},j,k}^{x,n+\frac{1}{2}} = \alpha_{i-\frac{1}{2},j,k} \; u^{MAC}_{i-\frac{1}{2},j,k} \; s_{i-\frac{1}{2},j,k}^{n+\frac{1}{2}}

on all x-faces,

.. math::

   F_{i,j-\frac{1}{2},k}^{y,n+\frac{1}{2}} = \alpha_{i,j-\frac{1}{2},k} \; v^{MAC}_{i,j-\frac{1}{2},k} \; s_{i,j-\frac{1}{2},k}^{n+\frac{1}{2}}

on all y-faces,

.. math::

   F_{i,j,k-\frac{1}{2}}^{z,n+\frac{1}{2}} = \alpha_{i,j,k-\frac{1}{2}} \; w^{MAC}_{i,j,k-\frac{1}{2}}\; s_{i,j,k-\frac{1}{2}}^{n+\frac{1}{2}}

on all z-faces.
Here :math:`\alpha_{i-\frac{1}{2},j,k}` is the area fraction of the lower x-face of cell-centered element :math:`(i,j,k)`, etc.


Extensive
~~~~~~~~~

Extensive fluxes are computed as

.. math::

   F_{i-\frac{1}{2},j,k}^{x,n+\frac{1}{2}} = \area_{i-\frac{1}{2},j,k} \; u^{MAC}_{i-\frac{1}{2},j,k} \; s_{i-\frac{1}{2},j,k}^{n+\frac{1}{2}}

on all x-faces,

.. math::

   F_{i,j-\frac{1}{2},k}^{y,n+\frac{1}{2}} = \area_{i,j-\frac{1}{2},k} \; v^{MAC}_{i,j-\frac{1}{2},k} \; s_{i,j-\frac{1}{2},k}^{n+\frac{1}{2}}

on all y-faces,

.. math::

   F_{i,j,k-\frac{1}{2}}^{z,n+\frac{1}{2}} = \area_{i,j,k-\frac{1}{2}} \; w^{MAC}_{i,j,k-\frac{1}{2}}\; s_{i,j,k-\frac{1}{2}}^{n+\frac{1}{2}}

on all z-faces.

where :math:`\area_{i-\frac{1}{2},j,k}` is the area of the lower x-face of cell-centered element :math:`(i,j,k)`, etc.

|

For EB, we simply scale area the by the area fraction in the above equations. For example, we use
:math:`\alpha_{i-\frac{1}{2},j,k} \area_{i-\frac{1}{2},j,k}` in place of :math:`\area_{i-\frac{1}{2},j,k}`, etc.


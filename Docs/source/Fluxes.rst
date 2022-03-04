.. include:: CustomCommands.rst

.. _fluxes:


Computing Fluxes
================

Doxygen links
`ComputeFluxes <https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceHydroUtils.html#ab70f040557a658e70ba076c9d105bab7>`_
and
`EB_ComputeFluxes <https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceHydroUtils.html#ab70f040557a658e70ba076c9d105bab7>`_ .

AMReX-Hydro has the option to compute intesive or extensive, i.e. area-weighted, fluxes.
Extensive fluxes are always used for problems using R-Z geometry. We first give the formulas for the EB-regular case,
and then those used when embedded boundaries are present.

Note that for problems with embedded boundaries, AMReX will carry area fractions.
Here we use :math:`af_{i-\frac{1}{2},j,k}` to denote the area fraction of the lower x-face of cell-centered element :math:`(i,j,k)`,
so that if :math:`\Delta_y,` and :math:`\Delta_z` are the cell sizes in the y- and z- directions, respectively, then
the total area of the :math:`(i-\half,j,k)` face is given by

.. math::

   a_{i-\frac{1}{2},j,k} = af_{i-\frac{1}{2},j,k} \Delta_y \Delta_z

The areas of y- and z-faces are given by cyclicly permutating  x, y, and z.

Intensive fluxes
~~~~~~~~~~~~~~~~

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

   F_{i-\frac{1}{2},j,k}^{x,n+\frac{1}{2}} = af_{i-\frac{1}{2},j,k} \; u^{MAC}_{i-\frac{1}{2},j,k} \; s_{i-\frac{1}{2},j,k}^{n+\frac{1}{2}}

on all x-faces,

.. math::

   F_{i,j-\frac{1}{2},k}^{y,n+\frac{1}{2}} = af_{i,j-\frac{1}{2},k} \; v^{MAC}_{i,j-\frac{1}{2},k} \; s_{i,j-\frac{1}{2},k}^{n+\frac{1}{2}}

on all y-faces,

.. math::

   F_{i,j,k-\frac{1}{2}}^{z,n+\frac{1}{2}} = af_{i,j,k-\frac{1}{2}} \; w^{MAC}_{i,j,k-\frac{1}{2}}\; s_{i,j,k-\frac{1}{2}}^{n+\frac{1}{2}}

on all z-faces.
Here :math:`af_{i-\frac{1}{2},j,k}` is the area fraction of the lower x-face of cell-centered element :math:`(i,j,k)`, etc.


Extensive Fluxes
~~~~~~~~~~~~~~~~

Extensive fluxes are computed as 

.. math::

   F_{i-\frac{1}{2},j,k}^{x,n+\frac{1}{2}} = a_{i-\frac{1}{2},j,k} \; u^{MAC}_{i-\frac{1}{2},j,k} \; s_{i-\frac{1}{2},j,k}^{n+\frac{1}{2}}

on all x-faces,

.. math::

   F_{i,j-\frac{1}{2},k}^{y,n+\frac{1}{2}} = a_{i,j-\frac{1}{2},k} \; v^{MAC}_{i,j-\frac{1}{2},k} \; s_{i,j-\frac{1}{2},k}^{n+\frac{1}{2}}

on all y-faces,

.. math::

   F_{i,j,k-\frac{1}{2}}^{z,n+\frac{1}{2}} = a_{i,j,k-\frac{1}{2}} \; w^{MAC}_{i,j,k-\frac{1}{2}}\; s_{i,j,k-\frac{1}{2}}^{n+\frac{1}{2}}

on all z-faces.

where :math:`a_{i-\frac{1}{2},j,k}` is the area of the lower x-face of cell-centered element :math:`(i,j,k)`, etc.

|

For EB, we simply use :math:`a_{i-\frac{1}{2},j,k} = af_{i-\frac{1}{2},j,k} \Delta_y \Delta_z`, etc., in the above equations,
where :math:`af_{i-\frac{1}{2},j,k}` is the area fraction of the lower x-face of cell-centered element :math:`(i,j,k)`, etc.


.. _advective_term:

Constructing the advective term
-------------------------------

AMReX-Hydro provides the option to compute the advective term either in
conservative (:math:`\grad \cdot \U s`) or convective form (:math:`\U \cdot \grad s`).
The exact equations differ slightly depending on whether the fluxes are
intensive (default for Cartesian geometries without embedded boundaries) or
extensive (EB or R-Z coordinates).

Intensive fluxes
~~~~~~~~~~~~~~~~
If the variable, :math:`s` is to be updated conservatively, we construct

.. math::

   \nabla \cdot \left({\bf u} s\right)^{n+\frac{1}{2}}
                             = & \frac{1}{dx} \left(F_{i+\frac{1}{2},j,k}^{x,n+\frac{1}{2}} -
                                  F_{i-\frac{1}{2},j,k}^{x,n+\frac{1}{2}}\right) + \\
                               & \frac{1}{dy} \left(F_{i,j+\frac{1}{2},k}^{y,n+\frac{1}{2}} -
                                  F_{i,j-\frac{1}{2},k}^{y,n+\frac{1}{2}}\right) + \\
                               & \frac{1}{dz} \left(F_{i,j,k+\frac{1}{2}}^{z,n+\frac{1}{2}} -
                                  F_{i,j,k-\frac{1}{2}}^{z,n+\frac{1}{2}}\right)

while if :math:`s` is to be updated in convective form, we construct

.. math::

   \left({\bf u}\cdot \nabla s\right)^{n+\frac{1}{2}} = \nabla \cdot \left({\bf u} s\right)^{n+\frac{1}{2}} - s_{i,j,k}^{n+\frac{1}{2}} \; \left(DU\right)^{MAC}

where

.. math::

   \left(DU\right)^{MAC} = \;
                   & \frac{1}{dx} \left(u^{MAC}_{i+\frac{1}{2},j,k} - u^{MAC}_{i-\frac{1}{2},j,k}\right) + \\
                   & \frac{1}{dy} \left(v^{MAC}_{i,j-\frac{1}{2},k} - v^{MAC}_{i,j-\frac{1}{2},k}\right) + \\
                   & \frac{1}{dz} \left(w^{MAC}_{i,j,k-\frac{1}{2}} - w^{MAC}_{i,j,k-\frac{1}{2}}\right)

and

.. math::

   s_{i,j,k}^{{n+\frac{1}{2}}} = \frac{1}{6} \left(
                    s_{i-\frac{1}{2},j,k}^{{n+\frac{1}{2}}} + s_{i+\frac{1}{2},j,k}^{{n+\frac{1}{2}}}
                +   s_{i,j-\frac{1}{2},k}^{{n+\frac{1}{2}}} + s_{i,j-\frac{1}{2},k}^{{n+\frac{1}{2}}}
                +   s_{i,j,k-\frac{1}{2}}^{{n+\frac{1}{2}}} + s_{i,j,k-\frac{1}{2}}^{{n+\frac{1}{2}}} \right)

|

The advective 

.. _EBadvective_term:

Extensive fluxes
~~~~~~~~~~~~~~~~

For problems with embedded boundaries or using R_Z coordinates, the fluxes are always extensive, and
so must be divided by the cell volume in constructing the advective term.
For R-Z systems, AMReX has functions for computing the position dependent cell volume and face areas
FIXME (need ref to AMReX here).
For EB sytstems, AMReX carries a volume fraction, :math:`Vf`, and area fractions, :math:`af`.
In this case, the total volume of the :math:`(i,j,k)`-th element is given by

.. math::

   V_{i,j,k} = Vf_(i,j,k) \Delta_x \Delta_y \Delta_z

and the total area of the lower x-face of cell-centered element :math:`(i,j,k)` is

.. math::

    a_{i-\frac{1}{2},j,k} = af_{i-\frac{1}{2},j,k} \Delta_y \Delta_z
   
where :math:`\Delta_x, \Delta_y,` and :math:`\Delta_z` are the cell sizes in the 3 dimensions,
and the y- and z-face areas are similarly formed with cyclical permutations of x, y, and z. 

If the variable, :math:`s` is to be updated conservatively, on all cells with :math:`Vf_{i,j,k} > 0` we construct

.. math::

    \nabla \cdot ({\bf u}s)^{n+\frac{1}{2}} = (
                           & ( F_{i+\frac{1}{2},j,k}^{{x,n+\frac{1}{2}}} -F_{i-\frac{1}{2},j,k}^{{x,n+\frac{1}{2}}}) + \\
                           & ( F_{i,j+\frac{1}{2},k}^{{y,n+\frac{1}{2}}} -F_{i,j-\frac{1}{2},k}^{{y,n+\frac{1}{2}}}) + \\
                           & ( F_{i,j,k+\frac{1}{2}}^{{z,n+\frac{1}{2}}} -F_{i,j,k-\frac{1}{2}}^{{z,n+\frac{1}{2}}}) ) / V_{i,j,k}

while if :math:`s` is to be updated in convective form, we construct

.. math::

   ({\bf u}\cdot \nabla s)^{n+\frac{1}{2}} = \nabla \cdot ({\bf u}s)^{n+\frac{1}{2}} - s_{i,j,k}^{{n+\frac{1}{2}}} (DU)^{MAC}

where

.. math::

   (DU)^{MAC}  = ( & (a_{i+\frac{1}{2},j,k} u^{MAC}_{i+\frac{1}{2},j,k}- a_{i-\frac{1}{2},j,k} u^{MAC}_{i-\frac{1}{2},j,k}) + \\
                   & (a_{i,j+\frac{1}{2},k} v^{MAC}_{i,j-\frac{1}{2},k}- a_{i,j-\frac{1}{2},k} v^{MAC}_{i,j-\frac{1}{2},k}) + \\
                   & (a_{i,j,k+\frac{1}{2}} w^{MAC}_{i,j,k-\frac{1}{2}}- a_{i,j,k-\frac{1}{2}} w^{MAC}_{i,j,k-\frac{1}{2}}) ) / V_{i,j,k}


.. math::

   s_{i,j,k}^{{n+\frac{1}{2}}} = (1/6) (
                    s_{i-\frac{1}{2},j,k}^{{n+\frac{1}{2}}} + s_{i+\frac{1}{2},j,k}^{{n+\frac{1}{2}}}
                +   s_{i,j-\frac{1}{2},k}^{{n+\frac{1}{2}}} + s_{i,j-\frac{1}{2},k}^{{n+\frac{1}{2}}}
                +   s_{i,j,k-\frac{1}{2}}^{{n+\frac{1}{2}}} + s_{i,j,k-\frac{1}{2}}^{{n+\frac{1}{2}}} )


where :math:`
etc., and, in terms of EB area fractions :math:`af`, is equivalent to :math:

|


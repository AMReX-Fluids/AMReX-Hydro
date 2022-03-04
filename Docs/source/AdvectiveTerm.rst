.. _advective_term:

Constructing the advective term
-------------------------------

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

.. _EBadvective_term:

Constructing the EB advective term
-------------------------------

If the variable, :math:`s` is to be updated conservatively, on all cells with :math:`V_{i,j,k} > 0` we construct

.. math::

    \nabla \cdot ({\bf u}s)^{n+\frac{1}{2}} = (
                           & ( F_{i+\frac{1}{2},j,k}^{{x,n+\frac{1}{2}}} -F_{i-\frac{1}{2},j,k}^{{x,n+\frac{1}{2}}}) + \\
                           & ( F_{i,j+\frac{1}{2},k}^{{y,n+\frac{1}{2}}} -F_{i,j-\frac{1}{2},k}^{{y,n+\frac{1}{2}}}) + \\
                           & ( F_{i,j,k+\frac{1}{2}}^{{z,n+\frac{1}{2}}} -F_{i,j,k-\frac{1}{2}}^{{z,n+\frac{1}{2}}}) ) / (V_{i,j,k} \Delta_x \Delta_y \Delta_z)

while if :math:`s` is to be updated in convective form, we construct

.. math::

   ({\bf u}\cdot \nabla s)^{n+\frac{1}{2}} = \nabla \cdot ({\bf u}s)^{n+\frac{1}{2}} - s_{i,j,k}^{{n+\frac{1}{2}}} (DU)^{MAC}

where

.. math::

   (DU)^{MAC}  = ( & (a_{i+\frac{1}{2},j,k} u^{MAC}_{i+\frac{1}{2},j,k}- a_{i-\frac{1}{2},j,k} u^{MAC}_{i-\frac{1}{2},j,k}) + \\
                   & (a_{i,j+\frac{1}{2},k} v^{MAC}_{i,j-\frac{1}{2},k}- a_{i,j-\frac{1}{2},k} v^{MAC}_{i,j-\frac{1}{2},k}) + \\
                   & (a_{i,j,k+\frac{1}{2}} w^{MAC}_{i,j,k-\frac{1}{2}}- a_{i,j,k-\frac{1}{2}} w^{MAC}_{i,j,k-\frac{1}{2}}) ) / V_{i,j,k}

and

.. math::

   s_{i,j,k}^{{n+\frac{1}{2}}} = (1/6) (
                    s_{i-\frac{1}{2},j,k}^{{n+\frac{1}{2}}} + s_{i+\frac{1}{2},j,k}^{{n+\frac{1}{2}}}
                +   s_{i,j-\frac{1}{2},k}^{{n+\frac{1}{2}}} + s_{i,j-\frac{1}{2},k}^{{n+\frac{1}{2}}}
                +   s_{i,j,k-\frac{1}{2}}^{{n+\frac{1}{2}}} + s_{i,j,k-\frac{1}{2}}^{{n+\frac{1}{2}}} )

|
|
|


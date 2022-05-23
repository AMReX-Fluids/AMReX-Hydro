.. include:: CustomCommands.rst

.. _advective_term:

Constructing the advective term
-------------------------------

AMReX-Hydro provides the option to compute the advective term either in
conservative (:math:`\grad \cdot \U s`) or convective form (:math:`\U \cdot \grad s`).
The exact equations differ slightly depending on whether the fluxes are
intensive or extensive (required for R-Z coordinates).

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


Extensive fluxes
~~~~~~~~~~~~~~~~

Given extensive fluxes we must divide by the cell volume in constructing the advective term.
Here, we give the form for EB and note that the EB-regular case is obtained by setting both :math:`\alpha = 1` and
:math:`\kappa=1`.
If the variable, :math:`s` is to be updated conservatively, on all cells with :math:`\kappa_{i,j,k} > 0` we construct

.. math::

    \nabla \cdot ({\bf u}s)^{n+\frac{1}{2}} = \frac{1}{ \kappa_{i,j,k} \vol_{i,j,k}}[
                           & ( F_{i+\frac{1}{2},j,k}^{{x,n+\frac{1}{2}}} -F_{i-\frac{1}{2},j,k}^{{x,n+\frac{1}{2}}}) + \\
                           & ( F_{i,j+\frac{1}{2},k}^{{y,n+\frac{1}{2}}} -F_{i,j-\frac{1}{2},k}^{{y,n+\frac{1}{2}}}) + \\
                           & ( F_{i,j,k+\frac{1}{2}}^{{z,n+\frac{1}{2}}} -F_{i,j,k-\frac{1}{2}}^{{z,n+\frac{1}{2}}}) ]

while if :math:`s` is to be updated in convective form, we construct

.. math::

   ({\bf u}\cdot \nabla s)^{n+\frac{1}{2}} = \nabla \cdot ({\bf u}s)^{n+\frac{1}{2}} - s_{i,j,k}^{{n+\frac{1}{2}}} (DU)^{MAC}

where

.. math::

   (DU)^{MAC}  = \frac{1}{\kappa_{i,j,k} \vol_{i,j,k}} [ & (\alpha_{i+\frac{1}{2},j,k} \area_{i+\frac{1}{2},j,k} u^{MAC}_{i+\frac{1}{2},j,k}-
   \alpha_{i-\frac{1}{2},j,k} \area_{i-\frac{1}{2},j,k} u^{MAC}_{i-\frac{1}{2},j,k}) + \\
                   & (\alpha_{i,j+\frac{1}{2},k} \area_{i,j+\frac{1}{2},k} v^{MAC}_{i,j-\frac{1}{2},k}-
           \alpha_{i,j-\frac{1}{2},k} \area_{i,j-\frac{1}{2},k} v^{MAC}_{i,j-\frac{1}{2},k}) + \\
                   & (\alpha_{i,j,k+\frac{1}{2}} \area_{i,j,k+\frac{1}{2}} w^{MAC}_{i,j,k-\frac{1}{2}}-
           \alpha_{i,j,k-\frac{1}{2}} \area_{i,j,k-\frac{1}{2}} w^{MAC}_{i,j,k-\frac{1}{2}}) ]

.. math::

   s_{i,j,k}^{{n+\frac{1}{2}}} = (1/6) (
                    s_{i-\frac{1}{2},j,k}^{{n+\frac{1}{2}}} + s_{i+\frac{1}{2},j,k}^{{n+\frac{1}{2}}}
                +   s_{i,j-\frac{1}{2},k}^{{n+\frac{1}{2}}} + s_{i,j-\frac{1}{2},k}^{{n+\frac{1}{2}}}
                +   s_{i,j,k-\frac{1}{2}}^{{n+\frac{1}{2}}} + s_{i,j,k-\frac{1}{2}}^{{n+\frac{1}{2}}} )





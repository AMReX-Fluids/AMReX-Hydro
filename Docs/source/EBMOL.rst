EBMOL
=====

The procedure for computing MAC velocities and edge states with EB-aware MOL 
does not involve any time derivatives. All slope computations use
second-order limited slopes as described in
`Slopes`_.

.. _`Slopes`: https://amrex-codes.github.io/amrex/hydro_html/Slopes.html

.. note::
   
   Note: if a cell and all of its immediate neighbors have volume fraction of 1 (i.e. they 
   are not cut cells), the EBMOL methodology will return exactly the same answer (to machine
   precision) as the MOL methodology.

We define :math:`\varepsilon = 1.e-8` in **Utils / hydro_constants.H**

Notation
--------

For every cut cell we define :math:`a_x`, :math:`a_y,` and :math:`a_z` to be the area fractions of the faces 
and :math:`V` is the volume fraction of the cell.  All area and volume fractions are greater than or equal to zero 
and less than or equal to 1.

Pre-MAC (`ExtrapVelToFaces`_)
---------------------------------------

.. _`ExtrapVelToFaces`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceEBMOL.html#a7add53a153ade9c5cb83e79a61ad1929

For computing the pre-MAC edge states to be MAC-projected, we define on every x-face with :math:`a_x > 0` :

.. math::

   \begin{aligned}
   u_L &=& u_{i-1,j,k} + \delta x \; {u^x}_{i-1,j,k} + \delta y \; {u^y}_{i-1,j,k} + \delta z \; {u^z}_{i-1,j,k} , \\
   u_R &=& u_{i,j,k}   - \delta x \; {u^x}_{i,j,k}   - \delta y \; {u^y}_{i,j,k}   - \delta z \; {u^z}_{i,j,k} ,\end{aligned}

where we calculate :math:`u^x`, :math:`u^y` and :math:`u^z` simultaneously using a least squares approach,
as described in `Slopes`_,
and :math:`\delta_x`, :math:`\delta_y`, and :math:`\delta_z` are the components of the distance vector from 
the cell centroid to the face centroid of the face at :math:`(i-\frac{1}{2},j,k).`

At each face we then upwind based on :math:`u_L` and :math:`u_R`

.. math::

   u_{i-\frac{1}{2},j,k} = 
   \begin{cases}
   0, & \mathrm{if} \; u_L < 0 \;\; \mathrm{and} \;\; u_R > 0 \; \mathrm{else} \\
   u_L, & \mathrm{if} \; u_L + u_R \ge  \varepsilon  \; \mathrm{else} \\
   u_R, & \mathrm{if} \; u_L + u_R \le  -\varepsilon  \; \mathrm{else} \\
   0
   \end{cases}

We similarly compute :math:`v_{i,j-\frac{1}{2},k}` on y-faces and
:math:`w_{i,j,k-\frac{1}{2}}` on z-faces.

Boundary conditions (`SetXEdgeBCs`_, `SetYEdgeBCs`_, `SetZEdgeBCs`_) 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _`SetXEdgeBCs`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceHydroBC.html#ab90f8ce229a7ebbc521dc27d65f2db9a
.. _`SetYEdgeBCs`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceHydroBC.html#a6865c2cfd50cc95f9b69ded1e8ac78ab
.. _`SetZEdgeBCs`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceHydroBC.html#a19ddc5ac50e9a6b9a98bc17f3815a62e

Domain boundary conditions affect the above in three ways.

(1) First, they potentially impact the slope computation in cells
adjacent to the domain boundary (see `Slopes`_).

(2) Second, if the face is on a domain boundary and the boundary
condition type is extdir, we set both :math:`u_L` and :math:`u_R` to the
boundary value. If the boundary condition type is foextrap, hoextrap, or 
reflecteven on the low side of the domain, 
we set :math:`u_L = u_R.` (If on the high side then we
set :math:`u_R = u_L.`) If the boundary condition type is reflectodd , we set
:math:`u_L = u_R = 0.`

(3) In addition, if the domain boundary condition on the low side is foextrap
or hoextrap, we set :math:`u_L = u_R = \min (u_R, 0).` If the domain boundary
condition on the high side is foextrap or hoextrap, we set
:math:`u_L = u_R = \max (u_L, 0).` This has the effect of not allowing
the velocity at an outflow face to flow back into the domain.

Note that the boundary conditions are imposed before the upwinding
described above.

Post-MAC (`ComputeEdgeState`_)
------------------------------

.. _`ComputeEdgeState`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceEBMOL.html#a94df1b279b45eac5141dfe0dff0a79bc

Once we have the MAC-projected velocities, we predict all quantities to faces with non-zero area fractions as above:

.. math::

   \begin{aligned}
   s_L &=& s_{i-1,j,k} + \delta x \; {s^x}_{i-1,j,k} + \delta y \; {s^y}_{i-1,j,k} + \delta z \; {s^z}_{i-1,j,k} , \\
   s_R &=& s_{i,j,k}   - \delta x \; {s^x}_{i,j,k}   - \delta y \; {s^y}_{i,j,k}   - \delta z \; {s^z}_{i,j,k} ,\end{aligned}

where we calculate :math:`s^x`, :math:`s^y` and :math:`s^z` simultaneously using a least squares approach,
as described in `Slopes`_,
and :math:`\delta_x`, :math:`\delta_y`, and :math:`\delta_z` are the components of the distance vector from 
the cell centroid to the face centroid of the face at :math:`(i-\frac{1}{2},j,k).`

The domain boundary conditions affect the solution as described above in
(1) and (2) for the pre-MAC step. We do not impose the
no-outflow-at-inflow condition quite as described in (3); instead we
impose that if, on the low side, :math:`u^{MAC}\ge 0` (i.e the flow is
coming in at an outflow face) and :math:`s` is the x-velocity, then
:math:`s_L = s_R = \min(s_R,0).` On the high side, if
:math:`u^{MAC}<= 0` on the domain face, then
:math:`s_L = s_R = \max(s_L,0).` This enforces that if :math:`u^{MAC}`
on an outflow face is inflowing, the normal velocity component must be
outflowing or zero.

At each face we then upwind based on :math:`u^{MAC}_{i-\frac{1}{2},j,k}`

.. math::

   s_{i-\frac{1}{2},j,k} = 
   \begin{cases}
   s_L, & \mathrm{if} \; u^{MAC}_{i-\frac{1}{2},j,k}\; \ge  \; \varepsilon  \; \mathrm{else} \\
   s_R, & \mathrm{if} \; u^{MAC}_{i-\frac{1}{2},j,k}\; \le  \; -\varepsilon  \; \mathrm{else} \\
   \frac{1}{2}(s_L + s_R), 
   \end{cases}

Computing the Fluxes (`ComputeFluxes`_)
---------------------------------------

.. _`ComputeFluxes`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceHydroUtils.html#ab70f040557a658e70ba076c9d105bab7

The fluxes are computed from the edge states above by defining, e.g.,

.. math::

   F_{i-\frac{1}{2},j,k}^{x,n} = a_{i-\frac{1}{2},j,k} \; u^{MAC}_{i-\frac{1}{2},j,k}\; s_{i-\frac{1}{2},j,k}^{n}

on all x-faces with non-zero area fraction,

.. math::

   F_{i,j-\frac{1}{2},k}^{y,n} = a_{i,j-\frac{1}{2},k} \; v^{MAC}_{i,j-\frac{1}{2},k}\; s_{i,j-\frac{1}{2},k}^{n}

on all y-faces with non-zero area fraction, and

.. math::

   F_{i,j,k-\frac{1}{2}}^{z,n} = a_{i,j,k-\frac{1}{2}} \; w^{MAC}_{i,j,k-\frac{1}{2}}\; s_{i,j,k-\frac{1}{2}}^{n}

on all z-face with non-zero area fraction

Constructing the update
-----------------------

If the variable, :math:`s` is to be updated conservatively, on all cells with :math:`V_{i,j,k} > 0` we construct

.. math::

   \nabla \cdot ({\bf u}s)^n  = (
                              & (F_{i+\frac{1}{2},j,k}^{x,n} - F_{i-\frac{1}{2},j,k}^{x,n}) + \\
                              & (F_{i,j+\frac{1}{2},k}^{y,n} - F_{i,j-\frac{1}{2},k}^{y,n}) + \\
                              & (F_{i,j,k+\frac{1}{2}}^{z,n} - F_{i,j,k-\frac{1}{2}}^{z,n}) ) / V_{i,j,k}

while if :math:`s` is to be updated in convective form, we construct

.. math::

   ({\bf u}\cdot \nabla s)^n = \nabla \cdot ({\bf u}s)^n - s_{i,j,k}^n (DU)^{MAC}

where

.. math::

   (DU)^{MAC}  = (
                 & (a_{i+\frac{1}{2},j,k} u^{MAC}_{i+\frac{1}{2},j,k}- a_{i-\frac{1}{2},j,k} u^{MAC}_{i-\frac{1}{2},j,k}) + \\
                 & (a_{i,j+\frac{1}{2},k} v^{MAC}_{i,j-\frac{1}{2},k}- a_{i,j-\frac{1}{2},k} v^{MAC}_{i,j-\frac{1}{2},k}) + \\
                 & (a_{i,j,k+\frac{1}{2}} w^{MAC}_{i,j,k-\frac{1}{2}}- a_{i,j,k-\frac{1}{2}} w^{MAC}_{i,j,k-\frac{1}{2}}) ) / V_{i,j,k}

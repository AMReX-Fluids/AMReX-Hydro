Method-of-Lines
~~~~~~~~~~~~~~~

MOL
===

The procedure for computing MAC velocities and edge states with MOL does
not involve any time derivatives. All slope computations use
second-order limited slopes as described in
`[sec:slopes] <#sec:slopes>`__.

We define :math:`\varepsilon = 1.e-8` in **Utils / hydro_constants.H**

Pre-MAC (**MOL::ExtrapVelToFacesBox** )
---------------------------------------

For computing the pre-MAC edge states to be MAC-projected, we define on
every x-face:

.. math::

   \begin{aligned}
   u_L &=& u_{i-1,j,k} + \half \Delta u_{\i-1,j,k}^x, \\
   u_R &=& u_{i,j,k}   - \half \Delta u_{i,j,k}^x,\end{aligned}

At each face we then upwind based on :math:`u_L` and :math:`u_R`

.. math::

   u_{\imh,j,k} = 
   \begin{cases}
   0, & \mathrm{if} \; u_L < 0 \;\; \mathrm{and} \;\; u_R > 0 \; \mathrm{else} \\
   u_L, & \mathrm{if} \; u_L + u_R \ge  \varepsilon  \; \mathrm{else} \\
   u_R, & \mathrm{if} \; u_L + u_R \le  -\varepsilon  \; \mathrm{else} \\
   0
   \end{cases}

We similarly compute :math:`v_{i,\jmh,k}` on y-faces and
:math:`w_{i,j,\kmh}` on z-faces.

Effect of boundary conditions (**SetXEdgeBCs** in **Utils / hydro_bcs_K.H** )
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Domain boundary conditions affect the above in two ways.

(1) First, they potentially impact the slope computation in cells
adjacent to the domain boundary (see `[sec:slopes] <#sec:slopes>`__).

(2) Second, if the face is on a domain boundary and the boundary
condition type is , we set both :math:`u_L` and :math:`u_R` to the
boundary value. If the boundary condition type is , or on the low side
of the domain, we set :math:`u_L = u_R.` (If on the high side then we
set :math:`u_R = u_L.`) If the boundary condition type is , we set
:math:`u_L = u_R = 0.`

(3) In addition, if the domain boundary condition on the low side is or
, we set :math:`u_L = u_R = \min (u_R, 0).` If the domain boundary
condition on the high side is or , we set
:math:`u_L = u_R = \max (u_L, 0).` This has the effect of not allowing
the velocity at an outflow face to flow back into the domain.

Note that the boundary conditions are imposed before the upwinding
described above.

MAC projection
--------------

We then perform a MAC projection on the face-centered velocities to enforce that they satisfy 

.. math:: \nabla \cdot (U^{MAC})  = S

where S is a specified divergence constraint (0 for incompressible flow).

We do this by solving 

.. math:: \nabla \cdot \frac{1}{\rho} \nabla \phi^{MAC} = \nabla \cdot \left(U^{pred} \right)

then defining

.. math:: U^{MAC} = U^{pred} - \frac{1}{\rho} \nabla \phi^{MAC}

Post-MAC
--------

Once we have the MAC-projected velocities, we project all quantities to
faces as above:

.. math::

   \begin{aligned}
   s_L &=& s_{i-1,j,k} + \half \Delta s_{i-1,j,k}^x, \\
   s_R &=& s_{i,j,k}   - \half \Delta s_{i,j,k}^x,\end{aligned}

The domain boundary conditions affect the solution as described above in
(1) and (2) for the pre-MAC step. We do not impose the
no-outflow-at-inflow condition quite as described in (3); instead we
impose that if, on the low side, :math:`\umac \ge 0` (i.e the flow is
coming in at an outflow face) and :math:`s` is the x-velocity, then
:math:`s_L = s_R = \min(s_R,0).` On the high side, if :math:`\umac <= 0`
on the domain face, then :math:`s_L = s_R = \max(s_L,0).` This enforces
that if :math:`\umac` on an outflow face is inflowing, the normal
velocity component must be outflowing or zero.

At each face we then upwind based on :math:`\umaclo`

.. math::

   s_{\imh,j,k} = 
   \begin{cases}
   s_L, & \mathrm{if} \; \umaclo \; \ge  \; \varepsilon  \; \mathrm{else} \\
   s_R, & \mathrm{if} \; \umaclo \; \le  \; -\varepsilon  \; \mathrm{else} \\
   \half (s_L + s_R), 
   \end{cases}

Constructing the update
=======================

If the variable, :math:`s` is to be updated conservatively, we construct

.. math::

   \begin{aligned}
   \nabla \cdot (\ub s) &=& (\umachi \; s_{\iph,j,k} - \umaclo \; s_{\imh,j,k}) \\
                        &+& (\vmaclo \; s_{i,\jph,k} - \vmaclo \; s_{i,\jmh,k} ) \\
                        &+& (\wmaclo \; s_{i,j,\kph} - \wmaclo \; s_{i,j,\kmh}) \end{aligned}

while if :math:`s` is to be updated in convective form, we construct

.. math::

   \begin{aligned}
   (\ub \cdot \nabla s) &=& (\umachi \; s_{\iph,j,k} - \umaclo \; s_{\imh,j,k}) \\
                       &+& (\vmaclo \; s_{i,\jph,k} - \vmaclo \; s_{i,\jmh,k} ) \\
                       &+& (\wmaclo \; s_{i,j,\kph} - \wmaclo \; s_{i,j,\kmh}) \\
                       &-& s_{i,j,k} \; (DU)^{MAC}\end{aligned}

where

.. math::

   \begin{aligned}
   (DU)^{MAC}  &=& (\umachi - \umaclo ) \\
               &+& (\vmaclo - \vmaclo ) \\
               &+& (\wmaclo - \wmaclo ) \\\end{aligned}

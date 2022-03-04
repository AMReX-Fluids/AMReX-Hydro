.. _mol:

MOL
-----

The procedure for computing MAC velocities and edge states with MOL does
not involve any time derivatives. All slope computations use
second-order limited slopes as described in :ref:`slopes`.

Domain boundary conditions are described in the :ref:`bcs` section.
Note that the boundary conditions are imposed before the upwinding descirbed below.

We define :math:`\varepsilon = 1.e-8` in **Utils / hydro_constants.H**

These alogrithms are applied in the MOL namespace. For API documentation, see
`Doxygen: MOL Namespace`_.

.. _`Doxygen: MOL Namespace`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceMOL.html


Pre-MAC (`ExtrapVelToFaces`_)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _`ExtrapVelToFaces`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceMOL.html#acdde2acf756048b8ef0bca332e4bf748

For computing the pre-MAC edge states to be MAC-projected, we define on
every x-face:

.. math::

   \begin{aligned}
   u_L &=& u_{i-1,j,k} + \frac{\Delta x}{2} {u_x}_{i-1,j,k}, \\
   u_R &=& u_{i,j,k}   - \frac{\Delta x}{2} {u_x}_{i,j,k}, \end{aligned}

where :math:`u^x` are the (limited) slopes in the x-direction.

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
v

Post-MAC (`ComputeEdgeState`_)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _`ComputeEdgeState`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceMOL.html#acdde2acf756048b8ef0bca332e4bf748

Once we have the MAC-projected velocities, we project all quantities to
faces as above:

.. math::

   \begin{aligned}
   s_L &=& s^{i-1,j,k} + \frac{\Delta x}{2} {s^x}_{i-1,j,k}, \\
   s_R &=& s^{i,j,k}   - \frac{\Delta x}{2} {s^x}_{i,j,k},   \end{aligned}

where :math:`s^x` are the (limited) slopes in the x-direction.

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



Fluxes and Convective Term
~~~~~~~~~~~~~~~~~~~~~~~~~~

The values on cell faces (or edge states) can be used to construct fluxes and then a convective term.
Details on this are in the :ref:`fluxes` and :ref:`advective_term` sections. 


.. _ebmol:

EBMOL
-----

AMReX-Hydro has also implemented an embedded boundary (EB) aware version of the MOL algorithm
discussed above.
The procedure for computing MAC velocities and edge states with EB-aware MOL
does not involve any time derivatives.

All slope computations use
second-order limited slopes as described in :ref:`EBslopes`.

Domain boundary conditions are described in the :ref:`bcs` section.
Note that the boundary conditions are imposed before the upwinding descirbed below.


.. note::

   Note: if a cell and all of its immediate neighbors have volume fraction of 1 (i.e. they
   are not cut or covered cells), the EBMOL methodology will return exactly the same answer (to machine
   precision) as the MOL methodology.

We define :math:`\varepsilon = 1.e-8` in **Utils / hydro_constants.H**

Notation
~~~~~~~~

For every cut cell we define :math:`a_x`, :math:`a_y,` and :math:`a_z` to be the area fractions of the faces
and :math:`V` is the volume fraction of the cell.  All area and volume fractions are greater than or equal to zero
and less than or equal to 1.

Pre-MAC (`ExtrapVelToFaces`_)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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


Post-MAC (`ComputeEdgeState`_)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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


Fluxes and Convective Term
~~~~~~~~~~~~~~~~~~~~~~~~~~

The values on cell faces (or edge states) can be used to construct fluxes and then a convective term.
Details on this are in the :ref:`EBfluxes` and :ref:`EBadvective_term` sections. 

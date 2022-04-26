.. _mol:

Method of Lines (MOL)
---------------------

The procedure for computing MAC velocities and edge states with MOL involves extrapolation in space only,
and does not involve any time derivatives. All slope computations use
second-order limited slopes as described in :ref:`slopes`.

These alogrithms are applied in the MOL namespace. For API documentation, see
`Doxygen: MOL Namespace`_.

.. _`Doxygen: MOL Namespace`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceMOL.html


Pre-MAC (API ref. `MOL::ExtrapVelToFaces <https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceMOL.html#acdde2acf756048b8ef0bca332e4bf748>`_)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For computing the pre-MAC edge states to be MAC-projected, we define on
every x-face:

.. math::

   \begin{aligned}
   u_L &=& u_{i-1,j,k} + \frac{\Delta x}{2} {u^x}_{i-1,j,k}, \\
   u_R &=& u_{i,j,k}   - \frac{\Delta x}{2} {u^x}_{i,j,k}, \end{aligned}

where :math:`u^x` are the (limited) slopes in the x-direction.

Boundary conditions are applied (as decribed in :ref:`bcs`).
Then, at each face we upwind based on :math:`u_L` and :math:`u_R`

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


Post-MAC (API ref. `MOL::ComputeEdgeState <https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceMOL.html#acdde2acf756048b8ef0bca332e4bf748>`_)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once we have the MAC-projected velocities, we extrapolate all quantities to
faces as above:

.. math::

   \begin{aligned}
   s_L &=& s^{i-1,j,k} + \frac{\Delta x}{2} {s^x}_{i-1,j,k}, \\
   s_R &=& s^{i,j,k}   - \frac{\Delta x}{2} {s^x}_{i,j,k},   \end{aligned}

where :math:`s^x` are the (limited) slopes in the x-direction.

Boundary conditions are applied (as decribed in :ref:`bcs`).
Then, at each face, we upwind based on :math:`u^{MAC}_{i-\frac{1}{2},j,k}`

.. math::

   s_{i-\frac{1}{2},j,k} =
   \begin{cases}
   s_L, & \mathrm{if} \; u^{MAC}_{i-\frac{1}{2},j,k}\; \ge  \; \varepsilon  \; \mathrm{else} \\
   s_R, & \mathrm{if} \; u^{MAC}_{i-\frac{1}{2},j,k}\; \le  \; -\varepsilon  \; \mathrm{else} \\
   \frac{1}{2}(s_L + s_R),
   \end{cases}



.. _ebmol:

Method of Lines with Embedded Boundaries (EBMOL)
------------------------------------------------

AMReX-Hydro has also implemented an embedded boundary (EB) aware version of the MOL algorithm
discussed above.
All slope computations use second-order limited slopes as described in :ref:`EBslopes`.


Pre-MAC (API ref. `EBMOL::ExtrapVelToFaces <https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceEBMOL.html#a7add53a153ade9c5cb83e79a61ad1929>`_)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For computing the pre-MAC edge states to be MAC-projected, we define on every x-face with non-zero area fraction:

.. math::

   \begin{aligned}
   u_L &=& u_{i-1,j,k} + \delta_x \; {u^x}_{i-1,j,k} + \delta_y \; {u^y}_{i-1,j,k} + \delta z \; {u^z}_{i-1,j,k} , \\
   u_R &=& u_{i,j,k}   - \delta_x \; {u^x}_{i,j,k}   - \delta_y \; {u^y}_{i,j,k}   - \delta z \; {u^z}_{i,j,k} ,\end{aligned}

where we calculate :math:`u^x`, :math:`u^y` and :math:`u^z` as described in :ref:`EBslopes`,
and :math:`\delta_x`, :math:`\delta_y`, and :math:`\delta_z` are the components of the distance vector from
the cell centroid to the face centroid of the face at :math:`(i-\frac{1}{2},j,k).`

Boundary conditions are applied (as decribed in :ref:`bcs`).
Then, at each face we upwind based on :math:`u_L` and :math:`u_R`

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


Post-MAC (API ref. `EBMOL::ComputeEdgeState <https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceEBMOL.html#a94df1b279b45eac5141dfe0dff0a79bc>`_)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once we have the MAC-projected velocities, we predict all quantities to faces with non-zero area fractions as above:

.. math::

   \begin{aligned}
   s_L &=& s_{i-1,j,k} + \delta_x \; {s^x}_{i-1,j,k} + \delta_y \; {s^y}_{i-1,j,k} + \delta z \; {s^z}_{i-1,j,k} , \\
   s_R &=& s_{i,j,k}   - \delta_x \; {s^x}_{i,j,k}   - \delta_y \; {s^y}_{i,j,k}   - \delta z \; {s^z}_{i,j,k} ,\end{aligned}

where we calculate :math:`s^x`, :math:`s^y` and :math:`s^z` as described in :ref:`EBslopes`,
and :math:`\delta_x`, :math:`\delta_y`, and :math:`\delta_z` are the components of the distance vector from
the cell centroid to the face centroid of the face at :math:`(i-\frac{1}{2},j,k).`

Boundary conditions are applied (as decribed in :ref:`bcs`).
Then, at each face we then upwind based on :math:`u^{MAC}_{i-\frac{1}{2},j,k}`

.. math::

   s_{i-\frac{1}{2},j,k} =
   \begin{cases}
   s_L, & \mathrm{if} \; u^{MAC}_{i-\frac{1}{2},j,k}\; \ge  \; \varepsilon  \; \mathrm{else} \\
   s_R, & \mathrm{if} \; u^{MAC}_{i-\frac{1}{2},j,k}\; \le  \; -\varepsilon  \; \mathrm{else} \\
   \frac{1}{2}(s_L + s_R),
   \end{cases}


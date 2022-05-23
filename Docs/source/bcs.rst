.. _bcs:


Boundary conditions
-------------------

AMReX-Hydro uses underlying AMReX functionality in implementing boundary conditions
(see AMReX's documentation section :ref:`amrex:sec:basics:boundary`).
Physical boundary conditions, such as
inflow, outflow, slip/no-slip walls, etc., and are ultimately linked to
mathematical Dirichlet or Neumann conditions.
See ``amrex/Src/Base/AMReX_BC_TYPES.H`` for common physical and mathematical types.

Domain boundary conditions affect the pre-MAC extrapolated velocities in three ways.

#. Potential impact to the slope computation in cells
   adjacent to the domain boundary (see :ref:`slopes` section).

#. Direct enforcement of the boundary condition: If the face is on a domain boundary and the boundary
   condition type is

   * External Dirichlet (``extdir``): we set :math:`u_L` to the boundary value, and then
     for the normal component of the velocity only, we set :math:`u_R = u_L`
     This is done because for turbulent inflow, there can be times when the inflow face
     actually has outflowing velocity. In this case, we want to use the normal component as
     specified by the BC, but then allow that outflowing velocity to transport values that come
     from the interior.

   * First-order extrapolation (``foextrap``), higher order extrapolation (``hoextrap``), or
     even reflection about the boundary (``reflecteven``):

     + on the low side of the domain, we set :math:`u_L = u_R.`

     + on the high side, we set :math:`u_R = u_L.`

   * Odd reflection about the boundary (``reflectodd``) , we set :math:`u_L = u_R = 0.`

#. To prohibit back flow into the domain at an outflow face (``foextrap`` or ``hoextrap`` mathematical BCs):

   * on the low side, we set :math:`u_L = u_R = \min (u_R, 0).`

   * on the high side, we set :math:`u_L = u_R = \max (u_L, 0).`

.. What about Godunov trans term bcs???

For the post-MAC edge state,

#. Same as pre-MAC

#. Same as pre-MAC

#. Here, we do not impose the _`no-inflow-at-outflow` condition quite as described above;
   instead we enforce that if :math:`u^{MAC}` on an outflow face is inflowing,
   the normal velocity component must be outflowing or zero. We do this by imposing

   * on the low side, if :math:`u^{MAC}\ge 0` (i.e the flow is
     coming in at an outflow face) and :math:`s` is the x-velocity, then
     :math:`s_L = s_R = \min(s_R,0).`

   * on the high side, if :math:`u^{MAC}<= 0` on the domain face, then
     :math:`s_L = s_R = \max(s_L,0).`

.. note::
   Boundary conditions are imposed before the upwinding described in the :ref:`schemes` section.

API documentation can be found in the Doxygen Technical Reference,
functions `SetXEdgeBCs`_, `SetYEdgeBCs`_, `SetZEdgeBCs`_ .

.. _`SetXEdgeBCs`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceHydroBC.html#ab90f8ce229a7ebbc521dc27d65f2db9a
.. _`SetYEdgeBCs`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceHydroBC.html#a6865c2cfd50cc95f9b69ded1e8ac78ab
.. _`SetZEdgeBCs`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceHydroBC.html#a19ddc5ac50e9a6b9a98bc17f3815a62e


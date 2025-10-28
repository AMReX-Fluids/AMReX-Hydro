.. role:: cpp(code)
   :language: c++

.. _bcs:


Boundary conditions
====================

AMReX-Hydro uses underlying AMReX functionality in implementing boundary conditions
(see AMReX's documentation section :ref:`amrex:sec:basics:boundary`).
Physical boundary conditions, such as
inflow, outflow, slip/no-slip walls, etc. are ultimately linked to
mathematical Dirichlet or Neumann conditions.
See ``amrex/Src/Base/AMReX_BC_TYPES.H`` for common physical and mathematical types.


.. _mixedBC:

Mixed Boundary Conditions
--------------------------

.. note::

   This has only been tested in geometries where the Neumann and Dirchlet areas are separated
   by an embedded boundary, and is not guaranteed to work for other cases.

In addition to the bc types listed in AMReX, mixed boundary conditions can be implemented.
To this end, the (EB)Godunov advection routines support position-dependent boundary conditions.
The routines accept a boundary condition MultiFab (BC MF), or Array4. The BC MF is a cell-centered
(integer) :cpp:`iMultiFab` that carries an :cpp:`amrex::BCType` in the first ghost cell
and must fully specify the BC on all faces.
If a position-dependent BCs are passed in, they take precedent and single BC per face :cpp:`BCRecs` are
ignored.


The MacProjector supports mixed boundary conditions by making use of the underlying linear solver's
Robin BC (:ref:`amrex:sec:linearsolver:bc`) option. Robin boundary conditions are formulated as
:math:`a\phi + b\frac{\partial\phi}{\partial n} = f`.
:math:`a`, :math:`b`, and :math:`f` are each a cell-centered (real) :cpp:`MultiFab` that carries
the relevant values in the first ghost cell. All three must be specified on each level with a call to
:cpp:`setLevelBC` (see :ref:`amrex:sec:linearsolver:bc` for usage).

The NodalProjector can support mixed boundary conditions through the use of an overset mask
(see :ref:`amrex:sec:linearsolver:bc`).
The overset mask specifies a Dirichlet BC with 0 (meaning no solve is needed since the solution is known) or Neumann with 1 (meaning do the solve). Note this is an integer (not bool) MultiFab, so the values must be only either 0 or 1.



Advective BC details
--------------------

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

   * Direction-dependent (``direction_dependent``): if the face velocity is inflowing, this bc
     is identical to ``extdir``; if the face velocity is outflowing, this is equivalent to ``foextrap``.

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
functions `SetExtrapVelBCsLo`_ , `SetExtrapVelBCsHi`_ ,`SetEdgeBCsLo`_ , and `SetEdgeBCsHi`_ .

.. _`SetExtrapVelBCsLo`: https://amrex-fluids.github.io/amrex-hydro/Doxygen/html/namespaceHydroBC.html#ab90f8ce229a7ebbc521dc27d65f2db9a
.. _`SetExtrapVelBCsHi`: https://amrex-fluids.github.io/amrex-hydro/Doxygen/html/namespaceHydroBC.html#ab90f8ce229a7ebbc521dc27d65f2db9a
.. _`SetEdgeBCsLo`: https://amrex-fluids.github.io/amrex-hydro/Doxygen/html/namespaceHydroBC.html#ab90f8ce229a7ebbc521dc27d65f2db9a
.. _`SetEdgeBCsHi`: https://amrex-fluids.github.io/amrex-hydro/Doxygen/html/namespaceHydroBC.html#a6865c2cfd50cc95f9b69ded1e8ac78ab

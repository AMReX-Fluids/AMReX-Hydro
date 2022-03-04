.. _bcs:


Boundary conditions (`SetXEdgeBCs`_, `SetYEdgeBCs`_, `SetZEdgeBCs`_)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _`SetXEdgeBCs`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceHydroBC.html#ab90f8ce229a7ebbc521dc27d65f2db9a
.. _`SetYEdgeBCs`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceHydroBC.html#a6865c2cfd50cc95f9b69ded1e8ac78ab
.. _`SetZEdgeBCs`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceHydroBC.html#a19ddc5ac50e9a6b9a98bc17f3815a62e

Domain boundary conditions affect the pre-MAC extrapolated velocities in three ways.

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

For the post-MAC edge state, the domain boundary conditions affect the solution as described above in
(1) and (2) for the pre-MAC step. But, we do not impose the
no-outflow-at-inflow condition quite as described in (3); instead we
impose that if, on the low side, :math:`u^{MAC}\ge 0` (i.e the flow is
coming in at an outflow face) and :math:`s` is the x-velocity, then
:math:`s_L = s_R = \min(s_R,0).` On the high side, if
:math:`u^{MAC}<= 0` on the domain face, then
:math:`s_L = s_R = \max(s_L,0).` This enforces that if :math:`u^{MAC}`
on an outflow face is inflowing, the normal velocity component must be
outflowing or zero.


.. note::
   Boundary conditions are imposed before the upwinding.

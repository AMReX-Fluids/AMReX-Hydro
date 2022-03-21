.. include:: CustomCommands.rst

BDS
===

The Bell-Dawson-Shubin (BDS) algorithm is a higher order Godunov method for scalar
conservation laws in multiple dimensions. Satisfying the maximum principal for
constant coefficient linear advection, the BDS routine provides
accurate resolution for smooth problems while avoiding undershoot and overshoot
for non-smooth profiles. Additional details and comparisons to other
schemes can be found in the references. In this implementation, BDS closely follows the Godunov approach. The
difference appears in the computation of edge states from the MAC-projected velocities.


Pre-MAC (`ExtrapVelToFaces`_)
-----------------------------

.. _`ExtrapVelToFaces`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceGodunov.html#a1c1dcedd6781260bd8322588e1290d94

Mirrors the approach taken in the Godunov routine.


Post-MAC (`ComputeEdgeState`_)
------------------------------

.. _`ComputeEdgeState`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceGodunov.html#addea54945ce554f8b4e28dabc1c74222

Once we have the MAC-projected velocities, we project all quantities to
faces as above:

.. math::
   :label: eq3

   \tilde{s}_{i+\frac{1}{2},j,k}^{L,{n+\frac{1}{2}}} & \approx s_{i,j,k}^n + \frac{dx}{2} s_x + \frac{dt}{2} s_t \\
    & = s_{i,j,k}^n + \left( \frac{dx}{2} - s^n_{i,j,k} \frac{dt}{2} \right) (s_x^{n,lim})_{i,j,k} \\
    & + \frac{dt}{2} (-(\widehat{v s_y})_{i,j,k} - (\widehat{w s_z})_{i,j,k} + f_{x,i,j,k}^n)

extrapolated from :math:`(i,j,k)`, and

.. math::
   :label: eq4

    \tilde{s}_{i+\frac{1}{2},j,k}^{R,{n+\frac{1}{2}}} & \approx s_{i+1,j,k}^n - \frac{dx}{2} s_x + \frac{dt}{2} s_t \\
    & = s_{i+1,j,k}^n - \left( \frac{dx}{2} + s^n_{i+1,j,k} \frac{dt}{2} \right)(s^{n,lim}_x)_{i+1,j,k} \\
    & + \frac{dt}{2} (-(\widehat{v s_y})_{i+1,j,k} - (\widehat{w s_z})_{i+1,j,k} + f_{x,i+1,j,k}^n)

extrapolated from :math:`(i+1,j,k).` Here, :math:`f` is the sum of external forces, discussed later.

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

   s_{i-\frac{1}{2},j,k}^{n+\frac{1}{2}} =
   \begin{cases}
   s_L, & \mathrm{if} \; u^{MAC}_{i-\frac{1}{2},j,k}\; \ge  \; \varepsilon  \; \mathrm{else} \\
   s_R, & \mathrm{if} \; u^{MAC}_{i-\frac{1}{2},j,k}\; \le  \; -\varepsilon  \; \mathrm{else} \\
   \frac{1}{2}(s_L + s_R),
   \end{cases}

Computing the Fluxes (`ComputeFluxes`_)
---------------------------------------


Constructing the update
-----------------------

|

These alogrithms are applied in the BDS namespace. For API documentation, see
`Doxygen: BDS Namespace`_.

..
.. _`Doxygen: BDS Namespace`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceBDS.html




My Notes:
=========

- Seems like the Pre-MAC stage is the same. -- QUESTION1: Confirm?
- Check on the boundary conditions: -- QUESTION2: Confirm that this is everything.
  - Using foextrap and hoextrap, depending on what the boundary cells are labeled as.
  - Set node values equal to average of ghost cells (Since ghost cells contain boundary information)
  - One cell inward - revert to a 4-point average
  - Compute values on the faces if they are allowed to be changed, otherwise fix at ghost cell value.

- Post-MAC is what is described in the paper. -- This is computing the edge states, s.
  QUESTION3: Any place to add more information without including entire paper contents. What
  is the right amount of info?

Computing the Fluxes -- Same -- using same hydroUtils

Constructing the Update -- This is also the SAME as the Godunov -- Confirmed



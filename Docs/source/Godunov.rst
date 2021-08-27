.. math::

    \newcommand{\half}{\frac{1}{2}}
    \newcommand{\nph}{{n + \half}}
    \newcommand{\nmh}{{n - \frac{1}{2}}}
    \newcommand{\iphj}{{i+\frac{1}{2},j,k}}
    \newcommand{\ijph}{{i,j+\frac{1}{2}},k}
    \newcommand{\imhj}{{i-\frac{1}{2},j,k}}
    \newcommand{\ijmh}{{i,j-\frac{1}{2}},k}
    \newcommand{\ijkmh}{{i,j,k-\frac{1}{2}}}
    \newcommand{\ijkph}{{i,j,k+\frac{1}{2}}}
    \newcommand{\grad}{\nabla}
    \newcommand{\del}{\nabla}
    \newcommand{\AN}{[(U \cdot \nabla)U]^{n+\frac{1}{2}}}
    \newcommand{\npk}{{n + \frac{p+\half}{R}}}
    \newcommand{\nak}{{n + \frac{p}{R}}}
    \newcommand{\nmk}{{n + \frac{p-\half}{R}}}
    \newcommand{\iph}{i+\half}
    \newcommand{\imh}{i-\half}
    \newcommand{\ipmh}{i\pm\half}
    \newcommand{\jph}{j+\half}
    \newcommand{\jmh}{j-\half}
    \newcommand{\jpmh}{j\pm\half}
    \newcommand{\GMAC}{C \rightarrow E}
    \newcommand{\DMAC}{E \rightarrow C}
    \newcommand{\U}{\boldsymbol{U}}
    \newcommand{\F}{\boldsymbol{F}}

Godunov
=======

All slope computations use second-order limited slopes as described in
`Slopes`_.

.. _`Slopes`: https://amrex-codes.github.io/amrex/hydro_html/Slopes.html

We define :math:`\varepsilon = 1.e-8` in `hydro_constants.H`_

.. _`hydro_constants.H`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/group__Utilities.html#ga57d5ce9bc3bca16e249c611342f3c550

Pre-MAC (`ExtrapVelToFaces`_)
-----------------------------

.. _`ExtrapVelToFaces`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceGodunov.html#a1c1dcedd6781260bd8322588e1290d94

We extrapolate the normal velocities to cell faces using a second-order Taylor series expansion
in space and time. For each face with a non-zero area fraction, we extrapolate the normal velocity
component from the centroids of the cells on either side to the face centroid, creating left (L)
and right (R) states. For face :math:`(i+1/2,j,k)` this gives

.. math::
   :label: eq1

   \tilde{u}_{i+\frac{1}{2},j,k}^{L,{n+\frac{1}{2}}} & \approx u_{i,j,k}^n + \frac{dx}{2} u_x + \frac{dt}{2} u_t \\
    & = u_{i,j,k}^n + \left( \frac{dx}{2} - u^n_{i,j,k} \frac{dt}{2} \right) (u_x^{n,lim})_{i,j,k} \\
    & + \frac{dt}{2} (-(\widehat{v u_y})_{i,j,k} - (\widehat{w u_z})_{i,j,k} + f_{x,i,j,k}^n)


extrapolated from :math:`(i,j,k)`, and

.. math::
   :label: eq2

    \tilde{u}_{i+\frac{1}{2},j,k}^{R,{n+\frac{1}{2}}} & \approx u_{i+1,j,k}^n - \frac{dx}{2} u_x + \frac{dt}{2} u_t \\
    & = u_{i+1,j,k}^n - \left( \frac{dx}{2} + u^n_{i+1,j,k} \frac{dt}{2} \right)(u^{n,lim}_x)_{i+1,j,k} \\
    & + \frac{dt}{2} (-(\widehat{v u_y})_{i+1,j,k} - (\widehat{w u_z})_{i+1,j,k} + f_{x,i+1,j,k}^n)

extrapolated from :math:`(i+1,j,k).` Here, :math:`f` is the sum of external forces, discussed later.

In evaluating these terms the first derivatives normal to the face (in this
case :math:`u_x^{n,lim}`) are evaluated using a monotonicity-limited fourth-order
slope approximation. The limiting is done on each component of the velocity at time :math:`n` individually.

The transverse derivative terms (:math:`\widehat{v u_y}` and
:math:`\widehat{w u_z}` in this case)
are evaluated by first extrapolating all velocity components
to the transverse faces from the cell centers on either side,
then choosing between these states using the upwinding procedure
defined below.  In particular, in the :math:`y` direction we define

.. math::

    \hat{\boldsymbol{U}}^F_{i,j+\frac{1}{2},k} =  \boldsymbol{U}_{i,j,k}^n +
    \left( \frac{dy}{2} - \frac{dt}{2} v_{i,j,k}^n \right)
    (\boldsymbol{U}^{n,lim}_y)_{i,j,k}  \;\;\;

.. math::

    \hat{\boldsymbol{U}}^B_{i,j+\frac{1}{2},k} =  \boldsymbol{U}_{i,j+1,k}^n -
    \left( \frac{dy}{2} + \frac{dt}{2} v_{i,j+1,k}^n \right)
    (\boldsymbol{U}^{n,lim}_y)_{i,j+1,k} \;\;\;

Values are similarly traced from :math:`(i,j,k)` and :math:`(i,j,k+1)`
to the :math:`(i,j,k+\frac{1}{2})` faces to define
:math:`\hat{\boldsymbol{U}}^D_{i,j,k+\frac{1}{2}}` and
:math:`\hat{\boldsymbol{U}}^{U}_{i,j,k+\frac{1}{2}}`, respectively.

In this upwinding procedure we first define a normal advective
velocity on the face
(suppressing the :math:`({i,j+\frac{1}{2},k})` spatial indices on front and back
states here and in the next equation):

.. math::

    \widehat{v}^{adv}_{{i,j+\frac{1}{2},k}} = \left\{\begin{array}{lll}
     \widehat{v}^F & \mbox{if $\widehat{v}^F > 0, \;\; \widehat{v}^F + \widehat{v}^B
     > 0$} \\
     0   & \mbox{if $\widehat{v}^F \leq 0, \widehat{v}^B \geq  0$ or
    $\widehat{v}^F + \widehat{v}^B = 0$ } \\
     \widehat{v}^B & \mbox{if $\widehat{v}^B < 0, \;\; \widehat{v}^F + \widehat{v}^B
     < 0$ .} \end{array} \right.


We now upwind :math:`\widehat{\boldsymbol{U}}` based on :math:`\widehat{v}_{{i,j+\frac{1}{2},k}}^{adv}`:

.. math::

    \widehat{\boldsymbol{U}}_{{i,j+\frac{1}{2},k}} = \left\{\begin{array}{lll}
     \widehat{\boldsymbol{U}}^F & \mbox{if $\widehat{v}_{{i,j+\frac{1}{2},k}}^{adv} > 0$} \\
    \frac{1}{2} (\widehat{\boldsymbol{U}}^F + \widehat{\boldsymbol{U}}^B)  & \mbox{if $\widehat{v}_{{i,j+\frac{1}{2},k}}^{adv} = 0
    $} \\
     \widehat{\boldsymbol{U}}^B &
    \mbox{if $\widehat{v}_{{i,j+\frac{1}{2},k}}^{adv} < 0$} \end{array} \right.

After constructing :math:`\widehat{\boldsymbol{U}}_{{i,j-\frac{1}{2},k}}, \widehat{\boldsymbol{U}}_{i,j,k+\frac{1}{2}}`
and :math:`\widehat{\boldsymbol{U}}_{i,j,k-\frac{1}{2}}` in a similar manner,
we use these upwind values to form the transverse derivatives in
Eqs. :eq:`eq1` and :eq:`eq2` :

.. math::

    (\widehat{v u_y})_{i,j,k} = \frac{1}{2dy} ( \widehat{v}_{{i,j+\frac{1}{2},k}}^{adv} +
   \widehat{v}_{{i,j-\frac{1}{2},k}}^{adv} ) ( \widehat{u}_{{i,j+\frac{1}{2},k}} - \widehat{u}_{{i,j-\frac{1}{2},k}} )

.. math::
    (\widehat{w u_z})_{i,j,k} = \frac{1}{2dz} (\widehat{w}_{i,j,k+\frac{1}{2}}^{adv} +
   \widehat{w}_{i,j,k-\frac{1}{2}}^{adv} ) ( \widehat{u}_{i,j,k+\frac{1}{2}} - \widehat{u}_{i,j,k-\frac{1}{2}} )

The normal velocity at each face is then determined by an upwinding procedure
based on the states predicted from the cell centers on either side.  The
procedure is similar to that described above, i.e.
(suppressing the (:math:`i+\frac{1}{2},j,k`) indices)

.. math::

    \tilde{u}^{n+\frac{1}{2}}_{{i+\frac{1}{2},j,k}} = \left\{\begin{array}{lll}
    \tilde{u}^{L,n+\frac{1}{2}}
    & \mbox{if $\tilde{u}^{L,n+\frac{1}{2}} > 0$ and $ \tilde{u}^{L,n+\frac{1}{2}} +
    \tilde{u}^{R,n+\frac{1}{2}} > 0$} \\
    0 & \mbox{if $\tilde{u}^{L,n+\frac{1}{2}} \leq 0, \tilde{u}^{R,n+\frac{1}{2}} \geq  0$ or
    $\tilde{u}^{L,n+\frac{1}{2}} + \tilde{u}^{R,n+\frac{1}{2}} = 0$ } \\
    \tilde{u}^{R,n+\frac{1}{2}}
    & \mbox{if $\tilde{u}^{R,n+\frac{1}{2}} < 0$ and $\tilde{u}^{L,n+\frac{1}{2}}
    + \tilde{u}^{R,n+\frac{1}{2}} < 0$}
    \end{array} \right.

We follow a similar
procedure to construct :math:`\tilde{v}^{n+\frac{1}{2}}_{i,j+\frac{1}{2},k}`
and :math:`\tilde{w}^{n+\frac{1}{2}}_{i,j,k+\frac{1}{2}}`. We refer to this unique value of
normal velocity on each face as :math:`\boldsymbol{U}^{MAC,*}`.

Boundary conditions (`SetXEdgeBCs`_, `SetYEdgeBCs`_, `SetZEdgeBCs`_)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

.. _`ComputeFluxes`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceHydroUtils.html#ab70f040557a658e70ba076c9d105bab7

The fluxes are computed from the edge states above by defining, e.g.

.. math::

   F_{i-\frac{1}{2},j,k}^{x,n+\frac{1}{2}} = u^{MAC}_{i-\frac{1}{2},j,k}\; s_{i-\frac{1}{2},j,k}^{n+\frac{1}{2}}

on all x-faces,

.. math::

   F_{i,j-\frac{1}{2},k}^{y,n+\frac{1}{2}} = v^{MAC}_{i,j-\frac{1}{2},k}\; s_{i,j-\frac{1}{2},k}^{n+\frac{1}{2}}

on all y-faces, and

.. math::

   F_{i,j,k-\frac{1}{2}}^{z,n+\frac{1}{2}} = w^{MAC}_{i,j,k-\frac{1}{2}}\; s_{i,j,k-\frac{1}{2}}^{n+\frac{1}{2}}

on all z-faces.

Constructing the update
-----------------------

If the variable, :math:`s` is to be updated conservatively, we construct

.. math::

   \nabla \cdot ({\bf u} s)^{n+\frac{1}{2}}
                             = & (F_{i+\frac{1}{2},j,k}^{x,n+\frac{1}{2}} -
                                  F_{i-\frac{1}{2},j,k}^{x,n+\frac{1}{2}})+ \\
                               & (F_{i,j+\frac{1}{2},k}^{y,n+\frac{1}{2}} -
                                  F_{i,j-\frac{1}{2},k}^{y,n+\frac{1}{2}})+ \\
                               & (F_{i,j,k+\frac{1}{2}}^{z,n+\frac{1}{2}} -
                                  F_{i,j,k-\frac{1}{2}}^{z,n+\frac{1}{2}})

while if :math:`s` is to be updated in convective form, we construct

.. math::

   ({\bf u}\cdot \nabla s)^{n+\frac{1}{2}} = \nabla \cdot ({\bf u} s)^{n+\frac{1}{2}} - s_{i,j,k}^{n+\frac{1}{2}} \; (DU)^{MAC}

where

.. math::

   (DU)^{MAC} = \; & (u^{MAC}_{i+\frac{1}{2},j,k} - u^{MAC}_{i-\frac{1}{2},j,k}) + (v^{MAC}_{i,j-\frac{1}{2},k} - v^{MAC}_{i,j-\frac{1}{2},k}) + \\
                   & (w^{MAC}_{i,j,k-\frac{1}{2}} - w^{MAC}_{i,j,k-\frac{1}{2}})

and

.. math::

   s_{i,j,k}^{{n+\frac{1}{2}}} = (1/6) (
                    s_{i-\frac{1}{2},j,k}^{{n+\frac{1}{2}}} + s_{i+\frac{1}{2},j,k}^{{n+\frac{1}{2}}}
                +   s_{i,j-\frac{1}{2},k}^{{n+\frac{1}{2}}} + s_{i,j-\frac{1}{2},k}^{{n+\frac{1}{2}}}
                +   s_{i,j,k-\frac{1}{2}}^{{n+\frac{1}{2}}} + s_{i,j,k-\frac{1}{2}}^{{n+\frac{1}{2}}} )

|

These alogrithms are applied in the Godunov namespace. For API documentation, see
`Doxygen: Godunov Namespace`_.

.. _`Doxygen: Godunov Namespace`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceGodunov.html

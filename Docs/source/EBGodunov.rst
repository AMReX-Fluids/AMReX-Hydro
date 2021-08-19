EBGodunov
=========

All slope computations use fourth-order limited slopes as described in `Slopes`_ for cells for which
this calculation would not use any information in cut or covered cells; otherwise the slope computations
use a least squares approach also described in `Slopes`_ .

.. _`Slopes`: https://amrex-codes.github.io/amrex/hydro_html/Slopes.html

We define :math:`\varepsilon = 1.e-8` in **Utils / hydro_constants.H**

Notation
--------

For every cut cell we define :math:`a_x`, :math:`a_y,` and :math:`a_z` to be the area fractions of the faces
and :math:`V` is the volume fraction of the cell.  All area and volume fractions are greater than or equal to zero
and less than or equal to 1.

Pre-MAC (`ExtrapVelToFaces`_)
----------------------------

.. _`ExtrapVelToFaces`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceEBGodunov.html#abea06da38cd7e2c6a6ed94d761c4e996

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
slope approximation for cells where the stencil would not include any cut or covered cells.
The limiting is done on each component of the velocity at time :math:`n` individually. 

When one or more cells on either side is a cut cell, we instead use a least squares fit centered on :math:`(i,j,k)` that uses
all regular and cut-cell neighbors, compute slopes in all three coordinate directions. 
We then define the left and right states by extrapolating from the cell centroid to the 
face centroid using slopes in all three coordinate directions as necessary.

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

Post-MAC (`ComputeEdgestate`_)
------------------------------

.. _`ComputeEdgeState`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceEBGodunov.html#afb5b3b4bcea09a8aeeb568ddde3a46e4

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

   s_{i-\frac{1}{2},j,k}^{{n+\frac{1}{2}}} = 
   \begin{cases}
   s_L, & \mathrm{if} \; u^{MAC}_{i-\frac{1}{2},j,k}\; \ge  \; \varepsilon  \; \mathrm{else} \\
   s_R, & \mathrm{if} \; u^{MAC}_{i-\frac{1}{2},j,k}\; \le  \; -\varepsilon  \; \mathrm{else} \\
   \frac{1}{2}(s_L + s_R), 
   \end{cases}

Constructing the update
-----------------------

If the variable, :math:`s` is to be updated conservatively, on all cells with :math:`V_{i,j,k} > 0` we construct

.. math::

   \nabla \cdot ({\bf u}s) = ( 
                           & (a_{i+\frac{1}{2},j,k} \; u^{MAC}_{i+\frac{1}{2},j,k}\; s_{i+\frac{1}{2},j,k}^{{n+\frac{1}{2}}} 
                             - a_{i-\frac{1}{2},j,k} \; u^{MAC}_{i-\frac{1}{2},j,k}\; s_{i-\frac{1}{2},j,k}^{{n+\frac{1}{2}}}) + \\
                           & (a_{i,j+\frac{1}{2},k} \; v^{MAC}_{i,j-\frac{1}{2},k}\; s_{i,j+\frac{1}{2},k}^{{n+\frac{1}{2}}} 
                            - a_{i,j-\frac{1}{2},k} \; v^{MAC}_{i,j-\frac{1}{2},k}\; s_{i,j-\frac{1}{2},k}^{{n+\frac{1}{2}}}) + \\
                           & (a_{i,j,k+\frac{1}{2}} \; w^{MAC}_{i,j,k-\frac{1}{2}}\; s_{i,j,k+\frac{1}{2}}^{{n+\frac{1}{2}}} 
                            - a_{i,j,k-\frac{1}{2}} \; w^{MAC}_{i,j,k-\frac{1}{2}}\; s_{i,j,k-\frac{1}{2}}^{{n+\frac{1}{2}}}) ) / V_{i,j,k}

while if :math:`s` is to be updated in convective form, we construct

.. math::

   ({\bf u}\cdot \nabla s) = \nabla \cdot ({\bf u}s) - s_{i,j,k}^{{n+\frac{1}{2}}} (DU)^{MAC}

where

.. math::

   (DU)^{MAC}  = ( & (a_{i+\frac{1}{2},j,k} u^{MAC}_{i+\frac{1}{2},j,k}- a_{i-\frac{1}{2},j,k} u^{MAC}_{i-\frac{1}{2},j,k}) + \\
                   & (a_{i,j+\frac{1}{2},k} v^{MAC}_{i,j-\frac{1}{2},k}- a_{i,j-\frac{1}{2},k} v^{MAC}_{i,j-\frac{1}{2},k}) + \\
                   & (a_{i,j,k+\frac{1}{2}} w^{MAC}_{i,j,k-\frac{1}{2}}- a_{i,j,k-\frac{1}{2}} w^{MAC}_{i,j,k-\frac{1}{2}}) ) / V_{i,j,k} 

and

.. math::

   s_{i,j,k}^{{n+\frac{1}{2}}} = (1/6) ( 
                    s_{i-\frac{1}{2},j,k}^{{n+\frac{1}{2}}} + s_{i+\frac{1}{2},j,k}^{{n+\frac{1}{2}}}
                +   s_{i,j-\frac{1}{2},k}^{{n+\frac{1}{2}}} + s_{i,j-\frac{1}{2},k}^{{n+\frac{1}{2}}}
                +   s_{i,j,k-\frac{1}{2}}^{{n+\frac{1}{2}}} + s_{i,j,k-\frac{1}{2}}^{{n+\frac{1}{2}}} )
|
|
|

These alogrithms are applied in the EBGodunov namespace. For API documentation, see 
`Doxygen: EBGodunov Namespace`_.

.. _`Doxygen: EBGodunov Namespace`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceEBGodunov.html

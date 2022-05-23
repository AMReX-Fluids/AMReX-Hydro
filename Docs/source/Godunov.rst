.. include:: CustomCommands.rst

.. _godunov:

Godunov Methods
---------------

AMReX-Hydro's implementation uses dimenensionally unsplit algorithms with full corner coupling in 3D,
with the option to use either piecewise linear (PLM) :cite:`colella:1990, saltzman`
or piecewise parabolic (PPM) :cite:`ppm, millercolella:2002` reconstructions of the state.

These alogrithms are applied in the Godunov namespace. For API documentation, see
`Doxygen: Godunov Namespace`_.

.. _`Doxygen: Godunov Namespace`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceGodunov.html


.. _godunov-pre-mac:

Pre-MAC (API ref. `Godunov::ExtrapVelToFaces <https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceGodunov.html#a1c1dcedd6781260bd8322588e1290d94>`_)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We extrapolate the normal velocities to cell faces using a second-order Taylor series expansion
in space and time. For each face, we extrapolate the normal velocity
component from the centers of the cells on either side to the face, creating left (L)
and right (R) states. For face :math:`(i+1/2,j,k)` this gives

.. math::
   :label: eq1

   \tilde{u}_{i+\frac{1}{2},j,k}^{L,{n+\frac{1}{2}}} \approx & u_{i,j,k}^n + \frac{dx}{2} u_x + \frac{dt}{2} u_t \\
    = & u_{i,j,k}^n + \left( \frac{dx}{2} - u^n_{i,j,k} \frac{dt}{2} \right) (u_x^{n,lim})_{i,j,k} \\
    & + \frac{dt}{2} (-(\widehat{v u_y})_{i,j,k} - (\widehat{w u_z})_{i,j,k} + F_{x,i,j,k}^n)


extrapolated from :math:`(i,j,k)`, and

.. math::
   :label: eq2

    \tilde{u}_{i+\frac{1}{2},j,k}^{R,{n+\frac{1}{2}}} \approx & u_{i+1,j,k}^n - \frac{dx}{2} u_x + \frac{dt}{2} u_t \\
    = & u_{i+1,j,k}^n - \left( \frac{dx}{2} + u^n_{i+1,j,k} \frac{dt}{2} \right)(u^{n,lim}_x)_{i+1,j,k} \\
    & + \frac{dt}{2} (-(\widehat{v u_y})_{i+1,j,k} - (\widehat{w u_z})_{i+1,j,k} + F_{x,i+1,j,k}^n)

extrapolated from :math:`(i+1,j,k).` Here, :math:`\F` is the sum of forcing terms, which typically
might include viscous, gravitational, and the pressure gradient terms.

.. Given the above equation, the viscous terms and pressure gradient have been lumped into F ...

If the parameter ``use_ppm`` is false, the first derivatives normal to the face (in this
case :math:`u_x^{n,lim}`) are evaluated using a monotonicity-limited fourth-order
slope approximation as described in :ref:`slopes`.
Otherwise the Piecewise Parabolic Method  (PPM) :cite:`ppm, millercolella:2002` is used
to compute

.. math::

    \hat{u}_{i+\frac{1}{2},j,k}^{L} = & u_{i,j,k}^n + \left( \frac{dx}{2} - u^n_{i,j,k} \frac{dt}{2} \right) (u_x^{n,lim})_{i,j,k} \\
    \hat{u}_{i+\frac{1}{2},j,k}^{R} = & u_{i+1,j,k}^n - \left( \frac{dx}{2} + u^n_{i+1,j,k} \frac{dt}{2} \right)(u^{n,lim}_x)_{i+1,j,k}

The transverse derivative terms (:math:`\widehat{v u_y}` and :math:`\widehat{w u_z}` in this case)
are evaluated using a three step process: first extrapolating all velocity components
to the transverse faces from the cell centers on either side (either via PLM or PPM),
applying boundary conditions
and then choosing between these states using the upwinding procedure
defined below.  In particular, in the :math:`y` direction we define

.. States on y-faces only have y BCs enforced, etc. This means y-faces in the x boundary don't have any additional BCs enforced (beyond the fact that the cell-centered U used to compute them had BCs enforced and slopes takes certain BCs into account). However, the hi/lo edge states have BCs enforced before the final upwind. No backflow prevention for trans terms, beacuse they're transverse.

.. math::

    \widehat{\boldsymbol{U}}^F_{i,j+\frac{1}{2},k} = & \boldsymbol{U}_{i,j,k}^n +
    \left( \frac{dy}{2} - \frac{dt}{2} v_{i,j,k}^n \right)
    (\boldsymbol{U}^{n,lim}_y)_{i,j,k}  \;\;\;

    \widehat{\boldsymbol{U}}^B_{i,j+\frac{1}{2},k} = & \boldsymbol{U}_{i,j+1,k}^n -
    \left( \frac{dy}{2} + \frac{dt}{2} v_{i,j+1,k}^n \right)
    (\boldsymbol{U}^{n,lim}_y)_{i,j+1,k} \;\;\;

Values are similarly traced from :math:`(i,j,k)` and :math:`(i,j,k+1)`
to the :math:`(i,j,k+\frac{1}{2})` faces to define
:math:`\widehat{\boldsymbol{U}}^D_{i,j,k+\frac{1}{2}}` and
:math:`\widehat{\boldsymbol{U}}^{U}_{i,j,k+\frac{1}{2}}`, respectively.

If the parameter ``use_forces_in_trans`` is true, the forcing terms (:math:`\F` in
Eqs. :eq:`eq1` and :eq:`eq2`) are added to :math:`\widehat{\U}` now. Otherwise, they are included later in the the
formation of the edge state.

.. FIXME? It appears that for use_forces_in_trans, they're added both before upwinding to define the advective velocity AND before upwinding to define the transverse velocities. That doesn't quite seem right... Ooops, I was making a local var into a pointer. What the code does is create a local var (twice) and then add the forcing term only to the local var. So the effect is the results have the correct factor of the forcing terms in there (not 2x as I was worried about).

Next, boundary conditions are enforced on domain faces as decribed in :ref:`bcs` #2.
Note that this means face-based values lying within the physical boundary but not exactly on the
boundary face (e.g. values located on the y-faces of ghost cells abutting the x-boundary but not on the y-boundary)
do not have boundary conditions enforced at this point.

Now, we define a normal advective velocity on the face
(suppressing the :math:`({i,j+\frac{1}{2},k})` spatial indices on front and back
states here and in the next equation):

.. math::

    \widehat{v}^{adv}_{{i,j+\frac{1}{2},k}} = \left\{\begin{array}{lll}
     \widehat{v}^F & \mbox{if $\widehat{v}^F > 0, \;\; \widehat{v}^F + \widehat{v}^B
     > \varepsilon $} \\
     0   & \mbox{if $\widehat{v}^F \leq 0, \widehat{v}^B \geq  0$ or
    $ \lvert \widehat{v}^F + \widehat{v}^B \rvert < \varepsilon $ } \\
     \widehat{v}^B & \mbox{if $\widehat{v}^B < 0, \;\; \widehat{v}^F + \widehat{v}^B
     < \varepsilon $ .} \end{array} \right.


We now upwind :math:`\widehat{\boldsymbol{U}}` based on :math:`\widehat{v}_{{i,j+\frac{1}{2},k}}^{adv}`:

.. math::

    \widehat{\boldsymbol{U}}_{{i,j+\frac{1}{2},k}} = \left\{\begin{array}{lll}
     \widehat{\boldsymbol{U}}^F & \mbox{if $\widehat{v}_{{i,j+\frac{1}{2},k}}^{adv} > \varepsilon $} \\
     \widehat{\boldsymbol{U}}^B & \mbox{if $\widehat{v}_{{i,j+\frac{1}{2},k}}^{adv} < - \varepsilon $} \\
     \frac{1}{2} (\widehat{\boldsymbol{U}}^F + \widehat{\boldsymbol{U}}^B)  & \mbox{otherwise}
    \end{array} \right.

.. This is where the corner coupling happens in 3d, then there's another BC and upwind in the code that seems to be missing here... Why does the code break up the derivative in the way it does though?

In 3D, we complete the intermediate transverse-face centered states by accounting for transverse corner coupling
following :cite:`colella:1990, saltzman`. For example, for :math:`u` predicted to y-faces we modify
with a factor of the z-derivative at cell centers

.. math::

   \breve{u}_{i,j+\half,k}^{F} = & \hat{u}_{i,j+\frac{1}{2},k}^{F} - \frac{dt}{3} \left( \hat{w}^{adv} u_z \right)_{i,j,k} \\
             = & \hat{u}_{i,j+\frac{1}{2},k}^{F}
           - \frac{dt}{3} \left( \frac{\hat{u}_{i,j,k+\half} \hat{w}^{adv}_{i,j,k+\half} - \hat{u}_{i,j,k-\half} \hat{w}^{adv}_{i,j,k-\half} }{dz} \right) \\
               & + \frac{dt}{3} \left( \frac{\hat{w}^{adv}_{i,j,k+\half} - \hat{w}^{adv}_{i,j,k-\half}}{dz} \right) u_{i,j,k} \\
   \breve{u}_{i,j+\half,k}^{B} = & \hat{u}_{i,j+\frac{1}{2},k}^{B} - \frac{dt}{3} \left( \hat{w}^{adv} u_z \right)_{i,j+1,k} \\
             = & \hat{u}_{i,j+\frac{1}{2},k}^{B}
           - \frac{dt}{3} \left( \frac{\hat{u}_{i,j+1,k+\half} \hat{w}^{adv}_{i,j+1,k+\half} - \hat{u}_{i,j+1,k-\half} \hat{w}^{adv}_{i,j+1,k-\half}}{dz} \right) \\
               & + \frac{dt}{3} \left( \frac{\hat{w}^{adv}_{i,j+1,k+\half} - \hat{w}^{adv}_{i,j+1,k-\half}}{dz} \right) u_{i,j+1,k} \\

and then apply boundary conditions on domian faces before upwinding according to
:math:`\hat{v}_{i,j+\frac{1}{2},k}^{adv}` as was done above for :math:`\widehat{\U}`.
:math:`\widebreve{\boldsymbol{U}}_{{i,j-\frac{1}{2},k}}, \widebreve{\boldsymbol{U}}_{i,j,k+\frac{1}{2}}`
and :math:`\widebreve{\boldsymbol{U}}_{i,j,k-\frac{1}{2}}` are constructed in a similar manner.
For 2D, we take :math:`\widebreve{\U} = \widehat{\U}`.

We use these upwind values to form the transverse derivatives in Eqs. :eq:`eq1` and :eq:`eq2` :

.. math::
   :label: trans

    (\widehat{v u_y})_{i,j,k} = \frac{1}{2dy} ( \hat{v}_{{i,j+\frac{1}{2},k}}^{adv} +
   \hat{v}_{{i,j-\frac{1}{2},k}}^{adv} ) ( \breve{u}_{{i,j+\frac{1}{2},k}} - \breve{u}_{{i,j-\frac{1}{2},k}} )

.. math::
    (\widehat{w u_z})_{i,j,k} = \frac{1}{2dz} (\hat{w}_{i,j,k+\frac{1}{2}}^{adv} +
   \hat{w}_{i,j,k-\frac{1}{2}}^{adv} ) ( \breve{u}_{i,j,k+\frac{1}{2}} - \breve{u}_{i,j,k-\frac{1}{2}} )

We now have all the terms needed to form :math:`\tilde{u}`.
If ``use_forces_in_trans`` is false, the forcing terms were not included in the computation of the
transverse deriviates and are instead included at this point.
We apply boundary conditions on domain faces,
including preventing backflow (as decribed in :ref:`bcs` #2 & 3).

The normal velocity at each face is then determined by an upwinding procedure
based on the states predicted from the cell centers on either side.  The
procedure is similar to that described above, i.e.
(suppressing the (:math:`i+\frac{1}{2},j,k`) indices)


.. math::

    \tilde{u}^{n+\frac{1}{2}}_{{i+\frac{1}{2},j,k}} = \left\{\begin{array}{lll}
    \tilde{u}^{L,n+\frac{1}{2}}
    & \mbox{if $\tilde{u}^{L,n+\frac{1}{2}} > 0$ and $ \tilde{u}^{L,n+\frac{1}{2}} +
    \tilde{u}^{R,n+\frac{1}{2}} > \varepsilon $} \\
    0 & \mbox{if $\tilde{u}^{L,n+\frac{1}{2}} \leq 0, \tilde{u}^{R,n+\frac{1}{2}} \geq  0$ or
    $ | \tilde{u}^{L,n+\frac{1}{2}} + \tilde{u}^{R,n+\frac{1}{2}} | < \varepsilon $ } \\
    \tilde{u}^{R,n+\frac{1}{2}}
    & \mbox{if $\tilde{u}^{R,n+\frac{1}{2}} < 0$ and $\tilde{u}^{L,n+\frac{1}{2}}
    + \tilde{u}^{R,n+\frac{1}{2}} < \varepsilon $}
    \end{array} \right.

We follow a similar
procedure to construct :math:`\tilde{v}^{n+\frac{1}{2}}_{i,j+\frac{1}{2},k}`
and :math:`\tilde{w}^{n+\frac{1}{2}}_{i,j,k+\frac{1}{2}}`. We refer to this unique value of
normal velocity on each face as :math:`\boldsymbol{U}^{MAC,*}`.


Post-MAC (API ref. `Godnuov::ComputeEdgeState <https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceGodunov.html#addea54945ce554f8b4e28dabc1c74222>`_)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Here, the face-centered advective velocity field, which we will call :math:`\U^{MAC}`, is already known.
Typically, this is the MAC projection of :math:`\U^{MAC,*}`, which would have been computed via the Pre-MAC
procedure detailed above.
We extrapolate all quantities to faces as above:

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

extrapolated from :math:`(i+1,j,k).` Details of how these terms are computed are analogous to the Pre-MAC case.

At each face we then upwind based on :math:`u^{MAC}_{i-\frac{1}{2},j,k}`

.. math::

   \tilde{s}_{i-\frac{1}{2},j,k}^{n+\frac{1}{2}} =
   \begin{cases}
   \tilde{s}_L, & \mathrm{if} \; u^{MAC}_{i-\frac{1}{2},j,k}\; \ge  \; \varepsilon  \; \mathrm{else} \\
   \tilde{s}_R, & \mathrm{if} \; u^{MAC}_{i-\frac{1}{2},j,k}\; \le  \; -\varepsilon  \; \mathrm{else} \\
   \frac{1}{2}(\tilde{s}_L + \tilde{s}_R),
   \end{cases}

.. NOTE the pieces are put together in a different way for the transverse terms in order to match up with EB. Need to look at this in more detail...

.. _ebgodunov:

Godunov with Embedded Boundaries (EBGodunov)
---------------------------------------------

AMReX-Hydro contains an embedded boundary (EB) aware version of the Godunov algorithm
discussed above, although with fewer options available.
This EB implementation employs a piecewise linear method (PLM) :cite:`needref`, and
always includes any forcing terms *after* the computation of the transverse terms.
EBGodunov attempts to use fourth-order limited slopes wherever possible, as described in :ref:`EBslopes`.


.. _pre-mac:

Pre-MAC (API ref. `EBGodunov::ExtrapVelToFaces <https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceEBGodunov.html#abea06da38cd7e2c6a6ed94d761c4e996>`_)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We extrapolate the normal velocities to cell faces using a second-order Taylor series expansion
in space and time. For each face with a non-zero area fraction, we extrapolate the normal velocity
component from the centroids of the cells on either side to the face centroid, creating left (L)
and right (R) states. For face :math:`(i+1/2,j,k)` this gives

.. math::
   :label: eq1-ebg

   \tilde{u}_{i+\frac{1}{2},j,k}^{L,\frac{1}{2}} = \hat{u}_{i+\frac{1}{2},j,k}^{L} +
   \frac{dt}{2} \; (-(\widehat{v u_y})_{i,j,k} - (\widehat{w u_z})_{i,j,k} + F_{x,i,j,k}^n)

extrapolated from :math:`(i,j,k)`, where

.. math::
   :label: eq1-ebg2

   \hat{u}_{i+\frac{1}{2},j,k}^{L} = u_{i,j,k}^n +
   \left( \delta_x - \frac{dt}{2} u_{i,j,k}^n \right)
   \; {u^x}_{i,j,k} +  \delta_y \; {u^y}_{i,j,k} + \delta_z \; {u^z}_{i,j,k}

and

.. math::
   :label: eq2-ebg

   \tilde{u}_{i+\frac{1}{2},j,k}^{R,\frac{1}{2}} = \hat{u}_{i+\frac{1}{2},j,k}^{R} +
   \frac{dt}{2} (-(\widehat{v u_y})_{i+1,j,k} - (\widehat{w u_z})_{i+1,j,k} + F_{x,i+1,j,k}^n)

extrapolated from :math:`(i+1,j,k),` where

.. math::
   :label: eq2-ebg2

   \hat{u}_{i+\frac{1}{2},j,k}^{R} = u_{i+1,j,k}^n +
   \left(\delta_x  - \frac{dt}{2} u_{i,j,k}^n \right)
   \; {u^x}_{i+1,j,k} +  \delta_y \; {u^y}_{i+1,j,k} + \delta_z \; {u^z}_{i+1,j,k}

Here, :math:`F` is the sum of forcing terms;
:math:`\delta_x,` :math:`\delta_y` and :math:`\delta_z` are the components of the distance vector
from the cell centroid to the face centroid of the :math:`x`-face at :math:`(i-\half,j,k)`;
and the slopes :math:`(u^x,u^y,u^z)` are calculated as described in the :ref:`EBslopes` section.

The transverse derivative terms ( :math:`\widehat{v u_y}` and :math:`\widehat{w u_z}` in this case)
are evaluated by first extrapolating all velocity components
to the face centroids of the transverse faces from the cell centers on either side,
applying boundary conditions on domain faces, and then choosing between these states using the upwinding procedure
defined below.  In particular, in the :math:`y` direction we define
:math:`\widehat{\boldsymbol{U}}^F_{i,j+\frac{1}{2},k}` and
:math:`\widehat{\boldsymbol{U}}^B_{i,j+\frac{1}{2},k}`
analogously to how we defined
:math:`\hat{u}_{i+\frac{1}{2},j,k}^{R}` and :math:`\hat{u}_{i+\frac{1}{2},j,k}^{L}`,
but here on the y-faces and including all three velocity components.
Values are similarly traced from :math:`(i,j,k)` and :math:`(i,j,k+1)`
to the :math:`(i,j,k+\frac{1}{2})` faces to define
:math:`\widehat{\boldsymbol{U}}^D_{i,j,k+\frac{1}{2}}` and
:math:`\widehat{\boldsymbol{U}}^{U}_{i,j,k+\frac{1}{2}}`, respectively.

In this upwinding procedure we first define a normal advective
velocity on the face
(suppressing the :math:`({i,j+\frac{1}{2},k})` spatial indices on front and back
states here and in the next equation):

.. math::
    \widehat{v}^{adv}_{{i,j+\frac{1}{2},k}} = \left\{\begin{array}{lll}
     \widehat{v}^F & \mbox{if $\widehat{v}^F > 0, \;\; \widehat{v}^F + \widehat{v}^B
     > \varepsilon $} \\
     0   & \mbox{if $\widehat{v}^F \leq 0, \widehat{v}^B \geq  0$ or
    $ | \widehat{v}^F + \widehat{v}^B | < \varepsilon $ } \\
     \widehat{v}^B & \mbox{if $\widehat{v}^B < 0, \;\; \widehat{v}^F + \widehat{v}^B
     < \varepsilon $ .} \end{array} \right.


We now upwind :math:`\widehat{\boldsymbol{U}}` based on :math:`\widehat{v}_{{i,j+\frac{1}{2},k}}^{adv}`:

.. math::
    \widehat{\boldsymbol{U}}_{{i,j+\frac{1}{2},k}} = \left\{\begin{array}{lll}
     \widehat{\boldsymbol{U}}^F & \mbox{if $\widehat{v}_{{i,j+\frac{1}{2},k}}^{adv} > 0$} \\
    \frac{1}{2} (\widehat{\boldsymbol{U}}^F + \widehat{\boldsymbol{U}}^B)  & \mbox{if $\widehat{v}_{{i,j+\frac{1}{2},k}}^{adv} = 0
    $} \\
     \widehat{\boldsymbol{U}}^B &
    \mbox{if $\widehat{v}_{{i,j+\frac{1}{2},k}}^{adv} < 0$} \end{array} \right.

In 3D, we modify the intermediate transverse-face centered states to accounting for transverse corner coupling.
For example, for :math:`u` predicted to y-faces we add a factor of the z-derivative

.. math::

   \breve{u}_{i,j+\half,k}^{F} = & \hat{u}_{i,j+\frac{1}{2},k}^{F} -  \frac{dt}{3} \left( \hat{w}^{adv} u_z \right)_{i,j,k}  \\
             = & \hat{u}_{i,j+\frac{1}{2},k}^{F}
           - \frac{dt}{3}   \left[ \left( \frac{\alpha_{i,j,k+\half} \hat{u}_{i,j,k+\half} \hat{w}^{adv}_{i,j,k+\half} - \alpha_{i,j,k-\half} \hat{u}_{i,j,k-\half} \hat{w}^{adv}_{i,j,k-\half}}{dz\ \kappa_{i,j,k}} \right) \right. \\
                & - \left. \left( \frac{\alpha_{i,j,k+\half} \hat{w}^{adv}_{i,j,k+\half} - \alpha_{i,j,k-\half} \hat{w}^{adv}_{i,j,k-\half}}{dz\ \kappa_{i,j,k}} \right) u_{i,j,k} \right] \\
   \breve{u}_{i,j+\half,k}^{B} = & \hat{u}_{i,j+\frac{1}{2},k}^{B} - \frac{dt}{3} \left( \hat{w}^{adv} u_z \right)_{i,j+1,k} \\
             = & \hat{u}_{i,j+\frac{1}{2},k}^{B}
           - \frac{dt}{3} \left[ \left( \frac{\alpha_{i,j+1,k+\half} \hat{u}_{i,j+1,k+\half} \hat{w}^{adv}_{i,j+1,k+\half} - \alpha_{i,j+1,k-\half} \hat{u}_{i,j+1,k-\half} \hat{w}^{adv}_{i,j+1,k-\half}}{dz\ \kappa_{i,j+1,k}} \right) \right. \\
               & - \left. \left( \frac{\alpha_{i,j+1,k+\half} \hat{w}^{adv}_{i,j+1,k+\half} - \alpha_{i,j+1,k-\half} \hat{w}^{adv}_{i,j+1,k-\half}}{dz\ \kappa_{i,j+1,k}} \right) u_{i,j+1,k} \right]\\

where :math:`alpha` are cell face area fractions and :math:`kappa` is the cell volume fraction.
Then we apply boundary conditions on domain faces before upwinding according to
:math:`\hat{v}_{i,j+\frac{1}{2},k}^{adv}` as was done above for :math:`\widehat{\U}`.
:math:`\widebreve{\boldsymbol{U}}_{{i,j-\frac{1}{2},k}}, \widebreve{\boldsymbol{U}}_{i,j,k+\frac{1}{2}}`
and :math:`\widebreve{\boldsymbol{U}}_{i,j,k-\frac{1}{2}}` are constructed in a similar manner.
For 2D, we take :math:`\widebreve{\U} = \widehat{\U}`.

If any of the four faces that contribute to the transverse derivatives for a particular
cell have zero area, all of the transverse *and* forcing terms are identically set to 0.  For example,
when constructing :math:`\tilde{u}_{i+\half,j,k}^{L,\nph}`, if any of the areas
:math:`a_{i,\jph,k}, a_{i,\jmh,k}, a_{i,j,\kmh}` or :math:`a_{i,j,\kph}` are zero, then we simply define

.. math::
   :label: eq2-ebg3

   \tilde{u}_{i+\half,j,k}^{L,\nph} = \hat{u}_{i+\half,j,k}^{L}

and use this in the upwinding step given by Eq. :eq:`eb-finalUpwind`.

For cut faces (i.e. faces intersecting the EB), :math:`\widebreve{\U}` and :math:`\widehat{\U}_{ad}` are linearly extrapolated
from face centroids to face centers. Then the transverse derivatives are constructed from these face centered values
using the same formulas as for the non-EB case (Eq. :eq:`trans`).

The normal velocity at each face is then determined by an upwinding procedure
based on the states predicted from the cells on either side.  The
procedure is similar to that described above, i.e.
(suppressing the (:math:`i+\frac{1}{2},j,k`) indices)

.. math::
   :label: eb-finalUpwind

    \tilde{u}^{n+\frac{1}{2}}_{{i+\frac{1}{2},j,k}} = \left\{\begin{array}{lll}
    \tilde{u}^{L,n+\frac{1}{2}}
    & \mbox{if $\tilde{u}^{L,n+\frac{1}{2}} > 0$ and $ \tilde{u}^{L,n+\frac{1}{2}} +
    \tilde{u}^{R,n+\frac{1}{2}} > \varepsilon $} \\
    0 & \mbox{if $\tilde{u}^{L,n+\frac{1}{2}} \leq 0, \tilde{u}^{R,n+\frac{1}{2}} \geq  0$ or
    $ | \tilde{u}^{L,n+\frac{1}{2}} + \tilde{u}^{R,n+\frac{1}{2}} | < \varepsilon$ } \\
    \tilde{u}^{R,n+\frac{1}{2}}
    & \mbox{if $\tilde{u}^{R,n+\frac{1}{2}} < 0$ and $\tilde{u}^{L,n+\frac{1}{2}}
    + \tilde{u}^{R,n+\frac{1}{2}} < \varepsilon $}
    \end{array} \right.

We follow a similar
procedure to construct :math:`\tilde{v}^{n+\frac{1}{2}}_{i,j+\frac{1}{2},k}`
and :math:`\tilde{w}^{n+\frac{1}{2}}_{i,j,k+\frac{1}{2}}`. We refer to these unique values of
normal velocity on each face as :math:`\boldsymbol{U}^{MAC,*}`.


Post-MAC (API ref. `EBGondunov::ComputeEdgestate <https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceEBGodunov.html#afb5b3b4bcea09a8aeeb568ddde3a46e4>`_)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here, the face-centered advective velocity field, which we will call :math:`\U^{MAC}`, is already known.
Typically, this is the MAC projection of :math:`\U^{MAC,*}`, which would have been computed via the Pre-MAC
procedure detailed above.
Now we predict all quantities to faces. Let the
scalar :math:`s` represent any advected quantities as well as all three velocity
components.  We now extrapolate :math:`s` from cell centroids to face centroids
as described in Sec. :ref:`pre-mac`. For example, on face :math:`(i+1/2,j,k)` we define

.. math::
   :label: postebg-eq1

   \tilde{s}_{i+\half,j,k}^{L,\nph} = \hat{s}_{i+\half,j,k}^{L}
    + \frac{dt}{2} \; (-(\widehat{v s_y})_{i,j,k} - (\widehat{w s_z})_{i,j,k} + f_{x,i,j,k}^n)

extrapolated from :math:`(i,j,k)`, where

.. math::
   :label: postebg-eq2

   \hat{s}_{i+\half,j,k}^{L} = s_{i,j,k}^n +
    \left( \delta_x - \frac{dt}{2} u_{i,j,k}^n \right)
    \; {s^x}_{i,j,k} +  \delta_y \; {s^y}_{i,j,k} + \delta_z \; {s^z}_{i,j,k}

and

.. math::
    \tilde{s}_{i+\half,j,k}^{R,\nph} = \hat{s}_{i+\half,j,k}^{R}
    + \frac{dt}{2} (-(\widehat{v s_y})_{i+1,j,k} - (\widehat{w s_z})_{i+1,j,k} + f_{x,i+1,j,k}^n)

extrapolated from :math:`(i+1,j,k),` where

.. math::
   :label: postebg-eq3

   \hat{u}_{i+\half,j,k}^{R} = u_{i+1,j,k}^n +
        \left(\delta_x  - \frac{dt}{2} u_{i,j,k}^n \right)
     \; {s^x}_{i+1,j,k} +  \delta_y \; {s^y}_{i+1,j,k} + \delta_z \; {s^z}_{i+1,j,k}

Here again :math:`\delta_x,` :math:`\delta_y` and :math:`\delta_z` are the components of the distance
vector from the cell centroid to the face centroid of the :math:`x`-face at :math:`(i-\half,j,k)`,
and the slopes :math:`(u^x,u^y,u^z)` are calculated as decribded in the :ref:`EBslopes` section.

The transverse terms are computed exactly as described earlier for the Pre-MAC case, except for the upwinding process.
Where we previously used the predicted states themselves to upwind, we now use the component of the
advective velocity normal to the face in question.

We upwind :math:`\tilde{s}_{i+\half,j,k}^{L,\nph}` and :math:`\tilde{s}_{i+\half,j,k}^{L,\nph}` using the
normal component of :math:`\U^{MAC}` to define :math:`\tilde{s}_{i+\half,j,k}^{\nph}.`  Again, suppressing
the subscripts, we define

.. math::
   :label: postebg-eq5

   \tilde{s}^{\nph} = \left\{\begin{array}{lll}
     \tilde{s}^{L,\nph}              & \mbox{if $u^{MAC} > \varepsilon $}  \\
     \tilde{s}^{R,\nph}  & \mbox{if $u^{MAC} < \varepsilon $} \\
   \frac{1}{2} (\tilde{s}^{L,\nph} + \tilde{s}^{R,\nph}) & \mbox{otherwise}
   \end{array} \right.



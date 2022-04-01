.. include:: CustomCommands.rst

.. _godunov:

Godunov
--------

All slope computations use fourth-order limited slopes as described in
:ref:`slopes`.

These alogrithms are applied in the Godunov namespace. For API documentation, see
`Doxygen: Godunov Namespace`_.

.. _`Doxygen: Godunov Namespace`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceGodunov.html


Pre-MAC (`ExtrapVelToFaces`_)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _`ExtrapVelToFaces`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceGodunov.html#a1c1dcedd6781260bd8322588e1290d94

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

extrapolated from :math:`(i+1,j,k).` Here, :math:`\F` is the sum of forcing terms, discussed later.

.. ANd we actually need to dicuss the forcing later. Given the above equation, the viscous terms and pressure gradient have been lumped in here...
   
In evaluating these terms the first derivatives normal to the face (in this
case :math:`u_x^{n,lim}`) are evaluated using a monotonicity-limited fourth-order
slope approximation. The limiting is done on each component of the velocity at time :math:`n` individually.

The transverse derivative terms (:math:`\widehat{v u_y}` and
:math:`\widehat{w u_z}` in this case)
are evaluated by first extrapolating all velocity components
to the transverse faces from the cell centers on either side including
applying boundary conditions (as described in :ref:`bcs`),
and then choosing between these states using the upwinding procedure
defined below.  In particular, in the :math:`y` direction we define

.. No backflow prevention for trans terms, beacuse they're transverse. Do we need to be explicit?
   
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

.. FIXME? It appears that for use_forces_in_trans, they're added both before upwinding to define the advective velocity AND before upwinding to define the transverse velocities. That doesn't quite seem right...

In this upwinding procedure we first define a normal advective
velocity on the face
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

The parameter ``use_forces_in_trans`` determines whether the the forcing terms (:math:`\F` in 
Eqs. :eq:`eq1` and :eq:`eq2`) are included in the transverse terms now, before applying boundary conditions
and upwinding, or if the forcing terms are included later in the the formation of the edge state.
     
We now upwind :math:`\widehat{\boldsymbol{U}}` based on :math:`\widehat{v}_{{i,j+\frac{1}{2},k}}^{adv}`:

.. math::

    \widehat{\boldsymbol{U}}_{{i,j+\frac{1}{2},k}} = \left\{\begin{array}{lll}
     \widehat{\boldsymbol{U}}^F & \mbox{if $\widehat{v}_{{i,j+\frac{1}{2},k}}^{adv} > \varepsilon $} \\
     \widehat{\boldsymbol{U}}^B & \mbox{if $\widehat{v}_{{i,j+\frac{1}{2},k}}^{adv} < - \varepsilon $} \\
     \frac{1}{2} (\widehat{\boldsymbol{U}}^F + \widehat{\boldsymbol{U}}^B)  & \mbox{otherwise}  
    \end{array} \right.

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

If ``use_forces_in_trans`` is false, the forcing terms were not included in the computation of the
transverse deriviates and are instead included at this point, before applying boundary conditions.
   
The normal velocity at each face is then determined by an upwinding procedure
based on the states predicted from the cell centers on either side.  The
procedure is similar to that described above, i.e.
(suppressing the (:math:`i+\frac{1}{2},j,k`) indices)

.. Once we have the transverse terms, we can sum everything up to get predicted lo & hi values at edge. BCs are enforced and we prevent backflow (2&3). Then we upwind.
   

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


Post-MAC (`ComputeEdgeState`_)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _`ComputeEdgeState`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceGodunov.html#addea54945ce554f8b4e28dabc1c74222

Once we have the MAC-projected velocities, we extrapolate all quantities to
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

At each face we then upwind based on :math:`u^{MAC}_{i-\frac{1}{2},j,k}`

.. math::

   s_{i-\frac{1}{2},j,k}^{n+\frac{1}{2}} =
   \begin{cases}
   s_L, & \mathrm{if} \; u^{MAC}_{i-\frac{1}{2},j,k}\; \ge  \; \varepsilon  \; \mathrm{else} \\
   s_R, & \mathrm{if} \; u^{MAC}_{i-\frac{1}{2},j,k}\; \le  \; -\varepsilon  \; \mathrm{else} \\
   \frac{1}{2}(s_L + s_R),
   \end{cases}



.. _ebgodunov:

EBGodunov
---------

AMReX-Hydro has also implemented an embedded boundary (EB) aware version of the Godunov algorithm
discussed above.
EBGodunov attempts to use frourth-order limited slopes where possible, as described in :ref:`EBslopes`.


.. _pre-mac:

Pre-MAC (`ExtrapVelToFaces`_)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _`ExtrapVelToFaces`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceEBGodunov.html#abea06da38cd7e2c6a6ed94d761c4e996

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
   \left( \delta x - \frac{dt}{2} u_{i,j,k}^n \right)
   \; {u^x}_{i,j,k} +  \delta y \; {u^y}_{i,j,k} + \delta z \; {u^z}_{i,j,k}

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
   \; {u^x}_{i+1,j,k} +  \delta y \; {u^y}_{i+1,j,k} + \delta z \; {u^z}_{i+1,j,k}

Here, :math:`F` is the sum of external forces, discussed later.

Here the slopes :math:`(u^x,u^y,u^z)` are calculated using a least-squares fit to available data and
:math:`\delta_x,` :math:`\delta_y` and :math:`\delta_z` are the components of the distance vector
from the cell centroid to the face centroid of the :math:`x`-face at :math:`(i-\frac{1}{2},j,k)`.
These slopes are limited with a Barth-Jesperson type of limiter that enforces no new maxima or minima
when the state is predicted to the face centroids. (If sufficient data is available for cells
with unit volume fraction, this computation instead uses a standard second- or fourth-order
slope calculation with limiting as described in Colella (1985).)

We note that if any of the four faces that contribute to the transverse derivatives for a particular cell have zero area, all of the transverse *and* forcing terms are identically set to 0.  For example, when constructing :math:`\tilde{u}_{i+\half,j,k}^{L,\nph}`, if any of the areas :math:`a_{i,\jph,k}, a_{i,\jmh,k}, a_{i,j,\kmh}` or :math:`a_{i,j,\kph}` are zero, then we simply define

.. math::
   :label: eq2-ebg3

   \tilde{u}_{i+\half,j,k}^{L,\nph} = \hat{u}_{i+\half,j,k}^{L}

The transverse derivative terms ( :math:`\widehat{v u_y}` and :math:`\widehat{w u_z}` in this case)
are evaluated by first extrapolating all velocity components
to the face centroids of the transverse faces from the cell centers on either side,
then choosing between these states using the upwinding procedure
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

After constructing :math:`\widehat{\boldsymbol{U}}_{{i,j-\frac{1}{2},k}}, \widehat{\boldsymbol{U}}_{i,j,k+\frac{1}{2}}`
and :math:`\widehat{\boldsymbol{U}}_{i,j,k-\frac{1}{2}}` in a similar manner,
we use these upwind values to form the transverse derivatives in
Eqs. :eq:`eq1-ebg` and :eq:`eq2-ebg` :

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
    \tilde{u}^{R,n+\frac{1}{2}} > \varepsilon $} \\
    0 & \mbox{if $\tilde{u}^{L,n+\frac{1}{2}} \leq 0, \tilde{u}^{R,n+\frac{1}{2}} \geq  0$ or
    $ | \tilde{u}^{L,n+\frac{1}{2}} + \tilde{u}^{R,n+\frac{1}{2}} | < \varepsilon$ } \\
    \tilde{u}^{R,n+\frac{1}{2}}
    & \mbox{if $\tilde{u}^{R,n+\frac{1}{2}} < 0$ and $\tilde{u}^{L,n+\frac{1}{2}}
    + \tilde{u}^{R,n+\frac{1}{2}} < \varepsilon $}
    \end{array} \right.

We follow a similar
procedure to construct :math:`\tilde{v}^{n+\frac{1}{2}}_{i,j+\frac{1}{2},k}`
and :math:`\tilde{w}^{n+\frac{1}{2}}_{i,j,k+\frac{1}{2}}`. We refer to this unique value of
normal velocity on each face as :math:`\boldsymbol{U}^{MAC,*}`.


Post-MAC (`ComputeEdgestate`_)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _`ComputeEdgeState`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceEBGodunov.html#afb5b3b4bcea09a8aeeb568ddde3a46e4

Once we have the MAC-projected velocities, we project all quantities to faces. Let the scalar :math:`s` represent any advected quantities as well as all three velocity components.  We now extrapolate :math:`s` from cell centroids to face centroids as described in Sec. :ref:`pre-mac`. For example, on face :math:`(i+1/2,j,k)` we define

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

Here again the slopes :math:`(s^x,s^y,s^z)` are calculated using a least-squares fit to available data and
:math:`\delta_x,` :math:`\delta_y` and :math:`\delta_z` are the components of the distance vector from the cell centroid to the face centroid of the :math:`x`-face at :math:`(i-\half,j,k).`  The transverse terms are computed exactly as described earlier except for the upwinding process; where we previously used the predicted states themselves to upwind, we now use the component of :math:`\U^{MAC}` normal to the face in question.

We note again that if any of the four faces that contribute to the transverse derivatives for a particular cell have zero area, all of the transverse *and* forcing terms are identically set to 0.  For example, when constructing :math:`\tilde{s}_{i+\half,j,k}^{L,\nph}`, if any of the areas :math:`a_{i,\jph,k}, a_{i,\jmh,k}, a_{i,j,\kmh}` or :math:`a_{i,j,\kph}` are zero, then we simply define

.. math::
   :label: postebg-eq4

   \tilde{s}_{i+\half,j,k}^{L,\nph} = \hat{s}_{i+\half,j,k}^{L}

We upwind :math:`\tilde{s}_{i+\half,j,k}^{L,\nph}` and :math:`\tilde{s}_{i+\half,j,k}^{L,\nph}` using the normal component of :math:`\U^{MAC}` to define :math:`\tilde{s}_{i+\half,j,k}^{\nph}.`  Again, suppressing the subscripts, we define

.. math::
   :label: postebg-eq5

   \tilde{s}^{\nph} = \left\{\begin{array}{lll}
     \tilde{s}^{L,\nph}              & \mbox{if $u^{MAC} > \varepsilon $}  \\
     \tilde{s}^{R,\nph}  & \mbox{if $u^{MAC} < \varepsilon $} \\
   \frac{1}{2} (\tilde{s}^{L,\nph} + \tilde{s}^{R,\nph}) & \mbox{otherwise}  
   \end{array} \right.



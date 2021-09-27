.. include:: CustomCommands.rst

EBGodunov
=========

All slope computations use fourth-order limited slopes as described in `Slopes`_ for cells for which
this calculation would not use any information in cut or covered cells; otherwise the slope computations
use a least squares approach also described in `Slopes`_ .

.. _`Slopes`: https://amrex-codes.github.io/amrex/hydro_html/Slopes.html

We define :math:`\varepsilon = 1.e-8` in `hydro_constants.H`_

.. _`hydro_constants.H`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/group__Utilities.html#ga57d5ce9bc3bca16e249c611342f3c550

Notation
--------

For every cut cell we define :math:`a_x`, :math:`a_y,` and :math:`a_z` to be the area fractions of the faces
and :math:`V` is the volume fraction of the cell.  All area and volume fractions are greater than or equal to zero
and less than or equal to 1.

.. _pre-mac:

Pre-MAC (`ExtrapVelToFaces`_)
-----------------------------

.. _`ExtrapVelToFaces`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceEBGodunov.html#abea06da38cd7e2c6a6ed94d761c4e996

We extrapolate the normal velocities to cell faces using a second-order Taylor series expansion
in space and time. For each face with a non-zero area fraction, we extrapolate the normal velocity
component from the centroids of the cells on either side to the face centroid, creating left (L)
and right (R) states. For face :math:`(i+1/2,j,k)` this gives

.. math::
   :label: eq1-ebg

   \tilde{u}_{i+\frac{1}{2},j,k}^{L,\frac{1}{2}} = \hat{u}_{i+\frac{1}{2},j,k}^{L} +
   \frac{dt}{2} \; (-(\widehat{v u_y})_{i,j,k} - (\widehat{w u_z})_{i,j,k} + f_{x,i,j,k}^n)

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
   \frac{dt}{2} (-(\widehat{v u_y})_{i+1,j,k} - (\widehat{w u_z})_{i+1,j,k} + f_{x,i+1,j,k}^n)

extrapolated from :math:`(i+1,j,k),` where 

.. math::
   :label: eq2-ebg2

   \hat{u}_{i+\frac{1}{2},j,k}^{R} = u_{i+1,j,k}^n + 
   \left(\delta_x  - \frac{dt}{2} u_{i,j,k}^n \right) 
   \; {u^x}_{i+1,j,k} +  \delta y \; {u^y}_{i+1,j,k} + \delta z \; {u^z}_{i+1,j,k}

Here, :math:`f` is the sum of external forces, discussed later.

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
:math:`\widehat{\boldsymbol{U}}^T_{i,j+\frac{1}{2},k}`
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
     \tilde{s}^{L,\nph}              & \mbox{if $u^{MAC} > 0$}  \\
   \frac{1}{2} (\tilde{s}^{L,\nph} + \tilde{s}^{R,\nph}) & \mbox{if $u^{MAC} = 0$}  \\
     \tilde{s}^{R,\nph}  & \mbox{if $u^{MAC} < 0$} 
   \end{array} \right.

Computing the Fluxes (`ComputeFluxes`_)
---------------------------------------

.. _`ComputeFluxes`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceHydroUtils.html#ab70f040557a658e70ba076c9d105bab7

The fluxes are computed from the edge states above by defining, e.g.,

.. math::
   :label: fluxebg-eq1

   F_{i-\frac{1}{2},j,k}^{x,n+\frac{1}{2}} = a_{i-\frac{1}{2},j,k} \; u^{MAC}_{i-\frac{1}{2},j,k} \; s_{i-\frac{1}{2},j,k}^{n+\frac{1}{2}} \; \Delta_y \; \Delta_z

on all x-faces with non-zero area fraction,

.. math::
   :label: fluxebg-eq2

   F_{i,j-\frac{1}{2},k}^{y,n+\frac{1}{2}} = a_{i,j-\frac{1}{2},k} \; v^{MAC}_{i,j-\frac{1}{2},k} \; s_{i,j-\frac{1}{2},k}^{n+\frac{1}{2}} \; \Delta_x \; \Delta_z

on all y-faces with non-zero area fraction, and

.. math::
   :label: fluxebg-eq3

   F_{i,j,k-\frac{1}{2}}^{z,n+\frac{1}{2}} = a_{i,j,k-\frac{1}{2}} \; w^{MAC}_{i,j,k-\frac{1}{2}}\; s_{i,j,k-\frac{1}{2}}^{n+\frac{1}{2}} \; \Delta_x \; \Delta_y

on all z-faces with non-zero area fraction.

Here :math:`\Delta_x, \Delta_y,` and :math:`\Delta_z` are the cell sizes in the 3 directions.

Constructing the update
-----------------------

If the variable, :math:`s` is to be updated conservatively, on all cells with :math:`V_{i,j,k} > 0` we construct

.. math::

    \nabla \cdot ({\bf u}s)^{n+\frac{1}{2}} = (
                           & ( F_{i+\frac{1}{2},j,k}^{{x,n+\frac{1}{2}}} -F_{i-\frac{1}{2},j,k}^{{x,n+\frac{1}{2}}}) + \\
                           & ( F_{i,j+\frac{1}{2},k}^{{y,n+\frac{1}{2}}} -F_{i,j-\frac{1}{2},k}^{{y,n+\frac{1}{2}}}) + \\
                           & ( F_{i,j,k+\frac{1}{2}}^{{z,n+\frac{1}{2}}} -F_{i,j,k-\frac{1}{2}}^{{z,n+\frac{1}{2}}}) ) / (V_{i,j,k} \Delta_x \Delta_y \Delta_z)

while if :math:`s` is to be updated in convective form, we construct

.. math::

   ({\bf u}\cdot \nabla s)^{n+\frac{1}{2}} = \nabla \cdot ({\bf u}s)^{n+\frac{1}{2}} - s_{i,j,k}^{{n+\frac{1}{2}}} (DU)^{MAC}

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

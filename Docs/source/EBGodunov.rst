EBGodunov
~~~~~~~~~

To create the normal velocities on faces, we first extrapolate from the cell centers on each side using the
slopes as computed earlier, and upwind the face value to define  :math:`U^{pred}`.
To compute the x-velocity on the x-faces of regular (ie not cut) cells, we call

   .. code:: shell

            AMREX_CUDA_HOST_DEVICE_FOR_3D(ubx, i, j, k,
             {
                 // X-faces
                 Real upls     = ccvel_fab(i  ,j,k,0) - 0.5 * xslopes_fab(i  ,j,k,0);
                 Real umns     = ccvel_fab(i-1,j,k,0) + 0.5 * xslopes_fab(i-1,j,k,0);
                 if ( umns < 0.0 && upls > 0.0 ) {
                    umac_fab(i,j,k) = 0.0;
                 } else {
                    Real avg = 0.5 * ( upls + umns );
                    if ( std::abs(avg) <  small_vel) { umac_fab(i,j,k) = 0.0;
                    } else if (avg >= 0)             { umac_fab(i,j,k) = umns;
                    } else                           { umac_fab(i,j,k) = upls;
                    }
                 }
             });

For cut cells we test on whether the area fraction is non-zero:

   .. code:: shell

             AMREX_CUDA_HOST_DEVICE_FOR_3D(ubx, i, j, k,
             {
                 // X-faces
                 if (ax_fab(i,j,k) > 0.0)
                 {
                    Real upls     = ccvel_fab(i  ,j,k,0) - 0.5 * xslopes_fab(i  ,j,k,0);
                    Real umns     = ccvel_fab(i-1,j,k,0) + 0.5 * xslopes_fab(i-1,j,k,0);
                    if ( umns < 0.0 && upls > 0.0 ) {
                       umac_fab(i,j,k) = 0.0;
                    } else {
                       Real avg = 0.5 * ( upls + umns );
                       if ( std::abs(avg) <  small_vel) { umac_fab(i,j,k) = 0.0;
                       } else if (avg >= 0)             { umac_fab(i,j,k) = umns;
                       } else                           { umac_fab(i,j,k) = upls;
                       }
                    }
                 } else {
                       umac_fab(i,j,k) = huge_vel;
                 }
             });

We then perform a MAC projection on the face-centered velocities to enforce that they satisfy

.. math:: \nabla \cdot (U^{MAC})  = 0

We do this by solving

.. math:: \nabla \cdot \frac{1}{\rho} \nabla \phi^{MAC} = \nabla \cdot \left(U^{pred} \right)

then defining

.. math:: U^{MAC} = U^{pred} - \frac{1}{\rho} \nabla \phi^{MAC}

----------------

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
slope approximation. The limiting is done on each component of the velocity at time :math:`n` individually. When one of the cells on either side is a cut cell, we
instead use a least squares fit centered on :math:`(i,j,k)` that uses
all regular and cut-cell neighbors, compute slopes in all three coordinate directions. We then define the left and right states by extrapolating from the cell centroid to the face centroid using slopes in all three coordinate directions
as necessary.

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



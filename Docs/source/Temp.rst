Temp 
====


Temporary holding document.

Computing the Fluxes
--------------------

Now let :math:`S =\{\boldsymbol{U}_g,\rho,c\}`.
Time-centered values :math:`\tilde{S}^{n+\frac{1}{2}}` at each face
(i.e., :math:`\tilde{\rho}^{n+\frac{1}{2}}`, :math:`\tilde{c}^{n+\frac{1}{2}}`, and :math:`\boldsymbol{U}^{MAC,*}`
including the normal velocity component)
are determined by upwinding using :math:`\boldsymbol{U}^{MAC}`, as follows:


.. math::

    \tilde{S}^{L} \approx 
    & S_{i,j,k} + \frac{dx}{2} (S_x^{lim})_{i,j,k} - \frac{dt}{2} \left( u^{MAC}_{i-\frac{1}{2},j,k}(S_x^{lim})_{i,j,k} \right) \\
    & - \frac{dt}{2dx(V_{i,j,k})} (S_{i,j,k}) (au^{MAC}_x)_{i-\frac{1}{2},j,k} \\
    & - \frac{dt}{2dy(V_{i,j,k})} (aS_{x|y}v^{MAC})_{y,i,j-\frac{1}{2},k} \\
    & - \frac{dt}{2dz(V_{i,j,k})} (aS_{x|z}w^{MAC})_{z,i,j,k-\frac{1}{2}} \\
    & -\frac{dt}{2} (f_{x,i,j,k})



.. math::

    \tilde{S}^{R} \approx 
    & S_{i+1,j,k} + \frac{dx}{2} (S_x^{lim})_{i+1,j,k} - \frac{dt}{2} \left( u^{MAC}_{i+\frac{1}{2},j,k}(S_x^{lim})_{i+1,j,k} \right) \\
    & - \frac{dt}{2dx(V_{i+1,j,k})} (S_{i+1,j,k}) (au^{MAC}_x)_{i+\frac{1}{2},j,k} \\
    & - \frac{dt}{2dy(V_{i+1,j,k})} (aS_{x|y}v^{MAC})_{y,i+1,j-\frac{1}{2},k} \\
    & - \frac{dt}{2dz(V_{i+1,j,k})} (aS_{x|z}w^{MAC})_{z,i+1,j,k-\frac{1}{2}} \\
    & -\frac{dt}{2} (f_{x,i+1,j,k})


Here :math:`a` is the area fraction normal to the face of the cell, 
:math:`V` is the volume fraction, and :math:`S_{x|y}, S_{x|z}` are the 
transverse terms.


.. math::

    \tilde{S}_{i+\frac{1}{2},j,k} = \left\{\begin{array}{lll}
     \tilde{S}^L 
    & \mbox{if $u_{i+\frac{1}{2},j,k}^{MAC} > 0$} \\
    \frac{1}{2} (\tilde{S}^L + \tilde{S}^R)  
    & \mbox{if $u_{i+\frac{1}{2},j,k}^{MAC} = 0$} \\
     \tilde{S}^R
    & \mbox{if $u_{i+\frac{1}{2},j,k}^{MAC} < 0$} \end{array} \right.

We multiply :math:`\epsilon^n_g \tilde{S}^{n+\frac{1}{2}}` by :math:`\boldsymbol{U}_g^{MAC}` using the interpolated :math:`\epsilon^n_g`, 
to construct the fluxes for the momentum equation.

.. math::

  \boldsymbol{F}^{adv}_S = \epsilon_g^n \boldsymbol{U}_g^{MAC} \tilde{S}^{n+\frac{1}{2}}


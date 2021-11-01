.. This rst page deomonstrates how to set math definitions in a single
   file that can be used across multiple rst pages.

.. The file `MathDefs.rst` contains user defined math commands of the type:
      \newcommand{\half}{\frac{1}{2}}

.. To include these commands for use within the equations on the page, use
   the include directive. As below, this is, ".. include:: MathDefs.rst"
   It will be necessary to add this line for each page that uses the
   user definitions.

.. include:: MathDefs.rst

Blah Now let :math:`S =\{\U_g,\rho,c\}.`
Time-centered values $\tilde{S}^{\nph}$ at each face
(i.e.\, $\tilde{\rho}^{\nph}$, $\tilde{c}^{\nph}$, and $\U^{MAC,*}$
including the normal velocity component)
are determined by upwinding using \( \U^{MAC} \), as follows:

.. math::

    \tilde{S}^{L} \approx
    & S_{i,j,k} + \frac{dx}{2} (S_x^{lim})_{i,j,k} - \frac{dt}{2} \left( u^{MAC}_{\imhj}(S_x^{lim})_{i,j,k} \right) \\
    & - \frac{dt}{2dx(V_{i,j,k})} (S_{i,j,k}) (au^{MAC}_x)_{\imhj} \\
    & - \frac{dt}{2dy(V_{i,j,k})} (aS_{x|y}v^{MAC})_{y,i,j-\half,k} \\
    & - \frac{dt}{2dz(V_{i,j,k})} (aS_{x|z}w^{MAC})_{z,i,j,k-\half} \\
    & -\frac{dt}{2} (f_{x,i,j,k})



.. math::

    \tilde{S}^{R} \approx
    & S_{i+1,j,k} + \frac{dx}{2} (S_x^{lim})_{i+1,j,k} - \frac{dt}{2} \left( u^{MAC}_{\iphj}(S_x^{lim})_{i+1,j,k} \right) \\
    & - \frac{dt}{2dx(V_{i+1,j,k})} (S_{i+1,j,k}) (au^{MAC}_x)_{\iphj} \\
    & - \frac{dt}{2dy(V_{i+1,j,k})} (aS_{x|y}v^{MAC})_{y,i+1,j-\half,k} \\
    & - \frac{dt}{2dz(V_{i+1,j,k})} (aS_{x|z}w^{MAC})_{z,i+1,j,k-\half} \\
    & -\frac{dt}{2} (f_{x,i+1,j,k})



Here \( a \) is the areA fraction normal to the face of the cell, \( V \) is the volume fraction, and \( S_{x|y}, S_{x|z} \) are the transverse terms.



We multiply \( \epsilon^n_g \tilde{S}^{\nph} \) by \( \U_g^{MAC} \) using the interpolated \( \epsilon^n_g \), to construct the fluxes for the momentum equation.

.. math::

  \F^{adv}_S = \epsilon_g^n \U_g^{MAC} \tilde{S}^\nph

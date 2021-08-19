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

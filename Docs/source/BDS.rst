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


Pre-MAC (API ref. `ExtrapVelToFaces`_)
---------------------------------------

.. _`ExtrapVelToFaces`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceGodunov.html#a1c1dcedd6781260bd8322588e1290d94

The BDS routine follows the Godunov PLM method to extrapolate velocities to cell faces, see `ExtrapVelToFaces`_.

Post-MAC (API ref. `ComputeEdgeState`_)
----------------------------------------

.. _`ComputeEdgeState`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceBDS.html#

We advance the solution in time using a
three step procedure described below. In the notation,
:math:`s` is a scalar field of the form :math:`s=s(x,y,z,t)`
and :math:`{\bf u}=(u,v,w)` represents the :math:`u^{MAC}` projected velocity field in the Pre-MAC step.
:math:`s^n_{ijk}` represents the average value of :math:`s` over the cell with index :math:`(ijk)` at
time :math:`t^n`. At each face the normal velocity (e.g., :math:`u_{i+1/2,j,k}`) is assumed constant
over the time step.

- **Step 1**: Construct a limited piecewise trilinear (bilinear in 2D) representation of the solution in
  each grid cell of the form,

.. math::
    \begin{eqnarray}
    s_{ijk}(x,y,z) &=& s_{ijk} + s_{x,ijk}\cdot(x-x_i) + s_{y,ijk}\cdot(y-y_j) + s_{z,ijk}\cdot(z-z_k) \nonumber \\
    && + s_{xy,ijk}\cdot(x-x_i)(y-y_j) + s_{xz,ijk}\cdot(x-x_i)(z-z_k) \nonumber \\
    && + s_{yz,ijk}\cdot(y-y_j)(z-z_k) + s_{xyz,ijk}\cdot(x-x_i)(y-y_j)(z-z_k).
    \end{eqnarray}

- **Step 2**: Construct edge states :math:`s_{i+1/2,j,k}`, etc. by integrating limited
  piecewise trilinear (bilinear in 2D) profiles over the space-time region determined by the characteristic
  domain of dependence on the face.

- **Step 3**: Advance the solution in time using the conservative update equation,

.. math::
    \begin{eqnarray}
    s_{ijk}^{n+1} = s_{ijk}^n &&
    - \frac{\dt}{\Delta x}(u_{i+\half,j,k}s_{i+\half,j,k} - u_{i-\half,j,k}s_{i-\half,j,k}) \nonumber \\
    && - \frac{\dt}{\Delta y}(v_{i,j+\half,k}s_{i,j+\half,k} - v_{i,j-\half,k}s_{i,j-\half,k}) \nonumber \\
    && - \frac{\dt}{\Delta z}(w_{i,j,k+\half}s_{i,j,k+\half} - w_{i,j,k-\half}s_{i,j,k-\half}).
    \end{eqnarray}


Boundary Conditions
~~~~~~~~~~~~~~~~~~~

The BDS algorithm supports periodic and Dirichlet boundary conditions.
In addition we impose that if, on the low side, :math:`{\bf u}\ge 0` (i.e the flow is
coming in at an outflow face) and :math:`s` is the x-velocity, then
:math:`s_L = s_R = \min(s_R,0).` On the high side, if
:math:`{\bf u}<= 0` on the domain face, then
:math:`s_L = s_R = \max(s_L,0).` This enforces that if :math:`{\bf u}`
on an outflow face is inflowing, the normal velocity component must be
outflowing or zero.


Computing the Fluxes and Constructing the Update
------------------------------------------------

Fluxes are computed and the update is constructed in a manner similar to the
Godunov routine, see :ref:`Godunov<godunov_section>`.

|

Additional Documentation
------------------------

These algorithms are applied in the BDS namespace. For API documentation, see
`Doxygen: BDS Namespace`_.

.. _`Doxygen: BDS Namespace`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceBDS.html

For additional details, please make a request by submitting an issue through GitHub or refer to
the following paper:

- *A Three-Dimensional, Unsplut Godunov Method For Scalar Conservation Laws*,
  A. Nonaka, S. May, A. S. Almgren, and J. B. Bell,
  SIAM Journal of Scientific Computation, Vol. 33, No.4, pp. 2039-2062
  https://ccse.lbl.gov/Publications/nonaka/BDS_3d.pdf






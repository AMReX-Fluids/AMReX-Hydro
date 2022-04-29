.. include:: CustomCommands.rst

BDS Algorithm
-------------

The Bell-Dawson-Shubin (BDS) algorithm is a higher order Godunov method for scalar
conservation laws in multiple dimensions. Satisfying the maximum principal for
constant coefficient linear advection, the BDS routine provides
accurate resolution for smooth problems while avoiding undershoot and overshoot
for non-smooth profiles. Additional details and comparisons to other
schemes can be found in :cite:`BDS_3d` and references therein.

This implementation of BDS closely follows the Godunov approach and leverages some of the
same code.
The difference appears in the computation of edge states given face-centered velocities,
i.e. the Post-MAC computation.
Currently, periodic, Dirichlet, and outflow (extrapolation)
boundary conditions are supported. Embedded boundaries are not supported within BDS at
this time.
If additional functionality is desired, or if questions remain after reading this guide,
further help is available by submitting an issue through
`Github <https://github.com/AMReX-Codes/AMReX-Hydro/issues/new>`_

..
  These lines can be added when API docs are ready.
  These algorithms are applied in the BDS namespace. For API documentation, see
  `Doxygen: BDS Namespace`_.
   .. _`Doxygen: BDS Namespace`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceBDS.html


Pre-MAC
~~~~~~~

The BDS routine follows the Godunov PLM method to extrapolate velocities to cell faces,
see :ref:`Godunov Methods: Pre-MAC <godunov-pre-mac>`.


Post-MAC
~~~~~~~~

..
    These lines can be added back when the Doxygen for BDS.rst is ready
    (API ref. `BDS::ComputeEdgeState`_)
    .. _`BDS::ComputeEdgeState`: https://amrex-codes.github.io/amrex-hydro/Doxygen/html/namespaceBDS.html#

In the notation below,
:math:`s` is a scalar field of the form :math:`s=s(x,y,z,t)`
and :math:`{\U}=(u,v,w)` represents a known face-centered velocity field, typically the projected velocity field from the Pre-MAC step (:math:`\U^{MAC}`).
:math:`s^n_{ijk}` represents the average value of :math:`s` over the cell with index :math:`(ijk)` at
time :math:`t^n`. At each face the normal velocity (e.g., :math:`u_{i+1/2,j,k}`) is assumed constant
over the time step.

Obtaining the edge states is a two step process:

- **Step 1**: Construct a limited piecewise trilinear (bilinear in 2D) representation of the solution in
  each grid cell of the form,

.. math::
    \begin{eqnarray}
    s_{ijk}(x,y,z) &=& s_{ijk} + s_{x,ijk}\cdot(x-x_i) + s_{y,ijk}\cdot(y-y_j) + s_{z,ijk}\cdot(z-z_k) \nonumber \\
    && + s_{xy,ijk}\cdot(x-x_i)(y-y_j) + s_{xz,ijk}\cdot(x-x_i)(z-z_k) \nonumber \\
    && + s_{yz,ijk}\cdot(y-y_j)(z-z_k) + s_{xyz,ijk}\cdot(x-x_i)(y-y_j)(z-z_k).
    \end{eqnarray}

- **Step 2**: Construct edge states :math:`s_{i+1/2,j,k}`, etc. by integrating the limited
  piecewise trilinear (bilinear in 2D) profiles over the space-time region determined by the characteristic
  domain of dependence of the face.
  We enforce no inflow at an outflow face as described in the post-MAC :ref:`Boundary Conditions Section<no-inflow-at-outflow>`.




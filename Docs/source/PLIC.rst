.. include:: CustomCommands.rst

.. _plic:

Volume of fluid with PLIC
=========================

The ``PLIC`` directory provides geometric volume-of-fluid transport with a piecewise-linear
interface calculation, for any AMReX-based solver. The interface in a cell is the plane
:math:`\mathbf{m} \cdot \mathbf{x} = \alpha` in unit-cube cell coordinates, with the liquid on the
side :math:`\mathbf{m} \cdot \mathbf{x} \le \alpha`, and the plane is fixed by the volume fraction
:math:`F` of the cell through the analytic relations of Scardovelli & Zaleski (JCP 164, 2000).

Like Godunov, MOL and BDS, the library provides spatial operators and leaves the update of
:math:`F`, its ghost cells and the synchronization of the levels to the solver. It requires
Cartesian coordinates and does not support embedded boundaries; ``PLIC::CheckMesh`` aborts on
either. The reconstruction is selected by name, as ``advection_type`` is in
``HydroUtils::ComputeFluxesOnBoxFromState``.

What the solver provides
------------------------

:math:`F` needs ``PLIC::nGrowF(reconstruction)`` filled ghost cells. They should

- use ``amrex::cell_cons_interp`` at coarse/fine boundaries, which keeps :math:`F` in :math:`[0,1]`;
- be filled with the ``BCRec`` of :math:`F`, as for any other scalar, and then passed with that
  ``BCRec`` to ``PLIC::FillWallGhosts``.

``PLIC::FillWallGhosts`` reads the boundary condition of each domain face:

.. list-table::
   :header-rows: 1

   * - ``BCType``
     - Boundary
     - Ghost cells
   * - ``int_dir``
     - periodic
     - left alone
   * - ``foextrap``, ``hoextrap``
     - wall or outflow
     - overwritten: the planes of the interior interface are extended through the face
   * - ``ext_dir``
     - inflow
     - kept, as the volume fraction that flows in
   * - ``reflect_even``
     - symmetry
     - kept, as the mirror image: the interface meets the plane at 90 degrees
   * - ``direction_dependent``
     - inflow or outflow, face by face
     - kept
   * - ``reflect_odd`` or any other
     - --
     - aborts

At a wall the extension gives a free contact line: the interface meets the wall at whatever angle it
has, while extrapolated values would impose 90 degrees. At an outflow it continues the interface
smoothly. A ghost value depends on cells further along the wall than one box holds, so
``FillWallGhosts`` works on a whole level, and its result does not depend on the box layout. Near
the edges of a domain, where the stencils read the ghost cells of another face, it is only
approximate for a planar interface.

Projection
----------

The solver owns the projection, for which ``PLIC::SharpFaceDensity`` gives a face density that
balances gravity for a planar interface (Francois et al., JCP 213, 2006): with
:math:`\rho_f = \rho_g + (\rho_l - \rho_g)\, \phi_f`, where :math:`\phi_f` is the liquid fraction
of the segment joining the cell centers on either side of the face, :math:`\rho_f \mathbf{g}` is a
discrete gradient. Its inverse is the face coefficient of ``Hydro::MacProjector``.

Reconstructions
---------------

``PLIC::FitPlanes`` fits the plane of each mixed cell, from the normal of the selected
reconstruction and the volume fraction of the cell:

=============== ===================================================================================
``youngs``      Youngs' finite-difference normal.
``myc``         Mixed Youngs-Centered normal (Aulisa et al., JCP 225, 2007).
``myc_height``  The height-function normal along the dominant direction of the Mixed
                Youngs-Centered normal, which is exact for a planar interface, and the latter where
                the heights are not defined. The default.
=============== ===================================================================================

The geometric kernels in ``hydro_plic_K.H`` and ``hydro_plic_normals_K.H`` derive from Kynema-SGF
(formerly AMR-Wind), under its BSD 3-Clause License; the headers keep its copyright notice.

Adding a normal
---------------

A normal runs inside device kernels, so it is a template parameter rather than a virtual function.
Write a struct in ``hydro_plic_normals_K.H``, as described there, and give it a name in
``PLIC::dispatch_normal``, which every operator of the library selects its normal with.

Verification
------------

``Tests/PLIC`` checks the plane geometry, the normals, the plane fit, the wall ghost cells and the
face density to round-off. A run fails when its checks do not hold. See its README for the expected results.

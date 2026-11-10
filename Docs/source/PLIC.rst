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

Advection
---------

``PLIC::ComputeEdgeStateAndFluxes`` gives the fluxes and edge states of :math:`F` on the faces of a
box, with the unsplit geometric fluxes of Owkes & Desjardins (JCP 270, 2014). It is called in the
solver's own ``MFIter`` loop, as ``Godunov::ComputeEdgeState`` is, and its fluxes follow the same
convention as those of ``HydroUtils::ComputeFluxesOnBoxFromState``: a flux is the edge state times
the face velocity, times the face area if ``fluxes_are_area_weighted``, and the solver updates

.. math::

   F^{n+1} = F^n - \Delta t \, \nabla \cdot \mathbf{flux}.

The flux through a face is the liquid inside its flux volume, the region bounded by the face and its
corners traced back along the flow over :math:`\Delta t`, divided by :math:`\Delta t`. The update
then conserves the liquid volume to round-off and keeps :math:`F` in :math:`[0,1]` up to round-off,
without clipping. The edge state is the liquid fraction of the flux volume.

Each flux depends only on data around its face, so the boxes sharing a face, and periodic images,
get the same flux to the bit. With one time step on all levels, a solver

#. computes the fluxes of each box of each level;
#. averages the fine fluxes onto the coarse faces under them (``amrex::average_down_faces``);
#. updates :math:`F` on each level, averages it down, and fills its ghost cells.

A subcycling solver advances each level on its own, and refluxes with the same fluxes: with rates,
the time factor of the ``FluxRegister`` is :math:`\Delta t`, as for any other scalar.

The face velocities must be at :math:`t^{n+1/2}`, discretely divergence-free and the same on faces
shared by boxes (as after a MAC projection), with one filled ghost layer, and the face Courant number
below 1.

.. note::

   Keep the interface on the finest level, by tagging the cells for which ``PLIC::interface_band``
   holds: flux volumes do not tile across coarse/fine boundaries, which must stay in pure cells. With
   subcycling, the reflux corrections then fall in pure cells too.

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
face density to round-off. It is also a small AMR solver over the library, with a prescribed,
discretely divergence-free velocity from a vector potential, and inputs for Zalesak's disk, the
single vortex, the 3D deformation field and a contact line. A run fails when its checks do not hold.
See its README for the expected results.

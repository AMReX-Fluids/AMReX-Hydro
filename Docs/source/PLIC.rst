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
   * - ``foextrap``, ``hoextrap``, ``hoextrapcc``
     - wall or outflow
     - overwritten: the planes of the interior interface are extended through the face
   * - ``ext_dir``, ``ext_dir_cc``
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
smoothly. The planes are extended beyond the wall cells within ``PLIC::wall_nsearch`` cells, along
the face, of the interface, and beyond the others a ghost cell takes the value of its wall cell: at
a shallow angle the extended plane would cross the ghost layers further from the contact line, where
no reconstruction stencil of an interface cell reaches. A ghost value depends on cells further along
the wall than one box holds, so ``FillWallGhosts`` works on a whole level, and its result does not
depend on the box layout. Nor does it depend on the order of the faces: each reads no ghost cells of
another, and a ghost cell beyond several (at an edge or a corner) takes the mean of their values,
over those that extend a plane to it, which for a planar interface is the plane.

Projection
----------

The solver owns the projection, for which ``PLIC::SharpFaceDensity`` gives a face density that
balances gravity for a planar interface (Francois et al., JCP 213, 2006): with
:math:`\rho_f = \rho_g + (\rho_l - \rho_g)\, \phi_f`, where :math:`\phi_f` is the liquid fraction
of the segment joining the cell centers on either side of the face, :math:`\rho_f \mathbf{g}` is a
discrete gradient. Its inverse is the face coefficient of ``Hydro::MacProjector``.

Each phase density is a value or a cell-centered field (``PLIC::PhaseDensity``), for a low Mach solver
whose gas density follows its temperature and composition. A field is the density of that phase, not
the mixture density: the solver defines it in every cell, extended across the interface, with one
filled ghost layer, and a face takes the arithmetic mean of the two cells beside it. A phase density
that varies along gravity is then balanced exactly away from the interface, and to
:math:`|\nabla\rho|\, h / 8` on the faces the interface cuts. On several AMR
levels it is called on each level; the projection averages the coefficients of a fine level down
onto the coarse faces it covers.

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
then conserves the liquid volume to round-off and, as the flux volumes tile, keeps :math:`F` in
:math:`[0,1]` up to round-off, without clipping. The exception is a face whose traced face has almost
no extent normal to it: its flux volume is held near the face and falls short of the Courant number,
and its liquid is the Courant number times the liquid fraction of what remains. A cell full of one
phase stays full there, but next to the interface :math:`F` can leave :math:`[0,1]` by more than
round-off. The edge state is the liquid fraction of the flux volume, held to :math:`[0,1]`:
where shear folds the flux volume onto itself, so that its net volume nearly cancels, the flux can
differ from the edge state times the face velocity, and it is the flux that conserves :math:`F`.

Each flux depends only on data around its face, so the boxes sharing a face, and periodic images,
get the same flux to the bit. With one time step on all levels, a solver

#. computes the fluxes of each box of each level;
#. averages the fine fluxes onto the coarse faces under them (``amrex::average_down_faces``);
#. updates :math:`F` on each level, averages it down, and fills its ghost cells.

A subcycling solver advances each level on its own, and refluxes with the same fluxes: with rates,
the time factor of the ``FluxRegister`` is :math:`\Delta t`, as for any other scalar.

The face velocities must be at :math:`t^{n+1/2}`, discretely divergence-free and the same on faces
shared by boxes (as after a MAC projection), with one filled ghost layer, and the face Courant number
below 1. That is checked in every build, box by box. The check waits for the GPU, so that the boxes
of the loop do not overlap: on GPUs, check the level once with ``PLIC::CheckCourantNumber`` before
the loop and pass ``check_courant = false``. No ghost face beyond a domain face that is not periodic is read, whatever its boundary
condition (an inflow profile in those ghost faces is not seen): the faces on the domain boundary stand
in for them.

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
                the heights are not defined. The most accurate of the three.
=============== ===================================================================================

On a structure one cell thick (a film, a gap, a filament or a drop), where the neighbors of a mixed
cell on both sides along some line are pure in the same phase, the gradient of :math:`F` is null or
noise. There every reconstruction takes the principal direction of the Hessian of :math:`F`
instead, which lies across the structure.

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

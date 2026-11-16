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

:math:`F` needs ``PLIC::nGrowF(reconstruction)`` filled ghost cells. At coarse/fine boundaries, fill
them with ``amrex::cell_cons_interp``, which keeps :math:`F` in :math:`[0,1]`.

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

``Tests/PLIC`` checks the plane geometry, the normals and the plane fit to round-off. A run fails
when its checks do not hold. See its README for the expected results.

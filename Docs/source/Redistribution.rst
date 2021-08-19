Redistribution
==============

AMReX-Hydro provides support for both flux redistribution and state redistribution.

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   FluxRedistribution
   StateRedistribution


Finite Volume Discretizations
-----------------------------

Consider a system of PDEs to advance a conserved quantity :math:`U` with fluxes
:math:`F`:

.. math:: \frac{\partial U}{\partial t} + \nabla \cdot F = 0.
  :label: eqn::hypsys

A conservative, finite volume discretization starts with the divergence theorm

.. math:: \int_V \nabla \cdot F dV = \int_{\partial V} F \cdot n dA.

In an embedded boundary cell, the "conservative divergence" is discretized (as
:math:`D^c(F)`) as follows

.. math::
  :label: eqn::ebdiv

   D^c(F) = \frac{1}{\kappa h} \left( \sum^D_{d = 1}
     (F_{d, \mathrm{hi}} \, \alpha_{d, \mathrm{hi}} - F_{d, \mathrm{lo}}\, \alpha_{d, \mathrm{lo}})
     + F^{EB} \alpha^{EB} \right).

Geometry is discretely represented by volumes (:math:`V = \kappa h^d`) and
apertures (:math:`A= \alpha h^{d-1}`), where :math:`h` is the (uniform) mesh
spacing at that AMR level, :math:`\kappa` is the volume fraction and
:math:`\alpha` are the area fractions. Without multivalued cells the volume
fractions, area fractions and cell and face centroids (see
:numref:`fig::volume`) are the only geometric information needed to compute
second-order fluxes centered at the face centroids, and to infer the
connectivity of the cells. Cells are connected if adjacent on the Cartesian
mesh, and only via coordinate-aligned faces on the mesh. If an aperture,
:math:`\alpha = 0`, between two cells, they are not directly connected to each
other.

.. raw:: latex

   \begin{center}

.. |a| image:: ./EB/areas_and_volumes.png
       :width: 100%

.. |b| image:: ./EB/eb_fluxes.png
       :width: 100%

.. _fig::volume:

.. table:: Illustration of embedded boundary cutting a two-dimensional cell.
   :align: center

   +-----------------------------------------------------+------------------------------------------------------+
   |                        |a|                          |                        |b|                           |
   +-----------------------------------------------------+------------------------------------------------------+
   | | A typical two-dimensional uniform cell that is    | | Fluxes in a cut cell.                              |
   | | cut by the embedded boundary. The grey area       | |                                                    |
   | | represents the region excluded from the           | |                                                    |
   | | calculation. The portion of the cell faces        | |                                                    |
   | | faces (labelled with A) through which fluxes      | |                                                    |
   | | flow are the "uncovered" regions of the full      | |                                                    |
   | | cell faces. The volume (labelled V) is the        | |                                                    |
   | | uncovered region of the interior.                 | |                                                    |
   +-----------------------------------------------------+------------------------------------------------------+

.. raw:: latex

   \end{center}


Small Cells And Stability
-------------------------

In the context of time-explicit advance methods for, say hyperbolic
conservation laws, a naive discretization in time of :eq:`eqn::hypsys` using
:eq:`eqn::ebdiv`,

.. math:: U^{n+1} = U^{n} - \delta t D^c(F)

would have a time step constraint :math:`\delta t \sim h \kappa^{1/D}/V_m`,
which goes to zero as the size of the smallest volume fraction :math:`\kappa` in
the calculation. Since EB volume fractions can be arbitrarily small, this is an
unacceptable constraint. One way to remedy this is to create "non-conservative"
approximation to the divergence :math:`D^{nc}`, which at a cell :math:`{\bf i}`,
can be formed as an average of the conservative divergences in the neighborhood,
:math:`N_{\bf i}`, of :math:`{\bf i}`.

.. math:: D^{nc}(F)_{\bf i}= \frac{\sum_{{\bf j}\in N_{\bf i}}\kappa_{\bf j}D(F)_{\bf j}}{\sum_{{\bf j}\in N_{\bf i}}\kappa_{\bf j}}

Incorporating this form, the solution can be updated using a *hybrid
divergence*, :math:`D^H(F) = \kappa D^c(F) + (1-\kappa)D^{nc}`:

.. math:: U^{n+1,*} = U^n - \delta t D^H(F)

However, we would like our finite-volume scheme to strictly conserve the field
quantities over the domain. To enforce this, we calculate :math:`\delta M`, the
mass gained or lost by not using :math:`D^c` directly,

.. math:: \delta M_{\bf i}= \kappa (1-\kappa)(D^c(F)_{\bf i}- D^{nc}(F)_{\bf i})

This "excess material" (mass, if :math:`U=\rho`) can be *redistributed* in a
time-explicit fashion to neighboring cells, :math:`{\bf j}\in N_{\bf i}`:

.. math:: \delta M_{\bf i}= \sum_{{\bf j}\in N_{\bf i}} \delta M_{{\bf j}, {\bf i}}.

in order to preserve strict conservation over :math:`N_{\bf i}`.

Note that the physics at hand may impact the optimal choice of precisely how the
excess mass is distributed in this fashion. We introduce a weighting for
redistribution, :math:`W`,

.. math::
  :label: eqn::massweight

   \delta M_{{\bf j}, {\bf i}} =  \frac{\delta M_{\bf i}\kappa_{\bf j}
     W_{\bf j}}{\sum_{{\bf k}\in N_{\bf i}} \kappa_{\bf k}W_{\bf k}}

For all :math:`{\bf j}\in N_{\bf i}`,

.. math::

   U^{n+1}_{\bf j}= U^{n+1,*}_{\bf j}+
    \frac{\delta M_{\bf i}
     W_{\bf j}}{\sum_{{\bf k}\in N_{\bf i}} \kappa_{\bf k}W_{\bf k}}.

Typically, the redistribution neighborhood for each cell is one that can be
reached via a monotonic path in each coordinate direction of unit length (see,
e.g., :numref:`fig::redistribution`)

.. raw:: latex

   \begin{center}

.. _fig::redistribution:

.. figure:: ./EB/redist.png
   :width: 50.0%

   : Redistribution illustration. Excess mass due to using a hybrid divergence
   :math:`D^H` instead of the conservative divergence :math:`D^C` is
   distributed to neighbor cells.

.. raw:: latex

   \end{center}


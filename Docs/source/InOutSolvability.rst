.. include:: CustomCommands.rst

Enforcing inflow-outflow solvability
------------------------------------

This routine enforces solvability for inflow-outflow boundaries,
which have both inflow and outflow cells.
A flux correction factor, :math:`\alpha_\text{fcf}` is introduced
which scales the outflow velocities to match the inflow:

.. math::

    \sum_{\partial\Omega_\text{in}} {\bf u} \cdot {\bf \area} + \alpha_\text{fcf} \sum_{\partial\Omega_\text{out}} {\bf u} \cdot {\bf \area} = 0


The new flux-conserving velocities to be used for the MAC/nodal projections,
:math:`{\bf u}_\text{fc}`, are calculated from the correction factor as:

.. math::
    \alpha_\text{fcf} = \frac{-\sum_{\partial\Omega_\text{in}} {\bf u} \cdot {\bf \area}}{\sum_{\partial\Omega_\text{out}} {\bf u} \cdot {\bf \area}},

.. math::
    {\bf u}_\text{fc} = \alpha_\text{fcf} \cdot {\bf u}, \ \forall {\bf x} \in  \partial\Omega_\text{out}.

It must be noted that this routine currently only accounts for boundaries
with the math BC ``BCType::direction_dependent``, which is to be used for
an inflow-outflow boundary. It does not compute or correct the boundary velocities
over other math BC types such as those representing pure inflow or pure outflow.
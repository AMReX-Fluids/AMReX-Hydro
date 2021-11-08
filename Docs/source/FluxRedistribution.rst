

Flux Redistribution
===================

Consider a conservative update in the form:

.. math:: (\rho \phi)_t + \nabla \cdot ( \rho \phi u) = RHS

For each valid cell in the domain, compute the conservative divergence, :math:`(\nabla \cdot F)^c` ,
of the convective fluxes, :math:`F`

.. math:: (\nabla \cdot {F})^c_i = \dfrac{1}{\mathcal{V}_i} \sum_{f=1}^{N_f} ({F}_f\cdot{n}_f) A_f

Here :math:`N_f` is the number of faces of cell :math:`i`, :math:`\vec{n}_f` and :math:`A_f`
are the unit normal and area of the :math:`f` -th face respectively,
and :math:`\mathcal{V}_i` is the volume of cell :math:`i` given by

.. math:: \mathcal{V}_i = (\Delta x \Delta y \Delta z)\cdot \mathcal{K}_i

where :math:`\mathcal{K}_i` is the volume fraction of cell :math:`i` .

Now, a conservative update can be written as

.. math:: \frac{ \rho^{n+1} \phi^{n+1} - \rho^{n} \phi^{n} }{\Delta t} = - \nabla \cdot{F}^c

For each cell cut by the EB geometry, compute the non-conservative update, :math:`\nabla \cdot {F}^{nc}` ,

.. math:: \nabla\cdot{F}^{nc}_i = \dfrac{\sum\limits_{j\in N(i) } \mathcal{K}_j\nabla \cdot {F}^c_j} {\sum\limits_{j\in N(i) } {\mathcal{K}}_j}

where :math:`N(i)` is the index set of cell :math:`i` and its neighbors.

For each cell cut by the EB geometry, compute the convective update :math:`\nabla \cdot{F}^{EB}` follows:

.. math:: \nabla \cdot{F}^{EB}_i = \mathcal{K}_i\nabla \cdot{F}^{c}_i +(1-\mathcal{K}_i)\mathcal{F}^{nc}_i

For each cell cut by the EB geometry, redistribute its mass loss, :math:`\delta M_i` , to its neighbors:

.. math::  \nabla \cdot {F}^{EB}_j :=   \nabla \cdot {F}^{EB}_j + w_{ij}\delta M_i\, \qquad \forall j\in N(i)\setminus i

where the mass loss in cell :math:`i` , :math:`\delta M_i` , is given by

.. math:: \delta M_i =  \mathcal{K}_i(1- \mathcal{K}_i)[ \nabla \cdot {F}^c_i-  \nabla \cdot {F}^{nc}_i]

and the weights, :math:`w_{ij}` , are

.. math:: w_{ij} = \dfrac{1}{\sum\limits_{j\in N(i)\setminus i}  \mathcal{K}_j}

Note that :math:`\nabla \cdot{F}_i^{EB}` gives an update for :math:`\rho \phi` ; i.e.,

.. math:: \frac{(\rho \phi_i)^{n+1} - (\rho \phi_i)^{n} }{\Delta t} = - \nabla \cdot{F}^{EB}_i

Typically, the redistribution neighborhood for each cell is one that can be
reached via a monotonic path in each coordinate direction of unit length (see,
e.g., :numref:`fig::redistribution`)

.. raw:: latex

   \begin{center}

.. _fig::redistribution:

.. figure:: ./EB/redist.png
   :width: 50.0%

   : Redistribution illustration. Excess update distributed to neighbor cells.

.. raw:: latex

   \end{center}




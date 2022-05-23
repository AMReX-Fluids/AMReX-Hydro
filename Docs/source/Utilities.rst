.. include:: CustomCommands.rst

.. _utilities:

Helper functions
^^^^^^^^^^^^^^^^^


Notation
---------

- :math:`dx`, :math:`dy`, :math:`dz` : cell sizes in the x-, y- and z- directions, respectively.

- :math:`\vol_{i,j,k}` : Volume of cell-centered element :math:`(i,j,k)`.

- :math:`\area` : Area of a cell face. For example,
  :math:`\area_{i-\half ,j,k}` represents the area of the lower x-face of the :math:`(i,j,k)`-th cell,


For problems with embedded boundaries, we also define

- :math:`\kappa_{i,j,k}` : Volume fraction of cell :math:`(i,j,k)`. Data are in the range
  of :math:`[0,1]` with zero representing covered cells and one for regular
  cells, so that for regular cells :math:`\kappa_{i,j,k} \vol_{i,j,k} =  dx\ dy\ dz`

- :math:`\alpha` : Area fractions. For example,
  :math:`\alpha_{i-\half ,j,k}` corresponds to the the lower x-face of the :math:`(i,j,k)`-th cell.
  Data are in the range of :math:`[0,1]` with zero representing a covered face
  and one an un-cut face, so that for regular cells
  :math:`\alpha_{i-\half ,j,k} \area_{i-\half ,j,k} = dy\ dz`.

Note that EB is only an option for Cartesian geometries as of this writing.


.. include::   Slopes.rst

.. include::   bcs.rst

.. include::   Fluxes.rst

.. include::   AdvectiveTerm.rst


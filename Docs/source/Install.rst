


Installation Guide
==================

AMReX-Hydro is not a stand-alone code: its routines are compiled together with an application that
uses them, with either CMake or GNU Make. In both cases the AMReX-Hydro configuration (dimension,
precision, MPI, OpenMP, GPU backend and embedded boundaries) must match that of AMReX.

Building with CMake
-------------------

AMReX-Hydro defines a CMake object library target, ``amrex_hydro``. An application can add
AMReX-Hydro to its own CMake build with ``add_subdirectory``, or AMReX-Hydro can be configured and
built on its own against an installed AMReX, as follows.

Install AMReX
~~~~~~~~~~~~~

Build and install AMReX with CMake as described in `Building with CMake`_ in the AMReX
documentation. AMReX-Hydro builds its embedded boundary (EB) routines by default (``HYDRO_EB``,
see below), which requires AMReX to be built with EB support, ``-DAMReX_EB=YES``.

.. _`Building with CMake`: https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html#building-with-cmake

Clone and Configure AMReX-Hydro
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Download AMReX-Hydro from GitHub by using

   ::

      git clone https://github.com/AMReX-Fluids/AMReX-Hydro.git

2. In the ``AMReX-Hydro`` folder, create and enter a build directory,

   ::

      mkdir build
      cd build

3. Call ``cmake`` to configure AMReX-Hydro. Set ``AMReX_ROOT`` to the AMReX installation
   directory so that CMake can find ``AMReXConfig.cmake``, and set ``CMAKE_INSTALL_PREFIX`` to the
   directory to install AMReX-Hydro in. Without it, CMake installs in ``/usr/local``, which usually
   requires administrator rights. For example,

   ::

      cmake -DAMReX_ROOT=/fullpathto/amrex/installdir \
            -DCMAKE_INSTALL_PREFIX=/fullpathto/AMReX-Hydro/installdir \
            <other configuration options> ..

   .. warning::

      CMake requires an absolute path for the ``AMReX_ROOT`` variable. Therefore, a
      relative path such as ``-DAMReX_ROOT=../../amrex/installdir``
      will result in an error because CMake will be unable to find ``AMReXConfig.cmake``.

   The AMReX-Hydro configuration options are listed in the table below. CMake checks that the
   installed AMReX was built with the matching components.

Build and Install AMReX-Hydro
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

4. While in the AMReX-Hydro build directory ``AMReX-Hydro/build``, build and install AMReX-Hydro
   by entering

   ::

      make -j4 install

5. The installation directory now contains the library ``libamrex_hydro_api`` in ``lib``, the
   header files in ``include``, and the CMake package files in ``lib/cmake/AMReX-Hydro`` (``lib64``
   instead of ``lib`` on some systems). An application can then find AMReX-Hydro with
   ``find_package(AMReX-Hydro)``, setting ``AMReX-Hydro_ROOT`` to the installation directory.

CMake Options
~~~~~~~~~~~~~

+--------------------+---------------------------------------------------------------+---------------+
| Option             | Description                                                   | Default       |
+====================+===============================================================+===============+
| HYDRO_SPACEDIM     | Dimension of the AMReX build: 2 or 3                          | 3             |
+--------------------+---------------------------------------------------------------+---------------+
| HYDRO_PRECISION    | Floating-point precision of the AMReX build: DOUBLE or SINGLE | DOUBLE        |
+--------------------+---------------------------------------------------------------+---------------+
| HYDRO_MPI          | Enable MPI                                                    | YES           |
+--------------------+---------------------------------------------------------------+---------------+
| HYDRO_OMP          | Enable OpenMP                                                 | NO            |
+--------------------+---------------------------------------------------------------+---------------+
| HYDRO_GPU_BACKEND  | GPU backend: NONE, SYCL, CUDA or HIP                          | NONE          |
+--------------------+---------------------------------------------------------------+---------------+
| HYDRO_EB           | Build the EB routines (EBMOL and EBGodunov). Requires AMReX   | YES           |
|                    | built with EB support. (``HYDRO_NO_EB=YES`` is a deprecated   |               |
|                    | alias for ``HYDRO_EB=NO``.)                                   |               |
+--------------------+---------------------------------------------------------------+---------------+
| HYDRO_PROJECTIONS  | Build the MAC and nodal projections                           | YES           |
+--------------------+---------------------------------------------------------------+---------------+
| HYDRO_FFT          | Build the FFT-based MAC projection. Requires AMReX built      | NO            |
|                    | with ``AMReX_FFT=ON``                                         |               |
+--------------------+---------------------------------------------------------------+---------------+

Building with GNU Make
----------------------

Each AMReX-Hydro directory (``BDS``, ``EBGodunov``, ``EBMOL``, ``Godunov``, ``MOL``,
``Projections`` and ``Utils``) has a ``Make.package`` file. With AMReX's GNU Make system, include
the ``Make.package`` files of the directories you need and add those directories to the
include and source paths, e.g.

::

   include $(AMREX_HYDRO_HOME)/Utils/Make.package
   INCLUDE_LOCATIONS += $(AMREX_HYDRO_HOME)/Utils
   VPATH_LOCATIONS   += $(AMREX_HYDRO_HOME)/Utils

The ``EBMOL`` and ``EBGodunov`` sources are only compiled with ``USE_EB = TRUE``, and
``NO_EB = TRUE`` disables AMReX-Hydro's EB support in an EB-enabled build. The FFT-based MAC
projection is compiled with ``USE_FFT = TRUE``. The ``GNUmakefile`` files in the ``Tests``
directory are complete examples.

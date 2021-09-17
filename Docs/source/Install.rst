


Installation Guide
==================

Object files of AMReX-Hydro routines can be built with CMake. To build the objects for linking
with your projects see the steps below:


Clone, Configure and Install AMReX
-----------------------------------

1. Download AMReX from GitHub by using

  ::

    git clone https://github.com/AMReX-Codes/amrex.git

2. In the ``amrex`` folder, create and enter a build directory,

  ::

     mkdir build
     cd build

3. Call ``cmake`` to configure the AMReX installation. For compatibility, AMReX-Hydro
   requires AMReX be built with support for embedded boundaries by setting the
   flag ``-DAMReX_EB=YES`` in the AMReX configuration,

  ::

     cmake -DAMReX_EB=YES <other configuration options> ..

For additional options building AMReX please refer to `Building with CMake`_ in the AMReX
Source Documentation.

.. _`Building with CMake`: https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html#building-with-cmake


4. To install AMReX type,

   ::

     make -j4 install

   This will create a directory called ``installdir`` that contains the AMReX library
   and header files. The optional flag ``-j4`` speeds up the compilation by using 4 processes. You can adjust
   this number as you see fit, or omit it altogether.


Clone and Configure AMReX-Hydro
-------------------------------

5. Download AMReX-Hydro from GitHub by using

  ::

    git clone https://github.com/AMReX-Codes/AMReX-Hydro.git

6. In the ``AMReX-Hydro`` folder, create and enter a build directory,

  ::

     mkdir build
     cd build

7. Call ``cmake`` to configure the AMReX-Hydro installation. You will need
   to specify the location of the file ``AMReXConfig.cmake`` be setting the flag
   to the directory that contains it. In the default case this could be

  ::

     cmake -DAMReX_ROOT=/fullpathto/amrex/installdir <other configuration options> ..


.. warning::

   CMake requires an absolute path for the ``AMReX_ROOT`` variable. Therefore, a
   relative path such as

   ::

      -DAMReX_ROOT=../../amrex/installdir

   will result in error because CMake will be unable to find ``AMReXConfig.cmake``.

Make The Object Files
----------------------

8. While in the AMReX-Hydro build directory ``AMReX-Hydro/build``, create the object
   files by entering

  ::

     make -j4


9. The object files for AMReX-Hydro routines can now be found in their respective
   subdirectories at

  ::

      AMReX-Hydro/build/CMakeFiles/amrex_hydro.dir


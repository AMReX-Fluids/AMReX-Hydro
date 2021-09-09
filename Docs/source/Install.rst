


Installation Guide
==================

AMReX-Hydro can be built with CMake. To build the library for linking
with your projects see the steps below:


1. Clone the AMReX and AMReX-Hydro repos from Github.

  ::

     git clone
     git clone

2. Create and Enter a build directory.

  ::

     mkdir build
     cd build

3. Call ``cmake`` to configure the AMReX installation.

  ::

      cmake <configuration options> ../amrex
     
4. Call ``cmake`` to configure the AMReX-Hydro installation.


  :: 

     cmake ../AMReX-Hydro

5. Type ``make install`` to install as configured. 
   


6. After ``make install`` completes navigate to the directory
   ``build/Src``. There you should find the static libraries file ``libamrex.a``.
   Linking to this file when you compile your own code enables the AMReX features. 



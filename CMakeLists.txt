cmake_minimum_required(VERSION 3.14)

project( AMReX-Hydro
    # DESCRIPTION "A software framework for massively parallel, block-structured adaptive mesh refinement (AMR) applications"
#    HOMEPAGE_URL "https://amrex-codes.github.io/amrex/"
    LANGUAGES C CXX
    )

message(STATUS "CMake version: ${CMAKE_VERSION}")

#
# Check if CMAKE_BUILD_TYPE is given. If not, use default
#
if ( NOT CMAKE_BUILD_TYPE )
   set(CMAKE_CONFIGURATION_TYPES "Release;Debug;MinSizeRel;RelWithDebInfo")
   set(CMAKE_BUILD_TYPE Release
       CACHE STRING
       "Choose the build type, e.g. Release, Debug, or RelWithDebInfo." FORCE)
else ()
   message(STATUS "Build type set by user to '${CMAKE_BUILD_TYPE}'.")
endif()

#
# Options
#
set(HYDRO_SPACEDIM 3 CACHE STRING "Dimension of AMReX build: <2,3>")
option( HYDRO_EB     "Enable Embedded-Boundary support" YES)
option( HYDRO_OMP    "Enable OpenMP" NO )
option( HYDRO_MPI    "Enable MPI"   YES )


set(HYDRO_GPU_BACKEND_VALUES NONE SYCL CUDA HIP)
set(HYDRO_GPU_BACKEND NONE CACHE STRING "On-node, accelerated GPU backend: <NONE,SYCL,CUDA,HIP>")
set_property(CACHE HYDRO_GPU_BACKEND PROPERTY STRINGS ${HYDRO_GPU_BACKEND_VALUES})
if (NOT HYDRO_GPU_BACKEND IN_LIST HYDRO_GPU_BACKEND_VALUES)
   message(FATAL_ERROR "HYDRO_GPU_BACKEND=${HYDRO_GPU_BACKEND} is not allowed."
      " Must be one of ${HYDRO_GPU_BACKEND_VALUES}")
endif ()

if (HYDRO_GPU_BACKEND STREQUAL "CUDA")
   enable_language(CUDA)
endif ()

#
# Find AMReX
#
if (NOT TARGET AMReX::amrex)
   set(AMREX_REQUIRED_COMPONENTS ${HYDRO_SPACEDIM}D DOUBLE)
   if (HYDRO_MPI)
      list(APPEND AMREX_REQUIRED_COMPONENTS MPI)
   endif ()
   if (HYDRO_OMP)
      list(APPEND AMREX_REQUIRED_COMPONENTS OMP)
   endif ()
   if (NOT HYDRO_GPU_BACKEND STREQUAL "NONE")
      list(APPEND AMREX_REQUIRED_COMPONENTS ${HYDRO_GPU_BACKEND})
   endif ()
   if (HYDRO_EB)
      list(APPEND AMREX_REQUIRED_COMPONENTS EB)
   endif ()

   find_package(AMReX CONFIG REQUIRED ${AMREX_REQUIRED_COMPONENTS} )
endif ()

#
# Enable CUDA if requested
#
if (HYDRO_GPU_BACKEND STREQUAL "CUDA")
   # CMake 3.18+: CMAKE_CUDA_ARCHITECTURES
   # https://cmake.org/cmake/help/latest/policy/CMP0104.html
   if(POLICY CMP0104)
      cmake_policy(SET CMP0104 OLD)
   endif()

   include(AMReX_SetupCUDA)
   include(AMReXTargetHelpers)
endif ()


#
# Define the object library to compile
#
add_library(amrex_hydro OBJECT)
target_link_libraries(amrex_hydro PUBLIC AMReX::amrex)
add_subdirectory(Slopes)
add_subdirectory(Utils)
add_subdirectory(MOL)
add_subdirectory(Godunov)

if (HYDRO_EB)
   add_subdirectory(EBMOL)
   add_subdirectory(EBGodunov)
   add_subdirectory(Redistribution)
endif ()

if (HYDRO_GPU_BACKEND STREQUAL "CUDA")
   setup_target_for_cuda_compilation(amrex_hydro)
endif ()

if ( NOT CMAKE_CXX_FLAGS )
   target_link_libraries(amrex_hydro PUBLIC AMReX::Flags_CXX)
endif ()

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


# #
# # Enable CUDA if requested
# #
# if (AMReX_CUDA)
#     # CMake 3.18+: CMAKE_CUDA_ARCHITECTURES
#     # https://cmake.org/cmake/help/latest/policy/CMP0104.html
#     if(POLICY CMP0104)
#         cmake_policy(SET CMP0104 OLD)
#     endif()

#     enable_language(CUDA)
#     include(AMReX_SetupCUDA)
# endif ()
set(HYDRO_SPACEDIM 3 CACHE STRING "Dimension of AMReX build: <1,2,3>")



# For now just link against amrex
find_package(AMReX REQUIRED)


#
# Define the object library to compile
#
add_library(amrex_hydro OBJECT)
target_link_libraries(amrex_hydro PUBLIC AMReX::amrex)
add_subdirectory(Utils)
add_subdirectory(MOL)
add_subdirectory(Godunov)
add_subdirectory(EBMOL)
add_subdirectory(Redistribution)

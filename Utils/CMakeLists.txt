target_include_directories(
   amrex_hydro
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
   )

target_sources(
   amrex_hydro
   PRIVATE
   hydro_utils.H
   hydro_utils.cpp
   hydro_extrap_vel_to_faces.cpp
   hydro_compute_edgestate_and_flux.cpp
   hydro_utils.cpp
   hydro_constants.H
   hydro_bcs_K.H
   )

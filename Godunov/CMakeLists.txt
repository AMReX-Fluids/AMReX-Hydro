target_include_directories(
   amrex_hydro
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
   )

target_sources(
   amrex_hydro
   PRIVATE
   hydro_godunov.H
   hydro_godunov_K.H
   hydro_godunov_edge_state_${HYDRO_SPACEDIM}D.cpp
   hydro_godunov_extrap_vel_to_faces_${HYDRO_SPACEDIM}D.cpp
   hydro_godunov_plm.H
   hydro_godunov_plm.cpp
   hydro_godunov_ppm.H
   hydro_godunov_ppm.cpp
   )

if (HYDRO_SPACEDIM EQUAL 3)
   target_sources(
      amrex_hydro
      PRIVATE
      hydro_godunov_corner_couple.H
      )
endif ()

target_include_directories(
   amrex_hydro
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
   )

target_sources(
   amrex_hydro
   PRIVATE
   hydro_ebgodunov.H
   hydro_ebgodunov_extrap_vel_to_faces.cpp
   hydro_ebgodunov_extrap_vel_to_faces_${HYDRO_SPACEDIM}D.cpp
   hydro_ebgodunov_edge_state_${HYDRO_SPACEDIM}D.cpp
   hydro_ebgodunov_transverse_${HYDRO_SPACEDIM}D_K.H
   hydro_ebgodunov_edge_state_${HYDRO_SPACEDIM}D.cpp
   hydro_ebgodunov_plm.H
   hydro_ebgodunov_plm.cpp
   hydro_ebgodunov_plm_fpu.cpp
   )

if (HYDRO_SPACEDIM EQUAL 3)
   target_sources(
      amrex_hydro
      PRIVATE
      hydro_ebgodunov_corner_couple.H
      )
endif ()

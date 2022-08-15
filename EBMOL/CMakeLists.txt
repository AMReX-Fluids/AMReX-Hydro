target_include_directories(
   amrex_hydro
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
   )

target_sources(
   amrex_hydro
   PRIVATE
   hydro_ebmol.H
   hydro_ebmol_edge_state_K.H
   hydro_ebmol_edge_state.cpp
   hydro_ebmol_extrap_vel_to_faces_box.cpp
   hydro_ebmol_extrap_vel_to_faces.cpp
   )

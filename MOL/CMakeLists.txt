target_include_directories(
   amrex_hydro
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
   )

target_sources(
   amrex_hydro
   PRIVATE
   hydro_mol.H
   hydro_mol_edge_state_K.H
   hydro_mol_edge_state.cpp
   hydro_mol_extrap_vel_to_faces_box.cpp
   hydro_mol_extrap_vel_to_faces.cpp
   )

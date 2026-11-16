target_include_directories(
   amrex_hydro
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
   )

target_sources(
   amrex_hydro
   PRIVATE
   hydro_plic.H
   hydro_plic.cpp
   hydro_plic_K.H
   hydro_plic_normals_K.H
   hydro_plic_reconstruction.cpp
   )

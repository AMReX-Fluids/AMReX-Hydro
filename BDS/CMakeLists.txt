target_include_directories(
   amrex_hydro
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
   )

target_sources(
   amrex_hydro
   PRIVATE
   hydro_bds.H
   hydro_bds_edge_state_${HYDRO_SPACEDIM}D.cpp
   )


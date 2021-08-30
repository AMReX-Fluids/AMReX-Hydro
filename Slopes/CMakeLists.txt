target_include_directories(
   amrex_hydro
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
   )

target_sources(
   amrex_hydro
   PRIVATE
   hydro_slopes_K.H
   hydro_eb_slopes_${HYDRO_SPACEDIM}D_K.H
   )

target_include_directories(
   amrex_hydro
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
   )

target_sources(
   amrex_hydro
   PRIVATE
   hydro_MacProjector.cpp
   hydro_MacProjector.H
   hydro_NodalProjector.cpp
   hydro_NodalProjector.H
   )

if (HYDRO_FFT)
    target_sources(amrex_hydro PRIVATE
        hydro_FFTMacProjector.cpp
        hydro_FFTMacProjector.H)
endif()

AMREX_HOME ?= ../../../amrex
AMREX_HYDRO_HOME = ../..

USE_FFT  = TRUE
USE_MPI  = TRUE
USE_OMP  = FALSE

COMP = gnu

DIM = 3

DEBUG = FALSE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package

Pdirs := AmrCore
Pdirs += Base
Pdirs += Boundary
Pdirs += FFT
Pdirs += LinearSolvers/MLMG

Ppack	+= $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)
Ppack	+= $(AMREX_HYDRO_HOME)/Projections/Make.package

include $(Ppack)

Bdirs := AmrCore
Bdirs += Base
Bdirs += Boundary
Bdirs += FFT
Bdirs += LinearSolvers/MLMG

Blocs	:= $(foreach dir, $(Bdirs), $(AMREX_HOME)/Src/$(dir))
Blocs	+= $(AMREX_HYDRO_HOME)/Projections

INCLUDE_LOCATIONS += $(Blocs)
VPATH_LOCATIONS   += $(Blocs)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

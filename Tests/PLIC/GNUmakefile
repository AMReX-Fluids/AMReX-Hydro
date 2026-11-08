AMREX_HOME ?= ../../../amrex
AMREX_HYDRO_HOME ?= ../..

USE_MPI  = TRUE
USE_OMP  = FALSE

COMP = gnu

DIM = 2

DEBUG = FALSE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package

Pdirs := AmrCore
Pdirs += Base
Pdirs += Boundary

Ppack	+= $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)
Hdirs := Utils Godunov MOL BDS PLIC

Ppack	+= $(foreach dir, $(Hdirs), $(AMREX_HYDRO_HOME)/$(dir)/Make.package)

include $(Ppack)

Blocs	:= $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir))
Blocs	+= $(foreach dir, $(Hdirs), $(AMREX_HYDRO_HOME)/$(dir))

INCLUDE_LOCATIONS += $(Blocs)
VPATH_LOCATIONS   += $(Blocs)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

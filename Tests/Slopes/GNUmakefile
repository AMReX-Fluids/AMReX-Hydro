TEST = TRUE
USE_ASSERTION = TRUE

USE_EB = TRUE
DEBUG = FALSE
USE_MPI  = TRUE
USE_OMP  = FALSE

COMP = gnu

DIM = 3

AMREX_HOME ?= ../..

include $(AMREX_HOME)/Tools/GNUMake/Make.defs
include ./Make.package

Pdirs := Base Boundary AmrCore
Pdirs += EB

Ppack	+= $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)

include $(Ppack)

INCLUDE_LOCATIONS += ../../Slopes
VPATH_LOCATIONS   += ../../Slopes

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

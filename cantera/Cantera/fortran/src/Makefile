#/bin/sh
###############################################################
#  $Author: hkmoffa $
#  $Date: 2006/04/30 22:15:06 $
#  $Revision: 1.11 $
#
#  Copyright 2001 California Institute of Technology
#
###############################################################

.SUFFIXES :
.SUFFIXES :  .cpp .d .o .f90 .mod

CXX_FLAGS = -D_NO_MPI_VERSION -D_GNU_GCC_3 -D_NO_ICEMCFD_VERSION -I/nfs/kris/d1/people/charestm/Programs/CFFC/sparselib++/include -I/nfs/kris/d1/people/charestm/Programs/CFFC/sparselib++/iml/include -I/nfs/kris/d1/people/charestm/Programs/CFFC/sparselib++/mv/include -I/nfs/kris/d1/people/charestm/Programs/CFFC/bpkit/src  $(CXX_OPT)
FORT_FLAGS = -fno-second-underscore -I. -O3


CXX_OBJS = fct.o fctxml.o

INTERFACE_MODULE_OBJS = fct_interface.o fctxml_interface.o 

USER_MODULE_OBJS = cantera_xml.o cantera_thermo.o cantera_kinetics.o \
                   cantera_transport.o cantera_funcs.o cantera.o

MODULES = $(INTERFACE_MODULE_OBJS:_interface.o=.mod) $(USER_MODULE_OBJS)
OBJS = $(CXX_OBJS) $(USER_MODULE_OBJS)

DEPENDS = $(CXX_OBJS:.o=.d)
MODFILES = $(MODULES:.o=.mod)

# Fortran libraries
FORT_LIBS = 

# the C++ compiler
CXX = g++ -O

# the Fortran 90/95 compiler
F90 = /usr/bin/f95

# external libraries
EXT_LIBS =  -luser -loneD -lzeroD -ltransport -lcantera -lcvode -lctlapack -lctblas -lctmath -ltpx -lctf2c -lconverters 

# the directory where the Cantera libraries are located
CANTERA_LIBDIR=/nfs/kris/d1/people/charestm/Programs/CFFC/cantera/build/lib/i686-pc-linux-gnu

# ../../../lib
LIB_DEPS = $(CANTERA_LIBDIR)/libcantera.a \
           $(CANTERA_LIBDIR)/libtransport.a

# the directory where Cantera include files may be found.
CANTERA_INCDIR=../../src

# the directory where module .mod files should be put
MODULE_DIR = /nfs/kris/d1/people/charestm/Programs/CFFC/cantera/build/include/cantera

CXX_INCLUDES = -I$(CANTERA_INCDIR)  

# flags passed to the C++ compiler/linker for the linking step
LCXX_FLAGS = -L$(CANTERA_LIBDIR) -D_NO_MPI_VERSION -D_GNU_GCC_3 -D_NO_ICEMCFD_VERSION -I/nfs/kris/d1/people/charestm/Programs/CFFC/sparselib++/include -I/nfs/kris/d1/people/charestm/Programs/CFFC/sparselib++/iml/include -I/nfs/kris/d1/people/charestm/Programs/CFFC/sparselib++/mv/include -I/nfs/kris/d1/people/charestm/Programs/CFFC/bpkit/src 

# how to compile C++ source files to object files
%.o : %.cpp
	$(CXX) -c $< $(CXX_INCLUDES) $(CXX_FLAGS)

%.o : %.f90
	$(F90) -c $< $(FORT_FLAGS)

%.mod : %_interface.f90
	$(F90) -c $< $(FORT_FLAGS)

%.mod : %.f90
	$(F90) -c $< $(FORT_FLAGS)

LIB_NAME=libfct.a

FTLIB = /nfs/kris/d1/people/charestm/Programs/CFFC/cantera/build/lib/i686-pc-linux-gnu/$(LIB_NAME)

$(FTLIB): $(MODFILES) $(USER_MODULE_OBJS) $(CXX_OBJS) $(LIB_DEPS)
	$(RM)     $(FTLIB)
	ar ruv $(FTLIB) $(OBJS)

clean:
	$(RM) $(OBJS) $(FTLIB) $(MODFILES)
	(if test -d SunWS_cache ; then \
           $(RM) -rf SunWS_cache ; \
         fi )

install:
	/usr/bin/install -c $(FTLIB) /usr/local/cantera/lib
	/usr/bin/install -c $(MODFILES) /usr/local/cantera/include/cantera

win-install:
	/usr/bin/install -c ../../../lib/fct.lib /usr/local/cantera/lib

%.d: %.cpp
	g++ -MM $(CXX_INCLUDES) $*.cpp > $*.d

depends: $(DEPENDS)
	cat *.d > .depends
	$(RM) $(DEPENDS) 

cantera_thermo.o: fct.mod cantera_xml.mod cantera_thermo.f90
fct.mod: fct_interface.f90
fctxml.mod: fctxml_interface.f90

ifeq ($(wildcard .depends), .depends)
include .depends
endif

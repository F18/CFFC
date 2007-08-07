#!/bin/sh

############################################################################
#
#  Makefile to compile and link a C++ or Fortran application to 
#  Cantera.
#
#############################################################################

# the name of the executable program to be created
PROG_NAME = __PROGRAM__

# the object files to be linked together.
OBJS = __OBJS__

# additional flags to be passed to the linker. If your program
# requires other external libraries, put them here
LINK_OPTIONS = 


#############################################################################

# the Fortran compiler
FORT = gfortran

# Fortran compile flags  
FORT_FLAGS = -O  -fno-second-underscore 

# Fortran libraries
FORT_LIBS = 

# the C++ compiler
CXX = g++ -O

# C++ compile flags
CXX_FLAGS = -D_NO_MPI_VERSION -D_GNU_GCC_3 -D_NO_ICEMCFD_VERSION -I/nfs/kris/d1/people/charestm/Programs/CFFC/sparselib++/include -I/nfs/kris/d1/people/charestm/Programs/CFFC/sparselib++/iml/include -I/nfs/kris/d1/people/charestm/Programs/CFFC/sparselib++/mv/include -I/nfs/kris/d1/people/charestm/Programs/CFFC/bpkit/src 

# external libraries
EXT_LIBS =  -luser -loneD -lzeroD -ltransport -lcantera -lcvode -lctlapack -lctblas -lctmath -ltpx -lctf2c -lconverters 



#------  you probably don't have to change anything below this line -----


# the directory where the Cantera libraries are located
CANTERA_LIBDIR=@CANTERA_LIBDIR@

# required Cantera libraries
CANTERA_LIBS =  

# the directory where Cantera include files may be found.
CANTERA_INCDIR=@CANTERA_INCDIR@

# flags passed to the C++ compiler/linker for the linking step
LCXX_FLAGS = -L$(CANTERA_LIBDIR)  -L/nfs/kris/d1/people/charestm/Programs/CFFC/cantera/build/lib/i686-pc-linux-gnu -D_NO_MPI_VERSION -D_GNU_GCC_3 -D_NO_ICEMCFD_VERSION -I/nfs/kris/d1/people/charestm/Programs/CFFC/sparselib++/include -I/nfs/kris/d1/people/charestm/Programs/CFFC/sparselib++/iml/include -I/nfs/kris/d1/people/charestm/Programs/CFFC/sparselib++/mv/include -I/nfs/kris/d1/people/charestm/Programs/CFFC/bpkit/src 

# how to compile C++ source files to object files
.cpp.o:
	$(CXX) -c $< -I$(CANTERA_INCDIR) $(CXX_FLAGS)

# how to compile Fortran source files to object files
.f.o: 
	$(FORT) -c $< $(FORT_FLAGS)

PROGRAM = $(PROG_NAME)$(EXE_EXT)

DEPENDS = $(OBJS:.o=.d)

all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(CXX) -o $(PROGRAM) $(OBJS) $(LCXX_FLAGS) $(CANTERA_LIBS) $(LINK_OPTIONS) $(EXT_LIBS)  $(FORT_LIBS)

%.d:
	g++ -MM $*.cpp > $*.d

clean:
	$(RM) $(OBJS) $(PROGRAM)

depends: $(DEPENDS)
	cat *.d > .depends
	$(RM) $(DEPENDS) 

TAGS: 
	etags *.h *.cpp

ifeq ($(wildcard .depends), .depends)
include .depends
endif







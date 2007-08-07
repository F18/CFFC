#!/bin/sh


# the name of the executable program to be created
PROG_NAME = demo

# the object files to be linked together. 
OBJS = demo.o demo_ftnlib.o

# additional flags to be passed to the linker. If your program
# requires other external libraries, put them here
LINK_OPTIONS = 

#---------------------------------------------------------------------------
# You probably don't need to edit anything below.


# the Fortran compiler
FORT = gfortran

# Fortran compile flags  
FORT_FLAGS = -O  -fno-second-underscore 

# Fortran libraries
FORT_LIBS = 

# the C++ compiler
CXX = g++ -O

# C++ compile flags
CXX_FLAGS = -D_NO_MPI_VERSION -D_GNU_GCC_3 -D_NO_ICEMCFD_VERSION -I/home/groth/CFFC/sparselib++/include -I/home/groth/CFFC/sparselib++/iml/include -I/home/groth/CFFC/sparselib++/mv/include -I/home/groth/CFFC/bpkit/src 

# external libraries
EXT_LIBS =  -luser -loneD -lzeroD -ltransport -lcantera -lcvode -lctlapack -lctblas -lctmath -ltpx -lctf2c -lconverters -lctcxx

# the directory where the Cantera libraries are located
CANTERA_LIBDIR=/usr/local/cantera/lib

# the directory where Cantera include files may be found.
CANTERA_INCDIR=/usr/local/cantera/include

# flags passed to the C++ compiler/linker for the linking step
LCXX_FLAGS = -L$(CANTERA_LIBDIR) -D_NO_MPI_VERSION -D_GNU_GCC_3 -D_NO_ICEMCFD_VERSION -I/home/groth/CFFC/sparselib++/include -I/home/groth/CFFC/sparselib++/iml/include -I/home/groth/CFFC/sparselib++/mv/include -I/home/groth/CFFC/bpkit/src 

# how to compile C++ source files to object files
.cpp.o:
	$(CXX) -c $< -I$(CANTERA_INCDIR) $(CXX_FLAGS)

# how to compile Fortran source files to object files
.f.o: 
	$(FORT) -c $< $(FORT_FLAGS)

PROGRAM = $(PROG_NAME)$(EXE_EXT)

DEPENDS = $(OBJS:.o=.d)

all: isentropic ctlib 

isentropic: isentropic.o demo_ftnlib.o
	$(CXX) -o isentropic isentropic.o demo_ftnlib.o $(LCXX_FLAGS) $(CANTERA_LIBS) $(LINK_OPTIONS) $(EXT_LIBS)  $(FORT_LIBS)

ctlib: ctlib.o demo_ftnlib.o
	$(CXX) -o ctlib ctlib.o demo_ftnlib.o $(LCXX_FLAGS) $(CANTERA_LIBS) $(LINK_OPTIONS) $(EXT_LIBS)  $(FORT_LIBS)

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







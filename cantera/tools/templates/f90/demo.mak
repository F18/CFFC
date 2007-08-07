#!/bin/sh

# This Makefile builds a Fortran 90 application that uses Cantera.  By
# default, the main program file is 'demo.f90,' which prints out some
# properties of a reacting gas mixture. 

# To build program 'demo', simply type 'make', or 'make -f <this
# file>' if this file is named something other than 'Makefile.'  

# Once you have verified that the demo runs, edit this file to replace
# object file 'demo.o' with your own object file or files. 


#------------------------  edit this block ---------------------------------

# the name of the executable program to be created
PROG_NAME = demo

# the object files to be linked together. 
OBJS = demo.o 

# additional flags to be passed to the linker. If your program
# requires other external libraries, put them here
LINK_OPTIONS =  

#---------------------------------------------------------------------------
# You probably don't need to edit anything below.

# the C++ compiler
CXX = g++ -O

# C++ compile flags
CXX_FLAGS = -D_NO_MPI_VERSION -D_GNU_GCC_3 -D_NO_ICEMCFD_VERSION -I/nfs/kris/d1/people/charestm/Programs/CFFC/sparselib++/include -I/nfs/kris/d1/people/charestm/Programs/CFFC/sparselib++/iml/include -I/nfs/kris/d1/people/charestm/Programs/CFFC/sparselib++/mv/include -I/nfs/kris/d1/people/charestm/Programs/CFFC/bpkit/src 

# external libraries
EXT_LIBS = -lfct -lclib  -luser -loneD -lzeroD -ltransport -lcantera -lcvode -lctlapack -lctblas -lctmath -ltpx -lctf2c -lconverters -lctcxx -lstdc++ 

# the Fortran 90/95 compiler
F90 = /usr/bin/f95

# Fortran compile flags
FORT_FLAGS = -fno-second-underscore -I. -I/usr/local/cantera/include/cantera -O3

# the directory where the Cantera libraries are located
CANTERA_LIBDIR=/usr/local/cantera/lib

# the directory where Cantera include files may be found.
CANTERA_INCDIR=/usr/local/cantera/include

# the directory where Cantera Fortran 90 modules may be found.
CANTERA_MODULE_DIR=/usr/local/cantera/include/cantera

# flags passed to the C++ compiler/linker for the linking step
LCXXFLAGS = -L$(CANTERA_LIBDIR)  -L/nfs/kris/d1/people/charestm/Programs/CFFC/cantera/build/lib/i686-pc-linux-gnu

# how to compile C++ source files to object files
%.o : %.cpp
	$(CXX) -c $< -I$(CANTERA_INCDIR) $(CXX_FLAGS)

# how to compile Fortran 90/95 source files to object files
%.o : %.f90
	$(F90) -c $< -I$(CANTERA_MODULE_DIR) $(FORT_FLAGS)

PROGRAM = $(PROG_NAME)$(EXE_EXT)

DEPENDS = $(OBJS:.o=.d)

all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(F90) -o $(PROGRAM) $(OBJS) $(LCXXFLAGS) $(CANTERA_LIBS) $(LINK_OPTIONS) $(EXT_LIBS)  

%.d : %.cpp
	g++ -MM -I$(CANTERA_INCDIR) $*.cpp > $*.d

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







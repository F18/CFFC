#!/bin/sh

# This Makefile builds a Fortran 77 application that uses Cantera.  By
# default, the main program file is 'demo.f,' which prints out some
# properties of a reacting gas mixture. It uses the library
# 'demo_ftnlib.cpp,' which contains Fortran-callable functions that
# are implemented with C++ calls to Cantera.

# To build program 'demo', simply type 'make', or 'make -f <this
# file>' if this file is named something other than 'Makefile.'  

# Once you have verified that the demo runs, edit this file to replace
# object file 'demo.o' with your own object file or files. You may
# continue to use 'demo_ftnlib' if it serves your needs, or else
# replace it with a different interface library.


#------------------------  edit this block ---------------------------------

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
CXX_FLAGS = -D_NO_MPI_VERSION -D_GNU_GCC_3 -D_NO_ICEMCFD_VERSION -I/nfs/kris/d1/people/charestm/Programs/CFFC/sparselib++/include -I/nfs/kris/d1/people/charestm/Programs/CFFC/sparselib++/iml/include -I/nfs/kris/d1/people/charestm/Programs/CFFC/sparselib++/mv/include -I/nfs/kris/d1/people/charestm/Programs/CFFC/bpkit/src 

# external libraries
EXT_LIBS =  -luser -loneD -lzeroD -ltransport -lcantera -lcvode -lctlapack -lctblas -lctmath -ltpx -lctf2c -lconverters -lctcxx

# the directory where the Cantera libraries are located
CANTERA_LIBDIR=/usr/local/cantera/lib

# the directory where Cantera include files may be found.
CANTERA_INCDIR=/usr/local/cantera/include

# flags passed to the C++ compiler/linker for the linking step
LCXXFLAGS = -L$(CANTERA_LIBDIR)  -L/nfs/kris/d1/people/charestm/Programs/CFFC/cantera/build/lib/i686-pc-linux-gnu 

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
	$(CXX) -o $(PROGRAM) $(OBJS) $(LCXXFLAGS) $(CANTERA_LIBS) $(LINK_OPTIONS) $(EXT_LIBS)  $(FORT_LIBS)

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







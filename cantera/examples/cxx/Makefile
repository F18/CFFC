#!/bin/sh

############################################################################
#
#  Makefile to compile and link a C++ application to 
#  Cantera.
#
#############################################################################

# the name of the executable program to be created
PROG_NAME = cxx_examples

# the object files to be linked together. List those generated from Fortran
# and from C/C++ separately
OBJS = examples.o kinetics_example1.o kinetics_example2.o \
       kinetics_example3.o \
       equil_example1.o \
       transport_example1.o transport_example2.o \
       rxnpath_example1.o

# additional flags to be passed to the linker. If your program
# requires other external libraries, put them here
LINK_OPTIONS = 


#############################################################################

# the Fortran compiler
FORT = /usr/bin/f95

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

# Ending C++ linking libraries
LCXX_END_LIBS =  -lctf2c -lm


#------  you probably don't have to change anything below this line -----


# the directory where the Cantera libraries are located
CANTERA_LIBDIR=/home/groth/CFFC/cantera/build/lib/i686-pc-linux-gnu

# required Cantera libraries
CANTERA_LIBS = 

# the directory where Cantera include files may be found.
CANTERA_INC=-I/home/groth/CFFC/cantera/build/include

# flags passed to the C++ compiler/linker for the linking step
LCXX_FLAGS = -L$(CANTERA_LIBDIR)  -L/home/groth/CFFC/cantera/build/lib/i686-pc-linux-gnu -D_NO_MPI_VERSION -D_GNU_GCC_3 -D_NO_ICEMCFD_VERSION -I/home/groth/CFFC/sparselib++/include -I/home/groth/CFFC/sparselib++/iml/include -I/home/groth/CFFC/sparselib++/mv/include -I/home/groth/CFFC/bpkit/src 

# how to compile C++ source files to object files
.cpp.o:
	$(CXX) -c $< $(CANTERA_INC) $(CXX_FLAGS)

# how to compile Fortran source files to object files
.f.o: 
	$(FORT) -c $< $(FORT_FLAGS)

PROGRAM = ./$(PROG_NAME)$(EXE_EXT)

all: $(PROGRAM)

DEPENDS=$(OBJS:.o=.d)

%.d:
	g++ -MM -I$(CANTERA_INCDIR) $(CXX_FLAGS) $*.cpp > $*.d

$(PROGRAM): $(OBJS) $(CANTERA_LIBDIR)/libcantera.a $(CANTERA_LIBDIR)/libzeroD.a
	$(CXX) -o $(PROGRAM) $(OBJS) $(LCXX_FLAGS) $(CANTERA_LIBS) \
                   $(LINK_OPTIONS) $(EXT_LIBS)  \
                  $(LCXX_END_LIBS)
# $(FORT_LIBS) \

clean:
	$(RM) $(OBJS) $(PROGRAM)
	(if test -d SunWS_cache ; then \
           $(RM) -rf SunWS_cache ; \
         fi )


depends: $(DEPENDS)
	cat *.d > .depends
	$(RM) $(DEPENDS)

TAGS:
	etags *.h *.cpp

ifeq ($(wildcard .depends), .depends)
include .depends
endif


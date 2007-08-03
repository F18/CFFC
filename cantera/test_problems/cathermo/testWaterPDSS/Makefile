#!/bin/sh

############################################################################
#
#  Makefile to compile and link a C++ application to 
#  Cantera.
#
#############################################################################

# addition to suffixes
.SUFFIXES : .d

# the name of the executable program to be created
PROG_NAME = testWaterPDSS

# the object files to be linked together. List those generated from Fortran
# and from C/C++ separately
OBJS = testWaterPDSS.o

# Location of the current build. Will assume that tests are run
# in the source directory tree location
src_dir_tree = 1

# additional flags to be passed to the linker. If your program
# requires other external libraries, put them here
LINK_OPTIONS = 

#############################################################################

# Check to see whether we are in the msvc++ environment
os_is_win = 0

# Fortran libraries
FORT_LIBS = 

# the C++ compiler
CXX = g++ -O

# C++ compile flags
ifeq ($(src_dir_tree), 1)
CXX_FLAGS = -DSRCDIRTREE -D_NO_MPI_VERSION -D_GNU_GCC_3 -D_NO_ICEMCFD_VERSION -I/home/groth/CFFC/sparselib++/include -I/home/groth/CFFC/sparselib++/iml/include -I/home/groth/CFFC/sparselib++/mv/include -I/home/groth/CFFC/bpkit/src 
else
CXX_FLAGS = -D_NO_MPI_VERSION -D_GNU_GCC_3 -D_NO_ICEMCFD_VERSION -I/home/groth/CFFC/sparselib++/include -I/home/groth/CFFC/sparselib++/iml/include -I/home/groth/CFFC/sparselib++/mv/include -I/home/groth/CFFC/bpkit/src 
endif

# Ending C++ linking libraries
LCXX_END_LIBS =  -lctf2c -lm

# the directory where the Cantera libraries are located
CANTERA_LIBDIR=/home/groth/CFFC/cantera/build/lib/i686-pc-linux-gnu

# required Cantera libraries
CANTERA_LIBS =  -luser -loneD -lzeroD -ltransport -lcantera -lcvode -lctlapack -lctblas -lctmath -ltpx -lctf2c -lconverters -lctcxx

# the directory where Cantera include files may be found.
ifeq ($(src_dir_tree), 1)
CANTERA_INCDIR=../../../Cantera/src
INCLUDES=-I$(CANTERA_INCDIR) -I$(CANTERA_INCDIR)/thermo
else
CANTERA_INCDIR=/home/groth/CFFC/cantera/build/include/cantera
INCLUDES=-I$(CANTERA_INCDIR) -I$(CANTERA_INCDIR)/kernel
endif

# flags passed to the C++ compiler/linker for the linking step
LCXX_FLAGS = -L$(CANTERA_LIBDIR)  -L/home/groth/CFFC/cantera/build/lib/i686-pc-linux-gnu -D_NO_MPI_VERSION -D_GNU_GCC_3 -D_NO_ICEMCFD_VERSION -I/home/groth/CFFC/sparselib++/include -I/home/groth/CFFC/sparselib++/iml/include -I/home/groth/CFFC/sparselib++/mv/include -I/home/groth/CFFC/bpkit/src 

# How to compile C++ source files to object files
.cpp.o:
	$(CXX) -c $< $(INCLUDES) $(CXX_FLAGS)

# How to compile the dependency file
.cpp.d:
	g++ -MM $(INCLUDES) $(CXX_FLAGS) $*.cpp > $*.d

# List of dependency files to be created
DEPENDS=$(OBJS:.o=.d)

# Program Name
PROGRAM = $(PROG_NAME)$(EXE_EXT)

all: $(PROGRAM) .depends

$(PROGRAM): $(OBJS) $(CANTERA_LIBDIR)/libcantera.a \
	     $(CANTERA_LIBDIR)/libcaThermo.a
	$(CXX) -o $(PROGRAM) $(OBJS) $(LCXX_FLAGS) $(LINK_OPTIONS) \
                  $(CANTERA_LIBS)  $(FORT_LIBS) \
                  $(LCXX_END_LIBS)

# depends target -> forces recalculation of dependencies
depends:
	make .depends

.depends: $(DEPENDS)
	cat $(DEPENDS) > .depends

# Do the test -> For the windows vc++ environment, we have to skip checking on
#                whether the program is uptodate, because we don't utilize make
#                in that environment to build programs.
test:
ifeq ($(os_is_win), 1)
else
	@ make -s $(PROGRAM)
endif
	@ ./runtest

clean:
	$(RM) $(OBJS) $(PROGRAM) $(DEPENDS) .depends
	../../../bin/rm_cvsignore
	(if test -d SunWS_cache ; then \
           $(RM) -rf SunWS_cache ; \
         fi )


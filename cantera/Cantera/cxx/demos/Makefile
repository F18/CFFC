#!/bin/sh

############################################################################
#
#  Makefile to compile and link a C++ application to 
#  Cantera.
#
#############################################################################
.SUFFIXES :
.SUFFIXES : .cpp .d .o .dh .h .h.gch


# the name of the executable program to be created
PROG_NAME = demos

# the object files to be linked together. List those generated from Fortran
# and from C/C++ separately
OBJS = demos.o

# additional flags to be passed to the linker. If your program
# requires other external libraries, put them here
LINK_OPTIONS = 

#############################################################################
# True if we are in the source directory tree
srcdirtree=1

# Fortran libraries
FORT_LIBS = 

# the C++ compiler
CXX = g++ -O

# C++ compile flags
CXX_FLAGS = -D_NO_MPI_VERSION -D_GNU_GCC_3 -D_NO_ICEMCFD_VERSION -I/home/groth/CFFC/sparselib++/include -I/home/groth/CFFC/sparselib++/iml/include -I/home/groth/CFFC/sparselib++/mv/include -I/home/groth/CFFC/bpkit/src  

# external libraries
EXT_LIBS =  -luser -loneD -lzeroD -ltransport -lcantera -lcvode -lctlapack -lctblas -lctmath -ltpx -lctf2c -lconverters -ltpx -lctcxx

# Ending C++ linking libraries
LCXX_END_LIBS =  -lctf2c -lm

# the directory where the Cantera libraries are located
CANTERA_LIBDIR=/home/groth/CFFC/cantera/build/lib/i686-pc-linux-gnu

# required Cantera libraries
CANTERA_LIBS =  

# Depends
ifeq ($srcdirtree, 1)
  LOCAL_DEFNS =  -DSRCDIRTREE
else
  LOCAL_DEFNS = 
endif 

# the directory where Cantera include files may be found.
CANTERA_INCDIR=/home/groth/CFFC/cantera/build/include

# flags passed to the C++ compiler/linker for the linking step
LCXX_FLAGS = -L$(CANTERA_LIBDIR)  -L/home/groth/CFFC/cantera/build/lib/i686-pc-linux-gnu -D_NO_MPI_VERSION -D_GNU_GCC_3 -D_NO_ICEMCFD_VERSION -I/home/groth/CFFC/sparselib++/include -I/home/groth/CFFC/sparselib++/iml/include -I/home/groth/CFFC/sparselib++/mv/include -I/home/groth/CFFC/bpkit/src 

# how to compile C++ source files to object files
.cpp.o:
	$(CXX) -c $< -I$(CANTERA_INCDIR) $(CXX_FLAGS) $(LOCAL_DEFNS)

# how to create a dependency file
.cpp.d:
	g++ -MM -I$(CANTERA_INCDIR) $(CXX_FLAGS) $(LOCAL_DEFNS) $*.cpp > $*.d

PROGRAM = $(PROG_NAME)$(EXE_EXT)

DEPENDS = $(OBJS:.o=.d)

all: $(PROGRAM) .depends

$(PROGRAM): $(OBJS) $(CANTERA_LIBDIR)/libcantera.a
	$(CXX) -o $(PROGRAM) $(OBJS) $(LCXX_FLAGS) $(CANTERA_LIBS) \
                   $(LINK_OPTIONS) $(EXT_LIBS)  $(FORT_LIBS) \
                  $(LCXX_END_LIBS)

test:
	make $(PROGRAM)
	echo 100 | $(PROGRAM) > $(PROGRAM).out

install:
	/usr/bin/install -c -d       /usr/local/cantera/demos/cxx
	/usr/bin/install -c *.cpp    /usr/local/cantera/demos/cxx
	/usr/bin/install -c *.h      /usr/local/cantera/demos/cxx
	/usr/bin/install -c Makefile.install /usr/local/cantera/demos/cxx/Makefile

depends: $(DEPENDS)
	make .depends


.depends: $(DEPENDS)
	cat $(DEPENDS) > .depends

clean:
	$(RM) $(OBJS) $(PROGRAM) .depends *.d
	../../../bin/rm_cvsignore
	(if test -d SunWS_cache ; then \
           $(RM) -rf SunWS_cache ; \
         fi )

TAGS:
	etags *.h *.cpp

ifeq ($(wildcard .depends), .depends)
include .depends
endif


#/bin/sh

LIBDIR = /home/groth/CFFC/cantera/build/lib/i686-pc-linux-gnu
INCDIR = /home/groth/CFFC/cantera/build/include/cantera
BINDIR = /home/groth/CFFC/cantera/bin
build_ck = 1

LCXX_FLAGS = -L$(LIBDIR)   -L/home/groth/CFFC/cantera/build/lib/i686-pc-linux-gnu -D_NO_MPI_VERSION -D_GNU_GCC_3 -D_NO_ICEMCFD_VERSION -I/home/groth/CFFC/sparselib++/include -I/home/groth/CFFC/sparselib++/iml/include -I/home/groth/CFFC/sparselib++/mv/include -I/home/groth/CFFC/bpkit/src 
LOCAL_LIBS = -lcantera -lckreader -lcvode -lctmath @LAPACK_LIBRARY@ @BLAS_LIBRARY@ -lctcxx

LCXX_END_LIBS =  -lctf2c -lm

OBJS = mdp_allo.o csvdiff.o tok_input_util.o
DEPENDS = $(OBJS:.o=.d)

CXX_FLAGS=-D_NO_MPI_VERSION -D_GNU_GCC_3 -D_NO_ICEMCFD_VERSION -I/home/groth/CFFC/sparselib++/include -I/home/groth/CFFC/sparselib++/iml/include -I/home/groth/CFFC/sparselib++/mv/include -I/home/groth/CFFC/bpkit/src 
CONFIGDEFS=-DHAVE_CONFIG_H

# Optimization is turned off (-O0) because a strange interal compiler
# error results when using g++ 3.0 on Mac OS X. Optimization should
# not be needed here anyway.

.cpp.o:
	g++ -O -c $< -I$(INCDIR) -D_NO_MPI_VERSION -D_GNU_GCC_3 -D_NO_ICEMCFD_VERSION -I/home/groth/CFFC/sparselib++/include -I/home/groth/CFFC/sparselib++/iml/include -I/home/groth/CFFC/sparselib++/mv/include -I/home/groth/CFFC/bpkit/src  -O0

all: $(BINDIR)/csvdiff

$(BINDIR)/csvdiff: mdp_allo.o csvdiff.o tok_input_util.o
	g++ -O -o $(BINDIR)/csvdiff  mdp_allo.o csvdiff.o tok_input_util.o \
	 $(LCXX_FLAGS) $(LCXX_END_LIBS) 

clean:
	$(RM) *.o *.*~ csvdiff $(BINDIR)/csvdiff

%.d:
	g++ -O -MM -I$(INCDIR) $(CXX_FLAGS) $*.cpp > $*.d

depends: $(DEPENDS)
	cat *.d &> .depends
	$(RM) $(DEPENDS)

TAGS: 
	etags *.h *.cpp

ifeq ($(wildcard .depends), .depends)
include .depends
endif

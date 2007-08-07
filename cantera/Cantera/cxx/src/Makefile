#/bin/sh
###############################################################
#  $Author: hkmoffa $
#  $Date: 2006/04/30 22:10:37 $
#  $Revision: 1.5 $
#
#  Copyright 2001 California Institute of Technology
#
###############################################################


SUFFIXES=
SUFFIXES= .cpp .d .o

CXX_FLAGS = -D_NO_MPI_VERSION -D_GNU_GCC_3 -D_NO_ICEMCFD_VERSION -I/nfs/kris/d1/people/charestm/Programs/CFFC/sparselib++/include -I/nfs/kris/d1/people/charestm/Programs/CFFC/sparselib++/iml/include -I/nfs/kris/d1/people/charestm/Programs/CFFC/sparselib++/mv/include -I/nfs/kris/d1/people/charestm/Programs/CFFC/bpkit/src  $(CXX_OPT)

OBJS = cxxutils.o 

DEPENDS = $(OBJS:.o=.d)

# the C++ compiler
CXX = g++ -O

# the directory where Cantera include files may be found.
CANTERA_INCDIR=../../src

CXX_INCLUDES = -I$(CANTERA_INCDIR)  

# how to compile C++ source files to object files
.cpp.o:
	$(CXX) -c $< $(CXX_INCLUDES) $(CXX_FLAGS)

LIB_NAME=libctcxx
CXXLIB=/nfs/kris/d1/people/charestm/Programs/CFFC/cantera/build/lib/i686-pc-linux-gnu/$(LIB_NAME).a

all: $(CXXLIB)

$(CXXLIB): $(OBJS)
	ar ruv $(CXXLIB) $(OBJS) 

clean:
	$(RM) $(OBJS) $(CXXLIB)
	(if test -d SunWS_cache ; then \
           $(RM) -rf SunWS_cache ; \
         fi )

install:
	/usr/bin/install -c $(CXXLIB) /usr/local/cantera/lib/cantera

%.d:
	g++ -MM $(CXX_INCLUDES) $*.cpp > $*.d

depends: $(DEPENDS)
	cat *.d > .depends
	$(RM) $(DEPENDS) 

ifeq ($(wildcard .depends), .depends)
include .depends
endif

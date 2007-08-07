#/bin/sh
###############################################################
#  $Author: hkmoffa $
#  $Date: 2006/08/21 16:42:16 $
#  $Revision: 1.3 $
#
#  Copyright 2002 California Institute of Technology
#
###############################################################

###############################################################
# list your object files here
# This is probably the only thing you need to change in this file.

OBJS    = user.o

###############################################################



SUFFIXES=
SUFFIXES= .cpp .d .o

OBJDIR = .

CXX_FLAGS = -D_NO_MPI_VERSION -D_GNU_GCC_3 -D_NO_ICEMCFD_VERSION -I/home/groth/CFFC/sparselib++/include -I/home/groth/CFFC/sparselib++/iml/include -I/home/groth/CFFC/sparselib++/mv/include -I/home/groth/CFFC/bpkit/src  $(CXX_OPT)


CXX_INCLUDES = -I../src

LIB = /home/groth/CFFC/cantera/build/lib/i686-pc-linux-gnu/libuser.a

DEPENDS = $(OBJS:.o=.d)

%.d:
	g++ -MM $(CXX_INCLUDES) $*.cpp > $*.d

.cpp.o:
	g++ -O -c $< $(CXX_FLAGS) $(CXX_INCLUDES)

.f.o:
	gfortran -c $< $(F77_FLAGS) 

all lib: $(LIB)

$(LIB): $(OBJS)
	ar ruv $(LIB) $(OBJS) > /dev/null

clean:
	$(RM) *.o *~ $(LIB)

depends: $(DEPENDS)
	cat *.d > .depends
	$(RM) $(DEPENDS) 

install:


TAGS: 
	etags *.h *.cpp

ifeq ($(wildcard .depends), .depends)
include .depends
endif


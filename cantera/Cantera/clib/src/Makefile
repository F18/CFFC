#/bin/sh
###############################################################
#  $Author: hkmoffa $
#  $Date: 2006/04/30 22:09:32 $
#  $Revision: 1.11 $
#
#  Copyright 2001 California Institute of Technology
#
###############################################################


.SUFFIXES :
.SUFFIXES : .cpp .d .o

INSTALL_TSC = ../../../bin/install_tsc
do_ranlib = 1

CXX_FLAGS = -D_NO_MPI_VERSION -D_GNU_GCC_3 -D_NO_ICEMCFD_VERSION -I/home/groth/CFFC/sparselib++/include -I/home/groth/CFFC/sparselib++/iml/include -I/home/groth/CFFC/sparselib++/mv/include -I/home/groth/CFFC/bpkit/src  $(CXX_OPT)

OBJS = ct.o Storage.o ctsurf.o ctrpath.o \
       ctreactor.o ctfunc.o ctxml.o ctonedim.o ctmultiphase.o
CLIB_H= Cabinet.h    ctreactor.h  clib_defs.h  ctfunc.h \
        ctnum.h ctsurf.h  ct.h   ctrpath.h      Storage.h \
        ctbdry.h    ctonedim.h    ctxml.h   ctmultiphase.h  ctstagn.h   

DEPENDS = $(OBJS:.o=.d)

# Fortran libraries
FORT_LIBS = 

shared_ctlib = 0

# the C++ compiler
CXX = g++ -O

# external libraries
EXT_LIBS = -lzeroD -loneD  -luser -loneD -lzeroD -ltransport -lcantera -lcvode -lctlapack -lctblas -lctmath -ltpx -lctf2c -lconverters 

# the directory where the Cantera libraries are located
CANTERA_LIBDIR=/home/groth/CFFC/cantera/build/lib/i686-pc-linux-gnu
# ../../../lib
LIB_DEPS = $(CANTERA_LIBDIR)/libcantera.a \
           $(CANTERA_LIBDIR)/libzeroD.a \
           $(CANTERA_LIBDIR)/liboneD.a 

# the directory where Cantera include files may be found.
CANTERA_INCDIR=../../src

CXX_INCLUDES = -I$(CANTERA_INCDIR)  

# flags passed to the C++ compiler/linker for the linking step
LCXX_FLAGS = -L$(CANTERA_LIBDIR) -D_NO_MPI_VERSION -D_GNU_GCC_3 -D_NO_ICEMCFD_VERSION -I/home/groth/CFFC/sparselib++/include -I/home/groth/CFFC/sparselib++/iml/include -I/home/groth/CFFC/sparselib++/mv/include -I/home/groth/CFFC/bpkit/src 

# how to compile C++ source files to object files
.cpp.o:
	$(CXX) -c $< $(CXX_INCLUDES) $(CXX_FLAGS) -fPIC

LIB_NAME=libclib

ifeq ($(shared_ctlib),1)
CTLIB = /home/groth/CFFC/cantera/build/lib/i686-pc-linux-gnu/$(LIB_NAME).so
else
CTLIB = /home/groth/CFFC/cantera/build/lib/i686-pc-linux-gnu/$(LIB_NAME).a
endif

all: $(CTLIB)
 
$(CTLIB): $(OBJS) $(LIB_DEPS)
	$(RM)     $(CTLIB)
ifeq ($(shared_ctlib),1)
	$(CXX) -o $(CTLIB) $(OBJS) $(LCXX_FLAGS)\
                -shared $(LINK_OPTIONS) \
                $(EXT_LIBS)  $(FORT_LIBS)
else
	ar ruv $(CTLIB) $(OBJS)
ifeq ($(do_ranlib),1)
	ranlib $(CTLIB)
endif
endif

clean:
	$(RM) $(OBJS) $(CTLIB)
	(if test -d SunWS_cache ; then \
           $(RM) -rf SunWS_cache ; \
         fi )


install:
	/usr/bin/install -c -d /usr/local/cantera/lib
	/usr/bin/install -c -m 644 $(CTLIB) /usr/local/cantera/lib
ifeq ($(do_ranlib),1)
ifeq ($(shared_ctlib),0)
	ranlib /usr/local/cantera/lib/$(LIB_NAME).a
endif
endif


win-install:
	/usr/bin/install -c ../../../lib/clib.lib /usr/local/cantera/lib/cantera


%.d:
	g++ -MM $(CXX_INCLUDES) $*.cpp > $*.d

depends: $(DEPENDS)
	cat *.d > .depends
	$(RM) $(DEPENDS) 

ifeq ($(wildcard .depends), .depends)
include .depends
endif

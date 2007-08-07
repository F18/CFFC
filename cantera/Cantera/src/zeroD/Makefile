#/bin/sh
###############################################################
#  $Author: dggoodwin $
#  $Date: 2006/11/11 13:23:11 $
#  $Revision: 1.9 $
#
#  Copyright 2002 California Institute of Technology
#
###############################################################

.SUFFIXES :
.SUFFIXES : .cpp .d .o .h

INCDIR = ../../../build/include/cantera/kernel/zeroD
INSTALL_TSC = ../../../bin/install_tsc
do_ranlib = 1

CXX_FLAGS = -D_NO_MPI_VERSION -D_GNU_GCC_3 -D_NO_ICEMCFD_VERSION -I/nfs/kris/d1/people/charestm/Programs/CFFC/sparselib++/include -I/nfs/kris/d1/people/charestm/Programs/CFFC/sparselib++/iml/include -I/nfs/kris/d1/people/charestm/Programs/CFFC/sparselib++/mv/include -I/nfs/kris/d1/people/charestm/Programs/CFFC/bpkit/src  $(CXX_OPT)

# stirred reactors
OBJS    = Reactor.o ReactorBase.o FlowDevice.o Wall.o ReactorNet.o \
          FlowReactor.o ConstPressureReactor.o ReactorFactory.o
ZEROD_H = Reactor.h ReactorBase.h FlowDevice.h Wall.h ReactorNet.h \
          flowControllers.h PID_Controller.h Reservoir.h FlowReactor.h \
          ConstPressureReactor.h ReactorFactory.h

CXX_INCLUDES = -I..
ZEROD_LIB = /nfs/kris/d1/people/charestm/Programs/CFFC/cantera/build/lib/i686-pc-linux-gnu/libzeroD.a

DEPENDS = $(OBJS:.o=.d)

all: .depends $(ZEROD_LIB)
	@(/usr/bin/install -c -d $(INCDIR))
	@(for lh in $(ZEROD_H) ; do \
          $(INSTALL_TSC) "$${lh}" $(INCDIR) ; \
	done)

%.d:
	g++ -MM $(CXX_INCLUDES) $*.cpp > $*.d

.cpp.o:
	g++ -O -c $< -DHAVE_CONFIG_H $(CXX_FLAGS) $(CXX_INCLUDES)

.f.o:
	gfortran -c $< $(F77_FLAGS) 

$(ZEROD_LIB): $(OBJS) $(ZEROD_H)
	ar ruv $(ZEROD_LIB) $(OBJS) > /dev/null
ifeq ($(do_ranlib),1)
	ranlib $(ZEROD_LIB)
endif

clean:
	@(for lh in $(ZEROD_H) ; do \
            th=$(INCDIR)/"$${lh}" ; \
            if test -f "$${th}" ; then \
              $(RM) "$${th}"  ; \
              echo "$(RM) $${th}" ; \
            fi \
          done)
	@(if test -f $(ZEROD_LIB) ; then \
            $(RM) $(ZEROD_LIB) ; \
            echo "$(RM) $(ZEROD_LIB)" ; \
          fi)
	$(RM) *.o *~ .depends
	(if test -d SunWS_cache ; then \
           $(RM) -rf SunWS_cache ; \
         fi )

depends: $(DEPENDS)
	cat *.d > .depends
	$(RM) $(DEPENDS) 

.depends: $(DEPENDS) Makefile
	cat *.d > .depends

TAGS: 
	etags *.h *.cpp

ifeq ($(wildcard .depends), .depends)
include .depends
endif


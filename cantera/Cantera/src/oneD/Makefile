#/bin/sh
###############################################################
#  $Author: dggoodwin $
#  $Date: 2006/04/23 07:24:27 $
#  $Revision: 1.8 $
#
#  Copyright 2002 California Institute of Technology
#
###############################################################

.SUFFIXES :
.SUFFIXES : .cpp .d .o

INCDIR = ../../../build/include/cantera/kernel/oneD
INSTALL_TSC = ../../../bin/install_tsc                      
do_ranlib = 1
CXX_FLAGS = -D_NO_MPI_VERSION -D_GNU_GCC_3 -D_NO_ICEMCFD_VERSION -I/nfs/kris/d1/people/charestm/Programs/CFFC/sparselib++/include -I/nfs/kris/d1/people/charestm/Programs/CFFC/sparselib++/iml/include -I/nfs/kris/d1/people/charestm/Programs/CFFC/sparselib++/mv/include -I/nfs/kris/d1/people/charestm/Programs/CFFC/bpkit/src  $(CXX_OPT)
CXX_INCLUDES = -I.. 

# stirred reactors
OBJS    = MultiJac.o MultiNewton.o newton_utils.o OneDim.o\
          StFlow.o boundaries1D.o refine.o Sim1D.o Domain1D.o
ONED_H  = Inlet1D.h   MultiJac.h  Sim1D.h         StFlow.h \
          Surf1D.h    Domain1D.h  MultiNewton.h   OneDim.h \
          Resid1D.h   Solid1D.h   refine.h

ONED_LIB = /nfs/kris/d1/people/charestm/Programs/CFFC/cantera/build/lib/i686-pc-linux-gnu/liboneD.a

DEPENDS = $(OBJS:.o=.d)

%.d:
	g++ -MM $(CXX_INCLUDES) $*.cpp > $*.d

.cpp.o:
	g++ -O -c $< -DHAVE_CONFIG_H $(CXX_FLAGS) $(CXX_INCLUDES)

all: $(ONED_LIB)
	@(/usr/bin/install -c -d $(INCDIR))
	@(for lh in $(ONED_H) ; do \
          $(INSTALL_TSC) "$${lh}" $(INCDIR) ; \
        done)

$(ONED_LIB): $(OBJS)
	ar ruv $(ONED_LIB) $(OBJS) > /dev/null
ifeq ($(do_ranlib),1)
	ranlib $(ONED_LIB)
endif

clean:        
	@(for lh in $(ONED_H) ; do \
            th=$(INCDIR)/"$${lh}" ; \
            if test -f "$${th}" ; then \
              $(RM) "$${th}"  ; \
              echo "$(RM) $${th}" ; \
            fi \
          done)
	$(RM) *.o *~ $(ONED_LIB) .depends
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


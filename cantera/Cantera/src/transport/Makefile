#/bin/sh
###############################################################
#  $Author: dggoodwin $
#  $Date: 2006/04/28 17:22:24 $
#  $Revision: 1.10 $
#
#  Copyright 2002 California Institute of Technology
#
###############################################################

.SUFFIXES :
.SUFFIXES : .cpp .d .o .h

INCDIR = ../../../build/include/cantera/kernel/transport
INSTALL_TSC = ../../../bin/install_tsc
do_ranlib = 1

CXX_FLAGS = -D_NO_MPI_VERSION -D_GNU_GCC_3 -D_NO_ICEMCFD_VERSION -I/home/groth/CFFC/sparselib++/include -I/home/groth/CFFC/sparselib++/iml/include -I/home/groth/CFFC/sparselib++/mv/include -I/home/groth/CFFC/bpkit/src  $(CXX_OPT)

# Transport Object Files
OBJS    = TransportFactory.o MultiTransport.o MixTransport.o MMCollisionInt.o \
          SolidTransport.o DustyGasTransport.o
TRAN_H  =  TransportFactory.h MultiTransport.h MixTransport.h \
           MMCollisionInt.h  SolidTransport.h DustyGasTransport.h \
           TransportBase.h L_matrix.h TransportParams.h

CXX_INCLUDES = -I.. 
LIB = /home/groth/CFFC/cantera/build/lib/i686-pc-linux-gnu/libtransport.a

DEPENDS = $(OBJS:.o=.d)

all: $(LIB)
	@(/usr/bin/install -c -d $(INCDIR))
	@(for lh in $(TRAN_H) ; do \
	  $(INSTALL_TSC) "$${lh}" $(INCDIR) ; \
	done)
%.d:
	g++ -MM $(CXX_INCLUDES) $*.cpp > $*.d

.cpp.o:
	g++ -O -c $< $(CXX_FLAGS) $(CXX_INCLUDES)

$(LIB): $(OBJS) $(TRAN_H)
	ar ruv $(LIB) $(OBJS) > /dev/null
ifeq ($(do_ranlib),1)
	ranlib $(LIB)
endif

clean:
	@(for lh in $(TRAN_H) ; do \
            th=$(INCDIR)/"$${lh}" ; \
	    if test -f "$${th}" ; then \
	      $(RM) "$${th}"  ; \
              echo "$(RM) $${th}" ; \
            fi \
	  done)
	@(if test -f $(LIB) ; then \
	    $(RM) $(LIB) ; \
            echo "$(RM) $(LIB)" ; \
          fi)
	$(RM) *.o *~  .depends
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


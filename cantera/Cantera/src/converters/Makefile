#/bin/sh
###############################################################
#  $Author: hkmoffa $
#  $Date: 2006/04/30 22:00:53 $
#  $Revision: 1.10 $
#
#  Copyright 2002 California Institute of Technology
#
###############################################################

.SUFFIXES :
.SUFFIXES : .cpp .d .o .h

INCDIR = ../../../build/include/cantera/kernel/converters
INSTALL_TSC = ../../../bin/install_tsc
do_ranlib = 1

CXX_FLAGS = -D_NO_MPI_VERSION -D_GNU_GCC_3 -D_NO_ICEMCFD_VERSION -I/home/groth/CFFC/sparselib++/include -I/home/groth/CFFC/sparselib++/iml/include -I/home/groth/CFFC/sparselib++/mv/include -I/home/groth/CFFC/bpkit/src  $(CXX_OPT) 

OBJS = atomicWeightDB.o CKParser.o CKReader.o Reaction.o ckr_utils.o \
       thermoFunctions.o writelog.o ck2ct.o 
# ck2ctml.o
CONV_H = CKReader.h      thermoFunctions.h \
         Element.h    Reaction.h      CKParser.h    \
         ckr_utils.h  RxnSpecies.h    writelog.h \
         ck2ct.h      ckr_defs.h      Constituent.h   \
         Species.h

CXX_INCLUDES = -I. -I..
CONV_LIB = /home/groth/CFFC/cantera/build/lib/i686-pc-linux-gnu/libconverters.a

DEPENDS = $(OBJS:.o=.d)

all: $(CONV_LIB)
	@(/usr/bin/install -c -d $(INCDIR))
	@(for lh in $(CONV_H) ; do \
          $(INSTALL_TSC) "$${lh}" $(INCDIR) ; \
        done)

%.d:
	g++ -MM $(CXX_INCLUDES) $*.cpp > $*.d

.cpp.o:
	g++ -O -c $< $(CXX_FLAGS) $(CXX_INCLUDES)

$(CONV_LIB): $(OBJS)
	ar ruv $(CONV_LIB) $(OBJS) 
ifeq ($(do_ranlib),1)
	ranlib $(CONV_LIB)
endif

clean:
	@(for lh in $(CONV_H) ; do \
            th=$(INCDIR)/"$${lh}" ; \
            if test -f "$${th}" ; then \
              $(RM) "$${th}"  ; \
              echo "$(RM) $${th}" ; \
            fi \
          done)
	$(RM) *.o *~ $(CONV_LIB) .depends
	(if test -d SunWS_cache ; then \
           $(RM) -rf SunWS_cache ; \
         fi )

CKParser.o: CKParser.cpp
	g++ -O -c CKParser.cpp $(CXX_FLAGS) $(CXX_INCLUDES) -O0

depends: $(DEPENDS)
	cat *.d > .depends
	$(RM) $(DEPENDS) 

TAGS: 
	etags *.h *.cpp

ifeq ($(wildcard .depends), .depends)
include .depends
endif


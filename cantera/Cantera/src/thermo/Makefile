#/bin/sh
###############################################################
#  $Author: hkmoffa $
#  $Date: 2006/08/15 21:59:24 $
#  $Revision: 1.11 $
#
#  Copyright 2002 California Institute of Technology
#
###############################################################

.SUFFIXES :
.SUFFIXES : .cpp .d .o .h

INCDIR = ../../../build/include/cantera/kernel/thermo
INSTALL_TSC = ../../../bin/install_tsc
do_ranlib = 1
do_electro = 0
do_issp    = 0

LOCAL_DEFS = -DDEBUG_MODE
CXX_FLAGS = -D_NO_MPI_VERSION -D_GNU_GCC_3 -D_NO_ICEMCFD_VERSION -I/home/groth/CFFC/sparselib++/include -I/home/groth/CFFC/sparselib++/iml/include -I/home/groth/CFFC/sparselib++/mv/include -I/home/groth/CFFC/bpkit/src  $(LOCAL_DEFS) $(CXX_OPT)

# Extended Cantera Thermodynamics Object Files

ifeq ($(do_electro),1)
ELECTRO_OBJ = SingleSpeciesTP.o StoichSubstanceSSTP.o \
              MolalityVPSSTP.o  VPStandardStateTP.o \
              IdealSolidSolnPhase.o IdealMolalSoln.o \
              WaterPropsIAPWSphi.o WaterPropsIAPWS.o WaterProps.o \
              PDSS.o WaterPDSS.o WaterTP.o \
	      HMWSoln.o HMWSoln_input.o DebyeHuckel.o

ELECTRO_H =  SingleSpeciesTP.h StoichSubstanceSSTP.h \
               MolalityVPSSTP.h  VPStandardStateTP.h \
               IdealSolidSolnPhase.h IdealMolalSoln.h \
               WaterPropsIAPWSphi.h WaterPropsIAPWS.h WaterProps.h \
               PDSS.h WaterPDSS.h WaterTP.h HMWSoln.h electrolytes.h \
	       DebyeHuckel.h
endif
ifeq ($(do_issp),1)
ISSP_OBJ = IdealSolidSolnPhase.o
ISSP_H   = IdealSolidSolnPhase.h
endif

CATHERMO_OBJ = $(ELECTRO_OBJ) $(ISSP_OBJ) 

CATHERMO_H   = $(ELECTRO_H) $(ISSP_H) 


CXX_INCLUDES = -I.. 
LIB = /home/groth/CFFC/cantera/build/lib/i686-pc-linux-gnu/libcaThermo.a

DEPENDS = $(CATHERMO_OBJ:.o=.d)

all: $(LIB) .depends
	@(/usr/bin/install -c -d $(INCDIR))
	@(for lh in $(CATHERMO_H) ; do \
	  $(INSTALL_TSC) "$${lh}" $(INCDIR) ; \
	done)
%.d:
	g++ -MM $(CXX_INCLUDES) $*.cpp > $*.d

.cpp.o:
	g++ -O -c $< $(CXX_FLAGS) $(CXX_INCLUDES)

$(LIB): $(CATHERMO_OBJ) $(CATHERMO_H)
	ar ruv $(LIB) $(CATHERMO_OBJ) > /dev/null
ifeq ($(do_ranlib),1)
	ranlib $(LIB)
endif

clean:
	@(for lh in dummy.h $(CATHERMO_H) ; do \
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
	$(RM) *.o *~  .depends *.d
	(if test -d SunWS_cache ; then \
           $(RM) -rf SunWS_cache ; \
         fi )

depends:
	make .depends

.depends: $(DEPENDS)
	cat $(DEPENDS) > .depends

TAGS: 
	etags *.h *.cpp

ifeq ($(wildcard .depends), .depends)
include .depends
endif


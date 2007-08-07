#/bin/sh
###############################################################
#  $Author: dggoodwin $
#  $Date: 2006/11/11 13:23:11 $
#  $Revision: 1.49 $
#
#  Copyright 2001 California Institute of Technology
#
###############################################################

.SUFFIXES :
.SUFFIXES : .cpp .d .o .dh .h .h.gch

INCDIR      = ../../build/include/cantera/kernel
INSTALL_TSC = ../../bin/install_tsc
CANTERA_LIB = /nfs/kris/d1/people/charestm/Programs/CFFC/cantera/build/lib/i686-pc-linux-gnu/libcantera.a
#LOCAL_DEFNS=-DDEBUG_PATHS -DDEBUG_MULTIPHASE_EQUIL
CXX_FLAGS   = -D_NO_MPI_VERSION -D_GNU_GCC_3 -D_NO_ICEMCFD_VERSION -I/nfs/kris/d1/people/charestm/Programs/CFFC/sparselib++/include -I/nfs/kris/d1/people/charestm/Programs/CFFC/sparselib++/iml/include -I/nfs/kris/d1/people/charestm/Programs/CFFC/sparselib++/mv/include -I/nfs/kris/d1/people/charestm/Programs/CFFC/bpkit/src  $(CXX_OPT) $(LOCAL_DEFNS)
EXT         = ../../ext
do_ranlib = 1
USE_SUNDIALS = 0
SUNDIALS_INC = 

#----------------------
#  kernel components
#----------------------

# basic components always needed
BASE_OBJ   = State.o Elements.o Constituents.o stringUtils.o misc.o \
             importCTML.o plots.o \
             xml.o Phase.o DenseMatrix.o ctml.o funcs.o \
             phasereport.o ct2ctml.o

BASE_H     = State.h Elements.h Constituents.h stringUtils.h global.h \
             importCTML.h plots.h xml.h Phase.h DenseMatrix.h ctml.h \
             funcs.h ct_defs.h Array.h PropertyCalculator.h \
             ctexceptions.h ctlapack.h units.h mix_defs.h\
             utilities.h vec_functions.h XML_Writer.h diagnostics.h \
             config.h logger.h
BASE       = $(BASE_OBJ)

# thermodynamic properties
THERMO_OBJ = ThermoPhase.o IdealGasPhase.o ConstDensityThermo.o \
             SpeciesThermoFactory.o ConstCpPoly.o \
             Mu0Poly.o GeneralSpeciesThermo.o \
             ThermoFactory.o  StoichSubstance.o PureFluidPhase.o LatticeSolidPhase.o LatticePhase.o

THERMO_H   = ThermoPhase.h IdealGasPhase.h ConstDensityThermo.h \
             SpeciesThermoFactory.h ThermoFactory.h \
             NasaPoly1.h NasaPoly2.h NasaThermo.h \
             ShomateThermo.h ShomatePoly.h ConstCpPoly.h\
             SimpleThermo.h SpeciesThermoMgr.h \
             SpeciesThermoInterpType.h \
             GeneralSpeciesThermo.h Mu0Poly.h \
             speciesThermoTypes.h SpeciesThermo.h SurfPhase.h \
             EdgePhase.h polyfit.h  Func1.h \
             FuncEval.h StoichManager.h  MetalPhase.h StoichSubstance.h PureFluidPhase.h LatticeSolidPhase.h LatticePhase.h
THERMO     = $(BASE_OBJ) $(THERMO_OBJ) $(THERMO_H)

# homogeneous kinetics
KINETICS_OBJ=GRI_30_Kinetics.o KineticsFactory.o GasKinetics.o \
             FalloffFactory.o ReactionStoichMgr.o Kinetics.o
KINETICS_H = GRI_30_Kinetics.h KineticsFactory.h Kinetics.h GasKinetics.h \
             FalloffFactory.h  ReactionStoichMgr.h reaction_defs.h \
             FalloffMgr.h ThirdBodyMgr.h RateCoeffMgr.h ReactionData.h \
             RxnRates.h Enhanced3BConc.h
KINETICS   = $(THERMO) $(KINETICS_OBJ) $(KINETICS_H)


# heterogeneous kinetics
HETEROKIN_OBJ=InterfaceKinetics.o ImplicitSurfChem.o SurfPhase.o \
              EdgeKinetics.o
HETEROKIN_H  =InterfaceKinetics.h ImplicitSurfChem.h SurfPhase.h \
              EdgeKinetics.h
HETEROKIN    = $(THERMO) $(HETEROKIN_OBJ)

# support for importing from Chemkin-compatible reaction mechanisms
CK         = $(KINETICS)

# chemical equilibrium
EQUIL_OBJ  = ChemEquil.o MultiPhaseEquil.o sort.o MultiPhase.o equilibrate.o
EQUIL_H    = ChemEquil.h MultiPhaseEquil.h MultiPhase.h Crystal.h sort.h equil.h
EQUIL      = $(THERMO) $(EQUIL_OBJ)

# reaction path analysis
RPATH_OBJ  = Group.o ReactionPath.o
RPATH_H    = Group.h ReactionPath.h
RPATH      = $(RPATH_OBJ)


# solvers
SOLVERS_OBJ = ODE_integrators.o  BandMatrix.o DAE_solvers.o
SOLVERS_H   = BandMatrix.h Integrator.h DAE_Solver.h ResidEval.h
SOLVERS    = $(SOLVERS_OBJ)

# 1D flow models
FLOW1D     =  $(KINETICS) $(SOLVERS)

ALL_OBJ = $(BASE_OBJ) $(THERMO_OBJ) $(KINETICS_OBJ) $(HETEROKIN_OBJ) \
          $(EQUIL_OBJ) $(RPATH_OBJ) $(SOLVERS_OBJ)


PCH = 
# config.h ct_defs.h utilities.h ThermoPhase.h Kinetics.h ReactionData.h \
# RateCoeffMgr.h ReactionStoichMgr.h

PCHGCH = $(PCH:.h=.h.gch)

all:  .depends $(PCHGCH) base thermo kinetics ck trprops equil reactor solvers flow1D rpath tpx $(CANTERA_LIB) cathermo

%.h.gch : %.h
ifeq (no,yes)
	g++ -O $*.h $(CXX_FLAGS) 
else
	@echo 'skipping precompiling header file $*.h'
endif

base: $(BASE)
	@(for lh in $(BASE_H) ; do \
             $(INSTALL_TSC) "$${lh}" $(INCDIR) ; \
          done)

thermo: $(THERMO)
	@(for lh in $(THERMO_H) ; do \
             $(INSTALL_TSC) "$${lh}" $(INCDIR) ; \
          done)

cathermo:
ifeq (, 1)
	cd thermo; make all
endif

kinetics: $(KINETICS) $(HETEROKIN)
	@(for lh in $(KINETICS_H) ; do \
             $(INSTALL_TSC) "$${lh}" $(INCDIR) ; \
          done)
	@(for lh in $(HETEROKIN_H) ; do \
             $(INSTALL_TSC) "$${lh}" $(INCDIR) ; \
          done)

ck: 
	cd converters; make

equil: $(EQUIL)
	@(for lh in $(EQUIL_H) ; do \
             $(INSTALL_TSC) "$${lh}" $(INCDIR) ; \
          done)

flow1d: $(FLOW1D)

trprops:
	cd transport; make

reactor:
	cd zeroD; make

rpath: $(RPATH)
	@(for lh in $(RPATH_H) ; do \
             $(INSTALL_TSC) "$${lh}" $(INCDIR) ; \
          done)

solvers: $(SOLVERS)
	@(for lh in $(SOLVERS_H) ; do \
             $(INSTALL_TSC) "$${lh}" $(INCDIR) ; \
          done)

tpx: 
	cd $(EXT)/tpx; make

flow1D:
	cd oneD;  make

CXX_LIBS = 
CXX_INCLUDES =  -I.
CANTERA_LIB = /nfs/kris/d1/people/charestm/Programs/CFFC/cantera/build/lib/i686-pc-linux-gnu/libcantera.a

DEPENDS = $(ALL_OBJ:.o=.d)
SRCS = $(ALL_OBJ:.o=.cpp)
ALL_H = $(BASE_H) $(THERMO_H) $(KINETICS_H) $(HETEROKIN_H) \
        $(EQUIL_H) $(SOLVERS_H) $(RPATH_H)

.cpp.d:
	g++ -MM $(CXX_INCLUDES) $*.cpp > $*.d

.cpp.o:
	g++ -O -c $< $(CXX_INCLUDES) $(CXX_FLAGS) 

ODE_integrators.o:
	g++ -O -c ODE_integrators.cpp $(CXX_INCLUDES) $(SUNDIALS_INC) \
            $(CXX_FLAGS) 
ODE_integrators.d: ODE_integrators.cpp Makefile
	g++ -MM ODE_integrators.cpp $(CXX_INCLUDES) $(SUNDIALS_INC) \
            $(CXX_FLAGS) > ODE_integrators.d 

DAE_solvers.o: DAE_solvers.cpp
	g++ -O -c DAE_solvers.cpp $(CXX_INCLUDES) $(SUNDIALS_INC) $(CXX_FLAGS)
DAE_solvers.d: DAE_solvers.cpp Makefile
	g++ -MM DAE_solvers.cpp $(CXX_INCLUDES) $(SUNDIALS_INC) \
            $(CXX_FLAGS) > DAE_solvers.d  

lib:  $(OBJ_LIB)
	$(RM)     $(CANTERA_LIB)
	ar ruv $(CANTERA_LIB) *.o > /dev/null
ifeq ($(do_ranlib),1)
	ranlib $(CANTERA_LIB)
endif

	cd ../.. ; make hdr-collect

$(CANTERA_LIB): $(BASE_OBJ) $(THERMO_OBJ) $(KINETICS_OBJ) $(HETEROKIN_OBJ) $(EQUIL_OBJ) $(SOLVERS_OBJ) $(RPATH_OBJ)
	$(RM)     $(CANTERA_LIB)
	ar ruv $(CANTERA_LIB) $(BASE_OBJ) $(THERMO_OBJ) $(KINETICS_OBJ) $(HETEROKIN_OBJ) $(EQUIL_OBJ) $(SOLVERS_OBJ) $(RPATH_OBJ) 
ifeq ($(do_ranlib),1)
	ranlib $(CANTERA_LIB)
endif

clean:
	$(RM) *.o *.gch *~ $(CANTERA_LIB) *.d
	@(for lh in $(ALL_H) ; do \
             th=$(INCDIR)/$$lh ; \
             if test -f $$th ; then \
               $(RM) $$th ; \
             fi \
          done)
	$(RM) *.o *.gch *~ $(CANTERA_LIB) .depends
	(if test -d SunWS_cache ; then \
           $(RM) -rf SunWS_cache ; \
         fi )
ifeq (1, 1)
	cd zeroD; make clean
endif
	cd oneD;  make clean
	cd converters; make clean
	cd transport; make clean
	cd thermo; make clean

depends: $(DEPENDS)
	cat *.d > .depends
	cd oneD; make depends
ifeq (1, 1)
	cd zeroD; make .depends
endif
	cd converters; make depends
	cd transport; make depends
ifeq (, 1)
	cd thermo; make depends
endif
	$(RM) $(DEPENDS) 

.depends: $(DEPENDS) Makefile
	cat *.d > .depends


TAGS: 
	etags *.h *.cpp

ifeq ($(wildcard .depends), .depends)
include .depends
endif

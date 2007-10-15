#
# Makefile for CFFC (3D version)
#
# Computational Framework for Fluids and Combustion (CFFC)
#

#
# Set CFFC directory location which can be overridden by 
# setting an environement variable CFFC_Path
CFFC_Path:= $(shell cd ../; pwd)

#
# Set CFFC source tree
CFFC_SRC_TREE = 3D

#
# Set standard directories
CFFC_BIN_DIR = $(CFFC_Path)/bin

#
# Platform compile options 
include makefile.def

#
# Assign sub directories
DEPEND_DIR = Depend/
AMR_DIR = AMR/
CFD_DIR = CFD/
EULER_DIR = Euler/
FANS_DIR = FANS/
GRID_DIR = Grid/
HEXA_DIR = HexaBlock/
ICEMCFD_DIR = ICEMCFD/
MATH_DIR = Math/
MPI_DIR = MPI/
MULTIGRID_DIR = FASMultigrid/
NAVIERSTOKES_DIR = NavierStokes/
NKS_DIR = NewtonKrylovSchwarz/
PHYSICS_DIR = Physics/
EOS_DIR = $(PHYSICS_DIR)EquationOfState/
NASADATA_DIR = $(PHYSICS_DIR)NASAData/
REACTIONS_DIR = Reactions/
SYSTEM_DIR = System/
TURBULENCE_DIR = TurbulenceModelling/

#
# Set default f90 module directory
MODULE_DIR = F90_Modules

###############################################################################
# C++, C, and Fortran Source Files
###############################################################################
# Math:
SRC_MATH = $(MATH_DIR)Vector2D.cc \
           $(MATH_DIR)Tensor2D.cc \
           $(MATH_DIR)Third_order_tensor2D.cc \
	   $(MATH_DIR)Tensor3D.cc \
           $(MATH_DIR)Spline2D.cc \
           $(MATH_DIR)Matrix.cc \
           $(MATH_DIR)LinearSystems.cc \
	   $(MATH_DIR)LinkedList.cc \
	   $(MATH_DIR)Polygon.cc \
	   $(MATH_DIR)Polyfit/Polyfit.cc \
	   $(MATH_DIR)Polyfit/dp1vlu.f \
	   $(MATH_DIR)Polyfit/dpcoef.f \
	   $(MATH_DIR)Polyfit/dpolft.f \
	   $(MATH_DIR)NumericalLibrary.cc
SRC_MATH_F90 = $(MATH_DIR)complexify.f90
SRC_MATH += $(SRC_MATH_F90)
SRC_MATH += Utilities/EpsilonTol.cc

# Add Math Testing If Required:
SRC_MATH_TESTS = TestFunctions/TestFunctions_1D.cc \
	         TestFunctions/TestFunctions_2D.cc \
	         TestFunctions/TestFunctions_3D.cc \
	         test_LeastSquaresSolvers.cc \
	         test_NumericalLibrary.cc
SRC_MATH_TESTS := $(addprefix $(MATH_DIR)UnitTests/, $(SRC_MATH_TESTS))
ifeq ($(TUT_TESTING),ON) 
  SRC_MATH += $(SRC_MATH_TESTS)
endif

# CFD and MPI:
SRC_CFD = $(CFD_DIR)CFD.cc \
          $(CFD_DIR)CFDInput.cc
SRC_MPI	= $(MPI_DIR)MPI.cc

# Grid:
SRC_GRID = $(GRID_DIR)Cell2D.cc \
           $(GRID_DIR)Grid2DQuadBlock.cc \
           $(GRID_DIR)Grid2DQuadMultiBlock.cc \
           $(GRID_DIR)Grid3DInput.cc \
           $(GRID_DIR)Grid3DHexaBlock.cc \
           $(GRID_DIR)Grid3DHexaMultiBlock.cc
SRC_ICEMCFD = $(ICEMCFD_DIR)ICEMCFD.cc \
	      $(ICEMCFD_DIR)ICEMCFD_ctype.c

# AMR:
ifeq ($(MPI_VERSION),NoMPI)
  SRC_AMR = $(AMR_DIR)AMRInput.cc \
            $(AMR_DIR)AdaptiveBlock3D.cc \
	    $(AMR_DIR)AdaptiveBlock3D_NoMPI.cc \
            $(AMR_DIR)Octree.cc
else
  SRC_AMR = $(AMR_DIR)AMRInput.cc \
            $(AMR_DIR)AdaptiveBlock3D.cc \
	    $(AMR_DIR)AdaptiveBlock3D_MPI.cc \
            $(AMR_DIR)Octree.cc
endif

# Multigrid:
SRC_MULTIGRID = $(MULTIGRID_DIR)FASMultigridInput.cc

# NKS:
SRC_NKS = $(NKS_DIR)NKSInput.cc

# Euler3D:
SRC_EULER_POLYTROPIC = $(EULER_DIR)Euler3DPolytropicState.cc

SRC_EULER_THERMPERFECT = $(EULER_DIR)Euler3DThermallyPerfectState.cc \
                         $(EULER_DIR)Euler3DThermallyPerfectHexaBlock.cc

# Add Euler3D testing if required:
SRC_EULER_POLYTROPIC_TESTS = test_Euler3DPolytropicState.cc
SRC_EULER_POLYTROPIC_TESTS := $(addprefix $(EULER_DIR)UnitTests/, $(SRC_EULER_POLYTROPIC_TESTS))
ifeq ($(TUT_TESTING),ON) 
  SRC_EULER_POLYTROPIC += $(SRC_EULER_POLYTROPIC_TESTS)
endif

SRC_EULER = $(SRC_EULER_POLYTROPIC) $(SRC_EULER_THERMPERFECT)

SRC_EULER_POLYTROPIC_EXE = $(EULER_DIR)Euler3DPolytropic.cc
SRC_EULER_THERMPERFECT_EXE = $(EULER_DIR)Euler3DThermallyPerfect.cc

# NavierStokes3D:
SRC_NAVIERSTOKES_POLYTROPIC =
SRC_NAVIERSTOKES_THERMPERFECT = $(NAVIERSTOKES_DIR)NavierStokes3DThermallyPerfectState.cc \
                                $(NAVIERSTOKES_DIR)NavierStokes3DThermallyPerfectHexaBlock.cc

SRC_NAVIERSTOKES = $(SRC_NAVIERSTOKES_POLYTROPIC) $(SRC_NAVIERSTOKES_THERMPERFECT)

SRC_NAVIERSTOKES_POLYTROPIC_EXE = $(NAVIERSTOKES_DIR)NavierStokes3DPolytropic.cc
SRC_NAVIERSTOKES_THERMPERFECT_EXE = $(NAVIERSTOKES_DIR)NavierStokes3DThermallyPerfect.cc

# FANS3D:
SRC_FANS_POLYTROPIC =
SRC_FANS_THERMPERFECT = $(FANS_DIR)FANS3DThermallyPerfectState.cc \
                        $(FANS_DIR)FANS3DThermallyPerfectHexaBlock.cc

SRC_FANS = $(SRC_FANS_POLYTROPIC) $(SRC_FANS_THERMPERFECT)

SRC_FANS_POLYTROPIC_EXE = $(FANS_DIR)FANS3DPolytropic.cc
SRC_FANS_THERMPERFECT_EXE = $(FANS_DIR)FANS3DThermallyPerfect.cc

# Physics, Turbulence, Thermodynamics, and Chemistry:
SRC_PHYSICS = $(PHYSICS_DIR)SpeciesInput.cc
SRC_TURBULENCE = $(TURBULENCE_DIR)TurbulenceModellingInput.cc \
                 $(TURBULENCE_DIR)TurbulenceModelling.cc 
SRC_EOS = $(EOS_DIR)tgas.f \
          $(EOS_DIR)ugas.f
SRC_NASADATA = $(NASADATA_DIR)NASARP1311data.cc
SRC_REACTIONS = $(REACTIONS_DIR)Reactions.cc \
	   	$(REACTIONS_DIR)15step_GRI211.f \
	  	$(REACTIONS_DIR)15step_GRI3.f
SRC_REACTIONS_F90 = $(REACTIONS_DIR)c_15step_GRI211.f90 \
                    $(REACTIONS_DIR)c_15step_GRI3.f90
SRC_REACTIONS += $(SRC_REACTIONS_F90)

# TUT testing:
SRC_TUT = TestRun.cc \
          TestData.cc \
	  BasicTests/Sample_Test.cc \
          BasicTests/test_TestData.cc
SRC_TUT := $(addprefix $(UNIT_TESTING_DIR), $(SRC_TUT))

# System:
SRC_SYSTEM = $(SYSTEM_DIR)System_Linux.cc

# Base and CFFC sources:
SRC_BASE = $(SRC_MATH) $(SRC_CFD) $(SRC_MPI) $(SRC_SYSTEM) \
           $(SRC_AMR) $(SRC_GRID) $(SRC_ICEMCFD) \
           $(SRC_MULTIGRID) $(SRC_NKS)
ifeq ($(TUT_TESTING),ON) 
  SRC_BASE += $(SRC_TUT)
endif
SRC_CFFC = $(SRC_BASE) \
           $(SRC_EULER) $(SRC_NAVIERSTOKES) $(SRC_FANS) \
           $(SRC_PHYSICS) $(SRC_TURBULENCE) $(SRC_EOS) $(SRC_NASADATA) \
           $(SRC_REACTIONS)

# Executable sources:
SRC_MAIN = cffc3D.cc

###############################################################################
# Compiled object code files
###############################################################################
# Math:
OBJ_MATH  = $(patsubst %.cc,%.o,$(filter %.cc,$(SRC_MATH)))
OBJ_MATH += $(patsubst %.f,%.o,$(filter %.f,$(SRC_MATH)))
OBJ_MATH += $(patsubst %.f90,%.o,$(filter %.f90,$(SRC_MATH)))
OBJ_MATH_F90 = $(SRC_MATH_F90:.f90=.o)
MOD_MATH_F90 = $(MODULE_DIR)/complexify.mod
OBJ_MATH_TESTS := $(SRC_MATH_TESTS:.cc=.o)

# CFD and MPI:
OBJ_CFD = $(SRC_CFD:.cc=.o)
OBJ_MPI = $(SRC_MPI:.cc=.o)

# Grid:
OBJ_GRID = $(SRC_GRID:.cc=.o)
OBJ_ICEMCFD = $(patsubst %.cc,%.o,$(filter %.cc,$(SRC_ICEMCFD)))
OBJ_ICEMCFD += $(patsubst %.c,%.o,$(filter %.c,$(SRC_ICEMCFD)))

# AMR:
OBJ_AMR = $(SRC_AMR:.cc=.o)

# Multigrid:
OBJ_MULTIGRID = $(SRC_MULTIGRID:.cc=.o)

# NKS:
OBJ_NKS = $(SRC_NKS:.cc=.o)

# Euler3D:
OBJ_EULER_POLYTROPIC = $(SRC_EULER_POLYTROPIC:.cc=.o)
OBJ_EULER_THERMPERFECT = $(SRC_EULER_THERMPERFECT:.cc=.o)
OBJ_EULER = $(SRC_EULER:.cc=.o)
OBJ_EULER_POLYTROPIC_EXE = $(SRC_EULER_POLYTROPIC_EXE:.cc=.o)
OBJ_EULER_THERMPERFECT_EXE = $(SRC_EULER_THERMPERFECT_EXE:.cc=.o)

# NavierStokes3D:
OBJ_NAVIERSTOKES_POLYTROPIC = $(SRC_NAVIERSTOKES_POLYTROPIC:.cc=.o)
OBJ_NAVIERSTOKES_THERMPERFECT = $(SRC_NAVIERSTOKES_THERMPERFECT:.cc=.o)
OBJ_NAVIERSTOKES = $(SRC_NAVIERSTOKES:.cc=.o)
OBJ_NAVIERSTOKES_POLYTROPIC_EXE = $(SRC_NAVIERSTOKES_POLYTROPIC_EXE:.cc=.o)
OBJ_NAVIERSTOKES_THERMPERFECT_EXE = $(SRC_NAVIERSTOKES_THERMPERFECT_EXE:.cc=.o)

# FANS3D:
OBJ_FANS_POLYTROPIC = $(SRC_FANS_POLYTROPIC:.cc=.o)
OBJ_FANS_THERMPERFECT = $(SRC_FANS_THERMPERFECT:.cc=.o)
OBJ_FANS = $(SRC_FANS:.cc=.o)
OBJ_FANS_POLYTROPIC_EXE = $(SRC_FANS_POLYTROPIC_EXE:.cc=.o)
OBJ_FANS_THERMPERFECT_EXE = $(SRC_FANS_THERMPERFECT_EXE:.cc=.o)

# Physics, Turbulence, Thermodynamics, and Chemistry:
OBJ_PHYSICS = $(SRC_PHYSICS:.cc=.o)
OBJ_TURBULENCE = $(SRC_TURBULENCE:.cc=.o)
OBJ_EOS = $(SRC_EOS:.f=.o)
OBJ_NASADATA = $(SRC_NASADATA:.cc=.o)
OBJ_REACTIONS = $(patsubst %.cc,%.o,$(filter %.cc,$(SRC_REACTIONS)))
OBJ_REACTIONS += $(patsubst %.f,%.o,$(filter %.f,$(SRC_REACTIONS)))
OBJ_REACTIONS += $(patsubst %.f90,%.o,$(filter %.f90,$(SRC_REACTIONS)))
OBJ_REACTIONS_F90 = $(SRC_REACTIONS_F90:.f90=.o)

# TUT testing:
OBJ_TUT := $(SRC_TUT:.cc=.o)

# System:
OBJ_SYSTEM = $(SRC_SYSTEM:.cc=.o)

# Base and CFFC objects:
OBJ_BASE = $(OBJ_MATH) $(OBJ_CFD) $(OBJ_MPI) $(OBJ_SYSTEM) \
           $(OBJ_AMR) $(OBJ_GRID) $(OBJ_ICEMCFD) \
           $(OBJ_MULTIGRID) $(OBJ_NKS)
ifeq ($(TUT_TESTING),ON)
  OBJ_BASE += $(OBJ_TUT)
endif
OBJ_CFFC = $(OBJ_BASE) \
           $(OBJ_EULER) $(OBJ_NAVIERSTOKES) $(OBJ_FANS) \
           $(OBJ_PHYSICS) $(OBJ_TURBULENCE) $(OBJ_EOS) $(OBJ_NASADATA) \
           $(OBJ_REACTIONS)

# Executable objects:
OBJ_MAIN = $(SRC_MAIN:.cc=.o)

###############################################################################
# Define the executables
###############################################################################
EXE_MAIN = cffc3D
EXE_EULER_POLYTROPIC = euler3Dpolytropic
EXE_EULER_THERMPERFECT = euler3Dthermallyperfect
EXE_NAVIERSTOKES_THERMPERFECT = navierstokes3Dthermallyperfect
EXE_FANS_THERMPERFECT = fans3Dthermallyperfect

###############################################################################
# Define dependencies
###############################################################################
help:
	@echo "+-----------------------------------------------------------------------------------+"
	@echo "|                                                                                   |"
	@echo "|                           makefile for CFFC                                       |"
	@echo "|                                                                                   |"
	@echo "| Usage: make all                  compile all executable programs                  |"
	@echo "|        make cffc3D               compile cffc3D program                           |"
	@echo "|        make euler3Dpolytropic    compile euler3Dpolytropic stand alone program    |"
	@echo "|        make euler3Dthermallyperfect                                               |"
	@echo "|             compile euler3Dthermallyperfect stand alone program                   |"
	@echo "|        make navierstokes3Dthermallyperfect                                        |"
	@echo "|             compile navierstokes3Dthermallyperfect stand alone program            |"
	@echo "|        make fans3Dthermallyperfect                                                |"
	@echo "|             compile fans3Dthermallyperfect stand alone program                    |"
	@echo "|        make install              copies cffc3D executable to bin directory        |"
	@echo "|        make install-all          copies all executables to bin directory          |"
	@echo "|        make uninstall            removes cffc3D executable from bin directory     |"
	@echo "|        make uinstall-all         removes all executables from bin directory       |"
	@echo "|        make documentation        create code documentaion using doxygen & latex   |"
	@echo "|        make libraries-all        create all supporting libraries                  |"
	@echo "|        make sparselib            create SparseLib++ library                       |"
	@echo "|        make bpkit                create BPKIT library                             |"
	@echo "|        make blkclib              create block connectivity library                |"
	@echo "|        make cantera              create Cantera library                           |"
	@echo "|        make fftw                 create fftw library                              |"
	@echo "|        make clean                clean up source directories                      |"
	@echo "|        make superclean           clean up documentation and library directories   |"
	@echo "|        make help                 display makefile usage information               |"
	@echo "|                                                                                   |"
	@echo "| Package Options: MPI_VERSION = NoMPI (default), MPICH, MPICH2, MPT (sgi)          |"
	@echo "|                  ICEMCFD_VERSION = NoICEMCFD (default),V41,V42,V43,V50,V10,V11    |"
	@echo "|                  GCC_VERSION = V296, V3+ (default), ICC                           |"
	@echo "|                  CANTERA_VERSION = NoCANTERA (default), V70                       |"
	@echo "|                  FFTW_VERSION = V32 (default), NoFFTW                             |"
	@echo "|                  TUT_TESTING = ON (default), OFF                                  |"
	@echo "|                  PROFILING = OFF (default), ON                                    |"
	@echo "|                                                                                   |"
	@echo "| Ensure that the source directory is clean and the ICEMCFD, SparseLib++, BPKIT,    |"
	@echo "| FFTW, and Cantera libraries are up to date before re-making CFFC from scratch.    |"
	@echo "|                                                                                   |"
	@echo "| An environment variable, CFFC_Path, can to be set to the location                 |"
	@echo "| of the CFFC source directory.  This can be done as follows:                       |"
	@echo "|                                                                                   |"
	@echo "|  bash:   export CFFC_Path=/your/path/here/CFFC                                    |"
	@echo "|  tcsh:   setenv CFFC_Path /your/path/here/CFFC                                    |"
	@echo "|                                                                                   |"
	@echo "+-----------------------------------------------------------------------------------+"

all: $(EXE_MAIN) $(EXE_EULER_POLYTROPIC) $(EXE_EULER_THERMPERFECT) $(EXE_NAVIERSTOKES_THERMPERFECT) $(EXE_FANS_THERMPERFECT)
	@echo ' '	
	@echo All programs have been brought up to date.
	@echo '-------------------------------------------'
	@echo ' '

libraries-all: sparselib bpkit blkclib fftw
	@echo ' '	
	@echo All libraries have been brought up to date.
	@echo '-------------------------------------------'
	@echo ' '

sparselib:
	@echo ' '
	@echo Building $(SPARSELIB) library for a $(PLATFORM) platform.
	@echo '--------------------------------------------------------'
	@echo ' '
	cd $(SPARSELIB_DIR); make CFFC_SRC_TREE=$(CFFC_SRC_TREE) wipe;
	cd $(SPARSELIB_DIR); make CFFC_SRC_TREE=$(CFFC_SRC_TREE) clean;
	cd $(SPARSELIB_DIR); make CFFC_SRC_TREE=$(CFFC_SRC_TREE) sp;
	@echo ' '
	@echo $(SPARSELIB) library has been brought up to date.
	@echo '------------------------------------------------'
	@echo ' '

bpkit:
	@echo ' '
	@echo Building $(BPKIT) library for a $(PLATFORM) platform.
	@echo '----------------------------------------------------'
	@echo ' '
	cd $(BPKIT_DIR)/src; make CFFC_SRC_TREE=$(CFFC_SRC_TREE) clean;
	cd $(BPKIT_DIR)/src; make CFFC_SRC_TREE=$(CFFC_SRC_TREE);
	@echo ' '
	@echo $(BPKIT) library has been brought up to date.
	@echo '--------------------------------------------'
	@echo ' '

blkclib:
	@echo ' '
	@echo Building $(BLKCLIB) library for a $(PLATFORM) platform.
	@echo '--------------------------------------------------------'
	@echo ' '
	cd $(BLKCLIB_DIR); ./configure;	
	cd $(BLKCLIB_DIR); $(MAKE) distclean;
	cd $(BLKCLIB_DIR); ./configure;
	cd $(BLKCLIB_DIR); $(MAKE);	
	@echo ' '
	@echo $(BLKCLIB) library has been brought up to date.
	@echo '------------------------------------------------'
	@echo ' '

cantera:
	@echo ' '
	@echo Building $(CANTERA) library for a $(PLATFORM) platform.
	@echo '----------------------------------------------------'
	@echo ' '
	cd $(CANTERA_DIR); ./preconfig CXX="$(CXX) $(CXXFLAGS)" CC="$(CC)" F77="$(FC)" \
          CXXFLAGS="$(CPPFLAGS)" CFLAGS="$(CFLAGS)" FFLAGS="$(FFLAGS)" \
	  LCXX_END_LIBS="$(CANTERA_LCXX_END_LIBS)" PYTHON_PACKAGE="$(CANTERA_PYTHON_PACKAGE)";
	cd $(CANTERA_DIR); make;
	@echo ' '
	@echo $(CANTERA) library has been brought up to date.
	@echo '--------------------------------------------'
	@echo ' '

fftw:
	@echo ' '
	@echo Building $(FFTW) library for a $(PLATFORM) platform.
	@echo '--------------------------------------------------------'
	@echo ' '
	cd $(FFTW_DIR); ./configure --prefix=$(FFTW_DIR);
	cd $(FFTW_DIR); make;
	cd $(FFTW_DIR); make install;
	@echo ' '
	@echo $(FFTW) library has been brought up to date.
	@echo '------------------------------------------------'
	@echo ' '

$(EXE_MAIN): $(OBJ_MAIN) $(OBJ_CFFC)
	@echo ' '
	@echo Building program $(EXE_MAIN) for a $(PLATFORM) platform.
	@echo '--------------------------------------------------------'
	@echo ' '
	$(LD) -o $(EXE_MAIN) $(OBJ_MAIN) $(OBJ_CFFC) $(LDFLAGS)
	@echo ' '
	@echo Program $(EXE_MAIN) has been brought up to date.
	@echo '------------------------------------------------'
	@echo ' '

$(EXE_EULER_POLYTROPIC): $(OBJ_EULER_POLYTROPIC_EXE) $(OBJ_BASE) $(OBJ_EULER_POLYTROPIC) $(OBJ_PHYSICS) $(OBJ_TURBULENCE) $(OBJ_REACTIONS)
	@echo ' '
	@echo Building program $(EXE_EULER_POLYTROPIC) for a $(PLATFORM) platform.
	@echo '-----------------------------------------------------------------'
	@echo ' '
	$(LD) -o $(EXE_EULER_POLYTROPIC) $(OBJ_EULER_POLYTROPIC_EXE) $(OBJ_BASE) $(OBJ_EULER_POLYTROPIC) $(OBJ_PHYSICS) $(OBJ_TURBULENCE)  $(OBJ_REACTIONS) $(LDFLAGS)
	@echo ' '
	@echo Program $(EXE_EULER_POLYTROPIC) has been brought up to date.
	@echo '---------------------------------------------------------'
	@echo ' '

$(EXE_EULER_THERMPERFECT): $(OBJ_EULER_THERMPERFECT_EXE) $(OBJ_BASE) $(OBJ_EULER_THERMPERFECT) $(OBJ_PHYSICS) $(OBJ_TURBULENCE) $(OBJ_NASADATA) $(OBJ_REACTIONS)
	@echo ' '
	@echo Building program $(EXE_EULER_THERMPERFECT) for a $(PLATFORM) platform.
	@echo '-----------------------------------------------------------------'
	@echo ' '
	$(LD) -o $(EXE_EULER_THERMPERFECT) $(OBJ_EULER_THERMPERFECT_EXE) $(OBJ_BASE) $(OBJ_EULER_THERMPERFECT) $(OBJ_PHYSICS) $(OBJ_TURBULENCE) $(OBJ_NASADATA) $(OBJ_REACTIONS) $(LDFLAGS)
	@echo ' '
	@echo Program $(EXE_EULER_THERMPERFECT) has been brought up to date.
	@echo '---------------------------------------------------------'
	@echo ' '

$(EXE_NAVIERSTOKES_THERMPERFECT): $(OBJ_NAVIERSTOKES_THERMPERFECT_EXE) $(OBJ_BASE) $(OBJ_NAVIERSTOKES_THERMPERFECT) $(OBJ_EULER_THERMPERFECT) $(OBJ_PHYSICS) $(OBJ_TURBULENCE) $(OBJ_NASADATA) $(OBJ_REACTIONS)
	@echo ' '
	@echo Building program $(EXE_NAVIERSTOKES_THERMPERFECT) for a $(PLATFORM) platform.
	@echo '-----------------------------------------------------------------'
	@echo ' '
	$(LD) -o $(EXE_NAVIERSTOKES_THERMPERFECT) $(OBJ_NAVIERSTOKES_THERMPERFECT_EXE) $(OBJ_BASE) $(OBJ_NAVIERSTOKES_THERMPERFECT) $(OBJ_EULER_THERMPERFECT) $(OBJ_PHYSICS) $(OBJ_TURBULENCE) $(OBJ_NASADATA) $(OBJ_REACTIONS) $(LDFLAGS)
	@echo ' '
	@echo Program $(EXE_NAVIERSTOKES_THERMPERFECT) has been brought up to date.
	@echo '---------------------------------------------------------'
	@echo ' '

$(EXE_FANS_THERMPERFECT): $(OBJ_FANS_THERMPERFECT_EXE) $(OBJ_BASE) $(OBJ_FANS_THERMPERFECT) $(OBJ_NAVIERSTOKES_THERMPERFECT) $(OBJ_EULER_THERMPERFECT) $(OBJ_PHYSICS) $(OBJ_TURBULENCE) $(OBJ_NASADATA) $(OBJ_REACTIONS)
	@echo ' '
	@echo Building program $(EXE_FANS_THERMPERFECT) for a $(PLATFORM) platform.
	@echo '-----------------------------------------------------------------'
	@echo ' '
	$(LD) -o $(EXE_FANS_THERMPERFECT) $(OBJ_FANS_THERMPERFECT_EXE) $(OBJ_BASE) $(OBJ_FANS_THERMPERFECT) $(OBJ_NAVIERSTOKES_THERMPERFECT) $(OBJ_EULER_THERMPERFECT) $(OBJ_PHYSICS) $(OBJ_TURBULENCE) $(OBJ_NASADATA) $(OBJ_REACTIONS) $(LDFLAGS)
	@echo ' '
	@echo Program $(EXE_FANS_THERMPERFECT) has been brought up to date.
	@echo '---------------------------------------------------------'
	@echo ' '

#
# Install/Uninstall targets
install: $(EXE_MAIN)
	@echo "Installing $(EXE_MAIN)"
	cp $(EXE_MAIN) $(CFFC_BIN_DIR)
	@echo " "
	@echo "Installation complete"
	@echo "---------------------"
	@echo " "

install-all: $(EXE_MAIN) $(EXE_EULER_POLYTROPIC) $(EXE_EULER_THERMPERFECT) \
             $(EXE_NAVIERSTOKES_THERMPERFECT) $(EXE_FANS_THERMPERFECT)
	@for f in $(EXE_MAIN) $(EXE_EULER_POLYTROPIC) $(EXE_EULER_THERMPERFECT) \
                  $(EXE_NAVIERSTOKES_THERMPERFECT) $(EXE_FANS_THERMPERFECT); \
	do \
	  echo "Installing $$f"; \
	  cp $$f $(CFFC_BIN_DIR); \
	done;
	@echo " "
	@echo "Installation complete"
	@echo "---------------------"
	@echo " "

uninstall:
	@echo "Uninstalling $(EXE_MAIN)"
	cd $(CFFC_BIN_DIR); rm -f $(EXE_MAIN)
	@echo " "
	@echo "Un-Install complete"
	@echo "---------------------"
	@echo " "

uninstall-all:
	@for f in $(EXE_MAIN) $(EXE_EULER_POLYTROPIC) $(EXE_EULER_THERMPERFECT) \
                  $(EXE_NAVIERSTOKES_THERMPERFECT) $(EXE_FANS_THERMPERFECT); \
	do \
	  echo "Uninstalling $$f"; \
	  cd $(CFFC_BIN_DIR) $(EXE_EULER_POLYTROPIC) $(EXE_EULER_THERMPERFECT) \
             $(EXE_NAVIERSTOKES_THERMPERFECT) $(EXE_FANS_THERMPERFECT); \
	  rm -f $$f; \
	done;
	@echo " "
	@echo "Un-Install complete"
	@echo "---------------------"
	@echo " "

#
# Create code documentation using doxygen and latex.
documentation:
	@echo ' '
	@echo Creating code documentation using doxygen.
	@echo '------------------------------------------'
	@echo ' '
	cp $(DOXYGEN_DIR_3D)/code_documentation.doxygen .
	doxygen code_documentation.doxygen
	rm code_documentation.doxygen
	@echo ' '
	@echo Code documentation created using doxygen.
	@echo '-----------------------------------------'
	@echo ' '
	@echo Creating user manuals using LaTeX.
	@echo '------------------------------------------'
	@echo ' '
	cd $(LATEX_DIR_3D)/users_guide_src; make.userguide;
	cd $(LATEX_DIR_3D)/programmers_guide_src; make.programmersguide; 
	cd $(LATEX_DIR_3D)/reference_manual_src; make.referencemanual;
	@echo ' '
	@echo User manuals created using LaTeX.
	@echo '-----------------------------------------'
	@echo ' '

# Clean up source directory
clean:
	@echo ' '
	@echo Cleaning source file directories.
	@echo '---------------------------------'
	@echo ' '
	rm -f *.depend *.o *.bak *.il *~
	rm -f $(MATH_DIR)*.o $(MATH_DIR)*.il $(MATH_DIR)*~
	rm -f $(MATH_DIR)Polyfit/*.o $(MATH_DIR)Polyfit/*.il $(MATH_DIR)Polyfit/*~
	rm -f $(MATH_DIR)UnitTests/*.o $(MATH_DIR)UnitTests/*.il $(MATH_DIR)UnitTests/*~
	rm -f $(MATH_DIR)UnitTests/TestFunctions/*.o $(MATH_DIR)UnitTests/TestFunctions/*.il $(MATH_DIR)UnitTests/TestFunctions/*~
	rm -f $(CFD_DIR)*.o $(CFD_DIR)*.il $(CFD_DIR)*~
	rm -f $(MPI_DIR)*.o $(MPI_DIR)*.il $(MPI_DIR)*~
	rm -f $(GRID_DIR)*.o $(GRID_DIR)*.il $(GRID_DIR)*~
	rm -f $(ICEMCFD_DIR)*.o $(ICEMCFD_DIR)*.il $(ICEMCFD_DIR)*~ $(ICEMCFD_DIR)*.depend $(ICEMCFD_DIR)*.bak
	rm -f $(AMR_DIR)*.o $(AMR_DIR)*.il $(AMR_DIR)*~
	rm -f $(EULER_DIR)*.o $(EULER_DIR)*.il $(EULER_DIR)*~ $(EULER_DIR)*.depend $(EULER_DIR)*.bak
	rm -f $(EULER_DIR)UnitTests/*.o $(EULER_DIR)UnitTests/*.il $(EULER_DIR)UnitTests/*~ $(EULER_DIR)UnitTests/*.bak
	rm -f $(NAVIERSTOKES_DIR)*.o $(NAVIERSTOKES_DIR)*.il $(NAVIERSTOKES_DIR)*~ $(NAVIERSTOKES_DIR)*.depend $(NAVIERSTOKES_DIR)*.bak
	rm -f $(FANS_DIR)*.o $(FANS_DIR)*.il $(FANS_DIR)*~ $(FANS_DIR)*.depend $(NAVIERSTOKES_DIR)*.bak
	rm -f $(HEXA_DIR)*.o $(HEXA_DIR)*.il $(HEXA_DIR)*~
	rm -f $(PHYSICS_DIR)*.o $(PHYSICS_DIR)*.il $(PHYSICS_DIR)*~
	rm -f $(TURBULENCE_DIR)*.o $(TURBULENCE_DIR)*.il $(TURBULENCE_DIR)*~
	rm -f $(EOS_DIR)*.o $(EOS_DIR)*.il $(EOS_DIR)*~
	rm -f $(NASADATA_DIR)*.o $(NASADATA_DIR)*.il $(NASADATA_DIR)*~
	rm -f $(REACTIONS_DIR)*.o $(REACTIONS_DIR)*.il $(REACTIONS_DIR)*~
	rm -f $(UNIT_TESTING_DIR)*.o $(UNIT_TESTING_DIR)*.il $(UNIT_TESTING_DIR)*~
	rm -f $(UNIT_TESTING_DIR)BasicTests/*.o $(UNIT_TESTING_DIR)BasicTests/*.il $(UNIT_TESTING_DIR)BasicTests/*~
	rm -f Utilities/*.o Utilities/.il Utilities/*~
	rm -f $(MULTIGRID_DIR)*.o $(MULTIGRID_DIR)*.il $(MULTIGRID_DIR)*~
	rm -f $(NKS_DIR)*.o $(NKS_DIR)*.il $(NKS_DIR)*~
	rm -f $(SYSTEM_DIR)*.o $(SYSTEM_DIR)*.il $(SYSTEM_DIR)*~
	rm -f $(MODULE_DIR)/*.mod $(MODULE_DIR)/*~
	rm -f $(EXE_MAIN)
	rm -f $(EXE_EULER_POLYTROPIC)
	rm -f $(EXE_EULER_THERMPERFECT)
	rm -f $(EXE_NAVIERSTOKES_THERMPERFECT)
	rm -f $(EXE_FANS_THERMPERFECT)
	rm -fr $(DEPEND_DIR)
	@echo ' '
	@echo All source file directories are now clean.
	@echo '------------------------------------------'
	@echo ' '

clean-extras:
	@echo ' '
	@echo Cleaning cantera and fftw directories.
	@echo '-----------------------------------------------'
	@echo ' '
	@if [ ! -f $(CANTERA_DIR)/Makefile ] ; then \
	  echo "Cantera already clean."; \
	else \
	  cd $(CANTERA_DIR); make -i clean; \
	fi
	@if [ ! -f $(FFTW_DIR)/Makefile ] ; then \
	  echo "fftw already clean."; \
	else \
	  cd $(FFTW_DIR); make clean; \
	  cd $(FFTW_DIR); make distclean; \
	fi

superclean: clean-extras
	@echo ' '
	@echo Cleaning documentation and library directories.
	@echo '-----------------------------------------------'
	@echo ' '
	rm -f $(DOXYGEN_DIR_3D)/html/*
	rm -f $(DOXYGEN_DIR_3D)/latex/*
	cd $(LATEX_DIR_3D)/users_guide_src; rm -f *~ *.ps *.pdf *.aux *.log *.dvi *.toc users_guide/*; 
	cd $(LATEX_DIR_3D)/programmers_guide_src; rm -f *~ *.ps *.pdf *.aux *.log *.dvi *.toc programmers_guide/*; 
	cd $(LATEX_DIR_3D)/reference_manual_src; rm -f *~ *.ps *.pdf *.aux *.log *.dvi *.toc reference_manual/*;
	cd $(SPARSELIB_DIR); make CFFC_SRC_TREE=$(CFFC_SRC_TREE) wipe;
	cd $(SPARSELIB_DIR); make CFFC_SRC_TREE=$(CFFC_SRC_TREE) clean;
	cd $(BPKIT_DIR)/src; make CFFC_SRC_TREE=$(CFFC_SRC_TREE) clean;
	cd $(BLKCLIB_DIR); make distclean;	
	@echo ' '
	@echo All documentation and library directories are now clean.
	@echo '-------------------------------------------------------'
	@echo ' '

#
# Dependencies for Fortran 90 objects/modules
$(OBJ_REACTIONS_F90): $(MOD_MATH_F90) $(SRC_REACTIONS_F90)

$(MOD_MATH_F90): $(MATH_DIR)complexify.f90 $(MATH_DIR)Complexify.h
	$(F90) $(F90FLAGS) $(MODFLAG) $(MODULE_DIR) -I$(MODULE_DIR) -c -o $(MATH_DIR)complexify.o $(MATH_DIR)complexify.f90

#
# Default dependencies for Fortran 90 objects
.SUFFIXES: .f90
.f90.o:
	$(F90) $(F90FLAGS) $(MODFLAG) $(MODULE_DIR) -I$(MODULE_DIR) -c -o $*.o $*.f90

###################################################################
#
# Include header file dependencies
#
#------------------------------------------------------------------
# For some commands (targets) there is no need to include the
# dependencies. Therefore a list of targets for which the dependencies
# DON'T need to be generated is provided.
# If you don't want to have dependencies generated for a particular target,
# just add the target name to 'NODEP_LIST' variable.
NODEP_LIST = help clean superclean clean-extras libraries-all \
             bpkit blkclib sparselib fftw cantera documentation \
             uninstall uninstall-all

#-------------------------------------------------
# Include header file dependencies when required
#-------------------------------------------------
ifdef MAKECMDGOALS
ifeq ($(strip $(filter $(MAKECMDGOALS), $(NODEP_LIST)) ), )
############# Include dependencies of files that have extensions different than " .cc " (e.g. *.c, *.f, *.f90)
	sinclude $(patsubst %.c,$(DEPEND_DIR)%.d, $(filter %.c, $(SRC_CFFC)))
	sinclude $(patsubst %.f,$(DEPEND_DIR)%.d, $(filter %.f, $(SRC_CFFC)))
	sinclude $(patsubst %.f90,$(DEPEND_DIR)%.d, $(filter %.f90, $(SRC_CFFC)))
############# Include the .d files coming from different executables
# This command includes most of the *.d files
	sinclude $(patsubst %.cc,$(DEPEND_DIR)%.d, $(filter %.cc, $(SRC_CFFC)))
# If you add a new executable, you must add a line here.
# all other files should be included in the previous rule.
	sinclude $(patsubst %.cc,$(DEPEND_DIR)%.d, $(filter %.cc, $(SRC_MAIN)))
	sinclude $(patsubst %.cc,$(DEPEND_DIR)%.d, $(filter %.cc, $(SRC_EULER_POLYTROPIC_EXE)))
	sinclude $(patsubst %.cc,$(DEPEND_DIR)%.d, $(filter %.cc, $(SRC_EULER_THERMPERFECT_EXE)))
	sinclude $(patsubst %.cc,$(DEPEND_DIR)%.d, $(filter %.cc, $(SRC_NAVIERSTOKES_THERMPERFECT_EXE)))
	sinclude $(patsubst %.cc,$(DEPEND_DIR)%.d, $(filter %.cc, $(SRC_FANS_THERMPERFECT_EXE)))
endif
endif

##################################################################################################
#------------------------------------- IMPLICIT RULES --------------------------------------------
##################################################################################################

# Rule for building automatic dependencies (#include <file_name>)
# These file are assigned '.d' extension
# See "http://www.ddj.com/linux-open-source/184406479" for a thorough explanation.
$(DEPEND_DIR)%.d : %.cc
	@mkdir -p $(dir $@)
	@makedepend -f- $(DEPENDSFLAGS) $(CPPFLAGS) $< 2>/dev/null \
        | sed 's,\($*\.o\)[ :]*\(.*\),$@ : $$\(wildcard \2\) =\1 : \2,g' | tr "=" "\n" > $@

$(DEPEND_DIR)%.d : %.c
	@mkdir -p $(dir $@)
	@makedepend -f- $(DEPENDSFLAGS) $(CPPFLAGS) $< 2>/dev/null \
        | sed 's,\($*\.o\)[ :]*\(.*\),$@ : $$\(wildcard \2\) \n\1 : \2,g' | tr "=" "\n" > $@

$(DEPEND_DIR)%.d : %.f90
	@mkdir -p $(dir $@)
	@makedepend -f- $(DEPENDSFLAGS) $(F90FLAGS) $< 2>/dev/null \
        | sed 's,\($*\.o\)[ :]*\(.*\),$@ : $$\(wildcard \2\) \n\1 : \2,g' | tr "=" "\n" > $@

$(DEPEND_DIR)%.d : %.f
	@mkdir -p $(dir $@)
	@makedepend -f- $(DEPENDSFLAGS) $(FFLAGS) $< 2>/dev/null \
        | sed 's,\($*\.o\)[ :]*\(.*\),$@ : $$\(wildcard \2\) \n\1 : \2,g' | tr "=" "\n" > $@

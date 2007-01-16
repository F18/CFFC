#
# Makefile for Reconstruction program
# Computer program for reconstructing functions
# based on the average values in the computational cells.
#

# Set Reconstruction directory location which 
# can be overridden by setting an environement 
# variable Reconstruction_Path
Reconstruction_Path= $(shell cd ../; pwd)
export Reconstruction_Path

# Platform compile options
include makefile.def

# Assign sub directories
TEST_DIR		= UnitTesting/Tests
CFD_DIR                 = CFD
GRID1D_DIR              = Grid/Grid1D
GRID2D_DIR              = Grid/Grid2D
GRID3D_DIR              = Grid/Grid3D
MATH_DIR                = Math
RECONSTRUCTION1D_DIR    = Reconstruction/Reconstruction1D
RECONSTRUCTION2D_DIR    = Reconstruction/Reconstruction2D
RECONSTRUCTION3D_DIR    = Reconstruction/Reconstruction3D
TESTFUNCTIONS_DIR       = TestFunctions
COMMON_DIR              = Common

# Source files
#
SRC_TUT			= Test_Run.cc
SRC_TUT                 := $(addprefix  $(TUT_DIR)/, $(SRC_TUT))
#
SRC_CFD                 = CFD.cc Euler2DState.cc Gaussian2DState.cc
SRC_CFD			:= $(addprefix  $(CFD_DIR)/, $(SRC_CFD))
#
SRC_MATH		= LinearSystems.cc Matrix.cc NumericalLibrary.cc
SRC_MATH		:= $(addprefix  $(MATH_DIR)/, $(SRC_MATH))
#
SRC_MATH_UnitTesting    = test_LinearSystems.cc test_NumericalLibrary.cc \
			  test_Vector2D.cc  test_HeaderData.cc \
			  test_ComputationalCell.cc

SRC_MATH_UnitTesting	:= $(addprefix  $(TEST_DIR)/, $(SRC_MATH_UnitTesting))
#
SRC_COMMON              = StateData.cc
SRC_COMMON		:= $(addprefix  $(COMMON_DIR)/, $(SRC_COMMON))
#
SRC_CONTAINERS_UnitTesting = test_PointWiseSolution.cc test_SubGridMesh.cc test_TaylorDerivatives.cc \
			     test_ReconstructionSolver.cc test_ComputationalDomain.cc
SRC_CONTAINERS_UnitTesting := $(addprefix  $(TEST_DIR)/, $(SRC_CONTAINERS_UnitTesting))
##################################################################################
#                          ONE-DIMENSION
##################################################################################
SRC_GRID1D              = Cell1D.cc
SRC_GRID1D		:= $(addprefix $(GRID1D_DIR)/, $(SRC_GRID1D))
#
SRC_1D_TEST_FUNCTIONS   = TestFunctions.cc
SRC_1D_TEST_FUNCTIONS   := $(addprefix $(TESTFUNCTIONS_DIR)/, $(SRC_1D_TEST_FUNCTIONS))
#
SRC_1D_UnitTesting	= test_Cell1D.cc test_Cell1DNonUniform.cc  \
		         test_CompCell1D.cc

SRC_1D_UnitTesting	:= $(addprefix $(TEST_DIR)/, $(SRC_1D_UnitTesting))
#
#*********************************
# SRC_RECONSTRUCTION 1D Main Files
#*********************************
SRC_RECONSTRUCTION1D	= CompCell1D.cc CompDomain1D.cc Reconstruct1DInput.cc \
			  Reconstruction1DSolver.cc 
SRC_RECONSTRUCTION1D	:= $(addprefix  $(RECONSTRUCTION1D_DIR)/, $(SRC_RECONSTRUCTION1D))
#*********************************
# SRC_RECONSTRUCTION 1D Adds-on Files
#*********************************
SRC_RECONSTRUCTION1D	:= $(SRC_RECONSTRUCTION1D) $(SRC_GRID1D) $(SRC_1D_TEST_FUNCTIONS)
#
#
##################################################################################
#                          TWO-DIMENSIONS
##################################################################################
SRC_GRID2D		= Cell2D.cc Grid2DQuadMultiBlock.cc  Grid2DQuadSingleBlock.cc Spline2D.cc \
			  QuadrilateralGrid.cc Grid2DQuad.cc
SRC_GRID2D		:= $(addprefix $(GRID2D_DIR)/, $(SRC_GRID2D))
#
SRC_2D_TEST_FUNCTIONS   = TestFunctions_2D.cc
SRC_2D_TEST_FUNCTIONS   := $(addprefix $(TESTFUNCTIONS_DIR)/, $(SRC_2D_TEST_FUNCTIONS))
#
SRC_2D_UnitTesting	= test_Cell2D.cc test_QuadrilateralGrid.cc
SRC_2D_UnitTesting	:= $(addprefix $(TEST_DIR)/, $(SRC_2D_UnitTesting))
#
#*********************************
# SRC_RECONSTRUCTION 2D Main Files
#*********************************
SRC_RECONSTRUCTION2D	= Reconstruct2DInput.cc Reconstruction2DSolver.cc
SRC_RECONSTRUCTION2D	:= $(addprefix  $(RECONSTRUCTION2D_DIR)/, $(SRC_RECONSTRUCTION2D))
#*********************************
# SRC_RECONSTRUCTION 2D Adds-on Files
#*********************************
SRC_RECONSTRUCTION2D	:= $(SRC_RECONSTRUCTION2D) $(SRC_GRID2D) $(SRC_2D_TEST_FUNCTIONS)
#
#
##################################################################################
#                          THREE-DIMENSIONS
##################################################################################
#
SRC_GRID3D		= 
SRC_GRID3D		:= $(addprefix $(GRID3D_DIR)/, $(SRC_GRID3D))
#
SRC_3D_TEST_FUNCTIONS   = TestFunctions_3D.cc
SRC_3D_TEST_FUNCTIONS   := $(addprefix $(TESTFUNCTIONS_DIR)/, $(SRC_3D_TEST_FUNCTIONS))
#
#*********************************
# SRC_RECONSTRUCTION 3D Main Files
#*********************************
SRC_RECONSTRUCTION3D	= Reconstruction3DSolver.cc
SRC_RECONSTRUCTION3D	:= $(addprefix  $(RECONSTRUCTION3D_DIR)/, $(SRC_RECONSTRUCTION3D))
#*********************************
# SRC_RECONSTRUCTION 2D Adds-on Files
#*********************************
SRC_RECONSTRUCTION3D	:= $(SRC_RECONSTRUCTION3D) $(SRC_GRID3D) $(SRC_3D_TEST_FUNCTIONS)
#
#
SRC_UnitTesting         = $(SRC_1D_UnitTesting) $(SRC_2D_UnitTesting) $(SRC_MATH_UnitTesting) \
			  $(SRC_CONTAINERS_UnitTesting)
SRC_RECONSTRUCTION	= main.cc Reconstruction/ReconstructionHelpers.cc \
			  $(SRC_RECONSTRUCTION1D) $(SRC_RECONSTRUCTION2D) $(SRC_RECONSTRUCTION3D) \
			  $(SRC_CFD) $(SRC_MATH) \
			  $(SRC_UnitTesting) $(SRC_TUT) $(SRC_COMMON)

# Object files
#
OBJ_RECONSTRUCTION   :=$(SRC_RECONSTRUCTION:.cc=.o)

# Define the executables
EXE_REC = reconstruct

# Define targets

# Phony targets
.PHONY: help all allnew sparselib bpkit clean show

help:
	@echo "+-------------------------------------------------------------------+"
	@echo "|                                                                   |"
	@echo "|           makefile for Space Reconstruction                       |"
	@echo "|                                                                   |"
	@echo "| Usage: make all                 compile all executable programs   |"
	@echo "|        make allnew              recompile all executable programs |"
	@echo "|                                 and the libraries.                |"
	@echo "|        make clean               delete all the intermediate files |"
	@echo "+-------------------------------------------------------------------+"

allnew:
	@$(MAKE) --no-print-directory sparselib
	@$(MAKE) --no-print-directory bpkit
	@$(MAKE) --no-print-directory $(EXE_REC)

all: $(EXE_REC)
	@echo ' '
	@echo All programs have been brought up to date.
	@echo '------------------------------------------'
	@echo ' '

$(EXE_REC) : $(OBJ_RECONSTRUCTION) $(EXE_REC).depend
	@echo ' '
	@echo Building program $(EXE_REC) for a $(PLATFORM) platform.
	@echo '--------------------------------------------------------'
	@echo ' '
	$(LD) $(OBJ_RECONSTRUCTION) $(LDFLAGS) -o $@
	@echo ' '
	@echo Program $(EXE_REC) has been brought up to date.
	@echo '------------------------------------------------'
	@echo ' '

$(EXE_REC).depend: $(SRC_RECONSTRUCTION)
	@echo ' '
	@echo Determining Reconstruction source file dependencies for a $(PLATFORM) platform.
	@echo '--------------------------------------------------------------------------------'
	@echo ' '
	touch $(EXE_REC).depend
	makedepend -f $(EXE_REC).depend $(DEPENDSFLAGS) $(CPPFLAGS) $(SRC_RECONSTRUCTION) >& /dev/null
	@echo ' '
	@echo Reconstruction source file dependencies have been brought up to date.
	@echo '----------------------------------------------------------------------'
	@echo ' '

# Build libraries
show:
	@echo $(SPARSELIB_DIR)

sparselib:
	@echo ' '
	@echo Building $(SPARSELIB) library for a $(PLATFORM) platform.
	@echo '--------------------------------------------------------'
	@echo ' '
	cd $(SPARSELIB_DIR); make clean;
	cd $(SPARSELIB_DIR); make sp;
	@echo ' '
	@echo $(SPARSELIB) library has been brought up to date.
	@echo '------------------------------------------------'
	@echo ' '

bpkit:
	@echo ' '
	@echo Building $(BPKIT) library for a $(PLATFORM) platform.
	@echo '----------------------------------------------------'
	@echo ' '
	cd $(BPKIT_DIR)/src; make clean;
	cd $(BPKIT_DIR)/src; make;
	@echo ' '
	@echo $(BPKIT) library has been brought up to date.
	@echo '--------------------------------------------'
	@echo ' '

# Clean directory
clean:
	@echo ' '
	@echo Cleaning source file directory.
	@echo '-------------------------------'
	@echo ' '
	@echo "Remove *.depend  *.bak  *.o  *.obj  *.s  *.rtl  core* *~ #*# $(EXE_REC)"
	@../miscellaneous/rmlots "*.depend"
	@../miscellaneous/rmlots "*.bak"
	@../miscellaneous/rmlots "*.o"
	@../miscellaneous/rmlots "*.obj"
	@../miscellaneous/rmlots "*.s"
	@../miscellaneous/rmlots "*.rtl"
	@../miscellaneous/rmlots "core*"
	@../miscellaneous/rmlots "*~"
	@rm -f $(EXE_REC) $(TUT_DIR)/*.o
	@echo ' '
	@echo Source file directory is now clean.
	@echo '-----------------------------------'
	@echo ' '

#Include header file dependencies
ifneq ($(MAKECMDGOALS),help)
 ifneq ($(MAKECMDGOALS),clean)
   include $(EXE_REC).depend
 endif
endif


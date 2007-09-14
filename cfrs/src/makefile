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

# File SVN
SVN_File = Common/SourceRevisionData.cc

# Assign sub directories
DEPENDDIR               = Depend/
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
SRC_COMMON              = StateData.cc SourceRevisionData.cc EpsilonTol.cc
SRC_COMMON		:= $(addprefix  $(COMMON_DIR)/, $(SRC_COMMON))
#
SRC_COMMON_UnitTesting  = test_SourceCode.cc BasicTestData.cc test_BasicTestData.cc

SRC_COMMON_UnitTesting	:= $(addprefix  $(TEST_DIR)/, $(SRC_COMMON_UnitTesting))
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
SRC_1D_UnitTesting	= test_Cell1D.cc test_Cell1DNonUniform.cc

SRC_1D_UnitTesting	:= $(addprefix $(TEST_DIR)/, $(SRC_1D_UnitTesting))
#
#*********************************
# SRC_RECONSTRUCTION 1D Main Files
#*********************************
SRC_RECONSTRUCTION1D	=  Reconstruct1DInput.cc \
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
			  $(SRC_CONTAINERS_UnitTesting) $(SRC_COMMON_UnitTesting)
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
	@$(MAKE) --no-print-directory clean
	@$(MAKE) --no-print-directory sparselib
	@$(MAKE) --no-print-directory bpkit
	@$(MAKE) --no-print-directory $(EXE_REC)

all: $(EXE_REC)
	@echo ' '
	@echo All programs have been brought up to date.
	@echo '------------------------------------------'
	@echo ' '

$(EXE_REC) : $(OBJ_RECONSTRUCTION)
	@echo ' '
	@echo Update compilation time and SVN revision data
	@echo '--------------------------------------------------------'
	@python ../miscellaneous/Update_SourceRevision_Data.py
        # Update $(SVN_File,.cc=.o)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(SVN_File) -c -o $(SVN_File:.cc=.o)
        # The updated SourceRevisionData is not required anymore ==> Revert SVN_File to the default (original)
	@svn revert $(SVN_File) >& /dev/null
	@echo ' '
	@echo Building program $(EXE_REC) for a $(PLATFORM) platform.
	@echo '--------------------------------------------------------'
	@echo ' '
	$(LD) $(OBJ_RECONSTRUCTION) $(LDFLAGS) -o $@
	@echo ' '
	@echo Program $(EXE_REC) has been brought up to date.
	@echo '------------------------------------------------'
	@echo ' '

# Build libraries
show:
	@echo $(SPARSELIB_DIR)
	@echo $(BPKIT_DIR)

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
	@echo "Remove *.bak  *.o  *.obj  *.s  *.rtl  core* *~ #*# $(EXE_REC)"
	@../miscellaneous/rmlots "*.bak"
	@../miscellaneous/rmlots "*.o"
	@../miscellaneous/rmlots "*.obj"
	@../miscellaneous/rmlots "*.s"
	@../miscellaneous/rmlots "*.rtl"
	@../miscellaneous/rmlots "core*"
	@../miscellaneous/rmlots "*~"
	@../miscellaneous/rmlots "#*#"
	@rm -f $(TUT_DIR)/*.o
	@echo "Remove the automatically created dependency files *.d"
	rm -fr $(DEPENDDIR)
	rm -f $(EXE_REC)
	@echo ' '
	@echo Source file directory is now clean.
	@echo '-----------------------------------'
	@echo ' '

clean-sparselib:
	@echo ' '
	@echo Cleaning $(SPARSELIB) library.
	@echo '--------------------------------------------------------'
	@echo ' '
	cd $(SPARSELIB_DIR); make clean;
	@echo ' '
	@echo $(SPARSELIB) library has been cleaned.
	@echo '------------------------------------------------'
	@echo ' '

clean-bpkit:
	@echo ' '
	@echo Cleaning $(BPKIT) library.
	@echo '----------------------------------------------------'
	@echo ' '
	cd $(BPKIT_DIR)/src; make clean;
	@echo ' '
	@echo $(BPKIT) library has been cleaned.
	@echo '--------------------------------------------'
	@echo ' '

clean-lib: clean-sparselib clean-bpkit clean-cgns
	@echo ' '
	@echo Libraries directories cleaned.
	@echo '----------------------------------------------------'

#
# Include header file dependencies
#
#------------------------------------------------------------------
# For some commands (targets) there is no need to include the 
# dependencies. Therefore a list of targets for which the dependencies
# DON'T need to be generated is provided.
# If you don't want to have dependencies generated, just add the target 
# to NODEP_LIST.
#
NODEP_LIST = help clean clean-lib clean-sparselib clean-bpkit \
	     allnew


#-----------------------------------
# Include header file dependencies
#-----------------------------------
ifeq ($(strip $(filter $(MAKECMDGOALS), $(NODEP_LIST)) ), )
############# Include dependencies of files that have extensions different than " .cc " (e.g. *.c)
       # sinclude $(patsubst %.c,$(DEPENDDIR)%.d, $(filter %.c, $(SRC_RECONSTRUCTION)))
############# Include the .d files comming from different executables
       # this command includes most of the *.d files
       sinclude $(patsubst %.cc,$(DEPENDDIR)%.d, $(filter %.cc, $(SRC_RECONSTRUCTION)))
endif

##################################################################################################
#------------------------------------- IMPLICIT RULES --------------------------------------------
##################################################################################################

# Rule for building automatic dependencies (#include <file_name>)
# These file are assigned '.d' extension
$(DEPENDDIR)%.d : %.cc
	@mkdir -p $(dir $@)
	@makedepend -f- $(DEPENDSFLAGS) $(CPPFLAGS) $< 2>/dev/null \
        | sed 's,\($*\.o\)[ :]*\(.*\),$@ : $$\(wildcard \2\) =\1 : \2,g' | tr "=" "\n" > $@ 

$(DEPENDDIR)%.d : %.c
	@mkdir -p $(dir $@)
	@makedepend -f- $(DEPENDSFLAGS) $(CPPFLAGS) $< 2>/dev/null \
        | sed 's,\($*\.o\)[ :]*\(.*\),$@ : $$\(wildcard \2\) \n\1 : \2,g' | tr "=" "\n" > $@


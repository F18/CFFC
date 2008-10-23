#
# Makefile for CFRS program
# 
# Computational Framework for Reconstruction Studies (CFRS)
#

# Set Reconstruction directory location which 
# can be overridden by setting an environement 
# variable Reconstruction_Path
Reconstruction_Path= $(shell cd ../; pwd)
export Reconstruction_Path

#
# Set CFFC directory location which can be overridden by 
# setting an environement variable CFFC_Path
CFFC_Path:= $(shell cd ../../; pwd)
export CFFC_Path

CFFC_2D_Path:= $(CFFC_Path)/src_2D
export CFFC_2D_Path

#
# Set CFFC source tree
CFFC_SRC_TREE = 2D
export CFFC_SRC_TREE

# Platform compile options
include makefile.def

# File GIT
GIT_File = Common/SourceRevisionData.cc

# Assign sub directories
DEPENDDIR               = Depend/
CFD_DIR                 = CFD
GRID1D_DIR              = Grid/Grid1D
GRID2D_DIR              = Grid/Grid2D
GRID3D_DIR              = Grid/Grid3D
MATH_DIR                = $(CFFC_2D_Path)/Math
MATH_TESTS_DIR          = Math/UnitTests
RECONSTRUCTION1D_DIR    = Reconstruction/Reconstruction1D
RECONSTRUCTION2D_DIR    = Reconstruction/Reconstruction2D
RECONSTRUCTION3D_DIR    = Reconstruction/Reconstruction3D
TESTFUNCTIONS_DIR       = TestFunctions
COMMON_DIR              = Common
CONTAINERS_TEST_DIR     = Common/Containers/UnitTests

###############################################################################
# C++, C, and Fortran Source Files
###############################################################################
# TUT testing:
SRC_TUT = TestRun.cc \
          TestData.cc \
	  BasicTests/Sample_Test.cc \
          BasicTests/test_TestData.cc
SRC_TUT := $(addprefix $(UNIT_TESTING_DIR), $(SRC_TUT))

# CFD:
SRC_CFD                 = CFD.cc Euler2DState.cc Gaussian2DState.cc
SRC_CFD			:= $(addprefix  $(CFD_DIR)/, $(SRC_CFD))

# Math: (files from CFFC2D)
SRC_MATH		= LinearSystems.cc Matrix.cc NumericalLibrary.cc
SRC_MATH		:= $(addprefix  $(MATH_DIR)/, $(SRC_MATH))
#--- SRC_MATH_TESTS ----
SRC_MATH_TESTS = test_LeastSquaresSolvers.cc \
	         test_NumericalLibrary.cc
#
SRC_MATH_TESTS	:= $(addprefix  $(MATH_TESTS_DIR)/, $(SRC_MATH_TESTS))
#--- Add SRC_MATH_TESTS to SRC_MATH ----
ifeq ($(TUT_TESTING),ON) 
  SRC_MATH += $(SRC_MATH_TESTS)
endif
#

# Common:
SRC_COMMON              = StateData.cc SourceRevisionData.cc
SRC_COMMON		:= $(addprefix  $(COMMON_DIR)/, $(SRC_COMMON))

# Common from CFFC_2D:
SRC_COMMON              += $(CFFC_2D_Path)/Utilities/EpsilonTol.cc

# Common from other directories:
SRC_COMMON              += Reconstruction/CENO_Tolerances.cc


#--- SRC_COMMON_TESTS ----
SRC_COMMON_TESTS  = 
#
#--- Add SRC_COMMON_TESTS to SRC_COMMON ----
ifeq ($(TUT_TESTING),ON) 
  SRC_COMMON += $(SRC_COMMON_TESTS)
endif
#

# Tests for Containers:
SRC_CONTAINERS_TESTS = $(CONTAINERS_TEST_DIR)/test_PointWiseSolution.cc \
		       $(CONTAINERS_TEST_DIR)/test_TaylorDerivatives.cc \
		       Reconstruction/UnitTests/test_ComputationalCell.cc \
		       Reconstruction/UnitTests/test_ComputationalDomain.cc

# Tests for SubGrid:
SRC_SUBGRID_TESTS := test_SubGridMesh.cc
SRC_SUBGRID_TESTS := $(addprefix Grid/SubGrid/UnitTests/, $(SRC_SUBGRID_TESTS))


##################################################################################
#                          ONE-DIMENSION
##################################################################################
# Grid:
SRC_GRID1D              = Cell1D.cc
SRC_GRID1D		:= $(addprefix $(GRID1D_DIR)/, $(SRC_GRID1D))

# Test functions:
SRC_1D_TEST_FUNCTIONS   = TestFunctions_1D.cc
SRC_1D_TEST_FUNCTIONS   := $(addprefix $(TESTFUNCTIONS_DIR)/, $(SRC_1D_TEST_FUNCTIONS))

# Reconstruction:
#*********************************
SRC_RECONSTRUCTION1D	=  Reconstruct1DInput.cc \
			  Reconstruction1DSolver.cc 
SRC_RECONSTRUCTION1D	:= $(addprefix  $(RECONSTRUCTION1D_DIR)/, $(SRC_RECONSTRUCTION1D))

#*********************************
# Final SRC_RECONSTRUCTION1D
#*********************************
SRC_RECONSTRUCTION1D	+= $(SRC_GRID1D) $(SRC_1D_TEST_FUNCTIONS)
#


##################################################################################
##################################################################################
#                          TWO-DIMENSIONS
##################################################################################
##################################################################################
# Grid:
SRC_GRID2D		= Cell2D.cc Grid2DQuadMultiBlock.cc  Grid2DQuadSingleBlock.cc Spline2D.cc \
			  QuadrilateralGrid.cc Grid2DQuad.cc
SRC_GRID2D		:= $(addprefix $(GRID2D_DIR)/, $(SRC_GRID2D))
#--- SRC_GRID2D_TESTS ----
SRC_GRID2D_TESTS  = test_Cell2D.cc  test_QuadrilateralGrid.cc
SRC_GRID2D_TESTS := $(addprefix Grid/Grid2D/UnitTests/, $(SRC_GRID2D_TESTS))

#
#--- Add SRC_GRID2D_TESTS to SRC_GRID2D ----
ifeq ($(TUT_TESTING),ON) 
  SRC_GRID2D += $(SRC_GRID2D_TESTS)

endif
#

# Test functions:
SRC_2D_TEST_FUNCTIONS   = TestFunctions_2D.cc
SRC_2D_TEST_FUNCTIONS   := $(addprefix $(TESTFUNCTIONS_DIR)/, $(SRC_2D_TEST_FUNCTIONS))

# Reconstruction:
#*********************************
SRC_RECONSTRUCTION2D	= Reconstruct2DInput.cc Reconstruction2DSolver.cc
SRC_RECONSTRUCTION2D	:= $(addprefix  $(RECONSTRUCTION2D_DIR)/, $(SRC_RECONSTRUCTION2D))

#*********************************
# Final SRC_RECONSTRUCTION2D
#*********************************
SRC_RECONSTRUCTION2D	+= $(SRC_GRID2D) $(SRC_2D_TEST_FUNCTIONS)
#

##################################################################################
##################################################################################
#                          THREE-DIMENSIONS
##################################################################################
##################################################################################
# Grid:
SRC_GRID3D		= Grid3DHexaBlock.cc
SRC_GRID3D		:= $(addprefix $(GRID3D_DIR)/, $(SRC_GRID3D))
#

# Test functions:
SRC_3D_TEST_FUNCTIONS   = TestFunctions_3D.cc
SRC_3D_TEST_FUNCTIONS   := $(addprefix $(TESTFUNCTIONS_DIR)/, $(SRC_3D_TEST_FUNCTIONS))


# Reconstruction:
#*********************************
SRC_RECONSTRUCTION3D	= Reconstruct3DInput.cc Reconstruction3DSolver.cc
SRC_RECONSTRUCTION3D	:= $(addprefix  $(RECONSTRUCTION3D_DIR)/, $(SRC_RECONSTRUCTION3D))

#*********************************
# Final SRC_RECONSTRUCTION3D
#*********************************
SRC_RECONSTRUCTION3D	+= $(SRC_GRID3D) $(SRC_3D_TEST_FUNCTIONS)
#

##################################################################################
##################################################################################
#
#                            SRC_RECONSTRUCTION
#
##################################################################################
##################################################################################
SRC_RECONSTRUCTION	= main.cc Reconstruction/ReconstructionHelpers.cc \
			  $(CFFC_2D_Path)/HighOrderReconstruction/ReconstructionHelpers.cc \
			  $(SRC_RECONSTRUCTION1D) $(SRC_RECONSTRUCTION2D) $(SRC_RECONSTRUCTION3D) \
			  $(SRC_CFD) $(SRC_MATH) \
			  $(SRC_TUT) $(SRC_COMMON)


#--- Add SRC_CONTAINERS_TESTS to SRC_RECONSTRUCTION ----
ifeq ($(TUT_TESTING),ON) 
  SRC_RECONSTRUCTION += $(SRC_CONTAINERS_TESTS)
  SRC_RECONSTRUCTION += $(SRC_SUBGRID_TESTS)
endif
#



#**********************  Object files ******************
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
	@echo Update compilation time and GIT revision data
	@echo '--------------------------------------------------------'
	@python ../miscellaneous/Update_SourceRevision_Data.py -f $(GIT_File)
        # Update $(GIT_File,.cc=.o)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(GIT_File) -c -o $(GIT_File:.cc=.o)
        # The updated SourceRevisionData is not required anymore ==> Revert GIT_File to the default (original)
	@python ../miscellaneous/Revert_Git_File.py -f $(GIT_File)
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
	@rm -f $(UNIT_TESTING_DIR)/*.o
	@echo Remove obj files from CFFC_2D
	@rm -f $(CFFC_2D_Path)/Utilities/EpsilonTol.o $(CFFC_2D_Path)/HighOrderReconstruction/ReconstructionHelpers.o
	@rm -f $(CFFC_2D_Path)/Math/LinearSystems.o $(CFFC_2D_Path)/Math/NumericalLibrary.o
	@rm -f $(CFFC_2D_Path)/Math/Matrix.o
	@rm -f $(CFFC_2D_Path)/Math/UnitTests/*.o $(CFFC_2D_Path)/Math/UnitTests/TestFunctions/*.o
	@rm -f $(CFFC_2D_Path)/UnitTesting/BasicTests/test_TestData.o $(CFFC_2D_Path)/UnitTesting/BasicTests/Sample_Test.o
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
ifdef MAKECMDGOALS
ifeq ($(strip $(filter $(MAKECMDGOALS), $(NODEP_LIST)) ), )
############# Include dependencies of files that have extensions different than " .cc " (e.g. *.c)
       # sinclude $(patsubst %.c,$(DEPENDDIR)%.d, $(filter %.c, $(SRC_RECONSTRUCTION)))
############# Include the .d files comming from different executables
       # this command includes most of the *.d files
       sinclude $(patsubst %.cc,$(DEPENDDIR)%.d, $(filter %.cc, $(SRC_RECONSTRUCTION)))
endif
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


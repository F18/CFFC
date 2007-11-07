#
# Makefile for CFFC (2D version)
#
# Computational Framework for Fluids and Combustion (CFFC)

#
# Set CFFC directory location which can be overridden by 
# setting an environement variable CFFC_Path
CFFC_Path:= $(shell cd ../; pwd)

#
# Set CFFC source tree
CFFC_SRC_TREE = 2D

#
# Set standard directories
CFFC_BIN_DIR = $(CFFC_Path)/bin

#
# Platform compile options 
include makefile.def

#
# Assign sub directories
DEPEND_DIR = Depend/
ADVECTDIFFUSE_DIR = AdvectDiffuse2D/
AMR_DIR = AMR/
CFD_DIR = CFD/
CHEM_DIR = Chem2D/
FLAME_DIR = Flame2D/
DUSTY_DIR = Dusty2D/
ELECTROSTATIC_DIR = Electrostatic2D/
EOS_DIR = EquationOfState/
EULER_1D_DIR = Euler1D/
EULER_2D_DIR = Euler2D/
LEVERMORE_1D_DIR = Levermore1D/
GRID_DIR = Grid/
HIGHTEMP_DIR = HighTemp2D/
ICEMCFD_DIR = ICEM/
INTERFACE_DIR = Interface2D/
ION_DIR = Ion5Moment2D/
LESPREMIXED_DIR = LESPremixed2D/
LEVELSET_DIR = LevelSet2D/
RTE_DIR = Rte2D/
MATH_DIR = Math/
MPI_DIR = MPI/
MULTIGRID_DIR = FASMultigrid2D/
NAVIERSTOKES_DIR = NavierStokes2D/
NKS_DIR = NewtonKrylovSchwarz2D/
PHYSICS_DIR = Physics/
SOLVERS1D_DIR = Solvers1D/
SYSTEM_DIR = System/
TURBULENT_DIR = Turbulent2D/
GAUSSIAN_DIR = Gaussian2D/
NASADATA_DIR = $(PHYSICS_DIR)NASAData/
SNBCK_DIR = $(PHYSICS_DIR)SNBCK/
REACTIONS_DIR = Reactions/
HIGH_ORDER_RECONSTRUCTION_DIR = HighOrderReconstruction/
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
SRC_CFD = $(CFD_DIR)CFD.cc
SRC_MPI = $(MPI_DIR)MPI.cc

# Grid:
SRC_GRID_1D = $(GRID_DIR)Cell1D.cc

SRC_GRID = $(GRID_DIR)Cell1D.cc \
           $(GRID_DIR)Cell2D.cc \
           $(GRID_DIR)Grid2DQuadSingleBlock.cc \
           $(GRID_DIR)Grid2DQuadMultiBlock.cc \
           $(GRID_DIR)Grid3DHexaSingleBlock.cc \
           $(GRID_DIR)Grid3DHexaMultiBlock.cc \
           $(GRID_DIR)HO_Cell2D.cc \
           $(GRID_DIR)HO_Spline2D.cc \
           $(GRID_DIR)HO_Spline2DInterval.cc \
	   $(GRID_DIR)HO_Grid2DQuad.cc \
           $(GRID_DIR)HO_Grid2DQuadMultiBlock.cc \
	   $(GRID_DIR)HO_Grid2DQuad_ExecutionMode.cc \
	   $(GRID_DIR)Tecplot_ExecutionMode.cc
SRC_ICEMCFD = $(ICEMCFD_DIR)ICEMCFD.cc \
	      $(ICEMCFD_DIR)ICEMCFD_ctype.c

# Add Grid Testing If Required:
SRC_GRID_TESTS := test_Grid2DQuadIntegration.cc \
		  test_HO_Grid2DQuad.cc \
		  test_HO_Grid2DQuadMultiBlock.cc \
		  test_HO_Spline2D.cc \
		  test_HO_Spline2D_Integration.cc \
		  test_HO_Spline2DInterval.cc \
		  test_HO_Cell2D.cc \
		  HO_Grid2DQuadMultiBlock_InputForTesting.cc
SRC_GRID_TESTS := $(addprefix $(GRID_DIR)UnitTests/,$(SRC_GRID_TESTS))
ifeq ($(TUT_TESTING),ON) 
  SRC_GRID += $(SRC_GRID_TESTS)
endif

# AMR:
ifeq ($(MPI_VERSION),NoMPI)
  SRC_AMR = $(AMR_DIR)AdaptiveBlock.cc \
	    $(AMR_DIR)AdaptiveBlock_NoMPI.cc \
            $(AMR_DIR)QuadTree.cc
else
  SRC_AMR = $(AMR_DIR)AdaptiveBlock.cc \
	    $(AMR_DIR)AdaptiveBlock_MPI.cc \
            $(AMR_DIR)QuadTree.cc
endif

# HighOrderReconstruction common for both 1D and 2D:
SRC_HIGH_ORDER_RECONSTRUCTION_COMMON = ReconstructionHelpers.cc \
                                       CENO_Tolerances.cc \
				       CENO_ExecutionMode.cc 
SRC_HIGH_ORDER_RECONSTRUCTION_COMMON := $(addprefix  $(HIGH_ORDER_RECONSTRUCTION_DIR), $(SRC_HIGH_ORDER_RECONSTRUCTION_COMMON))

# HighOrderReconstruction only for 1D:
SRC_HIGH_ORDER_RECONSTRUCTION_1D = ENO_ExecutionMode.cc \
				   AccuracyAssessment1D.cc
SRC_HIGH_ORDER_RECONSTRUCTION_1D := $(addprefix  $(HIGH_ORDER_RECONSTRUCTION_DIR), $(SRC_HIGH_ORDER_RECONSTRUCTION_1D))
# Add High-Order Reconstruction Testing If Required:
SRC_HIGH_ORDER_RECONSTRUCTION_TESTS_1D := test_HighOrder1D_Double.cc
SRC_HIGH_ORDER_RECONSTRUCTION_TESTS_1D := $(addprefix $(HIGH_ORDER_RECONSTRUCTION_DIR)UnitTests/,$(SRC_HIGH_ORDER_RECONSTRUCTION_TESTS_1D))
ifeq ($(TUT_TESTING),ON) 
  SRC_HIGH_ORDER_RECONSTRUCTION_1D += $(SRC_HIGH_ORDER_RECONSTRUCTION_TESTS_1D)
endif

# HighOrderReconstruction only for 2D:
SRC_HIGH_ORDER_RECONSTRUCTION_2D = AccuracyAssessment2DMultiBlock.cc \
				   HighOrder2D_Input.cc
SRC_HIGH_ORDER_RECONSTRUCTION_2D := $(addprefix  $(HIGH_ORDER_RECONSTRUCTION_DIR), $(SRC_HIGH_ORDER_RECONSTRUCTION_2D))

# Add High-Order Reconstruction Testing If Required:
SRC_HIGH_ORDER_RECONSTRUCTION_TESTS_2D := test_TaylorDerivatives2D_Double.cc \
				          test_HighOrder2D_Double.cc
SRC_HIGH_ORDER_RECONSTRUCTION_TESTS_2D := $(addprefix $(HIGH_ORDER_RECONSTRUCTION_DIR)UnitTests/,$(SRC_HIGH_ORDER_RECONSTRUCTION_TESTS_2D))
ifeq ($(TUT_TESTING),ON) 
  SRC_HIGH_ORDER_RECONSTRUCTION_2D += $(SRC_HIGH_ORDER_RECONSTRUCTION_TESTS_2D)
endif

# 1D:
SRC_CFD1D_INPUT = CFD1DInput.cc
SRC_CFD1D_INPUT := $(addprefix $(CFD_DIR), $(SRC_CFD1D_INPUT))

SRC_SCALAR = $(SOLVERS1D_DIR)Scalar1D.cc \
             $(SOLVERS1D_DIR)Scalar1DSolvers.cc
SRC_HEAT = $(SOLVERS1D_DIR)Heat1D.cc \
           $(SOLVERS1D_DIR)Heat1DSolvers.cc
SRC_HYPERHEAT = $(SOLVERS1D_DIR)HyperHeat1DState.cc \
                $(SOLVERS1D_DIR)HyperHeat1D.cc \
                $(SOLVERS1D_DIR)HyperHeat1DSolvers.cc
SRC_MHD1D = $(SOLVERS1D_DIR)MHD1DState.cc \
            $(SOLVERS1D_DIR)MHD1D.cc \
            $(SOLVERS1D_DIR)MHD1DSolvers.cc

# Euler 1D:
SRC_EULER1D = Euler1DState.cc \
	      Euler1D.cc \
	      Euler1DSolvers.cc \
	      ExactSolutions/ExactSolutions.cc
SRC_EULER1D := $(addprefix $(EULER_1D_DIR), $(SRC_EULER1D))

# Add Euler 1D Testing If Required:
SRC_EULER1D_TESTS = test_Euler1D_Simulations.cc \
		    test_HighOrder1D_EulerState.cc \
		    test_HighOrder1D_Euler1DUniformMesh.cc
SRC_EULER1D_TESTS := $(addprefix $(EULER_1D_DIR)UnitTests/, $(SRC_EULER1D_TESTS))
ifeq ($(TUT_TESTING),ON) 
  SRC_EULER1D += $(SRC_EULER1D_TESTS)
endif

SRC_EULER1D_EXE = $(EULER_1D_DIR)Euler1DMain.cc

SRC_EULER1D_SOLVER = $(SRC_MATH) $(SRC_CFD) $(SRC_CFD1D_INPUT) $(SRC_GRID_1D) \
                     $(SRC_SYSTEM) $(SRC_HIGH_ORDER_RECONSTRUCTION_COMMON) $(SRC_HIGH_ORDER_RECONSTRUCTION_1D) \
                     $(SRC_EULER1D) $(SRC_EULER1D_EXE) $(SRC_MPI)

# Solvers 1D:
SRC_SOLVERS1D = $(SRC_MATH) $(SRC_CFD) $(SRC_CFD1D_INPUT) $(SRC_GRID_1D) \
                $(SRC_SYSTEM) $(SRC_HIGH_ORDER_RECONSTRUCTION_COMMON) $(SRC_HIGH_ORDER_RECONSTRUCTION_1D) \
		$(SRC_SCALAR) $(SRC_HEAT) $(SRC_HYPERHEAT) \
                $(SRC_MHD1D) $(SRC_EULER1D) $(SRC_MPI)

SRC_SOLVERS1D_EXE = $(SOLVERS1D_DIR)Solvers1D.cc
SRC_SOLVERS1D_SOLVER = $(SRC_SOLVERS1D_EXE) $(SRC_SOLVERS1D)

# Levermore 1D:
SRC_LEVERMORE1D = #none yet
SRC_LEVERMORE1D := $(addprefix $(LEVERMORE_1D_DIR), $(SRC_LEVERMORE1D))

# Add Levermore 1D Testing If Required:
SRC_LEVERMORE1D_TESTS = test_Levermore1DVector.cc
SRC_LEVERMORE1D_TESTS := $(addprefix $(LEVERMORE_1D_DIR)UnitTests/, $(SRC_LEVERMORE1D_TESTS))
ifeq ($(TUT_TESTING),ON)
  SRC_LEVERMORE1D += $(SRC_LEVERMORE1D_TESTS)
endif

SRC_LEVERMORE1D_EXE = $(LEVERMORE_1D_DIR)Levermore1DMain.cc

# Advectdiffuse2D:
SRC_ADVECTDIFFUSE = $(ADVECTDIFFUSE_DIR)AdvectDiffuse2DState.cc \
	            $(ADVECTDIFFUSE_DIR)AdvectDiffuse2DInput.cc \
                    $(ADVECTDIFFUSE_DIR)AdvectDiffuse2DQuad.cc \
                    $(ADVECTDIFFUSE_DIR)AdvectDiffuse2DQuad_HighOrder.cc \
	            $(ADVECTDIFFUSE_DIR)AdvectDiffuse2DQuadSingleBlock.cc \
                    $(ADVECTDIFFUSE_DIR)AdvectDiffuse2DQuadMultiBlock.cc \
                    $(ADVECTDIFFUSE_DIR)AdvectDiffuse2DQuadGrid.cc \
                    $(ADVECTDIFFUSE_DIR)AdvectDiffuse2DQuadSolvers.cc \
		    $(ADVECTDIFFUSE_DIR)AdvectDiffuse2DExactSolutions.cc \
		    $(ADVECTDIFFUSE_DIR)ExactSolutions.cc \
                    $(ADVECTDIFFUSE_DIR)VelocityFields.cc \
                    $(ADVECTDIFFUSE_DIR)DiffusionFields.cc \
                    $(ADVECTDIFFUSE_DIR)SourceTermFields.cc \
                    $(ADVECTDIFFUSE_DIR)SourceFields.cc \
                    $(ADVECTDIFFUSE_DIR)AdvectDiffuse2DInflowField.cc \
                    $(ADVECTDIFFUSE_DIR)InflowFields.cc
SRC_ADVECTDIFFUSE_EXE = $(ADVECTDIFFUSE_DIR)AdvectDiffuse2D.cc

# Add AdvectDiffuse 2D Testing If Required:
SRC_ADVECTDIFFUSE_TESTS = test_AdvectDiffuse2DState.cc \
                          test_SourceTermFields.cc \
                          test_AdvectDiffuse2DQuad.cc \
			  test_AdvectDiffuse2DExactSolutions.cc \
			  test_AdvectDiffuse2D_Simulations.cc
SRC_ADVECTDIFFUSE_TESTS := $(addprefix $(ADVECTDIFFUSE_DIR)UnitTests/, $(SRC_ADVECTDIFFUSE_TESTS))
ifeq ($(TUT_TESTING),ON) 
  SRC_ADVECTDIFFUSE += $(SRC_ADVECTDIFFUSE_TESTS)
endif

# Euler2D:
SRC_EULER = $(EULER_2D_DIR)Euler2DState.cc \
            $(EULER_2D_DIR)Euler2DInput.cc \
            $(EULER_2D_DIR)Euler2DQuad.cc \
            $(EULER_2D_DIR)Euler2DQuadSingleBlock.cc \
            $(EULER_2D_DIR)Euler2DQuadMultiBlock.cc \
            $(EULER_2D_DIR)Euler2DQuadGrid.cc \
            $(EULER_2D_DIR)Euler2DQuadSolvers.cc 
SRC_EULER_EXE = $(EULER_2D_DIR)Euler2D.cc

# NavierStokes2D:
SRC_NAVIERSTOKES = $(NAVIERSTOKES_DIR)NavierStokes2DState.cc \
                   $(NAVIERSTOKES_DIR)NavierStokes2DInput.cc \
                   $(NAVIERSTOKES_DIR)NavierStokes2DQuad.cc \
                   $(NAVIERSTOKES_DIR)NavierStokes2DQuadSingleBlock.cc \
                   $(NAVIERSTOKES_DIR)NavierStokes2DQuadMultiBlock.cc \
                   $(NAVIERSTOKES_DIR)NavierStokes2DQuadIOSingleBlock.cc \
                   $(NAVIERSTOKES_DIR)NavierStokes2DQuadIOMultiBlock.cc \
                   $(NAVIERSTOKES_DIR)NavierStokes2DQuadTurbulenceSingleBlock.cc \
                   $(NAVIERSTOKES_DIR)NavierStokes2DQuadTurbulenceMultiBlock.cc \
                   $(NAVIERSTOKES_DIR)NavierStokes2DdRdU.cc \
                   $(NAVIERSTOKES_DIR)NavierStokes2DQuadGrid.cc \
                   $(NAVIERSTOKES_DIR)NavierStokes2DQuadSolvers.cc
SRC_NAVIERSTOKES_EXE = $(NAVIERSTOKES_DIR)NavierStokes2D.cc

# Chem2D:
SRC_CHEM = $(CHEM_DIR)Chem2DState.cc \
	   $(CHEM_DIR)Chem2DTools.cc \
	   $(CHEM_DIR)Chem2DInput.cc \
	   $(CHEM_DIR)Chem2DQuad.cc \
	   $(CHEM_DIR)Chem2DQuadGrid.cc \
	   $(CHEM_DIR)dRdU.cc \
	   $(CHEM_DIR)Chem2DQuadMultiBlock.cc \
	   $(CHEM_DIR)Chem2DQuadSingleBlock.cc \
	   $(CHEM_DIR)Chem2DQuadSolvers.cc
SRC_CHEM_EXE = $(CHEM_DIR)Chem2D.cc  

# Flame2D:
SRC_FLAME = $(FLAME_DIR)Flame2DState.cc \
	    $(FLAME_DIR)Mixture.cc \
	    $(FLAME_DIR)Flame2DTools.cc \
	    $(FLAME_DIR)Flame2DInput.cc \
	    $(FLAME_DIR)Flame2DQuad.cc \
	    $(FLAME_DIR)Flame2DQuadGrid.cc \
	    $(FLAME_DIR)Flame2DdRdU.cc \
	    $(FLAME_DIR)Flame2DQuadMultiBlock.cc \
	    $(FLAME_DIR)Flame2DQuadSingleBlock.cc \
	    $(FLAME_DIR)Flame2DQuadSolvers.cc \
	    $(FLAME_DIR)Soot2DState.cc
SRC_FLAME_EXE = $(FLAME_DIR)Flame2D.cc  

# HighTemp2D:
SRC_HIGHTEMP = $(HIGHTEMP_DIR)HighTemp2DState.cc \
               $(HIGHTEMP_DIR)HighTemp2DInput.cc  \
               $(HIGHTEMP_DIR)HighTemp2DQuad.cc  \
               $(HIGHTEMP_DIR)HighTemp2DQuadGrid.cc \
               $(HIGHTEMP_DIR)HighTemp2DQuadMultiBlock.cc  \
               $(HIGHTEMP_DIR)HighTemp2DQuadSingleBlock.cc  \
               $(HIGHTEMP_DIR)HighTemp2DQuadIOMultiBlock.cc  \
               $(HIGHTEMP_DIR)HighTemp2DQuadIOSingleBlock.cc  \
               $(HIGHTEMP_DIR)HighTemp2DQuadTurbulenceMultiBlock.cc  \
               $(HIGHTEMP_DIR)HighTemp2DQuadTurbulenceSingleBlock.cc  \
               $(HIGHTEMP_DIR)HighTemp2DdRdU.cc  \
               $(HIGHTEMP_DIR)HighTemp2DQuadSolvers.cc
SRC_HIGHTEMP_EXE = $(HIGHTEMP_DIR)HighTemp2D.cc

# Dusty2D:
SRC_DUSTY = $(DUSTY_DIR)Particle2DComponents.cc \
            $(DUSTY_DIR)Dusty2DState.cc \
            $(DUSTY_DIR)Dusty2DInput.cc \
            $(DUSTY_DIR)Dusty2DQuad.cc \
            $(DUSTY_DIR)Dusty2DQuadSingleBlock.cc \
            $(DUSTY_DIR)Dusty2DQuadMultiBlock.cc \
            $(DUSTY_DIR)Dusty2DQuadIOSingleBlock.cc \
            $(DUSTY_DIR)Dusty2DQuadIOMultiBlock.cc \
            $(DUSTY_DIR)Dusty2DQuadTurbulenceSingleBlock.cc \
            $(DUSTY_DIR)Dusty2DQuadTurbulenceMultiBlock.cc \
            $(DUSTY_DIR)Dusty2DQuadGrid.cc \
            $(DUSTY_DIR)Dusty2DQuadSolvers.cc
SRC_DUSTY_EXE = $(DUSTY_DIR)Dusty2D.cc  

# LESPremixed2D:
SRC_LESPREMIXED = $(LESPREMIXED_DIR)LESPremixed2DState.cc \
	          $(LESPREMIXED_DIR)LESPremixed2DInput.cc \
	          $(LESPREMIXED_DIR)LESPremixed2DQuad.cc \
	          $(LESPREMIXED_DIR)LESPremixed2DQuadSingleBlock.cc \
	          $(LESPREMIXED_DIR)LESPremixed2DQuadMultiBlock.cc \
	          $(LESPREMIXED_DIR)LESPremixed2DQuadGrid.cc \
	          $(LESPREMIXED_DIR)LESPremixed2DQuadSolvers.cc \
	          $(LESPREMIXED_DIR)LESPremixed2DTools.cc \
	          $(LESPREMIXED_DIR)LESPremixed2DTurbInit.cc \
	          $(LESPREMIXED_DIR)LESPremixed2DdRdU.cc \
	          $(LESPREMIXED_DIR)PowerLaw.cc \
	          $(LESPREMIXED_DIR)Scalars.cc 
SRC_LESPREMIXED_EXE = $(LESPREMIXED_DIR)LESPremixed2D.cc

# Gaussian2D:
SRC_GAUSSIAN = $(GAUSSIAN_DIR)Gaussian2DState.cc \
               $(GAUSSIAN_DIR)Gaussian2DInput.cc \
               $(GAUSSIAN_DIR)Gaussian2DCartesian.cc \
               $(GAUSSIAN_DIR)Gaussian2DCartesianSolvers.cc \
               $(GAUSSIAN_DIR)Gaussian2DQuadGrid.cc \
               $(GAUSSIAN_DIR)Gaussian2DQuadSingleBlock.cc \
               $(GAUSSIAN_DIR)Gaussian2DQuadMultiBlock.cc \
               $(GAUSSIAN_DIR)Gaussian2DQuadSolvers.cc \
               $(GAUSSIAN_DIR)Gaussian2DQuadIOSingleBlock.cc \
               $(GAUSSIAN_DIR)Gaussian2DQuadIOMultiBlock.cc
SRC_GAUSSIAN_EXE = $(GAUSSIAN_DIR)Gaussian2D.cc

# Ion5Moment2D:
SRC_ION = $(ION_DIR)Ion5Moment2DState.cc \
          $(ION_DIR)Ion5Moment2DInput.cc \
          $(ION_DIR)Ion5Moment2DQuad.cc \
          $(ION_DIR)Ion5Moment2DQuadSingleBlock.cc \
          $(ION_DIR)Ion5Moment2DQuadMultiBlock.cc \
          $(ION_DIR)Ion5Moment2DQuadGrid.cc \
          $(ION_DIR)Ion5Moment2DQuadSolvers.cc
SRC_ION_EXE = $(ION_DIR)Ion5Moment2D.cc

# Electrostatic2D:
SRC_ELECTROSTATIC = $(ELECTROSTATIC_DIR)Electrostatic2DState.cc \
	            $(ELECTROSTATIC_DIR)Electrostatic2DInput.cc \
                    $(ELECTROSTATIC_DIR)Electrostatic2DQuad.cc \
                    $(ELECTROSTATIC_DIR)Electrostatic2DQuadSingleBlock.cc \
                    $(ELECTROSTATIC_DIR)Electrostatic2DQuadMultiBlock.cc \
                    $(ELECTROSTATIC_DIR)Electrostatic2DQuadIOSingleBlock.cc \
                    $(ELECTROSTATIC_DIR)Electrostatic2DQuadIOMultiBlock.cc \
                    $(ELECTROSTATIC_DIR)Electrostatic2DQuadGrid.cc \
                    $(ELECTROSTATIC_DIR)Electrostatic2DQuadSolvers.cc
SRC_ELECTROSTATIC_EXE = $(ELECTROSTATIC_DIR)Electrostatic2D.cc   

# Rte2D:
SRC_RTE = $(RTE_DIR)Rte2DState.cc \
	  $(RTE_DIR)Medium2DState.cc \
          $(RTE_DIR)Rte2DInput.cc \
          $(RTE_DIR)Rte2DQuad.cc \
          $(RTE_DIR)Rte2DQuadSingleBlock.cc \
          $(RTE_DIR)Rte2DQuadMultiBlock.cc \
          $(RTE_DIR)Rte2DQuadGrid.cc \
          $(RTE_DIR)Rte2DQuadSolvers.cc \
          $(RTE_DIR)Rte2DTools.cc
SRC_RTE_EXE = $(RTE_DIR)Rte2D.cc

# EmbeddedBoundaries2D:
SRC_INTERFACE = $(INTERFACE_DIR)Interface2D.cc
SRC_EMBEDDEDBOUNDARIES_EXE = $(INTERFACE_DIR)EmbeddedBoundaries2D.cc

# LevelSet2D:
SRC_LEVELSET = $(LEVELSET_DIR)LevelSet2DInput.cc \
               $(LEVELSET_DIR)LevelSet2DQuad.cc \
               $(LEVELSET_DIR)LevelSet2DQuadSingleBlock.cc \
               $(LEVELSET_DIR)LevelSet2DQuadMultiBlock.cc \
               $(LEVELSET_DIR)LevelSet2DQuadIOSingleBlock.cc \
               $(LEVELSET_DIR)LevelSet2DQuadIOMultiBlock.cc \
               $(LEVELSET_DIR)LevelSet2DQuadInterfaceSingleBlock.cc \
               $(LEVELSET_DIR)LevelSet2DQuadInterfaceMultiBlock.cc \
               $(LEVELSET_DIR)LevelSet2DQuadHamiltonJacobiSingleBlock.cc \
               $(LEVELSET_DIR)LevelSet2DQuadHamiltonJacobiMultiBlock.cc \
               $(LEVELSET_DIR)LevelSet2DQuadEikonalSingleBlock.cc \
               $(LEVELSET_DIR)LevelSet2DQuadEikonalMultiBlock.cc \
               $(LEVELSET_DIR)LevelSet2DQuadScalarExtensionSingleBlock.cc \
               $(LEVELSET_DIR)LevelSet2DQuadScalarExtensionMultiBlock.cc \
               $(LEVELSET_DIR)LevelSet2DQuadGrid.cc \
               $(LEVELSET_DIR)LevelSet2DQuadSolvers.cc
SRC_LEVELSET_EXE = $(LEVELSET_DIR)LevelSet2D.cc

# Physics, Turbulence, Thermodynamics, and Chemistry:
SRC_TURBULENT = $(TURBULENT_DIR)SFSModelling.cc
SRC_EOS = $(EOS_DIR)tgas.f \
          $(EOS_DIR)ugas.f
SRC_NASADATA = $(NASADATA_DIR)NASARP1311dataclass.cc
SRC_REACTIONS = $(REACTIONS_DIR)Reactions.cc \
	   	$(REACTIONS_DIR)15step_GRI211.f \
	  	$(REACTIONS_DIR)15step_GRI3.f
SRC_REACTIONS_F90 = $(REACTIONS_DIR)c_15step_GRI211.f90 \
                    $(REACTIONS_DIR)c_15step_GRI3.f90
SRC_REACTIONS += $(SRC_REACTIONS_F90)
SRC_SNBCK = $(SNBCK_DIR)SNBCK.cc \
	    $(SNBCK_DIR)PlanckMean.cc

# TUT testing:
SRC_TUT = TestRun.cc \
          TestData.cc \
	  BasicTests/Sample_Test.cc \
          BasicTests/test_TestData.cc
SRC_TUT := $(addprefix $(UNIT_TESTING_DIR), $(SRC_TUT))

# System:
SRC_SYSTEM = $(SYSTEM_DIR)System_Linux.cc

# Base and CFFC sources:
SRC_BASE = $(SRC_MATH) $(SRC_CFD) $(SRC_MPI) $(SRC_GRID) $(SRC_AMR) \
           $(SRC_HIGH_ORDER_RECONSTRUCTION_COMMON) $(SRC_HIGH_ORDER_RECONSTRUCTION_2D) $(SRC_SYSTEM) 
ifeq ($(TUT_TESTING),ON) 
  SRC_BASE += $(SRC_TUT)
  SRC_SOLVERS1D += $(SRC_TUT)
  SRC_EULER1D_SOLVER += $(SRC_TUT)
endif
SRC_CFFC = $(SRC_BASE) \
           $(SRC_CFD1D_INPUT) $(SRC_SCALAR) $(SRC_HEAT) $(SRC_HYPERHEAT) $(SRC_EULER1D) $(SRC_MHD1D) \
           $(SRC_ADVECTDIFFUSE) $(SRC_EULER) $(SRC_NAVIERSTOKES) $(SRC_CHEM) $(SRC_HIGHTEMP) \
           $(SRC_DUSTY) $(SRC_LESPREMIXED) $(SRC_GAUSSIAN) $(SRC_ION) $(SRC_ELECTROSTATIC) \
           $(SRC_RTE) $(SRC_INTERFACE) $(SRC_LEVELSET) \
           $(SRC_TURBULENT) $(SRC_EOS) $(SRC_NASADATA) $(SRC_REACTIONS) $(SRC_SNBCK) $(SRC_ICEMCFD) \
           $(SRC_HIGH_ORDER_RECONSTRUCTION_1D)

# Executable sources:
SRC_MAIN = cffc2D.cc
SRC_ICEMCFD_UTIL1 = $(ICEMCFD_DIR)icemcfd2tecplot.cc 
SRC_ICEMCFD_UTIL2 = $(ICEMCFD_DIR)icemcfd2grid2dquad.cc

# Levermore 1D:
SRC_LEVERMORE1D_SOLVER = $(SRC_MATH) $(SRC_CFD) $(SRC_GRID_1D) \
                         $(SRC_SYSTEM) $(SRC_LEVERMORE1D) $(SRC_LEVERMORE1D_EXE)

# Add Unit Testing to Levermore 1D:
ifeq ($(TUT_TESTING),ON)
  SRC_LEVERMORE1D_SOLVER += $(SRC_TUT)
endif

#################################################################################
# Compiled object code files (automatically generated with pattern substitution)
#################################################################################
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

# 1D:
OBJ_SCALAR = $(SRC_SCALAR:.cc=.o)
OBJ_HEAT = $(SRC_HEAT:.cc=.o)
OBJ_HYPERHEAT = $(SRC_HYPERHEAT:.cc=.o)
OBJ_EULER1D = $(SRC_EULER1D:.cc=.o)
OBJ_LEVERMORE1D = $(SRC_LEVERMORE1D:.cc=.o)
OBJ_MHD1D = $(SRC_MHD1D:.cc=.o)
OBJ_CFD1D_INPUT := $(SRC_CFD1D_INPUT:.cc=.o)

# HighOrderReconstruction: (substitute *.o for *.cc)
OBJ_HIGH_ORDER_RECONSTRUCTION_COMMON := $(SRC_HIGH_ORDER_RECONSTRUCTION_COMMON:.cc=.o)
OBJ_HIGH_ORDER_RECONSTRUCTION_1D := $(SRC_HIGH_ORDER_RECONSTRUCTION_1D:.cc=.o)
OBJ_HIGH_ORDER_RECONSTRUCTION_2D := $(SRC_HIGH_ORDER_RECONSTRUCTION_2D:.cc=.o)

# Euler1D solver objects:
OBJ_EULER1D_SOLVER := $(patsubst %.cc,%.o, $(filter %.cc, $(SRC_EULER1D_SOLVER))) \
		      $(patsubst %.c,%.o, $(filter %.c, $(SRC_EULER1D_SOLVER))) \
		      $(patsubst %.f,%.o, $(filter %.f, $(SRC_EULER1D_SOLVER)))

# Levermore1D solver objects:
OBJ_LEVERMORE1D_SOLVER := $(patsubst %.cc,%.o, $(filter %.cc, $(SRC_LEVERMORE1D_SOLVER))) \
		          $(patsubst %.c,%.o, $(filter %.c, $(SRC_LEVERMORE1D_SOLVER))) \
		          $(patsubst %.f,%.o, $(filter %.f, $(SRC_LEVERMORE1D_SOLVER)))

# Solvers1D solver objects:
OBJ_SOLVERS1D_SOLVER := $(patsubst %.cc,%.o, $(filter %.cc, $(SRC_SOLVERS1D_SOLVER))) \
		        $(patsubst %.c,%.o, $(filter %.c, $(SRC_SOLVERS1D_SOLVER))) \
		        $(patsubst %.f,%.o, $(filter %.f, $(SRC_SOLVERS1D_SOLVER)))

# Advectdiffuse2D:
OBJ_ADVECTDIFFUSE = $(SRC_ADVECTDIFFUSE:.cc=.o)
OBJ_ADVECTDIFFUSE_EXE = $(SRC_ADVECTDIFFUSE_EXE:.cc=.o)

# Euler2D:
OBJ_EULER = $(SRC_EULER:.cc=.o)
OBJ_EULER_EXE = $(SRC_EULER_EXE:.cc=.o)

# NavierStokes2D:
OBJ_NAVIERSTOKES = $(SRC_NAVIERSTOKES:.cc=.o)
OBJ_NAVIERSTOKES_EXE = $(SRC_NAVIERSTOKES_EXE:.cc=.o)

# Chem2D:
OBJ_CHEM = $(SRC_CHEM:.cc=.o)
OBJ_CHEM_EXE = $(SRC_CHEM_EXE:.cc=.o)

# Flame2D:
OBJ_FLAME = $(SRC_FLAME:.cc=.o)
OBJ_FLAME_EXE = $(SRC_FLAME_EXE:.cc=.o)

# HighTemp2D:
OBJ_HIGHTEMP = $(SRC_HIGHTEMP:.cc=.o)
OBJ_HIGHTEMP_EXE = $(SRC_HIGHTEMP_EXE:.cc=.o)

# Dusty2D:
OBJ_DUSTY = $(SRC_DUSTY:.cc=.o)
OBJ_DUSTY_EXE = $(SRC_DUSTY_EXE:.cc=.o)

# LESPremixed2D:
OBJ_LESPREMIXED = $(SRC_LESPREMIXED:.cc=.o)
OBJ_LESPREMIXED_EXE = $(SRC_LESPREMIXED_EXE:.cc=.o)

# Gaussian2D:
OBJ_GAUSSIAN = $(SRC_GAUSSIAN:.cc=.o)
OBJ_GAUSSIAN_EXE = $(SRC_GAUSSIAN_EXE:.cc=.o)

# Ion5Moment2D:
OBJ_ION = $(SRC_ION:.cc=.o)
OBJ_ION_EXE = $(SRC_ION_EXE:.cc=.o)

# Electrostatic2D:
OBJ_ELECTROSTATIC = $(SRC_ELECTROSTATIC:.cc=.o)
OBJ_ELECTROSTATIC_EXE = $(SRC_ELECTROSTATIC_EXE:.cc=.o)

# Rte2D:
OBJ_RTE = $(SRC_RTE:.cc=.o)
OBJ_RTE_EXE = $(SRC_RTE_EXE:.cc=.o)
ifeq ($(CHEM2D_WITH_RTE),ON)
  OBJ_CHEM += $(OBJ_RTE)
endif
ifeq ($(FLAME2D_WITH_RTE),ON)
  OBJ_FLAME += $(OBJ_RTE)
endif

# EmbeddedBoundaries2D:
OBJ_INTERFACE = $(SRC_INTERFACE:.cc=.o)
OBJ_EMBEDDEDBOUNDARIES_EXE = $(SRC_EMBEDDEDBOUNDARIES_EXE:.cc=.o)

# LevelSet2D:
OBJ_LEVELSET = $(SRC_LEVELSET:.cc=.o)
OBJ_LEVELSET_EXE = $(SRC_LEVELSET_EXE:.cc=.o)

# Physics, Turbulence, Thermodynamics, and Chemistry:
OBJ_TURBULENT = $(SRC_TURBULENT:.cc=.o)
OBJ_EOS = $(SRC_EOS:.f=.o)
OBJ_NASADATA = $(SRC_NASADATA:.cc=.o)
OBJ_REACTIONS = $(patsubst %.cc,%.o,$(filter %.cc,$(SRC_REACTIONS)))
OBJ_REACTIONS += $(patsubst %.f,%.o,$(filter %.f,$(SRC_REACTIONS)))
OBJ_REACTIONS += $(patsubst %.f90,%.o,$(filter %.f90,$(SRC_REACTIONS)))
OBJ_REACTIONS_F90 = $(SRC_REACTIONS_F90:.f90=.o)
OBJ_SNBCK = $(SRC_SNBCK:.cc=.o)

# TUT testing:
OBJ_TUT := $(SRC_TUT:.cc=.o)

# System:
OBJ_SYSTEM = $(SRC_SYSTEM:.cc=.o)

# Base and CFFC objects:
OBJ_BASE = $(patsubst %.cc,%.o,$(filter %.cc,$(SRC_BASE)))
OBJ_BASE += $(patsubst %.c,%.o,$(filter %.c,$(SRC_BASE)))
OBJ_BASE += $(patsubst %.f90,%.o,$(filter %.f90,$(SRC_BASE)))
OBJ_BASE += $(patsubst %.f,%.o,$(filter %.f,$(SRC_BASE)))

OBJ_CFFC = $(patsubst %.cc,%.o,$(filter %.cc,$(SRC_CFFC)))
OBJ_CFFC += $(patsubst %.c,%.o,$(filter %.c,$(SRC_CFFC)))
OBJ_CFFC += $(patsubst %.f90,%.o,$(filter %.f90,$(SRC_CFFC)))
OBJ_CFFC += $(patsubst %.f,%.o,$(filter %.f,$(SRC_CFFC)))

# Executable objects:
OBJ_MAIN = $(SRC_MAIN:.cc=.o)
OBJ_ICEMCFD_UTIL1 = $(SRC_ICEMCFD_UTIL1:.cc=.o) \
                    $(OBJ_MATH) $(OBJ_CFD) $(OBJ_MPI) $(OBJ_GRID) $(OBJ_ICEMCFD) \
		    $(OBJ_HIGH_ORDER_RECONSTRUCTION_2D) $(OBJ_HIGH_ORDER_RECONSTRUCTION_COMMON)
OBJ_ICEMCFD_UTIL2 = $(SRC_ICEMCFD_UTIL2:.cc=.o) \
                    $(OBJ_MATH) $(OBJ_CFD) $(OBJ_MPI) $(OBJ_GRID) $(OBJ_ICEMCFD) \
		    $(OBJ_HIGH_ORDER_RECONSTRUCTION_2D) $(OBJ_HIGH_ORDER_RECONSTRUCTION_COMMON)
ifeq ($(TUT_TESTING),ON) 
  OBJ_ICEMCFD_UTIL1 += $(OBJ_TUT)
  OBJ_ICEMCFD_UTIL2 += $(OBJ_TUT)
endif

###############################################################################
# Define the executables
###############################################################################
EXE_MAIN = cffc2D
EXE_ICEMCFD_UTIL1 = icemcfd2tecplot
EXE_ICEMCFD_UTIL2 = icemcfd2grid2dquad
EXE_EULER1D = euler1D
EXE_SOLVERS1D = solvers1D
EXE_LEVERMORE1D = levermore1D
EXE_ADVECTDIFFUSE2D = advectdiffuse2D
EXE_EULER2D = euler2D
EXE_NAVIERSTOKES2D = navierstokes2D
EXE_CHEM2D = chem2D
EXE_FLAME2D = flame2D
EXE_HIGHTEMP2D = hightemp2D
EXE_DUSTY2D = dusty2D
EXE_LESPREMIXED2D = lespremixed2D
EXE_GAUSSIAN2D = gaussian2D
EXE_ION5MOMENT2D = ion5moment2D
EXE_ELECTROSTATIC2D = electrostatic2D
EXE_RTE2D = rte2D
EXE_LEVELSET2D = levelset2D
EXE_EMBEDDEDBOUNDARIES2D = embeddedboundaries2D

###############################################################################
# Define dependencies
###############################################################################
help:
	@echo "+-----------------------------------------------------------------------------------+"
	@echo "|                                                                                   |"
	@echo "|                           makefile for CFFC                                       |"
	@echo "|                                                                                   |"
	@echo "| Usage: make all                  compile all executable programs                  |"
	@echo "|        make cffc2D               compile cffc2D program                           |"
	@echo "|        make install              copies cffc2D executable to bin directory        |"
	@echo "|        make install-all          copies all executables to bin directory          |"
	@echo "|        make uninstall            removes cffc2D executable from bin directory     |"
	@echo "|        make uinstall-all         removes all executables from bin directory       |"
	@echo "|        make documentation        create code documentaion using doxygen & latex   |"
	@echo "|        make icemcfd2tecplot      compile icemcfd2tecplot program                  |"
	@echo "|        make icemcfd2grid2dquad   compile icemcfd2grid2dquad program               |"
	@echo "|        make euler1D              compile euler1D stand alone program              |"
	@echo "|        make euler1D-new          compile euler1D program from scratch             |"
	@echo "|        make levermore1D          compile levermore1D stand alone program          |"
	@echo "|        make solvers1D            compile solvers1D stand alone program            |"
	@echo "|        make solvers1D-new        compile solvers1D program from scratch           |"
	@echo "|        make advectdiffuse2D      compile advectdiffuse2D stand alone program      |"
	@echo "|        make euler2D              compile euler2D stand alone program              |"
	@echo "|        make navierstokes2D       compile navierstokes2D stand alone program       |"
	@echo "|        make chem2D               compile chem2D stand alone program               |"
	@echo "|        make flame2D              compile flame2D stand alone program              |"
	@echo "|        make hightemp2D           compile hightemp2D stand alone program           |"
	@echo "|        make dusty2D              compile dusty2D stand alone program              |"
	@echo "|        make lespremixed2D        compile lespremixed2D stand alone program        |"
	@echo "|        make gaussian2D           compile gaussian2D stand alone program           |"
	@echo "|        make ion5moment2D         compile ion5moment2D stand alone program         |"
	@echo "|        make electrostatic2D      compile electrostatic2D stand alone program      |"
	@echo "|        make rte2D                compile rte2D stand alone program                |"
	@echo "|        make levelset2D           compile levelset2D stand alone program           |"
	@echo "|        make embeddedboundaries2D compile embeddedboundaries2D stand alone program |"
	@echo "|        make libraries-all        create all supporting libraries                  |"
	@echo "|        make sparselib            create SparseLib++ library                       |"
	@echo "|        make bpkit                create BPKIT library                             |"
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
	@echo "|                  DO_CHECKS = ON (default), OFF                                    |"
	@echo "|                  NO_ASSERT = OFF (default), ON                                    |"
	@echo "|                                                                                   |"
	@echo "| Ensure that the source directory is clean and the ICEMCFD, SparseLib++, BPKIT,    |"
	@echo "| FFTW, and Cantera libraries are up to date before re-making CFFC from scratch.    |"
	@echo "|                                                                                   |"
	@echo "| An environment variable, CFFC_Path, needs to be set to the location               |"
	@echo "| of the CFFC source directory.  This can be done as follows:                       |"
	@echo "|                                                                                   |"
	@echo "|  bash:   export CFFC_Path=/your/path/here/CFFC                                    |"
	@echo "|  tcsh:   setenv CFFC_Path /your/path/here/CFFC                                    |"
	@echo "|                                                                                   |"
	@echo "+-----------------------------------------------------------------------------------+"

all: $(EXE_MAIN) $(EXE_ICEMCFD_UTIL1) $(EXE_ICEMCFD_UTIL2) \
     $(EXE_EULER1D) $(EXE_SOLVERS1D) \
     $(EXE_ADVECTDIFFUSE2D) $(EXE_EULER2D) $(EXE_NAVIERSTOKES2D) $(EXE_CHEM2D) \
     $(EXE_HIGHTEMP2D) $(EXE_DUSTY2D) $(EXE_LESPREMIXED2D) $(EXE_GAUSSIAN2D) \
     $(EXE_ION5MOMENT2D) $(EXE_ELECTROSTATIC2D) $(EXE_RTE2D) $(EXE_LEVELSET2D) \
     $(EXE_EMBEDDEDBOUNDARIES2D) $(EXE_LEVERMORE1D)
	@echo ' '
	@echo All programs have been brought up to date.
	@echo '------------------------------------------'
	@echo ' '

libraries-all: sparselib bpkit fftw
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

$(EXE_EULER1D): $(OBJ_EULER1D_SOLVER)
	@echo ' '
	@echo Building program $(EXE_EULER1D) for a $(PLATFORM) platform.
	@echo '-----------------------------------------------------------'
	@echo ' '
	$(LD) -o $(EXE_EULER1D) $(OBJ_EULER1D_SOLVER) $(LDFLAGS)
	@echo ' '
	@echo Program $(EXE_EULER1D) has been brought up to date.
	@echo '---------------------------------------------------'
	@echo ' '

$(EXE_EULER1D)-new:
	@$(MAKE) --no-print-directory clean
	@$(MAKE) --no-print-directory libraries-all 
	@$(MAKE) --no-print-directory $(EXE_EULER1D)


$(EXE_SOLVERS1D): $(OBJ_SOLVERS1D_SOLVER)
	@echo ' '
	@echo Building program $(EXE_SOLVERS1D) for a $(PLATFORM) platform.
	@echo '----------------------------------------------------------'
	@echo ' '
	$(LD) -o $(EXE_SOLVERS1D) $(OBJ_SOLVERS1D_SOLVER) $(LDFLAGS)
	@echo ' '
	@echo Program $(EXE_SOLVERS1D) has been brought up to date.
	@echo '---------------------------------------------------'
	@echo ' '

$(EXE_SOLVERS1D)-new:
	@$(MAKE) --no-print-directory clean
	@$(MAKE) --no-print-directory libraries-all 
	@$(MAKE) --no-print-directory $(EXE_SOLVERS1D)


$(EXE_LEVERMORE1D): $(OBJ_LEVERMORE1D_SOLVER)
	@echo ' '
	@echo Building program $(EXE_LEVERMORE1D) for a $(PLATFORM) platform.
	@echo '-----------------------------------------------------------'
	@echo ' '
	$(LD) -o $(EXE_LEVERMORE1D) $(OBJ_LEVERMORE1D_SOLVER) $(LDFLAGS)
	@echo ' '
	@echo Program $(EXE_LEVERMORE1D) has been brought up to date.
	@echo '---------------------------------------------------'
	@echo ' '

$(EXE_ADVECTDIFFUSE2D): $(OBJ_ADVECTDIFFUSE_EXE) $(OBJ_BASE) $(OBJ_ADVECTDIFFUSE) $(OBJ_ICEMCFD)
	@echo ' '
	@echo Building program $(EXE_ADVECTDIFFUSE2D) for a $(PLATFORM) platform.
	@echo '----------------------------------------------------------'
	@echo ' '
	$(LD) -o $(EXE_ADVECTDIFFUSE2D) $(OBJ_ADVECTDIFFUSE_EXE) $(OBJ_BASE) $(OBJ_ADVECTDIFFUSE) $(OBJ_ICEMCFD) $(LDFLAGS)
	@echo ' '
	@echo Program $(EXE_ADVECTDIFFUSE2D) has been brought up to date.
	@echo '--------------------------------------------------'
	@echo ' '

$(EXE_EULER2D): $(OBJ_EULER_EXE) $(OBJ_BASE) $(OBJ_EULER) $(OBJ_INTERFACE) $(OBJ_ICEMCFD)
	@echo ' '
	@echo Building program $(EXE_EULER2D) for a $(PLATFORM) platform.
	@echo '-----------------------------------------------------------'
	@echo ' '
	$(LD) -o $(EXE_EULER2D) $(OBJ_EULER_EXE) $(OBJ_BASE) $(OBJ_EULER) $(OBJ_INTERFACE) $(OBJ_ICEMCFD) $(LDFLAGS)
	@echo ' '
	@echo Program $(EXE_EULER2D) has been brought up to date.
	@echo '---------------------------------------------------'
	@echo ' '

$(EXE_NAVIERSTOKES2D): $(OBJ_NAVIERSTOKES_EXE) $(OBJ_BASE) $(OBJ_NAVIERSTOKES) $(OBJ_INTERFACE) $(OBJ_ICEMCFD)
	@echo ' '
	@echo Building program $(EXE_NAVIERSTOKES2D) for a $(PLATFORM) platform.
	@echo '------------------------------------------------------------------'
	@echo ' '
	$(LD) -o $(EXE_NAVIERSTOKES2D) $(OBJ_NAVIERSTOKES_EXE) $(OBJ_BASE) $(OBJ_NAVIERSTOKES) $(OBJ_INTERFACE) $(OBJ_ICEMCFD) $(LDFLAGS)
	@echo ' '
	@echo Program $(EXE_NAVIERSTOKES2D) has been brought up to date.
	@echo '----------------------------------------------------------'
	@echo ' '

$(EXE_CHEM2D): $(OBJ_CHEM_EXE) $(OBJ_BASE) $(OBJ_NASADATA) $(OBJ_REACTIONS) $(OBJ_SNBCK) $(OBJ_CHEM) $(OBJ_ICEMCFD)
	@echo ' '
	@echo Building program $(EXE_CHEM2D) for a $(PLATFORM) platform.
	@echo '----------------------------------------------------------'
	@echo ' '
	$(LD) -o $(EXE_CHEM2D) $(OBJ_CHEM_EXE) $(OBJ_BASE) $(OBJ_NASADATA) $(OBJ_REACTIONS) $(OBJ_SNBCK) $(OBJ_CHEM) $(OBJ_ICEMCFD) $(LDFLAGS)
	@echo ' '
	@echo Program $(EXE_CHEM2D) has been brought up to date.
	@echo '--------------------------------------------------'
	@echo ' '

$(EXE_FLAME2D): $(OBJ_FLAME_EXE) $(OBJ_BASE) $(OBJ_NASADATA) $(OBJ_SNBCK) $(OBJ_FLAME) $(OBJ_ICEMCFD)
	@echo ' '
	@echo Building program $(EXE_FLAME2D) for a $(PLATFORM) platform.
	@echo '----------------------------------------------------------'
	@echo ' '
	$(LD) -o $(EXE_FLAME2D) $(OBJ_FLAME_EXE) $(OBJ_BASE) $(OBJ_NASADATA) $(OBJ_SNBCK) $(OBJ_FLAME) $(OBJ_ICEMCFD) $(LDFLAGS)
	@echo ' '
	@echo Program $(EXE_FLAME2D) has been brought up to date.
	@echo '--------------------------------------------------'
	@echo ' '

$(EXE_HIGHTEMP2D): $(OBJ_HIGHTEMP_EXE) $(OBJ_BASE) $(OBJ_INTERFACE) $(OBJ_HIGHTEMP) $(OBJ_ICEMCFD) $(OBJ_EOS)
	@echo ' '
	@echo Building program $(EXE_HIGHTEMP2D) for a $(PLATFORM) platform.
	@echo '----------------------------------------------------------'
	@echo ' '
	$(LD) -o $(EXE_HIGHTEMP2D) $(OBJ_HIGHTEMP_EXE) $(OBJ_BASE) $(OBJ_INTERFACE) $(OBJ_HIGHTEMP) $(OBJ_ICEMCFD) $(LDFLAGS) $(OBJ_EOS)
	@echo ' '
	@echo Program $(EXE_HIGHTEMP2D) has been brought up to date.
	@echo '--------------------------------------------------'
	@echo ' '

$(EXE_DUSTY2D): $(OBJ_DUSTY_EXE) $(OBJ_BASE) $(OBJ_DUSTY) $(OBJ_INTERFACE) $(OBJ_ELECTROSTATIC) $(OBJ_ICEMCFD)
	@echo ' '
	@echo Building program $(EXE_DUSTY2D) for a $(PLATFORM) platform.
	@echo '-----------------------------------------------------------'
	@echo ' '
	$(LD) -o $(EXE_DUSTY2D) $(OBJ_DUSTY_EXE) $(OBJ_BASE) $(OBJ_DUSTY) $(OBJ_INTERFACE) $(OBJ_ELECTROSTATIC) $(OBJ_ICEMCFD) $(LDFLAGS)
	@echo ' '
	@echo Program $(EXE_DUSTY2D) has been brought up to date.
	@echo '---------------------------------------------------'
	@echo ' '

$(EXE_LESPREMIXED2D): $(OBJ_LESPREMIXED_EXE) $(OBJ_BASE) $(OBJ_NASADATA) $(OBJ_REACTIONS) $(OBJ_LESPREMIXED) $(OBJ_TURBULENT) $(OBJ_ICEMCFD) 
	@echo ' '
	@echo Building program $(EXE_LESPREMIXED2D) for a $(PLATFORM) platform.
	@echo '----------------------------------------------------------'
	@echo ' '
	$(LD) -o $(EXE_LESPREMIXED2D) $(OBJ_LESPREMIXED_EXE) $(OBJ_BASE) $(OBJ_NASADATA) $(OBJ_REACTIONS) $(OBJ_LESPREMIXED) $(OBJ_TURBULENT) $(OBJ_ICEMCFD) $(LDFLAGS)
	@echo ' '
	@echo Program $(EXE_LESPREMIXED2D) has been brought up to date.
	@echo '--------------------------------------------------'
	@echo ' '

$(EXE_GAUSSIAN2D): $(OBJ_GAUSSIAN_EXE) $(OBJ_BASE) $(OBJ_GAUSSIAN) $(OBJ_INTERFACE) $(OBJ_ICEMCFD)
	@echo ' '
	@echo Building program $(EXE_GAUSSIAN2D) for a $(PLATFORM) platform.
	@echo '--------------------------------------------------------'
	@echo ' '
	$(LD) -o $(EXE_GAUSSIAN2D) $(OBJ_GAUSSIAN_EXE) $(OBJ_BASE) $(OBJ_GAUSSIAN) $(OBJ_INTERFACE) $(OBJ_ICEMCFD) $(LDFLAGS)
	@echo ' '
	@echo Program $(EXE_GAUSSIAN2D) has been brought up to date.
	@echo '------------------------------------------------'
	@echo ' '

$(EXE_ION5MOMENT2D): $(OBJ_ION_EXE) $(OBJ_BASE) $(OBJ_ION) $(OBJ_EULER) $(OBJ_INTERFACE) $(OBJ_ICEMCFD)
	@echo ' '
	@echo Building program $(EXE_ION5MOMENT2D) for a $(PLATFORM) platform.
	@echo '----------------------------------------------------------------'
	@echo ' '
	$(LD) -o $(EXE_ION5MOMENT2D) $(OBJ_ION_EXE) $(OBJ_BASE) $(OBJ_ION) $(OBJ_EULER) $(OBJ_INTERFACE) $(OBJ_ICEMCFD) $(LDFLAGS)
	@echo ' '
	@echo Program $(EXE_ION5MOMENT2D) has been brought up to date.
	@echo '--------------------------------------------------------'
	@echo ' '

$(EXE_ELECTROSTATIC2D): $(OBJ_ELECTROSTATIC_EXE) $(OBJ_BASE) $(OBJ_ELECTROSTATIC) $(OBJ_ICEMCFD)
	@echo ' '
	@echo Building program $(EXE_ELECTROSTATIC2D) for a $(PLATFORM) platform.
	@echo '-----------------------------------------------------------'
	@echo ' '
	$(LD) -o $(EXE_ELECTROSTATIC2D) $(OBJ_ELECTROSTATIC_EXE) $(OBJ_BASE) $(OBJ_ELECTROSTATIC) $(OBJ_ICEMCFD) $(LDFLAGS)
	@echo ' '
	@echo Program $(EXE_ELECTROSTATIC2D) has been brought up to date.
	@echo '---------------------------------------------------'
	@echo ' '

$(EXE_RTE2D): $(OBJ_RTE_EXE) $(OBJ_BASE) $(OBJ_SNBCK) $(OBJ_RTE) $(OBJ_ICEMCFD)
	@echo ' '
	@echo Building program $(EXE_RTE2D) for a $(PLATFORM) platform.
	@echo '-----------------------------------------------------------'
	@echo ' '
	$(LD) -o $(EXE_RTE2D) $(OBJ_RTE_EXE) $(OBJ_BASE) $(OBJ_SNBCK) $(OBJ_RTE) $(OBJ_ICEMCFD) $(LDFLAGS)
	@echo ' '
	@echo Program $(EXE_RTE2D) has been brought up to date.
	@echo '---------------------------------------------------'
	@echo ' '

$(EXE_LEVELSET2D): $(OBJ_LEVELSET_EXE) $(OBJ_BASE) $(OBJ_LEVELSET) $(OBJ_INTERFACE) $(OBJ_ICEMCFD)
	@echo ' '
	@echo Building program $(EXE_LEVELSET2D) for a $(PLATFORM) platform.
	@echo '------------------------------------------------------------------------'
	@echo ' '
	$(LD) -o $(EXE_LEVELSET2D) $(OBJ_LEVELSET_EXE) $(OBJ_BASE) $(OBJ_LEVELSET) $(OBJ_INTERFACE) $(OBJ_ICEMCFD) $(LDFLAGS)
	@echo ' '
	@echo Program $(EXE_LEVELSET2D) has been brought up to date.
	@echo '------------------------------------------------------------------------'
	@echo ' '


$(EXE_EMBEDDEDBOUNDARIES2D): $(OBJ_EMBEDDEDBOUNDARIES_EXE) $(OBJ_BASE) $(OBJ_EMBEDDEDBOUNDARIES) $(OBJ_INTERFACE) $(OBJ_LEVELSET) $(OBJ_EULER) $(OBJ_NAVIERSTOKES) $(OBJ_DUSTY) $(OBJ_ELECTROSTATIC) $(OBJ_ICEMCFD) $(OBJ_GAUSSIAN)
	@echo ' '
	@echo Building program $(EXE_EMBEDDEDBOUNDARIES2D) for a $(PLATFORM) platform.
	@echo '------------------------------------------------------------------'
	@echo ' '
	$(LD) -o $(EXE_EMBEDDEDBOUNDARIES2D) $(OBJ_EMBEDDEDBOUNDARIES_EXE) $(OBJ_BASE) $(OBJ_EMBEDDEDBOUNDARIES) $(OBJ_INTERFACE) $(OBJ_LEVELSET) $(OBJ_EULER) $(OBJ_NAVIERSTOKES) $(OBJ_DUSTY) $(OBJ_ELECTROSTATIC) $(OBJ_ICEMCFD) $(OBJ_GAUSSIAN) $(LDFLAGS)
	@echo ' '
	@echo Program $(EXE_EMBEDDEDBOUNDARIES2D) has been brought up to date.
	@echo '----------------------------------------------------------'
	@echo ' '

$(EXE_ICEMCFD_UTIL1): $(OBJ_ICEMCFD_UTIL1)
	@echo ' '
	@echo Building program $(EXE_ICEMCFD_UTIL1) for a $(PLATFORM) platform.
	@echo '-----------------------------------------------------------------'
	@echo ' '
	$(LD) -o $(EXE_ICEMCFD_UTIL1) $(OBJ_ICEMCFD_UTIL1) $(LDFLAGS)
	@echo ' '
	@echo Program $(EXE_ICEMCFD_UTIL1) has been brought up to date.
	@echo '---------------------------------------------------------'
	@echo ' '

$(EXE_ICEMCFD_UTIL2): $(OBJ_ICEMCFD_UTIL2)
	@echo ' '
	@echo Building program $(EXE_ICEMCFD_UTIL2) for a $(PLATFORM) platform.
	@echo '-----------------------------------------------------------------'
	@echo ' '
	$(LD) -o $(EXE_ICEMCFD_UTIL2) $(OBJ_ICEMCFD_UTIL2) $(LDFLAGS)
	@echo ' '
	@echo Program $(EXE_ICEMCFD_UTIL2) has been brought up to date.
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

install-all: $(EXE_MAIN) $(EXE_ICEMCFD_UTIL1) $(EXE_ICEMCFD_UTIL2) $(EXE_EULER1D) $(EXE_SOLVERS1D) \
        $(EXE_ADVECTDIFFUSE2D) $(EXE_EULER2D) $(EXE_NAVIERSTOKES2D) $(EXE_CHEM2D) \
        $(EXE_HIGHTEMP2D) $(EXE_DUSTY2D) $(EXE_LESPREMIXED2D) $(EXE_GAUSSIAN2D) \
        $(EXE_ION5MOMENT2D) $(EXE_ELECTROSTATIC2D) $(EXE_RTE2D) $(EXE_LEVELSET2D) \
        $(EXE_EMBEDDEDBOUNDARIES2D) $(EXE_LEVERMORE1D)
	@for f in $(EXE_MAIN) $(EXE_ICEMCFD_UTIL1) $(EXE_ICEMCFD_UTIL2) $(EXE_EULER1D) $(EXE_SOLVERS1D) \
        $(EXE_ADVECTDIFFUSE2D) $(EXE_EULER2D) $(EXE_NAVIERSTOKES2D) $(EXE_CHEM2D) \
        $(EXE_HIGHTEMP2D) $(EXE_DUSTY2D) $(EXE_LESPREMIXED2D) $(EXE_GAUSSIAN2D) \
        $(EXE_ION5MOMENT2D) $(EXE_ELECTROSTATIC2D) $(EXE_RTE2D) $(EXE_LEVELSET2D) \
        $(EXE_EMBEDDEDBOUNDARIES2D) $(EXE_LEVERMORE1D); \
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
	@for f in $(EXE_MAIN) $(EXE_ICEMCFD_UTIL1) $(EXE_ICEMCFD_UTIL2) $(EXE_EULER1D) $(EXE_SOLVERS1D) \
        $(EXE_ADVECTDIFFUSE2D) $(EXE_EULER2D) $(EXE_NAVIERSTOKES2D) $(EXE_CHEM2D) \
        $(EXE_HIGHTEMP2D) $(EXE_DUSTY2D) $(EXE_LESPREMIXED2D) $(EXE_GAUSSIAN2D) \
        $(EXE_ION5MOMENT2D) $(EXE_ELECTROSTATIC2D) $(EXE_RTE2D) $(EXE_LEVELSET2D) \
        $(EXE_EMBEDDEDBOUNDARIES2D) $(EXE_LEVERMORE1D); \
	do \
	  echo "Uninstalling $$f"; \
	  cd $(CFFC_BIN_DIR); \
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
	cp $(DOXYGEN_DIR_2D)/code_documentation.doxygen .
	doxygen code_documentation.doxygen
	rm code_documentation.doxygen
	@echo ' '
	@echo Code documentation created using doxygen.
	@echo '-----------------------------------------'
	@echo ' '
	@echo Creating user manuals using LaTeX.
	@echo '------------------------------------------'
	@echo ' '
	cd $(LATEX_DIR_2D)/users_guide_src; make.userguide;
	cd $(LATEX_DIR_2D)/programmers_guide_src; make.programmersguide; 
	cd $(LATEX_DIR_2D)/reference_manual_src; make.referencemanual;
	@echo ' '
	@echo User manuals created using LaTeX.
	@echo '-----------------------------------------'
	@echo ' '

#
# Clean up source, documentation, and libraray directories
clean:
	@echo ' '
	@echo Cleaning source file directories.
	@echo '---------------------------------'
	@echo ' '
	rm -f *.depend *.o *.bak *.il *~
	rm -f journal.tut log.tut
	rm -f gmon.out tecplot.phy
	rm -f $(MATH_DIR)*.o $(MATH_DIR)*.il $(MATH_DIR)*~
	rm -f $(MATH_DIR)Polyfit/*.o $(MATH_DIR)Polyfit/*.il $(MATH_DIR)Polyfit/*~
	rm -f $(MATH_DIR)UnitTests/*.o $(MATH_DIR)UnitTests/*.il $(MATH_DIR)UnitTests/*~
	rm -f $(MATH_DIR)UnitTests/TestFunctions/*.o $(MATH_DIR)UnitTests/TestFunctions/*.il $(MATH_DIR)UnitTests/TestFunctions/*~
	rm -f $(CFD_DIR)*.o $(CFD_DIR)*.il $(CFD_DIR)*~
	rm -f $(MPI_DIR)*.o $(MPI_DIR)*.il $(MPI_DIR)*~
	rm -f $(GRID_DIR)*.o $(GRID_DIR)*.il $(GRID_DIR)*~
	rm -f $(ICEMCFD_DIR)*.o $(ICEMCFD_DIR)*.il $(ICEMCFD_DIR)*~ $(ICEMCFD_DIR)*.depend $(ICEMCFD_DIR)*.bak
	rm -f $(GRID_DIR)UnitTests/*.o $(GRID_DIR)UnitTests/*.il $(GRID_DIR)UnitTests/*~ $(GRID_DIR)UnitTests/*.bak
	rm -f $(AMR_DIR)*.o $(AMR_DIR)*.il $(AMR_DIR)*~
	rm -f $(SYSTEM_DIR)*.o $(SYSTEM_DIR)*.il $(SYSTEM_DIR)*~
	rm -f $(EULER_1D_DIR)*.o $(EULER_1D_DIR)*.il $(EULER_1D_DIR)*~ $(EULER_1D_DIR)*.bak
	rm -f $(LEVERMORE_1D_DIR)*.o $(LEVERMORE_1D_DIR)*.il $(LEVERMORE_1D_DIR)*~ $(LEVERMORE_1D_DIR)*.bak
	rm -f $(LEVERMORE_1D_DIR)UnitTests/*.o $(LEVERMORE_1D_DIR)UnitTests/*.il $(LEVERMORE_1D_DIR)UnitTests/*~ $(LEVERMORE_1D_DIR)UnitTests/*.bak
	rm -f $(EULER_1D_DIR)ExactSolutions/*.o $(EULER_1D_DIR)ExactSolutions/*.il $(EULER_1D_DIR)ExactSolutions/*~ $(EULER_1D_DIR)ExactSolutions/*.bak
	rm -f $(EULER_1D_DIR)UnitTests/*.o $(EULER_1D_DIR)UnitTests/*.il $(EULER_1D_DIR)UnitTests/*~ $(EULER_1D_DIR)UnitTests/*.bak
	rm -f $(SOLVERS1D_DIR)*.o $(SOLVERS1D_DIR)*.il $(SOLVERS1D_DIR)*~ $(SOLVERS1D_DIR)*.depend
	rm -f $(ADVECTDIFFUSE_DIR)*.o $(ADVECTDIFFUSE_DIR)*.il $(ADVECTDIFFUSE_DIR)*~ $(ADVECTDIFFUSE_DIR)*.depend $(ADVECTDIFFUSE_DIR)*.bak
	rm -f $(ADVECTDIFFUSE_DIR)UnitTests/*.o $(ADVECTDIFFUSE_DIR)UnitTests/*.il $(ADVECTDIFFUSE_DIR)UnitTests/*~ $(ADVECTDIFFUSE_DIR)UnitTests/*.bak
	rm -f $(EULER_2D_DIR)*.o $(EULER_2D_DIR)*.il $(EULER_2D_DIR)*~ $(EULER_2D_DIR)*.depend $(EULER_2D_DIR)*.bak
	rm -f $(NAVIERSTOKES_DIR)*.o $(NAVIERSTOKES_DIR)*.il $(NAVIERSTOKES_DIR)*~ $(NAVIERSTOKES_DIR)*.depend $(NAVIERSTOKES_DIR)*.bak
	rm -f $(CHEM_DIR)*.o $(CHEM_DIR)*.il $(CHEM_DIR)*~ $(CHEM_DIR)*.depend $(CHEM_DIR)*.bak
	rm -f $(FLAME_DIR)*.o $(FLAME_DIR)*.il $(FLAME_DIR)*~ $(FLAME_DIR)*.depend $(FLAME_DIR)*.bak
	rm -f $(HIGHTEMP_DIR)*.o $(HIGHTEMP_DIR)*.il $(HIGHTEMP_DIR)*~ $(HIGHTEMP_DIR)*.depend $(HIGHTEMP_DIR)*.bak
	rm -f $(DUSTY_DIR)*.o $(DUSTY_DIR)*.il $(DUSTY_DIR)*~ $(DUSTY_DIR)*.depend $(DUSTY_DIR)*.bak
	rm -f $(LESPREMIXED_DIR)*.o $(LESPREMIXED_DIR)*.il $(LESPREMIXED_DIR)*~ $(LESPREMIXED_DIR)*.depend $(LESPREMIXED_DIR)*.bak
	rm -f $(GAUSSIAN_DIR)*.o $(GAUSSIAN_DIR)*.il $(GAUSSIAN_DIR)*~ $(GAUSSIAN_DIR)*.depend $(GAUSSIAN_DIR)*.bak
	rm -f $(ION_DIR)*.o $(ION_DIR)*.il $(ION_DIR)*~ $(ION_DIR)*.depend $(ION_DIR)*.bak
	rm -f $(ELECTROSTATIC_DIR)*.o $(ELECTROSTATIC_DIR)*.il $(ELECTROSTATIC_DIR)*~ $(ELECTROSTATIC_DIR)*.depend $(ELECTROSTATIC_DIR)*.bak
	rm -f $(RTE_DIR)*.o $(RTE_DIR)*.il $(RTE_DIR)*~ $(RTE_DIR)*.depend $(RTE_DIR)*.bak
	rm -f $(INTERFACE_DIR)*.o $(INTERFACE_DIR)*.il $(INTERFACE_DIR)*~ $(INTERFACE_DIR)*.depend $(INTERFACE_DIR)*.bak
	rm -f $(LEVELSET_DIR)*.o $(LEVELSET_DIR)*.il $(LEVELSET_DIR)*~ $(LEVELSET_DIR)*.depend $(LEVELSET_DIR)*.bak
	rm -f $(PHYSICS_DIR)*.o $(PHYSICS_DIR)*.il $(PHYSICS_DIR)*~
	rm -f $(TURBULENT_DIR)*.o $(TURBULENT_DIR)*.il $(TURBULENT_DIR)*~
	rm -f $(EOS_DIR)*.o $(EOS_DIR)*.il $(EOS_DIR)*~
	rm -f $(NASADATA_DIR)*.o $(NASADATA_DIR)*.il $(NASADATA_DIR)*~
	rm -f $(REACTIONS_DIR)*.o $(REACTIONS_DIR)*.il $(REACTIONS_DIR)*~
	rm -f $(SNBCK_DIR)*.o $(SNBCK_DIR)*.il $(SNBCK_DIR)*~
	rm -f $(HIGH_ORDER_RECONSTRUCTION_DIR)*.o $(HIGH_ORDER_RECONSTRUCTION_DIR)*.il $(HIGH_ORDER_RECONSTRUCTION_DIR)*~
	rm -f $(HIGH_ORDER_RECONSTRUCTION_DIR)UnitTests/*.o $(HIGH_ORDER_RECONSTRUCTION_DIR)UnitTests/*.il $(HIGH_ORDER_RECONSTRUCTION_DIR)UnitTests/*~
	rm -f $(UNIT_TESTING_DIR)*.o $(UNIT_TESTING_DIR)*.il $(UNIT_TESTING_DIR)*~
	rm -f $(UNIT_TESTING_DIR)BasicTests/*.o $(UNIT_TESTING_DIR)BasicTests/*.il $(UNIT_TESTING_DIR)BasicTests/*~
	rm -f Utilities/*.o Utilities/.il Utilities/*~
	rm -f $(MULTIGRID_DIR)*~
	rm -f $(NKS_DIR)*~
	rm -f $(MODULE_DIR)/*.mod $(MODULE_DIR)/*~
	rm -f $(EXE_MAIN)
	rm -f $(EXE_ICEMCFD_UTIL1) $(EXE_ICEMCFD_UTIL2)
	rm -f $(EXE_EULER1D)
	rm -f $(EXE_LEVERMORE1D)
	rm -f $(EXE_SOLVERS1D)
	rm -f $(EXE_ADVECTDIFFUSE2D)
	rm -f $(EXE_EULER2D)
	rm -f $(EXE_NAVIERSTOKES2D)
	rm -f $(EXE_CHEM2D)
	rm -f $(EXE_FLAME2D)
	rm -f $(EXE_HIGHTEMP2D)
	rm -f $(EXE_DUSTY2D)
	rm -f $(EXE_LESPREMIXED2D)
	rm -f $(EXE_GAUSSIAN2D)
	rm -f $(EXE_ION5MOMENT2D)
	rm -f $(EXE_ELECTROSTATIC2D)
	rm -f $(EXE_RTE2D)
	rm -f $(EXE_LEVELSET2D)
	rm -f $(EXE_EMBEDDEDBOUNDARIES2D)
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
	rm -f $(DOXYGEN_DIR_2D)/html/*
	rm -f $(DOXYGEN_DIR_2D)/latex/*
	cd $(LATEX_DIR_2D)/users_guide_src; rm -f *~ *.ps *.pdf *.aux *.log *.dvi *.toc users_guide/*; 
	cd $(LATEX_DIR_2D)/programmers_guide_src; rm -f *~ *.ps *.pdf *.aux *.log *.dvi *.toc programmers_guide/*; 
	cd $(LATEX_DIR_2D)/reference_manual_src; rm -f *~ *.ps *.pdf *.aux *.log *.dvi *.toc reference_manual/*;
	cd $(SPARSELIB_DIR); make CFFC_SRC_TREE=$(CFFC_SRC_TREE) wipe;
	cd $(SPARSELIB_DIR); make CFFC_SRC_TREE=$(CFFC_SRC_TREE) clean;
	cd $(BPKIT_DIR)/src; make CFFC_SRC_TREE=$(CFFC_SRC_TREE) clean;
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
             bpkit sparselib fftw cantera documentation \
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
	sinclude $(patsubst %.cc,$(DEPEND_DIR)%.d, $(filter %.cc, $(SRC_ICEMCFD_UTIL1)))
	sinclude $(patsubst %.cc,$(DEPEND_DIR)%.d, $(filter %.cc, $(SRC_ICEMCFD_UTIL2)))
	sinclude $(patsubst %.cc,$(DEPEND_DIR)%.d, $(filter %.cc, $(SRC_SOLVERS1D_EXE)))
	sinclude $(patsubst %.cc,$(DEPEND_DIR)%.d, $(filter %.cc, $(SRC_EULER1D_EXE)))
	sinclude $(patsubst %.cc,$(DEPEND_DIR)%.d, $(filter %.cc, $(SRC_LEVERMORE1D_EXE)))
	sinclude $(patsubst %.cc,$(DEPEND_DIR)%.d, $(filter %.cc, $(SRC_ADVECTDIFFUSE_EXE)))
	sinclude $(patsubst %.cc,$(DEPEND_DIR)%.d, $(filter %.cc, $(SRC_EULER_EXE)))
	sinclude $(patsubst %.cc,$(DEPEND_DIR)%.d, $(filter %.cc, $(SRC_NAVIERSTOKES_EXE)))
	sinclude $(patsubst %.cc,$(DEPEND_DIR)%.d, $(filter %.cc, $(SRC_CHEM_EXE)))
	sinclude $(patsubst %.cc,$(DEPEND_DIR)%.d, $(filter %.cc, $(SRC_FLAME_EXE)))
	sinclude $(patsubst %.cc,$(DEPEND_DIR)%.d, $(filter %.cc, $(SRC_HIGHTEMP_EXE)))
	sinclude $(patsubst %.cc,$(DEPEND_DIR)%.d, $(filter %.cc, $(SRC_DUSTY_EXE)))
	sinclude $(patsubst %.cc,$(DEPEND_DIR)%.d, $(filter %.cc, $(SRC_LESPREMIXED_EXE)))
	sinclude $(patsubst %.cc,$(DEPEND_DIR)%.d, $(filter %.cc, $(SRC_GAUSSIAN_EXE)))
	sinclude $(patsubst %.cc,$(DEPEND_DIR)%.d, $(filter %.cc, $(SRC_ION_EXE)))
	sinclude $(patsubst %.cc,$(DEPEND_DIR)%.d, $(filter %.cc, $(SRC_ELECTROSTATIC_EXE)))
	sinclude $(patsubst %.cc,$(DEPEND_DIR)%.d, $(filter %.cc, $(SRC_RTE_EXE)))
	sinclude $(patsubst %.cc,$(DEPEND_DIR)%.d, $(filter %.cc, $(SRC_LEVELSET_EXE)))
	sinclude $(patsubst %.cc,$(DEPEND_DIR)%.d, $(filter %.cc, $(SRC_EMBEDDEDBOUNDARIES_EXE)))
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

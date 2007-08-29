#
# Makefile for CFFC (2D version)
#
# Computational Framework for Fluids and Combustion (CFFC)
#

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
ADVECTDIFFUSE_DIR = AdvectDiffuse2D/
AMR_DIR = AMR/
CFD_DIR = CFD/
CHEM_DIR = Chem2D/
DUSTY_DIR = Dusty2D/
ELECTROSTATIC_DIR = Electrostatic2D/
EOS_DIR = EquationOfState/
EULER_DIR = Euler2D/
GRID_DIR = Grid/
HIGHTEMP_DIR = HighTemp2D/
ICEM_DIR = ICEM/
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
NASARP1311_DIR = $(PHYSICS_DIR)NASAData/
REACTIONS_DIR = Reactions/

#
# C++ Source files
SRC_MATH = $(MATH_DIR)Vector2D.cc \
           $(MATH_DIR)Tensor2D.cc \
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
	   $(MATH_DIR)complexify.f90
SRC_CFD = $(CFD_DIR)CFD.cc
SRC_MPI = $(MPI_DIR)MPI.cc
SRC_SYSTEM = $(SYSTEM_DIR)System_Linux.cc
SRC_GRID = $(GRID_DIR)Cell1D.cc \
           $(GRID_DIR)Cell2D.cc \
           $(GRID_DIR)Grid2DQuadSingleBlock.cc \
           $(GRID_DIR)Grid2DQuadMultiBlock.cc \
           $(GRID_DIR)Grid3DHexaSingleBlock.cc \
           $(GRID_DIR)Grid3DHexaMultiBlock.cc
ifeq ($(MPI_VERSION),NoMPI)
  SRC_AMR = $(AMR_DIR)AdaptiveBlock.cc \
	    $(AMR_DIR)AdaptiveBlock_NoMPI.cc \
            $(AMR_DIR)QuadTree.cc
else
  SRC_AMR = $(AMR_DIR)AdaptiveBlock.cc \
	    $(AMR_DIR)AdaptiveBlock_MPI.cc \
            $(AMR_DIR)QuadTree.cc
endif
SRC_SCALAR = $(SOLVERS1D_DIR)Scalar1D.cc \
             $(SOLVERS1D_DIR)Scalar1DSolvers.cc
SRC_HEAT = $(SOLVERS1D_DIR)Heat1D.cc \
           $(SOLVERS1D_DIR)Heat1DSolvers.cc
SRC_HYPERHEAT = $(SOLVERS1D_DIR)HyperHeat1DState.cc \
                $(SOLVERS1D_DIR)HyperHeat1D.cc \
                $(SOLVERS1D_DIR)HyperHeat1DSolvers.cc
SRC_EULER1D = $(SOLVERS1D_DIR)Euler1DState.cc \
              $(SOLVERS1D_DIR)Euler1D.cc \
              $(SOLVERS1D_DIR)Euler1DSolvers.cc
SRC_MHD1D = $(SOLVERS1D_DIR)MHD1DState.cc \
            $(SOLVERS1D_DIR)MHD1D.cc \
            $(SOLVERS1D_DIR)MHD1DSolvers.cc
SRC_ADVECTDIFFUSE = $(ADVECTDIFFUSE_DIR)AdvectDiffuse2DState.cc \
	            $(ADVECTDIFFUSE_DIR)AdvectDiffuse2DInput.cc \
                    $(ADVECTDIFFUSE_DIR)AdvectDiffuse2DQuad.cc \
	            $(ADVECTDIFFUSE_DIR)AdvectDiffuse2DQuadSingleBlock.cc \
                    $(ADVECTDIFFUSE_DIR)AdvectDiffuse2DQuadMultiBlock.cc \
                    $(ADVECTDIFFUSE_DIR)AdvectDiffuse2DQuadGrid.cc \
                    $(ADVECTDIFFUSE_DIR)AdvectDiffuse2DQuadSolvers.cc \
                    $(ADVECTDIFFUSE_DIR)AdvectDiffuse2DQuad_NKS.cc
SRC_ADVECTDIFFUSE_EXE = $(ADVECTDIFFUSE_DIR)AdvectDiffuse2D.cc
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
SRC_EOS = $(EOS_DIR)tgas.f \
          $(EOS_DIR)ugas.f
SRC_GAUSSIAN = $(GAUSSIAN_DIR)Gaussian2DState.cc $(GAUSSIAN_DIR)Gaussian2DInput.cc \
               $(GAUSSIAN_DIR)Gaussian2DCartesian.cc $(GAUSSIAN_DIR)Gaussian2DCartesianSolvers.cc \
               $(GAUSSIAN_DIR)Gaussian2DQuadGrid.cc $(GAUSSIAN_DIR)Gaussian2DQuadSingleBlock.cc \
               $(GAUSSIAN_DIR)Gaussian2DQuadMultiBlock.cc $(GAUSSIAN_DIR)Gaussian2DQuadSolvers.cc \
               $(GAUSSIAN_DIR)Gaussian2DQuadIOSingleBlock.cc $(GAUSSIAN_DIR)Gaussian2DQuadIOMultiBlock.cc
SRC_GAUSSIAN_EXE = $(GAUSSIAN_DIR)Gaussian2D.cc
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
SRC_EULER = $(EULER_DIR)Euler2DState.cc \
            $(EULER_DIR)Euler2DInput.cc \
            $(EULER_DIR)Euler2DQuad.cc \
            $(EULER_DIR)Euler2DQuadSingleBlock.cc \
            $(EULER_DIR)Euler2DQuadMultiBlock.cc \
            $(EULER_DIR)Euler2DQuadGrid.cc \
            $(EULER_DIR)Euler2DQuadSolvers.cc 
SRC_EULER_EXE = $(EULER_DIR)Euler2D.cc
SRC_INTERFACE = $(INTERFACE_DIR)Interface2D.cc
SRC_EMBEDDEDBOUNDARIES_EXE = $(INTERFACE_DIR)EmbeddedBoundaries2D.cc
SRC_ION = $(ION_DIR)Ion5Moment2DState.cc \
          $(ION_DIR)Ion5Moment2DInput.cc \
          $(ION_DIR)Ion5Moment2DQuad.cc \
          $(ION_DIR)Ion5Moment2DQuadSingleBlock.cc \
          $(ION_DIR)Ion5Moment2DQuadMultiBlock.cc \
          $(ION_DIR)Ion5Moment2DQuadGrid.cc \
          $(ION_DIR)Ion5Moment2DQuadSolvers.cc
SRC_ION_EXE = $(ION_DIR)Ion5Moment2D.cc
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
SRC_RTE = $(RTE_DIR)Rte2DState.cc \
          $(RTE_DIR)Rte2DInput.cc \
          $(RTE_DIR)Rte2DQuad.cc \
          $(RTE_DIR)Rte2DQuadSingleBlock.cc \
          $(RTE_DIR)Rte2DQuadMultiBlock.cc \
          $(RTE_DIR)Rte2DQuadGrid.cc \
          $(RTE_DIR)Rte2DQuadSolvers.cc \
          $(RTE_DIR)Rte2DTools.cc \
          $(RTE_DIR)SNBCK.cc 
SRC_RTE_EXE = $(RTE_DIR)Rte2D.cc
SRC_NAVIERSTOKES = $(NAVIERSTOKES_DIR)NavierStokes2DState.cc \
                   $(NAVIERSTOKES_DIR)NavierStokes2DInput.cc \
                   $(NAVIERSTOKES_DIR)NavierStokes2DQuad.cc \
                   $(NAVIERSTOKES_DIR)NavierStokes2DQuadSingleBlock.cc \
                   $(NAVIERSTOKES_DIR)NavierStokes2DQuadMultiBlock.cc \
                   $(NAVIERSTOKES_DIR)NavierStokes2DQuadIOSingleBlock.cc \
                   $(NAVIERSTOKES_DIR)NavierStokes2DQuadIOMultiBlock.cc \
                   $(NAVIERSTOKES_DIR)NavierStokes2DQuadTurbulenceSingleBlock.cc \
                   $(NAVIERSTOKES_DIR)NavierStokes2DQuadTurbulenceMultiBlock.cc \
                   $(NAVIERSTOKES_DIR)NavierStokes2DQuadGrid.cc \
                   $(NAVIERSTOKES_DIR)NavierStokes2DQuadSolvers.cc
SRC_NAVIERSTOKES_EXE = $(NAVIERSTOKES_DIR)NavierStokes2D.cc
SRC_NASARP1311 = $(NASARP1311_DIR)NASARP1311dataclass.cc   
SRC_TURBULENT = $(TURBULENT_DIR)SFSModelling.cc
SRC_REACTIONS = $(REACTIONS_DIR)Reactions.cc \
	   	$(REACTIONS_DIR)15step_GRI211.f \
	  	$(REACTIONS_DIR)15step_GRI3.f \
	   	$(REACTIONS_DIR)c_15step_GRI211.f90 \
	   	$(REACTIONS_DIR)c_15step_GRI3.f90

SRC_ICEMCFD = $(ICEM_DIR)ICEMCFD.cc \
	      $(ICEM_DIR)ICEMCFD_ctype.c
SRC_BASE =  $(SRC_MATH) $(SRC_CFD) $(SRC_MPI) $(SRC_GRID) $(SRC_AMR) \
            $(SRC_SYSTEM)
SRC_CFFC = $(SRC_BASE) $(SRC_SCALAR) $(SRC_HEAT) $(SRC_HYPERHEAT) $(SRC_EULER1D) $(SRC_MHD1D) \
           $(SRC_ADVECTDIFFUSE) $(SRC_CHEM) $(SRC_DUSTY) $(SRC_EULER) $(SRC_HIGHTEMP) \
	   $(SRC_ION) $(SRC_ELECTROSTATIC) $(SRC_INTERFACE) $(SRC_LEVELSET) $(SRC_NAVIERSTOKES) \
	   $(SRC_GAUSSIAN) $(SRC_LESPREMIXED) $(SRC_RTE2D) $(SRC_TURBULENT) $(SRC_EOS) $(SRC_ICEMCFD)
SRC_MAIN = cffc2D.cc
SRC_ICEMCFD_UTIL1 = $(ICEM_DIR)icemcfd2tecplot.cc 
SRC_ICEMCFD_UTIL2 = $(ICEM_DIR)icemcfd2grid2dquad.cc

#
# C++ Header files
INC_MATH  = $(MATH_DIR)Math.h \
            $(MATH_DIR)Vector2D.h \
            $(MATH_DIR)Vector3D.h \
            $(MATH_DIR)Tensor2D.h \
            $(MATH_DIR)Tensor3D.h \
            $(MATH_DIR)SplineFit.h \
            $(MATH_DIR)Spline2D.h \
            $(MATH_DIR)Matrix.h \
            $(MATH_DIR)LinearSystems.h \
	    $(MATH_DIR)LinkedList.h \
	    $(MATH_DIR)Polygon.h \
	    $(MATH_DIR)Simpson.h \
	    $(MATH_DIR)Quadrature.h \
	    $(MATH_DIR)Polyfit/Polyfit.h \
            $(MATH_DIR)Complexify.h
INC_SYSTEM = $(SYSTEM_DIR)System_Linux.h
INC_CFD = $(CFD_DIR)CFD.h
INC_MPI = $(MPI_DIR)MPI.h
INC_PHYSICS = $(PHYSICS_DIR)GasConstants.h \
	      $(PHYSICS_DIR)SolidConstants.h \
	      $(PHYSICS_DIR)Species.h 
INC_GRID = $(GRID_DIR)Cell1D.h \
	   $(GRID_DIR)Cell2D.h \
	   $(GRID_DIR)Cell3D.h \
           $(GRID_DIR)Grid2DQuad.h \
           $(GRID_DIR)Grid3DHexa.h
INC_AMR = $(AMR_DIR)AdaptiveBlock.h \
          $(AMR_DIR)AdaptiveBlock2D_MessagePassing.h \
          $(AMR_DIR)QuadTree.h \
          $(AMR_DIR)QuadTree_Morton.h \
	  $(AMR_DIR)AMR.h
INC_SCALAR = $(SOLVERS1D_DIR)Scalar1D.h
INC_HEAT = $(SOLVERS1D_DIR)Heat1D.h
INC_HYPERHEAT = $(SOLVERS1D_DIR)HyperHeat1DState.h \
                $(SOLVERS1D_DIR)HyperHeat1D.h
INC_EULER1D = $(SOLVERS1D_DIR)Euler1DState.h \
              $(SOLVERS1D_DIR)Euler1D.h
INC_MHD1D = $(SOLVERS1D_DIR)MHD1DState.h \
            $(SOLVERS1D_DIR)MHD1D.h
INC_ADVECTDIFFUSE = $(ADVECTDIFFUSE_DIR)AdvectDiffuse2DState.h \
                    $(ADVECTDIFFUSE_DIR)AdvectDiffuse2DInput.h \
	            $(ADVECTDIFFUSE_DIR)AdvectDiffuse2DQuad.h \
		    $(ADVECTDIFFUSE_DIR)AdvectDiffuse2DQuad_NKS.h \
	            $(ADVECTDIFFUSE_DIR)AdvectDiffuse2DQuad_GMRES.h \
	            $(ADVECTDIFFUSE_DIR)AdvectDiffuse2DQuadMultigrid.h
INC_CHEM = $(CHEM_DIR)Chem2DInput.h \
	   $(CHEM_DIR)Chem2DQuad.h \
	   $(CHEM_DIR)Chem2DState.h \
	   $(CHEM_DIR)Chem2DQuadMultigrid.h 
INC_HIGHTEMP = $(HIGHTEMP_DIR)HighTemp2DState.h \
               $(HIGHTEMP_DIR)HighTemp2DInput.h  \
               $(HIGHTEMP_DIR)HighTemp2DQuad.h  \
               $(HIGHTEMP_DIR)HighTemp2DdRdU.h  \
               $(HIGHTEMP_DIR)HighTemp2DQuadMultigrid.h  \
               $(HIGHTEMP_DIR)HighTemp2DQuadNKS.h
INC_EOS = $(EOS_DIR)EquationOfState.h
INC_GAUSSIAN = $(GAUSSIAN_DIR)Gaussian2DState.h \
	       $(GAUSSIAN_DIR)Gaussian2DInput.h \
               $(GAUSSIAN_DIR)Gaussian2DCartesian.h \
	       $(GAUSSIAN_DIR)Gaussian2DQuad.h
INC_DUSTY = $(DUSTY_DIR)Particle2DState.h \
	    $(DUSTY_DIR)Particle2DComponents.h \
	    $(DUSTY_DIR)Dusty2DState.h \
	    $(DUSTY_DIR)Dusty2DInput.h \
	    $(DUSTY_DIR)Dusty2DQuad.h \
	    $(DUSTY_DIR)Dusty2DQuadMultigrid.h
INC_ELECTROSTATIC = $(ELECTROSTATIC_DIR)Electrostatic2DState.h \
	            $(ELECTROSTATIC_DIR)Electrostatic2DInput.h \
	            $(ELECTROSTATIC_DIR)Electrostatic2DQuad.h
INC_EULER = $(EULER_DIR)Euler2DState.h \
	    $(EULER_DIR)Euler2DInput.h \
	    $(EULER_DIR)Euler2DQuad.h 
INC_INTERFACE = $(INTERFACE_DIR)Interface2D.h
INC_EMBEDDEDBOUNDARIES = $(INTERFACE_DIR)EmbeddedBoundaries2D_Input.h \
	                 $(INTERFACE_DIR)EmbeddedBoundaries2D.h \
	                 $(INTERFACE_DIR)EmbeddedBoundaries2D_Solver.h \
	                 $(INTERFACE_DIR)EmbeddedBoundaries2D_FASMultigrid.h \
	                 $(INTERFACE_DIR)EmbeddedBoundaries2D_Euler.h \
	                 $(INTERFACE_DIR)EmbeddedBoundaries2D_NavierStokes.h \
	                 $(INTERFACE_DIR)EmbeddedBoundaries2D_Dusty.h
INC_ION = $(ION_DIR)Ion5Moment2DState.h \
          $(ION_DIR)Ion5Moment2DInput.h \
          $(ION_DIR)Ion5Moment2DQuad.h
INC_LEVELSET = $(LEVELSET_DIR)LevelSet2DState.h \
               $(LEVELSET_DIR)LevelSet2DInput.h \
               $(LEVELSET_DIR)LevelSet2DQuad.h
INC_LESPREMIXED = $(LESPREMIXED_DIR)LESPremixed2DState.h \
	          $(LESPREMIXED_DIR)LESPremixed2DInput.h \
	          $(LESPREMIXED_DIR)LESPremixed2DQuad.h \
	          $(LESPREMIXED_DIR)LESPremixed2DdRdU.h \
	          $(LESPREMIXED_DIR)LESPremixed2DQuadMultigrid.h \
	          $(LESPREMIXED_DIR)LESPremixed2DQuadNKS.h \
	          $(LESPREMIXED_DIR)LESPremixed2DTurbInit.h \
	          $(LESPREMIXED_DIR)LESPremixed2DReactions.h \
	          $(LESPREMIXED_DIR)PowerLaw.h \
	          $(LESPREMIXED_DIR)Scalars.h 
INC_NAVIERSTOKES = $(NAVIERSTOKES_DIR)NavierStokes2DState.h \
	           $(NAVIERSTOKES_DIR)NavierStokes2DInput.h \
	           $(NAVIERSTOKES_DIR)NavierStokes2DQuad.h \
	           $(NAVIERSTOKES_DIR)NavierStokes2DQuadMultigrid.h
INC_NASARP1311 = $(NASARP1311_DIR)NASARP1311data.h \
		 $(NASARP1311_DIR)LennardJones.h
INC_TURBULENT = $(TURBULENT_DIR)Turbulent2DWallData.h \
	        $(TURBULENT_DIR)TurbulenceModelling.h \
	        $(TURBULENT_DIR)SFSModelling.h
INC_RTE = $(RTE_DIR)Rte2DState.h \
	  $(RTE_DIR)Rte2DInput.h \
	  $(RTE_DIR)Rte2DQuad.h \
	  $(RTE_DIR)Rte2DQuadMultigrid.h \
	  $(RTE_DIR)Rte2DQuadNKS.h \
	  $(RTE_DIR)Rte2DQuadAMR.h \
	  $(RTE_DIR)Rte2DTools.h \
	  $(RTE_DIR)SNBCK.h \
	  $(RTE_DIR)Planck.h 
INC_MULTIGRID = $(MULTIGRID_DIR)FASMultigrid2D.h \
	        $(MULTIGRID_DIR)FASMultigrid2DInput.h
INC_NKS = $(NKS_DIR)NKS2D.h \
	  $(NKS_DIR)NKSInput2D.h \
	  $(NKS_DIR)GMRES2D.h \
	  $(NKS_DIR)Block_Preconditioner2D.h
INC_ICEMCFD = $(ICEM_DIR)ICEMCFD.h
INC_NASA_ROTORS = $(GRID_DIR)NASARotor37.h \
	          $(GRID_DIR)NASARotor67.h
INC_ELLIPTIC = $(CFD_DIR)EllipticOperatorAnalysis2D.h
INC_REACTIONS = $(REACTIONS_DIR)Reactions.h \
		$(REACTIONS_DIR)15step.h \

#
# Compiled object code files
OBJ_MATH = $(MATH_DIR)Vector2D.o \
           $(MATH_DIR)Tensor2D.o \
           $(MATH_DIR)Tensor3D.o \
           $(MATH_DIR)Spline2D.o \
           $(MATH_DIR)Matrix.o \
           $(MATH_DIR)LinearSystems.o \
	   $(MATH_DIR)LinkedList.o \
	   $(MATH_DIR)Polygon.o \
	   $(MATH_DIR)Polyfit/Polyfit.o \
	   $(MATH_DIR)Polyfit/dp1vlu.o \
	   $(MATH_DIR)Polyfit/dpcoef.o \
	   $(MATH_DIR)Polyfit/dpolft.o \
	   $(MATH_DIR)complexify.o 
OBJ_CFD = $(CFD_DIR)CFD.o
OBJ_MPI = $(MPI_DIR)MPI.o
OBJ_SYSTEM = $(SYSTEM_DIR)System_Linux.o
OBJ_GRID = $(GRID_DIR)Cell1D.o \
           $(GRID_DIR)Cell2D.o \
           $(GRID_DIR)Grid2DQuadSingleBlock.o \
           $(GRID_DIR)Grid2DQuadMultiBlock.o \
           $(GRID_DIR)Grid3DHexaSingleBlock.o \
           $(GRID_DIR)Grid3DHexaMultiBlock.o
ifeq ($(MPI_VERSION),NoMPI)
  OBJ_AMR = $(AMR_DIR)AdaptiveBlock.o \
	    $(AMR_DIR)AdaptiveBlock_NoMPI.o \
            $(AMR_DIR)QuadTree.o
else
  OBJ_AMR = $(AMR_DIR)AdaptiveBlock.o \
	    $(AMR_DIR)AdaptiveBlock_MPI.o \
            $(AMR_DIR)QuadTree.o
endif
OBJ_SCALAR = $(SOLVERS1D_DIR)Scalar1D.o \
             $(SOLVERS1D_DIR)Scalar1DSolvers.o
OBJ_HEAT = $(SOLVERS1D_DIR)Heat1D.o \
           $(SOLVERS1D_DIR)Heat1DSolvers.o
OBJ_HYPERHEAT = $(SOLVERS1D_DIR)HyperHeat1DState.o \
                $(SOLVERS1D_DIR)HyperHeat1D.o \
                $(SOLVERS1D_DIR)HyperHeat1DSolvers.o
OBJ_EULER1D = $(SOLVERS1D_DIR)Euler1DState.o \
              $(SOLVERS1D_DIR)Euler1D.o \
              $(SOLVERS1D_DIR)Euler1DSolvers.o
OBJ_MHD1D = $(SOLVERS1D_DIR)MHD1DState.o \
            $(SOLVERS1D_DIR)MHD1D.o \
            $(SOLVERS1D_DIR)MHD1DSolvers.o
OBJ_ADVECTDIFFUSE = $(ADVECTDIFFUSE_DIR)AdvectDiffuse2DState.o \
	            $(ADVECTDIFFUSE_DIR)AdvectDiffuse2DInput.o \
                    $(ADVECTDIFFUSE_DIR)AdvectDiffuse2DQuad.o \
                    $(ADVECTDIFFUSE_DIR)AdvectDiffuse2DQuadSingleBlock.o \
	            $(ADVECTDIFFUSE_DIR)AdvectDiffuse2DQuadMultiBlock.o \
                    $(ADVECTDIFFUSE_DIR)AdvectDiffuse2DQuadGrid.o \
	            $(ADVECTDIFFUSE_DIR)AdvectDiffuse2DQuadSolvers.o \
		    $(ADVECTDIFFUSE_DIR)AdvectDiffuse2DQuad_NKS.o
OBJ_ADVECTDIFFUSE_EXE = $(ADVECTDIFFUSE_DIR)AdvectDiffuse2D.o
OBJ_CHEM = $(CHEM_DIR)Chem2DState.o \
	   $(CHEM_DIR)dRdU.o \
	   $(CHEM_DIR)Chem2DTools.o \
	   $(CHEM_DIR)Chem2DInput.o \
	   $(CHEM_DIR)Chem2DQuad.o \
	   $(CHEM_DIR)Chem2DQuadGrid.o \
	   $(CHEM_DIR)Chem2DQuadMultiBlock.o \
	   $(CHEM_DIR)Chem2DQuadSingleBlock.o \
	   $(CHEM_DIR)Chem2DQuadSolvers.o 
OBJ_CHEM_EXE = $(CHEM_DIR)Chem2D.o
OBJ_HIGHTEMP = $(HIGHTEMP_DIR)HighTemp2DState.o \
               $(HIGHTEMP_DIR)HighTemp2DInput.o  \
               $(HIGHTEMP_DIR)HighTemp2DQuad.o  \
               $(HIGHTEMP_DIR)HighTemp2DQuadGrid.o \
               $(HIGHTEMP_DIR)HighTemp2DQuadMultiBlock.o  \
               $(HIGHTEMP_DIR)HighTemp2DQuadSingleBlock.o  \
               $(HIGHTEMP_DIR)HighTemp2DQuadIOMultiBlock.o  \
               $(HIGHTEMP_DIR)HighTemp2DQuadIOSingleBlock.o  \
               $(HIGHTEMP_DIR)HighTemp2DQuadTurbulenceMultiBlock.o  \
               $(HIGHTEMP_DIR)HighTemp2DQuadTurbulenceSingleBlock.o  \
               $(HIGHTEMP_DIR)HighTemp2DdRdU.o  \
               $(HIGHTEMP_DIR)HighTemp2DQuadSolvers.o
OBJ_HIGHTEMP_EXE = $(HIGHTEMP_DIR)HighTemp2D.o
OBJ_EOS = $(EOS_DIR)tgas.o \
          $(EOS_DIR)ugas.o
OBJ_GAUSSIAN = $(GAUSSIAN_DIR)Gaussian2DState.o $(GAUSSIAN_DIR)Gaussian2DInput.o \
               $(GAUSSIAN_DIR)Gaussian2DCartesian.o $(GAUSSIAN_DIR)Gaussian2DCartesianSolvers.o \
               $(GAUSSIAN_DIR)Gaussian2DQuadGrid.o $(GAUSSIAN_DIR)Gaussian2DQuadSingleBlock.o \
               $(GAUSSIAN_DIR)Gaussian2DQuadMultiBlock.o $(GAUSSIAN_DIR)Gaussian2DQuadSolvers.o \
               $(GAUSSIAN_DIR)Gaussian2DQuadIOSingleBlock.o $(GAUSSIAN_DIR)Gaussian2DQuadIOMultiBlock.o
OBJ_GAUSSIAN_EXE = $(GAUSSIAN_DIR)Gaussian2D.o
OBJ_DUSTY = $(DUSTY_DIR)Particle2DComponents.o \
            $(DUSTY_DIR)Dusty2DState.o \
            $(DUSTY_DIR)Dusty2DInput.o \
            $(DUSTY_DIR)Dusty2DQuad.o \
            $(DUSTY_DIR)Dusty2DQuadSingleBlock.o \
            $(DUSTY_DIR)Dusty2DQuadMultiBlock.o \
            $(DUSTY_DIR)Dusty2DQuadIOSingleBlock.o \
            $(DUSTY_DIR)Dusty2DQuadIOMultiBlock.o \
            $(DUSTY_DIR)Dusty2DQuadTurbulenceSingleBlock.o \
            $(DUSTY_DIR)Dusty2DQuadTurbulenceMultiBlock.o \
            $(DUSTY_DIR)Dusty2DQuadGrid.o \
            $(DUSTY_DIR)Dusty2DQuadSolvers.o
OBJ_DUSTY_EXE = $(DUSTY_DIR)Dusty2D.o
OBJ_ELECTROSTATIC = $(ELECTROSTATIC_DIR)Electrostatic2DState.o \
	            $(ELECTROSTATIC_DIR)Electrostatic2DInput.o \
                    $(ELECTROSTATIC_DIR)Electrostatic2DQuad.o \
                    $(ELECTROSTATIC_DIR)Electrostatic2DQuadSingleBlock.o \
                    $(ELECTROSTATIC_DIR)Electrostatic2DQuadMultiBlock.o \
                    $(ELECTROSTATIC_DIR)Electrostatic2DQuadIOSingleBlock.o \
                    $(ELECTROSTATIC_DIR)Electrostatic2DQuadIOMultiBlock.o \
                    $(ELECTROSTATIC_DIR)Electrostatic2DQuadGrid.o \
                    $(ELECTROSTATIC_DIR)Electrostatic2DQuadSolvers.o
OBJ_ELECTROSTATIC_EXE = $(ELECTROSTATIC_DIR)Electrostatic2D.o
OBJ_EULER = $(EULER_DIR)Euler2DState.o \
            $(EULER_DIR)Euler2DInput.o \
            $(EULER_DIR)Euler2DQuad.o \
            $(EULER_DIR)Euler2DQuadSingleBlock.o \
            $(EULER_DIR)Euler2DQuadMultiBlock.o \
            $(EULER_DIR)Euler2DQuadGrid.o \
            $(EULER_DIR)Euler2DQuadSolvers.o 
OBJ_EULER_EXE = $(EULER_DIR)Euler2D.o
OBJ_INTERFACE = $(INTERFACE_DIR)Interface2D.o
OBJ_EMBEDDEDBOUNDARIES_EXE = $(INTERFACE_DIR)EmbeddedBoundaries2D.o
OBJ_ION = $(ION_DIR)Ion5Moment2DState.o \
          $(ION_DIR)Ion5Moment2DInput.o \
          $(ION_DIR)Ion5Moment2DQuad.o \
          $(ION_DIR)Ion5Moment2DQuadSingleBlock.o \
          $(ION_DIR)Ion5Moment2DQuadMultiBlock.o \
          $(ION_DIR)Ion5Moment2DQuadGrid.o \
          $(ION_DIR)Ion5Moment2DQuadSolvers.o
OBJ_ION_EXE = $(ION_DIR)Ion5Moment2D.o
OBJ_LESPREMIXED = $(LESPREMIXED_DIR)LESPremixed2DState.o \
	          $(LESPREMIXED_DIR)LESPremixed2DInput.o \
	          $(LESPREMIXED_DIR)LESPremixed2DQuad.o \
	          $(LESPREMIXED_DIR)LESPremixed2DQuadSingleBlock.o \
	          $(LESPREMIXED_DIR)LESPremixed2DQuadMultiBlock.o \
	          $(LESPREMIXED_DIR)LESPremixed2DQuadGrid.o \
	          $(LESPREMIXED_DIR)LESPremixed2DQuadSolvers.o \
	          $(LESPREMIXED_DIR)LESPremixed2DTools.o \
	          $(LESPREMIXED_DIR)LESPremixed2DTurbInit.o \
	          $(LESPREMIXED_DIR)LESPremixed2DdRdU.o \
	          $(LESPREMIXED_DIR)PowerLaw.o \
	          $(LESPREMIXED_DIR)Scalars.o 
OBJ_LESPREMIXED_EXE = $(LESPREMIXED_DIR)LESPremixed2D.o
OBJ_LEVELSET = $(LEVELSET_DIR)LevelSet2DInput.o \
               $(LEVELSET_DIR)LevelSet2DQuad.o \
               $(LEVELSET_DIR)LevelSet2DQuadSingleBlock.o \
               $(LEVELSET_DIR)LevelSet2DQuadMultiBlock.o \
               $(LEVELSET_DIR)LevelSet2DQuadIOSingleBlock.o \
               $(LEVELSET_DIR)LevelSet2DQuadIOMultiBlock.o \
               $(LEVELSET_DIR)LevelSet2DQuadInterfaceSingleBlock.o \
               $(LEVELSET_DIR)LevelSet2DQuadInterfaceMultiBlock.o \
               $(LEVELSET_DIR)LevelSet2DQuadHamiltonJacobiSingleBlock.o \
               $(LEVELSET_DIR)LevelSet2DQuadHamiltonJacobiMultiBlock.o \
               $(LEVELSET_DIR)LevelSet2DQuadEikonalSingleBlock.o \
               $(LEVELSET_DIR)LevelSet2DQuadEikonalMultiBlock.o \
               $(LEVELSET_DIR)LevelSet2DQuadScalarExtensionSingleBlock.o \
               $(LEVELSET_DIR)LevelSet2DQuadScalarExtensionMultiBlock.o \
               $(LEVELSET_DIR)LevelSet2DQuadGrid.o \
               $(LEVELSET_DIR)LevelSet2DQuadSolvers.o
OBJ_LEVELSET_EXE = $(LEVELSET_DIR)LevelSet2D.o
OBJ_NAVIERSTOKES = $(NAVIERSTOKES_DIR)NavierStokes2DState.o \
                   $(NAVIERSTOKES_DIR)NavierStokes2DInput.o \
                   $(NAVIERSTOKES_DIR)NavierStokes2DQuad.o \
                   $(NAVIERSTOKES_DIR)NavierStokes2DQuadSingleBlock.o \
                   $(NAVIERSTOKES_DIR)NavierStokes2DQuadMultiBlock.o \
                   $(NAVIERSTOKES_DIR)NavierStokes2DQuadIOSingleBlock.o \
                   $(NAVIERSTOKES_DIR)NavierStokes2DQuadIOMultiBlock.o \
                   $(NAVIERSTOKES_DIR)NavierStokes2DQuadTurbulenceSingleBlock.o \
                   $(NAVIERSTOKES_DIR)NavierStokes2DQuadTurbulenceMultiBlock.o \
                   $(NAVIERSTOKES_DIR)NavierStokes2DQuadGrid.o \
                   $(NAVIERSTOKES_DIR)NavierStokes2DQuadSolvers.o
OBJ_NAVIERSTOKES_EXE = $(NAVIERSTOKES_DIR)NavierStokes2D.o 
OBJ_RTE = $(RTE_DIR)Rte2DState.o \
          $(RTE_DIR)Rte2DInput.o \
          $(RTE_DIR)Rte2DQuad.o \
          $(RTE_DIR)Rte2DQuadSingleBlock.o \
          $(RTE_DIR)Rte2DQuadMultiBlock.o \
          $(RTE_DIR)Rte2DQuadGrid.o \
          $(RTE_DIR)Rte2DQuadSolvers.o \
	  $(RTE_DIR)Rte2DTools.o \
          $(RTE_DIR)SNBCK.o
OBJ_RTE_EXE = $(RTE_DIR)Rte2D.o
OBJ_NASARP1311 = $(NASARP1311_DIR)NASARP1311dataclass.o  
OBJ_TURBULENT = $(TURBULENT_DIR)SFSModelling.o
OBJ_ICEMCFD = $(ICEM_DIR)ICEMCFD.o \
	      $(ICEM_DIR)ICEMCFD_ctype.o
OBJ_REACTIONS = $(REACTIONS_DIR)Reactions.o \
		$(REACTIONS_DIR)15step_GRI211.o \
		$(REACTIONS_DIR)15step_GRI3.o \
		$(REACTIONS_DIR)c_15step_GRI211.o \
		$(REACTIONS_DIR)c_15step_GRI3.o
OBJ_BASE = $(OBJ_MATH) $(OBJ_CFD) $(OBJ_MPI) $(OBJ_GRID) $(OBJ_AMR) \
	   $(OBJ_SYSTEM)
OBJ_CFFC = $(OBJ_BASE) $(OBJ_SCALAR) $(OBJ_HEAT) $(OBJ_HYPERHEAT) $(OBJ_EULER1D) $(OBJ_MHD1D) \
	   $(OBJ_ADVECTDIFFUSE) $(OBJ_CHEM) $(OBJ_DUSTY) $(OBJ_EULER) $(OBJ_HIGHTEMP) \
	   $(OBJ_ION) $(OBJ_ELECTROSTATIC) $(OBJ_INTERFACE) $(OBJ_LESPREMIXED) $(OBJ_LEVELSET) \
	   $(OBJ_NAVIERSTOKES) $(OBJ_GAUSSIAN) $(OBJ_EOS) $(OBJ_ICEMCFD) $(OBJ_TURBULENT) \
	   $(OBJ_NASARP1311) $(OBJ_REACTIONS)
OBJ_MAIN = cffc2D.o
OBJ_ICEMCFD_UTIL1 = $(ICEM_DIR)icemcfd2tecplot.o \
                    $(OBJ_MATH) $(OBJ_CFD) $(OBJ_MPI) $(OBJ_GRID) $(OBJ_ICEMCFD)
OBJ_ICEMCFD_UTIL2 = $(ICEM_DIR)icemcfd2grid2dquad.o \
                    $(OBJ_MATH) $(OBJ_CFD) $(OBJ_MPI) $(OBJ_GRID) $(OBJ_ICEMCFD) 

#
# Define the executables
EXE_MAIN = cffc2D
EXE_ICEMCFD_UTIL1 = icemcfd2tecplot
EXE_ICEMCFD_UTIL2 = icemcfd2grid2dquad
EXE_ADVECTDIFFUSE2D = advectdiffuse2D
EXE_CHEM2D = chem2D
EXE_HIGHTEMP2D = hightemp2D
EXE_DUSTY2D = dusty2D
EXE_ELECTROSTATIC2D = electrostatic2D
EXE_EULER2D = euler2D
EXE_ION5MOMENT2D = ion5moment2D
EXE_LESPREMIXED2D = lespremixed2D
EXE_LEVELSET2D = levelset2D
EXE_NAVIERSTOKES2D = navierstokes2D
EXE_RTE2D = rte2D
EXE_EMBEDDEDBOUNDARIES2D = embeddedboundaries2D
EXE_GAUSSIAN2D = gaussian2D

#
# implicit rule for f90
%.o : %.f90
	$(FC) -c $(FFLAGS) -o $@ $<

#
# Define dependencies
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
	@echo "|        make advectdiffuse2D      compile advectdiffuse2D stand alone program      |"
	@echo "|        make chem2D               compile chem2D stand alone program               |"
	@echo "|        make dusty2D              compile dusty2D stand alone program              |"
	@echo "|        make electrostatic2D      compile electrostatic2D stand alone program      |"
	@echo "|        make euler2D              compile euler2D stand alone program              |"
	@echo "|        make hightemp2D           compile hightemp2D stand alone program           |"
	@echo "|        make ion5moment2D         compile ion5moment2D stand alone program         |"
	@echo "|        make lespremixed2D        compile lespremixed2D stand alone program        |"
	@echo "|        make levelset2D           compile levelset2D stand alone program           |"
	@echo "|        make navierstokes2D       compile navierstokes2D stand alone program       |"
	@echo "|        make embeddedboundaries2D compile embeddedboundaries2D stand alone program |"
	@echo "|        make gaussian2D           compile gaussian2D stand alone program           |"
	@echo "|        make rte2D                compile rte2D stand alone program                |"
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
	@echo "|                  PROFILING = OFF (default), ON                                    |"
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

all: $(EXE_MAIN) $(EXE_ICEMCFD_UTIL1) $(EXE_ICEMCFD_UTIL2) $(EXE_CHEM2D) $(EXE_ADVECTDIFFUSE2D) \
     $(EXE_DUSTY2D) $(EXE_EULER2D) $(EXE_ION5MOMENT2D) $(EXE_LEVELSET2D) $(EXE_NAVIERSTOKES2D) \
     $(EXE_EMBEDDEDBOUNDARIES2D) $(EXE_GAUSSIAN2D) $(EXE_LESPREMIXED2D) $(EXE_HIGHTEMP2D) $(EXE_RTE2D)
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

$(EXE_MAIN): CFFC.depend cffc2D.depend $(OBJ_MAIN) $(OBJ_CFFC)
	@echo ' '
	@echo Building program $(EXE_MAIN) for a $(PLATFORM) platform.
	@echo '--------------------------------------------------------'
	@echo ' '
	$(LD) -o $(EXE_MAIN) $(OBJ_MAIN) $(OBJ_CFFC) $(LDFLAGS)
	@echo ' '
	@echo Program $(EXE_MAIN) has been brought up to date.
	@echo '------------------------------------------------'
	@echo ' '

CFFC.depend: $(SRC_CFFC)
	@echo ' '
	@echo Determining CFFC source file dependencies for a $(PLATFORM) platform.
	@echo '--------------------------------------------------------------------------------'
	@echo ' '
	touch CFFC.depend
	makedepend -fCFFC.depend $(DEPENDSFLAGS) $(CPPFLAGS) $(SRC_CFFC) >& /dev/null
	@echo ' '
	@echo CFFC source file dependencies have been brought up to date.
	@echo '-----------------------------------------------------------------'
	@echo ' '

cffc2D.depend: $(SRC_MAIN)
	@echo ' '
	@echo Determining cffc2D source file dependencies for a $(PLATFORM) platform.
	@echo '--------=--------------------------------------------------------------'
	@echo ' '
	touch cffc2D.depend
	makedepend -fcffc2D.depend $(DEPENDSFLAGS) $(CPPFLAGS) $(SRC_MAIN) >& /dev/null
	@echo ' '
	@echo cffc2D source file dependencies have been brought up to date.
	@echo '-----------------------------------------------------------------'
	@echo ' '

$(EXE_ADVECTDIFFUSE2D): $(ADVECTDIFFUSE_DIR)AdvectDiffuse2D.depend $(OBJ_ADVECTDIFFUSE_EXE) $(OBJ_BASE) $(OBJ_ADVECTDIFFUSE) $(OBJ_ICEMCFD)
	@echo ' '
	@echo Building program $(EXE_ADVECTDIFFUSE2D) for a $(PLATFORM) platform.
	@echo '----------------------------------------------------------'
	@echo ' '
	$(LD) -o $(EXE_ADVECTDIFFUSE2D) $(OBJ_ADVECTDIFFUSE_EXE) $(OBJ_BASE) $(OBJ_ADVECTDIFFUSE) $(OBJ_ICEMCFD) $(LDFLAGS)
	@echo ' '
	@echo Program $(EXE_ADVECTDIFFUSE2D) has been brought up to date.
	@echo '--------------------------------------------------'
	@echo ' '

$(ADVECTDIFFUSE_DIR)AdvectDiffuse2D.depend: $(SRC_ADVECTDIFFUSE_EXE) 
	@echo ' '
	@echo Determining AdvectDiffuse2D source file dependencies for a $(PLATFORM) platform.
	@echo '-----------------------------------------------------------------------'
	@echo ' '
	touch $(ADVECTDIFFUSE_DIR)AdvectDiffuse2D.depend
	makedepend -f$(ADVECTDIFFUSE_DIR)AdvectDiffuse2D.depend $(DEPENDSFLAGS) $(CPPFLAGS) $(SRC_ADVECTDIFFUSE_EXE) $(SRC_ADVECTDIFFUSE) >& /dev/null
	@echo ' '
	@echo AdvectDiffuse2D source file dependencies have been brought up to date.
	@echo '-------------------------------------------------------------'
	@echo ' '

$(EXE_CHEM2D): $(CHEM_DIR)Chem2D.depend $(OBJ_CHEM_EXE) $(OBJ_BASE) $(OBJ_NASARP1311) $(OBJ_REACTIONS) $(OBJ_CHEM) $(OBJ_ICEMCFD)
	@echo ' '
	@echo Building program $(EXE_CHEM2D) for a $(PLATFORM) platform.
	@echo '----------------------------------------------------------'
	@echo ' '
	$(LD) -o $(EXE_CHEM2D) $(OBJ_CHEM_EXE) $(OBJ_BASE) $(OBJ_NASARP1311) $(OBJ_REACTIONS) $(OBJ_CHEM) $(OBJ_ICEMCFD) $(LDFLAGS)
	@echo ' '
	@echo Program $(EXE_CHEM2D) has been brought up to date.
	@echo '--------------------------------------------------'
	@echo ' '

$(CHEM_DIR)Chem2D.depend: $(SRC_CHEM_EXE) 
	@echo ' '
	@echo Determining Chem2D source file dependencies for a $(PLATFORM) platform.
	@echo '-----------------------------------------------------------------------'
	@echo ' '
	touch $(CHEM_DIR)Chem2D.depend
	makedepend -f$(CHEM_DIR)Chem2D.depend $(DEPENDSFLAGS) $(CPPFLAGS) $(SRC_CHEM_EXE) $(SRC_CHEM) >& /dev/null
	@echo ' '
	@echo Chem2D source file dependencies have been brought up to date.
	@echo '-------------------------------------------------------------'
	@echo ' '

$(EXE_HIGHTEMP2D): $(HIGHTEMP_DIR)HighTemp2D.depend $(OBJ_HIGHTEMP_EXE) $(OBJ_BASE) $(OBJ_INTERFACE) $(OBJ_HIGHTEMP) $(OBJ_ICEMCFD) $(OBJ_EOS)
	@echo ' '
	@echo Building program $(EXE_HIGHTEMP2D) for a $(PLATFORM) platform.
	@echo '----------------------------------------------------------'
	@echo ' '
	$(LD) -o $(EXE_HIGHTEMP2D) $(OBJ_HIGHTEMP_EXE) $(OBJ_BASE) $(OBJ_INTERFACE) $(OBJ_HIGHTEMP) $(OBJ_ICEMCFD) $(LDFLAGS) $(OBJ_EOS)
	@echo ' '
	@echo Program $(EXE_HIGHTEMP2D) has been brought up to date.
	@echo '--------------------------------------------------'
	@echo ' '

$(HIGHTEMP_DIR)HighTemp2D.depend: $(SRC_HIGHTEMP_EXE) 
	@echo ' '
	@echo Determining HighTemp2D source file dependencies for a $(PLATFORM) platform.
	@echo '-----------------------------------------------------------------------'
	@echo ' '
	touch $(HIGHTEMP_DIR)HighTemp2D.depend
	makedepend -f$(HIGHTEMP_DIR)HighTemp2D.depend $(DEPENDSFLAGS) $(CPPFLAGS) $(SRC_HIGHTEMP_EXE) >& /dev/null
	@echo ' '
	@echo HighTemp2D source file dependencies have been brought up to date.
	@echo '-------------------------------------------------------------'
	@echo ' '

$(EXE_DUSTY2D): $(DUSTY_DIR)Dusty2D.depend $(OBJ_DUSTY_EXE) $(OBJ_BASE) $(OBJ_DUSTY) $(OBJ_INTERFACE) $(OBJ_ELECTROSTATIC) $(OBJ_ICEMCFD)
	@echo ' '
	@echo Building program $(EXE_DUSTY2D) for a $(PLATFORM) platform.
	@echo '-----------------------------------------------------------'
	@echo ' '
	$(LD) -o $(EXE_DUSTY2D) $(OBJ_DUSTY_EXE) $(OBJ_BASE) $(OBJ_DUSTY) $(OBJ_INTERFACE) $(OBJ_ELECTROSTATIC) $(OBJ_ICEMCFD) $(LDFLAGS)
	@echo ' '
	@echo Program $(EXE_DUSTY2D) has been brought up to date.
	@echo '---------------------------------------------------'
	@echo ' '

$(DUSTY_DIR)Dusty2D.depend: $(SRC_DUSTY_EXE)
	@echo ' '
	@echo Determining $(DUSTY_DIR)Dusty2D source file dependencies for a $(PLATFORM) platform.
	@echo '------------------------------------------------------------------------'
	@echo ' '
	touch $(DUSTY_DIR)Dusty2D.depend
	makedepend -f$(DUSTY_DIR)Dusty2D.depend $(DEPENDSFLAGS) $(CPPFLAGS) $(SRC_DUSTY_EXE) $(SRC_DUSTY) >& /dev/null
	@echo ' '
	@echo Dusty2D source file dependencies have been brought up to date.
	@echo '--------------------------------------------------------------'
	@echo ' '

$(EXE_EULER2D): $(EULER_DIR)Euler2D.depend $(OBJ_EULER_EXE) $(OBJ_BASE) $(OBJ_EULER) $(OBJ_INTERFACE) $(OBJ_ICEMCFD)
	@echo ' '
	@echo Building program $(EXE_EULER2D) for a $(PLATFORM) platform.
	@echo '-----------------------------------------------------------'
	@echo ' '
	$(LD) -o $(EXE_EULER2D) $(OBJ_EULER_EXE) $(OBJ_BASE) $(OBJ_EULER) $(OBJ_INTERFACE) $(OBJ_ICEMCFD) $(LDFLAGS)
	@echo ' '
	@echo Program $(EXE_EULER2D) has been brought up to date.
	@echo '---------------------------------------------------'
	@echo ' '

$(EULER_DIR)Euler2D.depend: $(SRC_EULER_EXE)
	@echo ' '
	@echo Determining Euler2D source file dependencies for a $(PLATFORM) platform.
	@echo '------------------------------------------------------------------------'
	@echo ' '
	touch $(EULER_DIR)Euler2D.depend
	makedepend -f$(EULER_DIR)Euler2D.depend $(DEPENDSFLAGS) $(CPPFLAGS) $(SRC_EULER_EXE) $(SRC_EULER) >& /dev/null
	@echo ' '
	@echo Euler2D source file dependencies have been brought up to date.
	@echo '--------------------------------------------------------------'
	@echo ' '

$(EXE_ELECTROSTATIC2D): $(ELECTROSTATIC_DIR)Electrostatic2D.depend $(OBJ_ELECTROSTATIC_EXE) $(OBJ_BASE) $(OBJ_ELECTROSTATIC) $(OBJ_ICEMCFD)
	@echo ' '
	@echo Building program $(EXE_ELECTROSTATIC2D) for a $(PLATFORM) platform.
	@echo '-----------------------------------------------------------'
	@echo ' '
	$(LD) -o $(EXE_ELECTROSTATIC2D) $(OBJ_ELECTROSTATIC_EXE) $(OBJ_BASE) $(OBJ_ELECTROSTATIC) $(OBJ_ICEMCFD) $(LDFLAGS)
	@echo ' '
	@echo Program $(EXE_ELECTROSTATIC2D) has been brought up to date.
	@echo '---------------------------------------------------'
	@echo ' '

$(ELECTROSTATIC_DIR)Electrostatic2D.depend: $(SRC_ELECTROSTATIC_EXE)
	@echo ' '
	@echo Determining Electrostatic2D source file dependencies for a $(PLATFORM) platform.
	@echo '------------------------------------------------------------------------'
	@echo ' '
	touch $(ELECTROSTATIC_DIR)Electrostatic2D.depend
	makedepend -f$(ELECTROSTATIC_DIR)Electrostatic2D.depend $(DEPENDSFLAGS) $(CPPFLAGS) $(SRC_ELECTROSTATIC_EXE) $(SRC_ELECTROSTATIC) >& /dev/null
	@echo ' '
	@echo Electrostatic2D source file dependencies have been brought up to date.
	@echo '--------------------------------------------------------------'
	@echo ' '

$(EXE_ION5MOMENT2D): $(ION_DIR)Ion5Moment2D.depend $(OBJ_ION_EXE) $(OBJ_BASE) $(OBJ_ION) $(OBJ_EULER) $(OBJ_INTERFACE) $(OBJ_ICEMCFD)
	@echo ' '
	@echo Building program $(EXE_ION5MOMENT2D) for a $(PLATFORM) platform.
	@echo '----------------------------------------------------------------'
	@echo ' '
	$(LD) -o $(EXE_ION5MOMENT2D) $(OBJ_ION_EXE) $(OBJ_BASE) $(OBJ_ION) $(OBJ_EULER) $(OBJ_INTERFACE) $(OBJ_ICEMCFD) $(LDFLAGS)
	@echo ' '
	@echo Program $(EXE_ION5MOMENT2D) has been brought up to date.
	@echo '--------------------------------------------------------'
	@echo ' '

$(ION_DIR)Ion5Moment2D.depend: $(SRC_ION_EXE) 
	@echo ' '
	@echo Determining Ion5Moment2D source file dependencies for a $(PLATFORM) platform.
	@echo '-----------------------------------------------------------------------------'
	@echo ' '
	touch $(ION_DIR)Ion5Moment2D.depend
	makedepend -f$(ION_DIR)Ion5Moment2D.depend $(DEPENDSFLAGS) $(CPPFLAGS) $(SRC_ION_EXE) $(SRC_ION) >& /dev/null
	@echo ' '
	@echo Ion5Moment2D source file dependencies have been brought up to date.
	@echo '-------------------------------------------------------------------'
	@echo ' '

$(EXE_LESPREMIXED2D): $(LESPREMIXED_DIR)LESPremixed2D.depend $(OBJ_LESPREMIXED_EXE) $(OBJ_BASE) $(OBJ_NASARP1311) $(OBJ_REACTIONS) $(OBJ_LESPREMIXED) $(OBJ_TURBULENT) $(OBJ_ICEMCFD) 
	@echo ' '
	@echo Building program $(EXE_LESPREMIXED2D) for a $(PLATFORM) platform.
	@echo '----------------------------------------------------------'
	@echo ' '
	$(LD) -o $(EXE_LESPREMIXED2D) $(OBJ_LESPREMIXED_EXE) $(OBJ_BASE) $(OBJ_NASARP1311) $(OBJ_REACTIONS) $(OBJ_LESPREMIXED) $(OBJ_TURBULENT) $(OBJ_ICEMCFD) $(LDFLAGS)
	@echo ' '
	@echo Program $(EXE_LESPREMIXED2D) has been brought up to date.
	@echo '--------------------------------------------------'
	@echo ' '

$(LESPREMIXED_DIR)LESPremixed2D.depend: $(SRC_LESPREMIXED_EXE) 
	@echo ' '
	@echo Determining LESPremixed2D source file dependencies for a $(PLATFORM) platform.
	@echo '-----------------------------------------------------------------------'
	@echo ' '
	touch $(LESPREMIXED_DIR)LESPremixed2D.depend
	makedepend -f$(LESPREMIXED_DIR)LESPremixed2D.depend $(DEPENDSFLAGS) $(CPPFLAGS) $(SRC_LESPREMIXED_EXE) >& /dev/null
	@echo ' '
	@echo LESPremixed2D source file dependencies have been brought up to date.
	@echo '-------------------------------------------------------------'
	@echo ' '

$(EXE_LEVELSET2D): $(LEVELSET_DIR)LevelSet2D.depend $(OBJ_LEVELSET_EXE) $(OBJ_BASE) $(OBJ_LEVELSET) $(OBJ_INTERFACE)
	@echo ' '
	@echo Building program $(EXE_LEVELSET2D) for a $(PLATFORM) platform.
	@echo '------------------------------------------------------------------------'
	@echo ' '
	$(LD) -o $(EXE_LEVELSET2D) $(OBJ_LEVELSET_EXE) $(OBJ_BASE) $(OBJ_LEVELSET) $(OBJ_INTERFACE) $(LDFLAGS)
	@echo ' '
	@echo Program $(EXE_LEVELSET2D) has been brought up to date.
	@echo '------------------------------------------------------------------------'
	@echo ' '

$(LEVELSET_DIR)LevelSet2D.depend: $(SRC_LEVELSET_EXE) 
	@echo ' '
	@echo Determining LevelSet2D source file dependencies for a $(PLATFORM) platform.
	@echo '------------------------------------------------------------------------'
	@echo ' '
	touch $(LEVELSET_DIR)LevelSet2D.depend
	makedepend -f$(LEVELSET_DIR)LevelSet2D.depend $(DEPENDSFLAGS) $(CPPFLAGS) $(SRC_LEVELSET_EXE) $(SRC_LEVELSET) >& /dev/null
	@echo ' '
	@echo LevelSet2D source file dependencies have been brought up to date.
	@echo '------------------------------------------------------------------------'
	@echo ' '

$(EXE_NAVIERSTOKES2D): $(NAVIERSTOKES_DIR)NavierStokes2D.depend $(OBJ_NAVIERSTOKES_EXE) $(OBJ_BASE) $(OBJ_NAVIERSTOKES) $(OBJ_INTERFACE) $(OBJ_ICEMCFD)
	@echo ' '
	@echo Building program $(EXE_NAVIERSTOKES2D) for a $(PLATFORM) platform.
	@echo '------------------------------------------------------------------'
	@echo ' '
	$(LD) -o $(EXE_NAVIERSTOKES2D) $(OBJ_NAVIERSTOKES_EXE) $(OBJ_BASE) $(OBJ_NAVIERSTOKES) $(OBJ_INTERFACE) $(OBJ_ICEMCFD) $(LDFLAGS)
	@echo ' '
	@echo Program $(EXE_NAVIERSTOKES2D) has been brought up to date.
	@echo '----------------------------------------------------------'
	@echo ' '

$(NAVIERSTOKES_DIR)NavierStokes2D.depend: $(SRC_NAVIERSTOKES_EXE)
	@echo ' '
	@echo Determining NavierStokes2D source file dependencies for a $(PLATFORM) platform.
	@echo '-------------------------------------------------------------------------------'
	@echo ' '
	touch $(NAVIERSTOKES_DIR)NavierStokes2D.depend
	makedepend -f$(NAVIERSTOKES_DIR)NavierStokes2D.depend $(DEPENDSFLAGS) $(CPPFLAGS) $(SRC_NAVIERSTOKES_EXE) $(SRC_NAVIERSTOKES) >& /dev/null
	@echo ' '
	@echo NavierStokes2D source file dependencies have been brought up to date.
	@echo '---------------------------------------------------------------------'
	@echo ' '

$(EXE_EMBEDDEDBOUNDARIES2D): $(INTERFACE_DIR)EmbeddedBoundaries2D.depend $(OBJ_EMBEDDEDBOUNDARIES_EXE) $(OBJ_BASE) $(OBJ_EMBEDDEDBOUNDARIES) $(OBJ_INTERFACE) $(OBJ_LEVELSET) $(OBJ_EULER) $(OBJ_NAVIERSTOKES) $(OBJ_DUSTY) $(OBJ_ELECTROSTATIC) $(OBJ_ICEMCFD) $(OBJ_GAUSSIAN)
	@echo ' '
	@echo Building program $(EXE_EMBEDDEDBOUNDARIES2D) for a $(PLATFORM) platform.
	@echo '------------------------------------------------------------------'
	@echo ' '
	$(LD) -o $(EXE_EMBEDDEDBOUNDARIES2D) $(OBJ_EMBEDDEDBOUNDARIES_EXE) $(OBJ_BASE) $(OBJ_EMBEDDEDBOUNDARIES) $(OBJ_INTERFACE) $(OBJ_LEVELSET) $(OBJ_EULER) $(OBJ_NAVIERSTOKES) $(OBJ_DUSTY) $(OBJ_ELECTROSTATIC) $(OBJ_ICEMCFD) $(OBJ_GAUSSIAN) $(LDFLAGS)
	@echo ' '
	@echo Program $(EXE_EMBEDDEDBOUNDARIES2D) has been brought up to date.
	@echo '----------------------------------------------------------'
	@echo ' '

$(INTERFACE_DIR)EmbeddedBoundaries2D.depend: $(SRC_EMBEDDEDBOUNDARIES_EXE)
	@echo ' '
	@echo Determining EmbeddedBoundaries2D source file dependencies for a $(PLATFORM) platform.
	@echo '-------------------------------------------------------------------------------'
	@echo ' '
	touch $(INTERFACE_DIR)EmbeddedBoundaries2D.depend
	makedepend -f$(INTERFACE_DIR)EmbeddedBoundaries2D.depend $(DEPENDSFLAGS) $(CPPFLAGS) $(SRC_EMBEDDEDBOUNDARIES_EXE) $(SRC_EMBEDDEDBOUNDARIES) >& /dev/null
	@echo ' '
	@echo EmbeddedBoundaries2D source file dependencies have been brought up to date.
	@echo '---------------------------------------------------------------------'
	@echo ' '

$(EXE_GAUSSIAN2D): $(GAUSSIAN_DIR)Gaussian2D.depend $(OBJ_GAUSSIAN_EXE) $(OBJ_BASE) $(OBJ_GAUSSIAN) $(OBJ_INTERFACE) $(OBJ_ICEMCFD)
	@echo ' '
	@echo Building program $(EXE_GAUSSIAN2D) for a $(PLATFORM) platform.
	@echo '--------------------------------------------------------'
	@echo ' '
	$(LD) -o $(EXE_GAUSSIAN2D) $(OBJ_GAUSSIAN_EXE) $(OBJ_BASE) $(OBJ_GAUSSIAN) $(OBJ_INTERFACE) $(OBJ_ICEMCFD) $(LDFLAGS)
	@echo ' '
	@echo Program $(EXE_GAUSSIAN2D) has been brought up to date.
	@echo '------------------------------------------------'
	@echo ' '

$(GAUSSIAN_DIR)Gaussian2D.depend: $(SRC_GAUSSIAN_EXE) 
	@echo ' '
	@echo Determining Gaussian2D source file dependencies for a $(PLATFORM) platform.
	@echo '-----------------------------------------------------------------------'
	@echo ' '
	touch $(GAUSSIAN_DIR)Gaussian2D.depend
	makedepend -f$(GAUSSIAN_DIR)Gaussian2D.depend $(DEPENDSFLAGS) $(CPPFLAGS) $(SRC_GAUSSIAN_EXE) $(SRC_GAUSSIAN) >& /dev/null
	@echo ' '
	@echo Gaussian2D source file dependencies have been brought up to date.
	@echo '-------------------------------------------------------------'
	@echo ' '

$(EXE_RTE2D): $(RTE_DIR)Rte2D.depend $(OBJ_RTE_EXE) $(OBJ_BASE) $(OBJ_RTE) $(OBJ_ICEMCFD)
	@echo ' '
	@echo Building program $(EXE_RTE2D) for a $(PLATFORM) platform.
	@echo '-----------------------------------------------------------'
	@echo ' '
	$(LD) -o $(EXE_RTE2D) $(OBJ_RTE_EXE) $(OBJ_BASE) $(OBJ_RTE) $(OBJ_ICEMCFD) $(LDFLAGS)
	@echo ' '
	@echo Program $(EXE_RTE2D) has been brought up to date.
	@echo '---------------------------------------------------'
	@echo ' '

$(RTE_DIR)Rte2D.depend: $(SRC_RTE_EXE)
	@echo ' '
	@echo Determining Euler2D source file dependencies for a $(PLATFORM) platform.
	@echo '------------------------------------------------------------------------'
	@echo ' '
	touch $(RTE_DIR)Rte2D.depend
	makedepend -f$(RTE_DIR)Rte2D.depend $(DEPENDSFLAGS) $(CPPFLAGS) $(SRC_RTE_EXE) $(SRC_RTE) >& /dev/null
	@echo ' '
	@echo Rte2D source file dependencies have been brought up to date.
	@echo '--------------------------------------------------------------'
	@echo ' '


$(EXE_ICEMCFD_UTIL1): CFFC.depend $(ICEM_DIR)icemcfd2tecplot.depend $(OBJ_ICEMCFD_UTIL1)
	@echo ' '
	@echo Building program $(EXE_ICEMCFD_UTIL1) for a $(PLATFORM) platform.
	@echo '-----------------------------------------------------------------'
	@echo ' '
	$(LD) -o $(EXE_ICEMCFD_UTIL1) $(OBJ_ICEMCFD_UTIL1) $(LDFLAGS)
	@echo ' '
	@echo Program $(EXE_ICEMCFD_UTIL1) has been brought up to date.
	@echo '---------------------------------------------------------'
	@echo ' '

$(EXE_ICEMCFD_UTIL2): CFFC.depend $(ICEM_DIR)icemcfd2grid2dquad.depend $(OBJ_ICEMCFD_UTIL2)
	@echo ' '
	@echo Building program $(EXE_ICEMCFD_UTIL2) for a $(PLATFORM) platform.
	@echo '-----------------------------------------------------------------'
	@echo ' '
	$(LD) -o $(EXE_ICEMCFD_UTIL2) $(OBJ_ICEMCFD_UTIL2) $(LDFLAGS)
	@echo ' '
	@echo Program $(EXE_ICEMCFD_UTIL2) has been brought up to date.
	@echo '---------------------------------------------------------'
	@echo ' '

$(ICEM_DIR)icemcfd2tecplot.depend: $(SRC_ICEMCFD_UTIL1)
	@echo ' '
	@echo Determining ICEM CFD utility program source file dependencies for a $(PLATFORM) platform.
	@echo '-----------------------------------------------------------------------------------------'
	@echo ' '
	touch $(ICEM_DIR)icemcfd2tecplot.depend
	makedepend -f$(ICEM_DIR)icemcfd2tecplot.depend $(DEPENDSFLAGS) $(CPPFLAGS) $(SRC_ICEMCFD_UTIL1) >& /dev/null
	@echo ' '
	@echo ICEM CFD utility program source file dependencies have been brought up to date.
	@echo '-------------------------------------------------------------------------------'
	@echo ' '

$(ICEM_DIR)icemcfd2grid2dquad.depend: $(SRC_ICEMCFD_UTIL2)
	@echo ' '
	@echo Determining ICEM CFD utility program source file dependencies for a $(PLATFORM) platform.
	@echo '-----------------------------------------------------------------------------------------'
	@echo ' '
	touch $(ICEM_DIR)icemcfd2grid2dquad.depend
	makedepend -f$(ICEM_DIR)icemcfd2grid2dquad.depend $(DEPENDSFLAGS) $(CPPFLAGS) $(SRC_ICEMCFD_UTIL2) >& /dev/null
	@echo ' '
	@echo ICEM CFD utility program source file dependencies have been brought up to date.
	@echo '-------------------------------------------------------------------------------'
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

install-all: $(EXE_MAIN) $(EXE_ICEMCFD_UTIL1) $(EXE_ICEMCFD_UTIL2) \
        $(EXE_ADVECTDIFFUSE2D) $(EXE_CHEM2D) $(EXE_DUSTY2D) $(EXE_EULER2D) \
        $(EXE_ION5MOMENT2D) $(EXE_LESPREMIXED2D) $(EXE_LEVELSET2D) $(EXE_NAVIERSTOKES2D) \
        $(EXE_EMBEDDEDBOUNDARIES2D) $(EXE_GAUSSIAN2D) $(EXE_HIGHTEMP2D)
	@for f in $(EXE_MAIN) $(EXE_ICEMCFD_UTIL1) $(EXE_ICEMCFD_UTIL2) \
	$(EXE_ADVECTDIFFUSE2D) $(EXE_CHEM2D) $(EXE_DUSTY2D) $(EXE_EULER2D) \
        $(EXE_ION5MOMENT2D) $(EXE_LEVELSET2D) $(EXE_NAVIERSTOKES2D) \
        $(EXE_EMBEDDEDBOUNDARIES2D) $(EXE_GAUSSIAN2D) $(EXE_LESPREMIXED2D) \
        $(EXE_HIGHTEMP2D) $(EXE_RTE2D); \
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
	@for f in $(EXE_MAIN) $(EXE_ICEMCFD_UTIL1) $(EXE_ICEMCFD_UTIL2) \
	$(EXE_ADVECTDIFFUSE2D) $(EXE_CHEM2D) $(EXE_DUSTY2D) $(EXE_EULER2D) \
        $(EXE_ION5MOMENT2D) $(EXE_LEVELSET2D) $(EXE_NAVIERSTOKES2D) \
        $(EXE_EMBEDDEDBOUNDARIES2D) $(EXE_LESPREMIXED2D) $(EXE_HIGHTEMP2D) \
        $(EXE_RTE2D); \
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
	rm -f $(SOLVERS1D_DIR)*.o $(SOLVERS1D_DIR)*.il $(SOLVERS1D_DIR)*~ 
	rm -f $(AMR_DIR)*.o $(AMR_DIR)*.il $(AMR_DIR)*~
	rm -f $(ADVECTDIFFUSE_DIR)*.o $(ADVECTDIFFUSE_DIR)*.il $(ADVECTDIFFUSE_DIR)*~ $(ADVECTDIFFUSE_DIR)*.depend $(ADVECTDIFFUSE_DIR)*.bak
	rm -f $(AMR_DIR)*.o $(AMR_DIR)*.il $(AMR_DIR)*~
	rm -f $(CFD_DIR)*.o $(CFD_DIR)*.il $(CFD_DIR)*~
	rm -f $(CHEM_DIR)*.o $(CHEM_DIR)*.il $(CHEM_DIR)*~ $(CHEM_DIR)*.depend $(CHEM_DIR)*.bak
	rm -f $(HIGHTEMP_DIR)*.o $(HIGHTEMP_DIR)*.il $(HIGHTEMP_DIR)*~ $(HIGHTEMP_DIR)*.depend $(HIGHTEMP_DIR)*.bak
	rm -f $(EOS_DIR)*.o $(EOS_DIR)*~
	rm -f $(DUSTY_DIR)*.o $(DUSTY_DIR)*.il $(DUSTY_DIR)*~ $(DUSTY_DIR)*.depend $(DUSTY_DIR)*.bak
	rm -f $(ELECTROSTATIC_DIR)*.o $(ELECTROSTATIC_DIR)*.il $(ELECTROSTATIC_DIR)*~ $(ELECTROSTATIC_DIR)*.depend $(ELECTROSTATIC_DIR)*.bak
	rm -f $(EULER_DIR)*.o $(EULER_DIR)*.il $(EULER_DIR)*~ $(EULER_DIR)*.depend $(EULER_DIR)*.bak
	rm -f $(RTE_DIR)*.o $(RTE_DIR)*.il $(RTE_DIR)*~ $(RTE_DIR)*.depend $(RTE_DIR)*.bak
	rm -f $(MULTIGRID_DIR)*~
	rm -f $(NKS_DIR)*~
	rm -f $(GRID_DIR)*.o $(GRID_DIR)*.il $(GRID_DIR)*~
	rm -f $(ICEM_DIR)*.o $(ICEM_DIR)*.il $(ICEM_DIR)*~ $(ICEM_DIR)*.depend $(ICEM_DIR)*.bak
	rm -f $(ION_DIR)*.o $(ION_DIR)*.il $(ION_DIR)*~ $(ION_DIR)*.depend $(ION_DIR)*.bak
	rm -f $(INTERFACE_DIR)*.o $(INTERFACE_DIR)*.il $(INTERFACE_DIR)*~ $(INTERFACE_DIR)*.depend $(INTERFACE_DIR)*.bak
	rm -f $(LEVELSET_DIR)*.o $(LEVELSET_DIR)*.il $(LEVELSET_DIR)*~ $(LEVELSET_DIR)*.depend $(LEVELSET_DIR)*.bak
	rm -f $(LESPREMIXED_DIR)*.o $(LESPREMIXED_DIR)*.il $(LESPREMIXED_DIR)*~ $(LESPREMIXED_DIR)*.depend $(LESPREMIXED_DIR)*.bak *.mod
	rm -f $(MATH_DIR)*.o $(MATH_DIR)*.il $(MATH_DIR)*~
	rm -f $(MPI_DIR)*.o $(MPI_DIR)*.il $(MPI_DIR)*~
	rm -f $(NAVIERSTOKES_DIR)*.o $(NAVIERSTOKES_DIR)*.il $(NAVIERSTOKES_DIR)*~ $(NAVIERSTOKES_DIR)*.depend $(NAVIERSTOKES_DIR)*.bak
	rm -f $(PHYSICS_DIR)*.o $(PHYSICS_DIR)*.il $(PHYSICS_DIR)*~
	rm -f $(SYSTEM_DIR)*.o $(SYSTEM_DIR)*.il $(SYSTEM_DIR)*~
	rm -f $(TURBULENT_DIR)*.o $(TURBULENT_DIR)*.il $(TURBULENT_DIR)*~
	rm -f $(GAUSSIAN_DIR)*.o $(GAUSSIAN_DIR)*.il $(GAUSSIAN_DIR)*~
	rm -f $(REACTIONS_DIR)*.o $(REACTIONS_DIR)*.il $(REACTIONS_DIR)*~
	rm -f $(NASARP1311_DIR)*.o $(NASARP1311_DIR)*.il $(NASARP1311_DIR)*~
	rm -f $(EXE_MAIN)
	rm -f $(EXE_ICEMCFD_UTIL1) $(EXE_ICEMCFD_UTIL2)
	rm -f $(EXE_ADVECTDIFFUSE2D)
	rm -f $(EXE_CHEM2D)
	rm -f $(EXE_HIGHTEMP2D)
	rm -f $(EXE_DUSTY2D)
	rm -f $(EXE_ELECTROSTATIC2D)
	rm -f $(EXE_EMBEDDEDBOUNDARIES2D)
	rm -f $(EXE_EULER2D)
	rm -f $(EXE_ION5MOMENT2D)
	rm -f $(EXE_LESPREMIXED2D)
	rm -f $(EXE_LEVELSET2D)
	rm -f $(EXE_RTE2D)
	rm -f $(EXE_NAVIERSTOKES2D)
	rm -f $(EXE_GAUSSIAN2D)
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
# Include header file dependencies
include CFFC.depend
include cffc2D.depend
include $(ADVECTDIFFUSE_DIR)AdvectDiffuse2D.depend
include $(CHEM_DIR)Chem2D.depend
include $(HIGHTEMP_DIR)HighTemp2D.depend
include $(DUSTY_DIR)Dusty2D.depend
include $(ELECTROSTATIC_DIR)Electrostatic2D.depend
include $(INTERFACE_DIR)EmbeddedBoundaries2D.depend
include $(EULER_DIR)Euler2D.depend
include $(ION_DIR)Ion5Moment2D.depend
include $(LESPREMIXED_DIR)LESPremixed2D.depend
include $(LEVELSET_DIR)LevelSet2D.depend
include $(NAVIERSTOKES_DIR)NavierStokes2D.depend
include $(RTE_DIR)Rte2D.depend
include $(ICEM_DIR)icemcfd2tecplot.depend
include $(ICEM_DIR)icemcfd2grid2dquad.depend
include $(GAUSSIAN_DIR)Gaussian2D.depend


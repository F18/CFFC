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
AMR_DIR = AMR/
CFD_DIR = CFD/
EULER_DIR = Euler/
FANS_DIR = FANS/
GRID_DIR = Grid/
HEXA_DIR = HexaBlock/
ICEM_DIR = ICEM/
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
	   $(MATH_DIR)Polyfit/dpolft.f
SRC_MATH_F90 = $(MATH_DIR)complexify.f90
SRC_MATH += $(SRC_MATH_F90)

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
SRC_ICEM = $(ICEM_DIR)ICEMCFD.cc \
	   $(ICEM_DIR)ICEMCFD_ctype.c

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
SRC_EULER = $(EULER_DIR)Euler3DPolytropicState.cc \
	    $(EULER_DIR)Euler3DThermallyPerfectState.cc \

# NavierStokes3D:
SRC_NAVIERSTOKES = $(NAVIERSTOKES_DIR)NavierStokes3DThermallyPerfectState.cc \
                   $(NAVIERSTOKES_DIR)NavierStokes3DThermallyPerfectHexaBlock.cc

# FANS3D:
SRC_FANS = $(FANS_DIR)FANS3DThermallyPerfectState.cc \
           $(FANS_DIR)FANS3DThermallyPerfectHexaBlock.cc

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
          BasicTests/test_TestData.cc \
	  ../Euler/test_Euler3DPolytropicState.cc
SRC_TUT := $(addprefix $(UNIT_TESTING_DIR), $(SRC_TUT))

# System:
SRC_SYSTEM = $(SYSTEM_DIR)System_Linux.cc

# Base and CFFC sources:
SRC_BASE = $(SRC_MATH) $(SRC_CFD) $(SRC_MPI) $(SRC_SYSTEM) \
           $(SRC_AMR) $(SRC_GRID) $(SRC_ICEM) \
           $(SRC_MULTIGRID) $(SRC_NKS)
ifeq ($(TUT_TESTING),ON) 
  SRC_BASE += $(SRC_TUT)
endif
SRC_CFFC = $(SRC_BASE) \
           $(SRC_EULER) $(SRC_NAVIERSTOKES) $(SRC_FANS) \
           $(SRC_PHYSICS) $(SRC_TURBULENCE) $(SRC_EOS) $(SRC_NASADATA) $(SRC_REACTIONS)

# Executable sources:
SRC_MAIN = cffc3D.cc

###############################################################################
# Compiled object code files
###############################################################################
# Math:
OBJ_MATH = $(MATH_DIR)Vector2D.o \
           $(MATH_DIR)Tensor2D.o \
           $(MATH_DIR)Third_order_tensor2D.o \
           $(MATH_DIR)Tensor3D.o \
           $(MATH_DIR)Spline2D.o \
           $(MATH_DIR)Matrix.o \
           $(MATH_DIR)LinearSystems.o \
	   $(MATH_DIR)LinkedList.o \
	   $(MATH_DIR)Polygon.o \
	   $(MATH_DIR)Polyfit/Polyfit.o \
	   $(MATH_DIR)Polyfit/dp1vlu.o \
	   $(MATH_DIR)Polyfit/dpcoef.o \
	   $(MATH_DIR)Polyfit/dpolft.o
OBJ_MATH_F90 = $(MATH_DIR)complexify.o 
OBJ_MATH += $(OBJ_MATH_F90)
MOD_MATH_F90 = $(MODULE_DIR)/complexify.mod

# CFD and MPI:
OBJ_CFD = $(CFD_DIR)CFD.o \
          $(CFD_DIR)CFDInput.o
OBJ_MPI	= $(MPI_DIR)MPI.o

# Grid:
OBJ_GRID = $(GRID_DIR)Cell2D.o \
           $(GRID_DIR)Grid2DQuadBlock.o \
           $(GRID_DIR)Grid2DQuadMultiBlock.o \
           $(GRID_DIR)Grid3DInput.o \
           $(GRID_DIR)Grid3DHexaBlock.o \
           $(GRID_DIR)Grid3DHexaMultiBlock.o
OBJ_ICEM = $(ICEM_DIR)ICEMCFD.o \
	   $(ICEM_DIR)ICEMCFD_ctype.o

# AMR:
ifeq ($(MPI_VERSION),NoMPI)
  OBJ_AMR = $(AMR_DIR)AMRInput.o \
            $(AMR_DIR)AdaptiveBlock3D.o \
	    $(AMR_DIR)AdaptiveBlock3D_NoMPI.o \
            $(AMR_DIR)Octree.o
else
  OBJ_AMR = $(AMR_DIR)AMRInput.o \
            $(AMR_DIR)AdaptiveBlock3D.o \
	    $(AMR_DIR)AdaptiveBlock3D_MPI.o \
            $(AMR_DIR)Octree.o
endif

# Multigrid:
OBJ_MULTIGRID = $(MULTIGRID_DIR)FASMultigridInput.o

# NKS:
OBJ_NKS = $(NKS_DIR)NKSInput.o

# Euler3D:
OBJ_EULER = $(EULER_DIR)Euler3DPolytropicState.o \
	    $(EULER_DIR)Euler3DThermallyPerfectState.o \


# NavierStokes3D:
OBJ_NAVIERSTOKES = $(NAVIERSTOKES_DIR)NavierStokes3DThermallyPerfectState.o \
                   $(NAVIERSTOKES_DIR)NavierStokes3DThermallyPerfectHexaBlock.o

# FANS3D:
OBJ_FANS = $(FANS_DIR)FANS3DThermallyPerfectState.o \
           $(FANS_DIR)FANS3DThermallyPerfectHexaBlock.o

# Physics, Turbulence, Thermodynamics, and Chemistry:
OBJ_PHYSICS = $(PHYSICS_DIR)SpeciesInput.o
OBJ_TURBULENCE = $(TURBULENCE_DIR)TurbulenceModellingInput.o \
                 $(TURBULENCE_DIR)TurbulenceModelling.o
OBJ_EOS = $(EOS_DIR)tgas.o \
          $(EOS_DIR)ugas.o
OBJ_NASADATA = $(NASADATA_DIR)NASARP1311data.o
OBJ_REACTIONS = $(REACTIONS_DIR)Reactions.o \
		$(REACTIONS_DIR)15step_GRI211.o \
		$(REACTIONS_DIR)15step_GRI3.o
OBJ_REACTIONS_F90 = $(REACTIONS_DIR)c_15step_GRI211.o \
		    $(REACTIONS_DIR)c_15step_GRI3.o
OBJ_REACTIONS += $(OBJ_REACTIONS_F90)

# TUT testing:
OBJ_TUT := $(SRC_TUT:.cc=.o)

# System:
OBJ_SYSTEM = $(SYSTEM_DIR)System_Linux.o

# Base and CFFC objects:
OBJ_BASE = $(OBJ_MATH) $(OBJ_CFD) $(OBJ_MPI) $(OBJ_SYSTEM) \
           $(OBJ_AMR) $(OBJ_GRID) $(OBJ_ICEM) \
           $(OBJ_MULTIGRID) $(OBJ_NKS)
ifeq ($(TUT_TESTING),ON)
  OBJ_BASE += $(OBJ_TUT)
endif
OBJ_CFFC = $(OBJ_BASE) \
           $(OBJ_EULER) $(OBJ_NAVIERSTOKES) $(OBJ_FANS) \
           $(OBJ_PHYSICS) $(OBJ_TURBULENCE) $(OBJ_EOS) $(OBJ_NASADATA) $(OBJ_REACTIONS)

# Executable objects:
OBJ_MAIN = cffc3D.o

###############################################################################
# Define the executables
###############################################################################
EXE_MAIN = cffc3D

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

all: $(EXE_MAIN)
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

$(EXE_MAIN): CFFC.depend cffc3D.depend $(OBJ_MAIN) $(OBJ_CFFC)
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
	@echo '-------------------------------------------------------------------------'
	@echo ' '
	touch CFFC.depend
	makedepend -fCFFC.depend $(DEPENDSFLAGS) $(CPPFLAGS) $(SRC_CFFC) >& /dev/null
	@echo ' '
	@echo CFFC source file dependencies have been brought up to date.
	@echo '-----------------------------------------------------------------'
	@echo ' '

cffc3D.depend: $(SRC_MAIN)
	@echo ' '
	@echo Determining cffc3D source file dependencies for a $(PLATFORM) platform.
	@echo '--------=--------------------------------------------------------------'
	@echo ' '
	touch cffc3D.depend
	makedepend -fcffc3D.depend $(DEPENDSFLAGS) $(CPPFLAGS) $(SRC_MAIN) >& /dev/null
	@echo ' '
	@echo cffc3D source file dependencies have been brought up to date.
	@echo '-----------------------------------------------------------------'
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

install-all: $(EXE_MAIN)
	@for f in $(EXE_MAIN); \
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
	@for f in $(EXE_MAIN); \
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
	rm -f $(CFD_DIR)*.o $(CFD_DIR)*.il $(CFD_DIR)*~
	rm -f $(MPI_DIR)*.o $(MPI_DIR)*.il $(MPI_DIR)*~
	rm -f $(GRID_DIR)*.o $(GRID_DIR)*.il $(GRID_DIR)*~
	rm -f $(ICEM_DIR)*.o $(ICEM_DIR)*.il $(ICEM_DIR)*~ $(ICEM_DIR)*.depend $(ICEM_DIR)*.bak
	rm -f $(AMR_DIR)*.o $(AMR_DIR)*.il $(AMR_DIR)*~
	rm -f $(EULER_DIR)*.o $(EULER_DIR)*.il $(EULER_DIR)*~ $(EULER_DIR)*.depend $(EULER_DIR)*.bak
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
	rm -f $(MULTIGRID_DIR)*.o $(MULTIGRID_DIR)*.il $(MULTIGRID_DIR)*~
	rm -f $(NKS_DIR)*.o $(NKS_DIR)*.il $(NKS_DIR)*~
	rm -f $(SYSTEM_DIR)*.o $(SYSTEM_DIR)*.il $(SYSTEM_DIR)*~
	rm -f $(MODULE_DIR)/*.mod $(MODULE_DIR)/*~
	rm -f $(EXE_MAIN)
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

#
# Include header file dependencies
include CFFC.depend
include cffc3D.depend


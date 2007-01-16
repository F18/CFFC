#####################################################
# Include Reconstruction platform and compiler flags
#####################################################
#
Reconstruction_Path:= $(shell cd ../../../; pwd)
include $(Reconstruction_Path)/src/makefile.def

#overide for linking
F77=g77
FFLAGS=

BPKIT_LIB = ../lib/libbpkit.a

OBJ = BpFortran.o \
      BpC.o \
      FGMRES.o \
      HBTMat.o \
      BlockMat.o \
      BRelax.o \
      BILUK.o \
      BTIF.o \
      BlockVec.o \
      CSRMat.o \
      CSRmath.o \
      CSRapinvs.o \
      CSRgmres.o \
      CSRilut.o \
      SparseUtil.o \
      DenseMat.o \
      BpResource.o \
      ilut.o \
      readhb.o \
      spmm.o \
      spsm.o

lib $(BPKIT_LIB): $(OBJ)
	ar rcv $(BPKIT_LIB) $(OBJ)
#	ranlib $(BPKIT_LIB)
# above ranlib needed on sun4

.SUFFIXES: .cc .f
.cc.o:
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $*.cc
.f.o:
	$(F77) $(FFLAGS) -c $*.f

clean:
	rm -f *.o $(BPKIT_LIB) *.il

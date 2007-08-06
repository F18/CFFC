#
# Include CFFC platform and compiler flags
#
CFFC_Path:= $(shell cd ../../; pwd)
CFFC_SRC_TREE = 2D

ifeq ($(CFFC_SRC_TREE),2D)
   include $(CFFC_Path)/src_2D/makefile.def
endif
ifeq ($(CFFC_SRC_TREE),3D)
   include $(CFFC_Path)/src_3D/makefile.def
endif

#
# Set Linux g++ flags
#
ifeq ($(PLATFORM),Linux_i386)
  ifeq ($(GCC_VERSION),V296)
    CXXFLAGS = -O3 '-D_GNU_GCC_296'
  endif
  ifeq ($(GCC_VERSION),V3+)
    CXXFLAGS = -O3 '-D_GNU_GCC_3'
  endif
endif

#
# Overide F77 flags for linking bpkit
#
#F77 = gfortran
#ifeq ($(PLATFORM),Linux_alpha)
#   F77 = g77
#endif
#FFLAGS = -O3

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
	@echo $(CFFC_Path)
	rm -f *.o $(BPKIT_LIB) *.il

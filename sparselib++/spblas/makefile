#
# Makefile for sparselib spblas
#

CFFC_Path:= $(shell cd ../../; pwd)
CFFC_SRC_TREE = 2D

include ../makefile.def

SRCS 	= spmm.cc spsm.cc
OBJS	= spmm.o  spsm.o

all: lib

lib: $(SPARSEBLAS_A)

$(SPARSEBLAS_A) : $(OBJS)
	@ ar rv $(SPARSEBLAS_A) $?
	@ case x$(HASRANLIB) in xt ) echo ranlib; ranlib $(SPARSEBLAS_A);; esac

clean: 
	/bin/rm -f $(OBJS) *.il

wipe:
	/bin/rm -f $(OBJS) $(SPARSEBLAS_A) *.il

#spblas.$(OBJ_EXT) : spblas.cc \
#	$(SPARSELIB_INCLUDE_DIR)/spblas.h

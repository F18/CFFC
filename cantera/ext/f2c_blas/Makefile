#/bin/sh
#/bin/sh
#
# $Source: /cvsroot/cantera/cantera/ext/f2c_blas/Makefile.in,v $
# $Author: dggoodwin $
# $Revision: 1.4 $
# $Date: 2005/06/18 17:01:33 $
# $License$
#

.SUFFIXES :
.SUFFIXES : .c .d .o

do_ranlib = 1

# the directory where the Cantera libraries are located
CANTERA_LIBDIR=/nfs/kris/d1/people/charestm/Programs/CFFC/cantera/build/lib/i686-pc-linux-gnu

# the directory where Cantera include files may be found.
CANTERA_INCDIR=/nfs/kris/d1/people/charestm/Programs/CFFC/cantera/build/include/cantera

# the C++ compiler
CXX = g++ -O

# the C compiler
CC = gcc

# C++ compile flags
CXX_FLAGS = -D_NO_MPI_VERSION -D_GNU_GCC_3 -D_NO_ICEMCFD_VERSION -I/nfs/kris/d1/people/charestm/Programs/CFFC/sparselib++/include -I/nfs/kris/d1/people/charestm/Programs/CFFC/sparselib++/iml/include -I/nfs/kris/d1/people/charestm/Programs/CFFC/sparselib++/mv/include -I/nfs/kris/d1/people/charestm/Programs/CFFC/bpkit/src  $(CXX_OPT)

# Local include files
CXX_INCLUDES=-I../f2c_libs

# How to compile the dependency file
.c.d:
	g++ -MM  $(CXX_FLAGS) $(CXX_INCLUDES) $*.c > $*.d

# How to compile a C file
.c.o:
	gcc -c $< -DHAVE_CONFIG_H $(CXX_FLAGS) $(CXX_INCLUDES)


# -----------------------------------------------

BLASLIB = /nfs/kris/d1/people/charestm/Programs/CFFC/cantera/build/lib/i686-pc-linux-gnu/libctblas.a

all: $(BLASLIB)

OBJS = \
dasum.o \
daxpy.o \
dcabs1.o \
dcopy.o \
ddot.o \
dgbmv.o \
dgemm.o \
dgemv.o \
dger.o \
dnrm2.o \
drot.o \
drotg.o \
drotm.o \
drotmg.o \
dsbmv.o \
dscal.o \
dsdot.o \
dspmv.o \
dspr.o \
dspr2.o \
dswap.o \
dsymm.o \
dsymv.o \
dsyr.o \
dsyr2.o \
dsyr2k.o \
dsyrk.o \
dtbmv.o \
dtbsv.o \
dtpmv.o \
dtpsv.o \
dtrmm.o \
dtrmv.o \
dtrsm.o \
dtrsv.o \
dzasum.o \
dznrm2.o \
idamax.o \
lsame.o \
xerbla.o 

SRCS = $(OBJS:.o=.cpp)

# List of dependency files to be created
DEPENDS=$(OBJS:.o=.d)

# rule to make library
$(BLASLIB):  $(OBJS)
	ar ruv $(BLASLIB) $(OBJS) > /dev/null
ifeq ($(do_ranlib),1)
	ranlib $(BLASLIB)
endif


# ------------------------------------------------
# Utility Targets

clean:
	$(RM) $(OBJS) $(BLASLIB) *.d .depends

# depends target
depends:
	$(RM) *.d .depends
	make .depends

.depends: $(DEPENDS)
	cat *.d > .depends

ifeq ($(wildcard .depends), .depends)
include .depends
endif


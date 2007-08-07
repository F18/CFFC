#/bin/sh
#
# $Source: /cvsroot/cantera/cantera/ext/f2c_math/Makefile.in,v $
# $Author: dggoodwin $
# $Revision: 1.3 $
# $Date: 2005/06/18 17:01:33 $
#

.SUFFIXES :
.SUFFIXES : .c .cpp .d .o

do_ranlib = 1
# the directory where the Cantera libraries are located
CANTERA_LIBDIR=/home/groth/CFFC/cantera/build/lib/i686-pc-linux-gnu

# the directory where Cantera include files may be found.
CANTERA_INCDIR=/home/groth/CFFC/cantera/build/include/cantera

# the C++ compiler
CXX = g++ -O

# the C compiler
CC = gcc

# C++ compile flags
CXX_FLAGS = -D_NO_MPI_VERSION -D_GNU_GCC_3 -D_NO_ICEMCFD_VERSION -I/home/groth/CFFC/sparselib++/include -I/home/groth/CFFC/sparselib++/iml/include -I/home/groth/CFFC/sparselib++/mv/include -I/home/groth/CFFC/bpkit/src  $(CXX_OPT)

# Local include files
CXX_INCLUDES=-I../f2c_libs

# How to compile the dependency file
.c.d:
	g++ -MM  $(CXX_FLAGS) $(CXX_INCLUDES) $*.c > $*.d

.cpp.d:
	g++ -MM  $(CXX_FLAGS) $(CXX_INCLUDES) $*.cpp > $*.d

# How to compile a C file
.c.o:
	gcc -c $< -DHAVE_CONFIG_H $(CXX_FLAGS) $(CXX_INCLUDES)

# How to compile a Cpp file
.cpp.o:
	g++ -O -c $< -DHAVE_CONFIG_H $(CXX_FLAGS) $(CXX_INCLUDES)


# -----------------------------------------------

LIB = /home/groth/CFFC/cantera/build/lib/i686-pc-linux-gnu/libctmath.a

all: $(LIB)

OBJS = \
mach.o \
ddaspk.o \
dgbefa.o \
dgbsl.o \
dgefa.o \
dgesl.o \
dp1vlu.o \
dpcoef.o \
dpolft.o \
fdump.o \
j4save.o \
pcoef.o \
polfit.o \
pvalue.o \
xercnt.o \
xerhlt.o \
xermsg.o \
xerprn.o \
xersve.o \
xgetua.o 

SRCS = $(OBJS:.o=.cpp)

# List of dependency files to be created
DEPENDS=$(OBJS:.o=.d)

# How to make the static library
$(LIB):  $(OBJS)
	ar ruv $(LIB) $(OBJS)
ifeq ($(do_ranlib),1)
	ranlib $(LIB)
endif

# ------------------------------------------------
# Utility Targets

clean:
	$(RM) $(OBJS) $(LIB) *.d .depends

# depends target
depends:
	$(RM) *.d .depends
	make .depends

.depends: $(DEPENDS)
	cat *.d > .depends

ifeq ($(wildcard .depends), .depends)
include .depends
endif


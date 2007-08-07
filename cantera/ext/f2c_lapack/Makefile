#/bin/sh
#
# $Source: /cvsroot/cantera/cantera/ext/f2c_lapack/Makefile.in,v $
# $Author: hkmoffa $
# $Revision: 1.3 $
# $Date: 2006/04/30 22:32:55 $
#

.SUFFIXES :
.SUFFIXES : .c .o .d

do_ranlib = 1

# the directory where the Cantera libraries are located
CANTERA_LIBDIR=/nfs/kris/d1/people/charestm/Programs/CFFC/cantera/build/lib/i686-pc-linux-gnu

# the C++ compiler
CXX = g++ -O

# the C compiler
CC = gcc

# C++ compile flags
CXX_FLAGS = -D_NO_MPI_VERSION -D_GNU_GCC_3 -D_NO_ICEMCFD_VERSION -I/nfs/kris/d1/people/charestm/Programs/CFFC/sparselib++/include -I/nfs/kris/d1/people/charestm/Programs/CFFC/sparselib++/iml/include -I/nfs/kris/d1/people/charestm/Programs/CFFC/sparselib++/mv/include -I/nfs/kris/d1/people/charestm/Programs/CFFC/bpkit/src  

# Local include files
CXX_INCLUDES=-I../f2c_libs

# How to compile the dependency file
.c.d:
	g++ -MM  $(CXX_FLAGS) $(CXX_INCLUDES) $*.c > $*.d

.c.o:
	gcc -c $< $(CXX_FLAGS) $(CXX_INCLUDES)

# -----------------------------------------------

LAPACKLIB = /nfs/kris/d1/people/charestm/Programs/CFFC/cantera/build/lib/i686-pc-linux-gnu/libctlapack.a

all: $(LAPACKLIB)

OBJS = \
dbdsqr.o \
dgbtrf.o \
dgbtf2.o \
dgbtrs.o \
dgbsv.o \
dgebd2.o \
dgebrd.o \
dgelq2.o \
dgelqf.o \
dgelss.o \
dgeqr2.o \
dgeqrf.o \
dgetf2.o \
dgetrf.o \
dgetri.o \
dgetrs.o \
dlabad.o \
dlabrd.o \
dlacpy.o \
dlamch.o \
dlange.o \
dlapy2.o \
dlarf.o \
dlarfb.o \
dlarfg.o \
dlarft.o \
dlartg.o \
dlas2.o \
dlascl.o \
dlaset.o \
dlasq1.o \
dlasq2.o \
dlasq3.o \
dlasq4.o \
dlasq5.o \
dlasq6.o \
dlasr.o \
dlasrt.o \
dlassq.o \
dlasv2.o \
dlaswp.o \
dorg2r.o \
dorgbr.o \
dorgl2.o \
dorglq.o \
dorgqr.o \
dorm2r.o \
dormbr.o \
dorml2.o \
dormlq.o \
dormqr.o \
drscl.o \
dtrtri.o \
dtrti2.o \
ieeeck.o \
ilaenv.o 


# list of source files
SRCS = $(OBJS:.o=.c)

# List of dependency files to be created
DEPENDS=$(OBJS:.o=.d)


$(LAPACKLIB):  $(OBJS)
	ar ruv $(LAPACKLIB) $(OBJS) > /dev/null
ifeq ($(do_ranlib),1)
	ranlib $(LAPACKLIB)
endif
# ------------------------------------------------
# Utility Targets

clean:
	$(RM) $(OBJS) $(LAPACKLIB) *.d .depends


# depends target
depends:
	$(RM) *.d .depends
	make .depends

.depends: $(DEPENDS)
	cat *.d > .depends
 
ifeq ($(wildcard .depends), .depends)
include .depends
endif






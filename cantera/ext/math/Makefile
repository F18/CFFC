#
# $Source: /cvsroot/cantera/cantera/ext/math/Makefile.in,v $
# $Author: hkmoffa $
# $Revision: 1.8 $
# $Date: 2006/03/03 20:31:49 $
#

.SUFFIXES :
.SUFFIXES : .f .d .o .cpp

do_ranlib = 1

LIB = /home/groth/CFFC/cantera/build/lib/i686-pc-linux-gnu/libctmath.a

all: $(LIB)


F_FLAGS = -O  -fno-second-underscore 

OBJS = \
mach.o \
ddaspk.o \
dgbfa.o \
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
all = $(LIB)

$(LIB):  $(OBJS)
	ar ruv $(LIB) $(OBJS) > /dev/null
ifeq ($(do_ranlib),1)
	ranlib $(LIB)
endif
%.o : %.cpp
	g++ -O -c $< -DHAVE_CONFIG_H  -D_NO_MPI_VERSION -D_GNU_GCC_3 -D_NO_ICEMCFD_VERSION -I/home/groth/CFFC/sparselib++/include -I/home/groth/CFFC/sparselib++/iml/include -I/home/groth/CFFC/sparselib++/mv/include -I/home/groth/CFFC/bpkit/src  $(INCLUDES)

%.o : %.f
	gfortran -c $< $(F_FLAGS)

clean:
	$(RM) $(OBJS) $(LIB)

depends:
	echo '...'

 




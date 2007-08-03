# $License$
#
# $Id: Makefile.in,v 1.6 2006/03/03 22:44:43 hkmoffa Exp $
#
#/bin/sh
.SUFFIXES :
.SUFFIXES : .f .d .o

LIB = /home/groth/CFFC/cantera/build/lib/i686-pc-linux-gnu/librecipes.a

do_ranlib = 1

F_FLAGS = -O  -fno-second-underscore 

OBJS = simp1.o simp2.o simp3.o simplx.o 
SRCS = $(OBJS:.o=.cpp)


all:  $(LIB)

$(LIB): $(OBJS)
	ar ruv $(LIB) $(OBJS) > /dev/null
ifeq ($(do_ranlib),1)
	ranlib $(LIB)
endif

%.o : %.f
	gfortran -c $< $(F_FLAGS)

clean:
	$(RM) $(OBJS) $(LIB)

depends: 




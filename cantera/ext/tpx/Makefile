#/bin/sh

.SUFFIXES :
.SUFFIXES : .cpp .d .o

do_ranlib = 1
CXX_FLAGS = -D_NO_MPI_VERSION -D_GNU_GCC_3 -D_NO_ICEMCFD_VERSION -I/home/groth/CFFC/sparselib++/include -I/home/groth/CFFC/sparselib++/iml/include -I/home/groth/CFFC/sparselib++/mv/include -I/home/groth/CFFC/bpkit/src  $(CXX_OPT)

#COBJS = tpx_files.o
COBJS = Methane.o Nitrogen.o Oxygen.o Water.o Hydrogen.o RedlichKwong.o \
        CarbonDioxide.o Heptane.o lk.o Sub.o utils.o HFC134a.o

FOBJS = 

CXX_LIBS = 
CXX_INCLUDES = -I../include 
TPLIB = /home/groth/CFFC/cantera/build/lib/i686-pc-linux-gnu/libtpx.a

DEPENDS = $(COBJS:.o=.d)

all:  $(TPLIB)

%.d:
	g++ -MM $(CXX_INCLUDES) $*.cpp > $*.d

.cpp.o:
	g++ -O -c $< -DHAVE_CONFIG_H $(CXX_FLAGS) $(CXX_INCLUDES)

.f.o:
	gfortran -c $< $(F77_FLAGS) 

$(TPLIB): $(COBJS) $(FOBJS)
	ar ruv $(TPLIB) $(COBJS) $(FOBJS) > /dev/null
ifeq ($(do_ranlib),1)
	ranlib $(TPLIB)
endif

clean:
	$(RM) $(COBJS) $(FOBJS) *~ *.d .depends

depends: $(DEPENDS)
	cat *.d > .depends
	$(RM) $(DEPENDS)

ifeq ($(wildcard .depends), .depends)
include .depends
endif


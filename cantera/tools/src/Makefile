#/bin/sh

LIBDIR = /home/groth/CFFC/cantera/build/lib/i686-pc-linux-gnu
INCDIR = -I/home/groth/CFFC/cantera/Cantera/src
BINDIR = /home/groth/CFFC/cantera/build/bin/i686-pc-linux-gnu
build_ck = 1

LCXX_FLAGS = -L$(LIBDIR)   -L/home/groth/CFFC/cantera/build/lib/i686-pc-linux-gnu -D_NO_MPI_VERSION -D_GNU_GCC_3 -D_NO_ICEMCFD_VERSION -I/home/groth/CFFC/sparselib++/include -I/home/groth/CFFC/sparselib++/iml/include -I/home/groth/CFFC/sparselib++/mv/include -I/home/groth/CFFC/bpkit/src 
LOCAL_LIBS = -lcantera -ltpx -lctcxx 
#LOCAL_LIBS = -lcantera -lcvode -lctmath -lctlapack -lctblas -lctcxx 

LCXX_END_LIBS =  -lctf2c -lm 

OBJS = ck2cti.o cti2ctml.o fixtext.o

DEPENDS = $(OBJS:.o=.d)

CONVLIB_DEP = /home/groth/CFFC/cantera/build/lib/i686-pc-linux-gnu/libconverters.a
CANTERALIB_DEP = /home/groth/CFFC/cantera/build/lib/i686-pc-linux-gnu/libcantera.a \
	         /home/groth/CFFC/cantera/build/lib/i686-pc-linux-gnu/libctcxx.a

progs = $(BINDIR)/cti2ctml $(BINDIR)/fixtext
ifeq ($(build_ck),1) 
exes =  $(progs) $(BINDIR)/ck2cti
else
exes = $(progs)
endif

.cpp.o:
	g++ -O -c $< -DHAVE_CONFIG_H  $(INCDIR) -D_NO_MPI_VERSION -D_GNU_GCC_3 -D_NO_ICEMCFD_VERSION -I/home/groth/CFFC/sparselib++/include -I/home/groth/CFFC/sparselib++/iml/include -I/home/groth/CFFC/sparselib++/mv/include -I/home/groth/CFFC/bpkit/src  $(CXX_FLAGS) 

all: $(exes)

ck: $(BINDIR)/ck2cti

$(BINDIR)/ck2cti: ck2cti.o $(CONVLIB_DEP) $(CANTERALIB_DEP)
	$(RM) $(BINDIR)/ck2cti
	g++ -O -o $(BINDIR)/ck2cti ck2cti.o $(LCXX_FLAGS) \
               -lconverters -lcantera -ltpx -lctcxx $(LCXX_END_LIBS)

$(BINDIR)/cti2ctml: cti2ctml.o $(CANTERALIB_DEP)
	$(RM) $(BINDIR)/cti2ctml
	g++ -O -o $(BINDIR)/cti2ctml cti2ctml.o $(LCXX_FLAGS) $(LOCAL_LIBS) \
               $(LCXX_END_LIBS)

$(BINDIR)/fixtext: fixtext.o 
	$(RM) $(BINDIR)/fixtext
	g++ -O -o $(BINDIR)/fixtext fixtext.o $(LCXX_FLAGS) \
               $(LCXX_END_LIBS)

$(CONVLIB_DEP):
	cd ../../Cantera/src/converters; make

clean:
	$(RM) *.o *.*~ $(BINDIR)/cti2ctml $(BINDIR)/ck2cti
	(if test -d SunWS_cache ; then \
           $(RM) -rf SunWS_cache ; \
        fi )

%.d:
	g++ -O -MM $(INCDIR) -D_NO_MPI_VERSION -D_GNU_GCC_3 -D_NO_ICEMCFD_VERSION -I/home/groth/CFFC/sparselib++/include -I/home/groth/CFFC/sparselib++/iml/include -I/home/groth/CFFC/sparselib++/mv/include -I/home/groth/CFFC/bpkit/src  $(CXX_FLAGS) $*.cpp > $*.d

depends: $(DEPENDS)
	cat *.d &> .depends
	$(RM) $(DEPENDS)

ifeq ((test -e .depends),0)
include .depends
endif

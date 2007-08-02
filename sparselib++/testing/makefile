# 
# Makefile for sparselib testing
#

CFFC_Path:= $(shell cd ../../; pwd)
CFFC_SRC_TREE = 2D

include ../makefile.def

SPARSELIB_LIB             = -L$(SPARSELIB_LIB_DIR) -lsparse -lspblas -lmv 

BLAS_LIB                  = -lblas -llapack -lg2c -lgcc

LDFLAGS2                  = $(BLAS_LIB) $(SPARSELIB_LIB) -lm

SRCS		= tsl.cc tpre.cc tspsm.cc tsp.cc
EXECS		= tsl    tpre    tspsm    tsp
OBJS		= tsl.o  tpre.o  tspsm.o  tsp.o

MYEXEC          = mytest
MYSRC           = mytest.cc
MYOBJ           = mytest.o

all: progs

progs: $(EXECS)

myprog: $(MYEXEC)

clean:
	/bin/rm -f $(OBJS) $(EXECS) out *.out

wipe: 
	/bin/rm -f $(OBJS) $(EXECS) out *.out

tsl.o : tsl.cc  
tpre.o : tpre.cc  
tspsm.o : tspsm.cc  
tsp.o : tsp.cc

mytest.o: mytest.cc

tsl : tsl.o
	$(CCC) $(CCCFLAGS) -o tsl tsl.o  $(LDFLAGS)

tpre : tpre.o
	$(CCC) $(CCCFLAGS) -o tpre tpre.o  $(LDFLAGS)

tspsm : tspsm.o
	$(CCC) $(CCCFLAGS) -o tspsm tspsm.o  $(LDFLAGS)

tsp : tsp.o
	$(CCC) $(CCCFLAGS) -o tsp tsp.o  $(LDFLAGS)

mytest : mytest.o
	$(CCC) $(CCCFLAGS) -o mytest mytest.o  $(LDFLAGS2)

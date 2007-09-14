CCC = g++
CC  = gcc
FC  = f77
INCDIRS = -I../src
CFLAGS  = $(INCDIRS) -O3
FFLAGS  = -O3
LDFLAGS = 

###### Solaris section
LALIB  = -L/project/darpa/lib/solaris -llapack -lblas
LALIB  = -llapack -lblas
FC_LIB = -L/opt/SUNWspro/SC4.2/lib -lF77 -lM77 -lsunmath

###### Sun 4 section
#LALIB  = -L/project/darpa/lib/sun4 -llapack -lblas
#FC_LIB = -L/opt/lang/SC2.0.1 -lF77 -lM77 -lansi /opt/lang/SC2.0.1/values-Xa.o
#FFLAGS = -O3 -dalign

###### IRIX 32-bit section
#LALIB  = -L/project/darpa/lib/irix -llapack -lblas
#FC_LIB = -L/usr/lib -lF77 -lI77 -lisam

###### DEC section
#LALIB  = /usr/opt/XMDLOA320/dxml/libdxml_ev4.a
#FC_LIB = -lc -lUfor -lfor -lFutil -lutil

###### Linux section
#LALIB  = -L/project/darpa2/lapack-linux-2.0.1-elf -llapack -lblas
#LALIB  = -llapack -lblas
#FC_LIB = -L/opt/gnu/lib -lf2c -lgcc
#FC_LIB = -L/usr/lib/gcc-lib/i386-redhat-linux/egcs-2.91.66 -lg2c -lgcc

###### Cray section
#CCC = CC
#CC = cc
#FC = cf77
#LALIB   =
#FC_LIB  =
#CFLAGS  = $(INCDIRS) -DCRAY
#FFLAGS  =
#LDFLAGS = -Wl"-i equiv"

LIBS = ../src/libbpkit.a $(LALIB) $(FC_LIB) -lm

# use C++ compiler to link, or else you must provide locations of C++ libraries

EXES = bp simple f1example f2example f3example c_example b22_test

all: $(EXES)

bp: bp.o ../src/libbpkit.a
	$(CCC) $(LDFLAGS) -o bp bp.o $(LIBS)

simple: simple.o ../src/libbpkit.a
	$(CCC) $(LDFLAGS) -o simple simple.o $(LIBS)

f1example: f1example.o ../src/libbpkit.a
	$(CCC) $(LDFLAGS) -o f1example f1example.o $(LIBS)

f2example: f2example.o ../src/libbpkit.a
	$(CCC) $(LDFLAGS) -o f2example f2example.o $(LIBS)

f3example: f3example.o ../src/libbpkit.a
	$(CCC) $(LDFLAGS) -o f3example f3example.o $(LIBS)

c_example: c_example.o ../src/libbpkit.a
	$(CCC) $(LDFLAGS) -o c_example c_example.o $(LIBS)

b22_test: b22_test.o B22.o ../src/libbpkit.a
	$(CCC) $(LDFLAGS) -o b22_test b22_test.o B22.o $(LIBS)

clean:
	rm -f *.o core $(EXES)

.SUFFIXES: .cc .c .f .o
.cc.o:
	$(CCC) $(CFLAGS) -c $<
.c.o:
	$(CC) $(CFLAGS) -c $<
.f.o:
	$(FC) $(FFLAGS) -c $<

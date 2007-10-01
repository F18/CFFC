#
# Makefile for tgas testing.
#
OBJ1 = tgas_test.o
OBJ2 = ugas_test.o
OBJ3 = tgas.o ugas.o
EXE1 = tgas_test
EXE2 = ugas_test
FORT = gfortran
BIND = gfortran
CPU = 
OPT = 
DEBUG =
#
# Define the dependencies for the main programs.
#
all : $(EXE1) $(EXE2)
	@echo Programs TGAS_TEST and UGAS_TEST have been brought up to date.

tgas : $(EXE1)
	@echo Program TGAS_TEST has been brought up to date.

ugas : $(EXE2)
	@echo Program UGAS_TEST has been brought up to date.
#
# Define the dependencies for the main executables.
#
$(EXE1) : $(OBJ1) $(OBJ3)
	$(BIND) -o $(EXE1) $(OBJ1) $(OBJ3)

$(EXE2) : $(OBJ2) $(OBJ3)
	$(BIND) -o $(EXE2) $(OBJ2) $(OBJ3)
#
# Define the dependencies for each object.
#
tgas_test.o : tgas_test.f
	$(FORT) -c $(CPU) $(OPT) $(DEBUG) tgas_test.f

ugas_test.o : ugas_test.f
	$(FORT) -c $(CPU) $(OPT) $(DEBUG) ugas_test.f

tgas.o : tgas.f
	$(FORT) -c $(CPU) $(OPT) $(DEBUG) tgas.f

ugas.o : ugas.f
	$(FORT) -c $(CPU) $(OPT) $(DEBUG) ugas.f

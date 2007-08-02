SOURCEDIR = /usr2/scratch/jai/Results/RinglebFlow/EmbeddedInterfaces/Data/src/

OBJ1 = order.o 

EXE = order

CC = g++
BIND = g++
CPU = 
OPT = 
DEBUG = -Wall
COMPUTER = 1
LIB = -lm

# Define the dependencies for the main programs.
all : $(EXE)
	@echo Program has been brought up to date.

# Define the dependencies for the main executables.
$(EXE) : $(OBJ1) 
	 $(BIND) -o $(EXE) $(OBJ1) $(LIB)

# Define the dependencies for each object.
order.o : order.cc
	$(CC) -c $(CPU) $(OPT) $(DEBUG) -D MACHINE=$(COMPUTER) order.cc

# Clean up the home directory.
clean:
	rm -f *.o
	rm -f *~

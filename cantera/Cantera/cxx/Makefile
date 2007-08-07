#/bin/sh
###############################################################
#  $Author: dggoodwin $
#  $Date: 2006/10/25 21:27:59 $
#  $Revision: 1.9 $
#
#  Copyright 2001 California Institute of Technology
#  See file License.txt for licensing information
#
###############################################################

build_f90=1
build_python=0

INCDIR =      ../../build/include/cantera
INSTALL_TSC = ../../../bin/install_tsc

CXX_H = Cantera.h  equilibrium.h  IncompressibleSolid.h \
        kinetics.h onedim.h  surface.h GRI30.h  integrators.h \
        Metal.h   PureFluid.h   transport.h Edge.h  \
        IdealGasMix.h  Interface.h    numerics.h  \
        reactionpaths.h  zerodim.h importPhase.h thermo.h

all:
	@(cd include ; \
          for lh in $(CXX_H) ; do \
          $(INSTALL_TSC) "$${lh}" "../"$(INCDIR) ; \
	done)
	cd src;   make
#	cd demos; make

clean:
	@(for lh in $(CXX_H) ; do \
            th=$(INCDIR)/"$${lh}" ; \
            if test -f "$${th}" ; then \
              $(RM) "$${th}"  ; \
              echo "$(RM) $${th}" ; \
            fi \
          done)
	cd src; make clean
	cd demos; make clean

depends:
	cd src; make depends
#	cd demos; make depends

install:
	cd src; make install
	cd demos; make install

# end of file


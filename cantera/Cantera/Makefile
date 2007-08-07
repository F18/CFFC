#/bin/sh
###############################################################
#  $Author: hkmoffa $
#  $Date: 2006/08/21 16:44:02 $
#  $Revision: 1.7 $
#
#  Copyright 2001 California Institute of Technology
#  See file License.txt for licensing information
#
###############################################################

build_f90=1
build_python=0
build_matlab = 0

all:
	cd src; make
	cd cxx; make
	cd clib/src; make
ifeq ($(build_f90),1)
	cd fortran/src; make
endif
	cd user; make

clean:
	cd src; make clean
	cd cxx; make clean
	cd clib/src; make clean
	cd python; make clean
ifeq ($(build_f90),1)
	cd fortran/src; make clean
endif
	cd user; make clean

depends:
	cd src; make depends
	cd cxx/src; make depends
	cd clib/src; make depends
ifeq ($(build_f90),1)
	cd fortran/src; make depends
endif
ifeq ($(build_matlab),1)
	cd matlab; make depends
endif
	cd user; make depends

install:
	cd src; make install
	cd cxx/src; make install
	cd clib/src; make install
ifeq ($(build_f90),1)
	cd fortran/src; make install
endif
	cd user; make install

# end of file


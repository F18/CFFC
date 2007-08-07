#/bin/sh
###############################################################
#  $Author: hkmoffa $
#  $Date: 2006/04/30 23:10:42 $
#  $Revision: 1.23 $
#
#  Copyright 2001 California Institute of Technology
#  See file License.txt for licensing information
#
###############################################################

have_python_site_package_topdir=0
python_site_package_topdir=/usr/local/cantera

CANTERA_LIBDIR= /nfs/kris/d1/people/charestm/Programs/CFFC/cantera/build/lib/i686-pc-linux-gnu
LIB_DEPS = $(CANTERA_LIBDIR)/libcantera.a $(CANTERA_LIBDIR)/libzeroD.a \
           $(CANTERA_LIBDIR)/liboneD.a \
           $(CANTERA_LIBDIR)/libtransport.a \
           $(CANTERA_LIBDIR)/libclib.a \
           $(CANTERA_LIBDIR)/libconverters.a 

WIN_LIB_DEPS = $(CANTERA_LIBDIR)/cantera.lib $(CANTERA_LIBDIR)/zeroD.lib \
           $(CANTERA_LIBDIR)/oneD.lib \
           $(CANTERA_LIBDIR)/transport.lib
SRCS = src/ctphase_methods.cpp \
      src/ctthermo_methods.cpp \
      src/ctkinetics_methods.cpp \
      src/cttransport_methods.cpp \
      src/ctxml_methods.cpp \
      src/ctfuncs.cpp \
      src/ctsurf_methods.cpp \
      src/ctbndry_methods.cpp \
      src/ctrpath_methods.cpp \
      src/ctreactor_methods.cpp \
      src/ctfunc_methods.cpp \
      src/ctonedim_methods.cpp \
      src/methods.h

os_is_win = 0
use_clib_dll = 0

all: _build 

win: _winbuild

#
# HKM -> Python seems to want to compile C++ code with the CC compiler.
#        The fix is to assign the CC compiler to the CXX compiler. 
#        the double quotes are needed if the CXX variable is multi-token.
#
_build: $(SRCS) $(LIB_DEPS) Makefile setup.py
	touch src/pycantera.cpp 
	(CXX="g++ -O"; export CXX; CC="g++ -O"; export CC; /usr/bin/python2 setup.py build)
	echo 'ok' > _build

#
# HKM -> If clib is built as a dll, it needs to be copied by hand to the
#        site package directory.
#
ifeq ($(use_clib_dll), 1)
CLIB_DLL=../../build/lib/i686-pc-win32/clib.dll
CLIB_EXP=../../build/lib/i686-pc-win32/clib.exp
else
CLIB_DLL=
CLIB_EXP=
endif

_winbuild: $(SRCS) $(WIN_LIB_DEPS) $(CLIB_DLL)
	touch src/pycantera.cpp 
	(/usr/bin/python2 setup.py build)
ifeq ($(use_clib_dll), 1)
ifeq ($(have_python_site_package_topdir), 1)
	(/usr/bin/install -c -m 755 $(CLIB_DLL) $(CLIB_EXP)  \
            $(python_site_package_topdir)/lib/site-packages/Cantera)
endif
endif
	echo 'ok' > _winbuild

minbuild:
	(/usr/bin/python2 setup.py build)

install:
ifeq (0,1)
	(/usr/bin/python2 setup.py install --prefix="/usr/local/cantera")
else
	(/usr/bin/python2 setup.py install)
endif

clean:
	/usr/bin/python2 setup.py clean
	rm -f _build; rm -f _winbuild
	cd src; rm -f *.o

depends:
	echo '-'

test:
	cd examples; /usr/bin/python2 reactor1.py &> reactor1.out
	cd examples; /usr/bin/python2 flame1.py  &> flame1.out
	cd examples; /usr/bin/python2 diamond.py &> diamond.out
	cd examples; /usr/bin/python2 critProperties.py  &> critProperties.out

# end of file




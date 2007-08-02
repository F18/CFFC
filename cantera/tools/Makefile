#
# $Revision: 1.11 $
# $Author: dggoodwin $
# $Date: 2004/08/05 14:26:42 $
#
build_ck=1

all: apps

install: templates-install apps-install man-install

apps:
	cd src; make
	cd testtools; make

apps-install:
	/usr/bin/install -c -d /usr/local/cantera/bin
ifeq ($(build_ck),1)
	/usr/bin/install -c /home/groth/CFFC/cantera/build/bin/i686-pc-linux-gnu/ck2cti /usr/local/cantera/bin
endif
	/usr/bin/install -c /home/groth/CFFC/cantera/build/bin/i686-pc-linux-gnu/cti2ctml /usr/local/cantera/bin

man-install:
	/usr/bin/install -c -d /usr/local/cantera/man1
	/usr/bin/install -c man/ck2cti.1 /usr/local/cantera/man1
	/usr/bin/install -c man/mixmaster.1 /usr/local/cantera/man1

uninstall:
	rm -r -f /usr/local/cantera/man1/ck2cti.1
	rm -f /usr/local/cantera/bin/ck2cti

clean:
	$(RM) *.*~ 
	cd src; make clean
	cd testtools; make clean

docs:
	cd doc; doxygen Cantera.cfg

depends:
	cd src; make depends
	cd testtools; make depends

templates-install:
	/usr/bin/install -c -d /usr/local/cantera/templates
	/usr/bin/install -c -d /usr/local/cantera/templates/f77
	/usr/bin/install -c -m 644 templates/f77/*.cpp /usr/local/cantera/templates/f77
	/usr/bin/install -c -m 644 templates/f77/*.f /usr/local/cantera/templates/f77
	/usr/bin/install -c -m 644 templates/f77/*.mak /usr/local/cantera/templates/f77
	/usr/bin/install -c -d /usr/local/cantera/templates/f90
	/usr/bin/install -c -m 644 templates/f90/*.f90 /usr/local/cantera/templates/f90
	/usr/bin/install -c -m 644 templates/f90/*.mak /usr/local/cantera/templates/f90
	/usr/bin/install -c -d /usr/local/cantera/templates/cxx
	/usr/bin/install -c -m 644 templates/cxx/*.cpp /usr/local/cantera/templates/cxx
	/usr/bin/install -c -m 644 templates/cxx/*.mak /usr/local/cantera/templates/cxx
	chown -R groth /usr/local/cantera/templates
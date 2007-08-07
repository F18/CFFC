#/bin/sh
###############################################################
#  $Author: dggoodwin $
#  $Date: 2005/06/18 17:01:00 $
#  $Revision: 1.26 $
#
#  Copyright 2001-2004 California Institute of Technology
#  See file License.txt for licensing information
#
###############################################################

LIBS =  -luser -loneD -lzeroD -ltransport -lcantera -lcvode -lctlapack -lctblas -lctmath -ltpx -lctf2c -lconverters  
SRCS = cantera/private/ctmethods.cpp \
       cantera/private/ctfunctions.cpp \
       cantera/private/xmlmethods.cpp \
       cantera/private/phasemethods.cpp \
       cantera/private/thermomethods.cpp \
       cantera/private/mixturemethods.cpp \
       cantera/private/kineticsmethods.cpp \
       cantera/private/transportmethods.cpp \
       cantera/private/reactormethods.cpp \
       cantera/private/reactornetmethods.cpp \
       cantera/private/wallmethods.cpp \
       cantera/private/flowdevicemethods.cpp \
       cantera/private/onedimmethods.cpp \
       cantera/private/surfmethods.cpp \
       cantera/private/funcmethods.cpp \
       cantera/private/write.cpp


CANTERA_LIBDIR=/nfs/kris/d1/people/charestm/Programs/CFFC/cantera/build/lib/i686-pc-linux-gnu
os_is_win=0

ifeq ($(os_is_win),0)
LIB_DEPS = $(CANTERA_LIBDIR)/libcantera.a $(CANTERA_LIBDIR)/libzeroD.a \
           $(CANTERA_LIBDIR)/liboneD.a \
           $(CANTERA_LIBDIR)/libtransport.a \
           $(CANTERA_LIBDIR)/libclib.a \
           $(CANTERA_LIBDIR)/libconverters.a 
else
LIB_DEPS = $(CANTERA_LIBDIR)/cantera.lib $(CANTERA_LIBDIR)/zeroD.lib \
           $(CANTERA_LIBDIR)/oneD.lib \
           $(CANTERA_LIBDIR)/transport.lib 
endif

DEPENDS = $(SRCS:.cpp=.d)

.cpp.d:
	g++ -MM $(CXX_INCLUDES) $*.cpp > $*.d

all: cantera/ctmethods.mexglx

cantera/ctmethods.mexglx: $(SRCS) $(LIBDEPS) Makefile
ifeq ($(os_is_win),0)
	/usr/bin/python2 setup_matlab.py /usr/local/cantera/bin /nfs/kris/d1/people/charestm/Programs/CFFC/cantera/build/lib/i686-pc-linux-gnu clib '$(LIBS)'
else
	/usr/bin/python2 setup_winmatlab.py
endif
	(matlab -nojvm -nosplash -r setup)
	rm -f setup.m 

install:
	/usr/bin/install -c -d /usr/local/cantera/matlab/toolbox/cantera/cantera
	/usr/bin/install -c -d /usr/local/cantera/matlab/toolbox/cantera/cantera-demos
	/usr/bin/install -c -d /usr/local/cantera/matlab/toolbox/cantera/cantera/private
	/usr/bin/install -c -d /usr/local/cantera/matlab/toolbox/cantera/cantera/@ThermoPhase/private
	/usr/bin/install -c -d /usr/local/cantera/matlab/toolbox/cantera/cantera/@Mixture/private
	/usr/bin/install -c -d /usr/local/cantera/matlab/toolbox/cantera/cantera/@Kinetics/private
	/usr/bin/install -c -d /usr/local/cantera/matlab/toolbox/cantera/cantera/@Transport/private
	/usr/bin/install -c -d /usr/local/cantera/matlab/toolbox/cantera/cantera/@Interface/private
	/usr/bin/install -c -d /usr/local/cantera/matlab/toolbox/cantera/cantera/@Solution
	/usr/bin/install -c -d /usr/local/cantera/matlab/toolbox/cantera/cantera/@XML_Node/private
	/usr/bin/install -c -d /usr/local/cantera/matlab/toolbox/cantera/cantera/@Reactor/private
	/usr/bin/install -c -d /usr/local/cantera/matlab/toolbox/cantera/cantera/@ReactorNet/private
	/usr/bin/install -c -d /usr/local/cantera/matlab/toolbox/cantera/cantera/@Wall/private
	/usr/bin/install -c -d /usr/local/cantera/matlab/toolbox/cantera/cantera/@FlowDevice/private
	/usr/bin/install -c -d /usr/local/cantera/matlab/toolbox/cantera/cantera/@Func/private
	/usr/bin/install -c -d /usr/local/cantera/matlab/toolbox/cantera/cantera/1D
	/usr/bin/install -c -d /usr/local/cantera/matlab/toolbox/cantera/cantera/1D/@Domain1D/private
	/usr/bin/install -c -d /usr/local/cantera/matlab/toolbox/cantera/cantera/1D/@Stack/private
	cd cantera; /usr/bin/install -c *.m *.mexglx /usr/local/cantera/matlab/toolbox/cantera/cantera
	cd cantera/private; /usr/bin/install -c *.m \
            /usr/local/cantera/matlab/toolbox/cantera/cantera/private
	cd cantera/examples; /usr/bin/install -c *.m /usr/local/cantera/matlab/toolbox/cantera/cantera-demos
	cd cantera/@ThermoPhase; /usr/bin/install -c \
            *.m /usr/local/cantera/matlab/toolbox/cantera/cantera/@ThermoPhase
	cd cantera/@ThermoPhase/private; /usr/bin/install -c *.m \
            /usr/local/cantera/matlab/toolbox/cantera/cantera/@ThermoPhase/private
	cd cantera/@Kinetics; /usr/bin/install -c *.m \
            /usr/local/cantera/matlab/toolbox/cantera/cantera/@Kinetics
	cd cantera/@Kinetics/private; /usr/bin/install -c *.m \
            /usr/local/cantera/matlab/toolbox/cantera/cantera/@Kinetics/private
	cd cantera/@Solution; /usr/bin/install -c *.m \
            /usr/local/cantera/matlab/toolbox/cantera/cantera/@Solution
	cd cantera/@Transport; /usr/bin/install -c *.m \
            /usr/local/cantera/matlab/toolbox/cantera/cantera/@Transport
	cd cantera/@Transport/private; /usr/bin/install -c *.m \
            /usr/local/cantera/matlab/toolbox/cantera/cantera/@Transport/private
	cd cantera/@Interface; /usr/bin/install -c \
            *.m /usr/local/cantera/matlab/toolbox/cantera/cantera/@Interface
	cd cantera/@Interface/private; /usr/bin/install -c *.m \
            /usr/local/cantera/matlab/toolbox/cantera/cantera/@Interface/private
	cd cantera/@Mixture; /usr/bin/install -c \
            *.m /usr/local/cantera/matlab/toolbox/cantera/cantera/@Mixture
	cd cantera/@Mixture/private; /usr/bin/install -c *.m \
            /usr/local/cantera/matlab/toolbox/cantera/cantera/@Mixture/private
	cd cantera/@XML_Node; /usr/bin/install -c *.m \
            /usr/local/cantera/matlab/toolbox/cantera/cantera/@XML_Node
	cd cantera/@Reactor; /usr/bin/install -c *.m \
            /usr/local/cantera/matlab/toolbox/cantera/cantera/@Reactor
	cd cantera/@Reactor/private; /usr/bin/install -c *.m \
            /usr/local/cantera/matlab/toolbox/cantera/cantera/@Reactor/private
	cd cantera/@ReactorNet; /usr/bin/install -c *.m \
            /usr/local/cantera/matlab/toolbox/cantera/cantera/@ReactorNet
	cd cantera/@ReactorNet/private; /usr/bin/install -c *.m \
            /usr/local/cantera/matlab/toolbox/cantera/cantera/@ReactorNet/private
	cd cantera/@Wall; /usr/bin/install -c *.m \
            /usr/local/cantera/matlab/toolbox/cantera/cantera/@Wall
	cd cantera/@Wall/private; /usr/bin/install -c *.m \
            /usr/local/cantera/matlab/toolbox/cantera/cantera/@Wall/private
	cd cantera/@FlowDevice; /usr/bin/install -c *.m \
            /usr/local/cantera/matlab/toolbox/cantera/cantera/@FlowDevice
	cd cantera/@FlowDevice/private; /usr/bin/install -c *.m \
            /usr/local/cantera/matlab/toolbox/cantera/cantera/@FlowDevice/private
	cd cantera/@Func; /usr/bin/install -c *.m \
            /usr/local/cantera/matlab/toolbox/cantera/cantera/@Func
	cd cantera/@Func/private; /usr/bin/install -c *.m \
            /usr/local/cantera/matlab/toolbox/cantera/cantera/@Func/private
	cd cantera/1D; /usr/bin/install -c *.m \
            /usr/local/cantera/matlab/toolbox/cantera/cantera/1D
	cd cantera/1D/@Domain1D; /usr/bin/install -c *.m \
            /usr/local/cantera/matlab/toolbox/cantera/cantera/1D/@Domain1D
	cd cantera/1D/@Domain1D/private; /usr/bin/install -c *.m \
            /usr/local/cantera/matlab/toolbox/cantera/cantera/1D/@Domain1D/private
	cd cantera/1D/@Stack; /usr/bin/install -c *.m \
            /usr/local/cantera/matlab/toolbox/cantera/cantera/1D/@Stack
	cd cantera/1D/@Stack/private; /usr/bin/install -c *.m \
            /usr/local/cantera/matlab/toolbox/cantera/cantera/1D/@Stack/private
	/usr/bin/install -c -d /usr/local/cantera/demos/matlab
	/usr/bin/install -c cantera/examples/*.m /usr/local/cantera/demos/matlab
	/usr/bin/install -c -d /usr/local/cantera/tutorials/matlab
	/usr/bin/install -c cantera/tutorial/*.m /usr/local/cantera/tutorials/matlab
	chown -R charestm /usr/local/cantera/demos/matlab
	chown -R charestm /usr/local/cantera/tutorials/matlab


clean:
	rm -f cantera/private/*.o
	rm -f cantera/ctmethods.mexglx

#depends: $(DEPENDS)
#	cat *.d > .depends
#	$(RM) $(DEPENDS) 

depends:
	echo "-"

TAGS: 
	etags *.h *.cpp

ifeq ($(wildcard .depends), .depends)
include .depends
endif

run-demo:
	(cd /usr/local/cantera; matlab -nojvm -nosplash -r cantera_demos)

test-demo:
	(cd /usr/local/cantera/demos/matlab; matlab -nojvm -nosplash -r test_demos)

# end of file




#
# $Revision: 1.3 $
# $Author: hkmoffa $
# $Date: 2006/10/20 21:56:39 $
#
#
test_python=0
test_ck=1
test_cathermo=

all:
	cd minDiamond;  make all
	cd negATest;  make all

test:
	cd minDiamond;  make -s test
	cd negATest;  make -s test

clean:
	$(RM) *.*~ 
	cd minDiamond;  make clean
	cd negATest;  make clean

depends:
	cd minDiamond;  make depends
	cd negATest;  make depends

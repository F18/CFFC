#
# $Revision: 1.14 $
# $Author: hkmoffa $
# $Date: 2006/10/20 21:44:28 $
#
#
test_python=0
test_ck=1
test_cathermo=

all:
	cd cxx_ex;       make all
	cd surfkin;      make all
	cd fracCoeff;    make all
	cd negATest;     make all
	cd diamondSurf;  make all
	cd ChemEquil_gri_matrix; make all
	cd ChemEquil_gri_pairs; make all
ifeq ($(test_python),1)
	cd min_python;   make all
endif
ifeq ($(test_python),2)
	cd min_python;   make all
endif
ifeq ($(test_ck),1)
	cd ck2cti_test;  make all
endif
ifeq ($(test_cathermo),1)
	cd cathermo;     make all
endif

test:
	cd cxx_ex;       make test
	cd surfkin;      make test
	cd fracCoeff;    make test
	cd negATest;     make test
	cd diamondSurf;  make test
	cd ChemEquil_gri_matrix; make test
	cd ChemEquil_gri_pairs; make test
ifeq ($(test_ck),1)
	cd ck2cti_test;  make test
endif
ifeq ($(test_python),1)
	cd min_python;   make test
endif
ifeq ($(test_python),2)
	cd min_python;   make test
	cd python;       make test
endif
ifeq ($(test_cathermo),1)
	cd cathermo;     make test
endif

clean:
	$(RM) *.*~ 
	cd cxx_ex;       make clean
	cd surfkin;      make clean
	cd fracCoeff;    make clean
	cd negATest;     make clean
	cd diamondSurf;  make clean
	cd ChemEquil_gri_matrix; make clean
	cd ChemEquil_gri_pairs; make clean
	cd ck2cti_test;  make clean
	cd min_python;   make clean
	cd python;       make clean
	cd cathermo;     make clean

depends:
	cd cxx_ex;       make depends
	cd surfkin;      make depends
	cd fracCoeff;    make depends
	cd negATest;     make depends
	cd diamondSurf;  make depends
	cd ChemEquil_gri_matrix; make depends
	cd ChemEquil_gri_pairs; make depends
ifeq ($(test_ck),1)
	cd ck2cti_test;  make depends
endif
ifeq ($(test_python),1)
	cd min_python;   make depends
endif
ifeq ($(test_python),2)
	cd min_python;   make depends
	cd python;       make depends
endif
ifeq ($(test_cathermo),1)
	cd cathermo;     make depends
endif

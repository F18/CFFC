#
# $Revision: 1.13 $
# $Author: hkmoffa $
# $Date: 2006/07/12 15:42:19 $
#
#
test_issp=0
test_electrolytes=0

all:
ifeq ($(test_issp),1)
	cd issp;  make all
endif
ifeq ($(test_electrolytes),1)
	cd ims; make all
	cd testIAPWS; make all
	cd testIAPWSPres; make all
	cd testIAPWSTripP; make all
	cd testWaterPDSS; make all
	cd testWaterTP; make all
	cd HMW_graph_GvT; make all
	cd HMW_graph_GvI; make all
	cd HMW_graph_HvT; make all
	cd HMW_graph_CpvT; make all
	cd HMW_graph_VvT; make all
	cd HMW_test_1; make all
	cd HMW_test_3; make all
	cd DH_graph_1; make all
endif

test:
ifeq ($(test_issp),1)
	cd issp;  make -s test
endif
ifeq ($(test_electrolytes),1)
	cd ims; make -s test
	cd testIAPWS; make -s test
	cd testIAPWSPres; make -s test
	cd testIAPWSTripP; make -s test
	cd testWaterPDSS; make -s test
	cd testWaterTP; make -s test
	cd HMW_graph_GvT; make -s test
	cd HMW_graph_GvI; make -s test
	cd HMW_graph_HvT; make -s test
	cd HMW_graph_CpvT; make -s test
	cd HMW_graph_VvT; make -s test
	cd HMW_test_1; make -s test
	cd HMW_test_3; make -s test
	cd DH_graph_1; make -s test
endif

clean:
	$(RM) *.*~ 
	cd issp; make clean
	cd ims;	make clean 
	cd testIAPWS; make clean
	cd testIAPWSPres; make clean
	cd testIAPWSTripP; make clean 
	cd testWaterPDSS; make clean
	cd testWaterTP; make clean
	cd HMW_graph_GvT; make clean
	cd HMW_graph_GvI; make clean
	cd HMW_graph_HvT; make clean
	cd HMW_graph_CpvT; make clean
	cd HMW_graph_VvT; make clean
	cd HMW_test_1; make clean
	cd HMW_test_3; make clean
	cd DH_graph_1; make clean

depends:
ifeq ($(test_issp),1)
	cd issp;make depends
endif
ifeq ($(test_electrolytes),1)
	cd ims; make depends
	cd testIAPWS; make depends
	cd testIAPWSPres; make depends
	cd testIAPWSTripP; make depends
	cd testWaterPDSS; make depends
	cd testWaterTP; make depends
	cd HMW_graph_GvT; make depends
	cd HMW_graph_GvI; make depends
	cd HMW_graph_HvT; make depends
	cd HMW_graph_CpvT; make depends
	cd HMW_graph_VvT; make depends
	cd HMW_test_1; make depends
	cd HMW_test_3; make depends
	cd DH_graph_1; make depends
endif

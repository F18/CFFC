

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#         ********   ***                                                */
#      *******  **  ***       ***      ***               SparseLib++    */
#       *****      ***     ******** ********                            */
#        *****    ***     ******** ********              R. Pozo        */
#   **  *******  ***   **   ***      ***                 K. Remington   */
#    ********   ********                                 A. Lumsdaine   */
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#
#               S p a r s e L i b  ++     m a k e f i l e 
#
#             Edit makefile.def to reflect your environment
#-------------------------------------------------------------------------

#
# Set CFFC source tree
CFFC_SRC_TREE = 2D

################################################################
#
# Subdirectories for building
#
################################################################

error:
	@echo "+-----------------------------------------------------------------+"
	@echo "|                                                                 |"
	@echo "|              SparseLib++ Sparse Matrix Library                  |"	
	@echo "|                                                                 |"
	@echo "| Usage: make all              install and test SparseLib++       |"
	@echo "|        make sp               install Sparselib++  library       |"
	@echo "|        make test             run test suite                     |"
	@echo "|        make clean            clean *.o and test executables     |"
	@echo "|                                                                 |"
	@echo "|  Make sure the system-specific makefile.def has been edited     |"
	@echo "|  to reflect your system configuration.                          |"
	@echo "|                                                                 |"
	@echo "+-----------------------------------------------------------------+"
all: sp test done


sp:
	cd ./src; make CFFC_SRC_TREE=$(CFFC_SRC_TREE);
	cd ./spblas; make CFFC_SRC_TREE=$(CFFC_SRC_TREE);
	cd ./mv/src; make CFFC_SRC_TREE=$(CFFC_SRC_TREE);

test:
	cd ./testing; make CFFC_SRC_TREE=$(CFFC_SRC_TREE);
	cd ./testing; sp_test > sp_test.out;
	@echo "  "
	@echo " +---------------------------------------------------------------+"
	@echo " |                                                               |"
	@echo " |     SparseLib++ test results are in ./testing/sp_test.out     |"
	@echo " |                                                               |"
	@echo " +---------------------------------------------------------------+"
	@echo "  "

clean:
	cd ./src; make CFFC_SRC_TREE=$(CFFC_SRC_TREE) clean;
	cd ./spblas; make CFFC_SRC_TREE=$(CFFC_SRC_TREE) clean;
	cd ./mv/src; make CFFC_SRC_TREE=$(CFFC_SRC_TREE) clean;
	cd ./testing; make CFFC_SRC_TREE=$(CFFC_SRC_TREE) clean;

wipe:
	cd ./src; make CFFC_SRC_TREE=$(CFFC_SRC_TREE) wipe;
	cd ./spblas; make CFFC_SRC_TREE=$(CFFC_SRC_TREE) wipe;
	cd ./mv/src; make CFFC_SRC_TREE=$(CFFC_SRC_TREE) wipe;
	cd ./testing; make CFFC_SRC_TREE=$(CFFC_SRC_TREE) wipe;

done:
	@echo "  "
	@echo " +---------------------------------------------------------------+"
	@echo " |                                                               |"
	@echo " |                   SparseLib++ installed.                      |"
	@echo " |                                                               |"
	@echo " +---------------------------------------------------------------+"
	@echo "  "
	@echo "  "

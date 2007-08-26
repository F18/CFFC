#!/bin/bash

echo "Running Rte2D tests"

#
# set some paths
RTE_EXE="valgrind --leak-check=full $CFFC_Path/src_2D/rte2D"
#RTE_EXE="$CFFC_Path/src_2D/rte2D"

#
# run rectangular coordinate system tests
echo 
echo "Rectangular Enclosure Tests"
echo "==========================="
echo
cd Rectangular_Enclosure
for input in $( ls *.in ); do
    echo
    echo @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    echo  @ $input
    echo @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    echo
    $RTE_EXE -f $input
done
cd ..

#
# run cylindrical coordinate system tests
echo 
echo "Cylindrical Enclosure Tests"
echo "==========================="
echo 
cd Cylindrical_Enclosure
for input in $( ls *.in ); do
    echo
    echo @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    echo @ $input
    echo @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    echo
    $RTE_EXE -f $input
done
cd ..

#
# run cylindrical coordinate system tests
echo 
echo "       SNBCK Tests         "
echo "==========================="
echo 
cd SNBCK
for input in $( ls *.in ); do
    echo
    echo @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    echo @ $input
    echo @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    echo
    $RTE_EXE -f $input
done
cd ..

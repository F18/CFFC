#/bin/sh
#
# $Source: /cvsroot/cantera/cantera/ext/Makefile.in,v $
# $Author: hkmoffa $
# $Revision: 1.8 $
# $Date: 2006/04/30 22:33:54 $
#
# Makefile for ext directory
#
BUILD_LAPACK=1
BUILD_BLAS=1
BUILD_WITH_F2C=1
BUILD_F2C_LIB=1
USE_SUNDIALS=0

LIBS = blas/libctblas.a lapack/libctlapack.a math/libctmath.a \
       cvode/libcvode.a tpx/libtpx.a
 
all: 
ifeq ($(BUILD_F2C_LIB),1)
	(if test -d "f2c_libs" ; then  \
            cd f2c_libs ; make ; \
         fi)
endif
ifeq ($(BUILD_WITH_F2C),1)
ifeq ($(BUILD_LAPACK),1)
	cd f2c_lapack; make
else
	cd f2c_lapack; make clean
endif
ifeq ($(BUILD_BLAS),1)
	cd f2c_blas; make
else
	cd f2c_blas; make clean
endif
#	cd f2c_recipes; make
	cd f2c_math; make
else
ifeq ($(BUILD_LAPACK),1)
	cd lapack; make
else
	cd lapack; make clean
endif
ifeq ($(BUILD_BLAS),1)
	cd blas; make
else
	cd blas; make clean
endif
#	cd recipes; make
	cd math; make
endif
ifeq ($(USE_SUNDIALS),0)
	cd cvode; make
endif
	cd tpx;   make

clean:
	cd lapack; make clean
	cd blas; make clean
	cd recipes; make clean
	cd cvode; make clean
	cd math; make clean
	cd tpx; make clean
	(if test -d "f2c_libs" ; then  \
            cd f2c_libs ; make clean ; \
         fi)
	(if test -d "f2c_lapack" ; then  \
            cd f2c_lapack ; make clean ; \
         fi)
	(if test -d "f2c_blas" ; then  \
            cd f2c_blas ; make clean ; \
         fi)
	(if test -d "f2c_math" ; then  \
            cd f2c_math ; make clean ; \
         fi)
#	(if test -d "f2c_recipes" ; then  \
#            cd f2c_recipes ; make clean ; \
#         fi)

depends:
ifeq ($(BUILD_WITH_F2C),1)
	(if test -d "f2c_libs" ; then  \
            cd f2c_libs ; make depends ; \
         fi)
	(if test -d "f2c_lapack" ; then  \
            cd f2c_lapack ; make depends ; \
         fi)
	(if test -d "f2c_blas" ; then  \
            cd f2c_blas ; make depends ; \
         fi)
	(if test -d "f2c_math" ; then  \
            cd f2c_math ; make depends ; \
         fi)
#	(if test -d "f2c_recipes" ; then  \
            cd f2c_recipes ; make depends ; \
         fi)
else
	cd lapack;   make depends
	cd blas;     make depends
#	cd recipes;  make depends
	cd math;     make depends
endif
	cd cvode;    make depends
	cd tpx;      make depends

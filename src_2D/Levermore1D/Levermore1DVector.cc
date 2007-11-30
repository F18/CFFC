/* Levermore1DVector.cc:  Source file for Levermore1D_Vector class. */

#include "Levermore1DVector.h"

int Levermore1D_Vector::length     = -1;  //default value of -1
                                          //will prohibit creation
                                          //of vectors prior to setting
                                          //the length.
int Levermore1D_Vector::length_set = 0;

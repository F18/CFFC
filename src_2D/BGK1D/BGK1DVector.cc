/* BGK1DVector.cc:  Source file for BGK1D_Vector class. */


#include"BGK1DVector.h"

int BGK1D_Vector::m_length = -1; //default value of -1
                                 //will prohibit creation
                                 //of vectors prior to setting
                                 //the length (which must be 
                                 //constant throughout the simulation).

int BGK1D_Vector::m_length_set = 0;


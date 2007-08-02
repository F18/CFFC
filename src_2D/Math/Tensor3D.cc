/* Tensor3D.cc:  Subroutines for 3D Tensor Class. */

/* Include 3D tensor header file. */

#ifndef _TENSOR3D_INCLUDED
#include "Tensor3D.h"
#endif // _TENSOR3D_INCLUDED

/*************************************************************
 * Tensor3D -- Create storage and initialize temp 3D vector. *
 *************************************************************/
Vector3D Tensor3D::temp_Vec = Vector3D_ZERO;

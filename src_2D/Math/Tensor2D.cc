/* Tensor2D.cc:  Subroutines for 2D Tensor Class. */

/* Include 2D tensor header file. */

#ifndef _TENSOR2D_INCLUDED
#include "Tensor2D.h"
#endif // _TENSOR2D_INCLUDED

/*************************************************************
 * Tensor2D -- Create storage and values for xz and yz.      *
 *************************************************************/
double Tensor2D::xz = ZERO;
double Tensor2D::yz = ZERO;

/*************************************************************
 * Tensor2D -- Create storage and initialize temp 3D vector. *
 *************************************************************/
Vector3D Tensor2D::temp_Vec = Vector3D_ZERO;

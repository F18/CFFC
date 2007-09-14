/* Matrix.cc:  Subroutines for various matrix classes. */

/* Include matrix header file. */

#ifndef _MATRIX_INCLUDED
#include "Matrix.h"
#endif // _MATRIX_INCLUDED

/*************************************************************
 * DenseMatrix -- Create storage for temp vector.            *
 *************************************************************/
RowVector    DenseMatrix::temp_RVec;

/*************************************************************
 * TriDiagonalMatrix -- Create storage for temp vector.      *
 *************************************************************/
RowVector    TriDiagonalMatrix::temp_RVec;


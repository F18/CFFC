/* Reconstruct1D.h:  Header file defining 
                     1D Reconstruction functions and classes. */

#ifndef _RECONSTRUCT_1D_INCLUDED
#define _RECONSTRUCT_1D_INCLUDED

/* Include header files. */
#include <cstring>
#include <cstdlib>
#include <cmath>

#include "Math/Matrix.h"
#include "CompDomain1D.h"

/**************************************************************************
 * Reconstruction 1D -- Solver.                                           *
 **************************************************************************/

double Compute_Geometric_Coeff(CompCell1D_NonUniform **Stencil,
			       int Index_Cell,
			       int Index_Reconstructed_Cell, int Order);

extern void Weighting_LS_Problem (DenseMatrix &A,
				  const RowVector &Delta_XC,
				  ColumnVector &Delta_U,
				  RowVector &W);

extern void Weighting_LS_Problem (DenseMatrix &, RowVector &,
				  const RowVector &, ColumnVector &,
				  const ColumnVector &, double &, int, int &);

extern void Weighting_LS_Problem (DenseMatrix &, RowVector &,
				  const RowVector &, ColumnVector &,
				  const ColumnVector &,
				  const double &, const double &, const double &,
				  int, int &);

extern void Carl_Weighting_LS_Problem (DenseMatrix &, RowVector &,
				  const RowVector &, ColumnVector &,
				  const ColumnVector &,
				  const double &, const double &,
				  int, int &);

extern void AnalyzeWeights (DenseMatrix &, RowVector &, int &);

extern void AnalyzeWeights (DenseMatrix &, RowVector &, const double , int &, const int &);

extern double DetermineMaxDeltaSolutionStencil (CompCell1D_NonUniform **, int , int);

#endif /* _RECONSTRUCT_1D_INCLUDED  */

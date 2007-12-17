/////////////////////////////////////////////////////////////////////
///
/// \file Flame2DdRdU.h
/// 
/// \author Marc R.J. Charest
/// 
/// \brief This header file contains the function prototypes for 
///        Jacobian formulations related to Flame2DState.
///        These Jacobians are used for:
///        -> Semi-Implicit (Source only)
///        -> Newton-Krylov GMRES Preconditioner
///           (1st and 2nd order stencil Inviscid, Viscous, and Source)
///
/////////////////////////////////////////////////////////////////////
#ifndef _FLAME2D_DRDU_INCLUDED
#define _FLAME2D_DRDU_INCLUDED

// inludes
#include "Flame2DQuad.h"

/////////////////////////////////////////////////////////////////////
/// Helper Functions
/////////////////////////////////////////////////////////////////////
void RotationMatrix(DenseMatrix &mat, const Vector2D &nface, const int A_matrix);

/////////////////////////////////////////////////////////////////////
/// Semi-Implicit Block Jacobi formulations
/////////////////////////////////////////////////////////////////////
void SemiImplicitBlockJacobi( DenseMatrix &dSdU,
			      Flame2D_Quad_Block &SolnBlk,
			      Flame2D_pState &Wo,
			      const int &ii, const int &jj);

void SemiImplicitBlockJacobi_dSdU( DenseMatrix &dSdU,
				   Flame2D_Quad_Block &SolnBlk,
				   Flame2D_pState &Wo,
				   const int &ii, const int &jj);
  
void SemiImplicitBlockJacobi_dSdW( DenseMatrix &dSdW,
				   Flame2D_Quad_Block &SolnBlk,
				   Flame2D_pState &Wo,
				   const int &ii, const int &jj);
  
  
//Calculate 2nd derivatives of dW/dx and dW/dy reconstructed using Green Gauss
void d_dWd_dW_Diamond(double &d_dWdx_dW, double &d_dWdy_dW, 
		      Flame2D_Quad_Block &SolnBlk, 
		      const double &LEFT, const double &RIGHT, 
		      const int &Orient_cell, const int &Orient_face,  
		      const int &i, const int &j);
  
//Calculate 2nd derivatives of dW/dx and dW/dy reconstructed using Green Gauss for cell center i,j
void d_dWd_dW_Center(double &d_dWdx_dW_C, double &d_dWdy_dW_C, 
		     Flame2D_Quad_Block &SolnBlk, 
		     const int &i, const int &j);
  
/////////////////////////////////////////////////////////////////////
/// Inviscid Flux Jacobian (FI)
/////////////////////////////////////////////////////////////////////
void dFIdW_Inviscid(DenseMatrix &dRdW, 
		    Flame2D_Quad_Block &SolnBlk, 
		    Flame2D_Input_Parameters &Input_Parameters,
		    const int &ii, const int &jj);

void dFIdW_Inviscid_HLLE(DenseMatrix &dRdW, 
			 Flame2D_Quad_Block &SolnBlk,
			 Flame2D_Input_Parameters &Input_Parameters, 
			 const int &ii, const int &jj, 
			 const int Orient);

void dFIdW_Inviscid_ROE(DenseMatrix& dRdW, Flame2D_Quad_Block &SolnBlk,  
			Flame2D_Input_Parameters &Input_Parameters,
			const int &ii, const int &jj, const int Orient);

void dFIdW_Inviscid_AUSM_plus_up(DenseMatrix& dRdW, 
				 Flame2D_Quad_Block &SolnBlk,  
				 Flame2D_Input_Parameters &Input_Parameters,
				 const int &ii, const int &jj, 
				 const int Orient);

/////////////////////////////////////////////////////////////////////
/// Viscous Flux Jacobian (GV)
/////////////////////////////////////////////////////////////////////
void dWfdWc_Diamond(DenseMatrix &dWfdWc_x,
		    DenseMatrix &dWfdWc_y, 
		    Flame2D_Quad_Block &SolnBlk,
		    const int &Orient_face, 
		    const int &i, const int &j, 
		    const int &Orient_cell);

void dFvdWf_Diamond(DenseMatrix &dFvdWf, 
		    DenseMatrix &dGvdWf, 
		    Flame2D_Quad_Block &SolnBlk,
		    const int &Orient, 
		    const int &ii, const int &jj);


#endif //_FLAME2D_DRDU_INCLUDED

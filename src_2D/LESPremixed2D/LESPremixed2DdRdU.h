#ifndef _LESPREMIXED2D_QUAD_INCLUDED
#include "LESPremixed2DQuad.h"
#endif 

/************************************************************************************
 *  LESPremixed2DdRdU.h and LESPremixed2DdRdU.cc contain Jacobian formulations      *
 *  related to LESPremixed2DState                                                   *
 *                                                                                  *
 *  These Jacobians are used for:                                                   *
 *     -> Semi-Implicit (Source only)                                               *
 *     -> Point-Implicit (Local Inviscid, Viscous, and Source                       *
 *     -> Newton-Krylov GMRES Preconditioner                                        *
 *        (1st and 2nd order stencil Inviscid, Viscous, and Source)                 *
 *                                                                                  * 
 ************************************************************************************/
extern void PointImplicitBlockJacobi( DenseMatrix &dRdU,
				      LESPremixed2D_Quad_Block &SolnBlk,
				      LESPremixed2D_Input_Parameters &Input_Parameters,
				      const int &ii, const int &jj);

extern void SemiImplicitBlockJacobi( DenseMatrix &dSdU,
				     LESPremixed2D_Quad_Block &SolnBlk,
				     const int &solver_type,
				     const int &ii, const int &jj);

extern void SemiImplicitBlockJacobi_dSdU( DenseMatrix &dSdU,
					  LESPremixed2D_Quad_Block &SolnBlk,
					  const int &solver_type,
					  const int &ii, const int &jj);

extern void SemiImplicitBlockJacobi_dSdW( DenseMatrix &dSdW,
					  LESPremixed2D_Quad_Block &SolnBlk,
					  const int &solver_type,
					  const int &ii, const int &jj);

extern DenseMatrix Rotation_Matrix_PM2D(Vector2D nface, int Size, int A_matrix); 

/*****************************************************************************/
//Inviscid Flux Jacobian (FI)
extern void dFIdW_Inviscid(DenseMatrix &dRdW, LESPremixed2D_Quad_Block &SolnBlk, 
			   LESPremixed2D_Input_Parameters &Input_Parameters,
			   const int &ii, const int &jj);

extern void dFIdW_Inviscid_HLLE(DenseMatrix &dRdW, LESPremixed2D_Quad_Block &SolnBlk,
				LESPremixed2D_Input_Parameters &Input_Parameters,
				const int &ii, const int &jj,const int Orient);

extern void dFIdW_Inviscid_ROE(DenseMatrix &dRdW, LESPremixed2D_Quad_Block &SolnBlk, 
			       LESPremixed2D_Input_Parameters &Input_Parameters,
			       const int &ii, const int &jj,const int Orient);

extern void dFIdW_Inviscid_ROE_FD(DenseMatrix &dRdW, LESPremixed2D_Quad_Block &SolnBlk, 
				  LESPremixed2D_Input_Parameters &Input_Parameters,
				  const int &ii, const int &jj,const int Orient);

//Needed by inviscid flux Jacobian -- reconstructed higher order left and right solution states
extern int Inviscid_Flux_Used_Reconstructed_LeftandRight_States(LESPremixed2D_pState &Wl, LESPremixed2D_pState &Wr, Vector2D &dX, 
                                                                LESPremixed2D_Quad_Block &SolnBlk, 
                                                                const int &Orient, const int &ii, const int &jj );

extern void dFIdW_Inviscid_AUSM_plus_up(DenseMatrix &dRdW, LESPremixed2D_Quad_Block &SolnBlk, 
					LESPremixed2D_Input_Parameters &Input_Parameters,
					const int &ii, const int &jj,const int Orient);


/*****************************************************************************/
//Viscous Flux Jacobian (GV)
extern void dGVdW_Viscous(DenseMatrix &dRdW, LESPremixed2D_Quad_Block &SolnBlk,
			  LESPremixed2D_Input_Parameters &Input_Parameters, 
			  const int &ii, const int &jj);

// Viscous Jacobian (including laminar and turbulent) -- Gradients are formulated by following a diamond path
extern  void dFvdWf_Diamond(DenseMatrix &dFvdWf, DenseMatrix &dGvdWf, LESPremixed2D_Quad_Block &SolnBlk, 
			    const int &Orient_face, const int &ii, const int &jj); 

//For diamond path, formulation of Viscous Jacobian needs cell face (Wf) to cell center (Wc) transformation Jacobians
extern  void dWfdWc_Diamond(DenseMatrix &dWfdWc_x,DenseMatrix &dWfdWc_y, LESPremixed2D_Quad_Block &SolnBlk, 
			    const int &Orient_face, const int &ii, const int &jj, const int& Orient_cell);

//Calculate 2nd derivatives of dW/dx and dW/dy reconstructed using Green Gauss
extern void d_dWd_dW_Diamond(double &d_dWdx_dW, double &d_dWdy_dW, LESPremixed2D_Quad_Block &SolnBlk, 
			     const double &LEFT, const double &RIGHT, const int &Orient_cell,
			     const int &Orient_face,  const int &i, const int &j);

//Calculate 2nd derivatives of dW/dx and dW/dy reconstructed using Green Gauss for cell center i,j
extern void d_dWd_dW_Center(double &d_dWdx_dW_C, double &d_dWdy_dW_C, LESPremixed2D_Quad_Block &SolnBlk, 
			    const int &i, const int &j);

/*****************************************************************************/
//Turbulence source Jacobian k-omega
extern int dS_tdW(DenseMatrix &dStdW,  LESPremixed2D_Quad_Block &SolnBlk,
		  double &d_dWdx_dW, double &d_dWdy_dW, const int &ii, const int &jj);

/****************************************************************/
extern int Local_NumVar(LESPremixed2D_Quad_Block &SolnBlk);


//ARE THESE EVER CALLED ???

// extern int Automatic_Wall_Treatment_Residual_Jacobian(LESPremixed2D_Quad_Block &SolnBlk, 
// 						      LESPremixed2D_Input_Parameters &Input_Parameters, 
// 						      int i, int j, 
// 						      DenseMatrix &dRdU);

// extern int Wall_Function_Residual_Jacobian(LESPremixed2D_Quad_Block &SolnBlk, 
// 					   LESPremixed2D_Input_Parameters &Input_Parameters, 
// 					   int i, int j, 
// 					   DenseMatrix &dRdU);

// extern void Low_Reynoldsnumber_Formulation_Residual_Jacobian(LESPremixed2D_Quad_Block &SolnBlk, 
// 							     LESPremixed2D_Input_Parameters &Input_Parameters,
// 							     int i, int j, 
// 							     DenseMatrix &dRdU);

// extern void BC_Residual_Jacobian(LESPremixed2D_Quad_Block &SolnBlk, 
// 				 LESPremixed2D_Input_Parameters &Input_Parameters, 
// 				 int i, int j, 
// 				 DenseMatrix &dRdU);

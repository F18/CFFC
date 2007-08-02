
#ifndef _CHEM2D_QUAD_INCLUDED
#include "Chem2DQuad.h"

#endif 

// The following functions are all used for the construction of Point Implicit Block Jacobi 
// Except the dSa_idU and dSa_VdU are implemented in State Class
// Constructing point implicity point Jacobi
// In blockmatrix P
// the diagonal entries are matrices
//     ---                   ---- 
//     |  ---                    |
//     | | C |                   |
//     |  --- .                  |
//     |         C               |
//     |                         |       
//     |            C            |
//     |                         |     // For laminar case n = 4, adn 6 for turbulent case
// P=  |                         |     C is a matrix that has dimension of (6+Num_Species) *(6+Num_Species). 
//     |                         |     C = C_inv  + C_vis + C_srcaxis_inv +S_srcaxis_vis + S_turbchem 
//     |                 C       |
//     |                         |
//     |                         | 
//     ----                   ----

/*****************************************************************************
 *   Viscous Jacobian Calculations                                           *
 *   Source Jacobians                                                        *                
 *   
 * DenseMatrix dFvdW_Viscous(..)  --- returns x viscous flux Jacobian        *
 * DenseMatrix dGvdW_Viscous(..)  --- returns y viscous flux Jacobian        *
 * DenseMatrix PointImplicitBlkJ_Viscous(...) --- returns C_vis              *
 * DenseMatrix PointImplicitBlkJ_Inviscid(...)--- returns C_inv              *
 * DenseMatrix PointImplicitBlkJ_Src_vis (...)--- returns C_srcaxis_vis      *
 * DenseMatrix PointImplicitBlkJ_Src_turbcomb(...)--- returns C_turbchem     *
 * DenseMatrix PointImplicitBlkJ_Diag(...)--- returns C                      *
 *
 * All the inviscid Jacobians are defined and implemented in Chem2DState.h 
 * and Chem2DState.h                                                         *
 *                                                                           *
 *                                                                           *
 *
 *        JAN. 16, 2003 Xinfeng Gao                                          *
 *                                                                           *
 *        Jun. 12, 2004 Xinfeng Gao                                          *
 *                                                                           *
 *****************************************************************************/
extern  DenseMatrix PointImplicitBlockJacobi(Chem2D_Quad_Block &SolnBlk,
					     Chem2D_Input_Parameters &Input_Parameters,
					     const int &ii, const int &jj);

extern  DenseMatrix SemiImplicitBlockJacobi(Chem2D_Quad_Block &SolnBlk,
					    const int &ii, const int &jj);
//Inviscid Flux Jacobian (FI ---symbole)
extern  DenseMatrix dFIdU_Inviscid(Chem2D_Quad_Block &SolnBlk, 
				   Chem2D_Input_Parameters &Input_Parameters,
				   const int &ii, const int &jj);
extern  DenseMatrix dFIdU_Inviscid_HLLE(Chem2D_Quad_Block &SolnBlk,
					const int &ii, const int &jj);
extern  DenseMatrix dFIdU_Inviscid_ROE(Chem2D_Quad_Block &SolnBlk, 
				       const int &ii, const int &jj);
/*****************************************************************************/
//Viscous Flux Jacobian (GV ---symbole)
extern  DenseMatrix dGVdW_Viscous(Chem2D_Quad_Block &SolnBlk,Chem2D_Input_Parameters &Input_Parameters, 
				  const int &ii, const int &jj);
/*****************************************************************************/
//Laminar Viscous Flux Jacobian  (Fv -- x direction; Gv -- y direction)
extern void dFvdW_Laminar(DenseMatrix &dFvdW,  Chem2D_Quad_Block &SolnBlk,
				 const int &ii, const int &jj); 
extern void dGvdW_Laminar(DenseMatrix &dGvdW,  Chem2D_Quad_Block &SolnBlk,
				 const int &ii, const int &jj); 

//The following two neigbour Jacobians are defined ... because ...
//when computing the viscous flux through cell faces by  averaging quantities of 
//e left (i,j) and its' neigbour cell  
extern void dFvdW_Laminar_Neigbour(DenseMatrix &dFvdW,Chem2D_Quad_Block &SolnBlk,
					  const string &Orient,const int &ii, const int &jj); 
					  
extern void dGvdW_Laminar_Neigbour(DenseMatrix &dGvdW,Chem2D_Quad_Block &SolnBlk,
					  const string &Orient, const int &ii, const int &jj); 

/*****************************************************************************/
extern void dFvdW_Turbulent(DenseMatrix &dFvdW, Chem2D_Quad_Block &SolnBlk,
				   const int &ii, const int &jj); 
extern void dGvdW_Turbulent(DenseMatrix &dGvdW,  Chem2D_Quad_Block &SolnBlk, 
				   const int &ii, const int &jj); 
//The following two neigbour Jacobians are defined ... because ...
//when computing the viscous flux through cell faces by  averaging quantities of the left (i,j) 
//d its' neigbour cell  
extern void dFvdW_Turbulent_Neigbour(DenseMatrix &dFvdW,  Chem2D_Quad_Block &SolnBlk,
					    const string &Orient, const int &ii, const int &jj); 
extern void dGvdW_Turbulent_Neigbour(DenseMatrix &dGvdW,  Chem2D_Quad_Block &SolnBlk,
					    const string &Orient, const int &ii, const int &jj); 

//Turbulence source Jacobian
extern void dStdW_Turbulence(DenseMatrix &dStdW,  Chem2D_Quad_Block &SolnBlk,
				    const int &ii, const int &jj); 
/********************************************************************clear*********/


extern void dFvdW_Laminar_Diamond(DenseMatrix &dFvdW,DenseMatrix &dWf_dWc, Chem2D_Quad_Block &SolnBlk,
					 const string &Orient, const int &ii, const int &jj); 
extern void dGvdW_Laminar_Diamond(DenseMatrix &dGvdW, DenseMatrix &dWf_dWc, Chem2D_Quad_Block &SolnBlk,
					 const string &Orient,	 const int &ii, const int &jj); 
extern void dFvdW_Turbulent_Diamond(DenseMatrix &dFvdW,DenseMatrix &dWf_dWc, Chem2D_Quad_Block &SolnBlk,
					   const string &Orient, const int &ii, const int &jj); 

extern void dGvdW_Turbulent_Diamond(DenseMatrix &dGvdW, DenseMatrix &dWf_dWc, Chem2D_Quad_Block &SolnBlk, 
					   const string &Orient,  const int &ii, const int &jj); 
//For diamond path, formulation of Viscous Jacobian Needs the transformation matrix
// Wf  --- primitive solution variables at cell face
// Wc  --- primitive solution variables at cell center
extern void dWfdWc(Chem2D_Quad_Block &SolnBlk, DenseMatrix dWf_dWc, const string &Orient,  const int &ii, const int &jj); 

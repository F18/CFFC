#ifndef _HIGHTEMP2D_QUAD_NKS_INCLUDED 
#define _HIGHTEMP2D_QUAD_NKS_INCLUDED 

/****************************************************************/
/****** NKS REQUIRED SPECIALIZATIONS & FUNCTIONS ****************/
/****************************************************************/  

#include "HighTemp2DQuad.h"
#include "../NewtonKrylovSchwarz2D/NKS2D.h"
#include "HighTemp2DdRdU.h"

/*! *****************************************************************************************
 *  Specialization of Newton_Update Function                                                *
 *                                                                                          *
 * This routine updates the previous solution data Uo with the deltaU from the GMRES        *                  
 * iterative solver.  U = Uo + GMRES.delatU                                                 *
 *                                                                                          *
 ********************************************************************************************/
template <>
int Newton_Update(HighTemp2D_Quad_Block *SolnBlk,
		  AdaptiveBlock2D_List &List_of_Local_Solution_Blocks,
		  HighTemp2D_Input_Parameters &Input_Parameters,
		  GMRES_RightPrecon_MatrixFree<HighTemp2D_pState,
				HighTemp2D_Quad_Block,
				HighTemp2D_Input_Parameters> &GMRES,
				double update_factor) 
{

	int Num_Var = SolnBlk[0].NumVar();  

	/* Update Solution. No updates to Ghost Cells, let the BC's take care of it */
	for ( int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks.Nblk ; ++Bcount ) {
		if (List_of_Local_Solution_Blocks.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
			for (int j = SolnBlk[Bcount].JCl; j <= SolnBlk[Bcount].JCu; j++){
				for (int i = SolnBlk[Bcount].ICl; i <= SolnBlk[Bcount].ICu; i++){

					// Update solutions in conservative variables:
					// U = Uo + deltaU = Uo + denormalized(x)
					for(int varindex =1; varindex <= Num_Var; varindex++){	
						SolnBlk[Bcount].U[i][j][varindex] = SolnBlk[Bcount].Uo[i][j][varindex] 
							+ update_factor * GMRES.deltaU(Bcount,i,j,varindex-1);
					} 	      

					// Apply update reduction while any one of the updated variables is negative. 
					// This seems to cause more problems than help.
					// Any ideas? Maybe use a global step reduction (broadcast with MPI)?
					if (SolnBlk[Bcount].U[i][j].rho   <= ZERO ||
							SolnBlk[Bcount].U[i][j].e() <= ZERO ||
							SolnBlk[Bcount].U[i][j].E   <= ZERO) {

						double update_reduction_factor = ONE;

						for (int n_update_reduction = 1; n_update_reduction <= 10; ++n_update_reduction) {		  
							update_reduction_factor = HALF*update_reduction_factor;		  		  
							for(int varindex = 1; varindex <= Num_Var; varindex++){		
								SolnBlk[Bcount].U[i][j][varindex] = SolnBlk[Bcount].Uo[i][j][varindex] 
									+ GMRES.deltaU(Bcount,i,j,varindex-1)*update_reduction_factor;
							}   

							cout<<"\n Applying Reduction to solution in NKS "<<n_update_reduction;

							if (SolnBlk[Bcount].U[i][j].rho   > ZERO &&
									SolnBlk[Bcount].U[i][j].E   > ZERO &&
									SolnBlk[Bcount].U[i][j].e() > ZERO ) break; 
						} 
					} 

					if (SolnBlk[Bcount].U[i][j].rho <= ZERO) { 
						cout << "\n NEGATIVE DENSITY : " << endl;
						cout << "SolnBlk["<<Bcount<<"].U["<<i<<"]["<<j<<"].d = " 
							<< SolnBlk[Bcount].U[i][j].rho 
							<< ": SolnBlk["<<Bcount<<"].Uo["<<i<<"]["<<j<<"].d = " 
							<< SolnBlk[Bcount].Uo[i][j].rho << endl;
						cout << "   G["<<Bcount<<"].x[G["<<Bcount<<"].index("<<i
							<<","<<j
							<<")] = " 
							<< GMRES.deltaU(Bcount,i,j,0) << endl;
					}	else if (SolnBlk[Bcount].U[i][j].e() <= ZERO) { 
						cout << "\n NEGATIVE INTERNAL ENERGY : " << endl;
						cout << "SolnBlk["<<Bcount<<"].U["<<i<<"]["<<j<<"].e() = " 
							<< SolnBlk[Bcount].U[i][j].e() << endl;
					}

					// Update solution in primitive variables.
					SolnBlk[Bcount].W[i][j] = W(SolnBlk[Bcount].U[i][j]);	  
				} 
			} 
		} 
	}   

	return 0; 
}

/*!**************************************************************
 * Specialization of Block_Preconditioner::Preconditioner_dFdU  *
 *                                                              *
 * Calculates the dFdU matrix used to generate the approximate  *               
 * Jacobian for the Block Preconditioner.                       *
 ****************************************************************/
template<> 
inline void Block_Preconditioner<HighTemp2D_pState,
			 HighTemp2D_Quad_Block,					    
			 HighTemp2D_Input_Parameters>::
Preconditioner_dFIdU(DenseMatrix &_dFdU, HighTemp2D_pState W)
{
  W.dFdU(_dFdU);
}

/*!**************************************************************
 * Specialization of Block_Preconditioner::Preconditioner_dFIdU *
 *                                                              *
 * Calculates the dFdU matrix used to generate the approximate  *               
 * Jacobian for the Block Preconditioner.                       *
 ****************************************************************/ 
template<> inline void Block_Preconditioner<HighTemp2D_pState,
					    HighTemp2D_Quad_Block,			    
					    HighTemp2D_Input_Parameters>::
Preconditioner_dFIdU_Roe(DenseMatrix &_dFIdU, int ii, int jj, int Orient)
{ 
  DenseMatrix dFI_dW(NUM_VAR_HIGHTEMP2D,NUM_VAR_HIGHTEMP2D,ZERO);
  
  dFIdW_Inviscid_ROE(dFI_dW, *SolnBlk, *Input_Parameters, ii, jj, Orient);  
 
  DenseMatrix dWdU(NUM_VAR_HIGHTEMP2D,NUM_VAR_HIGHTEMP2D,ZERO); 
  //transformation Jacobian  Wo == W here 
  SolnBlk->W[ii][jj].dWdU(dWdU);

  _dFIdU += dFI_dW*dWdU;
}

/*!**************************************************************
 * Specialization of Block_Preconditioner::                     *
 *                               normalize_Preconditioner_dFdU  *
 *                                                              *
 * Normalizes the dFdU matrix used to generate the approximate *               
 * Jacobian for the Block Preconditioner.                       *
 ****************************************************************/
template<> 
inline void Block_Preconditioner<HighTemp2D_pState,
			 HighTemp2D_Quad_Block,
			 HighTemp2D_Input_Parameters>::
normalize_Preconditioner_dFdU(DenseMatrix &dFdU) 
{
	double ao  = HighTemp2D_W_STDATM.a();

	// First row of dFdU is partial of (rho)(u) w.r.t. [U] 
	// and as such should be scaled as:
	dFdU(0,0) *= ONE/ao;
//dFdU(0,1) *= ONE;
//dFdU(0,2) *= ONE;
	dFdU(0,3) *= ao;
	dFdU(0,4) *= ONE/ao;
	dFdU(0,5) *= ONE/ao;

	// Next row is partial of (rho)(u)(u) w.r.t. [U]
	// so it is scaled as the first row, but with an additional
	// division by ao for each entry:
	dFdU(1,0) *= ONE/ao/ao;
	dFdU(1,1) *= ONE/ao;
	dFdU(1,2) *= ONE/ao;
//dFdU(1,3) *= ONE;
	dFdU(1,4) *= ONE/ao/ao;
	dFdU(1,5) *= ONE/ao/ao;

	// Next row is the same.
	dFdU(2,0) *= ONE/ao/ao;
	dFdU(2,1) *= ONE/ao;
	dFdU(2,2) *= ONE/ao;
//dFdU(2,3) *= ONE;
	dFdU(2,4) *= ONE/ao/ao;
	dFdU(2,5) *= ONE/ao/ao;

	// Add an additional /ao for energy flux:
	dFdU(3,0) *= ONE/ao/ao/ao;
	dFdU(3,1) *= ONE/ao/ao;
	dFdU(3,2) *= ONE/ao/ao;
	dFdU(3,3) *= ONE/ao;
	dFdU(3,4) *= ONE/ao/ao/ao;
	dFdU(3,5) *= ONE/ao/ao/ao;

	// We have scaled the two turbulence equations by rho. 
	// So their rows are like the density flux rows:
	dFdU(4,0) *= ONE/ao;
//dFdU(4,1) *= ONE;
//dFdU(4,2) *= ONE;
	dFdU(4,3) *= ao;
	dFdU(4,4) *= ONE/ao;
	dFdU(4,5) *= ONE/ao;

	dFdU(5,0) *= ONE/ao;
//dFdU(5,1) *= ONE;
//dFdU(5,2) *= ONE;
	dFdU(5,3) *= ao;
	dFdU(5,4) *= ONE/ao;
	dFdU(5,5) *= ONE/ao;

}

/************************************************************************/
/************ GMRES REQUIRED SPECIALIZATIONS & FUNCTIONS ****************/
/************************************************************************/   

/*!********************************************************
 * GMRES_Block::set_normalize_values for NKS/GMRES        * 
 *              *                                         *
 *   normalize_values[0] must be set to ao                *
 *   normalize_values[1-n] = values for index[1-n]        *
 *             where n = the number of solution variables *
 **********************************************************/
template<> 
inline void GMRES_Block<HighTemp2D_pState,
			 HighTemp2D_Quad_Block,
			 HighTemp2D_Input_Parameters>::
set_normalize_values(void)
{   
	double ao  = HighTemp2D_W_STDATM.a();
	double rho = HighTemp2D_W_STDATM.rho;

	normalize_valuesU[0] = rho;          //rho
	normalize_valuesU[1] = rho*ao;       //rho*u
	normalize_valuesU[2] = rho*ao;       //rho*v
	normalize_valuesU[3] = rho*ao*ao;    //rho*e
  normalize_valuesU[4] = rho;          //rho*k     
  normalize_valuesU[5] = rho;          //rho*omega 

	normalize_valuesR[0] = rho*ao;
	normalize_valuesR[1] = rho*ao*ao;
	normalize_valuesR[2] = rho*ao*ao;
	normalize_valuesR[3] = rho*ao*ao*ao;
  normalize_valuesR[4] = rho*ao;          //rho*k     
  normalize_valuesR[5] = rho*ao;          //rho*omega 
}

/*!**************************************************************
 * Specialization of Block_Preconditioner::Preconditioner_dFVdU *
 *                                                              *
 * Calculates the dRdU matrix used to generate the approximate  *               
 * Jacobian for the Block Preconditioner.                       *
 * Consult the comments above the function declaration in 
 * NewtonKrylovSchwarz2D/Block_Preconditioner2D.h
 * for an explanation of the input variables to this function
 ****************************************************************/
template<> 
inline void Block_Preconditioner<HighTemp2D_pState,
			 HighTemp2D_Quad_Block,					    
			 HighTemp2D_Input_Parameters>::
			 Preconditioner_dFVdU(DenseMatrix &dFvdU, int Rii, int Rjj, 
					 int Wii, int Wjj, int Orient_face, int Orient_cell) 
{

  double lface;
  Vector2D nface;

  switch(Orient_face){
  case NORTH:
    nface = SolnBlk->Grid.nfaceN(Rii, Rjj);
    lface = SolnBlk->Grid.lfaceN(Rii, Rjj);
    break;
  case SOUTH:
    nface = SolnBlk->Grid.nfaceS(Rii, Rjj);
    lface = SolnBlk->Grid.lfaceS(Rii, Rjj);
    break;
  case EAST:
    nface = SolnBlk->Grid.nfaceE(Rii, Rjj);
    lface = SolnBlk->Grid.lfaceE(Rii, Rjj);
    break;  
  case WEST:
    nface = SolnBlk->Grid.nfaceW(Rii, Rjj);
    lface = SolnBlk->Grid.lfaceW(Rii, Rjj);
    break;  
  }

	int n_vars_wo_vel = NUM_VAR_HIGHTEMP2D - 2;

	// [rho , u , v , p , k , omega , rhox , ux , uy , vx , vy , px , kx , omegax]
  int n_face_vars = (n_vars_wo_vel * 2) + 2 + 2*2;

  DenseMatrix dFvdWf(NUM_VAR_HIGHTEMP2D, n_face_vars, ZERO);
  DenseMatrix dGvdWf(NUM_VAR_HIGHTEMP2D, n_face_vars, ZERO);

  DenseMatrix dWfdWc_x(n_face_vars, NUM_VAR_HIGHTEMP2D, ZERO); 
  DenseMatrix dWfdWc_y(n_face_vars, NUM_VAR_HIGHTEMP2D, ZERO); 

  dFvdWf_Diamond(dFvdWf,   dGvdWf,   *SolnBlk, Orient_face, Rii, Rjj);
  dWfdWc_Diamond(dWfdWc_x, dWfdWc_y, *SolnBlk, Orient_face, Rii, Rjj, Orient_cell); 
  
  DenseMatrix dFvdWc = lface * (nface.x*(dFvdWf*dWfdWc_x) + nface.y*(dGvdWf*dWfdWc_y));
    
  DenseMatrix dWdU(NUM_VAR_HIGHTEMP2D, NUM_VAR_HIGHTEMP2D, ZERO); 

  SolnBlk->W[Wii][Wjj].dWdU(dWdU);

  dFvdU += dFvdWc*dWdU;
}

#endif // _HIGHTEMP2D_QUAD_NKS_INCLUDED 


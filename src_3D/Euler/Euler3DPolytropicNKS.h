#ifndef _EULER3D_POLYTROPIC_NKS_INLCUDED
#define _EULER3D_POLYTROPIC_NKS_INLCUDED

/****************************************************************/
/****** NKS REQUIRED SPECIALIZATIONS & FUNCTIONS ****************/
/****************************************************************/  

#ifndef _NKS_INCLUDED
#include "../NewtonKrylovSchwarz/NKS.h"
#endif // _NKS_INCLUDED

/*! *****************************************************************************************
 *  Specialization of Newton_Update Function                                                *
 *                                                                                          *
 * This routine updates the previous Solution Data Uo with the deltaU from the GMRES        *                  
 * iterative solver.  U = Uo + GMRES.delatU                                                 *
 *                                                                                          *
 ********************************************************************************************/
// template <>
// int Newton_Update(Euler2D_Quad_Block *SolnBlk,
// 		  AdaptiveBlock2D_List &List_of_Local_Solution_Blocks,
// 		  Euler2D_Input_Parameters &Input_Parameters,
// 		  GMRES_RightPrecon_MatrixFree<Euler2D_pState,Euler2D_Quad_Block,Euler2D_Input_Parameters> &GMRES,
// 		  double Relaxation_multiplier) {

//   int Num_Var = SolnBlk[0].NumVar();  
  
//   /* Update Solution. No updates to Ghost Cells, let the BC's take care of it */
//   for ( int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks.Nblk ; ++Bcount ) {
//     if (List_of_Local_Solution_Blocks.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
//       for (int j = SolnBlk[Bcount].JCl; j <= SolnBlk[Bcount].JCu; j++){
// 	for (int i = SolnBlk[Bcount].ICl; i <= SolnBlk[Bcount].ICu; i++){
	  
// 	  /* Update solutions in conversed variables  U = Uo + deltaU = Uo + denormalized(x) */	 
// 	  for(int varindex =1; varindex <= Num_Var; varindex++){	
// 	    SolnBlk[Bcount].U[i][j][varindex] = SolnBlk[Bcount].Uo[i][j][varindex] 
// 	      +  Relaxation_multiplier * GMRES.deltaU(Bcount,i,j,varindex-1);
// 	      } 	      
	  
// 	  /**************************************************************************/
// 	  /* Apply update reduction while any one of the updated variables is negative. */ 
// 	  if (SolnBlk[Bcount].U[i][j].d   <= ZERO ||
// 	      SolnBlk[Bcount].U[i][j].e() <= ZERO ||
// 	      SolnBlk[Bcount].U[i][j].E   <= ZERO) {    //THIS SEEMS TO CAUSE MORE PROBLEMS 
// 	                                                 //THAN HELP, ANY IDEAS, MAYBE MPI ACROSS ALL ????
// 	    double update_reduction_factor = ONE;
	    
// 	    for (int n_update_reduction = 1; n_update_reduction <= 10; ++n_update_reduction) {		  
// 	      update_reduction_factor = HALF*update_reduction_factor;		  		  
// 	      for(int varindex = 1; varindex <= Num_Var; varindex++){		
// 		SolnBlk[Bcount].U[i][j][varindex] = SolnBlk[Bcount].Uo[i][j][varindex] 
// 		  + GMRES.deltaU(Bcount,i,j,varindex-1)*update_reduction_factor;
// 	      }   
	      
// 	      cout<<"\n Applying Reduction to solution in NKS "<<n_update_reduction;
	      
// 	      if (SolnBlk[Bcount].U[i][j].d   > ZERO &&
// 		  SolnBlk[Bcount].U[i][j].E   > ZERO &&
// 		  SolnBlk[Bcount].U[i][j].e() > ZERO ) break; 
// 	    } 
// 	  } 
	  
// 	  /**************************************************************************/
// 	  /* Print error: Negative Density or Negative Pressure. */ 	  
// 	  if (SolnBlk[Bcount].U[i][j].d <= ZERO) { 
// 	    cout << "\n NEGATIVE DENSITY : " << endl;
// 	    cout << "SolnBlk["<<Bcount<<"].U["<<i<<"]["<<j<<"].d = " 
// 		 << SolnBlk[Bcount].U[i][j].d 
// 		 << ": SolnBlk["<<Bcount<<"].Uo["<<i<<"]["<<j<<"].d = " 
// 		 << SolnBlk[Bcount].Uo[i][j].d << endl;
// 	    cout << "   G["<<Bcount<<"].x[G["<<Bcount<<"].index("<<i
// 		 <<","<<j
// 		 <<")] = " 
// 		 << GMRES.deltaU(Bcount,i,j,0) << endl;
// 	  }	else if (SolnBlk[Bcount].U[i][j].e() <= ZERO) { 
// 	    cout << "\n NEGATIVE INTERNAL ENERGY : " << endl;
// 	    cout << "SolnBlk["<<Bcount<<"].U["<<i<<"]["<<j<<"].e() = " 
// 		 << SolnBlk[Bcount].U[i][j].e() << endl;
// 	  } /* endif */ 
	  
// 	  /**************************************************************************/
// 	  //Update solution in primitive variables.
// 	  SolnBlk[Bcount].W[i][j] = W(SolnBlk[Bcount].U[i][j]);	  
// 	} 
//       } 
//     } 
//   }   
//   return 0; 
// }

// /*!**************************************************************
//  * Specialization of Block_Preconditioner::Preconditioner_dFdU  *
//  *                                                              *
//  * Calculates the dFdU matrix used to generate the approximate  *               
//  * Jacobian for the Block Preconditioner.                       *
//  ****************************************************************/
// template<> inline void Block_Preconditioner<Euler2D_pState,
// 					    Euler2D_Quad_Block,					    
// 					    Euler2D_Input_Parameters>::
// Preconditioner_dFIdU(DenseMatrix &_dFdU, Euler2D_pState W)
// {
//   W.dFdU(_dFdU);
// }

/*!**************************************************************
 * Specialization of Block_Preconditioner::                     *
 *                               normalize_Preconditioner_dFdU  *
 *                                                              *
 * Normaliazes the dFdU matrix used to generate the approximate *               
 * Jacobian for the Block Preconditioner.                       *
 ****************************************************************/
template<> inline void Block_Preconditioner<Euler3D_Polytropic_pState, 
					    Euler3D_Polytropic_cState>::
normalize_Preconditioner_dFdU(DenseMatrix &dFdU) 
{
  double ao  = Euler3D_W_STDATM.a();

  dFdU(0,0) *= (ONE/ao);
  dFdU(0,4) *=  ao;

  dFdU(1,0) *= (ONE/sqr(ao));
  dFdU(1,1) *= (ONE/ao);
  dFdU(1,2) *= (ONE/ao);
  dFdU(1,3) *= (ONE/ao);

  dFdU(2,0) *= (ONE/sqr(ao));
  dFdU(2,1) *= (ONE/ao);
  dFdU(2,2) *= (ONE/ao);
  dFdU(2,3) *= (ONE/ao);

  dFdU(3,0) *= (ONE/sqr(ao));
  dFdU(3,1) *= (ONE/ao);
  dFdU(3,2) *= (ONE/ao);
  dFdU(3,3) *= (ONE/ao);

  dFdU(4,0) *= (ONE/cube(ao));
  dFdU(4,1) *= (ONE/sqr(ao));
  dFdU(4,2) *= (ONE/sqr(ao));
  dFdU(4,3) *= (ONE/sqr(ao));  
  dFdU(4,4) *= (ONE/ao); 
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
template<> inline void GMRES_Block<Euler3D_Polytropic_pState, 
				   Euler3D_Polytropic_cState>::
set_normalize_values(void)
{   
  double ao  = Euler3D_W_STDATM.a();
  double rho = Euler3D_W_STDATM.rho;

  normalize_valuesU[0] = rho;          //rho
  normalize_valuesU[1] = rho*ao;       //rho*u
  normalize_valuesU[2] = rho*ao;       //rho*v
  normalize_valuesU[3] = rho*ao;       //rho*w
  normalize_valuesU[4] = rho*ao*ao;    //rho*e

  normalize_valuesR[0] = rho*ao;
  normalize_valuesR[1] = rho*ao*ao;
  normalize_valuesR[2] = rho*ao*ao;
  normalize_valuesR[3] = rho*ao*ao;  
  normalize_valuesR[4] = rho*ao*ao*ao;

}

#endif // _EULER3D_POLYTROPIC_NKS_INLCUDED

#ifndef _EULER3D_THERMALLYPERFECT_NKS_INLCUDED
#define _EULER3D_THERMALLYPERFECT_NKS_INLCUDED

/****************************************************************/
/****** NKS REQUIRED SPECIALIZATIONS & FUNCTIONS ****************/
/****************************************************************/  

#ifndef _NKS_INCLUDED
#include "../NewtonKrylovSchwarz/NKS.h"
#endif // _NKS_INCLUDED


/*! *****************************************************************************************
 *  Euler3D_ThermallyPerfect Specialization of blocksize to use N-1 not N variables         *
 ********************************************************************************************/
template <> 
void Hexa_Newton_Krylov_Schwarz_Solver<Euler3D_ThermallyPerfect_pState, Euler3D_ThermallyPerfect_cState>::
set_blocksize()
{ 
  blocksize = Solution_Data->Local_Solution_Blocks.Soln_Blks[0].NumVar()-1; 
}


/*! *****************************************************************************************
 *  Euler3D_ThermallyPerfect_pState Newton_Update                                           *
 *                                                                                          *
 * This routine updates the previous Solution Data Uo with the deltaU from the GMRES        *                  
 * iterative solver.  U = Uo + GMRES.delatU                                                 *
 *                                                                                          *
 ********************************************************************************************/
template <>
int Hexa_Newton_Krylov_Schwarz_Solver<Euler3D_ThermallyPerfect_pState, Euler3D_ThermallyPerfect_cState>:: 
Newton_Update(){

  /* Update Solution. No updates to Ghost Cells, let the BC's take care of it */
    for ( int Bcount = 0 ; Bcount < Data->Local_Adaptive_Block_List.Nblk; ++Bcount ) {      
      if (Data->Local_Adaptive_Block_List.Block[Bcount].used) {

	for (int k = Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].KCl; 
	     k <= Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].KCu; k++){
	  for (int j = Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].JCl; 
	       j <= Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].JCu; j++){
	    for (int i = Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].ICl; 
		 i <= Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].ICu; i++){

	      //Update solutions in conversed variables  U = Uo + deltaU = Uo + denormalized(x)	 
	      for(int varindex =1; varindex <= Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].NumVar() ; varindex++){  
		Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].U[i][j][k][varindex] = 
		Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].Uo[i][j][k][varindex] 
		+  GMRES.deltaU(Bcount,i,j,k,varindex-1);
	      } 	      	  
 
	      //CHEM2D N-1 
	      Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].U[i][j][k][Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].NumVar()] =
		Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].U[i][j][k].rho*(ONE - 
	        Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].U[i][j][k].sum_species());	   
	
	      //NEED ERROR_CHECKS! FOR -ve species, energy, density 

	      //Update solution in primitive variables.
	      Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].W[i][j][k] = 
		Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].U[i][j][k].W(); 
	    } 
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
template<> inline void Block_Preconditioner<Euler3D_ThermallyPerfect_pState, 
					    Euler3D_ThermallyPerfect_cState>::
Preconditioner_dFIdU(DenseMatrix &_dFdU, Euler3D_ThermallyPerfect_pState W)
{
  W.dFIxdU(_dFdU);
}



/*!**************************************************************
 * Specialization of Block_Preconditioner::                     *
 *                               normalize_Preconditioner_dFdU  *
 *                                                              *
 * Normaliazes the dFdU matrix used to generate the approximate *               
 * Jacobian for the Block Preconditioner.                       *
 ****************************************************************/
template<> inline void Block_Preconditioner<Euler3D_ThermallyPerfect_pState, 
					    Euler3D_ThermallyPerfect_cState>::
normalize_Preconditioner_dFdU(DenseMatrix &dFdU) 
{
  

  Euler3D_ThermallyPerfect_pState W_STD_ATM;
  double ao  = W_STD_ATM.a();
  double rho = W_STD_ATM.rho;

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

  // cs's all have same normalization.
  for(int i= W_STD_ATM.NumVarSansSpecies(); 
      i<  W_STD_ATM.num_vars -1; i++){      //blocksize = W_STD_ATM.num_vars -1
    dFdU(0,i) *= (ONE/ao);
    dFdU(1,i) *= (ONE/(ao*ao));
    dFdU(2,i) *= (ONE/(ao*ao));
    dFdU(3,i) *= (ONE/(ao*ao));
    dFdU(4,i) *= (ONE/(ao*ao*ao));
    dFdU(i,0) *= (ONE/ao);
    dFdU(i,4) *= ao;
    for(int j = W_STD_ATM.NumVarSansSpecies() ; j< W_STD_ATM.num_vars-1; j++){   
      dFdU(i,j) *= (ONE/ao);            
    }
  } 
  
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
template<> inline void GMRES_Block<Euler3D_ThermallyPerfect_pState, 
				   Euler3D_ThermallyPerfect_cState>::
set_normalize_values(void)
{   
  Euler3D_ThermallyPerfect_pState W_STD_ATM;
  double ao  = W_STD_ATM.a();
  double rho = W_STD_ATM.rho;

  // Original Normalization from code  
  normalize_valuesU[0] = rho;          //rho
  normalize_valuesU[1] = rho*ao;       //rho*u
  normalize_valuesU[2] = rho*ao;       //rho*v 
  normalize_valuesU[3] = rho*ao;       //rho*w
  normalize_valuesU[4] = rho*ao*ao;    //rho*e
  for(int i=W_STD_ATM.NumVarSansSpecies();  i < blocksize; i++){
    normalize_valuesU[i] = rho;        //species mass fraction rho*cs
  }

  normalize_valuesR[0] = rho*ao;          //rho
  normalize_valuesR[1] = rho*ao*ao;       //rho*u
  normalize_valuesR[2] = rho*ao*ao;       //rho*v
  normalize_valuesR[3] = rho*ao*ao;       //rho*w
  normalize_valuesR[4] = rho*ao*ao*ao;    //rho*e
  for(int i=W_STD_ATM.NumVarSansSpecies(); i < blocksize; i++){
    normalize_valuesR[i] = rho*ao;        //species mass fraction rho*cs
  }

}

#endif // _EULER3D_THERMALLYPERFECT_NKS_INLCUDED

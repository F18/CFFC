#ifndef _RTE2D_NKS_INCLUDED 
#define _RTE2D_NKS_INCLUDED 

/****************************************************************/
/****** NKS REQUIRED SPECIALIZATIONS & FUNCTIONS ****************/
/****************************************************************/  
#include "Rte2DQuad.h"
#include "../NewtonKrylovSchwarz2D/NKS2D.h"


/*  BLOCK PRECONDITIONER FUNCTION PROTOTYPES */

void dFIdU_FD(DenseMatrix& dRdW, Rte2D_Quad_Block &SolnBlk,  
	      Rte2D_Input_Parameters &Input_Parameters,
	      const int ii, const int jj, const int Orient);

/****************************************************************/
/****** NKS REQUIRED SPECIALIZATIONS & FUNCTIONS ****************/
/****************************************************************/  



/*! *****************************************************************************************
 * Specialization of Newton_Update                                                          *
 *                                                                                          *
 * This routine updates the previous Solution Data Uo with the deltaU from the GMRES        *                  
 * iterative solver.  U = Uo + GMRES.delatU                                                 *
 *                                                                                          *
 ********************************************************************************************/
template <>
inline int Newton_Update(Rte2D_Quad_Block *SolnBlk,
			 AdaptiveBlock2D_List &List_of_Local_Solution_Blocks,
			 Rte2D_Input_Parameters &Input_Parameters,
			 GMRES_RightPrecon_MatrixFree<Rte2D_State,Rte2D_Quad_Block,Rte2D_Input_Parameters> &GMRES,
			 double Relaxation_multiplier) {

  int Num_Var = SolnBlk[0].NumVar();  
  int error_flag = 0;
 
  /* Update Solution. No updates to Ghost Cells, let the BC's take care of it */
  for ( int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks.Nblk ; ++Bcount ) {
    if (List_of_Local_Solution_Blocks.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
      for (int j = SolnBlk[Bcount].JCl; j <= SolnBlk[Bcount].JCu; j++){
	for (int i = SolnBlk[Bcount].ICl; i <= SolnBlk[Bcount].ICu; i++){
	  
	  /* Update solutions in conversed variables  U = Uo + deltaU = Uo + denormalized(x) */	 
	  for(int varindex =1; varindex <= Num_Var; varindex++){  
	    SolnBlk[Bcount].U[i][j][varindex] = SolnBlk[Bcount].Uo[i][j][varindex] 
	      +  Relaxation_multiplier * GMRES.deltaU(Bcount,i,j,varindex-1);
	  } 	      	  
	  // THIS FUNCTION HAS NO CHECKS FOR INVALID SOLUTIONS, 
	  // YOU PROBABLY WANT TO CREATE A SPECIALIZATION OF THIS FUNCTION SPECIFIC 
	  // FOR YOUR EQUATION SYSTEM see Euler2D, Chem2D, etc...
	  /***********************************************************************
	   *************************** RTE SPECIFIC ******************************/
 	  // Apply update reduction while any one of the updated variables is unphysical 
 	  if(SolnBlk[Bcount].U[i][j].Unphysical_Properties()){	   
 	    double update_reduction_factor = ONE;	    
 	    for (int n_update_reduction = 1; n_update_reduction <= 10; ++n_update_reduction) {		  
 	      update_reduction_factor = HALF*update_reduction_factor;		  		  
 	      for(int varindex = 1; varindex <= Num_Var; varindex++){		              
 		SolnBlk[Bcount].U[i][j][varindex] = SolnBlk[Bcount].Uo[i][j][varindex] 
 		  + GMRES.deltaU(Bcount,i,j,varindex-1)*update_reduction_factor;
 	      }   
 	      cout<<"\n Applying Reduction to solution in NKS "<<n_update_reduction;
 	      if( !SolnBlk[Bcount].U[i][j].Unphysical_Properties() )  break;	      
 	    } 
 	  }
	  // Error Check
 	  if(SolnBlk[Bcount].U[i][j].Unphysical_Properties()) error_flag = 1;
	  /*************************************************************************
	   *************************************************************************/
	} 
      } 
    } 
  } 
  
  return (error_flag); 
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
template<> inline void GMRES_Block<Rte2D_State,
				   Rte2D_Quad_Block,
				   Rte2D_Input_Parameters>::
set_normalize_values(void)
{   

  // NOT NORMALIZED
  for(int i=0; i < blocksize; i++) {
    normalize_valuesU[i] = BlackBody(1000.0);//ONE;
    normalize_valuesR[i] = FOUR*PI*BlackBody(1000.0);//ONE;
  }
}



/**************************************************************************
 * Routine: calculate_pertubed_residual                                   *
 **************************************************************************/
// Calculate Soln_ptr.U =  Soln_ptr.Uo + denormalize( epsilon * W(i) )
template <> inline void GMRES_Block<Rte2D_State,
                                    Rte2D_Quad_Block,
                                    Rte2D_Input_Parameters>::
calculate_perturbed_residual(const double &epsilon)
{    
  for (int j = JCl - Nghost ; j <= JCu + Nghost ; j++) {  //includes ghost cells 
    for (int i = ICl - Nghost ; i <= ICu + Nghost ; i++) {
      for(int varindex = 0; varindex < blocksize; varindex++){	
	SolnBlk->U[i][j][varindex+1] = SolnBlk->Uo[i][j][varindex+1] + 
	  denormalizeU( epsilon*W[search_directions*scalar_dim+index(i,j,varindex)], varindex); 
      }   
      /* Update primitive variables. */
      /***********************************************************************
       *************************** RTE SPECIFIC ******************************/
      // removed update primitive variables
      // SolnBlk->W[i][j] = SolnBlk->U[i][j].W();      
      /*************************************************************************
       *************************************************************************/
    }
  }  
}

// Copy forward difference & calculate backwards for 2nd order derivative 
template <> inline void GMRES_Block<Rte2D_State,
                                    Rte2D_Quad_Block,
                                    Rte2D_Input_Parameters>::
calculate_perturbed_residual_2nd(const double &epsilon)
{    
  for (int j = JCl - Nghost ; j <= JCu + Nghost ; j++) {  //includes ghost cells 
    for (int i = ICl - Nghost ; i <= ICu + Nghost ; i++) {
      //copy back R + epsilon * W(i)
      SolnBlk->dUdt[i][j][1] = SolnBlk->dUdt[i][j][0];

      for(int varindex = 0; varindex < blocksize; varindex++){	
	SolnBlk->U[i][j][varindex+1] = SolnBlk->Uo[i][j][varindex+1] - 
	  denormalizeU( epsilon*W[search_directions*scalar_dim+index(i,j,varindex)], varindex); 
      }   
      /* Update primitive variables. */
      /***********************************************************************
       *************************** RTE SPECIFIC ******************************/
      // removed update primitive variables
      // SolnBlk->W[i][j] = SolnBlk->U[i][j].W();      
      /*************************************************************************
       *************************************************************************/
    }
  }  
}

// Calculate Restart Soln_ptr.U =  Soln_ptr.Uo + denormalize( epsilon * x(i) )
template <> inline void GMRES_Block<Rte2D_State,
                                    Rte2D_Quad_Block,
                                    Rte2D_Input_Parameters>::
calculate_perturbed_residual_Restart(const double &epsilon)
{    
  for (int j = JCl - Nghost ; j <= JCu + Nghost ; j++) {
    for (int i = ICl - Nghost ; i <= ICu + Nghost ; i++) {
      for(int varindex = 0; varindex < blocksize; varindex++){	
	SolnBlk->U[i][j][varindex+1] = SolnBlk->Uo[i][j][varindex+1] + 
	  denormalizeU( epsilon*x[index(i,j,varindex)], varindex);
      }  
      /* Update primitive variables. */
      /***********************************************************************
       *************************** RTE SPECIFIC ******************************/
      // removed update primitive variables
      // SolnBlk->W[i][j] = SolnBlk->U[i][j].W();      
      /*************************************************************************
       *************************************************************************/
    }
  }  
}

// Copy forward difference & calculate backwards for 2nd order derivative  
template <> inline void GMRES_Block<Rte2D_State,
                                    Rte2D_Quad_Block,
                                    Rte2D_Input_Parameters>::
calculate_perturbed_residual_2nd_Restart(const double &epsilon)
{    
  for (int j = JCl - Nghost ; j <= JCu + Nghost ; j++) {
    for (int i = ICl - Nghost ; i <= ICu + Nghost ; i++) { 
      //copy back R + epsilon * W(i)
      SolnBlk->dUdt[i][j][1] = SolnBlk->dUdt[i][j][0];
      
      for(int varindex = 0; varindex < blocksize; varindex++){	
	SolnBlk->U[i][j][varindex+1] = SolnBlk->Uo[i][j][varindex+1] -
	  denormalizeU( epsilon*x[index(i,j,varindex)], varindex);
      }  
      /* Update primitive variables. */
      /***********************************************************************
       *************************** RTE SPECIFIC ******************************/
      // removed update primitive variables
      // SolnBlk->W[i][j] = SolnBlk->U[i][j].W();      
      /*************************************************************************
       *************************************************************************/
    }
  }  
}


/**************************************************************************
 * Routine: calculate_epsilon                                             *
 **************************************************************************/
//
// Restart vesion
//
// template <> 
// inline void GMRES_RightPrecon_MatrixFree<Rte2D_State,
// 					 Rte2D_Quad_Block,
// 					 Rte2D_Input_Parameters>::
// calculate_epsilon_restart(double &epsilon)
// {
//   double l2_norm_x(ZERO);
//   double l1_norm_u(ZERO);
//   int total_vars(0);
// 
//   //
//   // loop over each block
//   //
//   for ( int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) {
//     if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
// 
//       // add l2 norm of x component
//       l2_norm_x += sqr( G[Bcount].L2_Norm(G[Bcount].x) );
// 
//       // add l1 norm of U component
//       l1_norm_u += G[Bcount].L1_Norm_Unorm();
// 
//       // compute the total number of variables
//       total_vars += G[Bcount].scalar_dim;
//       
//     } // endif
//   }
// 
//   // compute the global norms
//   l2_norm_x = sqrt(CFFC_Summation_MPI(l2_norm_x));
//   l1_norm_u = CFFC_Summation_MPI(l1_norm_u) / CFFC_Summation_MPI(total_vars);
// 
//   // comput epsilon => (  (1/n)||U||  / ||V||_2 + 1 ) * b
//   epsilon = (l1_norm_u/l2_norm_x + 1.0) * Input_Parameters->NKS_IP.Epsilon_Naught;
// }
// 

//
// Non-Restart version
//
// template <> 
// inline void GMRES_RightPrecon_MatrixFree<Rte2D_State,
// 					 Rte2D_Quad_Block,
// 					 Rte2D_Input_Parameters>::
// calculate_epsilon(double &epsilon, const int &search_direction_counter)
// {
//   double l2_norm_z(ZERO); 
//   double l1_norm_u(ZERO); 
//   int total_vars(0);
// 
//   //
//   // loop over each block
//   //
//   for ( int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks->Nblk ; ++Bcount ) {
//     if ( List_of_Local_Solution_Blocks->Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
// 
//       // add l2 norm of z component
//       l2_norm_z += sqr(G[Bcount].L2_Norm(&(G[Bcount].W[(search_direction_counter)*
// 						       G[Bcount].scalar_dim])));
// 
//       // add l1 norm of U component
//       l1_norm_u += G[Bcount].L1_Norm_Unorm();
// 
//       // compute the total number of variables
//       total_vars += G[Bcount].scalar_dim;
//       
//     } // endif
//   }       
// 
//   // compute the global norms
//   l2_norm_z = sqrt(CFFC_Summation_MPI(l2_norm_z));
//   l1_norm_u = CFFC_Summation_MPI(l1_norm_u) / CFFC_Summation_MPI(total_vars);
// 
//   // comput epsilon => (  (1/n)||U||  / ||V||_2 + 1 ) * b
//   epsilon = (l1_norm_u/l2_norm_z + 1.0) * Input_Parameters->NKS_IP.Epsilon_Naught;
// }
// 

/****************************************************************/
/**** PRECONDTIONER REQUIRED SPECIALIZATIONS & FUNCTIONS ********/
/****************************************************************/  



/*!**************************************************************
 * Specialization of Block_Preconditioner::Preconditioner_dSdU  *
 *                                                              *
 * Calculates the dFdU matrix used to generate the approximate  *               
 * Jacobian for the Block Preconditioner.                       *
 ****************************************************************/
template<> inline void Block_Preconditioner<Rte2D_State,
					    Rte2D_Quad_Block,					    
					    Rte2D_Input_Parameters>::
Preconditioner_dSdU(int ii, int jj, DenseMatrix &dRdU){
  //Source Terms 

 
  //Add Jacobian for axisymmetric source terms
  if (SolnBlk->Axisymmetric) {
    SolnBlk->Uo[ii][jj].dSadU(dRdU,SolnBlk->Sp[ii][jj],
			      SolnBlk->Axisymmetric);
  }  

  //Add Jacobian for regular source terms
  SolnBlk->Uo[ii][jj].dSdU(dRdU, SolnBlk->M[ii][jj]);


}


/*!**************************************************************
 * Specialization of Block_Preconditioner::                     *
 *                               normalize_Preconditioner_dFdU  *
 *                                                              *
 * Normaliazes the dFdU matrix used to generate the approximate *               
 * Jacobian for the Block Preconditioner.                       *
 ****************************************************************/
template<> inline void Block_Preconditioner<Rte2D_State,
					    Rte2D_Quad_Block,
					    Rte2D_Input_Parameters>::
normalize_Preconditioner_dFdU(DenseMatrix &dFdU) { /*EMPTY - NOT NORMALIZED*/ }





/*!**************************************************************
 * Specialization of Block_Preconditioner::                     *
 *                           First_Order_Inviscid_Jacobian_HLLE *
 *                                                              *
 * Calculate First Order Local Jacobian Block(s) Coresponding   *
 * to Cell(i,j) using regular upwind flux function.  We are just*
 * using this function as a quick fix to fit in with the group  *
 * code.                                                        *
 ****************************************************************/
template<> inline void Block_Preconditioner<Rte2D_State,
					    Rte2D_Quad_Block,
					    Rte2D_Input_Parameters>::
First_Order_Inviscid_Jacobian_HLLE(const int &cell_index_i,
				   const int &cell_index_j, 
				   DenseMatrix* Jacobian){              
  
  //! Caculate normal vectors -> in Vector2D format. 
  Vector2D nface_N = SolnBlk->Grid.nfaceN(cell_index_i,cell_index_j-1);
  Vector2D nface_S = SolnBlk->Grid.nfaceS(cell_index_i,cell_index_j+1);
  Vector2D nface_E = SolnBlk->Grid.nfaceE(cell_index_i-1,cell_index_j);
  Vector2D nface_W = SolnBlk->Grid.nfaceW(cell_index_i+1,cell_index_j);

  //! Calculate dFdU using solutions in the rotated frame -> matrix in DenseMatrix format. 
  DenseMatrix dFdU_N(blocksize,blocksize,ZERO); 
  DenseMatrix dFdU_S(blocksize,blocksize,ZERO); 
  DenseMatrix dFdU_E(blocksize,blocksize,ZERO); 
  DenseMatrix dFdU_W(blocksize,blocksize,ZERO); 

  //! Calculate Jacobian matrix -> blocksizexblocksize matrix in DenseMatrix format
  //Solution Rotate provided in pState 
  dFndU( dFdU_N, SolnBlk->Uo[cell_index_i][cell_index_j], nface_N); 
  dFndU( dFdU_S, SolnBlk->Uo[cell_index_i][cell_index_j], nface_S);
  dFndU( dFdU_E, SolnBlk->Uo[cell_index_i][cell_index_j], nface_E);
  dFndU( dFdU_W, SolnBlk->Uo[cell_index_i][cell_index_j], nface_W);

  //! Calculate Jacobian matrix -> blocksizexblocksize matrix in DenseMatrix format
  //North
  Jacobian[NORTH] = SolnBlk->Grid.lfaceN(cell_index_i,cell_index_j-1) 
    * SolnBlk->SpN[cell_index_i][cell_index_j-1] * dFdU_N; 

  //South
  Jacobian[SOUTH] = SolnBlk->Grid.lfaceS(cell_index_i,cell_index_j+1) 
    * SolnBlk->SpS[cell_index_i][cell_index_j+1] * dFdU_S;

  //East
  Jacobian[EAST] = SolnBlk->Grid.lfaceE(cell_index_i-1,cell_index_j) 
    * SolnBlk->SpE[cell_index_i-1][cell_index_j] * dFdU_E;

  //West
  Jacobian[WEST] = SolnBlk->Grid.lfaceW(cell_index_i+1,cell_index_j) 
    * SolnBlk->SpW[cell_index_i+1][cell_index_j] * dFdU_W;

  //Center calculated from neighbours
  //! Using the fact that dF/dU(right) = - dF/dU(left) 
  Jacobian[CENTER] += (Jacobian[NORTH] + Jacobian[SOUTH] + Jacobian[EAST]  + Jacobian[WEST])
                     /(SolnBlk->Grid.Cell[cell_index_i][cell_index_j].A * SolnBlk->Sp[cell_index_i][cell_index_j]);

  Jacobian[NORTH] = -Jacobian[NORTH]/(SolnBlk->Grid.Cell[cell_index_i][cell_index_j-1].A*SolnBlk->Sp[cell_index_i][cell_index_j-1]);
  Jacobian[SOUTH] = -Jacobian[SOUTH]/(SolnBlk->Grid.Cell[cell_index_i][cell_index_j+1].A*SolnBlk->Sp[cell_index_i][cell_index_j+1]);
  Jacobian[EAST] = -Jacobian[EAST]/(SolnBlk->Grid.Cell[cell_index_i-1][cell_index_j].A*SolnBlk->Sp[cell_index_i-1][cell_index_j]);
  Jacobian[WEST] = -Jacobian[WEST]/(SolnBlk->Grid.Cell[cell_index_i+1][cell_index_j].A*SolnBlk->Sp[cell_index_i+1][cell_index_j]);

 
}




/*****************************************************************************
 * Specialization of Block_Preconditioner::                                  *
 *                           First_Order_Inviscid_Jacobian_Roe               *
 *                                                                           *
 *  Calculate First Order Local Jacobian Block(s) Coresponding to Cell(i,j)  *
 *  using Finite difference approximation.  The 'Roe' function is only used  *
 *  to fit in with the rest of the code.                                     *
 *****************************************************************************/
template<> inline void Block_Preconditioner<Rte2D_State,
					    Rte2D_Quad_Block,
					    Rte2D_Input_Parameters>::
First_Order_Inviscid_Jacobian_Roe(const int &cell_index_i,const int &cell_index_j, 
				   DenseMatrix* Jacobian){              
    
  //! Calculate Jacobian matrix -> blocksizexblocksize matrix in DenseMatrix format
  dFIdU_FD(Jacobian[NORTH],*SolnBlk,*Input_Parameters,cell_index_i,cell_index_j,NORTH);
  dFIdU_FD(Jacobian[SOUTH],*SolnBlk,*Input_Parameters,cell_index_i,cell_index_j,SOUTH); 
  dFIdU_FD(Jacobian[EAST],*SolnBlk,*Input_Parameters,cell_index_i,cell_index_j,EAST);        
  dFIdU_FD(Jacobian[WEST],*SolnBlk,*Input_Parameters,cell_index_i,cell_index_j,WEST); 

  //Center calculated from neighbours
  //! Using the fact that dF/dU(right) = - dF/dU(left) 
  Jacobian[CENTER] += (Jacobian[NORTH] + Jacobian[SOUTH] + Jacobian[EAST]  + Jacobian[WEST])
    /(SolnBlk->Grid.Cell[cell_index_i][cell_index_j].A * SolnBlk->Sp[cell_index_i][cell_index_j]);

  Jacobian[NORTH] = -Jacobian[NORTH]/(SolnBlk->Grid.Cell[cell_index_i][cell_index_j-1].A*SolnBlk->Sp[cell_index_i][cell_index_j-1]);
  Jacobian[SOUTH] = -Jacobian[SOUTH]/(SolnBlk->Grid.Cell[cell_index_i][cell_index_j+1].A*SolnBlk->Sp[cell_index_i][cell_index_j+1]);
  Jacobian[EAST] = -Jacobian[EAST]/(SolnBlk->Grid.Cell[cell_index_i-1][cell_index_j].A*SolnBlk->Sp[cell_index_i-1][cell_index_j]);
  Jacobian[WEST] = -Jacobian[WEST]/(SolnBlk->Grid.Cell[cell_index_i+1][cell_index_j].A*SolnBlk->Sp[cell_index_i+1][cell_index_j]);

}


/********************************************************
 * Routine: Flux Jacobian                               *
 *                                                      *
 *     Finite difference approximation                  *
 *                                                      *
 ********************************************************/
inline void dFIdU_FD(DenseMatrix& dRdU, Rte2D_Quad_Block &SolnBlk,  
		     Rte2D_Input_Parameters &Input_Parameters,
		     const int ii, const int jj, const int Orient)
{

  //
  // declares
  //
  int NUM_VAR_RTE2D = dRdU.get_n();   
  DenseMatrix dFidU(NUM_VAR_RTE2D, NUM_VAR_RTE2D,ZERO);
  int Ri, Rj;                 // left face location index
  Vector2D nface;             // face normal
  double lface;               // face length
  Rte2D_State Ul, Ur;         // unperturbed left and right state
  Rte2D_State UA, UB;         // perturbed right states
  Rte2D_State FluxA, FluxB;   // perturbed fluxes 
  Ul.Zero();  Ur.Zero();      // zero states
  double perturb = 5e-6;      // perturbation parameter

  
  //
  // Determine left and right states based on orientation
  //
  switch (Orient) {

    // NORTH
  case NORTH:
    Ri = ii; Rj=jj-1;
    nface = SolnBlk.Grid.nfaceN(Ri, Rj);
    lface = SolnBlk.Grid.lfaceN(Ri, Rj)*SolnBlk.SpN[Ri][Rj];
    Ur = SolnBlk.Uo[ii][jj];
    Ul = SolnBlk.Uo[Ri][Rj];
    break;

    // SOUTH
  case SOUTH:
    Ri = ii; Rj=jj+1;
    nface = SolnBlk.Grid.nfaceS(Ri, Rj);
    lface = SolnBlk.Grid.lfaceS(Ri, Rj)*SolnBlk.SpS[Ri][Rj];
    Ur = SolnBlk.Uo[ii][jj];
    Ul = SolnBlk.Uo[Ri][Rj];
    break;

    // EAST
  case EAST:
    Ri = ii-1; Rj=jj;
    nface = SolnBlk.Grid.nfaceE(Ri, Rj);     
    lface = SolnBlk.Grid.lfaceE(Ri, Rj)*SolnBlk.SpE[Ri][Rj];
    Ur = SolnBlk.Uo[ii][jj];
    Ul = SolnBlk.Uo[Ri][Rj];
    break;

    // WEST
  case WEST:
    Ri = ii+1; Rj=jj;
    nface = SolnBlk.Grid.nfaceW(Ri, Rj);
    lface = SolnBlk.Grid.lfaceW(Ri, Rj)*SolnBlk.SpW[Ri][Rj];
    Ur = SolnBlk.Uo[ii][jj];
    Ul = SolnBlk.Uo[Ri][Rj];
    break;

  } // endswitch 


  //
  // compute the derivatives using finite differences
  //
 for(int jcol=0; jcol<(NUM_VAR_RTE2D); jcol++){
    UA = Ur;
    UB = Ur;
    UA[jcol+1] += perturb; 	 
    UB[jcol+1] -= perturb; 
    
    FluxA = Flux_n(Ul,UA, nface);       
    FluxB = Flux_n(Ul,UB, nface);
    
    for(int irow=0; irow<(NUM_VAR_RTE2D); irow++){
      dFidU(irow,jcol) = (FluxA[irow+1] - FluxB[irow+1])/(TWO*perturb);           
    }
  } 

  // copy over
  dRdU += lface*dFidU;

}


#endif // _RTE2D_NKS_INCLUDED 

#ifndef _ADVECTDIFFUSE2D_IMPLICIT_SPECIALIZATION_INCLUDED 
#define _ADVECTDIFFUSE2D_IMPLICIT_SPECIALIZATION_INCLUDED 

/****************************************************************/
/****** NKS REQUIRED SPECIALIZATIONS & FUNCTIONS ****************/
/****************************************************************/  
#include "AdvectDiffuse2DQuad.h"
#include "../NewtonKrylovSchwarz2D/NKS2D.h"


/*! *****************************************************************************************
 *  Specialization of Newton_Update Function                                                *
 *                                                                                          *
 * This routine updates the previous Solution Data Uo with the deltaU from the GMRES        *                  
 * iterative solver.  U = Uo + GMRES.delatU                                                 *
 *                                                                                          *
 ********************************************************************************************/
template <>
int Newton_Update(AdvectDiffuse2D_Quad_Block *SolnBlk,
		  AdaptiveBlock2D_List &List_of_Local_Solution_Blocks,
		  AdvectDiffuse2D_Input_Parameters &Input_Parameters,
		  GMRES_RightPrecon_MatrixFree<AdvectDiffuse2D_State,
		  AdvectDiffuse2D_Quad_Block,
		  AdvectDiffuse2D_Input_Parameters> &GMRES,
		  double Relaxation_multiplier) {

  int Num_Var = SolnBlk[0].NumVar();  
  
  /* Update Solution. No updates to Ghost Cells, let the BC's take care of it */
  for ( int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks.Nblk ; ++Bcount ) {
    if (List_of_Local_Solution_Blocks.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
      for (int j = SolnBlk[Bcount].JCl; j <= SolnBlk[Bcount].JCu; j++){
	for (int i = SolnBlk[Bcount].ICl; i <= SolnBlk[Bcount].ICu; i++){
	  
	  /* Update solutions in conversed variables  U = Uo + deltaU = Uo + denormalized(x) */	 
	  for(int varindex =1; varindex <= Num_Var; varindex++){	
	    SolnBlk[Bcount].U[i][j][varindex] = ( SolnBlk[Bcount].Uo[i][j][varindex] 
						  +  Relaxation_multiplier * GMRES.deltaU(Bcount,i,j,varindex-1) );
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
template<> inline void Block_Preconditioner<AdvectDiffuse2D_State,
					    AdvectDiffuse2D_Quad_Block,					    
					    AdvectDiffuse2D_Input_Parameters>::
Preconditioner_dFIdU(DenseMatrix &_dFdU, AdvectDiffuse2D_State U)
{
  // Must be revisited
}

/*!**************************************************************
 * Specialization of Block_Preconditioner::                     *
 *                               normalize_Preconditioner_dFdU  *
 *                                                              *
 * Normaliazes the dFdU matrix used to generate the approximate *               
 * Jacobian for the Block Preconditioner.                       *
 ****************************************************************/
template<> inline void Block_Preconditioner<AdvectDiffuse2D_State,
					    AdvectDiffuse2D_Quad_Block,					    
					    AdvectDiffuse2D_Input_Parameters>::
normalize_Preconditioner_dFdU(DenseMatrix &dFdU) 
{
  // Must be revisited
}

/*!**************************************************************
 * Specialization of Block_Preconditioner::                     *
 *                           Update_Jacobian_and_Preconditioner *
 *                                                              *
 *  Update BlockMat with approximation to the Jacobian          *
 * \todo Must be revisited                                      *
 ****************************************************************/
template <>
void Block_Preconditioner<AdvectDiffuse2D_State,
			  AdvectDiffuse2D_Quad_Block,					    
			  AdvectDiffuse2D_Input_Parameters>::
Update_Jacobian_and_Preconditioner(const double &DTS_dTime)
{

  //!Local Variables and Temporary Storage
  int block_mat_size = SolnBlk->NCi*SolnBlk->NCj; 
  DenseMatrix *Jacobian_Data = new DenseMatrix[Jacobian_stencil_size];        //TEMP VAR  
  for(int i=0; i<Jacobian_stencil_size; i++) { Jacobian_Data[i] = DenseMatrix(blocksize,blocksize,ZERO); }
  int *block_i = new int[Jacobian_stencil_size];                                //TEMP VAR
  int *block_j = new int[Jacobian_stencil_size];                                  //TEMP VAR
  
  //! Initially assume no overlap
  int JCl_overlap = 0, JCu_overlap = 0, ICu_overlap =0, ICl_overlap = 0;
  
  //! If overlap determine which block boundaries are internal, ie. BC_NONE
  if(Input_Parameters->NKS_IP.GMRES_Overlap){	
    if (SolnBlk->Grid.BCtypeS[SolnBlk->ICl] == BC_NONE)  JCl_overlap = Input_Parameters->NKS_IP.GMRES_Overlap;
    if (SolnBlk->Grid.BCtypeN[SolnBlk->ICu] == BC_NONE)  JCu_overlap = Input_Parameters->NKS_IP.GMRES_Overlap;
    if (SolnBlk->Grid.BCtypeE[SolnBlk->JCu] == BC_NONE)  ICu_overlap = Input_Parameters->NKS_IP.GMRES_Overlap;
    if (SolnBlk->Grid.BCtypeW[SolnBlk->JCl] == BC_NONE)  ICl_overlap = Input_Parameters->NKS_IP.GMRES_Overlap;
  }

  //*********************************************************************************//
  /*! Calculate Jacobians for each cell and Update Global Jacobian Block Matrix
   * loop through all non-ghost cells, including overlap cells.  Ghost cells already
   * set to Zero or Identity by initialization.
   **********************************************************************************/
  for(int i= SolnBlk->ICl - ICl_overlap; i<= SolnBlk->ICu + ICu_overlap; i++){    
    for(int j= SolnBlk->JCl - JCl_overlap; j<= SolnBlk->JCu + ICu_overlap; j++){  
         
      //--------------------------------------------------------------------------//
      //! Calculate Local Approximate Jacobian                        
      switch(Input_Parameters->NKS_IP.Jacobian_Order){
      case SOURCE_TERMS_ONLY :
 	Implicit_Euler(i,j, Jacobian_Data,DTS_dTime);
	Preconditioner_dSdU(i,j,Jacobian_Data[CENTER]);
	break;
      }

      //--------------------------------------------------------------------------//
      //! Get Block Matrix locations that have components from a given Cell(i,j)
      Get_Block_Index(i,j, block_i, block_j);
      
      //fudge for iGhost Cell reset to zero   //jGhost cell already zero                
      if(block_i[NORTH] < TWO*SolnBlk->NCi) Jacobian_Data[NORTH].zero();
      if(block_i[SOUTH] > block_mat_size - TWO*SolnBlk->NCi) Jacobian_Data[SOUTH].zero();

      if(Jacobian_stencil_size == 9){		
	if(block_i[NORTH_EAST] < TWO*SolnBlk->NCi)    Jacobian_Data[NORTH_EAST].zero();
	if(block_i[NORTH_WEST] < TWO*SolnBlk->NCi)    Jacobian_Data[NORTH_WEST].zero(); 
	if(block_i[SOUTH_EAST] > block_mat_size - TWO*SolnBlk->NCi)    Jacobian_Data[SOUTH_EAST].zero();
	if(block_i[SOUTH_WEST] > block_mat_size - TWO*SolnBlk->NCi)    Jacobian_Data[SOUTH_WEST].zero();
      }
      //--------------------------------------------------------------------------//
      //! Update BlockMat with Local Approximate Jacobians 
      for( int block = 0; block < Jacobian_stencil_size; block++){
	// Normalize
	normalize_Preconditioner_dFdU(Jacobian_Data[block]);

	//can be sped up by more intelligent logic in bkpkit (BlockMat.cc  "setblock")
	Block_Jacobian_approx.setblock( block_i[block], block_j[block], DenseMatrix_to_DenseMat(Jacobian_Data[block]));

	Jacobian_Data[block].zero(); //Just in case to avoid +=/-= issues
      }     
    }
  }

  //Local Memory cleanup
  delete[] Jacobian_Data; delete[] block_i; delete[] block_j;

  //Setup appropriate Preconditioner after Jacobian has been formed/Updated
  Setup_Preconditioner();
}

/************************************************************************/
/************ GMRES REQUIRED SPECIALIZATIONS & FUNCTIONS ****************/
/************************************************************************/   

/*!********************************************************
 * Specialization of GMRES_Block::                        *
 *                   set_normalize_values for NKS/GMRES   * 
 *                                                        *
 *   normalize_values[0] must be set to ao                *
 *   normalize_values[1-n] = values for index[1-n]        *
 *             where n = the number of solution variables *
 **********************************************************/
template<> inline void GMRES_Block<AdvectDiffuse2D_State,
				   AdvectDiffuse2D_Quad_Block,					    
				   AdvectDiffuse2D_Input_Parameters>::
set_normalize_values(void)
{   
  // Must be revisited
}

/*!******************************************************************************
 * Specialization of GMRES_Block::calculate_perturbed_residual_Restart          * 
 *                                                                              *
 * Calculate Restart Soln_ptr.U =  Soln_ptr.Uo + denormalize( epsilon * x(i) )  *
 ********************************************************************************/
template <> inline 
void GMRES_Block<AdvectDiffuse2D_State,
		 AdvectDiffuse2D_Quad_Block,					    
		 AdvectDiffuse2D_Input_Parameters>::
calculate_perturbed_residual(const double &epsilon, const int &order) {    

  int i,j, varindex;
  for ( j = JCl - Nghost ; j <= JCu + Nghost ; j++) {
    for ( i = ICl - Nghost ; i <= ICu + Nghost ; i++) {
      if(order == SECOND_ORDER){
	//store R(Uo + perturb) 
	SolnBlk->dUdt[i][j][1] = SolnBlk->dUdt[i][j][0];
      }
     
      for( varindex = 0; varindex < blocksize; varindex++){
	SolnBlk->U[i][j][varindex+1] = perturbed_resiudal(epsilon,order,i,j,varindex);	
      }  

      //NO W for AdvectDiffuse2D
    }
  }  
}


#endif // _ADVECTDIFFUSE2D_IMPLICIT_SPECIALIZATION_INCLUDED

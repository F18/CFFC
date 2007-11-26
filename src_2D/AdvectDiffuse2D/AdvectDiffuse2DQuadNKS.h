#ifndef _ADVECTDIFFUSE2D_IMPLICIT_SPECIALIZATION_INCLUDED 
#define _ADVECTDIFFUSE2D_IMPLICIT_SPECIALIZATION_INCLUDED 

/****************************************************************/
/****** NKS REQUIRED SPECIALIZATIONS & FUNCTIONS ****************/
/****************************************************************/  
#include "New_AdvectDiffuse2DQuad.h"
#include "../NewtonKrylovSchwarz2D/NKS2D.h"


/*! *****************************************************************************************
 *  Specialization of Newton_Update Function                                                *
 *                                                                                          *
 * This routine updates the previous Solution Data Uo with the deltaU from the GMRES        *                  
 * iterative solver.  U = Uo + GMRES.delatU                                                 *
 *                                                                                          *
 ********************************************************************************************/
template <>
int Newton_Update(AdvectDiffuse2D_Quad_Block_New *SolnBlk,
		  AdaptiveBlock2D_List &List_of_Local_Solution_Blocks,
		  AdvectDiffuse2D_Input_Parameters &Input_Parameters,
		  GMRES_RightPrecon_MatrixFree<AdvectDiffuse2D_State_New,
		  AdvectDiffuse2D_Quad_Block_New,
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
template<> inline void Block_Preconditioner<AdvectDiffuse2D_State_New,
					    AdvectDiffuse2D_Quad_Block_New,					    
					    AdvectDiffuse2D_Input_Parameters>::
Preconditioner_dFIdU(DenseMatrix &_dFdU, AdvectDiffuse2D_State_New U)
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
template<> inline void Block_Preconditioner<AdvectDiffuse2D_State_New,
					    AdvectDiffuse2D_Quad_Block_New,					    
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
void Block_Preconditioner<AdvectDiffuse2D_State_New,
			  AdvectDiffuse2D_Quad_Block_New,					    
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
template<> inline void GMRES_Block<AdvectDiffuse2D_State_New,
				   AdvectDiffuse2D_Quad_Block_New,					    
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
void GMRES_Block<AdvectDiffuse2D_State_New,
		 AdvectDiffuse2D_Quad_Block_New,					    
		 AdvectDiffuse2D_Input_Parameters>::
calculate_perturbed_residual_Restart(const double &epsilon) {    
  int i,j, varindex;

  for ( j = JCl - Nghost ; j <= JCu + Nghost ; j++) {
    for ( i = ICl - Nghost ; i <= ICu + Nghost ; i++) {
      for( varindex = 0; varindex < blocksize; varindex++){	
	SolnBlk->U[i][j][varindex+1] = SolnBlk->Uo[i][j][varindex+1] + 
	  denormalizeU( epsilon*x[index(i,j,varindex)], varindex);
      }  
    }
  }  
}

/*!******************************************************************************
 * Specialization of GMRES_Block::calculate_perturbed_residual_2nd_Restart      * 
 *                                                                              *
 * Copy forward difference & calculate backwards for 2nd order derivative       *
 ********************************************************************************/
template <> inline 
void GMRES_Block<AdvectDiffuse2D_State_New,
		 AdvectDiffuse2D_Quad_Block_New,					    
		 AdvectDiffuse2D_Input_Parameters>::
calculate_perturbed_residual_2nd_Restart(const double &epsilon) {    
  for (int j = JCl - Nghost ; j <= JCu + Nghost ; j++) {
    for (int i = ICl - Nghost ; i <= ICu + Nghost ; i++) { 
      //copy back R + epsilon * W(i)
      SolnBlk->dUdt[i][j][1] = SolnBlk->dUdt[i][j][0];
      
      for(int varindex = 0; varindex < blocksize; varindex++){	
	SolnBlk->U[i][j][varindex+1] = SolnBlk->Uo[i][j][varindex+1] -
	  denormalizeU( epsilon*x[index(i,j,varindex)], varindex);
      }  
    }
  }  
}

/*!******************************************************************************
 * Specialization of GMRES_Block::calculate_perturbed_residual                  * 
 *                                                                              *
 * Calculate Soln_ptr.U =  Soln_ptr.Uo + denormalize( epsilon * W(i) )          *
 ********************************************************************************/
template <> inline 
void GMRES_Block<AdvectDiffuse2D_State_New,
		 AdvectDiffuse2D_Quad_Block_New,					    
		 AdvectDiffuse2D_Input_Parameters>::
calculate_perturbed_residual(const double &epsilon)
{    
  for (int j = JCl - Nghost ; j <= JCu + Nghost ; j++) {  //includes ghost cells 
    for (int i = ICl - Nghost ; i <= ICu + Nghost ; i++) {
      for(int varindex = 0; varindex < blocksize; varindex++){	
	SolnBlk->U[i][j][varindex+1] = SolnBlk->Uo[i][j][varindex+1] + 
	  denormalizeU( epsilon*W[search_directions*scalar_dim+index(i,j,varindex)], varindex); 
      }   
    }
  }  
}

/*!******************************************************************************
 * Specialization of GMRES_Block::calculate_perturbed_residual_2nd              *
 *                                                                              *
 * Copy forward difference & calculate backwards for 2nd order derivative       *
 ********************************************************************************/
template <> inline 
void GMRES_Block<AdvectDiffuse2D_State_New,
		 AdvectDiffuse2D_Quad_Block_New,					    
		 AdvectDiffuse2D_Input_Parameters>::
calculate_perturbed_residual_2nd(const double &epsilon) {    
  for (int j = JCl - Nghost ; j <= JCu + Nghost ; j++) {  //includes ghost cells 
    for (int i = ICl - Nghost ; i <= ICu + Nghost ; i++) {
      //copy back R + epsilon * W(i)
      SolnBlk->dUdt[i][j][1] = SolnBlk->dUdt[i][j][0];

      for(int varindex = 0; varindex < blocksize; varindex++){	
	SolnBlk->U[i][j][varindex+1] = SolnBlk->Uo[i][j][varindex+1] - 
	  denormalizeU( epsilon*W[search_directions*scalar_dim+index(i,j,varindex)], varindex); 
      }   
    }
  }  
}

/*******************************************************************************
 * Specialization of GMRES_Block::LoadSendBuffer_C2F                           *
 *                                                                             *
 *  Loads send message buffer for  coarse to fine block message passing.       *
 *                                                                             *
 * \todo Must be revisited                                                     *
 *******************************************************************************/
template <> inline 
int GMRES_Block<AdvectDiffuse2D_State_New,
		AdvectDiffuse2D_Quad_Block_New,					    
		AdvectDiffuse2D_Input_Parameters>::
LoadSendBuffer_C2F(double *buffer,
		   int &buffer_count,
		   const int buffer_size,
		   const int i_min, 
		   const int i_max,
		   const int i_inc,
		   const int j_min, 
		   const int j_max,
		   const int j_inc,
		   const int face,
		   const int sector) {

   return(1);
}

/**************************************************************************
 * Specialization of GMRES_Block::SubcellReconstruction --                *
 *              Performs the subcell reconstruction of solution state     *
 *              within a given cell (i,j) of the computational mesh for   *
 *              the specified quadrilateral solution block.               *
 **************************************************************************/
template <>
void GMRES_Block<AdvectDiffuse2D_State_New,
		 AdvectDiffuse2D_Quad_Block_New,					    
		 AdvectDiffuse2D_Input_Parameters>::
SubcellReconstruction(const int i, 
 		      const int j,
 		      const int Limiter) {
  
  int n, n2, n_pts, i_index[8], j_index[8], k;
  double u0, u0Min, u0Max, uQuad[4], phi_n;
  double DxDx_ave, DxDy_ave, DyDy_ave;
  Vector2D dX;
  AdvectDiffuse2D_State_New U0, DU, DUDx_ave, DUDy_ave, W_VACUUM;
  W_VACUUM.Vacuum();
  
  /* Carry out the limited solution reconstruction in
     each cell of the computational mesh. */

  //FOR VISCOUS -> CHANGED TO USE ALL 8
  if (i == ICl-Nghost || i == ICu+Nghost ||
      j == JCl-Nghost || j == JCu+Nghost) {
    n_pts = 0;
  } else {
    n_pts = 8;
    i_index[0] = i-1; j_index[0] = j-1;
    i_index[1] = i  ; j_index[1] = j-1;
    i_index[2] = i+1; j_index[2] = j-1;
    i_index[3] = i-1; j_index[3] = j  ;
    i_index[4] = i+1; j_index[4] = j  ;
    i_index[5] = i-1; j_index[5] = j+1;
    i_index[6] = i  ; j_index[6] = j+1;
    i_index[7] = i+1; j_index[7] = j+1;
  } /* endif */

  if (n_pts > 0) {
      DUDx_ave = W_VACUUM;
      DUDy_ave = W_VACUUM;
      DxDx_ave = ZERO;
      DxDy_ave = ZERO;
      DyDy_ave = ZERO;

      for ( k = 0 ; k < blocksize; ++ k) {
	if (vector_switch) {
	  U0[k+1] = W[(search_directions)*scalar_dim + index(i,j,k)];
	} else {
	  U0[k+1] = x[index(i,j,k)];
	}
      } 

      for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
	dX = SolnBlk->Grid.Cell[ i_index[n2] ][ j_index[n2] ].Xc - SolnBlk->Grid.Cell[i][j].Xc;
	for ( k = 0 ; k < blocksize; ++ k) {
	  if (vector_switch) {
	    DU[k+1] = W[(search_directions)*scalar_dim + index(i_index[n2] , j_index[n2] , k)] -  U0[k+1];
	  } else {
	    DU[k+1] = x[index( i_index[n2] , j_index[n2] , k)] -  U0[k+1];
	  } /* endif */
	} /* endfor */
	
	DUDx_ave += DU*dX.x;
	DUDy_ave += DU*dX.y;
	DxDx_ave += dX.x*dX.x;
	DxDy_ave += dX.x*dX.y;
	DyDy_ave += dX.y*dX.y;
      } /* endfor */
  					    
      DUDx_ave = DUDx_ave/double(n_pts);
      DUDy_ave = DUDy_ave/double(n_pts);
      DxDx_ave = DxDx_ave/double(n_pts);
      DxDy_ave = DxDy_ave/double(n_pts);
      DyDy_ave = DyDy_ave/double(n_pts);

      SolnBlk->dUdx[i][j] = (DUDx_ave*DyDy_ave-DUDy_ave*DxDy_ave)/
                            (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
      SolnBlk->dUdy[i][j] = (DUDy_ave*DxDx_ave-DUDx_ave*DxDy_ave)/
                            (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
  
      // Calculate slope limiters. 
      if (!SolnBlk->Freeze_Limiter) {
         for ( n = 1 ; n <= blocksize ; ++n ) {
	    u0 = U0[n];
            u0Min = U0[n];
            u0Max = u0Min;
            for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
               if (vector_switch) {
                  u0Min = min(u0Min, W[(search_directions)*scalar_dim + index(i_index[n2] , j_index[n2] ,n-1)]);
                  u0Max = max(u0Max, W[(search_directions)*scalar_dim + index(i_index[n2] , j_index[n2] ,n-1)]);
               } else {
                  u0Min = min(u0Min, x[index(i_index[n2] , j_index[n2] , n-1)]);
                  u0Max = max(u0Max, x[index(i_index[n2] , j_index[n2] , n-1)]);
               } /* endif */
            } /* endfor */
    
            dX = SolnBlk->Grid.xfaceE(i, j)-SolnBlk->Grid.Cell[i][j].Xc;
            uQuad[0] = u0 + 
                       SolnBlk->dUdx[i][j][n]*dX.x +
                       SolnBlk->dUdy[i][j][n]*dX.y ;
            dX = SolnBlk->Grid.xfaceW(i, j)-SolnBlk->Grid.Cell[i][j].Xc;
            uQuad[1] = u0 + 
                       SolnBlk->dUdx[i][j][n]*dX.x +
                       SolnBlk->dUdy[i][j][n]*dX.y ;
            dX = SolnBlk->Grid.xfaceN(i, j)-SolnBlk->Grid.Cell[i][j].Xc;
            uQuad[2] = u0 + 
                       SolnBlk->dUdx[i][j][n]*dX.x +
                       SolnBlk->dUdy[i][j][n]*dX.y ;
            dX = SolnBlk->Grid.xfaceS(i, j)-SolnBlk->Grid.Cell[i][j].Xc;
            uQuad[3] = u0 + 
                       SolnBlk->dUdx[i][j][n]*dX.x +
                       SolnBlk->dUdy[i][j][n]*dX.y ;
    
            switch(Limiter) {
              case LIMITER_ONE :
                phi_n = ONE;
                break;
              case LIMITER_ZERO :
                phi_n = ZERO;
                break;
              case LIMITER_BARTH_JESPERSEN :
                phi_n = Limiter_BarthJespersen(uQuad, u0, 
                                               u0Min, u0Max, 4);
                break;
              case LIMITER_VENKATAKRISHNAN :
                phi_n = Limiter_Venkatakrishnan(uQuad, u0, 
                                                u0Min, u0Max, 4);
                break;
              case LIMITER_VANLEER :
                phi_n = Limiter_VanLeer(uQuad, u0, 
                                        u0Min, u0Max, 4);
                break;
              case LIMITER_VANALBADA :
                phi_n = Limiter_VanAlbada(uQuad, u0, 
                                          u0Min, u0Max, 4);
                break;
              default:
                phi_n = Limiter_BarthJespersen(uQuad, u0, 
                                               u0Min, u0Max, 4);
                break;
            } /* endswitch */

	    SolnBlk->phi[i][j][n] = phi_n;

         } /* endfor */
      } /* endif */
  } else {
      SolnBlk->dUdx[i][j] = W_VACUUM;
      SolnBlk->dUdy[i][j] = W_VACUUM; 
      SolnBlk->phi[i][j]  = W_VACUUM;
  } /* endif */

}

#endif // _ADVECTDIFFUSE2D_IMPLICIT_SPECIALIZATION_INCLUDED

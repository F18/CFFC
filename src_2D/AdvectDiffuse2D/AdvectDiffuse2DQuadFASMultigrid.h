#ifndef _ADVECTDIFFUSE2D_FASMULTIGRID_SPECIALIZATION_INCLUDED 
#define _ADVECTDIFFUSE2D_FASMULTIGRID_SPECIALIZATION_INCLUDED 

#include "AdvectDiffuse2DQuad.h"
#include "../FASMultigrid2D/FASMultigrid2D.h"

/*! *****************************************************************************************
 *  Specialization of Update_Primitive_Variables                                            *
 *                                                                                          *
 ********************************************************************************************/
template<> inline
void FAS_Multigrid2D_Solver<AdvectDiffuse2D_State,
			    AdvectDiffuse2D_Quad_Block,
			    AdvectDiffuse2D_Input_Parameters>::
Update_Primitive_Variables(const int &Level) {
  // do nothing (there are no conserved and primitive variables in advection diffusion)
}

/*! *********************************************************************
 *  Specialization of Restrict_Boundary_Ref_States (for Multigrid)      *
 *                                                                      *
 *  This routine restricts the boundary refeference states (UoN/S/E/W)  *
 *  from Level_Fine to Level_Coarse for all solution blocks.  The       *
 *  values on the coarse grid level are overwritten by this routine.    *
 *  The restriction operator used is area weighted average.             *
 *                                                                      *
 ***********************************************************************/
template <>
void FAS_Multigrid2D_Solver<AdvectDiffuse2D_State,
			    AdvectDiffuse2D_Quad_Block,
			    AdvectDiffuse2D_Input_Parameters>::
Restrict_Boundary_Ref_States(const int &Level_Fine) {

  int i_fine, j_fine, ICl, ICu, JCl, JCu, Nghost;
  int Level_Coarse = Level_Fine + 1;
  int nb, i_coarse;
  
  // Loop through each local solution block.
  for (nb = 0; nb < List_of_Local_Solution_Blocks[Level_Fine].Nblk; nb++) {
    if (List_of_Local_Solution_Blocks[Level_Fine].Block[nb].used == ADAPTIVEBLOCK2D_USED) { 

      // Get the number of ghost cells.
      Nghost = Local_SolnBlk[Level_Fine][nb].Nghost;

      // Store the lower and upper j-indices on the fine grid.
      JCl = Local_SolnBlk[Level_Fine][nb].JCl;
      JCu = Local_SolnBlk[Level_Fine][nb].JCu;

      // Loop through the i-direction cells of the coarse grid.
      for (i_coarse = Local_SolnBlk[Level_Coarse][nb].ICl; i_coarse <= Local_SolnBlk[Level_Coarse][nb].ICu; i_coarse++) {
	// Determine the i-index of the corresponding SW corner
	// fine cell.
	i_fine = 2*(i_coarse-Nghost)+Nghost;
	// Restrict the north boundary reference state.
	Local_SolnBlk[Level_Coarse][nb].UoN[i_coarse] = 
	  (Local_SolnBlk[Level_Fine][nb].UoN[i_fine]*
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][JCu].A +
	   Local_SolnBlk[Level_Fine][nb].UoN[i_fine+1]*
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][JCu].A)/
	  (Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][JCu].A +
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][JCu].A);
	// Restrict the south boundary reference state.
	Local_SolnBlk[Level_Coarse][nb].UoS[i_coarse] = 
	  (Local_SolnBlk[Level_Fine][nb].UoS[i_fine]*
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][JCl].A +
	   Local_SolnBlk[Level_Fine][nb].UoS[i_fine+1]*
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][JCl].A)/
	  (Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][JCl].A +
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][JCl].A);
      }

      // Store the lower and upper i-indices on the fine grid.
      ICl = Local_SolnBlk[Level_Fine][nb].ICl;
      ICu = Local_SolnBlk[Level_Fine][nb].ICu;

      // Loop through the j-direction cells of the coarse grid.
      for (int j_coarse = Local_SolnBlk[Level_Coarse][nb].JCl; j_coarse <= Local_SolnBlk[Level_Coarse][nb].JCu; j_coarse++) {
	// Determine the j-index of the corresponding SW corner
	// fine cell.
	j_fine = 2*(j_coarse-Nghost)+Nghost;
	// Restrict the west boundary reference state.
	Local_SolnBlk[Level_Coarse][nb].UoW[j_coarse] = 
	  (Local_SolnBlk[Level_Fine][nb].UoW[j_fine]*
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[ICl][j_fine].A +
	   Local_SolnBlk[Level_Fine][nb].UoW[j_fine+1]*
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[ICl][j_fine+1].A)/
	  (Local_SolnBlk[Level_Fine][nb].Grid.Cell[ICl][j_fine].A +
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[ICl][j_fine+1].A);
	// Restrict the east boundary reference state.
	Local_SolnBlk[Level_Coarse][nb].UoE[j_coarse] = 
	  (Local_SolnBlk[Level_Fine][nb].UoE[j_fine]*
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[ICu][j_fine].A +
	   Local_SolnBlk[Level_Fine][nb].UoE[j_fine+1]*
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[ICu][j_fine+1].A)/
	  (Local_SolnBlk[Level_Fine][nb].Grid.Cell[ICu][j_fine].A +
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[ICu][j_fine+1].A);
      }

    }
  }

}


/*! *********************************************************************
 * Specialization: Additional_Solution_Block_Setup                      *
 *                                                                      *
 * Perform some additional setup on the local solution block before     *
 * beginning Multi Grid computions.                                     *
 *                                                                      *
 ***********************************************************************/
template <> 
inline void FAS_Multigrid2D_Solver<AdvectDiffuse2D_State,
				   AdvectDiffuse2D_Quad_Block,
				   AdvectDiffuse2D_Input_Parameters>::
Additional_Solution_Block_Setup(AdvectDiffuse2D_Quad_Block &SolnBlk,
				AdvectDiffuse2D_Quad_Block &SolnBlk_FinestLevel)
{

  // High-order related variables
  int i;
  vector<int> ReconstructionOrders;

  // High-order variables and their reconstruction order
  for (i = 0; i < SolnBlk_FinestLevel.NumberOfHighOrderObjects(); ++i){
    ReconstructionOrders.push_back(SolnBlk_FinestLevel.HighOrderVariable(i).RecOrder());
  }

  // allocate memory for high-order variables based on the setup of the FinestLevel
  SolnBlk.allocate_HighOrder(SolnBlk_FinestLevel.NumberOfHighOrderObjects(),
			     ReconstructionOrders);
  // allocate memory for high-order boundary conditions if necessary
  SolnBlk.allocate_HighOrder_BoundaryConditions();

}



#endif	// _ADVECTDIFFUSE2D_FASMULTIGRID_SPECIALIZATION_INCLUDED 

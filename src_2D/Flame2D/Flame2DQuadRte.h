/**************************************************************** 
   Flame2DQuadRte.h

   This header file defines the specializations required to 
   use the Radiation solver Rte2D with Flame2D.  The solution 
   of the radiation transfer equation is necessary to compute
   the radiation heat flux term in the energy equation.  It is
   solved sequentially.
     
*****************************************************************/
#ifndef _FLAME2D_RTE_INCLUDED 
#define _FLAME2D_RTE_INCLUDED 

/********************************************************
 * Required includes                                    *
 ********************************************************/
#include "Flame2DQuad.h"
#include "../Rte2D/Rte2DQuadSolvers.h"

/****************************************************************/
/****** RTE REQUIRED SPECIALIZATIONS & FUNCTIONS ****************/
/****************************************************************/  

/********************************************************
 * Copy_SRC_Solution_Vars                               *
 *                                                      *
 * Copy over the necessary Flame2D solution variables to *
 * the Rte2D solution block.  This is necessary for     *
 * solution of the radiation transfer equation.  Some   *
 * processing is required to compute the absorbsion     *
 * coefficients based on the Flame2D state species       *
 * concentration and the blackbody emmissive power.     *
 *                                                      *
 * NOTE: This function is a member of 'Rte2DSolver',    *
 * so all objects such as 'Soln_Block_List' and         *
 * 'SolnBlk' refer to the RTE2D solution state. Those   *
 * marked with a prefix 'Flame2D_' refer to Flame2D's.    *
 *                                                      *
 * NOTE: The meshes and quadtree should be identical,   *
 * so the block and cell index should also match.       *
 *                                                      *
 ********************************************************/
template <>
inline void Rte2DSolver::Copy_SRC_Solution_Vars(Flame2D_Quad_Block *Flame2D_SolnBlk) 
{ 
  // defines
  double xCO, xH2O, xCO2, xO2;
  static const double fsoot = ZERO;  // no soot for now
  
  // make some pointers to help readability
  Rte2D_Quad_Block *Rte2D_SolnBlk = Local_SolnBlk;
  AdaptiveBlock2D_List *Rte2D_Soln_Block_List = &List_of_Local_Solution_Blocks;

  //
  // Loop over each used block
  //
  for ( int n=0 ; n<Rte2D_Soln_Block_List->Nblk; n++ ) {
    if (Rte2D_Soln_Block_List->Block[n].used == ADAPTIVEBLOCK2D_USED) {

      //
      // loop over all the cells in the grid
      //
      for ( int i=Rte2D_SolnBlk[n].ICl-Rte2D_SolnBlk[n].Nghost; 
	    i<=Rte2D_SolnBlk[n].ICu+Rte2D_SolnBlk[n].Nghost; 
	    i++ ) {
	for ( int j=Rte2D_SolnBlk[n].JCl-Rte2D_SolnBlk[n].Nghost; 
	      j<=Rte2D_SolnBlk[n].JCu+Rte2D_SolnBlk[n].Nghost; 
	      j++ ) {

	  // first, get radiating species concentrations
	  Flame2D_SolnBlk[n].W[i][j].MoleFracOfRadSpec( xCO,  xH2O, xCO2, xO2 );

	  // compute the new state
 	  Rte2D_SolnBlk[n].M[i][j].ComputeNewState( Flame2D_SolnBlk[n].W[i][j].p,
						    Flame2D_SolnBlk[n].W[i][j].T(),
						    xCO,
						    xH2O,
						    xCO2,
						    xO2,
						    fsoot );

	} // endfor - i
      } // endfor - j

    } // endif
  } // endfor - n

}

/********************************************************
 * Copy_Rte2D_Solution_Vars                             *
 *                                                      *
 * Copy over the computed Rte2D solution variables back *
 * to the Flame2D solution block.  
 *                                                      *
 * NOTE: This function is a member of 'Rte2DSolver',    *
 * so all objects such as 'Soln_Block_List' and         *
 * 'SolnBlk' refer to the RTE2D solution state. Those   *
 * marked with a prefix 'Flame2D_' refer to Flame2D's.    *
 *                                                      *
 * NOTE: The meshes and quadtree should be identical,   *
 * so the block and cell index should also match.       *
 *                                                      *
 ********************************************************/
template <>
inline void Rte2DSolver::Copy_Rte2D_Solution_Vars(Flame2D_Quad_Block *Flame2D_SolnBlk)
{ 

  // make some pointers to help readability
  Rte2D_Quad_Block *Rte2D_SolnBlk = Local_SolnBlk;
  AdaptiveBlock2D_List *Rte2D_Soln_Block_List = &List_of_Local_Solution_Blocks;

  //
  // Loop over each used block
  //
  for ( int n=0 ; n<Rte2D_Soln_Block_List->Nblk; n++ ) {
    if (Rte2D_Soln_Block_List->Block[n].used == ADAPTIVEBLOCK2D_USED) {

      //
      // loop over all the cells in the grid
      //
      for ( int i=Rte2D_SolnBlk[n].ICl-Rte2D_SolnBlk[n].Nghost; 
	    i<=Rte2D_SolnBlk[n].ICu+Rte2D_SolnBlk[n].Nghost; 
	    i++ ) {
	for ( int j=Rte2D_SolnBlk[n].JCl-Rte2D_SolnBlk[n].Nghost; 
	      j<=Rte2D_SolnBlk[n].JCu+Rte2D_SolnBlk[n].Nghost; 
	      j++ ) {

	  // compute the divergence of the radiant heat flux
 	  Flame2D_SolnBlk[n].Srad[i][j] = 
	    Rte2D_SolnBlk[n].U[i][j].Sr( Rte2D_SolnBlk[n].M[i][j] );
						     
	} // endfor - i
      } // endfor - j

    } // endif
  } // endfor - n

};

#endif // _FLAME2D_RTE_INCLUDED 







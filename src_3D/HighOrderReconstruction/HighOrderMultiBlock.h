/*!\file HighOrder_MultiBlock.h
  \brief Header file implementing HighOrder_MultiBlock class. */

#ifndef _HIGHORDER_MULTIBLOCK_INCLUDED
#define _HIGHORDER_MULTIBLOCK_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#ifndef _HIGHORDER_INCLUDED
#include "HighOrder.h"	      // Include high-order class header file
#endif //_HIGHORDER_INCLUDED

#ifndef _HIGHORDER_INPUT_INCLUDED
#include "HighOrderInput.h"   // Include high-order input header file
#endif //_HIGHORDER_INPUT_INCLUDED

/*!
 * \class HighOrder_MultiBlock
 *
 * \brief Collection of templated functions that operate
 *        on high-order variables of 3D arrays of 
 *        quadrilateral solution blocks.
 * \nosubgrouping
 *********************************************************/
class HighOrder_MultiBlock {
public:

  //! Create initial high-order variables
  template<typename HEXA_BLOCK>
  static int Create_Initial_HighOrder_Variables(HEXA_BLOCK *Solution_Block,
						AdaptiveBlock3D_List &LocalSolnBlockList){

    int i;

    //  if (HighOrder_Input::NumberOfHighOrderReconstructions != 0){

    for (i = 0 ; i <= LocalSolnBlockList.Nblk-1 ; ++i ) {
      if (LocalSolnBlockList.Block[i].used == ADAPTIVEBLOCK3D_USED) {
	// allocate high-order objects
	Solution_Block[i].allocate_HighOrder();
	// allocate high-order boundary conditions
	// Solution_Block[i].allocate_HighOrder_BoundaryConditions(); //RR: allocate_HighOrder_Boundary_Conditions()
     
      } /* endif */
   
    }  /* endfor */
 
    //  } //endif (NumberOfHighOrderReconstructions)

    return(0);
  }

//  //! @brief Perform high-order reconstruction for all solution variables using the specified high-order object
//  template <class Quad_Soln_Block, class Input_Parameters_Type>
//  static void HighOrder_Reconstruction(Quad_Soln_Block *Soln_ptr,
//				       AdaptiveBlock3D_List &LocalSolnBlockList,
//				       Input_Parameters_Type & IP,
//				       const unsigned short int Pos,
//				       const typename Quad_Soln_Block::Soln_State &
//				       (Quad_Soln_Block::*ReconstructedSoln)(const int &,const int &) const = 
//				       &Quad_Soln_Block::CellSolution);

private:
  HighOrder_MultiBlock(void);   //!< Private default constructor
  HighOrder_MultiBlock(const HighOrder_MultiBlock&); //!< Private copy constructor
  HighOrder_MultiBlock& operator=(const HighOrder_MultiBlock&); //!< Private assignment operator

};

/*!
 * Create as many high-order block structures as specified
 * in the input parameters and initialize them based on the 
 * specified orders of reconstruction.
 */
//template<typename HEXA_BLOCK>
//static int Create_Initial_HighOrder_Variables(HEXA_BLOCK *Solution_Block,
//				       AdaptiveBlock3D_List &LocalSolnBlockList){
//
//  int i;
//
//  //  if (HighOrder_Input::NumberOfHighOrderReconstructions != 0){
//
//  for (i = 0 ; i <= LocalSolnBlockList.Nblk-1 ; ++i ) {
//    if (LocalSolnBlockList.Block[i].used == ADAPTIVEBLOCK3D_USED) {
//      // allocate high-order objects
//      Solution_Block[i].allocate_HighOrder();
//      // allocate high-order boundary conditions
//      // Solution_Block[i].allocate_HighOrder_BoundaryConditions(); //RR: allocate_HighOrder_Boundary_Conditions()
//     
//    } /* endif */
//   
//  }  /* endfor */
// 
//  //  } //endif (NumberOfHighOrderReconstructions)
//
//  return(0);
//}


/*!
 * Perform high-order reconstruction for 
 * all used solution blocks in the provided list.
 * 
 * \param Pos the index of the particular high-order object
 *            in the high-order variable array of the solution block.
 * \param ReconstructedSoln member function of Quad_Soln_Block which returns the solution.
 */

// --> RR: To be included later (reconstruct final solution)

//template <class Quad_Soln_Block, class Input_Parameters_Type>
//void HighOrder_MultiBlock::HighOrder_Reconstruction(Quad_Soln_Block *Soln_ptr,
//						      AdaptiveBlock3D_List &LocalSolnBlockList,
//						      Input_Parameters_Type & IP,
//						      const unsigned short int Pos,
//						      const typename Quad_Soln_Block::Soln_State &
//						      (Quad_Soln_Block::*ReconstructedSoln)(const int &,
//											    const int &) const){
//
//  for (int i = 0 ; i <= LocalSolnBlockList.Nblk-1 ; ++i ) {
//    if (LocalSolnBlockList.Block[i].used == ADAPTIVEBLOCK3D_USED) {
//      // Reconstruct the solution in the current solution block.
//      Soln_ptr[i].HighOrderVariable(Pos).ComputeHighOrderSolutionReconstruction(Soln_ptr[i],
//										IP.Limiter(),
//										ReconstructedSoln);
//    } /* endif */
//  }  /* endfor */
//
//}

#endif

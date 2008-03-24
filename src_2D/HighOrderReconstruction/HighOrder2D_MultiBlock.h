/*!\file HighOrder2D_MultiBlock.h
  \brief Header file implementing HighOrder2D_MultiBlock class. */

#ifndef _HIGHORDER_2D_MULTIBLOCK_INCLUDED
#define _HIGHORDER_2D_MULTIBLOCK_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "HighOrder2D.h"	// Include 2D high-order class header file
#include "HighOrder2D_Input.h"	// Include 2D high-order input header file

/*!
 * \class HighOrder2D_MultiBlock
 *
 * \brief Collection of templated functions that operate
 *        on high-order variables of 2D arrays of 
 *        quadrilateral solution blocks.
 * \nosubgrouping
 *********************************************************/
class HighOrder2D_MultiBlock {
public:

  //! Create initial high-order variables
  template <class Quad_Soln_Block>
  static void Create_Initial_HighOrder_Variables(Quad_Soln_Block *Soln_ptr,
						 AdaptiveBlock2D_List &LocalSolnBlockList);

private:
  HighOrder2D_MultiBlock(void);   //!< Private default constructor
  HighOrder2D_MultiBlock(const HighOrder2D_MultiBlock&); //!< Private copy constructor
  HighOrder2D_MultiBlock& operator=(const HighOrder2D_MultiBlock&); //!< Private assignment operator

};

/*!
 * Create as many high-order block structures as specified
 * in the input parameters and initialize them based on the 
 * specified orders of reconstruction.
 */
template <class Quad_Soln_Block>
void HighOrder2D_MultiBlock::Create_Initial_HighOrder_Variables(Quad_Soln_Block *Soln_ptr,
								AdaptiveBlock2D_List &LocalSolnBlockList){
  int i;

  if (HighOrder2D_Input::NumberOfHighOrderReconstructions != 0){

    for (i = 0 ; i <= LocalSolnBlockList.Nblk-1 ; ++i ) {
      if (LocalSolnBlockList.Block[i].used == ADAPTIVEBLOCK2D_USED) {
	// add initialization at block level
      } /* endif */
    }  /* endfor */

  } //endif (NumberOfHighOrderReconstructions)
  
}


#endif

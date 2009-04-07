/*!\file HighOrder.h
  \brief Header file defining the templated HighOrder class.*/

#ifndef _HIGHORDER_INCLUDED
#define _HIGHORDER_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
using std::ostream;
using std::istream;
using std::vector;

/* Include CFFC header files */
#include "../Math/Matrix.h"
#include "../Math/NumericalLibrary.h"
#include "TDContainer.h"
#include "CENO_ExecutionMode.h"
#include "CENO_Tolerances.h"
#include "../Grid/Grid3DHexaBlock.h"
//#include "ReconstructionHelpers.h"
//#include "Cauchy_BoundaryConditions.h" // --> RR: replace later
//#include "HighOrder2D_BlockBoundary.h" // --> RR: replace later

/*********************************
 * Declare the HighOrder   class *
 ********************************/
template <class SOLN_STATE> 
class HighOrder;

/************************************************
 *     Friend Functions : HighOrder             *
 ************************************************/
template<class SOLN_STATE>
ostream & operator<< (ostream & os, const HighOrder<SOLN_STATE> & Obj);

template<class SOLN_STATE>
istream & operator>> (istream & os, HighOrder<SOLN_STATE> & Obj);

//!  \class HighOrder
//   -----------------------------------------------------------------
/*!  
 *   \brief Templated class for high-order variables in 3D.
 *
 *///-----------------------------------------------------------------
template<class SOLN_STATE>
class HighOrder{

public:

private:

};


/* ---------------------------------------------------------------------------------------------- 
 * =============== INCLUDE THE IMPLEMENTATION OF HIGH-ORDER RECONSTRUCTIONS ==================
 * ---------------------------------------------------------------------------------------------*/
//#include "HighOrder_Reconstructions.h" --> RR: include later

#endif // _HIGHORDER_INCLUDED

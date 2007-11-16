/*!\file Grid2DQuadIntegration.h
  \brief Header file defining 2D quadrilateral mesh integration class. */

#ifndef _GRID2D_QUAD_INTEGRATION_INCLUDED
#define _GRID2D_QUAD_INTEGRATION_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../Math/NumericalLibrary.h"


template<class Grid2DQuadType>
class Grid2DQuadIntegration{
public:

  //! Constructor with Grid
  Grid2DQuadIntegration(Grid2DQuadType * AssociatedGrid);

  //! Destructor
  ~Grid2DQuadIntegration(void){ };

  //! Compute the integral of a general function over the domain of cell (ii,jj)
  template<typename FO, class ReturnType>
  ReturnType IntegrateFunctionOverCell(const int ii, const int jj, const FO FuncObj,
				       int digits, ReturnType _dummy_param);
  
  //! Compute the integral of a polynomial function over the domain of cell (ii,jj)
  template<typename FO, class ReturnType>
  ReturnType IntegratePolynomialOverCell(const int ii, const int jj, const FO FuncObj,
					 ReturnType _dummy_param);
  
  
private:
  Grid2DQuadType *Grid;		//!< pointer to the grid associated to this object

  Grid2DQuadIntegration(void);	//!< Private default constructor
};


// Constructor with the Grid that is going to be associated with this object
template<class Grid2DQuadType>
Grid2DQuadIntegration<Grid2DQuadType>::Grid2DQuadIntegration(Grid2DQuadType * AssociatedGrid){
  Grid = AssociatedGrid;
}



#endif

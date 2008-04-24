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

  //! Re-associate geometry pointer
  void AssociateGeometry(Grid2DQuadType * AssociatedGrid){ Grid = AssociatedGrid; }

  //! Access to grid
  Grid2DQuadType * getGrid(void) const { return Grid; }

  //! Compute the integral of a general function over the domain of cell (ii,jj)
  template<typename FO, class ReturnType>
  ReturnType IntegrateFunctionOverCell(const int &ii, const int &jj, const FO FuncObj,
				       int digits, ReturnType _dummy_param) const;

  //! Compute the integral of a general function over the potentially curved domain of cell (ii,jj)
  template<typename FO, class ReturnType>
  ReturnType IntegrateFunctionOverCell(const int &ii, const int &jj, const FO FuncObj,
				       const FO ContourIntegrand, int digits,
				       ReturnType _dummy_param) const;
  
  //! Compute the integral of a polynomial function over the domain of cell (ii,jj)
  template<typename FO, class ReturnType>
  ReturnType IntegratePolynomialOverCell(const int &ii, const int &jj, const FO FuncObj,
					 ReturnType _dummy_param) const;

private:
  Grid2DQuadType *Grid;	        //!< pointer to the grid associated to this object

  Grid2DQuadIntegration(void);	//!< Private default constructor  
};


// Constructor with the Grid that is going to be associated with this object
template<class Grid2DQuadType> inline
Grid2DQuadIntegration<Grid2DQuadType>::Grid2DQuadIntegration(Grid2DQuadType * AssociatedGrid){
  Grid = AssociatedGrid;
}

/*!
 * Integrate a general function (i.e. any function or pointer function)
 * over the domain of a cell (ii,jj), the boundaries of     
 * which are considered to be straight
 *
 * \param ii the i-index of the cell over which the integration is performed
 * \param jj the j-index of the cell over which the integration is performed
 * \param FuncObj the function to be integrated
 * \param digits the number of exact digits with which the result is computed (i.e. the accuracy of the calculation)
 * \param _dummy_param a parameter used only to determine the return type of the function FuncObj
 */
template<class Grid2DQuadType>
template<typename FO, class ReturnType> inline
ReturnType Grid2DQuadIntegration<Grid2DQuadType>::IntegrateFunctionOverCell(const int &ii, const int &jj,
									    const FO FuncObj, int digits, 
									    ReturnType _dummy_param) const {
  return QuadrilateralQuadrature(FuncObj,
				 Grid->nodeSW(ii,jj),
				 Grid->nodeNW(ii,jj),
				 Grid->nodeNE(ii,jj),
				 Grid->nodeSE(ii,jj),
				 digits,
				 _dummy_param);
}

/*!
 * Integrate a polynomial function (i.e. any function or pointer function)
 * over the domain of a cell (ii,jj), the boundaries of which are considered to be straight edges.
 * This subroutine is suitable for integrating exactly & efficiently polynomials of up to order 18.
 *
 * \param ii the i-index of the cell over which the integration is performed
 * \param jj the j-index of the cell over which the integration is performed
 * \param FuncObj the function to be integrated
 * \param digits the number of exact digits with which the result is computed (i.e. the accuracy of the calculation)
 * \param _dummy_param a parameter used only to determine the return type of the function FuncObj
 */
template<class Grid2DQuadType>
template<typename FO, class ReturnType> inline
ReturnType Grid2DQuadIntegration<Grid2DQuadType>::IntegratePolynomialOverCell(const int &ii, const int &jj,
									      const FO FuncObj,
									      ReturnType _dummy_param) const {
  return GaussLobattoQuadrilateralQuadrature(FuncObj,
					     Grid->nodeSW(ii,jj),
					     Grid->nodeNW(ii,jj),
					     Grid->nodeNE(ii,jj),
					     Grid->nodeSE(ii,jj),
					     _dummy_param);
}

/*!
 * Integrate a general function (i.e. any function or pointer function)
 * over the domain of a cell (ii,jj). The boundaries of the cell can    
 * be straight or curved.
 * If all boundaries are straight then the function is integrated with
 * the QuadrilateralQuadrature routine.
 * If curved boundaries are detected then Gauss integration along the
 * boundary contour is applied using the ContourIntegrand function.
 *
 * \param ii the i-index of the cell over which the integration is performed
 * \param jj the j-index of the cell over which the integration is performed
 * \param FuncObj the function to be integrated
 * \param ContourIntegrand the integrand with respect to x of FuncObj
 * \param digits the number of exact digits with which the result is computed (i.e. the accuracy of the calculation)
 * \param _dummy_param a parameter used only to determine the return type of the function FuncObj
 */
template<class Grid2DQuadType>
template<typename FO, class ReturnType> inline
ReturnType Grid2DQuadIntegration<Grid2DQuadType>::IntegrateFunctionOverCell(const int &ii, const int &jj, const FO FuncObj,
									    const FO ContourIntegrand, int digits,
									    ReturnType _dummy_param) const{
  
  // == check if high-order boundary treatment is required
  if ( Grid->IsHighOrderBoundary() ){

    // === Decide whether to use contour integration or surface integration

#if 0
    if ( (ii > Grid->ICl) && (ii < Grid->ICu) && (jj > Grid->JCl) && (jj < Grid->JCu) ){
      // This is an interior cell unaffected by curved boundaries
      return IntegrateFunctionOverCell(ii,jj,FuncObj,digits,_dummy_param);

    } else if ( (ii < Grid->ICl) || (ii > Grid->ICu) || (jj < Grid->JCl) || (jj > Grid->JCu) ){
      // This is a ghost cell unaffected by curved boundaries
      return IntegrateFunctionOverCell(ii,jj,FuncObj,digits,_dummy_param);

    } else if ( ( (ii == Grid->ICl  && Grid->BndWestSplineInfo != NULL)  || 
		( (ii == Grid->ICu  && Grid->BndEastSplineInfo != NULL)  || 
		( (jj == Grid->JCl  && Grid->BndSouthSplineInfo != NULL)  || 
		( (jj == Grid->JCu  && Grid->BndNorthSplineInfo != NULL)   ){

      // This cell needs Gauss contour integration
      
    } else {
      // This is a cell 
    }
#endif

  } else {
    // all cells are treated with low-order accuracy so use the integration over quadrilaterals.
    return IntegrateFunctionOverCell(ii,jj,FuncObj,digits,_dummy_param);
  } // endif
}

#endif

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
#include "../Math/RandomGenerator.h"


template<class Grid2DQuadType>
class Grid2DQuadIntegration{
public:

  //! Constructor with Grid
  Grid2DQuadIntegration(Grid2DQuadType * AssociatedGrid);

  //! Destructor
  ~Grid2DQuadIntegration(void){ Deallocate(); }

  //! Re-associate geometry pointer
  void AssociateGeometry(Grid2DQuadType * AssociatedGrid){ Grid = AssociatedGrid; }

  //! Access to grid
  Grid2DQuadType * getGrid(void) const { return Grid; }

  //! Deallocate object memory
  void Deallocate(void);

  //! Access to type of cell faces
  const bool* getCellFacesInfo(void){ return CellFacesInfo; }
  bool getWestFaceInfo(void){ return CellFacesInfo[0]; }
  bool getSouthFaceInfo(void){ return CellFacesInfo[1]; }
  bool getEastFaceInfo(void){ return CellFacesInfo[2]; }
  bool getNorthFaceInfo(void){ return CellFacesInfo[3]; }
  bool getFaceBlockEdgeCorrelation(void) const { if (IsTheSameBlockEdge == 1) { return true;} else { return false;} }

  //! Get the type of this cell (i.e. the type of each edge)
  void AnalyseCellFaces(const int &ii, const int &jj);

  //! Print information about the cell type (this is valid only for interior cells)
  void PrintCellInfo(ostream & os) const;

  //! Compute the integral of a general function over the domain of cell (ii,jj)
  template<typename FO, class ReturnType>
  ReturnType IntegrateFunctionOverCell(const int &ii, const int &jj, const FO FuncObj,
				       int digits, ReturnType _dummy_param) const;

  //! Compute the integral of a general function over the potentially curved domain of cell (ii,jj)
  template<typename FO, class ReturnType>
  ReturnType IntegrateFunctionOverCell(const int &ii, const int &jj, const FO FuncObj,
				       FO ContourIntegrand, int digits,
				       ReturnType _dummy_param);

  //! Compute the integral of a general function over the domain of cell (ii,jj) which has curved faces
  template<typename FO, class ReturnType>
  ReturnType IntegrateFunctionOverCellUsingContourIntegrand(const int &ii, const int &jj,
							    const FO FuncObj, FO ContourIntegrand,
							    int digits, ReturnType _dummy_param) const {
    throw runtime_error("Grid2DQuadIntegration()::IntegrateFunctionOverCellUsingContourIntegrand() not implemented for this Grid type");
  }

  //! Compute the averages given by a function over the sectors of cell (ii,jj) (i.e. SW, NW, SE, NE) using contour integration
  template<typename FO, class ReturnType>
  void IntegrateFunctionOverCellSectorsUsingContourIntegrand(const int &ii, const int &jj,
							     FO ContourIntegrand,
							     const ReturnType & CoarseAvgSoln,
							     ReturnType & AvgSoln_SW, ReturnType & AvgSoln_NW,
							     ReturnType & AvgSoln_SE, ReturnType & AvgSoln_NE) {
    throw runtime_error("Grid2DQuadIntegration()::IntegrateFunctionOverCellSectorsUsingContourIntegrand() not implemented for this Grid type");
  }

  //! Compute the integral of a general function over the domain of cell (ii,jj) which has curved faces using Monte Carlo method
  template<typename FO, class ReturnType>
  ReturnType IntegrateFunctionOverCellUsingMonteCarloMethod(const int &ii, const int &jj,
							    FO FuncObj,
							    int digits, ReturnType _dummy_param) const {
    throw runtime_error("Grid2DQuadIntegration()::IntegrateFunctionOverCellUsingMonteCarloMethod() not implemented for this Grid type");
  }
  
  /*! Compute the integral of a general function over the domain of cell (ii,jj) which
    has curved faces using adaptive quadrilaterals method */
  template<typename FO, class ReturnType>
  ReturnType IntegrateFunctionOverCellUsingPolygonalAdaptiveQuadratures(const int &ii, const int &jj,
									FO FuncObj,
									int digits, ReturnType _dummy_param) const {
    throw runtime_error("Grid2DQuadIntegration()::IntegrateFunctionOverCellUsingPolygonalAdaptiveQuadratures() not implemented for this Grid type");
  }

  //! Approximate the curvilinear domain of cell (ii,jj) with a polygon
  void ConvertCurvedIntegrationDomainToPolygon(const int &ii, const int &jj,
					       Polygon & PolygonalDomain) const {
    throw runtime_error("Grid2DQuadIntegration()::ConvertCurvedIntegrationDomainToPolygon() not implemented for this Grid type");
  }

  //! Compute the integral of an arbitrary function using Gauss quadrature rule with the provided data.
  template<typename FO, class ReturnType>
  ReturnType CalculateFunctionIntegralWithGaussQuadratures(FO FuncObj, 
							   const Vector2D * GaussQuads, const double * GaussWeights,
							   const int & NumGQPs, const double & DeltaY,
							   ReturnType _dummy_param) const;
  
  //! Compute the integral of a polynomial function over the domain of cell (ii,jj)
  template<typename FO, class ReturnType>
  ReturnType IntegratePolynomialOverCell(const int &ii, const int &jj, const FO FuncObj,
					 ReturnType _dummy_param) const;

  //! @brief Compute the integral of a polynomial function over an arbitrary quadrilateral domain
  template<typename FO, typename Node2DType, class ReturnType>
  ReturnType IntegratePolynomialOverQuadDomain(const FO FuncObj, 
					       const Node2DType &SW, const Node2DType &NW,
					       const Node2DType &NE, const Node2DType &SE,
					       ReturnType _dummy_param) const;

  //! @brief Compute the integral of the product between a scalar function and the normal vector
  template<typename FO, class ReturnType>
  void IntegrateFunctionProjectionAlongBoundarySpline(const int & BOUNDARY, const FO FuncObj,
						      ReturnType & ResultXdir, ReturnType & ResultYdir,
						      double & WettedSurface) const;

  //! @brief Similar to IntegrateFunctionProjectionAlongBoundarySpline() but for piecewise function (i.e. different for each cell)
  template<typename FO, class ReturnType>
  void IntegratePiecewiseFunctionProjectionAlongBoundarySpline(const int & BOUNDARY, FO FuncObj,
							       ReturnType & ResultXdir, ReturnType & ResultYdir,
							       double & WettedSurface) const;

  //! @brief Compute the integral of the wall shear stress along a boundary (i.e. using information from each cell)
  template<typename FO, typename FO_TestValidityDomain>
  void IntegratePiecewiseWallShearStressAlongBoundarySpline(const int & BOUNDARY, FO FuncObj,
							    double & ResultXdir, double & ResultYdir,
							    double & WettedSurface,
							    FO_TestValidityDomain & ValidateDomain) const;

private:
  Grid2DQuadType *Grid;	        //!< pointer to the grid associated to this object

  Grid2DQuadIntegration(void);	//!< Private default constructor  

  bool *CellFacesInfo;	        //!< Array to track the type of each cell face
  int IsTheSameBlockEdge;       /*!< Variable to mark the association between the curved cell face and the block side. 
				  (e.g a West cell face corresponds to a West or East block side)
				  It basically makes the distinction between interior cells and ghost cells. 
				  This variable takes value +1 for interior cells and -1 for ghost cells. */
  bool AtLeastOneCurvedFace;	//!< Flag to indicate whether the cell does have at least one curved face

  //! Define pointers to cell curved boundaries.
  typename Grid2DQuadType::BndSplineType *CellBndWest, *CellBndSouth, *CellBndEast, *CellBndNorth;


  /*!
   * \class IntegrandFunctionOverPolygonalDefinitionDomain
   *
   * Functor class which enforces an arbitrary function to 
   * be defined only on the inside domain of an arbitrary polygon.
   * Outside of the polygon, the function value is zero.
   */
  template<class FunctionType, class ReturnType>
  class IntegrandFunctionOverPolygonalDefinitionDomain {
  public:
    //! Main constructor
    IntegrandFunctionOverPolygonalDefinitionDomain(FunctionType & Func, 
						   Polygon &_DefinitionDomain_):
      Ptr_F(&Func), DefinitionDomain(&_DefinitionDomain_)
    { };

    //! Evaluate functor at a given location (x,y)
    ReturnType operator()(const double &x, const double &y);

  private:
    FunctionType* Ptr_F;
    Polygon* DefinitionDomain;

    //! Private default constructor
    IntegrandFunctionOverPolygonalDefinitionDomain(void){};
  }; 

};


// Constructor with the Grid that is going to be associated with this object
template<class Grid2DQuadType> inline
Grid2DQuadIntegration<Grid2DQuadType>::Grid2DQuadIntegration(Grid2DQuadType * AssociatedGrid):
  CellFacesInfo(NULL),
  CellBndWest(NULL), CellBndSouth(NULL),
  CellBndEast(NULL), CellBndNorth(NULL)
{
  Grid = AssociatedGrid;
  // The 4 faces are in counterclockwise order W,S,E and N.
  CellFacesInfo = new bool [4];
  CellFacesInfo[0] = CellFacesInfo[1] = CellFacesInfo[2] = CellFacesInfo[3] = false;
}

// Deallocate memory
template<class Grid2DQuadType>
void Grid2DQuadIntegration<Grid2DQuadType>::Deallocate(void){
  if (CellFacesInfo != NULL){
    delete [] CellFacesInfo;
    CellFacesInfo = NULL;
  }
  CellBndWest = NULL; CellBndSouth = NULL;
  CellBndEast = NULL; CellBndNorth = NULL; 
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
 * Integrate a polynomial function (i.e. any function or pointer function)
 * over an arbitrary quadrilateral domain with straight edges.
 * This subroutine is suitable for integrating exactly & efficiently polynomials of up to order 18.
 *
 * \param FuncObj the function to be integrated
 * \param SW the South-West quadrilateral vertex
 * \param SE the South-East quadrilateral vertex
 * \param NE the North-East quadrilateral vertex
 * \param NW the North-West quadrilateral vertex
 * \param _dummy_param a parameter used only to determine the return type of the function FuncObj
 */
template<class Grid2DQuadType>
template<typename FO, typename Node2DType, class ReturnType> inline
ReturnType Grid2DQuadIntegration<Grid2DQuadType>::IntegratePolynomialOverQuadDomain(const FO FuncObj, 
										    const Node2DType &SW, const Node2DType &NW,
										    const Node2DType &NE, const Node2DType &SE,
										    ReturnType _dummy_param) const {

  return GaussLobattoQuadrilateralQuadrature(FuncObj,
					     SW, NW, NE, SE,
					     _dummy_param);
}

/*! 
 * Analyse the type of the cell and store 
 * the information in the designated variables.
 * For the time being, only interior cells are
 * flagged as curved if there is a curved boundary.
 * That is, the adjacent ghost cell is flagged with straight edges.
 */
template<class Grid2DQuadType>
void Grid2DQuadIntegration<Grid2DQuadType>::AnalyseCellFaces(const int &ii, const int &jj){
  
  // Reset variable
  AtLeastOneCurvedFace = false;
  CellBndWest = NULL; CellBndSouth = NULL;
  CellBndEast = NULL; CellBndNorth = NULL; 

  // == check if high-order boundary treatment is required
  if ( Grid->IsHighOrderBoundary() ){

    // Consider to treat faces as straight edges
    CellFacesInfo[0] = CellFacesInfo[1] = CellFacesInfo[2] = CellFacesInfo[3] = false;

    if (ii >= Grid->ICl && ii <= Grid->ICu && jj >= Grid->JCl && jj <= Grid->JCu ){
      // Analyze interior cells

      // Set correlation between face and block side
      IsTheSameBlockEdge = 1;

      // === Set the face types of this interior cell

      if (ii == Grid->ICl  && Grid->BndWestSplineInfo != NULL){
	// the West face is curved
	CellFacesInfo[0] = true;
	AtLeastOneCurvedFace = true;
	CellBndWest = &Grid->BndWestSpline;
      }

      if (jj == Grid->JCl  && Grid->BndSouthSplineInfo != NULL){
	// the South face is curved
	CellFacesInfo[1] = true;
	AtLeastOneCurvedFace = true;
	CellBndSouth = &Grid->BndSouthSpline;
      }

      if (ii == Grid->ICu  && Grid->BndEastSplineInfo != NULL){
	// the East face is curved
	CellFacesInfo[2] = true;
	AtLeastOneCurvedFace = true;
	CellBndEast = &Grid->BndEastSpline;
      }

      if (jj == Grid->JCu  && Grid->BndNorthSplineInfo != NULL){
	// the North face is curved
	CellFacesInfo[3] = true;
	AtLeastOneCurvedFace = true;
	CellBndNorth = &Grid->BndNorthSpline;
      }

    } else {
      // Analyze ghost cells
      // Attention: For a ghost cell the east face is analyzed with the west boundary and vice-versa.
      //            Same thing for the south and north faces.

      // Set correlation between face and block side
      IsTheSameBlockEdge = -1;

      // === Set the face types of this ghost cell

      // === Check West face ===
      if (ii == Grid->ICl){
	if (jj < Grid->JCl && Grid->ExtendSouth_BndWestSplineInfo != NULL){
	  // the South extension of West is curved
	  CellFacesInfo[0] = true;
	  AtLeastOneCurvedFace = true;
	  CellBndWest = &Grid->ExtendSouth_BndWestSpline;
	} else if (jj > Grid->JCu && Grid->ExtendNorth_BndWestSplineInfo != NULL){
	  // the North extension of West is curved
	  CellFacesInfo[0] = true;
	  AtLeastOneCurvedFace = true;
	  CellBndWest = &Grid->ExtendNorth_BndWestSpline;
	}
      } else if (ii == Grid->ICu+1){
	if (jj < Grid->JCl && Grid->ExtendSouth_BndEastSplineInfo != NULL){
	  // the South extension of East is curved
	  CellFacesInfo[0] = true;
	  AtLeastOneCurvedFace = true;
	  CellBndWest = &Grid->ExtendSouth_BndEastSpline;	  
	} else if (jj >= Grid->JCl && jj <= Grid->JCu && Grid->BndEastSplineInfo != NULL){
	  // the West face is curved
	  CellFacesInfo[0] = true;
	  AtLeastOneCurvedFace = true;
	  CellBndWest = &Grid->BndEastSpline;
	} else if (jj > Grid->JCu && Grid->ExtendNorth_BndEastSplineInfo != NULL){
	  // the North extension of East is curved
	  CellFacesInfo[0] = true;
	  AtLeastOneCurvedFace = true;
	  CellBndWest = &Grid->ExtendNorth_BndEastSpline;
	}
      }


      // === Check South face ===
      if (jj == Grid->JCl){
	if (ii < Grid->ICl && Grid->ExtendWest_BndSouthSplineInfo != NULL){
	  CellFacesInfo[1] = true;
	  AtLeastOneCurvedFace = true;
	  CellBndSouth = &Grid->ExtendWest_BndSouthSpline;
	} else if (ii > Grid->ICu && Grid->ExtendEast_BndSouthSplineInfo != NULL){
	  CellFacesInfo[1] = true;
	  AtLeastOneCurvedFace = true;
	  CellBndSouth = &Grid->ExtendEast_BndSouthSpline;
	}
      } else if (jj == Grid->JCu+1){
	if (ii < Grid->ICl && Grid->ExtendWest_BndNorthSplineInfo != NULL){
	  // the West extension of North is curved
	  CellFacesInfo[1] = true;
	  AtLeastOneCurvedFace = true;
	  CellBndSouth = &Grid->ExtendWest_BndNorthSpline;
	} else if ( ii >= Grid->ICl && ii <= Grid->ICu && Grid->BndNorthSplineInfo != NULL){
	  // the South face is curved
	  CellFacesInfo[1] = true;
	  AtLeastOneCurvedFace = true;
	  CellBndSouth = &Grid->BndNorthSpline;
	} else if (ii > Grid->ICu && Grid->ExtendEast_BndNorthSplineInfo != NULL) {
	  // the East extension of North is curved
	  CellFacesInfo[1] = true;
	  AtLeastOneCurvedFace = true;
	  CellBndSouth = &Grid->ExtendEast_BndNorthSpline;
	}
      }


      // === Check East face ===
      if (ii == Grid->ICl-1){
	if (jj < Grid->JCl && Grid->ExtendSouth_BndWestSplineInfo != NULL){
	  // the South extension of West is curved
	  CellFacesInfo[2] = true;
	  AtLeastOneCurvedFace = true;
	  CellBndEast = &Grid->ExtendSouth_BndWestSpline;
	} else if (jj >= Grid->JCl && jj <= Grid->JCu && Grid->BndWestSplineInfo != NULL){
	  // the East face is curved
	  CellFacesInfo[2] = true;
	  AtLeastOneCurvedFace = true;
	  CellBndEast = &Grid->BndWestSpline;
	} else if (jj > Grid->JCu && Grid->ExtendNorth_BndWestSplineInfo != NULL){
	  // the North extension of West is curved
	  CellFacesInfo[2] = true;
	  AtLeastOneCurvedFace = true;
	  CellBndEast = &Grid->ExtendNorth_BndWestSpline;
	}
      } else if (ii == Grid->ICu){
	if (jj < Grid->JCl && Grid->ExtendSouth_BndEastSplineInfo != NULL){
	  // the South extension of East is curved
	  CellFacesInfo[2] = true;
	  AtLeastOneCurvedFace = true;
	  CellBndEast = &Grid->ExtendSouth_BndEastSpline;
	} else if (jj > Grid->JCu && Grid->ExtendNorth_BndEastSplineInfo != NULL){
	  // the North extension of East is curved
	  CellFacesInfo[2] = true;
	  AtLeastOneCurvedFace = true;
	  CellBndEast = &Grid->ExtendNorth_BndEastSpline;	  
	}
      }


      // === Check North face ===
      if (jj == Grid->JCl-1){
	if (ii < Grid->ICl && Grid->ExtendWest_BndSouthSplineInfo != NULL){
	  // the West extension of South is curved
	  CellFacesInfo[3] = true;
	  AtLeastOneCurvedFace = true;
	  CellBndNorth = &Grid->ExtendWest_BndSouthSpline;
	} else if (ii >= Grid->ICl && ii <= Grid->ICu && Grid->BndSouthSplineInfo != NULL){
	  // the North face is curved
	  CellFacesInfo[3] = true;
	  AtLeastOneCurvedFace = true;
	  CellBndNorth = &Grid->BndSouthSpline;
	} else if (ii > Grid->ICu && Grid->ExtendEast_BndSouthSplineInfo != NULL){
	  // the East extension of South is curved
	  CellFacesInfo[3] = true;
	  AtLeastOneCurvedFace = true;
	  CellBndNorth = &Grid->ExtendEast_BndSouthSpline;
	}
      } else if (jj == Grid->JCu){
	if (ii < Grid->ICl && Grid->ExtendWest_BndNorthSplineInfo != NULL){
	  // the West extension of North is curved
	  CellFacesInfo[3] = true;
	  AtLeastOneCurvedFace = true;
	  CellBndNorth = &Grid->ExtendWest_BndNorthSpline;	  
	} else if (ii > Grid->ICu && Grid->ExtendEast_BndNorthSplineInfo != NULL){
	  // the East extension of North is curved
	  CellFacesInfo[3] = true;
	  AtLeastOneCurvedFace = true;
	  CellBndNorth = &Grid->ExtendEast_BndNorthSpline;
	}
      }


    }

  } else {
    // Set correlation between face and block side
    IsTheSameBlockEdge = 1;
    CellFacesInfo[0] = CellFacesInfo[1] = CellFacesInfo[2] = CellFacesInfo[3] = false;
  }
}

/*!
 * Print information to the output stream
 */
template<class Grid2DQuadType>
void Grid2DQuadIntegration<Grid2DQuadType>::PrintCellInfo(ostream & os) const {

  if (AtLeastOneCurvedFace){
    if (CellFacesInfo[0]){
      os << "\n West cell face is curved";
    }
    if (CellFacesInfo[1]){
      os << "\n South cell face is curved";
    }
    if (CellFacesInfo[2]){
      os << "\n East cell face is curved";
    }
    if (CellFacesInfo[3]){
      os << "\n North cell face is curved";
    }

  } else {
    os << "\n All faces of this cell are straight\n";
  }
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
 * \param ContourIntegrand the primitive with respect to x of FuncObj
 * \param digits the number of exact digits with which the result is computed (i.e. the accuracy of the calculation)
 * \param _dummy_param a parameter used only to determine the return type of the function FuncObj
 */
template<class Grid2DQuadType>
template<typename FO, class ReturnType> inline
ReturnType Grid2DQuadIntegration<Grid2DQuadType>::IntegrateFunctionOverCell(const int &ii, const int &jj, const FO FuncObj,
									    FO ContourIntegrand, int digits,
									    ReturnType _dummy_param) {

  try {
    // === Analyse cell faces
    AnalyseCellFaces(ii,jj);
  
    // === Decide whether to use contour integration or surface integration
    if (AtLeastOneCurvedFace){
      // This cell needs Gauss contour integration (i.e. at least one of the faces is curved)
      return IntegrateFunctionOverCellUsingContourIntegrand(ii,jj,
							    FuncObj,
							    ContourIntegrand,
							    digits,_dummy_param);
    } else {
      // This is a cell unaffected by curved boundaries or for which 
      // the presence of curved boundaries is ignored.
      return IntegrateFunctionOverCell(ii,jj,FuncObj,digits,_dummy_param);
    }

  } catch(std::runtime_error){
    
    if ( Grid->IsPolygonalAdaptiveQuadraturesAllowed() ){
      // Use adaptive quadrilateral quadratures to integrate the function (i.e. at least one of the faces is curved)
      return IntegrateFunctionOverCellUsingPolygonalAdaptiveQuadratures(ii,jj,
									FuncObj,
									digits,_dummy_param);

    } else if ( Grid->IsMonteCarloIntegrationAllowed() ){
      // Use Monte Carlo method to integrate the function (i.e. at least one of the faces is curved)
      return IntegrateFunctionOverCellUsingMonteCarloMethod(ii,jj,
							    FuncObj,
							    digits,_dummy_param);

    } else if ( Grid->IsIntegrationAlongCurvedBoundariesToleratedToGeometricInaccuracies() ){
      // Force the integration by treating the cell as being unaffected
      // by the presence of curved boundaries.
      return IntegrateFunctionOverCell(ii,jj,FuncObj,digits,_dummy_param);      

    } else {
      // Be intolerant to imperfect geometry treatment or other errors.
      throw;
    }
  }
  
}

/*!
 * Integrate a general integrand with the x-dependency integrated
 * along a straight segment line with respect to the y-coordinate.
 *
 * \param FuncObj the primitive with respect to x of the integrand
 * \param GaussQuads the array of Gauss quadrature points used for integration
 * \param GaussWeights the array of Gauss weights for the Gauss quadrature points
 * \param NumGQPs the number of GQPs used for integration (i.e. useful array size)
 * \param DeltaY the y-difference of the end points.
 * \param _dummy_param a parameter used only to determine the return type of the function FuncObj
 *
 * \note The relationship is customized for line segments! That's why DeltaY is used instead of dYdS and Length.
 */
template<class Grid2DQuadType>
template<typename FO, class ReturnType> inline
ReturnType Grid2DQuadIntegration<Grid2DQuadType>::
CalculateFunctionIntegralWithGaussQuadratures(FO FuncObj, 
					      const Vector2D * GaussQuads, const double * GaussWeights,
					      const int & NumGQPs, const double & DeltaY,
					      ReturnType _dummy_param) const {

  ReturnType Result(0);

  for (int n=0; n<NumGQPs; ++n){
    Result += GaussWeights[n] * FuncObj(GaussQuads[n].x, GaussQuads[n].y);
  }

  return Result * DeltaY;
}

/*! 
 * Compute the integral of the product between a scalar function and the normal vector.
 * Use Gauss-quadrature rule to calculate the integral.
 */
template<class Grid2DQuadType>
template<typename FO, class ReturnType>
void Grid2DQuadIntegration<Grid2DQuadType>::
IntegrateFunctionProjectionAlongBoundarySpline(const int & BOUNDARY, const FO FuncObj,
					       ReturnType & ResultXdir, ReturnType & ResultYdir,
					       double & WettedSurface) const {

  ReturnType TempXdir(0), TempYdir(0), FuncVal;
  
  // Set number of Gauss quadrature points per face used to compute the integral
  // to be the same as that specified in Spline2DInterval_HO
  int NumGQP(Grid2DQuadType::BndSplineIntervalType::get_NumGQPoints_ContourIntegral());

  Vector2D *GaussQuadPoints = new Vector2D [NumGQP]; // the GQPs at which the function is evaluated
  double * GaussQuadWeights = new double [NumGQP];   // the Gauss integration weights for each Gauss quadrature
  Vector2D Normal;
  typename Grid2DQuadType::BndSplineType SplineCopy;
  int iCell, jCell, GQPoint;


  /* Set the GaussQuadWeights. */
  GaussQuadratureData::getGaussQuadWeights(GaussQuadWeights, NumGQP);


  switch(BOUNDARY){
  case NORTH:			// North Bnd

    // Ensure correct boundary spline definition if high-order boundary representation required
    if ( Grid->BndNorthSplineInfo != NULL ){
      // Copy the spline
      SplineCopy = Grid->BndNorthSpline;

      // Check if the North Spline is defined such that the normals at the GaussQuadratures point outside of the domain 
      // (i.e The spline pathlength increases from INu to INl)
      if ( SplineCopy.getS(Grid->Node[Grid->INl][Grid->JNu]) < SplineCopy.getS(Grid->Node[Grid->INu][Grid->JNu]) ){
	// Change the direction of increasing the pathlength
	SplineCopy.Reverse_Spline();
      }

    } // endif

    // Calculate the contribution to the integration result of each cell on the North block boundary
    for (iCell=Grid->ICl; iCell<=Grid->ICu; ++iCell){

      // == Check for the representation of the geometric boundary (i.e. high-order or low-order)
      if ( Grid->BndNorthSplineInfo != NULL){
	/* High-order boundary representation is required.
	   Use all geometric information from the correspondent BndSplineInfo */

	// Add the contribution of this spline segment to the final results
	Grid->BndNorthSplineInfo[iCell].IntegrateFunctionProjectionOnNormalDirections(SplineCopy,
										      FuncObj,
										      ResultXdir, ResultYdir,
										      WettedSurface);
      } else {
	/* Low-order boundary representation is required. */

	// Reset temporal results
	TempXdir = ReturnType(0);
	TempYdir = ReturnType(0);

	// Determine the location of the Gauss Quadrature Points
	Grid->getGaussQuadPointsFaceN(iCell,Grid->JCu,GaussQuadPoints,NumGQP);
	// Determine normal
	Normal = Grid->nfaceN(iCell,Grid->JCu);
	
	for (GQPoint = 0; GQPoint < NumGQP; ++GQPoint) { // for each Gauss Quadrature point
	  
	  // Evaluate weighted function at the current GQP
	  FuncVal = GaussQuadWeights[GQPoint] * FuncObj(GaussQuadPoints[GQPoint].x,
							GaussQuadPoints[GQPoint].y);

	  // Update the X projection
	  TempXdir += FuncVal* Normal.x;
	  
	  // Update the Y projection
	  TempYdir += FuncVal* Normal.y;
	}

	// Update final result with the contribution of the current spline segment
	ResultXdir += TempXdir * Grid->lfaceN(iCell,Grid->JCu);
	ResultYdir += TempYdir * Grid->lfaceN(iCell,Grid->JCu);
	WettedSurface += Grid->lfaceN(iCell,Grid->JCu);
      }
    }
    break;

  case SOUTH:			// South Bnd
    // Ensure correct boundary spline definition if high-order boundary representation required
    if ( Grid->BndSouthSplineInfo != NULL ){
      // Copy the spline
      SplineCopy = Grid->BndSouthSpline;

      // Check if the South Spline is defined such that the normals at the GaussQuadratures point outside of the domain 
      // (i.e The spline pathlength increases from INl to INu)
      if ( SplineCopy.getS(Grid->Node[Grid->INl][Grid->JNl]) > SplineCopy.getS(Grid->Node[Grid->INu][Grid->JNl]) ){
	// Change the direction of increasing the pathlength
	SplineCopy.Reverse_Spline();
      }

    } // endif

    // Calculate the contribution to the integration result of each cell on the South block boundary
    for (iCell=Grid->ICl; iCell<=Grid->ICu; ++iCell){

      // == Check for the representation of the geometric boundary (i.e. high-order or low-order)
      if ( Grid->BndSouthSplineInfo != NULL){
	/* High-order boundary representation is required.
	   Use all geometric information from the correspondent BndSplineInfo */

	// Add the contribution of this spline segment to the final results
	Grid->BndSouthSplineInfo[iCell].IntegrateFunctionProjectionOnNormalDirections(SplineCopy,
										      FuncObj,
										      ResultXdir, ResultYdir,
										      WettedSurface);
      } else {
	/* Low-order boundary representation is required. */

	// Reset temporal results
	TempXdir = ReturnType(0);
	TempYdir = ReturnType(0);

	// Determine the location of the Gauss Quadrature Points
	Grid->getGaussQuadPointsFaceS(iCell,Grid->JCl,GaussQuadPoints,NumGQP);
	// Determine normal
	Normal = Grid->nfaceS(iCell,Grid->JCl);
	
	for (GQPoint = 0; GQPoint < NumGQP; ++GQPoint) { // for each Gauss Quadrature point
	  
	  // Evaluate weighted function at the current GQP
	  FuncVal = GaussQuadWeights[GQPoint] * FuncObj(GaussQuadPoints[GQPoint].x,
							GaussQuadPoints[GQPoint].y);

	  // Update the X projection
	  TempXdir += FuncVal* Normal.x;
	  
	  // Update the Y projection
	  TempYdir += FuncVal* Normal.y;
	}

	// Update final result with the contribution of the current spline segment
	ResultXdir += TempXdir * Grid->lfaceS(iCell,Grid->JCl);
	ResultYdir += TempYdir * Grid->lfaceS(iCell,Grid->JCl);
	WettedSurface += Grid->lfaceS(iCell,Grid->JCl);
      }
    }
    break;

  case EAST:			// East Bnd
    // Ensure correct boundary spline definition if high-order boundary representation required
    if ( Grid->BndEastSplineInfo != NULL ){
      // Copy the spline
      SplineCopy = Grid->BndEastSpline;

      // Check if the East Spline is defined such that the normals at the GaussQuadratures point outside of the domain 
      // (i.e The spline pathlength increases from JNl to JNu)
      if ( SplineCopy.getS(Grid->Node[Grid->INu][Grid->JNl]) > SplineCopy.getS(Grid->Node[Grid->INu][Grid->JNu]) ){
	// Change the direction of increasing the pathlength
	SplineCopy.Reverse_Spline();
      }

    } // endif

    // Calculate the contribution to the integration result of each cell on the East block boundary
    for (jCell=Grid->JCl; jCell<=Grid->JCu; ++jCell){

      // == Check for the representation of the geometric boundary (i.e. high-order or low-order)
      if ( Grid->BndEastSplineInfo != NULL){
	/* High-order boundary representation is required.
	   Use all geometric information from the correspondent BndSplineInfo */

	// Add the contribution of this spline segment to the final results
	Grid->BndEastSplineInfo[jCell].IntegrateFunctionProjectionOnNormalDirections(SplineCopy,
										     FuncObj,
										     ResultXdir, ResultYdir,
										     WettedSurface);
      } else {
	/* Low-order boundary representation is required. */

	// Reset temporal results
	TempXdir = ReturnType(0);
	TempYdir = ReturnType(0);

	// Determine the location of the Gauss Quadrature Points
	Grid->getGaussQuadPointsFaceE(Grid->ICu,jCell,GaussQuadPoints,NumGQP);
	// Determine normal
	Normal = Grid->nfaceE(Grid->ICu,jCell);
	
	for (GQPoint = 0; GQPoint < NumGQP; ++GQPoint) { // for each Gauss Quadrature point
	  
	  // Evaluate weighted function at the current GQP
	  FuncVal = GaussQuadWeights[GQPoint] * FuncObj(GaussQuadPoints[GQPoint].x,
							GaussQuadPoints[GQPoint].y);

	  // Update the X projection
	  TempXdir += FuncVal* Normal.x;
	  
	  // Update the Y projection
	  TempYdir += FuncVal* Normal.y;
	}

	// Update final result with the contribution of the current spline segment
	ResultXdir += TempXdir * Grid->lfaceE(Grid->ICu,jCell);
	ResultYdir += TempYdir * Grid->lfaceE(Grid->ICu,jCell);
	WettedSurface += Grid->lfaceE(Grid->ICu,jCell);
      }
    }
    break;

  case WEST:			// West Bnd
    // Ensure correct boundary spline definition if high-order boundary representation required
    if ( Grid->BndWestSplineInfo != NULL ){
      // Copy the spline
      SplineCopy = Grid->BndWestSpline;

      // Check if the West Spline is defined such that the normals at the GaussQuadratures point outside of the domain 
      // (i.e The spline pathlength increases from JNu to JNl)
      if ( SplineCopy.getS(Grid->Node[Grid->INl][Grid->JNl]) < SplineCopy.getS(Grid->Node[Grid->INl][Grid->JNu]) ){
	// Change the direction of increasing the pathlength
	SplineCopy.Reverse_Spline();
      }

    } // endif

    // Calculate the contribution to the integration result of each cell on the West block boundary
    for (jCell=Grid->JCl; jCell<=Grid->JCu; ++jCell){

      // == Check for the representation of the geometric boundary (i.e. high-order or low-order)
      if ( Grid->BndWestSplineInfo != NULL){
	/* High-order boundary representation is required.
	   Use all geometric information from the correspondent BndSplineInfo */

	// Add the contribution of this spline segment to the final results
	Grid->BndWestSplineInfo[jCell].IntegrateFunctionProjectionOnNormalDirections(SplineCopy,
										     FuncObj,
										     ResultXdir, ResultYdir,
										     WettedSurface);
      } else {
	/* Low-order boundary representation is required. */

	// Reset temporal results
	TempXdir = ReturnType(0);
	TempYdir = ReturnType(0);
	
	// Determine the location of the Gauss Quadrature Points
	Grid->getGaussQuadPointsFaceW(Grid->ICl,jCell,GaussQuadPoints,NumGQP);
	// Determine normal
	Normal = Grid->nfaceW(Grid->ICl,jCell);
	
	for (GQPoint = 0; GQPoint < NumGQP; ++GQPoint) { // for each Gauss Quadrature point
	  
	  // Evaluate weighted function at the current GQP
	  FuncVal = GaussQuadWeights[GQPoint] * FuncObj(GaussQuadPoints[GQPoint].x,
							GaussQuadPoints[GQPoint].y);

	  // Update the X projection
	  TempXdir += FuncVal* Normal.x;
	  
	  // Update the Y projection
	  TempYdir += FuncVal* Normal.y;
	}

	// Update final result with the contribution of the current spline segment
	ResultXdir += TempXdir * Grid->lfaceW(Grid->ICl,jCell);
	ResultYdir += TempYdir * Grid->lfaceW(Grid->ICl,jCell);
	WettedSurface += Grid->lfaceW(Grid->ICl,jCell);
      }
    }
    break;
  }

}


/*! 
 * Compute the integral of the product between a piecewise scalar function and the normal vector.
 * Use Gauss-quadrature rule to calculate the integral.
 *
 * \note FuncObj is an object that has a member function called NewIndexes()!
 */
template<class Grid2DQuadType>
template<typename FO, class ReturnType>
void Grid2DQuadIntegration<Grid2DQuadType>::
IntegratePiecewiseFunctionProjectionAlongBoundarySpline(const int & BOUNDARY, FO FuncObj,
							ReturnType & ResultXdir, ReturnType & ResultYdir,
							double & WettedSurface) const {

  ReturnType TempXdir(0), TempYdir(0), FuncVal;
  
  // Set number of Gauss quadrature points per face used to compute the integral
  // to be the same as that specified in Spline2DInterval_HO
  int NumGQP(Grid2DQuadType::BndSplineIntervalType::get_NumGQPoints_ContourIntegral());

  Vector2D *GaussQuadPoints = new Vector2D [NumGQP]; // the GQPs at which the function is evaluated
  double * GaussQuadWeights = new double [NumGQP];   // the Gauss integration weights for each Gauss quadrature
  Vector2D Normal;
  typename Grid2DQuadType::BndSplineType SplineCopy;
  int iCell, jCell, GQPoint;


  /* Set the GaussQuadWeights. */
  GaussQuadratureData::getGaussQuadWeights(GaussQuadWeights, NumGQP);


  switch(BOUNDARY){
  case NORTH:			// North Bnd

    // Ensure correct boundary spline definition if high-order boundary representation required
    if ( Grid->BndNorthSplineInfo != NULL ){
      // Copy the spline
      SplineCopy = Grid->BndNorthSpline;

      // Check if the North Spline is defined such that the normals at the GaussQuadratures point outside of the domain 
      // (i.e The spline pathlength increases from INu to INl)
      if ( SplineCopy.getS(Grid->Node[Grid->INl][Grid->JNu]) < SplineCopy.getS(Grid->Node[Grid->INu][Grid->JNu]) ){
	// Change the direction of increasing the pathlength
	SplineCopy.Reverse_Spline();
      }

    } // endif

    // Calculate the contribution to the integration result of each cell on the North block boundary
    for (iCell=Grid->ICl; iCell<=Grid->ICu; ++iCell){

      // Change indexes of FuncObj
      FuncObj.NewIndexes(iCell,Grid->JCu);

      // == Check for the representation of the geometric boundary (i.e. high-order or low-order)
      if ( Grid->BndNorthSplineInfo != NULL){
	/* High-order boundary representation is required.
	   Use all geometric information from the correspondent BndSplineInfo */

	// Add the contribution of this spline segment to the final results
	Grid->BndNorthSplineInfo[iCell].IntegrateFunctionProjectionOnNormalDirections(SplineCopy,
										      FuncObj,
										      ResultXdir, ResultYdir,
										      WettedSurface);
      } else {
	/* Low-order boundary representation is required. */

	// Reset temporal results
	TempXdir = ReturnType(0);
	TempYdir = ReturnType(0);

	// Determine the location of the Gauss Quadrature Points
	Grid->getGaussQuadPointsFaceN(iCell,Grid->JCu,GaussQuadPoints,NumGQP);
	// Determine normal
	Normal = Grid->nfaceN(iCell,Grid->JCu);
	
	for (GQPoint = 0; GQPoint < NumGQP; ++GQPoint) { // for each Gauss Quadrature point
	  
	  // Evaluate weighted function at the current GQP
	  FuncVal = GaussQuadWeights[GQPoint] * FuncObj(GaussQuadPoints[GQPoint].x,
							GaussQuadPoints[GQPoint].y);

	  // Update the X projection
	  TempXdir += FuncVal* Normal.x;
	  
	  // Update the Y projection
	  TempYdir += FuncVal* Normal.y;
	}

	// Update final result with the contribution of the current spline segment
	ResultXdir += TempXdir * Grid->lfaceN(iCell,Grid->JCu);
	ResultYdir += TempYdir * Grid->lfaceN(iCell,Grid->JCu);
	WettedSurface += Grid->lfaceN(iCell,Grid->JCu);
      }
    }
    break;

  case SOUTH:			// South Bnd
    // Ensure correct boundary spline definition if high-order boundary representation required
    if ( Grid->BndSouthSplineInfo != NULL ){
      // Copy the spline
      SplineCopy = Grid->BndSouthSpline;

      // Check if the South Spline is defined such that the normals at the GaussQuadratures point outside of the domain 
      // (i.e The spline pathlength increases from INl to INu)
      if ( SplineCopy.getS(Grid->Node[Grid->INl][Grid->JNl]) > SplineCopy.getS(Grid->Node[Grid->INu][Grid->JNl]) ){
	// Change the direction of increasing the pathlength
	SplineCopy.Reverse_Spline();
      }

    } // endif

    // Calculate the contribution to the integration result of each cell on the South block boundary
    for (iCell=Grid->ICl; iCell<=Grid->ICu; ++iCell){

      // Change indexes of FuncObj
      FuncObj.NewIndexes(iCell,Grid->JCl);

      // == Check for the representation of the geometric boundary (i.e. high-order or low-order)
      if ( Grid->BndSouthSplineInfo != NULL){
	/* High-order boundary representation is required.
	   Use all geometric information from the correspondent BndSplineInfo */

	// Add the contribution of this spline segment to the final results
	Grid->BndSouthSplineInfo[iCell].IntegrateFunctionProjectionOnNormalDirections(SplineCopy,
										      FuncObj,
										      ResultXdir, ResultYdir,
										      WettedSurface);
      } else {
	/* Low-order boundary representation is required. */

	// Reset temporal results
	TempXdir = ReturnType(0);
	TempYdir = ReturnType(0);

	// Determine the location of the Gauss Quadrature Points
	Grid->getGaussQuadPointsFaceS(iCell,Grid->JCl,GaussQuadPoints,NumGQP);
	// Determine normal
	Normal = Grid->nfaceS(iCell,Grid->JCl);
	
	for (GQPoint = 0; GQPoint < NumGQP; ++GQPoint) { // for each Gauss Quadrature point
	  
	  // Evaluate weighted function at the current GQP
	  FuncVal = GaussQuadWeights[GQPoint] * FuncObj(GaussQuadPoints[GQPoint].x,
							GaussQuadPoints[GQPoint].y);

	  // Update the X projection
	  TempXdir += FuncVal* Normal.x;
	  
	  // Update the Y projection
	  TempYdir += FuncVal* Normal.y;
	}

	// Update final result with the contribution of the current spline segment
	ResultXdir += TempXdir * Grid->lfaceS(iCell,Grid->JCl);
	ResultYdir += TempYdir * Grid->lfaceS(iCell,Grid->JCl);
	WettedSurface += Grid->lfaceS(iCell,Grid->JCl);
      }
    }
    break;

  case EAST:			// East Bnd
    // Ensure correct boundary spline definition if high-order boundary representation required
    if ( Grid->BndEastSplineInfo != NULL ){
      // Copy the spline
      SplineCopy = Grid->BndEastSpline;

      // Check if the East Spline is defined such that the normals at the GaussQuadratures point outside of the domain 
      // (i.e The spline pathlength increases from JNl to JNu)
      if ( SplineCopy.getS(Grid->Node[Grid->INu][Grid->JNl]) > SplineCopy.getS(Grid->Node[Grid->INu][Grid->JNu]) ){
	// Change the direction of increasing the pathlength
	SplineCopy.Reverse_Spline();
      }

    } // endif

    // Calculate the contribution to the integration result of each cell on the East block boundary
    for (jCell=Grid->JCl; jCell<=Grid->JCu; ++jCell){

      // Change indexes of FuncObj
      FuncObj.NewIndexes(Grid->ICu,jCell);

      // == Check for the representation of the geometric boundary (i.e. high-order or low-order)
      if ( Grid->BndEastSplineInfo != NULL){
	/* High-order boundary representation is required.
	   Use all geometric information from the correspondent BndSplineInfo */

	// Add the contribution of this spline segment to the final results
	Grid->BndEastSplineInfo[jCell].IntegrateFunctionProjectionOnNormalDirections(SplineCopy,
										     FuncObj,
										     ResultXdir, ResultYdir,
										     WettedSurface);
      } else {
	/* Low-order boundary representation is required. */

	// Reset temporal results
	TempXdir = ReturnType(0);
	TempYdir = ReturnType(0);

	// Determine the location of the Gauss Quadrature Points
	Grid->getGaussQuadPointsFaceE(Grid->ICu,jCell,GaussQuadPoints,NumGQP);
	// Determine normal
	Normal = Grid->nfaceE(Grid->ICu,jCell);
	
	for (GQPoint = 0; GQPoint < NumGQP; ++GQPoint) { // for each Gauss Quadrature point
	  
	  // Evaluate weighted function at the current GQP
	  FuncVal = GaussQuadWeights[GQPoint] * FuncObj(GaussQuadPoints[GQPoint].x,
							GaussQuadPoints[GQPoint].y);

	  // Update the X projection
	  TempXdir += FuncVal* Normal.x;
	  
	  // Update the Y projection
	  TempYdir += FuncVal* Normal.y;
	}

	// Update final result with the contribution of the current spline segment
	ResultXdir += TempXdir * Grid->lfaceE(Grid->ICu,jCell);
	ResultYdir += TempYdir * Grid->lfaceE(Grid->ICu,jCell);
	WettedSurface += Grid->lfaceE(Grid->ICu,jCell);
      }
    }
    break;

  case WEST:			// West Bnd
    // Ensure correct boundary spline definition if high-order boundary representation required
    if ( Grid->BndWestSplineInfo != NULL ){
      // Copy the spline
      SplineCopy = Grid->BndWestSpline;

      // Check if the West Spline is defined such that the normals at the GaussQuadratures point outside of the domain 
      // (i.e The spline pathlength increases from JNu to JNl)
      if ( SplineCopy.getS(Grid->Node[Grid->INl][Grid->JNl]) < SplineCopy.getS(Grid->Node[Grid->INl][Grid->JNu]) ){
	// Change the direction of increasing the pathlength
	SplineCopy.Reverse_Spline();
      }

    } // endif

    // Calculate the contribution to the integration result of each cell on the West block boundary
    for (jCell=Grid->JCl; jCell<=Grid->JCu; ++jCell){

      // Change indexes of FuncObj
      FuncObj.NewIndexes(Grid->ICl,jCell);

      // == Check for the representation of the geometric boundary (i.e. high-order or low-order)
      if ( Grid->BndWestSplineInfo != NULL){
	/* High-order boundary representation is required.
	   Use all geometric information from the correspondent BndSplineInfo */

	// Add the contribution of this spline segment to the final results
	Grid->BndWestSplineInfo[jCell].IntegrateFunctionProjectionOnNormalDirections(SplineCopy,
										     FuncObj,
										     ResultXdir, ResultYdir,
										     WettedSurface);
      } else {
	/* Low-order boundary representation is required. */

	// Reset temporal results
	TempXdir = ReturnType(0);
	TempYdir = ReturnType(0);
	
	// Determine the location of the Gauss Quadrature Points
	Grid->getGaussQuadPointsFaceW(Grid->ICl,jCell,GaussQuadPoints,NumGQP);
	// Determine normal
	Normal = Grid->nfaceW(Grid->ICl,jCell);
	
	for (GQPoint = 0; GQPoint < NumGQP; ++GQPoint) { // for each Gauss Quadrature point
	  
	  // Evaluate weighted function at the current GQP
	  FuncVal = GaussQuadWeights[GQPoint] * FuncObj(GaussQuadPoints[GQPoint].x,
							GaussQuadPoints[GQPoint].y);

	  // Update the X projection
	  TempXdir += FuncVal* Normal.x;
	  
	  // Update the Y projection
	  TempYdir += FuncVal* Normal.y;
	}

	// Update final result with the contribution of the current spline segment
	ResultXdir += TempXdir * Grid->lfaceW(Grid->ICl,jCell);
	ResultYdir += TempYdir * Grid->lfaceW(Grid->ICl,jCell);
	WettedSurface += Grid->lfaceW(Grid->ICl,jCell);
      }
    }
    break;
  }

  // Delete memory
  delete [] GaussQuadPoints;
  delete [] GaussQuadWeights;

}


/*! 
 * Compute the integral of the wall shear stress along a boundary 
 * (i.e. using information from each cell).
 * Use Gauss-quadrature rule to calculate the integral.
 *
 * \note FuncObj is an object that has a member function called NewIndexes()!
 */
template<class Grid2DQuadType>
template<typename FO, typename FO_TestValidityDomain>
void Grid2DQuadIntegration<Grid2DQuadType>::
IntegratePiecewiseWallShearStressAlongBoundarySpline(const int & BOUNDARY, FO FuncObj,
						     double & ResultXdir, double & ResultYdir,
						     double & WettedSurface,
						     FO_TestValidityDomain & ValidateDomain) const {

  double TempXdir(0), TempYdir(0), FuncVal;
  
  // Set number of Gauss quadrature points per face used to compute the integral
  // to be the same as that specified in Spline2DInterval_HO
  int NumGQP(Grid2DQuadType::BndSplineIntervalType::get_NumGQPoints_ContourIntegral());

  Vector2D *GaussQuadPoints = new Vector2D [NumGQP]; // the GQPs at which the function is evaluated
  double * GaussQuadWeights = new double [NumGQP];   // the Gauss integration weights for each Gauss quadrature
  Vector2D Normal, Tangent;
  typename Grid2DQuadType::BndSplineType SplineCopy;
  int iCell, jCell, GQPoint;


  /* Set the GaussQuadWeights. */
  GaussQuadratureData::getGaussQuadWeights(GaussQuadWeights, NumGQP);


  switch(BOUNDARY){
  case NORTH:			// North Bnd

    // Ensure correct boundary spline definition if high-order boundary representation required
    if ( Grid->BndNorthSplineInfo != NULL ){
      // Copy the spline
      SplineCopy = Grid->BndNorthSpline;

      // Check if the North Spline is defined such that the normals at the GaussQuadratures point outside of the domain 
      // (i.e The spline pathlength increases from INu to INl)
      if ( SplineCopy.getS(Grid->Node[Grid->INl][Grid->JNu]) < SplineCopy.getS(Grid->Node[Grid->INu][Grid->JNu]) ){
	// Change the direction of increasing the pathlength
	SplineCopy.Reverse_Spline();
      }

    } // endif

    // Calculate the contribution to the integration result of each cell on the North block boundary
    for (iCell=Grid->ICl; iCell<=Grid->ICu; ++iCell){

      // == Check if the current cell position passes the test function. The cell centroid is used for the test!
      if ( ValidateDomain.PointInIntegrationDomainTest(Grid->CellCentroid(iCell,Grid->JCu)) ){

	// Change indexes of FuncObj
	FuncObj.NewIndexes(iCell,Grid->JCu);

	// == Check for the representation of the geometric boundary (i.e. high-order or low-order)
	if ( Grid->BndNorthSplineInfo != NULL){
	  /* High-order boundary representation is required.
	     Use all geometric information from the correspondent BndSplineInfo */

	  // Add the contribution of this spline segment to the final results
	  Grid->BndNorthSplineInfo[iCell].IntegrateWallShearStressContributions(SplineCopy,
										FuncObj,
										ResultXdir, ResultYdir,
										WettedSurface);
	} else {
	  /* Low-order boundary representation is required. */

	  // Reset temporal results
	  TempXdir = 0;
	  TempYdir = 0;

	  // Determine the location of the Gauss Quadrature Points
	  Grid->getGaussQuadPointsFaceN(iCell,Grid->JCu,GaussQuadPoints,NumGQP);
	  // Determine normal direction pointing outward from the body
	  Normal = -Grid->nfaceN(iCell,Grid->JCu);
	  // Determine the tangent
	  Grid->getTangent(Tangent,Normal);
	
	  for (GQPoint = 0; GQPoint < NumGQP; ++GQPoint) { // for each Gauss Quadrature point
	  
	    // Evaluate weighted function at the current GQP
	    FuncVal = GaussQuadWeights[GQPoint] * FuncObj(GaussQuadPoints[GQPoint],
							  Normal);

	    // Update the X projection
	    TempXdir += FuncVal* Tangent.x;
	  
	    // Update the Y projection
	    TempYdir += FuncVal* Tangent.y;
	  }

	  // Update final result with the contribution of the current spline segment
	  ResultXdir += TempXdir * Grid->lfaceN(iCell,Grid->JCu);
	  ResultYdir += TempYdir * Grid->lfaceN(iCell,Grid->JCu);
	  WettedSurface += Grid->lfaceN(iCell,Grid->JCu);
	}
      }	// endif

    }
    break;

  case SOUTH:			// South Bnd

    // Ensure correct boundary spline definition if high-order boundary representation required
    if ( Grid->BndSouthSplineInfo != NULL ){
      // Copy the spline
      SplineCopy = Grid->BndSouthSpline;

      // Check if the South Spline is defined such that the normals at the GaussQuadratures point outside of the domain 
      // (i.e The spline pathlength increases from INl to INu)
      if ( SplineCopy.getS(Grid->Node[Grid->INl][Grid->JNl]) > SplineCopy.getS(Grid->Node[Grid->INu][Grid->JNl]) ){
	// Change the direction of increasing the pathlength
	SplineCopy.Reverse_Spline();
      }

    } // endif

    // Calculate the contribution to the integration result of each cell on the South block boundary
    for (iCell=Grid->ICl; iCell<=Grid->ICu; ++iCell){

      // == Check if the current cell position passes the test function. The cell centroid is used for the test!
      if ( ValidateDomain.PointInIntegrationDomainTest(Grid->CellCentroid(iCell,Grid->JCl)) ){

	// Change indexes of FuncObj
	FuncObj.NewIndexes(iCell,Grid->JCl);

	// == Check for the representation of the geometric boundary (i.e. high-order or low-order)
	if ( Grid->BndSouthSplineInfo != NULL){
	  /* High-order boundary representation is required.
	     Use all geometric information from the correspondent BndSplineInfo */

	  // Add the contribution of this spline segment to the final results
	  Grid->BndSouthSplineInfo[iCell].IntegrateWallShearStressContributions(SplineCopy,
										FuncObj,
										ResultXdir, ResultYdir,
										WettedSurface);
	} else {
	  /* Low-order boundary representation is required. */

	  // Reset temporal results
	  TempXdir = 0;
	  TempYdir = 0;

	  // Determine the location of the Gauss Quadrature Points
	  Grid->getGaussQuadPointsFaceS(iCell,Grid->JCl,GaussQuadPoints,NumGQP);
	  // Determine normal direction pointing outward from the body
	  Normal = -Grid->nfaceS(iCell,Grid->JCl);
	  // Determine the tangent
	  Grid->getTangent(Tangent,Normal);
	
	  for (GQPoint = 0; GQPoint < NumGQP; ++GQPoint) { // for each Gauss Quadrature point
	  
	    // Evaluate weighted function at the current GQP
	    FuncVal = GaussQuadWeights[GQPoint] * FuncObj(GaussQuadPoints[GQPoint],
							  Normal);

	    // Update the X projection
	    TempXdir += FuncVal* Tangent.x; 
	  
	    // Update the Y projection
	    TempYdir += FuncVal* Tangent.y;
	  }

	  // Update final result with the contribution of the current spline segment
	  ResultXdir += TempXdir * Grid->lfaceS(iCell,Grid->JCl);
	  ResultYdir += TempYdir * Grid->lfaceS(iCell,Grid->JCl);
	  WettedSurface += Grid->lfaceS(iCell,Grid->JCl);
	}
      }	// endif

    }
    break;

  case EAST:			// East Bnd
    // Ensure correct boundary spline definition if high-order boundary representation required
    if ( Grid->BndEastSplineInfo != NULL ){
      // Copy the spline
      SplineCopy = Grid->BndEastSpline;

      // Check if the East Spline is defined such that the normals at the GaussQuadratures point outside of the domain 
      // (i.e The spline pathlength increases from JNl to JNu)
      if ( SplineCopy.getS(Grid->Node[Grid->INu][Grid->JNl]) > SplineCopy.getS(Grid->Node[Grid->INu][Grid->JNu]) ){
	// Change the direction of increasing the pathlength
	SplineCopy.Reverse_Spline();
      }

    } // endif

    // Calculate the contribution to the integration result of each cell on the East block boundary
    for (jCell=Grid->JCl; jCell<=Grid->JCu; ++jCell){

      // == Check if the current cell position passes the test function. The cell centroid is used for the test!
      if ( ValidateDomain.PointInIntegrationDomainTest(Grid->CellCentroid(Grid->ICu,jCell)) ){

	// Change indexes of FuncObj
	FuncObj.NewIndexes(Grid->ICu,jCell);

	// == Check for the representation of the geometric boundary (i.e. high-order or low-order)
	if ( Grid->BndEastSplineInfo != NULL){
	  /* High-order boundary representation is required.
	     Use all geometric information from the correspondent BndSplineInfo */

	  // Add the contribution of this spline segment to the final results
	  Grid->BndEastSplineInfo[jCell].IntegrateWallShearStressContributions(SplineCopy,
									       FuncObj,
									       ResultXdir, ResultYdir,
									       WettedSurface);
	} else {
	  /* Low-order boundary representation is required. */

	  // Reset temporal results
	  TempXdir = 0;
	  TempYdir = 0;

	  // Determine the location of the Gauss Quadrature Points
	  Grid->getGaussQuadPointsFaceE(Grid->ICu,jCell,GaussQuadPoints,NumGQP);
	  // Determine normal direction pointing outward from the body
	  Normal = -Grid->nfaceE(Grid->ICu,jCell);
	  // Determine the tangent
	  Grid->getTangent(Tangent,Normal);
	
	  for (GQPoint = 0; GQPoint < NumGQP; ++GQPoint) { // for each Gauss Quadrature point
	  
	    // Evaluate weighted function at the current GQP
	    FuncVal = GaussQuadWeights[GQPoint] * FuncObj(GaussQuadPoints[GQPoint],
							  Normal);

	    // Update the X projection
	    TempXdir += FuncVal* Tangent.x;
	  
	    // Update the Y projection
	    TempYdir += FuncVal* Tangent.y;
	  }

	  // Update final result with the contribution of the current spline segment
	  ResultXdir += TempXdir * Grid->lfaceE(Grid->ICu,jCell);
	  ResultYdir += TempYdir * Grid->lfaceE(Grid->ICu,jCell);
	  WettedSurface += Grid->lfaceE(Grid->ICu,jCell);
	}
      }	// endif

    }
    break;

  case WEST:			// West Bnd
    // Ensure correct boundary spline definition if high-order boundary representation required
    if ( Grid->BndWestSplineInfo != NULL ){
      // Copy the spline
      SplineCopy = Grid->BndWestSpline;

      // Check if the West Spline is defined such that the normals at the GaussQuadratures point outside of the domain 
      // (i.e The spline pathlength increases from JNu to JNl)
      if ( SplineCopy.getS(Grid->Node[Grid->INl][Grid->JNl]) < SplineCopy.getS(Grid->Node[Grid->INl][Grid->JNu]) ){
	// Change the direction of increasing the pathlength
	SplineCopy.Reverse_Spline();
      }

    } // endif

    // Calculate the contribution to the integration result of each cell on the West block boundary
    for (jCell=Grid->JCl; jCell<=Grid->JCu; ++jCell){

      // == Check if the current cell position passes the test function. The cell centroid is used for the test!
      if ( ValidateDomain.PointInIntegrationDomainTest(Grid->CellCentroid(Grid->ICl,jCell)) ){

	// Change indexes of FuncObj
	FuncObj.NewIndexes(Grid->ICl,jCell);

	// == Check for the representation of the geometric boundary (i.e. high-order or low-order)
	if ( Grid->BndWestSplineInfo != NULL){
	  /* High-order boundary representation is required.
	     Use all geometric information from the correspondent BndSplineInfo */

	  // Add the contribution of this spline segment to the final results
	  Grid->BndWestSplineInfo[jCell].IntegrateWallShearStressContributions(SplineCopy,
									       FuncObj,
									       ResultXdir, ResultYdir,
									       WettedSurface);
	} else {
	  /* Low-order boundary representation is required. */

	  // Reset temporal results
	  TempXdir = 0;
	  TempYdir = 0;
	
	  // Determine the location of the Gauss Quadrature Points
	  Grid->getGaussQuadPointsFaceW(Grid->ICl,jCell,GaussQuadPoints,NumGQP);
	  // Determine normal direction pointing outward from the body
	  Normal = -Grid->nfaceW(Grid->ICl,jCell);
	  // Determine the tangent
	  Grid->getTangent(Tangent,Normal);
	
	  for (GQPoint = 0; GQPoint < NumGQP; ++GQPoint) { // for each Gauss Quadrature point
	  
	    // Evaluate weighted function at the current GQP
	    FuncVal = GaussQuadWeights[GQPoint] * FuncObj(GaussQuadPoints[GQPoint],
							  Normal);

	    // Update the X projection
	    TempXdir += FuncVal* Tangent.x;
	  
	    // Update the Y projection
	    TempYdir += FuncVal* Tangent.y;
	  }

	  // Update final result with the contribution of the current spline segment
	  ResultXdir += TempXdir * Grid->lfaceW(Grid->ICl,jCell);
	  ResultYdir += TempYdir * Grid->lfaceW(Grid->ICl,jCell);
	  WettedSurface += Grid->lfaceW(Grid->ICl,jCell);
	}
      }	// endif

    }
    break;
  }

  // Delete memory
  delete [] GaussQuadPoints;
  delete [] GaussQuadWeights;

}


/**************************************************************************************
 * Implement member functions of IntegrandFunctionOverPolygonalDefinitionDomain class *
 **************************************************************************************/
//! Main constructor
template<class Grid2DQuadType>
template<class FunctionType, class ReturnType>
ReturnType Grid2DQuadIntegration<Grid2DQuadType>::
IntegrandFunctionOverPolygonalDefinitionDomain<FunctionType,ReturnType>::operator()(const double &x, const double &y){

  if (DefinitionDomain->IsPointInPolygon(Vector2D(x,y))){
    return (*Ptr_F)(x,y);
  } else {
    return ReturnType(0.0);
  }

}

#endif

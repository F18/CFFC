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

  //! Access to type of cell faces
  const vector<bool> getCellFacesInfo(void){ return CellFacesInfo; }
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


private:
  Grid2DQuadType *Grid;	        //!< pointer to the grid associated to this object

  Grid2DQuadIntegration(void);	//!< Private default constructor  

  vector<bool> CellFacesInfo;	//!< Array to track the type of each cell face
  int IsTheSameBlockEdge;       /*!< Variable to mark the association between the curved cell face and the block side. 
				  (e.g a West cell face corresponds to a West or East block side)
				  It basically makes the distinction between interior cells and ghost cells. 
				  This variable takes value +1 for interior cells and -1 for ghost cells. */
  bool AtLeastOneCurvedFace;	//!< Flag to indicate whether the cell does have at least one curved face
};


// Constructor with the Grid that is going to be associated with this object
template<class Grid2DQuadType> inline
Grid2DQuadIntegration<Grid2DQuadType>::Grid2DQuadIntegration(Grid2DQuadType * AssociatedGrid){
  Grid = AssociatedGrid;
  CellFacesInfo.assign(4, false); // The 4 faces are in counterclockwise order W,S,E and N.
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

  // == check if high-order boundary treatment is required
  if ( Grid->IsHighOrderBoundary() ){

    if (ii >= Grid->ICl && ii <= Grid->ICu && jj >= Grid->JCl && jj <= Grid->JCu ){
      // Analyze interior cells

      // Set correlation between face and block side
      IsTheSameBlockEdge = 1;

      // === Set the face types of this interior cell

      if (ii == Grid->ICl  && Grid->BndWestSplineInfo != NULL){
	// the West face is curved
	CellFacesInfo[0] = true;
	AtLeastOneCurvedFace = true;
      } else {
	// the West face is treated as straight edge
	CellFacesInfo[0] = false;
      }

      if (jj == Grid->JCl  && Grid->BndSouthSplineInfo != NULL){
	// the South face is curved
	CellFacesInfo[1] = true;
	AtLeastOneCurvedFace = true;
      } else {
	// the South face is treated as straight edge
	CellFacesInfo[1] = false;      
      }

      if (ii == Grid->ICu  && Grid->BndEastSplineInfo != NULL){
	// the East face is curved
	CellFacesInfo[2] = true;
	AtLeastOneCurvedFace = true;
      } else {
	// the East face is treated as straight edge
	CellFacesInfo[2] = false;
      }

      if (jj == Grid->JCu  && Grid->BndNorthSplineInfo != NULL){
	// the North face is curved
	CellFacesInfo[3] = true;
	AtLeastOneCurvedFace = true;      
      } else {
	// the North face is treated as straight edge
	CellFacesInfo[3] = false;
      }

    } else {
      // Analyze ghost cells
      // Attention: For a ghost cell the east face is analyzed with the west boundary and vice-versa.
      //            Same thing for the south and north faces.

      // Set correlation between face and block side
      IsTheSameBlockEdge = -1;

      // === Set the face types of this ghost cell

      if (ii == Grid->ICu+1 && jj >= Grid->JCl && jj <= Grid->JCu && Grid->BndEastSplineInfo != NULL){
	// the West face is curved
	CellFacesInfo[0] = true;
	AtLeastOneCurvedFace = true;
      } else {
	// the West face is treated as straight edge
	CellFacesInfo[0] = false;
      }

      if (jj == Grid->JCu+1 && ii >= Grid->ICl && ii <= Grid->ICu && Grid->BndNorthSplineInfo != NULL){
	// the South face is curved
	CellFacesInfo[1] = true;
	AtLeastOneCurvedFace = true;
      } else {
	// the South face is treated as straight edge
	CellFacesInfo[1] = false;
      }

      if (ii == Grid->ICl-1 && jj >= Grid->JCl && jj <= Grid->JCu && Grid->BndWestSplineInfo != NULL){
	// the East face is curved
	CellFacesInfo[2] = true;
	AtLeastOneCurvedFace = true;
      } else {
	// the East face is treated as straight edge
	CellFacesInfo[2] = false;
      }

      if (jj == Grid->JCl-1 && ii >= Grid->ICl && ii <= Grid->ICu && Grid->BndSouthSplineInfo != NULL){
	// the North face is curved
	CellFacesInfo[3] = true;
	AtLeastOneCurvedFace = true;
      } else {
	// the North face is treated as straight edge
	CellFacesInfo[3] = false;
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

#endif

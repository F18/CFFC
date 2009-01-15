/*!\file HO_Grid2DQuad_Specializations.h
  \brief Header file defining specializations for the 2D quadrilateral mesh class Grid2D_Quad_Block_HO.
  \note To use the specializations defined in this file include this file at the end of 'HO_Grid2DQuad.h'!
*/

#ifndef _HO_GRID2D_QUAD_BLOCK_SPECIALIZATIONS_INCLUDED 
#define _HO_GRID2D_QUAD_BLOCK_SPECIALIZATIONS_INCLUDED 

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
// None


/*!
 * Specialization of IntegrateFunctionOverCellUsingContourIntegrand()
 * from Grid2DQuadIntegration class.
 * Integrate a general function (i.e. any function or pointer function)
 * over the domain of a cell (ii,jj) using a contour integration.
 * The contour integrand of the function is used for this purpose.
 *
 * \param ii the i-index of the cell over which the integration is performed
 * \param jj the j-index of the cell over which the integration is performed
 * \param FuncObj the function to be integrated
 * \param ContourIntegrand the integrand with respect to x of FuncObj
 * \param digits the number of exact digits with which the result is computed (i.e. the accuracy of the calculation)
 * \param _dummy_param a parameter used only to determine the return type of the function FuncObj
 *
 * The dependent variable for the contour integration is the y-coordinate.
 */
template<>
template<typename FO, class ReturnType> inline
ReturnType Grid2DQuadIntegration<Grid2D_Quad_Block_HO>::
IntegrateFunctionOverCellUsingContourIntegrand(const int &ii, const int &jj,
					       const FO FuncObj, FO ContourIntegrand,
					       int digits, ReturnType _dummy_param) const {

  ostringstream ostm;

  // Check if the extension splines are not defined as curved.
  // This routine is not implemented to integrate correctly the function over ghost cells in these cases.
  if (Grid->IsNorthExtendWestBoundaryCurved() || 
      Grid->IsSouthExtendWestBoundaryCurved() ||
      Grid->IsNorthExtendEastBoundaryCurved() ||
      Grid->IsSouthExtendEastBoundaryCurved() ||
      Grid->IsEastExtendSouthBoundaryCurved() ||
      Grid->IsWestExtendSouthBoundaryCurved() ||
      Grid->IsEastExtendNorthBoundaryCurved() ||
      Grid->IsWestExtendNorthBoundaryCurved() ){

    // Build error message
    ostm << "Grid2DQuadIntegration<Grid2D_Quad_Block_HO>::IntegrateFunctionOverCellUsingContourIntegrand() ERROR!\n"
	 << "This routine doesn't know how to handle correctly integration over domains with curved extension splines.\n"
	 << "Options: Add this logic or use Monte Carlo integration instead of Gauss quadrature integration\n";

    throw runtime_error(ostm.str());
  }

  // SET VARIABLES USED IN THE INTEGRATION PROCESS

  ReturnType IntResult(0.0);		// the integration result
  double DeltaY;
  int NumGQP(Spline2DInterval_HO::get_NumGQPoints_ContourIntegral());   // Get the number of GQPs per each straight edge
  Vector2D *GaussQuadPoints = new Vector2D [NumGQP];   // the GQPs used to calculate the integral along a line segment
  double * GaussQuadWeights = new double [NumGQP];     // the Gauss integration weights for each Gauss quadrature
  int n;

  /* Set the GaussQuadWeights. */
  GaussQuadratureData::getGaussQuadWeights(GaussQuadWeights, NumGQP);

  // Integrate along West face
  if ( CellFacesInfo[0] ){
    // == This face is curved ==
    if ( getFaceBlockEdgeCorrelation() ){
      // this is an interior cell near West boundary
      IntResult +=  Grid->BndWestSplineInfo[jj].IntegrateFunctionWithRespectToY(ContourIntegrand,IntResult);
    } else {
      // this is a ghost cell near East boundary
      IntResult -=  Grid->BndEastSplineInfo[jj].IntegrateFunctionWithRespectToY(ContourIntegrand,IntResult);
    }
  } else {
    // == This face is straight ==

    // get the Gauss quadrature points
    Grid->getGaussQuadPointsFaceW(ii,jj,GaussQuadPoints,NumGQP);

    // get the derivative of the y-coordinate with respect to the path coordinate
    DeltaY = (Grid->nodeSW(ii,jj).y() - Grid->nodeNW(ii,jj).y());

    // integrate along the line segment
    IntResult += CalculateFunctionIntegralWithGaussQuadratures(ContourIntegrand, GaussQuadPoints, GaussQuadWeights,
							       NumGQP, DeltaY, IntResult);
  }

  // Integrate along South boundary
  if ( CellFacesInfo[1] ){
    // == This face is curved ==
    if ( getFaceBlockEdgeCorrelation() ){
      // this is an interior cell near South boundary
      IntResult +=  Grid->BndSouthSplineInfo[ii].IntegrateFunctionWithRespectToY(ContourIntegrand,IntResult);
    } else {
      // this is a ghost cell near North boundary
      IntResult -=  Grid->BndNorthSplineInfo[ii].IntegrateFunctionWithRespectToY(ContourIntegrand,IntResult);
    }
  } else {
    // == This face is straight ==

    // get the Gauss quadrature points
    Grid->getGaussQuadPointsFaceS(ii,jj,GaussQuadPoints,NumGQP);

    // get the derivative of the y-coordinate with respect to the path coordinate
    DeltaY = Grid->nodeSE(ii,jj).y() - Grid->nodeSW(ii,jj).y();

    // integrate along the line segment
    IntResult += CalculateFunctionIntegralWithGaussQuadratures(ContourIntegrand, GaussQuadPoints, GaussQuadWeights,
							       NumGQP, DeltaY, IntResult);
  }

  // Integrate along East boundary
  if ( CellFacesInfo[2] ){
    // == This face is curved ==
    if ( getFaceBlockEdgeCorrelation() ){
      // this is an interior cell near East boundary
      IntResult +=  Grid->BndEastSplineInfo[jj].IntegrateFunctionWithRespectToY(ContourIntegrand,IntResult);
    } else {
      // this is a ghost cell near West boundary
      IntResult -=  Grid->BndWestSplineInfo[jj].IntegrateFunctionWithRespectToY(ContourIntegrand,IntResult);
    }
  } else {
    // == This face is straight ==

    // get the Gauss quadrature points
    Grid->getGaussQuadPointsFaceE(ii,jj,GaussQuadPoints,NumGQP);

    // get the derivative of the y-coordinate with respect to the path coordinate
    DeltaY = Grid->nodeNE(ii,jj).y() - Grid->nodeSE(ii,jj).y();

    // integrate along the line segment
    IntResult += CalculateFunctionIntegralWithGaussQuadratures(ContourIntegrand, GaussQuadPoints, GaussQuadWeights,
							       NumGQP, DeltaY, IntResult);
  }

  // Integrate along North boundary
  if ( CellFacesInfo[3] ){
    // == This face is curved ==
    if ( getFaceBlockEdgeCorrelation() ){
      // this is an interior cell near North boundary
      IntResult +=  Grid->BndNorthSplineInfo[ii].IntegrateFunctionWithRespectToY(ContourIntegrand,IntResult);
    } else {
      // this is a ghost cell near South boundary
      IntResult -=  Grid->BndSouthSplineInfo[ii].IntegrateFunctionWithRespectToY(ContourIntegrand,IntResult);
    }
  } else {
    // == This face is straight ==

    // get the Gauss quadrature points
    Grid->getGaussQuadPointsFaceN(ii,jj,GaussQuadPoints,NumGQP);

    // get the derivative of the y-coordinate with respect to the path coordinate
    DeltaY = Grid->nodeNW(ii,jj).y() - Grid->nodeNE(ii,jj).y();

    // integrate along the line segment
    IntResult += CalculateFunctionIntegralWithGaussQuadratures(ContourIntegrand, GaussQuadPoints, GaussQuadWeights,
							       NumGQP, DeltaY, IntResult);
  }

  // Deallocate memory
  delete [] GaussQuadPoints;
  delete [] GaussQuadWeights;

  return IntResult;
}


#endif

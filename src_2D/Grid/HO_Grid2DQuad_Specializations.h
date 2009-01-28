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

/*!
 * Approximate the curvilinear domain of cell (ii,jj) with a polygon.
 */
template<> inline
void Grid2DQuadIntegration<Grid2D_Quad_Block_HO>::ConvertCurvedIntegrationDomainToPolygon(const int &ii, const int &jj,
											  Polygon & PolygonalDomain) const {

  LinkedList<Vector2D> PolygonalDomainVertexes;

  // Build the polygonal integration domain. 
  // Get the vertexes in counterclockwise direction.
  // Each cell edge adds its own start point to avoid duplication.
  
  // === Check West face ===
  if (CellBndWest != NULL){
    // The cell is bounded by a curved spline to the West
    CellBndWest->GenerateLineSegmentDiscretization(Grid->Node[ii][jj+1].X, Grid->Node[ii][jj].X,
						   PolygonalDomainVertexes);
    // Remove the last element to avoid duplication
    PolygonalDomainVertexes.removeLast();
  } else {
    // The cell has a straight edge to the West
    PolygonalDomainVertexes.add(Grid->Node[ii][jj+1].X);
  }

  // === Check South face ===
  if (CellBndSouth != NULL){
    // The cell is bounded by a curved spline to the South
    CellBndSouth->GenerateLineSegmentDiscretization(Grid->Node[ii][jj].X, Grid->Node[ii+1][jj].X,
						    PolygonalDomainVertexes);
    // Remove the last element to avoid duplication
    PolygonalDomainVertexes.removeLast();    
  } else {
    // The cell has a straight edge to the South
    PolygonalDomainVertexes.add(Grid->Node[ii][jj].X);
  }

  // === Check East face ===
  if (CellBndEast != NULL){
    // The cell is bounded by a curved spline to the East
    CellBndEast->GenerateLineSegmentDiscretization(Grid->Node[ii+1][jj].X, Grid->Node[ii+1][jj+1].X,
						   PolygonalDomainVertexes);
    // Remove the last element to avoid duplication
    PolygonalDomainVertexes.removeLast();    
  } else {
    // The cell has a straight edge to the East
    PolygonalDomainVertexes.add(Grid->Node[ii+1][jj].X);
  }

  // === Check North face ===
  if (CellBndNorth != NULL){
    // The cell is bounded by a curved spline to the North
    CellBndNorth->GenerateLineSegmentDiscretization(Grid->Node[ii+1][jj+1].X, Grid->Node[ii][jj+1].X,
						    PolygonalDomainVertexes);
    // Remove the last element to avoid duplication
    PolygonalDomainVertexes.removeLast();    
  } else {
    // The cell has a straight edge to the North
    PolygonalDomainVertexes.add(Grid->Node[ii+1][jj+1].X);
  }
  
  // Generate the polygonal integration domain which approximates the curvilinear one
  PolygonalDomain.convert(PolygonalDomainVertexes);

}


/*!
 * Specialization of IntegrateFunctionOverCellUsingMonteCarloMethod()
 * from Grid2DQuadIntegration class.
 * Integrate a general function (i.e. any function or pointer function)
 * over the domain of a cell (ii,jj) using a Monte Carlo integration method.
 *
 * \param ii the i-index of the cell over which the integration is performed
 * \param jj the j-index of the cell over which the integration is performed
 * \param FuncObj the function to be integrated
 * \param digits the number of exact digits with which the result is computed (i.e. the accuracy of the calculation)
 * \param _dummy_param a parameter used only to determine the return type of the function FuncObj
 *
 * /note Depending on the The exact digits might not be obtained.
 * /note The routine AnalyseCellFaces(ii,jj) must have been called before this routine.
 */
template<>
template<typename FO, class ReturnType>
ReturnType Grid2DQuadIntegration<Grid2D_Quad_Block_HO>::
IntegrateFunctionOverCellUsingMonteCarloMethod(const int &ii, const int &jj,
					       FO FuncObj,
					       int digits, ReturnType _dummy_param) const {

  // SET VARIABLES USED IN THE INTEGRATION PROCESS
  ReturnType SumFunc(0.0);
  Vector2D Xlo, Xhi;
  double x,y;  
  Polygon DefinitionDomain;
  RandomGen RandGen(17);                  // Initialize random number generator
  double DeltaX, DeltaY, IntegrationArea;
  int NumSamples(NumericalLibrary_Execution_Mode::Number_Monte_Carlo_Samples);
  int UsefulNumSamples(0);

  // Determine the definition domain of the integrand
  ConvertCurvedIntegrationDomainToPolygon(ii, jj, DefinitionDomain);

  // Determine the polygon bounding box
  DefinitionDomain.BoundingBoxCoordinates(Xlo, Xhi);
  DefinitionDomain.area();

  // Determine integration domain that includes the definition domain
  DeltaX = Xhi.x - Xlo.x;
  DeltaY = Xhi.y - Xlo.y;
  //  IntegrationArea = DeltaX * DeltaY;
  IntegrationArea = DefinitionDomain.A;

  // Apply Monte Carlo method
  for (int i=1; i<=NumSamples; ++i){

    // Get the coordinates of a random point in the DefinitionDomain bounding box
    x = Xlo.x + RandGen.doub()*DeltaX;
    y = Xlo.y + RandGen.doub()*DeltaY;
    
    if (DefinitionDomain.IsPointInPolygon(Vector2D(x,y))){
      // This point is inside the DefinitionDomain
      SumFunc += FuncObj(x,y); 
      ++UsefulNumSamples;
    }
  }

  return (IntegrationArea*SumFunc)/UsefulNumSamples;

}

/*!
 * Specialization of IntegrateFunctionOverCellUsingPolygonalAdaptiveQuadratures()
 * from Grid2DQuadIntegration class.
 * Integrate a general function (i.e. any function or pointer function)
 * over the domain of a cell (ii,jj) using an adaptive Gaussian quadrature integration method.
 *
 * \param ii the i-index of the cell over which the integration is performed
 * \param jj the j-index of the cell over which the integration is performed
 * \param FuncObj the function to be integrated
 * \param digits the number of exact digits with which the result is computed (i.e. the accuracy of the calculation)
 * \param _dummy_param a parameter used only to determine the return type of the function FuncObj
 *
 * /note The routine AnalyseCellFaces(ii,jj) must have been called before this routine.
 */
template<>
template<typename FO, class ReturnType>
ReturnType Grid2DQuadIntegration<Grid2D_Quad_Block_HO>::
IntegrateFunctionOverCellUsingPolygonalAdaptiveQuadratures(const int &ii, const int &jj,
							   FO FuncObj,
							   int digits, ReturnType _dummy_param) const {

  // SET VARIABLES USED IN THE INTEGRATION PROCESS
  ReturnType IntQuad;
  Vector2D Xlo, Xhi;
  Polygon DefinitionDomain;
  int BackUpMinRefLevels(NumericalLibrary_Execution_Mode::Adaptive_Integration_Minimum_Refinement_Levels);
  

  // Determine the definition domain of the integrand
  ConvertCurvedIntegrationDomainToPolygon(ii, jj, DefinitionDomain);

  // Determine the polygon bounding box
  DefinitionDomain.BoundingBoxCoordinates(Xlo, Xhi);

  // Build a local discontinuous integrand (i.e. outside of the DefinitionDomain it returns zero)
  IntegrandFunctionOverPolygonalDefinitionDomain<FO,ReturnType> LocalFO(FuncObj,
  									DefinitionDomain);

  // Ensure that a more suitable minimum number of adaptive refinement integration levels is required
  NumericalLibrary_Execution_Mode::
    Adaptive_Integration_Minimum_Refinement_Levels = Grid->Polygonal_Adaptive_Quadrature_Integration_Minimum_Levels;

  // Calculate integral
  IntQuad = QuadrilateralQuadrature(LocalFO,
				    typename Grid2D_Quad_Block_HO::NodeType(Xlo),
				    typename Grid2D_Quad_Block_HO::NodeType(Xlo.x, Xhi.y),
				    typename Grid2D_Quad_Block_HO::NodeType(Xhi),
				    typename Grid2D_Quad_Block_HO::NodeType(Xhi.x, Xlo.y),
				    digits,
				    _dummy_param);

  // Set back the minimum number of adaptive refinement integration levels
  NumericalLibrary_Execution_Mode::Adaptive_Integration_Minimum_Refinement_Levels = BackUpMinRefLevels;

  return IntQuad;
}

#endif

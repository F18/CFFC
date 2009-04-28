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
 * Compute the average solution for each sector of cell (ii,jj)
 * (i.e. South-West, South-East, North-West, North-East) resulting
 * from the integration of the provided function. 
 * These sectors are defined by the mesh refinement routine.
 *
 * \param [in] ii i-index of the cell which provides the geometry
 * \param [in] jj j-index of the cell which provides the geometry
 * \param [in] ContourIntegrand the function used for contour integration. 
 *             It is assumed to have been already integrated with respect to x.
 * \param [in] CoarseAvgSoln the average solution in (ii,jj) cell (precomputed already)
 * \param [out] AvgSoln_SW the average solution corresponding to the SW sector
 * \param [out] AvgSoln_NW the average solution corresponding to the NW sector
 * \param [out] AvgSoln_SE the average solution corresponding to the SE sector
 * \param [out] AvgSoln_NE the average solution corresponding to the NE sector
 *
 */
template<>
template<typename FO, class ReturnType> inline
void Grid2DQuadIntegration<Grid2D_Quad_Block_HO>::
IntegrateFunctionOverCellSectorsUsingContourIntegrand(const int &ii, const int &jj,
						      FO ContourIntegrand,
						      const ReturnType & CoarseAvgSoln,
						      ReturnType & AvgSoln_SW, ReturnType & AvgSoln_NW,
						      ReturnType & AvgSoln_SE, ReturnType & AvgSoln_NE) {

  Grid2D_Quad_Block_HO::NodeType MidN,MidS,MidE,MidW,CC; //< Nodes of the sectors in addition to the those of the current cell

  // SET VARIABLES USED IN THE INTEGRATION PROCESS

  double DeltaY;
  int NumGQP(3);   // Set the number of GQPs per each straight edge. Three GQPs are enough for polynomials up to quartic!
  Vector2D *GaussQuadPoints = new Vector2D [NumGQP];   // the GQPs used to calculate the integral along a line segment
  double * GaussQuadWeights = new double [NumGQP];     // the Gauss integration weights for each Gauss quadrature
  Spline2D_HO SplineCopy;

  /* Set the GaussQuadWeights. */
  GaussQuadratureData::getGaussQuadWeights(GaussQuadWeights, NumGQP);

  // Create spline interval variables in case there are curved faces (only those curved faces will use these variables)
  typename Grid2D_Quad_Block_HO::BndSplineIntervalType SplineIntervalInfo_1, SplineIntervalInfo_2;
  
  // Set storage for the contour integral along each face
  ReturnType faceN_NW_Int, faceN_NE_Int,
    faceW_SW_Int, faceW_NW_Int,
    faceS_SW_Int, faceS_SE_Int,
    faceE_SE_Int, faceE_NE_Int,
    faceW_CC_Int, faceS_CC_Int,
    faceE_CC_Int, faceN_CC_Int;
  double faceN_NW_area, faceN_NE_area,
    faceW_SW_area, faceW_NW_area,
    faceS_SW_area, faceS_SE_area,
    faceE_SE_area, faceE_NE_area,
    faceW_CC_area, faceS_CC_area,
    faceE_CC_area, faceN_CC_area,
    area_SW, area_NW, area_NE, area_SE;


  // === Analyse cell faces
  AnalyseCellFaces(ii,jj);

  // === Determine the nodes of the sectors in addition to the nodes of the current cell.
  MidW.setloc(Grid->getMidPointFaceW(ii,jj));
  MidN.setloc(Grid->getMidPointFaceN(ii,jj));
  MidE.setloc(Grid->getMidPointFaceE(ii,jj));
  MidS.setloc(Grid->getMidPointFaceS(ii,jj));
  CC.setloc(Grid->getNodeAverage(ii,jj));

  /* === Calculate the integral along each integration segment. 
     The arrow indicates the orientation considered for each segment.
     Special care must be paid to the sign of the interior integrals when they are summed up! */
  
  // Integrate Along Line Segment CC->MidW
  Grid->getLineSegmentGaussIntegrationData(CC.X, MidW.X, GaussQuadPoints, NumGQP, DeltaY);
  faceW_CC_Int = CalculateFunctionIntegralWithGaussQuadratures(ContourIntegrand, GaussQuadPoints, GaussQuadWeights,
							       NumGQP, DeltaY, faceW_CC_Int);
  faceW_CC_area = ZeroLineIntegration(CC, MidW);
  // Integrate Along Line Segment CC->MidS
  Grid->getLineSegmentGaussIntegrationData(CC.X, MidS.X, GaussQuadPoints, NumGQP, DeltaY);
  faceS_CC_Int = CalculateFunctionIntegralWithGaussQuadratures(ContourIntegrand, GaussQuadPoints, GaussQuadWeights,
							       NumGQP, DeltaY, faceS_CC_Int);
  faceS_CC_area = ZeroLineIntegration(CC, MidS);
  // Integrate Along Line Segment CC->MidE
  Grid->getLineSegmentGaussIntegrationData(CC.X, MidE.X, GaussQuadPoints, NumGQP, DeltaY);
  faceE_CC_Int = CalculateFunctionIntegralWithGaussQuadratures(ContourIntegrand, GaussQuadPoints, GaussQuadWeights,
							       NumGQP, DeltaY, faceE_CC_Int);
  faceE_CC_area = ZeroLineIntegration(CC, MidE);
  // Integrate Along Line Segment CC->MidN
  Grid->getLineSegmentGaussIntegrationData(CC.X, MidN.X, GaussQuadPoints, NumGQP, DeltaY);
  faceN_CC_Int = CalculateFunctionIntegralWithGaussQuadratures(ContourIntegrand, GaussQuadPoints, GaussQuadWeights,
							       NumGQP, DeltaY, faceN_CC_Int);
  faceN_CC_area = ZeroLineIntegration(CC, MidN);

  // Integrate along West boundary
  if ( CellFacesInfo[0] ){
    // == This face is curved ==

    // Check if the West Spline is defined such that the normals at the GaussQuadratures point outside of the domain 
    // (i.e The spline pathlength increases from JNu to JNl)
    if ( CellBndWest->getS(Grid->nodeSW(ii,jj)) > CellBndWest->getS(Grid->nodeNW(ii,jj)) ){
      // Generate spline intervals. The number of Gauss points is set by the NUMBER_OF_GQP_CONTOURINT variable and not by the value 1!
      SplineIntervalInfo_1.InitializeInterval(*CellBndWest,MidW,Grid->nodeNW(ii,jj),1);
      SplineIntervalInfo_2.InitializeInterval(*CellBndWest,Grid->nodeSW(ii,jj),MidW,1);
    } else {
      // Copy the spline
      SplineCopy = *CellBndWest;
      // Change the direction of increasing the pathlength
      SplineCopy.Reverse_Spline();
      // Generate spline intervals
      SplineIntervalInfo_1.InitializeInterval(SplineCopy,MidW,Grid->nodeNW(ii,jj),1);
      SplineIntervalInfo_2.InitializeInterval(SplineCopy,Grid->nodeSW(ii,jj),MidW,1);
    }//endif

    // Integrate Along Line Segment nodeNW->MidW
    faceW_NW_Int  = SplineIntervalInfo_1.IntegrateFunctionWithRespectToY(ContourIntegrand,faceW_NW_Int);
    faceW_NW_area = SplineIntervalInfo_1.AreaContribution();
    
    // Integrate Along Line Segment MidW->nodeSW
    faceW_SW_Int  = SplineIntervalInfo_2.IntegrateFunctionWithRespectToY(ContourIntegrand,faceW_SW_Int);
    faceW_SW_area = SplineIntervalInfo_2.AreaContribution();

  } else {
    // == This face is straight ==

    // Integrate Along Line Segment nodeNW->MidW
    Grid->getLineSegmentGaussIntegrationData(Grid->nodeNW(ii,jj).X, MidW.X, GaussQuadPoints, NumGQP, DeltaY);
    faceW_NW_Int = CalculateFunctionIntegralWithGaussQuadratures(ContourIntegrand, GaussQuadPoints, GaussQuadWeights,
								 NumGQP, DeltaY, faceW_NW_Int);
    faceW_NW_area = ZeroLineIntegration(Grid->nodeNW(ii,jj), MidW);
    
    // Integrate Along Line Segment MidW->nodeSW
    Grid->getLineSegmentGaussIntegrationData(MidW.X, Grid->nodeSW(ii,jj).X, GaussQuadPoints, NumGQP, DeltaY);
    faceW_SW_Int = CalculateFunctionIntegralWithGaussQuadratures(ContourIntegrand, GaussQuadPoints, GaussQuadWeights,
								 NumGQP, DeltaY, faceW_SW_Int);
    faceW_SW_area = ZeroLineIntegration(MidW, Grid->nodeSW(ii,jj));
  }

  // Integrate along South boundary
  if ( CellFacesInfo[1] ){
    // == This face is curved ==

    // Check if the South Spline is defined such that the normals at the GaussQuadratures point outside of the domain 
    // (i.e The spline pathlength increases from INl to INu)
    if ( CellBndSouth->getS(Grid->nodeSW(ii,jj)) < CellBndSouth->getS(Grid->nodeSE(ii,jj)) ){
      // Generate spline intervals
      SplineIntervalInfo_1.InitializeInterval(*CellBndSouth,Grid->nodeSW(ii,jj),MidS,1);
      SplineIntervalInfo_2.InitializeInterval(*CellBndSouth,MidS,Grid->nodeSE(ii,jj),1);
    } else {
      // Copy the spline
      SplineCopy = *CellBndSouth;
      // Change the direction of increasing the pathlength
      SplineCopy.Reverse_Spline();
      // Generate spline intervals
      SplineIntervalInfo_1.InitializeInterval(SplineCopy,Grid->nodeSW(ii,jj),MidS,1);
      SplineIntervalInfo_2.InitializeInterval(SplineCopy,MidS,Grid->nodeSE(ii,jj),1);
    }//endif

    // Integrate Along Line Segment nodeSW->MidS
    faceS_SW_Int  = SplineIntervalInfo_1.IntegrateFunctionWithRespectToY(ContourIntegrand,faceS_SW_Int);
    faceS_SW_area = SplineIntervalInfo_1.AreaContribution();
    
    // Integrate Along Line Segment MidS->nodeSE (Currently not used. It's substituted by the conservation law!!!)
    //  faceS_SE_Int  = SplineIntervalInfo_2.IntegrateFunctionWithRespectToY(ContourIntegrand,faceS_SE_Int);
    faceS_SE_area = SplineIntervalInfo_2.AreaContribution();

  } else {
    // == This face is straight ==

    // Integrate Along Line Segment nodeSW->MidS
    Grid->getLineSegmentGaussIntegrationData(Grid->nodeSW(ii,jj).X, MidS.X, GaussQuadPoints, NumGQP, DeltaY);
    faceS_SW_Int = CalculateFunctionIntegralWithGaussQuadratures(ContourIntegrand, GaussQuadPoints, GaussQuadWeights,
								 NumGQP, DeltaY, faceS_SW_Int);
    faceS_SW_area = ZeroLineIntegration(Grid->nodeSW(ii,jj), MidS);
    
    // Integrate Along Line Segment MidS->nodeSE (Currently not used. It's substituted by the conservation law!!!)
    //  Grid->getLineSegmentGaussIntegrationData(MidS.X, Grid->nodeSE(ii,jj).X, GaussQuadPoints, NumGQP, DeltaY);
    //  faceS_SE_Int = CalculateFunctionIntegralWithGaussQuadratures(ContourIntegrand, GaussQuadPoints, GaussQuadWeights,
    //  								 NumGQP, DeltaY, faceS_SE_Int);
    faceS_SE_area = ZeroLineIntegration(MidS, Grid->nodeSE(ii,jj));
  }

  // Integrate along East boundary
  if ( CellFacesInfo[2] ){
    // == This face is curved ==

    // Check if the East Spline is defined such that the normals at the GaussQuadratures point outside of the domain 
    // (i.e The spline pathlength increases from JNl to JNu)
    if ( CellBndEast->getS(Grid->nodeSE(ii,jj)) < CellBndEast->getS(Grid->nodeNE(ii,jj)) ){
      // Generate spline intervals
      SplineIntervalInfo_1.InitializeInterval(*CellBndEast,Grid->nodeSE(ii,jj),MidE,1);
      SplineIntervalInfo_2.InitializeInterval(*CellBndEast,MidE,Grid->nodeNE(ii,jj),1);
    } else {
      // Copy the spline
      SplineCopy = *CellBndEast;
      // Change the direction of increasing the pathlength
      SplineCopy.Reverse_Spline();
      // Generate spline intervals
      SplineIntervalInfo_1.InitializeInterval(SplineCopy,Grid->nodeSE(ii,jj),MidE,1);
      SplineIntervalInfo_2.InitializeInterval(SplineCopy,MidE,Grid->nodeNE(ii,jj),1);
    }//endif

    // Integrate Along Line Segment nodeSE->MidE (Currently not used. It's substituted by the conservation law!!!)
    //  faceE_SE_Int  = SplineIntervalInfo_1.IntegrateFunctionWithRespectToY(ContourIntegrand,faceE_SE_Int); 
    faceE_SE_area = SplineIntervalInfo_1.AreaContribution();

    // Integrate Along Line Segment MidE->nodeNE
    faceE_NE_Int  = SplineIntervalInfo_2.IntegrateFunctionWithRespectToY(ContourIntegrand,faceE_NE_Int);
    faceE_NE_area = SplineIntervalInfo_2.AreaContribution();

  } else {
    // == This face is straight ==

    // Integrate Along Line Segment nodeSE->MidE (Currently not used. It's substituted by the conservation law!!!)
    //  Grid->getLineSegmentGaussIntegrationData(Grid->nodeSE(ii,jj).X, MidE.X, GaussQuadPoints, NumGQP, DeltaY);
    //  faceE_SE_Int = CalculateFunctionIntegralWithGaussQuadratures(ContourIntegrand, GaussQuadPoints, GaussQuadWeights,
    //  								 NumGQP, DeltaY, faceE_SE_Int);
    faceE_SE_area = ZeroLineIntegration(Grid->nodeSE(ii,jj), MidE);

    // Integrate Along Line Segment MidE->nodeNE
    Grid->getLineSegmentGaussIntegrationData(MidE.X, Grid->nodeNE(ii,jj).X, GaussQuadPoints, NumGQP, DeltaY);
    faceE_NE_Int = CalculateFunctionIntegralWithGaussQuadratures(ContourIntegrand, GaussQuadPoints, GaussQuadWeights,
								 NumGQP, DeltaY, faceE_NE_Int);
    faceE_NE_area = ZeroLineIntegration(MidE, Grid->nodeNE(ii,jj));
  }

  // Integrate along North boundary
  if ( CellFacesInfo[3] ){
    // == This face is curved ==

    // Check if the North Spline is defined such that the normals at the GaussQuadratures point outside of the domain 
    // (i.e The spline pathlength increases from INu to INl)
    if ( CellBndNorth->getS(Grid->nodeNW(ii,jj)) > CellBndNorth->getS(Grid->nodeNE(ii,jj)) ){
      // Generate spline intervals
      SplineIntervalInfo_1.InitializeInterval(*CellBndNorth,Grid->nodeNW(ii,jj),MidN,1);
      SplineIntervalInfo_2.InitializeInterval(*CellBndNorth,MidN,Grid->nodeNE(ii,jj),1);
    } else {
      // Copy the spline
      SplineCopy = *CellBndNorth;
      // Change the direction of increasing the pathlength
      SplineCopy.Reverse_Spline();
      // Generate spline intervals
      SplineIntervalInfo_1.InitializeInterval(SplineCopy,Grid->nodeNW(ii,jj),MidN,1);
      SplineIntervalInfo_2.InitializeInterval(SplineCopy,MidN,Grid->nodeNE(ii,jj),1);
    }//endif

    // Integrate Along Line Segment nodeNE->MidN
    faceN_NE_Int  = SplineIntervalInfo_2.IntegrateFunctionWithRespectToY(ContourIntegrand,faceN_NE_Int);
    faceN_NE_area = SplineIntervalInfo_2.AreaContribution();

    // Integrate Along Line Segment MidN->nodeNW
    faceN_NW_Int  = SplineIntervalInfo_1.IntegrateFunctionWithRespectToY(ContourIntegrand,faceN_NW_Int);
    faceN_NW_area = SplineIntervalInfo_1.AreaContribution();

  } else {
    // == This face is straight ==

    // Integrate Along Line Segment nodeNE->MidN
    Grid->getLineSegmentGaussIntegrationData(Grid->nodeNE(ii,jj).X, MidN.X, GaussQuadPoints, NumGQP, DeltaY);
    faceN_NE_Int = CalculateFunctionIntegralWithGaussQuadratures(ContourIntegrand, GaussQuadPoints, GaussQuadWeights,
								 NumGQP, DeltaY, faceN_NE_Int);
    faceN_NE_area = ZeroLineIntegration(Grid->nodeNE(ii,jj), MidN);

    // Integrate Along Line Segment MidN->nodeNW
    Grid->getLineSegmentGaussIntegrationData(MidN.X, Grid->nodeNW(ii,jj).X, GaussQuadPoints, NumGQP, DeltaY);
    faceN_NW_Int = CalculateFunctionIntegralWithGaussQuadratures(ContourIntegrand, GaussQuadPoints, GaussQuadWeights,
								 NumGQP, DeltaY, faceN_NW_Int);
    faceN_NW_area = ZeroLineIntegration(MidN, Grid->nodeNW(ii,jj));
  }

  // === Calculate area for each cell sector
  area_SW = faceW_CC_area + faceW_SW_area + faceS_SW_area - faceS_CC_area;
  area_SE = faceS_CC_area + faceS_SE_area + faceE_SE_area - faceE_CC_area;
  area_NE = faceE_CC_area + faceE_NE_area + faceN_NE_area - faceN_CC_area;
  area_NW = faceN_CC_area + faceN_NW_area + faceW_NW_area - faceW_CC_area;
  
  // === Calculate integrals for each cell sector
  AvgSoln_SW  = (faceW_CC_Int + faceW_SW_Int + faceS_SW_Int - faceS_CC_Int);
  AvgSoln_NE  = (faceE_CC_Int + faceE_NE_Int + faceN_NE_Int - faceN_CC_Int);
  AvgSoln_NW  = (faceN_CC_Int + faceN_NW_Int + faceW_NW_Int - faceW_CC_Int);  
  
  // === Calculate final average solutions

  /* Calculate AvgSoln_SE based on conservation equation instead of 
     AvgSoln_SE  = (faceS_CC_Int + faceS_SE_Int + faceE_SE_Int - faceE_CC_Int)/area_SE;
     This approach is conservative and will also be accurate if high-order treatment of boundaries is performed.
     There is no particular reason why SE sector has been choosen instead of a different cell sector.
     A better approach is to use the sector with the most curved boundaries at the expense of a more complicated algorithm.
  */
  AvgSoln_SE = (CoarseAvgSoln*Grid->CellArea(ii,jj) - (AvgSoln_SW + AvgSoln_NE + AvgSoln_NW) )/area_SE;
  AvgSoln_SW /= area_SW;
  AvgSoln_NE /= area_NE;
  AvgSoln_NW /= area_NW;

  // Deallocate memory
  delete [] GaussQuadPoints;
  delete [] GaussQuadWeights;
  
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

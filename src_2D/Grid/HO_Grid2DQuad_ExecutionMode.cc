/*!\file HO_Grid2DQuad_ExecutionMode.cc
  \brief Initialize the flags that control the execution of high-order 2D grid. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "HO_Grid2DQuad_ExecutionMode.h"
#include "../MPI/MPI.h"

short HO_Grid2D_Execution_Mode::USE_HIGH_ORDER_GEOMETRIC_BOUNDARY_REPRESENTATION = OFF; // High-order boundaries are not used
short HO_Grid2D_Execution_Mode::METHOD_TO_INTEGRATE_ALONG_CURVED_EDGES = MIXED_QUADS_AND_LINES_METHOD; // Mixed method
short HO_Grid2D_Execution_Mode::NUMBER_OF_POINTS_FOR_GAUSS_QUADRATURE_INTEGRATION = 3; // use 3-point Gauss integration
short HO_Grid2D_Execution_Mode::EXERCISE_SPECIAL_CARE_TO_ROUNDOFF_ERRORS_FOR_THIS_MESH = OFF; // use Global Coordinate System


//! Set all flags to default values
// add all flag default values to this function
void HO_Grid2D_Execution_Mode::SetDefaults(void){

  USE_HIGH_ORDER_GEOMETRIC_BOUNDARY_REPRESENTATION = OFF; // High-order boundaries are not used
  METHOD_TO_INTEGRATE_ALONG_CURVED_EDGES = MIXED_QUADS_AND_LINES_METHOD; // Gauss quadrature integration mode
  NUMBER_OF_POINTS_FOR_GAUSS_QUADRATURE_INTEGRATION = 3; // use 3-point Gauss integration
  EXERCISE_SPECIAL_CARE_TO_ROUNDOFF_ERRORS_FOR_THIS_MESH = OFF; // use Global Coordinate System
}

//! Print the current execution mode
//  at the output stream
// \param [in] out_file the output stream
void HO_Grid2D_Execution_Mode::Print_Info(std::ostream & out_file){

  // output boundary representation mode
  if (USE_HIGH_ORDER_GEOMETRIC_BOUNDARY_REPRESENTATION == OFF){
    out_file << "\n  -> Boundary Accuracy: " << "2nd-Order";
  } else {
    out_file << "\n  -> Boundary Accuracy: " << "High-order";

    // output  mode
    if (METHOD_TO_INTEGRATE_ALONG_CURVED_EDGES == GAUSS_QUADRATURES_METHOD ){
      out_file << "\n     -> Integrate Along Curved Edges With: " << "Gauss quadrature";
    } else if (METHOD_TO_INTEGRATE_ALONG_CURVED_EDGES == ADAPTIVE_LINE_SEGMENTS_METHOD) {
      out_file << "\n     -> Integrate Along Curved Edges With: " << "Adaptive line segments";
    } else if (METHOD_TO_INTEGRATE_ALONG_CURVED_EDGES == MIXED_QUADS_AND_LINES_METHOD) {
      out_file << "\n     -> Integrate Along Curved Edges With: " << "Combination of Line Segments and Quads";
    } // endif

    // output  mode
    if (METHOD_TO_INTEGRATE_ALONG_CURVED_EDGES == GAUSS_QUADRATURES_METHOD || 
	METHOD_TO_INTEGRATE_ALONG_CURVED_EDGES == MIXED_QUADS_AND_LINES_METHOD){
      switch(NUMBER_OF_POINTS_FOR_GAUSS_QUADRATURE_INTEGRATION){
      case 3:
	out_file << "\n     -> Gauss points: " << "3";
	break;
      case 5:
	out_file << "\n     -> Gauss points: " << "5";
	break;
      }	// endswitch
    } // endif

  } // endif

  if (EXERCISE_SPECIAL_CARE_TO_ROUNDOFF_ERRORS_FOR_THIS_MESH == ON){
    out_file << "\n     -> Geometric properties: " << "Special attention to Round Offs";
  }
}

/*!
 * Broadcast the HO_Grid2D_Execution_Mode variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 */
void HO_Grid2D_Execution_Mode::Broadcast(void){
#ifdef _MPI_VERSION
  
  MPI::COMM_WORLD.Bcast(&USE_HIGH_ORDER_GEOMETRIC_BOUNDARY_REPRESENTATION,
 			1, 
 			MPI::SHORT, 0);
  MPI::COMM_WORLD.Bcast(&METHOD_TO_INTEGRATE_ALONG_CURVED_EDGES,
 			1, 
 			MPI::SHORT, 0);
  MPI::COMM_WORLD.Bcast(&NUMBER_OF_POINTS_FOR_GAUSS_QUADRATURE_INTEGRATION,
 			1, 
 			MPI::SHORT, 0);
  MPI::COMM_WORLD.Bcast(&EXERCISE_SPECIAL_CARE_TO_ROUNDOFF_ERRORS_FOR_THIS_MESH,
 			1, 
 			MPI::SHORT, 0);

  // Set the affected switches on all CPUs properly.
  Set_Affected_Switches();

#endif
}

/*!
 * Set all affected switches based on the parameters set
 * in HO_Grid2D_Execution_Mode class.
 */
void HO_Grid2D_Execution_Mode::Set_Affected_Switches(void){

  if (USE_HIGH_ORDER_GEOMETRIC_BOUNDARY_REPRESENTATION == ON){
    // Set the affected switch
    Grid2D_Quad_Block_HO::setHighOrderBoundaryRepresentation();
  } else {
    // Set the affected switch
    Grid2D_Quad_Block_HO::setLowOrderBoundaryRepresentation();
  }

  if (METHOD_TO_INTEGRATE_ALONG_CURVED_EDGES == GAUSS_QUADRATURES_METHOD){
    // Set the affected switch
    Grid2D_Quad_Block_HO::setContourIntegrationBasedOnGaussQuadratures();
  } else if (METHOD_TO_INTEGRATE_ALONG_CURVED_EDGES == ADAPTIVE_LINE_SEGMENTS_METHOD){
    // Set the affected switch
    Grid2D_Quad_Block_HO::setContourIntegrationBasedOnLinearSegments();
  } else if (METHOD_TO_INTEGRATE_ALONG_CURVED_EDGES == MIXED_QUADS_AND_LINES_METHOD){
    // Set the affected switch
    Grid2D_Quad_Block_HO::setMixedContourIntegration();
  }

  // Set the affected switch
  switch(NUMBER_OF_POINTS_FOR_GAUSS_QUADRATURE_INTEGRATION){
  case 3:
    Spline2DInterval_HO::setThreePointGaussQuadContourIntegration();
    break;
  case 5:
    Spline2DInterval_HO::setFivePointGaussQuadContourIntegration();
    break;
  }

  if (EXERCISE_SPECIAL_CARE_TO_ROUNDOFF_ERRORS_FOR_THIS_MESH == ON){
    // Set the affected switch
    Grid2D_Quad_Block_HO::setTreatMeshWithExtraCareForNumericalError();
  } else {
    // Set the affected switch
    Grid2D_Quad_Block_HO::setNoSpecialTreatmentForNumericalError();
  }
}

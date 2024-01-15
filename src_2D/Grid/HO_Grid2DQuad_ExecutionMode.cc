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
short HO_Grid2D_Execution_Mode::CHECK_FOR_INCORRECT_QUADRILATERALS = ON; // check for incorrect quads
short HO_Grid2D_Execution_Mode::REPORT_INCORRECT_QUADRILATERALS_BUT_CONTINUE_EXECUTION = OFF; // stop on detecting incorrect quads
short HO_Grid2D_Execution_Mode::USE_BROADCAST_MORE_THAN_RECOMPUTING = ON; // broadcast the majority of geometric properties
short HO_Grid2D_Execution_Mode::TOLERATE_INACCURATE_INTEGRATION_NEAR_CURVED_BOUNDARIES = OFF; // tolerate inaccurate integration
short HO_Grid2D_Execution_Mode::SMOOTH_QUAD_BLOCK_FLAG = ON; // smooth the grid
short HO_Grid2D_Execution_Mode::POLYGONAL_ADAPTIVE_QUADRATURE_INTEGRATION_FLAG = ON; // use polygonal adaptive integration
short HO_Grid2D_Execution_Mode::MONTE_CARLO_INTEGRATION_FLAG = OFF; // don't use Monte Carlo integration unless user specifies
short HO_Grid2D_Execution_Mode::POLYGONAL_ADAPTIVE_QUADRATURE_INTEGRATION_MINIMUM_LEVELS = 2;
short HO_Grid2D_Execution_Mode::ENFORCE_NONREFLECTED_SOUTH_BOUNDARY_GHOST_CELLS = OFF; // leave the BCs to take care of

// Block boundary flux calculation method
short HO_Grid2D_Execution_Mode::CUSTOMIZE_FLUX_CALCULATION_METHOD_AT_BOUNDARIES = OFF; // use the default settings
short HO_Grid2D_Execution_Mode::LOOPOVER_FLUX_CALCULATION_METHOD_AT_BOUNDARIES = OFF; // don't loop over the flux setting
short HO_Grid2D_Execution_Mode::WEST_RECONSTRUCTION_BASED_FLUX = OFF; // set to SolveRiemannProblem to compute the flux
short HO_Grid2D_Execution_Mode::SOUTH_RECONSTRUCTION_BASED_FLUX = OFF; // set to SolveRiemannProblem to compute the flux
short HO_Grid2D_Execution_Mode::NORTH_RECONSTRUCTION_BASED_FLUX = OFF; // set to SolveRiemannProblem to compute the flux
short HO_Grid2D_Execution_Mode::EAST_RECONSTRUCTION_BASED_FLUX = OFF; // set to SolveRiemannProblem to compute the flux


//! Set all flags to default values
// add all flag default values to this function
void HO_Grid2D_Execution_Mode::SetDefaults(void){

  USE_HIGH_ORDER_GEOMETRIC_BOUNDARY_REPRESENTATION = OFF; // High-order boundaries are not used
  METHOD_TO_INTEGRATE_ALONG_CURVED_EDGES = MIXED_QUADS_AND_LINES_METHOD; // Gauss quadrature integration mode
  NUMBER_OF_POINTS_FOR_GAUSS_QUADRATURE_INTEGRATION = 3; // use 3-point Gauss integration
  EXERCISE_SPECIAL_CARE_TO_ROUNDOFF_ERRORS_FOR_THIS_MESH = OFF; // use Global Coordinate System
  CHECK_FOR_INCORRECT_QUADRILATERALS = ON; // check for incorrect quads
  REPORT_INCORRECT_QUADRILATERALS_BUT_CONTINUE_EXECUTION = OFF; // stop on detecting incorrect quads
  USE_BROADCAST_MORE_THAN_RECOMPUTING = ON; // broadcast the majority of geometric properties
  CUSTOMIZE_FLUX_CALCULATION_METHOD_AT_BOUNDARIES = OFF; // use the default settings
  LOOPOVER_FLUX_CALCULATION_METHOD_AT_BOUNDARIES = OFF; // don't loop over the flux setting
  WEST_RECONSTRUCTION_BASED_FLUX = OFF; // set to SolveRiemannProblem to compute the flux
  SOUTH_RECONSTRUCTION_BASED_FLUX = OFF; // set to SolveRiemannProblem to compute the flux
  NORTH_RECONSTRUCTION_BASED_FLUX = OFF; // set to SolveRiemannProblem to compute the flux
  EAST_RECONSTRUCTION_BASED_FLUX = OFF; // set to SolveRiemannProblem to compute the flux
  TOLERATE_INACCURATE_INTEGRATION_NEAR_CURVED_BOUNDARIES = OFF; // tolerate inaccurate integration
  SMOOTH_QUAD_BLOCK_FLAG = ON; // smooth the grid
  POLYGONAL_ADAPTIVE_QUADRATURE_INTEGRATION_FLAG = ON; // use polygonal adaptive integration
  MONTE_CARLO_INTEGRATION_FLAG = OFF; // don't use Monte Carlo integration unless user specifies
  POLYGONAL_ADAPTIVE_QUADRATURE_INTEGRATION_MINIMUM_LEVELS = 2;
  ENFORCE_NONREFLECTED_SOUTH_BOUNDARY_GHOST_CELLS = OFF; // leave the boundary condition to take care of

  // Reset solid body counter in Spline2D_HO class
  Spline2D_HO::ResetCounter();
}

//! Print the current execution mode
//  at the output stream
// \param [in] out_file the output stream
void HO_Grid2D_Execution_Mode::Print_Info(std::ostream & out_file){

  // output flux calculation method through block boundaries
  if (CUSTOMIZE_FLUX_CALCULATION_METHOD_AT_BOUNDARIES == OFF){
    out_file << "\n  -> Flux Calculation Method: " << "Solve a 'Riemann' problem everywhere";

  } else {
    out_file << "\n  -> Flux Calculation Method: " << "Customized for each boundary";

    if (WEST_RECONSTRUCTION_BASED_FLUX == ON){
      out_file << "\n     -> West Flux Method: " << "Based on constrained reconstruction";
    } else {
      out_file << "\n     -> West Flux Method: " << "Solve 'Riemann' problem";
    }

    if (SOUTH_RECONSTRUCTION_BASED_FLUX == ON){
      out_file << "\n     -> South Flux Method: " << "Based on constrained reconstruction";
    } else {	     	  
      out_file << "\n     -> South Flux Method: " << "Solve a 'Riemann' problem";
    }

    if (EAST_RECONSTRUCTION_BASED_FLUX == ON){
      out_file << "\n     -> East Flux Method: " << "Based on constrained reconstruction";
    } else {	     	  
      out_file << "\n     -> East Flux Method: " << "Solve a 'Riemann' problem";
    }

    if (NORTH_RECONSTRUCTION_BASED_FLUX == ON){
      out_file << "\n     -> North Flux Method: " << "Based on constrained reconstruction";
    } else {	     	  
      out_file << "\n     -> North Flux Method: " << "Solve a 'Riemann' problem";
    }
  }

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

    if (POLYGONAL_ADAPTIVE_QUADRATURE_INTEGRATION_FLAG){
      out_file << "\n     -> Integration Along Curved Edges: " << "Polygonal adaptive quadrature"
	       << "\n     -> Minimum refinement levels: " << POLYGONAL_ADAPTIVE_QUADRATURE_INTEGRATION_MINIMUM_LEVELS;
    } else if (MONTE_CARLO_INTEGRATION_FLAG){
      out_file << "\n     -> Integration Along Curved Edges: " 
	       << "Monte Carlo, #Samples(" 
	       << NumericalLibrary_Execution_Mode::Number_Monte_Carlo_Samples 
	       << ")";
    } else if (TOLERATE_INACCURATE_INTEGRATION_NEAR_CURVED_BOUNDARIES ){
      out_file << "\n     -> Integration Along Curved Edges: " << "Force with straight edges";
    } // endif   

  } // endif

  if (EXERCISE_SPECIAL_CARE_TO_ROUNDOFF_ERRORS_FOR_THIS_MESH == ON){
    out_file << "\n     -> Geometric properties: " << "Special attention to Round Offs";
  }

#ifdef _MPI_VERSION
  if (USE_BROADCAST_MORE_THAN_RECOMPUTING == ON){
    out_file << "\n     -> MPI transfer of geometric properties: Broadcast";
  } else {
    out_file << "\n     -> MPI transfer of geometric properties: Recalculation";
  }
#endif

  if (CHECK_FOR_INCORRECT_QUADRILATERALS == OFF){
    out_file << "\n     -> Warning: " << "Mesh validity checking was completely deactivated!!!";
  } else if (REPORT_INCORRECT_QUADRILATERALS_BUT_CONTINUE_EXECUTION == ON){
    out_file << "\n     -> Warning: " << "Program won't stop the execution if invalid meshes are detected!!!";
  }

  if (SMOOTH_QUAD_BLOCK_FLAG == ON){
    out_file << "\n     -> Smooth Quad Block: Yes";
  } else {
    out_file << "\n     -> Smooth Quad Block: No";
  }

  if (ENFORCE_NONREFLECTED_SOUTH_BOUNDARY_GHOST_CELLS == ON){
    out_file << "\n     -> Ghost cells South boundary: Enforce non-reflected geometry";
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
  MPI::COMM_WORLD.Bcast(&CHECK_FOR_INCORRECT_QUADRILATERALS,
 			1, 
 			MPI::SHORT, 0);
  MPI::COMM_WORLD.Bcast(&REPORT_INCORRECT_QUADRILATERALS_BUT_CONTINUE_EXECUTION,
 			1, 
 			MPI::SHORT, 0);
  MPI::COMM_WORLD.Bcast(&USE_BROADCAST_MORE_THAN_RECOMPUTING,
 			1, 
 			MPI::SHORT, 0);
  MPI::COMM_WORLD.Bcast(&TOLERATE_INACCURATE_INTEGRATION_NEAR_CURVED_BOUNDARIES,
 			1, 
 			MPI::SHORT, 0);

  // broadcast flux calculation methods
  MPI::COMM_WORLD.Bcast(&CUSTOMIZE_FLUX_CALCULATION_METHOD_AT_BOUNDARIES,
 			1, 
 			MPI::SHORT, 0);
  MPI::COMM_WORLD.Bcast(&WEST_RECONSTRUCTION_BASED_FLUX,
 			1, 
 			MPI::SHORT, 0);
  MPI::COMM_WORLD.Bcast(&EAST_RECONSTRUCTION_BASED_FLUX,
 			1, 
 			MPI::SHORT, 0);
  MPI::COMM_WORLD.Bcast(&SOUTH_RECONSTRUCTION_BASED_FLUX,
 			1, 
 			MPI::SHORT, 0);
  MPI::COMM_WORLD.Bcast(&NORTH_RECONSTRUCTION_BASED_FLUX,
 			1, 
 			MPI::SHORT, 0);


  MPI::COMM_WORLD.Bcast(&SMOOTH_QUAD_BLOCK_FLAG,
 			1, 
 			MPI::SHORT, 0);

  MPI::COMM_WORLD.Bcast(&POLYGONAL_ADAPTIVE_QUADRATURE_INTEGRATION_FLAG,
 			1, 
 			MPI::SHORT, 0);

  MPI::COMM_WORLD.Bcast(&MONTE_CARLO_INTEGRATION_FLAG,
 			1, 
 			MPI::SHORT, 0);

  MPI::COMM_WORLD.Bcast(&POLYGONAL_ADAPTIVE_QUADRATURE_INTEGRATION_MINIMUM_LEVELS,
 			1, 
 			MPI::SHORT, 0);

  MPI::COMM_WORLD.Bcast(&ENFORCE_NONREFLECTED_SOUTH_BOUNDARY_GHOST_CELLS,
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

  if (TOLERATE_INACCURATE_INTEGRATION_NEAR_CURVED_BOUNDARIES == ON){
    // Set the affected switch
    Grid2D_Quad_Block_HO::setTolerateInaccurateIntegrationNearCurvedBoundaries();
  } else {
    // Set the affected switch
    Grid2D_Quad_Block_HO::setNoGeometricInaccuraciesForIntegrationNearCurvedBoundaries();
  }

  if (SMOOTH_QUAD_BLOCK_FLAG == ON){
    // Set the affected switch
    Grid2D_Quad_Block_HO::setGridSmoothing();
  } else {
    // Set the affected switch
    Grid2D_Quad_Block_HO::setNoGridSmoothing();
  }

  if (POLYGONAL_ADAPTIVE_QUADRATURE_INTEGRATION_FLAG){
    // Set the affected switches
    Grid2D_Quad_Block_HO::setPolygonalAdaptiveQuadratureIntegrationON();
    Grid2D_Quad_Block_HO::setMonteCarloIntegrationOFF();
    Grid2D_Quad_Block_HO::Polygonal_Adaptive_Quadrature_Integration_Minimum_Levels = 
      POLYGONAL_ADAPTIVE_QUADRATURE_INTEGRATION_MINIMUM_LEVELS;

  } else if (MONTE_CARLO_INTEGRATION_FLAG){
    // Set the affected switches
    Grid2D_Quad_Block_HO::setPolygonalAdaptiveQuadratureIntegrationOFF();
    Grid2D_Quad_Block_HO::setMonteCarloIntegrationON();

  } else {
    // Set the affected switches
    Grid2D_Quad_Block_HO::setPolygonalAdaptiveQuadratureIntegrationOFF();
    Grid2D_Quad_Block_HO::setMonteCarloIntegrationOFF();
  }

  if (ENFORCE_NONREFLECTED_SOUTH_BOUNDARY_GHOST_CELLS == ON){
    // Set the affected switch
    Grid2D_Quad_Block_HO::setNonReflectedGhostCellsNearSouthSolidBoundary();
  } else {
    // Set the affected switch
    Grid2D_Quad_Block_HO::setReflectedGhostCellsNearSouthSolidBoundary();
  }
}

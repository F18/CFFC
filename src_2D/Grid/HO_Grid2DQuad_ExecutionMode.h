/*!\file HO_Grid2DQuad_ExecutionMode.h
  \brief Definition of flags that control the execution of high-order 2D grid. */

#ifndef _HO_GRID2DQUAD_EXECUTIONMODE_INCLUDED
#define _HO_GRID2DQUAD_EXECUTIONMODE_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../Utilities/Utilities.h"
#include "../CFD/CFD.h"
#include "HO_Grid2DQuad.h"

/*!
 * \class HO_Grid2D_Execution_Mode
 *
 * @brief Definition of flags that control the execution of high-order 2D grid
 *
 */
class HO_Grid2D_Execution_Mode{
  
public:

  //! Structure that defines the available methods for integrating along curved boundaries.
  enum CurvilinearIntegrationMethod { GAUSS_QUADRATURES_METHOD, ADAPTIVE_LINE_SEGMENTS_METHOD, MIXED_QUADS_AND_LINES_METHOD };

  // set all flags to default values
  static void SetDefaults(void);

  /* This flag controls whether cell curved edges are 
     represented as straight lines or high-order polynomials.
     Turn ON if you want to use the high-order representation.
     Turn OFF if you want to approximate the edges with straight lines. (default)
     ATTENTION: Some settings/memory allocations are done only if this flag is ON.
     ----------------------------------------------------------------------------------------  */
  static short USE_HIGH_ORDER_GEOMETRIC_BOUNDARY_REPRESENTATION;

  
  /* This flag sets the method used to integrate along curved boundaries.
     There are three methods currently available:
         -> integration with 3-point or 5-point Gauss quadratures. (less expensive in general)
	 -> adaptive integration based on approximating the curve edges with line segments. (tipical more accurate)
         -> a combination of the above methods. This is a tradeoff between speed and accuracy.
	    The accuracy of the geometric moments is strictly dependent on the accuracy with which
            the cell area and the centroid are computed. Therefore, this method computes the area 
            and the centroid with the adaptive line segments and the higher-order geometric moments
            with the Gauss quadratures.
     Set to GAUSS_QUADRATURES_METHOD if you want to use Gauss quadratures.
     Set to ADAPTIVE_LINE_SEGMENTS_METHOD if you want to use adaptive integration based on line segments.
     Set to MIXED_QUADS_AND_LINES_METHOD if you want to use the combination of the above methods. (default)
     ---------------------------------------------------------------------------------------- */
  static short METHOD_TO_INTEGRATE_ALONG_CURVED_EDGES;


  /* This flag sets the number of Gauss points used to integrate along curved boundaries,
     if integration with Gauss quadratures is selected.
     Use 3 points along each smooth spline segment. (default)
     The other option currently available is 5 points.
     ----------------------------------------------------------------------------------------  */
  static short NUMBER_OF_POINTS_FOR_GAUSS_QUADRATURE_INTEGRATION;  


  /* This flag sets the designated mesh switch to inform the algorithm
     for computing cell geometric properties about the possibility of large
     round off errors due to the presence of large mesh dimensions.
     The numerical algorithm will try to avoid this problems by computing
     the geometric properties in a local coordinate system correspondent 
     to each computational cell.
     Turn ON if you want to use the local coordinate system (LCS). 
     Turn OFF if you want to have all calculations done in a single coordinate system. (default)
     ----------------------------------------------------------------------------------------  */
  static short EXERCISE_SPECIAL_CARE_TO_ROUNDOFF_ERRORS_FOR_THIS_MESH;  

  /* This flag is used to control the behaviour of the routine that 
     checks the validity of the mesh (i.e. presence of crossed or degenerated quadrilateral cells).
     Turn ON if you want to check for incorrect quadrilaterals. (default)
     Turn OFF if you don't want this subroutine to check for incorrect quadrilaterals.
     Note: It this flag is OFF, the subroutines is not even going to report the incorrect quads!
     ----------------------------------------------------------------------------------------  */
  static short CHECK_FOR_INCORRECT_QUADRILATERALS;
 
  /* This flag is used to control the behaviour of the routine that 
     checks the validity of the mesh.
     The default behavior of this subroutine when an invalid quadrilateral cell is
     encountered is to report the cell and to stop the execution. 
     There are moments when such a behavior is not desired (e.g. debugging), and for this
     situations this flag can be used to modify the default behavior.
     Turn ON if you want the subroutine to continue the execution and to report (i.e. output on screen)
     all incorrect quadrilaterals.
     Turn OFF if you want the mesh checking routine to report and exit at the first encountered
     incorrect quadrilateral. (default)
     ----------------------------------------------------------------------------------------  */
  static short REPORT_INCORRECT_QUADRILATERALS_BUT_CONTINUE_EXECUTION;

  template<class Input_Parameters_Type>
  static void Parse_Next_Input_Control_Parameter(Input_Parameters_Type & IP, int & i_command);

  static void Print_Info(std::ostream & out_file);

  static void Broadcast(void);

  static void Set_Affected_Switches(void);

protected:
  HO_Grid2D_Execution_Mode(void);   //!< Private default constructor
  HO_Grid2D_Execution_Mode(const HO_Grid2D_Execution_Mode&); //!< Private copy constructor
  HO_Grid2D_Execution_Mode& operator=(const HO_Grid2D_Execution_Mode&); //!< Private assignment operator

};

//! Parse the input control parameters for 
//  settings related to HO_Grid2D_Execution_Mode class
template<class Input_Parameters_Type> inline
void HO_Grid2D_Execution_Mode::Parse_Next_Input_Control_Parameter(Input_Parameters_Type & IP, int & i_command){

  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }

  char buffer[256];

  // Try to match the next control parameter
  if (strcmp(IP.Next_Control_Parameter, "High_Order_Boundary") == 0) {
    IP.Get_Next_Input_Control_Parameter();
    if ( strcmp(IP.Next_Control_Parameter, "ON") == 0 || strcmp(IP.Next_Control_Parameter,"On") == 0 ){
      USE_HIGH_ORDER_GEOMETRIC_BOUNDARY_REPRESENTATION = ON;
      // Set the affected switch
      Grid2D_Quad_Block_HO::setHighOrderBoundaryRepresentation();

    } else if ( strcmp(IP.Next_Control_Parameter, "OFF") == 0 || strcmp(IP.Next_Control_Parameter,"Off") == 0) {
      USE_HIGH_ORDER_GEOMETRIC_BOUNDARY_REPRESENTATION = OFF;
      // Set the affected switch
      Grid2D_Quad_Block_HO::setLowOrderBoundaryRepresentation();

    } else {
      i_command = INVALID_INPUT_VALUE;
      return;
    }
    i_command = 0;

  } else if (strcmp(IP.Next_Control_Parameter, "Curved_Boundary_Integration") == 0) {
    IP.Get_Next_Input_Control_Parameter();
    if ( strcmp(IP.Next_Control_Parameter, "Gauss_Quad") == 0 ){
      METHOD_TO_INTEGRATE_ALONG_CURVED_EDGES = GAUSS_QUADRATURES_METHOD;
      // Set the affected switches
      Grid2D_Quad_Block_HO::setContourIntegrationBasedOnGaussQuadratures();

    } else if ( strcmp(IP.Next_Control_Parameter, "Line_Segments") == 0 ) {
      METHOD_TO_INTEGRATE_ALONG_CURVED_EDGES = ADAPTIVE_LINE_SEGMENTS_METHOD;
      // Set the affected switches
      Grid2D_Quad_Block_HO::setContourIntegrationBasedOnLinearSegments();

    } else if ( strcmp(IP.Next_Control_Parameter, "Mixed_Gauss_Quad") == 0 ) {
      METHOD_TO_INTEGRATE_ALONG_CURVED_EDGES = MIXED_QUADS_AND_LINES_METHOD;
      // Set the affected switches
      Grid2D_Quad_Block_HO::setMixedContourIntegration();
 
    } else {
      i_command = INVALID_INPUT_VALUE;
      return;
    }
    i_command = 0;

  } else if (strcmp(IP.Next_Control_Parameter, "Gauss_Points_For_Boundary_Integration") == 0) {
    i_command = 0;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> NUMBER_OF_POINTS_FOR_GAUSS_QUADRATURE_INTEGRATION;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if ((NUMBER_OF_POINTS_FOR_GAUSS_QUADRATURE_INTEGRATION != 3) && 
	(NUMBER_OF_POINTS_FOR_GAUSS_QUADRATURE_INTEGRATION != 5) ){
      i_command = INVALID_INPUT_VALUE;
    } else {
      // Set the affected switch
      switch(NUMBER_OF_POINTS_FOR_GAUSS_QUADRATURE_INTEGRATION){
      case 3:
	Spline2DInterval_HO::setThreePointGaussQuadContourIntegration();
	break;
      case 5:
	Spline2DInterval_HO::setFivePointGaussQuadContourIntegration();
	break;
      }
    }

  } else if (strcmp(IP.Next_Control_Parameter, "Large_Domain_Mesh") == 0) {
    IP.Get_Next_Input_Control_Parameter();
    if ( strcmp(IP.Next_Control_Parameter, "ON") == 0 || strcmp(IP.Next_Control_Parameter, "On") == 0 ){
      EXERCISE_SPECIAL_CARE_TO_ROUNDOFF_ERRORS_FOR_THIS_MESH = ON;
      // Set the affected switch
      Grid2D_Quad_Block_HO::setTreatMeshWithExtraCareForNumericalError();

    } else if ( strcmp(IP.Next_Control_Parameter, "OFF") == 0 || strcmp(IP.Next_Control_Parameter, "Off") == 0 ) {
      EXERCISE_SPECIAL_CARE_TO_ROUNDOFF_ERRORS_FOR_THIS_MESH = OFF;
      // Set the affected switch
      Grid2D_Quad_Block_HO::setNoSpecialTreatmentForNumericalError();
 
    } else {
      i_command = INVALID_INPUT_VALUE;
      return;
    }
    i_command = 0;

  } else if (strcmp(IP.Next_Control_Parameter, "Check_Mesh_Validity") == 0) {
    IP.Get_Next_Input_Control_Parameter();
    if ( strcmp(IP.Next_Control_Parameter, "ON") == 0 || strcmp(IP.Next_Control_Parameter, "On") == 0 ){
      CHECK_FOR_INCORRECT_QUADRILATERALS = ON;
    } else if ( strcmp(IP.Next_Control_Parameter, "OFF") == 0 || strcmp(IP.Next_Control_Parameter, "Off") == 0 ) {
      CHECK_FOR_INCORRECT_QUADRILATERALS = OFF;
    } else {
      i_command = INVALID_INPUT_VALUE;
      return;
    }
    i_command = 0;

  } else if (strcmp(IP.Next_Control_Parameter, "If_Invalid_Mesh_Detected") == 0) {
    IP.Get_Next_Input_Control_Parameter();
    if ( strcmp(IP.Next_Control_Parameter, "Stop") == 0 ){
      REPORT_INCORRECT_QUADRILATERALS_BUT_CONTINUE_EXECUTION = OFF;
    } else if ( strcmp(IP.Next_Control_Parameter, "Continue") == 0 ) {
      REPORT_INCORRECT_QUADRILATERALS_BUT_CONTINUE_EXECUTION = ON;
    } else {
      i_command = INVALID_INPUT_VALUE;
      return;
    }
    i_command = 0;

  } else {
    i_command = INVALID_INPUT_CODE;
  } // endif

}

#endif

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
 * Class: HO_Grid2D_Execution_Mode
 *
 * @brief Definition of flags that control the execution of high-order 2D grid
 *
 */
class HO_Grid2D_Execution_Mode{
  
public:

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
     There are two methods currently available:
         -> integration with Gauss quadratures. (less expensive)
	 -> adaptive integration based on approximating the curve edges with line segments. (more accurate)
     Turn ON if you want to use Gauss quadratures. (default)
     Turn OFF if you want to use adaptive integration based on line segments.
     ---------------------------------------------------------------------------------------- */
  static short USE_GAUSS_QUADRATURES_TO_INTEGRATE_ALONG_CURVED_EDGES;


  /* This flag sets the number of Gauss points used to integrate along curved boundaries,
     if integration with Gauss quadratures is selected.
     Use 3 points along each smooth spline segment. (default)
     The other option currently available is 5 points.
     ----------------------------------------------------------------------------------------  */
  static short NUMBER_OF_POINTS_FOR_GAUSS_QUADRATURE_INTEGRATION;  

  
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
      USE_GAUSS_QUADRATURES_TO_INTEGRATE_ALONG_CURVED_EDGES = ON;
      // Set the affected switch
      Grid2D_Quad_Block_HO::setContourIntegrationBasedOnGaussQuadratures();

    } else if ( strcmp(IP.Next_Control_Parameter, "Line_Segments") == 0 ) {
      USE_GAUSS_QUADRATURES_TO_INTEGRATE_ALONG_CURVED_EDGES = OFF;
      // Set the affected switch
      Grid2D_Quad_Block_HO::setContourIntegrationBasedOnLinearSegments();
 
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

  } else {
    i_command = INVALID_INPUT_CODE;
  } // endif

}

#endif

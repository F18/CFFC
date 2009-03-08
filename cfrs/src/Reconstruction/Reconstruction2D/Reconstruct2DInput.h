/* Reconstruct2DInput.h:  Header file defining 
                          2D Reconstruction Input Parameter Class. */

#ifndef _RECONSTRUCT2D_INPUT_INCLUDED
#define _RECONSTRUCT2D_INPUT_INCLUDED

/* Include header files. */

#include <fstream>
#include "CFD/CFD.h"
#include "Grid/Grid2D/Cell2D.h"
#include "TestFunctions/TestFunctions_2D.h"
#include "include/TypeDefinition.h"
#include "../CENO_CFRS_ExecutionMode.h" // Include high-order CENO execution mode header file

using namespace std;

/* Define the structures and classes. */

#define	INPUT_PARAMETER_LENGTH_RECONSTRUCT2D    128

typedef double (* TestFunction2D) (double,double);
typedef double (* IntegralTestFunction2D) (double,double,double,double);


/**********************************************************
 * Class:  Reconstruct2D_Input_Parameters               *
 **********************************************************/
class Reconstruct2D_Input_Parameters{
  public:
  static const SpaceType IP_SpaceDimension = TwoD; 

  // Input file name:
  char Input_File_Name[INPUT_PARAMETER_LENGTH_RECONSTRUCT2D];

  // Input file stream:
  ifstream Input_File;

  // Input file line number:
  int Line_Number;
  int Message_Number;

  // Parameters characteristic for the method involved for solving the 
  // linear system of equations
  char Method_Used[INPUT_PARAMETER_LENGTH_RECONSTRUCT2D];
  int  Method;
  char Geometric_Weighting[INPUT_PARAMETER_LENGTH_RECONSTRUCT2D];
  int  geom_weighting;
  char Data_Dependent_Weighting[INPUT_PARAMETER_LENGTH_RECONSTRUCT2D];
  int  data_depend_weighting;
  char Pseudo_Inverse[INPUT_PARAMETER_LENGTH_RECONSTRUCT2D];

  // Function file name to be reconstructed and related input parameters:
  char Function_Type[INPUT_PARAMETER_LENGTH_RECONSTRUCT2D];
  int Reconstruction_Order;
  int i_Function;
  TestFunction2D TestF;
  IntegralTestFunction2D IntTestF;
  char Integration_Type[INPUT_PARAMETER_LENGTH_RECONSTRUCT2D];

  // Initial condition type indicator and related input parameters:
  char ICs_Type[INPUT_PARAMETER_LENGTH_RECONSTRUCT2D];
  int i_ICs;

  // Grid type indicator and related input parameters:
  char Grid_Type[INPUT_PARAMETER_LENGTH_RECONSTRUCT2D];
  char NACA_Aerofoil_Type[INPUT_PARAMETER_LENGTH_RECONSTRUCT2D];
  int i_Grid;
  char CellNumber_or_DeltaCell[2];
  int  Number_of_Cells_Idir, Number_of_Cells_Jdir,
       Number_of_Blocks_Idir, Number_of_Blocks_Jdir;
  int Number_of_Ghost_Cells;
  double Delta_X_Cell;
  double Delta_Y_Cell;

  double Box_Width, Box_Height, Plate_Length, 
         Pipe_Length, Pipe_Radius, 
         Blunt_Body_Radius, Blunt_Body_Mach_Number,
         Grain_Length, Grain_Radius, Grain_To_Throat_Length,
         Nozzle_Length, Nozzle_Radius_Exit, Nozzle_Radius_Throat,
         Cylinder_Radius, Ellipse_Length_X_Axis, 
         Ellipse_Length_Y_Axis, Chord_Length, Orifice_Radius,
         Wedge_Angle, Wedge_Length;

  double X_Scale, X_Rotate;
  Vector2D X_Shift;
  char Stretching_Function_I[INPUT_PARAMETER_LENGTH_RECONSTRUCT2D]; /* Grid stretching function */
  char Stretching_Function_J[INPUT_PARAMETER_LENGTH_RECONSTRUCT2D]; /* Grid stretching function */
  int Stretch_I, Stretch_J;	
  double Beta_I, Beta_J, Tau_I, Tau_J; /* Grid stretching parameters */
  int NumOfIter_UnsmoothMesh;
  double CharacteristicLength;
  double CutoffKnobHighOrderReconstruction;
  double CENO_Cutoff;

  char **ICEMCFD_FileNames;

  //Subgrid parameters 
  char SubGridPoints_or_DeltaSubGrid[2];
  int  Number_of_SubGrid_Points_Idir;
  int  Number_of_SubGrid_Points_Jdir;
  double Delta_X_of_SubGrid;
  double Delta_Y_of_SubGrid;

  double X_min, X_Length;
  double Y_min, Y_Length;

  // Limiter type indicator and related input parameters:
  char Limiter_Type[INPUT_PARAMETER_LENGTH_RECONSTRUCT2D];
  int i_Limiter;

  // Output file name:
  char Output_File_Name[INPUT_PARAMETER_LENGTH_RECONSTRUCT2D];

  // Multi-block mesh definition input file names:
  char Grid_File_Name[INPUT_PARAMETER_LENGTH_RECONSTRUCT2D];
  char Grid_Definition_File_Name[INPUT_PARAMETER_LENGTH_RECONSTRUCT2D];

  // Gnuplot file name:
  char Gnuplot_File_Name[INPUT_PARAMETER_LENGTH_RECONSTRUCT2D];

  // Next_Control_Parameter:
  char Next_Control_Parameter[INPUT_PARAMETER_LENGTH_RECONSTRUCT2D];

  // Output format type indicator:
  char Output_Format_Type[INPUT_PARAMETER_LENGTH_RECONSTRUCT2D];
  int i_Output_Format;

  // Access functions
  const int & NumSubGridI(void) const { return Number_of_SubGrid_Points_Idir; }
  const int & NumSubGridJ(void) const { return Number_of_SubGrid_Points_Jdir; }
  const int & RecOrder(void) const { return Reconstruction_Order; }
  const int & iCell(void) const { return Number_of_Cells_Idir;}
  const int & jCell(void) const { return Number_of_Cells_Jdir;}
  int kCell(void) const { return 0;}
  const int & Nghost(void) const { return Number_of_Ghost_Cells;}
  const TestFunction2D & TestFunction(void) const {return TestF;}
  double & CutoffKnob(void) { return CutoffKnobHighOrderReconstruction; }
  const double & CutoffKnob(void) const { return CutoffKnobHighOrderReconstruction; }
  int & ReconstructionMethod(void) {return Method;}
  const int & ReconstructionMethod(void) const {return Method;}
  const int & Limiter(void) const {return i_Limiter;}
  double & FitTolerance(void) {return CENO_Cutoff;}
  const double & FitTolerance(void) const {return CENO_Cutoff;}

  /* Input-output operators. */

  friend ostream &operator << (ostream &out_file,
		               const Reconstruct2D_Input_Parameters &IP);
  friend istream &operator >> (istream &in_file,
			       Reconstruct2D_Input_Parameters &IP);

};

/*************************************************************
 * Reconstruct2D_Input_Parameters -- External subroutines.         *
 *************************************************************/

extern void Open_Input_File(Reconstruct2D_Input_Parameters &IP);

extern void Close_Input_File(Reconstruct2D_Input_Parameters &IP);

extern void Set_Default_Input_Parameters(Reconstruct2D_Input_Parameters &IP);

extern void Get_Next_Input_Control_Parameter(Reconstruct2D_Input_Parameters &IP);

extern int Parse_Next_Input_Control_Parameter(Reconstruct2D_Input_Parameters &IP);

extern int Process_Input_Control_Parameter_File(Reconstruct2D_Input_Parameters
						&Input_Parameters,
                                                char *Input_File_Name_ptr,
                                                int &Command_Flag);

#endif /* _RECONSTRUCT2D_INPUT_INCLUDED  */

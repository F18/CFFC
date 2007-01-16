/* Reconstruct1DInput.h:  Header file defining 
                          1D Reconstruction Input Parameter Class. */

#ifndef _RECONSTRUCT1D_INPUT_INCLUDED
#define _RECONSTRUCT1D_INPUT_INCLUDED

/* Include header files. */

#include <fstream>
using namespace std;

#ifndef _CFD_INCLUDED
#include "CFD/CFD.h"
#endif // _CFD_INCLUDED

#ifndef _CELL1D_INCLUDED
#include "Grid/Grid1D/Cell1D.h"
#endif // _CELL1D_INCLUDED

#ifndef _TESTFUNCTIONS_INCLUDED
#include "TestFunctions/TestFunctions.h"
#endif //_TESTFUNCTIONS_INCLUDED

#ifndef _TYPEDEFINITION_INCLUDED
#include "include/TypeDefinition.h"
#endif // _TYPEDEFINITION_INCLUDED


/* Define the structures and classes. */

#define	INPUT_PARAMETER_LENGTH_RECONSTRUCT1D    128

typedef double (* TestFunction1D) (double);
typedef double (* IntegralTestFunction1D) (double,double);

/**********************************************************
 * Class:  Reconstruct1D_Input_Parameters               *
 **********************************************************/
class Reconstruct1D_Input_Parameters{
  public:
  static const SpaceType IP_SpaceDimension = OneD; 

  // Input file name:
  char Input_File_Name[INPUT_PARAMETER_LENGTH_RECONSTRUCT1D];

  // Input file stream:
  ifstream Input_File;

  // Input file line number:
  int Line_Number;
  int Message_Number;

  // linear system of equations
  // Characteristic parameters for the reconstruction method
  char Method_Used[INPUT_PARAMETER_LENGTH_RECONSTRUCT1D];
  int  Method;
  char Geometric_Weighting[INPUT_PARAMETER_LENGTH_RECONSTRUCT1D];
  bool geom_weighting;
  char Data_Dependent_Weighting[INPUT_PARAMETER_LENGTH_RECONSTRUCT1D];
  bool data_depend_weighting;
  char Use_Residual_DI[INPUT_PARAMETER_LENGTH_RECONSTRUCT1D];
  bool use_residual_DI;
  double CutoffKnobHighOrderReconstruction;
  double CENO_Cutoff;

  // Function file name to be reconstructed and related input parameters:
  char Function_Type[INPUT_PARAMETER_LENGTH_RECONSTRUCT1D];
  int Reconstruction_Order;
  int i_Function;
  TestFunction1D TestF;
  IntegralTestFunction1D IntTestF;
  char Integration_Type[INPUT_PARAMETER_LENGTH_RECONSTRUCT1D];

  // Grid type indicator and related input parameters:
  char Grid_Type[INPUT_PARAMETER_LENGTH_RECONSTRUCT1D];
  int  i_Grid;
  char CellNumber_or_DeltaCell[2];
  int  Number_of_Cells_Idir, Number_of_Ghost_Cells;
  double Delta_Cell;
  //Subgrid parameters 
  char SubGridPoints_or_DeltaSubGrid[2];
  int Number_of_SubGrid_Points;
  double Delta_of_SubGrid;
  
  double X_min;
  double X_Scale, X_Shift;
  double X_max;
  double CharacteristicLength;

  // Limiter type indicator and related input parameters:
  char Limiter_Type[INPUT_PARAMETER_LENGTH_RECONSTRUCT1D];
  int i_Limiter;
 
  // Output file name:
  char Output_File_Name[INPUT_PARAMETER_LENGTH_RECONSTRUCT1D];

  // Multi-block mesh definition input file names:
  char Grid_File_Name[INPUT_PARAMETER_LENGTH_RECONSTRUCT1D];
  char Grid_Definition_File_Name[INPUT_PARAMETER_LENGTH_RECONSTRUCT1D];

  // Gnuplot file name:
  char Gnuplot_File_Name[INPUT_PARAMETER_LENGTH_RECONSTRUCT1D];

  // Next_Control_Parameter:
  char Next_Control_Parameter[INPUT_PARAMETER_LENGTH_RECONSTRUCT1D];

  // Output format type indicator:
  char Output_Format_Type[INPUT_PARAMETER_LENGTH_RECONSTRUCT1D];
  int i_Output_Format;

  // Access functions
  const int & NumSubGridI(void) const { return Number_of_SubGrid_Points; }
  const int & RecOrder(void) const { return Reconstruction_Order; }
  const int & iCell(void) const { return Number_of_Cells_Idir;}
  int jCell(void) const { return 0;}
  int kCell(void) const { return 0;}
  const int & Nghost(void) const { return Number_of_Ghost_Cells;}
  const TestFunction1D & TestFunction(void) const {return TestF;}
  double & CutoffKnob(void) { return CutoffKnobHighOrderReconstruction; }
  const double & CutoffKnob(void) const { return CutoffKnobHighOrderReconstruction; }
  int & ReconstructionMethod(void) {return Method;}
  const int & ReconstructionMethod(void) const {return Method;}
  const double & MinX(void) const {return X_min;}
  double LengthX(void) const {return X_max - X_min;}
  const double & ScaleX(void) const {return X_Scale;}
  const double & ShiftX(void) const {return X_Shift;}
  const double & MaxX(void) const {return X_max;}
  const int & Limiter(void) const {return i_Limiter;}
  double & FitTolerance(void) {return CENO_Cutoff;}
  const double & FitTolerance(void) const {return CENO_Cutoff;}

  /* Input-output operators. */

  friend ostream &operator << (ostream &out_file,
		               const Reconstruct1D_Input_Parameters &IP);
  friend istream &operator >> (istream &in_file,
			       Reconstruct1D_Input_Parameters &IP);

};

/*************************************************************
 * Reconstruct1D_Input_Parameters -- External subroutines.         *
 *************************************************************/

extern void Open_Input_File(Reconstruct1D_Input_Parameters &IP);

extern void Close_Input_File(Reconstruct1D_Input_Parameters &IP);

extern void Set_Default_Input_Parameters(Reconstruct1D_Input_Parameters &IP);

extern void Get_Next_Input_Control_Parameter(Reconstruct1D_Input_Parameters &IP);

extern int Parse_Next_Input_Control_Parameter(Reconstruct1D_Input_Parameters &IP);

extern int Process_Input_Control_Parameter_File(Reconstruct1D_Input_Parameters &Input_Parameters,
                                                char *Input_File_Name_ptr,
                                                int &Command_Flag);

#endif /* _RECONSTRUCT1D_INPUT_INCLUDED  */

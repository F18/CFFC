/* Reconstruct3DInput.h:  Header file defining 
                          3D Reconstruction Input Parameter Class. */

#ifndef _RECONSTRUCT3D_INPUT_INCLUDED
#define _RECONSTRUCT3D_INPUT_INCLUDED

/* Include header files. */

#include <fstream>
#include "../../CFD/CFD.h"
#include "../../Grid/Grid3D/Cell3D.h"
#include "../../TestFunctions/TestFunctions_3D.h"
#include "../../include/TypeDefinition.h"
#include "../CENO_CFRS_ExecutionMode.h" // Include high-order CENO execution mode header file

using namespace std;

/* Define the structures and classes. */

#define	INPUT_PARAMETER_LENGTH_RECONSTRUCT3D    128

typedef double (* TestFunction3D) (double,double,double);
typedef double (* IntegralTestFunction3D) (double,double,double,double,double,double);


/**********************************************************
 * Class:  Reconstruct3D_Input_Parameters                 *
 **********************************************************/
class Reconstruct3D_Input_Parameters{
  public:
  static const SpaceType IP_SpaceDimension = ThreeD; 

  // Input file name:
  char Input_File_Name[INPUT_PARAMETER_LENGTH_RECONSTRUCT3D];

  // Input file stream:
  ifstream Input_File;

  // Input file line number:
  int Line_Number;
  int Message_Number;

  // Parameters characteristic for the method involved for solving the 
  // linear system of equations
  char Method_Used[INPUT_PARAMETER_LENGTH_RECONSTRUCT3D];
  int  Method;
  char Geometric_Weighting[INPUT_PARAMETER_LENGTH_RECONSTRUCT3D];
  int  geom_weighting;
  char Data_Dependent_Weighting[INPUT_PARAMETER_LENGTH_RECONSTRUCT3D];
  int  data_depend_weighting;
  char Pseudo_Inverse[INPUT_PARAMETER_LENGTH_RECONSTRUCT3D];
  
  char Reduce_Order[INPUT_PARAMETER_LENGTH_RECONSTRUCT3D];

  // Function file name to be reconstructed and related input parameters:
  char Function_Type[INPUT_PARAMETER_LENGTH_RECONSTRUCT3D];
  int Reconstruction_Order;
  int i_Function;
  TestFunction3D TestF;
  IntegralTestFunction3D IntTestF;
  char Integration_Type[INPUT_PARAMETER_LENGTH_RECONSTRUCT3D];

  // Initial condition type indicator and related input parameters:
  char ICs_Type[INPUT_PARAMETER_LENGTH_RECONSTRUCT3D];
  int i_ICs;

  // Grid type indicator and related input parameters:
  char Grid_Type[INPUT_PARAMETER_LENGTH_RECONSTRUCT3D];
  char NACA_Aerofoil_Type[INPUT_PARAMETER_LENGTH_RECONSTRUCT3D];
  int i_Grid;
  char CellNumber_or_DeltaCell[2];
  int  NCells_Idir, NCells_Jdir, NCells_Kdir,
       NBlk_Idir, NBlk_Jdir, NBlk_Kdir;
  int Nghost_Cells;
  double Delta_X_Cell;
  double Delta_Y_Cell;
  double Delta_Z_Cell;

  double Box_Width, Box_Height, Box_Length, Plate_Length, 
         Pipe_Length, Pipe_Radius, 
         Blunt_Body_Radius, Blunt_Body_Mach_Number,
         Grain_Length, Grain_Radius, Grain_To_Throat_Length,
         Nozzle_Length, Nozzle_Radius_Exit, Nozzle_Radius_Throat,
         Cylinder_Radius, Ellipse_Length_X_Axis, 
         Ellipse_Length_Y_Axis, Chord_Length, Orifice_Radius,
         Wedge_Angle, Wedge_Length;

  double X_Scale, X_Rotate;
  Vector3D X_Shift;
  char Stretching_Function_I[INPUT_PARAMETER_LENGTH_RECONSTRUCT3D]; /* Grid stretching function */
  char Stretching_Function_J[INPUT_PARAMETER_LENGTH_RECONSTRUCT3D]; /* Grid stretching function */
  char Stretching_Function_K[INPUT_PARAMETER_LENGTH_RECONSTRUCT3D]; /* Grid stretching function */
  int Stretch_I, Stretch_J, Stretch_K;
  double Beta_I, Beta_J, Beta_K, Tau_I, Tau_J, Tau_K; /* Grid stretching parameters */
  int NumOfIter_UnsmoothMesh;
  double CharacteristicLength;
  double CutoffKnobHighOrderReconstruction;
  double CENO_Cutoff;

  int Stretching_Type_Idir, Stretching_Type_Jdir, Stretching_Type_Kdir;
  double Stretching_Factor_Idir, Stretching_Factor_Jdir, Stretching_Factor_Kdir;

  char **ICEMCFD_FileNames;

  //Subgrid parameters 
  char SubGridPoints_or_DeltaSubGrid[2];
  int  NSubGrid_Points_Idir;
  int  NSubGrid_Points_Jdir;
  int  NSubGrid_Points_Kdir;
  double Delta_X_of_SubGrid;
  double Delta_Y_of_SubGrid;
  double Delta_Z_of_SubGrid;

  double X_min, X_Length;
  double Y_min, Y_Length;
  double Z_min, Z_Length;

  // Limiter type indicator and related input parameters:
  char Limiter_Type[INPUT_PARAMETER_LENGTH_RECONSTRUCT3D];
  int i_Limiter;

  // Output file name:
  char Output_File_Name[INPUT_PARAMETER_LENGTH_RECONSTRUCT3D];

  // Multi-block mesh definition input file names:
  char Grid_File_Name[INPUT_PARAMETER_LENGTH_RECONSTRUCT3D];
  char Grid_Definition_File_Name[INPUT_PARAMETER_LENGTH_RECONSTRUCT3D];

  // Gnuplot file name:
  char Gnuplot_File_Name[INPUT_PARAMETER_LENGTH_RECONSTRUCT3D];

  // Next_Control_Parameter:
  char Next_Control_Parameter[INPUT_PARAMETER_LENGTH_RECONSTRUCT3D];

  // Output format type indicator:
  char Output_Format_Type[INPUT_PARAMETER_LENGTH_RECONSTRUCT3D];
  int i_Output_Format;

  // Access functions
  const int & NumSubGridI(void) const { return NSubGrid_Points_Idir; }
  const int & NumSubGridJ(void) const { return NSubGrid_Points_Jdir; }
  const int & NumSubGridK(void) const { return NSubGrid_Points_Kdir; }
  const int & RecOrder(void) const { return Reconstruction_Order; }
  const int & iCell(void) const { return NCells_Idir;}
  const int & jCell(void) const { return NCells_Jdir;}
  const int & kCell(void) const { return NCells_Kdir;}
  const int & Nghost(void) const { return Nghost_Cells;}
  const TestFunction3D & TestFunction(void) const {return TestF;}
  double & CutoffKnob(void) { return CutoffKnobHighOrderReconstruction; }
  const double & CutoffKnob(void) const { return CutoffKnobHighOrderReconstruction; }
  int & ReconstructionMethod(void) {return Method;}
  const int & ReconstructionMethod(void) const {return Method;}
  const int & Limiter(void) const {return i_Limiter;}
  double & FitTolerance(void) {return CENO_Cutoff;}
  const double & FitTolerance(void) const {return CENO_Cutoff;}

  /* Input-output operators. */

  friend ostream &operator << (ostream &out_file, const Reconstruct3D_Input_Parameters &IP);
  friend istream &operator >> (istream &in_file, Reconstruct3D_Input_Parameters &IP);

};

/*************************************************************
 * Reconstruct3D_Input_Parameters -- External subroutines.         *
 *************************************************************/

extern void Open_Input_File(Reconstruct3D_Input_Parameters &IP);

extern void Close_Input_File(Reconstruct3D_Input_Parameters &IP);

extern void Set_Default_Input_Parameters(Reconstruct3D_Input_Parameters &IP);

extern void Get_Next_Input_Control_Parameter(Reconstruct3D_Input_Parameters &IP);

extern int Parse_Next_Input_Control_Parameter(Reconstruct3D_Input_Parameters &IP);

extern int Process_Input_Control_Parameter_File(Reconstruct3D_Input_Parameters
						&Input_Parameters,
                                                char *Input_File_Name_ptr,
                                                int &Command_Flag);

#endif /* _RECONSTRUCT3D_INPUT_INCLUDED  */

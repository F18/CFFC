/*!\file HO_Grid2DQuadMultiBlock_InputForTesting.h
  \brief Header file defining 2D Multi-block Quadrilateral Grid Input Class for testing purposes. */

#ifndef _GRID2DQUADMULTIBLOCK_INPUTFORTESTING_INCLUDED
#define _GRID2DQUADMULTIBLOCK_INPUTFORTESTING_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../../ICEM/ICEMCFD.h"    // Include ICEMCFD input header file.
#include "../../Utilities/TypeDefinition.h" // Include TypeDefinition header file.

/* Define the structures and classes. */

#define	INPUT_PARAMETER_LENGTH_GRID2DTESTING    128
#define	LONG_INPUT_PARAMETER_LENGTH_GRID2DTESTING    256

// Enviroment flag for CFFC root directory path
#define PATHVAR_GRID2DTESTING "CFFC_Path"

/*! 
 * \class Grid2DTesting_Input_Parameters
 *
 * @brief Definition and manipulation of grid 2D input variables.
 *        This class is used only for testing purposed.
 */
class Grid2DTesting_Input_Parameters{
private:
public:
  //! @name Input file parameters.
  //@{
  //! CFFC root directory path:
  char CFFC_Path[INPUT_PARAMETER_LENGTH_GRID2DTESTING];
  //! Input file name:
  char Input_File_Name[LONG_INPUT_PARAMETER_LENGTH_GRID2DTESTING];
  //! Input file stream:
  ifstream Input_File;
  //! Input file line number:
  int Line_Number;
  //@}

  //! @name Reconstruction type indicator and related input parameters:
  //@{
  char Reconstruction_Type[INPUT_PARAMETER_LENGTH_GRID2DTESTING]; /*!< Spatial reconstruction type. */
  int i_Reconstruction;		/*!< Index to store the reconstruction type. */
  int i_ReconstructionMethod;	/*!< Index to store the reconstruction method. */
  int Space_Accuracy;		/*!< Parameter to show the order of accuracy in space. */
  int IncludeHighOrderBoundariesRepresentation;	/*!< Flag for including or excluding high-order BCs. */
  char Viscous_Reconstruction_Type[INPUT_PARAMETER_LENGTH_GRID2DTESTING];
  int i_Viscous_Reconstruction; /*!< Index to store the method for computing face gradients for viscous fluxes. */
  //@}

  //! @name Flow geometry (planar or axisymmetric):
  //@{
  char Flow_Geometry_Type[INPUT_PARAMETER_LENGTH_GRID2DTESTING];
  int Axisymmetric;
  //@}

  //! @name Grid type indicator and related input parameters:
  //@{
  char Grid_Type[INPUT_PARAMETER_LENGTH_GRID2DTESTING];
  char NACA_Aerofoil_Type[INPUT_PARAMETER_LENGTH_GRID2DTESTING];
  int i_Grid;
  int Number_of_Cells_Idir, Number_of_Cells_Jdir, Number_of_Ghost_Cells,
      Number_of_Blocks_Idir, Number_of_Blocks_Jdir;
  double Box_Width, Box_Height, Plate_Length, 
         Pipe_Length, Pipe_Radius, 
         Blunt_Body_Radius, Blunt_Body_Mach_Number,
         Chamber_Length, Chamber_Radius, Chamber_To_Throat_Length,
         Nozzle_Length, Nozzle_Radius_Exit, Nozzle_Radius_Throat, Grain_Radius,
         Cylinder_Radius, Cylinder_Radius2, Ellipse_Length_X_Axis, 
         Ellipse_Length_Y_Axis, Chord_Length, Orifice_Radius,
         Annulus_Theta_Start, Annulus_Theta_End;
  int Nozzle_Type;
  Vector2D VertexSW, VertexSE, VertexNE, VertexNW;
  double X_Scale, X_Rotate;
  Vector2D X_Shift;
  char **ICEMCFD_FileNames;
  int IterationsOfInteriorNodesDisturbances; /*<! Number of iterations run to move (disturb) the interior nodes.
						 (create an unsmooth interior mesh). */
  int Num_Of_Spline_Control_Points;  /*<! Number of points used to define the spline curve. */
  //@}

  //! @name Mesh stretching factor.
  //@{
  int i_Mesh_Stretching;
  int Mesh_Stretching_Type_Idir;
  int Mesh_Stretching_Type_Jdir;
  double Mesh_Stretching_Factor_Idir;
  double Mesh_Stretching_Factor_Jdir;
  //@}

  //! @name Boundary conditions:
  //@{
  char Boundary_Conditions_Specified[INPUT_PARAMETER_LENGTH_GRID2DTESTING];
  int BCs_Specified;
  char BC_North_Type[INPUT_PARAMETER_LENGTH_GRID2DTESTING];
  char BC_South_Type[INPUT_PARAMETER_LENGTH_GRID2DTESTING];
  char BC_East_Type[INPUT_PARAMETER_LENGTH_GRID2DTESTING];
  char BC_West_Type[INPUT_PARAMETER_LENGTH_GRID2DTESTING];
  int BC_North, BC_South, BC_East, BC_West;
  //@}

  //! @name AMR input parameters:
  //@{
  //! Unsteady AMR flag.
  int AMR;
  //! Unsteady AMR frequency.
  int AMR_Frequency;
  //! \brief Unsteady AMR global reference --> set the time reference for unsteady AMR 
  //
  //!  --> Set to OFF, the reference is the start of the current simulation/restart <br>
  //!  --> Set to ON,  the reference is the absolute beginning of the simulation
  int AMR_Global_Reference;
  //! Number of initial mesh refinements.
  int Number_of_Initial_Mesh_Refinements;
  //! Number of uniform mesh refinements.
  int Number_of_Uniform_Mesh_Refinements;
  //! Number of boundary mesh refinements.
  int Number_of_Boundary_Mesh_Refinements;
  //! Maximum number of refinement levels.
  int Maximum_Refinement_Level;
  //! Minimum number of refinement levels.
  int Minimum_Refinement_Level;
  //! Threshold for refinement of mesh.
  double Threshold_for_Refinement;
  //! Threshold for coarsening of mesh.
  double Threshold_for_Coarsening;
  //! Smooth quad block flag:
  int i_Smooth_Quad_Block;
  //@}

  //! @name Output parameters:
  //@{
  //! Output file name:
  char Output_File_Name[LONG_INPUT_PARAMETER_LENGTH_GRID2DTESTING];
  //! Multi-block mesh definition input file names:
  char Grid_File_Name[LONG_INPUT_PARAMETER_LENGTH_GRID2DTESTING];
  char Grid_Definition_File_Name[LONG_INPUT_PARAMETER_LENGTH_GRID2DTESTING];
  //! Restart file name:
  char Restart_File_Name[LONG_INPUT_PARAMETER_LENGTH_GRID2DTESTING];
  //! Gnuplot file name:
  char Gnuplot_File_Name[LONG_INPUT_PARAMETER_LENGTH_GRID2DTESTING];
  //! Next_Control_Parameter:
  char Next_Control_Parameter[INPUT_PARAMETER_LENGTH_GRID2DTESTING];
  //! Output format type indicator:
  char Output_Format_Type[INPUT_PARAMETER_LENGTH_GRID2DTESTING];
  int i_Output_Format;
  int Restart_Solution_Save_Frequency;
  //! Output progress frequency:
  int Output_Progress_Frequency;
  //! Batch mode or verbose
  short verbose_flag;
  //@}

  //! @name Multi-block solution-adaption and parallel domain decomposition input parameters:
  //@{
  int Number_of_Processors, Number_of_Blocks_Per_Processor;
  //@}

  //! @name Default Constructor
  //@{
  Grid2DTesting_Input_Parameters(void);
  //@}

  //! @name Destructor
  //@{
  ~Grid2DTesting_Input_Parameters(void);
  //@}

  //! @name Obtain the CFFC root directory path:
  //@{
  void get_cffc_path();
  //@}

  //! @name Reconstruction related member functions:
  //@{
  int ReconstructionOrder(void) {return (Space_Accuracy-1);} //!< return order of reconstruction based on Space_Accuracy
  /*!
   * \brief Return the number of required ghost cells for the specified ReconstructionMethod
   */
  int Nghost(void) const;
  //@}

  //! @name Access fields:
  //@{
  short & Verbose(void) {return verbose_flag;}
  const short & Verbose(void) const {return verbose_flag;}
  void Verbose(const int & batch_flag){ (batch_flag != 0) ? verbose_flag=OFF: verbose_flag=ON; }
  bool OutputBoundaryReferenceState(const int & BCtype) const;
  //@}

  //! @name Operating functions:
  //@{
  int Parse_Input_File(char *Input_File_Name_ptr); //!< \brief Parse input file
  void Get_Next_Input_Control_Parameter(void);    //!< \brief Read the next input control parameter
  void Open_Input_File(void);
  void Close_Input_File(void);
  void Set_Default_Input_Parameters(void);
  int Parse_Next_Input_Control_Parameter(void);
  int Check_Input_Parameters(void);
  int Process_Input_Control_Parameter_File(char *Input_File_Name_ptr,
					   int &Command_Flag);
  //@}

  //! @name Input-output operators:
  //@{
  friend ostream &operator << (ostream &out_file,
		               const Grid2DTesting_Input_Parameters &IP);
  friend istream &operator >> (istream &in_file,
			       Grid2DTesting_Input_Parameters &IP);
  //@}

};

/************************************************************************
 * Grid2DTesting_Input_Parameters::Grid2DTesting_Input_Parameters() *
 * -->  Default Constructor                                             *
 ***********************************************************************/
inline Grid2DTesting_Input_Parameters::Grid2DTesting_Input_Parameters(void): verbose_flag(OFF){
  ICEMCFD_FileNames = NULL;
}

/************************************************************************
 * Grid2DTesting_Input_Parameters::~Grid2DTesting_Input_Parameters()*
 * -->  Default Constructor                                             *
 ***********************************************************************/
inline Grid2DTesting_Input_Parameters::~Grid2DTesting_Input_Parameters(void){
  if (ICEMCFD_FileNames != NULL){
    for (int i = 0 ; i < 3 ; ++i){
      delete [] ICEMCFD_FileNames[i]; ICEMCFD_FileNames[i] = NULL;
    }
    delete [] ICEMCFD_FileNames; ICEMCFD_FileNames = NULL;
  }
}

/*********************************************************************
 * Grid2DTesting_Input_Parameters::get_cffc_path -- Get CFFC path. *
 *********************************************************************/
inline void Grid2DTesting_Input_Parameters::get_cffc_path(){
  char *string_ptr;
  // Check to see if environment varible exists.
  if (getenv(PATHVAR_GRID2DTESTING) == NULL) {
    //Set default path
     string_ptr = "CFFC";
     strcpy(CFFC_Path, string_ptr);
  } else {
     //Set path specified by environment variable
     strcpy(CFFC_Path, getenv(PATHVAR_GRID2DTESTING));
  }
}

#endif /* _GRID2DQUADMULTIBLOCK_INPUTFORTESTING_INCLUDED  */

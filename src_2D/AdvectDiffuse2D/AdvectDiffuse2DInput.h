/*!\file AdvectDiffuse2DInput.h
  \brief Header file defining 2D Advection Diffusion Equation Input Parameter Class. */

#ifndef _ADVECTDIFFUSE2D_INPUT_INCLUDED
#define _ADVECTDIFFUSE2D_INPUT_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "AdvectDiffuse2DState.h" // Include 2D advection diffusion equation solution state header file
#include "../FASMultigrid2D/FASMultigrid2DInput.h" // Include multigrid input header file.
#include "../ICEM/ICEMCFD.h"      // Include ICEMCFD input header file.
#include "../Utilities/TypeDefinition.h" // Include TypeDefinition header file.
#include "AdvectDiffuse2DParameterFields.h" /* Include 2D advection diffusion parameter fields */
#include "../NewtonKrylovSchwarz2D/NKSInput2D.h" /* Include file for NKS */
#include "../HighOrderReconstruction/HighOrder2D_Input.h" /* Include file for high-order */

/* Define the structures and classes. */

#define	INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D    128
#define	LONG_INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D    256

// Enviroment flag for CFFC root directory path
#define PATHVAR_ADVECTDIFFUSE2D "CFFC_Path"

/*! 
 * \class AdvectDiffuse2D_Input_Parameters
 *
 * @brief Definition and manipulation of 2D Advection-Diffusion input variables.
 *
 */
class AdvectDiffuse2D_Input_Parameters{
private:
public:
  //! @name Input file parameters.
  //@{
  //! CFFC root directory path:
  char CFFC_Path[INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
  //! Input file name:
  char Input_File_Name[LONG_INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
  //! Input file stream:
  ifstream Input_File;
  //! Input file line number:
  int Line_Number;
  //@}

  //! @name Time integration type indicator and related input parameters:
  //@{
  char Time_Integration_Type[INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
  int i_Time_Integration;
  int Time_Accurate, Local_Time_Stepping, 
      Maximum_Number_of_Time_Steps, N_Stage;
  double CFL_Number, Time_Max;
  //! Residual variable:
  int i_Residual_Variable;
  //@}

  //! @name Implicit residual smoothing control parameters:
  //@{
  int Residual_Smoothing;
  double Residual_Smoothing_Epsilon;
  int Residual_Smoothing_Gauss_Seidel_Iterations;
  //@}

  //! @name Multigrid related input parameters:
  //@{
  Multigrid_Input_Parameters Multigrid_IP;
  //@}

  //@{ @name Newton-Krylov-Schwarz related input parameters:
  NKS_Input_Parameters  NKS_IP;
  int Solver_Type;
  //@}

  //! @name Reconstruction type indicator and related input parameters:
  //@{
  char Reconstruction_Type[INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D]; /*!< Spatial reconstruction type. */
  int i_Reconstruction;		/*!< Index to store the reconstruction type. */
  int i_ReconstructionMethod;	/*!< Index to store the reconstruction method. */
  int Space_Accuracy;		/*!< Parameter to show the order of accuracy in space. */
  int IncludeHighOrderBoundariesRepresentation;	/*!< Flag for including or excluding high-order BCs. */
  char Viscous_Reconstruction_Type[INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
  int i_Viscous_Reconstruction; /*!< Index to store the method for computing face gradients for viscous fluxes. */
  //@}

  //! @name Limiter type indicator and related input parameters:
  //@{
  char Limiter_Type[INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
  int i_Limiter;
  int  Freeze_Limiter;
  double Freeze_Limiter_Residual_Level;
  //@}

  //! @name Flux function type:
  //@{
  int i_Flux_Function;
  //@}

  //! @name Initial condition type indicator and related input parameters:
  //@{
  char ICs_Type[INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
  int i_ICs;
  AdvectDiffuse2D_State Uo, U1, U2;
  AdvectDiffuse2D_ExactSolutions *ExactSoln; /*!< Pointer to the exact solution */
  AdvectDiffuse2D_State RefU;		/*!< Reference state, used by CENO to normalize the
					   variables in the computation of the smoothness indicator. */
  unsigned int Exact_Integration_Digits;    //!< Number of exact digits with which the some integrations are carried out
  //@}

  //! @name Source term parameters:
  //@{
  SourceTermFields *SourceTerm;        /*!< Pointer to the source term field. */
  //@}

  //! @name Flow geometry (planar or axisymmetric):
  //@{
  char Flow_Geometry_Type[INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
  int Axisymmetric;
  //@}

  //! @name Grid type indicator and related input parameters:
  //@{
  char Grid_Type[INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
  char NACA_Aerofoil_Type[INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
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
  char Boundary_Conditions_Specified[INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
  int BCs_Specified;
  char BC_North_Type[INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
  char BC_South_Type[INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
  char BC_East_Type[INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
  char BC_West_Type[INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
  int BC_North, BC_South, BC_East, BC_West;
  //! Reference states for north and south boundary conditons
  AdvectDiffuse2D_State Ref_State_BC_North, Ref_State_BC_South;
  //! Reference states for east and west boundary conditons
  AdvectDiffuse2D_State Ref_State_BC_East, Ref_State_BC_West; 
  AdvectDiffuse2D_InflowField *Inflow;    /*!< Pointer to the inflow field */
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
  char Output_File_Name[LONG_INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
  //! Multi-block mesh definition input file names:
  char Grid_File_Name[LONG_INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
  char Grid_Definition_File_Name[LONG_INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
  //! Restart file name:
  char Restart_File_Name[LONG_INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
  //! Gnuplot file name:
  char Gnuplot_File_Name[LONG_INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
  //! Next_Control_Parameter:
  char Next_Control_Parameter[INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
  //! Output format type indicator:
  char Output_Format_Type[INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
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
  AdvectDiffuse2D_Input_Parameters(void);
  //@}

  //! @name Destructor
  //@{
  ~AdvectDiffuse2D_Input_Parameters(void);
  //@}

  //! @name Obtain the CFFC root directory path:
  //@{
  void get_cffc_path();
  //@}

  //! @name Reconstruction related member functions:
  //@{
  /*! Return order of reconstruction based on Space_Accuracy.
    To obtain a certain global space accuracy (i.e. convective, diffusive and source term)
    a piecewise polynomial reconstruction of the same order as the space accuracy must be performed. */
  int ReconstructionOrder(void) {return (Space_Accuracy);}
  int & Limiter(void) {return i_Limiter;}                    //!< write/read selected limiter
  const int & Limiter(void) const {return i_Limiter;}        //!< return selected limiter (read only)
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

  //! @name Accuracy assessment parameters:
  //@{
  unsigned int Accuracy_Assessment_Mode;
  unsigned int Accuracy_Assessment_Exact_Digits;
  unsigned int Accuracy_Assessment_Parameter;
  //@}

  //! @name Operating functions:
  //@{
  int Parse_Input_File(char *Input_File_Name_ptr); //!< \brief Parse input file
  void Get_Next_Input_Control_Parameter(void);    //!< \brief Read the next input control parameter
  //@}

  //! @name Input-output operators:
  //@{
  friend ostream &operator << (ostream &out_file,
		               const AdvectDiffuse2D_Input_Parameters &IP);
  friend istream &operator >> (istream &in_file,
			       AdvectDiffuse2D_Input_Parameters &IP);
  //@}

};

/************************************************************************
 * AdvectDiffuse2D_Input_Parameters::AdvectDiffuse2D_Input_Parameters() *
 * -->  Default Constructor                                             *
 ***********************************************************************/
inline AdvectDiffuse2D_Input_Parameters::AdvectDiffuse2D_Input_Parameters(void): verbose_flag(OFF){
  ICEMCFD_FileNames = NULL;

  // Get access to the SourceTermField
  SourceTerm = &SourceTermFields::getInstance();

  // Get access to the AdvectDiffuse2D_ExactSolutions object
  ExactSoln = &AdvectDiffuse2D_ExactSolutions::getInstance();

  // Get access to the AdvectDiffuse2D_InflowField object
  Inflow = &AdvectDiffuse2D_InflowField::getInstance();
}

/************************************************************************
 * AdvectDiffuse2D_Input_Parameters::~AdvectDiffuse2D_Input_Parameters()*
 * -->  Default Constructor                                             *
 ***********************************************************************/
inline AdvectDiffuse2D_Input_Parameters::~AdvectDiffuse2D_Input_Parameters(void){
  if (ICEMCFD_FileNames != NULL){
    for (int i = 0 ; i < 3 ; ++i){
      delete [] ICEMCFD_FileNames[i]; ICEMCFD_FileNames[i] = NULL;
    }
    delete [] ICEMCFD_FileNames; ICEMCFD_FileNames = NULL;
  }
}

/*********************************************************************
 * AdvectDiffuse2D_Input_Parameters::get_cffc_path -- Get CFFC path. *
 *********************************************************************/
inline void AdvectDiffuse2D_Input_Parameters::get_cffc_path(){
  char *string_ptr;
  // Check to see if environment varible exists.
  if (getenv(PATHVAR_ADVECTDIFFUSE2D) == NULL) {
    //Set default path
     string_ptr = "CFFC";
     strcpy(CFFC_Path, string_ptr);
  } else {
     //Set path specified by environment variable
     strcpy(CFFC_Path, getenv(PATHVAR_ADVECTDIFFUSE2D));
  }
}

/*************************************************************
 * AdvectDiffuse2D_Input_Parameters -- External subroutines. *
 *************************************************************/

extern void Open_Input_File(AdvectDiffuse2D_Input_Parameters &IP);

extern void Close_Input_File(AdvectDiffuse2D_Input_Parameters &IP);

extern void Set_Default_Input_Parameters(AdvectDiffuse2D_Input_Parameters &IP);

extern void Broadcast_Input_Parameters(AdvectDiffuse2D_Input_Parameters &IP);

#ifdef _MPI_VERSION
extern void Broadcast_Input_Parameters(AdvectDiffuse2D_Input_Parameters &IP,
                                       MPI::Intracomm &Communicator, 
                                       const int Source_CPU);
#endif

extern void Get_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters &IP);

extern int Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters &IP);

extern int Check_Input_Parameters(AdvectDiffuse2D_Input_Parameters &IP);

extern int Process_Input_Control_Parameter_File(AdvectDiffuse2D_Input_Parameters &Input_Parameters,
                                                char *Input_File_Name_ptr,
                                                int &Command_Flag);

#endif /* _ADVECTDIFFUSE2D_INPUT_INCLUDED  */

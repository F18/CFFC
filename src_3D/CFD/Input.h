/* Input.h:  Header file defining general CFD 
             solution input parameter class. */

#ifndef _INPUT_INCLUDED
#define _INPUT_INCLUDED

/* Include required CFFC header files. */

#ifndef _CFD_INCLUDED
#include "../CFD/CFD.h"
#endif // _CFD_INCLUDED

#ifndef _MPI_INCLUDED
#include "../MPI/MPI.h"
#endif // _MPI_INCLUDED

#ifndef _GRID3D_INPUT_INCLUDED
#include "../Grid/Grid3DInput.h"
#endif // _GRID3D_INPUT_INCLUDED

/* Define the structures and classes. */

#define	INPUT_PARAMETER_LENGTH    128

//Enviroment Flag 
#define PATHVAR "CFFC_Path"

/********************************************************
 * Class:  Input_Parameters                             *
 ********************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
class Input_Parameters {
  private:
  public:
   
  char CFFC_Path[INPUT_PARAMETER_LENGTH];
 
  // Input file name:
  char Input_File_Name[INPUT_PARAMETER_LENGTH];

  // Input file stream:
  ifstream Input_File;

  // Input file line number, the index of solution parameters:
  int Line_Number;
  int Solution_Parameters_Index;
  
  // grid input parameters
  Grid3D_Input_Parameters IP_Grid;
  
  // Time integration type indicator and related input parameters:
  char Time_Integration_Type[INPUT_PARAMETER_LENGTH];
  int i_Time_Integration;
  int Time_Accurate, Local_Time_Stepping, 
      Maximum_Number_of_Time_Steps, N_Stage;
  int p_Norms_Specified_Parameter; // 0 density, 1, momentum, 2 k,...
  double CFL_Number, Time_Max;

  // Implicit residual smoothing control parameters:
  int Residual_Smoothing;
  double Residual_Smoothing_Epsilon;
  int Residual_Smoothing_Gauss_Seidel_Iterations;

  // Reconstruction type indicator and related input parameters:
  char Reconstruction_Type[INPUT_PARAMETER_LENGTH];
  int i_Reconstruction;

  // Limiter type indicator and related input parameters:
  char Limiter_Type[INPUT_PARAMETER_LENGTH];
  int i_Limiter;
  int  Freeze_Limiter;
  double Freeze_Limiter_Residual_Level;

  // Flux_Function_Type and related input parameters
  char Flux_Function_Type[INPUT_PARAMETER_LENGTH];
  int i_Flux_Function;

  // Initial condition type indicator and related input parameters:
  char ICs_Type[INPUT_PARAMETER_LENGTH];
  int i_ICs;
  int i_Grid_Level;
  
  //Morton Ordering Parameters
  int Morton;
  int Morton_Reordering_Frequency;

  /**********  SPECIFIC ***************/
  string react_name;
  char React_Name[INPUT_PARAMETER_LENGTH];
  char **Multispecies;
  string *multispecies; //array of species names
  double *mass_fractions; //array of mass fractions
  int num_species;
 
  
  //Pa, K, m/s, degress from horizontal
  double Pressure, Temperature, Flow_Angle;
  double Mach_Number_Reference,Mach_Number, Re_lid;
  int Preconditioning;

  //BC Moving wall velocity
  Vector3D Moving_wall_velocity;
  Vector3D Pressure_Gradient;
 
  // Flow type indicator and related input parameters:
  char Flow_Type[INPUT_PARAMETER_LENGTH];
  int i_Flow_Type;

  // Flow geometry (planar or axisymmetric):
  char Flow_Geometry_Type[INPUT_PARAMETER_LENGTH];
  int Axisymmetric; //0 no, 1 yes
  int Gravity;      //0 no, 1 yes
  double Global_Schmidt;  //depricated, use each individual Schmidt's
  double *Schmidt;  //individual for each species
  int Wall_Boundary_Treatments; 
  //0, 1,2 , automatic, wall function, low_Reynolds number
 
  int Stretch_Level; //1, 2, higher stretching
  double Reynolds_Number;
  double Kinematic_Viscosity_Wall;
  double Eddy_Viscosity_Limit_Coefficient;
  
  //Debug Level 0 for none, 1,2,3... level of verboseness
  int debug_level;

  double Box_Length, Box_Width, Box_Height, Plate_Length, 
         Pipe_Length, Pipe_Radius, 
         Blunt_Body_Radius, Blunt_Body_Mach_Number,
         Grain_Length, Grain_Radius, Grain_To_Throat_Length,
         Nozzle_Length, Nozzle_Radius_Exit, Nozzle_Radius_Throat,
         Cylinder_Radius, Ellipse_Length_X_Axis, 
         Ellipse_Length_Y_Axis, Chord_Length, Orifice_Radius,
         Wedge_Angle, Wedge_Length;
  double Length_Shroud,Radius_Shroud,Length_BluffBody,Radius_BluffBody,
         Radius_Orifice, Length_Inlet_Pipe, Radius_Inlet_Pipe, Length_Combustor_Tube, 
         Radius_Combustor_Tube;
  double BluffBody_Coflow_Air_Velocity, BluffBody_Coflow_Fuel_Velocity;
  int BluffBody_Data_Usage; // 0 no, 1 yes, 
  
  double X_Scale, X_Rotate;
  
  // AMR
  int AMR, AMR_Frequency;

  // Number of initial mesh refinements.
  int Number_of_Initial_Mesh_Refinements;

  // Number of uniform mesh refinements.
  int Number_of_Uniform_Mesh_Refinements;

  // Number of boundary mesh refinements.
  int Number_of_Boundary_Mesh_Refinements;

  // Maximum and minimum number of refinement levels.
  int Maximum_Refinement_Level, Minimum_Refinement_Level;

  // Thresholds for refinement and coarsening of mesh.
  double Threshold_for_Refinement, Threshold_for_Coarsening;

  // Output file name:
  char Output_File_Name[INPUT_PARAMETER_LENGTH];

  // Multi-block mesh definition input file names:
  char Grid_File_Name[INPUT_PARAMETER_LENGTH];
  char Grid_Definition_File_Name[INPUT_PARAMETER_LENGTH];

  // Restart file name:
  char Restart_File_Name[INPUT_PARAMETER_LENGTH];

  // Gnuplot file name:
  char Gnuplot_File_Name[INPUT_PARAMETER_LENGTH];

  // Next_Control_Parameter:
  char Next_Control_Parameter[INPUT_PARAMETER_LENGTH];

  // Output format type indicator:
  char Output_Format_Type[INPUT_PARAMETER_LENGTH];
  int i_Output_Format;
  int Restart_Solution_Save_Frequency;
  int Time_Accurate_Plot_Freq;

  // Multi-block solution-adaption and parallel domain
  // decomposition input parameters

  int Number_of_Processors, Number_of_Blocks_Per_Processor;

  void Allocate();
  void Deallocate();
  
  SOLN_pSTATE Wo;
  SOLN_cSTATE Uo;
 
  // constructor ...set some default parameters
  Input_Parameters(void){
   
     Multispecies=NULL;
     multispecies=NULL;
     mass_fractions=NULL;
     
     Set_Default_Input_Parameters();
  }

  // member functions
  void Get_CFFC_Path(void);
  void Open_Input_File(void);
  void Close_Input_File(void);
  void Set_Default_Input_Parameters(void);
  void Broadcast_Input_Parameters(void);
  void Get_Next_Input_Control_Parameter(void);
  int Parse_Next_Input_Control_Parameter(void);
  int Process_Input_Control_Parameter_File(char *Input_File_Name_ptr,
                                           int &Command_Flag);

};

/* Input-output operators. */
template<class SOLN_pSTATE, class SOLN_cSTATE>
istream &operator >> (istream &in_file,
                             Input_Parameters<SOLN_pSTATE, 
                             SOLN_cSTATE> &IP);

template<class SOLN_pSTATE, class SOLN_cSTATE>
ostream &operator << (ostream &out_file,
                               const Input_Parameters<SOLN_pSTATE,
                               SOLN_cSTATE> &IP);

/*************************************************************
 * Input_Parameters -- Memory Management                     *
 *************************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
void Input_Parameters<SOLN_pSTATE, SOLN_cSTATE>::Allocate() {
   Multispecies = new char*[num_species];
   for (int i = 0; i < num_species; i++) {
      Multispecies[i] = new char[INPUT_PARAMETER_LENGTH];
   } 
   multispecies = new string[num_species]; 
   mass_fractions = new double[num_species];
   Schmidt = new double[num_species];
}

template<class SOLN_pSTATE, class SOLN_cSTATE>
void Input_Parameters<SOLN_pSTATE, SOLN_cSTATE>::Deallocate() {
   if(Multispecies != NULL ){
      for (int i = 0; i < num_species; i++) {
         delete[] Multispecies[i]; Multispecies[i]=NULL;
      }
      delete[] Multispecies; Multispecies=NULL;
   }
  
   if(multispecies != NULL){
      delete[] multispecies; multispecies=NULL;
   }
   
   if(mass_fractions != NULL){
      delete[] mass_fractions; mass_fractions=NULL;
      if( Schmidt != NULL) delete[] Schmidt; Schmidt = NULL;
   }  
   
} 

/********************************************************
 * Routine: Get_CFFC_Path                               *
 *                                                      *
 * Obtains the location of CFFC source directory        *
 *                                                      *
 ********************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
void Input_Parameters<SOLN_pSTATE, SOLN_cSTATE>::Get_CFFC_Path(void){
  char *string_ptr;
 
  // Check to see if environment varible exists.
  if (getenv(PATHVAR) == NULL) {
    //Set default path
     string_ptr = "CFFC";
     strcpy(CFFC_Path, string_ptr);
  } else {
     //Set path specified by environment variable
     strcpy(CFFC_Path, getenv(PATHVAR));
  }

}

/********************************************************
 * Routine: Open_Input_File                             *
 *                                                      *
 * Opens the appropriate input data file.               *
 *                                                      *
 ********************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
void Input_Parameters<SOLN_pSTATE, SOLN_cSTATE>::Open_Input_File(void) {

    Input_File.open(Input_File_Name, ios::in);
    if (!Input_File.bad()) {
       Line_Number = 0;
       Input_File.setf(ios::skipws);
    } /* endif */

}

/********************************************************
 * Routine: Close_Input_File                            *
 *                                                      *
 * Closes the appropriate input data file.              *
 *                                                      *
 ********************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
void Input_Parameters<SOLN_pSTATE, SOLN_cSTATE>::Close_Input_File(void) {

    Input_File.unsetf(ios::skipws);
    Input_File.close();

}

/********************************************************
 * Routine: Set_Default_Input_Parameters                *
 *                                                      *
 * Assigns default values to the input parameters.      *
 *                                                      *
 ********************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
void Input_Parameters<SOLN_pSTATE, SOLN_cSTATE>::
Set_Default_Input_Parameters(void) {

   
   int i;
   char *string_ptr;
   
   /* CFFC Directory Path */
   Get_CFFC_Path();

   string_ptr = "Euler3D.in";
   strcpy(Input_File_Name, string_ptr);
 
   string_ptr = "Explicit_Euler";
   strcpy(Time_Integration_Type, string_ptr);
   i_Time_Integration = TIME_STEPPING_EXPLICIT_EULER;
   Time_Accurate = 0;
   Local_Time_Stepping = SCALAR_LOCAL_TIME_STEPPING;
   Preconditioning = 0; //default off
   Maximum_Number_of_Time_Steps = 100;
   N_Stage = 1;
   Time_Max = ZERO;
   p_Norms_Specified_Parameter = 1;
   Residual_Smoothing = 0;
   Residual_Smoothing_Epsilon = ZERO;
   Residual_Smoothing_Gauss_Seidel_Iterations = 2;
   
   string_ptr = "Barth_Jespersen";
   strcpy(Limiter_Type, string_ptr);
   i_Limiter = LIMITER_BARTH_JESPERSEN;
   
   string_ptr = "HLLE";
   strcpy(Flux_Function_Type, string_ptr);
   i_Flux_Function = FLUX_FUNCTION_HLLE;
   
   string_ptr = "Uniform";
   strcpy(ICs_Type, string_ptr);
   i_ICs = IC_UNIFORM;
   
   debug_level = 0;  //default no debug information
   Mach_Number = ZERO;
   Mach_Number_Reference = ONE;  //for precondtioning
   Flow_Angle = ZERO; 
   Global_Schmidt = ONE;
   Reynolds_Number = 500000.0;
   Kinematic_Viscosity_Wall = 1.4590E-5;
   i_Grid_Level = 0;
   
   react_name ="NO_REACTIONS";   
   //Use air with 79% N2, and 21% 02 by volume.(ie. mol)
   num_species = 2;
   Allocate();
   multispecies[0] = "N2"; 
   multispecies[1] = "O2"; 
   mass_fractions[0] = 0.765; 
   mass_fractions[1] = 0.235;
   Schmidt[0] = Global_Schmidt;
   Schmidt[1] = Global_Schmidt;
    
   Moving_wall_velocity.zero(); 
   Re_lid = 100.0;
   Pressure_Gradient.zero(); 
   string_ptr = "Inviscid";
   strcpy(Flow_Type, string_ptr);
   i_Flow_Type = FLOWTYPE_INVISCID;
   
   string_ptr = "Planar";
   strcpy(Flow_Geometry_Type, string_ptr);
   Axisymmetric = 0;
   Gravity = 0;  //default sans gravity
   Wall_Boundary_Treatments = 0; 
   Stretch_Level = 1; 
   BluffBody_Data_Usage = 0; 
     
   Plate_Length = ONE;
   Pipe_Length = ONE;
   Pipe_Radius = 0.1;
   Blunt_Body_Radius = ONE;
   Blunt_Body_Mach_Number = TWO;
   Grain_Length = 0.835;
   Grain_Radius = 0.020;
   Grain_To_Throat_Length = 0.05;
   Nozzle_Length = 0.150;
   Nozzle_Radius_Exit = 0.030;
   Nozzle_Radius_Throat = 0.010;
   Cylinder_Radius = ONE;
   Ellipse_Length_X_Axis = TWO;
   Ellipse_Length_Y_Axis = HALF;
   Chord_Length = ONE;
   Orifice_Radius = ONE;
   Wedge_Angle = 25.0;
   Wedge_Length = HALF;

   X_Scale = ONE;
   X_Rotate = ZERO;
   Length_Shroud = 0.1;
   Radius_Shroud = 0.043;
   Length_BluffBody = 0.04;
   Radius_BluffBody = 0.02 ;
   Radius_Orifice = 0.001;
   Radius_Inlet_Pipe = 0.0508;
   Radius_Combustor_Tube = 0.0762;
   Length_Inlet_Pipe = 0.127;
   Length_Combustor_Tube = 0.508;
    
   BluffBody_Coflow_Air_Velocity = 20.0;
   BluffBody_Coflow_Fuel_Velocity = 61.0; 
   
   Morton = 0;
   Morton_Reordering_Frequency = 1000;
   
   AMR = 0;
   AMR_Frequency = 100;
   Number_of_Initial_Mesh_Refinements = 0;
   Number_of_Uniform_Mesh_Refinements = 0;
   Number_of_Boundary_Mesh_Refinements = 0;
   Maximum_Refinement_Level = 100;
   Minimum_Refinement_Level = 1;

   Threshold_for_Refinement = 0.50;
   Threshold_for_Coarsening = 0.10;

   string_ptr = "outputfile.dat";
   strcpy(Output_File_Name, string_ptr);
   
   string_ptr = "gridfile.grid";
   strcpy(Grid_File_Name, string_ptr);
   string_ptr = "gridfile.griddef";
   strcpy(Grid_Definition_File_Name, string_ptr);

   string_ptr = "restartfile.soln";
   strcpy(Restart_File_Name, string_ptr);
   
   string_ptr = "gnuplotfile.gplt";
   strcpy(Gnuplot_File_Name, string_ptr);

   string_ptr = "Tecplot";
   strcpy(Output_Format_Type, string_ptr);
   i_Output_Format = IO_TECPLOT;
   Restart_Solution_Save_Frequency = 1000;
   
   Time_Accurate_Plot_Freq = 0;
   string_ptr = " ";
   strcpy(Next_Control_Parameter, string_ptr);
   
   Line_Number = 0;
   
   Number_of_Processors = CFFC_MPI::Number_of_Processors;
   //  Number_of_Processors = 1;
   Number_of_Blocks_Per_Processor = 100;      

   Freeze_Limiter = 0;
   Freeze_Limiter_Residual_Level = 1e-4;
   
   //define multispecies with no reactions.
   react_name ="NO_REACTIONS";   
   //Use air with 79% N2, and 21% 02 by volume.(ie. mol)
   num_species = 2;
   Allocate();
   multispecies[0] = "N2"; 
   multispecies[1] = "O2"; 
   mass_fractions[0] = 0.765; 
   mass_fractions[1] = 0.235;
   Schmidt[0] = Global_Schmidt;
   Schmidt[1] = Global_Schmidt;
   
   Wo.React.set_reactions(react_name);
   Wo.React.set_species(multispecies,num_species);
 
   
   //Get Species parameters and set default initial values
   SOLN_pSTATE::set_species_data(Wo,
                                 num_species,
                                 multispecies,
                                 CFFC_Path,
                                 debug_level,
                                 Mach_Number_Reference,
                                 Schmidt);
   SOLN_cSTATE::set_species_data(Uo,
                                 num_species,
                                 multispecies,
                                 CFFC_Path,
                                 debug_level,
                                 Mach_Number_Reference,
                                 Schmidt);

   //Air at STD_ATM
   Pressure = Wo.p;
   Temperature = Wo.T(); 
   Wo.v.x =Mach_Number*Wo.a()*cos(TWO*PI*Flow_Angle/360.00);
   Wo.v.y =Mach_Number*Wo.a()*sin(TWO*PI*Flow_Angle/360.00);
   Wo.v.z =Mach_Number*Wo.a()*sin(TWO*PI*Flow_Angle/360.00);
   
   Wo.set_initial_values(mass_fractions);
   Uo.set_initial_values(mass_fractions);
 
   Uo = Wo.U();  
 
     
}

/********************************************************
 * Routine: Get_Next_Input_Control_Parameter            *
 *                                                      *
 * Get the next input control parameter from the input  *
 * file.                                                *
 *                                                      *
 ********************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
   void Input_Parameters<SOLN_pSTATE, SOLN_cSTATE>::
   Get_Next_Input_Control_Parameter(void) {

    int i;
    char buffer[256];

    Line_Number = Line_Number + 1;
    Input_File.getline(buffer, sizeof(buffer));
    i = 0;
    if (buffer[0] != '#') {
 
       while (1) {
          if (buffer[i] == ' ' || buffer[i] == '=' ) break;
          i = i + 1;
          if (i > strlen(buffer) ) break;
       } /* endwhile */
       buffer[i] = '\0';
    } /* endif */
    strcpy(Next_Control_Parameter, buffer);

    //    cout<<"\n "<<Next_Control_Parameter<<endl; cout.flush();

}

/********************************************************
 * Routine: Parse_Next_Input_Control_Parameter          *
 *                                                      *
 * Parses and executes the next input control parameter *
 * from the input file.                                 *
 *                                                      *
 ********************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
int Input_Parameters<SOLN_pSTATE, SOLN_cSTATE>::Parse_Next_Input_Control_Parameter(void) {

    int i_command;
    char buffer[256];
   
    i_command = 0;
   
    if (strcmp(Next_Control_Parameter, "CFFC_Path") == 0) {
      i_command = 1;
      Get_Next_Input_Control_Parameter();
      strcpy(CFFC_Path, Next_Control_Parameter);

    } else if (strcmp(Next_Control_Parameter, "Time_Integration_Type") == 0) {
      i_command = 2;
      Get_Next_Input_Control_Parameter();
      strcpy(Time_Integration_Type, Next_Control_Parameter);
       if (strcmp(Time_Integration_Type, "Explicit_Euler") == 0) {
           i_Time_Integration = TIME_STEPPING_EXPLICIT_EULER;
           N_Stage = 1;
       } else if (strcmp(Time_Integration_Type, "Explicit_Predictor_Corrector") == 0) {
           i_Time_Integration = TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR;
           N_Stage = 2;
       } else if (strcmp(Time_Integration_Type, "Explicit_Runge_Kutta") == 0) {
           i_Time_Integration = TIME_STEPPING_EXPLICIT_RUNGE_KUTTA;
           N_Stage = 4;
       } else if (strcmp(Time_Integration_Type, "Multistage_Optimal_Smoothing") == 0) {
           i_Time_Integration = TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING;
           N_Stage = 4;
       } else {
          i_Time_Integration = TIME_STEPPING_EXPLICIT_EULER;
          N_Stage = 1;
       } /* endif */
       
    } else if (strcmp(Next_Control_Parameter, "Reconstruction_Type") == 0) {
       i_command = 3;
       Get_Next_Input_Control_Parameter();
       strcpy(Reconstruction_Type, 
              Next_Control_Parameter);
       if (strcmp(Reconstruction_Type, "Green_Gauss") == 0) {
          i_Reconstruction = RECONSTRUCTION_GREEN_GAUSS;
       } else if (strcmp(Reconstruction_Type, "Least_Squares") == 0) {
          i_Reconstruction = RECONSTRUCTION_LEAST_SQUARES;
       } else {
          i_Reconstruction = RECONSTRUCTION_GREEN_GAUSS;
       } /* endif */

    } else if (strcmp(Next_Control_Parameter, "Limiter_Type") == 0) {
       i_command = 4;
       Get_Next_Input_Control_Parameter();
       strcpy(Limiter_Type, 
              Next_Control_Parameter);
       if (strcmp(Limiter_Type, "One") == 0) {
          i_Limiter = LIMITER_ONE;
       } else if (strcmp(Limiter_Type, "Zero") == 0) {
          i_Limiter = LIMITER_ZERO;
       } else if (strcmp(Limiter_Type, "VanLeer") == 0) {
          i_Limiter = LIMITER_VANLEER;
       } else if (strcmp(Limiter_Type, "VanAlbada") == 0) {
          i_Limiter = LIMITER_VANALBADA;
       } else if (strcmp(Limiter_Type, "Barth_Jespersen") == 0) {
          i_Limiter = LIMITER_BARTH_JESPERSEN;
       } else if (strcmp(Limiter_Type, "Venkatakrishnan") == 0) {
          i_Limiter = LIMITER_VENKATAKRISHNAN;
       } else {
          i_Limiter = LIMITER_VANLEER ;
       } /* endif */

    } else if (strcmp(Next_Control_Parameter, "Flux_Function_Type") == 0) {
       i_command = 5;
       Get_Next_Input_Control_Parameter();
       strcpy(Flux_Function_Type, 
              Next_Control_Parameter);
       if (strcmp(Flux_Function_Type, "Godunov") == 0) {
          i_Flux_Function = FLUX_FUNCTION_GODUNOV;
       } else if (strcmp(Flux_Function_Type, "Roe") == 0) {
          i_Flux_Function = FLUX_FUNCTION_ROE;
       } else if (strcmp(Flux_Function_Type, "Rusanov") == 0) {
          i_Flux_Function = FLUX_FUNCTION_RUSANOV;
       } else if (strcmp(Flux_Function_Type, "HLLE") == 0) {
          i_Flux_Function = FLUX_FUNCTION_HLLE;
       } else if (strcmp(Flux_Function_Type, "Linde") == 0) {
          i_Flux_Function = FLUX_FUNCTION_LINDE;
       } else if (strcmp(Flux_Function_Type, "HLLC") == 0) {
          i_Flux_Function = FLUX_FUNCTION_HLLC;
       } else {
          i_Flux_Function = FLUX_FUNCTION_ROE;
       } /* endif */

    } else if (strcmp(Next_Control_Parameter, "ICs_Type") == 0) {
       i_command = 6;
       Get_Next_Input_Control_Parameter();
       strcpy(ICs_Type, 
              Next_Control_Parameter);
       if (strcmp(ICs_Type, "Constant") == 0) {
          i_ICs = IC_CONSTANT;
       } else if (strcmp(ICs_Type, "Uniform") == 0) {
          i_ICs = IC_UNIFORM;
       } else if (strcmp(ICs_Type, "Sod") == 0) {
          i_ICs = IC_SOD;
       } else if (strcmp(ICs_Type, "Sod_Xdir") == 0) {
          i_ICs = IC_SOD_XDIR;
       } else if (strcmp(ICs_Type, "Sod_Ydir") == 0) {
          i_ICs = IC_SOD_YDIR;
       } else if (strcmp(ICs_Type, "Groth") == 0) {
          i_ICs = IC_GROTH;
       } else if (strcmp(ICs_Type, "Groth_Xdir") == 0) {
          i_ICs = IC_GROTH_XDIR;
       } else if (strcmp(ICs_Type, "Groth_Ydir") == 0) {
          i_ICs = IC_GROTH_YDIR;
       } else if (strcmp(ICs_Type, "Einfeldt") == 0) {
          i_ICs = IC_EINFELDT;
       } else if (strcmp(ICs_Type, "Einfeldt_Xdir") == 0) {
          i_ICs = IC_EINFELDT_XDIR;
       } else if (strcmp(ICs_Type, "Einfeldt_Ydir") == 0) {
          i_ICs = IC_EINFELDT_YDIR;
       } else if (strcmp(ICs_Type, "Shock_Wave_Xdir") == 0) {
          i_ICs = IC_SHOCK_WAVE_XDIR;
       } else if (strcmp(ICs_Type, "Shock_Wave_Ydir") == 0) {
          i_ICs = IC_SHOCK_WAVE_YDIR;
       } else if (strcmp(ICs_Type, "Contact_Surface_Xdir") == 0) {
          i_ICs = IC_CONTACT_SURFACE_XDIR;
       } else if (strcmp(ICs_Type, "Contact_Surface_Ydir") == 0) {
          i_ICs = IC_CONTACT_SURFACE_YDIR;
       } else if (strcmp(ICs_Type, "Rarefaction_Wave_Xdir") == 0) {
          i_ICs = IC_RAREFACTION_WAVE_XDIR;
       } else if (strcmp(ICs_Type, "Rarefaction_Wave_Ydir") == 0) {
          i_ICs = IC_RAREFACTION_WAVE_YDIR;
       } else if (strcmp(ICs_Type, "ShockBox") == 0) {
          i_ICs = IC_SHOCK_BOX;
       } else if (strcmp(ICs_Type, "High_Pressure_Reservoir") == 0) {
          i_ICs = IC_HIGH_PRESSURE_RESERVOIR;
       } else if (strcmp(ICs_Type, "Low_Pressure_Reservoir") == 0) {
          i_ICs = IC_LOW_PRESSURE_RESERVOIR;
       } else if (strcmp(ICs_Type, "Riemann") == 0) {
          i_ICs = IC_RIEMANN;
       } else if (strcmp(ICs_Type, "Riemann_Xdir") == 0) {
          i_ICs = IC_RIEMANN_XDIR;
       } else if (strcmp(ICs_Type, "Riemann_Ydir") == 0) {
	 i_ICs = IC_RIEMANN_YDIR;  
       } else if (strcmp(ICs_Type, "Wedge_Flow") == 0) {
	 i_ICs = IC_WEDGE_FLOW;	 
       } else if (strcmp(ICs_Type, "Mix") == 0) {
	 i_ICs = IC_GAS_MIX;
       } else if (strcmp(ICs_Type, "Core_Flame") == 0 ){
	 i_ICs = IC_CHEM_CORE_FLAME ;
       } else if (strcmp(ICs_Type, "Inverse_Flame") == 0 ){
	 i_ICs = IC_CHEM_INVERSE_FLAME ; 
       } else if (strcmp(ICs_Type, "Pressure_Gradient_x") == 0 ){
	 i_ICs = IC_PRESSURE_GRADIENT_X;
       } else if (strcmp(ICs_Type, "Pressure_Gradient_y") == 0 ){
	 i_ICs = IC_PRESSURE_GRADIENT_Y;
       } else if (strcmp(ICs_Type, "Pressure_Gradient_z") == 0 ){
	 i_ICs = IC_PRESSURE_GRADIENT_Z;
       } else if (strcmp(ICs_Type, "Couette") == 0 ){
          i_ICs = IC_VISCOUS_COUETTE; 
       } else if (strcmp(ICs_Type, "Couette_Pressure_Gradient_x") == 0 ){
          i_ICs = IC_VISCOUS_COUETTE_PRESSURE_GRADIENT_X;
       } else if (strcmp(ICs_Type, "Couette_Pressure_Gradient_y") == 0 ){
          i_ICs = IC_VISCOUS_COUETTE_PRESSURE_GRADIENT_Y;
       } else if (strcmp(ICs_Type, "Couette_Pressure_Gradient_z") == 0 ){
          i_ICs = IC_VISCOUS_COUETTE_PRESSURE_GRADIENT_Z;
       } else if (strcmp(ICs_Type, "1DPremixedFlame") == 0 ){
          i_ICs = IC_CHEM_1DFLAME;  
       } else if (strcmp(ICs_Type, "Pipe_Flow") == 0) {
          i_ICs = IC_TURBULENT_PIPE_FLOW;
       } else if (strcmp(ICs_Type, "Coflow") == 0) {
          i_ICs = IC_TURBULENT_COFLOW;
       }else if (strcmp(ICs_Type, "Driven_Cavity_Flow") == 0) {
          i_ICs = IC_VISCOUS_DRIVEN_CAVITY_FLOW;
       } else if (strcmp(ICs_Type, "Turbulent_Dump_Combustor") == 0) {
          i_ICs = IC_TURBULENT_DUMP_COMBUSTOR;
       }else if (strcmp(ICs_Type, "Laminar_Channel_Flow") == 0) {
          i_ICs = IC_CHANNEL_FLOW;
       } else if (strcmp(ICs_Type, "Turbulent_Channel_Flow") == 0) {
          i_ICs = IC_CHANNEL_FLOW;
       }else if (strcmp(ICs_Type, "Turbulent_Diffusion_Flame") == 0) {
          i_ICs = IC_TURBULENT_DIFFUSION_FLAME;
       }else if (strcmp(ICs_Type, "Turbulent_Free_Jet_Flame") == 0) {
          i_ICs = IC_FREE_JET_FLAME;
       }else if (strcmp(ICs_Type, "Restart") == 0) {
          i_ICs = IC_RESTART;
       } else {
          i_ICs = IC_UNIFORM;
       } /* endif */

   } else if (strcmp(Next_Control_Parameter, "Grid_Type") == 0) {
       i_command = 7;
       Get_Next_Input_Control_Parameter();
       strcpy(IP_Grid.Grid_Type, 
              Next_Control_Parameter);
       if (strcmp(IP_Grid.Grid_Type, "Cartesian") == 0) {
          IP_Grid.i_Grid = GRID_CARTESIAN_UNIFORM;
          IP_Grid.Box_Length = ONE;
          IP_Grid.Box_Width = ONE;
          IP_Grid.Box_Height = ONE;
       } else if (strcmp(IP_Grid.Grid_Type, "Cube") == 0) {
          IP_Grid.i_Grid = GRID_CUBE;
          IP_Grid.Box_Length = ONE;
          IP_Grid.Box_Width = ONE;
          IP_Grid.Box_Height = ONE;
       } else if (strcmp(IP_Grid.Grid_Type, "Channel") == 0) {
          IP_Grid.i_Grid = GRID_CHANNEL;
          IP_Grid.Box_Length = 0.2;
          IP_Grid.Box_Width = 0.001;
          IP_Grid.Box_Height = 0.001;
       } else if (strcmp(IP_Grid.Grid_Type, "Turbulent_Channel") == 0) {
          IP_Grid.i_Grid = GRID_CHANNEL;
          if (Pressure_Gradient.x !=ZERO) {
	    //IP_Grid.geometry_index = 1;
             IP_Grid.Box_Length = 1.524;
             IP_Grid.Box_Width  = 0.127;
             IP_Grid.Box_Height = 0.127;
          }
          if (Pressure_Gradient.y !=ZERO) {
	    //IP_Grid.geometry_index = 2;
             IP_Grid.Box_Length = 0.127;
             IP_Grid.Box_Width  = 1.524;
             IP_Grid.Box_Height = 0.127;
          }
          if (Pressure_Gradient.z !=ZERO) {
	    //IP_Grid.geometry_index = 3;
             IP_Grid.Box_Length =  0.127;
             IP_Grid.Box_Width  =  0.127;
             IP_Grid.Box_Height =  1.524;
          }
       } else if (strcmp(IP_Grid.Grid_Type, "Couette") == 0) {
          IP_Grid.i_Grid = GRID_COUETTE;
          if (Moving_wall_velocity.x !=ZERO) {
	    //IP_Grid.geometry_index = 1;
             IP_Grid.Box_Length = 0.2;
             IP_Grid.Box_Width  = 0.001;
             IP_Grid.Box_Height = 0.001;
          }
          if (Moving_wall_velocity.y !=ZERO) {
            // IP_Grid.geometry_index = 2;
             IP_Grid.Box_Length = 0.001;
             IP_Grid.Box_Width  = 0.2;
             IP_Grid.Box_Height = 0.001;
          }
          if (Moving_wall_velocity.z !=ZERO) {
	    //IP_Grid.geometry_index = 3;
             IP_Grid.Box_Length = 0.001;
             IP_Grid.Box_Width  = 0.001;
             IP_Grid.Box_Height = 0.2;
          }
       }else {
          IP_Grid.i_Grid = GRID_CUBE;
          IP_Grid.Box_Length = ONE;
          IP_Grid.Box_Width = ONE;
          IP_Grid.Box_Height = ONE;
       } /* endif */

   } else if (strcmp(Next_Control_Parameter, "Output_File_Name") == 0) {
      i_command = 8;
      Get_Next_Input_Control_Parameter();
      strcpy(Output_File_Name, 
             Next_Control_Parameter);
      strcat(Output_File_Name, ".dat");
      strcpy(Grid_File_Name, 
             Next_Control_Parameter);
      strcat(Grid_File_Name, ".grid");
      strcpy(Grid_Definition_File_Name, 
             Next_Control_Parameter);
      strcat(Grid_Definition_File_Name, ".griddef");
      strcpy(Restart_File_Name, 
             Next_Control_Parameter);
      strcat(Restart_File_Name, ".soln");
      strcpy(Gnuplot_File_Name, 
             Next_Control_Parameter);
      strcat(Gnuplot_File_Name, ".gplt");
      
    } else if (strcmp(Next_Control_Parameter, "Grid_File_Name") == 0) {
       i_command = 9;
       Get_Next_Input_Control_Parameter( );
       strcpy(Grid_File_Name, 
              Next_Control_Parameter);
       strcat(Grid_File_Name, ".grid");
       strcpy(Grid_Definition_File_Name, 
              Next_Control_Parameter);
       strcat(Grid_Definition_File_Name, ".griddef");

    } else if (strcmp(Next_Control_Parameter, "Restart_File_Name") == 0) {
       i_command = 10;
       Get_Next_Input_Control_Parameter( );
       strcpy(Restart_File_Name, 
              Next_Control_Parameter);
       strcat(Restart_File_Name, ".soln");

    } else if (strcmp(Next_Control_Parameter, "Number_of_Cells_Idir") == 0) {
       i_command = 11;
       Line_Number = Line_Number + 1;
       Input_File >> IP_Grid.NCells_Idir;
       Input_File.getline(buffer, sizeof(buffer));
       if ( IP_Grid.NCells_Idir < 1) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(Next_Control_Parameter, "Number_of_Cells_Jdir") == 0) {
       i_command = 12;
       Line_Number = Line_Number + 1;
       Input_File >> IP_Grid.NCells_Jdir;
       Input_File.getline(buffer, sizeof(buffer));
       if (IP_Grid.NCells_Jdir <1) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(Next_Control_Parameter, "Number_of_Cells_Kdir") == 0) {
       i_command = 13;
       Line_Number = Line_Number + 1;
       Input_File >> IP_Grid.NCells_Kdir;
       Input_File.getline(buffer, sizeof(buffer));
       if (IP_Grid.NCells_Kdir < 1) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(Next_Control_Parameter, "Number_of_Blocks_Idir") == 0) {
       i_command = 14;
       Line_Number = Line_Number + 1;
       Input_File >> IP_Grid.NBlk_Idir;
       Input_File.getline(buffer, sizeof(buffer));
       if ( IP_Grid.NBlk_Idir < 1) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(Next_Control_Parameter, "Number_of_Blocks_Jdir") == 0) {
       i_command = 15;
       Line_Number = Line_Number + 1;
       Input_File >> IP_Grid.NBlk_Jdir;
       Input_File.getline(buffer, sizeof(buffer));
       if ( IP_Grid.NBlk_Jdir < 1) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(Next_Control_Parameter, "Number_of_Blocks_Kdir") == 0) {
       i_command = 16;
       Line_Number = Line_Number + 1;
       Input_File >> IP_Grid.NBlk_Kdir;
       Input_File.getline(buffer, sizeof(buffer));
       if ( IP_Grid.NBlk_Kdir < 1) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(Next_Control_Parameter, "Stretching_Factor_Idir") == 0) {
       i_command = 17;
       Line_Number = Line_Number + 1;
       Input_File >> IP_Grid.Stretching_Factor_Idir;
       Input_File.getline(buffer, sizeof(buffer));
       if ( IP_Grid.Stretching_Factor_Idir < ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(Next_Control_Parameter, "Stretching_Factor_Jdir") == 0) {
       i_command = 18;
       Line_Number = Line_Number + 1;
       Input_File >> IP_Grid.Stretching_Factor_Jdir;
       Input_File.getline(buffer, sizeof(buffer));
       if ( IP_Grid.Stretching_Factor_Jdir < ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(Next_Control_Parameter, "Stretching_Factor_Kdir") == 0) {
       i_command = 19;
       Line_Number = Line_Number + 1;
       Input_File >> IP_Grid.Stretching_Factor_Kdir;
       Input_File.getline(buffer, sizeof(buffer));
       if ( IP_Grid.Stretching_Factor_Kdir < ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(Next_Control_Parameter, "Stretching_Type_Idir") == 0) {
       i_command = 20;
       Line_Number = Line_Number + 1;
       Input_File >> IP_Grid.Stretching_Type_Idir;
       Input_File.getline(buffer, sizeof(buffer));
       if ( IP_Grid.Stretching_Type_Idir < 0) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(Next_Control_Parameter, "Stretching_Type_Jdir") == 0) {
       i_command = 21;
       Line_Number = Line_Number + 1;
       Input_File >> IP_Grid.Stretching_Type_Jdir;
       Input_File.getline(buffer, sizeof(buffer));
       if ( IP_Grid.Stretching_Type_Jdir < 0) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(Next_Control_Parameter, "Stretching_Type_Kdir") == 0) {
       i_command = 22;
       Line_Number = Line_Number + 1;
       Input_File >> IP_Grid.Stretching_Type_Kdir;
       Input_File.getline(buffer, sizeof(buffer));
       if ( IP_Grid.Stretching_Type_Kdir < 0) i_command = INVALID_INPUT_VALUE;

    }else if (strcmp(Next_Control_Parameter, "Time_Accurate") == 0) {
       i_command = 23;
       Line_Number = Line_Number + 1;
       Input_File >> Time_Accurate;
       Input_File.getline(buffer, sizeof(buffer));
       if (Time_Accurate != 0 &&
           Time_Accurate != 1) Time_Accurate = 0;
       if (Time_Accurate) {
          Local_Time_Stepping = GLOBAL_TIME_STEPPING;
       } /* endif */

    } else if (strcmp(Next_Control_Parameter, "Local_Time_Stepping") == 0) {
       i_command = 24;
       Line_Number = Line_Number + 1;
       Input_File >> Local_Time_Stepping;
       Input_File.getline(buffer, sizeof(buffer));
       if (Local_Time_Stepping != GLOBAL_TIME_STEPPING &&
           Local_Time_Stepping != SCALAR_LOCAL_TIME_STEPPING &&
	   Local_Time_Stepping != MATRIX_LOCAL_TIME_STEPPING &&	   
           Local_Time_Stepping != LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER &&
	   Local_Time_Stepping != SEMI_IMPLICIT_LOCAL_TIME_STEPPING &&
	   Local_Time_Stepping != SEMI_IMPLICIT_LOW_MACH_NUMBER_PRECONDITIONER) {
         Local_Time_Stepping = SCALAR_LOCAL_TIME_STEPPING;
	}
       if(Local_Time_Stepping == LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER ||
       	  Local_Time_Stepping == SEMI_IMPLICIT_LOW_MACH_NUMBER_PRECONDITIONER){
	 Preconditioning = 1;
       }       

    } else if (strcmp(Next_Control_Parameter, "Maximum_Number_of_Time_Steps") == 0) {
       i_command = 25;
       Line_Number = Line_Number + 1;
       Input_File >> Maximum_Number_of_Time_Steps;
       Input_File.getline(buffer, sizeof(buffer));
       if (Maximum_Number_of_Time_Steps < 0) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(Next_Control_Parameter, "N_Stage") == 0) {
       i_command = 26;
       Line_Number = Line_Number + 1;
       Input_File >> N_Stage;
       Input_File.getline(buffer, sizeof(buffer));
       if (N_Stage < 0) i_command = INVALID_INPUT_VALUE;

    }else if (strcmp(Next_Control_Parameter, "p_norms") == 0) {
       i_command = 27;
       Line_Number = Line_Number + 1;
       Input_File >> p_Norms_Specified_Parameter;
       Input_File.getline(buffer, sizeof(buffer));
       if (N_Stage < 0) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(Next_Control_Parameter, "CFL_Number") == 0) {
       i_command = 28;
       Line_Number = Line_Number + 1;
       Input_File >> CFL_Number;
       Input_File.getline(buffer, sizeof(buffer));
       if (CFL_Number <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(Next_Control_Parameter, "Box_Width") == 0) {
       i_command = 29;
       Line_Number = Line_Number + 1;
       Input_File >> Box_Width;
       Input_File.getline(buffer, sizeof(buffer));
       if (Box_Width <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(Next_Control_Parameter, "Box_Height") == 0) {
       i_command = 30;
       Line_Number = Line_Number + 1;
       Input_File >> Box_Height;
       Input_File.getline(buffer, sizeof(buffer));
       if (Box_Height <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(Next_Control_Parameter, "Plate_Length") == 0) {
       i_command = 31;
       Line_Number = Line_Number + 1;
       Input_File >> Plate_Length;
       Input_File.getline(buffer, sizeof(buffer));
       if (Plate_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(Next_Control_Parameter, "Pipe_Length") == 0) {
       i_command = 32;
       Line_Number = Line_Number + 1;
       Input_File >> Pipe_Length;
       Input_File.getline(buffer, sizeof(buffer));
       if (Pipe_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(Next_Control_Parameter, "Pipe_Radius") == 0) {
       i_command = 33;
       Line_Number = Line_Number + 1;
       Input_File >> Pipe_Radius;
       Input_File.getline(buffer, sizeof(buffer));
       if (Pipe_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(Next_Control_Parameter, "Blunt_Body_Radius") == 0) {
       i_command = 34;
       Line_Number = Line_Number + 1;
       Input_File >> Blunt_Body_Radius;
       Input_File.getline(buffer, sizeof(buffer));
       if (Blunt_Body_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(Next_Control_Parameter, "Blunt_Body_Mach_Number") == 0) {
       i_command = 35;
       Line_Number = Line_Number + 1;
       Input_File >> Blunt_Body_Mach_Number;
       Input_File.getline(buffer, sizeof(buffer));
       if (Blunt_Body_Mach_Number <= ONE) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(Next_Control_Parameter, "Cylinder_Radius") == 0) {
       i_command = 36;
       Line_Number = Line_Number + 1;
       Input_File >> Cylinder_Radius;
       Input_File.getline(buffer, sizeof(buffer));
       if (Cylinder_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(Next_Control_Parameter, "Ellipse_Length_X_Axis") == 0) {
       i_command = 37;
       Line_Number = Line_Number + 1;
       Input_File >> Ellipse_Length_X_Axis;
       Input_File.getline(buffer, sizeof(buffer));
       if (Ellipse_Length_X_Axis <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(Next_Control_Parameter, "Ellipse_Length_Y_Axis") == 0) {
       i_command = 38;
       Line_Number = Line_Number + 1;
       Input_File >> Ellipse_Length_Y_Axis;
       Input_File.getline(buffer, sizeof(buffer));
       if (Ellipse_Length_Y_Axis <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(Next_Control_Parameter, "Chord_Length") == 0) {
       i_command = 39;
       Line_Number = Line_Number + 1;
       Input_File >> Chord_Length;
       Input_File.getline(buffer, sizeof(buffer));
       if (Chord_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(Next_Control_Parameter, "Orifice_Radius") == 0) {
       i_command = 41;
       Line_Number = Line_Number + 1;
       Input_File >> Orifice_Radius;
       Input_File.getline(buffer, sizeof(buffer));
       if (Orifice_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(Next_Control_Parameter, "Wedge_Angle") == 0) {
       i_command = 42;
       Line_Number = Line_Number + 1;
       Input_File >> Wedge_Angle;
       Input_File.getline(buffer, sizeof(buffer));
       if (Wedge_Angle <= ZERO) i_command = INVALID_INPUT_VALUE;
    
    } else if (strcmp(Next_Control_Parameter, "Wedge_Length") == 0) {
       i_command = 43;
       Line_Number = Line_Number + 1;
       Input_File >> Wedge_Length;
       Input_File.getline(buffer, sizeof(buffer));
       if (Wedge_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(Next_Control_Parameter, "Grain_Length") == 0) {
       i_command = 44;
       Line_Number = Line_Number + 1;
       Input_File >> Grain_Length;
       Input_File.getline(buffer, sizeof(buffer));
       if (Grain_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(Next_Control_Parameter, "Grain_Radius") == 0) {
       i_command = 45;
       Line_Number = Line_Number + 1;
       Input_File >> Grain_Radius;
       Input_File.getline(buffer, sizeof(buffer));
       if (Grain_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(Next_Control_Parameter, "Grain_To_Throat_Length") == 0) {
       i_command = 46;
       Line_Number = Line_Number + 1;
       Input_File >> Grain_To_Throat_Length;
       Input_File.getline(buffer, sizeof(buffer));
       if (Grain_To_Throat_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(Next_Control_Parameter, "Nozzle_Length") == 0) {
       i_command = 47;
       Line_Number = Line_Number + 1;
       Input_File >> Nozzle_Length;
       Input_File.getline(buffer, sizeof(buffer));
       if (Nozzle_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(Next_Control_Parameter, "Nozzle_Radius_Exit") == 0) {
       i_command = 48;
       Line_Number = Line_Number + 1;
       Input_File >> Nozzle_Radius_Exit;
       Input_File.getline(buffer, sizeof(buffer));
       if (Nozzle_Radius_Exit <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(Next_Control_Parameter, "Nozzle_Radius_Throat") == 0) {
       i_command = 49;
       Line_Number = Line_Number + 1;
       Input_File >> Nozzle_Radius_Throat;
       Input_File.getline(buffer, sizeof(buffer));
       if (Nozzle_Radius_Throat <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(Next_Control_Parameter, "Length_BluffBody") == 0) {
       i_command = 50;
       Line_Number = Line_Number + 1;
       Input_File >> Length_BluffBody;
       Input_File.getline(buffer, sizeof(buffer));
       if (Length_BluffBody <ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(Next_Control_Parameter, "Radius_Orifice") == 0) {
       i_command = 51 ;
       Line_Number = Line_Number + 1;
       Input_File >> Radius_Orifice;
       Input_File.getline(buffer, sizeof(buffer));
       if (Radius_Orifice <ZERO) i_command = INVALID_INPUT_VALUE;   

    }else if (strcmp(Next_Control_Parameter, "BluffBody_Coflow_Air_Velocity") == 0) {
       i_command = 52;
       Line_Number = Line_Number + 1;
       Input_File >> BluffBody_Coflow_Air_Velocity;
       Input_File.getline(buffer, sizeof(buffer));
       if (BluffBody_Coflow_Air_Velocity <ZERO) i_command = INVALID_INPUT_VALUE;
       
    }else if (strcmp(Next_Control_Parameter, "BluffBody_Coflow_Fuel_Velocity") == 0) {
       i_command = 53;
       Line_Number = Line_Number + 1;
       Input_File >> BluffBody_Coflow_Fuel_Velocity;
       Input_File.getline(buffer, sizeof(buffer));
       if (BluffBody_Coflow_Fuel_Velocity <ZERO) i_command = INVALID_INPUT_VALUE;       
       //Get next line and read in Schmidt numbers else will use defaults
       Get_Next_Input_Control_Parameter( );
            
    }  else if (strcmp(Next_Control_Parameter, "Time_Max") == 0) {
       i_command = 54;
       Line_Number = Line_Number + 1;
       Input_File >> Time_Max;
       Input_File.getline(buffer, sizeof(buffer));
       Time_Max = Time_Max/THOUSAND;
       if (Time_Max < ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(Next_Control_Parameter, "Number_of_Blocks_Per_Processor") == 0) {
       i_command = 55;
       Line_Number = Line_Number + 1;
       Input_File >> Number_of_Blocks_Per_Processor;
       Input_File.getline(buffer, sizeof(buffer));
       if (Number_of_Blocks_Per_Processor < 1) i_command = INVALID_INPUT_VALUE;
 
    } else if (strcmp(Next_Control_Parameter, "Output_Format_Type") == 0) {
       i_command = 56;
       Get_Next_Input_Control_Parameter( );
       strcpy(Output_Format_Type, 
              Next_Control_Parameter);
       if (strcmp(Output_Format_Type, "Gnuplot") == 0  ||
           strcmp(Output_Format_Type, "gnuplot") == 0  ||
           strcmp(Output_Format_Type, "GNUPLOT") == 0) {
          i_Output_Format = IO_GNUPLOT;
       } else if (strcmp(Output_Format_Type, "Tecplot") == 0  ||
                  strcmp(Output_Format_Type, "tecplot") == 0  ||
                  strcmp(Output_Format_Type, "TECPLOT") == 0) {
          i_Output_Format = IO_TECPLOT;
       } else if (strcmp(Output_Format_Type, "Matlab") == 0  ||
                  strcmp(Output_Format_Type, "matlab") == 0  ||
                  strcmp(Output_Format_Type, "MATLAB") == 0) {
          i_Output_Format = IO_MATLAB;
       } else if (strcmp(Output_Format_Type, "Octave") == 0  ||
                  strcmp(Output_Format_Type, "octave") == 0  ||
                  strcmp(Output_Format_Type, "OCTAVE") == 0) {
          i_Output_Format = IO_OCTAVE;
       } else {
          i_Output_Format = IO_TECPLOT;
       } /* endif */

    } else if (strcmp (Next_Control_Parameter, "Morton") == 0) {
      i_command = 57;
      Line_Number = Line_Number + 1;
      Input_File >> Morton;
      Input_File.getline(buffer, sizeof(buffer));
      if (Morton < ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(Next_Control_Parameter, "Morton_Reordering_Frequency") == 0) {
      i_command = 58;
      Line_Number = Line_Number + 1;
      Input_File >> Morton_Reordering_Frequency;
      Input_File.getline(buffer, sizeof(buffer));
      if (Morton_Reordering_Frequency < ZERO) i_command = INVALID_INPUT_VALUE;

    }else if (strcmp(Next_Control_Parameter, "Flow_Type") == 0) {
       i_command = 59;
       Get_Next_Input_Control_Parameter( );
       strcpy(Flow_Type, 
              Next_Control_Parameter);
       if (strcmp(Flow_Geometry_Type, "Inviscid") == 0) {
          i_Flow_Type = FLOWTYPE_INVISCID;
       } else if (strcmp(Flow_Type, "Laminar") == 0) {
	 i_Flow_Type = FLOWTYPE_LAMINAR;
       } else if (strcmp(Flow_Type, "Turbulent-k-epsilon") == 0) {
	 i_Flow_Type = FLOWTYPE_TURBULENT_RANS_K_EPSILON;
       } else if (strcmp(Flow_Type, "Turbulent-k-omega") == 0) {
	 i_Flow_Type = FLOWTYPE_TURBULENT_RANS_K_OMEGA;
       } else if (strcmp(Flow_Type, "Turbulent-LES") == 0) {
	 i_Flow_Type = FLOWTYPE_TURBULENT_LES;
       } else if (strcmp(Flow_Type, "Turbulent-DES-k-omega") == 0) {
	 i_Flow_Type = FLOWTYPE_TURBULENT_DES_K_OMEGA;
       } else if (strcmp(Flow_Type, "Turbulent-DNS") == 0) {
	 i_Flow_Type = FLOWTYPE_TURBULENT_DNS;
       } else {
         i_Flow_Type = FLOWTYPE_INVISCID;
       } /* endif */

    } else if (strcmp(Next_Control_Parameter, "Flow_Geometry_Type") == 0) {
       i_command = 60;
       Get_Next_Input_Control_Parameter( );
       strcpy(Flow_Geometry_Type, 
              Next_Control_Parameter);
       if (strcmp(Flow_Geometry_Type, "Planar") == 0) {
          Axisymmetric = 0;
       } else if (strcmp(Flow_Geometry_Type, "Axisymmetric") == 0) {
	  Axisymmetric = 1;
       } else if (strcmp(Flow_Geometry_Type, "Axisymmetric-x") == 0) {
	  Axisymmetric = 2;
       } else if (strcmp(Flow_Geometry_Type, "Axisymmetric-y") == 0) {
	  Axisymmetric = 1;
       } else {
          Axisymmetric = 0;
       } /* endif */

    } else if (strcmp(Next_Control_Parameter, "Restart_Solution_Save_Frequency") == 0) {
       i_command = 61;
       Line_Number = Line_Number + 1;
       Input_File >> Restart_Solution_Save_Frequency;
       Input_File.getline(buffer, sizeof(buffer));
       if (Restart_Solution_Save_Frequency < 1) i_command = INVALID_INPUT_VALUE;

   } else if (strcmp(Next_Control_Parameter, "Residual_Smoothing_Epsilon") == 0) {
      i_command = 70;
      Line_Number = Line_Number + 1;
      Input_File >> Residual_Smoothing_Epsilon;
      Input_File.getline(buffer, sizeof(buffer));
      if (Residual_Smoothing_Epsilon <= ZERO) {
         Residual_Smoothing = 0;
         Residual_Smoothing_Epsilon = ZERO;
      } else {
         Residual_Smoothing = 1;
      } /* endif */

    } else if (strcmp(Next_Control_Parameter, "Residual_Smoothing_Gauss_Seidel_Iterations") == 0) {
      i_command = 71;
      Line_Number = Line_Number + 1;
      Input_File >> Residual_Smoothing_Gauss_Seidel_Iterations;
      Input_File.getline(buffer, sizeof(buffer));
      if (Residual_Smoothing_Gauss_Seidel_Iterations < 0) {
         Residual_Smoothing_Gauss_Seidel_Iterations = 0;
      } /* endif */

    } else if (strcmp(Next_Control_Parameter, "AMR") == 0) {
      i_command = 72;
      Line_Number = Line_Number + 1;
      Input_File >> AMR;
      Input_File.getline(buffer, sizeof(buffer));
      if (AMR < 0) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(Next_Control_Parameter, "AMR_Frequency") == 0) {
      i_command = 73;
      Line_Number = Line_Number + 1;
      Input_File >> AMR_Frequency;
      Input_File.getline(buffer, sizeof(buffer));
      if (AMR_Frequency < 0) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(Next_Control_Parameter, "Number_of_Initial_Mesh_Refinements") == 0) {
      i_command = 74;
      Line_Number = Line_Number + 1;
      Input_File >> Number_of_Initial_Mesh_Refinements;
      Input_File.getline(buffer, sizeof(buffer));
      if (Number_of_Initial_Mesh_Refinements < 0) Number_of_Initial_Mesh_Refinements = 0;

    } else if (strcmp(Next_Control_Parameter,"Number_of_Uniform_Mesh_Refinements") == 0) {
      i_command = 75;
      Line_Number = Line_Number + 1;
      Input_File >> Number_of_Uniform_Mesh_Refinements;
      Input_File.getline(buffer,sizeof(buffer));
      if (Number_of_Uniform_Mesh_Refinements < 0) Number_of_Uniform_Mesh_Refinements = 0;

    } else if (strcmp(Next_Control_Parameter,"Number_of_Boundary_Mesh_Refinements") == 0) {
      i_command = 76;
      Line_Number = Line_Number + 1;
      Input_File >> Number_of_Boundary_Mesh_Refinements;
      Input_File.getline(buffer,sizeof(buffer));
      if (Number_of_Boundary_Mesh_Refinements < 0) Number_of_Boundary_Mesh_Refinements = 0;

    } else if (strcmp(Next_Control_Parameter,"Maximum_Refinement_Level") == 0) {
      i_command = 77;
      Line_Number = Line_Number + 1;
      Input_File >> Maximum_Refinement_Level;
      Input_File.getline(buffer,sizeof(buffer));
      if (Maximum_Refinement_Level < 1) Maximum_Refinement_Level = 1;

    } else if (strcmp(Next_Control_Parameter,"Minimum_Refinement_Level") == 0) {
      i_command = 78;
      Line_Number = Line_Number + 1;
      Input_File >> Minimum_Refinement_Level;
      Input_File.getline(buffer,sizeof(buffer));
      if (Minimum_Refinement_Level < 1) Minimum_Refinement_Level = 1;

    } else if (strcmp(Next_Control_Parameter, "Threshold_for_Refinement") == 0) {
       i_command = 79;
       Line_Number = Line_Number + 1;
       Input_File >> Threshold_for_Refinement;
       Input_File.getline(buffer, sizeof(buffer));
       if (Threshold_for_Refinement <= ZERO ||
           Threshold_for_Refinement > ONE) Threshold_for_Refinement = 0.50;

    } else if (strcmp(Next_Control_Parameter, "Threshold_for_Coarsening") == 0) {
      i_command = 80;
      Line_Number = Line_Number + 1;
      Input_File >> Threshold_for_Coarsening;
      Input_File.getline(buffer, sizeof(buffer));
      if (Threshold_for_Coarsening < ZERO ||
	  Threshold_for_Coarsening >= ONE) Threshold_for_Coarsening = 0.10;
      
    } else if (strcmp(Next_Control_Parameter, "Residual_Smoothing_Epsilon") == 0) {
       i_command = 81;
       Line_Number = Line_Number + 1;
       Input_File >> Residual_Smoothing_Epsilon;
       Input_File.getline(buffer, sizeof(buffer));
       if (Residual_Smoothing_Epsilon <= ZERO) {
          Residual_Smoothing = 0;
          Residual_Smoothing_Epsilon = ZERO;
       } else {
          Residual_Smoothing = 1;
       } /* endif */
     
    } else if (strcmp(Next_Control_Parameter, "Gravity") == 0) {
      i_command = 511; 
      Gravity = 1;
     
    } else if (strcmp(Next_Control_Parameter, "Schmidt") == 0) {
      i_command = 512; 
      Line_Number = Line_Number + 1;
      Input_File >> Global_Schmidt;
      Input_File.getline(buffer, sizeof(buffer));
      if (Global_Schmidt < 0) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(Next_Control_Parameter, "Debug_Level") == 0) {
      i_command = 513; 
      Line_Number = Line_Number + 1;
      Input_File >> debug_level;
      Input_File.getline(buffer, sizeof(buffer));
      if (debug_level != 0 && debug_level != 1 ) i_command = INVALID_INPUT_VALUE; 

    } else if (strcmp(Next_Control_Parameter, "Moving_Wall_Velocity_x") == 0) {
       i_command = 614;
       Line_Number = Line_Number + 1;
       Input_File >> Moving_wall_velocity.x;
       Input_File.getline(buffer, sizeof(buffer));

    }  else if (strcmp(Next_Control_Parameter, "Moving_Wall_Velocity_y") == 0) {
       i_command = 615;
       Line_Number = Line_Number + 1;
       Input_File >> Moving_wall_velocity.y;
       Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(Next_Control_Parameter, "Moving_Wall_Velocity_z") == 0) {
       i_command = 616;
       Line_Number = Line_Number + 1;
       Input_File >> Moving_wall_velocity.z;
       Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(Next_Control_Parameter, "Time_Accurate_Plot_Frequency") == 0) {
       i_command = 515; 
       Line_Number = Line_Number + 1;
       Input_File >> Time_Accurate_Plot_Freq;
       Input_File.getline(buffer, sizeof(buffer));
       if (Time_Accurate_Plot_Freq < 0) i_command = INVALID_INPUT_VALUE;
     
    }else if (strcmp(Next_Control_Parameter, "Wall_Boundary_Treatments") == 0) {
       i_command = 516;
       Line_Number = Line_Number + 1;
       Input_File >> Wall_Boundary_Treatments;
       Input_File.getline(buffer, sizeof(buffer));
        
    }else if (strcmp(Next_Control_Parameter, "Stretch_Level") == 0) {
       i_command = 517;
       Line_Number = Line_Number + 1;
       Input_File >> Stretch_Level;
       Input_File.getline(buffer, sizeof(buffer));
       
    }else if (strcmp(Next_Control_Parameter, "Pressure_Gradient_xdir") == 0) {
       i_command = 818;
       Line_Number = Line_Number + 1;
       Input_File >> Pressure_Gradient.x;
       Input_File.getline(buffer, sizeof(buffer));
    }else if (strcmp(Next_Control_Parameter, "Pressure_Gradient_ydir") == 0) {
       i_command = 819;
       Line_Number = Line_Number + 1;
       Input_File >> Pressure_Gradient.y;
       Input_File.getline(buffer, sizeof(buffer));
    }else if (strcmp(Next_Control_Parameter, "Pressure_Gradient_zdir") == 0) {
       i_command = 820;
       Line_Number = Line_Number + 1;
       Input_File >> Pressure_Gradient.z;
       Input_File.getline(buffer, sizeof(buffer));
    }else if (strcmp(Next_Control_Parameter, "Reynolds_Number") == 0) {
       i_command = 519;
       Line_Number = Line_Number + 1;
       Input_File >> Reynolds_Number;
       Input_File.getline(buffer, sizeof(buffer));
    }else if (strcmp(Next_Control_Parameter, "BluffBody_Data_Usage") == 0) {
       i_command = 520;
       BluffBody_Data_Usage = 1;
       
    }else if (strcmp(Next_Control_Parameter, "X_Scale") == 0) {
       i_command = 526;
       Line_Number = Line_Number + 1;
       Input_File >> X_Scale;
       Input_File.getline(buffer, sizeof(buffer));
       
    } else if (strcmp(Next_Control_Parameter, "X_Rotate") == 0) {
       i_command = 527;
       Line_Number = Line_Number + 1;
       Input_File >> X_Rotate;
       Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(Next_Control_Parameter, "Freeze_Limiter") == 0) {
      // Freeze_Limiter:
      i_command = 608;
      Line_Number = Line_Number + 1;
      Input_File >> Freeze_Limiter;
      Input_File.getline(buffer, sizeof(buffer));
      if (Freeze_Limiter < ZERO) i_command = INVALID_INPUT_VALUE;
      
    } else if (strcmp(Next_Control_Parameter, "Freeze_Limiter_Residual_Level") == 0) {
      // Freeze_Limiter_Residual_Level:
      i_command = 609;
      Line_Number = Line_Number + 1;
      Input_File >> Freeze_Limiter_Residual_Level;
      Input_File.getline(buffer, sizeof(buffer));
      if (Freeze_Limiter_Residual_Level < 0) i_command = INVALID_INPUT_VALUE;
             
    } else if (strcmp(Next_Control_Parameter, "Execute") == 0) {
       i_command = EXECUTE_CODE;
            
    } else if (strcmp(Next_Control_Parameter, "Terminate") == 0) {
       i_command = TERMINATE_CODE;
       
    } else if (strcmp(Next_Control_Parameter, "Continue") == 0) {
       i_command = CONTINUE_CODE;
       
    } else if (strcmp(Next_Control_Parameter, "Write_Output") == 0) {
       i_command = WRITE_OUTPUT_CODE;
       
    } else if (strcmp(Next_Control_Parameter, "Write_Output_Cells") == 0) {
       i_command = WRITE_OUTPUT_CELLS_CODE;
       
    } else if (strcmp(Next_Control_Parameter, "Write_Output_Nodes") == 0) {
       i_command = WRITE_OUTPUT_NODES_CODE;

    } else if (strcmp(Next_Control_Parameter, "Write_Restart") == 0) {
       i_command = WRITE_RESTART_CODE;
       
    } else if (strcmp(Next_Control_Parameter, "Write_Output_RHS") == 0) {
       i_command = WRITE_OUTPUT_RHS_CODE;
       
    } else if (strcmp(Next_Control_Parameter, "Write_Output_Mesh") == 0) {
       i_command = WRITE_OUTPUT_GRID_CODE;

    } else if (strcmp(Next_Control_Parameter, "Perturbation") == 0) {
       i_command = WRITE_OUTPUT_PERTURB_CODE;

    } else if (strcmp(Next_Control_Parameter, "Write_Mesh_Definition") == 0) {
       i_command = WRITE_GRID_DEFINITION_CODE;
       
    } else if (strcmp(Next_Control_Parameter, "Write_Output_Mesh_Nodes") == 0) {
       i_command = WRITE_OUTPUT_GRID_NODES_CODE;

    } else if (strcmp(Next_Control_Parameter, "Write_Output_Mesh_Cells") == 0) {
       i_command = WRITE_OUTPUT_GRID_CELLS_CODE;

    }  else if (strcmp(Next_Control_Parameter, "Morton_Ordering") == 0) {
       i_command = MORTON_ORDERING_CODE;

    } else if (Next_Control_Parameter[0] == '#') {
       i_command = COMMENT_CODE;

    } else if (strcmp(Next_Control_Parameter, "Refine_Grid") == 0) {
       i_command = REFINE_GRID_CODE;
       
    } else if (strcmp(Next_Control_Parameter, "Reaction_Mechanism") == 0) {
       
       i_command = 89;
       Get_Next_Input_Control_Parameter( );
       Deallocate();  //DEALLOCATE BEFORE CHANGING num_species
       int flag =0;
      
       //convert IP to string & define setup which Reaction Mechanism
       react_name = Next_Control_Parameter;
       Wo.React.set_reactions(react_name);
       
       num_species = Wo.React.num_species;      
       Allocate();
       
       //Get species and load appropriate data
       for(int i=0; i<num_species; i++){
          multispecies[i] = Wo.React.species[i];
          Schmidt[i] = Global_Schmidt;
       } 
       //Get next line and read in Schmidt numbers else will use defaults
       Get_Next_Input_Control_Parameter( );
       if (strcmp(Next_Control_Parameter, "Schmidt_Numbers") == 0){
          for(int i=0; i<num_species; i++){
             Input_File >> Schmidt[i];	 	 
          }	 
          //fudge the line number and istream counters
          Input_File.getline(buffer, sizeof(buffer)); 
          Line_Number = Line_Number + 1 ;
          flag = 1;
       } else { //Set to one value ~1
          for(int i=0; i<num_species; i++){
             Schmidt[i] = Global_Schmidt;	 	 
          }
          //To fix up line numbers
          Line_Number = Line_Number - 1 ;
       }
       
       //Set appropriate species data
       SOLN_pSTATE::set_species_data(Wo,
                                     num_species,
                                     multispecies,
                                     CFFC_Path,
                                     debug_level,
                                     Mach_Number_Reference,
                                     Schmidt); 
       SOLN_cSTATE::set_species_data(Uo,
                                     num_species,
                                     multispecies,
                                     CFFC_Path,
                                     debug_level,
                                     Mach_Number_Reference,
                                     Schmidt);
           
       //Get next line and read in mass fractions or set defaults
       if(flag){
          Get_Next_Input_Control_Parameter( );
       }
       if (strcmp(Next_Control_Parameter, "Mass_Fractions") == 0){
         
	 //Get Initial Mass Fractions from user	
          double temp=0.0;
          for(int i=0; i<num_species; i++){
                      
             Input_File >> mass_fractions[i];
             temp += mass_fractions[i];
                      
          }
          
          //check to make sure it adds to 1
          if(temp < ONE-MICRO || temp > ONE+MICRO){ 
             cout<<"\n Mass Fractions summed to "<<temp<<". Should sum to 1\n";
             i_command = INVALID_INPUT_VALUE;
          }
	 
          //Set inital Values; 
          Wo.set_initial_values(mass_fractions);  
          Uo.set_initial_values(mass_fractions);  
	  
          Uo = Wo.U();
          
          //fudge the line number and istream counters
          Input_File.getline(buffer, sizeof(buffer));  
          Line_Number = Line_Number + 1; 
          
          //Spit out appropriate mass fractions and exit
       } 
       
    } else if (strcmp(Next_Control_Parameter, "User_Reaction_Mechanism") == 0) { 
       // this will be added but its not quite yet
       i_command=90;
       cout<<endl<<Next_Control_Parameter<<"\n not currently available in freeware version :)\n";
       i_command = INVALID_INPUT_VALUE;   
       
    } else if (strcmp(Next_Control_Parameter, "Species") == 0) { 
          
       i_command = 91;
       
       Deallocate();  //DEALLOCATE BEFORE CHANGING num_species
       
       // Non Reaction case so set NO_REACTIONS flag in reactions class
       react_name ="NO_REACTIONS";
       Wo.React.set_reactions(react_name);
       
       //read in the number of species (should be first in line) 
       Input_File>>num_species;   
             
       Allocate();
       
       //read in species names
       for(int i=0; i<num_species; i++){
          Input_File >> multispecies[i];
          Schmidt[i] = Global_Schmidt;
         
       }
       
       //copy names into Reaction class for storage
       Wo.React.set_species(multispecies,num_species);
       
       //Setup State class data and find species thermo and transport properties
       SOLN_pSTATE::set_species_data(Wo,
                                     num_species,
                                     multispecies,
                                     CFFC_Path,
                                     debug_level,
                                     Mach_Number_Reference,
                                     Schmidt); 
       SOLN_cSTATE::set_species_data(Uo,
                                     num_species,
                                     multispecies,
                                     CFFC_Path,
                                     debug_level,
                                     Mach_Number_Reference,
                                     Schmidt);
         
       //More Fudging of lines 
       Line_Number = Line_Number + 1 ;
       Input_File.getline(buffer, sizeof(buffer));  
       
       //Get next line and read in mass fractions or set defaults
       Get_Next_Input_Control_Parameter();
       
       if (strcmp(Next_Control_Parameter, "Mass_Fractions") == 0){
          //Get Initial Mass Fractions from user 
          double temp=0.0;
          for(int i=0; i<num_species; i++){
             Input_File >> mass_fractions[i];
             temp += mass_fractions[i];
          }
          //check to make sure it adds to 1
          if(temp < ONE-MICRO || temp > ONE+MICRO){ 
             cout<<"\n Mass Fractions summed to "<<temp<<". Should be sum to 1\n";
             i_command = INVALID_INPUT_VALUE;
          }
          //Set inital Values; 
          Wo.set_initial_values(mass_fractions);  
          Uo.set_initial_values(mass_fractions);
          Uo = Wo.U();
          
          //fudge the line number and istream counters
          Input_File.getline(buffer, sizeof(buffer));  
          Line_Number = Line_Number + 1; 
       } 
       //If no mass fraction data is set to defaults (all equal to 1/num_species)
       else{        
          Uo = Wo.U();
          Line_Number = Line_Number - 1 ;
       }
       
    } else if (strcmp(Next_Control_Parameter, "Temperature") == 0) {
       
       i_command = 92;
       Line_Number = Line_Number + 1;
       Input_File >> Temperature;
       Input_File.getline(buffer, sizeof(buffer));
       if (Temperature <= ZERO) {
          i_command = INVALID_INPUT_VALUE;
       } else {
          
          Wo.rho = Pressure/(Wo.Rtot()*Temperature); 	
          
       } /* endif */
            
    } else if (strcmp(Next_Control_Parameter, "Pressure") == 0) {
       i_command = 93;
       Line_Number = Line_Number + 1;
       Input_File >> Pressure;
       Input_File.getline(buffer, sizeof(buffer));
       Pressure = Pressure*THOUSAND;
       if (Pressure <= ZERO) {
          i_command = INVALID_INPUT_VALUE;
       } else {
          Wo.rho = Pressure/(Wo.Rtot()*Temperature); 
          Wo.p = Pressure;	
          //Wo.v.zero();
       } /* endif */
       
    } else if (strcmp(Next_Control_Parameter, "Schmidt") == 0) {
       i_command = 94; 
       Line_Number = Line_Number + 1;
       Input_File >> Global_Schmidt;
       Input_File.getline(buffer, sizeof(buffer));
       if (Global_Schmidt < 0) i_command = INVALID_INPUT_VALUE;
       
    } else {
       i_command = INVALID_INPUT_CODE;
       cout<<"\n ERROR ERROR WILL ROBINSON.... INVALID_INPUT_CODE \n";
    } /* endif */

      
     /* Return the parser command type indicator. */
    return (i_command);
    
}

/********************************************************
 * Routine: Process_Input_Control_Parameter_File        *
 *                                                      *
 * Reads, parses, and executes the list of input        *
 * control parameters from the standard input file.     *
 *                                                      *
 ********************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
   int Input_Parameters<SOLN_pSTATE, SOLN_cSTATE>::
   Process_Input_Control_Parameter_File(
      char *Input_File_Name_ptr,
      int &Command_Flag) {
   
   int error_flag, line_number;
   
   /* Assign initial value for error indicator flag. */
   error_flag = 0;
   
   /* Assign default values to the input parameters. */
   Set_Default_Input_Parameters();
   
   /* Copy input file name (a string) to appropriate input parameter variable. */
   if (Input_File_Name_ptr != NULL) strcpy(Input_File_Name, Input_File_Name_ptr);
   
   /* Open the input file containing the input parameters. */
   Open_Input_File();
   error_flag = Input_File.bad();
   
   if (error_flag) {
      cout << "\n ERROR: Unable to open input data file.\n";
      return (error_flag);
   } /* endif */
    
   
   /* Read and parse control parameters contained in
      the input file. */
   while (1) {
      
      Get_Next_Input_Control_Parameter();
     
       Command_Flag = Parse_Next_Input_Control_Parameter();

       line_number = Line_Number;
       
       if (Command_Flag == EXECUTE_CODE) {
          
          break;
       } else if (Command_Flag == TERMINATE_CODE) {
          
          break;
       } else if (Command_Flag == INVALID_INPUT_CODE ||
                  Command_Flag == INVALID_INPUT_VALUE) {
          line_number = -line_number;
          cout << "\n ERROR: Error reading data at line #"
               << -line_number  << " of input data file.\n";
          error_flag = line_number;
          break;
          
       } /* endif */
     
       
    } /* endwhile */

    /* Initial processing of input control parameters complete.  
       Return the error indicator flag. */
    
    //Load the C-Type strings from the C++ strings 
    strcpy(React_Name,react_name.c_str());
    
    for (int i = 0; i < num_species; i++) {
       strcpy(Multispecies[i],multispecies[i].c_str());
    } 
 
    //Proper temperature 
    Temperature = Wo.T();
   

    return (error_flag);

}

/********************************************************
 * Routine: Broadcast_Input_Parameters                  *
 *                                                      *
 * Broadcast the input parameters variables to all      *
 * processors involved in the calculation from the      *
 * primary processor using the MPI broadcast routine.   *
 *                                                      *
*********************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
                  void Input_Parameters<SOLN_pSTATE, 
                                        SOLN_cSTATE>::
                  Broadcast_Input_Parameters(void) {
   
#ifdef _MPI_VERSION

   // Input file parameters  
   MPI::COMM_WORLD.Bcast(CFFC_Path, 
			 INPUT_PARAMETER_LENGTH, 
			 MPI::CHAR, 0);
   MPI::COMM_WORLD.Bcast(Input_File_Name, 
                         INPUT_PARAMETER_LENGTH, 
                         MPI::CHAR, 0);   
   MPI::COMM_WORLD.Bcast(&(Line_Number), 
                         1, 
                         MPI::INT, 0);

   // Grid parameters
   MPI::COMM_WORLD.Bcast(IP_Grid.Grid_Type, 
                         INPUT_PARAMETER_LENGTH, 
                         MPI::CHAR, 0);
   MPI::COMM_WORLD.Bcast(&(IP_Grid.i_Grid), 
                         1, 
                         MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(&(IP_Grid.NBlk_Idir), 
                         1, 
                         MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(&(IP_Grid.NBlk_Jdir), 
                         1, 
                         MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(&(IP_Grid.NBlk_Kdir), 
                         1, 
                         MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(&(IP_Grid.NCells_Idir), 
                         1, 
                         MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(&(IP_Grid.NCells_Jdir), 
                         1, 
                         MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(&(IP_Grid.NCells_Kdir), 
                         1, 
                         MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(&(IP_Grid.Nghost), 
                         1, 
                         MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(&(IP_Grid.Box_Length), 
                         1, 
                         MPI::DOUBLE, 0);
   MPI::COMM_WORLD.Bcast(&(IP_Grid.Box_Width), 
                         1, 
                         MPI::DOUBLE, 0);
   MPI::COMM_WORLD.Bcast(&(IP_Grid.Box_Height), 
                         1, 
                         MPI::DOUBLE, 0);
   MPI::COMM_WORLD.Bcast(&(IP_Grid.Stretching_Type_Idir),
	         	 1,
			 MPI::INT,0);
   MPI::COMM_WORLD.Bcast(&(IP_Grid.Stretching_Type_Jdir),
			 1,
			 MPI::INT,0);
   MPI::COMM_WORLD.Bcast(&(IP_Grid.Stretching_Type_Kdir),
			 1,
			 MPI::INT,0);
   MPI::COMM_WORLD.Bcast(&(IP_Grid.Stretching_Factor_Idir),
			 1,
			 MPI::DOUBLE,0);
   MPI::COMM_WORLD.Bcast(&(IP_Grid.Stretching_Factor_Jdir),
			 1,
			 MPI::DOUBLE,0);
   MPI::COMM_WORLD.Bcast(&(IP_Grid.Stretching_Factor_Kdir),
			 1,
			 MPI::DOUBLE,0);

   // Solver parameters
   MPI::COMM_WORLD.Bcast(Time_Integration_Type, 
                         INPUT_PARAMETER_LENGTH, 
                         MPI::CHAR, 0);
   MPI::COMM_WORLD.Bcast(&(i_Time_Integration), 
                         1, 
                         MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(&(Time_Accurate), 
                         1, 
                         MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(&(Local_Time_Stepping), 
                         1, 
                         MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(&(Maximum_Number_of_Time_Steps), 
                         1, 
                         MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(&(N_Stage), 
                         1, 
                         MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(&(CFL_Number), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Time_Max), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Residual_Smoothing), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(p_Norms_Specified_Parameter), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Residual_Smoothing_Epsilon), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Residual_Smoothing_Gauss_Seidel_Iterations), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(Reconstruction_Type, 
                          INPUT_PARAMETER_LENGTH, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(i_Reconstruction), 
                          1, 
                          MPI::INT, 0);   
     MPI::COMM_WORLD.Bcast(Limiter_Type, 
                          INPUT_PARAMETER_LENGTH, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(i_Limiter), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(Flux_Function_Type, 
                          INPUT_PARAMETER_LENGTH, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(i_Flux_Function), 
                          1, 
                          MPI::INT, 0);

    // Flow solution parameters
    MPI::COMM_WORLD.Bcast(ICs_Type, 
                          INPUT_PARAMETER_LENGTH, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(i_ICs), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(i_Grid_Level), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Pressure), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Temperature), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Mach_Number), 
			  1, 
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Mach_Number_Reference), 
			  1, 
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Flow_Angle), 
                          1, 
                          MPI::DOUBLE, 0);   
    MPI::COMM_WORLD.Bcast(&(Re_lid),
			  1,
			  MPI::DOUBLE,0);
     MPI::COMM_WORLD.Bcast(Flow_Type, 
                          INPUT_PARAMETER_LENGTH, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(i_Flow_Type), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(Flow_Geometry_Type, 
                          INPUT_PARAMETER_LENGTH, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(Axisymmetric), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Wall_Boundary_Treatments), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Stretch_Level), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Global_Schmidt), 
                          1, 
                          MPI::INT, 0); 
    MPI::COMM_WORLD.Bcast(&(Gravity), 
                          1, 
			  MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(debug_level), 
                          1, 
			  MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Preconditioning), 
                          1, 
			  MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(BluffBody_Data_Usage), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Moving_wall_velocity.x), 
                          1, 
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Moving_wall_velocity.y), 
                          1, 
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Moving_wall_velocity.z), 
                          1, 
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Pressure_Gradient.x), 
                          1, 
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Pressure_Gradient.y), 
                          1, 
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Pressure_Gradient.z), 
                          1, 
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Reynolds_Number), 
                          1, 
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Kinematic_Viscosity_Wall), 
                          1, 
			  MPI::DOUBLE, 0);

    // More grid parameters
    MPI::COMM_WORLD.Bcast(&(Plate_Length), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Pipe_Length), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Pipe_Radius), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Blunt_Body_Radius), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Blunt_Body_Mach_Number), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Grain_Length), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Grain_Radius), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Grain_To_Throat_Length), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Nozzle_Length), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Nozzle_Radius_Exit), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Nozzle_Radius_Throat), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Length_Shroud), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Radius_Shroud), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Length_BluffBody), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Radius_BluffBody), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Radius_Orifice), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(BluffBody_Coflow_Air_Velocity), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(BluffBody_Coflow_Fuel_Velocity), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Radius_Inlet_Pipe), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Radius_Combustor_Tube), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Length_Inlet_Pipe), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Length_Combustor_Tube), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Cylinder_Radius), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Ellipse_Length_X_Axis), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Ellipse_Length_Y_Axis), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Chord_Length), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Orifice_Radius), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Wedge_Angle), 
			  1, 
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Wedge_Length), 
			  1, 
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(X_Scale), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(X_Rotate), 
                          1, 
                          MPI::DOUBLE, 0);
   
    // AMR & Refinement Parameters
    MPI::COMM_WORLD.Bcast(&(AMR), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(AMR_Frequency),
                          1,
                          MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(Number_of_Initial_Mesh_Refinements), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Number_of_Uniform_Mesh_Refinements),
                          1,
                          MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(Number_of_Boundary_Mesh_Refinements),
                          1,
                          MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(Maximum_Refinement_Level),
                          1,
                          MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(Minimum_Refinement_Level),
                          1,
                          MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(Threshold_for_Refinement), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Threshold_for_Coarsening), 
                          1, 
                          MPI::DOUBLE, 0);

    // Morton Ordering Parameters
    MPI::COMM_WORLD.Bcast(&(IP.Morton), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Morton_Reordering_Frequency),
                          1,
                          MPI::INT,0);

    // File Names
    MPI::COMM_WORLD.Bcast(Output_File_Name, 
                          INPUT_PARAMETER_LENGTH, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(Grid_File_Name, 
                          INPUT_PARAMETER_LENGTH, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(Grid_Definition_File_Name, 
                          INPUT_PARAMETER_LENGTH, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(Restart_File_Name, 
                          INPUT_PARAMETER_LENGTH, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(Gnuplot_File_Name, 
                          INPUT_PARAMETER_LENGTH, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(Output_Format_Type, 
                          INPUT_PARAMETER_LENGTH, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(i_Output_Format), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Restart_Solution_Save_Frequency), 
                          1, 
                          MPI::INT, 0); 

    if (!CFFC_Primary_MPI_Processor()) {
       Number_of_Processors = CFFC_MPI::Number_of_Processors;
    } /* endif */
    MPI::COMM_WORLD.Bcast(&(Number_of_Blocks_Per_Processor), 
                          1, 
                          MPI::INT, 0); 
    // Freeze_Limiter
    MPI::COMM_WORLD.Bcast(&(Freeze_Limiter), 
			  1, 
			  MPI::INT, 0);
    // Freeze_Limiter_Residual_Level
    MPI::COMM_WORLD.Bcast(&(Freeze_Limiter_Residual_Level),
                          1, 
                          MPI::DOUBLE, 0);

   //reaction name
   MPI::COMM_WORLD.Bcast(React_Name, 
                         INPUT_PARAMETER_LENGTH, 
                         MPI::CHAR, 0);
    
   //delete current dynamic memory before changing num_species
   if(!CFFC_Primary_MPI_Processor()) {   
      Deallocate();
   } 
   // Number of species
   MPI::COMM_WORLD.Bcast(&(num_species), 
                         1, 
                         MPI::INT, 0);
   // Set up new dynamic memory
   if(!CFFC_Primary_MPI_Processor()) {   
      Allocate();
   } 

   //species names & mass fractions
   for(int i =0; i < num_species; i++){
      MPI::COMM_WORLD.Bcast(&(mass_fractions[i]), 
			    1, 
			    MPI::DOUBLE, 0);
      MPI::COMM_WORLD.Bcast(&(Schmidt[i]), 
			    1, 
			    MPI::DOUBLE, 0);
      MPI::COMM_WORLD.Bcast(Multispecies[i], 
			    INPUT_PARAMETER_LENGTH, 
			    MPI::CHAR, 0);
    }
   //set recaction and species parameters
   if (!CFFC_Primary_MPI_Processor()) {      
      react_name = React_Name;
      for (int i = 0; i < num_species; i++) {
         multispecies[i] = Multispecies[i];  
      }    
      
      //load reaction names
      Wo.React.set_reactions(react_name);
      
      //Set species if non-reacting
      if( Wo.React.reactset_flag == NO_REACTIONS){
         Wo.React.set_species(multispecies,num_species);
      }  
      
      SOLN_pSTATE::set_species_data(Wo,
                                    num_species,
                                    multispecies,
                                    CFFC_Path,
                                    debug_level,
                                    Mach_Number_Reference,
                                    Schmidt); 
      SOLN_cSTATE::set_species_data(Uo,
                                    num_species,
                                    multispecies,
                                    CFFC_Path,
                                    debug_level,
                                    Mach_Number_Reference,
                                    Schmidt);
    
      Wo.set_initial_values(mass_fractions);
      Uo.set_initial_values(mass_fractions);
      
      //set proper Temp & Pressure instead of defaults
      Wo.rho = Pressure/(Wo.Rtot()*Temperature); 
      Wo.p = Pressure;	
      Wo.v.zero();

      Uo = Wo.U();
   } 
   
   if(!CFFC_Primary_MPI_Processor()) {   
      Wo.v.x = Mach_Number*Wo.a()*cos(TWO*PI*Flow_Angle/360.00);
      Wo.v.y = Mach_Number*Wo.a()*sin(TWO*PI*Flow_Angle/360.00);
      Wo.v.z = Mach_Number*Wo.a()*sin(TWO*PI*Flow_Angle/360.00);
    }

#endif

}

/*************************************************************
 * Input_Parameters -- Input-output operators.               *
 *************************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
   ostream &operator << (ostream &out_file,
                         const Input_Parameters<SOLN_pSTATE, SOLN_cSTATE> &IP) {
   
   out_file << setprecision(6);
   
   out_file << "\n  -> CFFC Path: " << IP.CFFC_Path;
   out_file << "\n\n Solving 3D MulitSpecies";
   if (IP.i_Flow_Type ==  FLOWTYPE_INVISCID){
      out_file<<" Euler (Inviscid) ";
   } else {
      out_file<<" Navier-Stokes (Viscous) ";
   }
   out_file<<"equations (IBVP/BVP) "; 
   if (IP.i_Flow_Type ==  FLOWTYPE_INVISCID) {
   } else if (IP.i_Flow_Type == FLOWTYPE_LAMINAR) {
      out_file << "\n  -> Laminar flow";
   } else if (IP.i_Flow_Type == FLOWTYPE_TURBULENT_RANS_K_EPSILON) {
      out_file << "\n  -> Turbulent flow: RANS with k-epsilon turbulence model";
      if (IP.Wall_Boundary_Treatments==0){
         out_file <<" (Automatic wall boundary treatments)";
      } else {
         if(IP.Wall_Boundary_Treatments==2){
            out_file <<" (Low-Reynolds-number formulation)";
         }else{
            out_file <<" (Wall-function formulation)";
         }
      }
   } else if (IP.i_Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
      out_file << "\n  -> Turbulent flow: RANS with k-oemga turbulence model";
      if (IP.Wall_Boundary_Treatments==0){
         out_file <<" (Automatic wall boundary treatments)";
      } else {
         if(IP.Wall_Boundary_Treatments==2){
            out_file <<" (Low-Reynolds-number formulation)";
         }else{
            out_file <<" (Wall-function formulation)";
         }
      }
      
   } else if (IP.i_Flow_Type == FLOWTYPE_TURBULENT_LES) {
      out_file << "\n  -> Turbulent flow: LES ";
   } else if (IP.i_Flow_Type == FLOWTYPE_TURBULENT_DES_K_OMEGA) {
      out_file << "\n  -> Turbulent flow: DES with k-oemga SGS turbulence model ";
   } else if (IP.i_Flow_Type == FLOWTYPE_TURBULENT_DNS) {
      out_file << "\n  -> Turbulent flow: DNS ";
   }

   out_file << "\n  -> Input File Name: " 
            << IP.Input_File_Name;
   if (IP.Time_Accurate) { 
      out_file << "\n  -> Time Accurate (Unsteady) Solution";
   } else {
      out_file << "\n  -> Time Invariant (Steady-State) Solution";
   }
   
   if(IP.Gravity) {
      out_file << "\n  -> With Gravity (-z)";
   }

   if(IP.debug_level){
      out_file << "\n  -> Debug level: "
	       << IP.debug_level;
   }
   out_file << "\n  -> Time Integration: " 
            << IP.Time_Integration_Type;
   out_file << "\n  -> Number of Stages in Multi-Stage Scheme: " 
            << IP.N_Stage;
   if (IP.Local_Time_Stepping == GLOBAL_TIME_STEPPING) {
      out_file << "\n  -> Global Time Stepping";
   } else if (IP.Local_Time_Stepping == SCALAR_LOCAL_TIME_STEPPING) {
      out_file << "\n  -> Scalar Local Time Stepping";
   } else if (IP.Local_Time_Stepping == MATRIX_LOCAL_TIME_STEPPING) {
      out_file << "\n  -> Matrix Local Time Stepping";
   } else if (IP.Local_Time_Stepping == LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER) {
      out_file << "\n  -> Low-Mach-Number Local Preconditioning (Weiss-Smith)";
   } else if (IP.Local_Time_Stepping == SEMI_IMPLICIT_LOCAL_TIME_STEPPING) {
      out_file << "\n  -> Semi-Implicit Local Time Stepping";
   } else if (IP.Local_Time_Stepping == SEMI_IMPLICIT_LOW_MACH_NUMBER_PRECONDITIONER) {
      out_file << "\n  -> Semi-Implicit Low-Mach-Number Local Preconditioned Time Stepping";
   }
   if (IP.Preconditioning == 1){
      out_file <<"\n  -> Mach Number Reference "<<IP.Mach_Number_Reference;
   }
   
   if (IP.Residual_Smoothing) {
      out_file << "\n  -> Residual Smoothing";
      out_file << "\n  -> Epsilon: " 
               << IP.Residual_Smoothing_Epsilon;
      out_file << "\n  -> Gauss_Seidel_Iterations: " 
               << IP.Residual_Smoothing_Gauss_Seidel_Iterations;
   }

   out_file << "\n  -> Reconstruction: " 
            << IP.Reconstruction_Type;
   out_file << "\n  -> Limiter: " 
            << IP.Limiter_Type;   
   if (IP.Limiter_Type != LIMITER_ZERO && IP.Freeze_Limiter) {
      out_file << "\n  -> Freeze Limiter when L2-norm of residual is < "
	       << IP.Freeze_Limiter_Residual_Level;
   } 
   out_file << "\n  -> Flux Function: " 
            << IP.Flux_Function_Type;
   
   out_file << "\n  -> Reaction Mechanism: " 
            << IP.react_name;
   out_file << "\n  -> Species: "<<IP.Wo.ns
            << "\n  -> Initial mass fractions: ";
   for(int i=0; i<IP.Wo.ns; i++){
      out_file  <<"c["<<IP.multispecies[i]<<"]= ";
      out_file  << IP.Wo.spec[i].c<<", ";
   } 
//    if(IP.i_Flow_Type != FLOWTYPE_INVISCID){
//       out_file << "\n  -> Schmidt Numbers for Viscous flow: ";
//       for(int i=0; i <IP.Wo.ns; i++){
//          out_file  <<"Sc["<<IP.multispecies[i]<<"]= ";
//          out_file << IP.Schmidt[i]<<", ";
//       }
//    }
   out_file << "\n  -> Initial Conditions: " 
            << IP.ICs_Type;
   switch(IP.i_ICs) {
   case IC_CONSTANT :
      out_file << "\n  -> Pressure (kPa): " 
               << IP.Pressure/THOUSAND;
      out_file << "\n  -> Temperature (K): " 
               << IP.Temperature;
      out_file << "\n  -> Mach Number: " 
               << IP.Mach_Number;
      out_file << "\n  -> Flow Angle: " 
               << IP.Flow_Angle;
      break;
   case IC_UNIFORM :
      out_file << "\n  -> Pressure (kPa): " 
               << IP.Pressure/THOUSAND;
      out_file << "\n  -> Temperature (K): " 
               << IP.Temperature;
      out_file << "\n  -> Mach Number: " 
               << IP.Mach_Number;
      out_file << "\n  -> Flow Angle (degrees): " 
               << IP.Flow_Angle;
      break;
   case IC_RIEMANN_YDIR :
      out_file << "\n  -> Pressure (kPa): " 
               << IP.Pressure/THOUSAND;
      out_file << "\n  -> Temperature (K): " 
               << IP.Temperature;
      out_file << "\n  -> Mach Number: " 
               << IP.Mach_Number;
      out_file << "\n  -> Flow Angle (degrees): " 
               << IP.Flow_Angle;
      break;
   case IC_GAS_MIX :
      break;
   case IC_CHEM_CORE_FLAME:
      break;
   case IC_CHEM_INVERSE_FLAME:
      break;
   case IC_CHEM_1DFLAME:
      break;
   case IC_PRESSURE_GRADIENT_X:
      break;
   case IC_PRESSURE_GRADIENT_Y:
      break;
   case IC_VISCOUS_COUETTE:
      break;
   case IC_VISCOUS_COUETTE_PRESSURE_GRADIENT: 
      break;
      
   default:
      break;
   } /* endswitch */
   out_file << "\n  -> Grid: "
            << IP.IP_Grid.Grid_Type;
   switch(IP.IP_Grid.i_Grid) {

   case GRID_CUBE :
      out_file << "\n  -> Length of Solution Domain (m): "
               << IP.IP_Grid.Box_Length;
      out_file << "\n  -> Width of Solution Domain (m): "
               << IP.IP_Grid.Box_Width;
      out_file << "\n  -> Height of Solution Domain (m): "
               << IP.IP_Grid.Box_Height;
      break;
   case GRID_COUETTE :
      out_file << "\n  -> Length of Solution Domain (m): "
               << IP.IP_Grid.Box_Length;
      out_file << "\n  -> Width of Solution Domain (m): "
               << IP.IP_Grid.Box_Width;
      out_file << "\n  -> Height of Solution Domain (m): "
               << IP.IP_Grid.Box_Height;
      out_file << "\n  -> Moving Wall Velocity in x direction (m/s): "
               << IP.Moving_wall_velocity.x;
      out_file << "\n  -> Moving Wall Velocity in y direction (m/s): "
               << IP.Moving_wall_velocity.y;
      out_file << "\n  -> Moving Wall Velocity in z direction (m/s): "
               << IP.Moving_wall_velocity.z;
      break;

   default:
      out_file << "\n  -> Length of Solution Domain (m): "
               << IP.IP_Grid.Box_Length;
      out_file << "\n  -> Width of Solution Domain (m): "
               << IP.IP_Grid.Box_Width;
      out_file << "\n  -> Height of Solution Domain (m): "
               << IP.IP_Grid.Box_Height;
      break;
   } /* endswitch */
 
   if (IP.Number_of_Initial_Mesh_Refinements > 0)
      out_file << "\n  -> Number of Initial Mesh Refinements : " 
               << IP.Number_of_Initial_Mesh_Refinements;
   if (IP.Number_of_Uniform_Mesh_Refinements > 0)
      out_file << "\n  -> Number of Uniform Mesh Refinements : " 
               << IP.Number_of_Uniform_Mesh_Refinements;
   if (IP.Number_of_Boundary_Mesh_Refinements > 0)
      out_file << "\n  -> Number of Boundary Mesh Refinements : " 
               << IP.Number_of_Boundary_Mesh_Refinements;
   out_file << "\n  -> Number of Blocks i-direction: "
            << IP.IP_Grid.NBlk_Idir;
   out_file << "\n  -> Number of Blocks j-direction: " 
            << IP.IP_Grid.NBlk_Jdir;
   out_file << "\n  -> Number of Blocks k-direction: "
            << IP.IP_Grid.NBlk_Kdir;
   out_file << "\n  -> Number of Cells i-direction: "
            << IP.IP_Grid.NCells_Idir;
   out_file << "\n  -> Number of Cells j-direction: " 
            << IP.IP_Grid.NCells_Jdir;
   out_file << "\n  -> Number of Cells k-direction: " 
            << IP.IP_Grid.NCells_Kdir;
   out_file << "\n  -> CFL Number: " 
            << IP.CFL_Number;
   out_file << "\n  -> Maximum Time (ms): " 
            << IP.Time_Max*THOUSAND;
   out_file << "\n  -> Maximum Number of Time Steps (Iterations): " 
            << IP.Maximum_Number_of_Time_Steps;
   out_file << "\n  -> Number of Processors: " 
            << IP.Number_of_Processors;
   out_file << "\n  -> Number of Blocks Per Processor: " 
            << IP.Number_of_Blocks_Per_Processor;
   out_file << "\n  -> Output File Name: " 
            << IP.Output_File_Name;
   out_file << "\n  -> Output Format: " 
            << IP.Output_Format_Type;
   out_file << "\n  -> Restart Solution Save Frequency: "
            << IP.Restart_Solution_Save_Frequency
            << " steps (iterations)"; 
   if (IP.Time_Accurate_Plot_Freq !=0 && IP.Time_Accurate){
      out_file << "\n  -> Time Accurate Solution Plot Frequency: "
	       << IP.Time_Accurate_Plot_Freq
	       << " steps (iterations)"; 
   }
   return (out_file);
}

template<class SOLN_pSTATE, class SOLN_cSTATE>
istream  &operator >> (istream &in_file,
                       Input_Parameters<SOLN_pSTATE, SOLN_cSTATE> &IP) {
   return (in_file);
}

#endif /* _INPUT_INCLUDED  */

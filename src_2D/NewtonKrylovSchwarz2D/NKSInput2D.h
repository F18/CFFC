/* NKS2DInput.h:  Header file declaring related Inputs and functions */

#ifndef _NKS2DINPUT_INCLUDED
#define _NKS2DINPUT_INCLUDED
   
enum Block_Preconditioners { Block_ILUK, 
			     Block_Jacobi, 
			     Something_else };

enum Jacobian_Orders { SOURCE_TERMS_ONLY,
                       FIRST_ORDER_INVISCID_HLLE,
		       FIRST_ORDER_INVISCID_ROE,
		       FIRST_ORDER_INVISCID_AUSMPLUSUP,
		       SECOND_ORDER_DIAMOND_WITH_HLLE,
		       SECOND_ORDER_DIAMOND_WITH_ROE,
		       SECOND_ORDER_DIAMOND_WITH_AUSMPLUSUP,
                       SECOND_ORDER_OTHER };

enum Frechet_Derivatives { FIRST_ORDER = 1,
			   SECOND_ORDER = 2,
                           FIRST_ORDER_RESTART = 3,
                           SECOND_ORDER_RESTART = 4 };

enum output_formats { OF_SCOTT, OF_ALISTAIR };

enum Freeze_Limiter_Immediately_values { FLI_NOT_USED, FLI_NO, FLI_YES };

/*! ***************************************************************
 * Class: NKS_Input_Parameters
 *
 * @brief Input Parameters for NKS.
 *
 * This class defines and handles the input variables related to the
 * Newton Krylov Schwarz (NKS) method for solving 
 * nonlinear partial differential equations.  Note that unsteady
 * calculations can be performed using dual-time-stepping.
 *
 * Class: NKS_Input_Parameters                  
 *****************************************************************/ 
class NKS_Input_Parameters{
 private:
 public:
  //@{ @name Newton Parameters
  int    Maximum_Number_of_NKS_Iterations;  //!<Outer
  double Overall_Tolerance;
  double Relaxation_multiplier;
  //@}

  //@{ @name Dual Time Stepping Parameters
  bool Dual_Time_Stepping;              //!< Dual-time-stepping flag (on or off).  
  int Physical_Time_Integration;        //!< Implicit Euler, BDF2, ESDIRK, etc.
  double Physical_Time_CFL_Number;
  double Physical_Time_Step;
  int Maximum_Number_of_DTS_Steps; 
  int Time_Accurate_Plot_Frequency;
  double Total_Physical_Time;              //!< Storage HACK for BCs

  //int Dual_Time_Preconditioning;  
  //@}

  //@{ @name Implicit Euler Parameters
  bool   Finite_Time_Step;   
  double Finite_Time_Step_Initial_CFL;
  // Some of these finite-time step parameters are only used by explicit
  // specializations of the Finite_Time_Step() function.
  double Finite_Time_Step_Final_CFL; 
  double Finite_Time_Step_Max_CFL;
  //@}

  //@{ @name GMRES parameters 
  int    Maximum_Number_of_GMRES_Iterations; //!<Inner
  int    GMRES_Restart;
  int    GMRES_Overlap;
  bool   Normalization;   
  
  // The linear (GMRES) L2-norm tolerance at each Newton iteration 
  // is calculated as:
  //    (i) (L)^(x) 
  // where:
  //     i is GMRES_Initial_Tolerance
  //     L is the L2-norm of the non-linear (Newton) residual ratio
  // and x is calculated such that:
  //    (i) (Overall_Tolerance)^(x) == GMRES_Final_Tolerance
  // 
  // To turn off this variable linear tolerance, simply set 
  // the final tolerance equal to the initial tolerance.
  double GMRES_Initial_Tolerance;
  double GMRES_Final_Tolerance;
  //@}

  //@{ @name Matrix Free
  bool   GMRES_CHECK;
  int    GMRES_Frechet_Derivative_Order; 
  bool   GMRES_Frechet_Derivative_Chem;
  double Epsilon_Naught;
  //@}

  //@{ @name Preconditioner 
  int    GMRES_Block_Preconditioner; 
  int    Jacobian_Order;  
  int    GMRES_ILUK_Level_of_Fill;  
  //@}

  //@{ @name Convergence Parameters
  int    DCS_Window;  //!< size of the window for the detect convergence stall algorithm.
  bool   Detect_Convergence_Stall;
  // Hack by Alistair to play with Newton convergence.
  int    Freeze_Limiter_Immediately; 
  //@}

  //@{ @name Output & Debugging 
  int    output_format;
  int    output_precision, output_width;  
  int    NKS_Write_Output_Cells_Freq; //!< set to zero to turn off
  //@}

  //@{ @name NKS Parameters
  int Min_Number_of_Newton_Steps_With_Zero_Limiter;            //!< force 1st order for "N" steps
  int Number_of_Newton_Steps_Requiring_Jacobian_Update;        //!< force Jacobian updates for "N" Newton steps       
  double Min_L2_Norm_Requiring_Jacobian_Update;                //!< force Jacobian update for L2 < "N"
  double Min_Finite_Time_Step_Norm_Ratio;                      //!< ramp over to full newton over 8 orders of L2 magnitude
  int Jacobian_Update_Frequency;                               //!< update every 'n' intervals
  //@}

  
  //@{ @name Default Constructor 
  NKS_Input_Parameters() {
    Maximum_Number_of_NKS_Iterations = 0;
    Overall_Tolerance = 1e-5;      
    Relaxation_multiplier = 1.0;

    Dual_Time_Stepping = false;         
    Physical_Time_Integration = TIME_STEPPING_IMPLICIT_EULER;
    Physical_Time_CFL_Number = 1.0 ;
    Physical_Time_Step = 0.0;
    Maximum_Number_of_DTS_Steps = 0;
    Time_Accurate_Plot_Frequency = 10000000;
    Total_Physical_Time = 0.0;

    Finite_Time_Step = true;
    Finite_Time_Step_Initial_CFL = 1.0; 
    Finite_Time_Step_Final_CFL = 1.0e12;
    Finite_Time_Step_Max_CFL = 1.0e12;
     
    Maximum_Number_of_GMRES_Iterations = 0;
    GMRES_Restart = 30;
    GMRES_Overlap = 0;  
    GMRES_Initial_Tolerance = 1e-5;
    GMRES_Final_Tolerance = 1e-5;
    Normalization = true;
 
    GMRES_CHECK = false;
    GMRES_Frechet_Derivative_Order = FIRST_ORDER;   
    GMRES_Frechet_Derivative_Chem = false;
    Epsilon_Naught = 1.0e-6; // 1e-8 original, 1e-10 fixed minimum, viscous stable @ 1e-6
    
    GMRES_Block_Preconditioner = Block_Jacobi;
    Jacobian_Order = FIRST_ORDER_INVISCID_HLLE;
    GMRES_ILUK_Level_of_Fill = 0;   
    
    Detect_Convergence_Stall = false;
    DCS_Window = 15;
    output_format = OF_SCOTT;
    output_precision = 2;
    output_width = output_precision + 9;
    Freeze_Limiter_Immediately = FLI_NOT_USED;
    NKS_Write_Output_Cells_Freq = 0; 
   
    Min_Number_of_Newton_Steps_With_Zero_Limiter = 0 ;           
    Number_of_Newton_Steps_Requiring_Jacobian_Update = 1; //defaults to every step 
    Min_L2_Norm_Requiring_Jacobian_Update = 1.0e-08 ;
    Min_Finite_Time_Step_Norm_Ratio = 1.0e-10; 
    Jacobian_Update_Frequency = 1;
  };
  //@}

  ostream& Output(ostream&) const;
  int Parse_Next_Input_Control_Parameter(char *code, char *value);
  void Memory_Estimates(const int &, const int &, const int &);
  
  // Using Default Destructor 
  //~NKS_Input_Parameters() {};

  // Broadcast 
  void Broadcast_Input_Parameters() {
#ifdef _MPI_VERSION
    //Newton 
    MPI::COMM_WORLD.Bcast(&(Maximum_Number_of_NKS_Iterations), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Overall_Tolerance), 
			  1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Relaxation_multiplier), 1, MPI::DOUBLE, 0);
 
    // Dual Time Stepping
    MPI::COMM_WORLD.Bcast(&(Dual_Time_Stepping), 1, MPI::INT,    0);  //BOOL
    MPI::COMM_WORLD.Bcast(&(Physical_Time_Integration), 1, MPI::INT,   0);
    MPI::COMM_WORLD.Bcast(&(Physical_Time_CFL_Number), 1, MPI::DOUBLE, 0);  
    MPI::COMM_WORLD.Bcast(&(Physical_Time_Step), 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Maximum_Number_of_DTS_Steps), 1, MPI::INT,   0);
    MPI::COMM_WORLD.Bcast(&(Time_Accurate_Plot_Frequency), 1, MPI::INT,   0);

    // Finite Time Step
    MPI::COMM_WORLD.Bcast(&(Finite_Time_Step), 
			  1, 
			  MPI::INT, 0);  //BOOL
    MPI::COMM_WORLD.Bcast(&(Finite_Time_Step_Initial_CFL), 
			  1, 
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Finite_Time_Step_Final_CFL), 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Finite_Time_Step_Max_CFL), 
			  1, 
			  MPI::DOUBLE, 0);

    // GMRES
    MPI::COMM_WORLD.Bcast(&(Maximum_Number_of_GMRES_Iterations), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(GMRES_Restart), 
                          1, 
                          MPI::INT, 0); 
    MPI::COMM_WORLD.Bcast(&(GMRES_Overlap), 
			  1, 
			  MPI::INT, 0); 
    MPI::COMM_WORLD.Bcast(&(GMRES_Initial_Tolerance), 1, MPI::DOUBLE, 0); 
    MPI::COMM_WORLD.Bcast(&(GMRES_Final_Tolerance), 1, MPI::DOUBLE, 0); 
    MPI::COMM_WORLD.Bcast(&(Normalization), 
			  1, 
			  MPI::INT, 0); //BOOL
    MPI::COMM_WORLD.Bcast(&(GMRES_CHECK),
			  1, 
			  MPI::INT, 0); //BOOL
    MPI::COMM_WORLD.Bcast(&(GMRES_Frechet_Derivative_Order),
			  1, 
			  MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(GMRES_Frechet_Derivative_Chem),
			  1, 
			  MPI::INT, 0);  //BOOL
    MPI::COMM_WORLD.Bcast(&(Epsilon_Naught), 1, MPI::DOUBLE, 0);

    // GMRES Preconditioner Type
    MPI::COMM_WORLD.Bcast(&(GMRES_Block_Preconditioner), 
			  1, 
			  MPI::INT, 0);
    //Preconditoner Jacobian Order
    MPI::COMM_WORLD.Bcast(&(Jacobian_Order), 
                          1, 
                          MPI::INT, 0);    
    // GMRES ILUK_Level_of_Fill
    MPI::COMM_WORLD.Bcast(&(GMRES_ILUK_Level_of_Fill), 
			  1, 
			  MPI::INT, 0);

    MPI::COMM_WORLD.Bcast(&(Detect_Convergence_Stall), 1, MPI::INT, 0); // bool
    MPI::COMM_WORLD.Bcast(&(DCS_Window), 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(output_format), 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(output_precision), 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(output_width), 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Freeze_Limiter_Immediately), 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(NKS_Write_Output_Cells_Freq), 1, MPI::INT, 0);  

    MPI::COMM_WORLD.Bcast(&(Min_Number_of_Newton_Steps_With_Zero_Limiter), 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Number_of_Newton_Steps_Requiring_Jacobian_Update), 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&( Min_L2_Norm_Requiring_Jacobian_Update), 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&( Min_Finite_Time_Step_Norm_Ratio), 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Jacobian_Update_Frequency), 1, MPI::INT, 0);
#endif
  };

#ifdef _MPI_VERSION
  void Broadcast_Input_Parameters(MPI::Intracomm &Communicator,
     			          const int Source_Rank) {
    //Newton 
    Communicator.Bcast(&(Maximum_Number_of_NKS_Iterations), 
                          1, 
                          MPI::INT, Source_Rank);
    Communicator.Bcast(&(Overall_Tolerance), 
			  1, 
                          MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(Relaxation_multiplier), 1, MPI::DOUBLE, Source_Rank);
 
     // Dual Time Stepping
    Communicator.Bcast(&(Dual_Time_Stepping), 1, MPI::INT,    Source_Rank);  //BOOL
    Communicator.Bcast(&(Physical_Time_Integration), 1, MPI::INT,   Source_Rank);
    Communicator.Bcast(&(Physical_Time_CFL_Number), 1, MPI::DOUBLE, Source_Rank); 
    Communicator.Bcast(&(Physical_Time_Step), 1, MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(Maximum_Number_of_DTS_Steps), 1, MPI::INT,   Source_Rank); 
    Communicator.Bcast(&(Time_Accurate_Plot_Frequency), 1, MPI::INT,   Source_Rank);

    // Finite Time Step
    Communicator.Bcast(&(Finite_Time_Step), 
			  1, 
			  MPI::INT, Source_Rank);  //BOOL
    Communicator.Bcast(&(Finite_Time_Step_Initial_CFL), 
			  1, 
			  MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(Finite_Time_Step_Final_CFL), 1, MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(Finite_Time_Step_Max_CFL), 
			  1, 
			  MPI::DOUBLE, Source_Rank);

    // GMRES
    Communicator.Bcast(&(Maximum_Number_of_GMRES_Iterations), 
                          1, 
                          MPI::INT, Source_Rank);
    Communicator.Bcast(&(GMRES_Restart), 
                          1, 
                          MPI::INT, Source_Rank); 
    Communicator.Bcast(&(GMRES_Overlap), 
			  1, 
			  MPI::INT, Source_Rank); 
    Communicator.Bcast(&(GMRES_Initial_Tolerance), 1, MPI::DOUBLE, Source_Rank); 
    Communicator.Bcast(&(GMRES_Final_Tolerance), 1, MPI::DOUBLE, Source_Rank); 
    Communicator.Bcast(&(Normalization), 
			  1, 
			  MPI::INT, Source_Rank); //BOOL
    Communicator.Bcast(&(GMRES_CHECK),
			  1, 
			  MPI::INT, Source_Rank); //BOOL
    Communicator.Bcast(&(GMRES_Frechet_Derivative_Order),
			  1, 
			  MPI::INT, Source_Rank); 
    Communicator.Bcast(&(GMRES_Frechet_Derivative_Chem),
			  1,  
		       MPI::INT, Source_Rank); //BOOL
    Communicator.Bcast(&(Epsilon_Naught), 1, MPI::DOUBLE, Source_Rank);

    // GMRES Preconditioner Type
    Communicator.Bcast(&(GMRES_Block_Preconditioner), 
			  1, 
			  MPI::INT, Source_Rank);
    //Preconditoner Jacobian Order
    Communicator.Bcast(&(Jacobian_Order), 
                          1, 
                          MPI::INT, Source_Rank);    
    // GMRES ILUK_Level_of_Fill
    Communicator.Bcast(&(GMRES_ILUK_Level_of_Fill), 
			  1, 
			  MPI::INT, Source_Rank);

    Communicator.Bcast(&(Detect_Convergence_Stall), 1, MPI::INT, Source_Rank); // bool
    Communicator.Bcast(&(DCS_Window), 1, MPI::INT, Source_Rank);
    Communicator.Bcast(&(output_format), 1, MPI::INT, Source_Rank);
    Communicator.Bcast(&(output_precision), 1, MPI::INT, Source_Rank);
    Communicator.Bcast(&(output_width), 1, MPI::INT, Source_Rank);
    Communicator.Bcast(&(Freeze_Limiter_Immediately), 1, MPI::INT, Source_Rank);
    Communicator.Bcast(&(NKS_Write_Output_Cells_Freq), 1, MPI::INT, Source_Rank);  

    Communicator.Bcast(&(Min_Number_of_Newton_Steps_With_Zero_Limiter), 1, MPI::INT, Source_Rank);
    Communicator.Bcast(&(Number_of_Newton_Steps_Requiring_Jacobian_Update), 1, MPI::INT, Source_Rank);
    Communicator.Bcast(&(Min_L2_Norm_Requiring_Jacobian_Update), 1, MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(Min_Finite_Time_Step_Norm_Ratio), 1, MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(Jacobian_Update_Frequency), 1, MPI::INT, Source_Rank);

  };
#endif

/***************************************************************
 * NKS_Input_Parameters -- Output operator.                    *
 ***************************************************************/
 friend ostream &operator << (ostream &fout, const NKS_Input_Parameters &IP) {
     return IP.Output(fout);
 }

};

// Returns:
//  - INVALID_INPUT_VALUE if code is valid but value is invalid
//  - INVALID_INPUT_CODE  if unknown code
inline int NKS_Input_Parameters::
Parse_Next_Input_Control_Parameter(char *code, char *value)
{
  int i_command = INVALID_INPUT_CODE;
  char *ptr = NULL;

  // NEWTON 
  if (strcmp(code, "Maximum_Number_of_NKS_Iterations") == 0) {
    i_command = 64;
    Maximum_Number_of_NKS_Iterations = static_cast<int>(strtol(value, &ptr, 10));
    if (*ptr != '\0') { i_command = INVALID_INPUT_VALUE; }
  
  } else if (strcmp(code, "NKS_Overall_Tolerance") == 0) {
    i_command = 61;
    Overall_Tolerance = strtod(value, &ptr);
    if (*ptr != '\0') { i_command = INVALID_INPUT_VALUE; }
  
  } else if (strcmp(code, "Relaxation_multiplier") == 0) {
    i_command = 62;
    Relaxation_multiplier = strtod(value, &ptr);
    if (*ptr != '\0') { i_command = INVALID_INPUT_VALUE; }

     
    //DUAL TIME STEPPING 
  } else if (strcmp(code, "NKS_Dual_Time_Stepping") == 0) {
    i_command = 63;
    if (strcmp(value, "ON") == 0 || strcmp(value, "true") == 0 || strcmp(value, "1") == 0) {
      Dual_Time_Stepping = true;
    } else {
      Dual_Time_Stepping = false;
    }
   
  } else if (strcmp(code, "NKS_Physical_Time_Integration") == 0) {
    i_command = 65;
    if (strcmp(value, "Implicit_Euler") == 0) {
      Physical_Time_Integration = TIME_STEPPING_IMPLICIT_EULER;
    } else if(strcmp(value, "Second_Order_Backwards") == 0) {
      Physical_Time_Integration = TIME_STEPPING_IMPLICIT_SECOND_ORDER_BACKWARD;
    } else {
      Physical_Time_Integration = TIME_STEPPING_IMPLICIT_EULER;
    }
    
  } else if (strcmp(code, "NKS_Physical_Time_CFL") == 0) {
    i_command = 65;
    Physical_Time_CFL_Number = strtod(value, &ptr);
    if (*ptr != '\0') { i_command = INVALID_INPUT_VALUE; }
   
  } else if (strcmp(code, "NKS_Physical_Time_Step") == 0) {
    i_command = 65;
    Physical_Time_Step = strtod(value, &ptr);
    if (*ptr != '\0') { i_command = INVALID_INPUT_VALUE; }

  } else if (strcmp(code, "Maximum_Number_of_NKS_DTS_Steps" ) == 0) {
    i_command = 64;
    Maximum_Number_of_DTS_Steps = static_cast<int>(strtol(value, &ptr, 10));
    if (*ptr != '\0') { i_command = INVALID_INPUT_VALUE; }

  } else if (strcmp(code, "DTS_Time_Accurate_Plot_Frequency" ) == 0) {
    i_command = 64;
    Time_Accurate_Plot_Frequency = static_cast<int>(strtol(value, &ptr, 10));
    if (*ptr != '\0') { i_command = INVALID_INPUT_VALUE; }

    // FINITE TIME
  } else if (strcmp(code, "NKS_Finite_Time_Step") == 0) {
    i_command = 66; 
    if (strcmp(value, "OFF") == 0 || strcmp(value, "0") == 0) {
      Finite_Time_Step = false;
    } else {
      Finite_Time_Step = true;
    }

  } else if (strcmp(code, "NKS_Finite_Time_Step_Initial_CFL") == 0) {
    i_command = 70;
    Finite_Time_Step_Initial_CFL = strtod(value, &ptr);
    if (*ptr != '\0') { i_command = INVALID_INPUT_VALUE; }

  } else if (strcmp(code, "NKS_Finite_Time_Step_Final_CFL") == 0) {
    i_command = 71;
    Finite_Time_Step_Final_CFL = strtod(value, &ptr);
    if (*ptr != '\0') { i_command = INVALID_INPUT_VALUE; }

  } else if (strcmp(code, "NKS_Finite_Time_Step_Max_CFL") == 0) {
    i_command = 71;
    Finite_Time_Step_Max_CFL = strtod(value, &ptr);
    if (*ptr != '\0') { i_command = INVALID_INPUT_VALUE; }

    // GMRES 
  } else if (strcmp(code, "Maximum_Number_of_GMRES_Iterations") == 0) {
    i_command = 65;
    Maximum_Number_of_GMRES_Iterations = static_cast<int>(strtol(value, &ptr, 10));
    if (*ptr != '\0') { i_command = INVALID_INPUT_VALUE; }

  } else if (strcmp(code, "GMRES_Restart") == 0) {
    i_command = 58;
    GMRES_Restart = static_cast<int>(strtol(value, &ptr, 10));
    if (*ptr != '\0') { i_command = INVALID_INPUT_VALUE; }

  } else if (strcmp(code, "GMRES_Overlap") == 0) {
    i_command = 59;
    GMRES_Overlap = static_cast<int>(strtol(value, &ptr, 10));
    if (*ptr != '\0') { i_command = INVALID_INPUT_VALUE; }

  } else if (strcmp(code, "GMRES_Normalization") == 0) {
    i_command = 67;
    if (strcmp(value, "OFF") == 0 || strcmp(value, "0") == 0) {
      Normalization = false;
    } else {
      Normalization = true;
    }

  } else if (strcmp(code, "GMRES_Tolerance") == 0) { // backward compatible
    i_command = 60;
    GMRES_Initial_Tolerance = strtod(value, &ptr);
    GMRES_Final_Tolerance = GMRES_Initial_Tolerance; 
    if (*ptr != '\0') { i_command = INVALID_INPUT_VALUE; } 

  } else if (strcmp(code, "GMRES_Initial_Tolerance") == 0) {
    i_command = 60;
    GMRES_Initial_Tolerance = strtod(value, &ptr);
    if (*ptr != '\0') { i_command = INVALID_INPUT_VALUE; }

  } else if (strcmp(code, "GMRES_Final_Tolerance") == 0) {
    i_command = 60;
    GMRES_Final_Tolerance = strtod(value, &ptr);
    if (*ptr != '\0') { i_command = INVALID_INPUT_VALUE; }

  } else if (strcmp(code, "GMRES_Check") == 0) {
    i_command = 67;
    if (strcmp(value, "ON") == 0 || strcmp(value, "1") == 0) {
      GMRES_CHECK = true;
    } else {
      GMRES_CHECK = false;
    }

  } else if (strcmp(code, "GMRES_Frechet_Derivative_Order") == 0) {
    i_command = 62;
    if (strcmp(value, "First_Order") == 0) {
      GMRES_Frechet_Derivative_Order = FIRST_ORDER;
    } else if(strcmp(value, "Second_Order") == 0) {
      GMRES_Frechet_Derivative_Order = SECOND_ORDER;
    } else if(strcmp(value, "Second_Order_Chem") == 0) {
      GMRES_Frechet_Derivative_Order = SECOND_ORDER;
      GMRES_Frechet_Derivative_Chem = true;
    } else {
      GMRES_Frechet_Derivative_Order = FIRST_ORDER;
    }

  } else if (strcmp(code, "GMRES_Epsilon_Naught") == 0) {
    i_command = 60;
    Epsilon_Naught = strtod(value, &ptr);
    if (*ptr != '\0') { i_command = INVALID_INPUT_VALUE; }

  } else if (strcmp(code, "GMRES_Block_Preconditioner") == 0) {
    i_command = 62;
    if (strcmp(value, "Diagonal") == 0) {
      GMRES_Block_Preconditioner = Block_Jacobi;
    } else if(strcmp(value, "ILUK") == 0) {
      GMRES_Block_Preconditioner = Block_ILUK;
    } else {
      GMRES_Block_Preconditioner = Block_ILUK;
    }

  } else if (strcmp(code, "GMRES_ILUK_Level_of_Fill") == 0) {
    i_command = 63;
    GMRES_ILUK_Level_of_Fill = static_cast<int>(strtol(value, &ptr, 10));
    if (*ptr != '\0') { i_command = INVALID_INPUT_VALUE; }

  } else if (strcmp(code, "Jacobian_Order") == 0) {
    i_command = 631;
    if (strcmp(value, "Source_Terms_Only") == 0) {
      Jacobian_Order = SOURCE_TERMS_ONLY;
    } else if (strcmp(value, "First_Order_Inviscid_HLLE") == 0) {
      Jacobian_Order = FIRST_ORDER_INVISCID_HLLE;
    } else if (strcmp(value, "First_Order_Inviscid_Roe") == 0) {
      Jacobian_Order = FIRST_ORDER_INVISCID_ROE; 
    } else if (strcmp(value, "First_Order_Inviscid_AUSM_plus_up") == 0) {
      Jacobian_Order = FIRST_ORDER_INVISCID_AUSMPLUSUP;   
    } else if (strcmp(value, "Second_Order_Diamond_Path_with_HLLE") == 0) {
      Jacobian_Order = SECOND_ORDER_DIAMOND_WITH_HLLE;
    } else if (strcmp(value, "Second_Order_Diamond_Path_with_Roe") == 0) {
      Jacobian_Order = SECOND_ORDER_DIAMOND_WITH_ROE;	
    } else if (strcmp(value, "Second_Order_Diamond_Path_with_AUSM_plus_up") == 0) {
      Jacobian_Order = SECOND_ORDER_DIAMOND_WITH_AUSMPLUSUP;
    } else {
      cout << "\n***\n\nWarning: ";
      cout << "Unknown value given for Jacobian_Order. Defaulting to 1st-order HLLE.";
      cout << "\n\n***\n";
      Jacobian_Order = FIRST_ORDER_INVISCID_HLLE;
    }

  } else if (strcmp(code, "NKS_Detect_Convergence_Stall") == 0) {
    i_command = 68;
    if (strlen(value) > 1) {
	for (unsigned ii = 0; ii < strlen(value); ii++) {
	   value[ii] = tolower(value[ii]);
	}
    }
    if (strcmp(value, "on") == 0 || strcmp(value, "true") == 0 || strcmp(value, "1") == 0) {
	Detect_Convergence_Stall = true;
    } else {
	Detect_Convergence_Stall = false;
    }

  } else if (strcmp(code, "NKS_DCS_Window") == 0) {
    i_command = 69;
    DCS_Window = static_cast<int>(strtol(value, &ptr, 10));
    if (*ptr != '\0') { i_command = INVALID_INPUT_VALUE; }

  } else if (strcmp(code, "NKS_Output_Format") == 0) {
    i_command = 73;
    for (unsigned ii = 0; ii < strlen(value); ii++) {
	value[ii] = tolower(value[ii]);
    }
    if (strcmp(value, "alistair") == 0 || strcmp(value, "tight") == 0) {
	output_format = OF_ALISTAIR;
    } else {
	output_format = OF_SCOTT;
    }

  } else if (strcmp(code, "NKS_Output_Precision") == 0) {
    i_command = 74;
    output_precision = static_cast<int>(strtol(value, &ptr, 10));
    if (*ptr != '\0') { i_command = INVALID_INPUT_VALUE; }

  } else if (strcmp(code, "NKS_Output_Width") == 0) {
    i_command = 75;
    output_width = static_cast<int>(strtol(value, &ptr, 10));
    if (*ptr != '\0') { i_command = INVALID_INPUT_VALUE; }

  } else if (strcmp(code, "NKS_Freeze_Limiter_Immediately") == 0) {
    i_command = static_cast<int>(strtol(value, &ptr, 10));
    if (*ptr != '\0') { i_command = INVALID_INPUT_VALUE; }
    switch (i_command) {
	case 1: Freeze_Limiter_Immediately = FLI_NO;  break;
	case 2: Freeze_Limiter_Immediately = FLI_YES; break;
	default: Freeze_Limiter_Immediately = FLI_NOT_USED; break;
    }
    if (i_command != INVALID_INPUT_VALUE) { i_command = 76; }

  } else if (strcmp(code, "NKS_Write_Output_Cells_Freq") == 0) {
    i_command = 76;
    NKS_Write_Output_Cells_Freq = static_cast<int>(strtol(value, &ptr, 10));
    if (*ptr != '\0') { i_command = INVALID_INPUT_VALUE; }

  } else if (strcmp(code, "NKS_Min_Number_of_Newton_Steps_With_Zero_Limiter") == 0) {
    i_command = 77;
    Min_Number_of_Newton_Steps_With_Zero_Limiter = static_cast<int>(strtol(value, &ptr, 10));
    if (*ptr != '\0') { i_command = INVALID_INPUT_VALUE; }

  } else if (strcmp(code, "NKS_Number_of_Newton_Steps_Requiring_Jacobian_Update") == 0) {
    i_command = 79;
    Number_of_Newton_Steps_Requiring_Jacobian_Update = static_cast<int>(strtol(value, &ptr, 10));
    if (*ptr != '\0') { i_command = INVALID_INPUT_VALUE; }	

  } else if (strcmp(code, "NKS_Min_L2_Norm_Requiring_Jacobian_Update") == 0) {
    i_command = 80;
    Min_L2_Norm_Requiring_Jacobian_Update = strtod(value, &ptr);
    if (*ptr != '\0') { i_command = INVALID_INPUT_VALUE; }
    
  } else if (strcmp(code, "NKS_Min_Finite_Time_Step_Norm_Ratio") == 0) {
    i_command = 81;
    Min_Finite_Time_Step_Norm_Ratio = strtod(value, &ptr);
    if (*ptr != '\0') { i_command = INVALID_INPUT_VALUE; }
 
  } else if (strcmp(code, "NKS_Jacobian_Update_Frequency") == 0) {
    i_command = 77;
    Jacobian_Update_Frequency = static_cast<int>(strtol(value, &ptr, 10));
    if (*ptr != '\0') { i_command = INVALID_INPUT_VALUE; }

  } else {
    i_command = INVALID_INPUT_CODE;
  }
  return i_command;
  
}

/***************************************************************
 * NKS_Input_Parameters -- Display Output Operator             *
 ***************************************************************/
inline ostream& NKS_Input_Parameters::Output(ostream& fout) const {
  //fout.setf(ios::scientific);
  fout.unsetf(ios::scientific);
  fout << " " << endl;
  for (int star=0;star<75;star++){fout <<"*";}
  fout << "\n********                   Newton-Krylov-Schwarz                 **********" << endl;   
  for (int star=0;star<75;star++){fout<<"*";}
  
  fout <<"\n Overall Tolerance     ====> " << Overall_Tolerance << endl;

  fout <<" Relaxation Multiplier ====> " << Relaxation_multiplier << endl;

  fout <<" Dual Time Stepping    ====> ";
  if (Dual_Time_Stepping) { 
     fout << "ON\n"; 
     if ( Physical_Time_Integration == TIME_STEPPING_IMPLICIT_EULER) {
       fout<<" DTS Time Integration  ====> Implicit Euler \n";
     } else if ( Physical_Time_Integration == TIME_STEPPING_IMPLICIT_SECOND_ORDER_BACKWARD) {
       fout<<" DTS Time Integration  ====> Second Order Backwards \n";
     } 
     if( Physical_Time_Step > ZERO){
       fout <<" DTS Fixed Time Step   ====> " << Physical_Time_Step << endl; 
     } else { 
       fout <<" DTS CFL Number        ====> " << Physical_Time_CFL_Number << endl;   
     }
    fout <<" DTS Max Steps         ====> " << Maximum_Number_of_DTS_Steps  << endl;    
    fout <<" DTS Plot_Frequency    ====> " << Time_Accurate_Plot_Frequency << endl;
  } else { 
     fout << "OFF\n"; 
  }
  
  if (Finite_Time_Step == ON) {     
    fout <<" Finite Time Step      ====> ON" << endl;
    fout <<" Initial_CFL           ====> " << Finite_Time_Step_Initial_CFL << endl;
    fout <<" Final_CFL             ====> " << Finite_Time_Step_Final_CFL << endl;
    fout <<" Max_CFL               ====> " << Finite_Time_Step_Max_CFL << endl;
  } else {
    fout <<" Finite Time Step      ====> OFF" << endl; 
  } /* endif */ 
  
  // GMRES 
  fout <<" Maximum GMRES Its.    ====> " << Maximum_Number_of_GMRES_Iterations<< endl;
  fout <<" GMRES Restart Its.    ====> " << GMRES_Restart << endl;
  fout <<" Level of Overlap      ====> " << GMRES_Overlap << endl;

  if (fabs(GMRES_Initial_Tolerance-GMRES_Final_Tolerance) < 1e-10) {
    fout <<" GMRES Tolerance       ====> " << GMRES_Initial_Tolerance << endl;
  } else {
    fout <<" GMRES Initial Tol     ====> " << GMRES_Initial_Tolerance << endl;
    fout <<" GMRES Final Tol       ====> " << GMRES_Final_Tolerance << endl;
  }

  fout <<" Update Jacobian Steps ====> "<<Number_of_Newton_Steps_Requiring_Jacobian_Update<<endl;

  if (Normalization == ON) {
    fout <<" Normalization         ====> ON" << endl;
  } else {
    fout <<" Normalization         ====> OFF" << endl; 
  } 
  
  if (GMRES_CHECK) {
    fout <<" GMRES Check           ====> ON" << endl;
  }

  if( GMRES_Frechet_Derivative_Order == FIRST_ORDER) {
    fout<<   " Frechet Derivative    ====> First Order "<< endl;
  } else if ( GMRES_Frechet_Derivative_Order == SECOND_ORDER) {
    if (GMRES_Frechet_Derivative_Chem) {
      fout<<   " Frechet Derivative    ====> Second Order Chem Mod "<< endl;
    } else {
      fout<<   " Frechet Derivative    ====> Second Order "<< endl;
    }
  }

  fout <<" Matrix Free Epsilon0  ====> " << Epsilon_Naught << endl;
  
  // Precondtioner
  fout<< " Approximate Jacobian  ====>";
  if(Jacobian_Order == SOURCE_TERMS_ONLY) {
    fout<<" Source Terms Only "<<endl;
  } else if(Jacobian_Order == FIRST_ORDER_INVISCID_HLLE){
   fout<<" First Order Inviscid HLLE"<<endl;
  } else if (Jacobian_Order == FIRST_ORDER_INVISCID_ROE){
    fout<<" First Order Inviscid Roe"<<endl;
  } else if (Jacobian_Order == FIRST_ORDER_INVISCID_AUSMPLUSUP){
    fout<<" First Order Inviscid AUSM plus up"<<endl;   
  } else if (Jacobian_Order == SECOND_ORDER_DIAMOND_WITH_HLLE){
    fout<<" Second Order Diamond Path with HLLE"<<endl;
  } else if (Jacobian_Order == SECOND_ORDER_DIAMOND_WITH_ROE){
    fout<<" Second Order Diamond Path with Roe"<<endl;   
  } else if (Jacobian_Order == SECOND_ORDER_DIAMOND_WITH_AUSMPLUSUP){
    fout<<" Second Order Diamond Path with AUSM plus up"<<endl;   
  }
  
  if (GMRES_Block_Preconditioner == Block_ILUK) {        
    fout << " Local Preconditioner  ====> ILU("<< GMRES_ILUK_Level_of_Fill <<")" << endl;
  } else if (GMRES_Block_Preconditioner == Block_Jacobi){   // Diagonal
    fout << " Local Preconditioner  ====> Diagonal" << endl; 
  } 

  if (Detect_Convergence_Stall) {
    fout << " Detect Conv. Stall    ====> ON";
    fout << " (Window: " << DCS_Window << ")" << endl;
  }
  switch (Freeze_Limiter_Immediately) {
    case FLI_NO:  fout <<" Freeze Lim Immediately ===> OFF" << endl; break;
    case FLI_YES: fout <<" Freeze Lim Immediately ===> ON"  << endl; break;
    default: break;
  }

  if (NKS_Write_Output_Cells_Freq > 0) {
    fout << " Write Output Freq     ====> " << NKS_Write_Output_Cells_Freq << endl; 
  }
  
  //End 
  for (int star=0;star<75;star++){fout<<"*";}
  fout << endl;

  return fout;
}

/***************************************************************
 * NKS_Input_Parameters -- Memory Estimates                    *
 * These are NOT!! right yet, especially the PRECON            *
 ***************************************************************/
inline void NKS_Input_Parameters::Memory_Estimates(const int &blocksize, 
						   const int &block_mat_size,
						   const int &used_blocks){

  cerr<<" \n NKS_Input_Parameters::Memory_Estimates not working yet ";
//   int INT = sizeof(int);  //4 bytes
//   int DOUBLE = sizeof(double); //8bytes
//   double MB = 1024.0*1024.0;

//   double GMRES = DOUBLE*( (blocksize+1) + 
// 			  (GMRES_Restart +1) +
// 			  (GMRES_Restart*2) +
// 			  (block_mat_size*blocksize) +
// 			  (GMRES_Restart +1)*(GMRES_Restart) +
// 			  (GMRES_Restart)*(block_mat_size*blocksize) +
// 			  (GMRES_Restart+1)*(block_mat_size*blocksize) );
  
//   int nnz = 5*block_mat_size;  //ONLY FOR 1st order, 13 for 2nd
  
//   int JACOBIAN = DOUBLE*(nnz*blocksize*blocksize) +
//                  INT * ( (block_mat_size +1)*3 + nnz );

				    				    
//   //ONLY FOR ILU
//   int upper = nnz*blocksize*blocksize / 2 - block_mat_size;
//   int lower = upper;

//   int PRECON = DOUBLE*( (block_mat_size*blocksize*blocksize) +
// 			( GMRES_ILUK_Level_of_Fill + 2)*(upper) +
// 			( GMRES_ILUK_Level_of_Fill + 2)*(upper) + block_mat_size ) +
//                INT * (  (block_mat_size +1)*2 +
// 			( GMRES_ILUK_Level_of_Fill + 2)*(upper) +
// 			( GMRES_ILUK_Level_of_Fill + 2)*(upper) + block_mat_size );

  
//   //Output
//   cout<<" NKS Memory Requirement Estimate (MB)";
//   cout<<"\n Preconditioner (ILU) = "<<PRECON/MB; 
//   cout<<"\n GMRES                = "<<GMRES/MB; 
//   cout<<"\n Jacobian             = "<<JACOBIAN/MB; 
//   cout<<"\n Total                = "<<(PRECON+JACOBIAN+GMRES)/MB; 
//   cout<<endl;
//   for (int star=0;star<75;star++) {cout<<"*";}
//   cout <<endl;

} 

#endif // _NKS2DINPUT_INCLUDED

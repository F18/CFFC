/* NKSInput.h:  Header file declaring NKS input class. */

#ifndef _NKSINPUT_INCLUDED
#define _NKSINPUT_INCLUDED

// Include required C++ libraries.

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cassert>
#include <cstdlib>
#include <cstring>

using namespace std;

// Include CFD and MPI header files.

#ifndef _CFD_INCLUDED
#include "../CFD/CFD.h"
#endif // _CFD_INCLUDED

#ifndef _MPI_INCLUDED
#include "../MPI/MPI.h"
#endif // _MPI_INCLUDED

#define NKS_INPUT_PARAMETER_LENGTH 256

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
			   SECOND_ORDER = 2 };

enum Freeze_Limiter_Immediately_values { FLI_NOT_USED, FLI_NO, FLI_YES };

/*!
 * Class: NKS_Input_Parameters
 *
 * @brief Input Parameters for NKS Algorithm.
 *
 * This class defines and handles the input variables related to the
 * Newton-Krylov-Schwarz (NKS) algorithm for solving a coupled system
 * of nonlinear algebraic equations resulting from the spatial 
 * discretization of the governing system of partial differential 
 * equations of interest.  Note that unsteady calculations using a
 * desired time marching scheme can be performed using a 
 * dual-time-stepping approach.
 *
 */
class NKS_Input_Parameters{
 private:
 public:
  //@{ @name Newton parameters:
  //! Maximum number of Newton (outer) iterations for calculation
  int    Maximum_Number_of_NKS_Iterations;
  double Overall_Tolerance;
  double Relaxation_multiplier;
  //@}

  //@{ @name Dual Time Stepping Parameters
  bool Dual_Time_Stepping;              //!< Dual-time-stepping flag (on or off).  
  int Physical_Time_Integration;        //!< Implicit Euler, BDF2, ESDIRK, etc.
  double Physical_Time_CFL_Number;
  double Physical_Time_Step;
  int Maximum_Number_of_DTS_Steps; 
  //@}
  
  //@{ @name Allistar flags
  //! True/false flag for application of convergence stall detection algorithm
  bool   Detect_Convergence_Stall;
  //! Size of the window for the detection of convergence stall 
  int    DCS_Window;
  int    Freeze_Limiter_Immediately;  //!< Hack by Alistair to play with Newton convergence.
  //@}

  //@{ @name Implicit Euler parameters:
  bool   Finite_Time_Step;   
  double Finite_Time_Step_Initial_CFL;
  double Finite_Time_Step_Final_CFL; 
  double Finite_Time_Step_Max_CFL;
  //@}

  //@{ @name GMRES parameters:
  //! Maximum number of GMRES (inner) iterations for the calculation
  int    Maximum_Number_of_GMRES_Iterations;
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

  //@{ @name Matrix Free parameters:
  bool   GMRES_CHECK;
  int    GMRES_Frechet_Derivative_Order;
  double Epsilon_Naught;
  //@}

  //@{ @name Preconditioner parameters:
  int    GMRES_Block_Preconditioner; 
  int    Jacobian_Order;  
  int    GMRES_ILUK_Level_of_Fill;   
  //@}
   
  //@{ @name NKS Output parameters:
  int    Output_Format;
  int    Output_Precision, Output_Width;
  int    NKS_Write_Output_Cells_Freq; // set to zero to turn off
  //@}

  //@{ @name NKS Parameters
  int Min_Number_of_Newton_Steps_With_Zero_Limiter;            //!< force 1st order for "N" steps
  double Min_Number_of_Newton_Steps_Requiring_Jacobian_Update; //!< force Jacobian updates for "N" Newton steps       
  double Min_L2_Norm_Requiring_Jacobian_Update;                //!< force Jacobian update for L2 < "N"
  double Min_Finite_Time_Step_Norm_Ratio;                      //!< ramp over to full newton over 8 orders of L2 magnitude
  //@}

  //@{ @name Constructors and desctructors:
  //! Constructor (assign default values)
  NKS_Input_Parameters() {
    Maximum_Number_of_NKS_Iterations = 0;
    Overall_Tolerance = 1e-5;      
    Relaxation_multiplier = 1.0;
    //
    Dual_Time_Stepping = false;         
    Physical_Time_Integration = TIME_STEPPING_IMPLICIT_EULER;
    Physical_Time_CFL_Number = 1.0 ;
    Physical_Time_Step = 0.0;
    Maximum_Number_of_DTS_Steps = 0;
    //
    Finite_Time_Step = true;
    Finite_Time_Step_Initial_CFL = 1.0; 
    Finite_Time_Step_Final_CFL = 1.0e12;
    Finite_Time_Step_Max_CFL = 1.0e12;
    // 
    Maximum_Number_of_GMRES_Iterations = 0;
    GMRES_Restart = 30;
    GMRES_Overlap = 0;  
    GMRES_Initial_Tolerance = 1e-5;
    GMRES_Final_Tolerance = 1e-5;
    Normalization = true;
    //
    GMRES_CHECK = false;
    GMRES_Frechet_Derivative_Order = FIRST_ORDER; 
    Epsilon_Naught = 1e-6; // 1e-8 original, 1e-10 fixed minimum, viscous stable @ 1e-6
    //
    GMRES_Block_Preconditioner = Block_Jacobi;
    Jacobian_Order = FIRST_ORDER_INVISCID_HLLE;
    GMRES_ILUK_Level_of_Fill = 0;   
    //
    Detect_Convergence_Stall = false;
    DCS_Window = 15;
    Output_Format = 0;
    Output_Precision = 2;
    Output_Width = Output_Precision + 9;
    Freeze_Limiter_Immediately = FLI_NOT_USED;
    NKS_Write_Output_Cells_Freq = 0;  

    Min_Number_of_Newton_Steps_With_Zero_Limiter = 0 ;           
    Min_Number_of_Newton_Steps_Requiring_Jacobian_Update = 100; 
    Min_L2_Norm_Requiring_Jacobian_Update = 1.0e-08 ;
    Min_Finite_Time_Step_Norm_Ratio = 1.0e-10; 
  }

  //! Destructor
  ~NKS_Input_Parameters(void){}
  //@}

  //@{ @name Other Member functions:
  //! Broadcast input parameters to all processors
  void Broadcast(void);
  //! Parse next input line
  int Parse_Next_Input_Control_Parameter(char *code, stringstream &value);
  //! Check validity of specified input parameters
  int Check_Inputs(void);
  //! Estimate memory usage
  void Memory_Estimates(const int &, const int &, const int &);
  //@}

  //@{ @name Input-output operators:
  friend ostream &operator << (ostream &out_file,
		               NKS_Input_Parameters &IP);
  friend istream &operator >> (istream &in_file,
			       NKS_Input_Parameters &IP);
  void Output(ostream &fout) const;
  //@}

};

#endif // _NKSINPUT_INCLUDED

/*! ************************ NKS2D.h *********************************************
 *          2D Newton-Krylov-Schwarz Parallel Implicit Solver                    *  
 *                                                                               *
 * These Templated Classes & Functions are designed to work with all (well most) *
 * of the Equation systems built using the CFFC format (ie. Euler2D,             *
 * Chem2D, Dusty2D etc).                                                         *
 *                                                                               *
 *   Files:  NKS2D.h                                                             *
 *           NKSInput2D.h                                                        *
 *           GMRES2D.h                                                           *
 *           Block_Preconditioner2D.h                                            *             
 *           YourNKS.h                                                           *
 *                                                                               *
 * Explicit Specializations of certain functions and class member fuctions are   * 
 * required to tailor the NKS & GMRES to each specific equation set and are      *
 * listed below.                                                                 *
 *                                                                               *
 *      -> Newton_Update                                                         *
 *      -> Block_Preconditioner::Preconditioner_dFIdU                            *
 *      -> Block_Preconditioner::Preconditioner_dFVdU                            *
 *      -> Block_Preconditioner::normalize_Preconditioner_dFdU                   *
 *      -> GMRES_Block::set_normalize_values                                     *
 *                                                                               *
 *********************************************************************************/

#ifndef _NKS2D_INCLUDED
#define _NKS2D_INCLUDED

#include <stdlib.h>
#include <math.h>

#ifndef _BLOCK_PRECONDITONER_INCLUDED 
#include "Block_Preconditioner2D.h"
#endif 

#ifndef _GMRES2D_INCLUDED
#include "GMRES2D.h"
#endif

#ifndef _NKS_DTS_INCLUDED
#include "DTS_NKS.h"
#endif

/************************************
 * NKS "Magic" Numbers              *           //MAKE THESE INPUT PARAMETERS ???
 ************************************/ 
#define MIN_NUMBER_OF_NEWTON_STEPS_WITH_ZERO_LIMITER          0         // force 1st order for "N" steps
#define MIN_NUMBER_OF_NEWTON_STEPS_REQUIRING_JACOBIAN_UPDATE  100       // force Jacobian updates for "N" Newton steps       
#define MIN_L2_NORM_REQUIRING_JACOBIAN_UPDATE                 1.0e-08   // force Jacobian update for L2 < "N"
#define MIN_FINITE_TIME_STEP_NORM_RATIO                       1.0e-10   // ramp over to full newton over 8 orders of L2 magnitude

#define NKS_EPS 1e-18 // in case someone specifies a zero tolerance.

// These are for the Detect_Convergence_Stall algorithm.
enum DCS_States { CLEAR, WAITING_TO_FREEZE, FROZEN, WAITING_TO_STOP };


template <typename SOLN_VAR_TYPE,typename SOLN_BLOCK_TYPE, typename INPUT_TYPE>
int Newton_Update(SOLN_BLOCK_TYPE *SolnBlk,
		  AdaptiveBlock2D_List &List_of_Local_Solution_Blocks,
		  INPUT_TYPE &Input_Parameters,
		  GMRES_RightPrecon_MatrixFree<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE> &GMRES,
		  double Relaxation_multiplier);


template <typename INPUT_TYPE>
double Finite_Time_Step(const INPUT_TYPE &Input_Parameters, 
			const double &L2norm_first,
			const double &L2norm_current,
			const double &L2norm_current_n,				
			const int &Number_of_Newton_Steps);

template <typename SOLN_BLOCK_TYPE>
int set_blocksize(SOLN_BLOCK_TYPE &SolnBlk){ return (SolnBlk.NumVar()); }


/*! *****************************************************************************************
 *   Routine: Newton_Krylov_Solver
 *                                                       
 * This routine updates the specified solution block    
 * using Newton-Krylov-Schwarz method.                  
 *                                                      
 ********************************************************************************************/
template <typename SOLN_VAR_TYPE,typename SOLN_BLOCK_TYPE, typename INPUT_TYPE>
int Newton_Krylov_Schwarz_Solver(CPUTime &processor_cpu_time,
 		                 ostream &residual_file,   
		                 int &number_of_explicit_time_steps,
		                 SOLN_BLOCK_TYPE *SolnBlk,
		                 AdaptiveBlock2D_List &List_of_Local_Solution_Blocks,
		                 INPUT_TYPE &Input_Parameters) {

  int error_flag = 0;
  int blocksize = set_blocksize<SOLN_BLOCK_TYPE>(SolnBlk[0]); //Number of equations used for calculation

  //Reset Solver Type flag
  Input_Parameters.Solver_Type = IMPLICIT;

  /**************************************************************************/  
  /****************** SETTING UP STORAGE AND DATASTRUCTURES *****************/
  /**************************************************************************/  

  /* Count number of used Blocks on this processor */
  int Used_blocks_count = 0; 
  for ( int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks.Nblk; ++Bcount ) {
    if (List_of_Local_Solution_Blocks.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
      Used_blocks_count++;
    } 
  } 
  //total_Used_blocks = CFFC_Summation_MPI(Used_blocks_count);  ??      

  //DUAL TIME STEPPING 
  DTS_NKS_Quad_Block<SOLN_BLOCK_TYPE,INPUT_TYPE> *DTS_SolnBlk 
    = new DTS_NKS_Quad_Block<SOLN_BLOCK_TYPE,INPUT_TYPE>[Used_blocks_count];
  if (Input_Parameters.NKS_IP.Dual_Time_Stepping) {          
    for( int i=0; i<Used_blocks_count; i++) {  
      DTS_SolnBlk[i].allocate(SolnBlk[i].NCi,SolnBlk[i].NCj,blocksize,Input_Parameters);      
    }
  } 
  
  //GMRES 
  GMRES_RightPrecon_MatrixFree<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE> GMRES_(SolnBlk,
										List_of_Local_Solution_Blocks,
										Input_Parameters,
										blocksize,
										DTS_SolnBlk);
  // Block Preconditoner
  Block_Preconditioner<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE> *Block_precon 
    = new Block_Preconditioner<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>[Used_blocks_count];

  for( int i=0; i<Used_blocks_count; i++) {                                    
    Block_precon[i].Create_Preconditioner(SolnBlk[i],Input_Parameters,blocksize);
  }  

  /**************************************************************************/  
  /********* FINISHED SETTING UP STORAGE AND DATASTRUCTURES *****************/
  /**************************************************************************/  

  if (CFFC_Primary_MPI_Processor()){
    cout << Input_Parameters.NKS_IP; 
    //Input_Parameters.NKS_IP.Memory_Estimates(blocksize,SolnBlk[0].NCi*SolnBlk[0].NCj,Used_blocks_count);
  }

  /**************************************************************************/  
  /********* Unsteady Time Accurate, Dual Time Stepping  ********************/
  /**************************************************************************/  
  if (Input_Parameters.NKS_IP.Dual_Time_Stepping) { 
  
    int DTS_Step(1);
    double physical_time(ZERO);
    int physical_time_param(TIME_STEPPING_IMPLICIT_EULER);

    // Outer Loop (Physical Time)      
    while ( (DTS_Step < Input_Parameters.NKS_IP.Maximum_Number_of_DTS_Steps) &&
	    (Input_Parameters.Time_Max > physical_time ) ) {

      /**************************************************************************/    
      // First Step needs to be done with Implicit Euler
      if (DTS_Step == 1 ) {
	physical_time_param = Input_Parameters.NKS_IP.Physical_Time_Integration;
	Input_Parameters.NKS_IP.Physical_Time_Integration = TIME_STEPPING_IMPLICIT_EULER; 
      }
      /**************************************************************************/
      
      /**************************************************************************/
      // Determine global time step                                                        //NEED TO FIX FOR RESTARTS
      double DTS_dTime  = CFL(SolnBlk, List_of_Local_Solution_Blocks, Input_Parameters);   //ASSUMING STARTING @ TIME=0.0; 
      DTS_dTime = Input_Parameters.NKS_IP.Physical_Time_CFL_Number*CFFC_Minimum_MPI(DTS_dTime); 
  
      //Last Time sized to get Time_Max
      if( physical_time + DTS_dTime > Input_Parameters.Time_Max){
	DTS_dTime = Input_Parameters.Time_Max - physical_time;
      }
      /**************************************************************************/

      /**************************************************************************/
      // Store Previous Solution & Physical Time Step 
      for (int Bcount = 0; Bcount < List_of_Local_Solution_Blocks.Nblk; Bcount++) {
	if (List_of_Local_Solution_Blocks.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	  DTS_SolnBlk[Bcount].Store_Previous(SolnBlk[Bcount]);
	  DTS_SolnBlk[Bcount].DTS_dTime = DTS_dTime;
	}
      } 
      /**************************************************************************/

      /**************************************************************************/
      // Solve using NKS 
      error_flag = Internal_Newton_Krylov_Schwarz_Solver(processor_cpu_time,
							 residual_file,   
							 number_of_explicit_time_steps,
							 SolnBlk,
							 List_of_Local_Solution_Blocks,
							 Input_Parameters,
							 GMRES_,
							 Block_precon,
							 DTS_SolnBlk,
							 DTS_Step);              
      /**************************************************************************/

      /**************************************************************************/			
      // After first step, reset to requested Time Integration Method
      if (DTS_Step == 1) {
	Input_Parameters.NKS_IP.Physical_Time_Integration = physical_time_param;       
      }
      /**************************************************************************/

      /**************************************************************************/
      // Update Physical Time
      physical_time +=  DTS_dTime;
      /**************************************************************************/

      /**************************************************************************/
      // Error Checking 
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) { break; } 
      /**************************************************************************/

      /**************************************************************************/
      //DTS Output
      if (CFFC_Primary_MPI_Processor()) {
	cout << "\n *** End of DTS Step " << DTS_Step; 
	cout << " Time Step: " << DTS_dTime << "s  Real Time: " << physical_time << "s **** \n";
      }
      /**************************************************************************/

      // Increment DTS Steps
      DTS_Step++;
      
    } // end DTS while
    /**************************************************************************/  
  
    // Update global time steps
    number_of_explicit_time_steps += DTS_Step-1;

    
    /**************************************************************************/  
    /********* Steady State non-time accurate  ********************************/
    /**************************************************************************/  
  } else { 
    int dummy_int(0);
    error_flag = Internal_Newton_Krylov_Schwarz_Solver(processor_cpu_time,
						       residual_file,   
						       number_of_explicit_time_steps,
						       SolnBlk,
						       List_of_Local_Solution_Blocks,
						       Input_Parameters,
						       GMRES_,
						       Block_precon,
						       DTS_SolnBlk,
						       dummy_int);
  }


  /**************************************************************************/  
  // Housekeeping 
  delete[] Block_precon; 
  delete[] DTS_SolnBlk;

  /**************************************************************************/  
  //Reset Solver Type flag
  Input_Parameters.Solver_Type = EXPLICIT;

  /**************************************************************************/  
  return CFFC_OR_MPI(error_flag); 

}

/********************************************************************************************
//!  Routine: Newton_Krylov_Solver
 *                                                       
 * This routine updates the specified solution block    
 * using Newton-Krylov-Schwarz method.                  
 *                                                      
 ********************************************************************************************/
template <typename SOLN_VAR_TYPE,typename SOLN_BLOCK_TYPE, typename INPUT_TYPE>
int Internal_Newton_Krylov_Schwarz_Solver(CPUTime &processor_cpu_time,
		                          ostream &residual_file,   
		                          int &number_of_explicit_time_steps,
		                          SOLN_BLOCK_TYPE *SolnBlk,
		                          AdaptiveBlock2D_List &List_of_Local_Solution_Blocks,
		                          INPUT_TYPE &Input_Parameters,
		                          GMRES_RightPrecon_MatrixFree<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE> &GMRES_,
		                          Block_Preconditioner<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE> *Block_precon,
					  DTS_NKS_Quad_Block<SOLN_BLOCK_TYPE,INPUT_TYPE> *DTS_SolnBlk,
		                          int DTS_Step){
		     
  double *L2norm_current = new double[SolnBlk[0].Number_of_Residual_Norms]; 
  double *L1norm_current = new double[SolnBlk[0].Number_of_Residual_Norms]; 
  double *Max_norm_current = new double[SolnBlk[0].Number_of_Residual_Norms];  

  // IN case of restart need to reset this to last value for relative tolerancing to work...
  double L1norm_first     = -1.0,
         L2norm_first     = -1.0,
         Max_norm_first   = -1.0;
  double L2norm_current_n = ONE,  //ratio of first to current norm
         L1norm_current_n = ONE,
         Max_norm_current_n = ONE;

  double dTime, CFL_current;

  CPUTime total_cpu_time, start_total_cpu_time;
  double preconditioner_cputime = 0.0, res_cputime = 0.0;
  int nsetup = 0, res_nevals = 0;

  int  error_flag = 0;
  bool limiter_check= true;
  bool tell_me_about_limiter_freeze = false, 
       tell_me_about_output_tecplot = false, 
       tell_me_about_output_gmres = false;

  int Num_Var = SolnBlk[0].NumVar();  //Number of equations in "p/c" State
  int blocksize = set_blocksize<SOLN_BLOCK_TYPE>(SolnBlk[0]); //Number of equations used for calculation

  double linear_tolerance = 0.0;

  // dcs stands for detect convergence stall
  double *dcs_data = NULL;
  double dcs_slope = 0.0, dcs_L2_at_detect = 0.0, dcs_minL2 = 0.0;
  int dcs_position = 0, dcs_counter = 0; // potentially some duplication here. Let's not be too tricky.
  int dcs_counter_at_detect = 0;
  int dcs_state = CLEAR, dcs_freeze_now = 0, dcs_stop_now = 0;

  switch (Input_Parameters.i_Limiter) {
    case LIMITER_ONE:       // fall through
    case LIMITER_ZERO:      // fall through
       // TODO: make this work with using first-order for startup only.
       dcs_state = FROZEN;
       break;
  }

  if (CFFC_Primary_MPI_Processor() && 
     Input_Parameters.NKS_IP.Detect_Convergence_Stall) {
     dcs_data = new double[Input_Parameters.NKS_IP.DCS_Window];
  }

  //what ???
  int tmpp = cout.precision();

  bool GMRES_Restarted = false, GMRES_Failed = false;
  int GMRES_Restarts = 0, GMRES_Failures = 0, GMRES_Iters = 0, All_Iters_index = 0;
  int *GMRES_All_Iters = new int[Input_Parameters.NKS_IP.Maximum_Number_of_NKS_Iterations];

  double GMRES_Forcing_Power = 0.0;
  {
    double i = Input_Parameters.NKS_IP.GMRES_Initial_Tolerance;
    double f = Input_Parameters.NKS_IP.GMRES_Final_Tolerance;
    double L = max(Input_Parameters.NKS_IP.Overall_Tolerance, NKS_EPS);
    // Find x such that:
    // (i) (L)^(x) == f
    GMRES_Forcing_Power = log(f/i) / log(L);
  }
 

  /**************************************************************************/
  /******************* BEGIN NEWTON-KRYLOV-SCHWARZ CALCULATION **************/
  /**************************************************************************/
  bool NKS_continue_flag = true;    

  int Number_of_Newton_Steps = 1;  // FOR RESTART, THIS SHOULD BE SET TO LAST NKS STEP
  int i_limiter = Input_Parameters.i_Limiter;
 
  if (CFFC_Primary_MPI_Processor() && Input_Parameters.NKS_IP.output_format == OF_ALISTAIR &&
      (!Input_Parameters.NKS_IP.Dual_Time_Stepping ||
       (Input_Parameters.NKS_IP.Dual_Time_Stepping && DTS_Step == 1))) {
     int ww = Input_Parameters.NKS_IP.output_width; 
     cout.precision(Input_Parameters.NKS_IP.output_precision); 
     cout << scientific;
     cout << left;
     cout << " | NKS Outer Iteration " << endl;
     cout << setw(5) << " | " << "  " 
          << " | Residual L2-Norm (Using Entry " << Input_Parameters.i_Residual_Variable 
	  << " of the State Class)" << endl;
     cout << setw(5) << " | " << "  " << setw(ww) 
          << " | " << " | Ratio of Current to Initial Residual L2-Norm" << endl;
     cout << setw(5) << " | " << "  " << setw(ww) 
          << " | " << setw(ww) << " | " << " | CFL (for Continuation Method)" << endl;
     cout << setw(5) << " | " << "  " << setw(ww) << " | " << setw(ww) << " | " 
          << setw(ww) << " | " << "    | Number of GMRES (\"Inner\") Iterations " << endl;
     cout << setw(5) << " | " << "  " << setw(ww) << " | " 
          << setw(ww) << " | " << setw(ww) << " | " << "    | " << setw(5) << " ";
     cout << " | GMRES End Relative Residual L2-Norm times 1000" << endl;
     cout << setw(5) << " | " << "  " << setw(ww) << " | " 
          << setw(ww) << " | " << setw(ww) << " | " << "    | " << setw(5) 
          << " " << setw(ww) << " | ";
     cout << " | GMRES Initial Residual L2-Norm" << endl;
     cout << setw(5) << " | " << "  " << setw(ww) << " | " 
          << setw(ww) << " | " << setw(ww) << " | " << "    | " 
          << setw(5) << " " << setw(ww) << " | " << setw(ww) << " | ";
     cout << " | An \"F\" here means GMRES did not converge" << endl;
     cout << setw(5) << " | " << "  " << setw(ww) << " | " << setw(ww) << " | " 
          << setw(ww) << " | " << "    | " << setw(5) << " " << setw(ww) 
          << " | " << setw(ww) << " | ";
     cout << " | An \"R\" here means GMRES restarted at least once" << endl;
     cout << right;
  }

  processor_cpu_time.update();
  start_total_cpu_time.cput = CFFC_Summation_MPI(processor_cpu_time.cput); 
 
  while (NKS_continue_flag && Number_of_Newton_Steps <= Input_Parameters.NKS_IP.Maximum_Number_of_NKS_Iterations) {
    /**************************************************************************/
    // Limiter Switch to use  first order for first "N" newton steps, then switch to requested method 
    if (Number_of_Newton_Steps <= MIN_NUMBER_OF_NEWTON_STEPS_WITH_ZERO_LIMITER) {
      if (CFFC_Primary_MPI_Processor()) {  cout<<"\n Setting Limiter to ZERO, ie. Using First Order"; }
      Input_Parameters.i_Limiter = LIMITER_ZERO;   
    } else {
      Input_Parameters.i_Limiter = i_limiter;
    } 
    /**************************************************************************/
    clock_t t0 = clock();
    
    // -R(U_n) or  -R(U_n)* for DTS
    /**************************************************************************/
    /* Calculate residual: dUdt[i][j][0]  from U,W */
    for ( int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks.Nblk ; Bcount++ ) {
      if (List_of_Local_Solution_Blocks.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	error_flag = dUdt_Residual_Evaluation_NKS<SOLN_BLOCK_TYPE,INPUT_TYPE>
	                             (SolnBlk[Bcount],DTS_SolnBlk[Bcount],Input_Parameters);
      } 
    } 
    /**************************************************************************/

    /**************************************************************************/
    /* Send boundary flux corrections at block interfaces with resolution changes. */
    error_flag = Send_Conservative_Flux_Corrections(SolnBlk, 
		                		    List_of_Local_Solution_Blocks,
						    Num_Var);
    if (error_flag) {
       cout << "\n NKS ERROR: flux correction message passing error on processor "
            << List_of_Local_Solution_Blocks.ThisCPU
            << ".\n";
       cout.flush();
    } /* endif */
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);
	  
    /* Apply boundary flux corrections to residual to ensure that method is conservative. */
    Apply_Boundary_Flux_Corrections(SolnBlk, 
		                    List_of_Local_Solution_Blocks);

    res_cputime += double(clock() - t0) / double(CLOCKS_PER_SEC);
    res_nevals++;
    /**************************************************************************/

    /**************************************************************************/
    /* Calculate 1-, 2-, and max-norms of density residual (dUdt) for all blocks. */   

    L1_Norm_Residual(SolnBlk, List_of_Local_Solution_Blocks,L1norm_current);	  
    L2_Norm_Residual(SolnBlk, List_of_Local_Solution_Blocks,L2norm_current);       	  
    Max_Norm_Residual(SolnBlk, List_of_Local_Solution_Blocks,Max_norm_current);

    for(int q=0; q < SolnBlk[0].Number_of_Residual_Norms; q++){
      L1norm_current[q] = CFFC_Summation_MPI(L1norm_current[q]);      // L1 norm for all processors.
      L2norm_current[q] = sqr(L2norm_current[q]);
      L2norm_current[q] = sqrt(CFFC_Summation_MPI(L2norm_current[q])); // L2 norm for all processors.
      Max_norm_current[q] = CFFC_Maximum_MPI(Max_norm_current[q]);     // Max norm for all processors.
    }

    processor_cpu_time.update();
    total_cpu_time.cput = CFFC_Summation_MPI(processor_cpu_time.cput); 

    // Output Solution Progress 
    if (CFFC_Primary_MPI_Processor() &&
	(!Input_Parameters.NKS_IP.Dual_Time_Stepping ||
	 (Input_Parameters.NKS_IP.Dual_Time_Stepping && Number_of_Newton_Steps == 1)) ) {  //Output only 1st Step when using DTS
	int real_NKS_Step = Input_Parameters.NKS_IP.Dual_Time_Stepping ? DTS_Step : (Number_of_Newton_Steps-1);

	// FIXME - should not depend on Scott/Alistair but SolnBlk type.
	if (Input_Parameters.NKS_IP.output_format == OF_ALISTAIR) {
	   Output_Progress_to_File(residual_file, 
                                   number_of_explicit_time_steps+real_NKS_Step,
				   ZERO, 
				   total_cpu_time, 
				   L1norm_current[SolnBlk[0].residual_variable-1],
				   L2norm_current[SolnBlk[0].residual_variable-1],
				   Max_norm_current[SolnBlk[0].residual_variable-1]);
        } else {
           Output_Progress_to_File(residual_file,
		 		   number_of_explicit_time_steps+real_NKS_Step,
			           ZERO,                //DTS physical_time*THOUSAND ????
				   total_cpu_time, 
			           L1norm_current,      //maybe switch to current_n so all scale from 1 ???
			           L2norm_current,
			           Max_norm_current,
			           SolnBlk[0].residual_variable-1,
			           SolnBlk[0].Number_of_Residual_Norms);
	}
    }


    // Should we not decrease the initial (not the max) residual norm?
    //
    //  -- Alistair Tue Sep 12 2006 
    //
    if (Number_of_Newton_Steps == 1 ) {
    	L2norm_first = max(L2norm_first,L2norm_current[SolnBlk[0].residual_variable-1]);  //another restart cludge
    	L1norm_first = L1norm_current[SolnBlk[0].residual_variable-1];
    	Max_norm_first = Max_norm_current[SolnBlk[0].residual_variable-1];
    } else {
    	L2norm_first = max(L2norm_first, L2norm_current[SolnBlk[0].residual_variable-1]);
    	L1norm_first = max(L1norm_first, L1norm_current[SolnBlk[0].residual_variable-1]);
    	Max_norm_first = max(Max_norm_first, Max_norm_current[SolnBlk[0].residual_variable-1]);   
    } 

//     if (L2norm_first <= 0.0) {
//       L2norm_first = L2norm_current[SolnBlk[0].residual_variable-1];
//       L1norm_first = L1norm_current[SolnBlk[0].residual_variable-1];
//       Max_norm_first = Max_norm_current[SolnBlk[0].residual_variable-1];
//     } 

    L2norm_current_n   = L2norm_current[SolnBlk[0].residual_variable-1] / L2norm_first; 
    L1norm_current_n   = L1norm_current[SolnBlk[0].residual_variable-1] / L1norm_first; 
    Max_norm_current_n = Max_norm_current[SolnBlk[0].residual_variable-1] / Max_norm_first;

    /**************************************************************************/
    /*************** Restart & Progress ***************************************/
    /**************************************************************************/

//     // NEEDS SOME MODS TO Write and Read Restart_Solution to work properly....
//     // Periodically save restart solution files  -> !!!!!! SHOULD USE Implicit_Restart_Solution_Save_Frequency ???
//     if ( (number_of_explicit_time_steps + Number_of_Newton_Steps - 1)%Input_Parameters.Restart_Solution_Save_Frequency == 0 ) { 
//       if(CFFC_Primary_MPI_Processor()) cout << "\n\n  Saving solution to restart data file(s) after"
// 			   << " n = " << number_of_explicit_time_steps << " steps (iterations). \n";      
//       error_flag = Write_QuadTree(QuadTree,  Input_Parameters);
//       error_flag = Write_Restart_Solution(SolnBlk, 
// 					  List_of_Local_Solution_Blocks, 
// 					  Input_Parameters,
// 					  number_of_explicit_time_steps + Number_of_Newton_Steps - 1 ,
// 					  L2norm_first,
// 					  processor_cpu_time);
//       if (error_flag) {
// 	cout << "\n NKS ERROR: Unable to open restart output data file(s) on processor "
// 	     << List_of_Local_Solution_Blocks.ThisCPU<< ".\n";
// 	cout.flush();
//       } 
//       error_flag = CFFC_OR_MPI(error_flag);
//       if (error_flag) return (error_flag);  cout.flush();
//     } 

    if (Input_Parameters.NKS_IP.NKS_Write_Output_Cells_Freq > 0) {
       if (Number_of_Newton_Steps == 1 ||
	   (Number_of_Newton_Steps-1) % Input_Parameters.NKS_IP.NKS_Write_Output_Cells_Freq == 0) {
	  if (CFFC_Primary_MPI_Processor()) { tell_me_about_output_tecplot = true; }
	  do { // allows for a break
	     // TODO: this could be a true time value for DTS.
	     double no_time_please = 0.0;
	     int real_NKS_Step = Input_Parameters.NKS_IP.Dual_Time_Stepping ?
				 DTS_Step * 100 + (Number_of_Newton_Steps-1) : Number_of_Newton_Steps-1;
 	     char prefix[96], output_file_name[96];
	     ofstream output_file;    
	     int idx = 0;

	     for (idx = 0; idx < strlen(Input_Parameters.Output_File_Name); idx++) {
		if (Input_Parameters.Output_File_Name[idx] == ' ' ||
		    Input_Parameters.Output_File_Name[idx] == '.') {
		  break;
		}
	     }
   	     strncpy(prefix, Input_Parameters.Output_File_Name, idx);
	     prefix[idx] = '\0';
	     sprintf(output_file_name, "%s_n1%.4d_cells_cpu%.6d.dat",
		     prefix, real_NKS_Step, List_of_Local_Solution_Blocks.ThisCPU);
	     output_file.open(output_file_name, ios::out);
	     if (!output_file.good()) break;

	     int i_output_title = 1;
	     for (int nb = 0; nb < List_of_Local_Solution_Blocks.Nblk; nb++) {
		if (List_of_Local_Solution_Blocks.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
		   Output_Cells_Tecplot(SolnBlk[nb],
					Input_Parameters,
					real_NKS_Step, 
					no_time_please,
					List_of_Local_Solution_Blocks.Block[nb].gblknum,
					i_output_title,
			 	        output_file);
		   i_output_title = 0;
		}
	     }
	     output_file.close();

	     sprintf(output_file_name, "%s_n1%.4d_cpu%.6d.dat",
		     prefix, real_NKS_Step, List_of_Local_Solution_Blocks.ThisCPU);
	     output_file.open(output_file_name, ios::out);
	     if (!output_file.good()) break;

	     i_output_title = 1;
	     for (int nb = 0; nb < List_of_Local_Solution_Blocks.Nblk; nb++) {
		if (List_of_Local_Solution_Blocks.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
		   Output_Tecplot(SolnBlk[nb],
			          Input_Parameters,
				  real_NKS_Step, 
				  no_time_please,
				  List_of_Local_Solution_Blocks.Block[nb].gblknum,
				  i_output_title,
				  output_file);
		   i_output_title = 0;
		}
	     }
	     output_file.close();

	  } while (0); 
       }
    } // end if (Input_Parameters.NKS_IP.NKS_Write_Output_Cells_Freq > 0) 

    if (CFFC_Primary_MPI_Processor() &&	Input_Parameters.NKS_IP.Detect_Convergence_Stall) {
      double L2norm_cq = L2norm_current[SolnBlk[0].residual_variable-1]; // looks pretty on the page
      dcs_data[dcs_position++] = log10(L2norm_cq);
      if (dcs_position == Input_Parameters.NKS_IP.DCS_Window) {
	 dcs_position = 0; 
      }
      dcs_counter++;

      // Right here, dcs_counter equals the number of L2 norms
      // that we have recorded. For example, on the first
      // pass, dcs_counter is now 1.
      // And clearly from: "dcs_data[dcs_position++] = ...",
      // dcs_position now needs to equal the entry for the
      // next run.
      if ((dcs_counter == 1) || (L2norm_cq < dcs_minL2)) {
	dcs_minL2 = L2norm_cq;
      }

      if (dcs_counter >= Input_Parameters.NKS_IP.DCS_Window) {
 	 // Before calling Least_Squares_Slope(), dcs_data[dcs_position] must be 
         // the oldest entry.
	 dcs_slope = linear_regression_slope(dcs_data, 
                                             Input_Parameters.NKS_IP.DCS_Window, 
                                             dcs_position);

      //  The following are four reasons why we should transition from
      //  waiting to freeze the limiters to actually freezing the
      //  limiters. It goes something like this: "We stalled, eh? Well
      //  wait until we come close to the minimum norm encountered so far. . . .
      //  Still waiting?  Okay fine, I'll settle for half a decade
      //  (3.16) away from the minimum. . . . Still waiting? Okay
      //  fine, as long as we are a decade below the stall point then
      //  freeze. . . . Still waiting? Well, I refuse to freeze unless
      //  we go below the stall point."
      //  
      //  The same idea applies between waiting to stop and actually stopping.
#define ST1(thisL2, minL2, detectL2, steps) \
        (  (thisL2) < (1.1*(minL2))  )
#define ST2(thisL2, minL2, detectL2, steps) \
	(  (steps) > 10  && (thisL2) < (3.16*(minL2))  )
#define ST3(thisL2, minL2, detectL2, steps) \
	(  (steps) > 15  && (thisL2) < (0.10*(detectL2))  )
#define ST4(thisL2, minL2, detectL2, steps) \
	(  (steps) > 20  && (thisL2) < (     (detectL2))  )

      switch (dcs_state) {
	case CLEAR: 
	  if (dcs_slope > 0.0) {
	    cout << " Convergence stall detected. Will freeze limiters when next appropriate.\n";
	    dcs_state++;
	    dcs_L2_at_detect = L2norm_cq;
	    dcs_counter_at_detect = dcs_counter;
	  }
	  break;
	case WAITING_TO_FREEZE:
	  if (ST1(L2norm_cq, dcs_minL2, dcs_L2_at_detect, dcs_counter - dcs_counter_at_detect) ||
	      ST2(L2norm_cq, dcs_minL2, dcs_L2_at_detect, dcs_counter - dcs_counter_at_detect) ||
	      ST3(L2norm_cq, dcs_minL2, dcs_L2_at_detect, dcs_counter - dcs_counter_at_detect) ||
	      ST4(L2norm_cq, dcs_minL2, dcs_L2_at_detect, dcs_counter - dcs_counter_at_detect)) {
 	    cout << " Freezing limiters to help with convergence stall (";
	    if (ST1(L2norm_cq, dcs_minL2, dcs_L2_at_detect, dcs_counter - dcs_counter_at_detect)) { 
               cout << "ST1";
	    } else if (ST2(L2norm_cq, dcs_minL2, dcs_L2_at_detect, dcs_counter - dcs_counter_at_detect)) { 
               cout << "ST2";
	    } else if (ST3(L2norm_cq, dcs_minL2, dcs_L2_at_detect, dcs_counter - dcs_counter_at_detect)) { 
               cout << "ST3";
	    } else if (ST4(L2norm_cq, dcs_minL2, dcs_L2_at_detect, dcs_counter - dcs_counter_at_detect)) { 
               cout << "ST4";
	    }
	    cout << ").\n";
	    dcs_state++;
	    dcs_counter = 0;
	    dcs_position = 0;
	    dcs_freeze_now = 1;
	  }
 	  break;
	case FROZEN:
	  if (dcs_slope > 0.0) {
	    cout <<" Convergence stall detected even with limiters frozen. Will stop when next appropriate.\n";
	    dcs_state++;
	    dcs_L2_at_detect = L2norm_cq;
	    dcs_counter_at_detect = dcs_counter;
	  }
	  break;
        case WAITING_TO_STOP:
	  if (ST1(L2norm_cq, dcs_minL2, dcs_L2_at_detect, dcs_counter - dcs_counter_at_detect) ||
	      ST2(L2norm_cq, dcs_minL2, dcs_L2_at_detect, dcs_counter - dcs_counter_at_detect) ||
	      ST3(L2norm_cq, dcs_minL2, dcs_L2_at_detect, dcs_counter - dcs_counter_at_detect) ||
	      ST4(L2norm_cq, dcs_minL2, dcs_L2_at_detect, dcs_counter - dcs_counter_at_detect)) {
	     cout << " Limiters are frozen and yet convergence has stalled nonetheless. Stopping. (";
	     if (ST1(L2norm_cq, dcs_minL2, dcs_L2_at_detect, dcs_counter - dcs_counter_at_detect)) { 
                cout << "ST1";
	     } else if (ST2(L2norm_cq, dcs_minL2, dcs_L2_at_detect, dcs_counter - dcs_counter_at_detect)) { 
                cout << "ST2";
	     } else if (ST3(L2norm_cq, dcs_minL2, dcs_L2_at_detect, dcs_counter - dcs_counter_at_detect)) { 
                cout << "ST3";
	     } else if (ST4(L2norm_cq, dcs_minL2, dcs_L2_at_detect, dcs_counter - dcs_counter_at_detect)) { 
                cout << "ST4";
	     }
	     cout << ".)\n";
	     dcs_stop_now = 1;
	  }
	  break;
        }

      }// if (dcs_counter >= window_size) 
    } // if (Input_Parameters.NKS_IP.Detect_Convergence_Stall)

    CFFC_Broadcast_MPI(&dcs_stop_now, 1);
    if (dcs_stop_now) { break; }
    CFFC_Broadcast_MPI(&dcs_freeze_now, 1);

    if (limiter_check && Number_of_Newton_Steps > 1) {
       // These two if-statements should really be one but 
       // are split up to make it more readable.
       if ((Input_Parameters.Freeze_Limiter && 
            L2norm_current[SolnBlk[0].residual_variable-1] <= 
            Input_Parameters.Freeze_Limiter_Residual_Level) ||
	   (Input_Parameters.NKS_IP.Freeze_Limiter_Immediately == FLI_YES) || 
	    (dcs_freeze_now)) {
	  // Should one use the absolute (L2norm_current) or relative
	  // (L2norm_current_n) for when to freeze the limiter? Relative
	  // seems desirable since it appears more problem independent.
	  // But what if you use another scheme for startup? Well then
	  // you could pass the initial L2 norm from that to this
	  // function. Okay, what if you start with full multigrid? Ahh,
	  // you think about it.
	  //   -- Alistair Wood Tue Sep 19 2006
	  //
	  // This point becomes somewhat moot with the introduction of 
	  // my convergence stall detection algorithm.
	  //  -- Alistair Wood Sun Apr 01 2007
	  if (CFFC_Primary_MPI_Processor()) {
	    switch (Input_Parameters.NKS_IP.output_format) {
	      case OF_SCOTT:
		cout << "\n\n ********** Apply Limiter Freezing ********** \n";
		break;
	      case OF_ALISTAIR:
		if (!dcs_freeze_now) { 
	 	   tell_me_about_limiter_freeze = true;
		}
		break;
	      default: 
                break;
	    }
	  } 
	  Freeze_Limiters(SolnBlk, List_of_Local_Solution_Blocks);	
	  limiter_check = false;	
       } 
    } // if (limiter_check && Number_of_Newton_Steps > 1)
    if (dcs_freeze_now) { dcs_freeze_now = 0; }

    /**************************************************************************/
    /***************** NEWTON STEP ********************************************/
    /**************************************************************************/
    if (L2norm_current_n > Input_Parameters.NKS_IP.Overall_Tolerance){  
       
      /**************************************************************************/     
      if (Input_Parameters.NKS_IP.Finite_Time_Step) {  
	// Apply finite time step, ie. Implicit Euler ( as deltaT -> inf. Newton's Method)
	CFL_current = Finite_Time_Step(Input_Parameters,L2norm_first,
				       L2norm_current[SolnBlk[0].residual_variable-1] ,
				       L2norm_current_n, Number_of_Newton_Steps);	   
      } else { 
	CFL_current = Input_Parameters.NKS_IP.Finite_Time_Step_Final_CFL;   	
      } 

      // Determine time steps
      dTime = CFL(SolnBlk, List_of_Local_Solution_Blocks, Input_Parameters);
      
      if (Input_Parameters.Local_Time_Stepping == GLOBAL_TIME_STEPPING) {
	dTime = CFFC_Minimum_MPI(dTime); 
	dTime *= CFL_current;           
	Set_Global_TimeStep(SolnBlk, List_of_Local_Solution_Blocks, dTime);      
      } else {
	for (int Bcount = 0; Bcount < List_of_Local_Solution_Blocks.Nblk; Bcount++) {
	  if (List_of_Local_Solution_Blocks.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	    for (int i = 0; i < SolnBlk[Bcount].NCi; i++) {
	      for (int j = 0; j < SolnBlk[Bcount].NCj; j++) {
		SolnBlk[Bcount].dt[i][j] *= CFL_current;
	      }
	    }
	  }
	}
      }
      /**************************************************************************/
      
      /**************************************************************************/
      // Print out Newton Step info at the beginning of the iteration. 
      if (CFFC_Primary_MPI_Processor() && 
	  (!Input_Parameters.NKS_IP.Dual_Time_Stepping || Number_of_Newton_Steps ==1 )) {      	
	switch (Input_Parameters.NKS_IP.output_format ) {
	  case OF_SCOTT: 
	    cout.precision(10);
	    cout << "\n Newton Step (Outer It.) = " << Number_of_Newton_Steps;
	    cout << " L2norm = " << L2norm_current[SolnBlk[0].residual_variable-1];
	    cout << " L2norm_ratio = " << L2norm_current_n;
	    cout << " CFL = " << CFL_current;
	    //if(Input_Parameters.Preconditioning) cout<<" Mref = "<<Mrefnew;
	    break;
	  case OF_ALISTAIR: {
	    cout.precision(Input_Parameters.NKS_IP.output_precision); 
	    int ww = Input_Parameters.NKS_IP.output_width; 
	    cout << setw(5) << Number_of_Newton_Steps << "  ";
	    cout << setw(ww) << L2norm_current[SolnBlk[0].residual_variable-1];
	    cout << setw(ww) << L2norm_current_n;
	    if (Input_Parameters.NKS_IP.Dual_Time_Stepping) {
	       cout << setw(ww) << " ";
	    } else {
	       cout << setw(ww) << CFL_current;
	    }
	    }
            break;
	  default: 
            break;
	}
      }
      /**************************************************************************/

      /**************************************************************************/
      // Store  Uo, dt for use in GMRES & Precondtioner Calculations
      for ( int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks.Nblk ; ++Bcount ) {
	if (List_of_Local_Solution_Blocks.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {	  
	  /* Copy solutions in conserved variables, U to Uo for all blocks for update procedure. */
	  for (int j = SolnBlk[Bcount].JCl-SolnBlk[Bcount].Nghost;  j <= SolnBlk[Bcount].JCu+SolnBlk[Bcount].Nghost; j++){
	    for (int i = SolnBlk[Bcount].ICl-SolnBlk[Bcount].Nghost;  i <= SolnBlk[Bcount].ICu+SolnBlk[Bcount].Nghost; i++){
	      SolnBlk[Bcount].Uo[i][j] = SolnBlk[Bcount].U[i][j];
	    } 
	  } 	  
	} 
      }  
      /**************************************************************************/

      //DTS RULES: update Preconditioner only if GMRES iterations increase between Newton Its.
      bool GMRES_Iters_increaseing(true);
      if( GMRES_Iters != 0){
	if(GMRES_Iters ==1) {
	  GMRES_Iters_increaseing = false;
	} else if (GMRES_All_Iters[All_Iters_index-2] > GMRES_Iters) {
	  GMRES_Iters_increaseing = true;
	} else {
	  GMRES_Iters_increaseing = false;
	}
      }

      /**************************************************************************/
      /************* PRECONDTIONER "BLOCK" JACOBIANS ****************************/      
      /**************************************************************************/
      // Create/Update Jacobian Matrix(s) using Uo = U  
      if ( ( !Input_Parameters.NKS_IP.Dual_Time_Stepping &&
	     (Number_of_Newton_Steps < MIN_NUMBER_OF_NEWTON_STEPS_REQUIRING_JACOBIAN_UPDATE || 
	      L2norm_current_n > MIN_L2_NORM_REQUIRING_JACOBIAN_UPDATE) ) ||                          
	   ( Input_Parameters.NKS_IP.Dual_Time_Stepping && GMRES_Iters_increaseing) ) {

        if (CFFC_Primary_MPI_Processor() ) {	
	    switch (Input_Parameters.NKS_IP.output_format) {
	      case OF_SCOTT: 
                cout << "\n Creating/Updating Jacobian Matrix"; 
                break;
	      case OF_ALISTAIR: 
                cout << "    | "; 
                break;
	      default: 
                break;
	    }
	}
	clock_t t0 = clock();
	//Update for each block.
	for ( int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks.Nblk ; ++Bcount ) {
	  if (List_of_Local_Solution_Blocks.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	    // Jacobian changed and Preconditioner updated 
	    Block_precon[Bcount].Update_Jacobian_and_Preconditioner(DTS_SolnBlk[Bcount].DTS_dTime);
	  } 
	} 
	preconditioner_cputime += double(clock() - t0) / double(CLOCKS_PER_SEC);
	nsetup++;
      } else {
	if (CFFC_Primary_MPI_Processor() && 
            Input_Parameters.NKS_IP.output_format == OF_ALISTAIR) {
	cout << " noJ| "; 
	}
      }
      /**************************************************************************/

      /**************************************************************************/
      /************* LINEAR SYSTEM SOLVE WITH GMRES  ****************************/      
      /**************************************************************************/
      /* Solve system with right-preconditioned matrix free GMRES */  
      linear_tolerance = min(0.9, Input_Parameters.NKS_IP.GMRES_Initial_Tolerance*
			     pow(L2norm_current_n, GMRES_Forcing_Power));

      error_flag = GMRES_.solve(Block_precon, linear_tolerance,
				&GMRES_Restarted, &GMRES_Failed, &GMRES_Iters,
				&res_cputime, &res_nevals);
      if (CFFC_Primary_MPI_Processor() && error_flag) cout << "\n NKS2D: Error in GMRES \n";
      if (GMRES_Restarted) { GMRES_Restarts++; }
      if (GMRES_Failed)    { GMRES_Failures++; }
      GMRES_All_Iters[All_Iters_index++] = GMRES_Iters;
      /**************************************************************************/      
 
      /**************************************************************************/      
      // Solution Update U = Uo + GMRES.deltaU
      error_flag = Newton_Update<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>
	           (SolnBlk,List_of_Local_Solution_Blocks,
                    Input_Parameters,GMRES_, 
                    Input_Parameters.NKS_IP.Relaxation_multiplier);
      if (CFFC_Primary_MPI_Processor() && error_flag) cout <<  "\n NKS2D: Error in Solution Update \n";
      /**************************************************************************/

      if (Input_Parameters.NKS_IP.NKS_Write_Output_Cells_Freq > 0) {
        if (Number_of_Newton_Steps == 1 || (Number_of_Newton_Steps-1) % Input_Parameters.NKS_IP.NKS_Write_Output_Cells_Freq == 0) {
 	   if (CFFC_Primary_MPI_Processor()) { tell_me_about_output_gmres = true; }
  	   GMRES_.Output_GMRES_vars_Tecplot(Number_of_Newton_Steps-1,
					    L2norm_current[SolnBlk[0].residual_variable-1], L2norm_current_n);
	}			
      } 

      /**************************************************************************/
      /* Exchange solution information between neighbouring blocks. */    
      error_flag = Send_All_Messages(SolnBlk,
				     List_of_Local_Solution_Blocks,
				     Num_Var, 
				     OFF);
      if (error_flag) {
         cout << "\n 2D_NKS ERROR: 2D message passing error on processor "
	      << List_of_Local_Solution_Blocks.ThisCPU
	      << ".\n";
         cout.flush();
      } /* endif */
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return (error_flag);
      /**************************************************************************/

      /**************************************************************************/  
      /* Apply boundary conditions for Newton step. */
      BCs(SolnBlk,List_of_Local_Solution_Blocks,Input_Parameters);
      /**************************************************************************/

      /**************************************************************************/
      Number_of_Newton_Steps++;      

      // END OF  if (L2norm_current_n > Input_Parameters.NKS_IP.Overall_Tolerance)
      /**************************************************************************/
    } else {       
      // L2norm_current_n < Input_Parameters.NKS_IP.Overall_Tolerance so set flag to "stop"     
      NKS_continue_flag = false;
    } /* endif */
    /**************************************************************************/
 
    if (CFFC_Primary_MPI_Processor() && 
        Input_Parameters.NKS_IP.output_format == OF_ALISTAIR) { 
       if (tell_me_about_limiter_freeze) {
	  cout << "  Limiter Frozen "; 
	  tell_me_about_limiter_freeze = false;
       } 
       if (tell_me_about_output_tecplot) {
          cout << " Output_Tecplot()";
          tell_me_about_output_tecplot = false;
       }
       if (tell_me_about_output_gmres) {
          cout << " Output_GMRES()";
          tell_me_about_output_gmres = false;
       }
       cout << endl; 
    }

  }  // END OF NEWTON ITERATION WHILE LOOP

  /**************************************************************************/  
  /********* FINISHED NEWTON KRYLOV SCHWARZ *********************************/
  /**************************************************************************/  
  
  /**************************************************************************/
  /* Reset limiter. */
  Input_Parameters.i_Limiter = i_limiter;
  /**************************************************************************/
 
  processor_cpu_time.update();
  total_cpu_time.cput = CFFC_Summation_MPI(processor_cpu_time.cput); 
  double total_preconditioner_cputime = CFFC_Summation_MPI(preconditioner_cputime);
  double total_res_cputime = CFFC_Summation_MPI(res_cputime);
  int quick_used = 0;
  for (int blk = 0; blk < List_of_Local_Solution_Blocks.Nblk; ++blk) {
     if (List_of_Local_Solution_Blocks.Block[blk].used == ADAPTIVEBLOCK2D_USED) {
        quick_used++;
     } 
  }
  int quick_total_used = CFFC_Summation_MPI(quick_used);

  /**************************************************************************/
  /* Output final 2-norm for all blocks */ 
  if (CFFC_Primary_MPI_Processor() ) {  //&& !Input_Parameters.NKS_IP.Dual_Time_Stepping
    int ttmmpp = cout.precision();
    int nns = Number_of_Newton_Steps-1;
    cout << "\n ";
    for (int star=0;star<75;star++){cout <<"*";}
    //cout.setf(ios::scientific);
    cout << "\n End of Newton Steps = " << nns  << " L2norm = "
         << L2norm_current[SolnBlk[0].residual_variable-1] 
         << " L2norm_ratio = " << L2norm_current_n << "\n";
    //cout.unsetf(ios::scientific);
    //cout.setf(ios::fixed);
    int GMRESt = 0;
    for (int iii = 0; iii < nns; iii++) { 
       GMRESt += GMRES_All_Iters[iii]; 
    }
    cout << " GMRES Total Iterations: " << GMRESt ;
    cout << "   GMRES Restarts: " << GMRES_Restarts ;
    cout << "   GMRES Failures: " << GMRES_Failures ;
    /////////////////////////////////////////////////////////////////
    if (Input_Parameters.NKS_IP.output_format == OF_ALISTAIR) {
      qsort(GMRES_All_Iters, nns, sizeof(int), compare_integers);
      cout << "\nGMRES Iterations:" << endl;
      for (int iii = 0; iii < nns; iii++) {	
         cout << " " << setw(3) << GMRES_All_Iters[iii];
	 if (((iii+1) % 20) == 0) { 
            cout << endl; 
         }
      }
      cout << endl;
      cout.precision(1);
      // Suppose, for example, that there are five entries in 
      // GMRES_All_Iters (GAI):
      //    20, 22, 24, 26, 28
      // Then we would want to say:
      // 80% <= 26 where the avg was 23
      // Here, then, "nns" is 5, n_consider is 4 (5*0.8) and the
      // for-loop for the average considers entries GAI[0] to GAI[3].
      // 
      // For non-integer cases, say 90% for this example, round 
      // (and not floor or ceil) to the nearest entry; here, 
      // n_consider would round up to 5 (round(5*0.9)).
      for (double iii = 60.0; iii < 101.0; iii += 10.0) {
         int n_consider = (int) round(nns*iii/100.0); 
         double aavg = 0.0;
         for (int jjj = 0; jjj < n_consider; jjj++) {
	    aavg += (double) GMRES_All_Iters[jjj];
	 }
         aavg /= (double) n_consider;
         cout << setw(3) << (int) iii << " percent <= " << setw(3) << GMRES_All_Iters[n_consider - 1];
         cout << " where the avg was " << setw(5) << aavg << endl;
      }
      cout << endl;
      double NKS_sec = total_cpu_time.sec() - start_total_cpu_time.sec();
      int aa = SolnBlk[0].ICu - SolnBlk[0].ICl + 1, bb = SolnBlk[0].JCu - SolnBlk[0].JCl + 1;
      cout.precision(1);
      int gqw = 8;
      double time_else = NKS_sec - total_res_cputime - total_preconditioner_cputime;
      cout << "Total NKS cputime:            " << setw(gqw) << NKS_sec << "  s (" 
           << NKS_sec / 60.0 << " minutes)\n";
      cout << "Time spent on residual evals: " << setw(gqw) << total_res_cputime 
           << "  s (" << total_res_cputime/NKS_sec*100.0 << " %)\n";
      cout << "Time spent on preconditioners:" << setw(gqw) << total_preconditioner_cputime 
           << "  s (" << total_preconditioner_cputime/NKS_sec*100.0 << " %)\n";
      cout << "Time spent elsewhere:         " << setw(gqw) << time_else << "  s (" 
           << time_else/NKS_sec*100.0 << " %)\n";
      cout << "\n";
      double qq = total_res_cputime/double(res_nevals)*1000.0;
      cout << "Time for one residual evaluation:" << setw(gqw) << qq 
           << "  ms (average over " << res_nevals << " calls)\n";
      qq /= double(quick_total_used);
      cout << "One residual per block:          " << setw(gqw) << qq 
           << "  ms (" << quick_total_used << " blocks)\n";
      qq *= 1000.0 / double(aa*bb);
      cout << "One residual per cell:           " << setw(gqw) << qq << "  us (" 
           << aa << "x" << bb << " = " << aa*bb << " cells per block)\n";
      cout << "\n";
      qq = total_preconditioner_cputime/double(nsetup)*1000.0;
      cout << "Time for one complete precond:   " << setw(gqw) << qq 
           << "  ms (average over " << nsetup << " setups)\n";
      qq /= double(quick_total_used);
      cout << "One complete precond per block:  " << setw(gqw) << qq 
           << "  ms (" << quick_total_used << " blocks)\n";
      qq *= 1000.0 / double(aa*bb);
      cout << "One complete precond per cell:   " << setw(gqw) << qq 
           << "  us (" << aa << "x" << bb << " = " << aa*bb << " cells per block)\n";
    } // if (Input_Parameters.NKS_IP.output_format == OF_ALISTAIR)

    cout <<"\n ";
    for (int star=0;star<75;star++){
       cout <<"*";
    }
    cout.precision(ttmmpp);
    cout.unsetf(ios::fixed);
  } /* endif */
  /**************************************************************************/

  /**************************************************************************/
  if (!Input_Parameters.NKS_IP.Dual_Time_Stepping) {
     number_of_explicit_time_steps += Number_of_Newton_Steps-1;
  } // Otherwise the external solver will update the number of explicit steps.
  /**************************************************************************/
    
  /**************************************************************************/
  // Housekeeping 
  delete[] L2norm_current; delete[] L1norm_current; delete[] Max_norm_current;
  delete[] GMRES_All_Iters;
  if (CFFC_Primary_MPI_Processor() && 
      Input_Parameters.NKS_IP.Detect_Convergence_Stall) {
    delete [] dcs_data;
  }
  cout.precision(tmpp);
  /**************************************************************************/
  return (error_flag);

} /* End of Newton_Krylov_Schwarz_Solver. */




/*! *****************************************************************************************
 * Generic Newton_Update                                                                    *
 *                                                                                          *
 * This routine updates the previous Solution Data Uo with the deltaU from the GMRES        *                  
 * iterative solver.  U = Uo + GMRES.delatU                                                 *
 *                                                                                          *
 ********************************************************************************************/
template <typename SOLN_VAR_TYPE,typename SOLN_BLOCK_TYPE, typename INPUT_TYPE>
int Newton_Update(SOLN_BLOCK_TYPE *SolnBlk,
		  AdaptiveBlock2D_List &List_of_Local_Solution_Blocks,
		  INPUT_TYPE &Input_Parameters,
		  GMRES_RightPrecon_MatrixFree<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE> &GMRES,
		  double Relaxation_multiplier) {

  int Num_Var = SolnBlk[0].NumVar();  
  
  /* Update Solution. No updates to Ghost Cells, let the BC's take care of it */
  for ( int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks.Nblk ; ++Bcount ) {
    if (List_of_Local_Solution_Blocks.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
      for (int j = SolnBlk[Bcount].JCl; j <= SolnBlk[Bcount].JCu; j++){
	for (int i = SolnBlk[Bcount].ICl; i <= SolnBlk[Bcount].ICu; i++){
	  
	  /* Update solutions in conversed variables  U = Uo + deltaU = Uo + denormalized(x) */	 
	  for(int varindex =1; varindex <= Num_Var; varindex++){  
	    SolnBlk[Bcount].U[i][j][varindex] = SolnBlk[Bcount].Uo[i][j][varindex] 
	      +  Relaxation_multiplier * GMRES.deltaU(Bcount,i,j,varindex-1);
	  } 	      	  
	  // THIS FUNCTION HAS NO CHECKS FOR INVALID SOLUTIONS, 
	  // YOU PROBABLY WANT TO CREATE A SPECIALIZATION OF THIS FUNCTION SPECIFIC 
	  // FOR YOUR EQUATION SYSTEM see Euler2D, Chem2D, etc...
	 
	  //Update solution in primitive variables.
	  SolnBlk[Bcount].W[i][j] = W(SolnBlk[Bcount].U[i][j]);	  
	} 
      } 
    } 
  } 
  
  cout<<"\n USING GENERIC NEWTON_UPDATE \n";

  return 0; 
}



/*! *****************************************************************************************
 * Generic Finite_Time_Step                                                                 *
 *                                                                                          *
 * This routine calculates the Finite Time Step using a basic SER approach for startup,     *                  
 * returning a "CFL" number to multiply by the stability determined dt.                     *
 *                                                                                          *
 ********************************************************************************************/
template <typename INPUT_TYPE>
double Finite_Time_Step(const INPUT_TYPE &Input_Parameters, 
			const double &L2norm_first,
			const double &L2norm_current,
			const double &L2norm_current_n,		
			const int &Number_of_Newton_Steps) {

  double CFL_current;

  //SER 
  // Knoll and Keyes (JCP 193 (2004) 357-397) write that SER means that the
  // change in time step is based on the change in the non-linear residual norm
  // over the last two Newton iterations, that is, the time step is not based
  // simply on the current norm residual ratio. 
  //   -- Alistair Wood. Wed Aug 08 2007.
  if (L2norm_current_n > MIN_FINITE_TIME_STEP_NORM_RATIO ) { 
    CFL_current = Input_Parameters.NKS_IP.Finite_Time_Step_Initial_CFL*
                  pow( max(ONE, ONE/L2norm_current_n),ONE );
                  //pow(min(ONE, max(ONE, ONE/L2norm_current_n)*MIN_FINITE_TIME_STEP_NORM_RATIO),ONE );     
  } else {
     CFL_current = Input_Parameters.NKS_IP.Finite_Time_Step_Initial_CFL/
                   MIN_FINITE_TIME_STEP_NORM_RATIO;
  } 
 
  return CFL_current;

}

#endif  /* _NKS2D_INCLUDED */

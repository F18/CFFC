/*! ********************** NKS2D.h ***********************************************
 *          2D Newton-Krylov-Schwarz Parallel Implicit Solver                    *  
 *                                                                               *
 * These Templated Classes & Functions are designed to work with all (well most) *
 * of the Equation systems built using the CFDkit+caboodle format (ie. Euler2D,  *
 * Chem2D, Dusty2D etc).                                                         *
 *                                                                               *
 *   Files:  NKS2D.h                                                             *
 *           NKSInput2D.h                                                        *
 *           GMRES2D.h                                                           *
 *           Block_Preconditioner2D.h                                            *             
 *           Your_Specializations.h                                              *
 *                                                                               *
 * Explicit Specializations of certain functions and class member fuctions are   * 
 * required to tailor the NKS & GMRES to each specific equation set and are      *
 * listed below.                                                                 *
 *                                                                               *
 *      -> Newton_Update                                                         *
 *      -> Block_Preconditioner::Preconditioner_dFdU                             *
 *      -> Block_Preconditioner::normalize_Preconditioner_dFdU                   *
 *      -> GMRES_Block::set_normalize_values                                     *
 *      -> GMRES_Block::GMRES_BCs                                                *
 *                                                                               *
 *********************************************************************************/

#ifndef _NKS2D_INCLUDED
#define _NKS2D_INCLUDED

#ifndef _BLOCK_PRECONDITONER_INCLUDED 
#include "Block_Preconditioner2D.h"
#endif 

#ifndef _GMRES2D_INCLUDED
#include "GMRES2D.h"
#endif 

/************************************
 * NKS "Magic" Numbers              *           //MAKE THESE INPUT PARAMETERS ???
 ************************************/ 
#define MIN_NUMBER_OF_NEWTON_STEPS_WITH_ZERO_LIMITER          0         // force 1st order for "N" steps
#define MIN_NUMBER_OF_NEWTON_STEPS_REQUIRING_JACOBIAN_UPDATE  100       // force Jacobian updates for "N" Newton steps       
#define MIN_L2_NORM_REQUIRING_JACOBIAN_UPDATE                 1.0e-08   // force Jacobian update for L2 < "N"
#define MIN_FINITE_TIME_STEP_NORM_RATIO                       1.0e-10   // ramp over to full newton over 8 orders of L2 magnitude

template <typename SOLN_VAR_TYPE,typename SOLN_BLOCK_TYPE, typename INPUT_TYPE>
int Newton_Update(SOLN_BLOCK_TYPE *SolnBlk,
		  AdaptiveBlock2D_List &List_of_Local_Solution_Blocks,
		  INPUT_TYPE &Input_Parameters,
		  GMRES_RightPrecon_MatrixFree<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE> &GMRES);

template <typename INPUT_TYPE>
double Finite_Time_Step(const INPUT_TYPE &Input_Parameters, 
			const double &L2norm_first,
			const double &L2norm_current,
			const double &L2norm_current_n,				
			const int &Number_of_Newton_Steps);

template <typename SOLN_BLOCK_TYPE>
int set_blocksize(SOLN_BLOCK_TYPE &SolnBlk){ return (SolnBlk.NumVar()); }

/********************************************************************************************
//!  Routine: Newton_Krylov_Solver
 *                                                       
 * This routine updates the specified solution block    
 * using Newton-Krylov-Schwarz method.                  
 *                                                      
 ********************************************************************************************/
template <typename SOLN_VAR_TYPE,typename SOLN_BLOCK_TYPE, typename INPUT_TYPE>
int Newton_Krylov_Schwarz_Solver(CPUTime &NKS_processor_cpu_time,
				 CPUTime &NKS_total_cpu_time,
				 ostream &residual_file,   
 				 int &number_of_explicit_time_steps,
				 SOLN_BLOCK_TYPE *SolnBlk,
				 AdaptiveBlock2D_List &List_of_Local_Solution_Blocks,
				 INPUT_TYPE &Input_Parameters) { 
 

  double *L2norm_current = new double[SolnBlk[0].Number_of_Residual_Norms]; 
  double *L1norm_current = new double[SolnBlk[0].Number_of_Residual_Norms]; 
  double *Max_norm_current = new double[SolnBlk[0].Number_of_Residual_Norms];  

  // IN case of restart need to reset this to last value for relative tolerancing to work...
  double L2norm_first     = ZERO, 
         L1norm_first     = ZERO,
         Max_norm_first   = ZERO;
  double L2norm_current_n = ONE,  //ratio of first to current norm
         L1norm_current_n = ONE,
         Max_norm_current_n = ONE;
  double dTime, CFL_current;

  int  error_flag = 0;
  bool limiter_check= true;
  int Used_blocks_count = 0; 

  int Num_Var = SolnBlk[0].NumVar();  //Number of equations in "p/c" State
  int blocksize = set_blocksize<SOLN_BLOCK_TYPE>(SolnBlk[0]); //Number of equations used for calculation

  //Reset Solver Type flag
  Input_Parameters.Solver_Type = IMPLICIT;

  cout.setf(ios::scientific);

  /**************************************************************************/  
  /****************** SETTING UP STORAGE AND DATASTRUCTURES *****************/
  /**************************************************************************/  
  /* Count number of used Blocks on this processor */
  for ( int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks.Nblk; ++Bcount ) {
    if (List_of_Local_Solution_Blocks.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
      Used_blocks_count++;
    } 
  } 
  
  //SETUP GMRES Blocks 
  GMRES_RightPrecon_MatrixFree<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE> GMRES_(SolnBlk,
										List_of_Local_Solution_Blocks,
										Input_Parameters,
										blocksize);

  //SETUP Preconditioner Blocks
  Block_Preconditioner<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE> *Block_precon 
    = new Block_Preconditioner<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>[Used_blocks_count];

  for( int i=0; i<Used_blocks_count; i++) {  
    Block_precon[i].Create_Preconditioner(SolnBlk[i], Input_Parameters,blocksize);
  }  
  //total_Used_blocks = CFDkit_Summation_MPI(Used_blocks_count);  maybe ??      

  /**************************************************************************/  
  /********* FINISHED SETTING UP STORAGE AND DATASTRUCTURES *****************/
  /**************************************************************************/  
  if (CFDkit_Primary_MPI_Processor()){
    Input_Parameters.NKS_IP.Output(); 
    //Input_Parameters.NKS_IP.Memory_Estimates(blocksize,SolnBlk[0].NCi*SolnBlk[0].NCj,Used_blocks_count);
  }
   
  /**************************************************************************/
  /******************* BEGIN NEWTON-KRYLOV-SCHWARZ CALCULATION **************/
  /**************************************************************************/
  bool NKS_continue_flag = true;    
//   bool Update_Jacobian_flag = true;  //what's the deal with this flag, I think it is basically useless....

  int Number_of_Newton_Steps = 1;  // FOR RESTART, THIS SHOULD BE SET TO LAST NKS STEP
  int i_limiter = Input_Parameters.i_Limiter;
 
  while ( NKS_continue_flag && Number_of_Newton_Steps <= Input_Parameters.NKS_IP.Maximum_Number_of_NKS_Iterations) {

    /**************************************************************************/
    // Limiter Switch to use  first order for first "N" newton steps, then switch to requested method 
    if (Number_of_Newton_Steps <= MIN_NUMBER_OF_NEWTON_STEPS_WITH_ZERO_LIMITER) {
      if (CFDkit_Primary_MPI_Processor()) {  cout<<"\n Setting Limiter to ZERO, ie. Using First Order"; }
      Input_Parameters.i_Limiter = LIMITER_ZERO;   
    } else {
      Input_Parameters.i_Limiter = i_limiter;
    } 
    /**************************************************************************/
    
//     //THIS IS CHEM2D SPECIFIC SHOULD BE MOVED TO Chem2DQuad_Specializations.h !!!!!
//     /************** Change Low-Mach Number Preconditioning Mref ***************/ 
//     if(Input_Parameters.Preconditioning && Mrefnew > Input_Parameters.Mach_Number_Reference_target ){      
//       if(L2norm_current_n < 0.01){
// 	//start at 1 and decrease to Mref set in input
// 	Mrefnew = Input_Parameters.Mach_Number_Reference - 
// 	  (Input_Parameters.Mach_Number_Reference - Input_Parameters.Mach_Number_Reference_target)*(log10(ONE/L2norm_current_n)/FOUR);     
// 	if( Mrefnew < Input_Parameters.Mach_Number_Reference_target) Mrefnew = Input_Parameters.Mach_Number_Reference_target;
// 	}
//       Change_Mref(SolnBlk, List_of_Local_Solution_Blocks,Mrefnew);
//     }	
//     /**************************************************************************/

    // -R(U_n)
    /**************************************************************************/
    /* Calculate residual: dudt[i][j][0]  from U,W */
    for ( int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks.Nblk ; Bcount++ ) {
      if (List_of_Local_Solution_Blocks.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	error_flag = dUdt_Residual_Evaluation(SolnBlk[Bcount],Input_Parameters);
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
    error_flag = CFDkit_OR_MPI(error_flag);
    if (error_flag) return (error_flag);
	  
    /* Apply boundary flux corrections to residual to ensure that method is conservative. */
    Apply_Boundary_Flux_Corrections(SolnBlk, 
		                    List_of_Local_Solution_Blocks);
    /**************************************************************************/

    /**************************************************************************/
    /* Calculate 1-, 2-, and max-norms of density residual (dUdt) for all blocks. */   

    L1_Norm_Residual(SolnBlk, List_of_Local_Solution_Blocks,L1norm_current);	  
    L2_Norm_Residual(SolnBlk, List_of_Local_Solution_Blocks,L2norm_current);       	  
    Max_Norm_Residual(SolnBlk, List_of_Local_Solution_Blocks,Max_norm_current);

    for(int q=0; q < SolnBlk[0].Number_of_Residual_Norms; q++){
      L1norm_current[q] = CFDkit_Summation_MPI(L1norm_current[q]);      // L1 norm for all processors.
      L2norm_current[q] = sqr(L2norm_current[q]);
      L2norm_current[q] = sqrt(CFDkit_Summation_MPI(L2norm_current[q])); // L2 norm for all processors.
      Max_norm_current[q] = CFDkit_Maximum_MPI(Max_norm_current[q]);     // Max norm for all processors.
    }

    if (Number_of_Newton_Steps == 1 ) {
      L2norm_first = max(L2norm_first,L2norm_current[SolnBlk[0].residual_variable-1]);  //another restart cludge
      L1norm_first = L1norm_current[SolnBlk[0].residual_variable-1];
      Max_norm_first = Max_norm_current[SolnBlk[0].residual_variable-1];
    } else {
      L2norm_first = max(L2norm_first, L2norm_current[SolnBlk[0].residual_variable-1]);
      L1norm_first = max(L1norm_first, L1norm_current[SolnBlk[0].residual_variable-1]);
      Max_norm_first = max(Max_norm_first, Max_norm_current[SolnBlk[0].residual_variable-1]);   
    } 

    L2norm_current_n   = L2norm_current[SolnBlk[0].residual_variable-1] / L2norm_first; 
    L1norm_current_n   = L1norm_current[SolnBlk[0].residual_variable-1] / L1norm_first; 
    Max_norm_current_n = Max_norm_current[SolnBlk[0].residual_variable-1] / Max_norm_first;

    /* Calculate ratio of initial and current 2-norms. */
    /**************************************************************************/

    /**************************************************************************/
    /*************** Restart & Progress ***************************************/
    /**************************************************************************/
    NKS_processor_cpu_time.update();
    // Total CPU time for all processors.
    NKS_total_cpu_time.cput = CFDkit_Summation_MPI(NKS_processor_cpu_time.cput); 

//     // NEEDS SOME MODS TO Write and Read Restart_Solution to work properly....
//     // Periodically save restart solution files  -> !!!!!! SHOULD USE Implicit_Restart_Solution_Save_Frequency ???
//     if ( (number_of_explicit_time_steps + Number_of_Newton_Steps - 1)%Input_Parameters.Restart_Solution_Save_Frequency == 0 ) {       
//       if(CFDkit_Primary_MPI_Processor()) cout << "\n\n  Saving solution to restart data file(s) after"
// 			   << " n = " << number_of_explicit_time_steps << " steps (iterations). \n";      
//       error_flag = Write_QuadTree(QuadTree,  Input_Parameters);
//       error_flag = Write_Restart_Solution(SolnBlk, 
// 					  List_of_Local_Solution_Blocks, 
// 					  Input_Parameters,
// 					  number_of_explicit_time_steps + Number_of_Newton_Steps - 1 ,
// 					  L2norm_first,
// 					  NKS_processor_cpu_time);        //Should add explicit processor time as well????
//       if (error_flag) {
// 	cout << "\n NKS ERROR: Unable to open restart output data file(s) on processor "
// 	     << List_of_Local_Solution_Blocks.ThisCPU<< ".\n";
// 	cout.flush();
//       } 
//       error_flag = CFDkit_OR_MPI(error_flag);
//       if (error_flag) return (error_flag);  cout.flush();
//     } 

    // Output progress information for the calculation
    if (CFDkit_Primary_MPI_Processor()){
      Output_Progress_to_File(residual_file,
			      number_of_explicit_time_steps+Number_of_Newton_Steps-1,  
			      ZERO,
			      NKS_total_cpu_time,     //Should add explicit processor time as well????
			      L1norm_current,         //maybe switch to current_n so all scale from 1 ???
			      L2norm_current,
			      Max_norm_current,
			      SolnBlk[0].residual_variable-1,
			      SolnBlk[0].Number_of_Residual_Norms);
    }
    /**************************************************************************/
 
    /**************************************************************************/
    // Freeze Limiters if Residual less than given value or # of orders of reduction.
    if (Input_Parameters.Freeze_Limiter && limiter_check){
      if (Number_of_Newton_Steps > 1 &&  L2norm_current[SolnBlk[0].residual_variable-1] 
	  <= Input_Parameters.Freeze_Limiter_Residual_Level)  {    // absolute
	//if (Number_of_Newton_Steps > 1 &&  L2norm_current_n <= Input_Parameters.Freeze_Limiter_Residual_Level)  {  // relative 
	if (CFDkit_Primary_MPI_Processor()) {
	  cout << "\n\n ********** Apply Limiter Freezing ********** \n";
        } 
	Freeze_Limiters(SolnBlk, List_of_Local_Solution_Blocks);	
	limiter_check = false;	
      } 
    } 
    /**************************************************************************/

    /**************************************************************************/
    /***************** NEWTON STEP ********************************************/
    /**************************************************************************/
    if (L2norm_current_n > Input_Parameters.NKS_IP.Overall_Tolerance){  
     
      /**************************************************************************/
      /************************** TIME STEP *************************************/ 
      /**************************************************************************/
      /* Calculate delta t = min(delta t) among processors. */
      dTime = CFL(SolnBlk, List_of_Local_Solution_Blocks, Input_Parameters);  //deltaTref
      dTime = CFDkit_Minimum_MPI(dTime); 
   
      // Apply finite time step, ie. Implicit Euler ( as deltaT -> inf. Newtons Method)
      if (Input_Parameters.NKS_IP.Finite_Time_Step) {  
	CFL_current = Finite_Time_Step(Input_Parameters,L2norm_first,
				       L2norm_current[SolnBlk[0].residual_variable-1] ,
				       L2norm_current_n, Number_of_Newton_Steps);	   
      } else { 
        CFL_current = Input_Parameters.NKS_IP.Finite_Time_Step_Initial_CFL;   	
      }  

      //Set all dt[i][j] for global time stepping if requested (shouldn't be as this is NOT Time accurate)
      if(Input_Parameters.Local_Time_Stepping == GLOBAL_TIME_STEPPING ){
	Set_Global_TimeStep(SolnBlk, List_of_Local_Solution_Blocks, dTime);      
      }
      /**************************************************************************/
   

      /**************************************************************************/
      // Print out Newtwon Step info at the beginning of the iteration. 
      cout.precision(10);
      if (CFDkit_Primary_MPI_Processor()) {      	
	cout << "\n Newton Step (Outer It.) = " << Number_of_Newton_Steps << " L2norm = "
	     << L2norm_current[SolnBlk[0].residual_variable-1] << " L2norm_ratio = " << L2norm_current_n 
	     << " CFL = " << CFL_current <<" min_deltat = "<<CFL_current*dTime;
	//if(Input_Parameters.Preconditioning) cout<<" Mref = "<<Mrefnew;
      } 
      /**************************************************************************/

   
      /**************************************************************************/
      // Store  Uo, dt for use in GMRES & Precondtioner Calculations
      for ( int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks.Nblk ; ++Bcount ) {
	if (List_of_Local_Solution_Blocks.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {	  
	  /* Copy solutions in conserved variables, U to Uo for all blocks for update procedure. */
	  for (int j = SolnBlk[Bcount].JCl-SolnBlk[Bcount].Nghost;  j <= SolnBlk[Bcount].JCu+SolnBlk[Bcount].Nghost; j++){
	    for (int i = SolnBlk[Bcount].ICl-SolnBlk[Bcount].Nghost;  i <= SolnBlk[Bcount].ICu+SolnBlk[Bcount].Nghost; i++){
	      for(int varindex =1; varindex <= Num_Var; varindex++){		
		SolnBlk[Bcount].Uo[i][j][varindex] = SolnBlk[Bcount].U[i][j][varindex];
	      }    	      
	      /* Set time step. */
              SolnBlk[Bcount].dt[i][j] = CFL_current*SolnBlk[Bcount].dt[i][j];
	    } 
	  } 	  
	} 
      }  
      /**************************************************************************/

 
      /**************************************************************************/
      /************* PRECONDTIONER "BLOCK" JACOBIANS ****************************/      
      /**************************************************************************/
      // Create/Update Jacobian Matrix(s) using Uo = U  
      if (Number_of_Newton_Steps < MIN_NUMBER_OF_NEWTON_STEPS_REQUIRING_JACOBIAN_UPDATE
	   || L2norm_current_n > MIN_L2_NORM_REQUIRING_JACOBIAN_UPDATE ) { // && Update_Jacobian_flag ) { 
	
	if (CFDkit_Primary_MPI_Processor()) cout << "\n Creating/Updating Jacobian Matrix";  
	
	//Update for each block.
	for ( int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks.Nblk ; ++Bcount ) {
	  if (List_of_Local_Solution_Blocks.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	    // Jacobian changed and Preconditioner updated 
	    Block_precon[Bcount].Update_Jacobian_and_Preconditioner();
	  } 
	} 
      }
 //      else { 
// 	// don't update the Jacobian any more, but this is a stupid use of a flag 
// 	// as the Number of Iterations of tolerance must be violated first anyway. ????
// 	Update_Jacobian_flag == false;
//       } /* endif */ 
      /**************************************************************************/


      /**************************************************************************/
      /************* LINEAR SYSTEM SOLVE WITH GMRES  ****************************/      
      /**************************************************************************/
      /* Solve system with right-preconditioned matrix free GMRES */  
      error_flag = GMRES_.solve(Block_precon, Number_of_Newton_Steps);
      if (CFDkit_Primary_MPI_Processor() && error_flag) cout << "\n NKS2D: Error in GMRES \n";
      /**************************************************************************/      
 

      /**************************************************************************/      
      // Solution Update U = Uo + GMRES.deltaU
      error_flag = Newton_Update<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE>
	(SolnBlk,List_of_Local_Solution_Blocks,Input_Parameters,GMRES_);      
      if (CFDkit_Primary_MPI_Processor() && error_flag) cout <<  "\n NKS2D: Error in Solution Update \n";
      /**************************************************************************/


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
      error_flag = CFDkit_OR_MPI(error_flag);
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
 

  }  // END OF NEWTON ITERATION WHILE LOOP
 
  /**************************************************************************/  
  /********* FINISHED NEWTON KRYLOV SCHWARZ *********************************/
  /**************************************************************************/  
  
  /**************************************************************************/
  /* Reset limiter. */
  Input_Parameters.i_Limiter = i_limiter;
  /**************************************************************************/

  /**************************************************************************/
  /* Output final 2-norm for all blocks */ 
  if (CFDkit_Primary_MPI_Processor()) {     
    cout << " " << endl;
    for (int star=0;star<75;star++){cout <<"*";}
    cout << "\nEnd of Newton Steps = " << Number_of_Newton_Steps-1  << " L2norm = "
	 << L2norm_current[SolnBlk[0].residual_variable-1] << " L2norm_ratio = " << L2norm_current_n << endl;
    for (int star=0;star<75;star++){cout <<"*";}
  } /* endif */
  /**************************************************************************/

  /**************************************************************************/
  number_of_explicit_time_steps = number_of_explicit_time_steps + Number_of_Newton_Steps-1;
  /**************************************************************************/
    
  /**************************************************************************/
  // Housekeeping 
  delete[] Block_precon; 
  delete[] L2norm_current; delete[] L1norm_current; delete[] Max_norm_current;
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
		  GMRES_RightPrecon_MatrixFree<SOLN_VAR_TYPE,SOLN_BLOCK_TYPE,INPUT_TYPE> &GMRES) {

  int Num_Var = SolnBlk[0].NumVar();  
  
  /* Update Solution. No updates to Ghost Cells, let the BC's take care of it */
  for ( int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks.Nblk ; ++Bcount ) {
    if (List_of_Local_Solution_Blocks.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
      for (int j = SolnBlk[Bcount].JCl; j <= SolnBlk[Bcount].JCu; j++){
	for (int i = SolnBlk[Bcount].ICl; i <= SolnBlk[Bcount].ICu; i++){
	  
	  /* Update solutions in conversed variables  U = Uo + deltaU = Uo + denormalized(x) */	 
	  for(int varindex =1; varindex <= Num_Var; varindex++){  
	    SolnBlk[Bcount].U[i][j][varindex] = SolnBlk[Bcount].Uo[i][j][varindex] 
	      +  GMRES.deltaU(Bcount,i,j,varindex-1);
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
  if (L2norm_current_n > MIN_FINITE_TIME_STEP_NORM_RATIO ) { 
    CFL_current = Input_Parameters.NKS_IP.Finite_Time_Step_Initial_CFL*
      pow( max(ONE, ONE/L2norm_current_n),ONE ); 
      //      pow(min(ONE, max(ONE, ONE/L2norm_current_n)*MIN_FINITE_TIME_STEP_NORM_RATIO),ONE );     
  } else {
     CFL_current = Input_Parameters.NKS_IP.Finite_Time_Step_Initial_CFL/MIN_FINITE_TIME_STEP_NORM_RATIO;
  } 
 
  return CFL_current;

}

#endif  /* _NKS2D_INCLUDED */

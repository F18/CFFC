#ifndef _NKS_INCLUDED
#define _NKS_INCLUDED
/*! ********************** NKS2D.h ***********************************************
 *          3D Newton-Krylov-Schwarz Parallel Implicit Solver                    *  
 *                                                                               *
 * These Templated Classes & Functions are designed to work with all (well most) *
 * of the Equation systems built using the CFCC format                           *
 *                                                                               *
 *   Files:  NKS.h                                                               *
 *           NKSInput.h                                                          *
 *           GMRES.h                                                             *
 *           Block_Preconditioner.h                                              *             
 *           Your_Specializations.h                                              *
 *                                                                               *
 * Explicit Specializations of certain functions and class member fuctions are   * 
 * required to tailor the NKS & GMRES to each specific equation set and are      *
 * listed below.                                                                 *
 *                                                                               *
 *      -> Hexa_Newton_Krylov_Schwarz_Solver::Newton_Update                      *
 *      -> Block_Preconditioner::Preconditioner_dFdU                             *
 *      -> Block_Preconditioner::normalize_Preconditioner_dFdU                   *
 *      -> GMRES_Block::set_normalize_values                                     *
 *                                                                               *
 *********************************************************************************/

#include "GMRES.h"

/*! **************************************************************
 * class: Hexa_Newton_Krylov_Schwarz_Solver                      *
 *                                                               *
 *                                                               *
 * Data Members:                                                 *
 *                                                               *
 *****************************************************************/ 
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
class Hexa_Newton_Krylov_Schwarz_Solver {
  private:

  int blocksize;
  //Overload this to change blocksize (ie Species N-1 etc. )
  void set_blocksize(void) { blocksize = Solution_Data->Local_Solution_Blocks.Soln_Blks[0].NumVar(); }
 
  double *L2norm_current, *L1norm_current, *Max_norm_current; 
  double L2norm_first, L1norm_first, Max_norm_first;  
  double L2norm_current_n, L1norm_current_n, Max_norm_current_n;

  CPUTime total_cpu_time, start_total_cpu_time; //better names specific to NKS..

  protected: 
  public:

  //Data
  HexaSolver_Data *Data;
  HexaSolver_Solution_Data<SOLN_pSTATE, SOLN_cSTATE> *Solution_Data;

  DTS_Hexa_Block<SOLN_pSTATE,SOLN_cSTATE> *DTS_SolnBlk; 
  GMRES_RightPrecon_MatrixFree<SOLN_pSTATE,SOLN_cSTATE> GMRES;
  Block_Preconditioner<SOLN_pSTATE,SOLN_cSTATE> *Block_precon; 

  // Constructor
  Hexa_Newton_Krylov_Schwarz_Solver(void): Data(NULL), Solution_Data(NULL), DTS_SolnBlk(NULL) {} 
  Hexa_Newton_Krylov_Schwarz_Solver(HexaSolver_Data &Data_ptr, 
				    HexaSolver_Solution_Data<SOLN_pSTATE, SOLN_cSTATE> &Solution_Data_ptr);

  void allocate();
  void deallocate(); 
  // Destructor
  ~Hexa_Newton_Krylov_Schwarz_Solver() {deallocate();}

  //Member Functions
  void Calculate_Norms(const int &);
  int Restart_NKS();
  int Newton_Update();      
  double Finite_Time_Step(const int &); 

  int Solve();
  int Steady_Solve(const double &,const int &);

};


/*! *****************************************************************************************
 *   Routine: 
 ********************************************************************************************/
template <typename SOLN_pSTATE, typename SOLN_cSTATE> 
Hexa_Newton_Krylov_Schwarz_Solver<SOLN_pSTATE,SOLN_cSTATE>:: 
Hexa_Newton_Krylov_Schwarz_Solver(HexaSolver_Data &Data_ptr, 
				  HexaSolver_Solution_Data<SOLN_pSTATE, SOLN_cSTATE> &Solution_Data_ptr){
  Data = &(Data_ptr);
  Solution_Data = &(Solution_Data_ptr);  
  set_blocksize();
  allocate();
}

/*! *****************************************************************************************
 *   Routine: 
 ********************************************************************************************/
template <typename SOLN_pSTATE, typename SOLN_cSTATE> 
void Hexa_Newton_Krylov_Schwarz_Solver<SOLN_pSTATE,SOLN_cSTATE>:: 
allocate(){

  /* Count number of used Blocks on this processor */
  int Used_blocks_count = 0; 
  for ( int Bcount = 0 ; Bcount < Data->Local_Adaptive_Block_List.Nblk; ++Bcount ) {
    if (Data->Local_Adaptive_Block_List.Block[Bcount].used == ADAPTIVEBLOCK3D_USED) {
      Used_blocks_count++;
    } 
  }  //total_Used_blocks = CFFC_Summation_MPI(Used_blocks_count);  ??      
 
  //DUAL TIME STEPPING 
  DTS_SolnBlk = new DTS_Hexa_Block<SOLN_pSTATE,SOLN_cSTATE>[Used_blocks_count];
  if ( Solution_Data->Input.NKS_IP.Dual_Time_Stepping ) {          
    for( int i=0; i< Used_blocks_count; i++) {  
      DTS_SolnBlk[i].allocate(Solution_Data->Local_Solution_Blocks.Soln_Blks[i].NCi,
			      Solution_Data->Local_Solution_Blocks.Soln_Blks[i].NCj,
			      Solution_Data->Local_Solution_Blocks.Soln_Blks[i].NCk);
    }
  } 

  //GMRES 
  GMRES.Setup(*Data, *Solution_Data, *DTS_SolnBlk, blocksize);

  // Block Preconditoner
  Block_precon = new Block_Preconditioner<SOLN_pSTATE,SOLN_cSTATE>[Used_blocks_count];
  for( int i=0; i<Used_blocks_count; i++) {                                    
    Block_precon[i].Create_Preconditioner(Solution_Data->Local_Solution_Blocks.Soln_Blks[i],
					  Solution_Data->Input,
					  blocksize);
  }  

  //Norms 
  L2norm_current = new double[Solution_Data->Input.Number_of_Residual_Norms]; 
  L1norm_current = new double[Solution_Data->Input.Number_of_Residual_Norms]; 
  Max_norm_current = new double[Solution_Data->Input.Number_of_Residual_Norms];  

}

/*! *****************************************************************************************
 *   Routine: 
 ********************************************************************************************/
template <typename SOLN_pSTATE, typename SOLN_cSTATE> 
void Hexa_Newton_Krylov_Schwarz_Solver<SOLN_pSTATE,SOLN_cSTATE>:: 
deallocate() { 
  if(Block_precon != NULL)     delete[] Block_precon; 
  if(DTS_SolnBlk != NULL )     delete[] DTS_SolnBlk; 
  if(L2norm_current != NULL)   delete[] L2norm_current; 
  if(L1norm_current != NULL)   delete[] L1norm_current; 
  if(Max_norm_current != NULL) delete[] Max_norm_current;
} 

/*! *****************************************************************************************
 *   Routine: 
 ********************************************************************************************/
template <typename SOLN_pSTATE, typename SOLN_cSTATE> 
void Hexa_Newton_Krylov_Schwarz_Solver<SOLN_pSTATE,SOLN_cSTATE>::
Calculate_Norms(const int &Number_of_Newton_Steps){

  // Calculate 1-, 2-, and max-norms of density residual (dUdt) for all blocks. 
  
  // NEED TO MOD THIS FOR  Number_of_Residual_Norms >1
  L1norm_current[Solution_Data->Input.p_Norm_Indicator-1] = Solution_Data->Local_Solution_Blocks.L1_Norm_Residual();       
  L2norm_current[Solution_Data->Input.p_Norm_Indicator-1] = Solution_Data->Local_Solution_Blocks.L2_Norm_Residual();
  Max_norm_current[Solution_Data->Input.p_Norm_Indicator-1] = Solution_Data->Local_Solution_Blocks.Max_Norm_Residual();
  
//   for(int q=0; q < Solution_Data->Input.Number_of_Residual_Norms; q++){
//     L1norm_current[q] = CFFC_Summation_MPI(L1norm_current[q]);      // L1 norm for all processors.
//     L2norm_current[q] = sqr(L2norm_current[q]);
//     L2norm_current[q] = sqrt(CFFC_Summation_MPI(L2norm_current[q])); // L2 norm for all processors.
//     Max_norm_current[q] = CFFC_Maximum_MPI(Max_norm_current[q]);     // Max norm for all processors.
//   }
  
  //Relative Norms used for CFL scaling during startup
  if (Number_of_Newton_Steps == 1 ) {
    L2norm_first = L2norm_current[Solution_Data->Input.p_Norm_Indicator-1];
      //max(L2norm_first,L2norm_current[Solution_Data->Input.p_Norm_Indicator-1]);  //another restart cludge
    L1norm_first = L1norm_current[Solution_Data->Input.p_Norm_Indicator-1];
    Max_norm_first = Max_norm_current[Solution_Data->Input.p_Norm_Indicator-1];
  } else {
    L2norm_first = max(L2norm_first, L2norm_current[Solution_Data->Input.p_Norm_Indicator-1]);
    L1norm_first = max(L1norm_first, L1norm_current[Solution_Data->Input.p_Norm_Indicator-1]);
    Max_norm_first = max(Max_norm_first, Max_norm_current[Solution_Data->Input.p_Norm_Indicator-1]);   
  } 
  
  // Calculate ratio of initial and current 2-norms. 
  L2norm_current_n   = L2norm_current[Solution_Data->Input.p_Norm_Indicator-1] / L2norm_first; 
  L1norm_current_n   = L1norm_current[Solution_Data->Input.p_Norm_Indicator-1] / L1norm_first; 
  Max_norm_current_n = Max_norm_current[Solution_Data->Input.p_Norm_Indicator-1] / Max_norm_first;

}

/*! *****************************************************************************************
 *   Routine: 
 ********************************************************************************************/
template <typename SOLN_pSTATE, typename SOLN_cSTATE> 
int Hexa_Newton_Krylov_Schwarz_Solver<SOLN_pSTATE,SOLN_cSTATE>::
Restart_NKS(){
  int error_flag(0);

//     // NEEDS SOME MODS TO Write and Read Restart_Solution to work properly....
//     // Periodically save restart solution files  -> !!!!!! SHOULD USE Implicit_Restart_Solution_Save_Frequency ???
//     if ( (number_of_explicit_time_steps + Number_of_Newton_Steps - 1)%Solution_Data->Input.Restart_Solution_Save_Frequency == 0 ) {       
//       if(CFFC_Primary_MPI_Processor()) cout << "\n\n  Saving solution to restart data file(s) after"
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
//       error_flag = CFFC_OR_MPI(error_flag);
//       if (error_flag) return (error_flag);  cout.flush();
//     } 

  return error_flag;
}

/*! *****************************************************************************************
 * Generic Newton_Update                                                                    *
 *                                                                                          *
 * This routine updates the previous Solution Data Uo with the deltaU from the GMRES        *                  
 * iterative solver.  U = Uo + GMRES.delatU                                                 *
 *                                                                                          *
 ********************************************************************************************/
template <typename SOLN_pSTATE, typename SOLN_cSTATE> 
int Hexa_Newton_Krylov_Schwarz_Solver<SOLN_pSTATE,SOLN_cSTATE>:: 
Newton_Update(){

  /* Update Solution. No updates to Ghost Cells, let the BC's take care of it */
    for ( int Bcount = 0 ; Bcount < Data->Local_Adaptive_Block_List.Nblk; ++Bcount ) {      
      if (Data->Local_Adaptive_Block_List.Block[Bcount].used) {

	for (int k = Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].KCl; 
	     k <= Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].KCu; k++){
	  for (int j = Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].JCl; 
	       j <= Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].JCu; j++){
	    for (int i = Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].ICl; 
		 i <= Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].ICu; i++){

	      cout<<"\n Uo["<<i<<","<<j<<","<<k<<"] = "<<Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].Uo[i][j][k];
	      cout<<"\n dU["<<i<<","<<j<<","<<k<<"] = ";

	      /* Update solutions in conversed variables  U = Uo + deltaU = Uo + denormalized(x) */	 
	      for(int varindex =1; varindex <= Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].NumVar() ; varindex++){  
	      Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].U[i][j][k][varindex] = 
		Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].Uo[i][j][k][varindex] 
		+  GMRES.deltaU(Bcount,i,j,k,varindex-1);
	      cout<<" "<<GMRES.deltaU(Bcount,i,j,k,varindex-1);
	      } 	      	  

	      cout<<"\n U["<<i<<","<<j<<","<<k<<"] = "<<Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].U[i][j][k];

	      // THIS FUNCTION HAS NO CHECKS FOR INVALID SOLUTIONS, 
	      // YOU PROBABLY WANT TO CREATE A SPECIALIZATION OF THIS FUNCTION SPECIFIC 
	      // FOR YOUR EQUATION SYSTEM 
	      
	      //Update solution in primitive variables.
	      Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].W[i][j][k] = 
		Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].U[i][j][k].W(); 
	    } 
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
template <typename SOLN_pSTATE, typename SOLN_cSTATE> 
double Hexa_Newton_Krylov_Schwarz_Solver<SOLN_pSTATE,SOLN_cSTATE>:: 
Finite_Time_Step(const int &Number_of_Newton_Steps){

  double CFL_current;

  //SER 
  if (L2norm_current_n > Solution_Data->Input.NKS_IP.Min_Finite_Time_Step_Norm_Ratio ) { 
    CFL_current = Solution_Data->Input.NKS_IP.Finite_Time_Step_Initial_CFL*
      pow( max(ONE, ONE/L2norm_current_n),ONE ); 
      //      pow(min(ONE, max(ONE, ONE/L2norm_current_n)*MIN_FINITE_TIME_STEP_NORM_RATIO),ONE );     
  } else {
    CFL_current = Solution_Data->Input.NKS_IP.Finite_Time_Step_Initial_CFL;
  } 
 
  return CFL_current;

}


/*! *****************************************************************************************
 *  Routine: Hexa_Newton_Krylov_Schwarz_Solver
 *                                                       
 *  This routine updates the specified solution block using Newton-Krylov-Schwarz method.                  
 *                                                      
 ********************************************************************************************/
template <typename SOLN_pSTATE, typename SOLN_cSTATE> 
int Hexa_Newton_Krylov_Schwarz_Solver<SOLN_pSTATE,SOLN_cSTATE>::
Solve(){

  int error_flag = 0;

  /**************************************************************************/  
  // NKS Specific Output 
  if (CFFC_Primary_MPI_Processor() && !Data->batch_flag){
    //cout << Solution_Data->Input.NKS_IP;   // Causes linking error ?? so call output directly
    Solution_Data->Input.NKS_IP.Output(cout);
    //Input_Parameters.NKS_IP.Memory_Estimates(blocksize,SolnBlk[0].NCi*SolnBlk[0].NCj,Used_blocks_count);
  }  

  /**************************************************************************/  
  /********* Unsteady Time Accurate, Dual Time Stepping  ********************/
  /**************************************************************************/  
  int DTS_Step(1);    
  if (Solution_Data->Input.NKS_IP.Dual_Time_Stepping) { 

    double DTS_dTime(ZERO);
    double physical_time = (Data->Time);   
    int physical_time_param(TIME_STEPPING_IMPLICIT_EULER);

    // Outer Loop (Physical Time)      
    while ( (DTS_Step < Solution_Data->Input.NKS_IP.Maximum_Number_of_DTS_Steps) &&
	    (Solution_Data->Input.Time_Max > physical_time ) ) {

      /**************************************************************************/    
      // First Step needs to be done with Implicit Euler
      if (DTS_Step == 1 ) {
	physical_time_param = Solution_Data->Input.NKS_IP.Physical_Time_Integration;
	Solution_Data->Input.NKS_IP.Physical_Time_Integration = TIME_STEPPING_IMPLICIT_EULER; 
      }
      /**************************************************************************/
      
      /**************************************************************************/
      // Determine global time step using a "fixed" time step or using a CFL                                                       
      if(Solution_Data->Input.NKS_IP.Physical_Time_Step > ZERO){
	DTS_dTime = Solution_Data->Input.NKS_IP.Physical_Time_Step;
      } else {
	//DTS_dTime = CFL(SolnBlk, List_of_Local_Solution_Blocks, Solution_Data->Input);  
	DTS_dTime = Solution_Data->Input.NKS_IP.Physical_Time_CFL_Number*CFFC_Minimum_MPI(DTS_dTime); 
	
	//Last Time sized to get Time_Max
	if( physical_time + DTS_dTime > Solution_Data->Input.Time_Max){
	  DTS_dTime = Solution_Data->Input.Time_Max - physical_time;
	}
      }
      /**************************************************************************/

      /**************************************************************************/
      // Store Previous Solution & Physical Time Step 
      for (int Bcount = 0; Bcount <  Data->Local_Adaptive_Block_List.Nblk; Bcount++) {
	if ( Data->Local_Adaptive_Block_List.Block[Bcount].used == ADAPTIVEBLOCK3D_USED) {
	  //DTS_SolnBlk[Bcount].Store_Previous(SolnBlk[Bcount]);
	  DTS_SolnBlk[Bcount].DTS_dTime = DTS_dTime;
	}
      } 
      /**************************************************************************/

      /**************************************************************************/
      // Solve using NKS 
      error_flag = Steady_Solve(physical_time, DTS_Step);
             
      /**************************************************************************/

      /**************************************************************************/			
      // After first step, reset to requested Time Integration Method
      if (DTS_Step == 1) {
	Solution_Data->Input.NKS_IP.Physical_Time_Integration = physical_time_param;       
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
    //number_of_explicit_time_steps += DTS_Step-1;
    
    /**************************************************************************/  
    /********* Steady State non-time accurate  ********************************/
    /**************************************************************************/  
  } else { 
    error_flag = Steady_Solve(Data->Time, DTS_Step);
  }

  /**************************************************************************/  
  return CFFC_OR_MPI(error_flag); 
}


/*! *****************************************************************************************
 *  Routine: Newton_Krylov_Solver
 *                                                       
 * This routine updates the specified solution block using Newton-Krylov-Schwarz method.                  
 *                                                      
 ********************************************************************************************/
template <typename SOLN_pSTATE, typename SOLN_cSTATE> 
int Hexa_Newton_Krylov_Schwarz_Solver<SOLN_pSTATE,SOLN_cSTATE>::
Steady_Solve(const double &physical_time,const int &DTS_Step){
  
  int error_flag(0);

  /* Local Variabls */ 
  double dTime, CFL_current;
  CPUTime total_cpu_time, start_total_cpu_time; //??
  bool limiter_check= true;

  // SHOULD BE IN GMRES !!!!!! ?????????
  //bool GMRES_Restarted = false, GMRES_Failed = false;
  //int GMRES_Restarts = 0, GMRES_Failures = 0, GMRES_Iters = 0, All_Iters_index = 0;

  int *GMRES_All_Iters = new int[Solution_Data->Input.NKS_IP.Maximum_Number_of_NKS_Iterations+1];
  GMRES_All_Iters[0] = 0;

  //Reset relative Norms ???;
  L2norm_first = L1norm_first = Max_norm_first = ZERO;  
  L2norm_current_n = L1norm_current_n = Max_norm_current_n = ONE;

  /**************************************************************************/
  /******************* BEGIN NEWTON-KRYLOV-SCHWARZ CALCULATION **************/
  /**************************************************************************/
  bool NKS_continue_flag(true);    
  int Number_of_Newton_Steps(1);                // FOR RESTART, THIS SHOULD BE SET TO LAST NKS STEP ??
  int i_limiter = Solution_Data->Input.i_Limiter;

  Data->processor_cpu_time.update();                                         //WORK THESE OUT FOR TOTAL, NKS, & DTS
  start_total_cpu_time.cput = CFFC_Summation_MPI(Data->processor_cpu_time.cput); 

  while ( NKS_continue_flag && Number_of_Newton_Steps <= Solution_Data->Input.NKS_IP.Maximum_Number_of_NKS_Iterations) {

    /**************************************************************************/
    // Limiter Switch to use  first order for first "N" newton steps, then switch to requested method 
    if (Number_of_Newton_Steps <= Solution_Data->Input.NKS_IP.Min_Number_of_Newton_Steps_With_Zero_Limiter) {
      if (CFFC_Primary_MPI_Processor()) {  cout<<"\n Setting Limiter to ZERO, ie. Using First Order"; }
      Solution_Data->Input.i_Limiter = LIMITER_ZERO;   
    } else {
      Solution_Data->Input.i_Limiter = i_limiter;
    } 
    /**************************************************************************/

    // CLOCK //
    clock_t t0 = clock();

    // -R(U_n) or -R(U_n)* for DTS
    /**************************************************************************/
    /* Calculate residual: dUdt[i][j][0]  from U,W */ 
    for ( int Bcount = 0 ; Bcount < Data->Local_Adaptive_Block_List.Nblk; ++Bcount ) {
      if (Data->Local_Adaptive_Block_List.Block[Bcount].used == ADAPTIVEBLOCK3D_USED) {
	error_flag = dUdt_Residual_Evaluation_DTS<SOLN_pSTATE,SOLN_cSTATE>(Solution_Data,
									   DTS_SolnBlk,
									   Bcount);
      } 
    } 
    /**************************************************************************/


    /**************************************************************************/
//     /* Send boundary flux corrections at block interfaces with resolution changes. */
//     error_flag = Send_Conservative_Flux_Corrections(SolnBlk, 
// 		                		    List_of_Local_Solution_Blocks,
// 						    Num_Var);
//     if (error_flag) {
//        cout << "\n NKS ERROR: flux correction message passing error on processor "
//             << List_of_Local_Solution_Blocks.ThisCPU
//             << ".\n";
//        cout.flush();
//     } /* endif */
//     error_flag = CFDkit_OR_MPI(error_flag);
//     if (error_flag) return (error_flag);
	  
//     /* Apply boundary flux corrections to residual to ensure that method is conservative. */
//     Apply_Boundary_Flux_Corrections(SolnBlk, 
// 		                    List_of_Local_Solution_Blocks);
    /**************************************************************************/


    /**************************************************************************/
    /* Calculate 1-, 2-, and max-norms of density residual (dUdt) for all blocks. */   
    Calculate_Norms(Number_of_Newton_Steps);
    /**************************************************************************/
    Data->processor_cpu_time.update();
    Data->total_cpu_time.cput = CFFC_Summation_MPI(Data->processor_cpu_time.cput); 

    // Output Solution Progress 
    if (CFFC_Primary_MPI_Processor() &&
	(!Solution_Data->Input.NKS_IP.Dual_Time_Stepping ||
	 (Solution_Data->Input.NKS_IP.Dual_Time_Stepping && Number_of_Newton_Steps == 1)) ) {  //Output only 1st Step when using DTS

      Output_Progress_to_File(Data->residual_file,
			      Data->number_of_explicit_time_steps+Number_of_Newton_Steps,  //DTS!!!
			      physical_time*THOUSAND,
			      Data->total_cpu_time, 
			      L1norm_current[Solution_Data->Input.p_Norm_Indicator-1],      //maybe switch to current_n so all scale from 1 ???
			      L2norm_current[Solution_Data->Input.p_Norm_Indicator-1],
			      Max_norm_current[Solution_Data->Input.p_Norm_Indicator-1]);  //mod to output "N" norms
			      
    }

    /************************** Restart ***************************************/
    error_flag = Restart_NKS();
    /**************************************************************************/
 
    /**************************************************************************/
    // Freeze Limiters if Residual less than given value or # of orders of reduction.
    if (Solution_Data->Input.Freeze_Limiter && limiter_check){
      if (Number_of_Newton_Steps > 1 &&  L2norm_current[Solution_Data->Input.p_Norm_Indicator] 
	  <= Solution_Data->Input.Freeze_Limiter_Residual_Level)  {    // absolute
	//if (Number_of_Newton_Steps > 1 &&  L2norm_current_n <= Solution_Data->Input.Freeze_Limiter_Residual_Level)  {  // relative 
	if (CFFC_Primary_MPI_Processor()) {
	  cout << "\n\n ********** Apply Limiter Freezing ********** \n";
        } 
	//Freeze_Limiters(SolnBlk, List_of_Local_Solution_Blocks);	
	limiter_check = false;	
      } 
    } 
    /**************************************************************************/

    /**************************************************************************/
    /***************** NEWTON STEP ********************************************/
    /**************************************************************************/
    if (L2norm_current_n > Solution_Data->Input.NKS_IP.Overall_Tolerance){  
     
      /************************** TIME STEP *************************************/ 
      if (Solution_Data->Input.NKS_IP.Finite_Time_Step) {  
	// Apply finite time step, ie. Implicit Euler ( as deltaT -> inf. Newton's Method)
	CFL_current = Finite_Time_Step(Number_of_Newton_Steps);	   
      } else { 
	CFL_current = Solution_Data->Input.NKS_IP.Finite_Time_Step_Final_CFL;   	
      } 

      // Determine time steps
      dTime = Solution_Data->Local_Solution_Blocks.CFL(Solution_Data->Input); 
      
      if (Solution_Data->Input.Local_Time_Stepping == GLOBAL_TIME_STEPPING ){
	dTime = CFFC_Minimum_MPI(dTime); 
	dTime *= CFL_current;           
	Solution_Data->Local_Solution_Blocks.Set_Global_TimeStep(dTime);      
      } else { 
	for ( int Bcount = 0 ; Bcount < Data->Local_Adaptive_Block_List.Nblk; ++Bcount ) {
	  if (Data->Local_Adaptive_Block_List.Block[Bcount].used == ADAPTIVEBLOCK3D_USED) {
	    for (int k = Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].KCl; 
		 k <= Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].KCu; k++){
	      for (int j = Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].JCl; 
		   j <= Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].JCu; j++){
		for (int i = Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].ICl; 
		     i <= Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].ICu; i++){
		  Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].dt[i][j][k] *= CFL_current;
		}
	      }
	    }
	  }
	}
      }
      /**************************************************************************/   

      /**************************************************************************/
      // Print out Newtwon Step info at the beginning of the iteration. 
      cout.precision(10);
      if (CFFC_Primary_MPI_Processor()) {      	
	cout << "\n Newton Step (Outer It.) = " << Number_of_Newton_Steps << " L2norm = "
	     << L2norm_current[Solution_Data->Input.p_Norm_Indicator-1] << " L2norm_ratio = " << L2norm_current_n 
	     << " CFL = " << CFL_current <<" min_deltat = "<<CFL_current*dTime;
      } 
      /**************************************************************************/

   
      // POSSIBLE MOVE THIS TO ITS OWN FUNCTION, LIKE NEWTON_UPDATE
      /**************************************************************************/
      // Store  Uo, dt for use in GMRES & Preconditioner Calculations
      for ( int Bcount = 0 ; Bcount < Data->Local_Adaptive_Block_List.Nblk; ++Bcount ) {
	if (Data->Local_Adaptive_Block_List.Block[Bcount].used == ADAPTIVEBLOCK3D_USED) {
	  /* Copy solutions in conserved variables, U to Uo for all blocks for update procedure. */
	  for (int k = Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].KCl 
		 - Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].Nghost;
	       k <= Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].KCu 
		 + Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].Nghost; k++){
	    for (int j = Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].JCl 
		   - Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].Nghost; 
		 j <= Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].JCu 
		   + Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].Nghost; j++){
	      for (int i = Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].ICl 
		     - Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].Nghost; 
		   i <= Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].ICu 
		     + Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].Nghost; i++){	
		Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].Uo[i][j][k]= 
		  Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].U[i][j][k];
	          	      		
		//Solution_Data->Local_Solution_Blocks.Soln_Blks[Bcount].dt[i][j][k] *= CFL_current;  //instead of above???
	      }
	    } 
	  } 	  
	} 
      }  
      /**************************************************************************/
  

      //DTS RULES: update Preconditioner only if GMRES iterations increase between Newton Its.
      bool GMRES_Iters_increaseing(true);
//       if( GMRES_All_Iters[Number_of_Newton_Steps-1] == 1){
// 	GMRES_Iters_increaseing = false;
//       } else if (GMRES_All_Iters[Number_of_Newton_Steps-2] > GMRES_All_Iters[Number_of_Newton_Steps-1]) {
// 	GMRES_Iters_increaseing = true;
//       } else {
// 	GMRES_Iters_increaseing = false;
//       }

      /**************************************************************************/
      /************* PRECONDTIONER "BLOCK" JACOBIANS ****************************/      
      /**************************************************************************/
      // Create/Update Jacobian Matrix(s) using Uo = U  
      if ( ( !Solution_Data->Input.NKS_IP.Dual_Time_Stepping &&
	     (Number_of_Newton_Steps < Solution_Data->Input.NKS_IP.Min_Number_of_Newton_Steps_Requiring_Jacobian_Update || 
	      L2norm_current_n > Solution_Data->Input.NKS_IP.Min_L2_Norm_Requiring_Jacobian_Update) ) || 
	   ( Solution_Data->Input.NKS_IP.Dual_Time_Stepping && GMRES_Iters_increaseing) ) {                       
	
	cout << "\n Creating/Updating Jacobian Matrix"; 
  
	//CLOCK
	clock_t t0 = clock();

	//Update for each block.
	for ( int Bcount = 0 ; Bcount < Data->Local_Adaptive_Block_List.Nblk; ++Bcount ) {
	  if (Solution_Data->Local_Solution_Blocks.Block_Used[Bcount] == ADAPTIVEBLOCK3D_USED){
	    // Jacobian changed and Preconditioner updated 
	    Block_precon[Bcount].Update_Jacobian_and_Preconditioner(DTS_SolnBlk[Bcount].DTS_dTime);
	  } 
	} 

// 	preconditioner_cputime += double(clock() - t0) / double(CLOCKS_PER_SEC);

      }
      /**************************************************************************/


      /**************************************************************************/
      /************* LINEAR SYSTEM SOLVE WITH GMRES  ****************************/      
      /**************************************************************************/
      /* Solve system with right-preconditioned matrix free GMRES */  
      error_flag = GMRES.solve(Block_precon);
      if (CFFC_Primary_MPI_Processor() && error_flag) cout << "\n NKS2D: Error in GMRES \n";

      GMRES_All_Iters[Number_of_Newton_Steps-1] = GMRES.GMRES_Iterations();
      /**************************************************************************/      
 

      /**************************************************************************/      
      // Solution Update U = Uo + GMRES.deltaU
      error_flag = Newton_Update();      
      if (CFFC_Primary_MPI_Processor() && error_flag) cout <<  "\n NKS2D: Error in Solution Update \n";
      /**************************************************************************/
      

//       /**************************************************************************/
//       /* Exchange solution information between neighbouring blocks. */    
//       error_flag = Send_All_Messages(SolnBlk,
// 				     List_of_Local_Solution_Blocks,
// 				     Num_Var, 
// 				     OFF);
//       if (error_flag) {
//          cout << "\n 2D_NKS ERROR: 2D message passing error on processor "
// 	      << List_of_Local_Solution_Blocks.ThisCPU
// 	      << ".\n";
//          cout.flush();
//       } /* endif */
//       error_flag = CFFC_OR_MPI(error_flag);
//       if (error_flag) return (error_flag);
//       /**************************************************************************/

      /**************************************************************************/  
      /* Apply boundary conditions for Newton step. */
      Solution_Data->Local_Solution_Blocks.BCs(Solution_Data->Input);
      /**************************************************************************/

      /**************************************************************************/
      Number_of_Newton_Steps++;      

      // END OF  if (L2norm_current_n > Solution_Data->Input.NKS_IP.Overall_Tolerance)
      /**************************************************************************/
    } else {       
      // L2norm_current_n < Solution_Data->Input.NKS_IP.Overall_Tolerance so set flag to "stop"     
      NKS_continue_flag = false;
    } /* endif */
    /**************************************************************************/
 
  }  // END OF NEWTON ITERATION WHILE LOOP
 
  /**************************************************************************/  
  /********* FINISHED NEWTON KRYLOV SCHWARZ *********************************/
  /**************************************************************************/  
  
  /**************************************************************************/
  /* Reset limiter. */
  Solution_Data->Input.i_Limiter = i_limiter;
  /**************************************************************************/

  /**************************************************************************/
  /* Output final 2-norm for all blocks */ 
  if (CFFC_Primary_MPI_Processor()) {     
    cout << " " << endl;
    for (int star=0;star<75;star++){cout <<"*";}
    cout << "\nEnd of Newton Steps = " << Number_of_Newton_Steps-1  << " L2norm = "
	 << L2norm_current[Solution_Data->Input.p_Norm_Indicator-1] << " L2norm_ratio = " << L2norm_current_n << endl;
    for (int star=0;star<75;star++){cout <<"*";}
  } /* endif */
  /**************************************************************************/

  
//   //For total its count used in Restart ??
//   /**************************************************************************/
//   number_of_explicit_time_steps = number_of_explicit_time_steps + Number_of_Newton_Steps-1;
//   /**************************************************************************/    


  delete[] GMRES_All_Iters;


  return error_flag;
} // End of Steady Newton_Krylov_Schwarz_Solver.




#endif  /* _NKS_INCLUDED */

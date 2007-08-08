#ifndef _RTE2D_IMPLICIT_SPECIALIZATION_INCLUDED 
#define _RTE2D_IMPLICIT_SPECIALIZATION_INCLUDED 

#ifndef _RTE2D_QUAD_INCLUDED
#include "Rte2DQuad.h"
#endif // _RTE2D_QUAD_INCLUDED

#ifndef _NKS_INCLUDED
#include "../NewtonKrylovSchwarz2D/NKS2D.h"
#endif // _NKS_INCLUDED


// Prototypes required for instantiation

/*  NKS FUNCTION PROTOTYPES */
template <>
int Newton_Krylov_Schwarz_Solver<Rte2D_State,
				 Rte2D_Quad_Block,
				 Rte2D_Input_Parameters>(CPUTime &NKS_processor_cpu_time,
							 CPUTime &NKS_total_cpu_time,
							 ostream &residual_file,   
							 int &number_of_explicit_time_steps,
							 Rte2D_Quad_Block *Soln_ptr,
							 AdaptiveBlock2D_List &List_of_Local_Solution_Blocks,
							 Rte2D_Input_Parameters &Input_Parameters);

template <>
int Newton_Update<Rte2D_State,
		  Rte2D_Quad_Block,
		  Rte2D_Input_Parameters>(Rte2D_Quad_Block *Soln_ptr,
					  AdaptiveBlock2D_List &List_of_Local_Solution_Blocks,
					  Rte2D_Input_Parameters &Input_Parameters,
					  GMRES_RightPrecon_MatrixFree<Rte2D_State,Rte2D_Quad_Block,Rte2D_Input_Parameters> &GMRES);


/*  GMRES FUNCTION PROTOTYPES */

template <> 
void GMRES_Block<Rte2D_State,Rte2D_Quad_Block,Rte2D_Input_Parameters>::
  Output_U(int what);

template<> 
void GMRES_Block<Rte2D_State,Rte2D_Quad_Block,Rte2D_Input_Parameters>::
  set_normalize_values(void);

template<> 
void GMRES_Block<Rte2D_State,Rte2D_Quad_Block,Rte2D_Input_Parameters>::
  GMRES_BCs(double *v, const double epsilon);

template <> 
void GMRES_Block<Rte2D_State,Rte2D_Quad_Block,Rte2D_Input_Parameters>::
  calculate_perturbed_residual(const double &epsilon);

template <> 
void GMRES_Block<Rte2D_State,Rte2D_Quad_Block,Rte2D_Input_Parameters>::
  calculate_perturbed_residual_Restart(const double &epsilon);

template <> 
int GMRES_Block<Rte2D_State,Rte2D_Quad_Block,Rte2D_Input_Parameters>::
  LoadSendBuffer_F2C(double *buffer,
		     int &buffer_count,
		     const int buffer_size,
		     const int i_min, 
		     const int i_max,
		     const int i_inc,
		     const int j_min, 
		     const int j_max,
		     const int j_inc);

template <> 
int GMRES_Block<Rte2D_State,Rte2D_Quad_Block,Rte2D_Input_Parameters>::
  LoadSendBuffer_C2F(double *buffer,
		     int &buffer_count,
		     const int buffer_size,
		     const int i_min, 
		     const int i_max,
		     const int i_inc,
		     const int j_min, 
		     const int j_max,
		     const int j_inc);

template <> 
void GMRES_Block<Rte2D_State,Rte2D_Quad_Block,Rte2D_Input_Parameters>::
  SubcellReconstruction(const int i, 
			const int j,
			const int Limiter);


/*  BLOCK PRECONDITIONER FUNCTION PROTOTYPES */

template <>
void Block_Preconditioner<Rte2D_State,Rte2D_Quad_Block,Rte2D_Input_Parameters>::
  Setup_Jacobian_approximation();

template <>
void Block_Preconditioner<Rte2D_State,Rte2D_Quad_Block,Rte2D_Input_Parameters>::
  Update_Jacobian_and_Preconditioner();

template <>
void Block_Preconditioner<Rte2D_State,Rte2D_Quad_Block,Rte2D_Input_Parameters>::
  Implicit_Euler(const int &cell_index_i,const int &cell_index_j, DenseMatrix* Jacobian);

template<> 
void Block_Preconditioner<Rte2D_State,Rte2D_Quad_Block,Rte2D_Input_Parameters>::
  First_Order_Inviscid_Jacobian_HLLE(const int &cell_index_i,const int &cell_index_j, 
				     DenseMatrix* Jacobian);            

template<> 
void Block_Preconditioner<Rte2D_State,Rte2D_Quad_Block,Rte2D_Input_Parameters>::
  First_Order_Inviscid_Jacobian_Roe(const int &cell_index_i,const int &cell_index_j, 
				    DenseMatrix* Jacobian);          



/****************************************************************/
/****** NKS REQUIRED SPECIALIZATIONS & FUNCTIONS ****************/
/****************************************************************/  



/********************************************************************************************
//!  Routine: Newton_Krylov_Solver
 *                                                       
 * This routine updates the specified solution block    
 * using Newton-Krylov-Schwarz method.                  
 *                                                      
 ********************************************************************************************/
template <>
int Newton_Krylov_Schwarz_Solver<Rte2D_State,
				 Rte2D_Quad_Block,
				 Rte2D_Input_Parameters>(CPUTime &NKS_processor_cpu_time,
							 CPUTime &NKS_total_cpu_time,
							 ostream &residual_file,   
							 int &number_of_explicit_time_steps,
							 Rte2D_Quad_Block *Soln_ptr,
							 AdaptiveBlock2D_List &List_of_Local_Solution_Blocks,
							 Rte2D_Input_Parameters &Input_Parameters) { 
  
  double L2norm_current   = ZERO, 
         L1norm_current   = ZERO, 
         Max_norm_current = ZERO; 
  // IN case of restart need to reset this to last value for relative tolerancing to work...
  double L2norm_first     = ZERO, 
         L1norm_first     = ZERO,
         Max_norm_first   = ZERO;
  double L2norm_current_n = ONE,  //ratio of first to current norm
         L1norm_current_n = ONE,
         Max_norm_current_n = ONE;
  //double norm_ratio       = ZERO;
  double dTime, CFL_current, Mrefnew(STARTING_MREF);

  int  error_flag = 0;
  bool limiter_check= true;
  int Used_blocks_count = 0; 

  int Num_Var = Soln_ptr[0].NumVar();  //Number of equations in "p/c" State
  int blocksize = set_blocksize<Rte2D_Quad_Block>(Soln_ptr[0]); //Number of equations used for calculation

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
  GMRES_RightPrecon_MatrixFree<Rte2D_State,Rte2D_Quad_Block,Rte2D_Input_Parameters> GMRES_(
							    Soln_ptr,
							    List_of_Local_Solution_Blocks,
							    Input_Parameters,
							    blocksize);

  //SETUP Preconditioner Blocks
  Block_Preconditioner<Rte2D_State,Rte2D_Quad_Block,Rte2D_Input_Parameters> *Block_precon 
    = new Block_Preconditioner<Rte2D_State,Rte2D_Quad_Block,Rte2D_Input_Parameters>[Used_blocks_count];

  for( int i=0; i<Used_blocks_count; i++) {  
    Block_precon[i].Create_Preconditioner(Soln_ptr[i], Input_Parameters,blocksize);
  }  
  //total_Used_blocks = CFDkit_Summation_MPI(Used_blocks_count);  maybe ??      

  /**************************************************************************/  
  /********* FINISHED SETTING UP STORAGE AND DATASTRUCTURES *****************/
  /**************************************************************************/  
  if (CFDkit_Primary_MPI_Processor()){
    Input_Parameters.NKS_IP.Output(); 
    //Input_Parameters.NKS_IP.Memory_Estimates(blocksize,Soln_ptr[0].NCi*Soln_ptr[0].NCj,Used_blocks_count);
  }
   
  /**************************************************************************/
  /******************* BEGIN NEWTON-KRYLOV-SCHWARZ CALCULATION **************/
  /**************************************************************************/
  bool NKS_continue_flag = true;    
  bool Update_Jacobian_flag = true;  //what's the deal with this flag, I think it is basically useless....

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
//     if(Input_Parameters.Preconditioning && Mrefnew > Input_Parameters.Mach_Number_Reference ){      
//       if(L2norm_current_n < 0.1){
// 	//start at 1 and decrease to Mref set in input
// 	Mrefnew = STARTING_MREF - (STARTING_MREF-Input_Parameters.Mach_Number_Reference)*(log10(ONE/L2norm_current_n)/TWO);     
// 	if( Mrefnew < Input_Parameters.Mach_Number_Reference) Mrefnew = Input_Parameters.Mach_Number_Reference;
// 	}
//       Change_Mref(Soln_ptr, List_of_Local_Solution_Blocks,Mrefnew);
//     }	
//     /**************************************************************************/

    // -R(U_n)
    /**************************************************************************/
    /* Calculate residual : dudt[i][j][0]  from U,W */
    for ( int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks.Nblk ; Bcount++ ) {
      if (List_of_Local_Solution_Blocks.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	dUdt_Residual_Evaluation(Soln_ptr[Bcount],Input_Parameters);
      } 
    } 
    /**************************************************************************/

    /**************************************************************************/
    /* Send boundary flux corrections at block interfaces with resolution changes. */
    error_flag = Send_Conservative_Flux_Corrections(Soln_ptr, 
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
    Apply_Boundary_Flux_Corrections(Soln_ptr, 
		                    List_of_Local_Solution_Blocks);
    /**************************************************************************/

    /**************************************************************************/
    /* Calculate 1-, 2-, and max-norms of density residual (dUdt) for all blocks. */   
    L2norm_current   = sqr(L2_Norm_Residual(Soln_ptr, List_of_Local_Solution_Blocks));
    L2norm_current   = sqrt(CFDkit_Summation_MPI(L2norm_current));
    L1norm_current   = L1_Norm_Residual(Soln_ptr, List_of_Local_Solution_Blocks);
    L1norm_current   = CFDkit_Summation_MPI(L1norm_current);
    Max_norm_current = Max_Norm_Residual(Soln_ptr, List_of_Local_Solution_Blocks);
    Max_norm_current = CFDkit_Summation_MPI(Max_norm_current);

    if (Number_of_Newton_Steps == 1 ) {
      L2norm_first = max(L2norm_first,L2norm_current);  //another restart cludge
      L1norm_first = L1norm_current;
      Max_norm_first = Max_norm_current;
    } else {
      L2norm_first = max(L2norm_first, L2norm_current);
      L1norm_first = max(L1norm_first, L1norm_current);
      Max_norm_first = max(Max_norm_first, Max_norm_current);   
    } 
    
    L2norm_current_n   = L2norm_current / L2norm_first; 
    L1norm_current_n   = L1norm_current / L1norm_first; 
    Max_norm_current_n = Max_norm_current / Max_norm_first;

    /* Calculate ratio of initial and current 2-norms. */
    //norm_ratio = L2norm_first / L2norm_current;
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
//       error_flag = Write_Restart_Solution(Soln_ptr, 
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
//       error_flag = CFDkit_OR_MPI(error_flag);4
//       if (error_flag) return (error_flag);  cout.flush();
//     } 

    // Output progress information for the calculation
    if (CFDkit_Primary_MPI_Processor()){
      Output_Progress_to_File(residual_file,
			      number_of_explicit_time_steps+Number_of_Newton_Steps-1,
			      ZERO,
			      NKS_total_cpu_time,     //Should add explicit processor time as well????
			      L1norm_current,
			      L2norm_current,
			      Max_norm_current);
    }
    /**************************************************************************/
 
    /**************************************************************************/
    // Freeze Limiters if Residual less than given value or # of orders of reduction.
    if (Input_Parameters.Freeze_Limiter && limiter_check){
//       if (Number_of_Newton_Steps > 1 &&  L2norm_current <= Input_Parameters.Freeze_Limiter_Residual_Level)  {
      if (Number_of_Newton_Steps > 1)  {
	if (CFDkit_Primary_MPI_Processor()) {
	  cout << "\n\n ********** Apply Limiter Freezing ********** \n";
        } 
	Freeze_Limiters(Soln_ptr, List_of_Local_Solution_Blocks);	
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
      dTime = CFL(Soln_ptr, List_of_Local_Solution_Blocks, Input_Parameters);
      dTime = CFDkit_Minimum_MPI(dTime); 

      //NOTE: SHOULD PULL THIS OUT AS A TEMPLATE FUNCTION 
      // SO IT CAN BE SPECIALIZED IF REQUIRED .

      // Apply finite time step, ie. Implicit Euler ( as deltaT -> inf. Newtons Method)
      if (Input_Parameters.NKS_IP.Finite_Time_Step) {  

	//SER with curious scaling
        if (L2norm_current_n > MIN_FINITE_TIME_STEP_NORM_RATIO ) { 
	  //Original works for Invisicd, Viscous OK
	  CFL_current = Input_Parameters.NKS_IP.Finite_Time_Step_Initial_CFL*
	    pow(min(ONE, max(ONE, ONE/L2norm_current_n)*MIN_FINITE_TIME_STEP_NORM_RATIO),ONE ); 

// 	  if( L2norm_current_n > 0.1){
// 	    CFL_current = Input_Parameters.NKS_IP.Finite_Time_Step_Initial_CFL*
// 	      pow(min(ONE, max(ONE, ONE/L2norm_current_n)*MIN_FINITE_TIME_STEP_NORM_RATIO),ONE ); 
// 	  } else if( L2norm_current_n > 0.001){
// 	    CFL_current = Input_Parameters.NKS_IP.Finite_Time_Step_Initial_CFL*
// 	      pow(min(ONE, max(ONE, ONE/L2norm_current_n)*MIN_FINITE_TIME_STEP_NORM_RATIO),0.95); //0.92
// 	  } else if( L2norm_current_n > 0.0001){
// 	    CFL_current = Input_Parameters.NKS_IP.Finite_Time_Step_Initial_CFL*
// 	      pow(min(ONE, max(ONE, ONE/L2norm_current_n)*MIN_FINITE_TIME_STEP_NORM_RATIO),0.8); 
// 	  } else {
// 	    CFL_current = Input_Parameters.NKS_IP.Finite_Time_Step_Initial_CFL*
// 	      pow(min(ONE, max(ONE, ONE/L2norm_current_n)*MIN_FINITE_TIME_STEP_NORM_RATIO),0.7); 
// 	  }

        } else {
	  CFL_current = Input_Parameters.NKS_IP.Finite_Time_Step_Initial_CFL;
        } 

      } else { 
        CFL_current = Input_Parameters.NKS_IP.Finite_Time_Step_Initial_CFL; 
      }      
      dTime = CFL_current*dTime;
     

      //Set all dt[i][j] for global time stepping if requested (shouldn't be as this is NOT Time accurate)
      if(Input_Parameters.Local_Time_Stepping == GLOBAL_TIME_STEPPING ){
	Set_Global_TimeStep(Soln_ptr, List_of_Local_Solution_Blocks, dTime);      
      }
      /**************************************************************************/
  

      /**************************************************************************/
      // Print out Newtwon Step info at the beginning of the iteration. 
      cout.precision(10);
      if (CFDkit_Primary_MPI_Processor()) {      	
	cout << "\n Newton Step (Outer It.) = " << Number_of_Newton_Steps << " L2norm = "
	     << L2norm_current << " L2norm_ratio = " << L2norm_current_n << "0 CFL = " << CFL_current <<" min_deltat = "<<dTime;
	// Commented out for Rte2D
// 	if(Input_Parameters.Preconditioning) cout<<" Mref = "<<Mrefnew;
      } 
      /**************************************************************************/

   
      /**************************************************************************/
      // Store  Uo, dt for use in GMRES & Precondtioner Calculations
      for ( int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks.Nblk ; ++Bcount ) {
	if (List_of_Local_Solution_Blocks.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {	  
	  /* Copy solutions in conserved variables, U to Uo for all blocks for update procedure. */
	  for (int j = Soln_ptr[Bcount].JCl-Soln_ptr[Bcount].Nghost;  j <= Soln_ptr[Bcount].JCu+Soln_ptr[Bcount].Nghost; j++){
	    for (int i = Soln_ptr[Bcount].ICl-Soln_ptr[Bcount].Nghost;  i <= Soln_ptr[Bcount].ICu+Soln_ptr[Bcount].Nghost; i++){
// 	      for(int varindex =1; varindex <= Num_Var; varindex++){		
// 		Soln_ptr[Bcount].Uo[i][j][varindex] = Soln_ptr[Bcount].U[i][j][varindex];
// 	      }    	   
	      Soln_ptr[Bcount].Uo[i][j] = Soln_ptr[Bcount].U[i][j];
	      /* Set time step. */
              Soln_ptr[Bcount].dt[i][j] = CFL_current*Soln_ptr[Bcount].dt[i][j];
	    } 
	  } 	  
	} 
      }  
      /**************************************************************************/


      /**************************************************************************/
      /************* PRECONDTIONER "BLOCK" JACOBIANS ****************************/      
      /**************************************************************************/
      // Create/Update Jacobian Matrix(s) using Uo = U  
      if ((Number_of_Newton_Steps < MIN_NUMBER_OF_NEWTON_STEPS_REQUIRING_JACOBIAN_UPDATE
	   || L2norm_current_n > MIN_L2_NORM_REQUIRING_JACOBIAN_UPDATE ) && Update_Jacobian_flag ) { 
	
	if (CFDkit_Primary_MPI_Processor()) cout << "\n Creating/Updating Jacobian Matrix";  
	
	//Update for each block.
	for ( int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks.Nblk ; ++Bcount ) {
	  if (List_of_Local_Solution_Blocks.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	    // Jacobian changed and Preconditioner updated 
	    Block_precon[Bcount].Update_Jacobian_and_Preconditioner();
	  } 
	} 
      } else { 
	// don't update the Jacobian any more, but this is a stupid use of a flag 
	// as the Number of Iterations of tolerance must be violated first anyway. ????
	Update_Jacobian_flag == false;
      } /* endif */ 
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
      error_flag = Newton_Update<Rte2D_State,Rte2D_Quad_Block,Rte2D_Input_Parameters>
	(Soln_ptr,List_of_Local_Solution_Blocks,Input_Parameters,GMRES_);      
      if (CFDkit_Primary_MPI_Processor() && error_flag) cout <<  "\n NKS2D: Error in Solution Update \n";
      /**************************************************************************/


      /**************************************************************************/
      /* Exchange solution information between neighbouring blocks. */    
      error_flag = Send_All_Messages(Soln_ptr,
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
      BCs(Soln_ptr,List_of_Local_Solution_Blocks);
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
	 << L2norm_current << " L2norm_ratio = " << L2norm_current_n << endl;
    for (int star=0;star<75;star++){cout <<"*";}
  } /* endif */
  /**************************************************************************/

  /**************************************************************************/
  number_of_explicit_time_steps = number_of_explicit_time_steps + Number_of_Newton_Steps-1;
  /**************************************************************************/
    
  /**************************************************************************/
  // Housekeeping 
  delete[] Block_precon; 
  /**************************************************************************/

  return (error_flag);

} /* End of Newton_Krylov_Schwarz_Solver. */




/*! *****************************************************************************************
 *  Specialization of Newton_Update Function                                                *
 *                                                                                          *
 * This routine updates the previous Solution Data Uo with the deltaU from the GMRES        *                  
 * iterative solver.  U = Uo + GMRES.delatU                                                 *
 *                                                                                          *
 ********************************************************************************************/
template <>
int Newton_Update<Rte2D_State,
		  Rte2D_Quad_Block,
		  Rte2D_Input_Parameters>(Rte2D_Quad_Block *Soln_ptr,
					  AdaptiveBlock2D_List &List_of_Local_Solution_Blocks,
					  Rte2D_Input_Parameters &Input_Parameters,
					  GMRES_RightPrecon_MatrixFree<Rte2D_State,Rte2D_Quad_Block,Rte2D_Input_Parameters> &GMRES) 
{

  int Num_Var = Soln_ptr[0].NumVar();  
  int NegValue;
  
  /* Update Solution. No updates to Ghost Cells, let the BC's take care of it */
  for ( int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks.Nblk ; ++Bcount ) {
    if (List_of_Local_Solution_Blocks.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
      for (int j = Soln_ptr[Bcount].JCl; j <= Soln_ptr[Bcount].JCu; j++){
	for (int i = Soln_ptr[Bcount].ICl; i <= Soln_ptr[Bcount].ICu; i++){
	  
	  /* Update solutions in conversed variables  U = Uo + deltaU = Uo + denormalized(x) */	 
	  for(int varindex =1; varindex <= Num_Var; varindex++){	
	    Soln_ptr[Bcount].U[i][j][varindex] = Soln_ptr[Bcount].Uo[i][j][varindex] 
	      +  GMRES.deltaU(Bcount,i,j,varindex-1);
	      } 	      
	  
// 	  /**************************************************************************/
// 	  /* Apply update reduction while any one of the updated variables is negative. */ 
// 	  NegValue = Soln_ptr[Bcount].U[i][j].NegIntensity();
// 	  if (NegValue) {    //THIS SEEMS TO CAUSE MORE PROBLEMS 
// 	                     //THAN HELP, ANY IDEAS, MAYBE MPI ACROSS ALL ????
// 	    double update_reduction_factor = ONE;
	    
// 	    for (int n_update_reduction = 1; n_update_reduction <= 10; ++n_update_reduction) {		  
// 	      update_reduction_factor = HALF*update_reduction_factor;		  		  
// 	      for(int varindex = 1; varindex <= Num_Var; varindex++){		
// 		Soln_ptr[Bcount].U[i][j][varindex] = Soln_ptr[Bcount].Uo[i][j][varindex] 
// 		  + GMRES.deltaU(Bcount,i,j,varindex-1)*update_reduction_factor;
// 	      }   
	      
// 	      cout<<"\n Applying Reduction to solution in NKS "<<n_update_reduction;
// 	      NegValue = Soln_ptr[Bcount].U[i][j].NegIntensity();
// 	      if ( !NegValue ) break; 
// 	    } 
// 	  } 
	  
// 	  /**************************************************************************/
// 	  /* Print error: Negative Density or Negative Pressure. */ 	  
// 	  if (NegValue) { 
// 	    cout << "\n NEGATIVE INTENSITY : " << endl;
// 	    cout << "Soln_ptr["<<Bcount<<"].U["<<i<<"]["<<j<<"].I = " 
// 		 << Soln_ptr[Bcount].U[i][j]
// 		 << ": Soln_ptr["<<Bcount<<"].Uo["<<i<<"]["<<j<<"].d = " 
// 		 << Soln_ptr[Bcount].Uo[i][j] << endl;
// 	    cout << "   G["<<Bcount<<"].x[G["<<Bcount<<"].index("<<i
// 		 <<","<<j
// 		 <<")] = " 
// 		 << GMRES.deltaU(Bcount,i,j,0) << endl;
// 	  } /* endif */ 
	  
// 	  /**************************************************************************/
	  //Update solution in primitive variables.
// 	  Soln_ptr[Bcount].W[i][j] = W(Soln_ptr[Bcount].U[i][j]);	  
	} 
      } 
    } 
  }   
  return 0; 
}





/************************************************************************/
/************ GMRES REQUIRED SPECIALIZATIONS & FUNCTIONS ****************/
/************************************************************************/   



/**************************************************************************
 * Routine: TESTING FUNCITON                                              *
 **************************************************************************/
template <> 
inline void GMRES_Block<Rte2D_State,
                        Rte2D_Quad_Block,
                        Rte2D_Input_Parameters>::
Output_U(int what)
{    
  for (int j = JCl - Nghost ; j <= JCu + Nghost ; j++) {  //includes ghost cells 
    for (int i = ICl - Nghost ; i <= ICu + Nghost ; i++) {
      if(what == 1)cout<<"\n U ("<<i<<","<<j<<") = "<<SolBlk_ptr->U[i][j];      
      if(what == 3)cout<<"\n dUdt ("<<i<<","<<j<<") = "<<SolBlk_ptr->dUdt[i][j][0];
    }
  }  
}



/*!********************************************************
 * GMRES_Block::set_normalize_values for NKS/GMRES        * 
 *              *                                         *
 *   normalize_values[0] must be set to ao                *
 *   normalize_values[1-n] = values for index[1-n]        *
 *             where n = the number of solution variables *
 **********************************************************/
template<> inline void GMRES_Block<Rte2D_State,
				   Rte2D_Quad_Block,
				   Rte2D_Input_Parameters>::
set_normalize_values(void)
{   

  double temp = ZERO;
  for (int j = JCl - Nghost ; j <= JCu + Nghost ; j++) {  
    for (int i = ICl - Nghost ; i <= ICu + Nghost ; i++) {
      temp += fabs(Ib(Input_Parameters->Temperature));
    }
  }  
  
  // NOT NORMALIZED
  for(int i=0; i < blocksize; i++) {
    normalize_valuesU[i] = temp/scalar_dim;
    normalize_valuesR[i] = temp/scalar_dim;
  }
}



/********************************************************
 * Routine: GMRES_BCs                                   *
 *   Specialization of GMRES_Block::GMRES_BCS for       *       
 *   Chem2D_Quad_Block                                  *
 *                                                      *
 *  NOTE: Remeber "v" is the z or x vector representing *
 *        deltaU ie.                                    *
 ********************************************************/
template<> void GMRES_Block<Rte2D_State,
                            Rte2D_Quad_Block,
                            Rte2D_Input_Parameters>::
GMRES_BCs(double *v, const double epsilon)
{ 
  // DO NOTHING
}


/**************************************************************************
 * Routine: GMRES_Block::calculate_pertubed_residual                      *
 **************************************************************************/
// Calculate Soln_ptr.U =  Soln_ptr.Uo + denormalize( epsilon * W(i) )
template <> 
inline void GMRES_Block<Rte2D_State,Rte2D_Quad_Block,Rte2D_Input_Parameters>::
calculate_perturbed_residual(const double &epsilon)
{    
  for (int j = JCl - Nghost ; j <= JCu + Nghost ; j++) {  //includes ghost cells 
    for (int i = ICl - Nghost ; i <= ICu + Nghost ; i++) {
      for(int varindex = 0; varindex < blocksize; varindex++){	
	SolBlk_ptr->U[i][j][varindex+1] = SolBlk_ptr->Uo[i][j][varindex+1] + 
	  denormalizeU( epsilon*W[search_directions*scalar_dim+index(i,j,varindex)], varindex);
      }   
      
      //Chem2D spec_check  to make sure  U0 + epsilon*W(i) > ZERO  
//       if(SolBlk_ptr->U[i][j].NegIntensity()){
// 	cerr << SolBlk_ptr->U[i][j].NegIntensity() << endl;
// 	cerr << "Rte2DQuad_Implicit_Specializations - calculate_perturbed_residual_Restart:"
// 	     << " Negative value encountered.\n"
// 	     << " cell = (" << i << ", " << j << ") " 
// 	     << " X = " << SolBlk_ptr->Grid.Cell[i][j].Xc << "\n U = " 
// 	     << SolBlk_ptr->U[i][j]
// 	     << endl;
// 	exit(1);
//       }

      /* Update primitive variables. */
//       SolBlk_ptr->W[i][j] = SolBlk_ptr->U[i][j].W();      
    }
  }  
}

// Calculate Soln_ptr.U =  Soln_ptr.Uo + denormalize( epsilon * x(i) )
template <> 
inline void GMRES_Block<Rte2D_State,Rte2D_Quad_Block,Rte2D_Input_Parameters>::
calculate_perturbed_residual_Restart(const double &epsilon)
{    
  for (int j = JCl - Nghost ; j <= JCu + Nghost ; j++) {
    for (int i = ICl - Nghost ; i <= ICu + Nghost ; i++) {
      for(int varindex = 0; varindex < blocksize; varindex++){	
	SolBlk_ptr->U[i][j][varindex+1] = SolBlk_ptr->Uo[i][j][varindex+1] + 
	  denormalizeU( epsilon*x[index(i,j,varindex)], varindex);
      }  

      //Rte2D spec_check or maybe a check on epsilon *W(i) 
//       if(SolBlk_ptr->U[i][j].NegIntensity()){
// 	cerr << SolBlk_ptr->U[i][j].NegIntensity() << endl;
// 	cerr << "Rte2DQuad_Implicit_Specializations - calculate_perturbed_residual_Restart:"
// 	     << " Negative value encountered.\n"
// 	     << " cell = (" << i << ", " << j << ") " 
// 	     << " X = " << SolBlk_ptr->Grid.Cell[i][j].Xc << "\n U = " 
// 	     << SolBlk_ptr->U[i][j]  
// 	     << endl;
// 	exit(1);
//       }
 
      /* Update primitive variables. */
//       SolBlk_ptr->W[i][j] = SolBlk_ptr->U[i][j].W();      
    }
  }  
}

/*******************************************************************************
 * GMRES_Block::LoadSendBuffer_F2C -- Loads send message buffer for            *
 *                                    fine to coarse block message             *
 *                                    passing.                                 *
 *******************************************************************************/
template <> 
inline int GMRES_Block<Rte2D_State,Rte2D_Quad_Block,Rte2D_Input_Parameters>::
LoadSendBuffer_F2C(double *buffer,
		   int &buffer_count,
		   const int buffer_size,
		   const int i_min, 
		   const int i_max,
		   const int i_inc,
		   const int j_min, 
		   const int j_max,
		   const int j_inc) {
  int i, j, k;
  for ( j  = j_min ; ((j_inc+2)/4) ? (j < j_max):(j > j_max) ; j += j_inc ) {
     for ( i = i_min ;  ((i_inc+2)/4) ? (i < i_max):(i > i_max) ; i += i_inc ) {
        for ( k = 1 ; k <= blocksize; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
	   if (vector_switch) {
              buffer[buffer_count] = (SolBlk_ptr->Grid.Cell[i  ][j  ].A*SolBlk_ptr->Sp[i  ][j  ]*
                                      W[(search_directions)*scalar_dim + index(i  ,j  ,k)]+
                                      SolBlk_ptr->Grid.Cell[i+1][j  ].A*SolBlk_ptr->Sp[i+1][j  ]*
                                      W[(search_directions)*scalar_dim + index(i+1,j  ,k)]+
                                      SolBlk_ptr->Grid.Cell[i  ][j+1].A*SolBlk_ptr->Sp[i  ][j+1]*
                                      W[(search_directions)*scalar_dim + index(i  ,j+1,k)]+
                                      SolBlk_ptr->Grid.Cell[i+1][j+1].A*SolBlk_ptr->Sp[i+1][j+1]*
                                      W[(search_directions)*scalar_dim + index(i+1,j+1,k)])/
                                     (SolBlk_ptr->Grid.Cell[i  ][j  ].A*SolBlk_ptr->Sp[i  ][j  ]+
                                      SolBlk_ptr->Grid.Cell[i+1][j  ].A*SolBlk_ptr->Sp[i+1][j  ]+
                                      SolBlk_ptr->Grid.Cell[i  ][j+1].A*SolBlk_ptr->Sp[i  ][j+1]+
                                      SolBlk_ptr->Grid.Cell[i+1][j+1].A*SolBlk_ptr->Sp[i+1][j+1]);
/*               buffer[buffer_count] = (SolBlk_ptr->Grid.Cell[i  ][j  ].A* */
/*                                       W[(search_directions)*scalar_dim + index(i  ,j  ,k)]+ */
/*                                       SolBlk_ptr->Grid.Cell[i+1][j  ].A* */
/*                                       W[(search_directions)*scalar_dim + index(i+1,j  ,k)]+ */
/*                                       SolBlk_ptr->Grid.Cell[i  ][j+1].A* */
/*                                       W[(search_directions)*scalar_dim + index(i  ,j+1,k)]+ */
/*                                       SolBlk_ptr->Grid.Cell[i+1][j+1].A* */
/*                                       W[(search_directions)*scalar_dim + index(i+1,j+1,k)]); */
	   } else {
              buffer[buffer_count] = (SolBlk_ptr->Grid.Cell[i  ][j  ].A*SolBlk_ptr->Sp[i  ][j  ]*x[index(i  ,j  ,k)]+
                                      SolBlk_ptr->Grid.Cell[i+1][j  ].A*SolBlk_ptr->Sp[i+1][j  ]*x[index(i+1,j  ,k)]+
                                      SolBlk_ptr->Grid.Cell[i  ][j+1].A*SolBlk_ptr->Sp[i  ][j+1]*x[index(i  ,j+1,k)]+
                                      SolBlk_ptr->Grid.Cell[i+1][j+1].A*SolBlk_ptr->Sp[i+1][j+1]*x[index(i+1,j+1,k)])/
                                     (SolBlk_ptr->Grid.Cell[i  ][j  ].A*SolBlk_ptr->Sp[i  ][j  ]+
                                      SolBlk_ptr->Grid.Cell[i+1][j  ].A*SolBlk_ptr->Sp[i+1][j  ]+
                                      SolBlk_ptr->Grid.Cell[i  ][j+1].A*SolBlk_ptr->Sp[i  ][j+1]+
                                      SolBlk_ptr->Grid.Cell[i+1][j+1].A*SolBlk_ptr->Sp[i+1][j+1]);
/*               buffer[buffer_count] = (SolBlk_ptr->Grid.Cell[i  ][j  ].A*x[index(i  ,j  ,k)]+ */
/*                                       SolBlk_ptr->Grid.Cell[i+1][j  ].A*x[index(i+1,j  ,k)]+ */
/*                                       SolBlk_ptr->Grid.Cell[i  ][j+1].A*x[index(i  ,j+1,k)]+ */
/*                                       SolBlk_ptr->Grid.Cell[i+1][j+1].A*x[index(i+1,j+1,k)]); */
	   } /* endif */
        } /* endfor */
     } /* endfor */
  } /* endfor */
  return(0);
}



/*******************************************************************************
 * GMRES_Block::LoadSendBuffer_C2F -- Loads send message buffer for            *
 *                                    coarse to fine block message             *
 *                                    passing.                                 *
 *******************************************************************************/
template <> 
inline int GMRES_Block<Rte2D_State,Rte2D_Quad_Block,Rte2D_Input_Parameters>::
LoadSendBuffer_C2F(double *buffer,
		   int &buffer_count,
		   const int buffer_size,
		   const int i_min, 
		   const int i_max,
		   const int i_inc,
		   const int j_min, 
		   const int j_max,
		   const int j_inc) {
  int i, j, k;
  Vector2D dX;
  Rte2D_State Wcoarse, Wfine;

  if (j_min == j_max) { // North or south boundary.
     // Four different orderings to consider depending on the value of i_inc & j_inc.
     if (j_inc > 0) {
        if (i_inc > 0) {
           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
              // Perform limited linear least squares reconstruction in cell (i, j_min).
              SubcellReconstruction(i, j_min, LIMITER_ZERO);
              // Evaluate SW sub (fine) cell values.
              for (k = 1 ; k <= blocksize; ++ k) {
   	         if (vector_switch) {
	            Wcoarse[k] = W[(search_directions)*scalar_dim + index(i,j_min,k-1)];
	         } else {
	            Wcoarse[k] = x[index(i,j_min,k-1)];
	         } /* endif */
	      } /* endfor */
              dX = (SolBlk_ptr->Grid.Node[i][j_min].X+
                    HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min].X)+
                    HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i][j_min+1].X)+
                    SolBlk_ptr->Grid.Cell[i][j_min].Xc)/FOUR -
                   SolBlk_ptr->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dUdx[i][j_min])*dX.x +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dUdy[i][j_min])*dX.y;
              for (k = 1 ; k <= blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
              // Evaluate SE sub (fine) cell values.
              dX = (HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min].X)+
                    SolBlk_ptr->Grid.Node[i+1][j_min].X+
                    SolBlk_ptr->Grid.Cell[i][j_min].Xc+
                    HALF*(SolBlk_ptr->Grid.Node[i+1][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X))/FOUR -
                   SolBlk_ptr->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dUdx[i][j_min])*dX.x +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dUdy[i][j_min])*dX.y;
              for ( k = 1 ; k <= blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
           } /* endfor */
           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
              // Evaluate NW sub (fine) cell values.
              for (k = 1 ; k <= blocksize; ++ k) {
   	         if (vector_switch) {
	            Wcoarse[k] = W[(search_directions)*scalar_dim + index(i,j_min,k-1)];
	         } else {
	            Wcoarse[k] = x[index(i,j_min,k-1)];
	         } /* endif */
	      } /* endfor */
              dX = (HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i][j_min+1].X)+
                    SolBlk_ptr->Grid.Cell[i][j_min].Xc+
                    SolBlk_ptr->Grid.Node[i][j_min+1].X+
                    HALF*(SolBlk_ptr->Grid.Node[i][j_min+1].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X))/FOUR -
                   SolBlk_ptr->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dUdx[i][j_min])*dX.x +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dUdy[i][j_min])*dX.y;
              for (k = 1 ; k <= blocksize; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
              // Evaluate NE sub (fine) cell values.
              dX = (SolBlk_ptr->Grid.Cell[i][j_min].Xc+
                    HALF*(SolBlk_ptr->Grid.Node[i+1][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X)+
                    HALF*(SolBlk_ptr->Grid.Node[i][j_min+1].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X)+
                    SolBlk_ptr->Grid.Node[i+1][j_min+1].X)/FOUR -
                   SolBlk_ptr->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dUdx[i][j_min])*dX.x +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dUdy[i][j_min])*dX.y;
              for (k = 1 ; k <= blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
           } /* endfor */
        } else {
           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
              // Perform limited linear least squares reconstruction in cell (i, j_min).
              SubcellReconstruction(i, j_min, LIMITER_ZERO);
              // Evaluate SE sub (fine) cell values.
              for (k = 1 ; k <= blocksize; ++ k) {
   	         if (vector_switch) {
	            Wcoarse[k] = W[(search_directions)*scalar_dim + index(i,j_min,k-1)];
	         } else {
	            Wcoarse[k] = x[index(i,j_min,k-1)];
	         } /* endif */
	      } /* endfor */
              dX = (HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min].X)+
                    SolBlk_ptr->Grid.Node[i+1][j_min].X+
                    SolBlk_ptr->Grid.Cell[i][j_min].Xc+
                    HALF*(SolBlk_ptr->Grid.Node[i+1][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X))/FOUR -
                   SolBlk_ptr->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dUdx[i][j_min])*dX.x +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dUdy[i][j_min])*dX.y;
              for (k = 1 ; k <= blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
              // Evaluate SW sub (fine) cell values.
              dX = (SolBlk_ptr->Grid.Node[i][j_min].X+
                    HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min].X)+
                    HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i][j_min+1].X)+
                    SolBlk_ptr->Grid.Cell[i][j_min].Xc)/FOUR -
                   SolBlk_ptr->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dUdx[i][j_min])*dX.x +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dUdy[i][j_min])*dX.y;
              for ( k = 1 ; k <= blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
           } /* endfor */
           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
              // Evaluate NE sub (fine) cell values.
              for (k = 1 ; k <= blocksize; ++ k) {
   	         if (vector_switch) {
	            Wcoarse[k] = W[(search_directions)*scalar_dim + index(i,j_min,k-1)];
	         } else {
	            Wcoarse[k] = x[index(i,j_min,k-1)];
	         } /* endif */
	      } /* endfor */
              dX = (SolBlk_ptr->Grid.Cell[i][j_min].Xc+
                    HALF*(SolBlk_ptr->Grid.Node[i+1][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X)+
                    HALF*(SolBlk_ptr->Grid.Node[i][j_min+1].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X)+
                    SolBlk_ptr->Grid.Node[i+1][j_min+1].X)/FOUR -
                   SolBlk_ptr->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dUdx[i][j_min])*dX.x +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dUdy[i][j_min])*dX.y;
              for (k = 1 ; k <= blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
              // Evaluate NW sub (fine) cell values.
              dX = (HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i][j_min+1].X)+
                    SolBlk_ptr->Grid.Cell[i][j_min].Xc+
                    SolBlk_ptr->Grid.Node[i][j_min+1].X+
                    HALF*(SolBlk_ptr->Grid.Node[i][j_min+1].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X))/FOUR -
                   SolBlk_ptr->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dUdx[i][j_min])*dX.x +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dUdy[i][j_min])*dX.y;
              for (k = 1 ; k <= blocksize; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
           } /* endfor */
        } /* endif */
     } else {
        if (i_inc > 0) {
           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
              // Perform limited linear least squares reconstruction in cell (i, j_min).
              SubcellReconstruction(i, j_min, LIMITER_ZERO);
              // Evaluate NW sub (fine) cell values.
              for (k = 1 ; k <= blocksize; ++ k) {
   	         if (vector_switch) {
	            Wcoarse[k] = W[(search_directions)*scalar_dim + index(i,j_min,k-1)];
	         } else {
	            Wcoarse[k] = x[index(i,j_min,k-1)];
	         } /* endif */
	      } /* endfor */
              dX = (HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i][j_min+1].X)+
                    SolBlk_ptr->Grid.Cell[i][j_min].Xc+
                    SolBlk_ptr->Grid.Node[i][j_min+1].X+
                    HALF*(SolBlk_ptr->Grid.Node[i][j_min+1].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X))/FOUR -
                   SolBlk_ptr->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dUdx[i][j_min])*dX.x +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dUdy[i][j_min])*dX.y;
              for (k = 1 ; k <= blocksize; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
              // Evaluate NE sub (fine) cell values.
              dX = (SolBlk_ptr->Grid.Cell[i][j_min].Xc+
                    HALF*(SolBlk_ptr->Grid.Node[i+1][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X)+
                    HALF*(SolBlk_ptr->Grid.Node[i][j_min+1].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X)+
                    SolBlk_ptr->Grid.Node[i+1][j_min+1].X)/FOUR -
                   SolBlk_ptr->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dUdx[i][j_min])*dX.x +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dUdy[i][j_min])*dX.y;
              for (k = 1 ; k <= blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
           } /* endfor */
           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
              // Evaluate SW sub (fine) cell values.
             for (k = 1 ; k <= blocksize; ++ k) {
   	         if (vector_switch) {
	            Wcoarse[k] = W[(search_directions)*scalar_dim + index(i,j_min,k-1)];
	         } else {
	            Wcoarse[k] = x[index(i,j_min,k-1)];
	         } /* endif */
	      } /* endfor */
              dX = (SolBlk_ptr->Grid.Node[i][j_min].X+
                    HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min].X)+
                    HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i][j_min+1].X)+
                    SolBlk_ptr->Grid.Cell[i][j_min].Xc)/FOUR -
                   SolBlk_ptr->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dUdx[i][j_min])*dX.x +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dUdy[i][j_min])*dX.y;
              for ( k = 1 ; k <= blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
              // Evaluate SE sub (fine) cell values.
              dX = (HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min].X)+
                    SolBlk_ptr->Grid.Node[i+1][j_min].X+
                    SolBlk_ptr->Grid.Cell[i][j_min].Xc+
                    HALF*(SolBlk_ptr->Grid.Node[i+1][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X))/FOUR -
                   SolBlk_ptr->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dUdx[i][j_min])*dX.x +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dUdy[i][j_min])*dX.y;
              for (k = 1 ; k <= blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
           } /* endfor */
        } else {
           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
              // Perform limited linear least squares reconstruction in cell (i, j_min).
              SubcellReconstruction(i, j_min, LIMITER_ZERO);
              // Evaluate NE sub (fine) cell values.
              for (k = 1 ; k <= blocksize; ++ k) {
   	         if (vector_switch) {
	            Wcoarse[k] = W[(search_directions)*scalar_dim + index(i,j_min,k-1)];
	         } else {
	            Wcoarse[k] = x[index(i,j_min,k-1)];
	         } /* endif */
	      } /* endfor */
              dX = (SolBlk_ptr->Grid.Cell[i][j_min].Xc+
                    HALF*(SolBlk_ptr->Grid.Node[i+1][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X)+
                    HALF*(SolBlk_ptr->Grid.Node[i][j_min+1].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X)+
                    SolBlk_ptr->Grid.Node[i+1][j_min+1].X)/FOUR -
                   SolBlk_ptr->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dUdx[i][j_min])*dX.x +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dUdy[i][j_min])*dX.y;

	      // Wfine = Wfine.U();  ??MISTAKE MADE IN KALVINS ORIGINAL ?????

             for (k = 1 ; k <= blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
              // Evaluate NW sub (fine) cell values.
              dX = (HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i][j_min+1].X)+
                    SolBlk_ptr->Grid.Cell[i][j_min].Xc+
                    SolBlk_ptr->Grid.Node[i][j_min+1].X+
                    HALF*(SolBlk_ptr->Grid.Node[i][j_min+1].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X))/FOUR -
                   SolBlk_ptr->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dUdx[i][j_min])*dX.x +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dUdy[i][j_min])*dX.y;
              for (k = 1 ; k <= blocksize; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
           } /* endfor */
           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
              // Evaluate SE sub (fine) cell values.
              for (k = 1 ; k <= blocksize; ++ k) {
   	         if (vector_switch) {
	            Wcoarse[k] = W[(search_directions)*scalar_dim + index(i,j_min,k-1)];
	         } else {
	            Wcoarse[k] = x[index(i,j_min,k-1)];
	         } /* endif */
	      } /* endfor */
              dX = (HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min].X)+
                    SolBlk_ptr->Grid.Node[i+1][j_min].X+
                    SolBlk_ptr->Grid.Cell[i][j_min].Xc+
                    HALF*(SolBlk_ptr->Grid.Node[i+1][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X))/FOUR -
                   SolBlk_ptr->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dUdx[i][j_min])*dX.x +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dUdy[i][j_min])*dX.y;
              for ( k = 1 ; k <= blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
              // Evaluate SW sub (fine) cell values.
              dX = (SolBlk_ptr->Grid.Node[i][j_min].X+
                    HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min].X)+
                    HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i][j_min+1].X)+
                    SolBlk_ptr->Grid.Cell[i][j_min].Xc)/FOUR -
                   SolBlk_ptr->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dUdx[i][j_min])*dX.x +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dUdy[i][j_min])*dX.y;
              for (k = 1 ; k <= blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
           } /* endfor */
        } /* endif */
     } /* endif */
  } else { // East or west boundary.
     // Four different orderings to consider depending on the value of i_inc & j_inc.
     if (j_inc > 0) {
        if (i_inc > 0) {
           for ( j = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
              // Perform limited linear least squares reconstruction in cell (i_min, j).
              SubcellReconstruction(i_min, j, LIMITER_ZERO);
              // Evaluate SW sub (fine) cell values.
              for (k = 1 ; k <= blocksize; ++ k) {
   	         if (vector_switch) {
	            Wcoarse[k] = W[(search_directions)*scalar_dim + index(i_min,j,k-1)];
	         } else {
	            Wcoarse[k] = x[index(i_min,j,k-1)];
	         } /* endif */
	      } /* endfor */
              dX = (SolBlk_ptr->Grid.Node[i_min][j].X+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min+1][j].X)+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min][j+1].X)+
                    SolBlk_ptr->Grid.Cell[i_min][j].Xc)/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dUdx[i_min][j])*dX.x +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dUdy[i_min][j])*dX.y;
              for (k = 1 ; k <= blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
              // Evaluate SE sub (fine) cell values.
              dX = (HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min+1][j].X)+
                    SolBlk_ptr->Grid.Node[i_min+1][j].X+
                    SolBlk_ptr->Grid.Cell[i_min][j].Xc+
                    HALF*(SolBlk_ptr->Grid.Node[i_min+1][j].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X))/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dUdx[i_min][j])*dX.x +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dUdy[i_min][j])*dX.y;
              for (k = 1 ; k <= blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
              // Evaluate NW sub (fine) cell values.
              dX = (HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min][j+1].X)+
                    SolBlk_ptr->Grid.Cell[i_min][j].Xc+
                    SolBlk_ptr->Grid.Node[i_min][j+1].X+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j+1].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X))/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dUdx[i_min][j])*dX.x +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dUdy[i_min][j])*dX.y;
              for (k = 1 ; k <= blocksize; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
              // Evaluate NE sub (fine) cell values.
              dX = (SolBlk_ptr->Grid.Cell[i_min][j].Xc+
                    HALF*(SolBlk_ptr->Grid.Node[i_min+1][j].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X)+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j+1].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X)+
                    SolBlk_ptr->Grid.Node[i_min+1][j+1].X)/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dUdx[i_min][j])*dX.x +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dUdy[i_min][j])*dX.y;
              for (k = 1 ; k <= blocksize; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
           } /* endfor */
        } else {
           for ( j = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
              // Perform limited linear least squares reconstruction in cell (i_min, j).
              SubcellReconstruction(i_min, j, LIMITER_ZERO);
              // Evaluate SE sub (fine) cell values.
              for (k = 1 ; k <= blocksize; ++ k) {
   	         if (vector_switch) {
	            Wcoarse[k] = W[(search_directions)*scalar_dim + index(i_min,j,k-1)];
	         } else {
	            Wcoarse[k] = x[index(i_min,j,k-1)];
	         } /* endif */
	      } /* endfor */
              dX = (HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min+1][j].X)+
                    SolBlk_ptr->Grid.Node[i_min+1][j].X+
                    SolBlk_ptr->Grid.Cell[i_min][j].Xc+
                    HALF*(SolBlk_ptr->Grid.Node[i_min+1][j].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X))/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dUdx[i_min][j])*dX.x +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dUdy[i_min][j])*dX.y;
              for (k = 1 ; k <= blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
              // Evaluate SW sub (fine) cell values.
              dX = (SolBlk_ptr->Grid.Node[i_min][j].X+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min+1][j].X)+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min][j+1].X)+
                    SolBlk_ptr->Grid.Cell[i_min][j].Xc)/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dUdx[i_min][j])*dX.x +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dUdy[i_min][j])*dX.y;
              for (k = 1 ; k <= blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
              // Evaluate NE sub (fine) cell values.
              dX = (SolBlk_ptr->Grid.Cell[i_min][j].Xc+
                    HALF*(SolBlk_ptr->Grid.Node[i_min+1][j].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X)+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j+1].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X)+
                    SolBlk_ptr->Grid.Node[i_min+1][j+1].X)/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dUdx[i_min][j])*dX.x +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dUdy[i_min][j])*dX.y;
              for (k = 1 ; k <= blocksize; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
              // Evaluate NW sub (fine) cell values.
              dX = (HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min][j+1].X)+
                    SolBlk_ptr->Grid.Cell[i_min][j].Xc+
                    SolBlk_ptr->Grid.Node[i_min][j+1].X+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j+1].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X))/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dUdx[i_min][j])*dX.x +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dUdy[i_min][j])*dX.y;
              for (k = 1 ; k <= blocksize; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
           } /* endfor */
        } /* endif */
     } else {
        if (i_inc > 0) {
           for ( j = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
              // Perform limited linear least squares reconstruction in cell (i_min, j).
              SubcellReconstruction(i_min, j, LIMITER_ZERO);
              // Evaluate NW sub (fine) cell values.
              for (k = 1 ; k <= blocksize; ++ k) {
   	         if (vector_switch) {
	            Wcoarse[k] = W[(search_directions)*scalar_dim + index(i_min,j,k-1)];
	         } else {
	            Wcoarse[k] = x[index(i_min,j,k-1)];
	         } /* endif */
	      } /* endfor */
              dX = (HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min][j+1].X)+
                    SolBlk_ptr->Grid.Cell[i_min][j].Xc+
                    SolBlk_ptr->Grid.Node[i_min][j+1].X+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j+1].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X))/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dUdx[i_min][j])*dX.x +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dUdy[i_min][j])*dX.y;
              for (k = 1 ; k <= blocksize; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
              // Evaluate NE sub (fine) cell values.
              dX = (SolBlk_ptr->Grid.Cell[i_min][j].Xc+
                    HALF*(SolBlk_ptr->Grid.Node[i_min+1][j].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X)+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j+1].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X)+
                    SolBlk_ptr->Grid.Node[i_min+1][j+1].X)/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dUdx[i_min][j])*dX.x +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dUdy[i_min][j])*dX.y;
              for (k = 1 ; k <= blocksize; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
              // Evaluate SW sub (fine) cell values.
              dX = (SolBlk_ptr->Grid.Node[i_min][j].X+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min+1][j].X)+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min][j+1].X)+
                    SolBlk_ptr->Grid.Cell[i_min][j].Xc)/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dUdx[i_min][j])*dX.x +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dUdy[i_min][j])*dX.y;
              for (k = 1 ; k <= blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
              // Evaluate SE sub (fine) cell values.
              dX = (HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min+1][j].X)+
                    SolBlk_ptr->Grid.Node[i_min+1][j].X+
                    SolBlk_ptr->Grid.Cell[i_min][j].Xc+
                    HALF*(SolBlk_ptr->Grid.Node[i_min+1][j].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X))/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dUdx[i_min][j])*dX.x +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dUdy[i_min][j])*dX.y;
              for (k = 1 ; k <= blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
           } /* endfor */
        } else {
           for ( j = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
              // Perform limited linear least squares reconstruction in cell (i_min, j).
              SubcellReconstruction(i_min, j, LIMITER_ZERO);
              // Evaluate NE sub (fine) cell values.
              for (k = 1 ; k <= blocksize; ++ k) {
   	         if (vector_switch) {
	            Wcoarse[k] = W[(search_directions)*scalar_dim + index(i_min,j,k-1)];
	         } else {
	            Wcoarse[k] = x[index(i_min,j,k-1)];
	         } /* endif */
	      } /* endfor */
              dX = (SolBlk_ptr->Grid.Cell[i_min][j].Xc+
                    HALF*(SolBlk_ptr->Grid.Node[i_min+1][j].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X)+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j+1].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X)+
                    SolBlk_ptr->Grid.Node[i_min+1][j+1].X)/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dUdx[i_min][j])*dX.x +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dUdy[i_min][j])*dX.y;
              for (k = 1 ; k <= blocksize; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
              // Evaluate NW sub (fine) cell values.
              dX = (HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min][j+1].X)+
                    SolBlk_ptr->Grid.Cell[i_min][j].Xc+
                    SolBlk_ptr->Grid.Node[i_min][j+1].X+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j+1].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X))/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dUdx[i_min][j])*dX.x +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dUdy[i_min][j])*dX.y;
              for (k = 1 ; k <= blocksize; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
              // Evaluate SE sub (fine) cell values.
              dX = (HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min+1][j].X)+
                    SolBlk_ptr->Grid.Node[i_min+1][j].X+
                    SolBlk_ptr->Grid.Cell[i_min][j].Xc+
                    HALF*(SolBlk_ptr->Grid.Node[i_min+1][j].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X))/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dUdx[i_min][j])*dX.x +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dUdy[i_min][j])*dX.y;
              for (k = 1 ; k <= blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
              // Evaluate SW sub (fine) cell values.
              dX = (SolBlk_ptr->Grid.Node[i_min][j].X+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min+1][j].X)+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min][j+1].X)+
                    SolBlk_ptr->Grid.Cell[i_min][j].Xc)/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dUdx[i_min][j])*dX.x +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dUdy[i_min][j])*dX.y;
              for (k = 1 ; k <= blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
           } /* endfor */
        } /* endif */
     } /* endif */
  } /* endif */
  return(0);
}



/**************************************************************************
 * GMRES_Block::SubcellReconstruction --                                  *
 *              Performs the subcell reconstruction of solution state     *
 *              within a given cell (i,j) of the computational mesh for   *
 *              the specified quadrilateral solution block.               *
 **************************************************************************/
template <> 
inline void GMRES_Block<Rte2D_State,Rte2D_Quad_Block,Rte2D_Input_Parameters>::
SubcellReconstruction(const int i, 
		      const int j,
		      const int Limiter) {
  
  int n, n2, n_pts, i_index[8], j_index[8], k;
  double u0, u0Min, u0Max, uQuad[4], phi_n;
  double DxDx_ave, DxDy_ave, DyDy_ave;
  Vector2D dX;
  Rte2D_State U0, DU, DUDx_ave, DUDy_ave, W_VACUUM;
  W_VACUUM.Zero();
  
  /* Carry out the limited solution reconstruction in
     each cell of the computational mesh. */

   if (i == ICl-Nghost || i == ICu+Nghost ||
      j == JCl-Nghost || j == JCu+Nghost) {
    n_pts = 0;
    } else if ((i == ICl-Nghost+1) && 
               (Grid.BCtypeW[j] != BC_NONE)) {
      if (j == JCl-Nghost+1 || j == JCu+Nghost-1) {
         n_pts = 0;
      } else if (SolBlk_ptr->Grid.BCtypeW[j] == BC_PERIODIC ||
                 SolBlk_ptr->Grid.BCtypeW[j] == BC_CONSTANT_EXTRAPOLATION ||
                 SolBlk_ptr->Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION ||
                 SolBlk_ptr->Grid.BCtypeW[j] == BC_GRAY_WALL ) {
         if (j == JCl) {
            n_pts = 5;
            i_index[0] = i-1; j_index[0] = j  ;
            i_index[1] = i+1; j_index[1] = j  ;
            i_index[2] = i-1; j_index[2] = j+1;
            i_index[3] = i  ; j_index[3] = j+1;
            i_index[4] = i+1; j_index[4] = j+1;
         } else if (j == JCu) {
            n_pts = 5;
            i_index[0] = i-1; j_index[0] = j-1;
            i_index[1] = i  ; j_index[1] = j-1;
            i_index[2] = i+1; j_index[2] = j-1;
            i_index[3] = i-1; j_index[3] = j  ;
            i_index[4] = i+1; j_index[4] = j  ;
         } else {
            n_pts = 8;
            i_index[0] = i-1; j_index[0] = j-1;
            i_index[1] = i  ; j_index[1] = j-1;
            i_index[2] = i+1; j_index[2] = j-1;
            i_index[3] = i-1; j_index[3] = j  ;
            i_index[4] = i+1; j_index[4] = j  ;
            i_index[5] = i-1; j_index[5] = j+1;
            i_index[6] = i  ; j_index[6] = j+1;
            i_index[7] = i+1; j_index[7] = j+1;
         } /* endif */
      } else {
         if (j == JCl) {
            n_pts = 3;
            i_index[0] = i+1; j_index[0] = j  ;
            i_index[1] = i  ; j_index[1] = j+1;
            i_index[2] = i+1; j_index[2] = j+1;
         } else if (j == JCu) {
            n_pts = 3;
            i_index[0] = i  ; j_index[0] = j-1;
            i_index[1] = i+1; j_index[1] = j-1;
            i_index[2] = i+1; j_index[2] = j  ;
         } else {
            n_pts = 5;
            i_index[0] = i  ; j_index[0] = j-1;
            i_index[1] = i+1; j_index[1] = j-1;
            i_index[2] = i+1; j_index[2] = j  ;
            i_index[3] = i  ; j_index[3] = j+1;
            i_index[4] = i+1; j_index[4] = j+1;
         } /* endif */
      } /* endif */           
    } else if ((i == ICu+Nghost-1) && 
               (SolBlk_ptr->Grid.BCtypeE[j] != BC_NONE)) {
      if (j == JCl-Nghost+1 || j == JCu+Nghost-1) {
         n_pts = 0;
      } else if (SolBlk_ptr->Grid.BCtypeE[j] == BC_PERIODIC ||
                 SolBlk_ptr->Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
                 SolBlk_ptr->Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
                 SolBlk_ptr->Grid.BCtypeE[j] == BC_GRAY_WALL ) {
         if (j == JCl) {
            n_pts = 5;
            i_index[0] = i-1; j_index[0] = j  ;
            i_index[1] = i+1; j_index[1] = j  ;
            i_index[2] = i-1; j_index[2] = j+1;
            i_index[3] = i  ; j_index[3] = j+1;
            i_index[4] = i+1; j_index[4] = j+1;
         } else if (j == JCu) {
            n_pts = 5;
            i_index[0] = i-1; j_index[0] = j-1;
            i_index[1] = i  ; j_index[1] = j-1;
            i_index[2] = i+1; j_index[2] = j-1;
            i_index[3] = i-1; j_index[3] = j  ;
            i_index[4] = i+1; j_index[4] = j  ;
         } else {
            n_pts = 8;
            i_index[0] = i-1; j_index[0] = j-1;
            i_index[1] = i  ; j_index[1] = j-1;
            i_index[2] = i+1; j_index[2] = j-1;
            i_index[3] = i-1; j_index[3] = j  ;
            i_index[4] = i+1; j_index[4] = j  ;
            i_index[5] = i-1; j_index[5] = j+1;
            i_index[6] = i  ; j_index[6] = j+1;
            i_index[7] = i+1; j_index[7] = j+1;
         } /* endif */
      } else {
         if (j == JCl) {
            n_pts = 3;
            i_index[0] = i-1; j_index[0] = j  ;
            i_index[1] = i-1; j_index[1] = j+1;
            i_index[2] = i  ; j_index[2] = j+1;
         } else if (j == JCu) {
            n_pts = 3;
            i_index[0] = i-1; j_index[0] = j-1;
            i_index[1] = i  ; j_index[1] = j-1;
            i_index[2] = i-1; j_index[2] = j  ;
         } else {
            n_pts = 5;
            i_index[0] = i-1; j_index[0] = j-1;
            i_index[1] = i  ; j_index[1] = j-1;
            i_index[2] = i-1; j_index[2] = j  ;
            i_index[3] = i-1; j_index[3] = j+1;
            i_index[4] = i  ; j_index[4] = j+1;
         } /* endif */
      } /* endif */
    } else if ((j == JCl-Nghost+1) && 
               (SolBlk_ptr->Grid.BCtypeS[i] != BC_NONE)) {
      if (i == ICl-Nghost+1 || i == ICu+Nghost-1) {
         n_pts = 0;
      } else if (SolBlk_ptr->Grid.BCtypeS[i] == BC_PERIODIC ||
                 SolBlk_ptr->Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
                 SolBlk_ptr->Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
                 SolBlk_ptr->Grid.BCtypeS[i] == BC_GRAY_WALL ) {
         if (i == ICl) {
            n_pts = 5;
            i_index[0] = i  ; j_index[0] = j-1;
            i_index[1] = i+1; j_index[1] = j-1;
            i_index[2] = i+1; j_index[2] = j  ;
            i_index[3] = i  ; j_index[3] = j+1;
            i_index[4] = i+1; j_index[4] = j+1;
         } else if (i == ICu) {
            n_pts = 5;
            i_index[0] = i-1; j_index[0] = j-1;
            i_index[1] = i  ; j_index[1] = j-1;
            i_index[2] = i-1; j_index[2] = j  ;
            i_index[3] = i-1; j_index[3] = j+1;
            i_index[4] = i  ; j_index[4] = j+1;
         } else {
            n_pts = 8;
            i_index[0] = i-1; j_index[0] = j-1;
            i_index[1] = i  ; j_index[1] = j-1;
            i_index[2] = i+1; j_index[2] = j-1;
            i_index[3] = i-1; j_index[3] = j  ;
            i_index[4] = i+1; j_index[4] = j  ;
            i_index[5] = i-1; j_index[5] = j+1;
            i_index[6] = i  ; j_index[6] = j+1;
            i_index[7] = i+1; j_index[7] = j+1;
         } /* endif */
      } else {
         if (i == ICl) {
            n_pts = 3;
            i_index[0] = i+1; j_index[0] = j  ;
            i_index[1] = i  ; j_index[1] = j+1;
            i_index[2] = i+1; j_index[2] = j+1;
         } else if (i == ICu) {
            n_pts = 3;
            i_index[0] = i-1; j_index[0] = j  ;
            i_index[1] = i-1; j_index[1] = j+1;
            i_index[2] = i  ; j_index[2] = j+1;
         } else {
            n_pts = 5;
            i_index[0] = i-1; j_index[0] = j  ;
            i_index[1] = i+1; j_index[1] = j  ;
            i_index[2] = i-1; j_index[2] = j+1;
            i_index[3] = i  ; j_index[3] = j+1;
            i_index[4] = i+1; j_index[4] = j+1;
         } /* endif */
      } /* endif */
    } else if ((j == JCu+Nghost-1) && 
               (SolBlk_ptr->Grid.BCtypeN[i] != BC_NONE)) {
      if (i == ICl-Nghost+1 || i == ICu+Nghost-1) {
         n_pts = 0;
      } else if (SolBlk_ptr->Grid.BCtypeN[i] == BC_PERIODIC ||
                 SolBlk_ptr->Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
                 SolBlk_ptr->Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
                 SolBlk_ptr->Grid.BCtypeN[i] == BC_GRAY_WALL ) {
         if (i == ICl) {
            n_pts = 5;
            i_index[0] = i  ; j_index[0] = j-1;
            i_index[1] = i+1; j_index[1] = j-1;
            i_index[2] = i+1; j_index[2] = j  ;
            i_index[3] = i  ; j_index[3] = j+1;
            i_index[4] = i+1; j_index[4] = j+1;
         } else if (i == ICu) {
            n_pts = 5;
            i_index[0] = i-1; j_index[0] = j-1;
            i_index[1] = i  ; j_index[1] = j-1;
            i_index[2] = i-1; j_index[2] = j  ;
            i_index[3] = i-1; j_index[3] = j+1;
            i_index[4] = i  ; j_index[4] = j+1;
         } else {
            n_pts = 8;
            i_index[0] = i-1; j_index[0] = j-1;
            i_index[1] = i  ; j_index[1] = j-1;
            i_index[2] = i+1; j_index[2] = j-1;
            i_index[3] = i-1; j_index[3] = j  ;
            i_index[4] = i+1; j_index[4] = j  ;
            i_index[5] = i-1; j_index[5] = j+1;
            i_index[6] = i  ; j_index[6] = j+1;
            i_index[7] = i+1; j_index[7] = j+1;
         } /* endif */
      } else {
         if (i == ICl) {
            n_pts = 3;
            i_index[0] = i  ; j_index[0] = j-1;
            i_index[1] = i+1; j_index[1] = j-1;
            i_index[2] = i+1; j_index[2] = j  ;
         } else if (i == ICu) {
            n_pts = 3;
            i_index[0] = i-1; j_index[0] = j-1;
            i_index[1] = i  ; j_index[1] = j-1;
            i_index[2] = i-1; j_index[2] = j  ;
         } else {
            n_pts = 5;
            i_index[0] = i-1; j_index[0] = j-1;
            i_index[1] = i  ; j_index[1] = j-1;
            i_index[2] = i+1; j_index[2] = j-1;
            i_index[3] = i-1; j_index[3] = j  ;
            i_index[4] = i+1; j_index[4] = j  ;
         } /* endif */
      } /* endif */
  } else {
    n_pts = 8;
    i_index[0] = i-1; j_index[0] = j-1;
    i_index[1] = i  ; j_index[1] = j-1;
    i_index[2] = i+1; j_index[2] = j-1;
    i_index[3] = i-1; j_index[3] = j  ;
    i_index[4] = i+1; j_index[4] = j  ;
    i_index[5] = i-1; j_index[5] = j+1;
    i_index[6] = i  ; j_index[6] = j+1;
    i_index[7] = i+1; j_index[7] = j+1;
  } /* endif */

  if (n_pts > 0) {
      DUDx_ave = W_VACUUM;
      DUDy_ave = W_VACUUM;
      DxDx_ave = ZERO;
      DxDy_ave = ZERO;
      DyDy_ave = ZERO;

      for ( k = 1 ; k <= blocksize; ++ k) {
         if (vector_switch) {
            U0[k] = W[(search_directions)*scalar_dim + index(i,j,k-1)];
         } else {
            U0[k] = x[index(i,j,k-1)];
         } /* endif */
      } /* endfor */

      for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
          dX = SolBlk_ptr->Grid.Cell[ i_index[n2] ][ j_index[n2] ].Xc - 
               SolBlk_ptr->Grid.Cell[i][j].Xc;
          for ( k = 1 ; k <= blocksize; ++ k) {
             if (vector_switch) {
                DU[k] = W[(search_directions)*scalar_dim + index(i_index[n2] , j_index[n2] , k-1)] - 
                        U0[k];
             } else {
                DU[k] = x[index( i_index[n2] , j_index[n2] , k-1)] - 
                        U0[k];
             } /* endif */
          } /* endfor */
          DUDx_ave += DU*dX.x;
          DUDy_ave += DU*dX.y;
          DxDx_ave += dX.x*dX.x;
          DxDy_ave += dX.x*dX.y;
          DyDy_ave += dX.y*dX.y;
      } /* endfor */
  					    
      DUDx_ave = DUDx_ave/double(n_pts);
      DUDy_ave = DUDy_ave/double(n_pts);
      DxDx_ave = DxDx_ave/double(n_pts);
      DxDy_ave = DxDy_ave/double(n_pts);
      DyDy_ave = DyDy_ave/double(n_pts);
      SolBlk_ptr->dUdx[i][j] = (DUDx_ave*DyDy_ave-DUDy_ave*DxDy_ave)/
                               (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
      SolBlk_ptr->dUdy[i][j] = (DUDy_ave*DxDx_ave-DUDx_ave*DxDy_ave)/
                               (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
  
      // Calculate slope limiters. 
      if (!SolBlk_ptr->Freeze_Limiter) {
         for ( n = 1 ; n <= blocksize ; ++n ) {
	    u0 = U0[n];
            u0Min = U0[n];
            u0Max = u0Min;
            for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
               if (vector_switch) {
                  u0Min = min(u0Min, W[(search_directions)*scalar_dim + index(i_index[n2] , j_index[n2] ,n-1)]);
                  u0Max = max(u0Max, W[(search_directions)*scalar_dim + index(i_index[n2] , j_index[n2] ,n-1)]);
               } else {
                  u0Min = min(u0Min, x[index(i_index[n2] , j_index[n2] , n-1)]);
                  u0Max = max(u0Max, x[index(i_index[n2] , j_index[n2] , n-1)]);
               } /* endif */
            } /* endfor */
    
            dX = SolBlk_ptr->Grid.xfaceE(i, j)-SolBlk_ptr->Grid.Cell[i][j].Xc;
            uQuad[0] = u0 + 
                       SolBlk_ptr->dUdx[i][j][n]*dX.x +
                       SolBlk_ptr->dUdy[i][j][n]*dX.y ;
            dX = SolBlk_ptr->Grid.xfaceW(i, j)-SolBlk_ptr->Grid.Cell[i][j].Xc;
            uQuad[1] = u0 + 
                       SolBlk_ptr->dUdx[i][j][n]*dX.x +
                       SolBlk_ptr->dUdy[i][j][n]*dX.y ;
            dX = SolBlk_ptr->Grid.xfaceN(i, j)-SolBlk_ptr->Grid.Cell[i][j].Xc;
            uQuad[2] = u0 + 
                       SolBlk_ptr->dUdx[i][j][n]*dX.x +
                       SolBlk_ptr->dUdy[i][j][n]*dX.y ;
            dX = SolBlk_ptr->Grid.xfaceS(i, j)-SolBlk_ptr->Grid.Cell[i][j].Xc;
            uQuad[3] = u0 + 
                       SolBlk_ptr->dUdx[i][j][n]*dX.x +
                       SolBlk_ptr->dUdy[i][j][n]*dX.y ;
    
            switch(Limiter) {
              case LIMITER_ONE :
                phi_n = ONE;
                break;
              case LIMITER_ZERO :
                phi_n = ZERO;
                break;
              case LIMITER_BARTH_JESPERSEN :
                phi_n = Limiter_BarthJespersen(uQuad, u0, 
                                               u0Min, u0Max, 4);
                break;
              case LIMITER_VENKATAKRISHNAN :
                phi_n = Limiter_Venkatakrishnan(uQuad, u0, 
                                                u0Min, u0Max, 4);
                break;
              case LIMITER_VANLEER :
                phi_n = Limiter_VanLeer(uQuad, u0, 
                                        u0Min, u0Max, 4);
                break;
              case LIMITER_VANALBADA :
                phi_n = Limiter_VanAlbada(uQuad, u0, 
                                          u0Min, u0Max, 4);
                break;
              default:
                phi_n = Limiter_BarthJespersen(uQuad, u0, 
                                               u0Min, u0Max, 4);
                break;
            } /* endswitch */

	    SolBlk_ptr->phi[i][j][n] = phi_n;

         } /* endfor */
      } /* endif */
  } else {
      SolBlk_ptr->dUdx[i][j] = W_VACUUM;
      SolBlk_ptr->dUdy[i][j] = W_VACUUM; 
      SolBlk_ptr->phi[i][j]  = W_VACUUM;
  } /* endif */

}





/****************************************************************/
/**** PRECONDTIONER REQUIRED SPECIALIZATIONS & FUNCTIONS ********/
/****************************************************************/  





// /*!**************************************************************
//  * Specialization of Block_Preconditioner::Preconditioner_dFdU  *
//  *                                                              *
//  * Calculates the dFdU matrix used to generate the approximate  *               
//  * Jacobian for the Block Preconditioner.                       *
//  ****************************************************************/
// template<> inline void Block_Preconditioner<Rte2D_State,
// 					    Rte2D_Quad_Block,	   
// 					    Rte2D_Input_Parameters>::
// Preconditioner_dFIdU(DenseMatrix &_dFdU, Rte2D_State U)
// {
//   dFdU(_dFdU,U);
// }


/*!**************************************************************
 * Specialization of Block_Preconditioner::Preconditioner_dSdU  *
 *                                                              *
 * Calculates the dFdU matrix used to generate the approximate  *               
 * Jacobian for the Block Preconditioner.                       *
 ****************************************************************/
template<> inline void Block_Preconditioner<Rte2D_State,
					    Rte2D_Quad_Block,					    
					    Rte2D_Input_Parameters>::
Preconditioner_dSdU(int ii, int jj, DenseMatrix &dRdU){
  //Source Terms 

 
  //Add Jacobian for axisymmetric source terms
  if (Soln_Block_ptr->Axisymmetric) {
    Soln_Block_ptr->Uo[ii][jj].dSadU(dRdU,Soln_Block_ptr->Sp[ii][jj]);
  }  

  //Add Jacobian for regular source terms
  Soln_Block_ptr->Uo[ii][jj].dSdU(dRdU);

//   cout << endl;
//   cout << endl << dRdU << endl;
//   cout << endl;

}


/*!**************************************************************
 * Specialization of Block_Preconditioner::                     *
 *                               normalize_Preconditioner_dFdU  *
 *                                                              *
 * Normaliazes the dFdU matrix used to generate the approximate *               
 * Jacobian for the Block Preconditioner.                       *
 ****************************************************************/
template<> inline void Block_Preconditioner<Rte2D_State,
					    Rte2D_Quad_Block,
					    Rte2D_Input_Parameters>::
normalize_Preconditioner_dFdU(DenseMatrix &dFdU) { /*EMPTY - NOT NORMALIZED*/ }


/********************************************************************************
 *  Jacobian (dR/dU) Approximation Stencil Allocation - only called on STARTUP   *
 *********************************************************************************/
template <>
void Block_Preconditioner<Rte2D_State,
			  Rte2D_Quad_Block,
			  Rte2D_Input_Parameters>::
Setup_Jacobian_approximation(){

  int block_mat_size = Soln_Block_ptr->NCi*Soln_Block_ptr->NCj;       //block matrix size based on icells*jcells
  int nnz;

  //! 1st order 2D Jacobian, Sparse Block Matrix Memory Allocation of non-zero block entries
  if( Jacobian_stencil_size == 5){        
    nnz = 5*block_mat_size - 2*(Soln_Block_ptr->NCi + Soln_Block_ptr->NCj); 
    //! 2nd order 2D Jacobian
  }  else if ( Jacobian_stencil_size == 9) {
    nnz = 9*block_mat_size - 6*(Soln_Block_ptr->NCi + Soln_Block_ptr->NCj) + 4;    
  }      

  //!Temporary Storage Arrays
  int *i_index = new int[nnz];      //location of dense local blocks, in global sparse block Matrix   
  int *j_index = new int[nnz];
  double *Data = new double [nnz*blocksize*blocksize];
  int *block_i = new int[Jacobian_stencil_size];
  
  //!Setup Stencil for approximate Jacobian
  int nnz_count = 0;   
  //! Loop through each row of Block Matrix
  for (int i=0; i< block_mat_size; i++){       
    Get_Block_Index(i,block_i);
    int stencil = 0;
    //! Determine which entries correspond to 1st, 2nd, etc. order stencil    
    //cout<<"\n i "<<i;
    while( stencil < Jacobian_stencil_size) { 
      int j = block_i[stencil];
      //cout<<" "<<j;
      if( j >= 0 && j < block_mat_size){	
	i_index[nnz_count] = i+1;  // bpkit assumes Fortran indexing ie. starting at 1 
	j_index[nnz_count] = j+1;
	for(int k=0; k<blocksize*blocksize; k++){ //initialize Block Matrix as identity matrix
	  if( i == j && k%(blocksize+1) == 0){
	    Data[nnz_count*blocksize*blocksize + k] = ONE;
	  } else {
	    Data[nnz_count*blocksize*blocksize + k] = ZERO;
	  }
	}	 
	nnz_count++;
      }
      stencil++;	
    }      
  }

  if( nnz != nnz_count){ cerr<<"\n Number of nonzero blocks mismatch error in approximate Jacobian formation "
			     <<nnz<<" != "<<nnz_count<<"\n"; exit(1); }

  //! Create Sparse bpkit "BlockMat" Block Matrix 
  Block_Jacobian_approx.setup(block_mat_size,nnz,i_index,j_index,Data,blocksize);                  

  
//   /****************************************************************************************/
//   // Low Mach # Preconditioing  Initialization
//   if(Input_Parameters->Preconditioning){
//     int JCl_overlap = 0, JCu_overlap = 0, ICu_overlap =0, ICl_overlap = 0;
//     if(Input_Parameters->NKS_IP.GMRES_Overlap){	
//       if (Soln_Block_ptr->Grid.BCtypeS[Soln_Block_ptr->ICl] == BC_NONE)  JCl_overlap = Input_Parameters->NKS_IP.GMRES_Overlap;
//       if (Soln_Block_ptr->Grid.BCtypeN[Soln_Block_ptr->ICu] == BC_NONE)  JCu_overlap = Input_Parameters->NKS_IP.GMRES_Overlap;
//       if (Soln_Block_ptr->Grid.BCtypeE[Soln_Block_ptr->JCu] == BC_NONE)  ICu_overlap = Input_Parameters->NKS_IP.GMRES_Overlap;
//       if (Soln_Block_ptr->Grid.BCtypeW[Soln_Block_ptr->JCl] == BC_NONE)  ICl_overlap = Input_Parameters->NKS_IP.GMRES_Overlap;
//     }
//     int isize = Soln_Block_ptr->NCi;
//     int jsize = Soln_Block_ptr->NCj;

//     Low_Mach_Number_Preconditioner = new  DenseMatrix*[isize]; 
//     for(int i= 0 ; i<isize; i++){  
//       Low_Mach_Number_Preconditioner[i] = new DenseMatrix[jsize];
//       for(int j= 0; j<jsize; j++){  
// 	Low_Mach_Number_Preconditioner[i][j] = DenseMatrix(blocksize,blocksize,ZERO);
//       }
//     }
//   }
//   /****************************************************************************************/

  
  //! Clean up local memory
  delete[] i_index; delete[] j_index; delete[] Data;  delete[] block_i;

}



/**********************************************************
 *  Update BlockMat with approximation to the Jacobian    *
 **********************************************************/
template <>
void Block_Preconditioner<Rte2D_State,
			  Rte2D_Quad_Block,
			  Rte2D_Input_Parameters>::
Update_Jacobian_and_Preconditioner()
{
  
  //!Local Variables and Temporary Storage
  int block_mat_size = Soln_Block_ptr->NCi*Soln_Block_ptr->NCj; 
  DenseMatrix *Jacobian_Data = new DenseMatrix[Jacobian_stencil_size];
  for(int i=0; i<Jacobian_stencil_size; i++) { Jacobian_Data[i] = DenseMatrix(blocksize,blocksize,ZERO); }
  int *block_i = new int[Jacobian_stencil_size]; 
  int *block_j = new int[Jacobian_stencil_size]; 

//   DenseMatrix *Jacobian_test = new DenseMatrix[Jacobian_stencil_size];
//   for(int i=0; i<Jacobian_stencil_size; i++) { Jacobian_test[i] = DenseMatrix(blocksize,blocksize,ZERO); }

  //! Initially assume no overlap
  int JCl_overlap = 0, JCu_overlap = 0, ICu_overlap =0, ICl_overlap = 0;
  
  //! If overlap determine which block boundaries are internal, ie. BC_NONE
  if(Input_Parameters->NKS_IP.GMRES_Overlap){	
    if (Soln_Block_ptr->Grid.BCtypeS[Soln_Block_ptr->ICl] == BC_NONE)  JCl_overlap = Input_Parameters->NKS_IP.GMRES_Overlap;
    if (Soln_Block_ptr->Grid.BCtypeN[Soln_Block_ptr->ICu] == BC_NONE)  JCu_overlap = Input_Parameters->NKS_IP.GMRES_Overlap;
    if (Soln_Block_ptr->Grid.BCtypeE[Soln_Block_ptr->JCu] == BC_NONE)  ICu_overlap = Input_Parameters->NKS_IP.GMRES_Overlap;
    if (Soln_Block_ptr->Grid.BCtypeW[Soln_Block_ptr->JCl] == BC_NONE)  ICl_overlap = Input_Parameters->NKS_IP.GMRES_Overlap;
  }

  //*********************************************************************************//
  /*! Calculate Jacobians for each cell and Update Global Jacobian Block Matrix
   * loop through all non-ghost cells, including overlap cells.  Ghost cells already
   * set to Zero or Identity by initialization.
   **********************************************************************************/
  for(int i= Soln_Block_ptr->ICl - ICl_overlap; i<= Soln_Block_ptr->ICu + ICu_overlap; i++){    
    for(int j= Soln_Block_ptr->JCl - JCl_overlap; j<= Soln_Block_ptr->JCu + ICu_overlap; j++){  
  
      //--------------------------------------------------------------------------//
      //! Calculate Local Approximate Jacobian                        
      switch(Input_Parameters->NKS_IP.Jacobian_Order){
      case SOURCE_TERMS_ONLY :
	Implicit_Euler(i,j, Jacobian_Data);
	Preconditioner_dSdU(i,j,Jacobian_Data[CENTER]);
// 	Jacobian_Data[CENTER].identity();
	break;
      case FIRST_ORDER_INVISCID_HLLE : 
	Implicit_Euler(i,j, Jacobian_Data);
 	First_Order_Inviscid_Jacobian_HLLE(i,j, Jacobian_Data); // Not really HLLE
	Preconditioner_dSdU(i,j,Jacobian_Data[CENTER]);   
	break;
      case FIRST_ORDER_INVISCID_ROE : 
	Implicit_Euler(i,j, Jacobian_Data);
 	First_Order_Inviscid_Jacobian_Roe(i,j, Jacobian_Data); // Not really Roe
	Preconditioner_dSdU(i,j,Jacobian_Data[CENTER]);   
	break;
      }
                  
      //--------------------------------------------------------------------------//
      //! Get Block Matrix locations that have components from a given Cell(i,j)
      Get_Block_Index(i,j, block_i, block_j);
      
      //fudge for iGhost Cell reset to zero   //jGhost cell already zero                
      if(block_i[NORTH] < TWO*Soln_Block_ptr->NCi) Jacobian_Data[NORTH].zero();
      if(block_i[SOUTH] > block_mat_size - TWO*Soln_Block_ptr->NCi) Jacobian_Data[SOUTH].zero();

      if(Jacobian_stencil_size == 9){		
	if(block_i[NE] < TWO*Soln_Block_ptr->NCi)    Jacobian_Data[NE].zero();
	if(block_i[NW] < TWO*Soln_Block_ptr->NCi)    Jacobian_Data[NW].zero(); 
	if(block_i[SE] > block_mat_size - TWO*Soln_Block_ptr->NCi)    Jacobian_Data[SE].zero();
	if(block_i[SW] > block_mat_size - TWO*Soln_Block_ptr->NCi)    Jacobian_Data[SW].zero();
      }




      //--------------------------------------------------------------------------//
      //! Update BlockMat with Local Approximate Jacobians 
      for( int block = 0; block < Jacobian_stencil_size; block++){
	// Normalize
	normalize_Preconditioner_dFdU(Jacobian_Data[block]);

	//can be speed up by more intelligent logic in bkpkit (BlockMat.cc  "setblock")
	Block_Jacobian_approx.setblock( block_i[block], block_j[block], DenseMatrix_to_DenseMat(Jacobian_Data[block]));

 	//DEBUG
	//if(block ==CENTER) cout<<"  Normalized Block "<<i<<","<<j<<" \n "<<Jacobian_Data[CENTER];

	Jacobian_Data[block].zero(); //Just in case to avoid +=/-= issues
      }     

    }
  }
  
  //DEBUG
  //  cout<<"\n Cell i,j "<<i<<","<<j<<" Blocks ";
  //  for(int k =0; k<Jacobian_stencil_size; k++) cout<<" "<<block_i[k]<<","<<block_j[k];  
  // cout<<"\n Precondtioner "<<Block_Jacobian_approx;

  //Local Memory cleanup
  delete[] Jacobian_Data; delete[] block_i; delete[] block_j;

  //Setup appropriate Preconditioner after Jacobian has been formed/Updated
  Setup_Preconditioner();
}



/*****************************************************************************
 *  Add finite time step to diagonal.                                        *
 *****************************************************************************/ 
template <>
inline void Block_Preconditioner<Rte2D_State,
				 Rte2D_Quad_Block,
				 Rte2D_Input_Parameters>::
Implicit_Euler(const int &cell_index_i,const int &cell_index_j, DenseMatrix* Jacobian){   
  //Low Mach # Preconditioning 
//   if(Input_Parameters->Preconditioning){       
//     double delta_n = min( TWO*(Soln_Block_ptr->Grid.Cell[cell_index_i][cell_index_j].A/
// 			       (Soln_Block_ptr->Grid.lfaceE(cell_index_i, cell_index_j)
// 				+ Soln_Block_ptr->Grid.lfaceW(cell_index_i, cell_index_j))),
// 			  TWO*(Soln_Block_ptr->Grid.Cell[cell_index_i][cell_index_j].A/
// 			       (Soln_Block_ptr->Grid.lfaceN(cell_index_i, cell_index_j)
// 				+Soln_Block_ptr->Grid.lfaceS(cell_index_i, cell_index_j))));         
    
//     Soln_Block_ptr->Uo[cell_index_i][cell_index_j].
//       Low_Mach_Number_Preconditioner(Low_Mach_Number_Preconditioner[cell_index_i][cell_index_j],
// 				     Soln_Block_ptr->Flow_Type,
// 				     delta_n);    

//     Jacobian[CENTER] -= Low_Mach_Number_Preconditioner[cell_index_i][cell_index_j]
//                           / (Soln_Block_ptr->dt[cell_index_i][cell_index_j]);

//     //store normalized
//     normalize_Preconditioner_dFdU( Low_Mach_Number_Preconditioner[cell_index_i][cell_index_j]);

//   } else { 

  // I/deltat
  DenseMatrix II(blocksize,blocksize);  
  II.identity();    
  Jacobian[CENTER] -= (II / (Soln_Block_ptr->dt[cell_index_i][cell_index_j]));

// }

}



/*!**************************************************************
 * Specialization of Block_Preconditioner::                     *
 *                           First_Order_Inviscid_Jacobian_HLLE *
 *                                                              *
 * Calculate First Order Local Jacobian Block(s) Coresponding   *
 * to Cell(i,j) using regular upwind flux function.  We are just*
 * using this function as a quick fix to fit in with the group  *
 * code.                                                        *
 ****************************************************************/
template<> inline void Block_Preconditioner<Rte2D_State,
					    Rte2D_Quad_Block,
					    Rte2D_Input_Parameters>::
First_Order_Inviscid_Jacobian_HLLE(const int &cell_index_i,const int &cell_index_j, 
				   DenseMatrix* Jacobian){              
  
  //! Caculate normal vectors -> in Vector2D format. 
  Vector2D nface_N = Soln_Block_ptr->Grid.nfaceN(cell_index_i,cell_index_j-1);
  Vector2D nface_S = Soln_Block_ptr->Grid.nfaceS(cell_index_i,cell_index_j+1);
  Vector2D nface_E = Soln_Block_ptr->Grid.nfaceE(cell_index_i-1,cell_index_j);
  Vector2D nface_W = Soln_Block_ptr->Grid.nfaceW(cell_index_i+1,cell_index_j);

  //! Calculate dFdU using solutions in the rotated frame -> matrix in DenseMatrix format. 
  DenseMatrix dFdU_N(blocksize,blocksize,ZERO); 
  DenseMatrix dFdU_S(blocksize,blocksize,ZERO); 
  DenseMatrix dFdU_E(blocksize,blocksize,ZERO); 
  DenseMatrix dFdU_W(blocksize,blocksize,ZERO); 

  //! Calculate Jacobian matrix -> blocksizexblocksize matrix in DenseMatrix format
  //Solution Rotate provided in pState 
  dFdU_n( dFdU_N, Soln_Block_ptr->U[cell_index_i][cell_index_j], nface_N); 
  dFdU_n( dFdU_S, Soln_Block_ptr->U[cell_index_i][cell_index_j], nface_S);
  dFdU_n( dFdU_E, Soln_Block_ptr->U[cell_index_i][cell_index_j], nface_E);
  dFdU_n( dFdU_W, Soln_Block_ptr->U[cell_index_i][cell_index_j], nface_W);

  //! Calculate Jacobian matrix -> blocksizexblocksize matrix in DenseMatrix format
  //North
  Jacobian[NORTH] = Soln_Block_ptr->Grid.lfaceN(cell_index_i,cell_index_j-1) 
    * Soln_Block_ptr->SpN[cell_index_i][cell_index_j-1] * dFdU_N; 

  //South
  Jacobian[SOUTH] = Soln_Block_ptr->Grid.lfaceS(cell_index_i,cell_index_j+1) 
    * Soln_Block_ptr->SpS[cell_index_i][cell_index_j+1] * dFdU_S;

  //East
  Jacobian[EAST] = Soln_Block_ptr->Grid.lfaceE(cell_index_i-1,cell_index_j) 
    * Soln_Block_ptr->SpE[cell_index_i-1][cell_index_j] * dFdU_E;

  //West
  Jacobian[WEST] = Soln_Block_ptr->Grid.lfaceW(cell_index_i+1,cell_index_j) 
    * Soln_Block_ptr->SpW[cell_index_i+1][cell_index_j] * dFdU_W;

  //Center calculated from neighbours
  //! Using the fact that dF/dU(right) = - dF/dU(left) 
  Jacobian[CENTER] += (Jacobian[NORTH] + Jacobian[SOUTH] + Jacobian[EAST]  + Jacobian[WEST])
                     /(Soln_Block_ptr->Grid.Cell[cell_index_i][cell_index_j].A * Soln_Block_ptr->Sp[cell_index_i][cell_index_j]);

  Jacobian[NORTH] = -Jacobian[NORTH]/(Soln_Block_ptr->Grid.Cell[cell_index_i][cell_index_j-1].A*Soln_Block_ptr->Sp[cell_index_i][cell_index_j-1]);
  Jacobian[SOUTH] = -Jacobian[SOUTH]/(Soln_Block_ptr->Grid.Cell[cell_index_i][cell_index_j+1].A*Soln_Block_ptr->Sp[cell_index_i][cell_index_j+1]);
  Jacobian[EAST] = -Jacobian[EAST]/(Soln_Block_ptr->Grid.Cell[cell_index_i-1][cell_index_j].A*Soln_Block_ptr->Sp[cell_index_i-1][cell_index_j]);
  Jacobian[WEST] = -Jacobian[WEST]/(Soln_Block_ptr->Grid.Cell[cell_index_i+1][cell_index_j].A*Soln_Block_ptr->Sp[cell_index_i+1][cell_index_j]);

 
}



/********************************************************
 * Routine: Flux Jacobian                               *
 *                                                      *
 *     Finite difference approximation                  *
 *                                                      *
 ********************************************************/
template<> inline void Block_Preconditioner<Rte2D_State,
					    Rte2D_Quad_Block,
					    Rte2D_Input_Parameters>::
Preconditioner_dFIdU_Roe(DenseMatrix& dRdU, const int ii, const int jj, const int Orient)
{

  int NUM_VAR_RTE2D = dRdU.get_n();   
  DenseMatrix dFidU(NUM_VAR_RTE2D, NUM_VAR_RTE2D,ZERO);
  
  int Ri, Rj;
  Vector2D nface;
  double lface;
  Rte2D_State Ul, Ur;   
  Rte2D_State UA, UB;   
  Rte2D_State FluxA, FluxB;   
  Ul.Zero();  Ur.Zero();
  double perturb = 5e-6;

  
  switch (Orient) {

    // NORTH
  case NORTH:
    Ri = ii; Rj=jj-1;
    nface = Soln_Block_ptr->Grid.nfaceN(Ri, Rj);
    lface = Soln_Block_ptr->Grid.lfaceN(Ri, Rj)*Soln_Block_ptr->SpN[Ri][Rj];
    Ur = Soln_Block_ptr->Uo[ii][jj];
    Ul = Soln_Block_ptr->Uo[Ri][Rj];
    break;

    // SOUTH
  case SOUTH:
    Ri = ii; Rj=jj+1;
    nface = Soln_Block_ptr->Grid.nfaceS(Ri, Rj);
    lface = Soln_Block_ptr->Grid.lfaceS(Ri, Rj)*Soln_Block_ptr->SpS[Ri][Rj];
    Ur = Soln_Block_ptr->Uo[ii][jj];
    Ul = Soln_Block_ptr->Uo[Ri][Rj];
    break;

    // EAST
  case EAST:
    Ri = ii-1; Rj=jj;
    nface = Soln_Block_ptr->Grid.nfaceE(Ri, Rj);     
    lface = Soln_Block_ptr->Grid.lfaceE(Ri, Rj)*Soln_Block_ptr->SpE[Ri][Rj];
    Ur = Soln_Block_ptr->Uo[ii][jj];
    Ul = Soln_Block_ptr->Uo[Ri][Rj];
    break;

    // WEST
  case WEST:
    Ri = ii+1; Rj=jj;
    nface = Soln_Block_ptr->Grid.nfaceW(Ri, Rj);
    lface = Soln_Block_ptr->Grid.lfaceW(Ri, Rj)*Soln_Block_ptr->SpW[Ri][Rj];
    Ur = Soln_Block_ptr->Uo[ii][jj];
    Ul = Soln_Block_ptr->Uo[Ri][Rj];
    break;

  } /* endswitch */ 


  // compute the derivatives
  for(int jcol=0; jcol<(NUM_VAR_RTE2D); jcol++){
    UA = Ur;
    UB = Ur;
    UA[jcol+1] += perturb; 	 
    UB[jcol+1] -= perturb; 
    
    FluxA = Flux_n(Ul,UA, nface);       
    FluxB = Flux_n(Ul,UB, nface);
    
    for(int irow=0; irow<(NUM_VAR_RTE2D); irow++){
      dFidU(irow,jcol) = (FluxA[irow+1] - FluxB[irow+1])/(TWO*perturb);           
    }
  } 

  // copy over

  dRdU += lface*dFidU;

}

/*****************************************************************************
 *  Calculate First Order Local Jacobian Block(s) Coresponding to Cell(i,j)  *
 *  using Roe                                                                *
 *****************************************************************************/
template<> inline void Block_Preconditioner<Rte2D_State,
					    Rte2D_Quad_Block,
					    Rte2D_Input_Parameters>::
First_Order_Inviscid_Jacobian_Roe(const int &cell_index_i,const int &cell_index_j, 
				   DenseMatrix* Jacobian){              
    
  //! Calculate Jacobian matrix -> blocksizexblocksize matrix in DenseMatrix format
  Preconditioner_dFIdU_Roe(Jacobian[NORTH],cell_index_i,cell_index_j,NORTH);
  Preconditioner_dFIdU_Roe(Jacobian[SOUTH],cell_index_i,cell_index_j,SOUTH); 
  Preconditioner_dFIdU_Roe(Jacobian[EAST],cell_index_i,cell_index_j,EAST);        
  Preconditioner_dFIdU_Roe(Jacobian[WEST],cell_index_i,cell_index_j,WEST); 

  //Center calculated from neighbours
  //! Using the fact that dF/dU(right) = - dF/dU(left) 
  Jacobian[CENTER] += (Jacobian[NORTH] + Jacobian[SOUTH] + Jacobian[EAST]  + Jacobian[WEST])
    /(Soln_Block_ptr->Grid.Cell[cell_index_i][cell_index_j].A * Soln_Block_ptr->Sp[cell_index_i][cell_index_j]);

  Jacobian[NORTH] = -Jacobian[NORTH]/(Soln_Block_ptr->Grid.Cell[cell_index_i][cell_index_j-1].A*Soln_Block_ptr->Sp[cell_index_i][cell_index_j-1]);
  Jacobian[SOUTH] = -Jacobian[SOUTH]/(Soln_Block_ptr->Grid.Cell[cell_index_i][cell_index_j+1].A*Soln_Block_ptr->Sp[cell_index_i][cell_index_j+1]);
  Jacobian[EAST] = -Jacobian[EAST]/(Soln_Block_ptr->Grid.Cell[cell_index_i-1][cell_index_j].A*Soln_Block_ptr->Sp[cell_index_i-1][cell_index_j]);
  Jacobian[WEST] = -Jacobian[WEST]/(Soln_Block_ptr->Grid.Cell[cell_index_i+1][cell_index_j].A*Soln_Block_ptr->Sp[cell_index_i+1][cell_index_j]);

}



#endif // _RTE2D_IMPLICIT_SPECIALIZATION_INCLUDED 

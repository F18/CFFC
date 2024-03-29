#ifndef _HEXA_EXPLICIT_SOLVER
#define _HEXA_EXPLICIT_SOLVER

/*! ******************************************************
 * Routine: Hexa_Explicit_Solver                        *
 *                                                      *
 ********************************************************/
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
int Hexa_MultiStage_Explicit_Solver(HexaSolver_Data &Data,
				    HexaSolver_Solution_Data<SOLN_pSTATE, SOLN_cSTATE> &Solution_Data) {
  
  Explicit_Filter_Commands::Initialize_Filters(Data,Solution_Data);


  int error_flag(0);
   
  double dTime;
  double residual_l2norm_first, residual_ratio;
  double residual_l1_norm, residual_l2_norm, residual_max_norm;

  int first_step(1);
  int limiter_freezing_off(ON);
  
  CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
  Data.processor_cpu_time.reset();
      
  /*! ***************** Explicit Iterations ******************
   * Perform required number of iterations (time steps).    *
   **********************************************************/
  
  if ((!Solution_Data.Input.Time_Accurate && 
       Solution_Data.Input.Maximum_Number_of_Time_Steps > 0 &&
       Data.number_of_explicit_time_steps < Solution_Data.Input.Maximum_Number_of_Time_Steps) ||
      (Solution_Data.Input.Time_Accurate && 
       Solution_Data.Input.Time_Max > Data.Time)) {
    
    if (!Data.batch_flag) cout << "\n\n Beginning explicit time-marching computations on " << Date_And_Time() << ".\n\n";
    
    while (1){ // Data.number_of_time_steps <= Solution_Data.Input.Maximum_Number_of_Time_Steps) {  // instead of while(1) ??
      
      /*! ***************** Time Step *******************/
      /* Determine local and global time steps. */
      dTime = Solution_Data.Local_Solution_Blocks.CFL(Solution_Data.Input);
      dTime = CFFC_Minimum_MPI(dTime); 

      // Find global minimum time step for all processors for Time Accurate      
      if (Solution_Data.Input.Time_Accurate) {
	
	if ((Solution_Data.Input.i_Time_Integration != TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING) &&
	    (Data.Time + Solution_Data.Input.CFL_Number*dTime > Solution_Data.Input.Time_Max)) {
	  dTime = (Solution_Data.Input.Time_Max-Data.Time)/Solution_Data.Input.CFL_Number;
	} else if (Data.Time + Solution_Data.Input.CFL_Number*dTime*
		   MultiStage_Optimally_Smoothing(Solution_Data.Input.N_Stage, 
						  Solution_Data.Input.N_Stage,
						  Solution_Data.Input.i_Limiter) > Solution_Data.Input.Time_Max) {
	  dTime = (Solution_Data.Input.Time_Max-Data.Time)/(Solution_Data.Input.CFL_Number*
                  MultiStage_Optimally_Smoothing(Solution_Data.Input.N_Stage, 
                                                 Solution_Data.Input.N_Stage,
                                                 Solution_Data.Input.i_Limiter));
	} /* endif */
	dTime = CFFC_Minimum_MPI(dTime);   
      } /* endif */
      
      if (!Solution_Data.Input.Local_Time_Stepping) { 
	Solution_Data.Local_Solution_Blocks.Set_Global_TimeStep(dTime);
      } /* endif */
      /**************************************************/
  
      /**************************************************/
      /* Update CPU time used for the calculation so far. */
      Data.processor_cpu_time.update();
      Data.total_cpu_time.cput = CFFC_Summation_MPI(Data.processor_cpu_time.cput); 
      /**************************************************/

      /***************** RESTART **************************/            
      /* Periodically save restart solution files. */
      if (!first_step && 
	  Data.number_of_explicit_time_steps-Solution_Data.Input.Restart_Solution_Save_Frequency*
	  (Data.number_of_explicit_time_steps/Solution_Data.Input.Restart_Solution_Save_Frequency) == 0 ) {
	if (!Data.batch_flag) {
	  cout << "\n\n  Saving solution to restart data file(s) after"
	       << " n = " << Data.number_of_explicit_time_steps 
	       << " steps (iterations).";
	  cout.flush();
	} /* endif */
	error_flag = Write_Octree(Data.Octree, 
                                  Solution_Data.Input);
	if (error_flag) {
	  cout << "\n ERROR: Unable to open octree data file "
	       << "on processor "<< CFFC_MPI::This_Processor_Number
	       << ".\n"; cout.flush();
	} /* endif */
	error_flag = CFFC_OR_MPI(error_flag);
	if (error_flag) return (error_flag);	

	error_flag = Solution_Data.Local_Solution_Blocks.Write_Restart_Solution(Solution_Data.Input,
										Data.Local_Adaptive_Block_List,
										Data.number_of_explicit_time_steps,
										Data.Time,
										Data.processor_cpu_time);
	if (error_flag) {
	  cout << "\n ERROR: Unable to open restart output data file(s) on processor "
	       << CFFC_MPI::This_Processor_Number
	       << ".\n";
	  cout.flush();
	} /* endif */
	error_flag = CFFC_OR_MPI(error_flag);
	if (error_flag) return (error_flag);
      } /* endif */     
      /**************************************************/	           
      
      /***************** NORMS **************************/      
      /* Determine the L1, L2, and max norms of the solution residual. */
      residual_l1_norm = Solution_Data.Local_Solution_Blocks.L1_Norm_Residual(Solution_Data.Input.Residual_Norm);
      residual_l1_norm = CFFC_Summation_MPI(residual_l1_norm); 
      residual_l2_norm = Solution_Data.Local_Solution_Blocks.L2_Norm_Residual(Solution_Data.Input.Residual_Norm);
      residual_l2_norm = CFFC_Summation_MPI(residual_l2_norm); 
      residual_max_norm = Solution_Data.Local_Solution_Blocks.Max_Norm_Residual(Solution_Data.Input.Residual_Norm);
      residual_max_norm = CFFC_Maximum_MPI(residual_max_norm);

      /* Output progress information for the calculation. */
      if (!Data.batch_flag) Output_Progress_L2norm(Data.number_of_explicit_time_steps,
						   Data.Time*THOUSAND,
						   Data.total_cpu_time,
						   residual_l2_norm,
						   first_step,
						   50);
      if (CFFC_Primary_MPI_Processor() && !first_step) {
	Output_Progress_to_File(Data.residual_file,
				Data.number_of_explicit_time_steps,
				Data.Time*THOUSAND,
				Data.total_cpu_time,
				residual_l1_norm,
				residual_l2_norm,
				residual_max_norm);
      } /* endif */

      if (!first_step) {
         error_flag = Output_Other_Solution_Progress_Specialization_Data(Data,
                                                                         Solution_Data);
         if (error_flag) return (error_flag);
      } /* endif */
      /****************************************************/

      /**************************************************
       * Check to see if calculations are complete.     *
       **************************************************/
      if (!Solution_Data.Input.Time_Accurate && 
          Data.number_of_explicit_time_steps >= Solution_Data.Input.Maximum_Number_of_Time_Steps) break;
      if (Solution_Data.Input.Time_Accurate && 
          Data.Time >= Solution_Data.Input.Time_Max) break;
      /**************************************************/

      /**************************************************
       * Freeze limiters as necessary.                  *
       **************************************************/
      if (!first_step && 
          Solution_Data.Input.Freeze_Limiter && 
          limiter_freezing_off &&           
	  residual_l2_norm <= Solution_Data.Input.Freeze_Limiter_Residual_Level) {
	//Freeze_Limiters(Local_SolnBlk, 
        //                List_of_Solution_Data.Local_Solution_Blocks);
	limiter_freezing_off = ON;	
      } /* endif */
      /**************************************************/   

        
    /*************************** SOLUTION FILTERING *****************************/            
    /* Periodically filter the solution to eliminate build-up of high frequency */
    if (Solution_Data.Input.ExplicitFilters_IP.Solution_Filtering_Frequency != OFF) {
                        
        if (!first_step && 
            Data.number_of_explicit_time_steps-Solution_Data.Input.ExplicitFilters_IP.Solution_Filtering_Frequency*
            (Data.number_of_explicit_time_steps/Solution_Data.Input.ExplicitFilters_IP.Solution_Filtering_Frequency) == 0 ) {
            if (!Data.batch_flag) {
                cout << "o";   // This symbol will notify a solution filtering operation
                cout.flush();
            } /* endif */
            
            error_flag = Explicit_Filter_Commands::Filter_Solution(Data,Solution_Data,Explicit_Filter_Constants::SECONDARY_FILTER);
            if (error_flag) {
                cout << "\n ERROR: Could not filter solution "
                << "on processor "
                << CFFC_MPI::This_Processor_Number
                << ".\n";
                cout.flush();
            } /* endif */
        }
    }
        
        
      /*********** BLOCK SOLUTION UPDATE *************************
       * Update solution for next time step using a multistage   *
       * time stepping scheme.                                   *
       ***********************************************************/
      for (int i_stage  = 1 ; i_stage <= Solution_Data.Input.N_Stage ; ++i_stage ) {
      
	/*******************************************************************/
	// 1. Send/Copy the information between blocks...
	error_flag = Send_Messages<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> >
                        (Solution_Data.Local_Solution_Blocks.Soln_Blks,
			 Data.Local_Adaptive_Block_List);
	if (error_flag) {
           cout << "\n ERROR: Message passing error during explicit time marching "
                << "on processor "
                << CFFC_MPI::This_Processor_Number
                << ".\n";
           cout.flush();
	} /* endif */
	error_flag = CFFC_OR_MPI(error_flag);        
	if (error_flag) return (error_flag);

	/************* BOUNDARY CONDITIONS *********************************/
	// 2. Apply boundary conditions for stage.
	Solution_Data.Local_Solution_Blocks.BCs(Solution_Data.Input);
	
	error_flag = Hexa_BCs_Specializations(Data, Solution_Data);
	if (error_flag) return (error_flag);
	 /*************** UPDATE CORNER GHOST CELLS SOLUTION *****/
        // For those 3 blocks abutting cases, this is to fix the gradient-reconstruction.
        Solution_Data.Local_Solution_Blocks.Update_Corner_Cells_for_3_Blks_Abutting(Data.Local_Adaptive_Block_List);

	 /*************** UPDATE SOLUTION **********************************/
	// 3. Determine solution residuals for stage.
	error_flag = Solution_Data.Local_Solution_Blocks.dUdt_Multistage_Explicit(Solution_Data.Input,
                                                                                  i_stage);        
	if (error_flag) {
	  cout << "\n ERROR: Solution residual error on processor"
	       << CFFC_MPI::This_Processor_Number
               << ".\n";
	  cout.flush();
	} /* endif */
	error_flag = CFFC_OR_MPI(error_flag);        
	if (error_flag) return (error_flag);
          	
	/*******************************************************************/
	// 4. Send boundary flux corrections at block interfaces with resolution changes.

	/*******************************************************************/
	// 5. Apply boundary flux corrections to ensure that method is conservative.

	/*******************************************************************/
	// 6. Smooth the solution residual using implicit residual smoothing. */

	/*******************************************************************/
    // 7. Explicit filtering of the solution residual.      

    Explicit_Filter_Commands::Filter_Residual(i_stage,Data,Solution_Data);

    /*******************************************************************/
	// 8. Update solution for stage.
	error_flag = 
           Solution_Data.Local_Solution_Blocks.Update_Solution_Multistage_Explicit(Solution_Data.Input, 
                                                                                           i_stage);
	if (error_flag) {	  
	  cout << "\n ERROR: Solution update error on processor "
	       << CFFC_MPI::This_Processor_Number
               << ".\n";
	  cout.flush();
	} /* endif */
	error_flag = CFFC_OR_MPI(error_flag);
	if (error_flag) return (error_flag);
                    
          
    /*******************************************************************/
    // 9. Filter Solution in case of Filtering the variables instead of Residuals        
    if (Solution_Data.Input.ExplicitFilters_IP.Filter_Method == Explicit_Filter_Constants::FILTER_VARIABLES) {
        error_flag = Explicit_Filter_Commands::Filter_Solution(Data,Solution_Data,Explicit_Filter_Constants::PRIMARY_FILTER);
        if (error_flag) {
            cout << "\n ERROR: Could not filter solution "
            << "on processor "
            << CFFC_MPI::This_Processor_Number
            << ".\n";
            cout.flush();
        } /* endif */
    }
        
          
      }  // END Multistage for loop
      /************************************************************************/

        
      /******************* UPDATE TIMER & COUNTER *****************************
       *    Update time and time step counter.                                *
       ************************************************************************/
      if (first_step) first_step = 0;
      Data.number_of_explicit_time_steps++;
      if (Solution_Data.Input.i_Time_Integration != TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING) {	
	Data.Time = Data.Time + 
                    Solution_Data.Input.CFL_Number*dTime;
      } else {
	Data.Time = Data.Time + 
                    Solution_Data.Input.CFL_Number*dTime*
                    MultiStage_Optimally_Smoothing(Solution_Data.Input.N_Stage, 
			 		           Solution_Data.Input.N_Stage,
						   Solution_Data.Input.i_Limiter);	
      } /* endif */
       
    } // END WHILE(1) LOOP 
      
      if (!Data.batch_flag) {
          cout 
          << "\n\n Explicit time-marching computations complete on " 
          << Date_And_Time() << " after"
          << " n = " << Data.number_of_explicit_time_steps 
          << " steps (iterations).";
          cout.flush();   
      }
    
  } // END ( Time or Steps) IF

  CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
 
  /************************************************************************************  
   END OF EXPLICIT TIME MARCHING SCHEME CALCULATION:
   Update ghostcell information and prescribe boundary conditions to ensure
   that the solution is consistent on each block. 
   *************************************************************************************/
  error_flag = Send_Messages<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> >
                  (Solution_Data.Local_Solution_Blocks.Soln_Blks,
                   Data.Local_Adaptive_Block_List);
  if (error_flag) {
     cout << "\n ERROR: Message passing error after completion of explicit time marching "
          << "on processor "
          << CFFC_MPI::This_Processor_Number
          << ".\n";
     cout.flush();
  } /* endif */
  error_flag = CFFC_OR_MPI(error_flag);        
  if (error_flag) return (error_flag);

  Solution_Data.Local_Solution_Blocks.BCs(Solution_Data.Input);

  error_flag = Hexa_BCs_Specializations(Data, Solution_Data);
  if (error_flag) return (error_flag);

  /*************** UPDATE CORNER GHOST CELLS SOLUTION *****/
  // For those 3 blocks abutting cases, this is to fix the gradient-reconstruction.
  Solution_Data.Local_Solution_Blocks.Update_Corner_Cells_for_3_Blks_Abutting(Data.Local_Adaptive_Block_List);

  return error_flag;

}

#endif //_HEXA_EXPLICIT_SOLVER

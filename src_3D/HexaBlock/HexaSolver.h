#ifndef _HEXA_SOLVER_INCLUDED
#define _HEXA_SOLVER_INCLUDED

/* Include required CFFC header files. */

#ifndef _HEXA_MULTIBLOCK_INCLUDED
#include "HexaMultiBlock.h"
#endif //_HEXA_MULTIBLOCK_INCLUDED

#ifndef _AMR_INCLUDED
#include "../AMR/AMR.h"
#endif // _AMR_INCLUDED

#ifndef _ADAPTIVEBLOCK3D_MESSAGEPASSING_INCLUDED
#include "../AMR/AdaptiveBlock3D_MessagePassing.h"
#endif // _ADAPTIVEBLOCK3D_MESSAGEPASSING_INCLUDED

/********************************************************
 * Routine: HexaSolver                                  *
 *                                                      *
 * Computes solutions to the governing PDEs on          *
 * multi-block AMR mesh composed of body-fitted         *
 * hexahedral solution blocks.                          *
 *                                                      *
 ********************************************************/

template<class SOLN_pSTATE, class SOLN_cSTATE>
int HexaSolver(char *Input_File_Name_ptr,
               int batch_flag){
   
/********************************************************  
 * Local variable declarations.                         *
 ********************************************************/
   
  /* Define the primary solution variables for the solution
     of the problem using a multi-block hexahedral mesh with
     adaptive mesh refinement (AMR). */ 

   Input_Parameters<SOLN_pSTATE, SOLN_cSTATE>              Input;
   Grid3D_Hexa_Multi_Block_List                            Initial_Mesh;
   Hexa_Multi_Block<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> > Local_Solution_Blocks; 
    
   AdaptiveBlock3D_List                                    Local_Adaptive_Block_List; 
   AdaptiveBlock3D_ResourceList                            Global_Adaptive_Block_List;
   Octree_DataStructure                                    Octree;

   /* Define residual file and cpu time variables. */

   ofstream residual_file;
   CPUTime processor_cpu_time, total_cpu_time;
   
   /* Define other local variables. */

   int number_of_time_steps, first_step,
       command_flag, error_flag, line_number, 
       i_stage, perform_explicit_time_marching, limiter_freezing_off;
   
   double Time, dTime;
   double residual_l2norm_first, residual_ratio;
   double residual_l1_norm, residual_l2_norm, residual_max_norm;
 
   /********************************************************  
    * Set default values for the input solution parameters *
    * and then read user specified input values from the   *
    * specified input parameter file.                      *
    ********************************************************/

   // The primary MPI processor processes the input parameter file.
   if (CFFC_Primary_MPI_Processor()) {
      if (!batch_flag) {
         cout << "\n Reading input data file `"
              << Input_File_Name_ptr << "'.";
         cout.flush();
      } /* endif */
   
      error_flag = Input.Process_Input_Control_Parameter_File(Input_File_Name_ptr,
                                                              command_flag);
   } else {
      error_flag = 0;
   } /* endif */

   // Broadcast input solution parameters to other MPI processors.
   CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
   CFFC_Broadcast_MPI(&error_flag, 1);
   if (error_flag != 0) return (error_flag);
   CFFC_Broadcast_MPI(&command_flag, 1);
   if (command_flag == TERMINATE_CODE) return (0);
   Input.Broadcast();

   /********************************************************  
    * Create initial mesh and allocate solution            *
    * variables for specified IBVP/BVP problem.            *
    ********************************************************/
 
   execute_new_calculation: ;
   CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.

   /* Create initial mesh.  Read mesh from grid definition or data files 
      as specified by input parameters. */

   // The primary MPI processor creates the initial mesh.
   if (CFFC_Primary_MPI_Processor()) {
      Initial_Mesh.Create_Grid(Input.Grid_IP);
      error_flag = 0;
      //Outputting solution input parameters
      if (!batch_flag && error_flag == 0) {
         cout << Input << "\n";
         cout.flush();
      } /* endif */
   } /* endif */

   CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
   CFFC_Broadcast_MPI(&error_flag, 1); // Broadcast mesh creation error flag.
   if (error_flag) return (error_flag);

   //Broadcast the mesh to other MPI processors.
   Initial_Mesh.Broadcast();

   /* Create (allocate) list of hexahedral solution blocks on each 
      processor. */

   if (!batch_flag) {
      cout << "\n Creating multi-block octree data structure and assigning"
           << "\n  solution blocks corresponding to initial mesh.";
      cout.flush();
   } /* endif */

/* // ************************************************************** */
/* // output the grid geometry... */
/*    char prefix[256], extension[256], output_file_name[256]; */
/*    char *output_file_name_ptr; */
/*    ofstream output_file; */
   
/*    /\* Determine prefix of output data file names. *\/ */
/*    int i = 0; */
   
/*    while (1) { */
/*       if (Input.Output_File_Name[i] == ' ' || */
/*           Input.Output_File_Name[i] == '.') break; */
/*       prefix[i]=Input.Output_File_Name[i]; */
/*       i = i + 1; */
/*       if (i > strlen(Input.Output_File_Name) ) break; */
/*    } /\* endwhile *\/ */
/*    prefix[i] = '\0'; */
/*    strcat(prefix, "_cells_cpu"); */
   
/*    sprintf(extension, "%.6d"); */
/*    strcat(extension, ".dat"); */
/*    strcpy(output_file_name, prefix); */
/*    strcat(output_file_name, extension); */
/*    output_file_name_ptr = output_file_name; */
       
/*    // Open the output data file. */
/*    output_file.open(output_file_name_ptr, ios::out); */
/*    Initial_Mesh.Output_Cells_Tecplot(output_file); */
/*    output_file.close(); */
   
/*    for(int ii= -1; ii<2; ii++) */
/*       for(int jj= -1; jj<2; jj++) */
/*          for(int kk= -1; kk<2; kk++){ */
/*             Initial_Mesh.Connectivity[1].neighT_info.compute_message_tag(ii, jj, kk); */
            
/*          } */

   
   Create_Initial_Solution_Blocks<SOLN_pSTATE, SOLN_cSTATE>(Initial_Mesh,
                                                            Local_Solution_Blocks,
                                                            Input,
                                                            Octree,
                                                            Global_Adaptive_Block_List,
                                                            Local_Adaptive_Block_List);
   
   /********************************************************  
    * Initialize solution variables.                       *
    ********************************************************/

   /* Set the initial time level. */
   Time = ZERO;
   number_of_time_steps = 0;

   /* Set the CPU time to zero. */
   processor_cpu_time.zero();
   total_cpu_time.zero();

   /* Initialize the conserved and primitive state
      solution variables. */
         
   if (!batch_flag) {
      cout << "\n Prescribing initial data.";
      cout.flush();
   } /* endif */

   if (Input.i_ICs == IC_RESTART) {

      error_flag = Wall_Distance(Local_Solution_Blocks.Soln_Blks,
                                 Octree, 
                                 Local_Adaptive_Block_List);
      if (error_flag) {
         cout << "\n  ERROR: Difficulty determining the wall distance "
              << "on processor "<< CFFC_MPI::This_Processor_Number
              << ".\n";
         cout.flush();
      } /* endif */
      if (!batch_flag){ 
         cout << "\n Reading solution from restart data files."; 
         cout.flush();
      }
      //error_flag = Read_Octree(Octree,
      //			 Input);
      if (error_flag) {
         cout << "\n ERROR: Unable to open Octree data file "
 	      << "on processor "
  	      << Local_Adaptive_Block_List.ThisCPU
 	      << ".\n";
         cout.flush();
      } // endif 
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return (error_flag);

      error_flag = Local_Solution_Blocks.Read_Restart_Solution(Input,  
                                                               Local_Adaptive_Block_List,
                                                               number_of_time_steps,
                                                               Time,
                                                               processor_cpu_time);
      if (error_flag) {
         cout << "\n  ERROR: Unable to open restart input data file(s) "
              << "on processor "<< CFFC_MPI::This_Processor_Number
              << ".\n";
         cout.flush();
      } /* endif */
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return (error_flag);

      // Ensure each processor has the correct time and time!!!
      number_of_time_steps = CFFC_Maximum_MPI(number_of_time_steps);
      Time = CFFC_Maximum_MPI(Time);
      processor_cpu_time.cput = CFFC_Maximum_MPI(processor_cpu_time.cput);
      Input.Maximum_Number_of_Time_Steps = CFFC_Maximum_MPI(Input.Maximum_Number_of_Time_Steps);
      
   } else {
      error_flag = Wall_Distance(Local_Solution_Blocks.Soln_Blks,
                                 Octree, 
                                 Local_Adaptive_Block_List);
      if (error_flag) {
         cout << "\n  ERROR: Difficulty determining the wall distance "
              << "on processor "<< CFFC_MPI::This_Processor_Number
              << ".\n";
         cout.flush();
      } /* endif */
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return (error_flag);

      Local_Solution_Blocks.ICs(Input);
   } /* endif */

   /* Send solution information between neighbouring blocks to complete
     prescription of initial data. */

   Send_All_Messages<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> >(Local_Solution_Blocks.Soln_Blks,
                                                            Local_Adaptive_Block_List,
                                                            Local_Solution_Blocks.Soln_Blks[0].NumVar(),
                                                            OFF);
   exit(1);
   /* Prescribe boundary data consistent with initial data. */

 
   

   Local_Solution_Blocks.BCs(Input);

   /********************************************************  
    * Solve IBVP or BVP for conservation form of PDEs      *
    * on multi-block body-fitted AMR hexadedral mesh.      *
    ********************************************************/

   continue_existing_calculation: ;
   CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
 
   /* Open residual file and reset the CPU time. */

   first_step = 1;
   limiter_freezing_off = ON;
   
   if (CFFC_Primary_MPI_Processor()) {
      error_flag = Open_Progress_File(residual_file,
                                      Input.Output_File_Name,
                                      number_of_time_steps);
      if (error_flag) {
         cout << "\n ERROR: Unable to open residual file for the calculation.\n";
         cout.flush();
      } /* endif */
   } /* endif */
   CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
   CFFC_Broadcast_MPI(&error_flag, 1);
   if (error_flag) return (error_flag);

   processor_cpu_time.reset();

   /* Perform required number of iterations (time steps). */

   if ((!Input.Time_Accurate &&
        Input.Maximum_Number_of_Time_Steps > 0 &&
        number_of_time_steps < Input.Maximum_Number_of_Time_Steps) ||
       (Input.Time_Accurate &&
        Input.Time_Max > Time)) {
     
      if (!batch_flag) {
         cout << "\n\n Beginning computations on " 
              << Date_And_Time() 
              << ".\n\n";
         cout.flush();
      } /* endif */

      if ((!Input.Time_Accurate &&
           Input.Maximum_Number_of_Time_Steps > 0 &&
           number_of_time_steps < Input.Maximum_Number_of_Time_Steps) ||
          (Input.Time_Accurate &&
           Input.Time_Max > Time)) {
         perform_explicit_time_marching = ON;
      } else {
         perform_explicit_time_marching = OFF;
      } /* endif */
   
      while (perform_explicit_time_marching) {
          /* Determine local and global time steps. */
          dTime =  Local_Solution_Blocks.CFL(Input);
          dTime = CFFC_Minimum_MPI(dTime); // Find global minimum time step for all processors.

          if (Input.Time_Accurate) {
             if ((Input.i_Time_Integration != TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING) &&
                 (Time + Input.CFL_Number*dTime > Input.Time_Max)) {
                dTime = (Input.Time_Max-Time)/Input.CFL_Number;
             } else if (Time + Input.CFL_Number*dTime*
                        MultiStage_Optimally_Smoothing(Input.N_Stage, 
                                                       Input.N_Stage,
                                                       Input.i_Limiter) > Input.Time_Max) {
                dTime = (Input.Time_Max-Time)/
                        (Input.CFL_Number*
                        MultiStage_Optimally_Smoothing(Input.N_Stage, 
                                                       Input.N_Stage,
                                                       Input.i_Limiter));
             } /* endif */
          } /* endif */

          // Set global time step.
          if (!Input.Local_Time_Stepping) { 
             Local_Solution_Blocks.Set_Global_TimeStep(dTime);
          } /* endif */
       
          /* Determine the L1, L2, and max norms of the solution residual. */

          // L1 norm for all processors.
          residual_l1_norm =  Local_Solution_Blocks.L1_Norm_Residual();
          residual_l1_norm = CFFC_Summation_MPI(residual_l1_norm); 

          // L2 norm for all processors.
          residual_l2_norm =  Local_Solution_Blocks.L2_Norm_Residual();
  	  residual_l2_norm = sqr(residual_l2_norm);
          residual_l2_norm = CFFC_Summation_MPI(residual_l2_norm); 
          residual_l2_norm = sqrt(residual_l2_norm);
          
 	  // Max norm for all processors.  
  	  residual_max_norm =  Local_Solution_Blocks.Max_Norm_Residual();
 	  residual_max_norm = CFFC_Maximum_MPI(residual_max_norm);
 	  
          /* Update CPU time used for the calculation so far. */
 	  processor_cpu_time.update();
          total_cpu_time.cput = CFFC_Summation_MPI(processor_cpu_time.cput); // Total CPU time for all processors.
 	
          /* Periodically save restart solution files. */
          if (!first_step &&
              number_of_time_steps-Input.Restart_Solution_Save_Frequency*
              (number_of_time_steps/Input.Restart_Solution_Save_Frequency) == 0 ) {
	     if (!batch_flag) {
               cout << "\n\n  Saving solution to restart data file(s) after"
                    << " n = " << number_of_time_steps 
                    << " steps (iterations).";
               cout.flush();
	     } /* endif  */
             
             //error_flag = Write_Octree(Octree, Input);
             error_flag = 0;
             if (error_flag) {
                cout << "\n ERROR: Unable to open octree data file "
                     << "on processor "<< CFFC_MPI::This_Processor_Number
                     << ".\n";
                cout.flush();
             } /* endif */
             error_flag = CFFC_OR_MPI(error_flag);
             if (error_flag) return (error_flag);

             error_flag = Local_Solution_Blocks.Write_Restart_Solution(Input,
                                                                       Local_Adaptive_Block_List,
                                                                       number_of_time_steps,
                                                                       Time,
                                                                       processor_cpu_time);
             if (error_flag) {
                cout << "\n ERROR: Unable to open restart output data file(s) on processor "
                     << CFFC_MPI::This_Processor_Number
                     << ".\n";
                cout.flush();
             } /* endif */
             error_flag = CFFC_OR_MPI(error_flag);
             if (error_flag) return (error_flag);
  	  } /* endif */
       
 	  /* Output progress information for the calculation. */

          if (!batch_flag) {
             Output_Progress_L2norm(number_of_time_steps,
                                    Time*THOUSAND,
                                    total_cpu_time,
                                    residual_l2_norm,
                                    first_step,
                                    50);
          } /* endif */

          if (CFFC_Primary_MPI_Processor() && !first_step) {
             Output_Progress_to_File(residual_file,
                                     number_of_time_steps,
                                     Time*THOUSAND,
                                     total_cpu_time,
                                     residual_l1_norm,
                                     residual_l2_norm,
                                     residual_max_norm);
          } /* endif */

         /* Check to see if calculations are complete. */

         if (!Input.Time_Accurate && number_of_time_steps >= Input.Maximum_Number_of_Time_Steps) break;
         if (Input.Time_Accurate && Time >= Input.Time_Max) break;
         
         /* Freeze limiters as necessary. */

         if (!first_step && Input.Freeze_Limiter && limiter_freezing_off &&           
             residual_l2_norm <= Input.Freeze_Limiter_Residual_Level) {
            //Freeze_Limiters(Local_Solution_Blocks, 
            // 		      Local_Adaptive_Block_List);
            limiter_freezing_off = ON;	
         } /* endif */
       
         /* Update solution for next time step using a multistage
            time stepping scheme. */
         for ( i_stage  = 1 ; i_stage <= Input.N_Stage ; ++i_stage ) {
       
            // 1. Send/Copy the information between blocks...
            Send_All_Messages<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> >(Local_Solution_Blocks.Soln_Blks,
                                                                     Local_Adaptive_Block_List,
                                                                     Local_Solution_Blocks.Soln_Blks[0].NumVar(),
                                                                     OFF);
          
            // 2. Apply boundary conditions for stage.
            Local_Solution_Blocks.BCs(Input);

            // 3. Determine solution residuals for stage.
            error_flag =  Local_Solution_Blocks.dUdt_Multistage_Explicit(Input, i_stage);
            if (error_flag) {
               cout <<"\n ERROR: solution residual error on processor"
                    << CFFC_MPI::This_Processor_Number<< ".\n";
               cout.flush();
            } /* endif */
            error_flag = CFFC_OR_MPI(error_flag);
            if (error_flag) return (error_flag);
          
            // 7. Update solution for stage.
            error_flag =  Local_Solution_Blocks.Update_Solution_Multistage_Explicit(Input, i_stage);
            if (error_flag) {
               cout << "\n ERROR: solution update error on processor "
                    << CFFC_MPI::This_Processor_Number<<".\n";
               cout.flush();
            } /* endif */
            error_flag = CFFC_OR_MPI(error_flag);
            if (error_flag) return (error_flag);
          
         } /* endfor */
        
         /* Update time and time step counter. */
         if (first_step) first_step = 0;
         number_of_time_steps = number_of_time_steps + 1;
         if (Input.i_Time_Integration != TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING) {
             Time = Time + Input.CFL_Number*dTime;
         } else {
            Time = Time + Input.CFL_Number*dTime*MultiStage_Optimally_Smoothing(Input.N_Stage, 
                                                                                Input.N_Stage,
                                                                                Input.i_Limiter);
         } /* endif */
       
      } /* endwhile */
      
      if (!batch_flag) {
         cout << "\n\n Computations complete on " 
              << Date_And_Time() << ".\n";
         cout.flush();
      } /* endif */
      
   } /* endif */
  
   CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
    
   /* Update ghostcell information and prescribe boundary conditions to ensure
      that the solution is consistent on each block. */

   Send_All_Messages<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> >(Local_Solution_Blocks.Soln_Blks,
                                                            Local_Adaptive_Block_List,
                                                            Local_Solution_Blocks.Soln_Blks[0].NumVar(),
                                                            OFF);

   Local_Solution_Blocks.BCs(Input);

   /* Close residual file. */

   if (CFFC_Primary_MPI_Processor()) error_flag = Close_Progress_File(residual_file);

   /********************************************************
    * Solution calculations complete.                      *
    * Write current solution to output and restart files   *
    * as required, reset solution parameters, and run      *
    * other cases as specified by input parameters.        *
    ********************************************************/ 

   postprocess_current_calculation: ;
   CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.

   while (1) {
      if (CFFC_Primary_MPI_Processor()) {    
         Input.Get_Next_Input_Control_Parameter(true);
         command_flag = Input.Parse_Next_Input_Control_Parameter();
         line_number = Input.Line_Number;
      } /* endif */
      CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
      Input.Broadcast();
      CFFC_Broadcast_MPI(&command_flag, 1);

      if (command_flag == EXECUTE_CODE) {
         // Deallocate memory for exiting solution blocks.
	 if (!batch_flag) {
            cout <<"\n Deallocating existing solution blocks.";
            cout.flush();
         } /* endif */

         // Destructor of the class: Local_Block_List is automatically
         // called when the scope of the calculation is done...
         // ....

         // Output input parameters for new caluculation.
         if (!batch_flag)  {
            cout << "\n\n Starting a new calculation.";
            cout << Input << "\n";
            cout.flush();
         } /* endif */
         
         // Execute new calculation.
         goto execute_new_calculation;
         
      } else if (command_flag == TERMINATE_CODE) {
         // Deallocate memory for existing solution blocks.
         // Memory for the list of solution blocks is automatically
         // deallocated at the end execution.

         // Close input data file.         
	 if (!batch_flag) {
            cout <<"\n Closing input data file.";
            cout.flush();
         } /* endif */
         if (CFFC_Primary_MPI_Processor()) Input.Close_Input_File();

         // Terminate calculation
         return (0);
          
      } else if (command_flag == CONTINUE_CODE) {
        // Reset maximum time step counter.
        Input.Maximum_Number_of_Time_Steps += number_of_time_steps;
        
        // Output input parameters for continuing calculation.
        if (!batch_flag) {
	   cout << "\n\n Continuing existing calculation.";
	   cout << Input << "\n";
	   cout.flush();
        } /* endif */
      
        // Continue existing calculation.
        goto continue_existing_calculation;

      } else if (command_flag == WRITE_OUTPUT_CODE) {
         // Output solution data.
	 if (!batch_flag) {
           cout << "\n Writing solution to output data file(s).";
           cout.flush();
           cout.flush();
 	 } /* endif */

         error_flag = Local_Solution_Blocks.Output_Tecplot(Input, 
                                                           Local_Adaptive_Block_List,
                                                           number_of_time_steps,
                                                           Time);
         if (error_flag) {
            cout << "\n  ERROR: Unable to open  node output " 
               "data file(s) on processor "
             <<CFFC_MPI::This_Processor_Number<< ".\n";
            cout.flush();
         } /* endif */
         error_flag = CFFC_OR_MPI(error_flag);
         if (error_flag) return (error_flag);   
         
      } else if (command_flag == WRITE_OUTPUT_CELLS_CODE) {
         // Output solution data.
	 if (!batch_flag) {
           cout <<"\n Writing cell-centered solution to output data file(s).";
           cout.flush();
         } /* endif */
         
         error_flag = Local_Solution_Blocks.Output_Cells_Tecplot(Input, 
                                                                 Local_Adaptive_Block_List,
                                                                 number_of_time_steps,
                                                                 Time);
         if (error_flag) {
            cout << "\n  ERROR: Unable to open  cell output data file(s) on processor "
                 << CFFC_MPI::This_Processor_Number
                 << ".\n";
            cout.flush();
         } /* endif */   
         error_flag = CFFC_OR_MPI(error_flag);
         if (error_flag) return (error_flag);   
         
      } else if (command_flag == WRITE_OUTPUT_NODES_CODE) {
         // Output solution data.
	 if (!batch_flag) {
            cout<<"\n Writing nodal solution to output data file(s).";
            cout.flush();
         } /* endif */
         
         error_flag = Local_Solution_Blocks.Output_Nodes_Tecplot(Input, 
                                                                 Local_Adaptive_Block_List,
                                                                 number_of_time_steps,
                                                                 Time);
         if (error_flag) {
            cout << "\n  ERROR: Unable to open  cell output data file(s) on processor "
                 << CFFC_MPI::This_Processor_Number
                 << ".\n";
            cout.flush();
         } /* endif */   
         error_flag = CFFC_OR_MPI(error_flag);
         if (error_flag) return (error_flag);   

      } else if (command_flag == WRITE_RESTART_CODE) {
         // Write restart files.
	 if (!batch_flag) {
            cout <<"\n Writing solution to restart data file(s).";
	    cout.flush();
         } /* endif */

         error_flag = Write_Octree(Octree, Input);
         if (error_flag) {
            cout << "\n  ERROR: Unable to open  octree data file on processor " 
                 << CFFC_MPI::This_Processor_Number
                 << ".\n";
            cout.flush();
         } /* endif */
         error_flag = CFFC_OR_MPI(error_flag);
         if (error_flag) return (error_flag);
         
         error_flag = Local_Solution_Blocks.Write_Restart_Solution(Input, 
                                                                   Local_Adaptive_Block_List,
                                                                   number_of_time_steps,
                                                                   Time,
                                                                   processor_cpu_time);
         if (error_flag) {
            cout << "\n  ERROR: Unable to open restart output data file(s) on processor " 
                 << CFFC_MPI::This_Processor_Number
                 << ".\n";
            cout.flush();
         } /* endif */
         error_flag = CFFC_OR_MPI(error_flag);
         if (error_flag) return (error_flag);

      }else if (command_flag == READ_RESTART_CODE){
        //Read restart files.
        if (!batch_flag){ 
            cout << "\n Reading solution from restart data files."; 
            cout.flush();
        }
        //error_flag = Read_Octree(Octree,
   	//		           Input);
        if (error_flag) {
           cout << "\n ERROR: Unable to open Octree data file "
   	        << "on processor "
  	        << Local_Adaptive_Block_List.ThisCPU
 	        << ".\n";
           cout.flush();
        } // endif 
        error_flag = CFFC_OR_MPI(error_flag);
        if (error_flag) return (error_flag);

        error_flag = Local_Solution_Blocks.Read_Restart_Solution(Input,  
                                                                 Local_Adaptive_Block_List,
                                                                 number_of_time_steps,
                                                                 Time,
                                                                 processor_cpu_time);
        if (error_flag) {
           cout << "\n  ERROR: Unable to open restart input data file(s) on processor "
                << CFFC_MPI::This_Processor_Number << ".\n";
           cout.flush();
        } /* endif */
        error_flag = CFFC_OR_MPI(error_flag);
        if (error_flag) return (error_flag);

      } else if (command_flag == MORTON_ORDERING_CODE) {
         if (!batch_flag){ 
             cout << "\n\n Applying Morton re-ordering algorithm.";
             cout.flush();
         }
         error_flag = Morton_ReOrdering_of_Solution_Blocks(Octree,
				       			   Local_Adaptive_Block_List, 
                                                           Local_Solution_Blocks, 
                                                           Input, 
                                                           number_of_time_steps, 
                                                           Time, 
                                                           processor_cpu_time);
         if (error_flag) {
            cout <<"\n ERROR: Morton re-ordering error on processor "
                 << Local_Adaptive_Block_List.ThisCPU
                 << ".\n";
            cout.flush();
            return (error_flag);
         } // endif 
         error_flag = CFFC_OR_MPI(error_flag);
         if (error_flag) return (error_flag);

         //Output space filling curve in Tecplot format
         if (!batch_flag) cout << "\n Outputting space filling curve showing block loading for CPUs.";
         Morton_SFC_Output_Tecplot3D(Input,
                                     Local_Solution_Blocks,
				     Local_Adaptive_Block_List);
    
      } else if (command_flag == INVALID_INPUT_CODE ||
                 command_flag == INVALID_INPUT_VALUE) {
         line_number = -line_number;
         cout << "\n  ERROR: Error reading  data at line #"
              << -line_number  << " of input data file.\n";
         cout.flush();
         return (line_number);
      } /* endif */
         
   } /* endwhile */
 
   /********************************************************  
    * End of all solver computations and I/O.              *
    ********************************************************/
  
   return (0);
   
}

#endif // _HEXA_SOLVER_INCLUDED

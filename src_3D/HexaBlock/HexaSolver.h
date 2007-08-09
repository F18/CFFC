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
 * Computes solutions to 3D Euler equations on 3D       *
 * Hexaheral block mesh                                 *
 *                                                      *
 ********************************************************/

template<class SOLN_pSTATE, class SOLN_cSTATE>
   int HexaSolver(char *Input_File_Name_ptr,
               int batch_flag){
   
/********************************************************  
 * Local variable declarations.                         *
 ********************************************************/
   /* Define residual file and cpu time variables. */
   ofstream residual_file;
   CPUTime processor_cpu_time, total_cpu_time;
   
   /* Other local solution variables. */
   int number_of_time_steps, first_step,
      command_flag, error_flag, line_number, 
      i_stage, perform_explicit_time_marching, limiter_freezing_off;
   
   double Time, dTime;
   double residual_l2norm_first, residual_ratio;
   double residual_l1_norm, residual_l2_norm, residual_max_norm;
   
 
   Input_Parameters<SOLN_pSTATE, SOLN_cSTATE> IPs;
   Grid3D_Hexa_Multi_Block        Initial_MultiBlock_Grid;
    
   AdaptiveBlock3DResourceList    List_of_Global_Solution_Blocks;
   AdaptiveBlock3D_List           List_of_Local_Solution_Blocks; 
   OcTreeBlock_DataStructure      OcTree;
      
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
      } /* endif */
   
      error_flag = IPs.Process_Input_Control_Parameter_File
                   (Input_File_Name_ptr,
                    command_flag);

      if (!batch_flag && error_flag == 0) {
         cout << IPs << "\n";
         cout.flush();
      } /* endif */
   } else {
      error_flag = 0;
   } /* endif */

   // Broadcast input solution parameters to other MPI processors.
   CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
   CFFC_Broadcast_MPI(&error_flag, 1);
   if (error_flag != 0) return (error_flag);
   CFFC_Broadcast_MPI(&command_flag, 1);
   if (command_flag == TERMINATE_CODE) return (0);
   IPs.Broadcast_Input_Parameters();

   /********************************************************  
    * Create initial mesh and allocate solution            *
    * variables for specified IBVP/BVP problem.            *
    ********************************************************/
   execute_new_calculation: ;
   // create grids on each processor and then 
   // delete the unused grid on each processor 
   
   Initial_MultiBlock_Grid.Create_Grid(IPs.IP_Grid);
   
   // create the solution block list on each processor
   Hexa_MultiBlock< Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> > Local_SolnBlk(IPs);
      
   if (!batch_flag) cout << "\n Creating multi-block octree data structure and assigning"
                         << "\n  solution blocks corresponding to initial mesh.";
   
   Create_Initial_Solution_Blocks<SOLN_pSTATE, SOLN_cSTATE>(
      Initial_MultiBlock_Grid,
      Local_SolnBlk,
      IPs,
      OcTree,
      List_of_Global_Solution_Blocks,
      List_of_Local_Solution_Blocks);
     
   // The Grid3D_Hexa_Multi_Block class is no longer required.  Destruction of the
   // class will only delete unused grid blocks (used grid blocks have been handed
   // off to the local solution blocks who are now responsible for deleting that
   // memory).
   //    delete Create_MultiBlock_Grid; Create_MultiBlock_Grid = NULL;
  
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
         
   if (!batch_flag) cout << "\n Prescribing initial data.";
   if (IPs.i_ICs == IC_RESTART) {
      
      Local_SolnBlk.Create_Wall_Data( );
      
      error_flag = Local_SolnBlk.Read_Restart_Solution
         (IPs,  List_of_Local_Solution_Blocks,
          number_of_time_steps,
          Time,
          processor_cpu_time);
      
      if (error_flag) {
         cout << "\n  ERROR: Unable to open restart" 
            " input data file(s) "
              << "on processor "<< CFFC_MPI::This_Processor_Number
              << ".\n";
         cout.flush();
      } /* endif */

      if (error_flag) return (error_flag);
      // Ensure each processor has the correct time and time!!!
      number_of_time_steps = CFFC_Maximum_MPI(number_of_time_steps);
      Time = CFFC_Maximum_MPI(Time);
      processor_cpu_time.cput = CFFC_Maximum_MPI(processor_cpu_time.cput);
      IPs.Maximum_Number_of_Time_Steps =
         CFFC_Maximum_MPI(IPs.Maximum_Number_of_Time_Steps);
      
   } else {
      
      Local_SolnBlk.Create_Wall_Data( );

      Local_SolnBlk.ICs(IPs);
      
   } /* endif */

   Send_All_Messages<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> >(Local_SolnBlk.Hexa_Block_List,
                           List_of_Local_Solution_Blocks,
                           Local_SolnBlk.Hexa_Block_List[0]->NumVar(),
                           OFF);

//   Local_SolnBlk.Copy_MultiBlk_Boundary_Info( List_of_Local_Solution_Blocks,IPs.IP_Grid);
   
   
   /* Prescribe boundary data consistent with initial data. */

    Local_SolnBlk.BCs(IPs);

   
   //    /********************************************************  
//    * Solve IBVP or BVP for conservation form of 3D Euler  *
//    * equations on multi-block solution-adaptive           *
//    * hexadedral mesh.                                     *
//    ********************************************************/
   continue_existing_calculation: ;
   CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
 
// /* Open residual file and reset the CPU time. */
   first_step = 1;
   limiter_freezing_off = ON;
   
   if (CFFC_Primary_MPI_Processor()) {
      
      error_flag = Open_Progress_File(residual_file,
                                      IPs.Output_File_Name,
                                      number_of_time_steps);
      if (error_flag) {
         cout << "\n ERROR: Unable to open residual file for the"
            "calculation.\n";
         cout.flush();
      } /* endif */
   }

   CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
   CFFC_Broadcast_MPI(&error_flag, 1);
   if (error_flag) return (error_flag);
   processor_cpu_time.reset();

  
//    /* Perform required number of iterations (time steps). */
   if ((!IPs.Time_Accurate &&
        IPs.Maximum_Number_of_Time_Steps > 0 &&
        number_of_time_steps < IPs.Maximum_Number_of_Time_Steps) ||
       (IPs.Time_Accurate &&
        IPs.Time_Max > Time)) {
     
      if (!batch_flag) cout << "\n\n Beginning computations on "
                            << Date_And_Time() << ".\n\n";
      
      if ((!IPs.Time_Accurate &&
           IPs.Maximum_Number_of_Time_Steps > 0 &&
           number_of_time_steps < IPs.Maximum_Number_of_Time_Steps) ||
          (IPs.Time_Accurate &&
           IPs.Time_Max > Time)) {
         perform_explicit_time_marching = ON;
      } else {
         perform_explicit_time_marching = OFF;
      } /* endif */
   
      while (perform_explicit_time_marching) {
        
          /* Determine local and global time steps. */
          dTime =  Local_SolnBlk.CFL(IPs);
          dTime = CFFC_Minimum_MPI(dTime); 
          // Find global minimum time step for all processors.
                
          if (IPs.Time_Accurate) {

             if ((IPs.i_Time_Integration != 
                  TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING) &&
                 (Time + IPs.CFL_Number*dTime > 
                  IPs.Time_Max)) {
                dTime = (IPs.Time_Max-Time)/IPs.CFL_Number;
             } else if (Time + IPs.CFL_Number*dTime*
                        MultiStage_Optimally_Smoothing(IPs.N_Stage, 
                                                       IPs.N_Stage,
                                                       IPs.i_Limiter) > 
                        IPs.Time_Max) {
                dTime = (IPs.Time_Max-Time)/
                   (IPs.CFL_Number*
                    MultiStage_Optimally_Smoothing(IPs.N_Stage, 
                                                   IPs.N_Stage,
                                                   IPs.i_Limiter));
             } /* endif */

          } /* endif */
          if (!IPs.Local_Time_Stepping) { // Set global time step.
                Local_SolnBlk.Set_Global_TimeStep(dTime);
             
              
          } /* endif */

       
          /* Determine the L1, L2, and max norms of the solution residual. */
          residual_l1_norm =  Local_SolnBlk.L1_Norm_Residual();
          residual_l1_norm = CFFC_Summation_MPI(residual_l1_norm); 
          // L1 norm for all processors.
          residual_l2_norm =  Local_SolnBlk.L2_Norm_Residual();
          residual_l2_norm = CFFC_Summation_MPI(residual_l2_norm); 
          
 	 residual_max_norm =  Local_SolnBlk.Max_Norm_Residual();
         
 	 residual_max_norm = CFFC_Maximum_MPI(residual_max_norm);
 	 // Max norm for all processors.
 	 /* Update CPU time used for the calculation so far. */
 	 processor_cpu_time.update();
         total_cpu_time.cput = CFFC_Summation_MPI(processor_cpu_time.cput); 
 	 // Total CPU time for all processors.
         
          /* Periodically save restart solution files. */
          if (!first_step &&
              number_of_time_steps-IPs.Restart_Solution_Save_Frequency*
              (number_of_time_steps/IPs.Restart_Solution_Save_Frequency)
              == 0 ) {
             if (!batch_flag) cout << "\n\n  Saving solution to " 
                                       "restart data file(s) after"
                                   << " n = " << number_of_time_steps 
                                   << " steps (iterations).";
             //    error_flag = Write_Octree(Octree,
             //                                       IPs);
             if (error_flag) {
                cout << "\n ERROR: Unable to open octree data file "
                     << "on processor "<< CFFC_MPI::This_Processor_Number
                     << ".\n";
                cout.flush();
             } /* endif */
            
             if (error_flag) return (error_flag);
             error_flag =  Local_SolnBlk.Write_Restart_Solution
                (IPs,
                 List_of_Local_Solution_Blocks,
                 number_of_time_steps,
                 Time,
                 processor_cpu_time);
          
             if (error_flag) {
                cout << "\n ERROR: Unable to open restart "
                   "output  data file(s) on processor "
                     << CFFC_MPI::This_Processor_Number
                     << ".\n";
                cout.flush();
             } /* endif */
           
            if (error_flag) return (error_flag);
            cout << "\n";
            cout.flush();
 	} /* endif */
       
 	/* Output progress information for the calculation. */
         if (!batch_flag) Output_Progress_L2norm(number_of_time_steps,
                                                 Time*THOUSAND,
                                                 total_cpu_time,
                                                 residual_l2_norm,
                                                 first_step,
                                                 50);
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
         if (!IPs.Time_Accurate &&
             number_of_time_steps >= 
             IPs.Maximum_Number_of_Time_Steps) break;
         if (IPs.Time_Accurate &&
             Time >= IPs.Time_Max) break;
         
         /* Freeze limiters as necessary. */
       if (!first_step &&
           IPs.Freeze_Limiter &&
           limiter_freezing_off &&           
           residual_l2_norm <= IPs.Freeze_Limiter_Residual_Level) {
          //   Freeze_Limiters(Local_SolnBlk, 
          // 			   List_of_Local_Solution_Blocks);
          limiter_freezing_off = ON;	
       } /* endif */
       
 
   
         /* Update solution for next time step using a multistage
          time stepping scheme. */
       for ( i_stage  = 1 ; i_stage <= IPs.N_Stage ; ++i_stage ) {
       
          // 1. Send/Copy the information between blocks...
          
//             Send_All_Boundary_Info(MultiBlock_Connectivity, 
//                                                    IPs);

          Send_All_Messages<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> >(Local_SolnBlk.Hexa_Block_List,
                                                                   List_of_Local_Solution_Blocks,
                                                                   Local_SolnBlk.Hexa_Block_List[0]->NumVar(),
                                                                   OFF);

          //    Local_SolnBlk.Copy_MultiBlk_Boundary_Info( List_of_Local_Solution_Blocks, IPs.IP_Grid);
          
          // 2. Apply boundary conditions for stage.

          Local_SolnBlk.BCs(IPs);

          // 3. Determine solution residuals for stage.
          error_flag =  Local_SolnBlk.dUdt_Multistage_Explicit
             (IPs, i_stage);
          
                 
          if (error_flag) {
             cout<<"\n ERROR: solution residual error on processor"
             <<CFFC_MPI::This_Processor_Number<< ".\n";
             cout.flush();
          } /* endif */
          
          if (error_flag) return (error_flag);
          
          // 7. Update solution for stage.
          error_flag =  Local_SolnBlk.Update_Solution_Multistage_Explicit
             (IPs, i_stage);
          
          if (error_flag) {
             
             cout<<"\n ERROR: solution update error on processor "
             <<CFFC_MPI::This_Processor_Number<<".\n";
             cout.flush();
          } /* endif */
          if (error_flag) return (error_flag);
          
          
       } /* endfor */
        
       /* Update time and time step counter. */
       if (first_step) first_step = 0;
       number_of_time_steps = number_of_time_steps + 1;
       if (IPs.i_Time_Integration != 
           TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING) {
          
          Time = Time + IPs.CFL_Number*dTime;

       } else {
          Time = Time + IPs.CFL_Number*dTime*
             MultiStage_Optimally_Smoothing(IPs.N_Stage, 
                                            IPs.N_Stage,
                                            IPs.i_Limiter);
          
       } /* endif */
       
      } /* endwhile */
      
      if (!batch_flag) cout << "\n\n Computations complete on " 
                            << Date_And_Time() << ".\n";
      
      
      
    	
   } /* endif */
  
//     Send_All_Boundary_Info(MultiBlock_Connectivity, IPs);
   // copying boundary information between blocks

   //   Local_SolnBlk.Copy_MultiBlk_Boundary_Info( List_of_Local_Solution_Blocks, IPs.IP_Grid); 
   

   Send_All_Messages<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> >(Local_SolnBlk.Hexa_Block_List,
                                                            List_of_Local_Solution_Blocks,
                                                            Local_SolnBlk.Hexa_Block_List[0]->NumVar(),
                                                            OFF);
   
   // specifying boundary condition
   Local_SolnBlk.BCs(IPs);

   /* Close residual file. */
   if (CFFC_Primary_MPI_Processor()) error_flag = Close_Progress_File(residual_file);

//    /********************************************************
//    * Solution calculations complete.                      *
//    * Write 3D Euler solution to output and restart files  *
//    * as required, reset solution parameters, and run      *
//    * other cases as specified by input parameters.        *
//    ********************************************************/ 
   
 
  postprocess_current_calculation: ;
   while (1) {
      
      if (CFFC_Primary_MPI_Processor()) {    
         IPs.Get_Next_Input_Control_Parameter();
         command_flag = IPs.Parse_Next_Input_Control_Parameter();
         line_number = IPs.Line_Number;
      }
      
      CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
      IPs.Broadcast_Input_Parameters();
      CFFC_Broadcast_MPI(&command_flag, 1);
      if (command_flag == EXECUTE_CODE) {
         
//          Deallocate memory for 3D Euler equation solution.
         if (!batch_flag) cout<<"\n Deallocating solution variables.";
         // Destructor of the class: Local_Block_List is automatically
         // called when the scope of the calculation is done...
         // ....
         // Output input parameters for new caluculation.
         if (!batch_flag)  {
            cout << "\n\n Starting a new calculation.";
            cout << IPs << "\n";
            cout.flush();
         } /* endif */
         
          // Execute new calculation.
         goto execute_new_calculation;
         
      } else if (command_flag == TERMINATE_CODE) {
         // Deallocate memory for 3D Euler equation solution.
         // Memory for the list of solution blocks is automatically
         // deallocated at the end execution.
         // Close input data file.
         
         
         if (!batch_flag) cout<<"\n Closing input data file.";

         IPs.Close_Input_File();
         cout.flush();

         // if (CFFC_Primary_MPI_Processor()) IPs.Close_Input_File();
         // Terminate calculation
         return (0);
               
         
      } else if (command_flag == WRITE_OUTPUT_CODE) {
         // Output solution data.
         if (!batch_flag) cout<<"\n Writing nodal solution to " 
                             "output data file(s).";
         error_flag =  Local_SolnBlk.Output_Tecplot
            (IPs, List_of_Local_Solution_Blocks,
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
         cout.flush();
         
      } else if (command_flag == WRITE_OUTPUT_CELLS_CODE) {
         // Output solution data.
         if (!batch_flag) cout<<"\n Writing cell-centered  solution to " 
                             "output data file(s).";
         
         error_flag =  Local_SolnBlk.Output_Cells_Tecplot
            (IPs, List_of_Local_Solution_Blocks,
             number_of_time_steps,
             Time);
         
         if (error_flag) {
            cout << "\n  ERROR: Unable to open  cell output "
               "data file(s) on processor "
                 <<CFFC_MPI::This_Processor_Number
                 << ".\n";
            cout.flush();
         } /* endif */   
         
         error_flag = CFFC_OR_MPI(error_flag);
         if (error_flag) return (error_flag);   
         cout.flush();
         
      } else if (command_flag == WRITE_OUTPUT_NODES_CODE) {
         // Output solution data.
         if (!batch_flag) cout<<"\n Writing nodal solution to " 
                                "output data file(s).";
         
         error_flag =  Local_SolnBlk.Output_Nodes_Tecplot
                       (IPs, List_of_Local_Solution_Blocks,
                        number_of_time_steps,
                        Time);
         
         if (error_flag) {
            cout << "\n  ERROR: Unable to open  cell output "
               "data file(s) on processor "
                 << CFFC_MPI::This_Processor_Number
                 << ".\n";
            cout.flush();
         } /* endif */   
         
         error_flag = CFFC_OR_MPI(error_flag);
         if (error_flag) return (error_flag);   
         cout.flush();

      } else if (command_flag == WRITE_RESTART_CODE) {
         // Write restart files.
         if (!batch_flag) cout<<"\n Writing  solution to restart "
                             "data file(s).";
         //  error_flag = Write_Octree(Octree,
         //                                    IPs);
         if (error_flag) {
            cout << "\n  ERROR: Unable to open  octree data file "
                 << "on processor " <<CFFC_MPI::This_Processor_Number
                 << ".\n";
            cout.flush();
         } /* endif */
         
         error_flag =   Local_SolnBlk.Write_Restart_Solution
            (IPs, List_of_Local_Solution_Blocks,
             number_of_time_steps,
             Time,
             processor_cpu_time);
         
         if (error_flag) {
            cout << "\n  ERROR: Unable to open  restart "
               "output data file(s) on processor " 
               //  <<CFFC_MPI::This_Processor_Number
                 << ".\n";
            cout.flush();
            
         } /* endif */
         error_flag = CFFC_OR_MPI(error_flag);
         if (error_flag) return (error_flag);
         
      } else if (command_flag == INVALID_INPUT_CODE ||
                 command_flag == INVALID_INPUT_VALUE) {
         line_number = -line_number;
         cout << "\n  ERROR: Error reading  data at line #"
              << -line_number  << " of input data file.\n";
         cout.flush();
         return (line_number);
      } /* endif */
         
     
   } /* endwhile */
   cout.flush();
 
   /********************************************************  
    * End of all Solver computations and I/O.       *
    ********************************************************/
  
   return (0);
   
}

#endif // _HEXA_SOLVER_INCLUDED

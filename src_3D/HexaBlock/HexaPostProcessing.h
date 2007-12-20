#ifndef _HEXA_POST_PROCESSING_INCLUDED
#define _HEXA_POST_PROCESSING_INCLUDED

/*! *****************************************************
 * Routine: Hexa_Post_Processing                        *
 *                                                      *
 * Post Processes Solver Solution Data base on flags    *
 * in the Input File.  Also allows for the command_flag *
 * to be set to:                                        *
 *                                                      *
 *    EXECUTE_CODE   -> begin new calculation           *
 *    CONTINUE_CODE  -> continue existing calculation   *
 *    TERMINATE_CODE -> Finialize & Exit program        *
 *                                                      *
 ********************************************************/
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
int Hexa_Post_Processing(HexaSolver_Data &Data,
		 	 HexaSolver_Solution_Data<SOLN_pSTATE, SOLN_cSTATE> &Solution_Data) {


  int error_flag(0);
  int line_number;

  error_flag = Hexa_Post_Processing_Specializations(Data,
	  	  	                            Solution_Data);
  if (error_flag) return error_flag;

  while (Solution_Data.command_flag != TERMINATE_CODE) {
    
    if (CFFC_Primary_MPI_Processor()) {    
      Solution_Data.Input.Get_Next_Input_Control_Parameter(true);
      Solution_Data.command_flag = Solution_Data.Input.Parse_Next_Input_Control_Parameter();
      line_number = Solution_Data.Input.Line_Number;
    } /* endif */
  
    CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
    Solution_Data.Input.Broadcast();
    CFFC_Broadcast_MPI(&Solution_Data.command_flag, 1);

    if (Solution_Data.command_flag == EXECUTE_CODE) {
      // Deallocate memory for exiting solution blocks.
      if (!Data.batch_flag) {
	cout <<"\n Deallocating existing solution blocks.";
	cout.flush();
      } /* endif */
      
         // Destructor of the class: Local_Block_List is automatically
         // called when the scope of the calculation is done...
         // ....

         // Output input parameters for new caluculation.
         if (!Data.batch_flag)  {
            cout << "\n\n Starting a new calculation.";
            cout << Solution_Data.Input << "\n";
            cout.flush();
         } /* endif */
         
         // Execute new calculation.
	 return error_flag;
         
      } else if (Solution_Data.command_flag == TERMINATE_CODE) {
         // Deallocate memory for existing solution blocks.
         // Memory for the list of solution blocks is automatically
         // deallocated at the end execution.

         // Close input data file.         
	 if (!Data.batch_flag) {
            cout <<"\n Closing input data file.";
            cout.flush();
         } /* endif */
         if (CFFC_Primary_MPI_Processor()) Solution_Data.Input.Close_Input_File();

         // Terminate calculation
         return error_flag;
          
    } else if (Solution_Data.command_flag == CONTINUE_CODE) {

      // Adjust maximum time step counters.                     
      Solution_Data.Input.Maximum_Number_of_Time_Steps += Data.number_of_explicit_time_steps;
      Solution_Data.Input.NKS_IP.Maximum_Number_of_NKS_Iterations += Data.number_of_implicit_time_steps;
   
      // Output input parameters for continuing calculation.
        if (!Data.batch_flag) {
	   cout << "\n\n Continuing existing calculation.";
	   cout << Solution_Data.Input << "\n";
	   cout.flush();
        } /* endif */
      
        // Continue existing calculation. 
	return error_flag;

      } else if (Solution_Data.command_flag == WRITE_OUTPUT_CODE) {
         // Output solution data.
	 if (!Data.batch_flag) {
           cout << "\n Writing solution to output data file(s).";
           cout.flush();
           cout.flush();
 	 } /* endif */

         error_flag = Solution_Data.Local_Solution_Blocks.Output_Tecplot(Solution_Data.Input, 
									 Data.Local_Adaptive_Block_List,
									 Data.total_number_of_time_steps(),
									 Data.Time);
       
       
         
         if (error_flag) {
            cout << "\n  ERROR: Unable to open  node output " 
                         "data file(s) on processor "
             <<CFFC_MPI::This_Processor_Number<< ".\n";
            cout.flush();
         } /* endif */
         error_flag = CFFC_OR_MPI(error_flag);
         if (error_flag) return (error_flag);   
         
      } else if (Solution_Data.command_flag == WRITE_OUTPUT_CELLS_CODE) {
         // Output solution data.
	 if (!Data.batch_flag) {
           cout <<"\n Writing cell-centered solution to output data file(s).";
           cout.flush();
         } /* endif */
         
         error_flag = Solution_Data.Local_Solution_Blocks.Output_Cells_Tecplot(Solution_Data.Input, 
									       Data.Local_Adaptive_Block_List,
									       Data.total_number_of_time_steps(),
									       Data.Time);
       
         if (error_flag) {
            cout << "\n  ERROR: Unable to open  cell output data file(s) on processor "
                 << CFFC_MPI::This_Processor_Number
                 << ".\n";
            cout.flush();
         } /* endif */   
         error_flag = CFFC_OR_MPI(error_flag);
         if (error_flag) return (error_flag);   
         
      } else if (Solution_Data.command_flag == WRITE_OUTPUT_NODES_CODE) {
         // Output solution data.
	 if (!Data.batch_flag) {
            cout<<"\n Writing nodal solution to output data file(s).";
            cout.flush();
         } /* endif */
         
         error_flag = Solution_Data.Local_Solution_Blocks.Output_Nodes_Tecplot(Solution_Data.Input, 
									       Data.Local_Adaptive_Block_List,
									       Data.total_number_of_time_steps(),
									       Data.Time);
         if (error_flag) {
            cout << "\n  ERROR: Unable to open  cell output data file(s) on processor "
                 << CFFC_MPI::This_Processor_Number
                 << ".\n";
            cout.flush();
         } /* endif */   
         error_flag = CFFC_OR_MPI(error_flag);
         if (error_flag) return (error_flag);   

      } else if (Solution_Data.command_flag == WRITE_RESTART_CODE) {
         // Write restart files.
	 if (!Data.batch_flag) {
            cout <<"\n Writing solution to restart data file(s).";
	    cout.flush();
         } /* endif */

         error_flag = Write_Octree(Data.Octree, Solution_Data.Input);
         if (error_flag) {
            cout << "\n  ERROR: Unable to open  octree data file on processor " 
                 << CFFC_MPI::This_Processor_Number
                 << ".\n";
            cout.flush();
         } /* endif */
         error_flag = CFFC_OR_MPI(error_flag);
         if (error_flag) return (error_flag);
         
         error_flag = Solution_Data.Local_Solution_Blocks.Write_Restart_Solution(Solution_Data.Input, 
										 Data.Local_Adaptive_Block_List,
										 Data.total_number_of_time_steps(),
										 Data.Time,
										 Data.processor_cpu_time);
         if (error_flag) {
            cout << "\n  ERROR: Unable to open restart output data file(s) on processor " 
                 << CFFC_MPI::This_Processor_Number
                 << ".\n";
            cout.flush();
         } /* endif */
         error_flag = CFFC_OR_MPI(error_flag);
         if (error_flag) return (error_flag);

      }else if (Solution_Data.command_flag == READ_RESTART_CODE){
        //Read restart files.
        if (!Data.batch_flag){ 
            cout << "\n Reading solution from restart data files."; 
            cout.flush();
        }
        //error_flag = Read_Octree(Octree,
   	//		           Solution_Data.Input);
        if (error_flag) {
           cout << "\n ERROR: Unable to open Octree data file "
   	        << "on processor "
  	        << Data.Local_Adaptive_Block_List.ThisCPU
 	        << ".\n";
           cout.flush();
        } // endif 
        error_flag = CFFC_OR_MPI(error_flag);
        if (error_flag) return (error_flag);

        error_flag = Solution_Data.Local_Solution_Blocks.Read_Restart_Solution(Solution_Data.Input,  
									       Data.Local_Adaptive_Block_List,
									       Data.number_of_explicit_time_steps,  //??? FIX        
									       Data.Time,
									       Data.processor_cpu_time);
        if (error_flag) {
           cout << "\n  ERROR: Unable to open restart input data file(s) on processor "
                << CFFC_MPI::This_Processor_Number << ".\n";
           cout.flush();
        } /* endif */
        error_flag = CFFC_OR_MPI(error_flag);
        if (error_flag) return (error_flag);

      } else if (Solution_Data.command_flag == MORTON_ORDERING_CODE) {
         if (!Data.batch_flag){ 
             cout << "\n\n Applying Morton re-ordering algorithm.";
             cout.flush();
         }
         error_flag = Morton_ReOrdering_of_Solution_Blocks(Data.Octree,
				       			   Data.Local_Adaptive_Block_List, 
                                                           Solution_Data.Local_Solution_Blocks, 
                                                           Solution_Data.Input, 
                                                           Data.number_of_explicit_time_steps,  //??? FIX
                                                           Data.Time, 
                                                           Data.processor_cpu_time);
         if (error_flag) {
            cout <<"\n ERROR: Morton re-ordering error on processor "
                 << Data.Local_Adaptive_Block_List.ThisCPU
                 << ".\n";
            cout.flush();
            return (error_flag);
         } // endif 
         error_flag = CFFC_OR_MPI(error_flag);
         if (error_flag) return (error_flag);

         //Output space filling curve in Tecplot format
         if (!Data.batch_flag) cout << "\n Outputting space filling curve showing block loading for CPUs.";
         Morton_SFC_Output_Tecplot3D(Solution_Data.Input,
                                     Solution_Data.Local_Solution_Blocks,
				     Data.Local_Adaptive_Block_List);
    
      } else if (Solution_Data.command_flag == INVALID_INPUT_CODE ||
                 Solution_Data.command_flag == INVALID_INPUT_VALUE) {
         line_number = -line_number;
         cout << "\n  ERROR: Error reading  data at line #"
              << -line_number  << " of input data file.\n";
         cout.flush();
         return (line_number);  
      } /* endif */
         
   } /* endwhile */
  
   return error_flag;

}


/*! *****************************************************
 * Routine: Hexa_Post_Processing_Specializations        *
 *                                                      *
 *                                                      *
 ********************************************************/
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
int Hexa_Post_Processing_Specializations(HexaSolver_Data &Data,
		 	                 HexaSolver_Solution_Data<SOLN_pSTATE, SOLN_cSTATE> &Solution_Data) {


   int error_flag(0);
  
   return error_flag;

}

#endif // _HEXA_POST_PROCESSING_INCLUDED

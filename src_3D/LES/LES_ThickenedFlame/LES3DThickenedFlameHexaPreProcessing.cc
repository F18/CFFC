
#ifndef _LES3DTF_INCLUDED
#include "LES3DThickenedFlame.h"
#endif // _LES3DTF_INCLUDED

/********************************************************
 * Routine: Pre_Processing_Specializations              *
 ********************************************************/
template<>
int Hexa_Pre_Processing_Specializations(HexaSolver_Data &Data,
                                        HexaSolver_Solution_Data<LES3DTF_pState, 
                                                                 LES3DTF_cState> &Solution_Data) {

  int error_flag(0);
  
  RandomFieldRogallo<LES3DTF_pState, LES3DTF_cState>   Velocity_Field_Type(Solution_Data.Input);

  if (Solution_Data.Input.Grid_IP.i_Grid != GRID_BUNSEN_BURNER  &&
      Solution_Data.Input.Grid_IP.i_Grid != GRID_BUNSEN_BOX) {

    Data.Velocity_Field.Create(Data.Initial_Mesh, 
			       Solution_Data.Input.Grid_IP);
    
    error_flag = Velocity_Field_Type.Create_Homogeneous_Turbulence_Velocity_Field(Data.Initial_Mesh, 
										  Solution_Data.Input.Grid_IP,
										  Data.batch_flag,
										  Data.Velocity_Field);
  } else {

    Data.Velocity_Field.Create(Data.Auxiliary_Mesh, 
			       Solution_Data.Input.Grid_IP);

    error_flag = Velocity_Field_Type.Create_Homogeneous_Turbulence_Velocity_Field(Data.Auxiliary_Mesh, 
										  Solution_Data.Input.Grid_IP,
										  Data.batch_flag,
										  Data.Velocity_Field);
  }

  if (error_flag) return error_flag;
  

  // If required, do the interpolation of the turbulent field

//   if (Solution_Data.Input.Grid_IP.i_Grid == GRID_BUNSEN_BURNER ||
//       Solution_Data.Input.Grid_IP.i_Grid == GRID_BUNSEN_BOX) {

//     Turbulent_Velocity_Field_Multi_Block_List  Interpolated_Velocity_Field;

//     Interpolated_Velocity_Field.Create(Data.Initial_Mesh,   
// 				       Solution_Data.Input.Grid_IP);
//     Data.Velocity_Field.Interpolate_Turbulent_Field(Data.Initial_Mesh, 
// 					       Interpolated_Velocity_Field);

//     Assign_Homogeneous_Turbulence_Velocity_Field(Solution_Data.Local_Solution_Blocks.Soln_Blks,
// 						 Data.Local_Adaptive_Block_List,
// 						 Interpolated_Velocity_Field);
//   } else {
//     Assign_Homogeneous_Turbulence_Velocity_Field(Solution_Data.Local_Solution_Blocks.Soln_Blks,
// 						 Data.Local_Adaptive_Block_List,
// 						 Data.Velocity_Field);
//   }


  if (Solution_Data.Input.i_ICs != IC_RESTART) {

    // assign turbulent velocity fluctuations to the mixture of fresh gases
    if (Solution_Data.Input.Grid_IP.i_Grid == GRID_BUNSEN_BURNER  ||
	Solution_Data.Input.Grid_IP.i_Grid == GRID_BUNSEN_BOX) {

      error_flag = IC_Assign_Turbulence_Fresh_Gas(Solution_Data.Local_Solution_Blocks.Soln_Blks,
						  Data.Local_Adaptive_Block_List,
						  Data.Velocity_Field,
						  Solution_Data.Input);
      if (error_flag) return error_flag;
    }

    // ICs specializations
  error_flag = Solution_Data.Local_Solution_Blocks.ICs_Specializations(Solution_Data.Input);

  if (error_flag) return error_flag;
  }

  return error_flag;

}

/********************************************************
 * Routine: Hexa_Post_Processing_Specializations        *
 ********************************************************/
template<>
int Hexa_Post_Processing_Specializations(HexaSolver_Data &Data,
		 	                 HexaSolver_Solution_Data<LES3DTF_pState, 
                                                                  LES3DTF_cState>&Solution_Data) {

   int error_flag(0);
   double u_ave, v_ave, w_ave, sqr_u;
   
   if (Solution_Data.Input.Species_IP.reaction_mechanism_name != "NO_REACTIONS") {

     Conditional_Averaging_of_Velocity_Field(Solution_Data.Local_Solution_Blocks.Soln_Blks,
					     Data.Local_Adaptive_Block_List,
					     u_ave,
					     v_ave,
					     w_ave);

     Conditional_Averaging_of_Solution(Solution_Data.Local_Solution_Blocks.Soln_Blks,
				       Data.Local_Adaptive_Block_List,
				       u_ave,
				       v_ave,
				       w_ave,
				       sqr_u);
   } else {

   Time_Averaging_of_Velocity_Field(Solution_Data.Local_Solution_Blocks.Soln_Blks,
                                    Data.Local_Adaptive_Block_List,
                                    u_ave,
                                    v_ave,
                                    w_ave);

   Time_Averaging_of_Solution(Solution_Data.Local_Solution_Blocks.Soln_Blks,
                              Data.Local_Adaptive_Block_List,
                              u_ave,
                              v_ave,
                              w_ave,
                              sqr_u);
   }

   return error_flag;

}


/*****************************************************************
 * Routine: Open_Other_Solution_Progress_Specialization_Files    *
 *****************************************************************/
template<>
int Open_Other_Solution_Progress_Specialization_Files(HexaSolver_Data &Data,
                                                      HexaSolver_Solution_Data<LES3DTF_pState,
						                               LES3DTF_cState> &Solution_Data) {
   int error_flag(0);

   error_flag = Open_Turbulence_Progress_File(Data.other_solution_progress_files[0],
   	 		                      Solution_Data.Input.Output_File_Name,
   				              Data.number_of_explicit_time_steps);

   return error_flag;

}

/*****************************************************************
 * Routine: Close_Other_Solution_Progress_Specialization_Files   *
 *****************************************************************/
template<>
int Close_Other_Solution_Progress_Specialization_Files(HexaSolver_Data &Data,
                                                       HexaSolver_Solution_Data<LES3DTF_pState,
						                                LES3DTF_cState> &Solution_Data) {

   int error_flag(0);

   error_flag = Close_Turbulence_Progress_File(Data.other_solution_progress_files[0]);

   return error_flag;

}

/*****************************************************************
 * Routine: Output_Other_Solution_Progress_Specialization_Data   *
 *****************************************************************/
template<>
int Output_Other_Solution_Progress_Specialization_Data(HexaSolver_Data &Data,
                                                       HexaSolver_Solution_Data<LES3DTF_pState,
						                                LES3DTF_cState> &Solution_Data) {

   int error_flag(0);
   double total_TKE, total_enstrophy, u_prime, Taylor_scale, viscosity, turb_burning_rate;

   // Calculate various turbulence quantities
   if (Solution_Data.Input.Species_IP.reaction_mechanism_name != "NO_REACTIONS") {

     total_TKE = Conditional_Total_TKE<Hexa_Block<LES3DTF_pState, LES3DTF_cState> >(Solution_Data.Local_Solution_Blocks.Soln_Blks,
										    Data.Local_Adaptive_Block_List);

     total_enstrophy = Conditional_Total_Enstrophy<Hexa_Block<LES3DTF_pState, LES3DTF_cState> >(Solution_Data.Local_Solution_Blocks.Soln_Blks,
												Data.Local_Adaptive_Block_List);

     u_prime = Conditional_u_rms<Hexa_Block<LES3DTF_pState, LES3DTF_cState> >(Solution_Data.Local_Solution_Blocks.Soln_Blks,
									      Data.Local_Adaptive_Block_List);

     Taylor_scale = Conditional_Taylor_Scale<Hexa_Block<LES3DTF_pState, LES3DTF_cState> >(Solution_Data.Local_Solution_Blocks.Soln_Blks,
											  Data.Local_Adaptive_Block_List);

     viscosity = Conditional_Average_viscosity<Hexa_Block<LES3DTF_pState, LES3DTF_cState> >(Solution_Data.Local_Solution_Blocks.Soln_Blks,
											    Data.Local_Adaptive_Block_List);
   } else {

   total_TKE = Total_TKE<Hexa_Block<LES3DTF_pState, LES3DTF_cState> >(Solution_Data.Local_Solution_Blocks.Soln_Blks,
   					                                Data.Local_Adaptive_Block_List);

   total_enstrophy = Total_Enstrophy<Hexa_Block<LES3DTF_pState, LES3DTF_cState> >(Solution_Data.Local_Solution_Blocks.Soln_Blks,
   			 			                                    Data.Local_Adaptive_Block_List);

   u_prime = u_rms<Hexa_Block<LES3DTF_pState, LES3DTF_cState> >(Solution_Data.Local_Solution_Blocks.Soln_Blks,
   					                          Data.Local_Adaptive_Block_List);

   Taylor_scale = Taylor_Scale<Hexa_Block<LES3DTF_pState, LES3DTF_cState> >(Solution_Data.Local_Solution_Blocks.Soln_Blks,
   					                                      Data.Local_Adaptive_Block_List);

   viscosity = Average_viscosity<Hexa_Block<LES3DTF_pState, LES3DTF_cState> >(Solution_Data.Local_Solution_Blocks.Soln_Blks,
   					                                        Data.Local_Adaptive_Block_List);
   }

//    turb_burning_rate = Turbulent_Burning_Rate<Hexa_Block<LES3DTF_pState, LES3DTF_cState> >(Solution_Data.Local_Solution_Blocks.Soln_Blks,
// 											   Data.Local_Adaptive_Block_List,
// 											   Solution_Data.Input.Grid_IP);

   // Output turbulence statistics data to turbulence progress variable file
//    if (CFFC_Primary_MPI_Processor()) {
//      error_flag = Output_Turbulence_Progress_to_File(Data.other_solution_progress_files[0],
//                                                      Data.number_of_explicit_time_steps,
// 			   		             Data.Time*THOUSAND,
//  					             Data.total_cpu_time,
// 					             total_TKE,
//                                                      u_prime,
// 					             total_enstrophy,
// 					             Taylor_scale,
// 					             viscosity,
// 					             turb_burning_rate);
//    } /* endif */
   error_flag = CFFC_OR_MPI(error_flag);

   // Return error flag
   return error_flag;

}

/********************************************************
 * Routine: BCs_Specializations                         *
 ********************************************************/
template<>
int Hexa_BCs_Specializations(HexaSolver_Data &Data,
			     HexaSolver_Solution_Data<LES3DTF_pState, 
			                              LES3DTF_cState> &Solution_Data) {
  int error_flag(0);
  
  error_flag = Inflow_Turbulence_XY_Plane(Solution_Data.Local_Solution_Blocks.Soln_Blks,
					  Data.Local_Adaptive_Block_List,
					  Data.Velocity_Field,
					  Solution_Data.Input,
					  Data.Time);

  return error_flag;

}

/********************************************************
 * Routine: Initialize_Solution_Blocks_Specializations  *
 ********************************************************/
template<>
int Initialize_Solution_Blocks_Specializations(HexaSolver_Data &Data,
					       HexaSolver_Solution_Data<LES3DTF_pState, 
			                                                LES3DTF_cState> &Solution_Data) {
  int error_flag(0);
  
  /* Create the auxiliary mesh on the primary MPI processor. */

  if (CFFC_Primary_MPI_Processor()) {
    
    if (Solution_Data.Input.Grid_IP.i_Grid == GRID_BUNSEN_BURNER) {

      Solution_Data.Input.Grid_IP.i_Grid = GRID_TURBULENCE_BOX;
      Data.Auxiliary_Mesh.Create_Grid(Solution_Data.Input.Grid_IP);
      Solution_Data.Input.Grid_IP.i_Grid = GRID_BUNSEN_BURNER;

    } else if (Solution_Data.Input.Grid_IP.i_Grid == GRID_BUNSEN_BOX) {

      Solution_Data.Input.Grid_IP.i_Grid = GRID_TURBULENCE_BOX;
      Data.Auxiliary_Mesh.Create_Grid(Solution_Data.Input.Grid_IP);
      Solution_Data.Input.Grid_IP.i_Grid = GRID_BUNSEN_BOX;
    }   

  } /* endif */

  CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
  
  /* Broadcast the auxiliary mesh to other MPI processors. */

  if (Solution_Data.Input.Grid_IP.i_Grid == GRID_BUNSEN_BURNER) {

    Solution_Data.Input.Grid_IP.i_Grid = GRID_TURBULENCE_BOX;
    Data.Auxiliary_Mesh.Broadcast();
    Solution_Data.Input.Grid_IP.i_Grid = GRID_BUNSEN_BURNER;

  } else if(Solution_Data.Input.Grid_IP.i_Grid == GRID_BUNSEN_BOX) {

    Solution_Data.Input.Grid_IP.i_Grid = GRID_TURBULENCE_BOX;
    Data.Auxiliary_Mesh.Broadcast();
    Solution_Data.Input.Grid_IP.i_Grid = GRID_BUNSEN_BOX;
  }                  
 
  return error_flag;

}

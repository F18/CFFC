
#ifndef _LES3DFsd_INCLUDED
#include "LES3DFsd.h"
#endif // _LES3DFsd_INCLUDED


/********************************************************
 * Routine: Pre_Processing_Specializations              *
 ********************************************************/
template<>
int Hexa_Pre_Processing_Specializations(HexaSolver_Data &Data,
                                        HexaSolver_Solution_Data<LES3DFsd_pState, 
                                                                 LES3DFsd_cState> &Solution_Data) {

  int error_flag(0);
  //-------------------------------------------------
  //      RESTART
  //-------------------------------------------------
  if (Solution_Data.Input.i_ICs == IC_RESTART) {
  
    if (Solution_Data.Input.Grid_IP.i_Grid == GRID_BUNSEN_BURNER  ||
	Solution_Data.Input.Grid_IP.i_Grid == GRID_BUNSEN_BOX) {

      // Read the turbulent velocity field data (primary processor)
      if ( CFFC_Primary_MPI_Processor() ) {
	Data.Velocity_Field.Read_Turbulent_Velocity_Field();
      }

    }

  //-------------------------------------------------
  //      NON-RESTART
  //-------------------------------------------------
  } else {

    // the primary processor creates the initial turbulence
    if ( CFFC_Primary_MPI_Processor() ) {
      
      RandomFieldRogallo<LES3DFsd_pState, LES3DFsd_cState>   Velocity_Field_Type(Solution_Data.Input); 

      // use initial mesh
      if (Solution_Data.Input.Grid_IP.i_Grid != GRID_BUNSEN_BURNER  &&
	  Solution_Data.Input.Grid_IP.i_Grid != GRID_BUNSEN_BOX) {

	Data.Velocity_Field.Create(Data.Initial_Mesh, 
				   Solution_Data.Input.Grid_IP);

	error_flag = Velocity_Field_Type.Create_Homogeneous_Turbulence_Velocity_Field(Data.Initial_Mesh, 
										      Solution_Data.Input.Grid_IP,
										      Data.batch_flag,
      										      Data.Velocity_Field);
      // use auxiliary mesh
      } else {
	
	Data.Velocity_Field.Create(Data.Auxiliary_Mesh, 
				   Solution_Data.Input.Grid_IP);

	error_flag = Velocity_Field_Type.Create_Homogeneous_Turbulence_Velocity_Field(Data.Auxiliary_Mesh, 
										      Solution_Data.Input.Grid_IP,
										      Data.batch_flag,
										      Data.Velocity_Field);
      }

      if (error_flag) return error_flag;

      // Write the turbulent velocity field data
      if (Solution_Data.Input.Grid_IP.i_Grid == GRID_BUNSEN_BURNER  ||
	  Solution_Data.Input.Grid_IP.i_Grid == GRID_BUNSEN_BOX) {
      
	Data.Velocity_Field.Write_Turbulent_Velocity_Field();
      }
    }    

  } /* endif */
  

  CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
  
  /* Broadcast the turbulent velocity field to other MPI processors. */
  
  Data.Velocity_Field.Broadcast();
 

  //-------------------------------------------------
  //  Assign turbulent velocity fluctuations
  //-------------------------------------------------
  if (Solution_Data.Input.i_ICs == IC_TURBULENT_BOX ||
      Solution_Data.Input.i_ICs == IC_TURBULENT_PREMIXED_FLAME) {

    // assign turbulent velocity fluctuations
    Assign_Homogeneous_Turbulence_Velocity_Field(Solution_Data.Local_Solution_Blocks.Soln_Blks,
						 Data.Local_Adaptive_Block_List,
						 Data.Velocity_Field);
  
  } else if (Solution_Data.Input.i_ICs == IC_TURBULENT_BUNSEN_BOX ||
	     Solution_Data.Input.i_ICs == IC_TURBULENT_BUNSEN_FLAME) {

    // assign turbulent velocity fluctuations to the mixture of fresh gases
    error_flag = IC_Assign_Turbulence_Fresh_Gas(Solution_Data.Local_Solution_Blocks.Soln_Blks,
						Data.Local_Adaptive_Block_List,
						Data.Velocity_Field,
						Solution_Data.Input);
      
    if (error_flag) return error_flag;
   
  } else {

  RandomFieldRogallo<LES3DFsd_pState, LES3DFsd_cState>   Velocity_Field_Type(Solution_Data.Input);
  Turbulent_Velocity_Field_Multi_Block_List  Velocity_Field;

  error_flag = Velocity_Field_Type.Create_Homogeneous_Turbulence_Velocity_Field(Data.Initial_Mesh, 
										Solution_Data.Input.Grid_IP,
										Data.batch_flag,
										Data.Velocity_Field);

  if (error_flag) return error_flag;

// If required, do the interpolation of the turbulent field
//     Turbulent_Velocity_Field_Multi_Block_List  Interpolated_Velocity_Field;

//     Interpolated_Velocity_Field.Create(Data.Initial_Mesh,   
// 				       Solution_Data.Input.Grid_IP);

//     Data.Velocity_Field.Interpolate_Turbulent_Field(Data.Initial_Mesh, 
// 						    Interpolated_Velocity_Field);

//     Assign_Homogeneous_Turbulence_Velocity_Field(Solution_Data.Local_Solution_Blocks.Soln_Blks,
// 						 Data.Local_Adaptive_Block_List,
// 						 Interpolated_Velocity_Field);

  }

  //-------------------------------------------------
  // ICs specializations
  //-------------------------------------------------
  error_flag = Solution_Data.Local_Solution_Blocks.ICs_Specializations(Solution_Data.Input);

  if (error_flag) return error_flag;    
    
  return error_flag;

}

/********************************************************
 * Routine: Hexa_Post_Processing_Specializations        *
 ********************************************************/
template<>
int Hexa_Post_Processing_Specializations(HexaSolver_Data &Data,
		 	                 HexaSolver_Solution_Data<LES3DFsd_pState, 
                                                                  LES3DFsd_cState>&Solution_Data) {

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

   if (Data.Time == 0.0  &&  Data.number_of_explicit_time_steps == 0) { 
     Max_and_Min_Cell_Volumes(Solution_Data.Local_Solution_Blocks.Soln_Blks,
			      Data.Local_Adaptive_Block_List);
   }
   
   return error_flag;

}


/*****************************************************************
 * Routine: Open_Other_Solution_Progress_Specialization_Files    *
 *****************************************************************/
template<>
int Open_Other_Solution_Progress_Specialization_Files(HexaSolver_Data &Data,
                                                      HexaSolver_Solution_Data<LES3DFsd_pState,
						                               LES3DFsd_cState> &Solution_Data) {
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
                                                       HexaSolver_Solution_Data<LES3DFsd_pState,
						                                LES3DFsd_cState> &Solution_Data) {

   int error_flag(0);

   error_flag = Close_Turbulence_Progress_File(Data.other_solution_progress_files[0]);

   return error_flag;

}

/*****************************************************************
 * Routine: Output_Other_Solution_Progress_Specialization_Data   *
 *****************************************************************/
template<>
int Output_Other_Solution_Progress_Specialization_Data(HexaSolver_Data &Data,
                                                       HexaSolver_Solution_Data<LES3DFsd_pState,
						                                LES3DFsd_cState> &Solution_Data) {

   int error_flag(0);
   double total_TKE, total_enstrophy, u_prime, Taylor_scale, viscosity, turb_burning_rate;

   // Calculate various turbulence quantities

   if (Solution_Data.Input.Species_IP.reaction_mechanism_name != "NO_REACTIONS") {

     total_TKE = Conditional_Total_TKE<Hexa_Block<LES3DFsd_pState, LES3DFsd_cState> >(Solution_Data.Local_Solution_Blocks.Soln_Blks,
										    Data.Local_Adaptive_Block_List);

     total_enstrophy = Conditional_Total_Enstrophy<Hexa_Block<LES3DFsd_pState, LES3DFsd_cState> >(Solution_Data.Local_Solution_Blocks.Soln_Blks,
												Data.Local_Adaptive_Block_List);

     u_prime = Conditional_u_rms<Hexa_Block<LES3DFsd_pState, LES3DFsd_cState> >(Solution_Data.Local_Solution_Blocks.Soln_Blks,
									      Data.Local_Adaptive_Block_List);

     Taylor_scale = Conditional_Taylor_Scale<Hexa_Block<LES3DFsd_pState, LES3DFsd_cState> >(Solution_Data.Local_Solution_Blocks.Soln_Blks,
											  Data.Local_Adaptive_Block_List);

     viscosity = Conditional_Average_viscosity<Hexa_Block<LES3DFsd_pState, LES3DFsd_cState> >(Solution_Data.Local_Solution_Blocks.Soln_Blks,
											    Data.Local_Adaptive_Block_List);
   } else {

     total_TKE = Total_TKE<Hexa_Block<LES3DFsd_pState, LES3DFsd_cState> >(Solution_Data.Local_Solution_Blocks.Soln_Blks,
  					                                Data.Local_Adaptive_Block_List);

     total_enstrophy = Total_Enstrophy<Hexa_Block<LES3DFsd_pState, LES3DFsd_cState> >(Solution_Data.Local_Solution_Blocks.Soln_Blks,
   			 			                                    Data.Local_Adaptive_Block_List);

     u_prime = u_rms<Hexa_Block<LES3DFsd_pState, LES3DFsd_cState> >(Solution_Data.Local_Solution_Blocks.Soln_Blks,
   					                          Data.Local_Adaptive_Block_List);

     Taylor_scale = Taylor_Scale<Hexa_Block<LES3DFsd_pState, LES3DFsd_cState> >(Solution_Data.Local_Solution_Blocks.Soln_Blks,
   					                                      Data.Local_Adaptive_Block_List);

     viscosity = Average_viscosity<Hexa_Block<LES3DFsd_pState, LES3DFsd_cState> >(Solution_Data.Local_Solution_Blocks.Soln_Blks,
   					                                        Data.Local_Adaptive_Block_List);
   }

   turb_burning_rate = Turbulent_Burning_Rate(Solution_Data.Local_Solution_Blocks.Soln_Blks,
					      Data.Local_Adaptive_Block_List,
					      Solution_Data.Input.Grid_IP);

   // Output turbulence statistics data to turbulence progress variable file
   if (CFFC_Primary_MPI_Processor()) {
     error_flag = Output_Turbulence_Progress_to_File(Data.other_solution_progress_files[0],
                                                     Data.number_of_explicit_time_steps,
			   		             Data.Time*THOUSAND,
 					             Data.total_cpu_time,
					             total_TKE,
                                                     u_prime,
					             total_enstrophy,
					             Taylor_scale,
					             viscosity,
					             turb_burning_rate);
   } /* endif */

   error_flag = CFFC_OR_MPI(error_flag);

   // Return error flag
   return error_flag;

}

/********************************************************
 * Routine: BCs_Specializations                         *
 ********************************************************/
template<>
int Hexa_BCs_Specializations(HexaSolver_Data &Data,
			     HexaSolver_Solution_Data<LES3DFsd_pState, 
			                              LES3DFsd_cState> &Solution_Data) {
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
					       HexaSolver_Solution_Data<LES3DFsd_pState, 
			                                                LES3DFsd_cState> &Solution_Data) {
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

/********************************************************
 *       Burning rate                                   *
 ********************************************************/
template<>
double Turbulent_Burning_Rate(Hexa_Block<LES3DFsd_pState, LES3DFsd_cState> *Solution_Block,
			      AdaptiveBlock3D_List &LocalSolnBlockList,
			      Grid3D_Input_Parameters &IPs) {

  double local_vol, Yf_u, rho_u, burning_rate(ZERO), iso_surface_area(ZERO);
  double flame_height = 0.005;
  Yf_u = 0.05518;//Fresh_Fuel_Mass_Fraction;
  rho_u = 1.13;//Fresh_Density;

  for (int p = 0 ; p <= LocalSolnBlockList.Nblk-1 ; p++ ) {
    if (LocalSolnBlockList.Block[p].used == ADAPTIVEBLOCK3D_USED) {
      for (int i = Solution_Block[p].ICl ; i <= Solution_Block[p].ICu ; i++) {
        for (int j = Solution_Block[p].JCl ; j <= Solution_Block[p].JCu ; j++) {
           for (int k = Solution_Block[p].KCl ; k <= Solution_Block[p].KCu ; k++) {
	     local_vol = Solution_Block[p].Grid.volume(i,j,k);
 	     burning_rate +=  Solution_Block[p].W[i][j][k].Fsd*Solution_Block[p].W[i][j][k].rho*local_vol; 
	     if (Solution_Block[p].W[i][j][k].C <= 0.5 && 
                 Solution_Block[p].Grid.Cell[i][j][k].Xc.z > flame_height) {
	         flame_height = Solution_Block[p].Grid.Cell[i][j][k].Xc.z;
/* 	       iso_surface_area += propagation_dir_area(Solution_Block[p], i, j, k); */
	     }
	   }
	}
      }
    }
  }
  burning_rate = CFFC_Summation_MPI(burning_rate);
  burning_rate = burning_rate*0.403/(PI*0.0056*(0.0056+2.0*flame_height));
//   iso_surface_area = CFFC_Summation_MPI(iso_surface_area);
  
//   double ref_area, Lx, Ly, Lz;

//   if ( IPs.i_Grid == GRID_PERIODIC_BOX_WITH_INFLOW  ||
//        IPs.i_Grid == GRID_PERIODIC_BOX) {
//     Ly = IPs.Box_Height;
//     Lz = IPs.Box_Length;
//     ref_area = Ly*Lz;
//   } else if ( IPs.i_Grid == GRID_BUNSEN_BOX ) {
//     Lx = IPs.Box_Width;
//     ref_area = (0.025 + 2.0*0.02)*Lx;
//   }

//   if ( iso_surface_area > ref_area ) {
//     burning_rate = -burning_rate/(rho_u*Yf_u*iso_surface_area);
//   } else {
//     burning_rate = -burning_rate/(rho_u*Yf_u*ref_area);
//   }


  return burning_rate;
}

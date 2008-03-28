
#ifndef _LES3DFSD_INCLUDED
#include "LES3DFsd.h"
#endif // _LES3DFSD_INCLUDED

/********************************************************
 * Routine: Pre_Processing_Specializations              *
 ********************************************************/
template<>
int Hexa_Pre_Processing_Specializations(HexaSolver_Data &Data,
                                        HexaSolver_Solution_Data<LES3DFsd_pState, 
                                                                 LES3DFsd_cState> &Solution_Data) {

  int error_flag(0);
  
  RandomFieldRogallo<LES3DFsd_pState, LES3DFsd_cState>   Velocity_Field_Type(Solution_Data.Input);
  Turbulent_Velocity_Field_Multi_Block_List  Velocity_Field;

  Velocity_Field.Create(Data.Initial_Mesh, 
                        Solution_Data.Input.Grid_IP);

  error_flag = Velocity_Field_Type.Create_Homogeneous_Turbulence_Velocity_Field(Data.Initial_Mesh, 
										Solution_Data.Input.Grid_IP,
										Data.batch_flag,
										Data.Velocity_Field);
  if (error_flag) return error_flag;

  Assign_Homogeneous_Turbulence_Velocity_Field(Solution_Data.Local_Solution_Blocks.Soln_Blks,
                                               Data.Local_Adaptive_Block_List,
                                               Velocity_Field);

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
 *       Burning rate                                   *
 ********************************************************/
template<>
double Turbulent_Burning_Rate(Hexa_Block<LES3DFsd_pState, LES3DFsd_cState> *Solution_Block,
			      AdaptiveBlock3D_List &LocalSolnBlockList,
			      Grid3D_Input_Parameters &IPs) {

  double local_vol, Yf_u, rho_u, Ly, Lz, burning_rate(ZERO);
  Yf_u = 0.05518;//Fresh_Fuel_Mass_Fraction;
  rho_u = 1.13;//Fresh_Density;
  Ly = IPs.Box_Width;
  Lz = IPs.Box_Height;
  for (int p = 0 ; p <= LocalSolnBlockList.Nblk-1 ; p++ ) {
    if (LocalSolnBlockList.Block[p].used == ADAPTIVEBLOCK3D_USED) {
      for (int i = Solution_Block[p].ICl ; i <= Solution_Block[p].ICu ; i++) {
        for (int j = Solution_Block[p].JCl ; j <= Solution_Block[p].JCu ; j++) {
           for (int k = Solution_Block[p].KCl ; k <= Solution_Block[p].KCu ; k++) {
	     local_vol = Solution_Block[p].Grid.volume(i,j,k);
 	     burning_rate +=  Solution_Block[p].W[i][j][k].Fsd*local_vol*Solution_Block[p].W[i][j][k].rho; 
	   }
	}
      }
    }
  }
  burning_rate = CFFC_Summation_MPI(burning_rate);
  burning_rate = burning_rate*0.3837/(PI*0.0056*0.0056);//(Ly*Lz); //laminar_flame_speed/Ly;//(rho_u*Ly);  //(rho_u*Yf_u*Ly);

  return burning_rate;
}

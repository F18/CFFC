
#ifndef _LES3D_POLYTROPIC_INCLUDED
#include "LES3DPolytropic.h"
#endif // _LES3D_POLYTROPIC_INCLUDED

/********************************************************
 * Routine: Pre_Processing_Specializations              *
 ********************************************************/
template<>
int Hexa_Pre_Processing_Specializations(HexaSolver_Data &Data,
                                        HexaSolver_Solution_Data<LES3D_Polytropic_pState, 
                                                                 LES3D_Polytropic_cState> &Solution_Data) {

    int error_flag(0);
    
    /* -------------------------------------- *
     *                RESTART                 *
     * ---------------------------------------*/
    if (Solution_Data.Input.i_ICs == IC_RESTART) {
        
        
        /* ----------- Get turbulence statistics before computations -------- */
        
        // Get average velocity
        double u_ave, v_ave, w_ave, sqr_u;
        Time_Averaging_of_Velocity_Field(Solution_Data.Local_Solution_Blocks.Soln_Blks,
                                         Data.Local_Adaptive_Block_List,
                                         u_ave,
                                         v_ave,
                                         w_ave);
        
        // Make local velocity field blocks containing fluctuations around average velocity
        Turbulent_Velocity_Field_Multi_Block_List  Local_Velocity_Field;
        Create_Local_Velocity_Field_Multi_Block_List(Local_Velocity_Field, 
                                                     Solution_Data.Local_Solution_Blocks.Soln_Blks,
                                                     Data.Local_Adaptive_Block_List);
        Get_Local_Homogeneous_Turbulence_Velocity_Field(Solution_Data.Local_Solution_Blocks.Soln_Blks,
                                                        Data.Local_Adaptive_Block_List,
                                                        Vector3D(u_ave,v_ave,w_ave),
                                                        Local_Velocity_Field);
        
        // Collect all local velocity field blocks and put them in a global list
        Turbulent_Velocity_Field_Multi_Block_List  Global_Velocity_Field;
        Global_Velocity_Field.Create(Data.Initial_Mesh, Solution_Data.Input.Grid_IP);
        Global_Velocity_Field.Assign_Local_Velocity_Field_Blocks(Local_Velocity_Field);
        Global_Velocity_Field.Collect_Blocks_from_all_processors(Data.Octree,
                                                                 Data.Local_Adaptive_Block_List,
                                                                 Data.Global_Adaptive_Block_List);
        
        // Get the kinetic energy spectrum
        RandomFieldRogallo<LES3D_Polytropic_pState, LES3D_Polytropic_cState>  Spectrum(Solution_Data.Input);
        Spectrum.Get_Energy_Spectrum(Data.Initial_Mesh, 
                                     Data.batch_flag,
                                     Global_Velocity_Field);
        
        if(Solution_Data.Input.Turbulence_IP.rescale_spectrum == ON) {
            
            // Rescale the kinetic energy spectrum and update solution blocks
            if(!Data.batch_flag)
                cout << "\n ==============================="
                     << "\n    Rescaling Energy Spectrum   "
                     << "\n ===============================";
            Spectrum.Rescale_Energy_Spectrum_to_Initial_Spectrum();
            Spectrum.FFT_spectral_to_physical();
            Spectrum.Export_to_Velocity_Field_Blocks(Data.Initial_Mesh,Global_Velocity_Field);
            Assign_Homogeneous_Turbulence_Velocity_Field(Solution_Data.Local_Solution_Blocks.Soln_Blks,
                                                         Data.Local_Adaptive_Block_List,
                                                         Vector3D(u_ave,v_ave,w_ave),
                                                         Global_Velocity_Field);
        }


    /* -------------------------------------- *
     *             NON - RESTART              *
     * ---------------------------------------*/
    } else {
        
        /* ----------- Add turbulence to initial solution blocks -------- */

        // Make velocity field blocks where turbulence will be written
        Turbulent_Velocity_Field_Multi_Block_List  Velocity_Field;
        Velocity_Field.Create(Data.Initial_Mesh, Solution_Data.Input.Grid_IP);
        
        // Create turbulence spectrum and store in the velocity field blocks
        RandomFieldRogallo<LES3D_Polytropic_pState, LES3D_Polytropic_cState> Spectrum(Solution_Data.Input);
        error_flag = Spectrum.Create_Homogeneous_Turbulence_Velocity_Field(Data.Initial_Mesh, 
                                                                                      Data.batch_flag,
                                                                                      Velocity_Field);
        if (error_flag) {
            if (!Data.batch_flag)
                cerr << "Could not create the turbulence spectrum";
             return error_flag;
        }
        
        // Add the turbulence to the solution blocks
        Assign_Homogeneous_Turbulence_Velocity_Field(Solution_Data.Local_Solution_Blocks.Soln_Blks,
                                                     Data.Local_Adaptive_Block_List,
                                                     Velocity_Field);
        
        //if (false) { 
        if (Solution_Data.Input.Turbulence_IP.i_filter_type != FILTER_TYPE_IMPLICIT && Solution_Data.Input.Turbulence_IP.Filter_Initial_Condition) {
            Solution_Data.Explicit_Filter.Initialize(Data,Solution_Data);
            error_flag = Solution_Data.Local_Solution_Blocks.Explicitly_Filter_Initial_Condition(Solution_Data.Explicit_Filter);
            error_flag = Solution_Data.Explicit_Filter.Write_to_file();
            if (CFFC_Primary_MPI_Processor()) {
                Solution_Data.Explicit_Filter.transfer_function();
            }
        }
        
        
        
        /* -------------------- ICs Specializations --------------------- */
        error_flag = Solution_Data.Local_Solution_Blocks.ICs_Specializations(Solution_Data.Input);
        if (error_flag) return error_flag;
        
        
        
        /* ----------- Get turbulence statistics before computations -------- */
        
        // Get average velocity
        double u_ave, v_ave, w_ave, sqr_u;
        Time_Averaging_of_Velocity_Field(Solution_Data.Local_Solution_Blocks.Soln_Blks,
                                         Data.Local_Adaptive_Block_List,
                                         u_ave,
                                         v_ave,
                                         w_ave);
        
        // Make local velocity field blocks containing fluctuations around average velocity
        Turbulent_Velocity_Field_Multi_Block_List  Local_Velocity_Field;
        Create_Local_Velocity_Field_Multi_Block_List(Local_Velocity_Field, 
                                                     Solution_Data.Local_Solution_Blocks.Soln_Blks,
                                                     Data.Local_Adaptive_Block_List);
        Get_Local_Homogeneous_Turbulence_Velocity_Field(Solution_Data.Local_Solution_Blocks.Soln_Blks,
                                                        Data.Local_Adaptive_Block_List,
                                                        Vector3D(u_ave,v_ave,w_ave),
                                                        Local_Velocity_Field);
        
        // Collect all local velocity field blocks and put them in a global list
        Turbulent_Velocity_Field_Multi_Block_List  Global_Velocity_Field;
        Global_Velocity_Field.Create(Data.Initial_Mesh, Solution_Data.Input.Grid_IP);
        Global_Velocity_Field.Assign_Local_Velocity_Field_Blocks(Local_Velocity_Field);
        Global_Velocity_Field.Collect_Blocks_from_all_processors(Data.Octree,
                                                                 Data.Local_Adaptive_Block_List,
                                                                 Data.Global_Adaptive_Block_List);
        
        // Get the kinetic energy spectrum
        Spectrum.Get_Energy_Spectrum(Data.Initial_Mesh, 
                                     Data.batch_flag,
                                     Global_Velocity_Field);
        
        
        
        
        
        
    }
    return error_flag;
}

/********************************************************
 * Routine: Hexa_Post_Processing_Specializations        *
 ********************************************************/
template<>
int Hexa_Post_Processing_Specializations(HexaSolver_Data &Data,
		 	                 HexaSolver_Solution_Data<LES3D_Polytropic_pState, 
                                                                  LES3D_Polytropic_cState>&Solution_Data) {

   int error_flag(0);

    
    /* ----------- Get turbulence statistics after computations -------- */
    
    // Get average velocity
    double u_ave, v_ave, w_ave, sqr_u;
    Time_Averaging_of_Velocity_Field(Solution_Data.Local_Solution_Blocks.Soln_Blks,
                                     Data.Local_Adaptive_Block_List,
                                     u_ave,
                                     v_ave,
                                     w_ave);
    
    // Make local velocity field blocks containing fluctuations around average velocity
    Turbulent_Velocity_Field_Multi_Block_List  Local_Velocity_Field;
    Create_Local_Velocity_Field_Multi_Block_List(Local_Velocity_Field, 
                                                 Solution_Data.Local_Solution_Blocks.Soln_Blks,
                                                 Data.Local_Adaptive_Block_List);
    Get_Local_Homogeneous_Turbulence_Velocity_Field(Solution_Data.Local_Solution_Blocks.Soln_Blks,
                                                    Data.Local_Adaptive_Block_List,
                                                    Vector3D(u_ave,v_ave,w_ave),
                                                    Local_Velocity_Field);
    
    // Collect all local velocity field blocks and put them in a global list
    Turbulent_Velocity_Field_Multi_Block_List  Global_Velocity_Field;
    Global_Velocity_Field.Create(Data.Initial_Mesh, Solution_Data.Input.Grid_IP);
    Global_Velocity_Field.Assign_Local_Velocity_Field_Blocks(Local_Velocity_Field);
    Global_Velocity_Field.Collect_Blocks_from_all_processors(Data.Octree,
                                                             Data.Local_Adaptive_Block_List,
                                                             Data.Global_Adaptive_Block_List);
    
    // Get the kinetic energy spectrum
    RandomFieldRogallo<LES3D_Polytropic_pState, LES3D_Polytropic_cState>  Spectrum(Solution_Data.Input);
    Spectrum.Get_Energy_Spectrum(Data.Initial_Mesh, 
                                 Data.batch_flag,
                                 Global_Velocity_Field);

    
    
    Time_Averaging_of_Solution(Solution_Data.Local_Solution_Blocks.Soln_Blks,
                               Data.Local_Adaptive_Block_List,
                               u_ave,
                               v_ave,
                               w_ave,
                               sqr_u,
                               Data.batch_flag);
    
    SpectralAnalysis<LES3D_Polytropic_pState,LES3D_Polytropic_cState> Spectral_Analysis(Data,Solution_Data);
    typedef double (LES3D_Polytropic_pState::*member_ptr);
    typedef double (LES3D_Polytropic_pState::*member_function_ptr)(void);
    member_ptr rho_member = &LES3D_Polytropic_pState::rho;
    member_ptr p_member = &LES3D_Polytropic_pState::p;
    member_function_ptr E_member = &LES3D_Polytropic_pState::E;
    member_function_ptr ek_member = &LES3D_Polytropic_pState::ek;


    Spectral_Analysis.Get_Spectrum(rho_member,"density");
    Spectral_Analysis.Get_Spectrum(p_member,"pressure");
    Spectral_Analysis.Get_Spectrum(E_member,"total_energy");
    Spectral_Analysis.Get_Spectrum(ek_member,"kinetic_energy");

   return error_flag;

}

/*****************************************************************
 * Routine: Open_Other_Solution_Progress_Specialization_Files    *
 *****************************************************************/
template<>
int Open_Other_Solution_Progress_Specialization_Files(HexaSolver_Data &Data,
                                                      HexaSolver_Solution_Data<LES3D_Polytropic_pState,
						                               LES3D_Polytropic_cState> &Solution_Data) {
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
                                                       HexaSolver_Solution_Data<LES3D_Polytropic_pState,
						                                LES3D_Polytropic_cState> &Solution_Data) {

   int error_flag(0);

   error_flag = Close_Turbulence_Progress_File(Data.other_solution_progress_files[0]);

   return error_flag;

}

/*****************************************************************
 * Routine: Output_Other_Solution_Progress_Specialization_Data   *
 *****************************************************************/
template<>
int Output_Other_Solution_Progress_Specialization_Data(HexaSolver_Data &Data,
                                                       HexaSolver_Solution_Data<LES3D_Polytropic_pState,
						                                LES3D_Polytropic_cState> &Solution_Data) {

   int error_flag(0);
   double total_TKE, total_enstrophy, u_prime, Taylor_scale, viscosity, k_SFS;

   // Calculate various turbulence quantities
   total_TKE = Total_TKE<Hexa_Block<LES3D_Polytropic_pState, LES3D_Polytropic_cState> >(Solution_Data.Local_Solution_Blocks.Soln_Blks,
                                                                                        Data.Local_Adaptive_Block_List);

   total_enstrophy = Total_Enstrophy<Hexa_Block<LES3D_Polytropic_pState, LES3D_Polytropic_cState> >(Solution_Data.Local_Solution_Blocks.Soln_Blks,
                                                                                                    Data.Local_Adaptive_Block_List);

   u_prime = u_rms<Hexa_Block<LES3D_Polytropic_pState, LES3D_Polytropic_cState> >(Solution_Data.Local_Solution_Blocks.Soln_Blks,
                                                                                  Data.Local_Adaptive_Block_List);

   Taylor_scale = Taylor_Scale<Hexa_Block<LES3D_Polytropic_pState, LES3D_Polytropic_cState> >(Solution_Data.Local_Solution_Blocks.Soln_Blks,
                                                                                              Data.Local_Adaptive_Block_List);

   viscosity = Average_viscosity<Hexa_Block<LES3D_Polytropic_pState, LES3D_Polytropic_cState> >(Solution_Data.Local_Solution_Blocks.Soln_Blks,
                                                                                                Data.Local_Adaptive_Block_List);

   k_SFS = SFS_TKE<Hexa_Block<LES3D_Polytropic_pState, LES3D_Polytropic_cState> >(Solution_Data.Local_Solution_Blocks.Soln_Blks,
                                                                                  Data.Local_Adaptive_Block_List);

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
                                                     k_SFS);
   } /* endif */
   error_flag = CFFC_OR_MPI(error_flag);

   // Return error flag
   return error_flag;

}

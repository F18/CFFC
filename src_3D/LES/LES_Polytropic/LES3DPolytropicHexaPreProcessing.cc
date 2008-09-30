
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
        if(Solution_Data.Input.Turbulence_IP.filter_solution_before_execution == ON) {
            if(!Data.batch_flag){
                cout 
                << "\n ==============================="
                << "\n    Filtering solution first    "
                << "\n ===============================";
            }
            
            Explicit_Filters<LES3D_Polytropic_pState,LES3D_Polytropic_cState> solution_filter;
            solution_filter.Initialize(Data,Solution_Data);
            solution_filter.reset();
            Explicit_Filter_Properties::FGR = Solution_Data.Input.Turbulence_IP.FGR_secondary;
            
            error_flag = Solution_Data.Local_Solution_Blocks.Explicitly_Filter_Solution(solution_filter);
            if (error_flag) {
                cout << "\n ERROR: Could not filter solution "
                << "on processor "
                << CFFC_MPI::This_Processor_Number
                << ".\n";
                cout.flush();
            } /* endif */
            
            Solution_Data.Explicit_Filter.reset();
            Explicit_Filter_Properties::FGR = Solution_Data.Input.Turbulence_IP.FGR;

        }
        
        
        if (Solution_Data.Input.i_Flow_Type==FLOWTYPE_TURBULENT_LES) {
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
            Local_Velocity_Field.Create_Local_Velocity_Field_Multi_Block_List(Data.Initial_Mesh, 
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
            
        }

    /* -------------------------------------- *
     *             NON - RESTART              *
     * ---------------------------------------*/
    } else {
        
        if (Solution_Data.Input.i_Flow_Type==FLOWTYPE_TURBULENT_LES) {
            
            /* ----------- Add turbulence to initial solution blocks -------- */
            if ( CFFC_Primary_MPI_Processor() ) {
                                
                // Make velocity field blocks where turbulence will be written
                Data.Velocity_Field.Create(Data.Initial_Mesh, Solution_Data.Input.Grid_IP);
                   
                bool Use_Auxiliary_Mesh = false;
                if(Solution_Data.Input.Grid_IP.Mesh_Stretching==ON || Solution_Data.Input.Grid_IP.Disturb_Interior_Nodes!=OFF)
                    Use_Auxiliary_Mesh = true;
                if (!Use_Auxiliary_Mesh) {
                    // Create turbulence spectrum and store in the velocity field blocks
                    RandomFieldRogallo<LES3D_Polytropic_pState, LES3D_Polytropic_cState> Spectrum(Solution_Data.Input);
                    error_flag = Spectrum.Create_Homogeneous_Turbulence_Velocity_Field(Data.Initial_Mesh, 
                                                                                       Data.batch_flag,
                                                                                       Data.Velocity_Field);
                } else {
                    // Create a uniform single block with same dimensions as Initial_Mesh
                    cout << "\n\n Creating Auxiliary mesh " << endl;
                    Data.Auxiliary_Mesh.Create_Uniform_Initial_Grid(Solution_Data.Input.Grid_IP,Data.Initial_Mesh);
                              
                    char grid_file_name[256];
                    strcpy(grid_file_name,"Auxiliary_Mesh.dat");
                    ofstream auxiliary_mesh_file;
                    auxiliary_mesh_file.open(grid_file_name, ios::out);
                    if (auxiliary_mesh_file.bad()) return (1);
                    Data.Auxiliary_Mesh.Output_Tecplot(auxiliary_mesh_file);
                    
                    // Make velocity field blocks where turbulence will be written
                    Turbulent_Velocity_Field_Multi_Block_List  Auxiliary_Velocity_Field;
                    Auxiliary_Velocity_Field.Create(Data.Auxiliary_Mesh, Solution_Data.Input.Grid_IP);
                    
                    // Create turbulence spectrum and store in the velocity field blocks
                    RandomFieldRogallo<LES3D_Polytropic_pState, LES3D_Polytropic_cState> Spectrum(Solution_Data.Input);
                    error_flag = Spectrum.Create_Homogeneous_Turbulence_Velocity_Field(Data.Auxiliary_Mesh, 
                                                                                       Data.batch_flag,
                                                                                       Auxiliary_Velocity_Field);
                    
                    // Interpolate the Auxiliary_Velocity_Field to the Data.Velocity_Field
                    error_flag = Auxiliary_Velocity_Field.Interpolate_Turbulent_Field(Data.Initial_Mesh, Data.Velocity_Field);
                    
                }
                if (error_flag) {
                    if (!Data.batch_flag)
                        cerr << "Could not create the turbulence spectrum";
                    return error_flag;
                }
            }
                        
            CFFC_Barrier_MPI();
            Data.Velocity_Field.Broadcast();
            
            // Add the turbulence to the solution blocks
            Assign_Homogeneous_Turbulence_Velocity_Field(Solution_Data.Local_Solution_Blocks.Soln_Blks,
                                                         Data.Local_Adaptive_Block_List,
                                                         Data.Velocity_Field);
            
            
            /* ---------------------- Explicitly filter the initial condition --------------------- */
            if (Solution_Data.Input.Turbulence_IP.i_filter_type != FILTER_TYPE_IMPLICIT) {
                // Initialize the filter
                Solution_Data.Explicit_Filter.Initialize(Data,Solution_Data);

                // output the filter transfer function
                if (CFFC_Primary_MPI_Processor())
                    Solution_Data.Explicit_Filter.transfer_function(FILTER_MIDDLE_CELL);
            
                if (Solution_Data.Input.Turbulence_IP.Filter_Initial_Condition) {
                    // filter the initial condition
                    if (CFFC_Primary_MPI_Processor() && !Data.batch_flag) {
                        cout << endl;
                        cout << " ------------------------------------------------" << endl;
                        cout << "    Explicitly filtering the initial condition   " << endl;
                        cout << " ------------------------------------------------" << endl;        
                    }
                    error_flag = Solution_Data.Local_Solution_Blocks.Explicitly_Filter_Solution(Solution_Data.Explicit_Filter);
                    if (CFFC_Primary_MPI_Processor() && !Data.batch_flag) {
                        cout << "    Finished explicitly filtering the initial condition" << endl;
                    }
                    // save filter to file so don't have to recompute.
                    if (Solution_Data.Input.Turbulence_IP.i_filter_type != FILTER_TYPE_RESTART)
                        error_flag = Solution_Data.Explicit_Filter.Write_to_file();
                }
            }

            
            /* -------------------- ICs Specializations --------------------- */
            error_flag = Solution_Data.Local_Solution_Blocks.ICs_Specializations(Solution_Data.Input);
            if (error_flag) return error_flag;
            
            
            
            /* ----------- Get turbulence statistics before computations -------- */
            if (CFFC_Primary_MPI_Processor() && !Data.batch_flag) {
                cout << "    Get turbulence statistics before computations" << endl;
            }
            // Get average velocity
            double u_ave, v_ave, w_ave, sqr_u;
            Time_Averaging_of_Velocity_Field(Solution_Data.Local_Solution_Blocks.Soln_Blks,
                                             Data.Local_Adaptive_Block_List,
                                             u_ave,
                                             v_ave,
                                             w_ave);
            
            // Make local velocity field blocks containing fluctuations around average velocity
            Turbulent_Velocity_Field_Multi_Block_List  Local_Velocity_Field;
            Local_Velocity_Field.Create_Local_Velocity_Field_Multi_Block_List(Data.Initial_Mesh, 
                                                                              Solution_Data.Local_Solution_Blocks.Soln_Blks,
                                                                              Data.Local_Adaptive_Block_List);
            Get_Local_Homogeneous_Turbulence_Velocity_Field(Solution_Data.Local_Solution_Blocks.Soln_Blks,
                                                            Data.Local_Adaptive_Block_List,
                                                            Vector3D(u_ave,v_ave,w_ave),
                                                            Local_Velocity_Field);
            
            // Collect all local velocity field blocks and put them in a global list
            Data.Velocity_Field.Assign_Local_Velocity_Field_Blocks(Local_Velocity_Field);
            Data.Velocity_Field.Collect_Blocks_from_all_processors(Data.Octree,
                                                                   Data.Local_Adaptive_Block_List,
                                                                   Data.Global_Adaptive_Block_List);

            if (CFFC_Primary_MPI_Processor()) {
                bool Use_Auxiliary_Mesh = false;
                if(Solution_Data.Input.Grid_IP.Mesh_Stretching==ON || Solution_Data.Input.Grid_IP.Disturb_Interior_Nodes!=OFF)
                    Use_Auxiliary_Mesh = true;
                
                if (!Use_Auxiliary_Mesh) {
                    // Get the kinetic energy spectrum
                    RandomFieldRogallo<LES3D_Polytropic_pState, LES3D_Polytropic_cState> Spectrum(Solution_Data.Input);
                    Spectrum.Get_Energy_Spectrum(Data.Initial_Mesh, 
                                                 Data.batch_flag,
                                                 Data.Velocity_Field);
                } else {
                    // Make velocity field blocks where turbulence will be written
                    Turbulent_Velocity_Field_Multi_Block_List  Auxiliary_Velocity_Field;
                    Auxiliary_Velocity_Field.Create(Data.Auxiliary_Mesh, Solution_Data.Input.Grid_IP);
                    
                    // Interpolate the Data.Velocity_Field to the Auxiliary_Velocity_Field.
                    error_flag = Data.Velocity_Field.Interpolate_Turbulent_Field(Data.Auxiliary_Mesh, Auxiliary_Velocity_Field);
                    if (error_flag)   return error_flag;
                    
                    // Get the kinetic energy spectrum
                    RandomFieldRogallo<LES3D_Polytropic_pState, LES3D_Polytropic_cState> Spectrum(Solution_Data.Input);
                    Spectrum.Get_Energy_Spectrum(Data.Auxiliary_Mesh, 
                                                 Data.batch_flag,
                                                 Auxiliary_Velocity_Field);
                    //Deallocate Auxiliary Mesh
                    Data.Auxiliary_Mesh.Deallocate();
                }
            }
            Data.Velocity_Field.Deallocate();
        }
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

    if (Solution_Data.Input.i_Flow_Type==FLOWTYPE_TURBULENT_LES) {
        
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
        Local_Velocity_Field.Create_Local_Velocity_Field_Multi_Block_List(Data.Initial_Mesh, 
                                                                          Solution_Data.Local_Solution_Blocks.Soln_Blks,
                                                                          Data.Local_Adaptive_Block_List);
        Get_Local_Homogeneous_Turbulence_Velocity_Field(Solution_Data.Local_Solution_Blocks.Soln_Blks,
                                                        Data.Local_Adaptive_Block_List,
                                                        Vector3D(u_ave,v_ave,w_ave),
                                                        Local_Velocity_Field);
        
        // Collect all local velocity field blocks and put them in a global list
        Data.Velocity_Field.Create(Data.Initial_Mesh, Solution_Data.Input.Grid_IP);
        Data.Velocity_Field.Assign_Local_Velocity_Field_Blocks(Local_Velocity_Field);
        Data.Velocity_Field.Collect_Blocks_from_all_processors(Data.Octree,
                                                               Data.Local_Adaptive_Block_List,
                                                               Data.Global_Adaptive_Block_List);
        
        if (CFFC_Primary_MPI_Processor()) {
            bool Use_Auxiliary_Mesh = false;
            if(Solution_Data.Input.Grid_IP.Mesh_Stretching==ON || Solution_Data.Input.Grid_IP.Disturb_Interior_Nodes!=OFF)
                Use_Auxiliary_Mesh = true;
            
            if (!Use_Auxiliary_Mesh) {
                // Get the kinetic energy spectrum
                RandomFieldRogallo<LES3D_Polytropic_pState, LES3D_Polytropic_cState> Spectrum(Solution_Data.Input);
                Spectrum.Get_Energy_Spectrum(Data.Initial_Mesh, 
                                             Data.batch_flag,
                                             Data.Velocity_Field);
            } else {
                // Create a uniform single block with same dimensions as Initial_Mesh
                Data.Auxiliary_Mesh.Create_Uniform_Initial_Grid(Solution_Data.Input.Grid_IP,Data.Initial_Mesh);
                
                // Make velocity field blocks where turbulence will be written
                Turbulent_Velocity_Field_Multi_Block_List  Auxiliary_Velocity_Field;
                Auxiliary_Velocity_Field.Create(Data.Auxiliary_Mesh, Solution_Data.Input.Grid_IP);
                // Interpolate the Data.Velocity_Field to the Auxiliary_Velocity_Field.
                error_flag = Data.Velocity_Field.Interpolate_Turbulent_Field(Data.Auxiliary_Mesh, Auxiliary_Velocity_Field);
                if (error_flag)   return error_flag;
                
                // Get the kinetic energy spectrum
                RandomFieldRogallo<LES3D_Polytropic_pState, LES3D_Polytropic_cState> Spectrum(Solution_Data.Input);
                Spectrum.Get_Energy_Spectrum(Data.Auxiliary_Mesh, 
                                             Data.batch_flag,
                                             Auxiliary_Velocity_Field);
                //Deallocate Auxiliary Mesh
                Data.Auxiliary_Mesh.Deallocate();
            }
        }

        
        
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

    }
    
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

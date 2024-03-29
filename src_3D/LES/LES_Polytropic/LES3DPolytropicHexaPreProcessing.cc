
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
        if(Solution_Data.Input.ExplicitFilters_IP.Filter_Solution_Before_Execution) {
            if(!Data.batch_flag && CFFC_Primary_MPI_Processor()){
                cout 
                << "\n ==============================="
                << "\n    Filtering solution first    "
                << "\n ===============================";
                cout.flush();
            }
            
            Explicit_Filter_Commands::Initialize_Filters(Data, Solution_Data);
            Explicit_Filter_Commands::Filter_Solution(Data,Solution_Data,Explicit_Filter_Constants::SECONDARY_FILTER);
            
            if (error_flag) {
                cout << "\n ERROR: Could not filter solution "
                << "on processor "
                << CFFC_MPI::This_Processor_Number
                << ".\n";
                cout.flush();
            } /* endif */
        }
        
        if (Solution_Data.Input.i_Flow_Type==FLOWTYPE_TURBULENT_LES) {
            
            /* -------------------- ICs Specializations to generate reconstruction --------------------- */
            error_flag = Solution_Data.Local_Solution_Blocks.ICs_Specializations(Solution_Data.Input);
            if (error_flag) return error_flag;
            
            
            
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
            Local_Velocity_Field.Deallocate();
            
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
            
            Data.Velocity_Field.Deallocate();

//                if(Solution_Data.Input.Turbulence_IP.rescale_spectrum == ON) {
//                    
//                    // Rescale the kinetic energy spectrum and update solution blocks
//                    if(!Data.batch_flag)
//                        cout << "\n ==============================="
//                        << "\n    Rescaling Energy Spectrum   "
//                        << "\n ===============================";
//                    Spectrum.Rescale_Energy_Spectrum_to_Initial_Spectrum();
//                    Spectrum.FFT_spectral_to_physical();
//                    Spectrum.Export_to_Velocity_Field_Blocks(Data.Initial_Mesh,Global_Velocity_Field);
//                    Assign_Homogeneous_Turbulence_Velocity_Field(Solution_Data.Local_Solution_Blocks.Soln_Blks,
//                                                                 Data.Local_Adaptive_Block_List,
//                                                                 Vector3D(u_ave,v_ave,w_ave),
//                                                                 Global_Velocity_Field);
//                }
                
            
            
            
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
                    if(!Data.batch_flag && CFFC_Primary_MPI_Processor()){
                        cout << "\n Creating Auxiliary mesh for Turbulence generation" << endl;
                    }
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
                // Initialize the filter
                Explicit_Filter_Commands::Initialize_Filters(Data, Solution_Data);

                // output the filter transfer function
                Explicit_Filter_Commands::Transfer_Function(Data,Solution_Data);
            
                if (Solution_Data.Input.ExplicitFilters_IP.Filter_Initial_Condition[Explicit_Filter_Constants::PRIMARY_FILTER]) {
                    // filter the initial condition
                    if (CFFC_Primary_MPI_Processor() && !Data.batch_flag) {
                        cout << "\n\n Explicitly filtering the initial condition with filter 1 " << endl;
                    }
                    error_flag = Explicit_Filter_Commands::Filter_Solution(Data,Solution_Data,Explicit_Filter_Constants::PRIMARY_FILTER);
                }
                if (Solution_Data.Input.ExplicitFilters_IP.Filter_Initial_Condition[Explicit_Filter_Constants::SECONDARY_FILTER]) {
                    // filter the initial condition
                    if (CFFC_Primary_MPI_Processor() && !Data.batch_flag) {
                        cout << "\n\n Explicitly filtering the initial condition with filter 2  " << endl;
                    }
                    error_flag = Explicit_Filter_Commands::Filter_Solution(Data,Solution_Data,Explicit_Filter_Constants::SECONDARY_FILTER);
                }
            

            
            /* -------------------- ICs Specializations --------------------- */
            error_flag = Solution_Data.Local_Solution_Blocks.ICs_Specializations(Solution_Data.Input);
            if (error_flag) return error_flag;
            
            
            
            /* ----------- Get turbulence statistics before computations -------- */
            if (CFFC_Primary_MPI_Processor() && !Data.batch_flag) {
                cout << "\n\n Get turbulence statistics before computations" << endl;
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
            Local_Velocity_Field.Deallocate();

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
        Local_Velocity_Field.Deallocate();
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
        Data.Velocity_Field.Deallocate();
        
        
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

    // WARNING!!! UNNECESSARY COMPUTATION IN CASE OF CENO!!!! JUST FOR OUTPUT PURPOSE FOR NOW!
    if (Solution_Data.Input.i_Reconstruction == RECONSTRUCTION_HIGH_ORDER) {
        // Call Linear_Reconstruction_LeastSquares to set dWdx, dWdy, dWdz needed in the subroutines below
        Hexa_Block<LES3D_Polytropic_pState,LES3D_Polytropic_cState> *Soln_Blks = Solution_Data.Local_Solution_Blocks.Soln_Blks;
        for (int nBlk = 0; nBlk < Solution_Data.Local_Solution_Blocks.Number_of_Soln_Blks; nBlk++ ) {
            if (Solution_Data.Local_Solution_Blocks.Block_Used[nBlk]) {
                Soln_Blks[nBlk].Linear_Reconstruction_LeastSquares(LIMITER_ZERO); 
            }
        }
    }
    
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


#ifndef _FILTER_CONTROLLER_INCLUDED
#define _FILTER_CONTROLLER_INCLUDED

#include "Explicit_Filter_Commands.h"

#ifndef _SPECTRAL_ANALYSIS_INCLUDED
#include "../TurbulenceModelling/SpectralAnalysis.h"
#endif // _SPECTRAL_ANALYSIS_INCLUDED

template <typename Soln_pState, typename Soln_cState>
class Filter_Controller {
    
public:

    HexaSolver_Data                                     Data;  
    HexaSolver_Solution_Data<Soln_pState, Soln_cState>  Solution_Data; 

    int Initialize(char *Input_File_Name_ptr, int batch_flag) {
        Data.batch_flag=batch_flag;
        
        int error_flag;

        error_flag = Solution_Data.Get_Input_Parameters(Input_File_Name_ptr, batch_flag);  
        if(error_flag) return(error_flag);

        CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization. 

        error_flag = Initialize_Solution_Blocks(Data,
                                                Solution_Data);

        error_flag = CFFC_OR_MPI(error_flag);
        if (error_flag) return (error_flag);

        error_flag = Initial_Conditions(Data,
                                        Solution_Data);      
        if (error_flag) return (error_flag);
        
        return(0);
    }

    void Finalize(void) {
        Hexa_Post_Processing(Data,Solution_Data);
    }
    
    
    void Explicit_Filter_Operations(void){
        
        Explicit_Filter_Commands::Initialize_Filters(Data,Solution_Data);
        
        Explicit_Filter_Commands::Transfer_Function(Data,Solution_Data);
                
        //explicit_filter.Set_Filter_Property("debug_flag",ON);
        //explicit_filter_2.Set_Filter_Property("debug_flag",ON);

        

        typedef double (Soln_pState::*member_ptr);
        member_ptr rho_member = &Soln_pState::rho;
        member_ptr p_member = &Soln_pState::p;

        SpectralAnalysis<Soln_pState,Soln_cState> initial_spectrum(Data,Solution_Data), final_spectrum(Data,Solution_Data);
        //initial_spectrum.Set_Spectrum(rho_member);
        initial_spectrum.Get_Spectrum(rho_member,"density_before");

        typedef Soln_cState *** (Hexa_Block<Soln_pState,Soln_cState>::*Soln_cState_3D_ptr_type);
        Soln_cState_3D_ptr_type U_ptr = &Hexa_Block<Soln_pState,Soln_cState>::U; 
        
        //explicit_filter.filter(U_ptr);
        Explicit_Filter_Commands::Filter(rho_member,Data,Solution_Data);

        final_spectrum.Get_Spectrum_With_Reference(rho_member,"density",initial_spectrum);
        //explicit_filter.Calculate_Commutation_Error(rho_member);
        //explicit_filter.Write_to_file();

    }
    
    
    int Solver(char *Input_File_Name_ptr, int batch_flag) {
        int error_flag;
        Data.batch_flag=batch_flag;

        /******************* INPUT PARAMETERS  **********************************
         Set default values for the input solution parameters and then read user 
         specified input values from the specified input parameter file.               
         *************************************************************************/    
        error_flag = Solution_Data.Get_Input_Parameters(Input_File_Name_ptr, batch_flag);  
        if(error_flag) return(error_flag);
        
        CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
        
        /******************* SOLVER LOOP ****************************************
         Loop until command_flag is set to TERMINATE_CODE, most likely by
         Input Parameters read in Post_Processing.
         *************************************************************************/  
        while (Solution_Data.command_flag != TERMINATE_CODE) {
            
            /******************* INITIAL GRID & SOLUTION BLOCKS *********************************
             Create initial mesh and allocate solution variables for specified IBVP/BVP problem. 
             *************************************************************************************/    
            
            // New Calculation (  != CONTINUE_CODE )
            if (Solution_Data.command_flag == EXECUTE_CODE) { 
                CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization. 
                
                error_flag = Initialize_Solution_Blocks(Data,
                                                        Solution_Data);
                
                error_flag = CFFC_OR_MPI(error_flag);
                if (error_flag) return (error_flag);
                
                error_flag = Initial_Conditions(Data,
                                                Solution_Data);      
                if (error_flag) return (error_flag);
            } /* endif */
            
            //Explicit_Filter_Operations();
            
            /********************** MAIN SOLVER ****************************************
             Solve IBVP or BVP for conservation form of equations on multi-block 
             solution-adaptive quadrilateral mesh.                                  
             ****************************************************************************/    
            CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.    
            if (!Data.batch_flag) cout << "\n\n Initialization of mesh and solution data complete on " 
                << Date_And_Time() << ".";
                                   

            
            /***************************** POST PROCESSSING *******************************
             Solution calculations complete. Write 3D solution to output and restart files  
             as required, reset solution parameters, and run other cases as specified 
             by input parameters.        
             *******************************************************************************/     
            CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.    
            if (!Data.batch_flag) cout << "\n";
            
            error_flag = Hexa_Post_Processing(Data,
                                              Solution_Data);   
            if (error_flag) return (error_flag);
            
        } //END while
        
    }
    
};

#endif

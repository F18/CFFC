
#ifndef _CFFC_FILTER_CONTROLLER_INCLUDED
#define _CFFC_FILTER_CONTROLLER_INCLUDED

#include "Explicit_Filter.h"

template <typename Soln_pState, typename Soln_cState>
class CFFC_Filter_Controller {
    
public:

    HexaSolver_Data                                     Data;  
    HexaSolver_Solution_Data<Soln_pState, Soln_cState>  Solution_Data; 

    int Initialize(char *Input_File_Name_ptr, int batch_flag) {
        Data.batch_flag=batch_flag;
        
        int error_flag;

        error_flag = Solution_Data.Get_Input_Parameters(Input_File_Name_ptr, batch_flag);  
        if(error_flag) return(error_flag);

        CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.

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
        Explicit_Filters<Soln_pState,Soln_cState> explicit_filter;
        explicit_filter.Initialize(Data,Solution_Data);
        explicit_filter.transfer_function(FILTER_INNER_CELL);
        //explicit_filter.test();

        typedef double (Soln_pState::*member_ptr);
        member_ptr rho_member = &Soln_pState::rho;
        member_ptr p_member = &Soln_pState::p;


        typedef Soln_cState *** (Hexa_Block<Soln_pState,Soln_cState>::*Soln_cState_3D_ptr_type);
        Soln_cState_3D_ptr_type U_ptr = &Hexa_Block<Soln_pState,Soln_cState>::U; 
        
        //explicit_filter.filter(U_ptr);
        //explicit_filter.filter(rho_member);
        explicit_filter.Calculate_Commutation_Error(rho_member);
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
            
            Explicit_Filter_Operations();
            
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

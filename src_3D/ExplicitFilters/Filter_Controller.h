
#ifndef _FILTER_CONTROLLER_INCLUDED
#define _FILTER_CONTROLLER_INCLUDED

#include "Filter_State.h"
#include "Explicit_Filter.h"

//typedef RowVector Filter_pState;
//typedef ColumnVector Filter_cState;


template <typename Soln_pState, typename Soln_cState>
class Filter_Controller {
public:
    Input_Parameters<Soln_pState,Soln_cState>               Input;
    Grid3D_Hexa_Multi_Block_List        Initial_Mesh;
    Explicit_Filters<Soln_pState,Soln_cState> explicit_filter;

    int command_flag;
    int Get_Input_Parameters(char *Input_File_Name, int batch_flag);
    int Initialize_Solution_Blocks(void);
    int Explicit_Filter_Operations(void);
    int batch_flag;
    
    Filter_State **** Fstate;
    double **** Fdouble;
    
    
};

template <typename Soln_pState, typename Soln_cState>
int Filter_Controller<Soln_pState,Soln_cState>::Get_Input_Parameters(char *Input_File_Name, int batchflag) {
    
    int error_flag(0);
    batch_flag = batchflag;
    // The primary MPI processor processes the input parameter file.
    if (CFFC_Primary_MPI_Processor()) {
        if (!batch_flag) {
            cout << "\n Reading input data file `"<< Input_File_Name << "'.";
            cout.flush();
        } /* endif */   
        // Reads Input & sets "command_flag"
        error_flag = Input.Process_Input_Control_Parameter_File(Input_File_Name,
                                                                command_flag);
    } /* endif */
    
    // Barrier for synchronization of processors. 
    CFFC_Barrier_MPI(); 
    
    // Broadcast error flag and command flag to other MPI processors.
    CFFC_Broadcast_MPI(&error_flag, 1);
    if (error_flag != 0) return (error_flag);
    CFFC_Broadcast_MPI(&command_flag, 1);
    if (command_flag == TERMINATE_CODE) return error_flag;
    
    // Broadcast input solution parameters to other MPI processors.
    Input.Broadcast();
    
    return error_flag;
    
}

template <typename Soln_pState, typename Soln_cState>
int Filter_Controller<Soln_pState,Soln_cState>::Initialize_Solution_Blocks(void) {
    
    //NOTE: None of the initial "GRID" functions pass back an error_flag...
    int error_flag(0);
    
    /* Create the initial mesh on the primary MPI processor. */
    
    if (CFFC_Primary_MPI_Processor()) {
        Initial_Mesh.Create_Grid(Input.Grid_IP);
        Initial_Mesh.Update_Cells();
        //Outputting solution input parameters
        if (!batch_flag) {
            cout << Input << "\n";
            cout.flush(); 
        } /* endif */
    } /* endif */
    
    CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
    
    /* Broadcast the mesh to other MPI processors. */
    
    Initial_Mesh.Broadcast();
    
    /* Initialize solution blocks specializations. */
    
    //error_flag = Initialize_Solution_Blocks_Specializations(Data, Solution_Data);
    if (error_flag) return (error_flag);
    
    /* Create (allocate) list of hexahedral solution blocks on each processor. */
    
    if (!batch_flag) {
        cout << "\n Creating multi-block octree data structure and assigning"
        << "\n  solution blocks corresponding to initial mesh.";
        cout.flush();
        cout << endl;
    } /* endif */
    
    
    Fstate = new Filter_State ***[Initial_Mesh.NBlk];
    for (int n=0; n<Initial_Mesh.NBlk; n++){
        Fstate[n] = new Filter_State **[Initial_Mesh.Grid_Blks[0].NCi];
        for (int i=0; i<Initial_Mesh.Grid_Blks[0].NCi; i++) {
            Fstate[n][i] = new Filter_State *[Initial_Mesh.Grid_Blks[0].NCj];
            for (int j=0; j<Initial_Mesh.Grid_Blks[0].NCj; j++) {
                Fstate[n][i][j] = new Filter_State [Initial_Mesh.Grid_Blks[0].NCk];
            }
        }
    }

    Fdouble = new double ***[Initial_Mesh.NBlk];
    for (int n=0; n<Initial_Mesh.NBlk; n++) {
        Fdouble[n] = new double **[Initial_Mesh.Grid_Blks[n].NCi];
        for (int i=0; i<Initial_Mesh.Grid_Blks[0].NCi; i++) {
            Fdouble[n][i] = new double *[Initial_Mesh.Grid_Blks[n].NCj];
            for (int j=0; j<Initial_Mesh.Grid_Blks[0].NCj; j++) {
                Fdouble[n][i][j] = new double [Initial_Mesh.Grid_Blks[n].NCk];
                for (int k=0; k<Initial_Mesh.Grid_Blks[n].NCk; k++) {
                    Fdouble[n][i][j][k] = i*j*k;
                }
                
            }
        }        
    }
       
    // FROM AMR 
    /*error_flag = Create_Initial_Solution_Blocks<FILTER_STATE, FILTER_STATE>(Data.Initial_Mesh,
                                                                          Solution_Data.Local_Solution_Blocks,
                                                                          Solution_Data.Input,
                                                                          Data.Octree,
                                                                          Data.Global_Adaptive_Block_List,
                                                                          Data.Local_Adaptive_Block_List);
    */
     error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);
    
    /* Initialization of solution blocks complete, return. */
    
    return error_flag;
    
}


template <typename Soln_pState, typename Soln_cState>
int Filter_Controller<Soln_pState,Soln_cState>::Explicit_Filter_Operations(void) {
    explicit_filter.Initialize(Initial_Mesh,Input,batch_flag);
    explicit_filter.transfer_function(FILTER_INNER_CELL);
    
    //explicit_filter.test();
    
    typedef double (Filter_State::*member_ptr);
    member_ptr filter_member = &Filter_State::member;
    
    
    
    //typedef Soln_cState *** (Hexa_Block<Soln_pState,Soln_cState>::*Soln_cState_3D_ptr_type);
    //Soln_cState_3D_ptr_type U_ptr = &Hexa_Block<Soln_pState,Soln_cState>::U; 
    
   //explicit_filter.filter(Fdouble);
   // explicit_filter.filter(Fstate);
   // explicit_filter.Calculate_Commutation_Error(Fdouble);
    return 0;
}



/********************************************************
 * Routine: Set_Reference_Solution_States               *
 *                                                      *
 * Assigns default values to the reference solution     *
 * states.                                              *
 *                                                      *
 ********************************************************/
template<>
inline void Input_Parameters<Filter_pState, Filter_cState>::
Set_Reference_Solution_States(void) {

    Wo.member = 1.0;
    Uo.member = 2.0;
}

/********************************************************
 * Routine: Read_Reference_Solution_States              *
 *                                                      *
 * Read in the reference solution states.               *
 *                                                      *
 ********************************************************/
template<>
inline void Input_Parameters<Filter_pState, Filter_cState>::
Read_Reference_Solution_States(istream &restart_file) {

    
}

/********************************************************
 * Routine: Write_Reference_Solution_States             *
 *                                                      *
 * Write out the reference solution states.             *
 *                                                      *
 ********************************************************/
template<>
inline void Input_Parameters<Filter_pState, Filter_cState>::
Write_Reference_Solution_States(ostream &restart_file) {
    
    
}

/********************************************************
 * Routine: Deallocate_Static                           *
 *                                                      *
 * Deallocate static data of the reference solution     * 
 * states.                                              *
 *                                                      *
 ********************************************************/
template<>
inline void Input_Parameters<Filter_pState, Filter_cState>::
Deallocate_Static(void) {

    
}

//template<>
//int Hexa_Block<Filter_pState, Filter_cState>::NumVar(void) {
//    return (W[0][0][0].size());
//}

#endif

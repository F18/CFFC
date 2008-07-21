/*
 *  Explicit_Filterss.h
 *  CFFC
 *
 *  Created by Willem Deconinck on 25/03/08.
 *
 */

#ifndef _EXPLICIT_FILTER_INCLUDED
#define _EXPLICIT_FILTER_INCLUDED

#define FILTER_CFFC_MODE 1
#define FILTER_DESIGN_MODE 2

#include "../Grid/Grid3DHexaMultiBlock.h"
#include "../CFD/Input.h"
#include "Explicit_Filter_Helpers.h"
#include "General_Filter.h"
#include "Haselbacher_Filter.h"
#include "Vasilyev_Filter.h"
#include "Derivative_Reconstruction.h"


/* ------------------------------------------------------------------------------------------------------------------------------ */
/**
 * CLASS : Explicit_Filters
 */

template <typename Soln_pState, typename Soln_cState>
class Explicit_Filters {
public:
    Explicit_Filter_Properties properties;
    Explicit_Filter_Adaptor<Soln_pState,Soln_cState> adaptor;

    
    General_Filter<Soln_pState,Soln_cState> *filter_ptr;
    
    int mode;
    /* Mode 1 : group code integration */
    HexaSolver_Data *Data_ptr;
    HexaSolver_Solution_Data<Soln_pState,Soln_cState> *Solution_Data_ptr;
    
    /* Mode 2 : transfer function design */
    Grid3D_Hexa_Multi_Block_List *Grid_List_ptr;

        
    int initialized;
    
    /* ----- constructor ----- */
    
    Explicit_Filters(void) { initialized = false; }
    ~Explicit_Filters() {
        if (initialized)
            delete filter_ptr;
    }

    /* Mode 1 */
    void Initialize(HexaSolver_Data &Data, HexaSolver_Solution_Data<Soln_pState,Soln_cState> &Solution_Data);    

    /* Mode 2 */
    void Initialize(Grid3D_Hexa_Multi_Block_List &Grid_List, Input_Parameters<Soln_pState,Soln_cState> &Input, int batch_flag); 
    

    void Set_Properties(Input_Parameters<Soln_pState,Soln_cState> &IPs, int batch_flag);    
    
    void Create_filter(void);
    void transfer_function();
    void transfer_function(int flag);
    void transfer_function(int i, int j, int k);
    
    void ShowProgress(std::string message, int numIn, int maximum, int mode);
    
    template <typename Filter_Variable_Type>
    void filter(Filter_Variable_Type filter_variable);
    
    template <typename Filter_Variable_Type>
    void filter(Filter_Variable_Type filter_variable, int extra_specification);
    
    void filter_Blocks(void);
    template <typename Filter_Variable_Type>
    void filter_Blocks_design_mode(Filter_Variable_Type filter_variable);

    RowVector ***Filtered;
    void allocate_Filtered(Grid3D_Hexa_Block &Grid_Blk);
    void deallocate_Filtered(Grid3D_Hexa_Block &Grid_Blk);
    RowVector ***Divergence;
    void allocate_Divergence(Grid3D_Hexa_Block &Grid_Blk);
    void deallocate_Divergence(Grid3D_Hexa_Block &Grid_Blk);
    RowVector ***Filtered_Divergence;
    void allocate_Filtered_Divergence(Grid3D_Hexa_Block &Grid_Blk);
    void deallocate_Filtered_Divergence(Grid3D_Hexa_Block &Grid_Blk);
    RowVector ***Divergenced_Filtered;
    void allocate_Divergenced_Filtered(Grid3D_Hexa_Block &Grid_Blk);
    void deallocate_Divergenced_Filtered(Grid3D_Hexa_Block &Grid_Blk);
    RowVector ***Commutation_Error;
    void allocate_Commutation_Error(Grid3D_Hexa_Block &Grid_Blk);
    void deallocate_Commutation_Error(Grid3D_Hexa_Block &Grid_Blk);
    
    template <typename Filter_Variable_Type>
    void Calculate_Commutation_Error(Filter_Variable_Type filter_variable);
    int Calculate_Commutation_Error_Blocks(void);
    void Calculate_Commutation_Error_Block(Grid3D_Hexa_Block &Grid_Blk);
    void Output_Commutation(Grid3D_Hexa_Block &Grid_Blkf,
                            RowVector ***filtered,
                            RowVector ***divergence,
                            RowVector ***filtered_divergence,
                            RowVector ***divergence_filtered,
                            RowVector ***commutation_error,
                            int Block_Number,
                            bool Output_Title,
                            ofstream &Out_File);
    
    RowVector maxnorm(Grid3D_Hexa_Block &Grid_Blk,
                      RowVector ***Rows);
    RowVector p_norm(Grid3D_Hexa_Block &Grid_Blk,
                     RowVector ***Rows,
                     int p);
    
    void reset(void);
    void test(void);
    
    int Read_from_file(void);
    int Write_to_file(void);

};

#define FILTER_CORNER_CELL  0
#define FILTER_FACE_CELL    1
#define FILTER_EDGE_CELL    2
#define FILTER_INNER_CELL   3
#define FILTER_MIDDLE_CELL  4


template<typename Soln_pState,typename Soln_cState>
void Explicit_Filters<Soln_pState,Soln_cState>::Initialize(HexaSolver_Data &Data, HexaSolver_Solution_Data<Soln_pState,Soln_cState> &Solution_Data) {
    mode = FILTER_CFFC_MODE;
    if (!initialized) {
        Data_ptr = &Data;
        Solution_Data_ptr = &Solution_Data;

        Set_Properties(Solution_Data.Input, Data.batch_flag);
        
        if (Solution_Data.Input.i_ICs == IC_RESTART && !properties.restarted)
            properties.filter_type = FILTER_TYPE_RESTART;
        
        Create_filter();
        Solution_Data.Input.Turbulence_IP.i_filter_type = properties.filter_type;
        properties.progress_mode = Solution_Data.Input.Progress_Mode; // Don't use terminal mode when outputting to file
        initialized = true;
    }
}

template<typename Soln_pState,typename Soln_cState>
void Explicit_Filters<Soln_pState,Soln_cState>::Initialize(Grid3D_Hexa_Multi_Block_List &Grid_List, Input_Parameters<Soln_pState,Soln_cState> &Input, int batch_flag) {
    mode = FILTER_DESIGN_MODE;
    if (!initialized) {
        Grid_List_ptr  = &Grid_List;
        Set_Properties(Input,batch_flag);
        
        if (Input.i_ICs == IC_RESTART && !properties.restarted)
            properties.filter_type = FILTER_TYPE_RESTART;
        
        Create_filter();
        Input.Turbulence_IP.i_filter_type = properties.filter_type;
        properties.progress_mode = Input.Progress_Mode; // Don't use terminal mode when outputting to file
        initialized = true;
    }
}


template<typename Soln_pState,typename Soln_cState>
void Explicit_Filters<Soln_pState,Soln_cState>::Set_Properties(Input_Parameters<Soln_pState,Soln_cState> &IPs, int batch_flag) {
    properties.output_file_name = IPs.Output_File_Name_Prefix;
    properties.FGR = IPs.Turbulence_IP.FGR;
    properties.commutation_order = IPs.Turbulence_IP.commutation_order;
    properties.number_of_rings = IPs.Turbulence_IP.number_of_rings;
    properties.filter_type = IPs.Turbulence_IP.i_filter_type;
    properties.target_filter_sharpness = IPs.Turbulence_IP.Target_Filter_Sharpness;
    properties.LS_constraints = IPs.Turbulence_IP.LS_constraints;
    properties.Derivative_constraints = IPs.Turbulence_IP.Derivative_constraints;
    properties.Filter_Width_strict = IPs.Turbulence_IP.Filter_Width_strict;
    properties.Memory_Efficient = IPs.Turbulence_IP.Filter_Memory_Efficient;
    properties.batch_flag = batch_flag;
}


template <typename Soln_pState, typename Soln_cState>
void Explicit_Filters<Soln_pState,Soln_cState>::Create_filter(void) {
    int error_flag;
    switch (properties.filter_type) {
        case FILTER_TYPE_HASELBACHER:
            filter_ptr = new Haselbacher_Filter<Soln_pState,Soln_cState>;
            break;
        case FILTER_TYPE_VASILYEV:
            filter_ptr = new Vasilyev_Filter<Soln_pState,Soln_cState>;
            break;
        case FILTER_TYPE_RESTART:
            //error_flag = Read_from_file();
            if (error_flag == 1) cerr << "could not read filter_input_file" << endl;
            break;
        default:
            cerr << "Filter not defined" << endl;
            break;
    }
}

template <typename Soln_pState, typename Soln_cState>
void Explicit_Filters<Soln_pState,Soln_cState>::transfer_function(int flag) {
    if (!initialized) {
        cout << "Explicit_Filters<Soln_pState,Soln_cState> not initialized, can not return transfer_function" << endl;
        return;
    }
    int Nghost;
    int NMi, NMj, NMk;
    if (mode == FILTER_CFFC_MODE ) {
        for (int nBlk = 0; nBlk <  Solution_Data_ptr->Local_Solution_Blocks.Number_of_Soln_Blks; ++nBlk ) {
            if (Solution_Data_ptr->Local_Solution_Blocks.Block_Used[nBlk] == HEXA_BLOCK_USED){
                Nghost = Solution_Data_ptr->Local_Solution_Blocks.Soln_Blks[nBlk].Grid.Nghost;
                NMi = Solution_Data_ptr->Local_Solution_Blocks.Soln_Blks[nBlk].Grid.NCi/2;
                NMj = Solution_Data_ptr->Local_Solution_Blocks.Soln_Blks[nBlk].Grid.NCj/2;
                NMk = Solution_Data_ptr->Local_Solution_Blocks.Soln_Blks[nBlk].Grid.NCk/2;
                break;
            }
        }
    } else if (mode == FILTER_DESIGN_MODE) {
        for (int nBlk = 0; nBlk < Grid_List_ptr->NBlk; ++nBlk ) {
            Nghost = Grid_List_ptr->Grid_Blks[nBlk].Nghost;
            NMi = Grid_List_ptr->Grid_Blks[nBlk].NCi/2;
            NMj = Grid_List_ptr->Grid_Blks[nBlk].NCj/2;
            NMk = Grid_List_ptr->Grid_Blks[nBlk].NCk/2;
            break;
        }
    }
    switch (flag) {
        case FILTER_CORNER_CELL:
            transfer_function(Nghost, Nghost, Nghost);      break;
        case FILTER_FACE_CELL:
            transfer_function(Nghost, Nghost, properties.number_of_rings);     break;
        case FILTER_EDGE_CELL:
            transfer_function(Nghost, properties.number_of_rings, properties.number_of_rings);    break;
        case FILTER_INNER_CELL:
            transfer_function(properties.number_of_rings, properties.number_of_rings, properties.number_of_rings);    break;
        case FILTER_MIDDLE_CELL:
            transfer_function(NMi,NMj,NMk);    break;
    }
}

template <typename Soln_pState, typename Soln_cState>
void Explicit_Filters<Soln_pState,Soln_cState>::transfer_function(int i, int j, int k) {
    if (!initialized) {
        cout << "Explicit_Filters<Soln_pState,Soln_cState> not initialized, can not return transfer_function" << endl;
        return;
    }
    int NBlk;
    if (mode == FILTER_CFFC_MODE) {
        NBlk = Solution_Data_ptr->Local_Solution_Blocks.Number_of_Soln_Blks;
    } else if (mode == FILTER_DESIGN_MODE) {
        NBlk = Grid_List_ptr->NBlk;
    }
    if (!initialized) {
        cout << "Explicit_Filters<Soln_pState,Soln_cState> not initialized, can not return transfer_function" << endl;
        return;
    }
    cout << "\n Calculating filter transfer function for cell (" << i << "," << j << "," << k << ")."<<endl;
    for (int nBlk = 0; nBlk < NBlk; ++nBlk ) {
        if (mode == FILTER_CFFC_MODE) {
            if (Solution_Data_ptr->Local_Solution_Blocks.Block_Used[nBlk] == HEXA_BLOCK_USED){
                filter_ptr->transfer_function(Solution_Data_ptr->Local_Solution_Blocks.Soln_Blks[nBlk].Grid,
                                              Solution_Data_ptr->Local_Solution_Blocks.Soln_Blks[nBlk].Grid.Cell[i][j][k]); 
            }
        } else if (mode == FILTER_DESIGN_MODE) {
            filter_ptr->transfer_function(Grid_List_ptr->Grid_Blks[nBlk],
                                          Grid_List_ptr->Grid_Blks[nBlk].Cell[i][j][k]);            
        }
        return; // This makes sure it is called only once
    }
    
}

template <typename Soln_pState, typename Soln_cState>
void Explicit_Filters<Soln_pState,Soln_cState>::transfer_function() {
    if (!initialized) {
        cout << "Explicit_Filters<Soln_pState,Soln_cState> not initialized, can not return transfer_function" << endl;
        return;
    }
    transfer_function(FILTER_INNER_CELL);
}

template <typename Soln_pState, typename Soln_cState>
template <typename Filter_Variable_Type>
void Explicit_Filters<Soln_pState,Soln_cState>::filter(Filter_Variable_Type filter_variable) {
    if (!initialized) {
        cout << "Explicit_Filter not initialized, can not filter" << endl;
        return;
    }
    // Let an adaptor deal with what will be filtered
    adaptor.Set_Adaptor(filter_variable);
    if (mode == FILTER_CFFC_MODE) {
        filter_Blocks();
    } else if (mode == FILTER_DESIGN_MODE) {
        filter_Blocks_design_mode(filter_variable);
    }
}


template <typename Soln_pState, typename Soln_cState>
template <typename Filter_Variable_Type>
void Explicit_Filters<Soln_pState,Soln_cState>::filter(Filter_Variable_Type filter_variable, int extra_specification) {
    if (!initialized) {
        cout << "Explicit_Filter not initialized, can not filter" << endl;
        return;
    }
    // Let an adaptor deal with what will be filtered
    adaptor.Set_Adaptor(filter_variable,extra_specification);
    filter_Blocks();
}

template<typename Soln_pState, typename Soln_cState>
void Explicit_Filters<Soln_pState,Soln_cState>::filter_Blocks(void) {
    if (!initialized) {
        cout << "Explicit_Filters not initialized, can not filter" << endl;
        return;
    }
    
    Hexa_Block<Soln_pState,Soln_cState> *Soln_Blks = Solution_Data_ptr->Local_Solution_Blocks.Soln_Blks;
    int NBlk = Solution_Data_ptr->Local_Solution_Blocks.Number_of_Soln_Blks;
    int number_of_cells = 0;
    int number_of_processed_cells = 0;
    for (int nBlk = 0; nBlk < NBlk; nBlk++ ) {
        if (Solution_Data_ptr->Local_Solution_Blocks.Block_Used[nBlk] == HEXA_BLOCK_USED){
            number_of_cells += (Soln_Blks[nBlk].ICu - Soln_Blks[nBlk].ICl + 1)
            * (Soln_Blks[nBlk].JCu - Soln_Blks[nBlk].JCl + 1)
            * (Soln_Blks[nBlk].KCu - Soln_Blks[nBlk].KCl + 1);
        }
    }
    
    for (int nBlk = 0; nBlk < NBlk; ++nBlk ) {
        if (Solution_Data_ptr->Local_Solution_Blocks.Block_Used[nBlk] == HEXA_BLOCK_USED){
            adaptor.Set_Solution_Block(Soln_Blks[nBlk]); // allows access from filter to solution block
            allocate_Filtered(Soln_Blks[nBlk].Grid);
            for(int i=Soln_Blks[nBlk].ICl; i<=Soln_Blks[nBlk].ICu; i++) {
                for (int j=Soln_Blks[nBlk].JCl; j<=Soln_Blks[nBlk].JCu; j++) {
                    for (int k=Soln_Blks[nBlk].KCl; k<=Soln_Blks[nBlk].KCu; k++) {
                        Filtered[i][j][k] = filter_ptr->filter(Soln_Blks[nBlk].Grid,
                                                               Soln_Blks[nBlk].Grid.Cell[i][j][k]);
                        number_of_processed_cells++;
                        ShowProgress(" Filtering  ", number_of_processed_cells, number_of_cells, properties.progress_mode);
                    }
                }
            }
            adaptor.Load_into_Solution_Block(Filtered);
            deallocate_Filtered(Soln_Blks[nBlk].Grid);
        }
    }
    properties.progress_mode = PROGRESS_MODE_SILENT;
}

template <typename Soln_pState, typename Soln_cState>
template <typename Filter_Variable_Type>
void Explicit_Filters<Soln_pState,Soln_cState>::filter_Blocks_design_mode(Filter_Variable_Type filter_variable) {
    if (!initialized) {
        cout << "Explicit_Filters not initialized, can not filter" << endl;
        return;
    }
    
    int NBlk = Grid_List_ptr->NBlk;
    int number_of_cells = 0;
    int number_of_processed_cells = 0;
    for (int nBlk = 0; nBlk < NBlk; nBlk++ ) {
        number_of_cells += (Grid_List_ptr->Grid_Blks[nBlk].ICu - Grid_List_ptr->Grid_Blks[nBlk].ICl + 1)
                         * (Grid_List_ptr->Grid_Blks[nBlk].JCu - Grid_List_ptr->Grid_Blks[nBlk].JCl + 1)
                         * (Grid_List_ptr->Grid_Blks[nBlk].KCu - Grid_List_ptr->Grid_Blks[nBlk].KCl + 1);
    }
    
    for (int nBlk = 0; nBlk < NBlk; ++nBlk ) {
        adaptor.Set_Solution_Block(filter_variable,Grid_List_ptr->Grid_Blks[nBlk],nBlk); // allows access from filter to solution block
        allocate_Filtered(Grid_List_ptr->Grid_Blks[nBlk]);
        for(int i=Grid_List_ptr->Grid_Blks[nBlk].ICl; i<=Grid_List_ptr->Grid_Blks[nBlk].ICu; i++) {
            for (int j=Grid_List_ptr->Grid_Blks[nBlk].JCl; j<=Grid_List_ptr->Grid_Blks[nBlk].JCu; j++) {
                for (int k=Grid_List_ptr->Grid_Blks[nBlk].KCl; k<=Grid_List_ptr->Grid_Blks[nBlk].KCu; k++) {
                    Filtered[i][j][k] = filter_ptr->filter(Grid_List_ptr->Grid_Blks[nBlk],
                                                           Grid_List_ptr->Grid_Blks[nBlk].Cell[i][j][k]);
                    number_of_processed_cells++;
                    ShowProgress(" Filtering  ", number_of_processed_cells, number_of_cells, properties.progress_mode);
                }
            }
        }
        adaptor.Load_into_Solution_Block(Filtered);
        deallocate_Filtered(Grid_List_ptr->Grid_Blks[nBlk]);
    }
    properties.progress_mode = PROGRESS_MODE_SILENT;
}


template <typename Soln_pState, typename Soln_cState>
void Explicit_Filters<Soln_pState,Soln_cState>::allocate_Filtered(Grid3D_Hexa_Block &Grid_Blk){
    Filtered = new RowVector **[Grid_Blk.NCi];
    for (int i=0; i<Grid_Blk.NCi; i++) {
        Filtered[i] = new RowVector *[Grid_Blk.NCj];
        for (int j=0; j<Grid_Blk.NCj; j++) {
            Filtered[i][j] = new RowVector [Grid_Blk.NCk];
        }
    }
}

template <typename Soln_pState, typename Soln_cState>
void Explicit_Filters<Soln_pState,Soln_cState>::deallocate_Filtered(Grid3D_Hexa_Block &Grid_Blk) {
    for (int i=0; i<Grid_Blk.NCi; i++) {
        for (int j=0; j<Grid_Blk.NCj; j++) {
            delete[] Filtered[i][j];   Filtered[i][j] = NULL;
        }
        delete[] Filtered[i];   Filtered[i] = NULL;
    }
    delete[] Filtered;   Filtered = NULL;
}

template <typename Soln_pState, typename Soln_cState>
void Explicit_Filters<Soln_pState,Soln_cState>::allocate_Divergence(Grid3D_Hexa_Block &Grid_Blk){
    Divergence = new RowVector **[Grid_Blk.NCi];
    for (int i=0; i<Grid_Blk.NCi; i++) {
        Divergence[i] = new RowVector *[Grid_Blk.NCj];
        for (int j=0; j<Grid_Blk.NCj; j++) {
            Divergence[i][j] = new RowVector [Grid_Blk.NCk];
        }
    }
}

template <typename Soln_pState, typename Soln_cState>
void Explicit_Filters<Soln_pState,Soln_cState>::deallocate_Divergence(Grid3D_Hexa_Block &Grid_Blk) {
    for (int i=0; i<Grid_Blk.NCi; i++) {
        for (int j=0; j<Grid_Blk.NCj; j++) {
            delete[] Divergence[i][j];   Divergence[i][j] = NULL;
        }
        delete[] Divergence[i];   Divergence[i] = NULL;
    }
    delete[] Divergence;   Divergence = NULL;
}

template <typename Soln_pState, typename Soln_cState>
void Explicit_Filters<Soln_pState,Soln_cState>::allocate_Filtered_Divergence(Grid3D_Hexa_Block &Grid_Blk){
    Filtered_Divergence = new RowVector **[Grid_Blk.NCi];
    for (int i=0; i<Grid_Blk.NCi; i++) {
        Filtered_Divergence[i] = new RowVector *[Grid_Blk.NCj];
        for (int j=0; j<Grid_Blk.NCj; j++) {
            Filtered_Divergence[i][j] = new RowVector [Grid_Blk.NCk];
        }
    }
}

template <typename Soln_pState, typename Soln_cState>
void Explicit_Filters<Soln_pState,Soln_cState>::deallocate_Filtered_Divergence(Grid3D_Hexa_Block &Grid_Blk) {
    for (int i=0; i<Grid_Blk.NCi; i++) {
        for (int j=0; j<Grid_Blk.NCj; j++) {
            delete[] Filtered_Divergence[i][j];   Filtered_Divergence[i][j] = NULL;
        }
        delete[] Filtered_Divergence[i];   Filtered_Divergence[i] = NULL;
    }
    delete[] Filtered_Divergence;   Filtered_Divergence = NULL;
}


template <typename Soln_pState, typename Soln_cState>
void Explicit_Filters<Soln_pState,Soln_cState>::allocate_Divergenced_Filtered(Grid3D_Hexa_Block &Grid_Blk){
    Divergenced_Filtered = new RowVector **[Grid_Blk.NCi];
    for (int i=0; i<Grid_Blk.NCi; i++) {
        Divergenced_Filtered[i] = new RowVector *[Grid_Blk.NCj];
        for (int j=0; j<Grid_Blk.NCj; j++) {
            Divergenced_Filtered[i][j] = new RowVector [Grid_Blk.NCk];
        }
    }
}

template <typename Soln_pState, typename Soln_cState>
void Explicit_Filters<Soln_pState,Soln_cState>::deallocate_Divergenced_Filtered(Grid3D_Hexa_Block &Grid_Blk) {
    for (int i=0; i<Grid_Blk.NCi; i++) {
        for (int j=0; j<Grid_Blk.NCj; j++) {
            delete[] Divergenced_Filtered[i][j];   Divergenced_Filtered[i][j] = NULL;
        }
        delete[] Divergenced_Filtered[i];   Divergenced_Filtered[i] = NULL;
    }
    delete[] Divergenced_Filtered;   Divergenced_Filtered = NULL;
}

template <typename Soln_pState, typename Soln_cState>
void Explicit_Filters<Soln_pState,Soln_cState>::allocate_Commutation_Error(Grid3D_Hexa_Block &Grid_Blk){
    Commutation_Error = new RowVector **[Grid_Blk.NCi];
    for (int i=0; i<Grid_Blk.NCi; i++) {
        Commutation_Error[i] = new RowVector *[Grid_Blk.NCj];
        for (int j=0; j<Grid_Blk.NCj; j++) {
            Commutation_Error[i][j] = new RowVector [Grid_Blk.NCk];
        }
    }
}

template <typename Soln_pState, typename Soln_cState>
void Explicit_Filters<Soln_pState,Soln_cState>::deallocate_Commutation_Error(Grid3D_Hexa_Block &Grid_Blk) {
    for (int i=0; i<Grid_Blk.NCi; i++) {
        for (int j=0; j<Grid_Blk.NCj; j++) {
            delete[] Commutation_Error[i][j];   Commutation_Error[i][j] = NULL;
        }
        delete[] Commutation_Error[i];   Commutation_Error[i] = NULL;
    }
    delete[] Commutation_Error;   Commutation_Error = NULL;
}


template <typename Soln_pState, typename Soln_cState>
template <typename Filter_Variable_Type>
void Explicit_Filters<Soln_pState,Soln_cState>::Calculate_Commutation_Error(Filter_Variable_Type filter_variable) {
    if (!initialized) {
        cout << "Explicit_Filter not initialized, can not filter" << endl;
        return;
    }
    // Let an adaptor deal with what will be filtered
    adaptor.Set_Adaptor(filter_variable);
    Calculate_Commutation_Error_Blocks();
//    if (mode == FILTER_CFFC_MODE) {
//        filter_Blocks();
//    } else if (mode == FILTER_DESIGN_MODE) {
//        filter_Blocks_design_mode(filter_variable);
//    }
}

template<typename Soln_pState, typename Soln_cState>
int Explicit_Filters<Soln_pState,Soln_cState>::Calculate_Commutation_Error_Blocks(void) {
    
    char *prefix, cpu_id[256], extension[256], out_file_name[256];
    char *out_file_name_ptr;
    

    prefix = properties.output_file_name;
    strcat(prefix, "_commutation_error");
    strcpy(extension, ".dat");
    
    sprintf(cpu_id,"_cpu_%d",CFFC_MPI::This_Processor_Number);
    
    strcpy(out_file_name, prefix);
    strcat(out_file_name, cpu_id);
    strcat(out_file_name, extension);
    out_file_name_ptr = out_file_name;
    
    ofstream out_file;
    out_file.open(out_file_name, ios::out);
    if (out_file.bad()) return (1);
    
    bool first_flag = true;
        /* For every local solution block */
    if (mode == FILTER_CFFC_MODE) {
        if (!Explicit_Filter_Properties::batch_flag)
            cout << "\n\n Calculating Commutation Error: " << endl;
        Hexa_Block<Soln_pState,Soln_cState> *Soln_Blks = Solution_Data_ptr->Local_Solution_Blocks.Soln_Blks;
        for (int nBlk = 0; nBlk < Solution_Data_ptr->Local_Solution_Blocks.Number_of_Soln_Blks ; nBlk++ ) {
            if (Solution_Data_ptr->Local_Solution_Blocks.Block_Used[nBlk] == HEXA_BLOCK_USED) {
                /* ----- allocations ----- */
                allocate_Filtered(Soln_Blks[nBlk].Grid);
                allocate_Divergence(Soln_Blks[nBlk].Grid);
                allocate_Filtered_Divergence(Soln_Blks[nBlk].Grid);
                allocate_Divergenced_Filtered(Soln_Blks[nBlk].Grid);
                allocate_Commutation_Error(Soln_Blks[nBlk].Grid);
                
                /* ------ calculations -------- */
                adaptor.Set_Solution_Block(Soln_Blks[nBlk]);
                Calculate_Commutation_Error_Block(Soln_Blks[nBlk].Grid);
                
                /* ----- output commutation error ----- */
                 Output_Commutation(Soln_Blks[nBlk].Grid,
                                    Filtered,
                                    Divergence,
                                    Filtered_Divergence,
                                    Divergenced_Filtered,
                                    Commutation_Error,
                                    Data_ptr->Local_Adaptive_Block_List.Block[nBlk].info.gblknum,
                                    first_flag,
                                    out_file);
                 first_flag = false;
                 
                
                /* ----- deallocations ----- */
                deallocate_Filtered(Soln_Blks[nBlk].Grid);
                deallocate_Divergence(Soln_Blks[nBlk].Grid);
                deallocate_Filtered_Divergence(Soln_Blks[nBlk].Grid);
                deallocate_Divergenced_Filtered(Soln_Blks[nBlk].Grid);
                deallocate_Commutation_Error(Soln_Blks[nBlk].Grid);
            }         
        } 
    }         
        
    out_file.close();
    return (0);
    
}


template<typename Soln_pState, typename Soln_cState>
void Explicit_Filters<Soln_pState,Soln_cState>::Calculate_Commutation_Error_Block(Grid3D_Hexa_Block &Grid_Blk) {
    
    int number_of_cells_first = 0;
    int number_of_cells_second = 0;
    int number_of_processed_cells = 0;
    /* For every local solution block */
    
    if (mode == FILTER_CFFC_MODE) {
        Hexa_Block<Soln_pState,Soln_cState> *Soln_Blks = Solution_Data_ptr->Local_Solution_Blocks.Soln_Blks;
        for (int nBlk = 0; nBlk < Solution_Data_ptr->Local_Solution_Blocks.Number_of_Soln_Blks; nBlk++ ) {
            if (Solution_Data_ptr->Local_Solution_Blocks.Block_Used[nBlk] == HEXA_BLOCK_USED) {
                number_of_cells_first += (Soln_Blks[nBlk].ICu - Soln_Blks[nBlk].ICl + 1)
                                       * (Soln_Blks[nBlk].JCu - Soln_Blks[nBlk].JCl + 1)
                                       * (Soln_Blks[nBlk].KCu - Soln_Blks[nBlk].KCl + 1);
                number_of_cells_second += ((Soln_Blks[nBlk].ICu-properties.number_of_rings) - (Soln_Blks[nBlk].ICl+properties.number_of_rings) + 1)
                                        * ((Soln_Blks[nBlk].JCu-properties.number_of_rings) - (Soln_Blks[nBlk].JCl+properties.number_of_rings) + 1)
                                        * ((Soln_Blks[nBlk].KCu-properties.number_of_rings) - (Soln_Blks[nBlk].KCl+properties.number_of_rings) + 1);
            }
        }
    } else if (mode == FILTER_DESIGN_MODE) {
        for (int nBlk = 0; nBlk < Grid_List_ptr->NBlk; nBlk++ ) {
            number_of_cells_first += (Grid_List_ptr->Grid_Blks[nBlk].ICu - Grid_List_ptr->Grid_Blks[nBlk].ICl + 1)
                                   * (Grid_List_ptr->Grid_Blks[nBlk].JCu - Grid_List_ptr->Grid_Blks[nBlk].JCl + 1)
                                   * (Grid_List_ptr->Grid_Blks[nBlk].KCu - Grid_List_ptr->Grid_Blks[nBlk].KCl + 1);
            number_of_cells_second += ((Grid_List_ptr->Grid_Blks[nBlk].ICu-properties.number_of_rings) - (Grid_List_ptr->Grid_Blks[nBlk].ICl+properties.number_of_rings) + 1)
                                    * ((Grid_List_ptr->Grid_Blks[nBlk].JCu-properties.number_of_rings) - (Grid_List_ptr->Grid_Blks[nBlk].JCl+properties.number_of_rings) + 1)
                                    * ((Grid_List_ptr->Grid_Blks[nBlk].KCu-properties.number_of_rings) - (Grid_List_ptr->Grid_Blks[nBlk].KCl+properties.number_of_rings) + 1);
        }
    }
    

    
    /* For every cell */
    
    int temporary_adaptor_type = adaptor.adaptor_type;
    
    Derivative_Reconstruction<Soln_pState,Soln_cState> derivative_reconstructor(properties.commutation_order+1,2);
    
    /* ----- Filter ----- */
    //cout << " -- Filter " << endl;
    number_of_processed_cells = 0;
    for(int i=Grid_Blk.ICl; i<=Grid_Blk.ICu; i++) {
        for (int j=Grid_Blk.JCl; j<=Grid_Blk.JCu; j++) {
            for (int k=Grid_Blk.KCl; k<=Grid_Blk.KCu; k++) {
                Filtered[i][j][k] = filter_ptr->filter(Grid_Blk,Grid_Blk.Cell[i][j][k]);
                number_of_processed_cells++;
                ShowProgress(" -- Filter ",number_of_processed_cells,number_of_cells_first,properties.progress_mode);
            }
        }
    }
    
    /* ----- Divergence ----- */
    //cout << " -- Divergence " << endl;
    number_of_processed_cells = 0;
    for(int i=Grid_Blk.ICl; i<=Grid_Blk.ICu; i++) {
        for (int j=Grid_Blk.JCl; j<=Grid_Blk.JCu; j++) {
            for (int k=Grid_Blk.KCl; k<=Grid_Blk.KCu; k++) {
                Divergence[i][j][k] = derivative_reconstructor.divergence(Grid_Blk,Grid_Blk.Cell[i][j][k]);
                number_of_processed_cells++;
                ShowProgress(" -- Divergence ",number_of_processed_cells,number_of_cells_first,properties.progress_mode);
            }
        }
    }
    
    /* ----- Filter of Divergence ----- */
    //cout << " -- Filter of Divergence " << endl;
    number_of_processed_cells = 0;
    adaptor.Set_Commutation_RowVector(Divergence);
    for(int i=Grid_Blk.ICl+properties.number_of_rings; i<=Grid_Blk.ICu-properties.number_of_rings; i++) {
        for (int j=Grid_Blk.JCl+properties.number_of_rings; j<=Grid_Blk.JCu-properties.number_of_rings; j++) {
            for (int k=Grid_Blk.KCl+properties.number_of_rings; k<=Grid_Blk.KCu-properties.number_of_rings; k++) {
                Filtered_Divergence[i][j][k] = filter_ptr->filter(Grid_Blk,Grid_Blk.Cell[i][j][k]);
                number_of_processed_cells++;
                ShowProgress(" -- Filter of Divergence ",number_of_processed_cells,number_of_cells_second,properties.progress_mode);
                
            }
        }
    }
    
    /* ----- Divergence of Filtered ----- */
    //cout << " -- Divergence of Filtered " << endl;
    number_of_processed_cells = 0;
    adaptor.Set_Commutation_RowVector(Filtered);
    for(int i=Grid_Blk.ICl+properties.number_of_rings; i<=Grid_Blk.ICu-properties.number_of_rings; i++) {
        for (int j=Grid_Blk.JCl+properties.number_of_rings; j<=Grid_Blk.JCu-properties.number_of_rings; j++) {
            for (int k=Grid_Blk.KCl+properties.number_of_rings; k<=Grid_Blk.KCu-properties.number_of_rings; k++) {
                Divergenced_Filtered[i][j][k] = derivative_reconstructor.divergence(Grid_Blk,Grid_Blk.Cell[i][j][k]);
                number_of_processed_cells++;
                ShowProgress(" -- Divergence of Filtered ",number_of_processed_cells,number_of_cells_second,properties.progress_mode);
                Commutation_Error[i][j][k] = (RowVector(Filtered_Divergence[i][j][k] - Divergenced_Filtered[i][j][k])).absolute_values();
            }
        }
    }
    
    adaptor.Set_Adaptor(temporary_adaptor_type);
    
    /* for multiblock should be moved to other location!!! */
    RowVector Commutation_Error_maxnorm = maxnorm(Grid_Blk, Commutation_Error);
    RowVector Commutation_Error_L2norm = p_norm(Grid_Blk, Commutation_Error, 2);
    
    if (!Explicit_Filter_Properties::batch_flag){
        cout << "Filter : " << filter_ptr->filter_name() << endl;
        cout << "Commutation Order : " << properties.commutation_order << endl;
        cout << "Grid : " << (Grid_Blk.ICu - Grid_Blk.ICl + 1) 
        << "x" << (Grid_Blk.JCu - Grid_Blk.JCl + 1)
        << "x" << (Grid_Blk.KCu - Grid_Blk.KCl + 1) << endl;
        cout << "Commutation_Error:" << endl;
        cout << "   max norm = " << Commutation_Error_maxnorm;
        cout << "   L2 norm  = " << Commutation_Error_L2norm; 
    }

}


template<typename Soln_pState, typename Soln_cState>
RowVector Explicit_Filters<Soln_pState,Soln_cState>::p_norm(Grid3D_Hexa_Block &Grid_Blk,
                                                      RowVector ***Rows,
                                                      int p) {
    int imin = Grid_Blk.ICl+properties.number_of_rings,
        imax = Grid_Blk.ICu-properties.number_of_rings,
        jmin = Grid_Blk.JCl+properties.number_of_rings,
        jmax = Grid_Blk.JCu-properties.number_of_rings,
        kmin = Grid_Blk.KCl+properties.number_of_rings,
        kmax = Grid_Blk.KCu-properties.number_of_rings;
    
    
    int N = Rows[imin][jmin][jmax].size();
    
    RowVector normRow(N);
    normRow.zero();
    
    for (int k = kmin ; k <= kmax ; ++k) {
        for (int j = jmin ; j <= jmax ; ++j) {
            for (int i = imin ; i <= imax ; ++i) {
                for (int n = 0; n < N ; ++n) {
                    normRow(n) += pow(Rows[i][j][k](n),double(p));
                }
            }
        }
    }
    
    for (int n = 0; n < N ; ++n) {
        normRow(n) = pow(normRow(n),ONE/double(p));
    }
    return normRow;
}

template<typename Soln_pState, typename Soln_cState>
RowVector Explicit_Filters<Soln_pState,Soln_cState>::maxnorm(Grid3D_Hexa_Block &Grid_Blk,
                                                             RowVector ***Rows) {
    int imin = Grid_Blk.ICl+properties.number_of_rings,
        imax = Grid_Blk.ICu-properties.number_of_rings,
        jmin = Grid_Blk.JCl+properties.number_of_rings,
        jmax = Grid_Blk.JCu-properties.number_of_rings,
        kmin = Grid_Blk.KCl+properties.number_of_rings,
        kmax = Grid_Blk.KCu-properties.number_of_rings;
    
    int N = Rows[imin][jmin][jmax].size();
    
    RowVector normRow(N);
    normRow.zero();
    
    for (int k = kmin ; k <= kmax ; ++k) {
        for (int j = jmin ; j <= jmax ; ++j) {
            for (int i = imin ; i <= imax ; ++i) {
                for (int n = 0; n < N ; ++n) {
                    normRow(n) = max(normRow(n),Rows[i][j][k](n));
                }
            }
        }
    }
    
    return normRow;
}







template<typename Soln_pState, typename Soln_cState>
void Explicit_Filters<Soln_pState,Soln_cState>::reset(void) {
    if (!initialized) {
        cout << "Explicit_Filter not initialized, can not reset" << endl;
        return;
    }
    if (mode == FILTER_CFFC_MODE) {
        int NBlk;
        NBlk = Solution_Data_ptr->Local_Solution_Blocks.Number_of_Soln_Blks;
        for (int nBlk = 0; nBlk < NBlk; ++nBlk ) {
            if (Solution_Data_ptr->Local_Solution_Blocks.Block_Used[nBlk] == HEXA_BLOCK_USED){
                filter_ptr->Reset_Filter_Weights(Solution_Data_ptr->Local_Solution_Blocks.Soln_Blks[nBlk].Grid); 
            }
        }
    } else if (mode == FILTER_DESIGN_MODE) {
        int NBlk;
        NBlk = Grid_List_ptr->NBlk;
        for (int nBlk = 0; nBlk < NBlk; ++nBlk ) {
            filter_ptr->Reset_Filter_Weights(Grid_List_ptr->Grid_Blks[nBlk]); 
        }
    }
}

template<typename Soln_pState, typename Soln_cState>
void Explicit_Filters<Soln_pState,Soln_cState>::test(void) {
    if (!initialized) {
        cout << "Explicit_Filter not initialized, can not test" << endl;
        return;
    }
    if (mode == FILTER_CFFC_MODE) {
        int NBlk;
        NBlk = Solution_Data_ptr->Local_Solution_Blocks.Number_of_Soln_Blks;
        for (int nBlk = 0; nBlk < NBlk; ++nBlk ) {
            if (Solution_Data_ptr->Local_Solution_Blocks.Block_Used[nBlk] == HEXA_BLOCK_USED){
                filter_ptr->filter_tests(Solution_Data_ptr->Local_Solution_Blocks.Soln_Blks[nBlk].Grid,
                                         Solution_Data_ptr->Local_Solution_Blocks.Soln_Blks[nBlk].Grid.Cell[12][12][12]); 
            }
        }
    } else if (mode == FILTER_DESIGN_MODE) {
        int NBlk;
        NBlk = Grid_List_ptr->NBlk;
        for (int nBlk = 0; nBlk < NBlk; ++nBlk ) {
            filter_ptr->filter_tests(Grid_List_ptr->Grid_Blks[nBlk],
                                     Grid_List_ptr->Grid_Blks[nBlk].Cell[12][12][12]); 
        }
    }
    
}


template<typename Soln_pState, typename Soln_cState>
int Explicit_Filters<Soln_pState,Soln_cState>::Write_to_file(void){
    if (properties.Memory_Efficient) {
        cerr <<"\n Memory Efficient mode is ON --> No filter weights to write to file. " << endl;
        return (0);
    }
    if (!Explicit_Filter_Properties::batch_flag)
        cout << "\n Writing explicit filter coefficients to file." << endl ;
    int i;
    char *prefix, cpu_id[256], extension[256], out_file_name[256];
    char *out_file_name_ptr;
    
    prefix = properties.output_file_name;
    strcat(prefix, "_explicit_filter");
    strcpy(extension, ".dat");
    
    sprintf(cpu_id,"_cpu_%d",CFFC_MPI::This_Processor_Number);
    
    strcpy(out_file_name, prefix);
    strcat(out_file_name, cpu_id);
    strcat(out_file_name, extension);
    out_file_name_ptr = out_file_name;
    
    ofstream out_file;
    out_file.open(out_file_name, ios::out);
    if (out_file.bad()) return (1);
    
    out_file << filter_ptr->filter_type() << "\n";                    
    
    if (mode == FILTER_CFFC_MODE) {
        int NBlk;
        NBlk = Solution_Data_ptr->Local_Solution_Blocks.Number_of_Soln_Blks;
        for (int nBlk = 0; nBlk < NBlk; ++nBlk ) {
            if (Solution_Data_ptr->Local_Solution_Blocks.Block_Used[nBlk] == HEXA_BLOCK_USED){
                filter_ptr->Write_to_file(Solution_Data_ptr->Local_Solution_Blocks.Soln_Blks[nBlk].Grid,out_file);
            }         
        } 
    } else if (mode == FILTER_DESIGN_MODE) {
        int NBlk;
        NBlk = Grid_List_ptr->NBlk;
        for (int nBlk = 0; nBlk < NBlk; ++nBlk ) {
            filter_ptr->Write_to_file(Grid_List_ptr->Grid_Blks[nBlk],out_file);
        }
    }
    
    out_file.close();
    return (0);
}




template<typename Soln_pState, typename Soln_cState>
int Explicit_Filters<Soln_pState,Soln_cState>::Read_from_file(void) {
    if (!Explicit_Filter_Properties::batch_flag)
        cout << "\n Reading explicit filter coefficients from file." << endl ;
    int i;
    char *prefix, cpu_id[256], extension[256], in_file_name[256];
    char *in_file_name_ptr;
    
    prefix = properties.output_file_name;
    strcat(prefix, "_explicit_filter");
    strcpy(extension, ".dat");
    
    sprintf(cpu_id,"_cpu_%d",CFFC_MPI::This_Processor_Number);
    
    strcpy(in_file_name, prefix);
    strcat(in_file_name, cpu_id);
    strcat(in_file_name, extension);
    in_file_name_ptr = in_file_name;
    
    
    ifstream in_file;
    in_file.open(in_file_name, ios::in);
    if (in_file.bad()) return (1);        
    in_file.setf(ios::skipws);
    
    in_file >> properties.filter_type;
    Create_filter();
    
    if (mode == FILTER_CFFC_MODE) {
        int NBlk;
        NBlk = Solution_Data_ptr->Local_Solution_Blocks.Number_of_Soln_Blks;
        for (int nBlk = 0; nBlk < NBlk; ++nBlk ) {
            if (Solution_Data_ptr->Local_Solution_Blocks.Block_Used[nBlk] == HEXA_BLOCK_USED){
                    filter_ptr->Read_from_file(Solution_Data_ptr->Local_Solution_Blocks.Soln_Blks[nBlk].Grid,in_file);
            }         
        } 
    } else if (mode == FILTER_DESIGN_MODE) {
        int NBlk;
        NBlk = Grid_List_ptr->NBlk;
        for (int nBlk = 0; nBlk < NBlk; ++nBlk ) {
            filter_ptr->Read_from_file(Grid_List_ptr->Grid_Blks[nBlk],in_file);
        }
    }
    
    in_file.unsetf(ios::skipws);
    in_file.close();
    properties.restarted = true;
    return (0);
}


template<typename Soln_pState, typename Soln_cState>
void Explicit_Filters<Soln_pState,Soln_cState>::Output_Commutation(Grid3D_Hexa_Block &Grid_Blk,
                                                             RowVector ***filtered,
                                                             RowVector ***divergence,
                                                             RowVector ***filtered_divergence,
                                                             RowVector ***divergence_filtered,
                                                             RowVector ***commutation_error,
                                                             int Block_Number,
                                                             bool Output_Title,
                                                             ofstream &Out_File){
    
    
    
    int nVar = 0;
    
    
    int imin = Grid_Blk.ICl+properties.number_of_rings,
    imax = Grid_Blk.ICu-properties.number_of_rings,
    jmin = Grid_Blk.JCl+properties.number_of_rings,
    jmax = Grid_Blk.JCu-properties.number_of_rings,
    kmin = Grid_Blk.KCl+properties.number_of_rings,
    kmax = Grid_Blk.KCu-properties.number_of_rings;
    
    /* Ensure boundary conditions are updated before
     evaluating solution at the nodes. */
    
    //  BCs(IPs);
    
    /* Output nodal solution data. */
    
    Out_File << setprecision(14);
    
    if (Output_Title) {
        
        Out_File << "TITLE = \"" << CFFC_Name() << ": Commutation error "
        << "\"" << "\n"
        << "VARIABLES = \"x\" \\ \n"
        << "\"y\" \\ \n"
        << "\"z\" \\ \n"
        << "\"unfiltered\" \\ \n"
        << "\"filtered\" \\ \n"
        << "\"divergence\" \\ \n"
        << "\"filtered of divergence\" \\ \n"
        << "\"divergence of filtered\" \\ \n"
        << "\"commutation error\" \\ \n";
        Out_File << "DATASETAUXDATA filter = \"" << filter_ptr->filter_name() << "\" \\ \n"
                 << "DATASETAUXDATA commutation_order = \"" << properties.commutation_order << "\" \\ \n"
                 << "DATASETAUXDATA filter_grid_ratio = \"" << properties.FGR << "\" \\ \n"
                 << "DATASETAUXDATA number_of_rings = \"" << properties.number_of_rings << "\" \\ \n";
        Out_File<< "ZONE T =  \"Block Number = " << Block_Number
        << "\" \\ \n"
        << "I = " << imax - imin + 1 << " \\ \n"
        << "J = " << jmax - jmin + 1 << " \\ \n"
        << "K = " << kmax - kmin + 1 << " \\ \n"
        << "AUXDATA Commutation_Error_max_norm = \"" << maxnorm(Grid_Blk,Commutation_Error)(nVar) << "\" \\ \n"
        << "AUXDATA Commutation_Error_L2_norm = \"" << p_norm(Grid_Blk,Commutation_Error,2)(nVar) << "\" \\ \n"
        << "DATAPACKING = POINT \n";
    } else {
        Out_File << "ZONE T =  \"Block Number = " << Block_Number
        << "\" \\ \n"
        << "I = " << imax - imin + 1 << " \\ \n"
        << "J = " << jmax - jmin + 1 << " \\ \n"
        << "K = " << kmax - kmin + 1 << " \\ \n"
        << "AUXDATA Commutation_Error_max_norm = \"" << maxnorm(Grid_Blk,Commutation_Error)(nVar) << "\" \\ \n"
        << "AUXDATA Commutation_Error_L2_norm = \"" << p_norm(Grid_Blk,Commutation_Error,2)(nVar) << "\" \\ \n"
        << "DATAPACKING = POINT \n";              
    } /* endif */
    
    for (int k = kmin ; k <= kmax ; ++k) {
        for (int j = jmin ; j <= jmax ; ++j) {
            for (int i = imin ; i <= imax ; ++i) {
                Out_File << " "  << Grid_Blk.Cell[i][j][k].Xc;
                Out_File.setf(ios::scientific);
                Out_File    << " " << adaptor.FilterVariable(i,j,k)(nVar)
                            << " " << filtered[i][j][k](nVar)
                            << " " << divergence[i][j][k](nVar)
                            << " " << filtered_divergence[i][j][k](nVar)
                            << " " << divergence_filtered[i][j][k](nVar)
                            << " " << commutation_error[i][j][k](nVar) << "\n";
                Out_File.unsetf(ios::scientific);
            } /* endfor */
        } /* endfor */
    } /* endfor */
    
    Out_File << setprecision(6);
    
}





template <typename Soln_pState, typename Soln_cState>
void Explicit_Filters<Soln_pState,Soln_cState>::ShowProgress(std::string message, int numIn, int maximum, int mode) {
    if (Explicit_Filter_Properties::batch_flag)
        mode = PROGRESS_MODE_SILENT;
    int first_index = 1;
    int last_index = maximum;
    int percent = int(100*numIn/double(last_index));
    if (mode == PROGRESS_MODE_TERMINAL) {
        char barspin[16] = {'\\','\\','\\','\\',
            '|', '|','|','|',
            '/','/', '/', '/',
        '-','-','-','-'};
        
        int whichOne;
        
        whichOne = numIn % 16;
        
        std::cout << '\r'
        << message << setw(3)  <<  percent << " %"
        << "  " << barspin[whichOne] << " ";
        std::cout.flush();
        if (percent == 100) {
            std::cout << '\r'
            << message << setw(3)  <<  percent << " %      " << std::endl;
        }
        
    } else if (mode == PROGRESS_MODE_FILE) {
        if (numIn == first_index) {
            std::cout << message << endl;
        }
        int previous_percent = int(100*(numIn-1)/double(last_index));
        if (percent != previous_percent || numIn == first_index) {
            std::cout << " " << percent << "%";
            std::cout.flush();
        }
        if (percent == 100) {
            std::cout << std::endl;
        }
        
    } else if (mode == PROGRESS_MODE_MESSAGE) {
        // no progress, just message
        if (numIn == first_index) {
            std::cout << message << endl;
        }
    } else {
        // nothing
    }
    return;
}


#endif

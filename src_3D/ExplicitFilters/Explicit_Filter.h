/*
 *  Explicit_Filterss.h
 *  CFFC
 *
 *  Created by Willem Deconinck on 25/03/08.
 *
 */

#ifndef _EXPLICIT_FILTER_INCLUDED
#define _EXPLICIT_FILTER_INCLUDED

#include "../Grid/Grid3DHexaMultiBlock.h"
#include "../CFD/Input.h"
#include "Explicit_Filter_Helpers.h"
#include "General_Filter.h"
#include "Haselbacher_Filter.h"
#include "Vasilyev_Filter.h"
//#include "Derivative_Reconstruction.h"
//#include "Tophat_Filter.h"
//#include "Gaussian_Filter.h"
#include "Finite_Difference.h"
#include "Explicit_Filter_Constants.h"


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
    
    Hexa_Block<Soln_pState,Soln_cState> *SolnBlk_ptr;
    Input_Parameters<Soln_pState,Soln_cState> *Input_ptr;
    
    
    int initialized;
    
    /* ----- constructor ----- */
    
    //tmp hack for debugging
    Explicit_Filters() : initialized(false) {    }
    
    Explicit_Filters(Hexa_Block<Soln_pState,Soln_cState> *SolnBlk) :
        initialized(false),
        Filter_Weights_Allocated(false),
        Derivative_Reconstruction_Weights_Allocated(false),
        SolnBlk_ptr(SolnBlk) { }
    
    /* ----- destructor ----- */
    ~Explicit_Filters() {
        if (initialized)
            delete filter_ptr;
    }

    void Initialize(Explicit_Filter_Constants::Filter_Number& filter_number, int batch_flag, Input_Parameters<Soln_pState,Soln_cState> &Input);
    void Initialize(int batch_flag, Input_Parameters<Soln_pState,Soln_cState> &Input);    
    void Initialize_Secondary(int batch_flag, Input_Parameters<Soln_pState,Soln_cState> &Input);    

    
    void Set_Properties(Input_Parameters<Soln_pState,Soln_cState> &IPs, int batch_flag);
    
    template <typename T>
    void Set_Filter_Property(string property_name, T property_value);  
    
    template <typename T>
    void Set_Operating_Property(string property_name, T property_value);   

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
    RowVector ***Truncation_Error;
    void allocate_Truncation_Error(Grid3D_Hexa_Block &Grid_Blk);
    void deallocate_Truncation_Error(Grid3D_Hexa_Block &Grid_Blk);
    
    RowVector      ***Filter_Weights;  // Weights used in discrete explicit filtering
    bool  ***Filter_Weights_Assigned;  // Stores if the filterweights have been allocated
    bool    Filter_Weights_Allocated;
    void Allocate_Filter_Weights(void);
    void Deallocate_Filter_Weights(void);
    
    DenseMatrix    ***Derivative_Reconstruction_Weights;  // Weights used in reconstructing derivatives (in LES_Filters.h)
    bool  ***Derivative_Reconstruction_Weights_Assigned;  // Stores if the Derivative_Reconstruction_Weights have been allocated
    bool    Derivative_Reconstruction_Weights_Allocated;    
    void Allocate_Derivative_Reconstruction_Weights(void);
    void Deallocate_Derivative_Reconstruction_Weights(void);
        
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
                            RowVector ***truncation_error,
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


template<typename Soln_pState,typename Soln_cState>
void Explicit_Filters<Soln_pState,Soln_cState>::Initialize(Explicit_Filter_Constants::Filter_Number& filter_number, int batch_flag, Input_Parameters<Soln_pState,Soln_cState> &Input) {

// todo    
}

template<typename Soln_pState,typename Soln_cState>
void Explicit_Filters<Soln_pState,Soln_cState>::Initialize(int batch_flag, Input_Parameters<Soln_pState,Soln_cState> &Input) {
    if (!initialized) {
        
        Input_ptr = &Input;
        
        Set_Properties(Input, batch_flag);
        
        /*if (Solution_Data.Input.i_ICs == IC_RESTART && !properties.restarted)
            properties.filter_type = Explicit_Filter_Constants::RESTART_FILTER; */
        
        Create_filter();
        initialized = true;
    }
}

template<typename Soln_pState,typename Soln_cState>
void Explicit_Filters<Soln_pState,Soln_cState>::Initialize_Secondary(int batch_flag, Input_Parameters<Soln_pState,Soln_cState> &Input) {
    if (!initialized) {
        
        Input_ptr = &Input;

        Set_Properties(Input, batch_flag);
        
        /*if (Solution_Data.Input.i_ICs == IC_RESTART && !properties.restarted)
         properties.filter_type = Explicit_Filter_Constants::RESTART_FILTER; */
        
        //properties.Set_Operating_Property("memory_efficient",ON);
        string output_file_name;
        properties.Get_Property(output_file_name,"output_file_name");
        properties.Set_Operating_Property("output_file_name",output_file_name+string("_secondary_filter"));
        properties.Set_Filter_Property("FGR",Input.Turbulence_IP.FGR_secondary);
        if (properties.Get_Property_int("use_fixed_filter_width")) {
            properties.Set_Filter_Property("fixed_filter_width",Input.Turbulence_IP.Filter_Width_secondary);
        }
        properties.Set_Operating_Property("progress_mode",Input.Progress_Mode);

        properties.Set_Filter_Property("filter_type",Input.Turbulence_IP.i_filter_type_secondary);
        Create_filter();
        initialized = true;
    }
}



template<typename Soln_pState,typename Soln_cState>
void Explicit_Filters<Soln_pState,Soln_cState>::Set_Properties(Input_Parameters<Soln_pState,Soln_cState> &Input, int batch_flag) {
    properties.Set_Properties(Input,batch_flag);
}

template<typename Soln_pState,typename Soln_cState>
template<typename T>
void Explicit_Filters<Soln_pState,Soln_cState>::Set_Filter_Property(string property_name, T property_value){
    properties.Set_Filter_Property(property_name,property_value);
}

template<typename Soln_pState,typename Soln_cState>
template<typename T>
void Explicit_Filters<Soln_pState,Soln_cState>::Set_Operating_Property(string property_name, T property_value){
    properties.Set_Operating_Property(property_name,property_value);
}

template <typename Soln_pState, typename Soln_cState>
void Explicit_Filters<Soln_pState,Soln_cState>::Create_filter(void) {
    int error_flag;
    switch (properties.Get_Property_int("filter_type")) {
        case Explicit_Filter_Constants::HASELBACHER_FILTER:
            filter_ptr = new Haselbacher_Filter<Soln_pState,Soln_cState>(properties);
            break;
        case Explicit_Filter_Constants::VASILYEV_FILTER:
            filter_ptr = new Vasilyev_Filter<Soln_pState,Soln_cState>(properties);
            break;
        case Explicit_Filter_Constants::RESTART_FILTER:
            error_flag = Read_from_file();
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
    int Nghost, number_of_rings;
    int NMi, NMj, NMk;

    Nghost = SolnBlk_ptr->Grid.Nghost;
    NMi    = SolnBlk_ptr->Grid.NCi/2;
    NMj    = SolnBlk_ptr->Grid.NCj/2;
    NMk    = SolnBlk_ptr->Grid.NCk/2;

    number_of_rings = properties.Get_Property_int("number_of_rings");
    switch (flag) {
        case Explicit_Filter_Constants::CORNER_CELL:
            transfer_function(Nghost, Nghost, Nghost);      break;
        case Explicit_Filter_Constants::FACE_CELL:
            transfer_function(Nghost, Nghost, number_of_rings);     break;
        case Explicit_Filter_Constants::EDGE_CELL:
            transfer_function(Nghost, number_of_rings, number_of_rings);    break;
        case Explicit_Filter_Constants::INNER_CELL:
            transfer_function(max(Nghost,number_of_rings), max(Nghost,number_of_rings), max(Nghost,number_of_rings));    break;
        case Explicit_Filter_Constants::MIDDLE_CELL:
            transfer_function(NMi,NMj,NMk);    break;
    }
}

template <typename Soln_pState, typename Soln_cState>
void Explicit_Filters<Soln_pState,Soln_cState>::transfer_function(int i, int j, int k) {
    if (!initialized) {
        cout << "Explicit_Filters<Soln_pState,Soln_cState> not initialized, can not return transfer_function" << endl;
        return;
    }

    filter_ptr->transfer_function(SolnBlk_ptr->Grid,
                                  SolnBlk_ptr->Grid.Cell[i][j][k]);     
}

template <typename Soln_pState, typename Soln_cState>
void Explicit_Filters<Soln_pState,Soln_cState>::transfer_function() {
    if (!initialized) {
        cout << "Explicit_Filters<Soln_pState,Soln_cState> not initialized, can not return transfer_function" << endl;
        return;
    }
    transfer_function(Explicit_Filter_Constants::INNER_CELL);
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
    filter_Blocks();

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
    int progress_mode = properties.Get_Property_int("progress_mode");
    int number_of_cells = 0;
    int number_of_processed_cells = 0;

    number_of_cells = (SolnBlk_ptr->ICu - SolnBlk_ptr->ICl + 1)
    * (SolnBlk_ptr->JCu - SolnBlk_ptr->JCl + 1)
    * (SolnBlk_ptr->KCu - SolnBlk_ptr->KCl + 1);

    adaptor.Set_Solution_Block(*SolnBlk_ptr); // allows access from filter to solution block
    allocate_Filtered(SolnBlk_ptr->Grid);
    for(int i=SolnBlk_ptr->ICl; i<=SolnBlk_ptr->ICu; i++) {
        for (int j=SolnBlk_ptr->JCl; j<=SolnBlk_ptr->JCu; j++) {
            for (int k=SolnBlk_ptr->KCl; k<=SolnBlk_ptr->KCu; k++) {
                Filtered[i][j][k] = filter_ptr->filter(*this,SolnBlk_ptr->Grid,
                                                       SolnBlk_ptr->Grid.Cell[i][j][k]);
                number_of_processed_cells++;
                ShowProgress(" Filtering  ", number_of_processed_cells, number_of_cells, progress_mode);
            }
        }
    }
    adaptor.Load_into_Solution_Block(Filtered);
    deallocate_Filtered(SolnBlk_ptr->Grid);
    
    if (progress_mode!=PROGRESS_MODE_SILENT){
        properties.Set_Operating_Property("progress_mode",PROGRESS_MODE_SILENT);
    }
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
void Explicit_Filters<Soln_pState,Soln_cState>::allocate_Truncation_Error(Grid3D_Hexa_Block &Grid_Blk){
    Truncation_Error = new RowVector **[Grid_Blk.NCi];
    for (int i=0; i<Grid_Blk.NCi; i++) {
        Truncation_Error[i] = new RowVector *[Grid_Blk.NCj];
        for (int j=0; j<Grid_Blk.NCj; j++) {
            Truncation_Error[i][j] = new RowVector [Grid_Blk.NCk];
        }
    }
}

template <typename Soln_pState, typename Soln_cState>
void Explicit_Filters<Soln_pState,Soln_cState>::deallocate_Truncation_Error(Grid3D_Hexa_Block &Grid_Blk) {
    for (int i=0; i<Grid_Blk.NCi; i++) {
        for (int j=0; j<Grid_Blk.NCj; j++) {
            delete[] Truncation_Error[i][j];   Truncation_Error[i][j] = NULL;
        }
        delete[] Truncation_Error[i];   Truncation_Error[i] = NULL;
    }
    delete[] Truncation_Error;   Truncation_Error = NULL;
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

}

template<typename Soln_pState, typename Soln_cState>
int Explicit_Filters<Soln_pState,Soln_cState>::Calculate_Commutation_Error_Blocks(void) {
    
    char prefix[256], cpu_id[256], extension[256], out_file_name[256];
    char *out_file_name_ptr;
    
    string prefix_string;
    properties.Get_Property(prefix_string,"output_file_name");
    strcpy(prefix,prefix_string.c_str());
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

    if (!properties.Get_Property_int("batch_flag"))
        cout << "\n\n Calculating Commutation Error: " << endl;

    /* ----- allocations ----- */
    allocate_Filtered(SolnBlk_ptr->Grid);
    allocate_Divergence(SolnBlk_ptr->Grid);
    allocate_Filtered_Divergence(SolnBlk_ptr->Grid);
    allocate_Divergenced_Filtered(SolnBlk_ptr->Grid);
    allocate_Commutation_Error(SolnBlk_ptr->Grid);
    allocate_Truncation_Error(SolnBlk_ptr->Grid);
    
    /* ------ calculations -------- */
    adaptor.Set_Solution_Block(*SolnBlk_ptr);
    adaptor.Set_Initial_Condition(*Input_ptr);
    Calculate_Commutation_Error_Block(SolnBlk_ptr->Grid);
    
    /* ----- output commutation error ----- */
    // if first_flag = true --> Make Tecplot header
    bool first_flag = true;
    int gblknum = 0;
    Output_Commutation(SolnBlk_ptr->Grid,
                       Filtered,
                       Divergence,
                       Filtered_Divergence,
                       Divergenced_Filtered,
                       Commutation_Error,
                       Truncation_Error,
                       gblknum,
                       first_flag,
                       out_file);
    
    /* ----- deallocations ----- */
    deallocate_Filtered(SolnBlk_ptr->Grid);
    deallocate_Divergence(SolnBlk_ptr->Grid);
    deallocate_Filtered_Divergence(SolnBlk_ptr->Grid);
    deallocate_Divergenced_Filtered(SolnBlk_ptr->Grid);
    deallocate_Commutation_Error(SolnBlk_ptr->Grid);
    deallocate_Truncation_Error(SolnBlk_ptr->Grid);
        
    out_file.close();
    return (0);
    
}


template<typename Soln_pState, typename Soln_cState>
void Explicit_Filters<Soln_pState,Soln_cState>::Calculate_Commutation_Error_Block(Grid3D_Hexa_Block &Grid_Blk) {
    
    //properties.number_of_rings_increased = properties.number_of_rings;
    //Derivative_Reconstruction<Soln_pState,Soln_cState> derivative_reconstructor(properties.commutation_order,properties.number_of_rings_increased);
    
    properties.Set_Operating_Property("derivative_accuracy",properties.Get_Property_int("finite_differencing_order"));
    Finite_Difference_Class<Soln_pState,Soln_cState> finite_differencer(properties.Get_Property_int("derivative_accuracy"));
    properties.Set_Operating_Property("number_of_rings_increased",finite_differencer.Get_central_rings());
    
    //assert(Grid_Blk.Nghost >= properties.number_of_rings+properties.number_of_rings_increased);
    int number_of_rings = properties.Get_Property_int("number_of_rings");
    int number_of_rings_increased = properties.Get_Property_int("number_of_rings_increased");
    int progress_mode = properties.Get_Property_int("progress_mode");
    int number_of_cells_first = 0;
    int number_of_cells_second = 0;
    int number_of_cells_third = 0;
    int number_of_processed_cells = 0;
    
    number_of_cells_first += (SolnBlk_ptr->ICu - SolnBlk_ptr->ICl + 1)
                           * (SolnBlk_ptr->JCu - SolnBlk_ptr->JCl + 1)
                           * (SolnBlk_ptr->KCu - SolnBlk_ptr->KCl + 1);
    number_of_cells_second += ((SolnBlk_ptr->ICu-number_of_rings) - (SolnBlk_ptr->ICl+number_of_rings) + 1)
                            * ((SolnBlk_ptr->JCu-number_of_rings) - (SolnBlk_ptr->JCl+number_of_rings) + 1)
                            * ((SolnBlk_ptr->KCu-number_of_rings) - (SolnBlk_ptr->KCl+number_of_rings) + 1);
    number_of_cells_third += ((SolnBlk_ptr->ICu-number_of_rings_increased) - (SolnBlk_ptr->ICl+number_of_rings_increased) + 1)
                           * ((SolnBlk_ptr->JCu-number_of_rings_increased) - (SolnBlk_ptr->JCl+number_of_rings_increased) + 1)
                           * ((SolnBlk_ptr->KCu-number_of_rings_increased) - (SolnBlk_ptr->KCl+number_of_rings_increased) + 1);
                
    int imin,imax,jmin,jmax,kmin,kmax;
    
    /* For every cell */
    
    int temporary_adaptor_type = adaptor.adaptor_type;
    
    /* ----- Filter ----- */
    //cout << " -- Filter " << endl;
    
    imin = Grid_Blk.ICl+number_of_rings;
    imax = Grid_Blk.ICu-number_of_rings;
    jmin = Grid_Blk.JCl+number_of_rings;
    jmax = Grid_Blk.JCu-number_of_rings;
    kmin = Grid_Blk.KCl+number_of_rings;
    kmax = Grid_Blk.KCu-number_of_rings;
    
    number_of_cells_first = (imax-imin+1)*(jmax-jmin+1)*(kmax-kmin+1);
    
    number_of_processed_cells = 0;
    for(int i=imin; i<=imax; i++) {
        for (int j=jmin; j<=jmax; j++) {
            for (int k=kmin; k<=kmax; k++) {
                Filtered[i][j][k] = filter_ptr->filter(Grid_Blk,Grid_Blk.Cell[i][j][k]);
                //Filtered[i][j][k] = filter_ptr->filter_1D(Grid_Blk,Grid_Blk.Cell[i][j][k],X_DIRECTION);
                number_of_processed_cells++;
                ShowProgress(" -- Filter ",number_of_processed_cells,number_of_cells_first,progress_mode);
            }
        }
    }
    
    /* ----- Divergence ----- */
    //cout << " -- Divergence " << endl;
    imin = Grid_Blk.ICl+number_of_rings_increased;
    imax = Grid_Blk.ICu-number_of_rings_increased;
    jmin = Grid_Blk.JCl+number_of_rings_increased;
    jmax = Grid_Blk.JCu-number_of_rings_increased;
    kmin = Grid_Blk.KCl+number_of_rings_increased;
    kmax = Grid_Blk.KCu-number_of_rings_increased;
    
    number_of_cells_first = (imax-imin+1)*(jmax-jmin+1)*(kmax-kmin+1);

    number_of_processed_cells = 0;
    for(int i=imin; i<=imax; i++) {
        for (int j=jmin; j<=jmax; j++) {
            for (int k=kmin; k<=kmax; k++) {
                //Divergence[i][j][k] = derivative_reconstructor.dfdx(Grid_Blk,Grid_Blk.Cell[i][j][k]);
                Divergence[i][j][k] = finite_differencer.Finite_Difference(Grid_Blk,Grid_Blk.Cell[i][j][k], DFDX);
                number_of_processed_cells++;
                ShowProgress(" -- Divergence ",number_of_processed_cells,number_of_cells_first,progress_mode);
                Truncation_Error[i][j][k] = (RowVector(adaptor.Exact_Derivative(*Input_ptr,
                                                                                Grid_Blk.Cell[i][j][k].Xc.x,
                                                                                Grid_Blk.Cell[i][j][k].Xc.y,
                                                                                Grid_Blk.Cell[i][j][k].Xc.z)
                                             - Divergence[i][j][k])).absolute_values();
            }
        }
    }
    
    /* ----- Filter of Divergence ----- */
    //cout << " -- Filter of Divergence " << endl;
    imin = Grid_Blk.ICl+(number_of_rings+number_of_rings_increased);
    imax = Grid_Blk.ICu-(number_of_rings+number_of_rings_increased);
    jmin = Grid_Blk.JCl+(number_of_rings+number_of_rings_increased);
    jmax = Grid_Blk.JCu-(number_of_rings+number_of_rings_increased);
    kmin = Grid_Blk.KCl+(number_of_rings+number_of_rings_increased);
    kmax = Grid_Blk.KCu-(number_of_rings+number_of_rings_increased);
    
    number_of_cells_second = (imax-imin+1)*(jmax-jmin+1)*(kmax-kmin+1);

    number_of_processed_cells = 0;
    adaptor.Set_Commutation_RowVector(Divergence);
    number_of_processed_cells = 0;
    for(int i=imin; i<=imax; i++) {
        for (int j=jmin; j<=jmax; j++) {
            for (int k=kmin; k<=kmax; k++) {
                Filtered_Divergence[i][j][k] = filter_ptr->filter(Grid_Blk,Grid_Blk.Cell[i][j][k]);
                //Filtered_Divergence[i][j][k] = filter_ptr->filter_1D(Grid_Blk,Grid_Blk.Cell[i][j][k], X_DIRECTION);
                number_of_processed_cells++;
                ShowProgress(" -- Filter of Divergence ",number_of_processed_cells,number_of_cells_second,progress_mode);
                
            }
        }
    }
    
    /* ----- Divergence of Filtered ----- */
    //cout << " -- Divergence of Filtered " << endl;
    number_of_cells_third = (imax-imin+1)*(jmax-jmin+1)*(kmax-kmin+1);

    number_of_processed_cells = 0;
    adaptor.Set_Commutation_RowVector(Filtered);
    number_of_processed_cells = 0;
    for(int i=imin; i<=imax; i++) {
        for (int j=jmin; j<=jmax; j++) {
            for (int k=kmin; k<=kmax; k++) {
                //Divergenced_Filtered[i][j][k] = derivative_reconstructor.dfdx(Grid_Blk,Grid_Blk.Cell[i][j][k]);
                Divergenced_Filtered[i][j][k] = finite_differencer.Finite_Difference(Grid_Blk,Grid_Blk.Cell[i][j][k], DFDX);
                number_of_processed_cells++;
                ShowProgress(" -- Divergence of Filtered ",number_of_processed_cells,number_of_cells_third,progress_mode);
                Commutation_Error[i][j][k] = (RowVector(Filtered_Divergence[i][j][k] - Divergenced_Filtered[i][j][k])).absolute_values();
            }
        }
    }
    
    adaptor.Set_Adaptor(temporary_adaptor_type);
    
    /* for multiblock should be moved to other location!!! */
    RowVector Commutation_Error_maxnorm = maxnorm(Grid_Blk, Commutation_Error);
    RowVector Commutation_Error_L1norm = p_norm(Grid_Blk, Commutation_Error, 1);
    RowVector Commutation_Error_L2norm = p_norm(Grid_Blk, Commutation_Error, 2);
    
    if (!properties.Get_Property_int("batch_flag")){
        cout << "Filter : " << filter_ptr->filter_name() << endl;
        cout << "Commutation Order : " << properties.Get_Property_int("commutation_order") << endl;
        cout << "Grid : " << (Grid_Blk.ICu - Grid_Blk.ICl + 1) 
        << "x" << (Grid_Blk.JCu - Grid_Blk.JCl + 1)
        << "x" << (Grid_Blk.KCu - Grid_Blk.KCl + 1) << endl;
        cout << "Commutation_Error:" << endl;
        cout << "   max norm = " << Commutation_Error_maxnorm;
        cout << "   L1 norm  = " << Commutation_Error_L1norm; 
        cout << "   L2 norm  = " << Commutation_Error_L2norm; 
        cout << endl;
        cout << "Truncation_Error:" << endl;
        cout << "   max norm = " << maxnorm(Grid_Blk, Truncation_Error);
        cout << "   L1 norm  = " << p_norm(Grid_Blk, Truncation_Error, 1); 
        cout << "   L2 norm  = " << p_norm(Grid_Blk, Truncation_Error, 2); 
    }
    
    
    std::stringstream filenamestream;
    filenamestream << "commutation_error_norms_"<< filter_ptr->filter_name() << "_" << properties.Get_Property_int("commutation_order") << ".dat";
    
    ofstream commutation_error_norms_file;
    commutation_error_norms_file.open(filenamestream.str().c_str(),ios::app);    // open file for appending
    assert (!commutation_error_norms_file.fail( )); 
    
    commutation_error_norms_file 
    << Grid_Blk.ICu-Grid_Blk.ICl+1
    << " " << Commutation_Error_maxnorm(0)
    << " " << Commutation_Error_L1norm(0)
    << " " << Commutation_Error_L2norm(0)
    << endl;
    commutation_error_norms_file.close();
    assert (!commutation_error_norms_file.fail( )); 
    
    std::stringstream filenamestream2;
    filenamestream2 << "truncation_error_norms_"<< finite_differencer.Get_order() << ".dat";
    
    ofstream truncation_error_norms_file;
    truncation_error_norms_file.open(filenamestream2.str().c_str(),ios::app);    // open file for appending
    assert (!truncation_error_norms_file.fail( )); 
    
    truncation_error_norms_file 
    << Grid_Blk.ICu-Grid_Blk.ICl+1
    << " " << maxnorm(Grid_Blk, Truncation_Error)(0)
    << " " << p_norm(Grid_Blk, Truncation_Error, 1)(0)
    << " " << p_norm(Grid_Blk, Truncation_Error, 2)(0)
    << endl;
    truncation_error_norms_file.close();
    assert (!truncation_error_norms_file.fail( )); 

    

}


template<typename Soln_pState, typename Soln_cState>
RowVector Explicit_Filters<Soln_pState,Soln_cState>::p_norm(Grid3D_Hexa_Block &Grid_Blk,
                                                      RowVector ***Rows,
                                                      int p) {
    int number_of_rings = properties.Get_Property_int("number_of_rings");
    int number_of_rings_increased = properties.Get_Property_int("number_of_rings_increased");

    int imin = Grid_Blk.ICl+(number_of_rings+number_of_rings_increased),
        imax = Grid_Blk.ICu-(number_of_rings+number_of_rings_increased),
        jmin = Grid_Blk.JCl+(number_of_rings+number_of_rings_increased),
        jmax = Grid_Blk.JCu-(number_of_rings+number_of_rings_increased),
        kmin = Grid_Blk.KCl+(number_of_rings+number_of_rings_increased),
        kmax = Grid_Blk.KCu-(number_of_rings+number_of_rings_increased);
    
    
    int N = Rows[imin][jmin][jmax].size();
    
    RowVector normRow(N);
    normRow.zero();
    
    double total_volume = 0;
    for (int k = kmin ; k <= kmax ; ++k) {
        for (int j = jmin ; j <= jmax ; ++j) {
            for (int i = imin ; i <= imax ; ++i) {
                for (int n = 0; n < N ; ++n) {
                    normRow(n) += pow(Rows[i][j][k](n),double(p))*Grid_Blk.Cell[i][j][k].V;
                }
                total_volume += Grid_Blk.Cell[i][j][k].V;
            }
        }
    }
    
    for (int n = 0; n < N ; ++n) {
        normRow(n) = pow(normRow(n)/total_volume,ONE/double(p));
    }
    return normRow;
}

template<typename Soln_pState, typename Soln_cState>
RowVector Explicit_Filters<Soln_pState,Soln_cState>::maxnorm(Grid3D_Hexa_Block &Grid_Blk,
                                                             RowVector ***Rows) {
    int number_of_rings = properties.Get_Property_int("number_of_rings");
    int number_of_rings_increased = properties.Get_Property_int("number_of_rings_increased");
    
    int 
    imin = Grid_Blk.ICl+(number_of_rings+number_of_rings_increased),
    imax = Grid_Blk.ICu-(number_of_rings+number_of_rings_increased),
    jmin = Grid_Blk.JCl+(number_of_rings+number_of_rings_increased),
    jmax = Grid_Blk.JCu-(number_of_rings+number_of_rings_increased),
    kmin = Grid_Blk.KCl+(number_of_rings+number_of_rings_increased),
    kmax = Grid_Blk.KCu-(number_of_rings+number_of_rings_increased);
    
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
    

    if (Filter_Weights_Allocated) {
        for (int i=0; i<SolnBlk_ptr->NCi; i++) {
            for (int j=0; j<SolnBlk_ptr->NCj; j++) {            
                for (int k=0; k<SolnBlk_ptr->NCk; k++) {
                    Filter_Weights_Assigned[i][j][k] = false;
                }
            }
        }
    }
}

template<typename Soln_pState, typename Soln_cState>
void Explicit_Filters<Soln_pState,Soln_cState>::test(void) {
    if (!initialized) {
        cout << "Explicit_Filter not initialized, can not test" << endl;
        return;
    }
    filter_ptr->filter_tests(SolnBlk_ptr->Grid,
                             
                            SolnBlk_ptr->Grid.Cell[12][12][12]);     
}


template<typename Soln_pState, typename Soln_cState>
int Explicit_Filters<Soln_pState,Soln_cState>::Write_to_file(void){
    if (properties.Get_Property_int("memory_efficient")) {
        cerr <<"\n Memory Efficient mode is ON --> No filter weights to write to file. " << endl;
        return (0);
    }
    if (!properties.Get_Property_int("batch_flag"))
        cout << "\n Writing explicit filter coefficients to file." << endl ;
    int i;
    char prefix[256], cpu_id[256], extension[256], out_file_name[256];
    char *out_file_name_ptr;
    
    strcpy(prefix,properties.Get_Property_string("output_file_name").c_str());
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
    
    for (int i=0; i<SolnBlk_ptr->NCi; i++) {
        for (int j=0; j<SolnBlk_ptr->NCj; j++) {                
            for (int k=0; k<SolnBlk_ptr->NCk; k++) {
                out_file << Filter_Weights_Assigned[i][j][k] << " ";
                if (Filter_Weights_Assigned[i][j][k]) {
                    Filter_Weights[i][j][k].write(out_file);
                } else {
                    out_file << "\n";
                }
            }
        }
    }
    out_file << "\n"; // extra line to separate Grid_Blks
    
    out_file.close();
    return (0);
}




template<typename Soln_pState, typename Soln_cState>
int Explicit_Filters<Soln_pState,Soln_cState>::Read_from_file(void) {
    if (!properties.Get_Property_int("batch_flag"))
        cout << "\n Reading explicit filter coefficients from file." << endl ;
    int i;
    char prefix[256];
    char cpu_id[256], extension[256], in_file_name[256];
    char *in_file_name_ptr;
    
    strcpy(prefix,properties.Get_Property_string("output_file_name").c_str());
    strcat(prefix,  "_explicit_filter");
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
    int filter_type;
    in_file >> filter_type;
    properties.Set_Operating_Property("filter_type",filter_type);
    Create_filter();
    
    bool Store_Filter_Weights = !properties.Get_Property_int("memory_efficient");

    if (Store_Filter_Weights) {
        Allocate_Filter_Weights();
        for (int i=0; i<SolnBlk_ptr->NCi; i++) {
            for (int j=0; j<SolnBlk_ptr->NCj; j++) {                
                for (int k=0; k<SolnBlk_ptr->NCk; k++) {
                    in_file.setf(ios::skipws);
                    in_file >> Filter_Weights_Assigned[i][j][k];
                    in_file.unsetf(ios::skipws);
                    if (Filter_Weights_Assigned[i][j][k]) {
                        Filter_Weights[i][j][k].read(in_file);
                    }
                }
            }
        }
        Filter_Weights_Allocated = true;
        
    }
    
    in_file.unsetf(ios::skipws);
    in_file.close();
    properties.Set_Operating_Property("restarted",true);
    return (0);
}


template<typename Soln_pState, typename Soln_cState>
void Explicit_Filters<Soln_pState,Soln_cState>::Output_Commutation(Grid3D_Hexa_Block &Grid_Blk,
                                                                   RowVector ***filtered,
                                                                   RowVector ***divergence,
                                                                   RowVector ***filtered_divergence,
                                                                   RowVector ***divergence_filtered,
                                                                   RowVector ***commutation_error,
                                                                   RowVector ***truncation_error,
                                                                   int Block_Number,
                                                                   bool Output_Title,
                                                                   ofstream &Out_File){
    
    
    
    int nVar = 0;
    
    int number_of_rings = properties.Get_Property_int("number_of_rings");
    int number_of_rings_increased = properties.Get_Property_int("number_of_rings_increased");

    
    int 
    imin = Grid_Blk.ICl+(number_of_rings+number_of_rings_increased),
    imax = Grid_Blk.ICu-(number_of_rings+number_of_rings_increased),
    jmin = Grid_Blk.JCl+(number_of_rings+number_of_rings_increased),
    jmax = Grid_Blk.JCu-(number_of_rings+number_of_rings_increased),
    kmin = Grid_Blk.KCl+(number_of_rings+number_of_rings_increased),
    kmax = Grid_Blk.KCu-(number_of_rings+number_of_rings_increased);
    
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
        << "\"exact divergence\" \\ \n"
        << "\"filtered of divergence\" \\ \n"
        << "\"divergence of filtered\" \\ \n"
        << "\"commutation error\" \\ \n"
        << "\"truncation error\" \\ \n";
        Out_File << "DATASETAUXDATA filter = \"" << filter_ptr->filter_name() << "\" \\ \n"
                 << "DATASETAUXDATA commutation_order = \"" << properties.Get_Property_int("commutation_order") << "\" \\ \n"
                 << "DATASETAUXDATA filter_grid_ratio = \"" << properties.Get_Property_double("FGR") << "\" \\ \n"
                 << "DATASETAUXDATA number_of_rings = \"" << properties.Get_Property_int("number_of_rings")  << "\" \\ \n";
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
                            << " " << adaptor.Exact_Derivative(*Input_ptr,
                                                               Grid_Blk.Cell[i][j][k].Xc.x,
                                                               Grid_Blk.Cell[i][j][k].Xc.y,
                                                               Grid_Blk.Cell[i][j][k].Xc.z)(nVar)
                            << " " << filtered_divergence[i][j][k](nVar)
                            << " " << divergence_filtered[i][j][k](nVar)
                            << " " << commutation_error[i][j][k](nVar)
                            << " " << truncation_error[i][j][k](nVar) << "\n";
                Out_File.unsetf(ios::scientific);
            } /* endfor */
        } /* endfor */
    } /* endfor */
    
    Out_File << setprecision(6);
    
}

template<typename Soln_pState, typename Soln_cState>
void Explicit_Filters<Soln_pState,Soln_cState>::Allocate_Filter_Weights(void) {
    Deallocate_Filter_Weights();
    Filter_Weights = new RowVector **[SolnBlk_ptr->NCi];
    Filter_Weights_Assigned = new bool **[SolnBlk_ptr->NCi];
    for (int i=0; i<SolnBlk_ptr->NCi; i++) {
        Filter_Weights[i] = new RowVector *[SolnBlk_ptr->NCj];
        Filter_Weights_Assigned[i] = new bool *[SolnBlk_ptr->NCj];
        for (int j=0; j<SolnBlk_ptr->NCj; j++) {
            Filter_Weights[i][j] = new RowVector [SolnBlk_ptr->NCk];
            Filter_Weights_Assigned[i][j] = new bool [SolnBlk_ptr->NCk];
            
            for (int k=0; k<SolnBlk_ptr->NCk; k++) {
                Filter_Weights_Assigned[i][j][k] = false;
            }
        }
    }
    Filter_Weights_Allocated = true;
}

template<typename Soln_pState, typename Soln_cState>
void Explicit_Filters<Soln_pState,Soln_cState>::Deallocate_Filter_Weights(void) {
    if (Filter_Weights_Allocated) {
        for (int i=0; i<SolnBlk_ptr->NCi; i++) {
            for (int j=0; j<SolnBlk_ptr->NCj; j++) {
                delete[] Filter_Weights[i][j];          Filter_Weights[i][j] = NULL;
                delete[] Filter_Weights_Assigned[i][j]; Filter_Weights_Assigned[i][j] = NULL;
            }
            delete[] Filter_Weights[i];             Filter_Weights[i] = NULL;
            delete[] Filter_Weights_Assigned[i];    Filter_Weights_Assigned[i] = NULL;
        }
        delete[] Filter_Weights;            Filter_Weights = NULL;
        delete[] Filter_Weights_Assigned;   Filter_Weights_Assigned = NULL;
    }
    Filter_Weights_Allocated = false;
}

template<typename Soln_pState, typename Soln_cState>
void Explicit_Filters<Soln_pState,Soln_cState>::Allocate_Derivative_Reconstruction_Weights(void) {
    Deallocate_Derivative_Reconstruction_Weights();
    Derivative_Reconstruction_Weights = new DenseMatrix **[SolnBlk_ptr->NCi];
    Derivative_Reconstruction_Weights_Assigned = new bool **[SolnBlk_ptr->NCi];
    for (int i=0; i<SolnBlk_ptr->NCi; i++) {
        Derivative_Reconstruction_Weights[i] = new DenseMatrix *[SolnBlk_ptr->NCj];
        Derivative_Reconstruction_Weights_Assigned[i] = new bool *[SolnBlk_ptr->NCj];
        for (int j=0; j<SolnBlk_ptr->NCj; j++) {
            Derivative_Reconstruction_Weights[i][j] = new DenseMatrix [SolnBlk_ptr->NCk];
            Derivative_Reconstruction_Weights_Assigned[i][j] = new bool [SolnBlk_ptr->NCk];
            
            for (int k=0; k<SolnBlk_ptr->NCk; k++) {
                Derivative_Reconstruction_Weights_Assigned[i][j][k] = false;
            }
        }
    }
    Derivative_Reconstruction_Weights_Allocated = true;
}

template<typename Soln_pState, typename Soln_cState>
void Explicit_Filters<Soln_pState,Soln_cState>::Deallocate_Derivative_Reconstruction_Weights(void) {
    if (Derivative_Reconstruction_Weights_Allocated) {
        for (int i=0; i<SolnBlk_ptr->NCi; i++) {
            for (int j=0; j<SolnBlk_ptr->NCj; j++) {
                delete[] Derivative_Reconstruction_Weights[i][j];          Derivative_Reconstruction_Weights[i][j] = NULL;
                delete[] Derivative_Reconstruction_Weights_Assigned[i][j]; Derivative_Reconstruction_Weights_Assigned[i][j] = NULL;
            }
            delete[] Derivative_Reconstruction_Weights[i];             Derivative_Reconstruction_Weights[i] = NULL;
            delete[] Derivative_Reconstruction_Weights_Assigned[i];    Derivative_Reconstruction_Weights_Assigned[i] = NULL;
        }
        delete[] Derivative_Reconstruction_Weights;            Derivative_Reconstruction_Weights = NULL;
        delete[] Derivative_Reconstruction_Weights_Assigned;   Derivative_Reconstruction_Weights_Assigned = NULL;
    }
    Derivative_Reconstruction_Weights_Allocated = false;
}



template <typename Soln_pState, typename Soln_cState>
void Explicit_Filters<Soln_pState,Soln_cState>::ShowProgress(std::string message, int numIn, int maximum, int mode) {
    if (properties.Get_Property_int("batch_flag"))
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

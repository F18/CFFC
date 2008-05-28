/*
 *  LES_Filters.h
 *  CFFC
 *
 *  Created by Willem Deconinck on 25/03/08.
 *
 */

#ifndef _LES_FILTERS_INCLUDED
#define _LES_FILTERS_INCLUDED


#include "../../HexaBlock/HexaMultiBlock.h"
#include "../../HexaBlock/HexaSolverClasses.h"
#include "General_Filter.h"
#include "Haselbacher_Filter.h"
#include "Vasilyev_Filter.h"
#include "Vasilyev_LS_Filter.h"


/* ------------------------------------------------------------------------------------------------------------------------------ */
/**
 * CLASS : LES_Filter
 */

#define SOLN_CSTATE_3D      0
#define SOLN_CSTATE_4D      1
#define SOLN_PSTATE_DOUBLE  2

template<typename Soln_pState, typename Soln_cState>
class LES_Filter {
public:
    typedef double (Soln_pState::*Soln_pState_member_ptr_type);
    typedef Soln_cState *** (Hexa_Block<Soln_pState,Soln_cState>::*Soln_cState_3D_ptr_type);
    typedef Soln_cState **** (Hexa_Block<Soln_pState,Soln_cState>::*Soln_cState_4D_ptr_type);

    static Soln_pState_member_ptr_type Soln_pState_member_ptr;
    static Soln_cState_3D_ptr_type Soln_cState_3D_ptr;
    static Soln_cState_4D_ptr_type Soln_cState_4D_ptr;

    static int dUdt_k_residual;
    static int filter_variable_type;
    
    static int commutation_order;
    static double FGR;
    static int number_of_rings;
    
    static bool restarted;
    
    int filter_type;
    
    const char *output_file_name;
    
    General_Filter<Soln_pState,Soln_cState> *filter_ptr;
    
    AdaptiveBlock3D_List *LocalSolnBlkList_ptr;                     // List with properties of SolnBlks
    Hexa_Block<Soln_pState,Soln_cState> *Solution_Blocks_ptr;       // array of SolnBlks
    
    /* ----- constructor ----- */
    LES_Filter(HexaSolver_Data &Data,
               HexaSolver_Solution_Data<Soln_pState,Soln_cState> &Solution_Data,
               int filter_flag) {
        
        Solution_Blocks_ptr  = Solution_Data.Local_Solution_Blocks.Soln_Blks;
        LocalSolnBlkList_ptr = &(Data.Local_Adaptive_Block_List);
        FILTER_ONLY_ONE_SOLNBLK = false;
        FGR = Solution_Data.Input.Turbulence_IP.FGR;
        commutation_order = Solution_Data.Input.Turbulence_IP.commutation_order;
        number_of_rings = Solution_Data.Input.Turbulence_IP.number_of_rings;
        output_file_name = Solution_Data.Input.Output_File_Name;
        filter_type = filter_flag;
        
        Create_filter();
    }
    
    LES_Filter(HexaSolver_Data &Data,
               HexaSolver_Solution_Data<Soln_pState,Soln_cState> &Solution_Data) {
        
        Solution_Blocks_ptr  = Solution_Data.Local_Solution_Blocks.Soln_Blks;
        LocalSolnBlkList_ptr = &(Data.Local_Adaptive_Block_List);
        FILTER_ONLY_ONE_SOLNBLK = false;
        FGR = Solution_Data.Input.Turbulence_IP.FGR;
        commutation_order = Solution_Data.Input.Turbulence_IP.commutation_order;
        number_of_rings = Solution_Data.Input.Turbulence_IP.number_of_rings;
        filter_type = Solution_Data.Input.Turbulence_IP.i_filter_type;
        output_file_name = Solution_Data.Input.Output_File_Name;

        if (Solution_Data.Input.i_ICs == IC_RESTART && !restarted)
            filter_type = FILTER_TYPE_RESTART;
        
        Create_filter();
        Solution_Data.Input.Turbulence_IP.i_filter_type = filter_type;
    }
    
    LES_Filter(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk,
               Input_Parameters<Soln_pState,Soln_cState> &IPs,
               int filter_flag) {
        FILTER_ONLY_ONE_SOLNBLK = true;
        Solution_Blocks_ptr  = &SolnBlk;
        FGR = IPs.Turbulence_IP.FGR;
        commutation_order = IPs.Turbulence_IP.commutation_order;
        number_of_rings = IPs.Turbulence_IP.number_of_rings;
        filter_type = filter_flag;
        
        Create_filter();
    }
    
    
    LES_Filter(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk,
               Input_Parameters<Soln_pState,Soln_cState> &IPs) {
        FILTER_ONLY_ONE_SOLNBLK = true;
        Solution_Blocks_ptr  = &SolnBlk;
        FGR = IPs.Turbulence_IP.FGR;
        commutation_order = IPs.Turbulence_IP.commutation_order;
        number_of_rings = IPs.Turbulence_IP.number_of_rings;
        filter_type = IPs.Turbulence_IP.i_filter_type;
        if (filter_type == FILTER_TYPE_RESTART) {
            cerr << "Cannot read explicit filter from file with this constructor";
        }

        Create_filter();
    }
    
    ~LES_Filter() {
        delete filter_ptr;
    }
    
    
    void Create_filter(void) {
        switch (filter_type) {
            case FILTER_TYPE_HASELBACHER:
                filter_ptr = new Haselbacher_Filter<Soln_pState,Soln_cState>;
                break;
            case FILTER_TYPE_VASILYEV:
                filter_ptr = new Vasilyev_Filter<Soln_pState,Soln_cState>;
                break;
            case FILTER_TYPE_VASILYEV_LS:
                filter_ptr = new Vasilyev_LS_Filter<Soln_pState,Soln_cState>;
                break;
            case FILTER_TYPE_RESTART:
                int error_flag = Read_from_file();
                if (error_flag == 1) cerr << "could not read filter_input_file" << endl;
                break;
            default:
                cerr << "Filter not defined" << endl;
                break;
        }
    }
    
    int Read_from_file(void) {
        cout << "\n Reading explicit filter coefficients from file." << endl ;
        int i;
        char prefix[256], cpu_id[256], extension[256], in_file_name[256];
        char *in_file_name_ptr;
        
        i = 0;
        while (1) {
            if (output_file_name[i] == ' ' ||
                output_file_name[i] == '.') break;
            prefix[i] = output_file_name[i];
            i = i + 1;
            if (i > strlen(output_file_name) ) break;
        } /* endwhile */
        prefix[i] = '\0';
        strcat(prefix, "_explicit_filter");
        strcpy(extension, ".dat");
        
        sprintf(cpu_id,"_cpu_%d",LocalSolnBlkList_ptr->ThisCPU);
        
        strcpy(in_file_name, prefix);
        strcat(in_file_name, cpu_id);
        strcat(in_file_name, extension);
        in_file_name_ptr = in_file_name;
        
        
        ifstream in_file;
        in_file.open(in_file_name, ios::in);
        if (in_file.bad()) return (1);        
        in_file.setf(ios::skipws);
                
        in_file >> filter_type;
        Create_filter();
    
        if (LocalSolnBlkList_ptr->Nused() >= 1) {
            for (int nBlk = 0; nBlk < LocalSolnBlkList_ptr->Nused(); nBlk++ ) {
                if (LocalSolnBlkList_ptr->Block[nBlk].used == ADAPTIVEBLOCK3D_USED) {
                    filter_ptr->Read_from_file(Solution_Blocks_ptr[nBlk],in_file);
                }         
            } 
        } 
        
        
        
        in_file.unsetf(ios::skipws);
        in_file.close();
        restarted = true;
        return (0);
    }
    
    int Write_to_file(void){
        
        int i;
        char prefix[256], cpu_id[256], extension[256], out_file_name[256];
        char *out_file_name_ptr;
                
        i = 0;
        while (1) {
            if (output_file_name[i] == ' ' ||
                output_file_name[i] == '.') break;
            prefix[i] = output_file_name[i];
            i = i + 1;
            if (i > strlen(output_file_name) ) break;
        } /* endwhile */
        prefix[i] = '\0';
        strcat(prefix, "_explicit_filter");
        strcpy(extension, ".dat");

        sprintf(cpu_id,"_cpu_%d",LocalSolnBlkList_ptr->ThisCPU);
        
        strcpy(out_file_name, prefix);
        strcat(out_file_name, cpu_id);
        strcat(out_file_name, extension);
        out_file_name_ptr = out_file_name;
        
        ofstream out_file;
        out_file.open(out_file_name, ios::out);
        if (out_file.bad()) return (1);
        
        out_file << filter_ptr->filter_type() << "\n";                    

        if (LocalSolnBlkList_ptr->Nused() >= 1) {
            for (int nBlk = 0; nBlk < LocalSolnBlkList_ptr->Nused(); nBlk++ ) {
                if (LocalSolnBlkList_ptr->Block[nBlk].used == ADAPTIVEBLOCK3D_USED) {
                    filter_ptr->Write_to_file(Solution_Blocks_ptr[nBlk],out_file);
                }         
            } 
        } 
        out_file.close();
        return (0);
    }
    
    void filter(double Soln_pState::*&member) {
        Soln_pState_member_ptr = member;
        filter_variable_type = SOLN_PSTATE_DOUBLE;
        filter_Blocks();
    }
    void filter(Soln_cState *** Hexa_Block<Soln_pState,Soln_cState>::*&member) {
        Soln_cState_3D_ptr = member;
        filter_variable_type = SOLN_CSTATE_3D;
        filter_Blocks();
    }
    void filter(Soln_cState **** Hexa_Block<Soln_pState,Soln_cState>::*&member,int k_residual) {
        Soln_cState_4D_ptr = member;
        dUdt_k_residual = k_residual;
        filter_variable_type = SOLN_CSTATE_4D;
        filter_Blocks();
    }
    
    static void what_to_filter(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, int i, int j, int k, RowVector &x);
    static void what_to_filter(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, int i, int j, int k, DenseMatrix &b, int row_index);

    void filter_Blocks(void);
    
    double filter_width;
    void transfer_function();
    
    void reset(void);
    void test(void);
    
    double maximum_wavenumber();
    bool FILTER_ONLY_ONE_SOLNBLK;
    RowVector ***Filtered;
    
    void allocate_Filtered(Hexa_Block<Soln_pState,Soln_cState> &Soln_Blk){
        Filtered = new RowVector **[Soln_Blk.NCi];
        for (int i=0; i<Soln_Blk.NCi; i++) {
            Filtered[i] = new RowVector *[Soln_Blk.NCj];
            for (int j=0; j<Soln_Blk.NCj; j++) {
                Filtered[i][j] = new RowVector [Soln_Blk.NCk];
            }
        }
    }
    
    void deallocate_Filtered(Hexa_Block<Soln_pState,Soln_cState> &Soln_Blk) {
        for (int i=0; i<Soln_Blk.NCi; i++) {
            for (int j=0; j<Soln_Blk.NCj; j++) {
                delete[] Filtered[i][j];   Filtered[i][j] = NULL;
            }
            delete[] Filtered[i];   Filtered[i] = NULL;
        }
        delete[] Filtered;   Filtered = NULL;
    }
    
};

template<typename Soln_pState, typename Soln_cState>
int LES_Filter<Soln_pState,Soln_cState>::commutation_order = 2;

template<typename Soln_pState, typename Soln_cState>
int LES_Filter<Soln_pState,Soln_cState>::number_of_rings = 2;

template<typename Soln_pState, typename Soln_cState>
double LES_Filter<Soln_pState,Soln_cState>::FGR = 2.0;

template<typename Soln_pState, typename Soln_cState>
int LES_Filter<Soln_pState,Soln_cState>::filter_variable_type = SOLN_CSTATE_4D;

template<typename Soln_pState, typename Soln_cState>
int LES_Filter<Soln_pState,Soln_cState>::dUdt_k_residual = 0;

template<typename Soln_pState, typename Soln_cState>
double Soln_pState::* LES_Filter<Soln_pState,Soln_cState>::Soln_pState_member_ptr = NULL;

template<typename Soln_pState, typename Soln_cState>
Soln_cState *** Hexa_Block<Soln_pState,Soln_cState>::* LES_Filter<Soln_pState,Soln_cState>::Soln_cState_3D_ptr = NULL;

template<typename Soln_pState, typename Soln_cState>
Soln_cState **** Hexa_Block<Soln_pState,Soln_cState>::* LES_Filter<Soln_pState,Soln_cState>::Soln_cState_4D_ptr = NULL;

template<typename Soln_pState, typename Soln_cState>
bool LES_Filter<Soln_pState,Soln_cState>::restarted = false;

template<typename Soln_pState, typename Soln_cState>
 void LES_Filter<Soln_pState,Soln_cState>::what_to_filter(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, int i, int j, int k , RowVector &x) {
     int numvar;
     switch (filter_variable_type) {
        case SOLN_PSTATE_DOUBLE:
            x = RowVector(1);
            x(0) = SolnBlk.W[i][j][k].*Soln_pState_member_ptr;
            break;
        case SOLN_CSTATE_3D:
            numvar = SolnBlk.NumVar();
            x = RowVector(numvar);
            for (int n=1; n<=numvar; n++) {
                x(n-1) = (SolnBlk.*Soln_cState_3D_ptr)[i][j][k][n];
            }
            break;
        case SOLN_CSTATE_4D:
            numvar = SolnBlk.NumVar();
            x = RowVector(numvar);
            for (int n=1; n<=numvar; n++) {
                x(n-1) = (SolnBlk.*Soln_cState_4D_ptr)[i][j][k][dUdt_k_residual][n];
            }
            break;
    }
}


template<typename Soln_pState, typename Soln_cState>
void LES_Filter<Soln_pState,Soln_cState>::what_to_filter(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, int i, int j, int k , DenseMatrix &b, int row_index) {
    int numvar;
    switch (filter_variable_type) {
        case SOLN_PSTATE_DOUBLE:
            b(row_index,0) = SolnBlk.W[i][j][k].*Soln_pState_member_ptr;
            break;
        case SOLN_CSTATE_3D:
            numvar = SolnBlk.NumVar();
            for (int n=1; n<=numvar; n++) {
                b(row_index,n-1) = (SolnBlk.*Soln_cState_3D_ptr)[i][j][k][n];
            }
            break;
            case SOLN_CSTATE_4D:
            numvar = SolnBlk.NumVar();
            for (int n=1; n<=numvar; n++) {
                b(row_index,n-1) = (SolnBlk.*Soln_cState_4D_ptr)[i][j][k][dUdt_k_residual][n];
            }
            break;
    }
}


template<typename Soln_pState, typename Soln_cState>
double LES_Filter<Soln_pState,Soln_cState>::maximum_wavenumber() {
    double max_cell_volume = Max_Cell_Volume(Solution_Blocks_ptr,*LocalSolnBlkList_ptr);
    double kmax = PI/pow(max_cell_volume,ONE/THREE);
    return kmax;
}


template<typename Soln_pState, typename Soln_cState>
void LES_Filter<Soln_pState,Soln_cState>::transfer_function() {
    if (FILTER_ONLY_ONE_SOLNBLK) {
        filter_ptr->transfer_function(*Solution_Blocks_ptr,Solution_Blocks_ptr->Grid.Cell[12][12][12]);            
    }
    else if (LocalSolnBlkList_ptr->Nused() >= 1) {
        for (int nBlk = 0; nBlk <= LocalSolnBlkList_ptr->Nused(); ++nBlk ) {
            if (LocalSolnBlkList_ptr->Block[nBlk].used == ADAPTIVEBLOCK3D_USED) {
                filter_ptr->transfer_function(Solution_Blocks_ptr[nBlk],Solution_Blocks_ptr[nBlk].Grid.Cell[12][12][12]);            
                return; // This makes sure it is called only once
            }
        }
    }
}

template<typename Soln_pState, typename Soln_cState>
void LES_Filter<Soln_pState,Soln_cState>::reset() {
    if (LocalSolnBlkList_ptr->Nused() >= 1) {
        for (int nBlk = 0; nBlk <= LocalSolnBlkList_ptr->Nused(); ++nBlk ) {
            if (LocalSolnBlkList_ptr->Block[nBlk].used == ADAPTIVEBLOCK3D_USED) {
                filter_ptr->Reset_Filter_Weights(Solution_Blocks_ptr[nBlk]);
                return;
            }
        }
    }
}

template<typename Soln_pState, typename Soln_cState>
void LES_Filter<Soln_pState,Soln_cState>::test() {
    double kmax = maximum_wavenumber();
    if (LocalSolnBlkList_ptr->Nused() >= 1) {
        for (int nBlk = 0; nBlk <= LocalSolnBlkList_ptr->Nused(); ++nBlk ) {
            if (LocalSolnBlkList_ptr->Block[nBlk].used == ADAPTIVEBLOCK3D_USED) {                
                filter_ptr->filter_tests(Solution_Blocks_ptr[nBlk],Solution_Blocks_ptr[nBlk].Grid.Cell[12][12][12]);
                return;
            }
        }
    }
}


template<typename Soln_pState, typename Soln_cState>
void LES_Filter<Soln_pState,Soln_cState>::filter_Blocks(void) {
    if (FILTER_ONLY_ONE_SOLNBLK) {
        allocate_Filtered(*Solution_Blocks_ptr);
        for(int i=Solution_Blocks_ptr->ICl; i<=Solution_Blocks_ptr->ICu; i++) {
            for (int j=Solution_Blocks_ptr->JCl; j<=Solution_Blocks_ptr->JCu; j++) {
                for (int k=Solution_Blocks_ptr->KCl; k<=Solution_Blocks_ptr->KCu; k++) {
                    Filtered[i][j][k] = filter_ptr->filter(*Solution_Blocks_ptr,Solution_Blocks_ptr->Grid.Cell[i][j][k]);
                }
            }
        }
        for(int i=Solution_Blocks_ptr->ICl; i<=Solution_Blocks_ptr->ICu; i++) {
            for (int j=Solution_Blocks_ptr->JCl; j<=Solution_Blocks_ptr->JCu; j++) {
                for (int k=Solution_Blocks_ptr->KCl; k<=Solution_Blocks_ptr->KCu; k++) {
                    switch (filter_variable_type) {
                        case SOLN_PSTATE_DOUBLE:
                            Solution_Blocks_ptr->W[i][j][k].*Soln_pState_member_ptr = Filtered[i][j][k](0);
                            break;
                        case SOLN_CSTATE_3D:
                            for (int n=1; n<=Solution_Blocks_ptr->NumVar(); n++)
                                (Solution_Blocks_ptr->*Soln_cState_3D_ptr)[i][j][k][n] = Filtered[i][j][k](n-1);
                            break;
                        case SOLN_CSTATE_4D:
                            for (int n=1; n<=Solution_Blocks_ptr->NumVar(); n++)
                                (Solution_Blocks_ptr->*Soln_cState_4D_ptr)[i][j][k][dUdt_k_residual][n] = Filtered[i][j][k](n-1);
                            break;
                    }
                }
            }
        }
        deallocate_Filtered(*Solution_Blocks_ptr);        
    }
    else {
        /* For every local solution block */
        if (LocalSolnBlkList_ptr->Nused() >= 1) {
            for (int nBlk = 0; nBlk < LocalSolnBlkList_ptr->Nused(); nBlk++ ) {
                if (LocalSolnBlkList_ptr->Block[nBlk].used == ADAPTIVEBLOCK3D_USED) {
                    /* For every cell */
                    allocate_Filtered(Solution_Blocks_ptr[nBlk]);
                    for(int i=Solution_Blocks_ptr[nBlk].ICl; i<=Solution_Blocks_ptr[nBlk].ICu; i++) {
                        for (int j=Solution_Blocks_ptr[nBlk].JCl; j<=Solution_Blocks_ptr[nBlk].JCu; j++) {
                            for (int k=Solution_Blocks_ptr[nBlk].KCl; k<=Solution_Blocks_ptr[nBlk].KCu; k++) {
                                Filtered[i][j][k] = filter_ptr->filter(Solution_Blocks_ptr[nBlk],Solution_Blocks_ptr[nBlk].Grid.Cell[i][j][k]);
                            }
                        }
                    }
                    for(int i=Solution_Blocks_ptr[nBlk].ICl; i<=Solution_Blocks_ptr[nBlk].ICu; i++) {
                        for (int j=Solution_Blocks_ptr[nBlk].JCl; j<=Solution_Blocks_ptr[nBlk].JCu; j++) {
                            for (int k=Solution_Blocks_ptr[nBlk].KCl; k<=Solution_Blocks_ptr[nBlk].KCu; k++) {
                                switch (filter_variable_type) {
                                    case SOLN_PSTATE_DOUBLE:
                                        Solution_Blocks_ptr[nBlk].W[i][j][k].*Soln_pState_member_ptr = Filtered[i][j][k](0);
                                        break;
                                    case SOLN_CSTATE_3D:
                                        for (int n=1; n<=Solution_Blocks_ptr->NumVar(); n++)
                                            (Solution_Blocks_ptr[nBlk].*Soln_cState_3D_ptr)[i][j][k][n] = Filtered[i][j][k](n-1);
                                        break;
                                    case SOLN_CSTATE_4D:
                                        for (int n=1; n<=Solution_Blocks_ptr[nBlk].NumVar(); n++)
                                            (Solution_Blocks_ptr[nBlk].*Soln_cState_4D_ptr)[i][j][k][dUdt_k_residual][n] = Filtered[i][j][k](n-1);
                                        break;
                                }
                            }
                        }
                    }
                    deallocate_Filtered(Solution_Blocks_ptr[nBlk]);
                }         
            } 
        }         
    }
}








/* ----------------------------------------------------------------------------------------------------------------------- 
typedef double (theClass::*function_with_one_argument) (const double &abs_wave_num) const;
_Member_Function_Wrapper_<theClass,function_with_one_argument, double> mapped_function(this, &theClass::Energy_Spectrum_Value);
double dummy;
double TKE_entire_range = AdaptiveGaussianQuadrature(mapped_function, 0.0, k_eta/2.0, dummy,5);
/**/



template<typename Soln_pState, typename Soln_cState>
inline double return_it(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, double (func)(Soln_pState &)){
    /* For static member functions or any function */
    double b;
    b = func(SolnBlk.W[2][2][2]);
    
    return b;
}

template<typename Soln_pState, typename Soln_cState>
inline double return_it_2(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, double (Soln_pState::*&func)(void)){
    /* for non static member functions */
    double b;
    b = (SolnBlk.W[2][2][2].*func)();
    
    return b;
}

template<typename Soln_pState, typename Soln_cState>
inline double return_it_2(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, double (Soln_pState::*&func)){
    /* for non static member functions */
    double b;
    b = (SolnBlk.W[2][2][2].*func);
    
    return b;
}



template<typename Soln_pState, typename Soln_cState, typename func_object>
inline double return_it_3(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, func_object func) {
    /* for non static member functions */
    double b;
    b = (SolnBlk.W[2][2][2].*func);
    
    return b;
}

template<typename Soln_pState, typename Soln_cState, typename func_object>
inline double return_it_4(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, func_object func) {
    /* for non static member functions */
    double b;
    cout << "func = " << (SolnBlk.*func)[2][2][2] << endl;
    b = (SolnBlk.*func)[2][2][2].rho;
    
    return b;
}


template<typename Soln_pState, typename Soln_cState>
inline double return_it_44(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Soln_pState *** Hexa_Block<Soln_pState,Soln_cState>::*&func) {
    /* for non static member functions */
    double b;
    cout << "44" << endl;
    cout << "func = " << (SolnBlk.*func)[2][2][2] << endl;
    b = (SolnBlk.*func)[2][2][2].rho;
    
    return b;
}


template<typename Soln_pState, typename Soln_cState, typename func_object>
inline double return_it_5(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, func_object func, int istage) {
    /* for non static member functions */
    double b;
    cout << "func = " << (SolnBlk.*func)[2][2][2][istage] << endl;
    b = (SolnBlk.*func)[2][2][2][istage].rho;
    
    return b;
}



#endif

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

#define LES_FILTER_HASELBACHER  1
#define LES_FILTER_VASILYEV     2



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
    
    AdaptiveBlock3D_List *LocalSolnBlkList_ptr;                     // List with properties of SolnBlks
    Hexa_Block<Soln_pState,Soln_cState> *Solution_Blocks_ptr;       // array of SolnBlks
    /* ----- constructor ----- */
    LES_Filter(HexaSolver_Data &Data,
               HexaSolver_Solution_Data<Soln_pState,Soln_cState> &Solution_Data,
               int filter_flag) {
        
        Solution_Blocks_ptr  = Solution_Data.Local_Solution_Blocks.Soln_Blks;
        LocalSolnBlkList_ptr = &(Data.Local_Adaptive_Block_List);
        FILTER_ONLY_ONE_SOLNBLK = false;

        switch (filter_flag) {
            case LES_FILTER_HASELBACHER:
                filter_ptr = new Haselbacher_Filter<Soln_pState,Soln_cState>;
                break;
            case LES_FILTER_VASILYEV:
                filter_ptr = new Vasilyev_Filter<Soln_pState,Soln_cState>;
                break;
        }
    }
    
    LES_Filter(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk,
               int filter_flag) {
        FILTER_ONLY_ONE_SOLNBLK = true;
        Solution_Blocks_ptr  = &SolnBlk;
        
        switch (filter_flag) {
            case LES_FILTER_HASELBACHER:
                filter_ptr = new Haselbacher_Filter<Soln_pState,Soln_cState>;
                break;
            case LES_FILTER_VASILYEV:
                filter_ptr = new Vasilyev_Filter<Soln_pState,Soln_cState>;
                break;
        }
    }

    ~LES_Filter() {
        delete filter_ptr;
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

    
    General_Filter<Soln_pState,Soln_cState> *filter_ptr;
    
    void filter_Blocks(void);
    
    double filter_width;
    void transfer_function();
    
    double gaussian();
    double tophat();

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
inline void LES_Filter<Soln_pState,Soln_cState>::what_to_filter(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, int i, int j, int k , RowVector &x) {
    switch (filter_variable_type) {
        case SOLN_PSTATE_DOUBLE:
            x = RowVector(1);
            x(0) = SolnBlk.W[i][j][k].*Soln_pState_member_ptr;
            break;
        case SOLN_CSTATE_3D:
            x = RowVector(SolnBlk.NumVar());
            for (int n=1; n<=SolnBlk.NumVar(); n++) {
                x(n-1) = (SolnBlk.*Soln_cState_3D_ptr)[i][j][k][n];
            }
            break;
        case SOLN_CSTATE_4D:
            x = RowVector(SolnBlk.NumVar());
            for (int n=1; n<=SolnBlk.NumVar(); n++) {
                x(n-1) = (SolnBlk.*Soln_cState_4D_ptr)[i][j][k][dUdt_k_residual][n];
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
    double kmax = maximum_wavenumber();
    if (LocalSolnBlkList_ptr->Nused() >= 1) {
        for (int nBlk = 0; nBlk <= LocalSolnBlkList_ptr->Nused(); ++nBlk ) {
            if (LocalSolnBlkList_ptr->Block[nBlk].used == ADAPTIVEBLOCK3D_USED) {
            
                filter_ptr->transfer_function(Solution_Blocks_ptr[nBlk],Solution_Blocks_ptr[nBlk].Grid.Cell[6][6][6],kmax);
            
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

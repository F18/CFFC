/*
 *  Explicit_Filter_Properties.h
 *  CFFC
 *
 *  Created by Willem Deconinck on 25/03/08.
 *
 */

#ifndef _EXPLICIT_FILTER_HELPERS_INCLUDED
#define _EXPLICIT_FILTER_HELPERS_INCLUDED

#ifndef _NEIGHBOURS_INCLUDED
#include "Neighbours.h"
#endif

#ifndef _FILTER_STATE_INCLUDED
#include "Filter_State.h"
#endif

#ifndef _UTILITIES_INCLUDED
#include "../Utilities/Utilities.h"
#endif


#include <map>


template <typename Soln_pState, typename Soln_cState> 
class Hexa_Block;

/* ------------------------------------------------------------------------------------------------------------------------------ */
/**
 * CLASS : Explicit_Filter_Properties
 */

class Explicit_Filter_Properties {
private:
    std::map<std::string,std::string> property_map;
    bool properties_changed;

    template <typename T>
    void Set_Property(const std::string property_name, const T& property_value) {
        property_map[property_name] = to_str(property_value);
    }
    
public:
    template <typename T>
    T Get_Property(const std::string property_name) {
        return from_str<T>(property_map[property_name]);
    }
    
    template <typename T>
    void Get_Property(const std::string property_name, T& property_value) {
        property_value = from_str<T>(property_map[property_name]);
    }
    
    void Get_Property(const std::string property_name, std::string& property_value) {
        property_value = property_map[property_name];
    }
    
    
public:
    bool Changed(void) {
        return properties_changed;
    }
    
    void Properties_Read(void) {
        properties_changed = false;
    }
    
    
    template <typename Soln_pState, typename Soln_cState>
    void Set_Properties(Explicit_Filter_Constants::Filter_Number filter_number, Input_Parameters<Soln_pState,Soln_cState> &IPs, int batch_flag){
        Set_Property("Filter_Number",filter_number);
        Set_Property("G_cutoff",exp(-sqr(PI)/(4.*6.)));
        Set_Property("commutation_order",IPs.ExplicitFilters_IP.Commutation_Order);
        Set_Property("FGR",IPs.ExplicitFilters_IP.FGR[filter_number]);
        Set_Property("number_of_rings",IPs.ExplicitFilters_IP.Number_Of_Rings[filter_number]);
        Set_Property("filter_type",IPs.ExplicitFilters_IP.Filter_Type[filter_number]);
        Set_Property("target_filter_sharpness",IPs.ExplicitFilters_IP.Target_Filter_Sharpness[filter_number]);
        Set_Property("LS_constraints",IPs.ExplicitFilters_IP.LS_Constraints[filter_number]);
        Set_Property("Derivative_constraints",IPs.ExplicitFilters_IP.Derivative_Constraints[filter_number]);
        Set_Property("Filter_Width_strict",IPs.ExplicitFilters_IP.Filter_Width_Strict[filter_number]);
        Set_Property("memory_efficient",IPs.ExplicitFilters_IP.Filter_Memory_Efficient[filter_number]);
        Set_Property("relaxation_factor",IPs.ExplicitFilters_IP.Relaxation_Factor[filter_number]);
        Set_Property("least_squares_filter_weighting",IPs.ExplicitFilters_IP.Least_Squares_Filter_Weighting[filter_number]);
        Set_Property("least_squares_filter_weighting_factor",IPs.ExplicitFilters_IP.Least_Squares_Filter_Weighting_Factor[filter_number]);
        Set_Property("use_fixed_filter_width",!IPs.ExplicitFilters_IP.Filter_Relative[filter_number]);
        Set_Property("fixed_filter_width",IPs.ExplicitFilters_IP.Filter_Width[filter_number]);
        Set_Property("filter_strength",IPs.ExplicitFilters_IP.Filter_Strength[filter_number]);
        Set_Property("generate_at_startup",IPs.ExplicitFilters_IP.Generate_At_Startup[filter_number]);
        if (IPs.ExplicitFilters_IP.Reconstruction_Type[filter_number] == Explicit_Filter_Constants::DEFAULT_RECONSTRUCTION) {
            if( IPs.i_Reconstruction == RECONSTRUCTION_HIGH_ORDER ) {
                Set_Property("reconstruction_type",Explicit_Filter_Constants::CENO_RECONSTRUCTION);
            }
            else {
                Set_Property("reconstruction_type",Explicit_Filter_Constants::STANDARD_RECONSTRUCTION);
            }
        }
        else {
            Set_Property("reconstruction_type",IPs.ExplicitFilters_IP.Reconstruction_Type[filter_number]);
        }

        Set_Property("batch_flag",batch_flag);
        Set_Property("debug_flag",OFF);
        Set_Property("derivative_accuracy",3);
        Set_Property("number_of_rings_increased",3);
        Set_Property("derivative_accuracy",3);
        Set_Property("output_file_name",std::string(IPs.Output_File_Name_Prefix)+"_filter["+to_str(filter_number+1)+"]");
        Set_Property("finite_differencing_order",IPs.ExplicitFilters_IP.Finite_Differencing_Order);
        Set_Property("restarted",false);
        if (IPs.ExplicitFilters_IP.Filter_Type[filter_number]!=Explicit_Filter_Constants::RESTART_FILTER || batch_flag) {
            Set_Property("progress_mode",IPs.Progress_Mode);
        } else {
            Set_Property("progress_mode",PROGRESS_MODE_SILENT);
        }
        properties_changed = true;
    }
    
    template<typename T>
    void Set_Filter_Property(const std::string property_name, const T& property_value) {
        Set_Property(property_name, property_value);
        properties_changed = true;
    }
    template<typename T>
    void Set_Operating_Property(const std::string property_name, const T& property_value) {
        Set_Property(property_name, property_value);
    }
    
    void Output_Properties(void){
        for (map<std::string,string>::iterator p=property_map.begin(); p!=property_map.end(); p++){
            cout << "\n" << p->first << " = " << p->second; cout.flush();
        }
        cout << "\n" << "properties_changed = " << properties_changed; cout.flush();
    }
};

/* ------------------------------------------------------------------------------------------------------------------------------ */
/**
 * CLASS : Explicit_Filter_Adaptor
 */

#define SOLN_PSTATE_3D       0
#define SOLN_CSTATE_3D       1
#define SOLN_CSTATE_4D       2
#define SOLN_PSTATE_DOUBLE   3
#define SOLN_CSTATE_DOUBLE   4
#define FILTER_STATE_DOUBLE  5
#define FILTER_DOUBLE        6
#define COMMUTATION_ROWVECTOR     7
#define SOLN_PSTATE_3D_INDEXED 8
#define SOLN_CSTATE_3D_INDEXED 9

template <typename Soln_pState, typename Soln_cState>
class Explicit_Filter_Adaptor {
    
private:
    // adaptors
    typedef double (Soln_pState::*Soln_pState_member_ptr_type);
    typedef double (Soln_cState::*Soln_cState_member_ptr_type);
    typedef Soln_pState ***  (Hexa_Block<Soln_pState,Soln_cState>::*Soln_pState_3D_ptr_type);
    typedef Soln_cState ***  (Hexa_Block<Soln_pState,Soln_cState>::*Soln_cState_3D_ptr_type);
    typedef Soln_cState **** (Hexa_Block<Soln_pState,Soln_cState>::*Soln_cState_4D_ptr_type);
    typedef double (Filter_State::*Filter_State_member_ptr_type);
    static Soln_pState_member_ptr_type Soln_pState_member_ptr;
    static Soln_cState_member_ptr_type Soln_cState_member_ptr;
    static Soln_pState_3D_ptr_type Soln_pState_3D_ptr;
    static Soln_cState_3D_ptr_type Soln_cState_3D_ptr;
    static Soln_cState_4D_ptr_type Soln_cState_4D_ptr;
    static int dUdt_k_residual;
    static Filter_State_member_ptr_type Filter_State_member_ptr;
    static int startIdx, endIdx;

    // solution block
    static Hexa_Block<Soln_pState,Soln_cState> *Soln_Blk_ptr;
    static Filter_State ****Filter_State_Blk_ptr;
    static double ****Filter_double_Blk_ptr;
    static RowVector ***Commutation_RowVector_ptr;
    static Grid3D_Hexa_Block *Grid_Blk_ptr;
    static int nBlk;
    
    
public:
    static int adaptor_type;
    
    /* -------- Define all the adaptors ---------- */
    void Set_Adaptor(double Soln_pState::*&member);
    void Set_Adaptor(double Soln_cState::*&member);
    void Set_Adaptor(Soln_pState *** Hexa_Block<Soln_pState,Soln_cState>::*&member);
    void Set_Adaptor(Soln_pState *** Hexa_Block<Soln_pState,Soln_cState>::*&member, int startIdx, int endIdx);
    void Set_Adaptor(Soln_cState *** Hexa_Block<Soln_pState,Soln_cState>::*&member);
    void Set_Adaptor(Soln_cState *** Hexa_Block<Soln_pState,Soln_cState>::*&member, int startIdx, int endIdx);
    void Set_Adaptor(Soln_cState **** Hexa_Block<Soln_pState,Soln_cState>::*&member,int k_residual);
    void Set_Adaptor(double Filter_State::*&member);
    void Set_Adaptor(Filter_State ****);
    void Set_Adaptor(double ****);
    void Set_Adaptor(RowVector ***);
    void Set_Adaptor(int type);

    
    /* ------------------- Set the working solution block ---------------- */
    void Set_Solution_Block(Hexa_Block<Soln_pState,Soln_cState> &Soln_Blk);
    
    template <typename Filter_Variable_Type>
    void Set_Solution_Block(Filter_Variable_Type filter_variable, Grid3D_Hexa_Block &Grid_Blk, int n);
 
private:
    template <typename Filter_Variable_Type>
    void do_Set_Solution_Block(Filter_Variable_Type filter_variable, Grid3D_Hexa_Block &Grid_Blk, int n);
    void do_Set_Solution_Block(Filter_State ****Filter_Blk, Grid3D_Hexa_Block &Grid_Blk, int n);
    void do_Set_Solution_Block(double ****Filter_double_Blk, Grid3D_Hexa_Block &Grid_Blk, int n);

    /* --------------------------- Returns --------------------------------- */
    
    static void Asserts(void);
    static void Initialize_Result(RowVector& Result, int numvar);      

public:
    static RowVector FilterVariable(int i, int j, int k);
    static RowVector FilterVariable(Cell3D &theCell);
    
    static void Multiply_Weights_With_Values(RowVector& Result, RowVector& Weights, std::vector<Cell3D*> &stencil);

    /* -------------------- Replacement ------------------------ */
    void Load_into_Solution_Block(RowVector ***Filtered);
    
    /* ----------------- Just for commutation error calculations ------------------ */
    void Set_Commutation_RowVector(RowVector ***rows);
    
    void Set_Initial_Condition(Input_Parameters<Soln_pState,Soln_cState> &Input){
        int n = 8;
        int kind = 1;
        double b = 2.0;

        for(int i=0; i<Soln_Blk_ptr->NCi; i++) {
            for (int j=0; j<Soln_Blk_ptr->NCj; j++) {
                for (int k=0; k<Soln_Blk_ptr->NCk; k++) {
                    double x = Soln_Blk_ptr->Grid.Cell[i][j][k].Xc.x;
                    double y = Soln_Blk_ptr->Grid.Cell[i][j][k].Xc.y;
                    double z = Soln_Blk_ptr->Grid.Cell[i][j][k].Xc.z;
                    //double Dx = Input.Grid_IP.Box_Width;
                    //double Dy = Input.Grid_IP.Box_Height;
                    //double Dz = Input.Grid_IP.Box_Length;
                    double r = sqrt(sqr(x) + sqr(y) + sqr(z));
                    //Soln_Blk_ptr->W[i][j][k].*Soln_pState_member_ptr = chebyshev_polynomial(kind,n,r/scaling);
                    Soln_Blk_ptr->W[i][j][k].*Soln_pState_member_ptr = cos(b*x);
                                                                                            
                }
            }
        }
    }
    
    RowVector Exact_Derivative(Input_Parameters<Soln_pState,Soln_cState> &Input, double x, double y, double z){
        int n = 8;
        int kind = 1;
        double b = 2.0;
        //double Dx = Input.Grid_IP.Box_Width;
//        double Dy = Input.Grid_IP.Box_Height;
//        double Dz = Input.Grid_IP.Box_Length;
//        double scaling = sqrt(3.);
        double r = sqrt(sqr(x) + sqr(y) + sqr(z));
        RowVector temp(1);
        //temp(0) = 1./scaling*chebyshev_polynomial_derivative(kind,n,r);
        temp(0) = -b*sin(b*x);
        return temp;
    }
    
};
    


template <typename Soln_pState, typename Soln_cState>
void Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::Set_Adaptor(double Soln_pState::*&member) {
    Soln_pState_member_ptr = member;
    adaptor_type = SOLN_PSTATE_DOUBLE;
}

template <typename Soln_pState, typename Soln_cState>
void Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::Set_Adaptor(double Soln_cState::*&member) {
    Soln_pState_member_ptr = member;
    adaptor_type = SOLN_CSTATE_DOUBLE;
}

template <typename Soln_pState, typename Soln_cState>
void Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::Set_Adaptor(Soln_pState *** Hexa_Block<Soln_pState,Soln_cState>::*&member) {
    Soln_pState_3D_ptr = member;
    adaptor_type = SOLN_PSTATE_3D;
}

template <typename Soln_pState, typename Soln_cState>
void Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::Set_Adaptor(Soln_cState *** Hexa_Block<Soln_pState,Soln_cState>::*&member) {
    Soln_cState_3D_ptr = member;
    adaptor_type = SOLN_CSTATE_3D;
}

template <typename Soln_pState, typename Soln_cState>
void Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::Set_Adaptor(Soln_pState *** Hexa_Block<Soln_pState,Soln_cState>::*&member, int startIndex, int endIndex) {
    /* startIndex = 1 --> first element   *
     * endIndex = NumVar --> last element */
    startIdx = startIndex;
    endIdx = endIndex;
    Soln_pState_3D_ptr = member;
    adaptor_type = SOLN_PSTATE_3D_INDEXED;
}

template <typename Soln_pState, typename Soln_cState>
void Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::Set_Adaptor(Soln_cState *** Hexa_Block<Soln_pState,Soln_cState>::*&member, int startIndex, int endIndex) {
    /* startIndex = 1 --> first element   *
     * endIndex = NumVar --> last element */
    startIdx = startIndex;
    endIdx = endIndex;
    Soln_cState_3D_ptr = member;
    adaptor_type = SOLN_CSTATE_3D_INDEXED;
}

template <typename Soln_pState, typename Soln_cState>
void Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::Set_Adaptor(Soln_cState **** Hexa_Block<Soln_pState,Soln_cState>::*&member,int k_residual) {
    Soln_cState_4D_ptr = member;
    dUdt_k_residual = k_residual;
    adaptor_type = SOLN_CSTATE_4D;
}

template <typename Soln_pState, typename Soln_cState>
void Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::Set_Adaptor(double Filter_State::*&member) {
    adaptor_type = -1; // not defined
}

template <typename Soln_pState, typename Soln_cState>
void Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::Set_Adaptor(Filter_State ****F) {
    adaptor_type = FILTER_STATE_DOUBLE;
}

template <typename Soln_pState, typename Soln_cState>
void Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::Set_Adaptor(double ****F) {
    adaptor_type = FILTER_DOUBLE;
}

template <typename Soln_pState, typename Soln_cState>
void Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::Set_Adaptor(int type) {
    adaptor_type = type;
}

/* ------------------- Set the working solution block ---------------- */
template <typename Soln_pState, typename Soln_cState>
void Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::Set_Solution_Block(Hexa_Block<Soln_pState,Soln_cState> &Soln_Blk) {
    Soln_Blk_ptr = &Soln_Blk;
}

template <typename Soln_pState, typename Soln_cState>
template <typename Filter_Variable_Type>
void Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::Set_Solution_Block(Filter_Variable_Type filter_variable, Grid3D_Hexa_Block &Grid_Blk, int n){
    // This is a workaround to have specialization of templated functions in templated classes
    this->do_Set_Solution_Block(filter_variable,Grid_Blk, n);
}

template <typename Soln_pState, typename Soln_cState>
template <typename Filter_Variable_Type>
void Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::do_Set_Solution_Block(Filter_Variable_Type filter_variable, Grid3D_Hexa_Block &Grid_Blk, int n){
    cout << "This should never be called" << endl;
    exit(1);
}

template <typename Soln_pState, typename Soln_cState>
void Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::do_Set_Solution_Block(Filter_State ****Filter_Blk, Grid3D_Hexa_Block &Grid_Blk, int n){
    Filter_State_Blk_ptr = Filter_Blk;
    Grid_Blk_ptr = &Grid_Blk;
    nBlk = n;
}

template <typename Soln_pState, typename Soln_cState>
void Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::do_Set_Solution_Block(double ****Filter_double_Blk, Grid3D_Hexa_Block &Grid_Blk, int n){
    Filter_double_Blk_ptr = Filter_double_Blk;
    Grid_Blk_ptr = &Grid_Blk;
    nBlk = n;
}

/* ----------------- Just for commutation error calculations ------------------ */

template <typename Soln_pState, typename Soln_cState>
void Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::Set_Commutation_RowVector(RowVector ***rows) {
    adaptor_type = COMMUTATION_ROWVECTOR;
    Commutation_RowVector_ptr = rows;
}

/* --------------------------- Returns --------------------------------- */

template <typename Soln_pState, typename Soln_cState>
void Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::Asserts(void) {
    switch (adaptor_type) {
        case SOLN_PSTATE_DOUBLE :
        case SOLN_CSTATE_DOUBLE :
        case SOLN_PSTATE_3D :
        case SOLN_CSTATE_3D :
        case SOLN_CSTATE_4D :
            assert(Soln_Blk_ptr != NULL);
            break;
        case FILTER_STATE_DOUBLE :
            assert(Filter_State_Blk_ptr != NULL);
            break;
        case FILTER_DOUBLE :
            assert(Filter_double_Blk_ptr != NULL);
            break;
        case COMMUTATION_ROWVECTOR :
            assert(Commutation_RowVector_ptr != NULL);
            break;
    }
}

template <typename Soln_pState, typename Soln_cState>
RowVector Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::FilterVariable(int i, int j, int k) {
    Asserts();
    switch (adaptor_type) {
        case SOLN_PSTATE_DOUBLE : {
            RowVector temp(1);
            temp(0) = Soln_Blk_ptr->W[i][j][k].*Soln_pState_member_ptr;
            return temp;
        }
        case SOLN_CSTATE_DOUBLE : {
            RowVector temp(1);
            temp(0) = Soln_Blk_ptr->U[i][j][k].*Soln_cState_member_ptr;
            return temp;
        }
        case SOLN_PSTATE_3D: {
            int numvar = Soln_Blk_ptr->NumVar();
            RowVector temp(numvar);
            for (int n=0; n<numvar; n++) {
                temp(n) = (Soln_Blk_ptr->*Soln_pState_3D_ptr)[i][j][k][n+1];
            }
            return temp;
        }
        case SOLN_CSTATE_3D: {
            int numvar = Soln_Blk_ptr->NumVar();
            RowVector temp(numvar);
            for (int n=0; n<numvar; n++) {
                temp(n) = (Soln_Blk_ptr->*Soln_cState_3D_ptr)[i][j][k][n+1];
            }
            return temp;
        }
        case SOLN_PSTATE_3D_INDEXED: {
            RowVector temp(endIdx-startIdx+1);
            for (int n=startIdx; n<=endIdx; n++) {
                temp(n-startIdx) = (Soln_Blk_ptr->*Soln_pState_3D_ptr)[i][j][k][n];
            }
            return temp;
        }
        case SOLN_CSTATE_3D_INDEXED: {
            RowVector temp(endIdx-startIdx+1);
            for (int n=startIdx; n<=endIdx; n++) {
                temp(n-startIdx) = (Soln_Blk_ptr->*Soln_cState_3D_ptr)[i][j][k][n];
            }
            return temp;
        }
        case SOLN_CSTATE_4D: {
            int numvar = Soln_Blk_ptr->NumVar();
            RowVector temp(numvar);
            for (int n=0; n<numvar; n++) {
                temp(n) = (Soln_Blk_ptr->*Soln_cState_4D_ptr)[i][j][k][dUdt_k_residual][n+1];
            }
            return temp;
        }
        case FILTER_STATE_DOUBLE: {
            RowVector temp(1);
            temp(0) = Filter_State_Blk_ptr[nBlk][i][j][k].member;
            return temp;
        }
        case FILTER_DOUBLE: {
            RowVector temp(1);
            temp(0) = Filter_double_Blk_ptr[nBlk][i][j][k];
            return temp;
        }
        case COMMUTATION_ROWVECTOR:
            return Commutation_RowVector_ptr[i][j][k];
    }
}
template <typename Soln_pState, typename Soln_cState>
RowVector Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::FilterVariable(Cell3D &theCell) {
    return FilterVariable(theCell.I,theCell.J,theCell.K);
}

template <typename Soln_pState, typename Soln_cState>
void Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::Multiply_Weights_With_Values(RowVector& Result, RowVector& Weights, std::vector<Cell3D*> &stencil) {
    Asserts();
    // the dimensions of the matrix
    int number_of_neighbours = stencil.size();
    int I, J, K, cell;

    switch (adaptor_type) {
        case SOLN_PSTATE_DOUBLE :
            Initialize_Result(Result,1);
            for (cell=0; cell<number_of_neighbours; cell++) {
                I = stencil[cell]->I;
                J = stencil[cell]->J;
                K = stencil[cell]->K;
                Result(0) += Weights[cell] * Soln_Blk_ptr->W[I][J][K].*Soln_pState_member_ptr;
            }
            break;
        case SOLN_CSTATE_DOUBLE :
            Initialize_Result(Result,1);
            for (cell=0; cell<number_of_neighbours; cell++) {
                I = stencil[cell]->I;
                J = stencil[cell]->J;
                K = stencil[cell]->K;
                Result(0) += Weights[cell] * Soln_Blk_ptr->U[I][J][K].*Soln_cState_member_ptr;
            }
            break;
        case SOLN_PSTATE_3D: {
            int numvar = Soln_Blk_ptr->NumVar();
            Initialize_Result(Result,numvar);
            for (cell=0; cell<number_of_neighbours; cell++) {
                I = stencil[cell]->I;
                J = stencil[cell]->J;
                K = stencil[cell]->K;
                for (int n=0; n<numvar; n++) {
                    Result(n) += Weights[cell] * (Soln_Blk_ptr->*Soln_pState_3D_ptr)[I][J][K][n+1];
                }
            }
            break;  
        }
        case SOLN_CSTATE_3D: {
            int numvar = Soln_Blk_ptr->NumVar();
            Initialize_Result(Result,numvar);
            for (cell=0; cell<number_of_neighbours; cell++) {
                I = stencil[cell]->I;
                J = stencil[cell]->J;
                K = stencil[cell]->K;
                for (int n=0; n<numvar; n++) {
                    Result(n) += Weights[cell] * (Soln_Blk_ptr->*Soln_cState_3D_ptr)[I][J][K][n+1];
                }
            }
            break;
        }
        case SOLN_PSTATE_3D_INDEXED: {
            Initialize_Result(Result,endIdx-startIdx+1);
            for (cell=0; cell<number_of_neighbours; cell++) {
                I = stencil[cell]->I;
                J = stencil[cell]->J;
                K = stencil[cell]->K;
                for (int n=startIdx; n<=endIdx; n++) {
                    Result(n-startIdx) += Weights[cell] * (Soln_Blk_ptr->*Soln_pState_3D_ptr)[I][J][K][n];
                }
            }
            break;
        }
        case SOLN_CSTATE_3D_INDEXED: {
            Initialize_Result(Result,endIdx-startIdx+1);
            for (cell=0; cell<number_of_neighbours; cell++) {
                I = stencil[cell]->I;
                J = stencil[cell]->J;
                K = stencil[cell]->K;
                for (int n=startIdx; n<=endIdx; n++) {
                    Result(n-startIdx) += Weights[cell] * (Soln_Blk_ptr->*Soln_cState_3D_ptr)[I][J][K][n];
                }
            }
            break;
        }
        case SOLN_CSTATE_4D: {
            int numvar = Soln_Blk_ptr->NumVar();
            Initialize_Result(Result,numvar-1);
            for (cell=0; cell<number_of_neighbours; cell++) {
                I = stencil[cell]->I;
                J = stencil[cell]->J;
                K = stencil[cell]->K;
                for (int n=0; n<numvar-1; n++) {  // don't filter rho
                    Result(n) += Weights[cell] * (Soln_Blk_ptr->*Soln_cState_4D_ptr)[I][J][K][dUdt_k_residual][n+2];
                }
            }
            break;
        }
        case FILTER_STATE_DOUBLE:
            Initialize_Result(Result,1);
            for (cell=0; cell<number_of_neighbours; cell++) {
                I = stencil[cell]->I;
                J = stencil[cell]->J;
                K = stencil[cell]->K;
                Result(0) += Weights[cell] * Filter_State_Blk_ptr[nBlk][I][J][K].member;
            }
            break;
        case FILTER_DOUBLE:
            Initialize_Result(Result,1);
            for (cell=0; cell<number_of_neighbours; cell++) {
                I = stencil[cell]->I;
                J = stencil[cell]->J;
                K = stencil[cell]->K;
                Result(0) += Weights[cell] * Filter_double_Blk_ptr[nBlk][I][J][K];
            }
            break;
        case COMMUTATION_ROWVECTOR:
            I = stencil[0]->I;
            J = stencil[0]->J;
            K = stencil[0]->K;
            int numvar = Commutation_RowVector_ptr[I][J][K].size();
            Initialize_Result(Result,numvar);
            for (cell=0; cell<number_of_neighbours; cell++) {
                I = stencil[cell]->I;
                J = stencil[cell]->J;
                K = stencil[cell]->K;
                Result += Weights[cell] * Commutation_RowVector_ptr[I][J][K];
            }
            break;
    }
}

template <typename Soln_pState, typename Soln_cState>
void Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::Initialize_Result(RowVector& Result, int numvar) {
    if (Result.size() != numvar){
        Result.newsize(numvar);
    }
    Result.zero();
}


/* ---------------- Replace unfiltered variables with filtered variables --------------- */
template <typename Soln_pState, typename Soln_cState>
void Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::Load_into_Solution_Block(RowVector ***Filtered) {
    if (adaptor_type != FILTER_STATE_DOUBLE && 
        adaptor_type != FILTER_DOUBLE) {
        for(int i=Soln_Blk_ptr->ICl; i<=Soln_Blk_ptr->ICu; i++) {
            for (int j=Soln_Blk_ptr->JCl; j<=Soln_Blk_ptr->JCu; j++) {
                for (int k=Soln_Blk_ptr->KCl; k<=Soln_Blk_ptr->KCu; k++) {
                    switch (adaptor_type) {
                        case SOLN_PSTATE_DOUBLE:
                            Soln_Blk_ptr->W[i][j][k].*Soln_pState_member_ptr = Filtered[i][j][k](0);
                            break;
                        case SOLN_CSTATE_DOUBLE:
                            Soln_Blk_ptr->U[i][j][k].*Soln_cState_member_ptr = Filtered[i][j][k](0);
                            break;
                        case SOLN_PSTATE_3D:
                            for (int n=1; n<=Soln_Blk_ptr->NumVar(); n++)
                                (Soln_Blk_ptr->*Soln_pState_3D_ptr)[i][j][k][n] = Filtered[i][j][k](n-1);
                            break;
                        case SOLN_CSTATE_3D:
                            for (int n=1; n<=Soln_Blk_ptr->NumVar(); n++)
                                (Soln_Blk_ptr->*Soln_cState_3D_ptr)[i][j][k][n] = Filtered[i][j][k](n-1);
                            break;
                        case SOLN_PSTATE_3D_INDEXED: 
                            for (int n=startIdx; n<=endIdx; n++) 
                                (Soln_Blk_ptr->*Soln_pState_3D_ptr)[i][j][k][n] = Filtered[i][j][k](n-startIdx);
                            break;
                        case SOLN_CSTATE_3D_INDEXED: 
                            for (int n=startIdx; n<=endIdx; n++) 
                                (Soln_Blk_ptr->*Soln_cState_3D_ptr)[i][j][k][n] = Filtered[i][j][k](n-startIdx);
                            break;
                        case SOLN_CSTATE_4D:
                            for (int n=2; n<=Soln_Blk_ptr->NumVar(); n++)  // don't fill rho
                                (Soln_Blk_ptr->*Soln_cState_4D_ptr)[i][j][k][dUdt_k_residual][n] = Filtered[i][j][k](n-2);
                            break;
                    }
                }
            }
        }        
    } else {
        for(int i=Grid_Blk_ptr->ICl; i<=Grid_Blk_ptr->ICu; i++) {
            for (int j=Grid_Blk_ptr->JCl; j<=Grid_Blk_ptr->JCu; j++) {
                for (int k=Grid_Blk_ptr->KCl; k<=Grid_Blk_ptr->KCu; k++) {
                    switch (adaptor_type) {
                        case FILTER_STATE_DOUBLE:
                            Filter_State_Blk_ptr[nBlk][i][j][k].member = Filtered[i][j][k](0);
                            break;
                        case FILTER_DOUBLE:
                            Filter_double_Blk_ptr[nBlk][i][j][k] = Filtered[i][j][k](0);
                            break;
                    }
                }
            }
        }        
    }
}
/* --------------------------- Statics ----------------------------------- */
template<typename Soln_pState, typename Soln_cState>
int Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::dUdt_k_residual = 0;

template<typename Soln_pState, typename Soln_cState>
int Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::startIdx = 1;

template<typename Soln_pState, typename Soln_cState>
int Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::endIdx = 1;

template<typename Soln_pState, typename Soln_cState>
double Soln_pState::* Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::Soln_pState_member_ptr = NULL;

template<typename Soln_pState, typename Soln_cState>
double Soln_cState::* Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::Soln_cState_member_ptr = NULL;

template<typename Soln_pState, typename Soln_cState>
Soln_pState *** Hexa_Block<Soln_pState,Soln_cState>::* Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::Soln_pState_3D_ptr = NULL;

template<typename Soln_pState, typename Soln_cState>
Soln_cState *** Hexa_Block<Soln_pState,Soln_cState>::* Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::Soln_cState_3D_ptr = NULL;

template<typename Soln_pState, typename Soln_cState>
Soln_cState **** Hexa_Block<Soln_pState,Soln_cState>::* Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::Soln_cState_4D_ptr = NULL;

template<typename Soln_pState, typename Soln_cState>
double Filter_State::* Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::Filter_State_member_ptr = NULL;

template<typename Soln_pState, typename Soln_cState>
Hexa_Block<Soln_pState,Soln_cState> * Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::Soln_Blk_ptr = NULL;

template<typename Soln_pState, typename Soln_cState>
Filter_State **** Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::Filter_State_Blk_ptr = NULL;

template<typename Soln_pState, typename Soln_cState>
double **** Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::Filter_double_Blk_ptr = NULL;

template<typename Soln_pState, typename Soln_cState>
RowVector *** Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::Commutation_RowVector_ptr = NULL;

template<typename Soln_pState, typename Soln_cState>
Grid3D_Hexa_Block * Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::Grid_Blk_ptr = NULL;

template<typename Soln_pState, typename Soln_cState>
int Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::nBlk = 0;


template<typename Soln_pState, typename Soln_cState>
int Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::adaptor_type = SOLN_PSTATE_DOUBLE;


#endif

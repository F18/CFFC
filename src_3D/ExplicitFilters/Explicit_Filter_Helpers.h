/*
 *  Explicit_Filter_Properties.h
 *  CFFC
 *
 *  Created by Willem Deconinck on 25/03/08.
 *
 */

#ifndef _EXPLICIT_FILTER_HELPERS_INCLUDED
#define _EXPLICIT_FILTER_HELPERS_INCLUDED

#ifndef _CFD_INCLUDED
#include "../CFD/CFD.h"
#endif

#ifndef _HEXA_BLOCK_INCLUDED
#include "../HexaBlock/HexaBlock.h"
#endif

#ifndef _NEIGHBOURS_INCLUDED
#include "Neighbours.h"
#endif

#ifndef _FILTER_STATE_INCLUDED
#include "Filter_State.h"
#endif

#include <map>

/* ------------------------------------------------------------------------------------------------------------------------------ */
/**
 * CLASS : Explicit_Filter_Properties
 */


class Explicit_Filter_Properties {
public:
    std::map<string,int> map_int;
    std::map<string,double> map_double;
    std::map<string,char*> map_char;
    
public:
    
    static bool   debug_flag;
    static int    commutation_order;
    static int    finite_differencing_order;
    static double FGR;
    static int    number_of_rings;
    static double target_filter_sharpness;
    static bool   LS_constraints;
    static bool   Filter_Width_strict;
    static int    Derivative_constraints;
    static bool   Memory_Efficient;
    static bool   restarted;
    static int    progress_mode;
    static int    filter_type;
    static char   *output_file_name;
    static int    batch_flag;
    static int    number_of_rings_increased;
    static int    derivative_accuracy;
    static double G_cutoff;
    static double relaxation_factor;
    static int    least_squares_filter_weighting;
    static double least_squares_filter_weighting_factor;
    static bool   uniform_grid;
    Explicit_Filter_Properties(void){
        Set_New_Property("commutation_order",2);
        Set_New_Property("debug_flag",OFF);
        Set_New_Property("FGR",2.0);
        Set_New_Property("number_of_rings",2);
        Set_New_Property("memory_efficient",OFF);
        Set_New_Property("finite_differencing_order",4);
        Set_New_Property("LS_constraints",ON);
        Set_New_Property("Derivative_constraints",DEFAULT);
        Set_New_Property("Filter_Width_strict",OFF);
        Set_New_Property("target_filter_sharpness",-1.0);
        Set_New_Property("restarted",false);
        Set_New_Property("progress_mode",PROGRESS_MODE_MESSAGE);
        Set_New_Property("filter_type",FILTER_TYPE_VASILYEV);
        Set_New_Property("output_file_name","output_file_name");
        Set_New_Property("batch_flag",0);
        Set_New_Property("number_of_rings_increased",3);
        Set_New_Property("derivative_accuracy",3);
        Set_New_Property("G_cutoff",exp(-sqr(PI)/(4.*6.)));
        Set_New_Property("least_squares_filter_weighting",ON);
        Set_New_Property("relaxation_factor",double(DEFAULT));
        Set_New_Property("least_squares_filter_weighting_factor",double(DEFAULT));
        Set_New_Property("uniform_grid",ON);
    }
    
//    ~Explicit_Filter_Properties(void){
//        map_int.clear();
//        map_double.clear();
//        map_char.clear();
////        map_int.erase(map_int.begin(),map_int.end());
////        map_double.erase(map_double.begin(),map_double.end());
////        map_char.erase(map_char.begin(),map_char.end());
//    }
    
    template <typename Soln_pState, typename Soln_cState>
    void Set_Properties(Input_Parameters<Soln_pState,Soln_cState> &IPs, int batch_flag){
        Set_Property("output_file_name",IPs.Output_File_Name_Prefix);
        Set_Property("FGR",IPs.Turbulence_IP.FGR);
        Set_Property("commutation_order",IPs.Turbulence_IP.commutation_order);
        Set_Property("finite_differencing_order",IPs.Turbulence_IP.finite_differencing_order);
        Set_Property("number_of_rings",IPs.Turbulence_IP.number_of_rings);
        Set_Property("filter_type",IPs.Turbulence_IP.i_filter_type);
        Set_Property("target_filter_sharpness",IPs.Turbulence_IP.Target_Filter_Sharpness);
        Set_Property("LS_constraints",IPs.Turbulence_IP.LS_constraints);
        Set_Property("Derivative_constraints",IPs.Turbulence_IP.Derivative_constraints);
        Set_Property("Filter_Width_strict",IPs.Turbulence_IP.Filter_Width_strict);
        Set_Property("memory_efficient",IPs.Turbulence_IP.Filter_Memory_Efficient);
        Set_Property("relaxation_factor",IPs.Turbulence_IP.relaxation_factor);
        Set_Property("least_squares_filter_weighting",IPs.Turbulence_IP.least_squares_filter_weighting);
        Set_Property("least_squares_filter_weighting_factor",IPs.Turbulence_IP.least_squares_filter_weighting_factor);
        Set_Property("uniform_grid",IPs.Turbulence_IP.uniform_grid);
        Set_Property("batch_flag",batch_flag);
        if (IPs.Turbulence_IP.i_filter_type!=FILTER_TYPE_RESTART) {
            Set_Property("progress_mode",IPs.Progress_Mode);
        } else {
            Set_Property("progress_mode",PROGRESS_MODE_SILENT);
        }
    }

private:
    void Set_New_Property(string property_name, int property_value) {
        map_int[property_name]=property_value;
    }
    void Set_New_Property(string property_name, bool property_value) {
        map_int[property_name]=property_value;
    }
    void Set_New_Property(string property_name, double property_value) {
        map_double[property_name]=property_value;
    }
    void Set_New_Property(string property_name, char* property_value) {
        map_char[property_name]=property_value;
    }
    
public:
    void Set_Property(string property_name, int property_value) {
        map<string,int>::iterator p;
        p = map_int.find(property_name);
        if(p!=map_int.end()) {
            map_int[property_name]=property_value;
        } else {
            cout << "\n ERROR: Property " << property_name << " does not exist in Explicit_Filter_Helpers.h"; cout.flush();
        }
    }
    void Set_Property(string property_name, bool property_value) {
        map<string,int>::iterator p;
        p = map_int.find(property_name);
        if(p!=map_int.end()) {
            map_int[property_name]=property_value;
        } else {
            cout << "\n ERROR: Property " << property_name << " does not exist in Explicit_Filter_Helpers.h"; cout.flush();
        }
    }
    void Set_Property(string property_name, double property_value) {
        map<string,double>::iterator p;
        p = map_double.find(property_name);
        if(p!=map_double.end()) {
            map_double[property_name]=property_value;
        } else {
            cout << "\n ERROR: Property " << property_name << " does not exist in Explicit_Filter_Helpers.h"; cout.flush();
        }
        
    }
    void Set_Property(string property_name, char* property_value) {
        map<string,char*>::iterator p;
        p = map_char.find(property_name);
        if(p!=map_char.end()) {
            map_char[property_name]=property_value;
        } else {
            cout << "\n ERROR: Property " << property_name << " does not exist in Explicit_Filter_Helpers.h"; cout.flush();
        }
    }
    
    int Get_Property_int(string property_name) {
        map<string,int>::iterator p;
        p = map_int.find(property_name);
        assert(p!=map_int.end());
        return p->second;
    }
    double Get_Property_double(string property_name) {
        map<string,double>::iterator p;
        p = map_double.find(property_name);
        assert(p!=map_double.end());
    return p->second;    
    }
    char* Get_Property_char(string property_name) {
        map<string,char*>::iterator p;
        p = map_char.find(property_name);
        assert(p!=map_char.end());
    return p->second;   
    }
    
    void Output_Properties(void){
        for (map<string,int>::iterator p=map_int.begin(); p!=map_int.end(); p++){
            cout << "\n" << p->first << " = " << p->second; cout.flush();
        }
        for (map<string,double>::iterator p=map_double.begin(); p!=map_double.end(); p++){
            cout << "\n" << p->first << " = " << p->second; cout.flush();
        }
        for (map<string,char*>::iterator p=map_char.begin(); p!=map_char.end(); p++){
            cout << "\n" << p->first << " = " << p->second; cout.flush();
        }
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
    void Set_Adaptor(Soln_cState *** Hexa_Block<Soln_pState,Soln_cState>::*&member);
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
public:
    static RowVector FilterVariable(int i, int j, int k);
    static RowVector FilterVariable(Cell3D &theCell);
    
    static void FillRowVector(RowVector &row, int i, int j, int k);
    static void FillRowVector(RowVector &row, Cell3D &theCell);
    
    static void FillMatrixRow(DenseMatrix &matrix, int row_index, int i, int j, int k);
    static void FillMatrixRow(DenseMatrix &matrix, int row_index, Cell3D &theCell);
    
    static void FillMatrix(DenseMatrix &matrix, Neighbours &theNeighbours);
        
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
    
    double chebyshev_polynomial(const int& r, const int& n, const double& x)
    {
        double result;
        int i;
        double a;
        double b;
        
        
        //
        // Prepare A and B
        //
        if( r==1 ) {
            a = 1;
            b = x;
        }
        else {
            a = 1;
            b = 2*x;
        }
        
        //
        // Special cases: N=0 or N=1
        //
        if( n==0 ) {
            result = a;
            return result;
        }
        if( n==1 ) {
            result = b;
            return result;
        }
        
        //
        // General case: N>=2
        //
        for(i = 2; i <= n; i++) {
            result = 2*x*b-a;
            a = b;
            b = result;
        }
        return result;
    }
    
    
    double chebyshev_polynomial_derivative(const int& r, const int& n, const double& x) {
        double result;
        int i;
        double T_a , T_b, T;
        double dT_a , dT_b, dT;
        
        
        //
        // Prepare A and B
        //
        if( r==1 ) {
            T_a = 1;
            T_b = x;
            dT_a = 0;
            dT_b = 1;
        } else {
            T_a = 1;
            T_b = 2*x;
            dT_a = 0;
            dT_b = 2;
        }
        
        //
        // Special cases: N=0 or N=1
        //
        if( n==0 ) {
            result = dT_a;
            return result;
        }
        if( n==1 ) {
            result = dT_b;
            return result;
        }
        
        //
        // General case: N>=2
        //
        for(i = 2; i <= n; i++) {
            dT = 2*x*dT_b - dT_a + 2*T_b;
            dT_a = dT_b;
            dT_b = dT;
            
            T = 2*x*T_b-T_a;
            T_a = T_b;
            T_b = T;
        }
        return dT;
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
    Soln_cState_3D_ptr = member;
    adaptor_type = SOLN_PSTATE_3D;
}

template <typename Soln_pState, typename Soln_cState>
void Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::Set_Adaptor(Soln_cState *** Hexa_Block<Soln_pState,Soln_cState>::*&member) {
    Soln_cState_3D_ptr = member;
    adaptor_type = SOLN_CSTATE_3D;
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
                temp(n) = (Soln_Blk_ptr->*Soln_cState_3D_ptr)[i][j][k][n+1];
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
void Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::FillRowVector(RowVector &row, int i, int j, int k) {
    Asserts();
    switch (adaptor_type) {
        case SOLN_PSTATE_DOUBLE :
            row(0) = Soln_Blk_ptr->W[i][j][k].*Soln_pState_member_ptr;
            break;
        case SOLN_CSTATE_DOUBLE :
            row(0) = Soln_Blk_ptr->U[i][j][k].*Soln_cState_member_ptr;
            break;
        case SOLN_PSTATE_3D: {
            int numvar = Soln_Blk_ptr->NumVar();
            for (int n=0; n<numvar; n++) {
                row(n) = (Soln_Blk_ptr->*Soln_cState_3D_ptr)[i][j][k][n+1];
            }
            break;  
        }
        case SOLN_CSTATE_3D: {
            int numvar = Soln_Blk_ptr->NumVar();
            for (int n=0; n<numvar; n++) {
                row(n) = (Soln_Blk_ptr->*Soln_cState_3D_ptr)[i][j][k][n+1];
            }
            break;
        }
        case SOLN_CSTATE_4D: {
            int numvar = Soln_Blk_ptr->NumVar();
            for (int n=0; n<numvar; n++) {
                row(n) = (Soln_Blk_ptr->*Soln_cState_4D_ptr)[i][j][k][dUdt_k_residual][n+1];
            }
            break;
        }
        case FILTER_STATE_DOUBLE:
            row(0) = Filter_State_Blk_ptr[nBlk][i][j][k].member;
            break;
        case FILTER_DOUBLE:
            row(0) = Filter_double_Blk_ptr[nBlk][i][j][k];
            break;
        case COMMUTATION_ROWVECTOR:
            row = Commutation_RowVector_ptr[i][j][k];
            break;
    }
}
template <typename Soln_pState, typename Soln_cState>
void Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::FillRowVector(RowVector &row, Cell3D &theCell) {
    FillRowVector(row,theCell.I,theCell.J,theCell.K);
}

template <typename Soln_pState, typename Soln_cState>
void Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::FillMatrixRow(DenseMatrix &matrix, int row_index, int i, int j, int k) {
    Asserts();
    switch (adaptor_type) {
        case SOLN_PSTATE_DOUBLE :
            matrix(row_index,0) = Soln_Blk_ptr->W[i][j][k].*Soln_pState_member_ptr;
            break;
        case SOLN_CSTATE_DOUBLE :
            matrix(row_index,0) = Soln_Blk_ptr->U[i][j][k].*Soln_cState_member_ptr;
            break;
        case SOLN_PSTATE_3D: {
            int numvar = Soln_Blk_ptr->NumVar();
            for (int n=0; n<numvar; n++) {
                matrix(row_index,n) = (Soln_Blk_ptr->*Soln_cState_3D_ptr)[i][j][k][n+1];
            }
            break;
        }
        case SOLN_CSTATE_3D: {
            int numvar = Soln_Blk_ptr->NumVar();
            for (int n=0; n<numvar; n++) {
                matrix(row_index,n) = (Soln_Blk_ptr->*Soln_cState_3D_ptr)[i][j][k][n+1];
            }
            break;
        }
        case SOLN_CSTATE_4D: {
            int numvar = Soln_Blk_ptr->NumVar();
            for (int n=0; n<numvar; n++) {
                matrix(row_index,n) = (Soln_Blk_ptr->*Soln_cState_4D_ptr)[i][j][k][dUdt_k_residual][n+1];
            }
            break;
        }
        case FILTER_STATE_DOUBLE:
            matrix(row_index,0) = Filter_State_Blk_ptr[nBlk][i][j][k].member;
            break;
        case FILTER_DOUBLE:
            matrix(row_index,0) = Filter_double_Blk_ptr[nBlk][i][j][k];
            break;
        case COMMUTATION_ROWVECTOR:
            matrix.assignRow(row_index,Commutation_RowVector_ptr[i][j][k]);
            break;
    }
}
template <typename Soln_pState, typename Soln_cState>
void Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::FillMatrixRow(DenseMatrix &matrix, int row_index, Cell3D &theCell) {
    FillMatrixRow(matrix,row_index,theCell.I,theCell.J,theCell.K);
}

template <typename Soln_pState, typename Soln_cState>
void Explicit_Filter_Adaptor<Soln_pState,Soln_cState>::FillMatrix(DenseMatrix &matrix, Neighbours &theNeighbours) {
    Asserts();
    
    // the dimensions of the matrix
    int number_of_neighbours = theNeighbours.number_of_neighbours;
    int number_of_variables;
    switch (adaptor_type) {
        case SOLN_PSTATE_DOUBLE :
        case SOLN_CSTATE_DOUBLE :
        case FILTER_STATE_DOUBLE :
        case FILTER_DOUBLE :
            number_of_variables = 1;
            break;
        case SOLN_PSTATE_3D:
        case SOLN_CSTATE_3D:
        case SOLN_CSTATE_4D:
            number_of_variables = Soln_Blk_ptr->NumVar();
            break;
        case COMMUTATION_ROWVECTOR: {
            int I = theNeighbours.neighbour[0].I;
            int J = theNeighbours.neighbour[0].J;
            int K = theNeighbours.neighbour[0].K;
            number_of_variables = Commutation_RowVector_ptr[I][J][K].size();
        }
    }
    
    // Set dimensions of the matrix
    if (matrix.get_n()!=number_of_neighbours || matrix.get_m()!=number_of_variables) {
        matrix.newsize(number_of_neighbours,number_of_variables);
    }
    
    // Fill the matrix
    for (int i=0; i<number_of_neighbours; i++){
        FillMatrixRow(matrix, i, theNeighbours.neighbour[i]);
    }
    
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
                        case SOLN_CSTATE_4D:
                            for (int n=1; n<=Soln_Blk_ptr->NumVar(); n++)
                                (Soln_Blk_ptr->*Soln_cState_4D_ptr)[i][j][k][dUdt_k_residual][n] = Filtered[i][j][k](n-1);
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

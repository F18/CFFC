/*
 *  Vasilyev_Filter.h
 *  CFFC
 *
 *  Created by Willem Deconinck on 27/04/08.
 *
 */

#ifndef _VASILYEV_FILTER_INCLUDED
#define _VASILYEV_FILTER_INCLUDED

#include "Discrete_Filter.h"


#define G_CONSTRAINT 0
#define DG_CONSTRAINT 1
#define LS_CONSTRAINT 2

#define LS_ABS  0
#define LS_REAL 1
#define LS_IMAG 2
#define LS_ALL  3

#define MAXIMUM_NUMBER_OF_EXTRA_CONSTRAINTS 20
/**
 * CLASS: Vasilyev_Filter_Constraints
 */
class Vasilyev_Filter_Constraints {
public:
    Vector3D k;
    int type;
    double target;
    int order_of_derivative;
    int LS_DOF;
    int LS_type;
};


template<typename Object_type, typename Member_Pointer,typename Solution_type>
class LeastSquares_coefficient_function_class{
public:
    /* constructor */
    LeastSquares_coefficient_function_class(Object_type *object,
                      Member_Pointer mem_func,
                      const int &l_, 
                      const int &m_, 
                      const int &direction_): Obj(object), Ptr(mem_func), l(l_), m(m_), direction(direction_) { }
        
    // "member function evaluation" with one parameter
    Solution_type operator() (const double &k){
        return (Obj->*Ptr)(l,m,k,direction);
    }
    
private:
    LeastSquares_coefficient_function_class();	/* make default constructor private */
    Object_type *Obj;		/*!< pointer to the object */
    Member_Pointer Ptr;		/*!< pointer to the class member function */
    const int l;
    const int m;
    const int direction;
};

template<typename Object_type, typename Member_Pointer,typename Solution_type>
class LeastSquares_RHS_function_class{
public:
    /* constructor */
    LeastSquares_RHS_function_class(Object_type *object,
                     Member_Pointer mem_func,
                     const int &m_,
                     const int &LS_DOF_,
                     const int direction_) : Obj(object), Ptr(mem_func), m(m_), LS_DOF(LS_DOF_), direction(direction_) {  }

    
    // "member function evaluation" with one parameter
    Solution_type operator() (const double &k){
        return (Obj->*Ptr)(m,k,LS_DOF,direction);
    }

private:
    LeastSquares_RHS_function_class();	/* make default constructor private */
    Object_type *Obj;		/*!< pointer to the object */
    Member_Pointer Ptr;		/*!< pointer to the class member function */
    const int m;
    const int LS_DOF;
    const int direction;
};

/**
 * CLASS: Vasilyev_Filter
 */
template <typename Soln_pState, typename Soln_cState>
class Vasilyev_Filter : public Discrete_Filter<Soln_pState,Soln_cState> {
public:
    using Discrete_Filter<Soln_pState,Soln_cState>::theCell;
    using Discrete_Filter<Soln_pState,Soln_cState>::theNeighbours;
    using Discrete_Filter<Soln_pState,Soln_cState>::number_of_rings;
    using Discrete_Filter<Soln_pState,Soln_cState>::commutation_order;
    using Discrete_Filter<Soln_pState,Soln_cState>::FGR;
    using Discrete_Filter<Soln_pState,Soln_cState>::target_filter_sharpness;
    using Discrete_Filter<Soln_pState,Soln_cState>::LS_constraints;
    using Discrete_Filter<Soln_pState,Soln_cState>::Derivative_constraints;
    using Discrete_Filter<Soln_pState,Soln_cState>::Filter_Width_strict;
    using Discrete_Filter<Soln_pState,Soln_cState>::Neighbouring_Values;
    using Discrete_Filter<Soln_pState,Soln_cState>::Set_Neighbouring_Values;
    using Discrete_Filter<Soln_pState,Soln_cState>::check_filter_moments;
    using Discrete_Filter<Soln_pState,Soln_cState>::check_filter_moments_1D;
    using Discrete_Filter<Soln_pState,Soln_cState>::G_cutoff;


    
    Vasilyev_Filter(void) : Discrete_Filter<Soln_pState,Soln_cState>() {
        number_of_constraints = 0;
        Output_Constraints = true;
        if (Explicit_Filter_Properties::batch_flag)
            Output_Constraints = false;
    }
    
    int number_of_constraints;
    
    Vasilyev_Filter_Constraints Constraints[MAXIMUM_NUMBER_OF_EXTRA_CONSTRAINTS];
    
    bool Output_Constraints;
    
    RowVector symmetric_weights;
    RowVector symmetric_weights_1D;

    void Get_Neighbours(Cell3D &theCell_) {
        theCell = theCell_;
        theNeighbours.GetNeighbours(theCell, number_of_rings,FILTER_TYPE_VASILYEV);
    }
    
    void Get_Neighbours_1D(Cell3D &theCell_, int direction) {
        theCell = theCell_;
        theNeighbours.GetNeighbours_1D(theCell, number_of_rings,FILTER_TYPE_VASILYEV,direction);
    }
    
    RowVector Get_Weights(Cell3D &theCell, Neighbours &theNeighbours);
    RowVector Get_Weights_1D(Cell3D &theCell, Neighbours &theNeighbours, int direction);
    RowVector Calculate_Weights_1D(Cell3D &theCell, Neighbours &theNeighbours, int direction);

    int number_of_combinations(void);
    int fac(int n);
    int number_of_combinations_with_moment(int k);
    
    
    int Set_basic_constraints(Neighbours &theNeighbours);
    void Add_extra_constraints(const int type, const int LS_DOF, const int LS_type);
    void Add_extra_constraints(const int type, const double target, const Vector3D &k);
    void Add_extra_constraints(const int type, const double target, const Vector3D &k, const int p);
    void Add_extra_constraints(const int type, const double target, const Vector3D &k, const int p, const int q, const int LS_type);

    double LeastSquares_coefficient_function(const int &l, const int &m, const double &k, const int &direction);
    double LeastSquares_coefficient_function_real(const int &l, const int &m, const double &k, const int &direction);
    double LeastSquares_coefficient_function_imag(const int &l, const int &m, const double &k, const int &direction);
    double LeastSquares_coefficient(const int &l, const int &m, const int &direction, const int &type);

    Complex Gx_func(const double &k, const int &direction, const DenseMatrix &A, const ColumnVector &z); 
    Complex G_target(const double &k, const int &LS_DOF, const int &direction);
    Complex G0_func(const double &k_1D, const int &direction, const Cell3D &theCell, const Neighbours &theNeigbhours_subset, const ColumnVector &w0);
    double LeastSquares_RHS_function(const int &m, const double &k, const int &LS_DOF, const int &direction);
    double LeastSquares_RHS_function_real(const int &m, const double &k, const int &LS_DOF, const int &direction);
    double LeastSquares_RHS_function_imag(const int &m, const double &k, const int &LS_DOF, const int &direction);
    double LeastSquares_RHS(const int &m, const int &LS_DOF, const int &direction, const int &type);
    double filter_moment(Cell3D &theCell, Neighbours &theNeighbours, RowVector &w, int q, int r, int s);
    double filter_moment_1D(Cell3D &theCell, Neighbours &theNeighbours, RowVector &w, int q, int direction);

    void filter_tests(Grid3D_Hexa_Block &Grid_Blk, Cell3D &theCell);

    int filter_type(void) { return FILTER_TYPE_VASILYEV; }
    string filter_name(void) { return "Vasilyev"; }

};


template<typename Soln_pState, typename Soln_cState>
void Vasilyev_Filter<Soln_pState,Soln_cState>::filter_tests(Grid3D_Hexa_Block &Grid_Blk, Cell3D &theCell) {
    
    theNeighbours.set_grid(Grid_Blk);
    Get_Neighbours_1D(theCell,X_DIRECTION);
    RowVector w = Get_Weights_1D(theCell, theNeighbours,X_DIRECTION);
    cout << "w = " << w << endl;
    
    
    
    check_filter_moments(Grid_Blk,theCell);
    check_filter_moments_1D(Grid_Blk,theCell,X_DIRECTION);

    
}


template <typename Soln_pState, typename Soln_cState>
inline RowVector Vasilyev_Filter<Soln_pState,Soln_cState>::Calculate_Weights_1D(Cell3D &theCell, Neighbours &theNeighbours, int direction) {
    
    
    /* ---------------- direction dependent configuration --------------- */
    int K,L,N;
    double (Vector3D::*k_member);
    switch(direction){
        case X_DIRECTION:
            K = theNeighbours.Ki;     L = theNeighbours.Li;     N = theNeighbours.Ni;
            k_member = &Vector3D::x;
            break;
        case Y_DIRECTION:
            K = theNeighbours.Kj;     L = theNeighbours.Lj;     N = theNeighbours.Nj;
            k_member = &Vector3D::y;
            break;
        case Z_DIRECTION:
            K = theNeighbours.Kk;     L = theNeighbours.Lk;     N = theNeighbours.Nk;
            k_member = &Vector3D::z;
            break;
    }
    
    
    /* ------------------------ moment constraints ------------------------ */
    
    DenseMatrix A(commutation_order,N);
    int index_l;
    for (int q=0; q<commutation_order; q++) {
        for (int l=-K; l<=L; l++) {
            index_l = K+l;
            A(q,index_l) = pow(double(l),double(q));
        }
    }
    ColumnVector b(commutation_order);
    b.zero();
    b(0) = ONE;
    
    /* ------------------------ extra constraints ------------------------- */
    int type;
    double target;
    double k;
    int p;
    int t;
    int LS_DOF;
    int LS_type;
    RowVector A_row(K+L+1);
    
    for (int i=0; i<number_of_constraints; i++) {
        
        type = Constraints[i].type;
        target = Constraints[i].target;
        k = Constraints[i].k.*k_member;
        p = Constraints[i].order_of_derivative;
        LS_DOF = Constraints[i].LS_DOF;
        LS_type = Constraints[i].LS_type;
        switch (type) {
            case G_CONSTRAINT:
                for (int l=-K; l<=L; l++) {
                    index_l = K+l;
                    A_row(index_l) = real( exp(-I*(k*double(l))) );
                }
                A.append(A_row);      b.append(target);
                break;
                
            case DG_CONSTRAINT:
                for (int l=-K; l<=L; l++) {
                    index_l = K+l;
                    if ( pow(-1.0,p) < ZERO )
                        A_row(index_l) = imag( pow(-I*(k*double(l)),p) * exp(-I*(k*double(l))));
                    else
                        A_row(index_l) = real( pow(-I*(k*double(l)),p) * exp(-I*(k*double(l))));
                }
                A.append(A_row);     b.append(target);
                break;
                
            case LS_CONSTRAINT:
                t = A.get_n()-1;       // the last row index
                t -= K;                // the next undefined coefficient
                for (int l=-K; l<=L; l++) {
                    index_l = K+l;
                    A_row(index_l) = LeastSquares_coefficient(l,t,direction,LS_type);
                }
                A.append(A_row);     b.append(LeastSquares_RHS(t,LS_DOF,direction,LS_type));
                
            break;
        }
    }
    
    /* -------------------- Solve the linear system ---------------------- */
    ColumnVector w(K+L+1);
    assert(A.get_n() == A.get_m());
    assert(A.get_n() == b.size());
    assert(b.size() == w.size());
    w = A.pseudo_inverse()*b;    
    
    return w;
}


template <typename Soln_pState, typename Soln_cState>
inline RowVector Vasilyev_Filter<Soln_pState,Soln_cState>::Get_Weights(Cell3D &theCell, Neighbours &theNeighbours) {
    // Set FGR to properties
    FGR = Explicit_Filter_Properties::FGR;    
    if (theNeighbours.symmetric_stencil && symmetric_weights.size()!=0) {
        return symmetric_weights;
    } else {
        int number_of_moment_combinations = number_of_combinations();
        int number_of_neighbours = theNeighbours.number_of_neighbours;
        
        complex<double> I(0.0,1.0);
        
        
        int Lr; // number of remaining constraints after setting basic constraints
        Lr = Set_basic_constraints(theNeighbours); 
        
        
        RowVector w_i = Calculate_Weights_1D(theCell,theNeighbours,X_DIRECTION);
        RowVector w_j = Calculate_Weights_1D(theCell,theNeighbours,Y_DIRECTION);
        RowVector w_k = Calculate_Weights_1D(theCell,theNeighbours,Z_DIRECTION);
        
        // Load coefficients into neighbour_weights
        RowVector W(number_of_neighbours);
        double denominator = 0;
        int index_l, index_m, index_n;
        for(int i=0; i<number_of_neighbours; i++) {
            index_l = theNeighbours.Ki - (theCell.I - theNeighbours.neighbour[i].I);
            index_m = theNeighbours.Kj - (theCell.J - theNeighbours.neighbour[i].J);
            index_n = theNeighbours.Kk - (theCell.K - theNeighbours.neighbour[i].K);
            W(i) = w_i(index_l)*w_j(index_m)*w_k(index_n) ;// / theNeighbours.neighbour[i].Jacobian;
            //denominator += W(i);
        }
        
        //W *= theCell.Jacobian;
        
        if (theNeighbours.symmetric_stencil) {
            symmetric_weights = W;
        }
        
        return W;
        
    }
    
}


template <typename Soln_pState, typename Soln_cState>
inline RowVector Vasilyev_Filter<Soln_pState,Soln_cState>::Get_Weights_1D(Cell3D &theCell, Neighbours &theNeighbours, int direction) {
    
    if (theNeighbours.symmetric_stencil && symmetric_weights_1D.size()!=0) {
        return symmetric_weights_1D;
    } else {
        
    int number_of_moment_combinations = number_of_combinations();
    int number_of_neighbours = theNeighbours.number_of_neighbours;
    
    complex<double> I(0.0,1.0);
    
    
    int Lr; // number of remaining constraints after setting basic constraints
    Lr = Set_basic_constraints(theNeighbours); 
    
    
    RowVector w = Calculate_Weights_1D(theCell,theNeighbours,direction);
    
#ifdef _GNUPLOT
    if (Explicit_Filter_Properties::debug_flag) {
        int n=200;
        double *Gr = new double [n];
        double *Gi = new double [n];
        double *Gt = new double [n];
        double *k2 = new double [n];
        double kmax = PI;
        for (int i=0; i<n; i++) {
            double k = i*kmax/(n-ONE);
            k2[i] = k/kmax;
            Gr[i] = real( G0_func(k, direction, theCell, theNeighbours, w) );
            Gi[i] = imag( G0_func(k, direction, theCell, theNeighbours, w) );
            Gt[i] = real( G_target(k, 0, direction) );
        }
        
        Gnuplot_Control h2;
        h2.gnuplot_init();
        h2.gnuplot_setstyle("lines") ;
        h2.gnuplot_cmd("set grid");
        h2.gnuplot_set_xlabel("k");
        h2.gnuplot_set_ylabel("G(k)");
        //h2.gnuplot_cmd("set yrange [-1:1]");
        h2.gnuplot_set_title("transfer function LS");
        h2.gnuplot_plot1d_var2(k2,Gr,n,"real");
        h2.gnuplot_plot1d_var2(k2,Gi,n,"imag");
        h2.gnuplot_plot1d_var2(k2,Gt,n,"target");
        
        delete[] Gr;
        delete[] Gi;
        delete[] Gt;
        delete[] k2;        
    }
#endif
    
    /* ----------------------- non-uniform grid -------------------------- *
    int J_index;
    switch (direction) {
        case X_DIRECTION:
            J_index = 0;
            break;
        case Y_DIRECTION:
            J_index = 1;
            break;
        case Z_DIRECTION:
            J_index = 2;
            break;
    }
    DenseMatrix Jacob(3,3);
    DenseMatrix Jacob_Cell(3,3);
    double denominator = 0;
    theNeighbours.Grid_ptr->Jacobian_Matrix(Jacob_Cell,theCell.I,theCell.J,theCell.K,commutation_order+1);
    for (int n=0; n<theNeighbours.number_of_neighbours; n++) {
        int I,J,K;
        I = theNeighbours.neighbour[n].I;
        J = theNeighbours.neighbour[n].J;
        K = theNeighbours.neighbour[n].K;
        
        theNeighbours.Grid_ptr->Jacobian_Matrix(Jacob,I,J,K,8);
        w(n) *= Jacob(J_index,J_index) ; // * (theNeighbours.neighbour[n].dXc.x) ;// /Jacob_Cell(J_index,J_index);
        denominator = Jacob_Cell(J_index,J_index);
    }
    w /= denominator; */

    if (theNeighbours.symmetric_stencil) {
        symmetric_weights_1D = w;
    }
    return w;
    }
};
template<typename Soln_pState, typename Soln_cState>
inline int Vasilyev_Filter<Soln_pState,Soln_cState>::Set_basic_constraints(Neighbours &theNeighbours) {


    number_of_constraints = 0;

    Vector3D kmax;
    kmax.x = PI;
    kmax.y = PI;
    kmax.z = PI;
    Vector3D k_FGR = kmax/FGR;   

    int type;
    double target;
    Vector3D k;
    
    int number_of_extra_constraints = 2*number_of_rings+1 - commutation_order;
    if (number_of_extra_constraints <= 0){
        if (Output_Constraints && CFFC_Primary_MPI_Processor()) {
            cout << "Can not set filter width with " << number_of_rings << " number of rings and commutation order " << commutation_order << " ." << endl;
        }
        while (number_of_extra_constraints <= 0) {
            commutation_order--;
            number_of_extra_constraints = 2*number_of_rings+1 - commutation_order;
        }
        if (Output_Constraints && CFFC_Primary_MPI_Processor()) {
            cout << "Commutation order lowered to " << commutation_order << " ." << endl;
        }
    }
    if (Output_Constraints && CFFC_Primary_MPI_Processor()) {
        cout << "Setting extra constraints on discrete filter:" << endl;
        cout << "   --> Number of rings = " << number_of_rings << endl;
        cout << "   --> Commutation order = " << commutation_order << endl;
    }
    int number_of_remaining_constraints = number_of_extra_constraints;    
     
     
    if (number_of_remaining_constraints > 0) {
        /* ------------- Grid cut off --------------- */
        type = G_CONSTRAINT;
        target = ZERO;
        
        k = kmax;
        Add_extra_constraints(type, target, k);
        if (Output_Constraints && CFFC_Primary_MPI_Processor()) {
            cout << "   --> Transfer function at grid cut off = " << target << endl;
        }
        number_of_remaining_constraints--;
    }
    
    
    if (number_of_remaining_constraints>0 && Filter_Width_strict) {
        /* ------------- Filter Grid Ratio --------------- */
        type = G_CONSTRAINT;
        target = G_cutoff;
        
        k = k_FGR;
        Add_extra_constraints(type, target, k);
        if (Output_Constraints && CFFC_Primary_MPI_Processor()) {
            cout << "   --> Filter grid ratio = " << FGR << endl;
        }
        number_of_remaining_constraints--;   
    }
    
    int p=1;
    int plast(0);
    if (Derivative_constraints == DEFAULT) {
        if (number_of_remaining_constraints > 2) {
            plast = 1;
        }
    } else {
        plast = Derivative_constraints;
    }

    while (number_of_remaining_constraints>0  && (p<=plast || !LS_constraints)) {
        /* -------------- Derivatives ----------------- */
        type = DG_CONSTRAINT;
        target = ZERO;
        k = kmax;
        Add_extra_constraints(type, target, k, p);
        if (Output_Constraints && CFFC_Primary_MPI_Processor()) {
            cout << "   --> Derivative " << p << " of transfer function at grid cut off = " << target << endl;
        }
        number_of_remaining_constraints--;
        p++;
    }
    
    
    int Lr = number_of_remaining_constraints;
    if (Output_Constraints && CFFC_Primary_MPI_Processor()) {
        cout << "   --> Least Squares constraints = " << Lr << endl;
    }
    
    int LS_type;
    while (number_of_remaining_constraints>0) {
        /* ------------ Least Squares Constraints ------------ */
        if (Lr == 1) {
            LS_type = LS_REAL;
        } else if (pow(double(number_of_remaining_constraints),-1.0) > ZERO) {
            LS_type = LS_REAL;
        } else {
            LS_type = LS_IMAG;
        }
        LS_type = LS_ALL;
        type = LS_CONSTRAINT;
        Add_extra_constraints(type,Lr,LS_type);
        number_of_remaining_constraints--;
    }
    
    Output_Constraints = false;
    
    if (number_of_remaining_constraints != 0){
        cerr << "Not enough constraints set in Vasilyev's filter" << endl; exit(0);
    }
    
    return number_of_remaining_constraints;
}
    

template<typename Soln_pState, typename Soln_cState>
inline double Vasilyev_Filter<Soln_pState,Soln_cState>::LeastSquares_coefficient_function(const int &l, const int &m, const double &k, const int &direction) {
    return (/**/real(exp(-I*k*double(l+m)))/**/+imag(exp(-I*k*double(l+m)))/**/);
}
template<typename Soln_pState, typename Soln_cState>
inline double Vasilyev_Filter<Soln_pState,Soln_cState>::LeastSquares_coefficient_function_real(const int &l, const int &m, const double &k, const int &direction) {
    return (/**/real(exp(-I*k*double(l+m)))/**+imag(exp(-I*k*double(l+m)))/**/);
}
template<typename Soln_pState, typename Soln_cState>
inline double Vasilyev_Filter<Soln_pState,Soln_cState>::LeastSquares_coefficient_function_imag(const int &l, const int &m, const double &k, const int &direction) {
    return (/**real(exp(-I*k*double(l+m)))/**/+imag(exp(-I*k*double(l+m)))/**/);
}


template<typename Soln_pState, typename Soln_cState>
inline double Vasilyev_Filter<Soln_pState,Soln_cState>::LeastSquares_coefficient(const int &l, const int &m, const int &direction, const int &type) {
    
    typedef double (Vasilyev_Filter<Soln_pState,Soln_cState>::*LeastSquares_coefficient_function_ptr_type) (const int &, const int &, const double &, const int &);
    LeastSquares_coefficient_function_class<Vasilyev_Filter<Soln_pState,Soln_cState>,LeastSquares_coefficient_function_ptr_type,double> C_LS (this, &Vasilyev_Filter<Soln_pState,Soln_cState>::LeastSquares_coefficient_function_real, l,  m,  direction);
    LeastSquares_coefficient_function_class<Vasilyev_Filter<Soln_pState,Soln_cState>,LeastSquares_coefficient_function_ptr_type,double> C_LS_real (this, &Vasilyev_Filter<Soln_pState,Soln_cState>::LeastSquares_coefficient_function_real, l,  m,  direction);
    LeastSquares_coefficient_function_class<Vasilyev_Filter<Soln_pState,Soln_cState>,LeastSquares_coefficient_function_ptr_type,double> C_LS_imag (this, &Vasilyev_Filter<Soln_pState,Soln_cState>::LeastSquares_coefficient_function_imag, l,  m,  direction);

    double dummy;
    switch (type) {
        case LS_ABS:
            return AdaptiveGaussianQuadrature(C_LS, ZERO, PI, dummy, numeric_limits<double>::digits10);
        case LS_REAL:
            return AdaptiveGaussianQuadrature(C_LS_real, ZERO, PI, dummy, numeric_limits<double>::digits10);
        case LS_IMAG:
            return AdaptiveGaussianQuadrature(C_LS_imag, ZERO, PI, dummy, numeric_limits<double>::digits10);
        case LS_ALL:
            return AdaptiveGaussianQuadrature(C_LS_real, ZERO, PI, dummy, numeric_limits<double>::digits10)
                 + AdaptiveGaussianQuadrature(C_LS_imag, ZERO, PI, dummy, numeric_limits<double>::digits10);
    }
}

template<typename Soln_pState, typename Soln_cState>
inline double Vasilyev_Filter<Soln_pState,Soln_cState>::LeastSquares_RHS(const int &m, const int &LS_DOF, const int &direction, const int &type) {
      
    typedef double (Vasilyev_Filter<Soln_pState,Soln_cState>::*LeastSquares_RHS_function_ptr_type) (const int &, const double &, const int &, const int &);
    LeastSquares_RHS_function_class<Vasilyev_Filter<Soln_pState,Soln_cState>,LeastSquares_RHS_function_ptr_type,double> R_LS (this, &Vasilyev_Filter<Soln_pState,Soln_cState>::LeastSquares_RHS_function_real,m,LS_DOF,direction);
    LeastSquares_RHS_function_class<Vasilyev_Filter<Soln_pState,Soln_cState>,LeastSquares_RHS_function_ptr_type,double> R_LS_real (this, &Vasilyev_Filter<Soln_pState,Soln_cState>::LeastSquares_RHS_function_real,m,LS_DOF,direction);
    LeastSquares_RHS_function_class<Vasilyev_Filter<Soln_pState,Soln_cState>,LeastSquares_RHS_function_ptr_type,double> R_LS_imag (this, &Vasilyev_Filter<Soln_pState,Soln_cState>::LeastSquares_RHS_function_imag,m,LS_DOF,direction);

    double dummy;
    switch (type) {
        case LS_ABS:
            return AdaptiveGaussianQuadrature(R_LS, ZERO, PI, dummy, numeric_limits<double>::digits10);
        case LS_REAL:
            return AdaptiveGaussianQuadrature(R_LS_real, ZERO, PI, dummy, numeric_limits<double>::digits10);
        case LS_IMAG:
            return AdaptiveGaussianQuadrature(R_LS_imag, ZERO, PI, dummy, numeric_limits<double>::digits10);
        case LS_ALL:
            return AdaptiveGaussianQuadrature(R_LS_real, ZERO, PI, dummy, numeric_limits<double>::digits10)
                 + AdaptiveGaussianQuadrature(R_LS_imag, ZERO, PI, dummy, numeric_limits<double>::digits10);
    }
}

template<typename Soln_pState, typename Soln_cState>
inline double Vasilyev_Filter<Soln_pState,Soln_cState>::LeastSquares_RHS_function(const int &m, const double &k, const int &LS_DOF, const int &direction) {
    return (/**/real(G_target(k,LS_DOF,direction) * exp(-I*(k*double(m))))/**/ + imag(G_target(k,LS_DOF,direction) * exp(-I*(k*double(m))))/**/);
}
template<typename Soln_pState, typename Soln_cState>
inline double Vasilyev_Filter<Soln_pState,Soln_cState>::LeastSquares_RHS_function_real(const int &m, const double &k, const int &LS_DOF, const int &direction) {
     return (/**/real(G_target(k,LS_DOF,direction) * exp(-I*(k*double(m))))/** + imag(G_target(k,LS_DOF,direction) * exp(-I*(k*double(m))))/**/);
}
template<typename Soln_pState, typename Soln_cState>
inline double Vasilyev_Filter<Soln_pState,Soln_cState>::LeastSquares_RHS_function_imag(const int &m, const double &k, const int &LS_DOF, const int &direction) {
    return (/**real(G_target(k,LS_DOF,direction) * exp(-I*(k*double(m))))/**/ + imag(G_target(k,LS_DOF,direction) * exp(-I*(k*double(m))))/**/);
}

template<typename Soln_pState, typename Soln_cState>
inline Complex Vasilyev_Filter<Soln_pState,Soln_cState>::G_target(const double &k, const int &LS_DOF, const int &direction) {
    double kmax = PI;

//    double s;
//    if (target_filter_sharpness >= 0) {
//        s = -target_filter_sharpness*FGR;
//    } else {
//        s = -FGR;
//    }
//    return Complex(HALF + HALF*tanh(s*(k-kmax/FGR)),ZERO);
//   
    int m;
    if (target_filter_sharpness > 0) {
        m = int(target_filter_sharpness);
    } else {
        m = int(ceil(commutation_order/2.0));
    }
    double a = -(TWO*m)*log(G_cutoff)*pow(PI/FGR,-TWO*m);
    return Complex(exp(-a/(TWO*m)*pow(k,TWO*m)),ZERO);
}




template<typename Soln_pState, typename Soln_cState>
inline Complex Vasilyev_Filter<Soln_pState,Soln_cState>::G0_func(const double &k, const int &direction, const Cell3D &theCell, const Neighbours &theNeigbhours_subset, const ColumnVector &w0) {
    
    Complex G(0,0);
    
    int l,m,n, index_l,index_m,index_n;
    int Ki,Li, Kj,Lj, Kk,Lk;
    Ki=theNeigbhours_subset.Ki;  Li=theNeigbhours_subset.Li;
    Kj=theNeigbhours_subset.Kj;  Lj=theNeigbhours_subset.Lj;
    Kk=theNeigbhours_subset.Kk;  Lk=theNeigbhours_subset.Lk;
    
    for (int i=0; i<theNeigbhours_subset.number_of_neighbours; i++) {
        
        l = theNeigbhours_subset.neighbour[i].I - theCell.I;
        m = theNeigbhours_subset.neighbour[i].J - theCell.J;
        n = theNeigbhours_subset.neighbour[i].K - theCell.K;
        
        
        switch (direction) {
            case X_DIRECTION:
                if (m==0 && n==0) {
                    index_l = Ki - (theCell.I - theNeigbhours_subset.neighbour[i].I);
                    G += w0(index_l)  *  ( exp(-I*k*double(l)) );
                }
                break;
            case Y_DIRECTION:
                if (l==0 && n==0) {
                    index_m = Kj - (theCell.J - theNeigbhours_subset.neighbour[i].J);
                    G += w0(index_m)  *  ( exp(-I*k*double(m)) );
                }
                break;
            case Z_DIRECTION:
                if (l==0 && m==0) {
                    index_n = Kk - (theCell.K - theNeigbhours_subset.neighbour[i].K);
                    G += w0(index_n)  *  ( exp(-I*k*double(n)) );
                }
                break;
        }        
    }  
    
    return G;             
}

template<typename Soln_pState, typename Soln_cState>
inline Complex Vasilyev_Filter<Soln_pState,Soln_cState>::Gx_func(const double &k, const int &direction, const DenseMatrix &A, const ColumnVector &z) {
    
    Complex G(0,0);

    int L,K;
    switch (direction) {
        case X_DIRECTION:
            L = theNeighbours.Li;   K = theNeighbours.Ki;
            break;
        case Y_DIRECTION:
            L = theNeighbours.Lj;   K = theNeighbours.Kj;
            break;
        case Z_DIRECTION:
            L = theNeighbours.Lk;   K = theNeighbours.Kk;
            break;           
    }    
    int Lr = z.size();
    for (int m=0; m<Lr; m++) {
        for (int l=-K; l<=L; l++) {
            G += A(l+K,m) * exp(-I*(k*double(l)));
        }
        G *= z(m);
    }
    return G;             
}


template<typename Soln_pState, typename Soln_cState>
inline void Vasilyev_Filter<Soln_pState,Soln_cState>::Add_extra_constraints(int type, int LS_DOF, int LS_type) {
    Add_extra_constraints(type, ZERO, Vector3D_ZERO, 0, LS_DOF, LS_type);
}

template<typename Soln_pState, typename Soln_cState>
inline void Vasilyev_Filter<Soln_pState,Soln_cState>::Add_extra_constraints(const int type, const double target, const Vector3D &k) {
    Add_extra_constraints(type, target, k, 0, 0, 0);
}

template<typename Soln_pState, typename Soln_cState>
inline void Vasilyev_Filter<Soln_pState,Soln_cState>::Add_extra_constraints(const int type, const double target, const Vector3D &k, const int p) {
    Add_extra_constraints(type, target, k, p, 0, 0);
}

#define NUMBER_OF_CONSTRAINTS_PARAMETERS    4
template<typename Soln_pState, typename Soln_cState>
inline void Vasilyev_Filter<Soln_pState,Soln_cState>::Add_extra_constraints(const int type, const double target, const Vector3D &k, const int p, const int q, const int LS_type) {
    assert(number_of_constraints<MAXIMUM_NUMBER_OF_EXTRA_CONSTRAINTS);
    Constraints[number_of_constraints].type=type;
    Constraints[number_of_constraints].target=target;
    Constraints[number_of_constraints].k=k;
    Constraints[number_of_constraints].order_of_derivative=p;
    Constraints[number_of_constraints].LS_DOF=q;    
    Constraints[number_of_constraints].LS_type=LS_type;    
    number_of_constraints++;
}


/* ------------ calculates the number of unknowns for the given order ----------- */
template <typename Soln_pState, typename Soln_cState>
inline int Vasilyev_Filter<Soln_pState,Soln_cState>::number_of_combinations(void) {
    int number = 0;
    for (int k=0; k< commutation_order; k++) {
        number += number_of_combinations_with_moment(k);
    }
    return number;
}

/* ----- faculty ----- */
template <typename Soln_pState, typename Soln_cState>
int Vasilyev_Filter<Soln_pState,Soln_cState>::fac(int n) {
    if (n<=1)
		return 1;
	else
		return (n*fac(n-1));
}

/* how many ways can one color k eggs with n colors?
 *  (n+k-1)!
 *  ---------    =  number of combinations with moment k with n dimensions
 *  k! (n-1)!
 */
template <typename Soln_pState, typename Soln_cState>
inline int Vasilyev_Filter<Soln_pState,Soln_cState>::number_of_combinations_with_moment(int k){
    int n=3; // 3 dimensions
    return (  fac(n+k-1)/(fac(k)*fac(n-1))  );
}

template <typename Soln_pState, typename Soln_cState>
double Vasilyev_Filter<Soln_pState,Soln_cState>::filter_moment(Cell3D &theCell, Neighbours &theNeighbours, RowVector &w, int q, int r, int s) {
    
    double M(0);
    for(int i=0; i<theNeighbours.number_of_neighbours; i++) {
        int l = theNeighbours.neighbour[i].I - theCell.I;
        int m = theNeighbours.neighbour[i].J - theCell.J;
        int n = theNeighbours.neighbour[i].K - theCell.K;
        M += w(i) * pow(double(l),double(q)) * pow(double(m),double(r)) * pow(double(n),double(s)) ;
    }
    return M;
}

template <typename Soln_pState, typename Soln_cState>
double Vasilyev_Filter<Soln_pState,Soln_cState>::filter_moment_1D(Cell3D &theCell, Neighbours &theNeighbours, RowVector &w, int q, int direction) {
    int K,L;
    switch(direction){
        case X_DIRECTION:
            K = theNeighbours.Ki;    L = theNeighbours.Li;
            break;
        case Y_DIRECTION:
            K = theNeighbours.Kj;    L = theNeighbours.Lj;
            break;
        case Z_DIRECTION:
            K = theNeighbours.Kk;    L = theNeighbours.Lk;
            break;
    }
    double M(0);

    for(int l=-K; l<=L; l++) {
        M += w(K+l) * pow(double(l),double(q)) ;
    }
    return M;
}

#endif

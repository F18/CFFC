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


template <typename Soln_pState, typename Soln_cState>
class LES_Filter;

#define G_CONSTRAINT 0
#define DG_CONSTRAINT 1
#define LS_CONSTRAINT 2

#define X_DIRECTION 1
#define Y_DIRECTION 2
#define Z_DIRECTION 3

#define LS_ABS  0
#define LS_REAL 1
#define LS_IMAG 2

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
class Vasilyev_Filter_LeastSquares_coefficient_function_class{
public:
    /* constructor */
    Vasilyev_Filter_LeastSquares_coefficient_function_class(Object_type *object,
                      Member_Pointer mem_func,
                      const double &delta_l_, 
                      const double &delta_m_, 
                      const int &direction_): Obj(object), Ptr(mem_func), delta_l(delta_l_), delta_m(delta_m_), direction(direction_) { }
        
    // "member function evaluation" with one parameter
    Solution_type operator() (const double &k){
        return (Obj->*Ptr)(delta_l,delta_m,k,direction);
    }
    
private:
    Vasilyev_Filter_LeastSquares_coefficient_function_class();	/* make default constructor private */
    Object_type *Obj;		/*!< pointer to the object */
    Member_Pointer Ptr;		/*!< pointer to the class member function */
    const double delta_l;
    const double delta_m;
    const int direction;
};

template<typename Object_type, typename Member_Pointer,typename Solution_type>
class Vasilyev_Filter_LeastSquares_RHS_function_class{
public:
    /* constructor */
    Vasilyev_Filter_LeastSquares_RHS_function_class(Object_type *object,
                     Member_Pointer mem_func,
                     const double &delta_m_,
                     const int &LS_DOF_,
                     const int direction_) : Obj(object), Ptr(mem_func), delta_m(delta_m_), LS_DOF(LS_DOF_), direction(direction_) {  }

    
    // "member function evaluation" with one parameter
    Solution_type operator() (const double &k){
        return (Obj->*Ptr)(delta_m,k,LS_DOF,direction);
    }
    
private:
    Vasilyev_Filter_LeastSquares_RHS_function_class();	/* make default constructor private */
    Object_type *Obj;		/*!< pointer to the object */
    Member_Pointer Ptr;		/*!< pointer to the class member function */
    const double delta_m;
    const int LS_DOF;
    const int direction;
};

/**
 * CLASS: Vasilyev_Filter
 */
template <typename Soln_pState, typename Soln_cState>
class Vasilyev_Filter : public Discrete_Filter<Soln_pState,Soln_cState> {
public:
    using Discrete_Filter<Soln_pState,Soln_cState>::theNeighbours;
    using Discrete_Filter<Soln_pState,Soln_cState>::number_of_rings;
    using Discrete_Filter<Soln_pState,Soln_cState>::commutation_order;
    using Discrete_Filter<Soln_pState,Soln_cState>::FGR;
    using Discrete_Filter<Soln_pState,Soln_cState>::target_filter_sharpness;
    using Discrete_Filter<Soln_pState,Soln_cState>::Neighbouring_Values;
    using Discrete_Filter<Soln_pState,Soln_cState>::Set_Neighbouring_Values;

    Gnuplot_Control LS_COEFF;
    bool with_coefficients;
    
    Vasilyev_Filter(void) : Discrete_Filter<Soln_pState,Soln_cState>() {
        number_of_constraints = 0;
        Output_Constraints = true;
        
        with_coefficients = false;
        
#ifdef _GNUPLOT
#undef _GNUPLOT
#endif
//#define _GNUPLOT
#ifdef _GNUPLOT
        LS_COEFF.gnuplot_init();
        LS_COEFF.gnuplot_setstyle("lines");
        if (with_coefficients)
            LS_COEFF.gnuplot_set_title("with_coefficients");
        else
            LS_COEFF.gnuplot_set_title("with_deltas");
#endif

    }
    
    int number_of_constraints;
    
    Vasilyev_Filter_Constraints Constraints[MAXIMUM_NUMBER_OF_EXTRA_CONSTRAINTS];
    
    bool Output_Constraints;

    void Get_Neighbours(Cell3D &theCell) {
        theNeighbours.GetNeighbours_Vasilyev(theCell, number_of_rings);
    }
    RowVector Get_Weights(Cell3D &theCell, Neighbours &theNeighbours);

    int number_of_combinations(void);
    int fac(int n);
    int number_of_combinations_with_moment(int k);
    
    
    int Set_basic_constraints(Neighbours &theNeighbours);
    void Add_extra_constraints(const int type, const int LS_DOF, const int LS_type);
    void Add_extra_constraints(const int type, const double target, const Vector3D &k);
    void Add_extra_constraints(const int type, const double target, const Vector3D &k, const int p);
    void Add_extra_constraints(const int type, const double target, const Vector3D &k, const int p, const int q, const int LS_type);


    double LeastSquares_coefficient_function(const double &delta_l, const double &delta_m, const double &k, const int &direction);
    double LeastSquares_coefficient_function_real(const double &delta_l, const double &delta_m, const double &k, const int &direction);
    double LeastSquares_coefficient_function_imag(const double &delta_l, const double &delta_m, const double &k, const int &direction);
    double LeastSquares_coefficient(const double &delta_l, const double &delta_m, const int &direction, const int &type);
    double LeastSquares_RHS_function(const double &delta_m, const double &k, const int &LS_DOF, const int &direction);
    double LeastSquares_RHS_function_real(const double &delta_m, const double &k, const int &LS_DOF, const int &direction);
    double LeastSquares_RHS_function_imag(const double &delta_m, const double &k, const int &LS_DOF, const int &direction);
    double LeastSquares_RHS(const double &delta_m, const int &LS_DOF, const int &direction, const int &type);
    
    Complex Gx_func(const double &k, const int &direction, const DenseMatrix &A, const ColumnVector &z); 
    Complex G_target(const double &k, const int &LS_DOF, const int &direction);
    Complex G0_func(const double &k_1D, const int &direction, const Cell3D &theCell, const Neighbours &theNeigbhours_subset, const ColumnVector &w0);
    
    void filter_tests(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell);

    int filter_type(void) { return FILTER_TYPE_VASILYEV; }
    string filter_name(void) { return "Vasilyev"; }

};


template<typename Soln_pState, typename Soln_cState>
void Vasilyev_Filter<Soln_pState,Soln_cState>::filter_tests(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell) {
    return;
}



template <typename Soln_pState, typename Soln_cState>
inline RowVector Vasilyev_Filter<Soln_pState,Soln_cState>::Get_Weights(Cell3D &theCell, Neighbours &theNeighbours) {
        
    int number_of_moment_combinations = number_of_combinations();
    int number_of_neighbours = theNeighbours.number_of_neighbours;

    complex<double> I(0.0,1.0);

    
    int Lr; // number of remaining constraints after setting basic constraints
    Lr = Set_basic_constraints(theNeighbours); 
    
    /* ---------------------- Moment constraints ------------------------ */
    int Ki,Li,Ni , Kj,Lj,Nj , Kk,Lk,Nk;
    Ki = theNeighbours.Ki;     Li = theNeighbours.Li;       Ni = theNeighbours.Ni;
    Kj = theNeighbours.Kj;     Lj = theNeighbours.Lj;       Nj = theNeighbours.Nj;
    Kk = theNeighbours.Kk;     Lk = theNeighbours.Lk;       Nk = theNeighbours.Nk;
    DenseMatrix Ax(commutation_order,Ni);
    DenseMatrix Ay(commutation_order,Nj);
    DenseMatrix Az(commutation_order,Nk);

    int index_l, index_m, index_n;
    for (int q=0; q<commutation_order; q++) {
        for (int i=0; i<Ni; i++) {
            int l = theNeighbours.neighbour_x[i].I - theCell.I;
            index_l = Ki+l;
            Ax(q,index_l) = pow(theNeighbours.neighbour_x[i].Xc.x-theCell.Xc.x,double(q));
        }
        
        for (int i=0; i<Nj; i++) {
            int m = theNeighbours.neighbour_y[i].J - theCell.J;
            index_m = Kj+m;
            Ay(q,index_m) = pow(theNeighbours.neighbour_y[i].Xc.y-theCell.Xc.y,double(q));
        }
        
        for (int i=0; i<Nk; i++) {
            int n = theNeighbours.neighbour_z[i].K - theCell.K;
            index_n = Kk+n;
            Az(q,index_n) = pow(theNeighbours.neighbour_z[i].Xc.z-theCell.Xc.z,double(q));                
        }
    }
    ColumnVector bx(commutation_order);
    bx.zero();
    bx(0) = ONE;
    ColumnVector by(commutation_order);
    by.zero();
    by(0) = ONE;
    ColumnVector bz(commutation_order);
    bz.zero();
    bz(0) = ONE;
    

    
    /* ------------------------ extra constraints ------------------------- */

    int type;
    double target;
    Vector3D k;
    int p;
    int tx,ty,tz;
    double delta_tx, delta_ty, delta_tz, delta_l, delta_m, delta_n;
    int LS_DOF;
    int LS_type;
    RowVector Ax_row(Ki+Li+1);
    RowVector Ay_row(Kj+Lj+1);
    RowVector Az_row(Kk+Lk+1);

    Vector3D Delta = theNeighbours.Delta;
    
    for (int i=0; i<number_of_constraints; i++) {
        
        type = Constraints[i].type;
        target = Constraints[i].target;
        k = Constraints[i].k;
        p = Constraints[i].order_of_derivative;
        LS_DOF = Constraints[i].LS_DOF;
        LS_type = Constraints[i].LS_type;
        switch (type) {
            case G_CONSTRAINT:
                for (int i=0; i<Ni; i++) {
                    int l = theNeighbours.neighbour_x[i].I - theCell.I;
                    index_l = Ki+l;
                    Ax_row(index_l) = real( exp(-I*(k.x*(theNeighbours.neighbour_x[i].Xc.x-theCell.Xc.x) ) ) );
                }
                
                for (int i=0; i<Nj; i++) {
                    int m = theNeighbours.neighbour_y[i].J - theCell.J;
                    index_m = Kj+m;
                    Ay_row(index_m) = real( exp(-I*(k.y*(theNeighbours.neighbour_y[i].Xc.y-theCell.Xc.y) ) ) );
                }
                
                for (int i=0; i<Nk; i++) {
                    int n = theNeighbours.neighbour_z[i].K - theCell.K;
                    index_n = Kk+n;
                    Az_row(index_n) = real( exp(-I*(k.z*(theNeighbours.neighbour_z[i].Xc.z-theCell.Xc.z) ) ) );
                }
                Ax.append(Ax_row);      bx.append(target);
                Ay.append(Ay_row);      by.append(target);
                Az.append(Az_row);      bz.append(target);
                
                break;
                
            case DG_CONSTRAINT:
                for (int i=0; i<Ni; i++) {
                    int l = theNeighbours.neighbour_x[i].I - theCell.I;
                    index_l = Ki+l;
                    if ( pow(-1.0,p) < ZERO )
                        Ax_row(index_l) = imag( pow(-I*(k.x*(theNeighbours.neighbour_x[i].Xc.x-theCell.Xc.x)),p) * exp(-I*(k.x*(theNeighbours.neighbour_x[i].Xc.x-theCell.Xc.x))));
                    else
                        Ax_row(index_l) = real( pow(-I*(k.x*(theNeighbours.neighbour_x[i].Xc.x-theCell.Xc.x)),p) * exp(-I*(k.x*(theNeighbours.neighbour_x[i].Xc.x-theCell.Xc.x))));
                }
                
                for (int i=0; i<Nj; i++) {
                    int m = theNeighbours.neighbour_y[i].J - theCell.J;
                    index_m = Kj+m;
                    if ( pow(-1.0,p) < ZERO )
                        Ay_row(index_m) = imag( pow(-I*(k.y*(theNeighbours.neighbour_y[i].Xc.y-theCell.Xc.y)),p) * exp(-I*(k.y*(theNeighbours.neighbour_y[i].Xc.y-theCell.Xc.y))));
                    else
                        Ay_row(index_m) = real( pow(-I*(k.y*(theNeighbours.neighbour_y[i].Xc.y-theCell.Xc.y)),p) * exp(-I*(k.y*(theNeighbours.neighbour_y[i].Xc.y-theCell.Xc.y))));
                }
                
                for (int i=0; i<Nk; i++) {
                    int n = theNeighbours.neighbour_z[i].K - theCell.K;
                    index_n = Kk+n;
                    if ( pow(-1.0,p) < ZERO )
                        Az_row(index_n) = imag( pow(-I*(k.z*(theNeighbours.neighbour_z[i].Xc.z-theCell.Xc.z)),p) * exp(-I*(k.z*(theNeighbours.neighbour_z[i].Xc.z-theCell.Xc.z))));
                    else
                        Az_row(index_n) = real( pow(-I*(k.z*(theNeighbours.neighbour_z[i].Xc.z-theCell.Xc.z)),p) * exp(-I*(k.z*(theNeighbours.neighbour_z[i].Xc.z-theCell.Xc.z))));
                }

                Ax.append(Ax_row);      bx.append(target);
                Ay.append(Ay_row);      by.append(target);
                Az.append(Az_row);      bz.append(target);
                break;
            
            case LS_CONSTRAINT:
                tx = Ax.get_n()-1;       // the last row index
                tx -= Ki;                // the next undefined coefficient
                ty = Ay.get_n()-1;       // the last row index
                ty -= Kj;                // the next undefined coefficient
                tz = Az.get_n()-1;       // the last row index
                tz -= Kk;                // the next undefined coefficient
                
                /* t is number of coefficient
                 * delta_t is the distance between theCell en neighbour for coefficient t
                 */
                for (int i=0; i<Ni; i++) {
                    int l = theNeighbours.neighbour_x[i].I - theCell.I;
                    index_l = Ki+l;
                    if(l == tx)
                        delta_tx = theNeighbours.neighbour_x[i].Xc.x-theCell.Xc.x;
                }
                for (int i=0; i<Nj; i++) {
                    int m = theNeighbours.neighbour_y[i].J - theCell.J;
                    index_m = Kj+m;
                    if(m == ty)
                        delta_ty = theNeighbours.neighbour_y[i].Xc.y-theCell.Xc.y;
                }
                for (int i=0; i<Nk; i++) {
                    int n = theNeighbours.neighbour_z[i].K - theCell.K;
                    index_n = Kk+n;
                    if(n == tz)
                        delta_tz = theNeighbours.neighbour_z[i].Xc.z-theCell.Xc.z;               
                }
                
                
                for (int i=0; i<Ni; i++) {
                    int l = theNeighbours.neighbour_x[i].I - theCell.I;
                    index_l = Ki+l;
                    delta_l = theNeighbours.neighbour[i].Xc.x-theCell.Xc.x;
                    if (with_coefficients)
                        Ax_row(index_l) = LeastSquares_coefficient(double(l),double(tx),X_DIRECTION,LS_type);
                    else
                        Ax_row(index_l) = LeastSquares_coefficient(delta_l,delta_tx,X_DIRECTION,LS_type);                    
                }
                for (int i=0; i<Nj; i++) {
                    int m = theNeighbours.neighbour_y[i].J - theCell.J;
                    index_m = Kj+m;
                    delta_m = theNeighbours.neighbour[i].Xc.y-theCell.Xc.y;
                    if (with_coefficients)
                        Ay_row(index_m) = LeastSquares_coefficient(double(m),double(ty),Y_DIRECTION,LS_type);
                    else
                        Ay_row(index_m) = LeastSquares_coefficient(delta_m,delta_ty,Y_DIRECTION,LS_type);
                }
                for (int i=0; i<Nk; i++) {
                    int n = theNeighbours.neighbour_z[i].K - theCell.K;
                    index_n = Kk+n;
                    delta_n = theNeighbours.neighbour[i].Xc.z-theCell.Xc.z;
                    if (with_coefficients)
                        Az_row(index_n) = LeastSquares_coefficient(double(n),double(tz),Z_DIRECTION,LS_type);
                    else
                        Az_row(index_n) = LeastSquares_coefficient(delta_n,delta_tz,Z_DIRECTION,LS_type);             
                }

                if (with_coefficients) {
                    Ax.append(Ax_row);      bx.append(LeastSquares_RHS(double(tx),LS_DOF,X_DIRECTION,LS_type));
                    Ay.append(Ay_row);      by.append(LeastSquares_RHS(double(ty),LS_DOF,Y_DIRECTION,LS_type));
                    Az.append(Az_row);      bz.append(LeastSquares_RHS(double(tz),LS_DOF,Z_DIRECTION,LS_type));
                } else {
                    Ax.append(Ax_row);      bx.append(LeastSquares_RHS(delta_tx,LS_DOF,X_DIRECTION,LS_type));
                    Ay.append(Ay_row);      by.append(LeastSquares_RHS(delta_ty,LS_DOF,Y_DIRECTION,LS_type));
                    Az.append(Az_row);      bz.append(LeastSquares_RHS(delta_tz,LS_DOF,Z_DIRECTION,LS_type));
                }
                break;
        }
    }
    
    /* -------------------- Solve 3 linear system ---------------------- *
     *               these are the basic filter weights                  *
     * ----------------------------------------------------------------- */
    ColumnVector w_i(Ki+Li+1);

    assert(Ax.get_n() == Ax.get_m());
    assert(Ax.get_n() == bx.size());
    assert(bx.size() == w_i.size());
    w_i = Ax.pseudo_inverse()*bx;
    
    ColumnVector w_j(Kj+Lj+1);
    assert(Ay.get_n() == Ay.get_m());
    assert(Ay.get_n() == by.size());
    assert(by.size() == w_j.size());
    w_j = Ay.pseudo_inverse()*by;
    
    ColumnVector w_k(Kk+Lk+1);
    assert(Az.get_n() == Az.get_m());
    assert(Az.get_n() == bz.size());
    assert(bz.size() == w_k.size());
    w_k = Az.pseudo_inverse()*bz;
    
    
#ifdef _GNUPLOT
#undef _GNUPLOT
#endif

//#define _GNUPLOT
#ifdef _GNUPLOT
    int n=200;
    double *Gr = new double [n];
    double *Gi = new double [n];
    double *Gt = new double [n];
    double *k2 = new double [n];
    double dx = theNeighbours.Delta.x;
    double kmax = PI/dx;
    for (int i=0; i<n; i++) {
        k2[i] = i*kmax/(n-ONE);
        Gr[i] = real( G0_func(k2[i], X_DIRECTION, theCell, theNeighbours, w_i) );
        Gi[i] = imag( G0_func(k2[i], X_DIRECTION, theCell, theNeighbours, w_i) );
        Gt[i] = real( G_target(k2[i], Lr, X_DIRECTION) );
    }
    
    Gnuplot_Control h2;
    h2.gnuplot_init();
    h2.gnuplot_setstyle("lines") ;
    h2.gnuplot_cmd("set grid");
    h2.gnuplot_set_xlabel("k");
    h2.gnuplot_set_ylabel("G(k)");
    //h2.gnuplot_cmd("set yrange [0:1]");
    h2.gnuplot_set_title("transfer function");
    h2.gnuplot_plot1d_var2(k2,Gr,n,"real");
    h2.gnuplot_plot1d_var2(k2,Gi,n,"imag");
    h2.gnuplot_plot1d_var2(k2,Gt,n,"target");
    
    delete[] Gr;
    delete[] Gi;
    delete[] Gt;
    delete[] k2;
#endif

    // Load coefficients into neighbour_weights
    RowVector W(number_of_neighbours);
    for(int i=0; i<number_of_neighbours; i++) {
        index_l = Ki - (theCell.I - theNeighbours.neighbour[i].I);
        index_m = Kj - (theCell.J - theNeighbours.neighbour[i].J);
        index_n = Kk - (theCell.K - theNeighbours.neighbour[i].K);
        W(i) = w_i(index_l)*w_j(index_m)*w_k(index_n);
    }
    
    return W;
}


template<typename Soln_pState, typename Soln_cState>
inline int Vasilyev_Filter<Soln_pState,Soln_cState>::Set_basic_constraints(Neighbours &theNeighbours) {


    number_of_constraints = 0;

    Vector3D Delta = theNeighbours.Delta;
    Vector3D kmax;
    kmax.x = PI/Delta.x;
    kmax.y = PI/Delta.y;
    kmax.z = PI/Delta.z;
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
    
    
    if (number_of_remaining_constraints <= 2 ) {
        /* ------------- Filter Grid Ratio --------------- */
        type = G_CONSTRAINT;
        target = HALF;
        
        k = k_FGR;
        Add_extra_constraints(type, target, k);
        if (Output_Constraints && CFFC_Primary_MPI_Processor()) {
            cout << "   --> Filter grid ratio = " << FGR << endl;
        }
        number_of_remaining_constraints--;   
    }
    
    int p=1;
    int plast = 0;
    if (number_of_remaining_constraints > 2) {
        plast = 1;
    }
    while (number_of_remaining_constraints>0  && p<=plast) {
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
        type = LS_CONSTRAINT;
        Add_extra_constraints(type,Lr,LS_type);
        number_of_remaining_constraints--;
    }
    
    Output_Constraints = false;
    
    return number_of_remaining_constraints;
}
    

template<typename Soln_pState, typename Soln_cState>
inline double Vasilyev_Filter<Soln_pState,Soln_cState>::LeastSquares_coefficient_function(const double &delta_l, const double &delta_m, const double &k, const int &direction) {
    double Delta;
    int L,K;
    switch (direction) {
        case X_DIRECTION:
            Delta = theNeighbours.Delta.x;
            break;
        case Y_DIRECTION:
            Delta = theNeighbours.Delta.y;
            break;
        case Z_DIRECTION:
            Delta = theNeighbours.Delta.z;
            break;           
    }
    
    if (with_coefficients)
        return (/**/real(exp(-I*((k*Delta)*(delta_l+delta_m)))) /**/ + imag(exp(-I*((k*Delta)*(delta_l+delta_m))))/**/);
    else
        return (/**/real(exp(-I*(k*(delta_l+delta_m)))) /**/         + imag(exp(-I*(k*(delta_l+delta_m))))/**/);
}
template<typename Soln_pState, typename Soln_cState>
inline double Vasilyev_Filter<Soln_pState,Soln_cState>::LeastSquares_coefficient_function_real(const double &delta_l, const double &delta_m, const double &k, const int &direction) {
    double Delta;
    int L,K;
    switch (direction) {
        case X_DIRECTION:
            Delta = theNeighbours.Delta.x;
            break;
        case Y_DIRECTION:
            Delta = theNeighbours.Delta.y;
            break;
        case Z_DIRECTION:
            Delta = theNeighbours.Delta.z;
            break;           
    }
    
    if (with_coefficients)
        return (/**/real(exp(-I*((k*Delta)*(delta_l+delta_m)))) /** + imag(exp(-I*((k*Delta)*(delta_l+delta_m))))/**/);
    else
        return (/**/real(exp(-I*(k*(delta_l+delta_m)))) /**         + imag(exp(-I*(k*(delta_l+delta_m))))/**/);
}
template<typename Soln_pState, typename Soln_cState>
inline double Vasilyev_Filter<Soln_pState,Soln_cState>::LeastSquares_coefficient_function_imag(const double &delta_l, const double &delta_m, const double &k, const int &direction) {
    double Delta;
    int L,K;
    switch (direction) {
        case X_DIRECTION:
            Delta = theNeighbours.Delta.x;
            break;
        case Y_DIRECTION:
            Delta = theNeighbours.Delta.y;
            break;
        case Z_DIRECTION:
            Delta = theNeighbours.Delta.z;
            break;           
    }
    
    if (with_coefficients)
        return (/**real(exp(-I*((k*Delta)*(delta_l+delta_m)))) /**/ + imag(exp(-I*((k*Delta)*(delta_l+delta_m))))/**/);
    else
        return (/**real(exp(-I*(k*(delta_l+delta_m)))) /**/         + imag(exp(-I*(k*(delta_l+delta_m))))/**/);
}

template<typename Soln_pState, typename Soln_cState>
inline double Vasilyev_Filter<Soln_pState,Soln_cState>::LeastSquares_coefficient(const double &delta_l, const double &delta_m, const int &direction, const int &type) {
    
    typedef double (Vasilyev_Filter<Soln_pState,Soln_cState>::*LeastSquares_coefficient_function_ptr_type) (const double &, const double &, const double &, const int &);
    Vasilyev_Filter_LeastSquares_coefficient_function_class<Vasilyev_Filter<Soln_pState,Soln_cState>,LeastSquares_coefficient_function_ptr_type,double> C_LS (this, &Vasilyev_Filter<Soln_pState,Soln_cState>::LeastSquares_coefficient_function, delta_l,  delta_m,  direction);
    Vasilyev_Filter_LeastSquares_coefficient_function_class<Vasilyev_Filter<Soln_pState,Soln_cState>,LeastSquares_coefficient_function_ptr_type,double> C_LS_real (this, &Vasilyev_Filter<Soln_pState,Soln_cState>::LeastSquares_coefficient_function_real, delta_l,  delta_m,  direction);
    Vasilyev_Filter_LeastSquares_coefficient_function_class<Vasilyev_Filter<Soln_pState,Soln_cState>,LeastSquares_coefficient_function_ptr_type,double> C_LS_imag (this, &Vasilyev_Filter<Soln_pState,Soln_cState>::LeastSquares_coefficient_function_imag, delta_l,  delta_m,  direction);
    
    double dummy;
    double kmax;
    switch (direction) {
        case X_DIRECTION:
            kmax = PI/theNeighbours.Delta.x;
#ifdef _GNUPLOT
            char title[50];
            sprintf(title,"LS_COEFF %lf",delta_l);
            LS_COEFF.gnuplot_plot1d_function(C_LS,ZERO,kmax,100,title);
#endif
            break;
        case Y_DIRECTION:
            kmax = PI/theNeighbours.Delta.y;  break;
        case Z_DIRECTION:
            kmax = PI/theNeighbours.Delta.z;  break;
    }
                                   
    switch (type) {
        case LS_ABS:
            return AdaptiveGaussianQuadrature(C_LS, ZERO, kmax, dummy, numeric_limits<double>::digits10);
        case LS_REAL:
            return AdaptiveGaussianQuadrature(C_LS_real, ZERO, kmax, dummy, numeric_limits<double>::digits10);
        case LS_IMAG:
            return AdaptiveGaussianQuadrature(C_LS_imag, ZERO, kmax, dummy, numeric_limits<double>::digits10);
    }
}

template<typename Soln_pState, typename Soln_cState>
inline double Vasilyev_Filter<Soln_pState,Soln_cState>::LeastSquares_RHS(const double &delta_m, const int &LS_DOF, const int &direction, const int &type) {
      
    typedef double (Vasilyev_Filter<Soln_pState,Soln_cState>::*LeastSquares_RHS_function_ptr_type) (const double &, const double &, const int &, const int &);
    Vasilyev_Filter_LeastSquares_RHS_function_class<Vasilyev_Filter<Soln_pState,Soln_cState>,LeastSquares_RHS_function_ptr_type,double> R_LS (this, &Vasilyev_Filter<Soln_pState,Soln_cState>::LeastSquares_RHS_function,delta_m,LS_DOF,direction);
    Vasilyev_Filter_LeastSquares_RHS_function_class<Vasilyev_Filter<Soln_pState,Soln_cState>,LeastSquares_RHS_function_ptr_type,double> R_LS_real (this, &Vasilyev_Filter<Soln_pState,Soln_cState>::LeastSquares_RHS_function_real,delta_m,LS_DOF,direction);
    Vasilyev_Filter_LeastSquares_RHS_function_class<Vasilyev_Filter<Soln_pState,Soln_cState>,LeastSquares_RHS_function_ptr_type,double> R_LS_imag (this, &Vasilyev_Filter<Soln_pState,Soln_cState>::LeastSquares_RHS_function_imag,delta_m,LS_DOF,direction);

    double dummy;
    double kmax;
    switch (direction) {
        case X_DIRECTION:
            kmax = PI/theNeighbours.Delta.x; 
#ifdef _GNUPLOT
            char title[50];
            sprintf(title,"RS_COEFF %lf",delta_m);
            Gnuplot_Control h1;
            h1.gnuplot_init();
            h1.gnuplot_setstyle("lines");
            if (with_coefficients)  
                h1.gnuplot_set_title("with coefficients");
            else
                h1.gnuplot_set_title("with deltas");
            h1.gnuplot_plot1d_function(R_LS,ZERO,kmax,100,title);
#endif
            break;
        case Y_DIRECTION:
            kmax = PI/theNeighbours.Delta.y;  break;
        case Z_DIRECTION:
            kmax = PI/theNeighbours.Delta.z;  break;
    }    
    switch (type) {
        case LS_ABS:
            return AdaptiveGaussianQuadrature(R_LS, ZERO, kmax, dummy, numeric_limits<double>::digits10);
        case LS_REAL:
            return AdaptiveGaussianQuadrature(R_LS_real, ZERO, kmax, dummy, numeric_limits<double>::digits10);
        case LS_IMAG:
            return AdaptiveGaussianQuadrature(R_LS_imag, ZERO, kmax, dummy, numeric_limits<double>::digits10);
    }
}

template<typename Soln_pState, typename Soln_cState>
inline double Vasilyev_Filter<Soln_pState,Soln_cState>::LeastSquares_RHS_function(const double &delta_m, const double &k, const int &LS_DOF, const int &direction) {
    double Delta;
    int L,K;
    switch (direction) {
        case X_DIRECTION:
            Delta = theNeighbours.Delta.x;
            break;
        case Y_DIRECTION:
            Delta = theNeighbours.Delta.y;
            break;
        case Z_DIRECTION:
            Delta = theNeighbours.Delta.z;
            break;           
    }
    if (with_coefficients)
        return (/**/real(G_target(k,LS_DOF,direction) * exp(-I*(k*Delta*delta_m)))/**/ + imag(G_target(k,LS_DOF,direction) * exp(-I*(k*Delta*delta_m)))/**/);
    else
        return (/**/real(G_target(k,LS_DOF,direction) * exp(-I*(k*delta_m))) /**/      + imag(G_target(k,LS_DOF,direction) * exp(-I*(k*delta_m))) /**/);
}
template<typename Soln_pState, typename Soln_cState>
inline double Vasilyev_Filter<Soln_pState,Soln_cState>::LeastSquares_RHS_function_real(const double &delta_m, const double &k, const int &LS_DOF, const int &direction) {
    double Delta;
    int L,K;
    switch (direction) {
        case X_DIRECTION:
            Delta = theNeighbours.Delta.x;
            break;
        case Y_DIRECTION:
            Delta = theNeighbours.Delta.y;
            break;
        case Z_DIRECTION:
            Delta = theNeighbours.Delta.z;
            break;           
    }
    if (with_coefficients)
        return (/**/real(G_target(k,LS_DOF,direction) * exp(-I*(k*Delta*delta_m)))/** + imag(G_target(k,LS_DOF,direction) * exp(-I*(k*Delta*delta_m)))/**/);
    else
        return (/**/real(G_target(k,LS_DOF,direction) * exp(-I*(k*delta_m))) /**      + imag(G_target(k,LS_DOF,direction) * exp(-I*(k*delta_m))) /**/);
}
template<typename Soln_pState, typename Soln_cState>
inline double Vasilyev_Filter<Soln_pState,Soln_cState>::LeastSquares_RHS_function_imag(const double &delta_m, const double &k, const int &LS_DOF, const int &direction) {
    double Delta;
    int L,K;
    switch (direction) {
        case X_DIRECTION:
            Delta = theNeighbours.Delta.x;
            break;
        case Y_DIRECTION:
            Delta = theNeighbours.Delta.y;
            break;
        case Z_DIRECTION:
            Delta = theNeighbours.Delta.z;
            break;           
    }
    if (with_coefficients)
        return (/**real(G_target(k,LS_DOF,direction) * exp(-I*(k*Delta*delta_m)))/**/ + imag(G_target(k,LS_DOF,direction) * exp(-I*(k*Delta*delta_m)))/**/);
    else
        return (/**real(G_target(k,LS_DOF,direction) * exp(-I*(k*delta_m))) /**/      + imag(G_target(k,LS_DOF,direction) * exp(-I*(k*delta_m))) /**/);
}

template<typename Soln_pState, typename Soln_cState>
inline Complex Vasilyev_Filter<Soln_pState,Soln_cState>::G_target(const double &k, const int &LS_DOF, const int &direction) {
    double kmax;
    switch (direction) {
        case X_DIRECTION:
            kmax = PI/theNeighbours.Delta.x;  break;
        case Y_DIRECTION:
            kmax = PI/theNeighbours.Delta.y;  break;
        case Z_DIRECTION:
            kmax = PI/theNeighbours.Delta.z;  break;
    }

    double s;
    if (target_filter_sharpness >= 0) {
        s = -target_filter_sharpness*FGR*PI/kmax;
    } else {
        s = -FGR * PI / kmax;
    }
    return Complex(HALF + HALF*tanh(s*(k-kmax/FGR)),ZERO);
}




template<typename Soln_pState, typename Soln_cState>
inline Complex Vasilyev_Filter<Soln_pState,Soln_cState>::G0_func(const double &k_1D, const int &direction, const Cell3D &theCell, const Neighbours &theNeigbhours_subset, const ColumnVector &w0) {
    
    Complex G(0,0);
    
    Vector3D X0(theCell.Xc);
    Vector3D X;
    Vector3D dX;
    Vector3D k;
    k.zero();
    switch (direction) {
        case X_DIRECTION:
            k.x = k_1D;     break;
        case Y_DIRECTION:
            k.y = k_1D;     break;
        case Z_DIRECTION:
            k.z = k_1D;     break;
    }
    
    int l,m,n, index_l,index_m,index_n;
    int Ki,Li, Kj,Lj, Kk,Lk;
    Ki=theNeigbhours_subset.Ki;  Li=theNeigbhours_subset.Li;
    Kj=theNeigbhours_subset.Kj;  Lj=theNeigbhours_subset.Lj;
    Kk=theNeigbhours_subset.Kk;  Lk=theNeigbhours_subset.Lk;
    
    for (int i=0; i<theNeigbhours_subset.number_of_neighbours; i++) {
        dX = theNeigbhours_subset.neighbour[i].Xc - X0;
        
        l = theNeigbhours_subset.neighbour[i].I - theCell.I;
        m = theNeigbhours_subset.neighbour[i].J - theCell.J;
        n = theNeigbhours_subset.neighbour[i].K - theCell.K;
        
        
        switch (direction) {
            case X_DIRECTION:
                if (m==0 && n==0) {
                    index_l = Ki - (theCell.I - theNeigbhours_subset.neighbour[i].I);
                    G += w0(index_l)  *  exp(-I*(k*dX));
                }
                break;
            case Y_DIRECTION:
                if (l==0 && n==0) {
                    index_m = Kj - (theCell.J - theNeigbhours_subset.neighbour[i].J);
                    G += w0(index_m)  *  exp(-I*(k*dX));
                }
                break;
            case Z_DIRECTION:
                if (l==0 && m==0) {
                    index_n = Kk - (theCell.K - theNeigbhours_subset.neighbour[i].K);
                    G += w0(index_n)  *  exp(-I*(k*dX));
                }
                break;
        }        
    }  
    
    return G;             
}

template<typename Soln_pState, typename Soln_cState>
inline Complex Vasilyev_Filter<Soln_pState,Soln_cState>::Gx_func(const double &k, const int &direction, const DenseMatrix &A, const ColumnVector &z) {
    
    Complex G(0,0);

    double Delta;
    int L,K;
    switch (direction) {
        case X_DIRECTION:
            L = theNeighbours.Li;   K = theNeighbours.Ki;
            Delta = theNeighbours.Delta.x;
            break;
        case Y_DIRECTION:
            L = theNeighbours.Lj;   K = theNeighbours.Kj;
            Delta = theNeighbours.Delta.y;
            break;
        case Z_DIRECTION:
            L = theNeighbours.Lk;   K = theNeighbours.Kk;
            Delta = theNeighbours.Delta.z;
            break;           
    }    
    int Lr = z.size();
    for (int m=0; m<Lr; m++) {
        for (int l=-K; l<=L; l++) {
            G += A(l+K,m) * exp(-I*(k*double(l)*Delta));
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

#endif
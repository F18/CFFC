/*
 *  Vasilyev_LS_Filter.h
 *  CFFC
 *
 *  Created by Willem Deconinck on 27/04/08.
 *
 */

#ifndef _VASILYEV_LS_FILTER_INCLUDED
#define _VASILYEV_LS_FILTER_INCLUDED


inline double stepfunction(const double &k) {
    double kmax;
    kmax = PI/(TWO*PI/32.0);
    
    double kFGR(kmax/2.0);
    if (k <= kFGR)
        return ONE;
    else 
        return ZERO;
}

#include "Discrete_Filter.h"


template <typename Soln_pState, typename Soln_cState>
class LES_Filter;

#define G_CONSTRAINT 0
#define DG_CONSTRAINT 1
#define LS_CONSTRAINT 2

#define X_DIRECTION 1
#define Y_DIRECTION 2
#define Z_DIRECTION 3

#define MAXIMUM_NUMBER_OF_EXTRA_CONSTRAINTS 20
/**
 * CLASS: Vasilyev_LS_Filter_Constraints
 */
class Vasilyev_LS_Filter_Constraints {
public:
    Vector3D k;
    int type;
    double target;
    int order_of_derivative;
    int LS_DOF;
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
 * CLASS: Vasilyev_LS_Filter
 */
template <typename Soln_pState, typename Soln_cState>
class Vasilyev_LS_Filter : public Discrete_Filter<Soln_pState,Soln_cState> {
public:
    using Discrete_Filter<Soln_pState,Soln_cState>::theCell;
    using Discrete_Filter<Soln_pState,Soln_cState>::theNeighbours;
    using Discrete_Filter<Soln_pState,Soln_cState>::number_of_rings;
    using Discrete_Filter<Soln_pState,Soln_cState>::commutation_order;
    using Discrete_Filter<Soln_pState,Soln_cState>::FGR;
    using Discrete_Filter<Soln_pState,Soln_cState>::target_filter_sharpness;
    using Discrete_Filter<Soln_pState,Soln_cState>::Neighbouring_Values;
    using Discrete_Filter<Soln_pState,Soln_cState>::Set_Neighbouring_Values;

    
    Vasilyev_LS_Filter(void) : Discrete_Filter<Soln_pState,Soln_cState>() {
        number_of_constraints = 0;
        Output_Constraints = true;
    }
    
    int number_of_constraints;
    
    Vasilyev_LS_Filter_Constraints Constraints[MAXIMUM_NUMBER_OF_EXTRA_CONSTRAINTS];
    
    bool Output_Constraints;

    void Get_Neighbours(Cell3D &theCell_) {
        theCell = theCell_;
        theNeighbours.GetNeighbours_Vasilyev(theCell, number_of_rings);
    }
    RowVector Get_Weights(Cell3D &theCell, Neighbours &theNeighbours);

    int number_of_combinations(void);
    int fac(int n);
    int number_of_combinations_with_moment(int k);
    
    
    int Set_basic_constraints(Neighbours &theNeighbours);
    void Add_extra_constraints(const int type, const int LS_DOF);
    void Add_extra_constraints(const int type, const double target, const Vector3D &k);
    void Add_extra_constraints(const int type, const double target, const Vector3D &k, const int p);
    void Add_extra_constraints(const int type, const double target, const Vector3D &k, const int p, const int q);


    double LeastSquares_coefficient_function(const int &l, const int &m, const double &k, const int &direction);
    double LeastSquares_coefficient(const int &l, const int &m, const int &direction);

    Complex Gx_func(const double &k, const int &direction, const DenseMatrix &A, const ColumnVector &z); 
    Complex G_target(const double &k, const int &LS_DOF, const int &direction);
    Complex G0_func(const double &k_1D, const int &direction, const Cell3D &theCell, const Neighbours &theNeigbhours_subset, const ColumnVector &w0);
    double LeastSquares_RHS_function(const int &m, const double &k, const int &LS_DOF, const int &direction);
    double LeastSquares_RHS(const int &m, const int &LS_DOF, const int &direction);

    
    void filter_tests(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell);

    int filter_type(void) { return FILTER_TYPE_VASILYEV_LS; }
    string filter_name(void) { return "Vasilyev Least Squares"; }

};


template<typename Soln_pState, typename Soln_cState>
void Vasilyev_LS_Filter<Soln_pState,Soln_cState>::filter_tests(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell) {
    return;
}



template <typename Soln_pState, typename Soln_cState>
inline RowVector Vasilyev_LS_Filter<Soln_pState,Soln_cState>::Get_Weights(Cell3D &theCell, Neighbours &theNeighbours) {
        
    int number_of_moment_combinations = number_of_combinations();
    int number_of_neighbours = theNeighbours.number_of_neighbours;

    complex<double> I(0.0,1.0);

    
    int Lr; // number of remaining constraints after setting basic constraints
    Lr = Set_basic_constraints(theNeighbours); 
    assert(Lr == 0);
    //assert(pow(-ONE,Lr) > 0);   // Lr must be even
    //cout << "Lr = " << Lr << endl;
    
    //Neighbours theNeighbours_subset(theNeighbours);
    //theNeighbours_subset.allocate(number_of_rings-Lr/2);
    //theNeighbours_subset.GetNeighbours_Vasilyev(theCell,number_of_rings-Lr/2);
    
    Vector3D X0(theCell.Xc);
    double x,y,z;
    
    /* ---------------------- X - direction filter ------------------------ */
    int Ki,Li;
    Ki = theNeighbours.Ki;     Li = theNeighbours.Li;
    
    DenseMatrix Ax(commutation_order,Ki+Li+1);
    int index_l;
    for (int q=0; q<commutation_order; q++) {
        for (int l=-Ki; l<=Li; l++) {
            index_l = Ki+l;
            x = theNeighbours.neighbour_x[index_l].Xc.x;
//            Ax(q,index_l) = pow(x-X0.x,q);
            Ax(q,index_l) = pow(double(l),double(q));
        }
    }
    ColumnVector bx(commutation_order);
    bx.zero();
    bx(0) = ONE;
    
    /* ---------------------- Y - direction filter ------------------------ */
    int Kj,Lj;
    Kj = theNeighbours.Kj;     Lj = theNeighbours.Lj;
    DenseMatrix Ay(commutation_order,Kj+Lj+1);
    int index_m;
    for (int q=0; q<commutation_order; q++) {
        for (int m=-Kj; m<=Lj; m++) {
            index_m = Kj+m;
            y = theNeighbours.neighbour_y[index_m].Xc.y;
//            Ay(q,index_m) = pow(y-X0.y,q);
            Ay(q,index_m) = pow(double(m),double(q));
        }
    }
    ColumnVector by(commutation_order);
    by.zero();
    by(0) = ONE;
    
    /* ---------------------- Z - direction filter ------------------------ */
    int Kk,Lk;
    Kk = theNeighbours.Kk;     Lk = theNeighbours.Lk;
    DenseMatrix Az(commutation_order,Kk+Lk+1);
    int index_n;
    for (int q=0; q<commutation_order; q++) {
        for (int n=-Kk; n<=Lk; n++) {
            index_n = Kk+n;
            z = theNeighbours.neighbour_z[index_n].Xc.z;
//            Az(q,index_n) = pow(z-X0.z,q);
            Az(q,index_n) = pow(double(n),double(q));
        }
    }
    ColumnVector bz(commutation_order);
    bz.zero();
    bz(0) = ONE;
    
    /* ------------------------ extra constraints ------------------------- */

    int type;
    double target;
    Vector3D k;
    int p;
    int t;
    int LS_DOF;
    RowVector Ax_row(Ki+Li+1);
    RowVector Ay_row(Kj+Lj+1);
    RowVector Az_row(Kk+Lk+1);

    Vector3D Delta = theNeighbours.Delta;
    double dx = Delta.x, dy = Delta.y, dz = Delta.z;
    
    for (int i=0; i<number_of_constraints; i++) {
        
        type = Constraints[i].type;
        target = Constraints[i].target;
        k = Constraints[i].k;
        p = Constraints[i].order_of_derivative;
        LS_DOF = Constraints[i].LS_DOF;
        Lr = LS_DOF;
        switch (type) {
            case G_CONSTRAINT:
                for (double l=-Ki; l<=Li; l++) {
                    index_l = Ki+l;
//                    x = theNeighbours.neighbour_x[index_l].Xc.x;
//                    Ax_row(index_l) = real( exp(-I*(k.x*(x-X0.x))) );
                    Ax_row(index_l) = real( exp(-I*(k.x*l*dx)) );
                }
                Ax.append(Ax_row);      bx.append(target);
                
                for (double m=-Kj; m<=Lj; m++) {
                    index_m = Kj+m;
//                    y = theNeighbours.neighbour_y[index_m].Xc.y;
//                    Ay_row(index_m) = real( exp(-I*(k.y*(y-X0.y))) );
                    Ay_row(index_m) = real( exp(-I*(k.y*m*dy)) );
                }
                Ay.append(Ay_row);      by.append(target);
                
                for (double n=-Kk; n<=Lk; n++) {
                    index_n = Kk+n;
//                    z = theNeighbours.neighbour_z[index_n].Xc.z;
//                    Az_row(index_n) = real( exp(-I*(k.z*(z-X0.z))) );
                    Az_row(index_n) = real( exp(-I*(k.z*n*dz)) );
                }
                Az.append(Az_row);      bz.append(target);
                
                break;
                
            case DG_CONSTRAINT:
                for (double l=-Ki; l<=Li; l++) {
                    index_l = Ki+l;
//                    x = theNeighbours.neighbour_x[index_l].Xc.x;
                    if ( pow(-1.0,p) < ZERO )
                        Ax_row(index_l) = imag( pow(-I*(k.x*l*dx),p) * exp(-I*(k.x*l*dx)));
//                        Ax_row(index_l) = imag( pow(-I*(k.x*(x-X0.x)),p) * exp(-I*(k.x*(x-X0.x))));
                    else
                        Ax_row(index_l) = real( pow(-I*(k.x*l*dx),p) * exp(-I*(k.x*l*dx)));
//                        Ax_row(index_l) = real( pow(-I*(k.x*(x-X0.x)),p) * exp(-I*(k.x*(x-X0.x))));
                }
                Ax.append(Ax_row);     bx.append(target);
                
                for (double m=-Kj; m<=Lj; m++) {
                    index_m = Kj+m;
//                    y = theNeighbours.neighbour_y[index_m].Xc.y;
                    if ( pow(-1.0,p) < ZERO )
                        Ay_row(index_m) = imag( pow(-I*(k.y*(y-X0.y)),p) * exp(-I*(k.y*m*dy)));
//                        Ay_row(index_m) = imag( pow(-I*(k.y*(y-X0.y)),p) * exp(-I*(k.y*(y-X0.y))));
                    else
                        Ay_row(index_m) = real( pow(-I*(k.y*(y-X0.y)),p) * exp(-I*(k.y*m*dy)));
//                        Ay_row(index_m) = real( pow(-I*(k.y*(y-X0.y)),p) * exp(-I*(k.y*(y-X0.y))));
                }
                Ay.append(Ay_row);     by.append(target);
                
                for (double n=-Kk; n<=Lk; n++) {
                    index_n = Kk+n;
//                    z = theNeighbours.neighbour_z[index_n].Xc.z;
                    if ( pow(-1.0,p) < ZERO )
                        Az_row(index_n) = imag( pow(-I*(k.z*(z-X0.z)),p) * exp(-I*(k.z*n*dz)));
//                        Az_row(index_n) = imag( pow(-I*(k.z*(z-X0.z)),p) * exp(-I*(k.z*(z-X0.z))));

                    else
                        Az_row(index_n) = real( pow(-I*(k.z*(z-X0.z)),p) * exp(-I*(k.z*n*dz)));
//                        Az_row(index_n) = real( pow(-I*(k.z*(z-X0.z)),p) * exp(-I*(k.z*(z-X0.z))));
                }
                Az.append(Az_row);     bz.append(target);
                
                break;
            
            case LS_CONSTRAINT:
                t = Ax.get_n()-1;       // the last row index
                t -= Ki;                // the next undefined coefficient
                for (double l=-Ki; l<=Li; l++) {
                    index_l = Ki+l;
                    Ax_row(index_l) = LeastSquares_coefficient(l,t,X_DIRECTION);
                }
                Ax.append(Ax_row);     bx.append(LeastSquares_RHS(t,LS_DOF,X_DIRECTION));
                
                t = Ay.get_n()-1;       // the last row index
                t -= Kj;                // the next undefined coefficient
                for (double m=-Kj; m<=Lj; m++) {
                    index_m = Kj+m;
                    Ay_row(index_m) = LeastSquares_coefficient(m,t,Y_DIRECTION);
                }
                Ay.append(Ay_row);     by.append(LeastSquares_RHS(t,LS_DOF,Y_DIRECTION));
                
                t = Az.get_n()-1;       // the last row index
                t -= Kk;                // the next undefined coefficient
                for (double n=-Kk; n<=Lk; n++) {
                    index_n = Kk+n;
                    Az_row(index_n) = LeastSquares_coefficient(n,t,Z_DIRECTION);
                }
                Az.append(Az_row);     bz.append(LeastSquares_RHS(t,LS_DOF,Z_DIRECTION));
                
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
    h2.gnuplot_cmd("set yrange [0:1]");
    h2.gnuplot_set_title("transfer function");
    h2.gnuplot_plot1d_var2(k2,Gr,n,"real");
    h2.gnuplot_plot1d_var2(k2,Gi,n,"imag");
    h2.gnuplot_plot1d_var2(k2,Gt,n,"target");
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
inline int Vasilyev_LS_Filter<Soln_pState,Soln_cState>::Set_basic_constraints(Neighbours &theNeighbours) {


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
    
    /* ------------- Filter Grid Ratio --------------- */
    type = G_CONSTRAINT;
    target = HALF;
    
    k = k_FGR;
    Add_extra_constraints(type, target, k);
    if (Output_Constraints && CFFC_Primary_MPI_Processor()) {
        cout << "   --> Filter grid ratio = " << FGR << endl;
    }
    number_of_remaining_constraints--;   
     
     
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
    
//    int p=1;
//    while (number_of_remaining_constraints>0 && p<=2) {
//        /* -------------- Derivatives ----------------- */
//        type = DG_CONSTRAINT;
//        target = ZERO;
//        k = kmax;
//        Add_extra_constraints(type, target, k, p);
//        if (Output_Constraints && CFFC_Primary_MPI_Processor()) {
//            cout << "   --> Derivative " << p << " of transfer function at grid cut off = " << target << endl;
//        }
//        number_of_remaining_constraints--;
//        p++;
//    }
    int Lr = number_of_remaining_constraints;
    cout << "   --> Least Squares constraints = " << Lr << endl;

    while (number_of_remaining_constraints>0) {
        /* ------------ Least Squares Constraints ------------ */
        type = LS_CONSTRAINT;
        Add_extra_constraints(type,Lr);
        number_of_remaining_constraints--;
    }
    
    return number_of_remaining_constraints;
}
    

template<typename Soln_pState, typename Soln_cState>
inline double Vasilyev_LS_Filter<Soln_pState,Soln_cState>::LeastSquares_coefficient_function(const int &l, const int &m, const double &k, const int &direction) {
    double Delta;
    //int K;
    switch (direction) {
        case X_DIRECTION:
            Delta = theNeighbours.Delta.x;
            //K = theNeighbours.Ki;
            //Delta = theNeighbours.neighbour_x[K+l+m].Xc.x - theCell.Xc.x;
            break;
        case Y_DIRECTION:
            Delta = theNeighbours.Delta.y;
            //K = theNeighbours.Kj;
            //Delta = theNeighbours.neighbour_y[K+l+m].Xc.y - theCell.Xc.y;
            break;
        case Z_DIRECTION:
            Delta = theNeighbours.Delta.z;
            //K = theNeighbours.Kk;
            //Delta = theNeighbours.neighbour_z[K+l+m].Xc.z - theCell.Xc.z;
            break;           
    }   
    //return (real(exp(-I*((k*Delta)*(double(l+m)))))+imag(exp(-I*((k*Delta)*(double(l+m))))));
    return (real(exp(-I*(k*Delta)))+imag(exp(-I*(k*Delta))));
}

template<typename Soln_pState, typename Soln_cState>
inline double Vasilyev_LS_Filter<Soln_pState,Soln_cState>::LeastSquares_coefficient(const int &l, const int &m, const int &direction) {
    
    typedef double (Vasilyev_LS_Filter<Soln_pState,Soln_cState>::*LeastSquares_coefficient_function_ptr_type) (const int &, const int &, const double &, const int &);
    LeastSquares_coefficient_function_class<Vasilyev_LS_Filter<Soln_pState,Soln_cState>,LeastSquares_coefficient_function_ptr_type,double> C_LS (this, &Vasilyev_LS_Filter<Soln_pState,Soln_cState>::LeastSquares_coefficient_function, l,  m,  direction);
    
    double dummy;
    double kmax;
    switch (direction) {
        case X_DIRECTION:
            kmax = PI/theNeighbours.Delta.x;  break;
        case Y_DIRECTION:
            kmax = PI/theNeighbours.Delta.y;  break;
        case Z_DIRECTION:
            kmax = PI/theNeighbours.Delta.z;  break;
    }
    return AdaptiveGaussianQuadrature(C_LS, ZERO, kmax, dummy, numeric_limits<double>::digits10);
}

template<typename Soln_pState, typename Soln_cState>
inline double Vasilyev_LS_Filter<Soln_pState,Soln_cState>::LeastSquares_RHS(const int &m, const int &LS_DOF, const int &direction) {
      
    typedef double (Vasilyev_LS_Filter<Soln_pState,Soln_cState>::*LeastSquares_RHS_function_ptr_type) (const int &, const double &, const int &, const int &);
    LeastSquares_RHS_function_class<Vasilyev_LS_Filter<Soln_pState,Soln_cState>,LeastSquares_RHS_function_ptr_type,double> R_LS (this, &Vasilyev_LS_Filter<Soln_pState,Soln_cState>::LeastSquares_RHS_function,m,LS_DOF,direction);
    
    double dummy;
    double kmax;
    switch (direction) {
        case X_DIRECTION:
            kmax = PI/theNeighbours.Delta.x;  break;
        case Y_DIRECTION:
            kmax = PI/theNeighbours.Delta.y;  break;
        case Z_DIRECTION:
            kmax = PI/theNeighbours.Delta.z;  break;
    }    
    return AdaptiveGaussianQuadrature(R_LS, ZERO, kmax, dummy, numeric_limits<double>::digits10);
}

template<typename Soln_pState, typename Soln_cState>
inline double Vasilyev_LS_Filter<Soln_pState,Soln_cState>::LeastSquares_RHS_function(const int &m, const double &k, const int &LS_DOF, const int &direction) {
    double Delta;
    //int K;
    switch (direction) {
        case X_DIRECTION:
            Delta = theNeighbours.Delta.x;
            //K = theNeighbours.Ki;
            //Delta = theNeighbours.neighbour_x[K+m].Xc.x - theCell.Xc.x;
            break;
        case Y_DIRECTION:
            Delta = theNeighbours.Delta.y;
            //K = theNeighbours.Kj;
            //Delta = theNeighbours.neighbour_y[K+m].Xc.y - theCell.Xc.y;
            break;
        case Z_DIRECTION:
            Delta = theNeighbours.Delta.z;
            //K = theNeighbours.Kk;
            //Delta = theNeighbours.neighbour_z[K+m].Xc.z - theCell.Xc.z;
            break;           
    }   
    //return (real(G_target(k,LS_DOF,direction) * exp(-I*(k*Delta))) + imag(G_target(k,LS_DOF,direction) * exp(-I*(k*Delta))));
    return (real(G_target(k,LS_DOF,direction) * exp(-I*(k*double(m)*Delta))) + imag(G_target(k,LS_DOF,direction) * exp(-I*(k*double(m)*Delta))));
}

template<typename Soln_pState, typename Soln_cState>
inline Complex Vasilyev_LS_Filter<Soln_pState,Soln_cState>::G_target(const double &k, const int &LS_DOF, const int &direction) {
    double kmax;
    switch (direction) {
        case X_DIRECTION:
            kmax = PI/theNeighbours.Delta.x;  break;
        case Y_DIRECTION:
            kmax = PI/theNeighbours.Delta.y;  break;
        case Z_DIRECTION:
            kmax = PI/theNeighbours.Delta.z;  break;
    }
    double d = fabs((kmax/FOUR)/(LS_DOF)  - (kmax/TWENTY)*(commutation_order));
    if (target_filter_sharpness >= 0) {
        d = (kmax * target_filter_sharpness )/TWO;
    }
    double kFGR(kmax/FGR);
    d = min(min(d,kFGR),kmax-kFGR);
    if (k <= (kFGR-d))
        return Complex(ONE,ZERO);
    else if (k >= (kFGR+d))
        return Complex(ZERO,ZERO);
    else 
        return (Complex(ONE - (k-(kFGR-d))/(TWO*d),ZERO));
}




template<typename Soln_pState, typename Soln_cState>
inline Complex Vasilyev_LS_Filter<Soln_pState,Soln_cState>::G0_func(const double &k_1D, const int &direction, const Cell3D &theCell, const Neighbours &theNeigbhours_subset, const ColumnVector &w0) {
    
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
                    G += w0(index_l)  *  real( exp(-I*(k*dX)) );
                }
                break;
            case Y_DIRECTION:
                if (l==0 && n==0) {
                    index_m = Kj - (theCell.J - theNeigbhours_subset.neighbour[i].J);
                    G += w0(index_m)  *  real( exp(-I*(k*dX)) );
                }
                break;
            case Z_DIRECTION:
                if (l==0 && m==0) {
                    index_n = Kk - (theCell.K - theNeigbhours_subset.neighbour[i].K);
                    G += w0(index_n)  *  real( exp(-I*(k*dX)) );
                }
                break;
        }        
    }  
    
    return G;             
}

template<typename Soln_pState, typename Soln_cState>
inline Complex Vasilyev_LS_Filter<Soln_pState,Soln_cState>::Gx_func(const double &k, const int &direction, const DenseMatrix &A, const ColumnVector &z) {
    
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
inline void Vasilyev_LS_Filter<Soln_pState,Soln_cState>::Add_extra_constraints(int type, int LS_DOF) {
    Add_extra_constraints(type, ZERO, Vector3D_ZERO, 0, LS_DOF);
}

template<typename Soln_pState, typename Soln_cState>
inline void Vasilyev_LS_Filter<Soln_pState,Soln_cState>::Add_extra_constraints(const int type, const double target, const Vector3D &k) {
    Add_extra_constraints(type, target, k, 0, 0);
}

template<typename Soln_pState, typename Soln_cState>
inline void Vasilyev_LS_Filter<Soln_pState,Soln_cState>::Add_extra_constraints(const int type, const double target, const Vector3D &k, const int p) {
    Add_extra_constraints(type, target, k, p, 0);
}

#define NUMBER_OF_CONSTRAINTS_PARAMETERS    4
template<typename Soln_pState, typename Soln_cState>
inline void Vasilyev_LS_Filter<Soln_pState,Soln_cState>::Add_extra_constraints(const int type, const double target, const Vector3D &k, const int p, const int q) {
    assert(number_of_constraints<MAXIMUM_NUMBER_OF_EXTRA_CONSTRAINTS);
    Constraints[number_of_constraints].type=type;
    Constraints[number_of_constraints].target=target;
    Constraints[number_of_constraints].k=k;
    Constraints[number_of_constraints].order_of_derivative=p;
    Constraints[number_of_constraints].LS_DOF=q;    
    number_of_constraints++;
}


/* ------------ calculates the number of unknowns for the given order ----------- */
template <typename Soln_pState, typename Soln_cState>
inline int Vasilyev_LS_Filter<Soln_pState,Soln_cState>::number_of_combinations(void) {
    int number = 0;
    for (int k=0; k< commutation_order; k++) {
        number += number_of_combinations_with_moment(k);
    }
    return number;
}

/* ----- faculty ----- */
template <typename Soln_pState, typename Soln_cState>
int Vasilyev_LS_Filter<Soln_pState,Soln_cState>::fac(int n) {
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
inline int Vasilyev_LS_Filter<Soln_pState,Soln_cState>::number_of_combinations_with_moment(int k){
    int n=3; // 3 dimensions
    return (  fac(n+k-1)/(fac(k)*fac(n-1))  );
}

#endif
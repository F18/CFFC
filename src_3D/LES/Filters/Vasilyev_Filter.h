/*
 *  Vasilyev_Filter.h
 *  CFFC
 *
 *  Created by Willem Deconinck on 27/04/08.
 *
 */

#ifndef _VASILYEV_FIlTER_INCLUDED
#define _VASILYEV_FILTER_INCLUDED


#include "General_Filter.h"

#include <complex>
#include <cmath> 

#include "../../Math/Math.h"
#include "../../Math/LinearSystems.h"
#include "../../Math/Matrix.h"

#ifndef _GNUPLOT
#define _GNUPLOT
#include "../../System/gnuplot.h"
#endif

#include "../../Utilities/Utilities.h"

#include <cstdlib>

#include "Discrete_Filter.h"

template <typename Soln_pState, typename Soln_cState>
class LES_Filter;

#define G_CONSTRAINT 0
#define DG_CONSTRAINT 1

#define MAXIMUM_NUMBER_OF_EXTRA_CONSTRAINTS 10
/**
 * CLASS: Vasilyev_Filter_Constraints
 */
class Vasilyev_Filter_Constraints {
public:
    Vector3D k;
    int type;
    double target;
    int order_of_derivative;
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
    using Discrete_Filter<Soln_pState,Soln_cState>::Neighbouring_Values;
    using Discrete_Filter<Soln_pState,Soln_cState>::Set_Neighbouring_Values;

    
    Vasilyev_Filter(void) : Discrete_Filter<Soln_pState,Soln_cState>() {
        number_of_constraints = 0;
        Output_Constraints = true;
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
    
    
    void Set_extra_constraints(Neighbours &theNeighbours);
    void Add_extra_constraints(int type, double target, Vector3D &k);
    void Add_extra_constraints(int type, double target, Vector3D &k, int p);
    
    
    void filter_tests(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell);

    
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

    /* ---------------------- X - direction filter ------------------------ */
    int Ki,Li;
    Ki = theNeighbours.Ki;     Li = theNeighbours.Li;
    DenseMatrix Ax(commutation_order,Ki+Li+1);
    int index_l;
    for (int q=0; q<commutation_order; q++) {
        for (int l=-Ki; l<=Li; l++) {
            index_l = Ki+l;
            Ax(q,index_l) = pow(double(l),q);
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
            Ay(q,index_m) = pow(double(m),q);
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
            Az(q,index_n) = pow(double(n),q);
        }
    }
    ColumnVector bz(commutation_order);
    bz.zero();
    bz(0) = ONE;
    
    /* ------------------------ extra constraints ------------------------- */
    Set_extra_constraints(theNeighbours);

    int type;
    double target;
    Vector3D k;
    int p;
    RowVector Ax_row(Ki+Li+1);
    RowVector Ay_row(Kj+Lj+1);
    RowVector Az_row(Kk+Lk+1);

    Vector3D Delta = theNeighbours.Delta;
    
    for (int i=0; i<number_of_constraints; i++) {
        
        type = Constraints[i].type;
        target = Constraints[i].target;
        k = Constraints[i].k;
        p = Constraints[i].order_of_derivative;
        
        switch (type) {
            case G_CONSTRAINT:
                for (double l=-Ki; l<=Li; l++) {
                    index_l = Ki+l;
                    Ax_row(index_l) = real( exp(-I*(k.x*l*Delta.x)) );
                }
                Ax.append(Ax_row);      bx.append(target);
                
                for (double m=-Kj; m<=Lj; m++) {
                    index_m = Kj+m;
                    Ay_row(index_m) = real( exp(-I*(k.y*m*Delta.y)) );
                }
                Ay.append(Ay_row);      by.append(target);
                
                for (double n=-Kk; n<=Lk; n++) {
                    index_n = Kk+n;
                    Az_row(index_n) = real( exp(-I*(k.z*n*Delta.z)) );
                }
                Az.append(Az_row);      bz.append(target);
                
                break;
                
            case DG_CONSTRAINT:
                for (double l=-Ki; l<=Li; l++) {
                    index_l = Ki+l;
                    if ( pow(-1.0,p) < ZERO )
                        Ax_row(index_l) = imag( pow(-I*(k.x*l*Delta.x),p) * exp(-I*(k.x*l*Delta.x)));
                    else
                        Ax_row(index_l) = real( pow(-I*(k.x*l*Delta.x),p) * exp(-I*(k.x*l*Delta.x)));
                }
                Ax.append(Ax_row);     bx.append(target);
                
                for (double m=-Kj; m<=Lj; m++) {
                    index_m = Kj+m;
                    if ( pow(-1.0,p) < ZERO )
                        Ay_row(index_m) = imag( pow(-I*(k.y*m*Delta.y),p) * exp(-I*(k.y*m*Delta.y)));
                    else
                        Ay_row(index_m) = real( pow(-I*(k.y*m*Delta.y),p) * exp(-I*(k.y*m*Delta.y)));
                }
                Ay.append(Ay_row);     by.append(target);
                
                for (double n=-Kk; n<=Lk; n++) {
                    index_n = Kk+n;
                    if ( pow(-1.0,p) < ZERO )
                        Az_row(index_n) = imag( pow(-I*(k.z*n*Delta.z),p) * exp(-I*(k.z*n*Delta.z)));
                    else
                        Az_row(index_n) = real( pow(-I*(k.z*n*Delta.z),p) * exp(-I*(k.z*n*Delta.z)));
                }
                Az.append(Az_row);     bz.append(target);
                
                break;
        }
    }
    
    /* -------------------- Solve 3 linear system ---------------------- */
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
inline void Vasilyev_Filter<Soln_pState,Soln_cState>::Set_extra_constraints(Neighbours &theNeighbours) {


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


    
    /* ------------- Filter Grid Ratio --------------- */
    type = G_CONSTRAINT;
    target = HALF;
    
    k = k_FGR;
    Add_extra_constraints(type, target, k);
    if (Output_Constraints && CFFC_Primary_MPI_Processor()) {
        cout << "   --> Filter grid ratio = " << FGR << endl;
    }
    
    int number_of_remaining_constraints = number_of_extra_constraints-1;    
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


    int p=1;
    while (number_of_remaining_constraints > 0) {
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
    Output_Constraints = false;
}

template<typename Soln_pState, typename Soln_cState>
inline void Vasilyev_Filter<Soln_pState,Soln_cState>::Add_extra_constraints(int type, double target, Vector3D &k) {
    Add_extra_constraints(type, target, k,  0);
}

#define NUMBER_OF_CONSTRAINTS_PARAMETERS    4
template<typename Soln_pState, typename Soln_cState>
inline void Vasilyev_Filter<Soln_pState,Soln_cState>::Add_extra_constraints(int type, double target, Vector3D &k, int p) {
    assert(number_of_constraints<MAXIMUM_NUMBER_OF_EXTRA_CONSTRAINTS);
    Constraints[number_of_constraints].type=type;
    Constraints[number_of_constraints].target=target;
    Constraints[number_of_constraints].k=k;
    Constraints[number_of_constraints].order_of_derivative=p;    
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
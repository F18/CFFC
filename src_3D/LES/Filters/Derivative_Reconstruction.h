/*
 *  Derivative_Reconstruction.h
 *  CFFC
 *
 *  Created by Willem Deconinck on 09/06/08.
 *
 */

#ifndef _DERIVATIVE_RECONSTRUCTION_INCLUDED
#define _DERIVATIVE_RECONSTRUCTION_INCLUDED


#include "../../Math/Math.h"
#include "../../Math/LinearSystems.h"
#include "../../Math/Matrix.h"
#include "../../Utilities/Utilities.h"


template <typename Soln_pState, typename Soln_cState>
class LES_Filter;

/**
 * CLASS: Derivative_Reconstruction
 */
template <typename Soln_pState, typename Soln_cState>
class Derivative_Reconstruction {
public:
    
    Derivative_Reconstruction(void) {
        order_of_accuracy = LES_Filter<Soln_pState,Soln_cState>::commutation_order+1;
        number_of_rings = LES_Filter<Soln_pState,Soln_cState>::number_of_rings;
        Store_Derivative_Reconstruction_Weights = !(LES_Filter<Soln_pState,Soln_cState>::Memory_Efficient);
        theNeighbours.allocate(number_of_rings);
        the_number_of_unknowns = number_of_unknowns();
    }
    Derivative_Reconstruction(int order_of_accuracy_, int number_of_rings_) {
        order_of_accuracy = order_of_accuracy_;
        number_of_rings = number_of_rings_;
        Store_Derivative_Reconstruction_Weights = !(LES_Filter<Soln_pState,Soln_cState>::Memory_Efficient);
        theNeighbours.allocate(number_of_rings);
        the_number_of_unknowns = number_of_unknowns();
    }
    
    Neighbours theNeighbours;
    int order_of_accuracy;
    int number_of_rings;
    bool Store_Derivative_Reconstruction_Weights;
        
    DenseMatrix Delta_Neighbouring_Values;
    void Set_Delta_Neighbouring_Values(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell, Neighbours &theNeighbours, DenseMatrix &Delta_Neighbouring_Values);
    
    DenseMatrix gradient(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell);
    RowVector divergence(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell);
    
    void Allocate_Derivative_Reconstruction_Weights(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk);

    void Get_Neighbours(Cell3D &theCell) ;
    DenseMatrix Get_Weights(Cell3D &theCell, Neighbours &theNeighbours);
    DenseMatrix Matrix_A(Cell3D &theCell, Neighbours &theNeighbours);
    DenseMatrix Matrix_b(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Neighbours &theNeighbours);
    
    /* --------- functions for Matrix A calculations --------- */ 
    int number_of_unknowns();
    int the_number_of_unknowns;
    int number_of_terms_of_degree(int order);
    int fac(int n);
    double trinomial_coefficient(int n1, int n2, int n3);
    
    
};


template<typename Soln_pState, typename Soln_cState>
inline void Derivative_Reconstruction<Soln_pState,Soln_cState>::Allocate_Derivative_Reconstruction_Weights(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk) {
    SolnBlk.Derivative_Reconstruction_Weights = new DenseMatrix **[SolnBlk.NCi];
    SolnBlk.Derivative_Reconstruction_Weights_Assigned = new bool **[SolnBlk.NCi];
    for (int i=0; i<SolnBlk.NCi; i++) {
        SolnBlk.Derivative_Reconstruction_Weights[i] = new DenseMatrix *[SolnBlk.NCj];
        SolnBlk.Derivative_Reconstruction_Weights_Assigned[i] = new bool *[SolnBlk.NCj];
        for (int j=0; j<SolnBlk.NCj; j++) {
            SolnBlk.Derivative_Reconstruction_Weights[i][j] = new DenseMatrix [SolnBlk.NCk];
            SolnBlk.Derivative_Reconstruction_Weights_Assigned[i][j] = new bool [SolnBlk.NCk];
            
            for (int k=0; k<SolnBlk.NCk; k++) {
                SolnBlk.Derivative_Reconstruction_Weights_Assigned[i][j][k] = false;
            }
        }
    }
}

template<typename Soln_pState, typename Soln_cState>
inline void Derivative_Reconstruction<Soln_pState,Soln_cState>::Get_Neighbours(Cell3D &theCell) {
    theNeighbours.GetNeighbours(theCell, number_of_rings);
}

template <typename Soln_pState, typename Soln_cState>
void Derivative_Reconstruction<Soln_pState,Soln_cState>::Set_Delta_Neighbouring_Values(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell, Neighbours &theNeighbours, DenseMatrix &Delta_Neighbouring_Values) {
    RowVector This_Cell_Values;
    LES_Filter<Soln_pState,Soln_cState>::what_to_filter(SolnBlk,theCell.I,theCell.J,theCell.K,This_Cell_Values);
    int the_number_of_neighbours = theNeighbours.number_of_neighbours;
    int the_number_of_variables = This_Cell_Values.size();
    
    if (Delta_Neighbouring_Values.get_n()!=the_number_of_neighbours || Delta_Neighbouring_Values.get_m()!=the_number_of_variables) {
        Delta_Neighbouring_Values.newsize(the_number_of_neighbours,the_number_of_variables);
    }
    
    RowVector One_Neighbour_Values;
    
    int I,J,K;
    for (int i=0; i<the_number_of_neighbours; i++) {
        
        
        
        I = theNeighbours.neighbour[i].I;
        J = theNeighbours.neighbour[i].J;
        K = theNeighbours.neighbour[i].K;
        
        LES_Filter<Soln_pState,Soln_cState>::what_to_filter(SolnBlk,I,J,K,One_Neighbour_Values);
        RowVector Delta(One_Neighbour_Values-This_Cell_Values);
        Delta_Neighbouring_Values.assignRow(i,Delta);
    }    
}



template<typename Soln_pState, typename Soln_cState>
inline DenseMatrix Derivative_Reconstruction<Soln_pState,Soln_cState>::gradient(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell) {
    
    
    if (!SolnBlk.Derivative_Reconstruction_Weights_Allocated && Store_Derivative_Reconstruction_Weights) {
        Allocate_Derivative_Reconstruction_Weights(SolnBlk);
        SolnBlk.Derivative_Reconstruction_Weights_Allocated = true;
    }
    
    theNeighbours.set_grid(SolnBlk.Grid);
    Get_Neighbours(theCell);
    
    int I(theCell.I);
    int J(theCell.J);
    int K(theCell.K);
    if (Store_Derivative_Reconstruction_Weights) {
        if (!SolnBlk.Derivative_Reconstruction_Weights_Assigned[I][J][K]) {
            SolnBlk.Derivative_Reconstruction_Weights[I][J][K] = Get_Weights(theCell,theNeighbours);
            SolnBlk.Derivative_Reconstruction_Weights_Assigned[I][J][K] = true;
            Set_Delta_Neighbouring_Values(SolnBlk,theCell,theNeighbours,Delta_Neighbouring_Values);
            return SolnBlk.Derivative_Reconstruction_Weights[I][J][K]*Delta_Neighbouring_Values;
        }
    } else {
        DenseMatrix W = Get_Weights(theCell,theNeighbours);
        Set_Delta_Neighbouring_Values(SolnBlk,theCell,theNeighbours,Delta_Neighbouring_Values);
        return Get_Weights(theCell,theNeighbours) * Delta_Neighbouring_Values;
    }
 }

template<typename Soln_pState, typename Soln_cState>
inline RowVector Derivative_Reconstruction<Soln_pState,Soln_cState>::divergence(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell) {
    
    DenseMatrix grad = gradient(SolnBlk,theCell);
    
    int n = grad.get_m();
    RowVector diver(n);
    for (int i=0; i<n; i++) {
        diver(i) = grad(0,i) + grad(1,i) + grad(2,i);
    }
    return diver ;
}



/* ------------ calculates the number of unknowns for the given order ----------- */
template <typename Soln_pState, typename Soln_cState>
inline int Derivative_Reconstruction<Soln_pState,Soln_cState>::number_of_unknowns(void) {
    int number = 0;
    for (int k=0; k<= order_of_accuracy; k++) {
        number += number_of_terms_of_degree(k);
    }
    return number;
}

/* ----- faculty ----- */
template <typename Soln_pState, typename Soln_cState>
int Derivative_Reconstruction<Soln_pState,Soln_cState>::fac(int n) {
    if (n<=1)
		return 1;
	else
		return (n*fac(n-1));
}

/* how many ways can one color k eggs with n colors?
 *  (n+k-1)!
 *  ---------    =  number of terms of degree k with n dimensions
 *  k! (n-1)!
 */
template <typename Soln_pState, typename Soln_cState>
inline int Derivative_Reconstruction<Soln_pState,Soln_cState>::number_of_terms_of_degree(int k){
    int n=3; // 3 dimensions
    return (  fac(n+k-1)/(fac(k)*fac(n-1))  );
}


template<typename Soln_pState, typename Soln_cState>
inline DenseMatrix Derivative_Reconstruction<Soln_pState,Soln_cState>::Get_Weights(Cell3D &theCell, Neighbours &theNeighbours) {

    DenseMatrix A = Matrix_A(theCell,theNeighbours);
    
    //DiagonalMatrix W = Matrix_W(theCell, theNeighbours, weight_factor);
    //DenseMatrix Z = (W*A).pseudo_inverse()*W;
    
    DenseMatrix Z = A.pseudo_inverse();
    
    DenseMatrix weights(3,theNeighbours.number_of_neighbours);
    assert(weights.get_m() == Z[0].size());
    assert(weights.get_m() == Z[1].size());
    assert(weights.get_m() == Z[2].size());

    weights.assignRow(0,Z[0]);
    weights.assignRow(1,Z[1]);
    weights.assignRow(2,Z[2]);
    
    return weights;
}


template<typename Soln_pState, typename Soln_cState>
inline DenseMatrix Derivative_Reconstruction<Soln_pState,Soln_cState>::Matrix_b(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Neighbours &theNeighbours) {
    
    
    RowVector RowVec;
    LES_Filter<Soln_pState,Soln_cState>::what_to_filter(SolnBlk,0,0,0,RowVec);
    int the_number_of_neighbours = theNeighbours.number_of_neighbours;
    
    DenseMatrix b(the_number_of_neighbours,RowVec.size());
    
    int I,J,K;
    for (int i=0; i<the_number_of_neighbours; i++) {
        I = theNeighbours.neighbour[i].I;
        J = theNeighbours.neighbour[i].J;
        K = theNeighbours.neighbour[i].K;
        LES_Filter<Soln_pState,Soln_cState>::what_to_filter(SolnBlk,I,J,K,b,i);
    }
    return b;
}

template<typename Soln_pState, typename Soln_cState>
inline DenseMatrix Derivative_Reconstruction<Soln_pState,Soln_cState>::Matrix_A(Cell3D &theCell, Neighbours &theNeighbours) {
    int the_number_of_neighbours = theNeighbours.number_of_neighbours;
    DenseMatrix A(the_number_of_neighbours,the_number_of_unknowns);
    int j=0;
    double dx, dy, dz, coefficient;
    for (int k=1; k<=order_of_accuracy; k++) {  // from k=1 and not k=0 because the cell value is not an unknown!!!
        
        // for all terms of degree k
        for (int n3=0; n3<=k; n3++) {
            for (int n2=0; n2<=k; n2++) {
                for (int n1=0; n1<=k; n1++) {
                    if (n1+n2+n3 == k) {
                        
                        // fill this column
                        coefficient = trinomial_coefficient(n1,n2,n3)/double(fac(k));
                        for (int i=0; i<the_number_of_neighbours; i++) {
                            dx = theNeighbours.neighbour[i].Xc.x - theCell.Xc.x;
                            dy = theNeighbours.neighbour[i].Xc.y - theCell.Xc.y;
                            dz = theNeighbours.neighbour[i].Xc.z - theCell.Xc.z;
                            A(i,j) = coefficient * pow(dx,n1)*pow(dy,n2)*pow(dz,n3);
                            
                        }
                        j++; // go to next column
                        
                    }
                }
            }
        } 
        
    }
    return A;    
}


/* ------------------------ trinomial coefficient -------------------------- */
/* 
 * coefficient of x^n1 y^n2 z^n3 in expansion of (x+y+z)^n (with n=n1+n2+n3)
 *  (n1 + n2 + n3)!             (               (n1 + n2 + n3 + ... + nk)!   )
 *  ---------------             ( multinomial:  ---------------------------  )
 *   n1!  n2!  n3!              (                 n1!  n2!  n3!  ...  nk!    )
 */
template<typename Soln_pState, typename Soln_cState>
inline double Derivative_Reconstruction<Soln_pState,Soln_cState>::trinomial_coefficient(int n1, int n2, int n3){
    return (  double(fac(n1+n2+n3))/double( fac(n1) * fac(n2) * fac(n3) ) );
}


#endif

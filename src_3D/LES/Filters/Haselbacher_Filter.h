/*
 *  Haselbacher_Filter.h
 *  CFFC
 *
 *  Created by Willem Deconinck on 25/04/08.
 *
 */

#ifndef _HASELBACHER_FIlTER_INCLUDED
#define _HASELBACHER_FILTER_INCLUDED


#include "Discrete_Filter.h"

template <typename Soln_pState, typename Soln_cState>
class LES_Filter;


/**
 * CLASS: Haselbacher_Filter
 */
template <typename Soln_pState, typename Soln_cState>
class Haselbacher_Filter : public Discrete_Filter<Soln_pState,Soln_cState> {
public:
    using Discrete_Filter<Soln_pState,Soln_cState>::theNeighbours;
    using Discrete_Filter<Soln_pState,Soln_cState>::number_of_rings;
    using Discrete_Filter<Soln_pState,Soln_cState>::commutation_order;
    using Discrete_Filter<Soln_pState,Soln_cState>::FGR;
    using Discrete_Filter<Soln_pState,Soln_cState>::Neighbouring_Values;
    using Discrete_Filter<Soln_pState,Soln_cState>::Set_Neighbouring_Values;
    using Discrete_Filter<Soln_pState,Soln_cState>::G_function;
    using Discrete_Filter<Soln_pState,Soln_cState>::dG_function;
    using Discrete_Filter<Soln_pState,Soln_cState>::G_function_1D;
    using Discrete_Filter<Soln_pState,Soln_cState>::G_function_1D_100;
    using Discrete_Filter<Soln_pState,Soln_cState>::G_function_1D_110;
    using Discrete_Filter<Soln_pState,Soln_cState>::Calculate_wavenumber_of_Gvalue;
    using Discrete_Filter<Soln_pState,Soln_cState>::Calculate_wavenumber_of_dGvalue;
    using Discrete_Filter<Soln_pState,Soln_cState>::Filter_Grid_Ratio;
    


    
    //! Constructor
    Haselbacher_Filter(void) : Discrete_Filter<Soln_pState,Soln_cState>() {
        Assigned_w_Uniform_Grid = false;
        UNIFORM_GRID = true;
        weight_factor = 1.5;
        the_number_of_unknowns = number_of_unknowns();
        relaxation_factor = 0.08;
        relaxation_applied = false;
    }
    //! Destructor
    ~Haselbacher_Filter(void) {
        deallocate_W_row();
    }

    void Get_Neighbours(Cell3D &theCell) {
        theNeighbours.GetNeighbours(theCell, number_of_rings);
        theNeighbours.append_theCell(theCell);
    }
    
private:
    
    double relaxation_factor;
    int relaxation_flag;
    int weight_flag;
    double weight_factor;
    bool relaxation_applied;
    
    /* ---------- non uniform grid -- save rows ------------- */
    RowVector ***W_row;
    bool *** Assigned_W_row;
    bool Allocated_W_row;
    void allocate_W_row(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk);
    void deallocate_W_row(void);
    void Assign_W_row(Cell3D &theCell, Neighbours &theNeighbours);
    int NCi, NCj, NCk;
    
    /* ------------------- uniform grid ---------------------- */
    RowVector w_Uniform_Grid;
    bool Assigned_w_Uniform_Grid;
    bool UNIFORM_GRID;
    void Assign_w_Uniform_Grid(Cell3D &theCell, Neighbours &theNeighbours);
    
    /* --------- functions for Matrix A calculations --------- */ 
    int number_of_unknowns();
    int number_of_unknowns(int &commutation_order);
    int the_number_of_unknowns;
    int number_of_terms_of_degree(int order);
    int fac(int n);
    double trinomial_coefficient(int n1, int n2, int n3);
    
    
    /* -------------- Least Squares Reconstruction ----------- */
    RowVector LeastSquaresReconstruction(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell, Neighbours &theNeighbours);
    RowVector LeastSquaresReconstruction_Uniform_Grid(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell, Neighbours &theNeighbours);
    RowVector LeastSquaresReconstruction_W_row(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell, Neighbours &theNeighbours);
    DenseMatrix Matrix_A(Cell3D &theCell, Neighbours &theNeighbours);
    DenseMatrix Matrix_A(Cell3D &theCell, Neighbours &theNeighbours, int &commutation_order);
    DenseMatrix Matrix_b(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Neighbours &theNeighbours);
    DiagonalMatrix Matrix_W(Cell3D &theCell, Neighbours &theNeighbours);
    DiagonalMatrix Matrix_W(Cell3D &theCell, Neighbours &theNeighbours, double weight_factor);
    void Weight_Matrix_A_and_b(DenseMatrix &A, DenseMatrix &b, Cell3D &theCell, Neighbours &theNeighbours);

    RowVector Get_Weights(Cell3D &theCell, Neighbours &theNeighbours);
    void Apply_relaxation(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &kmax, RowVector &w);
    void Undo_relaxation(Neighbours &theNeighbours, RowVector &w);
    
    /* --------------- Transfer function calculations ------------- */

    
  //  void Optimize_Filter(Cell3D &theCell, Neighbours &theNeighbours, double &kmax, double &FGR, int &commutation_order);
    double Calculate_weight_factor(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &kmax, double &FGR, int &commutation_order);
    int Calculate_weight_factor(void);
    double Filter_Grid_Ratio(Cell3D &theCell, Neighbours &theNeighbours, double &weight, DenseMatrix &A, Vector3D &kmax);
    double Filter_Grid_Ratio(int number_of_rings, double commutation_order, double weight);
    double Calculate_relaxation_factor(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &kmax, RowVector &w);
    double filter_quality(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell, double &kmax, int number_of_rings, int commutation_order, double weight_factor);
    double filter_quality(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &kmax, RowVector &w);
    double filter_sharpness(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &kmax, RowVector &w);
    double filter_smoothness(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &kmax, RowVector &w, double max_smoothness);
    double filter_uniformity(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &kmax, RowVector &w, double max_distortion);

    int Output_Filter_types(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell);
    int Output_Filter_types(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell,  int number_of_rings);
    int Output_Filter_types(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell, int number_of_rings, int commutation_order);


    void filter_tests(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell);

    int filter_type(void) { return FILTER_TYPE_HASELBACHER; }
    string filter_name(void) { return "Haselbacher"; }
    
    
};


template<typename Soln_pState, typename Soln_cState>
inline RowVector Haselbacher_Filter<Soln_pState,Soln_cState>::Get_Weights(Cell3D &theCell, Neighbours &theNeighbours) {
    theNeighbours.delete_theCell();
    //cout << "theCell = (" << theCell.I << "," << theCell.J << "," << theCell.K << ")" << endl;
    int error_flag;
    DenseMatrix A = Matrix_A(theCell,theNeighbours);
    error_flag = Calculate_weight_factor();
    if (error_flag) {
        cout << " error in calculating weight for Cell ("<<theCell.I<<","<<theCell.J<<","<<theCell.K<<")" << endl;
        exit(2);
    }

    DiagonalMatrix W = Matrix_W(theCell, theNeighbours, weight_factor);
    DenseMatrix Z = (W*A).pseudo_inverse()*W;
    RowVector w = Z[0];
    
    Vector3D kmax;
    kmax.x = PI/theNeighbours.Delta.x;
    kmax.y = PI/theNeighbours.Delta.y;
    kmax.z = PI/theNeighbours.Delta.z;
    Apply_relaxation(theCell, theNeighbours, kmax, w);
    return w;
}

template<typename Soln_pState, typename Soln_cState>
inline void Haselbacher_Filter<Soln_pState,Soln_cState>::Apply_relaxation(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &kmax, RowVector &w){
    relaxation_factor = Calculate_relaxation_factor(theCell,theNeighbours,kmax,w);
    w *= (ONE-relaxation_factor);
    w.append(relaxation_factor);
    theNeighbours.append_theCell(theCell);
}

template<typename Soln_pState, typename Soln_cState>
inline void Haselbacher_Filter<Soln_pState,Soln_cState>::Undo_relaxation(Neighbours &theNeighbours, RowVector &w){
    theNeighbours.delete_theCell();
//    w /= (ONE-relaxation_factor);
//    int m = w.size();
//    RowVector temp(m-1);
//    for (int i=0; i<m-1; i++) {
//        temp(i) = w(i);
//    }
//    w = temp;
}

//template<typename Soln_pState, typename Soln_cState>
//inline RowVector Haselbacher_Filter<Soln_pState,Soln_cState>::filter(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell) {
//    if (true) {
//        if (!SolnBlk.Filter_Weights_Allocated) {
//            Allocate_Filter_Weights(SolnBlk);
//            SolnBlk.Filter_Weights_Allocated = true;
//        }
//        
//        theNeighbours.set_grid(SolnBlk.Grid);
//        theNeighbours.GetNeighbours(theCell,number_of_rings);
//        
//        int I,J,K;
//        I = theCell.I;  J = theCell.J;  K = theCell.K;
//        if (!SolnBlk.Filter_Weights_Assigned[I][J][K]) {
//            SolnBlk.Filter_Weights[I][J][K] = Get_Weights(theCell,theNeighbours);
//            SolnBlk.Filter_Weights_Assigned[I][J][K] = true;
//        }
//        Set_Neighbouring_Values(SolnBlk,theNeighbours,Neighbouring_Values);
//        
//        return SolnBlk.Filter_Weights[I][J][K]*Neighbouring_Values;
//    }
//
//    else {
//        Neighbours theNeighbours(SolnBlk.Grid);
//        theNeighbours.GetNeighbours(theCell,number_of_rings);
//        
//        if (true) {
//            RowVector W_row = Get_Weights(theCell,theNeighbours);
//            DenseMatrix b = Matrix_b(SolnBlk,theNeighbours);
//            return (W_row*b);
//        } 
//        else if (UNIFORM_GRID) {   // SUPERFAST
//            if (!Assigned_w_Uniform_Grid) {
//                Assign_w_Uniform_Grid(theCell,theNeighbours);  
//            }
//            if(w_Uniform_Grid.size() != theNeighbours.number_of_neighbours) {
//                cout << "size(w) = "<< w_Uniform_Grid.size() << "    number_of_neighbours = " << theNeighbours.number_of_neighbours << endl;
//                Assign_w_Uniform_Grid(theCell,theNeighbours);
//            }
//            return LeastSquaresReconstruction_Uniform_Grid(SolnBlk, theCell, theNeighbours);
//        }
//        else if (true) {   // VERY SLOW INITIALIZATION BUT AFTER SUPERFAST
//            if (!Allocated_W_row)
//                allocate_W_row(SolnBlk);
//            if (!Assigned_W_row[theCell.I][theCell.J][theCell.K]) 
//                Assign_W_row(theCell,theNeighbours);
//            return LeastSquaresReconstruction_W_row(SolnBlk, theCell, theNeighbours);
//        }
//        else{
//            return LeastSquaresReconstruction(SolnBlk, theCell, theNeighbours);   // UNREALISTICALLY SLOW
//            
//        }
//        
//    }
//    
//}


template<typename Soln_pState, typename Soln_cState>
DenseMatrix Haselbacher_Filter<Soln_pState,Soln_cState>::Matrix_A(Cell3D &theCell, Neighbours &theNeighbours, int &commutation_order) {
    
    int the_number_of_neighbours = theNeighbours.number_of_neighbours;
    int the_number_of_unknowns = number_of_unknowns(commutation_order);
    DenseMatrix A(the_number_of_neighbours,the_number_of_unknowns);
    int j=0;
    double dx, dy, dz, coefficient;
    for (int k=0; k<=commutation_order; k++) {
        
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


template<typename Soln_pState, typename Soln_cState>
inline DenseMatrix Haselbacher_Filter<Soln_pState,Soln_cState>::Matrix_A(Cell3D &theCell, Neighbours &theNeighbours) {
    int the_number_of_neighbours = theNeighbours.number_of_neighbours;
    DenseMatrix A(the_number_of_neighbours,the_number_of_unknowns);
    int j=0;
    double dx, dy, dz, coefficient;
    for (int k=0; k<=commutation_order; k++) {
        
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



template<typename Soln_pState, typename Soln_cState>
DiagonalMatrix Haselbacher_Filter<Soln_pState,Soln_cState>::Matrix_W(Cell3D &theCell, Neighbours &theNeighbours, double weight_factor) {
    
    int the_number_of_neighbours = theNeighbours.number_of_neighbours;
    DiagonalMatrix W(the_number_of_neighbours);
    
    
    Vector3D Delta(theNeighbours.Delta);
    double dr;
    
    Delta *= weight_factor;
    
    for (int i=0; i<the_number_of_neighbours; i++) {
        dr=(theNeighbours.neighbour[i].Xc - theCell.Xc).sqr();
        W(i) = sqrt(SIX/(PI*Delta.sqr())*exp(- SIX * dr/Delta.sqr())) ;
    }
    
    return W;
}

template<typename Soln_pState, typename Soln_cState>
inline DiagonalMatrix Haselbacher_Filter<Soln_pState,Soln_cState>::Matrix_W(Cell3D &theCell, Neighbours &theNeighbours) {
    return Matrix_W(theCell, theNeighbours, weight_factor);
}



template<typename Soln_pState, typename Soln_cState>
inline DenseMatrix Haselbacher_Filter<Soln_pState,Soln_cState>::Matrix_b(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Neighbours &theNeighbours) {
    
    
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
void Haselbacher_Filter<Soln_pState,Soln_cState>::Weight_Matrix_A_and_b(DenseMatrix &A, DenseMatrix &b,Cell3D &theCell, Neighbours &theNeighbours) {
    DiagonalMatrix W = Matrix_W(theCell, theNeighbours);
    A = W*A;
    b = W*b;
}

/* 
 * A x = b  Least Squares  and filtered variable is first element of x
 */
template<typename Soln_pState, typename Soln_cState>
RowVector Haselbacher_Filter<Soln_pState,Soln_cState>::LeastSquaresReconstruction(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell, Neighbours &theNeighbours) {
    
    /* ---------------- Matrix A ---------------- */
    DenseMatrix A = Matrix_A(theCell,theNeighbours);
    
    /* ---------------- Matrix b ---------------- */
    DenseMatrix b = Matrix_b(SolnBlk,theNeighbours);
    
    /* --------------- Weighting ---------------- */
    Weight_Matrix_A_and_b(A, b, theCell, theNeighbours);
    
    /* ------------ LS procedure for Unknowns DenseMatrix x ---------- */
    DenseMatrix x(the_number_of_unknowns,b.get_n());
    int krank;
    ColumnVector Rnorm(the_number_of_unknowns);
    Solve_LS_Householder(A,b,x,krank,Rnorm);
    RowVector xrow = x[0];
    
    /* --------------- relaxation -------------- */
    RowVector x0;
    LES_Filter<Soln_pState,Soln_cState>::what_to_filter(SolnBlk,theCell.I,theCell.J,theCell.K,x0);
    double w0 = relaxation_factor;
    xrow = x0*w0 + xrow*(ONE-w0);
    
    /* ------- return filtered variables for theCell ------- */
    return xrow;
}


template<typename Soln_pState, typename Soln_cState>
inline RowVector Haselbacher_Filter<Soln_pState,Soln_cState>::LeastSquaresReconstruction_Uniform_Grid(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell, Neighbours &theNeighbours) {
    
    /* ---------------- Matrix b ---------------- */
    DenseMatrix b = Matrix_b(SolnBlk,theNeighbours);
    
    /* ----------- unknowns directly solved ------------ */
    RowVector x = w_Uniform_Grid * b;
    
    /* --------------- relaxation -------------- */
    RowVector x0;
    LES_Filter<Soln_pState,Soln_cState>::what_to_filter(SolnBlk,theCell.I,theCell.J,theCell.K,x0);
    double w0 = relaxation_factor;
    x = x0*w0 + x*(ONE-w0);
    
    /* ------- return filtered variables for theCell ------- */
    return x;
}


template<typename Soln_pState, typename Soln_cState>
inline RowVector Haselbacher_Filter<Soln_pState,Soln_cState>::LeastSquaresReconstruction_W_row(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell, Neighbours &theNeighbours) {
    
    /* ---------------- Matrix b ---------------- */
    DenseMatrix b = Matrix_b(SolnBlk,theNeighbours);
    
    /* ----------- unknowns directly solved ------------ */
    RowVector x = W_row[theCell.I][theCell.J][theCell.K] * b;
    
    /* --------------- relaxation -------------- */
    RowVector x0;
    LES_Filter<Soln_pState,Soln_cState>::what_to_filter(SolnBlk,theCell.I,theCell.J,theCell.K,x0);
    double w0 = relaxation_factor;
    x = x0*w0 + x*(ONE-w0);
    
    /* ------- return filtered variables for theCell ------- */
    return x;
}

template<typename Soln_pState, typename Soln_cState>
void Haselbacher_Filter<Soln_pState,Soln_cState>::Assign_w_Uniform_Grid(Cell3D &theCell, Neighbours &theNeighbours) {
    DenseMatrix A = Matrix_A(theCell,theNeighbours);
    DiagonalMatrix W = Matrix_W(theCell,theNeighbours,weight_factor);
    DenseMatrix Z = (W*A).pseudo_inverse()*W;
    w_Uniform_Grid = Z[0]; // RowVector of weights
    Assigned_w_Uniform_Grid = true;
    cout << "   w_Uniform_Grid assigned   ";
}

template<typename Soln_pState, typename Soln_cState>
inline void Haselbacher_Filter<Soln_pState,Soln_cState>::allocate_W_row(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk) {
    NCi = SolnBlk.NCi;
    NCj = SolnBlk.NCj;
    NCk = SolnBlk.NCk;
    W_row = new RowVector **[NCi];
    Assigned_W_row = new bool **[NCi];
    for (int i=0; i<NCi; i++) {
        W_row[i] = new RowVector *[NCj];
        Assigned_W_row[i] = new bool *[NCj];
        for (int j=0; j<NCj; j++) {
            W_row[i][j] = new RowVector [NCk];
            Assigned_W_row[i][j] = new bool [NCk];
            for (int k=0; k<NCk; k++) {
                Assigned_W_row[i][j][k] = false;
            }
        }
    }
    Allocated_W_row = true;
}

template<typename Soln_pState, typename Soln_cState>
inline void Haselbacher_Filter<Soln_pState,Soln_cState>::deallocate_W_row(void) {
    if (Allocated_W_row){
        for (int i=0; i<NCi; i++) {
            for (int j=0; j<NCj; j++) {
                delete[] W_row[i][j];   W_row[i][j] = NULL;
                delete[] Assigned_W_row[i][j];  Assigned_W_row[i][j] = NULL;
            }
            delete[] W_row[i];   W_row[i] = NULL;
            delete[] Assigned_W_row[i];   Assigned_W_row[i] = NULL;
        }
        delete[] W_row;   W_row = NULL;
        delete[] Assigned_W_row;   Assigned_W_row = NULL;
        Allocated_W_row = false;        
    }
}

template<typename Soln_pState, typename Soln_cState>
void Haselbacher_Filter<Soln_pState,Soln_cState>::Assign_W_row(Cell3D &theCell, Neighbours &theNeighbours) {
    DenseMatrix A = Matrix_A(theCell,theNeighbours);
    DiagonalMatrix W = Matrix_W(theCell,theNeighbours,weight_factor);
    DenseMatrix Z = (W*A).pseudo_inverse()*W;
    W_row[theCell.I][theCell.J][theCell.K] = Z[0]; // RowVector of weights
    Assigned_W_row[theCell.I][theCell.J][theCell.K] = true;
}


/* ------------ calculates the number of unknowns for the given order ----------- */
template <typename Soln_pState, typename Soln_cState>
inline int Haselbacher_Filter<Soln_pState,Soln_cState>::number_of_unknowns(int &commutation_order) {
    int number = 0;
    for (int k=0; k<= commutation_order; k++) {
        number += number_of_terms_of_degree(k);
    }
    return number;
}

/* ------------ calculates the number of unknowns for the given order ----------- */
template <typename Soln_pState, typename Soln_cState>
inline int Haselbacher_Filter<Soln_pState,Soln_cState>::number_of_unknowns(void) {
    int number = 0;
    for (int k=0; k<= commutation_order; k++) {
        number += number_of_terms_of_degree(k);
    }
    return number;
}

/* ----- faculty ----- */
template <typename Soln_pState, typename Soln_cState>
int Haselbacher_Filter<Soln_pState,Soln_cState>::fac(int n) {
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
inline int Haselbacher_Filter<Soln_pState,Soln_cState>::number_of_terms_of_degree(int k){
    int n=3; // 3 dimensions
    return (  fac(n+k-1)/(fac(k)*fac(n-1))  );
}




/* ------------------------ trinomial coefficient -------------------------- */
/* 
 * coefficient of x^n1 y^n2 z^n3 in expansion of (x+y+z)^n (with n=n1+n2+n3)
 *  (n1 + n2 + n3)!             (               (n1 + n2 + n3 + ... + nk)!   )
 *  ---------------             ( multinomial:  ---------------------------  )
 *   n1!  n2!  n3!              (                 n1!  n2!  n3!  ...  nk!    )
 */
template<typename Soln_pState, typename Soln_cState>
inline double Haselbacher_Filter<Soln_pState,Soln_cState>::trinomial_coefficient(int n1, int n2, int n3){
    return (  double(fac(n1+n2+n3))/double( fac(n1) * fac(n2) * fac(n3) ) );
}




//template<typename Soln_pState, typename Soln_cState>
//inline double Haselbacher_Filter<Soln_pState,Soln_cState>::filter_quality(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell, double &kmax, int number_of_rings, int commutation_order, double weight_factor) {
//    
//    Neighbours theNeighbours(SolnBlk.Grid);
//    theNeighbours.GetNeighbours(theCell,number_of_rings);
//    
//    DenseMatrix A = Matrix_A(theCell, theNeighbours, commutation_order);
//    DiagonalMatrix W = Matrix_W(theCell, theNeighbours, weight_factor);
//    DenseMatrix Z = (W*A).pseudo_inverse()*W;
//    RowVector w = Z[0];
//    double w0 = Calculate_relaxation_factor(theCell,theNeighbours,kmax,w);
//    
//    return filter_quality(theCell,theNeighbours,kmax,w,w0);
//}


template<typename Soln_pState, typename Soln_cState>
inline double Haselbacher_Filter<Soln_pState,Soln_cState>::filter_quality(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &kmax, RowVector &w) {
    
    double sharpness  = filter_sharpness(theCell,theNeighbours,kmax,w);
    double smoothness  = filter_smoothness(theCell,theNeighbours,kmax,w,0.25);
    double uniformity = filter_uniformity(theCell,theNeighbours,kmax,w,0.3); 
    
    double quality = sharpness * smoothness * uniformity;
    
    Print_3(sharpness, smoothness, uniformity);
    
    return quality;
}


template<typename Soln_pState, typename Soln_cState>
inline double Haselbacher_Filter<Soln_pState,Soln_cState>::filter_sharpness(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &kmax, RowVector &w) {
    
    Vector3D k_90 = Calculate_wavenumber_of_Gvalue(theCell,theNeighbours,kmax,w,0.9);
    Vector3D k_10 = Calculate_wavenumber_of_Gvalue(theCell,theNeighbours,kmax,w,0.1);
    double sharpness = max(ZERO,ONE - (k_10 - k_90).abs()/(kmax.abs()/(sqrt(THREE))));
    return sqrt(sharpness);
}

template<typename Soln_pState, typename Soln_cState>
inline double Haselbacher_Filter<Soln_pState,Soln_cState>::filter_smoothness(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &kmax, RowVector &w, double max_smoothness) {
    
    Vector3D k_10 = Calculate_wavenumber_of_Gvalue(theCell,theNeighbours,kmax,w,0.1);
    Vector3D k_extreme = Calculate_wavenumber_of_dGvalue(theCell,theNeighbours,k_10,kmax,w,ZERO);
    Complex G = G_function(theCell,theNeighbours,k_extreme,w);
    double smoothness = fabs(real(G)) > max_smoothness ? ZERO : fabs(fabs(real(G))-max_smoothness)/max_smoothness;
    return smoothness;
}

template<typename Soln_pState, typename Soln_cState>
inline double Haselbacher_Filter<Soln_pState,Soln_cState>::filter_uniformity(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &kmax, RowVector &w, double max_distortion) {
    
    Vector3D k100 , k110;
    k100.zero();    k100.x = kmax.x;
    k110.zero();    k110.x = kmax.x;    k110.y = kmax.y;
    
    Complex G_111 = G_function(theCell,theNeighbours,kmax,w);
    Complex G_100 = G_function(theCell,theNeighbours,k100,w);
    Complex G_110 = G_function(theCell,theNeighbours,k100,w);
    
    double uniformity_100 = fabs(real(G_100)-real(G_111)) > max_distortion ? ZERO : fabs(fabs(real(G_100)-real(G_111))-max_distortion)/max_distortion;
    double uniformity_110 = fabs(real(G_110)-real(G_111)) > max_distortion ? ZERO : fabs(fabs(real(G_110)-real(G_111))-max_distortion)/max_distortion;
    
    return min(uniformity_100,uniformity_110);
}






template<typename Soln_pState, typename Soln_cState>
double Haselbacher_Filter<Soln_pState,Soln_cState>::
Calculate_relaxation_factor(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &kmax, RowVector &w) {
    
    double G_max = real(G_function(theCell,theNeighbours,kmax,w));
    
    double a,b,m,s,fa,fb,fm,f,w0;
    int counter;
    a=-100, b=100;
    fa = real(a+(ONE-a)*G_max);
    fb = real(b+(ONE-b)*G_max);
    f = 100;
    counter = 0;
    while( fabs(f) >= 0.0001 ) {
        m = a + (b-a)/TWO;
        fm = m+(ONE-m)*G_max;
        
        if (fa<ZERO)    s = -ONE;
        else            s = ONE;
        w0 = m + (m-a)*s*fm/sqrt(fm*fm - fa*fb);
        
        f = (w0+(ONE-w0)*G_max);
        
        if (fa*fm < ZERO) { b = m ; fb = fm; } else { a = m ; fa = fm; }; 
        if (fa*f  < ZERO) { b = w0; fb = f ; } else { a = w0; fa = f; }; 
        
        
        counter++;
        if(counter >= 100) {
            cout << "max reached for relaxation_factor" << endl;
            return 0.0;
        }
    }
    return w0;
}


//! curve fits
template<typename Soln_pState, typename Soln_cState>
inline double Haselbacher_Filter<Soln_pState,Soln_cState>::Filter_Grid_Ratio(int number_of_rings, double comm_order, double weight) {
    
    if (number_of_rings == 2) {
        if (comm_order == 2 || comm_order == 3) {
            // y = (a + b*x + c*x^2) / (1.0 + d*x + e*x^2)
            double temp;
            double x_in = weight;
            temp = 0.0;
            
            // coefficients
            double a = 1.2618132070805663E+00;
            double b = -1.5652264079596057E+00;
            double c = 1.0271073495343523E+00;
            double d = -6.8595127739351558E-01;
            double e = 4.4017996869326392E-01;
            
            temp += (a + b * x_in + c * x_in * x_in) / (1.0 + d * x_in + e * x_in * x_in);
            return temp;
        } else if (comm_order == 4 || comm_order == 5) {   // actually the same as previous... due to cancellation of symmetric scheme
            // y = a( x0.5 ) + b( x2 ) + c( x-2.0 ) + d( x-1.5 ) + offset 
            double temp;
            temp = 0.0;
            double x_in = weight;
            
            // coefficients
            double a = -5.9700306784590824E-02;
            double b = 1.4906257654202673E-04;
            double c = 4.2964779587709421E+00;
            double d = -5.6358841379335569E+00;
            double e = 2.2338885178389103E+00;
            
            temp += a * pow(x_in, 0.5);
            temp += b * pow(x_in, 2.0);
            temp += c * pow(x_in, -2.0);
            temp += d * pow(x_in, -1.5);
            temp += e;
            return temp;
        } else {
            // z = a( ln(x)tanh(y) ) + b( x^(-0.5) y^(-1.5) ) + offset
            double FGR;
            FGR = 0.0;
            
            // coefficients
            double a = -1.3199931671292526E+00;
            double b = -3.0701467221618390E+00;
            double c = 3.8963966490315074E+00;
            
            FGR += a * log(comm_order) * tanh(weight);
            FGR += b * pow(comm_order, -0.5) * pow(weight, -1.5);
            FGR += c;
            return FGR;
        }
    } else if (number_of_rings == 3) {
        if (commutation_order == 2 || commutation_order == 3) {
            double temp;
            temp = 0.0;
            double x_in = weight;
            // coefficients
            double a = -7.4547319587168670E-01;
            double b = -5.1576731798758173E+00;
            double c = -8.7194396722989236E+00;
            double d = -6.2091839207957884E+00;
            double e = 1.5705916687174923E+01;
            
            temp += a * log(x_in);
            temp += b * exp(-1.0 * x_in);
            temp += c * tanh(x_in);
            temp += d * pow(x_in, -0.5);
            temp += e;
            return temp; 
        }
        else {
            cout << "commutation_order " << commutation_order << " with " << number_of_rings << " number of rings not supported, make curve fit." << endl;
        }
    } 
    
    
    

    
}



template<typename Soln_pState, typename Soln_cState>
int Haselbacher_Filter<Soln_pState,Soln_cState>::Calculate_weight_factor(void) {
        
    if (number_of_rings == 2 || number_of_rings == 3) {
        
        // Ridder's method
        double a,b,m,s,fa,fb,fm,fp,p;
        int counter;
        a=1.0, b=10.0;
        fa = Filter_Grid_Ratio(number_of_rings,commutation_order,a) - FGR;
        fb = Filter_Grid_Ratio(number_of_rings,commutation_order,b) - FGR;
        fp = 100;
        counter = 0;
        while( fabs(fp) >= 0.001 ) {
            m = a + (b-a)/TWO;
            fm = Filter_Grid_Ratio(number_of_rings,commutation_order,m) - FGR;
            
            if (fa<ZERO)    s = -ONE;
            else            s = ONE;
            p = m + (m-a)*s*fm/(sqrt(fm*fm - fa*fb)+PICO);
            
            fp = Filter_Grid_Ratio(number_of_rings,commutation_order,p) - FGR;
            
            if (fa*fm < ZERO) { b = m; fb = fm; } else { a = m; fa = fm; }; 
            if (fa*fp < ZERO) { b = p; fb = fp; } else { a = p; fa = fp; }; 
            
            
            counter++;
            if(counter >= 10) {
                cout << "max reached for FGR" << endl;
                break;
            }
            //cout << "FGR = " << fp + FGR_target << endl;
        }
        
        weight_factor = p;
        
        if (isnan(weight_factor) || isinf(weight_factor)) {
            return 1;
        }
        return 0;
    }
    else {
        cout << "no curve-fits available for " << number_of_rings << "in Haselbacher's discrete filter. Make curve fit" << endl;
        return 1;
    }
}


template<typename Soln_pState, typename Soln_cState>
double Haselbacher_Filter<Soln_pState,Soln_cState>::Calculate_weight_factor(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &kmax, double &FGR_target, int &comm_order) {
    
    double weight;
    
    if (number_of_rings == 2 || number_of_rings == 3) {
        double a,b,m,s,fa,fb,fm,fp,p;
        int counter;
        a=1.0, b=10.0;
        fa = Filter_Grid_Ratio(number_of_rings,comm_order,a) - FGR_target;
        fb = Filter_Grid_Ratio(number_of_rings,comm_order,b) - FGR_target;
        fp = 100;
        counter = 0;
        while( fabs(fp) >= 0.001 ) {
            m = a + (b-a)/TWO;
            fm = Filter_Grid_Ratio(number_of_rings,comm_order,m) - FGR_target;
            
            if (fa<ZERO)    s = -ONE;
            else            s = ONE;
            p = m + (m-a)*s*fm/(sqrt(fm*fm - fa*fb)+PICO);
            
            fp = Filter_Grid_Ratio(number_of_rings,comm_order,p) - FGR_target;
            
            if (fa*fm < ZERO) { b = m; fb = fm; } else { a = m; fa = fm; }; 
            if (fa*fp < ZERO) { b = p; fb = fp; } else { a = p; fa = fp; }; 
            
            
            counter++;
            if(counter >= 10) {
                cout << "max reached for FGR" << endl;
                break;
            }
            //cout << "FGR = " << fp + FGR_target << endl;
        }
        
        weight = p;
        
        if (isnan(weight) || isinf(weight)) {
            cout << " error in calculating weight for Cell ("<<theCell.I<<","<<theCell.J<<","<<theCell.K<<")" << endl;
        }
        
        
    } else {
        DenseMatrix A = Matrix_A(theCell, theNeighbours);

        double a,b,m,s,fa,fb,fm,fp,p;
        int counter;
        a=1.0, b=10.0;
        fa = Filter_Grid_Ratio(theCell,theNeighbours,a,A,kmax) - FGR_target;
        fb = Filter_Grid_Ratio(theCell,theNeighbours,b,A,kmax) - FGR_target;
        fp = 100;
        counter = 0;
        while( fabs(fp) >= 0.001 ) {
            m = a + (b-a)/TWO;
            fm = Filter_Grid_Ratio(theCell,theNeighbours,m,A,kmax) - FGR_target;
            
            if (fa<ZERO)    s = -ONE;
            else            s = ONE;
            p = m + (m-a)*s*fm/sqrt(fm*fm - fa*fb);
            
            fp = Filter_Grid_Ratio(theCell,theNeighbours,p,A,kmax) - FGR_target;
            
            if (fa*fm < ZERO) { b = m; fb = fm; } else { a = m; fa = fm; }; 
            if (fa*fp < ZERO) { b = p; fb = fp; } else { a = p; fa = fp; }; 
            
            
            counter++;
            if(counter >= 10) {
                cout << "max reached for FGR" << endl;
                break;
            }
            //cout << "FGR = " << fp + FGR_target << endl;
        }
        
        weight = p;
    }

    return weight;    
}



template<typename Soln_pState, typename Soln_cState>
double Haselbacher_Filter<Soln_pState,Soln_cState>::Filter_Grid_Ratio(Cell3D &theCell, Neighbours &theNeighbours, double &weight, DenseMatrix &A, Vector3D &kmax) {
    DiagonalMatrix W = Matrix_W(theCell, theNeighbours,weight);
    DenseMatrix Z = (W*A).pseudo_inverse()*W;
    RowVector w = Z[0];
    Apply_relaxation(theCell,theNeighbours,kmax,w);
    double FGR = Filter_Grid_Ratio(theCell,theNeighbours,w,kmax);
    Undo_relaxation(theNeighbours,w);
    return FGR;
}



template<typename Soln_pState, typename Soln_cState>
void Haselbacher_Filter<Soln_pState,Soln_cState>::filter_tests(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell) {

    //Output_Filter_types(SolnBlk,theCell,kmax);
    
    //    Output_Filter_types(SolnBlk,theCell,kmax,1);
    //    Output_Filter_types(SolnBlk,theCell,kmax,2);
    //    Output_Filter_types(SolnBlk,theCell,kmax,3);
    //    Output_Filter_types(SolnBlk,theCell,kmax,4);
    Output_Filter_types(SolnBlk,theCell,3,1);
    Output_Filter_types(SolnBlk,theCell,3,3);
    Output_Filter_types(SolnBlk,theCell,3,5);
    

}




template<typename Soln_pState, typename Soln_cState>
int Haselbacher_Filter<Soln_pState,Soln_cState>::
Output_Filter_types(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell) {
    
    
    int N_number_of_rings = 5;
    int N_weight_factor = 29;   double max_weight_factor = 15.0;
    int N_commutation_order = 6;
    
    double *weight_factors = new double[N_weight_factor];
    int *commutation_orders = new int[N_commutation_order];
    int *rings = new int[N_number_of_rings];
    
    double ***FGR = new double **[N_number_of_rings];
    double ***Q = new double **[N_number_of_rings];
    double ***smoothness = new double **[N_number_of_rings];
    double ***sharpness = new double **[N_number_of_rings];
    double ***uniformity = new double **[N_number_of_rings];
    double ***cost = new double **[N_number_of_rings];
    for (int i=0; i<N_number_of_rings; i++) {
        FGR[i] = new double *[N_commutation_order];
        Q[i] = new double *[N_commutation_order];
        smoothness[i] = new double *[N_commutation_order];
        sharpness[i] = new double *[N_commutation_order];
        uniformity[i] = new double *[N_commutation_order];
        cost[i] = new double *[N_commutation_order];
        for (int j=0; j<N_commutation_order; j++) {
            FGR[i][j] = new double [N_weight_factor];
            Q[i][j] = new double [N_weight_factor];
            smoothness[i][j] = new double [N_weight_factor];
            sharpness[i][j] = new double [N_weight_factor];
            uniformity[i][j] = new double [N_weight_factor];
            cost[i][j] = new double [N_weight_factor];
        }
    }
    
    Neighbours theNeighbours(SolnBlk.Grid);
    
    
    for (int i=0; i<N_number_of_rings; i++) {
        int number_of_rings = i+1;
        rings[i] = number_of_rings;
        
        theNeighbours.GetNeighbours(theCell,number_of_rings);
        
        Vector3D Delta = theNeighbours.Delta;
        Vector3D kmax(PI/Delta.x,PI/Delta.y,PI/Delta.z);
        
        
        for (int j=0; j<N_commutation_order; j++) {
            int commutation_order = j+1;
            commutation_orders[j] = commutation_order;
            
            DenseMatrix A = Matrix_A(theCell, theNeighbours, commutation_order);
            
            
            for (int k=0; k<N_weight_factor; k++) {
                double weight_factor = 1.0 + k*(max_weight_factor-1.0)/(N_weight_factor-1.0);
                weight_factors[k]=weight_factor;
                
                //Print_3(number_of_rings,commutation_order,weight_factor);
                DiagonalMatrix W = Matrix_W(theCell, theNeighbours,weight_factor);
                DenseMatrix Z = (W*A).pseudo_inverse()*W;
                RowVector w = Z[0];
                Apply_relaxation(theCell,theNeighbours,kmax,w);
                
                FGR[i][j][k] = Filter_Grid_Ratio(theCell,theNeighbours,w,kmax);
                Q[i][j][k] = filter_quality(theCell, theNeighbours, kmax, w);
                sharpness[i][j][k] = filter_sharpness(theCell, theNeighbours, kmax, w);
                smoothness[i][j][k] = filter_smoothness(theCell, theNeighbours, kmax, w, 0.25);
                uniformity[i][j][k] = filter_uniformity(theCell, theNeighbours, kmax, w, 0.3);
                cost[i][j][k] = theNeighbours.number_of_neighbours * number_of_unknowns(commutation_order);
                Undo_relaxation(theNeighbours,w);
            }
        }
    }
    
    
    
    
    char prefix[256], extension[256], output_file_name[256];
    char *output_file_name_ptr;
    ofstream output_file;    
    
    /* Determine prefix of output data file names. */
    
    //    i = 0;
    //    while (1) {
    //        if (Input.Output_File_Name[i] == ' ' ||
    //            Input.Output_File_Name[i] == '.') break;
    //        prefix[i]=Input.Output_File_Name[i];
    //        i = i + 1;
    //        if (i > strlen(Input.Output_File_Name) ) break;
    //    } /* endwhile */
    //    prefix[i] = '\0';
    sprintf(prefix,"filter_types");
    
    /* Determine output data file name for this processor. */
    
    sprintf(extension, ".dat");
    strcpy(output_file_name, prefix);
    strcat(output_file_name, extension);
    output_file_name_ptr = output_file_name;
    
    cout << "Writing file " << output_file_name << endl;
    
    /* Open the output data file. */
    
    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.fail()) return (1);
    
    /* Write the solution data for each solution block. */
    
    output_file << setprecision(14);
    
    output_file << "TITLE = \"" << "Filter types, "
    << "\"" << "\n"
    << "VARIABLES = "
    << "\"Number of rings\" \\ \n"
    << "\"Commutation order\" \\ \n"
    << "\"Weight factor\" \\ \n"
    << "\"FGR\" \\ \n"
    << "\"filter quality\" \\ \n"
    << "\"sharpness\" \\ \n"
    << "\"smoothness\" \\ \n"
    << "\"uniformity\" \\ \n"
    << "\"cost\" \\ \n";
    output_file<< "ZONE T =  \"Block Number = 0" 
    << "\" \\ \n"
    << "I = " << N_number_of_rings << " \\ \n"
    << "J = " << N_commutation_order << " \\ \n"
    << "K = " << N_weight_factor << " \\ \n"
    << "DATAPACKING = POINT \n";
    
    
    for (int k=0; k<N_weight_factor; k++) {
        for (int j=0; j<N_commutation_order; j++) {
            for (int i=0; i<N_number_of_rings ; i++) {
                output_file << " " << rings[i] << " " << commutation_orders[j] << " " << weight_factors[k] 
                << " " << FGR[i][j][k]
                << " " << Q[i][j][k]*100
                << " " << sharpness[i][j][k]*100
                << " " << smoothness[i][j][k]*100
                << " " << uniformity[i][j][k]*100
                << " " << cost[i][j][k]
                << "\n";
            } 
        } 
    }
    
    output_file << setprecision(6);
    
    /* deallocate */
    for (int i=0; i<N_number_of_rings; i++) {
        for (int j=0; j<N_commutation_order; j++) {
            delete[] FGR[i][j];
            delete[] Q[i][j];
            delete[] smoothness[i][j];
            delete[] sharpness[i][j];
            delete[] uniformity[i][j];
            delete[] cost[i][j];
        }
        delete[] FGR[i];
        delete[] Q[i];
        delete[] smoothness[i];
        delete[] sharpness[i];
        delete[] uniformity[i];
        delete[] cost[i];
    }
    delete[] FGR;
    delete[] Q;
    delete[] smoothness;
    delete[] sharpness;
    delete[] uniformity;
    delete[] cost;
    delete[] weight_factors;
    delete[] commutation_orders;
    delete[] rings;
    
    /* Close the output data file. */
    
    output_file.close();
    
    /* Writing of output data files complete.  Return zero value. */
    
    return(0);
    
}


template<typename Soln_pState, typename Soln_cState>
int Haselbacher_Filter<Soln_pState,Soln_cState>::
Output_Filter_types(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell, int number_of_rings) {
    
    
    int N_weight_factor = 29;   double max_weight_factor = 15.0;
    int N_commutation_order = 7;
    
    double *weight_factors = new double[N_weight_factor];
    int *commutation_orders = new int[N_commutation_order];
    
    double **FGR = new double *[N_commutation_order];
    double **Q = new double *[N_commutation_order];
    double **smoothness = new double *[N_commutation_order];
    double **sharpness = new double *[N_commutation_order];
    double **uniformity = new double *[N_commutation_order];
    double **cost = new double *[N_commutation_order];
    for (int i=0; i<N_commutation_order; i++) {
        FGR[i] = new double [N_weight_factor];
        Q[i] = new double [N_weight_factor];
        smoothness[i] = new double [N_weight_factor];
        sharpness[i] = new double [N_weight_factor];
        uniformity[i] = new double [N_weight_factor];
        cost[i] = new double [N_weight_factor];
    }
    
    Neighbours theNeighbours(SolnBlk.Grid);
    theNeighbours.GetNeighbours(theCell,number_of_rings);
        
    Vector3D Delta = theNeighbours.Delta;
    Vector3D kmax(PI/Delta.x,PI/Delta.y,PI/Delta.z);
    
        
        for (int j=0; j<N_commutation_order; j+=2) {
            int commutation_order = j+1;
            commutation_orders[j] = commutation_order;
            
            DenseMatrix A = Matrix_A(theCell, theNeighbours, commutation_order);
            
            
            for (int k=0; k<N_weight_factor; k++) {
                double weight_factor = 1.0 + k*(max_weight_factor-1.0)/(N_weight_factor-1.0);
                weight_factors[k]=weight_factor;
                
                //Print_3(number_of_rings,commutation_order,weight_factor);
                DiagonalMatrix W = Matrix_W(theCell, theNeighbours,weight_factor);
                DenseMatrix Z = (W*A).pseudo_inverse()*W;
                RowVector w = Z[0];
                Apply_relaxation(theCell,theNeighbours,kmax,w);

                FGR[j][k] = Filter_Grid_Ratio(theCell,theNeighbours,w,kmax);
                Q[j][k] = filter_quality(theCell, theNeighbours, kmax, w);
                sharpness[j][k] = filter_sharpness(theCell, theNeighbours, kmax, w);
                smoothness[j][k] = filter_smoothness(theCell, theNeighbours, kmax, w, 0.25);
                uniformity[j][k] = filter_uniformity(theCell, theNeighbours, kmax, w, 0.3);
                cost[j][k] = theNeighbours.number_of_neighbours * number_of_unknowns(commutation_order);
                Undo_relaxation(theNeighbours,w);
                // FGR = FGR(weight,commutation_order,number_of_rings)
                
                //cout << "weight = " << weight_factor << "    w0 = " << w0 << "    FGR = " << FGR[i][j][k] << "    FQ = " << Q[i][j][k] << endl;
                
            }
        }
    
    
    
    
    char prefix[256], extension[256], output_file_name[256];
    char *output_file_name_ptr;
    ofstream output_file;    
    
    /* Determine prefix of output data file names. */
    
    //    i = 0;
    //    while (1) {
    //        if (Input.Output_File_Name[i] == ' ' ||
    //            Input.Output_File_Name[i] == '.') break;
    //        prefix[i]=Input.Output_File_Name[i];
    //        i = i + 1;
    //        if (i > strlen(Input.Output_File_Name) ) break;
    //    } /* endwhile */
    //    prefix[i] = '\0';
    sprintf(prefix,"filter_types_%d_rings",number_of_rings);
    
    /* Determine output data file name for this processor. */
    
    sprintf(extension, ".dat");
    strcpy(output_file_name, prefix);
    strcat(output_file_name, extension);
    output_file_name_ptr = output_file_name;
    
    /* Open the output data file. */
    cout << "Writing file " << output_file_name << endl;

    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.fail()) return (1);
    
    /* Write the solution data for each solution block. */
    
    output_file << setprecision(14);
    
    output_file << "TITLE = \"" << "Filter types, "
    << "\"" << "\n"
    << "VARIABLES = "
    << "\"Commutation order\" \\ \n"
    << "\"Weight factor\" \\ \n"
    << "\"FGR\" \\ \n"
    << "\"filter quality\" \\ \n"
    << "\"sharpness\" \\ \n"
    << "\"smoothness\" \\ \n"
    << "\"uniformity\" \\ \n"
    << "\"cost\" \\ \n";
    output_file<< "ZONE T =  \"Block Number = 0" 
    << "\" \\ \n"
    << "I = " << ceil(N_commutation_order/2.0) << " \\ \n"
    << "J = " << N_weight_factor << " \\ \n"
    << "DATAPACKING = POINT \n";
    
    
    for (int k=0; k<N_weight_factor; k++) {
        for (int j=0; j<N_commutation_order; j+=2) {
                output_file << " " << commutation_orders[j] << " " << weight_factors[k] 
                << " " << FGR[j][k]
                << " " << Q[j][k]*100
                << " " << sharpness[j][k]*100
                << " " << smoothness[j][k]*100
                << " " << uniformity[j][k]*100
                << " " << cost[j][k]
                << "\n";
        } 
    }
    
    output_file << setprecision(6);
    
    /* deallocate */
        for (int j=0; j<N_commutation_order; j++) {
            delete[] FGR[j];
            delete[] Q[j];
            delete[] smoothness[j];
            delete[] sharpness[j];
            delete[] uniformity[j];
            delete[] cost[j];
        }

    
    delete[] FGR;
    delete[] Q;
    delete[] smoothness;
    delete[] sharpness;
    delete[] uniformity;
    delete[] cost;
    delete[] weight_factors;
    delete[] commutation_orders;
    
    /* Close the output data file. */
    
    output_file.close();
    
    /* Writing of output data files complete.  Return zero value. */
    
    return(0);
    
}

template<typename Soln_pState, typename Soln_cState>
int Haselbacher_Filter<Soln_pState,Soln_cState>::
Output_Filter_types(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell, int number_of_rings, int commutation_order) {
    
    
    int N_weight_factor = 29;   double max_weight_factor = 15.0;
    
    double *weight_factors = new double[N_weight_factor];
    
    double *FGR = new double [N_weight_factor];
    double *Q = new double [N_weight_factor];
    double *smoothness = new double [N_weight_factor];
    double *sharpness = new double [N_weight_factor];
    double *uniformity = new double [N_weight_factor];
    double *cost = new double [N_weight_factor];
    
    Neighbours theNeighbours(SolnBlk.Grid);
    theNeighbours.GetNeighbours(theCell, number_of_rings);
    
    Vector3D Delta = theNeighbours.Delta;
    Vector3D kmax(PI/Delta.x,PI/Delta.y,PI/Delta.z);
    

        DenseMatrix A = Matrix_A(theCell, theNeighbours, commutation_order);
        
        
        for (int k=0; k<N_weight_factor; k++) {
            double weight_factor = 1.0 + k*(max_weight_factor-1.0)/(N_weight_factor-1.0);
            weight_factors[k]=weight_factor;
            
            //Print_3(number_of_rings,commutation_order,weight_factor);
            DiagonalMatrix W = Matrix_W(theCell, theNeighbours,weight_factor);
            DenseMatrix Z = (W*A).pseudo_inverse()*W;
            RowVector w = Z[0];
            
            Apply_relaxation(theCell,theNeighbours,kmax,w);
            
            /* ------- calculate k_HALF where G(k_HALF) = 0.5 -------- */
            FGR[k] = Filter_Grid_Ratio(theCell,theNeighbours,w,kmax);
            Q[k] = filter_quality(theCell, theNeighbours, kmax, w);
            sharpness[k] = filter_sharpness(theCell, theNeighbours, kmax, w);
            smoothness[k] = filter_smoothness(theCell, theNeighbours, kmax, w, 0.25);
            uniformity[k] = filter_uniformity(theCell, theNeighbours, kmax, w, 0.3);
            cost[k] = theNeighbours.number_of_neighbours * number_of_unknowns(commutation_order);
            
            
            Undo_relaxation(theNeighbours, w);
            // FGR = FGR(weight,commutation_order,number_of_rings)
            
            //cout << "weight = " << weight_factor << "    w0 = " << w0 << "    FGR = " << FGR[i][j][k] << "    FQ = " << Q[i][j][k] << endl;
            
        }

    
    
    
    
    char prefix[256], extension[256], output_file_name[256];
    char *output_file_name_ptr;
    ofstream output_file;    
    
    /* Determine prefix of output data file names. */
    
    //    i = 0;
    //    while (1) {
    //        if (Input.Output_File_Name[i] == ' ' ||
    //            Input.Output_File_Name[i] == '.') break;
    //        prefix[i]=Input.Output_File_Name[i];
    //        i = i + 1;
    //        if (i > strlen(Input.Output_File_Name) ) break;
    //    } /* endwhile */
    //    prefix[i] = '\0';
    sprintf(prefix,"filter_types_%d_rings_order_%d",number_of_rings,commutation_order);
    
    /* Determine output data file name for this processor. */
    
    sprintf(extension, ".dat");
    strcpy(output_file_name, prefix);
    strcat(output_file_name, extension);
    output_file_name_ptr = output_file_name;
    
    /* Open the output data file. */
    cout << "Writing file " << output_file_name << endl;

    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.fail()) return (1);
    
    /* Write the solution data for each solution block. */
    
    output_file << setprecision(14);
    
    output_file << "TITLE = \"" << "Filter types, "
    << "\"" << "\n"
    << "VARIABLES = "
    << "\"Weight factor\" \\ \n"
    << "\"FGR\" \\ \n"
    << "\"filter quality\" \\ \n"
    << "\"sharpness\" \\ \n"
    << "\"smoothness\" \\ \n"
    << "\"uniformity\" \\ \n"
    << "\"cost\" \\ \n";
    output_file<< "ZONE T =  \"Block Number = 0" 
    << "\" \\ \n"
    << "I = " << N_weight_factor << " \\ \n"
    << "DATAPACKING = POINT \n";
    
    
    for (int k=0; k<N_weight_factor; k++) {
            output_file << " " << weight_factors[k] 
            << " " << FGR[k]
            << " " << Q[k]*100
            << " " << sharpness[k]*100
            << " " << smoothness[k]*100
            << " " << uniformity[k]*100
            << " " << cost[k]
            << "\n";
    }
    
    output_file << setprecision(6);
    
    
    delete[] FGR;
    delete[] Q;
    delete[] smoothness;
    delete[] sharpness;
    delete[] uniformity;
    delete[] cost;
    delete[] weight_factors;
    
    /* Close the output data file. */
    
    output_file.close();
    
    /* Writing of output data files complete.  Return zero value. */
    
    return(0);
    
}


#endif
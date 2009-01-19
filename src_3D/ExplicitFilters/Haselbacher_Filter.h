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

/**
 * CLASS: Haselbacher_Filter
 */

template<typename Soln_pState, typename Soln_cState>
class Haselbacher_Filter : public Discrete_Filter<Soln_pState,Soln_cState> {
public:
    
//    using Discrete_Filter<Soln_pState,Soln_cState>::theNeighbours;
//    using Discrete_Filter<Soln_pState,Soln_cState>::number_of_rings;
//    using Discrete_Filter<Soln_pState,Soln_cState>::commutation_order;
//    using Discrete_Filter<Soln_pState,Soln_cState>::FGR;
//    using Discrete_Filter<Soln_pState,Soln_cState>::Neighbouring_Values;
//    using Discrete_Filter<Soln_pState,Soln_cState>::Set_Neighbouring_Values;
//    using Discrete_Filter<Soln_pState,Soln_cState>::check_filter_moments;
//    

    using Discrete_Filter<Soln_pState,Soln_cState>::theCell;
    using Discrete_Filter<Soln_pState,Soln_cState>::theNeighbours;
    using Discrete_Filter<Soln_pState,Soln_cState>::number_of_rings;
    using Discrete_Filter<Soln_pState,Soln_cState>::commutation_order;
    using Discrete_Filter<Soln_pState,Soln_cState>::FGR;
    using Discrete_Filter<Soln_pState,Soln_cState>::debug_flag;
    using Discrete_Filter<Soln_pState,Soln_cState>::batch_flag;
    using Discrete_Filter<Soln_pState,Soln_cState>::Neighbouring_Values;
    using Discrete_Filter<Soln_pState,Soln_cState>::Set_Neighbouring_Values;
    using Discrete_Filter<Soln_pState,Soln_cState>::check_filter_moments;
    using Discrete_Filter<Soln_pState,Soln_cState>::check_filter_moments_1D;
    using Discrete_Filter<Soln_pState,Soln_cState>::G_cutoff;
    using Discrete_Filter<Soln_pState,Soln_cState>::use_fixed_filter_width;
    using Discrete_Filter<Soln_pState,Soln_cState>::fixed_filter_width;
    using Discrete_Filter<Soln_pState,Soln_cState>::G_function;
    using Discrete_Filter<Soln_pState,Soln_cState>::dG_function;
    using Discrete_Filter<Soln_pState,Soln_cState>::G_function_1D;
    using Discrete_Filter<Soln_pState,Soln_cState>::G_function_1D_100;
    using Discrete_Filter<Soln_pState,Soln_cState>::G_function_1D_110;
    using Discrete_Filter<Soln_pState,Soln_cState>::Calculate_wavenumber_of_Gvalue;
    using Discrete_Filter<Soln_pState,Soln_cState>::Calculate_wavenumber_of_dGvalue;
    using Discrete_Filter<Soln_pState,Soln_cState>::Filter_Grid_Ratio;    
    using Discrete_Filter<Soln_pState,Soln_cState>::properties;
    
    

    Haselbacher_Filter(Explicit_Filter_Properties &explicit_filter_properties) : Discrete_Filter<Soln_pState,Soln_cState>(explicit_filter_properties) {
        Read_Properties();
    }
    
    
    
    void Read_Properties(void) {
        if (properties->Changed()) {
            Discrete_Filter<Soln_pState,Soln_cState>::Read_Basic_Properties();
            properties->Get_Property(relaxation_factor,"relaxation_factor");
            properties->Get_Property(weighting,"least_squares_filter_weighting");
            properties->Get_Property(least_squares_filter_weighting_factor,"least_squares_filter_weighting_factor");
            properties->Get_Property(uniform_grid,"uniform_grid");
            properties->Get_Property(target_filter_sharpness,"target_filter_sharpness");
            properties->Properties_Read();
            cout << endl;
            
            Assigned_w_Uniform_Grid = false;
            weight_factor = 1.5;
            the_number_of_unknowns = number_of_unknowns();
            
        }
        
    }
    
    
    
    
    
    
    
    
    
    
//    //! Constructor
//    Haselbacher_Filter(void) : Discrete_Filter<Soln_pState,Soln_cState>() {
//        Assigned_w_Uniform_Grid = false;
//        weight_factor = 1.5;
//        the_number_of_unknowns = number_of_unknowns();
//        relaxation_factor = Explicit_Filter_Properties::relaxation_factor;
//        weighting = Explicit_Filter_Properties::least_squares_filter_weighting;
//    }
    //! Destructor default
    //~Haselbacher_Filter(void) {
    //}

    void Get_Neighbours(Cell3D &theCell) {
        theNeighbours.GetNeighbours(theCell, number_of_rings, FILTER_TYPE_HASELBACHER);
        theNeighbours.append_theCell(theCell);
    }
    void Get_Neighbours_1D(Cell3D &theCell, int direction) {
        theNeighbours.GetNeighbours_1D(theCell, number_of_rings, FILTER_TYPE_HASELBACHER, direction);
        theNeighbours.append_theCell(theCell);
    }
    
private:
    
    double relaxation_factor;
    int weighting;
    double weight_factor;
    int uniform_grid;
    double target_filter_sharpness;
    double least_squares_filter_weighting_factor;
    
    
    /* ------------------- uniform grid ---------------------- */
    RowVector w_Uniform_Grid;
    bool Assigned_w_Uniform_Grid;
    
    /* --------- functions for Matrix A calculations --------- */ 
    int number_of_unknowns();
    int number_of_unknowns(int &commutation_order);
    int the_number_of_unknowns;
    int number_of_terms_of_degree(int order);
    int fac(int n);
    double trinomial_coefficient(int n1, int n2, int n3);
    
    
    /* -------------- Least Squares Reconstruction ----------- */
    string Reconstruction_Equation(int &commutation_order);
    DenseMatrix Matrix_A(Cell3D &theCell, Neighbours &theNeighbours);
    DenseMatrix Matrix_A(Cell3D &theCell, Neighbours &theNeighbours, int &commutation_order);
    DiagonalMatrix Matrix_W(Cell3D &theCell, Neighbours &theNeighbours);
    DiagonalMatrix Matrix_W(Cell3D &theCell, Neighbours &theNeighbours, double weight_factor);

    RowVector Get_Weights(Cell3D &theCell, Neighbours &theNeighbours);
    RowVector Get_Weights_1D(Cell3D &theCell, Neighbours &theNeighbours,int direction);
    void Apply_relaxation(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &kmax, RowVector &w);
    void Undo_relaxation(Neighbours &theNeighbours, RowVector &w);
    
    /* --------------- Transfer function calculations ------------- */

    
  //  void Optimize_Filter(Cell3D &theCell, Neighbours &theNeighbours, double &kmax, double &FGR, int &commutation_order);
    double Calculate_weight_factor(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &kmax, double &FGR, int &commutation_order);
    int Calculate_weight_factor(void);
    double Filter_Grid_Ratio(Cell3D &theCell, Neighbours &theNeighbours, double &weight, DenseMatrix &A, Vector3D &kmax);
    double Filter_Grid_Ratio(int number_of_rings, double commutation_order, double weight);
    double Calculate_relaxation_factor(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &kmax, RowVector &w);
    double filter_quality(Grid3D_Hexa_Block &Grid_Blk, Cell3D &theCell, double &kmax, int number_of_rings, int commutation_order, double weight_factor);
    double filter_quality(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &kmax, RowVector &w);
    double filter_sharpness(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &kmax, RowVector &w);
    double filter_smoothness(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &kmax, RowVector &w, double max_smoothness);
    double filter_uniformity(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &kmax, RowVector &w, double max_distortion);

    int Output_Filter_types(Grid3D_Hexa_Block &Grid_Blk, Cell3D &theCell);
    int Output_Filter_types(Grid3D_Hexa_Block &Grid_Blk, Cell3D &theCell,  int number_of_rings);
    int Output_Filter_types(Grid3D_Hexa_Block &Grid_Blk, Cell3D &theCell, int number_of_rings, int commutation_order);

    double filter_moment(Cell3D &theCell, Neighbours &theNeighbours, RowVector &w, int q, int r, int s);
    double filter_moment_1D(Cell3D &theCell, Neighbours &theNeighbours, RowVector &w, int q, int direction);
    
    void filter_tests(Grid3D_Hexa_Block &Grid_Blk, Cell3D &theCell);

    int filter_type(void) { return FILTER_TYPE_HASELBACHER; }
    string filter_name(void) { return "Haselbacher"; }
    
    
};

/******************************************************************************************************************************/
template <typename Soln_pState, typename Soln_cState>
double Haselbacher_Filter<Soln_pState,Soln_cState>::filter_moment(Cell3D &theCell, Neighbours &theNeighbours, RowVector &w, int q, int r, int s) {
    
    double M(0);
    Vector3D X0(theCell.Xc);
    Vector3D dX;
    for(int n=0; n<theNeighbours.number_of_neighbours; n++) {
        dX = theNeighbours.neighbour[n].Xc - X0;
        M += w(n) * pow(dX.x,double(q)) * pow(dX.y,double(r)) * pow(dX.z,double(s)) ;
    }
    return M;
}

template <typename Soln_pState, typename Soln_cState>
double Haselbacher_Filter<Soln_pState,Soln_cState>::filter_moment_1D(Cell3D &theCell, Neighbours &theNeighbours, RowVector &w, int q, int direction) {
    double (Vector3D::*dX_member);
    switch(direction){
        case X_DIRECTION:
            dX_member = &Vector3D::x;
            break;
        case Y_DIRECTION:
            dX_member = &Vector3D::y;
            break;
        case Z_DIRECTION:
            dX_member = &Vector3D::z;
            break;
    }
    double M(0);
    Vector3D X0(theCell.Xc);
    Vector3D dX;
    for(int n=0; n<theNeighbours.number_of_neighbours; n++) {
        dX = theNeighbours.neighbour[n].Xc - X0;
        double wn = w(n);
        M += wn * pow(dX.*dX_member,double(q)) ;
    }
    return M;
}


template <typename Soln_pState, typename Soln_cState>
inline RowVector Haselbacher_Filter<Soln_pState,Soln_cState>::Get_Weights(Cell3D &theCell, Neighbours &theNeighbours) {
    
    Read_Properties();
    
    if (uniform_grid && theNeighbours.symmetric_stencil) {
        if (Assigned_w_Uniform_Grid)
            return w_Uniform_Grid;
    }
    
    
    theNeighbours.delete_theCell();
    int error_flag;
    DenseMatrix A = Matrix_A(theCell,theNeighbours);
    error_flag = Calculate_weight_factor();
    if (error_flag) {
        cout << " error in calculating weight for Cell ("<<theCell.I<<","<<theCell.J<<","<<theCell.K<<")" << endl;
        exit(1);
    }
    
    DiagonalMatrix W = Matrix_W(theCell, theNeighbours, weight_factor);
    DenseMatrix Z = (W*A).pseudo_inverse()*W;
    RowVector w = Z[0];
    
    Vector3D kmax;
    kmax.x = PI/theNeighbours.Delta.x;
    kmax.y = PI/theNeighbours.Delta.y;
    kmax.z = PI/theNeighbours.Delta.z;
    Apply_relaxation(theCell, theNeighbours, kmax, w);
    
    if (uniform_grid && theNeighbours.symmetric_stencil) {
        w_Uniform_Grid = w;
        Assigned_w_Uniform_Grid = true;
    }
    return w;
}


template <typename Soln_pState, typename Soln_cState>
inline RowVector Haselbacher_Filter<Soln_pState,Soln_cState>::Get_Weights_1D(Cell3D &theCell, Neighbours &theNeighbours,int direction){
    return RowVector(0);
}

template <typename Soln_pState, typename Soln_cState>
inline void Haselbacher_Filter<Soln_pState,Soln_cState>::Apply_relaxation(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &kmax, RowVector &w){
    double w00;
    if (relaxation_factor==DEFAULT) {
        w00 = Calculate_relaxation_factor(theCell,theNeighbours,kmax,w);
    } else {
        w00 = relaxation_factor;
    }
    w *= (ONE-w00);
    w.append(w00);
    theNeighbours.append_theCell(theCell);
}

template <typename Soln_pState, typename Soln_cState>
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



template <typename Soln_pState, typename Soln_cState>
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

template <typename Soln_pState, typename Soln_cState>
string Haselbacher_Filter<Soln_pState,Soln_cState>::Reconstruction_Equation(int &commutation_order) {
    
    int the_number_of_neighbours = theNeighbours.number_of_neighbours;
    int the_number_of_unknowns = number_of_unknowns(commutation_order);
    DenseMatrix A(the_number_of_neighbours,the_number_of_unknowns);
    int j=0;
    string dx="dx", dy="dy", dz="dz";
    
    string equation = "\\phi_i = \\overline{\\phi}_0";
    for (int k=1; k<=commutation_order; k++) {
        
        // for all terms of degree k
        for (int n3=0; n3<=k; n3++) {
            for (int n2=0; n2<=k; n2++) {
                for (int n1=0; n1<=k; n1++) {
                    if (n1+n2+n3 == k) {
                        
                        std::stringstream term, numerator, denominator;
                        if (k==1) {
                            term << " + " ;
                            numerator << "\\partial \\phi_0";
                        }
                        else {
                            term << " + \\frac{" << trinomial_coefficient(n1,n2,n3) << "}{" << k << "!} ";
                            numerator << "\\partial^" << k << "\\phi_0";
                        }
                        
                        if (n1>0) {
                            if (n1==1)                        
                                term << " \\Delta x_{0i}";
                            else
                                term << " \\Delta x_{0i}^" << n1;
                        }
                        if (n2>0) {
                            if (n2==1)                        
                                term << " \\Delta y_{0i}";
                            else
                                term << " \\Delta y_{0i}^" << n2;
                        }
                        if (n3>0) {
                            if (n3==1)                        
                                term << " \\Delta z_{0i}";
                            else
                                term << " \\Delta z_{0i}^" << n3;
                        }
                        
                        if (n1>0) {
                            if (n1==1)
                                denominator << " \\partial x";
                            else
                                denominator << " \\partial x^" << n1;
                        }
                        if (n2>0) {
                            if (n2==1)
                                denominator << " \\partial y";
                            else
                                denominator << " \\partial y^" << n2;
                        }
                        if (n3>0) {
                            if (n3==1)
                                denominator << " \\partial z";
                            else
                                denominator << " \\partial z^" << n3;
                        }
                        
                        term << "\\frac{" << numerator.str() << "}{" << denominator.str() << "}";
                        
                        equation += term.str();
                        term.str(std::string());
                        numerator.str(std::string());
                        denominator.str(std::string());

                    }
                }
            }
        } 
    }
    return equation;
}

template <typename Soln_pState, typename Soln_cState>
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



template <typename Soln_pState, typename Soln_cState>
inline DiagonalMatrix Haselbacher_Filter<Soln_pState,Soln_cState>::Matrix_W(Cell3D &theCell, Neighbours &theNeighbours, double weight_factor) {
    
    int the_number_of_neighbours = theNeighbours.number_of_neighbours;
    DiagonalMatrix W(the_number_of_neighbours);
    
    if (weighting==ON) {
        Vector3D Delta;
        Vector3D DR;
        if (use_fixed_filter_width) {
            Delta = Vector3D(fixed_filter_width,fixed_filter_width,fixed_filter_width) * weight_factor;
            for (int i=0; i<the_number_of_neighbours; i++) {
                DR = (theNeighbours.neighbour[i].Xc - theCell.Xc);
                W(i) = sqrt(SIX/(PI*Delta.sqr()))*exp(- (sqr(DR.x/Delta.x)+sqr(DR.y/Delta.y)+sqr(DR.z/Delta.z)) ) ;
            }
        } else {
            Delta = theNeighbours.Delta*weight_factor;
            for (int i=0; i<the_number_of_neighbours; i++) {
                DR = (theNeighbours.neighbour[i].Xc - theCell.Xc);
                W(i) = sqrt(SIX/(PI*Delta.sqr()))*exp(- (sqr(DR.x/Delta.x)+sqr(DR.y/Delta.y)+sqr(DR.z/Delta.z)) ) ;
            }            
        }
    } else {
        W.identity();
    }

    
    return W;
}

template <typename Soln_pState, typename Soln_cState>
inline DiagonalMatrix Haselbacher_Filter<Soln_pState,Soln_cState>::Matrix_W(Cell3D &theCell, Neighbours &theNeighbours) {
    return Matrix_W(theCell, theNeighbours, weight_factor);
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
template <typename Soln_pState, typename Soln_cState>
inline double Haselbacher_Filter<Soln_pState,Soln_cState>::trinomial_coefficient(int n1, int n2, int n3){
    return (  double(fac(n1+n2+n3))/double( fac(n1) * fac(n2) * fac(n3) ) );
}




//template <typename Soln_pState, typename Soln_cState>
//double Haselbacher_Filter<Soln_pState,Soln_cState>::filter_quality(Grid3D_Hexa_Block &Grid_Blk, Cell3D &theCell, double &kmax, int number_of_rings, int commutation_order, double weight_factor) {
//    
//    Neighbours theNeighbours(Grid_Blk.Grid);
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


template <typename Soln_pState, typename Soln_cState>
double Haselbacher_Filter<Soln_pState,Soln_cState>::filter_quality(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &kmax, RowVector &w) {
    
    double sharpness  = filter_sharpness(theCell,theNeighbours,kmax,w);
    double smoothness  = filter_smoothness(theCell,theNeighbours,kmax,w,0.25);
    double uniformity = filter_uniformity(theCell,theNeighbours,kmax,w,0.3); 
    
    double quality = sharpness * smoothness * uniformity;
    
    Print_3(sharpness, smoothness, uniformity);
    
    return quality;
}


template <typename Soln_pState, typename Soln_cState>
double Haselbacher_Filter<Soln_pState,Soln_cState>::filter_sharpness(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &kmax, RowVector &w) {
    
    Vector3D k_90 = Calculate_wavenumber_of_Gvalue(theCell,theNeighbours,kmax,w,0.9);
    Vector3D k_10 = Calculate_wavenumber_of_Gvalue(theCell,theNeighbours,kmax,w,0.1);
    double sharpness = max(ZERO,ONE - (k_10 - k_90).abs()/(kmax.abs()/(sqrt(THREE))));
    return sqrt(sharpness);
}

template <typename Soln_pState, typename Soln_cState>
double Haselbacher_Filter<Soln_pState,Soln_cState>::filter_smoothness(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &kmax, RowVector &w, double max_smoothness) {
    
    Vector3D k_10 = Calculate_wavenumber_of_Gvalue(theCell,theNeighbours,kmax,w,0.1);
    Vector3D k_extreme = Calculate_wavenumber_of_dGvalue(theCell,theNeighbours,k_10,kmax,w,ZERO);
    Complex G = G_function(theCell,theNeighbours,k_extreme,w);
    double smoothness = fabs(real(G)) > max_smoothness ? ZERO : fabs(fabs(real(G))-max_smoothness)/max_smoothness;
    return smoothness;
}

template <typename Soln_pState, typename Soln_cState>
double Haselbacher_Filter<Soln_pState,Soln_cState>::filter_uniformity(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &kmax, RowVector &w, double max_distortion) {
    
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






template <typename Soln_pState, typename Soln_cState>
double Haselbacher_Filter<Soln_pState,Soln_cState>::
Calculate_relaxation_factor(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &kmax, RowVector &w) {
    
    double G_max = real(G_function(theCell,theNeighbours,kmax,w));
    
    double a,b,m,s,fa,fb,fm,f,w0;
    int counter;
    a=-100, b=100;
    fa = (a+(ONE-a)*G_max);
    fb = (b+(ONE-b)*G_max);
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
            if (!batch_flag)
                cout << "max reached for relaxation_factor" << endl;
            return 0.0;
        }
    }
    return w0;
}


//! curve fits
template <typename Soln_pState, typename Soln_cState>
inline double Haselbacher_Filter<Soln_pState,Soln_cState>::Filter_Grid_Ratio(int number_of_rings, double comm_order, double weight) {
/* 
 * Curve fits made with zunzun.com website 
 * Generate Data FGR versus weight for the desired commutation error and number of neighbouring rings
 * through Output_Filter_types(Grid_Blk, theCell, number_of_rings, commutation_order)
 * in filter_tests(Grid_Blk, theCell)
 */
    double x_in = weight;
    double temp(0.0);
    double a,b,c,d,e; // coefficients

    switch(number_of_rings) {
        case 2:
            switch (commutation_order) {
                case 2:
                case 3:
                    a = 1.2618132070805663E+00;
                    b = -1.5652264079596057E+00;
                    c = 1.0271073495343523E+00;
                    d = -6.8595127739351558E-01;
                    e = 4.4017996869326392E-01;
                    temp += (a + b * x_in + c * x_in * x_in) / (1.0 + d * x_in + e * x_in * x_in);
                    return temp;
                case 4:
                case 5:
                    a = -5.9700306784590824E-02;
                    b = 1.4906257654202673E-04;
                    c = 4.2964779587709421E+00;
                    d = -5.6358841379335569E+00;
                    e = 2.2338885178389103E+00;
                    temp += a * pow(x_in, 0.5);
                    temp += b * pow(x_in, 2.0);
                    temp += c * pow(x_in, -2.0);
                    temp += d * pow(x_in, -1.5);
                    temp += e;
                    return temp;
                default:
                    cout << "commutation_order " << commutation_order << 
                    " with " << number_of_rings << " number of rings not supported, make curve fit." << endl;
                    return ZERO;
            }
        case 3:
            switch (commutation_order) {
                case 2:
                case 3:
                    a = -7.4547319587168670E-01;
                    b = -5.1576731798758173E+00;
                    c = -8.7194396722989236E+00;
                    d = -6.2091839207957884E+00;
                    e = 1.5705916687174923E+01;
                    temp += a * log(x_in);
                    temp += b * exp(-1.0 * x_in);
                    temp += c * tanh(x_in);
                    temp += d * pow(x_in, -0.5);
                    temp += e;
                    return temp; 
                case 4:
                case 5:
                    a = 5.8366080575745464E+01;
                    b = 5.9578033255065705E+01;
                    c = -1.2855074423429333E+00;
                    d = -1.5683703991817932E+01;
                    e = -8.7921184020434339E+01;
                    temp += a * atan(x_in);
                    temp += b * pow(x_in, -1.0);
                    temp += c * tanh(x_in);
                    temp += d * pow(x_in, -2.0);
                    temp += e;
                    return temp;
                default:
                    cout << "commutation_order " << commutation_order << 
                    " with " << number_of_rings << " number of rings not supported, make curve fit." << endl;
                    return ZERO;
            }
        case 4:
            switch (commutation_order) {
                case 2:
                case 3:
                    a = 1.3194950599947615E-02;
                    b = -9.9714656466242175E+00;
                    c = -9.7502277171658793E+00;
                    d = -7.8971169401914315E-02;
                    e = 1.8375899683424766E+01;
                    temp += a * pow(x_in, 2.0);
                    temp += b * tanh(x_in);
                    temp += c * pow(x_in, -0.5);
                    temp += d * pow(x_in, 1.5);
                    temp += e;
                    return temp;
                case 4:
                case 5:
                    a = -3.4725296572424318E+00;
                    b = 4.5446201107087303E+00;
                    c = -6.1712195629000455E+00;
                    d = 2.7754065845661348E-02;
                    e = 8.9774240210873870E+00;
                    temp += a * pow(x_in, 0.5);
                    temp += b * log(x_in);
                    temp += c * tanh(x_in);
                    temp += d * pow(x_in, 1.5);
                    temp += e;
                    return temp;
                case 6:
                case 7:
                    a = 1.0345864710917796E+02;
                    b = 1.0604007211535077E+02;
                    c = 4.4684065741245442E+00;
                    d = -2.8261109245951442E+01;
                    e = -1.5982055989255858E+02;
                    temp += a * atan(x_in);
                    temp += b * pow(x_in, -1.0);
                    temp += c * exp(-1.0 * x_in);
                    temp += d * pow(x_in, -2.0);
                    temp += e;
                    return temp;
                default:
                    cout << "commutation_order " << commutation_order << 
                    " with " << number_of_rings << " number of rings not supported, make curve fit." << endl;
                    return ZERO;
            }
        case 5:
            switch (commutation_order) {
                case 6:
                case 7:
                    a = -3.7686465580425743E-01;
                    b = 6.9777909729062601E-03;
                    c = 2.8315002288423621E+00;
                    d = -5.8241253402595605E+00;
                    e = 5.6635934860988550E+00;
                    temp += a * x_in;
                    temp += b * pow(x_in, 2.0);
                    temp += c * log(x_in);
                    temp += d * tanh(x_in);
                    temp += e;
                    return temp;
                default:
                    cout << "commutation_order " << commutation_order << 
                    " with " << number_of_rings << " number of rings not supported, make curve fit." << endl;
                    return ZERO;
            }
        case 6:
            switch (commutation_order) {
                case 6:
                case 7:
                    a = 2.4808172161362614E+00;
                    b = -5.9317013888515309E+00;
                    c = -3.4691174352134446E-01;
                    d = 8.0955778565039394E-01;
                    temp = exp(a + (b/x_in) + c*log(x_in)) + d;
                    return temp;
                default:
                    cout << "commutation_order " << commutation_order << 
                    " with " << number_of_rings << " number of rings not supported, make curve fit." << endl;
                    return ZERO;
            }
        default:
            cout << "commutation_order " << commutation_order << 
            " with " << number_of_rings << " number of rings not supported, make curve fit." << endl;
            return ZERO;
    }
    
    
    
    
    
}



template <typename Soln_pState, typename Soln_cState>
int Haselbacher_Filter<Soln_pState,Soln_cState>::Calculate_weight_factor(void) {
    
    if (weighting==OFF) {
        weight_factor = 1.0;
        return 0;
    }
    
    if(least_squares_filter_weighting_factor!=DEFAULT){
        weight_factor = least_squares_filter_weighting_factor;
        return 0;
    }
    
    if (number_of_rings == 2 ||
        number_of_rings == 3 ||
        number_of_rings == 4 || 
        number_of_rings == 5 ||
        number_of_rings == 6) {
        
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
                if (!batch_flag)
                    cout << "max reached for FGR" << endl;
                break;
            }
        }
        
        weight_factor = p;
        
        if (isnan(weight_factor) || isinf(weight_factor)) {
            weight_factor = 1.0;
            return 1;
        }
        return 0;
    }
    else {
        cout << "no curve-fits available for " << number_of_rings << " number of rings in Haselbacher's discrete filter. Make curve fit" << endl;
        return 1;
    }
}


template <typename Soln_pState, typename Soln_cState>
double Haselbacher_Filter<Soln_pState,Soln_cState>::Calculate_weight_factor(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &kmax, double &FGR_target, int &comm_order) {
    
    if (weighting==OFF) {
        return 1.0;
    }
    
    if(least_squares_filter_weighting_factor!=DEFAULT){
        return least_squares_filter_weighting_factor;
    }
    
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
                if (!batch_flag)
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
                if (!batch_flag)
                    cout << "max reached for FGR" << endl;
                break;
            }
        }
        
        weight = p;
    }
    
    return weight;    
}



template <typename Soln_pState, typename Soln_cState>
inline double Haselbacher_Filter<Soln_pState,Soln_cState>::Filter_Grid_Ratio(Cell3D &theCell, Neighbours &theNeighbours, double &weight, DenseMatrix &A, Vector3D &kmax) {
    DiagonalMatrix W = Matrix_W(theCell, theNeighbours,weight);
    DenseMatrix Z = (W*A).pseudo_inverse()*W;
    RowVector w = Z[0];
    Apply_relaxation(theCell,theNeighbours,kmax,w);
    double FGR = Filter_Grid_Ratio(theCell,theNeighbours,w,kmax);
    Undo_relaxation(theNeighbours,w);
    return FGR;
}



template <typename Soln_pState, typename Soln_cState>
void Haselbacher_Filter<Soln_pState,Soln_cState>::filter_tests(Grid3D_Hexa_Block &Grid_Blk, Cell3D &theCell) {
    
    check_filter_moments(Grid_Blk,theCell);
    cout << "\n" << Reconstruction_Equation(commutation_order) << endl;
    
    //Output_Filter_types(Grid_Blk,theCell);
    
    //    Output_Filter_types(Grid_Blk,theCell,kmax,1);
    //    Output_Filter_types(Grid_Blk,theCell,kmax,2);
    //    Output_Filter_types(Grid_Blk,theCell,kmax,3);
    //    Output_Filter_types(Grid_Blk,theCell,kmax,4);
    //Output_Filter_types(Grid_Blk,theCell,4,2);
    //Output_Filter_types(Grid_Blk,theCell,4,6); 
    //Output_Filter_types(Grid_Blk,theCell,5,6);    
    //Output_Filter_types(Grid_Blk,theCell,6,6);    

}




template <typename Soln_pState, typename Soln_cState>
int Haselbacher_Filter<Soln_pState,Soln_cState>::
Output_Filter_types(Grid3D_Hexa_Block &Grid_Blk, Cell3D &theCell) {
    
    
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
    
    Neighbours theNeighbours(Grid_Blk);
    
    
    for (int i=0; i<N_number_of_rings; i++) {
        int number_of_rings = i+1;
        rings[i] = number_of_rings;
        
        theNeighbours.GetNeighbours(theCell,number_of_rings,FILTER_TYPE_HASELBACHER);
        
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
    
    
    
    
    char *prefix, suffix[256], extension[256], output_file_name[256];
    char *output_file_name_ptr;
    ofstream output_file;    
    
    /* Determine output data file name for this processor. */
    prefix = properties->Get_Property_string("output_file_name");
    sprintf(suffix,"_filter_types");
    strcat(prefix,suffix);
        
    sprintf(extension, ".dat");
    strcpy(output_file_name, prefix);
    strcat(output_file_name, extension);
    output_file_name_ptr = output_file_name;
    
    if (!batch_flag)
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


template <typename Soln_pState, typename Soln_cState>
int Haselbacher_Filter<Soln_pState,Soln_cState>::
Output_Filter_types(Grid3D_Hexa_Block &Grid_Blk, Cell3D &theCell, int number_of_rings) {
    
    
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
    
    Neighbours theNeighbours(Grid_Blk);
    theNeighbours.GetNeighbours(theCell,number_of_rings,FILTER_TYPE_HASELBACHER);
    
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
    
    
    
    
    char prefix[256], suffix[256], extension[256], output_file_name[256];
    char *output_file_name_ptr;
    ofstream output_file;    
    
    /* Determine output data file name for this processor. */
    strcpy(prefix,properties->Get_Property_string("output_file_name").c_str());
    sprintf(suffix,"_filter_types_%d_rings",number_of_rings);
    strcat(prefix,suffix);
    
    sprintf(extension, ".dat");
    strcpy(output_file_name, prefix);
    strcat(output_file_name, extension);
    output_file_name_ptr = output_file_name;
    
    /* Open the output data file. */
    if (!batch_flag)
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

template <typename Soln_pState, typename Soln_cState>
int Haselbacher_Filter<Soln_pState,Soln_cState>::
Output_Filter_types(Grid3D_Hexa_Block &Grid_Blk, Cell3D &theCell, int number_of_rings, int commutation_order) {
    
    
    int N_weight_factor = 29;   double max_weight_factor = 15.0;
    
    double *weight_factors = new double[N_weight_factor];
    
    double *FGR = new double [N_weight_factor];
    double *Q = new double [N_weight_factor];
    double *smoothness = new double [N_weight_factor];
    double *sharpness = new double [N_weight_factor];
    double *uniformity = new double [N_weight_factor];
    double *cost = new double [N_weight_factor];
    
    Neighbours theNeighbours(Grid_Blk);
    theNeighbours.GetNeighbours(theCell, number_of_rings, FILTER_TYPE_HASELBACHER);
    
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
    
    
    
    
    
    char prefix[256], suffix[256], extension[256], output_file_name[256];
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
    
    /* Determine output data file name for this processor. */
    strcpy(prefix,properties->Get_Property_string("output_file_name").c_str());
    sprintf(suffix,"_filter_types_%d_rings_order_%d",number_of_rings,commutation_order);
    
    strcat(prefix,suffix);
    
    sprintf(extension, ".dat");
    strcpy(output_file_name, prefix);
    strcat(output_file_name, extension);
    output_file_name_ptr = output_file_name;
    
    /* Open the output data file. */
    if (!batch_flag)
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

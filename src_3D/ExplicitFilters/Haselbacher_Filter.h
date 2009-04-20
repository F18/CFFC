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
            properties->Get_Property("relaxation_factor",relaxation_factor);
            properties->Get_Property("least_squares_filter_weighting",weighting);
            properties->Get_Property("least_squares_filter_weighting_factor",least_squares_filter_weighting_factor);
            properties->Get_Property("target_filter_sharpness",target_filter_sharpness);
            properties->Get_Property("filter_strength",filter_strength);
            properties->Get_Property("high_order",reconstruction_type);
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
        theNeighbours.GetNeighbours(theCell, number_of_rings, Explicit_Filter_Constants::HASELBACHER_FILTER);
        theNeighbours.append_theCell(theCell);
    }
    void Get_Neighbours_1D(Cell3D &theCell, int direction) {
        theNeighbours.GetNeighbours_1D(theCell, number_of_rings, Explicit_Filter_Constants::HASELBACHER_FILTER, direction);
        theNeighbours.append_theCell(theCell);
    }
    
private:
    
    double relaxation_factor;
    int weighting;
    double weight_factor;
    double target_filter_sharpness;
    double least_squares_filter_weighting_factor;
    double filter_strength;
    int reconstruction_type;
    
    
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
    
    DenseMatrix ComputeCellReconstructionPseudoInverse(Cell3D &theCell, Neighbours &theNeighbours);
    RowVector ComputeCellWeights(Cell3D &theCell, Neighbours &theNeighbours, DenseMatrix &Cell_LHS_Inv);
    Vector3D Geometric_Weighting_Config();
    double Geometric_Weighting(Vector3D &DeltaCellCenters, Vector3D &Delta);        

    RowVector Get_Weights(Cell3D &theCell, Neighbours &theNeighbours);
    RowVector Get_Weights_1D(Cell3D &theCell, Neighbours &theNeighbours,int direction);
    void Apply_relaxation(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &kmax, RowVector &w);
    void Undo_relaxation(Neighbours &theNeighbours, RowVector &w);
    
    /* --------------- Transfer function calculations ------------- */

    
  //  void Optimize_Filter(Cell3D &theCell, Neighbours &theNeighbours, double &kmax, double &FGR, int &commutation_order);
    double Calculate_weight_factor(Cell3D &theCell, Neighbours &theNeighbours, Vector3D &kmax, double &FGR, int &commutation_order);
    int Calculate_weight_factor(void);
    double Get_Weight_Factor_from_FGR(double FGR_in) ;
    double Filter_Grid_Ratio(Cell3D &theCell, Neighbours &theNeighbours, double &weight, DenseMatrix &A, Vector3D &kmax);
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

    int filter_type(void) { return Explicit_Filter_Constants::HASELBACHER_FILTER; }
    string filter_name(void) { return "Haselbacher"; }
    
    
};

/******************************************************************************************************************************/
template <typename Soln_pState, typename Soln_cState>
double Haselbacher_Filter<Soln_pState,Soln_cState>::filter_moment(Cell3D &theCell, Neighbours &theNeighbours, RowVector &w, int q, int r, int s) {
    
    double M(0);
    Vector3D X0(theCell.Xc);
    Vector3D dX;
    for(int n=0; n<theNeighbours.number_of_neighbours; n++) {
        dX = theNeighbours.neighbour[n]->Xc - X0;
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
        dX = theNeighbours.neighbour[n]->Xc - X0;
        double wn = w(n);
        M += wn * pow(dX.*dX_member,double(q)) ;
    }
    return M;
}


template <typename Soln_pState, typename Soln_cState>
inline RowVector Haselbacher_Filter<Soln_pState,Soln_cState>::Get_Weights(Cell3D &theCell, Neighbours &theNeighbours) {
    
    int error_flag;

    Read_Properties();
    
    // If the grid is uniform, only one set of weights needs to be calculated
    if (theNeighbours.uniform_stencil && theNeighbours.symmetric_stencil) {
        if (Assigned_w_Uniform_Grid)
            return w_Uniform_Grid;
    }
    
    // The centre cell is not part of the calculation --> delete
    theNeighbours.delete_theCell();
    
    // Calculate the weight factor defining the filter width 
    error_flag = Calculate_weight_factor();
    if (error_flag) {
        cout << " error in calculating weight for Cell ("<<theCell.I<<","<<theCell.J<<","<<theCell.K<<")" << endl;
        exit(1);
    }
    
    RowVector w(theNeighbours.number_of_neighbours);
    

    if (reconstruction_type == RECONSTRUCTION_HIGH_ORDER) {
        DenseMatrix Cell_LHS_Inv = ComputeCellReconstructionPseudoInverse(theCell,theNeighbours);
        w = ComputeCellWeights(theCell, theNeighbours, Cell_LHS_Inv);
    }
    else {
        // The system matrix used in Least Squares
        DenseMatrix A = Matrix_A(theCell,theNeighbours);
        
        // The weights matrix used in Weighted Least Squares
        DiagonalMatrix W = Matrix_W(theCell, theNeighbours, weight_factor);
        
        //                              -1
        // (W A) x = W b  -->  x = (W A)   W  b  -->  x = Z b
        //
        DenseMatrix Z = (W*A).pseudo_inverse()*W;
        
        // weights are first row of pseudo_inverse
        w = Z[0];
    }
    
    // Add weight to centre cell, to scale transfer function
    Vector3D kmax;
    kmax.x = PI/theNeighbours.Delta.x;
    kmax.y = PI/theNeighbours.Delta.y;
    kmax.z = PI/theNeighbours.Delta.z;
    Apply_relaxation(theCell, theNeighbours, kmax, w);
    
    
    if (theNeighbours.uniform_stencil && theNeighbours.symmetric_stencil) {
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
                            dx = theNeighbours.neighbour[i]->Xc.x - theCell.Xc.x;
                            dy = theNeighbours.neighbour[i]->Xc.y - theCell.Xc.y;
                            dz = theNeighbours.neighbour[i]->Xc.z - theCell.Xc.z;
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
                            dx = theNeighbours.neighbour[i]->Xc.x - theCell.Xc.x;
                            dy = theNeighbours.neighbour[i]->Xc.y - theCell.Xc.y;
                            dz = theNeighbours.neighbour[i]->Xc.z - theCell.Xc.z;
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



/*! 
 * Compute the pseudo-inverse of the left-hand-side term in the 
 * unlimited k-exact high-order reconstruction for a specified
 * computational cell based on the information provided by the
 * associated grid.
 */
template <typename Soln_pState, typename Soln_cState>
inline DenseMatrix Haselbacher_Filter<Soln_pState,Soln_cState>::ComputeCellReconstructionPseudoInverse(Cell3D &theCell,
                                                                                                Neighbours &theNeighbours){
    
    // Copied and modified from HighOrderReconstructions.h
    
    // SET VARIABLES USED IN THE RECONSTRUCTION PROCESS
    
    int StencilSize(theNeighbours.number_of_neighbours);
    int IndexSumZ, IndexSumY, IndexSumX, P1, P2, P3;
    double CombP1X, CombP2Y, CombP3Z;
    double  PowDistanceZC, PowDistanceYC, PowDistanceXC;
    int cell, i;
    double MaxWeight(0.0);
    double IntSum1(0.0), IntSum2(0.0);
    DenseMatrix Cell_LHS_Inv;
    std::vector<double> GeomWeights;
    std::vector<Vector3D> DeltaCellCenters;
    Vector3D Delta;
    TaylorDerivativesContainer<Soln_pState> CellTaylorDeriv(commutation_order);
    
    // Ensure that the LHS matrix is formated correctly.
    // Memory shouldn't be allocated here, only the dimensions should be defined properly.
    Cell_LHS_Inv.newsize(StencilSize, CellTaylorDeriv.size());
    GeomWeights.resize(StencilSize);
    DeltaCellCenters.resize(StencilSize);
    
    
    // START:   Set the LHS of the linear system
    // ***************************************************
    
    // ==== Set the geometric weight associated with the reconstructed cell
    GeomWeights[0] = 1;
    
    Delta = Geometric_Weighting_Config();
    
    // Step1. Compute the normalized geometric weights
    for (cell=0; cell<StencilSize; ++cell){ //for each neighbour cell in the stencil
        
        /* Compute the X, Y, and Z components of the distance between
         the cell centers of the neighbour and the reconstructed cell */
        DeltaCellCenters[cell] = theNeighbours.neighbour[cell]->Xc - theCell.Xc;
        
        /* Compute the geometric weight based on the centroid distance */
        GeomWeights[cell] = Geometric_Weighting(DeltaCellCenters[cell], Delta);
        
        /* Compute the maximum geometric weight (this is used for normalization) */
        MaxWeight = max(MaxWeight, GeomWeights[cell]);
    }
    
    // Step2. Set the approximate equations
    for (cell=0 ; cell<StencilSize; ++cell){ //for each neighbour cell in the stencil
        
        // compute the normalized geometric weight
        GeomWeights[cell] /= MaxWeight;
        
        Cell_LHS_Inv(cell,0) = 1.0;
        Cell_LHS_Inv(cell,0) *= GeomWeights[cell];
        
        // *** SET the matrix of the linear system (LHS) ***
        /* compute for each derivative the corresponding entry in the matrix of the linear system */
        for (i=1; i<=CellTaylorDeriv.LastElem(); ++i){
            // build the row of the matrix
            P1 = CellTaylorDeriv(i).P1();  // identify P1
            P2 = CellTaylorDeriv(i).P2();  // identify P2
            P3 = CellTaylorDeriv(i).P3();  // identify P3
            
            //----------------------------------------------
            Cell_LHS_Inv(cell,i) = 0.0;  // set sumation variable to zero
            CombP3Z = 1.0;        // the binomial coefficient "nC k" for k=0 is 1
            PowDistanceZC = 1.0;  // initialize PowDistanceZC
            
            // Compute geometric integral over the neighbour's domain
            for (IndexSumZ = 0; IndexSumZ<=P3; ++IndexSumZ){
                CombP2Y = 1.0;        // the binomial coefficient "nC k" for k=0 is 1
                PowDistanceYC = 1.0;  // initialize PowDistanceYC
                IntSum2 = 0.0;         // reset internal summation variable
                
                for (IndexSumY = 0; IndexSumY<=P2; ++IndexSumY){
                    CombP1X = 1.0;       // the binomial coefficient "nC k" for k=0 is 1
                    PowDistanceXC = 1.0; // initialize PowDistanceXC
                    IntSum1 = 0.0;	     // reset internal sumation variable
                    
                    for (IndexSumX = 0; IndexSumX<=P1; ++IndexSumX){
                        IntSum1 += ( CombP1X*PowDistanceXC*
                                    theNeighbours.neighbour[cell]->GeomCoeffValue(P1-IndexSumX,P2-IndexSumY,P3-IndexSumZ) );
                        
                        // update the binomial coefficients
                        CombP1X = (P1-IndexSumX)*CombP1X/(IndexSumX+1); // the index is still the old one => expression for "nC k+1"
                        PowDistanceXC *= DeltaCellCenters[cell].x;      // Update PowDistanceXC
                    }//endfor
                    
                    IntSum2 += CombP2Y*PowDistanceYC*IntSum1;
                    CombP2Y = (P2-IndexSumY)*CombP2Y/(IndexSumY+1); // the index is still the old one => expression for "nC k+1"
                    PowDistanceYC *= DeltaCellCenters[cell].y;      // Update PowDistanceYC
                }//endfor
                
                Cell_LHS_Inv(cell,i) += CombP3Z*PowDistanceZC*IntSum2;  // update the external sum
                
                CombP3Z = (P3-IndexSumZ)*CombP3Z/(IndexSumZ+1); // the index is still the old one => expression for "nC k+1"
                PowDistanceZC *= DeltaCellCenters[cell].z;      // Update PowDistanceYC
            }//endfor
            
            // Don't!!!! subtract the corresponding geometric moment of cell (iCell,jCell) 
            //Cell_LHS_Inv(cell,i) -= theCell.GeomCoeffValue(P1,P2,P3);
            
            // apply geometric weighting
            Cell_LHS_Inv(cell,i) *= GeomWeights[cell];
            
        } // endfor (i)
    }//endfor (cell)
    
    // STOP:   Matrix of the linear system (LHS) built. 
    //         For kExact_Reconstruction away from some special curved boundaries 
    //         the same matrix is used for all variables (same geometry) and  
    //         at every time step as long as the mesh is the same.
    // **********************************************************************
        
//    cout << "Cell_LHS = \n" << Cell_LHS_Inv << endl;

    // Compute the pseudo-inverse and override the LHS term.
    // This operation will change the dimensions of the matrix.
    Cell_LHS_Inv.pseudo_inverse_override();
//    cout << "Cell_LHS_Inv = \n" << Cell_LHS_Inv << endl;

    for (cell=0; cell<StencilSize; ++cell) {
        for (i=0; i<=CellTaylorDeriv.LastElem(); ++i) {
            Cell_LHS_Inv(i,cell) *= GeomWeights[cell];
        }
    }
    return Cell_LHS_Inv;
}

/*! 
 * Compute the pseudo-inverse of the left-hand-side term in the 
 * unlimited k-exact high-order reconstruction for a specified
 * computational cell based on the information provided by the
 * associated grid.
 */
template <typename Soln_pState, typename Soln_cState>
inline RowVector Haselbacher_Filter<Soln_pState,Soln_cState>::ComputeCellWeights(Cell3D &theCell, Neighbours &theNeighbours, DenseMatrix &Cell_LHS_Inv){
    
    // Copied and modified from HighOrderReconstructions.h
    
    // SET VARIABLES USED IN THE RECONSTRUCTION PROCESS
    
    int cell, i, P1, P2, P3;
    int StencilSize(theNeighbours.number_of_neighbours);
    RowVector Cell_Weights(StencilSize);
    TaylorDerivativesContainer<Soln_pState> CellTaylorDeriv(commutation_order);
      
    for (cell=0 ; cell<StencilSize; ++cell){ //for each neighbour cell in the stencil
        
        Cell_Weights[cell] = 0.0;
        Cell_Weights[cell] += Cell_LHS_Inv(0,cell);
        
        /* compute for each derivative the corresponding correction */
        for (i=1; i<=CellTaylorDeriv.LastElem(); ++i){
            // Get P1 P2 P3
            P1 = CellTaylorDeriv(i).P1();  // identify P1
            P2 = CellTaylorDeriv(i).P2();  // identify P2
            P3 = CellTaylorDeriv(i).P3();  // identify P3

            // correction
            Cell_Weights[cell] += theCell.GeomCoeffValue(P1,P2,P3) * Cell_LHS_Inv(i,cell);
        } // endfor (i)
    }//endfor (cell)
    
    return Cell_Weights;
}

template <typename Soln_pState, typename Soln_cState>
Vector3D Haselbacher_Filter<Soln_pState,Soln_cState>::Geometric_Weighting_Config() {
    Vector3D Delta(ZERO,ZERO,ZERO);
    if (weighting==ON) {
        if (use_fixed_filter_width) {
            
            Vector3D weight_factors(Get_Weight_Factor_from_FGR(fixed_filter_width/theNeighbours.Delta.x),
                                    Get_Weight_Factor_from_FGR(fixed_filter_width/theNeighbours.Delta.y),
                                    Get_Weight_Factor_from_FGR(fixed_filter_width/theNeighbours.Delta.z));
            Delta.x = theNeighbours.Delta.x * weight_factors.x;
            Delta.y = theNeighbours.Delta.y * weight_factors.y;
            Delta.z = theNeighbours.Delta.z * weight_factors.z;
        } else {
            Delta = theNeighbours.Delta*weight_factor;
        }
    }
    return Delta;
}


template <typename Soln_pState, typename Soln_cState>
double Haselbacher_Filter<Soln_pState,Soln_cState>::Geometric_Weighting(Vector3D &DeltaCellCenters, Vector3D &Delta) {
    double w;
    if (weighting==ON) {
        if (use_fixed_filter_width) {
            
            Vector3D weight_factors(Get_Weight_Factor_from_FGR(fixed_filter_width/theNeighbours.Delta.x),
                                    Get_Weight_Factor_from_FGR(fixed_filter_width/theNeighbours.Delta.y),
                                    Get_Weight_Factor_from_FGR(fixed_filter_width/theNeighbours.Delta.z));
            Delta.x = theNeighbours.Delta.x * weight_factors.x;
            Delta.y = theNeighbours.Delta.y * weight_factors.y;
            Delta.z = theNeighbours.Delta.z * weight_factors.z;
            
            w = sqrt(SIX/(PI*Delta.sqr()))*exp(- (sqr(DeltaCellCenters.x/Delta.x)+sqr(DeltaCellCenters.y/Delta.y)+sqr(DeltaCellCenters.z/Delta.z)) ) ;
        } else {
            w = sqrt(SIX/(PI*Delta.sqr()))*exp(- (sqr(DeltaCellCenters.x/Delta.x)+sqr(DeltaCellCenters.y/Delta.y)+sqr(DeltaCellCenters.z/Delta.z)) ) ;
        }
    } else {
        w=ONE;
    }
    return w;
}


template <typename Soln_pState, typename Soln_cState>
inline DiagonalMatrix Haselbacher_Filter<Soln_pState,Soln_cState>::Matrix_W(Cell3D &theCell, Neighbours &theNeighbours, double weight_factor) {
    
    int the_number_of_neighbours = theNeighbours.number_of_neighbours;
    DiagonalMatrix W(the_number_of_neighbours);
    
    if (weighting==ON) {
        Vector3D Delta;
        Vector3D DR;
        if (use_fixed_filter_width) {
            
            Vector3D weight_factors(Get_Weight_Factor_from_FGR(fixed_filter_width/theNeighbours.Delta.x),
                                    Get_Weight_Factor_from_FGR(fixed_filter_width/theNeighbours.Delta.y),
                                    Get_Weight_Factor_from_FGR(fixed_filter_width/theNeighbours.Delta.z));
            Delta.x = theNeighbours.Delta.x * weight_factors.x;
            Delta.y = theNeighbours.Delta.y * weight_factors.y;
            Delta.z = theNeighbours.Delta.z * weight_factors.z;
            
            for (int i=0; i<the_number_of_neighbours; i++) {
                DR = (theNeighbours.neighbour[i]->Xc - theCell.Xc);
                W(i) = sqrt(SIX/(PI*Delta.sqr()))*exp(- (sqr(DR.x/Delta.x)+sqr(DR.y/Delta.y)+sqr(DR.z/Delta.z)) ) ;
            }
        } else {
            Delta = theNeighbours.Delta*weight_factor;
            for (int i=0; i<the_number_of_neighbours; i++) {
                DR = (theNeighbours.neighbour[i]->Xc - theCell.Xc);
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
    
    double offset_factor = 1. - filter_strength;
    
    double G_max = real(G_function(theCell,theNeighbours,kmax,w));
    
    double a,b,m,s,fa,fb,fm,f,w0;
    int counter;
    a=-100, b=100;
    fa = (a+(ONE-a)*G_max) - offset_factor;
    fb = (b+(ONE-b)*G_max) - offset_factor;
    f = 100;
    counter = 0;
    while( fabs(f) >= 0.0001 ) {
        m = a + (b-a)/TWO;
        fm = m+(ONE-m)*G_max - offset_factor;
        
        if (fa<ZERO)    s = -ONE;
        else            s = ONE;
        w0 = m + (m-a)*s*fm/sqrt(fm*fm - fa*fb);
        
        f = (w0+(ONE-w0)*G_max) - offset_factor;
        
        
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
        
        weight_factor = Get_Weight_Factor_from_FGR(FGR);
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
double Haselbacher_Filter<Soln_pState,Soln_cState>::Get_Weight_Factor_from_FGR(double FGR_in) {
    /* 
     * Curve fits made with zunzun.com website 
     * Generate Data weight versus FGR for the desired commutation error and number of neighbouring rings
     * through Output_Filter_types(Grid_Blk, theCell, number_of_rings, commutation_order)
     * in filter_tests(Grid_Blk, theCell)
     * fitting function type = NIST Hahn
     */
    
    
    bool not_supported=false;
	double temp = 0.0;
    double x_in = FGR_in;
	double a,b,c,d,e,f,g;
    
	// coefficients
    
    switch (commutation_order){
        case 2:
        case 3:
            switch (number_of_rings){
                case 2:
                    // 2RINGS, ORDER2,3
                    a = -3.0140751826231790E+03;
                    b = 6.2841545987077097E+03;
                    c = -3.1317049743541656E+03;
                    d = 4.5293107921096691E+02;
                    e = 1.8146241964340607E+03;
                    f = -1.2922496009634167E+03;
                    g = 2.2980281466493403E+02;
                    break;
                case 3:
                    // 3RINGS, ORDER2,3
                    a = 3.4357518949360310E+03;
                    b = -6.3060701011942265E+03;
                    c = 1.9297038560169781E+03;
                    d = -1.4555898528024727E+02;
                    e = -2.3180770229651507E+03;
                    f = 1.0842858962506270E+03;
                    g = -1.2615464091327252E+02;
                    break;
                case 4:
                    // 4RINGS, ORDER2,3
                    a = -4.2404649518127042E+02;
                    b = 7.0955051057927221E+02;
                    c = -1.0902741777903421E+02;
                    d = -1.6941336473353654E+00;
                    e = 3.1322423354893340E+02;
                    f = -1.0567847273001135E+02;
                    g = 8.6762233720351443E+00;
                    break;
                case 5:
                    // 5RINGS, ORDER2,3
                    a = -1.1985412579860616E+09;
                    b = 1.8125645149833450E+09;
                    c = 1.2927408063997599E+08;
                    d = -5.3017437452576041E+07;
                    e = 1.0805095114115376E+09;
                    f = -2.4789600499725601E+08;
                    g = 1.2229223336981714E+07;
                    break;
                case 6:
                    // 6RINGS, ORDER2,3
                    a = 3.2061922586124655E+15;
                    b = -4.4084233884181905E+15;
                    c = -1.4971205587097480E+15;
                    d = 2.3504575565694050E+14;
                    e = -3.5419648474916920E+15;
                    f = 5.7778551364900050E+14;
                    g = -1.4271454806817041E+13;
                    break;
                default:
                    not_supported=true;
                    break;
            }
            break;
        case 4:
        case 5:
            switch (number_of_rings) {
                case 2:
                    // 2RINGS, ORDER4,5
                    a = -9.1002187107878208E+16;
                    b = 2.2217694352337606E+17;
                    c = -1.2095552722303850E+17;
                    d = 1.7737038684551338E+16;
                    e = 7.6213148553424512E+16;
                    f = -6.7032261872652584E+16;
                    g = 1.4715767437761356E+16;
                    break;
                case 3:
                    // 3RINGS, ORDER4,5
                    a = 1.5460528732337981E+28;
                    b = -3.4038946747384971E+28;
                    c = 3.9250754524436671E+28;
                    d = -9.9245625774474450E+27;
                    e = 1.4216863523877390E+28;
                    f = -4.8861228717812966E+27;
                    g = -5.2003705205911889E+25;
                    break;
                case 4:
                    // 4RINGS, ORDER4,5
                    a = 3.7027165890867844E+07;
                    b = -7.4543227123703763E+07;
                    c = 2.0395741348987468E+07;
                    d = -9.7054022492336365E+05;
                    e = -2.6962161138882220E+07;
                    f = 1.2939798619173352E+07;
                    g = -1.5345582271938729E+06;
                    break;
                case 5:
                    // 5RINGS, ORDER4,5
                    a = 1.2989111154329037E+03;
                    b = -2.4201118542869854E+03;
                    c = 3.2421503374454130E+02;
                    d = 2.5420812434345123E+01;
                    e = -1.0018822551301112E+03;
                    f = 3.7079265551148973E+02;
                    g = -3.3155930667554379E+01;
                    break;
                case 6:
                    // 6RINGS, ORDER4,5
                    a = -3.5032203193030575E+04;
                    b = 6.0149046859355643E+04;
                    c = 2.9083464693638043E+03;
                    d = -2.0036206825518561E+03;
                    e = 3.0155713796316581E+04;
                    f = -8.2943845958678758E+03;
                    g = 5.0875356408944953E+02;
                    break;
                default:
                    not_supported=true;
                    break;
            }
            break;
        default:
            not_supported=true;
            break;
    }

    if (not_supported) {
        cout << "commutation_order " << commutation_order << 
        " with " << number_of_rings << " number of rings not supported, make curve fit." << endl;
        return ZERO;
    } else {
        temp += (a + b * x_in + c * x_in * x_in + d * x_in * x_in * x_in) / (1.0 + e * x_in + f * x_in * x_in + g * x_in * x_in * x_in);
        return temp;
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
    //Output_Filter_types(Grid_Blk,theCell,2,commutation_order);
//    Output_Filter_types(Grid_Blk,theCell,3,commutation_order);
//    Output_Filter_types(Grid_Blk,theCell,4,commutation_order); 
//    Output_Filter_types(Grid_Blk,theCell,5,commutation_order);    
//    Output_Filter_types(Grid_Blk,theCell,6,commutation_order);    

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
        
        theNeighbours.GetNeighbours(theCell,number_of_rings,Explicit_Filter_Constants::HASELBACHER_FILTER);
        
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
    theNeighbours.GetNeighbours(theCell,number_of_rings,Explicit_Filter_Constants::HASELBACHER_FILTER);
    
    Vector3D Delta = theNeighbours.Delta;
    Vector3D kmax(PI/Delta.x,PI/Delta.y,PI/Delta.z);
    
    
    for (int j=0; j<N_commutation_order; j+=2) {
        int commutation_order = j+1;
        commutation_orders[j] = commutation_order;
        
        DenseMatrix A = Matrix_A(theCell, theNeighbours, commutation_order);
        
        
        for (int k=0; k<N_weight_factor; k++) {
            double weight_factor = 1.0 + k*(max_weight_factor-1.0)/(N_weight_factor-1.0);
            weight_factors[k]=weight_factor;
            
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
    
    
    int N_weight_factor = 29;   
    double min_weight_factor = 0.5;
    double max_weight_factor = 15.0;
    
    double *weight_factors = new double[N_weight_factor];
    
    double *FGR = new double [N_weight_factor];
    double *Q = new double [N_weight_factor];
    double *smoothness = new double [N_weight_factor];
    double *sharpness = new double [N_weight_factor];
    double *uniformity = new double [N_weight_factor];
    double *cost = new double [N_weight_factor];
    
    Neighbours theNeighbours(Grid_Blk);
    theNeighbours.GetNeighbours(theCell, number_of_rings, Explicit_Filter_Constants::HASELBACHER_FILTER);
    
    Vector3D Delta = theNeighbours.Delta;
    Vector3D kmax(PI/Delta.x,PI/Delta.y,PI/Delta.z);
    
    
    DenseMatrix A = Matrix_A(theCell, theNeighbours, commutation_order);
    
    
    for (int k=0; k<N_weight_factor; k++) {
        double weight_factor = min_weight_factor + k*(max_weight_factor-min_weight_factor)/(N_weight_factor-1.0);
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
    }
    
    
    
    
    
    char prefix[256], suffix[256], extension[256], output_file_name[256];
    char *output_file_name_ptr;
    ofstream output_file;    
    
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
    << "\"FGR\" \\ \n"
    << "\"Weight factor\" \\ \n"
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
        output_file         << " " << FGR[k]
        << " " << weight_factors[k] 
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

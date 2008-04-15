/*
 *  LES_Filters.h
 *  CFFC
 *
 *  Created by Willem Deconinck on 25/03/08.
 *
 */

#ifndef _LES_FILTERS_INCLUDED_
#define _LES_FILTERS_INCLUDED_

/*
 *
 *  Several speedups possible:
 *
 *  LES_Filter:
 *  -----------
 *   --> instead of copying an entire Soln_Blk, just the Soln_States
 *
 *  Haselbacher_Filter:
 *  -------------------
 *   --> Case Grid = uniform:
 *        Matrix_A and Matrix_W is same for every cell
 *
 *
 */



#include "../LES_Polytropic/LES3DPolytropic.h"
#include "../../Math/math.h"
#include "../../Math/LinearSystems.h"
#include "../../Math/Matrix.h"
#include "../../TurbulenceModelling/TurbulenceModelling.h"
#include "../../TurbulenceModelling/TurbulenceAveraging.h"
#include <complex>
#define _GNUPLOT
#include "../../System/gnuplot.h"
#include "../../Utilities/Utilities.h"

#define LES_FILTER_HASELBACHER 1
#define MAX_NUMBER_OF_NEIGHBOURS_IN_LES_FILTER 100



// Complex type definition 
typedef  complex<double>  Complex;

/* ------------------------------------------------------------------------------------------------------------------------------ */
/**
 * CLASS : Neighbours
 */
class Neighbours{
public:
    Cell3D *neighbour;
    Cell3D *theCell;
    Grid3D_Hexa_Block *Grid_ptr;
    
    
    int number_of_neighbours;
    Neighbours(Grid3D_Hexa_Block &Grid){
        Grid_ptr = &Grid;
        neighbour = new Cell3D [MAX_NUMBER_OF_NEIGHBOURS_IN_LES_FILTER];
    }
    ~Neighbours(void) {
        delete[] neighbour;
    }
    void GetNeighbours(Cell3D &theCell);
    
    friend ostream &operator << (ostream &out_file, const Neighbours &neighbours);
//    friend istream &operator >> (istream &in_file, Neighbours &neighbours);
    
};


/**
 * Neighbours -- Input-output operators. 
 */
inline ostream &operator << (ostream &out_file,
                             const Neighbours &neighbours) {
    out_file << " " << neighbours.number_of_neighbours << "\n";
    for (int i=0; i<neighbours.number_of_neighbours; i++) {
        out_file << " " << neighbours.neighbour[i] << "\n";
    }
    return (out_file);
}

/* ------------------------------------------------------------------------------------------------------------------------------ */
/**
 * CLASS: General Filter
 * Abstract class
 */
template <typename Soln_pState, typename Soln_cState>
class General_Filter {
public:
    virtual Soln_pState filter(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell) = 0;
    virtual void transfer_function(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell, int kmax) =0;

};

template <typename Soln_pState, typename Soln_cState>
class Haselbacher_Filter;


/* ------------------------------------------------------------------------------------------------------------------------------ */
/**
 * CLASS : LES_Filter
 */
template<typename Soln_pState, typename Soln_cState>
class LES_Filter {
public:
    
    AdaptiveBlock3D_List *LocalSolnBlkList_ptr;                     // List with properties of SolnBlks
    Hexa_Block<Soln_pState,Soln_cState> *Solution_Blocks_ptr;       // array of SolnBlks
    /* ----- constructor ----- */
    LES_Filter(HexaSolver_Data &Data,
               HexaSolver_Solution_Data<Soln_pState,Soln_cState> &Solution_Data,
               int filter_flag) {
        
        Solution_Blocks_ptr  = Solution_Data.Local_Solution_Blocks.Soln_Blks;
        LocalSolnBlkList_ptr = &(Data.Local_Adaptive_Block_List);
        
        if (filter_flag == LES_FILTER_HASELBACHER) {
            filter_ptr = new Haselbacher_Filter<Soln_pState,Soln_cState>;            
        }
    }

    General_Filter<Soln_pState,Soln_cState> *filter_ptr;
    void filter();
    
    double filter_width;
    void transfer_function();
    
    double gaussian();
    double tophat();

    double maximum_wavenumber();
};


template<typename Soln_pState, typename Soln_cState>
double LES_Filter<Soln_pState,Soln_cState>::maximum_wavenumber() {
    double max_cell_volume = Max_Cell_Volume(Solution_Blocks_ptr,*LocalSolnBlkList_ptr);
    double kmax = PI/pow(max_cell_volume,ONE/THREE);
    return kmax;
}


template<typename Soln_pState, typename Soln_cState>
void LES_Filter<Soln_pState,Soln_cState>::transfer_function() {

    if (LocalSolnBlkList_ptr->Nused() >= 1) {
        for (int nBlk = 0; nBlk <= LocalSolnBlkList_ptr->Nused(); ++nBlk ) {
            if (LocalSolnBlkList_ptr->Block[nBlk].used == ADAPTIVEBLOCK3D_USED) {
            
                filter_ptr->transfer_function(Solution_Blocks_ptr[nBlk],Solution_Blocks_ptr[nBlk].Grid.Cell[6][6][6],maximum_wavenumber());
            
                return;
            }
        }
    }

}

template<typename Soln_pState, typename Soln_cState>
void LES_Filter<Soln_pState,Soln_cState>::filter() {
    /* For every local solution block */
    if (LocalSolnBlkList_ptr->Nused() >= 1) {
        Hexa_Block<Soln_pState,Soln_cState> Filtered_Solution_Block;
        for (int nBlk = 0; nBlk <= LocalSolnBlkList_ptr->Nused(); ++nBlk ) {
            if (LocalSolnBlkList_ptr->Block[nBlk].used == ADAPTIVEBLOCK3D_USED) {
                /* For every cell */
                Filtered_Solution_Block.Copy(Solution_Blocks_ptr[nBlk]);    // could be more efficient
                for(int i=Solution_Blocks_ptr[nBlk].Grid.ICl ; i<=Solution_Blocks_ptr[nBlk].Grid.ICu ; i++) {
                    for (int j=Solution_Blocks_ptr[nBlk].Grid.JCl ; j<=Solution_Blocks_ptr[nBlk].Grid.JCu ; j++) {
                        for (int k=Solution_Blocks_ptr[nBlk].Grid.KCl ; k<=Solution_Blocks_ptr[nBlk].Grid.KCu ; k++) {
                            Filtered_Solution_Block.W[i][j][k] = 
                                filter_ptr->filter(Solution_Blocks_ptr[nBlk],Solution_Blocks_ptr[nBlk].Grid.Cell[i][j][k]);
                            Filtered_Solution_Block.U[i][j][k] = Filtered_Solution_Block.W[i][j][k].U();
                        }
                    }
                }
                Solution_Blocks_ptr[nBlk].Copy(Filtered_Solution_Block);   // could be more efficient
            }         
        } 
    } 
}


/* ------------------------------------------------------------------------------------------------------------------------------- */
/**
 * CLASS: Haselbacher_Filter
 */
template <typename Soln_pState, typename Soln_cState>
class Haselbacher_Filter : public General_Filter<Soln_pState,Soln_cState> {
public:

    
    Haselbacher_Filter(void) {
        commutation_order = 2;
        filter_width = 0.1;
        weight_factor = 3.1;
        the_number_of_unknowns = number_of_unknowns();
        
        relaxation_factor = 0.47;

    }

    double relaxation_factor;
    double filter_width;
    int commutation_order;
    int relaxation_flag;
    int weight_flag;
    double weight_factor;
    int number_of_unknowns();
    int the_number_of_unknowns;
    int number_of_terms_of_degree(int order);
    int fac(int n);
    double trinomial_coefficient(int n1, int n2, int n3);
    
    Soln_pState filter(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell);
    
    void transfer_function(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell, int kmax);

    Soln_pState LeastSquaresReconstruction(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell, Neighbours &theNeighbours);
    DenseMatrix Matrix_A(Cell3D &theCell, Neighbours &theNeighbours);
    DenseMatrix Matrix_b(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Neighbours &theNeighbours);
    DiagonalMatrix Matrix_W(Cell3D &theCell, Neighbours &theNeighbours);
    DiagonalMatrix Matrix_W(Cell3D &theCell, Neighbours &theNeighbours, double weight_factor);

    
    void Weight_Matrix_A_and_b(DenseMatrix &A, DenseMatrix &b, Cell3D &theCell, Neighbours &theNeighbours);
    
    
};

template<typename Soln_pState, typename Soln_cState>
void Haselbacher_Filter<Soln_pState,Soln_cState>::transfer_function(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell, int kmax) {
    
    Neighbours theNeighbours(SolnBlk.Grid);
    theNeighbours.GetNeighbours(theCell);
    
    DenseMatrix A = Matrix_A(theCell, theNeighbours);

    
    double w0;
    
    cout << "kmax = " << kmax << endl;
    
    
    /* --------------- allocations ----------------- */
    int N=50;
    Cell3D  ***kCells   = new Cell3D  **[N];
    Complex ***G        = new Complex **[N];
    for (int i=0; i<N; i++) {
        kCells[i]   = new Cell3D  *[N];
        G[i]        = new Complex *[N];
        for (int j=0; j<N; j++) {
            kCells[i][j]    = new Cell3D  [N];
            G[i][j]         = new Complex [N];
        }
    }
    
    /* --------------- assign kCells ------------------- */
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
            for (int k=0; k<N; k++) {
                kCells[i][j][k].setindex(i,j,k);
                kCells[i][j][k].setloc(i*kmax/(N-ONE),j*kmax/(N-ONE),k*kmax/(N-ONE));
            }
        }
    }

    
    
    K_BinaryTree K_tree;
    K_container *K = new K_container [N*N*N];  // Container of information for a wavenumber
    int index;
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
            for (int k=0; k<N; k++) {
                index = k+N*(j+N*i);
                K[index].k = kCells[i][j][k].Xc.abs();
                K[index].index = index;
                K_tree.InsertNode(K[index]);
            }
        }
    }
    delete[] K;
    K = K_tree.asArray();
    int nK = K_tree.countNodes();
    
    

    

    
    
    double x,y,z,x0,y0,z0;
    x0 = theCell.Xc.x;
    y0 = theCell.Xc.y;
    z0 = theCell.Xc.z;
    Complex G_max = G[N-1][N-1][N-1];
    double a=0.0, b=1.0;
    int counter;
    double *kk = new double [nK];
    double *GR = new double [nK];
    double *GRs = new double [nK];
    Gnuplot_Control h1;
    h1.gnuplot_init();
    h1.gnuplot_setstyle("lines") ;
    h1.gnuplot_cmd("set grid");
    h1.gnuplot_cmd("set xrange [0:50]");
    h1.gnuplot_cmd("set yrange [-1:1]");
    h1.gnuplot_set_xlabel("k");
    h1.gnuplot_set_ylabel("G(k)");
    h1.gnuplot_set_title("Transfer function");
    
    for (double weight_factor=0.1; weight_factor<=6.0; weight_factor+=0.5) {
        
        DiagonalMatrix W = Matrix_W(theCell, theNeighbours,weight_factor);
        DenseMatrix Z = (W*A).pseudo_inverse()*W;
        RowVector w = Z[0];
        
        for (int i=0; i<N; i++) {
            for (int j=0; j<N; j++) {
                for (int k=0; k<N; k++) {
                    G[i][j][k] = 0;
                    for(int n=0; n<theNeighbours.number_of_neighbours; n++) {
                        x=theNeighbours.neighbour[n].Xc.x;
                        y=theNeighbours.neighbour[n].Xc.y;
                        z=theNeighbours.neighbour[n].Xc.z;
                        G[i][j][k] += w(n) * exp(I*(kCells[i][j][k].Xc.x*(x0-x) + kCells[i][j][k].Xc.y*(y0-y) + kCells[i][j][k].Xc.z*(z0-z)));
                    }
                }
            }
        }
        
        
        Complex G_max = G[N-1][N-1][N-1];
        a=0.0, b=1.0;
        counter = 0;
        w0 = (a+b)/TWO;
        while( fabs(real(w0+(ONE-w0)*G_max)) >= 0.01 ) {
            if( real(a+(ONE-a)*G_max) * real(w0+(ONE-w0)*G_max) <= ZERO)
                b = w0;
            else
                a = w0;
            w0 = (a+b)/TWO;
            counter++;
            if(counter >=1000) {
                cout << "max reached" << endl;
                break;
            }
        }
        cout << "w0 = " << w0 << endl;
        
        
        
        for (int i=0; i<N; i++) {
            for (int j=0; j<N; j++) {
                for (int k=0; k<N; k++) {
                    G[i][j][k] = w0 + (ONE - w0) * G[i][j][k];
                }
            }
        }
        
        
        
        
        int i,j,k;
        for (int ii=0; ii<nK; ii++) {
            K[ii].Ek = 0;
            K[ii].Ek_smooth = 0;
            for(int jj=0; jj<K[ii].N; jj++) {       // For all grid points with this value of abs_k
                index = K[ii].indexes[jj];      // index corresponding to (i,j,l)
                i = index/(N*N);
                j = (index - N*N*i)/N;
                k = index - N*(j+N*i);
                //            Print_3(i,j,k);
                K[ii].Ek += real(G[i][j][k])/double(K[ii].N);
                K[ii].Ek_smooth += imag(G[i][j][k])/double(K[ii].N);
            }
            kk[ii]=K[ii].k;
            GR[ii]=K[ii].Ek;
        }
        polyfit_smoothing(nK, kk, GR, 3, 2.5, true, GRs)  ;  
        /* ---------- output ----------- */
        
        char title[80];
        sprintf(title,"weight = %3.2f  w0 = %3.2f",weight_factor,w0);
        
        h1.gnuplot_plot1d_var2(kk,GRs,nK,title);
        
    }
    
//    dpoint *dpr = new dpoint[nK];
//    dpoint *dpi = new dpoint[nK];
//
//    for (int ii=0; ii<nK; ii++) {
//        dpr[ii].x = K[ii].k;
//        dpr[ii].y = K[ii].Ek;
//        dpi[ii].x = K[ii].k;
//        dpi[ii].y = K[ii].Ek_smooth;
//    }


    
    /* ----------- deallocations ------------ */
    for (int ii=0; ii<nK; ii++) {
        delete[] K[ii].indexes;
    }
    delete[] K;

    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
            delete[] kCells[i][j];
            delete[] G[i][j];
        }
        delete[] kCells[i];
        delete[] G[i];
    }
    delete[] kCells;
    delete[] G;
    //delete[] dpr;
    //delete[] dpi;

    delete[] kk;
    delete[] GR;
    delete[] GRs;
    
    return; 
}

template<typename Soln_pState, typename Soln_cState>
Soln_pState Haselbacher_Filter<Soln_pState,Soln_cState>::filter(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell) {

    Neighbours theNeighbours(SolnBlk.Grid);
    theNeighbours.GetNeighbours(theCell);
    
    return LeastSquaresReconstruction(SolnBlk, theCell, theNeighbours);
}



template<typename Soln_pState, typename Soln_cState>
DenseMatrix Haselbacher_Filter<Soln_pState,Soln_cState>::Matrix_A(Cell3D &theCell, Neighbours &theNeighbours) {

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
                            //A(i,j) = coefficient * pow(1.0,n1)*pow(1.0,n2)*pow(1.0,n3);

                        }
                        j++; // go to next column
                        
                    }
                }
            }
        } 
        
    }
    
    bool output_this = false;
    if (output_this == true) {
        cout << endl;
        cout << "neighbours = " << theNeighbours << endl;
        cout << "the_number_of_unknowns = " << the_number_of_unknowns << endl;
        cout << "the_number_of_neighbours = " << the_number_of_neighbours << endl;
        cout << "dx = " << dx << endl;
        cout << "dy = " << dy << endl;
        cout << "dz = " << dz << endl;
        cout << "matrix A = " << endl << A << endl;   
    }
    
    return A;
}



template<typename Soln_pState, typename Soln_cState>
DiagonalMatrix Haselbacher_Filter<Soln_pState,Soln_cState>::Matrix_W(Cell3D &theCell, Neighbours &theNeighbours, double weight_factor) {
    
    int the_number_of_neighbours = theNeighbours.number_of_neighbours;
    DiagonalMatrix W(the_number_of_neighbours);
    
    
    double Delta;
    double r;
    
    Delta = pow(theCell.V,1.0/3.0) * weight_factor;
    
    for (int i=0; i<the_number_of_neighbours; i++) {
        
        r=(theNeighbours.neighbour[i].Xc - theCell.Xc).abs();
        W(i) = sqrt(6.0/(PI*pow(Delta,2)) * exp(- 6.0* pow(r/Delta,2))) ;
        
        //        cout << "weight_factor = " << weight_factor << "   Delta = " << Delta << "    r = " << r << "    W(i) = " << W(i) << endl;
        
    }
    
    return W;
}

template<typename Soln_pState, typename Soln_cState>
DiagonalMatrix Haselbacher_Filter<Soln_pState,Soln_cState>::Matrix_W(Cell3D &theCell, Neighbours &theNeighbours) {
    return Matrix_W(theCell, theNeighbours, weight_factor);
}



template<typename Soln_pState, typename Soln_cState>
DenseMatrix Haselbacher_Filter<Soln_pState,Soln_cState>::Matrix_b(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Neighbours &theNeighbours) {
    
    int the_number_of_neighbours = theNeighbours.number_of_neighbours;
    DenseMatrix b(the_number_of_neighbours,SolnBlk.NumVar());
    int I,J,K;
    for (int i=0; i<the_number_of_neighbours; i++) {
        I = theNeighbours.neighbour[i].I;
        J = theNeighbours.neighbour[i].J;
        K = theNeighbours.neighbour[i].K;
        for(int j=1; j<=SolnBlk.NumVar(); j++)
            b(i,j-1) = SolnBlk.W[I][J][K][j];
    }
    return b;
}


//template<typename Soln_pState, typename Soln_cState>
//DenseMatrix Haselbacher_Filter<Soln_pState,Soln_cState>::Relaxation_Matrix_x(DenseMatrix &x, Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Neighbours &theNeighbours) {
//    
//    double w0 = relaxation_factor;
//    
//    w0 + (ONE - w0)
//    int the_number_of_neighbours = theNeighbours.number_of_neighbours;
//    DenseMatrix b(the_number_of_neighbours,SolnBlk.NumVar());
//    int I,J,K;
//    for (int i=0; i<the_number_of_neighbours; i++) {
//        I = theNeighbours.neighbour[i].I;
//        J = theNeighbours.neighbour[i].J;
//        K = theNeighbours.neighbour[i].K;
//        for(int j=1; j<=SolnBlk.NumVar(); j++)
//            b(i,j-1) = SolnBlk.W[I][J][K][j];
//    }
//    return b;
//}

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
Soln_pState Haselbacher_Filter<Soln_pState,Soln_cState>::LeastSquaresReconstruction(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell, Neighbours &theNeighbours) {

    /* ---------------- Matrix A ---------------- */
    DenseMatrix A = Matrix_A(theCell,theNeighbours);
    
    /* ---------------- Matrix b ---------------- */
    DenseMatrix b = Matrix_b(SolnBlk,theNeighbours);
        
    Weight_Matrix_A_and_b(A, b, theCell, theNeighbours);
    
    /* ------------ LS procedure for Unknowns DenseMatrix x ---------- */
    DenseMatrix x(the_number_of_unknowns,SolnBlk.NumVar());
    int krank;
    ColumnVector Rnorm(the_number_of_unknowns);
    Solve_LS_Householder(A,b,x,krank,Rnorm);
    
    
    /* --------------- relaxation -------------- */
    RowVector xrow = x[0];
    RowVector x0(SolnBlk.NumVar());
    int I,J,K;
    I = theCell.I;
    J = theCell.J;
    K = theCell.K;
    for (int n=1; n<=SolnBlk.NumVar(); n++) {
        x0(n-1) = SolnBlk.W[I][J][K][n];
    }
    double w0 = relaxation_factor;

    xrow = x0*w0 + xrow*(ONE-w0);
    //xrow += w0*x0;
    
    /* ---------- Filtered value is first element of x ----------- */
    Soln_pState temp;
    for (int j=1; j<SolnBlk.NumVar(); j++)
        temp[j] = xrow(j-1);
    
    
    /* NOTE: MAKE "b" AND "x" INTO A "ColumnVector" FOR Single VARIABLES
    
//    /* --------------- Vector b --------------- */
//    ColumnVector b(the_number_of_neighbours);
//    int I,J,K;
//    for (int i=0; i<the_number_of_neighbours; i++) {
//        I = theNeighbours.neighbour[i].I;
//        J = theNeighbours.neighbour[i].J;
//        K = theNeighbours.neighbour[i].K;
//        b(i) = SolnBlk.W[I][J][K].v.x;
//    }
//
//    /* ------------ LS procedure for Unknowns Vector x ---------- */
//    ColumnVector x(the_number_of_unknowns);
//    int krank;
//    double Rnorm;
//    Solve_LS_Householder(A,b,x,krank,Rnorm);
//
//    /* ---------- Filtered value is first element of x ----------- */
//    return x(0); 
    
    return temp;  

}



/* ------------ calculates the number of unknowns for the given order ----------- */
template <typename Soln_pState, typename Soln_cState>
int Haselbacher_Filter<Soln_pState,Soln_cState>::number_of_unknowns() {
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





/* ----------------------------------------------------------------------------------------------------------------------- 
typedef double (theClass::*function_with_one_argument) (const double &abs_wave_num) const;
_Member_Function_Wrapper_<theClass,function_with_one_argument, double> mapped_function(this, &theClass::Energy_Spectrum_Value);
double dummy;
double TKE_entire_range = AdaptiveGaussianQuadrature(mapped_function, 0.0, k_eta/2.0, dummy,5);
/**/



template<typename Soln_pState, typename Soln_cState>
inline double return_it(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, double (func)(Soln_pState &)){
    /* For static member functions or any function */
    double b;
    b = func(SolnBlk.W[2][2][2]);
    
    return b;
}

template<typename Soln_pState, typename Soln_cState>
inline double return_it_2(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, double (Soln_pState::*&func)(void)){
    /* for non static member functions */
    double b;
    b = (SolnBlk.W[2][2][2].*func)();
    
    return b;
}

template<typename Soln_pState, typename Soln_cState>
inline double return_it_2(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, double (Soln_pState::*&func)){
    /* for non static member functions */
    double b;
    b = (SolnBlk.W[2][2][2].*func);
    
    return b;
}



template<typename Soln_pState, typename Soln_cState, typename func_object>
inline double return_it_3(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, func_object func) {
    /* for non static member functions */
    double b;
    b = (SolnBlk.W[2][2][2].*func);
    
    return b;
}




#endif

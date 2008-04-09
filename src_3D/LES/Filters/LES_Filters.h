/*
 *  LES_Filters.h
 *  CFFC
 *
 *  Created by Willem Deconinck on 25/03/08.
 *
 */

#include "../../Math/math.h"
#include "../../Math/LinearSystems.h"
#include "../../Math/Matrix.h"
#include "../LES_Polytropic/LES3DPolytropic.h"




template<typename Soln_pState, typename Soln_cState>
class Haselbacher_Filter;

template<typename Soln_pState, typename Soln_cState>
class LES_Filter {
public:


    
    
    double filter_width;
    double transfer_function(/*node*/);
    Haselbacher_Filter<Soln_pState,Soln_cState> haselbacher;
    double gaussian();
    double tophat();


};

template <typename Soln_pState, typename Soln_cState>
class Haselbacher_Filter {
public:
    AdaptiveBlock3D_List *LocalSolnBlkList_ptr;  // List with properties of SolnBlks
    Hexa_Block<Soln_pState,Soln_cState> *Solution_Blocks_ptr;              // array of SolnBlks
    /* ----- constructor ----- */
    Haselbacher_Filter(HexaSolver_Data &Data,
                       HexaSolver_Solution_Data<Soln_pState,Soln_cState> &Solution_Data) {

        Solution_Blocks_ptr = &(Solution_Data.Local_Solution_Blocks.Soln_Blks);
        LocalSolnBlkList_ptr = &(Data.Local_Adaptive_Block_List);
        
    }
    
//    Solution_Data.Local_Solution_Blocks.Soln_Blks,
//    Data.Local_Adaptive_Block_List
    

    
    double filter_width;
    int commutation_order;
    int relaxation_flag;
    int weight_flag;
    double weight_factor;
    void filter();
    int number_of_unknowns();
    int number_of_terms_of_degree(int order);
    int fac(int n);
    int trinomial_coefficient(int n1, int n2, int n3);
    double LeastSquaresReconstruction(Cell3D Xo, Cell3D *neighbours);
};

template<typename Soln_pState, typename Soln_cState>
void Haselbacher_Filter<Soln_pState,Soln_cState>::filter() {
    /* ----- number of unknowns ----- */
    int n = number_of_unknowns(commutation_order);   
    cout << "number of unknowns = " << n << endl;

	int nx,ny;
    double **F, **A, *RHS, *unknowns, *W;
	double ** filtered_signal = DenseMatrix(nx,ny);
    Hexa_Block<Soln_pState,Soln_cState> *SolnBlk;
    
    double *var;

    /* For every local solution block */
    if (LocalSolnBlkList_ptr->Nused() >= 1) {
        for (int nBlk = 0; nBlk <= LocalSolnBlkList_ptr->Nused(); ++nBlk ) {
            if (LocalSolnBlkList_ptr->Block[nBlk].used == ADAPTIVEBLOCK3D_USED) {
                SolnBlk = &(Solution_Blocks_ptr[nBlk]);
                
                /* For every cell */
                for(int i=SolnBlk.ICl ; i<=SolnBlk.ICu ; i++) {
                    for (int j=SolnBlk.JCl ; j<=SolnBlk.JCu ; j++) {
                        for (int k=SolnBlk.KCl ; k<=SolnBlk.KCu ; k++) {
                            
                            //neighbours = Get_Neighbours(SolnBlk.Grid,SolnBlk.Grid.Cell[i][j][k]);
                            
                            //v_filt.x = Least_Squares_Reconstruction(SolnBlk,SolnBlk.Grid.Cell[i][j][k]);
                            

                        }
                    }
                }
                
            }         
        } /* endfor */
    } /* endif */
        
    
    
    
}


/* ------------ calculates the number of unknowns for the given order ----------- */
template<typename Soln_pState, typename Soln_cState>
int Haselbacher_Filter<Soln_pState,Soln_cState>::number_of_unknowns() {
    int number = 0;
    for (int k=0; k<= commutation_order; k++) {
        number += number_of_terms_of_degree(k);
    }
    return number;
}

/* ----- faculty ----- */
template<typename Soln_pState, typename Soln_cState>
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
template<typename Soln_pState, typename Soln_cState>
int Haselbacher_Filter<Soln_pState,Soln_cState>::number_of_terms_of_degree(int k){
    int n=3; // 3 dimensions
    return (  fac(n+k-1)/(fac(k)*fac(n-1))  );
}




/* ------------------------ trinomial coefficient -------------------------- */
/* 
 * coefficient of x^n1 y^n2 z^n3 in expansion of (x+y+z)^n (with n=n1+n2+n3)
 *  (n1 + n2 + n3)!             (               (n1 + n2 + n3 + ... + nk)!   )
 *  ---------------             ( multinomial:  ---------------------------  )
 *  n1! + n2! + n3!             (               n1! + n2! + n3! + ... + nk!  )
 */
template<typename Soln_pState, typename Soln_cState>
int Haselbacher_Filter<Soln_pState,Soln_cState>::trinomial_coefficient(int n1, int n2, int n3){
    return (  fac(n1+n2+n3)/( fac(n1) + fac(n2) + fac(n3) ) );
}

/* 
 * A x = b  Least Squares  and filtered variable is first element of x
 */
template<typename Soln_pState, typename Soln_cState>
double Haselbacher_Filter<Soln_pState,Soln_cState>::LeastSquaresReconstruction(Cell3D theCell, Cell3D *Cells) {
    int number_of_neighbours;
    
    
    /* ---------------- Matrix A ---------------- */
    DenseMatrix A(number_of_neighbours,number_of_unknowns);
    int j=0;
    int index;
    double dx, dy, dz, coefficient;
    for (int k=0; k<=commutation_order; k++) {
        
        // for all terms of degree k
        for (int n3=0; n3<=k; n3++) {
            for (int n2=0; n2<=k; n2++) {
                for (int n1=0; n1<=k; n1++) {
                    if (n1+n2+n3 == k) {
                        
                        // fill this column
                        coefficient = trinomial_coefficient(n1,n2,n3)/fac(k);
                        for (int i=0; i<number_of_neighbours; i++) {
                            //index = neighbours[i];

                            dx = Cells[index].Xc.x - theCell.Xc.x;
                            dy = Cells[index].Xc.y - theCell.Xc.y;
                            dz = Cells[index].Xc.z - theCell.Xc.z;
                            A(i,j) = coefficient * pow(dx,n1)*pow(dy,n2)*pow(dz,n3);
                        }
                        j++; // go to next column
                        
                    }
                }
            }
        } 
        
    }
    
    
    /* NOTE: MAKE "b" AND "x" INTO A "DenseMatrix" FOR MULTIPLE VARIABLES
    
    /* --------------- Vector b --------------- */
    ColumnVector b(number_of_neighbours);
    for (int i=0; i<number_of_neighbours; i++) {
        b(i) = Solution_Blocks_ptr.W[index].v.x;
    }
    
    /* ------------ LS procedure for Unknowns Vector x ---------- */
    ColumnVector x(number_of_unknowns);
    int krank;
    double Rnorm;
    Solve_LS_Householder(A,b,x,krank,Rnorm);
    
    /* ---------- Filtered value is first element of x ----------- */
    return x(0);  
}


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

#define LES_FILTER_HASELBACHER 1


// Abstract class
template <typename Soln_pState, typename Soln_cState>
class General_Filter {
public:
    virtual void filter(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell) = 0;

};

template <typename Soln_pState, typename Soln_cState>
class Haselbacher_Filter;

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
    double transfer_function(/*node*/);
    
    double gaussian();
    double tophat();


};


template<typename Soln_pState, typename Soln_cState>
void LES_Filter<Soln_pState,Soln_cState>::filter() {
    
    /* For every local solution block */
    if (LocalSolnBlkList_ptr->Nused() >= 1) {
        for (int nBlk = 0; nBlk <= LocalSolnBlkList_ptr->Nused(); ++nBlk ) {
            if (LocalSolnBlkList_ptr->Block[nBlk].used == ADAPTIVEBLOCK3D_USED) {
                
                /* For every cell */
                for(int i=Solution_Blocks_ptr[nBlk].ICl ; i<=Solution_Blocks_ptr[nBlk].ICu ; i++) {
                    for (int j=Solution_Blocks_ptr[nBlk].JCl ; j<=Solution_Blocks_ptr[nBlk].JCu ; j++) {
                        for (int k=Solution_Blocks_ptr[nBlk].KCl ; k<=Solution_Blocks_ptr[nBlk].KCu ; k++) {
                            
                            filter_ptr->filter(Solution_Blocks_ptr[nBlk],Solution_Blocks_ptr[nBlk].Grid.Cell[i][j][k]);
                           
                        }
                    }
                }
                
            }         
        } /* endfor */
    } /* endif */
}

template <typename Soln_pState, typename Soln_cState>
class Haselbacher_Filter : public General_Filter<Soln_pState,Soln_cState> {
public:

    
//    Solution_Data.Local_Solution_Blocks.Soln_Blks,
//    Data.Local_Adaptive_Block_List
    
    Haselbacher_Filter(void) {
        commutation_order = 2;
        filter_width = 0.1;
    }

    
    double filter_width;
    int commutation_order;
    int relaxation_flag;
    int weight_flag;
    double weight_factor;
    int number_of_unknowns();
    int number_of_terms_of_degree(int order);
    int fac(int n);
    int trinomial_coefficient(int n1, int n2, int n3);
    
    void filter(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell);

    double LeastSquaresReconstruction(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell, Cell3D *Cells);
};

template<typename Soln_pState, typename Soln_cState>
inline void Haselbacher_Filter<Soln_pState,Soln_cState>::filter(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell) {

    cout << "filtering cell " << theCell << endl;

    // Get neighbours
    
    
    //double v_filt_x = LeastSquaresReconstruction(SolnBlk, theCell, neighbouringCells);
    
    
}


/* 
 * A x = b  Least Squares  and filtered variable is first element of x
 */
template<typename Soln_pState, typename Soln_cState>
double Haselbacher_Filter<Soln_pState,Soln_cState>::LeastSquaresReconstruction(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell, Cell3D *Cells) {
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
        b(i) = SolnBlk.W[index].v.x;
    }
    
    /* ------------ LS procedure for Unknowns Vector x ---------- */
    ColumnVector x(number_of_unknowns);
    int krank;
    double Rnorm;
    Solve_LS_Householder(A,b,x,krank,Rnorm);
    
    /* ---------- Filtered value is first element of x ----------- */
    return x(0);  
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


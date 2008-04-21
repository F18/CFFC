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
#define MAX_NUMBER_OF_NEIGHBOURING_RINGS_IN_LES_FILTER 6


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
    
    static int max_number_of_neighbours;
    static int number_of_rings;
    bool Allocated;
    int number_of_neighbours;
    Vector3D Delta;
    Neighbours(Grid3D_Hexa_Block &Grid){
        Grid_ptr = &Grid;
        Allocated = true;
        neighbour = new Cell3D [points(MAX_NUMBER_OF_NEIGHBOURING_RINGS_IN_LES_FILTER)];
    }
    ~Neighbours(void) {
        delete[] neighbour;
    }
    void GetNeighbours(Cell3D &theCell, int number_of_rings);
    
    int pointsPerEdge(int ring);
    int pointsPerFace(int ring);
    int points(int number_of_rings);
    
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
    virtual void transfer_function(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell, double &kmax) =0;

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

    ~LES_Filter() {
        delete filter_ptr;
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
    double kmax = maximum_wavenumber();
    if (LocalSolnBlkList_ptr->Nused() >= 1) {
        for (int nBlk = 0; nBlk <= LocalSolnBlkList_ptr->Nused(); ++nBlk ) {
            if (LocalSolnBlkList_ptr->Block[nBlk].used == ADAPTIVEBLOCK3D_USED) {
            
                filter_ptr->transfer_function(Solution_Blocks_ptr[nBlk],Solution_Blocks_ptr[nBlk].Grid.Cell[8][8][8],kmax);
            
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
        for (int nBlk = 0; nBlk < LocalSolnBlkList_ptr->Nused(); nBlk++ ) {
            if (LocalSolnBlkList_ptr->Block[nBlk].used == ADAPTIVEBLOCK3D_USED) {
                /* For every cell */
                Filtered_Solution_Block.Copy(Solution_Blocks_ptr[nBlk]);    // could be more efficient
                for(int i=Solution_Blocks_ptr[nBlk].Grid.ICl ; i<=Solution_Blocks_ptr[nBlk].Grid.ICu ; i++) {
                    for (int j=Solution_Blocks_ptr[nBlk].Grid.JCl ; j<=Solution_Blocks_ptr[nBlk].Grid.JCu ; j++) {
                        for (int k=Solution_Blocks_ptr[nBlk].Grid.KCl ; k<=Solution_Blocks_ptr[nBlk].Grid.KCu ; k++) {
                            Filtered_Solution_Block.W[i][j][k] = 
                                filter_ptr->filter(Solution_Blocks_ptr[nBlk],Solution_Blocks_ptr[nBlk].Grid.Cell[i][j][k]);                            
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
        Assigned_w_Uniform_Grid = false;
        UNIFORM_GRID = true;
        number_of_rings = 2; //2;
        commutation_order = 3; // 2;
        filter_width = 0.1;
        weight_factor = 2.5;
        the_number_of_unknowns = number_of_unknowns();
        
        relaxation_factor = 0.08;

    }

    int number_of_rings;
    double relaxation_factor;
    double filter_width;
    int commutation_order;
    int relaxation_flag;
    int weight_flag;
    double weight_factor;
    RowVector w_Uniform_Grid;
    bool Assigned_w_Uniform_Grid;
    bool UNIFORM_GRID;
    void Assign_w_Uniform_Grid(Cell3D &theCell, Neighbours &theNeighbours);
    int number_of_unknowns();
    int number_of_unknowns(int &commutation_order);
    int the_number_of_unknowns;
    int number_of_terms_of_degree(int order);
    int fac(int n);
    double trinomial_coefficient(int n1, int n2, int n3);
    
    Soln_pState filter(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell);
    
    void transfer_function(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell, double &kmax);

    Soln_pState LeastSquaresReconstruction(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell, Neighbours &theNeighbours);
    Soln_pState LeastSquaresReconstruction_Uniform_Grid(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell, Neighbours &theNeighbours);
    DenseMatrix Matrix_A(Cell3D &theCell, Neighbours &theNeighbours);
    DenseMatrix Matrix_A(Cell3D &theCell, Neighbours &theNeighbours, int &commutation_order);
    DenseMatrix Matrix_b(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Neighbours &theNeighbours);
    DiagonalMatrix Matrix_W(Cell3D &theCell, Neighbours &theNeighbours);
    DiagonalMatrix Matrix_W(Cell3D &theCell, Neighbours &theNeighbours, double weight_factor);

    
    void Weight_Matrix_A_and_b(DenseMatrix &A, DenseMatrix &b, Cell3D &theCell, Neighbours &theNeighbours);
    
    Complex G_function(Cell3D &theCell, Neighbours &theNeighbours, Cell3D &kCell, RowVector &w);
    Complex G_function(Cell3D &theCell, Neighbours &theNeighbours, Cell3D &kCell, RowVector &w, double &w0);
    Complex G_function_1D(Cell3D &theCell, Neighbours &theNeighbours, double &wave_number, RowVector &w);
    Complex G_function_1D(Cell3D &theCell, Neighbours &theNeighbours, double &wave_number, RowVector &w, double &w0);
    Complex G_function_1D_100(Cell3D &theCell, Neighbours &theNeighbours, double &wave_number, RowVector &w);
    Complex G_function_1D_100(Cell3D &theCell, Neighbours &theNeighbours, double &wave_number, RowVector &w, double &w0);
    Complex G_function_1D_110(Cell3D &theCell, Neighbours &theNeighbours, double &wave_number, RowVector &w);
    Complex G_function_1D_110(Cell3D &theCell, Neighbours &theNeighbours, double &wave_number, RowVector &w, double &w0);

    Complex dG_function(Cell3D &theCell, Neighbours &theNeighbours, Cell3D &kCell, RowVector &w);
    Complex dG_function(Cell3D &theCell, Neighbours &theNeighbours, Cell3D &kCell, RowVector &w, double &w0);
    Complex dG_function_1D(Cell3D &theCell, Neighbours &theNeighbours, double &wave_number, RowVector &w);
    Complex dG_function_1D(Cell3D &theCell, Neighbours &theNeighbours, double &wave_number, RowVector &w, double &w0);

    double Calculate_relaxation_factor(Cell3D &theCell, Neighbours &theNeighbours, double &maximum_wave_number, RowVector &w);
    double Calculate_wavenumber_of_Gvalue(Cell3D &theCell, Neighbours &theNeighbours, double &kmax, RowVector &w, double &w0, double G_value);
    double Calculate_wavenumber_of_dGvalue(Cell3D &theCell, Neighbours &theNeighbours, double kmin, double &kmax, RowVector &w, double &w0, double G_value);
    double filter_quality(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell, double &kmax, int number_of_rings, int commutation_order, double weight_factor);
    double filter_quality(Cell3D &theCell, Neighbours &theNeighbours, double &kmax, RowVector &w, double &w0);
    double filter_sharpness(Cell3D &theCell, Neighbours &theNeighbours, double &kmax, RowVector &w, double &w0);
    double filter_smoothness(Cell3D &theCell, Neighbours &theNeighbours, double &kmax, RowVector &w, double &w0, double max_smoothness);
    double filter_uniformity(Cell3D &theCell, Neighbours &theNeighbours, double &kmax, RowVector &w, double &w0, double max_distortion);
    int Output_Tecplot(Cell3D  ***kCells, Complex ***G, int N);
    int Output_Filter_types(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell, double kmax);


};


template<typename Soln_pState, typename Soln_cState>
int Haselbacher_Filter<Soln_pState,Soln_cState>::
Output_Tecplot(Cell3D  ***kCells, Complex ***G, int N) {
    
    int i, i_output_title;
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
    sprintf(prefix,"transfer_function");
    
    /* Determine output data file name for this processor. */
    
    sprintf(extension, ".dat");
    strcpy(output_file_name, prefix);
    strcat(output_file_name, extension);
    output_file_name_ptr = output_file_name;
    
    /* Open the output data file. */
    
    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.fail()) return (1);
    
    /* Write the solution data for each solution block. */
    
    output_file << setprecision(14);
    
        output_file << "TITLE = \"" << "Transfer function, "
        << "\"" << "\n"
        << "VARIABLES = "
        << "\"kx\" \\ \n"
        << "\"ky\" \\ \n"
        << "\"kz\" \\ \n"
        //<< "\"G abs\" \\ \n"
        << "\"G real\" \\ \n"
        << "\"G imag\" \\ \n";
        
        output_file<< "ZONE T =  \"Block Number = 0" 
        << "\" \\ \n"
        << "I = " << N << " \\ \n"
        << "J = " << N << " \\ \n"
        << "K = " << N << " \\ \n"
        << "DATAPACKING = POINT \n";
    
    
    for (int k=0; k<N; k++) {
        for (int j=0; j<N; j++) {
            for (int i=0; i<N ; i++) {
                output_file << " " << kCells[i][j][k].Xc 
                            << " " << real(G[i][j][k]) 
                            << " " << imag(G[i][j][k])
                            << "\n";
            } 
        } 
    }
    
    output_file << setprecision(6);
    
    
    
    /* Close the output data file. */
    
    output_file.close();
    
    /* Writing of output data files complete.  Return zero value. */
    
    return(0);
    
}






template<typename Soln_pState, typename Soln_cState>
Complex Haselbacher_Filter<Soln_pState,Soln_cState>::G_function_1D_100(Cell3D &theCell, Neighbours &theNeighbours, double &wave_number, RowVector &w) {
    Cell3D kCell;
    kCell.Xc = Vector3D(wave_number,ZERO,ZERO);    
    return (G_function(theCell, theNeighbours, kCell, w));
}

template<typename Soln_pState, typename Soln_cState>
Complex Haselbacher_Filter<Soln_pState,Soln_cState>::G_function_1D_100(Cell3D &theCell, Neighbours &theNeighbours, double &wave_number, RowVector &w, double &w0) {
    Cell3D kCell;
    kCell.Xc = Vector3D(wave_number,ZERO,ZERO);    
    return (G_function(theCell, theNeighbours, kCell, w, w0));
}

template<typename Soln_pState, typename Soln_cState>
Complex Haselbacher_Filter<Soln_pState,Soln_cState>::G_function_1D_110(Cell3D &theCell, Neighbours &theNeighbours, double &wave_number, RowVector &w) {
    Cell3D kCell;
    kCell.Xc = Vector3D(wave_number/sqrt(TWO),wave_number/sqrt(TWO),ZERO);    
    return (G_function(theCell, theNeighbours, kCell, w));
}

template<typename Soln_pState, typename Soln_cState>
Complex Haselbacher_Filter<Soln_pState,Soln_cState>::G_function_1D_110(Cell3D &theCell, Neighbours &theNeighbours, double &wave_number, RowVector &w, double &w0) {
    Cell3D kCell;
    kCell.Xc = Vector3D(wave_number/sqrt(TWO),wave_number/sqrt(TWO),ZERO);    
    return (G_function(theCell, theNeighbours, kCell, w, w0));
}

template<typename Soln_pState, typename Soln_cState>
Complex Haselbacher_Filter<Soln_pState,Soln_cState>::G_function_1D(Cell3D &theCell, Neighbours &theNeighbours, double &wave_number, RowVector &w) {
    Cell3D kCell;
    kCell.Xc = Vector3D(wave_number/sqrt(THREE),wave_number/sqrt(THREE),wave_number/sqrt(THREE));    
    return (G_function(theCell, theNeighbours, kCell, w));
}

template<typename Soln_pState, typename Soln_cState>
Complex Haselbacher_Filter<Soln_pState,Soln_cState>::G_function_1D(Cell3D &theCell, Neighbours &theNeighbours, double &wave_number, RowVector &w, double &w0) {
    Cell3D kCell;
    kCell.Xc = Vector3D(wave_number/sqrt(THREE),wave_number/sqrt(THREE),wave_number/sqrt(THREE));    
    return (G_function(theCell, theNeighbours, kCell, w, w0));
}

template<typename Soln_pState, typename Soln_cState>
Complex Haselbacher_Filter<Soln_pState,Soln_cState>::G_function(Cell3D &theCell, Neighbours &theNeighbours, Cell3D &kCell, RowVector &w, double &w0) {
    Complex G(0,0);
    double x,y,z,x0,y0,z0;
    x0 = theCell.Xc.x;
    y0 = theCell.Xc.y;
    z0 = theCell.Xc.z;
    for(int n=0; n<theNeighbours.number_of_neighbours; n++) {
        x=theNeighbours.neighbour[n].Xc.x;
        y=theNeighbours.neighbour[n].Xc.y;
        z=theNeighbours.neighbour[n].Xc.z;
        G += w(n) * exp(I*(kCell.Xc.x*(x0-x) + kCell.Xc.y*(y0-y) + kCell.Xc.z*(z0-z)));
    }
    G = Complex(w0,w0) + (ONE-w0)*G;
    return G;
}

template<typename Soln_pState, typename Soln_cState>
Complex Haselbacher_Filter<Soln_pState,Soln_cState>::G_function(Cell3D &theCell, Neighbours &theNeighbours, Cell3D &kCell, RowVector &w) {
    Complex G(0,0);
    double x,y,z,x0,y0,z0;
    x0 = theCell.Xc.x;
    y0 = theCell.Xc.y;
    z0 = theCell.Xc.z;
    for(int n=0; n<theNeighbours.number_of_neighbours; n++) {
        x=theNeighbours.neighbour[n].Xc.x;
        y=theNeighbours.neighbour[n].Xc.y;
        z=theNeighbours.neighbour[n].Xc.z;
        G += w(n) * exp(I*(kCell.Xc.x*(x0-x) + kCell.Xc.y*(y0-y) + kCell.Xc.z*(z0-z)));
    }
    return G;
}

template<typename Soln_pState, typename Soln_cState>
Complex Haselbacher_Filter<Soln_pState,Soln_cState>::dG_function_1D(Cell3D &theCell, Neighbours &theNeighbours, double &wave_number, RowVector &w) {
    Cell3D kCell;
    kCell.Xc = Vector3D(wave_number/sqrt(THREE),wave_number/sqrt(THREE),wave_number/sqrt(THREE));    
    return (dG_function(theCell, theNeighbours, kCell, w));
}

template<typename Soln_pState, typename Soln_cState>
Complex Haselbacher_Filter<Soln_pState,Soln_cState>::dG_function_1D(Cell3D &theCell, Neighbours &theNeighbours, double &wave_number, RowVector &w, double &w0) {
    Cell3D kCell;
    kCell.Xc = Vector3D(wave_number/sqrt(THREE),wave_number/sqrt(THREE),wave_number/sqrt(THREE));    
    return (dG_function(theCell, theNeighbours, kCell, w, w0));
}

template<typename Soln_pState, typename Soln_cState>
Complex Haselbacher_Filter<Soln_pState,Soln_cState>::dG_function(Cell3D &theCell, Neighbours &theNeighbours, Cell3D &kCell, RowVector &w, double &w0) {
    Complex G(0,0);
    double x,y,z,x0,y0,z0;
    double abs_k = kCell.Xc.abs();
    if (abs_k == ZERO) return ZERO;
    x0 = theCell.Xc.x;
    y0 = theCell.Xc.y;
    z0 = theCell.Xc.z;
    for(int n=0; n<theNeighbours.number_of_neighbours; n++) {
        x=theNeighbours.neighbour[n].Xc.x;
        y=theNeighbours.neighbour[n].Xc.y;
        z=theNeighbours.neighbour[n].Xc.z;
        G += w(n) * I * (kCell.Xc.x/abs_k*(x0-x) + kCell.Xc.y/abs_k*(y0-y) + kCell.Xc.z/abs_k*(z0-z)) * exp(I*(kCell.Xc.x*(x0-x) + kCell.Xc.y*(y0-y) + kCell.Xc.z*(z0-z)));
    }
    G = (ONE-w0)*G;
    return G;
}

template<typename Soln_pState, typename Soln_cState>
Complex Haselbacher_Filter<Soln_pState,Soln_cState>::dG_function(Cell3D &theCell, Neighbours &theNeighbours, Cell3D &kCell, RowVector &w) {
    Complex G(0,0);
    double x,y,z,x0,y0,z0;
    double abs_k = kCell.Xc.abs();
    if (abs_k == ZERO) return ZERO;
    x0 = theCell.Xc.x;
    y0 = theCell.Xc.y;
    z0 = theCell.Xc.z;
    for(int n=0; n<theNeighbours.number_of_neighbours; n++) {
        x=theNeighbours.neighbour[n].Xc.x;
        y=theNeighbours.neighbour[n].Xc.y;
        z=theNeighbours.neighbour[n].Xc.z;
        G += w(n) * (kCell.Xc.x/abs_k*(x0-x) + kCell.Xc.y/abs_k*(y0-y) + kCell.Xc.z/abs_k*(z0-z)) * exp(I*(kCell.Xc.x*(x0-x) + kCell.Xc.y*(y0-y) + kCell.Xc.z*(z0-z)));
    }
    return G;
}


template<typename Soln_pState, typename Soln_cState>
double Haselbacher_Filter<Soln_pState,Soln_cState>::filter_quality(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell, double &kmax, int number_of_rings, int commutation_order, double weight_factor) {
        
    Neighbours theNeighbours(SolnBlk.Grid);
    theNeighbours.GetNeighbours(theCell,number_of_rings);
    
    DenseMatrix A = Matrix_A(theCell, theNeighbours, commutation_order);
    DiagonalMatrix W = Matrix_W(theCell, theNeighbours, weight_factor);
    DenseMatrix Z = (W*A).pseudo_inverse()*W;
    RowVector w = Z[0];
    double w0 = Calculate_relaxation_factor(theCell,theNeighbours,kmax,w);
    
    return filter_quality(theCell,theNeighbours,kmax,w,w0);
}


template<typename Soln_pState, typename Soln_cState>
double Haselbacher_Filter<Soln_pState,Soln_cState>::filter_quality(Cell3D &theCell, Neighbours &theNeighbours, double &kmax, RowVector &w, double &w0) {
    
    double sharpness  = filter_sharpness(theCell,theNeighbours,kmax,w,w0);
    double smoothness  = filter_smoothness(theCell,theNeighbours,kmax,w,w0,0.25);
    double uniformity = filter_uniformity(theCell,theNeighbours,kmax,w,w0,0.3); 
    
    double quality = sharpness * smoothness * uniformity;
        
    Print_3(sharpness, smoothness, uniformity);
    
    return quality;
}


template<typename Soln_pState, typename Soln_cState>
double Haselbacher_Filter<Soln_pState,Soln_cState>::filter_sharpness(Cell3D &theCell, Neighbours &theNeighbours, double &kmax, RowVector &w, double &w0) {
    
    double k_90 = Calculate_wavenumber_of_Gvalue(theCell,theNeighbours,kmax,w,w0,0.9);
    double k_10 = Calculate_wavenumber_of_Gvalue(theCell,theNeighbours,kmax,w,w0,0.1);
    double sharpness = max(ZERO,ONE - (k_10 - k_90)/kmax);
    return sqrt(sharpness);
}

template<typename Soln_pState, typename Soln_cState>
double Haselbacher_Filter<Soln_pState,Soln_cState>::filter_smoothness(Cell3D &theCell, Neighbours &theNeighbours, double &kmax, RowVector &w, double &w0, double max_smoothness) {
    
    double k_10 = Calculate_wavenumber_of_Gvalue(theCell,theNeighbours,kmax,w,w0,0.1);
    double k_extreme = Calculate_wavenumber_of_dGvalue(theCell,theNeighbours,k_10,kmax,w,w0,ZERO);
    Complex G = G_function_1D(theCell,theNeighbours,k_extreme,w,w0);
    double smoothness = fabs(real(G)) > max_smoothness ? ZERO : fabs(fabs(real(G))-max_smoothness)/max_smoothness;
    return smoothness;
}

template<typename Soln_pState, typename Soln_cState>
double Haselbacher_Filter<Soln_pState,Soln_cState>::filter_uniformity(Cell3D &theCell, Neighbours &theNeighbours, double &kmax, RowVector &w, double &w0, double max_distortion) {
        
    Complex G_111 = G_function_1D(theCell,theNeighbours,kmax,w,w0);
    Complex G_100 = G_function_1D_100(theCell,theNeighbours,kmax,w,w0);
    Complex G_110 = G_function_1D_110(theCell,theNeighbours,kmax,w,w0);

    double uniformity_100 = fabs(real(G_100)-real(G_111)) > max_distortion ? ZERO : fabs(fabs(real(G_100)-real(G_111))-max_distortion)/max_distortion;
    double uniformity_110 = fabs(real(G_110)-real(G_111)) > max_distortion ? ZERO : fabs(fabs(real(G_110)-real(G_111))-max_distortion)/max_distortion;
    
    return min(uniformity_100,uniformity_110);
}


/* ------- calculate relaxation factor such that G_max = 0 ------- */
template<typename Soln_pState, typename Soln_cState>
double Haselbacher_Filter<Soln_pState,Soln_cState>::
Calculate_relaxation_factor(Cell3D &theCell, Neighbours &theNeighbours, double &kmax, RowVector &w) {
    
    double maximum_wave_number = sqrt(THREE)*kmax;
    Complex G_max = G_function_1D(theCell,theNeighbours,maximum_wave_number,w);
    
    double a,b,w0;
    int counter;
    a=-1.0, b=1.0;
    counter = 0;
    w0 = (a+b)/TWO;
    if (real(a+(ONE-a)*G_max) * real(b+(ONE-b)*G_max) >= ZERO) {
        cout << "relaxation_factor : root bracketed? " << endl;
    }
    while( fabs(real(w0+(ONE-w0)*G_max)) >= 0.001 ) {
        if( real(a+(ONE-a)*G_max) * real(w0+(ONE-w0)*G_max) <= ZERO)
            b = w0;
        else
            a = w0;
        w0 = (a+b)/TWO;
        counter++;
        if(counter >= 1000) {
            cout << "max reached for relaxation_factor" << endl;
            return 0.0;
        }
    }
    return w0;
}





/* ------- calculate k_HALF where G(k_HALF) = 0.5 -------- */
template<typename Soln_pState, typename Soln_cState>
double Haselbacher_Filter<Soln_pState,Soln_cState>::
Calculate_wavenumber_of_Gvalue(Cell3D &theCell, Neighbours &theNeighbours, double &kmax, RowVector &w, double &w0, double G_value) {
    double a,b;
    int counter = 0;
    a=0.0, b=sqrt(THREE)*kmax;
    double k_value = (a+b)/TWO;
    
    Complex G = G_function_1D(theCell,theNeighbours,k_value,w,w0);
    if( (real(G_function_1D(theCell,theNeighbours,a,w,w0))-G_value) * (real(G_function_1D(theCell,theNeighbours,b,w,w0))-G_value) >= ZERO) {
        cout << "Gvalue : root bracketed?" << endl;
    }
    while( fabs(real(G)-G_value) >= 0.001 ) {
        
        if( (real(G_function_1D(theCell,theNeighbours,a,w,w0))-G_value) * (real(G)-G_value) <= ZERO)
            b = k_value;
        else
            a = k_value;
        k_value = (a+b)/TWO;
        G = G_function_1D(theCell,theNeighbours,k_value,w,w0);
        counter++;
        if(counter >= 50) {
            cout << "max reached for Gvalue" << endl;
            break;
        }
    }
    return k_value;
}


template<typename Soln_pState, typename Soln_cState>
double Haselbacher_Filter<Soln_pState,Soln_cState>::
Calculate_wavenumber_of_dGvalue(Cell3D &theCell, Neighbours &theNeighbours, double kmin, double &kmax, RowVector &w, double &w0, double dG_value) {
    double a,b;
    int counter = 0;
    a=kmin, b=sqrt(THREE)*kmax;
    double k_value = (a+b)/TWO;
    
    Complex dG = dG_function_1D(theCell,theNeighbours,k_value,w,w0);
    while( fabs(real(dG)-dG_value) >= 0.001 ) {
        
        if( (real(dG_function_1D(theCell,theNeighbours,b,w,w0))-dG_value) * (real(dG)-dG_value) >= ZERO)
            b = k_value;
        else
            a = k_value;
        k_value = (a+b)/TWO;
        dG = dG_function_1D(theCell,theNeighbours,k_value,w,w0);
        counter++;
        if(counter >= 50) {
            cout << "max reached for dGvalue" << endl;
            return (sqrt(THREE)*kmax);
            break;
        }
    }
    return k_value;
}



template<typename Soln_pState, typename Soln_cState>
void Haselbacher_Filter<Soln_pState,Soln_cState>::transfer_function(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell, double &kmax) {
    
    Neighbours theNeighbours(SolnBlk.Grid);
    theNeighbours.GetNeighbours(theCell,number_of_rings);
    
    DenseMatrix A = Matrix_A(theCell, theNeighbours);

    
    
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
    
    
    
    
    
    
    

    double *k_111 = new double [N];
    double *k_110 = new double [N];
    double *k_100 = new double [N];
    double *G_111 = new double [N];
    double *G_110 = new double [N];
    double *G_100 = new double [N];
    Gnuplot_Control h1, h2, h3;
    h1.gnuplot_init();
    h2.gnuplot_init();
    h3.gnuplot_init();

    h1.gnuplot_setstyle("lines") ;
    h2.gnuplot_setstyle("lines") ;
    h3.gnuplot_setstyle("lines") ;

    h1.gnuplot_cmd("set grid");
    h2.gnuplot_cmd("set grid");
    h3.gnuplot_cmd("set grid");

    h1.gnuplot_set_xlabel("k");
    h2.gnuplot_set_xlabel("k");
    h3.gnuplot_set_xlabel("k");

    h1.gnuplot_set_ylabel("G(k)");
    h2.gnuplot_set_ylabel("G(k)");
    h3.gnuplot_set_ylabel("G(k)");

    h1.gnuplot_set_title("Transfer function 111");
    h2.gnuplot_set_title("Transfer function 110");
    h3.gnuplot_set_title("Transfer function 100");

    
    
    for (double weight_factor=0.5; weight_factor<=10.0; weight_factor+=0.5) {
        
        DiagonalMatrix W = Matrix_W(theCell, theNeighbours,weight_factor);
        DenseMatrix Z = (W*A).pseudo_inverse()*W;
        RowVector w = Z[0];
        
        double w0 = Calculate_relaxation_factor(theCell,theNeighbours,kmax,w);
        
        /* ------- calculate k_HALF where G(k_HALF) = 0.5 -------- */
        double k_HALF = Calculate_wavenumber_of_Gvalue(theCell,theNeighbours,kmax,w,w0,HALF);

        double FQ = filter_quality(theCell, theNeighbours, kmax, w, w0);
        //double FQ = filter_quality(SolnBlk, theCell, kmax, number_of_rings, commutation_order, weight_factor);
        // FGR = FGR(weight,commutation_order,number_of_rings)
        cout << "weight = " << weight_factor << "    w0 = " << w0 << "   k_HALF = " << k_HALF << "    FGR = " << kmax/k_HALF << "    FQ = " << FQ << endl;
        
        
        

//        int i,j,k;
//        for (int ii=0; ii<nK; ii++) {
//            K[ii].Ek = 0;
//            K[ii].Ek_smooth = 0;
//            for(int jj=0; jj<K[ii].N; jj++) {       // For all grid points with this value of abs_k
//                index = K[ii].indexes[jj];      // index corresponding to (i,j,l)
//                i = index/(N*N);
//                j = (index - N*N*i)/N;
//                k = index - N*(j+N*i);
//                K[ii].Ek += real(G[i][j][k])/double(K[ii].N);
//                K[ii].Ek_smooth += imag(G[i][0][0])/double(K[ii].N);
//            }
//            kk[ii]=K[ii].k;
//            GR[ii]=K[ii].Ek;
//        }
        
        for (int i=0; i<N; i++) {
            k_111[i] = i*sqrt(THREE)*kmax/(N-1.0);
            k_110[i] = i*sqrt(TWO)*kmax/(N-1.0);
            k_100[i] = i*kmax/(N-1.0);
            G_111[i]=real(G_function_1D(theCell,theNeighbours,k_111[i],w,w0));
            G_110[i]=real(G_function_1D_110(theCell,theNeighbours,k_110[i],w,w0));
            G_100[i]=real(G_function_1D_100(theCell,theNeighbours,k_100[i],w,w0));
            k_111[i]/=kmax;
            k_110[i]/=kmax;
            k_100[i]/=kmax;
        }
        

        
        
        //polyfit_smoothing(nK, kk, GR, 3, 2.5, true, GRs)  ;  
        /* ---------- output ----------- */
        
        char title[80];
        sprintf(title,"W = %3.2f  FGR = %3.2f",weight_factor,kmax/k_HALF);
        
        
        h1.gnuplot_plot1d_var2(k_111,G_111,N,title);
        h2.gnuplot_plot1d_var2(k_110,G_110,N,title);
        h3.gnuplot_plot1d_var2(k_100,G_100,N,title);
    }
    
    //cout << "filterquality(weight=5) = " << filter_quality(SolnBlk, theCell, kmax, number_of_rings, commutation_order, 5.0) << endl;
    


    
    DiagonalMatrix W = Matrix_W(theCell, theNeighbours, weight_factor);
    DenseMatrix Z = (W*A).pseudo_inverse()*W;
    RowVector w = Z[0];
    
    double w0 = Calculate_relaxation_factor(theCell,theNeighbours,kmax,w);

    
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
            for (int k=0; k<N; k++) {
                G[i][j][k] = G_function(theCell,theNeighbours,kCells[i][j][k],w,w0);
            }
        }
    }

    
    Output_Tecplot(kCells,G,N);
    
    //Output_Filter_types(SolnBlk,theCell,kmax);
    
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

    delete[] k_111;
    delete[] k_110;
    delete[] k_100;

    delete[] G_111;
    delete[] G_110;
    delete[] G_100;

    
    return; 
}

template<typename Soln_pState, typename Soln_cState>
Soln_pState Haselbacher_Filter<Soln_pState,Soln_cState>::filter(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell) {

    
    Neighbours theNeighbours(SolnBlk.Grid);
    theNeighbours.GetNeighbours(theCell,number_of_rings);
    
    if (UNIFORM_GRID) {
        if (!Assigned_w_Uniform_Grid) {
            Assign_w_Uniform_Grid(theCell,theNeighbours);
        }
        if(w_Uniform_Grid.size() != theNeighbours.number_of_neighbours)
            cout << "size(w) = "<< w_Uniform_Grid.size() << "    number_of_neighbours = " << theNeighbours.number_of_neighbours << endl;

        return LeastSquaresReconstruction_Uniform_Grid(SolnBlk, theCell, theNeighbours);
    }
    else{
        return LeastSquaresReconstruction(SolnBlk, theCell, theNeighbours);

    }
}


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
    Vector3D dr;
    
    Delta *= weight_factor;
    
    for (int i=0; i<the_number_of_neighbours; i++) {
        dr=(theNeighbours.neighbour[i].Xc - theCell.Xc);
        W(i) = sqrt(SIX/(PI*Delta.sqr())*exp(- SIX * dr.sqr()/Delta.sqr())) ;
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
    
    /* ---------- Filtered value is first element of x ----------- */
    Soln_pState temp;
    for (int j=1; j<=SolnBlk.NumVar(); j++)
        temp[j] = xrow(j-1);
    
    return temp;  

}


template<typename Soln_pState, typename Soln_cState>
Soln_pState Haselbacher_Filter<Soln_pState,Soln_cState>::LeastSquaresReconstruction_Uniform_Grid(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell, Neighbours &theNeighbours) {
        
    /* ---------------- Matrix b ---------------- */
    DenseMatrix b = Matrix_b(SolnBlk,theNeighbours);
        
    /* ----------- unknowns directly solved ------------ */
    RowVector x = w_Uniform_Grid * b;
    
    /* --------------- relaxation -------------- */
    RowVector x0(SolnBlk.NumVar());
    int I,J,K;
    I = theCell.I;
    J = theCell.J;
    K = theCell.K;
    for (int n=1; n<=SolnBlk.NumVar(); n++) {
        x0(n-1) = SolnBlk.W[I][J][K][n];
    }
    double w0 = relaxation_factor;
    
    x = x0*w0 + x*(ONE-w0);

    /* ---------- Filtered value is first element of x ----------- */
    Soln_pState temp;
    for (int j=1; j<=SolnBlk.NumVar(); j++)
        temp[j] = x(j-1);
    
    return temp;  
    
}

template<typename Soln_pState, typename Soln_cState>
void Haselbacher_Filter<Soln_pState,Soln_cState>::Assign_w_Uniform_Grid(Cell3D &theCell, Neighbours &theNeighbours) {
    DenseMatrix A = Matrix_A(theCell,theNeighbours);
    DiagonalMatrix W = Matrix_W(theCell,theNeighbours,weight_factor);
    DenseMatrix Z = (W*A).pseudo_inverse()*W;
    w_Uniform_Grid = Z[0]; // RowVector of weights
    Assigned_w_Uniform_Grid = true;
}



/* ------------ calculates the number of unknowns for the given order ----------- */
template <typename Soln_pState, typename Soln_cState>
int Haselbacher_Filter<Soln_pState,Soln_cState>::number_of_unknowns(int &commutation_order) {
    int number = 0;
    for (int k=0; k<= commutation_order; k++) {
        number += number_of_terms_of_degree(k);
    }
    return number;
}

/* ------------ calculates the number of unknowns for the given order ----------- */
template <typename Soln_pState, typename Soln_cState>
int Haselbacher_Filter<Soln_pState,Soln_cState>::number_of_unknowns(void) {
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




template<typename Soln_pState, typename Soln_cState>
int Haselbacher_Filter<Soln_pState,Soln_cState>::
Output_Filter_types(Hexa_Block<Soln_pState,Soln_cState> &SolnBlk, Cell3D &theCell, double kmax) {
    
    
    int N_number_of_rings = 5;
    int N_weight_factor = 59;   double max_weight_factor = 30.0;
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
        
        
        
        for (int j=0; j<N_commutation_order; j++) {
            int commutation_order = j+1;
            commutation_orders[j] = commutation_order;
            
            DenseMatrix A = Matrix_A(theCell, theNeighbours, commutation_order);
            
            
            for (int k=0; k<N_weight_factor; k++) {
                double weight_factor = 1.0 + k*(max_weight_factor-1.0)/(N_weight_factor-1.0);
                weight_factors[k]=weight_factor;
                
                Print_3(number_of_rings,commutation_order,weight_factor);
                DiagonalMatrix W = Matrix_W(theCell, theNeighbours,weight_factor);
                DenseMatrix Z = (W*A).pseudo_inverse()*W;
                RowVector w = Z[0];
                
                double w0 = Calculate_relaxation_factor(theCell,theNeighbours,kmax,w);
                
                /* ------- calculate k_HALF where G(k_HALF) = 0.5 -------- */
                double k_HALF = Calculate_wavenumber_of_Gvalue(theCell,theNeighbours,kmax,w,w0,HALF);
                
                FGR[i][j][k] = kmax/k_HALF;
                Q[i][j][k] = filter_quality(theCell, theNeighbours, kmax, w, w0);
                sharpness[i][j][k] = filter_sharpness(theCell, theNeighbours, kmax, w, w0);
                smoothness[i][j][k] = filter_smoothness(theCell, theNeighbours, kmax, w, w0, 0.25);
                uniformity[i][j][k] = filter_uniformity(theCell, theNeighbours, kmax, w, w0, 0.3);
                cost[i][j][k] = theNeighbours.number_of_neighbours * number_of_unknowns(commutation_order);
                
                // FGR = FGR(weight,commutation_order,number_of_rings)
                
                //cout << "weight = " << weight_factor << "    w0 = " << w0 << "    FGR = " << FGR[i][j][k] << "    FQ = " << Q[i][j][k] << endl;
                
                
                
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

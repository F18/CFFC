/*
 *  Neighbours.h
 *  CFFC_3D_Xcode
 *
 *  Created by Willem Deconinck on 25/04/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _LES_FILTERS_NEIGHBOURS_INCLUDED
#define _LES_FILTERS_NEIGHBOURS_INCLUDED


#include <cstdlib>     // defines the drand48() function
#include <ctime>       // defines the time() function

#define MAX_NUMBER_OF_NEIGHBOURING_RINGS_IN_LES_FILTER 8

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
    bool theCell_included;
    
    
    Vector3D Delta;
    Neighbours(Grid3D_Hexa_Block &Grid){
        Grid_ptr = &Grid;
        allocate();
        theCell_included = false;
    }
    Neighbours(Neighbours &anotherNeighbours){
        Grid_ptr = anotherNeighbours.Grid_ptr;
        theCell = anotherNeighbours.theCell;
        Allocated = false;
        theCell_included = false;
    }
    
    Neighbours(void) {
        Allocated = false;
        theCell_included = false;
    }
    void allocate(void) {
        neighbour = new Cell3D [points(MAX_NUMBER_OF_NEIGHBOURING_RINGS_IN_LES_FILTER)+1];
        Allocated = true;
    }
    void allocate(int number_of_rings) {
        neighbour = new Cell3D [points(number_of_rings)+1];
        Allocated = true;
    }
    
    void set_grid(Grid3D_Hexa_Block &Grid) {
        Grid_ptr = &Grid;
    }
    
    ~Neighbours(void) {
        if (Allocated)
            delete[] neighbour;
    }
    void GetNeighbours(Cell3D &theCell, int number_of_rings);
    void append_theCell(Cell3D &theCell);
    void delete_theCell(void);

    void GetNeighbours_Vasilyev(Cell3D &theCell, int number_of_rings);
    void GetNeighbours_Vasilyev_no_ghostcells(Cell3D &theCell, int number_of_rings);

    
    int Ki, Kj, Kk, Li, Lj, Lk;
    
    
    Cell3D random_cell(Cell3D &theCell, int imin, int imax, int jmin, int jmax, int kmin, int kmax);
    void shuffle_Neighbours(int number_of_cells);

    // useless because can calculate box - 1 point
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







// This is the original one with rings
inline void Neighbours::GetNeighbours(Cell3D &theCell, int number_of_rings) {
    number_of_neighbours = 0;
    assert(number_of_rings <= MAX_NUMBER_OF_NEIGHBOURING_RINGS_IN_LES_FILTER);    
    
    int imin = max(theCell.I-number_of_rings,Grid_ptr->ICl-Grid_ptr->Nghost);
    int imax = min(theCell.I+number_of_rings,Grid_ptr->ICu+Grid_ptr->Nghost);
    int jmin = max(theCell.J-number_of_rings,Grid_ptr->JCl-Grid_ptr->Nghost);
    int jmax = min(theCell.J+number_of_rings,Grid_ptr->JCu+Grid_ptr->Nghost);
    int kmin = max(theCell.K-number_of_rings,Grid_ptr->KCl-Grid_ptr->Nghost);
    int kmax = min(theCell.K+number_of_rings,Grid_ptr->KCu+Grid_ptr->Nghost);
    
    /* ---- assymetric stencil with same number of cells ---- */
    if (imin == Grid_ptr->ICl-Grid_ptr->Nghost)    imax = imin + 2*number_of_rings;
    if (jmin == Grid_ptr->JCl-Grid_ptr->Nghost)    jmax = jmin + 2*number_of_rings;
    if (kmin == Grid_ptr->KCl-Grid_ptr->Nghost)    kmax = kmin + 2*number_of_rings;
    if (imax == Grid_ptr->ICu+Grid_ptr->Nghost)    imin = imax - 2*number_of_rings;
    if (jmax == Grid_ptr->JCu+Grid_ptr->Nghost)    jmin = jmax - 2*number_of_rings;
    if (kmax == Grid_ptr->KCu+Grid_ptr->Nghost)    kmin = kmax - 2*number_of_rings;
    
    
    Vector3D Xmin(MILLION,MILLION,MILLION), Xmax(-MILLION,-MILLION,-MILLION);
    
    for (int i=imin; i<=imax; i++) {
        for (int j=jmin; j<=jmax; j++) {
            for (int k=kmin; k<=kmax; k++) {
                if (theCell != Grid_ptr->Cell[i][j][k]) {
                    neighbour[number_of_neighbours] = Grid_ptr->Cell[i][j][k];
                    Xmin.x = min(Xmin.x,neighbour[number_of_neighbours].Xc.x);
                    Xmin.y = min(Xmin.y,neighbour[number_of_neighbours].Xc.y);
                    Xmin.z = min(Xmin.z,neighbour[number_of_neighbours].Xc.z);
                    Xmax.x = max(Xmax.x,neighbour[number_of_neighbours].Xc.x);
                    Xmax.y = max(Xmax.y,neighbour[number_of_neighbours].Xc.y);
                    Xmax.z = max(Xmax.z,neighbour[number_of_neighbours].Xc.z);
                    number_of_neighbours++;
                }
            }
        }
    }
    
    Delta = (Xmax-Xmin);            // --> averaged dx dy dz over number of neighbouring rings
    Delta.x /= double(imax-imin);
    Delta.y /= double(jmax-jmin);
    Delta.z /= double(kmax-kmin);
    
    theCell_included = false;
    
}

inline void Neighbours::append_theCell(Cell3D &theCell) {
    if (!theCell_included) {
        neighbour[number_of_neighbours] = theCell;
        number_of_neighbours++;
        theCell_included = true; 
    }
}

inline void Neighbours::delete_theCell(void) {
    if (theCell_included)   number_of_neighbours--;
    theCell_included = false;
}


inline void Neighbours::GetNeighbours_Vasilyev(Cell3D &theCell, int number_of_rings) {
    number_of_neighbours = 0;
    assert(number_of_rings <= MAX_NUMBER_OF_NEIGHBOURING_RINGS_IN_LES_FILTER);    
    
    int imin = max(theCell.I-number_of_rings,Grid_ptr->ICl-Grid_ptr->Nghost);
    int jmin = max(theCell.J-number_of_rings,Grid_ptr->JCl-Grid_ptr->Nghost);
    int kmin = max(theCell.K-number_of_rings,Grid_ptr->KCl-Grid_ptr->Nghost);
    int imax = min(theCell.I+number_of_rings,Grid_ptr->ICu+Grid_ptr->Nghost);
    int jmax = min(theCell.J+number_of_rings,Grid_ptr->JCu+Grid_ptr->Nghost);
    int kmax = min(theCell.K+number_of_rings,Grid_ptr->KCu+Grid_ptr->Nghost);
    
    /* ---- assymetric stencil with same number of cells ---- */
    if (imin == Grid_ptr->ICl-Grid_ptr->Nghost)    imax = imin + 2*number_of_rings;
    if (jmin == Grid_ptr->JCl-Grid_ptr->Nghost)    jmax = jmin + 2*number_of_rings;
    if (kmin == Grid_ptr->KCl-Grid_ptr->Nghost)    kmax = kmin + 2*number_of_rings;
    if (imax == Grid_ptr->ICu+Grid_ptr->Nghost)    imin = imax - 2*number_of_rings;
    if (jmax == Grid_ptr->JCu+Grid_ptr->Nghost)    jmin = jmax - 2*number_of_rings;
    if (kmax == Grid_ptr->KCu+Grid_ptr->Nghost)    kmin = kmax - 2*number_of_rings;
    
    Ki = theCell.I - imin;      Li = imax - theCell.I;
    Kj = theCell.J - jmin;      Lj = jmax - theCell.J;
    Kk = theCell.K - kmin;      Lk = kmax - theCell.K;
    
    Vector3D Xmin(MILLION,MILLION,MILLION), Xmax(-MILLION,-MILLION,-MILLION);

    for (int i=imin; i<=imax; i++) {
        for (int j=jmin; j<=jmax; j++) {
            for (int k=kmin; k<=kmax; k++) {   // includes theCell !!!
                neighbour[number_of_neighbours] = Grid_ptr->Cell[i][j][k];
                Xmin.x = min(Xmin.x,neighbour[number_of_neighbours].Xc.x);
                Xmin.y = min(Xmin.y,neighbour[number_of_neighbours].Xc.y);
                Xmin.z = min(Xmin.z,neighbour[number_of_neighbours].Xc.z);
                Xmax.x = max(Xmax.x,neighbour[number_of_neighbours].Xc.x);
                Xmax.y = max(Xmax.y,neighbour[number_of_neighbours].Xc.y);
                Xmax.z = max(Xmax.z,neighbour[number_of_neighbours].Xc.z);
                number_of_neighbours++;
            }
        }
    }
    
    Delta = (Xmax-Xmin);            // --> averaged dx dy dz over number of neighbouring rings
    Delta.x /= double(imax-imin);
    Delta.y /= double(jmax-jmin);
    Delta.z /= double(kmax-kmin);
    
}



inline void Neighbours::GetNeighbours_Vasilyev_no_ghostcells(Cell3D &theCell, int number_of_rings) {
    number_of_neighbours = 0;
    assert(number_of_rings <= MAX_NUMBER_OF_NEIGHBOURING_RINGS_IN_LES_FILTER);    
    
    int imin = max(theCell.I-number_of_rings,Grid_ptr->ICl);
    int jmin = max(theCell.J-number_of_rings,Grid_ptr->JCl);
    int kmin = max(theCell.K-number_of_rings,Grid_ptr->KCl);
    int imax = min(theCell.I+number_of_rings,Grid_ptr->ICu);
    int jmax = min(theCell.J+number_of_rings,Grid_ptr->JCu);
    int kmax = min(theCell.K+number_of_rings,Grid_ptr->KCu);
    
    /* ---- assymetric stencil with same number of cells ---- */
    if (imin == Grid_ptr->ICl)    imax = imin + 2*number_of_rings;
    if (jmin == Grid_ptr->JCl)    jmax = jmin + 2*number_of_rings;
    if (kmin == Grid_ptr->KCl)    kmax = kmin + 2*number_of_rings;
    if (imax == Grid_ptr->ICu)    imin = imax - 2*number_of_rings;
    if (jmax == Grid_ptr->JCu)    jmin = jmax - 2*number_of_rings;
    if (kmax == Grid_ptr->KCu)    kmin = kmax - 2*number_of_rings;
    
    Ki = theCell.I - imin;      Li = imax - theCell.I;
    Kj = theCell.J - jmin;      Lj = jmax - theCell.J;
    Kk = theCell.K - kmin;      Lk = kmax - theCell.K;
    
    for (int i=imin; i<=imax; i++) {
        for (int j=jmin; j<=jmax; j++) {
            for (int k=kmin; k<=kmax; k++) {   // includes theCell !!!
                neighbour[number_of_neighbours] = Grid_ptr->Cell[i][j][k];
                number_of_neighbours++;
            }
        }
    }
    
}




//#define ODD -1.0
//#define EVEN 1.0
//inline void Neighbours::GetNeighbours(Cell3D &theCell, int number_of_rings) {
//    number_of_neighbours = 0;
//    
//    number_of_rings = 5;
//
//    
//    assert(number_of_rings <= MAX_NUMBER_OF_NEIGHBOURING_RINGS_IN_LES_FILTER);    
//    
//    
//    
//    int imin = max(theCell.I-number_of_rings,Grid_ptr->ICl-Grid_ptr->Nghost);
//    int imax = min(theCell.I+number_of_rings,Grid_ptr->ICu+Grid_ptr->Nghost);
//    int jmin = max(theCell.J-number_of_rings,Grid_ptr->JCl-Grid_ptr->Nghost);
//    int jmax = min(theCell.J+number_of_rings,Grid_ptr->JCu+Grid_ptr->Nghost);
//    int kmin = max(theCell.K-number_of_rings,Grid_ptr->KCl-Grid_ptr->Nghost);
//    int kmax = min(theCell.K+number_of_rings,Grid_ptr->KCu+Grid_ptr->Nghost);
//    int I = theCell.I;
//    int J = theCell.J;
//    int K = theCell.K;
//    
//    Vector3D Xmin(MILLION,MILLION,MILLION), Xmax(-MILLION,-MILLION,-MILLION);
//    
//    int seed = 1; // = time(NULL);      // assigns the current time to the seed
//    srand48(seed);                      // changes the seed for drand48()
//
//    cout << "max_pts = " << points(number_of_rings) << endl;
//    
////    for (int i=0; i<number_of_cells; i++) {
////        neighbour[number_of_neighbours] = random_cell(theCell,imin,imax,jmin,jmax,kmin,kmax);
////        Xmin.x = min(Xmin.x,neighbour[number_of_neighbours].Xc.x);
////        Xmin.y = min(Xmin.y,neighbour[number_of_neighbours].Xc.y);
////        Xmin.z = min(Xmin.z,neighbour[number_of_neighbours].Xc.z);
////        Xmax.x = max(Xmax.x,neighbour[number_of_neighbours].Xc.x);
////        Xmax.y = max(Xmax.y,neighbour[number_of_neighbours].Xc.y);
////        Xmax.z = max(Xmax.z,neighbour[number_of_neighbours].Xc.z);
////        number_of_neighbours++;
////    }                    
//                    
//    for (int i=imin; i<=imax; i++) {
//        for (int j=jmin; j<=jmax; j++) {
//            for (int k=kmin; k<=kmax; k++) {
//                if (theCell != Grid_ptr->Cell[i][j][k]) {
//                    neighbour[number_of_neighbours] = Grid_ptr->Cell[i][j][k];
//                    Xmin.x = min(Xmin.x,neighbour[number_of_neighbours].Xc.x);
//                    Xmin.y = min(Xmin.y,neighbour[number_of_neighbours].Xc.y);
//                    Xmin.z = min(Xmin.z,neighbour[number_of_neighbours].Xc.z);
//                    Xmax.x = max(Xmax.x,neighbour[number_of_neighbours].Xc.x);
//                    Xmax.y = max(Xmax.y,neighbour[number_of_neighbours].Xc.y);
//                    Xmax.z = max(Xmax.z,neighbour[number_of_neighbours].Xc.z);
//                    number_of_neighbours++;
//                }
//            }
//        }
//    }
//    shuffle_Neighbours(135);
//                                     
//                    
//                    
//                    /*if ( (abs(I-i)==1)  && (abs(J-j)==1) && (abs(K-k)==1) ) {
//                        neighbour[number_of_neighbours] = Grid_ptr->Cell[i][j][k];
//                        Xmin.x = min(Xmin.x,neighbour[number_of_neighbours].Xc.x);
//                        Xmin.y = min(Xmin.y,neighbour[number_of_neighbours].Xc.y);
//                        Xmin.z = min(Xmin.z,neighbour[number_of_neighbours].Xc.z);
//                        Xmax.x = max(Xmax.x,neighbour[number_of_neighbours].Xc.x);
//                        Xmax.y = max(Xmax.y,neighbour[number_of_neighbours].Xc.y);
//                        Xmax.z = max(Xmax.z,neighbour[number_of_neighbours].Xc.z);
//                        number_of_neighbours++;
//                    }
//                    if ( (abs(I-i)==3)  && (abs(J-j)==3) && (abs(K-k)==3) ) {
//                        neighbour[number_of_neighbours] = Grid_ptr->Cell[i][j][k];
//                        Xmin.x = min(Xmin.x,neighbour[number_of_neighbours].Xc.x);
//                        Xmin.y = min(Xmin.y,neighbour[number_of_neighbours].Xc.y);
//                        Xmin.z = min(Xmin.z,neighbour[number_of_neighbours].Xc.z);
//                        Xmax.x = max(Xmax.x,neighbour[number_of_neighbours].Xc.x);
//                        Xmax.y = max(Xmax.y,neighbour[number_of_neighbours].Xc.y);
//                        Xmax.z = max(Xmax.z,neighbour[number_of_neighbours].Xc.z);
//                        number_of_neighbours++;
//                    }
//                        
//                    
//                    else if ( pow(-ONE,(i-I)+(j-J)+(k-K)) == EVEN) {
//                        neighbour[number_of_neighbours] = Grid_ptr->Cell[i][j][k];
//                        Xmin.x = min(Xmin.x,neighbour[number_of_neighbours].Xc.x);
//                        Xmin.y = min(Xmin.y,neighbour[number_of_neighbours].Xc.y);
//                        Xmin.z = min(Xmin.z,neighbour[number_of_neighbours].Xc.z);
//                        Xmax.x = max(Xmax.x,neighbour[number_of_neighbours].Xc.x);
//                        Xmax.y = max(Xmax.y,neighbour[number_of_neighbours].Xc.y);
//                        Xmax.z = max(Xmax.z,neighbour[number_of_neighbours].Xc.z);
//                        number_of_neighbours++;
//                    }*/
////                }
////            }
////        }
////    }
//    
//    Delta = (Xmax-Xmin);            // --> averaged dx dy dz over number of neighbouring rings
//    Delta.x /= double(imax-imin);
//    Delta.y /= double(jmax-jmin);
//    Delta.z /= double(kmax-kmin);
//}


inline void Neighbours::shuffle_Neighbours(int number_of_cells){
    for (int i=0; i<number_of_cells; i++) {
        int r = i + (rand() % (number_of_neighbours-i)); // Random remaining position.
        
        // switch neighbour i and r
        Cell3D temp = neighbour[i]; 
        neighbour[i] = neighbour[r]; 
        neighbour[r] = temp;
    }
    number_of_neighbours = number_of_cells;
}

inline Cell3D Neighbours::random_cell(Cell3D &theCell, int imin, int imax, int jmin, int jmax, int kmin, int kmax){
    Cell3D theRandomCell;
    
    bool repeat_condition = false;
    do {
        int i = imin + rand() % (imax-imin+1); 
        int j = jmin + rand() % (jmax-jmin+1);
        int k = kmin + rand() % (kmax-kmin+1);
        
        theRandomCell = Grid_ptr->Cell[i][j][k];
        
        for (int i=0; i<number_of_neighbours; i++) {
            if (theRandomCell == neighbour[i] || theRandomCell == theCell) {
                repeat_condition = true;
                cout << "                   already present: " << theRandomCell << endl;
                break;
            }
        }
    } while (repeat_condition);
    cout << "return : " << theRandomCell << endl;
    return  theRandomCell;
}


inline int Neighbours::pointsPerEdge(int ring) {
    return 1+(ring-1)*2;
}

inline int Neighbours::pointsPerFace(int ring) {
    int n=1;
    for (int i=1; i<=ring; i++) {
        n += (i-1)*8;
    }
    return n;
}

inline int Neighbours::points(int num_of_rings) {
    int n=0;
    for (int i=1; i<=num_of_rings; i++) {
        n += 8 + 12*pointsPerEdge(i) + 4*pointsPerFace(i);
    }
    return n;
}

#endif
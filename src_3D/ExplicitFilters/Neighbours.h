/*
 *  Neighbours.h
 *  CFFC_3D_Xcode
 *
 *  Created by Willem Deconinck on 25/04/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _LES_FILTERS_NEIGHBOURS_NEW_INCLUDED
#define _LES_FILTERS_NEIGHBOURS_NEW_INCLUDED


#include <cstdlib>     // defines the drand48() function
#include <ctime>       // defines the time() function
#include "../Grid/Cell3D.h"
#include "../Grid/Grid3DHexaBlock.h"

#define MAX_NUMBER_OF_NEIGHBOURING_RINGS_IN_LES_FILTER 8

/**
 * CLASS : Neighbours
 */
class Neighbours{
public:
    Cell3D *neighbour;
   // Cell3D *theCell;
    Grid3D_Hexa_Block *Grid_ptr;
    
    /* for 1D-filters */
    Cell3D *neighbour_x;
    Cell3D *neighbour_y;
    Cell3D *neighbour_z;
    bool vasilyev_neighbours_allocated;
    
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
        vasilyev_neighbours_allocated = false;
    }
    Neighbours(Neighbours &anotherNeighbours){
        Grid_ptr = anotherNeighbours.Grid_ptr;
        //theCell = anotherNeighbours.theCell;
        Allocated = false;
        theCell_included = false;
        vasilyev_neighbours_allocated = false;
    }
    
    Neighbours(void) {
        Allocated = false;
        theCell_included = false;
        vasilyev_neighbours_allocated = false;
    }
    void allocate(void) {
        allocate(MAX_NUMBER_OF_NEIGHBOURING_RINGS_IN_LES_FILTER);
    }
    void allocate(int number_of_rings) {
        deallocate();
        neighbour = new Cell3D [int(pow(TWO*number_of_rings+ONE,THREE))];
        Allocated = true;
    }
    
    void allocate_vasilyev(int Ni, int Nj, int Nk) {
        deallocate_vasilyev();
        neighbour_x = new Cell3D[Ni];
        neighbour_y = new Cell3D[Nj];
        neighbour_z = new Cell3D[Nk];
        vasilyev_neighbours_allocated = true;
    }
    
    void deallocate_vasilyev(void) {
        if (vasilyev_neighbours_allocated) {
            delete[] neighbour_x;   neighbour_x = NULL;
            delete[] neighbour_y;   neighbour_y = NULL;
            delete[] neighbour_z;   neighbour_z = NULL;
        }
        vasilyev_neighbours_allocated = false;
    }
    
    void set_grid(Grid3D_Hexa_Block &Grid) {
        Grid_ptr = &Grid;
    }
    
    void deallocate(void) {
        if (Allocated) {
            delete[] neighbour;
        }
        Allocated = false;
    }
    
    ~Neighbours(void) {
        deallocate();
        deallocate_vasilyev();
    }
    void GetNeighbours(Cell3D &theCell, int number_of_rings);
    void append_theCell(Cell3D &theCell);
    void delete_theCell(void);

    void GetNeighbours_Vasilyev(Cell3D &theCell, int number_of_rings);
    
    int Ki, Kj, Kk, Li, Lj, Lk, Ni, Nj, Nk;
    
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
    
    for (int i=imin; i<=imax; i++) {
        for (int j=jmin; j<=jmax; j++) {
            for (int k=kmin; k<=kmax; k++) {
                if (theCell != Grid_ptr->Cell[i][j][k]) {
                    neighbour[number_of_neighbours] = Grid_ptr->Cell[i][j][k];
                    number_of_neighbours++;
                }
            }
        }
    }

    Delta = theCell.dXc;
    
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
    
    Ki = theCell.I - imin;      Li = imax - theCell.I;      Ni = Ki+Li+1;
    Kj = theCell.J - jmin;      Lj = jmax - theCell.J;      Nj = Kj+Lj+1;
    Kk = theCell.K - kmin;      Lk = kmax - theCell.K;      Nk = Kk+Lk+1;

    allocate_vasilyev(Ni,Nj,Nk);
    
    for (int i=imin; i<=imax; i++) {
        for (int j=jmin; j<=jmax; j++) {
            for (int k=kmin; k<=kmax; k++) {   // includes theCell !!!
                neighbour[number_of_neighbours] = Grid_ptr->Cell[i][j][k];
                number_of_neighbours++;
            }
        }
    }
    
    for (int i=imin, n=0; i<=imax; i++, n++) {
        neighbour_x[n] = Grid_ptr->Cell[i][theCell.J][theCell.K];
    }
    for (int j=jmin, n=0; j<=jmax; j++, n++) {
        neighbour_y[n] = Grid_ptr->Cell[theCell.I][j][theCell.K];
    }
    for (int k=kmin, n=0; k<=kmax; k++, n++) { 
        neighbour_z[n] = Grid_ptr->Cell[theCell.I][theCell.J][k];
    }
        
    Delta = theCell.dXc;

}

#endif

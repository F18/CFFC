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
    
    static int max_number_of_neighbours;
    static int number_of_rings;
    bool Allocated;
    int number_of_neighbours;
    bool theCell_included;
    bool symmetric_stencil;
    
    
    Vector3D Delta;
    double Delta_1D;
    
    Neighbours(Grid3D_Hexa_Block &Grid){
        Grid_ptr = &Grid;
        Allocated = false;
        allocate();
        theCell_included = false;
    }
    
    Neighbours(Neighbours &anotherNeighbours){
        Grid_ptr = anotherNeighbours.Grid_ptr;
        //theCell = anotherNeighbours.theCell;
        Allocated = false;
        theCell_included = false;
    }
    
    Neighbours(void) {
        Allocated = false;
        theCell_included = false;
    }
    void allocate(void) {
        allocate(MAX_NUMBER_OF_NEIGHBOURING_RINGS_IN_LES_FILTER);
    }
    
    void allocate(int number_of_rings) {
        deallocate();
        neighbour = new Cell3D [int(pow(TWO*number_of_rings+ONE,THREE))];
        Allocated = true;
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
    }
    
    void GetNeighbours(Cell3D &theCell, int number_of_rings,int filter_type);
    void GetNeighbours_1D(Cell3D &theCell, int number_of_rings,int filter_type, int direction);
    
    void append_theCell(Cell3D &theCell);
    void delete_theCell(void);
    
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
inline void Neighbours::GetNeighbours(Cell3D &theCell, int number_of_rings,int filter_type) {
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
    
    Ki = theCell.I - imin;      Li = imax - theCell.I;      Ni = Ki+Li+1;
    Kj = theCell.J - jmin;      Lj = jmax - theCell.J;      Nj = Kj+Lj+1;
    Kk = theCell.K - kmin;      Lk = kmax - theCell.K;      Nk = Kk+Lk+1;
    
    if (Ki==Li && Kj==Lj && Kk==Lk){
        symmetric_stencil = true;
    } else {
        symmetric_stencil = false;
    }
    
    switch (filter_type) {
        case FILTER_TYPE_HASELBACHER:
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
            theCell_included = false;
            break;
        case FILTER_TYPE_VASILYEV:
            for (int i=imin; i<=imax; i++) {
                for (int j=jmin; j<=jmax; j++) {
                    for (int k=kmin; k<=kmax; k++) {   // includes theCell !!!
                        neighbour[number_of_neighbours] = Grid_ptr->Cell[i][j][k];
                        number_of_neighbours++;
                    }
                }
            }
            break;
    }
    
    Delta = theCell.dXc.vector_abs();
}



inline void Neighbours::GetNeighbours_1D(Cell3D &theCell, int number_of_rings,int filter_type, int direction) {
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
    
    Ki = theCell.I - imin;      Li = imax - theCell.I;      Ni = Ki+Li+1;
    Kj = theCell.J - jmin;      Lj = jmax - theCell.J;      Nj = Kj+Lj+1;
    Kk = theCell.K - kmin;      Lk = kmax - theCell.K;      Nk = Kk+Lk+1;
    
    
    switch (filter_type) {
        case FILTER_TYPE_HASELBACHER:
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
            switch (direction){
                case X_DIRECTION:
                    for (int i=imin; i<=imax; i++) {
                        if (theCell != Grid_ptr->Cell[i][theCell.J][theCell.K]) {
                            neighbour[number_of_neighbours] = Grid_ptr->Cell[i][theCell.J][theCell.K];
                            number_of_neighbours++;
                        }
                    }
                    break;
                case Y_DIRECTION:
                    for (int j=jmin; j<=jmax; j++) {
                        if (theCell != Grid_ptr->Cell[theCell.I][j][theCell.K]) {
                            neighbour[number_of_neighbours] = Grid_ptr->Cell[theCell.I][j][theCell.K];
                            number_of_neighbours++;
                        }
                    }
                    break;
                case Z_DIRECTION:
                    for (int k=kmin; k<=kmax; k++) { 
                        if (theCell != Grid_ptr->Cell[theCell.I][theCell.J][k]) {
                            neighbour[number_of_neighbours] = Grid_ptr->Cell[theCell.I][theCell.J][k];
                            number_of_neighbours++;
                        }
                    }
                    break;
                    
            }
            theCell_included = false;
            break;
            case FILTER_TYPE_VASILYEV:
            switch (direction){
                case X_DIRECTION:
                    for (int i=imin; i<=imax; i++) {
                            neighbour[number_of_neighbours] = Grid_ptr->Cell[i][theCell.J][theCell.K];
                            number_of_neighbours++;
                        
                    }
                    break;
                case Y_DIRECTION:
                    for (int j=jmin; j<=jmax; j++) {
                            neighbour[number_of_neighbours] = Grid_ptr->Cell[theCell.I][j][theCell.K];
                            number_of_neighbours++;
                    }
                    break;
                case Z_DIRECTION:
                    for (int k=kmin; k<=kmax; k++) { 
                            neighbour[number_of_neighbours] = Grid_ptr->Cell[theCell.I][theCell.J][k];
                            number_of_neighbours++;
                    }
                    break;
            }
            break;
    }
    
    switch (direction){
        case X_DIRECTION:
            Delta_1D = theCell.dXc.x;
            if (Ki==Li){
                symmetric_stencil = true;
            } else {
                symmetric_stencil = false;
            }
            break;
        case Y_DIRECTION:
            Delta_1D = theCell.dXc.y;
            if (Kj==Lj){
                symmetric_stencil = true;
            } else {
                symmetric_stencil = false;
            }
            break;
        case Z_DIRECTION:
            Delta_1D = theCell.dXc.z;
            if (Kk==Lk){
                symmetric_stencil = true;
            } else {
                symmetric_stencil = false;
            }
            break;
    }    
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


#endif

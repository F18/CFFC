/*
 *  LES_Filters.cpp
 *  CFFC_3D_Xcode
 *
 *  Created by Willem Deconinck on 25/03/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "LES_Filters.h"

#define MAX_NEIGHBOURS_IN_LES_FILTERS_NOT_DEFINED  -1
#define NUMBER_OF_RINGS_IN_LES_FILTERS_NOT_DEFINED -1


int Neighbours::max_number_of_neighbours = MAX_NEIGHBOURS_IN_LES_FILTERS_NOT_DEFINED;
//int Neighbours::number_of_rings = NUMBER_OF_RINGS_IN_LES_FILTERS_NOT_DEFINED;


void Neighbours::GetNeighbours(Cell3D &theCell, int number_of_rings) {
    number_of_neighbours = 0;
    assert(number_of_rings <= MAX_NUMBER_OF_NEIGHBOURING_RINGS_IN_LES_FILTER);
//    cout << "number_of_rings = " << number_of_rings << "   and num_of_rings = " << num_of_rings << endl;
//    
//    if (number_of_rings != num_of_rings) {
//        cout << "changing number_of_rings" << endl;
//        number_of_rings = num_of_rings;
//        max_number_of_neighbours = points(number_of_rings);
//        if(Allocated) delete[] neighbour;
//        neighbour = new Cell3D[max_number_of_neighbours];
//        Allocated = true;
//    }
//
//    cout << "number_of_rings = " << number_of_rings << endl;


    
    int imin = max(theCell.I-number_of_rings,Grid_ptr->ICl-Grid_ptr->Nghost);
    int imax = min(theCell.I+number_of_rings,Grid_ptr->ICu+Grid_ptr->Nghost);
    int jmin = max(theCell.J-number_of_rings,Grid_ptr->JCl-Grid_ptr->Nghost);
    int jmax = min(theCell.J+number_of_rings,Grid_ptr->JCu+Grid_ptr->Nghost);
    int kmin = max(theCell.K-number_of_rings,Grid_ptr->KCl-Grid_ptr->Nghost);
    int kmax = min(theCell.K+number_of_rings,Grid_ptr->KCu+Grid_ptr->Nghost);
    
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
}


int Neighbours::pointsPerEdge(int ring) {
    return 1+(ring-1)*2;
}

int Neighbours::pointsPerFace(int ring) {
    int n=1;
    for (int i=1; i<=ring; i++) {
        n += (i-1)*8;
    }
    return n;
}

int Neighbours::points(int num_of_rings) {
    int n=0;
    for (int i=1; i<=num_of_rings; i++) {
        n += 8 + 12*pointsPerEdge(i) + 4*pointsPerFace(i);
    }
    return n;
}

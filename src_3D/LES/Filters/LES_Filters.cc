/*
 *  LES_Filters.cpp
 *  CFFC_3D_Xcode
 *
 *  Created by Willem Deconinck on 25/03/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "LES_Filters.h"

void Neighbours::GetNeighbours(Cell3D &theCell) {
    number_of_neighbours = 0;
    
    int imin = max(theCell.I-1,Grid_ptr->ICl-Grid_ptr->Nghost);
    int imax = min(theCell.I+1,Grid_ptr->ICu+Grid_ptr->Nghost);
    int jmin = max(theCell.J-1,Grid_ptr->JCl-Grid_ptr->Nghost);
    int jmax = min(theCell.J+1,Grid_ptr->JCu+Grid_ptr->Nghost);
    int kmin = max(theCell.K-1,Grid_ptr->KCl-Grid_ptr->Nghost);
    int kmax = min(theCell.K+1,Grid_ptr->KCu+Grid_ptr->Nghost);
    
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
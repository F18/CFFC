/* Cell2D.cc:  Subroutines and data for 2D cells. */

/* Include 2D cell header file. */

#include "HO_Cell2D.h"

/********************************************************************
 * Cell2D_Cartesian -- Create storage and set cell lengths.         *
 ********************************************************************/
Vector2D Cell2D_Cartesian_HO::dx = Vector2D_ZERO;

Cell2D_Cartesian_HO Cell2D_Cartesian_HO::Cell2D_Cartesian_HO_ONE(ONE);

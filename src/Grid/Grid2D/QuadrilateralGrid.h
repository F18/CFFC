/* Euler2DQuad.h:  Header file defining 
                   2D Euler Quadrilateral Mesh Solution Classes. */

#ifndef _QUADRILATERAL_GRID_INCLUDED
#define _QUADRILATERAL_GRID_INCLUDED

/* Include 2D quadrilateral multiblock input header files. */

#include "Grid/Grid2D/Grid2DQuad.h"
#include "CFD/CFD.h"
#include "Reconstruction/Reconstruction2D/Reconstruct2DInput.h"

/**************************************************************************
 * Euler2D_Quad_Block -- Multiple Block External Subroutines for Mesh.    *
 **************************************************************************/

extern Grid2D_Quad_Block** Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                            Reconstruct2D_Input_Parameters &Input_Parameters);

extern Grid2D_Quad_Block** Broadcast_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                                      Reconstruct2D_Input_Parameters &Input_Parameters);

extern int Write_Multi_Block_Grid_Definition(Grid2D_Quad_Block **Grid_ptr,
                                             Reconstruct2D_Input_Parameters &Input_Parameters);

extern Grid2D_Quad_Block** Read_Multi_Block_Grid_Definition(Grid2D_Quad_Block **Grid_ptr,
                                                            Reconstruct2D_Input_Parameters &Input_Parameters);

extern int Write_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                  Reconstruct2D_Input_Parameters &Input_Parameters);

extern Grid2D_Quad_Block** Read_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                                 Reconstruct2D_Input_Parameters &Input_Parameters);

extern int Output_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                          Reconstruct2D_Input_Parameters &Input_Parameters);

extern int Output_Nodes_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                                Reconstruct2D_Input_Parameters &Input_Parameters);

extern int Output_Cells_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                                Reconstruct2D_Input_Parameters &Input_Parameters);

#endif /* _EULER2D_QUAD_INCLUDED  */

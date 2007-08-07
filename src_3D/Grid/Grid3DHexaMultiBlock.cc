/* Grid3DHexaMultiBlock.cc:  Member functions for 
                             3D hexahedral multiblock grid class. */

/* Include 3D hexahedral block grid header file. */

#ifndef GRID3D_HEXA_MULTIBLOCK_INCLUDED
#include "Grid3DHexaMultiBlock.h"
#endif //_GRID3D_HEXA_MULTIBLOCK_INCLUDED

/********************************************************
 * Routine: Allocate_Grid_Blocks                        *
 *                                                      *
 * Allocate memory for a 3D array of 3D hexahedral      *
 * multi-block grids.                                   *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block::Allocate(const int Ni, 
                                       const int Nj, 
                                       const int Nk) {

   assert( Ni >= 1 && Nj >= 1 && Nk >= 1 && !Allocated);

   NBlk_Idir = Ni; NBlk_Jdir = Nj; NBlk_Kdir = Nk; Allocated = 1;

   Grid_Blks = new Grid3D_Hexa_Block**[NBlk_Idir];
   for (int i = 0 ; i < NBlk_Idir ; ++i ) {
      Grid_Blks[i] = new Grid3D_Hexa_Block*[NBlk_Jdir];
      for (int j = 0 ; j <NBlk_Jdir ; ++j ) {
         Grid_Blks[i][j] = new Grid3D_Hexa_Block[NBlk_Kdir];
      }  /* endfor */
   }/* endfor */ 

   Allocated = 1;

}

/********************************************************
 * Routine: Deallocate_Multi_Block_Grid                 *
 *                                                      *
 * Deallocate memory for a 3D array of 3D hexahedral    *
 * multi-block grids.                                   *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block::Deallocate(void) {

   assert( NBlk_Idir >= 1 && NBlk_Jdir >= 1 && NBlk_Kdir >= 1 && Allocated); 

   for (int i = 0 ; i <= NBlk_Idir-1 ; ++i ) {
      for ( int j=0 ; j <= NBlk_Jdir-1 ; ++j ) {
	 for ( int k = NBlk_Kdir-1 ; k >= 0; --k ) {
            if (Grid_Blks[i][j][k].Used) Grid_Blks[i][j][k].deallocate();
	 }/* end for */
	 delete []Grid_Blks[i][j];
	 Grid_Blks[i][j]=NULL;
       }  /* endfor */
       
       delete []Grid_Blks[i];
       Grid_Blks[i] = NULL;
    }  /* endfor */

    delete []Grid_Blks;
    Grid_Blks = NULL;

    Allocated = 0;

}

/********************************************************
 * Routine: Copy                                        *
 *                                                      *
 * Make a copy of multiblock hexahedral grid Grid2.     *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block::Copy(Grid3D_Hexa_Multi_Block &Grid2) {

  if (Grid2.Allocated) {

    /* Ensure multiblock grid arrays have same dimensions. */

    if (Allocated && (NBlk_Idir != Grid2.NBlk_Idir ||
                      NBlk_Jdir != Grid2.NBlk_Jdir ||
                      NBlk_Kdir != Grid2.NBlk_Kdir) ) {
      Deallocate();
      Allocate(Grid2.NBlk_Idir, Grid2.NBlk_Jdir, Grid2.NBlk_Kdir);
    } else if (!Allocated) {
      Allocate(Grid2.NBlk_Idir, Grid2.NBlk_Jdir, Grid2.NBlk_Kdir);
    } /* endif */

    /* Copy each grid block. */

    for (int  k = 0 ; k < NBlk_Kdir ; ++k ) {
       for (int  j = 0 ; j < NBlk_Jdir ; ++j ) {
          for (int  i = 0 ; i < NBlk_Idir ; ++i ) {
	     if (Grid2.Grid_Blks[i][j][k].Used) 
                Grid_Blks[i][j][k].Copy(Grid2.Grid_Blks[i][j][k]);
          } /* endfor */
       } /* endfor */
    } /* endfor */

  } /* endif */

}

/********************************************************
 * Routine: Output                                      *
 *                                                      *
 * Writes a 3D array of 3D hexahedral multi-block       *
 * grids to the specified output stream for retrieval   *
 * and re-use purposes.                                 *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block::Output(ostream &Out_File) {
   
  Out_File << NBlk_Idir << " " 
           << NBlk_Jdir << " "
           << NBlk_Kdir << "\n";
  
  for (int  k = 0 ; k < NBlk_Kdir ; ++k ) {
     for (int  j = 0 ; j < NBlk_Jdir ; ++j ) {
       for (int  i = 0 ; i < NBlk_Idir ; ++i ) {
          Out_File << setprecision(14) << Grid_Blks[i][j][k] << setprecision(6);
       } /* endfor */
     } /* endfor */
  } /* endfor */

}

/********************************************************
 * Routine: Output_Tecplot                              *
 *                                                      *
 * Writes the nodes of a 3D array of 3D hexahedral      *
 * multi-block grids to the specified output stream in  *
 * a format suitable for plotting the grid with         *
 * TECPLOT.                                             *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block::Output_Tecplot(ostream &Out_File) {

    int block_number, i_output_title;
    block_number = 0;
    i_output_title = 1;
    
    for (int  k = 0 ; k < NBlk_Kdir ; ++k ) {
       for (int  j = 0; j < NBlk_Jdir ; ++j ) {
	  for (int i = 0 ; i < NBlk_Idir ; ++i ) {
             if (Grid_Blks[i][j][k].Used) {
                Grid_Blks[i][j][k].Output_Tecplot(block_number,
                                                  i_output_title,
                                                  Out_File);
	        block_number = block_number + 1;
                if (i_output_title) i_output_title = 0;
             } /* endif */
	  } /* endfor */
       } /* endfor */
    } /* endfor */

}

/********************************************************
 * Routine: Output_Nodes_Tecplot                        *
 *                                                      *
 * Writes the nodes of a 3D array of 3D hexahedral      *
 * multi-block grids to the specified output stream in  *
 * a format suitable for plotting the grid with         *
 * TECPLOT.  Include boundary nodes.                    *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block::Output_Nodes_Tecplot(ostream &Out_File) {

    int block_number, i_output_title;
    block_number = 0;
    i_output_title = 1;
  
    for (int  k = 0 ; k < NBlk_Kdir ; ++k ) {
       for (int  j = 0; j < NBlk_Jdir ; ++j ) {
          for (int i = 0 ; i < NBlk_Idir ; ++i ) {
	     if (Grid_Blks[i][j][k].Used) {
                Grid_Blks[i][j][k].Output_Nodes_Tecplot(block_number,
                                                        i_output_title,
                                                        Out_File);
                block_number = block_number + 1;
                if (i_output_title) i_output_title = 0;
	     } /* endif */
          } /* endfor */
       } /* endfor */
    } /* endfor */ 

}

/********************************************************
 * Routine: Output_Cells_Tecplot                        *
 *                                                      *
 * Writes the cells of a 3D array of 3D hexahedral      *
 * multi-block grids to the specified output stream in  *
 * a format suitable for plotting the grid with         *
 * TECPLOT.                                             *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block::Output_Cells_Tecplot(ostream &Out_File) {

    int block_number, i_output_title;
    block_number = 0;
    i_output_title = 1;

    for (int  k = 0 ; k <NBlk_Kdir ; ++k ) {
       for (int  j = 0; j <NBlk_Jdir ; ++j ) {
          for (int i = 0 ; i <NBlk_Idir ; ++i ) {
	     if (Grid_Blks[i][j][k].Used) {
                Grid_Blks[i][j][k].Output_Cells_Tecplot(block_number,
                                                        i_output_title,
                                                        Out_File);
                block_number = block_number + 1;
                if (i_output_title) i_output_title = 0;
	     } /* endif */
          } /* endfor */
       } /* endfor */
    }/* endfor */ 

}

/********************************************************
 * Routine: Output_Gnuplot                              *
 *                                                      *
 * Writes the nodes of a 3D array of 3D hexahedral      *
 * multi-block grids to the specified output stream in  *
 * a format suitable for plotting the grid with         *
 * GNUPLOT.                                             *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block::Output_Gnuplot(ostream &Out_File) {
  
    int  block_number, i_output_title;
    block_number = 0;
    i_output_title = 1;
    
    for (int k = 0 ; k < NBlk_Kdir ; ++k ) {
        for (int j = 0; j < NBlk_Jdir ; ++j ) {
	   for (int i =0 ; i < NBlk_Idir ; ++i ) {
	      if (Grid_Blks[i][j][k].Used) {
	          Grid_Blks[i][j][k].Output_Gnuplot(block_number,
			                            i_output_title,
			                            Out_File);
                 block_number = block_number + 1;
	         if (i_output_title) i_output_title = 0;
	      } /* endif */
	   } /* endfor */
        } /* endfor */
    } /* endfor */
    
}

/********************************************************
 * Routine: Create_Grid                                 *
 *                                                      *
 * Generates a 3D multiblock mesh depending on          *
 * input parameters.                                    *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block::Create_Grid(Grid3D_Input_Parameters &Input){

    switch(Input.i_Grid) {
      case GRID_CUBE :
        Create_Grid_Cube(Input);
        break;
      case GRID_CHANNEL :
        Create_Grid_Channel(Input);
        break;
      case GRID_COUETTE :
        Create_Grid_Couette(Input);
        break;
      case GRID_BLUFF_BODY_BURNER :
        Create_Grid_Bluff_Body_Burner(Input);
        break;
      default:
        Create_Grid_Cube(Input);
        break;
    } /* endswitch */

}

/********************************************************
 * Routine: Create_Grid_Cube                            *
 *                                                      *
 * Generates a 3D Cartesian multiblock mesh for a cube. *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block::Create_Grid_Cube(Grid3D_Input_Parameters &Input){

    Grid2D_Quad_Block **Grid2D_Box_XYplane;

    /* Allocate required memory. */

    Allocate(Input.NBlk_Idir, Input.NBlk_Jdir, Input.NBlk_Kdir);

    /* Creat 2D cross-section grids from which the 3D grid
       will be extruded. */
    
    Grid2D_Box_XYplane = Grid_Rectangular_Box(Grid2D_Box_XYplane,
                                              Input.NBlk_Idir, 
                                              Input.NBlk_Jdir,
                                              Input.Box_Width,
                                              Input.Box_Height,
					      ON,
					      Input.Stretching_Type_Idir,
					      Input.Stretching_Type_Jdir,
					      Input.Stretching_Factor_Idir,
					      Input.Stretching_Factor_Jdir,
                                              Input.NCells_Idir,
                                              Input.NCells_Jdir,
					      Input.Nghost);


    /* Create the mesh for each block representing
       the complete grid. */

    for (int kBlk = 0; kBlk <= Input.NBlk_Kdir-1; ++kBlk) {
       for (int jBlk = 0; jBlk <= Input.NBlk_Jdir-1; ++jBlk) {
          for (int iBlk = 0; iBlk <= Input.NBlk_Idir-1; ++iBlk) {

             /* Create the splines defining the north, south,
                east, and west boundaries of the grid. */

             Grid_Blks[iBlk][jBlk][kBlk].Extrude(Grid2D_Box_XYplane[iBlk][jBlk],
                                                 Input.NCells_Kdir,
			                         Input.Stretching_Type_Kdir,
						 Input.Stretching_Factor_Kdir,
                                                 -HALF*Input.Box_Length+
                                                 (double(kBlk)/double(Input.NBlk_Kdir))*Input.Box_Length,
                                                 -HALF*Input.Box_Length+
                                                 (double(kBlk+1)/double(Input.NBlk_Kdir))*Input.Box_Length);

             /* Assign top and bottom boundary conditions. */

             for (int j = Grid_Blks[iBlk][jBlk][kBlk].JCl-Grid_Blks[iBlk][jBlk][kBlk].Nghost; 
                      j <= Grid_Blks[iBlk][jBlk][kBlk].JCu+Grid_Blks[iBlk][jBlk][kBlk].Nghost; ++j) {
                for (int i = Grid_Blks[iBlk][jBlk][kBlk].ICl-Grid_Blks[iBlk][jBlk][kBlk].Nghost; 
                         i <= Grid_Blks[iBlk][jBlk][kBlk].ICu+Grid_Blks[iBlk][jBlk][kBlk].Nghost; ++i) {
		   if (kBlk == Input.NBlk_Kdir-1) {
                      Grid_Blks[iBlk][jBlk][kBlk].BCtypeT[i][j] = BC_REFLECTION;
                   } else {
                      Grid_Blks[iBlk][jBlk][kBlk].BCtypeT[i][j] = BC_NONE;
                   } /* endif */
                   if (kBlk == 0) {
                      Grid_Blks[iBlk][jBlk][kBlk].BCtypeB[i][j] = BC_REFLECTION;
                   } else {
                      Grid_Blks[iBlk][jBlk][kBlk].BCtypeB[i][j] = BC_NONE;
                   } /* endif */
                } /* endfor */
             } /* endfor */

	  } /* endfor */
       } /* endfor */
    } /* endfor */

    /* Deallocate 2D grid. */

    Grid2D_Box_XYplane = Deallocate_Multi_Block_Grid(Grid2D_Box_XYplane,
                                                     Input.NBlk_Idir, 
                                                     Input.NBlk_Jdir);

}

/********************************************************
 * Routine: Create_Grid_Channel                         *
 *                                                      *
 * Generates a 3D Cartesian multiblock mesh for a       *
 * channel flow.                                        *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block::Create_Grid_Channel(Grid3D_Input_Parameters &Input){

    Grid2D_Quad_Block **Grid2D_Box_XYplane;
   
}

/********************************************************
 * Routine: Create_Grid_Couette                         *
 *                                                      *
 * Generates a 3D Cartesian multiblock mesh for a       *
 * couette flow.                                        *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block::Create_Grid_Couette(Grid3D_Input_Parameters &Input){
      
}

/********************************************************
 * Routine: Create_Grid_Bluff_Body_Burner               *
 *                                                      *
 * Generates a 3D multiblock mesh for TNF bluff body    *
 * burner.                                              *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block::Create_Grid_Bluff_Body_Burner(Grid3D_Input_Parameters &Input){

}

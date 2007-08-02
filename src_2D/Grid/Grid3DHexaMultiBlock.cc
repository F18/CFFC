/* Grid3DHexaMultiBlock.cc:  Multi-block subroutines for 
                             3D hexahedral block grid class. */

/* Include 3D quadrilateral block grid type header file. */

#ifndef _GRID3D_HEXA_BLOCK_INCLUDED
#include "Grid3DHexa.h"
#endif // _GRID3D_HEXA_BLOCK_INCLUDED

/*************************************************************************
 * Grid3D_Hexa_Block -- External subroutines for 3D array of grid blocks.*
 *************************************************************************/

/********************************************************
 * Routine: Allocate_Multi_Block_Grid                   *
 *                                                      *
 * Allocate memory for a 3D array of 3D hexahedral      *
 * multi-block grids.                                   *
 *                                                      *
 ********************************************************/
Grid3D_Hexa_Block***  Allocate_Multi_Block_Grid(Grid3D_Hexa_Block ***Grid_ptr,
						const int Number_of_Blocks_Idir,
						const int Number_of_Blocks_Jdir,
						const int Number_of_Blocks_Kdir) {

    
    /* Allocate memory. */
    
    Grid_ptr = new Grid3D_Hexa_Block**[Number_of_Blocks_Idir];
    for (int i = 0 ; i <= Number_of_Blocks_Idir-1 ; ++i ) {
      Grid_ptr[i] = new Grid3D_Hexa_Block*[Number_of_Blocks_Jdir];
      for (int j = 0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
	Grid_ptr[i][j] = new Grid3D_Hexa_Block[Number_of_Blocks_Kdir];
      }  /* endfor */
    }/* endfor */
    
    /* Return memory location. */
    return(Grid_ptr);
    
}

/********************************************************
 * Routine: Deallocate_Multi_Block_Grid                 *
 *                                                      *
 * Deallocate memory for a 3D array of 3D quadrilateral *
 * multi-block grids.                                   *
 *                                                      *
 ********************************************************/
Grid3D_Hexa_Block*** Deallocate_Multi_Block_Grid(Grid3D_Hexa_Block ***Grid_ptr,
						 const int Number_of_Blocks_Idir,
						 const int Number_of_Blocks_Jdir,
						 const int Number_of_Blocks_Kdir) {

 
    /* Deallocate memory. */

    for (int i = 0 ; i <= Number_of_Blocks_Idir-1 ; ++i ) {
       for ( int j=0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
	 for ( int k = Number_of_Blocks_Kdir-1 ; k >= 0; --k ) {
          if (Grid_ptr[i][j][k].Node != NULL && Grid_ptr[i][j][k].Cell != NULL) { 
             Grid_ptr[i][j][k].deallocate();
          } else if (Grid_ptr[i][j][k].Node != NULL) {
             Grid_ptr[i][j][k].deallocateNodes();
          } else if (Grid_ptr[i][j][k].Cell != NULL) {
             Grid_ptr[i][j][k].deallocateCells();
          } /* endif */
	 }/* end for */
	 delete []Grid_ptr[i][j];
	 Grid_ptr[i][j]=NULL;
       }  /* endfor */
       
       delete []Grid_ptr[i];
       Grid_ptr[i] = NULL;
    }  /* endfor */
    delete []Grid_ptr;
    Grid_ptr = NULL;

    /* Return memory location. */

    return(Grid_ptr);

}



/********************************************************
 * Routine: Write_Multi_Block_Grid                      *
 *                                                      *
 * Writes a 3D array of 3D hexahedral    multi-block    *
 * grids to the specified output stream for retrieval   *
 * and re-use purposes.                                 *
 *                                                      *
 ********************************************************/
void Write_Multi_Block_Grid(Grid3D_Hexa_Block ***Grid_ptr,
			    const int Number_of_Blocks_Idir,
		            const int Number_of_Blocks_Jdir,
			    const int Number_of_Blocks_Kdir,
                            ostream &Out_File) {
 
  Out_File << Number_of_Blocks_Idir << " " 
	   << Number_of_Blocks_Jdir << " "
	   << Number_of_Blocks_Kdir << "\n";
  
  for (int  k = 0 ; k <= Number_of_Blocks_Kdir-1 ; ++k ) {
    for (int  j = 0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
      for (int  i = 0 ; i <= Number_of_Blocks_Idir-1 ; ++i ) {
	Write_Hexa_Block(Grid_ptr[i][j][k], Out_File);
      }  /* endfor */
    }  /* endfor */
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
void Output_Tecplot(Grid3D_Hexa_Block ***Grid_ptr,
		    const int Number_of_Blocks_Idir,
		    const int Number_of_Blocks_Jdir,
		    const int Number_of_Blocks_Kdir,
	            ostream &Out_File) {

    int block_number, i_output_title;

    block_number = 0;
    i_output_title = 1;
    
    for (int  k = 0 ; k <= Number_of_Blocks_Kdir-1 ; ++k ) {
      for (int  j = 0; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
	for (int i = 0 ; i <= Number_of_Blocks_Idir-1 ; ++i ) {
          if (Grid_ptr[i][j][k].Node != NULL) {
	    Output_Tecplot(Grid_ptr[i][j][k],
			   block_number,
			   i_output_title,
			   Out_File);
	    block_number = block_number + 1;
             if (i_output_title) i_output_title = 0;
          } /* endif */
	}  /* endfor */
      }  /* endfor */
    }/* endfor */

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
void Output_Nodes_Tecplot(Grid3D_Hexa_Block ***Grid_ptr,
                          const int Number_of_Blocks_Idir,
                          const int Number_of_Blocks_Jdir,
			  const int Number_of_Blocks_Kdir,
	                  ostream &Out_File) {
  int  block_number = 0;
  int  i_output_title = 1;
  
  for (int  k = 0 ; k <= Number_of_Blocks_Kdir-1 ; ++k ) {
    for (int  j = 0; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
      for (int i = 0 ; i <= Number_of_Blocks_Idir-1 ; ++i ) {
	if (Grid_ptr[i][j][k].Node != NULL) {
	  Output_Nodes_Tecplot(Grid_ptr[i][j][k],
			       block_number,
			       i_output_title,
			       Out_File);
	  block_number = block_number + 1;
	  if (i_output_title) i_output_title = 0;
	} /* endif */
      }  /* endfor */
    }  /* endfor */
  }/* endfor */ 
}

/********************************************************
 * Routine: Output_Cells_Tecplot                        *
 *                                                      *
 * Writes the cells of a 2D array of 2D quadrilateral   *
 * multi-block grids to the specified output stream in  *
 * a format suitable for plotting the grid with         *
 * TECPLOT.                                             *
 *                                                      *
 ********************************************************/
void Output_Cells_Tecplot(Grid3D_Hexa_Block ***Grid_ptr,
                          const int Number_of_Blocks_Idir,
                          const int Number_of_Blocks_Jdir,
			  const int Number_of_Blocks_Kdir,
 	                  ostream &Out_File) {

    int block_number, i_output_title;

    block_number = 0;
    i_output_title = 1;
  for (int  k = 0 ; k <= Number_of_Blocks_Kdir-1 ; ++k ) {
    for (int  j = 0; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
      for (int i = 0 ; i <= Number_of_Blocks_Idir-1 ; ++i ) {
	if (Grid_ptr[i][j][k].Node != NULL) {
	  Output_Cells_Tecplot(Grid_ptr[i][j][k],
			       block_number,
			       i_output_title,
			       Out_File);
	  block_number = block_number + 1;
	  if (i_output_title) i_output_title = 0;
	} /* endif */
      }  /* endfor */
    }  /* endfor */
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
void Output_Gnuplot(Grid3D_Hexa_Block ***Grid_ptr,
		    const int Number_of_Blocks_Idir,
		    const int Number_of_Blocks_Jdir,
		    const int Number_of_Blocks_Kdir,
	            ostream &Out_File) {
  
  int  block_number, i_output_title;

    block_number = 0;
    i_output_title = 1;
    
    for (int k = 0 ; k <= Number_of_Blocks_Kdir-1 ; ++k ) {
      for (int j = 0; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
	for(int i =0 ; i <= Number_of_Blocks_Kdir-1 ; ++i ) {
	  if (Grid_ptr[i][j][k].Node != NULL) {
	    Output_Gnuplot(Grid_ptr[i][j][k],
			  block_number,
			  i_output_title,
			  Out_File);
	   block_number = block_number + 1;
	   if (i_output_title) i_output_title = 0;
	  } /* endif */
	}  /* endfor */
      }  /* endfor */
    } /* endfor */
    
}


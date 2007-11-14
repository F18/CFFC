/* Grid3DHexaMultiBlock.cc:  Member functions for 
                             3D hexahedral multiblock grid classes. */

/* Include 3D hexahedral block grid header file. */

#ifndef GRID3D_HEXA_MULTIBLOCK_INCLUDED
#include "Grid3DHexaMultiBlock.h"
#endif //_GRID3D_HEXA_MULTIBLOCK_INCLUDED

/********************************************************
 * Routine: Allocate                                    *
 *                                                      *
 * Allocate memory for a 1D array of 3D hexahedral      *
 * multi-block grids.                                   *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block_List::Allocate(const int Ni, 
                                            const int Nj, 
                                            const int Nk) {

   if (Ni >= 1 && Nj >= 1 && Nk >= 1 && !Allocated) {
      NBlk_Idir = Ni; 
      NBlk_Jdir = Nj; 
      NBlk_Kdir = Nk; 
      NBlk = Ni*Nj*Nk;

      Grid_Blks = new Grid3D_Hexa_Block[NBlk];

      neighT = new int[NBlk];
      neighB = new int[NBlk];
      neighN = new int[NBlk];
      neighS = new int[NBlk];
      neighE = new int[NBlk];
      neighW = new int[NBlk]; 
      neighNW = new int[NBlk];
      neighNE = new int[NBlk];
      neighSE = new int[NBlk];
      neighSW = new int[NBlk];
      neighTN = new int[NBlk];
      neighTS = new int[NBlk];
      neighTE = new int[NBlk];
      neighTW = new int[NBlk];  
      neighTNW = new int[NBlk];
      neighTSW = new int[NBlk];
      neighTNE = new int[NBlk];
      neighTSE = new int[NBlk];  
      neighBN = new int[NBlk];
      neighBS = new int[NBlk];
      neighBE = new int[NBlk];
      neighBW = new int[NBlk];  
      neighBNW = new int[NBlk];
      neighBSW = new int[NBlk];
      neighBNE = new int[NBlk];
      neighBSE = new int[NBlk];

      neighT_ctm = new Mesh_Orientation_Matrix[NBlk];
      neighB_ctm = new Mesh_Orientation_Matrix[NBlk];
      neighN_ctm = new Mesh_Orientation_Matrix[NBlk];
      neighS_ctm = new Mesh_Orientation_Matrix[NBlk];
      neighE_ctm = new Mesh_Orientation_Matrix[NBlk];
      neighW_ctm = new Mesh_Orientation_Matrix[NBlk]; 
      neighNW_ctm = new Mesh_Orientation_Matrix[NBlk];
      neighNE_ctm = new Mesh_Orientation_Matrix[NBlk];
      neighSE_ctm = new Mesh_Orientation_Matrix[NBlk];
      neighSW_ctm = new Mesh_Orientation_Matrix[NBlk];
      neighTN_ctm = new Mesh_Orientation_Matrix[NBlk];
      neighTS_ctm = new Mesh_Orientation_Matrix[NBlk];
      neighTE_ctm = new Mesh_Orientation_Matrix[NBlk];
      neighTW_ctm = new Mesh_Orientation_Matrix[NBlk];  
      neighTNW_ctm = new Mesh_Orientation_Matrix[NBlk];
      neighTSW_ctm = new Mesh_Orientation_Matrix[NBlk];
      neighTNE_ctm = new Mesh_Orientation_Matrix[NBlk];
      neighTSE_ctm = new Mesh_Orientation_Matrix[NBlk];  
      neighBN_ctm = new Mesh_Orientation_Matrix[NBlk];
      neighBS_ctm = new Mesh_Orientation_Matrix[NBlk];
      neighBE_ctm = new Mesh_Orientation_Matrix[NBlk];
      neighBW_ctm = new Mesh_Orientation_Matrix[NBlk];  
      neighBNW_ctm = new Mesh_Orientation_Matrix[NBlk];
      neighBSW_ctm = new Mesh_Orientation_Matrix[NBlk];
      neighBNE_ctm = new Mesh_Orientation_Matrix[NBlk];
      neighBSE_ctm = new Mesh_Orientation_Matrix[NBlk];

      Allocated = 1;

      for (int  i = 0 ; i < NBlk ; ++i ) {
         neighT[i] = GRID3D_NO_NEIGHBOUR;
         neighB[i] = GRID3D_NO_NEIGHBOUR;
         neighN[i] = GRID3D_NO_NEIGHBOUR;
         neighS[i] = GRID3D_NO_NEIGHBOUR;
         neighE[i] = GRID3D_NO_NEIGHBOUR;
         neighW[i] = GRID3D_NO_NEIGHBOUR; 
         neighNW[i] = GRID3D_NO_NEIGHBOUR;
         neighNE[i] = GRID3D_NO_NEIGHBOUR;
         neighSE[i] = GRID3D_NO_NEIGHBOUR;
         neighSW[i] = GRID3D_NO_NEIGHBOUR;
         neighTN[i] = GRID3D_NO_NEIGHBOUR;
         neighTS[i] = GRID3D_NO_NEIGHBOUR;
         neighTE[i] = GRID3D_NO_NEIGHBOUR;
         neighTW[i] = GRID3D_NO_NEIGHBOUR;  
         neighTNW[i] = GRID3D_NO_NEIGHBOUR;
         neighTSW[i] = GRID3D_NO_NEIGHBOUR;
         neighTNE[i] = GRID3D_NO_NEIGHBOUR;
         neighTSE[i] = GRID3D_NO_NEIGHBOUR;  
         neighBN[i] = GRID3D_NO_NEIGHBOUR;
         neighBS[i] = GRID3D_NO_NEIGHBOUR;
         neighBE[i] = GRID3D_NO_NEIGHBOUR;
         neighBW[i] = GRID3D_NO_NEIGHBOUR;  
         neighBNW[i] = GRID3D_NO_NEIGHBOUR;
         neighBSW[i] = GRID3D_NO_NEIGHBOUR;
         neighBNE[i] = GRID3D_NO_NEIGHBOUR;
         neighBSE[i] = GRID3D_NO_NEIGHBOUR;
      } /* endfor */
   } /* endif */

}

/********************************************************
 * Routine: Allocate                                    *
 *                                                      *
 * Allocate memory for a 1D array of 3D hexahedral      *
 * multi-block grids.                                   *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block_List::Allocate(const int N) {

   if (N >= 1 && !Allocated) {
      NBlk_Idir = N; 
      NBlk_Jdir = 1; 
      NBlk_Kdir = 1; 
      NBlk = N;

      Grid_Blks = new Grid3D_Hexa_Block[NBlk];

      neighT = new int[NBlk];
      neighB = new int[NBlk];
      neighN = new int[NBlk];
      neighS = new int[NBlk];
      neighE = new int[NBlk];
      neighW = new int[NBlk]; 
      neighNW = new int[NBlk];
      neighNE = new int[NBlk];
      neighSE = new int[NBlk];
      neighSW = new int[NBlk];
      neighTN = new int[NBlk];
      neighTS = new int[NBlk];
      neighTE = new int[NBlk];
      neighTW = new int[NBlk];  
      neighTNW = new int[NBlk];
      neighTSW = new int[NBlk];
      neighTNE = new int[NBlk];
      neighTSE = new int[NBlk];  
      neighBN = new int[NBlk];
      neighBS = new int[NBlk];
      neighBE = new int[NBlk];
      neighBW = new int[NBlk];  
      neighBNW = new int[NBlk];
      neighBSW = new int[NBlk];
      neighBNE = new int[NBlk];
      neighBSE = new int[NBlk];

      neighT_ctm = new Mesh_Orientation_Matrix[NBlk];
      neighB_ctm = new Mesh_Orientation_Matrix[NBlk];
      neighN_ctm = new Mesh_Orientation_Matrix[NBlk];
      neighS_ctm = new Mesh_Orientation_Matrix[NBlk];
      neighE_ctm = new Mesh_Orientation_Matrix[NBlk];
      neighW_ctm = new Mesh_Orientation_Matrix[NBlk]; 
      neighNW_ctm = new Mesh_Orientation_Matrix[NBlk];
      neighNE_ctm = new Mesh_Orientation_Matrix[NBlk];
      neighSE_ctm = new Mesh_Orientation_Matrix[NBlk];
      neighSW_ctm = new Mesh_Orientation_Matrix[NBlk];
      neighTN_ctm = new Mesh_Orientation_Matrix[NBlk];
      neighTS_ctm = new Mesh_Orientation_Matrix[NBlk];
      neighTE_ctm = new Mesh_Orientation_Matrix[NBlk];
      neighTW_ctm = new Mesh_Orientation_Matrix[NBlk];  
      neighTNW_ctm = new Mesh_Orientation_Matrix[NBlk];
      neighTSW_ctm = new Mesh_Orientation_Matrix[NBlk];
      neighTNE_ctm = new Mesh_Orientation_Matrix[NBlk];
      neighTSE_ctm = new Mesh_Orientation_Matrix[NBlk];  
      neighBN_ctm = new Mesh_Orientation_Matrix[NBlk];
      neighBS_ctm = new Mesh_Orientation_Matrix[NBlk];
      neighBE_ctm = new Mesh_Orientation_Matrix[NBlk];
      neighBW_ctm = new Mesh_Orientation_Matrix[NBlk];  
      neighBNW_ctm = new Mesh_Orientation_Matrix[NBlk];
      neighBSW_ctm = new Mesh_Orientation_Matrix[NBlk];
      neighBNE_ctm = new Mesh_Orientation_Matrix[NBlk];
      neighBSE_ctm = new Mesh_Orientation_Matrix[NBlk];

      Allocated = 1;

      for (int  i = 0 ; i < NBlk ; ++i ) {
         neighT[i] = GRID3D_NO_NEIGHBOUR;
         neighB[i] = GRID3D_NO_NEIGHBOUR;
         neighN[i] = GRID3D_NO_NEIGHBOUR;
         neighS[i] = GRID3D_NO_NEIGHBOUR;
         neighE[i] = GRID3D_NO_NEIGHBOUR;
         neighW[i] = GRID3D_NO_NEIGHBOUR; 
         neighNW[i] = GRID3D_NO_NEIGHBOUR;
         neighNE[i] = GRID3D_NO_NEIGHBOUR;
         neighSE[i] = GRID3D_NO_NEIGHBOUR;
         neighSW[i] = GRID3D_NO_NEIGHBOUR;
         neighTN[i] = GRID3D_NO_NEIGHBOUR;
         neighTS[i] = GRID3D_NO_NEIGHBOUR;
         neighTE[i] = GRID3D_NO_NEIGHBOUR;
         neighTW[i] = GRID3D_NO_NEIGHBOUR;  
         neighTNW[i] = GRID3D_NO_NEIGHBOUR;
         neighTSW[i] = GRID3D_NO_NEIGHBOUR;
         neighTNE[i] = GRID3D_NO_NEIGHBOUR;
         neighTSE[i] = GRID3D_NO_NEIGHBOUR;  
         neighBN[i] = GRID3D_NO_NEIGHBOUR;
         neighBS[i] = GRID3D_NO_NEIGHBOUR;
         neighBE[i] = GRID3D_NO_NEIGHBOUR;
         neighBW[i] = GRID3D_NO_NEIGHBOUR;  
         neighBNW[i] = GRID3D_NO_NEIGHBOUR;
         neighBSW[i] = GRID3D_NO_NEIGHBOUR;
         neighBNE[i] = GRID3D_NO_NEIGHBOUR;
         neighBSE[i] = GRID3D_NO_NEIGHBOUR;
      } /* endfor */
   } /* endif */

}

/********************************************************
 * Routine: Deallocate                                  *
 *                                                      *
 * Deallocate memory for a 1D array of 3D hexahedral    *
 * multi-block grids.                                   *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block_List::Deallocate(void) {

   if (NBlk >= 1 && Allocated) {
       delete []Grid_Blks;
       Grid_Blks = NULL;

       delete []neighT;
       delete []neighB;
       delete []neighN;
       delete []neighS;
       delete []neighE;
       delete []neighW; 
       delete []neighNW;
       delete []neighNE;
       delete []neighSE;
       delete []neighSW;
       delete []neighTN;
       delete []neighTS;
       delete []neighTE;
       delete []neighTW;  
       delete []neighTNW;
       delete []neighTSW;
       delete []neighTNE;
       delete []neighTSE;  
       delete []neighBN;
       delete []neighBS;
       delete []neighBE;
       delete []neighBW;  
       delete []neighBNW;
       delete []neighBSW;
       delete []neighBNE;
       delete []neighBSE;
       neighT = NULL;
       neighB = NULL;
       neighN = NULL;
       neighS = NULL;
       neighE = NULL;
       neighW = NULL; 
       neighNW = NULL;
       neighNE = NULL;
       neighSE = NULL;
       neighSW = NULL;
       neighTN = NULL;
       neighTS = NULL;
       neighTE = NULL;
       neighTW = NULL;  
       neighTNW = NULL;
       neighTSW = NULL;
       neighTNE = NULL;
       neighTSE = NULL;  
       neighBN = NULL;
       neighBS = NULL;
       neighBE = NULL;
       neighBW = NULL;  
       neighBNW = NULL;
       neighBSW = NULL;
       neighBNE = NULL;
       neighBSE = NULL;

       delete []neighT_ctm;
       delete []neighB_ctm;
       delete []neighN_ctm;
       delete []neighS_ctm;
       delete []neighE_ctm;
       delete []neighW_ctm; 
       delete []neighNW_ctm;
       delete []neighNE_ctm;
       delete []neighSE_ctm;
       delete []neighSW_ctm;
       delete []neighTN_ctm;
       delete []neighTS_ctm;
       delete []neighTE_ctm;
       delete []neighTW_ctm;  
       delete []neighTNW_ctm;
       delete []neighTSW_ctm;
       delete []neighTNE_ctm;
       delete []neighTSE_ctm;  
       delete []neighBN_ctm;
       delete []neighBS_ctm;
       delete []neighBE_ctm;
       delete []neighBW_ctm;  
       delete []neighBNW_ctm;
       delete []neighBSW_ctm;
       delete []neighBNE_ctm;
       delete []neighBSE_ctm;
       neighT_ctm = NULL;
       neighB_ctm = NULL;
       neighN_ctm = NULL;
       neighS_ctm = NULL;
       neighE_ctm = NULL;
       neighW_ctm = NULL; 
       neighNW_ctm = NULL;
       neighNE_ctm = NULL;
       neighSE_ctm = NULL;
       neighSW_ctm = NULL;
       neighTN_ctm = NULL;
       neighTS_ctm = NULL;
       neighTE_ctm = NULL;
       neighTW_ctm = NULL;  
       neighTNW_ctm = NULL;
       neighTSW_ctm = NULL;
       neighTNE_ctm = NULL;
       neighTSE_ctm = NULL;  
       neighBN_ctm = NULL;
       neighBS_ctm = NULL;
       neighBE_ctm = NULL;
       neighBW_ctm = NULL;  
       neighBNW_ctm = NULL;
       neighBSW_ctm = NULL;
       neighBNE_ctm = NULL;
       neighBSE_ctm = NULL;

       NBlk_Idir = 0; 
       NBlk_Jdir = 0; 
       NBlk_Kdir = 0;
       NBlk = 0;

       Allocated = 0;
   } /* endif */

}

/********************************************************
 * Routine: Copy                                        *
 *                                                      *
 * Make a copy of multiblock hexahedral grid Grid2.     *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block_List::Copy(Grid3D_Hexa_Multi_Block_List &Grid2) {

  if (Grid2.Allocated) {

    /* Ensure multiblock grid arrays have same dimensions. */

    if (Allocated && (NBlk      != Grid2.NBlk      ||
                      NBlk_Idir != Grid2.NBlk_Idir ||
                      NBlk_Jdir != Grid2.NBlk_Jdir ||
                      NBlk_Kdir != Grid2.NBlk_Kdir) ) {
      Deallocate();
      Allocate(Grid2.NBlk_Idir, Grid2.NBlk_Jdir, Grid2.NBlk_Kdir);
    } else if (!Allocated) {
      Allocate(Grid2.NBlk_Idir, Grid2.NBlk_Jdir, Grid2.NBlk_Kdir);
    } /* endif */

    /* Copy each grid block. */

    for (int  i = 0 ; i < NBlk ; ++i ) {
       if (Grid2.Grid_Blks[i].Allocated) 
          Grid_Blks[i].Copy(Grid2.Grid_Blks[i]);
    } /* endfor */

  } /* endif */

}

/********************************************************
 * Routine: Broadcast                                   *
 *                                                      *
 * Broadcast multiblock hexahedral grid.                *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block_List::Broadcast(void) {

#ifdef _MPI_VERSION
  int n, ni, nj, nk, grid_allocated;

  /* Broadcast the number of grid blocks. */

  if (CFFC_Primary_MPI_Processor()) {
     n = NBlk;
     ni = NBlk_Idir;
     nj = NBlk_Jdir;
     nk = NBlk_Kdir;
     grid_allocated = Allocated;
  } /* endif */

  MPI::COMM_WORLD.Bcast(&n, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&ni, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&nj, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&nk, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&grid_allocated, 1, MPI::INT, 0);

  /* On non-primary MPI processors, allocate (re-allocate) 
     memory for the grid blocks as necessary. */

  if (!CFFC_Primary_MPI_Processor()) {
     if (grid_allocated && 
         (NBlk != n ||
          NBlk_Idir != ni || 
          NBlk_Jdir != nj || 
          NBlk_Kdir != nk) ) { 
        if (Allocated) { 
          Deallocate();
        } /* endif */
        Allocate(ni, nj, nk);
     } /* endif */
  } /* endif */

  /* Broadcast each of the blocks in the multiblock mesh. */

  if (Allocated) {
    for (int  i = 0 ; i < NBlk ; ++i ) {
       Grid_Blks[i].Broadcast();
    } /* endfor */
  } /* endif */
#endif

}

/********************************************************
 * Routine: Output                                      *
 *                                                      *
 * Writes a 1D array of 3D hexahedral multi-block       *
 * grids to the specified output stream for retrieval   *
 * and re-use purposes.                                 *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block_List::Output(ostream &Out_File) {
   
  Out_File << NBlk << " "
           << NBlk_Idir << " "
           << NBlk_Jdir << " "
           << NBlk_Kdir << "\n";
  
  for (int  i = 0 ; i < NBlk ; ++i ) {
     if (Grid_Blks[i].Allocated) {
        Out_File << setprecision(14) << Grid_Blks[i] << setprecision(6);
     } /* endif */
  } /* endfor */

}

/********************************************************
 * Routine: Output_Tecplot                              *
 *                                                      *
 * Writes the nodes of a 1D array of 3D hexahedral      *
 * multi-block grids to the specified output stream in  *
 * a format suitable for plotting the grid with         *
 * TECPLOT.                                             *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block_List::Output_Tecplot(ostream &Out_File) {

    int i_output_title;
    i_output_title = 1;
    
    for (int i = 0 ; i < NBlk ; ++i ) {
       if (Grid_Blks[i].Allocated) {
          Grid_Blks[i].Output_Tecplot(i,
                                      i_output_title,
                                      Out_File);
          if (i_output_title) i_output_title = 0;
       } /* endif */
    } /* endfor */

}

/********************************************************
 * Routine: Output_Nodes_Tecplot                        *
 *                                                      *
 * Writes the nodes of a 1D array of 3D hexahedral      *
 * multi-block grids to the specified output stream in  *
 * a format suitable for plotting the grid with         *
 * TECPLOT.  Include boundary nodes.                    *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block_List::Output_Nodes_Tecplot(ostream &Out_File) {

    int i_output_title;
    i_output_title = 1;
  
    for (int i = 0 ; i < NBlk ; ++i ) {
       if (Grid_Blks[i].Allocated) {
          Grid_Blks[i].Output_Nodes_Tecplot(i,
                                            i_output_title,
                                            Out_File);
          if (i_output_title) i_output_title = 0;
       } /* endif */
    } /* endfor */ 

}

/********************************************************
 * Routine: Output_Cells_Tecplot                        *
 *                                                      *
 * Writes the cells of a 1D array of 3D hexahedral      *
 * multi-block grids to the specified output stream in  *
 * a format suitable for plotting the grid with         *
 * TECPLOT.                                             *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block_List::Output_Cells_Tecplot(ostream &Out_File) {

    int i_output_title;
    i_output_title = 1;

    for (int i = 0 ; i < NBlk ; ++i ) {
       if (Grid_Blks[i].Allocated) {
           Grid_Blks[i].Output_Cells_Tecplot(i,
                                             i_output_title,
                                             Out_File);
           if (i_output_title) i_output_title = 0;
	} /* endif */
    }/* endfor */ 

}

/********************************************************
 * Routine: Output_Gnuplot                              *
 *                                                      *
 * Writes the nodes of a 1D array of 3D hexahedral      *
 * multi-block grids to the specified output stream in  *
 * a format suitable for plotting the grid with         *
 * GNUPLOT.                                             *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block_List::Output_Gnuplot(ostream &Out_File) {
  
    int i_output_title;
    i_output_title = 1;
    
    for (int i =0 ; i < NBlk ; ++i ) {
       if (Grid_Blks[i].Allocated) {
	 Grid_Blks[i].Output_Gnuplot(i,
                                     i_output_title,
		                     Out_File);
          if (i_output_title) i_output_title = 0;
       } /* endif */
    } /* endfor */
    
}

/********************************************************
 * Routine: Create_Grid                                 *
 *                                                      *
 * Generates a 3D multiblock mesh depending on          *
 * input parameters.                                    *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block_List::Create_Grid(Grid3D_Input_Parameters &Input) {

    switch(Input.i_Grid) {
      case GRID_CUBE :
        Create_Grid_Cube(Input);
        break;

      case GRID_CHANNEL_XDIR :
      case GRID_CHANNEL_YDIR:
      case GRID_CHANNEL_ZDIR:
        Create_Grid_Channel(Input);
        break;

      case GRID_COUETTE_XDIR :
      case GRID_COUETTE_YDIR:
      case GRID_COUETTE_ZDIR:
        Create_Grid_Couette(Input);
        break;

      case GRID_PIPE :
        Create_Grid_Pipe(Input);
        break;

      case GRID_BLUFF_BODY_BURNER :
        Create_Grid_Bluff_Body_Burner(Input);
        break;

      case GRID_ICEMCFD :
        //Create_Grid_ICEMCFD(Input);
        break;

      case GRID_TURBULENT_PREMIXED_FLAME :
	Create_Grid_Turbulent_Premixed_Flame(Input);
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
void Grid3D_Hexa_Multi_Block_List::Create_Grid_Cube(Grid3D_Input_Parameters &Input) {

    int count_blocks;
    int BC_top, BC_bottom;
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

    count_blocks = 0;

    for (int kBlk = 0; kBlk <= Input.NBlk_Kdir-1; ++kBlk) {
       for (int jBlk = 0; jBlk <= Input.NBlk_Jdir-1; ++jBlk) {
          for (int iBlk = 0; iBlk <= Input.NBlk_Idir-1; ++iBlk) {

             /* Extrude each of the grid blocks from the
                appropriate 2D grid in XY-plane. */

             Grid_Blks[count_blocks].Extrude(Grid2D_Box_XYplane[iBlk][jBlk],
                                             Input.NCells_Kdir,
			                     Input.Stretching_Type_Kdir,
				             Input.Stretching_Factor_Kdir,
                                             -HALF*Input.Box_Length+
                                             (double(kBlk)/double(Input.NBlk_Kdir))*Input.Box_Length,
                                             -HALF*Input.Box_Length+
                                             (double(kBlk+1)/double(Input.NBlk_Kdir))*Input.Box_Length);

             /* Assign top and bottom boundary conditions. */

	     if (kBlk == Input.NBlk_Kdir-1) {
                BC_top = BC_REFLECTION;
             } else {
                BC_top = BC_NONE;
             } /* endif */
             if (kBlk == 0) {
                BC_bottom = BC_REFLECTION;
             } else {
                BC_bottom = BC_NONE;
             } /* endif */

             Grid_Blks[count_blocks].Set_BCs_Zdir(BC_top, BC_bottom);

             /* Update block counter. */

             count_blocks ++;

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
void Grid3D_Hexa_Multi_Block_List::Create_Grid_Channel(Grid3D_Input_Parameters &Input) {

    int count_blocks;
    int BC_east, BC_west, BC_north, BC_south, BC_top, BC_bottom;
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

    count_blocks = 0;

    for (int kBlk = 0; kBlk <= Input.NBlk_Kdir-1; ++kBlk) {
       for (int jBlk = 0; jBlk <= Input.NBlk_Jdir-1; ++jBlk) {
          for (int iBlk = 0; iBlk <= Input.NBlk_Idir-1; ++iBlk) {

             /* Extrude each of the grid blocks from the
                appropriate 2D grid in XY-plane. */

             Grid_Blks[count_blocks].Extrude(Grid2D_Box_XYplane[iBlk][jBlk],
                                             Input.NCells_Kdir,
			                     Input.Stretching_Type_Kdir,
					     Input.Stretching_Factor_Kdir,
                                             (double(kBlk)/double(Input.NBlk_Kdir))*Input.Box_Length,
                                             (double(kBlk+1)/double(Input.NBlk_Kdir))*Input.Box_Length);

             /* Assign top and bottom boundary conditions. */

             if (Input.i_Grid == GRID_CHANNEL_ZDIR) {

   	        if (iBlk == Input.NBlk_Idir-1) {
                   BC_east = BC_CONSTANT_EXTRAPOLATION;
                } else {
                   BC_east = BC_NONE;
                } /* endif */
                if (iBlk == 0) {
                   BC_west = BC_CONSTANT_EXTRAPOLATION;
                } else {
                   BC_west = BC_NONE;
                } /* endif */

	        if (jBlk == Input.NBlk_Jdir-1) {
                   BC_north = BC_WALL_VISCOUS;
                } else {
                   BC_north = BC_NONE;
                } /* endif */
                if (jBlk == 0) {
                   BC_south = BC_WALL_VISCOUS;
                } else {
                   BC_south = BC_NONE;
                } /* endif */

	        if (kBlk == Input.NBlk_Kdir-1) {
                   BC_top = BC_CHANNEL_OUTFLOW;
                } else {
                   BC_top = BC_NONE;
                } /* endif */
                if (kBlk == 0) {
                   BC_bottom = BC_CHANNEL_INFLOW;
                } else {
                   BC_bottom = BC_NONE;
                } /* endif */

             } else if (Input.i_Grid == GRID_CHANNEL_XDIR) {

    	        if (iBlk == Input.NBlk_Idir-1) {
                   BC_east = BC_CHANNEL_OUTFLOW;
                } else {
                   BC_east = BC_NONE;
                } /* endif */
                if (iBlk == 0) {
                   BC_west = BC_CHANNEL_INFLOW;
                } else {
                   BC_west = BC_NONE;
                } /* endif */

	        if (jBlk == Input.NBlk_Jdir-1) {
                   BC_north = BC_WALL_VISCOUS;
                } else {
                   BC_north = BC_NONE;
                } /* endif */
                if (jBlk == 0) {
                   BC_south = BC_WALL_VISCOUS;
                } else {
                   BC_south = BC_NONE;
                } /* endif */

	        if (kBlk == Input.NBlk_Kdir-1) {
                   BC_top = BC_CONSTANT_EXTRAPOLATION;
                } else {
                   BC_top = BC_NONE;
                } /* endif */
                if (kBlk == 0) {
                   BC_bottom = BC_CONSTANT_EXTRAPOLATION;
                } else {
                   BC_bottom = BC_NONE;
                } /* endif */

             } else if (Input.i_Grid == GRID_CHANNEL_YDIR) {

   	        if (iBlk == Input.NBlk_Idir-1) {
                   BC_east = BC_CONSTANT_EXTRAPOLATION;
                } else {
                   BC_east = BC_NONE;
                } /* endif */
                if (iBlk == 0) {
                   BC_west = BC_CONSTANT_EXTRAPOLATION;
                } else {
                   BC_west = BC_NONE;
                } /* endif */

	        if (jBlk == Input.NBlk_Jdir-1) {
                   BC_north = BC_CHANNEL_OUTFLOW;
                } else {
                   BC_north = BC_NONE;
                } /* endif */
                if (jBlk == 0) {
                   BC_south = BC_CHANNEL_INFLOW;
                } else {
                   BC_south = BC_NONE;
                } /* endif */

	        if (kBlk == Input.NBlk_Kdir-1) {
                   BC_top = BC_WALL_VISCOUS;
                } else {
                   BC_top = BC_NONE;
                } /* endif */
                if (kBlk == 0) {
                   BC_bottom = BC_WALL_VISCOUS;
                } else {
                   BC_bottom = BC_NONE;
                } /* endif */
	     } /* endif */

             Grid_Blks[count_blocks].Set_BCs(BC_east, 
                                             BC_west, 
                                             BC_north, 
                                             BC_south, 
                                             BC_top, 
                                             BC_bottom);
 
             /* Update block counter. */

             count_blocks ++;

	  } /* endfor */
       } /* endfor */
    } /* endfor */

    /* Deallocate 2D grid. */

    Grid2D_Box_XYplane = Deallocate_Multi_Block_Grid(Grid2D_Box_XYplane,
                                                     Input.NBlk_Idir, 
                                                     Input.NBlk_Jdir);

}

/********************************************************
 * Routine: Create_Grid_Couette                         *
 *                                                      *
 * Generates a 3D Cartesian multiblock mesh for a       *
 * Couette flow.                                        *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block_List::Create_Grid_Couette(Grid3D_Input_Parameters &Input) {
      
    int count_blocks;
    int BC_east, BC_west, BC_north, BC_south, BC_top, BC_bottom;
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

    count_blocks = 0;

    for (int kBlk = 0; kBlk <= Input.NBlk_Kdir-1; ++kBlk) {
       for (int jBlk = 0; jBlk <= Input.NBlk_Jdir-1; ++jBlk) {
          for (int iBlk = 0; iBlk <= Input.NBlk_Idir-1; ++iBlk) {

             /* Extrude each of the grid blocks from the
                appropriate 2D grid in XY-plane. */

             Grid_Blks[count_blocks].Extrude(Grid2D_Box_XYplane[iBlk][jBlk],
                                             Input.NCells_Kdir,
			                     Input.Stretching_Type_Kdir,
					     Input.Stretching_Factor_Kdir,
                                             (double(kBlk)/double(Input.NBlk_Kdir))*Input.Box_Length,
                                             (double(kBlk+1)/double(Input.NBlk_Kdir))*Input.Box_Length);

             /* Assign top and bottom boundary conditions. */

             if (Input.i_Grid == GRID_COUETTE_ZDIR) {

   	        if (iBlk == Input.NBlk_Idir-1) {
                   BC_east = BC_CONSTANT_EXTRAPOLATION;
                } else {
                   BC_east = BC_NONE;
                } /* endif */
                if (iBlk == 0) {
                   BC_west = BC_CONSTANT_EXTRAPOLATION;
                } else {
                   BC_west = BC_NONE;
                } /* endif */

	        if (jBlk == Input.NBlk_Jdir-1) {
		   BC_north = BC_MOVING_WALL;
                } else {
                   BC_north = BC_NONE;
                } /* endif */
                if (jBlk == 0) {
                   BC_south = BC_NO_SLIP;
                } else {
                   BC_south = BC_NONE;
                } /* endif */

	        if (kBlk == Input.NBlk_Kdir-1) {
                   BC_top = BC_FIXED_PRESSURE;
                } else {
                   BC_top = BC_NONE;
                } /* endif */
                if (kBlk == 0) {
                   BC_bottom = BC_FIXED_PRESSURE;
                } else {
                   BC_bottom = BC_NONE;
                } /* endif */

             } else if (Input.i_Grid == GRID_COUETTE_XDIR) {

   	        if (iBlk == Input.NBlk_Idir-1) {
                   BC_east = BC_FIXED_PRESSURE;
                } else {
                   BC_east = BC_NONE;
                } /* endif */
                if (iBlk == 0) {
                   BC_west = BC_FIXED_PRESSURE;
                } else {
                   BC_west = BC_NONE;
                } /* endif */

	        if (jBlk == Input.NBlk_Jdir-1) {
                   BC_north = BC_MOVING_WALL;
                } else {
                   BC_north = BC_NONE;
                } /* endif */
                if (jBlk == 0) {
                   BC_south = BC_NO_SLIP;
                } else {
                   BC_south = BC_NONE;
                } /* endif */

	        if (kBlk == Input.NBlk_Kdir-1) {
                   BC_top = BC_CONSTANT_EXTRAPOLATION;
                } else {
                   BC_top = BC_NONE;
                } /* endif */
                if (kBlk == 0) {
                   BC_bottom = BC_CONSTANT_EXTRAPOLATION;
                } else {
                   BC_bottom = BC_NONE;
                } /* endif */

             } else if (Input.i_Grid == GRID_COUETTE_YDIR) {

   	        if (iBlk == Input.NBlk_Idir-1) {
                   BC_east = BC_CONSTANT_EXTRAPOLATION;
                } else {
                   BC_east = BC_NONE;
                } /* endif */
                if (iBlk == 0) {
                   BC_west = BC_CONSTANT_EXTRAPOLATION;
                } else {
                   BC_west = BC_NONE;
                } /* endif */

	        if (jBlk == Input.NBlk_Jdir-1) {
                   BC_north = BC_FIXED_PRESSURE;
                } else {
                   BC_north = BC_NONE;
                } /* endif */
                if (jBlk == 0) {
                   BC_south = BC_FIXED_PRESSURE;
                } else {
                   BC_south = BC_NONE;
                } /* endif */

	        if (kBlk == Input.NBlk_Kdir-1) {
                   BC_top = BC_MOVING_WALL;
                } else {
                   BC_top = BC_NONE;
                } /* endif */
                if (kBlk == 0) {
                   BC_bottom = BC_NO_SLIP;
                } else {
                   BC_bottom = BC_NONE;
                } /* endif */
	     } /* endif */

             Grid_Blks[count_blocks].Set_BCs(BC_east, 
                                             BC_west, 
                                             BC_north, 
                                             BC_south, 
                                             BC_top, 
                                             BC_bottom);
 
             /* Update block counter. */

             count_blocks ++;

	  } /* endfor */
       } /* endfor */
    } /* endfor */

    /* Deallocate 2D grid. */

    Grid2D_Box_XYplane = Deallocate_Multi_Block_Grid(Grid2D_Box_XYplane,
                                                     Input.NBlk_Idir, 
                                                     Input.NBlk_Jdir);

}

/********************************************************
 * Routine: Create_Grid_Pipe                            *
 *                                                      *
 * Generates a 3D multiblock mesh for pipe flows.       *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block_List::Create_Grid_Pipe(Grid3D_Input_Parameters &Input) {

    int numblk_idir_pipe, numblk_jdir_pipe;
    Grid2D_Quad_Block **Grid2D_Pipe_XYplane;

    /* Allocate required memory. */

    Input.NBlk_Idir = 5;
    Input.NBlk_Jdir = 1;
    Input.NBlk_Kdir = 1;

    Allocate(Input.NBlk_Idir, Input.NBlk_Jdir, Input.NBlk_Kdir);

    /* Creat 2D cross-section grids from which the 3D grid
       will be extruded. */

    Grid2D_Pipe_XYplane = Grid_Tube_2D(
                                   Grid2D_Pipe_XYplane,
                                   numblk_idir_pipe,
		                   numblk_jdir_pipe,
                                   Input.Pipe_Radius,
                                   Input.NCells_Idir,
                                   Input.NCells_Jdir,
				   Input.Nghost,
                                   STRETCHING_FCN_MAX_CLUSTERING,
                                   1.25);

    /* Create the mesh for each block representing
       the complete grid. */

    for (int iBlk = 0; iBlk <= Input.NBlk_Idir-1; ++iBlk) {

        /* Extrude each of the grid blocks from the
           appropriate 2D grid in XY-plane. */

        Grid_Blks[iBlk].Extrude(Grid2D_Pipe_XYplane[iBlk][0],
                                Input.NCells_Kdir,
       	  	                Input.Stretching_Type_Kdir,
			        Input.Stretching_Factor_Kdir,
                                ZERO,
                                Input.Pipe_Length);
        if (iBlk == 0) {
           Grid_Blks[iBlk].Set_BCs(BC_NONE,
         	                   BC_NONE,
                                   BC_NONE,
       	                           BC_NONE,
                                   BC_NONE,
                                   BC_DIRICHLET);
        } else {
           Grid_Blks[iBlk].Set_BCs(BC_NONE,
       	                           BC_NONE,
				   BC_WALL_VISCOUS,
       	                           BC_NONE,
                                   BC_NONE,
                                   BC_DIRICHLET);
        } /* endif */

    } /* endfor */

    /* Deallocate 2D grid. */

    Grid2D_Pipe_XYplane = Deallocate_Multi_Block_Grid(
                                   Grid2D_Pipe_XYplane,
                                   numblk_idir_pipe,
		                   numblk_jdir_pipe);

}

/********************************************************
 * Routine: Create_Grid_Bluff_Body_Burner               *
 *                                                      *
 * Generates a 3D multiblock mesh for TNF bluff body    *
 * burner.                                              *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block_List::Create_Grid_Bluff_Body_Burner(Grid3D_Input_Parameters &Input) {

    int numblk_idir_fuel, numblk_jdir_fuel,
        numblk_idir_bluffbody, numblk_jdir_bluffbody,
        numblk_idir_coflow, numblk_jdir_coflow;
    Grid2D_Quad_Block **Grid2D_Fuel_Line_XYplane,
                      **Grid2D_Bluff_Body_Inner_XYplane,
                      **Grid2D_Bluff_Body_Outer_XYplane,
                      **Grid2D_Coflow_XYplane;

    /* Allocate required memory. */

    Input.NBlk_Idir = 42;
    Input.NBlk_Jdir = 1;
    Input.NBlk_Kdir = 1;

    Allocate(Input.NBlk_Idir, Input.NBlk_Jdir, Input.NBlk_Kdir);

    /* Creat 2D cross-section grids from which the 3D grid
       will be extruded. */

    Grid2D_Fuel_Line_XYplane = Grid_Tube_2D(
                                   Grid2D_Fuel_Line_XYplane,
                                   numblk_idir_fuel,
		                   numblk_jdir_fuel,
                                   Input.Radius_Fuel_Line,
                                   Input.NCells_Idir,
                                   Input.NCells_Jdir,
				   Input.Nghost,
                                   STRETCHING_FCN_MAX_CLUSTERING,
                                   1.25);

    Grid2D_Bluff_Body_Inner_XYplane = Grid_Annulus_2D(
                                   Grid2D_Bluff_Body_Inner_XYplane,
                                   numblk_idir_bluffbody,
		                   numblk_jdir_bluffbody,
                                   Input.Radius_Fuel_Line,
 			           HALF*Input.Radius_Bluff_Body,
                                   Input.NCells_Idir,
                                   Input.NCells_Jdir,
				   Input.Nghost,
                                   STRETCHING_FCN_MIN_CLUSTERING,
                                   1.10);

    Grid2D_Bluff_Body_Outer_XYplane = Grid_Annulus_2D(
                                   Grid2D_Bluff_Body_Outer_XYplane,
                                   numblk_idir_bluffbody,
		                   numblk_jdir_bluffbody,
                                   HALF*Input.Radius_Bluff_Body,
 			           Input.Radius_Bluff_Body,
                                   Input.NCells_Idir,
                                   Input.NCells_Jdir,
				   Input.Nghost,
                                   STRETCHING_FCN_MAX_CLUSTERING,
                                   1.10);

    Grid2D_Coflow_XYplane = Grid_Annulus_2D(
                                   Grid2D_Coflow_XYplane,
                                   numblk_idir_coflow,
		                   numblk_jdir_coflow,
 			           Input.Radius_Bluff_Body,
                                   Input.Radius_Coflow_Inlet_Pipe,
                                   Input.NCells_Idir,
                                   Input.NCells_Jdir,
				   Input.Nghost,
                                   STRETCHING_FCN_MIN_CLUSTERING,
                                   1.10);

    /* Create the mesh for each block representing
       the complete grid. */

    for (int iBlk = 0; iBlk <= Input.NBlk_Idir-1; ++iBlk) {

        /* Extrude each of the grid blocks from the
           appropriate 2D grid in XY-plane. */

        if (iBlk <= 4) {
           Grid_Blks[iBlk].Extrude(Grid2D_Fuel_Line_XYplane[iBlk][0],
                                   Input.NCells_Kdir,
			           STRETCHING_FCN_MIN_CLUSTERING,
                                   1.25,
                                   ZERO,
                                   0.25*Input.Length_Combustor_Tube);
           if (iBlk == 0) {
	      Grid_Blks[iBlk].Set_BCs(BC_NONE,
		                      BC_NONE,
                                      BC_NONE,
		                      BC_NONE,
                                      BC_NONE,
                                      BC_DIRICHLET);
	   } else {
              Grid_Blks[iBlk].Set_BCs(BC_NONE,
		                      BC_NONE,
                                      BC_NONE,
		                      BC_NONE,
                                      BC_NONE,
                                      BC_DIRICHLET);
           }  /* endif */

        } else if (iBlk >= 5 && iBlk <= 9) {
           Grid_Blks[iBlk].Extrude(Grid2D_Fuel_Line_XYplane[iBlk-5][0],
                                   Input.NCells_Kdir,
			           STRETCHING_FCN_LINEAR,
                                   ZERO,
                                   0.25*Input.Length_Combustor_Tube,
                                   Input.Length_Combustor_Tube);
           if (iBlk-5 == 0) {
	      Grid_Blks[iBlk].Set_BCs(BC_NONE,
		                      BC_NONE,
                                      BC_NONE,
		                      BC_NONE,
                                      BC_FIXED_PRESSURE,
                                      BC_NONE);
	   } else {
              Grid_Blks[iBlk].Set_BCs(BC_NONE,
		                      BC_NONE,
                                      BC_NONE,
		                      BC_NONE,
                                      BC_FIXED_PRESSURE,
                                      BC_NONE);
           } /* endif */

        } else if (iBlk >= 10 && iBlk <= 13) {
           Grid_Blks[iBlk].Extrude(Grid2D_Bluff_Body_Inner_XYplane[iBlk-10][0],
                                   Input.NCells_Kdir,
			           STRETCHING_FCN_MIN_CLUSTERING,
                                   1.25,
                                   ZERO,
                                   0.25*Input.Length_Combustor_Tube);
	   Grid_Blks[iBlk].Set_BCs(BC_NONE,
		                   BC_NONE,
                                   BC_NONE,
		                   BC_NONE,
                                   BC_NONE,
                                   BC_WALL_VISCOUS);

        } else if (iBlk >= 14 && iBlk <= 17) {
           Grid_Blks[iBlk].Extrude(Grid2D_Bluff_Body_Inner_XYplane[iBlk-14][0],
                                   Input.NCells_Kdir,
			           STRETCHING_FCN_LINEAR,
                                   ZERO,
                                   0.25*Input.Length_Combustor_Tube,
                                   Input.Length_Combustor_Tube);
	   Grid_Blks[iBlk].Set_BCs(BC_NONE,
		                   BC_NONE,
                                   BC_NONE,
		                   BC_NONE,
                                   BC_FIXED_PRESSURE,
                                   BC_NONE);

        } else if (iBlk >= 18 && iBlk <= 21) {
           Grid_Blks[iBlk].Extrude(Grid2D_Bluff_Body_Outer_XYplane[iBlk-18][0],
                                   Input.NCells_Kdir,
			           STRETCHING_FCN_MIN_CLUSTERING,
                                   1.25,
                                   ZERO,
                                   0.25*Input.Length_Combustor_Tube);
	   Grid_Blks[iBlk].Set_BCs(BC_NONE,
		                   BC_NONE,
                                   BC_NONE,
		                   BC_NONE,
                                   BC_NONE,
                                   BC_WALL_VISCOUS);

        } else if (iBlk >= 22 && iBlk <= 25) {
           Grid_Blks[iBlk].Extrude(Grid2D_Bluff_Body_Outer_XYplane[iBlk-22][0],
                                   Input.NCells_Kdir,
			           STRETCHING_FCN_LINEAR,
                                   ONE,
                                   0.25*Input.Length_Combustor_Tube,
                                   Input.Length_Combustor_Tube);
	   Grid_Blks[iBlk].Set_BCs(BC_NONE,
		                   BC_NONE,
                                   BC_NONE,
		                   BC_NONE,
                                   BC_FIXED_PRESSURE,
                                   BC_NONE);

        } else if (iBlk >= 26 && iBlk <= 29) {
           Grid_Blks[iBlk].Extrude(Grid2D_Coflow_XYplane[iBlk-26][0],
                                   Input.NCells_Kdir,
			           STRETCHING_FCN_MIN_CLUSTERING,
                                   1.25,
                                   ZERO,
                                   0.25*Input.Length_Combustor_Tube);
	   Grid_Blks[iBlk].Set_BCs(BC_NONE,
		                   BC_NONE,
                                   BC_REFLECTION,
		                   BC_NONE,
                                   BC_NONE,
                                   BC_NONE);

        } else if (iBlk >= 30 && iBlk <= 33) {
           Grid_Blks[iBlk].Extrude(Grid2D_Coflow_XYplane[iBlk-30][0],
                                   Input.NCells_Kdir,
			           STRETCHING_FCN_LINEAR,
                                   ZERO,
                                   0.25*Input.Length_Combustor_Tube,
                                   Input.Length_Combustor_Tube);
	   Grid_Blks[iBlk].Set_BCs(BC_NONE,
		                   BC_NONE,
                                   BC_REFLECTION,
		                   BC_NONE,
                                   BC_FIXED_PRESSURE,
                                   BC_NONE);

        } else if (iBlk >= 34 && iBlk <= 37) {
           Grid_Blks[iBlk].Extrude(Grid2D_Coflow_XYplane[iBlk-34][0],
                                   Input.NCells_Kdir,
			           STRETCHING_FCN_MAX_CLUSTERING,
                                   1.25,
                                   -0.25*Input.Length_Coflow_Inlet_Pipe,
                                   ZERO);
	   Grid_Blks[iBlk].Set_BCs(BC_NONE,
		                   BC_NONE,
                                   BC_REFLECTION,
		                   BC_WALL_VISCOUS,
                                   BC_NONE,
                                   BC_NONE);

        } else if (iBlk >= 38 && iBlk <= 41) {
           Grid_Blks[iBlk].Extrude(Grid2D_Coflow_XYplane[iBlk-38][0],
                                   Input.NCells_Kdir,
			           STRETCHING_FCN_LINEAR,
                                   ONE,
                                   -Input.Length_Coflow_Inlet_Pipe,
                                   -0.25*Input.Length_Coflow_Inlet_Pipe);
	   Grid_Blks[iBlk].Set_BCs(BC_NONE,
		                   BC_NONE,
                                   BC_REFLECTION,
		                   BC_WALL_VISCOUS,
                                   BC_DIRICHLET,
                                   BC_NONE);

        } /* endif */

    } /* endfor */

    /* Deallocate 2D grids. */

    Grid2D_Fuel_Line_XYplane = Deallocate_Multi_Block_Grid(
                                   Grid2D_Fuel_Line_XYplane,
                                   numblk_idir_fuel,
		                   numblk_jdir_fuel);

    Grid2D_Bluff_Body_Inner_XYplane = Deallocate_Multi_Block_Grid(
                                   Grid2D_Bluff_Body_Inner_XYplane,
                                   numblk_idir_bluffbody,
		                   numblk_jdir_bluffbody);

    Grid2D_Bluff_Body_Outer_XYplane = Deallocate_Multi_Block_Grid(
                                   Grid2D_Bluff_Body_Outer_XYplane,
                                   numblk_idir_bluffbody,
		                   numblk_jdir_bluffbody);

    Grid2D_Coflow_XYplane = Deallocate_Multi_Block_Grid(
                                   Grid2D_Coflow_XYplane,
                                   numblk_idir_coflow,
		                   numblk_jdir_coflow);

}

/********************************************************
 * Routine: Create_Grid_ICEMCFD                         *
 *                                                      *
 * Read ICEMCFD Mesh                                    *
 *                                                      *
 ********************************************************/
// void Grid3D_Hexa_Multi_Block_List::Create_Grid_ICEMCFD(Grid3D_Input_Parameters &Input) {
    
//    if (Allocated) Deallocate();

//    Grid_Blks = Grid_ICEMCFD(Grid_Blks, 
//                             Input.ICEMCFD_FileNames, 
//                             NBlk_Idir, 
//                             NBlk_Jdir, 
//                             NBlk_Kdir);

//    int found = 0;
//    for (int kBlk = 0; kBlk <= NBlk_Kdir-1 && !found; ++kBlk) {
//       for (int jBlk = 0; jBlk <= NBlk_Jdir-1 && !found; ++jBlk) {
//          for (int iBlk = 0; iBlk <= NBlk_Idir-1 && !found; ++iBlk) {
//               if (Grid_Blks[iBlk][jBlk][kBlk].Allocated){
//                   Input.NCells_Idir = Grid_Blks[iBlk][jBlk][kBlk].NCi-2*Grid_Blks[iBlk][jBlk][kBlk].Nghost;
//                   Input.NCells_Jdir = Grid_Blks[iBlk][jBlk][kBlk].NCj-2*Grid_Blks[iBlk][jBlk][kBlk].Nghost;
//                   Input.NCells_Kdir = Grid_Blks[iBlk][jBlk][kBlk].NCk-2*Grid_Blks[iBlk][jBlk][kBlk].Nghost;
//                   found = 1;
//               }/* endif */   
//          } /* endfor */
//       } /* endfor */
//    } /* endfor */

//    Input.NBlk_Idir = NBlk_Idir; 
//    Input.NBlk_Jdir = NBlk_Jdir; 
//    Input.NBlk_Kdir = NBlk_Kdir;

//    assert(NBlk_Idir >= 1 && NBlk_Jdir >= 1 && NBlk_Kdir >= 1); 
//    Allocated = 1;

// }

/********************************************************
 * Routine: Create_Grid_Turbulent_Premixed_Flame        *
 *                                                      *
 * Generates a 3D Cartesian multiblock mesh for a       *
 *                                                      *
 * turbulent premixed flame box                         *
 ********************************************************/
void Grid3D_Hexa_Multi_Block_List::Create_Grid_Turbulent_Premixed_Flame(Grid3D_Input_Parameters &Input) {

    int count_blocks;
    int BC_east, BC_west, BC_north, BC_south, BC_top, BC_bottom;
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

    count_blocks = 0;

    for (int kBlk = 0; kBlk <= Input.NBlk_Kdir-1; ++kBlk) {
       for (int jBlk = 0; jBlk <= Input.NBlk_Jdir-1; ++jBlk) {
          for (int iBlk = 0; iBlk <= Input.NBlk_Idir-1; ++iBlk) {

             /* Extrude each of the grid blocks from the
                appropriate 2D grid in XY-plane. */

             Grid_Blks[count_blocks].Extrude(Grid2D_Box_XYplane[iBlk][jBlk],
                                             Input.NCells_Kdir,
			                     Input.Stretching_Type_Kdir,
				             Input.Stretching_Factor_Kdir,
                                             -HALF*Input.Box_Length+
                                             (double(kBlk)/double(Input.NBlk_Kdir))*Input.Box_Length,
                                             -HALF*Input.Box_Length+
                                             (double(kBlk+1)/double(Input.NBlk_Kdir))*Input.Box_Length);

             /* Assign boundary conditions. */

	     if (iBlk == Input.NBlk_Idir-1) {
	       BC_east = BC_OUTFLOW_SUBSONIC;
             } else {
                BC_east = BC_NONE;
             } /* endif */
             if (iBlk == 0) {
	       BC_west = BC_INFLOW_SUBSONIC;
             } else {
                BC_west = BC_NONE;
             } /* endif */

	     if (jBlk == Input.NBlk_Jdir-1) {
                BC_north = BC_PERIODIC;
             } else {
                BC_north = BC_NONE;
             } /* endif */
             if (jBlk == 0) {
                BC_south = BC_PERIODIC;
             } else {
                BC_south = BC_NONE;
             } /* endif */

	     if (kBlk == Input.NBlk_Kdir-1) {
                BC_top = BC_PERIODIC;
             } else {
                BC_top = BC_NONE;
             } /* endif */
             if (kBlk == 0) {
                BC_bottom = BC_PERIODIC;
             } else {
                BC_bottom = BC_NONE;
             } /* endif */

             Grid_Blks[count_blocks].Set_BCs(BC_east, 
                                             BC_west, 
                                             BC_north, 
                                             BC_south, 
                                             BC_top, 
                                             BC_bottom);

             /* Update block counter. */

             count_blocks ++;

	  } /* endfor */
       } /* endfor */
    } /* endfor */

    /* Deallocate 2D grid. */

    Grid2D_Box_XYplane = Deallocate_Multi_Block_Grid(Grid2D_Box_XYplane,
                                                     Input.NBlk_Idir, 
                                                     Input.NBlk_Jdir);

}

/********************************************************
 * Routine: Find_Neighbours                             *
 *                                                      *
 * Determine neighbouring grid blocks and the relative  *
 * orientation of the neighbouring blocks with respect  *
 * to each block in the multi-block list.               *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block_List::Find_Neighbours(void) {

   if (Allocated) {

   } /* endif */

}

/********************************************************
 * Routine: Allocate                                    *
 *                                                      *
 * Allocate memory for a 3D array of 3D hexahedral      *
 * multi-block grids.                                   *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block::Allocate(const int Ni, 
                                       const int Nj, 
                                       const int Nk) {

   if (Ni >= 1 && Nj >= 1 && Nk >= 1 && !Allocated) {
      NBlk_Idir = Ni; 
      NBlk_Jdir = Nj; 
      NBlk_Kdir = Nk;

      Grid_Blks = new Grid3D_Hexa_Block**[NBlk_Idir];
      for (int i = 0 ; i < NBlk_Idir ; ++i ) {
         Grid_Blks[i] = new Grid3D_Hexa_Block*[NBlk_Jdir];
         for (int j = 0 ; j <NBlk_Jdir ; ++j ) {
            Grid_Blks[i][j] = new Grid3D_Hexa_Block[NBlk_Kdir];
         }  /* endfor */
      }/* endfor */ 

      Allocated = 1;
   } /* endif */

}

/********************************************************
 * Routine: Deallocate                                  *
 *                                                      *
 * Deallocate memory for a 3D array of 3D hexahedral    *
 * multi-block grids.                                   *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block::Deallocate(void) {

   if (NBlk_Idir >= 1 && NBlk_Jdir >= 1 && NBlk_Kdir >= 1 && Allocated) {
      for (int i = 0 ; i <= NBlk_Idir-1 ; ++i ) {
         for ( int j=0 ; j <= NBlk_Jdir-1 ; ++j ) {
	    for ( int k = NBlk_Kdir-1 ; k >= 0; --k ) {
               if (Grid_Blks[i][j][k].Allocated) Grid_Blks[i][j][k].deallocate();
	    }/* end for */
	    delete []Grid_Blks[i][j];
	    Grid_Blks[i][j]=NULL;
          }  /* endfor */
       
          delete []Grid_Blks[i];
          Grid_Blks[i] = NULL;
       }  /* endfor */

       delete []Grid_Blks;
       Grid_Blks = NULL;

       NBlk_Idir = 0; 
       NBlk_Jdir = 0; 
       NBlk_Kdir = 0;

       Allocated = 0;
   } /* endif */

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
	     if (Grid2.Grid_Blks[i][j][k].Allocated) 
                Grid_Blks[i][j][k].Copy(Grid2.Grid_Blks[i][j][k]);
          } /* endfor */
       } /* endfor */
    } /* endfor */

  } /* endif */

}

/********************************************************
 * Routine: Broadcast                                   *
 *                                                      *
 * Broadcast multiblock hexahedral grid.                *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block::Broadcast(void) {

#ifdef _MPI_VERSION
  int ni, nj, nk, grid_allocated;

  /* Broadcast the number of grid blocks. */

  if (CFFC_Primary_MPI_Processor()) {
     ni = NBlk_Idir;
     nj = NBlk_Jdir;
     nk = NBlk_Kdir;
     grid_allocated = Allocated;
  } /* endif */

  MPI::COMM_WORLD.Bcast(&ni, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&nj, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&nk, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&grid_allocated, 1, MPI::INT, 0);

  /* On non-primary MPI processors, allocate (re-allocate) 
     memory for the grid blocks as necessary. */

  if (!CFFC_Primary_MPI_Processor()) {
     if (grid_allocated && 
         (NBlk_Idir != ni || 
          NBlk_Jdir != nj || 
          NBlk_Kdir != nk) ) { 
        if (Allocated) { 
          Deallocate();
        } /* endif */
        Allocate(ni, nj, nk);
     } /* endif */
  } /* endif */

  /* Broadcast each of the blocks in the multiblock mesh. */

  if (Allocated) {
    for (int  k = 0 ; k < NBlk_Kdir ; ++k ) {
       for (int  j = 0 ; j < NBlk_Jdir ; ++j ) {
          for (int  i = 0 ; i < NBlk_Idir ; ++i ) {
	     Grid_Blks[i][j][k].Broadcast();
          } /* endfor */
       } /* endfor */
    } /* endfor */
  } /* endif */
#endif

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
          if (Grid_Blks[i][j][k].Allocated) 
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
             if (Grid_Blks[i][j][k].Allocated) {
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
	     if (Grid_Blks[i][j][k].Allocated) {
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
	     if (Grid_Blks[i][j][k].Allocated) {
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
	      if (Grid_Blks[i][j][k].Allocated) {
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
void Grid3D_Hexa_Multi_Block::Create_Grid(Grid3D_Input_Parameters &Input) {

    switch(Input.i_Grid) {
      case GRID_CUBE :
        Create_Grid_Cube(Input);
        break;

      case GRID_CHANNEL_XDIR :
      case GRID_CHANNEL_YDIR:
      case GRID_CHANNEL_ZDIR:
        Create_Grid_Channel(Input);
        break;

      case GRID_COUETTE_XDIR :
      case GRID_COUETTE_YDIR:
      case GRID_COUETTE_ZDIR:
        Create_Grid_Couette(Input);
        break;

      case GRID_PIPE :
        Create_Grid_Pipe(Input);
        break;

      case GRID_BLUFF_BODY_BURNER :
        Create_Grid_Bluff_Body_Burner(Input);
        break;

      case GRID_ICEMCFD :
        Create_Grid_ICEMCFD(Input);
        break;

      case GRID_TURBULENT_PREMIXED_FLAME :
	Create_Grid_Turbulent_Premixed_Flame(Input);
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
void Grid3D_Hexa_Multi_Block::Create_Grid_Cube(Grid3D_Input_Parameters &Input) {

    int BC_top, BC_bottom;
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

             /* Extrude each of the grid blocks from the
                appropriate 2D grid in XY-plane. */

             Grid_Blks[iBlk][jBlk][kBlk].Extrude(Grid2D_Box_XYplane[iBlk][jBlk],
                                                 Input.NCells_Kdir,
			                         Input.Stretching_Type_Kdir,
						 Input.Stretching_Factor_Kdir,
                                                 -HALF*Input.Box_Length+
                                                 (double(kBlk)/double(Input.NBlk_Kdir))*Input.Box_Length,
                                                 -HALF*Input.Box_Length+
                                                 (double(kBlk+1)/double(Input.NBlk_Kdir))*Input.Box_Length);

             /* Assign top and bottom boundary conditions. */

	     if (kBlk == Input.NBlk_Kdir-1) {
                BC_top = BC_REFLECTION;
             } else {
                BC_top = BC_NONE;
             } /* endif */
             if (kBlk == 0) {
                BC_bottom = BC_REFLECTION;
             } else {
                BC_bottom = BC_NONE;
             } /* endif */

             Grid_Blks[iBlk][jBlk][kBlk].Set_BCs_Zdir(BC_top, BC_bottom);

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
void Grid3D_Hexa_Multi_Block::Create_Grid_Channel(Grid3D_Input_Parameters &Input) {

    int BC_east, BC_west, BC_north, BC_south, BC_top, BC_bottom;
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

             /* Extrude each of the grid blocks from the
                appropriate 2D grid in XY-plane. */

             Grid_Blks[iBlk][jBlk][kBlk].Extrude(Grid2D_Box_XYplane[iBlk][jBlk],
                                                 Input.NCells_Kdir,
			                         Input.Stretching_Type_Kdir,
						 Input.Stretching_Factor_Kdir,
                                                 (double(kBlk)/double(Input.NBlk_Kdir))*Input.Box_Length,
                                                 (double(kBlk+1)/double(Input.NBlk_Kdir))*Input.Box_Length);

             /* Assign top and bottom boundary conditions. */

             if (Input.i_Grid == GRID_CHANNEL_ZDIR) {

   	        if (iBlk == Input.NBlk_Idir-1) {
                   BC_east = BC_CONSTANT_EXTRAPOLATION;
                } else {
                   BC_east = BC_NONE;
                } /* endif */
                if (iBlk == 0) {
                   BC_west = BC_CONSTANT_EXTRAPOLATION;
                } else {
                   BC_west = BC_NONE;
                } /* endif */

	        if (jBlk == Input.NBlk_Jdir-1) {
                   BC_north = BC_WALL_VISCOUS;
                } else {
                   BC_north = BC_NONE;
                } /* endif */
                if (jBlk == 0) {
                   BC_south = BC_WALL_VISCOUS;
                } else {
                   BC_south = BC_NONE;
                } /* endif */

	        if (kBlk == Input.NBlk_Kdir-1) {
                   BC_top = BC_CHANNEL_OUTFLOW;
                } else {
                   BC_top = BC_NONE;
                } /* endif */
                if (kBlk == 0) {
                   BC_bottom = BC_CHANNEL_INFLOW;
                } else {
                   BC_bottom = BC_NONE;
                } /* endif */

             } else if (Input.i_Grid == GRID_CHANNEL_XDIR) {

    	        if (iBlk == Input.NBlk_Idir-1) {
                   BC_east = BC_CHANNEL_OUTFLOW;
                } else {
                   BC_east = BC_NONE;
                } /* endif */
                if (iBlk == 0) {
                   BC_west = BC_CHANNEL_INFLOW;
                } else {
                   BC_west = BC_NONE;
                } /* endif */

	        if (jBlk == Input.NBlk_Jdir-1) {
                   BC_north = BC_WALL_VISCOUS;
                } else {
                   BC_north = BC_NONE;
                } /* endif */
                if (jBlk == 0) {
                   BC_south = BC_WALL_VISCOUS;
                } else {
                   BC_south = BC_NONE;
                } /* endif */

	        if (kBlk == Input.NBlk_Kdir-1) {
                   BC_top = BC_CONSTANT_EXTRAPOLATION;
                } else {
                   BC_top = BC_NONE;
                } /* endif */
                if (kBlk == 0) {
                   BC_bottom = BC_CONSTANT_EXTRAPOLATION;
                } else {
                   BC_bottom = BC_NONE;
                } /* endif */

             } else if (Input.i_Grid == GRID_CHANNEL_YDIR) {

   	        if (iBlk == Input.NBlk_Idir-1) {
                   BC_east = BC_CONSTANT_EXTRAPOLATION;
                } else {
                   BC_east = BC_NONE;
                } /* endif */
                if (iBlk == 0) {
                   BC_west = BC_CONSTANT_EXTRAPOLATION;
                } else {
                   BC_west = BC_NONE;
                } /* endif */

	        if (jBlk == Input.NBlk_Jdir-1) {
                   BC_north = BC_CHANNEL_OUTFLOW;
                } else {
                   BC_north = BC_NONE;
                } /* endif */
                if (jBlk == 0) {
                   BC_south = BC_CHANNEL_INFLOW;
                } else {
                   BC_south = BC_NONE;
                } /* endif */

	        if (kBlk == Input.NBlk_Kdir-1) {
                   BC_top = BC_WALL_VISCOUS;
                } else {
                   BC_top = BC_NONE;
                } /* endif */
                if (kBlk == 0) {
                   BC_bottom = BC_WALL_VISCOUS;
                } else {
                   BC_bottom = BC_NONE;
                } /* endif */
	     } /* endif */

             Grid_Blks[iBlk][jBlk][kBlk].Set_BCs(BC_east, 
                                                 BC_west, 
                                                 BC_north, 
                                                 BC_south, 
                                                 BC_top, 
                                                 BC_bottom);
 
	  } /* endfor */
       } /* endfor */
    } /* endfor */

    /* Deallocate 2D grid. */

    Grid2D_Box_XYplane = Deallocate_Multi_Block_Grid(Grid2D_Box_XYplane,
                                                     Input.NBlk_Idir, 
                                                     Input.NBlk_Jdir);

}

/********************************************************
 * Routine: Create_Grid_Couette                         *
 *                                                      *
 * Generates a 3D Cartesian multiblock mesh for a       *
 * Couette flow.                                        *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block::Create_Grid_Couette(Grid3D_Input_Parameters &Input) {
      
    int BC_east, BC_west, BC_north, BC_south, BC_top, BC_bottom;
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

             /* Extrude each of the grid blocks from the
                appropriate 2D grid in XY-plane. */

             Grid_Blks[iBlk][jBlk][kBlk].Extrude(Grid2D_Box_XYplane[iBlk][jBlk],
                                                 Input.NCells_Kdir,
			                         Input.Stretching_Type_Kdir,
						 Input.Stretching_Factor_Kdir,
                                                 (double(kBlk)/double(Input.NBlk_Kdir))*Input.Box_Length,
                                                 (double(kBlk+1)/double(Input.NBlk_Kdir))*Input.Box_Length);

             /* Assign top and bottom boundary conditions. */

             if (Input.i_Grid == GRID_COUETTE_ZDIR) {

   	        if (iBlk == Input.NBlk_Idir-1) {
                   BC_east = BC_CONSTANT_EXTRAPOLATION;
                } else {
                   BC_east = BC_NONE;
                } /* endif */
                if (iBlk == 0) {
                   BC_west = BC_CONSTANT_EXTRAPOLATION;
                } else {
                   BC_west = BC_NONE;
                } /* endif */

	        if (jBlk == Input.NBlk_Jdir-1) {
		   BC_north = BC_MOVING_WALL;
                } else {
                   BC_north = BC_NONE;
                } /* endif */
                if (jBlk == 0) {
                   BC_south = BC_NO_SLIP;
                } else {
                   BC_south = BC_NONE;
                } /* endif */

	        if (kBlk == Input.NBlk_Kdir-1) {
                   BC_top = BC_FIXED_PRESSURE;
                } else {
                   BC_top = BC_NONE;
                } /* endif */
                if (kBlk == 0) {
                   BC_bottom = BC_FIXED_PRESSURE;
                } else {
                   BC_bottom = BC_NONE;
                } /* endif */

             } else if (Input.i_Grid == GRID_COUETTE_XDIR) {

   	        if (iBlk == Input.NBlk_Idir-1) {
                   BC_east = BC_FIXED_PRESSURE;
                } else {
                   BC_east = BC_NONE;
                } /* endif */
                if (iBlk == 0) {
                   BC_west = BC_FIXED_PRESSURE;
                } else {
                   BC_west = BC_NONE;
                } /* endif */

	        if (jBlk == Input.NBlk_Jdir-1) {
                   BC_north = BC_MOVING_WALL;
                } else {
                   BC_north = BC_NONE;
                } /* endif */
                if (jBlk == 0) {
                   BC_south = BC_NO_SLIP;
                } else {
                   BC_south = BC_NONE;
                } /* endif */

	        if (kBlk == Input.NBlk_Kdir-1) {
                   BC_top = BC_CONSTANT_EXTRAPOLATION;
                } else {
                   BC_top = BC_NONE;
                } /* endif */
                if (kBlk == 0) {
                   BC_bottom = BC_CONSTANT_EXTRAPOLATION;
                } else {
                   BC_bottom = BC_NONE;
                } /* endif */

             } else if (Input.i_Grid == GRID_COUETTE_YDIR) {

   	        if (iBlk == Input.NBlk_Idir-1) {
                   BC_east = BC_CONSTANT_EXTRAPOLATION;
                } else {
                   BC_east = BC_NONE;
                } /* endif */
                if (iBlk == 0) {
                   BC_west = BC_CONSTANT_EXTRAPOLATION;
                } else {
                   BC_west = BC_NONE;
                } /* endif */

	        if (jBlk == Input.NBlk_Jdir-1) {
                   BC_north = BC_FIXED_PRESSURE;
                } else {
                   BC_north = BC_NONE;
                } /* endif */
                if (jBlk == 0) {
                   BC_south = BC_FIXED_PRESSURE;
                } else {
                   BC_south = BC_NONE;
                } /* endif */

	        if (kBlk == Input.NBlk_Kdir-1) {
                   BC_top = BC_MOVING_WALL;
                } else {
                   BC_top = BC_NONE;
                } /* endif */
                if (kBlk == 0) {
                   BC_bottom = BC_NO_SLIP;
                } else {
                   BC_bottom = BC_NONE;
                } /* endif */
	     } /* endif */

             Grid_Blks[iBlk][jBlk][kBlk].Set_BCs(BC_east, 
                                                 BC_west, 
                                                 BC_north, 
                                                 BC_south, 
                                                 BC_top, 
                                                 BC_bottom);
 
	  } /* endfor */
       } /* endfor */
    } /* endfor */

    /* Deallocate 2D grid. */

    Grid2D_Box_XYplane = Deallocate_Multi_Block_Grid(Grid2D_Box_XYplane,
                                                     Input.NBlk_Idir, 
                                                     Input.NBlk_Jdir);

}

/********************************************************
 * Routine: Create_Grid_Pipe                            *
 *                                                      *
 * Generates a 3D multiblock mesh for pipe flows.       *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block::Create_Grid_Pipe(Grid3D_Input_Parameters &Input) {

    int numblk_idir_pipe, numblk_jdir_pipe;
    Grid2D_Quad_Block **Grid2D_Pipe_XYplane;

    /* Allocate required memory. */

    Input.NBlk_Idir = 5;
    Input.NBlk_Jdir = 1;
    Input.NBlk_Kdir = 1;

    Allocate(Input.NBlk_Idir, Input.NBlk_Jdir, Input.NBlk_Kdir);

    /* Creat 2D cross-section grids from which the 3D grid
       will be extruded. */

    Grid2D_Pipe_XYplane = Grid_Tube_2D(
                                   Grid2D_Pipe_XYplane,
                                   numblk_idir_pipe,
		                   numblk_jdir_pipe,
                                   Input.Pipe_Radius,
                                   Input.NCells_Idir,
                                   Input.NCells_Jdir,
				   Input.Nghost,
                                   STRETCHING_FCN_MAX_CLUSTERING,
                                   1.25);

    /* Create the mesh for each block representing
       the complete grid. */

    for (int iBlk = 0; iBlk <= Input.NBlk_Idir-1; ++iBlk) {

        /* Extrude each of the grid blocks from the
           appropriate 2D grid in XY-plane. */

        Grid_Blks[iBlk][0][0].Extrude(Grid2D_Pipe_XYplane[iBlk][0],
                                      Input.NCells_Kdir,
       	  	                      Input.Stretching_Type_Kdir,
				      Input.Stretching_Factor_Kdir,
                                      ZERO,
                                      Input.Pipe_Length);
        if (iBlk == 0) {
           Grid_Blks[iBlk][0][0].Set_BCs(BC_NONE,
         	                         BC_NONE,
                                         BC_NONE,
       	                                 BC_NONE,
                                         BC_NONE,
                                         BC_DIRICHLET);
        } else {
           Grid_Blks[iBlk][0][0].Set_BCs(BC_NONE,
       	                                 BC_NONE,
					 BC_WALL_VISCOUS,
       	                                 BC_NONE,
                                         BC_NONE,
                                         BC_DIRICHLET);
        } /* endif */

    } /* endfor */

    /* Deallocate 2D grid. */

    Grid2D_Pipe_XYplane = Deallocate_Multi_Block_Grid(
                                   Grid2D_Pipe_XYplane,
                                   numblk_idir_pipe,
		                   numblk_jdir_pipe);

}

/********************************************************
 * Routine: Create_Grid_Bluff_Body_Burner               *
 *                                                      *
 * Generates a 3D multiblock mesh for TNF bluff body    *
 * burner.                                              *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block::Create_Grid_Bluff_Body_Burner(Grid3D_Input_Parameters &Input) {

    int numblk_idir_fuel, numblk_jdir_fuel,
        numblk_idir_bluffbody, numblk_jdir_bluffbody,
        numblk_idir_coflow, numblk_jdir_coflow;
    Grid2D_Quad_Block **Grid2D_Fuel_Line_XYplane,
                      **Grid2D_Bluff_Body_Inner_XYplane,
                      **Grid2D_Bluff_Body_Outer_XYplane,
                      **Grid2D_Coflow_XYplane;

    /* Allocate required memory. */

    Input.NBlk_Idir = 42;
    Input.NBlk_Jdir = 1;
    Input.NBlk_Kdir = 1;

    Allocate(Input.NBlk_Idir, Input.NBlk_Jdir, Input.NBlk_Kdir);

    /* Creat 2D cross-section grids from which the 3D grid
       will be extruded. */

    Grid2D_Fuel_Line_XYplane = Grid_Tube_2D(
                                   Grid2D_Fuel_Line_XYplane,
                                   numblk_idir_fuel,
		                   numblk_jdir_fuel,
                                   Input.Radius_Fuel_Line,
                                   Input.NCells_Idir,
                                   Input.NCells_Jdir,
				   Input.Nghost,
                                   STRETCHING_FCN_MAX_CLUSTERING,
                                   1.25);

    Grid2D_Bluff_Body_Inner_XYplane = Grid_Annulus_2D(
                                   Grid2D_Bluff_Body_Inner_XYplane,
                                   numblk_idir_bluffbody,
		                   numblk_jdir_bluffbody,
                                   Input.Radius_Fuel_Line,
 			           HALF*Input.Radius_Bluff_Body,
                                   Input.NCells_Idir,
                                   Input.NCells_Jdir,
				   Input.Nghost,
                                   STRETCHING_FCN_MIN_CLUSTERING,
                                   1.10);

    Grid2D_Bluff_Body_Outer_XYplane = Grid_Annulus_2D(
                                   Grid2D_Bluff_Body_Outer_XYplane,
                                   numblk_idir_bluffbody,
		                   numblk_jdir_bluffbody,
                                   HALF*Input.Radius_Bluff_Body,
 			           Input.Radius_Bluff_Body,
                                   Input.NCells_Idir,
                                   Input.NCells_Jdir,
				   Input.Nghost,
                                   STRETCHING_FCN_MAX_CLUSTERING,
                                   1.10);

    Grid2D_Coflow_XYplane = Grid_Annulus_2D(
                                   Grid2D_Coflow_XYplane,
                                   numblk_idir_coflow,
		                   numblk_jdir_coflow,
 			           Input.Radius_Bluff_Body,
                                   Input.Radius_Coflow_Inlet_Pipe,
                                   Input.NCells_Idir,
                                   Input.NCells_Jdir,
				   Input.Nghost,
                                   STRETCHING_FCN_MIN_CLUSTERING,
                                   1.10);

    /* Create the mesh for each block representing
       the complete grid. */

    for (int iBlk = 0; iBlk <= Input.NBlk_Idir-1; ++iBlk) {

        /* Extrude each of the grid blocks from the
           appropriate 2D grid in XY-plane. */

        if (iBlk <= 4) {
           Grid_Blks[iBlk][0][0].Extrude(
                              Grid2D_Fuel_Line_XYplane[iBlk][0],
                              Input.NCells_Kdir,
			      STRETCHING_FCN_MIN_CLUSTERING,
                              1.25,
                              ZERO,
                              0.25*Input.Length_Combustor_Tube);
           if (iBlk == 0) {
	      Grid_Blks[iBlk][0][0].Set_BCs(BC_NONE,
		                            BC_NONE,
                                            BC_NONE,
		                            BC_NONE,
                                            BC_NONE,
                                            BC_DIRICHLET);
	   } else {
              Grid_Blks[iBlk][0][0].Set_BCs(BC_NONE,
		                            BC_NONE,
                                            BC_NONE,
		                            BC_NONE,
                                            BC_NONE,
                                            BC_DIRICHLET);
           }  /* endif */

        } else if (iBlk >= 5 && iBlk <= 9) {
           Grid_Blks[iBlk][0][0].Extrude(
                              Grid2D_Fuel_Line_XYplane[iBlk-5][0],
                              Input.NCells_Kdir,
			      STRETCHING_FCN_LINEAR,
                              ZERO,
                              0.25*Input.Length_Combustor_Tube,
                              Input.Length_Combustor_Tube);
           if (iBlk-5 == 0) {
	      Grid_Blks[iBlk][0][0].Set_BCs(BC_NONE,
		                            BC_NONE,
                                            BC_NONE,
		                            BC_NONE,
                                            BC_FIXED_PRESSURE,
                                            BC_NONE);
	   } else {
              Grid_Blks[iBlk][0][0].Set_BCs(BC_NONE,
		                            BC_NONE,
                                            BC_NONE,
		                            BC_NONE,
                                            BC_FIXED_PRESSURE,
                                            BC_NONE);
           } /* endif */

        } else if (iBlk >= 10 && iBlk <= 13) {
           Grid_Blks[iBlk][0][0].Extrude(
                              Grid2D_Bluff_Body_Inner_XYplane[iBlk-10][0],
                              Input.NCells_Kdir,
			      STRETCHING_FCN_MIN_CLUSTERING,
                              1.25,
                              ZERO,
                              0.25*Input.Length_Combustor_Tube);
	   Grid_Blks[iBlk][0][0].Set_BCs(BC_NONE,
		                         BC_NONE,
                                         BC_NONE,
		                         BC_NONE,
                                         BC_NONE,
                                         BC_WALL_VISCOUS);

        } else if (iBlk >= 14 && iBlk <= 17) {
           Grid_Blks[iBlk][0][0].Extrude(
                              Grid2D_Bluff_Body_Inner_XYplane[iBlk-14][0],
                              Input.NCells_Kdir,
			      STRETCHING_FCN_LINEAR,
                              ZERO,
                              0.25*Input.Length_Combustor_Tube,
                              Input.Length_Combustor_Tube);
	   Grid_Blks[iBlk][0][0].Set_BCs(BC_NONE,
		                         BC_NONE,
                                         BC_NONE,
		                         BC_NONE,
                                         BC_FIXED_PRESSURE,
                                         BC_NONE);

        } else if (iBlk >= 18 && iBlk <= 21) {
           Grid_Blks[iBlk][0][0].Extrude(
                              Grid2D_Bluff_Body_Outer_XYplane[iBlk-18][0],
                              Input.NCells_Kdir,
			      STRETCHING_FCN_MIN_CLUSTERING,
                              1.25,
                              ZERO,
                              0.25*Input.Length_Combustor_Tube);
	   Grid_Blks[iBlk][0][0].Set_BCs(BC_NONE,
		                         BC_NONE,
                                         BC_NONE,
		                         BC_NONE,
                                         BC_NONE,
                                         BC_WALL_VISCOUS);

        } else if (iBlk >= 22 && iBlk <= 25) {
           Grid_Blks[iBlk][0][0].Extrude(
                              Grid2D_Bluff_Body_Outer_XYplane[iBlk-22][0],
                              Input.NCells_Kdir,
			      STRETCHING_FCN_LINEAR,
                              ONE,
                              0.25*Input.Length_Combustor_Tube,
                              Input.Length_Combustor_Tube);
	   Grid_Blks[iBlk][0][0].Set_BCs(BC_NONE,
		                         BC_NONE,
                                         BC_NONE,
		                         BC_NONE,
                                         BC_FIXED_PRESSURE,
                                         BC_NONE);

        } else if (iBlk >= 26 && iBlk <= 29) {
           Grid_Blks[iBlk][0][0].Extrude(
                              Grid2D_Coflow_XYplane[iBlk-26][0],
                              Input.NCells_Kdir,
			      STRETCHING_FCN_MIN_CLUSTERING,
                              1.25,
                              ZERO,
                              0.25*Input.Length_Combustor_Tube);
	   Grid_Blks[iBlk][0][0].Set_BCs(BC_NONE,
		                         BC_NONE,
                                         BC_REFLECTION,
		                         BC_NONE,
                                         BC_NONE,
                                         BC_NONE);

        } else if (iBlk >= 30 && iBlk <= 33) {
           Grid_Blks[iBlk][0][0].Extrude(
                              Grid2D_Coflow_XYplane[iBlk-30][0],
                              Input.NCells_Kdir,
			      STRETCHING_FCN_LINEAR,
                              ZERO,
                              0.25*Input.Length_Combustor_Tube,
                              Input.Length_Combustor_Tube);
	   Grid_Blks[iBlk][0][0].Set_BCs(BC_NONE,
		                         BC_NONE,
                                         BC_REFLECTION,
		                         BC_NONE,
                                         BC_FIXED_PRESSURE,
                                         BC_NONE);

        } else if (iBlk >= 34 && iBlk <= 37) {
           Grid_Blks[iBlk][0][0].Extrude(
                              Grid2D_Coflow_XYplane[iBlk-34][0],
                              Input.NCells_Kdir,
			      STRETCHING_FCN_MAX_CLUSTERING,
                              1.25,
                              -0.25*Input.Length_Coflow_Inlet_Pipe,
                              ZERO);
	   Grid_Blks[iBlk][0][0].Set_BCs(BC_NONE,
		                         BC_NONE,
                                         BC_REFLECTION,
		                         BC_WALL_VISCOUS,
                                         BC_NONE,
                                         BC_NONE);

        } else if (iBlk >= 38 && iBlk <= 41) {
           Grid_Blks[iBlk][0][0].Extrude(
                              Grid2D_Coflow_XYplane[iBlk-38][0],
                              Input.NCells_Kdir,
			      STRETCHING_FCN_LINEAR,
                              ONE,
                              -Input.Length_Coflow_Inlet_Pipe,
                              -0.25*Input.Length_Coflow_Inlet_Pipe);
	   Grid_Blks[iBlk][0][0].Set_BCs(BC_NONE,
		                         BC_NONE,
                                         BC_REFLECTION,
		                         BC_WALL_VISCOUS,
                                         BC_DIRICHLET,
                                         BC_NONE);

        } /* endif */

    } /* endfor */

    /* Deallocate 2D grids. */

    Grid2D_Fuel_Line_XYplane = Deallocate_Multi_Block_Grid(
                                   Grid2D_Fuel_Line_XYplane,
                                   numblk_idir_fuel,
		                   numblk_jdir_fuel);

    Grid2D_Bluff_Body_Inner_XYplane = Deallocate_Multi_Block_Grid(
                                   Grid2D_Bluff_Body_Inner_XYplane,
                                   numblk_idir_bluffbody,
		                   numblk_jdir_bluffbody);

    Grid2D_Bluff_Body_Outer_XYplane = Deallocate_Multi_Block_Grid(
                                   Grid2D_Bluff_Body_Outer_XYplane,
                                   numblk_idir_bluffbody,
		                   numblk_jdir_bluffbody);

    Grid2D_Coflow_XYplane = Deallocate_Multi_Block_Grid(
                                   Grid2D_Coflow_XYplane,
                                   numblk_idir_coflow,
		                   numblk_jdir_coflow);

}

/********************************************************
 * Routine: Create_Grid_ICEMCFD                         *
 *                                                      *
 * Read ICEMCFD Mesh                                    *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block::Create_Grid_ICEMCFD(Grid3D_Input_Parameters &Input) {
    
   if (Allocated) Deallocate();

   Grid_Blks = Grid_ICEMCFD(Grid_Blks, 
                            Input.ICEMCFD_FileNames, 
                            NBlk_Idir, 
                            NBlk_Jdir, 
                            NBlk_Kdir);

   int found = 0;
   for (int kBlk = 0; kBlk <= NBlk_Kdir-1 && !found; ++kBlk) {
      for (int jBlk = 0; jBlk <= NBlk_Jdir-1 && !found; ++jBlk) {
         for (int iBlk = 0; iBlk <= NBlk_Idir-1 && !found; ++iBlk) {
              if (Grid_Blks[iBlk][jBlk][kBlk].Allocated){
                  Input.NCells_Idir = Grid_Blks[iBlk][jBlk][kBlk].NCi-2*Grid_Blks[iBlk][jBlk][kBlk].Nghost;
                  Input.NCells_Jdir = Grid_Blks[iBlk][jBlk][kBlk].NCj-2*Grid_Blks[iBlk][jBlk][kBlk].Nghost;
                  Input.NCells_Kdir = Grid_Blks[iBlk][jBlk][kBlk].NCk-2*Grid_Blks[iBlk][jBlk][kBlk].Nghost;
                  found = 1;
              }/* endif */   
         } /* endfor */
      } /* endfor */
   } /* endfor */

   Input.NBlk_Idir = NBlk_Idir; 
   Input.NBlk_Jdir = NBlk_Jdir; 
   Input.NBlk_Kdir = NBlk_Kdir;

   assert(NBlk_Idir >= 1 && NBlk_Jdir >= 1 && NBlk_Kdir >= 1); 
   Allocated = 1;

}         

/********************************************************
 * Routine: Create_Grid_Turbulent_Premixed_Flame        *
 *                                                      *
 * Generates a 3D Cartesian multiblock mesh for a       *
 *                                                      *
 * turbulent premixed flame box                         *
 ********************************************************/
void Grid3D_Hexa_Multi_Block::Create_Grid_Turbulent_Premixed_Flame(Grid3D_Input_Parameters &Input) {

    int count_blocks;
    int BC_east, BC_west, BC_north, BC_south, BC_top, BC_bottom;
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

    //    count_blocks = 0;

    for (int kBlk = 0; kBlk <= Input.NBlk_Kdir-1; ++kBlk) {
       for (int jBlk = 0; jBlk <= Input.NBlk_Jdir-1; ++jBlk) {
          for (int iBlk = 0; iBlk <= Input.NBlk_Idir-1; ++iBlk) {

             /* Extrude each of the grid blocks from the
                appropriate 2D grid in XY-plane. */

             Grid_Blks[iBlk][jBlk][kBlk].Extrude(Grid2D_Box_XYplane[iBlk][jBlk],
                                             Input.NCells_Kdir,
			                     Input.Stretching_Type_Kdir,
				             Input.Stretching_Factor_Kdir,
                                             -HALF*Input.Box_Length+
                                             (double(kBlk)/double(Input.NBlk_Kdir))*Input.Box_Length,
                                             -HALF*Input.Box_Length+
                                             (double(kBlk+1)/double(Input.NBlk_Kdir))*Input.Box_Length);

             /* Assign boundary conditions. */

	     if (iBlk == Input.NBlk_Idir-1) {
                BC_east = BC_INFLOW_SUBSONIC;
             } else {
                BC_east = BC_NONE;
             } /* endif */
             if (iBlk == 0) {
                BC_west = BC_OUTFLOW_SUBSONIC;
             } else {
                BC_west = BC_NONE;
             } /* endif */

	     if (jBlk == Input.NBlk_Jdir-1) {
                BC_north = BC_PERIODIC;
             } else {
                BC_north = BC_NONE;
             } /* endif */
             if (jBlk == 0) {
                BC_south = BC_PERIODIC;
             } else {
                BC_south = BC_NONE;
             } /* endif */

	     if (kBlk == Input.NBlk_Kdir-1) {
                BC_top = BC_PERIODIC;
             } else {
                BC_top = BC_NONE;
             } /* endif */
             if (kBlk == 0) {
                BC_bottom = BC_PERIODIC;
             } else {
                BC_bottom = BC_NONE;
             } /* endif */

             Grid_Blks[iBlk][jBlk][kBlk].Set_BCs(BC_east, 
                                             BC_west, 
                                             BC_north, 
                                             BC_south, 
                                             BC_top, 
                                             BC_bottom);

             /* Update block counter. */

	     //             count_blocks ++;

	  } /* endfor */
       } /* endfor */
    } /* endfor */

    /* Deallocate 2D grid. */

    Grid2D_Box_XYplane = Deallocate_Multi_Block_Grid(Grid2D_Box_XYplane,
                                                     Input.NBlk_Idir, 
                                                     Input.NBlk_Jdir);

}

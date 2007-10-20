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
      Connectivity = new Grid3D_Hexa_Multi_Block_Connectivity[NBlk];
      Allocated = 1;
         
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
      Connectivity = new Grid3D_Hexa_Multi_Block_Connectivity[NBlk];

      Allocated = 1;
     
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

       delete []Connectivity;
       Connectivity = NULL;

  
       NBlk_Idir = 0; 
       NBlk_Jdir = 0; 
       NBlk_Kdir = 0;
       NBlk = 0;

       Allocated = 0;
   } /* endif */

}



/********************************************************
 * Routine: Deallocate                                  *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block_Connectivity::Deallocate(void) {

   
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


       delete []neighNW_info; 
       delete []neighNE_info; 
       delete []neighSE_info;
       delete []neighSW_info; 
       delete []neighTN_info; 
       delete []neighTS_info; 
       delete []neighTE_info; 
       delete []neighTW_info; 
       delete []neighTNW_info;
       delete []neighTSW_info; 
       delete []neighTNE_info; 
       delete []neighTSE_info; 
       delete []neighBN_info; 
       delete []neighBS_info; 
       delete []neighBE_info; 
       delete []neighBW_info; 
       delete []neighBNW_info;
       delete []neighBSW_info;
       delete []neighBNE_info;
       delete []neighBSE_info;
    
       
}
/********************************************************
 * Routine: Broadcast                                   *
 *                                                      *
 * Broadcast block connectivity info.                   *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block_Connectivity::Broadcast(void) {
#ifdef _MPI_VERSION

  MPI::COMM_WORLD.Bcast(&num_neighT, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&num_neighB, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&num_neighS, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&num_neighW, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&num_neighE, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&num_neighN, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&num_neighTN, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&num_neighTS, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&num_neighTW, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&num_neighTE, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&num_neighBN, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&num_neighBS, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&num_neighBE, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&num_neighBW, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&num_neighSW, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&num_neighSE, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&num_neighNW, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&num_neighNE, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&num_neighTNW, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&num_neighTNE, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&num_neighTSW, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&num_neighTSE, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&num_neighBNW, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&num_neighBNE, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&num_neighBSW, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&num_neighBSE, 1, MPI::INT, 0);
  
  MPI::COMM_WORLD.Bcast(&neighT, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighB, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighS, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighW, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighE, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighN, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighTN[0], 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighTS[0], 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighTW[0], 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighTE[0], 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighBN[0], 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighBS[0], 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighBE[0], 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighBW[0], 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighSW[0], 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighSE[0], 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighNW[0], 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighNE[0], 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighTNW[0], 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighTNE[0], 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighTSW[0], 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighTSE[0], 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighBNW[0], 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighBNE[0], 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighBSW[0], 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighBSE[0], 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighTN[1], 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighTS[1], 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighTW[1], 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighTE[1], 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighBN[1], 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighBS[1], 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighBE[1], 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighBW[1], 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighSW[1], 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighSE[1], 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighNW[1], 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighNE[1], 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighTNW[1], 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighTNE[1], 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighTSW[1], 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighTSE[1], 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighBNW[1], 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighBNE[1], 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighBSW[1], 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&neighBSE[1], 1, MPI::INT, 0);

  neighT_info.broadcast();
  neighB_info.broadcast();
  neighN_info.broadcast();
  neighS_info.broadcast();
  neighE_info.broadcast();
  neighW_info.broadcast();
  
  for(int i_nb = 0; i_nb < BLOCK_ORIENTATION_MAX_NEIGHBOUR; ++i_nb){
     
     neighNW_info[i_nb].broadcast();
     neighNE_info[i_nb].broadcast();
     neighSW_info[i_nb].broadcast();
     neighSE_info[i_nb].broadcast();
     neighTN_info[i_nb].broadcast();
     neighTE_info[i_nb].broadcast();
     neighTW_info[i_nb].broadcast();
     neighTS_info[i_nb].broadcast();
     neighBN_info[i_nb].broadcast();
     neighBE_info[i_nb].broadcast();
     neighBW_info[i_nb].broadcast();
     neighBS_info[i_nb].broadcast();
     neighTNW_info[i_nb].broadcast();
     neighTNE_info[i_nb].broadcast();
     neighTSW_info[i_nb].broadcast();
     neighTSE_info[i_nb].broadcast();
     neighBNW_info[i_nb].broadcast();
     neighBNE_info[i_nb].broadcast();
     neighBSW_info[i_nb].broadcast();
     neighBSE_info[i_nb].broadcast();
   
  }
  
#endif  


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
  /* Broadcast block connectivity info in the multiblock mesh. */
//  CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
  if (Allocated ) {
     for (int  i = 0 ; i < NBlk ; ++i ) {
        Grid_Blks[i].Broadcast();
        Connectivity[i].Broadcast();

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

      default:
        Create_Grid_Cube(Input);
        break;
    } /* endswitch */

    // call the function Find_Neighbours to obtain the neighbour block information
    // and assign values to data members in  Grid3D_Hexa_Multi_Block_Connectivity. 
    
    Find_Neighbours(Input);

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
 * Routine: Find_Neighbours                             *
 *                                                      *
 * Determine neighbouring grid blocks and the relative  *
 * orientation of the neighbouring blocks with respect  *
 * to each block in the multi-block list.               *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block_List::Find_Neighbours(Grid3D_Input_Parameters &Input) {
   
   if (Allocated) {
      // add all the blocks in the database in order to compute the
      // connectivity between blocks
      // constructor for blockConn
      
      BlkC::BlockConnectivity blkConn(NBlk, BlkC::CellCenter,Input.Nghost, 6);
      //blockConn.dbg_save_input("blkc_bluffbody.dat");
      
      BlkC::VertexPack vp;
      int id,jd,kd;
      int id_n, jd_n, kd_n;
      Vector3D dX_i, dX_j, dX_k;

      for (int iblk = 0 ; iblk < NBlk ; ++iblk ){
         // a vertex pack (8 vertices and 24 vecotrs) for each block
         for ( int i = 0; i != 2; ++i ) 
            for ( int j = 0; j != 2; ++j ) 
               for ( int k = 0; k != 2; ++k ) {
                  BlkC::ILoc_t iloc = ( i ) ? BlkC::IMax : BlkC::IMin;
                  BlkC::JLoc_t jloc = ( j ) ? BlkC::JMax : BlkC::JMin;
                  BlkC::KLoc_t kloc = ( k ) ? BlkC::KMax : BlkC::KMin;
                  
                  if(iloc == BlkC::IMin){
                     id = Input.Nghost;
                  }else{
                     id = Input.NCells_Idir + Input.Nghost;
                  }
                  if(jloc == BlkC::JMin){
                     jd =  Input.Nghost;
                  }else{
                     jd = Input.NCells_Jdir + Input.Nghost;
                  }
                  if(kloc == BlkC::KMin){
                     kd =  Input.Nghost;
                  }else{
                     kd = Input.NCells_Kdir + Input.Nghost;
                  }
                  
                  // Set coordinate
                  vp.set_coord(iloc, jloc, kloc, 
                               Grid_Blks[iblk].Node[id][jd][kd].X.x,
                               Grid_Blks[iblk].Node[id][jd][kd].X.y,  
                               Grid_Blks[iblk].Node[id][jd][kd].X.z);
                  
                  // Set vectors along each edge
                  if (id == Input.Nghost){
                     id_n = id +1;
                  }else{
                     id_n = id - 1;
                  }
                  if (jd == Input.Nghost){
                     jd_n = jd +1;
                  }else{
                     jd_n = jd - 1;
                  }
                  if (kd == Input.Nghost){
                     kd_n = kd +1;
                  }else{
                     kd_n = kd - 1;
                  }
                  dX_i = Grid_Blks[iblk].Node[id_n][jd][kd].X - Grid_Blks[iblk].Node[id][jd][kd].X;
                  dX_j = Grid_Blks[iblk].Node[id][jd_n][kd].X - Grid_Blks[iblk].Node[id][jd][kd].X;
                  dX_k = Grid_Blks[iblk].Node[id][jd][kd_n].X - Grid_Blks[iblk].Node[id][jd][kd].X;
                                
                  vp.set_vector(iloc, jloc, kloc, 'i', dX_i.x, dX_i.y, dX_i.z);
                  vp.set_vector(iloc, jloc, kloc, 'j', dX_j.x, dX_j.y, dX_j.z);
                  vp.set_vector(iloc, jloc, kloc, 'k', dX_k.x, dX_k.y, dX_k.z);
               }
               
         blkConn.add_block(iblk, Input.NCells_Idir, Input.NCells_Jdir, Input.NCells_Kdir, vp);
                  
               
      } /* endfor */
      
      
      // get the number of blocks on each element (total 26 elements) for each block
      // and the neighbour information
      for (int iblk = 0 ; iblk <NBlk ; ++iblk ){
       
         // number of blocks on the top face
         Connectivity[iblk].num_neighT = blkConn.num_neighbour_block(
            iblk, BlkC::IAll, BlkC::JAll, BlkC::KMax);
         // information of neighbour block of the top face
         Connectivity[iblk].neighT_info.set_block_orientation_info(
            blkConn, iblk, BlkC::IAll, BlkC::JAll, BlkC::KMax, 
            Connectivity[iblk].neighT);

         Connectivity[iblk].num_neighB = blkConn.num_neighbour_block(
            iblk,  BlkC::IAll, BlkC::JAll, BlkC::KMin);
         // information of neighbour block of the bottom face
         Connectivity[iblk].neighB_info.set_block_orientation_info(
            blkConn, iblk,  BlkC::IAll, BlkC::JAll, BlkC::KMin, 
            Connectivity[iblk].neighB);

         Connectivity[iblk].num_neighN = blkConn.num_neighbour_block(
            iblk,  BlkC::IAll, BlkC::JMax, BlkC::KAll);
         // information of neighbour block of the north face
         Connectivity[iblk].neighN_info.set_block_orientation_info(
            blkConn, iblk,  BlkC::IAll, BlkC::JMax, BlkC::KAll, 
            Connectivity[iblk].neighN);

         Connectivity[iblk].num_neighS = blkConn.num_neighbour_block(
            iblk,   BlkC::IAll, BlkC::JMin, BlkC::KAll);
         // information of neighbour block of the south face
         Connectivity[iblk].neighS_info.set_block_orientation_info(
            blkConn, iblk,   BlkC::IAll, BlkC::JMin, BlkC::KAll, 
            Connectivity[iblk].neighS);

         Connectivity[iblk].num_neighE = blkConn.num_neighbour_block(
            iblk, BlkC::IMax, BlkC::JAll, BlkC::KAll);
         // information of neighbour block of the east face
         Connectivity[iblk].neighE_info.set_block_orientation_info(
            blkConn, iblk, BlkC::IMax, BlkC::JAll, BlkC::KAll, 
            Connectivity[iblk].neighE);
         
         Connectivity[iblk].num_neighW = blkConn.num_neighbour_block(
            iblk,  BlkC::IMin, BlkC::JAll, BlkC::KAll);
         // information of neighbour block of the west face
         Connectivity[iblk].neighW_info.set_block_orientation_info(
            blkConn, iblk,  BlkC::IMin, BlkC::JAll, BlkC::KAll, 
            Connectivity[iblk].neighW);
          

         Connectivity[iblk].num_neighTN = blkConn.num_neighbour_block(
            iblk, BlkC::IAll, BlkC::JMax, BlkC::KMax);
         // information of neighbour block of the topnorth edge
         Connectivity[iblk].neighTN_info[0].set_block_orientation_info(
            blkConn, iblk, BlkC::IAll, BlkC::JMax, BlkC::KMax, 
            Connectivity[iblk].neighTN[0]);

         Connectivity[iblk].num_neighTS = blkConn.num_neighbour_block(
            iblk, BlkC::IAll, BlkC::JMin, BlkC::KMax);
         // information of neighbour block of the topsouth edge
         Connectivity[iblk].neighTS_info[0].set_block_orientation_info(
            blkConn, iblk,  BlkC::IAll, BlkC::JMin, BlkC::KMax, 
            Connectivity[iblk].neighTS[0]);
         
         Connectivity[iblk].num_neighTW = blkConn.num_neighbour_block(
            iblk,  BlkC::IMin, BlkC::JAll, BlkC::KMax);
         // information of neighbour block of the topwest edge
         Connectivity[iblk].neighTW_info[0].set_block_orientation_info(
            blkConn, iblk,  BlkC::IMin, BlkC::JAll, BlkC::KMax, 
            Connectivity[iblk].neighTW[0]);

         Connectivity[iblk].num_neighTE = blkConn.num_neighbour_block(
            iblk,  BlkC::IMax, BlkC::JAll, BlkC::KMax);
         // information of neighbour block of the topeast edge
         Connectivity[iblk].neighTE_info[0].set_block_orientation_info(
            blkConn, iblk,  BlkC::IMax, BlkC::JAll, BlkC::KMax, 
            Connectivity[iblk].neighTE[0]);

         Connectivity[iblk].num_neighBN = blkConn.num_neighbour_block(
            iblk,   BlkC::IAll, BlkC::JMax, BlkC::KMin);
         // information of neighbour block of the bottomnorth edge 
         Connectivity[iblk].neighBN_info[0].set_block_orientation_info(
            blkConn, iblk,  BlkC::IAll, BlkC::JMax, BlkC::KMin, 
            Connectivity[iblk].neighBN[0]);
         
         Connectivity[iblk].num_neighBS = blkConn.num_neighbour_block(
            iblk,   BlkC::IAll, BlkC::JMin, BlkC::KMin);
         // information of neighbour block of the bottomsouth edge 
         Connectivity[iblk].neighBS_info[0].set_block_orientation_info(
            blkConn, iblk,  BlkC::IAll, BlkC::JMin, BlkC::KMin, 
            Connectivity[iblk].neighBS[0]);

         Connectivity[iblk].num_neighBW = blkConn.num_neighbour_block(
            iblk,   BlkC::IMin, BlkC::JAll, BlkC::KMin);
         // information of neighbour block of the bottomwest edge  
         Connectivity[iblk].neighBW_info[0].set_block_orientation_info(
            blkConn, iblk,   BlkC::IMin, BlkC::JAll, BlkC::KMin, 
            Connectivity[iblk].neighBW[0]);

         Connectivity[iblk].num_neighBE = blkConn.num_neighbour_block(
            iblk,  BlkC::IMax, BlkC::JAll, BlkC::KMin);
         // information of neighbour block of the bottomeast edge
         Connectivity[iblk].neighBE_info[0].set_block_orientation_info(
            blkConn, iblk,  BlkC::IMax, BlkC::JAll, BlkC::KMin, 
            Connectivity[iblk].neighBE[0]);
         
         Connectivity[iblk].num_neighNW = blkConn.num_neighbour_block(
            iblk, BlkC::IMin, BlkC::JMax, BlkC::KAll);
         // information of neighbour block of the northwest edge
         Connectivity[iblk].neighNW_info[0].set_block_orientation_info(
            blkConn, iblk, BlkC::IMin, BlkC::JMax, BlkC::KAll, 
            Connectivity[iblk].neighNW[0]);
         
         Connectivity[iblk].num_neighNE = blkConn.num_neighbour_block(
            iblk,  BlkC::IMax, BlkC::JMax, BlkC::KAll);
         // information of neighbour block of the northeast edge
         Connectivity[iblk].neighNE_info[0].set_block_orientation_info(
            blkConn, iblk, BlkC::IMax, BlkC::JMax, BlkC::KAll, 
            Connectivity[iblk].neighNE[0]);

         Connectivity[iblk].num_neighSE = blkConn.num_neighbour_block(
            iblk,  BlkC::IMax, BlkC::JMin, BlkC::KAll);
         // information of neighbour block of the southeast edge 
         Connectivity[iblk].neighSE_info[0].set_block_orientation_info(
            blkConn, iblk,  BlkC::IMax, BlkC::JMin, BlkC::KAll, 
            Connectivity[iblk].neighSE[0]);
                
         Connectivity[iblk].num_neighSW = blkConn.num_neighbour_block(
            iblk,   BlkC::IMin, BlkC::JMin, BlkC::KAll);
         // information of neighbour block of the southwest edge 
         Connectivity[iblk].neighSW_info[0].set_block_orientation_info(
            blkConn, iblk,  BlkC::IMin, BlkC::JMin, BlkC::KAll, 
            Connectivity[iblk].neighSW[0]);

         Connectivity[iblk].num_neighTNW = blkConn.num_neighbour_block(
            iblk,  BlkC::IMin, BlkC::JMax, BlkC::KMax);
         // information of neighbour block of the topnorthwest vertex
         Connectivity[iblk].neighTNW_info[0].set_block_orientation_info(
            blkConn, iblk,  BlkC::IMin, BlkC::JMax, BlkC::KMax, 
            Connectivity[iblk].neighTNW[0]);

         Connectivity[iblk].num_neighTSW = blkConn.num_neighbour_block(
            iblk, BlkC::IMin, BlkC::JMin, BlkC::KMax);
         // information of neighbour block of the topsouthwest vertex
         Connectivity[iblk].neighTSW_info[0].set_block_orientation_info(
            blkConn, iblk,  BlkC::IMin, BlkC::JMin, BlkC::KMax, 
            Connectivity[iblk].neighTSW[0]);
         
         Connectivity[iblk].num_neighTNE = blkConn.num_neighbour_block(
            iblk,  BlkC::IMax, BlkC::JMax, BlkC::KMax);
         // information of neighbour block of the topnortheast vertex
         Connectivity[iblk].neighTNE_info[0].set_block_orientation_info(
            blkConn, iblk,  BlkC::IMax, BlkC::JMax, BlkC::KMax, 
            Connectivity[iblk].neighTNE[0]);

         Connectivity[iblk].num_neighTSE = blkConn.num_neighbour_block(
            iblk, BlkC::IMax, BlkC::JMin, BlkC::KMax);
         // information of neighbour block of the topsoutheast vertex
         Connectivity[iblk].neighTSE_info[0].set_block_orientation_info(
            blkConn, iblk, BlkC::IMax, BlkC::JMin, BlkC::KMax, 
            Connectivity[iblk].neighTSE[0]);

        Connectivity[iblk].num_neighBNW = blkConn.num_neighbour_block(
            iblk,  BlkC::IMin, BlkC::JMax, BlkC::KMin);
         // information of neighbour block of the bottomnrothwest vertex
         Connectivity[iblk].neighBNW_info[0].set_block_orientation_info(
            blkConn, iblk,  BlkC::IMin, BlkC::JMax, BlkC::KMin, 
            Connectivity[iblk].neighBNW[0]);

         Connectivity[iblk].num_neighBSW = blkConn.num_neighbour_block(
            iblk,  BlkC::IMin, BlkC::JMin, BlkC::KMin);
         // information of neighbour block of the bottomsouthwest vertex
         Connectivity[iblk].neighBSW_info[0].set_block_orientation_info(
            blkConn, iblk,  BlkC::IMin, BlkC::JMin, BlkC::KMin, 
            Connectivity[iblk].neighBSW[0]);

         Connectivity[iblk].num_neighBNE = blkConn.num_neighbour_block(
            iblk,  BlkC::IMax, BlkC::JMax, BlkC::KMin);
         // information of neighbour block of the bottomnortheast vertex
         Connectivity[iblk].neighBNE_info[0].set_block_orientation_info(
            blkConn, iblk,  BlkC::IMax, BlkC::JMax, BlkC::KMin, 
            Connectivity[iblk].neighBNE[0]);

         Connectivity[iblk].num_neighBSE = blkConn.num_neighbour_block(
            iblk,  BlkC::IMax, BlkC::JMin, BlkC::KMin);
         // information of neighbour block of the bottomsoutheast vertex
         Connectivity[iblk].neighBSE_info[0].set_block_orientation_info(
            blkConn, iblk, BlkC::IMax, BlkC::JMin, BlkC::KMin, 
            Connectivity[iblk].neighBSE[0]);

         // for roots only 
         if(Connectivity[iblk].num_neighT>ONE || Connectivity[iblk].num_neighB>ONE ||
            Connectivity[iblk].num_neighN>ONE || Connectivity[iblk].num_neighS>ONE ||
            Connectivity[iblk].num_neighE>ONE || Connectivity[iblk].num_neighW>ONE || 
            Connectivity[iblk].num_neighTN>BLOCK_ORIENTATION_MAX_NEIGHBOUR ||
            Connectivity[iblk].num_neighTS>BLOCK_ORIENTATION_MAX_NEIGHBOUR||
            Connectivity[iblk].num_neighTE>BLOCK_ORIENTATION_MAX_NEIGHBOUR||
            Connectivity[iblk].num_neighTW>BLOCK_ORIENTATION_MAX_NEIGHBOUR||
            Connectivity[iblk].num_neighTNW>BLOCK_ORIENTATION_MAX_NEIGHBOUR||
            Connectivity[iblk].num_neighTNE>BLOCK_ORIENTATION_MAX_NEIGHBOUR||
            Connectivity[iblk].num_neighTSW>BLOCK_ORIENTATION_MAX_NEIGHBOUR||
            Connectivity[iblk].num_neighTSE>BLOCK_ORIENTATION_MAX_NEIGHBOUR||
            Connectivity[iblk].num_neighBN>BLOCK_ORIENTATION_MAX_NEIGHBOUR ||
            Connectivity[iblk].num_neighBS>BLOCK_ORIENTATION_MAX_NEIGHBOUR||
            Connectivity[iblk].num_neighBE>BLOCK_ORIENTATION_MAX_NEIGHBOUR||
            Connectivity[iblk].num_neighBW>BLOCK_ORIENTATION_MAX_NEIGHBOUR||
            Connectivity[iblk].num_neighBNW>BLOCK_ORIENTATION_MAX_NEIGHBOUR||
            Connectivity[iblk].num_neighBNE>BLOCK_ORIENTATION_MAX_NEIGHBOUR||
            Connectivity[iblk].num_neighBSW>BLOCK_ORIENTATION_MAX_NEIGHBOUR||
            Connectivity[iblk].num_neighBSE>BLOCK_ORIENTATION_MAX_NEIGHBOUR){
            
            cerr<<"\n Number of neighbour blocks on some boudnary elements is out of bounds. \n";
            exit(1);

         }//
      
      } /* endfor */
  
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

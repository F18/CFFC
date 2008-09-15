/* Grid3DHexaMultiBlock.cc:  Member functions for 
 3D hexahedral multiblock grid classes. */

/* Include 3D hexahedral block grid header file. */

#ifndef GRID3D_HEXA_MULTIBLOCK_INCLUDED
#include "Grid3DHexaMultiBlock.h"
#endif //_GRID3D_HEXA_MULTIBLOCK_INCLUDED

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
    
    MPI::COMM_WORLD.Bcast(neighTN, GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(neighTS, GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(neighTW, GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(neighTE, GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(neighBN, GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(neighBS, GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(neighBE, GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(neighBW, GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(neighSW, GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(neighSE, GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(neighNW, GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(neighNE, GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS, MPI::INT, 0);
    
    MPI::COMM_WORLD.Bcast(neighTNW, GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(neighTNE, GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(neighTSW, GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(neighTSE, GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(neighBNW, GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(neighBNE, GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(neighBSW, GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(neighBSE, GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS, MPI::INT, 0);
    
    neighT_info.broadcast();
    neighB_info.broadcast();
    neighN_info.broadcast();
    neighS_info.broadcast();
    neighE_info.broadcast();
    neighW_info.broadcast();
    
    for (int i_neigh = 0; i_neigh < GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS; ++i_neigh) {
        neighNW_info[i_neigh].broadcast();
        neighNE_info[i_neigh].broadcast();
        neighSW_info[i_neigh].broadcast();
        neighSE_info[i_neigh].broadcast();
        neighTN_info[i_neigh].broadcast();
        neighTE_info[i_neigh].broadcast();
        neighTW_info[i_neigh].broadcast();
        neighTS_info[i_neigh].broadcast();
        neighBN_info[i_neigh].broadcast();
        neighBE_info[i_neigh].broadcast();
        neighBW_info[i_neigh].broadcast();
        neighBS_info[i_neigh].broadcast();
        
        neighTNW_info[i_neigh].broadcast();
        neighTNE_info[i_neigh].broadcast();
        neighTSW_info[i_neigh].broadcast();
        neighTSE_info[i_neigh].broadcast();
        neighBNW_info[i_neigh].broadcast();
        neighBNE_info[i_neigh].broadcast();
        neighBSW_info[i_neigh].broadcast();
        neighBSE_info[i_neigh].broadcast();
    } /* endfor */
    
    be.broadcast();
#endif  
    
}

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
        
        IsAuxiliary = Grid2.IsAuxiliary;
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
    MPI::COMM_WORLD.Bcast(&IsAuxiliary, 1, MPI::INT, 0);

    
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
        } else {
            Deallocate();
        } /* endif */
    } /* endif */
    
    /* Broadcast each of the blocks in the multiblock mesh. */
    /* Broadcast block connectivity info in the multiblock mesh. */
    if (Allocated) {
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
    
    // Create various multiblock multiblock grid depending on input parameters
    switch(Input.i_Grid) {
        case GRID_CUBE :
            Create_Grid_Cube(Input);
            break;
            
        case GRID_PERIODIC_BOX :
        case GRID_PERIODIC_BOX_WITH_INFLOW :
            Create_Grid_Periodic_Box(Input);
            break;
            
        case GRID_TURBULENCE_BOX :
            Create_Grid_Turbulence_Box(Input);
            break;
            
        case GRID_BUNSEN_INFLOW :
            Create_Grid_Bunsen_Inflow(Input);
            break;
            
        case GRID_BUNSEN_BOX :
            Create_Grid_Bunsen_Box(Input);
            break;
            
        case GRID_BUNSEN_BURNER :
            Create_Grid_Bunsen_Burner(Input);
            break;
            
        case GRID_FLAT_PLATE :
            Create_Grid_Flat_Plate(Input);
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
            
        case GRID_BUMP_CHANNEL_FLOW :
            Create_Grid_Bump_Channel_Flow(Input);
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
    
    if (Input.Disturb_Interior_Nodes != 0 && !IsAuxiliary) {
        Disturb_Interior_Nodes(Input.Disturb_Interior_Nodes);
    }
    
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
                                              Input.Mesh_Stretching,
                                              Input.Stretching_Type_Idir,
                                              Input.Stretching_Type_Jdir,
                                              Input.Stretching_Factor_Idir,
                                              Input.Stretching_Factor_Jdir,
                                              Input.NCells_Idir,
                                              Input.NCells_Jdir,
                                              Input.Nghost);
    
    if (Input.Mesh_Smoothing!=OFF){
        /* This will make the mesh curved in case of mesh stretching */
        for (int iBlk=0; iBlk<Input.NBlk_Idir; iBlk++) {
            for (int jBlk=0; jBlk<Input.NBlk_Jdir; jBlk++) {
                Smooth_Quad_Block(Grid2D_Box_XYplane[iBlk][jBlk],
                                  Input.Mesh_Smoothing);
            }
        }
    }
    
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
                                                Input.Mesh_Stretching,
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
    
    
    /* Call the function Find_Neighbours to obtain the neighbour block information
     and assign values to data members in the grid block connectivity data structure. */
    
    Find_Neighbours(Input);
    
}


/********************************************************
 * Routine: Create_Grid_Single_Block_Periodic_Box       *
 *                                                      *
 * Generates a 3D Cartesian multiblock mesh for a cube. *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block_List::Create_Grid_Single_Block_Periodic_Box(Grid3D_Input_Parameters &Input) {
    
    int count_blocks;
    int BC_top, BC_bottom;
    Grid2D_Quad_Block **Grid2D_Box_XYplane;
    
    /* Allocate required memory. */
    
    assert(Input.NBlk_Idir==1 && Input.NBlk_Jdir==1 && Input.NBlk_Kdir==1);
    
    Allocate(Input.NBlk_Idir, Input.NBlk_Jdir, Input.NBlk_Kdir);
    
    /* Creat 2D cross-section grids from which the 3D grid
     will be extruded. */
    
    Grid2D_Box_XYplane = Grid_Rectangular_Box(Grid2D_Box_XYplane,
                                              Input.NBlk_Idir, 
                                              Input.NBlk_Jdir,
                                              Input.Box_Width,
                                              Input.Box_Height,
                                              Input.Mesh_Stretching,
                                              Input.Stretching_Type_Idir,
                                              Input.Stretching_Type_Jdir,
                                              Input.Stretching_Factor_Idir,
                                              Input.Stretching_Factor_Jdir,
                                              Input.NCells_Idir,
                                              Input.NCells_Jdir,
                                              Input.Nghost);
    
    if (Input.Mesh_Smoothing!=OFF){
        /* This will make the mesh curved in case of mesh stretching */
        for (int iBlk=0; iBlk<Input.NBlk_Idir; iBlk++) {
            for (int jBlk=0; jBlk<Input.NBlk_Jdir; jBlk++) {
                Smooth_Quad_Block(Grid2D_Box_XYplane[iBlk][jBlk],
                                  Input.Mesh_Smoothing);
            }
        }
    }
    
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
                                                Input.Mesh_Stretching,
                                                Input.Stretching_Type_Kdir,
                                                Input.Stretching_Factor_Kdir,
                                                -HALF*Input.Box_Length+
                                                (double(kBlk)/double(Input.NBlk_Kdir))*Input.Box_Length,
                                                -HALF*Input.Box_Length+
                                                (double(kBlk+1)/double(Input.NBlk_Kdir))*Input.Box_Length);
                
                /* Set all boundary conditions to periodic */
                
                Grid_Blks[count_blocks].Set_BCs(BC_PERIODIC, 
                                                BC_PERIODIC, 
                                                BC_PERIODIC, 
                                                BC_PERIODIC, 
                                                BC_PERIODIC, 
                                                BC_PERIODIC);
                                
                /* Update block counter. */
                
                count_blocks ++;
                
            } /* endfor */
        } /* endfor */
    } /* endfor */
    
    /* Deallocate 2D grid. */
    
    Grid2D_Box_XYplane = Deallocate_Multi_Block_Grid(Grid2D_Box_XYplane,
                                                     Input.NBlk_Idir, 
                                                     Input.NBlk_Jdir);
    
    
    /* Call the function Find_Neighbours to obtain the neighbour block information
     and assign values to data members in the grid block connectivity data structure. */
    
    Find_Neighbours(Input);
    
}

/********************************************************
 * Routine: Create_Grid_Periodic_Box                    *
 *                                                      *
 * Generates a 3D Cartesian multiblock mesh for a box   *
 * with periodic boundaries on all six boundaries.      *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block_List::Create_Grid_Periodic_Box(Grid3D_Input_Parameters &Input) {
    
    /* No block connectivity needed in case of 1 single block --> No Message Passing for ghost cells */
    if (Input.NBlk_Idir==1 && Input.NBlk_Jdir==1 && Input.NBlk_Kdir==1 && Input.i_Grid!=GRID_PERIODIC_BOX_WITH_INFLOW) {
        return Create_Grid_Single_Block_Periodic_Box(Input);
    }

    
    int nBlk;
    int opposite_nBlk, opposite_iBlk, opposite_jBlk, opposite_kBlk;
    Direction_Indices Dir_Index;
    int BC_top, BC_bottom;
    Grid2D_Quad_Block **Grid2D_Box_XYplane;
    
    /* Allocate required memory. */
    
    Allocate(Input.NBlk_Idir, Input.NBlk_Jdir, Input.NBlk_Kdir);
    
    /* Creat 2D cross-section grids from which the 3D grid
     will be extruded. */
    
    Grid2D_Box_XYplane = Grid_Periodic_Box(Grid2D_Box_XYplane,
                                           Input.NBlk_Idir, 
                                           Input.NBlk_Jdir,
                                           Input.Box_Width,
                                           Input.Box_Height,
                                           Input.Mesh_Stretching,
                                           Input.Stretching_Type_Idir,
                                           Input.Stretching_Type_Jdir,
                                           Input.Stretching_Factor_Idir,
                                           Input.Stretching_Factor_Jdir,
                                           Input.NCells_Idir,
                                           Input.NCells_Jdir,
                                           Input.Nghost);
    
    if (Input.Mesh_Smoothing!=OFF){
        /* This will make the mesh curved in case of mesh stretching */
        for (int iBlk=0; iBlk<Input.NBlk_Idir; iBlk++) {
            for (int jBlk=0; jBlk<Input.NBlk_Jdir; jBlk++) {
                Smooth_Quad_Block(Grid2D_Box_XYplane[iBlk][jBlk],
                                  Input.Mesh_Smoothing);
            }
        }
    }
    
    /* Create the mesh for each block representing
     the complete grid. */
    
    for (int kBlk = 0; kBlk <= Input.NBlk_Kdir-1; ++kBlk) {
        for (int jBlk = 0; jBlk <= Input.NBlk_Jdir-1; ++jBlk) {
            for (int iBlk = 0; iBlk <= Input.NBlk_Idir-1; ++iBlk) {
                
                /* Determine grid block number */
                
                nBlk = iBlk + 
                jBlk*Input.NBlk_Idir + 
                kBlk*Input.NBlk_Idir*Input.NBlk_Jdir;
                
                /* Extrude each of the grid blocks from the
                 appropriate 2D grid in XY-plane. */
                
                Grid_Blks[nBlk].Extrude(Grid2D_Box_XYplane[iBlk][jBlk],
                                        Input.NCells_Kdir,
                                        Input.Mesh_Stretching,
                                        Input.Stretching_Type_Kdir, 
                                        Input.Stretching_Factor_Kdir,
                                        -HALF*Input.Box_Length+
                                        (double(kBlk)/double(Input.NBlk_Kdir))*Input.Box_Length,
                                        -HALF*Input.Box_Length+
                                        (double(kBlk+1)/double(Input.NBlk_Kdir))*Input.Box_Length);
                
                /* Assign top and bottom boundary conditions. */
                
                if (kBlk == Input.NBlk_Kdir-1) {
                    BC_top = BC_NONE;
                } else {
                    BC_top = BC_NONE;
                } /* endif */
                if (kBlk == 0) {
                    BC_bottom = BC_NONE;
                } else {
                    BC_bottom = BC_NONE;
                } /* endif */
                
                Grid_Blks[nBlk].Set_BCs_Zdir(BC_top, BC_bottom);
                
            } /* endfor */
        } /* endfor */
    } /* endfor */
    
    /* Fix boundary conditions for case of inflow. */
    
    if (Input.i_Grid == GRID_PERIODIC_BOX_WITH_INFLOW) {
        for (int kBlk = 0; kBlk <= Input.NBlk_Kdir-1; ++kBlk) {
            for (int jBlk = 0; jBlk <= Input.NBlk_Jdir-1; ++jBlk) {
                for (int iBlk = 0; iBlk <= Input.NBlk_Idir-1; ++iBlk) {
                    nBlk = iBlk + 
                    jBlk*Input.NBlk_Idir + 
                    kBlk*Input.NBlk_Idir*Input.NBlk_Jdir;
                    if (iBlk == 0 && iBlk == Input.NBlk_Idir-1) {
                        Grid_Blks[nBlk].Set_BCs_Xdir(BC_OUTFLOW_SUBSONIC,
                                                     BC_INFLOW_SUBSONIC);
                    } else if (iBlk == 0) {
                        Grid_Blks[nBlk].Set_BCs_Xdir(BC_NONE,
                                                     BC_INFLOW_SUBSONIC);
                    } else if (iBlk == Input.NBlk_Idir-1) {
                        Grid_Blks[nBlk].Set_BCs_Xdir(BC_OUTFLOW_SUBSONIC,
                                                     BC_NONE);
                    } /* endif */
                } /* endfor */
            } /* endfor */
        } /* endfor */
    } /* endif */
    
    /* Deallocate 2D grid. */
    
    Grid2D_Box_XYplane = Deallocate_Multi_Block_Grid(Grid2D_Box_XYplane,
                                                     Input.NBlk_Idir, 
                                                     Input.NBlk_Jdir);
    
    /* Call the function Find_Neighbours to obtain the neighbour block information
     and assign values to data members in the grid block connectivity data structure. */
    
    Find_Neighbours(Input);
    
    /* At periodic boundaries, add additional grid block neighbour information. */
    
    for (int kBlk = 0; kBlk <= Input.NBlk_Kdir-1; ++kBlk) {
        for (int jBlk = 0; jBlk <= Input.NBlk_Jdir-1; ++jBlk) {
            for (int iBlk = 0; iBlk <= Input.NBlk_Idir-1; ++iBlk) {
                
                nBlk = iBlk + 
                jBlk*Input.NBlk_Idir + 
                kBlk*Input.NBlk_Idir*Input.NBlk_Jdir;
                
                for (int nDir = 0; nDir <= MAX_BOUNDARY_ELEMENTS_FOR_A_BLOCK-1; ++nDir) {
                    Dir_Index = Dir_Index.boundary_element_number_to_direction_indices(nDir);
                    
                    if (iBlk == 0 && Dir_Index.i < 0) {
                        opposite_iBlk = Input.NBlk_Idir-1;
                    } else if (iBlk == Input.NBlk_Idir-1 && Dir_Index.i > 0) {
                        opposite_iBlk = 0;
                    } else {
                        opposite_iBlk = iBlk + Dir_Index.i;
                    } /* endif */
                    
                    if (jBlk == 0 && Dir_Index.j < 0) {
                        opposite_jBlk = Input.NBlk_Jdir-1;
                    } else if (jBlk == Input.NBlk_Jdir-1 && Dir_Index.j > 0) {
                        opposite_jBlk = 0;
                    } else {
                        opposite_jBlk = jBlk + Dir_Index.j;
                    } /* endif */
                    
                    if (kBlk == 0 && Dir_Index.k < 0) {
                        opposite_kBlk = Input.NBlk_Kdir-1;
                    } else if (kBlk == Input.NBlk_Kdir-1 && Dir_Index.k > 0) {
                        opposite_kBlk = 0;
                    } else {
                        opposite_kBlk = kBlk + Dir_Index.k;
                    } /* endif */
                    
                    opposite_nBlk = opposite_iBlk + 
                    opposite_jBlk*Input.NBlk_Idir + 
                    opposite_kBlk*Input.NBlk_Idir*Input.NBlk_Jdir;
                    
                    switch (nDir) { 
                        case BE::BSW : // 0
                            if (Connectivity[nBlk].num_neighBSW == 0 && 
                                (Input.i_Grid != GRID_PERIODIC_BOX_WITH_INFLOW ||
                                 (Input.i_Grid == GRID_PERIODIC_BOX_WITH_INFLOW &&
                                  iBlk > 0))
                                ) {
                                Connectivity[nBlk].num_neighBSW = 1;
                                Connectivity[nBlk].neighBSW[0] = opposite_nBlk;
                                Connectivity[nBlk].neighBSW_info[0].ctm_offsets[0] = 1; 
                                Connectivity[nBlk].neighBSW_info[0].ctm_offsets[1] = 2;  
                                Connectivity[nBlk].neighBSW_info[0].ctm_offsets[2] = 3;
                                Connectivity[nBlk].neighBSW_info[0].ctm_offsets[3] = Input.NCells_Idir; 
                                Connectivity[nBlk].neighBSW_info[0].ctm_offsets[4] = Input.NCells_Jdir;  
                                Connectivity[nBlk].neighBSW_info[0].ctm_offsets[5] = Input.NCells_Kdir;  
                                Connectivity[nBlk].neighBSW_info[0].direction_me_to_neighbour[0] = -1; 
                                Connectivity[nBlk].neighBSW_info[0].direction_me_to_neighbour[1] = -1;
                                Connectivity[nBlk].neighBSW_info[0].direction_me_to_neighbour[2] = -1; 
                                Connectivity[nBlk].neighBSW_info[0].direction_neighbour_to_me[0] = 1; 
                                Connectivity[nBlk].neighBSW_info[0].direction_neighbour_to_me[1] = 1;
                                Connectivity[nBlk].neighBSW_info[0].direction_neighbour_to_me[2] = 1; 
                                Connectivity[nBlk].be.on_grid_boundary[BE::BSW] = 0;
                            } /* endif */
                            break;
                            
                            case BE::SW :  // 1
                            if (Connectivity[nBlk].num_neighSW == 0 && 
                                (Input.i_Grid != GRID_PERIODIC_BOX_WITH_INFLOW ||
                                 (Input.i_Grid == GRID_PERIODIC_BOX_WITH_INFLOW &&
                                  iBlk > 0))
                                ) {
                                Connectivity[nBlk].num_neighSW = 1;
                                Connectivity[nBlk].neighSW[0] = opposite_nBlk;
                                Connectivity[nBlk].neighSW_info[0].ctm_offsets[0] = 1; 
                                Connectivity[nBlk].neighSW_info[0].ctm_offsets[1] = 2;  
                                Connectivity[nBlk].neighSW_info[0].ctm_offsets[2] = 3;
                                Connectivity[nBlk].neighSW_info[0].ctm_offsets[3] = Input.NCells_Idir; 
                                Connectivity[nBlk].neighSW_info[0].ctm_offsets[4] = Input.NCells_Jdir;  
                                Connectivity[nBlk].neighSW_info[0].ctm_offsets[5] = 0;
                                Connectivity[nBlk].neighSW_info[0].direction_me_to_neighbour[0] = -1; 
                                Connectivity[nBlk].neighSW_info[0].direction_me_to_neighbour[1] = -1;
                                Connectivity[nBlk].neighSW_info[0].direction_me_to_neighbour[2] = 0; 
                                Connectivity[nBlk].neighSW_info[0].direction_neighbour_to_me[0] = 1; 
                                Connectivity[nBlk].neighSW_info[0].direction_neighbour_to_me[1] = 1;
                                Connectivity[nBlk].neighSW_info[0].direction_neighbour_to_me[2] = 0;
                                Connectivity[nBlk].be.on_grid_boundary[BE::SW] = 0;
                            } /* endif */
                            break;
                            
                            case BE::TSW : // 2
                            if (Connectivity[nBlk].num_neighTSW == 0 && 
                                (Input.i_Grid != GRID_PERIODIC_BOX_WITH_INFLOW ||
                                 (Input.i_Grid == GRID_PERIODIC_BOX_WITH_INFLOW &&
                                  iBlk > 0))
                                ) {
                                Connectivity[nBlk].num_neighTSW = 1;
                                Connectivity[nBlk].neighTSW[0] = opposite_nBlk;
                                Connectivity[nBlk].neighTSW_info[0].ctm_offsets[0] = 1; 
                                Connectivity[nBlk].neighTSW_info[0].ctm_offsets[1] = 2;  
                                Connectivity[nBlk].neighTSW_info[0].ctm_offsets[2] = 3;
                                Connectivity[nBlk].neighTSW_info[0].ctm_offsets[3] = Input.NCells_Idir; 
                                Connectivity[nBlk].neighTSW_info[0].ctm_offsets[4] = Input.NCells_Jdir;  
                                Connectivity[nBlk].neighTSW_info[0].ctm_offsets[5] = -Input.NCells_Kdir;
                                Connectivity[nBlk].neighTSW_info[0].direction_me_to_neighbour[0] = -1; 
                                Connectivity[nBlk].neighTSW_info[0].direction_me_to_neighbour[1] = -1;
                                Connectivity[nBlk].neighTSW_info[0].direction_me_to_neighbour[2] = 1; 
                                Connectivity[nBlk].neighTSW_info[0].direction_neighbour_to_me[0] = 1; 
                                Connectivity[nBlk].neighTSW_info[0].direction_neighbour_to_me[1] = 1;
                                Connectivity[nBlk].neighTSW_info[0].direction_neighbour_to_me[2] = -1;
                                Connectivity[nBlk].be.on_grid_boundary[BE::TSW] = 0;
                            } /* endif */
                            break;
                            
                            case BE::BW :  // 3
                            if (Connectivity[nBlk].num_neighBW == 0 && 
                                (Input.i_Grid != GRID_PERIODIC_BOX_WITH_INFLOW ||
                                 (Input.i_Grid == GRID_PERIODIC_BOX_WITH_INFLOW &&
                                  iBlk > 0))
                                ) {
                                Connectivity[nBlk].num_neighBW = 1;
                                Connectivity[nBlk].neighBW[0] = opposite_nBlk;
                                Connectivity[nBlk].neighBW_info[0].ctm_offsets[0] = 1; 
                                Connectivity[nBlk].neighBW_info[0].ctm_offsets[1] = 2;  
                                Connectivity[nBlk].neighBW_info[0].ctm_offsets[2] = 3;
                                Connectivity[nBlk].neighBW_info[0].ctm_offsets[3] = Input.NCells_Idir; 
                                Connectivity[nBlk].neighBW_info[0].ctm_offsets[4] = 0;  
                                Connectivity[nBlk].neighBW_info[0].ctm_offsets[5] = Input.NCells_Kdir;
                                Connectivity[nBlk].neighBW_info[0].direction_me_to_neighbour[0] = -1; 
                                Connectivity[nBlk].neighBW_info[0].direction_me_to_neighbour[1] = 0;
                                Connectivity[nBlk].neighBW_info[0].direction_me_to_neighbour[2] = -1; 
                                Connectivity[nBlk].neighBW_info[0].direction_neighbour_to_me[0] = 1; 
                                Connectivity[nBlk].neighBW_info[0].direction_neighbour_to_me[1] = 0;
                                Connectivity[nBlk].neighBW_info[0].direction_neighbour_to_me[2] = 1;
                                Connectivity[nBlk].be.on_grid_boundary[BE::BW] = 0;
                            } /* endif */
                            break;
                            
                            case BE::W :   // 4
                            if (Connectivity[nBlk].num_neighW == 0 && 
                                (Input.i_Grid != GRID_PERIODIC_BOX_WITH_INFLOW ||
                                 (Input.i_Grid == GRID_PERIODIC_BOX_WITH_INFLOW &&
                                  iBlk > 0))
                                ) {
                                Connectivity[nBlk].num_neighW = 1;
                                Connectivity[nBlk].neighW = opposite_nBlk;
                                Connectivity[nBlk].neighW_info.ctm_offsets[0] = 1; 
                                Connectivity[nBlk].neighW_info.ctm_offsets[1] = 2;  
                                Connectivity[nBlk].neighW_info.ctm_offsets[2] = 3;
                                Connectivity[nBlk].neighW_info.ctm_offsets[3] = Input.NCells_Idir; 
                                Connectivity[nBlk].neighW_info.ctm_offsets[4] = 0;  
                                Connectivity[nBlk].neighW_info.ctm_offsets[5] = 0;
                                Connectivity[nBlk].neighW_info.direction_me_to_neighbour[0] = -1; 
                                Connectivity[nBlk].neighW_info.direction_me_to_neighbour[1] = 0;
                                Connectivity[nBlk].neighW_info.direction_me_to_neighbour[2] = 0; 
                                Connectivity[nBlk].neighW_info.direction_neighbour_to_me[0] = 1; 
                                Connectivity[nBlk].neighW_info.direction_neighbour_to_me[1] = 0;
                                Connectivity[nBlk].neighW_info.direction_neighbour_to_me[2] = 0;
                                Connectivity[nBlk].be.on_grid_boundary[BE::W] = 0;
                            } /* endif */
                            break;
                            
                            case BE::TW :  // 5
                            if (Connectivity[nBlk].num_neighTW == 0 && 
                                (Input.i_Grid != GRID_PERIODIC_BOX_WITH_INFLOW ||
                                 (Input.i_Grid == GRID_PERIODIC_BOX_WITH_INFLOW &&
                                  iBlk > 0))
                                ) {
                                Connectivity[nBlk].num_neighTW = 1;
                                Connectivity[nBlk].neighTW[0] = opposite_nBlk;
                                Connectivity[nBlk].neighTW_info[0].ctm_offsets[0] = 1; 
                                Connectivity[nBlk].neighTW_info[0].ctm_offsets[1] = 2;  
                                Connectivity[nBlk].neighTW_info[0].ctm_offsets[2] = 3;
                                Connectivity[nBlk].neighTW_info[0].ctm_offsets[3] = Input.NCells_Idir; 
                                Connectivity[nBlk].neighTW_info[0].ctm_offsets[4] = 0;  
                                Connectivity[nBlk].neighTW_info[0].ctm_offsets[5] = -Input.NCells_Kdir;
                                Connectivity[nBlk].neighTW_info[0].direction_me_to_neighbour[0] = -1; 
                                Connectivity[nBlk].neighTW_info[0].direction_me_to_neighbour[1] = 0;
                                Connectivity[nBlk].neighTW_info[0].direction_me_to_neighbour[2] = 1; 
                                Connectivity[nBlk].neighTW_info[0].direction_neighbour_to_me[0] = 1; 
                                Connectivity[nBlk].neighTW_info[0].direction_neighbour_to_me[1] = 0;
                                Connectivity[nBlk].neighTW_info[0].direction_neighbour_to_me[2] = -1;
                                Connectivity[nBlk].be.on_grid_boundary[BE::TW] = 0;
                            } /* endif */
                            break;
                            
                            case BE::BNW : // 6
                            if (Connectivity[nBlk].num_neighBNW == 0 && 
                                (Input.i_Grid != GRID_PERIODIC_BOX_WITH_INFLOW ||
                                 (Input.i_Grid == GRID_PERIODIC_BOX_WITH_INFLOW &&
                                  iBlk > 0))
                                ) {
                                Connectivity[nBlk].num_neighBNW = 1;
                                Connectivity[nBlk].neighBNW[0] = opposite_nBlk;
                                Connectivity[nBlk].neighBNW_info[0].ctm_offsets[0] = 1; 
                                Connectivity[nBlk].neighBNW_info[0].ctm_offsets[1] = 2;  
                                Connectivity[nBlk].neighBNW_info[0].ctm_offsets[2] = 3;
                                Connectivity[nBlk].neighBNW_info[0].ctm_offsets[3] = Input.NCells_Idir; 
                                Connectivity[nBlk].neighBNW_info[0].ctm_offsets[4] = -Input.NCells_Jdir;  
                                Connectivity[nBlk].neighBNW_info[0].ctm_offsets[5] = Input.NCells_Kdir;  
                                Connectivity[nBlk].neighBNW_info[0].direction_me_to_neighbour[0] = -1; 
                                Connectivity[nBlk].neighBNW_info[0].direction_me_to_neighbour[1] = 1;
                                Connectivity[nBlk].neighBNW_info[0].direction_me_to_neighbour[2] = -1; 
                                Connectivity[nBlk].neighBNW_info[0].direction_neighbour_to_me[0] = 1; 
                                Connectivity[nBlk].neighBNW_info[0].direction_neighbour_to_me[1] = -1;
                                Connectivity[nBlk].neighBNW_info[0].direction_neighbour_to_me[2] = 1;
                                Connectivity[nBlk].be.on_grid_boundary[BE::BNW] = 0;
                            } /* endif */
                            break;
                            
                            case BE::NW :  // 7
                            if (Connectivity[nBlk].num_neighNW == 0 && 
                                (Input.i_Grid != GRID_PERIODIC_BOX_WITH_INFLOW ||
                                 (Input.i_Grid == GRID_PERIODIC_BOX_WITH_INFLOW &&
                                  iBlk > 0))
                                ) {
                                Connectivity[nBlk].num_neighNW = 1;
                                Connectivity[nBlk].neighNW[0] = opposite_nBlk;
                                Connectivity[nBlk].neighNW_info[0].ctm_offsets[0] = 1; 
                                Connectivity[nBlk].neighNW_info[0].ctm_offsets[1] = 2;  
                                Connectivity[nBlk].neighNW_info[0].ctm_offsets[2] = 3;
                                Connectivity[nBlk].neighNW_info[0].ctm_offsets[3] = Input.NCells_Idir; 
                                Connectivity[nBlk].neighNW_info[0].ctm_offsets[4] = -Input.NCells_Jdir;  
                                Connectivity[nBlk].neighNW_info[0].ctm_offsets[5] = 0;  
                                Connectivity[nBlk].neighNW_info[0].direction_me_to_neighbour[0] = -1; 
                                Connectivity[nBlk].neighNW_info[0].direction_me_to_neighbour[1] = 1;
                                Connectivity[nBlk].neighNW_info[0].direction_me_to_neighbour[2] = 0; 
                                Connectivity[nBlk].neighNW_info[0].direction_neighbour_to_me[0] = 1; 
                                Connectivity[nBlk].neighNW_info[0].direction_neighbour_to_me[1] = -1;
                                Connectivity[nBlk].neighNW_info[0].direction_neighbour_to_me[2] = 0;
                                Connectivity[nBlk].be.on_grid_boundary[BE::NW] = 0;
                            } /* endif */
                            break;
                            
                            case BE::TNW : // 8
                            if (Connectivity[nBlk].num_neighTNW == 0 && 
                                (Input.i_Grid != GRID_PERIODIC_BOX_WITH_INFLOW ||
                                 (Input.i_Grid == GRID_PERIODIC_BOX_WITH_INFLOW &&
                                  iBlk > 0))
                                ) {
                                Connectivity[nBlk].num_neighTNW = 1;
                                Connectivity[nBlk].neighTNW[0] = opposite_nBlk;
                                Connectivity[nBlk].neighTNW_info[0].ctm_offsets[0] = 1; 
                                Connectivity[nBlk].neighTNW_info[0].ctm_offsets[1] = 2;  
                                Connectivity[nBlk].neighTNW_info[0].ctm_offsets[2] = 3;
                                Connectivity[nBlk].neighTNW_info[0].ctm_offsets[3] = Input.NCells_Idir; 
                                Connectivity[nBlk].neighTNW_info[0].ctm_offsets[4] = -Input.NCells_Jdir;  
                                Connectivity[nBlk].neighTNW_info[0].ctm_offsets[5] = -Input.NCells_Kdir;  
                                Connectivity[nBlk].neighTNW_info[0].direction_me_to_neighbour[0] = -1; 
                                Connectivity[nBlk].neighTNW_info[0].direction_me_to_neighbour[1] = 1;
                                Connectivity[nBlk].neighTNW_info[0].direction_me_to_neighbour[2] = 1; 
                                Connectivity[nBlk].neighTNW_info[0].direction_neighbour_to_me[0] = 1; 
                                Connectivity[nBlk].neighTNW_info[0].direction_neighbour_to_me[1] = -1;
                                Connectivity[nBlk].neighTNW_info[0].direction_neighbour_to_me[2] = -1;
                                Connectivity[nBlk].be.on_grid_boundary[BE::TNW] = 0;
                            } /* endif */
                            break;
                            
                            case BE::BS :  // 9
                            if (Connectivity[nBlk].num_neighBS == 0) {
                                Connectivity[nBlk].num_neighBS = 1;
                                Connectivity[nBlk].neighBS[0] = opposite_nBlk;
                                Connectivity[nBlk].neighBS_info[0].ctm_offsets[0] = 1; 
                                Connectivity[nBlk].neighBS_info[0].ctm_offsets[1] = 2;  
                                Connectivity[nBlk].neighBS_info[0].ctm_offsets[2] = 3;
                                Connectivity[nBlk].neighBS_info[0].ctm_offsets[3] = 0; 
                                Connectivity[nBlk].neighBS_info[0].ctm_offsets[4] = Input.NCells_Jdir;  
                                Connectivity[nBlk].neighBS_info[0].ctm_offsets[5] = Input.NCells_Kdir;  
                                Connectivity[nBlk].neighBS_info[0].direction_me_to_neighbour[0] = 0; 
                                Connectivity[nBlk].neighBS_info[0].direction_me_to_neighbour[1] = -1;
                                Connectivity[nBlk].neighBS_info[0].direction_me_to_neighbour[2] = -1; 
                                Connectivity[nBlk].neighBS_info[0].direction_neighbour_to_me[0] = 0; 
                                Connectivity[nBlk].neighBS_info[0].direction_neighbour_to_me[1] = 1;
                                Connectivity[nBlk].neighBS_info[0].direction_neighbour_to_me[2] = 1;
                                Connectivity[nBlk].be.on_grid_boundary[BE::BS] = 0;
                            } /* endif */
                            break;
                            
                            case BE::S :   // 10
                            if (Connectivity[nBlk].num_neighS == 0) {
                                Connectivity[nBlk].num_neighS = 1;
                                Connectivity[nBlk].neighS = opposite_nBlk;
                                Connectivity[nBlk].neighS_info.ctm_offsets[0] = 1; 
                                Connectivity[nBlk].neighS_info.ctm_offsets[1] = 2;  
                                Connectivity[nBlk].neighS_info.ctm_offsets[2] = 3;
                                Connectivity[nBlk].neighS_info.ctm_offsets[3] = 0; 
                                Connectivity[nBlk].neighS_info.ctm_offsets[4] = Input.NCells_Jdir;  
                                Connectivity[nBlk].neighS_info.ctm_offsets[5] = 0;
                                Connectivity[nBlk].neighS_info.direction_me_to_neighbour[0] = 0; 
                                Connectivity[nBlk].neighS_info.direction_me_to_neighbour[1] = -1;
                                Connectivity[nBlk].neighS_info.direction_me_to_neighbour[2] = 0; 
                                Connectivity[nBlk].neighS_info.direction_neighbour_to_me[0] = 0; 
                                Connectivity[nBlk].neighS_info.direction_neighbour_to_me[1] = 1;
                                Connectivity[nBlk].neighS_info.direction_neighbour_to_me[2] = 0;
                                Connectivity[nBlk].be.on_grid_boundary[BE::S] = 0;
                            } /* endif */
                            break;
                            
                            case BE::TS :  // 11
                            if (Connectivity[nBlk].num_neighTS == 0) {
                                Connectivity[nBlk].num_neighTS = 1;
                                Connectivity[nBlk].neighTS[0] = opposite_nBlk;
                                Connectivity[nBlk].neighTS_info[0].ctm_offsets[0] = 1; 
                                Connectivity[nBlk].neighTS_info[0].ctm_offsets[1] = 2;  
                                Connectivity[nBlk].neighTS_info[0].ctm_offsets[2] = 3;
                                Connectivity[nBlk].neighTS_info[0].ctm_offsets[3] = 0; 
                                Connectivity[nBlk].neighTS_info[0].ctm_offsets[4] = Input.NCells_Jdir;  
                                Connectivity[nBlk].neighTS_info[0].ctm_offsets[5] = -Input.NCells_Kdir;  
                                Connectivity[nBlk].neighTS_info[0].direction_me_to_neighbour[0] = 0; 
                                Connectivity[nBlk].neighTS_info[0].direction_me_to_neighbour[1] = -1;
                                Connectivity[nBlk].neighTS_info[0].direction_me_to_neighbour[2] = 1; 
                                Connectivity[nBlk].neighTS_info[0].direction_neighbour_to_me[0] = 0; 
                                Connectivity[nBlk].neighTS_info[0].direction_neighbour_to_me[1] = 1;
                                Connectivity[nBlk].neighTS_info[0].direction_neighbour_to_me[2] = -1;
                                Connectivity[nBlk].be.on_grid_boundary[BE::TS] = 0;
                            } /* endif */
                            break;
                            
                            case BE::B :   // 12
                            if (Connectivity[nBlk].num_neighB == 0) {
                                Connectivity[nBlk].num_neighB = 1;
                                Connectivity[nBlk].neighB = opposite_nBlk;
                                Connectivity[nBlk].neighB_info.ctm_offsets[0] = 1; 
                                Connectivity[nBlk].neighB_info.ctm_offsets[1] = 2;  
                                Connectivity[nBlk].neighB_info.ctm_offsets[2] = 3;
                                Connectivity[nBlk].neighB_info.ctm_offsets[3] = 0; 
                                Connectivity[nBlk].neighB_info.ctm_offsets[4] = 0;  
                                Connectivity[nBlk].neighB_info.ctm_offsets[5] = Input.NCells_Kdir;
                                Connectivity[nBlk].neighB_info.direction_me_to_neighbour[0] = 0; 
                                Connectivity[nBlk].neighB_info.direction_me_to_neighbour[1] = 0;
                                Connectivity[nBlk].neighB_info.direction_me_to_neighbour[2] = -1; 
                                Connectivity[nBlk].neighB_info.direction_neighbour_to_me[0] = 0; 
                                Connectivity[nBlk].neighB_info.direction_neighbour_to_me[1] = 0;
                                Connectivity[nBlk].neighB_info.direction_neighbour_to_me[2] = 1;
                                Connectivity[nBlk].be.on_grid_boundary[BE::B] = 0;
                            } /* endif */
                            break;
                            
                            case BE::T :   // 14
                            if (Connectivity[nBlk].num_neighT == 0) {
                                Connectivity[nBlk].num_neighT = 1;
                                Connectivity[nBlk].neighT = opposite_nBlk;
                                Connectivity[nBlk].neighT_info.ctm_offsets[0] = 1; 
                                Connectivity[nBlk].neighT_info.ctm_offsets[1] = 2;  
                                Connectivity[nBlk].neighT_info.ctm_offsets[2] = 3;
                                Connectivity[nBlk].neighT_info.ctm_offsets[3] = 0; 
                                Connectivity[nBlk].neighT_info.ctm_offsets[4] = 0;  
                                Connectivity[nBlk].neighT_info.ctm_offsets[5] = -Input.NCells_Kdir;
                                Connectivity[nBlk].neighT_info.direction_me_to_neighbour[0] = 0; 
                                Connectivity[nBlk].neighT_info.direction_me_to_neighbour[1] = 0;
                                Connectivity[nBlk].neighT_info.direction_me_to_neighbour[2] = 1; 
                                Connectivity[nBlk].neighT_info.direction_neighbour_to_me[0] = 0; 
                                Connectivity[nBlk].neighT_info.direction_neighbour_to_me[1] = 0;
                                Connectivity[nBlk].neighT_info.direction_neighbour_to_me[2] = -1;
                                Connectivity[nBlk].be.on_grid_boundary[BE::T] = 0;
                            } /* endif */
                            break;
                            
                            case BE::BN :  // 15
                            if (Connectivity[nBlk].num_neighBN == 0) {
                                Connectivity[nBlk].num_neighBN = 1;
                                Connectivity[nBlk].neighBN[0] = opposite_nBlk;
                                Connectivity[nBlk].neighBN_info[0].ctm_offsets[0] = 1; 
                                Connectivity[nBlk].neighBN_info[0].ctm_offsets[1] = 2;  
                                Connectivity[nBlk].neighBN_info[0].ctm_offsets[2] = 3;
                                Connectivity[nBlk].neighBN_info[0].ctm_offsets[3] = 0; 
                                Connectivity[nBlk].neighBN_info[0].ctm_offsets[4] = -Input.NCells_Jdir;  
                                Connectivity[nBlk].neighBN_info[0].ctm_offsets[5] = Input.NCells_Kdir;  
                                Connectivity[nBlk].neighBN_info[0].direction_me_to_neighbour[0] = 0; 
                                Connectivity[nBlk].neighBN_info[0].direction_me_to_neighbour[1] = 1;
                                Connectivity[nBlk].neighBN_info[0].direction_me_to_neighbour[2] = -1; 
                                Connectivity[nBlk].neighBN_info[0].direction_neighbour_to_me[0] = 0; 
                                Connectivity[nBlk].neighBN_info[0].direction_neighbour_to_me[1] = -1;
                                Connectivity[nBlk].neighBN_info[0].direction_neighbour_to_me[2] = 1;
                                Connectivity[nBlk].be.on_grid_boundary[BE::BN] = 0;
                            } /* endif */
                            break;
                            
                            case BE::N :   // 16
                            if (Connectivity[nBlk].num_neighN == 0) {
                                Connectivity[nBlk].num_neighN = 1;
                                Connectivity[nBlk].neighN = opposite_nBlk;
                                Connectivity[nBlk].neighN_info.ctm_offsets[0] = 1; 
                                Connectivity[nBlk].neighN_info.ctm_offsets[1] = 2;  
                                Connectivity[nBlk].neighN_info.ctm_offsets[2] = 3;
                                Connectivity[nBlk].neighN_info.ctm_offsets[3] = 0; 
                                Connectivity[nBlk].neighN_info.ctm_offsets[4] = -Input.NCells_Jdir;  
                                Connectivity[nBlk].neighN_info.ctm_offsets[5] = 0;
                                Connectivity[nBlk].neighN_info.direction_me_to_neighbour[0] = 0; 
                                Connectivity[nBlk].neighN_info.direction_me_to_neighbour[1] = 1;
                                Connectivity[nBlk].neighN_info.direction_me_to_neighbour[2] = 0; 
                                Connectivity[nBlk].neighN_info.direction_neighbour_to_me[0] = 0; 
                                Connectivity[nBlk].neighN_info.direction_neighbour_to_me[1] = -1;
                                Connectivity[nBlk].neighN_info.direction_neighbour_to_me[2] = 0;
                                Connectivity[nBlk].be.on_grid_boundary[BE::N] = 0;
                            } /* endif */
                            break;
                            
                            case BE::TN :  // 17
                            if (Connectivity[nBlk].num_neighTN == 0) {
                                Connectivity[nBlk].num_neighTN = 1;
                                Connectivity[nBlk].neighTN[0] = opposite_nBlk;
                                Connectivity[nBlk].neighTN_info[0].ctm_offsets[0] = 1; 
                                Connectivity[nBlk].neighTN_info[0].ctm_offsets[1] = 2;  
                                Connectivity[nBlk].neighTN_info[0].ctm_offsets[2] = 3;
                                Connectivity[nBlk].neighTN_info[0].ctm_offsets[3] = 0; 
                                Connectivity[nBlk].neighTN_info[0].ctm_offsets[4] = -Input.NCells_Jdir;  
                                Connectivity[nBlk].neighTN_info[0].ctm_offsets[5] = -Input.NCells_Kdir;  
                                Connectivity[nBlk].neighTN_info[0].direction_me_to_neighbour[0] = 0; 
                                Connectivity[nBlk].neighTN_info[0].direction_me_to_neighbour[1] = 1;
                                Connectivity[nBlk].neighTN_info[0].direction_me_to_neighbour[2] = 1; 
                                Connectivity[nBlk].neighTN_info[0].direction_neighbour_to_me[0] = 0; 
                                Connectivity[nBlk].neighTN_info[0].direction_neighbour_to_me[1] = -1;
                                Connectivity[nBlk].neighTN_info[0].direction_neighbour_to_me[2] = -1;
                                Connectivity[nBlk].be.on_grid_boundary[BE::TN] = 0;
                            } /* endif */
                            break;
                            
                            case BE::BSE : // 18
                            if (Connectivity[nBlk].num_neighBSE == 0 && 
                                (Input.i_Grid != GRID_PERIODIC_BOX_WITH_INFLOW ||
                                 (Input.i_Grid == GRID_PERIODIC_BOX_WITH_INFLOW &&
                                  iBlk < Input.NBlk_Idir-1))
                                ) {
                                Connectivity[nBlk].num_neighBSE = 1;
                                Connectivity[nBlk].neighBSE[0] = opposite_nBlk;
                                Connectivity[nBlk].neighBSE_info[0].ctm_offsets[0] = 1; 
                                Connectivity[nBlk].neighBSE_info[0].ctm_offsets[1] = 2;  
                                Connectivity[nBlk].neighBSE_info[0].ctm_offsets[2] = 3;
                                Connectivity[nBlk].neighBSE_info[0].ctm_offsets[3] = Input.NCells_Idir; 
                                Connectivity[nBlk].neighBSE_info[0].ctm_offsets[4] = -Input.NCells_Jdir;  
                                Connectivity[nBlk].neighBSE_info[0].ctm_offsets[5] = -Input.NCells_Kdir;  
                                Connectivity[nBlk].neighBSE_info[0].direction_me_to_neighbour[0] = 1; 
                                Connectivity[nBlk].neighBSE_info[0].direction_me_to_neighbour[1] = -1;
                                Connectivity[nBlk].neighBSE_info[0].direction_me_to_neighbour[2] = -1; 
                                Connectivity[nBlk].neighBSE_info[0].direction_neighbour_to_me[0] = -1; 
                                Connectivity[nBlk].neighBSE_info[0].direction_neighbour_to_me[1] = 1;
                                Connectivity[nBlk].neighBSE_info[0].direction_neighbour_to_me[2] = 1;
                                Connectivity[nBlk].be.on_grid_boundary[BE::BSE] = 0;
                            } /* endif */
                            break;
                            
                            case BE::SE :  // 19
                            if (Connectivity[nBlk].num_neighSE == 0 && 
                                (Input.i_Grid != GRID_PERIODIC_BOX_WITH_INFLOW ||
                                 (Input.i_Grid == GRID_PERIODIC_BOX_WITH_INFLOW &&
                                  iBlk < Input.NBlk_Idir-1))
                                ) {
                                Connectivity[nBlk].num_neighSE = 1;
                                Connectivity[nBlk].neighSE[0] = opposite_nBlk;
                                Connectivity[nBlk].neighSE_info[0].ctm_offsets[0] = 1; 
                                Connectivity[nBlk].neighSE_info[0].ctm_offsets[1] = 2;  
                                Connectivity[nBlk].neighSE_info[0].ctm_offsets[2] = 3;
                                Connectivity[nBlk].neighSE_info[0].ctm_offsets[3] = Input.NCells_Idir; 
                                Connectivity[nBlk].neighSE_info[0].ctm_offsets[4] = -Input.NCells_Jdir;  
                                Connectivity[nBlk].neighSE_info[0].ctm_offsets[5] = 0;  
                                Connectivity[nBlk].neighSE_info[0].direction_me_to_neighbour[0] = 1; 
                                Connectivity[nBlk].neighSE_info[0].direction_me_to_neighbour[1] = -1;
                                Connectivity[nBlk].neighSE_info[0].direction_me_to_neighbour[2] = 0; 
                                Connectivity[nBlk].neighSE_info[0].direction_neighbour_to_me[0] = -1; 
                                Connectivity[nBlk].neighSE_info[0].direction_neighbour_to_me[1] = 1;
                                Connectivity[nBlk].neighSE_info[0].direction_neighbour_to_me[2] = 0;
                                Connectivity[nBlk].be.on_grid_boundary[BE::SE] = 0;
                            } /* endif */
                            break;
                            
                            case BE::TSE : // 20
                            if (Connectivity[nBlk].num_neighTSE == 0 && 
                                (Input.i_Grid != GRID_PERIODIC_BOX_WITH_INFLOW ||
                                 (Input.i_Grid == GRID_PERIODIC_BOX_WITH_INFLOW &&
                                  iBlk < Input.NBlk_Idir-1))
                                ) {
                                Connectivity[nBlk].num_neighTSE = 1;
                                Connectivity[nBlk].neighTSE[0] = opposite_nBlk;
                                Connectivity[nBlk].neighTSE_info[0].ctm_offsets[0] = 1; 
                                Connectivity[nBlk].neighTSE_info[0].ctm_offsets[1] = 2;  
                                Connectivity[nBlk].neighTSE_info[0].ctm_offsets[2] = 3;
                                Connectivity[nBlk].neighTSE_info[0].ctm_offsets[3] = -Input.NCells_Idir; 
                                Connectivity[nBlk].neighTSE_info[0].ctm_offsets[4] = Input.NCells_Jdir;  
                                Connectivity[nBlk].neighTSE_info[0].ctm_offsets[5] = -Input.NCells_Kdir;  
                                Connectivity[nBlk].neighTSE_info[0].direction_me_to_neighbour[0] = 1; 
                                Connectivity[nBlk].neighTSE_info[0].direction_me_to_neighbour[1] = -1;
                                Connectivity[nBlk].neighTSE_info[0].direction_me_to_neighbour[2] = 1; 
                                Connectivity[nBlk].neighTSE_info[0].direction_neighbour_to_me[0] = -1; 
                                Connectivity[nBlk].neighTSE_info[0].direction_neighbour_to_me[1] = 1;
                                Connectivity[nBlk].neighTSE_info[0].direction_neighbour_to_me[2] = -1;
                                Connectivity[nBlk].be.on_grid_boundary[BE::TSE] = 0;
                            } /* endif */
                            break;
                            
                            case BE::BE :  // 21
                            if (Connectivity[nBlk].num_neighBE == 0 && 
                                (Input.i_Grid != GRID_PERIODIC_BOX_WITH_INFLOW ||
                                 (Input.i_Grid == GRID_PERIODIC_BOX_WITH_INFLOW &&
                                  iBlk < Input.NBlk_Idir-1))
                                ) {
                                Connectivity[nBlk].num_neighBE = 1;
                                Connectivity[nBlk].neighBE[0] = opposite_nBlk;
                                Connectivity[nBlk].neighBE_info[0].ctm_offsets[0] = 1; 
                                Connectivity[nBlk].neighBE_info[0].ctm_offsets[1] = 2;  
                                Connectivity[nBlk].neighBE_info[0].ctm_offsets[2] = 3;
                                Connectivity[nBlk].neighBE_info[0].ctm_offsets[3] = -Input.NCells_Idir; 
                                Connectivity[nBlk].neighBE_info[0].ctm_offsets[4] = 0;  
                                Connectivity[nBlk].neighBE_info[0].ctm_offsets[5] = Input.NCells_Kdir;  
                                Connectivity[nBlk].neighBE_info[0].direction_me_to_neighbour[0] = 1; 
                                Connectivity[nBlk].neighBE_info[0].direction_me_to_neighbour[1] = 0;
                                Connectivity[nBlk].neighBE_info[0].direction_me_to_neighbour[2] = -1; 
                                Connectivity[nBlk].neighBE_info[0].direction_neighbour_to_me[0] = -1; 
                                Connectivity[nBlk].neighBE_info[0].direction_neighbour_to_me[1] = 0;
                                Connectivity[nBlk].neighBE_info[0].direction_neighbour_to_me[2] = 1;
                                Connectivity[nBlk].be.on_grid_boundary[BE::BE] = 0;
                            } /* endif */
                            break;
                            
                            case BE::E :   // 22
                            if (Connectivity[nBlk].num_neighE == 0 && 
                                (Input.i_Grid != GRID_PERIODIC_BOX_WITH_INFLOW ||
                                 (Input.i_Grid == GRID_PERIODIC_BOX_WITH_INFLOW &&
                                  iBlk < Input.NBlk_Idir-1))
                                ) {
                                Connectivity[nBlk].num_neighE = 1;
                                Connectivity[nBlk].neighE = opposite_nBlk;
                                Connectivity[nBlk].neighE_info.ctm_offsets[0] = 1; 
                                Connectivity[nBlk].neighE_info.ctm_offsets[1] = 2;  
                                Connectivity[nBlk].neighE_info.ctm_offsets[2] = 3;
                                Connectivity[nBlk].neighE_info.ctm_offsets[3] = -Input.NCells_Idir; 
                                Connectivity[nBlk].neighE_info.ctm_offsets[4] = 0;  
                                Connectivity[nBlk].neighE_info.ctm_offsets[5] = 0;
                                Connectivity[nBlk].neighE_info.direction_me_to_neighbour[0] = 1; 
                                Connectivity[nBlk].neighE_info.direction_me_to_neighbour[1] = 0;
                                Connectivity[nBlk].neighE_info.direction_me_to_neighbour[2] = 0; 
                                Connectivity[nBlk].neighE_info.direction_neighbour_to_me[0] = -1; 
                                Connectivity[nBlk].neighE_info.direction_neighbour_to_me[1] = 0;
                                Connectivity[nBlk].neighE_info.direction_neighbour_to_me[2] = 0;
                                Connectivity[nBlk].be.on_grid_boundary[BE::E] = 0;
                            } /* endif */
                            break;
                            
                            case BE::TE :  // 23
                            if (Connectivity[nBlk].num_neighTE == 0 && 
                                (Input.i_Grid != GRID_PERIODIC_BOX_WITH_INFLOW ||
                                 (Input.i_Grid == GRID_PERIODIC_BOX_WITH_INFLOW &&
                                  iBlk < Input.NBlk_Idir-1))
                                ) {
                                Connectivity[nBlk].num_neighTE = 1;
                                Connectivity[nBlk].neighTE[0] = opposite_nBlk;
                                Connectivity[nBlk].neighTE_info[0].ctm_offsets[0] = 1; 
                                Connectivity[nBlk].neighTE_info[0].ctm_offsets[1] = 2;  
                                Connectivity[nBlk].neighTE_info[0].ctm_offsets[2] = 3;
                                Connectivity[nBlk].neighTE_info[0].ctm_offsets[3] = -Input.NCells_Idir; 
                                Connectivity[nBlk].neighTE_info[0].ctm_offsets[4] = 0;  
                                Connectivity[nBlk].neighTE_info[0].ctm_offsets[5] = -Input.NCells_Kdir;  
                                Connectivity[nBlk].neighTE_info[0].direction_me_to_neighbour[0] = 1; 
                                Connectivity[nBlk].neighTE_info[0].direction_me_to_neighbour[1] = 0;
                                Connectivity[nBlk].neighTE_info[0].direction_me_to_neighbour[2] = 1; 
                                Connectivity[nBlk].neighTE_info[0].direction_neighbour_to_me[0] = -1; 
                                Connectivity[nBlk].neighTE_info[0].direction_neighbour_to_me[1] = 0;
                                Connectivity[nBlk].neighTE_info[0].direction_neighbour_to_me[2] = -1;
                                Connectivity[nBlk].be.on_grid_boundary[BE::TE] = 0;
                            } /* endif */
                            break;
                            
                            case BE::BNE : // 24
                            if (Connectivity[nBlk].num_neighBNE == 0 && 
                                (Input.i_Grid != GRID_PERIODIC_BOX_WITH_INFLOW ||
                                 (Input.i_Grid == GRID_PERIODIC_BOX_WITH_INFLOW &&
                                  iBlk < Input.NBlk_Idir-1))
                                ) {
                                Connectivity[nBlk].num_neighBNE = 1;
                                Connectivity[nBlk].neighBNE[0] = opposite_nBlk;
                                Connectivity[nBlk].neighBNE[0] = opposite_nBlk;
                                Connectivity[nBlk].neighBNE_info[0].ctm_offsets[0] = 1; 
                                Connectivity[nBlk].neighBNE_info[0].ctm_offsets[1] = 2;  
                                Connectivity[nBlk].neighBNE_info[0].ctm_offsets[2] = 3;
                                Connectivity[nBlk].neighBNE_info[0].ctm_offsets[3] = -Input.NCells_Idir; 
                                Connectivity[nBlk].neighBNE_info[0].ctm_offsets[4] = -Input.NCells_Jdir;  
                                Connectivity[nBlk].neighBNE_info[0].ctm_offsets[5] = Input.NCells_Kdir;  
                                Connectivity[nBlk].neighBNE_info[0].direction_me_to_neighbour[0] = 1; 
                                Connectivity[nBlk].neighBNE_info[0].direction_me_to_neighbour[1] = 1;
                                Connectivity[nBlk].neighBNE_info[0].direction_me_to_neighbour[2] = -1; 
                                Connectivity[nBlk].neighBNE_info[0].direction_neighbour_to_me[0] = -1; 
                                Connectivity[nBlk].neighBNE_info[0].direction_neighbour_to_me[1] = -1;
                                Connectivity[nBlk].neighBNE_info[0].direction_neighbour_to_me[2] = 1;
                                Connectivity[nBlk].be.on_grid_boundary[BE::BNE] = 0;
                            } /* endif */
                            break;
                            
                            case BE::NE :  // 25
                            if (Connectivity[nBlk].num_neighNE == 0 && 
                                (Input.i_Grid != GRID_PERIODIC_BOX_WITH_INFLOW ||
                                 (Input.i_Grid == GRID_PERIODIC_BOX_WITH_INFLOW &&
                                  iBlk < Input.NBlk_Idir-1))
                                ) {
                                Connectivity[nBlk].num_neighNE = 1;
                                Connectivity[nBlk].neighNE[0] = opposite_nBlk;
                                Connectivity[nBlk].neighNE_info[0].ctm_offsets[0] = 1; 
                                Connectivity[nBlk].neighNE_info[0].ctm_offsets[1] = 2;  
                                Connectivity[nBlk].neighNE_info[0].ctm_offsets[2] = 3;
                                Connectivity[nBlk].neighNE_info[0].ctm_offsets[3] = -Input.NCells_Idir; 
                                Connectivity[nBlk].neighNE_info[0].ctm_offsets[4] = -Input.NCells_Jdir;  
                                Connectivity[nBlk].neighNE_info[0].ctm_offsets[5] = 0;  
                                Connectivity[nBlk].neighNE_info[0].direction_me_to_neighbour[0] = 1; 
                                Connectivity[nBlk].neighNE_info[0].direction_me_to_neighbour[1] = 1;
                                Connectivity[nBlk].neighNE_info[0].direction_me_to_neighbour[2] = 0; 
                                Connectivity[nBlk].neighNE_info[0].direction_neighbour_to_me[0] = -1; 
                                Connectivity[nBlk].neighNE_info[0].direction_neighbour_to_me[1] = -1;
                                Connectivity[nBlk].neighNE_info[0].direction_neighbour_to_me[2] = 0;
                                Connectivity[nBlk].be.on_grid_boundary[BE::NE] = 0;
                            } /* endif */
                            break;
                            
                            case BE::TNE : // 26
                            if (Connectivity[nBlk].num_neighTNE == 0 && 
                                (Input.i_Grid != GRID_PERIODIC_BOX_WITH_INFLOW ||
                                 (Input.i_Grid == GRID_PERIODIC_BOX_WITH_INFLOW &&
                                  iBlk < Input.NBlk_Idir-1))
                                ) {
                                Connectivity[nBlk].num_neighTNE = 1;
                                Connectivity[nBlk].neighTNE[0] = opposite_nBlk;
                                Connectivity[nBlk].neighTNE_info[0].ctm_offsets[0] = 1; 
                                Connectivity[nBlk].neighTNE_info[0].ctm_offsets[1] = 2;  
                                Connectivity[nBlk].neighTNE_info[0].ctm_offsets[2] = 3;
                                Connectivity[nBlk].neighTNE_info[0].ctm_offsets[3] = -Input.NCells_Idir; 
                                Connectivity[nBlk].neighTNE_info[0].ctm_offsets[4] = -Input.NCells_Jdir;  
                                Connectivity[nBlk].neighTNE_info[0].ctm_offsets[5] = -Input.NCells_Kdir;  
                                Connectivity[nBlk].neighTNE_info[0].direction_me_to_neighbour[0] = 1; 
                                Connectivity[nBlk].neighTNE_info[0].direction_me_to_neighbour[1] = 1;
                                Connectivity[nBlk].neighTNE_info[0].direction_me_to_neighbour[2] = 1; 
                                Connectivity[nBlk].neighTNE_info[0].direction_neighbour_to_me[0] = -1; 
                                Connectivity[nBlk].neighTNE_info[0].direction_neighbour_to_me[1] = -1;
                                Connectivity[nBlk].neighTNE_info[0].direction_neighbour_to_me[2] = -1;
                                Connectivity[nBlk].be.on_grid_boundary[BE::TNE] = 0;
                            } /* endif */
                            break;
                            
                            default :
                            break;
                    } /* endswitch */
                } /* endif */
                
            } /* endfor */
        } /* endfor */
    } /* endfor */
    
}

/**************************************************************
 * Routine: Create_Grid_Turbulence_Box                        *
 *                                                            *
 * Generates a 3D Cartesian multiblock mesh for a             * 
 * turbulence box.                                            *
 *                                                            *
 **************************************************************/
void Grid3D_Hexa_Multi_Block_List::Create_Grid_Turbulence_Box(Grid3D_Input_Parameters &Input) {
    
    IsAuxiliary = true;
    int count_blocks;
    int BC_east, BC_west, BC_north, BC_south, BC_top, BC_bottom;
    Grid2D_Quad_Block **Grid2D_Box_XYplane;
    
    
    
    int NBlk_Idir, NBlk_Jdir, NBlk_Kdir;
    NBlk_Idir = Input.NBlk_Idir;
    NBlk_Jdir = Input.NBlk_Jdir;
    NBlk_Kdir = Input.NBlk_Kdir;
    
    Input.NBlk_Idir = 1;
    Input.NBlk_Jdir = 1;
    Input.NBlk_Kdir = 1;
    
    /* Allocate required memory. */
    
    Allocate(Input.NBlk_Idir, Input.NBlk_Jdir, Input.NBlk_Kdir);
    
    /* Creat 2D cross-section grids from which the 3D grid
     will be extruded. */
    
    Grid2D_Box_XYplane = Grid_Rectangular_Box(Grid2D_Box_XYplane,
                                              Input.NBlk_Idir, 
                                              Input.NBlk_Jdir,
                                              Input.Turbulence_Box_Width,
                                              Input.Turbulence_Box_Height,
                                              OFF,  // no mesh stretching!
                                              Input.Stretching_Type_Idir,
                                              Input.Stretching_Type_Jdir,
                                              Input.Stretching_Factor_Idir,
                                              Input.Stretching_Factor_Jdir,
                                              Input.NCells_Turbulence_Idir,
                                              Input.NCells_Turbulence_Jdir,
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
                                                Input.NCells_Turbulence_Kdir,
                                                Input.Stretching_Type_Kdir,
                                                Input.Stretching_Factor_Kdir,
                                                ZERO*Input.Turbulence_Box_Length+  
                                                (double(kBlk)/double(Input.NBlk_Kdir))*Input.Turbulence_Box_Length,
                                                ZERO*Input.Turbulence_Box_Length+  
                                                (double(kBlk+1)/double(Input.NBlk_Kdir))*Input.Turbulence_Box_Length);
                
                
                if (iBlk == Input.NBlk_Idir-1) {
                    BC_east = BC_REFLECTION;
                } else {
                    BC_east = BC_NONE;
                } /* endif */
                if (iBlk == 0) {
                    BC_west = BC_REFLECTION;
                } else {
                    BC_west = BC_NONE;
                } /* endif */
                
                if (jBlk == Input.NBlk_Jdir-1) {
                    BC_north = BC_REFLECTION;
                } else {
                    BC_north = BC_NONE;
                } /* endif */
                if (jBlk == 0) {
                    BC_south = BC_REFLECTION;
                } else {
                    BC_south = BC_NONE;
                } /* endif */
                
                /* Assign top and bottom boundary conditions. */
                
                if (kBlk == Input.NBlk_Kdir-1) {
                    BC_top = BC_FIXED_PRESSURE;
                } else {
                    BC_top = BC_NONE;
                } /* endif */
                if (kBlk == 0) {
                    BC_bottom = BC_DIRICHLET;
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
                
                count_blocks++;
                
            } /* endfor */
        } /* endfor */
    } /* endfor */
    
    /* Deallocate 2D grid. */
    
    Grid2D_Box_XYplane = Deallocate_Multi_Block_Grid(Grid2D_Box_XYplane,
                                                     Input.NBlk_Idir, 
                                                     Input.NBlk_Jdir);
    
    
    /* Call the function Find_Neighbours to obtain the neighbour block information
     and assign values to data members in the grid block connectivity data structure. */
    
    //Find_Neighbours(Input);
    Input.NBlk_Idir = NBlk_Idir;
    Input.NBlk_Jdir = NBlk_Jdir;
    Input.NBlk_Kdir = NBlk_Kdir;
    
}


void Grid3D_Hexa_Multi_Block_List::Create_Uniform_Initial_Grid(Grid3D_Input_Parameters &Input, Grid3D_Hexa_Multi_Block_List &Initial_Mesh) {
    /*
     * Now only deals with a carthesian box-like configuration, 
     * for use in creating an initial turbulence field on a non-uniform grid
     *
     */
    
    IsAuxiliary = true;
    
    /* ----------------------- save old values --------------------------- */
    int NBlk_Idir = Input.NBlk_Idir;
    int NBlk_Jdir = Input.NBlk_Jdir;
    int NBlk_Kdir = Input.NBlk_Kdir;
    
    int NCells_Idir = Input.NCells_Idir;
    int NCells_Jdir = Input.NCells_Jdir;
    int NCells_Kdir = Input.NCells_Kdir;
    
    int Mesh_Stretching = Input.Mesh_Stretching;

    
    /* --------------------- Generage new Input values -------------------- */
    // Set dimensions same as Initial_Mesh
    Input.Turbulence_Box_Width  = Input.Box_Width;
    Input.Turbulence_Box_Height  = Input.Box_Height;
    Input.Turbulence_Box_Length = Input.Box_Length;
    
    // Use only one block since the turbulence will be generated on one processor anyway
    Input.NBlk_Idir = 1;
    Input.NBlk_Jdir = 1;
    Input.NBlk_Kdir = 1;
    
    // Turn off Mesh stretching
    Input.Mesh_Stretching = OFF;
    
    // Find smallest spacing
    Vector3D Delta_min = Initial_Mesh.Delta_minimum();
    cout << "Delta_min = " << Delta_min << endl;
    
    // Set number of cells in each direction
    Input.NCells_Turbulence_Idir = int(ceil(Input.Turbulence_Box_Width/Delta_min.x));
    Input.NCells_Turbulence_Jdir = int(ceil(Input.Turbulence_Box_Height/Delta_min.y));
    Input.NCells_Turbulence_Kdir = int(ceil(Input.Turbulence_Box_Length/Delta_min.z));

    Input.NCells_Idir = Input.NCells_Turbulence_Idir;
    Input.NCells_Jdir = Input.NCells_Turbulence_Jdir;
    Input.NCells_Kdir = Input.NCells_Turbulence_Kdir;

    cout << endl;
    cout << "NCells_Turbulence_Idir = " << Input.NCells_Turbulence_Idir << endl;
    cout << "NCells_Turbulence_Jdir = " << Input.NCells_Turbulence_Jdir << endl;
    cout << "NCells_Turbulence_Kdir = " << Input.NCells_Turbulence_Kdir << endl;

    /* ------------------------- Create the new grid ------------------------- */
    // Create the uniform grid with same dimensions as original and smallest spacing
    Create_Grid(Input);
    
    
    
    
    /* --------------------------- restore old values -------------------------- */
    Input.NBlk_Idir = NBlk_Idir;
    Input.NBlk_Jdir = NBlk_Jdir;
    Input.NBlk_Kdir = NBlk_Kdir;
    
    Input.NCells_Idir = NCells_Idir;
    Input.NCells_Jdir = NCells_Jdir;
    Input.NCells_Kdir = NCells_Kdir;

    Input.Mesh_Stretching = Mesh_Stretching;
}


/********************************************************
 * Routine: Create_Grid_Bunsen_Inflow                   *
 *                                                      *
 * Generates a 3D Cartesian multiblock mesh for a       *
 * Bunsen inflow turbulence field.                      *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block_List::Create_Grid_Bunsen_Inflow(Grid3D_Input_Parameters &Input) {
    
    int nBlk, numblk_idir_fuel, numblk_jdir_fuel;
    int NBlk_Idir, NBlk_Jdir, NBlk_Kdir;
    
    Grid2D_Quad_Block **Grid2D_Fuel_Line_XYplane;
    
    NBlk_Idir = 5;
    NBlk_Jdir = 1;
    NBlk_Kdir = 1;
    
    /* Allocate required memory. */
    
    Allocate(NBlk_Idir, NBlk_Jdir, NBlk_Kdir);
    
    /* Creat 2D cross-section grids from which the 3D grid
     will be extruded. */
    
    Grid2D_Fuel_Line_XYplane = Grid_Tube_2D(Grid2D_Fuel_Line_XYplane,
                                            numblk_idir_fuel,
                                            numblk_jdir_fuel,
                                            Input.Radius_Bunsen_Burner_Fuel_Line,
                                            Input.NCells_Idir,
                                            Input.NCells_Jdir,
                                            Input.Nghost,
                                            STRETCHING_FCN_MAX_CLUSTERING,
                                            1.25);
    
    /* Create the mesh for each block representing
     the complete grid. */
    
    for (int iBlk = 0; iBlk <= NBlk_Idir-1; ++iBlk) {
        
        /* Extrude each of the grid blocks from the
         appropriate 2D grid in XY-plane. */
        
        Grid_Blks[iBlk].Extrude(Grid2D_Fuel_Line_XYplane[iBlk][0],
                                Input.NCells_Kdir,
                                STRETCHING_FCN_MIN_CLUSTERING,
                                1.25,
                                ZERO*Input.Height_Bunsen_Burner,
                                0.5*Input.Height_Bunsen_Burner);
        
        Grid_Blks[iBlk].Set_BCs(BC_NONE,
                                BC_NONE,
                                BC_NONE,
                                BC_NONE,
                                BC_NONE,
                                BC_NONE);
        
    } /* endfor */
    
    /* Deallocate 2D grid. */
    
    Grid2D_Fuel_Line_XYplane = Deallocate_Multi_Block_Grid(Grid2D_Fuel_Line_XYplane,
                                                           numblk_idir_fuel,
                                                           numblk_jdir_fuel);
    
    /* Call the function Find_Neighbours to obtain the neighbour block information
     and assign values to data members in the grid block connectivity data structure. */
    
    Find_Neighbours(Input);
    
}

/**************************************************************
 * Routine: Create_Grid_Bunsen_Box                            *
 *                                                            *
 * Generates a 3D Cartesian multiblock mesh for a Bunsen box. *
 *                                                            *
 **************************************************************/
void Grid3D_Hexa_Multi_Block_List::Create_Grid_Bunsen_Box(Grid3D_Input_Parameters &Input) {
    
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
                                              OFF,  // ON 
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
                                                ZERO*Input.Box_Length+
                                                (double(kBlk)/double(Input.NBlk_Kdir))*Input.Box_Length,
                                                ZERO*Input.Box_Length+
                                                (double(kBlk+1)/double(Input.NBlk_Kdir))*Input.Box_Length);
                
                
                if (iBlk == Input.NBlk_Idir-1) {
                    BC_east = BC_OUTFLOW_SUBSONIC; 
                } else {
                    BC_east = BC_NONE;
                } /* endif */
                if (iBlk == 0) {
                    BC_west = BC_OUTFLOW_SUBSONIC;
                } else {
                    BC_west = BC_NONE;
                } /* endif */
                
                if (jBlk == Input.NBlk_Jdir-1) {
                    BC_north = BC_OUTFLOW_SUBSONIC; 
                } else {
                    BC_north = BC_NONE;
                } /* endif */
                if (jBlk == 0) {
                    BC_south = BC_OUTFLOW_SUBSONIC; 
                } else {
                    BC_south = BC_NONE;
                } /* endif */
                
                /* Assign top and bottom boundary conditions. */
                
                if (kBlk == Input.NBlk_Kdir-1) {
                    BC_top = BC_OUTFLOW_SUBSONIC; 
                } else {
                    BC_top = BC_NONE;
                } /* endif */
                if (kBlk == 0) {
                    if (jBlk == 2  ||  jBlk == 3) {  
                        BC_bottom = BC_INFLOW_TURBULENCE; 
                    } else {
                        BC_bottom = BC_INFLOW_SUBSONIC;  
                    }
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
                
                count_blocks++;
                
            } /* endfor */
        } /* endfor */
    } /* endfor */
    
    /* Deallocate 2D grid. */
    
    Grid2D_Box_XYplane = Deallocate_Multi_Block_Grid(Grid2D_Box_XYplane,
                                                     Input.NBlk_Idir, 
                                                     Input.NBlk_Jdir);
    
    
    /* Call the function Find_Neighbours to obtain the neighbour block information
     and assign values to data members in the grid block connectivity data structure. */
    
    Find_Neighbours(Input);
    
}

/********************************************************
 * Routine: Create_Grid_Bunsen_Burner                   *
 *                                                      *
 * Generates a 3D Cartesian multiblock mesh for a       *
 * Bunsen burner.                                       *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block_List::Create_Grid_Bunsen_Burner(Grid3D_Input_Parameters &Input) {
    
    int BC_east, BC_west, BC_north, BC_south, BC_top, BC_bottom;
    int nBlk, numblk_idir_fuel, numblk_jdir_fuel,
    numblk_idir_bunsenburner, numblk_jdir_bunsenburner,
    numblk_idir_coflow, numblk_jdir_coflow;
    Grid2D_Quad_Block **Grid2D_Fuel_Line_XYplane,
    **Grid2D_Bunsen_Burner_Inner_XYplane,
    **Grid2D_Bunsen_Burner_Middle_XYplane,
    **Grid2D_Bunsen_Burner_Outer_XYplane;
    
    /* Allocate required memory. */
    
    Allocate(Input.NBlk_Idir, Input.NBlk_Jdir, Input.NBlk_Kdir);
    
    /* Creat 2D cross-section grids from which the 3D grid
     will be extruded. */
    
    Grid2D_Fuel_Line_XYplane = Grid_Tube_2D(Grid2D_Fuel_Line_XYplane,
                                            numblk_idir_fuel,
                                            numblk_jdir_fuel,
                                            Input.Radius_Bunsen_Burner_Fuel_Line,
                                            Input.NCells_Idir,
                                            Input.NCells_Jdir,
                                            Input.Nghost,
                                            STRETCHING_FCN_LINEAR,//MAX_CLUSTERING,
                                            ZERO);//1.25);
    
    Grid2D_Bunsen_Burner_Inner_XYplane = Grid_Annulus_2D(Grid2D_Bunsen_Burner_Inner_XYplane,
							 numblk_idir_bunsenburner,
							 numblk_jdir_bunsenburner,
							 Input.Radius_Bunsen_Burner_Fuel_Line,
							 Input.Radius_Bunsen_Burner,
// 							 HALF*(Input.Radius_Bunsen_Burner+Input.Radius_Bunsen_Burner_Fuel_Line),
							 Input.NCells_Idir,
							 Input.NCells_Jdir,
							 Input.Nghost,
							 STRETCHING_FCN_LINEAR,//MIN_CLUSTERING,
							 ZERO);//1.10);

//     Grid2D_Bunsen_Burner_Middle_XYplane = Grid_Annulus_2D(Grid2D_Bunsen_Burner_Middle_XYplane,
// 							 numblk_idir_bunsenburner,
// 							 numblk_jdir_bunsenburner,
// 							 0.1*Input.Radius_Bunsen_Burner+Input.Radius_Bunsen_Burner_Fuel_Line,
// 							 0.45*Input.Radius_Bunsen_Burner+Input.Radius_Bunsen_Burner_Fuel_Line,
// 							 HALF*Input.Radius_Bunsen_Burner,
// 							 Input.NCells_Idir,
// 							 Input.NCells_Jdir,
// 							 Input.Nghost,
// 							 STRETCHING_FCN_MIN_CLUSTERING,
// 							 1.10);

//     Grid2D_Bunsen_Burner_Outer_XYplane = Grid_Annulus_2D(Grid2D_Bunsen_Burner_Outer_XYplane,
// 							 numblk_idir_bunsenburner,
// 							 numblk_jdir_bunsenburner,
// 							 HALF*(Input.Radius_Bunsen_Burner+Input.Radius_Bunsen_Burner_Fuel_Line),
// 							 Input.Radius_Bunsen_Burner,
// 							 Input.NCells_Idir,
// 							 Input.NCells_Jdir,
// 							 Input.Nghost,
// 							 STRETCHING_FCN_LINEAR,//MAX_CLUSTERING,
// 							 ZERO);//1.10);

    /* Create the mesh for each block representing
     the complete grid. */
    
    for (int iBlk = 0; iBlk <= Input.NBlk_Idir-1; ++iBlk) {
        
        /* Extrude each of the grid blocks from the
         appropriate 2D grid in XY-plane. */
        
        /* fuel lines */

        // iBlk = 72;
        if (iBlk <= 4) {
             Grid_Blks[iBlk].Extrude(Grid2D_Fuel_Line_XYplane[iBlk][0],
				     Input.NCells_Kdir,
				     STRETCHING_FCN_LINEAR,
				     ZERO,
				     ZERO*Input.Height_Bunsen_Burner,
				     0.125*Input.Height_Bunsen_Burner);

	      Grid_Blks[iBlk].Set_BCs(BC_NONE,
		                      BC_NONE,
                                      BC_NONE,
		                      BC_NONE,
                                      BC_NONE,
                                      BC_INFLOW_TURBULENCE);

           } else if (iBlk >= 5 && iBlk <= 9) {
             Grid_Blks[iBlk].Extrude(Grid2D_Fuel_Line_XYplane[iBlk-5][0],
				     Input.NCells_Kdir,
				     STRETCHING_FCN_LINEAR,
				     ZERO,
				     0.125*Input.Height_Bunsen_Burner,
				     0.25*Input.Height_Bunsen_Burner);

               Grid_Blks[iBlk].Set_BCs(BC_NONE,
				       BC_NONE,
				       BC_NONE,
				       BC_NONE,
				       BC_NONE,
				       BC_NONE);

           } else if (iBlk >= 10 && iBlk <= 14) {
             Grid_Blks[iBlk].Extrude(Grid2D_Fuel_Line_XYplane[iBlk-10][0],
				     Input.NCells_Kdir,
				     STRETCHING_FCN_LINEAR,
				     ZERO,
				     0.25*Input.Height_Bunsen_Burner,
				     0.375*Input.Height_Bunsen_Burner);

               Grid_Blks[iBlk].Set_BCs(BC_NONE,
				       BC_NONE,
				       BC_NONE,
				       BC_NONE,
				       BC_NONE,
				       BC_NONE);

           } else if (iBlk >= 15 && iBlk <= 19) {
             Grid_Blks[iBlk].Extrude(Grid2D_Fuel_Line_XYplane[iBlk-15][0],
				     Input.NCells_Kdir,
				     STRETCHING_FCN_LINEAR,
				     ZERO,
				     0.375*Input.Height_Bunsen_Burner,
				     0.5*Input.Height_Bunsen_Burner);

               Grid_Blks[iBlk].Set_BCs(BC_NONE,
				       BC_NONE,
				       BC_NONE,
				       BC_NONE,
				       BC_NONE,
				       BC_NONE);

           } else if (iBlk >= 20 && iBlk <= 24) {
             Grid_Blks[iBlk].Extrude(Grid2D_Fuel_Line_XYplane[iBlk-20][0],
				     Input.NCells_Kdir,
				     STRETCHING_FCN_LINEAR,
				     ZERO,
				     0.5*Input.Height_Bunsen_Burner,
				     0.625*Input.Height_Bunsen_Burner);

               Grid_Blks[iBlk].Set_BCs(BC_NONE,
				       BC_NONE,
				       BC_NONE,
				       BC_NONE,
				       BC_NONE,
				       BC_NONE);

           } else if (iBlk >= 25 && iBlk <= 29) {
             Grid_Blks[iBlk].Extrude(Grid2D_Fuel_Line_XYplane[iBlk-25][0],
				     Input.NCells_Kdir,
				     STRETCHING_FCN_LINEAR,
				     ZERO,
				     0.625*Input.Height_Bunsen_Burner,
				     0.75*Input.Height_Bunsen_Burner);

               Grid_Blks[iBlk].Set_BCs(BC_NONE,
				       BC_NONE,
				       BC_NONE,
				       BC_NONE,
				       BC_NONE,
				       BC_NONE);

           } else if (iBlk >= 30 && iBlk <= 34) {
             Grid_Blks[iBlk].Extrude(Grid2D_Fuel_Line_XYplane[iBlk-30][0],
				     Input.NCells_Kdir,
				     STRETCHING_FCN_LINEAR,
				     ZERO,
				     0.75*Input.Height_Bunsen_Burner,
				     0.875*Input.Height_Bunsen_Burner);

               Grid_Blks[iBlk].Set_BCs(BC_NONE,
				       BC_NONE,
				       BC_NONE,
				       BC_NONE,
				       BC_NONE,
				       BC_NONE);

           } else if (iBlk >= 35 && iBlk <= 39) {
             Grid_Blks[iBlk].Extrude(Grid2D_Fuel_Line_XYplane[iBlk-35][0],
				     Input.NCells_Kdir,
				     STRETCHING_FCN_LINEAR,
				     ZERO,
				     0.875*Input.Height_Bunsen_Burner,
				     Input.Height_Bunsen_Burner);

               Grid_Blks[iBlk].Set_BCs(BC_NONE,
				       BC_NONE,
				       BC_NONE,
				       BC_NONE,
				       BC_OUTFLOW_SUBSONIC,
				       BC_NONE);


	      /* air line */

           } else if (iBlk >= 40 && iBlk <= 43) {
              Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Inner_XYplane[iBlk-40][0],
				      Input.NCells_Kdir,
				      STRETCHING_FCN_LINEAR,
				      ZERO,
				      ZERO*Input.Height_Bunsen_Burner,
				      0.125*Input.Height_Bunsen_Burner);

               Grid_Blks[iBlk].Set_BCs(BC_NONE,
				       BC_NONE,
				       BC_OUTFLOW_SUBSONIC,
				       BC_NONE,
				       BC_NONE,
				       BC_INFLOW_SUBSONIC);

           } else if (iBlk >= 44 && iBlk <= 47) {
              Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Inner_XYplane[iBlk-44][0],
				      Input.NCells_Kdir,
				      STRETCHING_FCN_LINEAR,
				      ZERO,
				      0.125*Input.Height_Bunsen_Burner,
				      0.25*Input.Height_Bunsen_Burner);

               Grid_Blks[iBlk].Set_BCs(BC_NONE,
				       BC_NONE,
				       BC_OUTFLOW_SUBSONIC,
				       BC_NONE,
				       BC_NONE,
				       BC_NONE);

           } else if (iBlk >= 48 && iBlk <= 51) {
              Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Inner_XYplane[iBlk-48][0],
				      Input.NCells_Kdir,
				      STRETCHING_FCN_LINEAR,
				      ZERO,
				      0.25*Input.Height_Bunsen_Burner,
				      0.375*Input.Height_Bunsen_Burner);

               Grid_Blks[iBlk].Set_BCs(BC_NONE,
				       BC_NONE,
				       BC_OUTFLOW_SUBSONIC,
				       BC_NONE,
				       BC_NONE,
				       BC_NONE);

           } else if (iBlk >= 52 && iBlk <= 55) {
              Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Inner_XYplane[iBlk-52][0],
				      Input.NCells_Kdir,
				      STRETCHING_FCN_LINEAR,
				      ZERO,
				      0.375*Input.Height_Bunsen_Burner,
				      0.5*Input.Height_Bunsen_Burner);

               Grid_Blks[iBlk].Set_BCs(BC_NONE,
				       BC_NONE,
				       BC_OUTFLOW_SUBSONIC,
				       BC_NONE,
				       BC_NONE,
				       BC_NONE);

           } else if (iBlk >= 56 && iBlk <= 59) {
              Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Inner_XYplane[iBlk-56][0],
				      Input.NCells_Kdir,
				      STRETCHING_FCN_LINEAR,
				      ZERO,
				      0.5*Input.Height_Bunsen_Burner,
				      0.625*Input.Height_Bunsen_Burner);

               Grid_Blks[iBlk].Set_BCs(BC_NONE,
				       BC_NONE,
				       BC_OUTFLOW_SUBSONIC,
				       BC_NONE,
				       BC_NONE,
				       BC_NONE);

           } else if (iBlk >= 60 && iBlk <= 63) {
              Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Inner_XYplane[iBlk-60][0],
				      Input.NCells_Kdir,
				      STRETCHING_FCN_LINEAR,
				      ZERO,
				      0.625*Input.Height_Bunsen_Burner,
				      0.75*Input.Height_Bunsen_Burner);

               Grid_Blks[iBlk].Set_BCs(BC_NONE,
				       BC_NONE,
				       BC_OUTFLOW_SUBSONIC,
				       BC_NONE,
				       BC_NONE,
				       BC_NONE);

           } else if (iBlk >= 64 && iBlk <= 67) {
              Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Inner_XYplane[iBlk-64][0],
				      Input.NCells_Kdir,
				      STRETCHING_FCN_LINEAR,
				      ZERO,
				      0.75*Input.Height_Bunsen_Burner,
				      0.875*Input.Height_Bunsen_Burner);

               Grid_Blks[iBlk].Set_BCs(BC_NONE,
				       BC_NONE,
				       BC_OUTFLOW_SUBSONIC,
				       BC_NONE,
				       BC_NONE,
				       BC_NONE);

           } else if (iBlk >= 68 && iBlk <= 71) {
              Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Inner_XYplane[iBlk-68][0],
				      Input.NCells_Kdir,
				      STRETCHING_FCN_LINEAR,
				      ZERO,
				      0.875*Input.Height_Bunsen_Burner,
				      Input.Height_Bunsen_Burner);

               Grid_Blks[iBlk].Set_BCs(BC_NONE,
				       BC_NONE,
				       BC_OUTFLOW_SUBSONIC,
				       BC_NONE,
				       BC_INFLOW_SUBSONIC,
				       BC_NONE);


            // iBlk = 104;
	    /* fuel lines */
//           if (iBlk <= 4) {
//              Grid_Blks[iBlk].Extrude(Grid2D_Fuel_Line_XYplane[iBlk][0],
// 				     Input.NCells_Kdir,
// 				     STRETCHING_FCN_LINEAR,//_FCN_MIN_CLUSTERING,
// 				     ZERO,//1.10,//1.25,
// 				     ZERO*Input.Height_Bunsen_Burner,
// 				     0.125*Input.Height_Bunsen_Burner);

// 	      Grid_Blks[iBlk].Set_BCs(BC_NONE,
// 		                      BC_NONE,
//                                       BC_NONE,
// 		                      BC_NONE,
//                                       BC_NONE,
//                                       BC_INFLOW_TURBULENCE);

//            } else if (iBlk >= 5 && iBlk <= 9) {
//              Grid_Blks[iBlk].Extrude(Grid2D_Fuel_Line_XYplane[iBlk-5][0],
// 				     Input.NCells_Kdir,
// 				     STRETCHING_FCN_LINEAR,//_FCN_MIN_CLUSTERING,
// 				     ZERO,//1.10,//1.25,
// 				     0.125*Input.Height_Bunsen_Burner,
// 				     0.25*Input.Height_Bunsen_Burner);

//                Grid_Blks[iBlk].Set_BCs(BC_NONE,
// 				       BC_NONE,
// 				       BC_NONE,
// 				       BC_NONE,
// 				       BC_NONE,
// 				       BC_NONE);

//            } else if (iBlk >= 10 && iBlk <= 14) {
//              Grid_Blks[iBlk].Extrude(Grid2D_Fuel_Line_XYplane[iBlk-10][0],
// 				     Input.NCells_Kdir,
// 				     STRETCHING_FCN_LINEAR,//_FCN_MIN_CLUSTERING,
// 				     ZERO,//1.10,//1.25,
// 				     0.25*Input.Height_Bunsen_Burner,
// 				     0.375*Input.Height_Bunsen_Burner);

//                Grid_Blks[iBlk].Set_BCs(BC_NONE,
// 				       BC_NONE,
// 				       BC_NONE,
// 				       BC_NONE,
// 				       BC_NONE,
// 				       BC_NONE);

//            } else if (iBlk >= 15 && iBlk <= 19) {
//              Grid_Blks[iBlk].Extrude(Grid2D_Fuel_Line_XYplane[iBlk-15][0],
// 				     Input.NCells_Kdir,
// 				     STRETCHING_FCN_LINEAR,//_FCN_MIN_CLUSTERING,
// 				     ZERO,//1.10,//1.25,
// 				     0.375*Input.Height_Bunsen_Burner,
// 				     0.5*Input.Height_Bunsen_Burner);

//                Grid_Blks[iBlk].Set_BCs(BC_NONE,
// 				       BC_NONE,
// 				       BC_NONE,
// 				       BC_NONE,
// 				       BC_NONE,
// 				       BC_NONE);

//            } else if (iBlk >= 20 && iBlk <= 24) {
//              Grid_Blks[iBlk].Extrude(Grid2D_Fuel_Line_XYplane[iBlk-20][0],
// 				     Input.NCells_Kdir,
// 				     STRETCHING_FCN_LINEAR,//_FCN_MIN_CLUSTERING,
// 				     ZERO,//1.10,//1.25,
// 				     0.5*Input.Height_Bunsen_Burner,
// 				     0.625*Input.Height_Bunsen_Burner);

//                Grid_Blks[iBlk].Set_BCs(BC_NONE,
// 				       BC_NONE,
// 				       BC_NONE,
// 				       BC_NONE,
// 				       BC_NONE,
// 				       BC_NONE);

//            } else if (iBlk >= 25 && iBlk <= 29) {
//              Grid_Blks[iBlk].Extrude(Grid2D_Fuel_Line_XYplane[iBlk-25][0],
// 				     Input.NCells_Kdir,
// 				     STRETCHING_FCN_LINEAR,//_FCN_MIN_CLUSTERING,
// 				     ZERO,//1.10,//1.25,
// 				     0.625*Input.Height_Bunsen_Burner,
// 				     0.75*Input.Height_Bunsen_Burner);

//                Grid_Blks[iBlk].Set_BCs(BC_NONE,
// 				       BC_NONE,
// 				       BC_NONE,
// 				       BC_NONE,
// 				       BC_NONE,
// 				       BC_NONE);

//            } else if (iBlk >= 30 && iBlk <= 34) {
//              Grid_Blks[iBlk].Extrude(Grid2D_Fuel_Line_XYplane[iBlk-30][0],
// 				     Input.NCells_Kdir,
// 				     STRETCHING_FCN_LINEAR,//_FCN_MIN_CLUSTERING,
// 				     ZERO,//1.10,//1.25,
// 				     0.75*Input.Height_Bunsen_Burner,
// 				     0.875*Input.Height_Bunsen_Burner);

//                Grid_Blks[iBlk].Set_BCs(BC_NONE,
// 				       BC_NONE,
// 				       BC_NONE,
// 				       BC_NONE,
// 				       BC_NONE,
// 				       BC_NONE);

//            } else if (iBlk >= 35 && iBlk <= 39) {
//              Grid_Blks[iBlk].Extrude(Grid2D_Fuel_Line_XYplane[iBlk-35][0],
// 				     Input.NCells_Kdir,
// 				     STRETCHING_FCN_LINEAR,//_FCN_MIN_CLUSTERING,
// 				     ZERO,//1.10,//1.25,
// 				     0.875*Input.Height_Bunsen_Burner,
// 				     Input.Height_Bunsen_Burner);

//                Grid_Blks[iBlk].Set_BCs(BC_NONE,
// 				       BC_NONE,
// 				       BC_NONE,
// 				       BC_NONE,
// 				       BC_OUTFLOW_SUBSONIC,
// 				       BC_NONE);

// 	      /* first ring of air line */

//            } else if (iBlk >= 40 && iBlk <= 43) {
//               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Inner_XYplane[iBlk-40][0],
// 				      Input.NCells_Kdir,
// 				      STRETCHING_FCN_LINEAR,//_FCN_MIN_CLUSTERING,
// 				      ZERO,//1.10,//1.25,
// 				      ZERO*Input.Height_Bunsen_Burner,
// 				      0.125*Input.Height_Bunsen_Burner);

//                Grid_Blks[iBlk].Set_BCs(BC_NONE,
// 				       BC_NONE,
// 				       BC_OUTFLOW_SUBSONIC,
// 				       BC_NONE,
// 				       BC_NONE,
// 				       BC_INFLOW_SUBSONIC);

//            } else if (iBlk >= 44 && iBlk <= 47) {
//               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Inner_XYplane[iBlk-44][0],
// 				      Input.NCells_Kdir,
// 				      STRETCHING_FCN_LINEAR,
// 				      ZERO,//1.10,//1.25,
// 				      0.125*Input.Height_Bunsen_Burner,
// 				      0.25*Input.Height_Bunsen_Burner);

//                Grid_Blks[iBlk].Set_BCs(BC_NONE,
// 				       BC_NONE,
// 				       BC_OUTFLOW_SUBSONIC,
// 				       BC_NONE,
// 				       BC_NONE,
// 				       BC_NONE);

//            } else if (iBlk >= 48 && iBlk <= 51) {
//               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Inner_XYplane[iBlk-48][0],
// 				      Input.NCells_Kdir,
// 				      STRETCHING_FCN_LINEAR,
// 				      ZERO,//1.10,//1.25,
// 				      0.25*Input.Height_Bunsen_Burner,
// 				      0.375*Input.Height_Bunsen_Burner);

//                Grid_Blks[iBlk].Set_BCs(BC_NONE,
// 				       BC_NONE,
// 				       BC_OUTFLOW_SUBSONIC,
// 				       BC_NONE,
// 				       BC_NONE,
// 				       BC_NONE);

//            } else if (iBlk >= 52 && iBlk <= 55) {
//               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Inner_XYplane[iBlk-52][0],
// 				      Input.NCells_Kdir,
// 				      STRETCHING_FCN_LINEAR,
// 				      ZERO,//1.10,//1.25,
// 				      0.375*Input.Height_Bunsen_Burner,
// 				      0.5*Input.Height_Bunsen_Burner);

//                Grid_Blks[iBlk].Set_BCs(BC_NONE,
// 				       BC_NONE,
// 				       BC_OUTFLOW_SUBSONIC,
// 				       BC_NONE,
// 				       BC_NONE,
// 				       BC_NONE);

//            } else if (iBlk >= 56 && iBlk <= 59) {
//               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Inner_XYplane[iBlk-56][0],
// 				      Input.NCells_Kdir,
// 				      STRETCHING_FCN_LINEAR,
// 				      ZERO,//1.10,//1.25,
// 				      0.5*Input.Height_Bunsen_Burner,
// 				      0.625*Input.Height_Bunsen_Burner);

//                Grid_Blks[iBlk].Set_BCs(BC_NONE,
// 				       BC_NONE,
// 				       BC_OUTFLOW_SUBSONIC,
// 				       BC_NONE,
// 				       BC_NONE,
// 				       BC_NONE);

//            } else if (iBlk >= 60 && iBlk <= 63) {
//               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Inner_XYplane[iBlk-60][0],
// 				      Input.NCells_Kdir,
// 				      STRETCHING_FCN_LINEAR,
// 				      ZERO,//1.10,//1.25,
// 				      0.625*Input.Height_Bunsen_Burner,
// 				      0.75*Input.Height_Bunsen_Burner);

//                Grid_Blks[iBlk].Set_BCs(BC_NONE,
// 				       BC_NONE,
// 				       BC_OUTFLOW_SUBSONIC,
// 				       BC_NONE,
// 				       BC_NONE,
// 				       BC_NONE);

//            } else if (iBlk >= 64 && iBlk <= 67) {
//               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Inner_XYplane[iBlk-64][0],
// 				      Input.NCells_Kdir,
// 				      STRETCHING_FCN_LINEAR,
// 				      ZERO,//1.10,//1.25,
// 				      0.75*Input.Height_Bunsen_Burner,
// 				      0.875*Input.Height_Bunsen_Burner);

//                Grid_Blks[iBlk].Set_BCs(BC_NONE,
// 				       BC_NONE,
// 				       BC_OUTFLOW_SUBSONIC,
// 				       BC_NONE,
// 				       BC_NONE,
// 				       BC_NONE);

//            } else if (iBlk >= 68 && iBlk <= 71) {
//               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Inner_XYplane[iBlk-68][0],
// 				      Input.NCells_Kdir,
// 				      STRETCHING_FCN_LINEAR,
// 				      ZERO,//1.10,//1.25,
// 				      0.875*Input.Height_Bunsen_Burner,
// 				      Input.Height_Bunsen_Burner);

//                Grid_Blks[iBlk].Set_BCs(BC_NONE,
// 				       BC_NONE,
// 				       BC_OUTFLOW_SUBSONIC,
// 				       BC_NONE,
// 				       BC_OUTFLOW_SUBSONIC,
// 				       BC_NONE);

// 	      /* second ring of air line */

//            } else if (iBlk >= 72 && iBlk <= 75) {
//               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Outer_XYplane[iBlk-72][0],
// 				      Input.NCells_Kdir,
// 				      STRETCHING_FCN_LINEAR,//_FCN_MIN_CLUSTERING,
// 				      ZERO,//1.10,//1.25,
// 				      ZERO*Input.Height_Bunsen_Burner,
// 				      0.125*Input.Height_Bunsen_Burner);

//                Grid_Blks[iBlk].Set_BCs(BC_NONE,
// 				       BC_NONE,
// 				       BC_OUTFLOW_SUBSONIC,
// 				       BC_NONE,
// 				       BC_NONE,
// 				       BC_INFLOW_SUBSONIC);

//            } else if (iBlk >= 76 && iBlk <= 79) {
//               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Outer_XYplane[iBlk-76][0],
// 				      Input.NCells_Kdir,
// 				      STRETCHING_FCN_LINEAR,
// 				      ZERO,//1.10,//1.25,
// 				      0.125*Input.Height_Bunsen_Burner,
// 				      0.25*Input.Height_Bunsen_Burner);

//                Grid_Blks[iBlk].Set_BCs(BC_NONE,
// 				       BC_NONE,
// 				       BC_OUTFLOW_SUBSONIC,
// 				       BC_NONE,
// 				       BC_NONE,
// 				       BC_NONE);

//            } else if (iBlk >= 80 && iBlk <= 83) {
//               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Outer_XYplane[iBlk-80][0],
// 				      Input.NCells_Kdir,
// 				      STRETCHING_FCN_LINEAR,
// 				      ZERO,//1.10,//1.25,
// 				      0.25*Input.Height_Bunsen_Burner,
// 				      0.375*Input.Height_Bunsen_Burner);

//                Grid_Blks[iBlk].Set_BCs(BC_NONE,
// 				       BC_NONE,
// 				       BC_OUTFLOW_SUBSONIC,
// 				       BC_NONE,
// 				       BC_NONE,
// 				       BC_NONE);

//            } else if (iBlk >= 84 && iBlk <= 87) {
//               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Outer_XYplane[iBlk-84][0],
// 				      Input.NCells_Kdir,
// 				      STRETCHING_FCN_LINEAR,
// 				      ZERO,//1.10,//1.25,
// 				      0.375*Input.Height_Bunsen_Burner,
// 				      0.5*Input.Height_Bunsen_Burner);

//                Grid_Blks[iBlk].Set_BCs(BC_NONE,
// 				       BC_NONE,
// 				       BC_OUTFLOW_SUBSONIC,
// 				       BC_NONE,
// 				       BC_NONE,
// 				       BC_NONE);

//            } else if (iBlk >= 88 && iBlk <= 91) {
//               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Outer_XYplane[iBlk-88][0],
// 				      Input.NCells_Kdir,
// 				      STRETCHING_FCN_LINEAR,
// 				      ZERO,//1.10,//1.25,
// 				      0.5*Input.Height_Bunsen_Burner,
// 				      0.625*Input.Height_Bunsen_Burner);

//                Grid_Blks[iBlk].Set_BCs(BC_NONE,
// 				       BC_NONE,
// 				       BC_OUTFLOW_SUBSONIC,
// 				       BC_NONE,
// 				       BC_NONE,
// 				       BC_NONE);

//            } else if (iBlk >= 92 && iBlk <= 95) {
//               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Outer_XYplane[iBlk-92][0],
// 				      Input.NCells_Kdir,
// 				      STRETCHING_FCN_LINEAR,
// 				      ZERO,//1.10,//1.25,
// 				      0.625*Input.Height_Bunsen_Burner,
// 				      0.75*Input.Height_Bunsen_Burner);

//                Grid_Blks[iBlk].Set_BCs(BC_NONE,
// 				       BC_NONE,
// 				       BC_OUTFLOW_SUBSONIC,
// 				       BC_NONE,
// 				       BC_NONE,
// 				       BC_NONE);

//            } else if (iBlk >= 96 && iBlk <= 99) {
//               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Outer_XYplane[iBlk-96][0],
// 				      Input.NCells_Kdir,
// 				      STRETCHING_FCN_LINEAR,
// 				      ZERO,//1.10,//1.25,
// 				      0.75*Input.Height_Bunsen_Burner,
// 				      0.875*Input.Height_Bunsen_Burner);

//                Grid_Blks[iBlk].Set_BCs(BC_NONE,
// 				       BC_NONE,
// 				       BC_OUTFLOW_SUBSONIC,
// 				       BC_NONE,
// 				       BC_NONE,
// 				       BC_NONE);

//            } else if (iBlk >= 100 && iBlk <= 103) {
//               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Outer_XYplane[iBlk-100][0],
// 				      Input.NCells_Kdir,
// 				      STRETCHING_FCN_LINEAR,
// 				      ZERO,//1.10,//1.25,
// 				      0.875*Input.Height_Bunsen_Burner,
// 				      Input.Height_Bunsen_Burner);

//                Grid_Blks[iBlk].Set_BCs(BC_NONE,
// 				       BC_NONE,
// 				       BC_OUTFLOW_SUBSONIC,
// 				       BC_NONE,
// 				       BC_OUTFLOW_SUBSONIC,
// 				       BC_NONE);

            // iBlk = 90;
            /* fuel lines */
            
            //           if (iBlk <= 4) {
            //              Grid_Blks[iBlk].Extrude(Grid2D_Fuel_Line_XYplane[iBlk][0],
            // 				     Input.NCells_Kdir,
            // 				     STRETCHING_FCN_LINEAR,//_FCN_MIN_CLUSTERING,
            // 				     ZERO,//1.10,//1.25,
            // 				     ZERO*Input.Height_Bunsen_Burner,
            // 				     0.1*Input.Height_Bunsen_Burner);
            
            // 	      Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 		                      BC_NONE,
            //                                       BC_NONE,
            // 		                      BC_NONE,
            //                                       BC_NONE,
            //                                       BC_INFLOW_TURBULENCE);
            
            //            } else if (iBlk >= 5 && iBlk <= 9) {
            //              Grid_Blks[iBlk].Extrude(Grid2D_Fuel_Line_XYplane[iBlk-5][0],
            // 				     Input.NCells_Kdir,
            // 				     STRETCHING_FCN_LINEAR,//_FCN_MIN_CLUSTERING,
            // 				     ZERO,//1.10,//1.25,
            // 				     0.1*Input.Height_Bunsen_Burner,
            // 				     0.2*Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE);
            
            //            } else if (iBlk >= 10 && iBlk <= 14) {
            //              Grid_Blks[iBlk].Extrude(Grid2D_Fuel_Line_XYplane[iBlk-10][0],
            // 				     Input.NCells_Kdir,
            // 				     STRETCHING_FCN_LINEAR,//_FCN_MIN_CLUSTERING,
            // 				     ZERO,//1.10,//1.25,
            // 				     0.2*Input.Height_Bunsen_Burner,
            // 				     0.3*Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE);
            
            //            } else if (iBlk >= 15 && iBlk <= 19) {
            //              Grid_Blks[iBlk].Extrude(Grid2D_Fuel_Line_XYplane[iBlk-15][0],
            // 				     Input.NCells_Kdir,
            // 				     STRETCHING_FCN_LINEAR,//_FCN_MIN_CLUSTERING,
            // 				     ZERO,//1.10,//1.25,
            // 				     0.3*Input.Height_Bunsen_Burner,
            // 				     0.4*Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE);
            
            //            } else if (iBlk >= 20 && iBlk <= 24) {
            //              Grid_Blks[iBlk].Extrude(Grid2D_Fuel_Line_XYplane[iBlk-20][0],
            // 				     Input.NCells_Kdir,
            // 				     STRETCHING_FCN_LINEAR,//_FCN_MIN_CLUSTERING,
            // 				     ZERO,//1.10,//1.25,
            // 				     0.4*Input.Height_Bunsen_Burner,
            // 				     0.5*Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE);
            
            //            } else if (iBlk >= 25 && iBlk <= 29) {
            //              Grid_Blks[iBlk].Extrude(Grid2D_Fuel_Line_XYplane[iBlk-25][0],
            // 				     Input.NCells_Kdir,
            // 				     STRETCHING_FCN_LINEAR,//_FCN_MIN_CLUSTERING,
            // 				     ZERO,//1.10,//1.25,
            // 				     0.5*Input.Height_Bunsen_Burner,
            // 				     0.6*Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE);
            
            //            } else if (iBlk >= 30 && iBlk <= 34) {
            //              Grid_Blks[iBlk].Extrude(Grid2D_Fuel_Line_XYplane[iBlk-30][0],
            // 				     Input.NCells_Kdir,
            // 				     STRETCHING_FCN_LINEAR,//_FCN_MIN_CLUSTERING,
            // 				     ZERO,//1.10,//1.25,
            // 				     0.6*Input.Height_Bunsen_Burner,
            // 				     0.7*Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE);
            
            //            } else if (iBlk >= 35 && iBlk <= 39) {
            //              Grid_Blks[iBlk].Extrude(Grid2D_Fuel_Line_XYplane[iBlk-35][0],
            // 				     Input.NCells_Kdir,
            // 				     STRETCHING_FCN_LINEAR,//_FCN_MIN_CLUSTERING,
            // 				     ZERO,//1.10,//1.25,
            // 				     0.7*Input.Height_Bunsen_Burner,
            // 				     0.8*Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE);
            
            //            } else if (iBlk >= 40 && iBlk <= 44) {
            //              Grid_Blks[iBlk].Extrude(Grid2D_Fuel_Line_XYplane[iBlk-40][0],
            // 				     Input.NCells_Kdir,
            // 				     STRETCHING_FCN_LINEAR,//_FCN_MIN_CLUSTERING,
            // 				     ZERO,//1.10,//1.25,
            // 				     0.8*Input.Height_Bunsen_Burner,
            // 				     0.9*Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE);
            
            //            } else if (iBlk >= 45 && iBlk <= 49) {
            //              Grid_Blks[iBlk].Extrude(Grid2D_Fuel_Line_XYplane[iBlk-45][0],
            // 				     Input.NCells_Kdir,
            // 				     STRETCHING_FCN_LINEAR,//_FCN_MIN_CLUSTERING,
            // 				     ZERO,//1.10,//1.25,
            // 				     0.9*Input.Height_Bunsen_Burner,
            // 				     Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_OUTFLOW_SUBSONIC,
            // 				       BC_NONE);
            
            // 	      /* air line */
            
            //            } else if (iBlk >= 50 && iBlk <= 53) {
            //               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Inner_XYplane[iBlk-50][0],
            // 				      Input.NCells_Kdir,
            // 				      STRETCHING_FCN_LINEAR,//_FCN_MIN_CLUSTERING,
            // 				      ZERO,//1.10,//1.25,
            // 				      ZERO*Input.Height_Bunsen_Burner,
            // 				      0.1*Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_OUTFLOW_SUBSONIC,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_INFLOW_SUBSONIC);
            
            //            } else if (iBlk >= 54 && iBlk <= 57) {
            //               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Inner_XYplane[iBlk-54][0],
            // 				      Input.NCells_Kdir,
            // 				      STRETCHING_FCN_LINEAR,
            // 				      ZERO,//1.10,//1.25,
            // 				      0.1*Input.Height_Bunsen_Burner,
            // 				      0.2*Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_OUTFLOW_SUBSONIC,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE);
            
            //            } else if (iBlk >= 58 && iBlk <= 61) {
            //               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Inner_XYplane[iBlk-58][0],
            // 				      Input.NCells_Kdir,
            // 				      STRETCHING_FCN_LINEAR,
            // 				      ZERO,//1.10,//1.25,
            // 				      0.2*Input.Height_Bunsen_Burner,
            // 				      0.3*Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_OUTFLOW_SUBSONIC,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE);
            
            //            } else if (iBlk >= 62 && iBlk <= 65) {
            //               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Inner_XYplane[iBlk-62][0],
            // 				      Input.NCells_Kdir,
            // 				      STRETCHING_FCN_LINEAR,
            // 				      ZERO,//1.10,//1.25,
            // 				      0.3*Input.Height_Bunsen_Burner,
            // 				      0.4*Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_OUTFLOW_SUBSONIC,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE);
            
            //            } else if (iBlk >= 66 && iBlk <= 69) {
            //               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Inner_XYplane[iBlk-66][0],
            // 				      Input.NCells_Kdir,
            // 				      STRETCHING_FCN_LINEAR,
            // 				      ZERO,//1.10,//1.25,
            // 				      0.4*Input.Height_Bunsen_Burner,
            // 				      0.5*Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_OUTFLOW_SUBSONIC,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE);
            
            //            } else if (iBlk >= 70 && iBlk <= 73) {
            //               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Inner_XYplane[iBlk-70][0],
            // 				      Input.NCells_Kdir,
            // 				      STRETCHING_FCN_LINEAR,
            // 				      ZERO,//1.10,//1.25,
            // 				      0.5*Input.Height_Bunsen_Burner,
            // 				      0.6*Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_OUTFLOW_SUBSONIC,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE);
            
            //            } else if (iBlk >= 74 && iBlk <= 77) {
            //               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Inner_XYplane[iBlk-74][0],
            // 				      Input.NCells_Kdir,
            // 				      STRETCHING_FCN_LINEAR,
            // 				      ZERO,//1.10,//1.25,
            // 				      0.6*Input.Height_Bunsen_Burner,
            // 				      0.7*Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_OUTFLOW_SUBSONIC,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE);
            
            //            } else if (iBlk >= 78 && iBlk <= 81) {
            //               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Inner_XYplane[iBlk-78][0],
            // 				      Input.NCells_Kdir,
            // 				      STRETCHING_FCN_LINEAR,
            // 				      ZERO,//1.10,//1.25,
            // 				      0.7*Input.Height_Bunsen_Burner,
            // 				      0.8*Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_OUTFLOW_SUBSONIC,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE);
            
            //            } else if (iBlk >= 82 && iBlk <= 85) {
            //               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Inner_XYplane[iBlk-82][0],
            // 				      Input.NCells_Kdir,
            // 				      STRETCHING_FCN_LINEAR,
            // 				      ZERO,//1.10,//1.25,
            // 				      0.8*Input.Height_Bunsen_Burner,
            // 				      0.9*Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_OUTFLOW_SUBSONIC,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE);
            
            //            } else if (iBlk >= 86 && iBlk <= 89) {
            //               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Inner_XYplane[iBlk-86][0],
            // 				      Input.NCells_Kdir,
            // 				      STRETCHING_FCN_LINEAR,
            // 				      ZERO,//1.10,//1.25,
            // 				      0.9*Input.Height_Bunsen_Burner,
            // 				      Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_OUTFLOW_SUBSONIC,
            // 				       BC_NONE,
            // 				       BC_INFLOW_SUBSONIC,
            // 				       BC_NONE);
            
            // iBlk = 52;
            /* fuel lines */
            
            //           if (iBlk <= 4) {
            //              Grid_Blks[iBlk].Extrude(Grid2D_Fuel_Line_XYplane[iBlk][0],
            // 				     Input.NCells_Kdir,
            // 				     STRETCHING_FCN_MIN_CLUSTERING,
            // 				     1.25,
            // 				     ZERO*Input.Height_Bunsen_Burner,
            // 				     0.25*Input.Height_Bunsen_Burner);
            
            // 	      Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 		                      BC_NONE,
            //                                       BC_NONE,
            // 		                      BC_NONE,
            //                                       BC_NONE,
            //                                       BC_INFLOW_TURBULENCE);
            
            //            } else if (iBlk >= 5 && iBlk <= 9) {
            //              Grid_Blks[iBlk].Extrude(Grid2D_Fuel_Line_XYplane[iBlk-5][0],
            // 				     Input.NCells_Kdir,
            // 				     STRETCHING_FCN_MIN_CLUSTERING,
            // 				     1.25,
            // 				     0.25*Input.Height_Bunsen_Burner,
            // 				     0.5*Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE);
            
            //            } else if (iBlk >= 10 && iBlk <= 14) {
            //              Grid_Blks[iBlk].Extrude(Grid2D_Fuel_Line_XYplane[iBlk-10][0],
            // 				     Input.NCells_Kdir,
            // 				     STRETCHING_FCN_MIN_CLUSTERING,
            // 				     1.25,
            // 				     0.5*Input.Height_Bunsen_Burner,
            // 				     0.75*Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE);
            
            //            } else if (iBlk >= 15 && iBlk <= 19) {
            //              Grid_Blks[iBlk].Extrude(Grid2D_Fuel_Line_XYplane[iBlk-15][0],
            // 				     Input.NCells_Kdir,
            // 				     STRETCHING_FCN_MIN_CLUSTERING,
            // 				     1.25,
            // 				     0.75*Input.Height_Bunsen_Burner,
            // 				     Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_FIXED_PRESSURE,
            // 				       BC_NONE);
            
            // 	      /* first ring for air */
            
            //            } else if (iBlk >= 20 && iBlk <= 23) {
            //               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Inner_XYplane[iBlk-20][0],
            // 				      Input.NCells_Kdir,
            // 				      STRETCHING_FCN_MIN_CLUSTERING,
            // 				      1.25,
            // 				      ZERO*Input.Height_Bunsen_Burner,
            // 				      0.25*Input.Height_Bunsen_Burner);
            
            // 	      Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 		                      BC_NONE,
            //                                       BC_NONE,
            // 		                      BC_NONE,
            //                                       BC_NONE,
            //                                       BC_DIRICHLET);
            
            // 	   } else if (iBlk >= 24 && iBlk <= 27) {
            //               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Inner_XYplane[iBlk-24][0],
            // 				      Input.NCells_Kdir,
            // 				      STRETCHING_FCN_LINEAR,
            // 				      1.25,
            // 				      0.25*Input.Height_Bunsen_Burner,
            // 				      0.5*Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE);
            
            // 	   } else if (iBlk >= 28 && iBlk <= 31) {
            //               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Inner_XYplane[iBlk-28][0],
            // 				      Input.NCells_Kdir,
            // 				      STRETCHING_FCN_LINEAR,
            // 				      1.25,
            // 				      0.5*Input.Height_Bunsen_Burner,
            // 				      0.75*Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE);
            
            // 	   } else if (iBlk >= 32 && iBlk <= 35) {
            //               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Inner_XYplane[iBlk-32][0],
            // 				      Input.NCells_Kdir,
            // 				      STRETCHING_FCN_LINEAR,
            // 				      1.25,
            // 				      0.75*Input.Height_Bunsen_Burner,
            // 				      Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_FIXED_PRESSURE,
            // 				       BC_NONE);
            
            // 	      /* second ring for air */
            
            //            } else if (iBlk >= 36 && iBlk <= 39) {
            //               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Outer_XYplane[iBlk-36][0],
            // 				      Input.NCells_Kdir,
            // 				      STRETCHING_FCN_MIN_CLUSTERING,
            // 				      1.25,
            // 				      ZERO*Input.Height_Bunsen_Burner,
            // 				      0.25*Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_REFLECTION,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_DIRICHLET);
            
            //            } else if (iBlk >= 40 && iBlk <= 43) {
            //               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Outer_XYplane[iBlk-40][0],
            // 				      Input.NCells_Kdir,
            // 				      STRETCHING_FCN_LINEAR,
            // 				      1.25,
            // 				      0.25*Input.Height_Bunsen_Burner,
            // 				      0.5*Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_REFLECTION,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE);
            
            //            } else if (iBlk >= 44 && iBlk <= 47) {
            //               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Outer_XYplane[iBlk-44][0],
            // 				      Input.NCells_Kdir,
            // 				      STRETCHING_FCN_LINEAR,
            // 				      1.25,
            // 				      0.5*Input.Height_Bunsen_Burner,
            // 				      0.75*Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_REFLECTION,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE);
            
            //            } else if (iBlk >= 48 && iBlk <= 51) {
            //               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Outer_XYplane[iBlk-48][0],
            // 				      Input.NCells_Kdir,
            // 				      STRETCHING_FCN_LINEAR,
            // 				      1.25,
            // 				      0.75*Input.Height_Bunsen_Burner,
            // 				      Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_REFLECTION,
            // 				       BC_NONE,
            // 				       BC_FIXED_PRESSURE,
            // 				       BC_NONE);
            
            //iBlk = 39;
            /* fuel lines */
            
            //           if (iBlk <= 4) {
            //              Grid_Blks[iBlk].Extrude(Grid2D_Fuel_Line_XYplane[iBlk][0],
            // 				     Input.NCells_Kdir,
            // 				     STRETCHING_FCN_MIN_CLUSTERING,
            // 				     1.25,
            // 				     ZERO*Input.Height_Bunsen_Burner,
            // 				     1.0/3.0*Input.Height_Bunsen_Burner);
            
            // 	      Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 		                      BC_NONE,
            //                                       BC_NONE,
            // 		                      BC_NONE,
            //                                       BC_NONE,
            //                                       BC_DIRICHLET);
            
            //            } else if (iBlk >= 5 && iBlk <= 9) {
            //              Grid_Blks[iBlk].Extrude(Grid2D_Fuel_Line_XYplane[iBlk-5][0],
            // 				     Input.NCells_Kdir,
            // 				     STRETCHING_FCN_MIN_CLUSTERING,
            // 				     1.25,
            // 				     1.0/3.0*Input.Height_Bunsen_Burner,
            // 				     2.0/3.0*Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE);
            
            //            } else if (iBlk >= 10 && iBlk <= 14) {
            //              Grid_Blks[iBlk].Extrude(Grid2D_Fuel_Line_XYplane[iBlk-10][0],
            // 				     Input.NCells_Kdir,
            // 				     STRETCHING_FCN_MIN_CLUSTERING,
            // 				     1.25,
            // 				     2.0/3.0*Input.Height_Bunsen_Burner,
            // 				     Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_FIXED_PRESSURE,
            // 				       BC_NONE);
            
            // 	      /* first ring for air */
            
            //            } else if (iBlk >= 15 && iBlk <= 18) {
            //               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Inner_XYplane[iBlk-15][0],
            // 				      Input.NCells_Kdir,
            // 				      STRETCHING_FCN_MIN_CLUSTERING,
            // 				      1.25,
            // 				      ZERO*Input.Height_Bunsen_Burner,
            // 				      1.0/3.0*Input.Height_Bunsen_Burner);
            
            // 	      Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 		                      BC_NONE,
            //                                       BC_NONE,
            // 		                      BC_NONE,
            //                                       BC_NONE,
            //                                       BC_DIRICHLET);
            
            // 	   } else if (iBlk >= 19 && iBlk <= 22) {
            //               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Inner_XYplane[iBlk-19][0],
            // 				      Input.NCells_Kdir,
            // 				      STRETCHING_FCN_LINEAR,
            // 				      1.25,
            // 				      1.0/3.0*Input.Height_Bunsen_Burner,
            // 				      2.0/3.0*Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE);
            
            // 	   } else if (iBlk >= 23 && iBlk <= 26) {
            //               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Inner_XYplane[iBlk-23][0],
            // 				      Input.NCells_Kdir,
            // 				      STRETCHING_FCN_LINEAR,
            // 				      1.25,
            // 				      2.0/3.0*Input.Height_Bunsen_Burner,
            // 				      Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_FIXED_PRESSURE,
            // 				       BC_NONE);
            
            // 	      /* second ring for air */
            
            //            } else if (iBlk >= 27 && iBlk <= 30) {
            //               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Outer_XYplane[iBlk-27][0],
            // 				      Input.NCells_Kdir,
            // 				      STRETCHING_FCN_MIN_CLUSTERING,
            // 				      1.25,
            // 				      ZERO*Input.Height_Bunsen_Burner,
            // 				      1.0/3.0*Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_REFLECTION,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_DIRICHLET);
            
            //            } else if (iBlk >= 31 && iBlk <= 34) {
            //               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Outer_XYplane[iBlk-31][0],
            // 				      Input.NCells_Kdir,
            // 				      STRETCHING_FCN_LINEAR,
            // 				      1.25,
            // 				      1.0/3.0*Input.Height_Bunsen_Burner,
            // 				      2.0/3.0*Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_REFLECTION,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE);
            
            //            } else if (iBlk >= 35 && iBlk <= 38) {
            //               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Outer_XYplane[iBlk-35][0],
            // 				      Input.NCells_Kdir,
            // 				      STRETCHING_FCN_LINEAR,
            // 				      1.25,
            // 				      2.0/3.0*Input.Height_Bunsen_Burner,
            // 				      Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_REFLECTION,
            // 				       BC_NONE,
            // 				       BC_FIXED_PRESSURE,
            // 				       BC_NONE);
            
            // iBlk = 68;
            /* fuel lines */
            
            //           if (iBlk <= 4) {
            //              Grid_Blks[iBlk].Extrude(Grid2D_Fuel_Line_XYplane[iBlk][0],
            // 				     Input.NCells_Kdir,
            // 				     STRETCHING_FCN_MIN_CLUSTERING,
            // 				     1.25,
            // 				     ZERO*Input.Height_Bunsen_Burner,
            // 				     0.25*Input.Height_Bunsen_Burner);
            
            // 	      Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 		                      BC_NONE,
            //                                       BC_NONE,
            // 		                      BC_NONE,
            //                                       BC_NONE,
            //                                       BC_INFLOW_TURBULENCE);
            
            // //            } else if (iBlk >= 5 && iBlk <= 8) {
            // //              Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Inner_XYplane[iBlk-5][0],
            // // 				     Input.NCells_Kdir,
            // // 				     STRETCHING_FCN_MIN_CLUSTERING,
            // // 				     1.25,
            // // 				     ZERO*Input.Height_Bunsen_Burner,
            // // 				     0.25*Input.Height_Bunsen_Burner);
            
            // //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // // 				       BC_NONE,
            // // 				       BC_NONE,
            // // 				       BC_NONE,
            // // 				       BC_NONE,
            // // 				       BC_DIRICHLET);
            
            
            
            //            } else if (iBlk >= 5 && iBlk <= 9) {
            //              Grid_Blks[iBlk].Extrude(Grid2D_Fuel_Line_XYplane[iBlk-5][0],
            // 				     Input.NCells_Kdir,
            // 				     STRETCHING_FCN_MIN_CLUSTERING,
            // 				     1.25,
            // 				     0.25*Input.Height_Bunsen_Burner,
            // 				     0.5*Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE);
            
            //            } else if (iBlk >= 10 && iBlk <= 14) {
            //              Grid_Blks[iBlk].Extrude(Grid2D_Fuel_Line_XYplane[iBlk-10][0],
            // 				     Input.NCells_Kdir,
            // 				     STRETCHING_FCN_MIN_CLUSTERING,
            // 				     1.25,
            // 				     0.5*Input.Height_Bunsen_Burner,
            // 				     0.75*Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE);
            
            //            } else if (iBlk >= 15 && iBlk <= 19) {
            //              Grid_Blks[iBlk].Extrude(Grid2D_Fuel_Line_XYplane[iBlk-15][0],
            // 				     Input.NCells_Kdir,
            // 				     STRETCHING_FCN_MIN_CLUSTERING,
            // 				     1.25,
            // 				     0.75*Input.Height_Bunsen_Burner,
            // 				     Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_FIXED_PRESSURE,
            // 				       BC_NONE);
            
            // 	      /* first ring for air */
            
            //            } else if (iBlk >= 20 && iBlk <= 23) {
            //               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Inner_XYplane[iBlk-20][0],
            // 				      Input.NCells_Kdir,
            // 				      STRETCHING_FCN_MIN_CLUSTERING,
            // 				      1.25,
            // 				      ZERO*Input.Height_Bunsen_Burner,
            // 				      0.25*Input.Height_Bunsen_Burner);
            
            // 	      Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 		                      BC_NONE,
            //                                       BC_NONE,
            // 		                      BC_NONE,
            //                                       BC_NONE,
            //                                       BC_DIRICHLET);
            
            // 	   } else if (iBlk >= 24 && iBlk <= 27) {
            //               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Inner_XYplane[iBlk-24][0],
            // 				      Input.NCells_Kdir,
            // 				      STRETCHING_FCN_LINEAR,
            // 				      1.25,
            // 				      0.25*Input.Height_Bunsen_Burner,
            // 				      0.5*Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE);
            
            // 	   } else if (iBlk >= 28 && iBlk <= 31) {
            //               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Inner_XYplane[iBlk-28][0],
            // 				      Input.NCells_Kdir,
            // 				      STRETCHING_FCN_LINEAR,
            // 				      1.25,
            // 				      0.5*Input.Height_Bunsen_Burner,
            // 				      0.75*Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE);
            
            // 	   } else if (iBlk >= 32 && iBlk <= 35) {
            //               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Inner_XYplane[iBlk-32][0],
            // 				      Input.NCells_Kdir,
            // 				      STRETCHING_FCN_LINEAR,
            // 				      1.25,
            // 				      0.75*Input.Height_Bunsen_Burner,
            // 				      Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_FIXED_PRESSURE,
            // 				       BC_NONE);
            
            // 	      /* second ring for air */
            
            //            } else if (iBlk >= 36 && iBlk <= 39) {
            //               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Middle_XYplane[iBlk-36][0],
            // 				      Input.NCells_Kdir,
            // 				      STRETCHING_FCN_MIN_CLUSTERING,
            // 				      1.25,
            // 				      ZERO*Input.Height_Bunsen_Burner,
            // 				      0.25*Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_DIRICHLET);
            
            //            } else if (iBlk >= 40 && iBlk <= 43) {
            //               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Middle_XYplane[iBlk-40][0],
            // 				      Input.NCells_Kdir,
            // 				      STRETCHING_FCN_LINEAR,
            // 				      1.25,
            // 				      0.25*Input.Height_Bunsen_Burner,
            // 				      0.5*Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE);
            
            //            } else if (iBlk >= 44 && iBlk <= 47) {
            //               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Middle_XYplane[iBlk-44][0],
            // 				      Input.NCells_Kdir,
            // 				      STRETCHING_FCN_LINEAR,
            // 				      1.25,
            // 				      0.5*Input.Height_Bunsen_Burner,
            // 				      0.75*Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE);
            
            //            } else if (iBlk >= 48 && iBlk <= 51) {
            //               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Middle_XYplane[iBlk-48][0],
            // 				      Input.NCells_Kdir,
            // 				      STRETCHING_FCN_LINEAR,
            // 				      1.25,
            // 				      0.75*Input.Height_Bunsen_Burner,
            // 				      Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_FIXED_PRESSURE,
            // 				       BC_NONE);
            // 	      /* third ring for air */
            
            //            } else if (iBlk >= 52 && iBlk <= 55) {
            //               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Outer_XYplane[iBlk-52][0],
            // 				      Input.NCells_Kdir,
            // 				      STRETCHING_FCN_MIN_CLUSTERING,
            // 				      1.25,
            // 				      ZERO*Input.Height_Bunsen_Burner,
            // 				      0.25*Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_FIXED_PRESSURE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_DIRICHLET);
            
            //            } else if (iBlk >= 56 && iBlk <= 59) {
            //               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Outer_XYplane[iBlk-56][0],
            // 				      Input.NCells_Kdir,
            // 				      STRETCHING_FCN_LINEAR,
            // 				      1.25,
            // 				      0.25*Input.Height_Bunsen_Burner,
            // 				      0.5*Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_FIXED_PRESSURE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE);
            
            //            } else if (iBlk >= 60 && iBlk <= 63) {
            //               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Outer_XYplane[iBlk-60][0],
            // 				      Input.NCells_Kdir,
            // 				      STRETCHING_FCN_LINEAR,
            // 				      1.25,
            // 				      0.5*Input.Height_Bunsen_Burner,
            // 				      0.75*Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_FIXED_PRESSURE,
            // 				       BC_NONE,
            // 				       BC_NONE,
            // 				       BC_NONE);
            
            //            } else if (iBlk >= 64 && iBlk <= 67) {
            //               Grid_Blks[iBlk].Extrude(Grid2D_Bunsen_Burner_Outer_XYplane[iBlk-64][0],
            // 				      Input.NCells_Kdir,
            // 				      STRETCHING_FCN_LINEAR,
            // 				      1.25,
            // 				      0.75*Input.Height_Bunsen_Burner,
            // 				      Input.Height_Bunsen_Burner);
            
            //                Grid_Blks[iBlk].Set_BCs(BC_NONE,
            // 				       BC_NONE,
            // 				       BC_FIXED_PRESSURE,
            // 				       BC_NONE,
            // 				       BC_FIXED_PRESSURE,
            // 				       BC_NONE);
        } /* endif */
	} /* endfor */
    
    /* Deallocate 2D grid. */
    
    Grid2D_Fuel_Line_XYplane = Deallocate_Multi_Block_Grid(Grid2D_Fuel_Line_XYplane,
                                                           numblk_idir_fuel,
                                                           numblk_jdir_fuel);
    
    Grid2D_Bunsen_Burner_Inner_XYplane = Deallocate_Multi_Block_Grid(Grid2D_Bunsen_Burner_Inner_XYplane,
                                                                     numblk_idir_bunsenburner,
                                                                     numblk_jdir_bunsenburner);
    
    //     Grid2D_Bunsen_Burner_Middle_XYplane = Deallocate_Multi_Block_Grid(Grid2D_Bunsen_Burner_Middle_XYplane,
    // 								      numblk_idir_bunsenburner,
    // 								      numblk_jdir_bunsenburner);
    
    //     Grid2D_Bunsen_Burner_Outer_XYplane = Deallocate_Multi_Block_Grid(Grid2D_Bunsen_Burner_Outer_XYplane,
    // 								     numblk_idir_bunsenburner,
    // 								     numblk_jdir_bunsenburner);
    
    /* Call the function Find_Neighbours to obtain the neighbour block information
     and assign values to data members in the grid block connectivity data structure. */
    
    Find_Neighbours(Input);
    
}


/**************************************************************
 * Routine: Create_Flat_Plate                                 *
 *                                                            *
 * Generates a 3D Cartesian multiblock mesh for a Flat Plate. *
 *                                                            *
 **************************************************************/
void Grid3D_Hexa_Multi_Block_List::Create_Grid_Flat_Plate(Grid3D_Input_Parameters &Input) {
    
    int count_blocks;
    int BC_top, BC_bottom;
    Grid2D_Quad_Block **Grid2D_Flat_Plate_XYplane;
    
    /* Allocate required memory. */
    
    Allocate(Input.NBlk_Idir, Input.NBlk_Jdir, Input.NBlk_Kdir);
    
    /* Creat 2D cross-section grids from which the 3D grid
     will be extruded. */
    
    Grid2D_Flat_Plate_XYplane =  Grid_Flat_Plate(Grid2D_Flat_Plate_XYplane,
                                                 Input.NBlk_Idir, 
                                                 Input.NBlk_Jdir,
                                                 Input.Plate_Length,
                                                 Input.Box_Height,
                                                 BC_ADIABATIC_WALL,
                                                 ON,
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
                
                Grid_Blks[count_blocks].Extrude(Grid2D_Flat_Plate_XYplane[iBlk][jBlk],
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
    
    Grid2D_Flat_Plate_XYplane = Deallocate_Multi_Block_Grid(Grid2D_Flat_Plate_XYplane,
                                                     Input.NBlk_Idir, 
                                                     Input.NBlk_Jdir);
    
    
    /* Call the function Find_Neighbours to obtain the neighbour block information
     and assign values to data members in the grid block connectivity data structure. */
    
    Find_Neighbours(Input);
    
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
                        BC_east = BC_PERIODIC;
                    } else {
                        BC_east = BC_NONE;
                    } /* endif */
                    if (iBlk == 0) {
                        BC_west = BC_PERIODIC;
                    } else {
                        BC_west = BC_NONE;
                    } /* endif */
                    
                    if (jBlk == Input.NBlk_Jdir-1) {
                        BC_north = BC_ADIABATIC_WALL;
                    } else {
                        BC_north = BC_NONE;
                    } /* endif */
                    if (jBlk == 0) {
                        BC_south = BC_ADIABATIC_WALL;
                    } else {
                        BC_south = BC_NONE;
                    } /* endif */
                    
                    if (kBlk == Input.NBlk_Kdir-1) {
                        BC_top = BC_OUTFLOW_SUBSONIC;
                    } else {
                        BC_top = BC_NONE;
                    } /* endif */
                    if (kBlk == 0) {
                        BC_bottom = BC_INFLOW_SUBSONIC;
                    } else {
                        BC_bottom = BC_NONE;
                    } /* endif */
                    
                } else if (Input.i_Grid == GRID_CHANNEL_XDIR) {
                    
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
                        BC_top = BC_ADIABATIC_WALL;
                    } else {
                        BC_top = BC_NONE;
                    } /* endif */
                    if (kBlk == 0) {
                        BC_bottom = BC_ADIABATIC_WALL;
                    } else {
                        BC_bottom = BC_NONE;
                    } /* endif */
                    
                } else if (Input.i_Grid == GRID_CHANNEL_YDIR) {
                    
                    if (iBlk == Input.NBlk_Idir-1) {
                        BC_east = BC_ADIABATIC_WALL;
                    } else {
                        BC_east = BC_NONE;
                    } /* endif */
                    if (iBlk == 0) {
                        BC_west = BC_ADIABATIC_WALL;
                    } else {
                        BC_west = BC_NONE;
                    } /* endif */
                    
                    if (jBlk == Input.NBlk_Jdir-1) {
                        BC_north = BC_OUTFLOW_SUBSONIC;
                    } else {
                        BC_north = BC_NONE;
                    } /* endif */
                    if (jBlk == 0) {
                        BC_south = BC_INFLOW_SUBSONIC;
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
    
    
    /* Call the function Find_Neighbours to obtain the neighbour block information
     and assign values to data members in the grid block connectivity data structure. */
    
    Find_Neighbours(Input);
    
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
                        BC_east = BC_MOVING_WALL;
                    } else {
                        BC_east = BC_NONE;
                    } /* endif */
                    if (iBlk == 0) {
                        BC_west = BC_NO_SLIP;
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
                        BC_top = BC_CONSTANT_EXTRAPOLATION;
                    } else {
                        BC_top = BC_NONE;
                    } /* endif */
                    if (kBlk == 0) {
                        BC_bottom = BC_CONSTANT_EXTRAPOLATION;
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
    
    
    /* Call the function Find_Neighbours to obtain the neighbour block information
     and assign values to data members in the grid block connectivity data structure. */
    
    Find_Neighbours(Input);
    
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
    
    Grid2D_Pipe_XYplane = Grid_Tube_2D(Grid2D_Pipe_XYplane,
                                       numblk_idir_pipe,
                                       numblk_jdir_pipe,
                                       Input.Pipe_Radius,
                                       Input.NCells_Idir,
                                       Input.NCells_Jdir,
                                       Input.Nghost,
                                       STRETCHING_FCN_MAX_CLUSTERING,
                                       Input.Stretching_Factor_Jdir);
    
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
                                    BC_CONSTANT_EXTRAPOLATION,
                                    BC_FIXED_PRESSURE);
        } else {
            Grid_Blks[iBlk].Set_BCs(BC_NONE,
                                    BC_NONE,
                                    BC_NO_SLIP,
                                    BC_NONE,
                                    BC_CONSTANT_EXTRAPOLATION,
                                    BC_FIXED_PRESSURE);
        } /* endif */
        
    } /* endfor */
    
    /* Deallocate 2D grid. */
    
    Grid2D_Pipe_XYplane = Deallocate_Multi_Block_Grid(Grid2D_Pipe_XYplane,
                                                      numblk_idir_pipe,
                                                      numblk_jdir_pipe);
    
    
    /* Call the function Find_Neighbours to obtain the neighbour block information
     and assign values to data members in the grid block connectivity data structure. */
    
    Find_Neighbours(Input);
    
}

/**********************************************************************
 * Routine: Grid_Bump_Channel_Flow                                    *
 *                                                                    *
 * Generates a single block quadilateral mesh with clustering for     *
 * predicting supersonic flow around a cirucular cylinder blunt body. *
 *                                                                    *
 **********************************************************************/
void Grid3D_Hexa_Multi_Block_List::Create_Grid_Bump_Channel_Flow(Grid3D_Input_Parameters &Input) {
    
    int count_blocks;
    
    Grid2D_Quad_Block **Grid2D_XYplane;
    int BC_top, BC_bottom;
    
    //Fixed at 8 blocks (4x2)
    int Number_of_Blocks_Idir, Number_of_Blocks_Jdir;
    int Smooth_Bump = 0; //don't use smooth bump.
    
    /* Allocate required memory. */
    
    Input.NBlk_Idir = 4;
    Input.NBlk_Jdir = 2;
    Input.NBlk_Kdir = 1;
    
    Allocate(Input.NBlk_Idir, Input.NBlk_Jdir, Input.NBlk_Kdir);
    
    /* Creat 2D cross-section grids from which the 3D grid
     will be extruded. */
    
    Grid2D_XYplane = Grid_Bump_Channel_Flow(Grid2D_XYplane,
                                            Number_of_Blocks_Idir,
                                            Number_of_Blocks_Jdir,
                                            Smooth_Bump,
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
                
                Grid_Blks[count_blocks].Extrude(Grid2D_XYplane[iBlk][jBlk],
                                                Input.NCells_Kdir,
                                                Input.Stretching_Type_Kdir,
                                                Input.Stretching_Factor_Kdir,			       
                                                (double(kBlk)/double(Input.NBlk_Kdir))*Input.Box_Length,
                                                (double(kBlk+1)/double(Input.NBlk_Kdir))*Input.Box_Length);
                
                /* Assign top and bottom boundary conditions. */
                
                if (kBlk == Input.NBlk_Kdir-1) {
                    BC_top = BC_REFLECTION;//BC_CONSTANT_EXTRAPOLATION;
                } else {
                    BC_top = BC_NONE;
                } /* endif */
                if (kBlk == 0) {
                    BC_bottom =BC_REFLECTION;//BC_CONSTANT_EXTRAPOLATION;
                } else {
                    BC_bottom = BC_NONE;
                } /* endif */
                
                Grid_Blks[count_blocks].Set_BCs_Zdir(BC_top, BC_bottom);
                
                /* Update block counter. */
                
                count_blocks ++;
                
            } 
        } 
    } 
    
    /* Deallocate 2D grid. */
    
    Grid2D_XYplane = Deallocate_Multi_Block_Grid(Grid2D_XYplane,
                                                 Number_of_Blocks_Idir,
                                                 Number_of_Blocks_Jdir);
    
    
    /* Call the function Find_Neighbours to obtain the neighbour block information
     and assign values to data members in the grid block connectivity data structure. */
    
    Find_Neighbours(Input);
    
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
    Grid2D_Fuel_Line_XYplane = Grid_Tube_2D(Grid2D_Fuel_Line_XYplane,
                                            numblk_idir_fuel,
                                            numblk_jdir_fuel,
                                            Input.Radius_Fuel_Line,
                                            Input.NCells_Idir,
                                            Input.NCells_Jdir,
                                            Input.Nghost,
                                            STRETCHING_FCN_MAX_CLUSTERING,
                                            1.25);
    
    Grid2D_Bluff_Body_Inner_XYplane = Grid_Annulus_2D(Grid2D_Bluff_Body_Inner_XYplane,
                                                      numblk_idir_bluffbody,
                                                      numblk_jdir_bluffbody,
                                                      Input.Radius_Fuel_Line,
                                                      HALF*Input.Radius_Bluff_Body,
                                                      Input.NCells_Idir,
                                                      Input.NCells_Jdir,
                                                      Input.Nghost,
                                                      STRETCHING_FCN_MIN_CLUSTERING,
                                                      1.10);
    
    Grid2D_Bluff_Body_Outer_XYplane = Grid_Annulus_2D(Grid2D_Bluff_Body_Outer_XYplane,
                                                      numblk_idir_bluffbody,
                                                      numblk_jdir_bluffbody,
                                                      HALF*Input.Radius_Bluff_Body,
                                                      Input.Radius_Bluff_Body,
                                                      Input.NCells_Idir,
                                                      Input.NCells_Jdir,
                                                      Input.Nghost,
                                                      STRETCHING_FCN_MAX_CLUSTERING,
                                                      1.10);
    
    Grid2D_Coflow_XYplane = Grid_Annulus_2D(Grid2D_Coflow_XYplane,
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
                                    ZERO+Input.Length_Coflow_Inlet_Pipe,
                                    0.25*Input.Length_Combustor_Tube+Input.Length_Coflow_Inlet_Pipe);
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
                                    0.25*Input.Length_Combustor_Tube+Input.Length_Coflow_Inlet_Pipe,
                                    Input.Length_Combustor_Tube+Input.Length_Coflow_Inlet_Pipe);
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
                                    ZERO+Input.Length_Coflow_Inlet_Pipe,
                                    0.25*Input.Length_Combustor_Tube+Input.Length_Coflow_Inlet_Pipe);
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
                                    0.25*Input.Length_Combustor_Tube+Input.Length_Coflow_Inlet_Pipe,
                                    Input.Length_Combustor_Tube+Input.Length_Coflow_Inlet_Pipe);
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
                                    ZERO+Input.Length_Coflow_Inlet_Pipe,
                                    0.25*Input.Length_Combustor_Tube+Input.Length_Coflow_Inlet_Pipe);
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
                                    0.25*Input.Length_Combustor_Tube+Input.Length_Coflow_Inlet_Pipe,
                                    Input.Length_Combustor_Tube+Input.Length_Coflow_Inlet_Pipe);
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
                                    ZERO+Input.Length_Coflow_Inlet_Pipe,
                                    0.25*Input.Length_Combustor_Tube+Input.Length_Coflow_Inlet_Pipe);
            Grid_Blks[iBlk].Set_BCs(BC_NONE,
                                    BC_NONE,
                                    BC_REFLECTION,//BC_WALL_VISCOUS,//BC_REFLECTION,
                                    BC_NONE,
                                    BC_NONE,
                                    BC_NONE);
            
        } else if (iBlk >= 30 && iBlk <= 33) {
            Grid_Blks[iBlk].Extrude(Grid2D_Coflow_XYplane[iBlk-30][0],
                                    Input.NCells_Kdir,
                                    STRETCHING_FCN_LINEAR,
                                    ZERO,
                                    0.25*Input.Length_Combustor_Tube+Input.Length_Coflow_Inlet_Pipe,
                                    Input.Length_Combustor_Tube+Input.Length_Coflow_Inlet_Pipe);
            Grid_Blks[iBlk].Set_BCs(BC_NONE,
                                    BC_NONE,
                                    BC_REFLECTION,//BC_WALL_VISCOUS,//BC_REFLECTION,
                                    BC_NONE,
                                    BC_FIXED_PRESSURE,
                                    BC_NONE);
            
        } else if (iBlk >= 34 && iBlk <= 37) {
            Grid_Blks[iBlk].Extrude(Grid2D_Coflow_XYplane[iBlk-34][0],
                                    Input.NCells_Kdir,
                                    STRETCHING_FCN_MAX_CLUSTERING,
                                    1.25,
                                    -0.25*Input.Length_Coflow_Inlet_Pipe+Input.Length_Coflow_Inlet_Pipe,
                                    ZERO+Input.Length_Coflow_Inlet_Pipe);
            Grid_Blks[iBlk].Set_BCs(BC_NONE,
                                    BC_NONE,
                                    BC_REFLECTION,//BC_WALL_VISCOUS,//BC_REFLECTION,
                                    BC_WALL_VISCOUS,
                                    BC_NONE,
                                    BC_NONE);
            
        } else if (iBlk >= 38 && iBlk <= 41) {
            Grid_Blks[iBlk].Extrude(Grid2D_Coflow_XYplane[iBlk-38][0],
                                    Input.NCells_Kdir,
                                    STRETCHING_FCN_LINEAR,
                                    ONE,
                                    -Input.Length_Coflow_Inlet_Pipe+Input.Length_Coflow_Inlet_Pipe,
                                    -0.25*Input.Length_Coflow_Inlet_Pipe+Input.Length_Coflow_Inlet_Pipe);
            Grid_Blks[iBlk].Set_BCs(BC_NONE,
                                    BC_NONE,
                                    BC_REFLECTION,//BC_WALL_VISCOUS,//BC_REFLECTION,
                                    BC_WALL_VISCOUS,
                                    BC_DIRICHLET,
                                    BC_NONE);
            
        } /* endif */
    } /* endfor */
    
    /* Deallocate 2D grids. */
    
    Grid2D_Fuel_Line_XYplane = Deallocate_Multi_Block_Grid(Grid2D_Fuel_Line_XYplane,
                                                           numblk_idir_fuel,
                                                           numblk_jdir_fuel);
    
    Grid2D_Bluff_Body_Inner_XYplane = Deallocate_Multi_Block_Grid(Grid2D_Bluff_Body_Inner_XYplane,
                                                                  numblk_idir_bluffbody,
                                                                  numblk_jdir_bluffbody);
    
    Grid2D_Bluff_Body_Outer_XYplane = Deallocate_Multi_Block_Grid(Grid2D_Bluff_Body_Outer_XYplane,
                                                                  numblk_idir_bluffbody,
                                                                  numblk_jdir_bluffbody);
    
    Grid2D_Coflow_XYplane = Deallocate_Multi_Block_Grid(Grid2D_Coflow_XYplane,
                                                        numblk_idir_coflow,
                                                        numblk_jdir_coflow);
    
    /* Call the function Find_Neighbours to obtain the neighbour block information
     and assign values to data members in the grid block connectivity data structure. */
    
    Find_Neighbours(Input);
    
}

/********************************************************
 * Routine: Create_Grid_ICEMCFD                         *
 *                                                      *
 * Read ICEMCFD Mesh                                    *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block_List::Create_Grid_ICEMCFD(Grid3D_Input_Parameters &Input) {
    
    if (Allocated) Deallocate();
    
    assert(NBlk_Idir >= 1 && NBlk_Jdir >= 1 && NBlk_Kdir >= 1); 
    Allocated = 1;
    
    
    /* Call the function Find_Neighbours to obtain the neighbour block information
     and assign values to data members in the grid block connectivity data structure. */
    
    Find_Neighbours(Input);
    
}


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
        
        /* Add all the blocks in the database in order to compute the
         connectivity between blocks. */
        
        BlkC::BlockConnectivity blkConn(NBlk, 
                                        BlkC::CellCenter,
                                        Input.Nghost, 6);
        
        BlkC::VertexPack vp;
        int id,jd,kd;
        int id_n, jd_n, kd_n;
        Vector3D dX_i, dX_j, dX_k;
        
        for (int iblk = 0 ; iblk < NBlk ; ++iblk) {
            // A vertex pack (8 vertices and 24 vectors) for each block
            for (int i = 0; i != 2; ++i ) {
                for (int j = 0; j != 2; ++j ) {
                    for (int k = 0; k != 2; ++k ) {
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
                        
                        // Set coordinate of vertex
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
                    } /* endfor */
                } /* endfor */
            } /* endfor */
            
            blkConn.add_block(iblk, 
                              Input.NCells_Idir, 
                              Input.NCells_Jdir, 
                              Input.NCells_Kdir, 
                              vp);
        } /* endfor */
        
        /* Get the number of blocks on each element (total 26 elements) for each block
         and the neighbour information. */
        
        for (int iblk = 0 ; iblk <NBlk ; ++iblk) {
            /* Top neigbour */ 
            // Number of blocks on the top face
            Connectivity[iblk].num_neighT = blkConn.num_neighbour_block(iblk, 
                                                                        BlkC::IAll, 
                                                                        BlkC::JAll, 
                                                                        BlkC::KMax);
            
            // Information of neighbour block of the top face
            if (Connectivity[iblk].num_neighT > 0) {
                Connectivity[iblk].neighT_info.set_block_orientation_info(blkConn, 
                                                                          iblk,
                                                                          0,
                                                                          BlkC::IAll, 
                                                                          BlkC::JAll, 
                                                                          BlkC::KMax, 
                                                                          Connectivity[iblk].neighT);
            } /* endif */
            if (blkConn.at_domain_extent(iblk, BlkC::IAll, BlkC::JAll, BlkC::KMax)) {
                Connectivity[iblk].be.on_grid_boundary[BE::T] = 1;
            } /* endif */
            
            /* Bottom neigbour */ 
            // Number of blocks on the bottom face
            Connectivity[iblk].num_neighB = blkConn.num_neighbour_block(iblk,  
                                                                        BlkC::IAll, 
                                                                        BlkC::JAll, 
                                                                        BlkC::KMin);
            // Information of neighbour block of the bottom face
            if (Connectivity[iblk].num_neighB > 0) {
                Connectivity[iblk].neighB_info.set_block_orientation_info(blkConn, 
                                                                          iblk,
                                                                          0,
                                                                          BlkC::IAll, 
                                                                          BlkC::JAll, 
                                                                          BlkC::KMin, 
                                                                          Connectivity[iblk].neighB);
            } /* endif */
            if (blkConn.at_domain_extent(iblk, BlkC::IAll, BlkC::JAll, BlkC::KMin)) {
                Connectivity[iblk].be.on_grid_boundary[BE::B] = 1;
            } /* endif */
            
            /* North neigbour */ 
            // Number of blocks on the north face
            Connectivity[iblk].num_neighN = blkConn.num_neighbour_block(iblk,
                                                                        BlkC::IAll, 
                                                                        BlkC::JMax, 
                                                                        BlkC::KAll);
            // Information of neighbour block of the north face
            if (Connectivity[iblk].num_neighN > 0) {
                Connectivity[iblk].neighN_info.set_block_orientation_info(blkConn, 
                                                                          iblk,
                                                                          0,  
                                                                          BlkC::IAll, 
                                                                          BlkC::JMax, 
                                                                          BlkC::KAll, 
                                                                          Connectivity[iblk].neighN);
            } /* endif */
            if (blkConn.at_domain_extent(iblk, BlkC::IAll, BlkC::JMax, BlkC::KAll)) {
                Connectivity[iblk].be.on_grid_boundary[BE::N] = 1;
            } /* endif */
            
            /* South neigbour */ 
            // Number of blocks on the south face
            Connectivity[iblk].num_neighS = blkConn.num_neighbour_block(iblk,
                                                                        BlkC::IAll, 
                                                                        BlkC::JMin, 
                                                                        BlkC::KAll);
            // Information of neighbour block of the south face
            if (Connectivity[iblk].num_neighS > 0) {
                Connectivity[iblk].neighS_info.set_block_orientation_info(blkConn, 
                                                                          iblk,
                                                                          0,   
                                                                          BlkC::IAll, 
                                                                          BlkC::JMin, 
                                                                          BlkC::KAll, 
                                                                          Connectivity[iblk].neighS);
            } /* endif */
            if (blkConn.at_domain_extent(iblk, BlkC::IAll, BlkC::JMin, BlkC::KAll)) {
                Connectivity[iblk].be.on_grid_boundary[BE::S] = 1;
            } /* endif */
            
            /* East neigbour */ 
            // Number of blocks on the east face
            Connectivity[iblk].num_neighE = blkConn.num_neighbour_block(iblk, 
                                                                        BlkC::IMax, 
                                                                        BlkC::JAll, 
                                                                        BlkC::KAll);
            // Information of neighbour block of the east face
            if (Connectivity[iblk].num_neighE > 0) {
                Connectivity[iblk].neighE_info.set_block_orientation_info(blkConn, 
                                                                          iblk,
                                                                          0, 
                                                                          BlkC::IMax, 
                                                                          BlkC::JAll, 
                                                                          BlkC::KAll, 
                                                                          Connectivity[iblk].neighE);
            } /* endif */
            if (blkConn.at_domain_extent(iblk, BlkC::IMax, BlkC::JAll, BlkC::KAll)) {
                Connectivity[iblk].be.on_grid_boundary[BE::E] = 1;
            } /* endif */
            
            /* West neigbour */ 
            // Number of blocks on the west face
            Connectivity[iblk].num_neighW = blkConn.num_neighbour_block(iblk,  
                                                                        BlkC::IMin, 
                                                                        BlkC::JAll, 
                                                                        BlkC::KAll);
            // Information of neighbour block of the west face
            if (Connectivity[iblk].num_neighW > 0) {
                Connectivity[iblk].neighW_info.set_block_orientation_info(blkConn, 
                                                                          iblk,
                                                                          0,  
                                                                          BlkC::IMin, 
                                                                          BlkC::JAll, 
                                                                          BlkC::KAll, 
                                                                          Connectivity[iblk].neighW);
            } /* endif */
            if (blkConn.at_domain_extent(iblk, BlkC::IMin, BlkC::JAll, BlkC::KAll)) {
                Connectivity[iblk].be.on_grid_boundary[BE::W] = 1;
            } /* endif */
            
            /* Top-North neighbour */ 
            // Number of blocks on the top-north edge
            Connectivity[iblk].num_neighTN = blkConn.num_neighbour_block(iblk, 
                                                                         BlkC::IAll, 
                                                                         BlkC::JMax, 
                                                                         BlkC::KMax);
            // Information of neighbour block of the top-north edge
            for (int i_neigh = 0 ; i_neigh < Connectivity[iblk].num_neighTN ; ++i_neigh) {
                Connectivity[iblk].neighTN_info[i_neigh].set_block_orientation_info(blkConn, 
                                                                                    iblk,
                                                                                    i_neigh, 
                                                                                    BlkC::IAll, 
                                                                                    BlkC::JMax, 
                                                                                    BlkC::KMax, 
                                                                                    Connectivity[iblk].neighTN[i_neigh]);
            } /* endfor */
            if (blkConn.at_domain_extent(iblk,BlkC::IAll, BlkC::JMax, BlkC::KMax)) {
                Connectivity[iblk].be.on_grid_boundary[BE::TN] = 1;
            } /* endif */
            
            /* Top-South neighbour */ 
            // Number of blocks on the top-south edge
            Connectivity[iblk].num_neighTS = blkConn.num_neighbour_block(iblk, 
                                                                         BlkC::IAll, 
                                                                         BlkC::JMin, 
                                                                         BlkC::KMax);
            // Information of neighbour block of the top-south edge
            for (int i_neigh = 0 ; i_neigh < Connectivity[iblk].num_neighTS ; ++i_neigh) {
                Connectivity[iblk].neighTS_info[i_neigh].set_block_orientation_info(blkConn, 
                                                                                    iblk,
                                                                                    i_neigh,
                                                                                    BlkC::IAll, 
                                                                                    BlkC::JMin, 
                                                                                    BlkC::KMax, 
                                                                                    Connectivity[iblk].neighTS[i_neigh]);
            } /* endfor */
            if (blkConn.at_domain_extent(iblk, BlkC::IAll, BlkC::JMin, BlkC::KMax)) {
                Connectivity[iblk].be.on_grid_boundary[BE::TS] = 1;
            } /* endif */
            
            /* Top-West neighbour */ 
            // Number of blocks on the top-west edge
            Connectivity[iblk].num_neighTW = blkConn.num_neighbour_block(iblk,  
                                                                         BlkC::IMin, 
                                                                         BlkC::JAll, 
                                                                         BlkC::KMax);
            // Information of neighbour block of the top-west edge
            for (int i_neigh = 0 ; i_neigh < Connectivity[iblk].num_neighTW ; ++i_neigh) {
                Connectivity[iblk].neighTW_info[i_neigh].set_block_orientation_info(blkConn, 
                                                                                    iblk,
                                                                                    i_neigh,  
                                                                                    BlkC::IMin, 
                                                                                    BlkC::JAll, 
                                                                                    BlkC::KMax, 
                                                                                    Connectivity[iblk].neighTW[i_neigh]);
            } /* endfor */
            if (blkConn.at_domain_extent(iblk, BlkC::IMin, BlkC::JAll, BlkC::KMax)) {
                Connectivity[iblk].be.on_grid_boundary[BE::TW] = 1;
            } /* endif */
            
            /* Top-East neighbour */ 
            // Number of blocks on the top-east edge
            Connectivity[iblk].num_neighTE = blkConn.num_neighbour_block(iblk,  
                                                                         BlkC::IMax, 
                                                                         BlkC::JAll, 
                                                                         BlkC::KMax);
            // Information of neighbour block of the top-east edge
            for (int i_neigh = 0 ; i_neigh < Connectivity[iblk].num_neighTE ; ++i_neigh) {
                Connectivity[iblk].neighTE_info[i_neigh].set_block_orientation_info(blkConn, 
                                                                                    iblk,
                                                                                    i_neigh,  
                                                                                    BlkC::IMax, 
                                                                                    BlkC::JAll, 
                                                                                    BlkC::KMax, 
                                                                                    Connectivity[iblk].neighTE[i_neigh]);
            } /* endfor */
            if (blkConn.at_domain_extent(iblk, BlkC::IMax, BlkC::JAll, BlkC::KMax)) {
                Connectivity[iblk].be.on_grid_boundary[BE::TE] = 1;
            } /* endif */
            
            /* Bottom-North neighbour */ 
            // Number of blocks on the bottom-north edge
            Connectivity[iblk].num_neighBN = blkConn.num_neighbour_block(iblk,
                                                                         BlkC::IAll, 
                                                                         BlkC::JMax, 
                                                                         BlkC::KMin);
            // Information of neighbour block of the bottom-north edge 
            for (int i_neigh = 0 ; i_neigh < Connectivity[iblk].num_neighBN ; ++i_neigh) {
                Connectivity[iblk].neighBN_info[i_neigh].set_block_orientation_info(blkConn, 
                                                                                    iblk,
                                                                                    i_neigh,  
                                                                                    BlkC::IAll, 
                                                                                    BlkC::JMax, 
                                                                                    BlkC::KMin, 
                                                                                    Connectivity[iblk].neighBN[i_neigh]);
            } /* endfor */
            if (blkConn.at_domain_extent(iblk,  BlkC::IAll, BlkC::JMax, BlkC::KMin)) {
                Connectivity[iblk].be.on_grid_boundary[BE::BN] = 1;
            } /* endif */
            
            /* Bottom-South neighbour */ 
            // Number of blocks on the bottom-south edge
            Connectivity[iblk].num_neighBS = blkConn.num_neighbour_block(iblk,
                                                                         BlkC::IAll, 
                                                                         BlkC::JMin, 
                                                                         BlkC::KMin);
            // Information of neighbour block of the bottom-south edge 
            for (int i_neigh = 0 ; i_neigh < Connectivity[iblk].num_neighBS ; ++i_neigh) {
                Connectivity[iblk].neighBS_info[i_neigh].set_block_orientation_info(blkConn, 
                                                                                    iblk,
                                                                                    i_neigh,  
                                                                                    BlkC::IAll, 
                                                                                    BlkC::JMin, 
                                                                                    BlkC::KMin, 
                                                                                    Connectivity[iblk].neighBS[i_neigh]);
            } /* endfor */
            if (blkConn.at_domain_extent(iblk, BlkC::IAll, BlkC::JMin, BlkC::KMin)) {
                Connectivity[iblk].be.on_grid_boundary[BE::BS] = 1;
            } /* endif */
            
            /* Bottom-West neighbour */ 
            // Number of blocks on the bottom-west edge
            Connectivity[iblk].num_neighBW = blkConn.num_neighbour_block(iblk,
                                                                         BlkC::IMin, 
                                                                         BlkC::JAll, 
                                                                         BlkC::KMin);
            // Information of neighbour block of the bottom-west edge  
            for (int i_neigh = 0 ; i_neigh < Connectivity[iblk].num_neighBW ; ++i_neigh) {
                Connectivity[iblk].neighBW_info[i_neigh].set_block_orientation_info(blkConn, 
                                                                                    iblk,
                                                                                    i_neigh,   
                                                                                    BlkC::IMin, 
                                                                                    BlkC::JAll, 
                                                                                    BlkC::KMin, 
                                                                                    Connectivity[iblk].neighBW[i_neigh]);
            } /* endfor */
            if (blkConn.at_domain_extent(iblk,  BlkC::IMin, BlkC::JAll, BlkC::KMin)) {
                Connectivity[iblk].be.on_grid_boundary[BE::BW] = 1;
            } /* endif */
            
            /* Bottom-East neighbour */ 
            // Number of blocks on the bottom-east edge
            Connectivity[iblk].num_neighBE = blkConn.num_neighbour_block(iblk,  
                                                                         BlkC::IMax, 
                                                                         BlkC::JAll, 
                                                                         BlkC::KMin);
            // Information of neighbour block of the bottomeast edge
            for (int i_neigh = 0 ; i_neigh < Connectivity[iblk].num_neighBE ; ++i_neigh) {
                Connectivity[iblk].neighBE_info[i_neigh].set_block_orientation_info(blkConn, 
                                                                                    iblk,
                                                                                    i_neigh,  
                                                                                    BlkC::IMax, 
                                                                                    BlkC::JAll, 
                                                                                    BlkC::KMin, 
                                                                                    Connectivity[iblk].neighBE[i_neigh]);
            } /* endfor */
            if (blkConn.at_domain_extent(iblk,  BlkC::IMax, BlkC::JAll, BlkC::KMin)) {
                Connectivity[iblk].be.on_grid_boundary[BE::BE] = 1;
            } /* endif */
            
            /* North-West neighbour */ 
            // Number of blocks on the north-west edge         
            Connectivity[iblk].num_neighNW = blkConn.num_neighbour_block(iblk, 
                                                                         BlkC::IMin, 
                                                                         BlkC::JMax, 
                                                                         BlkC::KAll);
            // Information of neighbour block of the north-west edge
            for (int i_neigh = 0 ; i_neigh < Connectivity[iblk].num_neighNW ; ++i_neigh) {
                Connectivity[iblk].neighNW_info[i_neigh].set_block_orientation_info(blkConn, 
                                                                                    iblk,
                                                                                    i_neigh, 
                                                                                    BlkC::IMin, 
                                                                                    BlkC::JMax, 
                                                                                    BlkC::KAll, 
                                                                                    Connectivity[iblk].neighNW[i_neigh]);
            } /* endfor */
            if (blkConn.at_domain_extent(iblk, BlkC::IMin, BlkC::JMax, BlkC::KAll)) {
                Connectivity[iblk].be.on_grid_boundary[BE::NW] = 1;
            } /* endif */
            
            /* North-East neighbour */ 
            // Number of blocks on the north-east edge 
            Connectivity[iblk].num_neighNE = blkConn.num_neighbour_block(iblk,  
                                                                         BlkC::IMax, 
                                                                         BlkC::JMax, 
                                                                         BlkC::KAll);
            // Information of neighbour block of the north-east edge
            for (int i_neigh = 0 ; i_neigh < Connectivity[iblk].num_neighNE ; ++i_neigh) {
                Connectivity[iblk].neighNE_info[i_neigh].set_block_orientation_info(blkConn, 
                                                                                    iblk,
                                                                                    i_neigh,
                                                                                    BlkC::IMax, 
                                                                                    BlkC::JMax, 
                                                                                    BlkC::KAll, 
                                                                                    Connectivity[iblk].neighNE[i_neigh]);
            } /* endfor */
            if (blkConn.at_domain_extent(iblk, BlkC::IMax, BlkC::JMax, BlkC::KAll)) {
                Connectivity[iblk].be.on_grid_boundary[BE::NE] = 1;
            } /* endif */
            
            
            /* South-East neighbour */ 
            // Number of blocks on the south-east edge 
            Connectivity[iblk].num_neighSE = blkConn.num_neighbour_block(iblk,  
                                                                         BlkC::IMax, 
                                                                         BlkC::JMin, 
                                                                         BlkC::KAll);
            // Information of neighbour block of the south-east edge 
            for (int i_neigh = 0 ; i_neigh < Connectivity[iblk].num_neighSE ; ++i_neigh) {
                Connectivity[iblk].neighSE_info[i_neigh].set_block_orientation_info(blkConn, 
                                                                                    iblk,
                                                                                    i_neigh,  
                                                                                    BlkC::IMax, 
                                                                                    BlkC::JMin, 
                                                                                    BlkC::KAll, 
                                                                                    Connectivity[iblk].neighSE[i_neigh]);
            } /* endfor */
            if (blkConn.at_domain_extent(iblk, BlkC::IMax, BlkC::JMin, BlkC::KAll)) {
                Connectivity[iblk].be.on_grid_boundary[BE::SE] = 1;
            } /* endif */
            
            /* South-West neighbour */ 
            // Number of blocks on the south-west edge                
            Connectivity[iblk].num_neighSW = blkConn.num_neighbour_block(iblk,
                                                                         BlkC::IMin,
                                                                         BlkC::JMin,
                                                                         BlkC::KAll);
            // Information of neighbour block of the south-west edge 
            for (int i_neigh = 0 ; i_neigh < Connectivity[iblk].num_neighSW ; ++i_neigh) {
                Connectivity[iblk].neighSW_info[i_neigh].set_block_orientation_info(blkConn,
                                                                                    iblk,
                                                                                    i_neigh,
                                                                                    BlkC::IMin,
                                                                                    BlkC::JMin,
                                                                                    BlkC::KAll, 
                                                                                    Connectivity[iblk].neighSW[i_neigh]);
            } /* endfor */
            if (blkConn.at_domain_extent(iblk,  BlkC::IMin, BlkC::JMin, BlkC::KAll)) {
                Connectivity[iblk].be.on_grid_boundary[BE::SW] = 1;
            } /* endif */
            
            /* Top-North-West neighbour */ 
            // Number of blocks on the top-north-west vertex
            Connectivity[iblk].num_neighTNW = blkConn.num_neighbour_block(iblk,
                                                                          BlkC::IMin,
                                                                          BlkC::JMax,
                                                                          BlkC::KMax);
            // Information of neighbour block of the top-north-west vertex
            for (int i_neigh = 0 ; i_neigh < Connectivity[iblk].num_neighTNW ; ++i_neigh) {
                Connectivity[iblk].neighTNW_info[i_neigh].set_block_orientation_info(blkConn,
                                                                                     iblk,
                                                                                     i_neigh,
                                                                                     BlkC::IMin,
                                                                                     BlkC::JMax,
                                                                                     BlkC::KMax, 
                                                                                     Connectivity[iblk].neighTNW[i_neigh]);
            } /* endfor */
            if (blkConn.at_domain_extent(iblk,  BlkC::IMin, BlkC::JMax, BlkC::KMax)) {
                Connectivity[iblk].be.on_grid_boundary[BE::TNW] = 1;
            } /* endif */
            
            /* Top-South-West neighbour */
            // Number of blocks on the top-south-west vertex
            Connectivity[iblk].num_neighTSW = blkConn.num_neighbour_block(iblk, 
                                                                          BlkC::IMin, 
                                                                          BlkC::JMin, 
                                                                          BlkC::KMax);
            // Information of neighbour block of the top-south-west vertex
            for (int i_neigh = 0 ; i_neigh < Connectivity[iblk].num_neighTSW ; ++i_neigh) {
                Connectivity[iblk].neighTSW_info[i_neigh].set_block_orientation_info(blkConn, 
                                                                                     iblk,
                                                                                     i_neigh,
                                                                                     BlkC::IMin,
                                                                                     BlkC::JMin,
                                                                                     BlkC::KMax, 
                                                                                     Connectivity[iblk].neighTSW[i_neigh]);
            } /* endfor */
            if (blkConn.at_domain_extent(iblk,  BlkC::IMin, BlkC::JMin, BlkC::KMax)) {
                Connectivity[iblk].be.on_grid_boundary[BE::TSW] = 1;
            } /* endif */
            
            /* Top-North-East neighbour */
            // Number of blocks on the top-north-east vertex
            Connectivity[iblk].num_neighTNE = blkConn.num_neighbour_block(iblk,  
                                                                          BlkC::IMax, 
                                                                          BlkC::JMax, 
                                                                          BlkC::KMax);
            // Information of neighbour block of the top-north-east vertex
            for (int i_neigh = 0 ; i_neigh < Connectivity[iblk].num_neighTNE ; ++i_neigh) {
                Connectivity[iblk].neighTNE_info[i_neigh].set_block_orientation_info(blkConn, 
                                                                                     iblk,
                                                                                     i_neigh,  
                                                                                     BlkC::IMax, 
                                                                                     BlkC::JMax, 
                                                                                     BlkC::KMax, 
                                                                                     Connectivity[iblk].neighTNE[i_neigh]);
            } /* endfor */
            if (blkConn.at_domain_extent(iblk,  BlkC::IMax, BlkC::JMax, BlkC::KMax)) {
                Connectivity[iblk].be.on_grid_boundary[BE::TNE] = 1;
            } /* endif */
            
            /* Top-South-East neighbour */
            // Number of blocks on the top-south-east vertex
            Connectivity[iblk].num_neighTSE = blkConn.num_neighbour_block(iblk, 
                                                                          BlkC::IMax, 
                                                                          BlkC::JMin, 
                                                                          BlkC::KMax);
            // Information of neighbour block of the top-south-east vertex
            for (int i_neigh = 0 ; i_neigh < Connectivity[iblk].num_neighTSE ; ++i_neigh) {
                Connectivity[iblk].neighTSE_info[i_neigh].set_block_orientation_info(blkConn,
                                                                                     iblk,
                                                                                     i_neigh,
                                                                                     BlkC::IMax,
                                                                                     BlkC::JMin,
                                                                                     BlkC::KMax, 
                                                                                     Connectivity[iblk].neighTSE[i_neigh]);
            } /* endfor */
            if (blkConn.at_domain_extent(iblk,   BlkC::IMax, BlkC::JMin, BlkC::KMax)) {
                Connectivity[iblk].be.on_grid_boundary[BE::TSE] = 1;
            } /* endif */
            
            /* Bottom-North-West neighbour */
            // Number of blocks on the bottom-north-west vertex
            Connectivity[iblk].num_neighBNW = blkConn.num_neighbour_block(iblk,
                                                                          BlkC::IMin,
                                                                          BlkC::JMax,
                                                                          BlkC::KMin);
            // Information of neighbour block of the bottom-north-west vertex
            for (int i_neigh = 0 ; i_neigh < Connectivity[iblk].num_neighBNW ; ++i_neigh) {
                Connectivity[iblk].neighBNW_info[i_neigh].set_block_orientation_info(blkConn,
                                                                                     iblk,
                                                                                     i_neigh,
                                                                                     BlkC::IMin,
                                                                                     BlkC::JMax,
                                                                                     BlkC::KMin, 
                                                                                     Connectivity[iblk].neighBNW[i_neigh]);
            } /* endfor */
            if (blkConn.at_domain_extent(iblk,   BlkC::IMin, BlkC::JMax, BlkC::KMin)) {
                Connectivity[iblk].be.on_grid_boundary[BE::BNW] = 1;
            } /* endif */
            
            /* Bottom-South-West neighbour */
            // Number of blocks on the bottom-south-west vertex
            Connectivity[iblk].num_neighBSW = blkConn.num_neighbour_block(iblk,
                                                                          BlkC::IMin, 
                                                                          BlkC::JMin, 
                                                                          BlkC::KMin);
            // Information of neighbour block of the bottom-south-west vertex
            for (int i_neigh = 0 ; i_neigh < Connectivity[iblk].num_neighBSW ; ++i_neigh) {
                Connectivity[iblk].neighBSW_info[i_neigh].set_block_orientation_info(blkConn, 
                                                                                     iblk,
                                                                                     i_neigh,  
                                                                                     BlkC::IMin, 
                                                                                     BlkC::JMin, 
                                                                                     BlkC::KMin, 
                                                                                     Connectivity[iblk].neighBSW[i_neigh]);
            } /* endfor */
            if (blkConn.at_domain_extent(iblk,   BlkC::IMin, BlkC::JMin, BlkC::KMin)) {
                Connectivity[iblk].be.on_grid_boundary[BE::BSW] = 1;
            } /* endif */
            
            /* Bottom-North-East neighbour */
            // Number of blocks on the bottom-north-east vertex
            Connectivity[iblk].num_neighBNE = blkConn.num_neighbour_block(iblk,  
                                                                          BlkC::IMax, 
                                                                          BlkC::JMax, 
                                                                          BlkC::KMin);
            // Information of neighbour block of the bottom-north-east vertex
            for (int i_neigh = 0 ; i_neigh < Connectivity[iblk].num_neighBNE ; ++i_neigh) {
                Connectivity[iblk].neighBNE_info[i_neigh].set_block_orientation_info(blkConn, 
                                                                                     iblk,
                                                                                     i_neigh,  
                                                                                     BlkC::IMax, 
                                                                                     BlkC::JMax, 
                                                                                     BlkC::KMin, 
                                                                                     Connectivity[iblk].neighBNE[i_neigh]);
            } /* endfor */
            if (blkConn.at_domain_extent(iblk,   BlkC::IMax, BlkC::JMax, BlkC::KMin)) {
                Connectivity[iblk].be.on_grid_boundary[BE::BNE] = 1;
            } /* endif */
            
            /* Bottom-South-East neighbour */
            // Number of blocks on the bottom-south-east vertex
            Connectivity[iblk].num_neighBSE = blkConn.num_neighbour_block(iblk,  
                                                                          BlkC::IMax, 
                                                                          BlkC::JMin, 
                                                                          BlkC::KMin);
            // Information of neighbour block of the bottom-south-east vertex
            for (int i_neigh = 0 ; i_neigh < Connectivity[iblk].num_neighBSE ; ++i_neigh) {
                Connectivity[iblk].neighBSE_info[i_neigh].set_block_orientation_info(blkConn, 
                                                                                     iblk,
                                                                                     i_neigh, 
                                                                                     BlkC::IMax, 
                                                                                     BlkC::JMin, 
                                                                                     BlkC::KMin, 
                                                                                     Connectivity[iblk].neighBSE[i_neigh]);
            } /* endfor */
            if (blkConn.at_domain_extent(iblk,  BlkC::IMax, BlkC::JMin, BlkC::KMin)) {
                Connectivity[iblk].be.on_grid_boundary[BE::BSE] = 1;
            } /* endif */
            
            if (Connectivity[iblk].num_neighT > 1 ||
                Connectivity[iblk].num_neighB > 1 ||
                Connectivity[iblk].num_neighN > 1 || 
                Connectivity[iblk].num_neighS > 1 ||
                Connectivity[iblk].num_neighE > 1 || 
                Connectivity[iblk].num_neighW > 1 || 
                Connectivity[iblk].num_neighTN > GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS ||
                Connectivity[iblk].num_neighTS > GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS ||
                Connectivity[iblk].num_neighTE > GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS ||
                Connectivity[iblk].num_neighTW > GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS ||
                Connectivity[iblk].num_neighBN > GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS ||
                Connectivity[iblk].num_neighBS > GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS ||
                Connectivity[iblk].num_neighBE > GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS ||
                Connectivity[iblk].num_neighBW > GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS ||
                Connectivity[iblk].num_neighTNW > GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS ||
                Connectivity[iblk].num_neighTNE > GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS ||
                Connectivity[iblk].num_neighTSW > GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS ||
                Connectivity[iblk].num_neighTSE > GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS ||
                Connectivity[iblk].num_neighBNW > GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS ||
                Connectivity[iblk].num_neighBNE > GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS ||
                Connectivity[iblk].num_neighBSW > GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS ||
                Connectivity[iblk].num_neighBSE > GRID3D_HEXA_MULTI_BLOCK_MAX_NEIGHBOURS) {
                cerr<<"\n Number of neighbour blocks on some boudnary elements is out of bounds. \n";
                exit(1);
            } /* endif */
        } /* endfor */
        
    } /* endif */
    
}

void Grid3D_Hexa_Multi_Block_List::Update_Cells(void) {
    for (int n=0; n<NBlk; n++) {
        Grid_Blks[n].Update_Cells();
    }
}

void Grid3D_Hexa_Multi_Block_List::Disturb_Interior_Nodes(const int Number_of_Iterations) {
    srand48(1); // make sure every generation will be the same
    for (int n=0; n<NBlk; n++) {
        Grid_Blks[n].Disturb_Interior_Nodes(Number_of_Iterations);
    }
    Update_Cells();
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
        } else {
            Deallocate();
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
    
    // Create various multiblock multiblock grid depending on input parameters
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
            
        case GRID_BUMP_CHANNEL_FLOW :
            Create_Grid_Bump_Channel_Flow(Input);
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
    
    Grid2D_Pipe_XYplane = Grid_Tube_2D(Grid2D_Pipe_XYplane,
                                       numblk_idir_pipe,
                                       numblk_jdir_pipe,
                                       Input.Pipe_Radius,
                                       Input.NCells_Idir,
                                       Input.NCells_Jdir,
                                       Input.Nghost,
                                       STRETCHING_FCN_MAX_CLUSTERING,
                                       1.00001);
    
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
                                          BC_FIXED_PRESSURE,
                                          BC_FIXED_PRESSURE);
        } else {
            Grid_Blks[iBlk][0][0].Set_BCs(BC_NONE,
                                          BC_NONE,
                                          BC_NO_SLIP,
                                          BC_NONE,
                                          BC_FIXED_PRESSURE,
                                          BC_FIXED_PRESSURE);
        } /* endif */
        
    } /* endfor */
    
    /* Deallocate 2D grid. */
    
    Grid2D_Pipe_XYplane = Deallocate_Multi_Block_Grid(Grid2D_Pipe_XYplane,
                                                      numblk_idir_pipe,
                                                      numblk_jdir_pipe);
    
}

/**********************************************************************
 * Routine: Grid_Bump_Channel_Flow                                    *
 *                                                                    *
 * Generates a single block quadilateral mesh with clustering for     *
 * predicting supersonic flow around a cirucular cylinder blunt body. *
 *                                                                    *
 **********************************************************************/
void Grid3D_Hexa_Multi_Block::Create_Grid_Bump_Channel_Flow(Grid3D_Input_Parameters &Input) {
    
    Grid2D_Quad_Block **Grid2D_XYplane;
    int BC_top, BC_bottom;
    
    //Fixed at 4x2
    int Number_of_Blocks_Idir, Number_of_Blocks_Jdir;
    int Smooth_Bump = 0; //don't use smooth bump.
    
    /* Allocate required memory. */
    
    Input.NBlk_Idir = 4;
    Input.NBlk_Jdir = 2;
    Input.NBlk_Kdir = 1;
    
    Allocate(Input.NBlk_Idir, Input.NBlk_Jdir, Input.NBlk_Kdir);
    
    /* Creat 2D cross-section grids from which the 3D grid
     will be extruded. */
    
    Grid2D_XYplane = Grid_Bump_Channel_Flow(Grid2D_XYplane,
                                            Number_of_Blocks_Idir,
                                            Number_of_Blocks_Jdir,
                                            Smooth_Bump,
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
                
                Grid_Blks[iBlk][jBlk][kBlk].Extrude(Grid2D_XYplane[iBlk][jBlk],
                                                    Input.NCells_Kdir,
                                                    Input.Stretching_Type_Kdir,
                                                    Input.Stretching_Factor_Kdir,			       
                                                    (double(kBlk)/double(Input.NBlk_Kdir))*Input.Box_Length,
                                                    (double(kBlk+1)/double(Input.NBlk_Kdir))*Input.Box_Length);
                
                /* Assign top and bottom boundary conditions. */
                
                if (kBlk == Input.NBlk_Kdir-1) {
                    BC_top = BC_REFLECTION;//BC_CONSTANT_EXTRAPOLATION;
                } else {
                    BC_top = BC_NONE;
                } /* endif */
                if (kBlk == 0) {
                    BC_bottom =BC_REFLECTION;// BC_CONSTANT_EXTRAPOLATION;
                } else {
                    BC_bottom = BC_NONE;
                } /* endif */
                
                Grid_Blks[iBlk][jBlk][kBlk].Set_BCs_Zdir(BC_top, BC_bottom);
                
            } 
        } 
    } 
    
    /* Deallocate 2D grid. */
    Grid2D_XYplane = Deallocate_Multi_Block_Grid(Grid2D_XYplane,
                                                 Number_of_Blocks_Idir,
                                                 Number_of_Blocks_Jdir);
    
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
                if (Grid_Blks[iBlk][jBlk][kBlk].Allocated) {
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


Vector3D Grid3D_Hexa_Multi_Block_List::Delta_minimum(void) {
    double dx(1e10), dy(1e10), dz(1e10);
    Vector3D Delta;
    for (int  i = 0 ; i < NBlk ; ++i ) {
        if (Grid_Blks[i].Allocated) {
            Delta = Grid_Blks[i].Delta_minimum();
            if (Delta.x < dx)   dx = Delta.x;
            if (Delta.y < dy)   dy = Delta.y;
            if (Delta.z < dz)   dz = Delta.z;
        }  
    }
    return Vector3D(dx,dy,dz);
}

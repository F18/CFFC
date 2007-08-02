/* Grid3DHexaMultiBlock.cc:  Multi-block subroutines for 
                             3D hexahedral block grid class. */

/* Include 3D quadrilateral block grid type header file. */

#ifndef GRID3D_HEXA_MULTIBLOCK_INCLUDED
#include "Grid3DHexaMultiBlock.h"
#endif //_GRID3D_HEXA_MULTIBLOCK_INCLUDED


/*************************************************************************
 * Grid3D_Hexa_Block -- External subroutines for 3D array of grid blocks.*
 *************************************************************************/


/********************************************************
 * Routine: Write_Multi_Block_Grid                      *
 *                                                      *
 * Writes a 3D array of 3D hexahedral    multi-block    *
 * grids to the specified output stream for retrieval   *
 * and re-use purposes.                                 *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block::Write_Multi_Block_Grid(ostream &Out_File) {
   
   Out_File << NBI << " " 
            << NBJ << " "
            << NBK << "\n";
  
  for (int  k = 0 ; k < NBK ; ++k ) {
    for (int  j = 0 ; j < NBJ ; ++j ) {
      for (int  i = 0 ; i < NBI ; ++i ) {
	Grid_ptr[i][j][k]->Write_Hexa_Block((*Grid_ptr[i][j][k]), Out_File);
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
void Grid3D_Hexa_Multi_Block::Output_Tecplot(ostream &Out_File) {

    int block_number, i_output_title;

    block_number = 0;
    i_output_title = 1;
    
    for (int  k = 0 ; k <NBK ; ++k ) {
      for (int  j = 0; j <NBJ ; ++j ) {
	for (int i = 0 ; i <NBI ; ++i ) {
           if (Grid_ptr[i][j][k]->Node != NULL) {
              Grid_ptr[i][j][k]->Output_Tecplot(block_number,
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
void Grid3D_Hexa_Multi_Block::Output_Nodes_Tecplot(ostream &Out_File) {
  int  block_number = 0;
  int  i_output_title = 1;
  
  for (int  k = 0 ; k <NBK ; ++k ) {
    for (int  j = 0; j <NBJ ; ++j ) {
      for (int i = 0 ; i <NBI ; ++i ) {
	if (Grid_ptr[i][j][k]->Node != NULL) {
           Grid_ptr[i][j][k]->Output_Nodes_Tecplot(
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
void Grid3D_Hexa_Multi_Block::Output_Cells_Tecplot(ostream &Out_File) {

    int block_number, i_output_title;

    block_number = 0;
    i_output_title = 1;
  for (int  k = 0 ; k <NBK ; ++k ) {
    for (int  j = 0; j <NBJ ; ++j ) {
      for (int i = 0 ; i <NBI ; ++i ) {
	if (Grid_ptr[i][j][k]->Node != NULL) {
           Grid_ptr[i][j][k]->Output_Cells_Tecplot(
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
void Grid3D_Hexa_Multi_Block::Output_Gnuplot(ostream &Out_File) {
  
  int  block_number, i_output_title;

    block_number = 0;
    i_output_title = 1;
    
    for (int k = 0 ; k <NBK ; ++k ) {
      for (int j = 0; j <NBJ ; ++j ) {
	for(int i =0 ; i <NBI ; ++i ) {
	  if (Grid_ptr[i][j][k]->Node != NULL) {
	    Grid_ptr[i][j][k]->Output_Gnuplot(block_number,
			  i_output_title,
			  Out_File);
	   block_number = block_number + 1;
	   if (i_output_title) i_output_title = 0;
	  } /* endif */
	}  /* endfor */
      }  /* endfor */
    } /* endfor */
    
}

/********************************************************
 * Routine: Grid_Cube                                   *
 *                                                      *
 * Generates a 3D Cartesian mesh for a  cube            *
 *                                                      *
 *                                                      *
 * Usage:   Grid_Cube(Grid_ptr,                         *
 *                                        nblk_i,       *
 *                                        nblk_j,       *
 *                                        nblk_k,       *
 *                                        TEN,          *
 *                                        TEN,          *
 *                                        FIVE,         *
 *   	                                  100,          *
 *   	                                  100,          *
 *  	                                  50);          *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block::Grid_Cube
(const Grid3D_Input_Parameters  &IP){

   int iBlk, jBlk, kBlk, n_cells_i, n_cells_j, n_cells_k;
   double length, width, height;
   double x_orig, y_orig, z_orig;
   int   BCtype_top,   BCtype_bottom,   BCtype_north,
      BCtype_south,   BCtype_west,   BCtype_east;
   
 
   /* Allocate memory for grid block. */
   if (NBI < 0) NBI = 1;
   if (NBJ < 0) NBJ = 1;
   if (NBK < 0) NBK = 1;
  
     Grid3D_Hexa_Block *Grid_ptr_orig_blk ;
      Grid_ptr_orig_blk = new Grid3D_Hexa_Block(
         NBI*IP.ICells,  NBJ*IP.JCells , NBK*IP.KCells, IP.Nghost);
      // generate the total mesh in the original geometry
      
      BCtype_west =  BC_REFLECTION;
      BCtype_south = BC_REFLECTION;
      BCtype_bottom  = BC_REFLECTION;
      BCtype_east = BC_REFLECTION;
      BCtype_north = BC_REFLECTION;
      BCtype_top  = BC_REFLECTION;
      
      Grid_ptr_orig_blk->Create_Hexa_Block
         (  IP.Box_Length,
            IP.Box_Width,
            IP.Box_Height,
            0.,
            0.,
            0.,
            IP.I_stretching_factor,
            IP.J_stretching_factor,
            IP.K_stretching_factor,
            BCtype_top ,
            BCtype_bottom,
            BCtype_north,
            BCtype_south,
            BCtype_west,
            BCtype_east,
            NBI*IP.ICells,
            NBJ*IP.JCells,
            NBK*IP.KCells,
            IP.Nghost);
      // update cell coordinates
      Grid_ptr_orig_blk->Update_Cells();

  
   for ( kBlk = 0; kBlk <NBK; ++kBlk ) 
      for ( jBlk = 0; jBlk <NBJ; ++jBlk ) 
         for ( iBlk = 0; iBlk <NBI; ++iBlk ) {

            x_orig = iBlk*length;
            y_orig = jBlk*width ;
            z_orig = kBlk*height;

            if (iBlk==0 ){
               BCtype_west = BC_REFLECTION;
	    
                             
            }else {
               BCtype_west = BC_NONE;
            }
            if (jBlk==0 ){
               BCtype_south = BC_REFLECTION;
               
            }else {
               BCtype_south = BC_NONE;
            }
            if (kBlk==0 ){
               BCtype_bottom  = BC_REFLECTION;
               
	       
            }else {
               BCtype_bottom = BC_NONE;
            }
            if (iBlk== NBI-1 ){
               BCtype_east = BC_REFLECTION;
               
            }else {
               BCtype_east = BC_NONE;
            }
            if (jBlk== NBJ -1 ){
               BCtype_north = BC_REFLECTION;
               
 	    }else {
               BCtype_north = BC_NONE;
            }
            if (kBlk== NBK-1 ){
               BCtype_top  = BC_REFLECTION;
               
            }else {
               BCtype_top = BC_NONE;
            }

            /* copy the cells and nodes information  
               of the hexa block to the subblocks. */
            
            for (int k =  Grid_ptr[iBlk][jBlk][kBlk]->KNl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                 k <=  Grid_ptr[iBlk][jBlk][kBlk]->KNu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++k) 
               for (int j =  Grid_ptr[iBlk][jBlk][kBlk]->JNl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                    j <=  Grid_ptr[iBlk][jBlk][kBlk]->JNu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++j) 
                  for (int i =  Grid_ptr[iBlk][jBlk][kBlk]->INl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                       i <=  Grid_ptr[iBlk][jBlk][kBlk]->INu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++i) {
                     
                     Grid_ptr[iBlk][jBlk][kBlk]->Node[i][j][k].X.x = 
                        Grid_ptr_orig_blk->Node[iBlk*IP.ICells + i][jBlk*IP.JCells+ j][kBlk*IP.KCells +k].X.x;
                     
                     Grid_ptr[iBlk][jBlk][kBlk]->Node[i][j][k].X.y =  
                        Grid_ptr_orig_blk->Node[iBlk*IP.ICells + i][jBlk*IP.JCells+ j][kBlk*IP.KCells +k].X.y;
                     
                     Grid_ptr[iBlk][jBlk][kBlk]->Node[i][j][k].X.z =  
                        Grid_ptr_orig_blk->Node[iBlk*IP.ICells + i][jBlk*IP.JCells+ j][kBlk*IP.KCells +k].X.z;
                     
                  }
  
            Grid_ptr[iBlk][jBlk][kBlk]->Update_Cells();
               
            //store the boundary conditions
            for (int  j =  Grid_ptr[iBlk][jBlk][kBlk]->JCl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                 j <=  Grid_ptr[iBlk][jBlk][kBlk]->JCu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++j) 
               for (int  i =  Grid_ptr[iBlk][jBlk][kBlk]->ICl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                    i <=  Grid_ptr[iBlk][jBlk][kBlk]->ICu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++i) {
                  Grid_ptr[iBlk][jBlk][kBlk]->BCtypeT[i][j] =  BCtype_top;
                  Grid_ptr[iBlk][jBlk][kBlk]->BCtypeB[i][j] =  BCtype_bottom;
            
               }
            for ( int k =  Grid_ptr[iBlk][jBlk][kBlk]->KCl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                  k <=  Grid_ptr[iBlk][jBlk][kBlk]->KCu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++k) 
               for ( int j =  Grid_ptr[iBlk][jBlk][kBlk]->JCl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                     j <=  Grid_ptr[iBlk][jBlk][kBlk]->JCu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++j){
                  Grid_ptr[iBlk][jBlk][kBlk]->BCtypeW[j][k] =  BCtype_west;
                  Grid_ptr[iBlk][jBlk][kBlk]->BCtypeE[j][k] =  BCtype_east;
                     
               }
            for (int  k =  Grid_ptr[iBlk][jBlk][kBlk]->KCl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                 k <=  Grid_ptr[iBlk][jBlk][kBlk]->KCu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++k) 
               for ( int i =  Grid_ptr[iBlk][jBlk][kBlk]->ICl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                     i <=  Grid_ptr[iBlk][jBlk][kBlk]->ICu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++i){
                  Grid_ptr[iBlk][jBlk][kBlk]->BCtypeN[i][k] =  BCtype_north;
                  Grid_ptr[iBlk][jBlk][kBlk]->BCtypeS[i][k] =  BCtype_south;
                  
               }
         }


   delete Grid_ptr_orig_blk;
            
         


      
}// cubic box grid

void Grid3D_Hexa_Multi_Block::Grid_Channel
(const Grid3D_Input_Parameters  &IP){

   int iBlk, jBlk, kBlk, n_cells_i, n_cells_j, n_cells_k;
   double x_orig, y_orig, z_orig;
   int   BCtype_top,   BCtype_bottom,   BCtype_north,
      BCtype_south,   BCtype_west,   BCtype_east;
   
   
   /* Allocate memory for grid block. */
   if (NBI < 0) NBI = 1;
   if (NBJ < 0) NBJ = 1;
   if (NBK < 0) NBK = 1;
   
   if(IP.geometry_index == 1){

      Grid3D_Hexa_Block *Grid_ptr_orig_blk ;
      Grid_ptr_orig_blk = new Grid3D_Hexa_Block(
         NBI*IP.ICells,  NBJ*IP.JCells , NBK*IP.KCells, IP.Nghost);
      // generate the total mesh in the original geometry
      
      BCtype_west =  BC_CHANNEL_INFLOW;
      BCtype_south =  BC_NO_SLIP;
      BCtype_bottom  = BC_PERIODIC;
      BCtype_east = BC_CHANNEL_OUTFLOW;
      BCtype_north =  BC_NO_SLIP;
      BCtype_top  = BC_PERIODIC;
      
      Grid_ptr_orig_blk->Create_Hexa_Block
         (  IP.Box_Length,
            IP.Box_Width,
            IP.Box_Height,
            0.0,
            - IP.Box_Width/2.0,
            - IP.Box_Height/2.0,          
            IP.I_stretching_factor,
            IP.J_stretching_factor,
            IP.K_stretching_factor,
            BCtype_top ,
            BCtype_bottom,
            BCtype_north,
            BCtype_south,
            BCtype_west,
            BCtype_east,
            NBI*IP.ICells,
            NBJ*IP.JCells,
            NBK*IP.KCells,
            IP.Nghost);
      // update cell coordinates
      Grid_ptr_orig_blk->Update_Cells();


      for ( kBlk = 0; kBlk <NBK; ++kBlk ) 
         for ( jBlk = 0; jBlk <NBJ; ++jBlk ) 
            for ( iBlk = 0; iBlk <NBI; ++iBlk ) {
                
               if (iBlk==0 ){
                  BCtype_west = BC_CHANNEL_INFLOW;
               }else {
                  BCtype_west = BC_NONE;
               }
               if (jBlk==0 ){
                  BCtype_south = BC_NO_SLIP; //BC_PERIODIC;
               }else {
                  BCtype_south = BC_NONE;
               }
               if (kBlk==0 ){
                  BCtype_bottom  = BC_PERIODIC;//BC_NO_SLIP;
               }else {
                  BCtype_bottom = BC_NONE;
               }
               if (iBlk== NBI-1 ){
                  BCtype_east = BC_CHANNEL_OUTFLOW;
               }else {
                  BCtype_east = BC_NONE;
               }
               if (jBlk== NBJ -1 ){
                  BCtype_north = BC_NO_SLIP; //BC_PERIODIC;
               }else {
                  BCtype_north = BC_NONE;
               }
               if (kBlk== NBK-1 ){
                  BCtype_top  = BC_PERIODIC;//BC_NO_SLIP;
               }else {
                  BCtype_top = BC_NONE;
               }

        /* copy the cells and nodes information  
               of the hexa block to the subblocks. */

               for (int k =  Grid_ptr[iBlk][jBlk][kBlk]->KNl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                    k <=  Grid_ptr[iBlk][jBlk][kBlk]->KNu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++k) 
                  for (int j =  Grid_ptr[iBlk][jBlk][kBlk]->JNl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                       j <=  Grid_ptr[iBlk][jBlk][kBlk]->JNu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++j) 
                     for (int i =  Grid_ptr[iBlk][jBlk][kBlk]->INl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                          i <=  Grid_ptr[iBlk][jBlk][kBlk]->INu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++i) {
                        
                        Grid_ptr[iBlk][jBlk][kBlk]->Node[i][j][k].X.x = 
                           Grid_ptr_orig_blk->Node[iBlk*IP.ICells + i][jBlk*IP.JCells+ j][kBlk*IP.KCells +k].X.x;
                        
                        Grid_ptr[iBlk][jBlk][kBlk]->Node[i][j][k].X.y =  
                           Grid_ptr_orig_blk->Node[iBlk*IP.ICells + i][jBlk*IP.JCells+ j][kBlk*IP.KCells +k].X.y;
                        
                        Grid_ptr[iBlk][jBlk][kBlk]->Node[i][j][k].X.z =  
                           Grid_ptr_orig_blk->Node[iBlk*IP.ICells + i][jBlk*IP.JCells+ j][kBlk*IP.KCells +k].X.z;
                        
                     }
  
               Grid_ptr[iBlk][jBlk][kBlk]->Update_Cells();
               
               //store the boundary conditions
               for (int  j =  Grid_ptr[iBlk][jBlk][kBlk]->JCl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                    j <=  Grid_ptr[iBlk][jBlk][kBlk]->JCu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++j) 
                  for (int  i =  Grid_ptr[iBlk][jBlk][kBlk]->ICl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                       i <=  Grid_ptr[iBlk][jBlk][kBlk]->ICu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++i) {
                     Grid_ptr[iBlk][jBlk][kBlk]->BCtypeT[i][j] =  BCtype_top;
                     Grid_ptr[iBlk][jBlk][kBlk]->BCtypeB[i][j] =  BCtype_bottom;
            
                  }
               for ( int k =  Grid_ptr[iBlk][jBlk][kBlk]->KCl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                     k <=  Grid_ptr[iBlk][jBlk][kBlk]->KCu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++k) 
                  for ( int j =  Grid_ptr[iBlk][jBlk][kBlk]->JCl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                        j <=  Grid_ptr[iBlk][jBlk][kBlk]->JCu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++j){
                     Grid_ptr[iBlk][jBlk][kBlk]->BCtypeW[j][k] =  BCtype_west;
                     Grid_ptr[iBlk][jBlk][kBlk]->BCtypeE[j][k] =  BCtype_east;
                     
                  }
               for (int  k =  Grid_ptr[iBlk][jBlk][kBlk]->KCl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                    k <=  Grid_ptr[iBlk][jBlk][kBlk]->KCu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++k) 
                  for ( int i =  Grid_ptr[iBlk][jBlk][kBlk]->ICl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                        i <=  Grid_ptr[iBlk][jBlk][kBlk]->ICu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++i){
                     Grid_ptr[iBlk][jBlk][kBlk]->BCtypeN[i][k] =  BCtype_north;
                     Grid_ptr[iBlk][jBlk][kBlk]->BCtypeS[i][k] =  BCtype_south;
                     
                  }
            }


      delete Grid_ptr_orig_blk;
          
   }
   
   if(IP.geometry_index == 2){

      Grid3D_Hexa_Block *Grid_ptr_orig_blk ;


      Grid_ptr_orig_blk = new Grid3D_Hexa_Block(
         NBI*IP.ICells,  NBJ*IP.JCells , NBK*IP.KCells, IP.Nghost);
      // generate the total mesh in the original geometry
    
      BCtype_west = BC_PERIODIC;
      BCtype_south = BC_CHANNEL_OUTFLOW;
      BCtype_bottom  = BC_NO_SLIP;
      BCtype_east = BC_PERIODIC;
      BCtype_north = BC_CHANNEL_INFLOW ;
      BCtype_top  = BC_NO_SLIP;
      
      Grid_ptr_orig_blk->Create_Hexa_Block
         (  IP.Box_Length,
            IP.Box_Width,
            IP.Box_Height,
            - IP.Box_Length/2.0,
            0.0,
            - IP.Box_Height/2.0,          
            IP.I_stretching_factor,
            IP.J_stretching_factor,
            IP.K_stretching_factor,
            BCtype_top ,
            BCtype_bottom,
            BCtype_north,
            BCtype_south,
            BCtype_west,
            BCtype_east,
            NBI*IP.ICells,
            NBJ*IP.JCells,
            NBK*IP.KCells,
            IP.Nghost);
      // update cell coordinates
      Grid_ptr_orig_blk->Update_Cells();

      
      // copying the nodes coordinates from the parent block to 
      // subblocks.     
      for ( kBlk = 0; kBlk <NBK; ++kBlk ) 
         for ( jBlk = 0; jBlk <NBJ; ++jBlk ) 
            for ( iBlk = 0; iBlk <NBI; ++iBlk ) {
              
               if (iBlk==0 ){
                  BCtype_west = BC_PERIODIC;
               }else {
                  BCtype_west = BC_NONE;
               }
               if (jBlk==0 ){
                  BCtype_south = BC_CHANNEL_OUTFLOW;
               }else {
                  BCtype_south = BC_NONE;
               }
               if (kBlk==0 ){
                  BCtype_bottom  = BC_NO_SLIP;
               }else {
                  BCtype_bottom = BC_NONE;
               }
               if (iBlk== NBI-1 ){
                  BCtype_east = BC_PERIODIC;
               }else {
                  BCtype_east = BC_NONE;
               }
               if (jBlk== NBJ -1 ){
                  BCtype_north = BC_CHANNEL_INFLOW ;
               }else {
                  BCtype_north = BC_NONE;
               }
               if (kBlk== NBK-1 ){
                  BCtype_top  = BC_NO_SLIP;
               }else {
                  BCtype_top = BC_NONE;
               }

            /* copy the cells and nodes information  
               of the hexa block to the subblocks. */

               for (int k =  Grid_ptr[iBlk][jBlk][kBlk]->KNl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                    k <=  Grid_ptr[iBlk][jBlk][kBlk]->KNu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++k) 
                  for (int j =  Grid_ptr[iBlk][jBlk][kBlk]->JNl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                       j <=  Grid_ptr[iBlk][jBlk][kBlk]->JNu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++j) 
                     for (int i =  Grid_ptr[iBlk][jBlk][kBlk]->INl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                          i <=  Grid_ptr[iBlk][jBlk][kBlk]->INu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++i) {
                        
                        Grid_ptr[iBlk][jBlk][kBlk]->Node[i][j][k].X.x = 
                           Grid_ptr_orig_blk->Node[iBlk*IP.ICells + i][jBlk*IP.JCells+ j][kBlk*IP.KCells +k].X.x;
                        
                        Grid_ptr[iBlk][jBlk][kBlk]->Node[i][j][k].X.y =  
                           Grid_ptr_orig_blk->Node[iBlk*IP.ICells + i][jBlk*IP.JCells+ j][kBlk*IP.KCells +k].X.y;
                        
                        Grid_ptr[iBlk][jBlk][kBlk]->Node[i][j][k].X.z =  
                           Grid_ptr_orig_blk->Node[iBlk*IP.ICells + i][jBlk*IP.JCells+ j][kBlk*IP.KCells +k].X.z;
                        
                     }
  
               Grid_ptr[iBlk][jBlk][kBlk]->Update_Cells();
               
               //store the boundary conditions
               for (int  j =  Grid_ptr[iBlk][jBlk][kBlk]->JCl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                    j <=  Grid_ptr[iBlk][jBlk][kBlk]->JCu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++j) 
                  for (int  i =  Grid_ptr[iBlk][jBlk][kBlk]->ICl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                       i <=  Grid_ptr[iBlk][jBlk][kBlk]->ICu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++i) {
                     Grid_ptr[iBlk][jBlk][kBlk]->BCtypeT[i][j] =  BCtype_top;
                     Grid_ptr[iBlk][jBlk][kBlk]->BCtypeB[i][j] =  BCtype_bottom;
            
                  }
               for ( int k =  Grid_ptr[iBlk][jBlk][kBlk]->KCl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                     k <=  Grid_ptr[iBlk][jBlk][kBlk]->KCu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++k) 
                  for ( int j =  Grid_ptr[iBlk][jBlk][kBlk]->JCl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                        j <=  Grid_ptr[iBlk][jBlk][kBlk]->JCu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++j){
                     Grid_ptr[iBlk][jBlk][kBlk]->BCtypeW[j][k] =  BCtype_west;
                     Grid_ptr[iBlk][jBlk][kBlk]->BCtypeE[j][k] =  BCtype_east;
                     
                  }
               for (int  k =  Grid_ptr[iBlk][jBlk][kBlk]->KCl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                    k <=  Grid_ptr[iBlk][jBlk][kBlk]->KCu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++k) 
                  for ( int i =  Grid_ptr[iBlk][jBlk][kBlk]->ICl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                        i <=  Grid_ptr[iBlk][jBlk][kBlk]->ICu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++i){
                     Grid_ptr[iBlk][jBlk][kBlk]->BCtypeN[i][k] =  BCtype_north;
                     Grid_ptr[iBlk][jBlk][kBlk]->BCtypeS[i][k] =  BCtype_south;
                     
                  }
            }


      delete Grid_ptr_orig_blk;
      
   }
   

   if(IP.geometry_index == 3){

     Grid3D_Hexa_Block *Grid_ptr_orig_blk ;


     Grid_ptr_orig_blk = new Grid3D_Hexa_Block(
        NBI*IP.ICells,  NBJ*IP.JCells , NBK*IP.KCells, IP.Nghost);
     // generate the total mesh in the original geometry
    
      BCtype_west = BC_PERIODIC;
      BCtype_south = BC_NO_SLIP;
      BCtype_bottom  = BC_CHANNEL_INFLOW;
      BCtype_east = BC_PERIODIC;
      BCtype_north = BC_NO_SLIP;
      BCtype_top  = BC_CHANNEL_OUTFLOW;
      
      Grid_ptr_orig_blk->Create_Hexa_Block
         (  IP.Box_Length,
            IP.Box_Width,
            IP.Box_Height,
            - IP.Box_Length/2.0,
            - IP.Box_Width/2.0,
            0.0,
            IP.I_stretching_factor,
            IP.J_stretching_factor,
            IP.K_stretching_factor,
            BCtype_top ,
            BCtype_bottom,
            BCtype_north,
            BCtype_south,
            BCtype_west,
            BCtype_east,
            NBI*IP.ICells,
            NBJ*IP.JCells,
            NBK*IP.KCells,
            IP.Nghost);
      // update cell coordinates
      Grid_ptr_orig_blk->Update_Cells();


      
      for ( kBlk = 0; kBlk <NBK; ++kBlk ) 
         for ( jBlk = 0; jBlk <NBJ; ++jBlk ) 
            for ( iBlk = 0; iBlk <NBI; ++iBlk ) {
                 if (iBlk==0 ){
                  BCtype_west = BC_PERIODIC;
               }else {
                  BCtype_west = BC_NONE;
               }
               if (jBlk==0 ){
                  BCtype_south = BC_NO_SLIP;
               }else {
                  BCtype_south = BC_NONE;
               }
               if (kBlk==0 ){
                  BCtype_bottom  = BC_CHANNEL_INFLOW;
               }else {
                  BCtype_bottom = BC_NONE;
               }
               if (iBlk== NBI-1 ){
                  BCtype_east = BC_PERIODIC;
               }else {
                  BCtype_east = BC_NONE;
               }
               if (jBlk== NBJ -1 ){
                  BCtype_north = BC_NO_SLIP;
               }else {
                  BCtype_north = BC_NONE;
               }
               if (kBlk== NBK-1 ){
                  BCtype_top  = BC_CHANNEL_OUTFLOW;
               }else {
                  BCtype_top = BC_NONE;
               }
              

            /* copy the cells and nodes information  
               of the hexa block to the subblocks. */

               for (int k =  Grid_ptr[iBlk][jBlk][kBlk]->KNl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                    k <=  Grid_ptr[iBlk][jBlk][kBlk]->KNu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++k) 
                  for (int j =  Grid_ptr[iBlk][jBlk][kBlk]->JNl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                       j <=  Grid_ptr[iBlk][jBlk][kBlk]->JNu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++j) 
                     for (int i =  Grid_ptr[iBlk][jBlk][kBlk]->INl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                          i <=  Grid_ptr[iBlk][jBlk][kBlk]->INu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++i) {
                        
                        Grid_ptr[iBlk][jBlk][kBlk]->Node[i][j][k].X.x = 
                           Grid_ptr_orig_blk->Node[iBlk*IP.ICells + i][jBlk*IP.JCells+ j][kBlk*IP.KCells +k].X.x;
                        
                        Grid_ptr[iBlk][jBlk][kBlk]->Node[i][j][k].X.y =  
                           Grid_ptr_orig_blk->Node[iBlk*IP.ICells + i][jBlk*IP.JCells+ j][kBlk*IP.KCells +k].X.y;
                        
                        Grid_ptr[iBlk][jBlk][kBlk]->Node[i][j][k].X.z =  
                           Grid_ptr_orig_blk->Node[iBlk*IP.ICells + i][jBlk*IP.JCells+ j][kBlk*IP.KCells +k].X.z;
                        
                     }
  
               Grid_ptr[iBlk][jBlk][kBlk]->Update_Cells();
               
               //store the boundary conditions
               for (int  j =  Grid_ptr[iBlk][jBlk][kBlk]->JCl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                    j <=  Grid_ptr[iBlk][jBlk][kBlk]->JCu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++j) 
                  for (int  i =  Grid_ptr[iBlk][jBlk][kBlk]->ICl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                       i <=  Grid_ptr[iBlk][jBlk][kBlk]->ICu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++i) {
                     Grid_ptr[iBlk][jBlk][kBlk]->BCtypeT[i][j] =  BCtype_top;
                     Grid_ptr[iBlk][jBlk][kBlk]->BCtypeB[i][j] =  BCtype_bottom;
            
                  }
               for ( int k =  Grid_ptr[iBlk][jBlk][kBlk]->KCl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                     k <=  Grid_ptr[iBlk][jBlk][kBlk]->KCu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++k) 
                  for ( int j =  Grid_ptr[iBlk][jBlk][kBlk]->JCl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                        j <=  Grid_ptr[iBlk][jBlk][kBlk]->JCu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++j){
                     Grid_ptr[iBlk][jBlk][kBlk]->BCtypeW[j][k] =  BCtype_west;
                     Grid_ptr[iBlk][jBlk][kBlk]->BCtypeE[j][k] =  BCtype_east;
                     
                  }
               for (int  k =  Grid_ptr[iBlk][jBlk][kBlk]->KCl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                    k <=  Grid_ptr[iBlk][jBlk][kBlk]->KCu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++k) 
                  for ( int i =  Grid_ptr[iBlk][jBlk][kBlk]->ICl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                        i <=  Grid_ptr[iBlk][jBlk][kBlk]->ICu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++i){
                     Grid_ptr[iBlk][jBlk][kBlk]->BCtypeN[i][k] =  BCtype_north;
                     Grid_ptr[iBlk][jBlk][kBlk]->BCtypeS[i][k] =  BCtype_south;
                     
                  }
            }


      delete Grid_ptr_orig_blk;

   }
   

   
}// channel grid
void Grid3D_Hexa_Multi_Block::Grid_Couette
(const Grid3D_Input_Parameters  &IP){
   
   
   int iBlk, jBlk, kBlk, n_cells_i, n_cells_j, n_cells_k;
   double length, width, height;
   double x_orig, y_orig, z_orig;
   int   BCtype_top,   BCtype_bottom,   BCtype_north,
         BCtype_south,   BCtype_west,   BCtype_east;

   /* Allocate memory for grid block. */
   if (NBI < 0) NBI = 1;
   if (NBJ < 0) NBJ = 1;
   if (NBK < 0) NBK = 1;

  if(IP.geometry_index == 1){

      Grid3D_Hexa_Block *Grid_ptr_orig_blk ;
      Grid_ptr_orig_blk = new Grid3D_Hexa_Block(
         NBI*IP.ICells,  NBJ*IP.JCells , NBK*IP.KCells, IP.Nghost);
      // generate the total mesh in the original geometry
      
      BCtype_west =  BC_FIXED_PRESSURE;
      BCtype_south =   BC_PERIODIC;
      BCtype_bottom  =  BC_NO_SLIP;
      BCtype_east = BC_FIXED_PRESSURE;
      BCtype_north =  BC_PERIODIC;
      BCtype_top  = BC_MOVING_WALL;
      
      Grid_ptr_orig_blk->Create_Hexa_Block
         (  IP.Box_Length,
            IP.Box_Width,
            IP.Box_Height,
            0.0,
            - IP.Box_Width/2.0,
            - IP.Box_Height/2.0,          
            IP.I_stretching_factor,
            IP.J_stretching_factor,
            IP.K_stretching_factor,
            BCtype_top ,
            BCtype_bottom,
            BCtype_north,
            BCtype_south,
            BCtype_west,
            BCtype_east,
            NBI*IP.ICells,
            NBJ*IP.JCells,
            NBK*IP.KCells,
            IP.Nghost);
      // update cell coordinates
      Grid_ptr_orig_blk->Update_Cells();


      for ( kBlk = 0; kBlk <NBK; ++kBlk ) 
         for ( jBlk = 0; jBlk <NBJ; ++jBlk ) 
            for ( iBlk = 0; iBlk <NBI; ++iBlk ) {
         
               if (iBlk==0 ){
                  BCtype_west = BC_FIXED_PRESSURE;
               }else {
                  BCtype_west = BC_NONE;
               }
               if (jBlk==0 ){
                  BCtype_south = BC_PERIODIC;
               }else {
                  BCtype_south = BC_NONE;
               }
               if (kBlk==0 ){
                  BCtype_bottom  = BC_NO_SLIP;
               }else {
                  BCtype_bottom = BC_NONE;
               }
               if (iBlk== NBI-1 ){
                  BCtype_east = BC_FIXED_PRESSURE;
               }else {
                  BCtype_east = BC_NONE;
               }
               if (jBlk== NBJ -1 ){
                  BCtype_north = BC_PERIODIC;
               }else {
                  BCtype_north = BC_NONE;
               }
               if (kBlk== NBK-1 ){
                  BCtype_top  = BC_MOVING_WALL;
               }else {
                  BCtype_top = BC_NONE;
               } 
              

        /* copy the cells and nodes information  
               of the hexa block to the subblocks. */

               for (int k =  Grid_ptr[iBlk][jBlk][kBlk]->KNl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                    k <=  Grid_ptr[iBlk][jBlk][kBlk]->KNu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++k) 
                  for (int j =  Grid_ptr[iBlk][jBlk][kBlk]->JNl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                       j <=  Grid_ptr[iBlk][jBlk][kBlk]->JNu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++j) 
                     for (int i =  Grid_ptr[iBlk][jBlk][kBlk]->INl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                          i <=  Grid_ptr[iBlk][jBlk][kBlk]->INu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++i) {
                        
                        Grid_ptr[iBlk][jBlk][kBlk]->Node[i][j][k].X.x = 
                           Grid_ptr_orig_blk->Node[iBlk*IP.ICells + i][jBlk*IP.JCells+ j][kBlk*IP.KCells +k].X.x;
                        
                        Grid_ptr[iBlk][jBlk][kBlk]->Node[i][j][k].X.y =  
                           Grid_ptr_orig_blk->Node[iBlk*IP.ICells + i][jBlk*IP.JCells+ j][kBlk*IP.KCells +k].X.y;
                        
                        Grid_ptr[iBlk][jBlk][kBlk]->Node[i][j][k].X.z =  
                           Grid_ptr_orig_blk->Node[iBlk*IP.ICells + i][jBlk*IP.JCells+ j][kBlk*IP.KCells +k].X.z;
                        
                     }
  
               Grid_ptr[iBlk][jBlk][kBlk]->Update_Cells();
               
               //store the boundary conditions
               for (int  j =  Grid_ptr[iBlk][jBlk][kBlk]->JCl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                    j <=  Grid_ptr[iBlk][jBlk][kBlk]->JCu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++j) 
                  for (int  i =  Grid_ptr[iBlk][jBlk][kBlk]->ICl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                       i <=  Grid_ptr[iBlk][jBlk][kBlk]->ICu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++i) {
                     Grid_ptr[iBlk][jBlk][kBlk]->BCtypeT[i][j] =  BCtype_top;
                     Grid_ptr[iBlk][jBlk][kBlk]->BCtypeB[i][j] =  BCtype_bottom;
            
                  }
               for ( int k =  Grid_ptr[iBlk][jBlk][kBlk]->KCl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                     k <=  Grid_ptr[iBlk][jBlk][kBlk]->KCu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++k) 
                  for ( int j =  Grid_ptr[iBlk][jBlk][kBlk]->JCl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                        j <=  Grid_ptr[iBlk][jBlk][kBlk]->JCu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++j){
                     Grid_ptr[iBlk][jBlk][kBlk]->BCtypeW[j][k] =  BCtype_west;
                     Grid_ptr[iBlk][jBlk][kBlk]->BCtypeE[j][k] =  BCtype_east;
                     
                  }
               for (int  k =  Grid_ptr[iBlk][jBlk][kBlk]->KCl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                    k <=  Grid_ptr[iBlk][jBlk][kBlk]->KCu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++k) 
                  for ( int i =  Grid_ptr[iBlk][jBlk][kBlk]->ICl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                        i <=  Grid_ptr[iBlk][jBlk][kBlk]->ICu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++i){
                     Grid_ptr[iBlk][jBlk][kBlk]->BCtypeN[i][k] =  BCtype_north;
                     Grid_ptr[iBlk][jBlk][kBlk]->BCtypeS[i][k] =  BCtype_south;
                     
                  }
            }


      delete Grid_ptr_orig_blk;
          
   }
   
   if(IP.geometry_index == 2){

      Grid3D_Hexa_Block *Grid_ptr_orig_blk ;


      Grid_ptr_orig_blk = new Grid3D_Hexa_Block(
         NBI*IP.ICells,  NBJ*IP.JCells , NBK*IP.KCells, IP.Nghost);
      // generate the total mesh in the original geometry
    
      BCtype_west = BC_PERIODIC;
      BCtype_south =  BC_FIXED_PRESSURE;
      BCtype_bottom  = BC_NO_SLIP;
      BCtype_east = BC_PERIODIC;
      BCtype_north = BC_FIXED_PRESSURE ;
      BCtype_top  =  BC_MOVING_WALL;
      
      Grid_ptr_orig_blk->Create_Hexa_Block
         (  IP.Box_Length,
            IP.Box_Width,
            IP.Box_Height,
            - IP.Box_Length/2.0,
            0.0,
            - IP.Box_Height/2.0,          
            IP.I_stretching_factor,
            IP.J_stretching_factor,
            IP.K_stretching_factor,
            BCtype_top ,
            BCtype_bottom,
            BCtype_north,
            BCtype_south,
            BCtype_west,
            BCtype_east,
            NBI*IP.ICells,
            NBJ*IP.JCells,
            NBK*IP.KCells,
            IP.Nghost);
      // update cell coordinates
      Grid_ptr_orig_blk->Update_Cells();

      
      // copying the nodes coordinates from the parent block to 
      // subblocks.     
      for ( kBlk = 0; kBlk <NBK; ++kBlk ) 
         for ( jBlk = 0; jBlk <NBJ; ++jBlk ) 
            for ( iBlk = 0; iBlk <NBI; ++iBlk ) {
               
                if (iBlk==0 ){
                  BCtype_west = BC_PERIODIC;
               }else {
                  BCtype_west = BC_NONE;
               }
               if (jBlk==0 ){
                  BCtype_south = BC_FIXED_PRESSURE;
               }else {
                  BCtype_south = BC_NONE;
               }
               if (kBlk==0 ){
                  BCtype_bottom  = BC_NO_SLIP;
               }else {
                  BCtype_bottom = BC_NONE;
               }
               if (iBlk== NBI-1 ){
                  BCtype_east = BC_PERIODIC;
               }else {
                  BCtype_east = BC_NONE;
               }
               if (jBlk== NBJ -1 ){
                  BCtype_north = BC_FIXED_PRESSURE ;
               }else {
                  BCtype_north = BC_NONE;
               }
               if (kBlk== NBK-1 ){
                  BCtype_top  = BC_MOVING_WALL;
               }else {
                  BCtype_top = BC_NONE;
               }


            /* copy the cells and nodes information  
               of the hexa block to the subblocks. */

               for (int k =  Grid_ptr[iBlk][jBlk][kBlk]->KNl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                    k <=  Grid_ptr[iBlk][jBlk][kBlk]->KNu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++k) 
                  for (int j =  Grid_ptr[iBlk][jBlk][kBlk]->JNl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                       j <=  Grid_ptr[iBlk][jBlk][kBlk]->JNu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++j) 
                     for (int i =  Grid_ptr[iBlk][jBlk][kBlk]->INl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                          i <=  Grid_ptr[iBlk][jBlk][kBlk]->INu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++i) {
                        
                        Grid_ptr[iBlk][jBlk][kBlk]->Node[i][j][k].X.x = 
                           Grid_ptr_orig_blk->Node[iBlk*IP.ICells + i][jBlk*IP.JCells+ j][kBlk*IP.KCells +k].X.x;
                        
                        Grid_ptr[iBlk][jBlk][kBlk]->Node[i][j][k].X.y =  
                           Grid_ptr_orig_blk->Node[iBlk*IP.ICells + i][jBlk*IP.JCells+ j][kBlk*IP.KCells +k].X.y;
                        
                        Grid_ptr[iBlk][jBlk][kBlk]->Node[i][j][k].X.z =  
                           Grid_ptr_orig_blk->Node[iBlk*IP.ICells + i][jBlk*IP.JCells+ j][kBlk*IP.KCells +k].X.z;
                        
                     }
  
               Grid_ptr[iBlk][jBlk][kBlk]->Update_Cells();
               
               //store the boundary conditions
               for (int  j =  Grid_ptr[iBlk][jBlk][kBlk]->JCl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                    j <=  Grid_ptr[iBlk][jBlk][kBlk]->JCu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++j) 
                  for (int  i =  Grid_ptr[iBlk][jBlk][kBlk]->ICl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                       i <=  Grid_ptr[iBlk][jBlk][kBlk]->ICu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++i) {
                     Grid_ptr[iBlk][jBlk][kBlk]->BCtypeT[i][j] =  BCtype_top;
                     Grid_ptr[iBlk][jBlk][kBlk]->BCtypeB[i][j] =  BCtype_bottom;
            
                  }
               for ( int k =  Grid_ptr[iBlk][jBlk][kBlk]->KCl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                     k <=  Grid_ptr[iBlk][jBlk][kBlk]->KCu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++k) 
                  for ( int j =  Grid_ptr[iBlk][jBlk][kBlk]->JCl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                        j <=  Grid_ptr[iBlk][jBlk][kBlk]->JCu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++j){
                     Grid_ptr[iBlk][jBlk][kBlk]->BCtypeW[j][k] =  BCtype_west;
                     Grid_ptr[iBlk][jBlk][kBlk]->BCtypeE[j][k] =  BCtype_east;
                     
                  }
               for (int  k =  Grid_ptr[iBlk][jBlk][kBlk]->KCl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                    k <=  Grid_ptr[iBlk][jBlk][kBlk]->KCu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++k) 
                  for ( int i =  Grid_ptr[iBlk][jBlk][kBlk]->ICl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                        i <=  Grid_ptr[iBlk][jBlk][kBlk]->ICu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++i){
                     Grid_ptr[iBlk][jBlk][kBlk]->BCtypeN[i][k] =  BCtype_north;
                     Grid_ptr[iBlk][jBlk][kBlk]->BCtypeS[i][k] =  BCtype_south;
                     
                  }
            }


      delete Grid_ptr_orig_blk;
      
   }
   

   if(IP.geometry_index == 3){

     Grid3D_Hexa_Block *Grid_ptr_orig_blk ;


     Grid_ptr_orig_blk = new Grid3D_Hexa_Block(
        NBI*IP.ICells,  NBJ*IP.JCells , NBK*IP.KCells, IP.Nghost);
     // generate the total mesh in the original geometry
    
      BCtype_west = BC_PERIODIC;
      BCtype_south = BC_NO_SLIP;
      BCtype_bottom  =  BC_FIXED_PRESSURE;
      BCtype_east = BC_PERIODIC;
      BCtype_north = BC_MOVING_WALL;
      BCtype_top  =  BC_FIXED_PRESSURE;
      
      Grid_ptr_orig_blk->Create_Hexa_Block
         (  IP.Box_Length,
            IP.Box_Width,
            IP.Box_Height,
            - IP.Box_Length/2.0,
            - IP.Box_Width/2.0,
            0.0,
            IP.I_stretching_factor,
            IP.J_stretching_factor,
            IP.K_stretching_factor,
            BCtype_top ,
            BCtype_bottom,
            BCtype_north,
            BCtype_south,
            BCtype_west,
            BCtype_east,
            NBI*IP.ICells,
            NBJ*IP.JCells,
            NBK*IP.KCells,
            IP.Nghost);
      // update cell coordinates
      Grid_ptr_orig_blk->Update_Cells();


      
      for ( kBlk = 0; kBlk <NBK; ++kBlk ) 
         for ( jBlk = 0; jBlk <NBJ; ++jBlk ) 
            for ( iBlk = 0; iBlk <NBI; ++iBlk ) {
               
                              
               if (iBlk==0 ){
                  BCtype_west = BC_PERIODIC;
               }else {
                  BCtype_west = BC_NONE;
               }
               if (jBlk==0 ){
                  BCtype_south = BC_NO_SLIP;
               }else {
                  BCtype_south = BC_NONE;
               }
               if (kBlk==0 ){
                  BCtype_bottom  =  BC_FIXED_PRESSURE;
               }else {
                  BCtype_bottom = BC_NONE;
               }
               if (iBlk== NBI-1 ){
                  BCtype_east = BC_PERIODIC;
               }else {
                  BCtype_east = BC_NONE;
               }
               if (jBlk== NBJ -1 ){
                  BCtype_north = BC_MOVING_WALL;
               }else {
                  BCtype_north = BC_NONE;
               }
               if (kBlk== NBK-1 ){
                  BCtype_top  =  BC_FIXED_PRESSURE;
               }else {
                  BCtype_top = BC_NONE;
               }

            /* copy the cells and nodes information  
               of the hexa block to the subblocks. */

               for (int k =  Grid_ptr[iBlk][jBlk][kBlk]->KNl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                    k <=  Grid_ptr[iBlk][jBlk][kBlk]->KNu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++k) 
                  for (int j =  Grid_ptr[iBlk][jBlk][kBlk]->JNl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                       j <=  Grid_ptr[iBlk][jBlk][kBlk]->JNu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++j) 
                     for (int i =  Grid_ptr[iBlk][jBlk][kBlk]->INl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                          i <=  Grid_ptr[iBlk][jBlk][kBlk]->INu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++i) {
                        
                        Grid_ptr[iBlk][jBlk][kBlk]->Node[i][j][k].X.x = 
                           Grid_ptr_orig_blk->Node[iBlk*IP.ICells + i][jBlk*IP.JCells+ j][kBlk*IP.KCells +k].X.x;
                        
                        Grid_ptr[iBlk][jBlk][kBlk]->Node[i][j][k].X.y =  
                           Grid_ptr_orig_blk->Node[iBlk*IP.ICells + i][jBlk*IP.JCells+ j][kBlk*IP.KCells +k].X.y;
                        
                        Grid_ptr[iBlk][jBlk][kBlk]->Node[i][j][k].X.z =  
                           Grid_ptr_orig_blk->Node[iBlk*IP.ICells + i][jBlk*IP.JCells+ j][kBlk*IP.KCells +k].X.z;
                        
                     }
  
               Grid_ptr[iBlk][jBlk][kBlk]->Update_Cells();
               
               //store the boundary conditions
               for (int  j =  Grid_ptr[iBlk][jBlk][kBlk]->JCl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                    j <=  Grid_ptr[iBlk][jBlk][kBlk]->JCu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++j) 
                  for (int  i =  Grid_ptr[iBlk][jBlk][kBlk]->ICl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                       i <=  Grid_ptr[iBlk][jBlk][kBlk]->ICu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++i) {
                     Grid_ptr[iBlk][jBlk][kBlk]->BCtypeT[i][j] =  BCtype_top;
                     Grid_ptr[iBlk][jBlk][kBlk]->BCtypeB[i][j] =  BCtype_bottom;
            
                  }
               for ( int k =  Grid_ptr[iBlk][jBlk][kBlk]->KCl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                     k <=  Grid_ptr[iBlk][jBlk][kBlk]->KCu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++k) 
                  for ( int j =  Grid_ptr[iBlk][jBlk][kBlk]->JCl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                        j <=  Grid_ptr[iBlk][jBlk][kBlk]->JCu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++j){
                     Grid_ptr[iBlk][jBlk][kBlk]->BCtypeW[j][k] =  BCtype_west;
                     Grid_ptr[iBlk][jBlk][kBlk]->BCtypeE[j][k] =  BCtype_east;
                     
                  }
               for (int  k =  Grid_ptr[iBlk][jBlk][kBlk]->KCl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                    k <=  Grid_ptr[iBlk][jBlk][kBlk]->KCu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++k) 
                  for ( int i =  Grid_ptr[iBlk][jBlk][kBlk]->ICl- Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; 
                        i <=  Grid_ptr[iBlk][jBlk][kBlk]->ICu+ Grid_ptr[iBlk][jBlk][kBlk]->Nghost ; ++i){
                     Grid_ptr[iBlk][jBlk][kBlk]->BCtypeN[i][k] =  BCtype_north;
                     Grid_ptr[iBlk][jBlk][kBlk]->BCtypeS[i][k] =  BCtype_south;
                     
                  }
            }


      delete Grid_ptr_orig_blk;

   }
   

  //  if(IP.geometry_index == 1){

//       for ( kBlk = 0; kBlk <NBK; ++kBlk ) {
//          for ( jBlk = 0; jBlk <NBJ; ++jBlk ) {
//             for ( iBlk = 0; iBlk <NBI; ++iBlk ) {
//                x_orig = iBlk*length;
//                y_orig = jBlk*width - IP.Box_Width/2.0;
//                z_orig = kBlk*height - IP.Box_Height/2.0;
               
//                if (iBlk==0 ){
//                   BCtype_west = BC_FIXED_PRESSURE;
//                }else {
//                   BCtype_west = BC_NONE;
//                }
//                if (jBlk==0 ){
//                   BCtype_south = BC_PERIODIC;
//                }else {
//                   BCtype_south = BC_NONE;
//                }
//                if (kBlk==0 ){
//                   BCtype_bottom  = BC_NO_SLIP;
//                }else {
//                   BCtype_bottom = BC_NONE;
//                }
//                if (iBlk== NBI-1 ){
//                   BCtype_east = BC_FIXED_PRESSURE;
//                }else {
//                   BCtype_east = BC_NONE;
//                }
//                if (jBlk== NBJ -1 ){
//                   BCtype_north = BC_PERIODIC;
//                }else {
//                   BCtype_north = BC_NONE;
//                }
//                if (kBlk== NBK-1 ){
//                   BCtype_top  = BC_MOVING_WALL;
//                }else {
//                   BCtype_top = BC_NONE;
//                }

//             /* Allocate (re-allocate) memory for the cells and nodes 
//                of the hexa block. */
//             Grid_ptr[iBlk][jBlk][kBlk]->Create_Hexa_Block
//                (length,
//                 width,
//                 height,
//                 x_orig,
//                 y_orig,
//                 z_orig,
//                 IP.I_stretching_factor,
//                 IP.J_stretching_factor,
//                 IP.K_stretching_factor,
//                 BCtype_top,
//                 BCtype_bottom,
//                 BCtype_north,
//                 BCtype_south,
//                 BCtype_west,
//                 BCtype_east,
//                 IP.ICells,
//                 IP.JCells,
//                 IP.KCells,
//                 IP.Nghost);
            
//             }
//          }
//       }
//    }
   
//    if(IP.geometry_index == 2){

//       for ( kBlk = 0; kBlk <NBK; ++kBlk ) {
//          for ( jBlk = 0; jBlk <NBJ; ++jBlk ) {
//             for ( iBlk = 0; iBlk <NBI; ++iBlk ) {
            
//                x_orig = iBlk*length - IP.Box_Length/2.0;
//                y_orig = jBlk*width ;
//                z_orig = kBlk*height - IP.Box_Height/2.0;
               
//                if (iBlk==0 ){
//                   BCtype_west = BC_PERIODIC;
//                }else {
//                   BCtype_west = BC_NONE;
//                }
//                if (jBlk==0 ){
//                   BCtype_south = BC_FIXED_PRESSURE;
//                }else {
//                   BCtype_south = BC_NONE;
//                }
//                if (kBlk==0 ){
//                   BCtype_bottom  = BC_NO_SLIP;
//                }else {
//                   BCtype_bottom = BC_NONE;
//                }
//                if (iBlk== NBI-1 ){
//                   BCtype_east = BC_PERIODIC;
//                }else {
//                   BCtype_east = BC_NONE;
//                }
//                if (jBlk== NBJ -1 ){
//                   BCtype_north = BC_FIXED_PRESSURE ;
//                }else {
//                   BCtype_north = BC_NONE;
//                }
//                if (kBlk== NBK-1 ){
//                   BCtype_top  = BC_MOVING_WALL;
//                }else {
//                   BCtype_top = BC_NONE;
//                }

//             /* Allocate (re-allocate) memory for the cells and nodes 
//                of the hexa block. */
//             Grid_ptr[iBlk][jBlk][kBlk]->Create_Hexa_Block
//                (length,
//                 width,
//                 height,
//                 x_orig,
//                 y_orig,
//                 z_orig,
//                 IP.I_stretching_factor,
//                 IP.J_stretching_factor,
//                 IP.K_stretching_factor,
//                 BCtype_top,
//                 BCtype_bottom,
//                 BCtype_north,
//                 BCtype_south,
//                 BCtype_west,
//                 BCtype_east,
//                 IP.ICells,
//                 IP.JCells,
//                 IP.KCells,
//                 IP.Nghost);
            
//             }
//          }
//       }
//    }

//    if(IP.geometry_index == 3){
      
//       for ( kBlk = 0; kBlk <NBK; ++kBlk ) {
//          for ( jBlk = 0; jBlk <NBJ; ++jBlk ) {
//             for ( iBlk = 0; iBlk <NBI; ++iBlk ) {
               
//                x_orig = iBlk*length - IP.Box_Length/2.0 ;
//                y_orig = jBlk*width -  IP.Box_Width/2.0 ;
//                z_orig = kBlk*height ;
                             
                              
//                if (iBlk==0 ){
//                   BCtype_west = BC_PERIODIC;
//                }else {
//                   BCtype_west = BC_NONE;
//                }
//                if (jBlk==0 ){
//                   BCtype_south = BC_NO_SLIP;
//                }else {
//                   BCtype_south = BC_NONE;
//                }
//                if (kBlk==0 ){
//                   BCtype_bottom  = BC_CHANNEL_INFLOW;
//                }else {
//                   BCtype_bottom = BC_NONE;
//                }
//                if (iBlk== NBI-1 ){
//                   BCtype_east = BC_PERIODIC;
//                }else {
//                   BCtype_east = BC_NONE;
//                }
//                if (jBlk== NBJ -1 ){
//                   BCtype_north = BC_MOVING_WALL;
//                }else {
//                   BCtype_north = BC_NONE;
//                }
//                if (kBlk== NBK-1 ){
//                   BCtype_top  = BC_CHANNEL_OUTFLOW;
//                }else {
//                   BCtype_top = BC_NONE;
//                }

//             /* Allocate (re-allocate) memory for the cells and nodes 
//                of the hexa block. */
//             Grid_ptr[iBlk][jBlk][kBlk]->Create_Hexa_Block
//                (length,
//                 width,
//                 height,
//                 x_orig,
//                 y_orig,
//                 z_orig,
//                 IP.I_stretching_factor,
//                 IP.J_stretching_factor,
//                 IP.K_stretching_factor,
//                 BCtype_top,
//                 BCtype_bottom,
//                 BCtype_north,
//                 BCtype_south,
//                 BCtype_west,
//                 BCtype_east,
//                 IP.ICells,
//                 IP.JCells,
//                 IP.KCells,
//                 IP.Nghost);
            
//             }
//          }
//       }
//    }
   
}// Couette

void Grid3D_Hexa_Multi_Block::Deallocate_Unused_Grid_Blocks(void) {
   
   
   for (int i = 0 ; i <NBI ; ++i ) 
      for (int j = 0 ; j <NBJ ; ++j ) 
         for (int k = 0 ; k <NBK ; ++k ){
            
            if(Grid_ptr[i][j][k]->Used ==0){
               delete Grid_ptr[i][j][k]; //delete the unused grid block
            }
            
         }
}//endofdeallocate

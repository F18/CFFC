/*!\file HO_Grid2DQuadMultiBlock.cc
   \brief Multi-block subroutines for high-order 2D quadrilateral multi-block grid class. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "HO_Grid2DQuadMultiBlock.h"  /* Include 2D quadrilateral multi-block grid type header file. */
#include "../Grid/NASARotor37.h"       /* Include NASA rotor 37 header files. */
#include "../Grid/NASARotor67.h"       /* Include NASA rotor 67 header files. */


// ===== Member functions =====

/*!
 * Default constructor.
 */
Grid2D_Quad_MultiBlock_HO::Grid2D_Quad_MultiBlock_HO(void)
  : Number_of_Blocks_Idir(0), Number_of_Blocks_Jdir(0),
    Grid_ptr(NULL)
{
  //
}

/*!
 * Copy constructor. It is declared private
 */
Grid2D_Quad_MultiBlock_HO::Grid2D_Quad_MultiBlock_HO(const Grid2D_Quad_MultiBlock_HO &G)
  : Number_of_Blocks_Idir(0), Number_of_Blocks_Jdir(0),
    Grid_ptr(NULL)
{
  int i,j;

  allocate(G.Blocks_Idir(), G.Blocks_Jdir());

  // Copy the individual blocks
  if (G.Grid_ptr != NULL) {
    for (j = 0; j<=Last_jBlock(); ++j){
      for (i = 0; i<=Last_iBlock(); ++i){
	Grid_ptr[i][j] = G(i,j);
      }	// endfor
    } // endfor
  } // endif
}

/*!
 * Allocate memory for the 2D quadrilateral multi-block grid.
 *
 * \param _Number_of_Blocks_Idir_ number of blocks in i-direction
 * \param _Number_of_Blocks_Jdir_ number of blocks in j-direction
 */
void Grid2D_Quad_MultiBlock_HO::allocate(const int & _Number_of_Blocks_Idir_,
					 const int & _Number_of_Blocks_Jdir_) {

  // Check conditions
  assert( _Number_of_Blocks_Idir_ >= 1 && _Number_of_Blocks_Jdir_ >= 1 );

  // Check if the new required memory has dimensions different than the currently allocated ones
  if ( _Number_of_Blocks_Idir_ != Number_of_Blocks_Idir || 
       _Number_of_Blocks_Jdir_ != Number_of_Blocks_Jdir ){

    // Free the memory if there is memory allocated
    deallocate();
    
    // Set mesh parameters
    Number_of_Blocks_Idir = _Number_of_Blocks_Idir_;
    Number_of_Blocks_Jdir = _Number_of_Blocks_Jdir_;

    /* Allocate memory. */
    Grid_ptr = new Grid2D_Quad_Block_HO*[Number_of_Blocks_Idir];
    for (int i = 0 ; i <= Last_iBlock() ; ++i ) {
      Grid_ptr[i] = new Grid2D_Quad_Block_HO[Number_of_Blocks_Jdir];
    } /* endfor */

  } // endif
}

/*
 * Deallocate memory.
 */
void Grid2D_Quad_MultiBlock_HO::deallocate(void) {

  if (Grid_ptr != NULL){
  
    // deallocate the j-direction blocks
    for (int i = 0 ; i <= Last_iBlock() ; ++i ) {
      delete []Grid_ptr[i];
      Grid_ptr[i] = NULL;
    }

    // deallocate the i-direction blocks
    delete []Grid_ptr;
    Grid_ptr = NULL;

    // Reset indexes
    Number_of_Blocks_Idir = 0;
    Number_of_Blocks_Jdir = 0;
  }
}

/*!
 * Assignment operator =
 */
Grid2D_Quad_MultiBlock_HO& Grid2D_Quad_MultiBlock_HO::operator=(const Grid2D_Quad_MultiBlock_HO &G) {

  // !!! If the LHS grid block already has objects assigned, these are going to be deleted.
  // Handle self-assignment:
  if (this == &G) return *this;

  int i,j;

  // re-allocate memory if there isn't enough
  allocate(G.Blocks_Idir(),G.Blocks_Jdir());

  // Copy the values from each single block.
  if (G.Grid_ptr != NULL) {
    for (j = 0; j<=Last_jBlock(); ++j){
      for (i = 0; i<=Last_iBlock(); ++i){
	Grid_ptr[i][j] = G(i,j);
      }	// endfor
    } // endfor
  } // endif
}


/*!
 * Broadcasts the multi-block grid to all processors 
 * involved in the calculation from the primary processor 
 * using the MPI broadcast routine.
 */
void Grid2D_Quad_MultiBlock_HO::Broadcast_Multi_Block_Grid(void) {

#ifdef _MPI_VERSION
  int i, j;
  int Num_iBlocks(0), Num_jBlocks(0);
  
  if (CFFC_Primary_MPI_Processor()) {
    // initialize the number of blocks with values from the primary CPU
    Num_iBlocks = Number_of_Blocks_Idir;
    Num_jBlocks = Number_of_Blocks_Jdir;
  }

  // Broadcast the number of blocks in both directions
  CFFC_Broadcast_MPI(&Num_iBlocks,1);
  CFFC_Broadcast_MPI(&Num_jBlocks,1);

  // Re-allocate memory on CPUs different than the primary one
  if (!CFFC_Primary_MPI_Processor()) {
    allocate(Num_iBlocks,Num_jBlocks);
  } /* endif */
  
  // Broadcast the individual grid blocks
  for ( j = 0 ; j <= Last_jBlock() ; ++j ) {
    for ( i = 0; i <= Last_iBlock() ; ++i ) {
      Grid_ptr[i][j].Broadcast_Quad_Block();
    }  /* endfor */
  }  /* endfor */
#endif
}



/*!
 * Writes definition file information for a 2D array of 
 * 2D quadrilateral multi-block grids to the specified  
 * output stream for retrieval and re-use purposes.     
 */
void Grid2D_Quad_MultiBlock_HO::Write_Multi_Block_Grid_Definition(ostream &Out_File) {

    int i, j;
 
    Out_File << Number_of_Blocks_Idir << " " 
             << Number_of_Blocks_Jdir << "\n";

    for ( j = 0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
       for ( i = 0; i <= Number_of_Blocks_Idir-1 ; ++i ) {
 	  Grid_ptr[i][j].Write_Quad_Block_Definition(Out_File);
       }  /* endfor */
    }  /* endfor */

}


/*!
 * Reads definition file information for a 2D array of  
 * 2D quadrilateral multi-block grids from the          
 * specified input stream.                              
 */
void Grid2D_Quad_MultiBlock_HO::Read_Multi_Block_Grid_Definition(istream &In_File) {
  
  int i, j;
  
  In_File.setf(ios::skipws);
  In_File >> i >> j;
  In_File.unsetf(ios::skipws);
  
  // allocate/re-allocate memory for the multi-block grid
  allocate(i, j);
  
  for ( j = 0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
    for ( i = 0; i <= Number_of_Blocks_Idir-1 ; ++i ) {
      Grid_ptr[i][j].Read_Quad_Block_Definition(In_File);
    }  /* endfor */
  }  /* endfor */
  
}


/*!
 * Writes a 2D array of 2D quadrilateral multi-block    
 * grids to the specified output stream for retrieval   
 * and re-use purposes.                                 
 *                                                      
 */
void Grid2D_Quad_MultiBlock_HO::Write_Multi_Block_Grid(ostream &Out_File) {

    int i, j;
 
    Out_File << Number_of_Blocks_Idir << " " 
             << Number_of_Blocks_Jdir << "\n";

    for ( j = 0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
       for ( i = 0; i <= Number_of_Blocks_Idir-1 ; ++i ) {
 	  Grid_ptr[i][j].Write_Quad_Block(Out_File);
       }  /* endfor */
    }  /* endfor */

}

/*
 * Writes a 2D array of 2D quadrilateral multi-block    
 * grids from the specified input stream.               
 */
void Grid2D_Quad_MultiBlock_HO::Read_Multi_Block_Grid(istream &In_File) {
    
  int i, j;
  
  In_File.setf(ios::skipws);
  In_File >> i >> j;
  In_File.unsetf(ios::skipws);

  // allocate/re-allocate memory for the multi-block grid
  allocate(i, j);

  for ( j = 0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
    for ( i = 0; i <= Number_of_Blocks_Idir-1 ; ++i ) {
      Grid_ptr[i][j].Read_Quad_Block(In_File);
    }  /* endfor */
  }  /* endfor */
}

/*!
 * Translate the nodes of a 2D array of 2D              
 * quadrilateral multi-block grids.                     
 */
void Grid2D_Quad_MultiBlock_HO::Translate_Multi_Block_Grid(const Vector2D &V) {

  int i, j;
  
  for ( j = 0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
    for ( i = 0; i <= Number_of_Blocks_Idir-1 ; ++i ) {
      if (Grid_ptr[i][j].Node != NULL) {
	Grid_ptr[i][j].Translate_Quad_Block(V);
      } /* endif */
    }  /* endfor */
  }  /* endfor */

}

/*!
 * Translate the nodes of a 2D array of 2D              
 * quadrilateral multi-block grids.   
 * This subroutine DOESN'T update the geometric properties
 * of the block cells.                  
 */
void Grid2D_Quad_MultiBlock_HO::Translate_Multi_Block_Grid_Without_Update(const Vector2D &V) {

  int i, j;
  
  for ( j = 0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
    for ( i = 0; i <= Number_of_Blocks_Idir-1 ; ++i ) {
      if (Grid_ptr[i][j].Node != NULL) {
	Grid_ptr[i][j].Translate_Quad_Block_Without_Update(V);
      } /* endif */
    }  /* endfor */
  }  /* endfor */

}

/*
 * Scales 2D array of 2D quadrilateral multi-block      
 * grids.
 */
void Grid2D_Quad_MultiBlock_HO::Scale_Multi_Block_Grid(const double &Scaling_Factor) {

  int i, j;
  
  for ( j = 0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
    for ( i = 0; i <= Number_of_Blocks_Idir-1 ; ++i ) {
      if (Grid_ptr[i][j].Node != NULL) {
	Grid_ptr[i][j].Scale_Quad_Block(Scaling_Factor);
      } /* endif */
    }  /* endfor */
  }  /* endfor */

}

/*
 * Scales 2D array of 2D quadrilateral multi-block      
 * grids.                               
 * This subroutine DOESN'T update the geometric properties
 * of the block cells.                                  
 */
void Grid2D_Quad_MultiBlock_HO::Scale_Multi_Block_Grid_Without_Update(const double &Scaling_Factor) {

  int i, j;
  
  for ( j = 0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
    for ( i = 0; i <= Number_of_Blocks_Idir-1 ; ++i ) {
      if (Grid_ptr[i][j].Node != NULL) {
	Grid_ptr[i][j].Scale_Quad_Block_Without_Update(Scaling_Factor);
      } /* endif */
    }  /* endfor */
  }  /* endfor */

}

/*
 * Rotates 2D array of 2D quadrilateral multi-block     
 * grids.                                               
 */
void Grid2D_Quad_MultiBlock_HO::Rotate_Multi_Block_Grid(const double &Angle) {

  int i, j;

  for ( j = 0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
    for ( i = 0; i <= Number_of_Blocks_Idir-1 ; ++i ) {
      if (Grid_ptr[i][j].Node != NULL) {
	Grid_ptr[i][j].Rotate_Quad_Block(Angle);
      } /* endif */
    }  /* endfor */
  }  /* endfor */
  
}

/*
 * Rotates 2D array of 2D quadrilateral multi-block     
 * grids.                                               
 * This subroutine DOESN'T update the geometric properties
 * of the block cells.                                                                 
 */
void Grid2D_Quad_MultiBlock_HO::Rotate_Multi_Block_Grid_Without_Update(const double &Angle) {

  int i, j;

  for ( j = 0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
    for ( i = 0; i <= Number_of_Blocks_Idir-1 ; ++i ) {
      if (Grid_ptr[i][j].Node != NULL) {
	Grid_ptr[i][j].Rotate_Quad_Block_Without_Update(Angle);
      } /* endif */
    }  /* endfor */
  }  /* endfor */
  
}

/*!
 * Reflects 2D array of 2D quadrilateral multi-block    
 * grids about y=0 axis.                                
 */
void Grid2D_Quad_MultiBlock_HO::Reflect_Multi_Block_Grid(void) {

  int i, j;
  
  for ( j = 0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
    for ( i = 0; i <= Number_of_Blocks_Idir-1 ; ++i ) {
      if (Grid_ptr[i][j].Node != NULL) {
	Grid_ptr[i][j].Reflect_Quad_Block();
      } /* endif */
    }  /* endfor */
  }  /* endfor */
}

/*!
 * Reflects 2D array of 2D quadrilateral multi-block    
 * grids about y=0 axis.                                
 * This subroutine DOESN'T update the geometric properties
 * of the block cells.
 */
void Grid2D_Quad_MultiBlock_HO::Reflect_Multi_Block_Grid_Without_Update(void) {

  int i, j;
  
  for ( j = 0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
    for ( i = 0; i <= Number_of_Blocks_Idir-1 ; ++i ) {
      if (Grid_ptr[i][j].Node != NULL) {
	Grid_ptr[i][j].Reflect_Quad_Block_Without_Update();
      } /* endif */
    }  /* endfor */
  }  /* endfor */
}

/*!
 * Check the validity of 2D array of 2D quadrilateral   
 * multi-block grids.  Returns a non-zero result if     
 * mesh is not valid.                                   
 */
int Grid2D_Quad_MultiBlock_HO::Check_Multi_Block_Grid(void) {

  int i, j;
  
  for ( j = 0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
    for ( i = 0; i <= Number_of_Blocks_Idir-1 ; ++i ) {
      if (Grid_ptr[i][j].Node != NULL) {
	if (Grid_ptr[i][j].Check_Quad_Block()) return(1);
      } /* endif */
    }  /* endfor */
  }  /* endfor */
  
  return(0);

}

/*!
 * Check the validity of 2D array of 2D quadrilateral   
 * multi-block grids.  Returns a non-zero result if     
 * mesh is not valid.
 * This subroutine checks the validity of the ghost cells too.                                   
 */
int Grid2D_Quad_MultiBlock_HO::Check_Multi_Block_Grid_Completely(void) {

  int i, j;
  
  for ( j = 0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
    for ( i = 0; i <= Number_of_Blocks_Idir-1 ; ++i ) {
      if (Grid_ptr[i][j].Node != NULL) {
	if (Grid_ptr[i][j].Check_Quad_Block_Completely()) return(1);
      } /* endif */
    }  /* endfor */
  }  /* endfor */
  
  return(0);

}

/*!
 * Disturb randomly the interior nodes of 2D array of 2D quadrilateral   
 * multi-block grids.
 * This subroutine DOESN'T update the geometric properties
 * of the block cells.
 */
void Grid2D_Quad_MultiBlock_HO::Disturb_Interior_Nodes_Without_Update(const int &Number_of_Iterations){

  int i, j;
  
  for ( j = 0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
    for ( i = 0; i <= Number_of_Blocks_Idir-1 ; ++i ) {
      if (Grid_ptr[i][j].Node != NULL) {
	Grid_ptr[i][j].Disturb_Interior_Nodes_Without_Update(Number_of_Iterations);
      } /* endif */
    }  /* endfor */
  }  /* endfor */
}

/*!
 * Disturb randomly the interior nodes of 2D array of 2D quadrilateral   
 * multi-block grids.
 */
void Grid2D_Quad_MultiBlock_HO::Disturb_Interior_Nodes(const int &Number_of_Iterations){

  int i, j;
  
  for ( j = 0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
    for ( i = 0; i <= Number_of_Blocks_Idir-1 ; ++i ) {
      if (Grid_ptr[i][j].Node != NULL) {
	Grid_ptr[i][j].Disturb_Interior_Nodes(Number_of_Iterations);
      } /* endif */
    }  /* endfor */
  }  /* endfor */
}

/*!
 * Setup the required flux calculation method based on the flags set in 
 * HO_Grid2D_Execution_Mode class for each 2D quadrilateral block of the multi-block grids.
 */
void Grid2D_Quad_MultiBlock_HO::SetFluxCalculationMethod(void){
  
  int iBlock, jBlock;

  if (HO_Grid2D_Execution_Mode::CUSTOMIZE_FLUX_CALCULATION_METHOD_AT_BOUNDARIES){
    
    for (iBlock = 0; iBlock <= Last_iBlock() ; ++iBlock){
      for (jBlock = 0; jBlock <= Last_jBlock() ; ++jBlock){

	// North boundary
	if ( HO_Grid2D_Execution_Mode::NORTH_RECONSTRUCTION_BASED_FLUX && 
	     Grid_ptr[iBlock][jBlock].BndNorthSpline.bc[0] != BC_NONE ){
	  Grid_ptr[iBlock][jBlock].BndNorthSpline.setFluxCalcMethod(ReconstructionBasedFlux);
	} else {
	  Grid_ptr[iBlock][jBlock].BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem);
	}

	// South boundary
	if( HO_Grid2D_Execution_Mode::SOUTH_RECONSTRUCTION_BASED_FLUX && 
	    Grid_ptr[iBlock][jBlock].BndSouthSpline.bc[0] != BC_NONE ){
	  Grid_ptr[iBlock][jBlock].BndSouthSpline.setFluxCalcMethod(ReconstructionBasedFlux);
	} else {
	  Grid_ptr[iBlock][jBlock].BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem);
	}

	// East boundary
	if( HO_Grid2D_Execution_Mode::EAST_RECONSTRUCTION_BASED_FLUX && 
	    Grid_ptr[iBlock][jBlock].BndEastSpline.bc[0] != BC_NONE){
	  Grid_ptr[iBlock][jBlock].BndEastSpline.setFluxCalcMethod(ReconstructionBasedFlux);
	} else {
	  Grid_ptr[iBlock][jBlock].BndEastSpline.setFluxCalcMethod(SolveRiemannProblem);
	}

	// West boundary
	if( HO_Grid2D_Execution_Mode::WEST_RECONSTRUCTION_BASED_FLUX && 
	    Grid_ptr[iBlock][jBlock].BndWestSpline.bc[0] != BC_NONE){
	  Grid_ptr[iBlock][jBlock].BndWestSpline.setFluxCalcMethod(ReconstructionBasedFlux);
	} else {
	  Grid_ptr[iBlock][jBlock].BndWestSpline.setFluxCalcMethod(SolveRiemannProblem);
	}
      }	// endfor
    } // endfor
    
  } // endif
}

/*!
 * Update the exterior nodes in all 
 * mesh blocks.
 */
void Grid2D_Quad_MultiBlock_HO::Update_All_Exterior_Nodes(void){
  int i, j;
  
  for ( j = 0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
    for ( i = 0; i <= Number_of_Blocks_Idir-1 ; ++i ) {
      if (Grid_ptr[i][j].Node != NULL) {
	Grid_ptr[i][j].Update_Exterior_Nodes();
      } /* endif */
    }  /* endfor */
  }  /* endfor */  
}

/*!
 * Update cell geometric properties in all
 * mesh blocks.
 */
void Grid2D_Quad_MultiBlock_HO::Update_All_Cells(void){
  int i, j;
  
  for ( j = 0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
    for ( i = 0; i <= Number_of_Blocks_Idir-1 ; ++i ) {
      if (Grid_ptr[i][j].Node != NULL) {
	Grid_ptr[i][j].Update_Cells();
      } /* endif */
    }  /* endfor */
  }  /* endfor */  
}

/*!
 * Writes the nodes of a 2D array of 2D quadrilateral   
 * multi-block grids to the specified output stream in  
 * a format suitable for plotting the grid with         
 * TECPLOT.                                                                                                    *
 */
void Grid2D_Quad_MultiBlock_HO::Output_Tecplot(ostream &Out_File) {

  int i, j, block_number, i_output_title;
  
  block_number = 0;
  i_output_title = 1;
  
  for ( j = 0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
    for ( i = 0; i <= Number_of_Blocks_Idir-1 ; ++i ) {
      if (Grid_ptr[i][j].Node != NULL) {
	Grid_ptr[i][j].Output_Tecplot(block_number,
				      i_output_title,
				      Out_File);
	block_number = block_number + 1;
	if (i_output_title) i_output_title = 0;
      } /* endif */
    }  /* endfor */
  }  /* endfor */

}

/*!
 * Writes the nodes of a 2D array of 2D quadrilateral   
 * multi-block grids to the specified output stream in  
 * a format suitable for plotting the grid with         
 * TECPLOT.  Include boundary nodes.                    
 */
void Grid2D_Quad_MultiBlock_HO::Output_Nodes_Tecplot(ostream &Out_File) {
  
  int i, j, block_number, i_output_title;
  
  block_number = 0;
  i_output_title = 1;
  
  for ( j = 0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
    for ( i = 0; i <= Number_of_Blocks_Idir-1 ; ++i ) {
      if (Grid_ptr[i][j].Node != NULL) {
	Grid_ptr[i][j].Output_Nodes_Tecplot(block_number,
					    i_output_title,
					    Out_File);
	block_number = block_number + 1;
	if (i_output_title) i_output_title = 0;
      } /* endif */
    }  /* endfor */
  }  /* endfor */
  
}

/*
 * Writes the cells of a 2D array of 2D quadrilateral   
 * multi-block grids to the specified output stream in  
 * a format suitable for plotting the grid with         
 * TECPLOT.                                             
 */
void Grid2D_Quad_MultiBlock_HO::Output_Cells_Tecplot(ostream &Out_File) {

  int i, j, block_number, i_output_title;
  
  block_number = 0;
  i_output_title = 1;
  
  for ( j = 0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
    for ( i = 0; i <= Number_of_Blocks_Idir-1 ; ++i ) {
      if (Grid_ptr[i][j].Node != NULL) {
	Grid_ptr[i][j].Output_Cells_Tecplot(block_number,
					    i_output_title,
					    Out_File);
	block_number = block_number + 1;
	if (i_output_title) i_output_title = 0;
      } /* endif */
    }  /* endfor */
  }  /* endfor */

}

/*!
 * Writes the nodes of a 2D array of 2D quadrilateral   
 * multi-block grids to the specified output stream in  
 * a format suitable for plotting the grid with         
 * GNUPLOT.                                             
 *                                                      
 */
void Grid2D_Quad_MultiBlock_HO::Output_Gnuplot(ostream &Out_File) {

  int i, j, block_number, i_output_title;
  
  block_number = 0;
  i_output_title = 1;
  
  for ( j = 0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
    for ( i = 0; i <= Number_of_Blocks_Idir-1 ; ++i ) {
      if (Grid_ptr[i][j].Node != NULL) {
	Grid_ptr[i][j].Output_Gnuplot(block_number,
				      i_output_title,
				      Out_File);
	block_number = block_number + 1;
	if (i_output_title) i_output_title = 0;
      } /* endif */
    }  /* endfor */
  }  /* endfor */
  
}


/*!
 * Generates a uniform 2D Cartesian mesh for a          
 * rectangular box shaped domain.                       
 *                                                      
 * Usage: Grid_ptr = Grid_Rectangular_Box(Grid_ptr,     
 *                                        nblk_i,       
 *                                        nblk_j,       
 *                                        TEN,          
 *                                        FIVE,         
 *   	                                  100,          
 *  	                                  50,           
 *                                        2);           
 *                                                      
 * This subroutine DOESN'T update the ghost cells or
 * the geometric properties of the grid cells.
 */
void Grid2D_Quad_MultiBlock_HO::Grid_Rectangular_Box_Without_Update(int &_Number_of_Blocks_Idir_,
								    int &_Number_of_Blocks_Jdir_,
								    const double &Width,
								    const double &Height,
								    const int Number_of_Cells_Idir,
								    const int Number_of_Cells_Jdir,
								    const int Number_of_Ghost_Cells,
								    const int Highest_Order_of_Reconstruction) {
  
  int iBlk, jBlk, n_cells_i, n_cells_j, 
    Stretch_I, Stretch_J,
    Orthogonal_North, Orthogonal_South,
    Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D xc_NW, xc_NE, xc_SE, xc_SW;
  Spline2D_HO Bnd_Spline_North, Bnd_Spline_South,
              Bnd_Spline_East, Bnd_Spline_West;
  
  /* Allocate memory for grid block. */
  
  if (_Number_of_Blocks_Idir_ < 0) _Number_of_Blocks_Idir_ = 1;
  if (_Number_of_Blocks_Jdir_ < 0) _Number_of_Blocks_Jdir_ = 1;
  allocate(_Number_of_Blocks_Idir_, _Number_of_Blocks_Jdir_);
  
  /* Create the mesh for each block representing
     the complete grid. */
  
  for ( jBlk = 0; jBlk <= Number_of_Blocks_Jdir-1; ++jBlk ) {
    for ( iBlk = 0; iBlk <= Number_of_Blocks_Idir-1; ++iBlk ) {
      
      /* Assign values to the locations of the corners
	 of the rectangular box shaped domain. */
      
      xc_NW = Vector2D(-HALF*Width+double(iBlk)*Width/double(Number_of_Blocks_Idir), 
		       -HALF*Height+double(jBlk+1)*Height/double(Number_of_Blocks_Jdir));
      xc_NE = Vector2D(-HALF*Width+double(iBlk+1)*Width/double(Number_of_Blocks_Idir), 
		       -HALF*Height+double(jBlk+1)*Height/double(Number_of_Blocks_Jdir));
      xc_SE = Vector2D(-HALF*Width+double(iBlk+1)*Width/double(Number_of_Blocks_Idir), 
		       -HALF*Height+double(jBlk)*Height/double(Number_of_Blocks_Jdir));
      xc_SW = Vector2D(-HALF*Width+double(iBlk)*Width/double(Number_of_Blocks_Idir), 
		       -HALF*Height+double(jBlk)*Height/double(Number_of_Blocks_Jdir));
      
      /* Create the splines defining the north, south,
	 east, and west boundaries of the rectangular box. */
      
      Bnd_Spline_North.Create_Spline_Line(xc_NW, xc_NE, 2);
      Bnd_Spline_South.Create_Spline_Line(xc_SW, xc_SE, 2);
      Bnd_Spline_East.Create_Spline_Line(xc_SE, xc_NE, 2);
      Bnd_Spline_West.Create_Spline_Line(xc_SW, xc_NW, 2);
      
      /* Set the boundary condition types for each of the
	 boundary splines. */

      if (jBlk == Number_of_Blocks_Jdir-1) {
	Bnd_Spline_North.setBCtype(BC_REFLECTION);
      } else {
	Bnd_Spline_North.setBCtype(BC_NONE);
      } /* endif */
      if (jBlk == 0) {
	Bnd_Spline_South.setBCtype(BC_REFLECTION);
      } else {
	Bnd_Spline_South.setBCtype(BC_NONE);
      } /* endif */
      if (iBlk == Number_of_Blocks_Idir-1) {
	Bnd_Spline_East.setBCtype(BC_REFLECTION);
      } else {
	Bnd_Spline_East.setBCtype(BC_NONE);
      } /* endif */
      if (iBlk == 0) {
	Bnd_Spline_West.setBCtype(BC_REFLECTION);
      } else {
	Bnd_Spline_West.setBCtype(BC_NONE);
      } /* endif */
      
      /* Assign values to the stretching function parameters
	 and boundary grid line orthogonality parameters. */
      
      Stretch_I = STRETCHING_FCN_LINEAR;
      Beta_I = ZERO; 
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_LINEAR;
      Beta_J = ZERO;
      Tau_J = ZERO;
      Orthogonal_North = 0;
      Orthogonal_South = 0;
      Orthogonal_East = 0;
      Orthogonal_West = 0;
      
      /* Create the 2D quadrilateral grid block representing
	 the mesh. */
      
      Grid_ptr[iBlk][jBlk].Create_Quad_Block_Without_Update(Bnd_Spline_North,
							    Bnd_Spline_South,
							    Bnd_Spline_East,
							    Bnd_Spline_West,
							    Number_of_Cells_Idir/Number_of_Blocks_Idir,
							    Number_of_Cells_Jdir/Number_of_Blocks_Jdir,
							    Number_of_Ghost_Cells,
							    Highest_Order_of_Reconstruction,
							    GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
							    Stretch_I,
							    Beta_I, 
							    Tau_I,
							    Stretch_J,
							    Beta_J,
							    Tau_J,
							    Orthogonal_North,
							    Orthogonal_South,
							    Orthogonal_East,
							    Orthogonal_West);
      
    } /* endfor */
  } /* endfor */

}

void Grid2D_Quad_MultiBlock_HO::Grid_Rectangular_Box_Without_Update(int &_Number_of_Blocks_Idir_,
								    int &_Number_of_Blocks_Jdir_,
								    const double &Width,
								    const double &Height,
								    const int Stretching_Flag,
								    const int Stretching_Type_Idir,
								    const int Stretching_Type_Jdir,
								    const double &Stretching_Factor_Idir,
								    const double &Stretching_Factor_Jdir,
								    const int Number_of_Cells_Idir,
								    const int Number_of_Cells_Jdir,
								    const int Number_of_Ghost_Cells,
								    const int Highest_Order_of_Reconstruction) {

  int iBlk, jBlk, n_cells_i, n_cells_j,
      Stretch_I, Stretch_J,
      Orthogonal_North, Orthogonal_South,
      Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D xc_NW, xc_NE, xc_SE, xc_SW;
  Spline2D_HO Bnd_Spline_North, Bnd_Spline_South,
              Bnd_Spline_East, Bnd_Spline_West;

  /* Allocate memory for grid block. */

  if (_Number_of_Blocks_Idir_ < 0) _Number_of_Blocks_Idir_ = 1;
  if (_Number_of_Blocks_Jdir_ < 0) _Number_of_Blocks_Jdir_ = 1;
  allocate(_Number_of_Blocks_Idir_, _Number_of_Blocks_Jdir_);
  /* Create the mesh for each block representing
     the complete grid. */

  for ( jBlk = 0; jBlk <= Number_of_Blocks_Jdir-1; ++jBlk ) {
    for ( iBlk = 0; iBlk <= Number_of_Blocks_Idir-1; ++iBlk ) {

      /* Assign values to the locations of the corners
	 of the rectangular box shaped domain. */

      xc_NW = Vector2D(-HALF*Width+double(iBlk)*Width/double(Number_of_Blocks_Idir), 
		       -HALF*Height+double(jBlk+1)*Height/double(Number_of_Blocks_Jdir));
      xc_NE = Vector2D(-HALF*Width+double(iBlk+1)*Width/double(Number_of_Blocks_Idir), 
		       -HALF*Height+double(jBlk+1)*Height/double(Number_of_Blocks_Jdir));
      xc_SE = Vector2D(-HALF*Width+double(iBlk+1)*Width/double(Number_of_Blocks_Idir), 
		       -HALF*Height+double(jBlk)*Height/double(Number_of_Blocks_Jdir));
      xc_SW = Vector2D(-HALF*Width+double(iBlk)*Width/double(Number_of_Blocks_Idir), 
		       -HALF*Height+double(jBlk)*Height/double(Number_of_Blocks_Jdir));

      /* Create the splines defining the north, south,
	 east, and west boundaries of the rectangular box. */
      
      Bnd_Spline_North.Create_Spline_Line(xc_NW, xc_NE, 2);
      Bnd_Spline_South.Create_Spline_Line(xc_SW, xc_SE, 2);
      Bnd_Spline_East.Create_Spline_Line(xc_SE, xc_NE, 2);
      Bnd_Spline_West.Create_Spline_Line(xc_SW, xc_NW, 2);

      /* Set the boundary condition types for each of the
	 boundary splines. */

      if (jBlk == Number_of_Blocks_Jdir-1) {
	Bnd_Spline_North.setBCtype(BC_REFLECTION);
      } else {
	Bnd_Spline_North.setBCtype(BC_NONE);
      } /* endif */
      if (jBlk == 0) {
	Bnd_Spline_South.setBCtype(BC_REFLECTION);
      } else {
	Bnd_Spline_South.setBCtype(BC_NONE);
      } /* endif */
      if (iBlk == Number_of_Blocks_Idir-1) {
	Bnd_Spline_East.setBCtype(BC_REFLECTION);
      } else {
	Bnd_Spline_East.setBCtype(BC_NONE);
      } /* endif */
      if (iBlk == 0) {
	Bnd_Spline_West.setBCtype(BC_REFLECTION);
      } else {
	Bnd_Spline_West.setBCtype(BC_NONE);
      } /* endif */

      /* Assign values to the stretching function parameters
	 and boundary grid line orthogonality parameters. */

      if (Stretching_Flag) {
	Stretch_I = Stretching_Type_Idir;
	Stretch_J = Stretching_Type_Jdir;
	Beta_I = Stretching_Factor_Idir;
	Beta_J = Stretching_Factor_Jdir;
      } else {
	Stretch_I = STRETCHING_FCN_LINEAR;
	Stretch_J = STRETCHING_FCN_LINEAR;
	Beta_I = ZERO; 
	Beta_J = ZERO;
      }
      Tau_I = ZERO;
      Tau_J = ZERO;
      Orthogonal_North = 0;
      Orthogonal_South = 0;
      Orthogonal_East = 0;
      Orthogonal_West = 0;

      /* Create the 2D quadrilateral grid block representing
	 the mesh. */

      Grid_ptr[iBlk][jBlk].Create_Quad_Block_Without_Update(Bnd_Spline_North,
							    Bnd_Spline_South,
							    Bnd_Spline_East,
							    Bnd_Spline_West,
							    Number_of_Cells_Idir/Number_of_Blocks_Idir,
							    Number_of_Cells_Jdir/Number_of_Blocks_Jdir,
							    Number_of_Ghost_Cells,
							    Highest_Order_of_Reconstruction,
							    GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
							    Stretch_I,
							    Beta_I, 
							    Tau_I,
							    Stretch_J,
							    Beta_J,
							    Tau_J,
							    Orthogonal_North,
							    Orthogonal_South,
							    Orthogonal_East,
							    Orthogonal_West);

    } /* endfor */
  } /* endfor */


}

/*!
 * Generates a uniform 2D mesh for a deformed
 * box shaped domain. 
 *
 * \param VertexSW the SW corner of the mesh
 * \param VertexSE the SE corner of the mesh
 * \param VertexNE the NE corner of the mesh
 * \param VertexNW the NW corner of the mesh
 *                                                      
 * This subroutine DOESN'T update the ghost cells or
 * the geometric properties of the grid cells.
 */
void Grid2D_Quad_MultiBlock_HO::Grid_Deformed_Box_Without_Update(int &_Number_of_Blocks_Idir_,
								 int &_Number_of_Blocks_Jdir_,
								 const Vector2D &VertexSW,
								 const Vector2D &VertexSE,
								 const Vector2D &VertexNE,
								 const Vector2D &VertexNW,     
								 const int Stretching_Flag,
								 const int Stretching_Type_Idir,
								 const int Stretching_Type_Jdir,
								 const double &Stretching_Factor_Idir,
								 const double &Stretching_Factor_Jdir,
								 const int Number_of_Cells_Idir,
								 const int Number_of_Cells_Jdir,
								 const int Number_of_Ghost_Cells,
								 const int Highest_Order_of_Reconstruction){

  int iBlk, jBlk, n_cells_i, n_cells_j,
      Stretch_I, Stretch_J,
      Orthogonal_North, Orthogonal_South,
      Orthogonal_East, Orthogonal_West;
  int INl, JNl;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D xc_NW, xc_NE, xc_SE, xc_SW;
  Spline2D_HO Bnd_Spline_North, Bnd_Spline_South,
              Bnd_Spline_East, Bnd_Spline_West;
  int InfoQuad, InfoBlocks;

  Grid2D_Quad_Block_HO CornersMultiBlockMesh;

  // Check if it's possible to form a valid mesh with the given vertices.
  
  InfoQuad = Find_Quadrilateral_Type(VertexSW,VertexSE,VertexNE,VertexNW);

  if (InfoQuad == 0 || InfoQuad == 4){
    // The vertices form a degenerated or crossed quadrilateral
    return;
  }

  /* Allocate memory for grid block. */

  if (_Number_of_Blocks_Idir_ < 0) _Number_of_Blocks_Idir_ = 1;
  if (_Number_of_Blocks_Jdir_ < 0) _Number_of_Blocks_Jdir_ = 1;

  /* Check consistency of the number of blocks. */
  if ( _Number_of_Blocks_Idir_ == 1 && _Number_of_Blocks_Jdir_ == 1 ){
    // 1x1 blocks
    InfoBlocks = 1;
  } else {
    // nxm blocks (n and m are even numbers)
    _Number_of_Blocks_Idir_ += _Number_of_Blocks_Idir_ % 2;
    _Number_of_Blocks_Jdir_ += _Number_of_Blocks_Jdir_ % 2;
    InfoBlocks = 2;
  }
  allocate(_Number_of_Blocks_Idir_, _Number_of_Blocks_Jdir_);


  /* Assign values to the stretching function parameters
     and boundary grid line orthogonality parameters. */
  
  if (Stretching_Flag) {
    Stretch_I = Stretching_Type_Idir;
    Stretch_J = Stretching_Type_Jdir;
    Beta_I = Stretching_Factor_Idir;
    Beta_J = Stretching_Factor_Jdir;
  } else {
    Stretch_I = STRETCHING_FCN_LINEAR;
    Stretch_J = STRETCHING_FCN_LINEAR;
    Beta_I = ZERO; 
    Beta_J = ZERO;
  }
  Tau_I = ZERO;
  Tau_J = ZERO;
  Orthogonal_North = 0;
  Orthogonal_South = 0;
  Orthogonal_East = 0;
  Orthogonal_West = 0;

  // Create the boundaries of the multiblock mesh.
  Bnd_Spline_North.Create_Spline_Line(VertexNW,VertexNE,2);
  Bnd_Spline_South.Create_Spline_Line(VertexSW,VertexSE,2);
  Bnd_Spline_East.Create_Spline_Line(VertexSE,VertexNE,2);
  Bnd_Spline_West.Create_Spline_Line(VertexSW,VertexNW,2);

  if (InfoBlocks == 1){
    /* Set the boundary condition types for each of the
       boundary splines. */
    Bnd_Spline_North.setBCtype(BC_REFLECTION);
    Bnd_Spline_South.setBCtype(BC_REFLECTION);
    Bnd_Spline_East.setBCtype(BC_REFLECTION);
    Bnd_Spline_West.setBCtype(BC_REFLECTION);
    
    /* Create the 2D quadrilateral grid block representing
       the mesh. */
    
    Grid_ptr[0][0].Create_Quad_Block_Without_Update(Bnd_Spline_North,
						    Bnd_Spline_South,
						    Bnd_Spline_East,
						    Bnd_Spline_West,
						    Number_of_Cells_Idir,
						    Number_of_Cells_Jdir,
						    Number_of_Ghost_Cells,
						    Highest_Order_of_Reconstruction,
						    GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
						    Stretch_I,
						    Beta_I, 
						    Tau_I,
						    Stretch_J,
						    Beta_J,
						    Tau_J,
						    Orthogonal_North,
						    Orthogonal_South,
						    Orthogonal_East,
						    Orthogonal_West);
    
  } else {

    // Obtain the corner locations of the multiblock mesh.
    CornersMultiBlockMesh.Create_Quad_Block_Without_Update(Bnd_Spline_North,
							   Bnd_Spline_South,
							   Bnd_Spline_East,
							   Bnd_Spline_West,
							   Number_of_Blocks_Idir,
							   Number_of_Blocks_Jdir,
							   2,
							   0,
							   GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
							   STRETCHING_FCN_LINEAR,
							   0, 
							   0,
							   STRETCHING_FCN_LINEAR,
							   0, 
							   0,
							   0,
							   0,
							   0,
							   0);

    INl = CornersMultiBlockMesh.INl;
    JNl = CornersMultiBlockMesh.JNl;

    /* Create the mesh for each block representing
       the complete grid. */
    
    for ( jBlk = 0; jBlk <= Number_of_Blocks_Jdir-1; ++jBlk ) {
      for ( iBlk = 0; iBlk <= Number_of_Blocks_Idir-1; ++iBlk ) {

	/* Assign values to the locations of the corners
	   of the rectangular box shaped domain. */

	xc_SW = CornersMultiBlockMesh.Node[iBlk+INl  ][jBlk+JNl  ].X;
	xc_SE = CornersMultiBlockMesh.Node[iBlk+INl+1][jBlk+JNl  ].X;
	xc_NE = CornersMultiBlockMesh.Node[iBlk+INl+1][jBlk+JNl+1].X;
	xc_NW = CornersMultiBlockMesh.Node[iBlk+INl  ][jBlk+JNl+1].X;

	/* Create the splines defining the north, south,
	   east, and west boundaries of the rectangular box. */
      
	Bnd_Spline_North.Create_Spline_Line(xc_NW, xc_NE, 2);
	Bnd_Spline_South.Create_Spline_Line(xc_SW, xc_SE, 2);
	Bnd_Spline_East.Create_Spline_Line(xc_SE, xc_NE, 2);
	Bnd_Spline_West.Create_Spline_Line(xc_SW, xc_NW, 2);

	/* Set the boundary condition types for each of the
	   boundary splines. */

	if (jBlk == Number_of_Blocks_Jdir-1) {
	  Bnd_Spline_North.setBCtype(BC_REFLECTION);
	} else {
	  Bnd_Spline_North.setBCtype(BC_NONE);
	} /* endif */
	if (jBlk == 0) {
	  Bnd_Spline_South.setBCtype(BC_REFLECTION);
	} else {
	  Bnd_Spline_South.setBCtype(BC_NONE);
	} /* endif */
	if (iBlk == Number_of_Blocks_Idir-1) {
	  Bnd_Spline_East.setBCtype(BC_REFLECTION);
	} else {
	  Bnd_Spline_East.setBCtype(BC_NONE);
	} /* endif */
	if (iBlk == 0) {
	  Bnd_Spline_West.setBCtype(BC_REFLECTION);
	} else {
	  Bnd_Spline_West.setBCtype(BC_NONE);
	} /* endif */


	/* Create the 2D quadrilateral grid block representing
	   the mesh. */

	Grid_ptr[iBlk][jBlk].Create_Quad_Block_Without_Update(Bnd_Spline_North,
							      Bnd_Spline_South,
							      Bnd_Spline_East,
							      Bnd_Spline_West,
							      Number_of_Cells_Idir/Number_of_Blocks_Idir,
							      Number_of_Cells_Jdir/Number_of_Blocks_Jdir,
							      Number_of_Ghost_Cells,
							      Highest_Order_of_Reconstruction,
							      GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
							      Stretch_I,
							      Beta_I, 
							      Tau_I,
							      Stretch_J,
							      Beta_J,
							      Tau_J,
							      Orthogonal_North,
							      Orthogonal_South,
							      Orthogonal_East,
							      Orthogonal_West);
      } /* endfor */
    } /* endfor */
  }// endif 

}

/*!
 * Generates a quadilateral mesh with clustering        
 * consisting of two grid blocks for predicting viscous 
 * flow and boundary layer development over a flat      
 * plate.                                               
 *                                                      
 * Usage: Grid_ptr = Grid_Flat_Plate(Grid_ptr,          
 *                                   nblk_i,            
 *                                   nblk_j,            
 *                                   TWO,               
 *         		             100,               
 *         		             100,               
 *         		             2);                
 *             
 * This subroutine DOESN'T update the ghost cells or
 * the geometric properties of the grid cells.                                         
 */
void Grid2D_Quad_MultiBlock_HO::Grid_Flat_Plate_Without_Update(int &_Number_of_Blocks_Idir_,
							       int &_Number_of_Blocks_Jdir_,
							       const double &Length,
							       const int Flat_Plate_BC_Type,
							       const int Stretching_Flag,
							       const double &Stretching_Factor_Idir,
							       const double &Stretching_Factor_Jdir,
							       const int Number_of_Cells_Idir,
							       const int Number_of_Cells_Jdir,
							       const int Number_of_Ghost_Cells,
							       const int Highest_Order_of_Reconstruction) {
  
  int n_cells_i, n_cells_j, 
      Stretch_I, Stretch_J,
      Orthogonal_North, Orthogonal_South,
      Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D xc_NW, xc_NE, xc_SE, xc_SW;
  Spline2D_HO Bnd_Spline_North, Bnd_Spline_South,
              Bnd_Spline_East, Bnd_Spline_West;

  // Allocate memory for grid blocks.  There are two grid blocks for 
  // this mesh.
  if (Flat_Plate_BC_Type == BC_BURNING_SURFACE){
    _Number_of_Blocks_Idir_ = 3;
  } else {
    _Number_of_Blocks_Idir_ = 2;
  }
  _Number_of_Blocks_Jdir_ = 1;
  allocate(_Number_of_Blocks_Idir_,_Number_of_Blocks_Jdir_);

  // Create the mesh for each block representing the complete grid.
  for (int iBlk = 0; iBlk < Number_of_Blocks_Idir; iBlk++) {

    // Assign values to the locations of the corners of the 
    // rectangular box shaped domain representing each of the 
    // blocks in the grid.
    if (iBlk == 0) {
      xc_NW = TWO*Vector2D(-Length,Length);
      xc_NE = TWO*Vector2D(ZERO,Length);
      xc_SE = TWO*Vector2D(ZERO,ZERO);
      xc_SW = TWO*Vector2D(-Length,ZERO);
    } else if (iBlk == 1) {
      xc_NW = TWO*Vector2D(ZERO,Length);
      xc_NE = TWO*Vector2D(Length,Length);
      xc_SE = TWO*Vector2D(Length,ZERO);
      xc_SW = TWO*Vector2D(ZERO,  ZERO);
    } else {
      xc_NW = TWO*Vector2D(Length,Length);
      xc_NE = TWO*Vector2D(TWO*Length,Length);
      xc_SE = TWO*Vector2D(TWO*Length,ZERO);
      xc_SW = TWO*Vector2D(Length,ZERO);
    }

    // Create the splines defining the north, south, east, and west 
    // boundaries of the rectangular boxes.
    Bnd_Spline_North.Create_Spline_Line(xc_NW,xc_NE,2);
    Bnd_Spline_South.Create_Spline_Line(xc_SW,xc_SE,2);
    Bnd_Spline_East.Create_Spline_Line(xc_SE,xc_NE,2);
    Bnd_Spline_West.Create_Spline_Line(xc_SW,xc_NW,2);

    // Set the boundary condition types for each of the boundary splines.
    if (iBlk == 0) {
      Bnd_Spline_North.setBCtype(BC_FIXED);
      Bnd_Spline_South.setBCtype(BC_REFLECTION);
      Bnd_Spline_East.setBCtype(BC_NONE);
      Bnd_Spline_West.setBCtype(BC_FIXED);
    } else if (iBlk == 1) {
      Bnd_Spline_North.setBCtype(BC_FIXED);
      if (Flat_Plate_BC_Type != BC_BURNING_SURFACE) {
	Bnd_Spline_South.setBCtype(Flat_Plate_BC_Type);
	Bnd_Spline_East.setBCtype(BC_CONSTANT_EXTRAPOLATION);
      } else {
	Bnd_Spline_South.setBCtype(BC_WALL_VISCOUS_ISOTHERMAL);
	Bnd_Spline_East.setBCtype(BC_NONE);
      }
      Bnd_Spline_West.setBCtype(BC_NONE);
    } else {
      Bnd_Spline_North.setBCtype(BC_FIXED);
      Bnd_Spline_South.setBCtype(Flat_Plate_BC_Type);
      Bnd_Spline_East.setBCtype(BC_CONSTANT_EXTRAPOLATION);
      Bnd_Spline_West.setBCtype(BC_NONE);
    }

    // Determine the number of cells for this block.
    n_cells_i = Number_of_Cells_Idir;///2;
    n_cells_j = Number_of_Cells_Jdir;

    // Assign values to the stretching function parameters and boundary
    // grid line orthogonality parameters.
    Stretch_I = STRETCHING_FCN_LINEAR;
    Beta_I = ZERO;
    Tau_I = ZERO;
    Stretch_J = STRETCHING_FCN_LINEAR;
    Beta_J = ZERO;
    Tau_J = ZERO;
    if (Stretching_Flag) {
      if (iBlk == 0) {
	Stretch_I = STRETCHING_FCN_MAX_CLUSTERING;
      } else if (iBlk == 1) {
	if (Flat_Plate_BC_Type != BC_BURNING_SURFACE) {
	  Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
	} else {
	  Stretch_I = STRETCHING_FCN_MINMAX_CLUSTERING;
	}
      } else {
	Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
      }
      Beta_I = Stretching_Factor_Idir;
      Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
      Beta_J = Stretching_Factor_Jdir;
    }
    Orthogonal_North = 0;
    Orthogonal_South = 0;
    Orthogonal_East = 0;
    Orthogonal_West = 0;

    // Create the 2D quadrilateral grid block.
    Grid_ptr[iBlk][0].Create_Quad_Block_Without_Update(Bnd_Spline_North,
						       Bnd_Spline_South,
						       Bnd_Spline_East,
						       Bnd_Spline_West,
						       n_cells_i,
						       n_cells_j,
						       Number_of_Ghost_Cells,
						       Highest_Order_of_Reconstruction,
						       GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
						       Stretch_I,
						       Beta_I,
						       Tau_I,
						       Stretch_J,
						       Beta_J,
						       Tau_J,
						       Orthogonal_North,
						       Orthogonal_South,
						       Orthogonal_East,
						       Orthogonal_West);

  }
}


/*!
 * Generates a quadilateral mesh with clustering        
 * consisting of two grid blocks for predicting viscous 
 * flow and boundary layer development over a flat      
 * plate.                                               
 *                                                      
 * Usage: Grid_ptr = Grid_Flat_Plate_NK(Grid_ptr,       
 *                                      nblk_i,         
 *                                      nblk_j,         
 *                                      TWO,            
 *         		                100,            
 *         		                100);           
 *             
 * This subroutine DOESN'T update the ghost cells or
 * the geometric properties of the grid cells.                                         
 */
void Grid2D_Quad_MultiBlock_HO::Grid_Flat_Plate_NK_Without_Update(int &_Number_of_Blocks_Idir_,
								  int &_Number_of_Blocks_Jdir_,
								  const double &Length,
								  const int Stretching_Flag,
								  const double &Stretching_Factor_Idir,
								  const double &Stretching_Factor_Jdir,
								  const int Number_of_Cells_Idir,
								  const int Number_of_Cells_Jdir,
								  const int Number_of_Ghost_Cells,
								  const int Highest_Order_of_Reconstruction) {
  
  int iBlk, n_cells_i, n_cells_j, 
    Stretch_I, Stretch_J,
    Orthogonal_North, Orthogonal_South,
    Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D xc_NW, xc_NE, xc_SE, xc_SW;
  Spline2D_HO Bnd_Spline_North, Bnd_Spline_South,
              Bnd_Spline_East, Bnd_Spline_West;
  
  /* Allocate memory for grid blocks.  There are three grid blocks for this mesh. */
  _Number_of_Blocks_Idir_ = 3;
  _Number_of_Blocks_Jdir_ = 1;
  allocate(_Number_of_Blocks_Idir_,_Number_of_Blocks_Jdir_);
  
  /* Create the mesh for each block representing  the complete grid. */
  for ( iBlk = 0; iBlk < Number_of_Blocks_Idir; ++iBlk ) {
    
    /* Assign values to the locations of the corners
       of the rectangular box shaped domain representing
       each of the blocks in the grid. */
    if (iBlk == 0) {
      xc_NW = FOUR*Vector2D(-Length,Length);
      xc_NE = FOUR*Vector2D(ZERO,Length);
      xc_SE = FOUR*Vector2D(ZERO,ZERO);
      xc_SW = FOUR*Vector2D(-Length,ZERO);
    } else if( iBlk == 1) {
	xc_NW = Vector2D(ZERO,FOUR*Length);
	xc_NE = Vector2D(Length,FOUR*Length);
	xc_SE = Vector2D(Length,ZERO);
	xc_SW = Vector2D(ZERO,  ZERO);
      } else {
	xc_NW = Vector2D(Length,FOUR*Length);
	xc_NE = Vector2D(FIVE*Length,FOUR*Length);
	xc_SE = Vector2D(FIVE*Length,ZERO);
	xc_SW = Vector2D(Length,ZERO);	
      }

      /* Create the splines defining the north, south,
	 east, and west boundaries of the rectangular boxes. */
      
      Bnd_Spline_North.Create_Spline_Line(xc_NW, xc_NE, 2);
      Bnd_Spline_South.Create_Spline_Line(xc_SW, xc_SE, 2);
      Bnd_Spline_East.Create_Spline_Line(xc_SE, xc_NE, 2);
      Bnd_Spline_West.Create_Spline_Line(xc_SW, xc_NW, 2);
      
      /* Set the boundary condition types for each of the
	 boundary splines. */
      
      if (iBlk == 0) {
	Bnd_Spline_North.setBCtype(BC_CHARACTERISTIC); //BC_CONSTANT_EXTRAPOLATION);
	Bnd_Spline_South.setBCtype(BC_REFLECTION);
	Bnd_Spline_East.setBCtype(BC_NONE);
	Bnd_Spline_West.setBCtype(BC_FIXED);
      } else if (iBlk== 1) {
	Bnd_Spline_North.setBCtype(BC_CHARACTERISTIC); //BC_CONSTANT_EXTRAPOLATION);
	Bnd_Spline_South.setBCtype(BC_WALL_VISCOUS_HEATFLUX);  //BC_FIXED_TEMP_WALL);
	Bnd_Spline_East.setBCtype(BC_NONE); 
	Bnd_Spline_West.setBCtype(BC_NONE);
      } else {
	Bnd_Spline_North.setBCtype(BC_CHARACTERISTIC); //BC_CONSTANT_EXTRAPOLATION);
	Bnd_Spline_South.setBCtype(BC_REFLECTION); 
	Bnd_Spline_East.setBCtype(BC_CHARACTERISTIC); //BC_CONSTANT_EXTRAPOLATION); 
	Bnd_Spline_West.setBCtype(BC_NONE);
      } 

      /* Determine the number of cells for this block. */
      n_cells_i = Number_of_Cells_Idir;
      n_cells_j = Number_of_Cells_Jdir;

        /* Assign values to the stretching function parameters
           and boundary grid line orthogonality parameters. */

	Stretch_I = STRETCHING_FCN_LINEAR;
	Beta_I = ZERO;
	Tau_I = ZERO;
	Stretch_J = STRETCHING_FCN_LINEAR;
	Beta_J = ZERO;
	Tau_J = ZERO;
	//     if (Stretching_Flag) {
	//       if (iBlk == 0) {
	// 	Stretch_I = STRETCHING_FCN_MAX_CLUSTERING;
	// 	Beta_I = ONE + (Stretching_Factor_Idir-ONE)/6.50;
	//       } else {
	// 	Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
	// 	Beta_I = Stretching_Factor_Idir;
	//       }
	//       Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
	//       Beta_J = Stretching_Factor_Jdir;
	//     }
	if (Stretching_Flag) {
	  if (iBlk == 0) Stretch_I = STRETCHING_FCN_MAX_CLUSTERING;
	  else Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
	  Beta_I = Stretching_Factor_Idir;
	  Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
	  Beta_J = Stretching_Factor_Jdir;
	}
	Orthogonal_North = 0;
	Orthogonal_South = 0;
	Orthogonal_East = 0;
	Orthogonal_West = 0;
	
        /* Create the 2D quadrilateral grid block. */

        Grid_ptr[iBlk][0].Create_Quad_Block_Without_Update(Bnd_Spline_North,
							   Bnd_Spline_South,
							   Bnd_Spline_East,
							   Bnd_Spline_West,
							   n_cells_i,
							   n_cells_j,
							   Number_of_Ghost_Cells,
							   Highest_Order_of_Reconstruction,
							   GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
							   Stretch_I,
							   Beta_I, 
							   Tau_I,
							   Stretch_J,
							   Beta_J,
							   Tau_J,
							   Orthogonal_North,
							   Orthogonal_South,
							   Orthogonal_East,
							   Orthogonal_West);
	
    } /* endfor */    


}

/*!
 * Generates a quadilateral mesh with clustering        
 * consisting of two grid blocks for predicting viscous 
 * flow and boundary layer development over a flat      
 * plate.                                               
 *                                                      
 * Usage: Grid_ptr = Grid_Flat_Plate3(Grid_ptr,         
 *                                   nblk_i,            
 *                                   nblk_j,            
 *                                   TWO,               
 *         		             100,               
 *         		             100,               
 *         		             2);                
 *                                                      
 * This subroutine DOESN'T update the ghost cells or
 * the geometric properties of the grid cells.                                         
 */
void Grid2D_Quad_MultiBlock_HO::Grid_Flat_Plate3_Without_Update(int &_Number_of_Blocks_Idir_,
								int &_Number_of_Blocks_Jdir_,
								const double &Length,
								const int &Stretching_Flag,
								const double &Stretching_Factor_Idir,
								const double &Stretching_Factor_Jdir,
								const int Number_of_Cells_Idir,
								const int Number_of_Cells_Jdir,
								const int Number_of_Ghost_Cells,
								const int Highest_Order_of_Reconstruction) {

  int n_cells_i, n_cells_j, 
      Stretch_I, Stretch_J,
      Orthogonal_North, Orthogonal_South,
      Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D xc_NW, xc_NE, xc_SE, xc_SW;
  Spline2D_HO Bnd_Spline_North, Bnd_Spline_South,
              Bnd_Spline_East, Bnd_Spline_West;

  // Allocate memory for grid blocks.  There are two grid blocks for 
  // this mesh.
  _Number_of_Blocks_Idir_ = 2;
  _Number_of_Blocks_Jdir_ = 1;
  allocate(_Number_of_Blocks_Idir_,_Number_of_Blocks_Jdir_);

  // Create the mesh for each block representing the complete grid.
  for (int iBlk = 0; iBlk < Number_of_Blocks_Idir; iBlk++) {

    // Assign values to the locations of the corners of the 
    // rectangular box shaped domain representing each of the 
    // blocks in the grid.
    if (iBlk == 0) {
      xc_NW = FOUR*Vector2D(-Length,Length);
      xc_NE = FOUR*Vector2D(ZERO,Length);
      xc_SE = FOUR*Vector2D(ZERO,ZERO);
      xc_SW = FOUR*Vector2D(-Length,ZERO);
    } else {
      xc_NW = Vector2D(ZERO,FOUR*Length);
      xc_NE = Vector2D(Length,FOUR*Length);
      xc_SE = Vector2D(Length,ZERO);
      xc_SW = Vector2D(ZERO,  ZERO);
    }

    // Create the splines defining the north, south, east, and west 
    // boundaries of the rectangular boxes.
    Bnd_Spline_North.Create_Spline_Line(xc_NW,xc_NE,2);
    Bnd_Spline_South.Create_Spline_Line(xc_SW,xc_SE,2);
    Bnd_Spline_East.Create_Spline_Line(xc_SE,xc_NE,2);
    Bnd_Spline_West.Create_Spline_Line(xc_SW,xc_NW,2);

    // Set the boundary condition types for each of the boundary splines.
    if (iBlk == 0) {
      Bnd_Spline_North.setBCtype(BC_FIXED);
      Bnd_Spline_South.setBCtype(BC_REFLECTION);
      Bnd_Spline_East.setBCtype(BC_NONE);
      Bnd_Spline_West.setBCtype(BC_FIXED);
    } else {
      Bnd_Spline_North.setBCtype(BC_FIXED);
      Bnd_Spline_South.setBCtype(BC_WALL_VISCOUS_HEATFLUX);
      Bnd_Spline_East.setBCtype(BC_CONSTANT_EXTRAPOLATION);
      Bnd_Spline_West.setBCtype(BC_NONE);
    }

    // Determine the number of cells for this block.
    n_cells_i = Number_of_Cells_Idir;///2;
    n_cells_j = Number_of_Cells_Jdir;

    // Assign values to the stretching function parameters and boundary
    // grid line orthogonality parameters.
    Stretch_I = STRETCHING_FCN_LINEAR;
    Beta_I = ZERO;
    Tau_I = ZERO;
    Stretch_J = STRETCHING_FCN_LINEAR;
    Beta_J = ZERO;
    Tau_J = ZERO;
    if (Stretching_Flag) {
      if (iBlk == 0) {
	Stretch_I = STRETCHING_FCN_MAX_CLUSTERING;
	Beta_I = ONE + (Stretching_Factor_Idir-ONE)/6.50;
      } else {
	Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_I = Stretching_Factor_Idir;
      }
      Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
      Beta_J = Stretching_Factor_Jdir;
    }
    Orthogonal_North = 0;
    Orthogonal_South = 0;
    Orthogonal_East = 0;
    Orthogonal_West = 0;

    // Create the 2D quadrilateral grid block.
    Grid_ptr[iBlk][0].Create_Quad_Block_Without_Update(Bnd_Spline_North,
						       Bnd_Spline_South,
						       Bnd_Spline_East,
						       Bnd_Spline_West,
						       n_cells_i,
						       n_cells_j,
						       Number_of_Ghost_Cells,
						       Highest_Order_of_Reconstruction,
						       GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
						       Stretch_I,
						       Beta_I,
						       Tau_I,
						       Stretch_J,
						       Beta_J,
						       Tau_J,
						       Orthogonal_North,
						       Orthogonal_South,
						       Orthogonal_East,
						       Orthogonal_West);
    
  }

}

/*!
 * Generates a quadilateral mesh with clustering        
 * consisting of two grid blocks for predicting viscous 
 * flow and boundary layer development over a flat      
 * plate.                                               
 *                                                      
 * Usage: Grid_ptr = Grid_Flat_Plate4(Grid_ptr,         
 *                                   nblk_i,            
 *                                   nblk_j,            
 *                                   TWO,               
 *         		             100,               
 *         		             100,               
 *         		             2);                
 * This subroutine DOESN'T update the ghost cells or
 * the geometric properties of the grid cells.
 */
void Grid2D_Quad_MultiBlock_HO::Grid_Flat_Plate4_Without_Update(int &_Number_of_Blocks_Idir_,
								int &_Number_of_Blocks_Jdir_,
								const double &Length,
								const int &Stretching_Flag,
								const double &Stretching_Factor_Idir,
								const double &Stretching_Factor_Jdir,
								const int Number_of_Cells_Idir,
								const int Number_of_Cells_Jdir,
								const int Number_of_Ghost_Cells,
								const int Highest_Order_of_Reconstruction) {
  
  int n_cells_i, n_cells_j, 
      Stretch_I, Stretch_J,
      Orthogonal_North, Orthogonal_South,
      Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D xc_NW, xc_NE, xc_SE, xc_SW;
  Spline2D_HO Bnd_Spline_North, Bnd_Spline_South,
              Bnd_Spline_East, Bnd_Spline_West;

  // Allocate memory for grid blocks.  There are two grid blocks for 
  // this mesh.
  _Number_of_Blocks_Idir_ = 3;
  _Number_of_Blocks_Jdir_ = 1;
  allocate(_Number_of_Blocks_Idir_,_Number_of_Blocks_Jdir_);

  // Create the mesh for each block representing the complete grid.
  for (int iBlk = 0; iBlk < Number_of_Blocks_Idir; iBlk++) {

    // Assign values to the locations of the corners of the 
    // rectangular box shaped domain representing each of the 
    // blocks in the grid.
    if (iBlk == 0) {
      xc_NW = Vector2D(-4.0*Length,4.0*Length);
      xc_NE = Vector2D(ZERO,4.0*Length);
      xc_SE = Vector2D(ZERO,ZERO);
      xc_SW = Vector2D(-4.0*Length,ZERO);
    } else if (iBlk == 1) {
      xc_NW = Vector2D(ZERO,4.0*Length);
      xc_NE = Vector2D(Length,4.0*Length);
      xc_SE = Vector2D(Length,ZERO);
      xc_SW = Vector2D(ZERO,ZERO);
    } else {
      xc_NW = Vector2D(Length,4.0*Length);
      xc_NE = Vector2D(5.0*Length,4.0*Length);
      xc_SE = Vector2D(5.0*Length,ZERO);
      xc_SW = Vector2D(Length,ZERO);
    }

    // Create the splines defining the north, south, east, and west 
    // boundaries of the rectangular boxes.
    Bnd_Spline_North.Create_Spline_Line(xc_NW,xc_NE,2);
    Bnd_Spline_South.Create_Spline_Line(xc_SW,xc_SE,2);
    Bnd_Spline_East.Create_Spline_Line(xc_SE,xc_NE,2);
    Bnd_Spline_West.Create_Spline_Line(xc_SW,xc_NW,2);

    // Set the boundary condition types for each of the boundary splines.
    if (iBlk == 0) {
      Bnd_Spline_North.setBCtype(BC_FIXED);//BC_CONSTANT_EXTRAPOLATION);
      Bnd_Spline_South.setBCtype(BC_REFLECTION);
      Bnd_Spline_East.setBCtype(BC_NONE);
      Bnd_Spline_West.setBCtype(BC_FIXED);
    } else if (iBlk == 1) {
      Bnd_Spline_North.setBCtype(BC_FIXED);//BC_CONSTANT_EXTRAPOLATION);
      Bnd_Spline_South.setBCtype(BC_WALL_VISCOUS_HEATFLUX);
      Bnd_Spline_East.setBCtype(BC_NONE);
      Bnd_Spline_West.setBCtype(BC_NONE);
    } else {
      Bnd_Spline_North.setBCtype(BC_FIXED);//BC_CONSTANT_EXTRAPOLATION);
      Bnd_Spline_South.setBCtype(BC_REFLECTION);
      Bnd_Spline_East.setBCtype(BC_FIXED);//BC_CONSTANT_EXTRAPOLATION);
      Bnd_Spline_West.setBCtype(BC_NONE);
    }

    // Determine the number of cells for this block.
    n_cells_i = Number_of_Cells_Idir;
    n_cells_j = Number_of_Cells_Jdir;

    // Assign values to the stretching function parameters and boundary
    // grid line orthogonality parameters.
    if (iBlk == 0) {
      Stretch_I = STRETCHING_FCN_MAX_CLUSTERING;
      Beta_I = ONE + (Stretching_Factor_Idir-ONE)/6.50;
    } else if (iBlk == 1) {
      Stretch_I = STRETCHING_FCN_MINMAX_CLUSTERING;
      Beta_I = Stretching_Factor_Idir;
    } else {
      Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
      Beta_I = ONE + (Stretching_Factor_Idir-ONE)/6.50;
    }
    Tau_I = ZERO;
    Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
    Beta_J = Stretching_Factor_Jdir;
    Tau_J = ZERO;
    Orthogonal_North = 0;
    Orthogonal_South = 0;
    Orthogonal_East = 0;
    Orthogonal_West = 0;

    // Create the 2D quadrilateral grid block.
    Grid_ptr[iBlk][0].Create_Quad_Block_Without_Update(Bnd_Spline_North,
						       Bnd_Spline_South,
						       Bnd_Spline_East,
						       Bnd_Spline_West,
						       n_cells_i,
						       n_cells_j,
						       Number_of_Ghost_Cells,
						       Highest_Order_of_Reconstruction,
						       GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
						       Stretch_I,
						       Beta_I,
						       Tau_I,
						       Stretch_J,
						       Beta_J,
						       Tau_J,
						       Orthogonal_North,
						       Orthogonal_South,
						       Orthogonal_East,
						       Orthogonal_West);
    
  }

}



/*!
 * Generates a quadilateral mesh with clustering        
 * consisting of two grid blocks for predicting viscous 
 * flow and boundary layer development over a flat      
 * plate.                                               
 *                                                      
 * Usage: Grid_ptr = Grid_Flat_Plate9(Grid_ptr,         
 *                                   nblk_i,            
 *                                   nblk_j,            
 *                                   TWO,               
 *         		             100,               
 *         		             100,               
 *         		             2);                
 * This subroutine DOESN'T update the ghost cells or
 * the geometric properties of the grid cells.
 */
void Grid2D_Quad_MultiBlock_HO::Grid_Flat_Plate9_Without_Update(int &_Number_of_Blocks_Idir_,
								int &_Number_of_Blocks_Jdir_,
								const double &Length,
								const int &Flat_Plate_BC_Type,
								const int &Stretching_Flag,
								const double &Stretching_Factor_Idir,
								const double &Stretching_Factor_Jdir,
								const int Number_of_Cells_Idir,
								const int Number_of_Cells_Jdir,
								const int Number_of_Ghost_Cells,
								const int Highest_Order_of_Reconstruction) {
  
  int n_cells_i, n_cells_j, 
    Stretch_I, Stretch_J,
    Orthogonal_North, Orthogonal_South,
    Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D xc_NW, xc_NE, xc_SE, xc_SW;
  Spline2D_HO Bnd_Spline_North, Bnd_Spline_South,
              Bnd_Spline_East, Bnd_Spline_West;

  // Allocate memory for grid blocks.  There are two grid blocks for 
  // this mesh.
  _Number_of_Blocks_Idir_ = 2;
  _Number_of_Blocks_Jdir_ = 1;
  allocate(_Number_of_Blocks_Idir_,_Number_of_Blocks_Jdir_);

  // Create the mesh for each block representing the complete grid.
  for (int iBlk = 0; iBlk < Number_of_Blocks_Idir; iBlk++) {

    // Assign values to the locations of the corners of the 
    // rectangular box shaped domain representing each of the 
    // blocks in the grid.
    if (iBlk == 0) {
      xc_NW = TWO*Vector2D(-Length,Length);
      xc_NE = TWO*Vector2D(ZERO,Length);
      xc_SE = TWO*Vector2D(ZERO,ZERO);
      xc_SW = TWO*Vector2D(-Length,ZERO);
    } else {
      xc_NW = TWO*Vector2D(ZERO,Length);
      xc_NE = TWO*Vector2D(Length,Length);
      xc_SE = TWO*Vector2D(Length,ZERO);
      xc_SW = TWO*Vector2D(ZERO,  ZERO);
    }

    // Create the splines defining the north, south, east, and west 
    // boundaries of the rectangular boxes.
    Bnd_Spline_North.Create_Spline_Line(xc_NW,xc_NE,2);
    Bnd_Spline_South.Create_Spline_Line(xc_SW,xc_SE,2);
    Bnd_Spline_East.Create_Spline_Line(xc_SE,xc_NE,2);
    Bnd_Spline_West.Create_Spline_Line(xc_SW,xc_NW,2);

    // Set the boundary condition types for each of the boundary splines.
    if (iBlk == 0) {
      Bnd_Spline_North.setBCtype(BC_FIXED);
      Bnd_Spline_South.setBCtype(BC_REFLECTION);
      Bnd_Spline_East.setBCtype(BC_NONE);
      Bnd_Spline_West.setBCtype(BC_FIXED);
    } else {
      Bnd_Spline_North.setBCtype(BC_FIXED);
      Bnd_Spline_South.setBCtype(Flat_Plate_BC_Type);
      Bnd_Spline_East.setBCtype(BC_CONSTANT_EXTRAPOLATION);
      Bnd_Spline_West.setBCtype(BC_NONE);
    }

    // Determine the number of cells for this block.
    n_cells_i = Number_of_Cells_Idir;///2;
    n_cells_j = Number_of_Cells_Jdir;

    // Assign values to the stretching function parameters and boundary
    // grid line orthogonality parameters.
    Stretch_I = STRETCHING_FCN_LINEAR;
    Beta_I = ZERO;
    Tau_I = ZERO;
    Stretch_J = STRETCHING_FCN_LINEAR;
    Beta_J = ZERO;
    Tau_J = ZERO;
//     if (Stretching_Flag) {
//       if (iBlk == 0) {
// 	Stretch_I = STRETCHING_FCN_MAX_CLUSTERING;
// 	Beta_I = ONE + (Stretching_Factor_Idir-ONE)/6.50;
//       } else {
// 	Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
// 	Beta_I = Stretching_Factor_Idir;
//       }
//       Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
//       Beta_J = Stretching_Factor_Jdir;
//     }
    if (Stretching_Flag) {
      if (iBlk == 0) Stretch_I = STRETCHING_FCN_MAX_CLUSTERING;
      else Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
      Beta_I = Stretching_Factor_Idir;
      Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
      Beta_J = Stretching_Factor_Jdir;
    }
    Orthogonal_North = 0;
    Orthogonal_South = 0;
    Orthogonal_East = 0;
    Orthogonal_West = 0;

    // Create the 2D quadrilateral grid block.
    Grid_ptr[iBlk][0].Create_Quad_Block_Without_Update(Bnd_Spline_North,
						       Bnd_Spline_South,
						       Bnd_Spline_East,
						       Bnd_Spline_West,
						       n_cells_i,
						       n_cells_j,
						       Number_of_Ghost_Cells,
						       Highest_Order_of_Reconstruction,
						       GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
						       Stretch_I,
						       Beta_I,
						       Tau_I,
						       Stretch_J,
						       Beta_J,
						       Tau_J,
						       Orthogonal_North,
						       Orthogonal_South,
						       Orthogonal_East,
						       Orthogonal_West);
    
  }


}

/*!
 * Generates a quadilateral mesh with clustering        
 * consisting of one block clustered in the middle      
 * for predicting 1D flame speeds of a particular       
 * Chemical Mechanism.                                  
 *                                                      
 * Usage: Grid_ptr = Grid_1D_Flame (Grid_ptr,           
 *                                   nblk_i,            
 *                                   nblk_j,            
 *                                   TWO,               
 *                                   0.2,               
 *         		             100,               
 *         		             10,                
 *                                   2);                
 * This subroutine DOESN'T update the ghost cells or
 * the geometric properties of the grid cells.                                                      
 */
void Grid2D_Quad_MultiBlock_HO::Grid_1D_Flame_Without_Update(int &_Number_of_Blocks_Idir_,
							     int &_Number_of_Blocks_Jdir_,
							     const double &Length,
							     const double &Heigth,
							     const int Number_of_Cells_Idir,
							     const int Number_of_Cells_Jdir,
							     const int Number_of_Ghost_Cells,
							     const int Highest_Order_of_Reconstruction) {
  
  int iBlk, n_cells_i, n_cells_j, 
    Stretch_I, Stretch_J,
    Orthogonal_North, Orthogonal_South,
    Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D xc_NW, xc_NE, xc_SE, xc_SW;
  Spline2D_HO Bnd_Spline_North, Bnd_Spline_South,
              Bnd_Spline_East, Bnd_Spline_West;
  
  /* Allocate memory for grid blocks.  There are two grid
     blocks for this mesh. */
  
  _Number_of_Blocks_Idir_ = 1;
  _Number_of_Blocks_Jdir_ = 1;
  allocate(_Number_of_Blocks_Idir_,_Number_of_Blocks_Jdir_);
  
    /* Create the mesh for each block representing
       the complete grid. */

  for ( iBlk = 0; iBlk <= Number_of_Blocks_Idir-1; ++iBlk ) {
      
    /* Assign values to the locations of the corners
       of the rectangular box shaped domain representing
       each of the blocks in the grid. */
    
    xc_NW = Vector2D(ZERO  , Heigth);
    xc_NE = Vector2D(Length/Number_of_Blocks_Idir, Heigth);
    xc_SE = Vector2D(Length/Number_of_Blocks_Idir, ZERO);
    xc_SW = Vector2D(ZERO  , ZERO);
    
    /* Create the splines defining the north, south,
       east, and west boundaries of the rectangular boxes. */
    
    Bnd_Spline_North.Create_Spline_Line(xc_NW, xc_NE, 2);
    Bnd_Spline_South.Create_Spline_Line(xc_SW, xc_SE, 2);
    Bnd_Spline_East.Create_Spline_Line(xc_SE, xc_NE, 2);
    Bnd_Spline_West.Create_Spline_Line(xc_SW, xc_NW, 2);
    
    /* Set the boundary condition types for each of the
       boundary splines. */

    //Bnd_Spline_West.setBCtype(BC_FIXED); 
    Bnd_Spline_West.setBCtype(BC_1DFLAME_INFLOW); 
    Bnd_Spline_North.setBCtype(BC_CONSTANT_EXTRAPOLATION); // BC_REFLECTION  
    Bnd_Spline_South.setBCtype(BC_CONSTANT_EXTRAPOLATION); // BC_REFLECTION
    Bnd_Spline_East.setBCtype(BC_1DFLAME_OUTFLOW); 

    /* Determine the number of cells for this block. */
    
    n_cells_i = Number_of_Cells_Idir; 
    n_cells_j = Number_of_Cells_Jdir;
    
    /* Assign values to the stretching function parameters
           and boundary grid line orthogonality parameters. */
    
    Stretch_J = STRETCHING_FCN_LINEAR;
    Beta_J = ZERO;
    Tau_J = ZERO;
    Orthogonal_North = 0;
    Orthogonal_South = 0;
    Orthogonal_East = 0;
    Orthogonal_West = 0;
    Stretch_I = STRETCHING_FCN_MIDPT_CLUSTERING;
    Beta_I=0.5;
    Tau_I= 8.0;
    //Stretch_I = STRETCHING_FCN_LINEAR;
    
    /* Create the 2D quadrilateral grid block. */
    
    Grid_ptr[iBlk][0].Create_Quad_Block_Without_Update(Bnd_Spline_North,
						       Bnd_Spline_South,
						       Bnd_Spline_East,
						       Bnd_Spline_West,
						       n_cells_i,
						       n_cells_j,
						       Number_of_Ghost_Cells,
						       Highest_Order_of_Reconstruction,
						       GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
						       Stretch_I,
						       Beta_I, 
						       Tau_I,
						       Stretch_J,
						       Beta_J,
						       Tau_J,
						       Orthogonal_North,
						       Orthogonal_South,
						       Orthogonal_East,
						       Orthogonal_West);
    
  }

}


/*!
 * Generates a quadilateral mesh with clustering        
 * along the centerline (West) and entry (south)        
 * for the predicition of 2D laminar diffusion flames                         
 *                                                      
 * Usage: Grid_ptr = Grid_2D_Laminar_Flame (Grid_ptr,   
 *                                   nblk_i,            
 *                                   nblk_j,            
 *                                   TWO,               
 *                                   0.2,               
 *         		             100,               
 *         		             10,                
 *                                   2);                
 *
 * This subroutine DOESN'T update the ghost cells or
 * the geometric properties of the grid cells.                                                      
 */
void Grid2D_Quad_MultiBlock_HO::Grid_2D_Laminar_Flame_Without_Update(int &_Number_of_Blocks_Idir_,
								     int &_Number_of_Blocks_Jdir_,
								     const double &Length,
								     const double &Heigth,
								     const int Number_of_Cells_Idir,
								     const int Number_of_Cells_Jdir, 
								     const int Number_of_Ghost_Cells,
								     const int Highest_Order_of_Reconstruction,
								     const int Flame_Type_Flag) {
  

  int  n_cells_i, n_cells_j, Stretch_I, Stretch_J,
    Orthogonal_North, Orthogonal_South,
    Orthogonal_East, Orthogonal_West,
    Number_of_Blocks_Fuel,  Number_of_Blocks_Gap,
    Number_of_Blocks_Air, Number_of_Blocks_Free;

  double Beta_I, Tau_I, Beta_J, Tau_J, Top, Bot, East,West,
    fuel_spacing,  tube_spacing, air_spacing;
  
  Vector2D xc_NW, xc_NE, xc_SE, xc_SW;
  Spline2D_HO Bnd_Spline_North, Bnd_Spline_South,
              Bnd_Spline_East, Bnd_Spline_West;

  /******************************************************/
  //Standard Core Flame  //Flame_Type_Flag == IC_RESTART
  if( (Flame_Type_Flag == IC_RESTART || Flame_Type_Flag == IC_CHEM_CORE_FLAME) && _Number_of_Blocks_Idir_ == 12 ){
    fuel_spacing = 0.002;                   //m 
    tube_spacing = fuel_spacing + 0.00038;  //m 
    air_spacing = 0.025 - tube_spacing;     //m 
    
    //I-direction (inlet) blocks 
    Number_of_Blocks_Fuel = 4;
    Number_of_Blocks_Gap  = 1; 
    Number_of_Blocks_Air = 7;
    Number_of_Blocks_Free = 0;
    
  } else if( (Flame_Type_Flag == IC_RESTART || Flame_Type_Flag == IC_CHEM_CORE_FLAME) && _Number_of_Blocks_Idir_ == 3 ){
    
    fuel_spacing = 0.002;                   //m 
    tube_spacing = fuel_spacing;            //m 
    air_spacing = 0.025 - tube_spacing;     //m 
    
    Number_of_Blocks_Fuel = 1;
    Number_of_Blocks_Gap  = 0;
    Number_of_Blocks_Air = 2;
    Number_of_Blocks_Free = 0;  
       
  } else if( (Flame_Type_Flag == IC_RESTART || Flame_Type_Flag == IC_CHEM_CORE_FLAME) && _Number_of_Blocks_Idir_ == 1 ){
    
    fuel_spacing = 0.002;                   //m 
    tube_spacing = fuel_spacing;            //m 
    air_spacing = 0.025 - tube_spacing;     //m 
    
    Number_of_Blocks_Fuel = 1;
    Number_of_Blocks_Gap  = 0;
    Number_of_Blocks_Air = 0;
    Number_of_Blocks_Free = 0;
    
  } else if( Flame_Type_Flag == IC_CHEM_INVERSE_FLAME && _Number_of_Blocks_Idir_ == 8){
    
    /******************************************************/
    //Inverse Flames 1
    fuel_spacing = 0.0055;                 //m 
    tube_spacing = fuel_spacing + 0.0007;  //m 
    air_spacing = 0.020 - tube_spacing;    //m 
    
    //I-direction (inlet) blocks 
    Number_of_Blocks_Fuel = 3; //air
    Number_of_Blocks_Gap  = 1; 
    Number_of_Blocks_Air = 4;  //fuel
    Number_of_Blocks_Free = 0;
    
    //       /******************************************************/
    //       //Inverse Flames 2
    //        fuel_spacing = 0.005;                 //m 
    //        tube_spacing = fuel_spacing;  //m 
    //        air_spacing = 0.05 - tube_spacing;    //m 
    
    //       //I-direction (inlet) blocks 
    //        Number_of_Blocks_Fuel = 2;
    //        Number_of_Blocks_Gap  = 0; 
    //        Number_of_Blocks_Air = 4;
    //        Number_of_Blocks_Free = 2;
    
  } else {
    cerr<<"\n Initial conditions not valid for 2D Laminar Flame Grid"; cout.flush();
    exit(1);
  }

  if( _Number_of_Blocks_Idir_ != Number_of_Blocks_Fuel + Number_of_Blocks_Gap + Number_of_Blocks_Air + Number_of_Blocks_Free ){
    cout<<"\n WARNING: Grid_2D_Laminar_Flame has a fixed initial number of Blocks in the x-direction to insure proper BC's, ";
    cout<<" currently it is set to "<< Number_of_Blocks_Fuel + Number_of_Blocks_Gap + Number_of_Blocks_Air + Number_of_Blocks_Free; 
     }

  _Number_of_Blocks_Idir_ = Number_of_Blocks_Fuel + Number_of_Blocks_Gap + Number_of_Blocks_Air + Number_of_Blocks_Free;

  /* Allocate memory for grid blocks. */
  allocate(_Number_of_Blocks_Idir_, _Number_of_Blocks_Jdir_);
  
  /* Create the mesh for each block representing the complete grid. */
  
  for ( int iBlk = 0; iBlk < Number_of_Blocks_Idir; iBlk++ ) {
    for ( int jBlk = 0; jBlk < Number_of_Blocks_Jdir; jBlk++ ) {
      
      /* Assign values to the locations of the corners  of the rectangular box shaped domain representing
	 each of the blocks in the grid. */
      
      Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
      Beta_J = 1.05; 
      Tau_J = ZERO;
      Orthogonal_East = 0;
      Orthogonal_West = 0;
      Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
      Beta_I = 1.05; 
      Tau_I = ZERO;
      Orthogonal_North = 0;
      Orthogonal_South = 0;

      //Stretching for J Blocks
      Top = StretchingFcn(double(jBlk + 1)/double(Number_of_Blocks_Jdir), Beta_J, Tau_J, Stretch_J);
      Bot = StretchingFcn(double(jBlk)/double(Number_of_Blocks_Jdir), Beta_J, Tau_J, Stretch_J);  

      /**************** INLET **********************/
      // Only one block 
      if( Number_of_Blocks_Idir == 1 ) {
	xc_NW = Vector2D( ZERO , Heigth);
	xc_NE = Vector2D( Length , Heigth); 
	xc_SE = Vector2D( Length , ZERO); 
	xc_SW = Vector2D( ZERO , ZERO); 
	 
	//Fuel Inlet
      } else if(iBlk < Number_of_Blocks_Fuel ) {
	xc_NW = Vector2D( double(iBlk) * fuel_spacing/double(Number_of_Blocks_Fuel), Top*Heigth);
	xc_NE = Vector2D( double(iBlk + 1) * fuel_spacing/double(Number_of_Blocks_Fuel), Top*Heigth);			    
	xc_SE = Vector2D( double(iBlk + 1) * fuel_spacing/double(Number_of_Blocks_Fuel), Bot*Heigth);			    
	xc_SW = Vector2D( double(iBlk) * fuel_spacing/double(Number_of_Blocks_Fuel), Bot*Heigth);
			    
	//Dead space
      } else if( iBlk == Number_of_Blocks_Gap + Number_of_Blocks_Fuel - 1 ) {     
	xc_NW = Vector2D( fuel_spacing, Top*Heigth);
	xc_NE = Vector2D( tube_spacing,Top*Heigth); 
	xc_SE = Vector2D( tube_spacing, Bot*Heigth); 
	xc_SW = Vector2D( fuel_spacing, Bot*Heigth); 

	//Air Inlet
      } else if( iBlk < Number_of_Blocks_Gap + Number_of_Blocks_Fuel + Number_of_Blocks_Air) {

	//Air Coflow Block Stretching I-direction
	West = StretchingFcn(double(iBlk - (Number_of_Blocks_Fuel + Number_of_Blocks_Gap))/double(Number_of_Blocks_Air), Beta_I, Tau_I, Stretch_I);
	East = StretchingFcn(double(iBlk- (Number_of_Blocks_Fuel + Number_of_Blocks_Gap) + 1)/double(Number_of_Blocks_Air), Beta_I, Tau_I, Stretch_I);  

	xc_NW = Vector2D( tube_spacing + West*air_spacing, Top*Heigth);			    
	xc_NE = Vector2D( tube_spacing + East*air_spacing, Top*Heigth);					 			   
	xc_SE = Vector2D( tube_spacing + East*air_spacing, Bot*Heigth);					   
	xc_SW = Vector2D( tube_spacing + West*air_spacing, Bot*Heigth);
			   	  
	//Quiesent Air
      } else {
	xc_NW = Vector2D( tube_spacing + air_spacing + double(iBlk - (Number_of_Blocks_Fuel + Number_of_Blocks_Gap + Number_of_Blocks_Air)) 
			  * (Length - air_spacing - tube_spacing)/double(Number_of_Blocks_Free), Top*Heigth);
	xc_NE = Vector2D( tube_spacing + air_spacing + double(iBlk - (Number_of_Blocks_Fuel + Number_of_Blocks_Gap + Number_of_Blocks_Air) + 1) 
			  * (Length - air_spacing - tube_spacing)/double(Number_of_Blocks_Free), Top*Heigth);
	xc_SE = Vector2D( tube_spacing + air_spacing + double(iBlk - (Number_of_Blocks_Fuel + Number_of_Blocks_Gap + Number_of_Blocks_Air) + 1) 
			  * (Length - air_spacing - tube_spacing)/double(Number_of_Blocks_Free), Bot*Heigth);	
	xc_SW = Vector2D( tube_spacing + air_spacing + double(iBlk - (Number_of_Blocks_Fuel + Number_of_Blocks_Gap + Number_of_Blocks_Air ))
			  * (Length - air_spacing - tube_spacing)/double(Number_of_Blocks_Free), Bot*Heigth);	
      }
	
      Bnd_Spline_North.Create_Spline_Line(xc_NW, xc_NE, 2);
      Bnd_Spline_South.Create_Spline_Line(xc_SW, xc_SE, 2);
      Bnd_Spline_East.Create_Spline_Line(xc_SE, xc_NE, 2);
      Bnd_Spline_West.Create_Spline_Line(xc_SW, xc_NW, 2);
	
      /* Set the boundary condition types for each of the  boundary splines. */
      if( Number_of_Blocks_Idir == 1 && Number_of_Blocks_Jdir ==1 ) {
	Bnd_Spline_West.setBCtype(BC_REFLECTION);   
	Bnd_Spline_East.setBCtype(BC_REFLECTION);          //BC_FREE_SLIP_ISOTHERMAL);  
	Bnd_Spline_South.setBCtype(BC_FIXED);
	Bnd_Spline_North.setBCtype(BC_2DFLAME_OUTFLOW);
      } else {
	if (iBlk == 0 ) {
	  Bnd_Spline_West.setBCtype(BC_REFLECTION);    //centerline
	  Bnd_Spline_East.setBCtype(BC_NONE);
	} else if (iBlk == Number_of_Blocks_Idir-1 )  {
	  Bnd_Spline_West.setBCtype(BC_NONE);      
	  Bnd_Spline_East.setBCtype(BC_REFLECTION); //BC_FREE_SLIP_ISOTHERMAL);  //FIXED   //farfield right  
	} else { 
	  Bnd_Spline_West.setBCtype(BC_NONE);      
	  Bnd_Spline_East.setBCtype(BC_NONE);
	}
	  
	if (jBlk == 0 && iBlk == Number_of_Blocks_Gap + Number_of_Blocks_Fuel - 1 && Number_of_Blocks_Gap!=0) {
	  Bnd_Spline_South.setBCtype(BC_WALL_VISCOUS_HEATFLUX);      //Gap wall
	  Bnd_Spline_North.setBCtype(BC_NONE);
	} else if (jBlk == 0 && iBlk < Number_of_Blocks_Gap + Number_of_Blocks_Fuel + Number_of_Blocks_Air ) {
	  Bnd_Spline_South.setBCtype(BC_2DFLAME_INFLOW);     //BC_FIXED);        //Bottom Inflow Left
	  Bnd_Spline_North.setBCtype(BC_NONE);
	} else if(jBlk == 0) {
	  Bnd_Spline_South.setBCtype(BC_2DFLAME_INFLOW);    //BC_FIXED); //BC_WALL_VISCOUS_HEATFLUX);     //Bottom Right 
	  Bnd_Spline_North.setBCtype(BC_NONE);
	} else if(jBlk == Number_of_Blocks_Jdir-1 ){  
	  Bnd_Spline_North.setBCtype(BC_2DFLAME_OUTFLOW);   //BC_CHARACTERISTIC); //Top outflow 
	  Bnd_Spline_South.setBCtype(BC_NONE);
	} else {
	  Bnd_Spline_South.setBCtype(BC_NONE);
	  Bnd_Spline_North.setBCtype(BC_NONE);
	}
	
	
	/* Assign values to the stretching function parameters
	   and boundary grid line orthogonality parameters. */
	Stretch_I = STRETCHING_FCN_LINEAR; 
	Beta_I = ZERO; 
	Tau_I = ZERO;
	Stretch_J = STRETCHING_FCN_LINEAR;
	Beta_J = ZERO;
	Tau_J = ZERO; 
      }
	
      /* Determine the number of cells for this block. */
      n_cells_i = Number_of_Cells_Idir/Number_of_Blocks_Idir;
      n_cells_j = Number_of_Cells_Jdir/Number_of_Blocks_Jdir;
	
      /* Create the 2D quadrilateral grid block. */
      Grid_ptr[iBlk][jBlk].Create_Quad_Block_Without_Update(Bnd_Spline_North,
							    Bnd_Spline_South,
							    Bnd_Spline_East,
							    Bnd_Spline_West,
							    n_cells_i,
							    n_cells_j,
							    Number_of_Ghost_Cells,
							    Highest_Order_of_Reconstruction,
							    GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
							    Stretch_I,
							    Beta_I, 
							    Tau_I,
							    Stretch_J,
							    Beta_J,
							    Tau_J,
							    Orthogonal_North,
							    Orthogonal_South,
							    Orthogonal_East,
							    Orthogonal_West);

    } 
  }

}


/*!
 * Generates a single block quadilateral mesh with      
 * clustering for predicting viscous flow and boundary  
 * layer development in a cylindrical duct or pipe.     
 *                                                      
 * Usage: Grid_ptr = Grid_Pipe(Grid_ptr,                
 *                             nblk_i,                  
 *                             nblk_j,                  
 *                             TEN,                     
 *                             HALF,                    
 *   	                       100,                     
 *  	                       50,                      
 *                             2);                      
 * This subroutine DOESN'T update the ghost cells or
 * the geometric properties of the grid cells.
 */
void Grid2D_Quad_MultiBlock_HO::Grid_Pipe_Without_Update(int &_Number_of_Blocks_Idir_,
							 int &_Number_of_Blocks_Jdir_,
							 const double &Length,
							 const double &Radius,
							 const int Stretching_Flag,
							 const double Stretching_Factor,
							 const int Number_of_Cells_Idir,
							 const int Number_of_Cells_Jdir,
							 const int Number_of_Ghost_Cells,
							 const int Highest_Order_of_Reconstruction) {
  
  int Stretch_I, Stretch_J,
    Orthogonal_North, Orthogonal_South,
    Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D xc_NW, xc_NE, xc_SE, xc_SW;
  Spline2D_HO Bnd_Spline_North, Bnd_Spline_South,
              Bnd_Spline_East, Bnd_Spline_West;
  
  /* Allocate memory for grid block. */
  
  _Number_of_Blocks_Idir_ = 1;
  _Number_of_Blocks_Jdir_ = 1;
  allocate(_Number_of_Blocks_Idir_,_Number_of_Blocks_Jdir_);
  
  /* Assign values to the locations of the corners
     of the rectangular box defining the pipe geometry. */
  xc_NW = Vector2D(ZERO , Radius);
  xc_NE = Vector2D(Length, Radius);
  xc_SE = Vector2D(Length, ZERO);
  xc_SW = Vector2D(ZERO , ZERO);
 
  /* Create the splines defining the north, south,
     east, and west boundaries of the grid. */
  
  Bnd_Spline_North.Create_Spline_Line(xc_NW, xc_NE, 2);
  Bnd_Spline_South.Create_Spline_Line(xc_SW, xc_SE, 2);
  Bnd_Spline_East.Create_Spline_Line(xc_SE, xc_NE, 2);
  Bnd_Spline_West.Create_Spline_Line(xc_SW, xc_NW, 2);
  
  /* Set the boundary condition types for each of the
     boundary splines. */
   
  Bnd_Spline_North.setBCtype(BC_WALL_VISCOUS_ISOTHERMAL);
  Bnd_Spline_South.setBCtype(BC_REFLECTION);
  Bnd_Spline_East.setBCtype(BC_FIXED);
  Bnd_Spline_West.setBCtype(BC_FIXED);

  /* Assign values to the stretching function parameters
     and boundary grid line orthogonality parameters. */
  Stretch_I = STRETCHING_FCN_LINEAR;
  Beta_I = ZERO; 
  Tau_I = ZERO;
  Stretch_J = STRETCHING_FCN_MAX_CLUSTERING;
  Beta_J = ZERO;
  Tau_J = ZERO; 
  if (Stretching_Flag) Beta_J = Stretching_Factor;
  else Stretch_J = STRETCHING_FCN_LINEAR;

  Orthogonal_North = 0;
  Orthogonal_South = 0;
  Orthogonal_East = 0;
  Orthogonal_West = 0;

  /* Create the 2D quadrilateral grid block representing
       the mesh. */
  Grid_ptr[0][0].Create_Quad_Block_Without_Update(Bnd_Spline_North,
						  Bnd_Spline_South,
						  Bnd_Spline_East,
						  Bnd_Spline_West,
						  Number_of_Cells_Idir,
						  Number_of_Cells_Jdir,
						  Number_of_Ghost_Cells,
						  Highest_Order_of_Reconstruction,
						  GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
						  Stretch_I,
						  Beta_I, 
						  Tau_I,
						  Stretch_J,
						  Beta_J,
						  Tau_J,
						  Orthogonal_North,
						  Orthogonal_South,
						  Orthogonal_East,
						  Orthogonal_West);
  
}

/*!
 * Generates a single block quadilateral mesh with      
 * clustering for predicting viscous flow and boundary  
 * layer development in a cylindrical duct or pipe.     
 *                                                      
 * Usage: Grid_ptr = Grid_Pipe(Grid_ptr,                
 *                             nblk_i,                  
 *                             nblk_j,                  
 *                             TEN,                     
 *                             HALF,                    
 *   	                       100,                     
 *  	                       50,                      
 *                             2);                      
 * This subroutine DOESN'T update the ghost cells or
 * the geometric properties of the grid cells.
 */
void Grid2D_Quad_MultiBlock_HO::Grid_Pipe_Without_Update(int &_Number_of_Blocks_Idir_,
							 int &_Number_of_Blocks_Jdir_,
							 const double &Length,
							 const double &Radius,
							 const int &Axisymmetric,
							 const int Number_of_Cells_Idir,
							 const int Number_of_Cells_Jdir,
							 const int Number_of_Ghost_Cells,
							 const int Highest_Order_of_Reconstruction) {
  
  int Stretch_I, Stretch_J,
      Orthogonal_North, Orthogonal_South,
      Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D xc_NW, xc_NE, xc_SE, xc_SW;
  Spline2D_HO Bnd_Spline_North, Bnd_Spline_South,
              Bnd_Spline_East, Bnd_Spline_West;

  /* Allocate memory for grid block. */
  
  _Number_of_Blocks_Idir_ = 1;
  _Number_of_Blocks_Jdir_ = 1;
  allocate(_Number_of_Blocks_Idir_, _Number_of_Blocks_Jdir_);

  /* Assign values to the locations of the corners
     of the rectangular box defining the pipe geometry. */

  if(Axisymmetric ==1){
    xc_NW = Vector2D(ZERO , Radius);
    xc_NE = Vector2D(Length, Radius);
    xc_SE = Vector2D(Length, ZERO);
    xc_SW = Vector2D(ZERO , ZERO);
  }
  if(Axisymmetric ==2){
    xc_NW = Vector2D(ZERO , Length);
    xc_NE = Vector2D(Radius,Length);
    xc_SE = Vector2D(Radius, ZERO);
    xc_SW = Vector2D(ZERO , ZERO); //axisymmetric 2
  }
  /* Create the splines defining the north, south,
     east, and west boundaries of the grid. */

  Bnd_Spline_North.Create_Spline_Line(xc_NW, xc_NE, 2);
  Bnd_Spline_South.Create_Spline_Line(xc_SW, xc_SE, 2);
  Bnd_Spline_East.Create_Spline_Line(xc_SE, xc_NE, 2);
  Bnd_Spline_West.Create_Spline_Line(xc_SW, xc_NW, 2);

  /* Set the boundary condition types for each of the
     boundary splines. */

  if (Axisymmetric == 1) {
    Bnd_Spline_North.setBCtype(BC_WALL_VISCOUS_ISOTHERMAL);
    Bnd_Spline_South.setBCtype(BC_REFLECTION);
    Bnd_Spline_East.setBCtype(BC_FIXED);
    Bnd_Spline_West.setBCtype(BC_FIXED);

    /* Assign values to the stretching function parameters
       and boundary grid line orthogonality parameters. */
      
    Stretch_I = STRETCHING_FCN_LINEAR;
    Beta_I = ZERO; 
    Tau_I = ZERO;
    Stretch_J = STRETCHING_FCN_MAX_CLUSTERING;
    Beta_J = 1.0001;
    Tau_J = ZERO;  //axisymmetric #1
  }

  if (Axisymmetric == 2) {
    Bnd_Spline_East.setBCtype(BC_WALL_VISCOUS_ISOTHERMAL);
    Bnd_Spline_West.setBCtype(BC_REFLECTION);
    Bnd_Spline_North.setBCtype(BC_FIXED);
    Bnd_Spline_South.setBCtype(BC_FIXED);
    /* Assign values to the stretching function parameters
       and boundary grid line orthogonality parameters. */
    Stretch_I =  STRETCHING_FCN_MAX_CLUSTERING;
    Beta_I = 1.0001;
    Tau_I = ZERO;
    Stretch_J = STRETCHING_FCN_LINEAR;
    Beta_J =  ZERO; 
    Tau_J = ZERO;  
  }

  Orthogonal_North = 0;
  Orthogonal_South = 0;
  Orthogonal_East = 0;
  Orthogonal_West = 0;

  /* Create the 2D quadrilateral grid block representing
     the mesh. */

  Grid_ptr[0][0].Create_Quad_Block_Without_Update(Bnd_Spline_North,
						  Bnd_Spline_South,
						  Bnd_Spline_East,
						  Bnd_Spline_West,
						  Number_of_Cells_Idir,
						  Number_of_Cells_Jdir,
						  Number_of_Ghost_Cells,
						  Highest_Order_of_Reconstruction,
						  GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
						  Stretch_I,
						  Beta_I, 
						  Tau_I,
						  Stretch_J,
						  Beta_J,
						  Tau_J,
						  Orthogonal_North,
						  Orthogonal_South,
						  Orthogonal_East,
						  Orthogonal_West);

}

/*!
 * Generates a single block quadilateral mesh with      
 * clustering for predicting supersonic flow around     
 * a cirucular cylinder blunt body.                     
 *                                                      
 * Usage: Grid_ptr = Grid_Blunt_Body(Grid_ptr,          
 *                                   nblk_i,            
 *                                   nblk_j,            
 *                                   ONE,               
 *                                   FOUR,              
 *   	                             150,               
 *  	                             50,                
 *                                   2);                
 *                                                      
 * This subroutine DOESN'T update the ghost cells or
 * the geometric properties of the grid cells.
 */
void  Grid2D_Quad_MultiBlock_HO::Grid_Blunt_Body_Without_Update(int &_Number_of_Blocks_Idir_,
								int &_Number_of_Blocks_Jdir_,
								const double &Radius,
								const double &Mach_Number,
								const int Number_of_Cells_Idir,
								const int Number_of_Cells_Jdir,
								const int Number_of_Ghost_Cells,
								const int Highest_Order_of_Reconstruction) {

  int Stretch_I, Stretch_J,
    Orthogonal_North, Orthogonal_South,
    Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D x1, x2;
  Spline2D_HO Bnd_Spline_North, Bnd_Spline_South,
              Bnd_Spline_East, Bnd_Spline_West;

  /* Allocate memory for grid block. */

  _Number_of_Blocks_Idir_ = 1;
  _Number_of_Blocks_Jdir_ = 1;
  allocate(_Number_of_Blocks_Idir_, _Number_of_Blocks_Jdir_);

  /* Create the splines defining the north, south,
     east, and west boundaries of the grid. */

  Bnd_Spline_North.Create_Spline_Bow_Shock(Radius,
					   Mach_Number,
					   1,
					   181);
  x1 = Vector2D(ZERO , ZERO);
  Bnd_Spline_South.Create_Spline_Circular_Arc(x1,
					      Radius,
					      180.00,
					      90.00,
					      181);
  x1 = Bnd_Spline_South.Xp[0];
  x2 = Bnd_Spline_North.Xp[0];
  Bnd_Spline_West.Create_Spline_Line(x1, x2, 2);
  x1 = Bnd_Spline_South.Xp[Bnd_Spline_South.np-1];
  x2 = Bnd_Spline_North.Xp[Bnd_Spline_North.np-1];
  Bnd_Spline_East.Create_Spline_Line(x1, x2, 2);

  /* Set the boundary condition types for each of the
     boundary splines. */

  Bnd_Spline_North.setBCtype(BC_FIXED);
  Bnd_Spline_South.setBCtype(BC_REFLECTION);
  Bnd_Spline_East.setBCtype(BC_CONSTANT_EXTRAPOLATION);
  Bnd_Spline_West.setBCtype(BC_REFLECTION);

  /* Assign values to the stretching function parameters
     and boundary grid line orthogonality parameters. */

  Stretch_I = STRETCHING_FCN_LINEAR;
  Beta_I = ZERO; 
  Tau_I = ZERO;
  Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
  Beta_J = 1.05;
  Tau_J = ZERO;
  Orthogonal_North = 0;
  Orthogonal_South = 1;
  Orthogonal_East = 1;
  Orthogonal_West = 1;

  /* Create the 2D quadrilateral grid block representing
     the mesh. */

  Grid_ptr[0][0].Create_Quad_Block_Without_Update(Bnd_Spline_North,
						  Bnd_Spline_South,
						  Bnd_Spline_East,
						  Bnd_Spline_West,
						  Number_of_Cells_Idir,
						  Number_of_Cells_Jdir,
						  Number_of_Ghost_Cells,
						  Highest_Order_of_Reconstruction,
						  GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
						  Stretch_I,
						  Beta_I, 
						  Tau_I,
						  Stretch_J,
						  Beta_J,
						  Tau_J,
						  Orthogonal_North,
						  Orthogonal_South,
						  Orthogonal_East,
						  Orthogonal_West);
    
  /* Smooth the 2D quadrilateral grid block. */

  Grid_ptr[0][0].Smooth_Quad_Block(min(250, 2*max(Number_of_Cells_Idir,Number_of_Cells_Jdir)));

}


/*!
 * Generates a quadilateral mesh with clustering        
 * consisting of two grid blocks for predicting the     
 * axisymmetric core flow in a cylindrical grain solid  
 * propellant rocket motor.                             
 *                                                      
 * Usage: Grid_ptr = Grid_Rocket_Motor(Grid_ptr,        
 *                                     nblk_i,          
 *                                     nblk_j,          
 *                                     0.8350,          
 *                                     0.0200,          
 *                                     0.0500,          
 *                                     0.1500,          
 *                                     0.0300,          
 *                                     0.0100,          
 *         		               100,             
 *         		               100,             
 *                                     2);              
 * This subroutine DOESN'T update the ghost cells or
 * the geometric properties of the grid cells.
 */
void Grid2D_Quad_MultiBlock_HO::Grid_Rocket_Motor_Without_Update(int &_Number_of_Blocks_Idir_,
								 int &_Number_of_Blocks_Jdir_,
								 const double &Length_Chamber,
								 const double &Radius_Chamber,
								 const double &Length_Chamber_To_Throat,
								 const double &Length_Nozzle,
								 const double &Radius_Nozzle_Exit,
								 const double &Radius_Nozzle_Throat,
								 const double &Radius_Grain,
								 const int &Nozzle_Type,
								 const int &Chamber_BC_Type,
								 const int &Stretching_Flag,
								 const int Stretching_Type_Jdir,
								 const double &Stretching_Factor_Idir,
								 const double &Stretching_Factor_Jdir,
								 const int Number_of_Cells_Idir,
								 const int Number_of_Cells_Jdir,
								 const int Number_of_Ghost_Cells,
								 const int Highest_Order_of_Reconstruction) {
  
  int error_flag, block_flag;
  int Stretch_I, Stretch_J,
    Orthogonal_North, Orthogonal_South,
    Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D xc_NW, xc_NE, xc_SE, xc_SW;
  Spline2D_HO Bnd_Spline_North, Bnd_Spline_South,
              Bnd_Spline_East, Bnd_Spline_West;
  double num1 = ZERO, num2 = ZERO;
  int iBlk_Chamber, iBlk_Nozzle;

  // Allocate memory for grid blocks.
  if (Radius_Grain > Radius_Chamber) return;
  if (Radius_Grain > ZERO) {
    iBlk_Nozzle = 2;
    _Number_of_Blocks_Jdir_ = 2;
  } else if (Radius_Grain < ZERO) {
    iBlk_Nozzle = 1;
    _Number_of_Blocks_Jdir_ = 2;
  } else {
    iBlk_Nozzle = 1;
    _Number_of_Blocks_Jdir_ = 1;
  }
  iBlk_Chamber = _Number_of_Blocks_Idir_ - 1;
  _Number_of_Blocks_Idir_ = iBlk_Chamber + iBlk_Nozzle;
  allocate(_Number_of_Blocks_Idir_,_Number_of_Blocks_Jdir_);

  // Create the mesh for each block representing the complete grid.
  for (int jBlk = 0; jBlk < Number_of_Blocks_Jdir; jBlk++) {
    for (int iBlk = 0; iBlk < Number_of_Blocks_Idir; iBlk++) {

      // Default the current block to off (a hole).
      block_flag = OFF;

      // Default the stretching parameters.
      Stretch_I = STRETCHING_FCN_LINEAR;
      Beta_I = ONE;
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_LINEAR;
      Beta_J = ONE;
      Tau_J = ZERO;
      Orthogonal_North = 1;
      Orthogonal_South = 1;
      Orthogonal_East = 1;
      Orthogonal_West = 1;

      // Create the splines defining the north, south, east, and west 
      // boundaries of the grid block.
      if (jBlk == 0 && iBlk < iBlk_Chamber) {
	num1 = Length_Chamber*double(iBlk)/double(Number_of_Blocks_Idir-iBlk_Nozzle);
	num2 = Length_Chamber*double(iBlk+1)/double(Number_of_Blocks_Idir-iBlk_Nozzle);

	//   	num1 = Length_Chamber*StretchingFcn(double(iBlk)/double(iBlk_Chamber),
	//             Stretching_Factor_Idir,ZERO,STRETCHING_FCN_MIN_CLUSTERING);
	//   	num2 = Length_Chamber*StretchingFcn(double(iBlk+1)/double(iBlk_Chamber),
	//             Stretching_Factor_Idir,ZERO,STRETCHING_FCN_MIN_CLUSTERING);

	xc_SE = Vector2D(-Length_Chamber+num2,ZERO);
	xc_SW = Vector2D(-Length_Chamber+num1,ZERO);
	if (iBlk_Nozzle == 1) {
	  xc_NW = Vector2D(-Length_Chamber+num1,Radius_Chamber);
	  xc_NE = Vector2D(-Length_Chamber+num2,Radius_Chamber);
	} else {
	  xc_NW = Vector2D(-Length_Chamber+num1,Radius_Grain);
	  xc_NE = Vector2D(-Length_Chamber+num2,Radius_Grain);
	}
	Bnd_Spline_North.Create_Spline_Line(xc_NW,xc_NE,2);
	Bnd_Spline_South.Create_Spline_Line(xc_SW,xc_SE,2);
	Bnd_Spline_East.Create_Spline_Line(xc_SE,xc_NE,2);
	Bnd_Spline_West.Create_Spline_Line(xc_SW,xc_NW,2);
	if (iBlk_Nozzle == 1 && Radius_Grain < ZERO) Bnd_Spline_North.setBCtype(BC_NONE); 
	else Bnd_Spline_North.setBCtype(BC_BURNING_SURFACE);
	Bnd_Spline_South.setBCtype(BC_REFLECTION);
	Bnd_Spline_East.setBCtype(BC_NONE);
	if (iBlk == 0) Bnd_Spline_West.setBCtype(Chamber_BC_Type);
	else Bnd_Spline_West.setBCtype(BC_NONE);
	if (Stretching_Flag) {
	  if (iBlk == 0) {
	    Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
	    Beta_I = Stretching_Factor_Idir;
	    Tau_I = ZERO;
	  } else if (iBlk == iBlk_Chamber-1) {
	    Stretch_I = STRETCHING_FCN_MAX_CLUSTERING;
	    Stretch_I = STRETCHING_FCN_LINEAR;
	    Beta_I = Stretching_Factor_Idir;
	    Tau_I = ZERO;
	  }
	}
	Stretch_J = Stretching_Type_Jdir;
	Beta_J = Stretching_Factor_Jdir;
	Tau_J = ZERO;
	block_flag = ON;

      } else if (jBlk == 1 && iBlk < iBlk_Chamber && Radius_Grain < ZERO) {
	num1 = Length_Chamber*double(iBlk)/double(Number_of_Blocks_Idir-iBlk_Nozzle);
	num2 = Length_Chamber*double(iBlk+1)/double(Number_of_Blocks_Idir-iBlk_Nozzle);

	//  	num1 = Length_Chamber*StretchingFcn(double(iBlk)/double(iBlk_Chamber),
	//             Stretching_Factor_Idir,ZERO,STRETCHING_FCN_MIN_CLUSTERING);
	//  	num2 = Length_Chamber*StretchingFcn(double(iBlk+1)/double(iBlk_Chamber),
	//             Stretching_Factor_Idir,ZERO,STRETCHING_FCN_MIN_CLUSTERING);

	xc_NW = Vector2D(-Length_Chamber+num1,Radius_Chamber-Radius_Grain);
	xc_NE = Vector2D(-Length_Chamber+num2,Radius_Chamber-Radius_Grain);
	xc_SE = Vector2D(-Length_Chamber+num2,Radius_Chamber);
	xc_SW = Vector2D(-Length_Chamber+num1,Radius_Chamber);
	Bnd_Spline_North.Create_Spline_Line(xc_NW,xc_NE,2);
	Bnd_Spline_South.Create_Spline_Line(xc_SW,xc_SE,2);
	Bnd_Spline_East.Create_Spline_Line(xc_SE,xc_NE,2);
	Bnd_Spline_West.Create_Spline_Line(xc_SW,xc_NW,2);
	Bnd_Spline_North.setBCtype(BC_BURNING_SURFACE);
	Bnd_Spline_South.setBCtype(BC_NONE);
	Bnd_Spline_East.setBCtype(BC_NONE);
	if (iBlk == iBlk_Chamber-1) Bnd_Spline_East.setBCtype(Chamber_BC_Type);
	else Bnd_Spline_East.setBCtype(BC_NONE);
	if (iBlk == 0) Bnd_Spline_West.setBCtype(Chamber_BC_Type);
	else Bnd_Spline_West.setBCtype(BC_NONE);
	if (Stretching_Flag) {
	  if (iBlk == 0) {
	    Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
	    Beta_I = Stretching_Factor_Idir;
	    Tau_I = ZERO;
	  } else if (iBlk == iBlk_Chamber-1) {
	    Stretch_I = STRETCHING_FCN_MAX_CLUSTERING;
	    Stretch_I = STRETCHING_FCN_LINEAR;
	    Beta_I = Stretching_Factor_Idir;
	    Tau_I = ZERO;
	  }
	}
	//Stretch_J = Stretching_Type_Jdir;
	//Beta_J = Stretching_Factor_Jdir;
	//Tau_J = ZERO;
	block_flag = ON;

      } else if (jBlk == 0 && iBlk >= iBlk_Chamber && iBlk_Nozzle == 1) {
	Bnd_Spline_North.Create_Spline_Area_Variation(ZERO,
						      Length_Chamber_To_Throat,
						      Length_Nozzle,
						      Radius_Chamber,
						      Radius_Nozzle_Throat,
						      Radius_Nozzle_Exit,
						      Nozzle_Type,
						      501);
	xc_NW = Vector2D(ZERO,Radius_Chamber);
	xc_NE = Bnd_Spline_North.Xp[500];
	xc_SE = Vector2D(Length_Nozzle,ZERO);
	xc_SW = Vector2D(ZERO,ZERO);
	Bnd_Spline_East.Create_Spline_Line(xc_SE,xc_NE,2);
	Bnd_Spline_West.Create_Spline_Line(xc_SW,xc_NW,2);
	Bnd_Spline_South.Create_Spline_Line(xc_SW,xc_SE,2);
	Bnd_Spline_North.setBCtype(Chamber_BC_Type);
	Bnd_Spline_South.setBCtype(BC_REFLECTION);
	Bnd_Spline_East.setBCtype(BC_CONSTANT_EXTRAPOLATION);
	Bnd_Spline_East.setBCtype(BC_CHARACTERISTIC);
	Bnd_Spline_West.setBCtype(BC_NONE);
	Stretch_I = STRETCHING_FCN_MIDPT_CLUSTERING;
	Beta_I = Length_Chamber_To_Throat/Length_Nozzle;
	Tau_I = 3.00;
	Stretch_J = Stretching_Type_Jdir;
	Beta_J = Stretching_Factor_Jdir;
	Tau_J = ZERO;
	block_flag = ON;

      } else if (jBlk == 0 && iBlk - iBlk_Chamber == 0 && iBlk_Nozzle == 2) {
	xc_SW = Vector2D(ZERO,ZERO);
	xc_SE = Vector2D(Length_Chamber_To_Throat,ZERO);
	xc_NE = Vector2D(Length_Chamber_To_Throat,HALF*Radius_Nozzle_Throat);
	xc_NW = Vector2D(ZERO,Radius_Grain);
	Bnd_Spline_East.Create_Spline_Line(xc_SE,xc_NE,2);
	Bnd_Spline_West.Create_Spline_Line(xc_SW,xc_NW,2);
	Bnd_Spline_South.Create_Spline_Line(xc_SW,xc_SE,2);
	Bnd_Spline_North.Create_Spline_Converging_Nozzle(ZERO,
							 Length_Chamber_To_Throat,
							 Radius_Grain,
							 HALF*Radius_Nozzle_Throat,
							 251);
	Bnd_Spline_North.setBCtype(BC_NONE);
	Bnd_Spline_South.setBCtype(BC_REFLECTION);
	Bnd_Spline_East.setBCtype(BC_NONE);
	Bnd_Spline_West.setBCtype(BC_NONE);
	if (Stretching_Flag) {
	  Stretch_I = STRETCHING_FCN_MINMAX_CLUSTERING;
	  Beta_I = Stretching_Factor_Idir;
	  Tau_I = ZERO;
	  Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
	  Beta_J = Stretching_Factor_Jdir;
	  Tau_J = ZERO;
	}
	block_flag = ON;

      } else if (jBlk == 0 && iBlk - iBlk_Chamber == 1 && iBlk_Nozzle == 2) {
	xc_SW = Vector2D(Length_Chamber_To_Throat,ZERO);
	xc_SE = Vector2D(Length_Nozzle,ZERO);
	xc_NE = Vector2D(Length_Nozzle,HALF*Radius_Nozzle_Exit);
	xc_NW = Vector2D(Length_Chamber_To_Throat,HALF*Radius_Nozzle_Throat);
	Bnd_Spline_East.Create_Spline_Line(xc_SE,xc_NE,2);
	Bnd_Spline_West.Create_Spline_Line(xc_SW,xc_NW,2);
	Bnd_Spline_South.Create_Spline_Line(xc_SW,xc_SE,2);
	Bnd_Spline_North.Create_Spline_Diverging_Nozzle(Length_Chamber_To_Throat,
							Length_Nozzle,
							HALF*Radius_Nozzle_Throat,
							HALF*Radius_Nozzle_Exit,
							251);
	Bnd_Spline_North.setBCtype(BC_NONE);
	Bnd_Spline_South.setBCtype(BC_REFLECTION);
	Bnd_Spline_East.setBCtype(BC_CONSTANT_EXTRAPOLATION);
	Bnd_Spline_East.setBCtype(BC_CHARACTERISTIC);
	Bnd_Spline_West.setBCtype(BC_NONE);
	if (Stretching_Flag) {
	  Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
	  Beta_I = Stretching_Factor_Idir;
	  Tau_I = ZERO;
	  Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
	  Beta_J = Stretching_Factor_Jdir;
	  Tau_J = ZERO;
	}
	block_flag = ON;

      } else if (jBlk == 1 && iBlk - iBlk_Chamber == 0 && iBlk_Nozzle == 2) {
	xc_SW = Vector2D(ZERO,Radius_Grain);
	xc_SE = Vector2D(Length_Chamber_To_Throat,HALF*Radius_Nozzle_Throat);
	xc_NE = Vector2D(Length_Chamber_To_Throat,Radius_Nozzle_Throat);
	xc_NW = Vector2D(ZERO,Radius_Chamber);
	Bnd_Spline_East.Create_Spline_Line(xc_SE,xc_NE,2);
	Bnd_Spline_West.Create_Spline_Line(xc_SW,xc_NW,2);
	Bnd_Spline_South.Create_Spline_Converging_Nozzle(ZERO,
							 Length_Chamber_To_Throat,
							 Radius_Grain,
							 HALF*Radius_Nozzle_Throat,
							 251);
	Bnd_Spline_North.Create_Spline_Converging_Nozzle(ZERO,
							 Length_Chamber_To_Throat,
							 Radius_Chamber,
							 Radius_Nozzle_Throat,
							 251);
	Bnd_Spline_North.setBCtype(Chamber_BC_Type);
	Bnd_Spline_South.setBCtype(BC_NONE);
	Bnd_Spline_East.setBCtype(BC_NONE);
	Bnd_Spline_West.setBCtype(BC_NONE);
	if (Stretching_Flag) {
	  Stretch_I = STRETCHING_FCN_MINMAX_CLUSTERING;
	  Beta_I = Stretching_Factor_Idir;
	  Tau_I = ZERO;
	  Stretch_J = STRETCHING_FCN_MAX_CLUSTERING;
	  Beta_J = Stretching_Factor_Jdir;
	  Tau_J = ZERO;
	}
	block_flag = ON;

      } else if (jBlk == 1 && iBlk - iBlk_Chamber == 1 && iBlk_Nozzle == 2) {
	xc_SW = Vector2D(Length_Chamber_To_Throat,HALF*Radius_Nozzle_Throat);
	xc_SE = Vector2D(Length_Nozzle,HALF*Radius_Nozzle_Exit);
	xc_NE = Vector2D(Length_Nozzle,Radius_Nozzle_Exit);
	xc_NW = Vector2D(Length_Chamber_To_Throat,Radius_Nozzle_Throat);
	Bnd_Spline_East.Create_Spline_Line(xc_SE,xc_NE,2);
	Bnd_Spline_West.Create_Spline_Line(xc_SW,xc_NW,2);
	Bnd_Spline_South.Create_Spline_Diverging_Nozzle(Length_Chamber_To_Throat,
							Length_Nozzle,
							HALF*Radius_Nozzle_Throat,
							HALF*Radius_Nozzle_Exit,
							251);
	Bnd_Spline_North.Create_Spline_Diverging_Nozzle(Length_Chamber_To_Throat,
							Length_Nozzle,
							Radius_Nozzle_Throat,
							Radius_Nozzle_Exit,
							251);
	Bnd_Spline_North.setBCtype(Chamber_BC_Type);
	Bnd_Spline_South.setBCtype(BC_NONE);
	Bnd_Spline_East.setBCtype(BC_CONSTANT_EXTRAPOLATION);
	Bnd_Spline_East.setBCtype(BC_CHARACTERISTIC);
	Bnd_Spline_West.setBCtype(BC_NONE);
	if (Stretching_Flag) {
	  Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
	  Beta_I = Stretching_Factor_Idir;
	  Tau_I = ZERO;
	  Stretch_J = STRETCHING_FCN_MAX_CLUSTERING;
	  Beta_J = Stretching_Factor_Jdir;
	  Tau_J = ZERO;
	}
	block_flag = ON;

      } else {
	block_flag = OFF;

      }

      if (block_flag) {
	Grid_ptr[iBlk][jBlk].Create_Quad_Block_Without_Update(Bnd_Spline_North,
							      Bnd_Spline_South,
							      Bnd_Spline_East,
							      Bnd_Spline_West,
							      Number_of_Cells_Idir,
							      Number_of_Cells_Jdir,
							      Number_of_Ghost_Cells,
							      Highest_Order_of_Reconstruction,
							      GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
							      Stretch_I,
							      Beta_I,
							      Tau_I,
							      Stretch_J,
							      Beta_J,
							      Tau_J,
							      Orthogonal_North,
							      Orthogonal_South,
							      Orthogonal_East,
							      Orthogonal_West);

 	// Smooth the 2D quadrilateral grid block.
	// 	Smooth_Quad_Block(Grid_ptr[iBlk][jBlk],
	// 			  min(250,2*max(Number_of_Cells_Idir,
	// 					Number_of_Cells_Jdir)));
	Grid_ptr[iBlk][jBlk].Smooth_Rocket_Motor(Length_Chamber,
						 Radius_Chamber,
						 Length_Chamber_To_Throat,
						 Length_Nozzle,
						 Radius_Nozzle_Exit,
						 Radius_Nozzle_Throat,
						 Radius_Grain,
						 Nozzle_Type,
						 Stretching_Factor_Idir,
						 Stretching_Factor_Jdir,
						 -1,
						 0,
						 0,0,
						 iBlk,jBlk,
						 Number_of_Blocks_Idir,
						 Number_of_Blocks_Jdir);

      }

    }
  }


}

/*!
 * Routine: Grid_Nozzleless_Rocket_Motor                              
 * This subroutine DOESN'T update the ghost cells or
 * the geometric properties of the grid cells.
 */
void Grid2D_Quad_MultiBlock_HO::Grid_Nozzleless_Rocket_Motor_Without_Update(int &_Number_of_Blocks_Idir_,
									    int &_Number_of_Blocks_Jdir_,
									    const double &Length_Chamber,
									    const double &Radius_Chamber,
									    const double &Length_Nozzle,
									    const double &Radius_Nozzle_Exit,
									    const int &Chamber_BC_Type,
									    const int &Stretching_Flag,
									    const int Stretching_Type_Jdir,
									    const double &Stretching_Factor_Idir,
									    const double &Stretching_Factor_Jdir,
									    const int Number_of_Cells_Idir,
									    const int Number_of_Cells_Jdir,
									    const int Number_of_Ghost_Cells,
									    const int Highest_Order_of_Reconstruction) {
  
  int Stretch_I, Stretch_J,
      Orthogonal_North, Orthogonal_South,
      Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D xc_NW, xc_NE, xc_SE, xc_SW;
  Spline2D_HO Bnd_Spline_North, Bnd_Spline_South,
              Bnd_Spline_East, Bnd_Spline_West;

  // Allocate memory for grid blocks.
  _Number_of_Blocks_Jdir_ = 1;
  allocate(_Number_of_Blocks_Idir_,_Number_of_Blocks_Jdir_);

  // Create the mesh for each block representing the complete grid.
  for (int jBlk = 0; jBlk < Number_of_Blocks_Jdir; jBlk++) {
    for (int iBlk = 0; iBlk < Number_of_Blocks_Idir; iBlk++) {
      
      // Default the stretching parameters.
      Stretch_I = STRETCHING_FCN_LINEAR;
      Beta_I = ONE;
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_LINEAR;
      Beta_J = ONE;
      Tau_J = ZERO;
      Orthogonal_North = 1;
      Orthogonal_South = 1;
      Orthogonal_East = 1;
      Orthogonal_West = 1;
      
      // Create the splines defining the north, south, east, and west 
      // boundaries of the grid block.
      if (iBlk < Number_of_Blocks_Idir-1) {//2) {
	// 	xc_SW = Vector2D((Length_Chamber-0.002)*double(iBlk)/double(Number_of_Blocks_Idir-2),ZERO);
	// 	xc_SE = Vector2D((Length_Chamber-0.002)*double(iBlk+1)/double(Number_of_Blocks_Idir-2),ZERO);
	// 	xc_NW = Vector2D((Length_Chamber-0.002)*double(iBlk)/double(Number_of_Blocks_Idir-2),Radius_Chamber);
	// 	xc_NE = Vector2D((Length_Chamber-0.002)*double(iBlk+1)/double(Number_of_Blocks_Idir-2),Radius_Chamber);
	xc_SW = Vector2D(Length_Chamber*double(iBlk)/double(Number_of_Blocks_Idir-1),ZERO);
	xc_SE = Vector2D(Length_Chamber*double(iBlk+1)/double(Number_of_Blocks_Idir-1),ZERO);
	xc_NW = Vector2D(Length_Chamber*double(iBlk)/double(Number_of_Blocks_Idir-1),Radius_Chamber);
	xc_NE = Vector2D(Length_Chamber*double(iBlk+1)/double(Number_of_Blocks_Idir-1),Radius_Chamber);
	Bnd_Spline_North.Create_Spline_Line(xc_NW,xc_NE,2);
	Bnd_Spline_South.Create_Spline_Line(xc_SW,xc_SE,2);
	Bnd_Spline_East.Create_Spline_Line(xc_SE,xc_NE,2);
	Bnd_Spline_West.Create_Spline_Line(xc_SW,xc_NW,2);
	Bnd_Spline_North.setBCtype(BC_MASS_INJECTION);
	Bnd_Spline_South.setBCtype(BC_REFLECTION);
	Bnd_Spline_East.setBCtype(BC_NONE);
	if (iBlk == 0) Bnd_Spline_West.setBCtype(Chamber_BC_Type);
	else Bnd_Spline_West.setBCtype(BC_NONE);
	//       } else if (iBlk == Number_of_Blocks_Idir-2) {
	// 	xc_SW = Vector2D(Length_Chamber-0.002,ZERO);
	// 	xc_SE = Vector2D(Length_Chamber,ZERO);
	// 	xc_NW = Vector2D(Length_Chamber-0.002,Radius_Chamber);
	// 	xc_NE = Vector2D(Length_Chamber,Radius_Chamber);
	// 	Create_Spline_Line(Bnd_Spline_North,xc_NW,xc_NE,2);
	// 	Create_Spline_Line(Bnd_Spline_South,xc_SW,xc_SE,2);
	// 	Create_Spline_Line(Bnd_Spline_East,xc_SE,xc_NE,2);
	// 	Create_Spline_Line(Bnd_Spline_West,xc_SW,xc_NW,2);
	// 	Bnd_Spline_North.setBCtype(Chamber_BC_Type);
	// 	Bnd_Spline_South.setBCtype(BC_REFLECTION);
	// 	Bnd_Spline_East.setBCtype(BC_NONE);
	// 	if (iBlk == 0) Bnd_Spline_West.setBCtype(Chamber_BC_Type);
	// 	else Bnd_Spline_West.setBCtype(BC_NONE);
      } else {
	xc_SW = Vector2D(Length_Chamber,ZERO);
	xc_SE = Vector2D(Length_Chamber+Length_Nozzle,ZERO);
	xc_NW = Vector2D(Length_Chamber,Radius_Chamber);
	xc_NE = Vector2D(Length_Chamber+Length_Nozzle,Radius_Nozzle_Exit);
	Bnd_Spline_North.Create_Spline_Line(xc_NW,xc_NE,2);
	Bnd_Spline_South.Create_Spline_Line(xc_SW,xc_SE,2);
	Bnd_Spline_East.Create_Spline_Line(xc_SE,xc_NE,2);
	Bnd_Spline_West.Create_Spline_Line(xc_SW,xc_NW,2);
	Bnd_Spline_North.setBCtype(Chamber_BC_Type);
	Bnd_Spline_South.setBCtype(BC_REFLECTION);
	//Bnd_Spline_East.setBCtype(BC_CHARACTERISTIC);
	//Bnd_Spline_East.setBCtype(BC_FIXED_PRESSURE);
	Bnd_Spline_East.setBCtype(BC_CONSTANT_EXTRAPOLATION);
	Bnd_Spline_West.setBCtype(BC_NONE);
      }
      if (Stretching_Flag) {
	Stretch_J = Stretching_Type_Jdir;
	Beta_J = Stretching_Factor_Jdir;
	Tau_J = ZERO;
	if (iBlk == 0) {
	  Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
	  Beta_I = Stretching_Factor_Idir;
	} else if (iBlk == Number_of_Blocks_Idir-2) {//3) {
	  Stretch_I = STRETCHING_FCN_MAX_CLUSTERING;
	  Beta_I = Stretching_Factor_Idir;
	} else if (iBlk == Number_of_Blocks_Idir-1) {
	  Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
	  Beta_I = Stretching_Factor_Idir;
	}
      }
      
      // Create the quadrilateral solution block.
      Grid_ptr[iBlk][jBlk].Create_Quad_Block_Without_Update(Bnd_Spline_North,
							    Bnd_Spline_South,
							    Bnd_Spline_East,
							    Bnd_Spline_West,
							    Number_of_Cells_Idir,
							    Number_of_Cells_Jdir,
							    Number_of_Ghost_Cells,
							    Highest_Order_of_Reconstruction,
							    GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
							    Stretch_I,
							    Beta_I,
							    Tau_I,
							    Stretch_J,
							    Beta_J,
							    Tau_J,
							    Orthogonal_North,
							    Orthogonal_South,
							    Orthogonal_East,
							    Orthogonal_West);

    }
  }

}

/*!
 * Generates a two-block grid for a rocket motor nozzle.              
 * This subroutine DOESN'T update the ghost cells or
 * the geometric properties of the grid cells.
 */
void Grid2D_Quad_MultiBlock_HO::Grid_Nozzle_Without_Update(int &_Number_of_Blocks_Idir_,
							   int &_Number_of_Blocks_Jdir_,
							   const double &Length_Nozzle,
							   const double &Radius_Chamber,
							   const double &Radius_Nozzle_Exit,
							   const double &Radius_Nozzle_Throat,
							   const int &Nozzle_Type,
							   const int &Stretching_Flag,
							   const int &Stretching_Type_Idir,
							   const int &Stretching_Type_Jdir,
							   const double &Stretching_Factor_Idir,
							   const double &Stretching_Factor_Jdir,
							   const int Number_of_Cells_Idir,
							   const int Number_of_Cells_Jdir,
							   const int Number_of_Ghost_Cells,
							   const int Highest_Order_of_Reconstruction) {
  
  int error_flag, block_flag;
  int Stretch_I, Stretch_J,
      Orthogonal_North, Orthogonal_South,
      Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D xc_NW, xc_NE, xc_SE, xc_SW;
  Spline2D_HO Bnd_Spline_North, Bnd_Spline_South,
              Bnd_Spline_East, Bnd_Spline_West;

  // Allocate memory for grid blocks.
  if (!Nozzle_Type) _Number_of_Blocks_Idir_ = 1;
  else _Number_of_Blocks_Idir_ = 2;
  _Number_of_Blocks_Jdir_ = 1;
  allocate(_Number_of_Blocks_Idir_,_Number_of_Blocks_Jdir_);

  // Create the mesh for each block representing the complete grid.
  for (int iBlk = 0; iBlk < Number_of_Blocks_Idir; iBlk++) {

    if (iBlk == 0) {
      // Converging section of the nozzle.
      xc_SW = Vector2D(ZERO,ZERO);
      xc_SE = Vector2D(HALF*Length_Nozzle,ZERO);
      xc_NE = Vector2D(HALF*Length_Nozzle,Radius_Nozzle_Throat);
      xc_NW = Vector2D(ZERO,Radius_Chamber);
      Bnd_Spline_North.Create_Spline_Converging_Nozzle(ZERO,
						       HALF*Length_Nozzle,
						       Radius_Chamber,
						       Radius_Nozzle_Throat,
						       251);
      Bnd_Spline_South.Create_Spline_Line(xc_SW,xc_SE,2);
      Bnd_Spline_East.Create_Spline_Line(xc_SE,xc_NE,2);
      Bnd_Spline_West.Create_Spline_Line(xc_SW,xc_NW,2);
      Bnd_Spline_North.setBCtype(BC_REFLECTION);
      Bnd_Spline_South.setBCtype(BC_REFLECTION);
      Bnd_Spline_East.setBCtype(BC_NONE);
      Bnd_Spline_West.setBCtype(BC_FIXED);
    } else {
      // Diverging section of the nozzle.
      xc_SW = Vector2D(HALF*Length_Nozzle,ZERO);
      xc_SE = Vector2D(Length_Nozzle,ZERO);
      xc_NE = Vector2D(Length_Nozzle,Radius_Nozzle_Exit);
      xc_NW = Vector2D(HALF*Length_Nozzle,Radius_Nozzle_Throat);
      Bnd_Spline_North.Create_Spline_Diverging_Nozzle(HALF*Length_Nozzle,
						      Length_Nozzle,
						      Radius_Nozzle_Throat,
						      Radius_Nozzle_Exit,
						      251);
      Bnd_Spline_South.Create_Spline_Line(xc_SW,xc_SE,2);
      Bnd_Spline_East.Create_Spline_Line(xc_SE,xc_NE,2);
      Bnd_Spline_West.Create_Spline_Line(xc_SW,xc_NW,2);
      Bnd_Spline_North.setBCtype(BC_REFLECTION);
      Bnd_Spline_South.setBCtype(BC_REFLECTION);
      Bnd_Spline_East.setBCtype(BC_FIXED);
      Bnd_Spline_West.setBCtype(BC_NONE);
    }

    // Assign values to the stretching function parameters and 
    // boundary grid line orthogonality parameters.
    Stretch_I = STRETCHING_FCN_LINEAR; Beta_I = ONE; Tau_I = ZERO;
    Stretch_J = STRETCHING_FCN_LINEAR; Beta_J = ONE; Tau_J = ZERO;
    if (Stretching_Flag) {
      if (Stretching_Type_Idir != STRETCHING_FCN_LINEAR) {
	if (iBlk == 0) Stretch_I = STRETCHING_FCN_MAX_CLUSTERING;
	else Stretch_I = STRETCHING_FCN_MAX_CLUSTERING;
	Beta_I = Stretching_Factor_Idir;
      }
      Stretch_J = Stretching_Type_Jdir;
      Beta_J = Stretching_Factor_Jdir;
    }
    Orthogonal_North = 1;
    Orthogonal_South = 1;
    Orthogonal_East = 1;
    Orthogonal_West = 1;

    // Create the 2D quadrilateral grid block.
    Grid_ptr[iBlk][0].Create_Quad_Block_Without_Update(Bnd_Spline_North,
						       Bnd_Spline_South,
						       Bnd_Spline_East,
						       Bnd_Spline_West,
						       Number_of_Cells_Idir,
						       Number_of_Cells_Jdir,
						       Number_of_Ghost_Cells,
						       Highest_Order_of_Reconstruction,
						       GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
						       Stretch_I,
						       Beta_I,
						       Tau_I,
						       Stretch_J,
						       Beta_J,
						       Tau_J,
						       Orthogonal_North,
						       Orthogonal_South,
						       Orthogonal_East,
						       Orthogonal_West);

    // Smooth the 2D quadrilateral grid block.
    Grid_ptr[iBlk][0].Smooth_Quad_Block(min(250,2*max(Number_of_Cells_Idir,
						      Number_of_Cells_Jdir)));
    
  }
}

/*!
 * Generates a double-block O-type grid for predicting  
 * flow past a circular cylinder.                       
 *                                                      
 * Usage: Grid_ptr = Grid_Circular_Cylinder(Grid_ptr,   
 *                                          nblk_i,     
 *                                          nblk_j,     
 *                                          THREE,      
 *   		                            100,        
 *  		                            50,         
 *                                          2);FIXME         
 * This subroutine DOESN'T update the ghost cells or
 * the geometric properties of the grid cells.
 */
void Grid2D_Quad_MultiBlock_HO::Grid_Circular_Cylinder_Without_Update(int &_Number_of_Blocks_Idir_,
								      int &_Number_of_Blocks_Jdir_,
								      const double &Radius,
								      const int Stretching_Type_Idir,
								      const int Stretching_Type_Jdir,
								      const double &Stretching_Factor_Idir,
								      const double &Stretching_Factor_Jdir,
								      const int Number_of_Cells_Idir,
								      const int Number_of_Cells_Jdir,
								      const int Number_of_Ghost_Cells,
								      const int Highest_Order_of_Reconstruction) {
  
  Grid_Circular_Cylinder_Without_Update(_Number_of_Blocks_Idir_,
					_Number_of_Blocks_Jdir_,
					Radius,
					32.00*Radius,
					Stretching_Type_Idir,
					Stretching_Type_Jdir,
					Stretching_Factor_Idir,
					Stretching_Factor_Jdir,
					Number_of_Cells_Idir,
					Number_of_Cells_Jdir,
					Number_of_Ghost_Cells,
					Highest_Order_of_Reconstruction);
  
}

void Grid2D_Quad_MultiBlock_HO::Grid_Circular_Cylinder_Without_Update(int &_Number_of_Blocks_Idir_,
								      int &_Number_of_Blocks_Jdir_,
								      const double &Inner_Radius,
								      const double &Outer_Radius,
								      const int Stretching_Type_Idir,
								      const int Stretching_Type_Jdir,
								      const double &Stretching_Factor_Idir,
								      const double &Stretching_Factor_Jdir,
								      const int Number_of_Cells_Idir,
								      const int Number_of_Cells_Jdir,
								      const int Number_of_Ghost_Cells,
								      const int Highest_Order_of_Reconstruction) {
  
  int iBlk, n_cells_i, n_cells_j, Stretch_I, Stretch_J,
    Orthogonal_North, Orthogonal_South,
    Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D x1, x2;
  Spline2D_HO Bnd_Spline_North, Bnd_Spline_South,
              Bnd_Spline_East, Bnd_Spline_West;
  
  /* Allocate memory for grid blocks.  There are two grid
     blocks for this mesh. */
  _Number_of_Blocks_Idir_ = 2;
  _Number_of_Blocks_Jdir_ = 1;
  allocate(_Number_of_Blocks_Idir_, _Number_of_Blocks_Jdir_);

  /* Create the mesh for each block representing
     the complete grid. */
  
  for ( iBlk = 0; iBlk <= Number_of_Blocks_Idir-1; ++iBlk ) {

    /* Create the splines defining the north, south,
       east, and west boundaries of the grid. */

    if (iBlk == 0) {
      x1 = Vector2D(ZERO,ZERO);
      Bnd_Spline_North.Create_Spline_Circular_Arc(x1,
						  Outer_Radius,
						  360.00,
						  180.00,
						  361);
      Bnd_Spline_South.Create_Spline_Circular_Arc(x1,
						  Inner_Radius,
						  360.00,
						  180.00,
						  361);
      x1 = Vector2D(Inner_Radius, ZERO);
      x2 = Vector2D(Outer_Radius, ZERO);
      Bnd_Spline_West.Create_Spline_Line(x1, x2, 2);
      x1 = Vector2D(-Inner_Radius, ZERO);
      x2 = Vector2D(-Outer_Radius, ZERO);
      Bnd_Spline_East.Create_Spline_Line(x1, x2, 2);
    } else {
      x1 = Vector2D(ZERO,ZERO);
      Bnd_Spline_North.Create_Spline_Circular_Arc(x1,
						  Outer_Radius,
						  180.00,
						  ZERO,
						  361);
      Bnd_Spline_South.Create_Spline_Circular_Arc(x1,
						  Inner_Radius,
						  180.00,
						  ZERO,
						  361);
      x1 = Vector2D(-Inner_Radius, ZERO);
      x2 = Vector2D(-Outer_Radius, ZERO);
      Bnd_Spline_West.Create_Spline_Line(x1, x2, 2);
      x1 = Vector2D(Inner_Radius, ZERO);
      x2 = Vector2D(Outer_Radius, ZERO);
      Bnd_Spline_East.Create_Spline_Line(x1, x2, 2);
    } /* endif */

    /* Set the boundary condition types for each of the
       boundary splines. */

    if (iBlk == 0) {
      Bnd_Spline_North.setBCtype(BC_FIXED);
      Bnd_Spline_South.setBCtype(BC_REFLECTION);
      Bnd_Spline_South.makeSplineSolidBoundary(); // create first solid body
      Bnd_Spline_East.setBCtype(BC_NONE);
      Bnd_Spline_West.setBCtype(BC_NONE);
    } else {      
      Bnd_Spline_North.setBCtype(BC_FIXED);
      Bnd_Spline_South.setBCtype(BC_REFLECTION);
      Bnd_Spline_South.makeSplineSolidBoundary(1); // make this spline part of the first solid body
      Bnd_Spline_East.setBCtype(BC_NONE);
      Bnd_Spline_West.setBCtype(BC_NONE);
    } /* endif */

    /* Determine the number of cells for this block. */

    n_cells_i = Number_of_Cells_Idir/2;
    n_cells_j = Number_of_Cells_Jdir;

    /* Assign values to the stretching function parameters
       and boundary grid line orthogonality parameters. */

    if (iBlk == 0) {
      Stretch_I = Stretching_Type_Idir;
      Beta_I = Stretching_Factor_Idir;
      Tau_I = ZERO;
      Stretch_J = Stretching_Type_Jdir;
      Beta_J = Stretching_Factor_Jdir;
      Tau_J = ZERO;
    } else {
      Stretch_I = Stretching_Type_Idir;
      Beta_I = Stretching_Factor_Idir;
      Tau_I = ZERO;
      Stretch_J = Stretching_Type_Jdir;
      Beta_J = Stretching_Factor_Jdir;
      Tau_J = ZERO;
    } /* endif */

    Orthogonal_North = 0;
    Orthogonal_South = 0;
    Orthogonal_East = 0;
    Orthogonal_West = 0;

    /* Create the 2D quadrilateral grid block. */

    Grid_ptr[iBlk][0].Create_Quad_Block_Without_Update(Bnd_Spline_North,
						       Bnd_Spline_South,
						       Bnd_Spline_East,
						       Bnd_Spline_West,
						       n_cells_i,
						       n_cells_j,
						       Number_of_Ghost_Cells,
						       Highest_Order_of_Reconstruction,
						       GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
						       Stretch_I,
						       Beta_I, 
						       Tau_I,
						       Stretch_J,
						       Beta_J,
						       Tau_J,
						       Orthogonal_North,
						       Orthogonal_South,
						       Orthogonal_East,
						       Orthogonal_West);

  } /* endfor */

}

/*!
 * Generates a single-block O-type grid for predicting     
 * flow through an annulus of given inner and outer radii. 
 * This subroutine DOESN'T update the ghost cells or
 * the geometric properties of the grid cells.
 *                                                      
 * \param ThetaStart the start angle of the annulus (east edge) measured counterclockwise
 * \param ThetaEnd   the end angle of the annulus (west edge) measured counterclockwise
 */
void Grid2D_Quad_MultiBlock_HO::Grid_Annulus_Without_Update(int &_Number_of_Blocks_Idir_,
							    int &_Number_of_Blocks_Jdir_,
							    const double &Inner_Radius,
							    const double &Outer_Radius,
							    const double &ThetaStart,
							    const double &ThetaEnd,
							    const int Stretching_Type_Idir,
							    const int Stretching_Type_Jdir,
							    const double &Stretching_Factor_Idir,
							    const double &Stretching_Factor_Jdir,
							    const int Number_of_Cells_Idir,
							    const int Number_of_Cells_Jdir,
							    const int Number_of_Ghost_Cells,
							    const int Highest_Order_of_Reconstruction) {
  
  int n_cells_i, n_cells_j, Stretch_I, Stretch_J,
    Orthogonal_North, Orthogonal_South,
    Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D x1, x2;
  Spline2D_HO Bnd_Spline_North, Bnd_Spline_South,
              Bnd_Spline_East, Bnd_Spline_West;

  /* Allocate memory for grid blocks.  There is one grid
     blocks for this mesh. */
  
  _Number_of_Blocks_Idir_ = 1;
  _Number_of_Blocks_Jdir_ = 1;
  allocate(_Number_of_Blocks_Idir_, _Number_of_Blocks_Jdir_);
  
  /* Create the mesh for the block. */
  
  /* Create the splines defining the north, south,
     east, and west boundaries of the grid. */
  
  x1 = Vector2D(ZERO,ZERO);
  Bnd_Spline_North.Create_Spline_Circular_Arc(x1,
					      Outer_Radius,
					      ThetaEnd,
					      ThetaStart,
					      int(2*fabs(ThetaEnd-ThetaStart)+1)); // consider 2 points per degree
  Bnd_Spline_South.Create_Spline_Circular_Arc(x1,
					      Inner_Radius,
					      ThetaEnd,
					      ThetaStart,
					      int(2*fabs(ThetaEnd-ThetaStart)+1)); // consider 2 points per degree

  Bnd_Spline_East.Create_Spline_Line_Polar_Coordinates(Inner_Radius, Outer_Radius, ThetaStart, 2);
  Bnd_Spline_West.Create_Spline_Line_Polar_Coordinates(Inner_Radius, Outer_Radius, ThetaEnd, 2);

  /* Set the boundary condition types for each of the
     boundary splines. */
    
  Bnd_Spline_North.setBCtype(BC_DIRICHLET);
  Bnd_Spline_South.setBCtype(BC_DIRICHLET);
  Bnd_Spline_East.setBCtype(BC_INFLOW);	 // inflow BC is treated as a Dirichlet boundary condition
  Bnd_Spline_West.setBCtype(BC_OUTFLOW); // outflow BC is treated as a Neumann boundary condition

  /* Determine the number of cells for this block. */
    
  n_cells_i = Number_of_Cells_Idir;
  n_cells_j = Number_of_Cells_Jdir;
    
  /* Assign values to the stretching function parameters
     and boundary grid line orthogonality parameters. */
    
  Stretch_I = Stretching_Type_Idir;
  Beta_I = Stretching_Factor_Idir;
  Tau_I = ZERO;
  Stretch_J = Stretching_Type_Jdir;
  Beta_J = Stretching_Factor_Jdir;
  Tau_J = ZERO;

  Orthogonal_North = 0;
  Orthogonal_South = 0;
  Orthogonal_East = 0;
  Orthogonal_West = 0;
    
  /* Create the 2D quadrilateral grid block. */

  Grid_ptr[0][0].Create_Quad_Block_Without_Update(Bnd_Spline_North,
						  Bnd_Spline_South,
						  Bnd_Spline_East,
						  Bnd_Spline_West,
						  n_cells_i,
						  n_cells_j,
						  Number_of_Ghost_Cells,
						  Highest_Order_of_Reconstruction,
						  GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
						  Stretch_I,
						  Beta_I, 
						  Tau_I,
						  Stretch_J,
						  Beta_J,
						  Tau_J,
						  Orthogonal_North,
						  Orthogonal_South,
						  Orthogonal_East,
						  Orthogonal_West);

}

/*!
 * Generates a double-block O-type grid for predicting  
 * flow past an ellipse.                                
 *                                                      
 * Usage: Grid_ptr = Grid_Ellipse(Grid_ptr,             
 *                                nblk_i,               
 *                                nblk_j,               
 *                                FOUR,                 
 *                                ONE,                  
 *   		                  100,                  
 *  		                  50,                   
 *  		                  2);                   
 * This subroutine DOESN'T update the ghost cells or
 * the geometric properties of the grid cells.
 */
void Grid2D_Quad_MultiBlock_HO::Grid_Ellipse_Without_Update(int &_Number_of_Blocks_Idir_,
							    int &_Number_of_Blocks_Jdir_,
							    const double &A,
							    const double &B,
							    const int Number_of_Cells_Idir,
							    const int Number_of_Cells_Jdir,
							    const int Number_of_Ghost_Cells,
							    const int Highest_Order_of_Reconstruction) {
  
  int iBlk, n_cells_i, n_cells_j, Stretch_I, Stretch_J,
    Orthogonal_North, Orthogonal_South,
    Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D x1, x2;
  Spline2D_HO Bnd_Spline_North, Bnd_Spline_South,
              Bnd_Spline_East, Bnd_Spline_West;

  /* Allocate memory for grid blocks.  There are two grid
     blocks for this mesh. */
  
  _Number_of_Blocks_Idir_ = 2;
  _Number_of_Blocks_Jdir_ = 1;
  allocate(_Number_of_Blocks_Idir_, _Number_of_Blocks_Jdir_);

  /* Create the mesh for each block representing
     the complete grid. */
  
  for ( iBlk = 0; iBlk <= Number_of_Blocks_Idir-1; ++iBlk ) {
    
    /* Create the splines defining the north, south,
       east, and west boundaries of the grid. */
    
    if (iBlk == 0) {
      x1 = Vector2D(ZERO , ZERO);
      Bnd_Spline_North.Create_Spline_Circular_Arc(x1,
						  32.00*A,
						  360.00,
						  180.00,
						  361);
      Bnd_Spline_South.Create_Spline_Ellipsoidal_Arc(x1,
						     A,
						     B,
						     360.00,
						     180.00,
						     361);
      x1 = Vector2D(A, ZERO);
      x2 = Vector2D(32.00*A, ZERO);
      Bnd_Spline_West.Create_Spline_Line(x1, x2, 2);
      x1 = Vector2D(-A, ZERO);
      x2 = Vector2D(-32.00*A, ZERO);
      Bnd_Spline_East.Create_Spline_Line(x1, x2, 2);
    } else {
      x1 = Vector2D(ZERO , ZERO);
      Bnd_Spline_North.Create_Spline_Circular_Arc(x1,
						  32.00*A,
						  180.00,
						  ZERO,
						  361);
      Bnd_Spline_South.Create_Spline_Ellipsoidal_Arc(x1,
						     A,
						     B,
						     180.00,
						     ZERO,
						     361);
      x1 = Vector2D(-A, ZERO);
      x2 = Vector2D(-32.00*A, ZERO);
      Bnd_Spline_West.Create_Spline_Line(x1, x2, 2);
      x1 = Vector2D(A, ZERO);
      x2 = Vector2D(32.00*A, ZERO);
      Bnd_Spline_East.Create_Spline_Line(x1, x2, 2);
    } /* endif */

    /* Set the boundary condition types for each of the
       boundary splines. */

    if (iBlk == 0) {
      Bnd_Spline_North.setBCtype(BC_FIXED);
      Bnd_Spline_South.setBCtype(BC_REFLECTION);
      Bnd_Spline_East.setBCtype(BC_NONE);
      Bnd_Spline_West.setBCtype(BC_NONE);
    } else {
      Bnd_Spline_North.setBCtype(BC_FIXED);
      Bnd_Spline_South.setBCtype(BC_REFLECTION);
      Bnd_Spline_East.setBCtype(BC_NONE);
      Bnd_Spline_West.setBCtype(BC_NONE);
    } /* endif */

    /* Determine the number of cells for this block. */

    n_cells_i = Number_of_Cells_Idir/2;
    n_cells_j = Number_of_Cells_Jdir;

    /* Assign values to the stretching function parameters
       and boundary grid line orthogonality parameters. */

    if (iBlk == 0) {
      Stretch_I = STRETCHING_FCN_MINMAX_CLUSTERING;
      Beta_I = 1.025; 
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
      Beta_J = 1.001;
      Tau_J = ZERO;
      Orthogonal_North = 1;
      Orthogonal_South = 1;
      Orthogonal_East = 1;
      Orthogonal_West = 1;
    } else {
      Stretch_I = STRETCHING_FCN_MINMAX_CLUSTERING;
      Beta_I = 1.025; 
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
      Beta_J = 1.001;
      Tau_J = ZERO;
      Orthogonal_North = 1;
      Orthogonal_South = 1;
      Orthogonal_East = 1;
      Orthogonal_West = 1;
    } /* endif */

    /* Create the 2D quadrilateral grid block. */

    Grid_ptr[iBlk][0].Create_Quad_Block_Without_Update(Bnd_Spline_North,
						       Bnd_Spline_South,
						       Bnd_Spline_East,
						       Bnd_Spline_West,
						       n_cells_i,
						       n_cells_j,
						       Number_of_Ghost_Cells,
						       Highest_Order_of_Reconstruction,
						       GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
						       Stretch_I,
						       Beta_I, 
						       Tau_I,
						       Stretch_J,
						       Beta_J,
						       Tau_J,
						       Orthogonal_North,
						       Orthogonal_South,
						       Orthogonal_East,
						       Orthogonal_West);

    /* Smooth the 2D quadrilateral grid block. */

    Grid_ptr[iBlk][0].Smooth_Quad_Block(min(250, 2*max(n_cells_i,n_cells_j)));

  } /* endfor */
}

/*!
 * Generates a C-type grid consisting of four           
 * quadrilateral grid blocks for predicting flow past   
 * NACA 4-digit and 5-digit aerofoils.                  
 *                                                      
 * Usage:  Grid_ptr = Grid_NACA_Aerofoil(Grid_ptr,      
 *                                       nblk_i,        
 *                                       nblk_j,        
 *                                       "4412",        
 *                                       ONE,           
 *		                         120,           
 *		                         50,            
 *                                       2);            
 * This subroutine DOESN'T update the ghost cells or
 * the geometric properties of the grid cells.
 */
void Grid2D_Quad_MultiBlock_HO::Grid_NACA_Aerofoil_Without_Update(int &_Number_of_Blocks_Idir_,
								  int &_Number_of_Blocks_Jdir_,
								  char *NACA_Aerofoil_Type_ptr,
								  const double &Chord_Length,
								  const int Number_of_Cells_Idir,
								  const int Number_of_Cells_Jdir,
								  const int Number_of_Ghost_Cells,
								  const int Highest_Order_of_Reconstruction) {
  
  int iBlk, n_cells_i, n_cells_j,
    Stretch_I, Stretch_J,
    Orthogonal_North, Orthogonal_South,
    Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D x1, x2;
  Spline2D_HO Bnd_Spline_North, Bnd_Spline_South,
              Bnd_Spline_East, Bnd_Spline_West,
              s_tmp1, s_tmp2;

  /* Allocate memory for grid blocks.  There are four grid
     blocks for this mesh. */

  _Number_of_Blocks_Idir_ = 4;
  _Number_of_Blocks_Jdir_ = 1;
  allocate(_Number_of_Blocks_Idir_, _Number_of_Blocks_Jdir_);

  /* Create the mesh for each block representing
     the complete grid. */

  for ( iBlk = 0; iBlk <= Number_of_Blocks_Idir-1; ++iBlk ) {

    /* Create the splines defining the north, south,
       east, and west boundaries of the grid. */

    if (iBlk == 0) {
      x1 = Vector2D(32.00*Chord_Length, -32.00*Chord_Length);
      x2 = Vector2D(Chord_Length, -32.00*Chord_Length);
      Bnd_Spline_North.Create_Spline_Line(x1, x2, 2);
      x1 = Vector2D(32.00*Chord_Length, ZERO);
      x2 = Vector2D(Chord_Length, ZERO);
      Bnd_Spline_South.Create_Spline_Line(x1, x2, 2);
      x1 = Vector2D(32.00*Chord_Length, ZERO);
      x2 = Vector2D(32.00*Chord_Length, -32.00*Chord_Length);
      Bnd_Spline_West.Create_Spline_Line(x1, x2, 2);
      x1 = Vector2D(Chord_Length, ZERO);
      x2 = Vector2D(Chord_Length, -32.00*Chord_Length);
      Bnd_Spline_East.Create_Spline_Line(x1, x2, 2);
    } else if (iBlk == 1) {
      x1 = Vector2D(Chord_Length, -32.00*Chord_Length);
      x2 = Vector2D(ZERO, -32.00*Chord_Length);
      s_tmp1.Create_Spline_Line(x1, x2, 2);
      x1 = Vector2D(ZERO , ZERO);
      s_tmp2.Create_Spline_Circular_Arc(x1,
					32.00*Chord_Length,
					270.00,
					180.00,
					181);
      Bnd_Spline_North = Concatenate_Splines(s_tmp1, s_tmp2);
      s_tmp1.deallocate();
      s_tmp2.deallocate();
      Bnd_Spline_South.Create_Spline_NACA_Aerofoil(NACA_Aerofoil_Type_ptr,
						   Chord_Length,
						   -1,
						   501);
      x1 = Vector2D(Chord_Length, ZERO);
      x2 = Vector2D(Chord_Length, -32.00*Chord_Length);
      Bnd_Spline_West.Create_Spline_Line(x1, x2, 2);
      x1 = Vector2D(ZERO, ZERO);
      x2 = Vector2D(-32.00*Chord_Length, ZERO);
      Bnd_Spline_East.Create_Spline_Line(x1, x2, 2);
    } else if (iBlk == 2) {
      x1 = Vector2D(ZERO , ZERO);
      s_tmp1.Create_Spline_Circular_Arc(x1,
					32.00*Chord_Length,
					180.00,
					90.00,
					181);
      x1 = Vector2D(ZERO, 32.00*Chord_Length);
      x2 = Vector2D(Chord_Length, 32.00*Chord_Length);
      s_tmp2.Create_Spline_Line(x1, x2, 2);
      Bnd_Spline_North = Concatenate_Splines(s_tmp1, s_tmp2);
      s_tmp1.deallocate();
      s_tmp2.deallocate();
      Bnd_Spline_South.Create_Spline_NACA_Aerofoil(NACA_Aerofoil_Type_ptr,
						   Chord_Length,
						   1,
						   501);
      x1 = Vector2D(ZERO, ZERO);
      x2 = Vector2D(-32.00*Chord_Length, ZERO);
      Bnd_Spline_West.Create_Spline_Line(x1, x2, 2);
      x1 = Vector2D(Chord_Length, ZERO);
      x2 = Vector2D(Chord_Length, 32.00*Chord_Length);
      Bnd_Spline_East.Create_Spline_Line(x1, x2, 2);
    } else {
      x1 = Vector2D(Chord_Length, 32.00*Chord_Length);
      x2 = Vector2D(32.00*Chord_Length, 32.00*Chord_Length);
      Bnd_Spline_North.Create_Spline_Line(x1, x2, 2);
      x1 = Vector2D(Chord_Length, ZERO);
      x2 = Vector2D(32.00*Chord_Length, ZERO);
      Bnd_Spline_South.Create_Spline_Line(x1, x2, 2);
      x1 = Vector2D(Chord_Length, ZERO);
      x2 = Vector2D(Chord_Length, 32.00*Chord_Length);
      Bnd_Spline_West.Create_Spline_Line(x1, x2, 2);
      x1 = Vector2D(32.00*Chord_Length, ZERO);
      x2 = Vector2D(32.00*Chord_Length, 32.00*Chord_Length);
      Bnd_Spline_East.Create_Spline_Line(x1, x2, 2);
    } /* endif */

    /* Set the boundary condition types for each of the
       boundary splines. */

    if (iBlk == 0) {
      Bnd_Spline_North.setBCtype(BC_FIXED);
      Bnd_Spline_South.setBCtype(BC_NONE);
      Bnd_Spline_East.setBCtype(BC_NONE);
      Bnd_Spline_West.setBCtype(BC_FIXED);
    } else if (iBlk == 1) {
      Bnd_Spline_North.setBCtype(BC_FIXED);
      Bnd_Spline_South.setBCtype(BC_WALL_VISCOUS_HEATFLUX);//BC_REFLECTION);
      Bnd_Spline_South.makeSplineSolidBoundary(); // create first solid body      
      Bnd_Spline_East.setBCtype(BC_NONE);
      Bnd_Spline_West.setBCtype(BC_NONE);
    } else if (iBlk == 2) {
      Bnd_Spline_North.setBCtype(BC_FIXED);
      Bnd_Spline_South.setBCtype(BC_WALL_VISCOUS_HEATFLUX);//BC_REFLECTION);
      Bnd_Spline_South.makeSplineSolidBoundary(1); // make this spline part of the first solid body
      Bnd_Spline_East.setBCtype(BC_NONE);
      Bnd_Spline_West.setBCtype(BC_NONE);
    } else {
      Bnd_Spline_North.setBCtype(BC_FIXED);
      Bnd_Spline_South.setBCtype(BC_NONE);
      Bnd_Spline_East.setBCtype(BC_FIXED);
      Bnd_Spline_West.setBCtype(BC_NONE);
    } /* endif */

    /* Determine the number of cells for this block. */

    if (iBlk == 1 || iBlk == 2) {
      n_cells_i = 2*(((Number_of_Cells_Idir/3)+1)/2);
      n_cells_j = Number_of_Cells_Jdir;
    } else {
      n_cells_i = 2*(((Number_of_Cells_Idir/6)+1)/2);
      n_cells_j = Number_of_Cells_Jdir;
    } /* endif */
    n_cells_i = Number_of_Cells_Idir;
    n_cells_j = Number_of_Cells_Jdir;

    /* Assign values to the stretching function parameters
       and boundary grid line orthogonality parameters. */

    if (iBlk == 0) {
      Stretch_I = STRETCHING_FCN_MAX_CLUSTERING;
      Beta_I = 1.000075; 
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
      Beta_J = 1.00015;
      Tau_J = ZERO;
      Orthogonal_North = 0;
      Orthogonal_South = 0;
      Orthogonal_East = 0;
      Orthogonal_West = 0;
    } else if (iBlk == 1) {
      Stretch_I = STRETCHING_FCN_MINMAX_CLUSTERING;
      Beta_I = 1.01; 
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
      Beta_J = 1.00015;
      Tau_J = ZERO;
      Orthogonal_North = 1;
      Orthogonal_South = 1;
      Orthogonal_East = 1;
      Orthogonal_West = 1;
    } else if (iBlk == 2) {
      Stretch_I = STRETCHING_FCN_MINMAX_CLUSTERING;
      Beta_I = 1.01; 
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
      Beta_J = 1.00015;
      Tau_J = ZERO;
      Orthogonal_North = 1;
      Orthogonal_South = 1;
      Orthogonal_East = 1;
      Orthogonal_West = 1;
    } else {
      Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
      Beta_I = 1.000075; 
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
      Beta_J = 1.00015;
      Tau_J = ZERO;
      Orthogonal_North = 0;
      Orthogonal_South = 0;
      Orthogonal_East = 0;
      Orthogonal_West = 0;
    } /* endif */

    /* Create the 2D quadrilateral grid block. */

    Grid_ptr[iBlk][0].Create_Quad_Block_Without_Update(Bnd_Spline_North,
						       Bnd_Spline_South,
						       Bnd_Spline_East,
						       Bnd_Spline_West,
						       n_cells_i,
						       n_cells_j,
						       Number_of_Ghost_Cells,
						       Highest_Order_of_Reconstruction,
						       GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
						       Stretch_I,
						       Beta_I, 
						       Tau_I,
						       Stretch_J,
						       Beta_J,
						       Tau_J,
						       Orthogonal_North,
						       Orthogonal_South,
						       Orthogonal_East,
						       Orthogonal_West);
	
    /* Smooth the 2D quadrilateral grid block. */

    if (iBlk == 1 || iBlk == 2) 
      Grid_ptr[iBlk][0].Smooth_Quad_Block(min(250, 2*max(n_cells_i,n_cells_j)));

    /* Force the mesh to be symmetric. */

    //          if (iBlk == 2) {
    //             Copy_Quad_Block(Grid_ptr[2][0],
    //                             Grid_ptr[1][0]);
    //             Copy_Spline(Grid_ptr[2][0].BndNorthSpline,
    //                         Grid_ptr[1][0].BndNorthSpline);
    //             Copy_Spline(Grid_ptr[2][0].BndSouthSpline,
    //                         Grid_ptr[1][0].BndSouthSpline);
    //             Copy_Spline(Grid_ptr[2][0].BndEastSpline,
    //                         Grid_ptr[1][0].BndEastSpline);
    //             Copy_Spline(Grid_ptr[2][0].BndWestSpline,
    //                         Grid_ptr[1][0].BndWestSpline);
    //             Reflect_Quad_Block(Grid_ptr[2][0]);
    //          } else if (iBlk == 3) {
    //             Copy_Quad_Block(Grid_ptr[3][0],
    //                             Grid_ptr[0][0]);
    //             Copy_Spline(Grid_ptr[3][0].BndNorthSpline,
    //                         Grid_ptr[0][0].BndNorthSpline);
    //             Copy_Spline(Grid_ptr[3][0].BndSouthSpline,
    //                         Grid_ptr[0][0].BndSouthSpline);
    //             Copy_Spline(Grid_ptr[3][0].BndEastSpline,
    //                         Grid_ptr[0][0].BndEastSpline);
    //             Copy_Spline(Grid_ptr[3][0].BndWestSpline,
    //                         Grid_ptr[0][0].BndWestSpline);
    //             Reflect_Quad_Block(Grid_ptr[3][0]);
    //          } /* endif */

  } /* endfor */

}


/*!
 * Generates a multi-block grid for predicting free-jet 
 * flows associated with expansion through an orifice.  
 *                                                      
 * Usage: Grid_ptr = Grid_Free_Jet(Grid_ptr,            
 *                                 nblk_i,              
 *                                 nblk_j,              
 *                                 THREE,               
 *   		                   100,                 
 *  		                   50,                  
 *                                 2);                  
 *                                                      
 * This subroutine DOESN'T update the ghost cells or
 * the geometric properties of the grid cells.
 */
void Grid2D_Quad_MultiBlock_HO::Grid_Free_Jet_Without_Update(int &_Number_of_Blocks_Idir_,
							     int &_Number_of_Blocks_Jdir_,
							     const double &Radius,
							     const int Number_of_Cells_Idir,
							     const int Number_of_Cells_Jdir,
							     const int Number_of_Ghost_Cells,
							     const int Highest_Order_of_Reconstruction) {
  
  int iBlk, jBlk, n_cells_i, n_cells_j, Stretch_I, Stretch_J,
    Orthogonal_North, Orthogonal_South,
    Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D x1, x2;
  Spline2D_HO Bnd_Spline_North, Bnd_Spline_South,
              Bnd_Spline_East, Bnd_Spline_West;
  
  double x_orifice_up, x_orifice_down, x_up, x_down, r_max;
  
  /* Allocate memory for grid blocks.  There are two grid
     blocks for this mesh. */

  _Number_of_Blocks_Idir_ = 3;
  _Number_of_Blocks_Jdir_ = 2;
  allocate(_Number_of_Blocks_Idir_, _Number_of_Blocks_Jdir_);

  /* Create the mesh for each block representing
     the complete grid. */

  for ( jBlk = 0; jBlk <= Number_of_Blocks_Jdir-1; ++jBlk ) {
    for ( iBlk = 0; iBlk <= Number_of_Blocks_Idir-1; ++iBlk ) {

      x_orifice_up = -EIGHT*Radius;
      x_orifice_down = ZERO;
      x_up = -24*Radius;
      x_down = 60*Radius;
      r_max = 40*Radius;

      if (iBlk == 0 && jBlk == 0) {
	x1 = Vector2D(x_up, Radius);
	x2 = Vector2D(x_orifice_up, Radius);
	Bnd_Spline_North.Create_Spline_Line(x1, x2, 2);
	x1 = Vector2D(x_up, ZERO);
	x2 = Vector2D(x_orifice_up, ZERO);
	Bnd_Spline_South.Create_Spline_Line(x1, x2, 2);
	x1 = Vector2D(x_up, ZERO);
	x2 = Vector2D(x_up, Radius);
	Bnd_Spline_West.Create_Spline_Line(x1, x2, 2);
	x1 = Vector2D(x_orifice_up, ZERO);
	x2 = Vector2D(x_orifice_up, Radius);
	Bnd_Spline_East.Create_Spline_Line(x1, x2, 2);

	Bnd_Spline_North.setBCtype(BC_NONE);
	Bnd_Spline_South.setBCtype(BC_REFLECTION);
	Bnd_Spline_East.setBCtype(BC_NONE);
	Bnd_Spline_West.setBCtype(BC_FIXED);

	n_cells_i = Number_of_Cells_Idir/3;
	n_cells_j = Number_of_Cells_Jdir/2;

	Stretch_I = STRETCHING_FCN_MAX_CLUSTERING;
	Beta_I = 1.05; 
	Tau_I = ZERO;
	Stretch_J = STRETCHING_FCN_MINMAX_CLUSTERING;
	Beta_J = 1.5;
	Tau_J = ZERO;
	Orthogonal_North = 0;
	Orthogonal_South = 0;
	Orthogonal_East = 0;
	Orthogonal_West = 0;

	Grid_ptr[iBlk][jBlk].Create_Quad_Block_Without_Update(Bnd_Spline_North,
							      Bnd_Spline_South,
							      Bnd_Spline_East,
							      Bnd_Spline_West,
							      n_cells_i,
							      n_cells_j,
							      Number_of_Ghost_Cells,
							      Highest_Order_of_Reconstruction,
							      GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
							      Stretch_I,
							      Beta_I, 
							      Tau_I,
							      Stretch_J,
							      Beta_J,
							      Tau_J,
							      Orthogonal_North,
							      Orthogonal_South,
							      Orthogonal_East,
							      Orthogonal_West);
	      
	Bnd_Spline_North.deallocate();
	Bnd_Spline_South.deallocate();
	Bnd_Spline_East.deallocate();
	Bnd_Spline_West.deallocate();

      } else if (iBlk == 1 && jBlk == 0) {
	x1 = Vector2D(x_orifice_up, Radius);
	x2 = Vector2D(x_orifice_down, Radius);
	Bnd_Spline_North.Create_Spline_Line(x1, x2, 2);
	x1 = Vector2D(x_orifice_up, ZERO);
	x2 = Vector2D(x_orifice_down, ZERO);
	Bnd_Spline_South.Create_Spline_Line(x1, x2, 2);
	x1 = Vector2D(x_orifice_up, ZERO);
	x2 = Vector2D(x_orifice_up, Radius);
	Bnd_Spline_West.Create_Spline_Line(x1, x2, 2);
	x1 = Vector2D(x_orifice_down, ZERO);
	x2 = Vector2D(x_orifice_down, Radius);
	Bnd_Spline_East.Create_Spline_Line(x1, x2, 2);

	Bnd_Spline_North.setBCtype(BC_REFLECTION);
	Bnd_Spline_South.setBCtype(BC_REFLECTION);
	Bnd_Spline_East.setBCtype(BC_NONE);
	Bnd_Spline_West.setBCtype(BC_NONE);

	n_cells_i = Number_of_Cells_Idir/3;
	n_cells_j = Number_of_Cells_Jdir/2;

	Stretch_I = STRETCHING_FCN_MINMAX_CLUSTERING;
	Beta_I = 1.1; 
	Tau_I = ZERO;
	Stretch_J = STRETCHING_FCN_MINMAX_CLUSTERING;
	Beta_J = 1.5;
	Tau_J = ZERO;
	Orthogonal_North = 0;
	Orthogonal_South = 0;
	Orthogonal_East = 0;
	Orthogonal_West = 0;

	Grid_ptr[iBlk][jBlk].Create_Quad_Block_Without_Update(Bnd_Spline_North,
							      Bnd_Spline_South,
							      Bnd_Spline_East,
							      Bnd_Spline_West,
							      n_cells_i,
							      n_cells_j,
							      Number_of_Ghost_Cells,
							      Highest_Order_of_Reconstruction,
							      GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
							      Stretch_I,
							      Beta_I, 
							      Tau_I,
							      Stretch_J,
							      Beta_J,
							      Tau_J,
							      Orthogonal_North,
							      Orthogonal_South,
							      Orthogonal_East,
							      Orthogonal_West);

	Bnd_Spline_North.deallocate();
	Bnd_Spline_South.deallocate();
	Bnd_Spline_East.deallocate();
	Bnd_Spline_West.deallocate();

      } else if (iBlk == 2 && jBlk == 0) {
	x1 = Vector2D(x_orifice_down, Radius);
	x2 = Vector2D(x_down, Radius);
	Bnd_Spline_North.Create_Spline_Line(x1, x2, 2);
	x1 = Vector2D(x_orifice_down, ZERO);
	x2 = Vector2D(x_down, ZERO);
	Bnd_Spline_South.Create_Spline_Line(x1, x2, 2);
	x1 = Vector2D(x_orifice_down, ZERO);
	x2 = Vector2D(x_orifice_down, Radius);
	Bnd_Spline_West.Create_Spline_Line(x1, x2, 2);
	x1 = Vector2D(x_down, ZERO);
	x2 = Vector2D(x_down, Radius);
	Bnd_Spline_East.Create_Spline_Line(x1, x2, 2);

	Bnd_Spline_North.setBCtype(BC_NONE);
	Bnd_Spline_South.setBCtype(BC_REFLECTION);
	Bnd_Spline_East.setBCtype(BC_FIXED);
	Bnd_Spline_West.setBCtype(BC_NONE);

	n_cells_i = Number_of_Cells_Idir/3;
	n_cells_j = Number_of_Cells_Jdir/2;

	Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_I = 1.01; 
	Tau_I = ZERO;
	Stretch_J = STRETCHING_FCN_MINMAX_CLUSTERING;
	Beta_J = 1.5;
	Tau_J = ZERO;
	Orthogonal_North = 0;
	Orthogonal_South = 0;
	Orthogonal_East = 0;
	Orthogonal_West = 0;

	Grid_ptr[iBlk][jBlk].Create_Quad_Block_Without_Update(Bnd_Spline_North,
							      Bnd_Spline_South,
							      Bnd_Spline_East,
							      Bnd_Spline_West,
							      n_cells_i,
							      n_cells_j,
							      Number_of_Ghost_Cells,
							      Highest_Order_of_Reconstruction,
							      GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
							      Stretch_I,
							      Beta_I, 
							      Tau_I,
							      Stretch_J,
							      Beta_J,
							      Tau_J,
							      Orthogonal_North,
							      Orthogonal_South,
							      Orthogonal_East,
							      Orthogonal_West);
	      
	Bnd_Spline_North.deallocate();
	Bnd_Spline_South.deallocate();
	Bnd_Spline_East.deallocate();
	Bnd_Spline_West.deallocate();

      } else if (iBlk == 0 && jBlk == 1) {
	x1 = Vector2D(x_up, r_max);
	x2 = Vector2D(x_orifice_up, r_max);
	Bnd_Spline_North.Create_Spline_Line(x1, x2, 2);
	x1 = Vector2D(x_up, Radius);
	x2 = Vector2D(x_orifice_up, Radius);
	Bnd_Spline_South.Create_Spline_Line(x1, x2, 2);
	x1 = Vector2D(x_up, Radius);
	x2 = Vector2D(x_up, r_max);
	Bnd_Spline_West.Create_Spline_Line(x1, x2, 2);
	x1 = Vector2D(x_orifice_up, Radius);
	x2 = Vector2D(x_orifice_up, r_max);
	Bnd_Spline_East.Create_Spline_Line(x1, x2, 2);

	Bnd_Spline_North.setBCtype(BC_FIXED);
	Bnd_Spline_South.setBCtype(BC_NONE);
	Bnd_Spline_East.setBCtype(BC_REFLECTION);
	Bnd_Spline_West.setBCtype(BC_FIXED);

	n_cells_i = Number_of_Cells_Idir/3;
	n_cells_j = Number_of_Cells_Jdir/2;

	Stretch_I = STRETCHING_FCN_MAX_CLUSTERING;
	Beta_I = 1.05; 
	Tau_I = ZERO;
	Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_J = 1.004;
	Tau_J = ZERO;
	Orthogonal_North = 0;
	Orthogonal_South = 0;
	Orthogonal_East = 0;
	Orthogonal_West = 0;

	Grid_ptr[iBlk][jBlk].Create_Quad_Block_Without_Update(Bnd_Spline_North,
							      Bnd_Spline_South,
							      Bnd_Spline_East,
							      Bnd_Spline_West,
							      n_cells_i,
							      n_cells_j,
							      Number_of_Ghost_Cells,
							      Highest_Order_of_Reconstruction,
							      GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
							      Stretch_I,
							      Beta_I, 
							      Tau_I,
							      Stretch_J,
							      Beta_J,
							      Tau_J,
							      Orthogonal_North,
							      Orthogonal_South,
							      Orthogonal_East,
							      Orthogonal_West);
	      
	Bnd_Spline_North.deallocate();
	Bnd_Spline_South.deallocate();
	Bnd_Spline_East.deallocate();
	Bnd_Spline_West.deallocate();

      } else if (iBlk == 1 && jBlk == 1) {

      } else if (iBlk == 2 && jBlk == 1) {
	x1 = Vector2D(x_orifice_down, r_max);
	x2 = Vector2D(x_down, r_max);
	Bnd_Spline_North.Create_Spline_Line(x1, x2, 2);
	x1 = Vector2D(x_orifice_down, Radius);
	x2 = Vector2D(x_down, Radius);
	Bnd_Spline_South.Create_Spline_Line(x1, x2, 2);
	x1 = Vector2D(x_orifice_down, Radius);
	x2 = Vector2D(x_orifice_down, r_max);
	Bnd_Spline_West.Create_Spline_Line(x1, x2, 2);
	x1 = Vector2D(x_down, Radius);
	x2 = Vector2D(x_down, r_max);
	Bnd_Spline_East.Create_Spline_Line(x1, x2, 2);

	Bnd_Spline_North.setBCtype(BC_FIXED);
	Bnd_Spline_South.setBCtype(BC_NONE);
	Bnd_Spline_East.setBCtype(BC_FIXED);
	Bnd_Spline_West.setBCtype(BC_REFLECTION);

	n_cells_i = Number_of_Cells_Idir/3;
	n_cells_j = Number_of_Cells_Jdir/2;

	Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_I = 1.01; 
	Tau_I = ZERO;
	Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_J = 1.004;
	Tau_J = ZERO;
	Orthogonal_North = 0;
	Orthogonal_South = 0;
	Orthogonal_East = 0;
	Orthogonal_West = 0;

	Grid_ptr[iBlk][jBlk].Create_Quad_Block_Without_Update(Bnd_Spline_North,
							      Bnd_Spline_South,
							      Bnd_Spline_East,
							      Bnd_Spline_West,
							      n_cells_i,
							      n_cells_j,
							      Number_of_Ghost_Cells,
							      Highest_Order_of_Reconstruction,
							      GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
							      Stretch_I,
							      Beta_I, 
							      Tau_I,
							      Stretch_J,
							      Beta_J,
							      Tau_J,
							      Orthogonal_North,
							      Orthogonal_South,
							      Orthogonal_East,
							      Orthogonal_West);
	      
	Bnd_Spline_North.deallocate();
	Bnd_Spline_South.deallocate();
	Bnd_Spline_East.deallocate();
	Bnd_Spline_West.deallocate();
      } /* endif */

    } /* endfor */
  } /* endfor */
}


/*!
 * Generates a quadilateral mesh with clustering consisting of two    
 * grid blocks for predicting a wedge flow for looking at shock wave  
 * interactions at a wedge.                                           
 *                                                                    
 * Usage: Grid_ptr = Grid_Wedge(Grid_ptr,                             
 *                              nblk_i,                               
 *                              nblk_j,                               
 *                              45.0,                                 
 *                              0.50,                                 
 *         		        100,                                  
 *         		        100,                                  
 *                              2);                                   
 * This subroutine DOESN'T update the ghost cells or
 * the geometric properties of the grid cells.
 */
void Grid2D_Quad_MultiBlock_HO::Grid_Wedge_Without_Update(int &_Number_of_Blocks_Idir_,
							  int &_Number_of_Blocks_Jdir_,
							  const double &Wedge_Angle,
							  const double &Wedge_Length,
							  const int &Wedge_BC_Type,
							  const int &Stretching_Flag,
							  const double &Stretching_Factor_Idir,
							  const double &Stretching_Factor_Jdir,
							  const int Number_of_Cells_Idir,
							  const int Number_of_Cells_Jdir,
							  const int Number_of_Ghost_Cells,
							  const int Highest_Order_of_Reconstruction) {
  
  int Stretch_I, Stretch_J,
    Orthogonal_North, Orthogonal_South,
    Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D xc_NW, xc_NE, xc_SE, xc_SW;
  Spline2D_HO Bnd_Spline_North, Bnd_Spline_South,
              Bnd_Spline_East, Bnd_Spline_West;

  // Allocate memory for grid blocks.  There are two grid blocks for 
  // this mesh.
  _Number_of_Blocks_Idir_ = 2;
  _Number_of_Blocks_Jdir_ = 1;
  allocate(_Number_of_Blocks_Idir_, _Number_of_Blocks_Jdir_);

  // Create the mesh for each block representing the complete grid.
  for (int iBlk = 0; iBlk < Number_of_Blocks_Idir; iBlk++) {

    if (iBlk == 0) {
      xc_NW = Vector2D(-Wedge_Length,Wedge_Length);
      xc_NE = Vector2D(ZERO,Wedge_Length);
      xc_SE = Vector2D(ZERO,ZERO);
      xc_SW = Vector2D(-Wedge_Length,ZERO);
    } else if (iBlk == 1) {
      xc_NW = Vector2D(ZERO,Wedge_Length);
      xc_NE = Vector2D(Wedge_Length*cos(PI*Wedge_Angle/180.0),Wedge_Length);
      xc_SE = Vector2D(Wedge_Length*cos(PI*Wedge_Angle/180.0),Wedge_Length*sin(PI*Wedge_Angle/180.0));
      xc_SW = Vector2D(ZERO,ZERO);
    }
    // Create the splines defining the north, south, east, and west 
    // boundaries of the rectangular box.
    Bnd_Spline_North.Create_Spline_Line(xc_NW, xc_NE, 2);
    Bnd_Spline_South.Create_Spline_Line(xc_SW, xc_SE, 2);
    Bnd_Spline_East.Create_Spline_Line(xc_SE, xc_NE, 2);
    Bnd_Spline_West.Create_Spline_Line(xc_SW, xc_NW, 2);
    // Set the boundary condition types for each of the boundary splines.
    Bnd_Spline_North.setBCtype(BC_REFLECTION);
    Bnd_Spline_South.setBCtype(Wedge_BC_Type);
    if (iBlk == 0) {
      Bnd_Spline_East.setBCtype(BC_NONE);
      Bnd_Spline_West.setBCtype(BC_FIXED);
    } else if (iBlk == 1) {
      Bnd_Spline_East.setBCtype(BC_CONSTANT_EXTRAPOLATION);
      Bnd_Spline_West.setBCtype(BC_NONE);
    }
    // Assign values to the stretching function parameters and
    // boundary grid line orthogonality parameters.
    if (Stretching_Flag) {
      if (iBlk == 0) Stretch_I = STRETCHING_FCN_MAX_CLUSTERING;
      else Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
      Beta_I = Stretching_Factor_Idir;
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
      Beta_J = Stretching_Factor_Jdir;
    } else {
      Stretch_I = STRETCHING_FCN_LINEAR;
      Beta_I = ONE;
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_LINEAR;
      Beta_J = ONE;
      Tau_J = ZERO;
    }
    Orthogonal_North = 0;
    Orthogonal_South = 0;
    Orthogonal_East = 0;
    Orthogonal_West = 0;
    // Create the 2D quadrilateral grid block.
    Grid_ptr[iBlk][0].Create_Quad_Block_Without_Update(Bnd_Spline_North,
						       Bnd_Spline_South,
						       Bnd_Spline_East,
						       Bnd_Spline_West,
						       Number_of_Cells_Idir,
						       Number_of_Cells_Jdir,
						       Number_of_Ghost_Cells,
						       Highest_Order_of_Reconstruction,
						       GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
						       Stretch_I,
						       Beta_I, 
						       Tau_I,
						       Stretch_J,
						       Beta_J,
						       Tau_J,
						       Orthogonal_North,
						       Orthogonal_South,
						       Orthogonal_East,
						       Orthogonal_West);
    Smooth_Quad_Block(Grid_ptr[iBlk][0], min(250, 2*max(Number_of_Cells_Idir,Number_of_Cells_Jdir)));
  }

}

/*!
 * Generates a single block quadilateral mesh with clustering for     
 * predicting supersonic flow around a cirucular cylinder blunt body. 
 *                                                                    
 * Usage: Grid_ptr = Grid_Unsteady_Blunt_Body(Grid_ptr,               
 *                                            nblk_i,                 
 *                                            nblk_j,                 
 *                                            ONE,                    
 *                                            FOUR,                   
 *   	                                      150,                    
 *  	                                      50,                     
 *                                            2);                    
 * This subroutine DOESN'T update the ghost cells or
 * the geometric properties of the grid cells.
 */
void Grid2D_Quad_MultiBlock_HO::Grid_Unsteady_Blunt_Body_Without_Update(int &_Number_of_Blocks_Idir_,
									int &_Number_of_Blocks_Jdir_,
									const double &Radius,
									const double &Mach_Number,
									const int Number_of_Cells_Idir,
									const int Number_of_Cells_Jdir,
									const int Number_of_Ghost_Cells,
									const int Highest_Order_of_Reconstruction) {
  
  int Stretch_I, Stretch_J,
    Orthogonal_North, Orthogonal_South,
    Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D xc_NW, xc_NE, xc_SE, xc_SW;
  Spline2D_HO Bnd_Spline_North, Bnd_Spline_South,
              Bnd_Spline_East,  Bnd_Spline_West, Bow_Spline;
  
  // Allocate memory for grid block.
  _Number_of_Blocks_Idir_ = 1;
  _Number_of_Blocks_Jdir_ = 3;
  allocate(_Number_of_Blocks_Idir_, _Number_of_Blocks_Jdir_);
  
  // Create the mesh for each block representing the complete grid.
  for (int jBlk = 0; jBlk < Number_of_Blocks_Jdir; jBlk++) {
    // Create the splines defining the north, south, east, and west 
    // boundaries of the rectangular box.
    if (jBlk == 2) {
      Bow_Spline.Create_Spline_Bow_Shock(Radius,
					 Mach_Number,
					 1,
					 181);
      xc_SW = Bow_Spline.Xp[0] - Vector2D(TWO*Radius,ZERO);
      xc_SE = Vector2D(xc_SW.x,Bow_Spline.Xp[180].y);
      xc_NW = Bow_Spline.Xp[0] - Vector2D(FOUR*Radius,ZERO);
      xc_NE = Vector2D(xc_NW.x,Bow_Spline.Xp[180].y);
      Bnd_Spline_North.Create_Spline_Line(xc_NW,xc_NE,2);
      Bnd_Spline_South.Create_Spline_Line(xc_SW,xc_SE,2);
      Bnd_Spline_East.Create_Spline_Line(xc_SE,xc_NE,2);
      Bnd_Spline_West.Create_Spline_Line(xc_SW,xc_NW,2);
    } else if (jBlk == 1) {
      Bnd_Spline_South.Create_Spline_Bow_Shock(Radius,
					       Mach_Number,
					       1,
					       181);
      xc_SW = Bnd_Spline_South.Xp[0];
      xc_SE = Bnd_Spline_South.Xp[180];
      xc_NW = Bnd_Spline_South.Xp[0] - Vector2D(TWO*Radius,ZERO);
      xc_NE = Vector2D(xc_NW.x,Bnd_Spline_South.Xp[180].y);
      Bnd_Spline_North.Create_Spline_Line(xc_NW,xc_NE,2);
      Bnd_Spline_East.Create_Spline_Line(xc_SE,xc_NE,2);
      Bnd_Spline_West.Create_Spline_Line(xc_SW,xc_NW,2);
    } else if (jBlk == 0) {
      Bnd_Spline_North.Create_Spline_Bow_Shock(Radius,
					       Mach_Number,
					       1,
					       181);
      Bnd_Spline_South.Create_Spline_Circular_Arc(Vector2D(ZERO,ZERO),
						  Radius,
						  180.00,
						  90.00,
						  181);
      xc_SW = Bnd_Spline_South.Xp[0];
      xc_NW = Bnd_Spline_North.Xp[0];
      Bnd_Spline_West.Create_Spline_Line(xc_SW, xc_NW, 2);
      xc_SE = Bnd_Spline_South.Xp[Bnd_Spline_South.np-1];
      xc_NE = Bnd_Spline_North.Xp[Bnd_Spline_North.np-1];
      Bnd_Spline_East.Create_Spline_Line(xc_SE, xc_NE, 2);
    }
    // Set the boundary condition types for each of the boundary splines.
    if (jBlk == 2) {
      Bnd_Spline_North.setBCtype(BC_FIXED);
      Bnd_Spline_South.setBCtype(BC_NONE);
      Bnd_Spline_East.setBCtype(BC_CONSTANT_EXTRAPOLATION);
      Bnd_Spline_West.setBCtype(BC_REFLECTION);
    } else if (jBlk == 1) {
      Bnd_Spline_North.setBCtype(BC_NONE);
      Bnd_Spline_South.setBCtype(BC_NONE);
      Bnd_Spline_East.setBCtype(BC_CONSTANT_EXTRAPOLATION);
      Bnd_Spline_West.setBCtype(BC_REFLECTION);
    } else if (jBlk == 0) {
      Bnd_Spline_North.setBCtype(BC_NONE);
      Bnd_Spline_South.setBCtype(BC_REFLECTION);
      Bnd_Spline_East.setBCtype(BC_CONSTANT_EXTRAPOLATION);
      Bnd_Spline_West.setBCtype(BC_REFLECTION);
    }
    // Assign values to the stretching function parameters and
    // boundary grid line orthogonality parameters.
    if (jBlk == 2 || jBlk == 1) {
      Stretch_I = STRETCHING_FCN_LINEAR;
      Beta_I = ONE; 
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_LINEAR;
      Beta_J = ONE;
      Tau_J = ZERO;
      Orthogonal_North = 0;
      Orthogonal_South = 0;
      Orthogonal_East = 0;
      Orthogonal_West = 0;
    } else if (jBlk == 0) {
      Stretch_I = STRETCHING_FCN_LINEAR;
      Beta_I = ZERO; 
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
      Beta_J = 1.05;
      Tau_J = ZERO;
      Orthogonal_North = 0;
      Orthogonal_South = 1;
      Orthogonal_East = 1;
      Orthogonal_West = 1;
    }
    // Create the 2D quadrilateral grid block.
    Grid_ptr[0][jBlk].Create_Quad_Block_Without_Update(Bnd_Spline_North,
						       Bnd_Spline_South,
						       Bnd_Spline_East,
						       Bnd_Spline_West,
						       Number_of_Cells_Idir,
						       Number_of_Cells_Jdir,
						       Number_of_Ghost_Cells,
						       Highest_Order_of_Reconstruction,
						       GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
						       Stretch_I,
						       Beta_I, 
						       Tau_I,
						       Stretch_J,
						       Beta_J,
						       Tau_J,
						       Orthogonal_North,
						       Orthogonal_South,
						       Orthogonal_East,
						       Orthogonal_West);
    // Smooth the 2D quadrilateral grid block.
    Grid_ptr[0][jBlk].Smooth_Quad_Block(min(250, 2*max(Number_of_Cells_Idir,Number_of_Cells_Jdir)));
  }
  
}


void Grid2D_Quad_MultiBlock_HO::Determine_Coordinates_Ringleb_Flow(const double & Streamline, const double & Isotachline,
								   double & xLoc, double & yLoc){
  
  const double g(1.40);
  double c, J, rho;

  const double Isotachline2(Isotachline*Isotachline);
  const double Streamline2(Streamline*Streamline);

  c = sqrt(ONE - ((g-ONE)/TWO)*Isotachline2);
  rho = pow(c,TWO/(g-ONE));
  J = ONE/c + ONE/(THREE*pow(c,THREE)) + ONE/(FIVE*pow(c,FIVE)) - HALF*log((ONE+c)/(ONE-c));

  xLoc = (HALF/rho)*(TWO/(Streamline2) - ONE/(Isotachline2)) - HALF*J;
  yLoc = (ONE/(Streamline*rho*Isotachline))*sqrt(ONE - (Isotachline2)/(Streamline2));
  
}

/*!
 * Generates a uniform 2D mesh for Ringleb's flow.                    
 *                                                                    
 * Usage: Grid_ptr = Grid_Ringleb_Flow(Grid_ptr,                      
 *                                     nblk_i,                        
 *                                     nblk_j,                        
 *   	                               ncells_i,                      
 *  	                               ncells_j,                      
 *                                     Nghost_Cells);                 
 * This subroutine DOESN'T update the ghost cells or
 * the geometric properties of the grid cells.
 */
void Grid2D_Quad_MultiBlock_HO::Grid_Ringleb_Flow_Without_Update(int &_Number_of_Blocks_Idir_,
								 int &_Number_of_Blocks_Jdir_,
								 const double &Inner_Streamline_Number,
								 const double &Outer_Streamline_Number,
								 const double &Isotach_Line,
								 const int Number_of_Cells_Idir,
								 const int Number_of_Cells_Jdir,
								 const int Number_of_Ghost_Cells,
								 const int Highest_Order_of_Reconstruction) {

  assert(Inner_Streamline_Number > Outer_Streamline_Number);
  assert(Inner_Streamline_Number < 5.0/3.0);
  assert(Outer_Streamline_Number > Isotach_Line); 

  int nk, nq,
    Stretch_I, Stretch_J,
    Orthogonal_North, Orthogonal_South,
    Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D xc_NW, xc_NE, xc_SE, xc_SW;
  Spline2D_HO Bnd_Spline_North, Bnd_Spline_South,
    Bnd_Spline_East, Bnd_Spline_West;

  double k, k1, k2;
  double q, q1, q2;
  double delk, delq;
  int i,j;


  // Allocate memory for grid block.
  _Number_of_Blocks_Idir_ = 1;
  _Number_of_Blocks_Jdir_ = 1;
  allocate(_Number_of_Blocks_Idir_, _Number_of_Blocks_Jdir_);

  // Create the mesh for each block representing the complete grid.
  
  // Set number of points in each direction
  nk = min(32, 3*Number_of_Cells_Idir); // 3 points for each cell
  nq = min(50, 3*Number_of_Cells_Jdir); // 3 points for each cell

  // West streamline
  k  = Inner_Streamline_Number;
  q1 = Inner_Streamline_Number;
  q2 = Isotach_Line;
  delq = q2 - q1;

  // Allocate memory
  Bnd_Spline_West.allocate(nq);	 Bnd_Spline_West.settype(SPLINE2D_QUINTIC);

  for (j = 0; j < nq; ++j ){
    // Determine isotach line
    q = q1 + delq*pow((0.5 - cos((j*PI)/(nq - 1))/2.0),1.3);
    
    // Determine (xLoc,yLoc) of the control point
    Determine_Coordinates_Ringleb_Flow(k,q,
				       Bnd_Spline_West.Xp[j].x,
				       Bnd_Spline_West.Xp[j].y);
    if (j == 0 || j == nq-1) {
      Bnd_Spline_West.tp[j] = SPLINE2D_POINT_SHARP_CORNER;
    } else {
      Bnd_Spline_West.tp[j] = SPLINE2D_POINT_NORMAL;
    }

  }

  // East streamline
  k  = Outer_Streamline_Number;
  q1 = Outer_Streamline_Number; 
  q2 = Isotach_Line;
  delq = q2 - q1;
  
  // Allocate memory
  Bnd_Spline_East.allocate(nq);	 Bnd_Spline_East.settype(SPLINE2D_QUINTIC);

  for (j = 0; j < nq; ++j ){
    // Determine isotach line
    q = q1 + delq*pow((0.5 - cos((j*PI)/(nq - 1))/2.0),1.3);

    // Determine (xLoc,yLoc) of the control point
    Determine_Coordinates_Ringleb_Flow(k,q,
				       Bnd_Spline_East.Xp[j].x,
				       Bnd_Spline_East.Xp[j].y);

    if (j == 0 || j == nq-1) {
      Bnd_Spline_East.tp[j] = SPLINE2D_POINT_SHARP_CORNER;
    } else {
      Bnd_Spline_East.tp[j] = SPLINE2D_POINT_NORMAL;
    }
  }


  // South isotach 
  // Line between the end-points of West and East streamline.
  Bnd_Spline_South.Create_Spline_Line(Bnd_Spline_West.Xp[0],
				      Bnd_Spline_East.Xp[0],
				      2);

  // North isotach
  q  = Isotach_Line;
  k1 = Outer_Streamline_Number; 
  k2 = Inner_Streamline_Number; 
  delk = k2 - k1;

  // Allocate memory
  Bnd_Spline_North.allocate(nk); Bnd_Spline_North.settype(SPLINE2D_QUINTIC);
  
  for (i = 0; i < nk; ++i ){
    // Determine streamline
    k = k1 + delk*(1. + sin((i/(nk - 1.) - 1.)*PI/2));

    // Determine (xLoc,yLoc) of the control point
    Determine_Coordinates_Ringleb_Flow(k,q,
				       Bnd_Spline_North.Xp[nk-1-i].x,
				       Bnd_Spline_North.Xp[nk-1-i].y);
   
    if (i == 0 || i == nk-1) {
      Bnd_Spline_North.tp[nk-1-i] = SPLINE2D_POINT_SHARP_CORNER;
    } else {
      Bnd_Spline_North.tp[nk-1-i] = SPLINE2D_POINT_NORMAL;
    }
  }

  Bnd_Spline_North.pathlength();
  Bnd_Spline_South.pathlength();
  Bnd_Spline_East.pathlength();
  Bnd_Spline_West.pathlength();

  // Set the boundary condition types for each of the boundary 
  // splines.
  Bnd_Spline_North.setBCtype(BC_RINGLEB_FLOW);
  Bnd_Spline_South.setBCtype(BC_RINGLEB_FLOW);
  Bnd_Spline_East.setBCtype(BC_RINGLEB_FLOW);
  Bnd_Spline_West.setBCtype(BC_RINGLEB_FLOW);

  // Assign values to the stretching function parameters and
  // boundary grid line orthogonality parameters.
  Stretch_I = STRETCHING_FCN_LINEAR;
  Beta_I = ZERO; 
  Tau_I = ZERO;
  Stretch_J = STRETCHING_FCN_LINEAR;
  Beta_J = ONE;
  Tau_J = THREE;
  Orthogonal_North = 1;
  Orthogonal_South = 1;
  Orthogonal_East = 1;
  Orthogonal_West = 1;

  // Create the 2D quadrilateral grid block representing the mesh.
  Grid_ptr[0][0].Create_Quad_Block_Without_Update(Bnd_Spline_North,
						  Bnd_Spline_South,
						  Bnd_Spline_East,
						  Bnd_Spline_West,
						  Number_of_Cells_Idir,
						  Number_of_Cells_Jdir,
						  Number_of_Ghost_Cells,
						  Highest_Order_of_Reconstruction,
						  GRID2D_QUAD_BLOCK_INIT_PROCEDURE_EAST_WEST,
						  Stretch_I,
						  Beta_I,
						  Tau_I,
						  Stretch_J,
						  Beta_J,
						  Tau_J,
						  Orthogonal_North,
						  Orthogonal_South,
						  Orthogonal_East,
						  Orthogonal_West);

  Grid_ptr[0][0].Smooth_Quad_Block(min(250,2*max(Number_of_Cells_Idir,Number_of_Cells_Jdir)));

}

/*!
 * Generates a single block quadilateral mesh with clustering for     
 * predicting supersonic flow around a cirucular cylinder blunt body. 
 *                                                                    
 * Usage: Grid_ptr = Grid_ump_Channel_Flow(Grid_ptr,                  
 *                                         nblk_i,                    
 *                                         nblk_j,                    
 *                                         ONE,                       
 *   	                                   150,                       
 *  	                                   50,                        
 *                                         2);                        
 * This subroutine DOESN'T update the ghost cells or
 * the geometric properties of the grid cells.
 */
void Grid2D_Quad_MultiBlock_HO::Grid_Bump_Channel_Flow_Without_Update(int &_Number_of_Blocks_Idir_,
								      int &_Number_of_Blocks_Jdir_,
								      const int Smooth_Bump,
								      const int Number_of_Cells_Idir,
								      const int Number_of_Cells_Jdir,
								      const int Number_of_Ghost_Cells,
								      const int Highest_Order_of_Reconstruction) {
  
  int Stretch_I, Stretch_J,
      Orthogonal_North, Orthogonal_South,
      Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D x1, x2;
  Spline2D_HO Bnd_Spline_North, Bnd_Spline_South,
           Bnd_Spline_East, Bnd_Spline_West;
  double R, Theta;

  // Allocate memory for grid block.
  _Number_of_Blocks_Idir_ = 4;
  _Number_of_Blocks_Jdir_ = 2;
  allocate(_Number_of_Blocks_Idir_, _Number_of_Blocks_Jdir_);

  // Create the mesh for each block representing the complete grid.
  for (int jBlk = 0; jBlk < Number_of_Blocks_Jdir; jBlk++) {
    for (int iBlk = 0; iBlk < Number_of_Blocks_Idir; iBlk++) {

      if (iBlk == 0 && jBlk == 0) {
	// Create the splines defining the north, south, east, and west 
	// boundaries of the grid.
	x1 = Vector2D(-1.0,1.0);
	x2 = Vector2D( 0.0,1.0);
	Bnd_Spline_North.Create_Spline_Line(x1,x2,2);
	x1 = Vector2D(-1.0,0.0);
	x2 = Vector2D( 0.0,0.0);
	Bnd_Spline_South.Create_Spline_Line(x1,x2,2);
	x1 = Vector2D( 0.0,0.0);
	x2 = Vector2D( 0.0,1.0);
	Bnd_Spline_East.Create_Spline_Line(x1,x2,2);
	x1 = Vector2D(-1.0,0.0);
	x2 = Vector2D(-1.0,1.0);
	Bnd_Spline_West.Create_Spline_Line(x1,x2,2);
      } else if (iBlk == 1 && jBlk == 0) {
	// Create the splines defining the north, south, east, and west 
	// boundaries of the grid.
	x1 = Vector2D( 0.0,1.0);
	x2 = Vector2D( 1.0,1.0);
	Bnd_Spline_North.Create_Spline_Line(x1,x2,2);
	if (!Smooth_Bump) {
	  // Non-smooth circular bump.
	  R = (0.25 + 0.042*0.042)/(TWO*0.042);
	  Theta = acos((R - 0.042)/R);
	  Bnd_Spline_South.Create_Spline_Circular_Arc(Vector2D_ZERO,
						      R,
						      ZERO,
						      -TWO*Theta*180.0/PI,
						      31);
	  Bnd_Spline_South.Rotate_Spline(HALF*PI+Theta);
	  Bnd_Spline_South.Translate_Spline(Vector2D(HALF,0.042-R));
	} else {
	  // Smooth sin^2 bump.
	  Bnd_Spline_South.allocate(31);
	  Bnd_Spline_South.settype(SPLINE2D_QUINTIC);
	  for (int i = 0; i < 31; i++) {
	    Bnd_Spline_South.Xp[i].x = double(i)/30.0;
	    Bnd_Spline_South.Xp[i].y = 0.042*sqr(sin(PI*double(i)/30.0));
	    Bnd_Spline_South.bc[i] = BC_REFLECTION;
	    Bnd_Spline_South.tp[i] = SPLINE2D_POINT_NORMAL;
	  }
	  Bnd_Spline_South.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
	  Bnd_Spline_South.tp[30] = SPLINE2D_POINT_SHARP_CORNER;
	  Bnd_Spline_South.pathlength();
	}
	x1 = Vector2D( 1.0,0.0);
	x2 = Vector2D( 1.0,1.0);
	Bnd_Spline_East.Create_Spline_Line(x1,x2,2);
	x1 = Vector2D( 0.0,0.0);
	x2 = Vector2D( 0.0,1.0);
	Bnd_Spline_West.Create_Spline_Line(x1,x2,2);
      } else if (iBlk == 2 && jBlk == 0) {
	// Create the splines defining the north, south, east, and west 
	// boundaries of the grid.
	x1 = Vector2D( 1.00,1.00);
	x2 = Vector2D( 2.75,1.00);
	Bnd_Spline_North.Create_Spline_Line(x1,x2,2);
	x1 = Vector2D( 1.00,0.00);
	x2 = Vector2D( 2.75,0.00);
	Bnd_Spline_South.Create_Spline_Line(x1,x2,2);
	x1 = Vector2D( 2.75,0.00);
	x2 = Vector2D( 2.75,1.00);
	Bnd_Spline_East.Create_Spline_Line(x1,x2,2);
	x1 = Vector2D( 1.00,0.00);
	x2 = Vector2D( 1.00,1.00);
	Bnd_Spline_West.Create_Spline_Line(x1,x2,2);
      } else if (iBlk == 3 && jBlk == 0) {
	// Create the splines defining the north, south, east, and west 
	// boundaries of the grid.
	x1 = Vector2D( 2.75,1.00);
	x2 = Vector2D( 4.50,1.00);
	Bnd_Spline_North.Create_Spline_Line(x1,x2,2);
	x1 = Vector2D( 2.75,0.00);
	x2 = Vector2D( 4.50,0.00);
	Bnd_Spline_South.Create_Spline_Line(x1,x2,2);
	x1 = Vector2D( 4.50,0.00);
	x2 = Vector2D( 4.50,1.00);
	Bnd_Spline_East.Create_Spline_Line(x1,x2,2);
	x1 = Vector2D( 2.75,0.00);
	x2 = Vector2D( 2.75,1.00);
	Bnd_Spline_West.Create_Spline_Line(x1,x2,2);
      } else if (iBlk == 0 && jBlk == 1) {
	// Create the splines defining the north, south, east, and west 
	// boundaries of the grid.
	x1 = Vector2D(-1.00,2.00);
	x2 = Vector2D( 0.00,2.00);
	Bnd_Spline_North.Create_Spline_Line(x1,x2,2);
	x1 = Vector2D(-1.00,1.00);
	x2 = Vector2D( 0.00,1.00);
	Bnd_Spline_South.Create_Spline_Line(x1,x2,2);
	x1 = Vector2D( 0.00,1.00);
	x2 = Vector2D( 0.00,2.00);
	Bnd_Spline_East.Create_Spline_Line(x1,x2,2);
	x1 = Vector2D(-1.00,1.00);
	x2 = Vector2D(-1.00,2.00);
	Bnd_Spline_West.Create_Spline_Line(x1,x2,2);
      } else if (iBlk == 1 && jBlk == 1) {
	// Create the splines defining the north, south, east, and west 
	// boundaries of the grid.
	x1 = Vector2D( 0.00,2.00);
	x2 = Vector2D( 1.00,2.00);
	Bnd_Spline_North.Create_Spline_Line(x1,x2,2);
	x1 = Vector2D( 0.00,1.00);
	x2 = Vector2D( 1.00,1.00);
	Bnd_Spline_South.Create_Spline_Line(x1,x2,2);
	x1 = Vector2D( 1.00,1.00);
	x2 = Vector2D( 1.00,2.00);
	Bnd_Spline_East.Create_Spline_Line(x1,x2,2);
	x1 = Vector2D( 0.00,1.00);
	x2 = Vector2D( 0.00,2.00);
	Bnd_Spline_West.Create_Spline_Line(x1,x2,2);
      } else if (iBlk == 2 && jBlk == 1) {
	// Create the splines defining the north, south, east, and west 
	// boundaries of the grid.
	x1 = Vector2D( 1.00,2.00);
	x2 = Vector2D( 2.75,2.00);
	Bnd_Spline_North.Create_Spline_Line(x1,x2,2);
	x1 = Vector2D( 1.00,1.00);
	x2 = Vector2D( 2.75,1.00);
	Bnd_Spline_South.Create_Spline_Line(x1,x2,2);
	x1 = Vector2D( 2.75,1.00);
	x2 = Vector2D( 2.75,2.00);
	Bnd_Spline_East.Create_Spline_Line(x1,x2,2);
	x1 = Vector2D( 1.00,1.00);
	x2 = Vector2D( 1.00,2.00);
	Bnd_Spline_West.Create_Spline_Line(x1,x2,2);
      } else if (iBlk == 3 && jBlk == 1) {
	// Create the splines defining the north, south, east, and west 
	// boundaries of the grid.
	x1 = Vector2D( 2.75,2.00);
	x2 = Vector2D( 4.50,2.00);
	Bnd_Spline_North.Create_Spline_Line(x1,x2,2);
	x1 = Vector2D( 2.75,1.00);
	x2 = Vector2D( 4.50,1.00);
	Bnd_Spline_South.Create_Spline_Line(x1,x2,2);
	x1 = Vector2D( 4.50,1.00);
	x2 = Vector2D( 4.50,2.00);
	Bnd_Spline_East.Create_Spline_Line(x1,x2,2);
	x1 = Vector2D( 2.75,1.00);
	x2 = Vector2D( 2.75,2.00);
	Bnd_Spline_West.Create_Spline_Line(x1,x2,2);
      }

      // Set the boundary condition types for each of the boundary splines.
      if (jBlk == 0) Bnd_Spline_South.setBCtype(BC_REFLECTION);
      else Bnd_Spline_South.setBCtype(BC_NONE);
      if (jBlk == 1) Bnd_Spline_North.setBCtype(BC_REFLECTION);
      else Bnd_Spline_North.setBCtype(BC_NONE);
      if (iBlk == 0) Bnd_Spline_West.setBCtype(BC_FIXED);
      else Bnd_Spline_West.setBCtype(BC_NONE);
      //if (iBlk == 3) Bnd_Spline_East.setBCtype(BC_FIXED);
      //if (iBlk == 3) Bnd_Spline_East.setBCtype(BC_CHARACTERISTIC);
      if (iBlk == 3) Bnd_Spline_East.setBCtype(BC_CONSTANT_EXTRAPOLATION);
      else Bnd_Spline_East.setBCtype(BC_NONE);

      // Assign values to the stretching function parameters and boundary 
      // grid line orthogonality parameters.
      Stretch_I = STRETCHING_FCN_LINEAR;
      Beta_I = ZERO; 
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_LINEAR;
      Beta_J = ZERO;
      Tau_J = ZERO;
      Orthogonal_North = 1;
      Orthogonal_South = 1;
      Orthogonal_East = 1;
      Orthogonal_West = 1;

      // Create the 2D quadrilateral grid block representing the mesh.
      Grid_ptr[iBlk][jBlk].Create_Quad_Block_Without_Update(Bnd_Spline_North,
							    Bnd_Spline_South,
							    Bnd_Spline_East,
							    Bnd_Spline_West,
							    Number_of_Cells_Idir,
							    Number_of_Cells_Jdir,
							    Number_of_Ghost_Cells,
							    Highest_Order_of_Reconstruction,
							    GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
							    Stretch_I,
							    Beta_I, 
							    Tau_I,
							    Stretch_J,
							    Beta_J,
							    Tau_J,
							    Orthogonal_North,
							    Orthogonal_South,
							    Orthogonal_East,
							    Orthogonal_West);

      // Smooth the 2D quadrilateral grid block.
      Grid_ptr[iBlk][jBlk].Smooth_Quad_Block(min(250,2*max(Number_of_Cells_Idir,Number_of_Cells_Jdir)));
      
    }
  }

}

/*!
 * Generates a quadilateral mesh with for predicting jet flow coming  
 * out through an orifice.                                            
 *                                                                    
 * Usage: Grid_ptr = Grid_Jet_Flow(Grid_ptr,                          
 *                                 nblk_i,                            
 *                                 nblk_j,                            
 *                                 plate_length,                      
 *				   Stretching_Type_Idir,              
 *				   Stretching_Type_Jdir,              
 *				   Stretching_Factor_Idir,              
 *				   Stretching_Factor_Jdir,            
 *                                 Number_of_Cells_Idir,              
 *                                 Number_of_Cells_Jdir,              
 *                                 2)                                 
 * This subroutine DOESN'T update the ghost cells or
 * the geometric properties of the grid cells.
 */
void Grid2D_Quad_MultiBlock_HO::Grid_Jet_Flow_Without_Update(int &_Number_of_Blocks_Idir_,
							     int &_Number_of_Blocks_Jdir_,
							     const double &Radius,
							     const double &Mach,
							     const int &Stretching_Type_Idir,
							     const int &Stretching_Type_Jdir,
							     const double &Stretching_Factor_Idir,
							     const double &Stretching_Factor_Jdir,
							     const int Number_of_Cells_Idir,
							     const int Number_of_Cells_Jdir,
							     const int Number_of_Ghost_Cells,
							     const int Highest_Order_of_Reconstruction){
  
  int block_flag,
      BCtypeN, BCtypeS, BCtypeE, BCtypeW,
      Stretch_I, Stretch_J,
      Orthogonal_North, Orthogonal_South,
      Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D xNW, xNE, xSE, xSW;
  Spline2D_HO Bnd_Spline_North, Bnd_Spline_South,
           Bnd_Spline_East, Bnd_Spline_West;
 
  // Allocate memory for grid blocks. There are two grid blocks for 
  // this mesh.

  allocate(_Number_of_Blocks_Idir_,_Number_of_Blocks_Jdir_);
  // Create the mesh for each block representing the complete grid.
  for (int jBlk = 0 ; jBlk < Number_of_Blocks_Jdir ; jBlk++) {
    for (int iBlk = 0; iBlk < Number_of_Blocks_Idir ; iBlk++) {
      // Preset the block indicator and turn off when required:
      block_flag = ON;

       // Preset all boundary conditions to none and reset when required:
      BCtypeN = BC_NONE;
      BCtypeS = BC_NONE;
      BCtypeE = BC_NONE;
      BCtypeW = BC_NONE;    
 
      // Preset all stretching parameters and reset when required:
      Stretch_I = STRETCHING_FCN_LINEAR;
      Beta_I = ZERO;
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_LINEAR;
      Beta_J = ZERO;
      Tau_J = ZERO;
      Orthogonal_North = 0;
      Orthogonal_South = 0;
      Orthogonal_East = 0;
      Orthogonal_West = 0;

      if (Mach < 1.0){
	double Width_Domain = Radius*25.0; // Width of the domain
	double Len = 65.0*Radius; // Length of the domain
	double height_Block = (Width_Domain - 2*Radius)/(pow(TWO,Number_of_Blocks_Jdir-2)-ONE);
	double height_Block_sup = (2.0*Width_Domain - 2*Radius)/(pow(TWO,Number_of_Blocks_Jdir-2)-ONE);
	double ele = Len/(pow(TWO,Number_of_Blocks_Idir)-ONE);
	if (jBlk==0){
	  xSW = Vector2D(ele*(pow(TWO,iBlk)-1.0), ZERO);
	  xSE = Vector2D(ele*(pow(TWO,iBlk+1.0)-1.0), ZERO);
	  xNW = Vector2D(ele*(pow(TWO,iBlk)-1.0), Radius);
	  xNE = Vector2D(ele*(pow(TWO,iBlk+1.0)-1.0), Radius);
	  Stretch_J = STRETCHING_FCN_MAX_CLUSTERING;
	  Beta_J = Stretching_Factor_Jdir;
	}else {
	  if (jBlk ==1){
	    xSW = Vector2D(ele*(pow(TWO,iBlk)-1.0), Radius);
	    xSE = Vector2D(ele*(pow(TWO,iBlk+1.0)-1.0), Radius);
	    xNW = Vector2D(ele*(pow(TWO,iBlk)-1.0), TWO*Radius);
	    xNE = Vector2D(ele*(pow(TWO,iBlk+1)-1.0), TWO*Radius);
	    Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
	    Beta_J = Stretching_Factor_Jdir; 
	  }else{
	    xSW = Vector2D(ele*(pow(TWO,iBlk)-1.0), TWO*Radius + height_Block*(pow(TWO,jBlk-2)-1.0) );
	    xSE = Vector2D(ele*(pow(TWO,iBlk+1.0)-1.0), TWO*Radius + height_Block*(pow(TWO,jBlk-2)-1.0));
	    xNW = Vector2D(ele*(pow(TWO,iBlk)-1.0), TWO*Radius + height_Block*(pow(TWO,jBlk-1)-1.0));
	    xNE = Vector2D(ele*(pow(TWO,iBlk+1)-1.0), TWO*Radius + height_Block*(pow(TWO,jBlk-1)-1.0));
	  }
	  
	}
	
	if (iBlk==0){
	  Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
	  Beta_I = Stretching_Factor_Idir;
	}
	if (jBlk==2){
	  Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
	  Beta_J = (Stretching_Factor_Jdir-1.0)/2.0+1.0; 
	}

	//Specifying boundary conditions
	// 	if (jBlk==0 && iBlk ==0)
	// 	  BCtypeW = BC_FIXED; 
	// 	if (jBlk >0 && iBlk ==0)
	// 	  BCtypeW = BC_WALL_VISCOUS_ISOTHERMAL; 
	// 	if (jBlk == 0)
	// 	  BCtypeS = BC_REFLECTION;
	// 	if (jBlk == Number_of_Blocks_Jdir-1)
	// 	  BCtypeN = BC_FIXED;
	// 	if (iBlk ==  Number_of_Blocks_Idir-1)
	// 	  BCtypeE = BC_CONSTANT_EXTRAPOLATION;
	//Boundary condition specified

      }else{
	double Width_Domain = Radius*25.0; // Width of the domain
	double Len = 100.0*Radius; // Length of the domain
	double height_Block = (Width_Domain - 2*Radius)/(pow(TWO,Number_of_Blocks_Jdir-2)-ONE);
	double height_Block_sup = (2.0*Width_Domain - 2*Radius)/(pow(TWO,Number_of_Blocks_Jdir-2)-ONE);
	double ele = Len/(pow(TWO,Number_of_Blocks_Idir)-ONE);
	
	if (jBlk==0){
	  xSW = Vector2D(ele*(pow(TWO,iBlk)-1.0), ZERO);
	  xSE = Vector2D(ele*(pow(TWO,iBlk+1.0)-1.0), ZERO);
	  xNW = Vector2D(ele*(pow(TWO,iBlk)-1.0), Radius);
	  xNE = Vector2D(ele*(pow(TWO,iBlk+1.0)-1.0), Radius);
	  Stretch_J = STRETCHING_FCN_MAX_CLUSTERING;
	  Beta_J = Stretching_Factor_Jdir;
	}else if (jBlk ==1){
	    xSW = Vector2D(ele*(pow(TWO,iBlk)-1.0), Radius);
	    xSE = Vector2D(ele*(pow(TWO,iBlk+1.0)-1.0), Radius);
	    xNW = Vector2D(ele*(pow(TWO,iBlk)-1.0), TWO*Radius);
	    xNE = Vector2D(ele*(pow(TWO,iBlk+1)-1.0), TWO*Radius);
	    Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
	    Beta_J = Stretching_Factor_Jdir; 
	}else if (jBlk == 2){
	  xSW = Vector2D(ele*(pow(TWO,iBlk)-1.0), TWO*Radius + height_Block*(pow(TWO,jBlk-2)-1.0) );
	  xSE = Vector2D(ele*(pow(TWO,iBlk+1.0)-1.0), TWO*Radius + height_Block*(pow(TWO,jBlk-2)-1.0));
	  xNW = Vector2D(ele*(pow(TWO,iBlk)-1.0), TWO*Radius +  height_Block*(pow(TWO,jBlk-1)-1.0) + ele*(pow(TWO,iBlk)-1.0)/Len*(height_Block_sup*(pow(TWO,jBlk-1)-1.0)-height_Block*(pow(TWO,jBlk-1)-1.0)));
	  xNE = Vector2D(ele*(pow(TWO,iBlk+1)-1.0), TWO*Radius + height_Block*(pow(TWO,jBlk-1)-1.0) + ele*(pow(TWO,iBlk+1.0)-1.0)/Len*(height_Block_sup*(pow(TWO,jBlk-1)-1.0)-height_Block*(pow(TWO,jBlk-1)-1.0)));
	}else if (jBlk == 3){
	  xSW = Vector2D(ele*(pow(TWO,iBlk)-1.0), TWO*Radius +  height_Block*(pow(TWO,jBlk-2)-1.0) + ele*(pow(TWO,iBlk)-1.0)/Len*(height_Block_sup*(pow(TWO,jBlk-2)-1.0)-height_Block*(pow(TWO,jBlk-2)-1.0)));
	  xSE = Vector2D(ele*(pow(TWO,iBlk+1)-1.0), TWO*Radius + height_Block*(pow(TWO,jBlk-2)-1.0) + ele*(pow(TWO,iBlk+1.0)-1.0)/Len*(height_Block_sup*(pow(TWO,jBlk-2)-1.0)-height_Block*(pow(TWO,jBlk-2)-1.0)));	 
	  xNW = Vector2D(ele*(pow(TWO,iBlk)-1.0), TWO*Radius +  height_Block*(pow(TWO,jBlk-1)-1.0) + ele*(pow(TWO,iBlk)-1.0)/Len*(height_Block_sup*(pow(TWO,jBlk-1)-1.0)-height_Block*(pow(TWO,jBlk-1)-1.0)));
	  xNE = Vector2D(ele*(pow(TWO,iBlk+1)-1.0), TWO*Radius + height_Block*(pow(TWO,jBlk-1)-1.0) + ele*(pow(TWO,iBlk+1.0)-1.0)/Len*(height_Block_sup*(pow(TWO,jBlk-1)-1.0)-height_Block*(pow(TWO,jBlk-1)-1.0)));
	}
	if (iBlk==0){
	  Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
	  Beta_I = Stretching_Factor_Idir; 
	}
	if (jBlk==2){
	  Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
	  Beta_J = (Stretching_Factor_Jdir-1.0)/2.0+1.0; 
	}
	
	//Specifying boundary conditions
	if (jBlk==0 && iBlk ==0)
	  BCtypeW = BC_FIXED; 
	if (jBlk!=0 && iBlk ==0)
	  BCtypeW = BC_REFLECTION; 
	if (jBlk == 0)
	  BCtypeS = BC_REFLECTION;
	if (jBlk == Number_of_Blocks_Jdir-1)
	  BCtypeN = BC_FIXED;
	if (iBlk == Number_of_Blocks_Idir-1)
	  BCtypeE = BC_CONSTANT_EXTRAPOLATION;
	//Boundary condition specified
	
      }
      // Determine the block indicator and all boundary splines and
      // reset all boundary conditions and stretching parameters for the
      // current block: 
      
      if (block_flag) {
	
	// Create the splines defining the north, south, east, and west 
	// boundaries of the rectangular boxes.
	Bnd_Spline_North.Create_Spline_Line(xNW,xNE,2);
	Bnd_Spline_South.Create_Spline_Line(xSW,xSE,2);
	Bnd_Spline_East.Create_Spline_Line(xSE,xNE,2);
	Bnd_Spline_West.Create_Spline_Line(xSW,xNW,2);
	
	// Set the boundary condition types for each of the boundary
	// splines:
	Bnd_Spline_North.setBCtype(BCtypeN);
	Bnd_Spline_South.setBCtype(BCtypeS);
	Bnd_Spline_East.setBCtype(BCtypeE);
	Bnd_Spline_West.setBCtype(BCtypeW);
	
	// Create the 2D quadrilateral grid block.
	Grid_ptr[iBlk][jBlk].Create_Quad_Block_Without_Update(Bnd_Spline_North,
							      Bnd_Spline_South,
							      Bnd_Spline_East,
							      Bnd_Spline_West,
							      Number_of_Cells_Idir,
							      Number_of_Cells_Jdir,
							      Number_of_Ghost_Cells,
							      Highest_Order_of_Reconstruction,
							      GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
							      Stretch_I,
							      Beta_I,
							      Tau_I,
							      Stretch_J,
							      Beta_J,
							      Tau_J,
							      Orthogonal_North,
							      Orthogonal_South,
							      Orthogonal_East,
							      Orthogonal_West);
	
      }
      
    }    
  }
  
}

/*!
 * Generates a quadilateral mesh with for predicting Mixing Layer flow
 * coming out through an orifice.                                     
 *                                                                    
 * Usage: Grid_ptr = Grid_Mixing_Layer(Grid_ptr,                      
 *                                 nblk_i,                            
 *                                 nblk_j,                            
 *                                 plate_length,                      
 *				   Stretching_Type_Idir,              
 *				   Stretching_Type_Jdir,              
 *				   Stretching_Factor_Idir,              
 *				   Stretching_Factor_Jdir,            
 *                                 Number_of_Cells_Idir,              
 *                                 2)                                 
 * This subroutine DOESN'T update the ghost cells or
 * the geometric properties of the grid cells.
 */
void Grid2D_Quad_MultiBlock_HO::Grid_Mixing_Layer_Without_Update(int &_Number_of_Blocks_Idir_,
								 int &_Number_of_Blocks_Jdir_,
								 const double &Length,
								 const double &Mach,
								 const int &Stretching_Type_Idir,
								 const int &Stretching_Type_Jdir,
								 const double &Stretching_Factor_Idir,
								 const double &Stretching_Factor_Jdir,
								 const int Number_of_Cells_Idir,
								 const int Number_of_Cells_Jdir,
								 const int Number_of_Ghost_Cells,
								 const int Highest_Order_of_Reconstruction){
  
  int block_flag,
      BCtypeN, BCtypeS, BCtypeE, BCtypeW,
      Stretch_I, Stretch_J,
      Orthogonal_North, Orthogonal_South,
      Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D xNW, xNE, xSE, xSW;
  Spline2D_HO Bnd_Spline_North, Bnd_Spline_South,
           Bnd_Spline_East, Bnd_Spline_West;

  double factor = ONE; // this factor increases the length of the domain calculated
  // by the given Reynolds number 'factor' times. 

  double Fside = TWO; // this factor decides the width of the domain.
  //'Fside = fside' implies width of the domain on either side of centerline
  //will be 'factor*Length/fside'

  // Allocate memory for grid blocks. There are two grid blocks for 
  // this mesh.

  double len = Length*factor; // final length of the domain.
  double width = len/Fside;  // width on each side of the centerline.

  len = 2.0;
  width = 1.0;

  allocate(_Number_of_Blocks_Idir_,_Number_of_Blocks_Jdir_);
  // Create the mesh for each block representing the complete grid.
  
  for (int jBlk = 0; jBlk < Number_of_Blocks_Jdir; jBlk++) {
    for (int iBlk = 0; iBlk < Number_of_Blocks_Idir; iBlk++) {
      
      // Preset the block indicator and turn off when required:
       block_flag = ON;

       // Preset all boundary conditions to none and reset when required:
      BCtypeN = BC_NONE;
      BCtypeS = BC_NONE;
      BCtypeE = BC_NONE;
      BCtypeW = BC_NONE;    
 
      // Preset all stretching parameters and reset when required:
      Stretch_I = STRETCHING_FCN_LINEAR;
      Beta_I = ZERO;
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_LINEAR;
      Beta_J = ZERO;
      Tau_J = ZERO;
      Orthogonal_North = 0;
      Orthogonal_South = 0;
      Orthogonal_East = 0;
      Orthogonal_West = 0;

      //       if ( jBlk < Number_of_Blocks_Jdir/TWO ){
      // 	xSW = Vector2D(iBlk*len/Number_of_Blocks_Idir,-(pow(TWO,Number_of_Blocks_Jdir/TWO - jBlk)-ONE)*width/(pow(TWO,Number_of_Blocks_Jdir/TWO)-ONE));
      // 	xSE = Vector2D((iBlk+ONE)*len/Number_of_Blocks_Idir,-(pow(TWO,Number_of_Blocks_Jdir/TWO - jBlk)-ONE)*width/(pow(TWO,Number_of_Blocks_Jdir/TWO)-ONE));
      // 	xNW = Vector2D(iBlk*len/Number_of_Blocks_Idir,-(pow(TWO,Number_of_Blocks_Jdir/TWO-jBlk-ONE ) -ONE)*width/(pow(TWO,Number_of_Blocks_Jdir/TWO)-ONE));
      // 	xNE=Vector2D((iBlk+ONE)*len/Number_of_Blocks_Idir,-(pow(TWO,Number_of_Blocks_Jdir/TWO-jBlk-ONE)-ONE)*width/(pow(TWO,Number_of_Blocks_Jdir/TWO)-ONE));
      // 	if ( jBlk == 0 )
      // 	  BCtypeS = BC_CONSTANT_EXTRAPOLATION;
      // 	if ( iBlk == 0 )
      // 	  BCtypeW = BC_FIXED;
      // 	if ( iBlk == Number_of_Blocks_Idir - 1 )
      // 	  BCtypeE = BC_FIXED_PRESSURE;
	
      // 	if (jBlk == Number_of_Blocks_Jdir / 2 - 1 ){
      // 	  Stretch_J = STRETCHING_FCN_MAX_CLUSTERING;
      // 	  Beta_J = Stretching_Factor_Jdir; 
      // 	}
	
      //       }else{
      // 	xSW = Vector2D(iBlk*len/Number_of_Blocks_Idir,(pow(TWO,jBlk - Number_of_Blocks_Jdir/TWO) -ONE)*width/(pow(TWO,Number_of_Blocks_Jdir/TWO)-ONE));
      // 	xSE = Vector2D((iBlk+ONE)*len/Number_of_Blocks_Idir,(pow(TWO,jBlk -Number_of_Blocks_Jdir/TWO) -ONE)*width/(pow(TWO,Number_of_Blocks_Jdir/TWO)-ONE));
      // 	xNW = Vector2D(iBlk*len/Number_of_Blocks_Idir,(pow(TWO,jBlk-Number_of_Blocks_Jdir/TWO + ONE )-ONE)*width/(pow(TWO,Number_of_Blocks_Jdir/TWO)-ONE) );
      // 	xNE =Vector2D((iBlk+ONE)*len/Number_of_Blocks_Idir,(pow(TWO,jBlk-Number_of_Blocks_Jdir/TWO+ONE)-ONE)*width/(pow(TWO,Number_of_Blocks_Jdir/TWO)-ONE));
      // 	if ( jBlk == Number_of_Blocks_Jdir - 1 )
      // 	  BCtypeN = BC_CONSTANT_EXTRAPOLATION;
      // 	if ( iBlk == 0 )
      // 	  BCtypeW = BC_FIXED;
      // 	if ( iBlk == Number_of_Blocks_Idir - 1 ){
      // 	  if (Mach < ONE)
      // 	    BCtypeE = BC_FIXED_PRESSURE;
      // 	  else
      // 	    BCtypeE = BC_CONSTANT_EXTRAPOLATION;
      // 	}
      // 	if (jBlk == Number_of_Blocks_Jdir / 2 ){
      // 	  Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
      // 	  Beta_J = Stretching_Factor_Jdir; 
      // 	}
      //       }
      
      double w_part = width/(Number_of_Blocks_Idir/TWO);
      double l_part = len/(Number_of_Blocks_Jdir/TWO);

      xSW = Vector2D((iBlk-Number_of_Blocks_Idir/TWO)*w_part,(jBlk-Number_of_Blocks_Jdir/TWO)*l_part);
      xSE = Vector2D((iBlk+1.0-Number_of_Blocks_Idir/TWO)*w_part,(jBlk-Number_of_Blocks_Jdir/TWO)*l_part);
      xNW = Vector2D((iBlk-Number_of_Blocks_Idir/TWO)*w_part,(jBlk+1.0-Number_of_Blocks_Jdir/TWO)*l_part);
      xNE = Vector2D((iBlk+1.0-Number_of_Blocks_Idir/TWO)*w_part,(jBlk+1.0-Number_of_Blocks_Jdir/TWO)*l_part);

      if (jBlk == 0)
	BCtypeS = BC_REFLECTION;
      
      if(jBlk == Number_of_Blocks_Jdir - 1)
	BCtypeN = BC_REFLECTION;
    
      if(iBlk == Number_of_Blocks_Idir - 1)
	BCtypeE = BC_CHARACTERISTIC;

      if (iBlk ==0){
	if(jBlk < Number_of_Blocks_Jdir/TWO){
	  BCtypeW = BC_INFLOW_SUBSONIC;
	}else{
	  BCtypeW = BC_CHARACTERISTIC;
	}
      }

      if (jBlk == Number_of_Blocks_Jdir/TWO-1){
	  Stretch_J = STRETCHING_FCN_MAX_CLUSTERING;
 	  Beta_J = Stretching_Factor_Jdir; 
	  if (iBlk < Number_of_Blocks_Idir/TWO){
	    BCtypeN = BC_WALL_VISCOUS_HEATFLUX;
	  }
      }

      if (jBlk == Number_of_Blocks_Jdir/TWO){
	  Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
 	  Beta_J = Stretching_Factor_Jdir; 
	  if (iBlk < Number_of_Blocks_Idir/TWO){
	    BCtypeS = BC_WALL_VISCOUS_HEATFLUX;
	  }
      }
      
      if (block_flag) {
	
	// Create the splines defining the north, south, east, and west 
	// boundaries of the rectangular boxes.
	Bnd_Spline_North.Create_Spline_Line(xNW,xNE,2);
	Bnd_Spline_South.Create_Spline_Line(xSW,xSE,2);
	Bnd_Spline_East.Create_Spline_Line(xSE,xNE,2);
	Bnd_Spline_West.Create_Spline_Line(xSW,xNW,2);

	// Set the boundary condition types for each of the boundary
	// splines:
	Bnd_Spline_North.setBCtype(BCtypeN);
	Bnd_Spline_South.setBCtype(BCtypeS);
	Bnd_Spline_East.setBCtype(BCtypeE);
	Bnd_Spline_West.setBCtype(BCtypeW);
    
	// Create the 2D quadrilateral grid block.
	Grid_ptr[iBlk][jBlk].Create_Quad_Block_Without_Update(Bnd_Spline_North,
							      Bnd_Spline_South,
							      Bnd_Spline_East,
							      Bnd_Spline_West,
							      Number_of_Cells_Idir,
							      Number_of_Cells_Jdir,
							      Number_of_Ghost_Cells,
							      Highest_Order_of_Reconstruction,
							      GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
							      Stretch_I,
							      Beta_I,
							      Tau_I,
							      Stretch_J,
							      Beta_J,
							      Tau_J,
							      Orthogonal_North,
							      Orthogonal_South,
							      Orthogonal_East,
							      Orthogonal_West);
	
      }

    }    
  }

}

/*!
 * Generates a quadilateral mesh with clustering consisting of five   
 * grid blocks for predicting flow over a backward facing step.       
 * This subroutine DOESN'T update the ghost cells or
 * the geometric properties of the grid cells.
 **/
void Grid2D_Quad_MultiBlock_HO::Grid_Backward_Facing_Step_Without_Update(int &_Number_of_Blocks_Idir_,
									 int &_Number_of_Blocks_Jdir_,
									 const double &Step_Height,
									 const double &Top_Wall_Deflection,
									 const double &Stretching_Factor_Idir,
									 const double &Stretching_Factor_Jdir,
									 const int Number_of_Cells_Idir,
									 const int Number_of_Cells_Jdir,
									 const int Number_of_Ghost_Cells,
									 const int Highest_Order_of_Reconstruction) {

  int block_flag,
      BCtypeN, BCtypeS, BCtypeE, BCtypeW,
      Stretch_I, Stretch_J,
      Orthogonal_North, Orthogonal_South,
      Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D xNW, xNE, xSE, xSW;
  Spline2D_HO Bnd_Spline_North, Bnd_Spline_South,
           Bnd_Spline_East, Bnd_Spline_West;
  double Inlet_Height, Inlet_Length, Outlet_Length;

  // Set geometric constants.
  Inlet_Height = 8.0*Step_Height;
  Inlet_Length = 20.0*Step_Height;
  Outlet_Length = 50.0*Step_Height;

  // Allocate memory for grid blocks.  There are two grid blocks for 
  // this mesh.
  _Number_of_Blocks_Idir_ = 7;
  _Number_of_Blocks_Jdir_ = 3;
  allocate(_Number_of_Blocks_Idir_,_Number_of_Blocks_Jdir_);

  // Create the mesh for each block representing the complete grid.
  for (int jBlk = 0; jBlk < Number_of_Blocks_Jdir; jBlk++) {
    for (int iBlk = 0; iBlk < Number_of_Blocks_Idir; iBlk++) {

      // Preset the block indicator and turn off when required:
      block_flag = ON;

      // Preset all boundary conditions to none and reset when required:
      BCtypeN = BC_NONE;
      BCtypeS = BC_NONE;
      BCtypeE = BC_NONE;
      BCtypeW = BC_NONE;

      // Preset all stretching parameters and reset when required:
      Stretch_I = STRETCHING_FCN_LINEAR;
      Beta_I = ZERO;
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_LINEAR;
      Beta_J = ZERO;
      Tau_J = ZERO;
      Orthogonal_North = 0;
      Orthogonal_South = 0;
      Orthogonal_East = 0;
      Orthogonal_West = 0;

      // Determine the block indicator and all boundary splines and
      // reset all boundary conditions and stretching parameters for the
      // current block: 
      if (iBlk == 0 && jBlk == 0) {
	block_flag = OFF;
      } else if (iBlk == 1 && jBlk == 0) {
	block_flag = OFF;
      } else if (iBlk == 2 && jBlk == 0) {
	xSW = Vector2D(ZERO,ZERO);
	xSE = Vector2D(0.2*Outlet_Length,ZERO);
	xNW = Vector2D(ZERO,Step_Height);
	xNE = Vector2D(0.2*Outlet_Length,Step_Height);
	BCtypeS = BC_WALL_VISCOUS_HEATFLUX;
	BCtypeW = BC_WALL_VISCOUS_HEATFLUX;
	Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_I = 1.10;
	Stretch_J = STRETCHING_FCN_MINMAX_CLUSTERING;
	Beta_J = 1.10;
      } else if (iBlk == 3 && jBlk == 0) {
	xSW = Vector2D(0.2*Outlet_Length,ZERO);
	xSE = Vector2D(0.4*Outlet_Length,ZERO);
	xNW = Vector2D(0.2*Outlet_Length,Step_Height);
	xNE = Vector2D(0.4*Outlet_Length,Step_Height);
	BCtypeS = BC_WALL_VISCOUS_HEATFLUX;
	Stretch_J = STRETCHING_FCN_MINMAX_CLUSTERING;
	Beta_J = 1.10;
      } else if (iBlk == 4 && jBlk == 0) {
	xSW = Vector2D(0.4*Outlet_Length,ZERO);
	xSE = Vector2D(0.6*Outlet_Length,ZERO);
	xNW = Vector2D(0.4*Outlet_Length,Step_Height);
	xNE = Vector2D(0.6*Outlet_Length,Step_Height);
	BCtypeS = BC_WALL_VISCOUS_HEATFLUX;
	Stretch_J = STRETCHING_FCN_MINMAX_CLUSTERING;
	Beta_J = 1.10;
      } else if (iBlk == 5 && jBlk == 0) {
	xSW = Vector2D(0.6*Outlet_Length,ZERO);
	xSE = Vector2D(0.8*Outlet_Length,ZERO);
	xNW = Vector2D(0.6*Outlet_Length,Step_Height);
	xNE = Vector2D(0.8*Outlet_Length,Step_Height);
	BCtypeS = BC_WALL_VISCOUS_HEATFLUX;
	Stretch_J = STRETCHING_FCN_MINMAX_CLUSTERING;
	Beta_J = 1.10;
      } else if (iBlk == 6 && jBlk == 0) {
	xSW = Vector2D(0.8*Outlet_Length,ZERO);
	xSE = Vector2D(Outlet_Length,ZERO);
	xNW = Vector2D(0.8*Outlet_Length,Step_Height);
	xNE = Vector2D(Outlet_Length,Step_Height);
	BCtypeS = BC_WALL_VISCOUS_HEATFLUX;
	BCtypeE = BC_OUTFLOW_SUBSONIC;
	Stretch_J = STRETCHING_FCN_MINMAX_CLUSTERING;
	Beta_J = 1.10;

      } else if (iBlk == 0 && jBlk == 1) {
	xSW = Vector2D(-Inlet_Length,Step_Height);
	xSE = Vector2D(-HALF*Inlet_Length,Step_Height);
	xNW = Vector2D(-Inlet_Length,Step_Height+HALF*Inlet_Height);
	xNE = Vector2D(-HALF*Inlet_Length,Step_Height+HALF*Inlet_Height);
	BCtypeS = BC_WALL_VISCOUS_HEATFLUX;
	BCtypeW = BC_INFLOW_SUBSONIC;
	Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_J = 1.10;
      } else if (iBlk == 1 && jBlk == 1) {
	xSW = Vector2D(-HALF*Inlet_Length,Step_Height);
	xSE = Vector2D(ZERO,Step_Height);
	xNW = Vector2D(-HALF*Inlet_Length,Step_Height+HALF*Inlet_Height);
	xNE = Vector2D(ZERO,Step_Height+HALF*Inlet_Height);
	BCtypeS = BC_WALL_VISCOUS_HEATFLUX;
	Stretch_I = STRETCHING_FCN_MAX_CLUSTERING;
	Beta_I = 1.10;
	Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_J = 1.10;
      } else if (iBlk == 2 && jBlk == 1) {
	xSW = Vector2D(ZERO,Step_Height);
	xSE = Vector2D(0.2*Outlet_Length,Step_Height);
	xNW = Vector2D(ZERO,Step_Height+HALF*Inlet_Height);
	xNE = Vector2D(0.2*Outlet_Length,Step_Height+HALF*Inlet_Height);
	Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_I = 1.10;
	Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_J = 1.10;
      } else if (iBlk == 3 && jBlk == 1) {
	xSW = Vector2D(0.2*Outlet_Length,Step_Height);
	xSE = Vector2D(0.4*Outlet_Length,Step_Height);
	xNW = Vector2D(0.2*Outlet_Length,Step_Height+HALF*Inlet_Height);
	xNE = Vector2D(0.4*Outlet_Length,Step_Height+HALF*Inlet_Height);
	Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_J = 1.10;
      } else if (iBlk == 4 && jBlk == 1) {
	xSW = Vector2D(0.4*Outlet_Length,Step_Height);
	xSE = Vector2D(0.6*Outlet_Length,Step_Height);
	xNW = Vector2D(0.4*Outlet_Length,Step_Height+HALF*Inlet_Height);
	xNE = Vector2D(0.6*Outlet_Length,Step_Height+HALF*Inlet_Height);
	Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_J = 1.10;
      } else if (iBlk == 5 && jBlk == 1) {
	xSW = Vector2D(0.6*Outlet_Length,Step_Height);
	xSE = Vector2D(0.8*Outlet_Length,Step_Height);
	xNW = Vector2D(0.6*Outlet_Length,Step_Height+HALF*Inlet_Height);
	xNE = Vector2D(0.8*Outlet_Length,Step_Height+HALF*Inlet_Height);
	Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_J = 1.10;
      } else if (iBlk == 6 && jBlk == 1) {
	xSW = Vector2D(0.8*Outlet_Length,Step_Height);
	xSE = Vector2D(Outlet_Length,Step_Height);
	xNW = Vector2D(0.8*Outlet_Length,Step_Height+HALF*Inlet_Height);
	xNE = Vector2D(Outlet_Length,Step_Height+HALF*Inlet_Height);
	BCtypeE = BC_OUTFLOW_SUBSONIC;
	Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_J = 1.10;

      } else if (iBlk == 0 && jBlk == 2) {
	xSW = Vector2D(-Inlet_Length,Step_Height+HALF*Inlet_Height);
	xSE = Vector2D(-HALF*Inlet_Length,Step_Height+HALF*Inlet_Height);
	xNW = Vector2D(-Inlet_Length,Step_Height+Inlet_Height);
	xNE = Vector2D(-HALF*Inlet_Length,Step_Height+Inlet_Height);
	BCtypeN = BC_WALL_VISCOUS_HEATFLUX;
	BCtypeW = BC_INFLOW_SUBSONIC;
	Stretch_J = STRETCHING_FCN_MAX_CLUSTERING;
	Beta_J = 1.10;
      } else if (iBlk == 1 && jBlk == 2) {
	xSW = Vector2D(-HALF*Inlet_Length,Step_Height+HALF*Inlet_Height);
	xSE = Vector2D(ZERO,Step_Height+HALF*Inlet_Height);
	xNW = Vector2D(-HALF*Inlet_Length,Step_Height+Inlet_Height);
	xNE = Vector2D(ZERO,Step_Height+Inlet_Height);
	BCtypeN = BC_WALL_VISCOUS_HEATFLUX;
	Stretch_I = STRETCHING_FCN_MAX_CLUSTERING;
	Beta_I = 1.10;
	Stretch_J = STRETCHING_FCN_MAX_CLUSTERING;
	Beta_J = 1.10;
      } else if (iBlk == 2 && jBlk == 2) {
	xSW = Vector2D(ZERO,Step_Height+HALF*Inlet_Height);
	xSE = Vector2D(0.2*Outlet_Length,Step_Height+HALF*Inlet_Height);
	xNW = Vector2D(ZERO,Step_Height+Inlet_Height);
	xNE = Vector2D(0.2*Outlet_Length,Step_Height+Inlet_Height+0.20*Outlet_Length*tan(Top_Wall_Deflection*PI/180.0));
	BCtypeN = BC_WALL_VISCOUS_HEATFLUX;
	Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_I = 1.10;
	Stretch_J = STRETCHING_FCN_MAX_CLUSTERING;
	Beta_J = 1.10;
      } else if (iBlk == 3 && jBlk == 2) {
	xSW = Vector2D(0.2*Outlet_Length,Step_Height+HALF*Inlet_Height);
	xSE = Vector2D(0.4*Outlet_Length,Step_Height+HALF*Inlet_Height);
	xNW = Vector2D(0.2*Outlet_Length,Step_Height+Inlet_Height+0.20*Outlet_Length*tan(Top_Wall_Deflection*PI/180.0));
	xNE = Vector2D(0.4*Outlet_Length,Step_Height+Inlet_Height+0.40*Outlet_Length*tan(Top_Wall_Deflection*PI/180.0));
	BCtypeN = BC_WALL_VISCOUS_HEATFLUX;
	Stretch_J = STRETCHING_FCN_MAX_CLUSTERING;
	Beta_J = 1.10;
      } else if (iBlk == 4 && jBlk == 2) {
	xSW = Vector2D(0.4*Outlet_Length,Step_Height+HALF*Inlet_Height);
	xSE = Vector2D(0.6*Outlet_Length,Step_Height+HALF*Inlet_Height);
	xNW = Vector2D(0.4*Outlet_Length,Step_Height+Inlet_Height+0.40*Outlet_Length*tan(Top_Wall_Deflection*PI/180.0));
	xNE = Vector2D(0.6*Outlet_Length,Step_Height+Inlet_Height+0.60*Outlet_Length*tan(Top_Wall_Deflection*PI/180.0));
	BCtypeN = BC_WALL_VISCOUS_HEATFLUX;
	Stretch_J = STRETCHING_FCN_MAX_CLUSTERING;
	Beta_J = 1.10;
      } else if (iBlk == 5 && jBlk == 2) {
	xSW = Vector2D(0.6*Outlet_Length,Step_Height+HALF*Inlet_Height);
	xSE = Vector2D(0.8*Outlet_Length,Step_Height+HALF*Inlet_Height);
	xNW = Vector2D(0.6*Outlet_Length,Step_Height+Inlet_Height+0.60*Outlet_Length*tan(Top_Wall_Deflection*PI/180.0));
	xNE = Vector2D(0.8*Outlet_Length,Step_Height+Inlet_Height+0.80*Outlet_Length*tan(Top_Wall_Deflection*PI/180.0));
	BCtypeN = BC_WALL_VISCOUS_HEATFLUX;
	Stretch_J = STRETCHING_FCN_MAX_CLUSTERING;
	Beta_J = 1.10;
      } else if (iBlk == 6 && jBlk == 2) {
	xSW = Vector2D(0.8*Outlet_Length,Step_Height+HALF*Inlet_Height);
	xSE = Vector2D(Outlet_Length,Step_Height+HALF*Inlet_Height);
	xNW = Vector2D(0.8*Outlet_Length,Step_Height+Inlet_Height+0.80*Outlet_Length*tan(Top_Wall_Deflection*PI/180.0));
	xNE = Vector2D(Outlet_Length,Step_Height+Inlet_Height+Outlet_Length*tan(Top_Wall_Deflection*PI/180.0));
	BCtypeN = BC_WALL_VISCOUS_HEATFLUX;
	BCtypeE = BC_OUTFLOW_SUBSONIC;
	Stretch_J = STRETCHING_FCN_MAX_CLUSTERING;
	Beta_J = 1.10;

      }

      if (block_flag) {

	// Create the splines defining the north, south, east, and west 
	// boundaries of the rectangular boxes.
	Bnd_Spline_North.Create_Spline_Line(xNW,xNE,2);
	Bnd_Spline_South.Create_Spline_Line(xSW,xSE,2);
	Bnd_Spline_East.Create_Spline_Line(xSE,xNE,2);
	Bnd_Spline_West.Create_Spline_Line(xSW,xNW,2);

	// Set the boundary condition types for each of the boundary
	// splines:
	Bnd_Spline_North.setBCtype(BCtypeN);
	Bnd_Spline_South.setBCtype(BCtypeS);
	Bnd_Spline_East.setBCtype(BCtypeE);
	Bnd_Spline_West.setBCtype(BCtypeW);
    
	// Create the 2D quadrilateral grid block.
	Grid_ptr[iBlk][jBlk].Create_Quad_Block_Without_Update(Bnd_Spline_North,
							      Bnd_Spline_South,
							      Bnd_Spline_East,
							      Bnd_Spline_West,
							      Number_of_Cells_Idir,
							      Number_of_Cells_Jdir,
							      Number_of_Ghost_Cells,
							      Highest_Order_of_Reconstruction,
							      GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
							      Stretch_I,
							      Beta_I,
							      Tau_I,
							      Stretch_J,
							      Beta_J,
							      Tau_J,
							      Orthogonal_North,
							      Orthogonal_South,
							      Orthogonal_East,
							      Orthogonal_West);

      }

    }    
  }

}

/*!
 * Generates a quadilateral mesh with clustering consisting of five   
 * grid blocks for predicting flow over a forward facing step.        
 * This subroutine DOESN'T update the ghost cells or
 * the geometric properties of the grid cells.
 */
void Grid2D_Quad_MultiBlock_HO::Grid_Forward_Facing_Step_Without_Update(int &_Number_of_Blocks_Idir_,
									int &_Number_of_Blocks_Jdir_,
									const double &Step_Height,
									const double &Channel_Gap,
									const double &Stretching_Factor_Idir,
									const double &Stretching_Factor_Jdir,
									const int Number_of_Cells_Idir,
									const int Number_of_Cells_Jdir,
									const int Number_of_Ghost_Cells,
									const int Highest_Order_of_Reconstruction) {
  
  int BCtypeN, BCtypeS, BCtypeE, BCtypeW;
  Vector2D xNW, xNE, xSE, xSW;
  Spline2D_HO Bnd_Spline_North, Bnd_Spline_South,
           Bnd_Spline_East, Bnd_Spline_West;

  int Stretch_I = STRETCHING_FCN_LINEAR,
      Stretch_J = STRETCHING_FCN_LINEAR,
      Orthogonal_North = 1,
      Orthogonal_South = 1,
      Orthogonal_East = 1,
      Orthogonal_West = 1;
  double Beta_I=0.0, Tau_I=0.0, Beta_J=0.0, Tau_J=0.0;

  _Number_of_Blocks_Idir_ = 5;
  _Number_of_Blocks_Jdir_ = 5;
  allocate(_Number_of_Blocks_Idir_, _Number_of_Blocks_Jdir_);

  for (int jBlk = 0; jBlk < Number_of_Blocks_Jdir; jBlk++) {
    for (int iBlk = 0; iBlk < Number_of_Blocks_Idir; iBlk++) {

      // These blocks are to the east of the step 
      // and so are not part of the domain.
      if (iBlk >= 1 && jBlk == 0) { continue; }

      // Determine the block indicator and all boundary splines and
      // reset all boundary conditions and stretching parameters for the
      // current block: 
      xSW = Vector2D(double(iBlk)*Channel_Gap,double(jBlk)*Step_Height);
      xSE = Vector2D(double(iBlk+1)*Channel_Gap,double(jBlk)*Step_Height);
      xNW = Vector2D(double(iBlk)*Channel_Gap,double(jBlk+1)*Step_Height);
      xNE = Vector2D(double(iBlk+1)*Channel_Gap,double(jBlk+1)*Step_Height);

      BCtypeN = BC_NONE;
      BCtypeS = BC_NONE;
      BCtypeE = BC_NONE;
      BCtypeW = BC_NONE;

      if (iBlk == 0) { BCtypeW = BC_FIXED; }

      if (jBlk == 4) { BCtypeN = BC_REFLECTION; }

      if (iBlk == 4) { BCtypeE = BC_CONSTANT_EXTRAPOLATION; }

      if ((iBlk >= 1 && jBlk == 1) ||
	  (iBlk == 0 && jBlk == 0)) {
	BCtypeS = BC_REFLECTION; 
      }

      // The step:
      if (iBlk == 0 && jBlk == 0) { BCtypeE = BC_REFLECTION; }

      Bnd_Spline_North.Create_Spline_Line(xNW,xNE,2);
      Bnd_Spline_South.Create_Spline_Line(xSW,xSE,2);
      Bnd_Spline_East.Create_Spline_Line(xSE,xNE,2);
      Bnd_Spline_West.Create_Spline_Line(xSW,xNW,2);

      Bnd_Spline_North.setBCtype(BCtypeN);
      Bnd_Spline_South.setBCtype(BCtypeS);
      Bnd_Spline_East.setBCtype(BCtypeE);
      Bnd_Spline_West.setBCtype(BCtypeW);
    
      // Create the 2D quadrilateral grid block.
      Grid_ptr[iBlk][jBlk].Create_Quad_Block_Without_Update(Bnd_Spline_North,
							    Bnd_Spline_South,
							    Bnd_Spline_East,
							    Bnd_Spline_West,
							    Number_of_Cells_Idir,
							    Number_of_Cells_Jdir,
							    Number_of_Ghost_Cells,
							    Highest_Order_of_Reconstruction,
							    GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
							    Stretch_I,
							    Beta_I,
							    Tau_I,
							    Stretch_J,
							    Beta_J,
							    Tau_J,
							    Orthogonal_North,
							    Orthogonal_South,
							    Orthogonal_East,
							    Orthogonal_West);
      
    } /* endfor */    
  } /* endfor */
}

/*!
 * This routine creates a mesh corresponding to a simplified          
 * desolvation chamber designed by MDS-SCIEX.                         
 * This subroutine DOESN'T update the ghost cells or
 * the geometric properties of the grid cells.
 */
void Grid2D_Quad_MultiBlock_HO::Grid_Desolvation_Chamber_Without_Update(const int &Chamber_BC_Type,
									int &_Number_of_Blocks_Idir_,
									int &_Number_of_Blocks_Jdir_,
									const int Number_of_Cells_Idir,
									const int Number_of_Cells_Jdir,
									const int Number_of_Ghost_Cells,
									const int Highest_Order_of_Reconstruction) {
  
  int block_flag,
      BCtypeN, BCtypeS, BCtypeE, BCtypeW,
      Stretch_I, Stretch_J,
      Orthogonal_North, Orthogonal_South,
      Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D xNW, xNE, xSE, xSW;
  Spline2D_HO Bnd_Spline_North, Bnd_Spline_South,
           Bnd_Spline_East,  Bnd_Spline_West, Bow_Spline;

  // Allocate memory for grid block.
  _Number_of_Blocks_Idir_ = 7;
  _Number_of_Blocks_Jdir_ = 6;
  allocate(_Number_of_Blocks_Idir_, _Number_of_Blocks_Jdir_);

  // Create the mesh for each block representing the complete grid.
  for (int jBlk = 0; jBlk < Number_of_Blocks_Jdir; jBlk++) {
    for (int iBlk = 0; iBlk < Number_of_Blocks_Idir; iBlk++) {

      // Preset the block indicator and turn off when required:
      block_flag = ON;

      // Preset all boundary conditions to none and reset when required:
      BCtypeN = BC_NONE;
      BCtypeS = BC_NONE;
      BCtypeE = BC_NONE;
      BCtypeW = BC_NONE;

      // Preset all stretching parameters and reset when required:
      Stretch_I = STRETCHING_FCN_LINEAR;
      Beta_I = ZERO;
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_LINEAR;
      Beta_J = ZERO;
      Tau_J = ZERO;
      Orthogonal_North = 0;
      Orthogonal_South = 0;
      Orthogonal_East = 0;
      Orthogonal_West = 0;

      // Determine the block indicator and all boundary splines and
      // reset all boundary conditions and stretching parameters for the
      // current block: 
      if (iBlk == 0 && jBlk == 0) {
	xSW = Vector2D(-0.005000,0.000000);
	xSE = Vector2D( 0.000000,0.000000);
	xNW = Vector2D(-0.005000,0.000250);
	xNE = Vector2D( 0.000000,0.000250);
	BCtypeS = Chamber_BC_Type;
	BCtypeW = Chamber_BC_Type;
      } else if (iBlk == 1 && jBlk == 0) {
	xSW = Vector2D( 0.000000,0.000000);
	xSE = Vector2D( 0.006500,0.000000);
	xNW = Vector2D( 0.000000,0.000250);
	xNE = Vector2D( 0.006500,0.000250);
	BCtypeS = Chamber_BC_Type;
      } else if (iBlk == 2 && jBlk == 0) {
	xSW = Vector2D( 0.006500,0.000000);
	xSE = Vector2D( 0.013000,0.000000);
	xNW = Vector2D( 0.006500,0.000250);
	xNE = Vector2D( 0.013000,0.000250);
	BCtypeS = Chamber_BC_Type;
	Stretch_I = STRETCHING_FCN_MAX_CLUSTERING;
	Beta_I = 1.01;
	Tau_I = ZERO;
      } else if (iBlk == 3 && jBlk == 0) {
	xSW = Vector2D( 0.013000,0.000000);
	xSE = Vector2D( 0.014000,0.000000);
	xNW = Vector2D( 0.013000,0.000250);
	xNE = Vector2D( 0.014000,0.000250);
	BCtypeS = Chamber_BC_Type;
	Stretch_I = STRETCHING_FCN_MINMAX_CLUSTERING;
	Beta_I = 1.05;
	Tau_I = ZERO;
      } else if (iBlk == 4 && jBlk == 0) {
	xSW = Vector2D( 0.014000,0.000000);
	xSE = Vector2D( 0.014250,0.000000);
	xNW = Vector2D( 0.014000,0.000250);
	xNE = Vector2D( 0.014250,0.000250);
	BCtypeN = Chamber_BC_Type;
	BCtypeS = Chamber_BC_Type;
      } else if (iBlk == 5 && jBlk == 0) {
	xSW = Vector2D( 0.014250,0.000000);
	xSE = Vector2D( 0.040000,0.000000);
	xNW = Vector2D( 0.014250,0.000250);
	xNE = Vector2D( 0.040000,0.000250);
	BCtypeS = Chamber_BC_Type;
	Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_I = 1.001;
	Tau_I = ZERO;
      } else if (iBlk == 6 && jBlk == 0) {
	xSW = Vector2D( 0.040000,0.000000);
	xSE = Vector2D( 0.100000,0.000000);
	xNW = Vector2D( 0.040000,0.000250);
	xNE = Vector2D( 0.100000,0.000250);
	BCtypeS = Chamber_BC_Type;
	BCtypeE = BC_CONSTANT_EXTRAPOLATION;

      } else if (iBlk == 0 && jBlk == 1) {
	xSW = Vector2D(-0.005000,0.000250);
	xSE = Vector2D( 0.000000,0.000250);
	xNW = Vector2D(-0.005000,0.000500);
	xNE = Vector2D( 0.000000,0.000500);
	BCtypeW = Chamber_BC_Type;
      } else if (iBlk == 1 && jBlk == 1) {
	xSW = Vector2D( 0.000000,0.000250);
	xSE = Vector2D( 0.006500,0.000250);
	xNW = Vector2D( 0.000000,0.000500);
	xNE = Vector2D( 0.006500,0.000500);
	BCtypeN = Chamber_BC_Type;
      } else if (iBlk == 2 && jBlk == 1) {
	xSW = Vector2D( 0.006500,0.000250);
	xSE = Vector2D( 0.013000,0.000250);
	xNW = Vector2D( 0.006500,0.000500);
	xNE = Vector2D( 0.013000,0.000500);
	BCtypeN = Chamber_BC_Type;
	Stretch_I = STRETCHING_FCN_MAX_CLUSTERING;
	Beta_I = 1.01;
	Tau_I = ZERO;
      } else if (iBlk == 3 && jBlk == 1) {
	xSW = Vector2D( 0.013000,0.000250);
	xSE = Vector2D( 0.014000,0.000250);
	xNW = Vector2D( 0.013000,0.000500);
	xNE = Vector2D( 0.014000,0.000500);
	BCtypeE = Chamber_BC_Type;
	Stretch_I = STRETCHING_FCN_MINMAX_CLUSTERING;
	Beta_I = 1.05;
	Tau_I = ZERO;
      } else if (iBlk == 4 && jBlk == 1) {
	block_flag = OFF;
      } else if (iBlk == 5 && jBlk == 1) {
	xSW = Vector2D( 0.014250,0.000250);
	xSE = Vector2D( 0.040000,0.000250);
	xNW = Vector2D( 0.014250,0.000500);
	xNE = Vector2D( 0.040000,0.000500);
	BCtypeW = Chamber_BC_Type;
	Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_I = 1.001;
	Tau_I = ZERO;
      } else if (iBlk == 6 && jBlk == 1) {
	xSW = Vector2D( 0.040000,0.000250);
	xSE = Vector2D( 0.100000,0.000250);
	xNW = Vector2D( 0.040000,0.000500);
	xNE = Vector2D( 0.100000,0.000500);
	BCtypeE = BC_CONSTANT_EXTRAPOLATION;

      } else if (iBlk == 0 && jBlk == 2) {
	xSW = Vector2D(-0.005000,0.000500);
	xSE = Vector2D( 0.000000,0.000500);
	xNW = Vector2D(-0.005000,0.000925);
	xNE = Vector2D( 0.000000,0.000925);
	BCtypeE = Chamber_BC_Type;
	BCtypeW = Chamber_BC_Type;
      } else if (iBlk == 1 && jBlk == 2) {
	block_flag = OFF;
      } else if (iBlk == 2 && jBlk == 2) {
	block_flag = OFF;
      } else if (iBlk == 3 && jBlk == 2) {
	xSW = Vector2D( 0.01300,0.0005000);
	xSE = Vector2D( 0.01400,0.0005000);
	xNW = Vector2D( 0.01300,0.0009250);
	xNE = Vector2D( 0.01400,0.0009250);
	BCtypeN = Chamber_BC_Type;
	BCtypeE = Chamber_BC_Type;
	BCtypeW = Chamber_BC_Type;
	Stretch_I = STRETCHING_FCN_MINMAX_CLUSTERING;
	Beta_I = 1.05;
	Tau_I = ZERO;
      } else if (iBlk == 4 && jBlk == 2) {
	block_flag = OFF;
      } else if (iBlk == 5 && jBlk == 2) {
	xSW = Vector2D( 0.014250,0.000500);
	xSE = Vector2D( 0.040000,0.000500);
	xNW = Vector2D( 0.014250,0.000925);
	xNE = Vector2D( 0.040000,0.000925);
	BCtypeW = Chamber_BC_Type;
	Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_I = 1.001;
	Tau_I = ZERO;
      } else if (iBlk == 6 && jBlk == 2) {
	xSW = Vector2D( 0.040000,0.000500);
	xSE = Vector2D( 0.100000,0.000500);
	xNW = Vector2D( 0.040000,0.000925);
	xNE = Vector2D( 0.100000,0.000925);
	BCtypeE = BC_CONSTANT_EXTRAPOLATION;

      } else if (iBlk == 0 && jBlk == 3) {
	xSW = Vector2D(-0.005000,0.000925);
	xSE = Vector2D( 0.000000,0.000925);
	xNW = Vector2D(-0.005000,0.004500);
	xNE = Vector2D( 0.000000,0.004500);
	BCtypeE = Chamber_BC_Type;
	BCtypeW = Chamber_BC_Type;
      } else if (iBlk == 1 && jBlk == 3) {
	block_flag = OFF;
      } else if (iBlk == 2 && jBlk == 3) {
	block_flag = OFF;
      } else if (iBlk == 3 && jBlk == 3) {
	block_flag = OFF;
      } else if (iBlk == 4 && jBlk == 3) {
	block_flag = OFF;
      } else if (iBlk == 5 && jBlk == 3) {
	xSW = Vector2D( 0.014250,0.000925);
	xSE = Vector2D( 0.040000,0.000925);
	xNW = Vector2D( 0.014250,0.004500);
	xNE = Vector2D( 0.040000,0.004500);
	BCtypeW = Chamber_BC_Type;
	Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_I = 1.001;
	Tau_I = ZERO;
      } else if (iBlk == 6 && jBlk == 3) {
	xSW = Vector2D( 0.040000,0.000925);
	xSE = Vector2D( 0.100000,0.000925);
	xNW = Vector2D( 0.040000,0.004500);
	xNE = Vector2D( 0.100000,0.004500);
	BCtypeE = BC_CONSTANT_EXTRAPOLATION;

      } else if (iBlk == 0 && jBlk == 4) {
	xSW = Vector2D(-0.005000,0.004500);
	xSE = Vector2D( 0.000000,0.004500);
	xNW = Vector2D(-0.005000,0.009000);
	xNE = Vector2D( 0.000000,0.009000);
	BCtypeN = Chamber_BC_Type;
	BCtypeE = Chamber_BC_Type;
	BCtypeW = Chamber_BC_Type;
      } else if (iBlk == 1 && jBlk == 4) {
	block_flag = OFF;
      } else if (iBlk == 2 && jBlk == 4) {
	block_flag = OFF;
      } else if (iBlk == 3 && jBlk == 4) {
	block_flag = OFF;
      } else if (iBlk == 4 && jBlk == 4) {
	block_flag = OFF;
      } else if (iBlk == 5 && jBlk == 4) {
	xSW = Vector2D( 0.014250,0.004500);
	xSE = Vector2D( 0.040000,0.004500);
	xNW = Vector2D( 0.014250,0.009000);
	xNE = Vector2D( 0.040000,0.009000);
	BCtypeW = Chamber_BC_Type;
	Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_I = 1.001;
	Tau_I = ZERO;
      } else if (iBlk == 6 && jBlk == 4) {
	xSW = Vector2D( 0.040000,0.004500);
	xSE = Vector2D( 0.100000,0.004500);
	xNW = Vector2D( 0.040000,0.009000);
	xNE = Vector2D( 0.100000,0.009000);
	BCtypeE = BC_CONSTANT_EXTRAPOLATION;

      } else if (iBlk == 0 && jBlk == 5) {
	block_flag = OFF;
      } else if (iBlk == 1 && jBlk == 5) {
	block_flag = OFF;
      } else if (iBlk == 2 && jBlk == 5) {
	block_flag = OFF;
      } else if (iBlk == 3 && jBlk == 5) {
	block_flag = OFF;
      } else if (iBlk == 4 && jBlk == 5) {
	block_flag = OFF;
      } else if (iBlk == 5 && jBlk == 5) {
	xSW = Vector2D( 0.014250,0.009000);
	xSE = Vector2D( 0.040000,0.009000);
	xNW = Vector2D( 0.014250,0.070000);
	xNE = Vector2D( 0.040000,0.070000);
	BCtypeN = Chamber_BC_Type;
	BCtypeW = Chamber_BC_Type;
	Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_I = 1.001;
	Tau_I = ZERO;
	Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_J = 1.015;
	Tau_J = ZERO;
      } else if (iBlk == 6 && jBlk == 5) {
	xSW = Vector2D( 0.040000,0.009000);
	xSE = Vector2D( 0.100000,0.009000);
	xNW = Vector2D( 0.040000,0.070000);
	xNE = Vector2D( 0.100000,0.070000);
	BCtypeN = Chamber_BC_Type;
	BCtypeE = BC_CONSTANT_EXTRAPOLATION;
	Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_J = 1.015;
	Tau_J = ZERO;

      }

      if (block_flag) {

	// Create the splines defining the north, south, east, and west 
	// boundaries of the current block:
	Bnd_Spline_North.Create_Spline_Line(xNW,xNE,2);
	Bnd_Spline_South.Create_Spline_Line(xSW,xSE,2);
	Bnd_Spline_East.Create_Spline_Line(xSE,xNE,2);
	Bnd_Spline_West.Create_Spline_Line(xSW,xNW,2);

	// Set the boundary condition types for each of the boundary
	// splines:
	Bnd_Spline_North.setBCtype(BCtypeN);
	Bnd_Spline_South.setBCtype(BCtypeS);
	Bnd_Spline_East.setBCtype(BCtypeE);
	Bnd_Spline_West.setBCtype(BCtypeW);

	// Create the 2D quadrilateral grid block.
	Grid_ptr[iBlk][jBlk].Create_Quad_Block_Without_Update(Bnd_Spline_North,
							      Bnd_Spline_South,
							      Bnd_Spline_East,
							      Bnd_Spline_West,
							      Number_of_Cells_Idir,
							      Number_of_Cells_Jdir,
							      Number_of_Ghost_Cells,
							      Highest_Order_of_Reconstruction,
							      GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
							      Stretch_I,
							      Beta_I,
							      Tau_I,
							      Stretch_J,
							      Beta_J,
							      Tau_J,
							      Orthogonal_North,
							      Orthogonal_South,
							      Orthogonal_East,
							      Orthogonal_West);
	
      }

    }
  }

}


/*!
 * This routine creates a mesh to be used in conjunction with and     
 * embedded interface representing a section of a NASA rotor 37 blade 
 * at a given percent span.                                           
 *
 * \verbatim                                                         
 *  +---------------+-------------------+---------------+             
 *  |               |                   |               |             
 *  |               |                   |               |             
 *  |  Block (0,1)  |    Block (1,1)    |  Block (2,1)  |             
 *  |               |                   |               |             
 *  |               |                   |               |             
 *  +---------------<<<<<<<<BLADE>>>>>>>>---------------+             
 *  |               |                   |               |             
 *  |               |                   |               |             
 *  |  Block (0,0)  |    Block (1,0)    |  Block (2,0)  |             
 *  |               |                   |               |             
 *  |               |                   |               |             
 *  +---------------+-------------------+---------------+             
 *
 * \endverbatim
 *
 * This subroutine DOESN'T update the ghost cells or
 * the geometric properties of the grid cells.
 */
void Grid2D_Quad_MultiBlock_HO::Grid_NASA_Rotor_37_Without_Update(int &_Number_of_Blocks_Idir_,
								  int &_Number_of_Blocks_Jdir_,
								  const double &Rotor_Percent_Span,
								  const int Number_of_Cells_Idir,
								  const int Number_of_Cells_Jdir,
								  const int Number_of_Ghost_Cells,
								  const int Highest_Order_of_Reconstruction) {
  
  int error_flag;
  Spline2D_HO upperB, lowerB, upperMiddleB, lowerMiddleB, camberTrail, camberLead, camberBlade, Rotor_Spline;
  int iTrail, iLead, pos, mIndex;
  Vector2D swap, leadV, trailV, x_temp;
  double zlU, zlL, zrU, zrL, mLead, mTrail, mLeadc, mTrailc, 
         mTop, mBot, m1, m2, z, dz, A, m, zLead, zTrail, dm, swapm;
  double beta_i, beta_j;

  // NASA Rotor 37 variables.
  char NASA_Rotor_Data_Directory[128];
  int Rotor_Flow_Type;
  NASARotor37 NASA_Rotor_37;

  // Initialize NASA rotor 37 class.
  strcpy(NASA_Rotor_Data_Directory,"CFFC/data/NASA_Rotors/R37/");
  Rotor_Flow_Type = 1;//PEAK_FLOW;
  NASA_Rotor_37.init(Rotor_Flow_Type, NASA_Rotor_Data_Directory);

  // Allocate memory for grid block.
  _Number_of_Blocks_Idir_ = 3;
  _Number_of_Blocks_Jdir_ = 2;
  allocate(_Number_of_Blocks_Idir_,_Number_of_Blocks_Jdir_);

  // 'Hard-wire' stretching parameters.
  if (Number_of_Cells_Idir/2 < 10) {
    beta_i = 1.05;
  } else if (Number_of_Cells_Idir/2 < 25) {
    beta_i = 1.01;
  } else if (Number_of_Cells_Idir/2 < 50) {
    beta_i = 1.005;
  } else if (Number_of_Cells_Idir/2 < 100) {
    beta_i = 1.0025;
  } else {
    beta_i = 1.001;
  } 
  if (Number_of_Cells_Jdir < 10) {
    beta_j = 1.05;
  } else if (Number_of_Cells_Jdir < 25) {
    beta_j = 1.01;
  } else if (Number_of_Cells_Jdir < 50) {
    beta_j = 1.005;
  } else if (Number_of_Cells_Jdir < 100) {
    beta_j = 1.0025;
  } else {
    beta_j = 1.001;
  }

  // Find the index for the leading and trailing edge.
  NASA_Rotor_37.findLT(Rotor_Percent_Span,iTrail,iLead);

  // Get cross section spline.
  Rotor_Spline = NASA_Rotor_37.getBladeCS(Rotor_Percent_Span);

  // Re-organize rotor spline points such that they ascend in the 
  // clockwise direction
  for (int i = 0; i < 100; i++) {
    swap = Rotor_Spline.Xp[2*100-2-i];
    Rotor_Spline.Xp[2*100-2-i] = Rotor_Spline.Xp[i];
    Rotor_Spline.Xp[i] = swap;
  }

  // Exchange iTrail and iLead accordingly.
  iTrail = Rotor_Spline.np-1 - iTrail;
  iLead = Rotor_Spline.np-1 - iLead;

  // Save leading and trailing edge vectors.
  leadV = Rotor_Spline.Xp[iLead];
  trailV = Rotor_Spline.Xp[iTrail];

  // Get camber line splines - 40 pts.
  camberTrail = NASA_Rotor_37.getCamberTrail(Rotor_Percent_Span,40,NASA_Rotor_37.z_d); 
  camberLead = NASA_Rotor_37.getCamberLead(Rotor_Percent_Span,40,NASA_Rotor_37.z_up);
  camberBlade = NASA_Rotor_37.getCamberBlade(Rotor_Percent_Span,40);
  
  // Get upper, lower, upperMiddle, and lowerMiddle boundaries.
  upperB = NASA_Rotor_37.getCamberAndShift(Rotor_Percent_Span,40,NASA_Rotor_37.z_up,NASA_Rotor_37.z_d,PI/NASA_Rotor_37.num_blades);
  lowerB = NASA_Rotor_37.getCamberAndShift(Rotor_Percent_Span,40,NASA_Rotor_37.z_up,NASA_Rotor_37.z_d,-PI/NASA_Rotor_37.num_blades);
  upperMiddleB = NASA_Rotor_37.getCamberAndShift(Rotor_Percent_Span,40,NASA_Rotor_37.z_up,NASA_Rotor_37.z_d,PI/(2*NASA_Rotor_37.num_blades));
  lowerMiddleB = NASA_Rotor_37.getCamberAndShift(Rotor_Percent_Span,40,NASA_Rotor_37.z_up,NASA_Rotor_37.z_d,-PI/(2*NASA_Rotor_37.num_blades));

  // Locate boundary corners to ensure that boundaries meeting at the leading
  // and trailing edges are separated by approximately 45deg
   
  // Calculate slope of leading and trailing camber lines.
  mLeadc = (camberLead.Xp[camberLead.np-1].y - camberLead.Xp[0].y)/
           (camberLead.Xp[camberLead.np-1].x - camberLead.Xp[0].x);

  mTrailc = (camberTrail.Xp[camberTrail.np-1].y-camberTrail.Xp[0].y)/
            (camberTrail.Xp[camberTrail.np-1].x-camberTrail.Xp[0].x);

  // Calculate slope of tangent to cross-section at the leading and 
  // trailing edge using a central-difference
  z = Rotor_Spline.Xp[iLead].x;
  dz = (Rotor_Spline.Xp[iLead].x - Rotor_Spline.Xp[iLead+1].x)/TEN;
  mLead = (-   getY(z+2*dz,Rotor_Spline)[0].y 
	   + 8*getY(z+  dz,Rotor_Spline)[0].y
	   - 8*getY(z-  dz,Rotor_Spline)[0].y + 
               getY(z-2*dz,Rotor_Spline)[0].y)/(12*dz);

  z = Rotor_Spline.Xp[iTrail].x;
  dz = (Rotor_Spline.Xp[iTrail].x - Rotor_Spline.Xp[iTrail-1].x)/TEN;
  mTrail = (-   getY(z+2*dz,Rotor_Spline)[1].y
	    + 8*getY(z+  dz,Rotor_Spline)[1].y 
	    - 8*getY(z-  dz,Rotor_Spline)[1].y
	    +   getY(z-2*dz,Rotor_Spline)[1].y)/(12*dz);

  // Get slopes of lines approx 45deg to camber line at the LEADING EDGE.
  A = (mLead+mLeadc)/(1-mLead*mLeadc);
  mTop = (-1+sqrt(1+A*A))/A;
  mBot = (-1-sqrt(1+A*A))/A;

  // Find point on upperMiddle spline - leading edge.
  dm = 1e15;
  mIndex = -1;
  for (int i = 0; i < upperMiddleB.np; i++) {
    if (fabs(upperMiddleB.Xp[i].x-leadV.x) > TOLER) {
      m = (upperMiddleB.Xp[i].y-leadV.y)/(upperMiddleB.Xp[i].x-leadV.x);
      if (fabs(m-mTop) < dm) {
	dm = fabs(m-mTop);
	mIndex = i;
      }
    }
  }
  assert(mIndex != -1);
  zlU = upperMiddleB.Xp[mIndex].x;  

  // Find point on lowerMiddle spline - leading edge.
  dm = 1e15;
  mIndex = -1;
  for (int i = 0; i < lowerMiddleB.np; i++) {
    if (fabs(lowerMiddleB.Xp[i].x-leadV.x) > TOLER) {
      m = (lowerMiddleB.Xp[i].y-leadV.y)/(lowerMiddleB.Xp[i].x-leadV.x);
      if (fabs(m-mBot) < dm) {
	dm=fabs(m-mBot);
	mIndex=i;
      }
    }
  }
  assert(mIndex != -1);
  zlL = lowerMiddleB.Xp[mIndex].x;  

  // Get slopes of lines approx 45deg to camber line at the TRAILING EDGE.
  A = (mTrail+mTrailc)/(1-mTrail*mTrailc);
  mBot = (-1+sqrt(1+A*A))/A;
  mTop = (-1-sqrt(1+A*A))/A;
  
  // Find point on upperMiddle spline - trailing edge.
  dm = 1e15;
  mIndex = -1;
  for (int i=0; i < upperMiddleB.np; i++) {
    if (fabs(upperMiddleB.Xp[i].x-trailV.x) > TOLER) {
      m = (upperMiddleB.Xp[i].y-trailV.y)/(upperMiddleB.Xp[i].x-trailV.x);
      if (fabs(m-mTop) < dm) {
	dm=fabs(m-mTop);
	mIndex=i;
      }
    }
  }
  assert(mIndex != -1);
  zrU = upperMiddleB.Xp[mIndex].x;  

  dm = 1e15;
  mIndex = -1;
  for (int i = 0; i < upperMiddleB.np; i++) {
    if (fabs(upperMiddleB.Xp[i].x-trailV.x) > TOLER) {
      m = (upperMiddleB.Xp[i].y-trailV.y)/(upperMiddleB.Xp[i].x-trailV.x);
      if (fabs(m-mBot) < dm) {
	dm = fabs(m-mBot);
	mIndex = i;
      }
    }
  }
  assert(mIndex != -1);
  
  // Determine which to use.
  if (upperMiddleB.Xp[mIndex].x>zrU) {
    zrU = upperMiddleB.Xp[mIndex].x;  
    swapm = mBot;
    mBot = mTop;
    mTop = swapm;
  }

  // Find point on lowerMiddle spline - trailing edge.
  dm = 1e15;
  mIndex = -1;
  for (int i=0; i < lowerMiddleB.np; i++) {
    if (fabs(lowerMiddleB.Xp[i].x-trailV.x) > TOLER) {
      m = (lowerMiddleB.Xp[i].y-trailV.y)/(lowerMiddleB.Xp[i].x-trailV.x);
      if (fabs(m-mBot) < dm) {
	dm = fabs(m-mBot);
	mIndex = i;
      }
    }
  }
  assert(mIndex != -1);
  zrL = lowerMiddleB.Xp[mIndex].x;

  /////////////////
  // BLOCK (1,1) //
  /////////////////
  // North spline.
  Grid_ptr[1][1].BndNorthSpline.allocate(70);
  Grid_ptr[1][1].BndNorthSpline.settype(SPLINE2D_QUINTIC);
  zLead = (zlL+zlU)/2;
  zTrail = (zrL+zrU)/2;
  dz = (zTrail-zLead)/69;
  for (int i = 0; i < 70; i++) {
    Grid_ptr[1][1].BndNorthSpline.Xp[i] = getY(zLead+i*dz,upperB)[0];
    Grid_ptr[1][1].BndNorthSpline.bc[i] = BC_NONE;
    if (i == 0 || i == Grid_ptr[1][1].BndNorthSpline.np-1) {
      Grid_ptr[1][1].BndNorthSpline.tp[i] = SPLINE2D_POINT_SHARP_CORNER;
    } else {
      Grid_ptr[1][1].BndNorthSpline.tp[i] = SPLINE2D_POINT_NORMAL;
    }
  }
  Grid_ptr[1][1].BndNorthSpline.pathlength();
  // South spline.
  Grid_ptr[1][1].BndSouthSpline.allocate(iTrail-iLead+1);
  Grid_ptr[1][1].BndSouthSpline.settype(SPLINE2D_QUINTIC);
  for (int i = iLead; i <= iTrail; i++)
    Grid_ptr[1][1].BndSouthSpline.Xp[i-iLead] = Rotor_Spline.Xp[i];
  Grid_ptr[1][1].BndSouthSpline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[1][1].BndSouthSpline.tp[Grid_ptr[1][1].BndSouthSpline.np-1] = SPLINE2D_POINT_SHARP_CORNER;
  for (int i = 1; i < Grid_ptr[1][1].BndSouthSpline.np-1; i++)
    Grid_ptr[1][1].BndSouthSpline.tp[i]=SPLINE2D_POINT_NORMAL;
  Grid_ptr[1][1].BndSouthSpline.pathlength();
  // East spline.
  Grid_ptr[1][1].BndEastSpline.allocate(5);
  Grid_ptr[1][1].BndEastSpline.settype(SPLINE2D_CUBIC);
  Grid_ptr[1][1].BndEastSpline.Xp[0] = Rotor_Spline.Xp[iTrail];
  x_temp = getY(zrU,upperMiddleB)[0];
  Grid_ptr[1][1].BndEastSpline.Xp[1] = Rotor_Spline.Xp[iTrail]+HALF*(x_temp-Rotor_Spline.Xp[iTrail]);
  Grid_ptr[1][1].BndEastSpline.Xp[2] = x_temp;
  Grid_ptr[1][1].BndEastSpline.Xp[3] = x_temp + HALF*(Grid_ptr[1][1].BndNorthSpline.Xp[Grid_ptr[1][1].BndNorthSpline.np-1] - x_temp);
  Grid_ptr[1][1].BndEastSpline.Xp[4] = Grid_ptr[1][1].BndNorthSpline.Xp[Grid_ptr[1][1].BndNorthSpline.np-1];
  Grid_ptr[1][1].BndEastSpline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[1][1].BndEastSpline.tp[1] = SPLINE2D_POINT_NORMAL;
  Grid_ptr[1][1].BndEastSpline.tp[2] = SPLINE2D_POINT_NORMAL;
  Grid_ptr[1][1].BndEastSpline.tp[3] = SPLINE2D_POINT_NORMAL;
  Grid_ptr[1][1].BndEastSpline.tp[4] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[1][1].BndEastSpline.pathlength();
  // West spline.
  Grid_ptr[1][1].BndWestSpline.allocate(6);
  Grid_ptr[1][1].BndWestSpline.settype(SPLINE2D_CUBIC);
  Grid_ptr[1][1].BndWestSpline.Xp[0] = Rotor_Spline.Xp[iLead];
  x_temp = getY(zlU,upperMiddleB)[0];
  Grid_ptr[1][1].BndWestSpline.Xp[1] = Rotor_Spline.Xp[iLead]+0.3333*(x_temp-Rotor_Spline.Xp[iLead]);
  Grid_ptr[1][1].BndWestSpline.Xp[2] = HALF*(Rotor_Spline.Xp[iLead]+HALF*(x_temp-Rotor_Spline.Xp[iLead])+
			               x_temp+0.20*(Grid_ptr[1][1].BndNorthSpline.Xp[0]-x_temp));
  Grid_ptr[1][1].BndWestSpline.Xp[3] = HALF*(Rotor_Spline.Xp[iLead]+0.6666*(x_temp-Rotor_Spline.Xp[iLead])+
			               x_temp+0.3333*(Grid_ptr[1][1].BndNorthSpline.Xp[0]-x_temp));
  Grid_ptr[1][1].BndWestSpline.Xp[4] = x_temp+HALF*(Grid_ptr[1][1].BndNorthSpline.Xp[0]-x_temp);
  Grid_ptr[1][1].BndWestSpline.Xp[5] = Grid_ptr[1][1].BndNorthSpline.Xp[0];
  Grid_ptr[1][1].BndWestSpline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[1][1].BndWestSpline.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[1][1].BndWestSpline.tp[2] = SPLINE2D_POINT_NORMAL;
  Grid_ptr[1][1].BndWestSpline.tp[3] = SPLINE2D_POINT_NORMAL;
  Grid_ptr[1][1].BndWestSpline.tp[4] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[1][1].BndWestSpline.tp[5] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[1][1].BndWestSpline.pathlength();
  // Create quad block (1,1).
  Grid_ptr[1][1].Create_Quad_Block_Without_Update(Grid_ptr[1][1].BndNorthSpline,
						  Grid_ptr[1][1].BndSouthSpline,
						  Grid_ptr[1][1].BndEastSpline,
						  Grid_ptr[1][1].BndWestSpline,
						  Number_of_Cells_Idir,
						  Number_of_Cells_Jdir,
						  Number_of_Ghost_Cells,
						  Highest_Order_of_Reconstruction,
						  GRID2D_QUAD_BLOCK_INIT_PROCEDURE_TRANS_FINITE_XY,
						  STRETCHING_FCN_MINMAX_CLUSTERING,
						  ONE + (beta_i-ONE)/ONE,
						  0,
						  STRETCHING_FCN_MIN_CLUSTERING,
						  beta_j,
						  0,
						  ORTHOGONAL,
						  ORTHOGONAL,
						  ORTHOGONAL,
						  ORTHOGONAL);
  // Smooth quad block (1,1).
  Grid_ptr[1][1].Smooth_Quad_Block(
				   min(250,2*max(Number_of_Cells_Idir,
						 Number_of_Cells_Jdir)));

  /////////////////
  // BLOCK (1,0) //
  /////////////////
  // North spline.
  Grid_ptr[1][0].BndNorthSpline.allocate(Rotor_Spline.np-iTrail+iLead);
  Grid_ptr[1][0].BndNorthSpline.settype(SPLINE2D_QUINTIC);
  for (int i = iLead; i >= 0; i--)
    Grid_ptr[1][0].BndNorthSpline.Xp[iLead-i] = Rotor_Spline.Xp[i];
  for (int i=Rotor_Spline.np-2; i >= iTrail; i--)
    Grid_ptr[1][0].BndNorthSpline.Xp[iLead+1+(Rotor_Spline.np-2)-i] = Rotor_Spline.Xp[i];
  for (int i = 0; i < Grid_ptr[1][0].BndNorthSpline.np; i++) {
    Grid_ptr[1][0].BndNorthSpline.bc[i] = BC_NONE;
    if (i == 0 || i == Grid_ptr[1][0].BndNorthSpline.np-1) {
      Grid_ptr[1][0].BndNorthSpline.tp[i] = SPLINE2D_POINT_SHARP_CORNER;
    } else {
      Grid_ptr[1][0].BndNorthSpline.tp[i] = SPLINE2D_POINT_NORMAL;
    }
  }
  Grid_ptr[1][0].BndNorthSpline.pathlength();  
  // South spline.
  Grid_ptr[1][0].BndSouthSpline.allocate(70);
  Grid_ptr[1][0].BndSouthSpline.settype(SPLINE2D_QUINTIC);
  zLead = (zlL+zlU)/2;;
  zTrail = (zrL+zrU)/2;;
  dz = (zTrail-zLead)/69;
  for (int i = 0; i < 70; i++) 
    Grid_ptr[1][0].BndSouthSpline.Xp[i] = getY(zLead+i*dz,lowerB)[0];
  for (int i = 0; i < Grid_ptr[1][0].BndSouthSpline.np; i++) {
    Grid_ptr[1][0].BndSouthSpline.bc[i] = BC_NONE;
    if (i == 0 || i == Grid_ptr[1][0].BndSouthSpline.np-1) {
      Grid_ptr[1][0].BndSouthSpline.tp[i] = SPLINE2D_POINT_SHARP_CORNER;
    } else {
      Grid_ptr[1][0].BndSouthSpline.tp[i] = SPLINE2D_POINT_NORMAL;
    }
  }
  Grid_ptr[1][0].BndSouthSpline.pathlength();
  // West spline.
  Grid_ptr[1][0].BndWestSpline.allocate(5);
  Grid_ptr[1][0].BndWestSpline.settype(SPLINE2D_CUBIC);
  Grid_ptr[1][0].BndWestSpline.Xp[0] = Grid_ptr[1][0].BndSouthSpline.Xp[0];
  x_temp = getY(zlL,lowerMiddleB)[0];
  Grid_ptr[1][0].BndWestSpline.Xp[1] = Grid_ptr[1][0].BndSouthSpline.Xp[0]+
                                       HALF*(x_temp-Grid_ptr[1][0].BndSouthSpline.Xp[0]);
  Grid_ptr[1][0].BndWestSpline.Xp[2] = x_temp;
  Grid_ptr[1][0].BndWestSpline.Xp[3] = x_temp+HALF*(Rotor_Spline.Xp[iLead]-x_temp);
  Grid_ptr[1][0].BndWestSpline.Xp[4] = Rotor_Spline.Xp[iLead];
  for (int i = 0; i < 5; i++) Grid_ptr[1][0].BndWestSpline.bc[i] = BC_NONE;
  Grid_ptr[1][0].BndWestSpline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[1][0].BndWestSpline.tp[1] = SPLINE2D_POINT_NORMAL;
  Grid_ptr[1][0].BndWestSpline.tp[2] = SPLINE2D_POINT_NORMAL;
  Grid_ptr[1][0].BndWestSpline.tp[3] = SPLINE2D_POINT_NORMAL;
  Grid_ptr[1][0].BndWestSpline.tp[4] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[1][0].BndWestSpline.pathlength();
  // East spline.
  Grid_ptr[1][0].BndEastSpline.allocate(6);
  Grid_ptr[1][0].BndEastSpline.settype(SPLINE2D_CUBIC);
  Grid_ptr[1][0].BndEastSpline.Xp[0] = Grid_ptr[1][0].BndSouthSpline.Xp[Grid_ptr[1][0].BndSouthSpline.np-1];
  x_temp = getY(zrL,lowerMiddleB)[0];
  Grid_ptr[1][0].BndEastSpline.Xp[1] = x_temp + HALF*(Grid_ptr[1][0].BndSouthSpline.Xp[Grid_ptr[1][0].BndSouthSpline.np-1]-x_temp);
  Grid_ptr[1][0].BndEastSpline.Xp[2] = HALF*(Rotor_Spline.Xp[iTrail]+0.6666*(x_temp-Rotor_Spline.Xp[iTrail])+
                                       x_temp+0.3333*(Grid_ptr[1][0].BndSouthSpline.Xp[Grid_ptr[1][0].BndSouthSpline.np-1]-x_temp));
  Grid_ptr[1][0].BndEastSpline.Xp[3] = 0.95*(Rotor_Spline.Xp[iTrail]+HALF*(x_temp-Rotor_Spline.Xp[iTrail]))+
                                       0.05*(x_temp+0.80*(Grid_ptr[1][0].BndSouthSpline.Xp[Grid_ptr[1][0].BndSouthSpline.np-1]-x_temp));
  Grid_ptr[1][0].BndEastSpline.Xp[4] = Rotor_Spline.Xp[iTrail]+0.3333*(x_temp-Rotor_Spline.Xp[iTrail]);
  Grid_ptr[1][0].BndEastSpline.Xp[5] = Rotor_Spline.Xp[iTrail];
  for (int i = 0; i < 6; i++) Grid_ptr[1][0].BndEastSpline.bc[i] = BC_NONE;
  Grid_ptr[1][0].BndEastSpline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[1][0].BndEastSpline.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[1][0].BndEastSpline.tp[2] = SPLINE2D_POINT_NORMAL;
  Grid_ptr[1][0].BndEastSpline.tp[3] = SPLINE2D_POINT_NORMAL;
  Grid_ptr[1][0].BndEastSpline.tp[4] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[1][0].BndEastSpline.tp[5] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[1][0].BndEastSpline.pathlength();
  // Create quad block (1,0).
  Grid_ptr[1][0].Create_Quad_Block_Without_Update(Grid_ptr[1][0].BndNorthSpline,
						  Grid_ptr[1][0].BndSouthSpline,
						  Grid_ptr[1][0].BndEastSpline,
						  Grid_ptr[1][0].BndWestSpline,
						  Number_of_Cells_Idir,
						  Number_of_Cells_Jdir,
						  Number_of_Ghost_Cells,
						  Highest_Order_of_Reconstruction,
						  GRID2D_QUAD_BLOCK_INIT_PROCEDURE_TRANS_FINITE_XY,
						  STRETCHING_FCN_MINMAX_CLUSTERING,
						  ONE + (beta_i-ONE)/ONE,
						  0,
						  STRETCHING_FCN_MAX_CLUSTERING,
						  beta_j,
						  0,
						  ORTHOGONAL,
						  ORTHOGONAL,
						  ORTHOGONAL,
						  ORTHOGONAL);
  // Smooth quad block (1,0).
  Grid_ptr[1][0].Smooth_Quad_Block(
				   min(250,2*max(Number_of_Cells_Idir,
						 Number_of_Cells_Jdir)));

  /////////////////
  // BLOCK (0,0) //
  /////////////////
  // North spline.
  Grid_ptr[0][0].BndNorthSpline.allocate(2);
  Grid_ptr[0][0].BndNorthSpline.settype(SPLINE2D_LINEAR);
  Grid_ptr[0][0].BndNorthSpline.Xp[0] = camberLead.Xp[0];
  Grid_ptr[0][0].BndNorthSpline.Xp[1] = Rotor_Spline.Xp[iLead];
  Grid_ptr[0][0].BndNorthSpline.bc[0] = BC_NONE;
  Grid_ptr[0][0].BndNorthSpline.bc[1] = BC_NONE;
  Grid_ptr[0][0].BndNorthSpline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[0][0].BndNorthSpline.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[0][0].BndNorthSpline.pathlength();
  // South spline.
  Grid_ptr[0][0].BndSouthSpline.allocate(30);
  Grid_ptr[0][0].BndSouthSpline.settype(SPLINE2D_QUINTIC);
  zLead = lowerB.Xp[0].x;
  zTrail = Grid_ptr[1][0].BndSouthSpline.Xp[0].x;
  dz = (zTrail-zLead)/29;
  for (int i = 0; i < 30; i++)
    Grid_ptr[0][0].BndSouthSpline.Xp[i] = getY(zLead+i*dz, lowerB)[0];
  Grid_ptr[0][0].BndSouthSpline.Xp[0] = lowerB.Xp[0];
  for (int i = 0; i < Grid_ptr[0][0].BndSouthSpline.np; i++) {
    Grid_ptr[0][0].BndSouthSpline.bc[i] = BC_NONE;
    if (i == 0 || i == Grid_ptr[0][0].BndSouthSpline.np-1) {
      Grid_ptr[0][0].BndSouthSpline.tp[i] = SPLINE2D_POINT_SHARP_CORNER;
    } else {
      Grid_ptr[0][0].BndSouthSpline.tp[i] = SPLINE2D_POINT_NORMAL;
    }
  }
  Grid_ptr[0][0].BndSouthSpline.pathlength();
  // East spline.
  Grid_ptr[0][0].BndEastSpline = Grid_ptr[1][0].BndWestSpline;
  // West spline.
  Grid_ptr[0][0].BndWestSpline.allocate(2);
  Grid_ptr[0][0].BndWestSpline.settype(SPLINE2D_LINEAR);
  Grid_ptr[0][0].BndWestSpline.Xp[0] = Grid_ptr[0][0].BndSouthSpline.Xp[0];
  Grid_ptr[0][0].BndWestSpline.Xp[1] = Grid_ptr[0][0].BndNorthSpline.Xp[0];
  Grid_ptr[0][0].BndWestSpline.bc[0] = BC_FIXED;
  Grid_ptr[0][0].BndWestSpline.bc[1] = BC_FIXED;
  Grid_ptr[0][0].BndWestSpline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[0][0].BndWestSpline.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[0][0].BndWestSpline.pathlength();
  // Create quad block (0,0).
  Grid_ptr[0][0].Create_Quad_Block_Without_Update(Grid_ptr[0][0].BndNorthSpline,
						  Grid_ptr[0][0].BndSouthSpline,
						  Grid_ptr[0][0].BndEastSpline,
						  Grid_ptr[0][0].BndWestSpline,
						  Number_of_Cells_Idir,
						  Number_of_Cells_Jdir,
						  Number_of_Ghost_Cells,
						  Highest_Order_of_Reconstruction,
						  GRID2D_QUAD_BLOCK_INIT_PROCEDURE_EAST_WEST,
						  STRETCHING_FCN_MAX_CLUSTERING,
						  beta_i,
						  0,
						  STRETCHING_FCN_MAX_CLUSTERING,
						  beta_j,
						  0,
						  ORTHOGONAL,
						  ORTHOGONAL,
						  ORTHOGONAL,
						  NOT_ORTHOGONAL);
  // Smooth quad block (0,0).
  Grid_ptr[0][0].Smooth_Quad_Block(min(250,2*max(Number_of_Cells_Idir,
						 Number_of_Cells_Jdir)));
  
  /////////////////
  // BLOCK (2,0) //
  /////////////////
  // North spline.
  Grid_ptr[2][0].BndNorthSpline.allocate(2);
  Grid_ptr[2][0].BndNorthSpline.settype(SPLINE2D_LINEAR);
  Grid_ptr[2][0].BndNorthSpline.Xp[0] = Rotor_Spline.Xp[iTrail];
  Grid_ptr[2][0].BndNorthSpline.Xp[1] = camberTrail.Xp[camberTrail.np-1];
  Grid_ptr[2][0].BndNorthSpline.bc[0] = BC_NONE;
  Grid_ptr[2][0].BndNorthSpline.bc[1] = BC_NONE;
  Grid_ptr[2][0].BndNorthSpline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[2][0].BndNorthSpline.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[2][0].BndNorthSpline.pathlength();
  // South spline.
  Grid_ptr[2][0].BndSouthSpline.allocate(30);
  Grid_ptr[2][0].BndSouthSpline.settype(SPLINE2D_QUINTIC);
  zLead = (zrL+zrU)/2;;
  zTrail = lowerB.Xp[lowerB.np-1].x;
  dz = (zTrail-zLead)/29;
  for (int i = 0; i < 30; i++) {
    Grid_ptr[2][0].BndSouthSpline.Xp[i].x = zLead+i*dz;
    if(Grid_ptr[2][0].BndSouthSpline.Xp[i].x>zTrail) Grid_ptr[2][0].BndSouthSpline.Xp[i].x = zTrail;
    Grid_ptr[2][0].BndSouthSpline.Xp[i].y = getY(Grid_ptr[2][0].BndSouthSpline.Xp[i].x,lowerB)[0].y;
  }
  Grid_ptr[2][0].BndSouthSpline.Xp[Grid_ptr[2][0].BndSouthSpline.np-1] = lowerB.Xp[lowerB.np-1];
  for (int i = 0; i < Grid_ptr[2][0].BndSouthSpline.np; i++) {
    Grid_ptr[2][0].BndSouthSpline.tp[i] = BC_NONE;
    if (i == 0 || i == Grid_ptr[2][0].BndSouthSpline.np-1) {
      Grid_ptr[2][0].BndSouthSpline.tp[i] = SPLINE2D_POINT_SHARP_CORNER;
    } else {
      Grid_ptr[2][0].BndSouthSpline.tp[i] = SPLINE2D_POINT_NORMAL;
    }
  }
  Grid_ptr[2][0].BndSouthSpline.pathlength();
  // East spline.
  Grid_ptr[2][0].BndEastSpline.allocate(2);
  Grid_ptr[2][0].BndEastSpline.settype(SPLINE2D_LINEAR);
  Grid_ptr[2][0].BndEastSpline.Xp[0] = Grid_ptr[2][0].BndSouthSpline.Xp[Grid_ptr[2][0].BndSouthSpline.np-1];
  Grid_ptr[2][0].BndEastSpline.Xp[1] = Grid_ptr[2][0].BndNorthSpline.Xp[Grid_ptr[2][0].BndNorthSpline.np-1];
  Grid_ptr[2][0].BndEastSpline.bc[0] = BC_CHARACTERISTIC;
  Grid_ptr[2][0].BndEastSpline.bc[1] = BC_CHARACTERISTIC;
  Grid_ptr[2][0].BndEastSpline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[2][0].BndEastSpline.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[2][0].BndEastSpline.pathlength();
  // West spline.
  Grid_ptr[2][0].BndWestSpline = Grid_ptr[1][0].BndEastSpline; 
  // Create quad block (2,0).
  Grid_ptr[2][0].Create_Quad_Block_Without_Update(Grid_ptr[2][0].BndNorthSpline,
						  Grid_ptr[2][0].BndSouthSpline,
						  Grid_ptr[2][0].BndEastSpline,
						  Grid_ptr[2][0].BndWestSpline,
						  Number_of_Cells_Idir,
						  Number_of_Cells_Jdir,
						  Number_of_Ghost_Cells,
						  Highest_Order_of_Reconstruction,
						  GRID2D_QUAD_BLOCK_INIT_PROCEDURE_EAST_WEST,
						  STRETCHING_FCN_MIN_CLUSTERING,
						  beta_i,
						  0,
						  STRETCHING_FCN_MAX_CLUSTERING,
						  beta_j,
						  0,
						  ORTHOGONAL,//NOT_
						  ORTHOGONAL,
						  NOT_ORTHOGONAL,
						  ORTHOGONAL);
  // Smooth quad block (2,0).
  Grid_ptr[2][0].Smooth_Quad_Block(min(250,2*max(Number_of_Cells_Idir,
						 Number_of_Cells_Jdir)));
  
  /////////////////
  // BLOCK (0,1) //
  /////////////////
  // North spline.
  Grid_ptr[0][1].BndNorthSpline.allocate(30);
  Grid_ptr[0][1].BndNorthSpline.settype(SPLINE2D_QUINTIC);
  zLead = upperB.Xp[0].x;
  zTrail = (zlL+zlU)/2;
  dz = (zTrail-zLead)/29;
  for (int i = 0; i < 30; i++)
    Grid_ptr[0][1].BndNorthSpline.Xp[i] = getY(zLead+i*dz,upperB)[0];
  Grid_ptr[0][1].BndNorthSpline.Xp[0] = upperB.Xp[0];
  for (int i = 0; i < Grid_ptr[0][1].BndNorthSpline.np; i++) {
    Grid_ptr[0][1].BndNorthSpline.bc[i] = BC_NONE;
    if (i == 0 || i == Grid_ptr[0][1].BndNorthSpline.np-1) {
      Grid_ptr[0][1].BndNorthSpline.tp[i] = SPLINE2D_POINT_SHARP_CORNER;
    } else {
      Grid_ptr[0][1].BndNorthSpline.tp[i] = SPLINE2D_POINT_NORMAL;
    }
  }
  Grid_ptr[0][1].BndNorthSpline.pathlength();
  // South spline.
  Grid_ptr[0][1].BndSouthSpline = Grid_ptr[0][0].BndNorthSpline;
  // East spline.
  Grid_ptr[0][1].BndEastSpline = Grid_ptr[1][1].BndWestSpline;
  // West spline.
  Grid_ptr[0][1].BndWestSpline.allocate(2);
  Grid_ptr[0][1].BndWestSpline.settype(SPLINE2D_LINEAR);
  Grid_ptr[0][1].BndWestSpline.Xp[0] = Grid_ptr[0][1].BndSouthSpline.Xp[0];
  Grid_ptr[0][1].BndWestSpline.Xp[1] = Grid_ptr[0][1].BndNorthSpline.Xp[0];
  Grid_ptr[0][1].BndWestSpline.bc[0] = BC_FIXED;
  Grid_ptr[0][1].BndWestSpline.bc[1] = BC_FIXED;
  Grid_ptr[0][1].BndWestSpline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[0][1].BndWestSpline.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[0][1].BndWestSpline.pathlength();
  // Create quad block (0,1).
  Grid_ptr[0][1].Create_Quad_Block_Without_Update(Grid_ptr[0][1].BndNorthSpline,
						  Grid_ptr[0][1].BndSouthSpline,
						  Grid_ptr[0][1].BndEastSpline,
						  Grid_ptr[0][1].BndWestSpline,
						  Number_of_Cells_Idir,
						  Number_of_Cells_Jdir,
						  Number_of_Ghost_Cells,
						  Highest_Order_of_Reconstruction,
						  GRID2D_QUAD_BLOCK_INIT_PROCEDURE_EAST_WEST,
						  STRETCHING_FCN_MAX_CLUSTERING,
						  beta_i,
						  0,
						  STRETCHING_FCN_MIN_CLUSTERING,
						  beta_j,
						  0,
						  ORTHOGONAL,
						  ORTHOGONAL,
						  ORTHOGONAL,
						  NOT_ORTHOGONAL);
  // Smooth quad block (0,1).
  Grid_ptr[0][1].Smooth_Quad_Block(min(250,2*max(Number_of_Cells_Idir,
						 Number_of_Cells_Jdir)));
  
  /////////////////
  // BLOCK (2,1) //
  /////////////////
  // North spline.
  Grid_ptr[2][1].BndNorthSpline.allocate(30);
  Grid_ptr[2][1].BndNorthSpline.settype(SPLINE2D_QUINTIC);
  zLead = (zrL+zrU)/2;
  zTrail = upperB.Xp[upperB.np-1].x;
  dz = (zTrail-zLead)/29;
  for (int i = 0; i < 30; i++) {
    Grid_ptr[2][1].BndNorthSpline.Xp[i].x = zLead+i*dz;
    if (Grid_ptr[2][1].BndNorthSpline.Xp[i].x > zTrail)
      Grid_ptr[2][1].BndNorthSpline.Xp[i].x = zTrail;
    Grid_ptr[2][1].BndNorthSpline.Xp[i].y = getY(Grid_ptr[2][1].BndNorthSpline.Xp[i].x,upperB)[0].y;
  }
  Grid_ptr[2][1].BndNorthSpline.Xp[Grid_ptr[2][1].BndNorthSpline.np-1] = upperB.Xp[upperB.np-1];
  for (int i = 0; i < Grid_ptr[2][1].BndNorthSpline.np; i++) {
    Grid_ptr[2][1].BndNorthSpline.bc[i] = BC_NONE;
    if (i == 0 || i == Grid_ptr[2][1].BndNorthSpline.np-1) {
      Grid_ptr[2][1].BndNorthSpline.tp[i] = SPLINE2D_POINT_SHARP_CORNER;
    } else {
      Grid_ptr[2][1].BndNorthSpline.tp[i] = SPLINE2D_POINT_NORMAL;
    }
  }
  Grid_ptr[2][1].BndNorthSpline.pathlength();
  // South spline.
  Grid_ptr[2][1].BndSouthSpline = Grid_ptr[2][0].BndNorthSpline;
  // East spline.
  Grid_ptr[2][1].BndEastSpline.allocate(2);
  Grid_ptr[2][1].BndEastSpline.settype(SPLINE2D_LINEAR);
  Grid_ptr[2][1].BndEastSpline.Xp[0] = Grid_ptr[2][1].BndSouthSpline.Xp[Grid_ptr[2][1].BndSouthSpline.np-1];
  Grid_ptr[2][1].BndEastSpline.Xp[1] = Grid_ptr[2][1].BndNorthSpline.Xp[Grid_ptr[2][1].BndNorthSpline.np-1];
  Grid_ptr[2][1].BndEastSpline.bc[0] = BC_CHARACTERISTIC;
  Grid_ptr[2][1].BndEastSpline.bc[1] = BC_CHARACTERISTIC;
  Grid_ptr[2][1].BndEastSpline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[2][1].BndEastSpline.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[2][1].BndEastSpline.pathlength();
  // West spline.
  Grid_ptr[2][1].BndWestSpline = Grid_ptr[1][1].BndEastSpline;
  // Create quad block (2,1).
  Grid_ptr[2][1].Create_Quad_Block_Without_Update(Grid_ptr[2][1].BndNorthSpline,
						  Grid_ptr[2][1].BndSouthSpline,
						  Grid_ptr[2][1].BndEastSpline,
						  Grid_ptr[2][1].BndWestSpline,
						  Number_of_Cells_Idir,
						  Number_of_Cells_Jdir,
						  Number_of_Ghost_Cells,
						  Highest_Order_of_Reconstruction,
						  GRID2D_QUAD_BLOCK_INIT_PROCEDURE_EAST_WEST,
						  STRETCHING_FCN_MIN_CLUSTERING,
						  beta_i,
						  0,
						  STRETCHING_FCN_MIN_CLUSTERING,
						  beta_j,
						  0,
						  ORTHOGONAL,
						  ORTHOGONAL,
						  NOT_ORTHOGONAL,
						  ORTHOGONAL);
  // Smooth quad block (2,1).
  Grid_ptr[2][1].Smooth_Quad_Block(min(250,2*max(Number_of_Cells_Idir,
						 Number_of_Cells_Jdir)));
  
}

/*!
 * This routine creates a mesh to be used in conjunction with and     
 * embedded interface representing a section of a NASA rotor 67 blade 
 * at a given percent span.                                           
 *
 * \verbatim                                                         
 *  +---------------+-------------------+---------------+             
 *  |               |                   |               |             
 *  |               |                   |               |             
 *  |  Block (0,1)  |    Block (1,1)    |  Block (2,1)  |             
 *  |               |                   |               |             
 *  |               |                   |               |             
 *  +---------------<<<<<<<<BLADE>>>>>>>>---------------+             
 *  |               |                   |               |             
 *  |               |                   |               |             
 *  |  Block (0,0)  |    Block (1,0)    |  Block (2,0)  |             
 *  |               |                   |               |             
 *  |               |                   |               |             
 *  +---------------+-------------------+---------------+             
 * \endverbatim
 *
 * This subroutine DOESN'T update the ghost cells or
 * the geometric properties of the grid cells.
 */
void Grid2D_Quad_MultiBlock_HO::Grid_NASA_Rotor_67_Without_Update(int &_Number_of_Blocks_Idir_,
								  int &_Number_of_Blocks_Jdir_,
								  const double &Rotor_Percent_Span,
								  const int Number_of_Cells_Idir,
								  const int Number_of_Cells_Jdir,
								  const int Number_of_Ghost_Cells,
								  const int Highest_Order_of_Reconstruction) {

  int error_flag;
  Spline2D_HO upperB, lowerB, upperMiddleB, lowerMiddleB, camberTrail, camberLead, camberBlade, Rotor_Spline;
  int iTrail, iLead, pos, mIndex;
  Vector2D swap, leadV, trailV, x_temp;
  double zlU, zlL, zrU, zrL, mLead, mTrail, mLeadc, mTrailc, 
         mTop, mBot, m1, m2, z, dz, A, m, zLead, zTrail, dm, swapm;
  Spline2D_HO BndNorthSpline, tempBndNorthSpline, BndSouthSpline, tempBndSouthSpline, BndEastSpline, BndWestSpline;
  int Stretch_I, Stretch_J,
      Orthogonal_North, Orthogonal_South,
      Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D xNW, xNE, xSE, xSW;
  Spline2D_HO Bnd_Spline_North, Bnd_Spline_South,
           Bnd_Spline_East,  Bnd_Spline_West, Bow_Spline;

  // NASA Rotor 67 variables.
  char NASA_Rotor_Data_Directory[128];
  int Rotor_Flow_Type;
  NASARotor67 NASA_Rotor_67;

  // Initialize NASA rotor 67 class.
  strcpy(NASA_Rotor_Data_Directory,"CFFC/data/NASA_Rotors/R67/");
  Rotor_Flow_Type = 1;//PEAK_FLOW;
  NASA_Rotor_67.init(Rotor_Flow_Type,NASA_Rotor_Data_Directory);

  // Allocate memory for grid block.
  _Number_of_Blocks_Idir_ = 3;
  _Number_of_Blocks_Jdir_ = 2;
  allocate(_Number_of_Blocks_Idir_,_Number_of_Blocks_Jdir_);

  // Find the index for the leading and trailing edge.
  NASA_Rotor_67.findLT(Rotor_Percent_Span,iTrail,iLead);

  // Get cross section spline.
  Rotor_Spline = NASA_Rotor_67.getBladeCS(Rotor_Percent_Span);

  // Re-organize rotor spline points such that they ascend in the 
  // clockwise direction
  for (int i = 0; i < 70; i++ ){
    swap = Rotor_Spline.Xp[2*70-2-i];
    Rotor_Spline.Xp[2*70-2-i] = Rotor_Spline.Xp[i];
    Rotor_Spline.Xp[i] = swap;
  }

  // Exchange iTrail and iLead accordingly.
  iTrail = Rotor_Spline.np-1 - iTrail;
  iLead = Rotor_Spline.np-1 - iLead;

  // Save leading and trailing edge vectors.
  leadV = Rotor_Spline.Xp[iLead];
  trailV = Rotor_Spline.Xp[iTrail];

  // Get camber line splines - 40 pts.
  camberTrail = NASA_Rotor_67.getCamberTrail(Rotor_Percent_Span,40,NASA_Rotor_67.z_d); 
  camberLead = NASA_Rotor_67.getCamberLead(Rotor_Percent_Span,40,NASA_Rotor_67.z_up);
  camberBlade = NASA_Rotor_67.getCamberBlade(Rotor_Percent_Span,40);
  
  // Get upper, lower, upperMiddle, and lowerMiddle boundaries.
  upperB = NASA_Rotor_67.getCamberAndShift(Rotor_Percent_Span,40,NASA_Rotor_67.z_up,NASA_Rotor_67.z_d,PI/NASA_Rotor_67.num_blades);
  lowerB = NASA_Rotor_67.getCamberAndShift(Rotor_Percent_Span,40,NASA_Rotor_67.z_up,NASA_Rotor_67.z_d,-PI/NASA_Rotor_67.num_blades);
  upperMiddleB = NASA_Rotor_67.getCamberAndShift(Rotor_Percent_Span,40,NASA_Rotor_67.z_up,NASA_Rotor_67.z_d,PI/(2*NASA_Rotor_67.num_blades));
  lowerMiddleB = NASA_Rotor_67.getCamberAndShift(Rotor_Percent_Span,40,NASA_Rotor_67.z_up,NASA_Rotor_67.z_d,-PI/(2*NASA_Rotor_67.num_blades));

  // Locate boundary corners to ensure that boundaries meeting at the leading
  // and trailing edges are separated by approximately 45deg
   
  // Calculate slope of leading and trailing camber lines.
  mLeadc = (camberLead.Xp[camberLead.np-1].y - camberLead.Xp[0].y)/
           (camberLead.Xp[camberLead.np-1].x - camberLead.Xp[0].x);

  mTrailc = (camberTrail.Xp[camberTrail.np-1].y-camberTrail.Xp[0].y)/
            (camberTrail.Xp[camberTrail.np-1].x-camberTrail.Xp[0].x);

  // Calculate slope of tangent to cross-section at the leading and 
  // trailing edge using a central-difference
  z = Rotor_Spline.Xp[iLead].x;
  dz = (Rotor_Spline.Xp[iLead].x - Rotor_Spline.Xp[iLead+1].x)/TEN;
  mLead = (-   getY(z+2*dz,Rotor_Spline)[0].y 
	   + 8*getY(z+  dz,Rotor_Spline)[0].y
	   - 8*getY(z-  dz,Rotor_Spline)[0].y + 
               getY(z-2*dz,Rotor_Spline)[0].y)/(12*dz);

  z = Rotor_Spline.Xp[iTrail].x;
  dz = (Rotor_Spline.Xp[iTrail].x - Rotor_Spline.Xp[iTrail-1].x)/TEN;
  mTrail = (-   getY(z+2*dz,Rotor_Spline)[1].y
	    + 8*getY(z+  dz,Rotor_Spline)[1].y 
	    - 8*getY(z-  dz,Rotor_Spline)[1].y
	    +   getY(z-2*dz,Rotor_Spline)[1].y)/(12*dz);

  // Get slopes of lines approx 45deg to camber line at the LEADING EDGE.
  A = (mLead+mLeadc)/(1-mLead*mLeadc);
  mTop = (-1+sqrt(1+A*A))/A;
  mBot = (-1-sqrt(1+A*A))/A;

  // Find point on upperMiddle spline - leading edge.
  dm = 1e15;
  mIndex = -1;
  for (int i = 0; i < upperMiddleB.np; i++) {
    if (fabs(upperMiddleB.Xp[i].x-leadV.x) > TOLER) {
      m = (upperMiddleB.Xp[i].y-leadV.y)/(upperMiddleB.Xp[i].x-leadV.x);
      if (fabs(m-mTop) < dm) {
	dm = fabs(m-mTop);
	mIndex = i;
      }
    }
  }
  assert(mIndex != -1);
  zlU = upperMiddleB.Xp[mIndex].x;  

  // Find point on lowerMiddle spline - leading edge.
  dm = 1e15;
  mIndex = -1;
  for (int i = 0; i < lowerMiddleB.np; i++) {
    if (fabs(lowerMiddleB.Xp[i].x-leadV.x) > TOLER) {
      m = (lowerMiddleB.Xp[i].y-leadV.y)/(lowerMiddleB.Xp[i].x-leadV.x);
      if (fabs(m-mBot) < dm) {
	dm=fabs(m-mBot);
	mIndex=i;
      }
    }
  }
  assert(mIndex != -1);
  zlL = lowerMiddleB.Xp[mIndex].x;  

  // Get slopes of lines approx 45deg to camber line at the TRAILING EDGE.
  A = (mTrail+mTrailc)/(1-mTrail*mTrailc);
  mBot = (-1+sqrt(1+A*A))/A;
  mTop = (-1-sqrt(1+A*A))/A;
  
  // Find point on upperMiddle spline - trailing edge.
  dm = 1e15;
  mIndex = -1;
  for (int i=0; i < upperMiddleB.np; i++) {
    if (fabs(upperMiddleB.Xp[i].x-trailV.x) > TOLER) {
      m = (upperMiddleB.Xp[i].y-trailV.y)/(upperMiddleB.Xp[i].x-trailV.x);
      if (fabs(m-mTop) < dm) {
	dm=fabs(m-mTop);
	mIndex=i;
      }
    }
  }
  assert(mIndex != -1);
  zrU = upperMiddleB.Xp[mIndex].x;  

  dm = 1e15;
  mIndex = -1;
  for (int i = 0; i < upperMiddleB.np; i++) {
    if (fabs(upperMiddleB.Xp[i].x-trailV.x) > TOLER) {
      m = (upperMiddleB.Xp[i].y-trailV.y)/(upperMiddleB.Xp[i].x-trailV.x);
      if (fabs(m-mBot) < dm) {
	dm = fabs(m-mBot);
	mIndex = i;
      }
    }
  }
  assert(mIndex != -1);
  
  // Determine which to use.
  if (upperMiddleB.Xp[mIndex].x>zrU) {
    zrU = upperMiddleB.Xp[mIndex].x;  
    swapm = mBot;
    mBot = mTop;
    mTop = swapm;
  }

  // Find point on lowerMiddle spline - trailing edge.
  dm = 1e15;
  mIndex = -1;
  for (int i=0; i < lowerMiddleB.np; i++) {
    if (fabs(lowerMiddleB.Xp[i].x-trailV.x) > TOLER) {
      m = (lowerMiddleB.Xp[i].y-trailV.y)/(lowerMiddleB.Xp[i].x-trailV.x);
      if (fabs(m-mBot) < dm) {
	dm = fabs(m-mBot);
	mIndex = i;
      }
    }
  }
  assert(mIndex != -1);
  zrL = lowerMiddleB.Xp[mIndex].x;

  // NORTH BOUNDARY SPLINE.
  // First segment.
  BndNorthSpline.allocate(30);
  BndNorthSpline.settype(SPLINE2D_QUINTIC);
  zLead = upperB.Xp[0].x;
  zTrail = (zlL+zlU)/2;
  dz = (zTrail-zLead)/29;
  for (int i = 0; i < 30; i++)
    BndNorthSpline.Xp[i] = getY(zLead+i*dz,upperB)[0];
  BndNorthSpline.Xp[0] = upperB.Xp[0];
  BndNorthSpline.pathlength();
  // Second segment.
  tempBndNorthSpline.allocate(70);
  tempBndNorthSpline.settype(SPLINE2D_QUINTIC);
  zLead = (zlL+zlU)/2;
  zTrail = (zrL+zrU)/2;
  dz = (zTrail-zLead)/69;
  for (int i = 0; i < 70; i++)
    tempBndNorthSpline.Xp[i] = getY(zLead+i*dz,upperB)[0];
  tempBndNorthSpline.pathlength();
  // Concatenate first and second segments.
  BndNorthSpline = Concatenate_Splines(BndNorthSpline,tempBndNorthSpline);
  // Final segment.
  tempBndNorthSpline.allocate(30);
  tempBndNorthSpline.settype(SPLINE2D_QUINTIC);
  zLead = (zrL+zrU)/2;
  zTrail = upperB.Xp[upperB.np-1].x;
  dz = (zTrail-zLead)/29;
  for (int i = 0; i < 30; i++) {
    tempBndNorthSpline.Xp[i].x = zLead+i*dz;
    if (tempBndNorthSpline.Xp[i].x > zTrail)
      tempBndNorthSpline.Xp[i].x = zTrail;
    tempBndNorthSpline.Xp[i].y = getY(tempBndNorthSpline.Xp[i].x,upperB)[0].y;
  }
  tempBndNorthSpline.Xp[tempBndNorthSpline.np-1] = upperB.Xp[upperB.np-1];
  tempBndNorthSpline.pathlength();
  // Concatenate first and final segments.
  BndNorthSpline = Concatenate_Splines(BndNorthSpline,tempBndNorthSpline);
  for (int i = 0; i < BndNorthSpline.np; i++) {
    BndNorthSpline.bc[i] = BC_NONE;
    if (i == 0 || i == BndNorthSpline.np-1) {
      BndNorthSpline.tp[i] = SPLINE2D_POINT_SHARP_CORNER;
    } else {
      BndNorthSpline.tp[i] = SPLINE2D_POINT_NORMAL;
    }
  }
  //    BndNorthSpline.allocate(2);
  //    BndNorthSpline.settype(SPLINE2D_LINEAR);
  //    BndNorthSpline.Xp[0] = upperB.Xp[0];
  //    BndNorthSpline.Xp[1] = upperB.Xp[lowerB.np-1];
  //    BndNorthSpline.bc[0] = BC_NONE;
  //    BndNorthSpline.bc[1] = BC_NONE;
  //    BndNorthSpline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
  //    BndNorthSpline.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
  //    BndNorthSpline.pathlength();

  // SOUTH BOUNDARY SPLINE.
  // Second segment.
  tempBndSouthSpline.allocate(70);
  tempBndSouthSpline.settype(SPLINE2D_QUINTIC);
  zLead = (zlL+zlU)/2;;
  zTrail = (zrL+zrU)/2;;
  dz = (zTrail-zLead)/69;
  for (int i = 0; i < tempBndSouthSpline.np; i++)
    tempBndSouthSpline.Xp[i] = getY(zLead+i*dz,lowerB)[0];
  tempBndSouthSpline.pathlength();
  // First segment.
  BndSouthSpline.allocate(30);
  BndSouthSpline.settype(SPLINE2D_QUINTIC);
  zLead = lowerB.Xp[0].x;
  zTrail = tempBndSouthSpline.Xp[0].x;
  dz = (zTrail-zLead)/29;
  for (int i = 0; i < 30; i++)
    BndSouthSpline.Xp[i] = getY(zLead+i*dz,lowerB)[0];
  BndSouthSpline.Xp[0] = lowerB.Xp[0];
  BndSouthSpline.pathlength();
  // Concatenate first and second segments.
  BndSouthSpline = Concatenate_Splines(BndSouthSpline,tempBndSouthSpline);
  // Final segment.
  tempBndSouthSpline.allocate(30);
  tempBndSouthSpline.settype(SPLINE2D_QUINTIC);
  zLead = (zrL+zrU)/2;;
  zTrail = lowerB.Xp[lowerB.np-1].x;
  dz = (zTrail-zLead)/29;
  for (int i = 0; i < 30; i++) {
    tempBndSouthSpline.Xp[i].x = zLead+i*dz;
    if(tempBndSouthSpline.Xp[i].x>zTrail) tempBndSouthSpline.Xp[i].x = zTrail;
    tempBndSouthSpline.Xp[i].y = getY(tempBndSouthSpline.Xp[i].x,lowerB)[0].y;
  }
  tempBndSouthSpline.Xp[tempBndSouthSpline.np-1] = lowerB.Xp[lowerB.np-1];
  tempBndSouthSpline.pathlength();
  // Concatenate first and final segments.
  BndSouthSpline = Concatenate_Splines(BndSouthSpline,tempBndSouthSpline);
  for (int i = 0; i < BndSouthSpline.np; i++) {
    BndSouthSpline.bc[i] = BC_NONE;
    if (i == 0 || i == BndSouthSpline.np-1) {
      BndSouthSpline.tp[i] = SPLINE2D_POINT_SHARP_CORNER;
    } else {
      BndSouthSpline.tp[i] = SPLINE2D_POINT_NORMAL;
    }
  }
  //    BndSouthSpline.allocate(2);
  //    BndSouthSpline.settype(SPLINE2D_LINEAR);
  //    BndSouthSpline.Xp[0] = lowerB.Xp[0];
  //    BndSouthSpline.Xp[1] = lowerB.Xp[lowerB.np-1];
  //    BndSouthSpline.bc[0] = BC_NONE;
  //    BndSouthSpline.bc[1] = BC_NONE;
  //    BndSouthSpline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
  //    BndSouthSpline.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
  //    BndSouthSpline.pathlength();

  // EAST BOUNDARY SPLINE.
  BndEastSpline.allocate(2);
  BndEastSpline.settype(SPLINE2D_LINEAR);
  BndEastSpline.Xp[0] = BndSouthSpline.Xp[BndSouthSpline.np-1];
  BndEastSpline.Xp[1] = BndNorthSpline.Xp[BndNorthSpline.np-1];
  BndEastSpline.bc[0] = BC_CHARACTERISTIC;
  BndEastSpline.bc[1] = BC_CHARACTERISTIC;
  BndEastSpline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
  BndEastSpline.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
  BndEastSpline.pathlength();

  // WEST BOUNDARY SPLINE.
  BndWestSpline.allocate(2);
  BndWestSpline.settype(SPLINE2D_LINEAR);
  BndWestSpline.Xp[0] = BndSouthSpline.Xp[0];
  BndWestSpline.Xp[1] = BndNorthSpline.Xp[0];
  BndWestSpline.bc[0] = BC_FIXED;
  BndWestSpline.bc[1] = BC_FIXED;
  BndWestSpline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
  BndWestSpline.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
  BndWestSpline.pathlength();

  // Assign values to the stretching function parameters and
  // boundary grid line orthogonality parameters.
  Stretch_I = STRETCHING_FCN_LINEAR;
  Beta_I = ZERO; 
  Tau_I = ZERO;
  Stretch_J = STRETCHING_FCN_LINEAR;
  Beta_J = ZERO;
  Tau_J = ZERO;
  Orthogonal_North = 0;
  Orthogonal_South = 0;
  Orthogonal_East = 0;
  Orthogonal_West = 0;

  // Create the mesh for each block representing the complete grid.
  for (int jBlk = 0; jBlk < Number_of_Blocks_Jdir; jBlk++) {
    for (int iBlk = 0; iBlk < Number_of_Blocks_Idir; iBlk++) {
      // Create the splines defining the north, south, east, and west 
      // boundaries of the block.
      xSW.x = BndWestSpline.Xp[0].x + ((double(iBlk)/double(Number_of_Blocks_Idir))*
				       (BndEastSpline.Xp[0].x - BndWestSpline.Xp[0].x));
      xSW.y = getY(xSW.x,BndSouthSpline)[0].y + ((double(jBlk)/double(Number_of_Blocks_Jdir))*
						 (getY(xSW.x,BndNorthSpline)[0].y - getY(xSW.x,BndSouthSpline)[0].y));
      xSE.x = BndWestSpline.Xp[0].x + ((double(iBlk+1)/double(Number_of_Blocks_Idir))*
				       (BndEastSpline.Xp[0].x - BndWestSpline.Xp[0].x));
      xSE.y = getY(xSE.x,BndSouthSpline)[0].y + ((double(jBlk)/double(Number_of_Blocks_Jdir))*
						 (getY(xSE.x,BndNorthSpline)[0].y - getY(xSE.x,BndSouthSpline)[0].y));
      xNW.x = xSW.x;
      xNW.y = getY(xSW.x,BndSouthSpline)[0].y + ((double(jBlk+1)/double(Number_of_Blocks_Jdir))*
						 (getY(xNW.x,BndNorthSpline)[0].y - getY(xNW.x,BndSouthSpline)[0].y));
      xNE.x = BndWestSpline.Xp[0].x + ((double(iBlk+1)/double(Number_of_Blocks_Idir))*
				       (BndEastSpline.Xp[0].x - BndWestSpline.Xp[0].x));
      xNE.y = getY(xNE.x,BndSouthSpline)[0].y + ((double(jBlk+1)/double(Number_of_Blocks_Jdir))*
						 (getY(xNE.x,BndNorthSpline)[0].y - getY(xNE.x,BndSouthSpline)[0].y));
      // North Spline.
      Bnd_Spline_North.allocate(10);
      Bnd_Spline_North.settype(SPLINE2D_QUINTIC);
      for (int i = 0; i < Bnd_Spline_North.np; i++) {
	Bnd_Spline_North.Xp[i].x = xNW.x + (double(i)/double(Bnd_Spline_North.np-1))*(xNE.x - xNW.x);
	Bnd_Spline_North.Xp[i].y = getY(Bnd_Spline_North.Xp[i].x,BndSouthSpline)[0].y + 
                                   (double(jBlk+1)/double(Number_of_Blocks_Jdir))*
                                   (getY(Bnd_Spline_North.Xp[i].x,BndNorthSpline)[0].y - 
				    getY(Bnd_Spline_North.Xp[i].x,BndSouthSpline)[0].y);
  	Bnd_Spline_North.bc[i] = BC_NONE;
	if (i == 0 || i == Bnd_Spline_North.np-1) {
	  Bnd_Spline_North.tp[i] = SPLINE2D_POINT_SHARP_CORNER;
	} else {
	  Bnd_Spline_North.tp[i] = SPLINE2D_POINT_NORMAL;
	}
      }
      Bnd_Spline_North.pathlength();
      // South Spline.
      Bnd_Spline_South.allocate(10);
      Bnd_Spline_South.settype(SPLINE2D_QUINTIC);
      for (int i = 0; i < Bnd_Spline_South.np; i++) {
	Bnd_Spline_South.Xp[i].x = xSW.x + (double(i)/double(Bnd_Spline_South.np-1))*(xSE.x - xSW.x);
	Bnd_Spline_South.Xp[i].y = getY(Bnd_Spline_South.Xp[i].x,BndSouthSpline)[0].y + 
                                   (double(jBlk)/double(Number_of_Blocks_Jdir))*
	                           (getY(Bnd_Spline_South.Xp[i].x,BndNorthSpline)[0].y - 
				    getY(Bnd_Spline_South.Xp[i].x,BndSouthSpline)[0].y);
	Bnd_Spline_South.bc[i] = BC_NONE;
	if (i == 0 || i == Bnd_Spline_South.np-1) {
	  Bnd_Spline_South.tp[i] = SPLINE2D_POINT_SHARP_CORNER;
	} else {
	  Bnd_Spline_South.tp[i] = SPLINE2D_POINT_NORMAL;
	}
      }
      Bnd_Spline_South.pathlength();
      // East Spline.
      Bnd_Spline_East.allocate(2);
      Bnd_Spline_East.settype(SPLINE2D_LINEAR);
      Bnd_Spline_East.Xp[0] = xSE;
      Bnd_Spline_East.Xp[1] = xNE;
      if (iBlk == Number_of_Blocks_Idir-1) {
	Bnd_Spline_East.bc[0] = BC_CHARACTERISTIC;
	Bnd_Spline_East.bc[1] = BC_CHARACTERISTIC;
      } else {
	Bnd_Spline_East.bc[0] = BC_NONE;
	Bnd_Spline_East.bc[1] = BC_NONE;
      }
      Bnd_Spline_East.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
      Bnd_Spline_East.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
      Bnd_Spline_East.pathlength();
      // West Spline.
      Bnd_Spline_West.allocate(2);
      Bnd_Spline_West.settype(SPLINE2D_LINEAR);
      Bnd_Spline_West.Xp[0] = xSW;
      Bnd_Spline_West.Xp[1] = xNW;
      if (iBlk == 0) {
	Bnd_Spline_West.bc[0] = BC_FIXED;
	Bnd_Spline_West.bc[1] = BC_FIXED;
      } else {
	Bnd_Spline_West.bc[0] = BC_NONE;
	Bnd_Spline_West.bc[1] = BC_NONE;
      }
      Bnd_Spline_West.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
      Bnd_Spline_West.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
      Bnd_Spline_West.pathlength();
      // Create the 2D quadrilateral grid block representing the mesh.
      Grid_ptr[iBlk][jBlk].Create_Quad_Block_Without_Update(Bnd_Spline_North,
							    Bnd_Spline_South,
							    Bnd_Spline_East,
							    Bnd_Spline_West,
							    Number_of_Cells_Idir,
							    Number_of_Cells_Jdir,
							    Number_of_Ghost_Cells,
							    Highest_Order_of_Reconstruction,
							    GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
							    Stretch_I,
							    Beta_I, 
							    Tau_I,
							    Stretch_J,
							    Beta_J,
							    Tau_J,
							    Orthogonal_North,
							    Orthogonal_South,
							    Orthogonal_East,
							    Orthogonal_West);
      // Smooth the 2D quadrilateral grid block.
      Grid_ptr[iBlk][jBlk].Smooth_Quad_Block(min(250,2*max(Number_of_Cells_Idir,
							   Number_of_Cells_Jdir)));
    }
  }

}


/*!
 * Generates a mesh for the driven cavity flow.                       
 *                                                                    
 * Usage: Grid_ptr = Grid_Driven_Cavity(Grid_ptr,                     
 *                                      nblk_i,                       
 *                                      nblk_j,                       
 *                                      TEN,                          
 *                                      FIVE,                         
 *   	                                100,                          
 *  	                                50,                           
 *                                      2);                           
 * This subroutine DOESN'T update the ghost cells or
 * the geometric properties of the grid cells.
 */
void Grid2D_Quad_MultiBlock_HO::Grid_Driven_Cavity_Flow_Without_Update(int &_Number_of_Blocks_Idir_,
								       int &_Number_of_Blocks_Jdir_,
								       const double &Width,
								       const double &Height,
								       const int &Stretching_Type_Idir,
								       const int &Stretching_Type_Jdir,
								       const double &Stretching_Factor_Idir,
								       const double &Stretching_Factor_Jdir,
								       const int Number_of_Cells_Idir,
								       const int Number_of_Cells_Jdir,
								       const int Number_of_Ghost_Cells,
								       const int Highest_Order_of_Reconstruction) {

  int n_cells_i, n_cells_j,
      Stretch_I, Stretch_J,
      Orthogonal_North, Orthogonal_South,
      Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D xc_NW, xc_NE, xc_SE, xc_SW;
  Spline2D_HO Bnd_Spline_North, Bnd_Spline_South,
           Bnd_Spline_East, Bnd_Spline_West;

  // Allocate memory for grid block.
  _Number_of_Blocks_Idir_ = 1;
  _Number_of_Blocks_Jdir_ = 1;
  allocate(_Number_of_Blocks_Idir_, _Number_of_Blocks_Jdir_);

  // Create the mesh for each block representing the complete grid.

  for (int jBlk = 0; jBlk < Number_of_Blocks_Jdir; jBlk++) {
    for (int iBlk = 0; iBlk < Number_of_Blocks_Idir; iBlk++) {

      // Assign values to the locations of the corners of the
      // rectangular box shaped domain.

      xc_NW = Vector2D(-HALF*Width+double(iBlk)*Width/double(Number_of_Blocks_Idir), 
		       -HALF*Height+double(jBlk+1)*Height/double(Number_of_Blocks_Jdir));
      xc_NE = Vector2D(-HALF*Width+double(iBlk+1)*Width/double(Number_of_Blocks_Idir), 
		       -HALF*Height+double(jBlk+1)*Height/double(Number_of_Blocks_Jdir));
      xc_SE = Vector2D(-HALF*Width+double(iBlk+1)*Width/double(Number_of_Blocks_Idir), 
		       -HALF*Height+double(jBlk)*Height/double(Number_of_Blocks_Jdir));
      xc_SW = Vector2D(-HALF*Width+double(iBlk)*Width/double(Number_of_Blocks_Idir), 
		       -HALF*Height+double(jBlk)*Height/double(Number_of_Blocks_Jdir));

      // Create the splines defining the north, south, east, and west
      // boundaries of the rectangular box.
      Bnd_Spline_North.Create_Spline_Line(xc_NW, xc_NE, 2);
      Bnd_Spline_South.Create_Spline_Line(xc_SW, xc_SE, 2);
      Bnd_Spline_East.Create_Spline_Line(xc_SE, xc_NE, 2);
      Bnd_Spline_West.Create_Spline_Line(xc_SW, xc_NW, 2);

      // Set the boundary condition types for each of the boundary
      // splines.
      if (jBlk == Number_of_Blocks_Jdir-1) {
	Bnd_Spline_North.setBCtype(BC_MOVING_WALL_ISOTHERMAL);
      } else {
	Bnd_Spline_North.setBCtype(BC_NONE);
      }
      if (jBlk == 0) {
	Bnd_Spline_South.setBCtype(BC_WALL_VISCOUS_ISOTHERMAL);
      } else {
	Bnd_Spline_South.setBCtype(BC_NONE);
      }
      if (iBlk == Number_of_Blocks_Idir-1) {
	Bnd_Spline_East.setBCtype(BC_WALL_VISCOUS_ISOTHERMAL);
      } else {
	Bnd_Spline_East.setBCtype(BC_NONE);
      }
      if (iBlk == 0) {
	Bnd_Spline_West.setBCtype(BC_WALL_VISCOUS_ISOTHERMAL);
      } else {
	Bnd_Spline_West.setBCtype(BC_NONE);
      }

      // Assign values to the stretching function parameters and
      // boundary grid line orthogonality parameters.
      Stretch_I = Stretching_Type_Idir;
      Stretch_J = Stretching_Type_Jdir;
      Beta_I = Stretching_Factor_Idir; Tau_I = ZERO;
      Beta_J = Stretching_Factor_Jdir; Tau_J = ZERO;
      Orthogonal_North = 0;
      Orthogonal_South = 0;
      Orthogonal_East = 0;
      Orthogonal_West = 0;

      // Create the 2D quadrilateral grid block representing the mesh.
      Grid_ptr[iBlk][jBlk].Create_Quad_Block_Without_Update(Bnd_Spline_North,
							    Bnd_Spline_South,
							    Bnd_Spline_East,
							    Bnd_Spline_West,
							    Number_of_Cells_Idir/Number_of_Blocks_Idir,
							    Number_of_Cells_Jdir/Number_of_Blocks_Jdir,
							    Number_of_Ghost_Cells,
							    Highest_Order_of_Reconstruction,
							    GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
							    Stretch_I,
							    Beta_I, 
							    Tau_I,
							    Stretch_J,
							    Beta_J,
							    Tau_J,
							    Orthogonal_North,
							    Orthogonal_South,
							    Orthogonal_East,
							    Orthogonal_West);
      
    }
  }

}


/*!
 * Generates a quadilateral mesh with clustering        
 * consisting of two grid blocks for predicting viscous 
 * flow and boundary layer development over a flat      
 * plate with adiabatic BCs.                            
 *                                                      
 * Usage: Grid_ptr = Grid_Adiabatic_Flat_Plate(Grid_ptr,
 *                                             nblk_i,  
 *                                             nblk_j,  
 *                                             TWO,     
 *         		                       100,     
 *         		                       100,     
 *                                             2);      
 * This subroutine DOESN'T update the ghost cells or
 * the geometric properties of the grid cells.
 */
void Grid2D_Quad_MultiBlock_HO::Grid_Adiabatic_Flat_Plate_Without_Update(int &_Number_of_Blocks_Idir_,
									 int &_Number_of_Blocks_Jdir_,
									 const double &Length,
									 const int Number_of_Cells_Idir,
									 const int Number_of_Cells_Jdir,
									 const int Number_of_Ghost_Cells,
									 const int Highest_Order_of_Reconstruction) {
  
  int iBlk, n_cells_i, n_cells_j, 
        Stretch_I, Stretch_J,
    Orthogonal_North, Orthogonal_South,
    Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D xc_NW, xc_NE, xc_SE, xc_SW;
  Spline2D_HO Bnd_Spline_North, Bnd_Spline_South,
    Bnd_Spline_East, Bnd_Spline_West;
  
  double distance = TEN;//HUNDRED;
  double height = TEN;//HUNDRED;

  /* Allocate memory for grid blocks.  There are two grid
     blocks for this mesh. */

  _Number_of_Blocks_Idir_ = 3;
  _Number_of_Blocks_Jdir_ = 1;
  allocate(_Number_of_Blocks_Idir_, _Number_of_Blocks_Jdir_);

  /* Create the mesh for each block representing
     the complete grid. */

  for ( iBlk = 0; iBlk <= Number_of_Blocks_Idir-1; ++iBlk ) {

    /* Assign values to the locations of the corners
       of the rectangular box shaped domain representing
       each of the blocks in the grid. */

    if (iBlk == 0) {
      xc_NW = Vector2D(-distance*Length, height*Length);
      xc_NE = Vector2D(ZERO       , height*Length);
      xc_SE = Vector2D(ZERO       , ZERO);
      xc_SW = Vector2D(-distance*Length, ZERO);
    } else if (iBlk == 1) {
      xc_NW = Vector2D(ZERO  , height*Length);
      xc_NE = Vector2D(Length, height*Length);
      xc_SE = Vector2D(Length, ZERO);
      xc_SW = Vector2D(ZERO  , ZERO);
    } else {
      xc_NW = Vector2D(Length, height*Length);
      xc_NE = Vector2D((distance+1.0)*Length, height*Length);
      xc_SE = Vector2D((distance+1.0)*Length, ZERO);
      xc_SW = Vector2D(Length, ZERO);
    }/* endif */
   
    /* Create the splines defining the north, south,
       east, and west boundaries of the rectangular boxes. */

    Bnd_Spline_North.Create_Spline_Line(xc_NW, xc_NE, 2);
    Bnd_Spline_South.Create_Spline_Line(xc_SW, xc_SE, 2);
    Bnd_Spline_East.Create_Spline_Line(xc_SE, xc_NE, 2);
    Bnd_Spline_West.Create_Spline_Line(xc_SW, xc_NW, 2);

    /* Set the boundary condition types for each of the
       boundary splines. */

    if (iBlk == 0) {
      Bnd_Spline_North.setBCtype(BC_CHARACTERISTIC);
      Bnd_Spline_South.setBCtype(BC_REFLECTION);
      Bnd_Spline_East.setBCtype(BC_NONE);
      Bnd_Spline_West.setBCtype(BC_CHARACTERISTIC_VELOCITY);
    } else if (iBlk == 1) {
      Bnd_Spline_North.setBCtype(BC_CHARACTERISTIC);
      Bnd_Spline_South.setBCtype(BC_ADIABATIC_WALL);
      Bnd_Spline_East.setBCtype(BC_NONE);
      Bnd_Spline_West.setBCtype(BC_NONE);
    } else {
      Bnd_Spline_North.setBCtype(BC_CHARACTERISTIC);
      Bnd_Spline_South.setBCtype(BC_REFLECTION);
      Bnd_Spline_East.setBCtype(BC_CHARACTERISTIC_VELOCITY);
      Bnd_Spline_West.setBCtype(BC_NONE);
    }/* endif */

    /* Determine the number of cells for this block. */

    n_cells_i = Number_of_Cells_Idir;// /2;
    n_cells_j = Number_of_Cells_Jdir;

    /* Assign values to the stretching function parameters
       and boundary grid line orthogonality parameters. */

    if (iBlk == 0) {
      Stretch_I = STRETCHING_FCN_MAX_CLUSTERING;
      Beta_I = 1.001; 
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
      Beta_J = 1.00001;
      Tau_J = ZERO;
      Orthogonal_North = 0;
      Orthogonal_South = 0;
      Orthogonal_East = 0;
      Orthogonal_West = 0;
    } else if (iBlk == 1) {
      Stretch_I = STRETCHING_FCN_MINMAX_CLUSTERING;
      Beta_I = 1.005; 
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
      Beta_J = 1.00001;
      Tau_J = ZERO;
      Orthogonal_North = 0;
      Orthogonal_South = 0;
      Orthogonal_East = 0;
      Orthogonal_West = 0;
    } else {
      Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
      Beta_I = 1.001; 
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
      Beta_J = 1.00001;
      Tau_J = ZERO;
      Orthogonal_North = 0;
      Orthogonal_South = 0;
      Orthogonal_East = 0;
      Orthogonal_West = 0;
    }/* endif */

    /* Create the 2D quadrilateral grid block. */

    Grid_ptr[iBlk][0].Create_Quad_Block_Without_Update(Bnd_Spline_North,
						       Bnd_Spline_South,
						       Bnd_Spline_East,
						       Bnd_Spline_West,
						       n_cells_i,
						       n_cells_j,
						       Number_of_Ghost_Cells,
						       Highest_Order_of_Reconstruction,
						       GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
						       Stretch_I,
						       Beta_I, 
						       Tau_I,
						       Stretch_J,
						       Beta_J,
						       Tau_J,
						       Orthogonal_North,
						       Orthogonal_South,
						       Orthogonal_East,
						       Orthogonal_West);
	
  } /* endfor */    

}

/*!
 * Generates a double-block O-type grid for predicting            
 * flow past a circular cylinder with adiabatic BCs.              
 *                                                                
 * Usage: Grid_ptr = Grid_Adiabatic_Circular_Cylinder(Grid_ptr,   
 *                                                    nblk_i,     
 *                                                    nblk_j,     
 *                                                    THREE,      
 *   		                                      100,        
 *  		                                      50);        
 * This subroutine DOESN'T update the ghost cells or
 * the geometric properties of the grid cells.
 */
void Grid2D_Quad_MultiBlock_HO::Grid_Adiabatic_Circular_Cylinder_Without_Update(int &_Number_of_Blocks_Idir_,
										int &_Number_of_Blocks_Jdir_,
										const double &Radius,
										const int Number_of_Cells_Idir,
										const int Number_of_Cells_Jdir,
										const int Number_of_Ghost_Cells,
										const int Highest_Order_of_Reconstruction) {
  
  Grid_Adiabatic_Circular_Cylinder(_Number_of_Blocks_Idir_,
				   _Number_of_Blocks_Jdir_,
				   Radius,
				   96.0*Radius,
				   Number_of_Cells_Idir,
				   Number_of_Cells_Jdir,
				   Number_of_Ghost_Cells,
				   Highest_Order_of_Reconstruction);
}


void Grid2D_Quad_MultiBlock_HO::Grid_Adiabatic_Circular_Cylinder_Without_Update(int &_Number_of_Blocks_Idir_,
										int &_Number_of_Blocks_Jdir_,
										const double &Inner_Radius,
										const double &Outer_Radius,
										const int Number_of_Cells_Idir,
										const int Number_of_Cells_Jdir,
										const int Number_of_Ghost_Cells,
										const int Highest_Order_of_Reconstruction) {
  
  
  int iBlk, n_cells_i, n_cells_j, Stretch_I, Stretch_J,
    Orthogonal_North, Orthogonal_South,
    Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D x1, x2;
  Spline2D_HO Bnd_Spline_North, Bnd_Spline_South,
    Bnd_Spline_East, Bnd_Spline_West;

  /* Allocate memory for grid blocks.  There are two grid
     blocks for this mesh. */

  _Number_of_Blocks_Idir_ = 2;
  _Number_of_Blocks_Jdir_ = 1;
  allocate(_Number_of_Blocks_Idir_,_Number_of_Blocks_Jdir_);

  /* Create the mesh for each block representing
     the complete grid. */

  for ( iBlk = 0; iBlk <= Number_of_Blocks_Idir-1; ++iBlk ) {

    /* Create the splines defining the north, south,
       east, and west boundaries of the grid. */

    if (iBlk == 0) {
      x1 = Vector2D(ZERO , ZERO);
      Bnd_Spline_North.Create_Spline_Circular_Arc(x1,
						  Outer_Radius,   //was 96*R
						  360.00,
						  180.00,
						  361);
      Bnd_Spline_South.Create_Spline_Circular_Arc(x1,
						  Inner_Radius,
						  360.00,
						  180.00,
						  361);
      x1 = Vector2D(Inner_Radius, ZERO);
      x2 = Vector2D(Outer_Radius, ZERO);
      Bnd_Spline_West.Create_Spline_Line(x1, x2, 2);
      x1 = Vector2D(-Inner_Radius, ZERO);
      x2 = Vector2D(-Outer_Radius, ZERO);
      Bnd_Spline_East.Create_Spline_Line(x1, x2, 2);
    } else {
      x1 = Vector2D(ZERO , ZERO);
      Bnd_Spline_North.Create_Spline_Circular_Arc(x1,
						  Outer_Radius,
						  180.00,
						  ZERO,
						  361);
      Bnd_Spline_South.Create_Spline_Circular_Arc(x1,
						  Inner_Radius,
						  180.00,
						  ZERO,
						  361);
      x1 = Vector2D(-Inner_Radius, ZERO);
      x2 = Vector2D(-Outer_Radius, ZERO);
      Bnd_Spline_West.Create_Spline_Line(x1, x2, 2);
      x1 = Vector2D(Inner_Radius, ZERO);
      x2 = Vector2D(Outer_Radius, ZERO);
      Bnd_Spline_East.Create_Spline_Line(x1, x2, 2);
    } /* endif */

    /* Set the boundary condition types for each of the
       boundary splines. */

    if (iBlk == 0) {
      Bnd_Spline_North.setBCtype(BC_FIXED);
      Bnd_Spline_South.setBCtype(BC_ADIABATIC_WALL);
      Bnd_Spline_East.setBCtype(BC_NONE);
      Bnd_Spline_West.setBCtype(BC_NONE);
    } else {
      Bnd_Spline_North.setBCtype(BC_FIXED);
      Bnd_Spline_South.setBCtype(BC_ADIABATIC_WALL);
      Bnd_Spline_East.setBCtype(BC_NONE);
      Bnd_Spline_West.setBCtype(BC_NONE);
    } /* endif */

    /* Determine the number of cells for this block. */

    n_cells_i = Number_of_Cells_Idir/2;
    n_cells_j = Number_of_Cells_Jdir;

    /* Assign values to the stretching function parameters
       and boundary grid line orthogonality parameters. */

    if (iBlk == 0) {
      Stretch_I = STRETCHING_FCN_LINEAR;//MINMAX_CLUSTERING;
      Beta_I = 1.001;
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
      Beta_J = 1.01;//1.000001;
      Tau_J = ZERO;
      Orthogonal_North = 0;
      Orthogonal_South = 0;
      Orthogonal_East = 0;
      Orthogonal_West = 0;
    } else {
      Stretch_I = STRETCHING_FCN_LINEAR;//MINMAX_CLUSTERING;
      Beta_I = 1.001; 
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
      Beta_J = 1.01;//1.000001;
      Tau_J = ZERO;
      Orthogonal_North = 0;
      Orthogonal_South = 0;
      Orthogonal_East = 0;
      Orthogonal_West = 0;
    } /* endif */

    /* Create the 2D quadrilateral grid block. */

    Grid_ptr[iBlk][0].Create_Quad_Block_Without_Update(Bnd_Spline_North,
						       Bnd_Spline_South,
						       Bnd_Spline_East,
						       Bnd_Spline_West,
						       n_cells_i,
						       n_cells_j,
						       Number_of_Ghost_Cells,
						       Highest_Order_of_Reconstruction,
						       GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
						       Stretch_I,
						       Beta_I, 
						       Tau_I,
						       Stretch_J,
						       Beta_J,
						       Tau_J,
						       Orthogonal_North,
						       Orthogonal_South,
						       Orthogonal_East,
						       Orthogonal_West);
	
  } /* endfor */
}

/*!
 * Generates a quadilateral mesh with clustering        
 * consisting of two grid blocks for predicting viscous 
 * flow between two adiabatic plates.                   
 *                                                      
 * Usage: Grid_ptr = Grid_Adiabatic_Couette(Grid_ptr,   
 *                                          nblk_i,     
 *                                          nblk_j,     
 *                                          TWO,        
 *         		                    100,        
 *         		                    100);       
 * This subroutine DOESN'T update the ghost cells or
 * the geometric properties of the grid cells.
 */
void Grid2D_Quad_MultiBlock_HO::Grid_Adiabatic_Couette_Without_Update(int &_Number_of_Blocks_Idir_,
								      int &_Number_of_Blocks_Jdir_,
								      const double &Separation,
								      const int Number_of_Cells_Idir,
								      const int Number_of_Cells_Jdir,
								      const int Number_of_Ghost_Cells,
								      const int Highest_Order_of_Reconstruction) {
  
  int jBlk, n_cells_i, n_cells_j, 
    Stretch_I, Stretch_J,
    Orthogonal_North, Orthogonal_South,
    Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D xc_NW, xc_NE, xc_SE, xc_SW;
  Spline2D_HO Bnd_Spline_North, Bnd_Spline_South,
    Bnd_Spline_East, Bnd_Spline_West;

  /* Allocate memory for grid blocks.  There are two grid
     blocks for this mesh. */

  _Number_of_Blocks_Idir_ = 1;
  _Number_of_Blocks_Jdir_ = 2;
  allocate(_Number_of_Blocks_Idir_, _Number_of_Blocks_Jdir_);

  /* Create the mesh for each block representing
     the complete grid. */

  for ( jBlk = 0; jBlk <= Number_of_Blocks_Jdir-1; ++jBlk ) {
      
    /* Assign values to the locations of the corners
       of the rectangular box shaped domain representing
       each of the blocks in the grid. */

    if (jBlk == 0) {
      xc_NW = Vector2D(-Separation/8.0, ZERO);
      xc_NE = Vector2D( Separation/8.0, ZERO);
      xc_SE = Vector2D( Separation/8.0, -Separation/TWO);
      xc_SW = Vector2D(-Separation/8.0, -Separation/TWO);
    } else {
      xc_NW = Vector2D(-Separation/8.0, Separation/TWO);
      xc_NE = Vector2D( Separation/8.0, Separation/TWO);
      xc_SE = Vector2D( Separation/8.0, ZERO);
      xc_SW = Vector2D(-Separation/8.0, ZERO);
    } /* endif */
   
    /* Create the splines defining the north, south,
       east, and west boundaries of the rectangular boxes. */

    Bnd_Spline_North.Create_Spline_Line(xc_NW, xc_NE, 200);
    Bnd_Spline_South.Create_Spline_Line(xc_SW, xc_SE, 200);
    Bnd_Spline_East.Create_Spline_Line(xc_SE, xc_NE, 200);
    Bnd_Spline_West.Create_Spline_Line(xc_SW, xc_NW, 200);

    /* Set the boundary condition types for each of the
       boundary splines. */

    if (jBlk == 0) {
      Bnd_Spline_North.setBCtype(BC_NONE);
      Bnd_Spline_South.setBCtype(BC_ADIABATIC_WALL);
      Bnd_Spline_East.setBCtype(BC_CHARACTERISTIC);
      Bnd_Spline_West.setBCtype(BC_CHARACTERISTIC);
    } else {
      Bnd_Spline_North.setBCtype(BC_ADIABATIC_WALL);
      Bnd_Spline_South.setBCtype(BC_NONE);
      Bnd_Spline_East.setBCtype(BC_CHARACTERISTIC);
      Bnd_Spline_West.setBCtype(BC_CHARACTERISTIC);
    } /* endif */

    /* Determine the number of cells for this block. */

    n_cells_i = Number_of_Cells_Idir;
    n_cells_j = Number_of_Cells_Jdir/2;

    /* Assign values to the stretching function parameters
       and boundary grid line orthogonality parameters. */

    if (jBlk == 0) {
      Stretch_I = STRETCHING_FCN_LINEAR;
      Beta_I = ZERO; 
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_LINEAR;//MIN_CLUSTERING;
      Beta_J = ZERO;//1.01;
      Tau_J = ZERO;
      Orthogonal_North = 0;
      Orthogonal_South = 0;
      Orthogonal_East = 0;
      Orthogonal_West = 0;
    } else {
      Stretch_I = STRETCHING_FCN_LINEAR;
      Beta_I = ZERO; 
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_LINEAR;//MAX_CLUSTERING;
      Beta_J = ZERO;//1.01;
      Tau_J = ZERO;
      Orthogonal_North = 0;
      Orthogonal_South = 0;
      Orthogonal_East = 0;
      Orthogonal_West = 0;
    } /* endif */

    /* Create the 2D quadrilateral grid block. */

    Grid_ptr[0][jBlk].Create_Quad_Block_Without_Update(Bnd_Spline_North,
						       Bnd_Spline_South,
						       Bnd_Spline_East,
						       Bnd_Spline_West,
						       n_cells_i,
						       n_cells_j,
						       Number_of_Ghost_Cells,
						       Highest_Order_of_Reconstruction,
						       GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
						       Stretch_I,
						       Beta_I, 
						       Tau_I,
						       Stretch_J,
						       Beta_J,
						       Tau_J,
						       Orthogonal_North,
						       Orthogonal_South,
						       Orthogonal_East,
						       Orthogonal_West);
	
  } /* endfor */    

}

/*!
 * Generates a single block quadilateral mesh with      
 * clustering for predicting radiation heat transfer in 
 * a cylindrical duct or pipe.                          
 *                                                      
 * Usage: Grid_ptr = Grid_Pipe(Grid_ptr,                
 *                             nblk_i,                  
 *                             nblk_j,                  
 *                             TEN,                     
 *                             HALF,                    
 *   	                       100,                     
 *  	                       50,                      
 *                             2);FIXME
 *                                                      
 * This subroutine DOESN'T update the ghost cells or
 * the geometric properties of the grid cells.
 */
void Grid2D_Quad_MultiBlock_HO::Grid_Cylindrical_Encl_Without_Update(int &_Number_of_Blocks_Idir_,
								     int &_Number_of_Blocks_Jdir_,
								     const double &Length,
								     const double &Radius,
								     const int &Axisymmetric,
								     const int Number_of_Cells_Idir,
								     const int Number_of_Cells_Jdir,
								     const int Number_of_Ghost_Cells,
								     const int Highest_Order_of_Reconstruction) {

  int Stretch_I, Stretch_J,
    Orthogonal_North, Orthogonal_South,
    Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D xc_NW, xc_NE, xc_SE, xc_SW;
  Spline2D_HO Bnd_Spline_North, Bnd_Spline_South,
    Bnd_Spline_East, Bnd_Spline_West;

  /* Allocate memory for grid block. */

  if (_Number_of_Blocks_Idir_ < 0) _Number_of_Blocks_Idir_ = 1;
  if (_Number_of_Blocks_Jdir_ < 0) _Number_of_Blocks_Jdir_ = 1;
  allocate(_Number_of_Blocks_Idir_, _Number_of_Blocks_Jdir_);

  /* Assign values to the locations of the corners
     of the rectangular box defining the pipe geometry. */

  if(Axisymmetric ==1){
    xc_NW = Vector2D(ZERO , Radius);
    xc_NE = Vector2D(Length, Radius);
    xc_SE = Vector2D(Length, ZERO);
    xc_SW = Vector2D(ZERO , ZERO);

  } else if(Axisymmetric ==2){
    xc_NW = Vector2D(ZERO , Length);
    xc_NE = Vector2D(Radius,Length);
    xc_SE = Vector2D(Radius, ZERO);
    xc_SW = Vector2D(ZERO , ZERO); //axisymmetric 2
  }
  /* Create the splines defining the north, south,
     east, and west boundaries of the grid. */

  Bnd_Spline_North.Create_Spline_Line(xc_NW, xc_NE, 2);
  Bnd_Spline_South.Create_Spline_Line(xc_SW, xc_SE, 2);
  Bnd_Spline_East.Create_Spline_Line(xc_SE, xc_NE, 2);
  Bnd_Spline_West.Create_Spline_Line(xc_SW, xc_NW, 2);

  /* Set the boundary condition types for each of the
     boundary splines. */

  if (Axisymmetric == 1) {
    Bnd_Spline_North.setBCtype(BC_GRAY_WALL);
    Bnd_Spline_South.setBCtype(BC_REFLECTION);
    Bnd_Spline_East.setBCtype(BC_GRAY_WALL);
    Bnd_Spline_West.setBCtype(BC_GRAY_WALL);

    /* Assign values to the stretching function parameters
       and boundary grid line orthogonality parameters. */
      
    Stretch_I = STRETCHING_FCN_LINEAR;
    Beta_I = ZERO; 
    Tau_I = ZERO;
    Stretch_J = STRETCHING_FCN_LINEAR;
    Beta_J = ZERO;
    Tau_J = ZERO;  //axisymmetric #1

  } else if (Axisymmetric == 2) {
    Bnd_Spline_East.setBCtype(BC_GRAY_WALL);
    Bnd_Spline_West.setBCtype(BC_REFLECTION);
    Bnd_Spline_North.setBCtype(BC_GRAY_WALL);
    Bnd_Spline_South.setBCtype(BC_GRAY_WALL);
    /* Assign values to the stretching function parameters
       and boundary grid line orthogonality parameters. */
    Stretch_I =  STRETCHING_FCN_LINEAR;
    Beta_I = ZERO;
    Tau_I = ZERO;
    Stretch_J = STRETCHING_FCN_LINEAR;
    Beta_J =  ZERO; 
    Tau_J = ZERO;  
  }

  Orthogonal_North = 0;
  Orthogonal_South = 0;
  Orthogonal_East = 0;
  Orthogonal_West = 0;

  /* Create the 2D quadrilateral grid block representing
     the mesh. */

  Grid_ptr[0][0].Create_Quad_Block_Without_Update(Bnd_Spline_North,
						  Bnd_Spline_South,
						  Bnd_Spline_East,
						  Bnd_Spline_West,
						  Number_of_Cells_Idir,
						  Number_of_Cells_Jdir,
						  Number_of_Ghost_Cells,
						  Highest_Order_of_Reconstruction,
						  GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
						  Stretch_I,
						  Beta_I, 
						  Tau_I,
						  Stretch_J,
						  Beta_J,
						  Tau_J,
						  Orthogonal_North,
						  Orthogonal_South,
						  Orthogonal_East,
						  Orthogonal_West);
    
}

/*!
 * Generates a uniform 2D Cartesian mesh for a          
 * rectangular box shaped domain.                       
 *                                                      
 * Usage: Grid_ptr = Grid_Rectangular_Box(Grid_ptr,     
 *                                        nblk_i,       
 *                                        nblk_j,       
 *                                        TEN,          
 *                                        FIVE,         
 *   	                                  100,          
 *  	                                  50,           
 *                                        2);           
 * This subroutine DOESN'T update the ghost cells or
 * the geometric properties of the grid cells.
 */
void Grid2D_Quad_MultiBlock_HO::Grid_Rectangular_Encl_Without_Update(int &_Number_of_Blocks_Idir_,
								     int &_Number_of_Blocks_Jdir_,
								     const double &Width,
								     const double &Height,
								     const int Number_of_Cells_Idir,
								     const int Number_of_Cells_Jdir,
								     const int Number_of_Ghost_Cells,
								     const int Highest_Order_of_Reconstruction) {
  
  int iBlk, jBlk, n_cells_i, n_cells_j, 
    Stretch_I, Stretch_J,
    Orthogonal_North, Orthogonal_South,
    Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D xc_NW, xc_NE, xc_SE, xc_SW;
  Spline2D_HO Bnd_Spline_North, Bnd_Spline_South,
    Bnd_Spline_East, Bnd_Spline_West;

  /* Allocate memory for grid block. */

  if (_Number_of_Blocks_Idir_ < 0) _Number_of_Blocks_Idir_ = 1;
  if (_Number_of_Blocks_Jdir_ < 0) _Number_of_Blocks_Jdir_ = 1;
  allocate(_Number_of_Blocks_Idir_,_Number_of_Blocks_Jdir_);

  /* Create the mesh for each block representing
     the complete grid. */

  for ( jBlk = 0; jBlk <= Number_of_Blocks_Jdir-1; ++jBlk ) {
    for ( iBlk = 0; iBlk <= Number_of_Blocks_Idir-1; ++iBlk ) {

      /* Assign values to the locations of the corners
	 of the rectangular box shaped domain. */

      xc_NW = Vector2D(-HALF*Width+double(iBlk)*Width/double(Number_of_Blocks_Idir), 
		       -HALF*Height+double(jBlk+1)*Height/double(Number_of_Blocks_Jdir));
      xc_NE = Vector2D(-HALF*Width+double(iBlk+1)*Width/double(Number_of_Blocks_Idir), 
		       -HALF*Height+double(jBlk+1)*Height/double(Number_of_Blocks_Jdir));
      xc_SE = Vector2D(-HALF*Width+double(iBlk+1)*Width/double(Number_of_Blocks_Idir), 
		       -HALF*Height+double(jBlk)*Height/double(Number_of_Blocks_Jdir));
      xc_SW = Vector2D(-HALF*Width+double(iBlk)*Width/double(Number_of_Blocks_Idir), 
		       -HALF*Height+double(jBlk)*Height/double(Number_of_Blocks_Jdir));

      /* Create the splines defining the north, south,
	 east, and west boundaries of the rectangular box. */

      Bnd_Spline_North.Create_Spline_Line(xc_NW, xc_NE, 2);
      Bnd_Spline_South.Create_Spline_Line(xc_SW, xc_SE, 2);
      Bnd_Spline_East.Create_Spline_Line(xc_SE, xc_NE, 2);
      Bnd_Spline_West.Create_Spline_Line(xc_SW, xc_NW, 2);

      /* Set the boundary condition types for each of the
	 boundary splines. */

      if (jBlk == Number_of_Blocks_Jdir-1) {
	Bnd_Spline_North.setBCtype(BC_GRAY_WALL);
      } else {
	Bnd_Spline_North.setBCtype(BC_NONE);
      } /* endif */
      if (jBlk == 0) {
	Bnd_Spline_South.setBCtype(BC_GRAY_WALL);
      } else {
	Bnd_Spline_South.setBCtype(BC_NONE);
      } /* endif */
      if (iBlk == Number_of_Blocks_Idir-1) {
	Bnd_Spline_East.setBCtype(BC_GRAY_WALL);
      } else {
	Bnd_Spline_East.setBCtype(BC_NONE);
      } /* endif */
      if (iBlk == 0) {
	Bnd_Spline_West.setBCtype(BC_GRAY_WALL);
      } else {
	Bnd_Spline_West.setBCtype(BC_NONE);
      } /* endif */

      /* Assign values to the stretching function parameters
	 and boundary grid line orthogonality parameters. */

      Stretch_I = STRETCHING_FCN_LINEAR;
      Beta_I = ZERO; 
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_LINEAR;
      Beta_J = ZERO;
      Tau_J = ZERO;
      Orthogonal_North = 0;
      Orthogonal_South = 0;
      Orthogonal_East = 0;
      Orthogonal_West = 0;

      /* Create the 2D quadrilateral grid block representing
	 the mesh. */

      Grid_ptr[iBlk][jBlk].Create_Quad_Block_Without_Update(Bnd_Spline_North,
							    Bnd_Spline_South,
							    Bnd_Spline_East,
							    Bnd_Spline_West,
							    Number_of_Cells_Idir/Number_of_Blocks_Idir,
							    Number_of_Cells_Jdir/Number_of_Blocks_Jdir,
							    Number_of_Ghost_Cells,
							    Highest_Order_of_Reconstruction,
							    GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
							    Stretch_I,
							    Beta_I, 
							    Tau_I,
							    Stretch_J,
							    Beta_J,
							    Tau_J,
							    Orthogonal_North,
							    Orthogonal_South,
							    Orthogonal_East,
							    Orthogonal_West);

    } /* endfor */
  } /* endfor */

}

/*!
 * Routine: Grid_Tube_2D                                
 * This subroutine DOESN'T update the ghost cells or
 * the geometric properties of the grid cells.
 */
void Grid2D_Quad_MultiBlock_HO::Grid_Tube_2D_Without_Update(int &_Number_of_Blocks_Idir_,
							    int &_Number_of_Blocks_Jdir_,
							    const double &Radius,
							    const int Number_of_Cells_Idir,
							    const int Number_of_Cells_Jdir,
							    const int Number_of_Ghost_Cells,
							    const int Highest_Order_of_Reconstruction,
							    const int i_Stretching_Radial_Dir,
							    const double &Stretching_Radial_Dir) {

  int iBlk, Stretch_I, Stretch_J,
    Orthogonal_North, Orthogonal_South,
    Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D x1, x2;
  Vector2D xc_NW, xc_NE, xc_SE, xc_SW;
  Spline2D_HO Bnd_Spline_North, Bnd_Spline_South,
    Bnd_Spline_East, Bnd_Spline_West;

  _Number_of_Blocks_Idir_ = 5;
  _Number_of_Blocks_Jdir_ = 1;
  allocate(_Number_of_Blocks_Idir_,_Number_of_Blocks_Jdir_);
    
  /* Create the mesh for each block representing
     the complete grid. */

  xc_NW = Vector2D(-0.35*Radius,  0.35*Radius);
  xc_NE = Vector2D( 0.35*Radius,  0.35*Radius);
  xc_SE = Vector2D( 0.35*Radius, -0.35*Radius);
  xc_SW = Vector2D(-0.35*Radius, -0.35*Radius);

  for ( iBlk = 0; iBlk <= Number_of_Blocks_Idir-1; ++iBlk ) {
    /* Create the splines defining the north, south,
       east, and west boundaries of the grid. */

    if (iBlk == 0) {
      Bnd_Spline_North.Create_Spline_Line(xc_NW, xc_NE, 2);
      Bnd_Spline_South.Create_Spline_Line(xc_SW, xc_SE, 2);
      Bnd_Spline_East.Create_Spline_Line(xc_SE, xc_NE, 2);
      Bnd_Spline_West.Create_Spline_Line(xc_SW, xc_NW, 2);

    } else if (iBlk == 1) {
      x1 = Vector2D(ZERO,ZERO);
      Bnd_Spline_North.Create_Spline_Circular_Arc(x1,
						  Radius,
						  135.00,
						  45.00,
						  181);
      Bnd_Spline_South.Create_Spline_Line(xc_NW, xc_NE, 2);
      x2 = Bnd_Spline_North.Xp[0];
      Bnd_Spline_West.Create_Spline_Line(xc_NW, x2, 2);
      x2 = Bnd_Spline_North.Xp[Bnd_Spline_North.np-1];
      Bnd_Spline_East.Create_Spline_Line(xc_NE, x2, 2);

    } else if (iBlk == 2) {
      x1 = Vector2D(ZERO,ZERO);
      Bnd_Spline_North.Create_Spline_Circular_Arc(x1,
						  Radius,
						  45.00,
						  -45.00,
						  181);
      Bnd_Spline_South.Create_Spline_Line(xc_NE, xc_SE, 2);
      x2 = Bnd_Spline_North.Xp[0];
      Bnd_Spline_West.Create_Spline_Line(xc_NE, x2, 2);
      x2 = Bnd_Spline_North.Xp[Bnd_Spline_North.np-1];
      Bnd_Spline_East.Create_Spline_Line(xc_SE, x2, 2);

    } else if (iBlk == 3) {
      x1 = Vector2D(ZERO,ZERO);
      Bnd_Spline_North.Create_Spline_Circular_Arc(x1,
						  Radius,
						  -45.00,
						  -135.00,
						  181);
      Bnd_Spline_South.Create_Spline_Line(xc_SE, xc_SW, 2);
      x2 = Bnd_Spline_North.Xp[0];
      Bnd_Spline_West.Create_Spline_Line(xc_SE, x2, 2);
      x2 = Bnd_Spline_North.Xp[Bnd_Spline_North.np-1];
      Bnd_Spline_East.Create_Spline_Line(xc_SW, x2, 2);

    } else if (iBlk == 4) {
      x1 = Vector2D(ZERO,ZERO);
      Bnd_Spline_North.Create_Spline_Circular_Arc(x1,
						  Radius,
						  -135.00,
						  -225.00,
						  181);
      Bnd_Spline_South.Create_Spline_Line(xc_SW, xc_NW, 2);
      x2 = Bnd_Spline_North.Xp[0];
      Bnd_Spline_West.Create_Spline_Line(xc_SW, x2, 2);
      x2 = Bnd_Spline_North.Xp[Bnd_Spline_North.np-1];
      Bnd_Spline_East.Create_Spline_Line(xc_NW, x2, 2);

    } /* endif */

    /* Set the boundary condition types for each of the
       boundary splines. */

    Bnd_Spline_North.setBCtype(BC_NONE);
    Bnd_Spline_South.setBCtype(BC_NONE);
    Bnd_Spline_East.setBCtype(BC_NONE);
    Bnd_Spline_West.setBCtype(BC_NONE);

    /* Assign values to the stretching function parameters
       and boundary grid line orthogonality parameters. */

    if (iBlk == 0) {
      Stretch_I = STRETCHING_FCN_LINEAR;
      Beta_I = ZERO;
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_LINEAR;
      Beta_J = ZERO;
      Tau_J = ZERO;

    } else {
      Stretch_I = STRETCHING_FCN_LINEAR;
      Beta_I = ZERO;
      Tau_I = ZERO;
      Stretch_J = i_Stretching_Radial_Dir;
      Beta_J = Stretching_Radial_Dir;
      Tau_J = ZERO;

    } /* endif */

    Orthogonal_North = 0;
    Orthogonal_South = 0;
    Orthogonal_East = 0;
    Orthogonal_West = 0;

    /* Create the 2D quadrilateral grid block. */

    Grid_ptr[iBlk][0].Create_Quad_Block_Without_Update(Bnd_Spline_North,
						       Bnd_Spline_South,
						       Bnd_Spline_East,
						       Bnd_Spline_West,
						       Number_of_Cells_Idir,
						       Number_of_Cells_Jdir,
						       Number_of_Ghost_Cells,
						       Highest_Order_of_Reconstruction,
						       GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
						       Stretch_I,
						       Beta_I, 
						       Tau_I,
						       Stretch_J,
						       Beta_J,
						       Tau_J,
						       Orthogonal_North,
						       Orthogonal_South,
						       Orthogonal_East,
						       Orthogonal_West);
	
  } /* endfor */

}
 
/*!
 * Routine: Grid_Annulus_2D                             
 * This subroutine DOESN'T update the ghost cells or
 * the geometric properties of the grid cells.
 */
void Grid2D_Quad_MultiBlock_HO::Grid_Annulus_2D_Without_Update(int &_Number_of_Blocks_Idir_,
							       int &_Number_of_Blocks_Jdir_,
							       const double &Radius_Inner,
							       const double &Radius_Outer,
							       const int Number_of_Cells_Idir,
							       const int Number_of_Cells_Jdir,
							       const int Number_of_Ghost_Cells,
							       const int Highest_Order_of_Reconstruction,
							       const int i_Stretching_Radial_Dir,
							       const double &Stretching_Radial_Dir) {
  
  int iBlk, Stretch_I, Stretch_J,
    Orthogonal_North, Orthogonal_South,
    Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D x1, x2;
  Spline2D_HO Bnd_Spline_North, Bnd_Spline_South,
    Bnd_Spline_East, Bnd_Spline_West;

  _Number_of_Blocks_Idir_ = 4;
  _Number_of_Blocks_Jdir_ = 1;
  allocate(_Number_of_Blocks_Idir_, _Number_of_Blocks_Jdir_);

  /* Create the mesh for each block representing
     the complete grid. */

  for ( iBlk = 0; iBlk <= Number_of_Blocks_Idir-1; ++iBlk ) {
    /* Create the splines defining the north, south,
       east, and west boundaries of the grid. */

    if (iBlk == 0) {
      x1 = Vector2D(ZERO,ZERO);
      Bnd_Spline_North.Create_Spline_Circular_Arc(x1,
						  Radius_Outer,
						  135.00,
						  45.00,
						  181);
      Bnd_Spline_South.Create_Spline_Circular_Arc(x1,
						  Radius_Inner,
						  135.00,
						  45.00,
						  181);
      x1 = Bnd_Spline_South.Xp[0];
      x2 = Bnd_Spline_North.Xp[0];
      Bnd_Spline_West.Create_Spline_Line(x1, x2, 2);
      x1 = Bnd_Spline_South.Xp[Bnd_Spline_South.np-1];
      x2 = Bnd_Spline_North.Xp[Bnd_Spline_North.np-1];
      Bnd_Spline_East.Create_Spline_Line(x1, x2, 2);

    } else if (iBlk == 1) {
      x1 = Vector2D(ZERO,ZERO);
      Bnd_Spline_North.Create_Spline_Circular_Arc(x1,
						  Radius_Outer,
						  45.00,
						  -45.00,
						  181);
      Bnd_Spline_South.Create_Spline_Circular_Arc(x1,
						  Radius_Inner,
						  45.00,
						  -45.00,
						  181);
      x1 = Bnd_Spline_South.Xp[0];
      x2 = Bnd_Spline_North.Xp[0];
      Bnd_Spline_West.Create_Spline_Line(x1, x2, 2);
      x1 = Bnd_Spline_South.Xp[Bnd_Spline_South.np-1];
      x2 = Bnd_Spline_North.Xp[Bnd_Spline_North.np-1];
      Bnd_Spline_East.Create_Spline_Line(x1, x2, 2);

    } else if (iBlk == 2) {
      x1 = Vector2D(ZERO,ZERO);
      Bnd_Spline_North.Create_Spline_Circular_Arc(x1,
						  Radius_Outer,
						  -45.00,
						  -135.00,
						  181);
      Bnd_Spline_South.Create_Spline_Circular_Arc(x1,
						  Radius_Inner,
						  -45.00,
						  -135.00,
						  181);
      x1 = Bnd_Spline_South.Xp[0];
      x2 = Bnd_Spline_North.Xp[0];
      Bnd_Spline_West.Create_Spline_Line(x1, x2, 2);
      x1 = Bnd_Spline_South.Xp[Bnd_Spline_South.np-1];
      x2 = Bnd_Spline_North.Xp[Bnd_Spline_North.np-1];
      Bnd_Spline_East.Create_Spline_Line(x1, x2, 2);

    } else if (iBlk == 3) {
      x1 = Vector2D(ZERO,ZERO);
      Bnd_Spline_North.Create_Spline_Circular_Arc(x1,
						  Radius_Outer,
						  -135.00,
						  -225.00,
						  181);
      Bnd_Spline_South.Create_Spline_Circular_Arc(x1,
						  Radius_Inner,
						  -135.00,
						  -225.00,
						  181);
      x1 = Bnd_Spline_South.Xp[0];
      x2 = Bnd_Spline_North.Xp[0];
      Bnd_Spline_West.Create_Spline_Line(x1, x2, 2);
      x1 = Bnd_Spline_South.Xp[Bnd_Spline_South.np-1];
      x2 = Bnd_Spline_North.Xp[Bnd_Spline_North.np-1];
      Bnd_Spline_East.Create_Spline_Line(x1, x2, 2);

    } /* endif */

    /* Set the boundary condition types for each of the
       boundary splines. */

    Bnd_Spline_North.setBCtype(BC_NONE);
    Bnd_Spline_South.setBCtype(BC_NONE);
    Bnd_Spline_East.setBCtype(BC_NONE);
    Bnd_Spline_West.setBCtype(BC_NONE);

    /* Assign values to the stretching function parameters
       and boundary grid line orthogonality parameters. */

    Stretch_I = STRETCHING_FCN_LINEAR;
    Beta_I = ZERO;
    Tau_I = ZERO;
    Stretch_J = i_Stretching_Radial_Dir;
    Beta_J = Stretching_Radial_Dir;
    Tau_J = ZERO;

    Orthogonal_North = 0;
    Orthogonal_South = 0;
    Orthogonal_East = 0;
    Orthogonal_West = 0;

    /* Create the 2D quadrilateral grid block. */

    Grid_ptr[iBlk][0].Create_Quad_Block_Without_Update(Bnd_Spline_North,
						       Bnd_Spline_South,
						       Bnd_Spline_East,
						       Bnd_Spline_West,
						       Number_of_Cells_Idir,
						       Number_of_Cells_Jdir,
						       Number_of_Ghost_Cells,
						       Highest_Order_of_Reconstruction,
						       GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
						       Stretch_I,
						       Beta_I, 
						       Tau_I,
						       Stretch_J,
						       Beta_J,
						       Tau_J,
						       Orthogonal_North,
						       Orthogonal_South,
						       Orthogonal_East,
						       Orthogonal_West);
	
  } /* endfor */

}

/*
 * Output the cell data to the specified output stream.
 */
void Grid2D_Quad_MultiBlock_HO::Output_Cells_Data(ostream &Out_File){

  int i, j, block_number;
  
  block_number = 0;
  
  for ( j = 0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
    for ( i = 0; i <= Number_of_Blocks_Idir-1 ; ++i ) {
      Grid_ptr[i][j].Output_Cells_Data(block_number,
				       Out_File);
      block_number = block_number + 1;
    }  /* endfor */
  }  /* endfor */
}

/*
 * Output the cell properties that are invariant to 
 * translation and rotation to the specified output stream.
 */
void Grid2D_Quad_MultiBlock_HO::Output_Cells_Translation_Rotation_Invariant_Properties(ostream &Out_File){

  int i, j, block_number;
  
  block_number = 0;
  
  for ( j = 0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
    for ( i = 0; i <= Number_of_Blocks_Idir-1 ; ++i ) {
      Grid_ptr[i][j].Output_Cells_Translation_Rotation_Invariant_Properties(block_number,
									    Out_File);
      block_number = block_number + 1;
    }  /* endfor */
  }  /* endfor */
}

/*!
 * Output operator
 */
ostream &operator << (ostream &Out_File, const Grid2D_Quad_MultiBlock_HO &G){
  
  int i, j;
 
  Out_File << G.Blocks_Idir() << " " 
	   << G.Blocks_Jdir() << "\n";

  for ( j = 0 ; j <= G.Last_jBlock() ; ++j ) {
    for ( i = 0; i <= G.Last_iBlock() ; ++i ) {
      Out_File << G.Grid_ptr[i][j];
    }  /* endfor */
  }  /* endfor */
  
  
  return Out_File;
}

/*!
 * Input operator
 */
istream &operator >> (istream &In_File, Grid2D_Quad_MultiBlock_HO &G){
  
  int i, j;
  
  In_File.setf(ios::skipws);
  In_File >> i >> j;
  In_File.unsetf(ios::skipws);
  
  // allocate/re-allocate memory for the multi-block grid
  G.allocate(i, j);
  
  for ( j = 0 ; j <= G.Last_jBlock() ; ++j ) {
    for ( i = 0; i <= G.Last_iBlock() ; ++i ) {
      In_File >> G.Grid_ptr[i][j];
    }  /* endfor */
  }  /* endfor */
}

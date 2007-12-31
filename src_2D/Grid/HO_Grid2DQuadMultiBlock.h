/* \file HO_Grid2DQuadMultiBlock.h
   \brief Header file defining high-order 2D quadrilateral multi-block grid type. */

#ifndef _HO_GRID2D_QUAD_MULTIBLOCK_INCLUDED
#define _HO_GRID2D_QUAD_MULTIBLOCK_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "HO_Grid2DQuad.h"	    // Include high-order quadrilateral block grid header file.


/* Define the high-order quadrilateral 2D grid multi-block class. */

/*!
 * \class Grid2D_Quad_MultiBlock_HO
 *
 * @brief Definition of the 2D quadrilateral multi-block mesh.
 *
 * \verbatim
 * Member functions
 *
 * \endverbatim
 */
class Grid2D_Quad_MultiBlock_HO{
public:

  //! @name 
  //@{ 
  Grid2D_Quad_Block_HO  **Grid_ptr; //!< 2D array of quadrilateral block grids.
  //@}
			
  //! @name Constructors, destructor and assignment operator
  //@{  
  //! Default constructor.
  Grid2D_Quad_MultiBlock_HO(void);

  //! Destructor.
  ~Grid2D_Quad_MultiBlock_HO(void){ deallocate(); }


  //! Assignment operator.
  Grid2D_Quad_MultiBlock_HO& operator= (const Grid2D_Quad_MultiBlock_HO &G);
  //@}

  //! @name Memory allocation and deallocation
  //@{
  //! Allocate memory for structured quadrilateral multi-block grid.
  void allocate(const int & _Number_of_Blocks_Idir_, const int & _Number_of_Blocks_Jdir_);

  //! Deallocate memory for structured quadrilateral grid block.
  void deallocate(void);
  //@}

  //! @name Indexes for the multi-block grid
  //@{
  //! Get number of blocks in i-direction
  const int & Blocks_Idir(void) const { return Number_of_Blocks_Idir;}
  //! Get number of blocks in j-direction
  const int & Blocks_Jdir(void) const { return Number_of_Blocks_Jdir;}
  //! Index of the last block in i-direction
  int Last_iBlock(void) const {return Number_of_Blocks_Idir-1;}
  //! Index of the last block in j-direction
  int Last_jBlock(void) const {return Number_of_Blocks_Jdir-1;}


  //! @name Bracket operator.
  //@{
  Grid2D_Quad_Block_HO & operator()(const int &iBlk, const int &jBlk){ return Grid_ptr[iBlk][jBlk];}
  const Grid2D_Quad_Block_HO & operator()(const int &iBlk, const int &jBlk) const { return Grid_ptr[iBlk][jBlk];}
  //@}

  //! Create quad block
  //@{
  //@}

  //! @name Broadcast functions (MPI)
  //@{
  void Broadcast_Multi_Block_Grid(void);
  
  //@}
  
  //! @name Smooth quad grid
  //@{
  //@}
    
  //!@name Set Boundary Conditions.
  //@{
  //@}
  
  //!@name Update geometry and geometric properties.
  //@{
  //@}
  
  //!@name Input/Output functions
  //@{
  void Write_Multi_Block_Grid_Definition(ostream &Out_File);
  void Read_Multi_Block_Grid_Definition(istream &In_File);
  void Write_Multi_Block_Grid(ostream &Out_File);
  void Read_Multi_Block_Grid(istream &In_File);
  //@}

  //!@name Copy block
  //@{
  //@}
  
  //!@name Block manipulation
  //@{
  void Translate_Multi_Block_Grid(const Vector2D &V);
  void Scale_Multi_Block_Grid(const double &Scaling_Factor);
  void Rotate_Multi_Block_Grid(const double &Angle);
  void Reflect_Multi_Block_Grid(void);
  int Check_Multi_Block_Grid(void);
  //@}
  
  //!@name Output functions for plotting.
  //@{
  void Output_Tecplot(ostream &Out_File);
  void Output_Nodes_Tecplot(ostream &Out_File);
  void Output_Cells_Tecplot(ostream &Out_File);
  void Output_Gnuplot(ostream &Out_File);
  //@}
  

  //!@name Uniform 2D Cartesian mesh for different shapes.
  //@{
  void Grid_Rectangular_Box(int &_Number_of_Blocks_Idir_,
			    int &_Number_of_Blocks_Jdir_,
			    const double &Width,
			    const double &Height,
			    const int Number_of_Cells_Idir,
			    const int Number_of_Cells_Jdir,
			    const int Number_of_Ghost_Cells);
  
  void Grid_Rectangular_Box(const double &Width,
			    const double &Height,
			    const int Stretching_Flag,
			    const int Stretching_Type_Idir,
			    const int Stretching_Type_Jdir,
			    const double &Stretching_Factor_Idir,
			    const double &Stretching_Factor_Jdir,
			    const int Number_of_Cells_Idir,
			    const int Number_of_Cells_Jdir,
			    const int Number_of_Ghost_Cells);

  void Grid_Flat_Plate(const double &Length,
		       const int Flat_Plate_BC_Type,
		       const int Stretching_Flag,
		       const double &Stretching_Factor_Idir,
		       const double &Stretching_Factor_Jdir,
		       const int Number_of_Cells_Idir,
		       const int Number_of_Cells_Jdir,
		       const int Number_of_Ghost_Cells);

  void Grid_Flat_Plate_NK(const double &Length,
			  const int Stretching_Flag,
			  const double &Stretching_Factor_Idir,
			  const double &Stretching_Factor_Jdir,
			  const int Number_of_Cells_Idir,
			  const int Number_of_Cells_Jdir,
			  const int Number_of_Ghost_Cells);

  void Grid_Flat_Plate3(const double &Length,
			const int &Stretching_Flag,
			const double &Stretching_Factor_Idir,
			const double &Stretching_Factor_Jdir,
			const int Number_of_Cells_Idir,
			const int Number_of_Cells_Jdir,
			const int Number_of_Ghost_Cells);

  void Grid_Flat_Plate4(const double &Length,
			const int &Stretching_Flag,
			const double &Stretching_Factor_Idir,
			const double &Stretching_Factor_Jdir,
			const int Number_of_Cells_Idir,
			const int Number_of_Cells_Jdir,
			const int Number_of_Ghost_Cells);
  void Grid_Flat_Plate9(const double &Length,
			const int &Flat_Plate_BC_Type,
			const int &Stretching_Flag,
			const double &Stretching_Factor_Idir,
			const double &Stretching_Factor_Jdir,
			const int Number_of_Cells_Idir,
			const int Number_of_Cells_Jdir,
			const int Number_of_Ghost_Cells);

  void Grid_1D_Flame(const double &Length,
		     const double &Heigth,
		     const int Number_of_Cells_Idir,
		     const int Number_of_Cells_Jdir,
		     const int Number_of_Ghost_Cells);
  
  void Grid_2D_Laminar_Flame(const double &Length,
			     const double &Heigth,
			     const int Number_of_Cells_Idir,
			     const int Number_of_Cells_Jdir, 
			     const int Number_of_Ghost_Cells,
			     const int Flame_Type_Flag);

  void Grid_Pipe(const double &Length,
		 const double &Radius,
		 const int Stretching_Flag,
		 const double Stretching_Factor,
		 const int Number_of_Cells_Idir,
		 const int Number_of_Cells_Jdir,
		 const int Number_of_Ghost_Cells);  

  void Grid_Pipe(const double &Length,
		 const double &Radius,
		 const int &Axisymmetric,
		 const int Number_of_Cells_Idir,
		 const int Number_of_Cells_Jdir,
		 const int Number_of_Ghost_Cells);
  
  void  Grid_Blunt_Body(const double &Radius,
			const double &Mach_Number,
			const int Number_of_Cells_Idir,
			const int Number_of_Cells_Jdir,
			const int Number_of_Ghost_Cells);

  void Grid_Rocket_Motor(const double &Length_Chamber,
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
			 const int Number_of_Ghost_Cells);
  
  void Grid_Nozzleless_Rocket_Motor(const double &Length_Chamber,
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
				    const int Number_of_Ghost_Cells);

  void Grid_Nozzle(const double &Length_Nozzle,
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
		   const int Number_of_Ghost_Cells);
  
  void Grid_Circular_Cylinder(const double &Radius,
			      const int Stretching_Type_Idir,
			      const int Stretching_Type_Jdir,
			      const double &Stretching_Factor_Idir,
			      const double &Stretching_Factor_Jdir,
			      const int Number_of_Cells_Idir,
			      const int Number_of_Cells_Jdir,
			      const int Number_of_Ghost_Cells);

  void Grid_Circular_Cylinder(const double &Inner_Radius,
			      const double &Outer_Radius,
			      const int Stretching_Type_Idir,
			      const int Stretching_Type_Jdir,
			      const double &Stretching_Factor_Idir,
			      const double &Stretching_Factor_Jdir,
			      const int Number_of_Cells_Idir,
			      const int Number_of_Cells_Jdir,
			      const int Number_of_Ghost_Cells)
  //@}
  
  //!@name AMR related functions
  //@{
  //@}
  
  //! @name Binary arithmetic operators.
  //@{
  //@}

  //! @name Friend binary arithmetic operators.
  //@{
  //@}

  //! @name Input-output operators.
  //@{
  friend ostream &operator << (ostream &out_file, const Grid2D_Quad_MultiBlock_HO &G){};
  friend istream &operator >> (istream &in_file, Grid2D_Quad_MultiBlock_HO &G){};
  //@}

private:
  //! @name Mesh indexes
  //@{ 
  int   Number_of_Blocks_Idir; //!< Number of blocks in i-direction (zeta-direction).
  int   Number_of_Blocks_Jdir; //!< Number of blocks in j-direction (eta-direction).
  //@}


  Grid2D_Quad_MultiBlock_HO(const Grid2D_Quad_MultiBlock_HO &G);     //! Private copy constructor.

};


#endif	// _HO_GRID2D_QUAD_MULTIBLOCK_INCLUDED

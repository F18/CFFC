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
  
  void Grid_Rectangular_Box(int &_Number_of_Blocks_Idir_,
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
			    const int Number_of_Ghost_Cells);

  void Grid_Flat_Plate(int &_Number_of_Blocks_Idir_,
		       int &_Number_of_Blocks_Jdir_,
		       const double &Length,
		       const int Flat_Plate_BC_Type,
		       const int Stretching_Flag,
		       const double &Stretching_Factor_Idir,
		       const double &Stretching_Factor_Jdir,
		       const int Number_of_Cells_Idir,
		       const int Number_of_Cells_Jdir,
		       const int Number_of_Ghost_Cells);

  void Grid_Flat_Plate_NK(int &_Number_of_Blocks_Idir_,
			  int &_Number_of_Blocks_Jdir_,
			  const double &Length,
			  const int Stretching_Flag,
			  const double &Stretching_Factor_Idir,
			  const double &Stretching_Factor_Jdir,
			  const int Number_of_Cells_Idir,
			  const int Number_of_Cells_Jdir,
			  const int Number_of_Ghost_Cells);

  void Grid_Flat_Plate3(int &_Number_of_Blocks_Idir_,
			int &_Number_of_Blocks_Jdir_,
			const double &Length,
			const int &Stretching_Flag,
			const double &Stretching_Factor_Idir,
			const double &Stretching_Factor_Jdir,
			const int Number_of_Cells_Idir,
			const int Number_of_Cells_Jdir,
			const int Number_of_Ghost_Cells);

  void Grid_Flat_Plate4(int &_Number_of_Blocks_Idir_,
			int &_Number_of_Blocks_Jdir_,
			const double &Length,
			const int &Stretching_Flag,
			const double &Stretching_Factor_Idir,
			const double &Stretching_Factor_Jdir,
			const int Number_of_Cells_Idir,
			const int Number_of_Cells_Jdir,
			const int Number_of_Ghost_Cells);
  void Grid_Flat_Plate9(int &_Number_of_Blocks_Idir_,
			int &_Number_of_Blocks_Jdir_,
			const double &Length,
			const int &Flat_Plate_BC_Type,
			const int &Stretching_Flag,
			const double &Stretching_Factor_Idir,
			const double &Stretching_Factor_Jdir,
			const int Number_of_Cells_Idir,
			const int Number_of_Cells_Jdir,
			const int Number_of_Ghost_Cells);

  void Grid_1D_Flame(int &_Number_of_Blocks_Idir_,
		     int &_Number_of_Blocks_Jdir_,
		     const double &Length,
		     const double &Heigth,
		     const int Number_of_Cells_Idir,
		     const int Number_of_Cells_Jdir,
		     const int Number_of_Ghost_Cells);
  
  void Grid_2D_Laminar_Flame(int &_Number_of_Blocks_Idir_,
			     int &_Number_of_Blocks_Jdir_,
			     const double &Length,
			     const double &Heigth,
			     const int Number_of_Cells_Idir,
			     const int Number_of_Cells_Jdir, 
			     const int Number_of_Ghost_Cells,
			     const int Flame_Type_Flag);

  void Grid_Pipe(int &_Number_of_Blocks_Idir_,
		 int &_Number_of_Blocks_Jdir_,
		 const double &Length,
		 const double &Radius,
		 const int Stretching_Flag,
		 const double Stretching_Factor,
		 const int Number_of_Cells_Idir,
		 const int Number_of_Cells_Jdir,
		 const int Number_of_Ghost_Cells);  

  void Grid_Pipe(int &_Number_of_Blocks_Idir_,
		 int &_Number_of_Blocks_Jdir_,
		 const double &Length,
		 const double &Radius,
		 const int &Axisymmetric,
		 const int Number_of_Cells_Idir,
		 const int Number_of_Cells_Jdir,
		 const int Number_of_Ghost_Cells);
  
  void  Grid_Blunt_Body(int &_Number_of_Blocks_Idir_,
			int &_Number_of_Blocks_Jdir_,
			const double &Radius,
			const double &Mach_Number,
			const int Number_of_Cells_Idir,
			const int Number_of_Cells_Jdir,
			const int Number_of_Ghost_Cells);

  void Grid_Rocket_Motor(int &_Number_of_Blocks_Idir_,
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
			 const int Number_of_Ghost_Cells);
  
  void Grid_Nozzleless_Rocket_Motor(int &_Number_of_Blocks_Idir_,
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
				    const int Number_of_Ghost_Cells);

  void Grid_Nozzle(int &_Number_of_Blocks_Idir_,
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
		   const int Number_of_Ghost_Cells);
  
  void Grid_Circular_Cylinder(int &_Number_of_Blocks_Idir_,
			      int &_Number_of_Blocks_Jdir_,
			      const double &Radius,
			      const int Stretching_Type_Idir,
			      const int Stretching_Type_Jdir,
			      const double &Stretching_Factor_Idir,
			      const double &Stretching_Factor_Jdir,
			      const int Number_of_Cells_Idir,
			      const int Number_of_Cells_Jdir,
			      const int Number_of_Ghost_Cells);

  void Grid_Circular_Cylinder(int &_Number_of_Blocks_Idir_,
			      int &_Number_of_Blocks_Jdir_,
			      const double &Inner_Radius,
			      const double &Outer_Radius,
			      const int Stretching_Type_Idir,
			      const int Stretching_Type_Jdir,
			      const double &Stretching_Factor_Idir,
			      const double &Stretching_Factor_Jdir,
			      const int Number_of_Cells_Idir,
			      const int Number_of_Cells_Jdir,
			      const int Number_of_Ghost_Cells);

  void Grid_Annulus(int &_Number_of_Blocks_Idir_,
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
		    const int Number_of_Ghost_Cells);

  void Grid_Ellipse(int &_Number_of_Blocks_Idir_,
		    int &_Number_of_Blocks_Jdir_,
		    const double &A,
		    const double &B,
		    const int Number_of_Cells_Idir,
		    const int Number_of_Cells_Jdir,
		    const int Number_of_Ghost_Cells);
  
  void Grid_NACA_Aerofoil(int &_Number_of_Blocks_Idir_,
			  int &_Number_of_Blocks_Jdir_,
			  char *NACA_Aerofoil_Type_ptr,
			  const double &Chord_Length,
			  const int Number_of_Cells_Idir,
			  const int Number_of_Cells_Jdir,
			  const int Number_of_Ghost_Cells);
  
  void Grid_Free_Jet(int &_Number_of_Blocks_Idir_,
		     int &_Number_of_Blocks_Jdir_,
		     const double &Radius,
		     const int Number_of_Cells_Idir,
		     const int Number_of_Cells_Jdir,
		     const int Number_of_Ghost_Cells) ;

  void Grid_Wedge(int &_Number_of_Blocks_Idir_,
		  int &_Number_of_Blocks_Jdir_,
		  const double &Wedge_Angle,
		  const double &Wedge_Length,
		  const int &Wedge_BC_Type,
		  const int &Stretching_Flag,
		  const double &Stretching_Factor_Idir,
		  const double &Stretching_Factor_Jdir,
		  const int Number_of_Cells_Idir,
		  const int Number_of_Cells_Jdir,
		  const int Number_of_Ghost_Cells);

  void Grid_Unsteady_Blunt_Body(int &_Number_of_Blocks_Idir_,
				int &_Number_of_Blocks_Jdir_,
				const double &Radius,
				const double &Mach_Number,
				const int Number_of_Cells_Idir,
				const int Number_of_Cells_Jdir,
				const int Number_of_Ghost_Cells);

  void Grid_Ringleb_Flow(int &_Number_of_Blocks_Idir_,
			 int &_Number_of_Blocks_Jdir_,
			 const double &Inner_Streamline_Number,
			 const double &Outer_Streamline_Number,
			 const double &Isotach_Line,
			 const int Number_of_Cells_Idir,
			 const int Number_of_Cells_Jdir,
			 const int Number_of_Ghost_Cells);
  void Grid_Bump_Channel_Flow(int &_Number_of_Blocks_Idir_,
			      int &_Number_of_Blocks_Jdir_,
			      const int Smooth_Bump,
			      const int Number_of_Cells_Idir,
			      const int Number_of_Cells_Jdir,
			      const int Number_of_Ghost_Cells) ;
  void Grid_Jet_Flow(int &_Number_of_Blocks_Idir_,
		     int &_Number_of_Blocks_Jdir_,
		     const double &Radius,
		     const double &Mach,
		     const int &Stretching_Type_Idir,
		     const int &Stretching_Type_Jdir,
		     const double &Stretching_Factor_Idir,
		     const double &Stretching_Factor_Jdir,
		     const int Number_of_Cells_Idir,
		     const int Number_of_Cells_Jdir,
		     const int Number_of_Ghost_Cells);

  void Grid_Mixing_Layer(int &_Number_of_Blocks_Idir_,
			 int &_Number_of_Blocks_Jdir_,
			 const double &Length,
			 const double &Mach,
			 const int &Stretching_Type_Idir,
			 const int &Stretching_Type_Jdir,
			 const double &Stretching_Factor_Idir,
			 const double &Stretching_Factor_Jdir,
			 const int Number_of_Cells_Idir,
			 const int Number_of_Cells_Jdir,
			 const int Number_of_Ghost_Cells);

  void Grid_Backward_Facing_Step(int &_Number_of_Blocks_Idir_,
				 int &_Number_of_Blocks_Jdir_,
				 const double &Step_Height,
				 const double &Top_Wall_Deflection,
				 const double &Stretching_Factor_Idir,
				 const double &Stretching_Factor_Jdir,
				 const int Number_of_Cells_Idir,
				 const int Number_of_Cells_Jdir,
				 const int Number_of_Ghost_Cells);

  void Grid_Forward_Facing_Step(int &_Number_of_Blocks_Idir_,
				int &_Number_of_Blocks_Jdir_,
				const double &Step_Height,
				const double &Channel_Gap,
				const double &Stretching_Factor_Idir,
				const double &Stretching_Factor_Jdir,
				const int Number_of_Cells_Idir,
				const int Number_of_Cells_Jdir,
				const int Number_of_Ghost_Cells);

  void Grid_Desolvation_Chamber(const int &Chamber_BC_Type,
				int &_Number_of_Blocks_Idir_,
				int &_Number_of_Blocks_Jdir_,
				const int Number_of_Cells_Idir,
				const int Number_of_Cells_Jdir,
				const int Number_of_Ghost_Cells) ;

  void Grid_NASA_Rotor_37(int &_Number_of_Blocks_Idir_,
			  int &_Number_of_Blocks_Jdir_,
			  const double &Rotor_Percent_Span,
			  const int Number_of_Cells_Idir,
			  const int Number_of_Cells_Jdir,
			  const int Number_of_Ghost_Cells);

  void Grid_NASA_Rotor_67(int &_Number_of_Blocks_Idir_,
			  int &_Number_of_Blocks_Jdir_,
			  const double &Rotor_Percent_Span,
			  const int Number_of_Cells_Idir,
			  const int Number_of_Cells_Jdir,
			  const int Number_of_Ghost_Cells);

  void Grid_Driven_Cavity_Flow(int &_Number_of_Blocks_Idir_,
			       int &_Number_of_Blocks_Jdir_,
			       const double &Width,
			       const double &Height,
			       const int &Stretching_Type_Idir,
			       const int &Stretching_Type_Jdir,
			       const double &Stretching_Factor_Idir,
			       const double &Stretching_Factor_Jdir,
			       const int Number_of_Cells_Idir,
			       const int Number_of_Cells_Jdir,
			       const int Number_of_Ghost_Cells);

  void Grid_Adiabatic_Flat_Plate(int &_Number_of_Blocks_Idir_,
				 int &_Number_of_Blocks_Jdir_,
				 const double &Length,
				 const int Number_of_Cells_Idir,
				 const int Number_of_Cells_Jdir,
				 const int Number_of_Ghost_Cells);

  void Grid_Adiabatic_Circular_Cylinder(int &_Number_of_Blocks_Idir_,
					int &_Number_of_Blocks_Jdir_,
					const double &Radius,
					const int Number_of_Cells_Idir,
					const int Number_of_Cells_Jdir,
					const int Number_of_Ghost_Cells);

  void Grid_Adiabatic_Circular_Cylinder(int &_Number_of_Blocks_Idir_,
					int &_Number_of_Blocks_Jdir_,
					const double &Inner_Radius,
					const double &Outer_Radius,
					const int Number_of_Cells_Idir,
					const int Number_of_Cells_Jdir,
					const int Number_of_Ghost_Cells);

  void Grid_Adiabatic_Couette(int &_Number_of_Blocks_Idir_,
			      int &_Number_of_Blocks_Jdir_,
			      const double &Separation,
			      const int Number_of_Cells_Idir,
			      const int Number_of_Cells_Jdir,
			      const int Number_of_Ghost_Cells);

  void Grid_Cylindrical_Encl(int &_Number_of_Blocks_Idir_,
			     int &_Number_of_Blocks_Jdir_,
			     const double &Length,
			     const double &Radius,
			     const int &Axisymmetric,
			     const int Number_of_Cells_Idir,
			     const int Number_of_Cells_Jdir,
			     const int Number_of_Ghost_Cells);

  void Grid_Rectangular_Encl(int &_Number_of_Blocks_Idir_,
			     int &_Number_of_Blocks_Jdir_,
			     const double &Width,
			     const double &Height,
			     const int Number_of_Cells_Idir,
			     const int Number_of_Cells_Jdir,
			     const int Number_of_Ghost_Cells);

  void Grid_Tube_2D(int &_Number_of_Blocks_Idir_,
		    int &_Number_of_Blocks_Jdir_,
		    const double &Radius,
		    const int Number_of_Cells_Idir,
		    const int Number_of_Cells_Jdir,
		    const int Number_of_Ghost_Cells,
		    const int i_Stretching_Radial_Dir,
		    const double &Stretching_Radial_Dir) ;

  void Grid_Annulus_2D(int &_Number_of_Blocks_Idir_,
		       int &_Number_of_Blocks_Jdir_,
		       const double &Radius_Inner,
		       const double &Radius_Outer,
		       const int Number_of_Cells_Idir,
		       const int Number_of_Cells_Jdir,
		       const int Number_of_Ghost_Cells,
		       const int i_Stretching_Radial_Dir,
		       const double &Stretching_Radial_Dir);
  
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
  friend ostream &operator << (ostream &Out_File, const Grid2D_Quad_MultiBlock_HO &G);
  friend istream &operator >> (istream &In_File, Grid2D_Quad_MultiBlock_HO &G);
  //@}

  //! @name Multi-grid operations based on input parameters
  //@{
  template<typename Input_Parameters_Type>
  int Multi_Block_Grid(Input_Parameters_Type &Input_Parameters);

  template<typename Input_Parameters_Type>
  int Write_Multi_Block_Grid_Definition_Using_IP(const Input_Parameters_Type &Input_Parameters);

  template<typename Input_Parameters_Type>
  void Read_Multi_Block_Grid_Definition_Using_IP(Input_Parameters_Type &Input_Parameters);

  template<typename Input_Parameters_Type>
  int Write_Multi_Block_Grid_Using_IP(const Input_Parameters_Type &Input_Parameters);

  template<typename Input_Parameters_Type>
  void Read_Multi_Block_Grid_Using_IP(Input_Parameters_Type &Input_Parameters);

  template<typename Input_Parameters_Type>
  int Output_Tecplot_Using_IP(const Input_Parameters_Type &Input_Parameters);

  template<typename Input_Parameters_Type>
  int Output_Nodes_Tecplot_Using_IP(const Input_Parameters_Type &Input_Parameters);

  template<typename Input_Parameters_Type>
  int Output_Cells_Tecplot_Using_IP(const Input_Parameters_Type &Input_Parameters);
  //@}

private:
  //! @name Mesh indexes
  //@{ 
  int   Number_of_Blocks_Idir; //!< Number of blocks in i-direction (zeta-direction).
  int   Number_of_Blocks_Jdir; //!< Number of blocks in j-direction (eta-direction).
  //@}


  Grid2D_Quad_MultiBlock_HO(const Grid2D_Quad_MultiBlock_HO &G);     //! Private copy constructor.

};


/*!
 * Generates multi-block quadilateral mesh.
 */
template<typename Input_Parameters_Type>
int Grid2D_Quad_MultiBlock_HO::Multi_Block_Grid(Input_Parameters_Type &Input_Parameters){

  int error_flag(0);
  int iBlk, jBlk;
  int HiBlk, HjBlk;

  /* Generate appropriate mesh. */

  switch(Input_Parameters.i_Grid) {
  case GRID_READ_FROM_DEFINITION_FILE :
    Read_Multi_Block_Grid_Definition_Using_IP(Input_Parameters);

    if (Grid_ptr == NULL) {
      cout << "\n " << CFFC_Name() 
	   << " AdvectDiffuse2D ERROR: Unable to open multi-block mesh definition file "
	   << Input_Parameters.Grid_Definition_File_Name << ".\n";
    } /* endif */
    break;
  case GRID_READ_FROM_GRID_DATA_FILE :
    Read_Multi_Block_Grid_Using_IP(Input_Parameters);
    if (Grid_ptr == NULL) {
      cout << "\n " << CFFC_Name() 
	   << " AdvectDiffuse2D ERROR: Unable to open multi-block mesh data file "
	   << Input_Parameters.Grid_File_Name << ".\n";
    } /* endif */
    break;
  case GRID_SQUARE :
    Grid_Rectangular_Box(Input_Parameters.Number_of_Blocks_Idir,
			 Input_Parameters.Number_of_Blocks_Jdir,
			 Input_Parameters.Box_Width,
			 Input_Parameters.Box_Width,
			 Input_Parameters.Number_of_Cells_Idir,
			 Input_Parameters.Number_of_Cells_Jdir,
			 Input_Parameters.Number_of_Ghost_Cells);
    break;
  case GRID_RECTANGULAR_BOX :
    if (!Input_Parameters.i_Mesh_Stretching) {
      Grid_Rectangular_Box(Input_Parameters.Number_of_Blocks_Idir,
			   Input_Parameters.Number_of_Blocks_Jdir,
			   Input_Parameters.Box_Width,
			   Input_Parameters.Box_Height,
			   Input_Parameters.Number_of_Cells_Idir,
			   Input_Parameters.Number_of_Cells_Jdir,
			   Input_Parameters.Number_of_Ghost_Cells);
    } else {
      Grid_Rectangular_Box(Input_Parameters.Number_of_Blocks_Idir,
			   Input_Parameters.Number_of_Blocks_Jdir,
			   Input_Parameters.Box_Width,
			   Input_Parameters.Box_Height,
			   Input_Parameters.i_Mesh_Stretching,
			   Input_Parameters.Mesh_Stretching_Type_Idir,
			   Input_Parameters.Mesh_Stretching_Type_Jdir,
			   Input_Parameters.Mesh_Stretching_Factor_Idir,
			   Input_Parameters.Mesh_Stretching_Factor_Jdir,
			   Input_Parameters.Number_of_Cells_Idir,
			   Input_Parameters.Number_of_Cells_Jdir,
			   Input_Parameters.Number_of_Ghost_Cells);
    }
    break;
  case GRID_PERIODIC_BOX :
    Grid_Rectangular_Box(Input_Parameters.Number_of_Blocks_Idir,
			 Input_Parameters.Number_of_Blocks_Jdir,
			 Input_Parameters.Box_Width,
			 Input_Parameters.Box_Height,
			 Input_Parameters.i_Mesh_Stretching,
			 Input_Parameters.Mesh_Stretching_Type_Idir,
			 Input_Parameters.Mesh_Stretching_Type_Jdir,
			 Input_Parameters.Mesh_Stretching_Factor_Idir,
			 Input_Parameters.Mesh_Stretching_Factor_Jdir,
			 Input_Parameters.Number_of_Cells_Idir,
			 Input_Parameters.Number_of_Cells_Jdir,
			 Input_Parameters.Number_of_Ghost_Cells);

    // Impose the proper boundary conditions for periodic grid
    for ( jBlk = 0; jBlk <= Input_Parameters.Number_of_Blocks_Jdir-1; ++jBlk ) {
      for ( iBlk = 0; iBlk <= Input_Parameters.Number_of_Blocks_Idir-1; ++iBlk ) {
	// North
	if (jBlk == Input_Parameters.Number_of_Blocks_Jdir-1) {
	  if (Input_Parameters.Number_of_Blocks_Jdir==1) {
	    Grid_ptr[iBlk][jBlk].BndNorthSpline.setBCtype(BC_PERIODIC);
	  } else {
	    Grid_ptr[iBlk][jBlk].BndNorthSpline.setBCtype(BC_NONE);
	  }
	} else {
	  Grid_ptr[iBlk][jBlk].BndNorthSpline.setBCtype(BC_NONE);
	} 
	// South
	if (jBlk == 0) {
	  if (Input_Parameters.Number_of_Blocks_Jdir==1) {  
	    Grid_ptr[iBlk][jBlk].BndSouthSpline.setBCtype(BC_PERIODIC);
	  } else {
	    Grid_ptr[iBlk][jBlk].BndSouthSpline.setBCtype(BC_NONE);  
	  }
	} else {
	  Grid_ptr[iBlk][jBlk].BndSouthSpline.setBCtype(BC_NONE);
	}
	// East
	if (iBlk == Input_Parameters.Number_of_Blocks_Idir-1) {
	  if (Input_Parameters.Number_of_Blocks_Idir==1) {  
	    Grid_ptr[iBlk][jBlk].BndEastSpline.setBCtype(BC_PERIODIC);
	  } else {
	    Grid_ptr[iBlk][jBlk].BndEastSpline.setBCtype(BC_NONE);
	  }
	} else {
	  Grid_ptr[iBlk][jBlk].BndEastSpline.setBCtype(BC_NONE);
	}
	// West
	if (iBlk == 0) {
	  if (Input_Parameters.Number_of_Blocks_Idir==1) {  
	    Grid_ptr[iBlk][jBlk].BndWestSpline.setBCtype(BC_PERIODIC);
	  } else {
	    Grid_ptr[iBlk][jBlk].BndWestSpline.setBCtype(BC_NONE);
	  }
	} else {
	  Grid_ptr[iBlk][jBlk].BndWestSpline.setBCtype(BC_NONE);
	}
	Set_BCs(Grid_ptr[iBlk][jBlk]);
	Update_Exterior_Nodes(Grid_ptr[iBlk][jBlk]);
	Update_Cells(Grid_ptr[iBlk][jBlk]);
      } 
    }  
    
    break;   
  case GRID_INTERIOR_INFLOW_OUTFLOW_BOX :
    // This grid has a "cut" along the positive x axis, starting from (0,0).
    // Along this cut an inflow BC is specified on the upper edge and an outflow BC
    // is imposed on the lower edge.
    
    // Ensure that the grid is generated with an even
    //  number of blocks in order to obtain symmetry.
    Input_Parameters.Number_of_Blocks_Idir += Input_Parameters.Number_of_Blocks_Idir % 2;
    Input_Parameters.Number_of_Blocks_Jdir += Input_Parameters.Number_of_Blocks_Jdir % 2;

    // Generate a rectangular box
    Grid_Rectangular_Box(Input_Parameters.Number_of_Blocks_Idir,
			 Input_Parameters.Number_of_Blocks_Jdir,
			 Input_Parameters.Box_Width,
			 Input_Parameters.Box_Height,
			 Input_Parameters.i_Mesh_Stretching,
			 Input_Parameters.Mesh_Stretching_Type_Idir,
			 Input_Parameters.Mesh_Stretching_Type_Jdir,
			 Input_Parameters.Mesh_Stretching_Factor_Idir,
			 Input_Parameters.Mesh_Stretching_Factor_Jdir,
			 Input_Parameters.Number_of_Cells_Idir,
			 Input_Parameters.Number_of_Cells_Jdir,
			 Input_Parameters.Number_of_Ghost_Cells);

    // Impose the proper boundary conditions for the current grid type.
    HiBlk = Input_Parameters.Number_of_Blocks_Idir/2;
    HjBlk = Input_Parameters.Number_of_Blocks_Jdir/2;

    for (iBlk = HiBlk; iBlk <= Input_Parameters.Number_of_Blocks_Idir-1; ++iBlk ){

      // Impose the inflow BC on the upper edge of the cut to all affected blocks
      Grid_ptr[iBlk][HjBlk  ].BndSouthSpline.setBCtype(BC_INFLOW);

      // Update BCs and geometry
      Set_BCs(Grid_ptr[iBlk][HjBlk ]);
      Update_Exterior_Nodes(Grid_ptr[iBlk][HjBlk ]);
      Update_Cells(Grid_ptr[iBlk][HjBlk ]);

      // Impose the outflow BC on the lower edge of the cut to all affected blocks
      Grid_ptr[iBlk][HjBlk-1].BndNorthSpline.setBCtype(BC_OUTFLOW);

      // Update BCs and geometry
      Set_BCs(Grid_ptr[iBlk][HjBlk-1]);
      Update_Exterior_Nodes(Grid_ptr[iBlk][HjBlk-1]);
      Update_Cells(Grid_ptr[iBlk][HjBlk-1]);
    }
    break;   
  case GRID_FLAT_PLATE :
    Grid_Flat_Plate(Input_Parameters.Number_of_Blocks_Idir,
		    Input_Parameters.Number_of_Blocks_Jdir,
		    Input_Parameters.Plate_Length,
		    BC_DIRICHLET,
		    Input_Parameters.i_Mesh_Stretching,
		    Input_Parameters.Mesh_Stretching_Factor_Idir,
		    Input_Parameters.Mesh_Stretching_Factor_Jdir,
		    Input_Parameters.Number_of_Cells_Idir,
		    Input_Parameters.Number_of_Cells_Jdir,
		    Input_Parameters.Number_of_Ghost_Cells);
    break;
  case GRID_PIPE :
    Grid_Pipe(Input_Parameters.Number_of_Blocks_Idir,
	      Input_Parameters.Number_of_Blocks_Jdir,
	      Input_Parameters.Pipe_Length,
	      Input_Parameters.Pipe_Radius,
	      Input_Parameters.i_Mesh_Stretching,
	      Input_Parameters.Mesh_Stretching_Factor_Jdir,
	      Input_Parameters.Number_of_Cells_Idir,
	      Input_Parameters.Number_of_Cells_Jdir,
	      Input_Parameters.Number_of_Ghost_Cells);
    break;
  case GRID_BLUNT_BODY :
    Grid_Blunt_Body(Input_Parameters.Number_of_Blocks_Idir,
		    Input_Parameters.Number_of_Blocks_Jdir,
		    Input_Parameters.Blunt_Body_Radius,
		    Input_Parameters.Blunt_Body_Mach_Number,
		    Input_Parameters.Number_of_Cells_Idir,
		    Input_Parameters.Number_of_Cells_Jdir,
		    Input_Parameters.Number_of_Ghost_Cells);
    break;
  case GRID_ROCKET_MOTOR :
    Grid_Rocket_Motor(Input_Parameters.Number_of_Blocks_Idir,
		      Input_Parameters.Number_of_Blocks_Jdir,
		      Input_Parameters.Chamber_Length,
		      Input_Parameters.Chamber_Radius,
		      Input_Parameters.Chamber_To_Throat_Length,
		      Input_Parameters.Nozzle_Length,
		      Input_Parameters.Nozzle_Radius_Exit,
		      Input_Parameters.Nozzle_Radius_Throat,
		      Input_Parameters.Grain_Radius,
		      Input_Parameters.Nozzle_Type,
		      Input_Parameters.BC_North,
		      Input_Parameters.i_Mesh_Stretching,
		      Input_Parameters.Mesh_Stretching_Type_Jdir,
		      Input_Parameters.Mesh_Stretching_Factor_Idir,
		      Input_Parameters.Mesh_Stretching_Factor_Jdir,
		      Input_Parameters.Number_of_Cells_Idir,
		      Input_Parameters.Number_of_Cells_Jdir,
		      Input_Parameters.Number_of_Ghost_Cells);
    break;
  case GRID_NOZZLELESS_ROCKET_MOTOR :
    Grid_Nozzleless_Rocket_Motor(Input_Parameters.Number_of_Blocks_Idir,
				 Input_Parameters.Number_of_Blocks_Jdir,
				 Input_Parameters.Chamber_Length,
				 Input_Parameters.Chamber_Radius,
				 Input_Parameters.Nozzle_Length,
				 Input_Parameters.Nozzle_Radius_Exit,
				 BC_REFLECTION,
				 Input_Parameters.i_Mesh_Stretching,
				 Input_Parameters.Mesh_Stretching_Type_Jdir,
				 Input_Parameters.Mesh_Stretching_Factor_Idir,
				 Input_Parameters.Mesh_Stretching_Factor_Jdir,
				 Input_Parameters.Number_of_Cells_Idir,
				 Input_Parameters.Number_of_Cells_Jdir,
				 Input_Parameters.Number_of_Ghost_Cells);
    break;
  case GRID_CIRCULAR_CYLINDER :
    Grid_Circular_Cylinder(Input_Parameters.Number_of_Blocks_Idir,
			   Input_Parameters.Number_of_Blocks_Jdir,
			   Input_Parameters.Cylinder_Radius,
			   Input_Parameters.Cylinder_Radius2,
			   Input_Parameters.Mesh_Stretching_Type_Idir,
			   Input_Parameters.Mesh_Stretching_Type_Jdir,
			   Input_Parameters.Mesh_Stretching_Factor_Idir,
			   Input_Parameters.Mesh_Stretching_Factor_Jdir,
			   Input_Parameters.Number_of_Cells_Idir,
			   Input_Parameters.Number_of_Cells_Jdir,
			   Input_Parameters.Number_of_Ghost_Cells);
    break;
  case GRID_ANNULUS :
    Grid_Annulus(Input_Parameters.Number_of_Blocks_Idir,
		 Input_Parameters.Number_of_Blocks_Jdir,
		 Input_Parameters.Cylinder_Radius,
		 Input_Parameters.Cylinder_Radius2,
		 Input_Parameters.Annulus_Theta_Start,
		 Input_Parameters.Annulus_Theta_End,
		 Input_Parameters.Mesh_Stretching_Type_Idir,
		 Input_Parameters.Mesh_Stretching_Type_Jdir,
		 Input_Parameters.Mesh_Stretching_Factor_Idir,
		 Input_Parameters.Mesh_Stretching_Factor_Jdir,
		 Input_Parameters.Number_of_Cells_Idir,
		 Input_Parameters.Number_of_Cells_Jdir,
		 Input_Parameters.Number_of_Ghost_Cells);
    break;
  case GRID_ELLIPSE :
    Grid_Ellipse(Input_Parameters.Number_of_Blocks_Idir,
		 Input_Parameters.Number_of_Blocks_Jdir,
		 Input_Parameters.Ellipse_Length_X_Axis,
		 Input_Parameters.Ellipse_Length_Y_Axis,
		 Input_Parameters.Number_of_Cells_Idir,
		 Input_Parameters.Number_of_Cells_Jdir,
		 Input_Parameters.Number_of_Ghost_Cells);
    break;
  case GRID_NACA_AEROFOIL :
    Grid_NACA_Aerofoil(Input_Parameters.Number_of_Blocks_Idir,
		       Input_Parameters.Number_of_Blocks_Jdir,
		       Input_Parameters.NACA_Aerofoil_Type,
		       Input_Parameters.Chord_Length,
		       Input_Parameters.Number_of_Cells_Idir,
		       Input_Parameters.Number_of_Cells_Jdir,
		       Input_Parameters.Number_of_Ghost_Cells);
    break;
  case GRID_FREE_JET :
    Grid_Free_Jet(Input_Parameters.Number_of_Blocks_Idir,
		  Input_Parameters.Number_of_Blocks_Jdir,
		  Input_Parameters.Orifice_Radius,
		  Input_Parameters.Number_of_Cells_Idir,
		  Input_Parameters.Number_of_Cells_Jdir,
		  Input_Parameters.Number_of_Ghost_Cells);
    break;
  case GRID_ICEMCFD :
    ICEMCFD_Read(Input_Parameters.ICEMCFD_FileNames,
		 *this,
		 Input_Parameters.Number_of_Ghost_Cells,
		 &Input_Parameters.Number_of_Blocks_Idir,
		 &Input_Parameters.Number_of_Blocks_Jdir);
    break;
  default:
    Grid_Rectangular_Box(Input_Parameters.Number_of_Blocks_Idir,
			 Input_Parameters.Number_of_Blocks_Jdir,
			 Input_Parameters.Box_Width,
			 Input_Parameters.Box_Height,
			 Input_Parameters.Number_of_Cells_Idir,
			 Input_Parameters.Number_of_Cells_Jdir,
			 Input_Parameters.Number_of_Ghost_Cells);
    break;
  } /* endswitch */

    /* Reset boundary conditions if required. */

  if (Input_Parameters.BCs_Specified) {
    for (jBlk = 0; jBlk < Input_Parameters.Number_of_Blocks_Jdir; ++jBlk) {
      for ( iBlk = 0; iBlk < Input_Parameters.Number_of_Blocks_Idir; ++iBlk) {
	if (jBlk == Input_Parameters.Number_of_Blocks_Jdir-1)
	  Grid_ptr[iBlk][jBlk].BndNorthSpline.setBCtype(Input_Parameters.BC_North);
	if (jBlk == 0)
	  Grid_ptr[iBlk][jBlk].BndSouthSpline.setBCtype(Input_Parameters.BC_South);
	if (iBlk == Input_Parameters.Number_of_Blocks_Idir-1)
	  Grid_ptr[iBlk][jBlk].BndEastSpline.setBCtype(Input_Parameters.BC_East);
	if (iBlk == 0)
	  Grid_ptr[iBlk][jBlk].BndWestSpline.setBCtype(Input_Parameters.BC_West);
	Set_BCs(Grid_ptr[iBlk][jBlk]);
	Update_Exterior_Nodes(Grid_ptr[iBlk][jBlk]);
	Update_Cells(Grid_ptr[iBlk][jBlk]);
      }
    }
  }
  
  /* First translate quadrilateral mesh as specified by input parameters. */

  if (abs(Input_Parameters.X_Shift) > TOLER) {
    Translate_Multi_Block_Grid(Input_Parameters.X_Shift);
  } /* endif */

    /* Next scale quadrilateral mesh as specified by input parameters. */

  if (fabs(Input_Parameters.X_Scale-ONE) > TOLER) {
    Scale_Multi_Block_Grid(Input_Parameters.X_Scale);
  } /* endif */

    /* Finally rotate quadrilateral mesh as specified by input parameters. */

  if (fabs(Input_Parameters.X_Rotate) > TOLER) {
    Rotate_Multi_Block_Grid(TWO*PI*Input_Parameters.X_Rotate/360.00);
  } /* endif */

    /* Check the validity of the generated mesh. */
  if (Grid_ptr == NULL){
    error_flag = 1;
  } else if (Check_Multi_Block_Grid()){
    error_flag = 1;
  } else {
    error_flag = 0;
  }

  return error_flag;
}

/*!
 * Writes a grid definition file for a multi-block      
 * quadilateral mesh in a format suitable for retrieval 
 * and re-use purposes.  Returns a non-zero value if    
 * unable to write the grid definition file.            
 */
template<typename Input_Parameters_Type>
int Grid2D_Quad_MultiBlock_HO::Write_Multi_Block_Grid_Definition_Using_IP(const Input_Parameters_Type &Input_Parameters){

  char *mesh_definition_file_name_ptr;
  ofstream mesh_definition_file;

  /* Open the grid definition file. */

  mesh_definition_file_name_ptr = Input_Parameters.Grid_Definition_File_Name;
  mesh_definition_file.open(mesh_definition_file_name_ptr, ios::out);
  if (mesh_definition_file.fail()) return (1);

  /* Write grid type information. */

  mesh_definition_file << Input_Parameters.Grid_Type << "\n" 
		       << Input_Parameters.i_Grid << "\n";

  /* Write the grid definition information for each 
     quadrilateral grid block. */

  Write_Multi_Block_Grid_Definition(mesh_definition_file);

  /* Close the grid definition file. */

  mesh_definition_file.close();

  /* Writing of grid definition file complete.  Return zero value. */

  return(0);
}

/*!
 * Reads definition file information for a 2D array of
 * 2D quadrilateral multi-block grids from the        
 * specified input stream.                            
 */
template<typename Input_Parameters_Type>
void Grid2D_Quad_MultiBlock_HO::Read_Multi_Block_Grid_Definition_Using_IP(Input_Parameters_Type &Input_Parameters){
  
  char buffer[256];
  char *mesh_definition_file_name_ptr;
  ifstream mesh_definition_file;
  
  /* Open the grid definition file. */
  mesh_definition_file_name_ptr = Input_Parameters.Grid_Definition_File_Name;
  mesh_definition_file.open(mesh_definition_file_name_ptr, ios::in);
  if (mesh_definition_file.fail()){
    throw runtime_error("Grid2D_Quad_MultiBlock_HO::Read_Multi_Block_Grid_Definition() ERROR! Cannot open mesh definition file!");
  }
  
  /* Read grid type information. */
  mesh_definition_file.getline(buffer, sizeof(buffer));
  strcpy(Input_Parameters.Grid_Type, buffer);
  mesh_definition_file >> Input_Parameters.i_Grid;

  /* Read the grid definition information for each 
     quadrilateral grid block. */

  Read_Multi_Block_Grid_Definition(mesh_definition_file);

  /* Close the grid definition file. */
  mesh_definition_file.close();
}

/*!
 * Writes multi-block quadilateral mesh to a grid data 
 * file in a format suitable for retrieval and re-use  
 * purposes.  Returns a non-zero value if unable to    
 * write the grid data file.                           
 */
template<typename Input_Parameters_Type>
int Grid2D_Quad_MultiBlock_HO::Write_Multi_Block_Grid_Using_IP(const Input_Parameters_Type &Input_Parameters){

  char *mesh_file_name_ptr;
  ofstream mesh_file;

  /* Open the grid data output file. */
  mesh_file_name_ptr = Input_Parameters.Grid_File_Name;
  mesh_file.open(mesh_file_name_ptr, ios::out);
  if (mesh_file.fail()) return (1);

  /* Write grid type information. */
  mesh_file << Input_Parameters.Grid_Type << "\n" 
	    << Input_Parameters.i_Grid << "\n";

  /* Write the grid data for each quadrilateral grid block. */
  Write_Multi_Block_Grid(mesh_file);

  /* Close the grid data output file. */
  mesh_file.close();
}

/*!
 * Reads multi-block quadilateral mesh from
 * a grid data file.
 */
template<typename Input_Parameters_Type>
void Grid2D_Quad_MultiBlock_HO::Read_Multi_Block_Grid_Using_IP(Input_Parameters_Type &Input_Parameters){

  char buffer[256];
  char *mesh_file_name_ptr;
  ifstream mesh_file;

  /* Open the grid data input file. */
  mesh_file_name_ptr = Input_Parameters.Grid_File_Name;
  mesh_file.open(mesh_file_name_ptr, ios::in);
  if (mesh_file.fail()){
    throw runtime_error("Grid2D_Quad_MultiBlock_HO::Read_Multi_Block_Grid() ERROR! Cannot open grid data file!");
  }
  
  /* Read grid type information. */
  mesh_file.getline(buffer, sizeof(buffer));
  strcpy(Input_Parameters.Grid_Type, buffer);
  mesh_file >> Input_Parameters.i_Grid;

  /* Read the grid data for each quadrilateral grid block. */
  Read_Multi_Block_Grid(mesh_file);

  /* Close the grid data input file. */
  mesh_file.close();
}

/*!
 * Writes the nodes of a multi-block quadilateral mesh
 * in a format suitable for plotting with TECPLOT.    
 * Returns a non-zero value if unable to write the    
 * TECPLOT file.                                      
 */
template<typename Input_Parameters_Type>
int Grid2D_Quad_MultiBlock_HO::Output_Tecplot_Using_IP(const Input_Parameters_Type &Input_Parameters){

  int i, j;
  char prefix[256], extension[256], mesh_file_name[256];
  char *mesh_file_name_ptr;
  ofstream mesh_file;  
  
  /* Determine prefix of grid data output file name. */

  i = 0;
  while (1) {
    if (Input_Parameters.Output_File_Name[i] == ' ' ||
	Input_Parameters.Output_File_Name[i] == '.') break;
    prefix[i]=Input_Parameters.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(Input_Parameters.Output_File_Name) ) break;
  } /* endwhile */
  prefix[i] = '\0';
  strcat(prefix, "_mesh");

  /* Determine grid data output file name. */

  strcpy(extension, ".dat");
  strcpy(mesh_file_name, prefix);
  strcat(mesh_file_name, extension);
  mesh_file_name_ptr = mesh_file_name;

  /* Open the grid data output file. */

  mesh_file.open(mesh_file_name_ptr, ios::out);
  if (mesh_file.fail()) return (1);

  /* Write the node locations for each quadrilateral grid block. */
  Output_Tecplot(mesh_file);

  /* Close the grid data output file. */
  mesh_file.close();

}

/*!
 * Writes the nodes of a multi-block quadilateral mesh
 * in a format suitable for plotting with TECPLOT.    
 * Includes boundary nodes.  Returns a non-zero value 
 * if unable to write the TECPLOT file.               
 */
template<typename Input_Parameters_Type>
int Grid2D_Quad_MultiBlock_HO::Output_Nodes_Tecplot_Using_IP(const Input_Parameters_Type &Input_Parameters){

  int i, j;
  char prefix[256], extension[256], mesh_file_name[256];
  char *mesh_file_name_ptr;
  ofstream mesh_file;  

  /* Determine prefix of grid data output file name. */

  i = 0;
  while (1) {
    if (Input_Parameters.Output_File_Name[i] == ' ' ||
	Input_Parameters.Output_File_Name[i] == '.') break;
    prefix[i]=Input_Parameters.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(Input_Parameters.Output_File_Name) ) break;
  } /* endwhile */
  prefix[i] = '\0';
  strcat(prefix, "_mesh_nodes");

  /* Determine grid data output file name. */

  strcpy(extension, ".dat");
  strcpy(mesh_file_name, prefix);
  strcat(mesh_file_name, extension);
  mesh_file_name_ptr = mesh_file_name;

  /* Open the grid data output file. */

  mesh_file.open(mesh_file_name_ptr, ios::out);
  if (mesh_file.fail()) return (1);

  /* Write the node locations for each quadrilateral grid block. */

  Output_Nodes_Tecplot(mesh_file);

  /* Close the grid data output file. */

  mesh_file.close();

}

/*!
 * Writes the cells of a multi-block quadilateral mesh
 * in a format suitable for plotting with TECPLOT.    
 * Returns a non-zero value if unable to write the    
 * TECPLOT file.                                      
 */
template<typename Input_Parameters_Type>
int Grid2D_Quad_MultiBlock_HO::Output_Cells_Tecplot_Using_IP(const Input_Parameters_Type &Input_Parameters){

  int i, j;
  char prefix[256], extension[256], mesh_file_name[256];
  char *mesh_file_name_ptr;
  ofstream mesh_file;  

  /* Determine prefix of grid data output file name. */
  
  i = 0;
  while (1) {
    if (Input_Parameters.Output_File_Name[i] == ' ' ||
	Input_Parameters.Output_File_Name[i] == '.') break;
    prefix[i]=Input_Parameters.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(Input_Parameters.Output_File_Name) ) break;
  } /* endwhile */
  prefix[i] = '\0';
  strcat(prefix, "_mesh_cells");

  /* Determine grid data output file name. */

  strcpy(extension, ".dat");
  strcpy(mesh_file_name, prefix);
  strcat(mesh_file_name, extension);
  mesh_file_name_ptr = mesh_file_name;

  /* Open the grid data output file. */

  mesh_file.open(mesh_file_name_ptr, ios::out);
  if (mesh_file.fail()) return (1);

  /* Write the node locations for each quadrilateral grid block. */
  Output_Cells_Tecplot(mesh_file);

  /* Close the grid data output file. */
  mesh_file.close();

}

#endif	// _HO_GRID2D_QUAD_MULTIBLOCK_INCLUDED

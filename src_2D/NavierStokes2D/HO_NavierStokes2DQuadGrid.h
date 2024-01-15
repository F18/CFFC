/*!\file HO_NavierStokes2DQuadGrid.h
  \brief Header file defining specializations of 2D quadrilateral mesh related routines for the 2D Navier-Stokes solver.
*/

#ifndef _HO_NAVIERSTOKES2D_QUAD_GRID_INCLUDED 
#define _HO_NAVIERSTOKES2D_QUAD_GRID_INCLUDED 

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../Grid/HO_Grid2DQuadMultiBlock.h" /* Include 2D quadrilateral multiblock grid header file */
#include "NavierStokes2DInput.h"	     /* Include Navier-Stokes 2D input parameters header file */


////////////////////////////////////////////////////////
// Specialized routines for NavierStokes2D solver     //
////////////////////////////////////////////////////////

/*!
 * Specialization of Additional_Setup_Read_Multi_Block_Grid_Definition()
 * from Grid2D_Quad_MultiBlock_HO class.
 */
template<> inline
void Grid2D_Quad_MultiBlock_HO::
Additional_Setup_Read_Multi_Block_Grid_Definition(NavierStokes2D_Input_Parameters &Input_Parameters){
  
  // Handle NASA rotors 37 and 67
  if (Input_Parameters.i_Grid == GRID_NASA_ROTOR_37) {
    Input_Parameters.NASA_Rotor37.init(Input_Parameters.Rotor_Flow_Type,
				       Input_Parameters.NASA_Rotor37_Data_Directory);
    Input_Parameters.Wo.set_gas("AIR");
    Input_Parameters.NASA_Rotor37.getPstateREL_up(Input_Parameters.Wo,Input_Parameters.Rotor_Percent_Span);
    Input_Parameters.Pressure = Input_Parameters.Wo.p;
    Input_Parameters.Temperature = Input_Parameters.Wo.T();
    Input_Parameters.Mach_Number = 
      Input_Parameters.NASA_Rotor37.getMachREL_up(Input_Parameters.Rotor_Percent_Span);
    Input_Parameters.Flow_Angle = atan2(Input_Parameters.Wo.v.y, Input_Parameters.Wo.v.x); 
    if (Input_Parameters.Flow_Angle < ZERO) 
      Input_Parameters.Flow_Angle = TWO*PI + Input_Parameters.Flow_Angle;
    Input_Parameters.Flow_Angle = 180.00*Input_Parameters.Flow_Angle/PI;
  } else if (Input_Parameters.i_Grid == GRID_NASA_ROTOR_67) {
    Input_Parameters.NASA_Rotor67.init(Input_Parameters.Rotor_Flow_Type,
				       Input_Parameters.NASA_Rotor67_Data_Directory);
    Input_Parameters.Wo.set_gas("AIR");
    Input_Parameters.NASA_Rotor67.getPstateREL_up(Input_Parameters.Wo,Input_Parameters.Rotor_Percent_Span);
    Input_Parameters.Pressure = Input_Parameters.Wo.p;
    Input_Parameters.Temperature = Input_Parameters.Wo.T();
    Input_Parameters.Mach_Number = 
      Input_Parameters.NASA_Rotor67.getMachREL_up(Input_Parameters.Rotor_Percent_Span);
    Input_Parameters.Flow_Angle = atan2(Input_Parameters.Wo.v.y, Input_Parameters.Wo.v.x); 
    if (Input_Parameters.Flow_Angle < ZERO) 
      Input_Parameters.Flow_Angle = TWO*PI + Input_Parameters.Flow_Angle;
    Input_Parameters.Flow_Angle = 180.00*Input_Parameters.Flow_Angle/PI;
  } /* endif */
  
}


/*!
 * Specialization of Additional_Setup_Read_Multi_Block_Grid_Data()
 * from Grid2D_Quad_MultiBlock_HO class.
 */
template<> inline
void Grid2D_Quad_MultiBlock_HO::
Additional_Setup_Read_Multi_Block_Grid_Data(NavierStokes2D_Input_Parameters &Input_Parameters){

  // Handle NASA rotors 37 and 67
  if (Input_Parameters.i_Grid == GRID_NASA_ROTOR_37) {
    Input_Parameters.NASA_Rotor37.init(Input_Parameters.Rotor_Flow_Type,
				       Input_Parameters.NASA_Rotor37_Data_Directory);
    Input_Parameters.Wo.set_gas("AIR");
    Input_Parameters.NASA_Rotor37.getPstateREL_up(Input_Parameters.Wo,Input_Parameters.Rotor_Percent_Span);
  } else if (Input_Parameters.i_Grid == GRID_NASA_ROTOR_67) {
    Input_Parameters.NASA_Rotor67.init(Input_Parameters.Rotor_Flow_Type,
				       Input_Parameters.NASA_Rotor67_Data_Directory);
    Input_Parameters.Wo.set_gas("AIR");
    Input_Parameters.NASA_Rotor67.getPstateREL_up(Input_Parameters.Wo,Input_Parameters.Rotor_Percent_Span);
  } /* endif */
  
}

/*!
 * Specialization of Additional_Multi_Block_Grid_Setup()
 * from Grid2D_Quad_MultiBlock_HO class.
 */
template<> inline
void Grid2D_Quad_MultiBlock_HO::Additional_Multi_Block_Grid_Setup(NavierStokes2D_Input_Parameters &Input_Parameters){

  switch(Input_Parameters.i_Grid) {
  case GRID_NASA_ROTOR_37 :
    Input_Parameters.Number_of_Blocks_Idir = 3;
    Input_Parameters.Number_of_Blocks_Jdir = 2;
    if (!Input_Parameters.Interface_IP.Component_List.Ni) {
      Input_Parameters.NASA_Rotor37.genMeshH_3x2_AUTO(*this,
						      Input_Parameters.Rotor_Percent_Span,
						      Input_Parameters.NASA_Rotor37.z_up,
						      Input_Parameters.NASA_Rotor37.z_d,
						      Input_Parameters.Number_of_Cells_Idir,
						      Input_Parameters.Number_of_Cells_Jdir,
						      Input_Parameters.Number_of_Ghost_Cells,
						      HighOrder2D_Input::MaximumReconstructionOrder(),
						      min(250, 8*max(Input_Parameters.Number_of_Cells_Idir,
								     Input_Parameters.Number_of_Cells_Jdir)));
    } else {
      Grid_NASA_Rotor_37_Without_Update(Input_Parameters.Number_of_Blocks_Idir,
					Input_Parameters.Number_of_Blocks_Jdir,
					Input_Parameters.Rotor_Percent_Span,
					Input_Parameters.Number_of_Cells_Idir, 
					Input_Parameters.Number_of_Cells_Jdir,
					Input_Parameters.Number_of_Ghost_Cells,
					HighOrder2D_Input::MaximumReconstructionOrder());
    }
    break;
  case GRID_NASA_ROTOR_67 :
    Input_Parameters.Number_of_Blocks_Idir = 3;
    Input_Parameters.Number_of_Blocks_Jdir = 2;
    if (!Input_Parameters.Interface_IP.Component_List.Ni) {
      Input_Parameters.NASA_Rotor67.genMeshH_3x2_AUTO(*this,
						      Input_Parameters.Rotor_Percent_Span,
						      Input_Parameters.NASA_Rotor67.z_up,
						      Input_Parameters.NASA_Rotor67.z_d,
						      Input_Parameters.Number_of_Cells_Idir,
						      Input_Parameters.Number_of_Cells_Jdir,
						      Input_Parameters.Number_of_Ghost_Cells,
						      HighOrder2D_Input::MaximumReconstructionOrder(),
						      min(250, 8*max(Input_Parameters.Number_of_Cells_Idir,
								     Input_Parameters.Number_of_Cells_Jdir)));
    } else {
      Grid_NASA_Rotor_67_Without_Update(Input_Parameters.Number_of_Blocks_Idir,
					Input_Parameters.Number_of_Blocks_Jdir,
					Input_Parameters.Rotor_Percent_Span,
					Input_Parameters.Number_of_Cells_Idir,
					Input_Parameters.Number_of_Cells_Jdir,
					Input_Parameters.Number_of_Ghost_Cells,
					HighOrder2D_Input::MaximumReconstructionOrder());
    }
    break;
  case GRID_DRIVEN_CAVITY_FLOW :
    Grid_Driven_Cavity_Flow_Without_Update(Input_Parameters.Number_of_Blocks_Idir,
					   Input_Parameters.Number_of_Blocks_Jdir,
					   Input_Parameters.Box_Width,
					   Input_Parameters.Box_Height,
					   Input_Parameters.Mesh_Stretching_Type_Idir,
					   Input_Parameters.Mesh_Stretching_Type_Jdir,
					   Input_Parameters.Mesh_Stretching_Factor_Idir,
					   Input_Parameters.Mesh_Stretching_Factor_Jdir,
					   Input_Parameters.Number_of_Cells_Idir,
					   Input_Parameters.Number_of_Cells_Jdir,
					   Input_Parameters.Number_of_Ghost_Cells,
					   HighOrder2D_Input::MaximumReconstructionOrder());
    break;
  case GRID_JET_FLOW:
    Grid_Jet_Flow_Without_Update(Input_Parameters.Number_of_Blocks_Idir,
				 Input_Parameters.Number_of_Blocks_Jdir,
				 Input_Parameters.Pipe_Radius,
				 Input_Parameters.Mach_Number,
				 Input_Parameters.Mesh_Stretching_Type_Idir,
				 Input_Parameters.Mesh_Stretching_Type_Jdir,
				 Input_Parameters.Mesh_Stretching_Factor_Idir,
				 Input_Parameters.Mesh_Stretching_Factor_Jdir,
				 Input_Parameters.Number_of_Cells_Idir,
				 Input_Parameters.Number_of_Cells_Jdir,
				 Input_Parameters.Number_of_Ghost_Cells,
				 HighOrder2D_Input::MaximumReconstructionOrder());
    break; 
  case GRID_BACKWARD_FACING_STEP :
    Grid_Backward_Facing_Step(Input_Parameters.Number_of_Blocks_Idir,
			      Input_Parameters.Number_of_Blocks_Jdir,
			      Input_Parameters.Step_Height,
			      Input_Parameters.Top_Wall_Deflection,
			      Input_Parameters.Mesh_Stretching_Factor_Idir,
			      Input_Parameters.Mesh_Stretching_Factor_Jdir,
			      Input_Parameters.Number_of_Cells_Idir,
			      Input_Parameters.Number_of_Cells_Jdir,
			      Input_Parameters.Number_of_Ghost_Cells,
			      HighOrder2D_Input::MaximumReconstructionOrder());
    break;
  case GRID_MIXING_LAYER:
    Grid_Mixing_Layer_Without_Update(Input_Parameters.Number_of_Blocks_Idir,
				     Input_Parameters.Number_of_Blocks_Jdir,
				     Input_Parameters.Box_Width,
				     Input_Parameters.Mach_Number,
				     Input_Parameters.Mesh_Stretching_Type_Idir,
				     Input_Parameters.Mesh_Stretching_Type_Jdir,
				     Input_Parameters.Mesh_Stretching_Factor_Idir,
				     Input_Parameters.Mesh_Stretching_Factor_Jdir,
				     Input_Parameters.Number_of_Cells_Idir,
				     Input_Parameters.Number_of_Cells_Jdir,
				     Input_Parameters.Number_of_Ghost_Cells,
				     HighOrder2D_Input::MaximumReconstructionOrder());
    break;
  case GRID_FLAT_PLATE :
    Grid_Flat_Plate_Without_Update(Input_Parameters.Number_of_Blocks_Idir,
				   Input_Parameters.Number_of_Blocks_Jdir,
				   Input_Parameters.Plate_Length,
				   Input_Parameters.Flat_Plate_BC_Type,
				   Input_Parameters.i_Mesh_Stretching,
				   Input_Parameters.Mesh_Stretching_Factor_Idir,
				   Input_Parameters.Mesh_Stretching_Factor_Jdir,
				   Input_Parameters.Number_of_Cells_Idir,
				   Input_Parameters.Number_of_Cells_Jdir,
				   Input_Parameters.Number_of_Ghost_Cells,
				   HighOrder2D_Input::MaximumReconstructionOrder());
    break;
  default:
    Grid_Rectangular_Box_Without_Update(Input_Parameters.Number_of_Blocks_Idir,
					Input_Parameters.Number_of_Blocks_Jdir,
					Input_Parameters.Box_Width,
					Input_Parameters.Box_Height,
					Input_Parameters.Number_of_Cells_Idir,
					Input_Parameters.Number_of_Cells_Jdir,
					Input_Parameters.Number_of_Ghost_Cells,
					HighOrder2D_Input::MaximumReconstructionOrder());
  }
}

#endif // _HO_EULER2D_QUAD_GRID_INCLUDED 

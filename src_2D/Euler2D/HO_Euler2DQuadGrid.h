/*!\file HO_Euler2DQuadGrid.h
  \brief Header file defining specializations of 2D quadrilateral mesh related routines for the 2D Euler solver.
*/

#ifndef _HO_EULER2D_QUAD_GRID_INCLUDED 
#define _HO_EULER2D_QUAD_GRID_INCLUDED 

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../Grid/HO_Grid2DQuadMultiBlock.h" /* Include 2D quadrilateral multiblock grid header file */
#include "Euler2DInput.h"		     /* Include Euler 2D input parameters header file */


/////////////////////////////////////////////////
// Specialized routines for Euler2D solver     //
/////////////////////////////////////////////////

/*!
 * Specialization of Additional_Setup_Read_Multi_Block_Grid_Definition()
 * from Grid2D_Quad_MultiBlock_HO class.
 */
template<> inline
void Grid2D_Quad_MultiBlock_HO::Additional_Setup_Read_Multi_Block_Grid_Definition(Euler2D_Input_Parameters &Input_Parameters){
  
  // Handle NASA rotors 37 and 67
  if (Input_Parameters.i_Grid == GRID_NASA_ROTOR_37) {
    Input_Parameters.NASA_Rotor37.init(Input_Parameters.Rotor_Flow_Type,
				       Input_Parameters.NASA_Rotor37_Data_Directory);
    Input_Parameters.Wo.setgas("AIR");
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
    Input_Parameters.Wo.setgas("AIR");
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
void Grid2D_Quad_MultiBlock_HO::Additional_Setup_Read_Multi_Block_Grid_Data(Euler2D_Input_Parameters &Input_Parameters){

  // Handle NASA rotors 37 and 67
  if (Input_Parameters.i_Grid == GRID_NASA_ROTOR_37) {
    Input_Parameters.NASA_Rotor37.init(Input_Parameters.Rotor_Flow_Type,
				       Input_Parameters.NASA_Rotor37_Data_Directory);
    Input_Parameters.Wo.setgas("AIR");
    Input_Parameters.NASA_Rotor37.getPstateREL_up(Input_Parameters.Wo,Input_Parameters.Rotor_Percent_Span);
  } else if (Input_Parameters.i_Grid == GRID_NASA_ROTOR_67) {
    Input_Parameters.NASA_Rotor67.init(Input_Parameters.Rotor_Flow_Type,
				       Input_Parameters.NASA_Rotor67_Data_Directory);
    Input_Parameters.Wo.setgas("AIR");
    Input_Parameters.NASA_Rotor67.getPstateREL_up(Input_Parameters.Wo,Input_Parameters.Rotor_Percent_Span);
  } /* endif */
  
}

/*!
 * Specialization of Additional_Multi_Block_Grid_Setup()
 * from Grid2D_Quad_MultiBlock_HO class.
 */
template<> inline
void Grid2D_Quad_MultiBlock_HO::Additional_Multi_Block_Grid_Setup(Euler2D_Input_Parameters &Input_Parameters){

  switch(Input_Parameters.i_Grid) {
  case GRID_NASA_ROTOR_37 :
    //     Input_Parameters.Number_of_Blocks_Idir = 3;
    //     Input_Parameters.Number_of_Blocks_Jdir = 2;
    //     if (!Input_Parameters.Interface_IP.Component_List.Ni) {
    //       Grid_ptr = Input_Parameters.NASA_Rotor37.genMeshH_3x2_AUTO(Grid_ptr,
    // 								 Input_Parameters.Rotor_Percent_Span,
    // 								 Input_Parameters.NASA_Rotor37.z_up,
    // 								 Input_Parameters.NASA_Rotor37.z_d,
    // 								 Input_Parameters.Number_of_Cells_Idir,
    // 								 Input_Parameters.Number_of_Cells_Jdir,
    // 								 Input_Parameters.Number_of_Ghost_Cells,
    // 								 min(250, 8*max(Input_Parameters.Number_of_Cells_Idir,
    // 										Input_Parameters.Number_of_Cells_Jdir)));
    //     } else {
    //       Grid_ptr = Grid_NASA_Rotor_37(Grid_ptr,
    // 				    Input_Parameters.Number_of_Blocks_Idir,
    // 				    Input_Parameters.Number_of_Blocks_Jdir,
    // 				    Input_Parameters.Rotor_Percent_Span,
    // 				    Input_Parameters.Number_of_Cells_Idir, 
    // 				    Input_Parameters.Number_of_Cells_Jdir,
    // 				    Input_Parameters.Number_of_Ghost_Cells);
    //     }
    break;
  case GRID_NASA_ROTOR_67 :
    //     Input_Parameters.Number_of_Blocks_Idir = 3;
    //     Input_Parameters.Number_of_Blocks_Jdir = 2;
    //     if (!Input_Parameters.Interface_IP.Component_List.Ni) {
    //       Grid_ptr = Input_Parameters.NASA_Rotor67.genMeshH_3x2_AUTO(Grid_ptr,
    // 								 Input_Parameters.Rotor_Percent_Span,
    // 								 Input_Parameters.NASA_Rotor67.z_up,
    // 								 Input_Parameters.NASA_Rotor67.z_d,
    // 								 Input_Parameters.Number_of_Cells_Idir,
    // 								 Input_Parameters.Number_of_Cells_Jdir,
    // 								 Input_Parameters.Number_of_Ghost_Cells,
    // 								 min(250, 8*max(Input_Parameters.Number_of_Cells_Idir,
    // 										Input_Parameters.Number_of_Cells_Jdir)));
    //     } else {
    //       Grid_ptr = Grid_NASA_Rotor_67(Grid_ptr,
    // 				    Input_Parameters.Number_of_Blocks_Idir,
    // 				    Input_Parameters.Number_of_Blocks_Jdir,
    // 				    Input_Parameters.Rotor_Percent_Span,
    // 				    Input_Parameters.Number_of_Cells_Idir,
    // 				    Input_Parameters.Number_of_Cells_Jdir,
    // 				    Input_Parameters.Number_of_Ghost_Cells);
    //     }
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

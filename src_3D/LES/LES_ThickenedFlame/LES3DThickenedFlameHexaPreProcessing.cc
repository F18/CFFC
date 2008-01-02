
#ifndef _LES3DTF_INCLUDED
#include "LES3DThickenedFlame.h"
#endif // _LES3DTF_INCLUDED

/********************************************************
 * Routine: Pre_Processing_Specializations              *
 ********************************************************/
template<>
int Hexa_Pre_Processing_Specializations(HexaSolver_Data &Data,
                                        HexaSolver_Solution_Data<LES3DTF_pState, 
                                                                 LES3DTF_cState> &Solution_Data) {

  int error_flag(0);
  
  RandomFieldRogallo<LES3DTF_pState, LES3DTF_cState>   Velocity_Field_Type(SPECTRUM_HAWORTH_POINSOT);
  Turbulent_Velocity_Field_Multi_Block_List  Velocity_Field;

  Velocity_Field.Create(Data.Initial_Mesh, 
                        Solution_Data.Input.Grid_IP);

  error_flag = Velocity_Field_Type.Create_Homogeneous_Turbulence_Velocity_Field(Data.Initial_Mesh, 
										Velocity_Field,
                                                                                Solution_Data.Input.Grid_IP);
  if (error_flag) return error_flag;

  Assign_Homogeneous_Turbulence_Velocity_Field(Solution_Data.Local_Solution_Blocks.Soln_Blks,
                                               Data.Local_Adaptive_Block_List,
                                               Velocity_Field);

  error_flag = Solution_Data.Local_Solution_Blocks.ICs_Specializations(Solution_Data.Input);
  if (error_flag) return error_flag;

  return error_flag;

}

/********************************************************
 * Routine: Hexa_Post_Processing_Specializations        *
 ********************************************************/
template<>
int Hexa_Post_Processing_Specializations(HexaSolver_Data &Data,
		 	                 HexaSolver_Solution_Data<LES3DTF_pState, 
                                                                  LES3DTF_cState>&Solution_Data) {

   int error_flag(0);
   double u_ave, v_ave, w_ave, sqr_u;
   
   Time_Averaging_of_Velocity_Field(Solution_Data.Local_Solution_Blocks.Soln_Blks,
                                    Data.Local_Adaptive_Block_List,
                                    u_ave,
                                    v_ave,
                                    w_ave);

   Time_Averaging_of_Solution(Solution_Data.Local_Solution_Blocks.Soln_Blks,
                              Data.Local_Adaptive_Block_List,
                              u_ave,
                              v_ave,
                              w_ave,
                              sqr_u);

   Time_Averaging_of_Turbulent_Burning_Rate(Solution_Data.Local_Solution_Blocks.Soln_Blks,
                                            Data.Local_Adaptive_Block_List,
                                            Solution_Data.Input.Grid_IP);

   return error_flag;

}

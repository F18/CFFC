
#ifndef _LES3DFSD_INCLUDED
#include "LES3DFsd.h"
#endif // _LES3DFSD_INCLUDED

/********************************************************
 * Routine: Pre_Processing_Specializations              *
 ********************************************************/
template<>
int Hexa_Pre_Processing_Specializations(HexaSolver_Data &Data,
                                        HexaSolver_Solution_Data<LES3DFsd_pState, 
                                                                 LES3DFsd_cState> &Solution_Data) {

  int error_flag(0);
  
  RandomFieldRogallo<LES3DFsd_pState, LES3DFsd_cState>Create_Turbulence(SPECTRUM_HAWORTH_POINSOT);

  error_flag = Create_Turbulence.Generate_Velocity_Fluctuations(Data.Initial_Mesh, Solution_Data.Input.Grid_IP);

  return error_flag;

}

/********************************************************
 * Routine: Hexa_Post_Processing_Specializations        *
 ********************************************************/
template<>
int Hexa_Post_Processing_Specializations(HexaSolver_Data &Data,
		 	                 HexaSolver_Solution_Data<LES3DFsd_pState, 
                                                                  LES3DFsd_cState>&Solution_Data) {

   int error_flag(0);
   double u_ave, v_ave, w_ave, sqr_u;
   
   Time_Averaging_of_Velocity_Field(Solution_Data.Local_Solution_Blocks.Soln_Blks,
                                    Data.Local_Adaptive_Block_List,
                                    Solution_Data.Input.Grid_IP,
                                    u_ave,
                                    v_ave,
                                    w_ave);

   Time_Averaging_of_Solution(Solution_Data.Local_Solution_Blocks.Soln_Blks,
                              Data.Local_Adaptive_Block_List,
                              Solution_Data.Input.Grid_IP,
                              u_ave,
                              v_ave,
                              w_ave,
                              sqr_u);

   Time_Averaging_of_Turbulent_Burning_Rate(Solution_Data.Local_Solution_Blocks.Soln_Blks,
                                            Data.Local_Adaptive_Block_List,
                                            Solution_Data.Input.Grid_IP);

   return error_flag;

}

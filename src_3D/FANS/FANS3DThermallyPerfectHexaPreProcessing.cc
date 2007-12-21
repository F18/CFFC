
#ifndef _FANS3D_THERMALLYPERFECT_INCLUDED
#include "FANS3DThermallyPerfect.h"
#endif // _FANS3D_THERMALLYPERFECT_INCLUDED

/********************************************************
 * Routine: Pre_Processing_Specializations              *
 ********************************************************/
template<>
int Hexa_Pre_Processing_Specializations(HexaSolver_Data &Data,
                                        HexaSolver_Solution_Data<FANS3D_ThermallyPerfect_KOmega_pState, 
                                                                 FANS3D_ThermallyPerfect_KOmega_cState> &Solution_Data) {

   int error_flag(0);
  
   // Initialization of the flow field with available 2D numerical solution
   // for some cases.
 
   if ((Solution_Data.Local_Solution_Blocks.Soln_Blks[0].Flow_Type != FLOWTYPE_TURBULENT_RANS_K_OMEGA)) return 0;
   
   if  ((Solution_Data.Input.i_ICs != IC_TURBULENT_DIFFUSION_FLAME  && Solution_Data.Input.i_ICs != IC_TURBULENT_COFLOW)) return (0);
   
   FlowField_2D Numflowfield2D;

   Numflowfield2D.Allocate(Solution_Data.Input.i_ICs);

   if (Solution_Data.Input.i_ICs == IC_TURBULENT_DIFFUSION_FLAME){
      
      if(CFFC_Primary_MPI_Processor()) Numflowfield2D.read_numerical_solution_hotflow(Solution_Data.Input.CFFC_Path);
   }


   if (Solution_Data.Input.i_ICs == IC_TURBULENT_COFLOW){
      if(CFFC_Primary_MPI_Processor()) Numflowfield2D.read_numerical_solution_coldflow(Solution_Data.Input.CFFC_Path);
   }
   
#ifdef _MPI_VERSION  
   MPI::COMM_WORLD.Bcast(Numflowfield2D.data.begin(), Numflowfield2D.data.size(), MPI::DOUBLE, 0);
#endif
   
   Solution_Data.Local_Solution_Blocks.Interpolate_2Dto3D(Numflowfield2D);
     

   Numflowfield2D.Deallocate(Solution_Data.Input.i_ICs);
   
   return error_flag;
   
}

/********************************************************
 * Routine: Hexa_Post_Processing_Specializations        *
 ********************************************************/
template<>
int Hexa_Post_Processing_Specializations(HexaSolver_Data &Data,
		 	                 HexaSolver_Solution_Data<FANS3D_ThermallyPerfect_KOmega_pState, 
                                                                  FANS3D_ThermallyPerfect_KOmega_cState>&Solution_Data) {

   int error_flag(0);

   return error_flag;

}

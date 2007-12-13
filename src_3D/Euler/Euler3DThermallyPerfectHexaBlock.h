/* EulerThermallyPerfectHexaBlock.h 

   Header file defining 3D Hexahedral mesh solution classes
   for the Euler equations for a thermally perfect
   gaseous mixture. */

#ifndef _EULER3D_THERMALLYPERFECT_HEXA_BLOCK_INCLUDED
#define _EULER3D_THERMALLYPERFECT_HEXA_BLOCK_INCLUDED

/* Include required CFFC header files. */

#ifndef _HEXA_BLOCK_INCLUDED
#include "../HexaBlock/HexaBlock.h"
#endif // HEXA_BLOCK_INCLUDED

#ifndef _EULER3D_THERMALLYPERFECT_STATE_INCLUDED
#include "Euler3DThermallyPerfectState.h"
#endif // _EULER3D_THERMALLYPERFECT_STATE_INCLUDED   

/* Define required specializations. */

template<>
void Hexa_Block<Euler3D_ThermallyPerfect_pState,
	        Euler3D_ThermallyPerfect_cState>::
Output_Tecplot(Input_Parameters<Euler3D_ThermallyPerfect_pState,
      		                Euler3D_ThermallyPerfect_cState> &IPs,
               const int Number_of_Time_Steps,
	       const double &Time,  
               const int Block_Number,
	       const int Output_Title,
	       ostream &Out_File);
				
template<>
void Hexa_Block<Euler3D_ThermallyPerfect_pState, 
                Euler3D_ThermallyPerfect_cState>::
Output_Cells_Tecplot(Input_Parameters<Euler3D_ThermallyPerfect_pState,
      		                      Euler3D_ThermallyPerfect_cState> &IPs, 
                     const int Number_of_Time_Steps,
                     const double &Time,
                     const int Block_Number,
                     const int Output_Title,
                     ostream &Out_File);

template<>
void Hexa_Block<Euler3D_ThermallyPerfect_pState, 
                Euler3D_ThermallyPerfect_cState>::
Output_Nodes_Tecplot(Input_Parameters<Euler3D_ThermallyPerfect_pState,
      		                      Euler3D_ThermallyPerfect_cState> &IPs, 
                     const int Number_of_Time_Steps,
                     const double &Time,
                     const int Block_Number,
                     const int Output_Title,
                     ostream &Out_File);

template<>
int Hexa_Block<Euler3D_ThermallyPerfect_pState, 
               Euler3D_ThermallyPerfect_cState>::
Update_Solution_Multistage_Explicit(const int i_stage,
                                    Input_Parameters<Euler3D_ThermallyPerfect_pState,
                                                     Euler3D_ThermallyPerfect_cState> &IPs);

#endif // _EULER3D_THERMALLYPERFECT_HEXA_BLOCK_INCLUDED

/* NavierStokesThermallyPerfectHexaBlock.h 

   Header file defining 3D Hexahedral mesh solution classes
   for the Navier-Stokes equations for a thermally perfect
   gaseous mixture. */

#ifndef _NAVIERSTOKES3D_THERMALLYPERFECT_HEXA_BLOCK_INCLUDED
#define _NAVIERSTOKES3D_THERMALLYPERFECT_HEXA_BLOCK_INCLUDED

/* Include required CFFC header files. */

#ifndef _HEXA_BLOCK_INCLUDED
#include "../HexaBlock/HexaBlock.h"
#endif // HEXA_BLOCK_INCLUDED

#ifndef _NAVIERSTOKES3D_THERMALLYPERFECT_STATE_INCLUDED
#include "NavierStokes3DThermallyPerfectState.h"
#endif // _NAVIERSTOKES3D_THERMALLYPERFECT_STATE_INCLUDED   

/* Define required specializations. */

template<>
void Hexa_Block<NavierStokes3D_ThermallyPerfect_pState,
	        NavierStokes3D_ThermallyPerfect_cState>::
Output_Tecplot(Input_Parameters<NavierStokes3D_ThermallyPerfect_pState,
      		                NavierStokes3D_ThermallyPerfect_cState> &IPs,
               const int Number_of_Time_Steps,
	       const double &Time,  
               const int Block_Number,
	       const int Output_Title,
	       ostream &Out_File);
				
template<>
void Hexa_Block<NavierStokes3D_ThermallyPerfect_pState, 
                NavierStokes3D_ThermallyPerfect_cState>::
Output_Cells_Tecplot(Input_Parameters<NavierStokes3D_ThermallyPerfect_pState,
      		                      NavierStokes3D_ThermallyPerfect_cState> &IPs, 
                     const int Number_of_Time_Steps,
                     const double &Time,
                     const int Block_Number,
                     const int Output_Title,
                     ostream &Out_File);

template<>
void Hexa_Block<NavierStokes3D_ThermallyPerfect_pState, 
                NavierStokes3D_ThermallyPerfect_cState>::
Output_Nodes_Tecplot(Input_Parameters<NavierStokes3D_ThermallyPerfect_pState,
      		                      NavierStokes3D_ThermallyPerfect_cState> &IPs, 
                     const int Number_of_Time_Steps,
                     const double &Time,
                     const int Block_Number,
                     const int Output_Title,
                     ostream &Out_File);

template<>
int Hexa_Block<NavierStokes3D_ThermallyPerfect_pState, 
               NavierStokes3D_ThermallyPerfect_cState>::
ICs(Input_Parameters<NavierStokes3D_ThermallyPerfect_pState, 
                     NavierStokes3D_ThermallyPerfect_cState> &IPs);

template<>
double Hexa_Block<NavierStokes3D_ThermallyPerfect_pState,
                  NavierStokes3D_ThermallyPerfect_cState>::
CFL(Input_Parameters<NavierStokes3D_ThermallyPerfect_pState,
                     NavierStokes3D_ThermallyPerfect_cState> &IPs);

template<>
int Hexa_Block<NavierStokes3D_ThermallyPerfect_pState, 
               NavierStokes3D_ThermallyPerfect_cState>::
dUdt_Multistage_Explicit(const int i_stage,
                         Input_Parameters<NavierStokes3D_ThermallyPerfect_pState, 
                                          NavierStokes3D_ThermallyPerfect_cState> &IPs);

template<>
int Hexa_Block<NavierStokes3D_ThermallyPerfect_pState, 
               NavierStokes3D_ThermallyPerfect_cState>::
Update_Solution_Multistage_Explicit(const int i_stage,
                                    Input_Parameters<NavierStokes3D_ThermallyPerfect_pState,
                                                     NavierStokes3D_ThermallyPerfect_cState> &IPs);

#endif // _NAVIERSTOKES3D_THERMALLYPERFECT_HEXA_BLOCK_INCLUDED

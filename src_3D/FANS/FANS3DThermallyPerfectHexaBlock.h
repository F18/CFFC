/* FANS3DHexaBlock.h:  Header file defining 
                       3D Hexahedral mesh solution classes
                       for the Favre-Averaged Navier-Stokes 
                       equations. */

#ifndef _FANS3D_THERMALLYPERFECT_HEXA_BLOCK_INCLUDED
#define _FANS3D_THERMALLYPERFECT_HEXA_BLOCK_INCLUDED

/* Include required CFFC header files. */

#ifndef _HEXA_BLOCK_INCLUDED
#include "../HexaBlock/HexaBlock.h"
#endif // _HEXA_BLOCK_INCLUDED

#ifndef _FANS3D_THERMALLYPERFECT_STATE_INCLUDED
#include "FANS3DThermallyPerfectState.h"
#endif // _FANS3D_THERMALLYPERFECT_STATE_INCLUDED   

#ifndef _BLUFFBODY_DATABASE_INCLUDED
#include "BluffBodyBurner.h"
#endif // _BLUFFBODY_DATABASE_INCLUDE

/* Include 2D to 3D solution interpolation header file. */

#ifndef _INTERPOLATION2DTO3D_INCLUDED
#include "Interpolation2Dto3D.h"
#endif// _INTERPOLATION2DTO3D_INCLUDED

/* Define required specializations. */

template<>
int Hexa_Block<FANS3D_ThermallyPerfect_KOmega_pState,
	        FANS3D_ThermallyPerfect_KOmega_cState>::
Update_Corner_Cells_for_3_Blks_Abutting(const int i_elem, 
                                        const int j_elem, 
                                        const int k_elem, 
                                        const int numNeigh,
                                        const int be);

template<>
void Hexa_Block<FANS3D_ThermallyPerfect_KOmega_pState,
	        FANS3D_ThermallyPerfect_KOmega_cState>::
Output_Tecplot(Input_Parameters<FANS3D_ThermallyPerfect_KOmega_pState,
      		                FANS3D_ThermallyPerfect_KOmega_cState> &IPs,
               const int Number_of_Time_Steps,
	       const double &Time,  
               const int Block_Number,
	       const int Output_Title,
	       ostream &Out_File);
				
template<>
void Hexa_Block<FANS3D_ThermallyPerfect_KOmega_pState, 
                FANS3D_ThermallyPerfect_KOmega_cState>::
Output_Cells_Tecplot(Input_Parameters<FANS3D_ThermallyPerfect_KOmega_pState,
      		                      FANS3D_ThermallyPerfect_KOmega_cState> &IPs, 
                     const int Number_of_Time_Steps,
                     const double &Time,
                     const int Block_Number,
                     const int Output_Title,
                     ostream &Out_File);

template<>
void Hexa_Block<FANS3D_ThermallyPerfect_KOmega_pState, 
                FANS3D_ThermallyPerfect_KOmega_cState>::
Output_Nodes_Tecplot(Input_Parameters<FANS3D_ThermallyPerfect_KOmega_pState,
      		                      FANS3D_ThermallyPerfect_KOmega_cState> &IPs, 
                     const int Number_of_Time_Steps,
                     const double &Time,
                     const int Block_Number,
                     const int Output_Title,
                     ostream &Out_File);

template<>
int Hexa_Block<FANS3D_ThermallyPerfect_KOmega_pState, 
               FANS3D_ThermallyPerfect_KOmega_cState>::
ICs(Input_Parameters<FANS3D_ThermallyPerfect_KOmega_pState, 
                     FANS3D_ThermallyPerfect_KOmega_cState> &IPs);

template<>
int Hexa_Block<FANS3D_ThermallyPerfect_KOmega_pState,
               FANS3D_ThermallyPerfect_KOmega_cState>::
Interpolate_2Dto3D(const FlowField_2D &Numflowfield2D);

template<>
void Hexa_Block<FANS3D_ThermallyPerfect_KOmega_pState, 
                FANS3D_ThermallyPerfect_KOmega_cState>::
BCs(Input_Parameters<FANS3D_ThermallyPerfect_KOmega_pState, 
                     FANS3D_ThermallyPerfect_KOmega_cState> &IPs);

template<>
double Hexa_Block<FANS3D_ThermallyPerfect_KOmega_pState,
                  FANS3D_ThermallyPerfect_KOmega_cState>::
CFL(Input_Parameters<FANS3D_ThermallyPerfect_KOmega_pState,
                     FANS3D_ThermallyPerfect_KOmega_cState> &IPs);

template<>
int Hexa_Block<FANS3D_ThermallyPerfect_KOmega_pState,
               FANS3D_ThermallyPerfect_KOmega_cState>::
dUdt_Multistage_Explicit(const int i_stage,
                         Input_Parameters<FANS3D_ThermallyPerfect_KOmega_pState,
                                          FANS3D_ThermallyPerfect_KOmega_cState> &IPs);

template<>
int Hexa_Block<FANS3D_ThermallyPerfect_KOmega_pState, 
               FANS3D_ThermallyPerfect_KOmega_cState>::
Update_Solution_Multistage_Explicit(const int i_stage,
                                    Input_Parameters<FANS3D_ThermallyPerfect_KOmega_pState,
                                                     FANS3D_ThermallyPerfect_KOmega_cState> &IPs);

template<>
double Hexa_Block<FANS3D_ThermallyPerfect_KOmega_pState, 
                  FANS3D_ThermallyPerfect_KOmega_cState>::
Wall_Friction_Velocity(const int i, 
                       const int j, 
                       const int k);

template<> 
int Hexa_Block<FANS3D_ThermallyPerfect_KOmega_pState,  
               FANS3D_ThermallyPerfect_KOmega_cState>::Wall_Shear(void);

#endif // _FANS3D_THERMALLYPERFECT_HEXA_BLOCK_INCLUDED

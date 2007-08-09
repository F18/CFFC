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

//specilizations 
/********************************************************
 * Routine: dUdt_Multistage_Explicit                    *
 *                                                      *
 * This routine determines the solution residuals for a *
 * given stage of a variety of multi-stage explicit     *
 * time integration schemes for a given solution block. *
 *                                                      *
 ********************************************************/

template<>
int Hexa_Block<FANS3D_ThermallyPerfect_KOmega_pState,
               FANS3D_ThermallyPerfect_KOmega_cState>::dUdt_Multistage_Explicit(
                  const int i_stage,
                  Input_Parameters<FANS3D_ThermallyPerfect_KOmega_pState,
                  FANS3D_ThermallyPerfect_KOmega_cState> &IPs);

template<>
double Hexa_Block<FANS3D_ThermallyPerfect_KOmega_pState,
                  FANS3D_ThermallyPerfect_KOmega_cState>::CFL(
                     Input_Parameters<FANS3D_ThermallyPerfect_KOmega_pState,
                     FANS3D_ThermallyPerfect_KOmega_cState> &IPs);

template<>
double Hexa_Block<FANS3D_ThermallyPerfect_KOmega_pState, 
                     FANS3D_ThermallyPerfect_KOmega_cState>::Wall_Friction_Velocity(
                        int i, int j, int k);

template< >
void Hexa_Block<FANS3D_ThermallyPerfect_KOmega_pState,
	        FANS3D_ThermallyPerfect_KOmega_cState>::
Output_Tecplot(Input_Parameters<FANS3D_ThermallyPerfect_KOmega_pState,
      		                FANS3D_ThermallyPerfect_KOmega_cState> &IPs,
               const int Number_of_Time_Steps,
	       const double &Time,  
               const int Block_Number,
	       const int Output_Title,
	       ostream &Out_File);
				
template< >
void Hexa_Block<FANS3D_ThermallyPerfect_KOmega_pState, 
                FANS3D_ThermallyPerfect_KOmega_cState>::
Output_Cells_Tecplot(Input_Parameters<FANS3D_ThermallyPerfect_KOmega_pState,
      		                      FANS3D_ThermallyPerfect_KOmega_cState> &IPs, 
                     const int Number_of_Time_Steps,
                     const double &Time,
                     const int Block_Number,
                     const int Output_Title,
                     ostream &Out_File);

template< >
void Hexa_Block<FANS3D_ThermallyPerfect_KOmega_pState, 
                FANS3D_ThermallyPerfect_KOmega_cState>::
Output_Nodes_Tecplot(Input_Parameters<FANS3D_ThermallyPerfect_KOmega_pState,
      		                      FANS3D_ThermallyPerfect_KOmega_cState> &IPs, 
                     const int Number_of_Time_Steps,
                     const double &Time,
                     const int Block_Number,
                     const int Output_Title,
                     ostream &Out_File);

template< >
int Hexa_Block<FANS3D_ThermallyPerfect_KOmega_pState, 
               FANS3D_ThermallyPerfect_KOmega_cState>::
ICs(const int i_ICtype,
    Input_Parameters<FANS3D_ThermallyPerfect_KOmega_pState, 
    FANS3D_ThermallyPerfect_KOmega_cState> &IPs);

template< >
void Hexa_Block<FANS3D_ThermallyPerfect_KOmega_pState, 
                FANS3D_ThermallyPerfect_KOmega_cState>::BCs(
                   Input_Parameters<FANS3D_ThermallyPerfect_KOmega_pState, 
                   FANS3D_ThermallyPerfect_KOmega_cState> &IPs);
template< > 
int Hexa_Block<FANS3D_ThermallyPerfect_KOmega_pState,  
               FANS3D_ThermallyPerfect_KOmega_cState>::Wall_Shear(void);

template<>
int Hexa_Block<FANS3D_ThermallyPerfect_KOmega_pState, 
               FANS3D_ThermallyPerfect_KOmega_cState>::
Update_Solution_Multistage_Explicit(
   const int i_stage,
   Input_Parameters<FANS3D_ThermallyPerfect_KOmega_pState,
   FANS3D_ThermallyPerfect_KOmega_cState> &IPs);


#endif // _FANS3D_THERMALLYPERFECT_HEXA_BLOCK_INCLUDED

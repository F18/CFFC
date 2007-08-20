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
#endif // NAVIERSTOKES3D_THERMALLYPERFECT_STATE_INCLUDED   

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
int Hexa_Block<NavierStokes3D_ThermallyPerfect_pState, 
               NavierStokes3D_ThermallyPerfect_cState>::
dUdt_Multistage_Explicit(const int i_stage,
                         Input_Parameters<NavierStokes3D_ThermallyPerfect_pState, 
                                          NavierStokes3D_ThermallyPerfect_cState> &IPs);


#endif // _NAVIERSTOKES3D_THERMALLYPERFECT_HEXA_BLOCK_INCLUDED

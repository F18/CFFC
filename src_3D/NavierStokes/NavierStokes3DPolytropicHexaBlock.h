/* NavierStokesPolytropicHexaBlock.h 

   Header file defining 3D Hexahedral mesh solution classes
   for the Navier-Stokes equations for a Polytropic
   gaseous mixture. */

#ifndef _NAVIERSTOKES3D_POLYTROPIC_HEXA_BLOCK_INCLUDED
#define _NAVIERSTOKES3D_POLYTROPIC_HEXA_BLOCK_INCLUDED

/* Include required CFFC header files. */

#ifndef _HEXA_BLOCK_INCLUDED
#include "../HexaBlock/HexaBlock.h"
#endif // HEXA_BLOCK_INCLUDED

#ifndef _NAVIERSTOKES3D_POLYTROPIC_STATE_INCLUDED
#include "NavierStokes3DPolytropicState.h"
#endif // _NAVIERSTOKES3D_POLYTROPIC_STATE_INCLUDED   

/* Define required specializations. */

template<>
int Hexa_Block<NavierStokes3D_Polytropic_pState, 
               NavierStokes3D_Polytropic_cState>::
ICs(Input_Parameters<NavierStokes3D_Polytropic_pState, 
                     NavierStokes3D_Polytropic_cState> &IPs);

template<>
double Hexa_Block<NavierStokes3D_Polytropic_pState,
                  NavierStokes3D_Polytropic_cState>::
CFL(Input_Parameters<NavierStokes3D_Polytropic_pState,
                     NavierStokes3D_Polytropic_cState> &IPs);

template<>
int Hexa_Block<NavierStokes3D_Polytropic_pState, 
               NavierStokes3D_Polytropic_cState>::
dUdt_Multistage_Explicit(const int i_stage,
                         Input_Parameters<NavierStokes3D_Polytropic_pState, 
                                          NavierStokes3D_Polytropic_cState> &IPs);

template<>
int Hexa_Block<NavierStokes3D_Polytropic_pState, 
               NavierStokes3D_Polytropic_cState>::
Update_Solution_Multistage_Explicit(const int i_stage,
                                    Input_Parameters<NavierStokes3D_Polytropic_pState,
                                                     NavierStokes3D_Polytropic_cState> &IPs);

#endif // _NAVIERSTOKES3D_POLYTROPIC_HEXA_BLOCK_INCLUDED

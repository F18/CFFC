/* NavierStokes3DPolytropic.h:  

   Header file for the solution of the Navier-Stokes equations 
   governing polytropic non-reactive and combusting
   gaseous mixtures. */

#ifndef _NAVIERSTOKES3D_POLYTROPIC_INCLUDED
#define _NAVIERSTOKES3D_POLYTROPIC_INCLUDED

/* Include required CFFC header files. */

#ifndef _HEXA_SOLVER_INCLUDED
#include "../HexaBlock/HexaSolver.h"
#endif // _HEXA_SOLVER_INCLUDED

#ifndef _NAVIERSTOKES3D_POLYTROPIC_STATE_INCLUDED
#include "NavierStokes3DPolytropicState.h"
#endif // _NAVIERSTOKES3D_POLYTROPIC_STATE_INCLUDED

#ifndef _NAVIERSTOKES3D_POLYTROPIC_INPUT_INCLUDED
#include "NavierStokes3DPolytropicInput.h"
#endif // _NAVIERSTOKES3D_POLYTROPIC_INPUT_INCLUDED

#ifndef _NAVIERSTOKES3D_POLYTROPIC_HEXA_BLOCK_INCLUDED
#include "NavierStokes3DPolytropicHexaBlock.h"
#endif // _NAVIERSTOKES3D_POLYTROPIC_HEXA_BLOCK_INCLUDED

#ifndef _MPI_INCLUDED
#include "../MPI/MPI.h"
#endif // _MPI_INCLUDED

#ifndef _ICEMCFD_INCLUDED
#include "../ICEM/ICEMCFD.h"
#endif // _ICEMCFD_INCLUDED

#endif // _NAVIERSTOKES3D_POLYTROPIC_INCLUDED


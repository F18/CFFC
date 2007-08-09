/* NavierStokes3DThermallyPerfect.h:  

   Header file for the 3D Navier-Stokes equations 
   for thermally perfect gaseous mixtures. */

#ifndef _NAVIERSTOKES3D_THERMALLYPERFECT_INCLUDED
#define _NAVIERSTOKES3D_THERMALLYPERFECT_INCLUDED

/* Include required CFFC header files. */

#ifndef _HEXA_SOLVER_INCLUDED
#include "../HexaBlock/HexaSolver.h"
#endif // _HEXA_SOLVER_INCLUDED

#ifndef _NAVIERSTOKES3D_THERMALLYPERFECT_STATE_INCLUDED
#include "NavierStokes3DThermallyPerfectState.h"
#endif // _NAVIERSTOKES3D_THERMALLYPERFECT_STATE_INCLUDED

#ifndef _NAVIERSTOKES3_THERMALLYPERFECT__HEXA_BLOCK_INCLUDED
#include "NavierStokes3DThermallyPerfectHexaBlock.h"
#endif // _NAVIERSTOKES3D_THERMALLYPERFECT_HEXA_BLOCK_INCLUDED

#ifndef _MPI_INCLUDED
#include "../MPI/MPI.h"
#endif // _MPI_INCLUDED

#ifndef _ICEMCFD_INCLUDED
#include "../ICEM/ICEMCFD.h"
#endif // _ICEMCFD_INCLUDED

#endif // _NAVIERSTOKES3D_THERMALLYPERFECT_INCLUDED


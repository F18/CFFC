/* FANS3DThermallyPerfect.h:  

   Header file for the 3D Favre-Averaged Navier-Stokes (FANS) equations 
   for thermally perfect gaseous mixtures. */

#ifndef _FANS3D_THERMALLYPERFECT_INCLUDED
#define _FANS3D_THERMALLYPERFECT_INCLUDED

/* Include required CFFC header files. */

#ifndef _HEXA_SOLVER_INCLUDED
#include "../HexaBlock/HexaSolver.h"
#endif // _HEXA_SOLVER_INCLUDED

#ifndef _FANS3D_THERMALLYPERFECT_STATE_INCLUDED
#include "FANS3DThermallyPerfectState.h"
#endif // _FANS3D_THERMALLYPERFECT_STATE_INCLUDED

#ifndef _FANS3D_THERMALLYPERFECT_HEXA_BLOCK_INCLUDED
#include "FANS3DThermallyPerfectHexaBlock.h"
#endif // _FANS3D_THERMALLYPERFECT_HEXA_BLOCK_INCLUDED

#ifndef _FANS3D_THERMALLYPERFECT_HEXA_MULTIBLOCK_INCLUDED
#include "FANS3DThermallyPerfectHexaMultiBlock.h"
#endif // _FANS3D_THERMALLYPERFECT_HEXA_MULTIBLOCK_INCLUDED

#ifndef _MPI_INCLUDED
#include "../MPI/MPI.h"
#endif // _MPI_INCLUDED

#ifndef _ICEMCFD_INCLUDED
#include "../ICEM/ICEMCFD.h"
#endif // _ICEMCFD_INCLUDED

#endif // _FANS3D_THERMALLYPERFECT_INCLUDED


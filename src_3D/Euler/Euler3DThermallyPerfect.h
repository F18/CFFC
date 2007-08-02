/* Euler3DThermallyPerfect.h:  

   Header file for the 3D Euler equations 
   for thermally perfect gaseous mixtures. */

#ifndef _EULER3D_THERMALLYPERFECT_INCLUDED
#define _EULER3D_THERMALLYPERFECT_INCLUDED

/* Include required CFFC header files. */

#ifndef _HEXA_SOLVER_INCLUDED
#include "../HexaBlock/HexaSolver.h"
#endif // _HEXA_SOLVER_INCLUDED

#ifndef _EULER3D_THERMALLYPERFECT_STATE_INCLUDED
#include "Euler3DThermallyPerfectState.h"
#endif // _EULER3D_THERMALLYPERFECT_STATE_INCLUDED

#ifndef _MPI_INCLUDED
#include "../MPI/MPI.h"
#endif // _MPI_INCLUDED

#ifndef _ICEMCFD_INCLUDED
#include "../ICEM/ICEMCFD.h"
#endif // _ICEMCFD_INCLUDED

#endif // _EULER3D_THERMALLYPERFECT_INCLUDED


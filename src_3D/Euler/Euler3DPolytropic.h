/* Euler3DPolytropic.h:  

   Header file for the solution of the Euler equations 
   governing polytropic (thermally and calorically 
   perfect) gases. */

#ifndef _EULER3D_POLYTROPIC_INCLUDED
#define _EULER3D_POLYTROPIC_INCLUDED

/* Include required CFFC header files. */

#ifndef _HEXA_SOLVER_INCLUDED
#include "../HexaBlock/HexaSolver.h"
#endif // _HEXA_SOLVER_INCLUDED

#ifndef _EULER3D_POLYTROPIC_STATE_INCLUDED
#include "Euler3DPolytropicState.h"
#endif // _EULER3D_POLYTROPIC_STATE_INCLUDED

#ifndef _EULER3D_POLYTROPIC_INPUT_INCLUDED
#include "Euler3DPolytropicInput.h"
#endif // _EULER3D_POLYTROPIC_INPUT_INCLUDED

#ifndef _EULER3D_POLYTROPIC_NKS_INLCUDED
#include "Euler3DPolytropicNKS.h"
#endif //_EULER3D_POLYTROPIC_NKS_INLCUDED

#ifndef _MPI_INCLUDED
#include "../MPI/MPI.h"
#endif // _MPI_INCLUDED

#ifndef _ICEMCFD_INCLUDED
#include "../ICEM/ICEMCFD.h"
#endif // _ICEMCFD_INCLUDED

#endif // _EULER3D_POLYTROPIC_INCLUDED


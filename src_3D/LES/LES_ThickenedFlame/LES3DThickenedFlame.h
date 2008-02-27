/* LES3DThickenedFlame.h:  

   Header file for the 3D LES Navier-Stokes (LES) equations 
   governing thermally perfect reactive combusting gaseous 
   mixtures. */

#ifndef _LES3DTF_INCLUDED
#define _LES3DTF_INCLUDED

/* Include required CFFC header files. */

#ifndef _HEXA_SOLVER_INCLUDED
#include "../../HexaBlock/HexaSolver.h"
#endif // _HEXA_SOLVER_INCLUDED

#ifndef _LES3DTF_STATE_INCLUDED
#include "LES3DThickenedFlameState.h"
#endif // _LES3DTF_STATE_INCLUDED

#ifndef _LES3DTF_INPUT_INCLUDED
#include "LES3DThickenedFlameInput.h"
#endif // _LES3DTF_INPUT_INCLUDED

#ifndef _LES3DTF_HEXA_BLOCK_INCLUDED
#include "LES3DThickenedFlameHexaBlock.h"
#endif // _LES3DTF_HEXA_BLOCK_INCLUDED

#ifndef _LES3DTF_HEXA_PREPROCESSING_INCLUDED
#include "LES3DThickenedFlameHexaPreProcessing.h"
#endif // _LES3DTF_HEXA_PREPROCESSING_INCLUDED

#ifndef _MPI_INCLUDED
#include "../../MPI/MPI.h"
#endif // _MPI_INCLUDED

#ifndef _ICEMCFD_INCLUDED
#include "../../ICEM/ICEMCFD.h"
#endif // _ICEMCFD_INCLUDED

#endif // _LES3DTF_INCLUDED


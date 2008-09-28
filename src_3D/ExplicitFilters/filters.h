/* filters.h:  
 
 Header file for the solution of the Euler equations 
 governing polytropic (thermally and calorically 
 perfect) gases. */

#ifndef _FILTERS_INCLUDED
#define _FILTERS_INCLUDED

/* Include required CFFC header files. */

#ifndef _CFD_INCLUDED
#include "../CFD/CFD.h"
#endif // _CFD_INCLUDED

#ifndef _INPUT_INCLUDED
#include "../CFD/Input.h"
#endif // _INPUT_INCLUDED

#ifndef _HEXA_SOLVER_INCLUDED
#include "../HexaBlock/HexaSolver.h"
#endif // _HEXA_SOLVER_INCLUDED

#ifndef _EULER3D_POLYTROPIC_INCLUDED
#include "../Euler/Euler3DPolytropic.h"
#endif // _EULER3D_POLYTROPIC_STATE_INCLUDED

//#ifndef _EULER3D_POLYTROPIC_INPUT_INCLUDED
//#include "Euler3DPolytropicInput.h"
//#endif // _EULER3D_POLYTROPIC_INPUT_INCLUDED

//#ifndef _EULER3D_POLYTROPIC_NKS_INLCUDED
//#include "Euler3DPolytropicNKS.h"
//#endif //_EULER3D_POLYTROPIC_NKS_INLCUDED

#ifndef _MPI_INCLUDED
#include "../MPI/MPI.h"
#endif // _MPI_INCLUDED

#ifndef _ICEMCFD_INCLUDED
#include "../ICEM/ICEMCFD.h"
#endif // _ICEMCFD_INCLUDED

#ifndef _SPECTRAL_ANALYSIS_INCLUDED
#include "../TurbulenceModelling/SpectralAnalysis.h"
#endif // _SPECTRAL_ANALYSIS_INCLUDED

#ifndef _FILTER_CONTROLLER_INCLUDED
#include "Filter_Controller.h"
#endif // _FILTER_CONTROLLER_INCLUDED

#ifndef _CFFC_FILTER_CONTROLLER_INCLUDED
#include "CFFC_Filter_Controller.h"
#endif // _FILTER_CONTROLLER_INCLUDED

#ifndef _EXPLICIT_FILTER_INCLUDED
#include "Explicit_Filter.h"
#endif

#ifndef _EXLICIT_FILTER_TRANSLATOR
#include "Explicit_Filter_Translator.h"
#endif

#endif // _FILTERS_INCLUDED


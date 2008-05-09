/**************LES3DPolytropicHexaPreProcessing.h*************/

#ifndef _LES3D_POLYTROPIC_HEXA_PREPROCESSING_INCLUDED
#define _LES3D_POLYTROPIC_HEXA_PREPROCESSING_INCLUDED 

#ifndef _LES3D_POLYTROPIC_INCLUDED
#include "LES3DPolytropic.h"
#endif // _LES3D_POLYTROPIC_INCLUDED

/* Define the specializations. */

//! Pre_Processing_Specializations
template<>
int Hexa_Pre_Processing_Specializations(HexaSolver_Data &Data,
                                        HexaSolver_Solution_Data<LES3D_Polytropic_pState,
                                                                 LES3D_Polytropic_cState> &Solution_Data);

//! Hexa_Post_Processing_Specializations
template<>
int Hexa_Post_Processing_Specializations(HexaSolver_Data &Data,
                                         HexaSolver_Solution_Data<LES3D_Polytropic_pState,
                                                                  LES3D_Polytropic_cState> &Solution_Data);

//! Open_Other_Solution_Progress_Specialization_Files
template<>
int Open_Other_Solution_Progress_Specialization_Files(HexaSolver_Data &Data,
                                                      HexaSolver_Solution_Data<LES3D_Polytropic_pState,
                                                                               LES3D_Polytropic_cState> &Solution_Data);

//! Close_Other_Solution_Progress_Specialization_Files
template<>
int Close_Other_Solution_Progress_Specialization_Files(HexaSolver_Data &Data,
                                                       HexaSolver_Solution_Data<LES3D_Polytropic_pState,
                                                                                LES3D_Polytropic_cState> &Solution_Data);

//! Output_Other_Solution_Progress_Specialization_Data
template<>
int Output_Other_Solution_Progress_Specialization_Data(HexaSolver_Data &Data,
                                                       HexaSolver_Solution_Data<LES3D_Polytropic_pState,
                                                                                LES3D_Polytropic_cState> &Solution_Data);

#endif // _LES3D_POLYTROPIC_HEXA_PREPROCESSING_INCLUDED

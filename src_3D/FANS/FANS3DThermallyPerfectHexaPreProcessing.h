/**************FANS3DThermallyPerfectHexaPreProcessing.h*************/

#ifndef _FANS3D_THERMALLYPERFECT_HEXA_PREPROCESSING_INCLUDED
#define _FANS3D_THERMALLYPERFECT_HEXA_PREPROCESSING_INCLUDED

/* Define the specializations. */

//! Pre_Processing_Specializations
template<>
int Hexa_Pre_Processing_Specializations(HexaSolver_Data &Data,
                                        HexaSolver_Solution_Data<FANS3D_ThermallyPerfect_KOmega_pState,
                                                                 FANS3D_ThermallyPerfect_KOmega_cState> &Solution_Data);

//! Hexa_Post_Processing_Specializations
template<>
int Hexa_Post_Processing_Specializations(HexaSolver_Data &Data,
		 	                 HexaSolver_Solution_Data<FANS3D_ThermallyPerfect_KOmega_pState,
                                                                  FANS3D_ThermallyPerfect_KOmega_cState> &Solution_Data);

#endif // _FANS3D_THERMALLYPERFECT_HEXA_PREPROCESSING_INCLUDED

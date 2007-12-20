/**************LES3DFsdHexaPreProcessing.h*************/

#ifndef _LES3DFSD_HEXA_PREPROCESSING_INCLUDED
#define _LES3DFSD_HEXA_PREPROCESSING_INCLUDED

/* Define the specializations. */

//! Pre_Processing_Specializations
template<>
int Hexa_Pre_Processing_Specializations(HexaSolver_Data &Data,
                                        HexaSolver_Solution_Data<LES3DFsd_pState,
                                                                 LES3DFsd_cState> &Solution_Data);

//! Hexa_Post_Processing_Specializations
template<>
int Hexa_Post_Processing_Specializations(HexaSolver_Data &Data,
		 	                 HexaSolver_Solution_Data<LES3DFsd_pState,
                                                                  LES3DFsd_cState> &Solution_Data);

#endif // _LES3DFSD_HEXA_PREPROCESSING_INCLUDED

/**************LES3DThickenedFlameHexaPreProcessing.h*************/

#ifndef _LES3DTF_HEXA_PREPROCESSING_INCLUDED
#define _LES3DTF_HEXA_PREPROCESSING_INCLUDED


/* Define the specializations. */


//! Initialize_Solution_Blocks_Specializations
template<>
int Initialize_Solution_Blocks_Specializations(HexaSolver_Data &Data,
					       HexaSolver_Solution_Data<LES3DTF_pState, 
					                                LES3DTF_cState> &Solution_Data);

//! Pre_Processing_Specializations
template<>
int Hexa_Pre_Processing_Specializations(HexaSolver_Data &Data,
                                        HexaSolver_Solution_Data<LES3DTF_pState,
                                                                 LES3DTF_cState> &Solution_Data);

//! BCs_Specializations
template<>
int Hexa_BCs_Specializations(HexaSolver_Data &Data,
			     HexaSolver_Solution_Data<LES3DTF_pState, 
			                              LES3DTF_cState> &Solution_Data);

//! Hexa_Post_Processing_Specializations
template<>
int Hexa_Post_Processing_Specializations(HexaSolver_Data &Data,
		 	                 HexaSolver_Solution_Data<LES3DTF_pState,
                                                                  LES3DTF_cState> &Solution_Data);
//! Open_Other_Solution_Progress_Specialization_Files
template<>
int Open_Other_Solution_Progress_Specialization_Files(HexaSolver_Data &Data,
                                                      HexaSolver_Solution_Data<LES3DTF_pState,
                                                                               LES3DTF_cState> &Solution_Data);

//! Close_Other_Solution_Progress_Specialization_Files
template<>
int Close_Other_Solution_Progress_Specialization_Files(HexaSolver_Data &Data,
                                                       HexaSolver_Solution_Data<LES3DTF_pState,
						                                LES3DTF_cState> &Solution_Data);

//! Output_Other_Solution_Progress_Specialization_Data
template<>
int Output_Other_Solution_Progress_Specialization_Data(HexaSolver_Data &Data,
                                                       HexaSolver_Solution_Data<LES3DTF_pState,
						                                LES3DTF_cState> &Solution_Data);

//! Burning rate
template<>
double Turbulent_Burning_Rate(Hexa_Block<LES3DTF_pState, LES3DTF_cState> *Solution_Block,
			      AdaptiveBlock3D_List &LocalSolnBlockList,
			      Grid3D_Input_Parameters &IPs);


#endif // _LES3DTF_HEXA_PREPROCESSING_INCLUDED

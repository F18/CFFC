#ifndef _ADVECTDIFFUSE2D_FASMULTIGRID_SPECIALIZATION_INCLUDED 
#define _ADVECTDIFFUSE2D_FASMULTIGRID_SPECIALIZATION_INCLUDED 

#include "New_AdvectDiffuse2DQuad.h"
#include "../FASMultigrid2D/FASMultigrid2D.h"

/*! *****************************************************************************************
 *  Specialization of Update_Primitive_Variables                                            *
 *                                                                                          *
 ********************************************************************************************/
template<> inline
void FAS_Multigrid2D_Solver<AdvectDiffuse2D_State_New,
			    AdvectDiffuse2D_Quad_Block_New,
			    AdvectDiffuse2D_Input_Parameters>::Update_Primitive_Variables(const int &Level) {
  // do nothing (there are no conserved and primitive variables in advection diffusion)
}



#endif	// _ADVECTDIFFUSE2D_FASMULTIGRID_SPECIALIZATION_INCLUDED 

#ifndef _ADVECTDIFFUSE2D_QUAD_NKS_INCLUDED
#define _ADVECTDIFFUSE2D_QUAD_NKS_INCLUDED

/* Include header files. */

#ifndef _ADVECTDIFFUSE2D_QUAD_INCLUDED
#include "AdvectDiffuse2DQuad.h"
#endif // _ADVECTDIFFUSE2D_QUAD_INCLUDED

#ifndef _ADVECTDIFFUSE2D_QUAD_GMRES_INCLUDED
#include "AdvectDiffuse2DQuad_GMRES.h"
#endif // _ADVECTDIFFUSE2D_QUAD_GMRES_INCLUDED

/********************************************************
 * AdvectDiffuse2DQuad_NKS -- External subroutines.     *
 ********************************************************/
extern CompRow_Mat_double Jacobian(AdvectDiffuse2D_Quad_Block &SolnBlk, 
  				   double tau, 
				   int JaS, 
				   int JdS, 
				   int JsS, 
				   int overlap);

extern double Ja(AdvectDiffuse2D_Quad_Block &SolnBlk, 
		 int ncol, 
		 int nrow, 
		 int JaS, 
		 int overlap);

extern double Jd(AdvectDiffuse2D_Quad_Block &SolnBlk, 
		 int ncol, 
		 int nrow, 
		 int JdS, 
		 int JaS, 
		 int overlap);

extern double Js(AdvectDiffuse2D_Quad_Block &SolnBlk, 
		 int ncol, 
		 int nrow, 
		 int JsS, 
		 int JdS, 
		 int JaS, 
		 int overlap);

extern int Newton_Krylov_Schwarz_Solver(ostream &Progress_File,
                                        int &Number_of_Startup_Interations,
					AdvectDiffuse2D_Quad_Block *Soln_ptr, 
				        AdaptiveBlock2D_List &Soln_Block_List, 
  				        AdvectDiffuse2D_Input_Parameters &Input_Parameters);

#endif  /* _ADVECTDIFFUSE2D_QUAD_NKS_INCLUDED */



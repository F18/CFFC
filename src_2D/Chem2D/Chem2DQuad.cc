/**********************************************************************
 * Chem2DQuad.cc: Subroutines for the 2D chem quadrilateral solution  *
 *                block class.                                        *
 **********************************************************************/

// Include 2D Chem quadrilateral solution block class.

#ifndef _CHEM2D_QUAD_INCLUDED
#include "Chem2DQuad.h"
#endif // _CHEM2D_QUAD_INCLUDED

/**********************************************************************
 * Chem2D_Quad_Block -- Create storage for the total number of        *
 *                      variables.                                    *
 **********************************************************************/
int Chem2D_Quad_Block::residual_variable = 1;
int Chem2D_Quad_Block::Number_of_Residual_Norms = 4;


/********************************************************
 * Routine: FixSpecGrad                                 *
 *                                                      *
 * Recompute the gradient of the last species to ensure *
 * the gradients sum to zero.                           *
 *                                                      *
 ********************************************************/
void Chem2D_Quad_Block::FixSpecGrad(const int i,const int j, 
				    const bool&visc, const int nsp) {

  double sumX(0.0), sumY(0.0), 
         sumX_faceN(0.0), sumY_faceN(0.0),
         sumX_faceS(0.0), sumY_faceS(0.0),
         sumX_faceE(0.0), sumY_faceE(0.0),
         sumX_faceW(0.0), sumY_faceW(0.0);

  // For inviscid case
  if (!visc) {
    for (int k=0; k<nsp-1; k++) {
      sumX += dWdx[i][j].spec[k].c;
      sumY += dWdy[i][j].spec[k].c;
    }
    dWdx[i][j].spec[nsp-1].c = -sumX;
    dWdy[i][j].spec[nsp-1].c = -sumY;

  // For viscous case
  } else {
    for (int k=0; k<nsp-1; k++) {
      sumX += dWdx[i][j].spec[k].c;
      sumY += dWdy[i][j].spec[k].c;
      sumX_faceN += dWdx_faceN[i][j].spec[k].c;
      sumX_faceS += dWdx_faceS[i][j].spec[k].c;
      sumX_faceE += dWdx_faceE[i][j].spec[k].c;
      sumX_faceW += dWdx_faceW[i][j].spec[k].c;
      sumY_faceN += dWdy_faceN[i][j].spec[k].c;
      sumY_faceS += dWdy_faceS[i][j].spec[k].c;
      sumY_faceE += dWdy_faceE[i][j].spec[k].c;
      sumY_faceW += dWdy_faceW[i][j].spec[k].c;
    }
    dWdx[i][j].spec[nsp-1].c = -sumX;
    dWdy[i][j].spec[nsp-1].c = -sumY;
    dWdx_faceN[i][j].spec[nsp-1].c = -sumX_faceN;
    dWdx_faceS[i][j].spec[nsp-1].c = -sumX_faceS;
    dWdx_faceE[i][j].spec[nsp-1].c = -sumX_faceE;
    dWdx_faceW[i][j].spec[nsp-1].c = -sumX_faceW;
    dWdy_faceN[i][j].spec[nsp-1].c = -sumY_faceN;
    dWdy_faceS[i][j].spec[nsp-1].c = -sumY_faceS;
    dWdy_faceE[i][j].spec[nsp-1].c = -sumY_faceE;
    dWdy_faceW[i][j].spec[nsp-1].c = -sumY_faceW;
  }
}

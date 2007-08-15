#ifndef _HIGHTEMP2D_DRDU_INCLUDED
#define _HIGHTEMP2D_DRDU_INCLUDED

#include "HighTemp2DQuad.h"
#include "../Math/Matrix.h"

DenseMatrix Rotation_Matrix_NS2D(const Vector2D &nface, int Size, int A_matrix);

void dFIdW_Inviscid_ROE(DenseMatrix& dRdW, const HighTemp2D_Quad_Block &SolnBlk,  
			const HighTemp2D_Input_Parameters &Input_Parameters,
			int ii, int jj, int Orient);

int Inviscid_Flux_Used_Reconstructed_LeftandRight_States(
		HighTemp2D_pState &Wl, HighTemp2D_pState &Wr, 
		const HighTemp2D_Quad_Block &SolnBlk, 
		int Orient, int i, int j);

void dFvdWf_Diamond(DenseMatrix &dFvdWf, DenseMatrix &dGvdWf, 
		const HighTemp2D_Quad_Block &SolnBlk,
		int Orient_face, int Rii, int Rjj);

void dWfdWc_Diamond(DenseMatrix &dWfdWc_x, DenseMatrix &dWfdWc_y, 
		const HighTemp2D_Quad_Block &SolnBlk,
		int Orient_face, int Rii, int Rjj, int Orient_cell);

#endif // _HIGHTEMP2D_DRDU_INCLUDED


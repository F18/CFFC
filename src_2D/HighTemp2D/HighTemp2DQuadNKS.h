#ifndef _HIGHTEMP2D_QUAD_NKS_INCLUDED 
#define _HIGHTEMP2D_QUAD_NKS_INCLUDED 

/****************************************************************/
/****** NKS REQUIRED SPECIALIZATIONS & FUNCTIONS ****************/
/****************************************************************/  

#include "HighTemp2DQuad.h"
#include "../NewtonKrylovSchwarz2D/NKS2D.h"
#include "HighTemp2DdRdU.h"

/********************************************************************************************
 *  Specialization of Newton_Update Function                                                *
 *                                                                                          *
 * This routine updates the previous solution data Uo with the deltaU from the GMRES        *
 * iterative solver.  U = Uo + GMRES.delatU                                                 *
 *                                                                                          *
 ********************************************************************************************/
template <>
int Newton_Update(HighTemp2D_Quad_Block *SolnBlk,
		  AdaptiveBlock2D_List &List_of_Local_Solution_Blocks,
		  HighTemp2D_Input_Parameters &Input_Parameters,
		  GMRES_RightPrecon_MatrixFree<HighTemp2D_pState,
				               HighTemp2D_Quad_Block,
				               HighTemp2D_Input_Parameters> &GMRES,
		  double Relaxation_multiplier) 
{

	int Num_Var = SolnBlk[0].NumVar();  

	/* Update Solution. No updates to Ghost Cells, let the BC's take care of it */
	for ( int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks.Nblk ; ++Bcount ) {
		if (List_of_Local_Solution_Blocks.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
			for (int j = SolnBlk[Bcount].JCl; j <= SolnBlk[Bcount].JCu; j++){
				for (int i = SolnBlk[Bcount].ICl; i <= SolnBlk[Bcount].ICu; i++){

					// Update solutions in conservative variables:
					// U = Uo + deltaU = Uo + denormalized(x)
					for(int varindex =1; varindex <= Num_Var; varindex++){	
						SolnBlk[Bcount].U[i][j][varindex] = SolnBlk[Bcount].Uo[i][j][varindex] 
							+ Relaxation_multiplier * GMRES.deltaU(Bcount,i,j,varindex-1);
					} 	      

					// Apply update reduction while any one of the updated variables is negative. 
					// This seems to cause more problems than help.
					// Any ideas? Maybe use a global step reduction (broadcast with MPI)?
					if (SolnBlk[Bcount].U[i][j].rho   <= ZERO ||
							SolnBlk[Bcount].U[i][j].e() <= ZERO ||
							SolnBlk[Bcount].U[i][j].E   <= ZERO) {

						double update_reduction_factor = ONE;

						for (int n_update_reduction = 1; n_update_reduction <= 10; ++n_update_reduction) {		  
							update_reduction_factor = HALF*update_reduction_factor;		  		  
							for(int varindex = 1; varindex <= Num_Var; varindex++){		
								SolnBlk[Bcount].U[i][j][varindex] = SolnBlk[Bcount].Uo[i][j][varindex] 
									+ GMRES.deltaU(Bcount,i,j,varindex-1)*update_reduction_factor;
							}   

							cout<<"\n Applying Reduction to solution in NKS "<<n_update_reduction;

							if (SolnBlk[Bcount].U[i][j].rho   > ZERO &&
									SolnBlk[Bcount].U[i][j].E   > ZERO &&
									SolnBlk[Bcount].U[i][j].e() > ZERO ) break; 
						} 
					} 

					if (SolnBlk[Bcount].U[i][j].rho <= ZERO) { 
						cout << "\n NEGATIVE DENSITY : " << endl;
						cout << "SolnBlk["<<Bcount<<"].U["<<i<<"]["<<j<<"].d = " 
							<< SolnBlk[Bcount].U[i][j].rho 
							<< ": SolnBlk["<<Bcount<<"].Uo["<<i<<"]["<<j<<"].d = " 
							<< SolnBlk[Bcount].Uo[i][j].rho << endl;
						cout << "   G["<<Bcount<<"].x[G["<<Bcount<<"].index("<<i
							<<","<<j
							<<")] = " 
							<< GMRES.deltaU(Bcount,i,j,0) << endl;
					}	else if (SolnBlk[Bcount].U[i][j].e() <= ZERO) { 
						cout << "\n NEGATIVE INTERNAL ENERGY : " << endl;
						cout << "SolnBlk["<<Bcount<<"].U["<<i<<"]["<<j<<"].e() = " 
							<< SolnBlk[Bcount].U[i][j].e() << endl;
					}

					// Update solution in primitive variables.
					SolnBlk[Bcount].W[i][j] = W(SolnBlk[Bcount].U[i][j]);	  
				} 
			} 
		} 
	}   

	return 0; 
}

// HighTemp2D Specific Finite_Time_Step
template <>
double Finite_Time_Step(const HighTemp2D_Input_Parameters &Input_Parameters, 
			const double &L2norm_first,
			const double &L2norm_current,
			const double &L2norm_current_n,		
			const int &Number_of_Newton_Steps) {

	double i = Input_Parameters.NKS_IP.Finite_Time_Step_Initial_CFL;
	double f = Input_Parameters.NKS_IP.Finite_Time_Step_Final_CFL;
	double L = max(Input_Parameters.NKS_IP.Overall_Tolerance, NKS_EPS);

	// Find x such that:
	// (i) (L)^(x) == f
	double CFL_Power = log(f/i) / log(L);

	// If L2norm_current_n were a straight line on a semi-log plot 
	//  then CFL would  also be a straight line on a semi-log plot.
	return Input_Parameters.NKS_IP.Finite_Time_Step_Initial_CFL *
	       pow(L2norm_current_n, CFL_Power);
}

/*!**************************************************************
 * Specialization of Block_Preconditioner::Preconditioner_dFdU  *
 *                                                              *
 * Calculates the dFdU matrix used to generate the approximate  *               
 * Jacobian for the Block Preconditioner.                       *
 ****************************************************************/
template<> 
inline void Block_Preconditioner<HighTemp2D_pState,
			         HighTemp2D_Quad_Block,					    
			         HighTemp2D_Input_Parameters>::
Preconditioner_dFIdU(DenseMatrix &_dFdU, HighTemp2D_pState W)
{
  W.dFdU(_dFdU);
}

/*!**************************************************************
 * Specialization of Block_Preconditioner::Preconditioner_dFIdU *
 *                                                              *
 * Calculates the dFdU matrix used to generate the approximate  *               
 * Jacobian for the Block Preconditioner.                       *
 ****************************************************************/ 
template<> inline void Block_Preconditioner<HighTemp2D_pState,
					    HighTemp2D_Quad_Block,			    
					    HighTemp2D_Input_Parameters>::
Preconditioner_dFIdU_Roe(DenseMatrix &_dFIdU, int ii, int jj, int Orient)
{ 
  DenseMatrix dFI_dW(NUM_VAR_HIGHTEMP2D,NUM_VAR_HIGHTEMP2D,ZERO);
  
  dFIdW_Inviscid_ROE(dFI_dW, *SolnBlk, *Input_Parameters, ii, jj, Orient);  
 
  DenseMatrix dWdU(NUM_VAR_HIGHTEMP2D,NUM_VAR_HIGHTEMP2D,ZERO); 
  //transformation Jacobian  Wo == W here 
  SolnBlk->W[ii][jj].dWdU(dWdU);

  _dFIdU += dFI_dW*dWdU;
}

/*!**************************************************************
 * Specialization of Block_Preconditioner::                     *
 *                               normalize_Preconditioner_dFdU  *
 *                                                              *
 * Normalizes the dFdU matrix used to generate the approximate *               
 * Jacobian for the Block Preconditioner.                       *
 ****************************************************************/
template<> 
inline void Block_Preconditioner<HighTemp2D_pState,
			         HighTemp2D_Quad_Block,
			         HighTemp2D_Input_Parameters>::
normalize_Preconditioner_dFdU(DenseMatrix &dFdU) 
{
	double ao  = HighTemp2D_W_STDATM.a();
	
  // First row of dFdU is partial of (rho)(u) w.r.t. [U] 
  // and as such should be scaled as:
  dFdU(0,0) *= ONE/ao;
//dFdU(0,1) *= ONE;
//dFdU(0,2) *= ONE;
  dFdU(0,3) *= ao;
  for (int i = 4; i < blocksize; i++) {
    dFdU(0,i) *= ONE/ao;
  }

  // Next row is partial of (rho)(u)(u) w.r.t. [U]
  // so it is scaled as the first row, but with an additional
  // division by ao for each entry:
  dFdU(1,0) *= ONE/ao/ao;
  dFdU(1,1) *= ONE/ao;
  dFdU(1,2) *= ONE/ao;
//dFdU(1,3) *= ONE;
  for (int i = 4; i < blocksize; i++) {
    dFdU(1,i) *= ONE/ao/ao;
  }

  // Next row is the same.
  dFdU(2,0) *= ONE/ao/ao;
  dFdU(2,1) *= ONE/ao;
  dFdU(2,2) *= ONE/ao;
//dFdU(2,3) *= ONE;
  for (int i = 4; i < blocksize; i++) {
    dFdU(2,i) *= ONE/ao/ao;
  }

  // Add an additional /ao for energy flux:
  dFdU(3,0) *= ONE/ao/ao/ao;
  dFdU(3,1) *= ONE/ao/ao;
  dFdU(3,2) *= ONE/ao/ao;
  dFdU(3,3) *= ONE/ao;
  for (int i = 4; i < blocksize; i++) {
    dFdU(3,i) *= ONE/ao/ao/ao;
  }

  // and then for the turbulence modelling: 
  for (int i = 4; i < blocksize; i++) {
    dFdU(i,0) *= ONE/ao;
  //dFdU(i,1) *= ONE;
  //dFdU(i,2) *= ONE;
    dFdU(i,3) *= ao;
    for (int j = 4; j < blocksize; j++) {
      dFdU(i,j) *= ONE/ao;
      dFdU(i,j) *= ONE/ao;
    }
  }
}

/************************************************************************/
/************ GMRES REQUIRED SPECIALIZATIONS & FUNCTIONS ****************/
/************************************************************************/   

/*!********************************************************
 * GMRES_Block::set_normalize_values for NKS/GMRES        * 
 *              *                                         *
 *   normalize_values[0] must be set to ao                *
 *   normalize_values[1-n] = values for index[1-n]        *
 *             where n = the number of solution variables *
 **********************************************************/
template<> 
inline void GMRES_Block<HighTemp2D_pState,
			HighTemp2D_Quad_Block,
			HighTemp2D_Input_Parameters>::
set_normalize_values(void)
{   
	double ao  = HighTemp2D_W_STDATM.a();
	double rho = HighTemp2D_W_STDATM.rho;

	normalize_valuesU[0] = rho;          //rho
	normalize_valuesU[1] = rho*ao;       //rho*u
	normalize_valuesU[2] = rho*ao;       //rho*v
	normalize_valuesU[3] = rho*ao*ao;    //rho*e
  // turbulence modelling variables all scale by the same value.
  for (int i = 4; i < blocksize; i++) {
    normalize_valuesU[i] = rho;
  }

	normalize_valuesR[0] = rho*ao;
	normalize_valuesR[1] = rho*ao*ao;
	normalize_valuesR[2] = rho*ao*ao;
	normalize_valuesR[3] = rho*ao*ao*ao;
  for (int i = 4; i < blocksize; i++) {
    normalize_valuesR[i] = rho*ao;
  }
}

/*!**************************************************************
 * Specialization of Block_Preconditioner::Preconditioner_dFVdU *
 *                                                              *
 * Calculates the dRdU matrix used to generate the approximate  *               
 * Jacobian for the Block Preconditioner.                       *
 * Consult the comments above the function declaration in 
 * NewtonKrylovSchwarz2D/Block_Preconditioner2D.h
 * for an explanation of the input variables to this function
 ****************************************************************/
template<> 
inline void Block_Preconditioner<HighTemp2D_pState,
			         HighTemp2D_Quad_Block,					    
			         HighTemp2D_Input_Parameters>::
Preconditioner_dFVdU(DenseMatrix &dFvdU, int Rii, int Rjj, 
		     int Wii, int Wjj, int Orient_face, int Orient_cell) {

  double lface;
  Vector2D nface;

  switch(Orient_face){
  case NORTH:
    nface = SolnBlk->Grid.nfaceN(Rii, Rjj);
    lface = SolnBlk->Grid.lfaceN(Rii, Rjj);
    break;
  case SOUTH:
    nface = SolnBlk->Grid.nfaceS(Rii, Rjj);
    lface = SolnBlk->Grid.lfaceS(Rii, Rjj);
    break;
  case EAST:
    nface = SolnBlk->Grid.nfaceE(Rii, Rjj);
    lface = SolnBlk->Grid.lfaceE(Rii, Rjj);
    break;  
  case WEST:
    nface = SolnBlk->Grid.nfaceW(Rii, Rjj);
    lface = SolnBlk->Grid.lfaceW(Rii, Rjj);
    break;  
  }

  int spatial_dimension = 2;
  int n_vars_wo_vel = NUM_VAR_HIGHTEMP2D - spatial_dimension;

  // Wfx := [rho, rhox, u, ux, uy, v, vx, vy, p, px, k, kx, etc ...]
  int n_face_vars = (n_vars_wo_vel * 2) + spatial_dimension*(spatial_dimension+1);

  DenseMatrix dFvdWf(NUM_VAR_HIGHTEMP2D, n_face_vars, ZERO);
  DenseMatrix dGvdWf(NUM_VAR_HIGHTEMP2D, n_face_vars, ZERO);

  DenseMatrix dWfdWc_x(n_face_vars, NUM_VAR_HIGHTEMP2D, ZERO); 
  DenseMatrix dWfdWc_y(n_face_vars, NUM_VAR_HIGHTEMP2D, ZERO); 

  dFvdWf_Diamond(dFvdWf,   dGvdWf,   *SolnBlk, Orient_face, Rii, Rjj);
  dWfdWc_Diamond(dWfdWc_x, dWfdWc_y, *SolnBlk, Orient_face, Rii, Rjj, Orient_cell); 
  
  DenseMatrix dFvdWc = lface * (nface.x*(dFvdWf*dWfdWc_x) + nface.y*(dGvdWf*dWfdWc_y));
    
  DenseMatrix dWdU(NUM_VAR_HIGHTEMP2D, NUM_VAR_HIGHTEMP2D, ZERO); 

  SolnBlk->W[Wii][Wjj].dWdU(dWdU);

  dFvdU += dFvdWc*dWdU;
}

template <>
inline void Block_Preconditioner<HighTemp2D_pState,
			         HighTemp2D_Quad_Block,					    
			         HighTemp2D_Input_Parameters>::
First_Order_Inviscid_Jacobian_HLLE(const int &ci, const int &cj, DenseMatrix* J_column_data)
{

	Grid2D_Quad_Block *G = &(SolnBlk->Grid);
	double this_cellA = G->area(ci, cj); // likely need it more than once.
	HighTemp2D_pState** W = SolnBlk->W; 
	DenseMatrix I(blocksize, blocksize); I.identity();

	DenseMatrix dFdUtmp(blocksize,blocksize,ZERO);
	int nci = 0, ncj = 0;
	double len = 0.0, wsm = 0.0, wsp = 0.0;
	Vector2D nrml;
	
	for (int Jndx = 1; Jndx <= 4; Jndx++) { 

		// Stand in the neighbouring cell. 
		// Look towards the centre cell. 
		// You are now facing in the direction given below.

		switch (Jndx) {
			case NORTH: nci = ci  ; ncj = cj-1; break;
			case WEST:  nci = ci+1; ncj = cj  ; break;
			case SOUTH: nci = ci  ; ncj = cj+1; break;
			case EAST:  nci = ci-1; ncj = cj  ; break;
		}

		switch (Jndx) {
			case NORTH: nrml = G->nfaceS(ci,cj); len = G->lfaceS(ci,cj); break;
			case WEST:  nrml = G->nfaceE(ci,cj); len = G->lfaceE(ci,cj); break;
			case SOUTH: nrml = G->nfaceN(ci,cj); len = G->lfaceN(ci,cj); break;
			case EAST:  nrml = G->nfaceW(ci,cj); len = G->lfaceW(ci,cj); break;
		}

		// HLLE wavespeeds are assumed constant wrt U for this approximate Jacobian.
		// Note that wavespeed_fcn is a function pointer.
		Vector2D lms = HLLE_wavespeeds(W[ci][cj], W[nci][ncj], nrml);
		wsm = lms.x; wsp = lms.y;

		if (wsp <= ZERO) { // The flux at this side does not depend on this cell.
			J_column_data[Jndx].zero();
			continue; 
		} 

		DenseMatrix RotMat(     Rotation_Matrix(nrml, 1) );
		DenseMatrix RotMat_inv( Rotation_Matrix(nrml, 0) );

		//  We want (an approximation for) the matrix on the LHS of:
		//  
		//  	[ - I/h + dR/dU ] dU = - R
		//  
		//  where 
		//  
		//  	R = U' = - (len/area) F 
		//  
		//  after the dot product with n=(1,0). Since the method is
		//  conservative, F calculated from the perspective of one cell
		//  is exactly the negative of F caculated at the neighbouring
		//  cell. 

		dFdUtmp.zero(); // Ideally this would not be necessary.
		Preconditioner_dFIdU(dFdUtmp, Rotate(W[ci][cj], nrml));
		J_column_data[Jndx] = RotMat_inv * dFdUtmp * RotMat; // is this efficient?

		// if wsm is nearly zero, then treat it as zero 
		// to avoid round-off error in the else clause.
		if (wsm > -1.0e-2) { 

			J_column_data[Jndx] *= len / this_cellA; 
			J_column_data[0] -= J_column_data[Jndx];
			J_column_data[Jndx] *= this_cellA / G->area(nci,ncj);

		} else if (wsp > ZERO) { // could just be an "else" statement but let's not be too slick.

			// Remember that we are determining only one column of the (approximate)
			// Jacobian. Thus the derivative of the dependence of R(U) for this cell
			// on the cell to the "right" will be calculated when the cell to the
			// "right" is processed.

			J_column_data[Jndx] *= 1.0 / wsm;
			J_column_data[Jndx] -= I; // The identity matrix is invariant under rotation.
			J_column_data[Jndx] *= len * wsp * wsm / this_cellA / (wsp - wsm);
			J_column_data[0] -= J_column_data[Jndx];
			J_column_data[Jndx] *= this_cellA / G->area(nci,ncj);
		}
	} // for (int Jndx = 1; Jndx <= 4; Jndx++) 

}

#endif // _HIGHTEMP2D_QUAD_NKS_INCLUDED 


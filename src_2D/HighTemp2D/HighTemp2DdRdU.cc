#ifdef ALI_CHECK_HIGHTEMP
#define COMPILING_DRDU // not the cleanest...
#endif

#include "HighTemp2DdRdU.h"

/*********************************************************
 * Routine: Rotation_Matrix_NS2D                           *
 *                                                       *
 * This function returns either the rotation matrix, A,  *
 * or the inverse of A.                                  *
 *                                                       *
 * Note: A_matrix = 1 for returning A.                   *
 *                = 0 for returning inverse of A.        *
 *                                                       *
 *********************************************************/
// This rotation matrix is specific to HighTemp since it makes an
// assumption about the position of the velocity vector. But doesn't almost
// every other physical class (Euler, Chem, etc) also have velocity right after
// density? But to make a "physical" super-class does seem complicated.
DenseMatrix Rotation_Matrix_NS2D(const Vector2D &nface, int Size, int A_matrix) 
{
  double cos_angle = nface.x; 
  double sin_angle = nface.y;
    
  DenseMatrix mat(Size,Size);
  mat.identity();
 
  if (A_matrix) {
    // Rotation Matrix, A 
    mat(1,1) = cos_angle;
    mat(1,2) = sin_angle;
    mat(2,1) = -sin_angle;
    mat(2,2) = cos_angle;
  } else {
    // Rotation Matrix, Inv A 
    mat(1,1) = cos_angle;
    mat(1,2) = -sin_angle;
    mat(2,1) = sin_angle;
    mat(2,2) = cos_angle;
  } /* endif */

  return mat;
} 

/********************************************************
 * Routine: Inviscid Flux Jacobian using Roe            *
 *                                                      *
 * This routine returns the inviscid components of      *    
 * Jacobian matrix for the specified local solution     *
 * block calculated analytically.                       *
 *                                                      *
 * dF/dW_R                                              *
 ********************************************************/
void dFIdW_Inviscid_ROE(DenseMatrix& dRdW, const HighTemp2D_Quad_Block &SolnBlk,  
			const HighTemp2D_Input_Parameters &Input_Parameters,
			int ii, int jj, int Orient) {
   
  int overlap = Input_Parameters.NKS_IP.GMRES_Overlap;
  int Ri, Rj;

  if (ii < SolnBlk.ICl -overlap || ii > SolnBlk.ICu + overlap ||
      jj < SolnBlk.JCl -overlap || jj > SolnBlk.JCu + overlap) {
     // GHOST CELL so do nothing
     cout<<"\n Hey I am not suppose to be here! \n"; exit(1);
  } /* endif */

  DenseMatrix dFidW(NUM_VAR_HIGHTEMP2D, NUM_VAR_HIGHTEMP2D, ZERO);
  
  Vector2D nface,DX; double lface;   
  HighTemp2D_pState Wa, wavespeeds, Left, Right, Wl, Wr;
  double c_Avg, dpdrho_Avg, dpde_Avg;
  Left.Vacuum(); Right.Vacuum(); Wl.Vacuum(); Wr.Vacuum();
  			
  if (Orient == NORTH) {
     Ri = ii; Rj=jj-1;
     nface = SolnBlk.Grid.nfaceN(Ri, Rj);
     lface = SolnBlk.Grid.lfaceN(Ri, Rj);
  } else if (Orient == SOUTH) {
     Ri = ii; Rj=jj+1;
     nface = SolnBlk.Grid.nfaceS(Ri, Rj);
     lface = SolnBlk.Grid.lfaceS(Ri, Rj);
  } else if (Orient == EAST) { 
     Ri = ii-1; Rj=jj;
     nface = SolnBlk.Grid.nfaceE(Ri, Rj);     
     lface = SolnBlk.Grid.lfaceE(Ri, Rj);
  } else if (Orient == WEST) { 
     Ri = ii+1; Rj=jj;
     nface = SolnBlk.Grid.nfaceW(Ri, Rj);
     lface = SolnBlk.Grid.lfaceW(Ri, Rj);
  } /* endif */
   
  DenseMatrix A( Rotation_Matrix_NS2D(nface,NUM_VAR_HIGHTEMP2D, 1));
  DenseMatrix AI( Rotation_Matrix_NS2D(nface,NUM_VAR_HIGHTEMP2D, 0));
  
  //Left and Right States 
  Inviscid_Flux_Used_Reconstructed_LeftandRight_States(Wl, Wr, SolnBlk, Orient, ii, jj);
  Left.Rotate(Wl, nface);
  Right.Rotate(Wr, nface);
   
  //Determin Roe Averaged State
  Wa = RoeAverage(Left, Right);
	 
#ifdef ALI_CHECK_HIGHTEMP
  if (ii == 5 && jj == 5) {
     ali_dump_diffs_global = true;
  } else {
     ali_dump_diffs_global = false;
  }
  ali_dump_diffs_cpu = CFFC_MPI::This_Processor_Number;
#endif

  // Jacobian dF/dW 
  HighTemp2D_pState Wtemp;
  Wtemp.Rotate(SolnBlk.Uo[ii][jj].W(), nface);
  Wtemp.dFdW(dFidW);
  dFidW = HALF*dFidW;
	 
  // Determine Wave Speeds
  if (Wa.eos_type == EOS_TGAS) {
      //wavespeeds = HartenFixAbs(Wa.lambda_x(),
      //  		       Left.lambda_x(),
      //		       Right.lambda_x());
     Wa.RoeAverage_SoundSpeed(c_Avg, dpdrho_Avg, dpde_Avg, Wa, Left, Right);
     wavespeeds = HartenFixAbs(Wa.lambda_x(c_Avg),
     	  		       Left.lambda_x(),
     			       Right.lambda_x());
  } else {
     wavespeeds = HartenFixAbs(Wa.lambda_x(),
	  		       Left.lambda_x(),
			       Right.lambda_x());
  } /* endif */
		 
  //Loop through each wavespeed and each element of Jacobian(i,j)        
  for (int i=1; i <= NUM_VAR_HIGHTEMP2D; i++) {		   
     for (int irow =0; irow< NUM_VAR_HIGHTEMP2D; irow++){
	for (int jcol =0; jcol< NUM_VAR_HIGHTEMP2D; jcol++){
            if (Wa.eos_type == EOS_TGAS) { 
	       //dFidW(irow, jcol) -= HALF*wavespeeds[i]*Wa.lp_x(i)[jcol+1]*
	       //                     Wa.rc_x(i)[irow+1];
               dFidW(irow, jcol) -= HALF*wavespeeds[i]*Wa.lp_x(i, c_Avg)[jcol+1]*
	                            Wa.rc_x(i, dpde_Avg, dpdrho_Avg, c_Avg)[irow+1];
            } else {
	       dFidW(irow, jcol) -= HALF*wavespeeds[i]*Wa.lp_x(i)[jcol+1]*
                                    Wa.rc_x(i)[irow+1];
            } /* endif */
	} /* endfor */
     } /* endfor */
  } /* endfor */

  //Rotate back 
  dRdW += lface*AI*dFidW*A;
        
}

//Needed by inviscid flux Jacobian -- reconstructed higher order left and right solution states
int Inviscid_Flux_Used_Reconstructed_LeftandRight_States(
		HighTemp2D_pState &Wl, HighTemp2D_pState &Wr, 
		const HighTemp2D_Quad_Block &SolnBlk, 
		int Orient, int i, int j)
{
	switch(Orient) {
		case WEST:
			//Wl = SolnBlk.W[i][j];  
			Wl = SolnBlk.Uo[i][j].W(); 
			//	   + ((SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x + (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y); 
			//DX = SolnBlk.Grid.xfaceW(i+1, j)-SolnBlk.Grid.Cell[i+1][j].Xc;      

			/* Reconstruct left and right solution states at the east-west faces */
			if (i == SolnBlk.ICl && 
					(SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION ||
					 SolnBlk.Grid.BCtypeW[j] == BC_BURNING_SURFACE ||
					 SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX ||
					 SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
					 SolnBlk.Grid.BCtypeW[j] == BC_MOVING_WALL_HEATFLUX ||
					 SolnBlk.Grid.BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL ||
					 SolnBlk.Grid.BCtypeW[j] == BC_RINGLEB_FLOW ||
					 SolnBlk.Grid.BCtypeW[j] == BC_CHARACTERISTIC)) {

				if (SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION) {
					Wr = Reflect(Wl,SolnBlk.Grid.nfaceW(i,j));

				} else if (SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX) {
					Wr = WallViscousHeatFlux(Wl,SolnBlk.Grid.nfaceW(i,j));

				} else if (SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL) {
					Wr = WallViscousIsothermal(Wl,SolnBlk.Grid.nfaceW(i,j),SolnBlk.Twall);

				} else if (SolnBlk.Grid.BCtypeW[j] == BC_MOVING_WALL) {
					Wr = MovingWallHeatFlux(Wl,SolnBlk.Grid.nfaceW(i,j),SolnBlk.Vwall.x);

				} else if (SolnBlk.Grid.BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL) {
					Wr = MovingWallIsothermal(Wl,SolnBlk.Grid.nfaceW(i,j),SolnBlk.Vwall.x,SolnBlk.Twall);

				} else if (SolnBlk.Grid.BCtypeW[j] == BC_RINGLEB_FLOW) {
					Wr = RinglebFlow(Wl,SolnBlk.Grid.xfaceW(i,j));

				} else {
					Wr = BC_Characteristic_Pressure(Wl,SolnBlk.WoW[j],SolnBlk.Grid.nfaceW(i,j));
				}

			} else {
				//Wr = SolnBlk.W[i-1][j];
				Wr = SolnBlk.Uo[i-1][j].W();
			}     
			break;

		case EAST:
			//Wl = SolnBlk.W[i][j]; 
			Wl = SolnBlk.Uo[i][j].W(); 
			//     + ((SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x + (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y);
			//DX = SolnBlk.Grid.xfaceE(i, j)-SolnBlk.Grid.Cell[i][j].Xc;

			if (i == SolnBlk.ICu && 
					(SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION ||
					 SolnBlk.Grid.BCtypeE[j] == BC_BURNING_SURFACE ||
					 SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX ||
					 SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
					 SolnBlk.Grid.BCtypeE[j] == BC_MOVING_WALL_HEATFLUX ||
					 SolnBlk.Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL ||
					 SolnBlk.Grid.BCtypeE[j] == BC_RINGLEB_FLOW ||
					 SolnBlk.Grid.BCtypeE[j] == BC_CHARACTERISTIC)) {

				if (SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION) {
					Wr = Reflect(Wl,SolnBlk.Grid.nfaceE(i,j));

				} else if (SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX) {
					Wr = WallViscousHeatFlux(Wl,SolnBlk.Grid.nfaceE(i,j));

				} else if (SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL) {
					Wr = WallViscousIsothermal(Wl,SolnBlk.Grid.nfaceE(i,j),SolnBlk.Twall);

				} else if (SolnBlk.Grid.BCtypeE[j] == BC_MOVING_WALL_HEATFLUX) {
					Wr = MovingWallHeatFlux(Wl,SolnBlk.Grid.nfaceE(i,j),SolnBlk.Vwall.x);

				} else if (SolnBlk.Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL) {
					Wr = MovingWallIsothermal(Wl,SolnBlk.Grid.nfaceE(i,j),SolnBlk.Vwall.x,SolnBlk.Twall);

				} else if (SolnBlk.Grid.BCtypeE[j] == BC_RINGLEB_FLOW) {
					Wr = RinglebFlow(Wl,SolnBlk.Grid.xfaceE(i,j));

				} else {
					Wr = BC_Characteristic_Pressure(Wl,SolnBlk.WoE[j],SolnBlk.Grid.nfaceE(i,j));
				}

			} else { 
				//Wr = SolnBlk.W[i+1][j];
				Wr = SolnBlk.Uo[i+1][j].W();
			}       
			break;

		case SOUTH:

			// Wl = SolnBlk.W[i][j]; 
			Wl = SolnBlk.Uo[i][j].W(); 
			//+ (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x + (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y);
			// DX = SolnBlk.Grid.xfaceS(i, j)-SolnBlk.Grid.Cell[i][j].Xc;

			/* Reconstruct left and right solution states at the north-south faces */
			if (j == SolnBlk.JCl && 
					(SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION ||
					 SolnBlk.Grid.BCtypeS[i] == BC_BURNING_SURFACE ||
					 SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX ||
					 SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
					 SolnBlk.Grid.BCtypeS[i] == BC_MOVING_WALL_HEATFLUX ||
					 SolnBlk.Grid.BCtypeS[i] == BC_MOVING_WALL_ISOTHERMAL ||
					 SolnBlk.Grid.BCtypeS[i] == BC_RINGLEB_FLOW ||
					 SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC)) {

				if (SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION) {
					Wr = Reflect(Wl,SolnBlk.Grid.nfaceS(i,j));

				} else if (SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX) {
					Wr = WallViscousHeatFlux(Wl,SolnBlk.Grid.nfaceS(i,j));

				} else if (SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL) {
					Wr = WallViscousIsothermal(Wl,SolnBlk.Grid.nfaceS(i,j),SolnBlk.Twall);

				} else if (SolnBlk.Grid.BCtypeS[i] == BC_MOVING_WALL_HEATFLUX) {
					Wr = MovingWallHeatFlux(Wl,SolnBlk.Grid.nfaceS(i,j),SolnBlk.Vwall.x);

				} else if (SolnBlk.Grid.BCtypeS[i] == BC_MOVING_WALL_ISOTHERMAL) {
					Wr = MovingWallIsothermal(Wl,SolnBlk.Grid.nfaceS(i,j),SolnBlk.Vwall.x,SolnBlk.Twall);

				} else if (SolnBlk.Grid.BCtypeS[i] == BC_RINGLEB_FLOW) {
					Wr = RinglebFlow(Wl,SolnBlk.Grid.xfaceS(i,j));

				} else if (SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC) {
					Wr = BC_Characteristic_Pressure(Wl,SolnBlk.WoS[i],SolnBlk.Grid.nfaceS(i,j));
				}

			} else {
				//Wr = SolnBlk.W[i][j-1];
				Wr = SolnBlk.Uo[i][j-1].W();  
			}
			break;

		case NORTH:

			//Wl = SolnBlk.W[i][j]; 
			Wl = SolnBlk.Uo[i][j].W(); 
			//+ (SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x + (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y);
			// DX = SolnBlk.Grid.xfaceN(i, j)-SolnBlk.Grid.Cell[i][j].Xc;

			if (j == SolnBlk.JCu && 	
					(SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION ||
					 SolnBlk.Grid.BCtypeN[i] == BC_BURNING_SURFACE ||
					 SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX ||
					 SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
					 SolnBlk.Grid.BCtypeN[i] == BC_MOVING_WALL_HEATFLUX ||
					 SolnBlk.Grid.BCtypeN[i] == BC_MOVING_WALL_ISOTHERMAL ||
					 SolnBlk.Grid.BCtypeN[i] == BC_RINGLEB_FLOW ||
					 SolnBlk.Grid.BCtypeN[i] == BC_CHARACTERISTIC)) {

				if (SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION) {
					Wr = Reflect(Wl,SolnBlk.Grid.nfaceN(i,j));

				} else if (SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX) {
					Wr = WallViscousHeatFlux(Wl,SolnBlk.Grid.nfaceN(i,j));

				} else if (SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL) {
					Wr = WallViscousIsothermal(Wl,SolnBlk.Grid.nfaceN(i,j),SolnBlk.Twall);

				} else if (SolnBlk.Grid.BCtypeN[i] == BC_MOVING_WALL_HEATFLUX) {
					Wr = MovingWallHeatFlux(Wl,SolnBlk.Grid.nfaceN(i,j),SolnBlk.Vwall.x);

				} else if (SolnBlk.Grid.BCtypeN[i] == BC_MOVING_WALL_ISOTHERMAL) {
					Wr = MovingWallIsothermal(Wl,SolnBlk.Grid.nfaceN(i,j),SolnBlk.Vwall.x,SolnBlk.Twall);

				} else if (SolnBlk.Grid.BCtypeN[i] == BC_RINGLEB_FLOW) {
					Wr = RinglebFlow(Wl,SolnBlk.Grid.xfaceN(i,j));

				} else {
					Wr = BC_Characteristic_Pressure(Wl,SolnBlk.WoN[i],SolnBlk.Grid.nfaceN(i,j));
				}
				
			} else {      
				//Wr = SolnBlk.W[i][j+1]; 
				Wr = SolnBlk.Uo[i][j+1].W(); 
			}    
			break;
	}

	return 0;
}

#ifdef ALI_CHECK_HIGHTEMP
bool dump_diffs = true;
#endif

//  dFvdWf_Diamond:
//
//  This function calculates the Jacobian of the viscous
//  components of the residual with respect to the primitive
//  variables. 
//
//  Specifically, this returns in dFvdWf and dGvdWf the Jacobian
//  of the viscous part of the equation written for the
//  Orient_face face of cell (Rii, Rjj) with respect to the
//  primitive variables at that face. Since the equations are
//  non-linear this will depend on the values at the face. That
//  is, the structure of the matrices is independent of the
//  method used to determine the face gradients but the values do
//  depend on this method.

#define USE_INTERPOLATION_FOR_WFACE_IN_PRECONDITIONER 1

void dFvdWf_Diamond(DenseMatrix &dFvdWf, DenseMatrix &dGvdWf, 
		const HighTemp2D_Quad_Block &SolnBlk,
		int Orient_face, int Rii, int Rjj){

	// face_grad_work. Put these asserts elsewhere?
	assert(SolnBlk.face_grad_arrays_allocated);
	assert(SolnBlk.dWdx_faceN != NULL);

	HighTemp2D_pState QuadraturePoint_W;
	double rhox = 0, px = 0, ux = 0, vx = 0, kx = 0, omegax = 0;
	double rhoy = 0, py = 0, uy = 0, vy = 0, ky = 0, omegay = 0;

#if USE_INTERPOLATION_FOR_WFACE_IN_PRECONDITIONER
	int error_flag = 0;
#endif

	switch(Orient_face){
		case NORTH:
#if USE_INTERPOLATION_FOR_WFACE_IN_PRECONDITIONER
			// find QuadraturePoint_W:
			error_flag = Bilinear_Interpolation(
					SolnBlk.WwNE(Rii, Rjj), SolnBlk.Grid.nodeNE(Rii, Rjj).X,
					SolnBlk.W[Rii][Rjj],    SolnBlk.Grid.Cell[Rii][Rjj].Xc,
					SolnBlk.WwNW(Rii, Rjj), SolnBlk.Grid.nodeNW(Rii, Rjj).X,
					SolnBlk.W[Rii][Rjj+1],  SolnBlk.Grid.Cell[Rii][Rjj+1].Xc,
			    SolnBlk.Grid.xfaceN(Rii, Rjj), QuadraturePoint_W);
			// end if error_flag? Well that's a good question.
#else
			QuadraturePoint_W = HALF*(SolnBlk.WwNW(Rii, Rjj) + SolnBlk.WwNE(Rii, Rjj));
#endif

			rhox = SolnBlk.dWdx_faceN[Rii][Rjj].rho;
			rhoy = SolnBlk.dWdy_faceN[Rii][Rjj].rho;
			px = SolnBlk.dWdx_faceN[Rii][Rjj].p;
			py = SolnBlk.dWdy_faceN[Rii][Rjj].p;
			ux = SolnBlk.dWdx_faceN[Rii][Rjj].v.x;
			uy = SolnBlk.dWdy_faceN[Rii][Rjj].v.x;
			vx = SolnBlk.dWdx_faceN[Rii][Rjj].v.y;
			vy = SolnBlk.dWdy_faceN[Rii][Rjj].v.y;
			kx = SolnBlk.dWdx_faceN[Rii][Rjj].k;
			ky = SolnBlk.dWdy_faceN[Rii][Rjj].k;
			omegax = SolnBlk.dWdx_faceN[Rii][Rjj].omega;
			omegay = SolnBlk.dWdy_faceN[Rii][Rjj].omega;

			break;
		case EAST:
#if USE_INTERPOLATION_FOR_WFACE_IN_PRECONDITIONER
			// find QuadraturePoint_W:
			error_flag = Bilinear_Interpolation(
					SolnBlk.W[Rii+1][Rjj],  SolnBlk.Grid.Cell[Rii+1][Rjj].Xc,
					SolnBlk.WwSE(Rii, Rjj), SolnBlk.Grid.nodeSE(Rii, Rjj).X,
					SolnBlk.W[Rii][Rjj],    SolnBlk.Grid.Cell[Rii][Rjj].Xc,
					SolnBlk.WwNE(Rii, Rjj), SolnBlk.Grid.nodeNE(Rii, Rjj).X,
					SolnBlk.Grid.xfaceE(Rii, Rjj),
					QuadraturePoint_W);
#else
			QuadraturePoint_W = HALF*(SolnBlk.WwNE(Rii, Rjj) + SolnBlk.WwSE(Rii, Rjj) );     
#endif

			rhox = SolnBlk.dWdx_faceE[Rii][Rjj].rho; 
			rhoy = SolnBlk.dWdy_faceE[Rii][Rjj].rho;
			px = SolnBlk.dWdx_faceE[Rii][Rjj].p;
			py = SolnBlk.dWdy_faceE[Rii][Rjj].p;
			ux = SolnBlk.dWdx_faceE[Rii][Rjj].v.x;
			uy = SolnBlk.dWdy_faceE[Rii][Rjj].v.x;
			vx = SolnBlk.dWdx_faceE[Rii][Rjj].v.y;
			vy = SolnBlk.dWdy_faceE[Rii][Rjj].v.y;
			kx = SolnBlk.dWdx_faceE[Rii][Rjj].k;
			ky = SolnBlk.dWdy_faceE[Rii][Rjj].k;
			omegax = SolnBlk.dWdx_faceE[Rii][Rjj].omega;
			omegay = SolnBlk.dWdy_faceE[Rii][Rjj].omega;

			break;
		case SOUTH:      
#if USE_INTERPOLATION_FOR_WFACE_IN_PRECONDITIONER
			// find QuadraturePoint_W:
			error_flag = Bilinear_Interpolation(
					SolnBlk.WwSE(Rii, Rjj), SolnBlk.Grid.nodeSE(Rii, Rjj).X,
					SolnBlk.W[Rii][Rjj-1],  SolnBlk.Grid.Cell[Rii][Rjj-1].Xc,
					SolnBlk.WwSW(Rii, Rjj), SolnBlk.Grid.nodeSW(Rii, Rjj).X,
					SolnBlk.W[Rii][Rjj],    SolnBlk.Grid.Cell[Rii][Rjj].Xc,
					SolnBlk.Grid.xfaceS(Rii, Rjj),
					QuadraturePoint_W);
#else
			QuadraturePoint_W = HALF*(SolnBlk.WwSE(Rii, Rjj) + SolnBlk.WwSW(Rii, Rjj) );
#endif

			rhox = SolnBlk.dWdx_faceN[Rii][Rjj-1].rho; 
			rhoy = SolnBlk.dWdy_faceN[Rii][Rjj-1].rho;
			px = SolnBlk.dWdx_faceN[Rii][Rjj-1].p;
			py = SolnBlk.dWdy_faceN[Rii][Rjj-1].p;
			ux = SolnBlk.dWdx_faceN[Rii][Rjj-1].v.x;
			uy = SolnBlk.dWdy_faceN[Rii][Rjj-1].v.x;
			vx = SolnBlk.dWdx_faceN[Rii][Rjj-1].v.y;
			vy = SolnBlk.dWdy_faceN[Rii][Rjj-1].v.y;
			kx = SolnBlk.dWdx_faceN[Rii][Rjj-1].k;
			ky = SolnBlk.dWdy_faceN[Rii][Rjj-1].k;
			omegax = SolnBlk.dWdx_faceN[Rii][Rjj-1].omega;
			omegay = SolnBlk.dWdy_faceN[Rii][Rjj-1].omega;  

			break;
		case WEST:       
#if USE_INTERPOLATION_FOR_WFACE_IN_PRECONDITIONER
			// find QuadraturePoint_W:
			error_flag = Bilinear_Interpolation(
					SolnBlk.W[Rii][Rjj],    SolnBlk.Grid.Cell[Rii][Rjj].Xc,
					SolnBlk.WwSW(Rii, Rjj), SolnBlk.Grid.nodeSW(Rii, Rjj).X,
					SolnBlk.W[Rii-1][Rjj],  SolnBlk.Grid.Cell[Rii-1][Rjj].Xc,
					SolnBlk.WwNW(Rii, Rjj), SolnBlk.Grid.nodeNW(Rii, Rjj).X,
					SolnBlk.Grid.xfaceW(Rii, Rjj),
					QuadraturePoint_W);
#else
			QuadraturePoint_W = HALF*(SolnBlk.WwNW(Rii, Rjj) + SolnBlk.WwSW(Rii, Rjj));
#endif

			rhox = SolnBlk.dWdx_faceE[Rii-1][Rjj].rho; 
			rhoy = SolnBlk.dWdy_faceE[Rii-1][Rjj].rho;
			px = SolnBlk.dWdx_faceE[Rii-1][Rjj].p;
			py = SolnBlk.dWdy_faceE[Rii-1][Rjj].p;
			ux = SolnBlk.dWdx_faceE[Rii-1][Rjj].v.x;
			uy = SolnBlk.dWdy_faceE[Rii-1][Rjj].v.x;
			vx = SolnBlk.dWdx_faceE[Rii-1][Rjj].v.y;
			vy = SolnBlk.dWdy_faceE[Rii-1][Rjj].v.y;
			kx = SolnBlk.dWdx_faceE[Rii-1][Rjj].k;
			ky = SolnBlk.dWdy_faceE[Rii-1][Rjj].k;
			omegax = SolnBlk.dWdx_faceE[Rii-1][Rjj].omega;
			omegay = SolnBlk.dWdy_faceE[Rii-1][Rjj].omega;

			break;
	}

	double rho = QuadraturePoint_W.rho;
	double u = QuadraturePoint_W.v.x;
	double v = QuadraturePoint_W.v.y;
	double p = QuadraturePoint_W.p;

	double dTdp=0.0, dTdrho=0.0, ddTdp=0.0, ddTdrho=0.0, ddTdpdrho=0.0;

	if (HighTemp2D_pState::eos_type == EOS_TGAS) {
		dTdp      = QuadraturePoint_W.dTdp();
		dTdrho    = QuadraturePoint_W.dTdrho();
		ddTdp     = QuadraturePoint_W.ddTdp();
		ddTdrho   = QuadraturePoint_W.ddTdrho();
		ddTdpdrho = QuadraturePoint_W.ddTdpdrho();
	}

	// The transport coefficients are not constant!  A simplifying assumption is
	// that they are constant with respect to the variables for this
	 // preconditioner Jacobian but they are still a function of temperature.
   double mu = QuadraturePoint_W.mu();
   double kappa = QuadraturePoint_W.kappa();
	double R = 0.0;
	if (HighTemp2D_pState::eos_type == EOS_IDEAL) { R = HighTemp2D_cState::R; }

   /*********************** X - DIRECTION **************************************/
/* 
                                         dFvdWf

           [   /4 ux   2 vy\
Fv_temp := [0, |---- - ----| mu, (uy + vx) mu,
           [   \ 3      3  /
           [

      /4 ux   2 vy\                       kappa (-px rho + p rhox)      ]
    u |---- - ----| mu + v (uy + vx) mu - ------------------------, 0, 0]
      \ 3      3  /                                   2                 ]
                                                   rho  R               ]

        Wf_temp := [rho, u, v, p, k, omega, rhox, ux, uy, vx, vy, px, kx, omegax]

    [0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0]

    [                            4 mu             2 mu            ]
    [0 , 0 , 0 , 0 , 0 , 0 , 0 , ---- , 0 , 0 , - ---- , 0 , 0 , 0]
    [                             3                3              ]

    [0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , mu , mu , 0 , 0 , 0 , 0]

    [kappa (-px rho + 2 p rhox)     2 (-2 ux + vy) mu                    kappa rhox
    [-------------------------- , - ----------------- , (uy + vx) mu , - ---------- , 0
    [             3                         3                                 2
    [          rho  R                                                      rho  R

            kappa p   4 u mu                   2 u mu   kappa        ]
    , 0 , - ------- , ------ , v mu , v mu , - ------ , ----- , 0 , 0]
               2        3                        3      rho R        ]
            rho  R                                                   ]

    [0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0]

    [0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0]
*/

	// We are only called once on a given dFvdWf and moreover it is initialized
	// to zero. So you might be tempted to replace the "+=" with a "=",
	// especially for entries which are equal (say (3,8) and (3,9)). However, to
	// be safe (in case, for example, someone switched the order of this and the
	// turbulent terms) keep the "+=".

	//  The momentum viscous fluxes written in terms of the primitive
	//  variables do not require a gas equation of state. As such, the
	//  first three rows of dFvdWf do not depend on the equation of state.

	dFvdWf(1,7) += 4.0/3.0*mu;
	dFvdWf(1,10) += -2.0/3.0*mu;
	dFvdWf(2,8) += mu;
	dFvdWf(2,9) += mu;

	switch (HighTemp2D_pState::eos_type) {

		case EOS_IDEAL:

			dFvdWf(3,0) += kappa*(-px*rho+2.0*p*rhox)/(rho*rho*rho)/R;
			dFvdWf(3,1) += -2.0/3.0*(-2.0*ux+vy)*mu;
			dFvdWf(3,2) += (uy+vx)*mu;
			dFvdWf(3,3) += -kappa*rhox/(rho*rho)/R;
			dFvdWf(3,6) += -kappa*p/(rho*rho)/R;
			dFvdWf(3,7) += 4.0/3.0*u*mu;
			dFvdWf(3,8) += v*mu;
			dFvdWf(3,9) += v*mu;
			dFvdWf(3,10) += -2.0/3.0*u*mu;
			dFvdWf(3,11) += kappa/rho/R;
			break;

		case EOS_TGAS:

			dFvdWf(3,0) += kappa*(ddTdpdrho*px+ddTdrho*rhox);
			dFvdWf(3,1) += 2.0/3.0*(2.0*ux-vy)*mu;
			dFvdWf(3,2) += (uy+vx)*mu;
			dFvdWf(3,3) += kappa*(ddTdp*px+ddTdpdrho*rhox);
			dFvdWf(3,6) += kappa*dTdrho;
			dFvdWf(3,7) += 4.0/3.0*u*mu;
			dFvdWf(3,8) += v*mu;
			dFvdWf(3,9) += v*mu;
			dFvdWf(3,10) += -2.0/3.0*u*mu;
			dFvdWf(3,11) += kappa*dTdp;
			break;
	}

   /*********************** Y - DIRECTION **************************************/
/*
                                         dGvdWf

           [                 /  2 ux   4 vy\
Gv_temp := [0, (uy + vx) mu, |- ---- + ----| mu,
           [                 \   3      3  /
           [

                       /  2 ux   4 vy\      kappa (-py rho + p rhoy)      ]
    u (uy + vx) mu + v |- ---- + ----| mu - ------------------------, 0, 0]
                       \   3      3  /                  2                 ]
                                                     rho  R               ]

        Wf_temp := [rho, u, v, p, k, omega, rhoy, ux, uy, vx, vy, py, ky, omegay]

    [0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0]

    [0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , mu , mu , 0 , 0 , 0 , 0]

    [                              2 mu           4 mu            ]
    [0 , 0 , 0 , 0 , 0 , 0 , 0 , - ---- , 0 , 0 , ---- , 0 , 0 , 0]
    [                               3              3              ]

    [kappa (-py rho + 2 p rhoy)                  2 (-ux + 2 vy) mu     kappa rhoy
    [-------------------------- , (uy + vx) mu , ----------------- , - ---------- , 0 ,
    [             3                                      3                  2
    [          rho  R                                                    rho  R

          kappa p     2 v mu                 4 v mu   kappa        ]
    0 , - ------- , - ------ , u mu , u mu , ------ , ----- , 0 , 0]
             2          3                      3      rho R        ]
          rho  R                                                   ]

    [0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0]

    [0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0]
*/
	dGvdWf(1,8) += mu;
	dGvdWf(1,9) += mu;
	dGvdWf(2,7) += -2.0/3.0*mu;
	dGvdWf(2,10) += 4.0/3.0*mu;

	switch (HighTemp2D_pState::eos_type) {

		case EOS_IDEAL:

			dGvdWf(3,0) += kappa*(-py*rho+2.0*p*rhoy)/(rho*rho*rho)/R;
			dGvdWf(3,1) += (uy+vx)*mu;
			dGvdWf(3,2) += 2.0/3.0*(-ux+2.0*vy)*mu;
			dGvdWf(3,3) += -kappa*rhoy/(rho*rho)/R;
			dGvdWf(3,6) += -kappa*p/(rho*rho)/R;
			dGvdWf(3,7) += -2.0/3.0*v*mu;
			dGvdWf(3,8) += u*mu;
			dGvdWf(3,9) += u*mu;
			dGvdWf(3,10) += 4.0/3.0*v*mu;
			dGvdWf(3,11) += kappa/rho/R;
			break;

		case EOS_TGAS:

			dGvdWf(3,0) += kappa*(ddTdpdrho*py+ddTdrho*rhoy);
			dGvdWf(3,1) += (uy+vx)*mu;
			dGvdWf(3,2) += -2.0/3.0*(ux-2.0*vy)*mu;
			dGvdWf(3,3) += kappa*(ddTdp*py+ddTdpdrho*rhoy);
			dGvdWf(3,6) += kappa*dTdrho;
			dGvdWf(3,7) += -2.0/3.0*v*mu;
			dGvdWf(3,8) += u*mu;
			dGvdWf(3,9) += u*mu;
			dGvdWf(3,10) += 4.0/3.0*v*mu;
			dGvdWf(3,11) += kappa*dTdp;
			break;
	 }

	// This function ends here. The rest is debugging for my sanity.
	//   -- Alistair Wood Mon Apr 09 2007 

#ifdef ALI_CHECK_HIGHTEMP 
	if (dump_diffs && (Rii == 5) && (Rjj == 5)) {
		
		dump_diffs = false;

		char fname[100];
		sprintf(fname, "viscous_dump_cpu%.3d.txt", CFFC_MPI::This_Processor_Number);
		ofstream fout;
		fout.open(fname, ios::app);

		double R = HighTemp2D_cState::R;
		
		double ival[12], hval[12];
		for (int qi = 0; qi < 12; qi++) { ival[qi] = 0.0; hval[qi] = 0.0; }

		ival[0] = kappa*(-px*rho+2.0*p*rhox)/(rho*rho*rho)/R;
		ival[1] = -2.0/3.0*(-2.0*ux+vy)*mu;
		ival[2] = (uy+vx)*mu;
		ival[3] = -kappa*rhox/(rho*rho)/R;
		ival[6] = -kappa*p/(rho*rho)/R;
		ival[7] = 4.0/3.0*u*mu;
		ival[8] = v*mu;
		ival[9] = v*mu;
		ival[10] = -2.0/3.0*u*mu;
		ival[11] = kappa/rho/R;

		hval[0] = kappa*(ddTdpdrho*px+ddTdrho*rhox);
		hval[1] = 2.0/3.0*(2.0*ux-vy)*mu;
		hval[2] = (uy+vx)*mu;
		hval[3] = kappa*(ddTdp*px+ddTdpdrho*rhox);
		hval[6] = kappa*dTdrho;
		hval[7] = 4.0/3.0*u*mu;
		hval[8] = v*mu;
		hval[9] = v*mu;
		hval[10] = -2.0/3.0*u*mu;
		hval[11] = kappa*dTdp;
		
#define QTOL 0.00000001
		
		fout << "\n x-dir\n";
		fout << "entry      ideal      hightemp    diff normalized by ideal value\n";
		for (int qi = 0; qi < 12; qi++) {
			double norm_diff = (hval[qi] - ival[qi]) / max(fabs(ival[qi]), QTOL);
			if (fabs(norm_diff) < QTOL) { continue; }
			fout << setw(4) << qi;
			fout << fixed << setprecision(4) << setw(13) << ival[qi];
			fout << fixed << setprecision(4) << setw(13) << hval[qi];
			fout << scientific << setprecision(2) << setw(12) << norm_diff;
			fout << endl;
		}
		
		for (int qi = 0; qi < 12; qi++) { ival[qi] = 0.0; hval[qi] = 0.0; }

		ival[0] = kappa*(-py*rho+2.0*p*rhoy)/(rho*rho*rho)/R;
		ival[1] = (uy+vx)*mu;
		ival[2] = 2.0/3.0*(-ux+2.0*vy)*mu;
		ival[3] = -kappa*rhoy/(rho*rho)/R;
		ival[6] = -kappa*p/(rho*rho)/R;
		ival[7] = -2.0/3.0*v*mu;
		ival[8] = u*mu;
		ival[9] = u*mu;
		ival[10] = 4.0/3.0*v*mu;
		ival[11] = kappa/rho/R;

		hval[0] = kappa*(ddTdpdrho*py+ddTdrho*rhoy);
		hval[1] = (uy+vx)*mu;
		hval[2] = -2.0/3.0*(ux-2.0*vy)*mu;
		hval[3] = kappa*(ddTdp*py+ddTdpdrho*rhoy);
		hval[6] = kappa*dTdrho;
		hval[7] = -2.0/3.0*v*mu;
		hval[8] = u*mu;
		hval[9] = u*mu;
		hval[10] = 4.0/3.0*v*mu;
		hval[11] = kappa*dTdp;
		
		fout << "\n\n\n y-dir\n";
		fout << "entry      ideal      hightemp    diff normalized by ideal value\n";
		for (int qi = 0; qi < 12; qi++) {
			double norm_diff = (hval[qi] - ival[qi]) / max(fabs(ival[qi]), QTOL);
			if (fabs(norm_diff) < QTOL) { continue; }
			fout << setw(4) << qi;
			fout << fixed << setprecision(4) << setw(13) << ival[qi];
			fout << fixed << setprecision(4) << setw(13) << hval[qi];
			fout << scientific << setprecision(2) << setw(12) << norm_diff;
			fout << endl;
		}

		fout.close();
	} // if (dump_diffs && (Rii == 5) && (Rjj == 5)) 
#endif // #ifdef ALI_CHECK_HIGHTEMP
 
}


/********************************************************
 * Routine: dWfdWc_Diamond                              *
 *                                                      *
 * This routine calculates the transformation matrix    *
 * to convert from the Cell Face to Cell Center used    *
 * for constructing the Viscous Jacobians based         *
 * on a diamond path recontstruction.                   *  
 *                                                      *
 ********************************************************/
/*

In evaluating the dFvdWf Jacobian above we used the following order of
variables in the face variable array:

Wfx:
[rho , u , v , p , k , omega , rhox , ux , uy , vx , vy , px , kx , omegax]
Wfy:
[rho , u , v , p , k , omega , rhoy , ux , uy , vx , vy , py , ky , omegay]

Now we need to determine the Jacobian of these vectors with respect to
the cell-centered variables. Now we can express the above as:

Wfx:
[intp rho , intp u , intp v , intp p , intp k , intp omega , dwdx rho ,
u dwdx , u dwdy , v dwdx , v dwdy , dwdx p , dwdx k , dwdx omega]

Wfy:
[intp rho , intp u , intp v , intp p , intp k , intp omega , dwdy rho ,
u dwdx , u dwdy , v dwdx , v dwdy , dwdy p , dwdy k , dwdy omega]

where now rho, u, v and so on now refer to the cell-centered
variables. This states that all the variables at the face are
geometric weights (by some scalar intp) of the variables at the cell
centre (and the same weight for all variables, yes?). Similarly for
the variable gradients. 

This results in a banded structure for the face to centre Jacobian.
For example, the first three rows of both Jacobians are:

[intp      0        0        0         0         0   ]
[                                                    ]
[  0      intp      0        0         0         0   ]
[                                                    ]
[  0        0      intp      0         0         0   ]
*/
																																													 
void dWfdWc_Diamond(DenseMatrix &dWfdWc_x, DenseMatrix &dWfdWc_y, 
		const HighTemp2D_Quad_Block &SolnBlk,
		int Orient_face, int Rii, int Rjj, int Orient_cell) {

	double intp = 0.0, LL = 0.0, RR = 0.0, dwdx = 0.0, dwdy = 0.0;

	node_weights(&LL, &RR,
			Orient_face, Rii, Rjj, Orient_cell, 
			SolnBlk.Grid);

	intp = 0.5 * (LL+RR);
	for (int i = 0; i < 6; i++) {
		dWfdWc_x(i, i) = intp;
		dWfdWc_y(i, i) = intp;
	} 

	// Find dwdx and dwdy, the geometric weights used to obtain the gradients 
	// of the face variables from the cell-centred variables.
	diamond_gradient_weights(&dwdx, &dwdy,
			LL, RR, 
			Orient_face, Rii, Rjj, Orient_cell, 
			SolnBlk.Grid);

	// The x-direction Jacobian should only use dwdx ...
	dWfdWc_x( 6,0) = dwdx;
	dWfdWc_x( 7,1) = dwdx;        
	dWfdWc_x( 9,2) = dwdx;
	dWfdWc_x(11,3) = dwdx;
	dWfdWc_x(12,4) = dwdx;
	dWfdWc_x(13,5) = dwdx;

	// ... except for the velocity bits.
	dWfdWc_x( 8,1) = dwdy;
	dWfdWc_x(10,2) = dwdy;

	// The y-direction Jacobian should only use dwdy ...
	dWfdWc_y( 6,0) = dwdy;   
	dWfdWc_y( 8,1) = dwdy;
	dWfdWc_y(10,2) = dwdy;
	dWfdWc_y(11,3) = dwdy;
	dWfdWc_y(12,4) = dwdy;
	dWfdWc_y(13,5) = dwdy;

	// ... except for the velocity bits.
	dWfdWc_y( 7,1) = dwdx;
	dWfdWc_y( 9,2) = dwdx;
}


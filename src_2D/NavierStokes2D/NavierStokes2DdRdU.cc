#include "NavierStokes2DdRdU.h"

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
DenseMatrix Rotation_Matrix_NS2D(const Vector2D &nface, int Size, int A_matrix) {
  double cos_angle = nface.x; 
  double sin_angle = nface.y;

  DenseMatrix mat(Size,Size);
  mat.identity();

  if (A_matrix) {
    mat(1,1) = cos_angle;
    mat(1,2) = sin_angle;
    mat(2,1) = -sin_angle;
    mat(2,2) = cos_angle;
  } else {
    mat(1,1) = cos_angle;
    mat(1,2) = -sin_angle;
    mat(2,1) = sin_angle;
    mat(2,2) = cos_angle;
  }

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
void dFIdW_Inviscid_ROE(DenseMatrix& dRdW, const NavierStokes2D_Quad_Block &SolnBlk,  
    const NavierStokes2D_Input_Parameters &Input_Parameters,
    int ii, int jj, int Orient) {

  int overlap = Input_Parameters.NKS_IP.GMRES_Overlap;
  int Ri, Rj;

  if (ii < SolnBlk.ICl -overlap || ii > SolnBlk.ICu + overlap ||
      jj < SolnBlk.JCl -overlap || jj > SolnBlk.JCu + overlap) {

    // GHOST CELL so do nothing
    cout << "\ndFIdW_Inviscid_ROE: Not suppose to be here.\n"; 
    exit(1);
  }

  DenseMatrix dFidW(NUM_VAR_NAVIERSTOKES2D, NUM_VAR_NAVIERSTOKES2D, ZERO);

  Vector2D nface,DX; double lface;   
  NavierStokes2D_pState Wa, wavespeeds, Left, Right, Wl, Wr;   
  Left.Vacuum();   Right.Vacuum();  Wl.Vacuum();  Wr.Vacuum();

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
  }

  DenseMatrix A(Rotation_Matrix_NS2D(nface,NUM_VAR_NAVIERSTOKES2D, 1));
  DenseMatrix AI(Rotation_Matrix_NS2D(nface,NUM_VAR_NAVIERSTOKES2D, 0));

  // Left and Right States 
  Inviscid_Flux_Used_Reconstructed_LeftandRight_States(
      Wl, Wr, SolnBlk, Orient, ii, jj);
  Left = Rotate(Wl, nface);
  Right = Rotate(Wr, nface);

  // Determine Roe-Averaged State
  Wa = RoeAverage(Left,Right);

  // Jacobian dF/dW         
  NavierStokes2D_pState Wtemp;
  Wtemp = Rotate(SolnBlk.Uo[ii][jj].W(), nface);
  Wtemp.dFdW(dFidW);
  dFidW = HALF*dFidW;

  // Determine Wave Speeds
  wavespeeds = HartenFixAbs( Wa.lambda_x(),
      Left.lambda_x(),
      Right.lambda_x());         

  // Loop through each wavespeed and each element of Jacobian(i,j)        
  for (int i = 1; i <= NUM_VAR_NAVIERSTOKES2D; i++) {       
    for (int irow = 0; irow < NUM_VAR_NAVIERSTOKES2D; irow++) {
      for( int jcol = 0; jcol < NUM_VAR_NAVIERSTOKES2D; jcol++) {

        dFidW(irow, jcol) -= HALF*wavespeeds[i]*Wa.lp_x(i)[jcol+1]*Wa.rc_x(i)[irow+1];  

      }
    }
  } 

  // Rotate back 
  dRdW += lface*AI*dFidW*A;

}

// Needed by inviscid flux Jacobian - reconstructed higher-order left and right solution states
int Inviscid_Flux_Used_Reconstructed_LeftandRight_States(
    NavierStokes2D_pState &Wl, NavierStokes2D_pState &Wr, 
    const NavierStokes2D_Quad_Block &SolnBlk, 
    int Orient, int i, int j) {
  switch(Orient) {
    case WEST:
      //Wl = SolnBlk.W[i][j];  
      Wl = SolnBlk.Uo[i][j].W(); 
      //     + ((SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x + (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y); 
      //DX = SolnBlk.Grid.xfaceW(i+1, j)-SolnBlk.Grid.Cell[i+1][j].Xc;      

      // Reconstruct left and right solution states at the east-west faces 
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

        } else if (SolnBlk.Grid.BCtypeW[j] == BC_BURNING_SURFACE) {
          Wr = BurningSurface(Wl,SolnBlk.Grid.nfaceW(i,j));

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

        } else if (SolnBlk.Grid.BCtypeE[j] == BC_BURNING_SURFACE) {
          Wr = BurningSurface(Wl,SolnBlk.Grid.nfaceE(i,j));

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

      // Reconstruct left and right solution states at the north-south faces 
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

        } else if (SolnBlk.Grid.BCtypeS[i] == BC_BURNING_SURFACE) {
          Wr = BurningSurface(Wl,SolnBlk.Grid.nfaceS(i,j));

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

        } else if (SolnBlk.Grid.BCtypeN[i] == BC_BURNING_SURFACE) {
          Wr = BurningSurface(Wl,SolnBlk.Grid.nfaceN(i,j));

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

//  dFvdWf_Diamond: calculates the Jacobian of the viscous
//  components of the residual with respect to the primitive
//  variables. 
//
//  Not yet working for turbulence modeling.
#define USE_INTERPOLATION_FOR_WFACE_IN_PRECONDITIONER 1
void dFvdWf_Diamond(DenseMatrix &dFvdWf, DenseMatrix &dGvdWf, 
    NavierStokes2D_Quad_Block &SolnBlk,
    int Orient_face, int Rii, int Rjj) {

  // assert(SolnBlk.face_grad_arrays_allocated);
  // assert(SolnBlk.dWdx_faceN != NULL);

  NavierStokes2D_pState QuadraturePoint_W;
  double rhox = ZERO, px = ZERO, ux = ZERO, vx = ZERO, kx = ZERO, omegax = ZERO;
  double rhoy = ZERO, py = ZERO, uy = ZERO, vy = ZERO, ky = ZERO, omegay = ZERO;

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

  double mu = QuadraturePoint_W.mu();
  double kappa = QuadraturePoint_W.kappa();
  double R = NavierStokes2D_cState::R;

  // X - DIRECTION 
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

  dFvdWf(1,7) += 4.0/3.0*mu;
  dFvdWf(1,10) += -2.0/3.0*mu;
  dFvdWf(2,8) += mu;
  dFvdWf(2,9) += mu;
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

   // Y - DIRECTION 
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
}


/********************************************************
 * Routine: dWfdWc_Diamond                              *
 *                                                      *
 * This routine calculates the transformation matrix    *
 * to convert from the cell face to cell center used    *
 * for constructing the viscous Jacobians based         *
 * on a diamond-path reconstruction.                    *  
 *                                                      *
 ********************************************************/
void dWfdWc_Diamond(DenseMatrix &dWfdWc_x, DenseMatrix &dWfdWc_y, 
    const NavierStokes2D_Quad_Block &SolnBlk,
    int Orient_face, int Rii, int Rjj, int Orient_cell) {

  double intp = ZERO, LL = ZERO, RR = ZERO, dwdx = ZERO, dwdy = ZERO;

  // Determine the influence of the cell on the face of interest.
  SolnBlk.Grid.dFacedC(LL, RR, Rii, Rjj, Orient_face, Orient_cell);

  intp = 0.5 * (LL+RR);
  for (int i = 0; i < 6; i++) {
    dWfdWc_x(i, i) = intp;
    dWfdWc_y(i, i) = intp;
  } 

  // Find dwdx and dwdy, the geometric weights used to obtain the 
  // face solution gradients based on the cell-centred variables.
  SolnBlk.Grid.dDiamondPathdC(dwdx, dwdy, Rii, Rjj, Orient_face, LL, RR, Orient_cell);

  dWfdWc_x( 6,0) = dwdx;
  dWfdWc_x( 7,1) = dwdx;        
  dWfdWc_x( 9,2) = dwdx;
  dWfdWc_x(11,3) = dwdx;
  dWfdWc_x(12,4) = dwdx;
  dWfdWc_x(13,5) = dwdx;

  dWfdWc_x( 8,1) = dwdy;
  dWfdWc_x(10,2) = dwdy;

  dWfdWc_y( 6,0) = dwdy;   
  dWfdWc_y( 8,1) = dwdy;
  dWfdWc_y(10,2) = dwdy;
  dWfdWc_y(11,3) = dwdy;
  dWfdWc_y(12,4) = dwdy;
  dWfdWc_y(13,5) = dwdy;

  dWfdWc_y( 7,1) = dwdx;
  dWfdWc_y( 9,2) = dwdx;
}



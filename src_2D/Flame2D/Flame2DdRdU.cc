/////////////////////////////////////////////////////////////////////
///
/// \file Flame2DdRdU.cc
/// 
/// \author Marc R.J. Charest
/// 
/// \brief This file contains the function definitions for 
///        Jacobian formulations related to Flame2DState.
///        These Jacobians are used for:
///        -> Semi-Implicit (Source only)
///        -> Newton-Krylov GMRES Preconditioner
///           (1st and 2nd order stencil Inviscid, Viscous, and Source)
///
/////////////////////////////////////////////////////////////////////
#include "Flame2DdRdU.h"

// All functions should use Uo values as semi-implicit can be used in 
// a multistage scheme and this ensures that all Jacobians are taken as dR/dUo.
// For NKS, these functions are only called after an update, so U and Uo
// are equivalent.  To speed up calculations for semi-implicit methods, 
// the conserved Uo is converted to Wo and passed down through all the
// SemiImplicitBlockJacobi() functions.  However, the dFIdW() functions
// are now based off of W as they are only used by NKS.  If one were
// to use a point-implicit method which would require dFIdW, which
// can be used in a multistage scheme, this would no longer be consistent.


/////////////////////////////////////////////////////////////////////
/// Semi-Implicit Block Jacobi formulations
/////////////////////////////////////////////////////////////////////

/********************************************************
 * Routine: SemiImplicitBlockJacobi                     *
 *                                                      *
 * This routine adds all the appropriate source         *    
 * Jacobians based on flow type, for in semi-implicit   *
 * and implicit calculations.                           *
 *                                                      * 
 ********************************************************/ 
void SemiImplicitBlockJacobi(DenseMatrix &dSdU,
			     Flame2D_Quad_Block &SolnBlk,
			     Flame2D_pState &Wo,
			     const int &ii, const int &jj){ 
  
  if (SolnBlk.Axisymmetric && SolnBlk.Flow_Type != FLOWTYPE_INVISCID) { 

    // declares
    int n = dSdU.get_n(); 
    static DenseMatrix dRdW(n,n); dRdW.zero();
    static DenseMatrix dWdU(n,n); dWdU.zero();
    
    // Add Source Jacobians (viscous axisymmetric, turbulence)
    SemiImplicitBlockJacobi_dSdW(dRdW,SolnBlk,Wo,ii,jj);
    
    // Transformation Jacobian 
    Wo.dWdU(dWdU);  // <- need to load dih/dic temporary array (done outside)
    dSdU += dRdW*dWdU;
  }

  // Add Source Jacobians (inviscid axisymmetric, chemistry, gravity)
  SemiImplicitBlockJacobi_dSdU(dSdU,SolnBlk,Wo,ii,jj);                 

}


/********************************************************
 * Routine: SemiImplicitBlockJacobi_dSdW                *
 *                                                      *
 * This routine adds all the appropriate source         *    
 * Jacobians based on flow type, for in semi-implicit   *
 * and implicit calculations.                           *
 *                                                      * 
 ********************************************************/ 
void SemiImplicitBlockJacobi_dSdW(DenseMatrix &dSdW,
				  Flame2D_Quad_Block &SolnBlk,
				  Flame2D_pState &Wo,
				  const int &ii, const int &jj){ 
  
  //
  // Viscous Axisymmetric source term jacobian    
  //
  if(SolnBlk.Axisymmetric && SolnBlk.Flow_Type != FLOWTYPE_INVISCID){     
    
    //Cacluate 2nd derivatives  
    double d_dWdx_dW_C,d_dWdy_dW_C;
    d_dWd_dW_Center(d_dWdx_dW_C,d_dWdy_dW_C,SolnBlk,ii, jj);  
    
    // Compute the jacobian
    Wo.dSa_vdW(dSdW,
	       SolnBlk.dWdx[ii][jj],
	       SolnBlk.dWdy[ii][jj],
	       SolnBlk.Grid.Cell[ii][jj].Xc,
	       SolnBlk.Axisymmetric,
	       d_dWdx_dW_C,d_dWdy_dW_C);
  }
  
}


/********************************************************
 * Routine: SemiImplicitBlockJacobi_dSdU                *
 *                                                      *
 * This routine adds all the appropriate source         *    
 * Jacobians based on flow type, for in semi-implicit   *
 * and implicit calculations.                           *
 *                                                      * 
 ********************************************************/ 
void SemiImplicitBlockJacobi_dSdU(DenseMatrix &dSdU,
				  Flame2D_Quad_Block &SolnBlk,
				  Flame2D_pState &Wo,
				  const int &ii, const int &jj){
  
  //Add Jacobian for inviscid axisymmetric source terms
  if (SolnBlk.Axisymmetric) {
    Wo.dSa_idU(dSdU,
	       SolnBlk.Grid.Cell[ii][jj].Xc, 
	       SolnBlk.Axisymmetric);
  }
  
  //Add Jacobian for finite-rate chemistry source terms  
  if (Flame2D_pState::isReacting()) Wo.dSwdU(dSdU);  

  //Add Jacobian for gravitational source terms
  if (SolnBlk.Gravity){
    Wo.dSgdU(dSdU);
  } 
  
}



/////////////////////////////////////////////////////////////////////
/// Inviscid Flux Jacobian (FI)
/////////////////////////////////////////////////////////////////////

/*********************************************************
 * Routine: Rotation_Matrix2                             *
 *                                                       *
 * This function returns either the rotation matrix, A,  *
 * or the inverse of A.                                  *
 *                                                       *
 * Note: A_matrix = 1 for returning A.                   *
 *                = 0 for returning inverse of A.        *
 *                                                       *
 *********************************************************/
/*The rotation matrix is used for the inviscid flux calculations */
void RotationMatrix(DenseMatrix &mat, const Vector2D &nface, const int A_matrix) 
{
  const double &cos_angle = nface.x; 
  const double &sin_angle = nface.y;
    
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
void dFIdW_Inviscid_ROE(DenseMatrix& dRdW, 
			Flame2D_Quad_Block &SolnBlk,  
			Flame2D_Input_Parameters &Input_Parameters,
			const int &ii, const int &jj, 
			const int Orient){
   
  // declares
  const int overlap = Input_Parameters.NKS_IP.GMRES_Overlap;
  const int NUM_VAR = dRdW.get_n();
  static Vector2D nface;
  double lface;
  static DenseMatrix dFidW(NUM_VAR, NUM_VAR);
  static DenseMatrix A(NUM_VAR, NUM_VAR);
  static DenseMatrix AI(NUM_VAR, NUM_VAR);
  Flame2D_pState *Wl, *Wr;
  static Flame2D_pState Wa;
  static Flame2D_State lambdas_a, lambdas_l, lambdas_r, wavespeeds;
  dFidW.zero();

  // GHOST CELL so do nothing
  if (ii < SolnBlk.ICl -overlap || ii > SolnBlk.ICu + overlap ||
      jj < SolnBlk.JCl -overlap || jj > SolnBlk.JCu + overlap) {
    cout<<"\n Hey I am not suppose to be here! \n"; 
    exit(-1);
  }

  // get orientation
  int Ri, Rj;
  if (Orient == NORTH) {
    Ri = ii; Rj=jj-1;
    nface = SolnBlk.Grid.nfaceN(Ri, Rj);
    lface = SolnBlk.Grid.lfaceN(Ri, Rj);
    Wl = &SolnBlk.W[Ri  ][Rj  ]; 
    Wr = &SolnBlk.W[Ri  ][Rj+1]; 
  } else if (Orient == SOUTH) {
    Ri = ii; Rj=jj+1;
    nface = SolnBlk.Grid.nfaceS(Ri, Rj);
    lface = SolnBlk.Grid.lfaceS(Ri, Rj);
    Wl = &SolnBlk.W[Ri  ][Rj  ]; 
    Wr = &SolnBlk.W[Ri  ][Rj-1];  
  } else if (Orient == EAST) { 
    Ri = ii-1; Rj=jj;
    nface = SolnBlk.Grid.nfaceE(Ri, Rj);     
    lface = SolnBlk.Grid.lfaceE(Ri, Rj);
    Wl = &SolnBlk.W[Ri  ][Rj  ]; 
    Wr = &SolnBlk.W[Ri+1][Rj  ];
  } else if (Orient == WEST) { 
    Ri = ii+1; Rj=jj;
    nface = SolnBlk.Grid.nfaceW(Ri, Rj);
    lface = SolnBlk.Grid.lfaceW(Ri, Rj);
    Wl = &SolnBlk.W[Ri  ][Rj  ]; 
    Wr = &SolnBlk.W[Ri-1][Rj  ];
  }
    
  // Compute Rotation matrices
  RotationMatrix(A, nface, 1);
  RotationMatrix(AI, nface, 0);

  //---------------------------------------------------------------
  // Compute the Jacobian
  //---------------------------------------------------------------
  // temporary storage
  double ur( ((const Flame2D_pState*)Wr)->vx() ), 
    vr( ((const Flame2D_pState*)Wr)->vy() ), 
    ul( ((const Flame2D_pState*)Wl)->vx() ), 
    vl( ((const Flame2D_pState*)Wl)->vy() );

  //Rotate left and right states in place 
  Wl->Rotate(nface);
  Wr->Rotate(nface);
     
  //Determin Roe Averaged State
  Wa.RoeAverage( *Wl, *Wr );

  // Jacobian dF/dW
  // Note: Wr will always be i,j
  Wr->dFIdW(dFidW, HALF);

  //---------------------------------------------------------------
  // Regular Roe (no preconditioning)
  //---------------------------------------------------------------
  if(!Input_Parameters.Preconditioning){

    // Determine Wave Speeds
    Wl->lambda_x(lambdas_l);
    Wr->lambda_x(lambdas_r);
    Wa.lambda_x(lambdas_a);
    wavespeeds.HartenFix_Abs( lambdas_a,
			      lambdas_l,
			      lambdas_r );

    // Compute each element of Jacobian(i,j)
    // dFidW(irow, jcol) -= HALF*wavespeeds[i]*lp[jcol+1]*rc[irow+1];
    Wa.Flux_Dissipation_Jac(dFidW, wavespeeds, HALF);

      
    //---------------------------------------------------------------
    // LOW MACH NUMBER PRECONDITIONING
    //---------------------------------------------------------------
  } else if(Input_Parameters.Preconditioning){
       
    //preconditioner spacing -> THIS MAY NOT BE CONSISTENT !!!!!!!!!!!
    double deltax( SolnBlk.delta_n(ii, jj) );

    // Roeaverage state preconditioned velocity
    double MR2a( Wa.Mr2(SolnBlk.Flow_Type,deltax) );
      
    // Determine Preconditioned Wave Speeds
    Wl->lambda_preconditioned_x(lambdas_l, Wl->Mr2(SolnBlk.Flow_Type,deltax)); 
    Wr->lambda_preconditioned_x(lambdas_r, Wr->Mr2(SolnBlk.Flow_Type,deltax));
    Wa.lambda_preconditioned_x(lambdas_a, MR2a);
    wavespeeds.HartenFix_Abs( lambdas_a, 
			      lambdas_l, 
			      lambdas_r );       
       
    // Calculate the preconditioned upwind dissipation flux.
    // dFidW(irow, jcol) -= HALF*wavespeeds[i]*lp[jcol+1]*rc[irow+1];
    static DenseMatrix Flux_dissipation(NUM_VAR,NUM_VAR);
    Flux_dissipation.zero();
    Wa.Flux_Dissipation_Jac_precon(Flux_dissipation, MR2a, wavespeeds, HALF);

    // Evaluate the low-Mach-number local preconditioner for the Roe-averaged state.
    static DenseMatrix P( NUM_VAR, NUM_VAR, ZERO);          
    Wa.Low_Mach_Number_Preconditioner(P,SolnBlk.Flow_Type,deltax);
       
    // Add preconditioned dissipation to Inviscid Jacobian
    dFidW += P*Flux_dissipation;
       
  }     
     
  //---------------------------------------------------------------
  // Add to dRdW
  //---------------------------------------------------------------
  //Rotate back 
  dFidW *= lface;
  dRdW += AI*dFidW*A;
    
  // Rotate Back -> avoid roundoff by setting the exact values
  Wl->setVelocity( ul, vl );
  Wr->setVelocity( ur, vr );
      
}

/********************************************************
 * Routine: Inviscid Roe Flux Jacobian                  *
 *                                                      *
 *     Finite Differences                               *
 *                                                      *
 ********************************************************/
void dFIdW_Inviscid_ROE_FD(DenseMatrix& dRdW, 
			   Flame2D_Quad_Block &SolnBlk,  
			   Flame2D_Input_Parameters &Input_Parameters,
			   const int &ii, const int &jj, 
			   const int Orient){
   
  // declares
  const int overlap = Input_Parameters.NKS_IP.GMRES_Overlap;
  const int NUM_VAR = dRdW.get_n();   
  static Vector2D nface;
  double lface;
  static DenseMatrix dFidW(NUM_VAR, NUM_VAR);
  static DenseMatrix A(NUM_VAR, NUM_VAR);
  static DenseMatrix AI(NUM_VAR, NUM_VAR);
  Flame2D_pState *Wl, *Wr;
  dFidW.zero();

  // GHOST CELL so do nothing
  if (ii < SolnBlk.ICl -overlap || ii > SolnBlk.ICu + overlap ||
      jj < SolnBlk.JCl -overlap || jj > SolnBlk.JCu + overlap) {
     cout<<"\n Hey I am not suppose to be here! \n";
     exit(1);
  }

  // get orientation
  int Ri, Rj;
  if (Orient == NORTH) {
    Ri = ii; Rj=jj-1;
    nface = SolnBlk.Grid.nfaceN(Ri, Rj);
    lface = SolnBlk.Grid.lfaceN(Ri, Rj);
    Wl = &SolnBlk.W[Ri  ][Rj  ]; 
    Wr = &SolnBlk.W[Ri  ][Rj+1]; 
  } else if (Orient == SOUTH) {
    Ri = ii; Rj=jj+1;
    nface = SolnBlk.Grid.nfaceS(Ri, Rj);
    lface = SolnBlk.Grid.lfaceS(Ri, Rj);
    Wl = &SolnBlk.W[Ri  ][Rj  ]; 
    Wr = &SolnBlk.W[Ri  ][Rj-1];  
  } else if (Orient == EAST) { 
    Ri = ii-1; Rj=jj;
    nface = SolnBlk.Grid.nfaceE(Ri, Rj);     
    lface = SolnBlk.Grid.lfaceE(Ri, Rj);
    Wl = &SolnBlk.W[Ri  ][Rj  ]; 
    Wr = &SolnBlk.W[Ri+1][Rj  ];
  } else if (Orient == WEST) { 
    Ri = ii+1; Rj=jj;
    nface = SolnBlk.Grid.nfaceW(Ri, Rj);
    lface = SolnBlk.Grid.lfaceW(Ri, Rj);
    Wl = &SolnBlk.W[Ri  ][Rj  ]; 
    Wr = &SolnBlk.W[Ri-1][Rj  ];
  }
    
  // Compute Rotation matrices
  RotationMatrix(A, nface, 1);
  RotationMatrix(AI, nface, 0);

  // temporary storage
  double ur( ((const Flame2D_pState*)Wr)->vx() ), 
    vr( ((const Flame2D_pState*)Wr)->vy() ), 
    ul( ((const Flame2D_pState*)Wl)->vx() ), 
    vl( ((const Flame2D_pState*)Wl)->vy() );

  //Rotate left and right states in place 
  Wl->Rotate(nface);
  Wr->Rotate(nface);

  //---------------------------------------------------------------
  // Compute the Jacobian
  //---------------------------------------------------------------
  //Jacobain using Finite Differences
  static Flame2D_State FluxA, FluxB; 
  static Flame2D_pState WA, WB;
  const double perturb = 5e-6;
  double a;
     
  //For Preconditioning
  double delta_n = SolnBlk.delta_n(ii, jj);
  
  for(int jcol=0; jcol<(NUM_VAR); jcol++){
    WA = *Wr;
    WB = *Wr;
    
    if( jcol <NUM_FLAME2D_VAR_SANS_SPECIES) {
      WA[jcol+1] += perturb*max(ONE,(*Wr)[jcol+1]); 	 
      WB[jcol+1] -= perturb*max(ONE,(*Wr)[jcol+1]); 
    } else {
      a =  perturb*max(ONE,(*Wr)[jcol+1]); 
      WA[jcol+1] += a;
      WA[NUM_VAR+1] -= a;      
      WB[jcol+1] -= a;
      WB[NUM_VAR+1] += a;
    }
    
    FluxA.FluxRoe_x(*Wl,WA, Input_Parameters.Preconditioning,SolnBlk.Flow_Type,delta_n);       
    FluxB.FluxRoe_x(*Wl,WB, Input_Parameters.Preconditioning,SolnBlk.Flow_Type,delta_n);
    
    for(int irow=0; irow<(NUM_VAR); irow++){
      dFidW(irow,jcol) = (FluxA[irow+1] - FluxB[irow+1])/(TWO*perturb*max(ONE,(*Wr)[jcol+1]));
    }
  } 

  //---------------------------------------------------------------
  // Add to dRdW
  //---------------------------------------------------------------
  //Rotate back 
  dFidW *= lface;
  dRdW += AI*dFidW*A;
    
  // Rotate Back -> avoid roundoff by setting the exact values
  Wl->setVelocity( ul, vl );
  Wr->setVelocity( ur, vr );

}

/********************************************************
 * Routine: Inviscid Flux Jacobian using AUSM_plus_up   *
 *                                                      *
 * This routine returns the inviscid components of      *    
 * Jacobian matrix for the specified local solution     *
 * block calculated analytically.                       *
 *                                                      *
 * dF/dW_R                                              *
 ********************************************************/
void dFIdW_Inviscid_AUSM_plus_up(DenseMatrix& dRdW, 
				 Flame2D_Quad_Block &SolnBlk,  
				 Flame2D_Input_Parameters &Input_Parameters,
				 const int &ii, const int &jj, 
				 const int Orient){
   
  // declares
  const int overlap = Input_Parameters.NKS_IP.GMRES_Overlap;
  const int NUM_VAR = dRdW.get_n();
  static Vector2D nface;
  double lface;
  static DenseMatrix dFidW(NUM_VAR, NUM_VAR);
  static DenseMatrix A(NUM_VAR, NUM_VAR);
  static DenseMatrix AI(NUM_VAR, NUM_VAR);
  dFidW.zero();

  // GHOST CELL so do nothing
  if (ii < SolnBlk.ICl -overlap || ii > SolnBlk.ICu + overlap ||
      jj < SolnBlk.JCl -overlap || jj > SolnBlk.JCu + overlap) {
    cout<<"\n Hey I am not suppose to be here! \n";
    exit(-1);
  }

  // get orientation
  int Ri, Rj;
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
     
  // Compute Rotation matrices
  RotationMatrix(A, nface, 1);
  RotationMatrix(AI, nface, 0);
          
  //---------------------------------------------------------------
  // Compute the Jacobian
  //---------------------------------------------------------------

  // temporarily store velocity
  double u( ((const Flame2D_pState&)SolnBlk.W[ii][jj]).vx() ), 
    v( ((const Flame2D_pState&)SolnBlk.W[ii][jj]).vy() );

  // rotate in place
  SolnBlk.W[ii][jj].Rotate(nface);

  // Compute Jac
  SolnBlk.W[ii][jj].dFIdW(dFidW, HALF);

  // rotate back with exact values
  SolnBlk.W[ii][jj].setVelocity(u, v);

  //---------------------------------------------------------------
  // FILL IN AUSM SPECIFIC STUFF HERE
  //---------------------------------------------------------------
  cerr << "dFIdW_Inviscid_AUSM_plus_up(): Not finished yet lazy ass!!!!";
  exit(-1);

  //---------------------------------------------------------------
  // Add to dRdW
  //---------------------------------------------------------------
  //Rotate back 
  dFidW *= lface;
  dRdW += AI*dFidW*A;
        
}

/////////////////////////////////////////////////////////////////////
/// Viscous Flux Jacobian (Fv)
/////////////////////////////////////////////////////////////////////

/********************************************************
 * Routine: dFvdWf_Diamond                              *
 *                                                      *
 * This routine calculates the Viscous components of    *
 * the residual with respect to the primitive variables *   
 * in based on Green Gauss Diamond Path Reconstruction. *
 *                                                      *
 ********************************************************/
void dFvdWf_Diamond(DenseMatrix &dFvdWf, 
		    DenseMatrix &dGvdWf, 
		    Flame2D_Quad_Block &SolnBlk,
		    const int &Orient, 
		    const int &ii, const int &jj){
   
  // declares
  static Flame2D_pState QuadraturePoint_W;
  const Flame2D_State *QuadraturePoint_dWdx, *QuadraturePoint_dWdy;
  double radius(0.0);
  
  //---------------------------------------------------------------
  // Compute Gradients
  //---------------------------------------------------------------
  switch(Orient){

    /****************************** NORTH ******************************/
  case NORTH:
    // QuadraturePoint_W = HALF*(SolnBlk.UnoNW(ii, jj).W() +SolnBlk.UnoNE(ii, jj).W() );
    QuadraturePoint_W.Average( SolnBlk.Wnd[ii][jj+1], SolnBlk.Wnd[ii+1][jj+1] );
    QuadraturePoint_dWdx = &SolnBlk.dWdx_faceN[ii][jj];
    QuadraturePoint_dWdy = &SolnBlk.dWdy_faceN[ii][jj];

    if (SolnBlk.Axisymmetric == AXISYMMETRIC_X) {
      radius = (SolnBlk.Grid.xfaceN(ii,jj).x < MICRO) ? MICRO : SolnBlk.Grid.xfaceN(ii,jj).x;     
    } else if (SolnBlk.Axisymmetric == AXISYMMETRIC_Y) {    
      radius = (SolnBlk.Grid.xfaceN(ii,jj).y < MICRO) ? MICRO : SolnBlk.Grid.xfaceN(ii,jj).y;            
    } 

    break;

    /****************************** EAST ******************************/
  case EAST:
    // QuadraturePoint_W = HALF*(SolnBlk.UnoNE(ii, jj).W() + SolnBlk.UnoSE(ii, jj).W() );
    QuadraturePoint_W.Average( SolnBlk.Wnd[ii+1][jj+1], SolnBlk.Wnd[ii+1][jj] );
    QuadraturePoint_dWdx = &SolnBlk.dWdx_faceE[ii][jj];
    QuadraturePoint_dWdy = &SolnBlk.dWdy_faceE[ii][jj];

    if (SolnBlk.Axisymmetric == AXISYMMETRIC_X) {
      radius = ( SolnBlk.Grid.xfaceE(ii,jj).x < MICRO) ? MICRO : SolnBlk.Grid.xfaceE(ii,jj).x;     
    } else if (SolnBlk.Axisymmetric == AXISYMMETRIC_Y) {    
      radius = ( SolnBlk.Grid.xfaceE(ii,jj).y < MICRO) ? MICRO : SolnBlk.Grid.xfaceE(ii,jj).y;            
    } 

    break;

    /****************************** SOUTH ******************************/
  case SOUTH:      
    // QuadraturePoint_W = HALF*(SolnBlk.UnoSE(ii, jj).W() +SolnBlk.UnoSW(ii, jj).W() );
    QuadraturePoint_W.Average( SolnBlk.Wnd[ii+1][jj], SolnBlk.Wnd[ii][jj] );
    QuadraturePoint_dWdx = &SolnBlk.dWdx_faceS[ii][jj];
    QuadraturePoint_dWdy = &SolnBlk.dWdy_faceS[ii][jj];
     
    if (SolnBlk.Axisymmetric == AXISYMMETRIC_X) {
      radius = (SolnBlk.Grid.xfaceS(ii,jj).x < MICRO) ? MICRO : SolnBlk.Grid.xfaceS(ii,jj).x;     
    } else if (SolnBlk.Axisymmetric == AXISYMMETRIC_Y) {    
      radius = (SolnBlk.Grid.xfaceS(ii,jj).y < MICRO) ? MICRO : SolnBlk.Grid.xfaceS(ii,jj).y;            
    } 
    break;

    /****************************** WEST ******************************/
  case WEST:       
    // QuadraturePoint_W = HALF*(SolnBlk.UnoNW(ii, jj).W() +SolnBlk.UnoSW(ii, jj).W());
    QuadraturePoint_W.Average( SolnBlk.Wnd[ii][jj+1], SolnBlk.Wnd[ii][jj] );
    QuadraturePoint_dWdx = &SolnBlk.dWdx_faceW[ii][jj];
    QuadraturePoint_dWdy = &SolnBlk.dWdy_faceW[ii][jj];
     
    if (SolnBlk.Axisymmetric == AXISYMMETRIC_X) {
      radius = (SolnBlk.Grid.xfaceW(ii,jj).x < MICRO) ? MICRO : SolnBlk.Grid.xfaceW(ii,jj).x;     
    } else if (SolnBlk.Axisymmetric == AXISYMMETRIC_Y) {    
      radius = (SolnBlk.Grid.xfaceW(ii,jj).y < MICRO) ? MICRO : SolnBlk.Grid.xfaceW(ii,jj).y;            
    } 
    break;
  }

  //---------------------------------------------------------------
  // Compute Jacobian terms
  //---------------------------------------------------------------
  QuadraturePoint_W.dFvdWf_dGvdWf( dFvdWf, dGvdWf, 
				   *QuadraturePoint_dWdx, 
				   *QuadraturePoint_dWdy, 
				   SolnBlk.Axisymmetric, radius );
    
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
void dWfdWc_Diamond(DenseMatrix &dWfdWc_x,
		    DenseMatrix &dWfdWc_y, 
		    Flame2D_Quad_Block &SolnBlk,
		    const int &Orient_face, 
		    const int &i, const int &j, 
		    const int &Orient_cell){
   
  //---------------------------------------------------------------
  // All these confusing relationships are based on an outward
  // facing normal from the Orient_face ie  NORTH face ,
  // left is NW and right is NE.  
  //---------------------------------------------------------------
  double LL(ZERO),RR(ZERO); 
  switch (Orient_face){
  case NORTH:
    if(Orient_cell == CENTER){    
      LL = SolnBlk.dWn_dWc(i,j+1, NORTH_WEST);
      RR = SolnBlk.dWn_dWc(i+1, j+1, NORTH_EAST);
    } else if(Orient_cell == NORTH){     
      LL = SolnBlk.dWn_dWc(i,j+1, SOUTH_WEST);
      RR = SolnBlk.dWn_dWc(i+1, j+1, SOUTH_EAST);
    } else if(Orient_cell == EAST){
      LL = ZERO;
      RR = SolnBlk.dWn_dWc(i+1,j+1,NORTH_WEST);                  
    } else if(Orient_cell == WEST){ 
      LL = SolnBlk.dWn_dWc(i,j+1,NORTH_EAST);      
      RR = ZERO;                   
    } else if(Orient_cell == NORTH_WEST){
      LL = SolnBlk.dWn_dWc(i,j+1,SOUTH_EAST);      
      RR = ZERO;   
    } else if(Orient_cell == NORTH_EAST){
      LL = ZERO;
      RR = SolnBlk.dWn_dWc(i+1,j+1,SOUTH_WEST);    
    }
    break;
      
  case EAST:
    if(Orient_cell == CENTER){  
      LL = SolnBlk.dWn_dWc(i+1, j+1, NORTH_EAST);
      RR = SolnBlk.dWn_dWc(i+1, j, SOUTH_EAST); 
    } else if(Orient_cell == EAST){    
      LL = SolnBlk.dWn_dWc(i+1, j+1, NORTH_WEST);
      RR = SolnBlk.dWn_dWc(i+1, j, SOUTH_WEST); 
    } else if(Orient_cell == NORTH){     
      LL = SolnBlk.dWn_dWc(i+1,j+1,SOUTH_EAST);
      RR = ZERO;
    } else if(Orient_cell == SOUTH){           
      LL = ZERO;
      RR = SolnBlk.dWn_dWc(i+1,j,NORTH_EAST);       
    } else if(Orient_cell == NORTH_EAST){     
      LL = SolnBlk.dWn_dWc(i+1,j+1,SOUTH_WEST);
      RR = ZERO;
    } else if(Orient_cell == SOUTH_EAST){           
      LL = ZERO;
      RR = SolnBlk.dWn_dWc(i+1,j,NORTH_WEST);       
    }
    break;
  
  case SOUTH:
    if(Orient_cell == CENTER){  
      LL = SolnBlk.dWn_dWc(i+1,j,SOUTH_EAST);
      RR = SolnBlk.dWn_dWc(i,j, SOUTH_WEST);
    } else if(Orient_cell == SOUTH){     
      LL = SolnBlk.dWn_dWc(i+1,j, NORTH_EAST);
      RR = SolnBlk.dWn_dWc(i,j, NORTH_WEST);
    } else if(Orient_cell == EAST){ 
      LL = SolnBlk.dWn_dWc(i+1,j, SOUTH_WEST);      
      RR = ZERO;
    } else if(Orient_cell == WEST){    
      LL = ZERO; 
      RR = SolnBlk.dWn_dWc(i,j, SOUTH_EAST);                              
    } else if(Orient_cell == SOUTH_EAST){ 
      LL = SolnBlk.dWn_dWc(i+1,j, NORTH_WEST);      
      RR = ZERO;
    } else if(Orient_cell == SOUTH_WEST){    
      LL = ZERO; 
      RR = SolnBlk.dWn_dWc(i,j, NORTH_EAST);      
    }
    break;
   
  case WEST:
    if(Orient_cell == CENTER){  
      LL = SolnBlk.dWn_dWc(i,j, SOUTH_WEST);
      RR = SolnBlk.dWn_dWc(i,j+1,NORTH_WEST);  
    } else if(Orient_cell == WEST){  
      LL = SolnBlk.dWn_dWc(i,j, SOUTH_EAST);
      RR = SolnBlk.dWn_dWc(i,j+1,NORTH_EAST);  
    } else if(Orient_cell == SOUTH){           
      LL = SolnBlk.dWn_dWc(i,j,NORTH_WEST);
      RR = ZERO;
    } else if(Orient_cell == NORTH){   
      LL = ZERO;
      RR = SolnBlk.dWn_dWc(i,j+1,SOUTH_WEST);
    } else if(Orient_cell == SOUTH_WEST){           
      LL = SolnBlk.dWn_dWc(i,j,NORTH_EAST);
      RR = ZERO;
    } else if(Orient_cell == NORTH_WEST){   
      LL = ZERO;
      RR = SolnBlk.dWn_dWc(i,j+1,SOUTH_EAST);
    }
    break;     
  }

  //get 2nd derivatives
  double d_dWdx_dW(ZERO), d_dWdy_dW(ZERO);
  d_dWd_dW_Diamond(d_dWdx_dW,d_dWdy_dW,SolnBlk,LL,RR,Orient_cell,Orient_face,i,j);
 
  /*********************** X - DIRECTION **************************************/
  for(int nn=0; nn<6; nn++){
    dWfdWc_x(nn,nn) = HALF*(LL+RR);
    dWfdWc_y(nn,nn) = dWfdWc_x(nn,nn);
  } 
 
  dWfdWc_x(6,0) = d_dWdx_dW;
  dWfdWc_x(7,1) = d_dWdx_dW;        //NOTE 7,1 & 8,1 same for X and Y 
  dWfdWc_x(8,1) = d_dWdy_dW;
  dWfdWc_x(9,2) =  dWfdWc_x(7,1);
  dWfdWc_x(10,2) =  dWfdWc_x(8,1);
  dWfdWc_x(11,3) = d_dWdx_dW;
  dWfdWc_x(12,4) = d_dWdx_dW;
  dWfdWc_x(13,5) = d_dWdx_dW;
   
  const int ns = Flame2D_State::NumSpecies() - Flame2D_State::NSm1;
  for(int Num=0; Num<(ns); Num++){
    dWfdWc_x(14+Num,6+Num) =  d_dWdx_dW;
    dWfdWc_y(14+Num,6+Num) =  d_dWdy_dW;
  }

  /*********************** Y - DIRECTION **************************************/
  dWfdWc_y(6,0) = d_dWdy_dW;   
  dWfdWc_y(7,1) = d_dWdx_dW;
  dWfdWc_y(8,1) = d_dWdy_dW;
  dWfdWc_y(9,2) =  dWfdWc_y(7,1);
  dWfdWc_y(10,2) =  dWfdWc_y(8,1);
  dWfdWc_y(11,3) = d_dWdy_dW;
  dWfdWc_y(12,4) = d_dWdy_dW;
  dWfdWc_y(13,5) = d_dWdy_dW;
     
}

/////////////////////////////////////////////////////////////////////
/// 2nd deriavaites
/////////////////////////////////////////////////////////////////////

/********************************************************
 * Routine: d_dWd_dW_Diamond                            *
 *                                                      *
 * This routine calculates the 2nd deriavaites          *
 * associated with diamond path and bilinear            *
 * interpolation.                                       *
 *                                                      *
 ********************************************************/
void d_dWd_dW_Diamond(double &d_dWdx_dW, double &d_dWdy_dW, 
		      Flame2D_Quad_Block &SolnBlk, 
		      const double &LEFT, const double &RIGHT, 
		      const int &Orient_cell, const int &Orient_face,  
		      const int &i, const int &j){

  //  double area[4];
  double AREA;
  Vector2D norm[4];
  double  dWnNWdWc, dWnNEdWc,  dWnSWdWc, dWnSEdWc;
 
  switch(Orient_face){
    /*************** NORTH ****************************/
  case NORTH: 
    dWnNWdWc = LEFT;
    dWnNEdWc = RIGHT;

    //  normal vector of the SE side of a diamond 
    norm[0].x = SolnBlk.Grid.nodeNE(i,j).X.y - SolnBlk.Grid.Cell[i][j].Xc.y;
    norm[0].y = -( SolnBlk.Grid.nodeNE(i,j).X.x - SolnBlk.Grid.Cell[i][j].Xc.x);
    //  normal vector of the NE side of a diamond 
    norm[1].x = SolnBlk.Grid.Cell[i][j+1].Xc.y - SolnBlk.Grid.nodeNE(i,j).X.y;
    norm[1].y = -(SolnBlk.Grid.Cell[i][j+1].Xc.x - SolnBlk.Grid.nodeNE(i,j).X.x);
    //  normal vector of the NW side of a diamond 
    norm[2].x =  SolnBlk.Grid.nodeNW(i,j).X.y - SolnBlk.Grid.Cell[i][j+1].Xc.y;
    norm[2].y = -(SolnBlk.Grid.nodeNW(i,j).X.x - SolnBlk.Grid.Cell[i][j+1].Xc.x);
    //  normal vector of the SW side of a diamond 
    norm[3].x = SolnBlk.Grid.Cell[i][j].Xc.y - SolnBlk.Grid.nodeNW(i,j).X.y ;
    norm[3].y = -( SolnBlk.Grid.Cell[i][j].Xc.x - SolnBlk.Grid.nodeNW(i,j).X.x);
    
    AREA =  HALF*(fabs((SolnBlk.Grid.nodeNE(i,j).X-SolnBlk.Grid.Cell[i][j].Xc)^
		       (SolnBlk.Grid.nodeNW(i,j).X-SolnBlk.Grid.Cell[i][j].Xc)) +
		  fabs((SolnBlk.Grid.nodeNW(i,j).X-SolnBlk.Grid.Cell[i][j+1].Xc)^
		       (SolnBlk.Grid.nodeNE(i,j).X-SolnBlk.Grid.Cell[i][j+1].Xc)));
    
    if(Orient_cell == CENTER){
      d_dWdx_dW = HALF*((ONE+dWnNEdWc)* norm[0].x+dWnNEdWc* norm[1].x+ dWnNWdWc* norm[2].x+ (ONE+dWnNWdWc)* norm[3].x)/AREA;
      d_dWdy_dW = HALF*((ONE+dWnNEdWc)* norm[0].y+dWnNEdWc* norm[1].y+ dWnNWdWc* norm[2].y+ (ONE+dWnNWdWc)* norm[3].y)/AREA;  
    } else if( Orient_cell == NORTH) {
      d_dWdx_dW = HALF*(dWnNEdWc*norm[0].x + (ONE+dWnNEdWc)* norm[1].x + (ONE+dWnNWdWc)* norm[2].x + dWnNWdWc*norm[3].x)/AREA;
      d_dWdy_dW = HALF*(dWnNEdWc*norm[0].y + (ONE+dWnNEdWc)* norm[1].y + (ONE+dWnNWdWc)* norm[2].y + dWnNWdWc*norm[3].y)/AREA;  
    } else if( Orient_cell == EAST || Orient_cell == NORTH_EAST) {
      d_dWdx_dW = HALF*( dWnNEdWc*(norm[0].x+norm[1].x))/AREA;
      d_dWdy_dW = HALF*( dWnNEdWc*(norm[0].y+norm[1].y))/AREA;  
    } else if( Orient_cell == WEST || Orient_cell == NORTH_WEST) {
      d_dWdx_dW = HALF*( dWnNWdWc*(norm[2].x+norm[3].x))/AREA;
      d_dWdy_dW = HALF*( dWnNWdWc*(norm[2].y+norm[3].y))/AREA;  
    }

    break;
    
    /*************** EAST ****************************/
  case EAST:
    dWnNEdWc = LEFT;
    dWnSEdWc = RIGHT; 

    //  normal vector of the SE side of a diamond 
    norm[0].x =  SolnBlk.Grid.Cell[i+1][j].Xc.y - SolnBlk.Grid.nodeSE(i,j).X.y;
    norm[0].y = -(SolnBlk.Grid.Cell[i+1][j].Xc.x - SolnBlk.Grid.nodeSE(i,j).X.x);
    //  normal vector of the NE side of a diamond 
    norm[1].x = SolnBlk.Grid.nodeNE(i,j).X.y -  SolnBlk.Grid.Cell[i+1][j].Xc.y ;
    norm[1].y = -(SolnBlk.Grid.nodeNE(i,j).X.x -  SolnBlk.Grid.Cell[i+1][j].Xc.x );
    //  normal vector of the NW side of a diamond 
    norm[2].x =   SolnBlk.Grid.Cell[i][j].Xc.y - SolnBlk.Grid.nodeNE(i,j).X.y ;
    norm[2].y = -(SolnBlk.Grid.Cell[i][j].Xc.x - SolnBlk.Grid.nodeNE(i,j).X.x);
    //  normal vector of the SW side of a diamond 
    norm[3].x = SolnBlk.Grid.nodeSE(i,j).X.y - SolnBlk.Grid.Cell[i][j].Xc.y;
    norm[3].y = -(SolnBlk.Grid.nodeSE(i,j).X.x - SolnBlk.Grid.Cell[i][j].Xc.x);
    
    AREA =  HALF*(fabs((SolnBlk.Grid.nodeNE(i,j).X-SolnBlk.Grid.Cell[i+1][j].Xc)^
		       (SolnBlk.Grid.nodeSE(i,j).X-SolnBlk.Grid.Cell[i+1][j].Xc)) +
		  fabs((SolnBlk.Grid.nodeSE(i,j).X-SolnBlk.Grid.Cell[i][j].Xc)^
		       (SolnBlk.Grid.nodeNE(i,j).X-SolnBlk.Grid.Cell[i][j].Xc)));
   
    if(Orient_cell == CENTER){
      d_dWdx_dW = HALF*(dWnSEdWc* norm[0].x+ dWnNEdWc* norm[1].x+ (ONE+ dWnNEdWc)* norm[2].x+ (ONE+dWnSEdWc)* norm[3].x)/AREA;
      d_dWdy_dW = HALF*(dWnSEdWc* norm[0].y+ dWnNEdWc* norm[1].y+ (ONE+ dWnNEdWc)* norm[2].y+ (ONE+dWnSEdWc)* norm[3].y)/AREA;  
    } else if( Orient_cell == EAST) {
      d_dWdx_dW = HALF*( (ONE+dWnSEdWc)* norm[0].x + (ONE+dWnNEdWc)*norm[1].x + dWnNEdWc* norm[2].x + dWnSEdWc*norm[3].x)/AREA;
      d_dWdy_dW = HALF*( (ONE+dWnSEdWc)* norm[0].y + (ONE+dWnNEdWc)*norm[1].y + dWnNEdWc* norm[2].y + dWnSEdWc*norm[3].y)/AREA;  
    } else if( Orient_cell == NORTH || Orient_cell == NORTH_EAST) {
      d_dWdx_dW = HALF*( dWnNEdWc*(norm[1].x+norm[2].x))/AREA;
      d_dWdy_dW = HALF*( dWnNEdWc*(norm[1].y+norm[2].y))/AREA;  
    } else if( Orient_cell == SOUTH || Orient_cell == SOUTH_EAST) {
      d_dWdx_dW = HALF*( dWnSEdWc*(norm[0].x+norm[3].x))/AREA;
      d_dWdy_dW = HALF*( dWnSEdWc*(norm[0].y+norm[3].y))/AREA;  
    }
    break;

    /*************** SOUTH ****************************/
  case SOUTH:
    dWnSEdWc = LEFT;
    dWnSWdWc = RIGHT;

    //  normal vector of the SE side of a diamond 
    norm[0].x =  SolnBlk.Grid.nodeSE(i,j).X.y - SolnBlk.Grid.Cell[i][j-1].Xc.y;
    norm[0].y = -(SolnBlk.Grid.nodeSE(i,j).X.x - SolnBlk.Grid.Cell[i][j-1].Xc.x  );
    //  normal vector of the NE side of a diamond 
    norm[1].x = SolnBlk.Grid.Cell[i][j].Xc.y -  SolnBlk.Grid.nodeSE(i,j).X.y;
    norm[1].y = -(SolnBlk.Grid.Cell[i][j].Xc.x -  SolnBlk.Grid.nodeSE(i,j).X.x);
    //  normal vector of the NW side of a diamond 
    norm[2].x =   SolnBlk.Grid.nodeSW(i,j).X.y - SolnBlk.Grid.Cell[i][j].Xc.y ;
    norm[2].y = -(SolnBlk.Grid.nodeSW(i,j).X.x - SolnBlk.Grid.Cell[i][j].Xc.x);
    //  normal vector of the SW side of a diamond 
    norm[3].x =  SolnBlk.Grid.Cell[i][j-1].Xc.y - SolnBlk.Grid.nodeSW(i,j).X.y;
    norm[3].y = -(SolnBlk.Grid.Cell[i][j-1].Xc.x- SolnBlk.Grid.nodeSW(i,j).X.x);
    
    AREA =  HALF*(fabs((SolnBlk.Grid.nodeSE(i,j).X-SolnBlk.Grid.Cell[i][j-1].Xc)^
		       (SolnBlk.Grid.nodeSW(i,j).X-SolnBlk.Grid.Cell[i][j-1].Xc)) +
		  fabs((SolnBlk.Grid.nodeSE(i,j).X-SolnBlk.Grid.Cell[i][j].Xc)^
		       (SolnBlk.Grid.nodeSW(i,j).X-SolnBlk.Grid.Cell[i][j].Xc)));
    if(Orient_cell == CENTER){
      d_dWdx_dW = HALF*(dWnSEdWc* norm[0].x+ (ONE+dWnSEdWc)* norm[1].x+ (ONE+ dWnSWdWc)* norm[2].x+ (dWnSWdWc)* norm[3].x)/AREA;
      d_dWdy_dW = HALF*(dWnSEdWc* norm[0].y+ (ONE+dWnSEdWc)* norm[1].y+ (ONE+ dWnSWdWc)* norm[2].y+ (dWnSWdWc)* norm[3].y)/AREA;  
    } else if( Orient_cell == SOUTH) {
      d_dWdx_dW = HALF*( (ONE+dWnSEdWc)*norm[0].x + dWnSEdWc*norm[1].x + dWnSWdWc*norm[2].x + (ONE+dWnSWdWc)* norm[3].x)/AREA;
      d_dWdy_dW = HALF*( (ONE+dWnSEdWc)*norm[0].y + dWnSEdWc*norm[1].y + dWnSWdWc*norm[2].y + (ONE+dWnSWdWc)* norm[3].y)/AREA;  
    } else if( Orient_cell == EAST || Orient_cell == SOUTH_EAST) {
      d_dWdx_dW = HALF*( dWnSEdWc*(norm[0].x+norm[1].x))/AREA;
      d_dWdy_dW = HALF*( dWnSEdWc*(norm[0].y+norm[1].y))/AREA;  
    } else if( Orient_cell == WEST || Orient_cell == SOUTH_WEST) {
      d_dWdx_dW = HALF*( dWnSWdWc*(norm[2].x+norm[3].x))/AREA;
      d_dWdy_dW = HALF*( dWnSWdWc*(norm[2].y+norm[3].y))/AREA;  
    }
    break;
    
    /*************** WEST ****************************/
  case WEST:
    dWnSWdWc = LEFT;
    dWnNWdWc = RIGHT;
    
    //  normal vector of the SE side of a diamond 
    norm[0].x =   SolnBlk.Grid.Cell[i][j].Xc.y - SolnBlk.Grid.nodeSW(i,j).X.y;
    norm[0].y = -(SolnBlk.Grid.Cell[i][j].Xc.x - SolnBlk.Grid.nodeSW(i,j).X.x);
    //  normal vector of the NE side of a diamond 
    norm[1].x =  SolnBlk.Grid.nodeNW(i,j).X.y - SolnBlk.Grid.Cell[i][j].Xc.y;
    norm[1].y = -(SolnBlk.Grid.nodeNW(i,j).X.x - SolnBlk.Grid.Cell[i][j].Xc.x);
    //  normal vector of the NW side of a diamond 
    norm[2].x =  SolnBlk.Grid.Cell[i-1][j].Xc.y - SolnBlk.Grid.nodeNW(i,j).X.y ;
    norm[2].y = -( SolnBlk.Grid.Cell[i-1][j].Xc.x - SolnBlk.Grid.nodeNW(i,j).X.x);
    //  normal vector of the SW side of a diamond 
    norm[3].x =  SolnBlk.Grid.nodeSW(i,j).X.y - SolnBlk.Grid.Cell[i-1][j].Xc.y;
    norm[3].y = -(SolnBlk.Grid.nodeSW(i,j).X.x - SolnBlk.Grid.Cell[i-1][j].Xc.x);
    
    AREA =  HALF*(fabs((SolnBlk.Grid.nodeNW(i,j).X-SolnBlk.Grid.Cell[i][j].Xc)^
		       (SolnBlk.Grid.nodeSW(i,j).X-SolnBlk.Grid.Cell[i][j].Xc)) +
		  fabs((SolnBlk.Grid.nodeNW(i,j).X-SolnBlk.Grid.Cell[i-1][j].Xc)^
		       (SolnBlk.Grid.nodeSW(i,j).X-SolnBlk.Grid.Cell[i-1][j].Xc)));
    
    if(Orient_cell == CENTER){
      d_dWdx_dW = HALF*((ONE+dWnSWdWc)* norm[0].x+ (ONE+dWnNWdWc)* norm[1].x+ dWnNWdWc* norm[2].x+ dWnSWdWc* norm[3].x)/AREA;
      d_dWdy_dW = HALF*((ONE+dWnSWdWc)* norm[0].y+ (ONE+dWnNWdWc)* norm[1].y+ dWnNWdWc* norm[2].y+ dWnSWdWc* norm[3].y)/AREA;  
    }  else if( Orient_cell == WEST) {
      d_dWdx_dW = HALF*(dWnSWdWc*norm[0].x + dWnNWdWc*norm[1].x + (ONE+dWnNWdWc)* norm[2].x+ (ONE+dWnSWdWc)* norm[3].x)/AREA;
      d_dWdy_dW = HALF*(dWnSWdWc*norm[0].y + dWnNWdWc*norm[1].y + (ONE+dWnNWdWc)* norm[2].y+ (ONE+dWnSWdWc)* norm[3].y)/AREA;    
    } else if( Orient_cell == NORTH || Orient_cell == NORTH_WEST) {
      d_dWdx_dW = HALF*( dWnNWdWc*(norm[1].x+norm[2].x))/AREA;
      d_dWdy_dW = HALF*( dWnNWdWc*(norm[1].y+norm[2].y))/AREA;  
    } else if( Orient_cell == SOUTH || Orient_cell == SOUTH_WEST) {
      d_dWdx_dW = HALF*( dWnSWdWc*(norm[0].x+norm[3].x))/AREA;
      d_dWdy_dW = HALF*( dWnSWdWc*(norm[0].y+norm[3].y))/AREA;      
    }
    break;

  }

}

/********************************************************
 * Routine: d_dWd_dW_Center                             *
 *                                                      *
 * This routine calculates the 2nd deriavaites          *
 * associated with diamond path and bilinear            *
 * interpolation.                                       *
 *                                                      *
 ********************************************************/
void d_dWd_dW_Center(double &d_dWdx_dW_C, double &d_dWdy_dW_C, 
		     Flame2D_Quad_Block &SolnBlk, 
		     const int &i, const int &j){

  double area[4], d_dWdx_dW[4], d_dWdy_dW[4];

  // area weighted gradients at cell centers, 4 inside triangles
  area[0] = HALF*(SolnBlk.Grid.Node[i+1][j+1].X - SolnBlk.Grid.Node[i][j+1].X )^
    (SolnBlk.Grid.xfaceN(i,j)- SolnBlk.Grid.Cell[i][j].Xc);
  area[1] = HALF*(SolnBlk.Grid.Node[i+1][j+1].X - SolnBlk.Grid.Node[i+1][j].X )^
    (SolnBlk.Grid.Cell[i][j].Xc - SolnBlk.Grid.xfaceE(i,j)); 
  area[2] = HALF*(SolnBlk.Grid.Node[i+1][j].X - SolnBlk.Grid.Node[i][j].X )^
    (SolnBlk.Grid.Cell[i][j].Xc - SolnBlk.Grid.xfaceS(i,j) );  
  area[3] = HALF*(SolnBlk.Grid.Node[i][j+1].X - SolnBlk.Grid.Node[i][j].X )^
    ( SolnBlk.Grid.xfaceW(i, j) - SolnBlk.Grid.Cell[i][j].Xc );
  

  //NORTH
  d_dWd_dW_Diamond(d_dWdx_dW[0] ,d_dWdy_dW[0], SolnBlk,
		   SolnBlk.dWn_dWc(i,j+1, NORTH_WEST), SolnBlk.dWn_dWc(i+1, j+1, NORTH_EAST), 
		   CENTER, NORTH, i, j);  
  //EAST
  d_dWd_dW_Diamond(d_dWdx_dW[1] ,d_dWdy_dW[1], SolnBlk,
		   SolnBlk.dWn_dWc(i+1,j+1, NORTH_EAST), SolnBlk.dWn_dWc(i+1, j, SOUTH_EAST), 
		   CENTER, EAST, i, j);
  //SOUTH
  d_dWd_dW_Diamond(d_dWdx_dW[2] ,d_dWdy_dW[2], SolnBlk,
		   SolnBlk.dWn_dWc(i+1,j, SOUTH_EAST), SolnBlk.dWn_dWc(i, j, SOUTH_WEST), 
		   CENTER, SOUTH, i, j);
  //WEST
  d_dWd_dW_Diamond(d_dWdx_dW[3] ,d_dWdy_dW[3], SolnBlk,
		   SolnBlk.dWn_dWc(i,j, SOUTH_WEST), SolnBlk.dWn_dWc(i, j+1, NORTH_WEST), 
		   CENTER, WEST, i, j);
  
  //2nd derivative's at cell center 
  d_dWdx_dW_C = (d_dWdx_dW[0]*area[0] + d_dWdx_dW[1]*area[1] +
		 d_dWdx_dW[2]*area[2] + d_dWdx_dW[3]*area[3])/SolnBlk.Grid.Cell[i][j].A; 
  
  d_dWdy_dW_C = (d_dWdy_dW[0]*area[0]+d_dWdy_dW[1]*area[1] +
		 d_dWdy_dW[2]*area[2] + d_dWdy_dW[3]*area[3])/SolnBlk.Grid.Cell[i][j].A; 
      

}


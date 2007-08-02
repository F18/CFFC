#ifndef _CHEM2D_dRdU_INCLUDED
#include "dRdU.h"
#endif 


DenseMatrix PointImplicitBlockJacobi(Chem2D_Quad_Block &SolnBlk,
				     Chem2D_Input_Parameters &Input_Parameters,
				     const int &ii, const int &jj){
  
  int NUM_VAR_CHEM2D =  SolnBlk.NumVar() -1; 
  DenseMatrix PIBJ(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D);
  PIBJ.zero();
  //Inviscid
  DenseMatrix dFIdU(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D);
  dFIdU.zero();
  
  dFIdU = dFIdU_Inviscid(SolnBlk, Input_Parameters, ii, jj);
  PIBJ =  dFIdU;
 
  if (SolnBlk.Flow_Type != FLOWTYPE_INVISCID) {
    //Viscous
    DenseMatrix dGVdU(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D);
    DenseMatrix dGVdW(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D);
    dGVdW.zero();
    dGVdU.zero();
    DenseMatrix dWdQ(NUM_VAR_CHEM2D,NUM_VAR_CHEM2D); 
    dWdQ.zero();
    //transformation Jacobian
    SolnBlk.W[ii][jj].dWdU(dWdQ);
    //Inviscid + Viscous + Source (axisymmetric + turbulence + chemistry)     
    
    dGVdW = dGVdW_Viscous(SolnBlk,Input_Parameters, ii, jj);
    dGVdU = dGVdW *dWdQ;
    PIBJ  += dGVdU;
 }
  //  Source Jacobians (axisymmetric, turbulence and source)

  DenseMatrix SIBJ(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D);
  SIBJ.zero();
  SIBJ = SemiImplicitBlockJacobi(SolnBlk,ii,jj);
  PIBJ += SIBJ ;

  return PIBJ;
}
//SemiImplicit Block Jacobi
DenseMatrix SemiImplicitBlockJacobi(Chem2D_Quad_Block &SolnBlk,
				    const int &ii, const int &jj){
  
  int NUM_VAR_CHEM2D =  SolnBlk.NumVar() -1; 
  DenseMatrix SIBJ(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D);
  SIBJ.zero();

  //transformation matrix
  DenseMatrix dWdQ(NUM_VAR_CHEM2D,NUM_VAR_CHEM2D);  
  dWdQ.zero();
  SolnBlk.W[ii][jj].dWdU(dWdQ);
 
  //Add Jacobian for inviscid axisymmetric source terms
  if (SolnBlk.Axisymmetric) {
      SolnBlk.W[ii][jj].dSa_idU(SIBJ,SolnBlk.Grid.Cell[ii][jj].Xc,SolnBlk.Axisymmetric);
      // Add Jacobian for viscous axisymmetric source terms 
//     if(SolnBlk.Flow_Type != FLOWTYPE_INVISCID){
//        DenseMatrix dSa_VdU(NUM_VAR_CHEM2D,NUM_VAR_CHEM2D);
//        dSa_VdU.zero();
//        SolnBlk.W[ii][jj].dSa_vdU(dSa_VdU, dWdQ, SolnBlk.dWdx[ii][jj],SolnBlk.dWdy[ii][jj],
// 	 		      SolnBlk.Grid.Cell[ii][jj].Xc,
// 		 	      SolnBlk.Flow_Type, SolnBlk.Axisymmetric, SolnBlk.d_dWdx_dW[ii][jj][0], 
// 			      SolnBlk.d_dWdy_dW[ii][jj][0]);
//        SIBJ += dSa_VdU;
//     }
  }

  //Xinfeng: turbulence source Jacobian has a temperature outofrange problem
  //so commented out for debugging  ... temperarily
  //turbulence
  if((SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ) ||
     (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_EPSILON )){
     
    DenseMatrix dStdW(NUM_VAR_CHEM2D,NUM_VAR_CHEM2D); 
    dStdW.zero();
    DenseMatrix dStdU(NUM_VAR_CHEM2D,NUM_VAR_CHEM2D); 
    dStdU.zero();
    dStdW_Turbulence(dStdW, SolnBlk, ii, jj);
    dStdU = dStdW *dWdQ; 
    SIBJ +=dStdU;
  }
  //Add Jacobian for finite-rate chemistry source terms
  if (SolnBlk.W[ii][jj].React.reactset_flag != NO_REACTIONS){
    DenseMatrix dSwdU(NUM_VAR_CHEM2D,NUM_VAR_CHEM2D); 
    dSwdU.zero(); 
    SolnBlk.W[ii][jj].dSwdU(dSwdU);
    SIBJ += dSwdU;
  }  
  //Add Jacobian for gravitational source terms
  if (SolnBlk.Gravity){
    DenseMatrix dSgdU(NUM_VAR_CHEM2D,NUM_VAR_CHEM2D); 
    dSgdU.zero();
    SolnBlk.W[ii][jj].dSgdU(dSgdU);
    SIBJ += dSgdU;
  }
  
  return SIBJ;
  
}
/********************************************************
 * Routine: PointImplicitBlkJ  Inviscid Flux Jacobian   *
 *                                                      *
 * This routine returns the inviscid components of       *    
 * Point Implicit Block Jacobian matrix for the         *
 * specified local solution block.                      *
 *                                                      *
 ********************************************************/
// Based on HLLE || ROE Flux Function
DenseMatrix dFIdU_Inviscid(Chem2D_Quad_Block &SolnBlk, Chem2D_Input_Parameters &Input_Parameters,const int &ii, const int &jj){
  
  int NUM_VAR_CHEM2D =  SolnBlk.NumVar() -1; 
  DenseMatrix dFIdU(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D);
  dFIdU.zero();

   if (Input_Parameters.i_Flux_Function == FLUX_FUNCTION_HLLE ){
    dFIdU = dFIdU_Inviscid_HLLE(SolnBlk,ii, jj);
    
  }
  if (Input_Parameters.i_Flux_Function == FLUX_FUNCTION_ROE ){
    dFIdU = dFIdU_Inviscid_ROE(SolnBlk,ii, jj);
    
  }
  return dFIdU;
  
}
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
 DenseMatrix Rotation_Matrix2(Vector2D nface, int Size,  int A_matrix) 
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

} /* End of Rotation_Matrix. */

//Inviscid Flux Jacobian based on HLLE Flux Function
DenseMatrix dFIdU_Inviscid_HLLE(Chem2D_Quad_Block &SolnBlk, const int &ii, const int &jj){

  int NUM_VAR_CHEM2D =  SolnBlk.NumVar()-1; 
  DenseMatrix dFidU(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D),
              II(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D), 
              A_N(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D), 
              AI_N(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D),
              A_S(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D), 
              AI_S(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D),
              A_E(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D), 
              AI_E(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D),
              A_W(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D), 
              AI_W(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D),
              dFidU_N(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D),
              dFidU_S(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D),
              dFidU_E(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D),
	      dFidU_W(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D);

  Vector2D lambdas_N, lambdas_S, lambdas_E, lambdas_W;
  Vector2D nface_N, nface_S, nface_E, nface_W;
  
  double alpha_N, gamma_N;
  double alpha_S, gamma_S;
  double alpha_W, gamma_W;
  double alpha_E, gamma_E;

  if (ii < SolnBlk.ICl || ii > SolnBlk.ICu ||
      jj < SolnBlk.JCl || jj > SolnBlk.JCu) {
    // GHOST CELL
    dFidU.zero();
    return dFidU;
    
  } else {
//     // NON-GHOST CELL.
  
    dFidU.zero();
    II.identity();

    nface_N = SolnBlk.Grid.nfaceN(ii, jj);
    nface_S = SolnBlk.Grid.nfaceS(ii, jj);
    nface_E = SolnBlk.Grid.nfaceE(ii, jj);
    nface_W = SolnBlk.Grid.nfaceW(ii, jj);

     
    lambdas_N = HLLE_wavespeeds(W(SolnBlk.Uo[ii][jj]), 
				W(SolnBlk.Uo[ii][jj+1]), 
				nface_N);
    lambdas_S = HLLE_wavespeeds(W(SolnBlk.Uo[ii][jj]),
				W(SolnBlk.Uo[ii][jj-1]), 
				nface_S);
    lambdas_E = HLLE_wavespeeds(W(SolnBlk.Uo[ii][jj]), 
				W(SolnBlk.Uo[ii+1][jj]), 
				nface_E);
    lambdas_W = HLLE_wavespeeds(W(SolnBlk.Uo[ii][jj]),
				W(SolnBlk.Uo[ii-1][jj]),
				nface_W);
    
    dFidU_N.zero();
    dFIdU(dFidU_N, Rotate(W(SolnBlk.Uo[ii][jj]), nface_N));
    A_N  = Rotation_Matrix2(nface_N,NUM_VAR_CHEM2D, 1);
    AI_N = Rotation_Matrix2(nface_N,NUM_VAR_CHEM2D, 0);
    if (lambdas_N.x >= ZERO) {
   
       dFidU = (-SolnBlk.Grid.lfaceN(ii, jj)/SolnBlk.Grid.area(ii, jj))*AI_N*dFidU_N*A_N;
    } else if (lambdas_N.y <= ZERO) {
       // Do nothing.
    } else {
       alpha_N = lambdas_N.y/(lambdas_N.y-lambdas_N.x);
       gamma_N = (lambdas_N.x*lambdas_N.y)/(lambdas_N.y-lambdas_N.x);
     
       dFidU = (-SolnBlk.Grid.lfaceN(ii, jj)/SolnBlk.Grid.area(ii, jj))*AI_N*(alpha_N*dFidU_N-gamma_N*II)*A_N;
    
    } /* endif */

    dFidU_S.zero();
    dFIdU(dFidU_S, Rotate(W(SolnBlk.Uo[ii][jj]), nface_S));
    A_S  = Rotation_Matrix2(nface_S,NUM_VAR_CHEM2D, 1);
    AI_S = Rotation_Matrix2(nface_S,NUM_VAR_CHEM2D, 0);
    if (lambdas_S.x >= ZERO) {
       dFidU += (-SolnBlk.Grid.lfaceS(ii, jj)/SolnBlk.Grid.area(ii, jj))*AI_S*dFidU_S*A_S;
    } else if (lambdas_S.y <= ZERO) {
       // Do nothing.
    } else {
      alpha_S = lambdas_S.y/(lambdas_S.y-lambdas_S.x);
      gamma_S = (lambdas_S.x*lambdas_S.y)/(lambdas_S.y-lambdas_S.x);
      dFidU += (-SolnBlk.Grid.lfaceS(ii, jj)/SolnBlk.Grid.area(ii, jj))*AI_S*(alpha_S*dFidU_S-gamma_S*II)*A_S;
      
    } /* endif */

    dFidU_E.zero();
    dFIdU(dFidU_E, Rotate(W(SolnBlk.Uo[ii][jj]), nface_E));
    A_E  = Rotation_Matrix2(nface_E,NUM_VAR_CHEM2D, 1);
    AI_E = Rotation_Matrix2(nface_E,NUM_VAR_CHEM2D, 0);
    if (lambdas_E.x >= ZERO) {
       dFidU += (-SolnBlk.Grid.lfaceE(ii, jj)/SolnBlk.Grid.area(ii, jj))*AI_E*dFidU_E*A_E;
    } else if (lambdas_E.y <= ZERO) {
       // Do nothing.
    } else {
       alpha_E = lambdas_E.y/(lambdas_E.y-lambdas_E.x);
       gamma_E = (lambdas_E.x*lambdas_E.y)/(lambdas_E.y-lambdas_E.x);
       dFidU += (-SolnBlk.Grid.lfaceE(ii, jj)/SolnBlk.Grid.area(ii, jj))*AI_E*(alpha_E*dFidU_E-gamma_E*II)*A_E;
    } /* endif */

    dFidU_W.zero();
    dFIdU(dFidU_W, Rotate(W(SolnBlk.Uo[ii][jj]), nface_W));
    A_W  = Rotation_Matrix2(nface_W,NUM_VAR_CHEM2D, 1);
    AI_W = Rotation_Matrix2(nface_W,NUM_VAR_CHEM2D, 0);
    if (lambdas_W.x >= ZERO) {
       dFidU += (-SolnBlk.Grid.lfaceW(ii, jj)/SolnBlk.Grid.area(ii, jj))*AI_W*dFidU_W*A_W;
    } else if (lambdas_W.y <= ZERO) {
       // Do nothing.
    } else {
       alpha_W = lambdas_W.y/(lambdas_W.y-lambdas_W.x);
       gamma_W = (lambdas_W.x*lambdas_W.y)/(lambdas_W.y-lambdas_W.x);
       dFidU += (-SolnBlk.Grid.lfaceW(ii, jj)/SolnBlk.Grid.area(ii, jj))*AI_W*(alpha_W*dFidU_W-gamma_W*II)*A_W;
    } /* endif */
   }// endif

    
  return (dFidU);

} /*********/
//Inviscid Flux Jacobian based on ROE Flux Function
DenseMatrix dFIdU_Inviscid_ROE(Chem2D_Quad_Block &SolnBlk, const int &ii, const int &jj){
  
  int NUM_VAR_CHEM2D =  SolnBlk.NumVar()-1 ; 
  DenseMatrix dFidU(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D),
              II(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D), 
              A_N(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D), 
              AI_N(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D),
              A_S(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D), 
              AI_S(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D),
              A_E(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D), 
              AI_E(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D),
              A_W(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D), 
              AI_W(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D),
              dFidU_N(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D),
              dFidU_S(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D),
              dFidU_E(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D),
	      dFidU_W(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D);

 
  Vector2D nface_N, nface_S, nface_E, nface_W;

  Chem2D_pState Wa_N, Wa_S, Wa_E, Wa_W;
  Chem2D_pState lambdas_N, lambdas_S, lambdas_E, lambdas_W;
  Chem2D_pState Left, Right;

  double alpha_N, gamma_N;
  double alpha_S, gamma_S;
  double alpha_W, gamma_W;
  double alpha_E, gamma_E;

  dFidU.zero();
  if (ii < SolnBlk.ICl || ii > SolnBlk.ICu ||
      jj < SolnBlk.JCl || jj > SolnBlk.JCu) {
    // GHOST CELL
  
    return dFidU;
    
  } else {
//     // NON-GHOST CELL.
    
 
    II.identity();

    nface_N = SolnBlk.Grid.nfaceN(ii, jj);
    nface_S = SolnBlk.Grid.nfaceS(ii, jj);
    nface_E = SolnBlk.Grid.nfaceE(ii, jj);
    nface_W = SolnBlk.Grid.nfaceW(ii, jj);
  
    Wa_N = RoeAverage(W(SolnBlk.Uo[ii][jj]), W(SolnBlk.Uo[ii][jj+1]));
    Wa_S = RoeAverage(W(SolnBlk.Uo[ii][jj]), W(SolnBlk.Uo[ii][jj-1]));
    Wa_E = RoeAverage(W(SolnBlk.Uo[ii][jj]), W(SolnBlk.Uo[ii+1][jj]));
    Wa_W = RoeAverage(W(SolnBlk.Uo[ii][jj]), W(SolnBlk.Uo[ii-1][jj]));
  
    lambdas_N = Wa_N.lambda_x();  
    lambdas_S = Wa_S.lambda_x(); 
    lambdas_E = Wa_E.lambda_x(); 
    lambdas_W = Wa_W.lambda_x(); 
    
    dFidU_N.zero();
    dFIdU(dFidU_N, Rotate(W(SolnBlk.Uo[ii][jj]), nface_N));
    A_N  = Rotation_Matrix2(nface_N,NUM_VAR_CHEM2D, 1);
    AI_N = Rotation_Matrix2(nface_N,NUM_VAR_CHEM2D, 0);
  
    for(int index=1; index<=NUM_VAR_CHEM2D; index++)
      {
	II(index,index) =  lambdas_N[index];
      }
    // ROE FLUX FUNCTION
    dFidU =HALF*(-SolnBlk.Grid.lfaceN(ii, jj)/SolnBlk.Grid.area(ii, jj))*AI_N*(dFidU_N + II)*A_N;

    dFidU_S.zero();
    dFIdU(dFidU_S, Rotate(W(SolnBlk.Uo[ii][jj]), nface_S));
    A_S  = Rotation_Matrix2(nface_S,NUM_VAR_CHEM2D, 1);
    AI_S = Rotation_Matrix2(nface_S,NUM_VAR_CHEM2D, 0);
  
    dFidU += HALF*(-SolnBlk.Grid.lfaceS(ii, jj)/SolnBlk.Grid.area(ii, jj))*AI_S*(dFidU_S + II)*A_S;
    
    for(int index=1; index<=NUM_VAR_CHEM2D; index++)
      {
	II(index,index) =  lambdas_S[index];
      }
    
    dFidU_E.zero();
    dFIdU(dFidU_E, Rotate(W(SolnBlk.Uo[ii][jj]), nface_E));
    A_E  = Rotation_Matrix2(nface_E,NUM_VAR_CHEM2D, 1);
    AI_E = Rotation_Matrix2(nface_E,NUM_VAR_CHEM2D, 0);
   
    for(int index=1; index<=NUM_VAR_CHEM2D; index++)
      {
	II(index,index) =  lambdas_E[index];
      }
    
    dFidU += HALF*(-SolnBlk.Grid.lfaceE(ii, jj)/SolnBlk.Grid.area(ii, jj))*AI_E*(dFidU_E+ II )*A_E;
    dFidU_W.zero();
    dFIdU(dFidU_W, Rotate(W(SolnBlk.Uo[ii][jj]), nface_W));
    A_W  = Rotation_Matrix2(nface_W,NUM_VAR_CHEM2D, 1);
    AI_W = Rotation_Matrix2(nface_W,NUM_VAR_CHEM2D, 0);
  
    for(int index=1; index<=NUM_VAR_CHEM2D; index++)
      {
	II(index,index) =  lambdas_W[index];
      }
    
    dFidU += HALF*(-SolnBlk.Grid.lfaceW(ii, jj)/SolnBlk.Grid.area(ii, jj))*AI_W*(dFidU_W + II  )*A_W;
   
  }// endif

    
  return (dFidU);

} /*********/
/*Viscous Flux Jacobians  */
/*Viscous Flux Jacobians  */
void dFvdW_Laminar(DenseMatrix &dFvdW,  Chem2D_Quad_Block &SolnBlk,const int &ii, const int &jj){

  //planar flow

  double  d_dWdx_dW=0;
  double  d_dWdy_dW=0;

  d_dWdx_dW = SolnBlk.d_dWdx_dW[ii][jj][0];
  d_dWdy_dW = SolnBlk.d_dWdy_dW[ii][jj][0];
    
  double t1,t2,t5,t8,t9,t11,t13, t16,t18;
  double t20,t21, t22, t27, t35,t38,t46;
    
  double kappa, Cp, mu;
  double rho, U, V;
  double Temp =SolnBlk.W[ii][jj].T();
    
  kappa =SolnBlk.W[ii][jj].kappa();
  Cp =SolnBlk.W[ii][jj].Cp();
  mu  =SolnBlk.W[ii][jj].mu();
    
  rho = SolnBlk.W[ii][jj].rho;
  U =  SolnBlk.W[ii][jj].v.x;
  V =  SolnBlk.W[ii][jj].v.y;

  t1 = mu;
  t2 = d_dWdx_dW;
  t5 = d_dWdy_dW;
  t8 = SolnBlk.W[ii][jj].dmudT();
  t9 = SolnBlk.dWdx[ii][jj].v.x;
  t11 = SolnBlk.dWdy[ii][jj].v.y;
  t13 = 2.0/3.0*t9-t11/3.0;
  t16 = d_dWdy_dW;
  t18 = d_dWdx_dW;
  t20 = SolnBlk.dWdy[ii][jj].v.x;
  t21 = SolnBlk.dWdx[ii][jj].v.y;
  t22 = t20+t21;
  t35 = U*t1;
  t38 = V*t1;
  t46 = d_dWdx_dW;
    
  dFvdW(1,1) = 4.0/3.0*t1*t2;
  dFvdW(1,2) = -2.0/3.0*t1*t5;
  dFvdW(1,3) = 2.0*t8*t13;
    
  dFvdW(2,1) = t1*t16;
  dFvdW(2,2) = t1*t18;
  dFvdW(2,3) = t8*t22;

  double Sum_q = 0.0;
  double Sum_dq = 0.0;
  int ns = SolnBlk.W[0][0].ns;
    
  for(int Num = 0; Num<ns; Num++)
    {
      //for each species	// h*Dm*gradc  
      Sum_q +=  (SolnBlk.W[ii][jj].specdata[Num].Enthalpy(Temp)+SolnBlk.W[ii][jj].specdata[Num].Heatofform())*(SolnBlk.W[ii][jj].spec[Num].diffusion_coef)*SolnBlk.dWdx[ii][jj].spec[Num].c;
      //dhdT *(Dm)*gradc
      Sum_dq +=  SolnBlk.W[ii][jj].specdata[Num].Enthalpy_prime(Temp)*(SolnBlk.W[ii][jj].spec[Num].diffusion_coef)*SolnBlk.dWdx[ii][jj].spec[Num].c;
      
    }  
  dFvdW(3,0) = Sum_q ;
  dFvdW(3,1) = 2.0*t1*t13+4.0/3.0*t35*t2+t38*t16;
  dFvdW(3,2) = -2.0/3.0*t35*t5+t1*t22+t38*t18;
  dFvdW(3,3) = kappa*t46+rho*Sum_dq+2.0*U*t8*t13+V*t8*t22;
      
  //multispecies
  int NUM_VAR = NUM_CHEM2D_VAR_SANS_SPECIES; 
  for(int Num = 0; Num<(ns-1); Num++)
    {//h*rho*(Dm)*d_dWdx_dW
      dFvdW(3, NUM_VAR+Num) = (SolnBlk.W[ii][jj].specdata[Num].Enthalpy(Temp)+SolnBlk.W[ii][jj].specdata[Num].Heatofform())*rho*(SolnBlk.U[ii][jj].rhospec[Num].diffusion_coef)*d_dWdx_dW;
      //(Dm)*gradc
      dFvdW(NUM_VAR+Num, 0) =(SolnBlk.U[ii][jj].rhospec[Num].diffusion_coef)*SolnBlk.dWdx[ii][jj].spec[Num].c;
      //rho*(Dm)*dcxdc
      dFvdW(NUM_VAR+Num,NUM_VAR+Num) = rho*(SolnBlk.U[ii][jj].rhospec[Num].diffusion_coef)*d_dWdx_dW;
    }

  //The following entries are different for axisymmetric flow and planar flow
  // For those terms involves fluid molecular stess tau_xx and tau_yy
  // They are some terms in the scond row  and fourth of Jacobian matrix
  if(SolnBlk.Axisymmetric){ 
    double r;
    if (SolnBlk.Axisymmetric == 1) {
      r = SolnBlk.Grid.Cell[ii][jj].Xc.y;
      
      t13 -= (V/r)/3.0;
    
      dFvdW(1,2) = 2.0/3.0*t1*(-t5-ONE/r);
      dFvdW(1,3) = 2.0*t8*t13;
      dFvdW(3,2) = 2.0/3.0*t35*(-t5-ONE/r)+t1*t22+t38*t18;
      dFvdW(3,3) = kappa*t46+rho*Sum_dq +2.0*U*t8*t13+V*t8*t22;
    }
    if (SolnBlk.Axisymmetric == 2) {
      r = SolnBlk.Grid.Cell[ii][jj].Xc.x;
      double term1 = 2.0/3.0*t2 - 1.0/r*1.0/3.0;
      
      t13 -= (U/r)/3.0;

      dFvdW(1,1) = 4.0/3.0*t1*term1;
      dFvdW(1,3) = 2.0*t8*t13;
      dFvdW(3,1) = 2.0*t1*t13+2.0*t35*term1+t38*t16;
      dFvdW(3,3) = kappa*t46+rho*Sum_dq +2.0*U*t8*t13+V*t8*t22;
    }
  }//end of axisymmetric case


}//Laminar viscous Jacobian (X)
 
void dGvdW_Laminar(DenseMatrix &dGvdW,  Chem2D_Quad_Block &SolnBlk,const int &ii, const int &jj){

  //planar flow

  double  d_dWdx_dW=0;
  double  d_dWdy_dW=0;

  d_dWdx_dW = SolnBlk.d_dWdx_dW[ii][jj][0];
  d_dWdy_dW = SolnBlk.d_dWdy_dW[ii][jj][0];
    
  double t1,t2,t5,t8,t9,t11,t13, t16,t18;
  double t20,t21, t22, t27, t35,t38,t46;
    
  double kappa, Cp, mu;
  double rho, U, V;
  double Temp =SolnBlk.W[ii][jj].T();
    
  kappa =SolnBlk.W[ii][jj].kappa();
  Cp =SolnBlk.W[ii][jj].Cp();
  mu  =SolnBlk.W[ii][jj].mu();
    
  rho = SolnBlk.W[ii][jj].rho;
  U =  SolnBlk.W[ii][jj].v.x;
  V =  SolnBlk.W[ii][jj].v.y;

  t1 = mu;
  t2 = d_dWdy_dW;
  t5 = d_dWdx_dW;
  t8 = SolnBlk.W[ii][jj].dmudT();
  t9 = SolnBlk.dWdy[ii][jj].v.x;
  t11 = SolnBlk.dWdx[ii][jj].v.y;
  t13 = 2.0/3.0*t9-t11/3.0;
  t16 = d_dWdx_dW;
  t18 = d_dWdy_dW;
  t20 = SolnBlk.dWdy[ii][jj].v.y;
  t21 = SolnBlk.dWdx[ii][jj].v.x;
  t22 = t20+t21;
  t35 = U*t1;
  t38 = V*t1;
  t46 = d_dWdy_dW;
    
  dGvdW(1,1) = t1*t2;
  dGvdW(1,2) = -2.0/3.0*t1*t5;
  dGvdW(1,3) = t8*t22;
    
  dGvdW(2,1) = -2.0/3.0*t1*t16;
  dGvdW(2,2) = 4.0/3.0*t1*t18;
  dGvdW(2,3) = 2.0*t8*t13;

  double Sum_q = 0.0;
  double Sum_dq = 0.0;
  int ns = SolnBlk.W[0][0].ns;
    
  for(int Num = 0; Num<ns; Num++)
    {
      //for each species	// h*Dm*gradc  
      Sum_q +=  (SolnBlk.W[ii][jj].specdata[Num].Enthalpy(Temp)+SolnBlk.W[ii][jj].specdata[Num].Heatofform())*(SolnBlk.W[ii][jj].spec[Num].diffusion_coef)*SolnBlk.dWdy[ii][jj].spec[Num].c;
      //dhdT *(Dm)*gradc
      Sum_dq +=  SolnBlk.W[ii][jj].specdata[Num].Enthalpy_prime(Temp)*(SolnBlk.W[ii][jj].spec[Num].diffusion_coef)*SolnBlk.dWdy[ii][jj].spec[Num].c;
      
    }  
  dGvdW(3,0) = Sum_q ;
  dGvdW(3,1) = t1*t22+t35*t2-2.0/3.0*t38*t16;
  dGvdW(3,2) = t35*t5+2.0*t1*t13+4.0/3.0*t38*t18;
  dGvdW(3,3) = kappa*t46+rho*Sum_dq+U*t8*t22+2.0*V*t8*t13;
      
  //multispecies
  int NUM_VAR =NUM_CHEM2D_VAR_SANS_SPECIES; 
  for(int Num = 0; Num<(ns-1); Num++)
    {//h*rho*(Dm)*d_dWdy_dW
      dGvdW(3, NUM_VAR+Num) = (SolnBlk.W[ii][jj].specdata[Num].Enthalpy(Temp)+SolnBlk.W[ii][jj].specdata[Num].Heatofform())*rho*(SolnBlk.U[ii][jj].rhospec[Num].diffusion_coef)*d_dWdy_dW;
      //(Dm)*gradc
      dGvdW(NUM_VAR+Num, 0) =(SolnBlk.U[ii][jj].rhospec[Num].diffusion_coef)*SolnBlk.dWdy[ii][jj].spec[Num].c;
      //rho*(Dm)*dcydc
      dGvdW(NUM_VAR+Num,NUM_VAR+Num) = rho*(SolnBlk.U[ii][jj].rhospec[Num].diffusion_coef)*d_dWdy_dW;
    }

  //The following entries are different for axisymmetric flow and planar flow
  // For those terms involves fluid molecular stess tau_xx and tau_yy
  // They are some terms in the third row  and fourth of Jacobian matrix
  if(SolnBlk.Axisymmetric){ 
    double r;
    if (SolnBlk.Axisymmetric == 1) {
      r = SolnBlk.Grid.Cell[ii][jj].Xc.y;
      
      t13 -= (V/r)/3.0;
      t18 = 2.0/3.0*t18 -1.0/r*1.0/3.0;
    
      dGvdW(2,2) = 2.0*t1*t18;
      dGvdW(2,3) = 2.0*t8*t13;
      dGvdW(3,2) = t35*t5+2.0*t1*t13+2.0*t38*t18;
      dGvdW(3,3) = kappa*t46+rho*Sum_dq +2.0*U*t8*t22+2.0*V*t8*t13;
    }
    if (SolnBlk.Axisymmetric == 2) {
      r = SolnBlk.Grid.Cell[ii][jj].Xc.x;
      double term1 = -t16 - 1.0/r;
      
      t13 -= (U/r)/3.0;

      dGvdW(2,1) = -2.0/3.0*t1*term1;
      dGvdW(2,3) = 2.0*t8*t13;
      dGvdW(3,1) = t1*t22+t35*t2+2.0/3.0*t38*t13;
      dGvdW(3,3) = kappa*t46+rho*Sum_dq +U*t8*t22+2.0*V*t8*t13;
    }
  }//end of axisymmetric case

 

}//Laminar viscous Jacobian (Y)

void dFvdW_Laminar_Neigbour(DenseMatrix &dFvdW,Chem2D_Quad_Block &SolnBlk, 
				   const string &Orient, const int &ii, const int &jj){
  double d_dWdx_dW, d_dWdy_dW;
  int i, j;
  i = ii;
  j=jj;
  
  if( Orient == "NORTH"){
    j = jj+1;
    d_dWdx_dW = SolnBlk.d_dWdx_dW[ii][jj][1];
    d_dWdy_dW = SolnBlk.d_dWdy_dW[ii][jj][1];
      
  }else if( Orient == "SOUTH"){
    j = jj-1;
    d_dWdx_dW = SolnBlk.d_dWdx_dW[ii][jj][2];
    d_dWdy_dW = SolnBlk.d_dWdy_dW[ii][jj][2];
    
  }else if( Orient == "WEST"){
    i = ii-1;
    d_dWdx_dW = SolnBlk.d_dWdx_dW[ii][jj][3];
    d_dWdy_dW = SolnBlk.d_dWdy_dW[ii][jj][3];
    
  }else if( Orient == "EAST"){
    i = ii+1;
    d_dWdx_dW = SolnBlk.d_dWdx_dW[ii][jj][4];
    d_dWdy_dW = SolnBlk.d_dWdy_dW[ii][jj][4];
    
  }
 if (i < SolnBlk.ICl || i > SolnBlk.ICu ||
      j < SolnBlk.JCl || j > SolnBlk.JCu) {
    // GHOST CELL
  } else {
    //     // NON-GHOST CELL.
  double t1,t2,t5,t12, t15;
  double kappa = SolnBlk.W[i][j].kappa();
  double Temp =  SolnBlk.W[i][j].T();

  double rho,U,V;
  rho = SolnBlk.W[ii][jj].rho;
  U =  SolnBlk.W[ii][jj].v.x;
  V =  SolnBlk.W[ii][jj].v.y;
  
  t1 = SolnBlk.W[i][j].mu();
  t2 = d_dWdx_dW;
  t5 = d_dWdy_dW;
  t12 = U*t1;
  t15 = V*t1;
  
  dFvdW(1,1) = 4.0/3.0*t1*t2;
  dFvdW(1,2) = -2.0/3.0*t1*t5;
  dFvdW(2,1) = t1*t5;
  dFvdW(2,2) = t1*t2;
  
  dFvdW(3,1) = 4.0/3.0*t12*t2+t15*t5;
  dFvdW(3,2) = t15*t2-2.0/3.0*t12*t5;
  dFvdW(3,3) = kappa*d_dWdx_dW;
  
  int ns = SolnBlk.W[0][0].ns;
  int NUM_VAR =NUM_CHEM2D_VAR_SANS_SPECIES; 
  for(int Num = 0; Num<(ns-1); Num++)
    {
      dFvdW(3,NUM_VAR+Num) = (SolnBlk.W[ii][jj].specdata[Num].Enthalpy(Temp)+SolnBlk.W[ii][jj].specdata[Num].Heatofform())*rho*(SolnBlk.U[ii][jj].rhospec[Num].diffusion_coef)*d_dWdx_dW;
      dFvdW(NUM_VAR+Num,NUM_VAR+Num) = rho*(SolnBlk.U[ii][jj].rhospec[Num].diffusion_coef)*d_dWdx_dW;
      
    }
  }

}//end of neigbour cell for x viscous flux Jacobian
void dGvdW_Laminar_Neigbour(DenseMatrix &dGvdW,Chem2D_Quad_Block &SolnBlk, 
				   const string &Orient, const int &ii, const int &jj){
 double d_dWdx_dW, d_dWdy_dW;
  int i, j;
  i = ii;
  j=jj;

  if( Orient == "NORTH"){
    j = jj+1;
    d_dWdx_dW = SolnBlk.d_dWdx_dW[ii][jj][1];
    d_dWdy_dW = SolnBlk.d_dWdy_dW[ii][jj][1];

  }else if( Orient == "SOUTH"){
    j = jj-1;
    d_dWdx_dW = SolnBlk.d_dWdx_dW[ii][jj][2];
    d_dWdy_dW = SolnBlk.d_dWdy_dW[ii][jj][2];

  }else if( Orient == "WEST"){
    i = ii-1;
    d_dWdx_dW = SolnBlk.d_dWdx_dW[ii][jj][3];
    d_dWdy_dW = SolnBlk.d_dWdy_dW[ii][jj][3];

  }else if( Orient == "EAST"){
    i = ii+1;
    d_dWdx_dW = SolnBlk.d_dWdx_dW[ii][jj][4];
    d_dWdy_dW = SolnBlk.d_dWdy_dW[ii][jj][4];

  }
  if (i < SolnBlk.ICl || i > SolnBlk.ICu ||
      j < SolnBlk.JCl || j > SolnBlk.JCu) {
    // GHOST CELL
  } else {
    //     // NON-GHOST CELL.
    double t1,t2,t4, t12, t15;
    double kappa = SolnBlk.W[i][j].kappa();
    double Temp =  SolnBlk.W[i][j].T();
    
    double rho,U,V;
    rho = SolnBlk.W[ii][jj].rho;
    U =  SolnBlk.W[ii][jj].v.x;
    V =  SolnBlk.W[ii][jj].v.y;
    
    t1 = SolnBlk.W[i][j].mu();
    t2 = d_dWdy_dW;
    t4 = d_dWdx_dW;
    t12 = V*t1;
    t15 = U*t1;
   
    dGvdW(1, 1) = t1*t2;
    dGvdW(1, 2) = t1*t4;
    dGvdW(2, 1) = -2.0/3.0*t1*t4;
    dGvdW(2, 2) = 4.0/3.0*t1*t2;
    dGvdW(3, 1) = -2.0/3.0*t12*t4+t15*t2;
    dGvdW(3, 2) = t15*t4+4.0/3.0*t12*t2;
    dGvdW(3, 3) = kappa*d_dWdy_dW;
    
    int ns = SolnBlk.W[0][0].ns;
    int NUM_VAR = NUM_CHEM2D_VAR_SANS_SPECIES; 
    for(int Num = 0; Num<(ns-1); Num++)
      {
	dGvdW(3,NUM_VAR+Num) = (SolnBlk.W[ii][jj].specdata[Num].Enthalpy(Temp)+SolnBlk.W[ii][jj].specdata[Num].Heatofform())*rho*(SolnBlk.U[ii][jj].rhospec[Num].diffusion_coef)*d_dWdy_dW;
	dGvdW(NUM_VAR+Num,NUM_VAR+Num) = rho*(SolnBlk.U[ii][jj].rhospec[Num].diffusion_coef)*d_dWdy_dW;
	
      }
  }
 
}//end of neigbour cell for viscous Jacobian


//Turbulent Flux Jacobian
void dFvdW_Turbulent(DenseMatrix &dFvdW, Chem2D_Quad_Block &SolnBlk,
			    const int &ii, const int &jj){

  double  d_dWdx_dW, d_dWdy_dW;
  d_dWdx_dW = SolnBlk.d_dWdx_dW[ii][jj][0];
  d_dWdy_dW = SolnBlk.d_dWdy_dW[ii][jj][0];

  double t1,t2,t3,t5,t7,t11,t12,t13, t15,t18,t19,t23, t24, t27, t31,t32;
  double t33,t37,t38,t39,t41,t45,t46,t50, t56, t57;
  double t85,t86, t87, t90, t110, t140;
  double t149, t150, t153, t154, t156, t160, t164;
  double t194, t195, t198, t204;
  
  double kappa, Cp, mu, sigma, sigma_star;
  double rho, U, V, k, omega;
  double mu_turb, kappa_turb, Pr_turb, Dm_turb;
  double Temp =SolnBlk.W[ii][jj].T();
  double Rmix = SolnBlk.W[ii][jj].Rtot();  
   
  kappa =SolnBlk.W[ii][jj].kappa();
  Cp =SolnBlk.W[ii][jj].Cp();
  mu  =SolnBlk.W[ii][jj].mu();
    
  mu_turb = SolnBlk.W[ii][jj].eddy_viscosity();  
  Pr_turb = SolnBlk.W[ii][jj].Pr_turb(); 
  Dm_turb = SolnBlk.W[ii][jj].Dm_turb(); 
  
  sigma = SolnBlk.W[0][0].sigma;
  sigma_star = SolnBlk.W[0][0].sigma_star;  
  
  rho = SolnBlk.W[ii][jj].rho;
  U =  SolnBlk.W[ii][jj].v.x;
  V =  SolnBlk.W[ii][jj].v.y;
  k = SolnBlk.W[ii][jj].k;
  omega = SolnBlk.W[ii][jj].omega;
  
  t1 = 1/max(omega, TOLER);
  t2 = k*t1;
  t3 = SolnBlk.dWdx[ii][jj].v.x; 
  t5 =  SolnBlk.dWdy[ii][jj].v.y; 
  t7 = 2.0/3.0*t3-t5/3.0;
  t11 = 2.0*t2*t7-2.0/3.0*k;
  t12 = mu;
  t13 =  d_dWdx_dW;
  t15 = rho*k;
  t18 = t12*t13+t15*t1*t13;
  t19 =  d_dWdy_dW;
  t23 = -t12*t19-t15*t1*t19;
  t24 = SolnBlk.W[ii][jj].dmudT();
  t27 = rho*t1;
  t31 = 2.0*t27*t7-2.0/3.0*rho;
  t32 = omega*omega;
  t33 = 1/t32;
  t37 = SolnBlk.dWdy[ii][jj].v.x;
  t38 = SolnBlk.dWdx[ii][jj].v.y;
  t39 = t37+t38;
  t41 = d_dWdx_dW;
  t45 = t12*t41+t15*t1*t41;
  t46 = d_dWdx_dW;
  t50 = t12*t46+t15*t1*t46;
  t56 = Cp/Pr_turb;
  //t57 = dTdx
  t57 = (ONE/(SolnBlk.W[ii][jj].rho*Rmix)) * (SolnBlk.dWdx[ii][jj].p - 
	(SolnBlk.W[ii][jj].p/SolnBlk.W[ii][jj].rho)*SolnBlk.dWdx[ii][jj].rho);
  t85 = SolnBlk.dWdx[ii][jj].k;
  t86 = t1*t85;
  t87 = sigma_star*k*t86;
  t90 = t1*t39;
  t110 = d_dWdx_dW;
  t140 = t24*t85;
  t149 = sigma_star*rho;
  t150 = t149*t86;
  t153 = d_dWdx_dW;
  t154 = (t12+t149*t2)*t153;
  t156 = V*rho;
  t160 = k*t33;
  t164 = t149*t160*t85;
  t194 =  SolnBlk.dWdx[ii][jj].omega;
  t195 = t1*t194;
  t198 = sigma*rho;
  t204 =  d_dWdx_dW;

  dFvdW(1,0) = t11;
  dFvdW(1,1) = 4.0/3.0*t18;
  dFvdW(1,2) = 2.0/3.0*t23;
  dFvdW(1,3) = 2.0*t24*t7;
  dFvdW(1,4) = t31;
  dFvdW(1,5) = -2.0*t15*t33*t7;
  dFvdW(2,0) = t2*t39;
  dFvdW(2,1) = t45;
  dFvdW(2,2) = t50;
  dFvdW(2,3) = t24*t39;
  dFvdW(2,4) = t27*t39;
  dFvdW(2,5) = -t15*t33*t39;
    
  double Sum_q = 0.0;
  double Sum_dq = 0.0;
  int ns = SolnBlk.W[0][0].ns;
  
  for(int Num = 0; Num<(ns-1); Num++)
    {
      //for each species	// h*Dm*gradc  
      Sum_q +=  (SolnBlk.W[ii][jj].specdata[Num].Enthalpy(Temp)+SolnBlk.W[ii][jj].specdata[Num].Heatofform())*(SolnBlk.W[ii][jj].spec[Num].diffusion_coef+Dm_turb)*SolnBlk.dWdx[ii][jj].spec[Num].c;
      //dhdT *(Dm+Dmt)*gradc
      Sum_dq +=  SolnBlk.W[ii][jj].specdata[Num].Enthalpy_prime(Temp)*(SolnBlk.W[ii][jj].spec[Num].diffusion_coef+Dm_turb)*SolnBlk.dWdx[ii][jj].spec[Num].c;
      
    }
    
  dFvdW(3,0) = t56*t2*t57+ Sum_q+t87+U*t11+V*k*t90;
  dFvdW(3,1) = 2.0*t12*t7+2.0*t15*t1*t7-2.0/3.0*t15+4.0/3.0*U*t18+V*t45;
  dFvdW(3,2) = 2.0/3.0*U*t23+t12*t39+t15*t90+V*t50;
  dFvdW(3,3) = (kappa+t56*t15*t1)*t110+ rho*Sum_dq +t140+2.0*U*t24*t7+V*t24*t39;
  dFvdW(3,4) = t56*t27*t57+t150+t154+U*t31+t156*t90;
  dFvdW(3,5) = -t56*rho*t160*t57-t164-2.0*U*rho*t160*t7-t156*t160*t39;
      
  int NUM_VAR = NUM_CHEM2D_VAR_SANS_SPECIES; 
  for(int Num = 0; Num<(ns-1); Num++)
    {//h*rho*(D+Dt)*d_dWdx_dW
      dFvdW(3, NUM_VAR+Num) = (SolnBlk.W[ii][jj].specdata[Num].Enthalpy(Temp)+SolnBlk.W[ii][jj].specdata[Num].Heatofform())*rho*(SolnBlk.U[ii][jj].rhospec[Num].diffusion_coef+Dm_turb)*d_dWdx_dW;
      //(D+Dt)*gradc
      dFvdW(NUM_VAR+Num, 0) =(SolnBlk.U[ii][jj].rhospec[Num].diffusion_coef+Dm_turb)*SolnBlk.dWdx[ii][jj].spec[Num].c;
      //rho*(D+Dt)*dcxdc
      dFvdW(NUM_VAR+Num,NUM_VAR+Num) = rho*(SolnBlk.U[ii][jj].rhospec[Num].diffusion_coef+Dm_turb)*d_dWdx_dW;
    }// end the center ones
  
  dFvdW(4,0) = t87;
  dFvdW(4,3) = t140;
  dFvdW(4,4) = t150+t154;
  dFvdW(4,5) = -t164;
  dFvdW(5,0) = sigma*k*t195;
  dFvdW(5,3) = t24*t194;
  dFvdW(5,4) = t198*t195;
  dFvdW(5,5) = -t198*t160*t194+(t12+t198*t2)*t204;
    
  //The following entries are different for axisymmetric flow and planar flow
  // For matrix row 1 and row 3
  if (SolnBlk.Axisymmetric == 1) {

    double term1,term2,term3,term5, term7,term10;
    double term14,term15,term16,term18;
    double term21,term22,term23,term27,term28,term31,term35,term36, term37, term41, term42, term43;
    double term45, term49,term50, term54, term60, term61, term89,term90, term91,term94, term114; 
    double term144, term153, term154, term157, term158, term160, term164,term168 , term198, term199, term202, term208; 
    
    double r = SolnBlk.Grid.Cell[ii][jj].Xc.y;
    term1 = 1/max(omega,TOLER);
    term2 = k*term1;
    term3 = SolnBlk.dWdx[ii][jj].v.x; 
    term5 =  SolnBlk.dWdy[ii][jj].v.y; 
    term7 = 1/r;
    term10 = 2.0/3.0*term3-term5/3.0-V*term7/3.0;
    term14 = 2.0*term2*term10-2.0/3.0*k;
    term15 = mu;
    term16 = d_dWdx_dW;
    term18 = rho*k;
    term21 = term15*term16+term18*term1*term16;
    term22 = d_dWdy_dW;
    term23 = -term22-term7;
    term27 = term15*term23/3.0+term18*term1*term23/3.0;
    term28 = SolnBlk.W[ii][jj].dmudT();
    term31 = rho*term1;
    term35 = 2.0*term31*term10-2.0/3.0*rho;
    term36 = omega*omega;
    term37 = 1/term36;
    term41 = SolnBlk.dWdy[ii][jj].v.x; 
    term42 = SolnBlk.dWdx[ii][jj].v.y;
    term43 = term41+term42;
    term45 = d_dWdy_dW;
    term49 = term15*term45+term18*term1*term45;
    term50 = d_dWdx_dW;
    term54 = term15*term50+term18*term1*term50;
    term60 = Cp/Pr_turb;
    //term61 = dTdx
    term61 = (ONE/(SolnBlk.W[ii][jj].rho*Rmix)) * (SolnBlk.dWdx[ii][jj].p - 
	(SolnBlk.W[ii][jj].p/SolnBlk.W[ii][jj].rho)*SolnBlk.dWdx[ii][jj].rho);
    term89 =  SolnBlk.dWdx[ii][jj].k;
    term90 = term1*term89;
    term91 = sigma*k*term90;
    term94 = term1*term43;
    term114 = d_dWdx_dW;
    term144 = term28*term89;
    term153 = sigma*rho;
    term154 = term153*term90;
    term157 = d_dWdx_dW;
    term158 = (term15+term153*term2)*term157;
    term160 = V*rho;
    term164 = k*term37;
    term168 = term153*term164*term89;
    term198 = SolnBlk.dWdx[ii][jj].omega;
    term199 = term1*term198;
    term202 = sigma_star*rho;
    term208 = d_dWdx_dW;
	
    dFvdW(1,0) = term14;
    dFvdW(1,1) = 4.0/3.0*term21;
    dFvdW(1,2) = 2.0*term27;
    dFvdW(1,3) = 2.0*term28*term10;
    dFvdW(1,4) = term35;
    dFvdW(1,5) = -2.0*term18*term37*term10;
    
    dFvdW(3,0) = term60*term2*term61+ Sum_q+term91+U*term14+V*k*term94;
    dFvdW(3,1) = 2.0*term15*term10+2.0*term18*term1*term10-2.0/3.0*term18+4.0/3.0*U*term21+V*term49;
    dFvdW(3,2) = 2.0*U*term27+term15*term43+term18*term94+V*term54;
    dFvdW(3,3) = (kappa+term60*term18*term1)*term114+ rho*Sum_dq+term144+2.0*U*term28*term10+V*term28*term43;
    dFvdW(3,4) = term60*term31*term61+term154+term158+U*term35+term160*term94;
    dFvdW(3,5) = -term60*rho*term164*term61-term168-2.0*U*rho*term164*term10-term160*term164*term43;
    
  }

  if (SolnBlk.Axisymmetric == 2) {
    double tt1,tt2,tt3,tt5,tt7,tt10,tt14,tt15,tt16;
    double tt19,tt21,tt24,tt25,tt29,tt30,tt33;
    double tt37,tt38,tt39,tt43,tt44,tt45,tt47,tt51, tt52;
    double tt56,tt62,tt63,tt79,tt80, tt81,tt84;
    double tt104,tt120,tt129,tt130,tt133,tt134, tt136,tt140;
    double tt144,tt152,tt155,tt157,tt160,tt164, tt165,tt168,tt174;
    
    double r = SolnBlk.Grid.Cell[ii][jj].Xc.x;
    
    tt1 = 1/max(omega,TOLER);
    tt2 = k*tt1;
    tt3 =  SolnBlk.dWdx[ii][jj].v.x; 
    tt5 =  SolnBlk.dWdy[ii][jj].v.y;
    tt7 = 1/r;
    tt10 = 2.0/3.0*tt3-tt5/3.0-U*tt7/3.0;
    tt14 = 2.0*tt2*tt10-2.0/3.0*k;
    tt15 =  SolnBlk.W[ii][jj].mu();
    tt16 = d_dWdx_dW;
    tt19 = 2.0/3.0*tt16-tt7/3.0;
    tt21 = rho*k;
    tt24 = tt15*tt19+tt21*tt1*tt19;
    tt25 = d_dWdy_dW;
    tt29 = -tt15*tt25-tt21*tt1*tt25;
    tt30 = SolnBlk.W[ii][jj].dmudT();
    tt33 = rho*t1;
    tt37 = 2.0*tt33*tt10-2.0/3.0*rho;
    tt38 = omega*omega;
    tt39 = 1/tt38;
    tt43 = SolnBlk.dWdy[ii][jj].v.x;
    tt44 =SolnBlk.dWdx[ii][jj].v.y; 
    tt45 = tt43+tt44;
    tt47 = d_dWdy_dW;
    tt51 = tt15*tt47+tt21*tt1*tt47;
    tt52 = d_dWdx_dW;
    tt56 = tt15*tt52+tt21*tt1*tt52;
    tt62 = Cp/Pr_turb;
    tt79 = SolnBlk.dWdx[ii][jj].k;
    tt80 = tt1*tt79;
    tt81 = sigma*k*tt80;
    tt84 = tt1*tt45;
    tt104 = d_dWdx_dW;
    tt120 = tt30*tt79;
    tt129 = sigma*rho;
    tt130 = tt129*tt80;
    tt133 = d_dWdx_dW;
    tt134 = (tt15+tt129*tt2)*tt133;
    tt136 = V*rho;
    tt140 = k*tt39;
    tt144 = tt129*tt140*tt79;
    tt164 = SolnBlk.dWdx[ii][jj].omega;
    tt165 = tt1*tt164;
    tt168 = sigma_star*rho;
    tt174 = d_dWdx_dW;
    
    dFvdW(2,0) =tt2*tt45;
    dFvdW(2,1) =tt51;
    dFvdW(2,2) =tt56;
    dFvdW(2,3) =tt30*tt45;
    dFvdW(2,4) =tt33*tt45;
    dFvdW(2,5) =-tt21*tt39*tt45;
    dFvdW(3,0) =tt62*tt2*tt63+Sum_q+tt81+U*tt14+V*k*tt84;
    dFvdW(3,1) =2.0*tt15*tt10+2.0*tt21*tt1*tt10-2.0/3.0*tt21+2.0*U*tt24+V*tt51;
    dFvdW(3,2) =2.0/3.0*U*tt29+tt15*tt45+tt21*tt84+V*tt56;
    dFvdW(3,3) =(kappa+tt62*tt21*tt1)*tt104+rho*Sum_dq+tt120+2.0*U*tt30*tt10+V*tt30*tt45;
    dFvdW(3,4) =tt62*tt33*tt63+tt130+tt134+U*tt37+tt136*tt84;
    dFvdW(3,5) =-tt62*rho*tt140*tt63-tt144-2.0*U*rho*tt140*tt10-tt136*tt140*tt45;
  } 
 
  
}
void dGvdW_Turbulent(DenseMatrix &dGvdW, Chem2D_Quad_Block &SolnBlk, 
				   const int &ii, const int &jj){
  
  double  d_dWdx_dW,  d_dWdy_dW;
  
  d_dWdx_dW = SolnBlk.d_dWdx_dW[ii][jj][0];
  d_dWdy_dW = SolnBlk.d_dWdy_dW[ii][jj][0];
      
  double t1,t2,t3,t4,t5,t7,t8,t10;
  double t13,t14,t18,t19,t21,t23,t24,t27,t29;
  double t31,t35,t36,t40,t41,t45,t51,t56,t57;
  double t85, t86, t87,t89, t110;
  double t140, t149, t150, t153, t154, t155; 
  double t160, t164, t194, t195, t198, t204; 
      
  double kappa, Cp,mu;
  double rho, U, V, k, omega;
  double mu_turb, kappa_turb, Pr_turb, Dm_turb, sigma, sigma_star;
  double Rmix = SolnBlk.W[ii][jj].Rtot();    
  double Temp =SolnBlk.W[ii][jj].T();
  
  kappa =SolnBlk.W[ii][jj].kappa();
  Cp =SolnBlk.W[ii][jj].Cp();
  mu  =SolnBlk.W[ii][jj].mu();
  
  mu_turb = SolnBlk.W[ii][jj].eddy_viscosity();  
  Pr_turb = SolnBlk.W[ii][jj].Pr_turb(); 
  Dm_turb = SolnBlk.W[ii][jj].Dm_turb(); 
  
  sigma = SolnBlk.W[0][0].sigma;
  sigma_star = SolnBlk.W[0][0].sigma_star;  
  
  rho = SolnBlk.W[ii][jj].rho;
  U =  SolnBlk.W[ii][jj].v.x;
  V =  SolnBlk.W[ii][jj].v.y;
  k = SolnBlk.W[ii][jj].k;
  omega = SolnBlk.W[ii][jj].omega;
  
  t1 = 1/max(omega,TOLER);
  t2 = k*t1;
  t3 =  SolnBlk.dWdy[ii][jj].v.x;
  t4 = SolnBlk.dWdx[ii][jj].v.y;
  t5 = t3+t4;
  t7 = mu;
  t8 = d_dWdy_dW;
  t10 = rho*k;
  t13 = t7*t8+t10*t1*t8;
  t14 =  d_dWdx_dW;
  t18 = t7*t14+t10*t1*t14;
  t19 = SolnBlk.W[ii][jj].dmudT();
  t21 = rho*t1;
  t23 = omega*omega;
  t24 = 1/t23;
  t27 =SolnBlk.dWdy[ii][jj].v.y;
  t29 = SolnBlk.dWdx[ii][jj].v.x;
  t31 = 2.0/3.0*t27-t29/3.0;
  t35 = 2.0*t2*t31-2.0/3.0*k;
  t36 =  d_dWdx_dW;
  t40 = -t7*t36-t10*t1*t36;
  t41 = d_dWdy_dW;
  t45 = t7*t41+t10*t1*t41;
  t51 = 2.0*t21*t31-2.0/3.0*rho;
  t56 = Cp/Pr_turb;
  //t57 =  SolnBlk.dTdy(ii,jj);
  t57= (ONE/(SolnBlk.W[ii][jj].rho*Rmix)) * (SolnBlk.dWdy[ii][jj].p - 
                 (SolnBlk.W[ii][jj].p/SolnBlk.W[ii][jj].rho)*SolnBlk.dWdy[ii][jj].rho);
  
  t85 =SolnBlk.dWdy[ii][jj].k;
  t86 = t1*t85;
  t87 = sigma_star*k*t86;
  t89 = t1*t5;
  t110 =  d_dWdy_dW;

  t140 = t19*t85;
  t149 = sigma_star*rho;
  t150 = t149*t86;
  t153 =  d_dWdy_dW;
  t154 = (t7+t149*t2)*t153;
  t155 = U*rho;
  t160 = k*t24;
  t164 = t149*t160*t85;
      
  t194 =  SolnBlk.dWdy[ii][jj].omega;
  t195 = t1*t194;
  t198 = sigma*rho;
  t204 = d_dWdy_dW;
  
  dGvdW(1,0) = t2*t5;
  dGvdW(1,1) = t13;
  dGvdW(1,2) = t18;
  dGvdW(1,3) = t19*t5;
  dGvdW(1,4) = t21*t5;
  dGvdW(1,5) = -t10*t24*t5;
    
  dGvdW(2,0) = t35;
  dGvdW(2,1) = 2.0/3.0*t40;
  dGvdW(2,2) = 4.0/3.0*t45;
  dGvdW(2,3) = 2.0*t19*t31;
  dGvdW(2,4) = t51;
  dGvdW(2,5) = -2.0*t10*t24*t31;
  double Sum_q = 0.0;
  double Sum_dq = 0.0;
  int ns = SolnBlk.W[0][0].ns;
  // 6 represents rho, vr,vz, p, k, omega in 2D axisymmetric turbulent flows
  
  for(int Num = 0; Num<ns; Num++)
    {
      //for each species	// h*(Dm+Dt)*gradc  
      Sum_q +=  (SolnBlk.W[ii][jj].specdata[Num].Enthalpy(Temp)+SolnBlk.W[ii][jj].specdata[Num].Heatofform())*(SolnBlk.W[ii][jj].spec[Num].diffusion_coef+Dm_turb)*SolnBlk.dWdy[ii][jj].spec[Num].c;
      //dhdT *(Dm+Dmt)*gradc
      Sum_dq +=  SolnBlk.W[ii][jj].specdata[Num].Enthalpy_prime(Temp)*(SolnBlk.W[ii][jj].spec[Num].diffusion_coef+Dm_turb)*SolnBlk.dWdy[ii][jj].spec[Num].c;
      
    }
  
  dGvdW(3,0) = t56*t2*t57+Sum_q+t87+U*k*t89+V*t35;
  dGvdW(3,1) = t7*t5+t10*t89+U*t13+2.0/3.0*V*t40;
  dGvdW(3,2) = U*t18+2.0*t7*t31+2.0*t10*t1*t31-2.0/3.0*t10+4.0/3.0*V*t45;
  dGvdW(3,3) = (kappa+t56*t10*t1)*t110+ rho*Sum_dq+t140+U*t19*t5+2.0*V*t19*t31;
  dGvdW(3,4) = t56*t21*t57+t150+t154+t155*t89+V*t51;
  dGvdW(3,5) = -t56*rho*t160*t57-t164-t155*t160*t5-2.0*V*rho*t160*t31;
  
  int NUM_VAR = NUM_CHEM2D_VAR_SANS_SPECIES; 
  for(int Num = 0; Num<(ns-1); Num++)
    {//h*rho*(D+Dt)*d_dWdy_dW
      dGvdW(3, NUM_VAR+Num) = (SolnBlk.W[ii][jj].specdata[Num].Enthalpy(Temp)+SolnBlk.W[ii][jj].specdata[Num].Heatofform())*rho*(SolnBlk.U[ii][jj].rhospec[Num].diffusion_coef+Dm_turb)*d_dWdy_dW;
      //(D+Dt)*gradc
      dGvdW(NUM_VAR+Num, 0) =(SolnBlk.U[ii][jj].rhospec[Num].diffusion_coef+Dm_turb)*SolnBlk.dWdy[ii][jj].spec[Num].c;
      //rho*(D+Dt)*dcydc
      dGvdW(NUM_VAR+Num,NUM_VAR+Num) = rho*(SolnBlk.U[ii][jj].rhospec[Num].diffusion_coef+Dm_turb)*d_dWdy_dW;
    }// end the center ones
  
  dGvdW(4,0) = t87;
  dGvdW(4,3) = t140;
  dGvdW(4,4) = t150+t154;
  dGvdW(4,5) = -t164;
  
  dGvdW(5,0) = sigma*k*t195;
  dGvdW(5,3) = t19*t194;
  dGvdW(5,4) = t198*t195;
  dGvdW(5,5) = -t198*t160*t194+(t7+t198*t2)*t204;
  
  
  //The follwoing entries are different for axisymmetric flow and planar flow
  // For matrix row 2 and row 3
  if (SolnBlk.Axisymmetric == 1) {
    
    double term1,term2,term3,term4, term5,term7,term8,term10;
    double term13, term14,term18,term19;
    double term21,term22,term23,term24,term27,term29,term31,term34, term38, term39, term43, term44;
    double term47, term51,term57, term62, term63, term91, term92,term93, term95, term116; 
    double term146, term155, term156, term159, term160, term161, term166,term170 , term200, term201, term204, term210; 
    
    double r = SolnBlk.Grid.Cell[ii][jj].Xc.y;
    
    term1 = 1/max(omega,TOLER);
    term2 = k*term1;
    term3 =  SolnBlk.dWdy[ii][jj].v.x; 
    term4 =  SolnBlk.dWdx[ii][jj].v.y;
    term5 = term3+term4;
    term7 = mu;
    term8 = d_dWdy_dW;
    term10 = rho*k;
    term13 = term7*term8+term10*term1*term8;
    term14 = d_dWdx_dW;
    term18 = term7*term14+term10*term1*term14;
    term19 = SolnBlk.W[ii][jj].dmudT();
    term21 = rho*term1;
    term23 = omega*omega;
    term24 = 1/term23;
    term27 =  SolnBlk.dWdy[ii][jj].v.y;
    term29 = SolnBlk.dWdx[ii][jj].v.x; 
    term31 = 1/r;
    term34 = 2.0/3.0*term27-term29/3.0-V*term31/3.0;
    term38 = 2.0*term2*term34-2.0/3.0*k;
    term39 = d_dWdx_dW;
    term43 = -term7*term39-term10*term1*term39;
    term44 = d_dWdy_dW;
    term47 = 2.0/3.0*term44-term31/3.0;
    term51 = term7*term47+term10*term1*term47;
    term57 = 2.0*term21*term34-2.0/3.0*rho;
    term62 = Cp/Pr_turb;
    // term63 = SolnBlk.dTdy(ii,jj);
    term63 = (ONE/(SolnBlk.W[ii][jj].rho*Rmix)) * (SolnBlk.dWdy[ii][jj].p - 
	    (SolnBlk.W[ii][jj].p/SolnBlk.W[ii][jj].rho)*SolnBlk.dWdy[ii][jj].rho);
    term91 =SolnBlk.dWdy[ii][jj].k;
    term92 = term1*term91;
    term93 = sigma*k*term92;
    term95 = term1*term5;
    term116 = d_dWdy_dW;
    term146 = term19*term91;
    term155 = sigma*rho;
    term156 = term155*term92;
    term159 = d_dWdy_dW;
    term160 = (term7+term155*term2)*term159;
    term161 = U*rho;
    term166 = k*term24;
    term170 = term155*term166*term91;
    
    term200 =SolnBlk.dWdy[ii][jj].omega;
    term201 = term1*term200;
    term204 = sigma_star*rho;
    term210 = d_dWdy_dW;
    
    dGvdW(2,0) = term38;
    dGvdW(2,1) = 2.0/3.0*term43;
    dGvdW(2,2) = 2.0*term51;
    dGvdW(2,3) = 2.0*term19*term34;
    dGvdW(2,4) = term57;
    dGvdW(2,5) = -2.0*term10*term24*term34;
    
    dGvdW(3,0) = term62*term2*term63+Sum_q+term93+U*k*term95+V*term38;
    dGvdW(3,1) = term7*term5+term10*term95+U*term13+2.0/3.0*V*term43;
    dGvdW(3,2) = U*term18+2.0*term7*term34+2.0*term10*term1*term34-2.0/3.0*term10+2.0*V*term51;
    dGvdW(3,3) = (kappa+term62*term10*term1)*term116+ rho*Sum_dq+term146+U*term19*term5+2.0*V*term19*term34;
    dGvdW(3,4) = term62*term21*term63+term156+term160+term161*term95+V*term57;
    dGvdW(3,5) = -term62*rho*term166*term63-term170-term161*term166*term5-2.0*V*rho*term166*term34;
    
  }
 if (SolnBlk.Axisymmetric == 2) {

   double tt1,tt2,tt3,tt4,tt5,tt7,tt8, tt10,tt13,tt14,tt18;
   double tt19,tt21,tt23,tt24,tt27,tt29,tt31;
   double tt34,tt38,tt39,tt40,tt44,tt45,tt49,tt55, tt60;
   double tt61,tt77,tt78,tt79,tt81, tt102;
   double tt118,tt127,tt128,tt131,tt132,tt133, tt138,tt142;
   double tt150,tt153,tt155,tt158,tt162,tt163, tt166,tt172;
   
   
   double r = SolnBlk.Grid.Cell[ii][jj].Xc.x;
 
   tt1 = 1/max(omega,TOLER);
   tt2 = k*tt1;
   tt3 = SolnBlk.dWdy[ii][jj].v.x;
   tt4 = SolnBlk.dWdx[ii][jj].v.y; 
   tt5 = tt3+tt4;
   tt7 = SolnBlk.W[ii][jj].mu();
   tt8 = d_dWdy_dW;
   tt10 = rho*k;
   tt13 = tt7*tt8+tt10*tt1*tt8;
   tt14 = d_dWdx_dW;
   tt18 = tt7*tt14+tt10*tt1*tt14;
   tt19 =  SolnBlk.W[ii][jj].dmudT();
   tt21 = rho*tt1;
   tt23 = omega*omega;
   tt24 = 1/tt23;
   tt27 = SolnBlk.dWdy[ii][jj].v.y;
   tt29 = SolnBlk.dWdx[ii][jj].v.x;
   tt31 = 1/r;
   tt34 = 2.0/3.0*tt27-tt29/3.0-U*tt31/3.0;
   tt38 = 2.0*tt2*tt34-2.0/3.0*k;
   tt39 = d_dWdx_dW;
   tt40 = -tt39-tt31;
   tt44 = tt7*tt40/3.0+tt10*tt1*tt40/3.0;
   tt45 = d_dWdy_dW;
   tt49 = tt7*tt45+tt10*tt1*tt45;
   tt55 = 2.0*tt21*tt34-2.0/3.0*rho;
   tt60 = Cp/Pr_turb;
   // tt61 =SolnBlk.dTdy(ii,jj);
   tt61= (ONE/(SolnBlk.W[ii][jj].rho*Rmix)) * (SolnBlk.dWdy[ii][jj].p - 
                 (SolnBlk.W[ii][jj].p/SolnBlk.W[ii][jj].rho)*SolnBlk.dWdy[ii][jj].rho);
   tt77 = SolnBlk.dWdy[ii][jj].k;
   tt78 = tt1*tt77;
   tt79 = sigma*k*tt78;
   tt81 = tt1*tt5;
   tt102 = d_dWdy_dW;
   tt127 = sigma*rho;
   tt128 = tt127*tt78;
   tt131 = d_dWdy_dW;
   tt132 = (tt7+tt127*tt2)*tt131;
   tt133 = U*rho;
   tt138 = k*tt24;
   tt142 = tt127*tt138*tt77;
   tt162 = SolnBlk.dWdy[ii][jj].omega;
   tt163 = tt1*tt162;
   tt166 = sigma_star*rho;
   tt172 = d_dWdy_dW;
    
    
   dGvdW(2,0) = tt38;
   dGvdW(2,1) = 2.0*tt44;
   dGvdW(2,2) = 4.0/3.0*tt49;
   dGvdW(2,3) = 2.0*tt19*tt34;
   dGvdW(2,4) = tt55;
   dGvdW(2,5) = -2.0*tt10*tt24*tt34;
    
   dGvdW(3,0) = tt60*tt2*tt61+Sum_q+tt79+U*k*tt81+V*tt38;
   dGvdW(3,1) = tt7*tt5+tt10*tt81+U*tt13+2.0*V*tt44;
   dGvdW(3,2) = U*tt18+2.0*tt7*tt34+2.0*tt10*tt1*tt34-2.0/3.0*tt10+4.0/3.0*V*tt49;
   dGvdW(3,3) = (kappa+tt60*tt1*tt10)*tt102+ rho*Sum_dq+tt118+U*tt19*t5+2.0*V*tt19*tt34;
   dGvdW(3,4) = tt60*tt21*tt61+tt128+tt132+tt133*tt81+V*tt55;
   dGvdW(3,5) = -tt60*rho*tt138*tt61-tt142-tt133*tt138*tt5-2.0*V*rho*tt138*tt34;
    
     
 }
 
}
//The following two neigbour Jacobians are defined ... because ...
//when computing the viscous flux through cell faces by  averaging quantities of the left (i,j) and its' neigbour cell  
void dFvdW_Turbulent_Neigbour(DenseMatrix &dFvdW,  Chem2D_Quad_Block &SolnBlk, 
				     const string &Orient, const int &ii, const int &jj){
  double d_dWdx_dW, d_dWdy_dW;
  int i, j;
  i = ii;
  j=jj;
  
  if( Orient == "NORTH"){
    j = jj+1;
    d_dWdx_dW = SolnBlk.d_dWdx_dW[ii][jj][1];
    d_dWdy_dW = SolnBlk.d_dWdy_dW[ii][jj][1];
      
  }else if( Orient == "SOUTH"){
    j = jj-1;
    d_dWdx_dW = SolnBlk.d_dWdx_dW[ii][jj][2];
    d_dWdy_dW = SolnBlk.d_dWdy_dW[ii][jj][2];
    
  }else if( Orient == "WEST"){
    i = ii-1;
    d_dWdx_dW = SolnBlk.d_dWdx_dW[ii][jj][3];
    d_dWdy_dW = SolnBlk.d_dWdy_dW[ii][jj][3];
    
  }else if( Orient == "EAST"){
    i = ii+1;
    d_dWdx_dW = SolnBlk.d_dWdx_dW[ii][jj][4];
    d_dWdy_dW = SolnBlk.d_dWdy_dW[ii][jj][4];
    
  }
  if (i < SolnBlk.ICl || i > SolnBlk.ICu ||
      j < SolnBlk.JCl || j > SolnBlk.JCu) {
    // GHOST CELL
  } else {
    //     // NON-GHOST CELL.  
    double t1,t3,t4,t5,t6,t8, t10, t12,t16, t27, t30,t33;
    double t34,t35, t36, t38;
    double t40,t41, t43, t45, t46, t48;
    double t50,t51, t53, t62;   
    
    double kappa, Cp, mu;
    double rho, U, V, k, omega;
    double mu_turb, kappa_turb, Pr_turb, Dm_turb, sigma, sigma_star;
    
    double Temp =SolnBlk.W[i][j].T();
    
    kappa =SolnBlk.W[i][j].kappa();
    Cp =SolnBlk.W[i][j].Cp();
    mu  =SolnBlk.W[i][j].mu();
    
    mu_turb = SolnBlk.W[i][j].eddy_viscosity();  
    Pr_turb = SolnBlk.W[i][j].Pr_turb(); 
    Dm_turb = SolnBlk.W[i][j].Dm_turb(); 
    sigma = SolnBlk.W[0][0].sigma;
    sigma_star = SolnBlk.W[0][0].sigma_star;  
    
    rho = SolnBlk.W[i][j].rho;
    U =  SolnBlk.W[i][j].v.x;
    V =  SolnBlk.W[i][j].v.y;
    k = SolnBlk.W[i][j].k;
    omega = SolnBlk.W[i][j].omega;
    
    t1 = mu;
    t3 = 1/max(omega,TOLER);
    t4 = rho*k*t3;
    t5 = t1+t4;
    t6 = d_dWdx_dW;
    t8 =  d_dWdy_dW;
    t10 =  d_dWdy_dW;
    t12 = d_dWdx_dW;
    t16 = V*t5;
    t27 = d_dWdx_dW;
    t30 = k*t3;
    t33 =d_dWdx_dW;
    t34 = (t1+sigma*rho*t30)*t33;
    t62 =  d_dWdx_dW;
    
    dFvdW(1,1) = 4.0/3.0*t5*t6;
    dFvdW(1,2) = -2.0/3.0*t5*t8;
    
    dFvdW(2,1) = t5*t10;
    dFvdW(2,2) = t5*t12;
    
    dFvdW(3,1) = 4.0/3.0*U*t5*t6+t10*t16;
    dFvdW(3,2) = t16*t12-2.0/3.0*U*t5*t8;
    dFvdW(3,3) = (kappa+Cp/Pr_turb*t4)*t27;
    dFvdW(3,4) = t34;

    int ns = SolnBlk.W[0][0].ns;
    int NUM_VAR = NUM_CHEM2D_VAR_SANS_SPECIES; 
    for(int Num = 0; Num<(ns-1); Num++)
      {//h*rho*(D+Dt)*dcxdc;
	dFvdW(3,NUM_VAR+Num) =(SolnBlk.W[i][j].specdata[Num].Enthalpy(Temp)+SolnBlk.W[i][j].specdata[Num].Heatofform())*rho*(SolnBlk.U[ii][jj].rhospec[Num].diffusion_coef+Dm_turb)*d_dWdx_dW;
	//rho*(D+Dt)*dcxdc;
	dFvdW(NUM_VAR+Num,NUM_VAR+Num) = rho*(SolnBlk.U[i][j].rhospec[Num].diffusion_coef+Dm_turb)*d_dWdx_dW;
	
      }//endof the neigbours
 
  dFvdW(4,4) = t34;
  dFvdW(5,5) = (t1+sigma*rho*t30)*t62;
  
  }
  
}
 
void dGvdW_Turbulent_Neigbour(DenseMatrix &dGvdW,  Chem2D_Quad_Block &SolnBlk, 
				     const string &Orient, const int &ii, const int &jj){
  double d_dWdx_dW, d_dWdy_dW;
  int i, j;
  i = ii;
  j=jj;

  if( Orient == "NORTH"){
    j = jj+1;
    d_dWdx_dW = SolnBlk.d_dWdx_dW[ii][jj][1];
    d_dWdy_dW = SolnBlk.d_dWdy_dW[ii][jj][1];

  }else if( Orient == "SOUTH"){
    j = jj-1;
    d_dWdx_dW = SolnBlk.d_dWdx_dW[ii][jj][2];
    d_dWdy_dW = SolnBlk.d_dWdy_dW[ii][jj][2];

  }else if( Orient == "WEST"){
    i = ii-1;
    d_dWdx_dW = SolnBlk.d_dWdx_dW[ii][jj][3];
    d_dWdy_dW = SolnBlk.d_dWdy_dW[ii][jj][3];

  }else if( Orient == "EAST"){
    i = ii+1;
    d_dWdx_dW = SolnBlk.d_dWdx_dW[ii][jj][4];
    d_dWdy_dW = SolnBlk.d_dWdy_dW[ii][jj][4];

  }
  if (i < SolnBlk.ICl || i > SolnBlk.ICu ||
      j < SolnBlk.JCl || j > SolnBlk.JCu) {
    // GHOST CELL
  } else {
    //     // NON-GHOST CELL. 
    
    double t1,t3,t4,t5,t6,t8,t10, t12;
    double t16,t27,t30,t33;
    double t59,t61,t62,t71;
    
    
    double kappa, Cp, mu;
    double rho, U, V, k, omega;
    double mu_turb, kappa_turb, Pr_turb, Dm_turb, sigma, sigma_star;
    double Temp =SolnBlk.W[i][j].T();
    
    kappa =SolnBlk.W[i][j].kappa();
    Cp =SolnBlk.W[i][j].Cp();
    mu  =SolnBlk.W[i][j].mu();
    
    mu_turb = SolnBlk.W[i][j].eddy_viscosity();  
    Pr_turb = SolnBlk.W[i][j].Pr_turb(); 
    Dm_turb = SolnBlk.W[i][j].Dm_turb(); 
    
    sigma = SolnBlk.W[0][0].sigma;
    sigma_star = SolnBlk.W[0][0].sigma_star;  
    
    rho = SolnBlk.W[i][j].rho;
    U =  SolnBlk.W[i][j].v.x;
    V =  SolnBlk.W[i][j].v.y;
    k = SolnBlk.W[i][j].k;
    omega = SolnBlk.W[i][j].omega;
	
    
    t1 = mu;
    t3 = 1/max(omega,TOLER);
    t4 = rho*k*t3;
    t5 = t1+t4;
    t6 = d_dWdy_dW;
    t8 = d_dWdx_dW;
    t10 = d_dWdx_dW;
    t12 = d_dWdy_dW;
    t16 = U*t5;
    t27 = d_dWdy_dW;
    t30 = k*t3;
    t33 = d_dWdy_dW;
  
    t59 = rho;
    t61 = k;
    t62 = omega;
    t71 = d_dWdy_dW;
   
    dGvdW(1,1) = t5*t6;
    dGvdW(1,2) = t5*t8;
    dGvdW(2,1) = -2.0/3.0*t5*t10;
    dGvdW(2,2) = 4.0/3.0*t5*t12;
    dGvdW(3,1) = -2.0/3.0*V*t5*t10+t16*t6;
    dGvdW(3,2) = t16*t8+4.0/3.0*V*t5*t12;
    dGvdW(3,3) = (kappa+Cp/Pr_turb*t4)*t27;
    dGvdW(3,4) = (t1+sigma*rho*t30)*t33;
    int ns = SolnBlk.W[0][0].ns;
    int NUM_VAR = NUM_CHEM2D_VAR_SANS_SPECIES; 
    for(int Num = 0; Num<(ns-1); Num++)
      {//h*rho*(D+Dt)*dcydc;
	dGvdW(3,NUM_VAR+Num) =(SolnBlk.W[i][j].specdata[Num].Enthalpy(Temp)+SolnBlk.W[i][j].specdata[Num].Heatofform())*rho*(SolnBlk.U[ii][jj].rhospec[Num].diffusion_coef+Dm_turb)*d_dWdy_dW;
	//rho*(D+Dt)*dcydc;
	dGvdW(NUM_VAR+Num,NUM_VAR+Num) = rho*(SolnBlk.U[i][j].rhospec[Num].diffusion_coef+Dm_turb)*d_dWdy_dW;
	
      }//endof the neigbours
    
    dGvdW(4,4) = (t1+sigma*t59*t61/t62)*t33;
    
    dGvdW(5,5) = (t1+sigma*rho*t30)*t71;

  }
 
} 

//Turbulence source Jacobian
void dStdW_Turbulence(DenseMatrix &dStdW, Chem2D_Quad_Block &SolnBlk, const int &ii, const int &jj){


  double t1,t2,t3,t5,t7,t10;
  double t12,t13,t14,t15,t16,t17,t20;
  double t24,t28,t29,t30,t37,t38,t40,t41, t48;
  double t49,t50,t54,t64,t66, t67,t70, t72, t73, t77, t78; 
  double t81, t82, t86, t88, t92, t95 , t96 , t106, t112; 
  
  double alpha, beta, beta_star;
  double rho, U, V, k, omega;
  
  double d_dWdx_dW, d_dWdy_dW ;
 
  d_dWdx_dW = SolnBlk.d_dWdx_dW[ii][jj][0];
  d_dWdy_dW = SolnBlk.d_dWdy_dW[ii][jj][0];
  
   
  beta = SolnBlk.W[0][0].beta;
  beta_star = SolnBlk.W[0][0].beta_star;
  alpha =SolnBlk.W[0][0].alpha;
  
  t1 = 1/max(omega,TOLER);
  t2 = k*t1;
  t3 = SolnBlk.dWdx[ii][jj].v.x;
  t5 = SolnBlk.dWdy[ii][jj].v.y;
  t7 = 2.0/3.0*t3-t5/3.0;
  t10 = 2.0/3.0*k;
  t12 = (2.0*t2*t7-t10)*t3;
  t13 = SolnBlk.dWdy[ii][jj].v.x;
  t14 = SolnBlk.dWdx[ii][jj].v.y;
  t15 = t13+t14;
  t16 = t15*t15;
  t17 = t2*t16;
  t20 = 2.0/3.0*t5-t3/3.0;
  t24 = (2.0*t2*t20-t10)*t5;
  t28 = rho*k;
  t29 = d_dWdx_dW;
  t30 = t1*t29;
  t37 = 2.0/3.0*t28;
  t38 = 2.0*t28*t1*t7-t37;
  t40 = t1*t15;
  t41 = d_dWdy_dW;
  t48 = 4.0/3.0*t28*t30*t3+t38*t29+2.0*t28*t40*t41-2.0/3.0*t28*t30*t5;
  t49 = d_dWdy_dW;
  t50 = t1*t49;
  t54 =d_dWdx_dW;
  t64 = 2.0*t28*t1*t20-t37;
  t66 = -2.0/3.0*t28*t50*t3+2.0*t28*t40*t54+4.0/3.0*t28*t50*t5+t64*t49;
  t67 = rho*t1;
  t70 = 2.0/3.0*rho;
  t72 = (2.0*t67*t7-t70)*t3;
  t73 = t67*t16;
  t77 = (2.0*t67*t20-t70)*t5;
  t78 = beta_star*rho;
  t81 = max(omega,TOLER)*max(omega,TOLER);
  t82 = 1/t81;
  t86 = 2.0*t28*t82*t7*t3;
  t88 = t28*t82*t16;
  t92 = 2.0*t28*t82*t20*t5;
  t95 = alpha*omega;
  t96 = 1/max(k, TOLER);
  t106 =k*k;
  t112 = t38*t3+t28*t1*t16+t64*t5;
  
  dStdW(4,0) =  t12+t17+t24-beta_star*k*omega;
  dStdW(4,1) =  t48;
  dStdW(4,2) =  t66;
   
  dStdW(4,4) =  t72+t73+t77-t78*omega;
  dStdW(4,5) =  -t86-t88-t92-t78*k;
  
  dStdW(5,0) =  t95*t96*(t12+t17+t24)-beta*t81;
  dStdW(5,1) =  t95*t96*t48;
  dStdW(5,2) =  t95*t96*t66;
  
  dStdW(5,4) =  -t95/t106*t112+t95*t96*(t72+t73+t77);
  dStdW(5,5) =  alpha*t96*t112+t95*t96*(-t86-t88-t92)-2.0*beta*rho*omega;
//   //Axisymmetric casee
  if(SolnBlk.Axisymmetric){ 
    double r;
    if (SolnBlk.Axisymmetric == 1) {
      
      r = SolnBlk.Grid.Cell[ii][jj].Xc.y;
      double tt1,tt2,tt3,tt4,tt5,tt6, tt9,tt10;
      double tt16,tt23,tt24,tt25,tt26, tt27, tt30;
      double tt37,tt52,tt54,tt55,tt66,tt67,tt69;
      double tt78,tt79, tt80,tt88, tt96; 
      

      tt1 = 1/max(omega, TOLER);
      tt2 = k*tt1;
      tt3 = 1/r;
      tt4 = V*t3;
      tt5 = SolnBlk.dWdx[ii][jj].v.x;
      tt6 = tt4*tt5;
      tt9 = SolnBlk.dWdy[ii][jj].v.y;
      tt10 = tt4*tt9;
      tt16 = 2.0/3.0*tt4-tt5/3.0-tt9/3.0;
      tt23 = -2.0/3.0*tt2*tt6-2.0/3.0*tt2*tt10+(2.0*tt2*tt16-2.0/3.0*k)*V*tt3;
      tt24 = rho*k;
      tt25 = tt24*tt1;
      tt26 = d_dWdx_dW;
      tt27 = tt4*tt26;
      tt30 = tt1*tt3;
      tt37 = d_dWdy_dW;
      tt52 = 2.0*tt24*tt1*tt16-2.0/3.0*tt24;
      tt54 = -2.0/3.0*tt24*tt30*tt5-2.0/3.0*tt24*tt30*tt9-2.0/3.0*tt25*tt4*tt37+2.0*tt25*(
2.0/3.0*tt3-tt37/3.0)*V*tt3+tt52*tt3;
      tt55 = rho*tt1;
      tt66 = -2.0/3.0*tt55*tt6-2.0/3.0*tt55*tt10+(2.0*tt55*tt16-2.0/3.0*rho)*V*tt3;
      tt67 = max(omega, TOLER)*max(omega, TOLER);
      tt69 = tt24/tt67;
      tt78 = 2.0/3.0*tt69*tt6+2.0/3.0*tt69*tt10-2.0*tt69*tt16*V*tt3;
      tt79 = alpha*omega;
      tt80 = 1/max(k,TOLER);
      tt88 = max(k,TOLER)*max(k,TOLER);
      tt96 = -2.0/3.0*tt25*tt6-2.0/3.0*tt25*tt10+tt52*V*tt3;
     
      dStdW(4,0) +=tt23;
      dStdW(4,1) +=-4.0/3.0*tt25*tt27;
      dStdW(4,2) +=tt54;
      dStdW(4,4) +=tt66;
      dStdW(4,5) +=tt78;
   
      dStdW(5,0) +=tt79*tt80*tt23;
      dStdW(5,1) +=-4.0/3.0*alpha*rho*tt27;
      dStdW(5,2) +=tt79*tt80*tt54;
      dStdW(5,4) +=-tt79/tt88*tt96+tt79*tt80*tt66;
      dStdW(5,5) +=alpha*tt80*tt96+tt79*tt80*tt78;
      
    }//end of axisymmetric case (1) 
    if (SolnBlk.Axisymmetric == 2) {
      r = SolnBlk.Grid.Cell[ii][jj].Xc.x;
     
      double tt1,tt2,tt3,tt4,tt5,tt6, tt9,tt10;
      double tt16,tt23,tt24,tt25,tt29,tt30;
      double tt48,tt50,tt51,tt52,tt55,tt66,tt67,tt69;
      double tt78,tt79, tt80,tt88, tt96; 
     
      tt1 = 1/max(omega,TOLER);
      tt2 = k*tt1;
      tt3 = 1/r;
      tt4 = U*tt3;
      tt5 = SolnBlk.dWdx[ii][jj].v.x;
      tt6 = tt4*tt5;
      tt9 = SolnBlk.dWdy[ii][jj].v.y;
      tt10 = tt4*tt9;
      tt16 = 2.0/3.0*tt4-tt5/3.0-tt9/3.0;
      tt23 = -2.0/3.0*tt2*tt6-2.0/3.0*tt2*tt10+(2.0*tt2*tt16-2.0/3.0*k)*U*tt3;
      tt24 = rho*k;
      tt25 = tt1*tt3;
      tt29 = tt24*tt1;
      tt30 = d_dWdx_dW; 
      tt48 = 2.0*tt24*tt1*tt16-2.0/3.0*tt24;
      tt50 = -2.0/3.0*tt24*tt25*tt5-2.0/3.0*tt29*tt4*tt30-2.0/3.0*tt24*tt25*tt9+2.0*tt29*(
2.0/3.0*tt3-tt30/3.0)*U*tt3+tt48*tt3;
      tt51 = d_dWdy_dW;
      tt52 = tt4*tt51;
      tt55 = rho*tt1;
      tt66 = -2.0/3.0*tt55*tt6-2.0/3.0*tt55*tt10+(2.0*tt55*tt16-2.0/3.0*rho)*U*tt3;
      tt67 = max(omega, TOLER)*max(omega,TOLER);
      tt69 = tt24/tt67;
      tt78 = 2.0/3.0*tt69*tt6+2.0/3.0*tt69*tt10-2.0*tt69*tt16*U*tt3;
      tt79 = alpha*omega;
      tt80 = 1/max(k,TOLER);
      tt88 = max(k,TOLER)*max(k,TOLER);
      tt96 = -2.0/3.0*tt29*tt6-2.0/3.0*tt29*tt10+tt48*U*tt3;
     
      dStdW(4,0) += tt23;
      dStdW(4,1) += tt50;
      dStdW(4,2) += -4.0/3.0*tt29*tt52;
      dStdW(4,4) += tt66;
      dStdW(4,5) += tt78;
      dStdW(5,0) += tt79*tt80*tt23;
      dStdW(5,1) += tt79*tt80*tt50;
      dStdW(5,2) += -4.0/3.0*alpha*rho*tt52;
      dStdW(5,4) += -tt79/tt88*tt96+tt79*tt80*tt66;
      dStdW(5,5) == alpha*tt80*tt96+tt79*tt80*tt78;

    }//end of axisymmetric case (2)
  }

    
} 
// /********************************************************
//  * Routine: PointImplicitBlkJ_Viscous                   *
//  *                                                      *
//  * This routine returns the viscous components of       *    
//  * Point Implicit Block Jacobian matrix for the         *
//  * specified local solution block.                      *
//  *                                                      *
//  ********************************************************/
DenseMatrix dGVdW_Viscous(Chem2D_Quad_Block &SolnBlk, Chem2D_Input_Parameters &Input_Parameters,
			  const int &ii, const int &jj){
  
  int NUM_VAR_CHEM2D = SolnBlk.NumVar()-1; 

  
  DenseMatrix dFvdW_N(NUM_VAR_CHEM2D,NUM_VAR_CHEM2D),
              dFvdW_S(NUM_VAR_CHEM2D,NUM_VAR_CHEM2D),
              dFvdW_E(NUM_VAR_CHEM2D,NUM_VAR_CHEM2D),
              dFvdW_W(NUM_VAR_CHEM2D,NUM_VAR_CHEM2D);
  DenseMatrix dGvdW_N(NUM_VAR_CHEM2D,NUM_VAR_CHEM2D),
              dGvdW_S(NUM_VAR_CHEM2D,NUM_VAR_CHEM2D),
              dGvdW_E(NUM_VAR_CHEM2D,NUM_VAR_CHEM2D),
              dGvdW_W(NUM_VAR_CHEM2D,NUM_VAR_CHEM2D);

  DenseMatrix dFvdW_Center(NUM_VAR_CHEM2D,NUM_VAR_CHEM2D);
  DenseMatrix dGvdW_Center(NUM_VAR_CHEM2D,NUM_VAR_CHEM2D);
  DenseMatrix dGVdW(NUM_VAR_CHEM2D,NUM_VAR_CHEM2D);
  DenseMatrix dWf_dWc(NUM_VAR_CHEM2D,NUM_VAR_CHEM2D);
  Vector2D nface_N, nface_S, nface_E, nface_W;

  dWf_dWc.zero();
  string Orient;


  if (ii < SolnBlk.ICl || ii > SolnBlk.ICu ||
      jj < SolnBlk.JCl || jj > SolnBlk.JCu) {
    // GHOST CELL
    dGVdW.zero();
    return dGVdW;
    
  } else {
//     // NON-GHOST CELL.

    dGVdW.zero();

  nface_N = SolnBlk.Grid.nfaceN(ii, jj);
  nface_S = SolnBlk.Grid.nfaceS(ii, jj);
  nface_E = SolnBlk.Grid.nfaceE(ii, jj);
  nface_W = SolnBlk.Grid.nfaceW(ii, jj);

  // Compute the flux Jacobians at cell centers
  dFvdW_Center.zero();
  dGvdW_Center.zero();

  if (Input_Parameters.i_Viscous_Flux_Evaluation != VISCOUS_RECONSTRUCTION_DIAMOND_PATH) {

    if (SolnBlk.Flow_Type == FLOWTYPE_LAMINAR) {

      dFvdW_Laminar(dFvdW_Center, SolnBlk,  ii,jj);
      dGvdW_Laminar(dGvdW_Center, SolnBlk, ii,jj);
      //Viscous flux Jacobians at North face of cell (ii, jj)
      dFvdW_N.zero();
      dGvdW_N.zero();
      Orient =  "NORTH";
      //  cell (ii, jj+1)
      dFvdW_Laminar_Neigbour( dFvdW_N, SolnBlk,Orient,ii,jj);  
      dGvdW_Laminar_Neigbour( dGvdW_N, SolnBlk,Orient,ii,jj);
      dFvdW_N = HALF*(dFvdW_N+dFvdW_Center);
      dGvdW_N = HALF*(dGvdW_N+dGvdW_Center);

      // Viscous flux Jacobians at South face of cell (ii, jj)
      dFvdW_S.zero();
      dGvdW_S.zero();
      Orient =  "SOUTH";
      dFvdW_Laminar_Neigbour( dFvdW_S, SolnBlk,Orient,ii,jj);
      dGvdW_Laminar_Neigbour( dGvdW_S, SolnBlk,Orient,ii,jj);
      dFvdW_S = HALF*(dFvdW_S+dFvdW_Center);
      dGvdW_S = HALF*(dGvdW_S+dGvdW_Center);

      // Viscous flux Jacobians at West face of cell (ii, jj)
      dFvdW_W.zero();
      dGvdW_W.zero();
      Orient =  "WEST";
      dFvdW_Laminar_Neigbour( dFvdW_W, SolnBlk,Orient,ii,jj);
      dGvdW_Laminar_Neigbour( dGvdW_W, SolnBlk,Orient,ii,jj);
      dFvdW_W = HALF*(dFvdW_W+dFvdW_Center);
      dGvdW_W = HALF*(dGvdW_W+dGvdW_Center);
      // Viscous flux Jacobians at East face of cell (ii, jj)
      dFvdW_E.zero();
      dGvdW_E.zero();
      Orient =  "EAST";
      dFvdW_Laminar_Neigbour( dFvdW_E, SolnBlk,Orient,ii,jj);
      dGvdW_Laminar_Neigbour( dGvdW_E, SolnBlk,Orient,ii,jj);
      dFvdW_E = HALF*(dFvdW_E+dFvdW_Center);
      dGvdW_E = HALF*(dGvdW_E+dGvdW_Center);

    }

    if((SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ) ||
       (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_EPSILON )){
      //Xinfeng : differentiate diamond and average process
      dFvdW_Turbulent(dFvdW_Center, SolnBlk, ii,jj);
      dGvdW_Turbulent(dGvdW_Center, SolnBlk, ii,jj);
      //Viscous flux Jacobians at North face of cell (ii, jj)
      dFvdW_N.zero();
      dGvdW_N.zero();
      Orient =  "NORTH"; 
      //  cell (ii, jj+1)
      dFvdW_Turbulent_Neigbour( dFvdW_N, SolnBlk,Orient,ii,jj);  
      dGvdW_Turbulent_Neigbour( dGvdW_N, SolnBlk,Orient,ii,jj);
      dFvdW_N = HALF*(dFvdW_N+dFvdW_Center);
      dGvdW_N = HALF*(dGvdW_N+dGvdW_Center);
      
      // Viscous flux Jacobians at South face of cell (ii, jj)
      dFvdW_S.zero();
      dGvdW_S.zero();
      Orient =  "SOUTH";
      dFvdW_Turbulent_Neigbour( dFvdW_S, SolnBlk,Orient,ii,jj);
      dGvdW_Turbulent_Neigbour( dGvdW_S, SolnBlk,Orient,ii,jj);
      dFvdW_S = HALF*(dFvdW_S+dFvdW_Center);
      dGvdW_S = HALF*(dGvdW_S+dGvdW_Center);
    
      // Viscous flux Jacobians at West face of cell (ii, jj)
      dFvdW_W.zero();
      dGvdW_W.zero();
      Orient =  "WEST";
      dFvdW_Turbulent_Neigbour( dFvdW_W, SolnBlk,Orient,ii,jj);
      dGvdW_Turbulent_Neigbour( dGvdW_W, SolnBlk,Orient,ii,jj);
      dFvdW_W = HALF*(dFvdW_W+dFvdW_Center);
      dGvdW_W = HALF*(dGvdW_W+dGvdW_Center);
      // Viscous flux Jacobians at East face of cell (ii, jj)
      dFvdW_E.zero();
      dGvdW_E.zero();
      Orient =  "EAST";
      dFvdW_Turbulent_Neigbour( dFvdW_E, SolnBlk,Orient,ii,jj);
      dGvdW_Turbulent_Neigbour( dGvdW_E, SolnBlk,Orient,ii,jj);
      dFvdW_E = HALF*(dFvdW_E+dFvdW_Center);
      dGvdW_E = HALF*(dGvdW_E+dGvdW_Center);
      
    }//end of turbulent case
  }//end of the viscous formualtion not diamond path
  else{//Diamond path
    
    if(SolnBlk.Flow_Type == FLOWTYPE_LAMINAR){
      
      dFvdW_N.zero();
      dGvdW_N.zero();
      Orient =  "NORTH"; 
      dFvdW_Laminar_Diamond(dFvdW_N, dWf_dWc, SolnBlk, Orient,  ii,jj);
      dGvdW_Laminar_Diamond(dGvdW_N, dWf_dWc,SolnBlk, Orient, ii,jj);
    
      // Viscous flux Jacobians at South face of cell (ii, jj)
      dFvdW_S.zero();
      dGvdW_S.zero();
      Orient =  "SOUTH";
      dFvdW_Laminar_Diamond(dFvdW_S, dWf_dWc,SolnBlk, Orient,  ii,jj);
      dGvdW_Laminar_Diamond(dGvdW_S, dWf_dWc,SolnBlk, Orient, ii,jj);
      
      // Viscous flux Jacobians at West face of cell (ii, jj)
      dFvdW_W.zero();
      dGvdW_W.zero();
      Orient =  "WEST";
      dFvdW_Laminar_Diamond(dFvdW_W, dWf_dWc,SolnBlk, Orient,  ii,jj);
      dGvdW_Laminar_Diamond(dGvdW_W, dWf_dWc,SolnBlk, Orient, ii,jj);
 
      // Viscous flux Jacobians at East face of cell (ii, jj)
      dFvdW_E.zero();
      dGvdW_E.zero();
      Orient =  "EAST";
      dFvdW_Laminar_Diamond(dFvdW_E, dWf_dWc,SolnBlk, Orient,  ii,jj);
      dGvdW_Laminar_Diamond(dGvdW_E, dWf_dWc,SolnBlk, Orient, ii,jj);
 
      
    }

    if((SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ) ||
       (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_EPSILON )){
   
      //Viscous flux Jacobians at North face of cell (ii, jj)
      dFvdW_N.zero();
      dGvdW_N.zero();
      Orient =  "NORTH"; 
      dFvdW_Turbulent_Diamond(dFvdW_N, dWf_dWc,SolnBlk, Orient, ii,jj);
      dGvdW_Turbulent_Diamond(dGvdW_N, dWf_dWc,SolnBlk, Orient,ii,jj);
      // Viscous flux Jacobians at South face of cell (ii, jj)
      dFvdW_S.zero();
      dGvdW_S.zero();
      Orient =  "SOUTH";
      dFvdW_Turbulent_Diamond(dFvdW_S, dWf_dWc,SolnBlk, Orient, ii,jj);
      dGvdW_Turbulent_Diamond(dGvdW_S, dWf_dWc,SolnBlk, Orient,ii,jj);
     
      // Viscous flux Jacobians at West face of cell (ii, jj)
      dFvdW_W.zero();
      dGvdW_W.zero();
      Orient =  "WEST";
      dFvdW_Turbulent_Diamond(dFvdW_W, dWf_dWc,SolnBlk,Orient, ii,jj);
      dGvdW_Turbulent_Diamond(dGvdW_W, dWf_dWc,SolnBlk,Orient,ii,jj);
      // Viscous flux Jacobians at East face of cell (ii, jj)
      dFvdW_E.zero();
      dGvdW_E.zero();
      Orient =  "EAST";
      dFvdW_Turbulent_Diamond(dFvdW_E, dWf_dWc,SolnBlk, Orient, ii,jj);
      dGvdW_Turbulent_Diamond(dGvdW_E, dWf_dWc,SolnBlk, Orient,ii,jj);
         
    }
  }//DIAMOND PATH
  
  dGVdW = 1.0/ SolnBlk.Grid.Cell[ii][jj].A*(SolnBlk.Grid.lfaceN(ii, jj)* (nface_N.x* dFvdW_N + nface_N.y* dGvdW_N)
					    + SolnBlk.Grid.lfaceS(ii, jj)*(nface_S.x* dFvdW_S + nface_S.y* dGvdW_S)
					    + SolnBlk.Grid.lfaceW(ii, jj)*(nface_W.x* dFvdW_W + nface_W.y* dGvdW_W)
					    + SolnBlk.Grid.lfaceE(ii, jj)*(nface_E.x* dFvdW_E + nface_E.y* dGvdW_E)
					    );
  //  cout<< "This function is for computing the viscous components in point Implicit Block Jacobi"<<endl;
  return (dGVdW);
  }
  
}


//Laminar viscous Jacobian --Diamond path
void dFvdW_Laminar_Diamond(DenseMatrix &dFvdW, DenseMatrix &dWf_dWc, Chem2D_Quad_Block &SolnBlk,
				  const string &Orient,  const int &ii, const int &jj){
  //planar flow
  double  d_dWdx_dW=0;
  double  d_dWdy_dW=0;
    
  double t1,t2,t5,t8,t9,t11,t13, t16,t18;
  double t20,t21, t22, t27, t35,t38,t46;
    
  double kappa, Cp, mu;
  double rho, U, V;

  double  *Dm, *h, *dcdx, *dhdT;
  double dUdx,dUdy, dVdx,dVdy;
  double dmudT;
  int ns = SolnBlk.W[0][0].ns;
  Chem2D_pState  QuadraturePoint_W;
  Dm = new double [ns];
  h  = new double [ns];
  dcdx = new double [ns];
  dhdT = new double [ns];
  //needs to deallocate the momery
  if (Orient == "NORTH"){ 
    
    QuadraturePoint_W = HALF*(SolnBlk.WnNW(ii, jj) +SolnBlk.WnNE(ii, jj) );
    Cp = HALF*(SolnBlk.WnNW_Cp(ii, jj) +SolnBlk.WnNE_Cp(ii, jj));
    kappa = HALF*(SolnBlk.WnNW_kappa(ii, jj) +SolnBlk.WnNE_kappa(ii, jj)); 
    mu = HALF*(SolnBlk.WnNW_mu(ii, jj) +SolnBlk.WnNE_mu(ii, jj)); 
    dmudT = HALF*(SolnBlk.WnNW_dmudT(ii, jj) +SolnBlk.WnNE_dmudT(ii, jj)); 

    for(int Num = 0; Num<ns; Num++){
      h[Num] = HALF*(SolnBlk.WnNW_hi(ii, jj, Num) +SolnBlk.WnNE_hi(ii, jj, Num));
      Dm[Num] = HALF*(SolnBlk.WnNW_Dm(ii, jj, Num) +SolnBlk.WnNE_Dm(ii, jj, Num));
      dcdx[Num] = SolnBlk.dWdx_faceN[ii][jj].spec[Num].c;
      dhdT[Num] = HALF*(SolnBlk.WnNW_Cpi(ii, jj, Num) +SolnBlk.WnNE_Cpi(ii, jj, Num));
    }
    
    dUdx = SolnBlk.dWdx_faceN[ii][jj].v.x;
    dUdy = SolnBlk.dWdy_faceN[ii][jj].v.x;
    dVdx = SolnBlk.dWdx_faceN[ii][jj].v.y;
    dVdy = SolnBlk.dWdy_faceN[ii][jj].v.y;
    
    d_dWdx_dW =  SolnBlk.d_dWdx_dW[ii][jj][1];
    d_dWdy_dW =  SolnBlk.d_dWdy_dW[ii][jj][1];
  
    dWfdWc(SolnBlk, dWf_dWc, Orient,  ii, jj);
  
  }
  
  if (Orient == "EAST"){ 
   
    QuadraturePoint_W = HALF*(SolnBlk.WnNE(ii, jj) +SolnBlk.WnSE(ii, jj) );
    Cp = HALF*(SolnBlk.WnNE_Cp(ii, jj) +SolnBlk.WnSE_Cp(ii, jj));
    kappa = HALF*(SolnBlk.WnNE_kappa(ii, jj) +SolnBlk.WnSE_kappa(ii, jj)); 
    mu = HALF*(SolnBlk.WnNE_mu(ii, jj) +SolnBlk.WnSE_mu(ii, jj)); 
    dmudT = HALF*(SolnBlk.WnNE_dmudT(ii, jj) +SolnBlk.WnSE_dmudT(ii, jj)); 
       
    for(int Num = 0; Num<ns; Num++){
      h[Num] = HALF*(SolnBlk.WnNE_hi(ii, jj, Num) +SolnBlk.WnSE_hi(ii, jj, Num));
      Dm[Num] = HALF*(SolnBlk.WnNE_Dm(ii, jj, Num) +SolnBlk.WnSE_Dm(ii, jj, Num));
      dcdx[Num] = SolnBlk.dWdx_faceE[ii][jj].spec[Num].c;
      dhdT[Num] = HALF*(SolnBlk.WnNE_Cpi(ii, jj, Num) +SolnBlk.WnSE_Cpi(ii, jj, Num));
    }

    dUdx = SolnBlk.dWdx_faceE[ii][jj].v.x;
    dUdy = SolnBlk.dWdy_faceE[ii][jj].v.x;
    dVdx = SolnBlk.dWdx_faceE[ii][jj].v.y;
    dVdy = SolnBlk.dWdy_faceE[ii][jj].v.y;
    
    d_dWdx_dW =  SolnBlk.d_dWdx_dW[ii][jj][2];
    d_dWdy_dW =  SolnBlk.d_dWdy_dW[ii][jj][2];
    dWfdWc(SolnBlk, dWf_dWc, Orient,  ii, jj);
  
  }
    if (Orient == "SOUTH"){ 
      
      QuadraturePoint_W = HALF*(SolnBlk.WnSE(ii, jj) +SolnBlk.WnSW(ii, jj) );
      Cp = HALF*(SolnBlk.WnSW_Cp(ii, jj) +SolnBlk.WnSE_Cp(ii, jj));
      kappa = HALF*(SolnBlk.WnSW_kappa(ii, jj) +SolnBlk.WnSE_kappa(ii, jj)); 
      mu = HALF*(SolnBlk.WnSW_mu(ii, jj) +SolnBlk.WnSE_mu(ii, jj)); 
      dmudT = HALF*(SolnBlk.WnSW_dmudT(ii, jj) +SolnBlk.WnSE_dmudT(ii, jj)); 

      for(int Num = 0; Num<ns; Num++){
	h[Num] = HALF*(SolnBlk.WnSW_hi(ii, jj, Num) +SolnBlk.WnSE_hi(ii, jj, Num));
	Dm[Num] = HALF*(SolnBlk.WnSW_Dm(ii, jj, Num) +SolnBlk.WnSE_Dm(ii, jj, Num));
	dcdx[Num] = SolnBlk.dWdx_faceS[ii][jj].spec[Num].c;
	dhdT[Num] = HALF*(SolnBlk.WnSW_Cpi(ii, jj, Num) +SolnBlk.WnSE_Cpi(ii, jj, Num));
	
      }
      
      dUdx = SolnBlk.dWdx_faceS[ii][jj].v.x;
      dUdy = SolnBlk.dWdy_faceS[ii][jj].v.x;
      dVdx = SolnBlk.dWdx_faceS[ii][jj].v.y;
      dVdy = SolnBlk.dWdy_faceS[ii][jj].v.y;
      
      d_dWdx_dW =  SolnBlk.d_dWdx_dW[ii][jj][4];
      d_dWdy_dW =  SolnBlk.d_dWdy_dW[ii][jj][4];
      dWfdWc(SolnBlk, dWf_dWc, Orient,  ii, jj);
    }

    if (Orient == "WEST"){ 
      
      QuadraturePoint_W = HALF*(SolnBlk.WnNW(ii, jj) +SolnBlk.WnSW(ii, jj));
      Cp = HALF*(SolnBlk.WnNW_Cp(ii, jj) +SolnBlk.WnSW_Cp(ii, jj));
      kappa = HALF*(SolnBlk.WnNW_kappa(ii, jj) +SolnBlk.WnSW_kappa(ii, jj)); 
      mu = HALF*(SolnBlk.WnNW_mu(ii, jj) +SolnBlk.WnSW_mu(ii, jj)); 
      dmudT = HALF*(SolnBlk.WnNW_dmudT(ii, jj) +SolnBlk.WnSW_dmudT(ii, jj)); 
    
      for(int Num = 0; Num<ns; Num++){
	h[Num] = HALF*(SolnBlk.WnNW_hi(ii, jj, Num) +SolnBlk.WnSW_hi(ii, jj, Num));
	Dm[Num] = HALF*(SolnBlk.WnNW_Dm(ii, jj, Num) +SolnBlk.WnSW_Dm(ii, jj, Num));
	dcdx[Num] = SolnBlk.dWdx_faceW[ii][jj].spec[Num].c;
	dhdT[Num] = HALF*(SolnBlk.WnNW_Cpi(ii, jj, Num) +SolnBlk.WnSW_Cpi(ii, jj, Num));

      }

      dUdx = SolnBlk.dWdx_faceW[ii][jj].v.x;
      dUdy = SolnBlk.dWdy_faceW[ii][jj].v.x;
      dVdx = SolnBlk.dWdx_faceW[ii][jj].v.y;
      dVdy = SolnBlk.dWdy_faceW[ii][jj].v.y;
      
      d_dWdx_dW =  SolnBlk.d_dWdx_dW[ii][jj][3];
      d_dWdy_dW =  SolnBlk.d_dWdy_dW[ii][jj][3];
      dWfdWc(SolnBlk, dWf_dWc, Orient,  ii, jj);
      
    }
 
    rho = QuadraturePoint_W.rho;
    U = QuadraturePoint_W.v.x;
    V = QuadraturePoint_W.v.y;

    t1 = mu;
    t2 = d_dWdx_dW;
    t5 = d_dWdy_dW;
    t8 = dmudT;
    t9 = dUdx;
    t11 = dVdy;
    t13 = 2.0/3.0*t9-t11/3.0;
    t16 = d_dWdy_dW;
    t18 = d_dWdx_dW;
    t20 = dUdy;
    t21 = dVdx;
    t22 = t20+t21;
    t35 = U*t1;
    t38 = V*t1;
    t46 = d_dWdx_dW;
    
    dFvdW(1,1) = 4.0/3.0*t1*t2;
    dFvdW(1,2) = -2.0/3.0*t1*t5;
    dFvdW(1,3) = 2.0*t8*t13;
    
    dFvdW(2,1) = t1*t16;
    dFvdW(2,2) = t1*t18;
    dFvdW(2,3) = t8*t22;

    double Sum_q = 0.0;
    double Sum_dq = 0.0;
    for(int Num = 0; Num<ns; Num++)
      {
	//for each species	// h*Dm*gradc  
	Sum_q +=  h[Num]*Dm[Num]*dcdx[Num];
	//dhdT *(Dm)*gradc
	Sum_dq +=   dhdT[Num]*Dm[Num]*dcdx[Num];
	
      }  
    dFvdW(3,0) = Sum_q ;
    dFvdW(3,1) = 2.0*t1*t13+4.0/3.0*t35*t2+t38*t16;
    dFvdW(3,2) = -2.0/3.0*t35*t5+t1*t22+t38*t18;
    dFvdW(3,3) = kappa*t46+rho*Sum_dq+2.0*U*t8*t13+V*t8*t22;
      
  //multispecies
  int NUM_VAR = NUM_CHEM2D_VAR_SANS_SPECIES; 
  for(int Num = 0; Num<(ns-1); Num++)
    {//h*rho*(Dm)*d_dWdx_dW
      dFvdW(3, NUM_VAR+Num) = h[Num]*rho*Dm[Num]*d_dWdx_dW;
      //(Dm)*gradc
      dFvdW(NUM_VAR+Num, 0) = Dm[Num]*dcdx[Num];
      //rho*(Dm)*dcxdc
      dFvdW(NUM_VAR+Num,NUM_VAR+Num) = rho*Dm[Num]*d_dWdx_dW;
    }

  //The following entries are different for axisymmetric flow and planar flow
  // For those terms involves fluid molecular stess tau_xx and tau_yy
  // They are some terms in the scond row  and fourth of Jacobian matrix
  if(SolnBlk.Axisymmetric){ 
    double r;
    if (SolnBlk.Axisymmetric == 1) {
      r = SolnBlk.Grid.Cell[ii][jj].Xc.y;
      
      t13 -= (V/r)/3.0;
    
      dFvdW(1,2) = 2.0/3.0*t1*(-t5-ONE/r);
      dFvdW(1,3) = 2.0*t8*t13;
      dFvdW(3,2) = 2.0/3.0*t35*(-t5-ONE/r)+t1*t22+t38*t18;
      dFvdW(3,3) = kappa*t46+rho*Sum_dq +2.0*U*t8*t13+V*t8*t22;
    }
    if (SolnBlk.Axisymmetric == 2) {
      r = SolnBlk.Grid.Cell[ii][jj].Xc.x;
      double term1 = 2.0/3.0*t2 - 1.0/r*1.0/3.0;
      
      t13 -= (U/r)/3.0;

      dFvdW(1,1) = 4.0/3.0*t1*term1;
      dFvdW(1,3) = 2.0*t8*t13;
      dFvdW(3,1) = 2.0*t1*t13+2.0*t35*term1+t38*t16;
      dFvdW(3,3) = kappa*t46+rho*Sum_dq +2.0*U*t8*t13+V*t8*t22;
    }
  }//end of axisymmetric case
  // Here the use of dFvdW is only convenient, but
  // dGvdW is the Jacobian obtained (with respect to the primitive solution variable on the cell edges )
  // to obtain the Jacobian (with respect to the primitive soltuion variable at the cell center)
  // dGvdW *dWfdWc 
 
 
  dFvdW = dFvdW * dWf_dWc;
  

}//Laminar viscous Jacobian (X)

void dGvdW_Laminar_Diamond(DenseMatrix &dGvdW,  DenseMatrix &dWf_dWc, Chem2D_Quad_Block &SolnBlk,
				  const string &Orient,
				  const int &ii, const int &jj){

  //planar flow
  double  d_dWdx_dW=0;
  double  d_dWdy_dW=0;
  
  double t1,t2,t5,t8,t9,t11,t13, t16,t18;
  double t20,t21, t22, t27, t35,t38,t46;
  
  double kappa, mu;
  double rho, U, V;
  double dmudT; 
  double  *Dm, *h, *dcdy, *dhdT;
  double dUdx,dUdy, dVdx,dVdy;
  int ns = SolnBlk.W[0][0].ns;
  Chem2D_pState  QuadraturePoint_W;
  Dm = new double [ns];
  h  = new double [ns];
  dcdy = new double [ns];
  dhdT = new double [ns];
  //needs to deallocate the momery
 
  if (Orient == "NORTH"){ 
     
    QuadraturePoint_W = HALF*(SolnBlk.WnNW(ii, jj) +SolnBlk.WnNE(ii, jj) );
    kappa = HALF*(SolnBlk.WnNW_kappa(ii, jj) +SolnBlk.WnNE_kappa(ii, jj)); 
    mu = HALF*(SolnBlk.WnNW_mu(ii, jj) +SolnBlk.WnNE_mu(ii, jj)); 
    dmudT = HALF*(SolnBlk.WnNW_dmudT(ii, jj) +SolnBlk.WnNE_dmudT(ii, jj)); 
    for(int Num = 0; Num<ns; Num++){
      h[Num] = HALF*(SolnBlk.WnNW_hi(ii, jj, Num) +SolnBlk.WnNE_hi(ii, jj, Num));
      Dm[Num] = HALF*(SolnBlk.WnNW_Dm(ii, jj, Num) +SolnBlk.WnNE_Dm(ii, jj, Num));
      dcdy[Num] = SolnBlk.dWdy_faceN[ii][jj].spec[Num].c;
      dhdT[Num] = HALF*(SolnBlk.WnNW_Cpi(ii, jj, Num) +SolnBlk.WnNE_Cpi(ii, jj, Num));
    }
    
    dUdx = SolnBlk.dWdx_faceN[ii][jj].v.x;
    dUdy = SolnBlk.dWdy_faceN[ii][jj].v.x;
    dVdx = SolnBlk.dWdx_faceN[ii][jj].v.y;
    dVdy = SolnBlk.dWdy_faceN[ii][jj].v.y;
    
    d_dWdx_dW =  SolnBlk.d_dWdx_dW[ii][jj][1];
    d_dWdy_dW =  SolnBlk.d_dWdy_dW[ii][jj][1];
    dWfdWc(SolnBlk, dWf_dWc, Orient,  ii, jj);
     
    }
    
  if (Orient == "EAST"){ 
 
    QuadraturePoint_W = HALF*(SolnBlk.WnNE(ii, jj) +SolnBlk.WnSE(ii, jj) );
    kappa = HALF*(SolnBlk.WnNE_kappa(ii, jj) +SolnBlk.WnSE_kappa(ii, jj)); 
    mu = HALF*(SolnBlk.WnNE_mu(ii, jj) +SolnBlk.WnSE_mu(ii, jj)); 
    dmudT = HALF*(SolnBlk.WnNE_dmudT(ii, jj) +SolnBlk.WnSE_dmudT(ii, jj)); 
    
    for(int Num = 0; Num<ns; Num++){
      h[Num] = HALF*(SolnBlk.WnNE_hi(ii, jj, Num) +SolnBlk.WnSE_hi(ii, jj, Num));
      Dm[Num] = HALF*(SolnBlk.WnNE_Dm(ii, jj, Num) +SolnBlk.WnSE_Dm(ii, jj, Num));
      dcdy[Num] = SolnBlk.dWdy_faceE[ii][jj].spec[Num].c;
      dhdT[Num] = HALF*(SolnBlk.WnNE_Cpi(ii, jj, Num) +SolnBlk.WnSE_Cpi(ii, jj, Num));
    }
    
    dUdx = SolnBlk.dWdx_faceE[ii][jj].v.x;
    dUdy = SolnBlk.dWdy_faceE[ii][jj].v.x;
    dVdx = SolnBlk.dWdx_faceE[ii][jj].v.y;
    dVdy = SolnBlk.dWdy_faceE[ii][jj].v.y;
    
    d_dWdx_dW =  SolnBlk.d_dWdx_dW[ii][jj][2];
    d_dWdy_dW =  SolnBlk.d_dWdy_dW[ii][jj][2];
  
    dWfdWc(SolnBlk, dWf_dWc, Orient,  ii, jj);
       
  }
  if (Orient == "SOUTH"){ 
    
    QuadraturePoint_W = HALF*(SolnBlk.WnSE(ii, jj) +SolnBlk.WnSW(ii, jj) );
    kappa = HALF*(SolnBlk.WnSW_kappa(ii, jj) +SolnBlk.WnSE_kappa(ii, jj)); 
    mu = HALF*(SolnBlk.WnSW_mu(ii, jj) +SolnBlk.WnSE_mu(ii, jj)); 
    dmudT = HALF*(SolnBlk.WnSW_dmudT(ii, jj) +SolnBlk.WnSE_dmudT(ii, jj)); 
    
    for(int Num = 0; Num<ns; Num++){
      h[Num] = HALF*(SolnBlk.WnSW_hi(ii, jj, Num) +SolnBlk.WnSE_hi(ii, jj, Num));
      Dm[Num] = HALF*(SolnBlk.WnSW_Dm(ii, jj, Num) +SolnBlk.WnSE_Dm(ii, jj, Num));
      dcdy[Num] = SolnBlk.dWdy_faceS[ii][jj].spec[Num].c;
      dhdT[Num] = HALF*(SolnBlk.WnSW_Cpi(ii, jj, Num) +SolnBlk.WnSE_Cpi(ii, jj, Num));
      
    }
    
    dUdx = SolnBlk.dWdx_faceS[ii][jj].v.x;
    dUdy = SolnBlk.dWdy_faceS[ii][jj].v.x;
    dVdx = SolnBlk.dWdx_faceS[ii][jj].v.y;
    dVdy = SolnBlk.dWdy_faceS[ii][jj].v.y;
   
    d_dWdx_dW =  SolnBlk.d_dWdx_dW[ii][jj][4];
    d_dWdy_dW =  SolnBlk.d_dWdy_dW[ii][jj][4];
   
    dWfdWc(SolnBlk, dWf_dWc, Orient,  ii, jj);
  }

  if (Orient == "WEST"){ 
    
    QuadraturePoint_W = HALF*(SolnBlk.WnNW(ii, jj) +SolnBlk.WnSW(ii, jj));
    kappa = HALF*(SolnBlk.WnNW_kappa(ii, jj) +SolnBlk.WnSW_kappa(ii, jj)); 
    mu = HALF*(SolnBlk.WnNW_mu(ii, jj) +SolnBlk.WnSW_mu(ii, jj)); 
    dmudT = HALF*(SolnBlk.WnNW_dmudT(ii, jj) +SolnBlk.WnSW_dmudT(ii, jj)); 
      
    for(int Num = 0; Num<ns; Num++){
      h[Num] = HALF*(SolnBlk.WnNW_hi(ii, jj, Num) +SolnBlk.WnSW_hi(ii, jj, Num));
      Dm[Num] = HALF*(SolnBlk.WnNW_Dm(ii, jj, Num) +SolnBlk.WnSW_Dm(ii, jj, Num));
      dcdy[Num] = SolnBlk.dWdy_faceW[ii][jj].spec[Num].c;
      dhdT[Num] = HALF*(SolnBlk.WnNW_Cpi(ii, jj, Num) +SolnBlk.WnSW_Cpi(ii, jj, Num));

    }
    dUdx = SolnBlk.dWdx_faceW[ii][jj].v.x;
    dUdy = SolnBlk.dWdy_faceW[ii][jj].v.x;
    dVdx = SolnBlk.dWdx_faceW[ii][jj].v.y;
    dVdy = SolnBlk.dWdy_faceW[ii][jj].v.y;
    d_dWdx_dW =  SolnBlk.d_dWdx_dW[ii][jj][3];
    d_dWdy_dW =  SolnBlk.d_dWdy_dW[ii][jj][3];
    dWfdWc(SolnBlk, dWf_dWc, Orient,  ii, jj);
    
  }
  rho = QuadraturePoint_W.rho;
  U = QuadraturePoint_W.v.x;
  V = QuadraturePoint_W.v.y;
  t1 = mu;
  t2 = d_dWdy_dW;
  t5 = d_dWdx_dW;
  t8 = dmudT;
  t9 = dUdy;
  t11 = dVdx;
  t13 = 2.0/3.0*t9-t11/3.0;
  t16 = d_dWdx_dW;
  t18 = d_dWdy_dW;
  t20 = dVdy;
  t21 = dUdx;
  t22 = t20+t21;
  t35 = U*t1;
  t38 = V*t1;
  t46 = d_dWdy_dW;
  
  dGvdW(1,1) = t1*t2;
  dGvdW(1,2) = -2.0/3.0*t1*t5;
  dGvdW(1,3) = t8*t22;
    
  dGvdW(2,1) = -2.0/3.0*t1*t16;
  dGvdW(2,2) = 4.0/3.0*t1*t18;
  dGvdW(2,3) = 2.0*t8*t13;

  double Sum_q = 0.0;
  double Sum_dq = 0.0;
 
    
  for(int Num = 0; Num<ns; Num++)
    {
      //for each species	// h*Dm*gradc  
      Sum_q +=  h[Num]*Dm[Num]*dcdy[Num];
      //dhdT *(Dm)*gradc
      Sum_dq += dhdT[Num]*Dm[Num]*dcdy[Num];
      
    }  
  dGvdW(3,0) = Sum_q ;
  dGvdW(3,1) = t1*t22+t35*t2-2.0/3.0*t38*t16;
  dGvdW(3,2) = t35*t5+2.0*t1*t13+4.0/3.0*t38*t18;
  dGvdW(3,3) = kappa*t46+rho*Sum_dq+U*t8*t22+2.0*V*t8*t13;
      
  //multispecies
  int NUM_VAR =NUM_CHEM2D_VAR_SANS_SPECIES; 
  for(int Num = 0; Num<(ns-1); Num++)
    {//h*rho*(Dm)*d_dWdy_dW
      dGvdW(3, NUM_VAR+Num) = h[Num]*Dm[Num]*rho*d_dWdy_dW;
      //(Dm)*gradc
      dGvdW(NUM_VAR+Num, 0) = Dm[Num]*dcdy[Num];
      //rho*(Dm)*dcydc
      dGvdW(NUM_VAR+Num,NUM_VAR+Num) = rho*Dm[Num]*d_dWdy_dW;
    }

  //The following entries are different for axisymmetric flow and planar flow
  // For those terms involves fluid molecular stess tau_xx and tau_yy
  // They are some terms in the third row  and fourth of Jacobian matrix
  if(SolnBlk.Axisymmetric){ 
    double r;
    if (SolnBlk.Axisymmetric == 1) {
      r = SolnBlk.Grid.Cell[ii][jj].Xc.y;
      
      t13 -= (V/r)/3.0;
      t18 = 2.0/3.0*t18 -1.0/r*1.0/3.0;
    
      dGvdW(2,2) = 2.0*t1*t18;
      dGvdW(2,3) = 2.0*t8*t13;
      dGvdW(3,2) = t35*t5+2.0*t1*t13+2.0*t38*t18;
      dGvdW(3,3) = kappa*t46+rho*Sum_dq +2.0*U*t8*t22+2.0*V*t8*t13;
    }
    if (SolnBlk.Axisymmetric == 2) {
      r = SolnBlk.Grid.Cell[ii][jj].Xc.x;
      double term1 = -t16 - 1.0/r;
      
      t13 -= (U/r)/3.0;

      dGvdW(2,1) = -2.0/3.0*t1*term1;
      dGvdW(2,3) = 2.0*t8*t13;
      dGvdW(3,1) = t1*t22+t35*t2+2.0/3.0*t38*t13;
      dGvdW(3,3) = kappa*t46+rho*Sum_dq +U*t8*t22+2.0*V*t8*t13;
    }
  }//end of axisymmetric case
  dGvdW = dGvdW * dWf_dWc;


}//Laminar viscous Jacobian (Y)
//Turbulent Flux Jacobian
void dFvdW_Turbulent_Diamond(DenseMatrix &dFvdW,  DenseMatrix &dWf_dWc, Chem2D_Quad_Block &SolnBlk, 
				    const string &Orient, const int &ii, const int &jj){
 

  double  d_dWdx_dW=0;
  double  d_dWdy_dW=0;
  
  double t1,t2,t3,t5,t7,t11,t12,t13, t15,t18,t19,t23, t24, t27, t31,t32;
  double t33,t37,t38,t39,t41,t45,t46,t50, t56, t57;
  double t85,t86, t87, t90, t110, t140;
  double t149, t150, t153, t154, t156, t160, t164;
  double t194, t195, t198, t204;
 
  double kappa, Cp, mu, sigma, sigma_star;
  double rho, U, V, k, omega;
  double mu_turb, kappa_turb, Pr_turb, Dm_turb;
  double Temp =SolnBlk.W[ii][jj].T();

  double dmudT; 
  double  *Dm, *h, *dcdx, *dhdT;
  double dUdx,dUdy, dVdx,dVdy;
  double dkdx, dkdy, domegadx, domegady;
  double drhodx;
  int ns = SolnBlk.W[0][0].ns;
  double Rmix = SolnBlk.W[ii][jj].Rtot();    
  Chem2D_pState  QuadraturePoint_W;
  Dm = new double [ns];
  h  = new double [ns];
  dcdx = new double [ns];
  dhdT = new double [ns];

  Pr_turb = SolnBlk.W[ii][jj].Pr_turb(); 
  Dm_turb = SolnBlk.W[ii][jj].Dm_turb(); 
  
  //needs to deallocate the momery
 
  if (Orient == "NORTH"){ 
     
    QuadraturePoint_W = HALF*(SolnBlk.WnNW(ii, jj) +SolnBlk.WnNE(ii, jj) );
    kappa = HALF*(SolnBlk.WnNW_kappa(ii, jj) +SolnBlk.WnNE_kappa(ii, jj)); 
    mu = HALF*(SolnBlk.WnNW_mu(ii, jj) +SolnBlk.WnNE_mu(ii, jj)); 
    dmudT = HALF*(SolnBlk.WnNW_dmudT(ii, jj) +SolnBlk.WnNE_dmudT(ii, jj)); 
 
    for(int Num = 0; Num<ns; Num++){
      h[Num] = HALF*(SolnBlk.WnNW_hi(ii, jj, Num) +SolnBlk.WnNE_hi(ii, jj, Num));
      Dm[Num] = HALF*(SolnBlk.WnNW_Dm(ii, jj, Num) +SolnBlk.WnNE_Dm(ii, jj, Num));
      dcdx[Num] = SolnBlk.dWdx_faceN[ii][jj].spec[Num].c;
      dhdT[Num] = HALF*(SolnBlk.WnNW_Cpi(ii, jj, Num) +SolnBlk.WnNE_Cpi(ii, jj, Num));
    }
    drhodx = SolnBlk.dWdx_faceN[ii][jj].rho;
    dUdx = SolnBlk.dWdx_faceN[ii][jj].v.x;
    dUdy = SolnBlk.dWdy_faceN[ii][jj].v.x;
    dVdx = SolnBlk.dWdx_faceN[ii][jj].v.y;
    dVdy = SolnBlk.dWdy_faceN[ii][jj].v.y;
    dkdx = SolnBlk.dWdx_faceN[ii][jj].k;
    dkdy = SolnBlk.dWdy_faceN[ii][jj].k;
    domegadx = SolnBlk.dWdx_faceN[ii][jj].omega;
    domegady = SolnBlk.dWdy_faceN[ii][jj].omega;
   
 
    d_dWdx_dW =  SolnBlk.d_dWdx_dW[ii][jj][1];
    d_dWdy_dW =  SolnBlk.d_dWdy_dW[ii][jj][1];

    dWfdWc(SolnBlk, dWf_dWc , Orient,  ii, jj);  
  
    }
    
  if (Orient == "EAST"){ 
 
    QuadraturePoint_W = HALF*(SolnBlk.WnNE(ii, jj) +SolnBlk.WnSE(ii, jj) );
    kappa = HALF*(SolnBlk.WnNE_kappa(ii, jj) +SolnBlk.WnSE_kappa(ii, jj)); 
    mu = HALF*(SolnBlk.WnNE_mu(ii, jj) +SolnBlk.WnSE_mu(ii, jj)); 
    dmudT = HALF*(SolnBlk.WnNE_dmudT(ii, jj) +SolnBlk.WnSE_dmudT(ii, jj)); 
  
    for(int Num = 0; Num<ns; Num++){
      h[Num] = HALF*(SolnBlk.WnNE_hi(ii, jj, Num) +SolnBlk.WnSE_hi(ii, jj, Num));
      Dm[Num] = HALF*(SolnBlk.WnNE_Dm(ii, jj, Num) +SolnBlk.WnSE_Dm(ii, jj, Num));
      dcdx[Num] = SolnBlk.dWdx_faceE[ii][jj].spec[Num].c;
      dhdT[Num] = HALF*(SolnBlk.WnNE_Cpi(ii, jj, Num) +SolnBlk.WnSE_Cpi(ii, jj, Num));
    }
    drhodx = SolnBlk.dWdx_faceE[ii][jj].rho;
    dUdx = SolnBlk.dWdx_faceE[ii][jj].v.x;
    dUdy = SolnBlk.dWdy_faceE[ii][jj].v.x;
    dVdx = SolnBlk.dWdx_faceE[ii][jj].v.y;
    dVdy = SolnBlk.dWdy_faceE[ii][jj].v.y;
    dkdx = SolnBlk.dWdx_faceE[ii][jj].k;
    dkdy = SolnBlk.dWdy_faceE[ii][jj].k;
    domegadx = SolnBlk.dWdx_faceE[ii][jj].omega;
    domegady = SolnBlk.dWdy_faceE[ii][jj].omega;
   
    d_dWdx_dW =  SolnBlk.d_dWdx_dW[ii][jj][2];
    d_dWdy_dW =  SolnBlk.d_dWdy_dW[ii][jj][2];

    dWfdWc(SolnBlk, dWf_dWc , Orient,  ii, jj);  
         
  }
  if (Orient == "SOUTH"){ 
    
    QuadraturePoint_W = HALF*(SolnBlk.WnSE(ii, jj) +SolnBlk.WnSW(ii, jj) );
    kappa = HALF*(SolnBlk.WnSW_kappa(ii, jj) +SolnBlk.WnSE_kappa(ii, jj)); 
    mu = HALF*(SolnBlk.WnSW_mu(ii, jj) +SolnBlk.WnSE_mu(ii, jj)); 
    dmudT = HALF*(SolnBlk.WnSW_dmudT(ii, jj) +SolnBlk.WnSE_dmudT(ii, jj)); 
    
    for(int Num = 0; Num<ns; Num++){
      h[Num] = HALF*(SolnBlk.WnSW_hi(ii, jj, Num) +SolnBlk.WnSE_hi(ii, jj, Num));
      Dm[Num] = HALF*(SolnBlk.WnSW_Dm(ii, jj, Num) +SolnBlk.WnSE_Dm(ii, jj, Num));
      dcdx[Num] = SolnBlk.dWdx_faceS[ii][jj].spec[Num].c;
      dhdT[Num] = HALF*(SolnBlk.WnSW_Cpi(ii, jj, Num) +SolnBlk.WnSE_Cpi(ii, jj, Num));
      
    }
    drhodx = SolnBlk.dWdx_faceS[ii][jj].rho;
    dUdx = SolnBlk.dWdx_faceS[ii][jj].v.x;
    dUdy = SolnBlk.dWdy_faceS[ii][jj].v.x;
    dVdx = SolnBlk.dWdx_faceS[ii][jj].v.y;
    dVdy = SolnBlk.dWdy_faceS[ii][jj].v.y;
    dkdx = SolnBlk.dWdx_faceS[ii][jj].k;
    dkdy = SolnBlk.dWdy_faceS[ii][jj].k;
    domegadx = SolnBlk.dWdx_faceS[ii][jj].omega;
    domegady = SolnBlk.dWdy_faceS[ii][jj].omega;

    d_dWdx_dW =  SolnBlk.d_dWdx_dW[ii][jj][4];
    d_dWdy_dW =  SolnBlk.d_dWdy_dW[ii][jj][4];
  
    dWfdWc(SolnBlk, dWf_dWc , Orient,  ii, jj); 
  
  }

  if (Orient == "WEST"){ 
    
    QuadraturePoint_W = HALF*(SolnBlk.WnNW(ii, jj) +SolnBlk.WnSW(ii, jj));
    kappa = HALF*(SolnBlk.WnNW_kappa(ii, jj) +SolnBlk.WnSW_kappa(ii, jj)); 
    mu = HALF*(SolnBlk.WnNW_mu(ii, jj) +SolnBlk.WnSW_mu(ii, jj)); 
    dmudT = HALF*(SolnBlk.WnNW_dmudT(ii, jj) +SolnBlk.WnSW_dmudT(ii, jj)); 
      
    for(int Num = 0; Num<ns; Num++){
      h[Num] = HALF*(SolnBlk.WnNW_hi(ii, jj, Num) +SolnBlk.WnSW_hi(ii, jj, Num));
      Dm[Num] = HALF*(SolnBlk.WnNW_Dm(ii, jj, Num) +SolnBlk.WnSW_Dm(ii, jj, Num));
      dcdx[Num] = SolnBlk.dWdx_faceW[ii][jj].spec[Num].c;
      dhdT[Num] = HALF*(SolnBlk.WnNW_Cpi(ii, jj, Num) +SolnBlk.WnSW_Cpi(ii, jj, Num));

    }
    drhodx = SolnBlk.dWdx_faceW[ii][jj].rho;
    dUdx = SolnBlk.dWdx_faceW[ii][jj].v.x;
    dUdy = SolnBlk.dWdy_faceW[ii][jj].v.x;
    dVdx = SolnBlk.dWdx_faceW[ii][jj].v.y;
    dVdy = SolnBlk.dWdy_faceW[ii][jj].v.y;
    dkdx = SolnBlk.dWdx_faceW[ii][jj].k;
    dkdy = SolnBlk.dWdy_faceW[ii][jj].k;
    domegadx = SolnBlk.dWdx_faceW[ii][jj].omega;
    domegady = SolnBlk.dWdy_faceW[ii][jj].omega;

    d_dWdx_dW =  SolnBlk.d_dWdx_dW[ii][jj][3];
    d_dWdy_dW =  SolnBlk.d_dWdy_dW[ii][jj][3];
    dWfdWc(SolnBlk, dWf_dWc , Orient,  ii, jj);  
  }
   
  mu_turb =  QuadraturePoint_W.rho* QuadraturePoint_W.k/max(QuadraturePoint_W.omega, TOLER);
  sigma = SolnBlk.W[0][0].sigma;
  sigma_star = SolnBlk.W[0][0].sigma_star;  
  
  rho = QuadraturePoint_W.rho;
  U = QuadraturePoint_W.v.x;
  V = QuadraturePoint_W.v.y;
  k= QuadraturePoint_W.k;
  omega= QuadraturePoint_W.omega;
  
  t1 = 1/max(omega, TOLER);
  t2 = k*t1;
  t3 = dUdx; 
  t5 = dVdy;  
  t7 = 2.0/3.0*t3-t5/3.0;
  t11 = 2.0*t2*t7-2.0/3.0*k;
  t12 = mu;
  t13 =  d_dWdx_dW;
  t15 = rho*k;
  t18 = t12*t13+t15*t1*t13;
  t19 =  d_dWdy_dW;
  t23 = -t12*t19-t15*t1*t19;
  t24 = dmudT;
  t27 = rho*t1;
  t31 = 2.0*t27*t7-2.0/3.0*rho;
  t32 = omega*omega;
  t33 = 1/t32;
  t37 = dUdy; 
  t38 = dVdx; 
  t39 = t37+t38;
  t41 = d_dWdx_dW;
  t45 = t12*t41+t15*t1*t41;
  t46 = d_dWdx_dW;
  t50 = t12*t46+t15*t1*t46;
  t56 = Cp/Pr_turb;
  //t57 = SolnBlk.dTdx(ii,jj);
  t57 = (ONE/(rho*Rmix)) * (QuadraturePoint_W.p - 
	(QuadraturePoint_W.p/rho)*drhodx);
  t85 = dkdx;
  t86 = t1*t85;
  t87 = sigma_star*k*t86;
  t90 = t1*t39;
  t110 = d_dWdx_dW;
  t140 = t24*t85;
  t149 = sigma_star*rho;
  t150 = t149*t86;
  t153 = d_dWdx_dW;
  t154 = (t12+t149*t2)*t153;
  t156 = V*rho;
  t160 = k*t33;
  t164 = t149*t160*t85;
  t194 =  domegadx;
  t195 = t1*t194;
  t198 = sigma*rho;
  t204 =  d_dWdx_dW;

  dFvdW(1,0) = t11;
  dFvdW(1,1) = 4.0/3.0*t18;
  dFvdW(1,2) = 2.0/3.0*t23;
  dFvdW(1,3) = 2.0*t24*t7;
  dFvdW(1,4) = t31;
  dFvdW(1,5) = -2.0*t15*t33*t7;
  dFvdW(2,0) = t2*t39;
  dFvdW(2,1) = t45;
  dFvdW(2,2) = t50;
  dFvdW(2,3) = t24*t39;
  dFvdW(2,4) = t27*t39;
  dFvdW(2,5) = -t15*t33*t39;
  
  double Sum_q = 0.0;
  double Sum_dq = 0.0;
  for(int Num = 0; Num<(ns); Num++)
    {
      //for each species	// h*(Dm+Dm_turb)*gradc  
      Sum_q +=  h[Num]*(Dm[Num]+Dm_turb)*dcdx[Num];
      //dhdT *(Dm)*gradc
      Sum_dq +=   dhdT[Num]*(Dm[Num]+Dm_turb)*dcdx[Num];
      
    }   

    
  dFvdW(3,0) = t56*t2*t57+ Sum_q+t87+U*t11+V*k*t90;
  dFvdW(3,1) = 2.0*t12*t7+2.0*t15*t1*t7-2.0/3.0*t15+4.0/3.0*U*t18+V*t45;
  dFvdW(3,2) = 2.0/3.0*U*t23+t12*t39+t15*t90+V*t50;
  dFvdW(3,3) = (kappa+t56*t15*t1)*t110+ rho*Sum_dq +t140+2.0*U*t24*t7+V*t24*t39;
  dFvdW(3,4) = t56*t27*t57+t150+t154+U*t31+t156*t90;
  dFvdW(3,5) = -t56*rho*t160*t57-t164-2.0*U*rho*t160*t7-t156*t160*t39;
   
   //multispecies
  int NUM_VAR = NUM_CHEM2D_VAR_SANS_SPECIES; 
  for(int Num = 0; Num<(ns-1); Num++)
    {//h*rho*(Dm+Dmt)*d_dWdx_dW
      dFvdW(3, NUM_VAR+Num) = h[Num]*rho*(Dm[Num]+Dm_turb)*d_dWdx_dW;
      //(Dm+Dmt)*gradc
      dFvdW(NUM_VAR+Num, 0) = (Dm[Num]+Dm_turb)*dcdx[Num];
      //rho*(Dm+Dmt)*dcxdc
      dFvdW(NUM_VAR+Num,NUM_VAR+Num) = rho*(Dm[Num]+Dm_turb)*d_dWdx_dW;
    }
  
  dFvdW(4,0) = t87;
  dFvdW(4,3) = t140;
  dFvdW(4,4) = t150+t154;
  dFvdW(4,5) = -t164;
  dFvdW(5,0) = sigma*k*t195;
  dFvdW(5,3) = t24*t194;
  dFvdW(5,4) = t198*t195;
  dFvdW(5,5) = -t198*t160*t194+(t12+t198*t2)*t204;
    
  //The following entries are different for axisymmetric flow and planar flow
  // For matrix row 1 and row 3
  if (SolnBlk.Axisymmetric == 1) {

    double term1,term2,term3,term5, term7,term10;
    double term14,term15,term16,term18;
    double term21,term22,term23,term27,term28,term31,term35,term36, term37, term41, term42, term43;
    double term45, term49,term50, term54, term60, term61, term89,term90, term91,term94, term114; 
    double term144, term153, term154, term157, term158, term160, term164,term168 , term198, term199, term202, term208; 
    
    double r = SolnBlk.Grid.Cell[ii][jj].Xc.y;
    term1 = 1/max(omega,TOLER);
    term2 = k*term1;
    term3 = SolnBlk.dWdx[ii][jj].v.x; 
    term5 =  SolnBlk.dWdy[ii][jj].v.y; 
    term7 = 1/r;
    term10 = 2.0/3.0*term3-term5/3.0-V*term7/3.0;
    term14 = 2.0*term2*term10-2.0/3.0*k;
    term15 = mu;
    term16 = d_dWdx_dW;
    term18 = rho*k;
    term21 = term15*term16+term18*term1*term16;
    term22 = d_dWdy_dW;
    term23 = -term22-term7;
    term27 = term15*term23/3.0+term18*term1*term23/3.0;
    term28 = dmudT;
    term31 = rho*term1;
    term35 = 2.0*term31*term10-2.0/3.0*rho;
    term36 = omega*omega;
    term37 = 1/term36;
    term41 = dUdy; 
    term42 = dVdx;;
    term43 = term41+term42;
    term45 = d_dWdy_dW;
    term49 = term15*term45+term18*term1*term45;
    term50 = d_dWdx_dW;
    term54 = term15*term50+term18*term1*term50;
    term60 = Cp/Pr_turb;
    //  term61 = SolnBlk.dTdx(ii,jj);
    term61 = (ONE/(rho*Rmix)) * ( QuadraturePoint_W.p - 
	     ( QuadraturePoint_W.p/rho)*drhodx);  
    term89 = dkdx;
    term90 = term1*term89;
    term91 = sigma*k*term90;
    term94 = term1*term43;
    term114 = d_dWdx_dW;
    term144 = term28*term89;
    term153 = sigma*rho;
    term154 = term153*term90;
    term157 = d_dWdx_dW;
    term158 = (term15+term153*term2)*term157;
    term160 = V*rho;
    term164 = k*term37;
    term168 = term153*term164*term89;
    term198 = domegadx;
    term199 = term1*term198;
    term202 = sigma_star*rho;
    term208 = d_dWdx_dW;
	
    dFvdW(1,0) = term14;
    dFvdW(1,1) = 4.0/3.0*term21;
    dFvdW(1,2) = 2.0*term27;
    dFvdW(1,3) = 2.0*term28*term10;
    dFvdW(1,4) = term35;
    dFvdW(1,5) = -2.0*term18*term37*term10;
    
    dFvdW(3,0) = term60*term2*term61+ Sum_q+term91+U*term14+V*k*term94;
    dFvdW(3,1) = 2.0*term15*term10+2.0*term18*term1*term10-2.0/3.0*term18+4.0/3.0*U*term21+V*term49;
    dFvdW(3,2) = 2.0*U*term27+term15*term43+term18*term94+V*term54;
    dFvdW(3,3) = (kappa+term60*term18*term1)*term114+ rho*Sum_dq+term144+2.0*U*term28*term10+V*term28*term43;
    dFvdW(3,4) = term60*term31*term61+term154+term158+U*term35+term160*term94;
    dFvdW(3,5) = -term60*rho*term164*term61-term168-2.0*U*rho*term164*term10-term160*term164*term43;
    
  }

  if (SolnBlk.Axisymmetric == 2) {
    double tt1,tt2,tt3,tt5,tt7,tt10,tt14,tt15,tt16;
    double tt19,tt21,tt24,tt25,tt29,tt30,tt33;
    double tt37,tt38,tt39,tt43,tt44,tt45,tt47,tt51, tt52;
    double tt56,tt62,tt63,tt79,tt80, tt81,tt84;
    double tt104,tt120,tt129,tt130,tt133,tt134, tt136,tt140;
    double tt144,tt152,tt155,tt157,tt160,tt164, tt165,tt168,tt174;
    
    double r = SolnBlk.Grid.Cell[ii][jj].Xc.x;
    
    tt1 = 1/max(omega,TOLER);
    tt2 = k*tt1;
    tt3 =  dUdx; 
    tt5 =  dVdy;
    tt7 = 1/r;
    tt10 = 2.0/3.0*tt3-tt5/3.0-U*tt7/3.0;
    tt14 = 2.0*tt2*tt10-2.0/3.0*k;
    tt15 = mu;
    tt16 = d_dWdx_dW;
    tt19 = 2.0/3.0*tt16-tt7/3.0;
    tt21 = rho*k;
    tt24 = tt15*tt19+tt21*tt1*tt19;
    tt25 = d_dWdy_dW;
    tt29 = -tt15*tt25-tt21*tt1*tt25;
    tt30 = dmudT;
    tt33 = rho*t1;
    tt37 = 2.0*tt33*tt10-2.0/3.0*rho;
    tt38 = omega*omega;
    tt39 = 1/tt38;
    tt43 = dUdy;
    tt44 = dVdx; 
    tt45 = tt43+tt44;
    tt47 = d_dWdy_dW;
    tt51 = tt15*tt47+tt21*tt1*tt47;
    tt52 = d_dWdx_dW;
    tt56 = tt15*tt52+tt21*tt1*tt52;
    tt62 = Cp/Pr_turb;
    tt79 = dkdx;
    tt80 = tt1*tt79;
    tt81 = sigma*k*tt80;
    tt84 = tt1*tt45;
    tt104 = d_dWdx_dW;
    tt120 = tt30*tt79;
    tt129 = sigma*rho;
    tt130 = tt129*tt80;
    tt133 = d_dWdx_dW;
    tt134 = (tt15+tt129*tt2)*tt133;
    tt136 = V*rho;
    tt140 = k*tt39;
    tt144 = tt129*tt140*tt79;
    tt164 = domegadx;
    tt165 = tt1*tt164;
    tt168 = sigma_star*rho;
    tt174 = d_dWdx_dW;
    
    dFvdW(2,0) =tt2*tt45;
    dFvdW(2,1) =tt51;
    dFvdW(2,2) =tt56;
    dFvdW(2,3) =tt30*tt45;
    dFvdW(2,4) =tt33*tt45;
    dFvdW(2,5) =-tt21*tt39*tt45;
    dFvdW(3,0) =tt62*tt2*tt63+Sum_q+tt81+U*tt14+V*k*tt84;
    dFvdW(3,1) =2.0*tt15*tt10+2.0*tt21*tt1*tt10-2.0/3.0*tt21+2.0*U*tt24+V*tt51;
    dFvdW(3,2) =2.0/3.0*U*tt29+tt15*tt45+tt21*tt84+V*tt56;
    dFvdW(3,3) =(kappa+tt62*tt21*tt1)*tt104+rho*Sum_dq+tt120+2.0*U*tt30*tt10+V*tt30*tt45;
    dFvdW(3,4) =tt62*tt33*tt63+tt130+tt134+U*tt37+tt136*tt84;
    dFvdW(3,5) =-tt62*rho*tt140*tt63-tt144-2.0*U*rho*tt140*tt10-tt136*tt140*tt45;
  } 
 
   dFvdW = dFvdW * dWf_dWc;
  
  
}
void dGvdW_Turbulent_Diamond(DenseMatrix &dGvdW,  DenseMatrix &dWf_dWc, Chem2D_Quad_Block &SolnBlk, 
					   const string &Orient, const int &ii, const int &jj){
  
  //planar flow
  double  d_dWdx_dW=0;
  double  d_dWdy_dW=0;
       
  double t1,t2,t3,t4,t5,t7,t8,t10;
  double t13,t14,t18,t19,t21,t23,t24,t27,t29;
  double t31,t35,t36,t40,t41,t45,t51,t56,t57;
  double t85, t86, t87,t89, t110;
  double t140, t149, t150, t153, t154, t155; 
  double t160, t164, t194, t195, t198, t204; 
  
  double kappa, Cp, mu, mu_turb, Dm_turb, Pr_turb;
  double sigma, sigma_star;
  double rho, U, V, k, omega;
  double dmudT; 
  double  *Dm, *h, *dcdy, *dhdT;
  double dUdx,dUdy, dVdx,dVdy;
  double dkdx,dkdy, domegadx,domegady;
  double drhody;
  int ns = SolnBlk.W[0][0].ns;
  double Rmix = SolnBlk.W[ii][jj].Rtot();    
  Chem2D_pState  QuadraturePoint_W;
  Dm = new double [ns];
  h  = new double [ns];
  dcdy = new double [ns];
  dhdT = new double [ns];
  //needs to deallocate the momery
 
  if (Orient == "NORTH"){ 
     
    QuadraturePoint_W = HALF*(SolnBlk.WnNW(ii, jj) +SolnBlk.WnNE(ii, jj) );
    kappa = HALF*(SolnBlk.WnNW_kappa(ii, jj) +SolnBlk.WnNE_kappa(ii, jj)); 
    mu = HALF*(SolnBlk.WnNW_mu(ii, jj) +SolnBlk.WnNE_mu(ii, jj)); 
    dmudT = HALF*(SolnBlk.WnNW_dmudT(ii, jj) +SolnBlk.WnNE_dmudT(ii, jj)); 
    for(int Num = 0; Num<ns; Num++){
      h[Num] = HALF*(SolnBlk.WnNW_hi(ii, jj, Num) +SolnBlk.WnNE_hi(ii, jj, Num));
      Dm[Num] = HALF*(SolnBlk.WnNW_Dm(ii, jj, Num) +SolnBlk.WnNE_Dm(ii, jj, Num));
      dcdy[Num] = SolnBlk.dWdy_faceN[ii][jj].spec[Num].c;
      dhdT[Num] = HALF*(SolnBlk.WnNW_Cpi(ii, jj, Num) +SolnBlk.WnNE_Cpi(ii, jj, Num));
    }
    drhody = SolnBlk.dWdx_faceN[ii][jj].rho;
    dUdx = SolnBlk.dWdx_faceN[ii][jj].v.x;
    dUdy = SolnBlk.dWdy_faceN[ii][jj].v.x;
    dVdx = SolnBlk.dWdx_faceN[ii][jj].v.y;
    dVdy = SolnBlk.dWdy_faceN[ii][jj].v.y;
 
    dkdx = SolnBlk.dWdx_faceN[ii][jj].k;
    dkdy = SolnBlk.dWdy_faceN[ii][jj].k;
    domegadx = SolnBlk.dWdx_faceN[ii][jj].omega;
    domegady = SolnBlk.dWdy_faceN[ii][jj].omega;
    
    d_dWdx_dW =  SolnBlk.d_dWdx_dW[ii][jj][1];
    d_dWdy_dW =  SolnBlk.d_dWdy_dW[ii][jj][1];

    dWfdWc(SolnBlk, dWf_dWc,Orient,  ii, jj);  
   
    }
    
  if (Orient == "EAST"){ 
 
    QuadraturePoint_W = HALF*(SolnBlk.WnNE(ii, jj) +SolnBlk.WnSE(ii, jj) );
    kappa = HALF*(SolnBlk.WnNE_kappa(ii, jj) +SolnBlk.WnSE_kappa(ii, jj)); 
    mu = HALF*(SolnBlk.WnNE_mu(ii, jj) +SolnBlk.WnSE_mu(ii, jj)); 
    dmudT = HALF*(SolnBlk.WnNE_dmudT(ii, jj) +SolnBlk.WnSE_dmudT(ii, jj)); 
    
    for(int Num = 0; Num<ns; Num++){
      h[Num] = HALF*(SolnBlk.WnNE_hi(ii, jj, Num) +SolnBlk.WnSE_hi(ii, jj, Num));
      Dm[Num] = HALF*(SolnBlk.WnNE_Dm(ii, jj, Num) +SolnBlk.WnSE_Dm(ii, jj, Num));
      dcdy[Num] = SolnBlk.dWdy_faceE[ii][jj].spec[Num].c;
      dhdT[Num] = HALF*(SolnBlk.WnNE_Cpi(ii, jj, Num) +SolnBlk.WnSE_Cpi(ii, jj, Num));
    }
    drhody = SolnBlk.dWdx_faceE[ii][jj].rho;
    dUdx = SolnBlk.dWdx_faceE[ii][jj].v.x;
    dUdy = SolnBlk.dWdy_faceE[ii][jj].v.x;
    dVdx = SolnBlk.dWdx_faceE[ii][jj].v.y;
    dVdy = SolnBlk.dWdy_faceE[ii][jj].v.y;
   
    dkdx = SolnBlk.dWdx_faceE[ii][jj].k;
    dkdy = SolnBlk.dWdy_faceE[ii][jj].k;
    domegadx = SolnBlk.dWdx_faceE[ii][jj].omega;
    domegady = SolnBlk.dWdy_faceE[ii][jj].omega;
     
    d_dWdx_dW =  SolnBlk.d_dWdx_dW[ii][jj][2];
    d_dWdy_dW =  SolnBlk.d_dWdy_dW[ii][jj][2];
    dWfdWc(SolnBlk, dWf_dWc,Orient,  ii, jj);
   
  }
  if (Orient == "SOUTH"){ 
    
    QuadraturePoint_W = HALF*(SolnBlk.WnSE(ii, jj) +SolnBlk.WnSW(ii, jj) );
    kappa = HALF*(SolnBlk.WnSW_kappa(ii, jj) +SolnBlk.WnSE_kappa(ii, jj)); 
    mu = HALF*(SolnBlk.WnSW_mu(ii, jj) +SolnBlk.WnSE_mu(ii, jj)); 
    dmudT = HALF*(SolnBlk.WnSW_dmudT(ii, jj) +SolnBlk.WnSE_dmudT(ii, jj)); 
    
    for(int Num = 0; Num<ns; Num++){
      h[Num] = HALF*(SolnBlk.WnSW_hi(ii, jj, Num) +SolnBlk.WnSE_hi(ii, jj, Num));
      Dm[Num] = HALF*(SolnBlk.WnSW_Dm(ii, jj, Num) +SolnBlk.WnSE_Dm(ii, jj, Num));
      dcdy[Num] = SolnBlk.dWdy_faceS[ii][jj].spec[Num].c;
      dhdT[Num] = HALF*(SolnBlk.WnSW_Cpi(ii, jj, Num) +SolnBlk.WnSE_Cpi(ii, jj, Num));
      
    }
    drhody = SolnBlk.dWdx_faceS[ii][jj].rho;
    dUdx = SolnBlk.dWdx_faceS[ii][jj].v.x;
    dUdy = SolnBlk.dWdy_faceS[ii][jj].v.x;
    dVdx = SolnBlk.dWdx_faceS[ii][jj].v.y;
    dVdy = SolnBlk.dWdy_faceS[ii][jj].v.y;
   
    dkdx = SolnBlk.dWdx_faceS[ii][jj].k;
    dkdy = SolnBlk.dWdy_faceS[ii][jj].k;
    domegadx = SolnBlk.dWdx_faceS[ii][jj].omega;
    domegady = SolnBlk.dWdy_faceS[ii][jj].omega;
    
    d_dWdx_dW =  SolnBlk.d_dWdx_dW[ii][jj][4];
    d_dWdy_dW =  SolnBlk.d_dWdy_dW[ii][jj][4];
    dWfdWc(SolnBlk, dWf_dWc,Orient,  ii, jj);
   
  }

  if (Orient == "WEST"){ 
    
    QuadraturePoint_W = HALF*(SolnBlk.WnNW(ii, jj) +SolnBlk.WnSW(ii, jj));
    kappa = HALF*(SolnBlk.WnNW_kappa(ii, jj) +SolnBlk.WnSW_kappa(ii, jj)); 
    mu = HALF*(SolnBlk.WnNW_mu(ii, jj) +SolnBlk.WnSW_mu(ii, jj)); 
    dmudT = HALF*(SolnBlk.WnNW_dmudT(ii, jj) +SolnBlk.WnSW_dmudT(ii, jj)); 
      
    for(int Num = 0; Num<ns; Num++){
      h[Num] = HALF*(SolnBlk.WnNW_hi(ii, jj, Num) +SolnBlk.WnSW_hi(ii, jj, Num));
      Dm[Num] = HALF*(SolnBlk.WnNW_Dm(ii, jj, Num) +SolnBlk.WnSW_Dm(ii, jj, Num));
      dcdy[Num] = SolnBlk.dWdy_faceW[ii][jj].spec[Num].c;
      dhdT[Num] = HALF*(SolnBlk.WnNW_Cpi(ii, jj, Num) +SolnBlk.WnSW_Cpi(ii, jj, Num));

    }
    drhody = SolnBlk.dWdx_faceW[ii][jj].rho;
    dUdx = SolnBlk.dWdx_faceW[ii][jj].v.x;
    dUdy = SolnBlk.dWdy_faceW[ii][jj].v.x;
    dVdx = SolnBlk.dWdx_faceW[ii][jj].v.y;
    dVdy = SolnBlk.dWdy_faceW[ii][jj].v.y;

    dkdx = SolnBlk.dWdx_faceW[ii][jj].k;
    dkdy = SolnBlk.dWdy_faceW[ii][jj].k;
    domegadx = SolnBlk.dWdx_faceW[ii][jj].omega;
    domegady = SolnBlk.dWdy_faceW[ii][jj].omega;
    
    d_dWdx_dW =  SolnBlk.d_dWdx_dW[ii][jj][3];
    d_dWdy_dW =  SolnBlk.d_dWdy_dW[ii][jj][3];
    dWfdWc(SolnBlk, dWf_dWc,Orient,  ii, jj);
    
  }

  
  mu_turb =QuadraturePoint_W.rho*QuadraturePoint_W.k/max(QuadraturePoint_W.omega, TOLER);
  Pr_turb = SolnBlk.W[ii][jj].Pr_turb(); 
  Dm_turb = SolnBlk.W[ii][jj].Dm_turb(); 
  sigma = SolnBlk.W[0][0].sigma;
  sigma_star = SolnBlk.W[0][0].sigma_star;  
  
  rho = QuadraturePoint_W.rho;
  U =  QuadraturePoint_W.v.x;
  V =  QuadraturePoint_W.v.y;
  k = QuadraturePoint_W.k;
  omega = QuadraturePoint_W.omega;
  
  t1 = 1/max(omega,TOLER);
  t2 = k*t1;
  t3 = dUdy; 
  t4 = dVdx;
  t5 = t3+t4;
  t7 = mu;
  t8 = d_dWdy_dW;
  t10 = rho*k;
  t13 = t7*t8+t10*t1*t8;
  t14 =  d_dWdx_dW;
  t18 = t7*t14+t10*t1*t14;
  t19 = dmudT;
  t21 = rho*t1;
  t23 = omega*omega;
  t24 = 1/t23;
  t27 = dVdy;
  t29 = dUdx;
  t31 = 2.0/3.0*t27-t29/3.0;
  t35 = 2.0*t2*t31-2.0/3.0*k;
  t36 =  d_dWdx_dW;
  t40 = -t7*t36-t10*t1*t36;
  t41 = d_dWdy_dW;
  t45 = t7*t41+t10*t1*t41;
  t51 = 2.0*t21*t31-2.0/3.0*rho;
  t56 = Cp/Pr_turb;
  //t57 =  SolnBlk.dTdy(ii,jj);
  t57= (ONE/(rho*Rmix)) * (QuadraturePoint_W.p - 
                 (QuadraturePoint_W.p/rho)*drhody);
  t85 = dkdy;
  t86 = t1*t85;
  t87 = sigma_star*k*t86;
  t89 = t1*t5;
  t110 =  d_dWdy_dW;

  t140 = t19*t85;
  t149 = sigma_star*rho;
  t150 = t149*t86;
  t153 =  d_dWdy_dW;
  t154 = (t7+t149*t2)*t153;
  t155 = U*rho;
  t160 = k*t24;
  t164 = t149*t160*t85;
      
  t194 =  domegady;
  t195 = t1*t194;
  t198 = sigma*rho;
  t204 = d_dWdy_dW;
  
  dGvdW(1,0) = t2*t5;
  dGvdW(1,1) = t13;
  dGvdW(1,2) = t18;
  dGvdW(1,3) = t19*t5;
  dGvdW(1,4) = t21*t5;
  dGvdW(1,5) = -t10*t24*t5;
    
  dGvdW(2,0) = t35;
  dGvdW(2,1) = 2.0/3.0*t40;
  dGvdW(2,2) = 4.0/3.0*t45;
  dGvdW(2,3) = 2.0*t19*t31;
  dGvdW(2,4) = t51;
  dGvdW(2,5) = -2.0*t10*t24*t31;
 
  double Sum_q = 0.0;
  double Sum_dq = 0.0;
 
  for(int Num = 0; Num<ns; Num++)
    {
      //for each species	// h*(Dm+Dmt)*gradc  
      Sum_q +=  h[Num]*(Dm[Num]+Dm_turb)*dcdy[Num];
      //dhdT *(Dm)*gradc
      Sum_dq += dhdT[Num]*(Dm[Num]+Dm_turb)*dcdy[Num];
      
    }   
   
  dGvdW(3,0) = t56*t2*t57+Sum_q+t87+U*k*t89+V*t35;
  dGvdW(3,1) = t7*t5+t10*t89+U*t13+2.0/3.0*V*t40;
  dGvdW(3,2) = U*t18+2.0*t7*t31+2.0*t10*t1*t31-2.0/3.0*t10+4.0/3.0*V*t45;
  dGvdW(3,3) = (kappa+t56*t10*t1)*t110+ rho*Sum_dq+t140+U*t19*t5+2.0*V*t19*t31;
  dGvdW(3,4) = t56*t21*t57+t150+t154+t155*t89+V*t51;
  dGvdW(3,5) = -t56*rho*t160*t57-t164-t155*t160*t5-2.0*V*rho*t160*t31;
  

  //multispecies
  int NUM_VAR =NUM_CHEM2D_VAR_SANS_SPECIES; 
  for(int Num = 0; Num<(ns-1); Num++)
    {//h*rho*(Dm+Dmt)*d_dWdy_dW
      dGvdW(3, NUM_VAR+Num) = h[Num]*(Dm[Num]+Dm_turb)*rho*d_dWdy_dW;
      //(Dm+Dmt)*gradc
      dGvdW(NUM_VAR+Num, 0) = (Dm[Num]+Dm_turb)*dcdy[Num];
      //rho*(Dm+Dmt)*dcydc
      dGvdW(NUM_VAR+Num,NUM_VAR+Num) = rho*(Dm[Num]+Dm_turb)*d_dWdy_dW;
    }

  dGvdW(4,0) = t87;
  dGvdW(4,3) = t140;
  dGvdW(4,4) = t150+t154;
  dGvdW(4,5) = -t164;
  
  dGvdW(5,0) = sigma*k*t195;
  dGvdW(5,3) = t19*t194;
  dGvdW(5,4) = t198*t195;
  dGvdW(5,5) = -t198*t160*t194+(t7+t198*t2)*t204;
  
  
  //The follwoing entries are different for axisymmetric flow and planar flow
  // For matrix row 2 and row 3
  if (SolnBlk.Axisymmetric == 1) {
    
    double term1,term2,term3,term4, term5,term7,term8,term10;
    double term13, term14,term18,term19;
    double term21,term22,term23,term24,term27,term29,term31,term34, term38, term39, term43, term44;
    double term47, term51,term57, term62, term63, term91, term92,term93, term95, term116; 
    double term146, term155, term156, term159, term160, term161, term166,term170 , term200, term201, term204, term210; 
    
    double r = SolnBlk.Grid.Cell[ii][jj].Xc.y;
    
    term1 = 1/max(omega,TOLER);
    term2 = k*term1;
    term3 = dUdy; 
    term4 = dVdx;
    term5 = term3+term4;
    term7 = mu;
    term8 = d_dWdy_dW;
    term10 = rho*k;
    term13 = term7*term8+term10*term1*term8;
    term14 = d_dWdx_dW;
    term18 = term7*term14+term10*term1*term14;
    term19 = dmudT;
    term21 = rho*term1;
    term23 = omega*omega;
    term24 = 1/term23;
    term27 = dVdy;
    term29 = dUdx; 
    term31 = 1/r;
    term34 = 2.0/3.0*term27-term29/3.0-V*term31/3.0;
    term38 = 2.0*term2*term34-2.0/3.0*k;
    term39 = d_dWdx_dW;
    term43 = -term7*term39-term10*term1*term39;
    term44 = d_dWdy_dW;
    term47 = 2.0/3.0*term44-term31/3.0;
    term51 = term7*term47+term10*term1*term47;
    term57 = 2.0*term21*term34-2.0/3.0*rho;
    term62 = Cp/Pr_turb;
    //  term63 = SolnBlk.dTdy(ii,jj);
    term63=(ONE/(rho*Rmix)) * (QuadraturePoint_W.p - 
                 (QuadraturePoint_W.p/rho)*drhody); 
    term91 = dkdy;
    term92 = term1*term91;
    term93 = sigma*k*term92;
    term95 = term1*term5;
    term116 = d_dWdy_dW;
    
    term146 = term19*term91;
    term155 = sigma*rho;
    term156 = term155*term92;
    term159 = d_dWdy_dW;
    term160 = (term7+term155*term2)*term159;
    term161 = U*rho;
    term166 = k*term24;
    term170 = term155*term166*term91;
    
    term200 = domegady;
    term201 = term1*term200;
    term204 = sigma_star*rho;
    term210 = d_dWdy_dW;
    
    dGvdW(2,0) = term38;
    dGvdW(2,1) = 2.0/3.0*term43;
    dGvdW(2,2) = 2.0*term51;
    dGvdW(2,3) = 2.0*term19*term34;
    dGvdW(2,4) = term57;
    dGvdW(2,5) = -2.0*term10*term24*term34;
    
    dGvdW(3,0) = term62*term2*term63+Sum_q+term93+U*k*term95+V*term38;
    dGvdW(3,1) = term7*term5+term10*term95+U*term13+2.0/3.0*V*term43;
    dGvdW(3,2) = U*term18+2.0*term7*term34+2.0*term10*term1*term34-2.0/3.0*term10+2.0*V*term51;
    dGvdW(3,3) = (kappa+term62*term10*term1)*term116+ rho*Sum_dq+term146+U*term19*term5+2.0*V*term19*term34;
    dGvdW(3,4) = term62*term21*term63+term156+term160+term161*term95+V*term57;
    dGvdW(3,5) = -term62*rho*term166*term63-term170-term161*term166*term5-2.0*V*rho*term166*term34;
    
  }
 if (SolnBlk.Axisymmetric == 2) {

   double tt1,tt2,tt3,tt4,tt5,tt7,tt8, tt10,tt13,tt14,tt18;
   double tt19,tt21,tt23,tt24,tt27,tt29,tt31;
   double tt34,tt38,tt39,tt40,tt44,tt45,tt49,tt55, tt60;
   double tt61,tt77,tt78,tt79,tt81, tt102;
   double tt118,tt127,tt128,tt131,tt132,tt133, tt138,tt142;
   double tt150,tt153,tt155,tt158,tt162,tt163, tt166,tt172;
   
   
   double r = SolnBlk.Grid.Cell[ii][jj].Xc.x;
 
   tt1 = 1/max(omega,TOLER);
   tt2 = k*tt1;
   tt3 = dUdy;
   tt4 = dVdx; 
   tt5 = tt3+tt4;
   tt7 = mu;
   tt8 = d_dWdy_dW;
   tt10 = rho*k;
   tt13 = tt7*tt8+tt10*tt1*tt8;
   tt14 = d_dWdx_dW;
   tt18 = tt7*tt14+tt10*tt1*tt14;
   tt19 = dmudT;
   tt21 = rho*tt1;
   tt23 = omega*omega;
   tt24 = 1/tt23;
   tt27 = dVdy;
   tt29 = dUdx;
   tt31 = 1/r;
   tt34 = 2.0/3.0*tt27-tt29/3.0-U*tt31/3.0;
   tt38 = 2.0*tt2*tt34-2.0/3.0*k;
   tt39 = d_dWdx_dW;
   tt40 = -tt39-tt31;
   tt44 = tt7*tt40/3.0+tt10*tt1*tt40/3.0;
   tt45 = d_dWdy_dW;
   tt49 = tt7*tt45+tt10*tt1*tt45;
   tt55 = 2.0*tt21*tt34-2.0/3.0*rho;
   tt60 = Cp/Pr_turb;
   //tt61 =SolnBlk.dTdy(ii,jj);
   tt61= (ONE/(rho*Rmix)) * (QuadraturePoint_W.p - 
                 (QuadraturePoint_W.p/rho)*drhody); 
   tt77 = dkdy;
   tt78 = tt1*tt77;
   tt79 = sigma*k*tt78;
   tt81 = tt1*tt5;
   tt102 = d_dWdy_dW;
   tt127 = sigma*rho;
   tt128 = tt127*tt78;
   tt131 = d_dWdy_dW;
   tt132 = (tt7+tt127*tt2)*tt131;
   tt133 = U*rho;
   tt138 = k*tt24;
   tt142 = tt127*tt138*tt77;
   tt162 = domegady;
   tt163 = tt1*tt162;
   tt166 = sigma_star*rho;
   tt172 = d_dWdy_dW;
    
    
   dGvdW(2,0) = tt38;
   dGvdW(2,1) = 2.0*tt44;
   dGvdW(2,2) = 4.0/3.0*tt49;
   dGvdW(2,3) = 2.0*tt19*tt34;
   dGvdW(2,4) = tt55;
   dGvdW(2,5) = -2.0*tt10*tt24*tt34;
    
   dGvdW(3,0) = tt60*tt2*tt61+Sum_q+tt79+U*k*tt81+V*tt38;
   dGvdW(3,1) = tt7*tt5+tt10*tt81+U*tt13+2.0*V*tt44;
   dGvdW(3,2) = U*tt18+2.0*tt7*tt34+2.0*tt10*tt1*tt34-2.0/3.0*tt10+4.0/3.0*V*tt49;
   dGvdW(3,3) = (kappa+tt60*tt1*tt10)*tt102+ rho*Sum_dq+tt118+U*tt19*t5+2.0*V*tt19*tt34;
   dGvdW(3,4) = tt60*tt21*tt61+tt128+tt132+tt133*tt81+V*tt55;
   dGvdW(3,5) = -tt60*rho*tt138*tt61-tt142-tt133*tt138*tt5-2.0*V*rho*tt138*tt34;
     
 }
 // Here the use of dGvdW is only convenient, but
 // dGvdW is the Jacobian obtained (with respect to the primitive solution variable on the cell edges )
 // to obtain the Jacobian (with respect to the primitive soltuion variable at the cell center)
 // dGvdW *dWfdWc 
 dGvdW = dGvdW * dWf_dWc;

 
}

//For diamond path, formulation of Viscous Jacobian Needs the transformation matrix
// Wf  --- primitive solution variables at cell face
// Wc  --- primitive solution variables at cell center
void dWfdWc(Chem2D_Quad_Block &SolnBlk, DenseMatrix dWf_dWc, const string &Orient,  const int &ii, const int &jj){

  int NUM_VAR_CHEM2D =  SolnBlk.NumVar()-1; 
  string Left;
  string Right;
  if (Orient == "NORTH"){
    Left = "NW";
    Right = "NE";
  }
  if (Orient == "SOUTH"){
    Left = "SW";
    Right = "SE";
  }
 if (Orient == "WEST"){
    Left = "NW";
    Right = "SW";
  }
 if (Orient == "EAST"){
    Left = "NE";
    Right = "SE";
  }

  for(int nn=0; nn<(NUM_VAR_CHEM2D); nn++){
    dWf_dWc(nn,nn) = HALF*(SolnBlk.dWn_dWc(ii,jj, Left) +SolnBlk.dWn_dWc(ii,jj, Right) );
  } 
  if(SolnBlk.Flow_Type == FLOWTYPE_LAMINAR){
    dWf_dWc(4,4) = ZERO;
    dWf_dWc(5,5) = ZERO;
  } 
  
 
}

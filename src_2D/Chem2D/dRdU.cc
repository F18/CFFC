#ifndef _CHEM2D_dRdU_INCLUDED
#include "dRdU.h"
#endif 

// All functions setup to use "Uo" values as point implicit 
// can be used in a multistage scheme and this insures 
// that all Jacobians are taken as dR/dUo.

/********************************************************
 * Routine: PointImplicitBlockJacobi                    *
 *                                                      *
 * This routine adds all the appropriate Jacobians      *    
 * based on flow type, for use as a Point Implicit      *
 * Block Jacobi precondtioner.                          *
 *                                                      * 
 ********************************************************/ 
void PointImplicitBlockJacobi(DenseMatrix &dRdU,
			      Chem2D_Quad_Block &SolnBlk,
			      Chem2D_Input_Parameters &Input_Parameters,
			      const int &ii, const int &jj){
   
  int NUM_VAR_CHEM2D =  SolnBlk.NumVar()-1; 
  static DenseMatrix dRdW(NUM_VAR_CHEM2D,NUM_VAR_CHEM2D,ZERO);  
  dRdW.zero();

  //Inviscid dRdU
  dFIdW_Inviscid(dRdW, SolnBlk, Input_Parameters, ii,jj);
  
  //Viscous dRdU
  if (SolnBlk.Flow_Type != FLOWTYPE_INVISCID) {    
    dGVdW_Viscous(dRdW, SolnBlk,Input_Parameters, ii, jj); 
  }

  // Add Source Jacobians (axisymmetric, turbulence)
  SemiImplicitBlockJacobi_dSdW(dRdW,SolnBlk,EXPLICIT,ii,jj);                          

  static DenseMatrix dWdU(NUM_VAR_CHEM2D,NUM_VAR_CHEM2D,ZERO); 
  dWdU.zero();
  // Transformation Jacobian 
  SolnBlk.Uo[ii][jj].W().dWdU(dWdU, SolnBlk.Flow_Type); 
  dRdU += dRdW*dWdU;

  // Add Source Jacobians (axisymmetric, chemistry, gravity)
  SemiImplicitBlockJacobi_dSdU(dRdU,SolnBlk,EXPLICIT,ii,jj);                 

}

/********************************************************
 * Routine: SemiImplicitBlockJacobi                     *
 *                                                      *
 * This routine adds all the appropriate source         *    
 * Jacobians based on flow type, for in semi-implicit   *
 * and implicit calculations.                           *
 *                                                      * 
 ********************************************************/ 
void SemiImplicitBlockJacobi(DenseMatrix &dSdU,
			     Chem2D_Quad_Block &SolnBlk,
			     const int &solver_type,
			     const int &ii, const int &jj){ 
  
  if( (SolnBlk.Axisymmetric && SolnBlk.Flow_Type != FLOWTYPE_INVISCID) ||
      SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) { 

    int NUM_VAR_CHEM2D =  SolnBlk.NumVar()-1; 
    static DenseMatrix dRdW(NUM_VAR_CHEM2D,NUM_VAR_CHEM2D,ZERO);  
    dRdW.zero();

    // Add Source Jacobians (viscous axisymmetric, turbulence)
    SemiImplicitBlockJacobi_dSdW(dRdW,SolnBlk,EXPLICIT,ii,jj);                          
    
    static DenseMatrix dWdU(NUM_VAR_CHEM2D,NUM_VAR_CHEM2D,ZERO);
    dWdU.zero();
    
    // Transformation Jacobian 
    SolnBlk.W[ii][jj].dWdU(dWdU, SolnBlk.Flow_Type);   
    dSdU += dRdW*dWdU;
  }

  // Add Source Jacobians (inviscid axisymmetric, chemistry, gravity)
  SemiImplicitBlockJacobi_dSdU(dSdU,SolnBlk,EXPLICIT,ii,jj);                 

}

void SemiImplicitBlockJacobi_dSdW(DenseMatrix &dSdW,
				  Chem2D_Quad_Block &SolnBlk,
				  const int &solver_type,
				  const int &ii, const int &jj){ 
  
  //Cacluate 2nd derivatives  
  double d_dWdx_dW_C,d_dWdy_dW_C;
  d_dWd_dW_Center(d_dWdx_dW_C,d_dWdy_dW_C,SolnBlk,ii, jj);  
  
  // Viscous Axisymmetric source term jacobian    
  if(SolnBlk.Axisymmetric && SolnBlk.Flow_Type != FLOWTYPE_INVISCID){         
    SolnBlk.W[ii][jj].dSa_vdW(dSdW,
			      SolnBlk.dWdx[ii][jj],
			      SolnBlk.dWdy[ii][jj],
			      SolnBlk.Grid.Cell[ii][jj].Xc,
			      SolnBlk.Flow_Type, SolnBlk.Axisymmetric,
			      d_dWdx_dW_C,d_dWdy_dW_C);
  }
  
  // Add Jacobian for turbulence
  if((SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ) ||
     (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_EPSILON )){    
    dS_tdW(dSdW,SolnBlk, d_dWdx_dW_C, d_dWdy_dW_C, ii,jj);
  }
  
}

void SemiImplicitBlockJacobi_dSdU(DenseMatrix &dSdU,
				  Chem2D_Quad_Block &SolnBlk,
				  const int &solver_type,
				  const int &ii, const int &jj){
  
   //Add Jacobian for inviscid axisymmetric source terms
  if (SolnBlk.Axisymmetric) {
    // Inviscid source Jacobian
    SolnBlk.W[ii][jj].dSa_idU(dSdU,
			      SolnBlk.Grid.Cell[ii][jj].Xc, 
			      SolnBlk.Flow_Type,
			      SolnBlk.Axisymmetric);
  }
  
  //Add Jacobian for finite-rate chemistry source terms  
  if (SolnBlk.W[ii][jj].React.reactset_flag != NO_REACTIONS){    
    SolnBlk.W[ii][jj].dSwdU(dSdU, SolnBlk.Flow_Type,solver_type);    
  }  

  //Add Jacobian for gravitational source terms
  if (SolnBlk.Gravity){
    SolnBlk.W[ii][jj].dSgdU(dSdU);
  } 
  
}

/********************************************************
 * Routine: PointImplicitBlkJ Inviscid Flux Jacobian    *
 *                                                      *
 * This routine returns the inviscid components of      *    
 * Point Implicit Block Jacobian matrix for the         *
 * specified local solution block.                      *
 *                                                      *
 ********************************************************/
// Based on HLLE || ROE Flux Function
void dFIdW_Inviscid(DenseMatrix &dRdW, Chem2D_Quad_Block &SolnBlk, Chem2D_Input_Parameters 
		    &Input_Parameters,const int &ii, const int &jj){
  

  if (Input_Parameters.i_Flux_Function == FLUX_FUNCTION_HLLE ){    
    dFIdW_Inviscid_HLLE(dRdW,SolnBlk,Input_Parameters,ii, jj, NORTH);
    dFIdW_Inviscid_HLLE(dRdW,SolnBlk,Input_Parameters,ii, jj, SOUTH);
    dFIdW_Inviscid_HLLE(dRdW,SolnBlk,Input_Parameters,ii, jj, EAST);
    dFIdW_Inviscid_HLLE(dRdW,SolnBlk,Input_Parameters,ii, jj, WEST);
    dRdW = dRdW/SolnBlk.Grid.Cell[ii][jj].A;
    
  } else if (Input_Parameters.i_Flux_Function == FLUX_FUNCTION_ROE ){
    
    dFIdW_Inviscid_ROE(dRdW,SolnBlk,Input_Parameters,ii, jj, NORTH);
    dFIdW_Inviscid_ROE(dRdW,SolnBlk,Input_Parameters,ii, jj, SOUTH);
    dFIdW_Inviscid_ROE(dRdW,SolnBlk,Input_Parameters,ii, jj, EAST);
    dFIdW_Inviscid_ROE(dRdW,SolnBlk,Input_Parameters,ii, jj, WEST);
    dRdW = dRdW/SolnBlk.Grid.Cell[ii][jj].A;
    
  } else {
    cerr<<"\n NOT A VALID FLUX FUNCTION FOR USE WITH Point Implicit \n";
  }
  
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
void Rotation_Matrix2(DenseMatrix &mat, Vector2D nface,  int A_matrix) 
{

  double cos_angle = nface.x; 
  double sin_angle = nface.y;    
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
 * Routine: Inviscid Flux Jacobian using HLLE           *
 *                                                      *
 * This routine returns the inviscid components of      *    
 * Point Implicit Block Jacobian matrix for the         *
 * specified local solution block.                      *
 *                                                      *
 ********************************************************/
void dFIdW_Inviscid_HLLE(DenseMatrix &dRdW, Chem2D_Quad_Block &SolnBlk,
			 Chem2D_Input_Parameters &Input_Parameters, 
			 const int &ii, const int &jj, const int Orient){

  int NUM_VAR_CHEM2D =  SolnBlk.NumVar(); 
  int overlap = Input_Parameters.NKS_IP.GMRES_Overlap;

   static DenseMatrix dFidW(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D,ZERO);
   dFidW.zero();
   Vector2D nface, lambdas;   
    
   if (ii < SolnBlk.ICl -overlap || ii > SolnBlk.ICu +overlap ||
       jj < SolnBlk.JCl -overlap || jj > SolnBlk.JCu +overlap) {
     //Ghost cell so shouldn't be here
     exit(1);
   } else if (Orient == NORTH) {
     nface = SolnBlk.Grid.nfaceN(ii, jj-1);
     lambdas = HLLE_wavespeeds(SolnBlk.Uo[ii][jj-1].W(), 
			       SolnBlk.Uo[ii][jj].W(), nface);
   } else if (Orient == SOUTH) {
     nface = SolnBlk.Grid.nfaceS(ii, jj+1);
     lambdas = HLLE_wavespeeds(SolnBlk.Uo[ii][jj+1].W(), 
			       SolnBlk.Uo[ii][jj].W(), nface);
   } else if (Orient == EAST) {
     nface = SolnBlk.Grid.nfaceE(ii-1, jj);     
     lambdas = HLLE_wavespeeds(SolnBlk.Uo[ii-1][jj].W(), 
			       SolnBlk.Uo[ii][jj].W(), nface);
   } else if (Orient == WEST) {
     nface = SolnBlk.Grid.nfaceW(ii+1, jj);
     lambdas = HLLE_wavespeeds(SolnBlk.Uo[ii+1][jj].W(), 
			       SolnBlk.Uo[ii][jj].W(), nface);
   }

   static DenseMatrix A(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D,ZERO); 
   static DenseMatrix AI(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D,ZERO); 
   Rotation_Matrix2(A,nface, 1);
   Rotation_Matrix2(AI,nface, 0);
   static DenseMatrix II(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D,ZERO); 
   II.identity();
   
   //Weightings
   double gamma = (lambdas.x*lambdas.y)/(lambdas.y-lambdas.x);
   double beta  = - lambdas.x/(lambdas.y-lambdas.x);

   //Jacobian
   dFIdW(dFidW, Rotate(SolnBlk.Uo[ii][jj].W(), nface) , SolnBlk.Flow_Type);     

   if (Orient == NORTH) {                 
     dRdW += SolnBlk.Grid.lfaceN(ii, jj-1)*(AI*(beta*dFidW + gamma*II)*A); 
   } else if (Orient == SOUTH) {     
     dRdW += SolnBlk.Grid.lfaceS(ii, jj+1)*(AI*(beta*dFidW + gamma*II)*A); 
   } else if (Orient == EAST) {                
     dRdW += SolnBlk.Grid.lfaceE(ii-1, jj)*(AI*(beta*dFidW + gamma*II)*A);      
   } else if (Orient == WEST) {     
     dRdW += SolnBlk.Grid.lfaceW(ii+1, jj)*(AI*(beta*dFidW + gamma*II)*A); 
   } else {
     cerr<<" NOT A VALID ORIENTATION "; exit(1);
   }

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
void dFIdW_Inviscid_ROE(DenseMatrix& dRdW, Chem2D_Quad_Block &SolnBlk,  
			Chem2D_Input_Parameters &Input_Parameters,
			const int &ii, const int &jj, const int Orient){
   
  int overlap = Input_Parameters.NKS_IP.GMRES_Overlap;
  int Ri, Rj;

  if (ii < SolnBlk.ICl -overlap || ii > SolnBlk.ICu + overlap ||
      jj < SolnBlk.JCl -overlap || jj > SolnBlk.JCu + overlap) {

     // GHOST CELL so do nothing
     cout<<"\n Hey I am not suppose to be here! \n"; exit(1);

  } else {     
    int NUM_VAR_CHEM2D = dRdW.get_n();  //  SolnBlk.NumVar()-1;    
    static DenseMatrix dFidW(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D,ZERO);
    dFidW.zero();

    Vector2D nface,DX; double lface;   
    Chem2D_pState Wa, wavespeeds, Left, Right, Wl, Wr;   
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
     
     static DenseMatrix A(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D,ZERO); 
     static DenseMatrix AI(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D,ZERO); 
     Rotation_Matrix2(A,nface, 1);
     Rotation_Matrix2(AI,nface, 0);
     
     //Left and Right States                                                       ///Ri ,Rj fixed not ii, jj 
     Inviscid_Flux_Used_Reconstructed_LeftandRight_States(Wl, Wr, DX, SolnBlk, Orient, Ri, Rj );   
     Left  = Rotate(Wl, nface);
     Right = Rotate(Wr, nface);
     
     //Determin Roe Averaged State
     Wa = RoeAverage(Left,Right);       

     // Jacobian dF/dW         
     dFIdW(dFidW, Rotate(SolnBlk.Uo[ii][jj].W(), nface) , SolnBlk.Flow_Type);       
     //dFidW = HALF*dFidW;
     dFidW *= HALF; 
     
     /***************************** Regular Roe (no preconditioning) *************************************/
     if(!Input_Parameters.Preconditioning){
       // Determine Wave Speeds
       wavespeeds = HartenFixAbs( Wa.lambda_x(),
				  Left.lambda_x(),
				  Right.lambda_x());         
       
       //Loop through each wavespeed and each element of Jacobian(i,j)        
       for (int i=1; i <= NUM_VAR_CHEM2D; i++) {		   
	 for(int irow =0; irow< NUM_VAR_CHEM2D; irow++){
	   for(int jcol =0; jcol< NUM_VAR_CHEM2D; jcol++){

	     dFidW(irow, jcol) -= HALF*wavespeeds[i]*Wa.lp_x(i)[jcol+1]*Wa.rc_x(i)[irow+1];  
		   
             //     2nd Order terms	   
	     //        if(irow ==jcol){  // now the rotated Wc is used in Wr = Wc + phi*grad(Wc), chain rule...
	     //    	 dFidW(irow,jcol) -= SolnBlk.phi[ii][jj][irow+1]*(SolnBlk.d_dWdx_dW[ii][jj][0]*DX.x 
	     //                          + SolnBlk.d_dWdy_dW[ii][jj][0]*DX.y);
	     //        }
	   }
	 }
       } 
       
       /****************************** LOW MACH NUMBER PRECONDITIONING ************************************/
     } else if(Input_Parameters.Preconditioning){
       
       //THIS MAY NOT BE CONSISTENT !!!!!!!!!!!
       double deltax = min(TWO*(SolnBlk.Grid.Cell[ii][jj].A/(SolnBlk.Grid.lfaceE(ii, jj)+SolnBlk.Grid.lfaceW(ii, jj))),
			   TWO*(SolnBlk.Grid.Cell[ii][jj].A/(SolnBlk.Grid.lfaceN(ii, jj)+SolnBlk.Grid.lfaceS(ii, jj))));
       
       double MR2a = Wa.Mr2(SolnBlk.Flow_Type,deltax);  
       // Determine Preconditioned Wave Speeds                                                                   
       wavespeeds = HartenFixAbs( Wa.lambda_preconditioned_x(MR2a),
				  Left.lambda_preconditioned_x(Left.Mr2(SolnBlk.Flow_Type,deltax)),
				  Right.lambda_preconditioned_x(Right.Mr2(SolnBlk.Flow_Type,deltax)));
       
       
       //Calculate the preconditioned upwind dissipation flux.
       static DenseMatrix Flux_dissipation(NUM_VAR_CHEM2D,NUM_VAR_CHEM2D,ZERO); 
       Flux_dissipation.zero();
       for (int i=1; i <=  NUM_VAR_CHEM2D; i++) {		   
	 for(int irow =0; irow < NUM_VAR_CHEM2D; irow++){
	   for(int jcol =0; jcol < NUM_VAR_CHEM2D; jcol++){	   
	     Flux_dissipation(irow, jcol) -= HALF*wavespeeds[i]*Wa.lp_x_precon(i,MR2a)[jcol+1]*Wa.rc_x_precon(i,MR2a)[irow+1];   
	   }
	 }
       }
       
       // Evaluate the low-Mach-number local preconditioner for the Roe-averaged state.
       static DenseMatrix P( NUM_VAR_CHEM2D, NUM_VAR_CHEM2D,ZERO);    
       Wa.Low_Mach_Number_Preconditioner(P,SolnBlk.Flow_Type,deltax);
       
       // Add preconditioned dissipation to Inviscid Jacobian
       dFidW += P*Flux_dissipation;
       
     }     
     
     //Rotate back 
     dRdW += lface*AI*dFidW*A;
        
  } 
}

/********************************************************
 * Routine: Inviscid Roe Flux Jacobian                  *
 *                                                      *
 *     Finite Differences                               *
 *                                                      *
 ********************************************************/
void dFIdW_Inviscid_ROE_FD(DenseMatrix& dRdW, Chem2D_Quad_Block &SolnBlk,  
			   Chem2D_Input_Parameters &Input_Parameters,
			   const int &ii, const int &jj, const int Orient){
   
  int overlap = Input_Parameters.NKS_IP.GMRES_Overlap;
  int Ri, Rj;

  if (ii < SolnBlk.ICl -overlap || ii > SolnBlk.ICu + overlap ||
      jj < SolnBlk.JCl -overlap || jj > SolnBlk.JCu + overlap) {

     // GHOST CELL so do nothing
     cout<<"\n Hey I am not suppose to be here! \n"; exit(1);

  } else {     
     int NUM_VAR_CHEM2D = dRdW.get_n();   
     static DenseMatrix dFidW(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D,ZERO);
     dFidW.zero();

     Vector2D nface, DX;   double lface;
     Chem2D_pState Wa, wavespeeds, Left, Right, Wl, Wr;   
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
     
     static DenseMatrix A(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D,ZERO); 
     static DenseMatrix AI(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D,ZERO); 
     Rotation_Matrix2(A,nface, 1);
     Rotation_Matrix2(AI,nface, 0);
          
     //Left and Right States 
     Inviscid_Flux_Used_Reconstructed_LeftandRight_States(Wl, Wr, DX, SolnBlk, Orient, Ri, Rj );   
     Left  = Rotate(Wl, nface);
     Right = Rotate(Wr, nface);
     
     /**********************************************/
     //Jacobain using Finite Differences
     Chem2D_cState FluxA, FluxB; 
     Chem2D_pState WA, WB;
  
     double perturb = 5e-6;
     double a;
     
     //For Preconditioning
     double delta_n = min(TWO*(SolnBlk.Grid.Cell[ii][jj].A/
			      (SolnBlk.Grid.lfaceE(ii, jj)+SolnBlk.Grid.lfaceW(ii, jj))),
			 TWO*(SolnBlk.Grid.Cell[ii][jj].A/
			      (SolnBlk.Grid.lfaceN(ii, jj)+SolnBlk.Grid.lfaceS(ii, jj))));
 
     for(int jcol=0; jcol<(NUM_VAR_CHEM2D); jcol++){
       WA = Right;
       WB = Right;
       
       if( jcol <NUM_CHEM2D_VAR_SANS_SPECIES) {
	 WA[jcol+1] += perturb*max(ONE,Right[jcol+1]); 	 
	 WB[jcol+1] -= perturb*max(ONE,Right[jcol+1]); 
       } else {
	 a =  perturb*max(ONE,Right[jcol+1]); 
	 WA[jcol+1] += a;
	 WA[NUM_VAR_CHEM2D+1] -= a;      
	 WB[jcol+1] -= a;
	 WB[NUM_VAR_CHEM2D+1] += a;
       }

       FluxA = FluxRoe_x(Left,WA, Input_Parameters.Preconditioning,SolnBlk.Flow_Type,delta_n);       
       FluxB = FluxRoe_x(Left,WB, Input_Parameters.Preconditioning,SolnBlk.Flow_Type,delta_n);
       
       for(int irow=0; irow<(NUM_VAR_CHEM2D); irow++){
	 dFidW(irow,jcol) = (FluxA[irow+1] - FluxB[irow+1])/(TWO*perturb*max(ONE,Right[jcol+1]));           
       }
     } 
     /**********************************************/
        
     //Rotate back 
     dRdW += lface*AI*dFidW*A;

  } 
}

//ARE the BC checks required if everything calculated using Uo ???????

//Needed by inviscid flux Jacobian -- reconstructed higher order left and right solution states
int Inviscid_Flux_Used_Reconstructed_LeftandRight_States(Chem2D_pState &Wl, Chem2D_pState &Wr, 
                                                         Vector2D &DX, Chem2D_Quad_Block &SolnBlk, 
                                                         const int &Orient, const int &i, const int &j ){  
  switch(Orient) {
  case WEST:
    //Wl = SolnBlk.W[i][j];  
     Wl = SolnBlk.Uo[i][j].W(); 
     //	   + ((SolnBlk.phi[i][j]^SolnBlk.dWdx[i][j])*dX.x + (SolnBlk.phi[i][j]^SolnBlk.dWdy[i][j])*dX.y); 
     //DX = SolnBlk.Grid.xfaceW(i+1, j)-SolnBlk.Grid.Cell[i+1][j].Xc;      
 
     /* Reconstruct left and right solution states at the east-west faces */
     if (i == SolnBlk.ICl && 
	 (SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION ||
	  SolnBlk.Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
	  SolnBlk.Grid.BCtypeW[j] == BC_FREE_SLIP_ISOTHERMAL ||
	  SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
	  SolnBlk.Grid.BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL ||

	  SolnBlk.Grid.BCtypeW[j] == BC_1DFLAME_INFLOW ||
	  SolnBlk.Grid.BCtypeW[j] == BC_1DFLAME_OUTFLOW ||
	  SolnBlk.Grid.BCtypeW[j] == BC_2DFLAME_INFLOW ||
	  SolnBlk.Grid.BCtypeW[j] == BC_2DFLAME_OUTFLOW )) {
       
       if (SolnBlk.Grid.BCtypeW[j] == BC_REFLECTION) {
	 Wr = Reflect(Wl, SolnBlk.Grid.nfaceW(i, j));
       } else if (SolnBlk.Grid.BCtypeW[j] == BC_FREE_SLIP_ISOTHERMAL) {
	 Wr = Free_Slip(Wl,SolnBlk.WoW[j], SolnBlk.Grid.nfaceW(i, j),FIXED_TEMPERATURE_WALL);	      
       } else if (SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL) {
	 Wr = No_Slip(Wl,SolnBlk.WoW[j], SolnBlk.Grid.nfaceW(i, j),FIXED_TEMPERATURE_WALL);
       } else if (SolnBlk.Grid.BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL) {
	 Wr = Moving_Wall(Wl,SolnBlk.WoW[j], SolnBlk.Grid.nfaceW(i, j),
			  SolnBlk.Moving_wall_velocity,FIXED_TEMPERATURE_WALL);            
       } else if (SolnBlk.Grid.BCtypeW[j] == BC_1DFLAME_INFLOW){
	 Wr = BC_1DFlame_Inflow(Wl, 
			      SolnBlk.WoW[j],
			      SolnBlk.W[SolnBlk.ICu][j],
			      SolnBlk.Grid.nfaceW(i, j));
            
       } else if (SolnBlk.Grid.BCtypeW[j] == BC_1DFLAME_OUTFLOW){
	 Wr = BC_1DFlame_Outflow(Wl, 
			       SolnBlk.WoW[j], 
			       SolnBlk.W[SolnBlk.ICu][j],             
			       SolnBlk.Grid.nfaceW(i, j));
       } else { 
	 Wr = BC_Characteristic_Pressure(Wl, SolnBlk.WoW[j], 
					 SolnBlk.Grid.nfaceW(i, j));
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
	 SolnBlk.Grid.BCtypeE[j] == BC_CHARACTERISTIC ||
	 SolnBlk.Grid.BCtypeE[j] == BC_FREE_SLIP_ISOTHERMAL ||
	 SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
	 SolnBlk.Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL || 
	 SolnBlk.Grid.BCtypeE[j] == BC_1DFLAME_INFLOW ||
	 SolnBlk.Grid.BCtypeE[j] == BC_1DFLAME_OUTFLOW )) {
      
      if (SolnBlk.Grid.BCtypeE[j] == BC_REFLECTION) {
	Wr = Reflect(Wl, SolnBlk.Grid.nfaceE(i, j));   
      } else if (SolnBlk.Grid.BCtypeE[j] == BC_FREE_SLIP_ISOTHERMAL) {
	Wr = Free_Slip(Wl, SolnBlk.WoE[j],SolnBlk.Grid.nfaceE(i, j),FIXED_TEMPERATURE_WALL);
      } else if (SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL) {
	Wr = No_Slip(Wl, SolnBlk.WoE[j],SolnBlk.Grid.nfaceE(i, j),FIXED_TEMPERATURE_WALL);
      } else if (SolnBlk.Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL) {
	Wr = Moving_Wall(Wl,SolnBlk.WoE[j], SolnBlk.Grid.nfaceE(i, j),
			 SolnBlk.Moving_wall_velocity,FIXED_TEMPERATURE_WALL);
      } else if (SolnBlk.Grid.BCtypeE[j] == BC_1DFLAME_INFLOW){
	Wr = BC_1DFlame_Inflow(Wl, 
			     SolnBlk.WoE[j],
			     SolnBlk.W[SolnBlk.ICl][j],
			     SolnBlk.Grid.nfaceE(i, j));
      } else if (SolnBlk.Grid.BCtypeE[j] == BC_1DFLAME_OUTFLOW){
	Wr = BC_1DFlame_Outflow(Wl, 
			      SolnBlk.WoE[j],
			      SolnBlk.W[SolnBlk.ICl][j],
			      SolnBlk.Grid.nfaceE(i, j));
      } else {
	Wr = BC_Characteristic_Pressure(Wl, SolnBlk.WoE[j], 
					SolnBlk.Grid.nfaceE(i, j));
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
	 SolnBlk.Grid.BCtypeS[i] == BC_CHARACTERISTIC ||
	 SolnBlk.Grid.BCtypeS[i] == BC_FREE_SLIP_ISOTHERMAL ||
	 SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
	 SolnBlk.Grid.BCtypeS[i] == BC_MOVING_WALL_ISOTHERMAL ||
	 SolnBlk.Grid.BCtypeS[i] == BC_1DFLAME_OUTFLOW )) {
      
      if (SolnBlk.Grid.BCtypeS[i] == BC_REFLECTION) {
	Wr = Reflect(Wl, SolnBlk.Grid.nfaceS(i, j));	
      } else if (SolnBlk.Grid.BCtypeS[i] == BC_FREE_SLIP_ISOTHERMAL) {
	Wr = Free_Slip(Wl,SolnBlk.WoS[i], SolnBlk.Grid.nfaceS(i, j),FIXED_TEMPERATURE_WALL);
      } else if (SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL) {
	Wr = No_Slip(Wl,SolnBlk.WoS[i], SolnBlk.Grid.nfaceS(i, j),FIXED_TEMPERATURE_WALL);
      } else if (SolnBlk.Grid.BCtypeS[i] == BC_MOVING_WALL_ISOTHERMAL) {
	Wr = Moving_Wall(Wl,SolnBlk.WoS[i], SolnBlk.Grid.nfaceS(i, j),
			 SolnBlk.Moving_wall_velocity,FIXED_TEMPERATURE_WALL);
      } else if (SolnBlk.Grid.BCtypeS[i] == BC_1DFLAME_OUTFLOW){
	Wr = BC_1DFlame_Outflow(Wl, 
			      SolnBlk.WoS[i], 
			      SolnBlk.W[i][SolnBlk.JCu],
			      SolnBlk.Grid.nfaceS(i, j+1));
      } else { 
	Wr = BC_Characteristic_Pressure(Wl, 
					SolnBlk.WoS[i], 
					SolnBlk.Grid.nfaceS(i, j));
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
	 SolnBlk.Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
	 SolnBlk.Grid.BCtypeN[i] == BC_FREE_SLIP_ISOTHERMAL ||
	 SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
	 SolnBlk.Grid.BCtypeN[i] == BC_MOVING_WALL_ISOTHERMAL ||
	 SolnBlk.Grid.BCtypeN[i] == BC_1DFLAME_OUTFLOW ||
	 SolnBlk.Grid.BCtypeN[i] == BC_2DFLAME_OUTFLOW )) {                       
      
      if (SolnBlk.Grid.BCtypeN[i] == BC_REFLECTION) {
	Wr = Reflect(Wl, SolnBlk.Grid.nfaceN(i, j));	
      } else if (SolnBlk.Grid.BCtypeN[i] == BC_FREE_SLIP_ISOTHERMAL) {
	Wr = Free_Slip(Wl, SolnBlk.WoN[i], SolnBlk.Grid.nfaceN(i, j),FIXED_TEMPERATURE_WALL);
      } else if (SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL) {
	Wr = No_Slip(Wl, SolnBlk.WoN[i], SolnBlk.Grid.nfaceN(i, j),FIXED_TEMPERATURE_WALL);
      } else if (SolnBlk.Grid.BCtypeN[i] == BC_MOVING_WALL_ISOTHERMAL) { 	
	Wr = Moving_Wall(Wl,SolnBlk.WoN[i],  SolnBlk.Grid.nfaceN(i, j),
			 SolnBlk.Moving_wall_velocity,FIXED_TEMPERATURE_WALL);
      } else if (SolnBlk.Grid.BCtypeN[i] == BC_1DFLAME_OUTFLOW){
	Wr = BC_1DFlame_Outflow(Wl, 
			      SolnBlk.WoN[i], 
			      SolnBlk.W[i][SolnBlk.JCl],
			      SolnBlk.Grid.nfaceW(i, j));
      } else if (SolnBlk.Grid.BCtypeN[i] == BC_2DFLAME_OUTFLOW){
	Wr = BC_2DFlame_Outflow(Wl, 
				SolnBlk.WoN[i], 			      
				SolnBlk.Grid.nfaceW(i, j));
      } else {
	Wr = BC_Characteristic_Pressure(Wl, 
					SolnBlk.WoN[i], 
					SolnBlk.Grid.nfaceN(i, j));
      } 
    } else {      
      //Wr = SolnBlk.W[i][j+1]; 
      Wr = SolnBlk.Uo[i][j+1].W(); 
    }    
    break;
  }

  return 0;
   
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
void dFIdW_Inviscid_AUSM_plus_up(DenseMatrix& dRdW, Chem2D_Quad_Block &SolnBlk,  
				 Chem2D_Input_Parameters &Input_Parameters,
				 const int &ii, const int &jj, const int Orient){
   
  int overlap = Input_Parameters.NKS_IP.GMRES_Overlap;
  int Ri, Rj;

  if (ii < SolnBlk.ICl -overlap || ii > SolnBlk.ICu + overlap ||
      jj < SolnBlk.JCl -overlap || jj > SolnBlk.JCu + overlap) {
     // GHOST CELL so do nothing
     cout<<"\n Hey I am not suppose to be here! \n"; exit(1);

  } else {     
    int NUM_VAR_CHEM2D = dRdW.get_n();  //  SolnBlk.NumVar();    
    static DenseMatrix dFidW(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D,ZERO);
    dFidW.zero();

    Vector2D nface,DX; double lface;   
    Chem2D_pState Wa, wavespeeds, Left, Right, Wl, Wr;   
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

     static DenseMatrix A(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D,ZERO); 
     static DenseMatrix AI(NUM_VAR_CHEM2D, NUM_VAR_CHEM2D,ZERO); 
     Rotation_Matrix2(A,nface, 1);
     Rotation_Matrix2(AI,nface, 0);
               
     //********** FILL IN AUSM SPECIFIC STUFF HERE ********************//
     Left  = Rotate(Wl, nface);
     Right = Rotate(Wr, nface);
     

     //***************************************************************//

     // Jacobian dF/dW         
     dFIdW(dFidW, Rotate(SolnBlk.Uo[ii][jj].W(), nface) , SolnBlk.Flow_Type);       
     dFidW = HALF*dFidW;
          
     //Rotate back 
     dRdW += lface*AI*dFidW*A;
        
  } 
}


/********************************************************
 * Routine: dGVdW_Viscous                               *
 *                                                      *
 * This routine calculates the Viscous components of    *
 * the residual with respect to the primitive variables *   
 * ie dG/dW                                             *
 *                                                      *
 ********************************************************/
void dGVdW_Viscous(DenseMatrix &dRdW, Chem2D_Quad_Block &SolnBlk, 
		   Chem2D_Input_Parameters &Input_Parameters,
		   const int &ii, const int &jj){
   
   int NUM_VAR_CHEM2D = SolnBlk.NumVar()-1; 
   int ns = SolnBlk.W[ii][jj].ns-1;

   static DenseMatrix dFvdWf(NUM_VAR_CHEM2D, 14+(ns),ZERO);
   static DenseMatrix dWfdWx(14+(ns), NUM_VAR_CHEM2D,ZERO);  
   static DenseMatrix dGvdWf(NUM_VAR_CHEM2D, 14+(ns),ZERO);
   static DenseMatrix dWfdWy(14+(ns), NUM_VAR_CHEM2D,ZERO);
   static DenseMatrix dGVdW(NUM_VAR_CHEM2D,NUM_VAR_CHEM2D,ZERO);

   Vector2D nface;

   if (ii < SolnBlk.ICl || ii > SolnBlk.ICu ||
       jj < SolnBlk.JCl || jj > SolnBlk.JCu) {
     // GHOST CELL so do nothing to dRdW     
   } else {
     // NON-GHOST CELL.
           
     //Viscous flux Jacobians at North face of cell (ii, jj)
     nface = SolnBlk.Grid.nfaceN(ii, jj); 
     dFvdWf.zero(); dWfdWx.zero(); dGvdWf.zero(); dWfdWy.zero();   
     dFvdWf_Diamond(dFvdWf,dGvdWf, SolnBlk, NORTH, ii, jj);
     dWfdWc_Diamond(dWfdWx,dWfdWy, SolnBlk, NORTH, ii, jj, CENTER);     
     dGVdW = SolnBlk.Grid.lfaceN(ii, jj)* (nface.x*(dFvdWf*dWfdWx) + nface.y*(dGvdWf*dWfdWy));
     
     //Viscous flux Jacobians at South face of cell (ii, jj)
     nface = SolnBlk.Grid.nfaceS(ii, jj);
     dFvdWf.zero(); dGvdWf.zero(); dWfdWx.zero(); dWfdWy.zero();
     dFvdWf_Diamond(dFvdWf,dGvdWf, SolnBlk, SOUTH, ii, jj);
     dWfdWc_Diamond(dWfdWx,dWfdWy, SolnBlk, SOUTH, ii, jj,CENTER); 
     dGVdW += SolnBlk.Grid.lfaceS(ii, jj)* (nface.x*(dFvdWf*dWfdWx) + nface.y*(dGvdWf*dWfdWy));
     
     // Viscous flux Jacobians at West face of cell (ii, jj)
     nface = SolnBlk.Grid.nfaceW(ii, jj);
     dFvdWf.zero(); dGvdWf.zero(); dWfdWx.zero(); dWfdWy.zero();
     dFvdWf_Diamond(dFvdWf,dGvdWf, SolnBlk, WEST, ii, jj);
     dWfdWc_Diamond(dWfdWx,dWfdWy, SolnBlk, WEST, ii, jj, CENTER);
     dGVdW += SolnBlk.Grid.lfaceW(ii, jj)* (nface.x*(dFvdWf*dWfdWx) + nface.y*(dGvdWf*dWfdWy));
     
     // Viscous flux Jacobians at East face of cell (ii, jj)
     nface = SolnBlk.Grid.nfaceE(ii, jj);
     dFvdWf.zero(); dGvdWf.zero(); dWfdWx.zero(); dWfdWy.zero(); 
     dFvdWf_Diamond(dFvdWf,dGvdWf, SolnBlk, EAST, ii, jj);
     dWfdWc_Diamond(dWfdWx,dWfdWy, SolnBlk, EAST, ii, jj, CENTER);     
     dGVdW += SolnBlk.Grid.lfaceE(ii, jj)* (nface.x*(dFvdWf*dWfdWx) + nface.y*(dGvdWf*dWfdWy));
     
     dRdW += dGVdW/SolnBlk.Grid.Cell[ii][jj].A;
   }  
}

/********************************************************
 * Routine: dFvdWf_Diamond                              *
 *                                                      *
 * This routine calculates the Viscous components of    *
 * the residual with respect to the primitive variables *   
 * in based on Green Gauss Diamond Path Reconstruction. *
 *                                                      *
 ********************************************************/
void dFvdWf_Diamond(DenseMatrix &dFvdWf, DenseMatrix &dGvdWf, 
		    Chem2D_Quad_Block &SolnBlk,
		    const int &Orient, const int &ii, const int &jj){
   
  
   double kappa, Cp, mu, mu_t, kappa_t,Dm_t,Pr_t, Sc_t;
   double sigma, sigma_star, Rmix;
   double rho, U, V, p, k, omega;
   double  *h, *dcdx, *dcdy, *dhdT;
   double dUdx,dUdy, dVdx,dVdy;
   double dkdx, dkdy, domegadx, domegady;
   double drhodx, drhody, dpdx, dpdy, Temp;
   double radius;
   
   Chem2D_pState QuadraturePoint_W;
   QuadraturePoint_W.Vacuum();

   int ns_values = SolnBlk.W[ii][jj].ns;
   int ns_species = SolnBlk.W[ii][jj].ns-1;
   h  = new double [ns_values];
   dcdx = new double [ns_values]; 
   dcdy = new double [ns_values];
   dhdT = new double [ns_values];
  
   switch(Orient){
     /****************************** NORTH ******************************/
   case NORTH:
     QuadraturePoint_W = HALF*(SolnBlk.UnoNW(ii, jj).W() +SolnBlk.UnoNE(ii, jj).W() );
     Temp =  QuadraturePoint_W.T();
     for(int Num = 0; Num<ns_values; Num++){
       h[Num] =  QuadraturePoint_W.specdata[Num].Enthalpy(Temp)+QuadraturePoint_W.specdata[Num].Heatofform();
       dcdx[Num] = SolnBlk.dWdx_faceN[ii][jj].spec[Num].c;
       dcdy[Num] = SolnBlk.dWdy_faceN[ii][jj].spec[Num].c;
       dhdT[Num] = QuadraturePoint_W.specdata[Num].Enthalpy_prime(Temp); //HeatCapacity_p(Temp)
     }
     
     drhodx = SolnBlk.dWdx_faceN[ii][jj].rho;
     drhody = SolnBlk.dWdy_faceN[ii][jj].rho;
     dpdx = SolnBlk.dWdx_faceN[ii][jj].p;
     dpdy = SolnBlk.dWdy_faceN[ii][jj].p;
     dUdx = SolnBlk.dWdx_faceN[ii][jj].v.x;
     dUdy = SolnBlk.dWdy_faceN[ii][jj].v.x;
     dVdx = SolnBlk.dWdx_faceN[ii][jj].v.y;
     dVdy = SolnBlk.dWdy_faceN[ii][jj].v.y;
     dkdx = SolnBlk.dWdx_faceN[ii][jj].k;
     dkdy = SolnBlk.dWdy_faceN[ii][jj].k;
     domegadx = SolnBlk.dWdx_faceN[ii][jj].omega;
     domegady = SolnBlk.dWdy_faceN[ii][jj].omega;

     if (SolnBlk.Axisymmetric == AXISYMMETRIC_X) {
       radius = (SolnBlk.Grid.xfaceN(ii,jj).x < MICRO) ? MICRO : SolnBlk.Grid.xfaceN(ii,jj).x;     
     } else if (SolnBlk.Axisymmetric == AXISYMMETRIC_Y) {    
       radius = (SolnBlk.Grid.xfaceN(ii,jj).y < MICRO) ? MICRO : SolnBlk.Grid.xfaceN(ii,jj).y;            
     } 

     break;
     /****************************** EAST ******************************/
   case EAST:
     QuadraturePoint_W = HALF*(SolnBlk.UnoNE(ii, jj).W() + SolnBlk.UnoSE(ii, jj).W() );     
     Temp =  QuadraturePoint_W.T();   
     for(int Num = 0; Num<ns_values; Num++){
       h[Num] =  QuadraturePoint_W.specdata[Num].Enthalpy(Temp)+QuadraturePoint_W.specdata[Num].Heatofform();
       dcdx[Num] = SolnBlk.dWdx_faceE[ii][jj].spec[Num].c; 
       dcdy[Num] = SolnBlk.dWdy_faceE[ii][jj].spec[Num].c;
       dhdT[Num] = QuadraturePoint_W.specdata[Num].Enthalpy_prime(Temp);
     } 
     
     drhodx = SolnBlk.dWdx_faceE[ii][jj].rho; 
     drhody = SolnBlk.dWdy_faceE[ii][jj].rho;
     dpdx = SolnBlk.dWdx_faceE[ii][jj].p;
     dpdy = SolnBlk.dWdy_faceE[ii][jj].p;
     dUdx = SolnBlk.dWdx_faceE[ii][jj].v.x;
     dUdy = SolnBlk.dWdy_faceE[ii][jj].v.x;
     dVdx = SolnBlk.dWdx_faceE[ii][jj].v.y;
     dVdy = SolnBlk.dWdy_faceE[ii][jj].v.y;
     dkdx = SolnBlk.dWdx_faceE[ii][jj].k;
     dkdy = SolnBlk.dWdy_faceE[ii][jj].k;
     domegadx = SolnBlk.dWdx_faceE[ii][jj].omega;
     domegady = SolnBlk.dWdy_faceE[ii][jj].omega;

     if (SolnBlk.Axisymmetric == AXISYMMETRIC_X) {
       radius = ( SolnBlk.Grid.xfaceE(ii,jj).x < MICRO) ? MICRO : SolnBlk.Grid.xfaceE(ii,jj).x;     
     } else if (SolnBlk.Axisymmetric == AXISYMMETRIC_Y) {    
       radius = ( SolnBlk.Grid.xfaceE(ii,jj).y < MICRO) ? MICRO : SolnBlk.Grid.xfaceE(ii,jj).y;            
     } 

     break;
    /****************************** SOUTH ******************************/
   case SOUTH:      
     QuadraturePoint_W = HALF*(SolnBlk.UnoSE(ii, jj).W() +SolnBlk.UnoSW(ii, jj).W() );
    Temp =  QuadraturePoint_W.T();     
     for(int Num = 0; Num<ns_values; Num++){
       h[Num] =  QuadraturePoint_W.specdata[Num].Enthalpy(Temp)+QuadraturePoint_W.specdata[Num].Heatofform();
       dcdx[Num] = SolnBlk.dWdx_faceS[ii][jj].spec[Num].c;  
       dcdy[Num] = SolnBlk.dWdy_faceS[ii][jj].spec[Num].c;
       dhdT[Num] = QuadraturePoint_W.specdata[Num].Enthalpy_prime(Temp);
     }

     drhodx = SolnBlk.dWdx_faceS[ii][jj].rho; 
     drhody = SolnBlk.dWdy_faceS[ii][jj].rho;
     dpdx = SolnBlk.dWdx_faceS[ii][jj].p;
     dpdy = SolnBlk.dWdy_faceS[ii][jj].p;
     dUdx = SolnBlk.dWdx_faceS[ii][jj].v.x;
     dUdy = SolnBlk.dWdy_faceS[ii][jj].v.x;
     dVdx = SolnBlk.dWdx_faceS[ii][jj].v.y;
     dVdy = SolnBlk.dWdy_faceS[ii][jj].v.y;
     dkdx = SolnBlk.dWdx_faceS[ii][jj].k;
     dkdy = SolnBlk.dWdy_faceS[ii][jj].k;
     domegadx = SolnBlk.dWdx_faceS[ii][jj].omega;
     domegady = SolnBlk.dWdy_faceS[ii][jj].omega;  
     
     if (SolnBlk.Axisymmetric == AXISYMMETRIC_X) {
       radius = (SolnBlk.Grid.xfaceS(ii,jj).x < MICRO) ? MICRO : SolnBlk.Grid.xfaceS(ii,jj).x;     
     } else if (SolnBlk.Axisymmetric == AXISYMMETRIC_Y) {    
       radius = (SolnBlk.Grid.xfaceS(ii,jj).y < MICRO) ? MICRO : SolnBlk.Grid.xfaceS(ii,jj).y;            
     } 
     break;
     /****************************** WEST ******************************/
   case WEST:       
     QuadraturePoint_W = HALF*(SolnBlk.UnoNW(ii, jj).W() +SolnBlk.UnoSW(ii, jj).W());
     Temp =  QuadraturePoint_W.T();
     for(int Num = 0; Num<ns_values; Num++){
       h[Num] =  QuadraturePoint_W.specdata[Num].Enthalpy(Temp)+QuadraturePoint_W.specdata[Num].Heatofform();
       dcdx[Num] = SolnBlk.dWdx_faceW[ii][jj].spec[Num].c;
       dcdy[Num] = SolnBlk.dWdy_faceW[ii][jj].spec[Num].c;
       dhdT[Num] = QuadraturePoint_W.specdata[Num].Enthalpy_prime(Temp);
     } 
    
     drhodx = SolnBlk.dWdx_faceW[ii][jj].rho; 
     drhody = SolnBlk.dWdy_faceW[ii][jj].rho;
     dpdx = SolnBlk.dWdx_faceW[ii][jj].p;
     dpdy = SolnBlk.dWdy_faceW[ii][jj].p;
     dUdx = SolnBlk.dWdx_faceW[ii][jj].v.x;
     dUdy = SolnBlk.dWdy_faceW[ii][jj].v.x;
     dVdx = SolnBlk.dWdx_faceW[ii][jj].v.y;
     dVdy = SolnBlk.dWdy_faceW[ii][jj].v.y;
     dkdx = SolnBlk.dWdx_faceW[ii][jj].k;
     dkdy = SolnBlk.dWdy_faceW[ii][jj].k;
     domegadx = SolnBlk.dWdx_faceW[ii][jj].omega;
     domegady = SolnBlk.dWdy_faceW[ii][jj].omega;
     
     if (SolnBlk.Axisymmetric == AXISYMMETRIC_X) {
       radius = (SolnBlk.Grid.xfaceW(ii,jj).x < MICRO) ? MICRO : SolnBlk.Grid.xfaceW(ii,jj).x;     
     } else if (SolnBlk.Axisymmetric == AXISYMMETRIC_Y) {    
       radius = (SolnBlk.Grid.xfaceW(ii,jj).y < MICRO) ? MICRO : SolnBlk.Grid.xfaceW(ii,jj).y;            
     } 
     break;
   }
   
   /////////////////////////////////////////////////////////////////////////////////

   rho = QuadraturePoint_W.rho;
   p = QuadraturePoint_W.p;
   kappa = QuadraturePoint_W.kappa();
   Cp =   QuadraturePoint_W.Cp();
   mu = QuadraturePoint_W.mu();
   Rmix =  QuadraturePoint_W.Rtot();
   
   /*********************** X - DIRECTION **************************************/

   dFvdWf(1, 7) += FOUR/THREE*mu;
   dFvdWf(1, 10) -= TWO/THREE*mu;
   dFvdWf(2, 8) += mu;
   dFvdWf(2, 9) += mu;  
   
   double Sum_q(ZERO);
   double Sum_dq(ZERO);
   
   for(int Num = 0; Num<ns_values; Num++){
      Sum_q +=  dhdT[Num]/(rho*Rmix)*mu/QuadraturePoint_W.Schmidt[Num]*dcdx[Num];
      Sum_dq -= dhdT[Num]*(mu/QuadraturePoint_W.Schmidt[Num])*dcdx[Num]*p/(rho*rho*Rmix);
   } 

   dFvdWf(3,0) += kappa*(-dpdx/(rho*rho*Rmix) +TWO*p*drhodx/(rho*rho*rho*Rmix))+Sum_dq;
   dFvdWf(3,1) += TWO*mu*(TWO/THREE*dUdx - dVdy/THREE);
   dFvdWf(3,2) += mu*(dUdy+dVdx);
   dFvdWf(3,3) += -drhodx/(rho*rho*Rmix)*kappa+Sum_q;

   dFvdWf(3,6) -= p/(rho*rho*Rmix)*kappa;
   dFvdWf(3,7) += FOUR/THREE*QuadraturePoint_W.v.x*mu;
   dFvdWf(3,8) += QuadraturePoint_W.v.y*mu;
   dFvdWf(3,9) = dFvdWf(3,8);
   dFvdWf(3,10) -= TWO/THREE*QuadraturePoint_W.v.x*mu;
   dFvdWf(3,11) = kappa/(rho*Rmix); 

   // Axisymmetric
   if(SolnBlk.Axisymmetric == AXISYMMETRIC_Y){
     dFvdWf(1,2) -=  TWO/THREE*mu/radius;
     dFvdWf(3,1) -=  TWO/THREE*mu*QuadraturePoint_W.v.y/radius;
     dFvdWf(3,2) -=  TWO/THREE*mu*QuadraturePoint_W.v.x/radius;
   }
   if(SolnBlk.Axisymmetric == AXISYMMETRIC_X){    
     dFvdWf(1,1) -=  TWO/THREE*mu/radius;
     dFvdWf(3,1) -=  FOUR/THREE*mu*QuadraturePoint_W.v.x/radius;
   }
 
   //multispecies
   for(int Num = 0; Num<(ns_species); Num++){
     dFvdWf(3,14+Num) = mu/QuadraturePoint_W.Schmidt[Num]*h[Num];  //+  H3???
     dFvdWf(6+Num, 14+Num) = mu/QuadraturePoint_W.Schmidt[Num];
   }

   //Turbulence
   if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
       SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_EPSILON) {
      double k = QuadraturePoint_W.k;
      double omega = QuadraturePoint_W.omega;
      double taux = (TWO/THREE*dUdx - dVdy/THREE);
      double tauy = dUdy+dVdx;
      kappa_t = QuadraturePoint_W.Kappa_turb();
      Dm_t = QuadraturePoint_W.Dm_turb();
      mu_t =  QuadraturePoint_W.eddy_viscosity();
      Pr_t =  QuadraturePoint_W.Pr_turb();
      Sc_t =  QuadraturePoint_W.Sc_turb();
      double sigma_s = QuadraturePoint_W.sigma_star;
      double sigma = QuadraturePoint_W.sigma;
     
      U= QuadraturePoint_W.v.x;
      V= QuadraturePoint_W.v.y;
      
      dFvdWf(1, 0) += TWO*k/max(omega,TOLER)*taux;
      dFvdWf(1, 4) += TWO*rho/max(omega,TOLER)*taux;
      dFvdWf(1, 5) += -TWO*rho*k/max(omega*omega,TOLER)*taux;
      dFvdWf(1, 7) += FOUR/THREE*mu_t;
      dFvdWf(1, 10) += -TWO/THREE*mu_t;

      dFvdWf(2, 0) += k/max(omega, TOLER)*tauy;
      dFvdWf(2, 4) += rho/max(omega, TOLER)*tauy; 
      dFvdWf(2, 5) += -mu_t/max(omega,TOLER)*tauy;
      dFvdWf(2, 8) += mu_t;
      dFvdWf(2, 9) += mu_t;

      double Sum_dh = 0.0;
      double Sum_h = 0.0;
    
       for(int Num = 0; Num<ns_values; Num++){
         Sum_h +=  h[Num]*dcdx[Num]/Sc_t;
         Sum_dh -= dhdT[Num]*dcdx[Num]/Sc_t;
      }
      
      dFvdWf(3,0) += k/max(omega,TOLER)*Cp/Pr_t*p*drhodx/(Rmix*rho*rho) 
	+ Temp/rho*k/max(omega,TOLER)*Sum_dh +k/max(omega,TOLER)*Sum_h 
         +TWO*QuadraturePoint_W.v.x*k/max(omega,TOLER)*taux+
         QuadraturePoint_W.v.y*k/max(omega,TOLER)*tauy+k/max(omega,TOLER)*sigma_s*dkdx;
      dFvdWf(3,1) += TWO*mu_t*taux;
      dFvdWf(3,2) += mu_t*tauy;
      dFvdWf(3,3) += -k*Cp*drhodx/(max(omega,TOLER)*Pr_t*Rmix*rho)- k/max(omega,TOLER)*Sum_dh/(Rmix);
      dFvdWf(3,4) += Cp*(dpdx - p*drhodx/rho)/(max(omega, TOLER) *Pr_t*Rmix) 
	+ rho/max(omega,TOLER)*Sum_h+TWO*rho*QuadraturePoint_W.v.x*taux/max(omega,TOLER)
         +QuadraturePoint_W.v.y*rho*tauy/max(omega,TOLER)+rho*sigma_s*dkdx/max(omega,TOLER);
      
      dFvdWf(3,5) += -k*Cp*(dpdx-p*drhodx/rho)/(max(omega*omega, TOLER)*Pr_t*Rmix)
	- Sum_h*mu_t/max(omega,TOLER)-TWO*QuadraturePoint_W.v.x*mu_t*taux/max(omega,TOLER)
         -QuadraturePoint_W.v.y*mu_t*tauy/max(omega,TOLER)-mu_t*sigma_s*dkdx/max(omega,TOLER);
      dFvdWf(3,6) += -k*Cp*Temp/(max(omega,TOLER)*Pr_t);
      
      dFvdWf(3,7) += FOUR/THREE*QuadraturePoint_W.v.x*mu_t;
      dFvdWf(3,8) += QuadraturePoint_W.v.y*mu_t;
      dFvdWf(3,9) += QuadraturePoint_W.v.y*mu_t;
      dFvdWf(3,10) += -TWO/THREE*QuadraturePoint_W.v.x*mu_t;
      dFvdWf(3,11) += k*Cp/(max(omega,TOLER)*Rmix*Pr_t);
      dFvdWf(3,12) += mu+mu_t*sigma_s;
      
      for(int Num = 0; Num<(ns_species); Num++){ 
         dFvdWf(3,14+Num) += h[Num]*mu_t/Sc_t;
         dFvdWf(6+Num,0) += dcdx[Num]*k/max(omega,TOLER)/Sc_t;
         dFvdWf(6+Num,4) += dcdx[Num]*rho/max(omega,TOLER)/Sc_t;
         dFvdWf(6+Num,5) += -mu_t*dcdx[Num]/max(omega,TOLER)/Sc_t;     
         dFvdWf(6+Num,14+Num) += mu_t/Sc_t;
      } 

      dFvdWf(4,0) += k*sigma_s*dkdx/max(omega,TOLER);
      dFvdWf(4,4) += rho*sigma_s*dkdx/max(omega,TOLER);
      dFvdWf(4,5) -=mu_t*sigma_s*dkdx/max(omega,TOLER); 
      dFvdWf(4,12) += mu+mu_t*sigma_s;
      
      dFvdWf(5,0) += k*sigma*domegadx/max(omega,TOLER);
      dFvdWf(5,4) += rho*sigma*domegadx/max(omega,TOLER);
      dFvdWf(5,5) += -mu_t*sigma*domegadx/max(omega,TOLER);
      dFvdWf(5,13) += mu+mu_t*sigma;

      if(SolnBlk.Axisymmetric == AXISYMMETRIC_Y){            
	dFvdWf(1,0) -= TWO/THREE*V*k/(max(omega,TOLER)*radius);
	dFvdWf(1,2) -= TWO/THREE*mu_t/radius;
	dFvdWf(1,4) -= TWO/THREE*rho*V/(max(omega,TOLER)*radius);
	dFvdWf(1,5) += TWO/THREE*mu_t*V/(max(omega,TOLER)*radius);	
	dFvdWf(3,0) -= TWO/THREE*U*k*V/(max(omega,TOLER)*radius);
	dFvdWf(3,1) -= TWO/THREE*mu_t*V/radius;
	dFvdWf(3,2) -= TWO/THREE*U*mu_t/radius;
	dFvdWf(3,4) -= TWO/THREE*U*V*rho/(max(omega,TOLER)*radius);
	dFvdWf(3,5) += TWO/THREE*U*V*mu_t/(max(omega,TOLER)*radius);                 
      }//endofaxisymmetric                 
   }//endof turbulence
   
   /*********************** Y - DIRECTION **************************************/

   dGvdWf(1, 8) += mu;
   dGvdWf(1, 9) += mu; 
   dGvdWf(2, 7) -= TWO/THREE*mu;
   dGvdWf(2, 10) += FOUR/THREE*mu;

   Sum_q = ZERO;  Sum_dq = ZERO;

   for(int Num = 0; Num<ns_values; Num++){
      Sum_q +=  dhdT[Num]/(rho*Rmix)*mu/QuadraturePoint_W.Schmidt[Num]*dcdy[Num];
      Sum_dq -= dhdT[Num]*(mu/QuadraturePoint_W.Schmidt[Num])*dcdy[Num]*p/(rho*rho*Rmix);
   } 
   
   dGvdWf(3,0) += kappa*(-dpdy/(rho*rho*Rmix) +TWO*p*drhody/(rho*rho*rho*Rmix))+Sum_dq; 
   dGvdWf(3,1) += mu*(dUdy +dVdx);
   dGvdWf(3,2) += TWO*mu*(TWO/THREE*dVdy-dUdx/THREE);
   dGvdWf(3,3) += -drhody/(rho*rho*Rmix)*kappa+Sum_q;
   dGvdWf(3,6) -= p/(rho*rho*Rmix)*kappa;
 
   dGvdWf(3,7) -= TWO/THREE*QuadraturePoint_W.v.y*mu;
   dGvdWf(3,8) += QuadraturePoint_W.v.x*mu;
   dGvdWf(3,9) = dGvdWf(3,8);
   dGvdWf(3,10) += FOUR/THREE*QuadraturePoint_W.v.y*mu;
   dGvdWf(3,11) = kappa/(rho*Rmix); 

   //Axisymmetric 
   if(SolnBlk.Axisymmetric == AXISYMMETRIC_Y){    
     dGvdWf(2,2) -=  TWO/THREE*mu/radius;
     dGvdWf(3,2) -=  FOUR/THREE*mu*QuadraturePoint_W.v.y/radius;
   }
   if(SolnBlk.Axisymmetric == AXISYMMETRIC_X){
     dGvdWf(2,1) -=  TWO/THREE*mu/radius;
     dGvdWf(3,1) -=  TWO/THREE*mu*QuadraturePoint_W.v.y/radius;
     dGvdWf(3,2) -=  TWO/THREE*mu*QuadraturePoint_W.v.x/radius;
   }

   //multispecies
   for(int Num = 0; Num<(ns_species); Num++){
      dGvdWf(3,14+Num) = mu/QuadraturePoint_W.Schmidt[Num]*h[Num];
      dGvdWf(6+Num, 14+Num) =  mu/QuadraturePoint_W.Schmidt[Num];
   }

   if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
       SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_EPSILON) {
      double k = QuadraturePoint_W.k;
      double omega = QuadraturePoint_W.omega;
      double taux = (TWO/THREE*dVdy - dUdx/THREE);
      double tauy = dUdy+dVdx;
      kappa_t = QuadraturePoint_W.Kappa_turb();
      Dm_t = QuadraturePoint_W.Dm_turb();
      mu_t =  QuadraturePoint_W.eddy_viscosity();
      Pr_t =  QuadraturePoint_W.Pr_turb();
      Sc_t =  QuadraturePoint_W.Sc_turb();
      double sigma_s = QuadraturePoint_W.sigma_star;
      double sigma = QuadraturePoint_W.sigma;
      
      U= QuadraturePoint_W.v.x;
      V= QuadraturePoint_W.v.y;

      dGvdWf(1, 0) += k/max(omega, TOLER)*tauy;
      dGvdWf(1, 4) += rho/max(omega, TOLER)*tauy;
      dGvdWf(1, 5) += -mu_t/max(omega,TOLER)*tauy;
      dGvdWf(1, 8) += mu_t;
      dGvdWf(1, 9) += mu_t;

      dGvdWf(2, 0) += TWO*k/max(omega,TOLER)*taux;
      dGvdWf(2, 4) += TWO*rho/max(omega,TOLER)*taux;
      dGvdWf(2, 5) += -TWO*rho*k/max(omega*omega,TOLER)*taux;
      dGvdWf(2, 7) += -TWO/THREE*mu_t;
      dGvdWf(2, 10) += FOUR/THREE*mu_t;
 
      double Sum_dh = 0.0;
      double Sum_h = 0.0;
      for(int Num = 0; Num<ns_values; Num++){
         Sum_dh  -=  dhdT[Num]*dcdy[Num]/Sc_t;
         Sum_h +=  h[Num]*dcdy[Num]/Sc_t;
      }

      dGvdWf(3,0) += k/max(omega,TOLER)*Cp/Pr_t*drhody*Temp/rho 
	+ Temp/rho*k/max(omega,TOLER)* Sum_dh+k/max(omega,TOLER)*Sum_h
	+TWO*QuadraturePoint_W.v.y*k/max(omega,TOLER)*taux+
	QuadraturePoint_W.v.x*k/max(omega,TOLER)*tauy+k/max(omega,TOLER)*sigma_s*dkdy;
      dGvdWf(3,1) += mu_t*tauy;
      dGvdWf(3,2) += TWO*mu_t*taux;
      dGvdWf(3,3) += -k/max(omega,TOLER)*Cp/Pr_t*drhody/(rho*Rmix) -k/max(omega,TOLER)*ONE/(Rmix)*Sum_dh;
      dGvdWf(3,4) += Cp*(dpdy - p/rho*drhody)/(max(omega, TOLER)*Pr_t*Rmix) + 
	rho*Sum_h/max(omega,TOLER)+TWO*rho*QuadraturePoint_W.v.y*taux/max(omega,TOLER) 
	+ QuadraturePoint_W.v.x*rho*tauy/max(omega,TOLER)+rho*sigma_s*dkdy/max(omega,TOLER);      
      dGvdWf(3,5) += -k*Cp*(dpdy-p/rho*drhody)/(max(omega*omega,TOLER)*Pr_t*Rmix) - 
	Sum_h*mu_t/max(omega,TOLER)-TWO*QuadraturePoint_W.v.y*mu_t*taux/max(omega,TOLER)
         -QuadraturePoint_W.v.x*mu_t*tauy/max(omega,TOLER)-mu_t*sigma_s*dkdy/max(omega,TOLER);
      dGvdWf(3,6) += -k/max(omega,TOLER)*Cp/Pr_t*Temp;
      
      dGvdWf(3,7) += -TWO/THREE*QuadraturePoint_W.v.y*mu_t;
      dGvdWf(3,8) += QuadraturePoint_W.v.x*mu_t;
      dGvdWf(3,9) +=QuadraturePoint_W.v.x*mu_t;
      dGvdWf(3,10) += FOUR/THREE*QuadraturePoint_W.v.y*mu_t;
      dGvdWf(3,11) += k*Cp/(max(omega,TOLER)*Pr_t*Rmix);
      dGvdWf(3,12) += mu+mu_t*sigma_s;
      
      for(int Num = 0; Num<(ns_species); Num++){ 
         dGvdWf(3,14+Num) += h[Num]*mu_t/Sc_t;
         dGvdWf(6+Num,0) +=k/max(omega,TOLER)*dcdy[Num]/Sc_t;  
         dGvdWf(6+Num,4) += dcdy[Num]*rho/max(omega,TOLER)/Sc_t;
         dGvdWf(6+Num,5) += -dcdy[Num]*mu_t/max(omega,TOLER)/Sc_t;     
         dGvdWf(6+Num,14+Num) += mu_t/Sc_t;
      } 

      dGvdWf(4,0) += k*sigma_s*dkdy/max(omega,TOLER);
      dGvdWf(4,4) += rho*sigma_s*dkdy/max(omega,TOLER);
      dGvdWf(4,5) -=  mu_t*sigma_s*dkdy/max(omega,TOLER); 
      dGvdWf(4,12) +=mu+mu_t*sigma_s;
      
      dGvdWf(5,0) += k*sigma*domegady/max(omega,TOLER);
      dGvdWf(5,4) += rho*sigma*domegady/max(omega,TOLER);
      dGvdWf(5,5) += -mu_t*sigma*domegady/max(omega,TOLER);
      dGvdWf(5,13) += mu+mu_t*sigma;

      if(SolnBlk.Axisymmetric == AXISYMMETRIC_Y){        
	dGvdWf(2,0) -=  TWO/THREE*V*k/(max(omega,TOLER)*radius);
	dGvdWf(2,2) -=  TWO/THREE*mu_t/radius;
	dGvdWf(2,4) -=  TWO/THREE*V*rho/(max(omega,TOLER)*radius);
	dGvdWf(2,5) +=  TWO/THREE*V*mu_t/(max(omega,TOLER)*radius);	
	dGvdWf(3,0) -=  TWO/THREE*V*V*k/(max(omega,TOLER)*radius);
	dGvdWf(3,2) -=  FOUR/THREE*mu_t*V/radius;
	dGvdWf(3,4) -=  TWO/THREE*V*V*rho/(max(omega,TOLER)*radius);
	dGvdWf(3,5) +=  TWO/THREE*V*V*mu_t/(max(omega,TOLER)*radius);         
      }//endofaxisymmetric      
   }//endof turbulence


   // Memory cleanup
   delete []h;   h = NULL;
   delete []dcdx;   dcdx = NULL;
   delete []dcdy;   dcdy = NULL;
   delete []dhdT;   dhdT = NULL;   
   
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
void dWfdWc_Diamond(DenseMatrix &dWfdWc_x,DenseMatrix &dWfdWc_y, Chem2D_Quad_Block &SolnBlk,
		    const int &Orient_face, const int &i, const int &j, const int &Orient_cell){

   int ns = SolnBlk.W[i][j].ns-1;
     
   int Left, Right; 
   double LL(ZERO),RR(ZERO); 
   double d_dWdx_dW(ZERO), d_dWdy_dW(ZERO);
   
   // All these confusing relationships are based on an outward
   // facing normal from the Orient_face ie  NORTH face ,
   // left is NW and right is NE.  

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

/********************************************************
 * Routine: d_dWd_dW_Diamond                            *
 *                                                      *
 * This routine calculates the 2nd deriavaites          *
 * associated with diamond path and bilinear            *
 * interpolation.                                       *
 *                                                      *
 ********************************************************/
void d_dWd_dW_Diamond(double &d_dWdx_dW, double &d_dWdy_dW, Chem2D_Quad_Block &SolnBlk, 
		      const double &LEFT, const double &RIGHT, const int &Orient_cell,
		      const int &Orient_face,  const int &i, const int &j){

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

/*************** CENTER ****************************/
void d_dWd_dW_Center(double &d_dWdx_dW_C, double &d_dWdy_dW_C, 
		     Chem2D_Quad_Block &SolnBlk, 
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

/********************************************************
 * Routine:  dS_tdW                                     *
 *                                                      *
 * This routine calculates the Turbulence (k-omeaga)    *
 * source Jacobian.                                     *
 *                                                      *
 ********************************************************/
int dS_tdW(DenseMatrix &dStdW,  Chem2D_Quad_Block &SolnBlk,
	   double &d_dWdx_dW, double &d_dWdy_dW,
	   const int &ii, const int &jj){
  
   double rho, U,V, k, omega;
   double dUdx, dUdy, dVdx, dVdy;
   double mu_t, alpha, beta_star, beta;
   
   beta_star = SolnBlk.W[ii][jj].beta_star;
   beta = SolnBlk.W[ii][jj].beta;
   alpha = SolnBlk.W[ii][jj].alpha;
   mu_t = SolnBlk.W[ii][jj].eddy_viscosity();
   
   rho = SolnBlk.W[ii][jj].rho;
   U = SolnBlk.W[ii][jj].v.x;
   V = SolnBlk.W[ii][jj].v.y;
   k = SolnBlk.W[ii][jj].k;
   omega = SolnBlk.W[ii][jj].omega;
   
   dUdx =  SolnBlk.dWdx[ii][jj].v.x;
   dUdy =  SolnBlk.dWdy[ii][jj].v.x;
   dVdx =  SolnBlk.dWdx[ii][jj].v.y;
   dVdy =  SolnBlk.dWdy[ii][jj].v.y;

//    d_dWdx_dW =  SolnBlk.d_dWdx_dW[ii][jj][CENTER];
//    d_dWdy_dW =  SolnBlk.d_dWdy_dW[ii][jj][CENTER];

   double tau_x = TWO/THREE*dUdx - dVdy/THREE;
   double tau_y = TWO/THREE*dVdy-dUdx/THREE;
   double tau_z = dUdy+dVdx;
   
   dStdW(4,0) += (TWO*k/max(omega,TOLER)*tau_x - TWO/THREE*k)*dUdx +k/max(omega,TOLER)
      *tau_z*tau_z + (TWO*k/max(omega,TOLER)*tau_y - TWO/THREE*k)*dVdy - beta_star*k*omega;
   
   dStdW(4,1) += FOUR/THREE*mu_t*d_dWdx_dW*dUdx+(TWO*mu_t*tau_x - TWO/THREE*rho*k)*d_dWdx_dW
      + TWO*mu_t*tau_z*d_dWdy_dW - TWO/THREE*mu_t*d_dWdx_dW*dVdy;
   
   dStdW(4,2) += -TWO/THREE*mu_t*d_dWdy_dW*dUdx+TWO*mu_t*tau_z*d_dWdx_dW+FOUR/THREE*
      mu_t*d_dWdy_dW*dVdy+(TWO*mu_t*tau_y-TWO/THREE*rho*k)*d_dWdy_dW;
   
   dStdW(4,4) += (TWO*rho/max(omega,TOLER)*tau_x - TWO/THREE*rho)*dUdx +rho/max(omega,TOLER)*
      tau_z*tau_z + (TWO*rho/max(omega,TOLER)*tau_y - TWO/THREE*rho)*dVdy - beta_star*rho*omega;

   dStdW(4,5) += -TWO*mu_t/max(omega,TOLER)*tau_x*dUdx - mu_t/max(omega,TOLER)*tau_z*tau_z -
      TWO*mu_t/max(omega,TOLER)*tau_y*dVdy - beta_star*rho*k;
  
   dStdW(5,0) +=alpha*omega/max(k,TOLER)*((TWO*k/max(omega,TOLER)*tau_x - TWO/THREE*k)*dUdx +k/max(omega,TOLER)*tau_z*tau_z
                                         + (TWO*k/max(omega,TOLER)*tau_y - TWO/THREE*k)*dVdy) - beta*omega*omega;
   
   dStdW(5,1) +=alpha*omega/max(k,TOLER)*(FOUR/THREE*mu_t*d_dWdx_dW*dUdx+( TWO*mu_t*tau_x - TWO/THREE*rho*k)*d_dWdx_dW +
                                         TWO*mu_t*tau_z*d_dWdy_dW - TWO/THREE*mu_t*d_dWdx_dW*dVdy);
   
   dStdW(5,2) +=alpha*omega/max(k,TOLER)*(-TWO/THREE*mu_t*d_dWdy_dW*dUdx + TWO*mu_t*tau_z*d_dWdx_dW + 
      FOUR/THREE*mu_t*d_dWdy_dW*dVdy + (TWO*mu_t*tau_y - TWO/THREE*rho*k)*d_dWdy_dW  );
   
   dStdW(5,4) += -alpha*omega/max(k*k,TOLER)*((TWO*mu_t*tau_x - TWO/THREE*rho*k )*dUdx + mu_t*tau_z*tau_z 
                                             +(TWO*mu_t*tau_y - TWO/THREE*rho*k)*dVdy) +
      alpha*omega/max(k,TOLER)*((TWO*rho/max(omega,TOLER)*tau_x - TWO/THREE*rho)*dUdx+rho/max(omega,TOLER)*tau_z*tau_z 
       + (TWO*rho/max(omega,TOLER)*tau_y-TWO/THREE*rho)*dVdy);
   
   dStdW(5,5) +=  alpha/max(k,TOLER)*((TWO*mu_t*tau_x - TWO/THREE*rho*k)*dUdx + mu_t*tau_z*tau_z+(TWO*mu_t*tau_y-TWO*rho*k)*dVdy) +
      alpha/max(k,TOLER)*(-TWO*mu_t*tau_x*dUdx - mu_t*tau_z*tau_z - TWO*mu_t*tau_y*dVdy) - TWO*beta*rho*omega;


   if(SolnBlk.Axisymmetric ==AXISYMMETRIC_Y){
      double radius  = SolnBlk.Grid.Cell[ii][jj].Xc.y;
      if(radius !=ZERO){
         double tau_t = TWO/THREE*V/radius - dUdx/THREE - dVdy/THREE;
        
         dStdW(4,0) += -TWO/THREE*V*k/max(omega,TOLER)/radius*(dUdx +dVdy) + (TWO*k*tau_t/max(omega,TOLER) - TWO/THREE*k)*V/radius;
         dStdW(4,1) -=FOUR/THREE*mu_t*V*d_dWdx_dW/radius;
         dStdW(4,2) += - TWO/THREE*mu_t*(dUdx+dVdy)/radius - TWO/THREE*mu_t*V*d_dWdy_dW/radius+TWO*mu_t*(TWO/THREE/radius - d_dWdy_dW/THREE)*V/radius
	   + (TWO*mu_t*tau_t - TWO/THREE*rho*k)/radius;
         dStdW(4,4) += -TWO/THREE*rho*V/max(omega,TOLER)/radius*(dUdx+dVdy) + V*(TWO*rho/max(omega,TOLER)*tau_t - TWO/THREE*rho)/radius;
         dStdW(4,5) += TWO/THREE*mu_t*V/max(omega,TOLER)/radius*(dUdx+dVdy) - TWO*mu_t*V*tau_t/max(omega,TOLER)/radius;         
         dStdW(5,0) += alpha*omega/max(TOLER,k)*(-TWO/THREE*k*V*(dUdx+dVdy)/max(omega,TOLER)/radius + (TWO*k*tau_t/max(omega,TOLER) - TWO/THREE*k)*V/radius);
         dStdW(5,1) -= FOUR/THREE*alpha*rho*V*d_dWdx_dW/radius;
         dStdW(5,2) += alpha*omega/max(k,TOLER)*(-TWO/THREE*mu_t*(dUdx+dVdy)/radius-TWO/THREE*mu_t*V*d_dWdy_dW/radius
                                                 +TWO*mu_t*(TWO/THREE/radius - d_dWdy_dW/THREE)*V/radius + (TWO*mu_t*tau_t - TWO/THREE*rho*k)/radius);
	 dStdW(5,4) += - alpha*omega/max(k*k,TOLER)*(-TWO/THREE*mu_t*V*(dUdx + dVdy)/radius + (TWO*mu_t*tau_t - TWO/THREE*rho*k)*V/radius)+
	   alpha*omega/max(TOLER,k)*(-TWO/THREE*V*rho*(dUdx+dVdy)/max(omega,TOLER)/radius +(TWO*rho*tau_t/max(omega,TOLER) - TWO/THREE*rho)*V/radius );         
         dStdW(5,5) += alpha/max(k,TOLER)*(-TWO/THREE*mu_t*V*(dUdx+dVdy)/radius + (TWO*mu_t*tau_t - TWO/THREE*rho*k)*V/radius )+
            alpha/max(k,TOLER)*(TWO/THREE*mu_t*V*(dUdx+dVdy)/radius - TWO*mu_t*V*tau_t/radius);
      }      
   }
   
   return (0);
}


// ARE THES FUNCTIONS EVER CALLED ????


// int Automatic_Wall_Treatment_Residual_Jacobian(Chem2D_Quad_Block &SolnBlk, Chem2D_Input_Parameters &Input_Parameters, int i, int j, DenseMatrix &dRdU){
   
//    int NUM_VAR = SolnBlk.NumVar()-1;  

//     //First cells off wall  
//    if(((i==SolnBlk.ICl) && ((SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS ||
//                              SolnBlk.Grid.BCtypeW[j] == BC_NO_SLIP  ||
//                              SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
//                              SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX ||
//                              SolnBlk.Grid.BCtypeW[j] == BC_ADIABATIC_WALL))  ) 
//       ||((i==SolnBlk.ICu) &&(SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS ||
//                              SolnBlk.Grid.BCtypeE[j] == BC_NO_SLIP  ||
//                              SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
//                              SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX ||
//                              SolnBlk.Grid.BCtypeE[j] == BC_ADIABATIC_WALL))
//       || ((j == SolnBlk.JCl) &&((SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS ||
//                                  SolnBlk.Grid.BCtypeS[i] == BC_NO_SLIP  ||
//                                  SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
//                                  SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX ||
//                                  SolnBlk.Grid.BCtypeS[i] == BC_ADIABATIC_WALL))  )
//       || ((j ==SolnBlk.JCu) &&((SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS ||
//                                 SolnBlk.Grid.BCtypeN[i] == BC_NO_SLIP  ||
//                                 SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
//                                 SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX ||
//                                 SolnBlk.Grid.BCtypeN[i] == BC_ADIABATIC_WALL)))){

//       if(SolnBlk.Wall[i][j].yplus <=SolnBlk.W[i][j].y_sublayer){
//          for(int jcol=0; jcol<NUM_VAR; jcol++){
//             if(jcol!=5){
//                dRdU(5,jcol) = 0.0;
//             }
//          }
         
            
//       }
      
//       /* y+ > 2.5 && y+<30  apply  blending formulation */
//       if(SolnBlk.Wall[i][j].yplus >SolnBlk.W[i][j].y_sublayer
//          && SolnBlk.Wall[i][j].yplus <SolnBlk.W[i][j].Yplus_l){
//          for(int jcol=0; jcol<NUM_VAR; jcol++){
//             if(jcol!=4){
//                dRdU(4,jcol) = 0.0;
//             }
//          }
//          for(int jcol=0; jcol<NUM_VAR; jcol++){
//             if(jcol!=5){
//                dRdU(5,jcol) = 0.0;
//             }
//          }
//       }      
//       /* y+ >= 30  apply wall function */
//       if(SolnBlk.Wall[i][j].yplus >= SolnBlk.W[i][j].Yplus_l){
//          //cout<<"-----------#3"<<endl;
//          // Set k
//          for(int jcol=0; jcol<NUM_VAR; jcol++){
//             if(jcol!=4){
//                dRdU(4,jcol) = 0.0;
//             }
//          }
//          for(int jcol=0; jcol<NUM_VAR; jcol++){
//             if(jcol!=5){
//                dRdU(5,jcol) = 0.0;
//             }
                    
//          }
//       }
      
//     // first cells off walls
//    }else if ((SolnBlk.Wall[i][j].yplus >SolnBlk.W[i][j].y_sublayer
//               && SolnBlk.Wall[i][j].yplus <SolnBlk.W[i][j].Yplus_l) &&(
//                  (((i-1 ==SolnBlk.ICl) && ((SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS ||
//                                             SolnBlk.Grid.BCtypeW[j] == BC_NO_SLIP  ||
//                                             SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
//                                             SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX ||
//                                             SolnBlk.Grid.BCtypeW[j] == BC_ADIABATIC_WALL))  ) 
//                   ||((i+1==SolnBlk.ICu) &&(SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS ||
//                                            SolnBlk.Grid.BCtypeE[j] == BC_NO_SLIP  ||
//                                            SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
//                                            SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX ||
//                                            SolnBlk.Grid.BCtypeE[j] == BC_ADIABATIC_WALL))
//                   || ((j-1 == SolnBlk.JCl) &&((SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS ||
//                                                SolnBlk.Grid.BCtypeS[i] == BC_NO_SLIP  ||
//                                                SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
//                                                SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX ||
//                                                SolnBlk.Grid.BCtypeS[i] == BC_ADIABATIC_WALL))  )
//                   || ((j+1 ==SolnBlk.JCu) &&((SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS ||
//                                               SolnBlk.Grid.BCtypeN[i] == BC_NO_SLIP  ||
//                                               SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
//                                               SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX ||
//                                               SolnBlk.Grid.BCtypeN[i] == BC_ADIABATIC_WALL)))) ||
//                  (((i-2 ==SolnBlk.ICl) && ((SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS ||
//                                             SolnBlk.Grid.BCtypeW[j] == BC_NO_SLIP  ||
//                                             SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
//                                             SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX ||
//                                             SolnBlk.Grid.BCtypeW[j] == BC_ADIABATIC_WALL))  ) 
//                   ||((i+2==SolnBlk.ICu) &&(SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS ||
//                                            SolnBlk.Grid.BCtypeE[j] == BC_NO_SLIP  ||
//                                            SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
//                                            SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX ||
//                                            SolnBlk.Grid.BCtypeE[j] == BC_ADIABATIC_WALL))
//                   || ((j-2 == SolnBlk.JCl) &&((SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS ||
//                                                SolnBlk.Grid.BCtypeS[i] == BC_NO_SLIP  ||
//                                                SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
//                                                SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX ||
//                                                SolnBlk.Grid.BCtypeS[i] == BC_ADIABATIC_WALL))  )
//                   || ((j+2 ==SolnBlk.JCu) &&((SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS ||
//                                               SolnBlk.Grid.BCtypeN[i] == BC_NO_SLIP  ||
//                                               SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
//                                               SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX ||
//                                               SolnBlk.Grid.BCtypeN[i] == BC_ADIABATIC_WALL)))))) {
//       for(int jcol=0; jcol<NUM_VAR; jcol++){
//          if(jcol!=4){
//             dRdU(4,jcol) = 0.0;
//          }
//       }
//       for(int jcol=0; jcol<NUM_VAR; jcol++){
//          if(jcol!=5){
//             dRdU(5,jcol) = 0.0;
//          }         
//       }
      
//    } else if (SolnBlk.Wall[i][j].yplus <=SolnBlk.W[i][j].y_sublayer){
   
//       for(int jcol=0; jcol<NUM_VAR; jcol++){
//          if(jcol!=5){
//             dRdU(5,jcol) = 0.0;
//          }        
//       }     
//    }
   
//    return 0;

// }

// int Wall_Function_Residual_Jacobian(Chem2D_Quad_Block &SolnBlk, 
//                          Chem2D_Input_Parameters &Input_Parameters, 
//                          int i, int j, DenseMatrix &dRdU){
   
//   int NUM_VAR = SolnBlk.NumVar()-1;

//   if (SolnBlk.Wall[i][j].yplus <= SolnBlk.W[i][j].Yplus_u) {                    
//     for(int jcol=0; jcol<NUM_VAR; jcol++){
//       if(jcol!=4){
// 	dRdU(4,jcol) = 0.0;
//       }
//     }
//     for(int jcol=0; jcol<NUM_VAR; jcol++){
//       if(jcol!=5){
// 	dRdU(5,jcol) = 0.0;
//       }      
//     }    
//    } /* endif */
  
//   return 0;   
// }


// void Low_Reynoldsnumber_Formulation_Residual_Jacobian(Chem2D_Quad_Block &SolnBlk, 
// 						      Chem2D_Input_Parameters &Input_Parameters, 
// 						      int i, int j, DenseMatrix & dRdU){
  
//   int NUM_VAR = SolnBlk.NumVar()-1;
  
//   if (( SolnBlk.Wall[i][j].yplus <=SolnBlk.W[i][j].y_sublayer)
//       ||((i==SolnBlk.ICl) && ( (SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS ||
// 				SolnBlk.Grid.BCtypeW[j] == BC_NO_SLIP  ||
// 				SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
// 				SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX ||
// 				SolnBlk.Grid.BCtypeW[j] == BC_ADIABATIC_WALL))  ) 
//       ||  ((i==SolnBlk.ICu) &&(SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS ||
// 			       SolnBlk.Grid.BCtypeE[j] == BC_NO_SLIP  ||
// 			       SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
// 			       SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX ||
// 			       SolnBlk.Grid.BCtypeE[j] == BC_ADIABATIC_WALL))
//       || ((j == SolnBlk.JCl) &&((SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS ||
//                                  SolnBlk.Grid.BCtypeS[i] == BC_NO_SLIP  ||
//                                  SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
//                                  SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX ||
//                                  SolnBlk.Grid.BCtypeS[i] == BC_ADIABATIC_WALL))  )
//       || ((j ==SolnBlk.JCu) &&((SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS ||
// 				SolnBlk.Grid.BCtypeN[i] == BC_NO_SLIP  ||
// 				SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
// 				SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX ||
// 				SolnBlk.Grid.BCtypeN[i] == BC_ADIABATIC_WALL)))){
    
//     for(int jcol=0; jcol<NUM_VAR; jcol++){
//       if(jcol!=5){
// 	dRdU(5,jcol) = 0.0;
//       }      
//     }    
//   }   

// }


// void BC_Residual_Jacobian(Chem2D_Quad_Block &SolnBlk, Chem2D_Input_Parameters &Input_Parameters, int i, int j, DenseMatrix &dRdU){
   
//   if (((i==SolnBlk.ICl) && ( (SolnBlk.Grid.BCtypeW[j] != BC_NONE ))) 
//       ||  ((i==SolnBlk.ICu) &&(SolnBlk.Grid.BCtypeE[j] != BC_NONE))
//       || ((j == SolnBlk.JCl) &&((SolnBlk.Grid.BCtypeS[i] != BC_NONE )))
//       || ((j ==SolnBlk.JCu) &&((SolnBlk.Grid.BCtypeN[i] != BC_NONE )))){
          
//     dRdU.zero();
//     int NUM_VAR_CHEM2D = SolnBlk.NumVar();
    
//     for(int irow=0; irow<(NUM_VAR_CHEM2D-1); irow++)
//       for(int jcol=0; jcol<(NUM_VAR_CHEM2D-1); jcol++){
// 	if(irow==jcol){
// 	  dRdU(irow, jcol) = -ONE/SolnBlk.dt[i][j];
// 	}
//       }    
//   }
// }



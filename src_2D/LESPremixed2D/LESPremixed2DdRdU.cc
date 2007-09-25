#ifndef _LESPREMIXED2D_dRdU_INCLUDED
#include "LESPremixed2DdRdU.h"
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
			      LESPremixed2D_Quad_Block &SolnBlk,
			      LESPremixed2D_Input_Parameters &Input_Parameters,
			      const int &ii, const int &jj){
   
#ifdef THICKENED_FLAME_ON
  int NUM_VAR_LESPREMIXED2D =  SolnBlk.NumVar()-3; 
#else
  int NUM_VAR_LESPREMIXED2D =  SolnBlk.NumVar()-1; 
#endif
  DenseMatrix dRdW(NUM_VAR_LESPREMIXED2D,NUM_VAR_LESPREMIXED2D,ZERO);  

  //Inviscid dRdU
  dFIdW_Inviscid(dRdW, SolnBlk, Input_Parameters, ii,jj);
  
  //Viscous dRdU
  if (SolnBlk.Flow_Type != FLOWTYPE_INVISCID) {    
    dGVdW_Viscous(dRdW, SolnBlk,Input_Parameters, ii, jj); 
  }

  // Add Source Jacobians (axisymmetric, turbulence)
  SemiImplicitBlockJacobi_dSdW(dRdW,SolnBlk,EXPLICIT,ii,jj);                          

  DenseMatrix dWdU(NUM_VAR_LESPREMIXED2D,NUM_VAR_LESPREMIXED2D,ZERO);     
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
			     LESPremixed2D_Quad_Block &SolnBlk,
			     const int &solver_type,
			     const int &ii, const int &jj){ 
  
  if( (SolnBlk.Axisymmetric && SolnBlk.Flow_Type != FLOWTYPE_INVISCID) ||
       SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_LES_TF_K ||
       SolnBlk.Flow_Type == FLOWTYPE_LAMINAR_C || 
       SolnBlk.Flow_Type == FLOWTYPE_LAMINAR_C_ALGEBRAIC || 
       SolnBlk.Flow_Type == FLOWTYPE_LAMINAR_C_FSD || 
       SolnBlk.Flow_Type == FLOWTYPE_LAMINAR_NGT_C_FSD || 
       SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_LES_C || 
       SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_LES_C_ALGEBRAIC || 
       SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_SMAGORINSKY || 
       SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_CHARLETTE || 
       SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_LES_NGT_C_FSD_SMAGORINSKY || 
       SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_K ||
       SolnBlk.Flow_Type == FLOWTYPE_FROZEN_TURBULENT_LES_C_FSD) { 

#ifdef THICKENED_FLAME_ON
    int NUM_VAR_LESPREMIXED2D =  SolnBlk.NumVar()-3; 
#else
    int NUM_VAR_LESPREMIXED2D =  SolnBlk.NumVar()-1;
#endif

   if ( SolnBlk.Flow_Type == FLOWTYPE_LAMINAR_C || 
	SolnBlk.Flow_Type == FLOWTYPE_LAMINAR_C_ALGEBRAIC || 
	SolnBlk.Flow_Type == FLOWTYPE_LAMINAR_C_FSD || 
	SolnBlk.Flow_Type == FLOWTYPE_LAMINAR_NGT_C_FSD || 
	SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_LES_C || 
	SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_LES_C_ALGEBRAIC || 
	SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_SMAGORINSKY || 
	SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_CHARLETTE || 
	SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_LES_NGT_C_FSD_SMAGORINSKY || 
	SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_K ||
	SolnBlk.Flow_Type == FLOWTYPE_FROZEN_TURBULENT_LES_C_FSD ){
      NUM_VAR_LESPREMIXED2D =  SolnBlk.NumVar()-SolnBlk.W[0][0].ns;
    }

    DenseMatrix dRdW(NUM_VAR_LESPREMIXED2D,NUM_VAR_LESPREMIXED2D,ZERO);  
    
    // Add Source Jacobians (viscous axisymmetric, turbulence)
    SemiImplicitBlockJacobi_dSdW(dRdW,SolnBlk,EXPLICIT,ii,jj);                          
    
    DenseMatrix dWdU(NUM_VAR_LESPREMIXED2D,NUM_VAR_LESPREMIXED2D,ZERO);     
    // Transformation Jacobian 
    SolnBlk.Uo[ii][jj].W().dWdU(dWdU, SolnBlk.Flow_Type); 
    dSdU += dRdW*dWdU;
  }

  // Add Source Jacobians (inviscid axisymmetric, chemistry, gravity)
  SemiImplicitBlockJacobi_dSdU(dSdU,SolnBlk,EXPLICIT,ii,jj);                 

}

void SemiImplicitBlockJacobi_dSdW(DenseMatrix &dSdW,
				  LESPremixed2D_Quad_Block &SolnBlk,
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


  /////////////////////////////////////////////////////////////////
  //
  // dS_tdW should be added here for the k-equation and Wen's stuff
  //
  ////////////////////////////////////////////////////////////////

   if ( SolnBlk.Flow_Type == FLOWTYPE_LAMINAR_C || 
	SolnBlk.Flow_Type == FLOWTYPE_LAMINAR_C_ALGEBRAIC || 
	SolnBlk.Flow_Type == FLOWTYPE_LAMINAR_C_FSD || 
	SolnBlk.Flow_Type == FLOWTYPE_LAMINAR_NGT_C_FSD || 
	SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_LES_C || 
	SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_LES_C_ALGEBRAIC || 
	SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_SMAGORINSKY || 
	SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_CHARLETTE || 
	SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_LES_NGT_C_FSD_SMAGORINSKY || 
	SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_K ||
	SolnBlk.Flow_Type == FLOWTYPE_FROZEN_TURBULENT_LES_C_FSD ){
     dS_tdW(dSdW,SolnBlk, d_dWdx_dW_C, d_dWdy_dW_C, ii,jj);
   }
  
  // Add Jacobian for turbulence k-omega
  // if((SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ) ||
//      (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_EPSILON )){    
//     dS_tdW(dSdW,SolnBlk, d_dWdx_dW_C, d_dWdy_dW_C, ii,jj);
//   }
  
}

void SemiImplicitBlockJacobi_dSdU(DenseMatrix &dSdU,
				  LESPremixed2D_Quad_Block &SolnBlk,
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
  if (SolnBlk.W[ii][jj].React.reactset_flag != NO_REACTIONS &&
      SolnBlk.Flow_Type != FLOWTYPE_LAMINAR_C &&
      SolnBlk.Flow_Type != FLOWTYPE_LAMINAR_C_ALGEBRAIC && 
      SolnBlk.Flow_Type != FLOWTYPE_LAMINAR_C_FSD &&
      SolnBlk.Flow_Type != FLOWTYPE_LAMINAR_NGT_C_FSD &&
      SolnBlk.Flow_Type != FLOWTYPE_TURBULENT_LES_C && 
      SolnBlk.Flow_Type != FLOWTYPE_TURBULENT_LES_C_ALGEBRAIC && 
      SolnBlk.Flow_Type != FLOWTYPE_TURBULENT_LES_C_FSD_SMAGORINSKY &&
      SolnBlk.Flow_Type != FLOWTYPE_TURBULENT_LES_C_FSD_CHARLETTE &&
      SolnBlk.Flow_Type != FLOWTYPE_TURBULENT_LES_NGT_C_FSD_SMAGORINSKY &&
      SolnBlk.Flow_Type != FLOWTYPE_TURBULENT_LES_C_FSD_K &&
      SolnBlk.Flow_Type != FLOWTYPE_FROZEN_TURBULENT_LES_C_FSD){    
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
void dFIdW_Inviscid(DenseMatrix &dRdW, LESPremixed2D_Quad_Block &SolnBlk, LESPremixed2D_Input_Parameters 
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
 * Routine: Rotation_matrix2                             *
 *                                                       *
 * This function returns either the rotation matrix, A,  *
 * or the inverse of A.                                  *
 *                                                       *
 * Note: A_matrix = 1 for returning A.                   *
 *                = 0 for returning inverse of A.        *
 *                                                       *
 *********************************************************/
/*The rotation matrix is used for the inviscid flux calculations */
DenseMatrix Rotation_matrix2(Vector2D nface, int Size,  int A_matrix) 
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
 * Routine: Inviscid Flux Jacobian using HLLE           *
 *                                                      *
 * This routine returns the inviscid components of      *    
 * Point Implicit Block Jacobian matrix for the         *
 * specified local solution block.                      *
 *                                                      *
 ********************************************************/
void dFIdW_Inviscid_HLLE(DenseMatrix &dRdW, LESPremixed2D_Quad_Block &SolnBlk,
			 LESPremixed2D_Input_Parameters &Input_Parameters, 
			 const int &ii, const int &jj, const int Orient){

#ifdef THICKENED_FLAME_ON
  int NUM_VAR_LESPREMIXED2D =  SolnBlk.NumVar()-2; 
#else
  int NUM_VAR_LESPREMIXED2D =  SolnBlk.NumVar(); 
#endif
  int overlap = Input_Parameters.NKS_IP.GMRES_Overlap;

   DenseMatrix dFidW(NUM_VAR_LESPREMIXED2D, NUM_VAR_LESPREMIXED2D,ZERO);    
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

   DenseMatrix A( Rotation_matrix2(nface,NUM_VAR_LESPREMIXED2D, 1));
   DenseMatrix AI( Rotation_matrix2(nface,NUM_VAR_LESPREMIXED2D, 0));
   DenseMatrix II(NUM_VAR_LESPREMIXED2D, NUM_VAR_LESPREMIXED2D,ZERO); 
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
void dFIdW_Inviscid_ROE(DenseMatrix& dRdW, LESPremixed2D_Quad_Block &SolnBlk,  
			LESPremixed2D_Input_Parameters &Input_Parameters,
			const int &ii, const int &jj, const int Orient){
   
  int overlap = Input_Parameters.NKS_IP.GMRES_Overlap;
  int Ri, Rj;

  if (ii < SolnBlk.ICl -overlap || ii > SolnBlk.ICu + overlap ||
      jj < SolnBlk.JCl -overlap || jj > SolnBlk.JCu + overlap) {

     // GHOST CELL so do nothing
     cout<<"\n Hey I am not suppose to be here! \n"; exit(1);

  } else {     
    int NUM_VAR_LESPREMIXED2D = dRdW.get_n();  //  SolnBlk.NumVar()-1;    
    DenseMatrix dFidW(NUM_VAR_LESPREMIXED2D, NUM_VAR_LESPREMIXED2D,ZERO);
    
    Vector2D nface,DX; double lface;   
    LESPremixed2D_pState Wa, wavespeeds, Left, Right, Wl, Wr;   
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
     
     DenseMatrix A( Rotation_matrix2(nface,NUM_VAR_LESPREMIXED2D, 1));
     DenseMatrix AI( Rotation_matrix2(nface,NUM_VAR_LESPREMIXED2D, 0));
     
     //Left and Right States                                                       ///Ri ,Rj fixed not ii, jj 
     Inviscid_Flux_Used_Reconstructed_LeftandRight_States(Wl, Wr, DX, SolnBlk, Orient, Ri, Rj );   
     Left  = Rotate(Wl, nface);
     Right = Rotate(Wr, nface);
     
     //Determin Roe Averaged State
     Wa = RoeAverage(Left,Right);       

     // Jacobian dF/dW         
     dFIdW(dFidW, Rotate(SolnBlk.Uo[ii][jj].W(), nface) , SolnBlk.Flow_Type);       
     dFidW = HALF*dFidW;
     
     /***************************** Regular Roe (no preconditioning) *************************************/
     if(!Input_Parameters.Preconditioning){
       // Determine Wave Speeds
       wavespeeds = HartenFixAbs( Wa.lambda_x(),
				  Left.lambda_x(),
				  Right.lambda_x());         
       
       //Loop through each wavespeed and each element of Jacobian(i,j)        
       for (int i=1; i <= NUM_VAR_LESPREMIXED2D; i++) {		   
	 for(int irow =0; irow< NUM_VAR_LESPREMIXED2D; irow++){
	   for(int jcol =0; jcol< NUM_VAR_LESPREMIXED2D; jcol++){

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
				  Wl.lambda_preconditioned_x(Wl.Mr2(SolnBlk.Flow_Type,deltax)),
				  Wr.lambda_preconditioned_x(Wr.Mr2(SolnBlk.Flow_Type,deltax)));
       
       
       //Calculate the preconditioned upwind dissipation flux.
       DenseMatrix Flux_dissipation(NUM_VAR_LESPREMIXED2D,NUM_VAR_LESPREMIXED2D,ZERO);        
       for (int i=1; i <=  NUM_VAR_LESPREMIXED2D; i++) {		   
	 for(int irow =0; irow < NUM_VAR_LESPREMIXED2D; irow++){
	   for(int jcol =0; jcol < NUM_VAR_LESPREMIXED2D; jcol++){	   
	     Flux_dissipation(irow, jcol) -= HALF*wavespeeds[i]*Wa.lp_x_precon(i,MR2a)[jcol+1]*Wa.rc_x_precon(i,MR2a)[irow+1];   
	   }
	 }
       }
       
       // Evaluate the low-Mach-number local preconditioner for the Roe-averaged state.
       DenseMatrix P( NUM_VAR_LESPREMIXED2D, NUM_VAR_LESPREMIXED2D,ZERO);          
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
void dFIdW_Inviscid_ROE_FD(DenseMatrix& dRdW, LESPremixed2D_Quad_Block &SolnBlk,  
			   LESPremixed2D_Input_Parameters &Input_Parameters,
			   const int &ii, const int &jj, const int Orient){
   
  int overlap = Input_Parameters.NKS_IP.GMRES_Overlap;
  int Ri, Rj;

  if (ii < SolnBlk.ICl -overlap || ii > SolnBlk.ICu + overlap ||
      jj < SolnBlk.JCl -overlap || jj > SolnBlk.JCu + overlap) {

     // GHOST CELL so do nothing
     cout<<"\n Hey I am not suppose to be here! \n"; exit(1);

  } else {     
     int NUM_VAR_LESPREMIXED2D = dRdW.get_n();   
     DenseMatrix dFidW(NUM_VAR_LESPREMIXED2D, NUM_VAR_LESPREMIXED2D,ZERO);
     
     Vector2D nface, DX;   double lface;
     LESPremixed2D_pState Wa, wavespeeds, Left, Right, Wl, Wr;   
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
     
     DenseMatrix A( Rotation_matrix2(nface,NUM_VAR_LESPREMIXED2D, 1));
     DenseMatrix AI( Rotation_matrix2(nface,NUM_VAR_LESPREMIXED2D, 0));
     
     //Left and Right States 
     Inviscid_Flux_Used_Reconstructed_LeftandRight_States(Wl, Wr, DX, SolnBlk, Orient, Ri, Rj );   
     Left  = Rotate(Wl, nface);
     Right = Rotate(Wr, nface);
     
     /**********************************************/
     //Jacobian using Finite Differences
     LESPremixed2D_cState FluxA, FluxB; 
     LESPremixed2D_pState WA, WB;
  
     double perturb = 5e-6;
     double a;
     
     //For Preconditioning
     double delta_n = min(TWO*(SolnBlk.Grid.Cell[ii][jj].A/
			      (SolnBlk.Grid.lfaceE(ii, jj)+SolnBlk.Grid.lfaceW(ii, jj))),
			 TWO*(SolnBlk.Grid.Cell[ii][jj].A/
			      (SolnBlk.Grid.lfaceN(ii, jj)+SolnBlk.Grid.lfaceS(ii, jj))));
 
     for(int jcol=0; jcol<(NUM_VAR_LESPREMIXED2D); jcol++){
       WA = Right;
       WB = Right;
       
       if( jcol <NUM_LESPREMIXED2D_VAR_SANS_SPECIES) {
	 WA[jcol+1] += perturb*max(ONE,Right[jcol+1]); 	 
	 WB[jcol+1] -= perturb*max(ONE,Right[jcol+1]); 
       } else {
	 a =  perturb*max(ONE,Right[jcol+1]); 
	 WA[jcol+1] += a;
	 WA[NUM_VAR_LESPREMIXED2D+1] -= a;      
	 WB[jcol+1] -= a;
	 WB[NUM_VAR_LESPREMIXED2D+1] += a;
       }

       FluxA = FluxRoe_x(Left,WA, Input_Parameters.Preconditioning,SolnBlk.Flow_Type,delta_n);       
       FluxB = FluxRoe_x(Left,WB, Input_Parameters.Preconditioning,SolnBlk.Flow_Type,delta_n);
       
       for(int irow=0; irow<(NUM_VAR_LESPREMIXED2D); irow++){
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
int Inviscid_Flux_Used_Reconstructed_LeftandRight_States(LESPremixed2D_pState &Wl, LESPremixed2D_pState &Wr, 
                                                         Vector2D &DX, LESPremixed2D_Quad_Block &SolnBlk, 
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
void dFIdW_Inviscid_AUSM_plus_up(DenseMatrix& dRdW, LESPremixed2D_Quad_Block &SolnBlk,  
				 LESPremixed2D_Input_Parameters &Input_Parameters,
				 const int &ii, const int &jj, const int Orient){
   
  int overlap = Input_Parameters.NKS_IP.GMRES_Overlap;
  int Ri, Rj;

  if (ii < SolnBlk.ICl -overlap || ii > SolnBlk.ICu + overlap ||
      jj < SolnBlk.JCl -overlap || jj > SolnBlk.JCu + overlap) {
     // GHOST CELL so do nothing
     cout<<"\n Hey I am not suppose to be here! \n"; exit(1);

  } else {     
    int NUM_VAR_LESPREMIXED2D = dRdW.get_n();  //  SolnBlk.NumVar();    
    DenseMatrix dFidW(NUM_VAR_LESPREMIXED2D, NUM_VAR_LESPREMIXED2D,ZERO);
    
    Vector2D nface,DX; double lface;   
    LESPremixed2D_pState Wa, wavespeeds, Left, Right, Wl, Wr;   
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
     
     DenseMatrix A( Rotation_matrix2(nface,NUM_VAR_LESPREMIXED2D, 1));
     DenseMatrix AI( Rotation_matrix2(nface,NUM_VAR_LESPREMIXED2D, 0));
          
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
void dGVdW_Viscous(DenseMatrix &dRdW, LESPremixed2D_Quad_Block &SolnBlk, 
		   LESPremixed2D_Input_Parameters &Input_Parameters,
		   const int &ii, const int &jj){

#ifdef THICKENED_FLAME_ON   
   int NUM_VAR_LESPREMIXED2D = SolnBlk.NumVar()-3; 
#else
   int NUM_VAR_LESPREMIXED2D = SolnBlk.NumVar()-1;
#endif
   int ns = SolnBlk.W[ii][jj].ns-1;

   DenseMatrix dFvdWf(NUM_VAR_LESPREMIXED2D, 14+(ns),ZERO);
   DenseMatrix dWfdWx(14+(ns), NUM_VAR_LESPREMIXED2D,ZERO);  
   DenseMatrix dGvdWf(NUM_VAR_LESPREMIXED2D, 14+(ns),ZERO);
   DenseMatrix dWfdWy(14+(ns), NUM_VAR_LESPREMIXED2D,ZERO);
   
   DenseMatrix dGVdW(NUM_VAR_LESPREMIXED2D,NUM_VAR_LESPREMIXED2D,ZERO);
   Vector2D nface;

   if (ii < SolnBlk.ICl || ii > SolnBlk.ICu ||
       jj < SolnBlk.JCl || jj > SolnBlk.JCu) {
     // GHOST CELL so do nothing to dRdW     
   } else {
     // NON-GHOST CELL.
           
     //Viscous flux Jacobians at North face of cell (ii, jj)
     nface = SolnBlk.Grid.nfaceN(ii, jj);    
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
		    LESPremixed2D_Quad_Block &SolnBlk,
		    const int &Orient, const int &ii, const int &jj){
   
  
   double kappa, Cp, mu, mu_t, kappa_t,Dm_t,Pr_t, Sc_t;
   double sigma, sigma_star, Rmix;
   double rho, U, V, p;  //k, omega;
   double  *h, *dcdx, *dcdy, *dhdT;
   double dUdx,dUdy, dVdx,dVdy;
   //double dkdx, dkdy, domegadx, domegady;
   double drhodx, drhody, dpdx, dpdy, Temp;
   double radius;
   
   LESPremixed2D_pState QuadraturePoint_W;
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
     // dkdx = SolnBlk.dWdx_faceN[ii][jj].k;
//      dkdy = SolnBlk.dWdy_faceN[ii][jj].k;
//      domegadx = SolnBlk.dWdx_faceN[ii][jj].omega;
//      domegady = SolnBlk.dWdy_faceN[ii][jj].omega;

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
     // dkdx = SolnBlk.dWdx_faceE[ii][jj].k;
//      dkdy = SolnBlk.dWdy_faceE[ii][jj].k;
//      domegadx = SolnBlk.dWdx_faceE[ii][jj].omega;
//      domegady = SolnBlk.dWdy_faceE[ii][jj].omega;

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
     // dkdx = SolnBlk.dWdx_faceS[ii][jj].k;
//      dkdy = SolnBlk.dWdy_faceS[ii][jj].k;
//      domegadx = SolnBlk.dWdx_faceS[ii][jj].omega;
//      domegady = SolnBlk.dWdy_faceS[ii][jj].omega;  
     
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
     // dkdx = SolnBlk.dWdx_faceW[ii][jj].k;
//      dkdy = SolnBlk.dWdy_faceW[ii][jj].k;
//      domegadx = SolnBlk.dWdx_faceW[ii][jj].omega;
//      domegady = SolnBlk.dWdy_faceW[ii][jj].omega;
     
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
//    if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
//        SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_EPSILON) {
//       double k = QuadraturePoint_W.k;
//       double omega = QuadraturePoint_W.omega;
//       double taux = (TWO/THREE*dUdx - dVdy/THREE);
//       double tauy = dUdy+dVdx;
//       kappa_t = QuadraturePoint_W.Kappa_turb();
//       Dm_t = QuadraturePoint_W.Dm_turb();
//       mu_t =  QuadraturePoint_W.eddy_viscosity();
//       Pr_t =  QuadraturePoint_W.Pr_turb();
//       Sc_t =  QuadraturePoint_W.Sc_turb();
//       double sigma_s = QuadraturePoint_W.sigma_star;
//       double sigma = QuadraturePoint_W.sigma;
     
//       U= QuadraturePoint_W.v.x;
//       V= QuadraturePoint_W.v.y;
      
//       dFvdWf(1, 0) += TWO*k/max(omega,TOLER)*taux;
//       dFvdWf(1, 4) += TWO*rho/max(omega,TOLER)*taux;
//       dFvdWf(1, 5) += -TWO*rho*k/max(omega*omega,TOLER)*taux;
//       dFvdWf(1, 7) += FOUR/THREE*mu_t;
//       dFvdWf(1, 10) += -TWO/THREE*mu_t;

//       dFvdWf(2, 0) += k/max(omega, TOLER)*tauy;
//       dFvdWf(2, 4) += rho/max(omega, TOLER)*tauy; 
//       dFvdWf(2, 5) += -mu_t/max(omega,TOLER)*tauy;
//       dFvdWf(2, 8) += mu_t;
//       dFvdWf(2, 9) += mu_t;

//       double Sum_dh = 0.0;
//       double Sum_h = 0.0;
    
//        for(int Num = 0; Num<ns_values; Num++){
//          Sum_h +=  h[Num]*dcdx[Num]/Sc_t;
//          Sum_dh -= dhdT[Num]*dcdx[Num]/Sc_t;
//       }
      
//       dFvdWf(3,0) += k/max(omega,TOLER)*Cp/Pr_t*p*drhodx/(Rmix*rho*rho) 
// 	+ Temp/rho*k/max(omega,TOLER)*Sum_dh +k/max(omega,TOLER)*Sum_h 
//          +TWO*QuadraturePoint_W.v.x*k/max(omega,TOLER)*taux+
//          QuadraturePoint_W.v.y*k/max(omega,TOLER)*tauy+k/max(omega,TOLER)*sigma_s*dkdx;
//       dFvdWf(3,1) += TWO*mu_t*taux;
//       dFvdWf(3,2) += mu_t*tauy;
//       dFvdWf(3,3) += -k*Cp*drhodx/(max(omega,TOLER)*Pr_t*Rmix*rho)- k/max(omega,TOLER)*Sum_dh/(Rmix);
//       dFvdWf(3,4) += Cp*(dpdx - p*drhodx/rho)/(max(omega, TOLER) *Pr_t*Rmix) 
// 	+ rho/max(omega,TOLER)*Sum_h+TWO*rho*QuadraturePoint_W.v.x*taux/max(omega,TOLER)
//          +QuadraturePoint_W.v.y*rho*tauy/max(omega,TOLER)+rho*sigma_s*dkdx/max(omega,TOLER);
      
//       dFvdWf(3,5) += -k*Cp*(dpdx-p*drhodx/rho)/(max(omega*omega, TOLER)*Pr_t*Rmix)
// 	- Sum_h*mu_t/max(omega,TOLER)-TWO*QuadraturePoint_W.v.x*mu_t*taux/max(omega,TOLER)
//          -QuadraturePoint_W.v.y*mu_t*tauy/max(omega,TOLER)-mu_t*sigma_s*dkdx/max(omega,TOLER);
//       dFvdWf(3,6) += -k*Cp*Temp/(max(omega,TOLER)*Pr_t);
      
//       dFvdWf(3,7) += FOUR/THREE*QuadraturePoint_W.v.x*mu_t;
//       dFvdWf(3,8) += QuadraturePoint_W.v.y*mu_t;
//       dFvdWf(3,9) += QuadraturePoint_W.v.y*mu_t;
//       dFvdWf(3,10) += -TWO/THREE*QuadraturePoint_W.v.x*mu_t;
//       dFvdWf(3,11) += k*Cp/(max(omega,TOLER)*Rmix*Pr_t);
//       dFvdWf(3,12) += mu+mu_t*sigma_s;
      
//       for(int Num = 0; Num<(ns_species); Num++){ 
//          dFvdWf(3,14+Num) += h[Num]*mu_t/Sc_t;
//          dFvdWf(6+Num,0) += dcdx[Num]*k/max(omega,TOLER)/Sc_t;
//          dFvdWf(6+Num,4) += dcdx[Num]*rho/max(omega,TOLER)/Sc_t;
//          dFvdWf(6+Num,5) += -mu_t*dcdx[Num]/max(omega,TOLER)/Sc_t;     
//          dFvdWf(6+Num,14+Num) += mu_t/Sc_t;
//       } 

//       dFvdWf(4,0) += k*sigma_s*dkdx/max(omega,TOLER);
//       dFvdWf(4,4) += rho*sigma_s*dkdx/max(omega,TOLER);
//       dFvdWf(4,5) -=mu_t*sigma_s*dkdx/max(omega,TOLER); 
//       dFvdWf(4,12) += mu+mu_t*sigma_s;
      
//       dFvdWf(5,0) += k*sigma*domegadx/max(omega,TOLER);
//       dFvdWf(5,4) += rho*sigma*domegadx/max(omega,TOLER);
//       dFvdWf(5,5) += -mu_t*sigma*domegadx/max(omega,TOLER);
//       dFvdWf(5,13) += mu+mu_t*sigma;

//       if(SolnBlk.Axisymmetric == AXISYMMETRIC_Y){            
// 	dFvdWf(1,0) -= TWO/THREE*V*k/(max(omega,TOLER)*radius);
// 	dFvdWf(1,2) -= TWO/THREE*mu_t/radius;
// 	dFvdWf(1,4) -= TWO/THREE*rho*V/(max(omega,TOLER)*radius);
// 	dFvdWf(1,5) += TWO/THREE*mu_t*V/(max(omega,TOLER)*radius);	
// 	dFvdWf(3,0) -= TWO/THREE*U*k*V/(max(omega,TOLER)*radius);
// 	dFvdWf(3,1) -= TWO/THREE*mu_t*V/radius;
// 	dFvdWf(3,2) -= TWO/THREE*U*mu_t/radius;
// 	dFvdWf(3,4) -= TWO/THREE*U*V*rho/(max(omega,TOLER)*radius);
// 	dFvdWf(3,5) += TWO/THREE*U*V*mu_t/(max(omega,TOLER)*radius);                 
//       }//endofaxisymmetric                 

//    }//endof turbulence


   
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


//    if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
//        SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_EPSILON) {
//       double k = QuadraturePoint_W.k;
//       double omega = QuadraturePoint_W.omega;
//       double taux = (TWO/THREE*dVdy - dUdx/THREE);
//       double tauy = dUdy+dVdx;
//       kappa_t = QuadraturePoint_W.Kappa_turb();
//       Dm_t = QuadraturePoint_W.Dm_turb();
//       mu_t =  QuadraturePoint_W.eddy_viscosity();
//       Pr_t =  QuadraturePoint_W.Pr_turb();
//       Sc_t =  QuadraturePoint_W.Sc_turb();
//       double sigma_s = QuadraturePoint_W.sigma_star;
//       double sigma = QuadraturePoint_W.sigma;
      
//       U= QuadraturePoint_W.v.x;
//       V= QuadraturePoint_W.v.y;

//       dGvdWf(1, 0) += k/max(omega, TOLER)*tauy;
//       dGvdWf(1, 4) += rho/max(omega, TOLER)*tauy;
//       dGvdWf(1, 5) += -mu_t/max(omega,TOLER)*tauy;
//       dGvdWf(1, 8) += mu_t;
//       dGvdWf(1, 9) += mu_t;

//       dGvdWf(2, 0) += TWO*k/max(omega,TOLER)*taux;
//       dGvdWf(2, 4) += TWO*rho/max(omega,TOLER)*taux;
//       dGvdWf(2, 5) += -TWO*rho*k/max(omega*omega,TOLER)*taux;
//       dGvdWf(2, 7) += -TWO/THREE*mu_t;
//       dGvdWf(2, 10) += FOUR/THREE*mu_t;
 
//       double Sum_dh = 0.0;
//       double Sum_h = 0.0;
//       for(int Num = 0; Num<ns_values; Num++){
//          Sum_dh  -=  dhdT[Num]*dcdy[Num]/Sc_t;
//          Sum_h +=  h[Num]*dcdy[Num]/Sc_t;
//       }

//       dGvdWf(3,0) += k/max(omega,TOLER)*Cp/Pr_t*drhody*Temp/rho 
// 	+ Temp/rho*k/max(omega,TOLER)* Sum_dh+k/max(omega,TOLER)*Sum_h
// 	+TWO*QuadraturePoint_W.v.y*k/max(omega,TOLER)*taux+
// 	QuadraturePoint_W.v.x*k/max(omega,TOLER)*tauy+k/max(omega,TOLER)*sigma_s*dkdy;
//       dGvdWf(3,1) += mu_t*tauy;
//       dGvdWf(3,2) += TWO*mu_t*taux;
//       dGvdWf(3,3) += -k/max(omega,TOLER)*Cp/Pr_t*drhody/(rho*Rmix) -k/max(omega,TOLER)*ONE/(Rmix)*Sum_dh;
//       dGvdWf(3,4) += Cp*(dpdy - p/rho*drhody)/(max(omega, TOLER)*Pr_t*Rmix) + 
// 	rho*Sum_h/max(omega,TOLER)+TWO*rho*QuadraturePoint_W.v.y*taux/max(omega,TOLER) 
// 	+ QuadraturePoint_W.v.x*rho*tauy/max(omega,TOLER)+rho*sigma_s*dkdy/max(omega,TOLER);      
//       dGvdWf(3,5) += -k*Cp*(dpdy-p/rho*drhody)/(max(omega*omega,TOLER)*Pr_t*Rmix) - 
// 	Sum_h*mu_t/max(omega,TOLER)-TWO*QuadraturePoint_W.v.y*mu_t*taux/max(omega,TOLER)
//          -QuadraturePoint_W.v.x*mu_t*tauy/max(omega,TOLER)-mu_t*sigma_s*dkdy/max(omega,TOLER);
//       dGvdWf(3,6) += -k/max(omega,TOLER)*Cp/Pr_t*Temp;
      
//       dGvdWf(3,7) += -TWO/THREE*QuadraturePoint_W.v.y*mu_t;
//       dGvdWf(3,8) += QuadraturePoint_W.v.x*mu_t;
//       dGvdWf(3,9) +=QuadraturePoint_W.v.x*mu_t;
//       dGvdWf(3,10) += FOUR/THREE*QuadraturePoint_W.v.y*mu_t;
//       dGvdWf(3,11) += k*Cp/(max(omega,TOLER)*Pr_t*Rmix);
//       dGvdWf(3,12) += mu+mu_t*sigma_s;
      
//       for(int Num = 0; Num<(ns_species); Num++){ 
//          dGvdWf(3,14+Num) += h[Num]*mu_t/Sc_t;
//          dGvdWf(6+Num,0) +=k/max(omega,TOLER)*dcdy[Num]/Sc_t;  
//          dGvdWf(6+Num,4) += dcdy[Num]*rho/max(omega,TOLER)/Sc_t;
//          dGvdWf(6+Num,5) += -dcdy[Num]*mu_t/max(omega,TOLER)/Sc_t;     
//          dGvdWf(6+Num,14+Num) += mu_t/Sc_t;
//       } 

//       dGvdWf(4,0) += k*sigma_s*dkdy/max(omega,TOLER);
//       dGvdWf(4,4) += rho*sigma_s*dkdy/max(omega,TOLER);
//       dGvdWf(4,5) -=  mu_t*sigma_s*dkdy/max(omega,TOLER); 
//       dGvdWf(4,12) +=mu+mu_t*sigma_s;
      
//       dGvdWf(5,0) += k*sigma*domegady/max(omega,TOLER);
//       dGvdWf(5,4) += rho*sigma*domegady/max(omega,TOLER);
//       dGvdWf(5,5) += -mu_t*sigma*domegady/max(omega,TOLER);
//       dGvdWf(5,13) += mu+mu_t*sigma;

//       if(SolnBlk.Axisymmetric == AXISYMMETRIC_Y){        
// 	dGvdWf(2,0) -=  TWO/THREE*V*k/(max(omega,TOLER)*radius);
// 	dGvdWf(2,2) -=  TWO/THREE*mu_t/radius;
// 	dGvdWf(2,4) -=  TWO/THREE*V*rho/(max(omega,TOLER)*radius);
// 	dGvdWf(2,5) +=  TWO/THREE*V*mu_t/(max(omega,TOLER)*radius);	
// 	dGvdWf(3,0) -=  TWO/THREE*V*V*k/(max(omega,TOLER)*radius);
// 	dGvdWf(3,2) -=  FOUR/THREE*mu_t*V/radius;
// 	dGvdWf(3,4) -=  TWO/THREE*V*V*rho/(max(omega,TOLER)*radius);
// 	dGvdWf(3,5) +=  TWO/THREE*V*V*mu_t/(max(omega,TOLER)*radius);         
//       }//endofaxisymmetric      
//    }//endof turbulence


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
 * on a diamond path reconstruction.                    *
 *                                                      *
 ********************************************************/
void dWfdWc_Diamond(DenseMatrix &dWfdWc_x,DenseMatrix &dWfdWc_y, LESPremixed2D_Quad_Block &SolnBlk,
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
void d_dWd_dW_Diamond(double &d_dWdx_dW, double &d_dWdy_dW, LESPremixed2D_Quad_Block &SolnBlk, 
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
		     LESPremixed2D_Quad_Block &SolnBlk, 
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
 * This routine calculates the Turbulence (k-omega)     *
 * source Jacobian.                                     *
 *                                                      *
 ********************************************************/
int dS_tdW(DenseMatrix &dStdW,  LESPremixed2D_Quad_Block &SolnBlk,
	   double &d_dWdx_dW_C, double &d_dWdy_dW_C,
	   const int &ii, const int &jj){

  dStdW.zero();

  if ( SolnBlk.W[ii][jj].scalar[0]<0.99 && SolnBlk.W[ii][jj].scalar[0]>0.01 &&  
       SolnBlk.dWdx[ii][jj].scalar[0] !=ZERO && SolnBlk.dWdy[ii][jj].scalar[0] != ZERO ) {

	    Tensor2D strain_rate;
	    strain_rate = SolnBlk.W[ii][jj].Strain_Rate(SolnBlk.dWdx[ii][jj], SolnBlk.dWdy[ii][jj], 
						      SolnBlk.Flow_Type, SolnBlk.Axisymmetric, 
						      SolnBlk.Grid.Cell[ii][jj].Xc);  
  if ( SolnBlk.Flow_Type == FLOWTYPE_LAMINAR_C ){

    double t3,t4,t5,t6,t7,t9,t12,t13,t16,t18;
    double tau_fsd = SolnBlk.W[ii][jj].HeatRelease_Parameter();

      t3 = SolnBlk.W[ii][jj].reactants_den*SolnBlk.W[ii][jj].laminar_speed*(1.0+tau_fsd);
      t4 = SolnBlk.dWdx[ii][jj].scalar[0];//cx(c);
      t5 = t4*t4;
      t6 = SolnBlk.dWdy[ii][jj].scalar[0];//cy(c);
      t7 = t6*t6;
      t9 = sqrt(t5+t7);
      t12 = 1.0+tau_fsd*SolnBlk.W[ii][jj].scalar[0];//c;
      t13 = t12*t12;
      t16 = d_dWdx_dW_C;//diff(cx(c),c);
      t18 = d_dWdy_dW_C;//diff(cy(c),c);

      dStdW(4,4) = t3/t9/t13*(t4*t16+t6*t18)-2.0*t3*t9/t13/t12*tau_fsd;
  }

  if  ( SolnBlk.Flow_Type == FLOWTYPE_LAMINAR_C_ALGEBRAIC || 
        SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_LES_C_ALGEBRAIC ) { 

    double t2,t3,t4,t5,t6,t8;
    double t12,t13,t14,t15,t17;
    double t29,t30,t32,t33,t35,t39,t42;
    double tau_fsd = SolnBlk.W[ii][jj].HeatRelease_Parameter();
    double kappa_fsd = SolnBlk.W[ii][jj].Efficiency_Function_Fsd(SolnBlk.dWdx[ii][jj],SolnBlk.dWdy[ii][jj],SolnBlk.Flow_Type,strain_rate); 
    double k_fsd = SolnBlk.W[ii][jj].SFS_Kinetic_Energy_Fsd(SolnBlk.dWdx[ii][jj],SolnBlk.dWdy[ii][jj],SolnBlk.Flow_Type,strain_rate);

      t2 = 1.0+tau_fsd;
      t3 = SolnBlk.dWdx[ii][jj].scalar[0];//cx(c);
      t4 = t3*t3;
      t5 = SolnBlk.dWdy[ii][jj].scalar[0];//cy(c);
      t6 = t5*t5;
      t8 = sqrt(t4+t6);
      t12 = 1.0+tau_fsd*SolnBlk.W[ii][jj].scalar[0];//c;
      t13 = t12*t12;
      t14 = 1/t13;
      t15 = d_dWdx_dW_C;//diff(cx(c),c);
      t17 = d_dWdy_dW_C;//diff(cy(c),c);
      t29 = sqrt(k_fsd);
      t30 = kappa_fsd*t29;
      t32 = 1/t12;
      t33 = t2*SolnBlk.W[ii][jj].scalar[0];//c;
      t35 = 1.0-t33*t32;
      t39 = 1/SolnBlk.W[ii][jj].laminar_speed/SolnBlk.W[ii][jj].filter_width;
      t42 = t30*t33;
      dStdW(4,4) = SolnBlk.W[ii][jj].reactants_den*SolnBlk.W[ii][jj].laminar_speed*(t2/t8*t14*(t3*t15+t5*t17)-2.0*t2*t8/
t13/t12*tau_fsd+t30*t2*t32*t35*t39-t42*t14*t35*t39*tau_fsd+t42*t32*(-t2*t32+t33
*t14*tau_fsd)*t39);
  }

  if ( SolnBlk.Flow_Type == FLOWTYPE_LAMINAR_C_FSD ){
    
    double t3,t4,t5,t6,t7,t9,t12,t13,t16,t18;
    double tau_fsd = SolnBlk.W[ii][jj].HeatRelease_Parameter();

      t3 = SolnBlk.W[ii][jj].reactants_den*SolnBlk.W[ii][jj].laminar_speed*(1.0+tau_fsd);
      t4 = SolnBlk.dWdx[ii][jj].scalar[0];//cx(c);
      t5 = t4*t4;
      t6 = SolnBlk.dWdy[ii][jj].scalar[0];//cy(c);
      t7 = t6*t6;
      t9 = sqrt(t5+t7);
      t12 = 1.0+tau_fsd*SolnBlk.W[ii][jj].scalar[0];//c;
      t13 = t12*t12;
      t16 = d_dWdx_dW_C;//diff(cx(c),c);
      t18 = d_dWdy_dW_C;//diff(cy(c),c);

      dStdW(4,4) = t3/t9/t13*(t4*t16+t6*t18)-2.0*t3*t9/t13/t12*tau_fsd;
  
  }

  if ( SolnBlk.Flow_Type == FLOWTYPE_LAMINAR_NGT_C_FSD ){
      
    //counter-gradient in C
    double t1,t3,t4,t6,t7;
    double t11,t12,t13,t14,t19;
    double t20,t21,t27,t28;
    double t34,t35,t36,t37,t39;
    double t40,t42,t43,t44,t45,t46,t49;
    double t50,t52,t56,t57,t58,t59;
    double t60,t64,t65;
    double t73,t77,t79;
    double t85,t88;
    double t95,t96;
    double t100,t123,t127,t134,t163,t167;

    double tau_fsd = SolnBlk.W[ii][jj].HeatRelease_Parameter();
    double lam_speed_fsd = SolnBlk.W[ii][jj].laminar_speed;
    double rho_r = SolnBlk.W[ii][jj].reactants_den;

      t1 = rho_r*lam_speed_fsd;
      t3 = tau_fsd*lam_speed_fsd;
      t4 = 1.0-SolnBlk.W[ii][jj].scalar[0];//c;
      t6 = d_dWdx_dW_C;//diff(rhox(rho),rho);
      t7 = d_dWdy_dW_C;//diff(rhoy(rho),rho);
      t11 = 1.0-2.0*SolnBlk.W[ii][jj].scalar[0];//c;
      t12 = SolnBlk.dWdx[ii][jj].scalar[0];//cx(c);
      t13 = SolnBlk.dWdy[ii][jj].scalar[0];//cy(c);
      t14 = t12+t13;
      t19 = SolnBlk.dWdx[ii][jj].rho;//rhox(rho);
      t20 = SolnBlk.dWdy[ii][jj].rho;//rhoy(rho);
      t21 = t19+t20;
      t27 = d_dWdx_dW_C;//diff(cx(c),c);
      t28 = d_dWdy_dW_C;//diff(cy(c),c);
      t34 = t12*t12;
      t35 = t13*t13;
      t36 = t34+t35;
      t37 = 1/t36;
      t39 = 1.0-t34*t37;
      t40 = SolnBlk.dWdx[ii][jj].v.x;//Ux(U);
      t42 = t12*t37;
      t43 = SolnBlk.dWdy[ii][jj].v.x;//Uy(U);
      t44 = SolnBlk.dWdx[ii][jj].v.y;//Vx(V);
      t45 = t43+t44;
      t46 = t13*t45;
      t49 = 1.0-t35*t37;
      t50 = SolnBlk.dWdy[ii][jj].v.y;//Vy(V);
      t52 = t39*t40-t42*t46+t49*t50;
      t56 = lam_speed_fsd*(1.0+tau_fsd*SolnBlk.W[ii][jj].scalar[0]);//c);
      t57 = sqrt(t36);
      t58 = 1/t57;
      t59 = t12*t58;
      t60 = SolnBlk.dWdx[ii][jj].scalar[1];//Fsdx(Fsd);
      t64 = t13*t58;
      t65 = SolnBlk.dWdy[ii][jj].scalar[1];//Fsdy(Fsd);
      t73 = -t34*t58-t35*t58;
      t77 = d_dWdx_dW_C;//diff(Ux(U),U);
      t79 = d_dWdy_dW_C;//diff(Uy(U),U);
      t85 = d_dWdx_dW_C;//diff(Vx(V),V);
      t88 = d_dWdy_dW_C;//diff(Vy(V),V);
      t95 = t36*t36;
      t96 = 1/t95;
      t100 = t12*t27+t13*t28;
      t123 = t60*SolnBlk.W[ii][jj].rho+t19*SolnBlk.W[ii][jj].scalar[1];
      t127 = t65*SolnBlk.W[ii][jj].rho+t20*SolnBlk.W[ii][jj].scalar[1];
      t134 = 1/t57/t36;
      t163 = d_dWdx_dW_C;//diff(Fsdx(Fsd),Fsd);
      t167 = d_dWdy_dW_C;//diff(Fsdy(Fsd),Fsd);
      dStdW(4,0) = t1*SolnBlk.W[ii][jj].scalar[1]-t3*(SolnBlk.W[ii][jj].scalar[0]*t4*(t6+t7)+t11*t14);
      dStdW(4,4) = -t3*(t4*t21-SolnBlk.W[ii][jj].scalar[0]*t21-2.0*SolnBlk.W[ii][jj].rho*t14+SolnBlk.W[ii][jj].rho*t11*(t27+t28));
      dStdW(4,5) = t1*SolnBlk.W[ii][jj].rho;
      dStdW(5,0) = t52*SolnBlk.W[ii][jj].scalar[1]-t56*(-t59*(t60+t6*SolnBlk.W[ii][jj].scalar[1])-t64*(t65+t7*SolnBlk.W[ii][jj].scalar[1]))-t3*SolnBlk.W[ii][jj].scalar[1]*t73
;
      dStdW(5,1) = (t39*t77-t42*t13*t79)*SolnBlk.W[ii][jj].scalar[1]*SolnBlk.W[ii][jj].rho;
      dStdW(5,2) = (-t42*t13*t85+t49*t88)*SolnBlk.W[ii][jj].scalar[1]*SolnBlk.W[ii][jj].rho;
      dStdW(5,4) = ((-2.0*t42*t27+2.0*t34*t96*t100)*t40-t27*t37*t46+2.0*t12*
t96*t46*t100-t42*t28*t45+(-2.0*t13*t37*t28+2.0*t35*t96*t100)*t50)*SolnBlk.W[ii][jj].scalar[1]*SolnBlk.W[ii][jj].rho-t3*(-
t59*t123-t64*t127)-t56*(-t27*t58*t123+t12*t134*t123*t100-t28*t58*t127+t13*t134*
t127*t100)-t3*SolnBlk.W[ii][jj].rho*SolnBlk.W[ii][jj].scalar[1]*(-2.0*t59*t27+t34*t134*t100-2.0*t64*t28+t35*t134*t100);
      dStdW(5,5) = t52*SolnBlk.W[ii][jj].rho-t56*(-t59*(t163*SolnBlk.W[ii][jj].rho+t19)-t64*(t167*SolnBlk.W[ii][jj].rho+t20))-t3*SolnBlk.W[ii][jj].rho
*t73;
  }

  if ( SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_SMAGORINSKY ||
       SolnBlk.Flow_Type == FLOWTYPE_FROZEN_TURBULENT_LES_C_FSD ) {

    double t1,t4,t5,t6,t7,t8,t9;
    double t10,t12,t14,t15,t17,t18,t19;
    double t20,t21,t25,t26,t28;
    double t32,t33,t34,t35,t36,t37;
    double t41,t42,t43,t49;
    double t52,t55,t56,t57;
    double t60,t61,t63,t64,t69;
    double t71,t77;
    double t80,t85,t86,t88,t89;
    double t92,t94,t95,t98;
    double t101,t122,t124,t127,t129;
    double t136,t163,t165,t171,t175;

    double tau_fsd = SolnBlk.W[ii][jj].HeatRelease_Parameter();
    double k_fsd = SolnBlk.W[ii][jj].SFS_Kinetic_Energy_Fsd(SolnBlk.dWdx[ii][jj],SolnBlk.dWdy[ii][jj],SolnBlk.Flow_Type,strain_rate);
    double filter = SolnBlk.W[ii][jj].filter_width;
    double kappa_fsd = SolnBlk.W[ii][jj].Efficiency_Function_Fsd(SolnBlk.dWdx[ii][jj],SolnBlk.dWdy[ii][jj],SolnBlk.Flow_Type,strain_rate); 
    double beta_fsd=1.0;

      t1 = SolnBlk.W[ii][jj].reactants_den*SolnBlk.W[ii][jj].laminar_speed;
      t4 = SolnBlk.dWdx[ii][jj].scalar[0];//cx(c);
      t5 = t4*t4;
      t6 = SolnBlk.dWdy[ii][jj].scalar[0];//cy(c);
      t7 = t6*t6;
      t8 = t5+t7;
      t9 = 1/t8;
      t10 = t5*t9;
      t12 = t7*t9;
      t14 = 2.0/3.0-2.0/3.0*t10+t12/3.0;
      t15 = SolnBlk.dWdx[ii][jj].v.x;//Ux(U);
      t17 = t4*t9;
      t18 = SolnBlk.dWdy[ii][jj].v.x;//Uy(U);
      t19 = SolnBlk.dWdx[ii][jj].v.y;//Vx(V);
      t20 = t18+t19;
      t21 = t6*t20;
      t25 = 2.0/3.0-2.0/3.0*t12+t10/3.0;
      t26 = SolnBlk.dWdy[ii][jj].v.y;//Vy(V);
      t28 = t14*t15-t17*t21+t25*t26;
      t32 = SolnBlk.W[ii][jj].laminar_speed*(1.0+tau_fsd*SolnBlk.W[ii][jj].scalar[0]);//c);
      t33 = sqrt(t8);
      t34 = 1/t33;
      t35 = t4*t34;
      t36 = SolnBlk.dWdx[ii][jj].scalar[1];//Fsdx(Fsd);
      t37 = d_dWdx_dW_C;//diff(rhox(rho),rho);
      t41 = t6*t34;
      t42 = SolnBlk.dWdy[ii][jj].scalar[1];//Fsdy(Fsd);
      t43 = d_dWdy_dW_C;//diff(rhoy(rho),rho);
      t49 = SolnBlk.W[ii][jj].laminar_speed*tau_fsd;
      t52 = -t5*t34-t7*t34;
      t55 = sqrt(k_fsd);
      t56 = kappa_fsd*t55;
      t57 = 1/filter;
      t60 = beta_fsd*SolnBlk.W[ii][jj].laminar_speed;
      t61 = SolnBlk.W[ii][jj].scalar[1]*SolnBlk.W[ii][jj].scalar[1];//Fsd*Fsd;
      t63 = 1.0-SolnBlk.W[ii][jj].scalar[0];//c;
      t64 = 1/t63;
      t69 = d_dWdx_dW_C;//diff(Ux(U),U);
      t71 = d_dWdy_dW_C;//diff(Uy(U),U);
      t77 = d_dWdx_dW_C;//diff(Vx(V),V);
      t80 = d_dWdy_dW_C;//diff(Vy(V),V);
      t85 = d_dWdx_dW_C;//diff(cx(c),c);
      t86 = t17*t85;
      t88 = t8*t8;
      t89 = 1/t88;
      t92 = d_dWdy_dW_C;//diff(cy(c),c);
      t94 = t4*t85+t6*t92;
      t95 = 2.0*t5*t89*t94;
      t98 = t6*t9*t92;
      t101 = 2.0*t7*t89*t94;
      t122 = SolnBlk.dWdx[ii][jj].rho;//rhox(rho);
      t124 = t36*SolnBlk.W[ii][jj].rho+t122*SolnBlk.W[ii][jj].scalar[1];
      t127 = SolnBlk.dWdy[ii][jj].rho;//rhoy(rho);
      t129 = t42*SolnBlk.W[ii][jj].rho+t127*SolnBlk.W[ii][jj].scalar[1];
      t136 = 1/t33/t8;
      t163 = SolnBlk.W[ii][jj].rho*SolnBlk.W[ii][jj].rho;
      t165 = t63*t63;
      t171 = d_dWdx_dW_C;//diff(Fsdx(Fsd),Fsd);
      t175 = d_dWdy_dW_C;//diff(Fsdy(Fsd),Fsd);
      dStdW(4,0) = t1*SolnBlk.W[ii][jj].scalar[1];
      dStdW(4,5) = t1*SolnBlk.W[ii][jj].rho;
      dStdW(5,0) = t28*SolnBlk.W[ii][jj].scalar[1]-t32*(-t35*(t36+t37*SolnBlk.W[ii][jj].scalar[1])-t41*(t42+t43*SolnBlk.W[ii][jj].scalar[1]))-t49*SolnBlk.W[ii][jj].scalar[1]*t52+t56*SolnBlk.W[ii][jj].scalar[1]*t57-2.0*t60*t61*SolnBlk.W[ii][jj].rho*t64;
      dStdW(5,1) = (t14*t69-t17*t6*t71)*SolnBlk.W[ii][jj].scalar[1]*SolnBlk.W[ii][jj].rho;
      dStdW(5,2) = (-t17*t6*t77+t25*t80)*SolnBlk.W[ii][jj].scalar[1]*SolnBlk.W[ii][jj].rho;
      dStdW(5,4) = ((-4.0/3.0*t86+2.0/3.0*t95+2.0/3.0*t98-t101/3.0)*t15-t85*t9
*t21+2.0*t4*t89*t21*t94-t17*t92*t20+(-4.0/3.0*t98+2.0/3.0*t101+2.0/3.0*t86-t95/3.0)*t26)*SolnBlk.W[ii][jj].scalar[1]*SolnBlk.W[ii][jj].rho-t49*(-t35*t124-t41*t129)-t32*(-t85*t34*t124+t4*t136*t124*t94-t92*t34*t129+t6*t136*t129*t94)-t49*SolnBlk.W[ii][jj].rho*SolnBlk.W[ii][jj].scalar[1]*(-2.0*t35*t85+t5*t136*t94-2.0*t41*t92+t7*t136*t94)-t60*t61*t163/t165;
      dStdW(5,5) = t28*SolnBlk.W[ii][jj].rho-t32*(-t35*(t171*SolnBlk.W[ii][jj].rho+t122)-t41*(t175*SolnBlk.W[ii][jj].rho+t127))-t49*SolnBlk.W[ii][jj].rho*t52+t56*SolnBlk.W[ii][jj].rho*t57-2.0*t60*SolnBlk.W[ii][jj].scalar[1]*t163*t64;

  }

  if ( SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_CHARLETTE  ) {

    double t1,t4,t5,t6,t7,t8,t9;
    double t10,t12,t14,t15,t17,t18,t19;
    double t20,t21,t25,t26,t28;
    double t31,t32,t33,t34,t35,t37,t38;
    double t42,t45,t49;
    double t50,t51,t55,t57;
    double t63,t66;
    double t71,t72,t74,t75,t78;
    double t80,t81,t84,t87;
    double t107,t109,t116;
    double t144,t145,t147,t148;
    double t152,t161,t163,t164,t165,t167,t168,t169;
    double t170,t173,t175,t176;
    double t184,t196,t199;

    double tau_fsd = SolnBlk.W[ii][jj].HeatRelease_Parameter();
    double k_fsd = SolnBlk.W[ii][jj].SFS_Kinetic_Energy_Fsd(SolnBlk.dWdx[ii][jj],SolnBlk.dWdy[ii][jj],SolnBlk.Flow_Type,strain_rate);
    double kappa_fsd = SolnBlk.W[ii][jj].Efficiency_Function_Fsd(SolnBlk.dWdx[ii][jj],SolnBlk.dWdy[ii][jj],SolnBlk.Flow_Type,strain_rate); 
    double beta_fsd=1.0;

      t1 = SolnBlk.W[ii][jj].reactants_den*SolnBlk.W[ii][jj].laminar_speed;
      t4 = SolnBlk.dWdx[ii][jj].scalar[0];//cx(c);
      t5 = t4*t4;
      t6 = SolnBlk.dWdy[ii][jj].scalar[0];//cy(c);
      t7 = t6*t6;
      t8 = t5+t7;
      t9 = 1/t8;
      t10 = t5*t9;
      t12 = t7*t9;
      t14 = 2.0/3.0-2.0/3.0*t10+t12/3.0;
      t15 = SolnBlk.dWdx[ii][jj].v.x;//Ux(U);
      t17 = t4*t9;
      t18 = SolnBlk.dWdy[ii][jj].v.x;//Uy(U);
      t19 = SolnBlk.dWdx[ii][jj].v.y;//Vx(V);
      t20 = t18+t19;
      t21 = t6*t20;
      t25 = 2.0/3.0-2.0/3.0*t12+t10/3.0;
      t26 = SolnBlk.dWdx[ii][jj].scalar[0];//Vy(V);
      t28 = t14*t15-t17*t21+t25*t26;
      t31 = 1.0+tau_fsd*SolnBlk.W[ii][jj].scalar[0];//c;
      t32 = sqrt(t8);
      t33 = 1/t32;
      t34 = t4*t33;
      t35 = SolnBlk.dWdx[ii][jj].scalar[1];//Fsdx(Fsd);
      t37 = t6*t33;
      t38 = SolnBlk.dWdy[ii][jj].scalar[1];//Fsdy(Fsd);
      t42 = tau_fsd*SolnBlk.W[ii][jj].scalar[1];//Fsd;
      t45 = -t5*t33-t7*t33;
      t49 = sqrt(k_fsd);
      t50 = kappa_fsd*t49;
      t51 = 1/SolnBlk.W[ii][jj].filter_width;
      t55 = d_dWdx_dW_C;//diff(Ux(U),U);
      t57 = d_dWdy_dW_C;//diff(Uy(U),U);
      t63 = d_dWdx_dW_C;//diff(Vx(V),V);
      t66 = d_dWdy_dW_C;//diff(Vy(V),V);
      t71 = d_dWdx_dW_C;//diff(cx(c),c);
      t72 = t17*t71;
      t74 = t8*t8;
      t75 = 1/t74;
      t78 = d_dWdx_dW_C;//diff(cy(c),c);
      t80 = t4*t71+t6*t78;
      t81 = 2.0*t5*t75*t80;
      t84 = t6*t9*t78;
      t87 = 2.0*t7*t75*t80;
      t107 = SolnBlk.W[ii][jj].rho*t35;
      t109 = SolnBlk.W[ii][jj].rho*t38;
      t116 = 1/t32/t8;
      t144 = beta_fsd*SolnBlk.W[ii][jj].laminar_speed;
      t145 = 1.0+tau_fsd;
      t147 = t31*t31;
      t148 = 1/t147;
      t152 = t145*t32;
      t161 = 1/t145;
      t163 = t161/SolnBlk.W[ii][jj].scalar[0];//c;
      t164 = t145*SolnBlk.W[ii][jj].scalar[0];//c;
      t165 = 1/t31;
      t167 = 1.0-t164*t165;
      t168 = 1/t167;
      t169 = t31*t168;
      t170 = t163*t169;
      t173 = SolnBlk.W[ii][jj].scalar[1]-t152*t148;
      t175 = t144*t173*SolnBlk.W[ii][jj].scalar[1];
      t176 = SolnBlk.W[ii][jj].scalar[1]*SolnBlk.W[ii][jj].scalar[0];
      t184 = t167*t167;
      t196 = d_dWdx_dW_C;//diff(Fsdx(Fsd),Fsd);
      t199 = d_dWdy_dW_C;//diff(Fsdy(Fsd),Fsd);
      dStdW(4,0) = t1*SolnBlk.W[ii][jj].scalar[1];
      dStdW(4,5) = t1*SolnBlk.W[ii][jj].rho;
      dStdW(5,0) = t28*SolnBlk.W[ii][jj].scalar[1]-SolnBlk.W[ii][jj].laminar_speed*(t31*(-t34*t35-t37*t38)+t42*t45)+t50*SolnBlk.W[ii][jj].scalar[1]*t51;
      dStdW(5,1) = (t14*t55-t17*t6*t57)*SolnBlk.W[ii][jj].scalar[1]*SolnBlk.W[ii][jj].rho;
      dStdW(5,2) = (-t17*t6*t63+t25*t66)*SolnBlk.W[ii][jj].scalar[1]*SolnBlk.W[ii][jj].rho;
      dStdW(5,4) = ((-4.0/3.0*t72+2.0/3.0*t81+2.0/3.0*t84-t87/3.0)*t15-t71*t9*t21+2.0*t4*t75*t21*t80-t17*t78*t20+(-4.0/3.0*t84+2.0/3.0*t87+2.0/3.0*t72-t81/3.0)*t26)*SolnBlk.W[ii][jj].scalar[1]*SolnBlk.W[ii][jj].rho-SolnBlk.W[ii][jj].laminar_speed*(tau_fsd*(-t34*t107-t37*t109)+t31*(-t71*t33*t107+t4*t116*t107*t80-t78*t33*t109+t6*t116*t109*t80)+t42*SolnBlk.W[ii][jj].rho*(-2.0*t34*t71+t5*t116*t80-2.0*t37*t78+t7*t116*t80))-t144*(-t145*t33*t148*t80+2.0*t152/t147/t31*tau_fsd)*SolnBlk.W[ii][jj].scalar[1]*t170+t175*t161/t176*t169-t175*t163*tau_fsd*t168+t175*t163*t31/t184*(-t145*t165+t164*t148*tau_fsd);
      dStdW(5,5) = t28*SolnBlk.W[ii][jj].rho-SolnBlk.W[ii][jj].laminar_speed*(t31*(-t34*SolnBlk.W[ii][jj].rho*t196-t37*SolnBlk.W[ii][jj].rho*t199)+tau_fsd*SolnBlk.W[ii][jj].rho*t45)+t50*SolnBlk.W[ii][jj].rho*t51-t144*SolnBlk.W[ii][jj].scalar[1]*t170-t144*t173*t170;
  }

  if ( SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_LES_NGT_C_FSD_SMAGORINSKY  ) {

    //counter-gradient in C
    double t1,t3,t4,t6,t7;
    double t11,t12,t13,t14,t19;
    double t20,t21,t27,t28;
    double t34,t35,t36,t37,t38;
    double t40,t42,t43,t44,t45,t46,t47,t48,t49;
    double t53,t54,t56;
    double t60,t61,t62,t63,t64,t68,t69;
    double t77;
    double t80,t81,t82,t85,t86,t88;
    double t93,t95;
    double t101,t104,t109;
    double t111,t112,t116,t117;
    double t120,t123,t145,t149,t156;
    double t183,t185,t191,t195;

    double tau_fsd = SolnBlk.W[ii][jj].HeatRelease_Parameter();
    double lam_speed_fsd = SolnBlk.W[ii][jj].laminar_speed;
    double rho_r = SolnBlk.W[ii][jj].reactants_den;
    double k_fsd = SolnBlk.W[ii][jj].SFS_Kinetic_Energy_Fsd(SolnBlk.dWdx[ii][jj],SolnBlk.dWdy[ii][jj],SolnBlk.Flow_Type,strain_rate);
    double filter = SolnBlk.W[ii][jj].filter_width;
    double kappa_fsd = SolnBlk.W[ii][jj].Efficiency_Function_Fsd(SolnBlk.dWdx[ii][jj],SolnBlk.dWdy[ii][jj],SolnBlk.Flow_Type,strain_rate); 
    double beta_fsd=1.0;

      t1 = rho_r*lam_speed_fsd;
      t3 = tau_fsd*lam_speed_fsd;
      t4 = 1.0-SolnBlk.W[ii][jj].scalar[0];//c;
      t6 = d_dWdx_dW_C;//diff(rhox(rho),rho);
      t7 = d_dWdy_dW_C;//diff(rhoy(rho),rho);
      t11 = 1.0-2.0*SolnBlk.W[ii][jj].scalar[0];//c;
      t12 = SolnBlk.dWdx[ii][jj].scalar[0];//cx(c);
      t13 = SolnBlk.dWdy[ii][jj].scalar[0];//cy(c);
      t14 = t12+t13;
      t19 = SolnBlk.dWdx[ii][jj].rho;//rhox(rho);
      t20 = SolnBlk.dWdy[ii][jj].rho;//rhoy(rho);
      t21 = t19+t20;
      t27 = d_dWdx_dW_C;//diff(cx(c),c);
      t28 = d_dWdy_dW_C;//diff(cy(c),c);
      t34 = t12*t12;
      t35 = t13*t13;
      t36 = t34+t35;
      t37 = 1/t36;
      t38 = t34*t37;
      t40 = t35*t37;
      t42 = 2.0/3.0-2.0/3.0*t38+t40/3.0;
      t43 = SolnBlk.dWdx[ii][jj].v.x;//Ux(U);
      t45 = t12*t37;
      t46 = SolnBlk.dWdy[ii][jj].v.x;//Uy(U);
      t47 = SolnBlk.dWdx[ii][jj].v.y;//Vx(V);
      t48 = t46+t47;
      t49 = t13*t48;
      t53 = 2.0/3.0-2.0/3.0*t40+t38/3.0;
      t54 = SolnBlk.dWdy[ii][jj].v.y;//Vy(V);
      t56 = t42*t43-t45*t49+t53*t54;
      t60 = lam_speed_fsd*(1.0+tau_fsd*SolnBlk.W[ii][jj].scalar[0]);//c);
      t61 = sqrt(t36);
      t62 = 1/t61;
      t63 = t12*t62;
      t64 = SolnBlk.dWdx[ii][jj].scalar[1];//Fsdx(Fsd);
      t68 = t13*t62;
      t69 = SolnBlk.dWdy[ii][jj].scalar[1];//Fsdy(Fsd);
      t77 = -t34*t62-t35*t62;
      t80 = sqrt(k_fsd);
      t81 = kappa_fsd*t80;
      t82 = 1/filter;
      t85 = beta_fsd*lam_speed_fsd;
      t86 = SolnBlk.W[ii][jj].scalar[1]*SolnBlk.W[ii][jj].scalar[1];
      t88 = 1/t4;
      t93 = d_dWdx_dW_C;//diff(Ux(U),U);
      t95 = d_dWdy_dW_C;//diff(Uy(U),U);
      t101 = d_dWdx_dW_C;//diff(Vx(V),V);
      t104 = d_dWdy_dW_C;//diff(Vy(V),V);
      t109 = t45*t27;
      t111 = t36*t36;
      t112 = 1/t111;
      t116 = t12*t27+t13*t28;
      t117 = 2.0*t34*t112*t116;
      t120 = t13*t37*t28;
      t123 = 2.0*t35*t112*t116;
      t145 = t64*SolnBlk.W[ii][jj].rho+t19*SolnBlk.W[ii][jj].scalar[1];
      t149 = t69*SolnBlk.W[ii][jj].rho+t20*SolnBlk.W[ii][jj].scalar[1];
      t156 = 1/t61/t36;
      t183 = SolnBlk.W[ii][jj].rho*SolnBlk.W[ii][jj].rho;
      t185 = t4*t4;
      t191 = d_dWdx_dW_C;//diff(Fsdx(Fsd),Fsd);
      t195 = d_dWdy_dW_C;//diff(Fsdy(Fsd),Fsd);
      dStdW(4,0) = t1*SolnBlk.W[ii][jj].scalar[1]-t3*(SolnBlk.W[ii][jj].scalar[0]*t4*(t6+t7)+t11*t14);
      dStdW(4,4) = -t3*(t4*t21-SolnBlk.W[ii][jj].scalar[0]*t21-2.0*SolnBlk.W[ii][jj].rho*t14+SolnBlk.W[ii][jj].rho*t11*(t27+t28));
      dStdW(4,5) = t1*SolnBlk.W[ii][jj].rho;
      dStdW(5,0) = t56*SolnBlk.W[ii][jj].scalar[1]-t60*(-t63*(t64+t6*SolnBlk.W[ii][jj].scalar[1])-t68*(t69+t7*SolnBlk.W[ii][jj].scalar[1]))-t3*SolnBlk.W[ii][jj].scalar[1]*t77
+t81*SolnBlk.W[ii][jj].scalar[1]*t82-2.0*t85*t86*SolnBlk.W[ii][jj].rho*t88;
      dStdW(5,1) = (t42*t93-t45*t13*t95)*SolnBlk.W[ii][jj].scalar[1]*SolnBlk.W[ii][jj].rho;
      dStdW(5,2) = (-t45*t13*t101+t53*t104)*SolnBlk.W[ii][jj].scalar[1]*SolnBlk.W[ii][jj].rho;
      dStdW(5,4) = ((-4.0/3.0*t109+2.0/3.0*t117+2.0/3.0*t120-t123/3.0)*t43-t27
*t37*t49+2.0*t12*t112*t49*t116-t45*t28*t48+(-4.0/3.0*t120+2.0/3.0*t123+2.0/3.0*
t109-t117/3.0)*t54)*SolnBlk.W[ii][jj].scalar[1]*SolnBlk.W[ii][jj].rho-t3*(-t63*t145-t68*t149)-t60*(-t27*t62*t145+t12*t156
*t145*t116-t28*t62*t149+t13*t156*t149*t116)-t3*SolnBlk.W[ii][jj].rho*SolnBlk.W[ii][jj].scalar[1]*(-2.0*t63*t27+t34*t156*
t116-2.0*t68*t28+t35*t156*t116)-t85*t86*t183/t185;
      dStdW(5,5) = t56*SolnBlk.W[ii][jj].rho-t60*(-t63*(t191*SolnBlk.W[ii][jj].rho+t19)-t68*(t195*SolnBlk.W[ii][jj].rho+t20))-t3*SolnBlk.W[ii][jj].rho
*t77+t81*SolnBlk.W[ii][jj].rho*t82-2.0*t85*SolnBlk.W[ii][jj].scalar[1]*t183*t88;
  }

 if ( SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_K ) {

    //Simplified M with rho
 
    double t1,t4,t5,t6,t7,t8,t9;
    double t10,t12,t14,t15,t17,t18,t19;
    double t20,t21,t25,t26,t28;
    double t32,t33,t34,t35,t36,t37;
    double t41,t42,t43,t49;
    double t52,t55,t56,t57;
    double t62,t63,t65,t66;
    double t71,t73,t79;
    double t82,t87,t88;
    double t90,t91,t94,t96,t97;
    double t100,t103,t124,t126,t129,t131,t138;
    double t151,t170,t171,t174,t180,t184;
    double t206,t211,t214,t216,t219,t225,t226,t236,t255,t258,t264;

    double tau_fsd = SolnBlk.W[ii][jj].HeatRelease_Parameter();
    //    double k_fsd = SolnBlk.W[ii][jj].SFS_Kinetic_Energy_Fsd(SolnBlk.dWdx[ii][jj],SolnBlk.dWdy[ii][jj],SolnBlk.Flow_Type,strain_rate);
    double filter = SolnBlk.W[ii][jj].filter_width;
    double kappa_fsd = SolnBlk.W[ii][jj].Efficiency_Function_Fsd(SolnBlk.dWdx[ii][jj],SolnBlk.dWdy[ii][jj],SolnBlk.Flow_Type,strain_rate); 
    double beta_fsd=1.0;
    double Cv=0.088;
    double Cs=0.931;

      t1 = SolnBlk.W[ii][jj].reactants_den*SolnBlk.W[ii][ii].laminar_speed;
      t4 = SolnBlk.dWdx[ii][jj].scalar[0];//cx(c);
      t5 = t4*t4;
      t6 = SolnBlk.dWdy[ii][jj].scalar[0];//cy(c);
      t7 = t6*t6;
      t8 = t5+t7;
      t9 = 1/t8;
      t10 = t5*t9;
      t12 = t7*t9;
      t14 = 2.0/3.0-2.0/3.0*t10+t12/3.0;
      t15 = SolnBlk.dWdx[ii][jj].v.x;//Ux(U);
      t17 = t4*t9;
      t18 = SolnBlk.dWdy[ii][jj].v.x;//Uy(U);
      t19 = SolnBlk.dWdx[ii][jj].v.y;//Vx(V);
      t20 = t18+t19;
      t21 = t6*t20;
      t25 = 2.0/3.0-2.0/3.0*t12+t10/3.0;
      t26 = SolnBlk.dWdy[ii][jj].v.y;//Vy(V);
      t28 = t14*t15-t17*t21+t25*t26;
      t32 = SolnBlk.W[ii][jj].laminar_speed*(1.0+tau_fsd*SolnBlk.W[ii][jj].scalar[0]);//c);
      t33 = sqrt(t8);
      t34 = 1/t33;
      t35 = t4*t34;
      t36 = SolnBlk.dWdx[ii][jj].scalar[1];//Fsdx(Fsd);
      t37 = d_dWdx_dW_C;//diff(rhox(rho),rho);
      t41 = t6*t34;
      t42 = SolnBlk.dWdy[ii][jj].scalar[1];//Fsdy(Fsd);
      t43 = d_dWdy_dW_C;//diff(rhoy(rho),rho);
      t49 = SolnBlk.W[ii][jj].laminar_speed*tau_fsd;
      t52 = -t5*t34-t7*t34;
      t55 = sqrt(SolnBlk.W[ii][jj].scalar[2]);//k);
      t56 = kappa_fsd*t55;
      t57 = 1/filter;
      t62 = (1.0-t10-t12)*beta_fsd*SolnBlk.W[ii][jj].laminar_speed;
      t63 = SolnBlk.W[ii][jj].scalar[1]*SolnBlk.W[ii][jj].scalar[1];//Fsd*Fsd;
      t65 = 1.0-SolnBlk.W[ii][jj].scalar[0];//c;
      t66 = 1/t65;
      t71 = d_dWdx_dW_C;//diff(Ux(U),U);
      t73 = d_dWdy_dW_C;//diff(Uy(U),U);
      t79 = d_dWdx_dW_C;//diff(Vx(V),V);
      t82 = d_dWdy_dW_C;//diff(Vy(V),V);
      t87 = d_dWdx_dW_C;//diff(cx(c),c);
      t88 = t17*t87;
      t90 = t8*t8;
      t91 = 1/t90;
      t94 = d_dWdy_dW_C;//diff(cy(c),c);
      t96 = t4*t87+t6*t94;
      t97 = 2.0*t5*t91*t96;
      t100 = t6*t9*t94;
      t103 = 2.0*t7*t91*t96;
      t124 = SolnBlk.dWdx[ii][jj].rho;//rhox(rho);
      t126 = t36*SolnBlk.W[ii][jj].rho+t124*SolnBlk.W[ii][jj].scalar[1];//Fsd;
      t129 = SolnBlk.dWdy[ii][jj].rho;//rhoy(rho);
      t131 = t42*SolnBlk.W[ii][jj].rho+t129*SolnBlk.W[ii][jj].scalar[1];//Fsd;
      t138 = 1/t33/t8;
      t151 = SolnBlk.W[ii][jj].scalar[1]*SolnBlk.W[ii][jj].rho;//Fsd*rho;
      t170 = SolnBlk.W[ii][jj].rho*SolnBlk.W[ii][jj].rho;//rho*rho;
      t171 = t63*t170;
      t174 = t65*t65;
      t180 = d_dWdx_dW_C;//diff(Fsdx(Fsd),Fsd);
      t184 = d_dWdy_dW_C;//diff(Fsdy(Fsd),Fsd);
      t206 = 0.6666666667*t15-0.3333333333*t26;
      t211 = 0.6666666667*t26-0.3333333333*t15;
      t214 = pow(SolnBlk.W[ii][jj].scalar[2],1.5);
      t216 = 1/filter;
      t219 = 2.0*Cv*filter*SolnBlk.W[ii][jj].scalar[2];//eddy_viscosity*k;
      t225 = 0.6666666667*SolnBlk.W[ii][jj].rho*SolnBlk.W[ii][jj].scalar[2];//rho*k;
      t226 = t219*t206-t225;
      t236 = t219*t211-t225;
      t255 = 0.6666666667*SolnBlk.W[ii][jj].rho;
      t258 = t20*t20;
      t264 = pow(SolnBlk.W[ii][jj].scalar[2],0.5);
      dStdW(4,0) = t1*SolnBlk.W[ii][jj].scalar[1];
      dStdW(4,5) = t1*SolnBlk.W[ii][jj].rho;
      dStdW(5,0) = t28*SolnBlk.W[ii][jj].scalar[1]-t32*(-t35*(t36+t37*SolnBlk.W[ii][jj].scalar[1])-t41*(t42+t43*SolnBlk.W[ii][jj].scalar[1]))-t49*SolnBlk.W[ii][jj].scalar[1]*t52+t56*SolnBlk.W[ii][jj].scalar[1]*t57-2.0*t62*t63*SolnBlk.W[ii][jj].rho*t66;
      dStdW(5,1) = (t14*t71-t17*t6*t73)*SolnBlk.W[ii][jj].scalar[1]*SolnBlk.W[ii][jj].rho;
      dStdW(5,2) = (-t17*t6*t79+t25*t82)*SolnBlk.W[ii][jj].scalar[1]*SolnBlk.W[ii][jj].rho;
      dStdW(5,4) = ((-4.0/3.0*t88+2.0/3.0*t97+2.0/3.0*t100-t103/3.0)*t15-t87*t9*t21+2.0*t4*t91*t21*t96-t17*t94*t20+(-4.0/3.0*t100+2.0/3.0*t103+2.0/3.0*t88                    -t97/3.0)*t26)*SolnBlk.W[ii][jj].scalar[1]*SolnBlk.W[ii][jj].rho-t49*(-t35*t126-t41*t131)-t32*(-t87*t34*t126+t4*t138*t126*t96-t94*t34*t131+t6*t138*t131*t96)-t49*t151*(-2.0*t35*
                    t87+t5*t138*t96-2.0*t41*t94+t7*t138*t96)-(-2.0*t88+t97-2.0*t100+t103)*beta_fsd*SolnBlk.W[ii][jj].laminar_speed*t171*t66-t62*t171/t174;
      dStdW(5,5) = t28*SolnBlk.W[ii][jj].rho-t32*(-t35*(t180*SolnBlk.W[ii][jj].rho+t124)-t41*(t184*SolnBlk.W[ii][jj].rho+t129))-t49*SolnBlk.W[ii][jj].rho*t52+t56*SolnBlk.W[ii][jj].rho*t57-2.0*t62*SolnBlk.W[ii][jj].scalar[1]*t170*t66;
      dStdW(5,6) = kappa_fsd/t55*t151*t57/2.0;
      dStdW(6,0) = 0.6666666667*SolnBlk.W[ii][jj].scalar[2]*t206+0.6666666667*SolnBlk.W[ii][jj].scalar[2]*t211-Cs*t214*t216;
      dStdW(6,1) = -0.6666666667*t219*t71*t206-0.6666666667*t226*t71-2.0*t219*t20*t73+0.3333333333*t219*t71*t211+0.3333333333*t236*t71;
      dStdW(6,2) = 0.3333333333*t219*t82*t206+0.3333333333*t226*t82-2.0*t219*t20*t79-0.6666666667*t219*t82*t211-0.6666666667*t236*t82;
      dStdW(6,6) = -(2.0*Cv*SolnBlk.W[ii][jj].scalar[2]*t206-t255)*t206-2.0*Cv*SolnBlk.W[ii][jj].scalar[2]*t258-(2.0*Cv*SolnBlk.W[ii][jj].scalar[2]*t211-t255)*t211-0.15E1*Cs*SolnBlk.W[ii][jj].rho*t264*t216;
  }
  }  
  
//    double rho, U,V, k, omega;
//    double dUdx, dUdy, dVdx, dVdy;
//    double mu_t, alpha, beta_star, beta;
   
//    beta_star = SolnBlk.W[ii][jj].beta_star;
//    beta = SolnBlk.W[ii][jj].beta;
//    alpha = SolnBlk.W[ii][jj].alpha;
//    mu_t = SolnBlk.W[ii][jj].eddy_viscosity();
   
//    rho = SolnBlk.W[ii][jj].rho;
//    U = SolnBlk.W[ii][jj].v.x;
//    V = SolnBlk.W[ii][jj].v.y;
//    k = SolnBlk.W[ii][jj].k;
//    omega = SolnBlk.W[ii][jj].omega;
   
//    dUdx =  SolnBlk.dWdx[ii][jj].v.x;
//    dUdy =  SolnBlk.dWdy[ii][jj].v.x;
//    dVdx =  SolnBlk.dWdx[ii][jj].v.y;
//    dVdy =  SolnBlk.dWdy[ii][jj].v.y;

// //    d_dWdx_dW =  SolnBlk.d_dWdx_dW[ii][jj][CENTER];
// //    d_dWdy_dW =  SolnBlk.d_dWdy_dW[ii][jj][CENTER];

//    double tau_x = TWO/THREE*dUdx - dVdy/THREE;
//    double tau_y = TWO/THREE*dVdy-dUdx/THREE;
//    double tau_z = dUdy+dVdx;
   
//    dStdW(4,0) += (TWO*k/max(omega,TOLER)*tau_x - TWO/THREE*k)*dUdx +k/max(omega,TOLER)
//       *tau_z*tau_z + (TWO*k/max(omega,TOLER)*tau_y - TWO/THREE*k)*dVdy - beta_star*k*omega;
   
//    dStdW(4,1) += FOUR/THREE*mu_t*d_dWdx_dW*dUdx+(TWO*mu_t*tau_x - TWO/THREE*rho*k)*d_dWdx_dW
//       + TWO*mu_t*tau_z*d_dWdy_dW - TWO/THREE*mu_t*d_dWdx_dW*dVdy;
   
//    dStdW(4,2) += -TWO/THREE*mu_t*d_dWdy_dW*dUdx+TWO*mu_t*tau_z*d_dWdx_dW+FOUR/THREE*
//       mu_t*d_dWdy_dW*dVdy+(TWO*mu_t*tau_y-TWO/THREE*rho*k)*d_dWdy_dW;
   
//    dStdW(4,4) += (TWO*rho/max(omega,TOLER)*tau_x - TWO/THREE*rho)*dUdx +rho/max(omega,TOLER)*
//       tau_z*tau_z + (TWO*rho/max(omega,TOLER)*tau_y - TWO/THREE*rho)*dVdy - beta_star*rho*omega;

//    dStdW(4,5) += -TWO*mu_t/max(omega,TOLER)*tau_x*dUdx - mu_t/max(omega,TOLER)*tau_z*tau_z -
//       TWO*mu_t/max(omega,TOLER)*tau_y*dVdy - beta_star*rho*k;
  
//    dStdW(5,0) +=alpha*omega/max(k,TOLER)*((TWO*k/max(omega,TOLER)*tau_x - TWO/THREE*k)*dUdx +k/max(omega,TOLER)*tau_z*tau_z
//                                          + (TWO*k/max(omega,TOLER)*tau_y - TWO/THREE*k)*dVdy) - beta*omega*omega;
   
//    dStdW(5,1) +=alpha*omega/max(k,TOLER)*(FOUR/THREE*mu_t*d_dWdx_dW*dUdx+( TWO*mu_t*tau_x - TWO/THREE*rho*k)*d_dWdx_dW +
//                                          TWO*mu_t*tau_z*d_dWdy_dW - TWO/THREE*mu_t*d_dWdx_dW*dVdy);
   
//    dStdW(5,2) +=alpha*omega/max(k,TOLER)*(-TWO/THREE*mu_t*d_dWdy_dW*dUdx + TWO*mu_t*tau_z*d_dWdx_dW + 
//       FOUR/THREE*mu_t*d_dWdy_dW*dVdy + (TWO*mu_t*tau_y - TWO/THREE*rho*k)*d_dWdy_dW  );
   
//    dStdW(5,4) += -alpha*omega/max(k*k,TOLER)*((TWO*mu_t*tau_x - TWO/THREE*rho*k )*dUdx + mu_t*tau_z*tau_z 
//                                              +(TWO*mu_t*tau_y - TWO/THREE*rho*k)*dVdy) +
//       alpha*omega/max(k,TOLER)*((TWO*rho/max(omega,TOLER)*tau_x - TWO/THREE*rho)*dUdx+rho/max(omega,TOLER)*tau_z*tau_z 
//        + (TWO*rho/max(omega,TOLER)*tau_y-TWO/THREE*rho)*dVdy);
   
//    dStdW(5,5) +=  alpha/max(k,TOLER)*((TWO*mu_t*tau_x - TWO/THREE*rho*k)*dUdx + mu_t*tau_z*tau_z+(TWO*mu_t*tau_y-TWO*rho*k)*dVdy) +
//       alpha/max(k,TOLER)*(-TWO*mu_t*tau_x*dUdx - mu_t*tau_z*tau_z - TWO*mu_t*tau_y*dVdy) - TWO*beta*rho*omega;


//    if(SolnBlk.Axisymmetric ==AXISYMMETRIC_Y){
//       double radius  = SolnBlk.Grid.Cell[ii][jj].Xc.y;
//       if(radius !=ZERO){
//          double tau_t = TWO/THREE*V/radius - dUdx/THREE - dVdy/THREE;
        
//          dStdW(4,0) += -TWO/THREE*V*k/max(omega,TOLER)/radius*(dUdx +dVdy) + (TWO*k*tau_t/max(omega,TOLER) - TWO/THREE*k)*V/radius;
//          dStdW(4,1) -=FOUR/THREE*mu_t*V*d_dWdx_dW/radius;
//          dStdW(4,2) += - TWO/THREE*mu_t*(dUdx+dVdy)/radius - TWO/THREE*mu_t*V*d_dWdy_dW/radius+TWO*mu_t*(TWO/THREE/radius - d_dWdy_dW/THREE)*V/radius
// 	   + (TWO*mu_t*tau_t - TWO/THREE*rho*k)/radius;
//          dStdW(4,4) += -TWO/THREE*rho*V/max(omega,TOLER)/radius*(dUdx+dVdy) + V*(TWO*rho/max(omega,TOLER)*tau_t - TWO/THREE*rho)/radius;
//          dStdW(4,5) += TWO/THREE*mu_t*V/max(omega,TOLER)/radius*(dUdx+dVdy) - TWO*mu_t*V*tau_t/max(omega,TOLER)/radius;         
//          dStdW(5,0) += alpha*omega/max(TOLER,k)*(-TWO/THREE*k*V*(dUdx+dVdy)/max(omega,TOLER)/radius + (TWO*k*tau_t/max(omega,TOLER) - TWO/THREE*k)*V/radius);
//          dStdW(5,1) -= FOUR/THREE*alpha*rho*V*d_dWdx_dW/radius;
//          dStdW(5,2) += alpha*omega/max(k,TOLER)*(-TWO/THREE*mu_t*(dUdx+dVdy)/radius-TWO/THREE*mu_t*V*d_dWdy_dW/radius
//                                                  +TWO*mu_t*(TWO/THREE/radius - d_dWdy_dW/THREE)*V/radius + (TWO*mu_t*tau_t - TWO/THREE*rho*k)/radius);
// 	 dStdW(5,4) += - alpha*omega/max(k*k,TOLER)*(-TWO/THREE*mu_t*V*(dUdx + dVdy)/radius + (TWO*mu_t*tau_t - TWO/THREE*rho*k)*V/radius)+
// 	   alpha*omega/max(TOLER,k)*(-TWO/THREE*V*rho*(dUdx+dVdy)/max(omega,TOLER)/radius +(TWO*rho*tau_t/max(omega,TOLER) - TWO/THREE*rho)*V/radius );         
//          dStdW(5,5) += alpha/max(k,TOLER)*(-TWO/THREE*mu_t*V*(dUdx+dVdy)/radius + (TWO*mu_t*tau_t - TWO/THREE*rho*k)*V/radius )+
//             alpha/max(k,TOLER)*(TWO/THREE*mu_t*V*(dUdx+dVdy)/radius - TWO*mu_t*V*tau_t/radius);
//       }      
//    }
   
    return (0);
 }





// ARE THES FUNCTIONS EVER CALLED ????


// int Automatic_Wall_Treatment_Residual_Jacobian(LESPremixed2D_Quad_Block &SolnBlk, LESPremixed2D_Input_Parameters &Input_Parameters, int i, int j, DenseMatrix &dRdU){
   
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

// int Wall_Function_Residual_Jacobian(LESPremixed2D_Quad_Block &SolnBlk, 
//                          LESPremixed2D_Input_Parameters &Input_Parameters, 
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


// void Low_Reynoldsnumber_Formulation_Residual_Jacobian(LESPremixed2D_Quad_Block &SolnBlk, 
// 						      LESPremixed2D_Input_Parameters &Input_Parameters, 
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


// void BC_Residual_Jacobian(LESPremixed2D_Quad_Block &SolnBlk, LESPremixed2D_Input_Parameters &Input_Parameters, int i, int j, DenseMatrix &dRdU){
   
//   if (((i==SolnBlk.ICl) && ( (SolnBlk.Grid.BCtypeW[j] != BC_NONE ))) 
//       ||  ((i==SolnBlk.ICu) &&(SolnBlk.Grid.BCtypeE[j] != BC_NONE))
//       || ((j == SolnBlk.JCl) &&((SolnBlk.Grid.BCtypeS[i] != BC_NONE )))
//       || ((j ==SolnBlk.JCu) &&((SolnBlk.Grid.BCtypeN[i] != BC_NONE )))){
          
//     dRdU.zero();
//     int NUM_VAR_LESPREMIXED2D = SolnBlk.NumVar();
    
//     for(int irow=0; irow<(NUM_VAR_LESPREMIXED2D-1); irow++)
//       for(int jcol=0; jcol<(NUM_VAR_LESPREMIXED2D-1); jcol++){
// 	if(irow==jcol){
// 	  dRdU(irow, jcol) = -ONE/SolnBlk.dt[i][j];
// 	}
//       }    
//   }
// }



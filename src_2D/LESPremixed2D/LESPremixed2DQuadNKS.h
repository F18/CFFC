#ifndef _LESPREMIXED2D_NKS_INCLUDED 
#define _LESPREMIXED2D_NKS_INCLUDED 

/****************************************************************/
/****** NKS REQUIRED SPECIALIZATIONS & FUNCTIONS ****************/
/****************************************************************/  

#include "LESPremixed2DQuad.h"
#include "LESPremixed2DdRdU.h"
#include "../NewtonKrylovSchwarz2D/NKS2D.h"


/*! *****************************************************************************************
 *  LESPremixed2D Specialization of blocksize to use N-1 not N variables                    *
 ********************************************************************************************/
template <> int set_blocksize(LESPremixed2D_Quad_Block &SolnBlk){ 
#ifdef THICKENED_FLAME_ON
  return (SolnBlk.NumVar()-3);
#else
  return (SolnBlk.NumVar()-1);
#endif 
}


/*! *****************************************************************************************
 *  Specialization of Newton_Update                                                         *
 *                                                                                          *
 * This routine updates the previous Solution Data Uo with the deltaU from the GMRES        *                  
 * iterative solver.  U = Uo + GMRES.delatU                                                 *
 *                                                                                          *
 ********************************************************************************************/
template <>
int Newton_Update(LESPremixed2D_Quad_Block *SolnBlk,
		  AdaptiveBlock2D_List &List_of_Local_Solution_Blocks,
		  LESPremixed2D_Input_Parameters &Input_Parameters,
		  GMRES_RightPrecon_MatrixFree<LESPremixed2D_pState,LESPremixed2D_Quad_Block,LESPremixed2D_Input_Parameters> &GMRES,
      double Relaxation_multiplier) {

#ifdef THICKENED_FLAME_ON
  int Num_Var = SolnBlk[0].NumVar()-2;
#else
  int Num_Var = SolnBlk[0].NumVar();
#endif  	
  int error_flag = 0;
  
   /* Update Solution. No updates to Ghost Cells, let the BC's take care of it */
  for ( int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks.Nblk ; ++Bcount ) {
    if (List_of_Local_Solution_Blocks.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
      for (int j = SolnBlk[Bcount].JCl; j <= SolnBlk[Bcount].JCu; j++){
	for (int i = SolnBlk[Bcount].ICl; i <= SolnBlk[Bcount].ICu; i++){
	  
 	  /* Update solutions in conserved variables  U = Uo + RELAXATION*deltaU = Uo + denormalized(x) */	 
 	  for(int varindex =1; varindex < Num_Var; varindex++){	                       
 	    SolnBlk[Bcount].U[i][j][varindex] = SolnBlk[Bcount].Uo[i][j][varindex]  +  
	      Relaxation_multiplier*GMRES.deltaU(Bcount,i,j,varindex-1);
 	  } 	      	  
 	  //LESPREMIXED2D N-1 
 	  SolnBlk[Bcount].U[i][j][Num_Var] = SolnBlk[Bcount].U[i][j].rho*(ONE - SolnBlk[Bcount].U[i][j].sum_species());	   
	  	  
  	  /**************************************************************************/
 	  // Apply update reduction while any one of the updated variables is unphysical 
 	  if(! SolnBlk[Bcount].U[i][j].Unphysical_Properties_Check( SolnBlk[Bcount].Flow_Type, 10)){	   
 	    double update_reduction_factor = ONE;	    
 	    for (int n_update_reduction = 1; n_update_reduction <= 10; ++n_update_reduction) {		  
 	      update_reduction_factor = HALF*update_reduction_factor;		  		  
 	      for(int varindex = 1; varindex <= Num_Var-1; varindex++){		              
 		SolnBlk[Bcount].U[i][j][varindex] = SolnBlk[Bcount].Uo[i][j][varindex] 
 		  + GMRES.deltaU(Bcount,i,j,varindex-1)*update_reduction_factor;
 	      }   
 	      SolnBlk[Bcount].U[i][j][Num_Var] = SolnBlk[Bcount].U[i][j].rho*(ONE - SolnBlk[Bcount].U[i][j].sum_species());
 	      cout<<"\n Applying Reduction to solution in NKS "<<n_update_reduction;
 	      if( SolnBlk[Bcount].U[i][j].Unphysical_Properties_Check( SolnBlk[Bcount].Flow_Type,n_update_reduction))  break;	      
 	    } 
 	  }
	 
 	  /**************************************************************************/
 	  // Error Check
 	  if(! SolnBlk[Bcount].U[i][j].Unphysical_Properties_Check(SolnBlk[Bcount].Flow_Type,10)) error_flag = 1;
	 	  
 	  /**************************************************************************/
 	  //Update solution in primitive variables.	 
 	  SolnBlk[Bcount].W[i][j] = W(SolnBlk[Bcount].U[i][j]); 
	    
	} 
      }
      
// #ifdef _NKS_VERBOSE  
//       if (CFFC_Primary_MPI_Processor()) { 
// 	/**************************************************************************/
// 	//FOR DEBUGGING OUTPUT GMRES deltaU set to "SCALED" residual
// 	double *norm = new double[Num_Var-1];
// 	for(int i= 0; i< Num_Var-1; i++) norm[i]=ZERO;
	
	for (int j = SolnBlk[Bcount].JCl-SolnBlk[Bcount].Nghost; j <= SolnBlk[Bcount].JCu+SolnBlk[Bcount].Nghost; j++){
	  for (int i = SolnBlk[Bcount].ICl-SolnBlk[Bcount].Nghost; i <= SolnBlk[Bcount].ICu+SolnBlk[Bcount].Nghost; i++){
	    for(int varindex =1; varindex < Num_Var; varindex++){	  	   
	      //SolnBlk[Bcount].dUdt[i][j][0][varindex] = GMRES.deltaU_test(Bcount,i,j,varindex-1);
	      SolnBlk[Bcount].dUdt[i][j][0][varindex] = GMRES.b_test(Bcount,i,j,varindex-1);
	      //norm[varindex-1] += sqr(SolnBlk[Bcount].dUdt[i][j][0][varindex]);
	      // norm[varindex-1] = max(norm[varindex-1],fabs(SolnBlk[Bcount].dUdt[i][j][0][varindex]));
	    }
	  } 
	}
// 	cout<<"\n *************** ";
// 	for(int i= 0; i<11; i++){
// 	  //cout<<"\n L2 norm of variable "<<i<<" = "<<sqrt(norm[i]);
// 	  cout<<"\n max norm of variable "<<i<<" = "<<norm[i];
// 	}
// 	cout<<"\n *************** ";
// 	delete[] norm;
	/**************************************************************************/
//       }
// #endif

    } 
  }   
  return error_flag; 
}

/*! *****************************************************************************************
 * LESPremixed2D Specific  Finite_Time_Step                                                 *
 *                                                                                          *
 * This routine calculates the Finite Time Step using a basic SER approach for startup,     *                  
 * returning a "CFL" number to multiply by the stability determined dt.                     *
 *                                                                                          *
 ********************************************************************************************/
template <> double Finite_Time_Step(const LESPremixed2D_Input_Parameters &Input_Parameters, 
				    const double &L2norm_first,
				    const double &L2norm_current,
				    const double &L2norm_current_n,		
				    const int &Number_of_Newton_Steps) {

  double CFL_current;

  //SER 
  if (L2norm_current_n > MIN_FINITE_TIME_STEP_NORM_RATIO ) { 
    //Original works for Invisicd, Viscous OK, Reacting seems to be more picky
    CFL_current = Input_Parameters.NKS_IP.Finite_Time_Step_Initial_CFL*
      pow( max(ONE, ONE/L2norm_current_n),ONE ); 
    
    //Possible to help curb jumps, !!!!!!!!!!!!!!!! however can cause oscillations in residual!!!!!!!!!!!!!!!!
    if( L2norm_first == L2norm_current && Number_of_Newton_Steps != 1 &&
	CFL_current <= Input_Parameters.NKS_IP.Finite_Time_Step_Initial_CFL){ //ie still ONE
      //CFL_current = 0.1*CFL_current;      
    }

  } else {
    CFL_current = Input_Parameters.NKS_IP.Finite_Time_Step_Initial_CFL/MIN_FINITE_TIME_STEP_NORM_RATIO;
  }  

  // CAN SET A MAXIMUM, SOMETIMES USEFUL WITH CONVERGENCE STALL & BC ISSUES
  if( CFL_current > Input_Parameters.NKS_IP.Finite_Time_Step_Max_CFL) 
    CFL_current = Input_Parameters.NKS_IP.Finite_Time_Step_Max_CFL;

  return CFL_current;

}


/*****************************************************************************
 * SPECIALIZATION of finite time step addition to diagonal.                  *
 *****************************************************************************/ 
template<> inline void Block_Preconditioner<LESPremixed2D_pState,
					    LESPremixed2D_Quad_Block,					    
					    LESPremixed2D_Input_Parameters>::
Implicit_Euler(const int &cell_index_i,const int &cell_index_j, DenseMatrix* Jacobian)
{   
  //Low Mach # Preconditioning 
  if(Input_Parameters->Preconditioning){

    DenseMatrix Low_Mach_Number_Preconditioner(blocksize,blocksize,ZERO);         
    
    double delta_n = min( TWO*(SolnBlk->Grid.Cell[cell_index_i][cell_index_j].A/
			       (SolnBlk->Grid.lfaceE(cell_index_i, cell_index_j)
				+ SolnBlk->Grid.lfaceW(cell_index_i, cell_index_j))),
			  TWO*(SolnBlk->Grid.Cell[cell_index_i][cell_index_j].A/
			       (SolnBlk->Grid.lfaceN(cell_index_i, cell_index_j)
				+SolnBlk->Grid.lfaceS(cell_index_i, cell_index_j))));         
    
    SolnBlk->Uo[cell_index_i][cell_index_j].
      Low_Mach_Number_Preconditioner(Low_Mach_Number_Preconditioner,
				     SolnBlk->Flow_Type,
				     delta_n);    

    Jacobian[CENTER] -= Low_Mach_Number_Preconditioner/(SolnBlk->dt[cell_index_i][cell_index_j]);

  } else { // I/deltat

    DenseMatrix II(blocksize,blocksize);  
    II.identity();    
    Jacobian[CENTER] -= (II / (SolnBlk->dt[cell_index_i][cell_index_j]));
  }

}


/*!**************************************************************
 * Specialization of Block_Preconditioner::Preconditioner_dFIdU *
 *                                                              *
 * Calculates the dFdU matrix used to generate the approximate  *               
 * Jacobian for the Block Preconditioner.                       *
 ****************************************************************/
template<> inline void Block_Preconditioner<LESPremixed2D_pState,
					    LESPremixed2D_Quad_Block,					    
					    LESPremixed2D_Input_Parameters>::
Preconditioner_dFIdU(DenseMatrix &_dFIdU, LESPremixed2D_pState W)
{  
  dFIdU(_dFIdU,W,SolnBlk->Flow_Type);
}


/*!**************************************************************
 * Specialization of Block_Preconditioner::Preconditioner_dFIdU *
 *                                                              *
 * Calculates the dFdU matrix used to generate the approximate  *               
 * Jacobian for the Block Preconditioner.                       *
 ****************************************************************/ 
template<> inline void Block_Preconditioner<LESPremixed2D_pState,
					    LESPremixed2D_Quad_Block,			    
					    LESPremixed2D_Input_Parameters>::
Preconditioner_dFIdU_Roe(DenseMatrix &_dFIdU, int ii, int jj, int Orient)
{ 
#ifdef THICKENED_FLAME_ON
  int NUM_VAR_LESPREMIXED2D = SolnBlk->NumVar()-3;   
#else
  int NUM_VAR_LESPREMIXED2D = SolnBlk->NumVar()-1;
#endif
  DenseMatrix dFI_dW(NUM_VAR_LESPREMIXED2D,NUM_VAR_LESPREMIXED2D,ZERO);
  
  dFIdW_Inviscid_ROE(dFI_dW, *SolnBlk,*Input_Parameters, ii,jj,Orient);  
  //dFIdW_Inviscid_ROE_FD(dFI_dW, *SolnBlk,*Input_Parameters, ii,jj,Orient);
 
  //transformation Jacobian 
  DenseMatrix dWdU(NUM_VAR_LESPREMIXED2D,NUM_VAR_LESPREMIXED2D,ZERO); 
  //transformation Jacobian  Wo == W here 
  SolnBlk->W[ii][jj].dWdU(dWdU, SolnBlk->Flow_Type);
  _dFIdU += dFI_dW*dWdU;

}

/*!**************************************************************
 * Specialization of Block_Preconditioner::Preconditioner_dFIdU *
 *                                                              *
 * Calculates the dFdU matrix used to generate the approximate  *               
 * Jacobian for the Block Preconditioner.                       *
 ****************************************************************/ 
template<> inline void Block_Preconditioner<LESPremixed2D_pState,
					    LESPremixed2D_Quad_Block,			    
					    LESPremixed2D_Input_Parameters>::
Preconditioner_dFIdU_AUSM_plus_up(DenseMatrix &_dFIdU, int ii, int jj, int Orient)
{   
#ifdef THICKENED_FLAME_ON
  int NUM_VAR_LESPREMIXED2D = SolnBlk->NumVar()-3;   
#else
  int NUM_VAR_LESPREMIXED2D = SolnBlk->NumVar()-1;
#endif
  DenseMatrix dFIdW(NUM_VAR_LESPREMIXED2D,NUM_VAR_LESPREMIXED2D,ZERO);
  
  dFIdW_Inviscid_AUSM_plus_up(dFIdW, *SolnBlk,*Input_Parameters, ii,jj,Orient);
  
  //transformation Jacobian 
  DenseMatrix dWdU(NUM_VAR_LESPREMIXED2D,NUM_VAR_LESPREMIXED2D,ZERO);   
  SolnBlk->W[ii][jj].dWdU(dWdU, SolnBlk->Flow_Type);
  _dFIdU += dFIdW*dWdU;

}


/*!**************************************************************
 * Specialization of Block_Preconditioner::Preconditioner_dFVdU *
 *                                                              *
 * Calculates the dFdU matrix used to generate the approximate  *               
 * Jacobian for the Block Preconditioner.                       *
 ****************************************************************/
template<> inline void Block_Preconditioner<LESPremixed2D_pState,
					    LESPremixed2D_Quad_Block,					    
					    LESPremixed2D_Input_Parameters>::
Preconditioner_dFVdU(DenseMatrix &dFvdU, const int Rii, const int Rjj, 
		     const int Wii, const int Wjj, const int Orient_face, const int Orient_cell) {

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

#ifdef THICKENED_FLAME_ON
  int NUM_VAR_LESPREMIXED2D = SolnBlk->NumVar()-3;  
#else
  int NUM_VAR_LESPREMIXED2D = SolnBlk->NumVar()-1;
#endif
  int ns = SolnBlk->W[Rii][Rjj].ns-1; 
  //int Matrix_size = 2*NUM_VAR_LESPREMIXED2D+2;      // 14+Ns ?????  for variable R,k,mu, ie functions of ci see dRdU.cc
  int Matrix_size = 14 + ns;

  DenseMatrix dFvdWf(NUM_VAR_LESPREMIXED2D, Matrix_size,ZERO);
  DenseMatrix dWfdWx(Matrix_size, NUM_VAR_LESPREMIXED2D,ZERO); 
  DenseMatrix dGvdWf(NUM_VAR_LESPREMIXED2D, Matrix_size,ZERO);
  DenseMatrix dWfdWy(Matrix_size, NUM_VAR_LESPREMIXED2D,ZERO);  
  DenseMatrix dGVdW(NUM_VAR_LESPREMIXED2D,NUM_VAR_LESPREMIXED2D,ZERO);

  dFvdWf_Diamond(dFvdWf,dGvdWf,*SolnBlk, Orient_face, Rii, Rjj);
  dWfdWc_Diamond(dWfdWx,dWfdWy,*SolnBlk, Orient_face, Rii, Rjj, Orient_cell); 
  
  dGVdW = lface * (nface.x*(dFvdWf*dWfdWx) + nface.y*(dGvdWf*dWfdWy));
    
  //transformation Jacobian
  DenseMatrix dWdU(NUM_VAR_LESPREMIXED2D,NUM_VAR_LESPREMIXED2D,ZERO); 

  //transformation Jacobian  Wo == W here 
  SolnBlk->W[Wii][Wjj].dWdU(dWdU, SolnBlk->Flow_Type);  
  dFvdU += dGVdW*dWdU;


}

/*!**************************************************************
 * Specialization of Block_Preconditioner::Preconditioner_dSdU  *
 *                                                              *
 * Calculates the dFdU matrix used to generate the approximate  *               
 * Jacobian for the Block Preconditioner.                       *
 ****************************************************************/
template<> inline void Block_Preconditioner<LESPremixed2D_pState,
					    LESPremixed2D_Quad_Block,					    
					    LESPremixed2D_Input_Parameters>::
Preconditioner_dSdU(int ii, int jj, DenseMatrix &dRdU){
  SemiImplicitBlockJacobi(dRdU,*SolnBlk,IMPLICIT,ii,jj);
}


/*!**************************************************************
 *  Specialization of Block_Preconditioner::                    *
 *                               normalize_Preconditioner_dFdU  *
 *                                                              *
 * Normalizes the dFdU matrix used to generate the approximate  *               
 * Jacobian for the Block Preconditioner.                       *
 ****************************************************************/
void normalize_preconditioner(DenseMatrix &dFdU) 
{ 
 
  LESPremixed2D_pState W_STD_ATM;
  double ao  = W_STD_ATM.a();
  double rho = W_STD_ATM.rho;

  // ORIGINAL 
  dFdU(0,0) *= (ONE/ao); 
  dFdU(0,3) *=  ao;  
  dFdU(1,0) *= (ONE/(ao*ao));
  dFdU(1,1) *= (ONE/ao);
  dFdU(1,2) *= (ONE/ao);
  dFdU(2,0) *= (ONE/(ao*ao));
  dFdU(2,1) *= (ONE/ao);
  dFdU(2,2) *= (ONE/ao);
  dFdU(3,0) *= (ONE/(ao*ao*ao));
  dFdU(3,1) *= (ONE/(ao*ao));
  dFdU(3,2) *= (ONE/(ao*ao));
  dFdU(3,3) *= (ONE/ao);

  //scalars and cs's all have same normalization.
  for(int i=NUM_LESPREMIXED2D_VAR_SANS_SPECIES; 
      i< NUM_LESPREMIXED2D_VAR_SANS_SPECIES + W_STD_ATM.nscal + W_STD_ATM.ns-1; i++){   		  
    dFdU(0,i) *= (ONE/ao);
    dFdU(1,i) *= (ONE/(ao*ao));
    dFdU(2,i) *= (ONE/(ao*ao));
    dFdU(3,i) *= (ONE/(ao*ao*ao));
    dFdU(i,0) *= (ONE/ao);
    dFdU(i,3) *= ao;
    for(int j=NUM_LESPREMIXED2D_VAR_SANS_SPECIES; j< NUM_LESPREMIXED2D_VAR_SANS_SPECIES + W_STD_ATM.nscal + W_STD_ATM.ns-1; j++){   
      dFdU(i,j) *= (ONE/ao);            
    }
  } 

}

template<> inline void Block_Preconditioner<LESPremixed2D_pState,
					    LESPremixed2D_Quad_Block,
					    LESPremixed2D_Input_Parameters>::
normalize_Preconditioner_dFdU(DenseMatrix &dFdU) 
{ 
  normalize_preconditioner(dFdU); 
}


/************************************************************************/
/************ GMRES REQUIRED SPECIALIZATIONS & FUNCTIONS ****************/
/************************************************************************/   

/**************************************************************************
 * Routine: calculate_pertubed_residual                                   *
 **************************************************************************/
// Calculate SolnBlk.U =  SolnBlk.Uo + denormalize( epsilon * W(i) )
template <>inline void GMRES_Block<LESPremixed2D_pState,
				   LESPremixed2D_Quad_Block,
				   LESPremixed2D_Input_Parameters>::
calculate_perturbed_residual(const double &epsilon)
{    
  for (int j = JCl - Nghost ; j <= JCu + Nghost ; j++) {  //includes ghost cells 
    for (int i = ICl - Nghost ; i <= ICu + Nghost ; i++) {         
      for(int varindex = 0; varindex < blocksize; varindex++){	
	SolnBlk->U[i][j][varindex+1] = SolnBlk->Uo[i][j][varindex+1] +
	  denormalizeU( epsilon*W[search_directions*scalar_dim+index(i,j,varindex)], varindex);	    	
      }     
      //LESPremixed2D spec_check to make sure species (Uo + epsilon*W(i)) > ZERO , ie physical for dUdt calc 
      if(!SolnBlk->U[i][j].negative_speccheck(10)) { cerr<<"\n FAILURE in calculate_perturbed_residual"; exit(1); }       
      /* Update primitive variables. */
      SolnBlk->W[i][j] = SolnBlk->U[i][j].W();      
    }
  }  
}

// Calculate SolnBlk.U =  SolnBlk.Uo - denormalize( epsilon * W(i) )
template <>inline void GMRES_Block<LESPremixed2D_pState,
				   LESPremixed2D_Quad_Block,
				   LESPremixed2D_Input_Parameters>::
calculate_perturbed_residual_2nd(const double &epsilon)
{    
  for (int j = JCl - Nghost ; j <= JCu + Nghost ; j++) {  //includes ghost cells 
    for (int i = ICl - Nghost ; i <= ICu + Nghost ; i++) {      
      
      //copy back R + epsilon * W(i)
      SolnBlk->dUdt[i][j][1] = SolnBlk->dUdt[i][j][0];

      for(int varindex = 0; varindex < blocksize; varindex++){	
	SolnBlk->U[i][j][varindex+1] = SolnBlk->Uo[i][j][varindex+1] -
	  denormalizeU( epsilon*W[search_directions*scalar_dim+index(i,j,varindex)], varindex);	    	
      }     
      //LESPremixed2D spec_check to make sure species (Uo + epsilon*W(i)) > ZERO , ie physical for dUdt calc 
      if(!SolnBlk->U[i][j].negative_speccheck(10)) { cerr<<"\n FAILURE in calculate_perturbed_residual"; exit(1); }       
      /* Update primitive variables. */
      SolnBlk->W[i][j] = SolnBlk->U[i][j].W();      
    }
  }  
}

// Calculate SolnBlk.U =  SolnBlk.Uo + denormalize( epsilon * x(i) )
template <>inline void GMRES_Block<LESPremixed2D_pState,
				   LESPremixed2D_Quad_Block,
				   LESPremixed2D_Input_Parameters>::
calculate_perturbed_residual_Restart(const double &epsilon)
{    
  for (int j = JCl - Nghost ; j <= JCu + Nghost ; j++) {
    for (int i = ICl - Nghost ; i <= ICu + Nghost ; i++) {
      for(int varindex = 0; varindex < blocksize; varindex++){	
	SolnBlk->U[i][j][varindex+1] = SolnBlk->Uo[i][j][varindex+1] + denormalizeU( epsilon*x[index(i,j,varindex)], varindex);
      }  
      if(!SolnBlk->U[i][j].negative_speccheck(10)) { cerr<<"\n FAILURE in calculate_perturbed_residual_Restart "; exit(1); }
      /* Update primitive variables. */
      SolnBlk->W[i][j] = SolnBlk->U[i][j].W();      
    }
  }  
}

// Calculate SolnBlk.U =  SolnBlk.Uo + denormalize( epsilon * x(i) )
template <>inline void GMRES_Block<LESPremixed2D_pState,
				   LESPremixed2D_Quad_Block,
				   LESPremixed2D_Input_Parameters>::
calculate_perturbed_residual_2nd_Restart(const double &epsilon)
{    
  for (int j = JCl - Nghost ; j <= JCu + Nghost ; j++) {
    for (int i = ICl - Nghost ; i <= ICu + Nghost ; i++) {   
      
      //copy back R + epsilon * W(i)
      SolnBlk->dUdt[i][j][1] = SolnBlk->dUdt[i][j][0];
      
      for(int varindex = 0; varindex < blocksize; varindex++){	
	SolnBlk->U[i][j][varindex+1] = SolnBlk->Uo[i][j][varindex+1] - denormalizeU( epsilon*x[index(i,j,varindex)], varindex);
      }  
      if(!SolnBlk->U[i][j].negative_speccheck(10)) { cerr<<"\n FAILURE in calculate_perturbed_residual_Restart "; exit(1); }
      /* Update primitive variables. */
      SolnBlk->W[i][j] = SolnBlk->U[i][j].W();      
    }
  }  
}


/********************************************************
 * Routine: calculate_Matrix_Free with Preconditioning  *
 ********************************************************/
template <>inline void GMRES_Block<LESPremixed2D_pState,
				   LESPremixed2D_Quad_Block,
				   LESPremixed2D_Input_Parameters>::
calculate_Matrix_Free(const double &epsilon)
{
  //Taking into acount NKS overlap
  int JCl_overlap = 0; int JCu_overlap = 0;
  int ICu_overlap = 0; int ICl_overlap = 0;		  
  double value = ONE; 

  if(overlap){	
    if ( SolnBlk->Grid.BCtypeS[ICl] == BC_NONE)  JCl_overlap = overlap; 
    if ( SolnBlk->Grid.BCtypeN[ICu] == BC_NONE)  JCu_overlap = overlap;
    if ( SolnBlk->Grid.BCtypeE[JCu] == BC_NONE)  ICu_overlap = overlap;
    if ( SolnBlk->Grid.BCtypeW[JCl] == BC_NONE)  ICl_overlap = overlap;
  }
	   
  DenseMatrix Precon(blocksize,blocksize,ZERO); //SHOULD MOVE THIS OUT OF HERE AND STORE SOMEWHERE TO AVOID RECREATING EACH TIME
  
  // Non-Overlap Ghost Cells R(U) already set to zero by dUdt calculation  
  /* V(i+1) = ( R(U+epsilon*W) - b) / epsilon - (gamma) z / h */
  for (int j = JCl - JCl_overlap; j <= JCu + JCu_overlap; j++) {
    for (int i = ICl - ICl_overlap; i <= ICu + ICu_overlap; i++) {
      
      // Low Mach # Preconditioner                   
      if(Input_Parameters->Preconditioning){ 
	double delta_n = min( TWO*(SolnBlk->Grid.Cell[i][j].A/(SolnBlk->Grid.lfaceE(i, j) + SolnBlk->Grid.lfaceW(i, j))),
			      TWO*(SolnBlk->Grid.Cell[i][j].A/(SolnBlk->Grid.lfaceN(i, j) + SolnBlk->Grid.lfaceS(i, j))));             
	SolnBlk->Uo[i][j].Low_Mach_Number_Preconditioner(Precon, SolnBlk->Flow_Type, delta_n);       
      }

      //Update V 
      for(int k =0; k < blocksize; k++){	
	int iter = index(i,j,k);		
	//Matrix Free V(i+1)  
	
	if( Input_Parameters->NKS_IP.GMRES_Frechet_Derivative_Order == FIRST_ORDER ){
	  //forwards differenceing R(U+epsilon) - R(U) / epsilon
	  V[(search_directions+1)*scalar_dim+iter] = (normalizeR(SolnBlk->dUdt[i][j][0][k+1],k) - b[iter]) / epsilon ;
	} else if ( Input_Parameters->NKS_IP.GMRES_Frechet_Derivative_Order == SECOND_ORDER ){
	  //2nd order R(U+epsilon) - R(U-epsilon) / 2*epsilon
	  V[(search_directions+1)*scalar_dim+iter] = normalizeR( SolnBlk->dUdt[i][j][1][k+1] - SolnBlk->dUdt[i][j][0][k+1],k)/(TWO*epsilon);
	}

	//Finite Time Stepping
	if(Input_Parameters->Preconditioning){   //gamma(nxn)*z(nx1)/h(1x)    
	  value = ZERO;
	  for(int l =0; l < blocksize; l++){
	    value += Precon(k,l) * denormalizeU(W[(search_directions)*scalar_dim + index(i,j,l)],l);
	  }
	  V[(search_directions+1)*scalar_dim+iter] -= normalizeR(value/(SolnBlk->dt[i][j]),k);
	  
	//No Preconditioner
	} else { // z/h
	  V[(search_directions+1)*scalar_dim+iter] -= normalizeUtoR(W[(search_directions)*scalar_dim + iter]/(SolnBlk->dt[i][j]),k);	
	}
      }  
    
    } 
  } 
}

/****************************************************************
 * Routine: calculate_Matrix_Free_Restart with Preconditioning  *
 ****************************************************************/
template <>inline void GMRES_Block<LESPremixed2D_pState,
				   LESPremixed2D_Quad_Block,
				   LESPremixed2D_Input_Parameters>::
calculate_Matrix_Free_Restart(const double &epsilon)
{
  //Taking into acount NKS overlap
  int JCl_overlap = 0; int JCu_overlap = 0;
  int ICu_overlap = 0; int ICl_overlap = 0;		  
  double value = ONE; 

  if(overlap){	
    if ( SolnBlk->Grid.BCtypeS[ICl] == BC_NONE)  JCl_overlap = overlap; 
    if ( SolnBlk->Grid.BCtypeN[ICu] == BC_NONE)  JCu_overlap = overlap;
    if ( SolnBlk->Grid.BCtypeE[JCu] == BC_NONE)  ICu_overlap = overlap;
    if ( SolnBlk->Grid.BCtypeW[JCl] == BC_NONE)  ICl_overlap = overlap;
  }
	   
  DenseMatrix Precon(blocksize,blocksize,ZERO);
  
  // Non-Overlap Ghost Cells R(U) already set to zero by dUdt calculation  
  /* V(i+1) = ( R(U+epsilon*x) - b) / epsilon - (gamma) x / h */
  for (int j = JCl - JCl_overlap; j <= JCu + JCu_overlap; j++) {
    for (int i = ICl - ICl_overlap; i <= ICu + ICu_overlap; i++) {
      
      // Low Mach # Preconditioner                   
      if(Input_Parameters->Preconditioning){ 
	double delta_n = min( TWO*(SolnBlk->Grid.Cell[i][j].A/(SolnBlk->Grid.lfaceE(i, j) + SolnBlk->Grid.lfaceW(i, j))),
			      TWO*(SolnBlk->Grid.Cell[i][j].A/(SolnBlk->Grid.lfaceN(i, j) + SolnBlk->Grid.lfaceS(i, j))));             
	SolnBlk->Uo[i][j].Low_Mach_Number_Preconditioner(Precon, SolnBlk->Flow_Type, delta_n);       
      }

      //Update V 
      for(int k =0; k < blocksize; k++){	
	int iter = index(i,j,k);		
	//Matrix Free V(i+1) 

	if( Input_Parameters->NKS_IP.GMRES_Frechet_Derivative_Order == FIRST_ORDER ){
	  V[iter] = (normalizeR(SolnBlk->dUdt[i][j][0][k+1],k) - b[iter]) / epsilon ;
	} else if ( Input_Parameters->NKS_IP.GMRES_Frechet_Derivative_Order == SECOND_ORDER ){
	  V[iter] = normalizeR( SolnBlk->dUdt[i][j][1][k+1] - SolnBlk->dUdt[i][j][0][k+1],k)/(TWO*epsilon);
	}

	//Finite Time Stepping
	if(Input_Parameters->Preconditioning){   //gamma(nxn)*x(nx1)/h(1x)    
	  value = ZERO;
	  for(int l =0; l < blocksize; l++){
	    value += Precon(k,l) * denormalizeU( x[index(i,j,l)],l);
	  }
	  V[iter] -= normalizeR(value/(SolnBlk->dt[i][j]),k);
	  
	//No Preconditioner
	} else { // z/h
	  V[iter] -= normalizeUtoR(x[iter]/(SolnBlk->dt[i][j]),k);	
	}
      }  
    
    } 
  } 
}


/*!********************************************************
 * GMRES_Block::set_normalize_valuesR for NKS/GMRES       * 
 *                                                        *
 *   normalize_values[0] must be set to ao                *
 *   normalize_values[1-n] = values for cstate index[1-n] *
 *             premultiplied by ao ie for "R"             *
 *             where n = the number of solution variables *
 **********************************************************/
template<> inline void GMRES_Block<LESPremixed2D_pState,
				   LESPremixed2D_Quad_Block,
				   LESPremixed2D_Input_Parameters>::
set_normalize_values(void)
{   

  LESPremixed2D_pState W_STD_ATM;
  double ao  = W_STD_ATM.a();
  double rho = W_STD_ATM.rho;

  // Original Normalization from code  
  normalize_valuesU[0] = rho;          //rho
  normalize_valuesU[1] = rho*ao;       //rho*u
  normalize_valuesU[2] = rho*ao;       //rho*v
  normalize_valuesU[3] = rho*ao*ao;    //rho*e

  for(int i=NUM_LESPREMIXED2D_VAR_SANS_SPECIES; i < NUM_LESPREMIXED2D_VAR_SANS_SPECIES + W_STD_ATM.nscal; i++){
    normalize_valuesU[i] = rho;        //rhoscalar[i] 
  }
  for(int i=NUM_LESPREMIXED2D_VAR_SANS_SPECIES + W_STD_ATM.nscal; i < blocksize; i++){
    normalize_valuesU[i] = rho;        //species mass fraction rho*cs
  }



  normalize_valuesR[0] = rho*ao;          //rho
  normalize_valuesR[1] = rho*ao*ao;       //rho*u
  normalize_valuesR[2] = rho*ao*ao;       //rho*v
  normalize_valuesR[3] = rho*ao*ao*ao;    //rho*e

  for(int i=NUM_LESPREMIXED2D_VAR_SANS_SPECIES; i < NUM_LESPREMIXED2D_VAR_SANS_SPECIES + W_STD_ATM.nscal; i++){
    normalize_valuesR[i] = rho*ao;        //rhoscalar[i] 
  }
  for(int i=NUM_LESPREMIXED2D_VAR_SANS_SPECIES + W_STD_ATM.nscal; i < blocksize; i++){
    normalize_valuesR[i] = rho*ao;        //species mass fraction rho*cs
  }
 
}

#endif // _LESPREMIXED2D_NKS_INCLUDED 







#ifndef _FLAME2D_NKS_INCLUDED 
#define _FLAME2D_NKS_INCLUDED 

/****************************************************************/
/****** NKS REQUIRED SPECIALIZATIONS & FUNCTIONS ****************/
/****************************************************************/  
#include "Flame2DQuad.h"
#include "Flame2DdRdU.h"
#include "../NewtonKrylovSchwarz2D/NKS2D.h"

/****************************************************************/



/*! *************************************************************
 *  Specialization of Newton_Update                             *
 *                                                              *
 * This routine updates the previous Solution Data Uo with the  *
 * deltaU from the GMRES                                        *
 * iterative solver.  U = Uo + GMRES.delatU                     *
 *                                                              *
 ****************************************************************/
template <>
int Newton_Update(Flame2D_Quad_Block *SolnBlk,
		  AdaptiveBlock2D_List &List_of_Local_Solution_Blocks,
		  Flame2D_Input_Parameters &Input_Parameters,
		  GMRES_RightPrecon_MatrixFree<Flame2D_pState,
		  Flame2D_Quad_Block,
		  Flame2D_Input_Parameters> &GMRES,
		  double Relaxation_multiplier) {

  // declares
  const int Num_Var = SolnBlk[0].NumVar();  	
  int error_flag = 0;
  bool isGoodState;
  
  //
  // Update Solution. No updates to Ghost Cells, let the BC's take care of it
  //
  for ( int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks.Nblk ; ++Bcount ) {
    if (List_of_Local_Solution_Blocks.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
      for (int j = SolnBlk[Bcount].JCl; j <= SolnBlk[Bcount].JCu; j++){
	for (int i = SolnBlk[Bcount].ICl; i <= SolnBlk[Bcount].ICu; i++){
	  
	  // Update solutions in conserved variables  
	  // U = Uo + RELAXATION*deltaU = Uo + denormalized(x)
	  for(int varindex =1; varindex <= Num_Var; varindex++){	                       
	    SolnBlk[Bcount].U[i][j][varindex] = SolnBlk[Bcount].Uo[i][j][varindex]  +  
	      Relaxation_multiplier*GMRES.deltaU(Bcount,i,j,varindex-1);
	  } 	      	  
	  
	  // check species
	  isGoodState = SolnBlk[Bcount].U[i][j].isPhysical(10);

	  /*********************************************************/

	  // Apply update reduction while any one of the updated variables is unphysical 
	  if(!isGoodState){	   
 	    
	    double update_reduction_factor = ONE;	    
	    for (int n_update_reduction = 1; n_update_reduction <= 10; ++n_update_reduction) {
	      update_reduction_factor = HALF*update_reduction_factor;		  		  
	      for(int varindex = 1; varindex <= Num_Var; varindex++){
		SolnBlk[Bcount].U[i][j][varindex] = SolnBlk[Bcount].Uo[i][j][varindex] 
		  + GMRES.deltaU(Bcount,i,j,varindex-1)*update_reduction_factor;
	      }   
	      cout<<"\n Applying Reduction to solution in NKS "<<n_update_reduction;

	      isGoodState = SolnBlk[Bcount].U[i][j].isPhysical(n_update_reduction);
	      if(isGoodState)  break;	      
	    } 

	  } // endif - isGoodState
	 
	    /*********************************************************/

	    // Error Check
	  if(! isGoodState) error_flag = 1;
	 	  
	  //Update solution in primitive variables.
	  SolnBlk[Bcount].W[i][j].setU( SolnBlk[Bcount].U[i][j] );
	    
	} 
      }
      
// #ifdef _NKS_VERBOSE  
//       if (CFFC_Primary_MPI_Processor()) { 
//      /**************************************************************************/
// 	//FOR DEBUGGING OUTPUT GMRES deltaU set to "SCALED" residual
// 	double *norm = new double[Num_Var-1];
// 	for(int i= 0; i< Num_Var-1; i++) norm[i]=ZERO;
	
      for (int j = SolnBlk[Bcount].JCl-SolnBlk[Bcount].Nghost; j <= SolnBlk[Bcount].JCu+SolnBlk[Bcount].Nghost; j++){
	for (int i = SolnBlk[Bcount].ICl-SolnBlk[Bcount].Nghost; i <= SolnBlk[Bcount].ICu+SolnBlk[Bcount].Nghost; i++){
	  for(int varindex =1; varindex <= Num_Var; varindex++){	  	   
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
//      /**************************************************************************/
//       }
// #endif

    } 
  }   
  return error_flag; 
}

/*! *************************************************************
 * Flame2D Specific  Finite_Time_Step                           *
 *                                                              *
 * This routine calculates the Finite Time Step using a basic   *
 * SER approach for startup,                                    *
 * returning a "CFL" number to multiply by the stability        *
 * determined dt.                                               *
 *                                                              *
 ****************************************************************/
template <> double Finite_Time_Step(const Flame2D_Input_Parameters &Input_Parameters, 
				    const double &L2norm_first,
				    const double &L2norm_current,
				    const double &L2norm_current_n,		
				    const int &Number_of_Newton_Steps) {

  double CFL_current;

  //SER 
  if (L2norm_current_n > Input_Parameters.NKS_IP.Min_Finite_Time_Step_Norm_Ratio ) { 
    //Original works for Invisicd, Viscous OK, Reacting seems to be more picky
    CFL_current = Input_Parameters.NKS_IP.Finite_Time_Step_Initial_CFL*
      pow( max(ONE, ONE/L2norm_current_n),ONE ); 
    
    //Possible to help curb jumps, !!!!!!!!!!!!!!!! however can cause oscillations in residual!!!!!!!!!!!!!!!!
    if( L2norm_first == L2norm_current && Number_of_Newton_Steps != 1 &&
	CFL_current <= Input_Parameters.NKS_IP.Finite_Time_Step_Initial_CFL){ //ie still ONE
      //CFL_current = 0.1*CFL_current;      
    }

  } else {
    CFL_current = Input_Parameters.NKS_IP.Finite_Time_Step_Initial_CFL/
      Input_Parameters.NKS_IP.Min_Finite_Time_Step_Norm_Ratio;
  }  

  // CAN SET A MAXIMUM, SOMETIMES USEFUL WITH CONVERGENCE STALL & BC ISSUES
  if( CFL_current > Input_Parameters.NKS_IP.Finite_Time_Step_Max_CFL) 
    CFL_current = Input_Parameters.NKS_IP.Finite_Time_Step_Max_CFL;

  return CFL_current;

}


/*! *************************************************************
 * SPECIALIZATION of finite time step addition to diagonal.     *
 ****************************************************************/ 
template<> inline void Block_Preconditioner<Flame2D_pState,
					    Flame2D_Quad_Block,
					    Flame2D_Input_Parameters>::
Implicit_Euler(const int &cell_index_i,
	       const int &cell_index_j, 
	       DenseMatrix* Jacobian, 
	       const double& DTS_dTime)
{   
  //
  //Low Mach # Preconditioning 
  //
  if(Input_Parameters->Preconditioning){

    static Flame2D_pState Wo;
    static DenseMatrix Low_Mach_Number_Preconditioner(blocksize,blocksize,ZERO);
    
    // spacing for preconditioner
    double delta_n( SolnBlk->delta_n( cell_index_i, cell_index_j ) );
    
    // get initial pState
    Wo.setU( SolnBlk->Uo[cell_index_i][cell_index_j] );
    Wo.Low_Mach_Number_Preconditioner(Low_Mach_Number_Preconditioner,
				      SolnBlk->Flow_Type,
				      delta_n);    

    Jacobian[CENTER] -= Low_Mach_Number_Preconditioner*
      LHS_Time<Flame2D_Input_Parameters>(*Input_Parameters, 
					 SolnBlk->dt[cell_index_i][cell_index_j],
					 DTS_dTime);


    //
    // I/deltat
    //
  } else {

    static DenseMatrix II(blocksize,blocksize);  
    II.identity();    
    Jacobian[CENTER] -= II*LHS_Time<Flame2D_Input_Parameters>(*Input_Parameters, 
							      SolnBlk->dt[cell_index_i][cell_index_j],
							      DTS_dTime);
  }

}

/*!**************************************************************
 * Specialization of Block_Preconditioner::Preconditioner_dFIdU *
 *                                                              *
 * Calculates the dFdU matrix used to generate the approximate  *               
 * Jacobian for the Block Preconditioner.                       *
 ****************************************************************/
template<> inline void Block_Preconditioner<Flame2D_pState,
					    Flame2D_Quad_Block,					    
					    Flame2D_Input_Parameters>::
Preconditioner_dFIdU(DenseMatrix &_dFIdU, Flame2D_pState W)
{  
  //W.dFIdU(_dFIdU);
  cerr << "Preconditioner_dFIdU(): Function note used";
  exit(-1);
}

/*! *************************************************************
 *  Calculate First Order Local Jacobian Block(s) Coresponding  *
 *  to Cell(i,j) using HLLE                                     *
 ****************************************************************/ 
template<> inline void Block_Preconditioner<Flame2D_pState,
					    Flame2D_Quad_Block,
					    Flame2D_Input_Parameters>::
First_Order_Inviscid_Jacobian_HLLE(const int &cell_index_i,
				   const int &cell_index_j, 
				   DenseMatrix* Jacobian){              

  //! Caculate normal vectors -> in Vector2D format. 
  static Vector2D nface_N; 
  static Vector2D nface_S; 
  static Vector2D nface_E; 
  static Vector2D nface_W; 
  nface_N = SolnBlk->Grid.nfaceN(cell_index_i,cell_index_j-1);
  nface_S = SolnBlk->Grid.nfaceS(cell_index_i,cell_index_j+1);
  nface_E = SolnBlk->Grid.nfaceE(cell_index_i-1,cell_index_j);
  nface_W = SolnBlk->Grid.nfaceW(cell_index_i+1,cell_index_j);

  //! Calculate wavespeeds using solutions in the rotated frame 
  //! -> in Vector2D format.
  static Vector2D lambdas_N;
  Flame2D_pState::HLLE_wavespeeds(SolnBlk->W[cell_index_i][cell_index_j-1], 
				  SolnBlk->W[cell_index_i][cell_index_j], 
				  nface_N, lambdas_N);
  static Vector2D lambdas_S;
  Flame2D_pState::HLLE_wavespeeds(SolnBlk->W[cell_index_i][cell_index_j+1], 
				  SolnBlk->W[cell_index_i][cell_index_j], 
				  nface_S, lambdas_S);  
  static Vector2D lambdas_E;
  Flame2D_pState::HLLE_wavespeeds(SolnBlk->W[cell_index_i-1][cell_index_j], 
				  SolnBlk->W[cell_index_i][cell_index_j], 
				  nface_E, lambdas_E);
  static Vector2D lambdas_W;
  Flame2D_pState::HLLE_wavespeeds(SolnBlk->W[cell_index_i+1][cell_index_j], 
				  SolnBlk->W[cell_index_i][cell_index_j], 
				  nface_W, lambdas_W);

  //checks necessary ????
  //if ((lambdas_W.y-lambdas_W.x) == ZERO) 
  //  cout << " WEST : HLLE_wavespeeds " << endl;
 
  //! Calculate constants gamma and beta -> scalar values. 
  double gamma_N( (lambdas_N.x*lambdas_N.y)/(lambdas_N.y-lambdas_N.x) );
  double beta_N( - lambdas_N.x/(lambdas_N.y-lambdas_N.x) );
  double gamma_S( (lambdas_S.x*lambdas_S.y)/(lambdas_S.y-lambdas_S.x) );
  double beta_S( - lambdas_S.x/(lambdas_S.y-lambdas_S.x) );
  double gamma_E( (lambdas_E.x*lambdas_E.y)/(lambdas_E.y-lambdas_E.x) );
  double beta_E( - lambdas_E.x/(lambdas_E.y-lambdas_E.x) );
  double gamma_W( (lambdas_W.x*lambdas_W.y)/(lambdas_W.y-lambdas_W.x) );
  double beta_W( - lambdas_W.x/(lambdas_W.y-lambdas_W.x) );

  //! Obtain rotation matrices with normal vector -> matrices in 
  //! DenseMatrix format. 
  static DenseMatrix A_N(blocksize,blocksize);
  static DenseMatrix AI_N(blocksize,blocksize);
  static DenseMatrix A_S(blocksize,blocksize);
  static DenseMatrix AI_S(blocksize,blocksize);
  static DenseMatrix A_E(blocksize,blocksize);
  static DenseMatrix AI_E(blocksize,blocksize);
  static DenseMatrix A_W(blocksize,blocksize);
  static DenseMatrix AI_W(blocksize,blocksize);
  RotationMatrix(A_N, nface_N, 1);
  RotationMatrix(AI_N, nface_N, 0);
  RotationMatrix(A_S, nface_S, 1);
  RotationMatrix(AI_S, nface_S, 0);
  RotationMatrix(A_E, nface_E, 1);
  RotationMatrix(AI_E, nface_E, 0);
  RotationMatrix(A_W, nface_W, 1);
  RotationMatrix(AI_W, nface_W, 0);

  //Solution Rotate provided in pState 
  // Rotate in place
  double 
    u( ((const Flame2D_pState&)SolnBlk->W[cell_index_i][cell_index_j]).vx() ), 
    v( ((const Flame2D_pState&)SolnBlk->W[cell_index_i][cell_index_j]).vy() );

  SolnBlk->W[cell_index_i][cell_index_j].Rotate(nface_N);
  SolnBlk->W[cell_index_i][cell_index_j].dFIdU(Jacobian[NORTH]); 
  SolnBlk->W[cell_index_i][cell_index_j].setVelocity(u, v);

  SolnBlk->W[cell_index_i][cell_index_j].Rotate(nface_S);
  SolnBlk->W[cell_index_i][cell_index_j].dFIdU(Jacobian[SOUTH]);
  SolnBlk->W[cell_index_i][cell_index_j].setVelocity(u, v);

  SolnBlk->W[cell_index_i][cell_index_j].Rotate(nface_E);
  SolnBlk->W[cell_index_i][cell_index_j].dFIdU(Jacobian[EAST]);
  SolnBlk->W[cell_index_i][cell_index_j].setVelocity(u, v);

  SolnBlk->W[cell_index_i][cell_index_j].Rotate(nface_W);
  SolnBlk->W[cell_index_i][cell_index_j].dFIdU(Jacobian[WEST]);
  SolnBlk->W[cell_index_i][cell_index_j].setVelocity(u, v);

  //! Calculate Jacobian matrix -> blocksizexblocksize matrix in 
  //! DenseMatrix format
  //  NOTE: NORTH, SOUTH, EAST, and WEST jacobians are zero at this point,
  //        CENTER already has terms in it
  static DenseMatrix II(blocksize,blocksize);  II.identity();    
  //North
  Jacobian[NORTH] *= SolnBlk->Grid.lfaceN(cell_index_i,cell_index_j-1) * beta_N;
  Jacobian[NORTH] += (SolnBlk->Grid.lfaceN(cell_index_i,cell_index_j-1) * gamma_N) * II;
  Jacobian[NORTH]  = AI_N * Jacobian[NORTH] * A_N;

  //South
  Jacobian[SOUTH] *= SolnBlk->Grid.lfaceS(cell_index_i,cell_index_j+1) * beta_S;
  Jacobian[SOUTH] += (SolnBlk->Grid.lfaceS(cell_index_i,cell_index_j+1) * gamma_S) * II;
  Jacobian[SOUTH]  = AI_S * Jacobian[SOUTH] * A_S;

  //East
  Jacobian[EAST] *= SolnBlk->Grid.lfaceE(cell_index_i-1,cell_index_j) * beta_E;
  Jacobian[EAST] += (SolnBlk->Grid.lfaceE(cell_index_i-1,cell_index_j) * gamma_E) * II;
  Jacobian[EAST]  = AI_E * Jacobian[EAST] * A_E;

  //West
  Jacobian[WEST] *= SolnBlk->Grid.lfaceW(cell_index_i+1,cell_index_j) * beta_W;
  Jacobian[WEST] += (SolnBlk->Grid.lfaceW(cell_index_i+1,cell_index_j) * gamma_W) * II;
  Jacobian[WEST]  = AI_W * Jacobian[WEST] * A_W;

  //Center calculated from neighbours
  //! Using the fact that dF/dU(right) = - dF/dU(left)
  for (int i=0; i<blocksize; i++)
    for (int j=0; j<blocksize; j++) {
      Jacobian[CENTER](i,j) += ( (Jacobian[NORTH](i,j) + Jacobian[SOUTH](i,j) + 
				  Jacobian[EAST](i,j)  + Jacobian[WEST](i,j)) /
				 SolnBlk->Grid.Cell[cell_index_i][cell_index_j].A );
      
      Jacobian[NORTH](i,j) /= -SolnBlk->Grid.Cell[cell_index_i][cell_index_j-1].A;
      Jacobian[SOUTH](i,j) /= -SolnBlk->Grid.Cell[cell_index_i][cell_index_j+1].A;
      Jacobian[EAST](i,j) /= -SolnBlk->Grid.Cell[cell_index_i-1][cell_index_j].A;
      Jacobian[WEST](i,j) /= -SolnBlk->Grid.Cell[cell_index_i+1][cell_index_j].A;
    }
}

/*!**************************************************************
 * Specialization of Block_Preconditioner::Preconditioner_dFIdU *
 *                                                              *
 * Calculates the dFdU matrix used to generate the approximate  *               
 * Jacobian for the Block Preconditioner.                       *
 ****************************************************************/ 
template<> inline void Block_Preconditioner<Flame2D_pState,
					    Flame2D_Quad_Block,
					    Flame2D_Input_Parameters>::
Preconditioner_dFIdU_Roe(DenseMatrix &_dFIdU, 
			 int ii, int jj, 
			 int Orient)
{ 
  // declares
  static DenseMatrix dFIdW(blocksize,blocksize);
  dFIdW.zero();
  
  //! Calculate dFdU using solutions in the rotated frame -> matrix in DenseMatrix format. 
  dFIdW_Inviscid_ROE(dFIdW, *SolnBlk,*Input_Parameters, ii,jj,Orient);  
  //dFIdW_Inviscid_ROE_FD(dFIdW, *SolnBlk,*Input_Parameters, ii,jj,Orient);
 
  //transformation Jacobian 
  // Non neeed to zero more than once, always writing to the same spot
  static DenseMatrix dWdU(blocksize,blocksize,ZERO); 
  
  //transformation Jacobian  Wo == W here 
  SolnBlk->W[ii][jj].dWdU(dWdU);
  _dFIdU += dFIdW*dWdU;

}


/*!**************************************************************
 *  Calculate First Order Local Jacobian Block(s) Coresponding  *
 *  to Cell(i,j)                                                *
 *  using Roe                                                   *
 ****************************************************************/
template<> inline void Block_Preconditioner<Flame2D_pState,
					    Flame2D_Quad_Block,
					    Flame2D_Input_Parameters>::
First_Order_Inviscid_Jacobian_Roe(const int &cell_index_i,
				  const int &cell_index_j, 
				  DenseMatrix* Jacobian){              
    
  //! Calculate Jacobian matrix -> blocksizexblocksize matrix in 
  //! DenseMatrix format
  //  NOTE: NORTH, SOUTH, EAST, and WEST jacobians are zero at this point,
  //        CENTER already has terms in it
  Preconditioner_dFIdU_Roe(Jacobian[NORTH],cell_index_i,cell_index_j,NORTH);
  Preconditioner_dFIdU_Roe(Jacobian[SOUTH],cell_index_i,cell_index_j,SOUTH); 
  Preconditioner_dFIdU_Roe(Jacobian[EAST],cell_index_i,cell_index_j,EAST);
  Preconditioner_dFIdU_Roe(Jacobian[WEST],cell_index_i,cell_index_j,WEST); 

  //Center calculated from neighbours
  //! Using the fact that dF/dU(right) = - dF/dU(left)
  for (int i=0; i<blocksize; i++)
    for (int j=0; j<blocksize; j++) {
      Jacobian[CENTER](i,j) += ( (Jacobian[NORTH](i,j) + Jacobian[SOUTH](i,j) + 
				  Jacobian[EAST](i,j)  + Jacobian[WEST](i,j)) /
				 SolnBlk->Grid.Cell[cell_index_i][cell_index_j].A );
      
      Jacobian[NORTH](i,j) /= -SolnBlk->Grid.Cell[cell_index_i][cell_index_j-1].A;
      Jacobian[SOUTH](i,j) /= -SolnBlk->Grid.Cell[cell_index_i][cell_index_j+1].A;
      Jacobian[EAST](i,j) /= -SolnBlk->Grid.Cell[cell_index_i-1][cell_index_j].A;
      Jacobian[WEST](i,j) /= -SolnBlk->Grid.Cell[cell_index_i+1][cell_index_j].A;
    }

}

/*!**************************************************************
 * Specialization of Block_Preconditioner::Preconditioner_dFIdU *
 *                                                              *
 * Calculates the dFdU matrix used to generate the approximate  *               
 * Jacobian for the Block Preconditioner.                       *
 ****************************************************************/ 
template<> inline void Block_Preconditioner<Flame2D_pState,
					    Flame2D_Quad_Block,
					    Flame2D_Input_Parameters>::
Preconditioner_dFIdU_AUSM_plus_up(DenseMatrix &_dFIdU, 
				  int ii, int jj, 
				  int Orient)
{   
  // declares
  static DenseMatrix dFIdW(blocksize,blocksize);
  dFIdW.zero();
  
  //! Calculate dFdU using solutions in the rotated frame 
  //! -> matrix in DenseMatrix format. 
  dFIdW_Inviscid_AUSM_plus_up(dFIdW, *SolnBlk,*Input_Parameters, ii,jj,Orient);
  
  //transformation Jacobian 
  // Non neeed to zero more than once, always writing to the same spot
  static DenseMatrix dWdU(blocksize,blocksize,ZERO);

  //transformation Jacobian  Wo == W here 
  SolnBlk->W[ii][jj].dWdU(dWdU);
  _dFIdU += dFIdW*dWdU;

}

/*!**************************************************************
 *  Calculate First Order Local Jacobian Block(s) Coresponding  *
 *  to Cell(i,j) using AUSM_plus_up                             *
 ****************************************************************/
template<> inline void Block_Preconditioner<Flame2D_pState,
					    Flame2D_Quad_Block,
					    Flame2D_Input_Parameters>::
First_Order_Inviscid_Jacobian_AUSM_plus_up(const int &cell_index_i,
					   const int &cell_index_j, 
					   DenseMatrix* Jacobian){              

  //! Calculate Jacobian matrix -> blocksizexblocksize matrix in 
  //! DenseMatrix format
  //  NOTE: NORTH, SOUTH, EAST, and WEST jacobians are zero at this point,
  //        CENTER already has terms in it
  Preconditioner_dFIdU_AUSM_plus_up(Jacobian[NORTH],
				    cell_index_i,cell_index_j,NORTH);
  Preconditioner_dFIdU_AUSM_plus_up(Jacobian[SOUTH],
				    cell_index_i,cell_index_j,SOUTH); 
  Preconditioner_dFIdU_AUSM_plus_up(Jacobian[EAST],
				    cell_index_i,cell_index_j,EAST);        
  Preconditioner_dFIdU_AUSM_plus_up(Jacobian[WEST],
				    cell_index_i,cell_index_j,WEST); 

  //Center calculated from neighbours
  //! Using the fact that dF/dU(right) = - dF/dU(left)
  for (int i=0; i<blocksize; i++)
    for (int j=0; j<blocksize; j++) {
      Jacobian[CENTER](i,j) += ( (Jacobian[NORTH](i,j) + Jacobian[SOUTH](i,j) + 
				  Jacobian[EAST](i,j)  + Jacobian[WEST](i,j)) /
				 SolnBlk->Grid.Cell[cell_index_i][cell_index_j].A );
      
      Jacobian[NORTH](i,j) /= -SolnBlk->Grid.Cell[cell_index_i][cell_index_j-1].A;
      Jacobian[SOUTH](i,j) /= -SolnBlk->Grid.Cell[cell_index_i][cell_index_j+1].A;
      Jacobian[EAST](i,j) /= -SolnBlk->Grid.Cell[cell_index_i-1][cell_index_j].A;
      Jacobian[WEST](i,j) /= -SolnBlk->Grid.Cell[cell_index_i+1][cell_index_j].A;
    }
  
}

/*!**************************************************************
 * Specialization of Block_Preconditioner::Preconditioner_dFVdU *
 *                                                              *
 * Calculates the dFdU matrix used to generate the approximate  *               
 * Jacobian for the Block Preconditioner.                       *
 ****************************************************************/
template<> inline void Block_Preconditioner<Flame2D_pState,
					    Flame2D_Quad_Block,					    
					    Flame2D_Input_Parameters>::
Preconditioner_dFVdU(DenseMatrix &dFvdU, 
		     const int Rii, const int Rjj, 
		     const int Wii, const int Wjj, 
		     const int Orient_face, 
		     const int Orient_cell) {

  //-----------------------------------------------------------------
  // declares
  //-----------------------------------------------------------------
  // temporaries
  double lface;
  static Vector2D nface;

  // solution size
  const int ns = SolnBlk->W[Rii][Rjj].NumSpecies();
  const int Matrix_size = 10 + ns; // rho, u, v, p, drho, du/dx, du/dy, 
                                   // dv/dx, dv/dy, dp/dx, dcn/dx

  // temporary matrices
  // Only need to zero these matrices once as we are always writing to the same spot
  static DenseMatrix dFvdWf(blocksize, Matrix_size,ZERO);
  static DenseMatrix dGvdWf(blocksize, Matrix_size,ZERO); 
  static DenseMatrix dWfdWx(Matrix_size, blocksize,ZERO);
  static DenseMatrix dWfdWy(Matrix_size, blocksize,ZERO);
  static DenseMatrix dGVdW(blocksize, blocksize, ZERO);

  //-----------------------------------------------------------------
  // Main computation
  //-----------------------------------------------------------------

  // determine orientation
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

  //compute jacobian terms
  dFvdWf_Diamond(dFvdWf,dGvdWf,*SolnBlk, Orient_face, Rii, Rjj);
  dWfdWc_Diamond(dWfdWx,dWfdWy,*SolnBlk, Orient_face, Rii, Rjj, 
		 Orient_cell, lface*nface.x, lface*nface.y ); 
  
  // build the jacobian
  // dGVdW = lface * (nface.x*(dFvdWf*dWfdWx) + nface.y*(dGvdWf*dWfdWy));
  // NOTE: lface*nface.x and lface*nface.y already multiplied in
  if ( fabs(nface.x)>TOLER && fabs(nface.y)<=TOLER ) {
    dGVdW = dFvdWf*dWfdWx;
  } else if ( fabs(nface.x)<=TOLER  && fabs(nface.y)>TOLER)  {
    dGVdW = dGvdWf*dWfdWy;
  } else {
    dGVdW = dFvdWf*dWfdWx;
    dGVdW += dGvdWf*dWfdWy;
  }
  
  //transformation Jacobian
  // Non neeed to zero more than once, always writing to the same spot
  static DenseMatrix dWdU(blocksize,blocksize,ZERO);

  //transformation Jacobian  Wo == W here 
  SolnBlk->W[Wii][Wjj].dWdU(dWdU);  
  dFvdU += dGVdW*dWdU;

}


/*!**************************************************************
 * Calculate Second Order Local Jacobian Block(s) Coresponding  *
 * to Cell(i,j)                                                 *
 ****************************************************************/
template<> inline void Block_Preconditioner<Flame2D_pState,
					    Flame2D_Quad_Block,
					    Flame2D_Input_Parameters>::
Second_Order_Viscous_Jacobian(const int &cell_index_i,
			      const int &cell_index_j, 
			      DenseMatrix* Jacobian){

  // A real cludge with all the DenseMatrices and recalculations, 
  //  but just to test, need to change for performance....
  static DenseMatrix JacobianN(blocksize,blocksize);
  static DenseMatrix JacobianS(blocksize,blocksize);
  static DenseMatrix JacobianE(blocksize,blocksize);
  static DenseMatrix JacobianW(blocksize,blocksize);

  //Also should rewrite to minimize dR/dU calls and just 
  //call dWdU, but this needs to be done in Preconditioner_dFVdU
  
  /***************** dR(i,j)/dU(i,j) ******************************/
  //CENTER
  JacobianN.zero(); JacobianS.zero(); JacobianE.zero(); JacobianW.zero();
  Preconditioner_dFVdU(JacobianN,cell_index_i,cell_index_j,
		       cell_index_i,cell_index_j,NORTH,CENTER);
  Preconditioner_dFVdU(JacobianS,cell_index_i,cell_index_j,
		       cell_index_i,cell_index_j,SOUTH,CENTER);
  Preconditioner_dFVdU(JacobianE,cell_index_i,cell_index_j,
		       cell_index_i,cell_index_j,EAST,CENTER);
  Preconditioner_dFVdU(JacobianW,cell_index_i,cell_index_j,
		       cell_index_i,cell_index_j,WEST,CENTER);

  for (int i=0; i<blocksize; i++)
    for (int j=0; j<blocksize; j++) {
      Jacobian[CENTER](i,j) += ( (JacobianN(i,j) + JacobianS(i,j) + 
				  JacobianE(i,j)  + JacobianW(i,j)) /
				 SolnBlk->Grid.Cell[cell_index_i][cell_index_j].A );
    }				 

  /***************** dR(i,j-1)/dU(i,j) ****************************/
  //NORTH                           
  JacobianN.zero();
  Preconditioner_dFVdU(JacobianN,cell_index_i,cell_index_j-1,
		       cell_index_i,cell_index_j,EAST,NORTH); 
  Preconditioner_dFVdU(JacobianN,cell_index_i,cell_index_j-1,
		       cell_index_i,cell_index_j,NORTH,NORTH); 
  Preconditioner_dFVdU(JacobianN,cell_index_i,cell_index_j-1,
		       cell_index_i,cell_index_j,WEST,NORTH); 
 
  /***************** dR(i,j+1)/dU(i,j) ****************************/
  //SOUTH 
  JacobianS.zero();
  Preconditioner_dFVdU(JacobianS,cell_index_i,cell_index_j+1,
		       cell_index_i,cell_index_j,EAST,SOUTH);
  Preconditioner_dFVdU(JacobianS,cell_index_i,cell_index_j+1,
		       cell_index_i,cell_index_j,SOUTH,SOUTH);
  Preconditioner_dFVdU(JacobianS,cell_index_i,cell_index_j+1,
		       cell_index_i,cell_index_j,WEST,SOUTH);
  
  /***************** dR(i-1,j)/dU(i,j) ****************************/
  //EAST 
  JacobianE.zero();
  Preconditioner_dFVdU(JacobianE,cell_index_i-1,cell_index_j,
		       cell_index_i,cell_index_j,NORTH,EAST);
  Preconditioner_dFVdU(JacobianE,cell_index_i-1,cell_index_j,
		       cell_index_i,cell_index_j,EAST,EAST);
  Preconditioner_dFVdU(JacobianE,cell_index_i-1,cell_index_j,
		       cell_index_i,cell_index_j,SOUTH,EAST);

  /***************** dR(i+1,j)/dU(i,j) ****************************/
  //WEST
  JacobianW.zero(); 
  Preconditioner_dFVdU(JacobianW,cell_index_i+1,cell_index_j,
		       cell_index_i,cell_index_j, NORTH,WEST);
  Preconditioner_dFVdU(JacobianW,cell_index_i+1,cell_index_j,
		       cell_index_i,cell_index_j, WEST,WEST);
  Preconditioner_dFVdU(JacobianW,cell_index_i+1,cell_index_j,
		       cell_index_i,cell_index_j, SOUTH,WEST);

  /****************************************************************/
  // CORNER matrices are still zero at this point

  /***************** dR(i+1,j+1)/dU(i,j) **************************/
  //SOUTHWEST
  Preconditioner_dFVdU(Jacobian[SOUTH_WEST],cell_index_i+1,cell_index_j+1,
		       cell_index_i,cell_index_j,SOUTH,SOUTH_WEST);
  Preconditioner_dFVdU(Jacobian[SOUTH_WEST],cell_index_i+1,cell_index_j+1,
		       cell_index_i,cell_index_j, WEST,SOUTH_WEST);
  
  /***************** dR(i-1,j+1)/dU(i,j) **************************/
  //SOUTHEAST
  Preconditioner_dFVdU(Jacobian[SOUTH_EAST],cell_index_i-1,cell_index_j+1,
		       cell_index_i,cell_index_j,SOUTH,SOUTH_EAST);
  Preconditioner_dFVdU(Jacobian[SOUTH_EAST],cell_index_i-1,cell_index_j+1,
		       cell_index_i,cell_index_j,EAST,SOUTH_EAST);

  /***************** dR(i+1,j-1)/dU(i,j) **************************/
  //NORTHWEST
  Preconditioner_dFVdU(Jacobian[NORTH_WEST],cell_index_i+1,cell_index_j-1,
		       cell_index_i,cell_index_j, NORTH,NORTH_WEST);
  Preconditioner_dFVdU(Jacobian[NORTH_WEST],cell_index_i+1,cell_index_j-1,
		       cell_index_i,cell_index_j, WEST,NORTH_WEST);

  /***************** dR(i-1,j-1)/dU(i,j) **************************/
  //NORTHEAST
  Preconditioner_dFVdU(Jacobian[NORTH_EAST],cell_index_i-1,cell_index_j-1,
		       cell_index_i,cell_index_j, NORTH,NORTH_EAST);
  Preconditioner_dFVdU(Jacobian[NORTH_EAST],cell_index_i-1,cell_index_j-1,
		       cell_index_i,cell_index_j, EAST,NORTH_EAST);  
  /****************************************************************/
  for (int i=0; i<blocksize; i++)
    for (int j=0; j<blocksize; j++) {
      Jacobian[NORTH](i,j) += 
	JacobianN(i,j)/SolnBlk->Grid.Cell[cell_index_i][cell_index_j-1].A;
      Jacobian[SOUTH](i,j) += 
	JacobianS(i,j)/SolnBlk->Grid.Cell[cell_index_i][cell_index_j+1].A;
      Jacobian[EAST](i,j) += 
	JacobianE(i,j)/SolnBlk->Grid.Cell[cell_index_i-1][cell_index_j].A;
      Jacobian[WEST](i,j) += 
	JacobianW(i,j)/SolnBlk->Grid.Cell[cell_index_i+1][cell_index_j].A; 

      Jacobian[SOUTH_WEST](i,j) /= 
	SolnBlk->Grid.Cell[cell_index_i+1][cell_index_j+1].A;
      Jacobian[SOUTH_EAST](i,j) /= 
	SolnBlk->Grid.Cell[cell_index_i-1][cell_index_j+1].A;
      Jacobian[NORTH_WEST](i,j) /= 
	SolnBlk->Grid.Cell[cell_index_i+1][cell_index_j-1].A;
      Jacobian[NORTH_EAST](i,j) /= 
	SolnBlk->Grid.Cell[cell_index_i-1][cell_index_j-1].A;
    }
  /****************************************************************/

}



/*!**************************************************************
 * Specialization of Block_Preconditioner::Preconditioner_dSdU  *
 *                                                              *
 * Calculates the dFdU matrix used to generate the approximate  *               
 * Jacobian for the Block Preconditioner.                       *
 ****************************************************************/
template<> inline void Block_Preconditioner<Flame2D_pState,
					    Flame2D_Quad_Block,					    
					    Flame2D_Input_Parameters>::
Preconditioner_dSdU(int ii, int jj, DenseMatrix &dRdU){

  // build dRdU
  SemiImplicitBlockJacobi(dRdU,*SolnBlk,SolnBlk->W[ii][jj],ii,jj);

}


/*!**************************************************************
 *  Specialization of Block_Preconditioner::                    *
 *                               normalize_Preconditioner_dFdU  *
 *                                                              *
 * Normaliazes the dFdU matrix used to generate the approximate *               
 * Jacobian for the Block Preconditioner.                       *
 ****************************************************************/
void normalize_Preconditioner(DenseMatrix &dFdU) 
{ 
  int nsp( Flame2D_pState::NumSpecies() );
  static const Flame2D_pState W_STD_ATM;
  static double ao( W_STD_ATM.a() );
  static double rho( W_STD_ATM.rho() );

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

  //k,omega, and cs's all have same normalization.
  for(int i=NUM_FLAME2D_VAR_SANS_SPECIES; 
      i< NUM_FLAME2D_VAR_SANS_SPECIES + nsp; i++){   		  
    dFdU(0,i) *= (ONE/ao);
    dFdU(1,i) *= (ONE/(ao*ao));
    dFdU(2,i) *= (ONE/(ao*ao));
    dFdU(3,i) *= (ONE/(ao*ao*ao));
    dFdU(i,0) *= (ONE/ao);
    dFdU(i,3) *= ao;
    for(int j=NUM_FLAME2D_VAR_SANS_SPECIES; 
	j< NUM_FLAME2D_VAR_SANS_SPECIES + nsp; j++){   
      dFdU(i,j) *= (ONE/ao);            
    }
  } 

}

template<> inline void Block_Preconditioner<Flame2D_pState,
					    Flame2D_Quad_Block,
					    Flame2D_Input_Parameters>::
normalize_Preconditioner_dFdU(DenseMatrix &dFdU) 
{ 
  normalize_Preconditioner(dFdU); 
}

/*!**************************************************************
 *  Specialization of Block_Preconditioner::                    *
 *                                     Pre_Precon_SolnBlk_Init  *
 *                                                              *
 * Update the Wnd array of nodal values before we build a       *               
 * Second_Order_Viscous_Jacobian.                               *
 ****************************************************************/
template<> 
inline void Block_Preconditioner<Flame2D_pState,
				 Flame2D_Quad_Block,
				 Flame2D_Input_Parameters>::
Pre_Precon_SolnBlk_Init(void) {
  if ( Input_Parameters->NKS_IP.Jacobian_Order == SECOND_ORDER_DIAMOND_WITH_HLLE ||
       Input_Parameters->NKS_IP.Jacobian_Order == SECOND_ORDER_DIAMOND_WITH_ROE ||
       Input_Parameters->NKS_IP.Jacobian_Order == SECOND_ORDER_DIAMOND_WITH_AUSMPLUSUP ) {
    SolnBlk->Update_Nodal_Values();
  }
}


/************************************************************************/
/************ GMRES REQUIRED SPECIALIZATIONS & FUNCTIONS ****************/
/************************************************************************/ 
/**************************************************************************
 * GMRES_Block::SubcellReconstruction --                                  *
 *              Performs the subcell reconstruction of solution state     *
 *              within a given cell (i,j) of the computational mesh for   *
 *              the specified quadrilateral solution block.               *
 **************************************************************************/
template <>inline void GMRES_Block<Flame2D_pState,
				   Flame2D_Quad_Block,
				   Flame2D_Input_Parameters>::
SubcellReconstruction(const int i, 
		      const int j,
		      const int Limiter) {

  int n, n2, n_pts, i_index[8], j_index[8], k;
  double u0, u0Min, u0Max, uQuad[4], phi_n;
  double DxDx_ave, DxDy_ave, DyDy_ave, DU;
  static Vector2D dX, dXe, dXw, dXn, dXs;
  static Flame2D_State U0, DUDx_ave, DUDy_ave;

  /* Carry out the limited solution reconstruction in
     each cell of the computational mesh. */

  //FOR VISCOUS -> CHANGED TO USE ALL 8
  if (i == ICl-Nghost || i == ICu+Nghost ||
      j == JCl-Nghost || j == JCu+Nghost) {
    n_pts = 0;
  } else {
    n_pts = 8;
    i_index[0] = i-1; j_index[0] = j-1;
    i_index[1] = i  ; j_index[1] = j-1;
    i_index[2] = i+1; j_index[2] = j-1;
    i_index[3] = i-1; j_index[3] = j  ;
    i_index[4] = i+1; j_index[4] = j  ;
    i_index[5] = i-1; j_index[5] = j+1;
    i_index[6] = i  ; j_index[6] = j+1;
    i_index[7] = i+1; j_index[7] = j+1;
  } /* endif */

  if (n_pts > 0) {
    DUDx_ave.Vacuum();
    DUDy_ave.Vacuum();
    DxDx_ave = ZERO;
    DxDy_ave = ZERO;
    DyDy_ave = ZERO;

    for ( k = 0 ; k < blocksize; ++ k) {
      if (vector_switch) {
	U0[k+1] = W[(search_directions)*scalar_dim + index(i,j,k)];
      } else {
	U0[k+1] = x[index(i,j,k)];
      }
    } 

    for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
      dX = SolnBlk->Grid.Cell[ i_index[n2] ][ j_index[n2] ].Xc;
      dX -= SolnBlk->Grid.Cell[i][j].Xc;
      DxDx_ave += dX.x*dX.x;
      DxDy_ave += dX.x*dX.y;
      DyDy_ave += dX.y*dX.y;
      for ( k = 0 ; k < blocksize; ++ k) {
	if (vector_switch) {
	  DU = W[(search_directions)*scalar_dim + index(i_index[n2] , j_index[n2] , k)] -  U0[k+1];
	} else {
	  DU = x[index( i_index[n2] , j_index[n2] , k)] -  U0[k+1];
	} /* endif */
	DUDx_ave[k+1] += DU*dX.x;
	DUDy_ave[k+1] += DU*dX.y;
      } /* endfor */

    } /* endfor */
  					    
      // don't need to do this, it will cancel out
      // DUDx_ave = DUDx_ave/double(n_pts);
      // DUDy_ave = DUDy_ave/double(n_pts);
      // DxDx_ave = DxDx_ave/double(n_pts);
      // DxDy_ave = DxDy_ave/double(n_pts);
      // DyDy_ave = DyDy_ave/double(n_pts);

    for (k=1; k<=blocksize; k++) {
      SolnBlk->dWdx[i][j][k] = ( (DUDx_ave[k]*DyDy_ave-DUDy_ave[k]*DxDy_ave)/
				 (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave) );
      SolnBlk->dWdy[i][j][k] = ( (DUDy_ave[k]*DxDx_ave-DUDx_ave[k]*DxDy_ave)/
				 (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave) );
    }

    // Calculate slope limiters. 
    if (!SolnBlk->Freeze_Limiter) {

      dXe = SolnBlk->Grid.xfaceE(i, j)-SolnBlk->Grid.Cell[i][j].Xc;
      dXw = SolnBlk->Grid.xfaceW(i, j)-SolnBlk->Grid.Cell[i][j].Xc;
      dXn = SolnBlk->Grid.xfaceN(i, j)-SolnBlk->Grid.Cell[i][j].Xc;
      dXs = SolnBlk->Grid.xfaceS(i, j)-SolnBlk->Grid.Cell[i][j].Xc;

      for ( n = 1 ; n <= blocksize ; ++n ) {
	u0 = U0[n];
	u0Min = U0[n];
	u0Max = u0Min;
	for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
	  if (vector_switch) {
	    u0Min = min(u0Min, W[(search_directions)*scalar_dim + index(i_index[n2] , j_index[n2] ,n-1)]);
	    u0Max = max(u0Max, W[(search_directions)*scalar_dim + index(i_index[n2] , j_index[n2] ,n-1)]);
	  } else {
	    u0Min = min(u0Min, x[index(i_index[n2] , j_index[n2] , n-1)]);
	    u0Max = max(u0Max, x[index(i_index[n2] , j_index[n2] , n-1)]);
	  } /* endif */
	} /* endfor */
    
	uQuad[0] = u0 + 
	  SolnBlk->dWdx[i][j][n]*dXe.x +
	  SolnBlk->dWdy[i][j][n]*dXe.y ;
	uQuad[1] = u0 + 
	  SolnBlk->dWdx[i][j][n]*dXw.x +
	  SolnBlk->dWdy[i][j][n]*dXw.y ;
	uQuad[2] = u0 + 
	  SolnBlk->dWdx[i][j][n]*dXn.x +
	  SolnBlk->dWdy[i][j][n]*dXn.y ;
	uQuad[3] = u0 + 
	  SolnBlk->dWdx[i][j][n]*dXs.x +
	  SolnBlk->dWdy[i][j][n]*dXs.y ;
    
	switch(Limiter) {
	case LIMITER_ONE :
	  phi_n = ONE;
	  break;
	case LIMITER_ZERO :
	  phi_n = ZERO;
	  break;
	case LIMITER_BARTH_JESPERSEN :
	  phi_n = Limiter_BarthJespersen(uQuad, u0, 
					 u0Min, u0Max, 4);
	  break;
	case LIMITER_VENKATAKRISHNAN :
	  phi_n = Limiter_Venkatakrishnan(uQuad, u0, 
					  u0Min, u0Max, 4);
	  break;
	case LIMITER_VANLEER :
	  phi_n = Limiter_VanLeer(uQuad, u0, 
				  u0Min, u0Max, 4);
	  break;
	case LIMITER_VANALBADA :
	  phi_n = Limiter_VanAlbada(uQuad, u0, 
				    u0Min, u0Max, 4);
	  break;
	default:
	  phi_n = Limiter_BarthJespersen(uQuad, u0, 
					 u0Min, u0Max, 4);
	  break;
	} /* endswitch */

	SolnBlk->phi[i][j][n] = phi_n;

      } /* endfor */
    } /* endif */
  } else {
    SolnBlk->dWdx[i][j].Vacuum();
    SolnBlk->dWdy[i][j].Vacuum();
    SolnBlk->phi[i][j].Vacuum();
  } /* endif */


}
  
/**************************************************************************
 * GMRES_Block::LoadSendBuffer_C2F -- Loads send message buffer for       *
 *                                    coarse to fine block message        *
 *                                    passing.                            *
 **************************************************************************/
template <> 
inline int GMRES_Block<Flame2D_pState,
		       Flame2D_Quad_Block,
		       Flame2D_Input_Parameters>::
LoadSendBuffer_C2F(double *buffer,
		   int &buffer_count,
		   const int buffer_size,
		   const int i_min, 
		   const int i_max,
		   const int i_inc,
		   const int j_min, 
		   const int j_max,
		   const int j_inc,
		   const int face,
		   const int sector) {
  int i, j, k;
  static Vector2D dX;
  static Flame2D_State Wcoarse, Wfine;
  int LIMITER = LIMITER_ZERO; //LIMITER_VENKATAKRISHNAN 

  if (j_min == j_max) { // North or south boundary.
    // Four different orderings to consider depending on the value of i_inc & j_inc.
    if (j_inc > 0) {             
      if (i_inc > 0) {

	/************************ CASE #1 ******************************/
	for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	  // Perform limited linear least squares reconstruction in cell (i, j_min).
	  SubcellReconstruction(i, j_min, LIMITER);                    
	  
	  // Evaluate SW sub (fine) cell values.
	  for (k = 0 ; k < blocksize; ++ k) {
	    if (vector_switch) {
	      Wcoarse[k+1] = W[(search_directions)*scalar_dim + index(i,j_min,k)];
	    } else {
	      Wcoarse[k+1] = x[index(i,j_min,k)];
	    } 
	  } 
	  dX = (SolnBlk->Grid.Node[i][j_min].X+
		HALF*(SolnBlk->Grid.Node[i][j_min].X+SolnBlk->Grid.Node[i+1][j_min].X)+
		HALF*(SolnBlk->Grid.Node[i][j_min].X+SolnBlk->Grid.Node[i][j_min+1].X)+
		SolnBlk->Grid.Cell[i][j_min].Xc)/FOUR -
	    SolnBlk->Grid.Cell[i][j_min].Xc;
	  Wfine = Wcoarse +
	    (SolnBlk->phi[i][j_min]^SolnBlk->dWdx[i][j_min])*dX.x +
	    (SolnBlk->phi[i][j_min]^SolnBlk->dWdy[i][j_min])*dX.y;

	  for (k = 0 ; k < blocksize; ++ k) {
	    buffer_count = buffer_count + 1;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Wfine[k+1];
	  } 
	  
	  // Evaluate SE sub (fine) cell values.
	  dX = (HALF*(SolnBlk->Grid.Node[i][j_min].X+SolnBlk->Grid.Node[i+1][j_min].X)+
		SolnBlk->Grid.Node[i+1][j_min].X + SolnBlk->Grid.Cell[i][j_min].Xc+
		HALF*(SolnBlk->Grid.Node[i+1][j_min].X+SolnBlk->Grid.Node[i+1][j_min+1].X))/FOUR -
	    SolnBlk->Grid.Cell[i][j_min].Xc;
	  Wfine = Wcoarse +
	    (SolnBlk->phi[i][j_min]^SolnBlk->dWdx[i][j_min])*dX.x +
	    (SolnBlk->phi[i][j_min]^SolnBlk->dWdy[i][j_min])*dX.y;
            
	  for ( k = 0 ; k < blocksize; ++ k) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Wfine[k+1];
	  }
	} 
	
	for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	  // Evaluate NW sub (fine) cell values.
	  for (k = 0 ; k < blocksize; ++ k) {
	    if (vector_switch) {
	      Wcoarse[k+1] = W[(search_directions)*scalar_dim + index(i,j_min,k)];
	    } else {
	      Wcoarse[k+1] = x[index(i,j_min,k)];
	    } /* endif */
	  } /* endfor */
	  dX = (HALF*(SolnBlk->Grid.Node[i][j_min].X+SolnBlk->Grid.Node[i][j_min+1].X)+
		SolnBlk->Grid.Cell[i][j_min].Xc+
		SolnBlk->Grid.Node[i][j_min+1].X+
		HALF*(SolnBlk->Grid.Node[i][j_min+1].X+SolnBlk->Grid.Node[i+1][j_min+1].X))/FOUR -
	    SolnBlk->Grid.Cell[i][j_min].Xc;
	  Wfine = Wcoarse +
	    (SolnBlk->phi[i][j_min]^SolnBlk->dWdx[i][j_min])*dX.x +
	    (SolnBlk->phi[i][j_min]^SolnBlk->dWdy[i][j_min])*dX.y;
          
	  for (k = 0 ; k < blocksize; ++ k) {
	    buffer_count = buffer_count + 1;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Wfine[k+1];
	  } /* endfor */
	  // Evaluate NE sub (fine) cell values.
	  dX = (SolnBlk->Grid.Cell[i][j_min].Xc+
		HALF*(SolnBlk->Grid.Node[i+1][j_min].X+SolnBlk->Grid.Node[i+1][j_min+1].X)+
		HALF*(SolnBlk->Grid.Node[i][j_min+1].X+SolnBlk->Grid.Node[i+1][j_min+1].X)+
		SolnBlk->Grid.Node[i+1][j_min+1].X)/FOUR -
	    SolnBlk->Grid.Cell[i][j_min].Xc;
	  Wfine = Wcoarse +
	    (SolnBlk->phi[i][j_min]^SolnBlk->dWdx[i][j_min])*dX.x +
	    (SolnBlk->phi[i][j_min]^SolnBlk->dWdy[i][j_min])*dX.y;
	  for (k = 0 ; k < blocksize; ++ k) {
	    buffer_count = buffer_count + 1;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Wfine[k+1];
	  } /* endfor */
	} /* endfor */
	
	/************************ CASE #2 ******************************/
      } else {
	for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	  // Perform limited linear least squares reconstruction in cell (i, j_min).
	  SubcellReconstruction(i, j_min, LIMITER);
	  // Evaluate SE sub (fine) cell values.
	  for (k = 0 ; k < blocksize; ++ k) {
	    if (vector_switch) {
	      Wcoarse[k+1] = W[(search_directions)*scalar_dim + index(i,j_min,k)];
	    } else {
	      Wcoarse[k+1] = x[index(i,j_min,k)];
	    } /* endif */
	  } /* endfor */
	  dX = (HALF*(SolnBlk->Grid.Node[i][j_min].X+SolnBlk->Grid.Node[i+1][j_min].X)+
		SolnBlk->Grid.Node[i+1][j_min].X+
		SolnBlk->Grid.Cell[i][j_min].Xc+
		HALF*(SolnBlk->Grid.Node[i+1][j_min].X+SolnBlk->Grid.Node[i+1][j_min+1].X))/FOUR -
	    SolnBlk->Grid.Cell[i][j_min].Xc;
	  Wfine = Wcoarse +
	    (SolnBlk->phi[i][j_min]^SolnBlk->dWdx[i][j_min])*dX.x +
	    (SolnBlk->phi[i][j_min]^SolnBlk->dWdy[i][j_min])*dX.y;
	  for (k = 0 ; k < blocksize; ++ k) {
	    buffer_count = buffer_count + 1;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Wfine[k+1];
	  } /* endfor */
	  // Evaluate SW sub (fine) cell values.
	  dX = (SolnBlk->Grid.Node[i][j_min].X+
		HALF*(SolnBlk->Grid.Node[i][j_min].X+SolnBlk->Grid.Node[i+1][j_min].X)+
		HALF*(SolnBlk->Grid.Node[i][j_min].X+SolnBlk->Grid.Node[i][j_min+1].X)+
		SolnBlk->Grid.Cell[i][j_min].Xc)/FOUR -
	    SolnBlk->Grid.Cell[i][j_min].Xc;
	  Wfine = Wcoarse +
	    (SolnBlk->phi[i][j_min]^SolnBlk->dWdx[i][j_min])*dX.x +
	    (SolnBlk->phi[i][j_min]^SolnBlk->dWdy[i][j_min])*dX.y;
	  for ( k = 0 ; k < blocksize; ++ k) {
	    buffer_count = buffer_count + 1;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Wfine[k+1];
	  } /* endfor */
	} /* endfor */
	for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	  // Evaluate NE sub (fine) cell values.
	  for (k = 0 ; k < blocksize; ++ k) {
	    if (vector_switch) {
	      Wcoarse[k+1] = W[(search_directions)*scalar_dim + index(i,j_min,k)];
	    } else {
	      Wcoarse[k+1] = x[index(i,j_min,k)];
	    } /* endif */
	  } /* endfor */
	  dX = (SolnBlk->Grid.Cell[i][j_min].Xc+
		HALF*(SolnBlk->Grid.Node[i+1][j_min].X+SolnBlk->Grid.Node[i+1][j_min+1].X)+
		HALF*(SolnBlk->Grid.Node[i][j_min+1].X+SolnBlk->Grid.Node[i+1][j_min+1].X)+
		SolnBlk->Grid.Node[i+1][j_min+1].X)/FOUR -
	    SolnBlk->Grid.Cell[i][j_min].Xc;
	  Wfine = Wcoarse +
	    (SolnBlk->phi[i][j_min]^SolnBlk->dWdx[i][j_min])*dX.x +
	    (SolnBlk->phi[i][j_min]^SolnBlk->dWdy[i][j_min])*dX.y;
	  for (k = 0 ; k < blocksize; ++ k) {
	    buffer_count = buffer_count + 1;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Wfine[k+1];
	  } /* endfor */
	  // Evaluate NW sub (fine) cell values.
	  dX = (HALF*(SolnBlk->Grid.Node[i][j_min].X+SolnBlk->Grid.Node[i][j_min+1].X)+
		SolnBlk->Grid.Cell[i][j_min].Xc+
		SolnBlk->Grid.Node[i][j_min+1].X+
		HALF*(SolnBlk->Grid.Node[i][j_min+1].X+SolnBlk->Grid.Node[i+1][j_min+1].X))/FOUR -
	    SolnBlk->Grid.Cell[i][j_min].Xc;
	  Wfine = Wcoarse +
	    (SolnBlk->phi[i][j_min]^SolnBlk->dWdx[i][j_min])*dX.x +
	    (SolnBlk->phi[i][j_min]^SolnBlk->dWdy[i][j_min])*dX.y;
	  for (k = 0 ; k < blocksize; ++ k) {
	    buffer_count = buffer_count + 1;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Wfine[k+1];
	  } /* endfor */
	} /* endfor */
      } /* endif */

	/************************ CASE #3 ******************************/
    } else {
      if (i_inc > 0) {
	for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	  // Perform limited linear least squares reconstruction in cell (i, j_min).
	  SubcellReconstruction(i, j_min, LIMITER);
	  // Evaluate NW sub (fine) cell values.
	  for (k = 0 ; k < blocksize; ++ k) {
	    if (vector_switch) {
	      Wcoarse[k+1] = W[(search_directions)*scalar_dim + index(i,j_min,k)];
	    } else {
	      Wcoarse[k+1] = x[index(i,j_min,k)];
	    } /* endif */
	  } /* endfor */
	  dX = (HALF*(SolnBlk->Grid.Node[i][j_min].X+SolnBlk->Grid.Node[i][j_min+1].X)+
		SolnBlk->Grid.Cell[i][j_min].Xc+
		SolnBlk->Grid.Node[i][j_min+1].X+
		HALF*(SolnBlk->Grid.Node[i][j_min+1].X+SolnBlk->Grid.Node[i+1][j_min+1].X))/FOUR -
	    SolnBlk->Grid.Cell[i][j_min].Xc;
	  Wfine = Wcoarse +
	    (SolnBlk->phi[i][j_min]^SolnBlk->dWdx[i][j_min])*dX.x +
	    (SolnBlk->phi[i][j_min]^SolnBlk->dWdy[i][j_min])*dX.y;
	  for (k = 0 ; k < blocksize; ++ k) {
	    buffer_count = buffer_count + 1;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Wfine[k+1];
	  } /* endfor */
	  // Evaluate NE sub (fine) cell values.
	  dX = (SolnBlk->Grid.Cell[i][j_min].Xc+
		HALF*(SolnBlk->Grid.Node[i+1][j_min].X+SolnBlk->Grid.Node[i+1][j_min+1].X)+
		HALF*(SolnBlk->Grid.Node[i][j_min+1].X+SolnBlk->Grid.Node[i+1][j_min+1].X)+
		SolnBlk->Grid.Node[i+1][j_min+1].X)/FOUR -
	    SolnBlk->Grid.Cell[i][j_min].Xc;
	  Wfine = Wcoarse +
	    (SolnBlk->phi[i][j_min]^SolnBlk->dWdx[i][j_min])*dX.x +
	    (SolnBlk->phi[i][j_min]^SolnBlk->dWdy[i][j_min])*dX.y;
	  for (k = 0 ; k < blocksize; ++ k) {
	    buffer_count = buffer_count + 1;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Wfine[k+1];
	  } /* endfor */
	} /* endfor */
	for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	  // Evaluate SW sub (fine) cell values.
	  for (k = 0 ; k < blocksize; ++ k) {
	    if (vector_switch) {
	      Wcoarse[k+1] = W[(search_directions)*scalar_dim + index(i,j_min,k)];
	    } else {
	      Wcoarse[k+1] = x[index(i,j_min,k)];
	    } /* endif */
	  } /* endfor */
	  dX = (SolnBlk->Grid.Node[i][j_min].X+
		HALF*(SolnBlk->Grid.Node[i][j_min].X+SolnBlk->Grid.Node[i+1][j_min].X)+
		HALF*(SolnBlk->Grid.Node[i][j_min].X+SolnBlk->Grid.Node[i][j_min+1].X)+
		SolnBlk->Grid.Cell[i][j_min].Xc)/FOUR -
	    SolnBlk->Grid.Cell[i][j_min].Xc;
	  Wfine = Wcoarse +
	    (SolnBlk->phi[i][j_min]^SolnBlk->dWdx[i][j_min])*dX.x +
	    (SolnBlk->phi[i][j_min]^SolnBlk->dWdy[i][j_min])*dX.y;
	  for ( k = 0 ; k < blocksize; ++ k) {
	    buffer_count = buffer_count + 1;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Wfine[k+1];
	  } /* endfor */
	  // Evaluate SE sub (fine) cell values.
	  dX = (HALF*(SolnBlk->Grid.Node[i][j_min].X+SolnBlk->Grid.Node[i+1][j_min].X)+
		SolnBlk->Grid.Node[i+1][j_min].X+
		SolnBlk->Grid.Cell[i][j_min].Xc+
		HALF*(SolnBlk->Grid.Node[i+1][j_min].X+SolnBlk->Grid.Node[i+1][j_min+1].X))/FOUR -
	    SolnBlk->Grid.Cell[i][j_min].Xc;
	  Wfine = Wcoarse +
	    (SolnBlk->phi[i][j_min]^SolnBlk->dWdx[i][j_min])*dX.x +
	    (SolnBlk->phi[i][j_min]^SolnBlk->dWdy[i][j_min])*dX.y;
	  for (k = 0 ; k < blocksize; ++ k) {
	    buffer_count = buffer_count + 1;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Wfine[k+1];
	  } /* endfor */
	} /* endfor */
	/************************ CASE #4 ******************************/
      } else {
	for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	  // Perform limited linear least squares reconstruction in cell (i, j_min).
	  SubcellReconstruction(i, j_min, LIMITER);
	  // Evaluate NE sub (fine) cell values.
	  for (k = 0 ; k < blocksize; ++ k) {
	    if (vector_switch) {
	      Wcoarse[k+1] = W[(search_directions)*scalar_dim + index(i,j_min,k)];
	    } else {
	      Wcoarse[k+1] = x[index(i,j_min,k)];
	    } /* endif */
	  } /* endfor */
	  dX = (SolnBlk->Grid.Cell[i][j_min].Xc+
		HALF*(SolnBlk->Grid.Node[i+1][j_min].X+SolnBlk->Grid.Node[i+1][j_min+1].X)+
		HALF*(SolnBlk->Grid.Node[i][j_min+1].X+SolnBlk->Grid.Node[i+1][j_min+1].X)+
		SolnBlk->Grid.Node[i+1][j_min+1].X)/FOUR -
	    SolnBlk->Grid.Cell[i][j_min].Xc;
	  Wfine = Wcoarse +
	    (SolnBlk->phi[i][j_min]^SolnBlk->dWdx[i][j_min])*dX.x +
	    (SolnBlk->phi[i][j_min]^SolnBlk->dWdy[i][j_min])*dX.y;

	  // Wfine = Wfine.U();  ??MISTAKE MADE IN KALVINS ORIGINAL ?????

	  for (k = 0 ; k < blocksize; ++ k) {
	    buffer_count = buffer_count + 1;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Wfine[k+1];
	  } /* endfor */
	  // Evaluate NW sub (fine) cell values.
	  dX = (HALF*(SolnBlk->Grid.Node[i][j_min].X+SolnBlk->Grid.Node[i][j_min+1].X)+
		SolnBlk->Grid.Cell[i][j_min].Xc+
		SolnBlk->Grid.Node[i][j_min+1].X+
		HALF*(SolnBlk->Grid.Node[i][j_min+1].X+SolnBlk->Grid.Node[i+1][j_min+1].X))/FOUR -
	    SolnBlk->Grid.Cell[i][j_min].Xc;
	  Wfine = Wcoarse +
	    (SolnBlk->phi[i][j_min]^SolnBlk->dWdx[i][j_min])*dX.x +
	    (SolnBlk->phi[i][j_min]^SolnBlk->dWdy[i][j_min])*dX.y;
	  for (k = 0 ; k < blocksize; ++ k) {
	    buffer_count = buffer_count + 1;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Wfine[k+1];
	  } /* endfor */
	} /* endfor */
	for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	  // Evaluate SE sub (fine) cell values.
	  for (k = 0 ; k < blocksize; ++ k) {
	    if (vector_switch) {
	      Wcoarse[k+1] = W[(search_directions)*scalar_dim + index(i,j_min,k)];
	    } else {
	      Wcoarse[k+1] = x[index(i,j_min,k)];
	    } /* endif */
	  } /* endfor */
	  dX = (HALF*(SolnBlk->Grid.Node[i][j_min].X+SolnBlk->Grid.Node[i+1][j_min].X)+
		SolnBlk->Grid.Node[i+1][j_min].X+
		SolnBlk->Grid.Cell[i][j_min].Xc+
		HALF*(SolnBlk->Grid.Node[i+1][j_min].X+SolnBlk->Grid.Node[i+1][j_min+1].X))/FOUR -
	    SolnBlk->Grid.Cell[i][j_min].Xc;
	  Wfine = Wcoarse +
	    (SolnBlk->phi[i][j_min]^SolnBlk->dWdx[i][j_min])*dX.x +
	    (SolnBlk->phi[i][j_min]^SolnBlk->dWdy[i][j_min])*dX.y;
	  for ( k = 0 ; k < blocksize; ++ k) {
	    buffer_count = buffer_count + 1;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Wfine[k+1];
	  } /* endfor */
	  // Evaluate SW sub (fine) cell values.
	  dX = (SolnBlk->Grid.Node[i][j_min].X+
		HALF*(SolnBlk->Grid.Node[i][j_min].X+SolnBlk->Grid.Node[i+1][j_min].X)+
		HALF*(SolnBlk->Grid.Node[i][j_min].X+SolnBlk->Grid.Node[i][j_min+1].X)+
		SolnBlk->Grid.Cell[i][j_min].Xc)/FOUR -
	    SolnBlk->Grid.Cell[i][j_min].Xc;
	  Wfine = Wcoarse +
	    (SolnBlk->phi[i][j_min]^SolnBlk->dWdx[i][j_min])*dX.x +
	    (SolnBlk->phi[i][j_min]^SolnBlk->dWdy[i][j_min])*dX.y;
	  for (k = 0 ; k < blocksize; ++ k) {
	    buffer_count = buffer_count + 1;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Wfine[k+1];
	  } /* endfor */
	} /* endfor */
      } /* endif */
    } /* endif */


    /************************ CASE #5 ******************************/
  } else { // East or west boundary.
    // Four different orderings to consider depending on the value of i_inc & j_inc.
    if (j_inc > 0) {
      if (i_inc > 0) {
	for ( j = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
	  // Perform limited linear least squares reconstruction in cell (i_min, j).
	  SubcellReconstruction(i_min, j, LIMITER);
	  // Evaluate SW sub (fine) cell values.
	  for (k = 0 ; k < blocksize; ++ k) {
	    if (vector_switch) {
	      Wcoarse[k+1] = W[(search_directions)*scalar_dim + index(i_min,j,k)];
	    } else {
	      Wcoarse[k+1] = x[index(i_min,j,k)];
	    } /* endif */
	  } /* endfor */
	  dX = (SolnBlk->Grid.Node[i_min][j].X+
		HALF*(SolnBlk->Grid.Node[i_min][j].X+SolnBlk->Grid.Node[i_min+1][j].X)+
		HALF*(SolnBlk->Grid.Node[i_min][j].X+SolnBlk->Grid.Node[i_min][j+1].X)+
		SolnBlk->Grid.Cell[i_min][j].Xc)/FOUR -
	    SolnBlk->Grid.Cell[i_min][j].Xc;
	  Wfine = Wcoarse +
	    (SolnBlk->phi[i_min][j]^SolnBlk->dWdx[i_min][j])*dX.x +
	    (SolnBlk->phi[i_min][j]^SolnBlk->dWdy[i_min][j])*dX.y;
	  for (k = 0 ; k < blocksize; ++ k) {
	    buffer_count = buffer_count + 1;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Wfine[k+1];
	  } /* endfor */
	  // Evaluate SE sub (fine) cell values.
	  dX = (HALF*(SolnBlk->Grid.Node[i_min][j].X+SolnBlk->Grid.Node[i_min+1][j].X)+
		SolnBlk->Grid.Node[i_min+1][j].X+
		SolnBlk->Grid.Cell[i_min][j].Xc+
		HALF*(SolnBlk->Grid.Node[i_min+1][j].X+SolnBlk->Grid.Node[i_min+1][j+1].X))/FOUR -
	    SolnBlk->Grid.Cell[i_min][j].Xc;
	  Wfine = Wcoarse +
	    (SolnBlk->phi[i_min][j]^SolnBlk->dWdx[i_min][j])*dX.x +
	    (SolnBlk->phi[i_min][j]^SolnBlk->dWdy[i_min][j])*dX.y;
	  for (k = 0 ; k < blocksize; ++ k) {
	    buffer_count = buffer_count + 1;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Wfine[k+1];
	  } /* endfor */
	  // Evaluate NW sub (fine) cell values.
	  dX = (HALF*(SolnBlk->Grid.Node[i_min][j].X+SolnBlk->Grid.Node[i_min][j+1].X)+
		SolnBlk->Grid.Cell[i_min][j].Xc+
		SolnBlk->Grid.Node[i_min][j+1].X+
		HALF*(SolnBlk->Grid.Node[i_min][j+1].X+SolnBlk->Grid.Node[i_min+1][j+1].X))/FOUR -
	    SolnBlk->Grid.Cell[i_min][j].Xc;
	  Wfine = Wcoarse +
	    (SolnBlk->phi[i_min][j]^SolnBlk->dWdx[i_min][j])*dX.x +
	    (SolnBlk->phi[i_min][j]^SolnBlk->dWdy[i_min][j])*dX.y;
	  for (k = 0 ; k < blocksize; ++ k) {
	    buffer_count = buffer_count + 1;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Wfine[k+1];
	  } /* endfor */
	  // Evaluate NE sub (fine) cell values.
	  dX = (SolnBlk->Grid.Cell[i_min][j].Xc+
		HALF*(SolnBlk->Grid.Node[i_min+1][j].X+SolnBlk->Grid.Node[i_min+1][j+1].X)+
		HALF*(SolnBlk->Grid.Node[i_min][j+1].X+SolnBlk->Grid.Node[i_min+1][j+1].X)+
		SolnBlk->Grid.Node[i_min+1][j+1].X)/FOUR -
	    SolnBlk->Grid.Cell[i_min][j].Xc;
	  Wfine = Wcoarse +
	    (SolnBlk->phi[i_min][j]^SolnBlk->dWdx[i_min][j])*dX.x +
	    (SolnBlk->phi[i_min][j]^SolnBlk->dWdy[i_min][j])*dX.y;
	  for (k = 0 ; k < blocksize; ++ k) {
	    buffer_count = buffer_count + 1;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Wfine[k+1];
	  } /* endfor */
	} /* endfor */
	/************************ CASE #6 ******************************/
      } else {
	for ( j = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
	  // Perform limited linear least squares reconstruction in cell (i_min, j).
	  SubcellReconstruction(i_min, j, LIMITER);
	  // Evaluate SE sub (fine) cell values.
	  for (k = 0 ; k < blocksize; ++ k) {
	    if (vector_switch) {
	      Wcoarse[k+1] = W[(search_directions)*scalar_dim + index(i_min,j,k)];
	    } else {
	      Wcoarse[k+1] = x[index(i_min,j,k)];
	    } /* endif */
	  } /* endfor */
	  dX = (HALF*(SolnBlk->Grid.Node[i_min][j].X+SolnBlk->Grid.Node[i_min+1][j].X)+
		SolnBlk->Grid.Node[i_min+1][j].X+
		SolnBlk->Grid.Cell[i_min][j].Xc+
		HALF*(SolnBlk->Grid.Node[i_min+1][j].X+SolnBlk->Grid.Node[i_min+1][j+1].X))/FOUR -
	    SolnBlk->Grid.Cell[i_min][j].Xc;
	  Wfine = Wcoarse +
	    (SolnBlk->phi[i_min][j]^SolnBlk->dWdx[i_min][j])*dX.x +
	    (SolnBlk->phi[i_min][j]^SolnBlk->dWdy[i_min][j])*dX.y;
	  for (k = 0 ; k < blocksize; ++ k) {
	    buffer_count = buffer_count + 1;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Wfine[k+1];
	  } /* endfor */
	  // Evaluate SW sub (fine) cell values.
	  dX = (SolnBlk->Grid.Node[i_min][j].X+
		HALF*(SolnBlk->Grid.Node[i_min][j].X+SolnBlk->Grid.Node[i_min+1][j].X)+
		HALF*(SolnBlk->Grid.Node[i_min][j].X+SolnBlk->Grid.Node[i_min][j+1].X)+
		SolnBlk->Grid.Cell[i_min][j].Xc)/FOUR -
	    SolnBlk->Grid.Cell[i_min][j].Xc;
	  Wfine = Wcoarse +
	    (SolnBlk->phi[i_min][j]^SolnBlk->dWdx[i_min][j])*dX.x +
	    (SolnBlk->phi[i_min][j]^SolnBlk->dWdy[i_min][j])*dX.y;
	  for (k = 0 ; k < blocksize; ++ k) {
	    buffer_count = buffer_count + 1;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Wfine[k+1];
	  } /* endfor */
	  // Evaluate NE sub (fine) cell values.
	  dX = (SolnBlk->Grid.Cell[i_min][j].Xc+
		HALF*(SolnBlk->Grid.Node[i_min+1][j].X+SolnBlk->Grid.Node[i_min+1][j+1].X)+
		HALF*(SolnBlk->Grid.Node[i_min][j+1].X+SolnBlk->Grid.Node[i_min+1][j+1].X)+
		SolnBlk->Grid.Node[i_min+1][j+1].X)/FOUR -
	    SolnBlk->Grid.Cell[i_min][j].Xc;
	  Wfine = Wcoarse +
	    (SolnBlk->phi[i_min][j]^SolnBlk->dWdx[i_min][j])*dX.x +
	    (SolnBlk->phi[i_min][j]^SolnBlk->dWdy[i_min][j])*dX.y;
	  for (k = 0 ; k < blocksize; ++ k) {
	    buffer_count = buffer_count + 1;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Wfine[k+1];
	  } /* endfor */
	  // Evaluate NW sub (fine) cell values.
	  dX = (HALF*(SolnBlk->Grid.Node[i_min][j].X+SolnBlk->Grid.Node[i_min][j+1].X)+
		SolnBlk->Grid.Cell[i_min][j].Xc+
		SolnBlk->Grid.Node[i_min][j+1].X+
		HALF*(SolnBlk->Grid.Node[i_min][j+1].X+SolnBlk->Grid.Node[i_min+1][j+1].X))/FOUR -
	    SolnBlk->Grid.Cell[i_min][j].Xc;
	  Wfine = Wcoarse +
	    (SolnBlk->phi[i_min][j]^SolnBlk->dWdx[i_min][j])*dX.x +
	    (SolnBlk->phi[i_min][j]^SolnBlk->dWdy[i_min][j])*dX.y;
	  for (k = 0 ; k < blocksize; ++ k) {
	    buffer_count = buffer_count + 1;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Wfine[k+1];
	  } /* endfor */
	} /* endfor */
      } /* endif */	
	/************************ CASE #7 ******************************/
    } else {
      if (i_inc > 0) {
	for ( j = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
	  // Perform limited linear least squares reconstruction in cell (i_min, j).
	  SubcellReconstruction(i_min, j, LIMITER);
	  // Evaluate NW sub (fine) cell values.
	  for (k = 0 ; k < blocksize; ++ k) {
	    if (vector_switch) {
	      Wcoarse[k+1] = W[(search_directions)*scalar_dim + index(i_min,j,k)];
	    } else {
	      Wcoarse[k+1] = x[index(i_min,j,k)];
	    } /* endif */
	  } /* endfor */
	  dX = (HALF*(SolnBlk->Grid.Node[i_min][j].X+SolnBlk->Grid.Node[i_min][j+1].X)+
		SolnBlk->Grid.Cell[i_min][j].Xc+
		SolnBlk->Grid.Node[i_min][j+1].X+
		HALF*(SolnBlk->Grid.Node[i_min][j+1].X+SolnBlk->Grid.Node[i_min+1][j+1].X))/FOUR -
	    SolnBlk->Grid.Cell[i_min][j].Xc;
	  Wfine = Wcoarse +
	    (SolnBlk->phi[i_min][j]^SolnBlk->dWdx[i_min][j])*dX.x +
	    (SolnBlk->phi[i_min][j]^SolnBlk->dWdy[i_min][j])*dX.y;
	  for (k = 0 ; k < blocksize; ++ k) {
	    buffer_count = buffer_count + 1;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Wfine[k+1];
	  } /* endfor */
	  // Evaluate NE sub (fine) cell values.
	  dX = (SolnBlk->Grid.Cell[i_min][j].Xc+
		HALF*(SolnBlk->Grid.Node[i_min+1][j].X+SolnBlk->Grid.Node[i_min+1][j+1].X)+
		HALF*(SolnBlk->Grid.Node[i_min][j+1].X+SolnBlk->Grid.Node[i_min+1][j+1].X)+
		SolnBlk->Grid.Node[i_min+1][j+1].X)/FOUR -
	    SolnBlk->Grid.Cell[i_min][j].Xc;
	  Wfine = Wcoarse +
	    (SolnBlk->phi[i_min][j]^SolnBlk->dWdx[i_min][j])*dX.x +
	    (SolnBlk->phi[i_min][j]^SolnBlk->dWdy[i_min][j])*dX.y;
	  for (k = 0 ; k < blocksize; ++ k) {
	    buffer_count = buffer_count + 1;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Wfine[k+1];
	  } /* endfor */
	  // Evaluate SW sub (fine) cell values.
	  dX = (SolnBlk->Grid.Node[i_min][j].X+
		HALF*(SolnBlk->Grid.Node[i_min][j].X+SolnBlk->Grid.Node[i_min+1][j].X)+
		HALF*(SolnBlk->Grid.Node[i_min][j].X+SolnBlk->Grid.Node[i_min][j+1].X)+
		SolnBlk->Grid.Cell[i_min][j].Xc)/FOUR -
	    SolnBlk->Grid.Cell[i_min][j].Xc;
	  Wfine = Wcoarse +
	    (SolnBlk->phi[i_min][j]^SolnBlk->dWdx[i_min][j])*dX.x +
	    (SolnBlk->phi[i_min][j]^SolnBlk->dWdy[i_min][j])*dX.y;
	  for (k = 0 ; k < blocksize; ++ k) {
	    buffer_count = buffer_count + 1;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Wfine[k+1];
	  } /* endfor */
	  // Evaluate SE sub (fine) cell values.
	  dX = (HALF*(SolnBlk->Grid.Node[i_min][j].X+SolnBlk->Grid.Node[i_min+1][j].X)+
		SolnBlk->Grid.Node[i_min+1][j].X+
		SolnBlk->Grid.Cell[i_min][j].Xc+
		HALF*(SolnBlk->Grid.Node[i_min+1][j].X+SolnBlk->Grid.Node[i_min+1][j+1].X))/FOUR -
	    SolnBlk->Grid.Cell[i_min][j].Xc;
	  Wfine = Wcoarse +
	    (SolnBlk->phi[i_min][j]^SolnBlk->dWdx[i_min][j])*dX.x +
	    (SolnBlk->phi[i_min][j]^SolnBlk->dWdy[i_min][j])*dX.y;
	  for (k = 0 ; k < blocksize; ++ k) {
	    buffer_count = buffer_count + 1;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Wfine[k+1];
	  } /* endfor */
	} /* endfor */\

	/************************ CASE #8 ******************************/
      } else {
	for ( j = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
	  // Perform limited linear least squares reconstruction in cell (i_min, j).
	  SubcellReconstruction(i_min, j, LIMITER);
	  // Evaluate NE sub (fine) cell values.
	  for (k = 0 ; k < blocksize; ++ k) {
	    if (vector_switch) {
	      Wcoarse[k+1] = W[(search_directions)*scalar_dim + index(i_min,j,k)];
	    } else {
	      Wcoarse[k+1] = x[index(i_min,j,k)];
	    } /* endif */
	  } /* endfor */
	  dX = (SolnBlk->Grid.Cell[i_min][j].Xc+
		HALF*(SolnBlk->Grid.Node[i_min+1][j].X+SolnBlk->Grid.Node[i_min+1][j+1].X)+
		HALF*(SolnBlk->Grid.Node[i_min][j+1].X+SolnBlk->Grid.Node[i_min+1][j+1].X)+
		SolnBlk->Grid.Node[i_min+1][j+1].X)/FOUR -
	    SolnBlk->Grid.Cell[i_min][j].Xc;
	  Wfine = Wcoarse +
	    (SolnBlk->phi[i_min][j]^SolnBlk->dWdx[i_min][j])*dX.x +
	    (SolnBlk->phi[i_min][j]^SolnBlk->dWdy[i_min][j])*dX.y;
	  for (k = 0 ; k < blocksize; ++ k) {
	    buffer_count = buffer_count + 1;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Wfine[k+1];
	  } /* endfor */
	  // Evaluate NW sub (fine) cell values.
	  dX = (HALF*(SolnBlk->Grid.Node[i_min][j].X+SolnBlk->Grid.Node[i_min][j+1].X)+
		SolnBlk->Grid.Cell[i_min][j].Xc+
		SolnBlk->Grid.Node[i_min][j+1].X+
		HALF*(SolnBlk->Grid.Node[i_min][j+1].X+SolnBlk->Grid.Node[i_min+1][j+1].X))/FOUR -
	    SolnBlk->Grid.Cell[i_min][j].Xc;
	  Wfine = Wcoarse +
	    (SolnBlk->phi[i_min][j]^SolnBlk->dWdx[i_min][j])*dX.x +
	    (SolnBlk->phi[i_min][j]^SolnBlk->dWdy[i_min][j])*dX.y;
	  for (k = 0 ; k < blocksize; ++ k) {
	    buffer_count = buffer_count + 1;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Wfine[k+1];
	  } /* endfor */
	  // Evaluate SE sub (fine) cell values.
	  dX = (HALF*(SolnBlk->Grid.Node[i_min][j].X+SolnBlk->Grid.Node[i_min+1][j].X)+
		SolnBlk->Grid.Node[i_min+1][j].X+
		SolnBlk->Grid.Cell[i_min][j].Xc+
		HALF*(SolnBlk->Grid.Node[i_min+1][j].X+SolnBlk->Grid.Node[i_min+1][j+1].X))/FOUR -
	    SolnBlk->Grid.Cell[i_min][j].Xc;
	  Wfine = Wcoarse +
	    (SolnBlk->phi[i_min][j]^SolnBlk->dWdx[i_min][j])*dX.x +
	    (SolnBlk->phi[i_min][j]^SolnBlk->dWdy[i_min][j])*dX.y;
	  for (k = 0 ; k < blocksize; ++ k) {
	    buffer_count = buffer_count + 1;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Wfine[k+1];
	  } /* endfor */
	  // Evaluate SW sub (fine) cell values.
	  dX = (SolnBlk->Grid.Node[i_min][j].X+
		HALF*(SolnBlk->Grid.Node[i_min][j].X+SolnBlk->Grid.Node[i_min+1][j].X)+
		HALF*(SolnBlk->Grid.Node[i_min][j].X+SolnBlk->Grid.Node[i_min][j+1].X)+
		SolnBlk->Grid.Cell[i_min][j].Xc)/FOUR -
	    SolnBlk->Grid.Cell[i_min][j].Xc;
	  Wfine = Wcoarse +
	    (SolnBlk->phi[i_min][j]^SolnBlk->dWdx[i_min][j])*dX.x +
	    (SolnBlk->phi[i_min][j]^SolnBlk->dWdy[i_min][j])*dX.y;
	  for (k = 0 ; k < blocksize; ++ k) {
	    buffer_count = buffer_count + 1;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Wfine[k+1];
	  } /* endfor */
	} /* endfor */
      } /* endif */
    } /* endif */
  } /* endif */
  return(0);

}

/**************************************************************************
 * Routine: calculate_pertubed_residual                                   *
 **************************************************************************/
//-------------------------------------------------------------------
// Calculate SolnBlk.U =  SolnBlk.Uo + denormalize( epsilon * W(i) )
//-------------------------------------------------------------------
template <>inline void GMRES_Block<Flame2D_pState,
				   Flame2D_Quad_Block,
				   Flame2D_Input_Parameters>::
calculate_perturbed_residual(const double &epsilon)
{    
  //
  // Loop over grid, including ghost cells
  // 
  for (int j = JCl - Nghost ; j <= JCu + Nghost ; j++) { 
    for (int i = ICl - Nghost ; i <= ICu + Nghost ; i++) {

      // compute perturbed state
      for(int varindex = 0; varindex < blocksize; varindex++){
	SolnBlk->U[i][j][varindex+1] = SolnBlk->Uo[i][j][varindex+1] +
	  denormalizeU( epsilon*W[search_directions*scalar_dim + 
				  index(i,j,varindex)], varindex);
      }     

      //Flame2D spec_check to make sure species (Uo + epsilon*W(i)) > ZERO
      // ie physical for dUdt calc 
      if(!SolnBlk->U[i][j].speciesOK(10)) { 
	cerr<<"\n FAILURE in calculate_perturbed_residual"; 
	exit(1); 
      }

      // Update primitive variables.
      SolnBlk->W[i][j].setU( SolnBlk->U[i][j] );
    }
  }  
}

//-------------------------------------------------------------------
// Calculate SolnBlk.U =  SolnBlk.Uo - denormalize( epsilon * W(i) )
//-------------------------------------------------------------------
template <>inline void GMRES_Block<Flame2D_pState,
				   Flame2D_Quad_Block,
				   Flame2D_Input_Parameters>::
calculate_perturbed_residual_2nd(const double &epsilon)
{    
  //
  // Loop over grid, including ghost cells
  // 
  for (int j = JCl - Nghost ; j <= JCu + Nghost ; j++) { 
    for (int i = ICl - Nghost ; i <= ICu + Nghost ; i++) {      
      
      //copy back R + epsilon * W(i)
      SolnBlk->dUdt[i][j][1] = SolnBlk->dUdt[i][j][0];

      // compute perturbed state
      for(int varindex = 0; varindex < blocksize; varindex++){	
	SolnBlk->U[i][j][varindex+1] = SolnBlk->Uo[i][j][varindex+1] -
	  denormalizeU( epsilon*W[search_directions*scalar_dim + 
				  index(i,j,varindex)], varindex);	    	
      }     

      //Flame2D spec_check to make sure species (Uo + epsilon*W(i)) > ZERO
      // ie physical for dUdt calc 
      if(!SolnBlk->U[i][j].speciesOK(10)) { 
	cerr<<"\n FAILURE in calculate_perturbed_residual"; 
	exit(1); 
      }

      // Update primitive variables.
      SolnBlk->W[i][j].setU( SolnBlk->U[i][j] );
    }
  }  
}

//-------------------------------------------------------------------
// Calculate SolnBlk.U =  SolnBlk.Uo + denormalize( epsilon * x(i) )
//-------------------------------------------------------------------
template <>inline void GMRES_Block<Flame2D_pState,
				   Flame2D_Quad_Block,
				   Flame2D_Input_Parameters>::
calculate_perturbed_residual_Restart(const double &epsilon)
{    
  //
  // Loop over grid, including ghost cells
  // 
  for (int j = JCl - Nghost ; j <= JCu + Nghost ; j++) {
    for (int i = ICl - Nghost ; i <= ICu + Nghost ; i++) {

      // compute perturbed state
      for(int varindex = 0; varindex < blocksize; varindex++){	
	SolnBlk->U[i][j][varindex+1] = SolnBlk->Uo[i][j][varindex+1] + 
	  denormalizeU( epsilon*x[index(i,j,varindex)], varindex);
      }  

      //Flame2D spec_check to make sure species (Uo + epsilon*W(i)) > ZERO
      // ie physical for dUdt calc 
      if(!SolnBlk->U[i][j].speciesOK(10)) { 
	cerr<<"\n FAILURE in calculate_perturbed_residual_Restart "; 
	exit(1); 
      }

      // Update primitive variables.
      SolnBlk->W[i][j].setU( SolnBlk->U[i][j] );
    }
  }  
}

//-------------------------------------------------------------------
// Calculate SolnBlk.U =  SolnBlk.Uo + denormalize( epsilon * x(i) )
//-------------------------------------------------------------------
template <>inline void GMRES_Block<Flame2D_pState,
				   Flame2D_Quad_Block,
				   Flame2D_Input_Parameters>::
calculate_perturbed_residual_2nd_Restart(const double &epsilon)
{    
  //
  // Loop over grid, including ghost cells
  // 
  for (int j = JCl - Nghost ; j <= JCu + Nghost ; j++) {
    for (int i = ICl - Nghost ; i <= ICu + Nghost ; i++) {   
      
      //copy back R + epsilon * W(i)
      SolnBlk->dUdt[i][j][1] = SolnBlk->dUdt[i][j][0];
      
      // compute perturbed state
      for(int varindex = 0; varindex < blocksize; varindex++){	
	SolnBlk->U[i][j][varindex+1] = SolnBlk->Uo[i][j][varindex+1] - 
	  denormalizeU( epsilon*x[index(i,j,varindex)], varindex);
      }  

      //Flame2D spec_check to make sure species (Uo + epsilon*W(i)) > ZERO
      // ie physical for dUdt calc 
      if(!SolnBlk->U[i][j].speciesOK(10)) { 
	cerr<<"\n FAILURE in calculate_perturbed_residual_Restart "; 
	exit(1); 
      }
      
      // Update primitive variables.
      SolnBlk->W[i][j].setU( SolnBlk->U[i][j] );
    }
  }  
}


/********************************************************
 * Routine: calculate_Matrix_Free with Preconditioning  *
 ********************************************************/
template <>inline void GMRES_Block<Flame2D_pState,
				   Flame2D_Quad_Block,
				   Flame2D_Input_Parameters>::
calculate_Matrix_Free(const double &epsilon)
{
  // declares
  static Flame2D_pState Wo;

  //Taking into acount NKS overlap
  int JCl_overlap( 0 ); int JCu_overlap( 0 );
  int ICu_overlap( 0 ); int ICl_overlap( 0 );
  double value( ONE );

  if(overlap){	
    if ( SolnBlk->Grid.BCtypeS[ICl] == BC_NONE)  JCl_overlap = overlap; 
    if ( SolnBlk->Grid.BCtypeN[ICu] == BC_NONE)  JCu_overlap = overlap;
    if ( SolnBlk->Grid.BCtypeE[JCu] == BC_NONE)  ICu_overlap = overlap;
    if ( SolnBlk->Grid.BCtypeW[JCl] == BC_NONE)  ICl_overlap = overlap;
  }
	   
  static DenseMatrix Precon(blocksize,blocksize,ZERO);
  
  // Non-Overlap Ghost Cells R(U) already set to zero by dUdt calculation  
  /* V(i+1) = ( R(U+epsilon*W) - b) / epsilon - (gamma) z / h */
  for (int j = JCl - JCl_overlap; j <= JCu + JCu_overlap; j++) {
    for (int i = ICl - ICl_overlap; i <= ICu + ICu_overlap; i++) {
      
      //-------------------------------------------------------------
      // Low Mach # Preconditioner                   
      //-------------------------------------------------------------
      if(Input_Parameters->Preconditioning){ 
	// preconditioner spacing
	double delta_n = SolnBlk->delta_n(i,j);
	// get initial pState
	Wo.setU( SolnBlk->Uo[i][j] );
	//build preconditioner
	Wo.Low_Mach_Number_Preconditioner(Precon, SolnBlk->Flow_Type, delta_n);       
      }

      //-------------------------------------------------------------
      //Update V 
      //-------------------------------------------------------------
      for(int k =0; k < blocksize; k++){	
	int iter = index(i,j,k);		
	//Matrix Free V(i+1)  
	
	if( Input_Parameters->NKS_IP.GMRES_Frechet_Derivative_Order == FIRST_ORDER ){
	  //forwards differenceing R(U+epsilon) - R(U) / epsilon
	  V[(search_directions+1)*scalar_dim+iter] = 
	    (normalizeR(SolnBlk->dUdt[i][j][0][k+1],k) - 
	     b[iter]) / epsilon ;
	} else if ( Input_Parameters->NKS_IP.GMRES_Frechet_Derivative_Order == SECOND_ORDER ){
	  //2nd order R(U+epsilon) - R(U-epsilon) / 2*epsilon
	  V[(search_directions+1)*scalar_dim+iter] = 
	    normalizeR( SolnBlk->dUdt[i][j][1][k+1] - 
			SolnBlk->dUdt[i][j][0][k+1],k)/(TWO*epsilon);
	}

	//
	//Finite Time Stepping
	//

	// Preconditioner
	if(Input_Parameters->Preconditioning){
	  //gamma(nxn)*z(nx1)/h(1x)
	  value = ZERO;
	  for(int l =0; l < blocksize; l++){
	    value += Precon(k,l) * denormalizeU(W[(search_directions)*scalar_dim + index(i,j,l)],l);
	  }
	  V[(search_directions+1)*scalar_dim+iter] -= 
	    normalizeR(value * LHS_Time<Flame2D_Input_Parameters>(*Input_Parameters,
								  SolnBlk->dt[i][j],
								  DTS_ptr->DTS_dTime),k);   
	  //No Preconditioner
	} else {
	  // z/h
	  V[(search_directions+1)*scalar_dim+iter] -= 
	    normalizeUtoR( W[(search_directions)*scalar_dim + iter] 
			   * LHS_Time<Flame2D_Input_Parameters>(*Input_Parameters,
								SolnBlk->dt[i][j],
								DTS_ptr->DTS_dTime),k);
	}       

// #ifdef _NKS_VERBOSE_NAN_CHECK
// 	// nan check most commonly caused by nans in dUdt !!!!
// 	if (V[(search_directions+1)*scalar_dim+iter] != V[(search_directions+1)*scalar_dim+iter] ){
// 	  cout<<"\n nan in V[ "<<(search_directions+1)*scalar_dim+iter<<"] at "<<i<<" "<<j<<" "<<k
// 	      <<" dUdt "<<  normalizeR(SolnBlk->dUdt[i][j][0][k+1],k) <<" b "<< b[iter]
// 	      <<" z "<<W[(search_directions)*scalar_dim + iter]<< " h "<<( SolnBlk->dt[i][j]*ao);
// 	}
// #endif
      } 
    } 
  } 
}

/****************************************************************
 * Routine: calculate_Matrix_Free_Restart with Preconditioning  *
 ****************************************************************/
template <>inline void GMRES_Block<Flame2D_pState,
				   Flame2D_Quad_Block,
				   Flame2D_Input_Parameters>::
calculate_Matrix_Free_Restart(const double &epsilon)
{
  // declares
  static Flame2D_pState Wo;

  //Taking into acount NKS overlap
  int JCl_overlap( 0 ); int JCu_overlap( 0 );
  int ICu_overlap( 0 ); int ICl_overlap( 0 );
  double value( ONE );

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
      
      //-------------------------------------------------------------
      // Low Mach # Preconditioner                   
      //-------------------------------------------------------------
      if(Input_Parameters->Preconditioning){ 
	// preconditioner spacing
	double delta_n( SolnBlk->delta_n(i, j) );
	// get initial pState
	Wo.setU( SolnBlk->Uo[i][j] );
	//build preconditioner
	Wo.Low_Mach_Number_Preconditioner(Precon, SolnBlk->Flow_Type, delta_n);       
      }

      //-------------------------------------------------------------
      //Update V 
      //-------------------------------------------------------------
      for(int k =0; k < blocksize; k++){	
	int iter = index(i,j,k);		
	//Matrix Free V(i+1) 

	if( Input_Parameters->NKS_IP.GMRES_Frechet_Derivative_Order == FIRST_ORDER ){
	  V[iter] = (normalizeR(SolnBlk->dUdt[i][j][0][k+1],k) - b[iter]) / epsilon ;
	} else if ( Input_Parameters->NKS_IP.GMRES_Frechet_Derivative_Order == SECOND_ORDER ){
	  V[iter] = normalizeR( SolnBlk->dUdt[i][j][1][k+1] - SolnBlk->dUdt[i][j][0][k+1],k)/(TWO*epsilon);
	}

	//
	//Finite Time Stepping
	//

	// Preconditioner
	if(Input_Parameters->Preconditioning){
	  //gamma(nxn)*x(nx1)/h(1x)    
	  value = ZERO;
	  for(int l =0; l < blocksize; l++){
	    value += Precon(k,l) * denormalizeU( x[index(i,j,l)],l);
	  }
	  V[iter] -= normalizeR(value * LHS_Time<Flame2D_Input_Parameters>(*Input_Parameters,
									   SolnBlk->dt[i][j],
									   DTS_ptr->DTS_dTime),k);
	  
	  //No Preconditioner
	} else { 
	  // z/h
	  V[iter] -= normalizeUtoR(x[iter] * LHS_Time<Flame2D_Input_Parameters>(*Input_Parameters,
										SolnBlk->dt[i][j],
										DTS_ptr->DTS_dTime),k);	
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
template<> inline void GMRES_Block<Flame2D_pState,
				   Flame2D_Quad_Block,
				   Flame2D_Input_Parameters>::
set_normalize_values(void)
{   

  static const Flame2D_pState W_STD_ATM;
  static double ao( W_STD_ATM.a() );
  static double rho( W_STD_ATM.rho() );

  // Original Normalization from code  
  normalize_valuesU[0] = rho;          //rho
  normalize_valuesU[1] = rho*ao;       //rho*u
  normalize_valuesU[2] = rho*ao;       //rho*v
  normalize_valuesU[3] = rho*ao*ao;    //rho*e
  for(int i=NUM_FLAME2D_VAR_SANS_SPECIES; i < blocksize; i++){
    normalize_valuesU[i] = rho;        //species mass fraction rho*cs
  }

  normalize_valuesR[0] = rho*ao;          //rho
  normalize_valuesR[1] = rho*ao*ao;       //rho*u
  normalize_valuesR[2] = rho*ao*ao;       //rho*v
  normalize_valuesR[3] = rho*ao*ao*ao;    //rho*e
  for(int i=NUM_FLAME2D_VAR_SANS_SPECIES; i < blocksize; i++){
    normalize_valuesR[i] = rho*ao;        //species mass fraction rho*cs
  }
 
}

#endif // _FLAME2D_NKS_INCLUDED 







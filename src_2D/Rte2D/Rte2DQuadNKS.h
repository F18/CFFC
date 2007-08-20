#ifndef _RTE2D_NKS_INCLUDED 
#define _RTE2D_NKS_INCLUDED 

/****************************************************************/
/****** NKS REQUIRED SPECIALIZATIONS & FUNCTIONS ****************/
/****************************************************************/  
#include "Rte2DQuad.h"
#include "../NewtonKrylovSchwarz2D/NKS2D.h"


// Prototypes required for instantiation


/*  GMRES FUNCTION PROTOTYPES */

template <>
void GMRES_Block<Rte2D_State,Rte2D_Quad_Block,Rte2D_Input_Parameters>::
  SubcellReconstruction(const int i,
			const int j,
			const int Limiter);



/*  BLOCK PRECONDITIONER FUNCTION PROTOTYPES */

void dFIdU_FD(DenseMatrix& dRdW, Rte2D_Quad_Block &SolnBlk,  
	      Rte2D_Input_Parameters &Input_Parameters,
	      const int ii, const int jj, const int Orient);

/****************************************************************/
/****** NKS REQUIRED SPECIALIZATIONS & FUNCTIONS ****************/
/****************************************************************/  



/*! *****************************************************************************************
 * Specialization of Newton_Update                                                          *
 *                                                                                          *
 * This routine updates the previous Solution Data Uo with the deltaU from the GMRES        *                  
 * iterative solver.  U = Uo + GMRES.delatU                                                 *
 *                                                                                          *
 ********************************************************************************************/
template <>
int Newton_Update(Rte2D_Quad_Block *SolnBlk,
		  AdaptiveBlock2D_List &List_of_Local_Solution_Blocks,
		  Rte2D_Input_Parameters &Input_Parameters,
		  GMRES_RightPrecon_MatrixFree<Rte2D_State,Rte2D_Quad_Block,Rte2D_Input_Parameters> &GMRES) {

  int Num_Var = SolnBlk[0].NumVar();  
  int error_flag = 0;
 
  /* Update Solution. No updates to Ghost Cells, let the BC's take care of it */
  for ( int Bcount = 0 ; Bcount < List_of_Local_Solution_Blocks.Nblk ; ++Bcount ) {
    if (List_of_Local_Solution_Blocks.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
      for (int j = SolnBlk[Bcount].JCl; j <= SolnBlk[Bcount].JCu; j++){
	for (int i = SolnBlk[Bcount].ICl; i <= SolnBlk[Bcount].ICu; i++){
	  
	  /* Update solutions in conversed variables  U = Uo + deltaU = Uo + denormalized(x) */	 
	  for(int varindex =1; varindex <= Num_Var; varindex++){  
	    SolnBlk[Bcount].U[i][j][varindex] = SolnBlk[Bcount].Uo[i][j][varindex] 
	      +  GMRES.deltaU(Bcount,i,j,varindex-1);
	  } 	      	  
	  // THIS FUNCTION HAS NO CHECKS FOR INVALID SOLUTIONS, 
	  // YOU PROBABLY WANT TO CREATE A SPECIALIZATION OF THIS FUNCTION SPECIFIC 
	  // FOR YOUR EQUATION SYSTEM see Euler2D, Chem2D, etc...
	  /***********************************************************************
	   *************************** RTE SPECIFIC ******************************/
 	  // Apply update reduction while any one of the updated variables is unphysical 
 	  if(! SolnBlk[Bcount].U[i][j].Unphysical_Properties()){	   
 	    double update_reduction_factor = ONE;	    
 	    for (int n_update_reduction = 1; n_update_reduction <= 10; ++n_update_reduction) {		  
 	      update_reduction_factor = HALF*update_reduction_factor;		  		  
 	      for(int varindex = 1; varindex <= Num_Var; varindex++){		              
 		SolnBlk[Bcount].U[i][j][varindex] = SolnBlk[Bcount].Uo[i][j][varindex] 
 		  + GMRES.deltaU(Bcount,i,j,varindex-1)*update_reduction_factor;
 	      }   
 	      cout<<"\n Applying Reduction to solution in NKS "<<n_update_reduction;
 	      if( SolnBlk[Bcount].U[i][j].Unphysical_Properties() )  break;	      
 	    } 
 	  }
	  // Error Check
 	  if(! SolnBlk[Bcount].U[i][j].Unphysical_Properties()) error_flag = 1;
	  /*************************************************************************
	   *************************************************************************/
	} 
      } 
    } 
  } 
  
  return (error_flag); 
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
template<> inline void GMRES_Block<Rte2D_State,
				   Rte2D_Quad_Block,
				   Rte2D_Input_Parameters>::
set_normalize_values(void)
{   

  // NOT NORMALIZED
  for(int i=0; i < blocksize; i++) {
    normalize_valuesU[i] = ONE;
    normalize_valuesR[i] = ONE;
  }
}




/**************************************************************************
 * Routine: TESTING FUNCITON                                              *
 **************************************************************************/
template <> 
inline void GMRES_Block<Rte2D_State,
                        Rte2D_Quad_Block,
                        Rte2D_Input_Parameters>::
Output_U(int what)
{    
  for (int j = JCl - Nghost ; j <= JCu + Nghost ; j++) {  //includes ghost cells 
    for (int i = ICl - Nghost ; i <= ICu + Nghost ; i++) {
      if(what == 1)cout<<"\n U ("<<i<<","<<j<<") = "<<SolnBlk->U[i][j];      
      /***********************************************************************
       *************************** RTE SPECIFIC ******************************/
      // removed output primitive variables
      // if(what == 2)cout<<"\n W ("<<i<<","<<j<<") = "<<SolnBlk->W[i][j];
      /*************************************************************************
       *************************************************************************/
      if(what == 3)cout<<"\n dUdt ("<<i<<","<<j<<") = "<<SolnBlk->dUdt[i][j][0];
    }
  }  
}



/**************************************************************************
 * Routine: calculate_pertubed_residual                                   *
 **************************************************************************/
// Calculate Soln_ptr.U =  Soln_ptr.Uo + denormalize( epsilon * W(i) )
template <> inline void GMRES_Block<Rte2D_State,
                                    Rte2D_Quad_Block,
                                    Rte2D_Input_Parameters>::
calculate_perturbed_residual(const double &epsilon)
{    
  for (int j = JCl - Nghost ; j <= JCu + Nghost ; j++) {  //includes ghost cells 
    for (int i = ICl - Nghost ; i <= ICu + Nghost ; i++) {
      for(int varindex = 0; varindex < blocksize; varindex++){	
	SolnBlk->U[i][j][varindex+1] = SolnBlk->Uo[i][j][varindex+1] + 
	  denormalizeU( epsilon*W[search_directions*scalar_dim+index(i,j,varindex)], varindex); 
      }   
      /* Update primitive variables. */
      /***********************************************************************
       *************************** RTE SPECIFIC ******************************/
      // removed update primitive variables
      // SolnBlk->W[i][j] = SolnBlk->U[i][j].W();      
      /*************************************************************************
       *************************************************************************/
    }
  }  
}

// Copy forward difference & calculate backwards for 2nd order derivative 
template <> inline void GMRES_Block<Rte2D_State,
                                    Rte2D_Quad_Block,
                                    Rte2D_Input_Parameters>::
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
      /* Update primitive variables. */
      /***********************************************************************
       *************************** RTE SPECIFIC ******************************/
      // removed update primitive variables
      // SolnBlk->W[i][j] = SolnBlk->U[i][j].W();      
      /*************************************************************************
       *************************************************************************/
    }
  }  
}

// Calculate Restart Soln_ptr.U =  Soln_ptr.Uo + denormalize( epsilon * x(i) )
template <> inline void GMRES_Block<Rte2D_State,
                                    Rte2D_Quad_Block,
                                    Rte2D_Input_Parameters>::
calculate_perturbed_residual_Restart(const double &epsilon)
{    
  for (int j = JCl - Nghost ; j <= JCu + Nghost ; j++) {
    for (int i = ICl - Nghost ; i <= ICu + Nghost ; i++) {
      for(int varindex = 0; varindex < blocksize; varindex++){	
	SolnBlk->U[i][j][varindex+1] = SolnBlk->Uo[i][j][varindex+1] + 
	  denormalizeU( epsilon*x[index(i,j,varindex)], varindex);
      }  
      /* Update primitive variables. */
      /***********************************************************************
       *************************** RTE SPECIFIC ******************************/
      // removed update primitive variables
      // SolnBlk->W[i][j] = SolnBlk->U[i][j].W();      
      /*************************************************************************
       *************************************************************************/
    }
  }  
}

// Copy forward difference & calculate backwards for 2nd order derivative  
template <> inline void GMRES_Block<Rte2D_State,
                                    Rte2D_Quad_Block,
                                    Rte2D_Input_Parameters>::
calculate_perturbed_residual_2nd_Restart(const double &epsilon)
{    
  for (int j = JCl - Nghost ; j <= JCu + Nghost ; j++) {
    for (int i = ICl - Nghost ; i <= ICu + Nghost ; i++) { 
      //copy back R + epsilon * W(i)
      SolnBlk->dUdt[i][j][1] = SolnBlk->dUdt[i][j][0];
      
      for(int varindex = 0; varindex < blocksize; varindex++){	
	SolnBlk->U[i][j][varindex+1] = SolnBlk->Uo[i][j][varindex+1] -
	  denormalizeU( epsilon*x[index(i,j,varindex)], varindex);
      }  
      /* Update primitive variables. */
      /***********************************************************************
       *************************** RTE SPECIFIC ******************************/
      // removed update primitive variables
      // SolnBlk->W[i][j] = SolnBlk->U[i][j].W();      
      /*************************************************************************
       *************************************************************************/
    }
  }  
}






/*******************************************************************************
 * GMRES_Block::LoadSendBuffer_F2C -- Loads send message buffer for            *
 *                                    fine to coarse block message             *
 *                                    passing.                                 *
 *******************************************************************************/
template <> 
inline int GMRES_Block<Rte2D_State,Rte2D_Quad_Block,Rte2D_Input_Parameters>::
LoadSendBuffer_F2C(double *buffer,
		   int &buffer_count,
		   const int buffer_size,
		   const int i_min, 
		   const int i_max,
		   const int i_inc,
		   const int j_min, 
		   const int j_max,
		   const int j_inc) {
  
  for (int j  = j_min ; ((j_inc+2)/4) ? (j < j_max):(j > j_max) ; j += j_inc ) {
    for (int i = i_min ;  ((i_inc+2)/4) ? (i < i_max):(i > i_max) ; i += i_inc ) {
      for (int k = 0 ; k < blocksize; ++ k) {
	buffer_count++;
	if (buffer_count >= buffer_size) return(1);
	if (vector_switch) {
	  /***********************************************************************
	   *************************** RTE SPECIFIC ******************************
	  // ORIGINAL : 
	  buffer[buffer_count] = (SolnBlk->Grid.Cell[i  ][j  ].A*W[(search_directions)*scalar_dim + index(i  ,j  ,k)]+
				  SolnBlk->Grid.Cell[i+1][j  ].A*W[(search_directions)*scalar_dim + index(i+1,j  ,k)]+
				  SolnBlk->Grid.Cell[i  ][j+1].A*W[(search_directions)*scalar_dim + index(i  ,j+1,k)]+
				  SolnBlk->Grid.Cell[i+1][j+1].A*W[(search_directions)*scalar_dim + index(i+1,j+1,k)])/
	                         (SolnBlk->Grid.Cell[i  ][j  ].A+
				  SolnBlk->Grid.Cell[i+1][j  ].A+
				  SolnBlk->Grid.Cell[i  ][j+1].A+
				  SolnBlk->Grid.Cell[i+1][j+1].A);
	  *************************************************************************/
	  // RTE2D: Added Quasi-3D scaling param 
	  buffer[buffer_count] = (SolnBlk->Grid.Cell[i  ][j  ].A*SolnBlk->Sp[i  ][j  ]*
				  W[(search_directions)*scalar_dim + index(i  ,j  ,k)]+
				  SolnBlk->Grid.Cell[i+1][j  ].A*SolnBlk->Sp[i+1][j  ]*
				  W[(search_directions)*scalar_dim + index(i+1,j  ,k)]+
				  SolnBlk->Grid.Cell[i  ][j+1].A*SolnBlk->Sp[i  ][j+1]*
				  W[(search_directions)*scalar_dim + index(i  ,j+1,k)]+
				  SolnBlk->Grid.Cell[i+1][j+1].A*SolnBlk->Sp[i+1][j+1]*
				  W[(search_directions)*scalar_dim + index(i+1,j+1,k)])/
	                         (SolnBlk->Grid.Cell[i  ][j  ].A*SolnBlk->Sp[i  ][j  ]+
				  SolnBlk->Grid.Cell[i+1][j  ].A*SolnBlk->Sp[i+1][j  ]+
				  SolnBlk->Grid.Cell[i  ][j+1].A*SolnBlk->Sp[i  ][j+1]+
				  SolnBlk->Grid.Cell[i+1][j+1].A*SolnBlk->Sp[i+1][j+1]);
	  /*************************************************************************
	   *************************************************************************/
	} else {
	  /***********************************************************************
	   *************************** RTE SPECIFIC ******************************
	  // ORIGINAL : 
	  buffer[buffer_count] = (SolnBlk->Grid.Cell[i  ][j  ].A*x[index(i  ,j  ,k)]+
				  SolnBlk->Grid.Cell[i+1][j  ].A*x[index(i+1,j  ,k)]+
				  SolnBlk->Grid.Cell[i  ][j+1].A*x[index(i  ,j+1,k)]+
				  SolnBlk->Grid.Cell[i+1][j+1].A*x[index(i+1,j+1,k)])/
	                         (SolnBlk->Grid.Cell[i  ][j  ].A+
				  SolnBlk->Grid.Cell[i+1][j  ].A+
				  SolnBlk->Grid.Cell[i  ][j+1].A+
				  SolnBlk->Grid.Cell[i+1][j+1].A);
	  *************************************************************************/
	  // RTE2D: Added Quasi-3D scaling param 
	  buffer[buffer_count] = (SolnBlk->Grid.Cell[i  ][j  ].A*SolnBlk->Sp[i  ][j  ]*x[index(i  ,j  ,k)]+
				  SolnBlk->Grid.Cell[i+1][j  ].A*SolnBlk->Sp[i+1][j  ]*x[index(i+1,j  ,k)]+
				  SolnBlk->Grid.Cell[i  ][j+1].A*SolnBlk->Sp[i  ][j+1]*x[index(i  ,j+1,k)]+
				  SolnBlk->Grid.Cell[i+1][j+1].A*SolnBlk->Sp[i+1][j+1]*x[index(i+1,j+1,k)])/
	                         (SolnBlk->Grid.Cell[i  ][j  ].A*SolnBlk->Sp[i  ][j  ]+
				  SolnBlk->Grid.Cell[i+1][j  ].A*SolnBlk->Sp[i+1][j  ]+
				  SolnBlk->Grid.Cell[i  ][j+1].A*SolnBlk->Sp[i  ][j+1]+
				  SolnBlk->Grid.Cell[i+1][j+1].A*SolnBlk->Sp[i+1][j+1]);
	  /*************************************************************************
	   *************************************************************************/
	} 
      } 
    } 
  } 
  return(0);
}




/*******************************************************************************
 * GMRES_Block::LoadSendBuffer_C2F -- Loads send message buffer for            *
 *                                    coarse to fine block message             *
 *                                    passing.                                 *
 *******************************************************************************/
template <> 
inline int GMRES_Block<Rte2D_State,Rte2D_Quad_Block,Rte2D_Input_Parameters>::
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
  Vector2D dX;
  Rte2D_State Wcoarse, Wfine;
  int LIMITER = LIMITER_ZERO; //LIMITER_VENKATAKRISHNAN 

      /***********************************************************************
       *************************** RTE SPECIFIC ******************************/
      // changed dWdx -> dUdx, dWdy -> dUdy since Rte2D has no primitive class
      /*************************************************************************
       *************************************************************************/


  if (j_min == j_max) { // North or south boundary.
     // Four different orderings to consider depending on the value of i_inc & j_inc.
    if (j_inc > 0) {             
      if (i_inc > 0) {

	/******************************* CASE #1 ***************************************/
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
	    (SolnBlk->phi[i][j_min]^SolnBlk->dUdx[i][j_min])*dX.x +
	    (SolnBlk->phi[i][j_min]^SolnBlk->dUdy[i][j_min])*dX.y;

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
	    (SolnBlk->phi[i][j_min]^SolnBlk->dUdx[i][j_min])*dX.x +
	    (SolnBlk->phi[i][j_min]^SolnBlk->dUdy[i][j_min])*dX.y;
            
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
                      (SolnBlk->phi[i][j_min]^SolnBlk->dUdx[i][j_min])*dX.x +
                      (SolnBlk->phi[i][j_min]^SolnBlk->dUdy[i][j_min])*dX.y;
          
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
                      (SolnBlk->phi[i][j_min]^SolnBlk->dUdx[i][j_min])*dX.x +
                      (SolnBlk->phi[i][j_min]^SolnBlk->dUdy[i][j_min])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
           } /* endfor */
	
	/******************************* CASE #2 ***************************************/
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
                      (SolnBlk->phi[i][j_min]^SolnBlk->dUdx[i][j_min])*dX.x +
                      (SolnBlk->phi[i][j_min]^SolnBlk->dUdy[i][j_min])*dX.y;
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
                      (SolnBlk->phi[i][j_min]^SolnBlk->dUdx[i][j_min])*dX.x +
                      (SolnBlk->phi[i][j_min]^SolnBlk->dUdy[i][j_min])*dX.y;
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
                      (SolnBlk->phi[i][j_min]^SolnBlk->dUdx[i][j_min])*dX.x +
                      (SolnBlk->phi[i][j_min]^SolnBlk->dUdy[i][j_min])*dX.y;
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
                      (SolnBlk->phi[i][j_min]^SolnBlk->dUdx[i][j_min])*dX.x +
                      (SolnBlk->phi[i][j_min]^SolnBlk->dUdy[i][j_min])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
           } /* endfor */
        } /* endif */

      /******************************* CASE #3 ***************************************/
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
                      (SolnBlk->phi[i][j_min]^SolnBlk->dUdx[i][j_min])*dX.x +
                      (SolnBlk->phi[i][j_min]^SolnBlk->dUdy[i][j_min])*dX.y;
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
                      (SolnBlk->phi[i][j_min]^SolnBlk->dUdx[i][j_min])*dX.x +
                      (SolnBlk->phi[i][j_min]^SolnBlk->dUdy[i][j_min])*dX.y;
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
                      (SolnBlk->phi[i][j_min]^SolnBlk->dUdx[i][j_min])*dX.x +
                      (SolnBlk->phi[i][j_min]^SolnBlk->dUdy[i][j_min])*dX.y;
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
                      (SolnBlk->phi[i][j_min]^SolnBlk->dUdx[i][j_min])*dX.x +
                      (SolnBlk->phi[i][j_min]^SolnBlk->dUdy[i][j_min])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
           } /* endfor */
	   /******************************* CASE #4 ***************************************/
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
                      (SolnBlk->phi[i][j_min]^SolnBlk->dUdx[i][j_min])*dX.x +
                      (SolnBlk->phi[i][j_min]^SolnBlk->dUdy[i][j_min])*dX.y;

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
                      (SolnBlk->phi[i][j_min]^SolnBlk->dUdx[i][j_min])*dX.x +
                      (SolnBlk->phi[i][j_min]^SolnBlk->dUdy[i][j_min])*dX.y;
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
                      (SolnBlk->phi[i][j_min]^SolnBlk->dUdx[i][j_min])*dX.x +
                      (SolnBlk->phi[i][j_min]^SolnBlk->dUdy[i][j_min])*dX.y;
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
                      (SolnBlk->phi[i][j_min]^SolnBlk->dUdx[i][j_min])*dX.x +
                      (SolnBlk->phi[i][j_min]^SolnBlk->dUdy[i][j_min])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
           } /* endfor */
        } /* endif */
     } /* endif */


    /******************************* CASE #5 ***************************************/
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
                      (SolnBlk->phi[i_min][j]^SolnBlk->dUdx[i_min][j])*dX.x +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dUdy[i_min][j])*dX.y;
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
                      (SolnBlk->phi[i_min][j]^SolnBlk->dUdx[i_min][j])*dX.x +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dUdy[i_min][j])*dX.y;
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
                      (SolnBlk->phi[i_min][j]^SolnBlk->dUdx[i_min][j])*dX.x +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dUdy[i_min][j])*dX.y;
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
                      (SolnBlk->phi[i_min][j]^SolnBlk->dUdx[i_min][j])*dX.x +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dUdy[i_min][j])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
           } /* endfor */
	   /******************************* CASE #6 ***************************************/
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
                      (SolnBlk->phi[i_min][j]^SolnBlk->dUdx[i_min][j])*dX.x +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dUdy[i_min][j])*dX.y;
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
                      (SolnBlk->phi[i_min][j]^SolnBlk->dUdx[i_min][j])*dX.x +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dUdy[i_min][j])*dX.y;
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
                      (SolnBlk->phi[i_min][j]^SolnBlk->dUdx[i_min][j])*dX.x +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dUdy[i_min][j])*dX.y;
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
                      (SolnBlk->phi[i_min][j]^SolnBlk->dUdx[i_min][j])*dX.x +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dUdy[i_min][j])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
           } /* endfor */
        } /* endif */	
	/******************************* CASE #7 ***************************************/
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
                      (SolnBlk->phi[i_min][j]^SolnBlk->dUdx[i_min][j])*dX.x +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dUdy[i_min][j])*dX.y;
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
                      (SolnBlk->phi[i_min][j]^SolnBlk->dUdx[i_min][j])*dX.x +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dUdy[i_min][j])*dX.y;
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
                      (SolnBlk->phi[i_min][j]^SolnBlk->dUdx[i_min][j])*dX.x +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dUdy[i_min][j])*dX.y;
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
                      (SolnBlk->phi[i_min][j]^SolnBlk->dUdx[i_min][j])*dX.x +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dUdy[i_min][j])*dX.y;
              for (k = 0 ; k < blocksize; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k+1];
              } /* endfor */
           } /* endfor */\

	   /******************************* CASE #8 ***************************************/
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
                      (SolnBlk->phi[i_min][j]^SolnBlk->dUdx[i_min][j])*dX.x +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dUdy[i_min][j])*dX.y;
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
                      (SolnBlk->phi[i_min][j]^SolnBlk->dUdx[i_min][j])*dX.x +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dUdy[i_min][j])*dX.y;
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
                      (SolnBlk->phi[i_min][j]^SolnBlk->dUdx[i_min][j])*dX.x +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dUdy[i_min][j])*dX.y;
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
                      (SolnBlk->phi[i_min][j]^SolnBlk->dUdx[i_min][j])*dX.x +
                      (SolnBlk->phi[i_min][j]^SolnBlk->dUdy[i_min][j])*dX.y;
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
 * GMRES_Block::SubcellReconstruction --                                  *
 *              Performs the subcell reconstruction of solution state     *
 *              within a given cell (i,j) of the computational mesh for   *
 *              the specified quadrilateral solution block.               *
 **************************************************************************/
template <> 
inline void GMRES_Block<Rte2D_State,Rte2D_Quad_Block,Rte2D_Input_Parameters>::
SubcellReconstruction(const int i, 
		      const int j,
		      const int Limiter) {
  
  /***********************************************************************
   *************************** RTE SPECIFIC ******************************/
  // changed dWdx -> dUdx, dWdy -> dUdy since Rte2D has no primitive class
  // also changec Vacuum() calls to Zero() calls.
  /*************************************************************************
   *************************************************************************/



  int n, n2, n_pts, i_index[8], j_index[8], k;
  double u0, u0Min, u0Max, uQuad[4], phi_n;
  double DxDx_ave, DxDy_ave, DyDy_ave;
  Vector2D dX;
  Rte2D_State U0, DU, DUDx_ave, DUDy_ave, W_VACUUM;
  W_VACUUM.Zero();
  
  /* Carry out the limited solution reconstruction in
     each cell of the computational mesh. */

   if (i == ICl-Nghost || i == ICu+Nghost ||
      j == JCl-Nghost || j == JCu+Nghost) {
    n_pts = 0;
    } else if ((i == ICl-Nghost+1) && 
               (Grid.BCtypeW[j] != BC_NONE)) {
      if (j == JCl-Nghost+1 || j == JCu+Nghost-1) {
         n_pts = 0;
      } else if (SolnBlk->Grid.BCtypeW[j] == BC_PERIODIC ||
                 SolnBlk->Grid.BCtypeW[j] == BC_CONSTANT_EXTRAPOLATION ||
                 SolnBlk->Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION ||
                 SolnBlk->Grid.BCtypeW[j] == BC_GRAY_WALL ) {
         if (j == JCl) {
            n_pts = 5;
            i_index[0] = i-1; j_index[0] = j  ;
            i_index[1] = i+1; j_index[1] = j  ;
            i_index[2] = i-1; j_index[2] = j+1;
            i_index[3] = i  ; j_index[3] = j+1;
            i_index[4] = i+1; j_index[4] = j+1;
         } else if (j == JCu) {
            n_pts = 5;
            i_index[0] = i-1; j_index[0] = j-1;
            i_index[1] = i  ; j_index[1] = j-1;
            i_index[2] = i+1; j_index[2] = j-1;
            i_index[3] = i-1; j_index[3] = j  ;
            i_index[4] = i+1; j_index[4] = j  ;
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
      } else {
         if (j == JCl) {
            n_pts = 3;
            i_index[0] = i+1; j_index[0] = j  ;
            i_index[1] = i  ; j_index[1] = j+1;
            i_index[2] = i+1; j_index[2] = j+1;
         } else if (j == JCu) {
            n_pts = 3;
            i_index[0] = i  ; j_index[0] = j-1;
            i_index[1] = i+1; j_index[1] = j-1;
            i_index[2] = i+1; j_index[2] = j  ;
         } else {
            n_pts = 5;
            i_index[0] = i  ; j_index[0] = j-1;
            i_index[1] = i+1; j_index[1] = j-1;
            i_index[2] = i+1; j_index[2] = j  ;
            i_index[3] = i  ; j_index[3] = j+1;
            i_index[4] = i+1; j_index[4] = j+1;
         } /* endif */
      } /* endif */           
    } else if ((i == ICu+Nghost-1) && 
               (SolnBlk->Grid.BCtypeE[j] != BC_NONE)) {
      if (j == JCl-Nghost+1 || j == JCu+Nghost-1) {
         n_pts = 0;
      } else if (SolnBlk->Grid.BCtypeE[j] == BC_PERIODIC ||
                 SolnBlk->Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
                 SolnBlk->Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
                 SolnBlk->Grid.BCtypeE[j] == BC_GRAY_WALL ) {
         if (j == JCl) {
            n_pts = 5;
            i_index[0] = i-1; j_index[0] = j  ;
            i_index[1] = i+1; j_index[1] = j  ;
            i_index[2] = i-1; j_index[2] = j+1;
            i_index[3] = i  ; j_index[3] = j+1;
            i_index[4] = i+1; j_index[4] = j+1;
         } else if (j == JCu) {
            n_pts = 5;
            i_index[0] = i-1; j_index[0] = j-1;
            i_index[1] = i  ; j_index[1] = j-1;
            i_index[2] = i+1; j_index[2] = j-1;
            i_index[3] = i-1; j_index[3] = j  ;
            i_index[4] = i+1; j_index[4] = j  ;
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
      } else {
         if (j == JCl) {
            n_pts = 3;
            i_index[0] = i-1; j_index[0] = j  ;
            i_index[1] = i-1; j_index[1] = j+1;
            i_index[2] = i  ; j_index[2] = j+1;
         } else if (j == JCu) {
            n_pts = 3;
            i_index[0] = i-1; j_index[0] = j-1;
            i_index[1] = i  ; j_index[1] = j-1;
            i_index[2] = i-1; j_index[2] = j  ;
         } else {
            n_pts = 5;
            i_index[0] = i-1; j_index[0] = j-1;
            i_index[1] = i  ; j_index[1] = j-1;
            i_index[2] = i-1; j_index[2] = j  ;
            i_index[3] = i-1; j_index[3] = j+1;
            i_index[4] = i  ; j_index[4] = j+1;
         } /* endif */
      } /* endif */
    } else if ((j == JCl-Nghost+1) && 
               (SolnBlk->Grid.BCtypeS[i] != BC_NONE)) {
      if (i == ICl-Nghost+1 || i == ICu+Nghost-1) {
         n_pts = 0;
      } else if (SolnBlk->Grid.BCtypeS[i] == BC_PERIODIC ||
                 SolnBlk->Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
                 SolnBlk->Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
                 SolnBlk->Grid.BCtypeS[i] == BC_GRAY_WALL ) {
         if (i == ICl) {
            n_pts = 5;
            i_index[0] = i  ; j_index[0] = j-1;
            i_index[1] = i+1; j_index[1] = j-1;
            i_index[2] = i+1; j_index[2] = j  ;
            i_index[3] = i  ; j_index[3] = j+1;
            i_index[4] = i+1; j_index[4] = j+1;
         } else if (i == ICu) {
            n_pts = 5;
            i_index[0] = i-1; j_index[0] = j-1;
            i_index[1] = i  ; j_index[1] = j-1;
            i_index[2] = i-1; j_index[2] = j  ;
            i_index[3] = i-1; j_index[3] = j+1;
            i_index[4] = i  ; j_index[4] = j+1;
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
      } else {
         if (i == ICl) {
            n_pts = 3;
            i_index[0] = i+1; j_index[0] = j  ;
            i_index[1] = i  ; j_index[1] = j+1;
            i_index[2] = i+1; j_index[2] = j+1;
         } else if (i == ICu) {
            n_pts = 3;
            i_index[0] = i-1; j_index[0] = j  ;
            i_index[1] = i-1; j_index[1] = j+1;
            i_index[2] = i  ; j_index[2] = j+1;
         } else {
            n_pts = 5;
            i_index[0] = i-1; j_index[0] = j  ;
            i_index[1] = i+1; j_index[1] = j  ;
            i_index[2] = i-1; j_index[2] = j+1;
            i_index[3] = i  ; j_index[3] = j+1;
            i_index[4] = i+1; j_index[4] = j+1;
         } /* endif */
      } /* endif */
    } else if ((j == JCu+Nghost-1) && 
               (SolnBlk->Grid.BCtypeN[i] != BC_NONE)) {
      if (i == ICl-Nghost+1 || i == ICu+Nghost-1) {
         n_pts = 0;
      } else if (SolnBlk->Grid.BCtypeN[i] == BC_PERIODIC ||
                 SolnBlk->Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
                 SolnBlk->Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
                 SolnBlk->Grid.BCtypeN[i] == BC_GRAY_WALL ) {
         if (i == ICl) {
            n_pts = 5;
            i_index[0] = i  ; j_index[0] = j-1;
            i_index[1] = i+1; j_index[1] = j-1;
            i_index[2] = i+1; j_index[2] = j  ;
            i_index[3] = i  ; j_index[3] = j+1;
            i_index[4] = i+1; j_index[4] = j+1;
         } else if (i == ICu) {
            n_pts = 5;
            i_index[0] = i-1; j_index[0] = j-1;
            i_index[1] = i  ; j_index[1] = j-1;
            i_index[2] = i-1; j_index[2] = j  ;
            i_index[3] = i-1; j_index[3] = j+1;
            i_index[4] = i  ; j_index[4] = j+1;
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
      } else {
         if (i == ICl) {
            n_pts = 3;
            i_index[0] = i  ; j_index[0] = j-1;
            i_index[1] = i+1; j_index[1] = j-1;
            i_index[2] = i+1; j_index[2] = j  ;
         } else if (i == ICu) {
            n_pts = 3;
            i_index[0] = i-1; j_index[0] = j-1;
            i_index[1] = i  ; j_index[1] = j-1;
            i_index[2] = i-1; j_index[2] = j  ;
         } else {
            n_pts = 5;
            i_index[0] = i-1; j_index[0] = j-1;
            i_index[1] = i  ; j_index[1] = j-1;
            i_index[2] = i+1; j_index[2] = j-1;
            i_index[3] = i-1; j_index[3] = j  ;
            i_index[4] = i+1; j_index[4] = j  ;
         } /* endif */
      } /* endif */
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
      DUDx_ave = W_VACUUM;
      DUDy_ave = W_VACUUM;
      DxDx_ave = ZERO;
      DxDy_ave = ZERO;
      DyDy_ave = ZERO;

      for ( k = 1 ; k <= blocksize; ++ k) {
         if (vector_switch) {
            U0[k] = W[(search_directions)*scalar_dim + index(i,j,k-1)];
         } else {
            U0[k] = x[index(i,j,k-1)];
         } /* endif */
      } /* endfor */

      for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
          dX = SolnBlk->Grid.Cell[ i_index[n2] ][ j_index[n2] ].Xc - 
               SolnBlk->Grid.Cell[i][j].Xc;
          for ( k = 1 ; k <= blocksize; ++ k) {
             if (vector_switch) {
                DU[k] = W[(search_directions)*scalar_dim + index(i_index[n2] , j_index[n2] , k-1)] - 
                        U0[k];
             } else {
                DU[k] = x[index( i_index[n2] , j_index[n2] , k-1)] - 
                        U0[k];
             } /* endif */
          } /* endfor */
          DUDx_ave += DU*dX.x;
          DUDy_ave += DU*dX.y;
          DxDx_ave += dX.x*dX.x;
          DxDy_ave += dX.x*dX.y;
          DyDy_ave += dX.y*dX.y;
      } /* endfor */
  					    
      DUDx_ave = DUDx_ave/double(n_pts);
      DUDy_ave = DUDy_ave/double(n_pts);
      DxDx_ave = DxDx_ave/double(n_pts);
      DxDy_ave = DxDy_ave/double(n_pts);
      DyDy_ave = DyDy_ave/double(n_pts);

      SolnBlk->dUdx[i][j] = (DUDx_ave*DyDy_ave-DUDy_ave*DxDy_ave)/
                               (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
      SolnBlk->dUdy[i][j] = (DUDy_ave*DxDx_ave-DUDx_ave*DxDy_ave)/
                               (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
  
      // Calculate slope limiters. 
      if (!SolnBlk->Freeze_Limiter) {
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
    
            dX = SolnBlk->Grid.xfaceE(i, j)-SolnBlk->Grid.Cell[i][j].Xc;
            uQuad[0] = u0 + 
                       SolnBlk->dUdx[i][j][n]*dX.x +
                       SolnBlk->dUdy[i][j][n]*dX.y ;
            dX = SolnBlk->Grid.xfaceW(i, j)-SolnBlk->Grid.Cell[i][j].Xc;
            uQuad[1] = u0 + 
                       SolnBlk->dUdx[i][j][n]*dX.x +
                       SolnBlk->dUdy[i][j][n]*dX.y ;
            dX = SolnBlk->Grid.xfaceN(i, j)-SolnBlk->Grid.Cell[i][j].Xc;
            uQuad[2] = u0 + 
                       SolnBlk->dUdx[i][j][n]*dX.x +
                       SolnBlk->dUdy[i][j][n]*dX.y ;
            dX = SolnBlk->Grid.xfaceS(i, j)-SolnBlk->Grid.Cell[i][j].Xc;
            uQuad[3] = u0 + 
                       SolnBlk->dUdx[i][j][n]*dX.x +
                       SolnBlk->dUdy[i][j][n]*dX.y ;
    
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
      SolnBlk->dUdx[i][j] = W_VACUUM;
      SolnBlk->dUdy[i][j] = W_VACUUM; 
      SolnBlk->phi[i][j]  = W_VACUUM;
  } /* endif */

}





/****************************************************************/
/**** PRECONDTIONER REQUIRED SPECIALIZATIONS & FUNCTIONS ********/
/****************************************************************/  



/*!**************************************************************
 * Specialization of Block_Preconditioner::Preconditioner_dSdU  *
 *                                                              *
 * Calculates the dFdU matrix used to generate the approximate  *               
 * Jacobian for the Block Preconditioner.                       *
 ****************************************************************/
template<> inline void Block_Preconditioner<Rte2D_State,
					    Rte2D_Quad_Block,					    
					    Rte2D_Input_Parameters>::
Preconditioner_dSdU(int ii, int jj, DenseMatrix &dRdU){
  //Source Terms 

 
  //Add Jacobian for axisymmetric source terms
  if (SolnBlk->Axisymmetric) {
    SolnBlk->Uo[ii][jj].dSadU(dRdU,SolnBlk->Sp[ii][jj],
			      SolnBlk->Axisymmetric);
  }  

  //Add Jacobian for regular source terms
  SolnBlk->Uo[ii][jj].dSdU(dRdU);


}


/*!**************************************************************
 * Specialization of Block_Preconditioner::                     *
 *                               normalize_Preconditioner_dFdU  *
 *                                                              *
 * Normaliazes the dFdU matrix used to generate the approximate *               
 * Jacobian for the Block Preconditioner.                       *
 ****************************************************************/
template<> inline void Block_Preconditioner<Rte2D_State,
					    Rte2D_Quad_Block,
					    Rte2D_Input_Parameters>::
normalize_Preconditioner_dFdU(DenseMatrix &dFdU) { /*EMPTY - NOT NORMALIZED*/ }





/*!**************************************************************
 * Specialization of Block_Preconditioner::                     *
 *                           First_Order_Inviscid_Jacobian_HLLE *
 *                                                              *
 * Calculate First Order Local Jacobian Block(s) Coresponding   *
 * to Cell(i,j) using regular upwind flux function.  We are just*
 * using this function as a quick fix to fit in with the group  *
 * code.                                                        *
 ****************************************************************/
template<> inline void Block_Preconditioner<Rte2D_State,
					    Rte2D_Quad_Block,
					    Rte2D_Input_Parameters>::
First_Order_Inviscid_Jacobian_HLLE(const int &cell_index_i,const int &cell_index_j, 
				   DenseMatrix* Jacobian){              
  
  //! Caculate normal vectors -> in Vector2D format. 
  Vector2D nface_N = SolnBlk->Grid.nfaceN(cell_index_i,cell_index_j-1);
  Vector2D nface_S = SolnBlk->Grid.nfaceS(cell_index_i,cell_index_j+1);
  Vector2D nface_E = SolnBlk->Grid.nfaceE(cell_index_i-1,cell_index_j);
  Vector2D nface_W = SolnBlk->Grid.nfaceW(cell_index_i+1,cell_index_j);

  //! Calculate dFdU using solutions in the rotated frame -> matrix in DenseMatrix format. 
  DenseMatrix dFdU_N(blocksize,blocksize,ZERO); 
  DenseMatrix dFdU_S(blocksize,blocksize,ZERO); 
  DenseMatrix dFdU_E(blocksize,blocksize,ZERO); 
  DenseMatrix dFdU_W(blocksize,blocksize,ZERO); 

  //! Calculate Jacobian matrix -> blocksizexblocksize matrix in DenseMatrix format
  //Solution Rotate provided in pState 
  dFndU( dFdU_N, SolnBlk->Uo[cell_index_i][cell_index_j], nface_N); 
  dFndU( dFdU_S, SolnBlk->Uo[cell_index_i][cell_index_j], nface_S);
  dFndU( dFdU_E, SolnBlk->Uo[cell_index_i][cell_index_j], nface_E);
  dFndU( dFdU_W, SolnBlk->Uo[cell_index_i][cell_index_j], nface_W);

  //! Calculate Jacobian matrix -> blocksizexblocksize matrix in DenseMatrix format
  //North
  Jacobian[NORTH] = SolnBlk->Grid.lfaceN(cell_index_i,cell_index_j-1) 
    * SolnBlk->SpN[cell_index_i][cell_index_j-1] * dFdU_N; 

  //South
  Jacobian[SOUTH] = SolnBlk->Grid.lfaceS(cell_index_i,cell_index_j+1) 
    * SolnBlk->SpS[cell_index_i][cell_index_j+1] * dFdU_S;

  //East
  Jacobian[EAST] = SolnBlk->Grid.lfaceE(cell_index_i-1,cell_index_j) 
    * SolnBlk->SpE[cell_index_i-1][cell_index_j] * dFdU_E;

  //West
  Jacobian[WEST] = SolnBlk->Grid.lfaceW(cell_index_i+1,cell_index_j) 
    * SolnBlk->SpW[cell_index_i+1][cell_index_j] * dFdU_W;

  //Center calculated from neighbours
  //! Using the fact that dF/dU(right) = - dF/dU(left) 
  Jacobian[CENTER] += (Jacobian[NORTH] + Jacobian[SOUTH] + Jacobian[EAST]  + Jacobian[WEST])
                     /(SolnBlk->Grid.Cell[cell_index_i][cell_index_j].A * SolnBlk->Sp[cell_index_i][cell_index_j]);

  Jacobian[NORTH] = -Jacobian[NORTH]/(SolnBlk->Grid.Cell[cell_index_i][cell_index_j-1].A*SolnBlk->Sp[cell_index_i][cell_index_j-1]);
  Jacobian[SOUTH] = -Jacobian[SOUTH]/(SolnBlk->Grid.Cell[cell_index_i][cell_index_j+1].A*SolnBlk->Sp[cell_index_i][cell_index_j+1]);
  Jacobian[EAST] = -Jacobian[EAST]/(SolnBlk->Grid.Cell[cell_index_i-1][cell_index_j].A*SolnBlk->Sp[cell_index_i-1][cell_index_j]);
  Jacobian[WEST] = -Jacobian[WEST]/(SolnBlk->Grid.Cell[cell_index_i+1][cell_index_j].A*SolnBlk->Sp[cell_index_i+1][cell_index_j]);

 
}




/*****************************************************************************
 * Specialization of Block_Preconditioner::                                  *
 *                           First_Order_Inviscid_Jacobian_Roe               *
 *                                                                           *
 *  Calculate First Order Local Jacobian Block(s) Coresponding to Cell(i,j)  *
 *  using Finite difference approximation.  The 'Roe' function is only used  *
 *  to fit in with the rest of the code.                                     *
 *****************************************************************************/
template<> inline void Block_Preconditioner<Rte2D_State,
					    Rte2D_Quad_Block,
					    Rte2D_Input_Parameters>::
First_Order_Inviscid_Jacobian_Roe(const int &cell_index_i,const int &cell_index_j, 
				   DenseMatrix* Jacobian){              
    
  //! Calculate Jacobian matrix -> blocksizexblocksize matrix in DenseMatrix format
  dFIdU_FD(Jacobian[NORTH],*SolnBlk,*Input_Parameters,cell_index_i,cell_index_j,NORTH);
  dFIdU_FD(Jacobian[SOUTH],*SolnBlk,*Input_Parameters,cell_index_i,cell_index_j,SOUTH); 
  dFIdU_FD(Jacobian[EAST],*SolnBlk,*Input_Parameters,cell_index_i,cell_index_j,EAST);        
  dFIdU_FD(Jacobian[WEST],*SolnBlk,*Input_Parameters,cell_index_i,cell_index_j,WEST); 

  //Center calculated from neighbours
  //! Using the fact that dF/dU(right) = - dF/dU(left) 
  Jacobian[CENTER] += (Jacobian[NORTH] + Jacobian[SOUTH] + Jacobian[EAST]  + Jacobian[WEST])
    /(SolnBlk->Grid.Cell[cell_index_i][cell_index_j].A * SolnBlk->Sp[cell_index_i][cell_index_j]);

  Jacobian[NORTH] = -Jacobian[NORTH]/(SolnBlk->Grid.Cell[cell_index_i][cell_index_j-1].A*SolnBlk->Sp[cell_index_i][cell_index_j-1]);
  Jacobian[SOUTH] = -Jacobian[SOUTH]/(SolnBlk->Grid.Cell[cell_index_i][cell_index_j+1].A*SolnBlk->Sp[cell_index_i][cell_index_j+1]);
  Jacobian[EAST] = -Jacobian[EAST]/(SolnBlk->Grid.Cell[cell_index_i-1][cell_index_j].A*SolnBlk->Sp[cell_index_i-1][cell_index_j]);
  Jacobian[WEST] = -Jacobian[WEST]/(SolnBlk->Grid.Cell[cell_index_i+1][cell_index_j].A*SolnBlk->Sp[cell_index_i+1][cell_index_j]);

}


/********************************************************
 * Routine: Flux Jacobian                               *
 *                                                      *
 *     Finite difference approximation                  *
 *                                                      *
 ********************************************************/
void dFIdU_FD(DenseMatrix& dRdU, Rte2D_Quad_Block &SolnBlk,  
	      Rte2D_Input_Parameters &Input_Parameters,
	      const int ii, const int jj, const int Orient)
{

  int NUM_VAR_RTE2D = dRdU.get_n();   
  DenseMatrix dFidU(NUM_VAR_RTE2D, NUM_VAR_RTE2D,ZERO);
  
  int Ri, Rj;
  Vector2D nface;
  double lface;
  Rte2D_State Ul, Ur;   
  Rte2D_State UA, UB;   
  Rte2D_State FluxA, FluxB;   
  Ul.Zero();  Ur.Zero();
  double perturb = 5e-6;

  
  switch (Orient) {

    // NORTH
  case NORTH:
    Ri = ii; Rj=jj-1;
    nface = SolnBlk.Grid.nfaceN(Ri, Rj);
    lface = SolnBlk.Grid.lfaceN(Ri, Rj)*SolnBlk.SpN[Ri][Rj];
    Ur = SolnBlk.Uo[ii][jj];
    Ul = SolnBlk.Uo[Ri][Rj];
    break;

    // SOUTH
  case SOUTH:
    Ri = ii; Rj=jj+1;
    nface = SolnBlk.Grid.nfaceS(Ri, Rj);
    lface = SolnBlk.Grid.lfaceS(Ri, Rj)*SolnBlk.SpS[Ri][Rj];
    Ur = SolnBlk.Uo[ii][jj];
    Ul = SolnBlk.Uo[Ri][Rj];
    break;

    // EAST
  case EAST:
    Ri = ii-1; Rj=jj;
    nface = SolnBlk.Grid.nfaceE(Ri, Rj);     
    lface = SolnBlk.Grid.lfaceE(Ri, Rj)*SolnBlk.SpE[Ri][Rj];
    Ur = SolnBlk.Uo[ii][jj];
    Ul = SolnBlk.Uo[Ri][Rj];
    break;

    // WEST
  case WEST:
    Ri = ii+1; Rj=jj;
    nface = SolnBlk.Grid.nfaceW(Ri, Rj);
    lface = SolnBlk.Grid.lfaceW(Ri, Rj)*SolnBlk.SpW[Ri][Rj];
    Ur = SolnBlk.Uo[ii][jj];
    Ul = SolnBlk.Uo[Ri][Rj];
    break;

  } /* endswitch */ 


  // compute the derivatives
  for(int jcol=0; jcol<(NUM_VAR_RTE2D); jcol++){
    UA = Ur;
    UB = Ur;
    UA[jcol+1] += perturb; 	 
    UB[jcol+1] -= perturb; 
    
    FluxA = Flux_n(Ul,UA, nface);       
    FluxB = Flux_n(Ul,UB, nface);
    
    for(int irow=0; irow<(NUM_VAR_RTE2D); irow++){
      dFidU(irow,jcol) = (FluxA[irow+1] - FluxB[irow+1])/(TWO*perturb);           
    }
  } 

  // copy over

  dRdU += lface*dFidU;

}


#endif // _RTE2D_NKS_INCLUDED 

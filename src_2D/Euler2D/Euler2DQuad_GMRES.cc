/* Include 2D Euler quadrilateral mesh GMRES header file. */

#ifndef _EULER2D_QUAD_GMRES_INCLUDED
#include "Euler2DQuad_GMRES.h"
#endif // _EULER2D_QUAD_GMRES_INCLUDED

#define ABS(x) ((x)<0 ? (-(x)) : (x))

void GeneratePlaneRotation(double &dx, double &dy, double &cs, double &sn) {
  if (dy == 0.0) {
    cs = 1.0;
    sn = 0.0;
  } else if (ABS(dy) > ABS(dx)) {
    double temp = dx / dy;
    sn = 1.0 / sqrt( 1.0 + temp*temp );
    cs = temp * sn;
  } else {
    double temp = dy / dx;
    cs = 1.0 / sqrt( 1.0 + temp*temp );
    sn = temp * cs;
  }
}

void ApplyPlaneRotation(double &dx, double &dy, double &cs, double &sn) {
  double temp  =  cs * dx + sn * dy;
  dy = -sn * dx + cs * dy;
  dx = temp;
}

GMRES_Block_RightPrecon_MatrixFree::GMRES_Block_RightPrecon_MatrixFree(int dim_, 
                                                                       int max_iter_, 
                                                                       double tol_) {
    // override defaults and resource values
    dim = dim_;
    max_iter = max_iter_;
    tol = tol_;
}

/****************************************************************************
 * GMRES_Block_RightPrecon_MatrixFree::solve -- apply right-preconditioned  *
 *                                              matrix-free GMRES.          *
 ****************************************************************************/
void GMRES_Block_RightPrecon_MatrixFree::solve(Euler2D_Quad_Block *Soln_ptr, 
					       AdaptiveBlock2D_List &Soln_Block_List, 
                                               Euler2D_Input_Parameters &Input_Parameters,
					       GMRES_Block *G,
					       BpPrecon *M,
					       const int &normalization,
                                               const int Number_of_Newton_Steps) {   
  int m1 = dim+1; // used inside H macro
  int i, j, k, i_first_time_through;
  int Bcount, NBLK = Soln_Block_List.Nblk;
  int xcount, ycount, num, temp_num;
  int index, Grid_i, Grid_j;
  int error_flag;
  
  double beta, resid0;
  double  epsilon;
  
  double total_norm_H = 0.0;
  double total_norm_z = 0.0;
  double total_norm_x = 0.0;
  
  integer inc = 1;  // vector stride is always 1
  doublereal temp;
  int Number_of_GMRES_Iterations = 0;
  iter = 0;

  /* Casting preconditioner pointer. */
  BILUK   *ILUK_Precon;
  BJacobi *Jacobi_Precon;
  if (G[0].P_Switch == 1) {
    ILUK_Precon  = (BILUK*)M;
  } else if (G[0].P_Switch == 2) {
    Jacobi_Precon = (BJacobi*)M;
  }
  
  integer * N = new integer[NBLK];
  
  for (Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
    if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	N[Bcount] = G[Bcount].scalar_dim;
    } /* endif */
  } /* endfor */
  
  for (Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
    if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
      
      if (normalization) normalize_RHS(G[Bcount].b, Soln_ptr[Bcount]);
      
      /* Require R(Uo), which is negative of the RHS stored in b vector. */
      for (num = 0; num < N[Bcount]; ++num)
	G[Bcount].Q2[num] = -1.0 * G[Bcount].b[num];
    
    } /* endif */
  } /* endfor */

  i_first_time_through = ON;
  do {
    // Compute initial residual and its norm
    /* Calculate global epsilon base on 2-norm of x. */
    total_norm_x = ZERO;
    for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
      if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	total_norm_x = total_norm_x + sqr(F77NAME(dnrm2)
					  (&N[Bcount], G[Bcount].x,&inc));
	
      } /* endif */
    } /* endfor */
    total_norm_x = sqrt(CFDkit_Summation_MPI(total_norm_x));
    
    if (!i_first_time_through) {
      epsilon = 1.0e-7/sqrt(total_norm_x);
//       cout << "\n Number_of_GMRES_Iterations = " << Number_of_GMRES_Iterations << " total_norm_x = " 
//            << total_norm_x << " epsilon = " << epsilon;
    } /* endif */
 
    if (i_first_time_through) { // FIRST TIME THROUGH - NO RESTART NEEDED
      for (Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
         if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
  	    G[Bcount].tmp = &(G[Bcount].V[(0)*N[Bcount]]);
	    for (num = 0; num < N[Bcount]; ++num) G[Bcount].tmp[num] = 0.0;
         } /* endif */
      } /* endfor */
    } else { // NOT FIRST TIME THROUGH - RESTART APPLIED
       /* BEGIN MATRIX-FREE FOR RESTART ************************************/
       for (Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
         if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {

 	   if (CFDkit_Primary_MPI_Processor() &&  Bcount == 0) { 
	     cout << "  **** RESTART ****" << endl; cout.flush();
	   } /* endif */
	  
	   // Q = 0.0;
	   for (num = 0; num < N[Bcount]; ++num) G[Bcount].Q[num] = 0.0;  
	   // Q = Q + epsilon * x;
	   F77NAME(daxpy)(&N[Bcount],&epsilon,G[Bcount].x,&inc,G[Bcount].Q,&inc);  
	  
	   if (normalization) denormalize_solution(G[Bcount].Q, Soln_ptr[Bcount]);
	  
	   /* Copy Q vector to Soln_ptr.U[i][j] (conserved variables)
	      to evaluate the residual. */
	   for (ycount = Soln_ptr[Bcount].JCl-2; 
	        ycount <= Soln_ptr[Bcount].JCu+2; ++ycount) {
	     for (xcount = Soln_ptr[Bcount].ICl-2; 
		  xcount <= Soln_ptr[Bcount].ICu+2; ++xcount) {
	      
	       Soln_ptr[Bcount].U[xcount][ycount].d    = 
		 G[Bcount].Q[G[Bcount].index(xcount,ycount)]
		 + G[Bcount].U[G[Bcount].index(xcount,ycount)];
	      
	       Soln_ptr[Bcount].U[xcount][ycount].dv.x = 
		 G[Bcount].Q[G[Bcount].index(xcount,ycount,1)]
		 + G[Bcount].U[G[Bcount].index(xcount,ycount,1)];
	      
	       Soln_ptr[Bcount].U[xcount][ycount]. dv.y = 
		 G[Bcount].Q[G[Bcount].index(xcount,ycount,2)]
		 + G[Bcount].U[G[Bcount].index(xcount,ycount,2)];
	      
	       Soln_ptr[Bcount].U[xcount][ycount].E    = 
		 G[Bcount].Q[G[Bcount].index(xcount,ycount,3)]
		 + G[Bcount].U[G[Bcount].index(xcount,ycount,3)];
	      
	       /* Update primitive variables. */
	       Soln_ptr[Bcount].W[xcount][ycount] = 
		 W(Soln_ptr[Bcount].U[xcount][ycount]);
	      
	     } /* endfor */
	   } /* endfor */

	   dUdt_Residual_Evaluation(Soln_ptr[Bcount],
		 		    Input_Parameters);
	  
         } /* endif */
       } /* endfor */

       // Send boundary flux corrections at block interfaces with resolution changes.
       error_flag = Send_Conservative_Flux_Corrections(Soln_ptr,
		       			               Soln_Block_List,
						       NUM_VAR_EULER2D);
       if (error_flag) {
	 cout << "\n Euler2D GMRES ERROR: Euler2D flux correction message passing error on processor "
	      << Soln_Block_List.ThisCPU
	      << ".\n";
	 cout.flush();
       } /* endif */
       error_flag = CFDkit_OR_MPI(error_flag);
	  
       // Apply boundary flux corrections to residual to ensure that method is conservative.
       Apply_Boundary_Flux_Corrections(Soln_ptr,
				       Soln_Block_List);

       for (Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
         if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
 	   G[Bcount].copy_residual(G[Bcount].Q1, Soln_ptr[Bcount]);
	  
	   if (normalization) normalize_RHS(G[Bcount].Q1, Soln_ptr[Bcount]);
	  
// 	   for (num = 0; num < N[Bcount]; ++num) G[Bcount].dt_vec[num] = ZERO; 	  
// 	   if (normalization) G[Bcount].dTime = G[Bcount].dTime * Euler2D_W_STDATM.a();

           /* Apply Finite Time Step. */
	   for (num = 0; num < N[Bcount]; ++num) { 
	       temp_num = num/G[Bcount].blocksize;
	       Grid_i = temp_num%G[Bcount].xpts;
	       Grid_j = int(temp_num/G[Bcount].xpts);

 	       if (normalization) {
                  G[Bcount].dt_vec[num] = G[Bcount].x[num] / 
                                          (Soln_ptr[Bcount].dt[Grid_i][Grid_j] * Euler2D_W_STDATM.a());
               } else {
                  G[Bcount].dt_vec[num] = G[Bcount].x[num] / 
                                          Soln_ptr[Bcount].dt[Grid_i][Grid_j];
               } /* endif */

           } /* endfor */
	  
	   G[Bcount].tmp = &(G[Bcount].V[(0)*N[Bcount]]);
	  
	   for (num = 0; num < N[Bcount]; ++num) {
	     if (num%G[Bcount].blocksize == 0) {
	       temp_num = num/G[Bcount].blocksize;
	       Grid_i = temp_num%G[Bcount].xpts;
	       Grid_j = int(temp_num/G[Bcount].xpts);
	      
	       if (Grid_j < G[Bcount].Nghost-G[Bcount].overlap || 
		   Grid_j > Soln_ptr[Bcount].JCu+G[Bcount].overlap  || 
		   Grid_i < G[Bcount].Nghost-G[Bcount].overlap  ||  
		   Grid_i > Soln_ptr[Bcount].ICu+G[Bcount].overlap ) {
		  /************** GHOST CELLS *************/
		  // V(i+1) = (Q1 - Q2) / epsilon;
		  for (j=num; j<num+G[Bcount].blocksize; j++) {
		    G[Bcount].tmp[j] = 
		      (G[Bcount].Q1[j]-G[Bcount].Q2[j]) / epsilon;
		  } /* endfor */
	       } else {
		  /************ INTERIOR CELLS ************/ 
		  // V(i+1) = (Q1 - Q2) / epsilon - z / h;
		  for (j=num; j<num+G[Bcount].blocksize; j++) {
		    G[Bcount].tmp[j] = 
		      (G[Bcount].Q1[j]-G[Bcount].Q2[j]) / epsilon 
		      - G[Bcount].dt_vec[j];
		  } /* endfor */
	       } /* endif */ 
	     } /* endif */ 
	   } /* endfor */
	  
         } /* endif */
       } /* endfor */
    
       /* END OF MATRIX-FREE ***********************************/

    } /* endif */ // NOT FIRST TIME THROUGH - RESTART APPLIED

    // CALCULATE NORM OF FIRST SEARCH VECTOR, V(0).
    beta = 0;
    for (Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
      if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	/* V(0) = V(0) - b */
	temp = - 1.0;
	F77NAME(daxpy)(&N[Bcount], &temp, G[Bcount].b, &inc, 
		       &(G[Bcount].V[(0)*N[Bcount]]), &inc);
	/* beta = norm(V(0)) */
	beta = beta + sqr(F77NAME(dnrm2)(&N[Bcount], &(G[Bcount].V[(0)*N[Bcount]]), 
					 &inc));
      }
    }
    beta = sqrt(CFDkit_Summation_MPI(beta));
    
    // RESSCALE V(0) USING NORM.
    for (Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
      if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	
	/* V(0) = -V(0)/beta */
	temp = -1.0/beta;
	F77NAME(dscal)(&N[Bcount], &temp, &(G[Bcount].V[(0)*N[Bcount]]), &inc);
	
	/* save the very first residual norm */
	if (Number_of_GMRES_Iterations == 0) {
	  resid0 = beta;
	}
	
	for (i = 1; i < dim+1; i++)
	  G[Bcount].s[i] = 0.0;
	
	G[Bcount].s[0] = beta;
	
      } /* endif */
    } /* endfor */
    
    i = -1;
    //--------------------------//
    // Begin Primary GMRES Loop //
    //--------------------------//
    do {
      
      i++;
      Number_of_GMRES_Iterations++;
      iter = Number_of_GMRES_Iterations;
      
      // Calculate z vector using preconditioner.
      for (Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
	if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	  
	  /* Set vector switch and search direction counter. */
	  G[Bcount].vector_switch     = 1;
	  G[Bcount].search_directions = i;
	  
	  /* z = Minv * V(i) -> stored in W(i). */
	  if (G[0].P_Switch == 1) {
	    ILUK_Precon[Bcount].apply(N[Bcount], 1, &(G[Bcount].V[(i)*N[Bcount]]), 
				      N[Bcount], &(G[Bcount].W[(i)*N[Bcount]]), 
				      N[Bcount]);
	  } else if (G[0].P_Switch == 2) {
	    Jacobi_Precon[Bcount].apply(N[Bcount], 1, &(G[Bcount].V[(i)*N[Bcount]]), 
					N[Bcount], &(G[Bcount].W[(i)*N[Bcount]]), 
					N[Bcount]);
	  } /* endif */
	  
	} /* endif */
      } /* endfor */   
      
      ///////////////////////////////////////////////////////////
      ////////////////* BEGIN MESSAGE PASSING *//////////////////
      ///////////////////////////////////////////////////////////
      
      /* MPI barrier to ensure processor synchronization. */
      CFDkit_Barrier_MPI();  
      
      /* Send solution information between neighbouring blocks.*/
      error_flag = Send_All_Messages(G, 
				     Soln_Block_List,
				     NUM_VAR_EULER2D, 
				     OFF);
      if (error_flag) {
	cout << "\n Euler2D ERROR: Euler2D message passing error on processor "
	     << Soln_Block_List.ThisCPU
	     << ".\n";
	cout.flush();
      } /* endif */
      error_flag = CFDkit_OR_MPI(error_flag);
      
      /////////////////////////////////////////////////////////////
      //////////////*  END OF MESSAGE PASSING  *///////////////////
      /////////////////////////////////////////////////////////////  
      
      for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) { 
	if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {   
	  
	  // Apply boundary conditions on z-vector;
	  GMRES_BCs(Soln_ptr[Bcount], &G[Bcount].W[(i)*N[Bcount]], G[Bcount], normalization);      
	  
	  /* FINITE TIME STEP VECTOR *******************************/
	  /* dtime = 0.0 */
	  for (j=0;j<N[Bcount];j++) G[Bcount].dt_vec[j] = 0.0;
	  
	  /* dtime = dtime + z */
	  temp = 1.0;
	  F77NAME(daxpy)(&N[Bcount],&temp,&(G[Bcount].W[(i)*N[Bcount]]),
			 &inc,G[Bcount].dt_vec,&inc);
	  
// 	  if (normalization) G[Bcount].dTime = G[Bcount].dTime * Euler2D_W_STDATM.a();
	  
          /* Apply Finite Time Step. */
 	  /* dtime = dtime / h */
	  for (num = 0; num < N[Bcount]; ++num) {
	      temp_num = num/G[Bcount].blocksize;
	      Grid_i = temp_num%G[Bcount].xpts;
	      Grid_j = int(temp_num/G[Bcount].xpts);

	      if (normalization) {
                 G[Bcount].dt_vec[num] = G[Bcount].dt_vec[num] / 
                                         (Soln_ptr[Bcount].dt[Grid_i][Grid_j]*Euler2D_W_STDATM.a());
              } else {
                 G[Bcount].dt_vec[num] = G[Bcount].dt_vec[num] / 
                                         Soln_ptr[Bcount].dt[Grid_i][Grid_j];
              } /* endif */
          } /* endfor */
	  /* END OF FINITE TIME STEP VECTOR ***********************/
	  
	} /* endif */
      } /* endfor */
      
      /* Calculate global epsilon base on 2-norm of z. */
      total_norm_z = 0.0;
      for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
	if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	  total_norm_z = total_norm_z + sqr(F77NAME(dnrm2)
					    (&N[Bcount], 
					     &(G[Bcount].W[(i)*N[Bcount]]), 
					     &inc));
	  
	} /* endif */
      } /* endfor */
      total_norm_z = sqrt(CFDkit_Summation_MPI(total_norm_z));     	    
      
      epsilon = 1.0e-07/sqrt(total_norm_z);
//        cout << "\n Number_of_GMRES_Iterations = " << Number_of_GMRES_Iterations << " total_norm_z = " 
//             << total_norm_z << " epsilon = " << epsilon;     

      /* BEGIN MATRIX-FREE FOR PRIMARY GMRES LOOP ************************************/

      for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) { 
	if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) { 
	  
	  /* Q = 0.0 */
	  for (j=0 ; j < N[Bcount] ; j++) G[Bcount].Q[j] = ZERO;
	  
	  /* Q = Q + epsilon * z = epsilon * W(i) */
	  F77NAME(daxpy)(&N[Bcount],&epsilon,&(G[Bcount].W[(i)*N[Bcount]]),
			 &inc,G[Bcount].Q,&inc);
	  
	  if (normalization) denormalize_solution(G[Bcount].Q, Soln_ptr[Bcount]);
	  
	  /* Copy Q vector to Soln_ptr.U[i][j] (conserved variables)
	     to evaluate the residual. */
	  for (ycount = Soln_ptr[Bcount].JCl-2; 
	       ycount <= Soln_ptr[Bcount].JCu+2; ++ycount) {
	    for (xcount = Soln_ptr[Bcount].ICl-2; 
		 xcount <= Soln_ptr[Bcount].ICu+2; ++xcount) {
	      
	      Soln_ptr[Bcount].U[xcount][ycount].d    = 
		G[Bcount].U[G[Bcount].index(xcount,ycount)]
		+ G[Bcount].Q[G[Bcount].index(xcount,ycount)];
	      
	      Soln_ptr[Bcount].U[xcount][ycount].dv.x = 
		G[Bcount].U[G[Bcount].index(xcount,ycount,1)]
		+ G[Bcount].Q[G[Bcount].index(xcount,ycount,1)];
	      
	      Soln_ptr[Bcount].U[xcount][ycount].dv.y = 
		G[Bcount].U[G[Bcount].index(xcount,ycount,2)]
		+ G[Bcount].Q[G[Bcount].index(xcount,ycount,2)];
	      
	      Soln_ptr[Bcount].U[xcount][ycount].E    = 
		G[Bcount].U[G[Bcount].index(xcount,ycount,3)]
		+ G[Bcount].Q[G[Bcount].index(xcount,ycount,3)];
	      
	      /* Update primitive variables. */
	      Soln_ptr[Bcount].W[xcount][ycount] = 
		W(Soln_ptr[Bcount].U[xcount][ycount]);
		
	    } /* endfor */
	  } /* endfor */
	  
	  dUdt_Residual_Evaluation(Soln_ptr[Bcount],
				   Input_Parameters);
	  
	 } /* endif */
      } /* endfor */

      // Send boundary flux corrections at block interfaces with resolution changes.
      error_flag = Send_Conservative_Flux_Corrections(Soln_ptr,
	              			              Soln_Block_List,
						      NUM_VAR_EULER2D);
      if (error_flag) {
	cout << "\n Euler2D GMRES ERROR: Euler2D flux correction message passing error on processor "
	     << Soln_Block_List.ThisCPU
	     << ".\n";
        cout.flush();
      } /* endif */
      error_flag = CFDkit_OR_MPI(error_flag);
	  
      // Apply boundary flux corrections to residual to ensure that method is conservative.
      Apply_Boundary_Flux_Corrections(Soln_ptr,
	       		              Soln_Block_List);

      for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) { 
	if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) { 
	  G[Bcount].copy_residual(G[Bcount].Q1, Soln_ptr[Bcount]);
	  
	  if (normalization) normalize_RHS(G[Bcount].Q1, Soln_ptr[Bcount]);
	  
	  /* tmp = V(i+1) */
	  G[Bcount].tmp = &(G[Bcount].V[(i+1)*N[Bcount]]);
	  
	  for (num = 0; num < N[Bcount]; ++num) {
	    if (num%G[Bcount].blocksize == 0) {
	      temp_num = num/G[Bcount].blocksize;
	      Grid_i = temp_num%G[Bcount].xpts;
	      Grid_j = int(temp_num/G[Bcount].xpts);
	      
	      if (Grid_j < G[Bcount].Nghost-G[Bcount].overlap || 
		  Grid_j > Soln_ptr[Bcount].JCu+G[Bcount].overlap  || 
		  Grid_i < G[Bcount].Nghost-G[Bcount].overlap  ||  
		  Grid_i > Soln_ptr[Bcount].ICu+G[Bcount].overlap ) 
		{ 
		  /************** GHOST CELLS *************/
		  /* V(i+1) = (Q1 - Q2) / epsilon */
		  for (j=num; j<num+G[Bcount].blocksize; j++) {
		    G[Bcount].tmp[j] = 
		      (G[Bcount].Q1[j]-G[Bcount].Q2[j]) / epsilon; 
		    
		  }
		} else {
		  /************ INTERIOR CELLS ************/ 
		  /* V(i+1) = (Q1 - Q2) / epsilon - z / h */
		  for (j=num; j<num+G[Bcount].blocksize; j++) {
		    G[Bcount].tmp[j] = 
		      (G[Bcount].Q1[j]-G[Bcount].Q2[j]) / epsilon 
		      - G[Bcount].dt_vec[j];
		    
		  } /* endfor */ 
		} /* endif */ 
	    } /* endif */ 
	  } /* endfor */
	  
	} /* endif */
      } /* endfor */
      
      /* END OF MATRIX-FREE ***********************************/

      for (k = 0; k <= i; k++) {
	total_norm_H = 0.0;
	for (Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
	  if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED){
	    
	    /* total_H = dot(W,V(k)) */
	    total_norm_H = total_norm_H + F77NAME(ddot) 
	      (&N[Bcount], 
	       &(G[Bcount].V[(i+1)*N[Bcount]]),
	       &inc, 
	       &(G[Bcount].V[(k)*N[Bcount]]), 
	       &inc);
	  }
	}
	total_norm_H = CFDkit_Summation_MPI(total_norm_H);
	
	for (Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
	  if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED){
	    
	    /* V(i+1) -= H(k, i) * V(k) */
	    G[Bcount].H[(i)*m1+(k)] = total_norm_H;
	    temp = -1.0 * G[Bcount].H[(i)*m1+(k)];
	    F77NAME(daxpy)(&N[Bcount], &temp, &(G[Bcount].V[(k)*N[Bcount]]), &inc,
			   &(G[Bcount].V[(i+1)*N[Bcount]]), &inc);
	  } /* endif */
	} /* endfor */
      } /* endfor */
      
      total_norm_H = 0.0;
      /* Calculate 2-norm of V(i+1) -> H(i+1,1) */
      for (Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
	if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED){
	  
	  total_norm_H = total_norm_H + sqr(F77NAME(dnrm2)
					    (&N[Bcount], &(G[Bcount].V[(i+1)*N[Bcount]]),
					     &inc));
	} /* endif */
      } /* endfor */
      total_norm_H = sqrt(CFDkit_Summation_MPI(total_norm_H)); 
      
      for (Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
	if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED){
	  
	  G[Bcount].H[(i)*m1+(i+1)] = total_norm_H;
	  temp = 1.0 / G[Bcount].H[(i)*m1+(i+1)];
	  
	  /* V(i+1) = V(i+1) / H(i+1, i) */
	  F77NAME(dscal)(&N[Bcount], &temp, &(G[Bcount].V[(i+1)*N[Bcount]]), &inc);
	  
	  } /* endif */
      } /* endfor */
      
      for (Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
	if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED){
	  
	  for (k = 0; k < i; k++) {
	    ApplyPlaneRotation(G[Bcount].H[(i)*m1+(k)], 
			       G[Bcount].H[(i)*m1+(k+1)], 
			       G[Bcount].cs[k], G[Bcount].sn[k]);
	  } /* endfor */
	  
	  GeneratePlaneRotation(G[Bcount].H[(i)*m1+(i)],
				G[Bcount].H[(i)*m1+(i+1)], 
				G[Bcount].cs[i], G[Bcount].sn[i]);
	  
	  
	  ApplyPlaneRotation(G[Bcount].H[(i)*m1+(i)],
			     G[Bcount].H[(i)*m1+(i+1)], 
			     G[Bcount].cs[i], G[Bcount].sn[i]);
	  
	  
	  ApplyPlaneRotation(G[Bcount].s[i], G[Bcount].s[i+1], 
			     G[Bcount].cs[i], G[Bcount].sn[i]);
	  
	} /* endif */
      } /* endfor */
      
      rel_resid = ABS(G[0].s[i+1]) / resid0;
      cout << "\n Number_of_GMRES_Iterations = " << Number_of_GMRES_Iterations << " resid0 = " 
           << resid0 << " resid = " << ABS(G[0].s[i+1]) << " relative_residual = " << rel_resid;
      if (rel_resid <= tol) break;
      
    } while (i+1 < dim && Number_of_GMRES_Iterations+1 <= max_iter);
    
    //------------------//
    // UPDATE SOLUTION. // 
    //------------------//
    
    for (Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
      if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED){  
	
	/* Set vector switch. */
	G[Bcount].vector_switch = 0;
	
	/* solve upper triangular system in place */	  
	for (j = i; j >= 0; j--) { 
	  G[Bcount].s[j] /= G[Bcount].H[(j)*m1+(j)];
	  for (k = j-1; k >= 0; k--)
	    G[Bcount].s[k] -= G[Bcount].H[(j)*m1+(k)] * G[Bcount].s[j];
	}
	
	/* update solution, x */ 
	for (j = 0; j <= i; j++) {
	  /* x = x + s[j] * W(j) */
	  F77NAME(daxpy)(&N[Bcount], &G[Bcount].s[j], &(G[Bcount].W[(j)*N[Bcount]]), 
			 &inc, G[Bcount].x, &inc); 
	} /* endfor */  
      } /* endif */
    } /* endfor */   
    
    ///////////////////////////////////////////////////////////
    ////////////////* BEGIN MESSAGE PASSING *//////////////////
    ///////////////////////////////////////////////////////////
    
    /* MPI barrier to ensure processor synchronization. */
    CFDkit_Barrier_MPI();  
    
    /* Send solution information between neighbouring blocks.*/
    
    error_flag = Send_All_Messages(G, 
				   Soln_Block_List,
				   NUM_VAR_EULER2D, 
				   OFF);
    if (error_flag) {
      cout << "\n Euler2D ERROR: Euler2D message passing error on processor "
	   << Soln_Block_List.ThisCPU
	   << ".\n";
      cout.flush();
    } /* endif */
    error_flag = CFDkit_OR_MPI(error_flag);
    
    /////////////////////////////////////////////////////////////
    //////////////*  END OF MESSAGE PASSING  *///////////////////
    ///////////////////////////////////////////////////////////// 
    
    // Apply boundary conditions on x-vector;
    for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
      if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	GMRES_BCs(Soln_ptr[Bcount], G[Bcount].x, G[Bcount], normalization);
      } /* endif */
    } /* endfor */       
    
    i_first_time_through = OFF;

  } while (rel_resid > tol && Number_of_GMRES_Iterations+1 <= max_iter);
  
  delete [] N;
  
} /* End of GMRES_Block_RightPrecon_MatrixFree::solve. */ 

/***********************************************************************
 * Iterative routine -- GMRES_Algorithm                                *
 *                      (Right-preconditioned matrix-free)             *
 *                                                                     * 
 *  GMRES solves the unsymmetric linear system Ax = b using the        *
 *  Generalized Minimum Residual method                                * 
 *                                                                     *
 *         GMRES_Algorithm(Soln_ptr,                                   *
 *                         Soln_Block_List,                            *
 *                         Input_Parameters,                           *
 *       		   G, P_Switch, restart, gmrestol);            *
 *                                                                     *
 * List of required inputs:                                            *
 *                                                                     *
 *         Soln_ptr --  solution block pointer                         *
 *  Soln_Block_List -- local solution block list                       *
 * Input_Parameters --  parameters describing solution                 *
 *                G --  GMRES object containing all function variables *
 *         P_Switch --  preconditioner switch                          *
 *          restart --  number of iterations for each restart          *    
 *         gmrestol --  the residual after the final iteration         *
 *                                                                     *
 * Note:  GMRES follows the algorithm described on p. 22 of the        *
 *        Block Preconditioning Toolkit (BPKIT) reference manual.      *
 *                                                                     *
 ***********************************************************************/
void GMRES_Algorithm(Euler2D_Quad_Block *Soln_ptr, 
		     AdaptiveBlock2D_List &Soln_Block_List, 
                     Euler2D_Input_Parameters &Input_Parameters,
		     GMRES_Block *G, double &gmrestol, 
		     BILUK *MBILUK, BJacobi *MBJacobi,
		     const int normalization,
                     const int Number_of_Newton_Steps)
{
  int i, j;
  int Bcount;
  int NBLK = Soln_Block_List.Nblk;
  
  for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
    if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
      /* Initialize all GMRES variables except b and U. */
      for (i=0;i<G[Bcount].m+1;++i)  G[Bcount].s[i] = 0.0;
      for (i=0;i<G[Bcount].m;++i)   G[Bcount].cs[i] = 0.0;
      for (i=0;i<G[Bcount].m;++i)   G[Bcount].sn[i] = 0.0;
     
      for (i=0;i<G[Bcount].scalar_dim;++i)  G[Bcount].Q1[i]     = 0.0;
      for (i=0;i<G[Bcount].scalar_dim;++i)  G[Bcount].Q2[i]     = 0.0;
      for (i=0;i<G[Bcount].scalar_dim;++i)  G[Bcount].dt_vec[i] = 0.0;
      for (i=0;i<G[Bcount].scalar_dim;++i)  G[Bcount].x[i]      = 0.0;
      for (i=0;i<G[Bcount].scalar_dim;++i)  G[Bcount].Q[i]      = 0.0;
      
      for (i=0;i<(G[Bcount].m*(G[Bcount].m+1));++i)  G[Bcount].H[i] = 0.0;
      for (i=0;i<G[Bcount].m* G[Bcount].scalar_dim;++i)    
	G[Bcount].W[i] = 0.0;
      for (i=0;i<((G[Bcount].m + 1) * G[Bcount].scalar_dim);++i) 
	G[Bcount].V[i] = 0.0;

    } /* endif */
  } /* endfor */

  GMRES_Block_RightPrecon_MatrixFree GMRES_(G[0].m, 
	        			    G[0].max_gmres_iter, 
					    gmrestol);

  if (G[0].P_Switch == 1) {
    GMRES_.solve(Soln_ptr, Soln_Block_List,
                 Input_Parameters, 
		 G, MBILUK, normalization,
                 Number_of_Newton_Steps);
  } else if (G[0].P_Switch == 2) { 
    GMRES_.solve(Soln_ptr, Soln_Block_List, 
		 Input_Parameters,
		 G, MBJacobi, normalization,
                 Number_of_Newton_Steps);
  }

  if (normalization) {
    for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
      if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	denormalize_solution(G[Bcount].x, Soln_ptr[Bcount]);
      } /* endif */
    } /* endfor */
  } /* endif */

  if (GMRES_.get_iter() == G[0].max_gmres_iter) {
    if (CFDkit_Primary_MPI_Processor()) {  
      cout << "\n   Euler2D - GMRES ERROR: Unable to reach the specified convergence" << endl;
      cout << "   tolerance in the maximum number of iterations...... " << endl;
      cout << "   final gmrestol = " 
	   << setw(10) << GMRES_.get_rel_resid_norm() << endl;
      cout.flush();
    }
  } else {
    if (CFDkit_Primary_MPI_Processor()) { 
      cout << "\n  GMRES iter = " << GMRES_.get_iter() 
	   << "  gmrestol = " << setw(10) << GMRES_.get_rel_resid_norm() 
	   << endl;
    }    
  } /* endif */
  
  for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
    if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
      G[Bcount].tmp = NULL;
    }
  }

} /* End of GMRES_Algorithm. */

/********************************************************
 * Routine: GMRES_BCs                                   *
 *                                                      *
 * Apply boundary conditions at boundaries of the       *
 * specified GMRES vector.                              *
 *                                                      *
 ********************************************************/
void GMRES_BCs(Euler2D_Quad_Block &SolnBlk,
	       double * v, 
	       GMRES_Block &G,
	       int normalization) 
{
    int i, j;

    /* For BC_REFLECTION only. ***********************/
    DenseMatrix dUdW_(G.blocksize,G.blocksize, 0.0);
    DenseMatrix dWdU_(G.blocksize,G.blocksize, 0.0); 
    ColumnVector cv(G.blocksize, 0.0);
    Euler2D_pState Wo;
    /*************************************************/

    for ( j = SolnBlk.JCl-SolnBlk.Nghost  ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
      if ( (j >= SolnBlk.JCl && j <= SolnBlk.JCu) ||
           (j < SolnBlk.JCl && 
            (SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_NONE ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_PERIODIC ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_CONSTANT_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_LINEAR_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_CHARACTERISTIC) ) ||
           (j > SolnBlk.JCu && 
            (SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_NONE ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_PERIODIC ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_CONSTANT_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_LINEAR_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_CHARACTERISTIC) ) ) {
        switch(SolnBlk.Grid.BCtypeW[j]) {
	case BC_NONE :
	  break;
	case BC_FIXED :
	  for (int numNG = 1; numNG <= SolnBlk.Nghost; numNG++) {
	    for (int num = 0; num < G.blocksize; num++) {
	      v[G.index(SolnBlk.ICl-numNG,j)+num] = 0.0;
	    }
	  }
	  break;
	case BC_CONSTANT_EXTRAPOLATION :
	  for (int numNG = 1; numNG <= SolnBlk.Nghost; numNG++) {
	    for (int num = 0; num < G.blocksize; num++) {
	      v[G.index(SolnBlk.ICl-numNG,j)+num] = v[G.index(SolnBlk.ICl,j)+num];
	    }
	  }
	  break;
	case BC_LINEAR_EXTRAPOLATION :
	  for (int numNG = 1; numNG <= SolnBlk.Nghost; numNG++) {
	    for (int num = 0; num < G.blocksize; num++) {
	      v[G.index(SolnBlk.ICl-numNG,j)+num] = v[G.index(SolnBlk.ICl,j)+num];
	    }
	  }
	  break;
	case BC_REFLECTION :
	  for (int numNG = 1; numNG <= SolnBlk.Nghost; numNG++) {  
	    /*************************************************************/
	    /* Loaction: [SolnBlk.ICl][j]   -> [SolnBlk.ICl-1][j] */
	    /* Loaction: [SolnBlk.ICl+1][j] -> [SolnBlk.ICl-2][j] */
	    /* Copy values from v-vector to temporary ColumnVector, cv. */
	    cv(0) = v[G.index(SolnBlk.ICl+numNG-1,j)];
	    cv(1) = v[G.index(SolnBlk.ICl+numNG-1,j)+1];
	    cv(2) = v[G.index(SolnBlk.ICl+numNG-1,j)+2];
	    cv(3) = v[G.index(SolnBlk.ICl+numNG-1,j)+3];
	    
	    /* Copy converted Uo values to temporary storage, Wo. */
	    Wo = W(SolnBlk.Uo[SolnBlk.ICl+1][j]);
	    
	    /* Calculate dWdU using Uo values. */
	    dWdU_.zero();
	    dWdU(dWdU_, Wo);
	    if (normalization) normalize_dWdU(dWdU_);

	    /* Multiple dWdU matrix by the ColumnVector. */
	    cv = dWdU_ * cv;
	    
	    /* Calculate rotated values for the ColumnVector. */
	    cv = GMRES_Reflect(cv, SolnBlk.Grid.nfaceW(SolnBlk.ICl,j),
			       G.blocksize);
	    
	    /* Calculate dUdW using Uo values. */
	    dUdW_.zero();
	    dUdW(dUdW_, Wo);
	    if (normalization) normalize_dUdW(dUdW_);

	    /* Multiple dUdW matrix by the ColumnVector. */
	    cv = dUdW_ * cv;
	    
	    /* Copy solutions back to the v-vector. */
	    v[G.index(SolnBlk.ICl-numNG,j)]   = cv(0);
	    v[G.index(SolnBlk.ICl-numNG,j)+1] = cv(1);
	    v[G.index(SolnBlk.ICl-numNG,j)+2] = cv(2);
	    v[G.index(SolnBlk.ICl-numNG,j)+3] = cv(3);
	    /*************************************************************/
	  }
	  break;
	case BC_BURNING_SURFACE :
	  break;
	case BC_PERIODIC :
	  for (int numNG = 1; numNG <= SolnBlk.Nghost; numNG++) {
	    for (int num = 0; num < G.blocksize; num++) {
	      v[G.index(SolnBlk.ICl-numNG,j)+num] = 
		v[G.index(SolnBlk.ICu-numNG,j)+num];
	    }
	  }
	  break;
	case BC_CHARACTERISTIC :
	  for (int numNG = 1; numNG <= SolnBlk.Nghost; numNG++) {
	    for (int num = 0; num < G.blocksize; num++) {
	      v[G.index(SolnBlk.ICl-numNG,j)+num] = 0.0;
	    }
	  }
// 	  for (int numNG = 1; numNG <= SolnBlk.Nghost; numNG++) {
// 	    for (int num = 0; num < G.blocksize; num++) {
// 	      v[G.index(SolnBlk.ICl-numNG,j)+num] = v[G.index(SolnBlk.ICl,j)+num];
// 	    }
// 	  }
	  break;
	default:
	  for (int numNG = 1; numNG <= SolnBlk.Nghost; numNG++) {
	    for (int num = 0; num < G.blocksize; num++) {
	      v[G.index(SolnBlk.ICl-numNG,j)+num] = v[G.index(SolnBlk.ICl,j)+num];
	    }
	  }
	  break;
        } /* endswitch */
      } /* endif */

      if ( (j >= SolnBlk.JCl && j <= SolnBlk.JCu) ||
           (j < SolnBlk.JCl && 
            (SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_NONE ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_PERIODIC ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_CONSTANT_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_LINEAR_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_CHARACTERISTIC) ) ||
           (j > SolnBlk.JCu && 
            (SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_NONE ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_PERIODIC ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_CONSTANT_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_LINEAR_EXTRAPOLATION ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_CHARACTERISTIC) ) ) {
        switch(SolnBlk.Grid.BCtypeE[j]) {
	case BC_NONE :
	  break;
	case BC_FIXED :
	  for (int numNG = 1; numNG <= SolnBlk.Nghost; numNG++) {	    
	    for (int num = 0; num < G.blocksize; num++) {
	      v[G.index(SolnBlk.ICu+numNG,j)+num] = 0.0;
	    }
	  }
	  break;
	case BC_CONSTANT_EXTRAPOLATION :
	  for (int numNG = 1; numNG <= SolnBlk.Nghost; numNG++) {	
	    for (int num = 0; num < G.blocksize; num++) {
	      v[G.index(SolnBlk.ICu+numNG,j)+num] = v[G.index(SolnBlk.ICu,j)+num];
	    }
	  }
	  break;
	case BC_LINEAR_EXTRAPOLATION :
	  for (int numNG = 1; numNG <= SolnBlk.Nghost; numNG++) {	
	    for (int num = 0; num < G.blocksize; num++) {
	      v[G.index(SolnBlk.ICu+numNG,j)+num] = v[G.index(SolnBlk.ICu,j)+num];
	    }
	  }
	  break;
	case BC_REFLECTION :
	  for (int numNG = 1; numNG <= SolnBlk.Nghost; numNG++) {	
	    /*************************************************************/
	    /* Loaction: [SolnBlk.ICu][j]   -> [SolnBlk.ICu+1][j] */
	    /* Loaction: [SolnBlk.ICu-1][j] -> [SolnBlk.ICu+2][j] */
	    /* Copy values from v-vector to temporary ColumnVector, cv. */
	    cv(0) = v[G.index(SolnBlk.ICu-numNG-1,j)];
	    cv(1) = v[G.index(SolnBlk.ICu-numNG-1,j)+1];
	    cv(2) = v[G.index(SolnBlk.ICu-numNG-1,j)+2];
	    cv(3) = v[G.index(SolnBlk.ICu-numNG-1,j)+3];
	    
	    /* Copy converted Uo values to temporary storage, Wo. */
	    Wo = W(SolnBlk.Uo[SolnBlk.ICu-1][j]);
	    
	    /* Calculate dWdU using Uo values. */
	    dWdU_.zero();
	    dWdU(dWdU_, Wo);
	    if (normalization) normalize_dWdU(dWdU_);

	    /* Multiple dWdU matrix by the ColumnVector. */
	    cv = dWdU_ * cv;
	    
	    /* Calculate rotated values for the ColumnVector. */
	    cv = GMRES_Reflect(cv, SolnBlk.Grid.nfaceE(SolnBlk.ICu,j),
			       G.blocksize);
	    
	    /* Calculate dUdW using Uo values. */
	    dUdW_.zero();
	    dUdW(dUdW_, Wo);
	    if (normalization) normalize_dUdW(dUdW_);

	    /* Multiple dUdW matrix by the ColumnVector. */
	    cv = dUdW_ * cv;
	    
	    /* Copy solutions back to the v-vector. */
	    v[G.index(SolnBlk.ICu+numNG,j)]   = cv(0);
	    v[G.index(SolnBlk.ICu+numNG,j)+1] = cv(1);
	    v[G.index(SolnBlk.ICu+numNG,j)+2] = cv(2);
	    v[G.index(SolnBlk.ICu+numNG,j)+3] = cv(3);
	    /*************************************************************/
	  }
	  break;
	case BC_BURNING_SURFACE :
	  break;
	case BC_PERIODIC :
	  for (int numNG = 1; numNG <= SolnBlk.Nghost; numNG++) {	
	    for (int num = 0; num < G.blocksize; num++) {
	      v[G.index(SolnBlk.ICu+numNG,j)+num] = 
		v[G.index(SolnBlk.ICl+numNG,j)+num];
	    }
	  }
	  break;
	case BC_CHARACTERISTIC :
	  for (int numNG = 1; numNG <= SolnBlk.Nghost; numNG++) {	    
	    for (int num = 0; num < G.blocksize; num++) {
	      v[G.index(SolnBlk.ICu+numNG,j)+num] = 0.0;
	    }
	  }
// 	  for (int numNG = 1; numNG <= SolnBlk.Nghost; numNG++) {	
// 	    for (int num = 0; num < G.blocksize; num++) {
// 	      v[G.index(SolnBlk.ICu+numNG,j)+num] = v[G.index(SolnBlk.ICu,j)+num];
// 	    }
// 	  }
	  break;
	default:
	  for (int numNG = 1; numNG <= SolnBlk.Nghost; numNG++) {	
	    for (int num = 0; num < G.blocksize; num++) {
	    v[G.index(SolnBlk.ICu+numNG,j)+num] = v[G.index(SolnBlk.ICu,j)+num];
	    }
	  }
	  break;
        } /* endswitch */
      } /* endif */
    } /* endfor */

    for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
        switch(SolnBlk.Grid.BCtypeS[i]) {
          case BC_NONE :
            break;
          case BC_FIXED :
	    for (int numNG = 1; numNG <= SolnBlk.Nghost; numNG++) {	
	      for (int num = 0; num < G.blocksize; num++) {
		v[G.index(i,SolnBlk.JCl-numNG)+num] = 0.0;
	      }
	    }
            break;
 	  case BC_CONSTANT_EXTRAPOLATION :
	    for (int numNG = 1; numNG <= SolnBlk.Nghost; numNG++) {	   
	      for (int num = 0; num < G.blocksize; num++) {
	        v[G.index(i,SolnBlk.JCl-numNG)+num] = v[G.index(i,SolnBlk.JCl)+num];
	      }
	    }
	    break;
	  case BC_LINEAR_EXTRAPOLATION :
	    for (int numNG = 1; numNG <= SolnBlk.Nghost; numNG++) {	   
	      for (int num = 0; num < G.blocksize; num++) {
	        v[G.index(i,SolnBlk.JCl-numNG)+num] = v[G.index(i,SolnBlk.JCl)+num];
	      }
	    }
	    break;
	  case BC_REFLECTION :
	    for (int numNG = 1; numNG <= SolnBlk.Nghost; numNG++) {	
	      /*************************************************************/
	      /* Loaction: [i][SolnBlk.JCl] -> [i][SolnBlk.JCl-1] */
	      /* Loaction: [i][SolnBlk.JCl+1] -> [i][SolnBlk.JCl-2] */
	      /* Copy values from v-vector to temporary ColumnVector, cv. */
	      cv(0) = v[G.index(i,SolnBlk.JCl+numNG-1)];
	      cv(1) = v[G.index(i,SolnBlk.JCl+numNG-1)+1];
	      cv(2) = v[G.index(i,SolnBlk.JCl+numNG-1)+2];
	      cv(3) = v[G.index(i,SolnBlk.JCl+numNG-1)+3];

	      /* Copy converted Uo values to temporary storage, Wo. */
	      Wo = W(SolnBlk.Uo[i][SolnBlk.JCl+1]);

  	      /* Calculate dWdU using Uo values. */
	      dWdU_.zero();
	      dWdU(dWdU_, Wo);
	      if (normalization) normalize_dWdU(dWdU_);

	      /* Multiple dWdU matrix by the ColumnVector. */
	      cv = dWdU_ * cv;

	      /* Calculate rotated values for the ColumnVector. */
	      cv = GMRES_Reflect(cv, SolnBlk.Grid.nfaceS(i,SolnBlk.JCl), 
		  	         G.blocksize);

	      /* Calculate dUdW using Uo values. */
	      dUdW_.zero();
	      dUdW(dUdW_, Wo);
	      if (normalization) normalize_dUdW(dUdW_);

	      /* Multiple dUdW matrix by the ColumnVector. */
	      cv = dUdW_ * cv;

	      /* Copy solutions back to the v-vector. */
	      v[G.index(i,SolnBlk.JCl-numNG)]   = cv(0);
	      v[G.index(i,SolnBlk.JCl-numNG)+1] = cv(1);
	      v[G.index(i,SolnBlk.JCl-numNG)+2] = cv(2);
	      v[G.index(i,SolnBlk.JCl-numNG)+3] = cv(3);
	      /*************************************************************/
	    }
	    break;
	  case BC_BURNING_SURFACE :
	    break;
	  case BC_PERIODIC :
	    for (int numNG = 1; numNG <= SolnBlk.Nghost; numNG++) {	
	      for (int num = 0; num < G.blocksize; num++) {
	        v[G.index(i,SolnBlk.JCl-numNG)+num] = 
		  v[G.index(i,SolnBlk.JCu-numNG)+num];
	      }
	    }
	    break;
	  case BC_CHARACTERISTIC :
	    for (int numNG = 1; numNG <= SolnBlk.Nghost; numNG++) {	
	      for (int num = 0; num < G.blocksize; num++) {
		v[G.index(i,SolnBlk.JCl-numNG)+num] = 0.0;
	      }
	    }
// 	    for (int numNG = 1; numNG <= SolnBlk.Nghost; numNG++) {	   
// 	      for (int num = 0; num < G.blocksize; num++) {
// 	        v[G.index(i,SolnBlk.JCl-numNG)+num] = v[G.index(i,SolnBlk.JCl)+num];
// 	      }
// 	    }
	    break;
	  default:
	    for (int numNG = 1; numNG <= SolnBlk.Nghost; numNG++) {	
	      for (int num = 0; num < G.blocksize; num++) {
	        v[G.index(i,SolnBlk.JCl-numNG)+num] = v[G.index(i,SolnBlk.JCl)+num];
	      }
	    }
	    break;
        } /* endswitch */
	
        switch(SolnBlk.Grid.BCtypeN[i]) {
  	  case BC_NONE :
	    break;
	  case BC_FIXED :
	    for (int numNG = 1; numNG <= SolnBlk.Nghost; numNG++) {	
	      for (int num = 0; num < G.blocksize; num++) {
	        v[G.index(i,SolnBlk.JCu+numNG)+num] = 0.0;
	      }
	    }
	    break;
          case BC_CONSTANT_EXTRAPOLATION :
	    for (int numNG = 1; numNG <= SolnBlk.Nghost; numNG++) {	
	      for (int num = 0; num < G.blocksize; num++) {
	        v[G.index(i,SolnBlk.JCu+numNG)+num] = v[G.index(i,SolnBlk.JCu)+num];
	      }
	    }
            break;
	  case BC_LINEAR_EXTRAPOLATION :
	    for (int numNG = 1; numNG <= SolnBlk.Nghost; numNG++) {	
	      for (int num = 0; num < G.blocksize; num++) {
	        v[G.index(i,SolnBlk.JCu+numNG)+num] = 0.0;
	      }
	    }
	    break;
	  case BC_REFLECTION :
	    for (int numNG = 1; numNG <= SolnBlk.Nghost; numNG++) {	
	      /*************************************************************/
	      /* Loaction: [i][SolnBlk.JCu]   -> [i][SolnBlk.JCu+1] */
	      /* Loaction: [i][SolnBlk.JCu-1] -> [i][SolnBlk.JCu+2] */
	      /* Copy values from v-vector to temporary ColumnVector, cv. */
	      cv(0) = v[G.index(i,SolnBlk.JCu-numNG-1)];
	      cv(1) = v[G.index(i,SolnBlk.JCu-numNG-1)+1];
	      cv(2) = v[G.index(i,SolnBlk.JCu-numNG-1)+2];
	      cv(3) = v[G.index(i,SolnBlk.JCu-numNG-1)+3];

	      /* Copy converted Uo values to temporary storage, Wo. */
	      Wo = W(SolnBlk.Uo[i][SolnBlk.JCu-1]);

	      /* Calculate dWdU using Uo values. */
	      dWdU_.zero();
	      dWdU(dWdU_, Wo);
	      if (normalization) normalize_dWdU(dWdU_);

	      /* Multiple dWdU matrix by the ColumnVector. */
	      cv = dWdU_ * cv;

	      /* Calculate rotated values for the ColumnVector. */
	      cv = GMRES_Reflect(cv, SolnBlk.Grid.nfaceN(i,SolnBlk.JCu), 
			       G.blocksize);

	      /* Calculate dUdW using Uo values. */
	      dUdW_.zero();
	      dUdW(dUdW_, Wo);
	      if (normalization) normalize_dUdW(dUdW_);

	      /* Multiple dUdW matrix by the ColumnVector. */
	      cv = dUdW_ * cv;

	      /* Copy solutions back to the v-vector. */
	      v[G.index(i,SolnBlk.JCu+numNG)]   = cv(0);
	      v[G.index(i,SolnBlk.JCu+numNG)+1] = cv(1);
	      v[G.index(i,SolnBlk.JCu+numNG)+2] = cv(2);
	      v[G.index(i,SolnBlk.JCu+numNG)+3] = cv(3);
	      /*************************************************************/
	    }
	    break;
	  case BC_BURNING_SURFACE :
	    break;
	  case BC_PERIODIC :
	    for (int numNG = 1; numNG <= SolnBlk.Nghost; numNG++) {	
	       for (int num = 0; num < G.blocksize; num++) {
	         v[G.index(i,SolnBlk.JCu+numNG)+num] = 
		  v[G.index(i,SolnBlk.JCl+numNG)+num];
	       }
	    }
	    break;
	  case BC_CHARACTERISTIC :
	    for (int numNG = 1; numNG <= SolnBlk.Nghost; numNG++) {	
	      for (int num = 0; num < G.blocksize; num++) {
	        v[G.index(i,SolnBlk.JCu+numNG)+num] = 0.0;
	      }
	    }
// 	    for (int numNG = 1; numNG <= SolnBlk.Nghost; numNG++) {	
// 	      for (int num = 0; num < G.blocksize; num++) {
// 	        v[G.index(i,SolnBlk.JCu+numNG)+num] = v[G.index(i,SolnBlk.JCu)+num];
// 	      }
// 	    }
	    break;
          default:
	    for (int numNG = 1; numNG <= SolnBlk.Nghost; numNG++) {	
	      for (int num = 0; num < G.blocksize; num++) {
	        v[G.index(i,SolnBlk.JCu+numNG)+num] = v[G.index(i,SolnBlk.JCu)+num];
	      }
	    }
            break;
        } /* endswitch */
    } /* endfor */
    
    /* BC fix for corner points with burning surfaces on either side. */

    if (SolnBlk.Grid.BCtypeW[SolnBlk.JCl] == BC_BURNING_SURFACE &&
        SolnBlk.Grid.BCtypeS[SolnBlk.ICl] == BC_BURNING_SURFACE) {
       SolnBlk.W[SolnBlk.ICl-1][SolnBlk.JCl-1] = HALF*(SolnBlk.W[SolnBlk.ICl-1][SolnBlk.JCl]+
                                                       SolnBlk.W[SolnBlk.ICl][SolnBlk.JCl-1]);
       SolnBlk.U[SolnBlk.ICl-1][SolnBlk.JCl-1] = U(SolnBlk.W[SolnBlk.ICl-1][SolnBlk.JCl-1]);
    } /* endif */

    if (SolnBlk.Grid.BCtypeW[SolnBlk.JCu] == BC_BURNING_SURFACE &&
        SolnBlk.Grid.BCtypeN[SolnBlk.ICl] == BC_BURNING_SURFACE) {
       SolnBlk.W[SolnBlk.ICl-1][SolnBlk.JCu+1] = HALF*(SolnBlk.W[SolnBlk.ICl-1][SolnBlk.JCu]+
                                                       SolnBlk.W[SolnBlk.ICl][SolnBlk.JCu+1]);
       SolnBlk.U[SolnBlk.ICl-1][SolnBlk.JCu+1] = U(SolnBlk.W[SolnBlk.ICl-1][SolnBlk.JCu+1]);
    } /* endif */

    if (SolnBlk.Grid.BCtypeE[SolnBlk.JCl] == BC_BURNING_SURFACE &&
        SolnBlk.Grid.BCtypeS[SolnBlk.ICu] == BC_BURNING_SURFACE) {
       SolnBlk.W[SolnBlk.ICu+1][SolnBlk.JCl-1] = HALF*(SolnBlk.W[SolnBlk.ICu+1][SolnBlk.JCl]+
                                                       SolnBlk.W[SolnBlk.ICu][SolnBlk.JCl-1]);
       SolnBlk.U[SolnBlk.ICu+1][SolnBlk.JCl-1] = U(SolnBlk.W[SolnBlk.ICu+1][SolnBlk.JCl-1]);
    } /* endif */

    if (SolnBlk.Grid.BCtypeE[SolnBlk.JCu] == BC_BURNING_SURFACE &&
        SolnBlk.Grid.BCtypeN[SolnBlk.ICu] == BC_BURNING_SURFACE) {
       SolnBlk.W[SolnBlk.ICu+1][SolnBlk.JCu+1] = HALF*(SolnBlk.W[SolnBlk.ICu+1][SolnBlk.JCu]+
                                                       SolnBlk.W[SolnBlk.ICu][SolnBlk.JCu+1]);
       SolnBlk.U[SolnBlk.ICu+1][SolnBlk.JCu+1] = U(SolnBlk.W[SolnBlk.ICu+1][SolnBlk.JCu+1]);
    } /* endif */

} /* End of GMRES_BCs. */

/********************************************************
 * Routine: GMRES_Reflect                               *
 *                                                      *
 * This function returns the reflected solution state   *
 * in a given direction given the specified gmres       *
 * vector and the unit normal vector in the             *
 * direction of interest.                               *
 *                                                      *
 ********************************************************/
ColumnVector GMRES_Reflect(const ColumnVector cv,
			   const Vector2D &norm_dir,
			   int blocksize) {
  
  double ur, vr, u, v;
  double cos_angle, sin_angle;
  
  ColumnVector result(blocksize,0.0);
  u = cv(1);
  v = cv(2);

  /* Determine the direction cosine's for the frame rotation. */
  cos_angle = norm_dir.x; 
  sin_angle = norm_dir.y;
  
  /* Apply the frame rotation and calculate the primitive
     solution state variables in the local rotated frame
     defined by the unit normal vector. */
  
  ur =   u*cos_angle + v*sin_angle;
  vr = - u*sin_angle + v*cos_angle;

  /* Reflect the normal velocity in the rotated frame. */
  
  ur = -ur;
  
  /* Rotate back to the original Cartesian reference frame. */
  
  result(0) = cv(0);
  result(1) = ur*cos_angle - vr*sin_angle;
  result(2) = ur*sin_angle + vr*cos_angle;
  result(3) = cv(3);

  /* Return the reflected state. */  
  return (result);
  
} /* End of GMRES_Reflect. */

/********************************************************
 * Routine: normalize_RHS                               *
 *                                                      *
 * This routine returns the normalized RHS.             *
 *                                                      *
 ********************************************************/
void normalize_RHS(double * RHS,
		   Euler2D_Quad_Block &SolnBlk) {

  double ao  = Euler2D_W_STDATM.a();
  double rho = Euler2D_W_STDATM.d;  

  int xpts = (SolnBlk.ICu-SolnBlk.ICl+1)+SolnBlk.Nghost * 2;

  /* Normalize the solution residuals to "RHS" vector. */
    int index = 0;
    for (int j = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; i++) {
	index = (j*xpts+i)*NUM_VAR_EULER2D;
	RHS[index] = RHS[index] / (rho * ao);      index++;      
	RHS[index] = RHS[index] / (rho * ao * ao); index++;
	RHS[index] = RHS[index] / (rho * ao * ao); index++;
	RHS[index] = RHS[index] / (rho * ao * ao * ao);

      } /* endfor */
    } /* endfor */

}/* End of normalize_RHS. */

/********************************************************
 * Routine: denormalize_solution                        *
 *                                                      *
 * This routine returns denormalized solution vector.   *
 *                                                      *
 ********************************************************/
void denormalize_solution(double * SolVec,
		          Euler2D_Quad_Block &SolnBlk) {

  double ao  = Euler2D_W_STDATM.a();
  double rho = Euler2D_W_STDATM.d;
  
  int xpts = (SolnBlk.ICu-SolnBlk.ICl+1) + SolnBlk.Nghost * 2;
  
  /* Denormalize the solution vector. */
    int index = 0;
    for (int j = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; j++) {
      for (int i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; i++) {
	index = (j*xpts+i)*NUM_VAR_EULER2D;
	SolVec[index] = SolVec[index] * rho;      index++;
	SolVec[index] = SolVec[index] * rho * ao; index++;
	SolVec[index] = SolVec[index] * rho * ao; index++;
	SolVec[index] = SolVec[index] * rho * ao * ao;

      } /* endfor */
    } /* endfor */

}/* End of denormalize_solution. */

/********************************************************
 * Routine: normalize_dUdW                              *
 *                                                      *
 * This function returns the dUdW matrix.               * 
 *                                                      *
 ********************************************************/
void normalize_dUdW(DenseMatrix &dUdW) {

double ao  = Euler2D_W_STDATM.ao();
double rho = Euler2D_W_STDATM.d;

dUdW(0,1) = dUdW(0,1) * (ao/rho);
dUdW(0,2) = dUdW(0,2) * (ao/rho);
dUdW(0,3) = dUdW(0,3) * sqr(ao);
dUdW(1,0) = dUdW(1,0) * (1.0/ao);
dUdW(1,1) = dUdW(1,1) * (1.0/rho);
dUdW(1,2) = dUdW(1,2) * (1.0/rho);
dUdW(1,3) = dUdW(1,3) * (ao);
dUdW(2,0) = dUdW(2,0) * (1.0/ao);
dUdW(2,1) = dUdW(2,1) * (1.0/rho);
dUdW(2,2) = dUdW(2,2) * (1.0/rho);
dUdW(2,3) = dUdW(2,3) * (ao);
dUdW(3,0) = dUdW(3,0) * (1.0/sqr(ao));
dUdW(3,1) = dUdW(3,1) * (1.0/(rho*ao));
dUdW(3,2) = dUdW(3,2) * (1.0/(rho*ao));

}/* End of normalize_dUdW. */

/********************************************************
 * Routine: normalize_dWdU                              *
 *                                                      *
 * This function returns the dWdU matrix.               * 
 *                                                      *
 ********************************************************/
void normalize_dWdU(DenseMatrix &dWdU) {

double ao  = Euler2D_W_STDATM.ao();
double rho = Euler2D_W_STDATM.d;

dWdU(0,1) = dWdU(0,1) * (ao);
dWdU(0,2) = dWdU(0,2) * (ao);
dWdU(0,3) = dWdU(0,3) * sqr(ao);
dWdU(1,0) = dWdU(1,0) * (rho/ao);
dWdU(1,1) = dWdU(1,1) * (rho);
dWdU(1,2) = dWdU(1,2) * (rho);
dWdU(1,3) = dWdU(1,3) * (rho*ao);
dWdU(2,0) = dWdU(2,0) * (rho/ao);
dWdU(2,1) = dWdU(2,1) * (rho);
dWdU(2,2) = dWdU(2,2) * (rho);
dWdU(2,3) = dWdU(2,3) * (rho*ao);
dWdU(3,0) = dWdU(3,0) * (1.0/sqr(ao));
dWdU(3,1) = dWdU(3,1) * (1.0/ao);
dWdU(3,2) = dWdU(3,2) * (1.0/ao);

}/* End of normalize_dWdU. */

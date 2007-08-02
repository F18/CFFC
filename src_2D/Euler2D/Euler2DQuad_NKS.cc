/* Euler2DQuad_NKS.cc:  Functions for NKS Solver to solve 2D Euler Equation. */

/* Include 2D Euler equation quadrilateral NKS header file. */

#ifndef _EULER2D_QUAD_NKS_INCLUDED 
#include "Euler2DQuad_NKS.h" 
#endif // _EULER2D_QUAD_NKS_INCLUDED 

/********************************************************
 * Routine: Newton_Krylov_Solver                        *
 *                                                      *
 * This routine updates the specified solution block    *
 * using Newton-Krylov-Schwarz method.                  *
 *                                                      *
 ********************************************************/ 
int Newton_Krylov_Schwarz_Solver(CPUTime &NKS_time,
				 ostream &Progress_File,
                                 int &Number_of_Startup_Interations,
				 Euler2D_Quad_Block *Soln_ptr,
				 AdaptiveBlock2D_List &Soln_Block_List,
				 Euler2D_Input_Parameters &Input_Parameters) { 

  /* Local variable definitions. */

  /*Overall Convergence tolerance. */
  double overall_tol  = Input_Parameters.Overall_Toler;              
  /* GMRES Convergence tolerance. */
  double gmrestol     = Input_Parameters.GMRES_Toler;
  /* GMRES Restart. */
  int restart         = Input_Parameters.GMRES_Restart;
  /* GMRES Level of overlap. */
  int overlap         = Input_Parameters.GMRES_Overlap;
  /* GMRES Level of fill. */
  int level_of_fill   = Input_Parameters.GMRES_ILUK_Level_of_Fill;
  /* Maximum Newton Iterations. */
  int max_newton_step = Input_Parameters.Maximum_Number_of_NKS_Iterations;
  /* Maximum GMRES Iterations. */
  int max_gmres_iter  = Input_Parameters.Maximum_Number_of_GMRES_Iterations;
  /* Switch for preconditioner. */
  int P_Switch        = Input_Parameters.GMRES_P_Switch;

  double L2norm_current   = 0.0, 
         L1norm_current   = 0.0, 
         Max_norm_current = 0.0; 
  double L2norm_first     = 0.0, 
         L1norm_first     = 0.0,
         Max_norm_first   = 0.0;
  double L2norm_current_n = 0.0, 
         L1norm_current_n = 0.0,
         Max_norm_current_n = 0.0;
  double norm_ratio       = 0.0;
  double dTime, CFL_current;

  int error_flag;
  int count = 0;
  int Number_of_Newton_Steps = 0;
  int stop  = 1;
  int flag  = ON;
  int NBLK  = Soln_Block_List.Nblk;
  int finite_time_step = Input_Parameters.Finite_Time_Step;
  int normalization    = Input_Parameters.Normalization;

  int    limiter_check    = ON;
  int    Freeze_Limiter = Input_Parameters.Freeze_Limiter;
  double Freeze_Limiter_Residual_Level = Input_Parameters.Freeze_Limiter_Residual_Level;
  int    i_limiter;

  int ni, i, j, num = 0;
  int xcount, ycount;
  int Bcount;
  int xpts, ypts, zpts = 1;
  int index, indexRHS;

  int n_update_reduction;
  double update_reduction_factor;

  if (CFDkit_Primary_MPI_Processor()) {
    cout << " " << endl;
    for (int star=0;star<75;star++){cout <<"*";}
    cout << "\n**************       BPKIT BPKIT BPKIT BPKIT BPKIT BPKIT   ****************" << endl;
    cout << "********                   Newton-Krylov-Schwarz                 **********" << endl;   
    for (int star=0;star<75;star++){cout<<"*";}
    cout << "\nLimiter Type          ====> ";
    if (Input_Parameters.i_Limiter == LIMITER_ZERO)   { 
      cout << "ZERO" << endl;
    } else if (Input_Parameters.i_Limiter == LIMITER_BARTH_JESPERSEN) { 
      cout <<"BARTH_JESPERSEN" << endl;
    } else if (Input_Parameters.i_Limiter == LIMITER_VENKATAKRISHNAN) { 
      cout << "VENKATAKRISHNAN" << endl;
    } else { 
      cout << "??" << endl;
    } /* endif */
    if (P_Switch == 1) {         // ILU(0) 
      cout << "Local Preconditioner  ====> ILU("<< level_of_fill <<")" << endl;
    } else if (P_Switch == 2){   // Diagonal
      cout << "Local Preconditioner  ====> Diagonal" << endl; 
    } /* endif */  
    cout <<     "Overall Tolerance     ====> " << overall_tol << endl;
    cout <<     "GMRES Tolerance       ====> " << gmrestol << endl;
    cout <<     "GMRES Restart         ====> " << restart << endl;
    cout <<     "Level of Overlap      ====> " << overlap << endl;
    cout <<     "Number of Processors  ====> " 
	 << Input_Parameters.Number_of_Processors << endl;
    cout <<     "Overall Block Size    ====> " 
	 << Input_Parameters.Number_of_Blocks_Idir << " x " 
	 << Input_Parameters.Number_of_Blocks_Jdir << endl;

    if (finite_time_step == ON) {
      cout << "Finite Time Step      ====> ON" << endl;
      cout << "Initial_CFL           ====> " 
	   << Input_Parameters.Finite_Time_Step_Initial_CFL << endl;
    } else {
      cout << "Finite Time Step      ====> OFF" << endl; 
    } /* endif */ 
    if (Freeze_Limiter == ON) {
      cout << "Freeze_Limiter        ====> ON" << endl;
      cout << "Residual_Level        ====> "
	   << Input_Parameters.Freeze_Limiter_Residual_Level << endl;
    } else {
      cout << "Freeze_Limiter      ====> OFF" << endl; 
    } /* endif */
    if (normalization == ON) {
      cout << "Normalization         ====> ON" << endl;
    } else {
      cout << "Normalization         ====> OFF" << endl; 
    } /* endif */ 
    for (int star=0;star<75;star++){cout<<"*";}
    cout << "  " << endl;

  } /* endif (CFDkit_Primary_MPI_Processor())  */ 

  /***********************************************************************/
  /* fill_pattern = 1 ->  5 coefficients(nearest neighboring cells only) */
  /*                2 -> 13 coefficients(both nearest and second nearest */
  /* 		                         neighboring cells)              */
  /*                                                                     */
  /* Number of Variables = NUM_VAR_EULER2D (4 variables in Euler         */
  /*                                        Equations)                   */                              
  /***********************************************************************/
  int fill_pattern = 1;
  int blocksize    = NUM_VAR_EULER2D;

  /* Create local variables for constructing Jacobian matrices. **************/
  DenseMatrix   dFdU(blocksize,blocksize, 0.0);
  DenseMatrix   zero(blocksize,blocksize, 0.0);
  Index M_ind_CENTER, M_ind_NORTH, M_ind_SOUTH, M_ind_EAST, M_ind_WEST;
  /**************************************************************************/
  
  /* Create GMRES vector for each block. ************************************/
  
  GMRES_Block * G = new GMRES_Block[NBLK];

  /* Initialize variables in GMRES_Block objects. */
  for ( Bcount = 0 ; Bcount < NBLK; ++Bcount ) {
    if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
      count++;
      xpts = (Soln_ptr[Bcount].ICu-Soln_ptr[Bcount].ICl+1) + 
	      Soln_ptr[Bcount].Nghost * 2;
      ypts = (Soln_ptr[Bcount].JCu-Soln_ptr[Bcount].JCl+1) + 
	      Soln_ptr[Bcount].Nghost * 2;
      G[Bcount].allocate(restart, blocksize, xpts, ypts);
      G[Bcount].xpts           = xpts;
      G[Bcount].ypts           = ypts;
      G[Bcount].m              = restart;
      G[Bcount].overlap        = overlap;
      G[Bcount].P_Switch       = P_Switch;
      G[Bcount].blocksize      = blocksize;
      G[Bcount].max_gmres_iter = max_gmres_iter;

      G[Bcount].NCi = Soln_ptr[Bcount].NCi;
      G[Bcount].ICl = Soln_ptr[Bcount].ICl;
      G[Bcount].ICu = Soln_ptr[Bcount].ICu;
      G[Bcount].NCj = Soln_ptr[Bcount].NCj;
      G[Bcount].JCl = Soln_ptr[Bcount].JCl;
      G[Bcount].JCu = Soln_ptr[Bcount].JCu;
      G[Bcount].Nghost = Soln_ptr[Bcount].Nghost;

      G[Bcount].SolBlk_ptr = &(Soln_ptr[Bcount]);

    } /* endif */
  } /* endfor */
 /**************************************************************************/

  /** Create BlockMatrix Objects. *******************************************/
  HBTMat * model = new HBTMat[count]; 
  HBTMat_Class hbtmat_object(blocksize, fill_pattern);

  for (int NUM=0; NUM<count ; NUM++)
    model[NUM] = 
      Create_Block_Matrix(G[NUM].xpts, G[NUM].ypts, zpts, hbtmat_object);

  BlockMatrix * blockmatrix = new BlockMatrix[count];
  for (i=0; i<count; i++)
    blockmatrix[i].create(G[i].xpts, G[i].ypts, zpts, 
			 model[i], blocksize, fill_pattern);
  /**************************************************************************/

  /** Create Preconditioner Objects. ****************************************/
  BILUK   * MBILUK   = new   BILUK[NBLK];
  BJacobi * MBJacobi = new BJacobi[NBLK];
  /**************************************************************************/
  
  /*******************************************/
  /* BEGIN NEWTON-KRYLOV-SCHWARZ CALCULATION */
  /*******************************************/
  Number_of_Newton_Steps = 1;
  i_limiter = Input_Parameters.i_Limiter;
  while (stop && Number_of_Newton_Steps <= max_newton_step) {

    /* Limiter Switch. */
    if (Number_of_Newton_Steps <= 5) {
      //Input_Parameters.i_Limiter = LIMITER_ZERO;
    } else {
      Input_Parameters.i_Limiter = i_limiter;
    } /* endif */

    /* Calculate residual : use dudt[i][j][0] for residual calculations. */
    for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
      if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	for (num=0;num<G[Bcount].scalar_dim;++num) G[Bcount].b[num] = 0.0;

	dUdt_Residual_Evaluation(Soln_ptr[Bcount],
				 Input_Parameters);
      } /* endif */
    } /* endfor */

    /* Send boundary flux corrections at block interfaces with resolution changes. */
    error_flag = Send_Conservative_Flux_Corrections(Soln_ptr, 
		                		    Soln_Block_List,
						    NUM_VAR_EULER2D);
    if (error_flag) {
       cout << "\n Euler2D NKS ERROR: AdvectDiffuse2D flux correction message passing error on processor "
            << Soln_Block_List.ThisCPU
            << ".\n";
       cout.flush();
    } /* endif */
    error_flag = CFDkit_OR_MPI(error_flag);
    if (error_flag) return (error_flag);
	  
    /* Apply boundary flux corrections to residual to ensure that method is conservative. */
    Apply_Boundary_Flux_Corrections(Soln_ptr, 
		                    Soln_Block_List);

    for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
      if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	G[Bcount].copy_residual(G[Bcount].b, Soln_ptr[Bcount]);
	
	/* RHS: (dR/dU)(delta U) = -R */
     	for (i=0; i<G[Bcount].scalar_dim; ++i)
	  G[Bcount].b[i] = -1.0 * G[Bcount].b[i];
	
      } /* endif */
    } /* endfor */

    /* Calculate 1-, 2-, and max-norms or residual for all blocks. */   
    L2norm_current   = sqr(L2_Norm_Residual(Soln_ptr, Soln_Block_List));
    L2norm_current   = sqrt(CFDkit_Summation_MPI(L2norm_current));
    L1norm_current   = L1_Norm_Residual(Soln_ptr, Soln_Block_List);
    L1norm_current   = CFDkit_Summation_MPI(L1norm_current);
    Max_norm_current = Max_Norm_Residual(Soln_ptr, Soln_Block_List);
    Max_norm_current = CFDkit_Summation_MPI(Max_norm_current);

    if (Number_of_Newton_Steps == 1) {
      L2norm_first = L2norm_current;
      L1norm_first = L1norm_current;
      Max_norm_first = Max_norm_current;
    } /* endif */ 
    
    L2norm_current_n   = L2norm_current / L2norm_first; 
    L1norm_current_n   = L1norm_current / L1norm_first; 
    Max_norm_current_n = Max_norm_current / Max_norm_first;

    /* Calculate ratio of initial and current 2-norms. */
    norm_ratio = L2norm_first / L2norm_current;

    /* Output info to file. */
    NKS_time.update();
    Output_Progress_to_File(Progress_File,
			    Number_of_Startup_Interations+Number_of_Newton_Steps-1,
			    ZERO,
			    NKS_time,
			    L1norm_current,
			    L2norm_current,
			    Max_norm_current);

    /* Apply Freeze Limiter */
    if (Freeze_Limiter) {
      if (L2norm_current_n <= Freeze_Limiter_Residual_Level && 
	  limiter_check == ON)  {
	if (CFDkit_Primary_MPI_Processor()) cout << "********** Apply Limiter Freezing **********" << endl;

	Freeze_Limiters(Soln_ptr, Soln_Block_List);
	
	limiter_check = OFF;	
      } /* endif */
    } /* endif */
    
    if (L2norm_current_n > overall_tol){  
      
      /* Calculate delta t = min(delta t) among processors. */
      dTime = CFL(Soln_ptr, Soln_Block_List, Input_Parameters);
      dTime = CFDkit_Minimum_MPI(dTime); 

      /* Apply finite time step. */
      if (finite_time_step) {
        if (norm_ratio <= 1.0e08) {
	   CFL_current = Input_Parameters.Finite_Time_Step_Initial_CFL*
                         pow(min(ONE, max(ONE, norm_ratio)/1.0e08), ONE);
        } else {
	   CFL_current = Input_Parameters.Finite_Time_Step_Initial_CFL;
        } /* endif */
	dTime = CFL_current*dTime;
      } else {
        CFL_current = Input_Parameters.Finite_Time_Step_Initial_CFL;
	dTime = CFL_current*dTime;
      }
      
//       if (finite_time_step) {
//         if (Number_of_Newton_Steps <= 100) {
//  	   CFL_current = Input_Parameters.Finite_Time_Step_Initial_CFL*
//                       pow((double(Number_of_Newton_Steps)/100.0), 4.0);
//         } else {
// 	   CFL_current = Input_Parameters.Finite_Time_Step_Initial_CFL;
//         } /* endif */
// 	   dTime = CFL_current*dTime;
//       } else {
//         CFL_current = Input_Parameters.Finite_Time_Step_Initial_CFL;
// 	   dTime = CFL_current*dTime;
//       } /* endif */
//       if (CFDkit_Primary_MPI_Processor()) 
//          cout << " CFL_current = " << CFL_current << " modified dTime = " << dTime << endl;        
//       Set_Global_TimeStep(Soln_ptr, Soln_Block_List, dTime);

      /* Print out info at the beginning of the iteration. */
      if (CFDkit_Primary_MPI_Processor()) {      
	cout << "\nBegin Newton Step (Outer Iterations) = "
	     << Number_of_Newton_Steps << " L2norm = "
	     << L2norm_current << " norm_ratio = " 
             << norm_ratio << " CFL = " << CFL_current << endl;
      } /* endif */

      /* Create/Update Jacobian Matrix - blockmatrix. */
      if ((Number_of_Newton_Steps < 10 || L2norm_current_n > 1.0e-02 ) && flag == ON) { 
	if (CFDkit_Primary_MPI_Processor()) cout << "  **** Create/Update Jacobian Matrix **** " << endl;  
	for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
	  if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	    
	  G[Bcount].dTime = dTime;
	  
	  for (j  = Soln_ptr[Bcount].JCl-Soln_ptr[Bcount].Nghost; 
	       j <= Soln_ptr[Bcount].JCu+Soln_ptr[Bcount].Nghost; ++j) {
	    for (i = Soln_ptr[Bcount].ICl-Soln_ptr[Bcount].Nghost; 
		 i <= Soln_ptr[Bcount].ICu+Soln_ptr[Bcount].Nghost; ++i) {
	      
	      /* Find the Jacobian Matrix location (x,y) based on (i,j) location from the mesh. */ 
	      M_ind_CENTER = blockmatrix[Bcount].Mindex_2D(i, j, CENTER);
	      M_ind_NORTH  = blockmatrix[Bcount].Mindex_2D(i, j, NORTH);
	      M_ind_SOUTH  = blockmatrix[Bcount].Mindex_2D(i, j, SOUTH);
	      M_ind_EAST   = blockmatrix[Bcount].Mindex_2D(i, j, EAST);
	      M_ind_WEST   = blockmatrix[Bcount].Mindex_2D(i, j, WEST);
	      
	      /* Calculate Coefficient X - SOUTH, and assign the values (4x4 block) into the corresponding location */
	      /* in the Jacobian matrix. */
	      if (M_ind_SOUTH.x >= 0 && M_ind_SOUTH.y >= 0) {
		dFdU = Jacobian_LocalBlock(i, j, Soln_ptr[Bcount], SOUTH, 
					   M_ind_SOUTH,
					   blockmatrix[Bcount].blocksize,
					   normalization);
		blockmatrix[Bcount].setblock(M_ind_SOUTH.x,M_ind_SOUTH.y, dFdU);
	      } /* endif */
	      
	      /* Calculate Coefficient B - WEST, and assign the values (4x4 block) into the corresponding location */
	      /* in the Jacobian matrix. */
	      if (M_ind_WEST.x >= 0 && M_ind_WEST.y >= 0) {
		dFdU = Jacobian_LocalBlock(i, j, Soln_ptr[Bcount], WEST, 
					   M_ind_WEST,
					   blockmatrix[Bcount].blocksize,
					   normalization);
		blockmatrix[Bcount].setblock(M_ind_WEST.x,M_ind_WEST.y, dFdU);
	      } /* endif */ 
	      
	      /* Calculate Coefficient C - CENTER, and assign the values (4x4 block) into the corresponding location */
	      /* in the Jacobian matrix. */
	      if (M_ind_CENTER.x >= 0 && M_ind_CENTER.y >= 0) {
		dFdU = Jacobian_LocalBlock(i, j, Soln_ptr[Bcount], CENTER, 
					   M_ind_CENTER,
					   blockmatrix[Bcount].blocksize, 
					   normalization, CFL_current*Soln_ptr[Bcount].dt[i][j]);
		blockmatrix[Bcount].setblock(M_ind_CENTER.x,M_ind_CENTER.y, dFdU);
	      } /* endif */
	      
	      /* Calculate Coefficient D - EAST, and assign the values (4x4 block) into the corresponding location */
	      /* in the Jacobian matrix. */
	      if (M_ind_EAST.x >= 0 && M_ind_EAST.y >= 0) {
		dFdU = Jacobian_LocalBlock(i, j, Soln_ptr[Bcount], EAST,
					   M_ind_EAST,
					   blockmatrix[Bcount].blocksize,
					   normalization);
		blockmatrix[Bcount].setblock(M_ind_EAST.x,M_ind_EAST.y, dFdU);
	      } /* endif */
	      
	      /* Calculate Coefficient Y - NORTH, and assign the values (4x4 block) into the corresponding location */
	      /* in the Jacobian matrix. */
	      if (M_ind_NORTH.x >= 0 && M_ind_NORTH.y >= 0) {
		dFdU = Jacobian_LocalBlock(i, j, Soln_ptr[Bcount], NORTH, 
					   M_ind_NORTH,
					   blockmatrix[Bcount].blocksize,
					   normalization);
		blockmatrix[Bcount].setblock(M_ind_NORTH.x,M_ind_NORTH.y, dFdU);
	      } /* endif */ 
	      
	    } /* endfor */
	  } /* endfor */
	  
	} /* endif */
      } /* endfor */
      } else {
	flag = OFF;
      } /* endif */ 

      for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
	if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {

	  /* Choose Preconditioner */
	  if (P_Switch == 1) {      
	    /* ILU(k) */
	    MBILUK[Bcount].localprecon(LP_INVERSE);
	    MBILUK[Bcount].setup(blockmatrix[Bcount], level_of_fill);

	  } else if (P_Switch == 2) {
	    /* Diagonal */
	    MBJacobi[Bcount].localprecon(LP_INVERSE);
	    MBJacobi[Bcount].setup(blockmatrix[Bcount]);

	  } /* endif */
	  
	} /* endif */
      } /* endfor */

      for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
	if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	  
	  /* Initialize G[Bcount].U[num] to zero for GMRES calculation. */
	  for (num=0;num<G[Bcount].scalar_dim;++num) G[Bcount].U[num] = 0.0;
	  
  	  /* Copy solutions in conserved variables, U to Uo for all blocks for update procedure. */
	  for (j = Soln_ptr[Bcount].JCl-Soln_ptr[Bcount].Nghost; 
	       j <= Soln_ptr[Bcount].JCu+Soln_ptr[Bcount].Nghost; j++){
	    for (i = Soln_ptr[Bcount].ICl-Soln_ptr[Bcount].Nghost; 
		 i <= Soln_ptr[Bcount].ICu+Soln_ptr[Bcount].Nghost; i++){
	      Soln_ptr[Bcount].Uo[i][j].d    = Soln_ptr[Bcount].U[i][j].d;
	      Soln_ptr[Bcount].Uo[i][j].dv.x = Soln_ptr[Bcount].U[i][j].dv.x;
	      Soln_ptr[Bcount].Uo[i][j].dv.y = Soln_ptr[Bcount].U[i][j].dv.y;
	      Soln_ptr[Bcount].Uo[i][j].E    = Soln_ptr[Bcount].U[i][j].E;

	      /* Copy Soln_ptr[Bcount].U to G[Bcount].U for GMRES calculation. */ 
	      G[Bcount].U[G[Bcount].index(i,j)]   = 
		Soln_ptr[Bcount].U[i][j].d;
	      G[Bcount].U[G[Bcount].index(i,j,1)] = 
		Soln_ptr[Bcount].U[i][j].dv.x;
	      G[Bcount].U[G[Bcount].index(i,j,2)] = 
		Soln_ptr[Bcount].U[i][j].dv.y;
	      G[Bcount].U[G[Bcount].index(i,j,3)] = 
		Soln_ptr[Bcount].U[i][j].E;

	      /* Set time step. */
              Soln_ptr[Bcount].dt[i][j] = CFL_current*Soln_ptr[Bcount].dt[i][j];

	    } /* endfor */
	  } /* endfor */
	  
	} /* endif */
      } /* endfor */

      /* Solve system with right-preconditioned matrix free GMRES */ 
      GMRES_Algorithm(Soln_ptr,
		      Soln_Block_List,
                      Input_Parameters,
		      G, gmrestol,
		      MBILUK, MBJacobi,
		      normalization,
                      Number_of_Newton_Steps);

      /* Update Solution. */
      int index = 0;
      for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
	if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	  for (j = Soln_ptr[Bcount].JCl-Soln_ptr[Bcount].Nghost; 
	       j <= Soln_ptr[Bcount].JCu+Soln_ptr[Bcount].Nghost; j++){
	    for (i = Soln_ptr[Bcount].ICl-Soln_ptr[Bcount].Nghost; 
		 i <= Soln_ptr[Bcount].ICu+Soln_ptr[Bcount].Nghost; i++){

	      /* Update solutions in conversed variables. */
	      Soln_ptr[Bcount].U[i][j].d    = Soln_ptr[Bcount].Uo[i][j].d
		+ G[Bcount].x[G[Bcount].index(i,j)];
	      Soln_ptr[Bcount].U[i][j].dv.x = Soln_ptr[Bcount].Uo[i][j].dv.x +
		G[Bcount].x[G[Bcount].index(i,j,1)];
	      Soln_ptr[Bcount].U[i][j].dv.y = Soln_ptr[Bcount].Uo[i][j].dv.y +
		G[Bcount].x[G[Bcount].index(i,j,2)];
	      Soln_ptr[Bcount].U[i][j].E    = Soln_ptr[Bcount].Uo[i][j].E
		+ G[Bcount].x[G[Bcount].index(i,j,3)];

	      /* Apply update reduction while any one of the updated variables is negative. */ 
	      if (Soln_ptr[Bcount].U[i][j].d   <= ZERO ||
                  Soln_ptr[Bcount].U[i][j].e() <= ZERO ||
                  Soln_ptr[Bcount].U[i][j].E   <= ZERO) { 
 	         update_reduction_factor = ONE;
	     
	         for (n_update_reduction = 1; n_update_reduction <= 10; ++n_update_reduction) {
	           update_reduction_factor = HALF*update_reduction_factor;
	           G[Bcount].x[G[Bcount].index(i,j)] *= update_reduction_factor;
	           G[Bcount].x[G[Bcount].index(i,j,1)] *= update_reduction_factor;
	           G[Bcount].x[G[Bcount].index(i,j,2)] *= update_reduction_factor;
	           G[Bcount].x[G[Bcount].index(i,j,3)] *= update_reduction_factor;
   	           Soln_ptr[Bcount].U[i][j].d    = Soln_ptr[Bcount].Uo[i][j].d +
		                                   G[Bcount].x[G[Bcount].index(i,j)];
	           Soln_ptr[Bcount].U[i][j].dv.x = Soln_ptr[Bcount].Uo[i][j].dv.x +
 		                                   G[Bcount].x[G[Bcount].index(i,j,1)];
	           Soln_ptr[Bcount].U[i][j].dv.y = Soln_ptr[Bcount].Uo[i][j].dv.y +
		                                   G[Bcount].x[G[Bcount].index(i,j,2)];
	           Soln_ptr[Bcount].U[i][j].E    = Soln_ptr[Bcount].Uo[i][j].E +
		                                   G[Bcount].x[G[Bcount].index(i,j,3)];
	           if (Soln_ptr[Bcount].U[i][j].d   > ZERO &&
		       Soln_ptr[Bcount].U[i][j].E   > ZERO &&
		       Soln_ptr[Bcount].U[i][j].e() > ZERO ) break;
	         } /* endfor */
	      } /* endif */

	      /* Print error: Negative Density or Negative Pressure. */ 
	      if (Soln_ptr[Bcount].U[i][j].d <= ZERO) { 
		cout << "\n NEGATIVE DENSITY : " << endl;
		cout << "Soln_ptr["<<Bcount<<"].U["<<i<<"]["<<j<<"].d = " 
		     << Soln_ptr[Bcount].U[i][j].d 
		     << ": Soln_ptr["<<Bcount<<"].Uo["<<i<<"]["<<j<<"].d = " 
		     << Soln_ptr[Bcount].Uo[i][j].d << endl;
		cout << "   G["<<Bcount<<"].x[G["<<Bcount<<"].index("<<i
		     <<","<<j
		     <<")] = " 
		     << G[Bcount].x[G[Bcount].index(i,j)] << endl;
	      }	else if (Soln_ptr[Bcount].U[i][j].e() <= ZERO) { 
		cout << "\n NEGATIVE INTERNAL ENERGY : " << endl;
		cout << "Soln_ptr["<<Bcount<<"].U["<<i<<"]["<<j<<"].e() = " 
		     << Soln_ptr[Bcount].U[i][j].e() << endl;
	      } /* endif */ 

	      /* Update solution in primitive variables. */
	      Soln_ptr[Bcount].W[i][j] = W(Soln_ptr[Bcount].U[i][j]);

	    } /* endfor */
	  } /* endfor */
	} /* endif */
      } /* endfor */

      /* Exchange solution information between neighbouring blocks. */
      error_flag = Send_All_Messages(Soln_ptr,
				     Soln_Block_List,
				     NUM_VAR_EULER2D, 
				     OFF);
      if (error_flag) {
         cout << "\n Euler2D_NKS ERROR: Euler2D message passing error on processor "
	      << Soln_Block_List.ThisCPU
	      << ".\n";
         cout.flush();
      } /* endif */
      error_flag = CFDkit_OR_MPI(error_flag);
      if (error_flag) return (error_flag);
	  
      /* Apply boundary conditions for Newton step. */
      BCs(Soln_ptr,
	  Soln_Block_List,
	  Input_Parameters);
      
      Number_of_Newton_Steps++;
    } else {
      stop = 0;
    } /* endif */
  } /* endwhile */
  
  /* Reset limiter. */
  Input_Parameters.i_Limiter = i_limiter;

  /* Calculate final residual with latest solutions for all blocks before termination. */
  for (Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
    if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
      
      dUdt_Residual_Evaluation(Soln_ptr[Bcount],
			       Input_Parameters);

      G[Bcount].copy_residual(G[Bcount].b, Soln_ptr[Bcount]);
      
      for (i=0; i<G[Bcount].scalar_dim; ++i) 
	G[Bcount].b[i] = -1.0 * G[Bcount].b[i];
      
    } /* endif */
  } /* endfor */

  /* Calculate 2-norm for all blocks */ 
  double L2norm, L2norm_n;   
  L2norm = sqr(L2_Norm_Residual(Soln_ptr, Soln_Block_List));
  L2norm = sqrt(CFDkit_Summation_MPI(L2norm));
  L2norm_n = L2norm / L2norm_first;

  if (CFDkit_Primary_MPI_Processor()) {     
    cout << " " << endl;
    for (int star=0;star<75;star++){cout <<"*";}
    cout.precision(20);
    cout << "\nEnd of Newton Steps = " << setw(3) << Number_of_Newton_Steps-1 
	 << "   L2norm = " << setw(12) << L2norm 
	 << endl;
    
    for (int star=0;star<75;star++){cout <<"*";}
  } /* endif */

  /* Update iteration counter. */
  Number_of_Startup_Interations = Number_of_Startup_Interations+Number_of_Newton_Steps-1;

  /* Deallocate memory. */
  delete [] MBILUK;
  delete [] MBJacobi;

  for ( Bcount = 0 ; Bcount < NBLK; ++Bcount ) {
    if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
      G[Bcount].deallocate();
    } /* endif */ 
  } /* endfor */
  delete [] G;

  hbtmat_object.deallocate_mem();
  delete [] blockmatrix;

  return 0;

} /* End of Newton_Krylov_Schwarz_Solver. */

/*********************************************************
 * Routine: Rotation_Matrix                              *
 *                                                       *
 * This function returns either the rotation matrix, A,  *
 * or the inverse of A.                                  *
 *                                                       *
 * Note: A_matrix = 1 for returning A.                   *
 *                = 0 for returning inverse of A.        *
 *                                                       *
 *********************************************************/
DenseMatrix Rotation_Matrix(Vector2D nface, int A_matrix) 
{
  double cos_angle = nface.x; 
  double sin_angle = nface.y;
    
  DenseMatrix mat(4,4);
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

/*********************************************************
 * Routine: Nface                                        *
 *                                                       *
 * This function returns the normal vector for the       *
 * specified face of the cell.                           *
 *                                                       *
 *********************************************************/
Vector2D Nface(Euler2D_Quad_Block SolnBlk, int xpt, int ypt, int location)
{
  Vector2D nface;

  switch(location) { 
  case NORTH:
    nface = SolnBlk.Grid.nfaceN(xpt,ypt);
    break;
  case SOUTH:
    nface = SolnBlk.Grid.nfaceS(xpt,ypt);
    break;
  case EAST:
    nface = SolnBlk.Grid.nfaceE(xpt,ypt);
    break;
  case WEST:
    nface = SolnBlk.Grid.nfaceW(xpt,ypt);
    break;
  default:
    cout << "Problem : Vector2D Nface......" << endl;
    break;
  } /* endif */

  return nface;

} /* End of Nface. */

/********************************************************
 * Routine: Jacobian_LocalBlock                         *
 *                                                      *
 * This routine returns Jacobian matrix for the         *
 * specified local solution block.                      *
 *                                                      *
 ********************************************************/
DenseMatrix Jacobian_LocalBlock(int x_pt, 
                                int y_pt, 
                                Euler2D_Quad_Block &SolnBlk, 
                                int location, 
                                Index M_ind, 
                                int blocksize, 
                                int normalization, 
                                double dTime) {

  /* Declaration. */
  Euler2D_pState W;
  DenseMatrix dFdU_temp(blocksize, blocksize, 0.0);
  DenseMatrix dFdU_(blocksize, blocksize, 0.0);
  DenseMatrix II(blocksize, blocksize, 0.0);
  DenseMatrix IM(blocksize, blocksize, 1.0);
  DenseMatrix A_N(blocksize, blocksize), AI_N(blocksize, blocksize);
  DenseMatrix A_S(blocksize, blocksize), AI_S(blocksize, blocksize);
  DenseMatrix A_E(blocksize, blocksize), AI_E(blocksize, blocksize);
  DenseMatrix A_W(blocksize, blocksize), AI_W(blocksize, blocksize);

  Vector2D lambdas_N, lambdas_S, lambdas_E, lambdas_W;
  Vector2D nface_N, nface_S, nface_E, nface_W;
  Index G_ind, G_ind_N, G_ind_S, G_ind_E, G_ind_W;

  double alpha_N, beta_N, gamma_N;
  double alpha_S, beta_S, gamma_S;
  double alpha_W, beta_W, gamma_W;
  double alpha_E, beta_E, gamma_E;
  double cos_angle, sin_angle;

  G_ind = Gindex_2D(x_pt, y_pt, location);
  int Fswitch = M_ind.y - M_ind.x; 
  int xpts = (SolnBlk.ICu-SolnBlk.ICl+1) + SolnBlk.Nghost * 2;
  int ypts = (SolnBlk.JCu-SolnBlk.JCl+1) + SolnBlk.Nghost * 2; 

  int a_xpts = SolnBlk.ICu-SolnBlk.ICl+1;
  int a_ypts = SolnBlk.JCu-SolnBlk.JCl+1;
  int overlap = 0;

  if ((G_ind.x >= a_xpts + SolnBlk.Nghost + overlap) ||
      (G_ind.x < SolnBlk.Nghost - overlap) || 
      (G_ind.y >= a_ypts + SolnBlk.Nghost + overlap) || 
      (G_ind.y < SolnBlk.Nghost - overlap)) 
    {
      /* GHOST CELL. */
      if (Fswitch == 0) {
	/* Identity Matrix -> Diagonal Blocks. */
	dFdU_.identity();
      } else {
	/* Zero -> Off-diagonal Blocks. */
	dFdU_.zero();
      } 
      return dFdU_;

    } else {

      /* NON-GHOST CELL. */
      II.identity();
      if (Fswitch == - xpts) {
	
	dFdU_.zero();
	
	/* Construct Coefficient X - SOUTH (4x4 block). */
	
	/* Obtain normal vector -> in Vector2D format. */
	nface_S = SolnBlk.Grid.nfaceS(x_pt, y_pt);

	/* Calculate wavespeeds using solutions in the rotated frame -> in Vector2D format. */
	lambdas_S = HLLE_wavespeeds(SolnBlk.W[x_pt][y_pt],
			            SolnBlk.W[G_ind.x][G_ind.y], 
			            nface_S);

	if ((lambdas_S.y-lambdas_S.x) == ZERO)
	  cout << " PROBLEM X - SOUTH :  HLLE_wavespeeds " << endl;

	/* Calculate constants gamma and beta -> scalar values. */
	gamma_S = (lambdas_S.x*lambdas_S.y)/(lambdas_S.y-lambdas_S.x);
	beta_S  = -1.0 * lambdas_S.x/(lambdas_S.y-lambdas_S.x);

	/* Obtain rotation matrices with normal vector -> 4x4 matrices in DenseMatrix format. */
	A_S  = Rotation_Matrix(nface_S, 1);
	AI_S = Rotation_Matrix(nface_S, 0);

	/* Calculate dFdU using solutions in the rotated frame -> 4x4 matrix in DenseMatrix format. */ 
	dFdU_temp.zero();
	dFdU(dFdU_temp, Rotate(SolnBlk.W[G_ind.x][G_ind.y], nface_S));

	/* Calculate Jacobian matrix -> 4x4 matrix in DenseMatrix format. */
	dFdU_ =  -1.0 * ((SolnBlk.Grid.lfaceS(x_pt, y_pt) * 
		AI_S * (beta_S * dFdU_temp + gamma_S * II) * A_S))
	        / SolnBlk.Grid.area(x_pt, y_pt);

	/* Apply normalization. */
	if (normalization == ON) normalize_Jacobian(dFdU_);

      } else if (Fswitch == -1) {

	dFdU_.zero();

	/* Construct Coefficient B - WEST (4x4 block). */

	/* Obtain normal vector -> in Vector2D format -> in Vector2D format. */
	nface_W = SolnBlk.Grid.nfaceW(x_pt, y_pt);

	/* Calculate wavespeeds using solutions in the rotated frame -> in Vector2D format. */
	lambdas_W = HLLE_wavespeeds(SolnBlk.W[x_pt][y_pt],
			            SolnBlk.W[G_ind.x][G_ind.y], 
			            nface_W);

	if ((lambdas_W.y-lambdas_W.x) == ZERO)
	  cout << " PROBLEM B - WEST : HLLE_wavespeeds " << endl;

	/* Calculate constants gamma and beta -> scalar values. */
	gamma_W = (lambdas_W.x*lambdas_W.y)/(lambdas_W.y-lambdas_W.x);
	beta_W  = -1.0 * lambdas_W.x/(lambdas_W.y-lambdas_W.x);

	/* Obtain rotation matrices with normal vector -> 4x4 matrices in DenseMatrix format. */
	A_W  = Rotation_Matrix(nface_W, 1);
	AI_W = Rotation_Matrix(nface_W, 0);
    
	/* Calculate dFdU using solutions in the rotated frame -> 4x4 matrix in DenseMatrix format. */
	dFdU_temp.zero();
	dFdU(dFdU_temp, Rotate(SolnBlk.W[G_ind.x][G_ind.y], nface_W));

	/* Calculate Jacobian matrix -> 4x4 matrix in DenseMatrix format. */
	dFdU_ =  -1.0 * (SolnBlk.Grid.lfaceW(x_pt, y_pt) * 
		 AI_W * (beta_W * dFdU_temp + gamma_W * II) * A_W)
	        / SolnBlk.Grid.area(x_pt, y_pt);

	/* Apply normalization. */
	if (normalization == ON) normalize_Jacobian(dFdU_);

      } else if (Fswitch == 0) {

	dFdU_.zero();
	II.identity();

	DenseMatrix dFdU_temp_N(blocksize, blocksize, 0.0),
	            dFdU_temp_S(blocksize, blocksize, 0.0),
	            dFdU_temp_E(blocksize, blocksize, 0.0),
	            dFdU_temp_W(blocksize, blocksize, 0.0);

	/* Construct Coefficient C - CENTER (4x4 block). */    		
	
	/* Obtain grid indices for all faces -> in Index format. */
	G_ind_N  = Gindex_2D(x_pt, y_pt, NORTH);
	G_ind_S  = Gindex_2D(x_pt, y_pt, SOUTH);
	G_ind_E  = Gindex_2D(x_pt, y_pt, EAST);
	G_ind_W  = Gindex_2D(x_pt, y_pt, WEST);
	
	/* Obtain all normal vectors -> in Vector2D format. */
	nface_N = SolnBlk.Grid.nfaceN(x_pt, y_pt);
	nface_S = SolnBlk.Grid.nfaceS(x_pt, y_pt);
	nface_E = SolnBlk.Grid.nfaceE(x_pt, y_pt);
	nface_W = SolnBlk.Grid.nfaceW(x_pt, y_pt);

	/* Obtain rotation matrices with differnt normal vectors -> 4x4 matrices in DenseMatrix format. */
	A_N  = Rotation_Matrix(nface_N, 1);
	AI_N = Rotation_Matrix(nface_N, 0);
	A_S  = Rotation_Matrix(nface_S, 1);
	AI_S = Rotation_Matrix(nface_S, 0);
	A_E  = Rotation_Matrix(nface_E, 1);
	AI_E = Rotation_Matrix(nface_E, 0);
	A_W  = Rotation_Matrix(nface_W, 1);
	AI_W = Rotation_Matrix(nface_W, 0);

	/* Calculate wavespeeds using solutions in the rotated frame -> 4x4 matrices in DenseMatrix format. */
	lambdas_N = HLLE_wavespeeds(SolnBlk.W[x_pt][y_pt], 
			            SolnBlk.W[G_ind_N.x][G_ind_N.y], 
			            nface_N);
	lambdas_S = HLLE_wavespeeds(SolnBlk.W[x_pt][y_pt],
			            SolnBlk.W[G_ind_S.x][G_ind_S.y], 
			            nface_S);
	lambdas_E = HLLE_wavespeeds(SolnBlk.W[x_pt][y_pt], 
			            SolnBlk.W[G_ind_E.x][G_ind_E.y], 
			            nface_E);
	lambdas_W = HLLE_wavespeeds(SolnBlk.W[x_pt][y_pt],
			            SolnBlk.W[G_ind_W.x][G_ind_W.y],
			            nface_W);

	if ((lambdas_N.y-lambdas_N.x) == ZERO)
	  cout << " PROBLEM C - CENTER->NORTH : HLLE_wavespeeds " << endl;
	if ((lambdas_S.y-lambdas_S.x) == ZERO)
	  cout << " PROBLEM C - CENTER->SOUTH : HLLE_wavespeeds " << endl;
	if ((lambdas_E.y-lambdas_E.x) == ZERO)
	  cout << " PROBLEM C - CENTER->EAST : HLLE_wavespeeds " << endl;
	if ((lambdas_W.y-lambdas_W.x) == ZERO)
	  cout << " PROBLEM C - CENTER->WEST : HLLE_wavespeeds " << endl;

	/* Calculate constants alpha and gamma -> scalar values. */
	alpha_W = lambdas_W.y/(lambdas_W.y-lambdas_W.x);
	alpha_S = lambdas_S.y/(lambdas_S.y-lambdas_S.x);
	alpha_N = lambdas_N.y/(lambdas_N.y-lambdas_N.x);
	alpha_E = lambdas_E.y/(lambdas_E.y-lambdas_E.x);

	gamma_N = (lambdas_N.x*lambdas_N.y)/(lambdas_N.y-lambdas_N.x);
	gamma_S = (lambdas_S.x*lambdas_S.y)/(lambdas_S.y-lambdas_S.x);
	gamma_E = (lambdas_E.x*lambdas_E.y)/(lambdas_E.y-lambdas_E.x);
	gamma_W = (lambdas_W.x*lambdas_W.y)/(lambdas_W.y-lambdas_W.x);

	/* Calculate four dFdU using solutions in the rotated frame -> 4x4 matrices in DenseMatrix format. */
	dFdU_temp_N.zero();
	dFdU(dFdU_temp_N, Rotate(SolnBlk.W[x_pt][y_pt], nface_N));
	dFdU_temp_S.zero();
	dFdU(dFdU_temp_S, Rotate(SolnBlk.W[x_pt][y_pt], nface_S));
	dFdU_temp_E.zero();
	dFdU(dFdU_temp_E, Rotate(SolnBlk.W[x_pt][y_pt], nface_E));
	dFdU_temp_W.zero();
	dFdU(dFdU_temp_W, Rotate(SolnBlk.W[x_pt][y_pt], nface_W));

	/* Calculate Jacobian matrix -> 4x4 matrix in DenseMatrix format. */
	dFdU_ = -1.0 * ((SolnBlk.Grid.lfaceN(x_pt, y_pt) * 
			 AI_N * (alpha_N * dFdU_temp_N -  gamma_N * II) * A_N)
		      + (SolnBlk.Grid.lfaceS(x_pt, y_pt) * 
		         AI_S * (alpha_S * dFdU_temp_S -  gamma_S * II) * A_S)
		      + (SolnBlk.Grid.lfaceE(x_pt, y_pt) * 
		         AI_E * (alpha_E * dFdU_temp_E -  gamma_E * II) * A_E)
		      + (SolnBlk.Grid.lfaceW(x_pt, y_pt) * 
		         AI_W * (alpha_W * dFdU_temp_W -  gamma_W * II) * A_W)
		       ) / SolnBlk.Grid.area(x_pt, y_pt);

	/* Apply Finite Time Step, dTime, on Diagonal Blocks Only -> 4x4 matrix in DenseMatrix format. */
	dFdU_ = dFdU_ - II * (1.0/dTime);

	/* Apply normalization. */
	if (normalization == ON) normalize_Jacobian(dFdU_);

      } else if (Fswitch == 1) {

 	dFdU_.zero();

	/* Construct Coefficient D - EAST (4x4 block). */

	/* Obtain normal vector -> in Vector2D format. */
	nface_E = SolnBlk.Grid.nfaceE(x_pt, y_pt);

	/* Calculate wavespeeds using solutions in the rotated frame -> 4x4 matrix in DenseMatrix format. */
	lambdas_E = HLLE_wavespeeds(SolnBlk.W[x_pt][y_pt], 
			            SolnBlk.W[G_ind.x][G_ind.y], nface_E);

	if ((lambdas_E.y-lambdas_E.x) == ZERO)
	  cout << " PROBLEM D - EAST : HLLE_wavespeeds " << endl;

	/* Calculate constants gamma and beta -> scalar values. */
	beta_E  = -1.0 * lambdas_E.x/(lambdas_E.y-lambdas_E.x);
	gamma_E = (lambdas_E.x*lambdas_E.y)/(lambdas_E.y-lambdas_E.x);

	/* Obtain rotation matrices with normal vector -> 4x4 matrices in DenseMatrix format. */
	A_E  = Rotation_Matrix(nface_E, 1);
	AI_E = Rotation_Matrix(nface_E, 0);
    
	/* Calculate dFdU using solutions in the rotated frame -> 4x4 matrix in DenseMatrix format. */
        dFdU_temp.zero();
	dFdU(dFdU_temp, Rotate(SolnBlk.W[G_ind.x][G_ind.y], nface_E));

	/* Calculate Jacobian matrix -> 4x4 matrix in DenseMatrix format. */
	dFdU_ = -1.0 * (SolnBlk.Grid.lfaceE(x_pt, y_pt) * 
			AI_E * (beta_E * dFdU_temp +  gamma_E * II) * A_E)
	        / SolnBlk.Grid.area(x_pt, y_pt);

	/* Apply normalization. */
	if (normalization == ON) normalize_Jacobian(dFdU_);

      } else if (Fswitch == xpts) {

	dFdU_.zero();

	/* Construct Coefficient Y - NORTH (4x4 block). */

	/* Obtain normal vector -> in Vector2D format. */
	nface_N = SolnBlk.Grid.nfaceN(x_pt, y_pt);

	/* Calculate wavespeeds using solutions in the rotated frame -> 4x4 matrix in DenseMatrix format. */
	lambdas_N = HLLE_wavespeeds(SolnBlk.W[x_pt][y_pt], 
			            SolnBlk.W[G_ind.x][G_ind.y], nface_N);

	if ((lambdas_N.y-lambdas_N.x) == ZERO)
	  cout << " PROBLEM N - NORTH : HLLE_wavespeeds " << endl;

	/* Calculate constants gamma and beta -> scalar values. */
	beta_N  = -1.0 * lambdas_N.x/(lambdas_N.y-lambdas_N.x);
	gamma_N = (lambdas_N.x*lambdas_N.y)/(lambdas_N.y-lambdas_N.x);

	/* Obtain rotation matrices with normal vector -> 4x4 matrices in DenseMatrix format. */
	A_N  = Rotation_Matrix(nface_N, 1);
	AI_N = Rotation_Matrix(nface_N, 0);

	/* Calculate dFdU using solutions in the rotated frame -> 4x4 matrix in DenseMatrix format.  */
        dFdU_temp.zero();
	dFdU(dFdU_temp, Rotate(SolnBlk.W[G_ind.x][G_ind.y], nface_N));

	/* Calculate Jacobian matrix -> 4x4 matrix in DenseMatrix format. */
	dFdU_ = -1.0 * (SolnBlk.Grid.lfaceN(x_pt, y_pt) * 
			AI_N * (beta_N * dFdU_temp +  gamma_N * II) * A_N)
	        / SolnBlk.Grid.area(x_pt, y_pt);
	
	/* Apply normalization. */
	if (normalization == ON) normalize_Jacobian(dFdU_);

      } else { 
	cout << "\n**Problem : Jacobian Calculation!! \n\n" << endl;
      } /* endif */

      return dFdU_;
    } /* endif */

} /* End of Jacobian_LocalBlock. */

/********************************************************
 * Routine: Create_Block_Matrix                         *
 *                                                      *
 * This routine returns the model structure in HBTMat   *
 * (Harwell Boeing) fomat for a specified blocksize.    *
 *                                                      *
 ********************************************************/ 
HBTMat Create_Block_Matrix(int nx, int ny, int nz, HBTMat_Class &HC) {
  HBTMat model; 
  if (nz == 1) {
    if (HC.fill_pattern == 1) {
      model = Create_BlockMat_Object_5 (nx, ny, HC); 
    } else { 
      model = Create_BlockMat_Object_13(nx, ny, HC);
    }  
  } else { 
    cout << " Euler3D : NOT READY !" << endl; 
  } 

  return model;

} /* End of Create_Block_Matrix. */

/********************************************************
 * Routine: Create_BlockMat_Object_5                    *
 *                                                      *
 * This routine returns a data structure of a           *
 * five-coefficient matrix in Harwell-Boeing matrix     *
 * format using the nearest neighboring cells.          *
 ********************************************************/ 
HBTMat Create_BlockMat_Object_5(int xpts, int ypts, HBTMat_Class &HC) {

  /* Local variable definitions. */
  int blocksize = HC.blocksize;

  double ****BSM = new double ***[xpts*ypts];
  for (int i = 0; i < xpts*ypts; i++) {
    BSM[i] = new double **[5];
  }
  for (int i = 0; i < xpts*ypts; i++) {
    for (int j = 0; j < 5; j++) {
      BSM[i][j] = new double *[blocksize];
    }
  }
  for (int i = 0; i < xpts*ypts; i++) {
    for (int j = 0; j < 5; j++) {
      for (int k = 0; k < blocksize; k++) {
	BSM[i][j][k] = new double [blocksize];
      }
    }
  }

  short int ** Sx    = new short int * [xpts*ypts];
  short int ** Sy    = new short int * [xpts*ypts];
  short int ** model = new short int * [xpts*ypts];
  for (int i = 0; i < xpts*ypts; i++) {
      Sx[i]= new short int [5];
      Sy[i]= new short int [5];
      model[i]= new short int [5];
  }

  int i,j,k;
  int ncol   = 0;
  int snz    = 0;
  int count  = 0;
  int rownum = 1;

  /* Initialization. */
 for (int i=0; i<xpts*ypts; i++){
    for (int j=0; j<5; j++){
     Sx[i][j] = 0;
     Sy[i][j] = 0;
     model[i][j] = 0;
    }
  }

  for (int i=0; i<xpts*ypts; i++){
    for (int j=0; j<5; j++){
      for (int n=0; n<blocksize; n++){
	for (int m=0; m<blocksize; m++){
	  BSM[i][j][n][m] = 0.0;
	}
      }
    }
  }

  /* Create Sx - a 2D vector with a size of (xpts*ypts x 4) modeling a 3-point central-differencing point operator */
  /* for a second dervivative in x-direction. */ 
  Sx[0][2]=-2;
  Sx[0][3]=1;
  for (int i=1; i<xpts*ypts-1; i++){
    Sx[i][1]=1;
    Sx[i][2]=-2;
    Sx[i][3]=1;}
  Sx[xpts*ypts-1][1]=1;
  Sx[xpts*ypts-1][2]=-2;

  for (int i=xpts; i<xpts*ypts; i+=xpts){
    Sx[i][1]=0;
    Sx[i-1][3]=0;}

  /* Create Sy - a 2D vector with a size of (xpts*ypts x 4) modeling a 3-point central-differencing point operator */
  /* for a second dervivative in y-direction. */   
  for (int i=0; i<xpts; i++){
    Sy[i][2]=-2;
    Sy[i][4]=1;}
  for (int i=xpts; i<xpts*ypts-xpts; i++){
    Sy[i][0]=1;
    Sy[i][2]=-2;
    Sy[i][4]=1;}
  for (int i=xpts*ypts-xpts; i<xpts*ypts; i++){
    Sy[i][0]=1;
    Sy[i][2]=-2;}

  /* Create SM = Sx+Sy. */
  for (int i=0; i<xpts*ypts; i++){
    for (int j=0; j<5; j++){
      model[i][j]=Sx[i][j]+Sy[i][j];
      if(model[i][j] != 0) snz++;
    }
  }
  
  /* Create HB_SM -  a 4D vector with a size of (xpts*ypts x 4 x 4 x 4), extending all entries in SM from a scalar to a 4x4 block. */
  for (int i=0; i<xpts*ypts; i++){
    for (int j=0; j<5; j++){
      if (model[i][j]!=0){
	BSM[i][j][0][0]=model[i][j];
	BSM[i][j][0][1]=model[i][j];
	BSM[i][j][0][2]=model[i][j];
	BSM[i][j][0][3]=model[i][j];
	BSM[i][j][1][0]=model[i][j];
	BSM[i][j][1][1]=model[i][j];
	BSM[i][j][1][2]=model[i][j];
	BSM[i][j][1][3]=model[i][j];
	BSM[i][j][2][0]=model[i][j];
	BSM[i][j][2][1]=model[i][j];
	BSM[i][j][2][2]=model[i][j];
	BSM[i][j][2][3]=model[i][j];
	BSM[i][j][3][0]=model[i][j];
	BSM[i][j][3][1]=model[i][j];
	BSM[i][j][3][2]=model[i][j];
	BSM[i][j][3][3]=model[i][j];
      }
    }
  }

  /* Create vectors: valv, ncolv, rowind, rhs, guess, exact for HBTMat. */
  HC.ncolv   = new int[snz * blocksize * blocksize];
  HC.nrowv   = new int[snz * blocksize * blocksize];
  HC.rowind  = new int[xpts * ypts * blocksize + 1];
  
  HC.valv  = new double[snz * blocksize * blocksize];
  HC.rhs   = new double[xpts * ypts * blocksize * blocksize];
  HC.guess = new double[xpts * ypts * blocksize * blocksize];
  HC.exact = new double[xpts * ypts * blocksize * blocksize];
  
  for (int i=0; i<xpts*ypts; i++){
    for (int p=0; p<blocksize; p++){    
      for (int j=0; j<5; j++){
	for (int q=0; q<blocksize; q++){
	  if (BSM[i][j][p][q]!=0){
	    if (j==0){
	      ncol = (i-xpts)*blocksize+q;}
	    else if (j==1){
	      ncol = (i-1)*blocksize+q;}
	    else if (j==2){
	      ncol = i*blocksize+q;}
	    else if (j==3){
	      ncol = (i+1)*blocksize+q;}
	    else if (j==4){
	      ncol = (i+xpts)*blocksize+q;}
	    else {cout << "Create_BlockMat_Object_5: PROBLEM!!!!!!" << endl;}
    
	    HC.ncolv[count] = ncol;
	    HC.nrowv[count] = i*blocksize+p;
	    HC.valv[count] = BSM[i][j][p][q];
	    count ++;
	  }}}}}

  HC.rowind[0]=0;
  for (int i=0; i<=count-2; i++){
    if (HC.nrowv[i] != HC.nrowv[i+1]) {
      HC.rowind[rownum] = i+1;
      rownum++;
    }
  }
  HC.rowind[rownum]=count;

  for (int i=0; i<xpts*ypts*blocksize; i++) {
    HC.rhs[i] = 0.0;
    HC.guess[i] = 0.0;
    HC.exact[i] = 0.0;
  }

  HBTMat example(xpts*ypts*blocksize, HC.valv, HC.ncolv, HC.rowind, 
		 HC.rhs, HC.guess, HC.exact);

  for (int i = 0; i < xpts*ypts; i++) {
      delete [] model[i];
  }
  delete [] model;

  for (int i = 0; i < xpts*ypts; i++) {
    for (int j = 0; j < 5; j++) {
      for (int k = 0; k < blocksize; k++) {
	delete [] BSM[i][j][k];
      }
    }
  }
  for (int i = 0; i < xpts*ypts; i++) {
    for (int j = 0; j < 5; j++) {
      delete [] BSM[i][j];
 
   }
  }
  for (int i = 0; i < xpts*ypts; i++) {
    delete [] BSM[i];
  }

  delete [] BSM;
  
  return example;

} /* End of Create_BlockMat_Object_5. */

/********************************************************
 * Routine: Create_BlockMat_Object_13                   *
 *                                                      *
 * This routine returns a data structure of a           *
 * thirteen-coefficient matrix in Harwell-              *
 * Boeing matrix format using the nearest and the       *
 * second nearest neighboring cells .                   *
 ********************************************************/ 
HBTMat Create_BlockMat_Object_13(int xpts, int ypts, HBTMat_Class &HC) {
  int blocksize = HC.blocksize;

  double ****BSM = new double ***[xpts*ypts];
  for (int i = 0; i < xpts*ypts; i++) {
    BSM[i] = new double **[13];
  } /* endfor */
  for (int i = 0; i < xpts*ypts; i++) {
    for (int j = 0; j < 13; j++) {
      BSM[i][j] = new double *[blocksize];
    } /* endfor */
  } /* endfor */
  for (int i = 0; i < xpts*ypts; i++) {
    for (int j = 0; j < 13; j++) {
      for (int k = 0; k < blocksize; k++) {
	BSM[i][j][k] = new double [blocksize];
      } /* endfor */
    } /* endfor */
  } /* endfor */

  short int ** Sx    = new short int * [xpts*ypts];
  short int ** Sy    = new short int * [xpts*ypts];
  short int ** model = new short int * [xpts*ypts];
  for (int i = 0; i < xpts*ypts; i++) {
      Sx[i]= new short int [13];
      Sy[i]= new short int [13];
      model[i]= new short int [13];
  } /* endfor */
  
  int i,j,k;
  int ncol   = 0;
  int snz    = 0;
  int count  = 0;
  int rownum = 1;

  // Initialization.
 for (int i=0; i<xpts*ypts; i++){
    for (int j=0; j<13; j++){
     Sx[i][j] = 0;
     Sy[i][j] = 0;
     model[i][j] = 0;
    } /* endfor */
  } /* endfor */

  for (int i=0; i<xpts*ypts; i++){
    for (int j=0; j<13; j++){
      for (int n=0; n<blocksize; n++){
	for (int m=0; m<blocksize; m++){
	  BSM[i][j][n][m] = 0.0;
	} /* endfor */
      } /* endfor */
    } /* endfor */
  } /* endfor */

  /**************************
   * Construct Model Matrix *
   **************************/
  for (i = 0; i < xpts * ypts; i++) {
    for (j = 0; j < 13 ; j++) {
      model[i][j] = 0;
    } /* endfor */
  } /* endfor */

  for ( i = 0; i < xpts*ypts; i++){
    if ( i >= xpts * 2) {
      model[i][0]=1; 
    } /* endif */
    if ( i >= xpts){
      if ( i % xpts != 0){
	model[i][1]=2;}
      model[i][2]=3;
      if ( i % xpts != xpts - 1){
	model[i][3]=4;}
    } /* endif */
    if (( i % xpts != 0) & ( i % xpts != 1)){
      model[i][4]=5;}
    if ( i % xpts != 0){
      model[i][5]=6;}
    model[i][6]=7;
    if ( i % xpts != xpts - 1){
      model[i][7]=8;}
    if (( i % xpts != xpts - 1) && ( i % xpts != xpts - 2)){
      model[i][8]=9;}
    if( i < xpts * ypts - xpts){
      if ( i % xpts != 0){
	model[i][9]=10;}
      model[i][10]=11;
      if ( i % xpts != xpts - 1){
	model[i][11]=12;}
    } /* endif */
    if ( i < xpts * ypts - xpts * 2){
      model[i][12]=13;}
  } /* endfor */

  //Create HB_SM
  for (int i=0; i<xpts*ypts; i++){
    for (int j=0; j<13; j++){
      if (model[i][j] != 0){
	snz++;
	BSM[i][j][0][0]=model[i][j];
	BSM[i][j][0][1]=model[i][j];
	BSM[i][j][0][2]=model[i][j];
	BSM[i][j][0][3]=model[i][j];
	BSM[i][j][1][0]=model[i][j];
	BSM[i][j][1][1]=model[i][j];
	BSM[i][j][1][2]=model[i][j];
	BSM[i][j][1][3]=model[i][j];
	BSM[i][j][2][0]=model[i][j];
	BSM[i][j][2][1]=model[i][j];
	BSM[i][j][2][2]=model[i][j];
	BSM[i][j][2][3]=model[i][j];
	BSM[i][j][3][0]=model[i][j];
	BSM[i][j][3][1]=model[i][j];
	BSM[i][j][3][2]=model[i][j];
	BSM[i][j][3][3]=model[i][j];
      } /* endif */
    } /* endfor */
  } /* endfor */

  // Create vectors: valv, ncolv, rowind, rhs, guess, exact for HBTMat.
  HC.ncolv   = new int[snz*blocksize];
  HC.nrowv   = new int[snz*blocksize];
  HC.rowind  = new int[xpts*ypts*blocksize+1];
  
  HC.valv  = new double[snz*blocksize];
  HC.rhs   = new double[xpts*ypts*blocksize];
  HC.guess = new double[xpts*ypts*blocksize];
  HC.exact = new double[xpts*ypts*blocksize];
  
  for (int i=0; i<xpts*ypts; i++){
    for (int p=0; p<blocksize; p++){    
      for (int j=0; j<13; j++){
	for (int q=0; q<blocksize; q++){
	  if (BSM[i][j][p][q] != 0){
	    if (j==0) {
	      ncol = (i-xpts*2)*blocksize+q;}
	    else if (j==1){
	      ncol = (i-xpts-1)*blocksize+q;}
	    else if (j==2){
	      ncol = (i-xpts)*blocksize+q;}
	    else if (j==3){
	      ncol = (i-xpts+1)*blocksize+q;}
	    else if (j==4){
	      ncol = (i-2)*blocksize+q;}
	    else if (j==5){
	      ncol = (i-1)*blocksize+q;}
	    else if (j==6){
	      ncol = i*blocksize+q;}
	    else if (j==7){
	      ncol = (i+1)*blocksize+q;}
	    else if (j==8){
	      ncol = (i+2)*blocksize+q;}
	    else if (j==9){
	      ncol = (i+xpts-1)*blocksize+q;}
	    else if (j==10){
	      ncol = (i+xpts)*blocksize+q;}
	    else if (j==11){
	      ncol = (i+xpts+1)*blocksize+q;}
	    else if (j==12){
	      ncol = (i+xpts*2)*blocksize+q;}
	    else {cout << " Create_BlockMat_Object_13: PROBLEM !!!!!!" << endl;
	    } /* endif */
	    HC.ncolv[count] = ncol;
	    HC.nrowv[count] = i*blocksize+p;
	    HC.valv[count] = BSM[i][j][p][q];
	    count ++;
	  } /* endif */
	} /* endfor */
      } /* endfor */
    } /* endfor */
  } /* endfor */

  cout << " count = " << count << endl;
  
  HC.rowind[0]=0;
  for (int i=0; i<=count-2; i++) {
    if (HC.nrowv[i] != HC.nrowv[i+1]) {
      HC.rowind[rownum] = i+1;
      rownum++;
    } /* endif */
  } /* endfor */
  HC.rowind[rownum]=count;

  for (int i=0; i<xpts*ypts*blocksize; i++) {
    HC.rhs[i] = 0.0;
    HC.guess[i] = 0.0;
    HC.exact[i] = 0.0;
  } /* endfor */

  for (int i=0; i<count; i++) {
    cout << " ncol["<<i<<"]= "<< HC.ncolv[i] 
	 << " nrow["<<i<<"]= "<< HC.nrowv[i] 
	 << " valv["<<i<<"]= "<< HC.valv[i] 
	 << endl;
  } /* endfor */
  for (int i=0; i<xpts*ypts*blocksize+1; i++) {
    cout << " rowind["<<i<<"]= "<< HC.rowind[i] << endl;
  } /* endfor */

  HBTMat example(xpts*ypts*blocksize, HC.valv, HC.ncolv, HC.rowind,
		 HC.rhs, HC.guess, HC.exact);

  for (int i = 0; i < xpts*ypts; i++) {
      delete [] model[i];
  } /* endfor */
  delete [] model;

  for (int i = 0; i < xpts*ypts; i++) {
    for (int j = 0; j < 13; j++) {
      for (int k = 0; k < blocksize; k++) {
	delete [] BSM[i][j][k];
      } /* endfor */
    } /* endfor */
  } /* endfor */
  for (int i = 0; i < xpts*ypts; i++) {
    for (int j = 0; j < 13; j++) {
      delete [] BSM[i][j];
    } /* endfor */
  } /* endfor */
  for (int i = 0; i < xpts*ypts; i++) {
    delete [] BSM[i];
  } /* endfor */

  delete [] BSM;
  
  return example;

} /* End of Create_BlockMat_Object_13. */

/********************************************************
 * Routine: normalize_Jacobian                          *
 *                                                      *
 * This routine returns the normalized Jacobian matrix. *
 *                                                      *
 ********************************************************/
void normalize_Jacobian(DenseMatrix &dFdU) {

  double ao  = Euler2D_W_STDATM.a();

  dFdU(0,0) = dFdU(0,0) * (ONE/ao);
  //dFdU(0,1) = dFdU(0,1);
  //dFdU(0,2) = dFdU(0,2);
  dFdU(0,3) = dFdU(0,3) * ao;

  dFdU(1,0) = dFdU(1,0) * (ONE/sqr(ao));
  dFdU(1,1) = dFdU(1,1) * (ONE/ao);
  dFdU(1,2) = dFdU(1,2) * (ONE/ao);
  //dFdU(1,3) = dFdU(1,3);

  dFdU(2,0) = dFdU(2,0) * (ONE/sqr(ao));
  dFdU(2,1) = dFdU(2,1) * (ONE/ao);
  dFdU(2,2) = dFdU(2,2) * (ONE/ao);
  //dFdU(2,3) = dFdU(2,3);

  dFdU(3,0) = dFdU(3,0) * (ONE/cube(ao));
  dFdU(3,1) = dFdU(3,1) * (ONE/sqr(ao));
  dFdU(3,2) = dFdU(3,2) * (ONE/sqr(ao));
  dFdU(3,3) = dFdU(3,3) * (ONE/ao);

}/* End of normalize_RHS. */

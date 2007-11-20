/* Gaussian2DQuadSingleBlock.cc:  Single-Block Versions of IO Subroutines for 2D Gaussian
                                  Multi-Block Quadrilateral Mesh 
                                  Solution Classes. */

/* Include 2D Gaussian quadrilateral mesh solution header file. */

#ifndef _GAUSSIAN2D_QUAD_INCLUDED
#include "Gaussian2DQuad.h"
#endif // _GAUSSIAN2D_QUAD_INCLUDED

/**************************************************************************
 * Gaussian2D_Quad_Block -- Single Block External Subroutines.            *
 **************************************************************************/

/********************************************************
 * Routine: Write_Solution_Block                        *
 *                                                      *
 * Writes the cell centred solution values of the       *
 * specified quadrilateral solution block to the        *
 * specified output stream for restart purposes.        *
 *                                                      *
 ********************************************************/
void Write_Solution_Block(Gaussian2D_Quad_Block &SolnBlk,
	                  ostream &Out_File) {

    Out_File << setprecision(14) << SolnBlk << setprecision(6);

}

/********************************************************
 * Routine: Read_Solution_Block                         *
 *                                                      *
 * Reads the cell centred solution values for the       *
 * specified quadrilateral solution block from the      *
 * specified input stream as required for restart       *
 * purposes.                                            *
 *                                                      *
 ********************************************************/
void Read_Solution_Block(Gaussian2D_Quad_Block &SolnBlk,
	                 istream &In_File) {

    In_File >> SolnBlk;

}


/********************************************************
 * Routine: Output_Tecplot                              *
 *                                                      *
 * Writes the solution values at the nodes of the       *
 * specified quadrilateral solution block to the        *
 * specified output stream suitable for plotting with   *
 * TECPLOT.                                             *
 *                                                      *
 ********************************************************/
void Output_Tecplot(Gaussian2D_Quad_Block &SolnBlk,
		    Gaussian2D_Input_Parameters &IP,
                    const int Number_of_Time_Steps,
                    const double &Time,
                    const int Block_Number,
                    const int Output_Title,
	            ostream &Out_File) {
  Output_Tecplot(SolnBlk, Number_of_Time_Steps, Time, Block_Number, Output_Title, Out_File);
}

void Output_Tecplot(Gaussian2D_Quad_Block &SolnBlk,
                    const int Number_of_Time_Steps,
                    const double &Time,
                    const int Block_Number,
                    const int Output_Title,
	            ostream &Out_File) {

    int i, j;
    Gaussian2D_pState W_node;
    Vector2D dX;

    /* Ensure boundary conditions are updated before
       evaluating solution at the nodes. Linear Extrapolation
       is used to ensure ouput is not affected by "strange"
       ghostcell information */

    //left and right
    for ( j = SolnBlk.Grid.JCl; j <= SolnBlk.Grid.JCu; j++){
      //left
      if(SolnBlk.Grid.BCtypeW[j] != BC_NONE) {
	Linear_Reconstruction_LeastSquares(SolnBlk,
					   SolnBlk.ICl, j, 
					   LIMITER_BARTH_JESPERSEN);
	dX = SolnBlk.Grid.Cell[SolnBlk.ICl-1][j].Xc -
	  SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc;
	SolnBlk.W[SolnBlk.ICl-1][j] = SolnBlk.W[SolnBlk.ICl][j] +
	  (SolnBlk.phi[SolnBlk.ICl][j]^SolnBlk.dWdx[SolnBlk.ICl][j])*dX.x +
	  (SolnBlk.phi[SolnBlk.ICl][j]^SolnBlk.dWdy[SolnBlk.ICl][j])*dX.y;
	SolnBlk.U[SolnBlk.ICl-1][j] = U(SolnBlk.W[SolnBlk.ICl-1][j]);
      }
    //right
      if(SolnBlk.Grid.BCtypeE[j] != BC_NONE) {
	Linear_Reconstruction_LeastSquares(SolnBlk,
					   SolnBlk.ICu, j,
					   LIMITER_BARTH_JESPERSEN);
	dX = SolnBlk.Grid.Cell[SolnBlk.ICu+1][j].Xc -
	  SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc;
	SolnBlk.W[SolnBlk.ICu+1][j] = SolnBlk.W[SolnBlk.ICu][j] +
	  (SolnBlk.phi[SolnBlk.ICu][j]^SolnBlk.dWdx[SolnBlk.ICu][j])*dX.x +
	  (SolnBlk.phi[SolnBlk.ICu][j]^SolnBlk.dWdy[SolnBlk.ICu][j])*dX.y;
	SolnBlk.U[SolnBlk.ICu+1][j] = U(SolnBlk.W[SolnBlk.ICu+1][j]);
      }
    }

    //bottom and top
    for ( i = SolnBlk.Grid.ICl; i <= SolnBlk.Grid.ICu; i++){
      //bottom
      if(SolnBlk.Grid.BCtypeS[i] != BC_NONE) {
	Linear_Reconstruction_LeastSquares(SolnBlk,
					   i, SolnBlk.JCl,
					   LIMITER_BARTH_JESPERSEN);
	dX = SolnBlk.Grid.Cell[i][SolnBlk.JCl-1].Xc -
	  SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc;
	SolnBlk.W[i][SolnBlk.JCl-1] = SolnBlk.W[i][SolnBlk.JCl] +
	  (SolnBlk.phi[i][SolnBlk.JCl]^SolnBlk.dWdx[i][SolnBlk.JCl])*dX.x +
	  (SolnBlk.phi[i][SolnBlk.JCl]^SolnBlk.dWdy[i][SolnBlk.JCl])*dX.y;
	SolnBlk.U[i][SolnBlk.JCl-1] = U(SolnBlk.W[i][SolnBlk.JCl-1]);
      }
      //top
      if(SolnBlk.Grid.BCtypeN[i] != BC_NONE) {
	Linear_Reconstruction_LeastSquares(SolnBlk,
					   i, SolnBlk.JCu,
					   LIMITER_BARTH_JESPERSEN);
	dX = SolnBlk.Grid.Cell[i][SolnBlk.JCu+1].Xc -
	  SolnBlk.Grid.Cell[i][SolnBlk.JCu].Xc;
	SolnBlk.W[i][SolnBlk.JCu+1] = SolnBlk.W[i][SolnBlk.JCu] +
	  (SolnBlk.phi[i][SolnBlk.JCu]^SolnBlk.dWdx[i][SolnBlk.JCu])*dX.x +
	  (SolnBlk.phi[i][SolnBlk.JCu]^SolnBlk.dWdy[i][SolnBlk.JCu])*dX.y;
	SolnBlk.U[i][SolnBlk.JCu+1] = U(SolnBlk.W[i][SolnBlk.JCu+1]);
      }
    }

    /* Output node solution data. */

    Out_File << setprecision(14);
    if (Output_Title) {
       Out_File << "TITLE = \"" << CFFC_Name() << ": 2D Gaussian Solution, "
                << "Time Step/Iteration Level = " << Number_of_Time_Steps
                << ", Time = " << Time
                << "\"" << "\n"
	        << "VARIABLES = \"x\" \\ \n"
                << "\"y\" \\ \n"
                << "\"rho\" \\ \n"
                << "\"u\" \\ \n"
                << "\"v\" \\ \n"
                << "\"pxx\" \\ \n"
                << "\"pxy\" \\ \n"
                << "\"pyy\" \\ \n"
                << "\"pzz\" \\ \n"
                << "\"erot\" \\ \n"
		<< "\"v_abs\" \\ \n"
		<< "\"M\" \\ \n"
                << "\"T\" \\ \n"
                << "\"DetP\" \n"
                << "ZONE T =  \"Block Number = " << Block_Number
                << "\" \\ \n"
                << "I = " << SolnBlk.Grid.INu - SolnBlk.Grid.INl + 1 << " \\ \n"
                << "J = " << SolnBlk.Grid.JNu - SolnBlk.Grid.JNl + 1 << " \\ \n"
                << "F = POINT \n";
    } else {
       Out_File << "ZONE T =  \"Block Number = " << Block_Number
                << "\" \\ \n"
                << "I = " << SolnBlk.Grid.INu - SolnBlk.Grid.INl + 1 << " \\ \n"
                << "J = " << SolnBlk.Grid.JNu - SolnBlk.Grid.JNl + 1 << " \\ \n"
                << "F = POINT \n";
    } /* endif */

    for ( j  = SolnBlk.Grid.JNl ; j <= SolnBlk.Grid.JNu ; ++j ) {
       for ( i = SolnBlk.Grid.INl ; i <= SolnBlk.Grid.INu ; ++i ) {
	   W_node = SolnBlk.Wn(i, j);
           Out_File << " "  << SolnBlk.Grid.Node[i][j].X << W_node;
           Out_File.setf(ios::scientific);
           Out_File << " " << sqrt(sqr(W_node.v.x)+sqr(W_node.v.y))
		    << " " << sqrt(sqr(W_node.v.x)+sqr(W_node.v.y))/W_node.sound()
		    << " " << W_node.T() <<  " " << W_node.DetP() << "\n";
           Out_File.unsetf(ios::scientific);
       } /* endfor */
    } /* endfor */
    Out_File << setprecision(6);
    
}

/********************************************************
 * Routine: Output_Cells_Tecplot                        *
 *                                                      *
 * Writes the cell centred solution values of the       *
 * specified quadrilateral solution block to the        *
 * specified output stream suitable for plotting with   *
 * TECPLOT.                                             *
 *                                                      *
 ********************************************************/
void Output_Cells_Tecplot(Gaussian2D_Quad_Block &SolnBlk,
			  Gaussian2D_Input_Parameters &IP,
                          const int Number_of_Time_Steps,
                          const double &Time,
                          const int Block_Number,
                          const int Output_Title,
	                  ostream &Out_File) {

    int i, j;

    /* Ensure boundary conditions are updated before
       outputting solution */

    BCs(SolnBlk,IP);

    Out_File << setprecision(14);
    if (Output_Title) {
       Out_File << "TITLE = \"" << CFFC_Name() << ": 2D Gaussian Solution, "
                << "Time Step/Iteration Level = " << Number_of_Time_Steps
                << ", Time = " << Time
                << "\"" << "\n"
   	        << "VARIABLES = \"x\" \\ \n"
                << "\"y\" \\ \n"
                << "\"rho\" \\ \n"
                << "\"u\" \\ \n"
                << "\"v\" \\ \n"
                << "\"pxx\" \\ \n"
                << "\"pxy\" \\ \n"
                << "\"pyy\" \\ \n"
                << "\"pzz\" \\ \n"
                << "\"erot\" \\ \n"
		<< "\"v_abs\" \\ \n"
                << "\"T\" \\ \n"
                << "\"DetP\" \\ \n"
                << "\"d_rho_u\" \\ \n"
                << "\"Mean Free Path\" \\ \n"
#ifdef _GAUSSIAN_HEAT_TRANSFER_
		<< "\"Qxxx\" \\ \n"
		<< "\"Qxxy\" \\ \n"
		<< "\"Qxyy\" \\ \n"
		<< "\"Qxzz\" \\ \n"
		<< "\"Qyyy\" \\ \n"
		<< "\"Qyzz\" \\ \n"
		<< "\"qx\" \\ \n"
		<< "\"qy\" \\ \n"
		<< "\"K\" \\ \n"
#endif
                << "\"drhodx\" \\ \n"
                << "\"dudx\" \\ \n"
                << "\"dvdx\" \\ \n"
                << "\"dpxxdx\" \\ \n"
                << "\"dpxydx\" \\ \n"
                << "\"dpyydx\" \\ \n"
                << "\"dpzzdx\" \\ \n"
                << "\"derotdx\" \\ \n"
                << "\"drhody\" \\ \n"
                << "\"dudy\" \\ \n"
                << "\"dvdy\" \\ \n"
                << "\"dpxxdy\" \\ \n"
                << "\"dpxydy\" \\ \n"
                << "\"dpyydy\" \\ \n"
                << "\"dpzzdy\" \\ \n"
                << "\"derotdy\" \\ \n"
                << "\"phi_d\" \\ \n"
                << "\"phi_vx\" \\ \n"
                << "\"phi_vy\" \\ \n"
                << "\"phi_pxx\" \\ \n"
                << "\"phi_pxy\" \\ \n"
                << "\"phi_pyy\" \\ \n"
                << "\"phi_pzz\" \\ \n"
                << "\"phi_erot\" \\ \n"
                << "\"abs_grad_rho\" \\ \n"
                << "\"abs_div_v\" \\ \n"
		<< "\"abs_curl_v_z\" \n"
		<< "\"cell_area\" \n"
		<< "\"i\" \n"
		<< "\"j\" \n"
                << "ZONE T =  \"Block Number = " << Block_Number
                << "\" \\ \n"
                << "I = " << SolnBlk.Grid.ICu - SolnBlk.Grid.ICl + 2*SolnBlk.Nghost + 1 << " \\ \n"
                << "J = " << SolnBlk.Grid.JCu - SolnBlk.Grid.JCl + 2*SolnBlk.Nghost + 1 << " \\ \n"
                << "F = POINT \n";
    } else {
       Out_File << "ZONE T =  \"Block Number = " << Block_Number
                << "\" \\ \n"
                << "I = " << SolnBlk.Grid.ICu - SolnBlk.Grid.ICl + 2*SolnBlk.Nghost + 1 << " \\ \n"
                << "J = " << SolnBlk.Grid.JCu - SolnBlk.Grid.JCl + 2*SolnBlk.Nghost + 1 << " \\ \n"
                << "F = POINT \n";
    } /* endif */

    for ( j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
       for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {

#ifdef _GAUSSIAN_HEAT_TRANSFER_
	   //note: this does not compute the heat terms in exactly the same way as in
	   //      the residual calculation.  Be sure you know what you're computing.
	   SolnBlk.W[i][j].ComputeHeatTerms(SolnBlk.dWdx[i][j],SolnBlk.dWdy[i][j],
					    Vector2D(0.0,0.0),0);
#endif
	 
           Out_File << " "  << SolnBlk.Grid.Cell[i][j].Xc
                    << SolnBlk.W[i][j];
           Out_File.setf(ios::scientific);
	   Out_File << " " << sqrt(sqr(SolnBlk.W[i][j].v.x)+sqr(SolnBlk.W[i][j].v.y))
		    << " " << SolnBlk.W[i][j].T() << " " << SolnBlk.W[i][j].DetP()
		    << " " << SolnBlk.dUdt[i][j][0][1]<< " " << SolnBlk.W[i][j].mfp() 
#ifdef _GAUSSIAN_HEAT_TRANSFER_
		    << " " << SolnBlk.W[i][j].q
		    << " " << (SolnBlk.W[i][j].q.xxx+SolnBlk.W[i][j].q.xyy+SolnBlk.W[i][j].q.xzz)/2.0
		    << " " << (SolnBlk.W[i][j].q.xxy+SolnBlk.W[i][j].q.yyy+SolnBlk.W[i][j].q.yzz)/2.0
		    << " " << SolnBlk.W[i][j].K()
#endif
		    << " " << SolnBlk.dWdx[i][j]
		    << " " << SolnBlk.dWdy[i][j]
		    << " " << SolnBlk.phi[i][j] 
		    << " " << sqrt(sqr(SolnBlk.dWdx[i][j].d) + sqr(SolnBlk.dWdy[i][j].d))
		    << " " << fabs(SolnBlk.dWdx[i][j].v.x + SolnBlk.dWdy[i][j].v.y)
		    << " " << fabs(SolnBlk.dWdx[i][j].v.y - SolnBlk.dWdy[i][j].v.x)
		    << " " << SolnBlk.Grid.Cell[i][j].A
		    << " " << i
		    << " " << j << "\n";
           Out_File.unsetf(ios::scientific);
       } /* endfor */
    } /* endfor */
    Out_File << setprecision(6);
    
}

/**********************************************************************
 * Routine: Output_Flat_Plate                                         *
 *                                                                    *
 * This routine outputs the non-dimensionalized computed flat plate   *
 * solution and the corresponding Blasius solution.                   *
 *                                                                    *
 **********************************************************************/
void Output_Flat_Plate(Gaussian2D_Quad_Block &SolnBlk,
		       const int Number_of_Time_Steps,
		       const double &Time,
		       const int Block_Number,
		       const int Output_Title,
		       ostream &Out_File,
		       const Gaussian2D_pState &Winf) {

  Gaussian2D_pState W, We;
  Vector2D X;
  double eta, f, fp, fpp, Rex, linf, Cf, Cfe;

  // Output node solution data.  
  Out_File << setprecision(14);
  if (Output_Title) {
    Out_File      << "TITLE = \"" << CFFC_Name() << ": 2D Gaussian Solution, "
                  << "Time Step/Iteration Level = " << Number_of_Time_Steps
                  << ", Time = " << Time
		  << "\"" << "\n"
		  << "VARIABLES = \"x\" \\ \n"
		  << "\"y\" \\ \n"
		  << "\"rho\" \\ \n"
		  << "\"u\" \\ \n"
		  << "\"v\" \\ \n"
		  << "\"pxx\" \\ \n"
		  << "\"pxy\" \\ \n"
		  << "\"pyy\" \\ \n"
		  << "\"pzz\" \\ \n"
		  << "\"erot\" \\ \n"
		  << "\"u_exact\" \\ \n"
		  << "\"v_exact\" \\ \n"
		  << "\"eta\" \\ \n"
		  << "\"f\" \\ \n"
		  << "\"fp\" \\ \n"
		  << "\"fpp\" \\ \n"
		  << "\"u/u0\" \\ \n"
		  << "\"u_exact/u0\" \\ \n"
                  << "ZONE T =  \"Block Number = " << Block_Number
		  << "\" \\ \n"
		  << "I = " << SolnBlk.Grid.ICu - SolnBlk.Grid.ICl + 1 << " \\ \n"
		  << "J = " << SolnBlk.Grid.JCu - SolnBlk.Grid.JCl + 1 << " \\ \n"
		  << "F = POINT \n";
  } else {
    Out_File      << "ZONE T =  \"Block Number = " << Block_Number
		  << "\" \\ \n"
		  << "I = " << SolnBlk.Grid.ICu - SolnBlk.Grid.ICl + 1 << " \\ \n"
		  << "J = " << SolnBlk.Grid.JCu - SolnBlk.Grid.JCl + 1 << " \\ \n"
		  << "F = POINT \n";
  }

  for (int j = SolnBlk.Grid.JCl; j <= SolnBlk.Grid.JCu; j++) {
    for (int i = SolnBlk.Grid.ICl; i <= SolnBlk.Grid.ICu; i++) {
      // Get cell position and solution data.
      X = SolnBlk.Grid.Cell[i][j].Xc;
      W = SolnBlk.W[i][j];
      // Get exact solution.
      We = FlatPlate(Winf,X,eta,f,fp,fpp);
      // Output data.
      Out_File      << " " << X
		    << W
		    << We.v;
      Out_File.setf(ios::scientific);
      Out_File	    << " " << eta
		    << " " << f
		    << " " << fp
		    << " " << fpp
		    << " " << W.v.x/Winf.v.x
		    << " " << We.v.x/Winf.v.x
		    << endl;
      Out_File.unsetf(ios::scientific);
    } /* endfor */
  } /* endfor */
  Out_File << setprecision(6);

}

/**********************************************************************
 * Routine: Output_Cylinder_Free_Molecular                            *
 *                                                                    *
 * This routine outputs the boudary pressure and shear profile for a  *
 * cylinder in free molecular flow.                                   *
 *                                                                    *
 **********************************************************************/
void Output_Cylinder_Free_Molecular(Gaussian2D_Quad_Block &SolnBlk,
				    const int Number_of_Time_Steps,
				    const double &Time,
				    const int Block_Number,
				    const int Output_Title,
				    ostream &Out_File,
				    const Gaussian2D_pState &Winf,
				    Vector2D *nodes,
				    const int number_of_nodes) {

  int i, j;
  double psi, cos_angle, sin_angle, Pn, tau, Pn_exact, tau_exact, cpsi;
  double T_ratio(1.0); //ratio of temperature between incoming and reflected particles
  Vector2D dX;
  Gaussian2D_pState W_node, W_exact;

  /* Ensure boundary conditions are updated before
     evaluating solution at the nodes. Linear Extrapolation
     is used to ensure ouput is not affected by "strange"
     ghostcell information */

  //left and right
  for ( j = SolnBlk.Grid.JCl; j <= SolnBlk.Grid.JCu; j++){
    //left
    if(SolnBlk.Grid.BCtypeW[j] != BC_NONE) {
      Linear_Reconstruction_LeastSquares(SolnBlk, 
					 SolnBlk.ICl, j, 
					 LIMITER_BARTH_JESPERSEN);
      dX = SolnBlk.Grid.Cell[SolnBlk.ICl-1][j].Xc -
	SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc;
      SolnBlk.W[SolnBlk.ICl-1][j] = SolnBlk.W[SolnBlk.ICl][j] + 
	(SolnBlk.phi[SolnBlk.ICl][j]^SolnBlk.dWdx[SolnBlk.ICl][j])*dX.x +
	(SolnBlk.phi[SolnBlk.ICl][j]^SolnBlk.dWdy[SolnBlk.ICl][j])*dX.y;
      SolnBlk.U[SolnBlk.ICl-1][j] = U(SolnBlk.W[SolnBlk.ICl-1][j]);
    }
    //right
    if(SolnBlk.Grid.BCtypeE[j] != BC_NONE) {
      Linear_Reconstruction_LeastSquares(SolnBlk, 
					 SolnBlk.ICu, j, 
					 LIMITER_BARTH_JESPERSEN);
      dX = SolnBlk.Grid.Cell[SolnBlk.ICu+1][j].Xc -
	SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc;
      SolnBlk.W[SolnBlk.ICu+1][j] = SolnBlk.W[SolnBlk.ICu][j] + 
	(SolnBlk.phi[SolnBlk.ICu][j]^SolnBlk.dWdx[SolnBlk.ICu][j])*dX.x +
	(SolnBlk.phi[SolnBlk.ICu][j]^SolnBlk.dWdy[SolnBlk.ICu][j])*dX.y;
      SolnBlk.U[SolnBlk.ICu+1][j] = U(SolnBlk.W[SolnBlk.ICu+1][j]);
    }
  }

  //bottom and top
  for ( i = SolnBlk.Grid.ICl; i <= SolnBlk.Grid.ICu; i++){
    //bottom
    if(SolnBlk.Grid.BCtypeS[i] != BC_NONE) {
      Linear_Reconstruction_LeastSquares(SolnBlk, 
					 i, SolnBlk.JCl, 
					 LIMITER_BARTH_JESPERSEN);
      dX = SolnBlk.Grid.Cell[i][SolnBlk.JCl-1].Xc -
	SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc;
      SolnBlk.W[i][SolnBlk.JCl-1] = SolnBlk.W[i][SolnBlk.JCl] + 
	(SolnBlk.phi[i][SolnBlk.JCl]^SolnBlk.dWdx[i][SolnBlk.JCl])*dX.x +
	(SolnBlk.phi[i][SolnBlk.JCl]^SolnBlk.dWdy[i][SolnBlk.JCl])*dX.y;
      SolnBlk.U[i][SolnBlk.JCl-1] = U(SolnBlk.W[i][SolnBlk.JCl-1]);
    }
    //top
    if(SolnBlk.Grid.BCtypeN[i] != BC_NONE) {
      Linear_Reconstruction_LeastSquares(SolnBlk, 
					 i, SolnBlk.JCu, 
					 LIMITER_BARTH_JESPERSEN);
      dX = SolnBlk.Grid.Cell[i][SolnBlk.JCu+1].Xc -
	SolnBlk.Grid.Cell[i][SolnBlk.JCu].Xc;
      SolnBlk.W[i][SolnBlk.JCu+1] = SolnBlk.W[i][SolnBlk.JCu] + 
	(SolnBlk.phi[i][SolnBlk.JCu]^SolnBlk.dWdx[i][SolnBlk.JCu])*dX.x +
	(SolnBlk.phi[i][SolnBlk.JCu]^SolnBlk.dWdy[i][SolnBlk.JCu])*dX.y;
      SolnBlk.U[i][SolnBlk.JCu+1] = U(SolnBlk.W[i][SolnBlk.JCu+1]);
    }
  }

  /* Output node solution data. */

  Out_File << setprecision(14);
  if (Output_Title) {
    Out_File << "TITLE = \"" << CFFC_Name() << ": 2D Gaussian Solution, "
	     << "Time Step/Iteration Level = " << Number_of_Time_Steps
	     << ", Time = " << Time
	     << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
	     << "\"y\" \\ \n"
	     << "\"rho\" \\ \n"
	     << "\"u\" \\ \n"
	     << "\"v\" \\ \n"
	     << "\"pxx\" \\ \n"
	     << "\"pxy\" \\ \n"
	     << "\"pyy\" \\ \n"
	     << "\"pzz\" \\ \n"
	     << "\"erot\" \\ \n"
	     << "\"rho_exact\" \\ \n"
	     << "\"u_exact\" \\ \n"
	     << "\"v_exact\" \\ \n"
	     << "\"pxx_exact\" \\ \n"
	     << "\"pxy_exact\" \\ \n"
	     << "\"pyy_exact\" \\ \n"
	     << "\"pzz_exact\" \\ \n"
	     << "\"erot_exact\" \\ \n"
	     << "ZONE T =  \"Block Number = " << Block_Number
	     << "\" \\ \n"
	     << "I = " << SolnBlk.Grid.INu - SolnBlk.Grid.INl + 1 << " \\ \n"
	     << "J = " << SolnBlk.Grid.JNu - SolnBlk.Grid.JNl + 1 << " \\ \n"
	     << "F = POINT \n";
  } else {
    Out_File << "ZONE T =  \"Block Number = " << Block_Number
	     << "\" \\ \n"
	     << "I = " << SolnBlk.Grid.INu - SolnBlk.Grid.INl + 1 << " \\ \n"
	     << "J = " << SolnBlk.Grid.JNu - SolnBlk.Grid.JNl + 1 << " \\ \n"
	     << "F = POINT \n";
  } /* endif */

  for ( j  = SolnBlk.Grid.JNl ; j <= SolnBlk.Grid.JNu ; ++j ) {
    for ( i = SolnBlk.Grid.INl ; i <= SolnBlk.Grid.INu ; ++i ) {
      W_node = SolnBlk.Wn(i, j);
      W_exact = Free_Molecular_Exact(SolnBlk.Grid.Node[i][j].X, 
				     Winf, 
				     nodes, 
				     number_of_nodes);

      Out_File << " "  << SolnBlk.Grid.Node[i][j].X << W_node << W_exact << endl;
    } /* endfor */
  } /* endfor */
  Out_File << setprecision(6);

}

/**********************************************************************
 * Routine: Append_nodes_to_send_buffer                               *
 *                                                                    *
 * This routine appends node location to a send buffer for the        *
 * free molecular solution procedure.                                 *
 *                                                                    *
 **********************************************************************/

void Append_nodes_to_send_buffer(Gaussian2D_Quad_Block &SolnBlk,
				 double *buffer1, double *buffer2, 
				 int &count) {

  for(int i = SolnBlk.Grid.ICl; i <= SolnBlk.Grid.ICu; i++) {

    if(SolnBlk.Grid.BCtypeS[i] == BC_ADIABATIC_WALL){
      buffer1[count] = SolnBlk.Grid.Node[i][SolnBlk.Grid.JNl].X.x;
      buffer2[count] = SolnBlk.Grid.Node[i][SolnBlk.Grid.JNl].X.y;
      count++;
    }
  }

  return;
}



/**********************************************************************
 * Routine: Output_Drag                                               *
 *                                                                    *
 * This routine outputs the Drag on a body (only south face is done!) *
 *                                                                    *
 **********************************************************************/
void Output_Drag(Gaussian2D_Quad_Block &SolnBlk,
		 double &drag, double &lift) {
  int i;
  double cos_angle, sin_angle;
  double normal_pressure, shear_pressure, cell_drag, cell_lift;
  Vector2D dX;
  Gaussian2D_pState W;

  //North and South walls

  for(i = SolnBlk.Grid.ICl; i <= SolnBlk.Grid.ICu; i++) {

    if(SolnBlk.Grid.BCtypeS[i] == BC_ADIABATIC_WALL ||
       SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
       SolnBlk.Grid.BCtypeS[i] == BC_TEMPERATURE_SLIP){

      Linear_Reconstruction_LeastSquares(SolnBlk,i,SolnBlk.Grid.JCl,
					 LIMITER_VENKATAKRISHNAN);

      dX = SolnBlk.Grid.xfaceS(i, SolnBlk.Grid.JCl)-SolnBlk.Grid.Cell[i][SolnBlk.Grid.JCl].Xc;
      W = SolnBlk.W[i][SolnBlk.Grid.JCl] +
        (SolnBlk.phi[i][SolnBlk.Grid.JCl]^SolnBlk.dWdx[i][SolnBlk.Grid.JCl])*dX.x +
	(SolnBlk.phi[i][SolnBlk.Grid.JCl]^SolnBlk.dWdy[i][SolnBlk.Grid.JCl])*dX.y;

      ////////////////////////////////////////////////////
      /*
      W = SolnBlk.W[i][SolnBlk.Grid.JCl];

      cout << endl << "Cell center = " <<SolnBlk.Grid.Cell[i][SolnBlk.Grid.JCl].Xc << endl;
      cout << "Face length = " << SolnBlk.Grid.lfaceS(i,SolnBlk.Grid.JCl) << endl;
      cout << "W = " << W << endl;
      */
      ////////////////////////////////////////////////////

      // Determine the direction cosine's for the frame
      // rotation.

      cos_angle = SolnBlk.Grid.nfaceS(i,SolnBlk.Grid.JCl).x; 
      sin_angle = SolnBlk.Grid.nfaceS(i,SolnBlk.Grid.JCl).y;

      // Apply the frame rotation and calculate the primitive
      // solution state variables in the local rotated frame
      // defined by the unit normal vector. 

      normal_pressure  = W.p.xx*cos_angle*cos_angle+W.p.yy*sin_angle*sin_angle
 	                 +2.0*W.p.xy*cos_angle*sin_angle;
      shear_pressure  = W.p.xy*(cos_angle*cos_angle-sin_angle*sin_angle)
	                -(W.p.xx-W.p.yy)*cos_angle*sin_angle;

      //shear_pressure = 0.0;

      cell_drag = (normal_pressure*cos_angle-shear_pressure*sin_angle)*
	          SolnBlk.Grid.lfaceS(i,SolnBlk.Grid.JCl);

      cell_lift = (normal_pressure*sin_angle+shear_pressure*cos_angle)*
	          SolnBlk.Grid.lfaceS(i,SolnBlk.Grid.JCl);
      ///////////////////////////////////////////////
      /*
      cout << "Pn = " << normal_pressure << endl;
      cout << "Pt = " << shear_pressure << endl;
      cout << "cell_drag = " << cell_drag << endl;
      cout << "cell_lift = " << cell_lift << endl;
      */
      //////////////////////////////////////////////

      drag += cell_drag;
      lift += cell_lift;
    }

    if(SolnBlk.Grid.BCtypeN[i] == BC_ADIABATIC_WALL ||
       SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
       SolnBlk.Grid.BCtypeN[i] == BC_TEMPERATURE_SLIP){
      //Not Done
    }
  }

  //East and West walls

  for(i = SolnBlk.Grid.JCl; i <= SolnBlk.Grid.JCu; i++) {

    if(SolnBlk.Grid.BCtypeE[i] == BC_ADIABATIC_WALL ||
       SolnBlk.Grid.BCtypeE[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
       SolnBlk.Grid.BCtypeE[i] == BC_TEMPERATURE_SLIP){
      //Not Done
    }

    if(SolnBlk.Grid.BCtypeW[i] == BC_ADIABATIC_WALL ||
       SolnBlk.Grid.BCtypeW[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
       SolnBlk.Grid.BCtypeW[i] == BC_TEMPERATURE_SLIP){
      //Not Done
    }
  }

}

/**********************************************************************
 * Routine: Output_Gradients_Tecplot                                  *
 *                                                                    *
 * Writes the cell centred primitive variable state gradients and     *
 * limiters of the specified quadrilateral solution block to the      *
 * specified output stream suitable for plotting with TECPLOT.        *
 *                                                                    *
 **********************************************************************/
void Output_Gradients_Tecplot(Gaussian2D_Quad_Block &SolnBlk,
			      const int Number_of_Time_Steps,
			      const double &Time,
			      const int Block_Number,
			      const int Output_Title,
			      ostream &Out_File) {

  cout << "This has not been implemented." << endl;

  //this is here only for compalibility with embedded boundaries
  return;
}

/********************************************************
 * Routine: Output_Shock_Structure                      *
 *                                                      *
 * Writes 1-D shock structure output.  This function    *
 * will write the solution for the first cell above the *
 * x axis.  It is therefore necessary that all blocks   *
 * be at the same refinement level.                     *
 *                                                      *
 ********************************************************/

void Output_Shock_Structure(Gaussian2D_Quad_Block &SolnBlk,
			    Gaussian2D_Input_Parameters &IP,
			    const int Number_of_Time_Steps,
			    const double &Time,
			    const int Block_Number,
			    int &Output_Title,
			    ostream &Out_File,
			    const double &y,
			    const double &y_tol) {
  int i, j, output_flag(0);

    /* Ensure boundary conditions are updated before
       outputting solution */

    BCs(SolnBlk,IP);

    Out_File << setprecision(14);

    // determine whether to output
    for ( j  = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
      if(fabs(SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc.y-y)<y_tol && 
	 SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc.y > 0.0) {
	output_flag = 1;
	break;
      }
    }

    if(output_flag) {
      if (Output_Title) {
	Out_File << "TITLE = \"" << CFFC_Name() << ": 2D Gaussian Solution, "
		 << "Time Step/Iteration Level = " << Number_of_Time_Steps
		 << ", Time = " << Time
		 << "\"" << "\n"
		 << "VARIABLES = \"x\" \\ \n"
		 << "\"y\" \\ \n"
		 << "\"Mean Free Path\" \\ \n"
		 << "\"x/mfp_ref\" \\ \n"
		 << "\"rho\" \\ \n"
		 << "\"u\" \\ \n"
		 << "\"v\" \\ \n"
		 << "\"pxx\" \\ \n"
		 << "\"pxy\" \\ \n"
		 << "\"pyy\" \\ \n"
		 << "\"pzz\" \\ \n"
		 << "\"erot\" \\ \n"
		 << "\"T\" \\ \n"
		 << "\"DetP\" \\ \n"
		 << "ZONE T =  \"Block Number = " << Block_Number
		 << "\" \\ \n"
		 << "I = " << SolnBlk.Grid.ICu - SolnBlk.Grid.ICl + 1 << " \\ \n"
		 << "F = POINT \n";
	Output_Title = 0;
      } else {
	Out_File << "ZONE T =  \"Block Number = " << Block_Number
		 << "\" \\ \n"
		 << "I = " << SolnBlk.Grid.ICu - SolnBlk.Grid.ICl +1 << " \\ \n"
		 << "F = POINT \n";
      } /* endif */

      for ( j  = SolnBlk.JCl ; j <= SolnBlk.JCu ; ++j ) {
	for ( i = SolnBlk.ICl ; i <= SolnBlk.ICu ; ++i ) {
	  if(fabs(SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc.y-y)<y_tol && 
	     SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc.y > 0.0) {
	    Out_File << " "  << SolnBlk.Grid.Cell[i][j].Xc;
	    Out_File.setf(ios::scientific);
	    Out_File << " "  << SolnBlk.W[i][j].mfp() 
		     << " "  << SolnBlk.Grid.Cell[i][j].Xc.x/SolnBlk.WoN[0].mfp(); //I don't use IP.Wo because it is not reset after a restart.
	    Out_File.unsetf(ios::scientific);
	    Out_File << SolnBlk.W[i][j];
	    Out_File.setf(ios::scientific);
	    Out_File << " " << SolnBlk.W[i][j].T() << " " << SolnBlk.W[i][j].DetP() << "\n";
	    Out_File.unsetf(ios::scientific);
	  }
	} /* endfor */
      } /* endfor */
    }
    Out_File << setprecision(6);
}

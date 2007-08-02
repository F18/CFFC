
/********************************************************
 ********************* CHEM2D Tools *********************
 ********************************************************/
#ifndef _CHEM2D_QUAD_INCLUDED
#include "Chem2DQuad.h"
#endif // _CHEM2D_QUAD_INCLUDED

/********************************************************
 * Routine: Open_Time_Accurate_File                     *
 *                                                      *
 * This routine opens a file for time accurate data     *
 *                                                      *
 *                                                      *
 ********************************************************/
int Open_Time_Accurate_File(ofstream &Time_Accurate_File,
			    char *File_Name,
			    const int Append_to_File,
			    const Chem2D_pState &Soln) {

    int i;
    char prefix[256], extension[256], 
      time_accurate_file_name[256];
    char *time_accurate_name_ptr;
    
    /* Determine the name of the file. */

    i = 0;
    while (1) {  //wow this is weird way to do this
       if (File_Name[i] == ' ' ||
           File_Name[i] == '.') break;
       prefix[i] = File_Name[i];
       i = i + 1;
       if (i > strlen(File_Name) ) break;
    } /* endwhile */
    prefix[i] = '\0';
    strcat(prefix, "_time_accurate");
    strcpy(extension, ".dat");
    strcpy(time_accurate_file_name, prefix);
    strcat(time_accurate_file_name, extension);

    time_accurate_name_ptr = time_accurate_file_name;

    /* Open the file. */

    if (Append_to_File) {
       Time_Accurate_File.open(time_accurate_name_ptr, ios::out|ios::app);
    } else {
       Time_Accurate_File.open(time_accurate_name_ptr, ios::out);
    } /* endif */
    if (Time_Accurate_File.bad()) return (1);

    if(!Append_to_File){
      /* Write the appropriate Tecplot header information. */
      Time_Accurate_File << "TITLE = \" Unsteady Mass Fraction plots "
			 << "\"" << "\n"  
			 << "VARIABLES = \"time\" \\ \n"
			 << "\"rho\" \\ \n"
			 << "\"u\" \\ \n"
			 << "\"v\" \\ \n"
			 << "\"p\" \\ \n"
			 << "\"k\" \\ \n"
			 << "\"omega\" \\ \n";
      //n species mass fractions names
      for(int i =0 ;i<Soln.ns ;i++){
	Time_Accurate_File <<"\"c"<<Soln.specdata[i].Speciesname()<<"\" \\ \n";
	//	cout<<"\"c"<<Soln.specdata[i].Speciesname()<<"\" \\ \n";
      }   
      //Viscous Terms 
      Time_Accurate_File << "\"qflux_x\" \\ \n"  
			 << "\"qflux_y\" \\ \n"   
			 << "\"Tau_xx\" \\ \n"  //rr -axisymmetric
			 << "\"Tau_xy\" \\ \n"  //rz
			 << "\"Tau_yy\" \\ \n"  //zz
			 << "\"Tau_zz\" \\ \n" //thetatheta 
			 << "\"theta_x\" \\ \n"  
			 << "\"theta_y\" \\ \n"   
			 << "\"lambda_xx\" \\ \n"   //rr -axisymmetric
			 << "\"lambda_xy\" \\ \n"   //rz
			 << "\"lambda_yy\" \\ \n"   //zz
			 << "\"lambda_zz\" \\ \n"	
			 << "\"T\" \\ \n"
			 << "\"e\" \\ \n"
			 << "\"es\" \\ \n"
			 << "\"e_ref\" \\ \n"; 
//       for(int i=0;i<Soln.React.num_reactions; i++){
// 	Time_Accurate_File <<"\"kf"<<i<<"\" \\ \n";
// 	Time_Accurate_File <<"\"kb"<<i<<"\" \\ \n";
//       }
    }
    return(0);

}

/********************************************************
 * Routine: Close_Time_Accurate                         *
 *                                                      *
 * This routine closes the time accurate data file      *
 * for a CFD  calculation.                              *
 *                                                      *
 ********************************************************/
int Close_Time_Accurate_File(ofstream &Time_Accurate_File) {
    Time_Accurate_File.close();
    return(0);
}

/********************************************************
 * Routine: Output_Progress_to_File                     *
 *                                                      *
 * This routine writes out progress information for     *
 * a CFD calculation to a progress file, including      *
 * iteration level, time, CPU time, and residual norms. *
 *                                                      *
 ********************************************************/
void Output_to_Time_Accurate_File(ostream &Time_Accurate_File,
			       const double &Time,
			       const Chem2D_pState &Soln){
  double Temp=Soln.T();

  Time_Accurate_File << setprecision(6);
  Time_Accurate_File << Time << " " <<Soln 
		     <<" "<<Temp
		     <<" "<<Soln.e()
		     <<" "<<Soln.es()
                     <<" "<<Soln.eref();
//   for(int i=0;i<Soln.React.num_reactions; i++){
//     double kf = Soln.React.reactions[i].kf(Temp);
//     Time_Accurate_File <<" "<<kf
//                        <<" "<<kf/Soln.React.reactions[i].keq(Soln,Temp);
//   }
  Time_Accurate_File <<endl;

//Time_Accurate_File.setf(ios::scientific); 
//Time_Accurate_File.unsetf(ios::scientific);
    Time_Accurate_File.flush();
}

// FROM DUSTY2D


/**********************************************************************
 * Routine: Output_Ringleb_Solution                                   *
 *                                                                    *
 * Writes the exact and computed Ringleb's flow solution values at    *
 * the nodes of the specified quadrilateral solution block to the     *
 * specified output stream suitable for plotting with TECPLOT.  The   *
 * error norms are also computed.                                     *
 *                                                                    *
 **********************************************************************/
void Output_Ringleb_Solution(Chem2D_Quad_Block &SolnBlk,
			     const int Block_Number,
			     const int Output_Title,
			     ostream &Out_File) {

  Chem2D_pState We, Ww;

  // Output node solution data.  
  Out_File << setprecision(14);
  if (Output_Title)
    Out_File << "TITLE = \"" << CFDkit_Name() << ": 2D Dusty Ringleb Flow Solution "
	     << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
	     << "\"y\" \\ \n"
	     << "\"rho\" \\ \n"
	     << "\"a\" \\ \n"
	     << "\"rho_exact\" \\ \n"
	     << "\"a_exact\" \\ \n"
	     << "\"rho_error\" \\ \n"
	     << "\"a_error\" \\ \n";
//    Out_File << "ZONE T =  \"Block Number = " << Block_Number
//  	   << "\" \\ \n"
//  	   << "I = " << SolnBlk.Grid.ICu - SolnBlk.Grid.ICl + 1 << " \\ \n"
//  	   << "J = " << SolnBlk.Grid.JCu - SolnBlk.Grid.JCl + 1 << " \\ \n"
//  	   << "F = POINT \\ \n";

//    Out_File.setf(ios::scientific);
//    for (int j = SolnBlk.Grid.JCl; j <= SolnBlk.Grid.JCu; j++) {
//      for (int i = SolnBlk.Grid.ICl; i <= SolnBlk.Grid.ICu; i++) {
//        if (SolnBlk.Grid.cellstatus[i][j] == CELL_STATUS_ACTIVE) {
//  	W = RinglebFlow(W,SolnBlk.Grid.Cell[i][j].Xc);
//        } else {
//  	W = SolnBlk.W[i][j];
//        }
//        Out_File << SolnBlk.Grid.Cell[i][j].Xc << " " 
//  	       << SolnBlk.W[i][j].rho << " " 
//  	       << SolnBlk.W[i][j].a() << " " 
//  	       << W.rho << " " 
//  	       << W.a()   << " " 
//  	       << fabs(SolnBlk.W[i][j].rho - W.rho) << " " 
//  	       << fabs(SolnBlk.W[i][j].a() - W.a()) << endl;
//      }
//    }
  Out_File << "ZONE T =  \"Block Number = " << Block_Number
	   << "\" \\ \n"
	   << "I = " << SolnBlk.Grid.INu - SolnBlk.Grid.INl + 1 << " \\ \n"
	   << "J = " << SolnBlk.Grid.JNu - SolnBlk.Grid.JNl + 1 << " \\ \n"
	   << "F = POINT \\ \n";

  Out_File.setf(ios::scientific);
  for (int j = SolnBlk.Grid.JNl; j <= SolnBlk.Grid.JNu; j++) {
    for (int i = SolnBlk.Grid.INl; i <= SolnBlk.Grid.INu; i++) {
      Ww = SolnBlk.Wn(i,j);
      We = RinglebFlow(We,SolnBlk.Grid.Node[i][j].X);
      Out_File << SolnBlk.Grid.Node[i][j].X << " " 
	       << Ww.rho << " " 
	       << Ww.a() << " " 
	       << We.rho << " " 
	       << We.a()   << " " 
	       << fabs(Ww.rho - We.rho) << " " 
	       << fabs(Ww.a() - We.a()) << endl;
    }
  }

}

/**********************************************************************
 * Routine: Output_Ringleb_Error                                      *
 *                                                                    *
 *                                                                    *
 **********************************************************************/
void Output_Ringleb_Error(Chem2D_Quad_Block &SolnBlk,
			  double &l1_norm,
			  double &l2_norm,
			  double &max_norm,
			  int &numberofactivecells) {

  Chem2D_pState W;

  for (int j = SolnBlk.Grid.JCl; j <= SolnBlk.Grid.JCu; j++) {
    for (int i = SolnBlk.Grid.ICl; i <= SolnBlk.Grid.ICu; i++) {
	W = RinglebFlow(W,SolnBlk.Grid.Cell[i][j].Xc);
	l1_norm += fabs(W.rho - SolnBlk.W[i][j].rho);
	l2_norm += sqr(W.rho - SolnBlk.W[i][j].rho);
	max_norm = max(max_norm,fabs(W.rho - SolnBlk.W[i][j].rho));
    }
  }
  
  l2_norm = sqrt(l2_norm);
//    if (numberofactivecells > 0) {
//      l1_norm = l1_norm/numberofactivecells;
//      l2_norm = sqrt(l2_norm/numberofactivecells);
//    }

//    Vector2D X, X1, X2, X3, X4;
//    double *w;
//    Chem2D_pState **f;
//    double *epsilon;
//    double *eta;

//    double epsilon1 = -ONE, epsilon2 = ONE, epsilon3 = -ONE, epsilon4 = ONE;
//    double eta1 = -ONE, eta2 = -ONE, eta3 =  ONE, eta4 = ONE;
//    double N1, N2, N3, N4;

//    l1_norm = ZERO;
//    l2_norm = ZERO;
//    max_norm = ZERO;
 
//    // Allocate memory.
//    w = new double[3];
//    f = new Chem2D_pState*[3];
//    epsilon = new double[3];
//    eta = new double[3];
//    for (int i = 0; i < 3; i++) f[i] = new Chem2D_pState[3];

//    // Assign values... the dumbest way possible.
//    w[0] = 5.0/9.0;
//    w[1] = 8.0/9.0;
//    w[2] = 5.0/9.0;
//    epsilon[0] = -sqrt(15.0)/5.0;
//    epsilon[1] = 0.0;
//    epsilon[2] = sqrt(15.0)/5.0;
//    eta[0] = -sqrt(15.0)/5.0;
//    eta[1] = 0.0;
//    eta[2] = sqrt(15.0)/5.0;

//    numberofactivecells = 0;

//    // Get the value of 'f' at each point then perform the integration and finally compute the norm.
//    for (int j = SolnBlk.Grid.JCl; j <= SolnBlk.Grid.JCu; j++) {
//      for (int i = SolnBlk.Grid.ICl; i <= SolnBlk.Grid.ICu; i++) {

//        // For all active cells (required for embedded boundary testing).
//        if (SolnBlk.Grid.cellstatus[i][j] == CELL_STATUS_ACTIVE) {

//  	// Get the value of 'f' at each point.
//  	for (int jj = 0; jj < 3; jj++) {
//  	  for (int ii = 0; ii < 3; ii++) {

//  	    // Save node values.
//  	    X1 = SolnBlk.Grid.Node[i  ][j  ].X;
//  	    X2 = SolnBlk.Grid.Node[i+1][j  ].X;
//  	    X3 = SolnBlk.Grid.Node[i+1][j+1].X;
//  	    X4 = SolnBlk.Grid.Node[i  ][j+1].X;

//  	    // Set basis functions.
//  	    N1 = 0.25*(1 + epsilon[ii]*epsilon1)*(1 + eta[jj]*eta1);
//  	    N2 = 0.25*(1 + epsilon[ii]*epsilon2)*(1 + eta[jj]*eta2);
//  	    N3 = 0.25*(1 + epsilon[ii]*epsilon3)*(1 + eta[jj]*eta3);
//  	    N4 = 0.25*(1 + epsilon[ii]*epsilon4)*(1 + eta[jj]*eta4);

//  	    // Get point X.
//  	    X = X1*N1 + X2*N2 + X3*N3 + X4*N4;

//  	    // Determine the value of 'f'.
//  	    f[ii][jj] = RinglebFlow(f[ii][jj],X);
	  
//  	  }
//  	}

//  	// Perform the integration to find the average value of the exact value of rho.
//  	W = Chem2D_W_VACUUM;
//  	for (int jj = 0; jj < 3; jj++)
//  	  for (int ii = 0; ii < 3; ii++)
//  	    W += w[ii]*w[jj]*f[ii][jj]/FOUR;

//  	// Compute the norms.
//  	l1_norm += fabs(W.rho - SolnBlk.W[i][j].rho);
//  	l2_norm += sqr(W.rho - SolnBlk.W[i][j].rho);
//  	max_norm = max(max_norm,fabs(W.rho - SolnBlk.W[i][j].rho));

//  	// Increment the number of active cells.
//  	numberofactivecells++;

//        }

//      }
//    }

//    //l1_norm = l1_norm/((SolnBlk.NCi-2*SolnBlk.Nghost)*(SolnBlk.NCj-2*SolnBlk.Nghost));
//    //l2_norm = sqrt(l2_norm/((SolnBlk.NCi-2*SolnBlk.Nghost)*(SolnBlk.NCj-2*SolnBlk.Nghost)));
//    l1_norm = l1_norm/numberofactivecells;
//    l2_norm = sqrt(l2_norm/numberofactivecells);

//    // Delete memory.
//    for (int i = 0; i < 3; i++) {
//      delete []f[i]; f[i] = NULL;
//    }
//    delete []w; w = NULL;
//    delete []f; f = NULL;
//    delete []epsilon; epsilon = NULL;
//    delete []eta; eta = NULL;

}

/**********************************************************************
 * Routine: Output_Viscous_Channel                                    *
 *                                                                    *
 * Writes the exact and computed Viscous Channel flow solution values *
 * at the nodes of the specified quadrilateral solution block to the  *
 * specified output stream suitable for plotting with TECPLOT.  The   *
 * error norms are also computed.                                     *
 *                                                                    *
 **********************************************************************/
void Output_Viscous_Channel(Chem2D_Quad_Block &SolnBlk,
			    const int Block_Number,
			    const int Output_Title,
			    ostream &Out_File,
			    double &l1_norm,
			    double &l2_norm,
			    double &max_norm,
			    double &Vwall,
			    const double dp) {

  Chem2D_pState We, W;

  // Output node solution data.  
  Out_File << setprecision(14);
  if (Output_Title)
    Out_File << "TITLE = \"" << CFDkit_Name() << ": 2D Chem2D Channel Flow Solution "
	     << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
	     << "\"y\" \\ \n"
	     << "\"rho\" \\ \n"
	     << "\"u\" \\ \n"
	     << "\"v\" \\ \n"
	     << "\"u_exact\" \\ \n"
	     << "\"v_exact\" \\ \n"
	     << "\"p\" \\ \n"
	     << "\"T\" \\ \n"
	     << "\"u_error\" \\ \n"
	     << "\"v_error\" \\ \n";
  Out_File << "ZONE T =  \"Block Number = " << Block_Number
	   << "\" \\ \n"
	   << "I = " << SolnBlk.Grid.NCi << " \\ \n"
	   << "J = " << SolnBlk.Grid.NCj << " \\ \n"
	   << "F = POINT \\ \n";

  Out_File.setf(ios::scientific);
  for (int j = 0; j < SolnBlk.Grid.NCj; j++) {
    for (int i = 0; i < SolnBlk.Grid.NCi; i++) {
      W = SolnBlk.W[i][j];
      if (i < SolnBlk.Nghost || i > SolnBlk.ICu ||
	  j < SolnBlk.Nghost || j > SolnBlk.JCu) {
	We = SolnBlk.W[i][j];
      } else {
	We = ViscousChannelFlow(W,SolnBlk.Grid.Cell[i][j].Xc,Vwall,dp);
	l1_norm += fabs(W.v.x - We.v.x);
	l2_norm += sqr(W.v.x - We.v.x);
	max_norm = max(max_norm,fabs(W.v.x - We.v.x));
      }
      Out_File << SolnBlk.Grid.Cell[i][j].Xc << " " 
	       << W.rho << " " 
	       << W.v.x << " " 
	       << W.v.y << " " 
	       << We.v.x << " " 
	       << We.v.y << " " 
	       << W.p   << " " 
	       << W.T() << " " 
	       << fabs(W.v.x - We.v.x) << " "
	       << fabs(W.v.y - We.v.y) << endl;
    }
  }
  l2_norm = sqrt(l2_norm);

}


/**********************************************************************
 * Routine: Output_Flat_Plate                                         *
 *                                                                    *
 * This routine outputs the non-dimensionalized computed flat plate   *
 * solution and the corresponding Blasius solution.                   *
 *                                                                    *
 **********************************************************************/
void Output_Flat_Plate(Chem2D_Quad_Block &SolnBlk,
		       const int Block_Number,
		       const int Output_Title_Soln,
		       ostream &Out_File_Soln,
		       const int Output_Title_Skin,
		       ostream &Out_File_Skin,
		       const Chem2D_pState &Winf,
		       double &l1_norm,
		       double &l2_norm,
		       double &max_norm){

  Chem2D_pState W, We;
  Vector2D X;
  double eta, f, fp, fpp, Rex, linf, Cf, Cfe, xpt;

  // Output node solution data.  
  Out_File_Soln << setprecision(14);
  if (Output_Title_Soln) {
    Out_File_Soln << "TITLE = \"" << CFDkit_Name() << ": 2D Chem2D Solution, "
		  << "\"" << "\n"
		  << "VARIABLES = \"x\" \\ \n"
		  << "\"y\" \\ \n"
		  << "\"rho\" \\ \n"
		  << "\"vx\" \\ \n"
		  << "\"vy\" \\ \n"
		  << "\"p\" \\ \n"
		  << "\"rho_e\" \\ \n"
		  << "\"vx_e\" \\ \n"
		  << "\"vy_e\" \\ \n"
		  << "\"p_e\" \\ \n"
		  << "\"eta\" \\ \n"
		  << "\"f\" \\ \n"
		  << "\"fp\" \\ \n"
		  << "\"fpp\" \\ \n"
		  << "\"vx/Vxinf\" \\ \n"
		  << "\"vx_e/Vxinf\" \\ \n"
		  << "\"phiv\" \\ \n"
		  << "\"phiv_e\" \\ \n";
  }


  for (int npts = 0; npts < 6; npts++) {
    Out_File_Soln << "ZONE T =  \"Block Number = " << Block_Number
		  << ", xpt = " << 2000.0*double(npts)
		  << "\" \\ \n"
		  << "I = " << 1 << " \\ \n"
		  << "J = " << SolnBlk.Grid.JNu - SolnBlk.Grid.JNl + 1 << " \\ \n"
		  << "F = POINT \n";
    for (int j = SolnBlk.Grid.JNl; j <= SolnBlk.Grid.JNu; j++) {
      for (int i = SolnBlk.Grid.INl; i <= SolnBlk.Grid.INu-1; i++) {
	// Get cell position and solution data.
	X = SolnBlk.Grid.Node[i][j].X;

	W = SolnBlk.Wn(i,j);

	if (fabs(X.x - 0.20*double(npts)) < fabs(X.x - SolnBlk.Grid.Node[i+1][j].X.x) &&
	    X.x <= 0.20*double(npts) ) {
	  We = FlatPlate(Winf,X,eta,f,fp,fpp);
	  linf = 10000.0*(Winf.mu()/Winf.rho)/Winf.v.x;
	  Rex = (Winf.v.x/(Winf.mu()/Winf.rho))*((X.x+NANO)/ONE);
	  // Output data.
	  Out_File_Soln.setf(ios::scientific);
	  Out_File_Soln << " " << X
			<< " " << W.rho 
			<< W.v
			<< " " << W.p
			<< " " << We.rho
			<< We.v
			<< " " << We.p
			<< " " << eta
			<< " " << f
			<< " " << fp
			<< " " << fpp
			<< " " << W.v.x/Winf.v.x
			<< " " << We.v.x/Winf.v.x
			<< " " << (W.v.y/Winf.v.x)*sqrt(Winf.v.x/((Winf.mu()/Winf.rho)*(X.x+NANO)))
			<< " " << HALF*(eta*fp-f)
			<< endl;
	}
      }
    }
  }

  // Calculate the norms of the skin friction coefficient.
  for (int j = SolnBlk.Grid.JCl; j <= SolnBlk.Grid.JCu; j++) {
    for (int i = SolnBlk.Grid.ICl; i <= SolnBlk.Grid.ICu; i++) {
      // Get cell position and solution data.
      X = SolnBlk.Grid.Cell[i][j].Xc;
      W = SolnBlk.W[i][j];
      // Determine the skin friction coefficient.
      linf = 10000.0*(Winf.mu()/Winf.rho)/Winf.v.x;
      Rex = (Winf.v.x/(Winf.mu()/Winf.rho))*((X.x+NANO)/ONE);
      if (X.x >= ZERO && j == 2 && SolnBlk.Grid.BCtypeS[i] == BC_NO_SLIP) {//BC_WALL_VISCOUS_HEATFLUX
	// Get exact skin friction coefficient.
	Cfe = TWO*0.32206/sqrt(Rex);
	// Get computed skin friction coefficient.
	Cf  = TWO*WallShearStress(W,X,
				  SolnBlk.Grid.nodeSW(i,j).X,
				  SolnBlk.Grid.nodeSE(i,j).X,
				  Vector2D(ONE,ZERO))/(Winf.rho*Winf.v.x*Winf.v.x);
	// Calculate error norms.
	l1_norm += fabs(Cf - Cfe);
	l2_norm += sqr(Cf - Cfe);
	max_norm = max(max_norm,fabs(Cf - Cfe));
      } else {
	Cf  = ZERO;
	Cfe = ZERO;
      }
    }
  }
  l2_norm = sqrt(l2_norm);


  if (Output_Title_Skin) {
      Out_File_Skin << "TITLE = \"" << CFDkit_Name() << ": 2D Chem2D Solution, "
		    << "\"" << "\n"
		    << "VARIABLES = \"x\" \\ \n"
		    << "\"Rex\" \\ \n"
		    << "\"Cf\" \\ \n"
		    << "\"Cf_e\" \\ \n";
      Out_File_Skin << "ZONE T =  \"Block Number = " << Block_Number
		    << "\" \\ \n"
		    << "I = " <<  16 << " \\ \n"
		    << "J = " << 1 << " \\ \n"
		    << "F = POINT \n";
    } else {
      Out_File_Skin << "ZONE T =  \"Block Number = " << Block_Number
		    << "\" \\ \n"
		    << "I = " <<  16<< " \\ \n"
		    << "J = " << 1 << " \\ \n"
		    << "F = POINT \n";
    }
    for (int j = SolnBlk.Grid.JCl; j <= SolnBlk.Grid.JCu; j++) {
      for (int i = SolnBlk.Grid.ICl; i <= SolnBlk.Grid.ICu; i++) {
	  // Get cell position and solution data.
	  X = SolnBlk.Grid.Cell[i][j].Xc;
	  W = SolnBlk.W[i][j];
	  // Determine the skin friction coefficient.
	  linf = 10000.0*(Winf.mu()/Winf.rho)/Winf.v.x;
	  Rex = (Winf.v.x/(Winf.mu()/Winf.rho))*((X.x+NANO)/ONE);
	  if (X.x >= ZERO && j == 2 && SolnBlk.Grid.BCtypeS[i] == BC_NO_SLIP){	     
	    Cfe = TWO*0.32206/sqrt(Rex);
	    Cf  = TWO*WallShearStress(W,X,
				      SolnBlk.Grid.nodeSW(i,j).X,
				      SolnBlk.Grid.nodeSE(i,j).X,
				      Vector2D(ONE,ZERO))/(Winf.rho*Winf.v.x*Winf.v.x);
	    // Output data.
	    Out_File_Skin.setf(ios::scientific);
	    Out_File_Skin << " " << X.x
			  << " " << Rex*linf
			  << " " << Cf
			  << " " << Cfe
			  << endl;
	  }		      
      }
    }
  
}
 


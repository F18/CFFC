/////////////////////////////////////////////////////////////////////
///
/// \file Flame2DTools.cc
/// 
/// \author Marc R.J. Charest
/// 
/// \brief This header file contains the function definition for 
///        the Flame2DTools.  Basically functions that do not fit
///        anywhere else.
///
/////////////////////////////////////////////////////////////////////
#include "Flame2DTools.h"


/////////////////////////////////////////////////////////////////////
/// Time Accurate Ouput
/////////////////////////////////////////////////////////////////////

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
			    const Flame2D_pState &Soln) {

  int i;
  char prefix[256], extension[256], 
    time_accurate_file_name[256];
  char *time_accurate_name_ptr;
    
  // Determine the name of the file.
  i = 0;
  while (1) {  //wow this is weird way to do this
    if (File_Name[i] == ' ' ||
	File_Name[i] == '.') break;
    prefix[i] = File_Name[i];
    i = i + 1;
    if (i > strlen(File_Name) ) break;
  } // endwhile
  prefix[i] = '\0';
  strcat(prefix, "_time_accurate");
  strcpy(extension, ".dat");
  strcpy(time_accurate_file_name, prefix);
  strcat(time_accurate_file_name, extension);

  time_accurate_name_ptr = time_accurate_file_name;

  // Open the file.
  if (Append_to_File) {
    Time_Accurate_File.open(time_accurate_name_ptr, ios::out|ios::app);
  } else {
    Time_Accurate_File.open(time_accurate_name_ptr, ios::out);
  } /* endif */
  if (Time_Accurate_File.fail()) return (1);

  //
  // Write the appropriate Tecplot header information.
  //
  if(!Append_to_File){
    Time_Accurate_File << "TITLE = \" Unsteady Mass Fraction plots "
		       << "\"" << "\n"  
		       << "VARIABLES = \"time\" \\ \n"
		       << "\"rho\" \\ \n"
		       << "\"u\" \\ \n"
		       << "\"v\" \\ \n"
		       << "\"p\" \\ \n";
    //n species mass fractions names
    for(int i =0; i<Flame2D_pState::NumSpecies(); i++){
      Time_Accurate_File <<"\"c"<<Flame2D_pState::speciesName(i)<<"\" \\ \n";
    }
    //Calculated values
    Time_Accurate_File << "\"T\" \\ \n"
		       << "\"R\" \\ \n"
		       << "\"viscosity\" \\ \n"
		       << "\"thermal conduct\" \\ \n"
		       << "\"Prandtl\" \\ \n"
		       << "\"rho*H\"  \\ \n"  
		       <<"\"h\" \\ \n"
		       <<"\"h_s\" \\ \n"
		       <<"\"rho*E\" \\ \n"
		       << "\"e\" \\  \n" 
		       << "\"e_s\" \\ \n";
    for(int i =0; i<Flame2D_pState::NumSpecies(); i++){
      Time_Accurate_File <<"\"omega_c"<<Flame2D_pState::speciesName(i)<<"\" \\ \n";
    }
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
 * Routine: Output_Time_Accurate_File                   *
 *                                                      *
 * This routine writes out progress information for     *
 * a CFD calculation to a progress file, including      *
 * iteration level, time, CPU time, and residual norms. *
 *                                                      *
 ********************************************************/
void Output_to_Time_Accurate_File(ostream &Time_Accurate_File,
				  const double &Time,
				  Flame2D_pState &Soln){

  // compute reaction rates
  Flame2D_State omega;
  omega.Vacuum();
  if (Flame2D_pState::isReacting()) Soln.Sw( omega );

  Soln.updateViscosity();

  Time_Accurate_File << setprecision(6);
  Time_Accurate_File.setf(ios::scientific);
  Time_Accurate_File << " " << Time 
		     << " " << Soln
		     << " " << Soln.T()
		     << " " << Soln.Rtot()
		     << " " << Soln.mu()
		     << " " << Soln.kappa()
		     << " " << Soln.Pr()
		     << " " << Soln.H() 
		     << " " << Soln.h() 
		     << " " << Soln.hs()
		     << " " << ((const Flame2D_pState&)Soln).E() 
		     << " " << Soln.e()
		     << " " << Soln.es();
  for(int k=0; k<Flame2D_pState::NumSpecies(); k++){
    Time_Accurate_File <<" "<<omega.rhoc(k);
  }
  Time_Accurate_File.unsetf(ios::scientific);
  Time_Accurate_File <<endl;
  Time_Accurate_File.flush();
}



/////////////////////////////////////////////////////////////////////
/// Viscous Channel
/////////////////////////////////////////////////////////////////////

/**********************************************************************
 * Routine: Output_Viscous_Channel                                    *
 *                                                                    *
 * Writes the exact and computed viscous channel flow solution values *
 * at the nodes for a 1D array of 2D quadrilateral multi-block        *
 * solution blocks to the specified output data file(s) in a format   *
 * suitable for plotting with TECPLOT.  Returns a non-zero value if   *
 * cannot write any of the TECPLOT solution files.  This routine will *
 * also compute the error norms of the computed solution (L1, L2, and *
 * max norms).                                                        *
 *                                                                    *
 **********************************************************************/
int Output_Viscous_Channel(Flame2D_Quad_Block *Soln_ptr,
			   AdaptiveBlock2D_List &Soln_Block_List,
			   Flame2D_Input_Parameters &IP) {

  int i, i_output_title;
  char prefix[256], extension[256], output_file_name[256];
  char *output_file_name_ptr;
  ofstream output_file;    
  double l1_norm, l2_norm, max_norm;
  double l1_temp, l2_temp, max_temp;
  int numberofactivecells, numberofactivecells_temp;

  // Initialize error variables.
  l1_norm = ZERO; l2_norm = ZERO; max_norm = ZERO;
  l1_temp = ZERO; l2_temp = ZERO; max_temp = ZERO; 
  numberofactivecells = 0; numberofactivecells_temp = 0;

  // Determine prefix of output data file names.
  i = 0;
  while (1) {
    if (IP.Output_File_Name[i] == ' ' ||
	IP.Output_File_Name[i] == '.') break;
    prefix[i] = IP.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(IP.Output_File_Name)) break;
  }
  prefix[i] = '\0';
  strcat(prefix,"_viscous_channel_cpu");
  
  // Determine output data file name for this processor.
  sprintf(extension,"%.6d",Soln_Block_List.ThisCPU);
  strcat(extension,".dat");
  strcpy(output_file_name,prefix);
  strcat(output_file_name,extension);
  output_file_name_ptr = output_file_name;
  
  // Open the output data file.
  output_file.open(output_file_name_ptr,ios::out);
  if (output_file.fail()) return 1;
  
  //double dpdx = -3177.7;  
  double dpdx = IP.Pressure_Gradient;

  // Write the solution data for each solution block.
  i_output_title = 1;
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      l1_temp = ZERO; l2_temp = ZERO; max_temp = ZERO;
      Output_Viscous_Channel(Soln_ptr[nb],nb,i_output_title,output_file,
			     l1_temp,l2_temp,max_temp,IP.Moving_wall_velocity,dpdx);
      l1_norm += l1_temp;
      l2_norm += sqr(l2_temp);
      max_norm = max(max_norm,max_temp);
      numberofactivecells_temp = (Soln_ptr[nb].Grid.NCi-2*Soln_ptr[nb].Nghost)*(Soln_ptr[nb].Grid.NCj-2*Soln_ptr[nb].Nghost);
      numberofactivecells += numberofactivecells_temp;
      if (i_output_title) i_output_title = 0;
    }
  }
  
  // Close the output data file.
  output_file.close();
  
#ifdef _MPI_VERSION
  l1_norm = CFFC_Summation_MPI(l1_norm);
  l2_norm = CFFC_Summation_MPI(l2_norm);
  max_norm = CFFC_Maximum_MPI(max_norm);
  numberofactivecells = CFFC_Summation_MPI(numberofactivecells);
#endif

  // Calculate the L1-norm and L2-norm for all blocks.
  if (Soln_Block_List.Block[0].used == ADAPTIVEBLOCK2D_USED) {
    l1_norm = l1_norm/double(numberofactivecells);
    l2_norm = sqrt(l2_norm/double(numberofactivecells));
  }

  if (CFFC_Primary_MPI_Processor()) {
    cout << endl
	 << endl
	 << " ==================================================================== "
	 << endl
	 << " Error norms for the viscous channel flow:"
	 << endl
	 << "   L1_Norm = " << l1_norm
	 << endl
	 << "   L2_Norm = " << l2_norm
	 << endl
	 << "   Max_Norm = " << max_norm
	 << endl
	 << "   Number of cells = " << numberofactivecells
	 << endl
	 << " ==================================================================== "
	 << endl;
  }

  // Writing of output data files complete.  Return zero value.
  return 0;

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
void Output_Viscous_Channel(Flame2D_Quad_Block &SolnBlk,
			    const int Block_Number,
			    const int Output_Title,
			    ostream &Out_File,
			    double &l1_norm,
			    double &l2_norm,
			    double &max_norm,
			    double &Vwall,
			    const double dpdx) {

  Flame2D_pState We, W;

  // Output node solution data.  
  Out_File << setprecision(14);
  if (Output_Title)
    Out_File << "TITLE = \"" << CFFC_Name() << ": 2D Chem2D Channel Flow Solution "
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
 
  //
  // Loop over the grid
  //
  for (int j = 0; j < SolnBlk.Grid.NCj; j++) {
    for (int i = 0; i < SolnBlk.Grid.NCi; i++) {

      // set the state
      W = SolnBlk.W[i][j];

      // Compute exact solution
      if (i < SolnBlk.Nghost || i > SolnBlk.ICu ||
	  j < SolnBlk.Nghost || j > SolnBlk.JCu) {
	We = SolnBlk.W[i][j];
      } else {
	We.ViscousChannelFlow(SolnBlk.Grid.Cell[i][j].Xc,Vwall,dpdx);
	// compute norms
	double diff = ( ((const Flame2D_pState&)W).vx() - 
			((const Flame2D_pState&)We).vx() );
	l1_norm += fabs( diff );
	l2_norm += sqr( diff );
	max_norm = max(max_norm,fabs( diff ));
      }

      // output results
      const Flame2D_pState& We_ptr = We;
      const Flame2D_pState& W_ptr = W;
      Out_File << SolnBlk.Grid.Cell[i][j].Xc << " " 
	       << W_ptr.rho() << " " 
	       << W_ptr.vx() << " " 
	       << W_ptr.vy() << " " 
	       << We_ptr.vx() << " " 
	       << We_ptr.vy() << " " 
	       << W_ptr.p()   << " " 
	       << W_ptr.T() << " " 
	       << ( W_ptr.vx() - We_ptr.vx() ) << " "
	       << ( W_ptr.vy() - We_ptr.vy() ) << endl;

    } // endfor - i 
  } // endfor - j
  l2_norm = sqrt(l2_norm);

}


/////////////////////////////////////////////////////////////////////
/// Flat Plate
/////////////////////////////////////////////////////////////////////

/**********************************************************************
 * Routine: Output_Flat_Plate                                         *
 *                                                                    *
 * This routine outputs the non-dimensionalized computed flat plate   *
 * solution and the corresponding Blasius solution.                   *
 *                                                                    *
 **********************************************************************/
int Output_Flat_Plate(Flame2D_Quad_Block *Soln_ptr,
		      AdaptiveBlock2D_List &Soln_Block_List,
		      Flame2D_Input_Parameters &IP) {

  int i, i_output_title_soln, i_output_title_skin;
  char prefix_soln[256], prefix_skin[256], extension_soln[256], extension_skin[256];
  char output_file_name_soln[256], output_file_name_skin[256];
  char *output_file_name_soln_ptr, *output_file_name_skin_ptr;
  ofstream output_file_soln, output_file_skin;
  double l1_norm, l2_norm, max_norm;
  double l1_temp, l2_temp, max_temp;
  
  // Initialize error variables.
  l1_norm = ZERO; l2_norm = ZERO; max_norm = ZERO;
  l1_temp = ZERO; l2_temp = ZERO; max_temp = ZERO; 

  // Determine prefix of output data file names.
  i = 0;
  while (1) {
    if (IP.Output_File_Name[i] == ' ' || IP.Output_File_Name[i] == '.') break;
    prefix_soln[i] = IP.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(IP.Output_File_Name)) break;
  }
  prefix_soln[i] = '\0';
  strcat(prefix_soln,"_flatplate_soln_cpu");
  i = 0;
  while (1) {
    if (IP.Output_File_Name[i] == ' ' || IP.Output_File_Name[i] == '.') break;
    prefix_skin[i] = IP.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(IP.Output_File_Name)) break;
  }
  prefix_skin[i] = '\0';
  strcat(prefix_skin,"_flatplate_skin_friction_cpu");
  
  // Determine output data file name for this processor.
  sprintf(extension_soln,"%.6d",Soln_Block_List.ThisCPU);
  strcat(extension_soln,".dat");
  strcpy(output_file_name_soln,prefix_soln);
  strcat(output_file_name_soln,extension_soln);
  output_file_name_soln_ptr = output_file_name_soln;
  sprintf(extension_skin,"%.6d",Soln_Block_List.ThisCPU);
  strcat(extension_skin,".dat");
  strcpy(output_file_name_skin,prefix_skin);
  strcat(output_file_name_skin,extension_skin);
  output_file_name_skin_ptr = output_file_name_skin;
  
  // Open the output data files.
  output_file_soln.open(output_file_name_soln_ptr,ios::out);
  if (output_file_soln.fail()) return 1;
  output_file_skin.open(output_file_name_skin_ptr,ios::out);
  if (output_file_skin.fail()) return 1;
  
  //
  // Write the solution data for each solution block.
  //
  i_output_title_soln = 1;
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {

      // zero
      l1_temp = ZERO; l2_temp = ZERO; max_temp = ZERO; 

      // only output blocks along plate
      Output_Flat_Plate(Soln_ptr[nb],nb,
			i_output_title_soln,output_file_soln,
			i_output_title_skin,output_file_skin,
			IP.Wo,
			l1_temp,l2_temp,max_temp);

      // compute norms
      l1_norm += l1_temp;
      l2_norm += sqr(l2_temp);
      max_norm = max(max_norm,max_temp);

      if (i_output_title_soln) i_output_title_soln = 0;

    }
  }
  
  // Close the output data file.
  output_file_soln.close();
  output_file_skin.close();

  // Calculate the L1-norm and L2-norm for all blocks.
#ifdef _MPI_VERSION
  l1_norm = CFFC_Summation_MPI(l1_norm);
  l2_norm = CFFC_Summation_MPI(l2_norm);
  max_norm = CFFC_Maximum_MPI(max_norm);
#endif
 
  l2_norm = sqrt(l2_norm);

  if (CFFC_Primary_MPI_Processor()) {
    cout << endl
	 << endl
	 << " ==================================================================== "
	 << endl
	 << " Error norms for the skin friction coefficient for the laminar flow"
	 << " over a flat plate:"
	 << endl
	 << "   L1_Norm = " << l1_norm
	 << endl
	 << "   L2_Norm = " << l2_norm
	 << endl
	 << "   Max_Norm = " << max_norm
	 << endl
	 << " ==================================================================== "
	 << endl;
  }

  // Writing of output data files complete.  Return zero value.
  return 0;

}


/**********************************************************************
 * Routine: Output_Flat_Plate                                         *
 *                                                                    *
 * This routine outputs the non-dimensionalized computed flat plate   *
 * solution and the corresponding Blasius solution.                   *
 *                                                                    *
 **********************************************************************/
void Output_Flat_Plate(Flame2D_Quad_Block &SolnBlk,
		       const int Block_Number,
		       const int Output_Title_Soln,
		       ostream &Out_File_Soln,
		       const int Output_Title_Skin,
		       ostream &Out_File_Skin,
		       Flame2D_pState &Winf,
		       double &l1_norm,
		       double &l2_norm,
		       double &max_norm){

  Flame2D_pState W, We;
  const Flame2D_pState &Winf_ptr = Winf;
  Vector2D X;
  double eta, f, fp, fpp, Rex, linf, Cf, Cfe, xpt;   //NEEDS LOGIC TO ONLY OUTPUT BLOCKS ALONG PLATE !!!!!!
  Winf.updateViscosity();

  //-------------------------------------------------------------------
  // Output node solution data.  
  //-------------------------------------------------------------------
  // output title
  Out_File_Soln << setprecision(14);
  if (Output_Title_Soln) {
    Out_File_Soln << "TITLE = \"" << CFFC_Name() << ": 2D Chem2D Solution, "
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
  
    //
    // loop over blocks
    //
    bool output = true;
    for (int j = SolnBlk.Grid.JNl; j <= SolnBlk.Grid.JNu; j++) {
      for (int i = SolnBlk.Grid.INl; i <= SolnBlk.Grid.INu-1; i++) {
	
	// Get cell position and solution data.
	X = SolnBlk.Grid.Node[i][j].X;
	W = SolnBlk.Wn(i,j);

	// is this a point in the boudary layer
	if (fabs(X.x - 0.20*double(npts)) < fabs(X.x - SolnBlk.Grid.Node[i+1][j].X.x) &&
	    X.x <= 0.20*double(npts) ) {

	  // output zone title
	  if (output) {
	    Out_File_Soln << "ZONE T =  \"Block Number = " << Block_Number
			  << ", xpt = " << 2000.0*double(npts)
			  << "\" \\ \n"
			  << "I = " << 1 << " \\ \n"
			  << "J = " << SolnBlk.Grid.JNu - SolnBlk.Grid.JNl + 1 << " \\ \n"
			  << "F = POINT \n";
	    output = false;
	  }

	  // compute flat plate solution
	  We.FlatPlate(Winf,X,eta,f,fp,fpp);
	  linf = 10000.0*(Winf.mu()/Winf_ptr.rho())/Winf_ptr.vx();
	  Rex = (Winf_ptr.vx()/(Winf.mu()/Winf_ptr.rho()))*((X.x+NANO)/ONE);
	  
	  // Output data.
	  const Flame2D_pState& We_ptr = We;
	  const Flame2D_pState& W_ptr = W;
	  Out_File_Soln.setf(ios::scientific);
	  Out_File_Soln << " " << X
			<< " " << W_ptr.rho() 
			<< W_ptr.v()
			<< " " << W_ptr.p()
			<< " " << We_ptr.rho()
			<< We_ptr.v()
			<< " " << We_ptr.p()
			<< " " << eta
			<< " " << f
			<< " " << fp
			<< " " << fpp
			<< " " << W_ptr.vx()/Winf_ptr.vx()
			<< " " << We_ptr.vx()/Winf_ptr.vx()
			<< " " << (W_ptr.vy()/Winf_ptr.vx())*sqrt(Winf_ptr.vx()/((Winf.mu()/Winf_ptr.rho())*(X.x+NANO)))
			<< " " << HALF*(eta*fp-f)
			<< endl;

	} // endif - relevant point

      } // endfor - i
    } // endfor - j

  } // endfor - points

    //-------------------------------------------------------------------
    // Output skin friction coefficient.
    //-------------------------------------------------------------------
    // only output for blocks along the plate
  if (SolnBlk.Grid.BCtypeS[SolnBlk.Nghost] == BC_WALL_VISCOUS_HEATFLUX) {
    
    // output title
    if (Output_Title_Skin) {
      Out_File_Skin << "TITLE = \"" << CFFC_Name() << ": 2D Chem2D Solution, "
		    << "\"" << "\n"
		    << "VARIABLES = \"x\" \\ \n"
		    << "\"Rex\" \\ \n"
		    << "\"Cf\" \\ \n"
		    << "\"Cf_e\" \\ \n";
    }
    Out_File_Skin << "ZONE T =  \"Block Number = " << Block_Number
		  << "\" \\ \n"
		  << "I = " << SolnBlk.Grid.INu - SolnBlk.Grid.INl << " \\ \n"
		  << "J = " << 1 << " \\ \n"
		  << "F = POINT \n";

    //
    // Loop over the grid
    //
    for (int j = SolnBlk.Grid.JCl; j <= SolnBlk.Grid.JCu; j++) {
      for (int i = SolnBlk.Grid.ICl; i <= SolnBlk.Grid.ICu; i++) {

	// Get cell position and solution data.
	X = SolnBlk.Grid.Cell[i][j].Xc;
	W = SolnBlk.W[i][j];

	// Determine the skin friction coefficient.
	linf = 10000.0*(Winf.mu()/Winf_ptr.rho())/Winf_ptr.vx();
	Rex = (Winf_ptr.vx()/(Winf.mu()/Winf_ptr.rho()))*((X.x+NANO)/ONE);

	if (X.x >= ZERO && j == 2 && SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX){
	     
	  // Get exact skin friction coefficient.
	  Cfe = TWO*0.32206/sqrt(Rex);
	  // Get computed skin friction coefficient.
	  Cf  = TWO*W.WallShearStress(X,
				      SolnBlk.Grid.nodeSW(i,j).X,
				      SolnBlk.Grid.nodeSE(i,j).X,
				      Vector2D(ONE,ZERO))/(Winf_ptr.rho()*Winf_ptr.vx()*Winf_ptr.vx());
	  // Output data.
	  Out_File_Skin.setf(ios::scientific);
	  Out_File_Skin << " " << X.x
			<< " " << Rex*linf
			<< " " << Cf
			<< " " << Cfe
			<< endl;
	  // Calculate error norms.
	  l1_norm += fabs(Cf - Cfe);
	  l2_norm += sqr(Cf - Cfe);
	  max_norm = max(max_norm,fabs(Cf - Cfe));
	  
	} // endif
	
      } // endfor - i
    } // endfor - j
    
    l2_norm = sqrt(l2_norm);

  } // endif - plate block
  
}
 

/////////////////////////////////////////////////////////////////////
/// Driven Cavity Flow
/////////////////////////////////////////////////////////////////////

/**********************************************************************
 * Routine: Output_Driven_Cavity_Flow                                 *
 *                                                                    *
 * This routine outputs a comparison of the computed solution for the *
 * driven cavity flow with the computations done by Ghia et al. (J.   *
 * Comp. Phys. Vol. 48 1982) for a 1D array of 2D quadrilateral       *
 * multi-block solution blocks to the specified output data file(s)   *
 * in a format suitable for plotting with tecplot.                    *
 *                                                                    *
 **********************************************************************/
int Output_Driven_Cavity_Flow(Flame2D_Quad_Block *Soln_ptr,
			      AdaptiveBlock2D_List &Soln_Block_List,
			      Flame2D_Input_Parameters &IP) {

  int i, i_output_title;
  char prefix_u[256], extension_u[256], output_file_name_u[256];
  char *output_file_name_ptr_u;
  ofstream output_file_u;
  char prefix_v[256], extension_v[256], output_file_name_v[256];
  char *output_file_name_ptr_v;
  ofstream output_file_v;

  // Determine prefix of output data file names.
  i = 0;
  while (1) {
    if (IP.Output_File_Name[i] == ' ' || IP.Output_File_Name[i] == '.') break;
    prefix_u[i] = IP.Output_File_Name[i];
    prefix_v[i] = IP.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(IP.Output_File_Name)) break;
  }
  prefix_u[i] = '\0';
  prefix_v[i] = '\0';
  strcat(prefix_u,"_driven_cavity_flow_u_cpu");
  strcat(prefix_v,"_driven_cavity_flow_v_cpu");

  // Determine output data file name for this processor.
  sprintf(extension_u,"%.6d",Soln_Block_List.ThisCPU);
  sprintf(extension_v,"%.6d",Soln_Block_List.ThisCPU);
  strcat(extension_u,".dat");
  strcat(extension_v,".dat");
  strcpy(output_file_name_u,prefix_u);
  strcpy(output_file_name_v,prefix_v);
  strcat(output_file_name_u,extension_u);
  strcat(output_file_name_v,extension_v);
  output_file_name_ptr_u = output_file_name_u;
  output_file_name_ptr_v = output_file_name_v;

  // Open the output data files.
  output_file_u.open(output_file_name_ptr_u,ios::out);
  if (output_file_u.fail()) return 1;
  output_file_v.open(output_file_name_ptr_v,ios::out);
  if (output_file_v.fail()) return 2;

  // Write the solution data for each solution block.
  i_output_title = 1;
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      Output_Driven_Cavity_Flow(Soln_ptr[nb],nb,
				i_output_title,
				output_file_u,
				output_file_v,
				IP.Re_lid,
				IP.Moving_wall_velocity,
				IP.Box_Width);
      if (i_output_title) i_output_title = 0;
    }
  }
  
  // Close the output data files.
  output_file_u.close();
  output_file_v.close();

  // Writing of output data files complete.  Return zero value.
  return 0;

}


/**********************************************************************
 * Routine: Output_Driven_Cavity_Flow                                 *
 *                                                                    *
 * This routine writes a comparison of the computed solution for the  *
 * driven cavity flow with the computations done by Ghia et al. (J.   *
 * Comp. Phys. Vol. 48 1982) for the specified quadrilateral solution *
 * block to the specified output stream suitable for plotting with    *
 * TECPLOT.                                                           *
 *                                                                    *
 **********************************************************************/
void Output_Driven_Cavity_Flow(Flame2D_Quad_Block &SolnBlk,
			       const int Block_Number,
			       const int Output_Title,
			       ostream &Out_File_u,
			       ostream &Out_File_v,
			       const double &Re,
			       const double &Vwall,
			       const double &length) {

  // Set the output precision and type.
  Out_File_u << setprecision(14);
  Out_File_u.setf(ios::scientific);
  Out_File_v << setprecision(14);
  Out_File_v.setf(ios::scientific);

  // Output node solution data.  
  if (Output_Title) {
    Out_File_u << "TITLE = \"" << CFFC_Name() << ": 2D Dusty Driven Cavity Flow Solution u-velocity Comparison "
	       << "\"" << "\n"
	       << "VARIABLES = \"x\" \\ \n"
	       << "\"u\" \\ \n";
    Out_File_v << "TITLE = \"" << CFFC_Name() << ": 2D Dusty Driven Cavity Flow Solution v-velocity Comparison "
	       << "\"" << "\n"
	       << "VARIABLES = \"x\" \\ \n"
	       << "\"v\" \\ \n";
  }

  // If the output title flag is on (first time the current processor
  // is writing data) and the current processor is the primary processor
  // then write the solution data of Ghia et al. to the output file.
  if (Output_Title && CFFC_Primary_MPI_Processor()) {
    // Output a zone corresponding to the u-velocity component at the
    // geometric centre of the cavity.
    Out_File_u << "ZONE T = \"u-velocity"
	       << "\" \\ \n"
	       << "I = " << 17 << " \\ \n"
	       << "J = " << 1 << " \\ \n"
	       << "F = POINT \\ \n";
    if (Re == 100.0) {
      Out_File_u << 0.0000 << " " <<  0.00000 << endl
		 << 0.0547 << " " << -0.03717 << endl
		 << 0.0625 << " " << -0.04192 << endl
		 << 0.0703 << " " << -0.04775 << endl
		 << 0.1016 << " " << -0.06434 << endl
		 << 0.1719 << " " << -0.10150 << endl
		 << 0.2813 << " " << -0.15662 << endl
		 << 0.4531 << " " << -0.21090 << endl
		 << 0.5000 << " " << -0.20581 << endl
		 << 0.6172 << " " << -0.13641 << endl
		 << 0.7344 << " " <<  0.00332 << endl
		 << 0.8516 << " " <<  0.23151 << endl
		 << 0.9531 << " " <<  0.68717 << endl
		 << 0.9609 << " " <<  0.73722 << endl
		 << 0.9688 << " " <<  0.78871 << endl
		 << 0.9766 << " " <<  0.84123 << endl
		 << 1.0000 << " " <<  1.00000 << endl;
    } else if (Re == 400.0) {
      Out_File_u << 0.0000 << " " <<  0.00000 << endl
		 << 0.0547 << " " << -0.08186 << endl
		 << 0.0625 << " " << -0.09266 << endl
		 << 0.0703 << " " << -0.10338 << endl
		 << 0.1016 << " " << -0.14612 << endl
		 << 0.1719 << " " << -0.24299 << endl
		 << 0.2813 << " " << -0.32726 << endl
		 << 0.4531 << " " << -0.17119 << endl
		 << 0.5000 << " " << -0.11477 << endl
		 << 0.6172 << " " <<  0.02135 << endl
		 << 0.7344 << " " <<  0.16256 << endl
		 << 0.8516 << " " <<  0.29093 << endl
		 << 0.9531 << " " <<  0.55892 << endl
		 << 0.9609 << " " <<  0.61756 << endl
		 << 0.9688 << " " <<  0.68439 << endl
		 << 0.9766 << " " <<  0.75837 << endl
		 << 1.0000 << " " <<  1.00000 << endl;
    }
    // Output a zone corresponding to the v-velocity component at the
    // geometric centre of the cavity.
    Out_File_v << "ZONE T = \"v-velocity"
	       << "\" \\ \n"
	       << "I = " << 17 << " \\ \n"
	       << "J = " << 1 << " \\ \n"
	       << "F = POINT \\ \n";
    if (Re == 100.0) {
      Out_File_v << 0.0000 << " " <<  0.00000 << endl
		 << 0.0625 << " " <<  0.09233 << endl
		 << 0.0703 << " " <<  0.10091 << endl
		 << 0.0781 << " " <<  0.10890 << endl
		 << 0.0938 << " " <<  0.12317 << endl
		 << 0.1563 << " " <<  0.16077 << endl
		 << 0.2266 << " " <<  0.17507 << endl
		 << 0.2344 << " " <<  0.17527 << endl
		 << 0.5000 << " " <<  0.05454 << endl
		 << 0.8047 << " " << -0.24533 << endl
		 << 0.8594 << " " << -0.22445 << endl
		 << 0.9063 << " " << -0.16914 << endl
		 << 0.9453 << " " << -0.10313 << endl
		 << 0.9531 << " " << -0.08864 << endl
		 << 0.9609 << " " << -0.07391 << endl
		 << 0.9688 << " " << -0.05906 << endl
		 << 1.0000 << " " <<  0.00000 << endl;
    } else if (Re == 400.0) {
      Out_File_v << 0.0000 << " " <<  0.00000 << endl
		 << 0.0625 << " " <<  0.18360 << endl
		 << 0.0703 << " " <<  0.19713 << endl
		 << 0.0781 << " " <<  0.20920 << endl
		 << 0.0938 << " " <<  0.22965 << endl
		 << 0.1563 << " " <<  0.28124 << endl
		 << 0.2266 << " " <<  0.30203 << endl
		 << 0.2344 << " " <<  0.30174 << endl
		 << 0.5000 << " " <<  0.05186 << endl
		 << 0.8047 << " " << -0.38598 << endl
		 << 0.8594 << " " << -0.44993 << endl
		 << 0.9063 << " " << -0.33827 << endl
		 << 0.9453 << " " << -0.22847 << endl
		 << 0.9531 << " " << -0.19254 << endl
		 << 0.9609 << " " << -0.15663 << endl
		 << 0.9688 << " " << -0.12146 << endl
		 << 1.0000 << " " <<  0.00000 << endl;
    }
  }

  if (fabs(SolnBlk.Grid.nodeSE(SolnBlk.ICu,SolnBlk.JCl).X.x) < TOLER*TOLER) {
    Out_File_u << "ZONE T =  \"u-velocity " << Block_Number
	       << "\" \\ \n"
	       << "I = " << (SolnBlk.JCu+1)-(SolnBlk.JCl-1) << " \\ \n"
	       << "J = " << 1 << " \\ \n"
	       << "F = POINT \\ \n";
    for (int j = SolnBlk.JCl; j <= SolnBlk.JCu+1; j++) {
      Out_File_u << SolnBlk.Grid.nodeSE(SolnBlk.ICu,j).X.y/length + HALF << " "
		 << SolnBlk.WnSE(SolnBlk.ICu,j).vx()/Vwall << endl;
    }
  }
  if (fabs(SolnBlk.Grid.nodeSW(SolnBlk.ICl,SolnBlk.JCl).X.y) < TOLER*TOLER) {
    Out_File_v << "ZONE T =  \"v-velocity " << Block_Number
	       << "\" \\ \n"
	       << "I = " << (SolnBlk.ICu+1)-(SolnBlk.ICl-1) << " \\ \n"
	       << "J = " << 1 << " \\ \n"
	       << "F = POINT \\ \n";
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu+1; i++) {
      Out_File_v << SolnBlk.Grid.nodeSW(i,SolnBlk.JCl).X.x/length + HALF << " "
		 << SolnBlk.WnSW(i,SolnBlk.JCl).vy()/Vwall << endl;
    }
  }

}

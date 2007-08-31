/**********************************************************************
 * NavierStokes2DQuadIOSingleBlock.cc: Single-block versions of input *
 *                                     and output subroutines for 2D  *
 *                                     Navier-Stokes multi-block      *
 *                                     quadrilateral mesh solution    *
 *                                     classes.                       *
 **********************************************************************/

// Include 2D NavierStokes quadrilateral mesh solution header file.

#ifndef _NAVIERSTOKES2D_QUAD_INCLUDED
#include "NavierStokes2DQuad.h"
#endif // _NAVIERSTOKES2D_QUAD_INCLUDED

/**********************************************************************
 * NavierStokes2D_Quad_Block -- IO Single Block External Subroutines. *
 **********************************************************************/

/**********************************************************************
 * Routine: Write_Solution_Block                                      *
 *                                                                    *
 * Writes the cell centred solution values of the specified           *
 * quadrilateral solution block to the specified output stream for    *
 * restart purposes.                                                  *
 *                                                                    *
 **********************************************************************/
void Write_Solution_Block(NavierStokes2D_Quad_Block &SolnBlk,
	                  ostream &Out_File) {

  Out_File << setprecision(14) << SolnBlk << setprecision(6);

}

/**********************************************************************
 * Routine: Read_Solution_Block                                       *
 *                                                                    *
 * Reads the cell centred solution values for the specified           *
 * quadrilateral solution block from the specified input stream as    *
 * required for restart purposes.                                     *
 *                                                                    *
 **********************************************************************/
void Read_Solution_Block(NavierStokes2D_Quad_Block &SolnBlk,
	                 istream &In_File) {

  In_File >> SolnBlk;

}

/**********************************************************************
 * Routine: Output_Tecplot                                            *
 *                                                                    *
 * Writes the solution values at the nodes of the specified           *
 * quadrilateral solution block to the specified output stream        *
 * suitable for plotting with TECPLOT.                                *
 *                                                                    *
 **********************************************************************/
void Output_Tecplot(NavierStokes2D_Quad_Block &SolnBlk,
		    NavierStokes2D_Input_Parameters &IP,
		    const int Number_of_Time_Steps,
		    const double &Time,
		    const int Block_Number,
		    const int Output_Title,
		    ostream &Out_File) {

  NavierStokes2D_pState W_node;

  // Ensure boundary conditions are updated before evaluating
  // solution at the nodes.
  BCs(SolnBlk,IP);

  // Output node solution data.  
  Out_File << setprecision(14);
  if (Output_Title) {
    Out_File << "TITLE = \"" << CFFC_Name() << ": 2D NavierStokes Solution, "
	     << "Time Step/Iteration Level = " << Number_of_Time_Steps
	     << ", Time = " << Time
	     << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
	     << "\"y\" \\ \n";
    Out_File << "\"rho\" \\ \n"
	     << "\"u\" \\ \n"
	     << "\"v\" \\ \n"
	     << "\"p\" \\ \n"
	     << "\"T\" \\ \n"
	     << "\"M\" \\ \n"
	     << "\"H\" \\ \n"
	     << "\"s\" \\ \n";
    if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
      Out_File << "\"k\" \\ \n"
	       << "\"omega\" \\ \n"
	       << "\"epsilon\" \\ \n"
	       << "\"ell\" \\ \n"
	       << "\"p_modified\" \\ \n";
    }
    Out_File << "\"Rex\" \\ \n";
  }

  Out_File << "ZONE T =  \"Block Number = " << Block_Number
	   << "\" \\ \n"
	   << "I = " << SolnBlk.Grid.INu - SolnBlk.Grid.INl + 1 << " \\ \n"
	   << "J = " << SolnBlk.Grid.JNu - SolnBlk.Grid.JNl + 1 << " \\ \n"
	   << "F = POINT \\ \n";

  for (int j = SolnBlk.Grid.JNl; j <= SolnBlk.Grid.JNu; j++) {
    for (int i = SolnBlk.Grid.INl; i <= SolnBlk.Grid.INu; i++) {
      Out_File.setf(ios::scientific);
      Out_File << " " << SolnBlk.Grid.Node[i][j].X;
      W_node = SolnBlk.Wn(i,j);
      Out_File << " " << W_node.rho 
	       << W_node.v
	       << " " << W_node.p
	       << " " << W_node.T()
	       << " " << W_node.v.abs()/W_node.a() 
	       << " " << W_node.H()
	       << " " << W_node.s();
      if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
	Out_File << " " << W_node.k
		 << " " << W_node.omega
		 << " " << W_node.epsilon()
		 << " " << W_node.ell()
		 << " " << W_node.pmodified();
      }
      Out_File << " " << IP.Wo.v.x/IP.Wo.nu()*SolnBlk.Grid.Node[i][j].X.x;
      Out_File << endl;
    }
  }

  Out_File << setprecision(6);

}

/**********************************************************************
 * Routine: Output_Cells_Tecplot                                      *
 *                                                                    *
 * Writes the cell centred solution values of the specified           *
 * quadrilateral solution block to the specified output stream        *
 * suitable for plotting with TECPLOT.                                *
 *                                                                    *
 **********************************************************************/
void Output_Cells_Tecplot(NavierStokes2D_Quad_Block &SolnBlk,
		          NavierStokes2D_Input_Parameters &IP,
                          const int Number_of_Time_Steps,
                          const double &Time,
                          const int Block_Number,
                          const int Output_Title,
	                  ostream &Out_File) {

  int nghost = 1;

  Out_File << setprecision(14);
  if (Output_Title) {
    Out_File << "TITLE = \"" << CFFC_Name() << ": 2D NavierStokes Solution, "
	     << "Time Step/Iteration Level = " << Number_of_Time_Steps
	     << ", Time = " << Time
	     << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
	     << "\"y\" \\ \n"
	     << "\"rho\" \\ \n"
	     << "\"u\" \\ \n"
	     << "\"v\" \\ \n"
	     << "\"p\" \\ \n"
	     << "\"T\" \\ \n"
	     << "\"M\" \\ \n"
	     << "\"H\" \\ \n"
	     << "\"s\" \\ \n";
    if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
      Out_File << "\"k\" \\ \n"
	       << "\"omega\" \\ \n"
	       << "\"epsilon\" \\ \n"
	       << "\"ell\" \\ \n"
	       << "\"p_modified\" \\ \n"
	       << "\"Mt\" \\ \n";
    }
    if (SolnBlk.Flow_Type) {
      Out_File << "\"tau_xx\" \\ \n"
	       << "\"tau_xy\" \\ \n"
	       << "\"tau_yy\" \\ \n";
      if (SolnBlk.Axisymmetric)
	Out_File << "\"tau_zz\" \\ \n";
      Out_File << "\"qx\" \\ \n"
	       << "\"qy\" \\ \n";
    }
    Out_File << "\"mu\" \\ \n";
    if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
      Out_File << "\"muT\" \\ \n"
	       << "\"ywall\" \\ \n"
	       << "\"Xwall.x\" \\ \n"
	       << "\"Xwall.y\" \\ \n"
	       << "\"nwall.x\" \\ \n"
	       << "\"nwall.y\" \\ \n"
	       << "\"BCwall\" \\ \n"
	       << "\"yplus\" \\ \n"
	       << "\"utau\" \\ \n"
	       << "\"tauw\" \\ \n"
	       << "\"vinj\" \\ \n"
	       << "\"vwplus\" \\ \n";
    }
    Out_File << "\"dt\" \\ \n";
    Out_File << "\"A\" \\ \n";
    Out_File << "\"AR\" \\ \n";
  }
  Out_File << "ZONE T =  \"Block Number = " << Block_Number
	   << "\" \\ \n"
	   << "I = " << SolnBlk.Grid.ICu - SolnBlk.Grid.ICl + 2*nghost + 1 << " \\ \n"
	   << "J = " << SolnBlk.Grid.JCu - SolnBlk.Grid.JCl + 2*nghost + 1 << " \\ \n"
	   << "F = POINT \n";

  for (int j = SolnBlk.JCl-nghost; j <= SolnBlk.JCu+nghost; j++) {
    for (int i = SolnBlk.ICl-nghost; i <= SolnBlk.ICu+nghost; i++) {
      if (SolnBlk.Flow_Type != FLOWTYPE_INVISCID) {
 	SolnBlk.W[i][j].ComputeViscousTerms(SolnBlk.dWdx[i][j],
 					    SolnBlk.dWdy[i][j],
 					    SolnBlk.Grid.Cell[i][j].Xc,
 					    SolnBlk.Axisymmetric);
 	SolnBlk.U[i][j].tau = SolnBlk.W[i][j].tau;
 	SolnBlk.U[i][j].q = SolnBlk.W[i][j].q;
      }
      Out_File.setf(ios::scientific);
      Out_File << " " << SolnBlk.Grid.Cell[i][j].Xc;
      Out_File << " " << SolnBlk.W[i][j].rho 
	       << SolnBlk.W[i][j].v
	       << " " << SolnBlk.W[i][j].p
	       << " " << SolnBlk.W[i][j].T() 
	       << " " << SolnBlk.W[i][j].v.abs()/SolnBlk.W[i][j].a() 
	       << " " << SolnBlk.W[i][j].H() 
	       << " " << SolnBlk.W[i][j].s();
      if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
	Out_File << " " << SolnBlk.W[i][j].k
		 << " " << SolnBlk.W[i][j].omega
 		 << " " << SolnBlk.W[i][j].epsilon()
		 << " " << SolnBlk.W[i][j].ell()
		 << " " << SolnBlk.W[i][j].pmodified()
		 << " " << SolnBlk.W[i][j].Mt();
      }
      if (SolnBlk.Flow_Type) {
	Out_File << " " << SolnBlk.W[i][j].tau.xx
		 << " " << SolnBlk.W[i][j].tau.xy
		 << " " << SolnBlk.W[i][j].tau.yy;
	if (SolnBlk.Axisymmetric)
	  Out_File << " " << SolnBlk.W[i][j].tau.zz;
	Out_File << " " << SolnBlk.W[i][j].q.x
		 << " " << SolnBlk.W[i][j].q.y;
      }
      Out_File.unsetf(ios::scientific);
      Out_File << " " << SolnBlk.W[i][j].mu();
      if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
	Out_File << " " << SolnBlk.W[i][j].muT()
		 << " " << SolnBlk.Wall[i][j].ywall
		 << SolnBlk.Wall[i][j].Xwall
		 << SolnBlk.Wall[i][j].nwall
		 << " " << SolnBlk.Wall[i][j].BCwall
		 << " " << SolnBlk.Wall[i][j].yplus
		 << " " << SolnBlk.Wall[i][j].utau
		 << " " << SolnBlk.Wall[i][j].tauw
		 << " " << SolnBlk.Wall[i][j].vinj
		 << " " << SolnBlk.Wall[i][j].vwplus;
      }
      Out_File << " " << SolnBlk.dt[i][j];
      Out_File << " " << SolnBlk.Grid.Cell[i][j].A;
      Out_File << " " << (abs(SolnBlk.Grid.xfaceE(i,j) - SolnBlk.Grid.xfaceW(i,j))/
			  abs(SolnBlk.Grid.xfaceN(i,j) - SolnBlk.Grid.xfaceS(i,j)));
      Out_File << endl;
    }
  }
  Out_File << setprecision(6);

}

void Output_Nozzleless_Tecplot(NavierStokes2D_Quad_Block &SolnBlk,
			       const int Number_of_Time_Steps,
			       const double &Time,
			       const int Block_Number,
			       const int Output_Title,
			       ostream &Out_File,
			       const double &po,
			       const double &rhoo) {

  int nghost = 1;

  Out_File << setprecision(14);
  if (Output_Title) {
    Out_File << "TITLE = \"" << CFFC_Name() << ": 2D NavierStokes Solution, "
	     << "Time Step/Iteration Level = " << Number_of_Time_Steps
	     << ", Time = " << Time
	     << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
	     << "\"y\" \\ \n"
	     << "\"rho\" \\ \n"
	     << "\"u\" \\ \n"
	     << "\"v\" \\ \n"
	     << "\"p\" \\ \n"
	     << "\"T\" \\ \n"
	     << "\"M\" \\ \n"
	     << "\"H\" \\ \n"
	     << "\"s\" \\ \n";
    if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
      Out_File << "\"k\" \\ \n"
	       << "\"omega\" \\ \n"
	       << "\"epsilon\" \\ \n"
	       << "\"ell\" \\ \n"
	       << "\"p_modified\" \\ \n"
	       << "\"Mt\" \\ \n";
    }
    if (SolnBlk.Flow_Type) {
      Out_File << "\"tau_xx\" \\ \n"
	       << "\"tau_xy\" \\ \n"
	       << "\"tau_yy\" \\ \n";
      if (SolnBlk.Axisymmetric)
	Out_File << "\"tau_zz\" \\ \n";
      Out_File << "\"qx\" \\ \n"
	       << "\"qy\" \\ \n";
    }
    Out_File << "\"mu\" \\ \n";
    if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
      Out_File << "\"muT\" \\ \n"
	       << "\"ywall\" \\ \n"
	       << "\"Xwall.x\" \\ \n"
	       << "\"Xwall.y\" \\ \n"
	       << "\"nwall.x\" \\ \n"
	       << "\"nwall.y\" \\ \n"
	       << "\"BCwall\" \\ \n"
	       << "\"yplus\" \\ \n"
	       << "\"utau\" \\ \n"
	       << "\"tauw\" \\ \n"
	       << "\"vinj\" \\ \n"
	       << "\"vwplus\" \\ \n";
    }
    Out_File << "\"dt\" \\ \n";
    Out_File << "\"A\" \\ \n";
    Out_File << "\"AR\" \\ \n";
    Out_File << "\"xbar\" \\ \n"
	     << "\"ybar\" \\ \n"
	     << "\"rhobar\" \\ \n"
	     << "\"pbar\" \\ \n"
	     << "\"mdot\" \\ \n";
  }
  Out_File << "ZONE T =  \"Block Number = " << Block_Number
	   << "\" \\ \n"
	   << "I = " << SolnBlk.Grid.ICu - SolnBlk.Grid.ICl + 2*nghost + 1 << " \\ \n"
	   << "J = " << SolnBlk.Grid.JCu - SolnBlk.Grid.JCl + 2*nghost + 1 << " \\ \n"
	   << "F = POINT \n";

  for (int j = SolnBlk.JCl-nghost; j <= SolnBlk.JCu+nghost; j++) {
    for (int i = SolnBlk.ICl-nghost; i <= SolnBlk.ICu+nghost; i++) {
      if (SolnBlk.Flow_Type != FLOWTYPE_INVISCID) {
 	SolnBlk.W[i][j].ComputeViscousTerms(SolnBlk.dWdx[i][j],
 					    SolnBlk.dWdy[i][j],
 					    SolnBlk.Grid.Cell[i][j].Xc,
 					    SolnBlk.Axisymmetric);
 	SolnBlk.U[i][j].tau = SolnBlk.W[i][j].tau;
 	SolnBlk.U[i][j].q = SolnBlk.W[i][j].q;
      }
      Out_File.setf(ios::scientific);
      Out_File << " " << SolnBlk.Grid.Cell[i][j].Xc;
      Out_File << " " << SolnBlk.W[i][j].rho 
	       << SolnBlk.W[i][j].v
	       << " " << SolnBlk.W[i][j].p
	       << " " << SolnBlk.W[i][j].T() 
	       << " " << SolnBlk.W[i][j].v.abs()/SolnBlk.W[i][j].a() 
	       << " " << SolnBlk.W[i][j].H() 
	       << " " << SolnBlk.W[i][j].s();
      if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
	Out_File << " " << SolnBlk.W[i][j].k
		 << " " << SolnBlk.W[i][j].omega
 		 << " " << SolnBlk.W[i][j].epsilon()
		 << " " << SolnBlk.W[i][j].ell()
		 << " " << SolnBlk.W[i][j].pmodified()
		 << " " << SolnBlk.W[i][j].Mt();
      }
      if (SolnBlk.Flow_Type) {
	Out_File << " " << SolnBlk.W[i][j].tau.xx
		 << " " << SolnBlk.W[i][j].tau.xy
		 << " " << SolnBlk.W[i][j].tau.yy;
	if (SolnBlk.Axisymmetric)
	  Out_File << " " << SolnBlk.W[i][j].tau.zz;
	Out_File << " " << SolnBlk.W[i][j].q.x
		 << " " << SolnBlk.W[i][j].q.y;
      }
      Out_File.unsetf(ios::scientific);
      Out_File << " " << SolnBlk.W[i][j].mu();
      if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
	Out_File << " " << SolnBlk.W[i][j].muT()
		 << " " << SolnBlk.Wall[i][j].ywall
		 << SolnBlk.Wall[i][j].Xwall
		 << SolnBlk.Wall[i][j].nwall
		 << " " << SolnBlk.Wall[i][j].BCwall
		 << " " << SolnBlk.Wall[i][j].yplus
		 << " " << SolnBlk.Wall[i][j].utau
		 << " " << SolnBlk.Wall[i][j].tauw
		 << " " << SolnBlk.Wall[i][j].vinj
		 << " " << SolnBlk.Wall[i][j].vwplus;
      }
      Out_File << " " << SolnBlk.dt[i][j];
      Out_File << " " << SolnBlk.Grid.Cell[i][j].A;
      Out_File << " " << (abs(SolnBlk.Grid.xfaceE(i,j) - SolnBlk.Grid.xfaceW(i,j))/
			  abs(SolnBlk.Grid.xfaceN(i,j) - SolnBlk.Grid.xfaceS(i,j)));
      Out_File << " " << SolnBlk.Grid.Cell[i][j].Xc.x/0.01
	       << " " << SolnBlk.Grid.Cell[i][j].Xc.y/0.01
	       << " " << SolnBlk.W[i][j].rho/rhoo
	       << " " << SolnBlk.W[i][j].p/po
	       << " " << SolnBlk.W[i][j].rho*SolnBlk.W[i][j].v.y;
      Out_File << endl;
    }
  }
  Out_File << setprecision(6);

}

/**********************************************************************
 * Routine: Output_Nodes_Tecplot                                      *
 *                                                                    *
 * Writes the node values of the specified quadrilateral solution     *
 * block to the specified output stream suitable for plotting with    *
 * TECPLOT.                                                           *
 *                                                                    *
 **********************************************************************/
void Output_Nodes_Tecplot(NavierStokes2D_Quad_Block &SolnBlk,
                          const int Number_of_Time_Steps,
                          const double &Time,
                          const int Block_Number,
                          const int Output_Title,
	                  ostream &Out_File) {

  Out_File << setprecision(14);
  if (Output_Title) {
    Out_File << "TITLE = \"" << CFFC_Name() << ": 2D NavierStokes Nodes, "
	     << "Time Step/Iteration Level = " << Number_of_Time_Steps
	     << ", Time = " << Time
	     << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
	     << "\"y\" \\ \n";
  }
  Out_File << "ZONE T =  \"Block Number = " << Block_Number
	   << "\" \\ \n"
	   << "I = " << SolnBlk.Grid.INu - SolnBlk.Grid.INl + 1 + 2*SolnBlk.Nghost << " \\ \n"
	   << "J = " << SolnBlk.Grid.JNu - SolnBlk.Grid.JNl + 1 + 2*SolnBlk.Nghost << " \\ \n"
	   << "F = POINT \\ \n";

  Out_File.setf(ios::scientific);
  for (int j = SolnBlk.Grid.JNl - SolnBlk.Nghost; j <= SolnBlk.Grid.JNu + SolnBlk.Nghost; j++) {
    for (int i = SolnBlk.Grid.INl - SolnBlk.Nghost; i <= SolnBlk.Grid.INu + SolnBlk.Nghost; i++) {
      Out_File << " " << SolnBlk.Grid.Node[i][j].X << endl;
    }
  }
  Out_File.unsetf(ios::scientific);
  Out_File << setprecision(6);

}

/**********************************************************************
 * Routine: Output_Gradients_Tecplot                                  *
 *                                                                    *
 * Writes the cell centred primitive variable state gradients and     *
 * limiters of the specified quadrilateral solution block to the      *
 * specified output stream suitable for plotting with TECPLOT.        *
 *                                                                    *
 **********************************************************************/
void Output_Gradients_Tecplot(NavierStokes2D_Quad_Block &SolnBlk,
			      const int Number_of_Time_Steps,
			      const double &Time,
			      const int Block_Number,
			      const int Output_Title,
			      ostream &Out_File) {

  Out_File << setprecision(14);
  if (Output_Title) {
    Out_File << "TITLE = \"" << CFFC_Name() << ": 2D NavierStokes Solution, "
	     << "Time Step/Iteration Level = " << Number_of_Time_Steps
	     << ", Time = " << Time
	     << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
 	     << "\"y\" \\ \n"
 	     << "\"drhodx\" \\ \n"
 	     << "\"drhody\" \\ \n"
 	     << "\"dudx\" \\ \n"
 	     << "\"dudy\" \\ \n"
 	     << "\"dvdx\" \\ \n"
 	     << "\"dvdy\" \\ \n"
 	     << "\"dpdx\" \\ \n"
 	     << "\"dpdy\" \\ \n";
    if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
      Out_File << "\"dkdx\" \\ \n"
	       << "\"dkdy\" \\ \n"
	       << "\"domegadx\" \\ \n"
	       << "\"domgeady\" \\ \n";
    }
    Out_File << "\"dilatation\" \\ \n"
	     << "\"divergence\" \\ \n"
	     << "\"vorticity\" \\ \n";
    Out_File << "\"phi_rho\" \\ \n"
 	     << "\"phi_u\" \\ \n"
 	     << "\"phi_v\" \\ \n"
 	     << "\"phi_p\" \\ \n";
    if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
      Out_File << "\"phi_k\" \\ \n"
	       << "\"phi_omega\" \\ \n";
    }
  }
  Out_File << "ZONE T =  \"Block Number = " << Block_Number
	   << "\" \\ \n"
	   << "I = " << SolnBlk.Grid.ICu - SolnBlk.Grid.ICl + 3 << " \\ \n"
	   << "J = " << SolnBlk.Grid.JCu - SolnBlk.Grid.JCl + 3 << " \\ \n"
	   << "F = POINT \n";

  for (int j = SolnBlk.JCl-1; j <= SolnBlk.JCu+1; j++) {
    for (int i = SolnBlk.ICl-1; i <= SolnBlk.ICu+1; i++) {
      Linear_Reconstruction_LeastSquares(SolnBlk,i,j,LIMITER_VENKATAKRISHNAN);
      Out_File.setf(ios::scientific);
      Out_File << " " << SolnBlk.Grid.Cell[i][j].Xc
	       << " " << SolnBlk.dWdx[i][j].rho
 	       << " " << SolnBlk.dWdy[i][j].rho
 	       << " " << SolnBlk.dWdx[i][j].v.x
 	       << " " << SolnBlk.dWdy[i][j].v.x
 	       << " " << SolnBlk.dWdx[i][j].v.y
 	       << " " << SolnBlk.dWdy[i][j].v.y
 	       << " " << SolnBlk.dWdx[i][j].p
 	       << " " << SolnBlk.dWdy[i][j].p;
      if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
	Out_File << " " << SolnBlk.dWdx[i][j].k
		 << " " << SolnBlk.dWdy[i][j].k
		 << " " << SolnBlk.dWdx[i][j].omega
		 << " " << SolnBlk.dWdy[i][j].omega;
      }
      Out_File << " " << SolnBlk.dWdx[i][j].rho + SolnBlk.dWdy[i][j].rho
	       << " " << SolnBlk.dWdx[i][j].v.x + SolnBlk.dWdy[i][j].v.y
	       << " " << SolnBlk.dWdx[i][j].v.y - SolnBlk.dWdy[i][j].v.x;
      Out_File << " " << SolnBlk.phi[i][j].rho
 	       << " " << SolnBlk.phi[i][j].v.x
 	       << " " << SolnBlk.phi[i][j].v.y
 	       << " " << SolnBlk.phi[i][j].p;
      if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
	Out_File << " " << SolnBlk.phi[i][j].k
		 << " " << SolnBlk.phi[i][j].omega;
      }
      Out_File << endl;
      Out_File.unsetf(ios::scientific);
    }
  }
  Out_File << setprecision(6);

}

/**********************************************************************
 * Routine: Output_Quasi3D_Tecplot                                    *
 *                                                                    *
 * Writes the solution values at the nodes of the specified active    *
 * quadrilateral solution block to the specified output stream        *
 * suitable for plotting with TECPLOT in a quasi-3D format.           *
 *                                                                    *
 **********************************************************************/
void Output_Quasi3D_Tecplot(NavierStokes2D_Quad_Block &SolnBlk,
			    NavierStokes2D_Input_Parameters &IP,
			    const int Number_of_Time_Steps,
			    const double &Time,
			    const int Block_Number,
			    const int Output_Title,
			    ostream &Out_File) {

//   NavierStokes2D_pState W_node;
//   int numberofcells = (SolnBlk.JCu-SolnBlk.JCl+1)*(SolnBlk.ICu-SolnBlk.ICl+1);
//   int nrr, numberofrotations = 360/15;
//   int numberofnodes = ((SolnBlk.Grid.INu - SolnBlk.Grid.INl + 1)*
// 		       (SolnBlk.Grid.JNu - SolnBlk.Grid.JNl + 1));

//   // Ensure boundary conditions are updated before evaluating
//   // solution at the nodes.
//   BCs(SolnBlk,IP);

//   // Output node solution data.  
//   Out_File << setprecision(14);
//   if (Output_Title) {
//     Out_File << "TITLE = \"" << CFFC_Name() << ": 2D NavierStokes Solution, "
// 	     << "Time Step/Iteration Level = " << Number_of_Time_Steps
// 	     << ", Time = " << Time
// 	     << "\"" << "\n"
// 	     << "VARIABLES = \"x\" \\ \n"
// 	     << "\"y\" \\ \n"
// 	     << "\"z\" \\ \n"
// 	     << "\"rho\" \\ \n"
// 	     << "\"u\" \\ \n"
// 	     << "\"v\" \\ \n"
// 	     << "\"w\" \\ \n"
// 	     << "\"p\" \\ \n"
// 	     << "\"T\" \\ \n"
// 	     << "\"M\" \\ \n"
// 	     << "\"H\" \\ \n"
// 	     << "\"s\" \\ \n";
//     if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
//       Out_File << "\"k\" \\ \n"
// 	       << "\"omega\" \\ \n"
// 	       << "\"epsilon\" \\ \n"
// 	       << "\"ell\" \\ \n"
// 	       << "\"p_modified\" \\ \n";
//     }
//   }

//   for (int nr = 0; nr < numberofrotations; nr++) {
//     Out_File << "ZONE T =  \"Block Number = " << Block_Number << nr
// 	     << "\" \\ \n"
// 	     << "N = " << 2*numberofnodes << " \\ \n"
// 	     << "E = " << numberofcells << " \\ \n"
// 	     << "F = FEPOINT \\ \n"
// 	     << "ET = BRICK \n";
//     for (int j = SolnBlk.Grid.JNl; j <= SolnBlk.Grid.JNu; j++) {
//       for (int i = SolnBlk.Grid.INl; i <= SolnBlk.Grid.INu; i++) {
// 	W_node = SolnBlk.Wn(i,j);
// 	Out_File.setf(ios::scientific);
// 	Out_File << " " << SolnBlk.Grid.Node[i][j].X.x
// 		 << " " << SolnBlk.Grid.Node[i][j].X.y*sin(TWO*PI*double(nr)/double(numberofrotations))
// 		 << " " << SolnBlk.Grid.Node[i][j].X.y*cos(TWO*PI*double(nr)/double(numberofrotations))
// 		 << " " << W_node.rho
// 		 << " " << W_node.v.x
// 		 << " " << W_node.v.y*sin(TWO*PI*double(nr)/double(numberofrotations))
// 		 << " " << W_node.v.y*cos(TWO*PI*double(nr)/double(numberofrotations))
// 		 << " " << W_node.p
// 		 << " " << W_node.T()
// 		 << " " << W_node.v.abs()/W_node.a() 
// 		 << " " << W_node.H()
// 		 << " " << W_node.s();
// 	if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
// 	  Out_File << " " << W_node.k
// 		   << " " << W_node.omega
// 		   << " " << W_node.epsilon()
// 		   << " " << W_node.ell()
// 		   << " " << W_node.pmodified();
// 	}
// 	Out_File << endl;
//       }
//     }
//     if (nr < numberofrotations-1) nrr = nr + 1;
//     else nrr = 0;
//     for (int j = SolnBlk.Grid.JNl; j <= SolnBlk.Grid.JNu; j++) {
//       for (int i = SolnBlk.Grid.INl; i <= SolnBlk.Grid.INu; i++) {
// 	W_node = SolnBlk.Wn(i,j);
// 	Out_File.setf(ios::scientific);
// 	Out_File << " " << SolnBlk.Grid.Node[i][j].X.x
// 		 << " " << SolnBlk.Grid.Node[i][j].X.y*sin(TWO*PI*double(nrr)/double(numberofrotations))
// 		 << " " << SolnBlk.Grid.Node[i][j].X.y*cos(TWO*PI*double(nrr)/double(numberofrotations))
// 		 << " " << W_node.rho
// 		 << " " << W_node.v.x
// 		 << " " << W_node.v.y*sin(TWO*PI*double(nrr)/double(numberofrotations))
// 		 << " " << W_node.v.y*cos(TWO*PI*double(nrr)/double(numberofrotations))
// 		 << " " << W_node.p
// 		 << " " << W_node.T()
// 		 << " " << W_node.v.abs()/W_node.a()
// 		 << " " << W_node.H()
// 		 << " " << W_node.s();
// 	if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
// 	  Out_File << " " << W_node.k
// 		   << " " << W_node.omega
// 		   << " " << W_node.epsilon()
// 		   << " " << W_node.ell()
// 		   << " " << W_node.pmodified();
// 	}
// 	Out_File << endl;
//       }
//     }

//     // Connectivity table.
//     for (int j = SolnBlk.Grid.JCl; j <= SolnBlk.Grid.JCu; j++) {
//       for (int i = SolnBlk.Grid.ICl; i <= SolnBlk.Grid.ICu; i++) {
// 	Out_File << (j-2)*(SolnBlk.Grid.INu-1) + i - 1                         << " "
// 		 << (j-2)*(SolnBlk.Grid.INu-1) + i                             << " "
// 		 << (j-1)*(SolnBlk.Grid.INu-2) + i + 1 + (j-2)                 << " "
// 		 << (j-1)*(SolnBlk.Grid.INu-2) + i     + (j-2)                 << " "
// 		 << (j-2)*(SolnBlk.Grid.INu-1) + i - 1         + numberofnodes << " "
// 		 << (j-2)*(SolnBlk.Grid.INu-1) + i             + numberofnodes << " "
// 		 << (j-1)*(SolnBlk.Grid.INu-2) + i + 1 + (j-2) + numberofnodes << " "
// 		 << (j-1)*(SolnBlk.Grid.INu-2) + i     + (j-2) + numberofnodes << endl;
//       }
//     }
//   }

//   Out_File << setprecision(6);

  NavierStokes2D_pState W_node;
  int numberofcells = (SolnBlk.JCu-SolnBlk.JCl+1)*(SolnBlk.ICu-SolnBlk.ICl+1);
  int nrr, numberofrotations = 360/15;
  int numberofnodes = ((SolnBlk.Grid.INu - SolnBlk.Grid.INl + 1)*
		       (SolnBlk.Grid.JNu - SolnBlk.Grid.JNl + 1));

  // Ensure boundary conditions are updated before evaluating
  // solution at the nodes.
  BCs(SolnBlk,IP);

  // Output node solution data.  
  Out_File << setprecision(14);
  if (Output_Title) {
    Out_File << "TITLE = \"" << CFFC_Name() << ": 2D NavierStokes Solution, "
	     << "Time Step/Iteration Level = " << Number_of_Time_Steps
	     << ", Time = " << Time
	     << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
	     << "\"y\" \\ \n"
	     << "\"z\" \\ \n"
	     << "\"rho\" \\ \n"
	     << "\"u\" \\ \n"
	     << "\"v\" \\ \n"
	     << "\"w\" \\ \n"
	     << "\"p\" \\ \n"
	     << "\"T\" \\ \n"
	     << "\"M\" \\ \n";
    if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
      Out_File << "\"k\" \\ \n"
	       << "\"omega\" \\ \n";
    }
  }

  Out_File << "ZONE T =  \"Block Number = " << Block_Number
	   << "\" \\ \n"
	   << "I = " << SolnBlk.Grid.INu - SolnBlk.Grid.INl + 1 << " \\ \n"
	   << "J = " << SolnBlk.Grid.JNu - SolnBlk.Grid.JNl + 1 << " \\ \n"
	   << "K = " << numberofrotations
	   << "DATAPACKIING = POINT \\ \n";
  for (int j = SolnBlk.Grid.JNl; j <= SolnBlk.Grid.JNu; j++) {
    for (int i = SolnBlk.Grid.INl; i <= SolnBlk.Grid.INu; i++) {
      W_node = SolnBlk.Wn(i,j);
      for (int nr = 0; nr < 3; nr++) {
	Out_File.setf(ios::scientific);
	Out_File << " " << SolnBlk.Grid.Node[i][j].X.x
		 << " " << SolnBlk.Grid.Node[i][j].X.y*sin(TWO*PI*double(nr)/double(numberofrotations))
		 << " " << SolnBlk.Grid.Node[i][j].X.y*cos(TWO*PI*double(nr)/double(numberofrotations))
		 << " " << W_node.rho
		 << " " << W_node.v.x
		 << " " << W_node.v.y*sin(TWO*PI*double(nr)/double(numberofrotations))
		 << " " << W_node.v.y*cos(TWO*PI*double(nr)/double(numberofrotations))
		 << " " << W_node.p
		 << " " << W_node.T()
	       << " " << W_node.M();
	if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
	  Out_File << " " << W_node.k
		   << " " << W_node.omega;
	}
	Out_File << endl;
      }
    }
  }

  Out_File << setprecision(6);

}

/**********************************************************************
 * Routine: Output_Ringleb_Tecplot                                    *
 *                                                                    *
 * Writes the exact and computed Ringleb's flow solution values at    *
 * the nodes of the specified quadrilateral solution block to the     *
 * specified output stream suitable for plotting with TECPLOT.  The   *
 * error norms are also computed.                                     *
 *                                                                    *
 **********************************************************************/
void Output_Ringleb_Tecplot(NavierStokes2D_Quad_Block &SolnBlk,
			    const int Block_Number,
			    const int Output_Title,
			    ostream &Out_File,
			    double &l1_norm,
			    double &l2_norm,
			    double &max_norm,
			    double &area) {

  NavierStokes2D_pState We;

  // Output node solution data.  
  Out_File << setprecision(14);
  if (Output_Title)
    Out_File << "TITLE = \"" << CFFC_Name() << ": 2D NavierStokes Ringleb Flow Solution "
	     << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
	     << "\"y\" \\ \n"
	     << "\"rho\" \\ \n"
	     << "\"a\" \\ \n"
	     << "\"rho_exact\" \\ \n"
	     << "\"a_exact\" \\ \n"
	     << "\"rho_error\" \\ \n"
	     << "\"a_error\" \\ \n";
  Out_File << "ZONE T =  \"Block Number = " << Block_Number << "\" \\ \n"
	   << "I = " << SolnBlk.Grid.ICu - SolnBlk.Grid.ICl + 1 << " \\ \n"
	   << "J = " << SolnBlk.Grid.JCu - SolnBlk.Grid.JCl + 1 << " \\ \n"
	   << "F = POINT \\ \n";
  Out_File.setf(ios::scientific);
  for (int j = SolnBlk.Grid.JCl; j <= SolnBlk.Grid.JCu; j++) {
    for (int i = SolnBlk.Grid.ICl; i <= SolnBlk.Grid.ICu; i++) {
      We = RinglebFlowAverageState(We,SolnBlk.Grid.nodeSW(i,j).X,
				   SolnBlk.Grid.nodeSE(i,j).X,
				   SolnBlk.Grid.nodeNE(i,j).X,
				   SolnBlk.Grid.nodeNW(i,j).X);
      Out_File << SolnBlk.Grid.Cell[i][j].Xc << " " 
	       << SolnBlk.W[i][j].rho << " " 
	       << SolnBlk.W[i][j].a() << " " 
	       << We.rho << " " 
	       << We.a()   << " " 
	       << fabs(SolnBlk.W[i][j].rho - We.rho) << " "
	       << fabs(SolnBlk.W[i][j].a() - We.a()) << endl;
    }
  }
  Out_File.unsetf(ios::scientific);

  // Determine the error norms.
  for (int j = SolnBlk.Grid.JCl; j <= SolnBlk.Grid.JCu; j++) {
    for (int i = SolnBlk.Grid.ICl; i <= SolnBlk.Grid.ICu; i++) {
      We = RinglebFlowAverageState(We,SolnBlk.Grid.nodeSW(i,j).X,
				   SolnBlk.Grid.nodeSE(i,j).X,
				   SolnBlk.Grid.nodeNE(i,j).X,
				   SolnBlk.Grid.nodeNW(i,j).X);
      l1_norm += fabs(We.rho - SolnBlk.W[i][j].rho)*SolnBlk.Grid.Cell[i][j].A;
      l2_norm += sqr(We.rho - SolnBlk.W[i][j].rho)*SolnBlk.Grid.Cell[i][j].A;
      max_norm = max(max_norm,fabs(We.rho - SolnBlk.W[i][j].rho));
      area += SolnBlk.Grid.Cell[i][j].A;
    }
  }

}

/**********************************************************************
 * Routine: Output_Viscous_Channel_Tecplot                            *
 *                                                                    *
 * Writes the exact and computed laminar channel flow solution        *
 * (Couette/Poiseuille solutions) at the nodes of the specified       *
 * quadrilateral solution block to the specified output stream        *
 * suitable for plotting with TECPLOT.  The error norms of the        *
 * u-velocity component are also computed (L1, L2, and max norms).    *
 *                                                                    *
 **********************************************************************/
void Output_Viscous_Channel_Tecplot(NavierStokes2D_Quad_Block &SolnBlk,
				    NavierStokes2D_Input_Parameters &IP,
				    const int Block_Number,
				    const int Output_Title,
				    ostream &Out_File,
				    double &l1_norm,
				    double &l2_norm,
				    double &max_norm) {

  NavierStokes2D_pState We, W;
  Vector2D dX;

  // Output node solution data.  
  Out_File << setprecision(14);
  if (Output_Title)
    Out_File << "TITLE = \"" << CFFC_Name() << ": 2D NavierStokes Viscous Channel Flow Solution "
	     << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
	     << "\"y\" \\ \n"
	     << "\"rho\" \\ \n"
	     << "\"u\" \\ \n"
	     << "\"v\" \\ \n"
	     << "\"p\" \\ \n"
	     << "\"T\" \\ \n"
	     << "\"rho_exact\" \\ \n"
	     << "\"u_exact\" \\ \n"
	     << "\"v_exact\" \\ \n"
	     << "\"p_exact\" \\ \n"
	     << "\"rho_error\" \\ \n"
	     << "\"u_error\" \\ \n"
	     << "\"v_error\" \\ \n"
	     << "\"p_error\" \\ \n";
  Out_File << "ZONE T =  \"Block Number = " << Block_Number
	   << "\" \\ \n"
	   << "I = " << SolnBlk.Grid.NCi << " \\ \n"
	   << "J = " << SolnBlk.Grid.NCj << " \\ \n"
	   << "F = POINT \\ \n";

  Out_File.setf(ios::scientific);
  for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
    for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
      // Get the computed solution.
      W = SolnBlk.W[i][j];
      // Determine the exact solution.  Apply the appropriate boundary
      // condition to the exact solution as necessary.
      if (j == SolnBlk.JCl-1 && i >= SolnBlk.ICl && i <= SolnBlk.ICu) {
	// Apply the adiabtic viscous wall boundary condition to the
	// exact solution at cell (i,JCl).
	We = ViscousChannelFlow(IP.Wo,SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc,
				SolnBlk.Vwall,IP.dp,IP.Box_Width,IP.Box_Height);
	We = WallViscousHeatFlux(We,SolnBlk.Grid.nfaceS(i,SolnBlk.JCl));
      } else if (j == SolnBlk.JCl-2 && i >= SolnBlk.ICl && i <= SolnBlk.ICu) {
	// Apply the adiabtic viscous wall boundary condition to the
	// exact solution at cell (i,JCl+1).
	We = ViscousChannelFlow(IP.Wo,SolnBlk.Grid.Cell[i][SolnBlk.JCl+1].Xc,
				SolnBlk.Vwall,IP.dp,IP.Box_Width,IP.Box_Height);
	We = WallViscousHeatFlux(We,SolnBlk.Grid.nfaceS(i,SolnBlk.JCl));
      } else if (j == SolnBlk.JCu+1 && i >= SolnBlk.ICl && i <= SolnBlk.ICu) {
	// Apply the adiabtic viscous wall boundary condition to the
	// exact solution at cell (i,JCu).
 	We = ViscousChannelFlow(IP.Wo,SolnBlk.Grid.Cell[i][SolnBlk.JCu].Xc,
 				SolnBlk.Vwall,IP.dp,IP.Box_Width,IP.Box_Height);
 	We = MovingWallHeatFlux(We,SolnBlk.Grid.nfaceN(i,SolnBlk.JCu),SolnBlk.Vwall.x);
      } else if (j == SolnBlk.JCu+2 && i >= SolnBlk.ICl && i <= SolnBlk.ICu) {
	// Apply the adiabtic viscous wall boundary condition to the
	// exact solution at cell (i,JCu-1).
 	We = ViscousChannelFlow(IP.Wo,SolnBlk.Grid.Cell[i][SolnBlk.JCu-1].Xc,
 				SolnBlk.Vwall,IP.dp,IP.Box_Width,IP.Box_Height);
 	We = MovingWallHeatFlux(We,SolnBlk.Grid.nfaceN(i,SolnBlk.JCu),SolnBlk.Vwall.x);
      } else if (i == SolnBlk.ICl-1 && j >= SolnBlk.JCl && j <= SolnBlk.JCu) {
	// Apply the subsonic inflow boundary condition to the exact
	// solution at cell (ICl,j).
 	We = ViscousChannelFlow(IP.Wo,SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc,
 				SolnBlk.Vwall,IP.dp,IP.Box_Width,IP.Box_Height);
	dX = SolnBlk.Grid.Cell[SolnBlk.ICl-1][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc;
	We.p = SolnBlk.W[SolnBlk.ICl][j].p + ((SolnBlk.WoE[j].p - SolnBlk.WoW[j].p)/
					      fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc.x - 
						   SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc.x))*dX.x;
	We = ViscousChannelFlowVelocity(We,SolnBlk.Grid.Cell[i][j].Xc,
					SolnBlk.Vwall,IP.dp,IP.Box_Width,IP.Box_Height);
      } else if (i == SolnBlk.ICl-2 && j >= SolnBlk.JCl && j <= SolnBlk.JCu) {
	// Apply the subsonic inflow boundary condition to the exact
	// solution at cell (ICl,j).
 	We = ViscousChannelFlow(IP.Wo,SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc,
 				SolnBlk.Vwall,IP.dp,IP.Box_Width,IP.Box_Height);
	dX = SolnBlk.Grid.Cell[SolnBlk.ICl-2][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc;
	We.p = SolnBlk.W[SolnBlk.ICl][j].p + ((SolnBlk.WoE[j].p - SolnBlk.WoW[j].p)/
					      fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc.x - 
						   SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc.x))*dX.x;
	We = ViscousChannelFlowVelocity(We,SolnBlk.Grid.Cell[i][j].Xc,
					SolnBlk.Vwall,IP.dp,IP.Box_Width,IP.Box_Height);
      } else if (i == SolnBlk.ICu+1 && j >= SolnBlk.JCl && j <= SolnBlk.JCu) {
	// Apply the subsonic outflow boundary condition to the exact
	// solution at cell (ICu,j).
 	We = ViscousChannelFlow(IP.Wo,SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc,
 				SolnBlk.Vwall,IP.dp,IP.Box_Width,IP.Box_Height);
	dX = SolnBlk.Grid.Cell[SolnBlk.ICu+1][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc;
	We.p = SolnBlk.W[SolnBlk.ICu][j].p + ((SolnBlk.WoE[j].p - SolnBlk.WoW[j].p)/
					      fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc.x - 
						   SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc.x))*dX.x;
	We = ViscousChannelFlowVelocity(We,SolnBlk.Grid.Cell[i][j].Xc,
					SolnBlk.Vwall,IP.dp,IP.Box_Width,IP.Box_Height);
      } else if (i == SolnBlk.ICu+2 && j >= SolnBlk.JCl && j <= SolnBlk.JCu) {
	// Apply the subsonic outflow boundary condition to the exact
	// solution at cell (ICu,j).
 	We = ViscousChannelFlow(IP.Wo,SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc,
 				SolnBlk.Vwall,IP.dp,IP.Box_Width,IP.Box_Height);
	dX = SolnBlk.Grid.Cell[SolnBlk.ICu+2][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc;
	We.p = SolnBlk.W[SolnBlk.ICu][j].p + ((SolnBlk.WoE[j].p - SolnBlk.WoW[j].p)/
					      fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc.x - 
						   SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc.x))*dX.x;
	We = ViscousChannelFlowVelocity(We,SolnBlk.Grid.Cell[i][j].Xc,
					SolnBlk.Vwall,IP.dp,IP.Box_Width,IP.Box_Height);
      } else if (i < SolnBlk.ICl && j == SolnBlk.JCl-1) {
	// Apply the subsonic inflow boundary condition to the exact
	// solution at cell (ICl,JCl) then apply the adiabatic viscous
	// wall boundary condition.
 	We = ViscousChannelFlow(IP.Wo,SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCl].Xc,
 				SolnBlk.Vwall,IP.dp,IP.Box_Width,IP.Box_Height);
	dX = SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc - SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCl].Xc;
	We.p = SolnBlk.W[SolnBlk.ICl][SolnBlk.JCl].p + ((SolnBlk.WoE[SolnBlk.JCl].p - SolnBlk.WoW[SolnBlk.JCl].p)/
							fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCl].Xc.x - 
							     SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCl].Xc.x))*dX.x;
	We = ViscousChannelFlowVelocity(We,SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc,
					SolnBlk.Vwall,IP.dp,IP.Box_Width,IP.Box_Height);
 	We = WallViscousHeatFlux(We,SolnBlk.Grid.nfaceS(i,SolnBlk.JCl));
      } else if (i < SolnBlk.ICl && j == SolnBlk.JCl-2) {
	// Apply the subsonic inflow boundary condition to the exact
	// solution at cell (ICl,JCl+1) then apply the adiabatic viscous
	// wall boundary condition.
 	We = ViscousChannelFlow(IP.Wo,SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCl+1].Xc,
 				SolnBlk.Vwall,IP.dp,IP.Box_Width,IP.Box_Height);
	dX = SolnBlk.Grid.Cell[i][SolnBlk.JCl+1].Xc - SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCl+1].Xc;
	We.p = SolnBlk.W[SolnBlk.ICl][SolnBlk.JCl+1].p + ((SolnBlk.WoE[SolnBlk.JCl+1].p - SolnBlk.WoW[SolnBlk.JCl+1].p)/
							  fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCl+1].Xc.x - 
							       SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCl+1].Xc.x))*dX.x;
	We = ViscousChannelFlowVelocity(We,SolnBlk.Grid.Cell[i][SolnBlk.JCl+1].Xc,
					SolnBlk.Vwall,IP.dp,IP.Box_Width,IP.Box_Height);
 	We = WallViscousHeatFlux(We,SolnBlk.Grid.nfaceS(i,SolnBlk.JCl));
      } else if (i > SolnBlk.ICu && j == SolnBlk.JCl-1) {
	// Apply the subsonic inflow boundary condition to the exact
	// solution at cell (ICu,JCl) then apply the adiabatic viscous
	// wall boundary condition.
 	We = ViscousChannelFlow(IP.Wo,SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCl].Xc,
 				SolnBlk.Vwall,IP.dp,IP.Box_Width,IP.Box_Height);
	dX = SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc - SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCl].Xc;
	We.p = SolnBlk.W[SolnBlk.ICu][SolnBlk.JCl].p + ((SolnBlk.WoE[SolnBlk.JCl].p - SolnBlk.WoW[SolnBlk.JCl].p)/
							fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCl].Xc.x - 
							     SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCl].Xc.x))*dX.x;
	We = ViscousChannelFlowVelocity(We,SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc,
					SolnBlk.Vwall,IP.dp,IP.Box_Width,IP.Box_Height);
 	We = WallViscousHeatFlux(We,SolnBlk.Grid.nfaceS(i,SolnBlk.JCl));
      } else if (i > SolnBlk.ICu && j == SolnBlk.JCl-2) {
	// Apply the subsonic inflow boundary condition to the exact
	// solution at cell (ICu,JCl+1) then apply the adiabatic viscous
	// wall boundary condition.
 	We = ViscousChannelFlow(IP.Wo,SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCl+1].Xc,
 				SolnBlk.Vwall,IP.dp,IP.Box_Width,IP.Box_Height);
	dX = SolnBlk.Grid.Cell[i][SolnBlk.JCl+1].Xc - SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCl+1].Xc;
	We.p = SolnBlk.W[SolnBlk.ICu][SolnBlk.JCl+1].p + ((SolnBlk.WoE[SolnBlk.JCl+1].p - SolnBlk.WoW[SolnBlk.JCl+1].p)/
							  fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCl+1].Xc.x - 
							       SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCl+1].Xc.x))*dX.x;
	We = ViscousChannelFlowVelocity(We,SolnBlk.Grid.Cell[i][SolnBlk.JCl+1].Xc,
					SolnBlk.Vwall,IP.dp,IP.Box_Width,IP.Box_Height);
 	We = WallViscousHeatFlux(We,SolnBlk.Grid.nfaceS(i,SolnBlk.JCl));
      } else if (i < SolnBlk.ICl && j == SolnBlk.JCu+1) {
	// Apply the subsonic inflow boundary condition to the exact
	// solution at cell (ICl,JCu) then apply the adiabatic viscous
	// wall boundary condition.
 	We = ViscousChannelFlow(IP.Wo,SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCu].Xc,
 				SolnBlk.Vwall,IP.dp,IP.Box_Width,IP.Box_Height);
	dX = SolnBlk.Grid.Cell[i][SolnBlk.JCu].Xc - SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCu].Xc;
	We.p = SolnBlk.W[SolnBlk.ICl][SolnBlk.JCu].p + ((SolnBlk.WoE[SolnBlk.JCu].p - SolnBlk.WoW[SolnBlk.JCu].p)/
							fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCu].Xc.x - 
							     SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCu].Xc.x))*dX.x;
	We = ViscousChannelFlowVelocity(We,SolnBlk.Grid.Cell[i][SolnBlk.JCu].Xc,
					SolnBlk.Vwall,IP.dp,IP.Box_Width,IP.Box_Height);
 	We = MovingWallHeatFlux(We,SolnBlk.Grid.nfaceS(i,SolnBlk.JCu),SolnBlk.Vwall.x);
      } else if (i < SolnBlk.ICl && j == SolnBlk.JCu+2) {
	// Apply the subsonic inflow boundary condition to the exact
	// solution at cell (ICl,JCu-1) then apply the adiabatic viscous
	// wall boundary condition.
 	We = ViscousChannelFlow(IP.Wo,SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCu-1].Xc,
 				SolnBlk.Vwall,IP.dp,IP.Box_Width,IP.Box_Height);
	dX = SolnBlk.Grid.Cell[i][SolnBlk.JCu-1].Xc - SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCu-1].Xc;
	We.p = SolnBlk.W[SolnBlk.ICl][SolnBlk.JCu-1].p + ((SolnBlk.WoE[SolnBlk.JCu-1].p - SolnBlk.WoW[SolnBlk.JCu-1].p)/
							  fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCu-1].Xc.x - 
							       SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCu-1].Xc.x))*dX.x;
	We = ViscousChannelFlowVelocity(We,SolnBlk.Grid.Cell[i][SolnBlk.JCu-1].Xc,
					SolnBlk.Vwall,IP.dp,IP.Box_Width,IP.Box_Height);
 	We = MovingWallHeatFlux(We,SolnBlk.Grid.nfaceS(i,SolnBlk.JCu),SolnBlk.Vwall.x);
      } else if (i > SolnBlk.ICu && j == SolnBlk.JCu+1) {
	// Apply the subsonic inflow boundary condition to the exact
	// solution at cell (ICu,JCu) then apply the adiabatic viscous
	// wall boundary condition.
 	We = ViscousChannelFlow(IP.Wo,SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCu].Xc,
 				SolnBlk.Vwall,IP.dp,IP.Box_Width,IP.Box_Height);
	dX = SolnBlk.Grid.Cell[i][SolnBlk.JCu].Xc - SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCu].Xc;
	We.p = SolnBlk.W[SolnBlk.ICu][SolnBlk.JCu].p + ((SolnBlk.WoE[SolnBlk.JCu].p - SolnBlk.WoW[SolnBlk.JCu].p)/
							fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCu].Xc.x - 
							     SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCu].Xc.x))*dX.x;
	We = ViscousChannelFlowVelocity(We,SolnBlk.Grid.Cell[i][SolnBlk.JCu].Xc,
					SolnBlk.Vwall,IP.dp,IP.Box_Width,IP.Box_Height);
 	We = MovingWallHeatFlux(We,SolnBlk.Grid.nfaceS(i,SolnBlk.JCu),SolnBlk.Vwall.x);
      } else if (i > SolnBlk.ICu && j == SolnBlk.JCu+2) {
	// Apply the subsonic inflow boundary condition to the exact
	// solution at cell (ICu,JCu-1) then apply the adiabatic viscous
	// wall boundary condition.
 	We = ViscousChannelFlow(IP.Wo,SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCu-1].Xc,
 				SolnBlk.Vwall,IP.dp,IP.Box_Width,IP.Box_Height);
	dX = SolnBlk.Grid.Cell[i][SolnBlk.JCu-1].Xc - SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCu-1].Xc;
	We.p = SolnBlk.W[SolnBlk.ICu][SolnBlk.JCu-1].p + ((SolnBlk.WoE[SolnBlk.JCu-1].p - SolnBlk.WoW[SolnBlk.JCu-1].p)/
							  fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCu-1].Xc.x - 
							       SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCu-1].Xc.x))*dX.x;
	We = ViscousChannelFlowVelocity(We,SolnBlk.Grid.Cell[i][SolnBlk.JCu-1].Xc,
					SolnBlk.Vwall,IP.dp,IP.Box_Width,IP.Box_Height);
 	We = MovingWallHeatFlux(We,SolnBlk.Grid.nfaceS(i,SolnBlk.JCu),SolnBlk.Vwall.x);
      } else {
	// Internal point.  Determine the exact solution and compute error norms.
 	We = ViscousChannelFlow(IP.Wo,SolnBlk.Grid.Cell[i][j].Xc,
 				SolnBlk.Vwall,IP.dp,IP.Box_Width,IP.Box_Height);
	l1_norm += fabs(W.v.x - We.v.x);
	l2_norm += sqr(W.v.x - We.v.x);
	max_norm = max(max_norm,fabs(W.v.x - We.v.x));
      }
      //if (i < SolnBlk.ICl || i > SolnBlk.ICu || j < SolnBlk.JCl || j > SolnBlk.JCu) We = SolnBlk.W[i][j];
      // Output the data.
      Out_File << SolnBlk.Grid.Cell[i][j].Xc << " " 
	       << W.rho << " " 
	       << W.v.x << " " 
	       << W.v.y << " " 
	       << W.p << " " 
	       << W.T() << " " 
	       << We.rho << " " 
	       << We.v.x << " " 
	       << We.v.y << " " 
	       << We.p << " " 
	       << fabs(W.rho - We.rho) << " "
	       << fabs(W.v.x - We.v.x) << " "
	       << fabs(W.v.y - We.v.y) << " "
	       << fabs(W.p - We.p) << " "
	       << endl;
    }
  }
  l2_norm = sqrt(l2_norm);

}

/**********************************************************************
 * Routine: Output_Viscous_Pipe_Tecplot                               *
 *                                                                    *
 * Writes the exact and computed laminar pipe flow solution (Hagen-   *
 * Poiseuille solution) at the nodes of the specified quadrilateral   *
 * solution block to the specified output stream suitable for         *
 * plotting with TECPLOT.  The error norms of the u-velocity          *
 * component are also computed (L1, L2, and max norms).               *
 *                                                                    *
 **********************************************************************/
void Output_Viscous_Pipe_Tecplot(NavierStokes2D_Quad_Block &SolnBlk,
				 NavierStokes2D_Input_Parameters &IP,
				 const int Block_Number,
				 const int Output_Title,
				 ostream &Out_File,
				 double &l1_norm,
				 double &l2_norm,
				 double &max_norm) {

  NavierStokes2D_pState We, W;
  Vector2D dX;

  // Output node solution data.  
  Out_File << setprecision(14);
  if (Output_Title)
    Out_File << "TITLE = \"" << CFFC_Name() << ": 2D NavierStokes Viscous Pipe Flow Solution "
	     << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
	     << "\"y\" \\ \n"
	     << "\"rho\" \\ \n"
	     << "\"u\" \\ \n"
	     << "\"v\" \\ \n"
	     << "\"p\" \\ \n"
	     << "\"T\" \\ \n"
	     << "\"rho_exact\" \\ \n"
	     << "\"u_exact\" \\ \n"
	     << "\"v_exact\" \\ \n"
	     << "\"p_exact\" \\ \n"
	     << "\"rho_error\" \\ \n"
	     << "\"u_error\" \\ \n"
	     << "\"v_error\" \\ \n"
	     << "\"p_error\" \\ \n";
  Out_File << "ZONE T =  \"Block Number = " << Block_Number
	   << "\" \\ \n"
	   << "I = " << SolnBlk.Grid.NCi << " \\ \n"
	   << "J = " << SolnBlk.Grid.NCj << " \\ \n"
	   << "F = POINT \\ \n";

  Out_File.setf(ios::scientific);
  for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
    for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
      // Get the computed solution.
      W = SolnBlk.W[i][j];
      // Determine the exact solution.  Apply the appropriate boundary
      // condition to the exact solution as necessary.
      if (j == SolnBlk.JCl-1 && i >= SolnBlk.ICl && i <= SolnBlk.ICu) {
	// Apply the reflection boundary condition to the exact solution
	// at cell (i,JCl).
	We = ViscousPipeFlow(IP.Wo,SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc,
			     IP.dp,IP.Pipe_Length,IP.Pipe_Radius);
	We = Reflect(We,SolnBlk.Grid.nfaceS(i,SolnBlk.JCl));
      } else if (j == SolnBlk.JCl-2 && i >= SolnBlk.ICl && i <= SolnBlk.ICu) {
	// Apply the reflection boundary condition to the exact solution
	// at cell (i,JCl+1).
	We = ViscousPipeFlow(IP.Wo,SolnBlk.Grid.Cell[i][SolnBlk.JCl+1].Xc,
			     IP.dp,IP.Pipe_Length,IP.Pipe_Radius);
	We = Reflect(We,SolnBlk.Grid.nfaceS(i,SolnBlk.JCl));
      } else if (j == SolnBlk.JCu+1 && i >= SolnBlk.ICl && i <= SolnBlk.ICu) {
	// Apply the adiabtic viscous wall boundary condition to the
	// exact solution at cell (i,JCu).
 	We = ViscousPipeFlow(IP.Wo,SolnBlk.Grid.Cell[i][SolnBlk.JCu].Xc,
			     IP.dp,IP.Pipe_Length,IP.Pipe_Radius);
 	We = WallViscousHeatFlux(We,SolnBlk.Grid.nfaceN(i,SolnBlk.JCu));
      } else if (j == SolnBlk.JCu+2 && i >= SolnBlk.ICl && i <= SolnBlk.ICu) {
	// Apply the adiabtic viscous wall boundary condition to the
	// exact solution at cell (i,JCu-1).
 	We = ViscousPipeFlow(IP.Wo,SolnBlk.Grid.Cell[i][SolnBlk.JCu-1].Xc,
			     IP.dp,IP.Pipe_Length,IP.Pipe_Radius);
 	We = WallViscousHeatFlux(We,SolnBlk.Grid.nfaceN(i,SolnBlk.JCu));
      } else if (i == SolnBlk.ICl-1 && j >= SolnBlk.JCl && j <= SolnBlk.JCu) {
	// Apply the subsonic inflow boundary condition to the exact
	// solution at cell (ICl,j).
 	We = ViscousPipeFlow(IP.Wo,SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc,
			     IP.dp,IP.Pipe_Length,IP.Pipe_Radius);
	dX = SolnBlk.Grid.Cell[SolnBlk.ICl-1][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc;
	We.p = SolnBlk.W[SolnBlk.ICl][j].p + ((SolnBlk.WoE[j].p - SolnBlk.WoW[j].p)/
					      fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc.x - 
						   SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc.x))*dX.x;
      } else if (i == SolnBlk.ICl-2 && j >= SolnBlk.JCl && j <= SolnBlk.JCu) {
	// Apply the subsonic inflow boundary condition to the exact
	// solution at cell (ICl,j).
 	We = ViscousPipeFlow(IP.Wo,SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc,
			     IP.dp,IP.Pipe_Length,IP.Pipe_Radius);
	dX = SolnBlk.Grid.Cell[SolnBlk.ICl-2][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc;
	We.p = SolnBlk.W[SolnBlk.ICl][j].p + ((SolnBlk.WoE[j].p - SolnBlk.WoW[j].p)/
					      fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc.x - 
						   SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc.x))*dX.x;
      } else if (i == SolnBlk.ICu+1 && j >= SolnBlk.JCl && j <= SolnBlk.JCu) {
	// Apply the subsonic outflow boundary condition to the exact
	// solution at cell (ICu,j).
 	We = ViscousPipeFlow(IP.Wo,SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc,
			     IP.dp,IP.Pipe_Length,IP.Pipe_Radius);
	dX = SolnBlk.Grid.Cell[SolnBlk.ICu+1][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc;
	We.p = SolnBlk.W[SolnBlk.ICu][j].p + ((SolnBlk.WoE[j].p - SolnBlk.WoW[j].p)/
					      fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc.x - 
						   SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc.x))*dX.x;
      } else if (i == SolnBlk.ICu+2 && j >= SolnBlk.JCl && j <= SolnBlk.JCu) {
	// Apply the subsonic outflow boundary condition to the exact
	// solution at cell (ICu,j).
 	We = ViscousPipeFlow(IP.Wo,SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc,
			     IP.dp,IP.Pipe_Length,IP.Pipe_Radius);
	dX = SolnBlk.Grid.Cell[SolnBlk.ICu+2][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc;
	We.p = SolnBlk.W[SolnBlk.ICu][j].p + ((SolnBlk.WoE[j].p - SolnBlk.WoW[j].p)/
					      fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc.x - 
						   SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc.x))*dX.x;
      } else if (i < SolnBlk.ICl && j == SolnBlk.JCl-1) {
	// Apply the subsonic inflow boundary condition to the exact
	// solution at cell (ICl,JCl) then apply the reflection
	// boundary condition.
 	We = ViscousPipeFlow(IP.Wo,SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCl].Xc,
			     IP.dp,IP.Pipe_Length,IP.Pipe_Radius);
	dX = SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc - SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCl].Xc;
	We.p = SolnBlk.W[SolnBlk.ICl][SolnBlk.JCl].p + ((SolnBlk.WoE[SolnBlk.JCl].p - SolnBlk.WoW[SolnBlk.JCl].p)/
							fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCl].Xc.x - 
							     SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCl].Xc.x))*dX.x;
 	We = Reflect(We,SolnBlk.Grid.nfaceS(i,SolnBlk.JCl));
      } else if (i < SolnBlk.ICl && j == SolnBlk.JCl-2) {
	// Apply the subsonic inflow boundary condition to the exact
	// solution at cell (ICl,JCl+1) then apply the reflection
	// boundary condition.
 	We = ViscousPipeFlow(IP.Wo,SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCl+1].Xc,
			     IP.dp,IP.Pipe_Length,IP.Pipe_Radius);
	dX = SolnBlk.Grid.Cell[i][SolnBlk.JCl+1].Xc - SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCl+1].Xc;
	We.p = SolnBlk.W[SolnBlk.ICl][SolnBlk.JCl+1].p + ((SolnBlk.WoE[SolnBlk.JCl+1].p - SolnBlk.WoW[SolnBlk.JCl+1].p)/
							  fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCl+1].Xc.x - 
							       SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCl+1].Xc.x))*dX.x;
 	We = Reflect(We,SolnBlk.Grid.nfaceS(i,SolnBlk.JCl));
      } else if (i > SolnBlk.ICu && j == SolnBlk.JCl-1) {
	// Apply the subsonic inflow boundary condition to the exact
	// solution at cell (ICu,JCl) then apply the adiabatic viscous
	// wall boundary condition.
 	We = ViscousPipeFlow(IP.Wo,SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCl].Xc,
			     IP.dp,IP.Pipe_Length,IP.Pipe_Radius);
	dX = SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc - SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCl].Xc;
	We.p = SolnBlk.W[SolnBlk.ICu][SolnBlk.JCl].p + ((SolnBlk.WoE[SolnBlk.JCl].p - SolnBlk.WoW[SolnBlk.JCl].p)/
							fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCl].Xc.x - 
							     SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCl].Xc.x))*dX.x;
 	We = Reflect(We,SolnBlk.Grid.nfaceS(i,SolnBlk.JCl));
      } else if (i > SolnBlk.ICu && j == SolnBlk.JCl-2) {
	// Apply the subsonic inflow boundary condition to the exact
	// solution at cell (ICu,JCl+1) then apply the adiabatic viscous
	// wall boundary condition.
 	We = ViscousPipeFlow(IP.Wo,SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCl+1].Xc,
			     IP.dp,IP.Pipe_Length,IP.Pipe_Radius);
	dX = SolnBlk.Grid.Cell[i][SolnBlk.JCl+1].Xc - SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCl+1].Xc;
	We.p = SolnBlk.W[SolnBlk.ICu][SolnBlk.JCl+1].p + ((SolnBlk.WoE[SolnBlk.JCl+1].p - SolnBlk.WoW[SolnBlk.JCl+1].p)/
							  fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCl+1].Xc.x - 
							       SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCl+1].Xc.x))*dX.x;
 	We = Reflect(We,SolnBlk.Grid.nfaceS(i,SolnBlk.JCl));
      } else if (i < SolnBlk.ICl && j == SolnBlk.JCu+1) {
	// Apply the subsonic inflow boundary condition to the exact
	// solution at cell (ICl,JCu) then apply the adiabatic viscous
	// wall boundary condition.
 	We = ViscousPipeFlow(IP.Wo,SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCu].Xc,
			     IP.dp,IP.Pipe_Length,IP.Pipe_Radius);
	dX = SolnBlk.Grid.Cell[i][SolnBlk.JCu].Xc - SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCu].Xc;
	We.p = SolnBlk.W[SolnBlk.ICl][SolnBlk.JCu].p + ((SolnBlk.WoE[SolnBlk.JCu].p - SolnBlk.WoW[SolnBlk.JCu].p)/
							fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCu].Xc.x - 
							     SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCu].Xc.x))*dX.x;
 	We = WallViscousHeatFlux(We,SolnBlk.Grid.nfaceS(i,SolnBlk.JCu));
      } else if (i < SolnBlk.ICl && j == SolnBlk.JCu+2) {
	// Apply the subsonic inflow boundary condition to the exact
	// solution at cell (ICl,JCu-1) then apply the adiabatic viscous
	// wall boundary condition.
 	We = ViscousPipeFlow(IP.Wo,SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCu-1].Xc,
			     IP.dp,IP.Pipe_Length,IP.Pipe_Radius);
	dX = SolnBlk.Grid.Cell[i][SolnBlk.JCu-1].Xc - SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCu-1].Xc;
	We.p = SolnBlk.W[SolnBlk.ICl][SolnBlk.JCu-1].p + ((SolnBlk.WoE[SolnBlk.JCu-1].p - SolnBlk.WoW[SolnBlk.JCu-1].p)/
							  fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCu-1].Xc.x - 
							       SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCu-1].Xc.x))*dX.x;
 	We = WallViscousHeatFlux(We,SolnBlk.Grid.nfaceS(i,SolnBlk.JCu));
      } else if (i > SolnBlk.ICu && j == SolnBlk.JCu+1) {
	// Apply the subsonic inflow boundary condition to the exact
	// solution at cell (ICu,JCu) then apply the adiabatic viscous
	// wall boundary condition.
 	We = ViscousPipeFlow(IP.Wo,SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCu].Xc,
			     IP.dp,IP.Pipe_Length,IP.Pipe_Radius);
	dX = SolnBlk.Grid.Cell[i][SolnBlk.JCu].Xc - SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCu].Xc;
	We.p = SolnBlk.W[SolnBlk.ICu][SolnBlk.JCu].p + ((SolnBlk.WoE[SolnBlk.JCu].p - SolnBlk.WoW[SolnBlk.JCu].p)/
							fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCu].Xc.x - 
							     SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCu].Xc.x))*dX.x;
 	We = WallViscousHeatFlux(We,SolnBlk.Grid.nfaceS(i,SolnBlk.JCu));
      } else if (i > SolnBlk.ICu && j == SolnBlk.JCu+2) {
	// Apply the subsonic inflow boundary condition to the exact
	// solution at cell (ICu,JCu-1) then apply the adiabatic viscous
	// wall boundary condition.
 	We = ViscousPipeFlow(IP.Wo,SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCu-1].Xc,
			     IP.dp,IP.Pipe_Length,IP.Pipe_Radius);
	dX = SolnBlk.Grid.Cell[i][SolnBlk.JCu-1].Xc - SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCu-1].Xc;
	We.p = SolnBlk.W[SolnBlk.ICu][SolnBlk.JCu-1].p + ((SolnBlk.WoE[SolnBlk.JCu-1].p - SolnBlk.WoW[SolnBlk.JCu-1].p)/
							  fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCu-1].Xc.x - 
							       SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCu-1].Xc.x))*dX.x;
 	We = WallViscousHeatFlux(We,SolnBlk.Grid.nfaceS(i,SolnBlk.JCu));
      } else {
	// Internal point.  Determine the exact solution and compute error norms.
 	We = ViscousPipeFlow(IP.Wo,SolnBlk.Grid.Cell[i][j].Xc,
 			     IP.dp,IP.Pipe_Length,IP.Pipe_Radius);
	l1_norm += fabs(W.v.x - We.v.x);
	l2_norm += sqr(W.v.x - We.v.x);
	max_norm = max(max_norm,fabs(W.v.x - We.v.x));
      }
      // Output the data.
      Out_File << SolnBlk.Grid.Cell[i][j].Xc << " "
	       << W.rho << " "
	       << W.v.x << " "
	       << W.v.y << " "
	       << W.p   << " "
	       << W.T() << " "
	       << We.rho << " "
	       << We.v.x << " "
	       << We.v.y << " "
	       << We.p << " "
	       << fabs(W.rho - We.rho) << " "
	       << fabs(W.v.x - We.v.x) << " "
	       << fabs(W.v.y - We.v.y) << " "
	       << fabs(W.p - We.p) << endl;
    }
  }
  l2_norm = sqrt(l2_norm);

}

/**********************************************************************
 * Routine: Output_Turbulent_Pipe_Tecplot                             *
 *                                                                    *
 * Writes the experiment and computed solution data for a turbulent   *
 * pipe flow for a 1D array of 2D quadrilateral multi-block solution  *
 * blocks to the specified output data file(s) in a format suitable   *
 * for plotting with TECPLOT.  Returns a non-zero value if cannot     *
 * write any of the TECPLOT solution files.  The experimental data    *
 * was reported by Laufer (19??).                                     *
 *                                                                    *
 **********************************************************************/
void Output_Turbulent_Pipe_Tecplot(NavierStokes2D_Quad_Block &SolnBlk,
				   const int Block_Number,
				   const int Output_Title,
				   const int Output_Data,
				   ostream &Out_File,
				   const double &Re,
				   const double &Pipe_Radius,
				   const int &variable_flag) {

  NavierStokes2D_pState We, W;

  // Output node solution data.  
  Out_File << setprecision(14);

  if (variable_flag == 1) {
    // Axial velocity output file.
    if (Output_Title) {
      Out_File << "TITLE = \"" << CFFC_Name() << ": 2D NavierStokes Turbulent Pipe Flow Axial Velocity "
	       << "\"" << "\n"
	       << "VARIABLES = \"r/R\" \\ \n"
	       << "\"r\" \\ \n"
	       << "\"u\" \\ \n";
      if (Output_Data) {
	Out_File << "ZONE T =  \"Experimental data\" \\ \n"
		 << "I = " << 13 << " \\ \n"
		 << "J = " << 1 << " \\ \n"
		 << "F = POINT \\ \n";
	if (fabs(Re-50000.0) < TOLER) {
	  Out_File << 0.000 << " " << 0.000000 << " " << 3.0480 << endl
		   << 0.100 << " " << 0.012300 << " " << 2.9830 << endl
		   << 0.200 << " " << 0.024600 << " " << 2.9400 << endl
		   << 0.300 << " " << 0.036900 << " " << 2.8800 << endl
		   << 0.400 << " " << 0.049200 << " " << 2.8000 << endl
		   << 0.500 << " " << 0.061500 << " " << 2.7200 << endl
		   << 0.600 << " " << 0.073800 << " " << 2.6000 << endl
		   << 0.700 << " " << 0.086100 << " " << 2.5200 << endl
		   << 0.850 << " " << 0.104550 << " " << 2.4000 << endl
		   << 0.900 << " " << 0.110700 << " " << 2.2000 << endl
		   << 0.930 << " " << 0.114390 << " " << 2.0800 << endl
		   << 0.980 << " " << 0.120540 << " " << 1.1600 << endl
		   << 1.000 << " " << 0.123000 << " " << 0.0000 << endl;
	} else if (fabs(Re-500000.0) < TOLER) {
	  Out_File << 0.000 << " " << 0.000000 << " " << 30.480 << endl
		   << 0.100 << " " << 0.012300 << " " << 30.360 << endl
		   << 0.200 << " " << 0.024600 << " " << 30.000 << endl
		   << 0.300 << " " << 0.036900 << " " << 29.640 << endl
		   << 0.400 << " " << 0.049200 << " " << 28.930 << endl
		   << 0.500 << " " << 0.061500 << " " << 28.570 << endl
		   << 0.600 << " " << 0.073800 << " " << 27.500 << endl
		   << 0.700 << " " << 0.086100 << " " << 26.790 << endl
		   << 0.770 << " " << 0.094710 << " " << 26.070 << endl
		   << 0.840 << " " << 0.103320 << " " << 25.000 << endl
		   << 0.900 << " " << 0.110700 << " " << 23.210 << endl
		   << 0.980 << " " << 0.120540 << " " << 19.640 << endl
		   << 1.000 << " " << 0.123000 << " " <<  0.000 << endl;
	}
      }
    }

  } else if (variable_flag == 2) {
    // Kinetic energy output file.
    if (Output_Title) {
      Out_File << "TITLE = \"" << CFFC_Name() << ": 2D NavierStokes Turbulent Pipe Flow Kinetic Energy "
	       << "\"" << "\n"
	       << "VARIABLES = \"r/R\" \\ \n"
	       << "\"r\" \\ \n"
	       << "\"k\" \\ \n";
      if (Output_Data) {
	if (fabs(Re-50000.0) < TOLER) {
	  Out_File << "ZONE T =  \"Experimental data\" \\ \n"
		   << "I = " << 16 << " \\ \n"
		   << "J = " << 1 << " \\ \n"
		   << "F = POINT \\ \n";
	  Out_File << 0.000 << " " << 0.000000 << " " << 0.01380 << endl
		   << 0.100 << " " << 0.012300 << " " << 0.01420 << endl
		   << 0.200 << " " << 0.024600 << " " << 0.01724 << endl
		   << 0.300 << " " << 0.036900 << " " << 0.01980 << endl
		   << 0.400 << " " << 0.049200 << " " << 0.02460 << endl	
		   << 0.500 << " " << 0.061500 << " " << 0.02930 << endl
		   << 0.600 << " " << 0.073800 << " " << 0.03450 << endl
		   << 0.770 << " " << 0.094710 << " " << 0.04300 << endl
		   << 0.850 << " " << 0.104550 << " " << 0.04660 << endl
		   << 0.930 << " " << 0.114390 << " " << 0.05340 << endl
		   << 0.950 << " " << 0.116850 << " " << 0.06200 << endl
		   << 0.975 << " " << 0.119925 << " " << 0.06550 << endl
		   << 0.983 << " " << 0.120909 << " " << 0.06720 << endl
		   << 0.987 << " " << 0.121400 << " " << 0.05300 << endl
		   << 0.992 << " " << 0.122000 << " " << 0.01640 << endl
		   << 1.000 << " " << 0.123000 << " " << 0.00000 << endl;
	} else if (fabs(Re-500000.0) < TOLER) {
	  Out_File << "ZONE T =  \"Experimental data\" \\ \n"
		   << "I = " << 13 << " \\ \n"
		   << "J = " << 1 << " \\ \n"
		   << "F = POINT \\ \n";
	  Out_File << 0.000 << " " << 0.000000 << " " << 0.97100 << endl
		   << 0.100 << " " << 0.012300 << " " << 0.94000 << endl
		   << 0.200 << " " << 0.024600 << " " << 1.11800 << endl
		   << 0.300 << " " << 0.036900 << " " << 1.35300 << endl
		   << 0.400 << " " << 0.049200 << " " << 1.61800 << endl
		   << 0.500 << " " << 0.061500 << " " << 1.94100 << endl
		   << 0.600 << " " << 0.073800 << " " << 2.35300 << endl
		   << 0.750 << " " << 0.092250 << " " << 2.82400 << endl
		   << 0.850 << " " << 0.104550 << " " << 3.52900 << endl
		   << 0.930 << " " << 0.114390 << " " << 4.26000 << endl
		   << 0.970 << " " << 0.119310 << " " << 4.64700 << endl
		   << 0.990 << " " << 0.121770 << " " << 4.75300 << endl
		   << 1.000 << " " << 0.123000 << " " << 5.58800 << endl;
	}
      }
    }

  } else if (variable_flag == 3) {
    // Shear stress output file.
    if (Output_Title) {
      Out_File << "TITLE = \"" << CFFC_Name() << ": 2D NavierStokes Turbulent Pipe Flow uv "
	       << "\"" << "\n"
	       << "VARIABLES = \"r/R\" \\ \n"
	       << "\"r\" \\ \n"
	       << "\"uv\" \\ \n";
      if (Output_Data) {
	if (fabs(Re-50000.0) < TOLER) {
	  Out_File << "ZONE T =  \"Experimental data\" \\ \n"
		   << "I = " << 11 << " \\ \n"
		   << "J = " << 1 << " \\ \n"
		   << "F = POINT \\ \n";
	  Out_File << 0.000 << " " << 0.00000 << " " << 0.00000 << endl
		   << 0.100 << " " << 0.01230 << " " << 0.00162 << endl
		   << 0.200 << " " << 0.02460 << " " << 0.00353 << endl
		   << 0.300 << " " << 0.03690 << " " << 0.00457 << endl
		   << 0.400 << " " << 0.04920 << " " << 0.00619 << endl	
		   << 0.500 << " " << 0.06150 << " " << 0.00773 << endl
		   << 0.600 << " " << 0.07380 << " " << 0.00854 << endl
		   << 0.720 << " " << 0.08856 << " " << 0.01000 << endl
		   << 0.850 << " " << 0.10455 << " " << 0.01240 << endl
		   << 0.930 << " " << 0.11440 << " " << 0.01330 << endl
		   << 1.000 << " " << 0.12300 << " " << 0.00000 << endl;
	} else if (fabs(Re-500000.0) < TOLER) {
	  Out_File << "ZONE T =  \"Experimental data\" \\ \n"
		   << "I = " << 10 << " \\ \n"
		   << "J = " << 1 << " \\ \n"
		   << "F = POINT \\ \n";
	  Out_File << 0.000 << " " << 0.00000 << " " << 0.00000 << endl
		   << 0.100 << " " << 0.01230 << " " << 0.15000 << endl
		   << 0.200 << " " << 0.02460 << " " << 0.25000 << endl
		   << 0.400 << " " << 0.04920 << " " << 0.50000 << endl
		   << 0.500 << " " << 0.06150 << " " << 0.52400 << endl
		   << 0.600 << " " << 0.07380 << " " << 0.63400 << endl
		   << 0.700 << " " << 0.08610 << " " << 0.87500 << endl
		   << 0.850 << " " << 0.10455 << " " << 1.08000 << endl
		   << 0.950 << " " << 0.11690 << " " << 1.19000 << endl
		   << 1.000 << " " << 0.12300 << " " << 0.00000 << endl;
	}
      }
    }
  }

  Out_File << "ZONE T = \"Block Number = " << Block_Number
	   << "\" \\ \n"
	   << "I = " << SolnBlk.Grid.INu - SolnBlk.Grid.INl + 1 << " \\ \n"
	   << "J = " << SolnBlk.Grid.JNu - SolnBlk.Grid.JNl + 1 << " \\ \n"
	   << "F = POINT \\ \n";
  Out_File.setf(ios::scientific);
  for (int j = SolnBlk.Grid.JNl; j <= SolnBlk.Grid.JNu; j++) {
    for (int i = SolnBlk.Grid.INl; i <= SolnBlk.Grid.INu; i++) {
      W = SolnBlk.Wn(i,j);
      Out_File << SolnBlk.Grid.Node[i][j].X.y/Pipe_Radius << " " 
	       << SolnBlk.Grid.Node[i][j].X.y << " ";
      if (variable_flag == 1) {
	Out_File << W.v.x;
      } else if (variable_flag == 2) {
	Out_File << W.k;
      } else if (variable_flag == 3) {
	Out_File << -SolnBlk.W[i][j].tau.xy;
      }
      Out_File << endl;
    }
  }

}

/**********************************************************************
 * Routine: Output_Flat_Plate_Tecplot                                 *
 *                                                                    *
 * This routine outputs the non-dimensionalized computed flat plate   *
 * solution and the corresponding Blasius solution.                   *
 *                                                                    *
 **********************************************************************/
void Output_Flat_Plate_Tecplot(NavierStokes2D_Quad_Block &SolnBlk,
			       const int Block_Number,
			       const int Output_Title_Soln,
			       ostream &Out_File_Soln,
			       const int Output_Title_Skin,
			       ostream &Out_File_Skin,
			       const NavierStokes2D_pState &Winf,
			       const double &plate_length,
			       double &l1_norm,
			       double &l2_norm,
			       double &max_norm,
			       double &area,
			       int &numberofactivecells,
			       double &l1_norm_cf,
			       double &l2_norm_cf,
			       double &max_norm_cf,
			       double &area_cf,
			       int &numberofactivecells_cf) {

  int skin_friction_flag = OFF;
  NavierStokes2D_pState W, We;
  Vector2D X;
  double eta, f, fp, fpp, Rex, Cf, Cf2, Cfe, phiv;

  Linear_Reconstruction_LeastSquares(SolnBlk,LIMITER_ONE);

  // Output node solution data.  
  Out_File_Soln << setprecision(14);
  if (Output_Title_Soln) {
    Out_File_Soln << "TITLE = \"" << CFFC_Name() << ": 2D NavierStokes Solution, "
		  << "\"" << "\n"
		  << "VARIABLES = \"x\" \\ \n"
		  << "\"y\" \\ \n"
		  << "\"rho\" \\ \n"
		  << "\"u\" \\ \n"
		  << "\"v\" \\ \n"
		  << "\"p\" \\ \n"
		  << "\"rho_e\" \\ \n"
		  << "\"u_e\" \\ \n"
		  << "\"v_e\" \\ \n"
		  << "\"p_e\" \\ \n"
		  << "\"eta\" \\ \n"
		  << "\"f\" \\ \n"
		  << "\"fp\" \\ \n"
		  << "\"fpp\" \\ \n"
		  << "\"u/Uinf\" \\ \n"
		  << "\"u_e/Uinf\" \\ \n"
		  << "\"(u - u_e)/Uinf\" \\ \n"
		  << "\"phiv\" \\ \n"
		  << "\"phiv_e\" \\ \n"
		  << "\"phiv - phiv_e\" \\ \n";
  }
  Out_File_Soln << "ZONE T =  \"Block Number = " << Block_Number << "\" \\ \n"
		<< "I = " << SolnBlk.Grid.INu - SolnBlk.Grid.INl + 1 << " \\ \n"
		<< "J = " << SolnBlk.Grid.JNu - SolnBlk.Grid.JNl + 1 << " \\ \n"
		<< "F = POINT \n";
  for (int j = SolnBlk.Grid.JNl; j <= SolnBlk.Grid.JNu; j++) {
    for (int i = SolnBlk.Grid.INl; i <= SolnBlk.Grid.INu; i++) {
      // Get cell position and solution data.
      X = SolnBlk.Grid.Node[i][j].X;
      W = SolnBlk.Wn(i,j);
      We = FlatPlate(Winf,X,plate_length,eta,f,fp,fpp);
      if (X.x > ZERO && X.x < plate_length) phiv = (W.v.y/Winf.v.x)*sqrt(Winf.v.x*X.x/Winf.nu());
      else phiv = ZERO;
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
  		    << " " << (W.v.x - We.v.x)/Winf.v.x
  		    << " " << phiv
  		    << " " << HALF*(eta*fp-f)
  		    << " " << phiv - HALF*(eta*fp-f)
		    << endl;
      Out_File_Soln.unsetf(ios::scientific);
    }
  }

  // Calculate the norms of the u-velocity component and the skin
  // friction coefficient.
  for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
      // Get cell position and solution data.
      X = SolnBlk.Grid.Cell[i][j].Xc;
      W = SolnBlk.W[i][j];
      // Determine the norms of the u-velocity component.
      if (X.x >= ZERO && X.x <= plate_length) {
	We = FlatPlate(Winf,X,plate_length,eta,f,fp,fpp);
	l1_norm += fabs(W.v.x - We.v.x)*SolnBlk.Grid.Cell[i][j].A;
	l2_norm += sqr(W.v.x - We.v.x)*SolnBlk.Grid.Cell[i][j].A;
	max_norm = max(max_norm,fabs(W.v.x - We.v.x));
	area += SolnBlk.Grid.Cell[i][j].A;
	numberofactivecells++;
      }
      // Determine the norms of the skin friction coefficient.
      if (X.x >= ZERO && X.x <= plate_length && j == 2 &&
	  (SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX ||
	   SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL)) {
	// Get exact skin friction coefficient.
	Rex = (Winf.v.x/Winf.nu())*X.x;
	if (SolnBlk.Flow_Type == FLOWTYPE_LAMINAR) {
	  Cfe = TWO*0.32206/max(sqrt(Rex),NANO);
	} else if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
	  Cfe = 0.37*pow(log(Rex),-2.584);
	}
	// Get computed skin friction coefficient.
	Cf  = TWO*WallShearStress(W,X,
				  SolnBlk.Grid.nodeSW(i,j).X,
				  SolnBlk.Grid.nodeSE(i,j).X,
				  -SolnBlk.Grid.nfaceS(i,j))/(Winf.rho*Winf.v.x*Winf.v.x);
	// Calculate error norms.
	l1_norm_cf += fabs(Cf - Cfe)*SolnBlk.Grid.Cell[i][j].A;
	l2_norm_cf += sqr(Cf - Cfe)*SolnBlk.Grid.Cell[i][j].A;
	max_norm_cf = max(max_norm_cf,fabs(Cf - Cfe));
	area_cf += SolnBlk.Grid.Cell[i][j].A;
	numberofactivecells_cf++;
	if (!skin_friction_flag) skin_friction_flag = ON;
      } else {
	Cf  = ZERO;
	Cfe = ZERO;
      }
    }
  }

  if (Output_Title_Skin) {
    Out_File_Skin << "TITLE = \"" << CFFC_Name() << ": 2D NavierStokes Solution, "
		  << "\"" << "\n"
		  << "VARIABLES = \"x\" \\ \n"
		  << "\"Rex\" \\ \n"
		  << "\"Cf\" \\ \n"
		  << "\"Cf2\" \\ \n"
		  << "\"Cf_e\" \\ \n"
		  << "\"Cf-Cf_e\" \\ \n"
		  << "\"Cf2-Cf_e\" \\ \n";
  }
  if (skin_friction_flag) {
    Out_File_Skin << "ZONE T =  \"Block Number = " << Block_Number
		  << "\" \\ \n"
		  << "I = " << numberofactivecells_cf << " \\ \n"
		  << "J = " << 1 << " \\ \n"
		  << "F = POINT \n";
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
      if (SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc.x >= ZERO &&
	  SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc.x <= plate_length &&
	  (SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX ||
	   SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL)) {
	Rex = (Winf.v.x/Winf.nu())*SolnBlk.Grid.xfaceS(i,SolnBlk.JCl).x;
	if (SolnBlk.Flow_Type == FLOWTYPE_LAMINAR) {
	  Cfe = TWO*0.32206/max(sqrt(Rex),NANO);
	} else if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
	  Cfe = 0.37*pow(log(Rex)/log(TEN),-2.584);
	}
	Cf  = TWO*WallShearStress(SolnBlk.W[i][SolnBlk.JCl],
				  SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc,
				  SolnBlk.Grid.nodeSW(i,SolnBlk.JCl).X,
				  SolnBlk.Grid.nodeSE(i,SolnBlk.JCl).X,
				  -SolnBlk.Grid.nfaceS(i,SolnBlk.JCl))/(Winf.rho*sqr(Winf.v.x));
	Cf2 = TWO*WallShearStress2(SolnBlk.Grid.xfaceS(i,SolnBlk.JCl),
				   SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc,
				   SolnBlk.W[i][SolnBlk.JCl],
				   SolnBlk.dWdx[i][SolnBlk.JCl],
				   SolnBlk.dWdy[i][SolnBlk.JCl],
				   -SolnBlk.Grid.nfaceS(i,SolnBlk.JCl))/(Winf.rho*sqr(Winf.v.x));
	if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
	  Cf2 = TWO*SolnBlk.Wall[i][SolnBlk.JCl].tauw/(Winf.rho*sqr(Winf.v.x));
	}
	Out_File_Skin.setf(ios::scientific);
	Out_File_Skin << " " << SolnBlk.Grid.xfaceS(i,SolnBlk.JCl).x
		      << " " << Rex
		      << " " << Cf
		      << " " << Cf2
		      << " " << Cfe
		      << " " << Cf - Cfe
		      << " " << Cf2 - Cfe
		      << endl;
      }
    }
  }

}

/**********************************************************************
 * Routine: Output_Driven_Cavity_Flow_Tecplot                         *
 *                                                                    *
 * This routine writes a comparison of the computed solution for the  *
 * driven cavity flow with the computations done by Ghia et al. (J.   *
 * Comp. Phys. Vol. 48 1982) for the specified quadrilateral solution *
 * block to the specified output stream suitable for plotting with    *
 * TECPLOT.                                                           *
 *                                                                    *
 **********************************************************************/
void Output_Driven_Cavity_Flow_Tecplot(NavierStokes2D_Quad_Block &SolnBlk,
				       const int Block_Number,
				       const int Output_Title,
				       ostream &Out_File_u,
				       ostream &Out_File_v,
				       const double &Re,
				       const Vector2D &Vwall,
				       const double &length) {

  // Set the output precision and type.
  Out_File_u << setprecision(14);
  Out_File_u.setf(ios::scientific);
  Out_File_v << setprecision(14);
  Out_File_v.setf(ios::scientific);

  // Output node solution data.  
  if (Output_Title) {
    Out_File_u << "TITLE = \"" << CFFC_Name() << ": 2D NavierStokes Driven Cavity Flow Solution u-velocity Comparison "
	       << "\"" << "\n"
	       << "VARIABLES = \"x\" \\ \n"
	       << "\"u\" \\ \n";
    Out_File_v << "TITLE = \"" << CFFC_Name() << ": 2D NavierStokes Driven Cavity Flow Solution v-velocity Comparison "
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
		 << SolnBlk.WnSE(SolnBlk.ICu,j).v.x/Vwall.x << endl;
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
		 << SolnBlk.WnSW(i,SolnBlk.JCl).v.y/Vwall.x << endl;
    }
  }

}

/**********************************************************************
 * Routine: Output_Backward_Facing_Step_Tecplot                       *
 *                                                                    *
 * This routine writes a comparison of the computed solution for the  *
 * backward facing step flow with the experimental data published by  *
 * Driver and Seegmiller (AIAA J. Vol. 23 No. 2 1985) for the         *
 * specified quadrilateral solution block to the specified output     *
 * stream in a format suitable for plotting with tecplot.             *
 *                                                                    *
 **********************************************************************/
void Output_Backward_Facing_Step_Tecplot(NavierStokes2D_Quad_Block &SolnBlk,
					 const int Block_Number,
					 const int Output_Title,
					 ostream &Out_File,
					 const double &step_height,
					 const double &top_wall_deflection) {

//   assert(step_height == 0.0127);
//   assert(top_wall_deflection == ZERO || top_wall_deflection == 6.0);

//   // Set the output precision and type.
//   Out_File << setprecision(14);
//   Out_File.setf(ios::scientific);

//   // Output node solution data.  
//   if (Output_Title) {
//     Out_File << "TITLE = \"" << CFFC_Name() << ": 2D NavierStokes Backward Facing Step \""
// 	     << endl << "VARIABLES = \"x\" \\"
// 	     << endl << "\"y\" \\"
// 	     << endl << "\"u\" \\"
// 	     << endl << "\"v\" \\"
// 	     << endl << "\"k\" \\"
// 	     << endl << "\"-uv\" \\"
// 	     << endl << "\"vuu+vvv\" \\";
//   }

//   // If the output title flag is on (first time the current processor
//   // is writing data) and the current processor is the primary processor
//   // then write the solution data of Driver and Seegmiller to the output
//   // file.
//   if (Output_Title && CFFC_Primary_MPI_Processor()) {

//     Out_File << endl << "ZONE T = \"Station " << 1 << "\" \\"
// 	     << endl << "I = " << 17 << " \\"
// 	     << endl << "J = " << 1 << " \\"
// 	     << endl << " F = POINT \\";
//     Out_File << endl << " " << -4.0 << " " << 1.10 << " " <<  0.657 << " " <<  0.000 << " " << 0.2158245 << " " <<  -2.70 << " " <<  0.08
// 	     << endl << " " << -4.0 << " " << 1.15 << " " <<  0.696 << " " <<  0.000 << " " << 0.2422080 << " " <<  -2.19 << " " <<  0.10
// 	     << endl << " " << -4.0 << " " << 1.20 << " " <<  0.719 << " " <<  0.001 << " " << 0.2584810 << " " <<  -2.23 << " " <<  0.03
// 	     << endl << " " << -4.0 << " " << 1.30 << " " <<  0.760 << " " <<  0.000 << " " << 0.2888000 << " " <<  -2.32 << " " << -0.03
// 	     << endl << " " << -4.0 << " " << 1.40 << " " <<  0.790 << " " <<  0.001 << " " << 0.3120505 << " " <<  -2.10 << " " <<  0.04
// 	     << endl << " " << -4.0 << " " << 1.50 << " " <<  0.818 << " " << -0.001 << " " << 0.3345625 << " " <<  -2.05 << " " <<  0.01
// 	     << endl << " " << -4.0 << " " << 1.70 << " " <<  0.870 << " " <<  0.000 << " " << 0.3784500 << " " <<  -1.64 << " " <<  0.00
// 	     << endl << " " << -4.0 << " " << 2.00 << " " <<  0.926 << " " <<  0.000 << " " << 0.4287380 << " " <<  -1.18 << " " <<  0.01
// 	     << endl << " " << -4.0 << " " << 2.40 << " " <<  0.982 << " " <<  0.003 << " " << 0.4821665 << " " <<  -0.63 << " " << -0.02
// 	     << endl << " " << -4.0 << " " << 2.80 << " " <<  1.003 << " " <<  0.007 << " " << 0.5030290 << " " <<  -0.32 << " " <<  0.01
// 	     << endl << " " << -4.0 << " " << 3.20 << " " <<  1.003 << " " <<  0.006 << " " << 0.5030225 << " " <<  -0.30 << " " <<  0.03
// 	     << endl << " " << -4.0 << " " << 3.60 << " " <<  1.000 << " " <<  0.006 << " " << 0.5000180 << " " <<  -0.28 << " " <<  0.01
// 	     << endl << " " << -4.0 << " " << 4.00 << " " <<  1.001 << " " <<  0.003 << " " << 0.5010050 << " " <<  -0.28 << " " <<  0.01
// 	     << endl << " " << -4.0 << " " << 5.00 << " " <<  1.001 << " " <<  0.000 << " " << 0.5010005 << " " <<  -0.27 << " " <<  0.01
// 	     << endl << " " << -4.0 << " " << 6.00 << " " <<  1.000 << " " << -0.004 << " " << 0.5000080 << " " <<  -0.26 << " " <<  0.01
// 	     << endl << " " << -4.0 << " " << 7.50 << " " <<  1.001 << " " << -0.006 << " " << 0.5010185 << " " <<  -0.28 << " " <<  0.01
// 	     << endl << " " << -4.0 << " " << 8.20 << " " <<  0.943 << " " <<  0.002 << " " << 0.4446265 << " " <<  -1.28 << " " << -0.25;
//     Out_File << endl << "ZONE T = \"Station " << 4 << "\" \\"
// 	     << endl << "I = " << 17 << " \\"
// 	     << endl << "J = " << 1 << " \\"
// 	     << endl << " F = POINT \\";
//     Out_File << endl << " " <<  0.0 << " " << 1.10 << " " <<  0.700 << " " << -0.017 << " " << 0.2451445 << " " <<  -3.54 << " " <<  0.00
// 	     << endl << " " <<  0.0 << " " << 1.15 << " " <<  0.729 << " " << -0.017 << " " << 0.2658650 << " " <<  -2.21 << " " << -0.01
// 	     << endl << " " <<  0.0 << " " << 1.20 << " " <<  0.757 << " " << -0.020 << " " << 0.2867245 << " " <<  -2.14 << " " <<  0.01
// 	     << endl << " " <<  0.0 << " " << 1.30 << " " <<  0.794 << " " << -0.020 << " " << 0.3154180 << " " <<  -2.01 << " " <<  0.04
// 	     << endl << " " <<  0.0 << " " << 1.40 << " " <<  0.825 << " " << -0.020 << " " << 0.3405125 << " " <<  -1.85 << " " <<  0.00
// 	     << endl << " " <<  0.0 << " " << 1.50 << " " <<  0.853 << " " << -0.022 << " " << 0.3640465 << " " <<  -1.80 << " " <<  0.02
// 	     << endl << " " <<  0.0 << " " << 1.70 << " " <<  0.895 << " " << -0.023 << " " << 0.4007770 << " " <<  -1.52 << " " <<  0.02
// 	     << endl << " " <<  0.0 << " " << 2.00 << " " <<  0.947 << " " << -0.022 << " " << 0.4486465 << " " <<  -1.16 << " " << -0.02
// 	     << endl << " " <<  0.0 << " " << 2.40 << " " <<  0.998 << " " << -0.020 << " " << 0.4982020 << " " <<  -0.69 << " " << -0.01
// 	     << endl << " " <<  0.0 << " " << 2.80 << " " <<  1.020 << " " << -0.017 << " " << 0.5203445 << " " <<  -0.43 << " " <<  0.01
// 	     << endl << " " <<  0.0 << " " << 3.20 << " " <<  1.019 << " " << -0.016 << " " << 0.5193085 << " " <<  -0.37 << " " <<  0.02
// 	     << endl << " " <<  0.0 << " " << 3.60 << " " <<  1.016 << " " << -0.016 << " " << 0.5162560 << " " <<  -0.37 << " " <<  0.02
// 	     << endl << " " <<  0.0 << " " << 4.00 << " " <<  1.014 << " " << -0.017 << " " << 0.5142425 << " " <<  -0.39 << " " <<  0.02
// 	     << endl << " " <<  0.0 << " " << 5.00 << " " <<  1.009 << " " << -0.022 << " " << 0.5092825 << " " <<  -0.40 << " " <<  0.01
// 	     << endl << " " <<  0.0 << " " << 6.00 << " " <<  1.007 << " " << -0.020 << " " << 0.5072245 << " " <<  -0.36 << " " <<  0.00
// 	     << endl << " " <<  0.0 << " " << 7.50 << " " <<  1.007 << " " << -0.022 << " " << 0.5072665 << " " <<  -0.39 << " " <<  0.02
// 	     << endl << " " <<  0.0 << " " << 8.20 << " " <<  0.940 << " " << -0.015 << " " << 0.4419125 << " " <<  -1.37 << " " << -0.29;
//     Out_File << endl << "ZONE T = \"Station " << 5 << "\" \\"
// 	     << endl << "I = " << 28 << " \\"
// 	     << endl << "J = " << 1 << " \\"
// 	     << endl << " F = POINT \\";
//     Out_File << endl << " " <<  1.0 << " " << 0.10 << " " << -0.008 << " " << -0.005 << " " << 4.450e-05 << " " <<  -2.94 << " " << -0.16
// 	     << endl << " " <<  1.0 << " " << 0.15 << " " << -0.017 << " " << -0.004 << " " << 0.0001525 << " " <<  -3.22 << " " << -0.08
// 	     << endl << " " <<  1.0 << " " << 0.20 << " " << -0.027 << " " << -0.001 << " " << 0.0003650 << " " <<  -3.29 << " " << -0.12
// 	     << endl << " " <<  1.0 << " " << 0.30 << " " << -0.035 << " " <<  0.004 << " " << 0.0006205 << " " <<  -3.49 << " " << -0.16
// 	     << endl << " " <<  1.0 << " " << 0.40 << " " << -0.035 << " " <<  0.014 << " " << 0.0007105 << " " <<  -3.67 << " " <<  0.07
// 	     << endl << " " <<  1.0 << " " << 0.50 << " " << -0.035 << " " <<  0.022 << " " << 0.0008545 << " " <<  -2.90 << " " << -0.05
// 	     << endl << " " <<  1.0 << " " << 0.60 << " " << -0.027 << " " <<  0.028 << " " << 0.0007565 << " " <<  -3.00 << " " << -0.05
// 	     << endl << " " <<  1.0 << " " << 0.70 << " " << -0.005 << " " <<  0.025 << " " << 0.0003250 << " " <<  -3.44 << " " << -0.19
// 	     << endl << " " <<  1.0 << " " << 0.80 << " " <<  0.081 << " " <<  0.006 << " " << 0.0032985 << " " <<  -7.38 << " " << -0.02
// 	     << endl << " " <<  1.0 << " " << 0.90 << " " <<  0.292 << " " << -0.027 << " " << 0.0429965 << " " << -11.77 << " " << -0.06
// 	     << endl << " " <<  1.0 << " " << 1.00 << " " <<  0.519 << " " << -0.033 << " " << 0.1352250 << " " << -10.41 << " " << -0.41
// 	     << endl << " " <<  1.0 << " " << 1.10 << " " <<  0.704 << " " << -0.023 << " " << 0.2480725 << " " <<  -3.79 << " " <<  0.13
// 	     << endl << " " <<  1.0 << " " << 1.15 << " " <<  0.740 << " " << -0.022 << " " << 0.2740420 << " " <<  -2.84 << " " <<  0.09
// 	     << endl << " " <<  1.0 << " " << 1.20 << " " <<  0.762 << " " << -0.023 << " " << 0.2905865 << " " <<  -2.51 << " " <<  0.13
// 	     << endl << " " <<  1.0 << " " << 1.30 << " " <<  0.793 << " " << -0.031 << " " << 0.3149050 << " " <<  -2.26 << " " <<  0.05
// 	     << endl << " " <<  1.0 << " " << 1.40 << " " <<  0.827 << " " << -0.027 << " " << 0.3423290 << " " <<  -2.05 << " " <<  0.06
// 	     << endl << " " <<  1.0 << " " << 1.50 << " " <<  0.851 << " " << -0.031 << " " << 0.3625810 << " " <<  -1.98 << " " <<  0.02
// 	     << endl << " " <<  1.0 << " " << 1.70 << " " <<  0.891 << " " << -0.034 << " " << 0.3975185 << " " <<  -1.57 << " " << -0.02
// 	     << endl << " " <<  1.0 << " " << 2.00 << " " <<  0.942 << " " << -0.033 << " " << 0.4442265 << " " <<  -1.21 << " " <<  0.00
// 	     << endl << " " <<  1.0 << " " << 2.40 << " " <<  0.990 << " " << -0.031 << " " << 0.4905305 << " " <<  -0.67 << " " <<  0.00
// 	     << endl << " " <<  1.0 << " " << 2.80 << " " <<  1.009 << " " << -0.028 << " " << 0.5094325 << " " <<  -0.36 << " " <<  0.02
// 	     << endl << " " <<  1.0 << " " << 3.20 << " " <<  1.003 << " " << -0.028 << " " << 0.5033965 << " " <<  -0.33 << " " <<  0.02
// 	     << endl << " " <<  1.0 << " " << 3.60 << " " <<  1.002 << " " << -0.026 << " " << 0.5023400 << " " <<  -0.33 << " " <<  0.01
// 	     << endl << " " <<  1.0 << " " << 4.00 << " " <<  1.000 << " " << -0.028 << " " << 0.5003920 << " " <<  -0.30 << " " <<  0.01
// 	     << endl << " " <<  1.0 << " " << 5.00 << " " <<  0.997 << " " << -0.027 << " " << 0.4973690 << " " <<  -0.29 << " " <<  0.01
// 	     << endl << " " <<  1.0 << " " << 6.00 << " " <<  0.994 << " " << -0.027 << " " << 0.4943825 << " " <<  -0.27 << " " <<  0.00
// 	     << endl << " " <<  1.0 << " " << 7.50 << " " <<  0.993 << " " << -0.025 << " " << 0.4933370 << " " <<  -0.30 << " " <<  0.01
// 	     << endl << " " <<  1.0 << " " << 8.20 << " " <<  0.919 << " " << -0.013 << " " << 0.4223650 << " " <<  -1.56 << " " << -0.33;
//     Out_File << endl << "ZONE T = \"Station " << 6 << "\" \\"
// 	     << endl << "I = " << 27 << " \\"
// 	     << endl << "J = " << 1 << " \\"
// 	     << endl << " F = POINT \\";
//     Out_File << endl << " " <<  1.5 << " " << 0.10 << " " << -0.069 << " " << -0.003 << " " << 0.0023850 << " " <<  -2.70 << " " << -0.18
// 	     << endl << " " <<  1.5 << " " << 0.15 << " " << -0.081 << " " <<  0.004 << " " << 0.0032885 << " " <<  -2.27 << " " << -0.10
// 	     << endl << " " <<  1.5 << " " << 0.20 << " " << -0.085 << " " <<  0.009 << " " << 0.0036530 << " " <<  -2.35 << " " << -0.09
// 	     << endl << " " <<  1.5 << " " << 0.30 << " " << -0.088 << " " <<  0.020 << " " << 0.0040720 << " " <<  -2.78 << " " << -0.05
// 	     << endl << " " <<  1.5 << " " << 0.40 << " " << -0.080 << " " <<  0.025 << " " << 0.0035125 << " " <<  -3.13 << " " << -0.21
// 	     << endl << " " <<  1.5 << " " << 0.50 << " " << -0.062 << " " <<  0.025 << " " << 0.0022345 << " " <<  -3.44 << " " << -0.11
// 	     << endl << " " <<  1.5 << " " << 0.60 << " " << -0.026 << " " <<  0.020 << " " << 0.0005380 << " " <<  -4.67 << " " << -0.48
// 	     << endl << " " <<  1.5 << " " << 0.70 << " " <<  0.055 << " " <<  0.000 << " " << 0.0015125 << " " <<  -7.49 << " " << -0.80
// 	     << endl << " " <<  1.5 << " " << 0.80 << " " <<  0.208 << " " << -0.034 << " " << 0.0222100 << " " << -10.56 << " " << -1.49
// 	     << endl << " " <<  1.5 << " " << 0.90 << " " <<  0.373 << " " << -0.046 << " " << 0.0706225 << " " << -11.22 << " " << -0.66
// 	     << endl << " " <<  1.5 << " " << 1.00 << " " <<  0.554 << " " << -0.047 << " " << 0.1545625 << " " <<  -7.99 << " " <<  0.16
// 	     << endl << " " <<  1.5 << " " << 1.10 << " " <<  0.709 << " " << -0.039 << " " << 0.2521010 << " " <<  -3.89 << " " <<  0.39
// 	     << endl << " " <<  1.5 << " " << 1.15 << " " <<  0.739 << " " << -0.037 << " " << 0.2737450 << " " <<  -3.12 << " " <<  0.29
// 	     << endl << " " <<  1.5 << " " << 1.20 << " " <<  0.768 << " " << -0.035 << " " << 0.2955245 << " " <<  -2.35 << " " <<  0.19
// 	     << endl << " " <<  1.5 << " " << 1.30 << " " <<  0.803 << " " << -0.037 << " " << 0.3230890 << " " <<  -2.03 << " " <<  0.09
// 	     << endl << " " <<  1.5 << " " << 1.40 << " " <<  0.832 << " " << -0.037 << " " << 0.3467965 << " " <<  -1.79 << " " <<  0.00
// 	     << endl << " " <<  1.5 << " " << 1.50 << " " <<  0.856 << " " << -0.036 << " " << 0.3670160 << " " <<  -1.61 << " " << -0.03
// 	     << endl << " " <<  1.5 << " " << 1.70 << " " <<  0.898 << " " << -0.036 << " " << 0.4038500 << " " <<  -1.35 << " " <<  0.02
// 	     << endl << " " <<  1.5 << " " << 2.00 << " " <<  0.949 << " " << -0.034 << " " << 0.4508785 << " " <<  -0.96 << " " << -0.02
// 	     << endl << " " <<  1.5 << " " << 2.40 << " " <<  0.997 << " " << -0.030 << " " << 0.4974545 << " " <<  -0.50 << " " << -0.01
// 	     << endl << " " <<  1.5 << " " << 2.80 << " " <<  1.012 << " " << -0.026 << " " << 0.5124100 << " " <<  -0.28 << " " <<  0.00
// 	     << endl << " " <<  1.5 << " " << 3.20 << " " <<  1.009 << " " << -0.025 << " " << 0.5093530 << " " <<  -0.26 << " " <<  0.00
// 	     << endl << " " <<  1.5 << " " << 3.60 << " " <<  1.001 << " " << -0.028 << " " << 0.5013925 << " " <<  -0.24 << " " <<  0.00
// 	     << endl << " " <<  1.5 << " " << 4.00 << " " <<  1.000 << " " << -0.028 << " " << 0.5003920 << " " <<  -0.24 << " " <<  0.00
// 	     << endl << " " <<  1.5 << " " << 5.00 << " " <<  0.996 << " " << -0.027 << " " << 0.4963725 << " " <<  -0.23 << " " <<  0.00
// 	     << endl << " " <<  1.5 << " " << 6.00 << " " <<  0.992 << " " << -0.025 << " " << 0.4923445 << " " <<  -0.23 << " " <<  0.00
// 	     << endl << " " <<  1.5 << " " << 7.50 << " " <<  0.991 << " " << -0.022 << " " << 0.4912825 << " " <<  -0.25 << " " <<  0.00;
//     Out_File << endl << "ZONE T = \"Station " << 7 << "\" \\"
// 	     << endl << "I = " << 28 << " \\"
// 	     << endl << "J = " << 1 << " \\"
// 	     << endl << " F = POINT \\";
//     Out_File << endl << " " <<  2.0 << " " << 0.10 << " " << -0.117 << " " <<  0.013 << " " << 0.0069290 << " " <<  -3.98 << " " << -0.22
// 	     << endl << " " <<  2.0 << " " << 0.15 << " " << -0.132 << " " <<  0.016 << " " << 0.0088400 << " " <<  -3.30 << " " << -0.05
// 	     << endl << " " <<  2.0 << " " << 0.20 << " " << -0.131 << " " <<  0.018 << " " << 0.0087425 << " " <<  -3.40 << " " << -0.06
// 	     << endl << " " <<  2.0 << " " << 0.30 << " " << -0.116 << " " <<  0.025 << " " << 0.0070405 << " " <<  -3.96 << " " << -0.04
// 	     << endl << " " <<  2.0 << " " << 0.40 << " " << -0.093 << " " <<  0.028 << " " << 0.0047165 << " " <<  -4.93 << " " << -0.10
// 	     << endl << " " <<  2.0 << " " << 0.50 << " " << -0.050 << " " <<  0.025 << " " << 0.0015625 << " " <<  -6.33 << " " << -0.09
// 	     << endl << " " <<  2.0 << " " << 0.60 << " " <<  0.027 << " " <<  0.009 << " " << 0.0004050 << " " <<  -8.68 << " " << -0.73
// 	     << endl << " " <<  2.0 << " " << 0.70 << " " <<  0.139 << " " << -0.017 << " " << 0.0098050 << " " << -10.48 << " " << -1.04
// 	     << endl << " " <<  2.0 << " " << 0.80 << " " <<  0.280 << " " << -0.040 << " " << 0.0400000 << " " << -11.68 << " " << -0.57
// 	     << endl << " " <<  2.0 << " " << 0.90 << " " <<  0.406 << " " << -0.048 << " " << 0.0835700 << " " << -11.55 << " " <<  0.11
// 	     << endl << " " <<  2.0 << " " << 1.00 << " " <<  0.564 << " " << -0.050 << " " << 0.1602980 << " " <<  -9.31 << " " << -0.59
// 	     << endl << " " <<  2.0 << " " << 1.10 << " " <<  0.701 << " " << -0.046 << " " << 0.2467585 << " " <<  -5.20 << " " <<  0.06
// 	     << endl << " " <<  2.0 << " " << 1.15 << " " <<  0.735 << " " << -0.045 << " " << 0.2711250 << " " <<  -4.04 << " " <<  0.12
// 	     << endl << " " <<  2.0 << " " << 1.20 << " " <<  0.769 << " " << -0.044 << " " << 0.2966485 << " " <<  -2.88 << " " <<  0.18
// 	     << endl << " " <<  2.0 << " " << 1.30 << " " <<  0.806 << " " << -0.044 << " " << 0.3257860 << " " <<  -2.30 << " " <<  0.07
// 	     << endl << " " <<  2.0 << " " << 1.40 << " " <<  0.835 << " " << -0.045 << " " << 0.3496250 << " " <<  -2.04 << " " <<  0.03
// 	     << endl << " " <<  2.0 << " " << 1.50 << " " <<  0.858 << " " << -0.045 << " " << 0.3690945 << " " <<  -1.79 << " " << -0.01
// 	     << endl << " " <<  2.0 << " " << 1.70 << " " <<  0.898 << " " << -0.046 << " " << 0.4042600 << " " <<  -1.50 << " " << -0.02
// 	     << endl << " " <<  2.0 << " " << 2.00 << " " <<  0.944 << " " << -0.045 << " " << 0.4465805 << " " <<  -1.11 << " " << -0.02
// 	     << endl << " " <<  2.0 << " " << 2.40 << " " <<  0.990 << " " << -0.043 << " " << 0.4909745 << " " <<  -0.62 << " " << -0.01
// 	     << endl << " " <<  2.0 << " " << 2.80 << " " <<  1.005 << " " << -0.038 << " " << 0.5057345 << " " <<  -0.35 << " " <<  0.00
// 	     << endl << " " <<  2.0 << " " << 3.20 << " " <<  1.000 << " " << -0.036 << " " << 0.5006480 << " " <<  -0.31 << " " << -0.01
// 	     << endl << " " <<  2.0 << " " << 3.60 << " " <<  0.996 << " " << -0.034 << " " << 0.4965860 << " " <<  -0.30 << " " <<  0.01
// 	     << endl << " " <<  2.0 << " " << 4.00 << " " <<  0.996 << " " << -0.033 << " " << 0.4965525 << " " <<  -0.29 << " " <<  0.00
// 	     << endl << " " <<  2.0 << " " << 5.00 << " " <<  0.992 << " " << -0.031 << " " << 0.4925125 << " " <<  -0.29 << " " <<  0.00
// 	     << endl << " " <<  2.0 << " " << 6.00 << " " <<  0.988 << " " << -0.029 << " " << 0.4884925 << " " <<  -0.29 << " " <<  0.00
// 	     << endl << " " <<  2.0 << " " << 7.50 << " " <<  0.987 << " " << -0.025 << " " << 0.4873970 << " " <<  -0.29 << " " <<  0.00
// 	     << endl << " " <<  2.0 << " " << 8.20 << " " <<  0.911 << " " << -0.014 << " " << 0.4150585 << " " <<  -1.61 << " " << -0.35;
//     Out_File << endl << "ZONE T = \"Station " << 8 << "\" \\"
// 	     << endl << "I = " << 28 << " \\"
// 	     << endl << "J = " << 1 << " \\"
// 	     << endl << " F = POINT \\";
//     Out_File << endl << " " <<  2.5 << " " << 0.10 << " " << -0.177 << " " <<  0.008 << " " << 0.0156965 << " " <<  -4.21 << " " <<  0.13
// 	     << endl << " " <<  2.5 << " " << 0.15 << " " << -0.180 << " " <<  0.008 << " " << 0.0162320 << " " <<  -4.18 << " " << -0.26
// 	     << endl << " " <<  2.5 << " " << 0.20 << " " << -0.171 << " " <<  0.012 << " " << 0.0146925 << " " <<  -4.73 << " " << -0.34
// 	     << endl << " " <<  2.5 << " " << 0.30 << " " << -0.141 << " " <<  0.016 << " " << 0.0100685 << " " <<  -6.31 << " " << -0.57
// 	     << endl << " " <<  2.5 << " " << 0.40 << " " << -0.095 << " " <<  0.013 << " " << 0.0045970 << " " <<  -7.92 << " " << -0.30
// 	     << endl << " " <<  2.5 << " " << 0.50 << " " << -0.012 << " " <<  0.002 << " " << 7.400e-05 << " " << -10.50 << " " << -0.76
// 	     << endl << " " <<  2.5 << " " << 0.60 << " " <<  0.093 << " " << -0.015 << " " << 0.0044370 << " " << -12.43 << " " << -2.12
// 	     << endl << " " <<  2.5 << " " << 0.70 << " " <<  0.226 << " " << -0.038 << " " << 0.0262600 << " " << -12.10 << " " << -1.84
// 	     << endl << " " <<  2.5 << " " << 0.80 << " " <<  0.341 << " " << -0.049 << " " << 0.0593410 << " " << -11.95 << " " << -0.01
// 	     << endl << " " <<  2.5 << " " << 0.90 << " " <<  0.461 << " " << -0.056 << " " << 0.1078285 << " " << -10.77 << " " <<  0.29
// 	     << endl << " " <<  2.5 << " " << 1.00 << " " <<  0.599 << " " << -0.056 << " " << 0.1809685 << " " <<  -8.36 << " " << -0.31
// 	     << endl << " " <<  2.5 << " " << 1.10 << " " <<  0.707 << " " << -0.052 << " " << 0.2512765 << " " <<  -5.32 << " " <<  0.13
// 	     << endl << " " <<  2.5 << " " << 1.15 << " " <<  0.740 << " " << -0.053 << " " << 0.2752045 << " " <<  -4.09 << " " <<  0.16
// 	     << endl << " " <<  2.5 << " " << 1.20 << " " <<  0.774 << " " << -0.054 << " " << 0.3009960 << " " <<  -2.87 << " " <<  0.21
// 	     << endl << " " <<  2.5 << " " << 1.30 << " " <<  0.811 << " " << -0.054 << " " << 0.3303185 << " " <<  -2.29 << " " <<  0.10
// 	     << endl << " " <<  2.5 << " " << 1.40 << " " <<  0.838 << " " << -0.053 << " " << 0.3525265 << " " <<  -1.95 << " " << -0.01
// 	     << endl << " " <<  2.5 << " " << 1.50 << " " <<  0.863 << " " << -0.055 << " " << 0.3738970 << " " <<  -1.73 << " " << -0.02
// 	     << endl << " " <<  2.5 << " " << 1.70 << " " <<  0.902 << " " << -0.055 << " " << 0.4083145 << " " <<  -1.40 << " " <<  0.02
// 	     << endl << " " <<  2.5 << " " << 2.00 << " " <<  0.948 << " " << -0.054 << " " << 0.4508100 << " " <<  -0.97 << " " << -0.02
// 	     << endl << " " <<  2.5 << " " << 2.40 << " " <<  0.992 << " " << -0.049 << " " << 0.4932325 << " " <<  -0.53 << " " << -0.01
// 	     << endl << " " <<  2.5 << " " << 2.80 << " " <<  1.006 << " " << -0.044 << " " << 0.5069860 << " " <<  -0.34 << " " <<  0.00
// 	     << endl << " " <<  2.5 << " " << 3.20 << " " <<  1.001 << " " << -0.040 << " " << 0.5018005 << " " <<  -0.29 << " " <<  0.01
// 	     << endl << " " <<  2.5 << " " << 3.60 << " " <<  1.000 << " " << -0.036 << " " << 0.5006480 << " " <<  -0.26 << " " <<  0.01
// 	     << endl << " " <<  2.5 << " " << 4.00 << " " <<  0.999 << " " << -0.032 << " " << 0.4995125 << " " <<  -0.27 << " " <<  0.01
// 	     << endl << " " <<  2.5 << " " << 5.00 << " " <<  0.995 << " " << -0.028 << " " << 0.4954045 << " " <<  -0.28 << " " <<  0.00
// 	     << endl << " " <<  2.5 << " " << 6.00 << " " <<  0.990 << " " << -0.024 << " " << 0.4903380 << " " <<  -0.29 << " " <<  0.00
// 	     << endl << " " <<  2.5 << " " << 7.50 << " " <<  0.992 << " " << -0.019 << " " << 0.4922125 << " " <<  -0.28 << " " <<  0.01
// 	     << endl << " " <<  2.5 << " " << 8.20 << " " <<  0.914 << " " << -0.010 << " " << 0.4177480 << " " <<  -1.46 << " " << -0.32;
//     Out_File << endl << "ZONE T = \"Station " << 9 << "\" \\"
// 	     << endl << "I = " << 28 << " \\"
// 	     << endl << "J = " << 1 << " \\"
// 	     << endl << " F = POINT \\";
//     Out_File << endl << " " <<  3.0 << " " << 0.10 << " " << -0.198 << " " << -0.003 << " " << 0.0196065 << " " <<  -6.10 << " " <<  0.12
// 	     << endl << " " <<  3.0 << " " << 0.15 << " " << -0.197 << " " <<  0.000 << " " << 0.0194045 << " " <<  -5.93 << " " << -0.47
// 	     << endl << " " <<  3.0 << " " << 0.20 << " " << -0.184 << " " << -0.001 << " " << 0.0169285 << " " <<  -7.16 << " " << -0.33
// 	     << endl << " " <<  3.0 << " " << 0.30 << " " << -0.143 << " " << -0.002 << " " << 0.0102265 << " " <<  -9.75 << " " << -0.09
// 	     << endl << " " <<  3.0 << " " << 0.40 << " " << -0.066 << " " << -0.014 << " " << 0.0022760 << " " << -13.31 << " " <<  0.47
// 	     << endl << " " <<  3.0 << " " << 0.50 << " " <<  0.038 << " " << -0.029 << " " << 0.0011425 << " " << -15.51 << " " << -2.08
// 	     << endl << " " <<  3.0 << " " << 0.60 << " " <<  0.162 << " " << -0.046 << " " << 0.0141800 << " " << -15.29 << " " << -3.61
// 	     << endl << " " <<  3.0 << " " << 0.70 << " " <<  0.282 << " " << -0.057 << " " << 0.0413865 << " " << -13.35 << " " << -1.94
// 	     << endl << " " <<  3.0 << " " << 0.80 << " " <<  0.388 << " " << -0.064 << " " << 0.0773200 << " " << -11.63 << " " << -0.18
// 	     << endl << " " <<  3.0 << " " << 0.90 << " " <<  0.506 << " " << -0.068 << " " << 0.1303300 << " " << -10.25 << " " <<  0.35
// 	     << endl << " " <<  3.0 << " " << 1.00 << " " <<  0.613 << " " << -0.064 << " " << 0.1899325 << " " <<  -8.06 << " " <<  0.10
// 	     << endl << " " <<  3.0 << " " << 1.10 << " " <<  0.716 << " " << -0.067 << " " << 0.2585725 << " " <<  -4.80 << " " <<  0.01
// 	     << endl << " " <<  3.0 << " " << 1.15 << " " <<  0.747 << " " << -0.067 << " " << 0.2812490 << " " <<  -3.94 << " " <<  0.06
// 	     << endl << " " <<  3.0 << " " << 1.20 << " " <<  0.777 << " " << -0.066 << " " << 0.3040425 << " " <<  -3.07 << " " <<  0.12
// 	     << endl << " " <<  3.0 << " " << 1.30 << " " <<  0.813 << " " << -0.064 << " " << 0.3325325 << " " <<  -2.28 << " " <<  0.13
// 	     << endl << " " <<  3.0 << " " << 1.40 << " " <<  0.841 << " " << -0.063 << " " << 0.3556250 << " " <<  -1.81 << " " <<  0.06
// 	     << endl << " " <<  3.0 << " " << 1.50 << " " <<  0.864 << " " << -0.063 << " " << 0.3752325 << " " <<  -1.55 << " " <<  0.01
// 	     << endl << " " <<  3.0 << " " << 1.70 << " " <<  0.903 << " " << -0.063 << " " << 0.4096890 << " " <<  -1.30 << " " << -0.01
// 	     << endl << " " <<  3.0 << " " << 2.00 << " " <<  0.948 << " " << -0.059 << " " << 0.4510925 << " " <<  -0.83 << " " << -0.01
// 	     << endl << " " <<  3.0 << " " << 2.40 << " " <<  0.989 << " " << -0.055 << " " << 0.4905730 << " " <<  -0.40 << " " << -0.01
// 	     << endl << " " <<  3.0 << " " << 2.80 << " " <<  0.998 << " " << -0.048 << " " << 0.4991540 << " " <<  -0.25 << " " <<  0.00
// 	     << endl << " " <<  3.0 << " " << 3.20 << " " <<  0.995 << " " << -0.044 << " " << 0.4959805 << " " <<  -0.23 << " " <<  0.00
// 	     << endl << " " <<  3.0 << " " << 3.60 << " " <<  0.994 << " " << -0.042 << " " << 0.4949000 << " " <<  -0.22 << " " <<  0.00
// 	     << endl << " " <<  3.0 << " " << 4.00 << " " <<  0.992 << " " << -0.039 << " " << 0.4927925 << " " <<  -0.21 << " " <<  0.00
// 	     << endl << " " <<  3.0 << " " << 5.00 << " " <<  0.988 << " " << -0.034 << " " << 0.4886500 << " " <<  -0.21 << " " <<  0.00
// 	     << endl << " " <<  3.0 << " " << 6.00 << " " <<  0.984 << " " << -0.030 << " " << 0.4845780 << " " <<  -0.22 << " " <<  0.00
// 	     << endl << " " <<  3.0 << " " << 7.50 << " " <<  0.983 << " " << -0.025 << " " << 0.4834570 << " " <<  -0.23 << " " <<  0.01
// 	     << endl << " " <<  3.0 << " " << 8.20 << " " <<  0.899 << " " << -0.015 << " " << 0.4042130 << " " <<  -1.55 << " " << -0.40;
//     Out_File << endl << "ZONE T = \"Station " << 10 << "\" \\"
// 	     << endl << "I = " << 28 << " \\"
// 	     << endl << "J = " << 1 << " \\"
// 	     << endl << " F = POINT \\";
//     Out_File << endl << " " <<  4.0 << " " << 0.10 << " " << -0.171 << " " << -0.014 << " " << 0.0147185 << " " <<  -7.67 << " " <<  0.02
// 	     << endl << " " <<  4.0 << " " << 0.15 << " " << -0.154 << " " << -0.020 << " " << 0.0120580 << " " << -10.10 << " " << -0.24
// 	     << endl << " " <<  4.0 << " " << 0.20 << " " << -0.134 << " " << -0.027 << " " << 0.0093425 << " " << -12.98 << " " <<  0.27
// 	     << endl << " " <<  4.0 << " " << 0.30 << " " << -0.057 << " " << -0.046 << " " << 0.0026825 << " " << -16.16 << " " <<  0.31
// 	     << endl << " " <<  4.0 << " " << 0.40 << " " <<  0.052 << " " << -0.068 << " " << 0.0036640 << " " << -18.95 << " " << -0.61
// 	     << endl << " " <<  4.0 << " " << 0.50 << " " <<  0.169 << " " << -0.082 << " " << 0.0176425 << " " << -18.32 << " " << -2.96
// 	     << endl << " " <<  4.0 << " " << 0.60 << " " <<  0.279 << " " << -0.093 << " " << 0.0432450 << " " << -15.90 << " " << -2.24
// 	     << endl << " " <<  4.0 << " " << 0.70 << " " <<  0.381 << " " << -0.094 << " " << 0.0769985 << " " << -13.03 << " " << -0.20
// 	     << endl << " " <<  4.0 << " " << 0.80 << " " <<  0.473 << " " << -0.097 << " " << 0.1165690 << " " << -10.98 << " " <<  0.51
// 	     << endl << " " <<  4.0 << " " << 0.90 << " " <<  0.568 << " " << -0.097 << " " << 0.1660165 << " " <<  -8.30 << " " <<  0.06
// 	     << endl << " " <<  4.0 << " " << 1.00 << " " <<  0.656 << " " << -0.092 << " " << 0.2194000 << " " <<  -6.33 << " " << -0.23
// 	     << endl << " " <<  4.0 << " " << 1.10 << " " <<  0.730 << " " << -0.089 << " " << 0.2704105 << " " <<  -4.20 << " " <<  0.10
// 	     << endl << " " <<  4.0 << " " << 1.15 << " " <<  0.754 << " " << -0.087 << " " << 0.2880425 << " " <<  -3.52 << " " <<  0.20
// 	     << endl << " " <<  4.0 << " " << 1.20 << " " <<  0.778 << " " << -0.084 << " " << 0.3061700 << " " <<  -2.83 << " " <<  0.29
// 	     << endl << " " <<  4.0 << " " << 1.30 << " " <<  0.811 << " " << -0.083 << " " << 0.3323050 << " " <<  -2.19 << " " <<  0.13
// 	     << endl << " " <<  4.0 << " " << 1.40 << " " <<  0.835 << " " << -0.081 << " " << 0.3518930 << " " <<  -1.80 << " " <<  0.09
// 	     << endl << " " <<  4.0 << " " << 1.50 << " " <<  0.856 << " " << -0.079 << " " << 0.3694885 << " " <<  -1.52 << " " <<  0.06
// 	     << endl << " " <<  4.0 << " " << 1.70 << " " <<  0.895 << " " << -0.076 << " " << 0.4034005 << " " <<  -1.25 << " " <<  0.03
// 	     << endl << " " <<  4.0 << " " << 2.00 << " " <<  0.940 << " " << -0.072 << " " << 0.4443920 << " " <<  -0.78 << " " <<  0.00
// 	     << endl << " " <<  4.0 << " " << 2.40 << " " <<  0.978 << " " << -0.063 << " " << 0.4802265 << " " <<  -0.39 << " " <<  0.01
// 	     << endl << " " <<  4.0 << " " << 2.80 << " " <<  0.986 << " " << -0.056 << " " << 0.4876660 << " " <<  -0.26 << " " <<  0.00
// 	     << endl << " " <<  4.0 << " " << 3.20 << " " <<  0.984 << " " << -0.050 << " " << 0.4853780 << " " <<  -0.23 << " " <<  0.00
// 	     << endl << " " <<  4.0 << " " << 3.60 << " " <<  0.985 << " " << -0.047 << " " << 0.4862170 << " " <<  -0.23 << " " <<  0.00
// 	     << endl << " " <<  4.0 << " " << 4.00 << " " <<  0.985 << " " << -0.043 << " " << 0.4860370 << " " <<  -0.21 << " " <<  0.00
// 	     << endl << " " <<  4.0 << " " << 5.00 << " " <<  0.983 << " " << -0.037 << " " << 0.4838290 << " " <<  -0.22 << " " <<  0.00
// 	     << endl << " " <<  4.0 << " " << 6.00 << " " <<  0.980 << " " << -0.031 << " " << 0.4806805 << " " <<  -0.23 << " " <<  0.00
// 	     << endl << " " <<  4.0 << " " << 7.50 << " " <<  0.980 << " " << -0.025 << " " << 0.4805125 << " " <<  -0.25 << " " <<  0.00
// 	     << endl << " " <<  4.0 << " " << 8.20 << " " <<  0.892 << " " << -0.013 << " " << 0.3979165 << " " <<  -1.59 << " " << -0.35;
//     Out_File << endl << "ZONE T = \"Station " << 11 << "\" \\"
// 	     << endl << "I = " << 28 << " \\"
// 	     << endl << "J = " << 1 << " \\"
// 	     << endl << " F = POINT \\";
//     Out_File << endl << " " <<  5.0 << " " << 0.10 << " " << -0.032 << " " << -0.046 << " " << 0.0015700 << " " << -13.16 << " " <<  1.15
// 	     << endl << " " <<  5.0 << " " << 0.15 << " " << -0.006 << " " << -0.053 << " " << 0.0014225 << " " << -15.57 << " " <<  0.42
// 	     << endl << " " <<  5.0 << " " << 0.20 << " " <<  0.026 << " " << -0.063 << " " << 0.0023225 << " " << -18.23 << " " <<  0.72
// 	     << endl << " " <<  5.0 << " " << 0.30 << " " <<  0.117 << " " << -0.088 << " " << 0.0107165 << " " << -20.49 << " " << -0.49
// 	     << endl << " " <<  5.0 << " " << 0.40 << " " <<  0.204 << " " << -0.097 << " " << 0.0255125 << " " << -19.81 << " " << -1.22
// 	     << endl << " " <<  5.0 << " " << 0.50 << " " <<  0.293 << " " << -0.103 << " " << 0.0482290 << " " << -17.30 << " " << -2.13
// 	     << endl << " " <<  5.0 << " " << 0.60 << " " <<  0.375 << " " << -0.106 << " " << 0.0759305 << " " << -15.07 << " " << -0.08
// 	     << endl << " " <<  5.0 << " " << 0.70 << " " <<  0.452 << " " << -0.104 << " " << 0.1075600 << " " << -13.26 << " " << -0.09
// 	     << endl << " " <<  5.0 << " " << 0.80 << " " <<  0.536 << " " << -0.103 << " " << 0.1489525 << " " << -10.57 << " " << -0.21
// 	     << endl << " " <<  5.0 << " " << 0.90 << " " <<  0.610 << " " << -0.097 << " " << 0.1907545 << " " <<  -8.40 << " " << -0.05
// 	     << endl << " " <<  5.0 << " " << 1.00 << " " <<  0.679 << " " << -0.096 << " " << 0.2351285 << " " <<  -6.14 << " " <<  0.60
// 	     << endl << " " <<  5.0 << " " << 1.10 << " " <<  0.734 << " " << -0.095 << " " << 0.2738905 << " " <<  -4.05 << " " <<  0.29
// 	     << endl << " " <<  5.0 << " " << 1.15 << " " <<  0.753 << " " << -0.093 << " " << 0.2878290 << " " <<  -3.50 << " " <<  0.30
// 	     << endl << " " <<  5.0 << " " << 1.20 << " " <<  0.771 << " " << -0.091 << " " << 0.3013610 << " " <<  -2.95 << " " <<  0.33
// 	     << endl << " " <<  5.0 << " " << 1.30 << " " <<  0.800 << " " << -0.090 << " " << 0.3240500 << " " <<  -2.12 << " " <<  0.18
// 	     << endl << " " <<  5.0 << " " << 1.40 << " " <<  0.826 << " " << -0.088 << " " << 0.3450100 << " " <<  -1.85 << " " <<  0.10
// 	     << endl << " " <<  5.0 << " " << 1.50 << " " <<  0.845 << " " << -0.086 << " " << 0.3607105 << " " <<  -1.61 << " " <<  0.07
// 	     << endl << " " <<  5.0 << " " << 1.70 << " " <<  0.880 << " " << -0.082 << " " << 0.3905620 << " " <<  -1.25 << " " <<  0.03
// 	     << endl << " " <<  5.0 << " " << 2.00 << " " <<  0.926 << " " << -0.076 << " " << 0.4316260 << " " <<  -0.83 << " " << -0.01
// 	     << endl << " " <<  5.0 << " " << 2.40 << " " <<  0.965 << " " << -0.068 << " " << 0.4679245 << " " <<  -0.44 << " " <<  0.02
// 	     << endl << " " <<  5.0 << " " << 2.80 << " " <<  0.973 << " " << -0.060 << " " << 0.4751645 << " " <<  -0.33 << " " <<  0.03
// 	     << endl << " " <<  5.0 << " " << 3.20 << " " <<  0.972 << " " << -0.054 << " " << 0.4738500 << " " <<  -0.29 << " " <<  0.02
// 	     << endl << " " <<  5.0 << " " << 3.60 << " " <<  0.973 << " " << -0.050 << " " << 0.4746145 << " " <<  -0.30 << " " <<  0.03
// 	     << endl << " " <<  5.0 << " " << 4.00 << " " <<  0.974 << " " << -0.046 << " " << 0.4753960 << " " <<  -0.29 << " " <<  0.03
// 	     << endl << " " <<  5.0 << " " << 5.00 << " " <<  0.974 << " " << -0.040 << " " << 0.4751380 << " " <<  -0.28 << " " <<  0.02
// 	     << endl << " " <<  5.0 << " " << 6.00 << " " <<  0.975 << " " << -0.034 << " " << 0.4758905 << " " <<  -0.28 << " " <<  0.01
// 	     << endl << " " <<  5.0 << " " << 7.50 << " " <<  0.975 << " " << -0.026 << " " << 0.4756505 << " " <<  -0.32 << " " <<  0.02
// 	     << endl << " " <<  5.0 << " " << 8.20 << " " <<  0.880 << " " << -0.018 << " " << 0.3873620 << " " <<  -1.91 << " " << -0.47;
//     Out_File << endl << "ZONE T = \"Station " << 12 << "\" \\"
// 	     << endl << "I = " << 28 << " \\"
// 	     << endl << "J = " << 1 << " \\"
// 	     << endl << " F = POINT \\";
//     Out_File << endl << " " <<  5.5 << " " << 0.10 << " " <<  0.053 << " " << -0.044 << " " << 0.0023725 << " " << -13.38 << " " << -0.71
// 	     << endl << " " <<  5.5 << " " << 0.15 << " " <<  0.064 << " " << -0.050 << " " << 0.0032980 << " " << -15.15 << " " << -0.22
// 	     << endl << " " <<  5.5 << " " << 0.20 << " " <<  0.100 << " " << -0.062 << " " << 0.0069220 << " " << -17.42 << " " << -1.16
// 	     << endl << " " <<  5.5 << " " << 0.30 << " " <<  0.178 << " " << -0.081 << " " << 0.0191225 << " " << -18.91 << " " << -0.33
// 	     << endl << " " <<  5.5 << " " << 0.40 << " " <<  0.254 << " " << -0.095 << " " << 0.0367705 << " " << -17.38 << " " << -0.85
// 	     << endl << " " <<  5.5 << " " << 0.50 << " " <<  0.327 << " " << -0.103 << " " << 0.0587690 << " " << -15.32 << " " << -0.53
// 	     << endl << " " <<  5.5 << " " << 0.60 << " " <<  0.406 << " " << -0.104 << " " << 0.0878260 << " " << -13.37 << " " <<  0.31
// 	     << endl << " " <<  5.5 << " " << 0.70 << " " <<  0.485 << " " << -0.104 << " " << 0.1230205 << " " << -11.04 << " " << -0.01
// 	     << endl << " " <<  5.5 << " " << 0.80 << " " <<  0.560 << " " << -0.101 << " " << 0.1619005 << " " <<  -8.80 << " " << -0.05
// 	     << endl << " " <<  5.5 << " " << 0.90 << " " <<  0.626 << " " << -0.097 << " " << 0.2006425 << " " <<  -6.91 << " " << -0.28
// 	     << endl << " " <<  5.5 << " " << 1.00 << " " <<  0.684 << " " << -0.092 << " " << 0.2381600 << " " <<  -5.15 << " " << -0.09
// 	     << endl << " " <<  5.5 << " " << 1.10 << " " <<  0.732 << " " << -0.089 << " " << 0.2718725 << " " <<  -3.66 << " " <<  0.03
// 	     << endl << " " <<  5.5 << " " << 1.15 << " " <<  0.748 << " " << -0.089 << " " << 0.2837125 << " " <<  -3.18 << " " <<  0.05
// 	     << endl << " " <<  5.5 << " " << 1.20 << " " <<  0.765 << " " << -0.089 << " " << 0.2965730 << " " <<  -2.70 << " " <<  0.06
// 	     << endl << " " <<  5.5 << " " << 1.30 << " " <<  0.792 << " " << -0.085 << " " << 0.3172445 << " " <<  -2.28 << " " <<  0.11
// 	     << endl << " " <<  5.5 << " " << 1.40 << " " <<  0.816 << " " << -0.084 << " " << 0.3364560 << " " <<  -1.99 << " " <<  0.05
// 	     << endl << " " <<  5.5 << " " << 1.50 << " " <<  0.836 << " " << -0.081 << " " << 0.3527285 << " " <<  -1.72 << " " <<  0.04
// 	     << endl << " " <<  5.5 << " " << 1.70 << " " <<  0.872 << " " << -0.079 << " " << 0.3833125 << " " <<  -1.30 << " " << -0.01
// 	     << endl << " " <<  5.5 << " " << 2.00 << " " <<  0.918 << " " << -0.074 << " " << 0.4241000 << " " <<  -0.83 << " " << -0.01
// 	     << endl << " " <<  5.5 << " " << 2.40 << " " <<  0.957 << " " << -0.066 << " " << 0.4601025 << " " <<  -0.42 << " " << -0.01
// 	     << endl << " " <<  5.5 << " " << 2.80 << " " <<  0.965 << " " << -0.058 << " " << 0.4672945 << " " <<  -0.29 << " " <<  0.00
// 	     << endl << " " <<  5.5 << " " << 3.20 << " " <<  0.966 << " " << -0.052 << " " << 0.4679300 << " " <<  -0.27 << " " <<  0.00
// 	     << endl << " " <<  5.5 << " " << 3.60 << " " <<  0.968 << " " << -0.049 << " " << 0.4697125 << " " <<  -0.26 << " " <<  0.00
// 	     << endl << " " <<  5.5 << " " << 4.00 << " " <<  0.969 << " " << -0.044 << " " << 0.4704485 << " " <<  -0.25 << " " <<  0.00
// 	     << endl << " " <<  5.5 << " " << 5.00 << " " <<  0.969 << " " << -0.038 << " " << 0.4702025 << " " <<  -0.26 << " " <<  0.00
// 	     << endl << " " <<  5.5 << " " << 6.00 << " " <<  0.970 << " " << -0.032 << " " << 0.4709620 << " " <<  -0.26 << " " <<  0.00
// 	     << endl << " " <<  5.5 << " " << 7.50 << " " <<  0.971 << " " << -0.025 << " " << 0.4717330 << " " <<  -0.29 << " " <<  0.00
// 	     << endl << " " <<  5.5 << " " << 8.20 << " " <<  0.874 << " " << -0.016 << " " << 0.3820660 << " " <<  -1.75 << " " << -0.43;
//     Out_File << endl << "ZONE T = \"Station " << 13 << "\" \\"
// 	     << endl << "I = " << 28 << " \\"
// 	     << endl << "J = " << 1 << " \\"
// 	     << endl << " F = POINT \\";
//     Out_File << endl << " " <<  6.0 << " " << 0.10 << " " << 0.147 << " " << -0.052 << " " << 0.0121565 << " " << -12.81 << " " << -1.40
// 	     << endl << " " <<  6.0 << " " << 0.15 << " " << 0.169 << " " << -0.061 << " " << 0.0161410 << " " << -15.02 << " " << -1.07
// 	     << endl << " " <<  6.0 << " " << 0.20 << " " << 0.188 << " " << -0.067 << " " << 0.0199165 << " " << -16.50 << " " << -1.12
// 	     << endl << " " <<  6.0 << " " << 0.30 << " " << 0.248 << " " << -0.084 << " " << 0.0342800 << " " << -16.98 << " " << -1.03
// 	     << endl << " " <<  6.0 << " " << 0.40 << " " << 0.306 << " " << -0.091 << " " << 0.0509585 << " " << -16.60 << " " << -0.18
// 	     << endl << " " <<  6.0 << " " << 0.50 << " " << 0.369 << " " << -0.097 << " " << 0.0727850 << " " << -15.89 << " " << -0.87
// 	     << endl << " " <<  6.0 << " " << 0.60 << " " << 0.435 << " " << -0.099 << " " << 0.0995130 << " " << -13.86 << " " << -0.38
// 	     << endl << " " <<  6.0 << " " << 0.70 << " " << 0.501 << " " << -0.095 << " " << 0.1300130 << " " << -12.00 << " " << -0.25
// 	     << endl << " " <<  6.0 << " " << 0.80 << " " << 0.569 << " " << -0.094 << " " << 0.1662985 << " " <<  -9.63 << " " << -0.36
// 	     << endl << " " <<  6.0 << " " << 0.90 << " " << 0.632 << " " << -0.091 << " " << 0.2038525 << " " <<  -7.62 << " " << -0.65
// 	     << endl << " " <<  6.0 << " " << 1.00 << " " << 0.686 << " " << -0.089 << " " << 0.2392585 << " " <<  -5.57 << " " << -0.52
// 	     << endl << " " <<  6.0 << " " << 1.10 << " " << 0.730 << " " << -0.086 << " " << 0.2701480 << " " <<  -3.95 << " " <<  0.15
// 	     << endl << " " <<  6.0 << " " << 1.15 << " " << 0.747 << " " << -0.086 << " " << 0.2827025 << " " <<  -3.35 << " " <<  0.10
// 	     << endl << " " <<  6.0 << " " << 1.20 << " " << 0.763 << " " << -0.085 << " " << 0.2946970 << " " <<  -2.75 << " " <<  0.05
// 	     << endl << " " <<  6.0 << " " << 1.30 << " " << 0.794 << " " << -0.083 << " " << 0.3186625 << " " <<  -2.19 << " " <<  0.21
// 	     << endl << " " <<  6.0 << " " << 1.40 << " " << 0.815 << " " << -0.082 << " " << 0.3354745 << " " <<  -1.79 << " " <<  0.13
// 	     << endl << " " <<  6.0 << " " << 1.50 << " " << 0.834 << " " << -0.081 << " " << 0.3510585 << " " <<  -1.53 << " " <<  0.04
// 	     << endl << " " <<  6.0 << " " << 1.70 << " " << 0.869 << " " << -0.077 << " " << 0.3805450 << " " <<  -1.24 << " " <<  0.03
// 	     << endl << " " <<  6.0 << " " << 2.00 << " " << 0.915 << " " << -0.074 << " " << 0.4213505 << " " <<  -0.80 << " " << -0.02
// 	     << endl << " " <<  6.0 << " " << 2.40 << " " << 0.949 << " " << -0.066 << " " << 0.4524785 << " " <<  -0.44 << " " <<  0.02
// 	     << endl << " " <<  6.0 << " " << 2.80 << " " << 0.954 << " " << -0.060 << " " << 0.4568580 << " " <<  -0.33 << " " <<  0.02
// 	     << endl << " " <<  6.0 << " " << 3.20 << " " << 0.956 << " " << -0.054 << " " << 0.4584260 << " " <<  -0.33 << " " <<  0.02
// 	     << endl << " " <<  6.0 << " " << 3.60 << " " << 0.959 << " " << -0.050 << " " << 0.4610905 << " " <<  -0.31 << " " <<  0.02
// 	     << endl << " " <<  6.0 << " " << 4.00 << " " << 0.961 << " " << -0.046 << " " << 0.4628185 << " " <<  -0.29 << " " <<  0.03
// 	     << endl << " " <<  6.0 << " " << 5.00 << " " << 0.963 << " " << -0.040 << " " << 0.4644845 << " " <<  -0.30 << " " <<  0.02
// 	     << endl << " " <<  6.0 << " " << 6.00 << " " << 0.965 << " " << -0.035 << " " << 0.4662250 << " " <<  -0.31 << " " <<  0.01
// 	     << endl << " " <<  6.0 << " " << 7.50 << " " << 0.969 << " " << -0.027 << " " << 0.4698450 << " " <<  -0.34 << " " <<  0.01
// 	     << endl << " " <<  6.0 << " " << 8.20 << " " << 0.870 << " " << -0.018 << " " << 0.3786120 << " " <<  -1.84 << " " << -0.31;
//     Out_File << endl << "ZONE T = \"Station " << 14 << "\" \\"
// 	     << endl << "I = " << 28 << " \\"
// 	     << endl << "J = " << 1 << " \\"
// 	     << endl << " F = POINT \\";
//     Out_File << endl << " " <<  6.5 << " " << 0.10 << " " << 0.188 << " " << -0.046 << " " << 0.0187300 << " " << -11.39 << " " << -0.53
// 	     << endl << " " <<  6.5 << " " << 0.15 << " " << 0.198 << " " << -0.044 << " " << 0.0205700 << " " << -12.71 << " " << -0.22
// 	     << endl << " " <<  6.5 << " " << 0.20 << " " << 0.220 << " " << -0.052 << " " << 0.0255520 << " " << -14.18 << " " << -0.38
// 	     << endl << " " <<  6.5 << " " << 0.30 << " " << 0.268 << " " << -0.065 << " " << 0.0380245 << " " << -16.27 << " " << -0.34
// 	     << endl << " " <<  6.5 << " " << 0.40 << " " << 0.321 << " " << -0.074 << " " << 0.0542585 << " " << -16.52 << " " << -0.53
// 	     << endl << " " <<  6.5 << " " << 0.50 << " " << 0.379 << " " << -0.080 << " " << 0.0750205 << " " << -15.60 << " " << -0.36
// 	     << endl << " " <<  6.5 << " " << 0.60 << " " << 0.440 << " " << -0.083 << " " << 0.1002445 << " " << -14.15 << " " << -0.36
// 	     << endl << " " <<  6.5 << " " << 0.70 << " " << 0.504 << " " << -0.083 << " " << 0.1304525 << " " << -12.63 << " " <<  0.32
// 	     << endl << " " <<  6.5 << " " << 0.80 << " " << 0.570 << " " << -0.083 << " " << 0.1658945 << " " << -10.54 << " " << -0.70
// 	     << endl << " " <<  6.5 << " " << 0.90 << " " << 0.630 << " " << -0.081 << " " << 0.2017305 << " " <<  -7.19 << " " <<  0.06
// 	     << endl << " " <<  6.5 << " " << 1.00 << " " << 0.682 << " " << -0.080 << " " << 0.2357620 << " " <<  -5.52 << " " <<  0.43
// 	     << endl << " " <<  6.5 << " " << 1.10 << " " << 0.726 << " " << -0.079 << " " << 0.2666585 << " " <<  -3.75 << " " <<  0.05
// 	     << endl << " " <<  6.5 << " " << 1.15 << " " << 0.743 << " " << -0.078 << " " << 0.2790665 << " " <<  -3.21 << " " <<  0.22
// 	     << endl << " " <<  6.5 << " " << 1.20 << " " << 0.761 << " " << -0.077 << " " << 0.2925250 << " " <<  -2.67 << " " <<  0.39
// 	     << endl << " " <<  6.5 << " " << 1.30 << " " << 0.787 << " " << -0.076 << " " << 0.3125725 << " " <<  -2.04 << " " <<  0.24
// 	     << endl << " " <<  6.5 << " " << 1.40 << " " << 0.811 << " " << -0.075 << " " << 0.3316730 << " " <<  -1.65 << " " <<  0.12
// 	     << endl << " " <<  6.5 << " " << 1.50 << " " << 0.831 << " " << -0.073 << " " << 0.3479450 << " " <<  -1.42 << " " <<  0.08
// 	     << endl << " " <<  6.5 << " " << 1.70 << " " << 0.865 << " " << -0.071 << " " << 0.3766330 << " " <<  -1.09 << " " <<  0.07
// 	     << endl << " " <<  6.5 << " " << 2.00 << " " << 0.908 << " " << -0.067 << " " << 0.4144765 << " " <<  -0.67 << " " <<  0.02
// 	     << endl << " " <<  6.5 << " " << 2.40 << " " << 0.941 << " " << -0.061 << " " << 0.4446010 << " " <<  -0.38 << " " <<  0.03
// 	     << endl << " " <<  6.5 << " " << 2.80 << " " << 0.949 << " " << -0.056 << " " << 0.4518685 << " " <<  -0.27 << " " <<  0.03
// 	     << endl << " " <<  6.5 << " " << 3.20 << " " << 0.952 << " " << -0.051 << " " << 0.4544525 << " " <<  -0.24 << " " <<  0.02
// 	     << endl << " " <<  6.5 << " " << 3.60 << " " << 0.955 << " " << -0.048 << " " << 0.4571645 << " " <<  -0.24 << " " <<  0.03
// 	     << endl << " " <<  6.5 << " " << 4.00 << " " << 0.956 << " " << -0.045 << " " << 0.4579805 << " " <<  -0.23 << " " <<  0.03
// 	     << endl << " " <<  6.5 << " " << 5.00 << " " << 0.960 << " " << -0.039 << " " << 0.4615605 << " " <<  -0.23 << " " <<  0.02
// 	     << endl << " " <<  6.5 << " " << 6.00 << " " << 0.963 << " " << -0.033 << " " << 0.4642290 << " " <<  -0.22 << " " <<  0.02
// 	     << endl << " " <<  6.5 << " " << 7.50 << " " << 0.963 << " " << -0.027 << " " << 0.4640490 << " " <<  -0.31 << " " <<  0.02
// 	     << endl << " " <<  6.5 << " " << 8.20 << " " << 0.864 << " " << -0.017 << " " << 0.3733925 << " " <<  -1.83 << " " << -0.42;
//     Out_File << endl << "ZONE T = \"Station " << 15 << "\" \\"
// 	     << endl << "I = " << 27 << " \\"
// 	     << endl << "J = " << 1 << " \\"
// 	     << endl << " F = POINT \\";
//     Out_File << endl << " " <<  7.0 << " " << 0.10 << " " << 0.251 << " " << -0.040 << " " << 0.0323005 << " " <<  -9.95 << " " << -0.60
// 	     << endl << " " <<  7.0 << " " << 0.15 << " " << 0.253 << " " << -0.042 << " " << 0.0328865 << " " << -11.21 << " " << -0.69
// 	     << endl << " " <<  7.0 << " " << 0.20 << " " << 0.274 << " " << -0.050 << " " << 0.0387880 << " " << -12.69 << " " << -1.28
// 	     << endl << " " <<  7.0 << " " << 0.30 << " " << 0.315 << " " << -0.061 << " " << 0.0514730 << " " << -14.39 << " " << -0.19
// 	     << endl << " " <<  7.0 << " " << 0.40 << " " << 0.362 << " " << -0.069 << " " << 0.0679025 << " " << -14.58 << " " << -0.11
// 	     << endl << " " <<  7.0 << " " << 0.50 << " " << 0.418 << " " << -0.075 << " " << 0.0901745 << " " << -13.64 << " " <<  0.30
// 	     << endl << " " <<  7.0 << " " << 0.60 << " " << 0.473 << " " << -0.079 << " " << 0.1149850 << " " << -12.45 << " " <<  0.17
// 	     << endl << " " <<  7.0 << " " << 0.70 << " " << 0.529 << " " << -0.079 << " " << 0.1430410 << " " << -10.77 << " " << -0.12
// 	     << endl << " " <<  7.0 << " " << 0.80 << " " << 0.585 << " " << -0.078 << " " << 0.1741545 << " " <<  -9.36 << " " << -0.74
// 	     << endl << " " <<  7.0 << " " << 0.90 << " " << 0.639 << " " << -0.076 << " " << 0.2070485 << " " <<  -7.27 << " " << -0.50
// 	     << endl << " " <<  7.0 << " " << 1.00 << " " << 0.688 << " " << -0.073 << " " << 0.2393365 << " " <<  -5.49 << " " << -0.09
// 	     << endl << " " <<  7.0 << " " << 1.10 << " " << 0.729 << " " << -0.074 << " " << 0.2684585 << " " <<  -4.06 << " " << -0.08
// 	     << endl << " " <<  7.0 << " " << 1.15 << " " << 0.744 << " " << -0.072 << " " << 0.2793600 << " " <<  -3.58 << " " <<  0.08
// 	     << endl << " " <<  7.0 << " " << 1.20 << " " << 0.760 << " " << -0.070 << " " << 0.2912500 << " " <<  -3.10 << " " <<  0.24
// 	     << endl << " " <<  7.0 << " " << 1.30 << " " << 0.788 << " " << -0.071 << " " << 0.3129925 << " " <<  -2.28 << " " <<  0.21
// 	     << endl << " " <<  7.0 << " " << 1.40 << " " << 0.809 << " " << -0.068 << " " << 0.3295525 << " " <<  -1.99 << " " <<  0.22
// 	     << endl << " " <<  7.0 << " " << 1.50 << " " << 0.831 << " " << -0.068 << " " << 0.3475925 << " " <<  -1.65 << " " <<  0.16
// 	     << endl << " " <<  7.0 << " " << 1.70 << " " << 0.867 << " " << -0.068 << " " << 0.3781565 << " " <<  -1.21 << " " <<  0.07
// 	     << endl << " " <<  7.0 << " " << 2.00 << " " << 0.910 << " " << -0.065 << " " << 0.4161625 << " " <<  -0.76 << " " <<  0.01
// 	     << endl << " " <<  7.0 << " " << 2.40 << " " << 0.942 << " " << -0.059 << " " << 0.4454225 << " " <<  -0.45 << " " <<  0.02
// 	     << endl << " " <<  7.0 << " " << 2.80 << " " << 0.950 << " " << -0.055 << " " << 0.4527625 << " " <<  -0.35 << " " <<  0.03
// 	     << endl << " " <<  7.0 << " " << 3.20 << " " << 0.953 << " " << -0.051 << " " << 0.4554050 << " " <<  -0.32 << " " <<  0.02
// 	     << endl << " " <<  7.0 << " " << 3.60 << " " << 0.956 << " " << -0.048 << " " << 0.4581200 << " " <<  -0.31 << " " <<  0.02
// 	     << endl << " " <<  7.0 << " " << 4.00 << " " << 0.958 << " " << -0.045 << " " << 0.4598945 << " " <<  -0.32 << " " <<  0.06
// 	     << endl << " " <<  7.0 << " " << 5.00 << " " << 0.962 << " " << -0.039 << " " << 0.4634825 << " " <<  -0.31 << " " <<  0.04
// 	     << endl << " " <<  7.0 << " " << 6.00 << " " << 0.965 << " " << -0.035 << " " << 0.4662250 << " " <<  -0.34 << " " <<  0.06
// 	     << endl << " " <<  7.0 << " " << 7.50 << " " << 0.966 << " " << -0.028 << " " << 0.4669700 << " " <<  -0.38 << " " <<  0.02;
//     Out_File << endl << "ZONE T = \"Station " << 16 << "\" \\"
// 	     << endl << "I = " << 28 << " \\"
// 	     << endl << "J = " << 1 << " \\"
// 	     << endl << " F = POINT \\";
//     Out_File << endl << " " <<  8.0 << " " << 0.10 << " " << 0.313 << " " << -0.037 << " " << 0.0496690 << " " <<  -8.23 << " " << -0.03
// 	     << endl << " " <<  8.0 << " " << 0.15 << " " << 0.327 << " " << -0.039 << " " << 0.0542250 << " " <<  -9.26 << " " << -0.92
// 	     << endl << " " <<  8.0 << " " << 0.20 << " " << 0.340 << " " << -0.044 << " " << 0.0587680 << " " << -10.46 << " " << -0.51
// 	     << endl << " " <<  8.0 << " " << 0.30 << " " << 0.377 << " " << -0.053 << " " << 0.0724690 << " " << -12.21 << " " <<  0.53
// 	     << endl << " " <<  8.0 << " " << 0.40 << " " << 0.415 << " " << -0.062 << " " << 0.0880345 << " " << -12.53 << " " << -0.28
// 	     << endl << " " <<  8.0 << " " << 0.50 << " " << 0.459 << " " << -0.068 << " " << 0.1076525 << " " << -12.25 << " " << -0.06
// 	     << endl << " " <<  8.0 << " " << 0.60 << " " << 0.511 << " " << -0.071 << " " << 0.1330810 << " " << -11.05 << " " << -0.42
// 	     << endl << " " <<  8.0 << " " << 0.70 << " " << 0.558 << " " << -0.072 << " " << 0.1582740 << " " <<  -9.46 << " " << -0.34
// 	     << endl << " " <<  8.0 << " " << 0.80 << " " << 0.609 << " " << -0.068 << " " << 0.1877525 << " " <<  -8.34 << " " << -0.40
// 	     << endl << " " <<  8.0 << " " << 0.90 << " " << 0.657 << " " << -0.067 << " " << 0.2180690 << " " <<  -6.60 << " " << -0.77
// 	     << endl << " " <<  8.0 << " " << 1.00 << " " << 0.695 << " " << -0.064 << " " << 0.2435605 << " " <<  -5.19 << " " << -0.40
// 	     << endl << " " <<  8.0 << " " << 1.10 << " " << 0.731 << " " << -0.062 << " " << 0.2691025 << " " <<  -4.06 << " " << -0.14
// 	     << endl << " " <<  8.0 << " " << 1.15 << " " << 0.748 << " " << -0.062 << " " << 0.2816740 << " " <<  -3.48 << " " <<  0.09
// 	     << endl << " " <<  8.0 << " " << 1.20 << " " << 0.764 << " " << -0.063 << " " << 0.2938325 << " " <<  -2.89 << " " <<  0.32
// 	     << endl << " " <<  8.0 << " " << 1.30 << " " << 0.791 << " " << -0.062 << " " << 0.3147625 << " " <<  -2.17 << " " <<  0.22
// 	     << endl << " " <<  8.0 << " " << 1.40 << " " << 0.813 << " " << -0.060 << " " << 0.3322845 << " " <<  -1.78 << " " <<  0.18
// 	     << endl << " " <<  8.0 << " " << 1.50 << " " << 0.832 << " " << -0.060 << " " << 0.3479120 << " " <<  -1.49 << " " <<  0.12
// 	     << endl << " " <<  8.0 << " " << 1.70 << " " << 0.868 << " " << -0.058 << " " << 0.3783940 << " " <<  -1.07 << " " <<  0.03
// 	     << endl << " " <<  8.0 << " " << 2.00 << " " << 0.909 << " " << -0.057 << " " << 0.4147650 << " " <<  -0.65 << " " <<  0.04
// 	     << endl << " " <<  8.0 << " " << 2.40 << " " << 0.935 << " " << -0.052 << " " << 0.4384645 << " " <<  -0.37 << " " <<  0.04
// 	     << endl << " " <<  8.0 << " " << 2.80 << " " << 0.943 << " " << -0.048 << " " << 0.4457765 << " " <<  -0.27 << " " <<  0.02
// 	     << endl << " " <<  8.0 << " " << 3.20 << " " << 0.946 << " " << -0.047 << " " << 0.4485625 << " " <<  -0.25 << " " <<  0.04
// 	     << endl << " " <<  8.0 << " " << 3.60 << " " << 0.949 << " " << -0.043 << " " << 0.4512250 << " " <<  -0.23 << " " <<  0.03
// 	     << endl << " " <<  8.0 << " " << 4.00 << " " << 0.950 << " " << -0.042 << " " << 0.4521320 << " " <<  -0.22 << " " <<  0.04
// 	     << endl << " " <<  8.0 << " " << 5.00 << " " << 0.954 << " " << -0.038 << " " << 0.4557800 << " " <<  -0.23 << " " <<  0.04
// 	     << endl << " " <<  8.0 << " " << 6.00 << " " << 0.958 << " " << -0.034 << " " << 0.4594600 << " " <<  -0.24 << " " <<  0.04
// 	     << endl << " " <<  8.0 << " " << 7.50 << " " << 0.959 << " " << -0.028 << " " << 0.4602325 << " " <<  -0.31 << " " <<  0.02
// 	     << endl << " " <<  8.0 << " " << 8.20 << " " << 0.855 << " " << -0.019 << " " << 0.3656930 << " " <<  -1.97 << " " << -0.41;
//     Out_File << endl << "ZONE T = \"Station " << 17 << "\" \\"
// 	     << endl << "I = " << 28 << " \\"
// 	     << endl << "J = " << 1 << " \\"
// 	     << endl << " F = POINT \\";
//     Out_File << endl << " " << 10.0 << " " << 0.10 << " " << 0.418 << " " << -0.029 << " " << 0.0877825 << " " << -5.63 << " " <<  0.27
// 	     << endl << " " << 10.0 << " " << 0.15 << " " << 0.433 << " " << -0.027 << " " << 0.0941090 << " " << -6.15 << " " << -0.36
// 	     << endl << " " << 10.0 << " " << 0.20 << " " << 0.447 << " " << -0.032 << " " << 0.1004165 << " " << -6.97 << " " << -0.16
// 	     << endl << " " << 10.0 << " " << 0.30 << " " << 0.471 << " " << -0.036 << " " << 0.1115685 << " " << -8.34 << " " <<  0.18
// 	     << endl << " " << 10.0 << " " << 0.40 << " " << 0.496 << " " << -0.041 << " " << 0.1238485 << " " << -9.08 << " " <<  0.09
// 	     << endl << " " << 10.0 << " " << 0.50 << " " << 0.526 << " " << -0.043 << " " << 0.1392625 << " " << -9.16 << " " <<  0.18
// 	     << endl << " " << 10.0 << " " << 0.60 << " " << 0.560 << " " << -0.045 << " " << 0.1578125 << " " << -9.05 << " " << -0.16
// 	     << endl << " " << 10.0 << " " << 0.70 << " " << 0.596 << " " << -0.047 << " " << 0.1787125 << " " << -8.44 << " " << -0.18
// 	     << endl << " " << 10.0 << " " << 0.80 << " " << 0.638 << " " << -0.047 << " " << 0.2046265 << " " << -7.75 << " " << -0.07
// 	     << endl << " " << 10.0 << " " << 0.90 << " " << 0.676 << " " << -0.046 << " " << 0.2295460 << " " << -6.44 << " " << -0.44
// 	     << endl << " " << 10.0 << " " << 1.00 << " " << 0.712 << " " << -0.047 << " " << 0.2545765 << " " << -5.35 << " " << -0.57
// 	     << endl << " " << 10.0 << " " << 1.10 << " " << 0.745 << " " << -0.046 << " " << 0.2785705 << " " << -4.32 << " " << -0.21
// 	     << endl << " " << 10.0 << " " << 1.15 << " " << 0.758 << " " << -0.045 << " " << 0.2882945 << " " << -3.92 << " " << -0.26
// 	     << endl << " " << 10.0 << " " << 1.20 << " " << 0.772 << " " << -0.044 << " " << 0.2989600 << " " << -3.52 << " " << -0.32
// 	     << endl << " " << 10.0 << " " << 1.30 << " " << 0.798 << " " << -0.042 << " " << 0.3192840 << " " << -2.80 << " " <<  0.10
// 	     << endl << " " << 10.0 << " " << 1.40 << " " << 0.819 << " " << -0.040 << " " << 0.3361805 << " " << -2.21 << " " <<  0.15
// 	     << endl << " " << 10.0 << " " << 1.50 << " " << 0.839 << " " << -0.039 << " " << 0.3527210 << " " << -1.80 << " " <<  0.12
// 	     << endl << " " << 10.0 << " " << 1.70 << " " << 0.871 << " " << -0.038 << " " << 0.3800425 << " " << -1.03 << " " <<  0.06
// 	     << endl << " " << 10.0 << " " << 2.00 << " " << 0.908 << " " << -0.034 << " " << 0.4128100 << " " << -0.66 << " " <<  0.06
// 	     << endl << " " << 10.0 << " " << 2.40 << " " << 0.929 << " " << -0.032 << " " << 0.4320325 << " " << -0.35 << " " <<  0.05
// 	     << endl << " " << 10.0 << " " << 2.80 << " " << 0.933 << " " << -0.031 << " " << 0.4357250 << " " << -0.27 << " " <<  0.04
// 	     << endl << " " << 10.0 << " " << 3.20 << " " << 0.937 << " " << -0.032 << " " << 0.4394965 << " " << -0.25 << " " <<  0.06
// 	     << endl << " " << 10.0 << " " << 3.60 << " " << 0.939 << " " << -0.030 << " " << 0.4413105 << " " << -0.24 << " " <<  0.05
// 	     << endl << " " << 10.0 << " " << 4.00 << " " << 0.940 << " " << -0.031 << " " << 0.4422805 << " " << -0.22 << " " <<  0.03
// 	     << endl << " " << 10.0 << " " << 5.00 << " " << 0.944 << " " << -0.029 << " " << 0.4459885 << " " << -0.22 << " " <<  0.03
// 	     << endl << " " << 10.0 << " " << 6.00 << " " << 0.947 << " " << -0.028 << " " << 0.4487965 << " " << -0.23 << " " <<  0.04
// 	     << endl << " " << 10.0 << " " << 7.50 << " " << 0.949 << " " << -0.024 << " " << 0.4505885 << " " << -0.32 << " " <<  0.02
// 	     << endl << " " << 10.0 << " " << 8.20 << " " << 0.827 << " " << -0.019 << " " << 0.3421450 << " " << -2.18 << " " << -0.41;
//     Out_File << endl << "ZONE T = \"Station " << 18 << "\" \\"
// 	     << endl << "I = " << 28 << " \\"
// 	     << endl << "J = " << 1 << " \\"
// 	     << endl << " F = POINT \\";
//     Out_File << endl << " " << 12.0 << " " << 0.10 << " " << 0.461 << " " << -0.015 << " " << 0.1063730 << " " << -5.97 << " " <<  0.24
// 	     << endl << " " << 12.0 << " " << 0.15 << " " << 0.475 << " " << -0.014 << " " << 0.1129105 << " " << -5.38 << " " <<  0.34
// 	     << endl << " " << 12.0 << " " << 0.20 << " " << 0.488 << " " << -0.015 << " " << 0.1191845 << " " << -5.93 << " " <<  0.26
// 	     << endl << " " << 12.0 << " " << 0.30 << " " << 0.511 << " " << -0.019 << " " << 0.1307410 << " " << -7.45 << " " <<  0.28
// 	     << endl << " " << 12.0 << " " << 0.40 << " " << 0.537 << " " << -0.021 << " " << 0.1444050 << " " << -8.20 << " " <<  0.21
// 	     << endl << " " << 12.0 << " " << 0.50 << " " << 0.563 << " " << -0.027 << " " << 0.1588490 << " " << -8.70 << " " <<  0.21
// 	     << endl << " " << 12.0 << " " << 0.60 << " " << 0.594 << " " << -0.027 << " " << 0.1767825 << " " << -9.11 << " " <<  0.33
// 	     << endl << " " << 12.0 << " " << 0.70 << " " << 0.625 << " " << -0.028 << " " << 0.1957045 << " " << -8.46 << " " << -0.40
// 	     << endl << " " << 12.0 << " " << 0.80 << " " << 0.657 << " " << -0.027 << " " << 0.2161890 << " " << -7.94 << " " << -0.31
// 	     << endl << " " << 12.0 << " " << 0.90 << " " << 0.691 << " " << -0.027 << " " << 0.2391050 << " " << -7.82 << " " << -0.88
// 	     << endl << " " << 12.0 << " " << 1.00 << " " << 0.721 << " " << -0.025 << " " << 0.2602330 << " " << -6.89 << " " << -0.76
// 	     << endl << " " << 12.0 << " " << 1.10 << " " << 0.748 << " " << -0.026 << " " << 0.2800900 << " " << -4.69 << " " << -0.19
// 	     << endl << " " << 12.0 << " " << 1.15 << " " << 0.763 << " " << -0.025 << " " << 0.2913970 << " " << -4.30 << " " << -0.22
// 	     << endl << " " << 12.0 << " " << 1.20 << " " << 0.779 << " " << -0.025 << " " << 0.3037330 << " " << -3.91 << " " << -0.26
// 	     << endl << " " << 12.0 << " " << 1.30 << " " << 0.804 << " " << -0.025 << " " << 0.3235205 << " " << -3.19 << " " << -0.17
// 	     << endl << " " << 12.0 << " " << 1.40 << " " << 0.826 << " " << -0.021 << " " << 0.3413585 << " " << -2.69 << " " << -0.02
// 	     << endl << " " << 12.0 << " " << 1.50 << " " << 0.845 << " " << -0.021 << " " << 0.3572330 << " " << -2.10 << " " <<  0.18
// 	     << endl << " " << 12.0 << " " << 1.70 << " " << 0.878 << " " << -0.020 << " " << 0.3856420 << " " << -1.23 << " " <<  0.09
// 	     << endl << " " << 12.0 << " " << 2.00 << " " << 0.913 << " " << -0.018 << " " << 0.4169465 << " " << -0.69 << " " <<  0.06
// 	     << endl << " " << 12.0 << " " << 2.40 << " " << 0.935 << " " << -0.017 << " " << 0.4372570 << " " << -0.39 << " " <<  0.06
// 	     << endl << " " << 12.0 << " " << 2.80 << " " << 0.939 << " " << -0.016 << " " << 0.4409885 << " " << -0.28 << " " <<  0.04
// 	     << endl << " " << 12.0 << " " << 3.20 << " " << 0.941 << " " << -0.017 << " " << 0.4428850 << " " << -0.25 << " " <<  0.04
// 	     << endl << " " << 12.0 << " " << 3.60 << " " << 0.943 << " " << -0.017 << " " << 0.4447690 << " " << -0.24 << " " <<  0.04
// 	     << endl << " " << 12.0 << " " << 4.00 << " " << 0.944 << " " << -0.018 << " " << 0.4457300 << " " << -0.23 << " " <<  0.04
// 	     << endl << " " << 12.0 << " " << 5.00 << " " << 0.947 << " " << -0.018 << " " << 0.4485665 << " " << -0.24 << " " <<  0.03
// 	     << endl << " " << 12.0 << " " << 6.00 << " " << 0.949 << " " << -0.018 << " " << 0.4504625 << " " << -0.24 << " " <<  0.03
// 	     << endl << " " << 12.0 << " " << 7.50 << " " << 0.947 << " " << -0.015 << " " << 0.4485170 << " " << -0.36 << " " <<  0.00
// 	     << endl << " " << 12.0 << " " << 8.20 << " " << 0.825 << " " << -0.007 << " " << 0.3403370 << " " << -2.11 << " " << -0.39;
//     Out_File << endl << "ZONE T = \"Station " << 19 << "\" \\"
// 	     << endl << "I = " << 28 << " \\"
// 	     << endl << "J = " << 1 << " \\"
// 	     << endl << " F = POINT \\";
//     Out_File << endl << " " << 14.0 << " " << 0.10 << " " << 0.505 << " " << -0.014 << " " << 0.1276105 << " " << -4.12 << " " <<  0.19
// 	     << endl << " " << 14.0 << " " << 0.15 << " " << 0.511 << " " << -0.015 << " " << 0.1306730 << " " << -5.67 << " " <<  0.16
// 	     << endl << " " << 14.0 << " " << 0.20 << " " << 0.529 << " " << -0.013 << " " << 0.1400050 << " " << -5.19 << " " <<  0.19
// 	     << endl << " " << 14.0 << " " << 0.30 << " " << 0.553 << " " << -0.014 << " " << 0.1530025 << " " << -6.13 << " " <<  0.00
// 	     << endl << " " << 14.0 << " " << 0.40 << " " << 0.570 << " " << -0.016 << " " << 0.1625780 << " " << -6.50 << " " <<  0.26
// 	     << endl << " " << 14.0 << " " << 0.50 << " " << 0.593 << " " << -0.021 << " " << 0.1760450 << " " << -6.79 << " " <<  0.29
// 	     << endl << " " << 14.0 << " " << 0.60 << " " << 0.618 << " " << -0.020 << " " << 0.1911620 << " " << -6.82 << " " <<  0.09
// 	     << endl << " " << 14.0 << " " << 0.70 << " " << 0.643 << " " << -0.021 << " " << 0.2069450 << " " << -6.90 << " " <<  0.10
// 	     << endl << " " << 14.0 << " " << 0.80 << " " << 0.671 << " " << -0.023 << " " << 0.2253850 << " " << -6.50 << " " <<  0.11
// 	     << endl << " " << 14.0 << " " << 0.90 << " " << 0.700 << " " << -0.023 << " " << 0.2452645 << " " << -5.85 << " " <<  0.07
// 	     << endl << " " << 14.0 << " " << 1.00 << " " << 0.726 << " " << -0.022 << " " << 0.2637800 << " " << -5.47 << " " <<  0.16
// 	     << endl << " " << 14.0 << " " << 1.10 << " " << 0.754 << " " << -0.021 << " " << 0.2844785 << " " << -4.81 << " " << -0.22
// 	     << endl << " " << 14.0 << " " << 1.15 << " " << 0.768 << " " << -0.021 << " " << 0.2951325 << " " << -4.37 << " " << -0.22
// 	     << endl << " " << 14.0 << " " << 1.20 << " " << 0.783 << " " << -0.021 << " " << 0.3067650 << " " << -3.93 << " " << -0.22
// 	     << endl << " " << 14.0 << " " << 1.30 << " " << 0.809 << " " << -0.021 << " " << 0.3274610 << " " << -3.46 << " " << -0.18
// 	     << endl << " " << 14.0 << " " << 1.40 << " " << 0.831 << " " << -0.019 << " " << 0.3454610 << " " << -2.92 << " " << -0.08
// 	     << endl << " " << 14.0 << " " << 1.50 << " " << 0.851 << " " << -0.018 << " " << 0.3622625 << " " << -2.28 << " " <<  0.05
// 	     << endl << " " << 14.0 << " " << 1.70 << " " << 0.885 << " " << -0.017 << " " << 0.3917570 << " " << -1.43 << " " <<  0.06
// 	     << endl << " " << 14.0 << " " << 2.00 << " " << 0.919 << " " << -0.014 << " " << 0.4223785 << " " << -0.74 << " " <<  0.04
// 	     << endl << " " << 14.0 << " " << 2.40 << " " << 0.938 << " " << -0.012 << " " << 0.4399940 << " " << -0.41 << " " <<  0.04
// 	     << endl << " " << 14.0 << " " << 2.80 << " " << 0.941 << " " << -0.011 << " " << 0.4428010 << " " << -0.31 << " " <<  0.04
// 	     << endl << " " << 14.0 << " " << 3.20 << " " << 0.943 << " " << -0.012 << " " << 0.4446965 << " " << -0.27 << " " <<  0.04
// 	     << endl << " " << 14.0 << " " << 3.60 << " " << 0.944 << " " << -0.012 << " " << 0.4456400 << " " << -0.25 << " " <<  0.03
// 	     << endl << " " << 14.0 << " " << 4.00 << " " << 0.945 << " " << -0.014 << " " << 0.4466105 << " " << -0.26 << " " <<  0.03
// 	     << endl << " " << 14.0 << " " << 5.00 << " " << 0.948 << " " << -0.013 << " " << 0.4494365 << " " << -0.27 << " " <<  0.03
// 	     << endl << " " << 14.0 << " " << 6.00 << " " << 0.951 << " " << -0.012 << " " << 0.4522725 << " " << -0.28 << " " <<  0.04
// 	     << endl << " " << 14.0 << " " << 7.50 << " " << 0.946 << " " << -0.012 << " " << 0.4475300 << " " << -0.49 << " " <<  0.00
// 	     << endl << " " << 14.0 << " " << 8.20 << " " << 0.824 << " " << -0.005 << " " << 0.3395005 << " " << -2.33 << " " << -0.33;
//     Out_File << endl << "ZONE T = \"Station " << 20 << "\" \\"
// 	     << endl << "I = " << 28 << " \\"
// 	     << endl << "J = " << 1 << " \\"
// 	     << endl << " F = POINT \\";
//     Out_File << endl << " " << 16.0 << " " << 0.10 << " " << 0.517 << " " << -0.011 << " " << 0.1337050 << " " << -3.74 << " " <<  0.27
// 	     << endl << " " << 16.0 << " " << 0.15 << " " << 0.538 << " " << -0.011 << " " << 0.1447825 << " " << -3.74 << " " <<  0.22
// 	     << endl << " " << 16.0 << " " << 0.20 << " " << 0.547 << " " << -0.010 << " " << 0.1496545 << " " << -4.27 << " " << -0.09
// 	     << endl << " " << 16.0 << " " << 0.30 << " " << 0.570 << " " << -0.012 << " " << 0.1625220 << " " << -4.98 << " " << -0.04
// 	     << endl << " " << 16.0 << " " << 0.40 << " " << 0.590 << " " << -0.015 << " " << 0.1741625 << " " << -5.52 << " " << -0.03
// 	     << endl << " " << 16.0 << " " << 0.50 << " " << 0.610 << " " << -0.016 << " " << 0.1861780 << " " << -5.93 << " " <<  0.14
// 	     << endl << " " << 16.0 << " " << 0.60 << " " << 0.631 << " " << -0.016 << " " << 0.1992085 << " " << -6.24 << " " <<  0.21
// 	     << endl << " " << 16.0 << " " << 0.70 << " " << 0.654 << " " << -0.019 << " " << 0.2140385 << " " << -5.96 << " " <<  0.22
// 	     << endl << " " << 16.0 << " " << 0.80 << " " << 0.674 << " " << -0.017 << " " << 0.2272825 << " " << -6.05 << " " << -0.08
// 	     << endl << " " << 16.0 << " " << 0.90 << " " << 0.699 << " " << -0.018 << " " << 0.2444625 << " " << -5.68 << " " <<  0.08
// 	     << endl << " " << 16.0 << " " << 1.00 << " " << 0.722 << " " << -0.019 << " " << 0.2608225 << " " << -5.17 << " " << -0.25
// 	     << endl << " " << 16.0 << " " << 1.10 << " " << 0.747 << " " << -0.018 << " " << 0.2791665 << " " << -4.66 << " " << -0.32
// 	     << endl << " " << 16.0 << " " << 1.15 << " " << 0.758 << " " << -0.017 << " " << 0.2874265 << " " << -4.46 << " " << -0.28
// 	     << endl << " " << 16.0 << " " << 1.20 << " " << 0.769 << " " << -0.016 << " " << 0.2958085 << " " << -4.27 << " " << -0.25
// 	     << endl << " " << 16.0 << " " << 1.30 << " " << 0.796 << " " << -0.017 << " " << 0.3169525 << " " << -3.62 << " " << -0.32
// 	     << endl << " " << 16.0 << " " << 1.40 << " " << 0.817 << " " << -0.015 << " " << 0.3338570 << " " << -3.20 << " " << -0.20
// 	     << endl << " " << 16.0 << " " << 1.50 << " " << 0.838 << " " << -0.015 << " " << 0.3512345 << " " << -2.58 << " " << -0.17
// 	     << endl << " " << 16.0 << " " << 1.70 << " " << 0.874 << " " << -0.013 << " " << 0.3820225 << " " << -1.65 << " " <<  0.05
// 	     << endl << " " << 16.0 << " " << 2.00 << " " << 0.911 << " " << -0.011 << " " << 0.4150210 << " " << -0.87 << " " <<  0.05
// 	     << endl << " " << 16.0 << " " << 2.40 << " " << 0.932 << " " << -0.008 << " " << 0.4343440 << " " << -0.47 << " " <<  0.07
// 	     << endl << " " << 16.0 << " " << 2.80 << " " << 0.935 << " " << -0.007 << " " << 0.4371370 << " " << -0.34 << " " <<  0.04
// 	     << endl << " " << 16.0 << " " << 3.20 << " " << 0.936 << " " << -0.008 << " " << 0.4380800 << " " << -0.28 << " " <<  0.04
// 	     << endl << " " << 16.0 << " " << 3.60 << " " << 0.937 << " " << -0.008 << " " << 0.4390165 << " " << -0.25 << " " <<  0.04
// 	     << endl << " " << 16.0 << " " << 4.00 << " " << 0.938 << " " << -0.010 << " " << 0.4399720 << " " << -0.23 << " " <<  0.03
// 	     << endl << " " << 16.0 << " " << 5.00 << " " << 0.939 << " " << -0.011 << " " << 0.4409210 << " " << -0.23 << " " <<  0.03
// 	     << endl << " " << 16.0 << " " << 6.00 << " " << 0.941 << " " << -0.011 << " " << 0.4428010 << " " << -0.22 << " " <<  0.04
// 	     << endl << " " << 16.0 << " " << 7.50 << " " << 0.933 << " " << -0.011 << " " << 0.4353050 << " " << -0.41 << " " << -0.03
// 	     << endl << " " << 16.0 << " " << 8.20 << " " << 0.798 << " " << -0.006 << " " << 0.3184200 << " " << -2.12 << " " << -0.44;
//     Out_File << endl << "ZONE T = \"Station " << 21 << "\" \\"
// 	     << endl << "I = " << 28 << " \\"
// 	     << endl << "J = " << 1 << " \\"
// 	     << endl << " F = POINT \\";
//     Out_File << endl << " " << 20.0 << " " << 0.10 << " " << 0.540 << " " << -0.012 << " " << 0.1458720 << " " << -3.05 << " " <<  0.15
// 	     << endl << " " << 20.0 << " " << 0.15 << " " << 0.566 << " " << -0.009 << " " << 0.1602185 << " " << -2.88 << " " << -0.05
// 	     << endl << " " << 20.0 << " " << 0.20 << " " << 0.582 << " " << -0.011 << " " << 0.1694225 << " " << -3.48 << " " << -0.07
// 	     << endl << " " << 20.0 << " " << 0.30 << " " << 0.604 << " " << -0.012 << " " << 0.1824800 << " " << -4.06 << " " <<  0.02
// 	     << endl << " " << 20.0 << " " << 0.40 << " " << 0.622 << " " << -0.012 << " " << 0.1935140 << " " << -4.58 << " " <<  0.18
// 	     << endl << " " << 20.0 << " " << 0.50 << " " << 0.638 << " " << -0.015 << " " << 0.2036345 << " " << -4.86 << " " <<  0.00
// 	     << endl << " " << 20.0 << " " << 0.60 << " " << 0.658 << " " << -0.015 << " " << 0.2165945 << " " << -5.12 << " " <<  0.13
// 	     << endl << " " << 20.0 << " " << 0.70 << " " << 0.675 << " " << -0.017 << " " << 0.2279570 << " " << -5.04 << " " << -0.07
// 	     << endl << " " << 20.0 << " " << 0.80 << " " << 0.693 << " " << -0.017 << " " << 0.2402690 << " " << -5.10 << " " <<  0.08
// 	     << endl << " " << 20.0 << " " << 0.90 << " " << 0.714 << " " << -0.017 << " " << 0.2550425 << " " << -5.14 << " " <<  0.08
// 	     << endl << " " << 20.0 << " " << 1.00 << " " << 0.736 << " " << -0.020 << " " << 0.2710480 << " " << -4.75 << " " << -0.02
// 	     << endl << " " << 20.0 << " " << 1.10 << " " << 0.754 << " " << -0.017 << " " << 0.2844025 << " " << -4.51 << " " << -0.03
// 	     << endl << " " << 20.0 << " " << 1.15 << " " << 0.764 << " " << -0.018 << " " << 0.2920100 << " " << -4.38 << " " << -0.07
// 	     << endl << " " << 20.0 << " " << 1.20 << " " << 0.773 << " " << -0.019 << " " << 0.2989450 << " " << -4.25 << " " << -0.11
// 	     << endl << " " << 20.0 << " " << 1.30 << " " << 0.796 << " " << -0.018 << " " << 0.3169700 << " " << -3.77 << " " << -0.05
// 	     << endl << " " << 20.0 << " " << 1.40 << " " << 0.817 << " " << -0.018 << " " << 0.3339065 << " " << -3.28 << " " << -0.06
// 	     << endl << " " << 20.0 << " " << 1.50 << " " << 0.833 << " " << -0.017 << " " << 0.3470890 << " " << -2.97 << " " << -0.19
// 	     << endl << " " << 20.0 << " " << 1.70 << " " << 0.866 << " " << -0.015 << " " << 0.3750905 << " " << -2.17 << " " << -0.06
// 	     << endl << " " << 20.0 << " " << 2.00 << " " << 0.908 << " " << -0.012 << " " << 0.4123040 << " " << -1.27 << " " <<  0.02
// 	     << endl << " " << 20.0 << " " << 2.40 << " " << 0.934 << " " << -0.009 << " " << 0.4362185 << " " << -0.75 << " " <<  0.06
// 	     << endl << " " << 20.0 << " " << 2.80 << " " << 0.939 << " " << -0.007 << " " << 0.4408850 << " " << -0.52 << " " <<  0.06
// 	     << endl << " " << 20.0 << " " << 3.20 << " " << 0.940 << " " << -0.008 << " " << 0.4418320 << " " << -0.43 << " " <<  0.04
// 	     << endl << " " << 20.0 << " " << 3.60 << " " << 0.940 << " " << -0.008 << " " << 0.4418320 << " " << -0.41 << " " <<  0.03
// 	     << endl << " " << 20.0 << " " << 4.00 << " " << 0.941 << " " << -0.009 << " " << 0.4427810 << " " << -0.38 << " " <<  0.03
// 	     << endl << " " << 20.0 << " " << 5.00 << " " << 0.942 << " " << -0.010 << " " << 0.4437320 << " " << -0.36 << " " <<  0.03
// 	     << endl << " " << 20.0 << " " << 6.00 << " " << 0.943 << " " << -0.011 << " " << 0.4446850 << " " << -0.33 << " " <<  0.03
// 	     << endl << " " << 20.0 << " " << 7.50 << " " << 0.931 << " " << -0.010 << " " << 0.4334305 << " " << -0.55 << " " << -0.04
// 	     << endl << " " << 20.0 << " " << 8.20 << " " << 0.800 << " " << -0.003 << " " << 0.3200045 << " " << -2.29 << " " << -0.42;
//     Out_File << endl << "ZONE T = \"Station " << 22 << "\" \\"
// 	     << endl << "I = " << 27 << " \\"
// 	     << endl << "J = " << 1 << " \\"
// 	     << endl << " F = POINT \\";
//     Out_File << endl << " " << 32.0 << " " << 0.10 << " " << 0.581 << " " << -0.005 << " " << 0.1687930 << " " << -2.95 << " " <<  0.20
// 	     << endl << " " << 32.0 << " " << 0.15 << " " << 0.604 << " " << -0.004 << " " << 0.1824160 << " " << -2.86 << " " <<  0.15
// 	     << endl << " " << 32.0 << " " << 0.20 << " " << 0.621 << " " << -0.004 << " " << 0.1928285 << " " << -3.02 << " " <<  0.16
// 	     << endl << " " << 32.0 << " " << 0.30 << " " << 0.645 << " " << -0.005 << " " << 0.2080250 << " " << -3.32 << " " <<  0.15
// 	     << endl << " " << 32.0 << " " << 0.40 << " " << 0.662 << " " << -0.004 << " " << 0.2191300 << " " << -3.35 << " " <<  0.06
// 	     << endl << " " << 32.0 << " " << 0.50 << " " << 0.679 << " " << -0.005 << " " << 0.2305330 << " " << -3.63 << " " <<  0.11
// 	     << endl << " " << 32.0 << " " << 0.60 << " " << 0.691 << " " << -0.006 << " " << 0.2387585 << " " << -3.90 << " " <<  0.10
// 	     << endl << " " << 32.0 << " " << 0.70 << " " << 0.705 << " " << -0.006 << " " << 0.2485305 << " " << -3.85 << " " <<  0.10
// 	     << endl << " " << 32.0 << " " << 0.80 << " " << 0.718 << " " << -0.008 << " " << 0.2577940 << " " << -3.96 << " " <<  0.10
// 	     << endl << " " << 32.0 << " " << 0.90 << " " << 0.734 << " " << -0.009 << " " << 0.2694185 << " " << -4.06 << " " <<  0.14
// 	     << endl << " " << 32.0 << " " << 1.00 << " " << 0.747 << " " << -0.010 << " " << 0.2790545 << " " << -4.15 << " " <<  0.02
// 	     << endl << " " << 32.0 << " " << 1.10 << " " << 0.761 << " " << -0.009 << " " << 0.2896010 << " " << -4.01 << " " <<  0.14
// 	     << endl << " " << 32.0 << " " << 1.15 << " " << 0.769 << " " << -0.011 << " " << 0.2957410 << " " << -3.90 << " " <<  0.12
// 	     << endl << " " << 32.0 << " " << 1.20 << " " << 0.778 << " " << -0.013 << " " << 0.3027265 << " " << -3.80 << " " <<  0.09
// 	     << endl << " " << 32.0 << " " << 1.30 << " " << 0.792 << " " << -0.013 << " " << 0.3137165 << " " << -3.73 << " " << -0.04
// 	     << endl << " " << 32.0 << " " << 1.40 << " " << 0.808 << " " << -0.013 << " " << 0.3265165 << " " << -3.53 << " " << -0.10
// 	     << endl << " " << 32.0 << " " << 1.50 << " " << 0.821 << " " << -0.013 << " " << 0.3371050 << " " << -3.37 << " " << -0.13
// 	     << endl << " " << 32.0 << " " << 1.70 << " " << 0.848 << " " << -0.012 << " " << 0.3596240 << " " << -2.86 << " " << -0.23
// 	     << endl << " " << 32.0 << " " << 2.00 << " " << 0.886 << " " << -0.011 << " " << 0.3925585 << " " << -2.01 << " " << -0.20
// 	     << endl << " " << 32.0 << " " << 2.40 << " " << 0.922 << " " << -0.008 << " " << 0.4250740 << " " << -1.10 << " " << -0.02
// 	     << endl << " " << 32.0 << " " << 2.80 << " " << 0.936 << " " << -0.006 << " " << 0.4380660 << " " << -0.64 << " " <<  0.04
// 	     << endl << " " << 32.0 << " " << 3.20 << " " << 0.941 << " " << -0.003 << " " << 0.4427450 << " " << -0.45 << " " <<  0.07
// 	     << endl << " " << 32.0 << " " << 3.60 << " " << 0.942 << " " << -0.004 << " " << 0.4436900 << " " << -0.37 << " " <<  0.03
// 	     << endl << " " << 32.0 << " " << 4.00 << " " << 0.943 << " " << -0.005 << " " << 0.4446370 << " " << -0.33 << " " <<  0.03
// 	     << endl << " " << 32.0 << " " << 5.00 << " " << 0.944 << " " << -0.006 << " " << 0.4455860 << " " << -0.32 << " " <<  0.04
// 	     << endl << " " << 32.0 << " " << 6.00 << " " << 0.944 << " " << -0.007 << " " << 0.4455925 << " " << -0.30 << " " <<  0.05
// 	     << endl << " " << 32.0 << " " << 7.50 << " " << 0.912 << " " << -0.005 << " " << 0.4158845 << " " << -0.87 << " " << -0.18;
//   }

// //   if (fabs(SolnBlk.Grid.nodeSE(SolnBlk.ICu,SolnBlk.JCl).X.x) < TOLER*TOLER) {
// //     Out_File << "ZONE T =  \"u-velocity " << Block_Number
// // 	       << "\" \\ \n"
// // 	       << "I = " << (SolnBlk.JCu+1)-(SolnBlk.JCl-1) << " \\ \n"
// // 	       << "J = " << 1 << " \\ \n"
// // 	       << "F = POINT \\ \n";
// //     for (int j = SolnBlk.JCl; j <= SolnBlk.JCu+1; j++) {
// //       Out_File << SolnBlk.Grid.nodeSE(SolnBlk.ICu,j).X.y/length + HALF << " "
// // 		 << SolnBlk.WnSE(SolnBlk.ICu,j).v.x/Vwall.x << endl;
// //     }
// //   }

}

/**********************************************************************
 * Dusty2DQuadIOSingleBlock.cc                                        *
 *                                                                    *
 * Single-block versions of input and output subroutines for 2D Dusty *
 * multi-block quadrilateral mesh solution classes.                   *
 *                                                                    *
 **********************************************************************/

// Include 2D Dusty quadrilateral mesh solution header file.

#ifndef _DUSTY2D_QUAD_INCLUDED
#include "Dusty2DQuad.h"
#endif // _DUSTY2D_QUAD_INCLUDED

/**********************************************************************
 * Dusty2D_Quad_Block -- IO Single Block External Subroutines.        *
 **********************************************************************/

/**********************************************************************
 * Routine: Write_Solution_Block                                      *
 *                                                                    *
 * Writes the cell centred solution values of the specified           *
 * quadrilateral solution block to the specified output stream for    *
 * restart purposes.                                                  *
 *                                                                    *
 **********************************************************************/
void Write_Solution_Block(Dusty2D_Quad_Block &SolnBlk,
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
void Read_Solution_Block(Dusty2D_Quad_Block &SolnBlk,
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
void Output_Tecplot(Dusty2D_Quad_Block &SolnBlk,
		    Dusty2D_Input_Parameters &IP,
		    const int Number_of_Time_Steps,
		    const double &Time,
		    const int Block_Number,
		    const int Output_Title,
		    ostream &Out_File) {

  Dusty2D_pState W_node;
  Electrostatic2DState E_node;
  Vector2D J_node;
  Particle2D_pState Wp;
  int NUM_CMP_PART = W_node.NUM_CMP_PART;

  // Ensure boundary conditions are updated before evaluating
  // solution at the nodes.
  BCs(SolnBlk,IP);

  // Output node solution data.  
  Out_File << setprecision(14);
  if (Output_Title) {
    Out_File << "TITLE = \"" << CFFC_Name() << ": 2D Dusty Solution, "
	     << "Time Step/Iteration Level = " << Number_of_Time_Steps
	     << ", Time = " << Time
	     << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
	     << "\"y\" \\ \n"
	     << "\"rho\" \\ \n"
	     << "\"vx\" \\ \n"
	     << "\"vy\" \\ \n"
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
	       << "\"muT/mu\" \\ \n";
    }
    if (SolnBlk.Particles == PARTICLE_PHASE_EULERIAN_FORMULATION) {
      Out_File << "\"sigma\" \\ \n"
	       << "\"ux\" \\ \n"
	       << "\"uy\" \\ \n"
	       << "\"Tp\" \\ \n"
	       << "\"chi\" \\ \n"
	       << "\"zetap\" \\ \n"
	       << "\"phip\" \\ \n";
      Out_File << "\"ln(sigmap)\" \\ \n";
      if (NUM_CMP_PART == PARTICLE2D_MULTI_VELOCITY_FORMULATION) {
	for (int npc = 0; npc < NUM_CMP_PART; npc++) {
	  Out_File << "\"sigma" << npc+1 << "\" \\ \n"
		   << "\"ux" << npc+1 << "\" \\ \n"
		   << "\"uy" << npc+1 << "\" \\ \n"
		   << "\"Tp" << npc+1 << "\" \\ \n";
	}
      }
    }
    if (SolnBlk.Electrostatic) {
      Out_File << "\"Ex\" \\ \n"
	       << "\"Ey\" \\ \n"
	       << "\"V\" \\ \n";
      if (SolnBlk.Particles == PARTICLE_PHASE_EULERIAN_FORMULATION) {
	Out_File << "\"Jx\" \\ \n"
		 << "\"Jy\" \\ \n"
		 << "\"J\" \\ \n";
	if (NUM_CMP_PART == PARTICLE2D_MULTI_VELOCITY_FORMULATION) {
	  for (int npc = 0; npc < NUM_CMP_PART; npc++) {
	    Out_File << "\"Jx" << npc+1 << "\" \\ \n"
		     << "\"Jy" << npc+1 << "\" \\ \n"
		     << "\"J"  << npc+1 << "\" \\ \n";
	  }
	}
      }
    }
  }

  Out_File << "ZONE T =  \"Block Number = " << Block_Number
	   << "\" \\ \n"
	   << "I = " << SolnBlk.Grid.INu - SolnBlk.Grid.INl + 1 << " \\ \n"
	   << "J = " << SolnBlk.Grid.JNu - SolnBlk.Grid.JNl + 1 << " \\ \n"
	   << "F = POINT \\ \n";

  for (int j = SolnBlk.Grid.JNl; j <= SolnBlk.Grid.JNu; j++) {
    for (int i = SolnBlk.Grid.INl; i <= SolnBlk.Grid.INu; i++) {
      W_node = SolnBlk.Wn(i,j);
      Out_File.setf(ios::scientific);
      Out_File << " " << SolnBlk.Grid.Node[i][j].X
	       << " " << W_node.rho 
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
		 << " " << W_node.pmodified()
		 << " " << W_node.muT()/W_node.mu();
      }
      if (SolnBlk.Particles == PARTICLE_PHASE_EULERIAN_FORMULATION) {
	if (NUM_CMP_PART == PARTICLE2D_SINGLE_VELOCITY_FORMULATION) {
	  Wp = W_node.Wp[0];
	} else {
	  Wp = Particle2D_pState(W_node.Wp.sigma(),W_node.Wp.u(),W_node.Wp.Tp());
	}
	if (Wp.sigma < 0.0001) Wp.Vacuum();
	Out_File << " " << Wp.sigma
		 << Wp.u
		 << " " << Wp.Tp
		 << " " << W_node.chi()
		 << " " << W_node.zetap()
		 << " " << W_node.phip();
	Out_File << " " << log(Wp.sigma + TOLER);
	if (NUM_CMP_PART == PARTICLE2D_MULTI_VELOCITY_FORMULATION) {
	  for (int npc = 0; npc < NUM_CMP_PART; npc++) {
	    Out_File << W_node.Wp[npc];
	  }
	}
      }
      if (SolnBlk.Electrostatic) {
	E_node = SolnBlk.En(i,j);
	Out_File << " " << E_node.E.x
		 << " " << E_node.E.y
		 << " " << E_node.V;
	if (SolnBlk.Particles == PARTICLE_PHASE_EULERIAN_FORMULATION) {
	  J_node = W_node.Wp.jc(W_node.qe,W_node.mp);
	  Out_File << " " << J_node.x
		   << " " << J_node.y
		   << " " << abs(J_node);
	  if (NUM_CMP_PART == PARTICLE2D_MULTI_VELOCITY_FORMULATION) {
	    for (int npc = 0; npc < NUM_CMP_PART; npc++) {
	      J_node = W_node.Wp[npc].jc(W_node.qe,W_node.mp);
	      Out_File << " " << J_node.x
		       << " " << J_node.y
		       << " " << abs(J_node);
	    }
	  }
	}
      }
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
void Output_Cells_Tecplot(Dusty2D_Quad_Block &SolnBlk,
		          Dusty2D_Input_Parameters &IP,
                          const int Number_of_Time_Steps,
                          const double &Time,
                          const int Block_Number,
                          const int Output_Title,
	                  ostream &Out_File) {

  Particle2D_pState Wp;
  int NUM_CMP_PART = SolnBlk.W[0][0].NUM_CMP_PART;

  Out_File << setprecision(14);
  if (Output_Title) {
    Out_File << "TITLE = \"" << CFFC_Name() << ": 2D Dusty Solution, "
	     << "Time Step/Iteration Level = " << Number_of_Time_Steps
	     << ", Time = " << Time
	     << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
	     << "\"y\" \\ \n"
	     << "\"rho\" \\ \n"
	     << "\"vx\" \\ \n"
	     << "\"vy\" \\ \n"
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
    if (SolnBlk.Particles == PARTICLE_PHASE_EULERIAN_FORMULATION) {
      Out_File << "\"sigma\" \\ \n"
	       << "\"ux\" \\ \n"
	       << "\"uy\" \\ \n"
	       << "\"Tp\" \\ \n"
	       << "\"chi\" \\ \n"
	       << "\"zetap\" \\ \n"
	       << "\"phip\" \\ \n";
      if (NUM_CMP_PART == PARTICLE2D_MULTI_VELOCITY_FORMULATION) {
	for (int npc = 0; npc < NUM_CMP_PART; npc++) {
	  Out_File << "\"sigma" << npc+1 << "\" \\ \n"
		   << "\"ux" << npc+1 << "\" \\ \n"
		   << "\"uy" << npc+1 << "\" \\ \n"
		   << "\"Tp" << npc+1 << "\" \\ \n";
	}
      }
    }
    if (SolnBlk.Electrostatic) {
      Out_File << "\"Ex\" \\ \n"
	       << "\"Ey\" \\ \n"
	       << "\"V\" \\ \n";
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
	   << "I = " << SolnBlk.Grid.ICu - SolnBlk.Grid.ICl + 2*SolnBlk.Nghost + 1 << " \\ \n"
	   << "J = " << SolnBlk.Grid.JCu - SolnBlk.Grid.JCl + 2*SolnBlk.Nghost + 1 << " \\ \n"
	   << "F = POINT \n";

  for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
    for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
      if (SolnBlk.Flow_Type |= FLOWTYPE_INVISCID) {
// 	// Determine the cell-centered gradient.
//  	Linear_Reconstruction_LeastSquares(SolnBlk,i,j,LIMITER_VENKATAKRISHNAN);//IP.i_Limiter);
// 	// Compute the viscous terms.
 	SolnBlk.W[i][j].ComputeViscousTerms(SolnBlk.dWdx[i][j],
 					    SolnBlk.dWdy[i][j],
 					    SolnBlk.Grid.Cell[i][j].Xc,
 					    SolnBlk.Axisymmetric,
 					    OFF);
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
      if (SolnBlk.Particles == PARTICLE_PHASE_EULERIAN_FORMULATION) {
	Wp = Particle2D_pState(SolnBlk.W[i][j].Wp.sigma(),
			       SolnBlk.W[i][j].Wp.u(),
			       SolnBlk.W[i][j].Wp.Tp());
	Out_File << " " << Wp.sigma
		 << Wp.u
		 << " " << Wp.Tp
		 << " " << SolnBlk.W[i][j].chi()
		 << " " << SolnBlk.W[i][j].zetap()
		 << " " << SolnBlk.W[i][j].phip();
	if (NUM_CMP_PART == PARTICLE2D_MULTI_VELOCITY_FORMULATION) {
	  for (int npc = 0; npc < NUM_CMP_PART; npc++) {
	    Out_File << SolnBlk.W[i][j].Wp[npc];
	  }
	}
      }
      Out_File.unsetf(ios::scientific);
      if (SolnBlk.Electrostatic) {
	Out_File.unsetf(ios::scientific);
	Out_File << SolnBlk.E[i][j].E;
	Out_File.setf(ios::scientific);
	Out_File << " " << SolnBlk.E[i][j].V;
	Out_File.unsetf(ios::scientific);
      }
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

//////////////////////////////////////////////////////////////////////

//   double mass = 0.0;
//   double momentum = 0.0;
//   double energy = 0.0;

//   for (int j = 3; j <= 3; j++) {
//     for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
//       mass += SolnBlk.U[i][j].rho + SolnBlk.U[i][j].Up[0].sigma();
//       momentum += SolnBlk.U[i][j].dv.x + SolnBlk.U[i][j].Up[0].du().x;
//       energy += SolnBlk.U[i][j].E + SolnBlk.U[i][j].Up[0].ep;
//     }
//   }

//   cout << endl << " Conservation properties:";
//   cout << endl << " -> Mass = " << mass;
//   cout << endl << " -> Momentum = " << momentum;
//   cout << endl << " -> Energy = " << energy;

//////////////////////////////////////////////////////////////////////

}

/**********************************************************************
 * Routine: Output_Nodes_Tecplot                                      *
 *                                                                    *
 * Writes the node values of the specified quadrilateral solution     *
 * block to the specified output stream suitable for plotting with    *
 * TECPLOT.                                                           *
 *                                                                    *
 **********************************************************************/
void Output_Nodes_Tecplot(Dusty2D_Quad_Block &SolnBlk,
                          const int Number_of_Time_Steps,
                          const double &Time,
			  const int Block_Number,
			  const int Output_Title,
			  ostream &Out_File) {

  // Output node solution data.  
  Out_File << setprecision(14);
  if (Output_Title) {
    Out_File << "TITLE = \"" << CFFC_Name() << ": 2D Dusty Node Values, "
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
void Output_Gradients_Tecplot(Dusty2D_Quad_Block &SolnBlk,
			      const int Number_of_Time_Steps,
			      const double &Time,
			      const int Block_Number,
			      const int Output_Title,
			      ostream &Out_File) {

  int NUM_CMP_PART = SolnBlk.W[0][0].NUM_CMP_PART;

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
    if (SolnBlk.Particles == PARTICLE_PHASE_EULERIAN_FORMULATION) {
      if (NUM_CMP_PART == PARTICLE2D_MULTI_VELOCITY_FORMULATION) {
	for (int npc = 0; npc < NUM_CMP_PART; npc++) {
	  Out_File << "\"dsigmax" << npc+1 << "\" \\ \n"
		   << "\"dsigmay" << npc+1 << "\" \\ \n"
		   << "\"duxx" << npc+1 << "\" \\ \n"
		   << "\"duxy" << npc+1 << "\" \\ \n"
		   << "\"duyx" << npc+1 << "\" \\ \n"
		   << "\"duyy" << npc+1 << "\" \\ \n"
		   << "\"dTpx" << npc+1 << "\" \\ \n"
		   << "\"dTpy" << npc+1 << "\" \\ \n";
	}
      } else {
	Out_File << "\"dsigmax\" \\ \n"
		 << "\"dsigmay\" \\ \n"
		 << "\"duxx\" \\ \n"
		 << "\"duxy\" \\ \n"
		 << "\"duyx\" \\ \n"
		 << "\"duyy\" \\ \n"
		 << "\"dTpx\" \\ \n"
		 << "\"dTpy\" \\ \n";
      }
    }
    Out_File << "\"phi_rho\" \\ \n"
 	     << "\"phi_u\" \\ \n"
 	     << "\"phi_v\" \\ \n"
 	     << "\"phi_p\" \\ \n";
    if (SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
      Out_File << "\"phi_k\" \\ \n"
	       << "\"phi_omega\" \\ \n";
    }
    if (SolnBlk.Particles == PARTICLE_PHASE_EULERIAN_FORMULATION) {
      if (NUM_CMP_PART == PARTICLE2D_MULTI_VELOCITY_FORMULATION) {
	for (int npc = 0; npc < NUM_CMP_PART; npc++) {
	  Out_File << "\"phi_sigma" << npc+1 << "\" \\ \n"
		   << "\"phi_ux" << npc+1 << "\" \\ \n"
		   << "\"phi_uy" << npc+1 << "\" \\ \n"
		   << "\"phi_Tp" << npc+1 << "\" \\ \n";
	}
      } else {
	Out_File << "\"phi_sigma\" \\ \n"
		 << "\"phi_ux\" \\ \n"
		 << "\"phi_uy\" \\ \n"
		 << "\"phi_Tp\" \\ \n";
      }
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
 	       << " " << SolnBlk.dWdy[i][j].p
	       << " " << SolnBlk.phi[i][j].rho
 	       << " " << SolnBlk.phi[i][j].v.x
 	       << " " << SolnBlk.phi[i][j].v.y
 	       << " " << SolnBlk.phi[i][j].p
	       << endl;
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
void Output_Quasi3D_Tecplot(Dusty2D_Quad_Block &SolnBlk,
			    Dusty2D_Input_Parameters &IP,
			    const int Number_of_Time_Steps,
			    const double &Time,
			    const int Block_Number,
			    const int Output_Title,
			    ostream &Out_File) {

  Dusty2D_pState W_node;
  Electrostatic2DState E_node;
  Particle2D_pState Wp;
  Vector2D J_node;
  int NUM_CMP_PART = W_node.NUM_CMP_PART;
  int nrr, numberofrotations = 360/15;
  int numberofnodes = ((SolnBlk.Grid.INu - SolnBlk.Grid.INl + 1)*
		       (SolnBlk.Grid.JNu - SolnBlk.Grid.JNl + 1));

  // Ensure boundary conditions are updated before evaluating
  // solution at the nodes.
  BCs(SolnBlk,IP);

  // Output node solution data.  
  Out_File << setprecision(14);
  if (Output_Title) {
    Out_File << "TITLE = \"" << CFFC_Name() << ": 2D Dusty Solution, "
	     << "Time Step/Iteration Level = " << Number_of_Time_Steps
	     << ", Time = " << Time
	     << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
	     << "\"y\" \\ \n"
	     << "\"z\" \\ \n"
	     << "\"rho\" \\ \n"
	     << "\"vx\" \\ \n"
	     << "\"vy\" \\ \n"
	     << "\"vz\" \\ \n"
	     << "\"p\" \\ \n"
	     << "\"T\" \\ \n"
	     << "\"M\" \\ \n"
	     << "\"H\" \\ \n"
	     << "\"s\" \\ \n";
    if (SolnBlk.Particles == PARTICLE_PHASE_EULERIAN_FORMULATION) {
      Out_File << "\"sigma\" \\ \n"
	       << "\"ux\" \\ \n"
	       << "\"uy\" \\ \n"
	       << "\"uz\" \\ \n"
	       << "\"Tp\" \\ \n";
      if (NUM_CMP_PART == PARTICLE2D_MULTI_VELOCITY_FORMULATION) {
	for (int npc = 0; npc < NUM_CMP_PART; npc++) {
	  Out_File << "\"sigma" << npc+1 << "\" \\ \n"
		   << "\"ux" << npc+1 << "\" \\ \n"
		   << "\"uy" << npc+1 << "\" \\ \n"
		   << "\"uz" << npc+1 << "\" \\ \n"
		   << "\"Tp" << npc+1 << "\" \\ \n";
	}
      }
    }
    if (SolnBlk.Electrostatic) {
      Out_File << "\"Ex\" \\ \n"
	       << "\"Ey\" \\ \n"
	       << "\"Ez\" \\ \n"
	       << "\"V\" \\ \n";
      if (SolnBlk.Particles == PARTICLE_PHASE_EULERIAN_FORMULATION) {
	Out_File << "\"Jx\" \\ \n"
		 << "\"Jy\" \\ \n"
		 << "\"Jz\" \\ \n"
		 << "\"J\" \\ \n";
	if (NUM_CMP_PART == PARTICLE2D_MULTI_VELOCITY_FORMULATION) {
	  for (int npc = 0; npc < NUM_CMP_PART; npc++) {
	    Out_File << "\"Jx" << npc+1 << "\" \\ \n"
		     << "\"Jy" << npc+1 << "\" \\ \n"
		     << "\"Jz" << npc+1 << "\" \\ \n"
		     << "\"J"  << npc+1 << "\" \\ \n";
	  }
	}
      }
    }
  }

  for (int nr = 0; nr < numberofrotations; nr++) {
    Out_File << "ZONE T =  \"Block Number = " << Block_Number << nr
	     << "\" \\ \n"
	     << "N = " << 2*numberofnodes << " \\ \n"
	     << "E = " << (SolnBlk.ICu-SolnBlk.ICl+1)*(SolnBlk.JCu-SolnBlk.JCl+1) << " \\ \n"
	     << "F = FEPOINT \\ \n"
	     << "ET = BRICK \n";
    for (int j = SolnBlk.Grid.JNl; j <= SolnBlk.Grid.JNu; j++) {
      for (int i = SolnBlk.Grid.INl; i <= SolnBlk.Grid.INu; i++) {
	W_node = SolnBlk.Wn(i,j);
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
		 << " " << W_node.v.abs()/W_node.a() 
		 << " " << W_node.H()
		 << " " << W_node.s();
	if (SolnBlk.Particles == PARTICLE_PHASE_EULERIAN_FORMULATION) {
	  Wp = Particle2D_pState(W_node.Wp.sigma(),
				 W_node.Wp.u(),
				 W_node.Wp.Tp());
	  Out_File << " " << Wp.sigma
		   << " " << Wp.u.x
		   << " " << Wp.u.y*sin(TWO*PI*double(nr)/double(numberofrotations))
		   << " " << Wp.u.y*cos(TWO*PI*double(nr)/double(numberofrotations))
		   << " " << Wp.Tp;
	  if (NUM_CMP_PART == PARTICLE2D_MULTI_VELOCITY_FORMULATION) {
	    for (int npc = 0; npc < NUM_CMP_PART; npc++) {
	      Out_File << " " << W_node.Wp[npc].sigma 
		       << " " << W_node.Wp[npc].u.x
		       << " " << W_node.Wp[npc].u.y*sin(TWO*PI*double(nr)/double(numberofrotations))
		       << " " << W_node.Wp[npc].u.y*cos(TWO*PI*double(nr)/double(numberofrotations))
		       << " " << W_node.Wp[npc].Tp;
	    }
	  }
	}
	if (SolnBlk.Electrostatic) {
	  E_node = SolnBlk.En(i,j);
	  Out_File << " " << E_node.E.x
		   << " " << E_node.E.y*sin(TWO*PI*double(nr)/double(numberofrotations))
		   << " " << E_node.E.y*cos(TWO*PI*double(nr)/double(numberofrotations))
		   << " " << E_node.V;
	  if (SolnBlk.Particles == PARTICLE_PHASE_EULERIAN_FORMULATION) {
	    J_node = W_node.Wp.jc(W_node.qe,W_node.mp);
	    Out_File << " " << J_node.x
		     << " " << J_node.y*sin(TWO*PI*double(nr)/double(numberofrotations))
		     << " " << J_node.y*cos(TWO*PI*double(nr)/double(numberofrotations))
		     << " " << abs(J_node);
	    if (NUM_CMP_PART == PARTICLE2D_MULTI_VELOCITY_FORMULATION) {
	      for (int npc = 0; npc < NUM_CMP_PART; npc++) {
		J_node = W_node.Wp[npc].jc(W_node.qe,W_node.mp);
		Out_File << " " << J_node.x
			 << " " << J_node.y*sin(TWO*PI*double(nr)/double(numberofrotations))
			 << " " << J_node.y*cos(TWO*PI*double(nr)/double(numberofrotations))
			 << " " << abs(J_node);
	      }
	    }
	  }
	}
	Out_File << endl;
      }
    }
    if (nr < numberofrotations-1) nrr = nr + 1;
    else nrr = 0;
    for (int j = SolnBlk.Grid.JNl; j <= SolnBlk.Grid.JNu; j++) {
      for (int i = SolnBlk.Grid.INl; i <= SolnBlk.Grid.INu; i++) {
	W_node = SolnBlk.Wn(i,j);
	Out_File.setf(ios::scientific);
	Out_File << " " << SolnBlk.Grid.Node[i][j].X.x
		 << " " << SolnBlk.Grid.Node[i][j].X.y*sin(TWO*PI*double(nrr)/double(numberofrotations))
		 << " " << SolnBlk.Grid.Node[i][j].X.y*cos(TWO*PI*double(nrr)/double(numberofrotations))
		 << " " << W_node.rho
		 << " " << W_node.v.x
		 << " " << W_node.v.y*sin(TWO*PI*double(nrr)/double(numberofrotations))
		 << " " << W_node.v.y*cos(TWO*PI*double(nrr)/double(numberofrotations))
		 << " " << W_node.p
		 << " " << W_node.T()
		 << " " << W_node.v.abs()/W_node.a()
		 << " " << W_node.H()
		 << " " << W_node.s();
	if (SolnBlk.Particles == PARTICLE_PHASE_EULERIAN_FORMULATION) {
	  Wp = Particle2D_pState(W_node.Wp.sigma(),
				 W_node.Wp.u(),
				 W_node.Wp.Tp());
	  Out_File << " " << Wp.sigma
		   << " " << Wp.u.x
		   << " " << Wp.u.y*sin(TWO*PI*double(nrr)/double(numberofrotations))
		   << " " << Wp.u.y*cos(TWO*PI*double(nrr)/double(numberofrotations))
		   << " " << Wp.Tp;
	  if (NUM_CMP_PART == PARTICLE2D_MULTI_VELOCITY_FORMULATION) {
	    for (int npc = 0; npc < NUM_CMP_PART; npc++) {
	      Out_File << " " << W_node.Wp[npc].sigma 
		       << " " << W_node.Wp[npc].u.x
		       << " " << W_node.Wp[npc].u.y*sin(TWO*PI*double(nrr)/double(numberofrotations))
		       << " " << W_node.Wp[npc].u.y*cos(TWO*PI*double(nrr)/double(numberofrotations))
		       << " " << W_node.Wp[npc].Tp;
	    }
	  }
	}
	if (SolnBlk.Electrostatic) {
	  E_node = SolnBlk.En(i,j);
	  Out_File << " " << E_node.E.x
		   << " " << E_node.E.y*sin(TWO*PI*double(nrr)/double(numberofrotations))
		   << " " << E_node.E.y*cos(TWO*PI*double(nrr)/double(numberofrotations))
		   << " " << E_node.V;
	  if (SolnBlk.Particles == PARTICLE_PHASE_EULERIAN_FORMULATION) {
	    J_node = W_node.Wp.jc(W_node.qe,W_node.mp);
	    Out_File << " " << J_node.x
		     << " " << J_node.y*sin(TWO*PI*double(nrr)/double(numberofrotations))
		     << " " << J_node.y*cos(TWO*PI*double(nrr)/double(numberofrotations))
		     << " " << abs(J_node);
	    if (NUM_CMP_PART == PARTICLE2D_MULTI_VELOCITY_FORMULATION) {
	      for (int npc = 0; npc < NUM_CMP_PART; npc++) {
		J_node = W_node.Wp[npc].jc(W_node.qe,W_node.mp);
		Out_File << " " << J_node.x
			 << " " << J_node.y*sin(TWO*PI*double(nrr)/double(numberofrotations))
			 << " " << J_node.y*cos(TWO*PI*double(nrr)/double(numberofrotations))
			 << " " << abs(J_node);
	      }
	    }
	  }
	}
	Out_File << endl;
      }
    }
    
  }

  Out_File << setprecision(6);

}

/**********************************************************************
 * Routine: Output_Ringleb                                            *
 *                                                                    *
 * Writes the exact and computed Ringleb's flow solution values at    *
 * the nodes of the specified quadrilateral solution block to the     *
 * specified output stream suitable for plotting with TECPLOT.  The   *
 * error norms are also computed.                                     *
 *                                                                    *
 **********************************************************************/
void Output_Ringleb(Dusty2D_Quad_Block &SolnBlk,
		    const int Block_Number,
		    const int Output_Title,
		    ostream &Out_File,
		    double &l1_norm,
		    double &l2_norm,
		    double &max_norm,
		    double &area) {

  Dusty2D_pState Wn, We;

  // Output node solution data.  
  Out_File << setprecision(14);
  if (Output_Title)
    Out_File << "TITLE = \"" << CFFC_Name() << ": 2D Dusty Ringleb Flow Solution "
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
 * Routine: Output_Viscous_Channel                                    *
 *                                                                    *
 * Writes the exact and computed laminar channel flow solution        *
 * (Couette/Poiseuille solutions) at the nodes of the specified       *
 * quadrilateral solution block to the specified output stream        *
 * suitable for plotting with TECPLOT.  The error norms of the        *
 * u-velocity component are also computed (L1, L2, and max norms).    *
 *                                                                    *
 **********************************************************************/
void Output_Viscous_Channel(Dusty2D_Quad_Block &SolnBlk,
			    const int Block_Number,
			    const int Output_Title,
			    ostream &Out_File,
			    double &l1_norm,
			    double &l2_norm,
			    double &max_norm,
			    const Vector2D Vwall,
			    const double dp,
			    const double length,
			    const double height) {

  Dusty2D_pState We, W;
  Vector2D dX;

  // Output node solution data.  
  Out_File << setprecision(14);
  if (Output_Title)
    Out_File << "TITLE = \"" << CFFC_Name() << ": 2D Dusty Viscous Channel Flow Solution "
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
	We = ViscousChannelFlow(SolnBlk.W[i][SolnBlk.JCl],SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc,
				Vwall,dp,length,height);
	We = WallViscousHeatFlux(We,SolnBlk.Grid.nfaceS(i,SolnBlk.JCl));
      } else if (j == SolnBlk.JCl-2 && i >= SolnBlk.ICl && i <= SolnBlk.ICu) {
	// Apply the adiabtic viscous wall boundary condition to the
	// exact solution at cell (i,JCl+1).
	We = ViscousChannelFlow(SolnBlk.W[i][SolnBlk.JCl+1],SolnBlk.Grid.Cell[i][SolnBlk.JCl+1].Xc,
				Vwall,dp,length,height);
	We = WallViscousHeatFlux(We,SolnBlk.Grid.nfaceS(i,SolnBlk.JCl));
      } else if (j == SolnBlk.JCu+1 && i >= SolnBlk.ICl && i <= SolnBlk.ICu) {
	// Apply the adiabtic viscous wall boundary condition to the
	// exact solution at cell (i,JCu).
 	We = ViscousChannelFlow(SolnBlk.W[i][SolnBlk.JCu],SolnBlk.Grid.Cell[i][SolnBlk.JCu].Xc,
 				Vwall,dp,length,height);
 	We = WallViscousHeatFlux(We,SolnBlk.Grid.nfaceN(i,SolnBlk.JCu));
      } else if (j == SolnBlk.JCu+2 && i >= SolnBlk.ICl && i <= SolnBlk.ICu) {
	// Apply the adiabtic viscous wall boundary condition to the
	// exact solution at cell (i,JCu-1).
 	We = ViscousChannelFlow(SolnBlk.W[i][SolnBlk.JCu-1],SolnBlk.Grid.Cell[i][SolnBlk.JCu-1].Xc,
 				Vwall,dp,length,height);
 	We = WallViscousHeatFlux(We,SolnBlk.Grid.nfaceN(i,SolnBlk.JCu));
      } else if (i == SolnBlk.ICl-1 && j >= SolnBlk.JCl && j <= SolnBlk.JCu) {
	// Apply the subsonic inflow boundary condition to the exact
	// solution at cell (ICl,j).
 	We = ViscousChannelFlow(SolnBlk.W[SolnBlk.ICl][j],SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc,
 				Vwall,dp,length,height);
	dX = SolnBlk.Grid.Cell[SolnBlk.ICl-1][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc;
	We.p = SolnBlk.W[SolnBlk.ICl][j].p + ((SolnBlk.WoE[j].p - SolnBlk.WoW[j].p)/
					      fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc.x - 
						   SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc.x))*dX.x;
	We = ViscousChannelFlowVelocity(We,SolnBlk.Grid.Cell[i][j].Xc,
					Vwall,dp,length,height);
      } else if (i == SolnBlk.ICl-2 && j >= SolnBlk.JCl && j <= SolnBlk.JCu) {
	// Apply the subsonic inflow boundary condition to the exact
	// solution at cell (ICl,j).
 	We = ViscousChannelFlow(SolnBlk.W[SolnBlk.ICl][j],SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc,
 				Vwall,dp,length,height);
	dX = SolnBlk.Grid.Cell[SolnBlk.ICl-2][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc;
	We.p = SolnBlk.W[SolnBlk.ICl][j].p + ((SolnBlk.WoE[j].p - SolnBlk.WoW[j].p)/
					      fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc.x - 
						   SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc.x))*dX.x;
	We = ViscousChannelFlowVelocity(We,SolnBlk.Grid.Cell[i][j].Xc,
					Vwall,dp,length,height);
      } else if (i == SolnBlk.ICu+1 && j >= SolnBlk.JCl && j <= SolnBlk.JCu) {
	// Apply the subsonic outflow boundary condition to the exact
	// solution at cell (ICu,j).
 	We = ViscousChannelFlow(SolnBlk.W[SolnBlk.ICu][j],SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc,
 				Vwall,dp,length,height);
	dX = SolnBlk.Grid.Cell[SolnBlk.ICu+1][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc;
	We.p = SolnBlk.W[SolnBlk.ICu][j].p + ((SolnBlk.WoE[j].p - SolnBlk.WoW[j].p)/
					      fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc.x - 
						   SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc.x))*dX.x;
	We = ViscousChannelFlowVelocity(We,SolnBlk.Grid.Cell[i][j].Xc,
					Vwall,dp,length,height);
      } else if (i == SolnBlk.ICu+2 && j >= SolnBlk.JCl && j <= SolnBlk.JCu) {
	// Apply the subsonic outflow boundary condition to the exact
	// solution at cell (ICu,j).
 	We = ViscousChannelFlow(SolnBlk.W[SolnBlk.ICu][j],SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc,
 				Vwall,dp,length,height);
	dX = SolnBlk.Grid.Cell[SolnBlk.ICu+2][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc;
	We.p = SolnBlk.W[SolnBlk.ICu][j].p + ((SolnBlk.WoE[j].p - SolnBlk.WoW[j].p)/
					      fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc.x - 
						   SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc.x))*dX.x;
	We = ViscousChannelFlowVelocity(We,SolnBlk.Grid.Cell[i][j].Xc,
					Vwall,dp,length,height);
      } else if (i < SolnBlk.ICl && j == SolnBlk.JCl-1) {
	// Apply the subsonic inflow boundary condition to the exact
	// solution at cell (ICl,JCl) then apply the adiabatic viscous
	// wall boundary condition.
 	We = ViscousChannelFlow(SolnBlk.W[SolnBlk.ICl][SolnBlk.JCl],SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCl].Xc,
 				Vwall,dp,length,height);
	dX = SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc - SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCl].Xc;
	We.p = SolnBlk.W[SolnBlk.ICl][SolnBlk.JCl].p + ((SolnBlk.WoE[SolnBlk.JCl].p - SolnBlk.WoW[SolnBlk.JCl].p)/
							fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCl].Xc.x - 
							     SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCl].Xc.x))*dX.x;
	We = ViscousChannelFlowVelocity(We,SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc,
					Vwall,dp,length,height);
 	We = WallViscousHeatFlux(We,SolnBlk.Grid.nfaceS(i,SolnBlk.JCl));
      } else if (i < SolnBlk.ICl && j == SolnBlk.JCl-2) {
	// Apply the subsonic inflow boundary condition to the exact
	// solution at cell (ICl,JCl+1) then apply the adiabatic viscous
	// wall boundary condition.
 	We = ViscousChannelFlow(SolnBlk.W[SolnBlk.ICl][SolnBlk.JCl+1],SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCl+1].Xc,
 				Vwall,dp,length,height);
	dX = SolnBlk.Grid.Cell[i][SolnBlk.JCl+1].Xc - SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCl+1].Xc;
	We.p = SolnBlk.W[SolnBlk.ICl][SolnBlk.JCl+1].p + ((SolnBlk.WoE[SolnBlk.JCl+1].p - SolnBlk.WoW[SolnBlk.JCl+1].p)/
							fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCl+1].Xc.x - 
							     SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCl+1].Xc.x))*dX.x;
	We = ViscousChannelFlowVelocity(We,SolnBlk.Grid.Cell[i][SolnBlk.JCl+1].Xc,
					Vwall,dp,length,height);
 	We = WallViscousHeatFlux(We,SolnBlk.Grid.nfaceS(i,SolnBlk.JCl));
      } else if (i > SolnBlk.ICu && j == SolnBlk.JCl-1) {
	// Apply the subsonic inflow boundary condition to the exact
	// solution at cell (ICu,JCl) then apply the adiabatic viscous
	// wall boundary condition.
 	We = ViscousChannelFlow(SolnBlk.W[SolnBlk.ICu][SolnBlk.JCl],SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCl].Xc,
 				Vwall,dp,length,height);
	dX = SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc - SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCl].Xc;
	We.p = SolnBlk.W[SolnBlk.ICu][SolnBlk.JCl].p + ((SolnBlk.WoE[SolnBlk.JCl].p - SolnBlk.WoW[SolnBlk.JCl].p)/
							fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCl].Xc.x - 
							     SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCl].Xc.x))*dX.x;
	We = ViscousChannelFlowVelocity(We,SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc,
					Vwall,dp,length,height);
 	We = WallViscousHeatFlux(We,SolnBlk.Grid.nfaceS(i,SolnBlk.JCl));
      } else if (i > SolnBlk.ICu && j == SolnBlk.JCl-2) {
	// Apply the subsonic inflow boundary condition to the exact
	// solution at cell (ICu,JCl+1) then apply the adiabatic viscous
	// wall boundary condition.
 	We = ViscousChannelFlow(SolnBlk.W[SolnBlk.ICu][SolnBlk.JCl+1],SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCl+1].Xc,
 				Vwall,dp,length,height);
	dX = SolnBlk.Grid.Cell[i][SolnBlk.JCl+1].Xc - SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCl+1].Xc;
	We.p = SolnBlk.W[SolnBlk.ICu][SolnBlk.JCl+1].p + ((SolnBlk.WoE[SolnBlk.JCl+1].p - SolnBlk.WoW[SolnBlk.JCl+1].p)/
							fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCl+1].Xc.x - 
							     SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCl+1].Xc.x))*dX.x;
	We = ViscousChannelFlowVelocity(We,SolnBlk.Grid.Cell[i][SolnBlk.JCl+1].Xc,
					Vwall,dp,length,height);
 	We = WallViscousHeatFlux(We,SolnBlk.Grid.nfaceS(i,SolnBlk.JCl));
      } else if (i < SolnBlk.ICl && j == SolnBlk.JCu+1) {
	// Apply the subsonic inflow boundary condition to the exact
	// solution at cell (ICl,JCu) then apply the adiabatic viscous
	// wall boundary condition.
 	We = ViscousChannelFlow(SolnBlk.W[SolnBlk.ICl][SolnBlk.JCu],SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCu].Xc,
 				Vwall,dp,length,height);
	dX = SolnBlk.Grid.Cell[i][SolnBlk.JCu].Xc - SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCu].Xc;
	We.p = SolnBlk.W[SolnBlk.ICl][SolnBlk.JCu].p + ((SolnBlk.WoE[SolnBlk.JCu].p - SolnBlk.WoW[SolnBlk.JCu].p)/
							fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCu].Xc.x - 
							     SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCu].Xc.x))*dX.x;
	We = ViscousChannelFlowVelocity(We,SolnBlk.Grid.Cell[i][SolnBlk.JCu].Xc,
					Vwall,dp,length,height);
 	We = WallViscousHeatFlux(We,SolnBlk.Grid.nfaceS(i,SolnBlk.JCu));
      } else if (i < SolnBlk.ICl && j == SolnBlk.JCu+2) {
	// Apply the subsonic inflow boundary condition to the exact
	// solution at cell (ICl,JCu-1) then apply the adiabatic viscous
	// wall boundary condition.
 	We = ViscousChannelFlow(SolnBlk.W[SolnBlk.ICl][SolnBlk.JCu-1],SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCu-1].Xc,
 				Vwall,dp,length,height);
	dX = SolnBlk.Grid.Cell[i][SolnBlk.JCu-1].Xc - SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCu-1].Xc;
	We.p = SolnBlk.W[SolnBlk.ICl][SolnBlk.JCu-1].p + ((SolnBlk.WoE[SolnBlk.JCu-1].p - SolnBlk.WoW[SolnBlk.JCu-1].p)/
							fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCu-1].Xc.x - 
							     SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCu-1].Xc.x))*dX.x;
	We = ViscousChannelFlowVelocity(We,SolnBlk.Grid.Cell[i][SolnBlk.JCu-1].Xc,
					Vwall,dp,length,height);
 	We = WallViscousHeatFlux(We,SolnBlk.Grid.nfaceS(i,SolnBlk.JCu));
      } else if (i > SolnBlk.ICu && j == SolnBlk.JCu+1) {
	// Apply the subsonic inflow boundary condition to the exact
	// solution at cell (ICu,JCu) then apply the adiabatic viscous
	// wall boundary condition.
 	We = ViscousChannelFlow(SolnBlk.W[SolnBlk.ICu][SolnBlk.JCu],SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCu].Xc,
 				Vwall,dp,length,height);
	dX = SolnBlk.Grid.Cell[i][SolnBlk.JCu].Xc - SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCu].Xc;
	We.p = SolnBlk.W[SolnBlk.ICu][SolnBlk.JCu].p + ((SolnBlk.WoE[SolnBlk.JCu].p - SolnBlk.WoW[SolnBlk.JCu].p)/
							fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCu].Xc.x - 
							     SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCu].Xc.x))*dX.x;
	We = ViscousChannelFlowVelocity(We,SolnBlk.Grid.Cell[i][SolnBlk.JCu].Xc,
					Vwall,dp,length,height);
 	We = WallViscousHeatFlux(We,SolnBlk.Grid.nfaceS(i,SolnBlk.JCu));
      } else if (i > SolnBlk.ICu && j == SolnBlk.JCu+2) {
	// Apply the subsonic inflow boundary condition to the exact
	// solution at cell (ICu,JCu-1) then apply the adiabatic viscous
	// wall boundary condition.
 	We = ViscousChannelFlow(SolnBlk.W[SolnBlk.ICu][SolnBlk.JCu-1],SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCu-1].Xc,
 				Vwall,dp,length,height);
	dX = SolnBlk.Grid.Cell[i][SolnBlk.JCu-1].Xc - SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCu-1].Xc;
	We.p = SolnBlk.W[SolnBlk.ICu][SolnBlk.JCu-1].p + ((SolnBlk.WoE[SolnBlk.JCu-1].p - SolnBlk.WoW[SolnBlk.JCu-1].p)/
							fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCu-1].Xc.x - 
							     SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCu-1].Xc.x))*dX.x;
	We = ViscousChannelFlowVelocity(We,SolnBlk.Grid.Cell[i][SolnBlk.JCu-1].Xc,
					Vwall,dp,length,height);
 	We = WallViscousHeatFlux(We,SolnBlk.Grid.nfaceS(i,SolnBlk.JCu));
      } else {
	// Internal point.  Determine the exact solution and compute error norms.
 	We = ViscousChannelFlow(SolnBlk.W[i][j],SolnBlk.Grid.Cell[i][j].Xc,
 				Vwall,dp,length,height);
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
 * Routine: Output_Viscous_Pipe                                       *
 *                                                                    *
 * Writes the exact and computed laminar pipe flow solution (Hagen-   *
 * Poiseuille solution) at the nodes of the specified quadrilateral   *
 * solution block to the specified output stream suitable for         *
 * plotting with TECPLOT.  The error norms of the u-velocity          *
 * component are also computed (L1, L2, and max norms).               *
 *                                                                    *
 **********************************************************************/
void Output_Viscous_Pipe(Dusty2D_Quad_Block &SolnBlk,
			 const int Block_Number,
			 const int Output_Title,
			 ostream &Out_File,
			 double &l1_norm,
			 double &l2_norm,
			 double &max_norm,
			 const double dp,
			 const double length,
			 const double radius) {

  Dusty2D_pState We, W;
  Vector2D dX;

  // Output node solution data.  
  Out_File << setprecision(14);
  if (Output_Title)
    Out_File << "TITLE = \"" << CFFC_Name() << ": 2D Dusty Viscous Pipe Flow Solution "
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
	We = ViscousPipeFlow(SolnBlk.W[i][SolnBlk.JCl],SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc,
			     dp,length,radius);
	We = Reflect(We,SolnBlk.Grid.nfaceS(i,SolnBlk.JCl));
      } else if (j == SolnBlk.JCl-2 && i >= SolnBlk.ICl && i <= SolnBlk.ICu) {
	// Apply the reflection boundary condition to the exact solution
	// at cell (i,JCl+1).
	We = ViscousPipeFlow(SolnBlk.W[i][SolnBlk.JCl+1],SolnBlk.Grid.Cell[i][SolnBlk.JCl+1].Xc,
			     dp,length,radius);
	We = Reflect(We,SolnBlk.Grid.nfaceS(i,SolnBlk.JCl));
      } else if (j == SolnBlk.JCu+1 && i >= SolnBlk.ICl && i <= SolnBlk.ICu) {
	// Apply the adiabtic viscous wall boundary condition to the
	// exact solution at cell (i,JCu).
 	We = ViscousPipeFlow(SolnBlk.W[i][SolnBlk.JCu],SolnBlk.Grid.Cell[i][SolnBlk.JCu].Xc,
			     dp,length,radius);
 	We = WallViscousHeatFlux(We,SolnBlk.Grid.nfaceN(i,SolnBlk.JCu));
      } else if (j == SolnBlk.JCu+2 && i >= SolnBlk.ICl && i <= SolnBlk.ICu) {
	// Apply the adiabtic viscous wall boundary condition to the
	// exact solution at cell (i,JCu-1).
 	We = ViscousPipeFlow(SolnBlk.W[i][SolnBlk.JCu-1],SolnBlk.Grid.Cell[i][SolnBlk.JCu-1].Xc,
			     dp,length,radius);
 	We = WallViscousHeatFlux(We,SolnBlk.Grid.nfaceN(i,SolnBlk.JCu));
      } else if (i == SolnBlk.ICl-1 && j >= SolnBlk.JCl && j <= SolnBlk.JCu) {
	// Apply the subsonic inflow boundary condition to the exact
	// solution at cell (ICl,j).
 	We = ViscousPipeFlow(SolnBlk.W[SolnBlk.ICl][j],SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc,
			     dp,length,radius);
	dX = SolnBlk.Grid.Cell[SolnBlk.ICl-1][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc;
	We.p = SolnBlk.W[SolnBlk.ICl][j].p + ((SolnBlk.WoE[j].p - SolnBlk.WoW[j].p)/
					      fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc.x - 
						   SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc.x))*dX.x;
      } else if (i == SolnBlk.ICl-2 && j >= SolnBlk.JCl && j <= SolnBlk.JCu) {
	// Apply the subsonic inflow boundary condition to the exact
	// solution at cell (ICl,j).
 	We = ViscousPipeFlow(SolnBlk.W[SolnBlk.ICl][j],SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc,
 				dp,length,radius);
	dX = SolnBlk.Grid.Cell[SolnBlk.ICl-2][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc;
	We.p = SolnBlk.W[SolnBlk.ICl][j].p + ((SolnBlk.WoE[j].p - SolnBlk.WoW[j].p)/
					      fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc.x - 
						   SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc.x))*dX.x;
      } else if (i == SolnBlk.ICu+1 && j >= SolnBlk.JCl && j <= SolnBlk.JCu) {
	// Apply the subsonic outflow boundary condition to the exact
	// solution at cell (ICu,j).
 	We = ViscousPipeFlow(SolnBlk.W[SolnBlk.ICu][j],SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc,
 				dp,length,radius);
	dX = SolnBlk.Grid.Cell[SolnBlk.ICu+1][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc;
	We.p = SolnBlk.W[SolnBlk.ICu][j].p + ((SolnBlk.WoE[j].p - SolnBlk.WoW[j].p)/
					      fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc.x - 
						   SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc.x))*dX.x;
      } else if (i == SolnBlk.ICu+2 && j >= SolnBlk.JCl && j <= SolnBlk.JCu) {
	// Apply the subsonic outflow boundary condition to the exact
	// solution at cell (ICu,j).
 	We = ViscousPipeFlow(SolnBlk.W[SolnBlk.ICu][j],SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc,
 				dp,length,radius);
	dX = SolnBlk.Grid.Cell[SolnBlk.ICu+2][j].Xc - SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc;
	We.p = SolnBlk.W[SolnBlk.ICu][j].p + ((SolnBlk.WoE[j].p - SolnBlk.WoW[j].p)/
					      fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc.x - 
						   SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc.x))*dX.x;
      } else if (i < SolnBlk.ICl && j == SolnBlk.JCl-1) {
	// Apply the subsonic inflow boundary condition to the exact
	// solution at cell (ICl,JCl) then apply the reflection
	// boundary condition.
 	We = ViscousPipeFlow(SolnBlk.W[SolnBlk.ICl][SolnBlk.JCl],SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCl].Xc,
 				dp,length,radius);
	dX = SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc - SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCl].Xc;
	We.p = SolnBlk.W[SolnBlk.ICl][SolnBlk.JCl].p + ((SolnBlk.WoE[SolnBlk.JCl].p - SolnBlk.WoW[SolnBlk.JCl].p)/
							fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCl].Xc.x - 
							     SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCl].Xc.x))*dX.x;
 	We = Reflect(We,SolnBlk.Grid.nfaceS(i,SolnBlk.JCl));
      } else if (i < SolnBlk.ICl && j == SolnBlk.JCl-2) {
	// Apply the subsonic inflow boundary condition to the exact
	// solution at cell (ICl,JCl+1) then apply the reflection
	// boundary condition.
 	We = ViscousPipeFlow(SolnBlk.W[SolnBlk.ICl][SolnBlk.JCl+1],SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCl+1].Xc,
 				dp,length,radius);
	dX = SolnBlk.Grid.Cell[i][SolnBlk.JCl+1].Xc - SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCl+1].Xc;
	We.p = SolnBlk.W[SolnBlk.ICl][SolnBlk.JCl+1].p + ((SolnBlk.WoE[SolnBlk.JCl+1].p - SolnBlk.WoW[SolnBlk.JCl+1].p)/
							fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCl+1].Xc.x - 
							     SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCl+1].Xc.x))*dX.x;
 	We = Reflect(We,SolnBlk.Grid.nfaceS(i,SolnBlk.JCl));
      } else if (i > SolnBlk.ICu && j == SolnBlk.JCl-1) {
	// Apply the subsonic inflow boundary condition to the exact
	// solution at cell (ICu,JCl) then apply the adiabatic viscous
	// wall boundary condition.
 	We = ViscousPipeFlow(SolnBlk.W[SolnBlk.ICu][SolnBlk.JCl],SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCl].Xc,
 				dp,length,radius);
	dX = SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc - SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCl].Xc;
	We.p = SolnBlk.W[SolnBlk.ICu][SolnBlk.JCl].p + ((SolnBlk.WoE[SolnBlk.JCl].p - SolnBlk.WoW[SolnBlk.JCl].p)/
							fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCl].Xc.x - 
							     SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCl].Xc.x))*dX.x;
 	We = Reflect(We,SolnBlk.Grid.nfaceS(i,SolnBlk.JCl));
      } else if (i > SolnBlk.ICu && j == SolnBlk.JCl-2) {
	// Apply the subsonic inflow boundary condition to the exact
	// solution at cell (ICu,JCl+1) then apply the adiabatic viscous
	// wall boundary condition.
 	We = ViscousPipeFlow(SolnBlk.W[SolnBlk.ICu][SolnBlk.JCl+1],SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCl+1].Xc,
 				dp,length,radius);
	dX = SolnBlk.Grid.Cell[i][SolnBlk.JCl+1].Xc - SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCl+1].Xc;
	We.p = SolnBlk.W[SolnBlk.ICu][SolnBlk.JCl+1].p + ((SolnBlk.WoE[SolnBlk.JCl+1].p - SolnBlk.WoW[SolnBlk.JCl+1].p)/
							fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCl+1].Xc.x - 
							     SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCl+1].Xc.x))*dX.x;
 	We = Reflect(We,SolnBlk.Grid.nfaceS(i,SolnBlk.JCl));
      } else if (i < SolnBlk.ICl && j == SolnBlk.JCu+1) {
	// Apply the subsonic inflow boundary condition to the exact
	// solution at cell (ICl,JCu) then apply the adiabatic viscous
	// wall boundary condition.
 	We = ViscousPipeFlow(SolnBlk.W[SolnBlk.ICl][SolnBlk.JCu],SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCu].Xc,
 				dp,length,radius);
	dX = SolnBlk.Grid.Cell[i][SolnBlk.JCu].Xc - SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCu].Xc;
	We.p = SolnBlk.W[SolnBlk.ICl][SolnBlk.JCu].p + ((SolnBlk.WoE[SolnBlk.JCu].p - SolnBlk.WoW[SolnBlk.JCu].p)/
							fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCu].Xc.x - 
							     SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCu].Xc.x))*dX.x;
 	We = WallViscousHeatFlux(We,SolnBlk.Grid.nfaceS(i,SolnBlk.JCu));
      } else if (i < SolnBlk.ICl && j == SolnBlk.JCu+2) {
	// Apply the subsonic inflow boundary condition to the exact
	// solution at cell (ICl,JCu-1) then apply the adiabatic viscous
	// wall boundary condition.
 	We = ViscousPipeFlow(SolnBlk.W[SolnBlk.ICl][SolnBlk.JCu-1],SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCu-1].Xc,
 				dp,length,radius);
	dX = SolnBlk.Grid.Cell[i][SolnBlk.JCu-1].Xc - SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCu-1].Xc;
	We.p = SolnBlk.W[SolnBlk.ICl][SolnBlk.JCu-1].p + ((SolnBlk.WoE[SolnBlk.JCu-1].p - SolnBlk.WoW[SolnBlk.JCu-1].p)/
							fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCu-1].Xc.x - 
							     SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCu-1].Xc.x))*dX.x;
 	We = WallViscousHeatFlux(We,SolnBlk.Grid.nfaceS(i,SolnBlk.JCu));
      } else if (i > SolnBlk.ICu && j == SolnBlk.JCu+1) {
	// Apply the subsonic inflow boundary condition to the exact
	// solution at cell (ICu,JCu) then apply the adiabatic viscous
	// wall boundary condition.
 	We = ViscousPipeFlow(SolnBlk.W[SolnBlk.ICu][SolnBlk.JCu],SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCu].Xc,
 				dp,length,radius);
	dX = SolnBlk.Grid.Cell[i][SolnBlk.JCu].Xc - SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCu].Xc;
	We.p = SolnBlk.W[SolnBlk.ICu][SolnBlk.JCu].p + ((SolnBlk.WoE[SolnBlk.JCu].p - SolnBlk.WoW[SolnBlk.JCu].p)/
							fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCu].Xc.x - 
							     SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCu].Xc.x))*dX.x;
 	We = WallViscousHeatFlux(We,SolnBlk.Grid.nfaceS(i,SolnBlk.JCu));
      } else if (i > SolnBlk.ICu && j == SolnBlk.JCu+2) {
	// Apply the subsonic inflow boundary condition to the exact
	// solution at cell (ICu,JCu-1) then apply the adiabatic viscous
	// wall boundary condition.
 	We = ViscousPipeFlow(SolnBlk.W[SolnBlk.ICu][SolnBlk.JCu-1],SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCu-1].Xc,
 				dp,length,radius);
	dX = SolnBlk.Grid.Cell[i][SolnBlk.JCu-1].Xc - SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCu-1].Xc;
	We.p = SolnBlk.W[SolnBlk.ICu][SolnBlk.JCu-1].p + ((SolnBlk.WoE[SolnBlk.JCu-1].p - SolnBlk.WoW[SolnBlk.JCu-1].p)/
							fabs(SolnBlk.Grid.Cell[SolnBlk.ICu][SolnBlk.JCu-1].Xc.x - 
							     SolnBlk.Grid.Cell[SolnBlk.ICl][SolnBlk.JCu-1].Xc.x))*dX.x;
 	We = WallViscousHeatFlux(We,SolnBlk.Grid.nfaceS(i,SolnBlk.JCu));
      } else {
	// Internal point.  Determine the exact solution and compute error norms.
 	We = ViscousPipeFlow(SolnBlk.W[i][j],SolnBlk.Grid.Cell[i][j].Xc,
 			     dp,length,radius);
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
void Output_Turbulent_Pipe_Tecplot(Dusty2D_Quad_Block &SolnBlk,
				   const int Block_Number,
				   const int Output_Title,
				   const int Output_Data,
				   ostream &Out_File,
				   const double &Re,
				   const double &Pipe_Radius,
				   const int &variable_flag) {

  Dusty2D_pState We, W;

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
 * Routine: Output_Flat_Plate                                         *
 *                                                                    *
 * This routine outputs the non-dimensionalized computed flat plate   *
 * solution and the corresponding Blasius solution.                   *
 *                                                                    *
 **********************************************************************/
void Output_Flat_Plate(Dusty2D_Quad_Block &SolnBlk,
		       const int Block_Number,
		       const int Output_Title_Soln,
		       ostream &Out_File_Soln,
		       const int Output_Title_Skin,
		       ostream &Out_File_Skin,
		       const Dusty2D_pState &Winf,
		       const double &plate_length,
		       double &l1_norm,
		       double &l2_norm,
		       double &max_norm,
		       double &area,
		       int &numberofcells,
		       double &l1_norm_cf,
		       double &l2_norm_cf,
		       double &max_norm_cf,
		       double &area_cf,
		       int &numberofcells_cf) {

  int skin_friction_flag = OFF;
  Dusty2D_pState W, We;
  Vector2D X;
  double eta, f, fp, fpp, Rex, linf, Cf, Cf2, Cfe, phiv;

  // Output node solution data.  
  Out_File_Soln << setprecision(14);
  if (Output_Title_Soln) {
    Out_File_Soln << "TITLE = \"" << CFFC_Name() << ": 2D Dusty Solution, "
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
		  << "\"(vx - vx_e)/Vxinf\" \\ \n"
		  << "\"phiv\" \\ \n"
		  << "\"phiv_e\" \\ \n"
		  << "\"phiv - phiv_e\" \\ \n";
  }
  Out_File_Soln << "ZONE T =  \"Block Number = " << Block_Number << "\" \\ \n"
		<< "I = " << SolnBlk.Grid.ICu - SolnBlk.Grid.ICl + 1 << " \\ \n"
		<< "J = " << SolnBlk.Grid.JCu - SolnBlk.Grid.JCl + 1 << " \\ \n"
		<< "F = POINT \n";
  for (int j = SolnBlk.Grid.JCl; j <= SolnBlk.Grid.JCu; j++) {
    for (int i = SolnBlk.Grid.ICl; i <= SolnBlk.Grid.ICu; i++) {
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
  for (int j = SolnBlk.Grid.JCl; j <= SolnBlk.Grid.JCu; j++) {
    for (int i = SolnBlk.Grid.ICl; i <= SolnBlk.Grid.ICu; i++) {
      // Get cell position and solution data.
      X = SolnBlk.Grid.Cell[i][j].Xc;
      W = SolnBlk.W[i][j];
      // Determine the norms of the u-velocity component.
      if (X.x >= ZERO && X.x <= plate_length) {
	We = FlatPlate(Winf,X,plate_length,eta,f,fp,fpp);
	l1_norm += fabs(W.v.x - We.v.x);
	l2_norm += sqr(W.v.x - We.v.x);
	max_norm = max(max_norm,fabs(W.v.x - We.v.x));
	area += SolnBlk.Grid.Cell[i][j].A;
	numberofcells++;
      }
      // Determine the norms of the skin friction coefficient.
      if (X.x >= ZERO && j == 2 &&
	  X.x <= plate_length && SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX) {
	// Get exact skin friction coefficient.
	Rex = (Winf.v.x/Winf.nu())*((X.x+NANO)/ONE);
	Cfe = TWO*0.32206/sqrt(Rex);
	// Get computed skin friction coefficient.
	Cf  = TWO*WallShearStress(W,X,
				  SolnBlk.Grid.nodeSW(i,j).X,
				  SolnBlk.Grid.nodeSE(i,j).X,
				  -SolnBlk.Grid.nfaceS(i,j))/(Winf.rho*Winf.v.x*Winf.v.x);
	// Calculate error norms.
	l1_norm_cf += fabs(Cf - Cfe);
	l2_norm_cf += sqr(Cf - Cfe);
	max_norm_cf = max(max_norm_cf,fabs(Cf - Cfe));
	area_cf += SolnBlk.Grid.Cell[i][j].A;
	numberofcells_cf++;
	if (!skin_friction_flag) skin_friction_flag = ON;
      } else {
	Cf  = ZERO;
	Cfe = ZERO;
      }
    }
  }

  if (Output_Title_Skin) {
    Out_File_Skin << "TITLE = \"" << CFFC_Name() << ": 2D Dusty Solution, "
		  << "\"" << "\n"
		  << "VARIABLES = \"x\" \\ \n"
		  << "\"Rex\" \\ \n"
		  << "\"Cf\" \\ \n"
		  << "\"Cf_e\" \\ \n"
		  << "\"Cf2\" \\ \n"
		  << "\"Cf-Cf_e\" \\ \n";
  }
  if (skin_friction_flag) {
    Out_File_Skin << "ZONE T =  \"Block Number = " << Block_Number
		  << "\" \\ \n"
		  << "I = " << SolnBlk.Grid.ICu-SolnBlk.Grid.ICl+2 << " \\ \n"
		  << "J = " << 1 << " \\ \n"
		  << "F = POINT \n";
    for (int i = SolnBlk.Grid.ICl-1; i <= SolnBlk.Grid.ICu+1; i++) {
      if (SolnBlk.Grid.Cell[i][SolnBlk.Grid.JCl].Xc.x >= ZERO &&
	  SolnBlk.Grid.Cell[i][SolnBlk.Grid.JCl].Xc.x <= plate_length &&
	  SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX) {
	Rex = (Winf.v.x/Winf.nu())*SolnBlk.Grid.Cell[i][SolnBlk.Grid.JCl].Xc.x;
	Cfe = TWO*0.33206/max(sqrt(Rex),NANO);
	Cf  = TWO*WallShearStress(SolnBlk.W[i][SolnBlk.Grid.JCl],
				  SolnBlk.Grid.Cell[i][SolnBlk.Grid.JCl].Xc,
				  SolnBlk.Grid.nodeSW(i,SolnBlk.Grid.JCl).X,
				  SolnBlk.Grid.nodeSE(i,SolnBlk.Grid.JCl).X,
				  -SolnBlk.Grid.nfaceS(i,SolnBlk.Grid.JCl))/(Winf.rho*sqr(Winf.v.x));
	Cf2 = TWO*ShearStress(SolnBlk.W[i][SolnBlk.Grid.JCl],
			      SolnBlk.dWdx[i][SolnBlk.Grid.JCl],
			      SolnBlk.dWdy[i][SolnBlk.Grid.JCl],
			      -SolnBlk.Grid.nfaceS(i,SolnBlk.Grid.JCl))/(Winf.rho*sqr(Winf.v.x));
	Out_File_Skin.setf(ios::scientific);
	Out_File_Skin << " " << SolnBlk.Grid.Cell[i][SolnBlk.Grid.JCl].Xc.x
		      << " " << Rex
		      << " " << Cf
		      << " " << Cfe
		      << " " << Cf - Cfe
		      << endl;
      }
    }
  }

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
void Output_Driven_Cavity_Flow(Dusty2D_Quad_Block &SolnBlk,
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

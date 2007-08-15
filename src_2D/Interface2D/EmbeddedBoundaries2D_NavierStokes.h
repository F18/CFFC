/**********************************************************************
 * EmbeddedBoundaries2D: Header file declaring 2D embedded boundary   *
 *                       classes and functions.                       *
 **********************************************************************/

#ifndef _EMBEDDEDBOUNDARIES2D_NAVIERSTOKES2D_INCLUDED
#define _EMBEDDEDBOUNDARIES2D_NAVIERSTOKES2D_INCLUDED

// Include 2D Navier-Stokes quadrilateral mesh solution header file.

#ifndef _NAVIERSTOKES2D_QUAD_INCLUDED
#include "../NavierStokes2D/NavierStokes2DQuad.h"
#endif // _NAVIERSTOKES2D_QUAD_INCLUDED

// Include Embedded Boundaries input header file.

#ifndef _EMBEDDEDBOUNDARIES2D_INCLUDED
#include "EmbeddedBoundaries2D.h"
#endif // _EMBEDDEDBOUNDARIES2D_INCLUDED

/**********************************************************************
 * EmbeddedBoundaries2D::Calculate_Refinement_Criteria --             *
 **********************************************************************/
template <> inline void EmbeddedBoundaries2D<NavierStokes2D_cState,
					     NavierStokes2D_pState,
					     NavierStokes2D_Quad_Block,
					     NavierStokes2D_Input_Parameters>::
Calculate_Refinement_Criteria(const int &nb, const int &i, const int &j,
			      double *refinement_criteria) {

  double grad_rho_x, grad_rho_y, grad_rho_abs, grad_rho_criteria,
         div_V, div_V_criteria, curl_V_z, curl_V_abs, curl_V_criteria;
  int refinement_criteria_number;

  // Refinement criteria that can be used:
  // (1) Refinement criteria based on the gradient of the density field;
  // (2) Refinement criteria based on the divergence of the velocity vector;
  // (3) Refinement criteria based on the curl of the velocity vector;

  // Reconstruct the solution within the cell.
  Linear_Least_Squares_Reconstruction(nb,i,j);

  // Evaluate refinement criteria #1 based on the gradient of the 
  // density field.
  if (IP->Refinement_Criteria_Gradient_Density) {
    grad_rho_x = Local_SolnBlk[nb].dWdx[i][j].rho;
    grad_rho_y = Local_SolnBlk[nb].dWdy[i][j].rho;
    grad_rho_abs = sqrt(sqr(grad_rho_x) + sqr(grad_rho_y));
    grad_rho_criteria = sqrt(Local_SolnBlk[nb].Grid.Cell[i][j].A)*grad_rho_abs/Local_SolnBlk[nb].W[i][j].rho;
  } else {
    grad_rho_criteria = ZERO;
  }

  // Evaluate refinement criteria #2 based on the divergence of the
  // velocity vector.
  if (IP->Refinement_Criteria_Divergence_Velocity) {
    div_V = Local_SolnBlk[nb].dWdx[i][j].v.x + Local_SolnBlk[nb].dWdy[i][j].v.y;
    div_V_criteria = sqrt(Local_SolnBlk[nb].Grid.Cell[i][j].A)*fabs(div_V)/Local_SolnBlk[nb].W[i][j].a();
  } else {
    div_V_criteria = ZERO;
  }

  // Evaluate refinement criteria #3 based on the curl of the
  // velocity vector.
  if (IP->Refinement_Criteria_Curl_Velocity) {
    curl_V_z = Local_SolnBlk[nb].dWdx[i][j].v.y - Local_SolnBlk[nb].dWdy[i][j].v.x; 
    //curl_V_z = Local_SolnBlk[nb].dWdy[i][j].v.x;
    curl_V_abs = sqrt(sqr(curl_V_z)); 
    curl_V_criteria = sqrt(Local_SolnBlk[nb].Grid.Cell[i][j].A)*curl_V_abs/Local_SolnBlk[nb].W[i][j].a();
  } else {
    curl_V_criteria = ZERO;
  }

  // Set the refinement criteria.
  refinement_criteria_number = 0;
  if (IP->Refinement_Criteria_Gradient_Density) {
    refinement_criteria[refinement_criteria_number] = max(refinement_criteria[refinement_criteria_number],
							  grad_rho_criteria);
    refinement_criteria_number++;
  }
  if (IP->Refinement_Criteria_Divergence_Velocity) {
    refinement_criteria[refinement_criteria_number] = max(refinement_criteria[refinement_criteria_number],
							  div_V_criteria);
    refinement_criteria_number++;
  }
  if (IP->Refinement_Criteria_Curl_Velocity) {
    refinement_criteria[refinement_criteria_number] = max(refinement_criteria[refinement_criteria_number],
							  curl_V_criteria);
    refinement_criteria_number++;
  }

}

/**********************************************************************
 * EmbeddedBoundaries2D::Output_Flat_Plate_Tecplot --                 *
 **********************************************************************/
template <> inline void EmbeddedBoundaries2D<NavierStokes2D_cState,
					     NavierStokes2D_pState,
					     NavierStokes2D_Quad_Block,
					     NavierStokes2D_Input_Parameters>::
Output_Flat_Plate_Tecplot(const int &nb,
			  const int &Output_Title_Soln,
			  ostream &Out_File_Soln,
			  const int &Output_Title_Skin,
			  ostream &Out_File_Skin,
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
  NavierStokes2D_pState W, We;
  Vector2D X;
  double eta, f, fp, fpp, Rex, linf, Cf, Cf2, Cfe, phiv;

  Linear_Least_Squares_Reconstruction(nb);

  // Output node solution data.  
  Out_File_Soln << setprecision(14);
  if (Output_Title_Soln) {
    Out_File_Soln << "TITLE = \"" << CFFC_Name() << ": 2D Flat Plate Solution, "
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
		  << "\"phiv - phiv_e\" \\ \n"
		  << "\"dudx\" \\ \n"
		  << "\"dvdx\" \\ \n"
		  << "\"dudy\" \\ \n"
		  << "\"dvdy\" \\ \n";
  }
  Out_File_Soln << "ZONE T =  \"Block Number = " << Local_Solution_Block_List->Block[nb].gblknum << "\" \\ \n"
		<< "I = " << Local_SolnBlk[nb].ICu - Local_SolnBlk[nb].ICl + 1 << " \\ \n"
		<< "J = " << Local_SolnBlk[nb].JCu - Local_SolnBlk[nb].JCl + 1 << " \\ \n"
		<< "F = POINT \n";
  for (int j = Local_SolnBlk[nb].JCl; j <= Local_SolnBlk[nb].JCu; j++) {
    for (int i = Local_SolnBlk[nb].ICl; i <= Local_SolnBlk[nb].ICu; i++) {
      X = Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
      W = Local_SolnBlk[nb].W[i][j];
      if (Mesh[nb].cell_status[i][j] == CELL_STATUS_ACTIVE) {
	We = FlatPlate(IP->Wo,X,IP->Plate_Length,eta,f,fp,fpp);
	if (X.x > ZERO && X.x < IP->Plate_Length) phiv = (W.v.y/IP->Wo.v.x)*sqrt(IP->Wo.v.x*X.x/IP->Wo.nu());
	else phiv = ZERO;
      } else {
	We = W; eta = ZERO; f = ZERO; fp = ZERO; fpp = ZERO; phiv = ZERO;
      }
      // Output data.
      Out_File_Soln.setf(ios::scientific);
      Out_File_Soln << " " << X
		    << " " << W[1] 
		    << W.v
		    << " " << W.p
		    << " " << We[1]
		    << We.v
		    << " " << We.p
		    << " " << eta
		    << " " << f
		    << " " << fp
		    << " " << fpp
		    << " " << W.v.x/IP->Wo.v.x
		    << " " << We.v.x/IP->Wo.v.x
		    << " " << (W.v.x - We.v.x)/IP->Wo.v.x
		    << " " << phiv
		    << " " << HALF*(eta*fp-f)
		    << " " << phiv - HALF*(eta*fp-f)
		    << " " << Local_SolnBlk[nb].dWdx[i][j].v.x
		    << " " << Local_SolnBlk[nb].dWdx[i][j].v.y
		    << " " << Local_SolnBlk[nb].dWdy[i][j].v.x
		    << " " << Local_SolnBlk[nb].dWdy[i][j].v.y
		    << endl;
      Out_File_Soln.unsetf(ios::scientific);
    }
  }

  // Calculate the norms of the u-velocity component and the skin
  // friction coefficient.
  for (int j = Local_SolnBlk[nb].JCl; j <= Local_SolnBlk[nb].JCu; j++) {
    for (int i = Local_SolnBlk[nb].ICl; i <= Local_SolnBlk[nb].ICu; i++) {
      if (Mesh[nb].cell_status[i][j] == CELL_STATUS_ACTIVE) {
	// Get cell position and solution data.
	X = Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	W = Local_SolnBlk[nb].W[i][j];
	// Determine the norms of the u-velocity component.
	if (X.x >= ZERO && X.x <= IP->Plate_Length) {
	  We = FlatPlate(IP->Wo,X,IP->Plate_Length,eta,f,fp,fpp);
	  l1_norm += fabs(W.v.x - We.v.x)*Local_SolnBlk[nb].Grid.Cell[i][j].A;
	  l2_norm += sqr(W.v.x - We.v.x)*Local_SolnBlk[nb].Grid.Cell[i][j].A;
	  max_norm = max(max_norm,fabs(W.v.x - We.v.x));
	  area += Local_SolnBlk[nb].Grid.Cell[i][j].A;
	  numberofcells++;
	}
	// Determine the norms of the skin friction coefficient.
	if (!Interface_Component_List.Ni && X.x >= ZERO && j == 2 &&
	    X.x <= IP->Plate_Length && Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX) {
	  // Get exact skin friction coefficient.
	  Rex = (IP->Wo.v.x/IP->Wo.nu())*((X.x+NANO)/ONE);
	  Cfe = TWO*0.32206/sqrt(Rex);
	  // Get computed skin friction coefficient.
	  Cf  = TWO*WallShearStress(W,X,
				    Local_SolnBlk[nb].Grid.nodeSW(i,j).X,
				    Local_SolnBlk[nb].Grid.nodeSE(i,j).X,
				    -Local_SolnBlk[nb].Grid.nfaceS(i,j))/(IP->Wo.rho*IP->Wo.v.x*IP->Wo.v.x);
	  // Calculate error norms.
	  l1_norm_cf += fabs(Cf - Cfe)*Local_SolnBlk[nb].Grid.Cell[i][j].A;
	  l2_norm_cf += sqr(Cf - Cfe)*Local_SolnBlk[nb].Grid.Cell[i][j].A;
	  max_norm_cf = max(max_norm_cf,fabs(Cf - Cfe));
	  area_cf += Local_SolnBlk[nb].Grid.Cell[i][j].A;
	  numberofcells_cf++;
	  if (!skin_friction_flag) skin_friction_flag = ON;
 	} else if (Interface_Component_List.Ni &&
		   (Mesh[nb].cell_status[i-1][j  ] != CELL_STATUS_ACTIVE ||
		    Mesh[nb].cell_status[i+1][j  ] != CELL_STATUS_ACTIVE ||
		    Mesh[nb].cell_status[i  ][j-1] != CELL_STATUS_ACTIVE ||
		    Mesh[nb].cell_status[i  ][j+1] != CELL_STATUS_ACTIVE) &&
		   X.x >= ZERO && X.x <= IP->Plate_Length) {
	  if (!skin_friction_flag) skin_friction_flag = ON;
 	  if (Mesh[nb].cell_status[i-1][j] != CELL_STATUS_ACTIVE &&
 	      Local_SolnBlk[nb].Grid.lfaceW(i,j) > TOLER) {
 	    Cf  = TWO*WallShearStress(W,X,
 				      Local_SolnBlk[nb].Grid.nodeNW(i,j).X,
 				      Local_SolnBlk[nb].Grid.nodeSW(i,j).X,
 				      -Local_SolnBlk[nb].Grid.nfaceW(i,j))/(IP->Wo.rho*sqr(IP->Wo.v.x));
	    l1_norm_cf += fabs(Cf - Cfe)*Local_SolnBlk[nb].Grid.Cell[i][j].A;
	    l2_norm_cf += sqr(Cf - Cfe)*Local_SolnBlk[nb].Grid.Cell[i][j].A;
	    max_norm_cf = max(max_norm_cf,fabs(Cf - Cfe));
	    area_cf += Local_SolnBlk[nb].Grid.Cell[i][j].A;
	    numberofcells_cf++;
 	  } else if (Mesh[nb].cell_status[i+1][j] != CELL_STATUS_ACTIVE &&
 		     Local_SolnBlk[nb].Grid.lfaceE(i,j) > TOLER) {
 	    Cf  = TWO*WallShearStress(W,X,
				      Local_SolnBlk[nb].Grid.nodeSE(i,j).X,
				      Local_SolnBlk[nb].Grid.nodeNE(i,j).X,
				      -Local_SolnBlk[nb].Grid.nfaceE(i,j))/(IP->Wo.rho*sqr(IP->Wo.v.x));
	    l1_norm_cf += fabs(Cf - Cfe)*Local_SolnBlk[nb].Grid.Cell[i][j].A;
	    l2_norm_cf += sqr(Cf - Cfe)*Local_SolnBlk[nb].Grid.Cell[i][j].A;
	    max_norm_cf = max(max_norm_cf,fabs(Cf - Cfe));
	    area_cf += Local_SolnBlk[nb].Grid.Cell[i][j].A;
	    numberofcells_cf++;
	  } else if (Mesh[nb].cell_status[i][j-1] != CELL_STATUS_ACTIVE &&
		     Local_SolnBlk[nb].Grid.lfaceS(i,j) > TOLER) {
	    Cf  = TWO*WallShearStress(W,X,
				      Local_SolnBlk[nb].Grid.nodeSW(i,j).X,
				      Local_SolnBlk[nb].Grid.nodeSE(i,j).X,
				      -Local_SolnBlk[nb].Grid.nfaceS(i,j))/(IP->Wo.rho*sqr(IP->Wo.v.x));
	    l1_norm_cf += fabs(Cf - Cfe)*Local_SolnBlk[nb].Grid.Cell[i][j].A;
	    l2_norm_cf += sqr(Cf - Cfe)*Local_SolnBlk[nb].Grid.Cell[i][j].A;
	    max_norm_cf = max(max_norm_cf,fabs(Cf - Cfe));
	    area_cf += Local_SolnBlk[nb].Grid.Cell[i][j].A;
	    numberofcells_cf++;
	  } else if (Mesh[nb].cell_status[i][j+1] != CELL_STATUS_ACTIVE &&
		     Local_SolnBlk[nb].Grid.lfaceN(i,j) > TOLER) {
	    Cf  = TWO*WallShearStress(W,X,
				      Local_SolnBlk[nb].Grid.nodeNE(i,j).X,
				      Local_SolnBlk[nb].Grid.nodeNW(i,j).X,
				      -Local_SolnBlk[nb].Grid.nfaceN(i,j))/(IP->Wo.rho*sqr(IP->Wo.v.x));
	    l1_norm_cf += fabs(Cf - Cfe)*Local_SolnBlk[nb].Grid.Cell[i][j].A;
	    l2_norm_cf += sqr(Cf - Cfe)*Local_SolnBlk[nb].Grid.Cell[i][j].A;
	    max_norm_cf = max(max_norm_cf,fabs(Cf - Cfe));
	    area_cf += Local_SolnBlk[nb].Grid.Cell[i][j].A;
	    numberofcells_cf++;
	  }
	} else {
	  Cf  = ZERO;
	  Cfe = ZERO;
	}
      }
    }
  }

  if (Output_Title_Skin) {
    Out_File_Skin << "TITLE = \"" << CFFC_Name() << ": Flat Plate Skin Friction Coefficient, "
		  << "\"" << "\n"
		  << "VARIABLES = \"x\" \\ \n"
		  << "\"Rex\" \\ \n"
		  << "\"Cf\" \\ \n"
		  << "\"Cf2\" \\ \n"
		  << "\"Cf_e\" \\ \n"
		  << "\"Cf-Cf_e\" \\ \n"
		  << "\"Cf2-Cf_e\" \\ \n"
		  << "\"A\" \\ \n";
  }
  if (skin_friction_flag) {
    if (!Interface_Component_List.Ni) {
      Out_File_Skin << "ZONE T =  \"Block Number = " << nb
		    << "\" \\ \n"
		    << "I = " << Local_SolnBlk[nb].ICu-Local_SolnBlk[nb].ICl+2 << " \\ \n"
		    << "J = " << 1 << " \\ \n"
		    << "F = POINT \n";
      for (int i = Local_SolnBlk[nb].ICl-1; i <= Local_SolnBlk[nb].ICu+1; i++) {
	if (Local_SolnBlk[nb].Grid.Cell[i][Local_SolnBlk[nb].JCl].Xc.x >= ZERO &&
	    Local_SolnBlk[nb].Grid.Cell[i][Local_SolnBlk[nb].JCl].Xc.x <= IP->Plate_Length &&
	    Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX) {
	  Rex = (IP->Wo.v.x/IP->Wo.nu())*Local_SolnBlk[nb].Grid.xfaceS(i,Local_SolnBlk[nb].JCl).x;
	  Cfe = TWO*0.33206/max(sqrt(Rex),NANO);
	  Cf  = TWO*WallShearStress(Local_SolnBlk[nb].W[i][Local_SolnBlk[nb].JCl],
				    Local_SolnBlk[nb].Grid.Cell[i][Local_SolnBlk[nb].JCl].Xc,
				    Local_SolnBlk[nb].Grid.nodeSW(i,Local_SolnBlk[nb].JCl).X,
				    Local_SolnBlk[nb].Grid.nodeSE(i,Local_SolnBlk[nb].JCl).X,
				    -Local_SolnBlk[nb].Grid.nfaceS(i,Local_SolnBlk[nb].JCl))/(IP->Wo.rho*sqr(IP->Wo.v.x));
	  Cf2 = TWO*WallShearStress2(Local_SolnBlk[nb].Grid.xfaceS(i,Local_SolnBlk[nb].JCl),
				     Local_SolnBlk[nb].Grid.Cell[i][Local_SolnBlk[nb].JCl].Xc,
				     Local_SolnBlk[nb].W[i][Local_SolnBlk[nb].JCl],
				     Local_SolnBlk[nb].dWdx[i][Local_SolnBlk[nb].JCl],
				     Local_SolnBlk[nb].dWdy[i][Local_SolnBlk[nb].JCl],
				     -Local_SolnBlk[nb].Grid.nfaceS(i,Local_SolnBlk[nb].JCl))/(IP->Wo.rho*sqr(IP->Wo.v.x));
	  Out_File_Skin.setf(ios::scientific);
	  Out_File_Skin << " " << Local_SolnBlk[nb].Grid.Cell[i][Local_SolnBlk[nb].JCl].Xc.x
			<< " " << Rex
			<< " " << Cf
			<< " " << Cf2
			<< " " << Cfe
			<< " " << Cf - Cfe
			<< " " << Cf2 - Cfe
			<< " " << Local_SolnBlk[nb].Grid.Cell[i][Local_SolnBlk[nb].JCl].A
			<< endl;
	}
      }
    } else if (Interface_Component_List.Ni && numberofcells_cf) {
      Out_File_Skin << "ZONE T =  \"Block Number = " << nb
		    << "\" \\ \n"
		    << "I = " << numberofcells_cf << " \\ \n"
		    << "J = " << 1 << " \\ \n"
		    << "F = POINT \n";
      for (int j = Local_SolnBlk[nb].JCl; j <= Local_SolnBlk[nb].JCu; j++) {
	for (int i = Local_SolnBlk[nb].ICl; i <= Local_SolnBlk[nb].ICu; i++) {
	  if (Mesh[nb].cell_status[i][j] == CELL_STATUS_ACTIVE &&
	      (Mesh[nb].cell_status[i-1][j] != CELL_STATUS_ACTIVE ||
	       Mesh[nb].cell_status[i+1][j] != CELL_STATUS_ACTIVE ||
	       Mesh[nb].cell_status[i][j-1] != CELL_STATUS_ACTIVE ||
	       Mesh[nb].cell_status[i][j+1] != CELL_STATUS_ACTIVE) &&
	      Local_SolnBlk[nb].Grid.Cell[i][j].Xc.x >= TOLER &&
	      Local_SolnBlk[nb].Grid.Cell[i][j].Xc.x <= IP->Plate_Length) {
	    if (Mesh[nb].cell_status[i-1][j] != CELL_STATUS_ACTIVE &&
		Local_SolnBlk[nb].Grid.lfaceW(i,j) > TOLER) {
	      Rex = (IP->Wo.v.x/IP->Wo.nu())*Local_SolnBlk[nb].Grid.xfaceW(i,j).x;
	      Cfe = TWO*0.33206/max(sqrt(Rex),NANO);
	      Cf  = TWO*WallShearStress(Local_SolnBlk[nb].W[i][j],
					Local_SolnBlk[nb].Grid.Cell[i][j].Xc,
					Local_SolnBlk[nb].Grid.nodeNW(i,j).X,
					Local_SolnBlk[nb].Grid.nodeSW(i,j).X,
					-Local_SolnBlk[nb].Grid.nfaceW(i,j))/(IP->Wo.rho*sqr(IP->Wo.v.x));
	      Cf2 = TWO*WallShearStress2(Local_SolnBlk[nb].Grid.xfaceW(i,j),
					 Local_SolnBlk[nb].Grid.Cell[i][j].Xc,
					 Local_SolnBlk[nb].W[i][j],
					 Local_SolnBlk[nb].dWdx[i][j],
					 Local_SolnBlk[nb].dWdy[i][j],
					 -Local_SolnBlk[nb].Grid.nfaceW(i,j))/(IP->Wo.rho*sqr(IP->Wo.v.x));
	      Out_File_Skin.setf(ios::scientific);
	      Out_File_Skin << " " << Local_SolnBlk[nb].Grid.xfaceW(i,j).x//Cell[i][j].Xc.x
			    << " " << Rex
			    << " " << Cf
			    << " " << Cf2
			    << " " << Cfe
			    << " " << Cf - Cfe
			    << " " << Cf2 - Cfe
			    << " " << Local_SolnBlk[nb].Grid.Cell[i][j].A
			    << endl;
	    } else if (Mesh[nb].cell_status[i+1][j] != CELL_STATUS_ACTIVE &&
		       Local_SolnBlk[nb].Grid.lfaceE(i,j) > TOLER) {
	      Rex = (IP->Wo.v.x/IP->Wo.nu())*Local_SolnBlk[nb].Grid.xfaceE(i,j).x;
	      Cfe = TWO*0.33206/max(sqrt(Rex),NANO);
	      Cf  = TWO*WallShearStress(Local_SolnBlk[nb].W[i][j],
					Local_SolnBlk[nb].Grid.Cell[i][j].Xc,
					Local_SolnBlk[nb].Grid.nodeSE(i,j).X,
					Local_SolnBlk[nb].Grid.nodeNE(i,j).X,
					-Local_SolnBlk[nb].Grid.nfaceE(i,j))/(IP->Wo.rho*sqr(IP->Wo.v.x));
	      Cf2 = TWO*WallShearStress2(Local_SolnBlk[nb].Grid.xfaceE(i,j),
					 Local_SolnBlk[nb].Grid.Cell[i][j].Xc,
					 Local_SolnBlk[nb].W[i][j],
					 Local_SolnBlk[nb].dWdx[i][j],
					 Local_SolnBlk[nb].dWdy[i][j],
					 -Local_SolnBlk[nb].Grid.nfaceE(i,j))/(IP->Wo.rho*sqr(IP->Wo.v.x));
	      Out_File_Skin.setf(ios::scientific);
	      Out_File_Skin << " " << Local_SolnBlk[nb].Grid.xfaceE(i,j).x//Cell[i][j].Xc.x
			    << " " << Rex
			    << " " << Cf
			    << " " << Cf2
			    << " " << Cfe
			    << " " << Cf - Cfe
			    << " " << Cf2 - Cfe
			    << " " << Local_SolnBlk[nb].Grid.Cell[i][j].A
			    << endl;
	    } else if (Mesh[nb].cell_status[i][j-1] != CELL_STATUS_ACTIVE &&
		       Local_SolnBlk[nb].Grid.lfaceS(i,j) > TOLER) {
	      Rex = (IP->Wo.v.x/IP->Wo.nu())*Local_SolnBlk[nb].Grid.xfaceS(i,j).x;
	      Cfe = TWO*0.33206/max(sqrt(Rex),NANO);
	      Cf  = TWO*WallShearStress(Local_SolnBlk[nb].W[i][j],
					Local_SolnBlk[nb].Grid.Cell[i][j].Xc,
					Local_SolnBlk[nb].Grid.nodeSW(i,j).X,
					Local_SolnBlk[nb].Grid.nodeSE(i,j).X,
					-Local_SolnBlk[nb].Grid.nfaceS(i,j))/(IP->Wo.rho*sqr(IP->Wo.v.x));
	      Cf2 = TWO*WallShearStress2(Local_SolnBlk[nb].Grid.xfaceS(i,j),
					 Local_SolnBlk[nb].Grid.Cell[i][j].Xc,
					 Local_SolnBlk[nb].W[i][j],
					 Local_SolnBlk[nb].dWdx[i][j],
					 Local_SolnBlk[nb].dWdy[i][j],
					 -Local_SolnBlk[nb].Grid.nfaceS(i,j))/(IP->Wo.rho*sqr(IP->Wo.v.x));
	      Out_File_Skin.setf(ios::scientific);
	      Out_File_Skin << " " << Local_SolnBlk[nb].Grid.xfaceS(i,j).x//Cell[i][j].Xc.x
			    << " " << Rex
			    << " " << Cf
			    << " " << Cf2
			    << " " << Cfe
			    << " " << Cf - Cfe
			    << " " << Cf2 - Cfe
			    << " " << Local_SolnBlk[nb].Grid.Cell[i][j].A
			    << endl;
	    } else if (Mesh[nb].cell_status[i][j+1] != CELL_STATUS_ACTIVE &&
		       Local_SolnBlk[nb].Grid.lfaceN(i,j) > TOLER) {
	      Rex = (IP->Wo.v.x/IP->Wo.nu())*Local_SolnBlk[nb].Grid.xfaceN(i,j).x;
	      Cfe = TWO*0.33206/max(sqrt(Rex),NANO);
	      Cf  = TWO*WallShearStress(Local_SolnBlk[nb].W[i][j],
					Local_SolnBlk[nb].Grid.Cell[i][j].Xc,
					Local_SolnBlk[nb].Grid.nodeNE(i,j).X,
					Local_SolnBlk[nb].Grid.nodeNW(i,j).X,
					-Local_SolnBlk[nb].Grid.nfaceN(i,j))/(IP->Wo.rho*sqr(IP->Wo.v.x));
	      Cf2 = TWO*WallShearStress2(Local_SolnBlk[nb].Grid.xfaceN(i,j),
					 Local_SolnBlk[nb].Grid.Cell[i][j].Xc,
					 Local_SolnBlk[nb].W[i][j],
					 Local_SolnBlk[nb].dWdx[i][j],
					 Local_SolnBlk[nb].dWdy[i][j],
					 -Local_SolnBlk[nb].Grid.nfaceN(i,j))/(IP->Wo.rho*sqr(IP->Wo.v.x));
	      Out_File_Skin.setf(ios::scientific);
	      Out_File_Skin << " " << Local_SolnBlk[nb].Grid.xfaceN(i,j).x//Cell[i][j].Xc.x
			    << " " << Rex
			    << " " << Cf
			    << " " << Cf2
			    << " " << Cfe
			    << " " << Cf - Cfe
			    << " " << Cf2 - Cfe
			    << " " << Local_SolnBlk[nb].Grid.Cell[i][j].A
			    << endl;
	    }
	  }
	}
      }
    }
  }

}

/**********************************************************************
 * EmbeddedBoundaries2D::BCs_Interface --                             *
 **********************************************************************/
template <> inline int EmbeddedBoundaries2D<NavierStokes2D_cState,
					    NavierStokes2D_pState,
					    NavierStokes2D_Quad_Block,
					    NavierStokes2D_Input_Parameters>::
BCs_Interface(const int &nb, const double &Time) {

  // Exit immediately if none of the embedded boundaries are present in
  // the current solution block.
  if (!Adjustment_Data[nb].Interface_Present[0]) return 0;

  int error_flag, Ni, neighbour_flag;
  int Interface_BC_Type;
  double length;//, Temperature;
  Vector2D V;

  for (int j = Local_SolnBlk[nb].JCl-Local_SolnBlk[nb].Nghost; j <= Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost; j++) {
    for (int i = Local_SolnBlk[nb].ICl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].ICu+Local_SolnBlk[nb].Nghost; i++) {

      // Only apply embedded interface boundary data if the current cell
      // is not active (internal to one of the embedded interfaces) and
      // one of its neighbours is active.  Note that the inactive cell
      // status number corresponds to the interface number (1..Ni) the
      // cell is internal to.
      if (Mesh[nb].cell_status[i][j] != CELL_STATUS_ACTIVE) {

	// Store the interface number.
	Ni = Mesh[nb].cell_status[i][j];
	if (Ni < 1 || Ni > Interface_Union_List.Ni) return 4501;

	// Reset neighbour flag, primitive variable solution state, and
	// the total face length.
	neighbour_flag = OFF;
	Local_SolnBlk[nb].W[i][j].Vacuum();
	length = ZERO;

	// Apply boundary condition for any active neighbour cells if
	// required.

	// NORTH face.
	if (j < Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost) {
	  if (Mesh[nb].cell_status[i][j+1] == CELL_STATUS_ACTIVE &&
	      Local_SolnBlk[nb].Grid.lfaceN(i,j) > TOLER) {
	    length += Local_SolnBlk[nb].Grid.lfaceN(i,j);
	    // Determine the boundary condition type.
	    Interface_BC_Type = Interface_Union_List[Ni].Determine_Interface_BC_Type(Local_SolnBlk[nb].Grid.xfaceN(i,j));
	    // Determine the boundary velocity.
	    V = Interface_Union_List[Ni].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceN(i,j),Time);
	    // Apply the boundary condition.
	    switch(Interface_BC_Type) {
	    case INTERFACE_BC_REFLECTION :
	      Local_SolnBlk[nb].W[i][j] += Reflect(Local_SolnBlk[nb].W[i][j+1],
						   -Local_SolnBlk[nb].Grid.nfaceN(i,j),
						   V)*Local_SolnBlk[nb].Grid.lfaceN(i,j);
	      break;
	    case INTERFACE_BC_BURNING_SURFACE :
	      Local_SolnBlk[nb].W[i][j] += BurningSurface(Local_SolnBlk[nb].W[i][j+1],
							  -Local_SolnBlk[nb].Grid.nfaceN(i,j))*Local_SolnBlk[nb].Grid.lfaceN(i,j);
	      break;
	    case INTERFACE_BC_WALL_VISCOUS_HEATFLUX :
	      Local_SolnBlk[nb].W[i][j] += WallViscousHeatFlux(Local_SolnBlk[nb].W[i][j+1],
							       -Local_SolnBlk[nb].Grid.nfaceN(i,j),
							       V)*Local_SolnBlk[nb].Grid.lfaceN(i,j);
	      break;
	    case INTERFACE_BC_WALL_VISCOUS_ISOTHERMAL :
	      // Determine the boundary temperature.
	      //Temperature = Interface_Union_List[Ni].Determine_Interface_Temperature(Local_SolnBlk[nb].Grid.xfaceN(i,j));
	      Local_SolnBlk[nb].W[i][j] += WallViscousIsothermal(Local_SolnBlk[nb].W[i][j+1],
								 -Local_SolnBlk[nb].Grid.nfaceN(i,j),
								 V,IP->Twall)*Local_SolnBlk[nb].Grid.lfaceN(i,j);
	      break;
	    case INTERFACE_BC_RINGLEB :
	      Local_SolnBlk[nb].W[i][j] += RinglebFlowAverageState(Local_SolnBlk[nb].W[i][j],
								   Local_SolnBlk[nb].Grid.nodeSW(i,j).X,
								   Local_SolnBlk[nb].Grid.nodeSE(i,j).X,
								   Local_SolnBlk[nb].Grid.nodeNE(i,j).X,
								   Local_SolnBlk[nb].Grid.nodeNW(i,j).X)*Local_SolnBlk[nb].Grid.lfaceN(i,j);
	      break;
	    };
	    // Update neighbour_flag.
	    neighbour_flag = ON;
	  }
	}

	// SOUTH face.
	if (j > Local_SolnBlk[nb].JCl-Local_SolnBlk[nb].Nghost) {
	  if (Mesh[nb].cell_status[i][j-1] == CELL_STATUS_ACTIVE &&
	      Local_SolnBlk[nb].Grid.lfaceS(i,j) > TOLER) {
	    length += Local_SolnBlk[nb].Grid.lfaceS(i,j);
	    // Determine the boundary condition type.
	    Interface_BC_Type = Interface_Union_List[Ni].Determine_Interface_BC_Type(Local_SolnBlk[nb].Grid.xfaceS(i,j));
	    // Determine the boundary velocity.
	    V = Interface_Union_List[Ni].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceS(i,j),Time);	    
	    // Apply the boundary condition.
	    switch(Interface_BC_Type) {
	    case INTERFACE_BC_REFLECTION :
	      Local_SolnBlk[nb].W[i][j] += Reflect(Local_SolnBlk[nb].W[i][j-1],
						   -Local_SolnBlk[nb].Grid.nfaceS(i,j),
						   V)*Local_SolnBlk[nb].Grid.lfaceS(i,j);
	      break;
	    case INTERFACE_BC_BURNING_SURFACE :
	      Local_SolnBlk[nb].W[i][j] += BurningSurface(Local_SolnBlk[nb].W[i][j-1],
							  -Local_SolnBlk[nb].Grid.nfaceS(i,j))*Local_SolnBlk[nb].Grid.lfaceS(i,j);
	      break;
	    case INTERFACE_BC_WALL_VISCOUS_HEATFLUX :
	      Local_SolnBlk[nb].W[i][j] += WallViscousHeatFlux(Local_SolnBlk[nb].W[i][j-1],
							       -Local_SolnBlk[nb].Grid.nfaceS(i,j),
							       V)*Local_SolnBlk[nb].Grid.lfaceS(i,j);
	      break;
	    case INTERFACE_BC_WALL_VISCOUS_ISOTHERMAL :
	      // Determine the boundary temperature.
	      //Temperature = Interface_Union_List[Ni].Determine_Interface_Temperature(Local_SolnBlk[nb].Grid.xfaceS(i,j));
	      Local_SolnBlk[nb].W[i][j] += WallViscousIsothermal(Local_SolnBlk[nb].W[i][j-1],
								 -Local_SolnBlk[nb].Grid.nfaceS(i,j),
								 V,IP->Twall)*Local_SolnBlk[nb].Grid.lfaceS(i,j);
	      break;
	    case INTERFACE_BC_RINGLEB :
	      Local_SolnBlk[nb].W[i][j] += RinglebFlowAverageState(Local_SolnBlk[nb].W[i][j],
								   Local_SolnBlk[nb].Grid.nodeSW(i,j).X,
								   Local_SolnBlk[nb].Grid.nodeSE(i,j).X,
								   Local_SolnBlk[nb].Grid.nodeNE(i,j).X,
								   Local_SolnBlk[nb].Grid.nodeNW(i,j).X)*Local_SolnBlk[nb].Grid.lfaceS(i,j);
	      break;
	    default :
	      cout << endl << " -> " << Local_SolnBlk[nb].Grid.Cell[i][j].Xc << " " << Interface_BC_Type;
	      return 777;
	    };
	    // Update neighbour_flag.
	    neighbour_flag = ON;
	  }
	}

	// EAST face.
	if (i < Local_SolnBlk[nb].ICu+Local_SolnBlk[nb].Nghost) {
	  if (Mesh[nb].cell_status[i+1][j] == CELL_STATUS_ACTIVE &&
	      Local_SolnBlk[nb].Grid.lfaceE(i,j) > TOLER) {
	    length += Local_SolnBlk[nb].Grid.lfaceE(i,j);
	    // Determine the boundary condition type.
	    Interface_BC_Type = Interface_Union_List[Ni].Determine_Interface_BC_Type(Local_SolnBlk[nb].Grid.xfaceE(i,j));
	    // Determine the boundary velocity.
	    V = Interface_Union_List[Ni].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceE(i,j),Time);
	    // Apply the boundary condition.
	    switch(Interface_BC_Type) {
	    case INTERFACE_BC_REFLECTION :
	      Local_SolnBlk[nb].W[i][j] += Reflect(Local_SolnBlk[nb].W[i+1][j],
						   -Local_SolnBlk[nb].Grid.nfaceE(i,j),
						   V)*Local_SolnBlk[nb].Grid.lfaceE(i,j);
	      break;
	    case INTERFACE_BC_BURNING_SURFACE :
	      Local_SolnBlk[nb].W[i][j] += BurningSurface(Local_SolnBlk[nb].W[i+1][j],
							  -Local_SolnBlk[nb].Grid.nfaceE(i,j))*Local_SolnBlk[nb].Grid.lfaceE(i,j);
	      break;
	    case INTERFACE_BC_WALL_VISCOUS_HEATFLUX :
	      Local_SolnBlk[nb].W[i][j] += WallViscousHeatFlux(Local_SolnBlk[nb].W[i+1][j],
							       -Local_SolnBlk[nb].Grid.nfaceE(i,j),
							       V)*Local_SolnBlk[nb].Grid.lfaceE(i,j);
	      break;
	    case INTERFACE_BC_WALL_VISCOUS_ISOTHERMAL :
	      // Determine the boundary temperature.
	      //Temperature = Interface_Union_List[Ni].Determine_Interface_Temperature(Local_SolnBlk[nb].Grid.xfaceE(i,j));
	      Local_SolnBlk[nb].W[i][j] += WallViscousIsothermal(Local_SolnBlk[nb].W[i+1][j],
								 -Local_SolnBlk[nb].Grid.nfaceE(i,j),
								 V,IP->Twall)*Local_SolnBlk[nb].Grid.lfaceE(i,j);
	      break;
	    case INTERFACE_BC_RINGLEB :
	      Local_SolnBlk[nb].W[i][j] += RinglebFlowAverageState(Local_SolnBlk[nb].W[i][j],
								   Local_SolnBlk[nb].Grid.nodeSW(i,j).X,
								   Local_SolnBlk[nb].Grid.nodeSE(i,j).X,
								   Local_SolnBlk[nb].Grid.nodeNE(i,j).X,
								   Local_SolnBlk[nb].Grid.nodeNW(i,j).X)*Local_SolnBlk[nb].Grid.lfaceE(i,j);
	      break;
	    };
	    // Update neighbour_flag.
	    neighbour_flag = ON;
	  }
	}

	// WEST face.
	if (i > Local_SolnBlk[nb].ICl-Local_SolnBlk[nb].Nghost) {
	  if (Mesh[nb].cell_status[i-1][j] == CELL_STATUS_ACTIVE &&
	      Local_SolnBlk[nb].Grid.lfaceW(i,j) > TOLER) {
	    length += Local_SolnBlk[nb].Grid.lfaceW(i,j);
	    // Determine the boundary condition type.
	    Interface_BC_Type = Interface_Union_List[Ni].Determine_Interface_BC_Type(Local_SolnBlk[nb].Grid.xfaceW(i,j));
	    // Determine the boundary velocity.
	    V = Interface_Union_List[Ni].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceW(i,j),Time);
	    // Apply the boundary condition.
	    switch(Interface_BC_Type) {
	    case INTERFACE_BC_REFLECTION :
	      Local_SolnBlk[nb].W[i][j] += Reflect(Local_SolnBlk[nb].W[i-1][j],
						   -Local_SolnBlk[nb].Grid.nfaceW(i,j),
						   V)*Local_SolnBlk[nb].Grid.lfaceW(i,j);
	      break;
	    case INTERFACE_BC_BURNING_SURFACE :
	      Local_SolnBlk[nb].W[i][j] += BurningSurface(Local_SolnBlk[nb].W[i-1][j],
							  -Local_SolnBlk[nb].Grid.nfaceW(i,j))*Local_SolnBlk[nb].Grid.lfaceW(i,j);
	      break;
	    case INTERFACE_BC_WALL_VISCOUS_HEATFLUX :
	      Local_SolnBlk[nb].W[i][j] += WallViscousHeatFlux(Local_SolnBlk[nb].W[i-1][j],
							       -Local_SolnBlk[nb].Grid.nfaceW(i,j),
							       V)*Local_SolnBlk[nb].Grid.lfaceW(i,j);
	      break;
	    case INTERFACE_BC_WALL_VISCOUS_ISOTHERMAL :
	      // Determine the boundary temperature.
	      //Temperature = Interface_Union_List[Ni].Determine_Interface_Temperature(Local_SolnBlk[nb].Grid.xfaceW(i,j));
	      Local_SolnBlk[nb].W[i][j] += WallViscousIsothermal(Local_SolnBlk[nb].W[i-1][j],
								 -Local_SolnBlk[nb].Grid.nfaceW(i,j),
								 V,IP->Twall)*Local_SolnBlk[nb].Grid.lfaceW(i,j);
	      break;
	    case INTERFACE_BC_RINGLEB :
	      Local_SolnBlk[nb].W[i][j] += RinglebFlowAverageState(Local_SolnBlk[nb].W[i][j],
								   Local_SolnBlk[nb].Grid.nodeSW(i,j).X,
								   Local_SolnBlk[nb].Grid.nodeSE(i,j).X,
								   Local_SolnBlk[nb].Grid.nodeNE(i,j).X,
								   Local_SolnBlk[nb].Grid.nodeNW(i,j).X)*Local_SolnBlk[nb].Grid.lfaceW(i,j);
	      break;
	    };
	    // Update neighbour_flag.
	    neighbour_flag = ON;
	  }
	}

	// Apply boundary condition for any active distance-2 neighbour
	// cells if required (if any exist and there were no active
	// neighbour cells.
	if (!neighbour_flag) {

	  // NORTH face.
	  if (j <= Local_SolnBlk[nb].JCu) {
	    if (Mesh[nb].cell_status[i][j+2] == CELL_STATUS_ACTIVE &&
		Local_SolnBlk[nb].Grid.lfaceN(i,j+1) > TOLER) {
	      length += Local_SolnBlk[nb].Grid.lfaceN(i,j+1);
	      // Determine the boundary condition type.
	      Interface_BC_Type = Interface_Union_List[Ni].Determine_Interface_BC_Type(Local_SolnBlk[nb].Grid.xfaceN(i,j+1));
	      // Determine the boundary velocity.
	      V = Interface_Union_List[Ni].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceN(i,j+1),Time);
	      // Apply the boundary condition.
	      switch(Interface_BC_Type) {
	      case INTERFACE_BC_REFLECTION :
		Local_SolnBlk[nb].W[i][j] += Reflect(Local_SolnBlk[nb].W[i][j+2],
						     -Local_SolnBlk[nb].Grid.nfaceN(i,j+1),
						     V)*Local_SolnBlk[nb].Grid.lfaceN(i,j+1);
		break;
	      case INTERFACE_BC_BURNING_SURFACE :
		Local_SolnBlk[nb].W[i][j] += BurningSurface(Local_SolnBlk[nb].W[i][j+2],
							    -Local_SolnBlk[nb].Grid.nfaceN(i,j+1))*Local_SolnBlk[nb].Grid.lfaceN(i,j+1);
		break;
	      case INTERFACE_BC_WALL_VISCOUS_HEATFLUX :
		Local_SolnBlk[nb].W[i][j] += WallViscousHeatFlux(Local_SolnBlk[nb].W[i][j+2],
								 -Local_SolnBlk[nb].Grid.nfaceN(i,j+1),
								 V)*Local_SolnBlk[nb].Grid.lfaceN(i,j+1);
		break;
	      case INTERFACE_BC_WALL_VISCOUS_ISOTHERMAL :
		// Determine the boundary temperature.
		//Temperature = Interface_Union_List[Ni].Determine_Interface_Temperature(Local_SolnBlk[nb].Grid.xfaceN(i,j+1));
		Local_SolnBlk[nb].W[i][j] += WallViscousIsothermal(Local_SolnBlk[nb].W[i][j+2],
								   -Local_SolnBlk[nb].Grid.nfaceN(i,j+1),
								   V,IP->Twall)*Local_SolnBlk[nb].Grid.lfaceN(i,j+1);
		break;
	      case INTERFACE_BC_RINGLEB :
		Local_SolnBlk[nb].W[i][j] += RinglebFlowAverageState(Local_SolnBlk[nb].W[i][j],
								     Local_SolnBlk[nb].Grid.nodeSW(i,j).X,
								     Local_SolnBlk[nb].Grid.nodeSE(i,j).X,
								     Local_SolnBlk[nb].Grid.nodeNE(i,j).X,
								     Local_SolnBlk[nb].Grid.nodeNW(i,j).X)*Local_SolnBlk[nb].Grid.lfaceN(i,j+1);
		break;
	      };
	    }
	  }

	  // SOUTH face.
	  if (j >= Local_SolnBlk[nb].JCl) {
	    if (Mesh[nb].cell_status[i][j-2] == CELL_STATUS_ACTIVE &&
		Local_SolnBlk[nb].Grid.lfaceS(i,j-1) > TOLER) {
	      length += Local_SolnBlk[nb].Grid.lfaceS(i,j-1);
	      // Determine the boundary condition type.
	      Interface_BC_Type = Interface_Union_List[Ni].Determine_Interface_BC_Type(Local_SolnBlk[nb].Grid.xfaceS(i,j-1));
	      // Determine the boundary velocity.
	      V = Interface_Union_List[Ni].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceS(i,j-1),Time);
	      // Apply the boundary condition.
	      switch(Interface_BC_Type) {
	      case INTERFACE_BC_REFLECTION :
		Local_SolnBlk[nb].W[i][j] += Reflect(Local_SolnBlk[nb].W[i][j-2],
						     -Local_SolnBlk[nb].Grid.nfaceS(i,j-1),
						     V)*Local_SolnBlk[nb].Grid.lfaceS(i,j-1);
		break;
	      case INTERFACE_BC_BURNING_SURFACE :
		Local_SolnBlk[nb].W[i][j] += BurningSurface(Local_SolnBlk[nb].W[i][j-2],
							    -Local_SolnBlk[nb].Grid.nfaceS(i,j-1))*Local_SolnBlk[nb].Grid.lfaceS(i,j-1);
		break;
	      case INTERFACE_BC_WALL_VISCOUS_HEATFLUX :
		Local_SolnBlk[nb].W[i][j] += WallViscousHeatFlux(Local_SolnBlk[nb].W[i][j-2],
								 -Local_SolnBlk[nb].Grid.nfaceS(i,j-1),
								 V)*Local_SolnBlk[nb].Grid.lfaceS(i,j-1);
		break;
	      case INTERFACE_BC_WALL_VISCOUS_ISOTHERMAL :
		// Determine the boundary temperature.
		//Temperature = Interface_Union_List[Ni].Determine_Interface_Temperature(Local_SolnBlk[nb].Grid.xfaceS(i,j-1));
		Local_SolnBlk[nb].W[i][j] += WallViscousIsothermal(Local_SolnBlk[nb].W[i][j-2],
								   -Local_SolnBlk[nb].Grid.nfaceS(i,j-1),
								   V,IP->Twall)*Local_SolnBlk[nb].Grid.lfaceS(i,j-1);
		break;
	      case INTERFACE_BC_RINGLEB :
		Local_SolnBlk[nb].W[i][j] += RinglebFlowAverageState(Local_SolnBlk[nb].W[i][j],
								     Local_SolnBlk[nb].Grid.nodeSW(i,j).X,
								     Local_SolnBlk[nb].Grid.nodeSE(i,j).X,
								     Local_SolnBlk[nb].Grid.nodeNE(i,j).X,
								     Local_SolnBlk[nb].Grid.nodeNW(i,j).X)*Local_SolnBlk[nb].Grid.lfaceS(i,j-1);
		break;
	      };
	    }
	  }

	  // EAST face.
	  if (i <= Local_SolnBlk[nb].ICu) {
	    if (Mesh[nb].cell_status[i+2][j] == CELL_STATUS_ACTIVE &&
		Local_SolnBlk[nb].Grid.lfaceE(i+1,j) > TOLER) {
	      length += Local_SolnBlk[nb].Grid.lfaceE(i+1,j);
	      // Determine the boundary condition type.
	      Interface_BC_Type = Interface_Union_List[Ni].Determine_Interface_BC_Type(Local_SolnBlk[nb].Grid.xfaceE(i+1,j));
	      // Determine the boundary velocity.
	      V = Interface_Union_List[Ni].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceE(i+1,j),Time);
	      // Apply the boundary condition.
	      switch(Interface_BC_Type) {
	      case INTERFACE_BC_REFLECTION :
		Local_SolnBlk[nb].W[i][j] += Reflect(Local_SolnBlk[nb].W[i+2][j],
						     -Local_SolnBlk[nb].Grid.nfaceE(i+1,j),
						     V)*Local_SolnBlk[nb].Grid.lfaceE(i+1,j);
		break;
	      case INTERFACE_BC_BURNING_SURFACE :
		Local_SolnBlk[nb].W[i][j] += BurningSurface(Local_SolnBlk[nb].W[i+2][j],
							    -Local_SolnBlk[nb].Grid.nfaceE(i+1,j))*Local_SolnBlk[nb].Grid.lfaceE(i+1,j);
		break;
	      case INTERFACE_BC_WALL_VISCOUS_HEATFLUX :
		Local_SolnBlk[nb].W[i][j] += WallViscousHeatFlux(Local_SolnBlk[nb].W[i+2][j],
								 -Local_SolnBlk[nb].Grid.nfaceE(i+1,j),
								 V)*Local_SolnBlk[nb].Grid.lfaceE(i+1,j);
		break;
	      case INTERFACE_BC_WALL_VISCOUS_ISOTHERMAL :
		// Determine the boundary temperature.
		//Temperature = Interface_Union_List[Ni].Determine_Interface_Temperature(Local_SolnBlk[nb].Grid.xfaceE(i+1,j));
		Local_SolnBlk[nb].W[i][j] += WallViscousIsothermal(Local_SolnBlk[nb].W[i+2][j],
								   -Local_SolnBlk[nb].Grid.nfaceE(i+1,j),
								   V,IP->Twall)*Local_SolnBlk[nb].Grid.lfaceE(i+1,j);
		break;
	      case INTERFACE_BC_RINGLEB :
		Local_SolnBlk[nb].W[i][j] += RinglebFlowAverageState(Local_SolnBlk[nb].W[i][j],
							   Local_SolnBlk[nb].Grid.nodeSW(i,j).X,
							   Local_SolnBlk[nb].Grid.nodeSE(i,j).X,
							   Local_SolnBlk[nb].Grid.nodeNE(i,j).X,
							   Local_SolnBlk[nb].Grid.nodeNW(i,j).X)*Local_SolnBlk[nb].Grid.lfaceE(i+1,j);
		break;
	      };
	    }
	  }

	  // WEST face.
	  if (i >= Local_SolnBlk[nb].ICl) {
	    if (Mesh[nb].cell_status[i-2][j] == CELL_STATUS_ACTIVE &&
		Local_SolnBlk[nb].Grid.lfaceW(i-1,j) > TOLER) {
	      length += Local_SolnBlk[nb].Grid.lfaceW(i-1,j);
	      // Determine the boundary condition type.
	      Interface_BC_Type = Interface_Union_List[Ni].Determine_Interface_BC_Type(Local_SolnBlk[nb].Grid.xfaceW(i-1,j));
	      // Determine the boundary velocity.
	      V = Interface_Union_List[Ni].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceW(i-1,j),Time);
	      // Apply the boundary condition.
	      switch(Interface_BC_Type) {
	      case INTERFACE_BC_REFLECTION :
		Local_SolnBlk[nb].W[i][j] += Reflect(Local_SolnBlk[nb].W[i-2][j],
						     -Local_SolnBlk[nb].Grid.nfaceW(i-1,j),
						     V)*Local_SolnBlk[nb].Grid.lfaceW(i-1,j);
		break;
	      case INTERFACE_BC_BURNING_SURFACE :
		Local_SolnBlk[nb].W[i][j] += BurningSurface(Local_SolnBlk[nb].W[i-2][j],
							    -Local_SolnBlk[nb].Grid.nfaceW(i-1,j))*Local_SolnBlk[nb].Grid.lfaceW(i-1,j);
		break;
	      case INTERFACE_BC_WALL_VISCOUS_HEATFLUX :
		Local_SolnBlk[nb].W[i][j] += WallViscousHeatFlux(Local_SolnBlk[nb].W[i-2][j],
								 -Local_SolnBlk[nb].Grid.nfaceW(i-1,j),
								 V)*Local_SolnBlk[nb].Grid.lfaceW(i-1,j);
		break;
	      case INTERFACE_BC_WALL_VISCOUS_ISOTHERMAL :
		// Determine the boundary temperature.
		//Temperature = Interface_Union_List[Ni].Determine_Interface_Temperature(Local_SolnBlk[nb].Grid.xfaceW(i-1,j));
		Local_SolnBlk[nb].W[i][j] += WallViscousIsothermal(Local_SolnBlk[nb].W[i-2][j],
								   -Local_SolnBlk[nb].Grid.nfaceW(i-1,j),
								   V,IP->Twall)*Local_SolnBlk[nb].Grid.lfaceW(i-1,j);
		break;
	      case INTERFACE_BC_RINGLEB :
		Local_SolnBlk[nb].W[i][j] += RinglebFlowAverageState(Local_SolnBlk[nb].W[i][j],
								     Local_SolnBlk[nb].Grid.nodeSW(i,j).X,
								     Local_SolnBlk[nb].Grid.nodeSE(i,j).X,
								     Local_SolnBlk[nb].Grid.nodeNE(i,j).X,
								     Local_SolnBlk[nb].Grid.nodeNW(i,j).X)*Local_SolnBlk[nb].Grid.lfaceW(i-1,j);
		break;
	      };
	    }
	  }

	}

	// Determine the final (length-averaged) boundary state.
	if (length > ZERO) { Local_SolnBlk[nb].W[i][j] = Local_SolnBlk[nb].W[i][j]/length; }
	else { Local_SolnBlk[nb].W[i][j] = IP->Wo; Local_SolnBlk[nb].W[i][j].v = Vector2D_ZERO; }

	// Set conserved variable solution state.
	Local_SolnBlk[nb].U[i][j] = U(Local_SolnBlk[nb].W[i][j]);

      }

    }
  }

  // Embedded boundary boundary conditions applied successfully.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D:: CFL --                                      *
 **********************************************************************/
template <> inline double EmbeddedBoundaries2D<NavierStokes2D_cState,
					       NavierStokes2D_pState,
					       NavierStokes2D_Quad_Block,
					       NavierStokes2D_Input_Parameters>::
CFL(const double &Time) {

  double dtMin = MILLION, a, nu, v, V, dt_inv, dt_vis;

  // Determine the allowable time step for each solution block.
  for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
    if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      for (int j = Local_SolnBlk[nb].JCl-Local_SolnBlk[nb].Nghost; j <= Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost; j++) {
	for (int i = Local_SolnBlk[nb].ICl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].ICu+Local_SolnBlk[nb].Nghost; i++) {
	  if (i < Local_SolnBlk[nb].ICl || i > Local_SolnBlk[nb].ICu || j < Local_SolnBlk[nb].JCl || j > Local_SolnBlk[nb].JCu) {
	    Local_SolnBlk[nb].dt[i][j] = ZERO;
	  } else if (Mesh[nb].cell_status[i][j] != CELL_STATUS_ACTIVE) {
	    Local_SolnBlk[nb].dt[i][j] = ZERO;
	  } else {

	    // Initialize the time-step.
	    dt_inv = ZERO; dt_vis = ZERO;

	    // Determine the speed of sound in the current cell
	    // for the inviscid time-step calculation:
	    a = Local_SolnBlk[nb].W[i][j].c();

	    // Determine the viscosity in the current cell for the
	    // viscous time-step calculation.
	    if (Local_SolnBlk[nb].Flow_Type) nu = max(Local_SolnBlk[nb].W[i][j].nu(),Local_SolnBlk[nb].W[i][j].nuT());

	    // Determine the NORTH face contribution if required.
	    if (Local_SolnBlk[nb].Grid.lfaceN(i,j) > TOLER) {
	      // Inviscid time-step.
	      if (Mesh[nb].cell_status[i][j+1] == CELL_STATUS_ACTIVE) V = ZERO;
	      else V = Interface_Union_List[Mesh[nb].cell_status[i][j+1]].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceN(i,j),Time)*Local_SolnBlk[nb].Grid.nfaceN(i,j);
	      v = Local_SolnBlk[nb].W[i][j].v*Local_SolnBlk[nb].Grid.nfaceN(i,j);
	      dt_inv += (a+fabs(v-V)+fabs(V))*Local_SolnBlk[nb].Grid.lfaceN(i,j);
	      // Viscous time-step.
	      if (Local_SolnBlk[nb].Flow_Type) dt_vis += HALF*sqr(Local_SolnBlk[nb].Grid.lfaceN(i,j))/nu;
	    }

	    // Determine the SOUTH face contribution if required.
	    if (Local_SolnBlk[nb].Grid.lfaceS(i,j) > TOLER) {
	      // Inviscid time-step.
	      if (Mesh[nb].cell_status[i][j-1] == CELL_STATUS_ACTIVE) V = ZERO;
	      else V = Interface_Union_List[Mesh[nb].cell_status[i][j-1]].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceS(i,j),Time)*Local_SolnBlk[nb].Grid.nfaceS(i,j);
	      v = Local_SolnBlk[nb].W[i][j].v*Local_SolnBlk[nb].Grid.nfaceS(i,j);
	      dt_inv += (a+fabs(v-V)+fabs(V))*Local_SolnBlk[nb].Grid.lfaceS(i,j);
	      // Viscous time-step.
	      if (Local_SolnBlk[nb].Flow_Type) dt_vis += HALF*sqr(Local_SolnBlk[nb].Grid.lfaceS(i,j))/nu;
	    }

	    // Determine the EAST face contribution if required.
	    if (Local_SolnBlk[nb].Grid.lfaceE(i,j) > TOLER) {
	      // Inviscid time-step.
	      if (Mesh[nb].cell_status[i+1][j] == CELL_STATUS_ACTIVE) V = ZERO;
	      else V = Interface_Union_List[Mesh[nb].cell_status[i+1][j]].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceE(i,j),Time)*Local_SolnBlk[nb].Grid.nfaceE(i,j);
	      v = Local_SolnBlk[nb].W[i][j].v*Local_SolnBlk[nb].Grid.nfaceE(i,j);
	      dt_inv += (a+fabs(v-V)+fabs(V))*Local_SolnBlk[nb].Grid.lfaceE(i,j);
	      // Viscous time-step.
	      if (Local_SolnBlk[nb].Flow_Type) dt_vis += HALF*sqr(Local_SolnBlk[nb].Grid.lfaceE(i,j))/nu;
	    }

	    // Determine the WEST face contribution if required.
	    if (Local_SolnBlk[nb].Grid.lfaceW(i,j) > TOLER) {
	      // Inviscid time-step.
	      if (Mesh[nb].cell_status[i-1][j] == CELL_STATUS_ACTIVE) V = ZERO;
	      else V = Interface_Union_List[Mesh[nb].cell_status[i-1][j]].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceW(i,j),Time)*Local_SolnBlk[nb].Grid.nfaceW(i,j);
	      v = Local_SolnBlk[nb].W[i][j].v*Local_SolnBlk[nb].Grid.nfaceW(i,j);
	      dt_inv += (a+fabs(v-V)+fabs(V))*Local_SolnBlk[nb].Grid.lfaceW(i,j);
	      // Viscous time-step.
	      if (Local_SolnBlk[nb].Flow_Type) dt_vis += HALF*sqr(Local_SolnBlk[nb].Grid.lfaceW(i,j))/nu;
	    }

	    // Determine the minimum time-step for the cell.
	    Local_SolnBlk[nb].dt[i][j] = Local_SolnBlk[nb].Grid.Cell[i][j].A/dt_inv;
	    if (Local_SolnBlk[nb].Flow_Type) Local_SolnBlk[nb].dt[i][j] = min(Local_SolnBlk[nb].dt[i][j],
									      Local_SolnBlk[nb].Grid.Cell[i][j].A/dt_vis);

	    // Determine the mimimum global time-step.
	    dtMin = min(dtMin,Local_SolnBlk[nb].dt[i][j]);

	  }
	}
      }

      for (int j = Local_SolnBlk[nb].JCl-Local_SolnBlk[nb].Nghost; j <= Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost; j++) {
	for (int i = Local_SolnBlk[nb].ICl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].ICu+Local_SolnBlk[nb].Nghost; i++) {
	  if ((i < Local_SolnBlk[nb].ICl || i > Local_SolnBlk[nb].ICu ||
	       j < Local_SolnBlk[nb].JCl || j > Local_SolnBlk[nb].JCu) ||
	      Mesh[nb].cell_status[i][j] != CELL_STATUS_ACTIVE) {
	    Local_SolnBlk[nb].dt[i][j] = dtMin;
	  }
	}
      }

    }
  }

  // Return the global time step.
  return dtMin;

}

/**********************************************************************
 * EmbeddedBoundaries2D::dUdt_Residual_Evaluation --                  *
 *                                                                    *
 * This routine evaluates the residual for the specified solution     *
 * block using a 2nd-order limited upwind finite-volume spatial       *
 * discretization scheme.  The residual is stored in dUdt[][][0].     *
 *                                                                    *
 **********************************************************************/
template <> int EmbeddedBoundaries2D<NavierStokes2D_cState,
				     NavierStokes2D_pState,
				     NavierStokes2D_Quad_Block,
				     NavierStokes2D_Input_Parameters>::
dUdt_Residual_Evaluation(const double &Time) {

  Vector2D dX;
  NavierStokes2D_pState Wl, Wr;
  NavierStokes2D_cState Flux;

  NavierStokes2D_pState Wu, Wd, dWdxl, dWdyl, dWdxr, dWdyr;
  Vector2D Xl, Xr, Xu, Xd;

  int Interface_BC_Type, Ni;
  Vector2D V, Vu, Vd;
  //double Temperature;

  double q1, q2, k;

  // Evaluate the time rate of change of the solution (i.e., the
  // solution residuals) using a second-order limited upwind scheme
  // with a variety of flux functions.

  // Evaluate the stage solution residual for each solution block.
  for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
    if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {

      // Perform the linear reconstruction within each cell of the
      // computational grid for this stage.
      Linear_Least_Squares_Reconstruction(nb);

      // Add i-direction (zeta-direction) fluxes.
      for (int j = Local_SolnBlk[nb].JCl-1; j <= Local_SolnBlk[nb].JCu+1; j++) {
	Local_SolnBlk[nb].dUdt[Local_SolnBlk[nb].ICl-1][j][0].Vacuum();

	for (int i = Local_SolnBlk[nb].ICl-1; i <= Local_SolnBlk[nb].ICu; i++) {

	  Local_SolnBlk[nb].dUdt[i+1][j][0].Vacuum();

	  if (j >= Local_SolnBlk[nb].JCl && j <= Local_SolnBlk[nb].JCu) {

	    // Evaluate the cell interface i-direction INVISCID fluxes.
	    if (Local_SolnBlk[nb].Grid.lfaceE(i,j) < TOLER) {

	      // EAST face of cell (i,j) has zero length.
	      Wl.Standard_Atmosphere(); Wr.Standard_Atmosphere(); V = Vector2D_ZERO;

	    } else if (Mesh[nb].cell_status[i  ][j] != CELL_STATUS_ACTIVE &&
		       Mesh[nb].cell_status[i+1][j] != CELL_STATUS_ACTIVE) {

	      // EAST face of cell (i,j) is inactive.
	      Wl.Standard_Atmosphere(); Wr.Standard_Atmosphere(); V = Vector2D_ZERO;

	    } else if (Mesh[nb].cell_status[i  ][j] != CELL_STATUS_ACTIVE &&
		       Mesh[nb].cell_status[i+1][j] == CELL_STATUS_ACTIVE) {

	      // WEST face of cell (i+1,j) corresponds to an embedded boundary.
	      dX = Local_SolnBlk[nb].Grid.xfaceW(i+1,j)-Local_SolnBlk[nb].Grid.Cell[i+1][j].Xc;
	      Wr = Local_SolnBlk[nb].W[i+1][j] + (Local_SolnBlk[nb].phi[i+1][j]^Local_SolnBlk[nb].dWdx[i+1][j])*dX.x +
		                                 (Local_SolnBlk[nb].phi[i+1][j]^Local_SolnBlk[nb].dWdy[i+1][j])*dX.y;
	      Ni = Mesh[nb].cell_status[i][j];
	      // Determine the boundary condition type.
	      Interface_BC_Type = Interface_Union_List[Ni].Determine_Interface_BC_Type(Local_SolnBlk[nb].Grid.xfaceW(i+1,j));
	      // Determine the boundary velocity.
	      V = Interface_Union_List[Ni].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceW(i+1,j),Time);
	      // Determine the left state by applying the appropriate
	      // boundary condition.
	      if (Interface_BC_Type == INTERFACE_BC_REFLECTION) {
		// WEST face of cell (i+1,j) is a REFLECTION boundary.
		Wl = Reflect(Wr,Local_SolnBlk[nb].Grid.nfaceW(i+1,j),V);
	      } else if (Interface_BC_Type == INTERFACE_BC_BURNING_SURFACE) {
		// WEST face of cell (i+1,j) is a BURNING_SURFACE boundary.
		Wl = BurningSurface(Wr,Local_SolnBlk[nb].Grid.nfaceW(i+1,j));
	      } else if (Interface_BC_Type == INTERFACE_BC_WALL_VISCOUS_HEATFLUX) {
		// WEST face of cell (i+1,j) is a WALL_VISCOUS_HEATFLUX boundary.
		Wl = WallViscousHeatFlux(Wr,Local_SolnBlk[nb].Grid.nfaceW(i+1,j),V);
	      } else if (Interface_BC_Type == INTERFACE_BC_WALL_VISCOUS_ISOTHERMAL) {
		// WEST face of cell (i+1,j) is a WALL_VISCOUS_ISOTHERMAL boundary.
		//Temperature = Interface_Union_List[Ni].Determine_Interface_Temperature(Local_SolnBlk[nb].Grid.xfaceW(i+1,j));
		Wl = WallViscousIsothermal(Wr,Local_SolnBlk[nb].Grid.nfaceW(i+1,j),V,IP->Twall);
	      } else if (Interface_BC_Type == INTERFACE_BC_RINGLEB) {
		// WEST face of cell (i+1,j) is a RINGLEB_FLOW boundary.
		Wl = RinglebFlow(Wr,Local_SolnBlk[nb].Grid.Cell[i][j].Xc,q1,k);
		Wl = RinglebFlow(Wr,Local_SolnBlk[nb].Grid.Cell[i+1][j].Xc,q2,k);
		Wl = RinglebFlow(Wr,Local_SolnBlk[nb].Grid.xfaceW(i+1,j));
		if (q1 < IP->Isotach_Line && q2 > IP->Isotach_Line)
		  Wl = BC_Characteristic_Pressure(Wr,Wl,Local_SolnBlk[nb].Grid.nfaceW(i+1,j));
	      }

	    } else if (Mesh[nb].cell_status[i  ][j] == CELL_STATUS_ACTIVE &&
		       Mesh[nb].cell_status[i+1][j] != CELL_STATUS_ACTIVE) {

	      // EAST face of cell (i,j) corresponds to an embedded boundary.
	      dX = Local_SolnBlk[nb].Grid.xfaceE(i,j) - Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	      Wl = Local_SolnBlk[nb].W[i][j] + (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdx[i][j])*dX.x +
	                                       (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdy[i][j])*dX.y;
	      Ni = Mesh[nb].cell_status[i+1][j];
	      // Determine the boundary condition type.
	      Interface_BC_Type = Interface_Union_List[Ni].Determine_Interface_BC_Type(Local_SolnBlk[nb].Grid.xfaceE(i,j));
	      // Determine the boundary velocity.
	      V = Interface_Union_List[Ni].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceE(i,j),Time);
	      // Determine the right state by applying the appropriate
	      // boundary condition.
	      if (Interface_BC_Type == INTERFACE_BC_REFLECTION) {
		// EAST face of cell (i,j) is a REFLECTION boundary.
		Wr = Reflect(Wl,Local_SolnBlk[nb].Grid.nfaceE(i,j),V);
	      } else if (Interface_BC_Type == INTERFACE_BC_BURNING_SURFACE) {
		// EAST face of cell (i,j) is a BURNING_SURFACE boundary.
		Wr = BurningSurface(Wl,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      } else if (Interface_BC_Type == INTERFACE_BC_WALL_VISCOUS_HEATFLUX) {
		// EAST face of cell (i,j) is a WALL_VISCOUS_HEATFLUX boundary.
		Wr = WallViscousHeatFlux(Wl,Local_SolnBlk[nb].Grid.nfaceE(i,j),V);
	      } else if (Interface_BC_Type == INTERFACE_BC_WALL_VISCOUS_ISOTHERMAL) {
		// EAST face of cell (i,j) is a WALL_VISCOUS_ISOTHERMAL boundary.
		//Temperature = Interface_Union_List[Ni].Determine_Interface_Temperature(Local_SolnBlk[nb].Grid.xfaceE(i,j));
		Wr = WallViscousIsothermal(Wl,Local_SolnBlk[nb].Grid.nfaceE(i,j),V,IP->Twall);
	      } else if (Interface_BC_Type == INTERFACE_BC_RINGLEB) {
		// EAST face of cell (i,j) is a RINGLEB FLOW boundary.
		Wr = RinglebFlow(Wr,Local_SolnBlk[nb].Grid.Cell[i+1][j].Xc,q1,k);
		Wr = RinglebFlow(Wr,Local_SolnBlk[nb].Grid.Cell[i][j].Xc,q2,k);
		Wr = RinglebFlow(Wr,Local_SolnBlk[nb].Grid.xfaceE(i,j));
		if (q1 < IP->Isotach_Line && q2 > IP->Isotach_Line)
		  Wr = BC_Characteristic_Pressure(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      }

	    } else if (Mesh[nb].cell_status[i  ][j] == CELL_STATUS_ACTIVE &&
		       Mesh[nb].cell_status[i+1][j] == CELL_STATUS_ACTIVE) {

	      V = Vector2D_ZERO;

	      if (i == Local_SolnBlk[nb].ICl-1 && 
		  (Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_REFLECTION ||
		   Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_BURNING_SURFACE ||
		   Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX ||
		   Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
		   Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_MOVING_WALL_HEATFLUX ||
		   Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL ||
		   Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_RINGLEB_FLOW ||
		   Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_CHARACTERISTIC)) {

		dX = Local_SolnBlk[nb].Grid.xfaceW(i+1,j) - Local_SolnBlk[nb].Grid.Cell[i+1][j].Xc;
		Wr = Local_SolnBlk[nb].W[i+1][j] + (Local_SolnBlk[nb].phi[i+1][j]^Local_SolnBlk[nb].dWdx[i+1][j])*dX.x +
	                                           (Local_SolnBlk[nb].phi[i+1][j]^Local_SolnBlk[nb].dWdy[i+1][j])*dX.y;

		// WEST face of cell (i+1,j) is a normal boundary.
		if (Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_REFLECTION) {
		  // WEST face of cell (i+1,j) is a REFLECTION boundary.
		  Wl = Reflect(Wr,Local_SolnBlk[nb].Grid.nfaceW(i+1,j));
		} else if (Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_BURNING_SURFACE) {
		  // WEST face of cell (i+1,j) is a BURNING_SURFACE boundary.
		  Wl = BurningSurface(Wr,Local_SolnBlk[nb].Grid.nfaceW(i+1,j));
		} else if (Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX) {
		  // WEST face of cell (i+1,j) is a WALL_VISCOUS_HEATFLUX boundary.
		  Wl = WallViscousHeatFlux(Wr,Local_SolnBlk[nb].Grid.nfaceW(i+1,j));
		} else if (Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL) {
		  // WEST face of cell (i+1,j) is a WALL_VISCOUS_ISOTHERMAL boundary.
		  Wl = WallViscousIsothermal(Wr,Local_SolnBlk[nb].Grid.nfaceW(i+1,j),Local_SolnBlk[nb].Twall);
		} else if (Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_MOVING_WALL) {
		  // WEST face of cell (i+1,j) is a MOVINGWALL_HEATFLUX boundary.
		  Wl = MovingWallHeatFlux(Wr,Local_SolnBlk[nb].Grid.nfaceW(i+1,j),Local_SolnBlk[nb].Vwall.x);
		} else if (Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL) {
		  // WEST face of cell (i+1,j) is a MOVINGWALL_ISOTHERMAL boundary.
		  Wl = MovingWallIsothermal(Wr,Local_SolnBlk[nb].Grid.nfaceW(i+1,j),Local_SolnBlk[nb].Vwall.x,Local_SolnBlk[nb].Twall);
		} else if (Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_RINGLEB_FLOW) {
		  // WEST face of cell (i+1,j) is a RINGLEB_FLOW boundary.
		  Wl = RinglebFlow(Wl,Local_SolnBlk[nb].Grid.xfaceW(i+1,j));
		} else {
		  // WEST face of cell (i+1,j) is a CHARACTERISTIC boundary.
		  Wl = BC_Characteristic_Pressure(Wr,Local_SolnBlk[nb].WoW[j],Local_SolnBlk[nb].Grid.nfaceW(i+1,j));
		}

	      } else if (i == Local_SolnBlk[nb].ICu &&
			 (Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_REFLECTION ||
			  Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_BURNING_SURFACE ||
			  Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX ||
			  Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
			  Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_MOVING_WALL_HEATFLUX ||
			  Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL ||
			  Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_RINGLEB_FLOW ||
			  Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_CHARACTERISTIC)) {

		dX = Local_SolnBlk[nb].Grid.xfaceE(i,j) - Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
		Wl = Local_SolnBlk[nb].W[i][j] + (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdx[i][j])*dX.x +
	                                         (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdy[i][j])*dX.y;

		// EAST face of cell (i,j) is a normal boundary.
		if (Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_REFLECTION) {
		  // EAST face of cell (i,j) is a REFLECTION boundary.
		  Wr = Reflect(Wl,Local_SolnBlk[nb].Grid.nfaceE(i,j));
		} else if (Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_BURNING_SURFACE) {
		  // EAST face of cell (i,j) is a BURNING_SURFACE boundary.
		  Wr = BurningSurface(Wl,Local_SolnBlk[nb].Grid.nfaceE(i,j));
		} else if (Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX) {
		  // EAST face of cell (i,j) is a WALL_VISCOUS_HEATFLUX boundary.
		  Wr = WallViscousHeatFlux(Wl,Local_SolnBlk[nb].Grid.nfaceE(i,j));
		} else if (Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL) {
		  // EAST face of cell (i,j) is a WALL_VISCOUS_ISOTHERMAL boundary.
		  Wr = WallViscousIsothermal(Wl,Local_SolnBlk[nb].Grid.nfaceE(i,j),Local_SolnBlk[nb].Twall);
		} else if (Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_MOVING_WALL_HEATFLUX) {
		  // EAST face of cell (i,j) is a MOVINGWALL_HEATFLUX boundary.
		  Wr = MovingWallHeatFlux(Wl,Local_SolnBlk[nb].Grid.nfaceE(i,j),Local_SolnBlk[nb].Vwall.x);
		} else if (Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL) {
		  // EAST face of cell (i,j) is a MOVINGWALL_ISOTHERMAL boundary.
		  Wr = MovingWallIsothermal(Wl,Local_SolnBlk[nb].Grid.nfaceE(i,j),Local_SolnBlk[nb].Vwall.x,Local_SolnBlk[nb].Twall);
		} else if (Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_RINGLEB_FLOW) {
		  // EAST face of cell (i,j) is a RINGLEB_FLOW boundary.
		  Wr = RinglebFlow(Wr,Local_SolnBlk[nb].Grid.xfaceE(i,j));
		} else {
		  // EAST face of cell (i,j) is a CHARACTERISTIC boundary.
		  Wr = BC_Characteristic_Pressure(Wl,Local_SolnBlk[nb].WoE[j],Local_SolnBlk[nb].Grid.nfaceE(i,j));
		}

	      } else {

		// EAST face is either a normal cell or possibly a FIXED, 
		// NONE or EXTRAPOLATION boundary.
		dX = Local_SolnBlk[nb].Grid.xfaceE(i  ,j) - Local_SolnBlk[nb].Grid.Cell[i  ][j].Xc;
		Wl = Local_SolnBlk[nb].W[i  ][j] + (Local_SolnBlk[nb].phi[i  ][j]^Local_SolnBlk[nb].dWdx[i  ][j])*dX.x +
	                                           (Local_SolnBlk[nb].phi[i  ][j]^Local_SolnBlk[nb].dWdy[i  ][j])*dX.y;
		dX = Local_SolnBlk[nb].Grid.xfaceW(i+1,j) - Local_SolnBlk[nb].Grid.Cell[i+1][j].Xc;
		Wr = Local_SolnBlk[nb].W[i+1][j] + (Local_SolnBlk[nb].phi[i+1][j]^Local_SolnBlk[nb].dWdx[i+1][j])*dX.x +
	                                           (Local_SolnBlk[nb].phi[i+1][j]^Local_SolnBlk[nb].dWdy[i+1][j])*dX.y;

	      }

	    }

	    // Determine EAST face INVISCID flux.
	    switch(IP->i_Flux_Function) {
	    case FLUX_FUNCTION_GODUNOV :
	      Flux = FluxGodunov_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_ROE :
	      Flux = FluxRoe_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_RUSANOV :
	      Flux = FluxRusanov_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_HLLE :
	      Flux = FluxHLLE_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_HLLL :
	      Flux = FluxHLLL_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_HLLC :
	      Flux = FluxHLLC_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_VANLEER :
	      Flux = FluxVanLeer_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_AUSM :
	      Flux = FluxAUSM_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_AUSMplus :
	      Flux = FluxAUSMplus_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_GODUNOV_MB :
	      Flux = FluxGodunov_MB_n(Wl,Wr,V,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_ROE_MB :
	      Flux = FluxRoe_MB_n(Wl,Wr,V,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_HLLE_MB :
	      Flux = FluxHLLE_MB_n(Wl,Wr,V,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_VANLEER_MB :
	      Flux = FluxVanLeer_MB_n(Wl,Wr,V,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    default:
	      Flux = FluxRoe_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    };

	    // Compute the cell centred stress tensor and heat flux vector if required.
 	    if (Local_SolnBlk[nb].Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
 		(Local_SolnBlk[nb].Flow_Type && Local_SolnBlk[nb].Axisymmetric)) {
 	      Local_SolnBlk[nb].W[i][j].ComputeViscousTerms(Local_SolnBlk[nb].dWdx[i][j],
 							    Local_SolnBlk[nb].dWdy[i][j],
 							    Local_SolnBlk[nb].Grid.Cell[i][j].Xc,
 							    Local_SolnBlk[nb].Axisymmetric);
 	      Local_SolnBlk[nb].U[i][j].tau = Local_SolnBlk[nb].W[i][j].tau;
 	      Local_SolnBlk[nb].U[i][j].q = Local_SolnBlk[nb].W[i][j].q;
 	    }

	    // Determine the i-direction viscous flux.

			do { // with while (0) to allow jumps to the end of the viscous flux calculation.
	    	if (Local_SolnBlk[nb].Flow_Type == 0) { break; }

				if (Local_SolnBlk[nb].Grid.lfaceE(i,j) < TOLER) { break; }

				if (Mesh[nb].cell_status[i  ][j] != CELL_STATUS_ACTIVE &&
						Mesh[nb].cell_status[i+1][j] != CELL_STATUS_ACTIVE) {
					break;
				}

				int adiabatic_flag = OFF;
				Vector2D Xface = Local_SolnBlk[nb].Grid.xfaceE(i, j);
				Vector2D norm_dir = Local_SolnBlk[nb].Grid.nfaceE(i, j);

	      if (Mesh[nb].cell_status[i  ][j] != CELL_STATUS_ACTIVE &&
				    Mesh[nb].cell_status[i+1][j] == CELL_STATUS_ACTIVE) {

					// WEST face of cell (i+1,j) corresponds to an embedded boundary.
					Xr = Local_SolnBlk[nb].Grid.Cell[i+1][j].Xc; 
					Wr = Local_SolnBlk[nb].W[i+1][j];

					Ni = Mesh[nb].cell_status[i][j];
					Interface_BC_Type = Interface_Union_List[Ni].Determine_Interface_BC_Type(Xface);
					
					// No viscous flux for a reflected boundary. Is this true?
					if (Interface_BC_Type == INTERFACE_BC_REFLECTION) { break; } 

					switch (Interface_BC_Type) {
						// The face velocity is the same regardless of the BC type 
						// and is set below.
						case INTERFACE_BC_WALL_VISCOUS_HEATFLUX:
							Wface.rho = Wr.rho; 
							Wface.p = Wr.p; Wface.k = Wr.k; Wface.omega = Wr.omega;
							break;
						case INTERFACE_BC_WALL_VISCOUS_ISOTHERMAL:
							Wface.rho = Wr.p/(Wr.R*Local_SolnBlk[nb].Twall); 
							Wface.p = Wr.p; Wface.k = Wr.k; Wface.omega = Wr.omega;
							break;
						case INTERFACE_BC_BURNING_SURFACE:
							Wface = BurningSurface(Wr, norm_dir);
							break;
					}

					Wface.v = Interface_Union_List[Ni].Determine_Interface_Velocity(Xface, Time);

					if (Interface_BC_Type == INTERFACE_BC_WALL_VISCOUS_HEATFLUX) {
						adiabatic_flag = ON;
					}

					switch(IP->i_Viscous_Reconstruction) {
						case VISCOUS_RECONSTRUCTION_CARTESIAN :
						case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
							{
							int stencil_flag = DIAMONDPATH_RIGHT_TRIANGLE;
							Xd = Local_SolnBlk[nb].Grid.Node[i+1][j  ].X;
							Wd = Wface;
							Wd.v = Interface_Union_List[Ni].Determine_Interface_Velocity(Xd, Time);

							Xu = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
							Wu = Wface;
							Wu.v = Interface_Union_List[Ni].Determine_Interface_Velocity(Xu, Time);

							DiamondPath_Find_dWdX(dWdx, dWdy,
									BlankVector, BlankState, 
									Xd, Wd,
									Xr, Wr,
									Xu, Wu,
									stencil_flag);
							}
							break;
						case VISCOUS_RECONSTRUCTION_HYBRID :
							dWdxr = Local_SolnBlk[nb].dWdx[i+1][j]; 
							dWdyr = Local_SolnBlk[nb].dWdy[i+1][j];

							Xl = Xface;
							Wl = Wface;
							dWdxl = dWdxr; dWdyl = dWdyr; // Hmm. 

							Hybrid_Find_dWdX(dWdx, dWdy, 
									Xl, Wl, dWdxl, dWdyl,
									Xr, Wr, dWdxr, dWdyr,
									norm_dir);
							break;
						case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
							dWdx = Local_SolnBlk[nb].dWdx[i+1][j]; 
							dWdy = Local_SolnBlk[nb].dWdy[i+1][j];
							break;
					}
				} else if (Mesh[nb].cell_status[i  ][j] == CELL_STATUS_ACTIVE &&
				           Mesh[nb].cell_status[i+1][j] != CELL_STATUS_ACTIVE) {

					// EAST face of cell (i,j) corresponds to an embedded boundary.
					Xl = Local_SolnBlk[nb].Grid.Cell[i][j].Xc; 
					Wl = Local_SolnBlk[nb].W[i][j];

					Ni = Mesh[nb].cell_status[i+1][j];
					Interface_BC_Type = Interface_Union_List[Ni].Determine_Interface_BC_Type(Xface);

					// No viscous flux for a reflected boundary. Is this true?
					if (Interface_BC_Type == INTERFACE_BC_REFLECTION) { break; } 

					switch (Interface_BC_Type) {
						// The face velocity is the same regardless of the BC type 
						// and is set below.
						case INTERFACE_BC_WALL_VISCOUS_HEATFLUX:
							Wface.rho = Wl.rho; 
							Wface.p = Wl.p; Wface.k = Wl.k; Wface.omega = Wl.omega;
							break;
						case INTERFACE_BC_WALL_VISCOUS_ISOTHERMAL:
							Wface.rho = Wl.p/(Wl.R*IP->Twall); 
							Wface.p = Wl.p; Wface.k = Wl.k; Wface.omega = Wl.omega;
							break;
						case INTERFACE_BC_BURNING_SURFACE:
							Wface = BurningSurface(Wl, norm_dir);
							break;
					}

					Wface.v = Interface_Union_List[Ni].Determine_Interface_Velocity(Xface, Time);

					if (Interface_BC_Type == INTERFACE_BC_WALL_VISCOUS_HEATFLUX) {
						adiabatic_flag = ON;
					}

					switch(IP->i_Viscous_Reconstruction) {
						case VISCOUS_RECONSTRUCTION_CARTESIAN :
						case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
							{
							int stencil_flag = DIAMONDPATH_LEFT_TRIANGLE;
							Xd = Local_SolnBlk[nb].Grid.Node[i+1][j  ].X;
							Wd = Wface;
							Wd.v = Interface_Union_List[Ni].Determine_Interface_Velocity(Xd, Time);

							Xu = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
							Wu = Wface;
							Wu.v = Interface_Union_List[Ni].Determine_Interface_Velocity(Xu, Time);

							DiamondPath_Find_dWdX(dWdx, dWdy,
									Xl, Wl,
									Xd, Wd,
									BlankVector, BlankState,
									Xu, Wu,
									stencil_flag);
							}
							break;
						case VISCOUS_RECONSTRUCTION_HYBRID :
							dWdxl = Local_SolnBlk[nb].dWdx[i][j]; 
							dWdyl = Local_SolnBlk[nb].dWdy[i][j]; 

							Xr = Xface;
							Wr = Wface;
							dWdxr = dWdxl; dWdyr = dWdyl; // Hmm.

							Hybrid_Find_dWdX(dWdx, dWdy, 
									Xl, Wl, dWdxl, dWdyl,
									Xr, Wr, dWdxr, dWdyr,
									norm_dir);
							break;
						case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
							dWdx = Local_SolnBlk[nb].dWdx[i][j]; 
							dWdy = Local_SolnBlk[nb].dWdy[i][j];
							break;
					}

				} else if (Mesh[nb].cell_status[i  ][j] == CELL_STATUS_ACTIVE &&
				           Mesh[nb].cell_status[i+1][j] == CELL_STATUS_ACTIVE) {

					if (i == Local_SolnBlk[nb].ICl-1 && 
							(Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX ||
							 Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
							 Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_MOVING_WALL_HEATFLUX ||
							 Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL ||
							 Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_BURNING_SURFACE)) {

						// WEST face of cell (i+1,j) is a normal boundary.

						Xr = Local_SolnBlk[nb].Grid.Cell[i+1][j].Xc; 
						Wr = Local_SolnBlk[nb].W[i+1][j];

						if (Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX ||
								Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_MOVING_WALL_HEATFLUX) {
							adiabatic_flag = ON;
						}

						switch (Local_SolnBlk[nb].Grid.BCtypeW[j]) {
							case BC_WALL_VISCOUS_HEATFLUX:
								Wface.rho = Wr.rho; Wface.v.x = ZERO; Wface.v.y = ZERO; 
								Wface.p = Wr.p; Wface.k = Wr.k; Wface.omega =	Wr.omega;
								break;
							case BC_WALL_VISCOUS_ISOTHERMAL:
								Wface.rho = Wr.p/(Wr.R*Local_SolnBlk[nb].Twall); Wface.v.x = ZERO; Wface.v.y = ZERO; 
								Wface.p = Wr.p; Wface.k = Wr.k; Wface.omega = Wr.omega;
								break;
							case BC_MOVING_WALL_HEATFLUX:
								Wface.rho = Wr.rho; Wface.v.x = Local_SolnBlk[nb].Vwall.x; Wface.v.y = Local_SolnBlk[nb].Vwall.y; 
								Wface.p = Wr.p; Wface.k = Wr.k; Wface.omrga = Wr.omega;
								break;
							case BC_MOVING_WALL_ISOTHERMAL:
								Wface.rho = Wr.p/(Wr.R*Local_SolnBlk[nb].Twall); Wface.v.x = Local_SolnBlk[nb].Vwall.x; Wface.v.y = Local_SolnBlk[nb].Vwall.y; 
								Wface.p = Wr.p; Wface.k = Wr.k; Wface.omega = Wr.omega;
								break;
							case BC_BURNING_SURFACE:
								Wface = BurningSurface(Wr, norm_dir);
								break;
						}

						switch(IP->i_Viscous_Reconstruction) {
							case VISCOUS_RECONSTRUCTION_CARTESIAN :
							case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
								{
								int stencil_flag = DIAMONDPATH_RIGHT_TRIANGLE;
								Xu = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
								Xd = Local_SolnBlk[nb].Grid.Node[i+1][j  ].X; 
								DiamondPath_Find_dWdX(dWdx, dWdy,
										BlankVector, BlankState, 
										Xd, Wface,
										Xr, Wr,
										Xu, Wface,
										stencil_flag);
								}
								break;
							case VISCOUS_RECONSTRUCTION_HYBRID :
								dWdxr = Local_SolnBlk[nb].dWdx[i+1][j];
								dWdyr = Local_SolnBlk[nb].dWdy[i+1][j];

								Xl = Xface;
								Wl = Wface;
								dWdxl = dWdxr; dWdyl = dWdyr; // Hmm.

								Hybrid_Find_dWdX(dWdx, dWdy, 
										Xl, Wl, dWdxl, dWdyl,
										Xr, Wr, dWdxr, dWdyr,
										norm_dir);
								break;
							case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
								dWdx = Local_SolnBlk[nb].dWdx[i+1][j];
								dWdy = Local_SolnBlk[nb].dWdy[i+1][j];
								break;
						}

					} else if (i == Local_SolnBlk[nb].ICu &&
							(Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX ||
							 Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
							 Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_MOVING_WALL_HEATFLUX ||
							 Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL ||
							 Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_BURNING_SURFACE)) {
						
						// EAST face of cell (i,j) is a normal boundary.
						
						Xl = Local_SolnBlk[nb].Grid.Cell[i][j].Xc; 
						Wl = Local_SolnBlk[nb].W[i][j];

						if (Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX ||
						    Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_MOVING_WALL_HEATFLUX) {
							adiabatic_flag = ON;
						}

						switch (Local_SolnBlk[nb].Grid.BCtypeE[j]) {
							case BC_WALL_VISCOUS_HEATFLUX:
								Wface.rho = Wl.rho; Wface.v.x = ZERO; Wface.v.y = ZERO; 
								Wface.p = Wl.p; Wface.k = Wl.k; Wface.omega = Wl.omega;
								break;
							case BC_WALL_VISCOUS_ISOTHERMAL:
								Wface.rho = Wl.p/(Wl.R*Local_SolnBlk[nb].Twall); Wface.v.x = ZERO; Wface.v.y = ZERO; 
								Wface.p = Wl.p; Wface.k = Wl.k; Wface.omega = Wl.omega;
								break;
							case BC_MOVING_WALL_HEATFLUX:
								Wface.rho = Wl.rho; Wface.v.x = Local_SolnBlk[nb].Vwall.x; Wface.v.y = Local_SolnBlk[nb].Vwall.y; 
								Wface.p = Wl.p; Wface.k = Wl.k; Wface.omega = Wl.omega;
								break;
							case BC_MOVING_WALL_ISOTHERMAL:
								Wface.rho = Wl.p/(Wl.R*Local_SolnBlk[nb].Twall); Wface.v.x = Local_SolnBlk[nb].Vwall.x; Wface.v.y = Local_SolnBlk[nb].Vwall.y; 
								Wface.p = Wl.p; Wface.k = Wl.k; Wface.omega = Wl.omega;
								break;
							case BC_MOVING_WALL_ISOTHERMAL:
								Wface = BurningSurface(Wl, norm_dir);
								break;
						}

						switch(IP->i_Viscous_Reconstruction) {
							case VISCOUS_RECONSTRUCTION_CARTESIAN :
							case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
								{
								int stencil_flag = DIAMONDPATH_LEFT_TRIANGLE;
								Xu = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
								Xd = Local_SolnBlk[nb].Grid.Node[i+1][j  ].X; 
								DiamondPath_Find_dWdX(dWdx, dWdy,
										Xl, Wl,
										Xd, Wface,
										BlankVector, BlankState,
										Xu, Wface,
										stencil_flag);
								}
								break;
							case VISCOUS_RECONSTRUCTION_HYBRID :
								dWdxl = Local_SolnBlk[nb].dWdx[i][j]; 
								dWdyl = Local_SolnBlk[nb].dWdy[i][j]; 

								Xr = Xface; 
								Wr = Wface;
								dWdxr = dWdxl; dWdyr = dWdyl; 

								Hybrid_Find_dWdX(dWdx, dWdy, 
										Xl, Wl, dWdxl, dWdyl,
										Xr, Wr, dWdxr, dWdyr,
										norm_dir);
								break;
							case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
								dWdx = Local_SolnBlk[nb].dWdx[i][j]; 
								dWdy = Local_SolnBlk[nb].dWdy[i][j]; 
								break;
						}

					} else {
						// EAST face is either a normal cell or possibly a non-viscous
						// boundary condition.
						Xl = Local_SolnBlk[nb].Grid.Cell[i  ][j].Xc; Wl = Local_SolnBlk[nb].W[i  ][j];
						Xr = Local_SolnBlk[nb].Grid.Cell[i+1][j].Xc; Wr = Local_SolnBlk[nb].W[i+1][j];

						if (IP->i_Viscous_Reconstruction == VISCOUS_RECONSTRUCTION_CARTESIAN ||
						    IP->i_Viscous_Reconstruction == VISCOUS_RECONSTRUCTION_DIAMOND_PATH) {

							int stencil_flag = DIAMONDPATH_QUADRILATERAL_RECONSTRUCTION;

							Xu = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X; 
							Xd = Local_SolnBlk[nb].Grid.Node[i+1][j  ].X; 

							// These two (WnNE and WnSE) are what make diamond path so expensive.
							Wu = Local_SolnBlk[nb].WnNE(i,j);
							Wd = Local_SolnBlk[nb].WnSE(i,j);

							DiamondPath_Find_dWdX(dWdx, dWdy,
									Xl, Wl, Xd, Wd, Xr, Wr, Xu, Wu, 
									stencil_flag);

							// Find Wface:
							error_flag = Bilinear_Interpolation(Wl, Xl, Wu, Xu, Wr, Xr, Wd, Xd,
									Xface, Wface);
							// and if error_flag?

						} else {
							Wface = HALF*(Wl + Wr);

							dWdxl = Local_SolnBlk[nb].dWdx[i][j]; dWdxr = Local_SolnBlk[nb].dWdy[i+1][j];
							dWdyl = Local_SolnBlk[nb].dWdy[i][j]; dWdyr = Local_SolnBlk[nb].dWdy[i+1][j];

							if (IP->i_Viscous_Reconstruction == VISCOUS_RECONSTRUCTION_HYBRID) {
								Hybrid_Find_dWdX(dWdx, dWdy, 
										Xl, Wl, dWdxl, dWdyl,
										Xr, Wr, dWdxr, dWdyr,
										norm_dir);
							} else if (IP->i_Viscous_Reconstruction == VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE) {
								dWdx = HALF*(dWdxl + dWdxr);
								dWdy = HALF*(dWdyl + dWdyr);
							}
						}
					}
				}
				// Compute the EAST face viscous flux.
				Flux -= ViscousFlux_n(Xface, Wface, dWdx, dWdy,
						norm_dir, Local_SolnBlk[nb].Axisymmetric, adiabatic_flag);
			} while (0); // end of viscous flux clause.

	    // Evaluate cell-averaged solution changes.
	    Local_SolnBlk[nb].dUdt[i  ][j][0] -= Flux*Local_SolnBlk[nb].Grid.lfaceE(i,j)/Local_SolnBlk[nb].Grid.Cell[i][j].A;
	    Local_SolnBlk[nb].dUdt[i+1][j][0] += Flux*Local_SolnBlk[nb].Grid.lfaceW(i+1,j)/Local_SolnBlk[nb].Grid.Cell[i+1][j].A;

	    // Include all required source terms.
	    Local_SolnBlk[nb].dUdt[i][j][0] += Local_SolnBlk[nb].W[i][j].S(Local_SolnBlk[nb].Grid.Cell[i][j].Xc,
									   Local_SolnBlk[nb].dWdx[i][j],
									   Local_SolnBlk[nb].dWdy[i][j],
									   Local_SolnBlk[nb].Axisymmetric);

	    // Include area change source term.
	    if (Interface_Union_List.Ni)
	      Local_SolnBlk[nb].dUdt[i][j][0] -= Local_SolnBlk[nb].U[i][j]*dAdt(nb,i,j,Time);

	    // Save west and east face boundary flux.
	    if (i == Local_SolnBlk[nb].ICl-1) {
	      Local_SolnBlk[nb].FluxW[j] = -Flux*Local_SolnBlk[nb].Grid.lfaceW(i+1,j);
	    } else if (i == Local_SolnBlk[nb].ICu) {
	      Local_SolnBlk[nb].FluxE[j] =  Flux*Local_SolnBlk[nb].Grid.lfaceE(i,j);
	    }

	  }
	}

	if (j > Local_SolnBlk[nb].JCl-1 && j < Local_SolnBlk[nb].JCu+1) {
	  Local_SolnBlk[nb].dUdt[Local_SolnBlk[nb].ICl-1][j][0].Vacuum();
	  Local_SolnBlk[nb].dUdt[Local_SolnBlk[nb].ICu+1][j][0].Vacuum();
	}

      }

      // Add j-direction (eta-direction) fluxes.
      for (int i = Local_SolnBlk[nb].ICl; i <= Local_SolnBlk[nb].ICu; i++) {
 	for (int j = Local_SolnBlk[nb].JCl-1; j <= Local_SolnBlk[nb].JCu; j++) {

 	  // Evaluate the cell interface j-direction INVISCID fluxes.
 	  if (Local_SolnBlk[nb].Grid.lfaceN(i,j) < TOLER) {
 	    // NORTH face of cell (i,j) has zero length.
 	    Wl.Standard_Atmosphere(); Wr.Standard_Atmosphere(); V = Vector2D_ZERO;

 	  } else if (Mesh[nb].cell_status[i][j  ] != CELL_STATUS_ACTIVE &&
 		     Mesh[nb].cell_status[i][j+1] != CELL_STATUS_ACTIVE) {
 	    // NORTH face of cell (i,j) is inactive.
 	    Wl.Standard_Atmosphere(); Wr.Standard_Atmosphere(); V = Vector2D_ZERO;

 	  } else if (Mesh[nb].cell_status[i][j  ] != CELL_STATUS_ACTIVE &&
 		     Mesh[nb].cell_status[i][j+1] == CELL_STATUS_ACTIVE) {
	    // SOUTH face of cell (i,j+1) corresponds to an embedded boundary.
	    dX = Local_SolnBlk[nb].Grid.xfaceS(i,j+1) - Local_SolnBlk[nb].Grid.Cell[i][j+1].Xc;
	    Wr = Local_SolnBlk[nb].W[i][j+1] + (Local_SolnBlk[nb].phi[i][j+1]^Local_SolnBlk[nb].dWdx[i][j+1])*dX.x +
	                                       (Local_SolnBlk[nb].phi[i][j+1]^Local_SolnBlk[nb].dWdy[i][j+1])*dX.y;
	    Ni = Mesh[nb].cell_status[i][j];
	    // Determine the boundary condition type.
	    Interface_BC_Type = Interface_Union_List[Ni].Determine_Interface_BC_Type(Local_SolnBlk[nb].Grid.xfaceS(i,j+1));
	    // Determine the boundary velocity.
	    V = Interface_Union_List[Ni].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceS(i,j+1),Time);
	    // Determine the left state by applying the appropriate
	    // boundary condition.
	    if (Interface_BC_Type == INTERFACE_BC_REFLECTION) {
	      // SOUTH face of cell (i,j+1) is a REFLECTION boundary.
	      Wl = Reflect(Wr,Local_SolnBlk[nb].Grid.nfaceS(i,j+1),V);
	    } else if (Interface_BC_Type == INTERFACE_BC_BURNING_SURFACE) {
	      // SOUTH face of cell (i,j+1) is a BURNING_SURFACE boundary.
	      Wl = BurningSurface(Wr,Local_SolnBlk[nb].Grid.nfaceS(i,j+1));
	    } else if (Interface_BC_Type == INTERFACE_BC_WALL_VISCOUS_HEATFLUX) {
	      // SOUTH face of cell (i,j+1) is a WALL_VISCOUS_HEATFLUX boundary.
	      Wl = WallViscousHeatFlux(Wr,Local_SolnBlk[nb].Grid.nfaceS(i,j+1),V);
	    } else if (Interface_BC_Type == INTERFACE_BC_WALL_VISCOUS_ISOTHERMAL) {
	      // SOUTH face of cell (i,j+1) is a WALL_VISCOUS_ISOTHERMAL boundary.
	      Wl = WallViscousIsothermal(Wr,Local_SolnBlk[nb].Grid.nfaceS(i,j+1),V,Local_SolnBlk[nb].Twall);
	    } else if (Interface_BC_Type == INTERFACE_BC_RINGLEB) {
	      // SOUTH face of cell (i,j+1) is a RINGLEB_FLOW boundary.
	      Wl = RinglebFlow(Wl,Local_SolnBlk[nb].Grid.Cell[i][j].Xc,q1,k);
	      Wl = RinglebFlow(Wl,Local_SolnBlk[nb].Grid.Cell[i][j+1].Xc,q2,k);
	      Wl = RinglebFlow(Wl,Local_SolnBlk[nb].Grid.xfaceS(i,j+1));
	      if (q1 < IP->Isotach_Line && q2 > IP->Isotach_Line)
		Wl = BC_Characteristic_Pressure(Wr,Wl,Local_SolnBlk[nb].Grid.nfaceS(i,j+1));
	    }

 	  } else if (Mesh[nb].cell_status[i][j  ] == CELL_STATUS_ACTIVE &&
 		     Mesh[nb].cell_status[i][j+1] != CELL_STATUS_ACTIVE) {
	    // NORTH face of cell (i,j) corresponds to an embedded boundary.
	    dX = Local_SolnBlk[nb].Grid.xfaceN(i,j) - Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	    Wl = Local_SolnBlk[nb].W[i][j] + (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdx[i][j])*dX.x +
	                                     (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdy[i][j])*dX.y;
	    Ni = Mesh[nb].cell_status[i][j+1];
	    // Determine the boundary condition type.
	    Interface_BC_Type = Interface_Union_List[Ni].Determine_Interface_BC_Type(Local_SolnBlk[nb].Grid.xfaceN(i,j));
	    // Determine the boundary velocity.
	    V = Interface_Union_List[Ni].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceN(i,j),Time);
	    // Determine the right state by applying the appropriate
	    // boundary condition.
	    if (Interface_BC_Type == INTERFACE_BC_REFLECTION) {
	      // NORTH face of cell (i,j) is a REFLECTION boundary.
	      Wr = Reflect(Wl,Local_SolnBlk[nb].Grid.nfaceN(i,j),V);
	    } else if (Interface_BC_Type == INTERFACE_BC_BURNING_SURFACE) {
	      // NORTH face of cell (i,j) is a BURNING_SURFACE boundary.
	      Wr = BurningSurface(Wl,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    } else if (Interface_BC_Type == INTERFACE_BC_WALL_VISCOUS_HEATFLUX) {
	      // NORTH face of cell (i,j) is a WALL_VISCOUS_HEATFLUX boundary.
	      Wr = WallViscousHeatFlux(Wl,Local_SolnBlk[nb].Grid.nfaceN(i,j),V);
	    } else if (Interface_BC_Type == INTERFACE_BC_WALL_VISCOUS_ISOTHERMAL) {
	      // NORTH face of cell (i,j) is a WALL_VISCOUS_ISOTHERMAL boundary.
	      Wr = WallViscousIsothermal(Wl,Local_SolnBlk[nb].Grid.nfaceN(i,j),V,Local_SolnBlk[nb].Twall);
	    } else if (Interface_BC_Type == INTERFACE_BC_RINGLEB) {
	      // NORTH face of cell (i,j) is a RINGLEB_FLOW boundary.
	      Wr = RinglebFlow(Wr,Local_SolnBlk[nb].Grid.Cell[i][j+1].Xc,q1,k);
	      Wr = RinglebFlow(Wr,Local_SolnBlk[nb].Grid.Cell[i][j].Xc,q2,k);
	      Wr = RinglebFlow(Wr,Local_SolnBlk[nb].Grid.xfaceN(i,j));
	      if (q1 < IP->Isotach_Line && q2 > IP->Isotach_Line)
		Wr = BC_Characteristic_Pressure(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    }

 	  } else if (Mesh[nb].cell_status[i][j  ] == CELL_STATUS_ACTIVE &&
 		     Mesh[nb].cell_status[i][j+1] == CELL_STATUS_ACTIVE) {

 	    V = Vector2D_ZERO;

 	    if (j == Local_SolnBlk[nb].JCl-1 && 
 		(Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_REFLECTION ||
 		 Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_BURNING_SURFACE ||
 		 Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX ||
 		 Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
 		 Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_MOVING_WALL_HEATFLUX ||
 		 Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_MOVING_WALL_ISOTHERMAL ||
 		 Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_RINGLEB_FLOW ||
 		 Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_CHARACTERISTIC)) {

	      dX = Local_SolnBlk[nb].Grid.xfaceS(i,j+1) - Local_SolnBlk[nb].Grid.Cell[i][j+1].Xc;
	      Wr = Local_SolnBlk[nb].W[i][j+1] + (Local_SolnBlk[nb].phi[i][j+1]^Local_SolnBlk[nb].dWdx[i][j+1])*dX.x +
   	                                         (Local_SolnBlk[nb].phi[i][j+1]^Local_SolnBlk[nb].dWdy[i][j+1])*dX.y;

	      // SOUTH face of cell (i,j+1) is a normal boundary.
	      if (Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_REFLECTION) {
		// SOUTH face of cell (i,j+1) is a REFLECTION boundary.
		Wl = Reflect(Wr,Local_SolnBlk[nb].Grid.nfaceS(i,j+1));
	      } else if (Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_BURNING_SURFACE) {
		// SOUTH face of cell (i,j+1) is a BURNING_SURFACE boundary.
		Wl = BurningSurface(Wr,Local_SolnBlk[nb].Grid.nfaceS(i,j+1));
	      } else if (Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX) {
		// SOUTH face of cell (i,j+1) is a WALL_VISCOUS_HEATFLUX boundary.
		Wl = WallViscousHeatFlux(Wr,Local_SolnBlk[nb].Grid.nfaceS(i,j+1));
	      } else if (Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL) {
		// SOUTH face of cell (i,j+1) is a WALL_VISCOUS_ISOTHERMAL boundary.
		Wl = WallViscousIsothermal(Wr,Local_SolnBlk[nb].Grid.nfaceS(i,j+1),Local_SolnBlk[nb].Twall);
	      } else if (Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_MOVING_WALL_HEATFLUX) {
		// SOUTH face of cell (i,j+1) is a MOVINGWALL_HEATFLUX boundary.
		Wl = MovingWallHeatFlux(Wr,Local_SolnBlk[nb].Grid.nfaceS(i,j+1),Local_SolnBlk[nb].Vwall.x);
	      } else if (Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_MOVING_WALL_ISOTHERMAL) {
		// SOUTH face of cell (i,j+1) is a MOVINGWALL_ISOTHERMAL boundary.
		Wl = MovingWallIsothermal(Wr,Local_SolnBlk[nb].Grid.nfaceS(i,j+1),Local_SolnBlk[nb].Vwall.x,Local_SolnBlk[nb].Twall);
	      } else if (Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_RINGLEB_FLOW) {
		// SOUTH face of cell (i,j+1) is a RINGLEB_FLOW boundary.
		Wl = RinglebFlow(Wl,Local_SolnBlk[nb].Grid.xfaceS(i,j+1));
		Wl = BC_Characteristic_Pressure(Wr,Wl,Local_SolnBlk[nb].Grid.nfaceS(i,j+1));
	      } else if (Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_CHARACTERISTIC) {
		// SOUTH face of cell (i,j+1) is a CHARACTERISTIC boundary.
		Wl = BC_Characteristic_Pressure(Wr,Local_SolnBlk[nb].WoS[i],Local_SolnBlk[nb].Grid.nfaceS(i,j+1));
	      }

	    } else if (j == Local_SolnBlk[nb].JCu && 
		       (Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_REFLECTION ||
			Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_BURNING_SURFACE ||
			Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX ||
			Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
			Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_MOVING_WALL_HEATFLUX ||
			Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_MOVING_WALL_ISOTHERMAL ||
			Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_RINGLEB_FLOW ||
			Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_CHARACTERISTIC)) {

	      dX = Local_SolnBlk[nb].Grid.xfaceN(i,j) - Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	      Wl = Local_SolnBlk[nb].W[i][j] + (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdx[i][j])*dX.x +
   	                                       (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdy[i][j])*dX.y;

	      // NORTH face of cell (i,j) is a normal boundary.
	      if (Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_REFLECTION) {
		// NORTH face of cell (i,j) is a REFLECTION boundary.
		Wr = Reflect(Wl,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	      } else if (Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_BURNING_SURFACE) {
		// NORTH face of cell (i,j) is a BURNING_SURFACE boundary.
		Wr = BurningSurface(Wl,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	      } else if (Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX) {
		// NORTH face of cell (i,j) is a WALL_VISCOUS_HEATFLUX boundary.
		Wr = WallViscousHeatFlux(Wl,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	      } else if (Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL) {
		// NORTH face of cell (i,j) is a WALL_VISCOUS_ISOTHERMAL boundary.
		Wr = WallViscousIsothermal(Wl,Local_SolnBlk[nb].Grid.nfaceN(i,j),Local_SolnBlk[nb].Twall);
	      } else if (Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_MOVING_WALL_HEATFLUX) {
		// NORTH face of cell (i,j) is a MOVINGWALL_HEATFLUX boundary.
		Wr = MovingWallHeatFlux(Wl,Local_SolnBlk[nb].Grid.nfaceN(i,j),Local_SolnBlk[nb].Vwall.x);
	      } else if (Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_MOVING_WALL_ISOTHERMAL) {
		// NORTH face of cell (i,j) is a MOVINGWALL_ISOTHERMAL boundary.
		Wr = MovingWallIsothermal(Wl,Local_SolnBlk[nb].Grid.nfaceN(i,j),Local_SolnBlk[nb].Vwall.x,Local_SolnBlk[nb].Twall);
	      } else if (Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_RINGLEB_FLOW) {
		// NORTH face of cell (i,j) is a RINGLEB_FLOW boundary.
		Wr = RinglebFlow(Wr,Local_SolnBlk[nb].Grid.xfaceN(i,j));
		Wr = BC_Characteristic_Pressure(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	      } else {
		// NORTH face of cell (i,j) is a CHARACTERISTIC boundary.
		Wr = BC_Characteristic_Pressure(Wl,Local_SolnBlk[nb].WoN[i],Local_SolnBlk[nb].Grid.nfaceN(i,j));
	      }

	    } else {

	      // NORTH face is either a normal cell or possibly a FIXED, 
	      // NONE or EXTRAPOLATION boundary.
	      dX = Local_SolnBlk[nb].Grid.xfaceN(i,j  ) - Local_SolnBlk[nb].Grid.Cell[i][j  ].Xc;
	      Wl = Local_SolnBlk[nb].W[i][j  ] + (Local_SolnBlk[nb].phi[i][j  ]^Local_SolnBlk[nb].dWdx[i][j  ])*dX.x +
                                                 (Local_SolnBlk[nb].phi[i][j  ]^Local_SolnBlk[nb].dWdy[i][j  ])*dX.y;
	      dX = Local_SolnBlk[nb].Grid.xfaceS(i,j+1) - Local_SolnBlk[nb].Grid.Cell[i][j+1].Xc;
	      Wr = Local_SolnBlk[nb].W[i][j+1] + (Local_SolnBlk[nb].phi[i][j+1]^Local_SolnBlk[nb].dWdx[i][j+1])*dX.x +
                                                 (Local_SolnBlk[nb].phi[i][j+1]^Local_SolnBlk[nb].dWdy[i][j+1])*dX.y;

	    }

	  }

	  // Determine NORTH face inviscid flux.
	  switch(IP->i_Flux_Function) {
	  case FLUX_FUNCTION_GODUNOV :
	    Flux = FluxGodunov_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_ROE :
	    Flux = FluxRoe_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_RUSANOV :
	    Flux = FluxRusanov_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_HLLE :
	    Flux = FluxHLLE_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_HLLL :
	    Flux = FluxHLLL_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_HLLC :
	    Flux = FluxHLLC_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_VANLEER :
	    Flux = FluxVanLeer_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_AUSM :
	    Flux = FluxAUSM_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_AUSMplus :
	    Flux = FluxAUSMplus_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_GODUNOV_MB :
	    Flux = FluxGodunov_MB_n(Wl,Wr,V,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_ROE_MB :
	    Flux = FluxRoe_MB_n(Wl,Wr,V,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_HLLE_MB :
	    Flux = FluxHLLE_MB_n(Wl,Wr,V,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_VANLEER_MB :
	    Flux = FluxVanLeer_MB_n(Wl,Wr,V,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  default:
	    Flux = FluxRoe_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  };

	    // Determine the j-direction viscous flux.

			do { // with while (0) to allow jumps to the end of the viscous flux calculation.
	    	if (Local_SolnBlk[nb].Flow_Type == 0) { break; }

				if (Local_SolnBlk[nb].Grid.lfaceN(i,j) < TOLER) { break; }

				if (Mesh[nb].cell_status[i][j  ] != CELL_STATUS_ACTIVE &&
						Mesh[nb].cell_status[i][j+1] != CELL_STATUS_ACTIVE) {
					break;
				}

				// Please rotate your head 90 degrees to the left before
				// attempting to read this code. For example, Xr is the
				// position of the cell to the north of the face in
				// question. This is necessary, for example, for
				// DiamondPath_Find_dWdX, which only understands "left" and
				// "right" triangles.

				int adiabatic_flag = OFF;
				Vector2D Xface = Local_SolnBlk[nb].Grid.xfaceN(i, j);
				Vector2D norm_dir = Local_SolnBlk[nb].Grid.nfaceN(i, j);

	      if (Mesh[nb].cell_status[i][j  ] != CELL_STATUS_ACTIVE &&
				    Mesh[nb].cell_status[i][j+1] == CELL_STATUS_ACTIVE) {

					// SOUTH face of cell (i,j+1) corresponds to an embedded boundary.
					Xr = Local_SolnBlk[nb].Grid.Cell[i][j+1].Xc; 
					Wr = Local_SolnBlk[nb].W[i][j+1];

					Ni = Mesh[nb].cell_status[i][j];
					Interface_BC_Type = Interface_Union_List[Ni].Determine_Interface_BC_Type(Xface);
					
					// No viscous flux for a reflected boundary. Is this true?
					if (Interface_BC_Type == INTERFACE_BC_REFLECTION) { break; } 

					switch (Interface_BC_Type) {
						// The face velocity is the same regardless of the BC type 
						// and is set below.
						case INTERFACE_BC_WALL_VISCOUS_HEATFLUX:
							Wface.rho = Wr.rho; 
							Wface.p = Wr.p; Wface.k = Wr.k; Wface.omega = Wr.omega;
							break;
						case INTERFACE_BC_WALL_VISCOUS_ISOTHERMAL:
							Wface.rho = Wr.p/(Wr.R*Local_SolnBlk[nb].Twall); 
							Wface.p = Wr.p; Wface.k = Wr.k; Wface.omega = Wr.omega;
							break;
						case INTERFACE_BC_BURNING_SURFACE:
							Wface = BurningSurface(Wr, norm_dir);
							break;
					}

					Wface.v = Interface_Union_List[Ni].Determine_Interface_Velocity(Xface, Time);

					if (Interface_BC_Type == INTERFACE_BC_WALL_VISCOUS_HEATFLUX) {
						adiabatic_flag = ON;
					}

					switch(IP->i_Viscous_Reconstruction) {
						case VISCOUS_RECONSTRUCTION_CARTESIAN :
						case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
							{
							int stencil_flag = DIAMONDPATH_RIGHT_TRIANGLE;
							Xd = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
							Wd = Wface;
							Wd.v = Interface_Union_List[Ni].Determine_Interface_Velocity(Xd, Time);

							Xu = Local_SolnBlk[nb].Grid.Node[i  ][j+1].X;
							Wu = Wface;
							Wu.v = Interface_Union_List[Ni].Determine_Interface_Velocity(Xu, Time);

							DiamondPath_Find_dWdX(dWdx, dWdy,
									BlankVector, BlankState, 
									Xd, Wd,
									Xr, Wr,
									Xu, Wu,
									stencil_flag);
							}
							break;
						case VISCOUS_RECONSTRUCTION_HYBRID :
							dWdxr = Local_SolnBlk[nb].dWdx[i][j+1]; 
							dWdyr = Local_SolnBlk[nb].dWdy[i][j+1];

							Xl = Xface;
							Wl = Wface;
							dWdxl = dWdxr; dWdyl = dWdyr; // Hmm. 

							Hybrid_Find_dWdX(dWdx, dWdy, 
									Xl, Wl, dWdxl, dWdyl,
									Xr, Wr, dWdxr, dWdyr,
									norm_dir);
							break;
						case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
							dWdx = Local_SolnBlk[nb].dWdx[i][j+1]; 
							dWdy = Local_SolnBlk[nb].dWdy[i][j+1];
							break;
					}
				} else if (Mesh[nb].cell_status[i][j  ] == CELL_STATUS_ACTIVE &&
				           Mesh[nb].cell_status[i][j+1] != CELL_STATUS_ACTIVE) {

					// NORTH face of cell (i,j) corresponds to an embedded boundary.
					Xl = Local_SolnBlk[nb].Grid.Cell[i][j].Xc; 
					Wl = Local_SolnBlk[nb].W[i][j];

					Ni = Mesh[nb].cell_status[i][j+1];
					Interface_BC_Type = Interface_Union_List[Ni].Determine_Interface_BC_Type(Xface);

					// No viscous flux for a reflected boundary. Is this true?
					if (Interface_BC_Type == INTERFACE_BC_REFLECTION) { break; } 

					switch (Interface_BC_Type) {
						// The face velocity is the same regardless of the BC type 
						// and is set below.
						case INTERFACE_BC_WALL_VISCOUS_HEATFLUX:
							Wface.rho = Wl.rho; 
							Wface.p = Wl.p; Wface.k = Wl.k; Wface.omega = Wl.omega;
							break;
						case INTERFACE_BC_WALL_VISCOUS_ISOTHERMAL:
							Wface.rho = Wl.p/(Wl.R*IP->Twall); 
							Wface.p = Wl.p; Wface.k = Wl.k; Wface.omega = Wl.omega;
							break;
						case INTERFACE_BC_BURNING_SURFACE:
							Wface = BurningSurface(Wl, norm_dir);
							break;
					}

					Wface.v = Interface_Union_List[Ni].Determine_Interface_Velocity(Xface, Time);

					if (Interface_BC_Type == INTERFACE_BC_WALL_VISCOUS_HEATFLUX) {
						adiabatic_flag = ON;
					}

					switch(IP->i_Viscous_Reconstruction) {
						case VISCOUS_RECONSTRUCTION_CARTESIAN :
						case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
							{
							int stencil_flag = DIAMONDPATH_LEFT_TRIANGLE;
							Xd = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
							Wd = Wface;
							Wd.v = Interface_Union_List[Ni].Determine_Interface_Velocity(Xd, Time);

							Xu = Local_SolnBlk[nb].Grid.Node[i  ][j+1].X;
							Wu = Wface;
							Wu.v = Interface_Union_List[Ni].Determine_Interface_Velocity(Xu, Time);

							DiamondPath_Find_dWdX(dWdx, dWdy,
									Xl, Wl,
									Xd, Wd,
									BlankVector, BlankState,
									Xu, Wu,
									stencil_flag);
							}
							break;
						case VISCOUS_RECONSTRUCTION_HYBRID :
							dWdxl = Local_SolnBlk[nb].dWdx[i][j]; 
							dWdyl = Local_SolnBlk[nb].dWdy[i][j]; 

							Xr = Xface;
							Wr = Wface;
							dWdxr = dWdxl; dWdyr = dWdyl; // Hmm.

							Hybrid_Find_dWdX(dWdx, dWdy, 
									Xl, Wl, dWdxl, dWdyl,
									Xr, Wr, dWdxr, dWdyr,
									norm_dir);
							break;
						case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
							dWdx = Local_SolnBlk[nb].dWdx[i][j]; 
							dWdy = Local_SolnBlk[nb].dWdy[i][j];
							break;
					}
				} else if (Mesh[nb].cell_status[i][j  ] == CELL_STATUS_ACTIVE &&
				           Mesh[nb].cell_status[i][j+1] == CELL_STATUS_ACTIVE) {

					if (j == Local_SolnBlk[nb].JCl-1 && 
							(Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX ||
							 Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
							 Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_MOVING_WALL_HEATFLUX ||
							 Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_MOVING_WALL_ISOTHERMAL ||
							 Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_BURNING_SURFACE)) {

						// SOUTH face of cell (i,j+1) is a normal boundary.

						Xr = Local_SolnBlk[nb].Grid.Cell[i][j+1].Xc; 
						Wr = Local_SolnBlk[nb].W[i][j+1];

						if (Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX ||
								Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_MOVING_WALL_HEATFLUX) {
							adiabatic_flag = ON;
						}

						switch (Local_SolnBlk[nb].Grid.BCtypeS[i]) {
							case BC_WALL_VISCOUS_HEATFLUX:
								Wface.rho = Wr.rho; Wface.v.x = ZERO; Wface.v.y = ZERO; 
								Wface.p = Wr.p; Wface.k = Wr.k; Wface.omega =	Wr.omega;
								break;
							case BC_WALL_VISCOUS_ISOTHERMAL:
								Wface.rho = Wr.p/(Wr.R*Local_SolnBlk[nb].Twall); Wface.v.x = ZERO; Wface.v.y = ZERO; 
								Wface.p = Wr.p; Wface.k = Wr.k; Wface.omega = Wr.omega;
								break;
							case BC_MOVING_WALL_HEATFLUX:
								Wface.rho = Wr.rho; Wface.v.x = Local_SolnBlk[nb].Vwall.x; Wface.v.y = Local_SolnBlk[nb].Vwall.y; 
								Wface.p = Wr.p; Wface.k = Wr.k; Wface.omrga = Wr.omega;
								break;
							case BC_MOVING_WALL_ISOTHERMAL:
								Wface.rho = Wr.p/(Wr.R*Local_SolnBlk[nb].Twall); Wface.v.x = Local_SolnBlk[nb].Vwall.x; Wface.v.y = Local_SolnBlk[nb].Vwall.y; 
								Wface.p = Wr.p; Wface.k = Wr.k; Wface.omega = Wr.omega;
								break;
							case BC_BURNING_SURFACE:
								Wface = BurningSurface(Wr, norm_dir);
								break;
						}

						switch(IP->i_Viscous_Reconstruction) {
							case VISCOUS_RECONSTRUCTION_CARTESIAN :
							case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
								{
								int stencil_flag = DIAMONDPATH_RIGHT_TRIANGLE;
								Xu = Local_SolnBlk[nb].Grid.Node[i  ][j+1].X;
								Xd = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X; 
								DiamondPath_Find_dWdX(dWdx, dWdy,
										BlankVector, BlankState, 
										Xd, Wface,
										Xr, Wr,
										Xu, Wface,
										stencil_flag);
								}
								break;
							case VISCOUS_RECONSTRUCTION_HYBRID :
								dWdxr = Local_SolnBlk[nb].dWdx[i][j+1];
								dWdyr = Local_SolnBlk[nb].dWdy[i][j+1];

								Xl = Xface;
								Wl = Wface;
								dWdxl = dWdxr; dWdyl = dWdyr; // Hmm.

								Hybrid_Find_dWdX(dWdx, dWdy, 
										Xl, Wl, dWdxl, dWdyl,
										Xr, Wr, dWdxr, dWdyr,
										norm_dir);
								break;
							case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
								dWdx = Local_SolnBlk[nb].dWdx[i][j+1];
								dWdy = Local_SolnBlk[nb].dWdy[i][j+1];
								break;
						}
					} else if (j == Local_SolnBlk[nb].JCu &&
							(Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX ||
							 Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
							 Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_MOVING_WALL_HEATFLUX ||
							 Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_MOVING_WALL_ISOTHERMAL ||
							 Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_BURNING_SURFACE)) {
						
						// NORTH face of cell (i,j) is a normal boundary.
						
						Xl = Local_SolnBlk[nb].Grid.Cell[i][j].Xc; 
						Wl = Local_SolnBlk[nb].W[i][j];

						if (Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX ||
						    Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_MOVING_WALL_HEATFLUX) {
							adiabatic_flag = ON;
						}

						switch (Local_SolnBlk[nb].Grid.BCtypeN[i]) {
							case BC_WALL_VISCOUS_HEATFLUX:
								Wface.rho = Wl.rho; Wface.v.x = ZERO; Wface.v.y = ZERO; 
								Wface.p = Wl.p; Wface.k = Wl.k; Wface.omega = Wl.omega;
								break;
							case BC_WALL_VISCOUS_ISOTHERMAL:
								Wface.rho = Wl.p/(Wl.R*Local_SolnBlk[nb].Twall); Wface.v.x = ZERO; Wface.v.y = ZERO; 
								Wface.p = Wl.p; Wface.k = Wl.k; Wface.omega = Wl.omega;
								break;
							case BC_MOVING_WALL_HEATFLUX:
								Wface.rho = Wl.rho; Wface.v.x = Local_SolnBlk[nb].Vwall.x; Wface.v.y = Local_SolnBlk[nb].Vwall.y; 
								Wface.p = Wl.p; Wface.k = Wl.k; Wface.omega = Wl.omega;
								break;
							case BC_MOVING_WALL_ISOTHERMAL:
								Wface.rho = Wl.p/(Wl.R*Local_SolnBlk[nb].Twall); Wface.v.x = Local_SolnBlk[nb].Vwall.x; Wface.v.y = Local_SolnBlk[nb].Vwall.y; 
								Wface.p = Wl.p; Wface.k = Wl.k; Wface.omega = Wl.omega;
								break;
							case BC_MOVING_WALL_ISOTHERMAL:
								Wface = BurningSurface(Wl, norm_dir);
								break;
						}

						switch(IP->i_Viscous_Reconstruction) {
							case VISCOUS_RECONSTRUCTION_CARTESIAN :
							case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
								{
								int stencil_flag = DIAMONDPATH_LEFT_TRIANGLE;
								Xu = Local_SolnBlk[nb].Grid.Node[i  ][j+1].X;
								Xd = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X; 
								DiamondPath_Find_dWdX(dWdx, dWdy,
										Xl, Wl,
										Xd, Wface,
										BlankVector, BlankState,
										Xu, Wface,
										stencil_flag);
								}
								break;
							case VISCOUS_RECONSTRUCTION_HYBRID :
								dWdxl = Local_SolnBlk[nb].dWdx[i][j]; 
								dWdyl = Local_SolnBlk[nb].dWdy[i][j]; 

								Xr = Xface; 
								Wr = Wface;
								dWdxr = dWdxl; dWdyr = dWdyl; 

								Hybrid_Find_dWdX(dWdx, dWdy, 
										Xl, Wl, dWdxl, dWdyl,
										Xr, Wr, dWdxr, dWdyr,
										norm_dir);
								break;
							case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
								dWdx = Local_SolnBlk[nb].dWdx[i][j]; 
								dWdy = Local_SolnBlk[nb].dWdy[i][j]; 
								break;
						}
					} else {
						// EAST face is either a normal cell or possibly a non-viscous
						// boundary condition.
						Xl = Local_SolnBlk[nb].Grid.Cell[i][j  ].Xc; Wl = Local_SolnBlk[nb].W[i][j  ];
						Xr = Local_SolnBlk[nb].Grid.Cell[i][j+1].Xc; Wr = Local_SolnBlk[nb].W[i][j+1];

						if (IP->i_Viscous_Reconstruction == VISCOUS_RECONSTRUCTION_CARTESIAN ||
						    IP->i_Viscous_Reconstruction == VISCOUS_RECONSTRUCTION_DIAMOND_PATH) {

							int stencil_flag = DIAMONDPATH_QUADRILATERAL_RECONSTRUCTION;

							Xu = Local_SolnBlk[nb].Grid.Node[i  ][j+1].X; 
							Xd = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X; 

							// These two (WnNW and WnNE) are what make diamond path so expensive.
							Wu = Local_SolnBlk[nb].WnNW(i,j);
							Wd = Local_SolnBlk[nb].WnNE(i,j);

							DiamondPath_Find_dWdX(dWdx, dWdy,
									Xl, Wl, Xd, Wd, Xr, Wr, Xu, Wu, 
									stencil_flag);

							// Find Wface:
							error_flag = Bilinear_Interpolation(Wl, Xl, Wu, Xu, Wr, Xr, Wd, Xd,
									Xface, Wface);
							// and if error_flag?

						} else {
							Wface = HALF*(Wl + Wr);

							dWdxl = Local_SolnBlk[nb].dWdx[i][j]; dWdxr = Local_SolnBlk[nb].dWdy[i][j+1];
							dWdyl = Local_SolnBlk[nb].dWdy[i][j]; dWdyr = Local_SolnBlk[nb].dWdy[i][j+1];

							if (IP->i_Viscous_Reconstruction == VISCOUS_RECONSTRUCTION_HYBRID) {
								Hybrid_Find_dWdX(dWdx, dWdy, 
										Xl, Wl, dWdxl, dWdyl,
										Xr, Wr, dWdxr, dWdyr,
										norm_dir);
							} else if (IP->i_Viscous_Reconstruction == VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE) {
								dWdx = HALF*(dWdxl + dWdxr);
								dWdy = HALF*(dWdyl + dWdyr);
							}
						}
					}
				}
				// Compute the NORTH face viscous flux.
				Flux -= ViscousFlux_n(Xface, Wface, dWdx, dWdy,
						norm_dir, Local_SolnBlk[nb].Axisymmetric, adiabatic_flag);
			} while (0); // end of viscous flux clause.

	  // Evaluate cell-averaged solution changes.
	  Local_SolnBlk[nb].dUdt[i][j  ][0] -= Flux*Local_SolnBlk[nb].Grid.lfaceN(i,j)/Local_SolnBlk[nb].Grid.Cell[i][j].A;
	  Local_SolnBlk[nb].dUdt[i][j+1][0] += Flux*Local_SolnBlk[nb].Grid.lfaceS(i,j+1)/Local_SolnBlk[nb].Grid.Cell[i][j+1].A;

	  // Save south and north face boundary flux.
	  if (j == Local_SolnBlk[nb].JCl-1) {
	    Local_SolnBlk[nb].FluxS[i] = -Flux*Local_SolnBlk[nb].Grid.lfaceS(i,j+1);
	  } else if (j == Local_SolnBlk[nb].JCu) {
	    Local_SolnBlk[nb].FluxN[i] = Flux*Local_SolnBlk[nb].Grid.lfaceN(i,j);
	  }

	}

	Local_SolnBlk[nb].dUdt[i][Local_SolnBlk[nb].JCl-1][0].Vacuum();
	Local_SolnBlk[nb].dUdt[i][Local_SolnBlk[nb].JCu+1][0].Vacuum();

      }

    }
  }

  // Residuals for each quadrilateral multi-block solution block
  // successfully calculated.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::dUdt_Multistage_Explicit --                  *
 *                                                                    *
 * This routine evaluates the stage solution residual for a 1D array  *
 * of 2D quadrilateral multi-block solution blocks.  A variety of     *
 * multistage explicit time integration and upwind finite-volume      *
 * spatial discretization procedures can be used depending on the     *
 * specified input values.                                            *
 *                                                                    *
 **********************************************************************/
template <> int EmbeddedBoundaries2D<NavierStokes2D_cState,
				     NavierStokes2D_pState,
				     NavierStokes2D_Quad_Block,
				     NavierStokes2D_Input_Parameters>::
dUdt_Multistage_Explicit(const int &i_stage, const double &Time) {

  int error_flag, k_residual;
  double omega;
  Vector2D dX;
  NavierStokes2D_pState Wl, Wr;
  NavierStokes2D_cState Flux;

  NavierStokes2D_pState Wu, Wd, dWdxl, dWdyl, dWdxr, dWdyr;
  Vector2D Xl, Xr, Xu, Xd;

  int Interface_BC_Type, Ni;
  Vector2D V;

  double q1, q2, k;

  // Evaluate the solution residual for stage i_stage of n_stage scheme.

  // Evaluate the time step fraction and residual storage location for
  // the stage.
  switch(IP->i_Time_Integration) {
  case TIME_STEPPING_EXPLICIT_EULER :
    omega = Runge_Kutta(i_stage,IP->N_Stage);
    k_residual = 0;
    break;
  case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
    omega = Runge_Kutta(i_stage,IP->N_Stage);
    k_residual = 0;
    break;
  case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA :
    omega = Runge_Kutta(i_stage,IP->N_Stage);
    k_residual = 0;
    if (IP->N_Stage == 4) {
      if (i_stage == 4) k_residual = 0;
      else k_residual = i_stage - 1;
    }
    break;
  case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
    omega = MultiStage_Optimally_Smoothing(i_stage,IP->N_Stage,IP->i_Limiter);
    k_residual = 0;
    break;
  default:
    omega = Runge_Kutta(i_stage,IP->N_Stage);
    k_residual = 0;
    break;
  };

  // Evaluate the time rate of change of the solution (i.e., the
  // solution residuals) using a second-order limited upwind scheme
  // with a variety of flux functions.

  // Evaluate the stage solution residual for each solution block.
  for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
    if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {

      // Perform the linear reconstruction within each cell of the
      // computational grid for this stage.
      Linear_Least_Squares_Reconstruction(nb);

      // Add i-direction (zeta-direction) fluxes.
      for (int j = Local_SolnBlk[nb].JCl-1; j <= Local_SolnBlk[nb].JCu+1; j++) {
	if (i_stage == 1) {
	  Local_SolnBlk[nb].Uo[Local_SolnBlk[nb].ICl-1][j] = Local_SolnBlk[nb].U[Local_SolnBlk[nb].ICl-1][j];
	  Local_SolnBlk[nb].dUdt[Local_SolnBlk[nb].ICl-1][j][k_residual].Vacuum();
	} else {
	  Local_SolnBlk[nb].dUdt[Local_SolnBlk[nb].ICl-1][j][k_residual].Vacuum();
	}

	for (int i = Local_SolnBlk[nb].ICl-1; i <= Local_SolnBlk[nb].ICu; i++) {

	  if (i_stage == 1) {
	    Local_SolnBlk[nb].Uo[i+1][j] = Local_SolnBlk[nb].U[i+1][j];
	    Local_SolnBlk[nb].dUdt[i+1][j][k_residual].Vacuum();
	  } else if (j > Local_SolnBlk[nb].JCl-1 && j < Local_SolnBlk[nb].JCu+1) {
	    switch(IP->i_Time_Integration) {
	    case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
	      //Local_SolnBlk[nb].dUdt[i+1][j][k_residual] = Local_SolnBlk[nb].dUdt[i+1][j][k_residual];
	      break;
	    case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA :
	      if (IP->N_Stage == 2) {
		//Local_SolnBlk[nb].dUdt[i+1][j][k_residual] = Local_SolnBlk[nb].dUdt[i+1][j][k_residual];
	      } else if (IP->N_Stage == 4 && i_stage == 4) {
		Local_SolnBlk[nb].dUdt[i+1][j][k_residual] = Local_SolnBlk[nb].dUdt[i+1][j][0] + 
                                                             TWO*Local_SolnBlk[nb].dUdt[i+1][j][1] +
                                                             TWO*Local_SolnBlk[nb].dUdt[i+1][j][2];
	      } else {
		Local_SolnBlk[nb].dUdt[i+1][j][k_residual].Vacuum();
	      }
	      break;
	    case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
	      Local_SolnBlk[nb].dUdt[i+1][j][k_residual].Vacuum();
	      break;
	    default:
	      Local_SolnBlk[nb].dUdt[i+1][j][k_residual].Vacuum();
	      break;
	    };
	  }

	  if (j >= Local_SolnBlk[nb].JCl && j <= Local_SolnBlk[nb].JCu) {

	    // Evaluate the cell interface i-direction INVISCID fluxes.
	    if (Local_SolnBlk[nb].Grid.lfaceE(i,j) < TOLER) {

	      // EAST face of cell (i,j) has zero length.
	      Wl.Standard_Atmosphere(); Wr.Standard_Atmosphere(); V = Vector2D_ZERO;

	    } else if (Mesh[nb].cell_status[i  ][j] != CELL_STATUS_ACTIVE &&
		       Mesh[nb].cell_status[i+1][j] != CELL_STATUS_ACTIVE) {

	      // EAST face of cell (i,j) is inactive.
	      Wl.Standard_Atmosphere(); Wr.Standard_Atmosphere(); V = Vector2D_ZERO;

	    } else if (Mesh[nb].cell_status[i  ][j] != CELL_STATUS_ACTIVE &&
		       Mesh[nb].cell_status[i+1][j] == CELL_STATUS_ACTIVE) {

	      // WEST face of cell (i+1,j) corresponds to an embedded boundary.
	      dX = Local_SolnBlk[nb].Grid.xfaceW(i+1,j)-Local_SolnBlk[nb].Grid.Cell[i+1][j].Xc;
	      Wr = Local_SolnBlk[nb].W[i+1][j] + (Local_SolnBlk[nb].phi[i+1][j]^Local_SolnBlk[nb].dWdx[i+1][j])*dX.x +
		                                 (Local_SolnBlk[nb].phi[i+1][j]^Local_SolnBlk[nb].dWdy[i+1][j])*dX.y;
	      Ni = Mesh[nb].cell_status[i][j];
	      // Determine the boundary condition type.
	      Interface_BC_Type = Interface_Union_List[Ni].Determine_Interface_BC_Type(Local_SolnBlk[nb].Grid.xfaceW(i+1,j));
	      // Determine the boundary velocity.
	      V = Interface_Union_List[Ni].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceW(i+1,j),Time);
	      // Determine the left state by applying the appropriate
	      // boundary condition.
	      if (Interface_BC_Type == INTERFACE_BC_REFLECTION) {
		// WEST face of cell (i+1,j) is a REFLECTION boundary.
		Wl = Reflect(Wr,Local_SolnBlk[nb].Grid.nfaceW(i+1,j),V);
	      } else if (Interface_BC_Type == INTERFACE_BC_BURNING_SURFACE) {
		// WEST face of cell (i+1,j) is a BURNING_SURFACE boundary.
		Wl = BurningSurface(Wr,Local_SolnBlk[nb].Grid.nfaceW(i+1,j));
	      } else if (Interface_BC_Type == INTERFACE_BC_WALL_VISCOUS_HEATFLUX) {
		// WEST face of cell (i+1,j) is a WALL_VISCOUS_HEATFLUX boundary.
		Wl = WallViscousHeatFlux(Wr,Local_SolnBlk[nb].Grid.nfaceW(i+1,j),V);
	      } else if (Interface_BC_Type == INTERFACE_BC_WALL_VISCOUS_ISOTHERMAL) {
		// WEST face of cell (i+1,j) is a WALL_VISCOUS_ISOTHERMAL boundary.
		//Temperature = Interface_Union_List[Ni].Determine_Interface_Temperature(Local_SolnBlk[nb].Grid.xfaceW(i+1,j));
		Wl = WallViscousIsothermal(Wr,Local_SolnBlk[nb].Grid.nfaceW(i+1,j),V,IP->Twall);
	      } else if (Interface_BC_Type == INTERFACE_BC_RINGLEB) {
		// WEST face of cell (i+1,j) is a RINGLEB_FLOW boundary.
		Wl = RinglebFlow(Wr,Local_SolnBlk[nb].Grid.Cell[i][j].Xc,q1,k);
		Wl = RinglebFlow(Wr,Local_SolnBlk[nb].Grid.Cell[i+1][j].Xc,q2,k);
		Wl = RinglebFlow(Wr,Local_SolnBlk[nb].Grid.xfaceW(i+1,j));
		if (q1 < IP->Isotach_Line && q2 > IP->Isotach_Line)
		  Wl = BC_Characteristic_Pressure(Wr,Wl,Local_SolnBlk[nb].Grid.nfaceW(i+1,j));
	      }

	    } else if (Mesh[nb].cell_status[i  ][j] == CELL_STATUS_ACTIVE &&
		       Mesh[nb].cell_status[i+1][j] != CELL_STATUS_ACTIVE) {

	      // EAST face of cell (i,j) corresponds to an embedded boundary.
	      dX = Local_SolnBlk[nb].Grid.xfaceE(i,j) - Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	      Wl = Local_SolnBlk[nb].W[i][j] + (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdx[i][j])*dX.x +
	                                       (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdy[i][j])*dX.y;
	      Ni = Mesh[nb].cell_status[i+1][j];
	      // Determine the boundary condition type.
	      Interface_BC_Type = Interface_Union_List[Ni].Determine_Interface_BC_Type(Local_SolnBlk[nb].Grid.xfaceE(i,j));
	      // Determine the boundary velocity.
	      V = Interface_Union_List[Ni].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceE(i,j),Time);
	      // Determine the right state by applying the appropriate
	      // boundary condition.
	      if (Interface_BC_Type == INTERFACE_BC_REFLECTION) {
		// EAST face of cell (i,j) is a REFLECTION boundary.
		Wr = Reflect(Wl,Local_SolnBlk[nb].Grid.nfaceE(i,j),V);
	      } else if (Interface_BC_Type == INTERFACE_BC_BURNING_SURFACE) {
		// EAST face of cell (i,j) is a BURNING_SURFACE boundary.
		Wr = BurningSurface(Wl,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      } else if (Interface_BC_Type == INTERFACE_BC_WALL_VISCOUS_HEATFLUX) {
		// EAST face of cell (i,j) is a WALL_VISCOUS_HEATFLUX boundary.
		Wr = WallViscousHeatFlux(Wl,Local_SolnBlk[nb].Grid.nfaceE(i,j),V);
	      } else if (Interface_BC_Type == INTERFACE_BC_WALL_VISCOUS_ISOTHERMAL) {
		// EAST face of cell (i,j) is a WALL_VISCOUS_ISOTHERMAL boundary.
		//Temperature = Interface_Union_List[Ni].Determine_Interface_Temperature(Local_SolnBlk[nb].Grid.xfaceE(i,j));
		Wr = WallViscousIsothermal(Wl,Local_SolnBlk[nb].Grid.nfaceE(i,j),V,IP->Twall);
	      } else if (Interface_BC_Type == INTERFACE_BC_RINGLEB) {
		// EAST face of cell (i,j) is a RINGLEB FLOW boundary.
		Wr = RinglebFlow(Wr,Local_SolnBlk[nb].Grid.Cell[i+1][j].Xc,q1,k);
		Wr = RinglebFlow(Wr,Local_SolnBlk[nb].Grid.Cell[i][j].Xc,q2,k);
		Wr = RinglebFlow(Wr,Local_SolnBlk[nb].Grid.xfaceE(i,j));
		if (q1 < IP->Isotach_Line && q2 > IP->Isotach_Line)
		  Wr = BC_Characteristic_Pressure(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      }

	    } else if (Mesh[nb].cell_status[i  ][j] == CELL_STATUS_ACTIVE &&
		       Mesh[nb].cell_status[i+1][j] == CELL_STATUS_ACTIVE) {

	      V = Vector2D_ZERO;

	      if (i == Local_SolnBlk[nb].ICl-1 && 
		  (Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_REFLECTION ||
		   Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_BURNING_SURFACE ||
		   Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX ||
		   Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
		   Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_MOVING_WALL_HEATFLUX ||
		   Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL ||
		   Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_RINGLEB_FLOW ||
		   Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_CHARACTERISTIC)) {

		dX = Local_SolnBlk[nb].Grid.xfaceW(i+1,j) - Local_SolnBlk[nb].Grid.Cell[i+1][j].Xc;
		Wr = Local_SolnBlk[nb].W[i+1][j] + (Local_SolnBlk[nb].phi[i+1][j]^Local_SolnBlk[nb].dWdx[i+1][j])*dX.x +
	                                           (Local_SolnBlk[nb].phi[i+1][j]^Local_SolnBlk[nb].dWdy[i+1][j])*dX.y;

		// WEST face of cell (i+1,j) is a normal boundary.
		if (Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_REFLECTION) {
		  // WEST face of cell (i+1,j) is a REFLECTION boundary.
		  Wl = Reflect(Wr,Local_SolnBlk[nb].Grid.nfaceW(i+1,j));
		} else if (Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_BURNING_SURFACE) {
		  // WEST face of cell (i+1,j) is a BURNING_SURFACE boundary.
		  Wl = BurningSurface(Wr,Local_SolnBlk[nb].Grid.nfaceW(i+1,j));
		} else if (Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX) {
		  // WEST face of cell (i+1,j) is a WALL_VISCOUS_HEATFLUX boundary.
		  Wl = WallViscousHeatFlux(Wr,Local_SolnBlk[nb].Grid.nfaceW(i+1,j));
		} else if (Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL) {
		  // WEST face of cell (i+1,j) is a WALL_VISCOUS_ISOTHERMAL boundary.
		  Wl = WallViscousIsothermal(Wr,Local_SolnBlk[nb].Grid.nfaceW(i+1,j),Local_SolnBlk[nb].Twall);
		} else if (Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_MOVING_WALL) {
		  // WEST face of cell (i+1,j) is a MOVINGWALL_HEATFLUX boundary.
		  Wl = MovingWallHeatFlux(Wr,Local_SolnBlk[nb].Grid.nfaceW(i+1,j),Local_SolnBlk[nb].Vwall.x);
		} else if (Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL) {
		  // WEST face of cell (i+1,j) is a MOVINGWALL_ISOTHERMAL boundary.
		  Wl = MovingWallIsothermal(Wr,Local_SolnBlk[nb].Grid.nfaceW(i+1,j),Local_SolnBlk[nb].Vwall.x,Local_SolnBlk[nb].Twall);
		} else if (Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_RINGLEB_FLOW) {
		  // WEST face of cell (i+1,j) is a RINGLEB_FLOW boundary.
		  Wl = RinglebFlow(Wl,Local_SolnBlk[nb].Grid.xfaceW(i+1,j));
		} else {
		  // WEST face of cell (i+1,j) is a CHARACTERISTIC boundary.
		  Wl = BC_Characteristic_Pressure(Wr,Local_SolnBlk[nb].WoW[j],Local_SolnBlk[nb].Grid.nfaceW(i+1,j));
		}

	      } else if (i == Local_SolnBlk[nb].ICu &&
			 (Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_REFLECTION ||
			  Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_BURNING_SURFACE ||
			  Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX ||
			  Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
			  Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_MOVING_WALL_HEATFLUX ||
			  Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL ||
			  Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_RINGLEB_FLOW ||
			  Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_CHARACTERISTIC)) {

		dX = Local_SolnBlk[nb].Grid.xfaceE(i,j) - Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
		Wl = Local_SolnBlk[nb].W[i][j] + (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdx[i][j])*dX.x +
	                                         (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdy[i][j])*dX.y;

		// EAST face of cell (i,j) is a normal boundary.
		if (Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_REFLECTION) {
		  // EAST face of cell (i,j) is a REFLECTION boundary.
		  Wr = Reflect(Wl,Local_SolnBlk[nb].Grid.nfaceE(i,j));
		} else if (Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_BURNING_SURFACE) {
		  // EAST face of cell (i,j) is a BURNING_SURFACE boundary.
		  Wr = BurningSurface(Wl,Local_SolnBlk[nb].Grid.nfaceE(i,j));
		} else if (Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX) {
		  // EAST face of cell (i,j) is a WALL_VISCOUS_HEATFLUX boundary.
		  Wr = WallViscousHeatFlux(Wl,Local_SolnBlk[nb].Grid.nfaceE(i,j));
		} else if (Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL) {
		  // EAST face of cell (i,j) is a WALL_VISCOUS_ISOTHERMAL boundary.
		  Wr = WallViscousIsothermal(Wl,Local_SolnBlk[nb].Grid.nfaceE(i,j),Local_SolnBlk[nb].Twall);
		} else if (Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_MOVING_WALL_HEATFLUX) {
		  // EAST face of cell (i,j) is a MOVINGWALL_HEATFLUX boundary.
		  Wr = MovingWallHeatFlux(Wl,Local_SolnBlk[nb].Grid.nfaceE(i,j),Local_SolnBlk[nb].Vwall.x);
		} else if (Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL) {
		  // EAST face of cell (i,j) is a MOVINGWALL_ISOTHERMAL boundary.
		  Wr = MovingWallIsothermal(Wl,Local_SolnBlk[nb].Grid.nfaceE(i,j),Local_SolnBlk[nb].Vwall.x,Local_SolnBlk[nb].Twall);
		} else if (Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_RINGLEB_FLOW) {
		  // EAST face of cell (i,j) is a RINGLEB_FLOW boundary.
		  Wr = RinglebFlow(Wr,Local_SolnBlk[nb].Grid.xfaceE(i,j));
		} else {
		  // EAST face of cell (i,j) is a CHARACTERISTIC boundary.
		  Wr = BC_Characteristic_Pressure(Wl,Local_SolnBlk[nb].WoE[j],Local_SolnBlk[nb].Grid.nfaceE(i,j));
		}

	      } else {

		// EAST face is either a normal cell or possibly a FIXED, 
		// NONE or EXTRAPOLATION boundary.
		dX = Local_SolnBlk[nb].Grid.xfaceE(i  ,j) - Local_SolnBlk[nb].Grid.Cell[i  ][j].Xc;
		Wl = Local_SolnBlk[nb].W[i  ][j] + (Local_SolnBlk[nb].phi[i  ][j]^Local_SolnBlk[nb].dWdx[i  ][j])*dX.x +
	                                           (Local_SolnBlk[nb].phi[i  ][j]^Local_SolnBlk[nb].dWdy[i  ][j])*dX.y;
		dX = Local_SolnBlk[nb].Grid.xfaceW(i+1,j) - Local_SolnBlk[nb].Grid.Cell[i+1][j].Xc;
		Wr = Local_SolnBlk[nb].W[i+1][j] + (Local_SolnBlk[nb].phi[i+1][j]^Local_SolnBlk[nb].dWdx[i+1][j])*dX.x +
	                                           (Local_SolnBlk[nb].phi[i+1][j]^Local_SolnBlk[nb].dWdy[i+1][j])*dX.y;

	      }

	    }

// 	    dout << endl << " EAST";
// 	    dout << endl << Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
// 	    dout << endl << Wl;
// 	    dout << endl << Wr;
// 	    dout << endl << FluxRoe_MB_n(Wl,Wr,V,Local_SolnBlk[nb].Grid.nfaceE(i,j));

	    // Determine EAST face INVISCID flux.
	    switch(IP->i_Flux_Function) {
	    case FLUX_FUNCTION_GODUNOV :
	      Flux = FluxGodunov_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_ROE :
	      Flux = FluxRoe_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_RUSANOV :
	      Flux = FluxRusanov_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_HLLE :
	      Flux = FluxHLLE_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_HLLL :
	      Flux = FluxHLLL_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_HLLC :
	      Flux = FluxHLLC_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_VANLEER :
	      Flux = FluxVanLeer_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_AUSM :
	      Flux = FluxAUSM_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_AUSMplus :
	      Flux = FluxAUSMplus_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_GODUNOV_MB :
	      Flux = FluxGodunov_MB_n(Wl,Wr,V,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_ROE_MB :
	      Flux = FluxRoe_MB_n(Wl,Wr,V,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_HLLE_MB :
	      Flux = FluxHLLE_MB_n(Wl,Wr,V,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_VANLEER_MB :
	      Flux = FluxVanLeer_MB_n(Wl,Wr,V,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    default:
	      Flux = FluxRoe_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    };

	    // Compute the cell centred stress tensor and heat flux vector if required.
 	    if (Local_SolnBlk[nb].Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
 		(Local_SolnBlk[nb].Flow_Type && Local_SolnBlk[nb].Axisymmetric)) {
 	      Local_SolnBlk[nb].W[i][j].ComputeViscousTerms(Local_SolnBlk[nb].dWdx[i][j],
 							    Local_SolnBlk[nb].dWdy[i][j],
 							    Local_SolnBlk[nb].Grid.Cell[i][j].Xc,
 							    Local_SolnBlk[nb].Axisymmetric);
 	      Local_SolnBlk[nb].U[i][j].tau = Local_SolnBlk[nb].W[i][j].tau;
 	      Local_SolnBlk[nb].U[i][j].q = Local_SolnBlk[nb].W[i][j].q;
 	    }

	    // Determine the i-direction viscous flux.

			do { // with while (0) to allow jumps to the end of the viscous flux calculation.
	    	if (Local_SolnBlk[nb].Flow_Type == 0) { break; }

				if (Local_SolnBlk[nb].Grid.lfaceE(i,j) < TOLER) { break; }

				if (Mesh[nb].cell_status[i  ][j] != CELL_STATUS_ACTIVE &&
						Mesh[nb].cell_status[i+1][j] != CELL_STATUS_ACTIVE) {
					break;
				}

				int adiabatic_flag = OFF;
				Vector2D Xface = Local_SolnBlk[nb].Grid.xfaceE(i, j);
				Vector2D norm_dir = Local_SolnBlk[nb].Grid.nfaceE(i, j);

	      if (Mesh[nb].cell_status[i  ][j] != CELL_STATUS_ACTIVE &&
				    Mesh[nb].cell_status[i+1][j] == CELL_STATUS_ACTIVE) {

					// WEST face of cell (i+1,j) corresponds to an embedded boundary.
					Xr = Local_SolnBlk[nb].Grid.Cell[i+1][j].Xc; 
					Wr = Local_SolnBlk[nb].W[i+1][j];

					Ni = Mesh[nb].cell_status[i][j];
					Interface_BC_Type = Interface_Union_List[Ni].Determine_Interface_BC_Type(Xface);
					
					// No viscous flux for a reflected boundary. Is this true?
					if (Interface_BC_Type == INTERFACE_BC_REFLECTION) { break; } 

					switch (Interface_BC_Type) {
						// The face velocity is the same regardless of the BC type 
						// and is set below.
						case INTERFACE_BC_WALL_VISCOUS_HEATFLUX:
							Wface.rho = Wr.rho; 
							Wface.p = Wr.p; Wface.k = Wr.k; Wface.omega = Wr.omega;
							break;
						case INTERFACE_BC_WALL_VISCOUS_ISOTHERMAL:
							Wface.rho = Wr.p/(Wr.R*Local_SolnBlk[nb].Twall); 
							Wface.p = Wr.p; Wface.k = Wr.k; Wface.omega = Wr.omega;
							break;
						case INTERFACE_BC_BURNING_SURFACE:
							Wface = BurningSurface(Wr, norm_dir);
							break;
					}

					Wface.v = Interface_Union_List[Ni].Determine_Interface_Velocity(Xface, Time);

					if (Interface_BC_Type == INTERFACE_BC_WALL_VISCOUS_HEATFLUX) {
						adiabatic_flag = ON;
					}

					switch(IP->i_Viscous_Reconstruction) {
						case VISCOUS_RECONSTRUCTION_CARTESIAN :
						case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
							{
							int stencil_flag = DIAMONDPATH_RIGHT_TRIANGLE;
							Xd = Local_SolnBlk[nb].Grid.Node[i+1][j  ].X;
							Wd = Wface;
							Wd.v = Interface_Union_List[Ni].Determine_Interface_Velocity(Xd, Time);

							Xu = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
							Wu = Wface;
							Wu.v = Interface_Union_List[Ni].Determine_Interface_Velocity(Xu, Time);

							DiamondPath_Find_dWdX(dWdx, dWdy,
									BlankVector, BlankState, 
									Xd, Wd,
									Xr, Wr,
									Xu, Wu,
									stencil_flag);
							}
							break;
						case VISCOUS_RECONSTRUCTION_HYBRID :
							dWdxr = Local_SolnBlk[nb].dWdx[i+1][j]; 
							dWdyr = Local_SolnBlk[nb].dWdy[i+1][j];

							Xl = Xface;
							Wl = Wface;
							dWdxl = dWdxr; dWdyl = dWdyr; // Hmm. 

							Hybrid_Find_dWdX(dWdx, dWdy, 
									Xl, Wl, dWdxl, dWdyl,
									Xr, Wr, dWdxr, dWdyr,
									norm_dir);
							break;
						case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
							dWdx = Local_SolnBlk[nb].dWdx[i+1][j]; 
							dWdy = Local_SolnBlk[nb].dWdy[i+1][j];
							break;
					}
// #ifdef _EB_PARALLEL_DEBUG_
// 		dout << endl << " ==========================================================";
// 		dout << endl << " WEST face of cell (" << i << "," << j << "," << Local_Solution_Block_List->Block[nb].gblknum << ") VISCOUS flux:";
// 		dout << endl << " Embbedded boundary with bc = " << Interface_BC_Type;
// 		dout << endl << " Diamond-path = " << viscous_bc_flag;
// 		dout << endl << " xfaceE =" << Local_SolnBlk[nb].Grid.xfaceE(i,j);
// 		dout << endl << " nfaceE =" << Local_SolnBlk[nb].Grid.nfaceE(i,j);
// 		dout << endl << " lfaceE = " << Local_SolnBlk[nb].Grid.lfaceE(i,j);
// 		dout << endl << " Xr =" << Xr;
// 		dout << endl << " Xu =" << Xu;
// 		dout << endl << " Xd =" << Xd;
// 		dout << endl << " Wr =" << Wr;
// 		dout << endl << " Wu =" << Wu;
// 		dout << endl << " Wd =" << Wd;
// 		dout << endl << " F =" << ViscousFluxDiamondPath_n(Local_SolnBlk[nb].Grid.xfaceE(i,j),
// 								   Xl,Wl,Xd,Wd,Xr,Wr,Xu,Wu,
// 								   Local_SolnBlk[nb].Grid.nfaceE(i,j),
// 								   Local_SolnBlk[nb].Axisymmetric,
// 								   viscous_bc_flag);
// #endif
				} else if (Mesh[nb].cell_status[i  ][j] == CELL_STATUS_ACTIVE &&
				           Mesh[nb].cell_status[i+1][j] != CELL_STATUS_ACTIVE) {

					// EAST face of cell (i,j) corresponds to an embedded boundary.
					Xl = Local_SolnBlk[nb].Grid.Cell[i][j].Xc; 
					Wl = Local_SolnBlk[nb].W[i][j];

					Ni = Mesh[nb].cell_status[i+1][j];
					Interface_BC_Type = Interface_Union_List[Ni].Determine_Interface_BC_Type(Xface);

					// No viscous flux for a reflected boundary. Is this true?
					if (Interface_BC_Type == INTERFACE_BC_REFLECTION) { break; } 

					switch (Interface_BC_Type) {
						// The face velocity is the same regardless of the BC type 
						// and is set below.
						case INTERFACE_BC_WALL_VISCOUS_HEATFLUX:
							Wface.rho = Wl.rho; 
							Wface.p = Wl.p; Wface.k = Wl.k; Wface.omega = Wl.omega;
							break;
						case INTERFACE_BC_WALL_VISCOUS_ISOTHERMAL:
							Wface.rho = Wl.p/(Wl.R*IP->Twall); 
							Wface.p = Wl.p; Wface.k = Wl.k; Wface.omega = Wl.omega;
							break;
						case INTERFACE_BC_BURNING_SURFACE:
							Wface = BurningSurface(Wl, norm_dir);
							break;
					}

					Wface.v = Interface_Union_List[Ni].Determine_Interface_Velocity(Xface, Time);

					if (Interface_BC_Type == INTERFACE_BC_WALL_VISCOUS_HEATFLUX) {
						adiabatic_flag = ON;
					}

					switch(IP->i_Viscous_Reconstruction) {
						case VISCOUS_RECONSTRUCTION_CARTESIAN :
						case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
							{
							int stencil_flag = DIAMONDPATH_LEFT_TRIANGLE;
							Xd = Local_SolnBlk[nb].Grid.Node[i+1][j  ].X;
							Wd = Wface;
							Wd.v = Interface_Union_List[Ni].Determine_Interface_Velocity(Xd, Time);

							Xu = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
							Wu = Wface;
							Wu.v = Interface_Union_List[Ni].Determine_Interface_Velocity(Xu, Time);

							DiamondPath_Find_dWdX(dWdx, dWdy,
									Xl, Wl,
									Xd, Wd,
									BlankVector, BlankState,
									Xu, Wu,
									stencil_flag);
							}
							break;
						case VISCOUS_RECONSTRUCTION_HYBRID :
							dWdxl = Local_SolnBlk[nb].dWdx[i][j]; 
							dWdyl = Local_SolnBlk[nb].dWdy[i][j]; 

							Xr = Xface;
							Wr = Wface;
							dWdxr = dWdxl; dWdyr = dWdyl; // Hmm.

							Hybrid_Find_dWdX(dWdx, dWdy, 
									Xl, Wl, dWdxl, dWdyl,
									Xr, Wr, dWdxr, dWdyr,
									norm_dir);
							break;
						case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
							dWdx = Local_SolnBlk[nb].dWdx[i][j]; 
							dWdy = Local_SolnBlk[nb].dWdy[i][j];
							break;
					}
// #ifdef _EB_PARALLEL_DEBUG_
// 		dout << endl << " ==========================================================";
// 		dout << endl << " EAST face of cell (" << i << "," << j << "," << Local_Solution_Block_List->Block[nb].gblknum << ") VISCOUS flux:";
// 		dout << endl << " Embbedded boundary with bc = " << Interface_BC_Type;
// 		dout << endl << " Diamond-path = " << viscous_bc_flag;
// 		dout << endl << " xfaceE =" << Local_SolnBlk[nb].Grid.xfaceE(i,j);
// 		dout << endl << " nfaceE =" << Local_SolnBlk[nb].Grid.nfaceE(i,j);
// 		dout << endl << " lfaceE = " << Local_SolnBlk[nb].Grid.lfaceE(i,j);
// 		dout << endl << " Xl =" << Xl;
// 		dout << endl << " Xu =" << Xu;
// 		dout << endl << " Xd =" << Xd;
// 		dout << endl << " Wl =" << Wl;
// 		dout << endl << " Wu =" << Wu;
// 		dout << endl << " Wd =" << Wd;
// 		dout << endl << " F =" << ViscousFluxDiamondPath_n(Local_SolnBlk[nb].Grid.xfaceE(i,j),
// 								   Xl,Wl,Xd,Wd,Xr,Wr,Xu,Wu,
// 								   Local_SolnBlk[nb].Grid.nfaceE(i,j),
// 								   Local_SolnBlk[nb].Axisymmetric,
// 								   viscous_bc_flag);
// #endif

				} else if (Mesh[nb].cell_status[i  ][j] == CELL_STATUS_ACTIVE &&
				           Mesh[nb].cell_status[i+1][j] == CELL_STATUS_ACTIVE) {

					if (i == Local_SolnBlk[nb].ICl-1 && 
							(Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX ||
							 Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
							 Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_MOVING_WALL_HEATFLUX ||
							 Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL ||
							 Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_BURNING_SURFACE)) {

						// WEST face of cell (i+1,j) is a normal boundary.

						Xr = Local_SolnBlk[nb].Grid.Cell[i+1][j].Xc; 
						Wr = Local_SolnBlk[nb].W[i+1][j];

						if (Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX ||
								Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_MOVING_WALL_HEATFLUX) {
							adiabatic_flag = ON;
						}

						switch (Local_SolnBlk[nb].Grid.BCtypeW[j]) {
							case BC_WALL_VISCOUS_HEATFLUX:
								Wface.rho = Wr.rho; Wface.v.x = ZERO; Wface.v.y = ZERO; 
								Wface.p = Wr.p; Wface.k = Wr.k; Wface.omega =	Wr.omega;
								break;
							case BC_WALL_VISCOUS_ISOTHERMAL:
								Wface.rho = Wr.p/(Wr.R*Local_SolnBlk[nb].Twall); Wface.v.x = ZERO; Wface.v.y = ZERO; 
								Wface.p = Wr.p; Wface.k = Wr.k; Wface.omega = Wr.omega;
								break;
							case BC_MOVING_WALL_HEATFLUX:
								Wface.rho = Wr.rho; Wface.v.x = Local_SolnBlk[nb].Vwall.x; Wface.v.y = Local_SolnBlk[nb].Vwall.y; 
								Wface.p = Wr.p; Wface.k = Wr.k; Wface.omrga = Wr.omega;
								break;
							case BC_MOVING_WALL_ISOTHERMAL:
								Wface.rho = Wr.p/(Wr.R*Local_SolnBlk[nb].Twall); Wface.v.x = Local_SolnBlk[nb].Vwall.x; Wface.v.y = Local_SolnBlk[nb].Vwall.y; 
								Wface.p = Wr.p; Wface.k = Wr.k; Wface.omega = Wr.omega;
								break;
							case BC_BURNING_SURFACE:
								Wface = BurningSurface(Wr, norm_dir);
								break;
						}

						switch(IP->i_Viscous_Reconstruction) {
							case VISCOUS_RECONSTRUCTION_CARTESIAN :
							case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
								{
								int stencil_flag = DIAMONDPATH_RIGHT_TRIANGLE;
								Xu = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
								Xd = Local_SolnBlk[nb].Grid.Node[i+1][j  ].X; 
								DiamondPath_Find_dWdX(dWdx, dWdy,
										BlankVector, BlankState, 
										Xd, Wface,
										Xr, Wr,
										Xu, Wface,
										stencil_flag);
								}
								break;
							case VISCOUS_RECONSTRUCTION_HYBRID :
								dWdxr = Local_SolnBlk[nb].dWdx[i+1][j];
								dWdyr = Local_SolnBlk[nb].dWdy[i+1][j];

								Xl = Xface;
								Wl = Wface;
								dWdxl = dWdxr; dWdyl = dWdyr; // Hmm.

								Hybrid_Find_dWdX(dWdx, dWdy, 
										Xl, Wl, dWdxl, dWdyl,
										Xr, Wr, dWdxr, dWdyr,
										norm_dir);
								break;
							case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
								dWdx = Local_SolnBlk[nb].dWdx[i+1][j];
								dWdy = Local_SolnBlk[nb].dWdy[i+1][j];
								break;
						}
// #ifdef _EB_PARALLEL_DEBUG_
// 		dout << endl << " ==========================================================";
// 		dout << endl << " WEST face of cell (" << i+1 << "," << j << "," << Local_Solution_Block_List->Block[nb].gblknum << ") VISCOUS flux:";
// 		dout << endl << " Normal boundary with bc = " << Local_SolnBlk[nb].Grid.BCtypeW[j];
// 		dout << endl << " Diamond-path = " << viscous_bc_flag;
// 		dout << endl << " xfaceW =" << Local_SolnBlk[nb].Grid.xfaceW(i+1,j);
// 		dout << endl << " nfaceW =" << Local_SolnBlk[nb].Grid.nfaceW(i+1,j);
// 		dout << endl << " lfaceW = " << Local_SolnBlk[nb].Grid.lfaceW(i+1,j);
// 		dout << endl << " Xr =" << Xr;
// 		dout << endl << " Xu =" << Xu;
// 		dout << endl << " Xd =" << Xd;
// 		dout << endl << " Wr =" << Wr;
// 		dout << endl << " Wu =" << Wu;
// 		dout << endl << " Wd =" << Wd;
// 		dout << endl << " F =" << ViscousFluxDiamondPath_n(Local_SolnBlk[nb].Grid.xfaceE(i,j),
// 								   Xl,Wl,Xd,Wd,Xr,Wr,Xu,Wu,
// 								   Local_SolnBlk[nb].Grid.nfaceE(i,j),
// 								   Local_SolnBlk[nb].Axisymmetric,
// 								   viscous_bc_flag);
// #endif

					} else if (i == Local_SolnBlk[nb].ICu &&
							(Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX ||
							 Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
							 Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_MOVING_WALL_HEATFLUX ||
							 Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL ||
							 Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_BURNING_SURFACE)) {
						
						// EAST face of cell (i,j) is a normal boundary.
						
						Xl = Local_SolnBlk[nb].Grid.Cell[i][j].Xc; 
						Wl = Local_SolnBlk[nb].W[i][j];

						if (Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX ||
						    Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_MOVING_WALL_HEATFLUX) {
							adiabatic_flag = ON;
						}

						switch (Local_SolnBlk[nb].Grid.BCtypeE[j]) {
							case BC_WALL_VISCOUS_HEATFLUX:
								Wface.rho = Wl.rho; Wface.v.x = ZERO; Wface.v.y = ZERO; 
								Wface.p = Wl.p; Wface.k = Wl.k; Wface.omega = Wl.omega;
								break;
							case BC_WALL_VISCOUS_ISOTHERMAL:
								Wface.rho = Wl.p/(Wl.R*Local_SolnBlk[nb].Twall); Wface.v.x = ZERO; Wface.v.y = ZERO; 
								Wface.p = Wl.p; Wface.k = Wl.k; Wface.omega = Wl.omega;
								break;
							case BC_MOVING_WALL_HEATFLUX:
								Wface.rho = Wl.rho; Wface.v.x = Local_SolnBlk[nb].Vwall.x; Wface.v.y = Local_SolnBlk[nb].Vwall.y; 
								Wface.p = Wl.p; Wface.k = Wl.k; Wface.omega = Wl.omega;
								break;
							case BC_MOVING_WALL_ISOTHERMAL:
								Wface.rho = Wl.p/(Wl.R*Local_SolnBlk[nb].Twall); Wface.v.x = Local_SolnBlk[nb].Vwall.x; Wface.v.y = Local_SolnBlk[nb].Vwall.y; 
								Wface.p = Wl.p; Wface.k = Wl.k; Wface.omega = Wl.omega;
								break;
							case BC_MOVING_WALL_ISOTHERMAL:
								Wface = BurningSurface(Wl, norm_dir);
								break;
						}

						switch(IP->i_Viscous_Reconstruction) {
							case VISCOUS_RECONSTRUCTION_CARTESIAN :
							case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
								{
								int stencil_flag = DIAMONDPATH_LEFT_TRIANGLE;
								Xu = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
								Xd = Local_SolnBlk[nb].Grid.Node[i+1][j  ].X; 
								DiamondPath_Find_dWdX(dWdx, dWdy,
										Xl, Wl,
										Xd, Wface,
										BlankVector, BlankState,
										Xu, Wface,
										stencil_flag);
								}
								break;
							case VISCOUS_RECONSTRUCTION_HYBRID :
								dWdxl = Local_SolnBlk[nb].dWdx[i][j]; 
								dWdyl = Local_SolnBlk[nb].dWdy[i][j]; 

								Xr = Xface; 
								Wr = Wface;
								dWdxr = dWdxl; dWdyr = dWdyl; 

								Hybrid_Find_dWdX(dWdx, dWdy, 
										Xl, Wl, dWdxl, dWdyl,
										Xr, Wr, dWdxr, dWdyr,
										norm_dir);
								break;
							case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
								dWdx = Local_SolnBlk[nb].dWdx[i][j]; 
								dWdy = Local_SolnBlk[nb].dWdy[i][j]; 
								break;
						}
// #ifdef _EB_PARALLEL_DEBUG_
// 		dout << endl << " ==========================================================";
// 		dout << endl << " EAST face of cell (" << i << "," << j << "," << Local_Solution_Block_List->Block[nb].gblknum << ") VISCOUS flux:";
// 		dout << endl << " Normal boundary with bc = " << Local_SolnBlk[nb].Grid.BCtypeE[j];
// 		dout << endl << " Diamond-path = " << viscous_bc_flag;
// 		dout << endl << " xfaceE =" << Local_SolnBlk[nb].Grid.xfaceE(i,j);
// 		dout << endl << " nfaceE =" << Local_SolnBlk[nb].Grid.nfaceE(i,j);
// 		dout << endl << " lfaceE = " << Local_SolnBlk[nb].Grid.lfaceE(i,j);
// 		dout << endl << " Xl =" << Xl;
// 		dout << endl << " Xu =" << Xu;
// 		dout << endl << " Xd =" << Xd;
// 		dout << endl << " Wl =" << Wl;
// 		dout << endl << " Wu =" << Wu;
// 		dout << endl << " Wd =" << Wd;
// 		dout << endl << " F =" << ViscousFluxDiamondPath_n(Local_SolnBlk[nb].Grid.xfaceE(i,j),
// 								   Xl,Wl,Xd,Wd,Xr,Wr,Xu,Wu,
// 								   Local_SolnBlk[nb].Grid.nfaceE(i,j),
// 								   Local_SolnBlk[nb].Axisymmetric,
// 								   viscous_bc_flag);
// #endif

					} else {
						// EAST face is either a normal cell or possibly a non-viscous
						// boundary condition.
						Xl = Local_SolnBlk[nb].Grid.Cell[i  ][j].Xc; Wl = Local_SolnBlk[nb].W[i  ][j];
						Xr = Local_SolnBlk[nb].Grid.Cell[i+1][j].Xc; Wr = Local_SolnBlk[nb].W[i+1][j];

						if (IP->i_Viscous_Reconstruction == VISCOUS_RECONSTRUCTION_CARTESIAN ||
						    IP->i_Viscous_Reconstruction == VISCOUS_RECONSTRUCTION_DIAMOND_PATH) {

							int stencil_flag = DIAMONDPATH_QUADRILATERAL_RECONSTRUCTION;

							Xu = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X; 
							Xd = Local_SolnBlk[nb].Grid.Node[i+1][j  ].X; 

							// These two (WnNE and WnSE) are what make diamond path so expensive.
							Wu = Local_SolnBlk[nb].WnNE(i,j);
							Wd = Local_SolnBlk[nb].WnSE(i,j);

							DiamondPath_Find_dWdX(dWdx, dWdy,
									Xl, Wl, Xd, Wd, Xr, Wr, Xu, Wu, 
									stencil_flag);

							// Find Wface:
							error_flag = Bilinear_Interpolation(Wl, Xl, Wu, Xu, Wr, Xr, Wd, Xd,
									Xface, Wface);
							// and if error_flag?

						} else {
							Wface = HALF*(Wl + Wr);

							dWdxl = Local_SolnBlk[nb].dWdx[i][j]; dWdxr = Local_SolnBlk[nb].dWdy[i+1][j];
							dWdyl = Local_SolnBlk[nb].dWdy[i][j]; dWdyr = Local_SolnBlk[nb].dWdy[i+1][j];

							if (IP->i_Viscous_Reconstruction == VISCOUS_RECONSTRUCTION_HYBRID) {
								Hybrid_Find_dWdX(dWdx, dWdy, 
										Xl, Wl, dWdxl, dWdyl,
										Xr, Wr, dWdxr, dWdyr,
										norm_dir);
							} else if (IP->i_Viscous_Reconstruction == VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE) {
								dWdx = HALF*(dWdxl + dWdxr);
								dWdy = HALF*(dWdyl + dWdyr);
							}
						}
// #ifdef _EB_PARALLEL_DEBUG_
// 		dout << endl << " ==========================================================";
// 		dout << endl << " EAST face of cell (" << i << "," << j <<  "," << Local_Solution_Block_List->Block[nb].gblknum << ") VISCOUS flux:";
// 		dout << endl << " Diamond-path = " << viscous_bc_flag;
// 		dout << endl << " xfaceE =" << Local_SolnBlk[nb].Grid.xfaceE(i,j);
// 		dout << endl << " nfaceE =" << Local_SolnBlk[nb].Grid.nfaceE(i,j);
// 		dout << endl << " lfaceE = " << Local_SolnBlk[nb].Grid.lfaceE(i,j);
// 		dout << endl << " Xl =" << Xl;
// 		dout << endl << " Xr =" << Xr;
// 		dout << endl << " Xu =" << Xu;
// 		dout << endl << " Xd =" << Xd;
// 		dout << endl << " Wl =" << Wl;
// 		dout << endl << " Wr =" << Wr;
// 		dout << endl << " Wu =" << Wu;
// 		dout << endl << " Wd =" << Wd;
// 		dout << endl << " F =" << ViscousFluxDiamondPath_n(Local_SolnBlk[nb].Grid.xfaceE(i,j),
// 								   Xl,Wl,Xd,Wd,Xr,Wr,Xu,Wu,
// 								   Local_SolnBlk[nb].Grid.nfaceE(i,j),
// 								   Local_SolnBlk[nb].Axisymmetric,
// 								   viscous_bc_flag);
// #endif

					}
				}
				// Compute the EAST face viscous flux.
				Flux -= ViscousFlux_n(Xface, Wface, dWdx, dWdy,
						norm_dir, Local_SolnBlk[nb].Axisymmetric, adiabatic_flag);
// 		dout << endl << " l:" << Xl << Wl;
// 		dout << endl << " r:" << Xr << Wr;
// 		dout << endl << " u:" << Xu << Wu;
// 		dout << endl << " d:" << Xd << Wd;
// 		dout << endl << " ndl = " << Vector2D((Xd.y-Xl.y),-(Xd.x-Xl.x));
// 		dout << endl << " nud = " << Vector2D((Xu.y-Xd.y),-(Xu.x-Xd.x));
// 		dout << endl << " nlu = " << Vector2D((Xl.y-Xu.y),-(Xl.x-Xu.x));
// 		dout << endl << " Al = " << HALF*((Xd-Xl)^(Xu-Xl));
// 		dout << endl << " nrd = " << Vector2D((Xr.y-Xd.y),-(Xr.x-Xd.x));
// 		dout << endl << " nur = " << Vector2D((Xu.y-Xr.y),-(Xu.x-Xr.x));
// 		dout << endl << " ndu = " << Vector2D((Xd.y-Xu.y),-(Xd.x-Xu.x));
// 		dout << endl << " Ar = " << HALF*((Xr-Xu)^(Xr-Xd));
// 		if (viscous_bc_flag == DIAMONDPATH_QUADRILATERAL_RECONSTRUCTION) {
// 		dout << endl << Bilinear_Interpolation(Wl,Xl,Wu,Xu,Wr,Xr,Wd,Xd,Local_SolnBlk[nb].Grid.xfaceE(i,j),dWdxl) << dWdxl;
// 		dout << endl << Bilinear_Interpolation(Wl,Xl,Wd,Xd,Wr,Xr,Wu,Xu,Local_SolnBlk[nb].Grid.xfaceE(i,j),dWdxl) << dWdxl;
// 		dout << endl << Bilinear_Interpolation_HC(Wl,Xl,Wu,Xu,Wr,Xr,Wd,Xd,Local_SolnBlk[nb].Grid.xfaceE(i,j),dWdxl) << dWdxl;
// 		dout << endl << Bilinear_Interpolation_HC(Wl,Xl,Wd,Xd,Wr,Xr,Wu,Xu,Local_SolnBlk[nb].Grid.xfaceE(i,j),dWdxl) << dWdxl;
// 		}
// 		dout << endl << " F =" << ViscousFluxDiamondPath_n(Local_SolnBlk[nb].Grid.xfaceE(i,j),
// 								   Xl,Wl,Xd,Wd,Xr,Wr,Xu,Wu,
// 								   Local_SolnBlk[nb].Grid.nfaceE(i,j),
// 								   Local_SolnBlk[nb].Axisymmetric,
// 								   viscous_bc_flag);
			} while (0); // end of viscous flux clause.

	    // Evaluate cell-averaged solution changes.
	    Local_SolnBlk[nb].dUdt[i  ][j][k_residual] -= (IP->CFL_Number*Local_SolnBlk[nb].dt[i][j])*
                                                          Flux*Local_SolnBlk[nb].Grid.lfaceE(i,j)/Local_SolnBlk[nb].Grid.Cell[i][j].A;
	    Local_SolnBlk[nb].dUdt[i+1][j][k_residual] += (IP->CFL_Number*Local_SolnBlk[nb].dt[i+1][j])*
	                                                  Flux*Local_SolnBlk[nb].Grid.lfaceW(i+1,j)/Local_SolnBlk[nb].Grid.Cell[i+1][j].A;

	    // Include all required source terms.
	    Local_SolnBlk[nb].dUdt[i][j][k_residual] += (IP->CFL_Number*Local_SolnBlk[nb].dt[i][j])*Local_SolnBlk[nb].W[i][j].S(Local_SolnBlk[nb].Grid.Cell[i][j].Xc,
																Local_SolnBlk[nb].dWdx[i][j],
																Local_SolnBlk[nb].dWdy[i][j],
																Local_SolnBlk[nb].Axisymmetric);

	    // Include area change source term.
	    if (Interface_Union_List.Ni)
	      Local_SolnBlk[nb].dUdt[i][j][k_residual] -= (IP->CFL_Number*Local_SolnBlk[nb].dt[i][j])*Local_SolnBlk[nb].U[i][j]*dAdt(nb,i,j,Time);

	    // Save west and east face boundary flux.
	    if (i == Local_SolnBlk[nb].ICl-1) {
	      Local_SolnBlk[nb].FluxW[j] = -Flux*Local_SolnBlk[nb].Grid.lfaceW(i+1,j);
	    } else if (i == Local_SolnBlk[nb].ICu) {
	      Local_SolnBlk[nb].FluxE[j] =  Flux*Local_SolnBlk[nb].Grid.lfaceE(i,j);
	    }

	  }
	}

	if (j > Local_SolnBlk[nb].JCl-1 && j < Local_SolnBlk[nb].JCu+1) {
	  Local_SolnBlk[nb].dUdt[Local_SolnBlk[nb].ICl-1][j][k_residual].Vacuum();
	  Local_SolnBlk[nb].dUdt[Local_SolnBlk[nb].ICu+1][j][k_residual].Vacuum();
	}

      }

      // Add j-direction (eta-direction) fluxes.
      for (int i = Local_SolnBlk[nb].ICl; i <= Local_SolnBlk[nb].ICu; i++) {
 	for (int j = Local_SolnBlk[nb].JCl-1; j <= Local_SolnBlk[nb].JCu; j++) {

 	  // Evaluate the cell interface j-direction INVISCID fluxes.
 	  if (Local_SolnBlk[nb].Grid.lfaceN(i,j) < TOLER) {
 	    // NORTH face of cell (i,j) has zero length.
 	    Wl.Standard_Atmosphere(); Wr.Standard_Atmosphere(); V = Vector2D_ZERO;

 	  } else if (Mesh[nb].cell_status[i][j  ] != CELL_STATUS_ACTIVE &&
 		     Mesh[nb].cell_status[i][j+1] != CELL_STATUS_ACTIVE) {
 	    // NORTH face of cell (i,j) is inactive.
 	    Wl.Standard_Atmosphere(); Wr.Standard_Atmosphere(); V = Vector2D_ZERO;

 	  } else if (Mesh[nb].cell_status[i][j  ] != CELL_STATUS_ACTIVE &&
 		     Mesh[nb].cell_status[i][j+1] == CELL_STATUS_ACTIVE) {
	    // SOUTH face of cell (i,j+1) corresponds to an embedded boundary.
	    dX = Local_SolnBlk[nb].Grid.xfaceS(i,j+1) - Local_SolnBlk[nb].Grid.Cell[i][j+1].Xc;
	    Wr = Local_SolnBlk[nb].W[i][j+1] + (Local_SolnBlk[nb].phi[i][j+1]^Local_SolnBlk[nb].dWdx[i][j+1])*dX.x +
	                                       (Local_SolnBlk[nb].phi[i][j+1]^Local_SolnBlk[nb].dWdy[i][j+1])*dX.y;
	    Ni = Mesh[nb].cell_status[i][j];
	    // Determine the boundary condition type.
	    Interface_BC_Type = Interface_Union_List[Ni].Determine_Interface_BC_Type(Local_SolnBlk[nb].Grid.xfaceS(i,j+1));
	    // Determine the boundary velocity.
	    V = Interface_Union_List[Ni].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceS(i,j+1),Time);
	    // Determine the left state by applying the appropriate
	    // boundary condition.
	    if (Interface_BC_Type == INTERFACE_BC_REFLECTION) {
	      // SOUTH face of cell (i,j+1) is a REFLECTION boundary.
	      Wl = Reflect(Wr,Local_SolnBlk[nb].Grid.nfaceS(i,j+1),V);
	    } else if (Interface_BC_Type == INTERFACE_BC_BURNING_SURFACE) {
	      // SOUTH face of cell (i,j+1) is a BURNING_SURFACE boundary.
	      Wl = BurningSurface(Wr,Local_SolnBlk[nb].Grid.nfaceS(i,j+1));
	    } else if (Interface_BC_Type == INTERFACE_BC_WALL_VISCOUS_HEATFLUX) {
	      // SOUTH face of cell (i,j+1) is a WALL_VISCOUS_HEATFLUX boundary.
	      Wl = WallViscousHeatFlux(Wr,Local_SolnBlk[nb].Grid.nfaceS(i,j+1),V);
	    } else if (Interface_BC_Type == INTERFACE_BC_WALL_VISCOUS_ISOTHERMAL) {
	      // SOUTH face of cell (i,j+1) is a WALL_VISCOUS_ISOTHERMAL boundary.
	      Wl = WallViscousIsothermal(Wr,Local_SolnBlk[nb].Grid.nfaceS(i,j+1),V,Local_SolnBlk[nb].Twall);
	    } else if (Interface_BC_Type == INTERFACE_BC_RINGLEB) {
	      // SOUTH face of cell (i,j+1) is a RINGLEB_FLOW boundary.
	      Wl = RinglebFlow(Wl,Local_SolnBlk[nb].Grid.Cell[i][j].Xc,q1,k);
	      Wl = RinglebFlow(Wl,Local_SolnBlk[nb].Grid.Cell[i][j+1].Xc,q2,k);
	      Wl = RinglebFlow(Wl,Local_SolnBlk[nb].Grid.xfaceS(i,j+1));
	      if (q1 < IP->Isotach_Line && q2 > IP->Isotach_Line)
		Wl = BC_Characteristic_Pressure(Wr,Wl,Local_SolnBlk[nb].Grid.nfaceS(i,j+1));
	    }
// #ifdef _EB_PARALLEL_DEBUG_
// 	    dout << endl << " ==========================================================";
// 	    dout << endl << " SOUTH face of cell (" << i << "," << j+1 <<  "," << Local_Solution_Block_List->Block[nb].gblknum << ") INVISCID flux:";
// 	    dout << endl << " Embedded boundary with bc = " << Interface_BC_Type;
// 	    dout << endl << " xfaceS =" << Local_SolnBlk[nb].Grid.xfaceS(i,j+1);
// 	    dout << endl << " nfaceS =" << Local_SolnBlk[nb].Grid.nfaceS(i,j+1);
// 	    dout << endl << " lfaceS = " << Local_SolnBlk[nb].Grid.lfaceS(i,j+1);
// 	    dout << endl << " Wl =" << Wl;
// 	    dout << endl << " Wr =" << Wr;
// 	    dout << endl << " F =" << FluxRoe_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
// #endif

 	  } else if (Mesh[nb].cell_status[i][j  ] == CELL_STATUS_ACTIVE &&
 		     Mesh[nb].cell_status[i][j+1] != CELL_STATUS_ACTIVE) {
	    // NORTH face of cell (i,j) corresponds to an embedded boundary.
	    dX = Local_SolnBlk[nb].Grid.xfaceN(i,j) - Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	    Wl = Local_SolnBlk[nb].W[i][j] + (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdx[i][j])*dX.x +
	                                     (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdy[i][j])*dX.y;
	    Ni = Mesh[nb].cell_status[i][j+1];
	    // Determine the boundary condition type.
	    Interface_BC_Type = Interface_Union_List[Ni].Determine_Interface_BC_Type(Local_SolnBlk[nb].Grid.xfaceN(i,j));
	    // Determine the boundary velocity.
	    V = Interface_Union_List[Ni].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceN(i,j),Time);
	    // Determine the right state by applying the appropriate
	    // boundary condition.
	    if (Interface_BC_Type == INTERFACE_BC_REFLECTION) {
	      // NORTH face of cell (i,j) is a REFLECTION boundary.
	      Wr = Reflect(Wl,Local_SolnBlk[nb].Grid.nfaceN(i,j),V);
	    } else if (Interface_BC_Type == INTERFACE_BC_BURNING_SURFACE) {
	      // NORTH face of cell (i,j) is a BURNING_SURFACE boundary.
	      Wr = BurningSurface(Wl,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    } else if (Interface_BC_Type == INTERFACE_BC_WALL_VISCOUS_HEATFLUX) {
	      // NORTH face of cell (i,j) is a WALL_VISCOUS_HEATFLUX boundary.
	      Wr = WallViscousHeatFlux(Wl,Local_SolnBlk[nb].Grid.nfaceN(i,j),V);
	    } else if (Interface_BC_Type == INTERFACE_BC_WALL_VISCOUS_ISOTHERMAL) {
	      // NORTH face of cell (i,j) is a WALL_VISCOUS_ISOTHERMAL boundary.
	      Wr = WallViscousIsothermal(Wl,Local_SolnBlk[nb].Grid.nfaceN(i,j),V,Local_SolnBlk[nb].Twall);
	    } else if (Interface_BC_Type == INTERFACE_BC_RINGLEB) {
	      // NORTH face of cell (i,j) is a RINGLEB_FLOW boundary.
	      Wr = RinglebFlow(Wr,Local_SolnBlk[nb].Grid.Cell[i][j+1].Xc,q1,k);
	      Wr = RinglebFlow(Wr,Local_SolnBlk[nb].Grid.Cell[i][j].Xc,q2,k);
	      Wr = RinglebFlow(Wr,Local_SolnBlk[nb].Grid.xfaceN(i,j));
	      if (q1 < IP->Isotach_Line && q2 > IP->Isotach_Line)
		Wr = BC_Characteristic_Pressure(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    }
// #ifdef _EB_PARALLEL_DEBUG_
// 	    dout << endl << " ==========================================================";
// 	    dout << endl << " NORTH face of cell (" << i << "," << j <<  "," << Local_Solution_Block_List->Block[nb].gblknum << ") INVISCID flux:";
// 	    dout << endl << " Embedded boundary with bc = " << Interface_BC_Type;
// 	    dout << endl << " xfaceN =" << Local_SolnBlk[nb].Grid.xfaceN(i,j);
// 	    dout << endl << " nfaceN =" << Local_SolnBlk[nb].Grid.nfaceN(i,j);
// 	    dout << endl << " lfaceN = " << Local_SolnBlk[nb].Grid.lfaceN(i,j);
// 	    dout << endl << " Wl =" << Wl;
// 	    dout << endl << " Wr =" << Wr;
// 	    dout << endl << " F =" << FluxRoe_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
// #endif

 	  } else if (Mesh[nb].cell_status[i][j  ] == CELL_STATUS_ACTIVE &&
 		     Mesh[nb].cell_status[i][j+1] == CELL_STATUS_ACTIVE) {

 	    V = Vector2D_ZERO;

 	    if (j == Local_SolnBlk[nb].JCl-1 && 
 		(Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_REFLECTION ||
 		 Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_BURNING_SURFACE ||
 		 Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX ||
 		 Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
 		 Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_MOVING_WALL_HEATFLUX ||
 		 Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_MOVING_WALL_ISOTHERMAL ||
 		 Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_RINGLEB_FLOW ||
 		 Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_CHARACTERISTIC)) {

	      dX = Local_SolnBlk[nb].Grid.xfaceS(i,j+1) - Local_SolnBlk[nb].Grid.Cell[i][j+1].Xc;
	      Wr = Local_SolnBlk[nb].W[i][j+1] + (Local_SolnBlk[nb].phi[i][j+1]^Local_SolnBlk[nb].dWdx[i][j+1])*dX.x +
   	                                         (Local_SolnBlk[nb].phi[i][j+1]^Local_SolnBlk[nb].dWdy[i][j+1])*dX.y;

	      // SOUTH face of cell (i,j+1) is a normal boundary.
	      if (Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_REFLECTION) {
		// SOUTH face of cell (i,j+1) is a REFLECTION boundary.
		Wl = Reflect(Wr,Local_SolnBlk[nb].Grid.nfaceS(i,j+1));
	      } else if (Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_BURNING_SURFACE) {
		// SOUTH face of cell (i,j+1) is a BURNING_SURFACE boundary.
		Wl = BurningSurface(Wr,Local_SolnBlk[nb].Grid.nfaceS(i,j+1));
	      } else if (Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX) {
		// SOUTH face of cell (i,j+1) is a WALL_VISCOUS_HEATFLUX boundary.
		Wl = WallViscousHeatFlux(Wr,Local_SolnBlk[nb].Grid.nfaceS(i,j+1));
	      } else if (Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL) {
		// SOUTH face of cell (i,j+1) is a WALL_VISCOUS_ISOTHERMAL boundary.
		Wl = WallViscousIsothermal(Wr,Local_SolnBlk[nb].Grid.nfaceS(i,j+1),Local_SolnBlk[nb].Twall);
	      } else if (Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_MOVING_WALL_HEATFLUX) {
		// SOUTH face of cell (i,j+1) is a MOVINGWALL_HEATFLUX boundary.
		Wl = MovingWallHeatFlux(Wr,Local_SolnBlk[nb].Grid.nfaceS(i,j+1),Local_SolnBlk[nb].Vwall.x);
	      } else if (Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_MOVING_WALL_ISOTHERMAL) {
		// SOUTH face of cell (i,j+1) is a MOVINGWALL_ISOTHERMAL boundary.
		Wl = MovingWallIsothermal(Wr,Local_SolnBlk[nb].Grid.nfaceS(i,j+1),Local_SolnBlk[nb].Vwall.x,Local_SolnBlk[nb].Twall);
	      } else if (Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_RINGLEB_FLOW) {
		// SOUTH face of cell (i,j+1) is a RINGLEB_FLOW boundary.
		Wl = RinglebFlow(Wl,Local_SolnBlk[nb].Grid.xfaceS(i,j+1));
		Wl = BC_Characteristic_Pressure(Wr,Wl,Local_SolnBlk[nb].Grid.nfaceS(i,j+1));
	      } else if (Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_CHARACTERISTIC) {
		// SOUTH face of cell (i,j+1) is a CHARACTERISTIC boundary.
		Wl = BC_Characteristic_Pressure(Wr,Local_SolnBlk[nb].WoS[i],Local_SolnBlk[nb].Grid.nfaceS(i,j+1));
	      }
// #ifdef _EB_PARALLEL_DEBUG_
// 	    dout << endl << " ==========================================================";
// 	    dout << endl << " SOUTH face of cell (" << i << "," << j+1 <<  "," << Local_Solution_Block_List->Block[nb].gblknum << ") INVISCID flux:";
// 	    dout << endl << " Normal boundary with bc = " << Local_SolnBlk[nb].Grid.BCtypeS[i];
// 	    dout << endl << " xfaceS =" << Local_SolnBlk[nb].Grid.xfaceS(i,j+1);
// 	    dout << endl << " nfaceS =" << Local_SolnBlk[nb].Grid.nfaceS(i,j+1);
// 	    dout << endl << " lfaceS = " << Local_SolnBlk[nb].Grid.lfaceS(i,j+1);
// 	    dout << endl << " Wl =" << Wl;
// 	    dout << endl << " Wr =" << Wr;
// 	    dout << endl << " F =" << FluxRoe_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
// #endif

	    } else if (j == Local_SolnBlk[nb].JCu && 
		       (Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_REFLECTION ||
			Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_BURNING_SURFACE ||
			Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX ||
			Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
			Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_MOVING_WALL_HEATFLUX ||
			Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_MOVING_WALL_ISOTHERMAL ||
			Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_RINGLEB_FLOW ||
			Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_CHARACTERISTIC)) {

	      dX = Local_SolnBlk[nb].Grid.xfaceN(i,j) - Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	      Wl = Local_SolnBlk[nb].W[i][j] + (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdx[i][j])*dX.x +
   	                                       (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdy[i][j])*dX.y;

	      // NORTH face of cell (i,j) is a normal boundary.
	      if (Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_REFLECTION) {
		// NORTH face of cell (i,j) is a REFLECTION boundary.
		Wr = Reflect(Wl,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	      } else if (Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_BURNING_SURFACE) {
		// NORTH face of cell (i,j) is a BURNING_SURFACE boundary.
		Wr = BurningSurface(Wl,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	      } else if (Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX) {
		// NORTH face of cell (i,j) is a WALL_VISCOUS_HEATFLUX boundary.
		Wr = WallViscousHeatFlux(Wl,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	      } else if (Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL) {
		// NORTH face of cell (i,j) is a WALL_VISCOUS_ISOTHERMAL boundary.
		Wr = WallViscousIsothermal(Wl,Local_SolnBlk[nb].Grid.nfaceN(i,j),Local_SolnBlk[nb].Twall);
	      } else if (Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_MOVING_WALL_HEATFLUX) {
		// NORTH face of cell (i,j) is a MOVINGWALL_HEATFLUX boundary.
		Wr = MovingWallHeatFlux(Wl,Local_SolnBlk[nb].Grid.nfaceN(i,j),Local_SolnBlk[nb].Vwall.x);
	      } else if (Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_MOVING_WALL_ISOTHERMAL) {
		// NORTH face of cell (i,j) is a MOVINGWALL_ISOTHERMAL boundary.
		Wr = MovingWallIsothermal(Wl,Local_SolnBlk[nb].Grid.nfaceN(i,j),Local_SolnBlk[nb].Vwall.x,Local_SolnBlk[nb].Twall);
	      } else if (Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_RINGLEB_FLOW) {
		// NORTH face of cell (i,j) is a RINGLEB_FLOW boundary.
		Wr = RinglebFlow(Wr,Local_SolnBlk[nb].Grid.xfaceN(i,j));
		Wr = BC_Characteristic_Pressure(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	      } else {
		// NORTH face of cell (i,j) is a CHARACTERISTIC boundary.
		Wr = BC_Characteristic_Pressure(Wl,Local_SolnBlk[nb].WoN[i],Local_SolnBlk[nb].Grid.nfaceN(i,j));
	      }
// #ifdef _EB_PARALLEL_DEBUG_
// 	    dout << endl << " ==========================================================";
// 	    dout << endl << " NORTH face of cell (" << i << "," << j <<  "," << Local_Solution_Block_List->Block[nb].gblknum << ") INVISCID flux:";
// 	    dout << endl << " Normal boundary with bc = " << Local_SolnBlk[nb].Grid.BCtypeN[i];
// 	    dout << endl << " xfaceN =" << Local_SolnBlk[nb].Grid.xfaceN(i,j);
// 	    dout << endl << " nfaceN =" << Local_SolnBlk[nb].Grid.nfaceN(i,j);
// 	    dout << endl << " lfaceN = " << Local_SolnBlk[nb].Grid.lfaceN(i,j);
// 	    dout << endl << " Wl =" << Wl;
// 	    dout << endl << " Wr =" << Wr;
// 	    dout << endl << " F =" << FluxRoe_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
// #endif

	    } else {

	      // NORTH face is either a normal cell or possibly a FIXED, 
	      // NONE or EXTRAPOLATION boundary.
	      dX = Local_SolnBlk[nb].Grid.xfaceN(i,j  ) - Local_SolnBlk[nb].Grid.Cell[i][j  ].Xc;
	      Wl = Local_SolnBlk[nb].W[i][j  ] + (Local_SolnBlk[nb].phi[i][j  ]^Local_SolnBlk[nb].dWdx[i][j  ])*dX.x +
                                                 (Local_SolnBlk[nb].phi[i][j  ]^Local_SolnBlk[nb].dWdy[i][j  ])*dX.y;
	      dX = Local_SolnBlk[nb].Grid.xfaceS(i,j+1) - Local_SolnBlk[nb].Grid.Cell[i][j+1].Xc;
	      Wr = Local_SolnBlk[nb].W[i][j+1] + (Local_SolnBlk[nb].phi[i][j+1]^Local_SolnBlk[nb].dWdx[i][j+1])*dX.x +
                                                 (Local_SolnBlk[nb].phi[i][j+1]^Local_SolnBlk[nb].dWdy[i][j+1])*dX.y;
// #ifdef _EB_PARALLEL_DEBUG_
// 	    dout << endl << " ==========================================================";
// 	    dout << endl << " NORTH face of cell (" << i << "," << j <<  "," << Local_Solution_Block_List->Block[nb].gblknum << ") INVISCID flux:";
// 	    dout << endl << " xfaceN =" << Local_SolnBlk[nb].Grid.xfaceN(i,j);
// 	    dout << endl << " nfaceN =" << Local_SolnBlk[nb].Grid.nfaceN(i,j);
// 	    dout << endl << " lfaceN = " << Local_SolnBlk[nb].Grid.lfaceN(i,j);
// 	    dout << endl << " Wl =" << Wl;
// 	    dout << endl << " Wr =" << Wr;
// 	    dout << endl << " F =" << FluxRoe_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
// #endif

	    }

	  }

// 	  dout << endl << " NORTH";
// 	  dout << endl << Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
// 	  dout << endl << Wl;
// 	  dout << endl << Wr;
// 	  dout << endl << FluxRoe_MB_n(Wl,Wr,V,Local_SolnBlk[nb].Grid.nfaceN(i,j));

	  // Determine NORTH face inviscid flux.
	  switch(IP->i_Flux_Function) {
	  case FLUX_FUNCTION_GODUNOV :
	    Flux = FluxGodunov_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_ROE :
	    Flux = FluxRoe_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_RUSANOV :
	    Flux = FluxRusanov_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_HLLE :
	    Flux = FluxHLLE_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_HLLL :
	    Flux = FluxHLLL_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_HLLC :
	    Flux = FluxHLLC_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_VANLEER :
	    Flux = FluxVanLeer_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_AUSM :
	    Flux = FluxAUSM_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_AUSMplus :
	    Flux = FluxAUSMplus_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_GODUNOV_MB :
	    Flux = FluxGodunov_MB_n(Wl,Wr,V,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_ROE_MB :
	    Flux = FluxRoe_MB_n(Wl,Wr,V,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_HLLE_MB :
	    Flux = FluxHLLE_MB_n(Wl,Wr,V,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_VANLEER_MB :
	    Flux = FluxVanLeer_MB_n(Wl,Wr,V,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  default:
	    Flux = FluxRoe_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  };

	    // Determine the j-direction viscous flux.

			do { // with while (0) to allow jumps to the end of the viscous flux calculation.
	    	if (Local_SolnBlk[nb].Flow_Type == 0) { break; }

				if (Local_SolnBlk[nb].Grid.lfaceN(i,j) < TOLER) { break; }

				if (Mesh[nb].cell_status[i][j  ] != CELL_STATUS_ACTIVE &&
						Mesh[nb].cell_status[i][j+1] != CELL_STATUS_ACTIVE) {
					break;
				}

				// Please rotate your head 90 degrees to the left before
				// attempting to read this code. For example, Xr is the
				// position of the cell to the north of the face in
				// question. This is necessary, for example, for
				// DiamondPath_Find_dWdX, which only understands "left" and
				// "right" triangles.

				int adiabatic_flag = OFF;
				Vector2D Xface = Local_SolnBlk[nb].Grid.xfaceN(i, j);
				Vector2D norm_dir = Local_SolnBlk[nb].Grid.nfaceN(i, j);

	      if (Mesh[nb].cell_status[i][j  ] != CELL_STATUS_ACTIVE &&
				    Mesh[nb].cell_status[i][j+1] == CELL_STATUS_ACTIVE) {

					// SOUTH face of cell (i,j+1) corresponds to an embedded boundary.
					Xr = Local_SolnBlk[nb].Grid.Cell[i][j+1].Xc; 
					Wr = Local_SolnBlk[nb].W[i][j+1];

					Ni = Mesh[nb].cell_status[i][j];
					Interface_BC_Type = Interface_Union_List[Ni].Determine_Interface_BC_Type(Xface);
					
					// No viscous flux for a reflected boundary. Is this true?
					if (Interface_BC_Type == INTERFACE_BC_REFLECTION) { break; } 

					switch (Interface_BC_Type) {
						// The face velocity is the same regardless of the BC type 
						// and is set below.
						case INTERFACE_BC_WALL_VISCOUS_HEATFLUX:
							Wface.rho = Wr.rho; 
							Wface.p = Wr.p; Wface.k = Wr.k; Wface.omega = Wr.omega;
							break;
						case INTERFACE_BC_WALL_VISCOUS_ISOTHERMAL:
							Wface.rho = Wr.p/(Wr.R*Local_SolnBlk[nb].Twall); 
							Wface.p = Wr.p; Wface.k = Wr.k; Wface.omega = Wr.omega;
							break;
						case INTERFACE_BC_BURNING_SURFACE:
							Wface = BurningSurface(Wr, norm_dir);
							break;
					}

					Wface.v = Interface_Union_List[Ni].Determine_Interface_Velocity(Xface, Time);

					if (Interface_BC_Type == INTERFACE_BC_WALL_VISCOUS_HEATFLUX) {
						adiabatic_flag = ON;
					}

					switch(IP->i_Viscous_Reconstruction) {
						case VISCOUS_RECONSTRUCTION_CARTESIAN :
						case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
							{
							int stencil_flag = DIAMONDPATH_RIGHT_TRIANGLE;
							Xd = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
							Wd = Wface;
							Wd.v = Interface_Union_List[Ni].Determine_Interface_Velocity(Xd, Time);

							Xu = Local_SolnBlk[nb].Grid.Node[i  ][j+1].X;
							Wu = Wface;
							Wu.v = Interface_Union_List[Ni].Determine_Interface_Velocity(Xu, Time);

							DiamondPath_Find_dWdX(dWdx, dWdy,
									BlankVector, BlankState, 
									Xd, Wd,
									Xr, Wr,
									Xu, Wu,
									stencil_flag);
							}
							break;
						case VISCOUS_RECONSTRUCTION_HYBRID :
							dWdxr = Local_SolnBlk[nb].dWdx[i][j+1]; 
							dWdyr = Local_SolnBlk[nb].dWdy[i][j+1];

							Xl = Xface;
							Wl = Wface;
							dWdxl = dWdxr; dWdyl = dWdyr; // Hmm. 

							Hybrid_Find_dWdX(dWdx, dWdy, 
									Xl, Wl, dWdxl, dWdyl,
									Xr, Wr, dWdxr, dWdyr,
									norm_dir);
							break;
						case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
							dWdx = Local_SolnBlk[nb].dWdx[i][j+1]; 
							dWdy = Local_SolnBlk[nb].dWdy[i][j+1];
							break;
					}
// #ifdef _EB_PARALLEL_DEBUG_
// 	      dout << endl << " ==========================================================";
// 	      dout << endl << " SOUTH face of cell (" << i << "," << j+1 << "," << Local_Solution_Block_List->Block[nb].gblknum << ") VISCOUS flux:";
// 	      dout << endl << " Embbedded boundary with bc = " << Interface_BC_Type;
// 	      dout << endl << " Diamond-path = " << viscous_bc_flag;
// 	      dout << endl << " xfaceS =" << Local_SolnBlk[nb].Grid.xfaceS(i,j+1);
// 	      dout << endl << " nfaceS =" << Local_SolnBlk[nb].Grid.nfaceS(i,j+1);
// 	      dout << endl << " lfaceS = " << Local_SolnBlk[nb].Grid.lfaceS(i,j+1);
// 	      dout << endl << " Xr =" << Xr;
// 	      dout << endl << " Xu =" << Xu;
// 	      dout << endl << " Xd =" << Xd;
// 	      dout << endl << " Wr =" << Wr;
// 	      dout << endl << " Wu =" << Wu;
// 	      dout << endl << " Wd =" << Wd;
// 	      dout << endl << " F =" << ViscousFluxDiamondPath_n(Local_SolnBlk[nb].Grid.xfaceN(i,j),
// 								 Xl,Wl,Xd,Wd,Xr,Wr,Xu,Wu,
// 								 Local_SolnBlk[nb].Grid.nfaceN(i,j),
// 								 Local_SolnBlk[nb].Axisymmetric,
// 								 viscous_bc_flag);
// #endif
				} else if (Mesh[nb].cell_status[i][j  ] == CELL_STATUS_ACTIVE &&
				           Mesh[nb].cell_status[i][j+1] != CELL_STATUS_ACTIVE) {

					// NORTH face of cell (i,j) corresponds to an embedded boundary.
					Xl = Local_SolnBlk[nb].Grid.Cell[i][j].Xc; 
					Wl = Local_SolnBlk[nb].W[i][j];

					Ni = Mesh[nb].cell_status[i][j+1];
					Interface_BC_Type = Interface_Union_List[Ni].Determine_Interface_BC_Type(Xface);

					// No viscous flux for a reflected boundary. Is this true?
					if (Interface_BC_Type == INTERFACE_BC_REFLECTION) { break; } 

					switch (Interface_BC_Type) {
						// The face velocity is the same regardless of the BC type 
						// and is set below.
						case INTERFACE_BC_WALL_VISCOUS_HEATFLUX:
							Wface.rho = Wl.rho; 
							Wface.p = Wl.p; Wface.k = Wl.k; Wface.omega = Wl.omega;
							break;
						case INTERFACE_BC_WALL_VISCOUS_ISOTHERMAL:
							Wface.rho = Wl.p/(Wl.R*IP->Twall); 
							Wface.p = Wl.p; Wface.k = Wl.k; Wface.omega = Wl.omega;
							break;
						case INTERFACE_BC_BURNING_SURFACE:
							Wface = BurningSurface(Wl, norm_dir);
							break;
					}

					Wface.v = Interface_Union_List[Ni].Determine_Interface_Velocity(Xface, Time);

					if (Interface_BC_Type == INTERFACE_BC_WALL_VISCOUS_HEATFLUX) {
						adiabatic_flag = ON;
					}

					switch(IP->i_Viscous_Reconstruction) {
						case VISCOUS_RECONSTRUCTION_CARTESIAN :
						case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
							{
							int stencil_flag = DIAMONDPATH_LEFT_TRIANGLE;
							Xd = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
							Wd = Wface;
							Wd.v = Interface_Union_List[Ni].Determine_Interface_Velocity(Xd, Time);

							Xu = Local_SolnBlk[nb].Grid.Node[i  ][j+1].X;
							Wu = Wface;
							Wu.v = Interface_Union_List[Ni].Determine_Interface_Velocity(Xu, Time);

							DiamondPath_Find_dWdX(dWdx, dWdy,
									Xl, Wl,
									Xd, Wd,
									BlankVector, BlankState,
									Xu, Wu,
									stencil_flag);
							}
							break;
						case VISCOUS_RECONSTRUCTION_HYBRID :
							dWdxl = Local_SolnBlk[nb].dWdx[i][j]; 
							dWdyl = Local_SolnBlk[nb].dWdy[i][j]; 

							Xr = Xface;
							Wr = Wface;
							dWdxr = dWdxl; dWdyr = dWdyl; // Hmm.

							Hybrid_Find_dWdX(dWdx, dWdy, 
									Xl, Wl, dWdxl, dWdyl,
									Xr, Wr, dWdxr, dWdyr,
									norm_dir);
							break;
						case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
							dWdx = Local_SolnBlk[nb].dWdx[i][j]; 
							dWdy = Local_SolnBlk[nb].dWdy[i][j];
							break;
					}
// #ifdef _EB_PARALLEL_DEBUG_
// 	      dout << endl << " ==========================================================";
// 	      dout << endl << " NORTH face of cell (" << i << "," << j << "," << Local_Solution_Block_List->Block[nb].gblknum << ") VISCOUS flux:";
// 	      dout << endl << " Embbedded boundary with bc = " << Interface_BC_Type;
// 	      dout << endl << " Diamond-path = " << viscous_bc_flag;
// 	      dout << endl << " xfaceN =" << Local_SolnBlk[nb].Grid.xfaceN(i,j);
// 	      dout << endl << " nfaceN =" << Local_SolnBlk[nb].Grid.nfaceN(i,j);
// 	      dout << endl << " lfaceN = " << Local_SolnBlk[nb].Grid.lfaceN(i,j);
// 	      dout << endl << " Xl =" << Xl;
// 	      dout << endl << " Xu =" << Xu;
// 	      dout << endl << " Xd =" << Xd;
// 	      dout << endl << " Wl =" << Wl;
// 	      dout << endl << " Wu =" << Wu;
// 	      dout << endl << " Wd =" << Wd;
// 	      dout << endl << " F =" << ViscousFluxDiamondPath_n(Local_SolnBlk[nb].Grid.xfaceN(i,j),
// 								 Xl,Wl,Xd,Wd,Xr,Wr,Xu,Wu,
// 								 Local_SolnBlk[nb].Grid.nfaceN(i,j),
// 								 Local_SolnBlk[nb].Axisymmetric,
// 								 viscous_bc_flag);
// #endif
				} else if (Mesh[nb].cell_status[i][j  ] == CELL_STATUS_ACTIVE &&
				           Mesh[nb].cell_status[i][j+1] == CELL_STATUS_ACTIVE) {

					if (j == Local_SolnBlk[nb].JCl-1 && 
							(Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX ||
							 Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
							 Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_MOVING_WALL_HEATFLUX ||
							 Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_MOVING_WALL_ISOTHERMAL ||
							 Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_BURNING_SURFACE)) {

						// SOUTH face of cell (i,j+1) is a normal boundary.

						Xr = Local_SolnBlk[nb].Grid.Cell[i][j+1].Xc; 
						Wr = Local_SolnBlk[nb].W[i][j+1];

						if (Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX ||
								Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_MOVING_WALL_HEATFLUX) {
							adiabatic_flag = ON;
						}

						switch (Local_SolnBlk[nb].Grid.BCtypeS[i]) {
							case BC_WALL_VISCOUS_HEATFLUX:
								Wface.rho = Wr.rho; Wface.v.x = ZERO; Wface.v.y = ZERO; 
								Wface.p = Wr.p; Wface.k = Wr.k; Wface.omega =	Wr.omega;
								break;
							case BC_WALL_VISCOUS_ISOTHERMAL:
								Wface.rho = Wr.p/(Wr.R*Local_SolnBlk[nb].Twall); Wface.v.x = ZERO; Wface.v.y = ZERO; 
								Wface.p = Wr.p; Wface.k = Wr.k; Wface.omega = Wr.omega;
								break;
							case BC_MOVING_WALL_HEATFLUX:
								Wface.rho = Wr.rho; Wface.v.x = Local_SolnBlk[nb].Vwall.x; Wface.v.y = Local_SolnBlk[nb].Vwall.y; 
								Wface.p = Wr.p; Wface.k = Wr.k; Wface.omrga = Wr.omega;
								break;
							case BC_MOVING_WALL_ISOTHERMAL:
								Wface.rho = Wr.p/(Wr.R*Local_SolnBlk[nb].Twall); Wface.v.x = Local_SolnBlk[nb].Vwall.x; Wface.v.y = Local_SolnBlk[nb].Vwall.y; 
								Wface.p = Wr.p; Wface.k = Wr.k; Wface.omega = Wr.omega;
								break;
							case BC_BURNING_SURFACE:
								Wface = BurningSurface(Wr, norm_dir);
								break;
						}

						switch(IP->i_Viscous_Reconstruction) {
							case VISCOUS_RECONSTRUCTION_CARTESIAN :
							case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
								{
								int stencil_flag = DIAMONDPATH_RIGHT_TRIANGLE;
								Xu = Local_SolnBlk[nb].Grid.Node[i  ][j+1].X;
								Xd = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X; 
								DiamondPath_Find_dWdX(dWdx, dWdy,
										BlankVector, BlankState, 
										Xd, Wface,
										Xr, Wr,
										Xu, Wface,
										stencil_flag);
								}
								break;
							case VISCOUS_RECONSTRUCTION_HYBRID :
								dWdxr = Local_SolnBlk[nb].dWdx[i][j+1];
								dWdyr = Local_SolnBlk[nb].dWdy[i][j+1];

								Xl = Xface;
								Wl = Wface;
								dWdxl = dWdxr; dWdyl = dWdyr; // Hmm.

								Hybrid_Find_dWdX(dWdx, dWdy, 
										Xl, Wl, dWdxl, dWdyl,
										Xr, Wr, dWdxr, dWdyr,
										norm_dir);
								break;
							case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
								dWdx = Local_SolnBlk[nb].dWdx[i][j+1];
								dWdy = Local_SolnBlk[nb].dWdy[i][j+1];
								break;
						}
// #ifdef _EB_PARALLEL_DEBUG_
// 		dout << endl << " ==========================================================";
// 		dout << endl << " SOUTH face of cell (" << i << "," << j+1 << "," << Local_Solution_Block_List->Block[nb].gblknum << ") VISCOUS flux:";
// 		dout << endl << " Normal boundary with bc = " << Local_SolnBlk[nb].Grid.BCtypeS[i];
// 		dout << endl << " Diamond-path = " << viscous_bc_flag;
// 		dout << endl << " xfaceS =" << Local_SolnBlk[nb].Grid.xfaceS(i,j+1);
// 		dout << endl << " nfaceS =" << Local_SolnBlk[nb].Grid.nfaceS(i,j+1);
// 		dout << endl << " lfaceS = " << Local_SolnBlk[nb].Grid.lfaceS(i,j+1);
// 		dout << endl << " Xr =" << Xr;
// 		dout << endl << " Xu =" << Xu;
// 		dout << endl << " Xd =" << Xd;
// 		dout << endl << " Wr =" << Wr;
// 		dout << endl << " Wu =" << Wu;
// 		dout << endl << " Wd =" << Wd;
// 		dout << endl << " F =" << ViscousFluxDiamondPath_n(Local_SolnBlk[nb].Grid.xfaceN(i,j),
// 								   Xl,Wl,Xd,Wd,Xr,Wr,Xu,Wu,
// 								   Local_SolnBlk[nb].Grid.nfaceN(i,j),
// 								   Local_SolnBlk[nb].Axisymmetric,
// 								   viscous_bc_flag);
// #endif
					} else if (j == Local_SolnBlk[nb].JCu &&
							(Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX ||
							 Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
							 Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_MOVING_WALL_HEATFLUX ||
							 Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_MOVING_WALL_ISOTHERMAL ||
							 Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_BURNING_SURFACE)) {
						
						// NORTH face of cell (i,j) is a normal boundary.
						
						Xl = Local_SolnBlk[nb].Grid.Cell[i][j].Xc; 
						Wl = Local_SolnBlk[nb].W[i][j];

						if (Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX ||
						    Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_MOVING_WALL_HEATFLUX) {
							adiabatic_flag = ON;
						}

						switch (Local_SolnBlk[nb].Grid.BCtypeN[i]) {
							case BC_WALL_VISCOUS_HEATFLUX:
								Wface.rho = Wl.rho; Wface.v.x = ZERO; Wface.v.y = ZERO; 
								Wface.p = Wl.p; Wface.k = Wl.k; Wface.omega = Wl.omega;
								break;
							case BC_WALL_VISCOUS_ISOTHERMAL:
								Wface.rho = Wl.p/(Wl.R*Local_SolnBlk[nb].Twall); Wface.v.x = ZERO; Wface.v.y = ZERO; 
								Wface.p = Wl.p; Wface.k = Wl.k; Wface.omega = Wl.omega;
								break;
							case BC_MOVING_WALL_HEATFLUX:
								Wface.rho = Wl.rho; Wface.v.x = Local_SolnBlk[nb].Vwall.x; Wface.v.y = Local_SolnBlk[nb].Vwall.y; 
								Wface.p = Wl.p; Wface.k = Wl.k; Wface.omega = Wl.omega;
								break;
							case BC_MOVING_WALL_ISOTHERMAL:
								Wface.rho = Wl.p/(Wl.R*Local_SolnBlk[nb].Twall); Wface.v.x = Local_SolnBlk[nb].Vwall.x; Wface.v.y = Local_SolnBlk[nb].Vwall.y; 
								Wface.p = Wl.p; Wface.k = Wl.k; Wface.omega = Wl.omega;
								break;
							case BC_MOVING_WALL_ISOTHERMAL:
								Wface = BurningSurface(Wl, norm_dir);
								break;
						}

						switch(IP->i_Viscous_Reconstruction) {
							case VISCOUS_RECONSTRUCTION_CARTESIAN :
							case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
								{
								int stencil_flag = DIAMONDPATH_LEFT_TRIANGLE;
								Xu = Local_SolnBlk[nb].Grid.Node[i  ][j+1].X;
								Xd = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X; 
								DiamondPath_Find_dWdX(dWdx, dWdy,
										Xl, Wl,
										Xd, Wface,
										BlankVector, BlankState,
										Xu, Wface,
										stencil_flag);
								}
								break;
							case VISCOUS_RECONSTRUCTION_HYBRID :
								dWdxl = Local_SolnBlk[nb].dWdx[i][j]; 
								dWdyl = Local_SolnBlk[nb].dWdy[i][j]; 

								Xr = Xface; 
								Wr = Wface;
								dWdxr = dWdxl; dWdyr = dWdyl; 

								Hybrid_Find_dWdX(dWdx, dWdy, 
										Xl, Wl, dWdxl, dWdyl,
										Xr, Wr, dWdxr, dWdyr,
										norm_dir);
								break;
							case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
								dWdx = Local_SolnBlk[nb].dWdx[i][j]; 
								dWdy = Local_SolnBlk[nb].dWdy[i][j]; 
								break;
						}
// #ifdef _EB_PARALLEL_DEBUG_
// 	      dout << endl << " ==========================================================";
// 	      dout << endl << " NORTH face of cell (" << i << "," << j << "," << Local_Solution_Block_List->Block[nb].gblknum << ") VISCOUS flux:";
// 	      dout << endl << " Embbedded boundary with bc = " << Interface_BC_Type;
// 	      dout << endl << " Diamond-path = " << viscous_bc_flag;
// 	      dout << endl << " xfaceN =" << Local_SolnBlk[nb].Grid.xfaceN(i,j);
// 	      dout << endl << " nfaceN =" << Local_SolnBlk[nb].Grid.nfaceN(i,j);
// 	      dout << endl << " lfaceN = " << Local_SolnBlk[nb].Grid.lfaceN(i,j);
// 	      dout << endl << " Xl =" << Xl;
// 	      dout << endl << " Xu =" << Xu;
// 	      dout << endl << " Xd =" << Xd;
// 	      dout << endl << " Wl =" << Wl;
// 	      dout << endl << " Wu =" << Wu;
// 	      dout << endl << " Wd =" << Wd;
// 	      dout << endl << " F =" << ViscousFluxDiamondPath_n(Local_SolnBlk[nb].Grid.xfaceN(i,j),
// 								 Xl,Wl,Xd,Wd,Xr,Wr,Xu,Wu,
// 								 Local_SolnBlk[nb].Grid.nfaceN(i,j),
// 								 Local_SolnBlk[nb].Axisymmetric,
// 								 viscous_bc_flag);
// #endif
					} else {
						// EAST face is either a normal cell or possibly a non-viscous
						// boundary condition.
						Xl = Local_SolnBlk[nb].Grid.Cell[i][j  ].Xc; Wl = Local_SolnBlk[nb].W[i][j  ];
						Xr = Local_SolnBlk[nb].Grid.Cell[i][j+1].Xc; Wr = Local_SolnBlk[nb].W[i][j+1];

						if (IP->i_Viscous_Reconstruction == VISCOUS_RECONSTRUCTION_CARTESIAN ||
						    IP->i_Viscous_Reconstruction == VISCOUS_RECONSTRUCTION_DIAMOND_PATH) {

							int stencil_flag = DIAMONDPATH_QUADRILATERAL_RECONSTRUCTION;

							Xu = Local_SolnBlk[nb].Grid.Node[i  ][j+1].X; 
							Xd = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X; 

							// These two (WnNW and WnNE) are what make diamond path so expensive.
							Wu = Local_SolnBlk[nb].WnNW(i,j);
							Wd = Local_SolnBlk[nb].WnNE(i,j);

							DiamondPath_Find_dWdX(dWdx, dWdy,
									Xl, Wl, Xd, Wd, Xr, Wr, Xu, Wu, 
									stencil_flag);

							// Find Wface:
							error_flag = Bilinear_Interpolation(Wl, Xl, Wu, Xu, Wr, Xr, Wd, Xd,
									Xface, Wface);
							// and if error_flag?

						} else {
							Wface = HALF*(Wl + Wr);

							dWdxl = Local_SolnBlk[nb].dWdx[i][j]; dWdxr = Local_SolnBlk[nb].dWdy[i][j+1];
							dWdyl = Local_SolnBlk[nb].dWdy[i][j]; dWdyr = Local_SolnBlk[nb].dWdy[i][j+1];

							if (IP->i_Viscous_Reconstruction == VISCOUS_RECONSTRUCTION_HYBRID) {
								Hybrid_Find_dWdX(dWdx, dWdy, 
										Xl, Wl, dWdxl, dWdyl,
										Xr, Wr, dWdxr, dWdyr,
										norm_dir);
							} else if (IP->i_Viscous_Reconstruction == VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE) {
								dWdx = HALF*(dWdxl + dWdxr);
								dWdy = HALF*(dWdyl + dWdyr);
							}
						}
// #ifdef _EB_PARALLEL_DEBUG_
// 		dout << endl << " ==========================================================";
// 		dout << endl << " NORTH face of cell (" << i << "," << j << "," << Local_Solution_Block_List->Block[nb].gblknum << ") VISCOUS flux:";
// 		dout << endl << " Diamond-path = " << viscous_bc_flag;
// 		dout << endl << " xfaceN =" << Local_SolnBlk[nb].Grid.xfaceN(i,j);
// 		dout << endl << " nfaceN =" << Local_SolnBlk[nb].Grid.nfaceN(i,j);
// 		dout << endl << " lfaceN = " << Local_SolnBlk[nb].Grid.lfaceN(i,j);
// 		dout << endl << " Xl =" << Xl;
// 		dout << endl << " Xr =" << Xr;
// 		dout << endl << " Xu =" << Xu;
// 		dout << endl << " Xd =" << Xd;
// 		dout << endl << " Wl =" << Wl;
// 		dout << endl << " Wr =" << Wr;
// 		dout << endl << " Wu =" << Wu;
// 		dout << endl << " Wd =" << Wd;
// 		dout << endl << " F =" << ViscousFluxDiamondPath_n(Local_SolnBlk[nb].Grid.xfaceN(i,j),
// 								   Xl,Wl,Xd,Wd,Xr,Wr,Xu,Wu,
// 								   Local_SolnBlk[nb].Grid.nfaceN(i,j),
// 								   Local_SolnBlk[nb].Axisymmetric,
// 								   viscous_bc_flag);
// #endif
					}
				}
				// Compute the NORTH face viscous flux.
				Flux -= ViscousFlux_n(Xface, Wface, dWdx, dWdy,
						norm_dir, Local_SolnBlk[nb].Axisymmetric, adiabatic_flag);
			} while (0); // end of viscous flux clause.

	  // Evaluate cell-averaged solution changes.
	  Local_SolnBlk[nb].dUdt[i][j  ][k_residual] -= (IP->CFL_Number*Local_SolnBlk[nb].dt[i][j  ])*
                                                        Flux*Local_SolnBlk[nb].Grid.lfaceN(i,j)/Local_SolnBlk[nb].Grid.Cell[i][j].A;
	  Local_SolnBlk[nb].dUdt[i][j+1][k_residual] += (IP->CFL_Number*Local_SolnBlk[nb].dt[i][j+1])*
                                                        Flux*Local_SolnBlk[nb].Grid.lfaceS(i,j+1)/Local_SolnBlk[nb].Grid.Cell[i][j+1].A;

	  // Save south and north face boundary flux.
	  if (j == Local_SolnBlk[nb].JCl-1) {
	    Local_SolnBlk[nb].FluxS[i] = -Flux*Local_SolnBlk[nb].Grid.lfaceS(i,j+1);
	  } else if (j == Local_SolnBlk[nb].JCu) {
	    Local_SolnBlk[nb].FluxN[i] = Flux*Local_SolnBlk[nb].Grid.lfaceN(i,j);
	  }

	}

	Local_SolnBlk[nb].dUdt[i][Local_SolnBlk[nb].JCl-1][k_residual].Vacuum();
	Local_SolnBlk[nb].dUdt[i][Local_SolnBlk[nb].JCu+1][k_residual].Vacuum();

      }

    }
  }

  // Residuals for each quadrilateral multi-block solution block
  // successfully calculated.
  return 0;

}

#endif // _EMBEDDEDBOUNDARIES2D_NAVIERSTOKES2D_INCLUDED

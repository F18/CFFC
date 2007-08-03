/**********************************************************************
 * EmbeddedBoundaries2D: Header file declaring 2D embedded boundary   *
 *                       classes and functions.                       *
 **********************************************************************/

#ifndef _EMBEDDEDBOUNDARIES2D_GAUSSIAN2D_INCLUDED
#define _EMBEDDEDBOUNDARIES2D_GAUSSIAN2D_INCLUDED

// Include 2D Gaussian quadrilateral mesh solution header file.

#ifndef _GAUSSIAN2D_QUAD_INCLUDED
#include "../Gaussian2D/Gaussian2DQuad.h"
#endif // _GAUSSIAN2D_QUAD_INCLUDED

// Include Embedded Boundaries input header file.

#ifndef _EMBEDDEDBOUNDARIES2D_INCLUDED
#include "EmbeddedBoundaries2D.h"
#endif // _EMBEDDEDBOUNDARIES2D_INCLUDED

/**********************************************************************
 * EmbeddedBoundaries2D::Calculate_Refinement_Criteria --             *
 **********************************************************************/
template <> inline void EmbeddedBoundaries2D<Gaussian2D_cState,
					     Gaussian2D_pState,
					     Gaussian2D_Quad_Block,
					     Gaussian2D_Input_Parameters>::
Calculate_Refinement_Criteria(const int &nb, const int &i, const int &j,
			      double *refinement_criteria) {

  double grad_rho_x, grad_rho_y, grad_rho_abs, grad_rho_criteria,
         div_V, div_V_criteria, curl_V_z, curl_V_abs, curl_V_criteria;
  int refinement_criteria_number;

  double a = 0.5*(Local_SolnBlk[nb].W[i][j].axx()+Local_SolnBlk[nb].W[i][j].ayy());

  // Refinement criteria that can be used:
  // (1) Refinement criteria based on the gradient of the density field;
  // (2) Refinement criteria based on the divergence of the velocity vector;
  // (3) Refinement criteria based on the curl of the velocity vector;

  // Reconstruct the solution within the cell.
  Linear_Least_Squares_Reconstruction(nb,i,j);

  // Evaluate refinement criteria #1 based on the gradient of the 
  // density field.
  if (IP->Refinement_Criteria_Gradient_Density) {
    grad_rho_x = Local_SolnBlk[nb].dWdx[i][j].d;
    grad_rho_y = Local_SolnBlk[nb].dWdy[i][j].d;
    grad_rho_abs = sqrt(sqr(grad_rho_x) + sqr(grad_rho_y));
    //grad_rho_criteria = sqrt(Local_SolnBlk[nb].Grid.Cell[i][j].A)*grad_rho_abs/Local_SolnBlk[nb].W[i][j].d;
    grad_rho_criteria = grad_rho_abs/Local_SolnBlk[nb].W[i][j].d;
  } else {
    grad_rho_criteria = ZERO;
  }

  // Evaluate refinement criteria #2 based on the divergence of the
  // velocity vector.
  if (IP->Refinement_Criteria_Divergence_Velocity) {
    div_V = Local_SolnBlk[nb].dWdx[i][j].v.x + Local_SolnBlk[nb].dWdy[i][j].v.y;
    //div_V_criteria = sqrt(Local_SolnBlk[nb].Grid.Cell[i][j].A)*fabs(div_V)/a;
    div_V_criteria = fabs(div_V)/a;
  } else {
    div_V_criteria = ZERO;
  }

  // Evaluate refinement criteria #3 based on the curl of the
  // velocity vector.
  if (IP->Refinement_Criteria_Curl_Velocity) {
    curl_V_z = Local_SolnBlk[nb].dWdx[i][j].v.y - Local_SolnBlk[nb].dWdy[i][j].v.x; 
    curl_V_abs = sqrt(sqr(curl_V_z)); 
    //curl_V_criteria = sqrt(Local_SolnBlk[nb].Grid.Cell[i][j].A)*curl_V_abs/a;
    curl_V_criteria = curl_V_abs/a;
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
template <> inline void EmbeddedBoundaries2D<Gaussian2D_cState,
					     Gaussian2D_pState,
					     Gaussian2D_Quad_Block,
					     Gaussian2D_Input_Parameters>::
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

  Gaussian2D_pState W, We;
  Vector2D X, dX;
  double eta, f, fp, fpp, Rex, linf, Cf, Cfe;

  Linear_Least_Squares_Reconstruction(nb);

  // Output node solution data.  
  Out_File_Soln << setprecision(14);
  if (Output_Title_Soln) {
    Out_File_Soln << "TITLE = \"" << CFDkit_Name() << ": 2D Flat Plate Gaussian Solution."
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
		  << "\"pxy_exact\" \\ \n"
		  << "\"eta\" \\ \n"
		  << "\"f\" \\ \n"
		  << "\"fp\" \\ \n"
		  << "\"fpp\" \\ \n"
		  << "\"u/u0\" \\ \n"
		  << "\"u_exact/u0\" \\ \n";
  } 
  Out_File_Soln << "ZONE T =  \"Block Number = " << Local_Solution_Block_List->Block[nb].gblknum
		<< "\" \\ \n"
		<< "I = " << Local_SolnBlk[nb].ICu - Local_SolnBlk[nb].ICl + 1 << " \\ \n"
		<< "J = " << Local_SolnBlk[nb].JCu - Local_SolnBlk[nb].JCl + 1 << " \\ \n"
		<< "F = POINT \n";
  

  for (int j = Local_SolnBlk[nb].JCl; j <= Local_SolnBlk[nb].JCu; j++) {
    for (int i = Local_SolnBlk[nb].ICl; i <= Local_SolnBlk[nb].ICu; i++) {
      // Get cell position and solution data.
      X = Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
      W = Local_SolnBlk[nb].W[i][j];
      // Get exact solution.
      if (Mesh[nb].cell_status[i][j] == CELL_STATUS_ACTIVE) {
      	We = FlatPlate(IP->Wo,X,eta,f,fp,fpp);
      } else {
	We.Standard_Atmosphere(); eta = ZERO; f = ZERO; fp = ZERO; fpp = ZERO;
      }

      // Output data.
      Out_File_Soln << " " << X
		    << W
		    << We.v;
      Out_File_Soln.setf(ios::scientific);
      Out_File_Soln << " " << We.p.xy
		    << " " << eta
		    << " " << f
		    << " " << fp
		    << " " << fpp
		    << " " << W.v.x/IP->Wo.v.x
		    << " " << We.v.x/IP->Wo.v.x
		    << endl;
      Out_File_Soln.unsetf(ios::scientific);
    } /* endfor */
  } /* endfor */
  Out_File_Soln << setprecision(6);

  // Output Skin frinction coefficient data

  if (Output_Title_Skin) {
    Out_File_Skin << "TITLE = \"" << CFDkit_Name() << ": Flat Plate Skin Friction Coefficient, "
		  << "\"" << "\n"
		  << "VARIABLES = \"x\" \\ \n"
		  << "\"y\" \\ \n"
		  << "\"Rex\" \\ \n"
		  << "\"Cf\" \\ \n"
		  << "\"Cf_e\" \\ \n"
		  << "\"Cf-Cf_e\" \\ \n"
		  << "\"pxx\" \\ \n"
		  << "\"pyy\" \\ \n"
		  << "\"total_p\" \\ \n"
		  << "\"pxy\" \\ \n"
		  << "\"Area\" \\ \n";
  }

  if (Interface_Component_List.Ni) {

    //count number of points that will be outputted (inefficient but works)
    for (int j = Local_SolnBlk[nb].JCl; j <= Local_SolnBlk[nb].JCu; j++) {
      for (int i = Local_SolnBlk[nb].ICl; i <= Local_SolnBlk[nb].ICu; i++) {
	if (Mesh[nb].cell_status[i][j] == CELL_STATUS_ACTIVE &&
	    (Mesh[nb].cell_status[i-1][j] != CELL_STATUS_ACTIVE ||
	     Mesh[nb].cell_status[i+1][j] != CELL_STATUS_ACTIVE ||
	     Mesh[nb].cell_status[i][j-1] != CELL_STATUS_ACTIVE ||
	     Mesh[nb].cell_status[i][j+1] != CELL_STATUS_ACTIVE) &&
	    Local_SolnBlk[nb].Grid.Cell[i][j].Xc.x >= TOLER &&
	    Local_SolnBlk[nb].Grid.Cell[i][j].Xc.x <= IP->Plate_Length) {
	  numberofcells_cf++;
	}
      }
    }

    //output 
    if(numberofcells_cf > 0) {
//      Out_File_Skin << "ZONE T =  \"Block Number = " << nb
//        	    << "\" \\ \n"
//    		    << "I = " << numberofcells_cf << " \\ \n"
//    		    << "F = POINT \n";
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
	      dX = Local_SolnBlk[nb].Grid.xfaceW(i,j)-Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	      W = Local_SolnBlk[nb].W[i][j] + (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdx[i  ][j])*dX.x +
		(Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdy[i][j])*dX.y;
	      Cf = -TWO*W.p.xy/(IP->Wo.d*sqr(IP->Wo.v.x));
	      Out_File_Skin.setf(ios::scientific);
	      Out_File_Skin << " " << Local_SolnBlk[nb].Grid.xfaceW(i,j).x
			    << " " << Local_SolnBlk[nb].Grid.xfaceW(i,j).y
			    << " " << Rex
			    << " " << Cf
			    << " " << Cfe
			    << " " << Cf - Cfe
			    << " " << W.p.xx
			    << " " << W.p.yy
			    << " " << W.pressure()
			    << " " << W.p.xy
			    << " " << Local_SolnBlk[nb].Grid.Cell[i][j].A
			    << endl;
	    } else if (Mesh[nb].cell_status[i+1][j] != CELL_STATUS_ACTIVE &&
		       Local_SolnBlk[nb].Grid.lfaceE(i,j) > TOLER) {
	      Rex = (IP->Wo.v.x/IP->Wo.nu())*Local_SolnBlk[nb].Grid.xfaceE(i,j).x;
	      Cfe = -TWO*0.33206/max(sqrt(Rex),NANO);
	      dX = Local_SolnBlk[nb].Grid.xfaceE(i,j)-Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	      W = Local_SolnBlk[nb].W[i][j] + (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdx[i  ][j])*dX.x +
		(Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdy[i][j])*dX.y;
	      Cf = TWO*W.p.xy/(IP->Wo.d*sqr(IP->Wo.v.x));
	      Out_File_Skin.setf(ios::scientific);
	      Out_File_Skin << " " << Local_SolnBlk[nb].Grid.xfaceE(i,j).x
			    << " " << Local_SolnBlk[nb].Grid.xfaceE(i,j).y
			    << " " << Rex
			    << " " << Cf
			    << " " << Cfe
			    << " " << Cf - Cfe
			    << " " << W.p.xx
			    << " " << W.p.yy
			    << " " << W.pressure()
			    << " " << W.p.xy
			    << " " << Local_SolnBlk[nb].Grid.Cell[i][j].A
			    << endl;
	    } else if (Mesh[nb].cell_status[i][j-1] != CELL_STATUS_ACTIVE &&
		       Local_SolnBlk[nb].Grid.lfaceS(i,j) > TOLER) {
	      Rex = (IP->Wo.v.x/IP->Wo.nu())*Local_SolnBlk[nb].Grid.xfaceS(i,j).x;
	      Cfe = TWO*0.33206/max(sqrt(Rex),NANO);
	      dX = Local_SolnBlk[nb].Grid.xfaceS(i,j)-Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	      W = Local_SolnBlk[nb].W[i][j] + (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdx[i  ][j])*dX.x +
		(Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdy[i][j])*dX.y;
	      Cf = -TWO*W.p.xy/(IP->Wo.d*sqr(IP->Wo.v.x));
	      Out_File_Skin.setf(ios::scientific);
	      Out_File_Skin << " " << Local_SolnBlk[nb].Grid.xfaceS(i,j).x
			    << " " << Local_SolnBlk[nb].Grid.xfaceS(i,j).y
			    << " " << Rex
			    << " " << Cf
			    << " " << Cfe
			    << " " << Cf - Cfe
			    << " " << W.p.xx
			    << " " << W.p.yy
			    << " " << W.pressure()
			    << " " << W.p.xy
			    << " " << Local_SolnBlk[nb].Grid.Cell[i][j].A
			    << endl;
	    } else if (Mesh[nb].cell_status[i][j+1] != CELL_STATUS_ACTIVE &&
		       Local_SolnBlk[nb].Grid.lfaceN(i,j) > TOLER) {
	      Rex = (IP->Wo.v.x/IP->Wo.nu())*Local_SolnBlk[nb].Grid.xfaceN(i,j).x;
	      Cfe = TWO*0.33206/max(sqrt(Rex),NANO);
	      dX = Local_SolnBlk[nb].Grid.xfaceN(i,j)-Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	      W = Local_SolnBlk[nb].W[i][j] + (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdx[i  ][j])*dX.x +
		(Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdy[i][j])*dX.y;
	      Cf = -TWO*W.p.xy/(IP->Wo.d*sqr(IP->Wo.v.x));
	      Out_File_Skin.setf(ios::scientific);
	      Out_File_Skin << " " << Local_SolnBlk[nb].Grid.xfaceN(i,j).x
			    << " " << Local_SolnBlk[nb].Grid.xfaceN(i,j).y
			    << " " << Rex
			    << " " << Cf
			    << " " << Cfe
			    << " " << Cf - Cfe
			    << " " << W.p.xx
			    << " " << W.p.yy
			    << " " << W.pressure()
			    << " " << W.p.xy
			    << " " << Local_SolnBlk[nb].Grid.Cell[i][j].A
			    << endl;
	    }
	  }
	}
      }
    }
  }

  return;
}

/**********************************************************************
 * EmbeddedBoundaries2D::Net_Force --                                 *
 **********************************************************************/
template <> inline Vector2D EmbeddedBoundaries2D<Gaussian2D_cState,
					     Gaussian2D_pState,
					     Gaussian2D_Quad_Block,
					     Gaussian2D_Input_Parameters>::
Net_Force(void) {

  double fx(0.0), fy(0.0);
  double average_edge, sin_angle, cos_angle;
  double normal_pressure, shear_pressure;
  Vector2D dX;
  Gaussian2D_pState W;

  for(int nb=0; nb < Local_Solution_Block_List->Nblk; nb++) {
    if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      Linear_Least_Squares_Reconstruction(nb);

      for (int j = Local_SolnBlk[nb].JCl; j <= Local_SolnBlk[nb].JCu; j++) {
	for (int i = Local_SolnBlk[nb].ICl; i <= Local_SolnBlk[nb].ICu; i++) {

	  if (Mesh[nb].cell_status[i][j] == CELL_STATUS_ACTIVE) {
	    average_edge = 0.25*(Local_SolnBlk[nb].Grid.lfaceS(i,j) +
				 Local_SolnBlk[nb].Grid.lfaceN(i,j) +
				 Local_SolnBlk[nb].Grid.lfaceW(i,j) +
				 Local_SolnBlk[nb].Grid.lfaceE(i,j));

	    //south edge
	    if(Mesh[nb].cell_status[i][j-1] != CELL_STATUS_ACTIVE &&
	       Local_SolnBlk[nb].Grid.lfaceS(i,j) > TOLER*average_edge) {
	      dX = Local_SolnBlk[nb].Grid.xfaceS(i,j)-Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	      W = Local_SolnBlk[nb].W[i][j] + (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdx[i][j])*dX.x +
		(Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdy[i][j])*dX.y;
	      cos_angle = Local_SolnBlk[nb].Grid.nfaceS(i,j).x; 
	      sin_angle = Local_SolnBlk[nb].Grid.nfaceS(i,j).y;
	      normal_pressure  = W.p.xx*cos_angle*cos_angle+W.p.yy*sin_angle*sin_angle
		+2.0*W.p.xy*cos_angle*sin_angle;
	      shear_pressure  = W.p.xy*(cos_angle*cos_angle-sin_angle*sin_angle)
		-(W.p.xx-W.p.yy)*cos_angle*sin_angle;
	      fx += (normal_pressure*cos_angle-shear_pressure*sin_angle)*
		Local_SolnBlk[nb].Grid.lfaceS(i,j);
	      fy += (normal_pressure*sin_angle+shear_pressure*cos_angle)*
		Local_SolnBlk[nb].Grid.lfaceS(i,j);
	    }
	    //north edge
	    if(Mesh[nb].cell_status[i][j+1] != CELL_STATUS_ACTIVE &&
	       Local_SolnBlk[nb].Grid.lfaceN(i,j) > TOLER*average_edge) {
	      dX = Local_SolnBlk[nb].Grid.xfaceN(i,j)-Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	      W = Local_SolnBlk[nb].W[i][j] + (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdx[i][j])*dX.x +
		(Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdy[i][j])*dX.y;
	      cos_angle = Local_SolnBlk[nb].Grid.nfaceN(i,j).x; 
	      sin_angle = Local_SolnBlk[nb].Grid.nfaceN(i,j).y;
	      normal_pressure  = W.p.xx*cos_angle*cos_angle+W.p.yy*sin_angle*sin_angle
		+2.0*W.p.xy*cos_angle*sin_angle;
	      shear_pressure  = W.p.xy*(cos_angle*cos_angle-sin_angle*sin_angle)
		-(W.p.xx-W.p.yy)*cos_angle*sin_angle;
	      fx += (normal_pressure*cos_angle-shear_pressure*sin_angle)*
		Local_SolnBlk[nb].Grid.lfaceN(i,j);
	      fy += (normal_pressure*sin_angle+shear_pressure*cos_angle)*
		Local_SolnBlk[nb].Grid.lfaceN(i,j);
	    }
	    //west edge
	    if(Mesh[nb].cell_status[i-1][j] != CELL_STATUS_ACTIVE &&
	       Local_SolnBlk[nb].Grid.lfaceW(i,j) > TOLER*average_edge) {
	      dX = Local_SolnBlk[nb].Grid.xfaceW(i,j)-Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	      W = Local_SolnBlk[nb].W[i][j] + (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdx[i][j])*dX.x +
	        (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdy[i][j])*dX.y;
	      cos_angle = Local_SolnBlk[nb].Grid.nfaceW(i,j).x; 
	      sin_angle = Local_SolnBlk[nb].Grid.nfaceW(i,j).y;
	      normal_pressure  = W.p.xx*cos_angle*cos_angle+W.p.yy*sin_angle*sin_angle
		+2.0*W.p.xy*cos_angle*sin_angle;
	      shear_pressure  = W.p.xy*(cos_angle*cos_angle-sin_angle*sin_angle)
		-(W.p.xx-W.p.yy)*cos_angle*sin_angle;
	      fx += (normal_pressure*cos_angle-shear_pressure*sin_angle)*
		Local_SolnBlk[nb].Grid.lfaceW(i,j);
	      fy += (normal_pressure*sin_angle+shear_pressure*cos_angle)*
		Local_SolnBlk[nb].Grid.lfaceW(i,j);
	    }
	    //east edge
	    if(Mesh[nb].cell_status[i+1][j] != CELL_STATUS_ACTIVE &&
	       Local_SolnBlk[nb].Grid.lfaceE(i,j) > TOLER*average_edge) {
	      dX = Local_SolnBlk[nb].Grid.xfaceE(i,j)-Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	      W = Local_SolnBlk[nb].W[i][j] + (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdx[i][j])*dX.x +
	        (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdy[i][j])*dX.y;
	      cos_angle = Local_SolnBlk[nb].Grid.nfaceE(i,j).x; 
	      sin_angle = Local_SolnBlk[nb].Grid.nfaceE(i,j).y;
	      normal_pressure  = W.p.xx*cos_angle*cos_angle+W.p.yy*sin_angle*sin_angle
		+2.0*W.p.xy*cos_angle*sin_angle;
	      shear_pressure  = W.p.xy*(cos_angle*cos_angle-sin_angle*sin_angle)
		-(W.p.xx-W.p.yy)*cos_angle*sin_angle;
	      fx += (normal_pressure*cos_angle-shear_pressure*sin_angle)*
		Local_SolnBlk[nb].Grid.lfaceE(i,j);
	      fy += (normal_pressure*sin_angle+shear_pressure*cos_angle)*
		Local_SolnBlk[nb].Grid.lfaceE(i,j);
	    }

	  }  //if cell is active

	}  //end i-direction
      }  //end j-direction
    } // if block is used
  } //end number of blocks


#ifdef _MPI_VERSION
  fx = CFDkit_Summation_MPI(fx);
  fy = CFDkit_Summation_MPI(fy);
#endif

  return Vector2D(fx,fy);
}

/**********************************************************************
 * EmbeddedBoundaries2D::Net_Pressure_Force --                        *
 **********************************************************************/
template <> inline Vector2D EmbeddedBoundaries2D<Gaussian2D_cState,
					     Gaussian2D_pState,
					     Gaussian2D_Quad_Block,
					     Gaussian2D_Input_Parameters>::
Net_Pressure_Force(void) {

  double fx(0.0), fy(0.0);
  double average_edge, sin_angle, cos_angle;
  double normal_pressure;
  Vector2D dX;
  Gaussian2D_pState W;

  for(int nb=0; nb < Local_Solution_Block_List->Nblk; nb++) {
    if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      Linear_Least_Squares_Reconstruction(nb);

      for (int j = Local_SolnBlk[nb].JCl; j <= Local_SolnBlk[nb].JCu; j++) {
	for (int i = Local_SolnBlk[nb].ICl; i <= Local_SolnBlk[nb].ICu; i++) {

	  if (Mesh[nb].cell_status[i][j] == CELL_STATUS_ACTIVE) {
	    average_edge = 0.25*(Local_SolnBlk[nb].Grid.lfaceS(i,j) +
				 Local_SolnBlk[nb].Grid.lfaceN(i,j) +
				 Local_SolnBlk[nb].Grid.lfaceW(i,j) +
				 Local_SolnBlk[nb].Grid.lfaceE(i,j));

	    //south edge
	    if(Mesh[nb].cell_status[i][j-1] != CELL_STATUS_ACTIVE &&
	       Local_SolnBlk[nb].Grid.lfaceS(i,j) > TOLER*average_edge) {
	      dX = Local_SolnBlk[nb].Grid.xfaceS(i,j)-Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	      W = Local_SolnBlk[nb].W[i][j] + (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdx[i][j])*dX.x +
		(Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdy[i][j])*dX.y;
	      cos_angle = Local_SolnBlk[nb].Grid.nfaceS(i,j).x; 
	      sin_angle = Local_SolnBlk[nb].Grid.nfaceS(i,j).y;
	      normal_pressure  = W.p.xx*cos_angle*cos_angle+W.p.yy*sin_angle*sin_angle
		+2.0*W.p.xy*cos_angle*sin_angle;
	      fx += (normal_pressure*cos_angle)*
		Local_SolnBlk[nb].Grid.lfaceS(i,j);
	      fy += (normal_pressure*sin_angle)*
		Local_SolnBlk[nb].Grid.lfaceS(i,j);
	    }
	    //north edge
	    if(Mesh[nb].cell_status[i][j+1] != CELL_STATUS_ACTIVE &&
	       Local_SolnBlk[nb].Grid.lfaceN(i,j) > TOLER*average_edge) {
	      dX = Local_SolnBlk[nb].Grid.xfaceN(i,j)-Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	      W = Local_SolnBlk[nb].W[i][j] + (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdx[i][j])*dX.x +
		(Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdy[i][j])*dX.y;
	      cos_angle = Local_SolnBlk[nb].Grid.nfaceN(i,j).x; 
	      sin_angle = Local_SolnBlk[nb].Grid.nfaceN(i,j).y;
	      normal_pressure  = W.p.xx*cos_angle*cos_angle+W.p.yy*sin_angle*sin_angle
		+2.0*W.p.xy*cos_angle*sin_angle;
	      fx += (normal_pressure*cos_angle)*
		Local_SolnBlk[nb].Grid.lfaceN(i,j);
	      fy += (normal_pressure*sin_angle)*
		Local_SolnBlk[nb].Grid.lfaceN(i,j);
	    }
	    //west edge
	    if(Mesh[nb].cell_status[i-1][j] != CELL_STATUS_ACTIVE &&
	       Local_SolnBlk[nb].Grid.lfaceW(i,j) > TOLER*average_edge) {
	      dX = Local_SolnBlk[nb].Grid.xfaceW(i,j)-Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	      W = Local_SolnBlk[nb].W[i][j] + (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdx[i][j])*dX.x +
	        (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdy[i][j])*dX.y;
	      cos_angle = Local_SolnBlk[nb].Grid.nfaceW(i,j).x; 
	      sin_angle = Local_SolnBlk[nb].Grid.nfaceW(i,j).y;
	      normal_pressure  = W.p.xx*cos_angle*cos_angle+W.p.yy*sin_angle*sin_angle
		+2.0*W.p.xy*cos_angle*sin_angle;
	      fx += (normal_pressure*cos_angle)*
		Local_SolnBlk[nb].Grid.lfaceW(i,j);
	      fy += (normal_pressure*sin_angle)*
		Local_SolnBlk[nb].Grid.lfaceW(i,j);
	    }
	    //east edge
	    if(Mesh[nb].cell_status[i+1][j] != CELL_STATUS_ACTIVE &&
	       Local_SolnBlk[nb].Grid.lfaceE(i,j) > TOLER*average_edge) {
	      dX = Local_SolnBlk[nb].Grid.xfaceE(i,j)-Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	      W = Local_SolnBlk[nb].W[i][j] + (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdx[i][j])*dX.x +
	        (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdy[i][j])*dX.y;
	      cos_angle = Local_SolnBlk[nb].Grid.nfaceE(i,j).x; 
	      sin_angle = Local_SolnBlk[nb].Grid.nfaceE(i,j).y;
	      normal_pressure  = W.p.xx*cos_angle*cos_angle+W.p.yy*sin_angle*sin_angle
		+2.0*W.p.xy*cos_angle*sin_angle;
	      fx += (normal_pressure*cos_angle)*
		Local_SolnBlk[nb].Grid.lfaceE(i,j);
	      fy += (normal_pressure*sin_angle)*
		Local_SolnBlk[nb].Grid.lfaceE(i,j);
	    }

	  }  //if cell is active

	}  //end i-direction
      }  //end j-direction
    } // if block is used
  } //end number of blocks


#ifdef _MPI_VERSION
  fx = CFDkit_Summation_MPI(fx);
  fy = CFDkit_Summation_MPI(fy);
#endif

  return Vector2D(fx,fy);
}

/**********************************************************************
 * EmbeddedBoundaries2D::Output_Cylinder_Drag --                      *
 **********************************************************************/
template <> inline int EmbeddedBoundaries2D<Gaussian2D_cState,
					     Gaussian2D_pState,
					     Gaussian2D_Quad_Block,
					     Gaussian2D_Input_Parameters>::
Output_Couette(void) {

  double u_top_wall_ave(0.0), u_bottom_wall_ave(0.0), pxy_ave(0.0);
  double u_top_length(0.0), u_bottom_length(0.0), pxy_area(0.0);
  Vector2D dX;
  Gaussian2D_pState W;
  double pxy_ref, average_edge;

  pxy_ref = IP->Wo.d*30.0*sqrt(2.0*BOLTZMANN*IP->Wo.T()/
			       (PI*IP->Wo.M/AVOGADRO/THOUSAND));

  for(int nb=0; nb < Local_Solution_Block_List->Nblk; nb++) {
    if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      Linear_Least_Squares_Reconstruction(nb);

      for (int j = Local_SolnBlk[nb].JCl; j <= Local_SolnBlk[nb].JCu; j++) {
	for (int i = Local_SolnBlk[nb].ICl; i <= Local_SolnBlk[nb].ICu; i++) {

	  if (Mesh[nb].cell_status[i][j] == CELL_STATUS_ACTIVE) {
	    average_edge = 0.25*(Local_SolnBlk[nb].Grid.lfaceS(i,j) +
				 Local_SolnBlk[nb].Grid.lfaceN(i,j) +
				 Local_SolnBlk[nb].Grid.lfaceW(i,j) +
				 Local_SolnBlk[nb].Grid.lfaceE(i,j));
	    pxy_ave  += Local_SolnBlk[nb].W[i][j].p.xy*Local_SolnBlk[nb].Grid.Cell[i][j].A;
	    pxy_area += Local_SolnBlk[nb].Grid.Cell[i][j].A;

	    //south edge
	    if(Mesh[nb].cell_status[i][j-1] != CELL_STATUS_ACTIVE &&
	       Local_SolnBlk[nb].Grid.lfaceS(i,j) > TOLER*average_edge) {
	      dX = Local_SolnBlk[nb].Grid.xfaceS(i,j)-Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	      W = Local_SolnBlk[nb].W[i][j] + (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdx[i][j])*dX.x +
		(Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdy[i][j])*dX.y;
	      if(Local_SolnBlk[nb].Grid.Cell[i][j].Xc.y > 0.0) {
		u_top_wall_ave += W.v.x * Local_SolnBlk[nb].Grid.lfaceS(i,j);
		u_top_length   += Local_SolnBlk[nb].Grid.lfaceS(i,j);
	      } else {
		u_bottom_wall_ave += W.v.x*Local_SolnBlk[nb].Grid.lfaceS(i,j);
		u_bottom_length   += Local_SolnBlk[nb].Grid.lfaceS(i,j);
	      }
	    }
	    //north edge
	    if(Mesh[nb].cell_status[i][j+1] != CELL_STATUS_ACTIVE &&
	       Local_SolnBlk[nb].Grid.lfaceN(i,j) > TOLER*average_edge) {
	      dX = Local_SolnBlk[nb].Grid.xfaceN(i,j)-Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	      W = Local_SolnBlk[nb].W[i][j] + (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdx[i][j])*dX.x +
		(Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdy[i][j])*dX.y;
	      if(Local_SolnBlk[nb].Grid.Cell[i][j].Xc.y > 0.0) {
		u_top_wall_ave += W.v.x * Local_SolnBlk[nb].Grid.lfaceN(i,j);
		u_top_length   += Local_SolnBlk[nb].Grid.lfaceN(i,j);
	      } else {
		u_bottom_wall_ave += W.v.x*Local_SolnBlk[nb].Grid.lfaceN(i,j);
		u_bottom_length   += Local_SolnBlk[nb].Grid.lfaceN(i,j);
	      }
	    }
	    //west edge
	    if(Mesh[nb].cell_status[i-1][j] != CELL_STATUS_ACTIVE &&
	       Local_SolnBlk[nb].Grid.lfaceW(i,j) > TOLER*average_edge) {
	      dX = Local_SolnBlk[nb].Grid.xfaceW(i,j)-Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	      W = Local_SolnBlk[nb].W[i][j] + (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdx[i][j])*dX.x +
	        (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdy[i][j])*dX.y;
	      if(Local_SolnBlk[nb].Grid.Cell[i][j].Xc.y > 0.0) {
		u_top_wall_ave += W.v.x * Local_SolnBlk[nb].Grid.lfaceW(i,j);
		u_top_length   += Local_SolnBlk[nb].Grid.lfaceW(i,j);
	      } else {
		u_bottom_wall_ave += W.v.x*Local_SolnBlk[nb].Grid.lfaceW(i,j);
		u_bottom_length   += Local_SolnBlk[nb].Grid.lfaceW(i,j);
	      }
	    }
	    //east edge
	    if(Mesh[nb].cell_status[i+1][j] != CELL_STATUS_ACTIVE &&
	       Local_SolnBlk[nb].Grid.lfaceE(i,j) > TOLER*average_edge) {
	      dX = Local_SolnBlk[nb].Grid.xfaceE(i,j)-Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	      W = Local_SolnBlk[nb].W[i][j] + (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdx[i][j])*dX.x +
	        (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdy[i][j])*dX.y;
	      if(Local_SolnBlk[nb].Grid.Cell[i][j].Xc.y > 0.0) {
		u_top_wall_ave += W.v.x * Local_SolnBlk[nb].Grid.lfaceE(i,j);
		u_top_length   += Local_SolnBlk[nb].Grid.lfaceE(i,j);
	      } else {
		u_bottom_wall_ave += W.v.x*Local_SolnBlk[nb].Grid.lfaceE(i,j);
		u_bottom_length   += Local_SolnBlk[nb].Grid.lfaceE(i,j);
	      }
	    }

	  }  //if cell is active

	}  //end i-direction
      }  //end j-direction
    } // if block is used
  } //end number of blocks


#ifdef _MPI_VERSION
  u_top_wall_ave = CFDkit_Summation_MPI(u_top_wall_ave);
  u_top_length   = CFDkit_Summation_MPI(u_top_length);
  u_bottom_wall_ave = CFDkit_Summation_MPI(u_top_wall_ave);
  u_bottom_length   = CFDkit_Summation_MPI(u_top_length);
  pxy_ave = CFDkit_Summation_MPI(pxy_ave);
  pxy_area   = CFDkit_Summation_MPI(pxy_area);
#endif
  if(CFDkit_Primary_MPI_Processor()) {

    u_top_wall_ave    = u_top_wall_ave/u_top_length;
    u_bottom_wall_ave = u_bottom_wall_ave/u_bottom_length;
    pxy_ave           = pxy_ave/pxy_area;

    cout << endl
	 << endl
	 << " ==================================================================== "
	 << endl << "                    Couette Info"
	 << endl << "   U_top    = " << u_top_wall_ave
	 << endl << "   U_bottom = " << u_bottom_wall_ave
	 << endl << "   Pxy      = " << pxy_ave << endl
	 << endl << "   U_top/U  = " << fabs(u_top_wall_ave/30.0) << "   (assumes 30 m/s)"
	 << endl << "   U_bot/U  = " << fabs(u_bottom_wall_ave/30.0) << "   (assumes 30 m/s)"
	 << endl << "   Pxy/ref  = " << fabs(pxy_ave/pxy_ref) << "   (assumes 30 m/s)"
	 << endl
	 << " ==================================================================== "
	 << endl;
  }
 
  return 0;
}

/**********************************************************************
 * EmbeddedBoundaries2D::Output_Cylinder_Drag --                      *
 **********************************************************************/
template <> inline int EmbeddedBoundaries2D<Gaussian2D_cState,
					     Gaussian2D_pState,
					     Gaussian2D_Quad_Block,
					     Gaussian2D_Input_Parameters>::
Output_Cylinder_Drag(void) {

  double lift, drag, Cl, Cd, speed, diam, S, Ma, Re, Kn;
  Vector2D force;

  force = Net_Force();

  lift = force.y;
  drag = force.x;

  if(CFDkit_Primary_MPI_Processor()) {

  diam = 2.0*Interface_Component_List[1].Length1;
  speed = sqrt(sqr(IP->Wo.v.x)+sqr(IP->Wo.v.y));
  Cd = drag/(0.5*IP->Wo.d*sqr(speed)*diam);
  Cl = lift/(0.5*IP->Wo.d*sqr(speed)*diam);

  S = speed/sqrt(2.0*AVOGADRO*BOLTZMANN*IP->Temperature/IP->Wo.M*THOUSAND);
  Ma = speed/IP->Wo.sound();
  Kn = IP->Wo.mfp()/diam;
  Re = IP->Wo.d*speed*diam/IP->Wo.viscosity();

    cout << endl
	 << endl
	 << " ==================================================================== "
	 << endl << "               Cylinder Drag and Lift Information"
	 << endl << "   Kn = " << Kn
	 << endl << "    S = " << S
	 << endl << "   Ma = " << Ma
	 << endl << "   Re = " << Re
	 << endl << endl << " drag = " << drag << " N/m "
	 << endl << " lift = " << lift << " N/m "
	 << endl << "   Cl = " << Cl
	 << endl << "   Cd = " << Cd
	 << endl
	 << " ==================================================================== "
	 << endl;
  }

  return 0;
}
/**********************************************************************
 * EmbeddedBoundaries2D::Output_Aerodynamic_Coefficients_Tecplot --   *
 **********************************************************************/
template <> inline int EmbeddedBoundaries2D<Gaussian2D_cState,
					     Gaussian2D_pState,
					     Gaussian2D_Quad_Block,
					     Gaussian2D_Input_Parameters>::
Output_Aerodynamic_Coefficients_Tecplot(const int &number_of_time_steps,
					const double &Time) {

  double lift, drag, Cl, Cd, Cn, speed, chord, S, Ma, Re, Kn, alpha;
  Vector2D force, pressure_force, nhat, ahat;

  force = Net_Force();
  pressure_force = Net_Pressure_Force();

  lift = force.y;
  drag = force.x;

  if(CFDkit_Primary_MPI_Processor()) {

  // Determine the axial and normal directions of the aerofoil:
  ahat = (Interface_Union_List[1].Spline.Xp[0] - Interface_Union_List[1].Spline.Xp[(Interface_Union_List[1].Spline.np-1)/2])/
         abs(Interface_Union_List[1].Spline.Xp[0] - Interface_Union_List[1].Spline.Xp[(Interface_Union_List[1].Spline.np-1)/2]);
  nhat = Vector2D(-ahat.y,ahat.x);

  alpha = acos(ahat.x);

  chord = Interface_Component_List[1].Length1;
  speed = sqrt(sqr(IP->Wo.v.x)+sqr(IP->Wo.v.y));
  Cd = drag/(0.5*IP->Wo.d*sqr(speed)*chord);
  Cl = lift/(0.5*IP->Wo.d*sqr(speed)*chord);
  Cn = dot(pressure_force,nhat)/(0.5*IP->Wo.d*sqr(speed)*chord);

  S = speed/sqrt(2.0*AVOGADRO*BOLTZMANN*IP->Temperature/IP->Wo.M*THOUSAND);
  Ma = speed/IP->Wo.sound();
  Kn = IP->Wo.mfp()/chord;
  Re = IP->Wo.d*speed*chord/IP->Wo.viscosity();

    cout << endl
	 << endl
	 << " ==================================================================== "
	 << endl << "               Aerodynamic Coefficient Information"
	 << endl << "                        alpha = " << alpha
	 << endl << "   Kn = " << Kn
	 << endl << "    S = " << S
	 << endl << "   Ma = " << Ma
	 << endl << "   Re = " << Re
	 << endl << endl << " drag = " << drag << " N/m "
	 << endl << " lift = " << lift << " N/m "
	 << endl << "   Cl = " << Cl
	 << endl << "   Cd = " << Cd
	 << endl << "   Cn = " << Cn
	 << endl
	 << " ==================================================================== "
	 << endl;
  }

  return 0;
}


/**********************************************************************
 * EmbeddedBoundaries2D::BCs_Interface --                             *
 **********************************************************************/
template <> inline int EmbeddedBoundaries2D<Gaussian2D_cState,
					    Gaussian2D_pState,
					    Gaussian2D_Quad_Block,
					    Gaussian2D_Input_Parameters>::
BCs_Interface(const int &nb, const double &Time) {

  // Exit immediately if none of the embedded boundaries are present in
  // the current solution block.
  if (!Adjustment_Data[nb].Interface_Present[0]) return 0;

  int error_flag, Ni, neighbour_flag;
  int Interface_BC_Type;
  double length, Temperature;
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
	      Local_SolnBlk[nb].W[i][j] += Reflect(Local_SolnBlk[nb].W[i][j+1],-Local_SolnBlk[nb].Grid.nfaceN(i,j),V)*Local_SolnBlk[nb].Grid.lfaceN(i,j);
	      break;
	    case INTERFACE_BC_ADIABATIC_WALL :
	      Local_SolnBlk[nb].W[i][j] += Adiabatic_Wall(Local_SolnBlk[nb].W[i][j+1],V,-Local_SolnBlk[nb].Grid.nfaceN(i,j))*Local_SolnBlk[nb].Grid.lfaceN(i,j);
	      break;
	    default :
	      cout << endl << " -> " << Local_SolnBlk[nb].Grid.Cell[i][j].Xc << " " << Interface_BC_Type;
	      return 777;
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
	      Local_SolnBlk[nb].W[i][j] += Reflect(Local_SolnBlk[nb].W[i][j-1],-Local_SolnBlk[nb].Grid.nfaceS(i,j),V)*Local_SolnBlk[nb].Grid.lfaceS(i,j);
	      break;
	    case INTERFACE_BC_ADIABATIC_WALL :
	      Local_SolnBlk[nb].W[i][j] += Adiabatic_Wall(Local_SolnBlk[nb].W[i][j-1],V,-Local_SolnBlk[nb].Grid.nfaceS(i,j))*Local_SolnBlk[nb].Grid.lfaceS(i,j);
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
	      Local_SolnBlk[nb].W[i][j] += Reflect(Local_SolnBlk[nb].W[i+1][j],-Local_SolnBlk[nb].Grid.nfaceE(i,j),V)*Local_SolnBlk[nb].Grid.lfaceE(i,j);
	      break;
	    case INTERFACE_BC_ADIABATIC_WALL :
	      Local_SolnBlk[nb].W[i][j] += Adiabatic_Wall(Local_SolnBlk[nb].W[i+1][j],V,-Local_SolnBlk[nb].Grid.nfaceE(i,j))*Local_SolnBlk[nb].Grid.lfaceE(i,j);
	      break;
	    default :
	      cout << endl << " -> " << Local_SolnBlk[nb].Grid.Cell[i][j].Xc << " " << Interface_BC_Type;
	      return 777;
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
	      Local_SolnBlk[nb].W[i][j] += Reflect(Local_SolnBlk[nb].W[i-1][j],-Local_SolnBlk[nb].Grid.nfaceW(i,j),V)*Local_SolnBlk[nb].Grid.lfaceW(i,j);
	      break;
	    case INTERFACE_BC_ADIABATIC_WALL :
	      Local_SolnBlk[nb].W[i][j] += Adiabatic_Wall(Local_SolnBlk[nb].W[i-1][j],V,-Local_SolnBlk[nb].Grid.nfaceW(i,j))*Local_SolnBlk[nb].Grid.lfaceW(i,j);
	      break;
	    default :
	      cout << endl << " -> " << Local_SolnBlk[nb].Grid.Cell[i][j].Xc << " " << Interface_BC_Type;
	      return 777;
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
		Local_SolnBlk[nb].W[i][j] += Reflect(Local_SolnBlk[nb].W[i][j+2],-Local_SolnBlk[nb].Grid.nfaceN(i,j+1),V)*Local_SolnBlk[nb].Grid.lfaceN(i,j+1);
		break;
	      case INTERFACE_BC_ADIABATIC_WALL :
		Local_SolnBlk[nb].W[i][j] += Adiabatic_Wall(Local_SolnBlk[nb].W[i][j+2],V,-Local_SolnBlk[nb].Grid.nfaceN(i,j+1))*Local_SolnBlk[nb].Grid.lfaceN(i,j+1);
		break;
	      default :
		cout << endl << " -> " << Local_SolnBlk[nb].Grid.Cell[i][j].Xc << " " << Interface_BC_Type;
		return 777;
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
		Local_SolnBlk[nb].W[i][j] += Reflect(Local_SolnBlk[nb].W[i][j-2],-Local_SolnBlk[nb].Grid.nfaceS(i,j-1),V)*Local_SolnBlk[nb].Grid.lfaceS(i,j-1);
		break;
	      case INTERFACE_BC_ADIABATIC_WALL :
		Local_SolnBlk[nb].W[i][j] += Adiabatic_Wall(Local_SolnBlk[nb].W[i][j-2],V,-Local_SolnBlk[nb].Grid.nfaceS(i,j-1))*Local_SolnBlk[nb].Grid.lfaceS(i,j-1);
		break;
	      default :
		cout << endl << " -> " << Local_SolnBlk[nb].Grid.Cell[i][j].Xc << " " << Interface_BC_Type;
		return 777;
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
		Local_SolnBlk[nb].W[i][j] += Reflect(Local_SolnBlk[nb].W[i+2][j],-Local_SolnBlk[nb].Grid.nfaceE(i+1,j),V)*Local_SolnBlk[nb].Grid.lfaceE(i+1,j);
		break;
	      case INTERFACE_BC_ADIABATIC_WALL :
		Local_SolnBlk[nb].W[i][j] += Adiabatic_Wall(Local_SolnBlk[nb].W[i+2][j],V,-Local_SolnBlk[nb].Grid.nfaceE(i+1,j))*Local_SolnBlk[nb].Grid.lfaceE(i+1,j);
		break;
	      default :
		cout << endl << " -> " << Local_SolnBlk[nb].Grid.Cell[i][j].Xc << " " << Interface_BC_Type;
		return 777;
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
		Local_SolnBlk[nb].W[i][j] += Reflect(Local_SolnBlk[nb].W[i-2][j],-Local_SolnBlk[nb].Grid.nfaceW(i-1,j),V)*Local_SolnBlk[nb].Grid.lfaceW(i-1,j);
		break;
	      case INTERFACE_BC_ADIABATIC_WALL :
		Local_SolnBlk[nb].W[i][j] += Adiabatic_Wall(Local_SolnBlk[nb].W[i-2][j],V,-Local_SolnBlk[nb].Grid.nfaceW(i-1,j))*Local_SolnBlk[nb].Grid.lfaceW(i-1,j);
		break;
	      default :
		cout << endl << " -> " << Local_SolnBlk[nb].Grid.Cell[i][j].Xc << " " << Interface_BC_Type;
		return 777;
	      };
	    }
	  }

	}

	// Determine the final (length-averaged) boundary state.
	if (length > ZERO) Local_SolnBlk[nb].W[i][j] = Local_SolnBlk[nb].W[i][j]/length;
	else Local_SolnBlk[nb].W[i][j] = IP->Wo;

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
template <> inline double EmbeddedBoundaries2D<Gaussian2D_cState,
					       Gaussian2D_pState,
					       Gaussian2D_Quad_Block,
					       Gaussian2D_Input_Parameters>::
CFL(const double &Time) {

  double dtMin = MILLION, a, v, V, dt_inv, rel_tol;

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

	    rel_tol = TOLER*rel_tol*(Local_SolnBlk[nb].Grid.lfaceN(i,j) +
				     Local_SolnBlk[nb].Grid.lfaceS(i,j) +
				     Local_SolnBlk[nb].Grid.lfaceE(i,j) +
				     Local_SolnBlk[nb].Grid.lfaceW(i,j));

	    // Initialize the time-step.
	    dt_inv = ZERO;

	    // Determine the speed of sound in the current cell
	    // for the inviscid time-step calculation:
	    //
	    // Here we use the average of axx and ayy, this is
	    // not correct, but it is easy to implement
	    // and doesn't really matter too much, maybe I'll change it later
	    // (Remmber that the fastest wavespeed it sqrt(3.0*p/rho) or sqrt(3.0)*a_ii )
	    //                               ~james

	    a = sqrt(3.0)*(0.5*(Local_SolnBlk[nb].W[i][j].axx()+Local_SolnBlk[nb].W[i][j].ayy()));

	    // Determine the NORTH face contribution if required.
	    if (Local_SolnBlk[nb].Grid.lfaceN(i,j) > rel_tol) {
	      // Inviscid time-step.
	      if (Mesh[nb].cell_status[i][j+1] == CELL_STATUS_ACTIVE) V = ZERO;
	      else V = Interface_Union_List[Mesh[nb].cell_status[i][j+1]].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceN(i,j),Time)*Local_SolnBlk[nb].Grid.nfaceN(i,j);
	      v = Local_SolnBlk[nb].W[i][j].v*Local_SolnBlk[nb].Grid.nfaceN(i,j);
	      dt_inv += (a+fabs(v-V)+fabs(V))*Local_SolnBlk[nb].Grid.lfaceN(i,j);
	    }

	    // Determine the SOUTH face contribution if required.
	    if (Local_SolnBlk[nb].Grid.lfaceS(i,j) > rel_tol) {
	      // Inviscid time-step.
	      if (Mesh[nb].cell_status[i][j-1] == CELL_STATUS_ACTIVE) V = ZERO;
	      else V = Interface_Union_List[Mesh[nb].cell_status[i][j-1]].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceS(i,j),Time)*Local_SolnBlk[nb].Grid.nfaceS(i,j);
	      v = Local_SolnBlk[nb].W[i][j].v*Local_SolnBlk[nb].Grid.nfaceS(i,j);
	      dt_inv += (a+fabs(v-V)+fabs(V))*Local_SolnBlk[nb].Grid.lfaceS(i,j);
	    }

	    // Determine the EAST face contribution if required.
	    if (Local_SolnBlk[nb].Grid.lfaceE(i,j) > rel_tol) {
	      // Inviscid time-step.
	      if (Mesh[nb].cell_status[i+1][j] == CELL_STATUS_ACTIVE) V = ZERO;
	      else V = Interface_Union_List[Mesh[nb].cell_status[i+1][j]].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceE(i,j),Time)*Local_SolnBlk[nb].Grid.nfaceE(i,j);
	      v = Local_SolnBlk[nb].W[i][j].v*Local_SolnBlk[nb].Grid.nfaceE(i,j);
	      dt_inv += (a+fabs(v-V)+fabs(V))*Local_SolnBlk[nb].Grid.lfaceE(i,j);
	    }

	    // Determine the WEST face contribution if required.
	    if (Local_SolnBlk[nb].Grid.lfaceW(i,j) > rel_tol) {
	      // Inviscid time-step.
	      if (Mesh[nb].cell_status[i-1][j] == CELL_STATUS_ACTIVE) V = ZERO;
	      else V = Interface_Union_List[Mesh[nb].cell_status[i-1][j]].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceW(i,j),Time)*Local_SolnBlk[nb].Grid.nfaceW(i,j);
	      v = Local_SolnBlk[nb].W[i][j].v*Local_SolnBlk[nb].Grid.nfaceW(i,j);
	      dt_inv += (a+fabs(v-V)+fabs(V))*Local_SolnBlk[nb].Grid.lfaceW(i,j);
	    }

	    // Determine the minimum time-step for the cell.
	    Local_SolnBlk[nb].dt[i][j] = Local_SolnBlk[nb].Grid.Cell[i][j].A/dt_inv;

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
template <> int EmbeddedBoundaries2D<Gaussian2D_cState,
				     Gaussian2D_pState,
				     Gaussian2D_Quad_Block,
				     Gaussian2D_Input_Parameters>::
dUdt_Residual_Evaluation(const double &Time) {

  Vector2D dX;
  Gaussian2D_pState Wl, Wr;
  Gaussian2D_cState Flux;

  int Interface_BC_Type, Ni;
  Vector2D V;

  double q1, q2, k, rel_tol;

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
	  rel_tol = TOLER*0.25*(Local_SolnBlk[nb].Grid.lfaceN(i,j) +
				Local_SolnBlk[nb].Grid.lfaceS(i,j) +
				Local_SolnBlk[nb].Grid.lfaceE(i,j) +
				Local_SolnBlk[nb].Grid.lfaceW(i,j));

	  Local_SolnBlk[nb].dUdt[i+1][j][0].Vacuum();

	  if (j >= Local_SolnBlk[nb].JCl && j <= Local_SolnBlk[nb].JCu) {

 	    // Evaluate the cell interface i-direction fluxes.
 	    if (Local_SolnBlk[nb].Grid.lfaceE(i,j) < rel_tol) {

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
	      } else if (Interface_BC_Type == INTERFACE_BC_ADIABATIC_WALL) {
		// WEST face of cell (i+1,j) is an ADIABATIC_WALL boundary.
		Wl = Adiabatic_Wall(Wr,V,Local_SolnBlk[nb].Grid.nfaceW(i+1,j));
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
	      } else if (Interface_BC_Type == INTERFACE_BC_ADIABATIC_WALL) {
		// EAST face of cell (i,j) is an ADIABATIC_WALL boundary.
		Wr = Adiabatic_Wall(Wl,V,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      }

 	    } else if (Mesh[nb].cell_status[i  ][j] == CELL_STATUS_ACTIVE &&
 		       Mesh[nb].cell_status[i+1][j] == CELL_STATUS_ACTIVE) {

 	      V = Vector2D_ZERO;

 	      if (i == Local_SolnBlk[nb].ICl-1 && 
 		  (Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_REFLECTION ||
 		   Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_ADIABATIC_WALL ||
 		   //Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_COUETTE ||
 		   Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
 		   Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_CHARACTERISTIC_VELOCITY)) {

		dX = Local_SolnBlk[nb].Grid.xfaceW(i+1,j) - Local_SolnBlk[nb].Grid.Cell[i+1][j].Xc;
		Wr = Local_SolnBlk[nb].W[i+1][j] + (Local_SolnBlk[nb].phi[i+1][j]^Local_SolnBlk[nb].dWdx[i+1][j])*dX.x +
	                                           (Local_SolnBlk[nb].phi[i+1][j]^Local_SolnBlk[nb].dWdy[i+1][j])*dX.y;

		// WEST face of cell (i+1,j) is a normal boundary.
		if (Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_REFLECTION) {
		  // WEST face of cell (i+1,j) is a REFLECTION boundary.
		  Wl = Reflect(Wr,Local_SolnBlk[nb].Grid.nfaceW(i+1,j));
		} else if (Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_ADIABATIC_WALL) {
		  // WEST face of cell (i+1,j) is an ADIABATIC_WALL boundary.
		  Wl = Adiabatic_Wall(Wr,Local_SolnBlk[nb].WoW[j],Local_SolnBlk[nb].Grid.nfaceW(i+1,j));
		//} else if (Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_COUETTE) {
		//  // WEST face of cell (i+1,j) is a COUETTE boundary.
		//  Wl = BC_Couette(Wr,Local_SolnBlk[nb].WoW[j]);
		} else if (Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_CHARACTERISTIC_VELOCITY) {
		  // WEST face of cell (i+1,j) is a CHARACTERISTIC_VELOCITY boundary.
		  Wl = BC_Characteristic_Velocity(Wr,Local_SolnBlk[nb].WoW[j],Local_SolnBlk[nb].Grid.nfaceW(i+1,j));
		} else {
		  // WEST face of cell (i+1,j) is a CHARACTERISTIC boundary.
		  Wl = BC_Characteristic_Pressure(Wr,Local_SolnBlk[nb].WoW[j],Local_SolnBlk[nb].Grid.nfaceW(i+1,j));
		}

	      } else if (i == Local_SolnBlk[nb].ICu &&
			 (Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_REFLECTION ||
			  Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_ADIABATIC_WALL ||
			  //Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_COUETTE ||
			  Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_CHARACTERISTIC ||
			  Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_CHARACTERISTIC_VELOCITY)) {

		dX = Local_SolnBlk[nb].Grid.xfaceE(i,j) - Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
		Wl = Local_SolnBlk[nb].W[i][j] + (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdx[i][j])*dX.x +
	                                         (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdy[i][j])*dX.y;

		// EAST face of cell (i,j) is a normal boundary.
		if (Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_REFLECTION) {
		  // EAST face of cell (i,j) is a REFLECTION boundary.
		  Wr = Reflect(Wl,Local_SolnBlk[nb].Grid.nfaceE(i,j));
		} else if (Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_ADIABATIC_WALL) {
		  // EAST face of cell (i,j) is an ADIABATIC_WALL boundary.
		  Wr = Adiabatic_Wall(Wl,Local_SolnBlk[nb].WoE[j],Local_SolnBlk[nb].Grid.nfaceE(i,j));
		//} else if (Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_COUETTE) {
		//  // EAST face of cell (i,j) is a COUETTE boundary.
		//  Wr = BC_Couette(Wl,Local_SolnBlk[nb].WoE[j]);
		} else if (Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_CHARACTERISTIC_VELOCITY) {
		  // EAST face of cell (i,j) is a CHARATERISTIC_VELOCITY boundary.
		  Wr = BC_Characteristic_Velocity(Wl,Local_SolnBlk[nb].WoE[j],Local_SolnBlk[nb].Grid.nfaceE(i,j));
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
	    case FLUX_FUNCTION_ROE :
	      Flux = FluxRoe_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_HLLE :
	      Flux = FluxHLLE_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_ROE_MB :
	      Flux = FluxRoe_MB_n(Wl,Wr,V,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    default:
	      Flux = FluxRoe_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    };

	    // Evaluate cell-averaged solution changes.
	    Local_SolnBlk[nb].dUdt[i  ][j][0] -= Flux*Local_SolnBlk[nb].Grid.lfaceE(i,j)/
	                                         Local_SolnBlk[nb].Grid.Cell[i][j].A;
	    Local_SolnBlk[nb].dUdt[i+1][j][0] += Flux*Local_SolnBlk[nb].Grid.lfaceW(i+1,j)/
	                                         Local_SolnBlk[nb].Grid.Cell[i+1][j].A;

	    // Include all axisymmetric source term if required.
	    if (Local_SolnBlk[nb].Axisymmetric)
	      Local_SolnBlk[nb].dUdt[i][j][0] += Local_SolnBlk[nb].W[i][j].S(Local_SolnBlk[nb].Grid.Cell[i][j].Xc);

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
	  rel_tol = TOLER*0.25*(Local_SolnBlk[nb].Grid.lfaceN(i,j) +
				Local_SolnBlk[nb].Grid.lfaceS(i,j) +
				Local_SolnBlk[nb].Grid.lfaceE(i,j) +
				Local_SolnBlk[nb].Grid.lfaceW(i,j));

 	  // Evaluate the cell interface j-direction INVISCID fluxes.
 	  if (Local_SolnBlk[nb].Grid.lfaceN(i,j) < rel_tol) {
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
	    } else if (Interface_BC_Type == INTERFACE_BC_ADIABATIC_WALL) {
	      // SOUTH face of cell (i,j+1) is an ADIABATIC_WALL boundary.
	      Wl = Adiabatic_Wall(Wr,V,Local_SolnBlk[nb].Grid.nfaceS(i,j+1));
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
	    } else if (Interface_BC_Type == INTERFACE_BC_ADIABATIC_WALL) {
	      // NORTH face of cell (i,j) is an ADIABATIC_WALL boundary.
	      Wr = Adiabatic_Wall(Wl,V,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    }

 	  } else if (Mesh[nb].cell_status[i][j  ] == CELL_STATUS_ACTIVE &&
 		     Mesh[nb].cell_status[i][j+1] == CELL_STATUS_ACTIVE) {

 	    V = Vector2D_ZERO;

 	    if (j == Local_SolnBlk[nb].JCl-1 && 
 		(Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_REFLECTION ||
 		 Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_ADIABATIC_WALL ||
 		 //Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_COUETTE ||
 		 Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_CHARACTERISTIC ||
 		 Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_CHARACTERISTIC_VELOCITY)) {

	      dX = Local_SolnBlk[nb].Grid.xfaceS(i,j+1) - Local_SolnBlk[nb].Grid.Cell[i][j+1].Xc;
	      Wr = Local_SolnBlk[nb].W[i][j+1] + (Local_SolnBlk[nb].phi[i][j+1]^Local_SolnBlk[nb].dWdx[i][j+1])*dX.x +
   	                                         (Local_SolnBlk[nb].phi[i][j+1]^Local_SolnBlk[nb].dWdy[i][j+1])*dX.y;
	      // SOUTH face of cell (i,j+1) is a normal boundary.
	      if (Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_REFLECTION) {
		// SOUTH face of cell (i,j+1) is a REFLECTION boundary.
		Wl = Reflect(Wr,Local_SolnBlk[nb].Grid.nfaceS(i,j+1));
	      } else if (Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_ADIABATIC_WALL) {
		// SOUTH face of cell (i,j+1) is an ADIABATIC_WALL boundary.
		Wl = Adiabatic_Wall(Wr,Local_SolnBlk[nb].WoS[i],Local_SolnBlk[nb].Grid.nfaceS(i,j+1));
	      //else if (Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_COUETTE) {
	      //// SOUTH face of cell (i,j+1) is a COUETTE boundary.
	      //Wl = BC_Couette(Wr,Local_SolnBlk[nb].WoS[i]);
	      } else if (Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_CHARACTERISTIC_VELOCITY) {
		// SOUTH face of cell (i,j+1) is a CHARACTERISTIC_VELOCITY boundary.
		Wl = BC_Characteristic_Velocity(Wr,Local_SolnBlk[nb].WoS[i],Local_SolnBlk[nb].Grid.nfaceS(i,j+1));
	      } else if (Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_CHARACTERISTIC) {
		// SOUTH face of cell (i,j+1) is a CHARACTERISTIC boundary.
		Wl = BC_Characteristic_Pressure(Wr,Local_SolnBlk[nb].WoS[i],Local_SolnBlk[nb].Grid.nfaceS(i,j+1));
	      }

	    } else if (j == Local_SolnBlk[nb].JCu && 
		       (Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_REFLECTION ||
			Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_ADIABATIC_WALL ||
			//Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_COUETTE ||
			Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
			Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_CHARACTERISTIC_VELOCITY)) {

	      dX = Local_SolnBlk[nb].Grid.xfaceN(i,j) - Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	      Wl = Local_SolnBlk[nb].W[i][j] + (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdx[i][j])*dX.x +
   	                                       (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdy[i][j])*dX.y;
	      // NORTH face of cell (i,j) is a normal boundary.
	      if (Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_REFLECTION) {
		// NORTH face of cell (i,j) is a REFLECTION boundary.
		Wr = Reflect(Wl,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	      } else if (Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_ADIABATIC_WALL) {
		// NORTH face of cell (i,j) is an ADIABATIC_WALL boundary.
		Wr = Adiabatic_Wall(Wl,Local_SolnBlk[nb].WoN[i],Local_SolnBlk[nb].Grid.nfaceN(i,j));
	      //else if (Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_COUETTE) {
	      //// NORTH face of cell (i,j) is a COUETTE boundary.
	      //Wr = BC_Couette(Wl,Local_SolnBlk[nb].WoN[i]);
	      } else if (Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_CHARACTERISTIC_VELOCITY) {
		// NORTH face of cell (i,j) is a CHARACTERISTIC_VELOCITY boundary.
		Wr = BC_Characteristic_Velocity(Wl,Local_SolnBlk[nb].WoN[i],Local_SolnBlk[nb].Grid.nfaceN(i,j));
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
	  case FLUX_FUNCTION_ROE :
	    Flux = FluxRoe_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_HLLE :
	    Flux = FluxHLLE_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_ROE_MB :
	    Flux = FluxRoe_MB_n(Wl,Wr,V,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  default:
	    Flux = FluxRoe_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  };

	  // Evaluate cell-averaged solution changes.
	  Local_SolnBlk[nb].dUdt[i][j  ][0] -= Flux*Local_SolnBlk[nb].Grid.lfaceN(i,j)/
	                                       Local_SolnBlk[nb].Grid.Cell[i][j].A;
	  Local_SolnBlk[nb].dUdt[i][j+1][0] += Flux*Local_SolnBlk[nb].Grid.lfaceS(i,j+1)/
	                                       Local_SolnBlk[nb].Grid.Cell[i][j+1].A;

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
template <> int EmbeddedBoundaries2D<Gaussian2D_cState,
				     Gaussian2D_pState,
				     Gaussian2D_Quad_Block,
				     Gaussian2D_Input_Parameters>::
dUdt_Multistage_Explicit(const int &i_stage,
			 const double &Time) {

  int error_flag, k_residual;
  double omega, rel_tol;
  Vector2D dX;
  Gaussian2D_pState Wl, Wr;
  Gaussian2D_cState Flux;

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
	  rel_tol = TOLER*0.25*(Local_SolnBlk[nb].Grid.lfaceN(i,j) +
				Local_SolnBlk[nb].Grid.lfaceS(i,j) +
				Local_SolnBlk[nb].Grid.lfaceE(i,j) +
				Local_SolnBlk[nb].Grid.lfaceW(i,j));

	  if (i_stage == 1) {
	    Local_SolnBlk[nb].Uo[i+1][j] = Local_SolnBlk[nb].U[i+1][j];
	    Local_SolnBlk[nb].dUdt[i+1][j][k_residual].Vacuum();
	  } else if (j > Local_SolnBlk[nb].JCl-1 && j < Local_SolnBlk[nb].JCu+1) {
	    switch(IP->i_Time_Integration) {
	    case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
	      //Local_SolnBlk[nb].dUdt[i+1][j][k_residual] = Local_SolnBlk[nb].dUdt[i+1][j][k_residual];
	      break;
	    default:
	      Local_SolnBlk[nb].dUdt[i+1][j][k_residual].Vacuum();
	      break;
	    };
	  }

	  if (j >= Local_SolnBlk[nb].JCl && j <= Local_SolnBlk[nb].JCu) {

 	    // Evaluate the cell interface i-direction INVISCID fluxes.
 	    if (Local_SolnBlk[nb].Grid.lfaceE(i,j) < rel_tol) {

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
	      } else if (Interface_BC_Type == INTERFACE_BC_ADIABATIC_WALL) {
		// WEST face of cell (i+1,j) is an ADIABATIC_WALL boundary.
		Wl = Adiabatic_Wall(Wr,V,Local_SolnBlk[nb].Grid.nfaceW(i+1,j));
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
	      } else if (Interface_BC_Type == INTERFACE_BC_ADIABATIC_WALL) {
		// EAST face of cell (i,j) is a ADIABATIC_WALL boundary.
		Wr = Adiabatic_Wall(Wl,V,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      }

 	    } else if (Mesh[nb].cell_status[i  ][j] == CELL_STATUS_ACTIVE &&
 		       Mesh[nb].cell_status[i+1][j] == CELL_STATUS_ACTIVE) {

 	      V = Vector2D_ZERO;

 	      if (i == Local_SolnBlk[nb].ICl-1 && 
 		  (Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_REFLECTION ||
 		   Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_ADIABATIC_WALL ||
 		   //Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_COUETTE ||
 		   Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
 		   Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_CHARACTERISTIC_VELOCITY)) {

		dX = Local_SolnBlk[nb].Grid.xfaceW(i+1,j) - Local_SolnBlk[nb].Grid.Cell[i+1][j].Xc;
		Wr = Local_SolnBlk[nb].W[i+1][j] + (Local_SolnBlk[nb].phi[i+1][j]^Local_SolnBlk[nb].dWdx[i+1][j])*dX.x +
	                                           (Local_SolnBlk[nb].phi[i+1][j]^Local_SolnBlk[nb].dWdy[i+1][j])*dX.y;

		// WEST face of cell (i+1,j) is a normal boundary.
		if (Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_REFLECTION) {
		  // WEST face of cell (i+1,j) is a REFLECTION boundary.
		  Wl = Reflect(Wr,Local_SolnBlk[nb].Grid.nfaceW(i+1,j));
		} else if (Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_ADIABATIC_WALL) {
		  // WEST face of cell (i+1,j) is an ADIABATIC_WALL boundary.
		  Wl = Adiabatic_Wall(Wr,Local_SolnBlk[nb].WoW[j],Local_SolnBlk[nb].Grid.nfaceW(i+1,j));
		//} else if (Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_COUETTE) {
		//  // WEST face of cell (i+1,j) is a COUETTE boundary.
		//  Wl = BC_Couette(Wr,Local_SolnBlk[nb].WoW[j]);
		} else if (Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_CHARACTERISTIC_VELOCITY) {
		  // WEST face of cell (i+1,j) is a CHARACTERISTIC_VELOCITY boundary.
		  Wl = BC_Characteristic_Velocity(Wr,Local_SolnBlk[nb].WoW[j],Local_SolnBlk[nb].Grid.nfaceW(i+1,j));
		} else {
		  // WEST face of cell (i+1,j) is a CHARACTERISTIC boundary.
		  Wl = BC_Characteristic_Pressure(Wr,Local_SolnBlk[nb].WoW[j],Local_SolnBlk[nb].Grid.nfaceW(i+1,j));
		}

	      } else if (i == Local_SolnBlk[nb].ICu &&
			 (Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_REFLECTION ||
			  Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_ADIABATIC_WALL ||
			  //Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_COUETTE ||
			  Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_CHARACTERISTIC ||
			  Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_CHARACTERISTIC_VELOCITY)) {

		dX = Local_SolnBlk[nb].Grid.xfaceE(i,j) - Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
		Wl = Local_SolnBlk[nb].W[i][j] + (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdx[i][j])*dX.x +
	                                         (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdy[i][j])*dX.y;

		// EAST face of cell (i,j) is a normal boundary.
		if (Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_REFLECTION) {
		  // EAST face of cell (i,j) is a REFLECTION boundary.
		  Wr = Reflect(Wl,Local_SolnBlk[nb].Grid.nfaceE(i,j));
		} else if (Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_ADIABATIC_WALL) {
		  // EAST face of cell (i,j) is an ADIABATIC_WALL boundary.
		  Wr = Adiabatic_Wall(Wl,Local_SolnBlk[nb].WoE[j],Local_SolnBlk[nb].Grid.nfaceE(i,j));
		//} else if (Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_COUETTE) {
		//  // EAST face of cell (i,j) is a COUETTE boundary.
		//  Wr = BC_Couette(Wl,Local_SolnBlk[nb].WoE[j]);
		} else if (Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_CHARACTERISTIC_VELOCITY) {
		  // EAST face of cell (i,j) is a CHARACTERISTIC_VELOCITY boundary.
		  Wr = BC_Characteristic_Velocity(Wl,Local_SolnBlk[nb].WoE[j],Local_SolnBlk[nb].Grid.nfaceE(i,j));
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
	    case FLUX_FUNCTION_ROE :
	      Flux = FluxRoe_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_HLLE :
	      Flux = FluxHLLE_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    case FLUX_FUNCTION_ROE_MB :
	      Flux = FluxRoe_MB_n(Wl,Wr,V,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    default:
	      Flux = FluxRoe_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceE(i,j));
	      break;
	    };

	    // Evaluate cell-averaged solution changes.
	    Local_SolnBlk[nb].dUdt[i  ][j][k_residual] -= (IP->CFL_Number*Local_SolnBlk[nb].dt[i][j])*
                                                          Flux*Local_SolnBlk[nb].Grid.lfaceE(i,j)/
	                                                  Local_SolnBlk[nb].Grid.Cell[i][j].A;
	    Local_SolnBlk[nb].dUdt[i+1][j][k_residual] += (IP->CFL_Number*Local_SolnBlk[nb].dt[i+1][j])*
	                                                  Flux*Local_SolnBlk[nb].Grid.lfaceW(i+1,j)/
	                                                  Local_SolnBlk[nb].Grid.Cell[i+1][j].A;

	    // Include all axisymmetric source term if required.
	    if (Local_SolnBlk[nb].Axisymmetric)
	      Local_SolnBlk[nb].dUdt[i][j][k_residual] += (IP->CFL_Number*Local_SolnBlk[nb].dt[i][j])*
		                                          Local_SolnBlk[nb].W[i][j].S(Local_SolnBlk[nb].Grid.Cell[i][j].Xc);

 	    // Include area change source term.
 	    if (Interface_Union_List.Ni)
 	      Local_SolnBlk[nb].dUdt[i][j][k_residual] -= (IP->CFL_Number*Local_SolnBlk[nb].dt[i][j])*
 		                                          Local_SolnBlk[nb].U[i][j]*dAdt(nb,i,j,Time);

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
	  rel_tol = TOLER*0.25*(Local_SolnBlk[nb].Grid.lfaceN(i,j) +
				Local_SolnBlk[nb].Grid.lfaceS(i,j) +
				Local_SolnBlk[nb].Grid.lfaceE(i,j) +
				Local_SolnBlk[nb].Grid.lfaceW(i,j));

 	  // Evaluate the cell interface j-direction INVISCID fluxes.
 	  if (Local_SolnBlk[nb].Grid.lfaceN(i,j) < rel_tol) {
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
	    } else if (Interface_BC_Type == INTERFACE_BC_ADIABATIC_WALL) {
	      // SOUTH face of cell (i,j+1) is an ADIABATIC_WALL boundary.
	      Wl = Adiabatic_Wall(Wr,V,Local_SolnBlk[nb].Grid.nfaceS(i,j+1));
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
	    } else if (Interface_BC_Type == INTERFACE_BC_ADIABATIC_WALL) {
	      // NORTH face of cell (i,j) is an ADIABATIC_WALL boundary.
	      Wr = Adiabatic_Wall(Wl,V,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    }

 	  } else if (Mesh[nb].cell_status[i][j  ] == CELL_STATUS_ACTIVE &&
 		     Mesh[nb].cell_status[i][j+1] == CELL_STATUS_ACTIVE) {

 	    V = Vector2D_ZERO;

 	    if (j == Local_SolnBlk[nb].JCl-1 && 
 		(Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_REFLECTION ||
 		 Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_ADIABATIC_WALL ||
 		 //Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_COUETTE ||
 		 Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_CHARACTERISTIC ||
 		 Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_CHARACTERISTIC_VELOCITY)) {

	      dX = Local_SolnBlk[nb].Grid.xfaceS(i,j+1) - Local_SolnBlk[nb].Grid.Cell[i][j+1].Xc;
	      Wr = Local_SolnBlk[nb].W[i][j+1] + (Local_SolnBlk[nb].phi[i][j+1]^Local_SolnBlk[nb].dWdx[i][j+1])*dX.x +
   	                                         (Local_SolnBlk[nb].phi[i][j+1]^Local_SolnBlk[nb].dWdy[i][j+1])*dX.y;

	      // SOUTH face of cell (i,j+1) is a normal boundary.
	      if (Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_REFLECTION) {
		// SOUTH face of cell (i,j+1) is a REFLECTION boundary.
		Wl = Reflect(Wr,Local_SolnBlk[nb].Grid.nfaceS(i,j+1));
	      } else if (Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_ADIABATIC_WALL) {
		// SOUTH face of cell (i,j+1) is an ADIABATIC_WALL boundary.
		Wl = Adiabatic_Wall(Wr,Local_SolnBlk[nb].WoS[i],Local_SolnBlk[nb].Grid.nfaceS(i,j+1));
	      //else if (Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_COUETTE) {
	      //// SOUTH face of cell (i,j+1) is a COUETTE boundary.
	      //Wl = BC_Couette(Wr,Local_SolnBlk[nb].WoS[i]);
	      } else if (Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_CHARACTERISTIC_VELOCITY) {
		// SOUTH face of cell (i,j+1) is a CHARACTERISTIC_VELOCITY boundary.
		Wl = BC_Characteristic_Velocity(Wr,Local_SolnBlk[nb].WoS[i],Local_SolnBlk[nb].Grid.nfaceS(i,j+1));
	      } else if (Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_CHARACTERISTIC) {
		// SOUTH face of cell (i,j+1) is a CHARACTERISTIC boundary.
		Wl = BC_Characteristic_Pressure(Wr,Local_SolnBlk[nb].WoS[i],Local_SolnBlk[nb].Grid.nfaceS(i,j+1));
	      }

	    } else if (j == Local_SolnBlk[nb].JCu && 
		       (Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_REFLECTION ||
			Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_ADIABATIC_WALL ||
			//Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_COUETTE ||
			Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
			Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_CHARACTERISTIC_VELOCITY)) {

	      dX = Local_SolnBlk[nb].Grid.xfaceN(i,j) - Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	      Wl = Local_SolnBlk[nb].W[i][j] + (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdx[i][j])*dX.x +
   	                                       (Local_SolnBlk[nb].phi[i][j]^Local_SolnBlk[nb].dWdy[i][j])*dX.y;

	      // NORTH face of cell (i,j) is a normal boundary.
	      if (Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_REFLECTION) {
		// NORTH face of cell (i,j) is a REFLECTION boundary.
		Wr = Reflect(Wl,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	      } else if (Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_ADIABATIC_WALL) {
		// NORTH face of cell (i,j) is an ADIABATIC_WALL boundary.
		Wr = Adiabatic_Wall(Wl,Local_SolnBlk[nb].WoN[i],Local_SolnBlk[nb].Grid.nfaceN(i,j));
	      //else if (Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_COUETTE) {
	      //// NORTH face of cell (i,j) is a COUETTE boundary.
	      //Wr = BC_Couette(Wl,Local_SolnBlk[nb].WoN[i]);
	      } else if (Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_CHARACTERISTIC_VELOCITY) {
		// NORTH face of cell (i,j) is a CHARACTERISTIC_VELOCITY boundary.
		Wr = BC_Characteristic_Velocity(Wl,Local_SolnBlk[nb].WoN[i],Local_SolnBlk[nb].Grid.nfaceN(i,j));
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
	  case FLUX_FUNCTION_ROE :
	    Flux = FluxRoe_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_HLLE :
	    Flux = FluxHLLE_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  case FLUX_FUNCTION_ROE_MB :
	      Flux = FluxRoe_MB_n(Wl,Wr,V,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  default:
	    Flux = FluxRoe_n(Wl,Wr,Local_SolnBlk[nb].Grid.nfaceN(i,j));
	    break;
	  };

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

#endif // _EMBEDDEDBOUNDARIES2D_GAUSSIAN2D_INCLUDED

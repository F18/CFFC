/*! \file Levermore1D.cc
  \brief Implementation of subroutines prototyped in Levermore1D.h file. */

/* Include 1D Levermore solution header file. */
#include "Levermore1D.h"

/******************************************************//**
 * Routine: Allocate
 *
 * Allocate memory for 1D Levermore equation solution.
 *
 ********************************************************/
Levermore1D_UniformMesh* Allocate(Levermore1D_UniformMesh *Soln_ptr,
				  const CFD1D_Input_Parameters &IP) {

  int NC;                       // number of cells in the computational domain
  int Nghost; 			// number of ghost cells

  /* Calculate the total number of computational cells */
  // IP.Nghost() : calculates the number of ghost cells based on the order and the method of reconstruction.
  Nghost = IP.Nghost();
  NC = IP.Number_of_Cells + 2 * Nghost;

  /* Allocate memory. */
  Soln_ptr = new Levermore1D_UniformMesh[NC];

  /* Set preliminary mesh parameters */
  for (int i=0; i<= NC-1; ++i){
    // store domain indexes in each cell (shouldn't these all be static?)
    Soln_ptr[i].Nghost = Nghost;
    Soln_ptr[i].ICl = Nghost;
    Soln_ptr[i].ICu = NC - 1 - Nghost;
  }//endfor

  /* Return memory location. */

  return(Soln_ptr);
}

/******************************************************//**
 * Routine: Deallocate
 *
 * Deallocate memory for 1D Levermore equation solution.
 *
 ********************************************************/
Levermore1D_UniformMesh* Deallocate(Levermore1D_UniformMesh *Soln_ptr) {

  /* Deallocate memory. */
  if(Soln_ptr != NULL) {
    delete []Soln_ptr;
    Soln_ptr = NULL;
  }
  /* Return memory location. */
  return(Soln_ptr);
}

/******************************************************//**
 * Routine: Output_Gnuplot
 *
 * Writes the solution to specified output stream
 * suitable for plotting with GNUPLOT.
 *
 ********************************************************/
void Output_Gnuplot(Levermore1D_UniformMesh *Soln,
                    const int Number_of_Cells,
                    const int Number_of_Time_Steps,
                    const double &Time,
	            ostream &out_file) {

  int i;
  int ICl, ICu;

  // Set the limits of the plotted domain
  ICl = Soln[0].ICl;
  ICu = Soln[0].ICu;

  out_file << "# " << CFFC_Name() << ": 1D Levermore Solution, "
	   << "Time Step/Iteration Level = "
	   << Number_of_Time_Steps
	   << ", Time = " << Time*THOUSAND << " (ms)\n"
	   << "# cell, x(m), dx(m), rho (kg/m^3), "
	   << "v (m/s), p (Pa), T (K)\n";

  for ( i = ICl ; i <= ICu ; ++i ) {
    out_file << " " << i << " " << Soln[i] << "\n";
  } /* endfor */

  out_file << "\n";

}

/******************************************************//**
 * Routine: Output_Tecplot
 *
 * Writes the solution to specified output stream
 * suitable for plotting with TECPLOT.
 *
 ********************************************************/
void Output_Tecplot(Levermore1D_UniformMesh *Soln,
                    const CFD1D_Input_Parameters &IP,
                    const int Number_of_Time_Steps,
                    const double &Time,
	            ostream &out_file) {
//  if (IP.i_ReconstructionMethod == RECONSTRUCTION_CENO){
//    Output_Tecplot_HighOrder(Soln,
//			     IP.Number_of_Cells,
//			     Number_of_Time_Steps,
//			     Time,
//			     out_file);
//  } else {
    Output_Tecplot(Soln,
		   IP.Number_of_Cells,
		   Number_of_Time_Steps,
		   Time,
		   out_file);
//  }
}

/******************************************************//**
 * Routine: Output_Tecplot
 *
 * Writes the solution to specified output stream
 * suitable for plotting with TECPLOT.
 *
 ********************************************************/
void Output_Tecplot(Levermore1D_UniformMesh *Soln,
                    const int Number_of_Cells,
                    const int Number_of_Time_Steps,
                    const double &Time,
	            ostream &out_file) {

  int ICl, ICu;

  // Set the limits of the plotted domain
  ICl = Soln[0].ICl;
  ICu = Soln[0].ICu;

  int i;

  out_file << "TITLE = \"" << CFFC_Name() << ": 1D Levermore Solution, "
	   << Levermore1D_Vector::get_length() << " moments, "
	   << "Time Step/Iteration Level = "
	   << Number_of_Time_Steps
	   << ", Time = " << Time*THOUSAND << " (ms)\"" << "\n"
	   << "VARIABLES = \"x\" \\ \n"
	   << "\"dx\" \n";
  for(i = 0; i < Levermore1D_Vector::get_length(); ++i) {
    out_file << "\"W" << i << "\" \\ \n";
  }
  for(i = 0; i < Levermore1D_Vector::get_length(); ++i) {
    out_file << "\"U" << i << "\" \\ \n";
  }
  for(i = 0; i < Levermore1D_Vector::get_length(); ++i) {
    out_file << "\"A" << i << "\" \\ \n";
  }
  out_file << "\"detector\" \n";
  out_file << "\"N_resyncs\" \n";
  out_file << "ZONE \n";

  for ( i = ICl ; i <= ICu ; ++i ) {
    out_file << " " << Soln[i] << " " << Soln[i].U << " " << Soln[i].A << " " 
	     << Soln[i].detector << " " << Soln[i].number_of_resyncs << "\n";
  } /* endfor */

  out_file << "\n";

}

/******************************************************//**
 * Routine: Grid
 *
 * Generates a uniform mesh and assign the locations of
 * the cell centers to appropriate solution variables.
 *
 ********************************************************/
void Grid(Levermore1D_UniformMesh *Soln,
          const double &xMin,
	  const double &xMax,
	  const int Number_of_Cells) {

  int i;
  double delta_x;

  int TC;

  TC = Number_of_Cells+2*Soln[0].Nghost; // total number of cells

  /* Determine the mesh spacing. */

  delta_x = (xMax - xMin)/double(Number_of_Cells);
  Soln[0].X.setsize(delta_x);  //static function, only needs to be called once.

  /* Create the cells. */

  Soln[0].X.x = xMin - (Soln[0].Nghost - HALF)*delta_x;
//  Soln[0].CellHighOrder().AssociateGeometry(Soln[0].X);   // Associate geometry with high-order solution variables


  for ( i = 1 ; i <= TC-1 ; ++i ) {
    // Initialize the coordinate of the centroids
    Soln[i].X.x =  Soln[0].X.x + double(i)*delta_x;

//    // Associate geometry with high-order solution variables
//    Soln[i].CellHighOrder().AssociateGeometry(Soln[i].X);
  } /* endfor */
}

/******************************************************//**
 * Routine: ICs
 *
 * Assigns initial conditions and data to the
 * solution variables.
 *
 ********************************************************/
void ICs(Levermore1D_UniformMesh *Soln,
         char *gas_ptr,
	 const int i_ICtype,
         const int Number_of_Cells,
	 CFD1D_Input_Parameters &IP) {

  int i;
  Levermore1D_pState  Wl, Wr;
  Levermore1D_cState  Ul, Ur;
  Levermore1D_weights Al, Ar;

  int ICl, ICu, TC;
  double xmin, xmax;
  double a, b, dx;

  ICl = Soln[0].ICl;
  ICu = Soln[0].ICu;
  TC = Number_of_Cells+2*Soln[0].Nghost; // total number of cells

  /* Assign the gas constants for the gas of interest. */

  Soln[0].A.setgas(gas_ptr);
//  Soln[0].U.setgas(gas_ptr);
//  Soln[0].CellSolutionCharactVar().setgas(gas_ptr);

  /* Assign the initial data for the IVP of interest. */

  switch(i_ICtype) {
  case IC_SOD :
    Wl = Levermore1D_pState(DENSITY_STDATM, ZERO, PRESSURE_STDATM);
    Wr = Levermore1D_pState(DENSITY_STDATM/EIGHT,
			     ZERO,
			     PRESSURE_STDATM/TEN);
    Ul = Levermore1D_cState(Wl);
    Ur = Levermore1D_cState(Wr);
    Al = Levermore1D_weights(Ul);
    Ar = Levermore1D_weights(Ur);

    cout << endl << endl << "Ul  = " << Ul << endl
	 << "Ual = " << Levermore1D_cState(Al,Wl[2])<< endl
	 << "Al  = " << Al << endl
	 << "Wl  = " << Wl << endl << endl;
    cout << "Ur  = " << Ur << endl
	 << "Uar = " << Levermore1D_cState(Ar,Wr[2])<< endl
	 << "Ar  = " << Ar << endl
	 << "Wr  = " << Wr << endl << endl;

    for ( i = 0 ; i <= TC-1 ; ++i ) {
      if (Soln[i].X.x <= ZERO) {
	Soln[i].set_state(Wl,Ul,Al);
      } else {
	Soln[i].set_state(Wr,Ur,Ar);
      } /* end if */
    } /* endfor */
    break;
  case IC_CONSTANT :
  case IC_UNIFORM :
  default:
    Wl = Levermore1D_pState(DENSITY_STDATM, ZERO, PRESSURE_STDATM);
    Ul = Levermore1D_cState(Wl);
    Al = Levermore1D_weights(Ul);
    for ( i = 0 ; i <= TC-1 ; ++i ) {
      Soln[i].set_state(Wl, Ul, Al);
    } /* endfor */
    break;
  } /* endswitch */

  //Calculate Hessians
  for ( i = 0 ; i <= TC-1 ; ++i ) {
    Soln[i].calculate_Hessians();
  }

  //set relaxation times
  Levermore1D_pState::set_relaxation_time(IP.relaxation_time);
  Levermore1D_cState::set_relaxation_time(IP.relaxation_time);
  Levermore1D_cState::m_resync_tol = IP.resync_tol;
}

/******************************************************//**
 * Routine: CFL
 *
 * Determines the allowable global and local time steps
 * (for explicit Euler time stepping scheme) according
 * to the Courant-Friedrichs-Lewy condition.
 *
 ********************************************************/
double CFL(Levermore1D_UniformMesh *Soln,
           const int Number_of_Cells) {

  int i;
  double dtMin;
  double l_max;
  /* Determine local and global time steps. */

  dtMin = MILLION;

  for ( i = Soln[0].ICl; i <= Soln[0].ICu ; ++i ) {
    l_max = max(fabs(Soln[i].lambda_max),fabs(Soln[i].lambda_min));
    Soln[i].dt = Soln[i].X.dx/l_max;
    dtMin = min(dtMin, Soln[i].dt);
  } /* endfor */

    /* Return the global time step. */

  return (dtMin);
}

/******************************************************//**
 * Routine: Linear_Reconstruction_MUSCL
 *
 * Peforms the reconstruction of a limited piecewise
 * linear solution state within each cell of the
 * computational mesh.  The MUSCL approach of Van Leer
 * (1978) is used.  Several slope limiters may be used.
 *
 ********************************************************/
void Linear_Reconstruction_MUSCL(Levermore1D_UniformMesh *Soln,
                                 const int Number_of_Cells,
                                 const int Limiter) {

    int i, n;
    double phi;
    Levermore1D_pState a, b;

    /* Carry out the limited solution reconstruction in
       each cell. */

    b = (Soln[1].W-Soln[0].W)/Soln[0].X.dx;
    for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
        Soln[i].dWdx = (Soln[i+1].W-Soln[i-1].W)*HALF/
                       Soln[i].X.dx;
        a = b;
        b = (Soln[i+1].W-Soln[i].W)/Soln[i].X.dx;

        for ( n = 1 ; n <= Levermore1D_Vector::get_length() ; ++n ) {
           switch(Limiter) {
             case LIMITER_ZERO :
	       phi = ZERO;
               break;
             case LIMITER_ONE :
	       phi = ONE;
               break;
             case LIMITER_MINMOD :
	       phi = minmod(a[n], b[n])/(Soln[i].dWdx[n]+sgn(Soln[i].dWdx[n])*TOLER*TOLER);
               break;
             case LIMITER_UMIST :
               phi = minmod(TWO*a[n],
                            TWO*b[n],
                            (a[n]+THREE*b[n])/FOUR,
                            (THREE*a[n]+b[n])/FOUR)/
                     (Soln[i].dWdx[n]+sgn(Soln[i].dWdx[n])*TOLER*TOLER);
               break;
             case LIMITER_DOUBLE_MINMOD :
               phi = minmod(HALF*(a[n]+b[n]),
	                    TWO*a[n],
	                    TWO*b[n])/(Soln[i].dWdx[n]+sgn(Soln[i].dWdx[n])*TOLER*TOLER);
               break;
             case LIMITER_SUPERBEE :
               phi = superbee(a[n], b[n])/
                     (Soln[i].dWdx[n]+sgn(Soln[i].dWdx[n])*TOLER*TOLER);
               break;
             case LIMITER_PHI :
               phi = philimiter(a[n], b[n], 1.50)/
                     (Soln[i].dWdx[n]+sgn(Soln[i].dWdx[n])*TOLER*TOLER);
               break;
             case LIMITER_VANLEER :
               phi = vanleer(a[n], b[n])/
                     (Soln[i].dWdx[n]+sgn(Soln[i].dWdx[n])*TOLER*TOLER);
               break;
             case LIMITER_VANALBADA :
               phi = vanalbada(a[n], b[n], 0.10)/
                     (Soln[i].dWdx[n]+sgn(Soln[i].dWdx[n])*TOLER*TOLER);
               break;
	     default:
               phi = philimiter(a[n], b[n], ONE)/
                     (Soln[i].dWdx[n]+sgn(Soln[i].dWdx[n])*TOLER*TOLER);
               break;
           } /* endswitch */

           Soln[i].phi[n] = phi;
        } /* endfor */

    } /* endfor */

    Soln[0].dWdx = Soln[1].phi^Soln[1].dWdx;
    Soln[0].phi.set_all(1.0);

    Soln[Number_of_Cells+1].dWdx = Soln[Number_of_Cells].phi^Soln[Number_of_Cells].dWdx;
    Soln[Number_of_Cells+1].phi.set_all(1.0);

}

/******************************************************//**
 * Routine: Linear_Reconstruction_GreenGauss
 *
 * Peforms the reconstruction of a limited piecewise
 * linear solution state within each cell of the
 * computational mesh.  A Green-Gauss approach is used
 * in the evaluation of the unlimited solution
 * gradients.  Several slope limiters may be used.
 *
 ********************************************************/
void Linear_Reconstruction_GreenGauss(Levermore1D_UniformMesh *Soln,
                                      const int Number_of_Cells,
                                      const int Limiter) {

    int i, n;
    double u0Min, u0Max, uQuad[2], phi;

    /* Carry out the limited solution reconstruction in
       each cell. */

    for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
        Soln[i].dWdx = (Soln[i+1].W-Soln[i-1].W)*HALF/
                       Soln[i].X.dx;

        for ( n = 1 ; n <= Levermore1D_Vector::get_length() ; ++n ) {
	   u0Min = min(Soln[i-1].W[n], Soln[i].W[n]);
	   u0Min = min(u0Min, Soln[i+1].W[n]);
	   u0Max = max(Soln[i-1].W[n], Soln[i].W[n]);
	   u0Max = max(u0Max, Soln[i+1].W[n]);
	   uQuad[0] = Soln[i].W[n] - HALF*Soln[i].dWdx[n]*Soln[i].X.dx;
	   uQuad[1] = Soln[i].W[n] + HALF*Soln[i].dWdx[n]*Soln[i].X.dx;

           switch(Limiter) {
             case LIMITER_ZERO :
	       phi = ZERO;
               break;
             case LIMITER_ONE :
	       phi = ONE;
               break;
             case LIMITER_BARTH_JESPERSEN :
               phi = Limiter_BarthJespersen(uQuad, Soln[i].W[n], u0Min, u0Max, 2);
               break;
             case LIMITER_VENKATAKRISHNAN :
               phi = Limiter_Venkatakrishnan(uQuad, Soln[i].W[n], u0Min, u0Max, 2);
               break;
             case LIMITER_VANLEER :
               phi = Limiter_VanLeer(uQuad, Soln[i].W[n], u0Min, u0Max, 2);
               break;
             case LIMITER_VANALBADA :
               phi = Limiter_VanAlbada(uQuad, Soln[i].W[n], u0Min, u0Max, 2);
               break;
	     default:
               phi = Limiter_BarthJespersen(uQuad, Soln[i].W[n], u0Min, u0Max, 2);
               break;
           } /* endswitch */

           Soln[i].phi[n] = phi;
        } /* endfor */

    } /* endfor */

    Soln[0].dWdx = Soln[1].phi^Soln[1].dWdx;
    Soln[0].phi.set_all(1.0);

    Soln[Number_of_Cells+1].dWdx = Soln[Number_of_Cells].phi^Soln[Number_of_Cells].dWdx;
    Soln[Number_of_Cells+1].phi.set_all(1.0);

}


/******************************************************//**
 * Routine: Linear_Reconstruction_LeastSquares
 *
 * Peforms the reconstruction of a limited piecewise
 * linear solution state within each cell of the
 * computational mesh.  A least squares approach is
 * used in the evaluation of the unlimited solution
 * gradients.  Several slope limiters may be used.
 *
 ********************************************************/
void Linear_Reconstruction_LeastSquares(Levermore1D_UniformMesh *Soln,
                                        const int Number_of_Cells,
                                        const int Limiter) {

    int i, n, n2, n_pts, index[2];
    double u0Min, u0Max, uQuad[2], phi;
    double Dx, DxDx_ave;
    Levermore1D_pState DU, DUDx_ave;

    /* Carry out the limited solution reconstruction in
       each cell. */

    for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
        n_pts = 2;
        index[0] = i-1;
        index[1] = i+1;

        DUDx_ave.zero();
        DxDx_ave = ZERO;

        for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
            Dx = Soln[ index[n2] ].X.x - Soln[i].X.x;
            DU = Soln[ index[n2] ].W - Soln[i].W;
            DUDx_ave += DU*Dx;
            DxDx_ave += Dx*Dx;
        } /* endfor */

        DUDx_ave = DUDx_ave/double(n_pts);
        DxDx_ave = DxDx_ave/double(n_pts);

        Soln[i].dWdx = DUDx_ave/DxDx_ave;

        for ( n = 1 ; n <= Levermore1D_Vector::get_length() ; ++n ) {
           u0Min = Soln[i].W[n];
           u0Max = u0Min;
           for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
              u0Min = min(u0Min, Soln[ index[n2] ].W[n]);
              u0Max = max(u0Max, Soln[ index[n2] ].W[n]);
           } /* endfor */

	   uQuad[0] = Soln[i].W[n] - HALF*Soln[i].dWdx[n]*Soln[i].X.dx;
	   uQuad[1] = Soln[i].W[n] + HALF*Soln[i].dWdx[n]*Soln[i].X.dx;

           switch(Limiter) {
             case LIMITER_ZERO :
	       phi = ZERO;
               break;
             case LIMITER_ONE :
	       phi = ONE;
               break;
             case LIMITER_BARTH_JESPERSEN :
               phi = Limiter_BarthJespersen(uQuad, Soln[i].W[n], u0Min, u0Max, 2);
               break;
             case LIMITER_VENKATAKRISHNAN :
               phi = Limiter_Venkatakrishnan(uQuad, Soln[i].W[n], u0Min, u0Max, 2);
               break;
             case LIMITER_VANLEER :
               phi = Limiter_VanLeer(uQuad, Soln[i].W[n], u0Min, u0Max, 2);
               break;
             case LIMITER_VANALBADA :
               phi = Limiter_VanAlbada(uQuad, Soln[i].W[n], u0Min, u0Max, 2);
               break;
	     default:
               phi = Limiter_BarthJespersen(uQuad, Soln[i].W[n], u0Min, u0Max, 2);
               break;
           } /* endswitch */

	   Soln[i].phi[n] = phi;
        } /* endfor */

    } /* endfor */

    Soln[0].dWdx = Soln[1].phi^Soln[1].dWdx;
    Soln[0].phi.set_all(1.0);

    Soln[Number_of_Cells+1].dWdx = Soln[Number_of_Cells].phi^Soln[Number_of_Cells].dWdx;
    Soln[Number_of_Cells+1].phi.set_all(1.0);

}

/******************************************************//**
 * Routine: dUdt_explicitEuler_upwind
 *
 * This routine updates the solution using a 1st-order
 * explicit Levermore time integration and 1st-order upwind
 * spatial discretization scheme in conjunction with
 * either the Godunov, Roe, Rusanov, HLLE, Linde, or
 * HLLC flux functions.
 *
 ********************************************************/
int dUdt_explicitEuler_upwind(Levermore1D_UniformMesh *Soln,
	                      const int Number_of_Cells,
			      double &dtMin,
			      const double &CFL_Number,
                              const int Flux_Function_Type,
			      const int Local_Time_Stepping) {
    int i, count(0);
    Levermore1D_cState Flux;
    Levermore1D_Vector temp;
    ColumnVector Update(Levermore1D_Vector::get_length());
    ColumnVector RHS(Levermore1D_Vector::get_length());
    DenseMatrix LHS(Levermore1D_Vector::get_length(),
		    Levermore1D_Vector::get_length());

    /* Evaluate the time rate of change of the solution
       (i.e., the solution residuals) using the first-order
       upwind scheme with a variety of flux functions. */

    Soln[0].dUdt.zero();
    for ( i = 0 ; i <= Number_of_Cells ; ++i ) {
        Soln[i+1].dUdt.zero();

        switch(Flux_Function_Type) {
          case FLUX_FUNCTION_HLLE :
            Flux = FluxHLLE(Soln[i].U,
			    Soln[i].A,
			    Soln[i].lambda_min,
			    Soln[i+1].U,
			    Soln[i+1].A,
			    Soln[i+1].lambda_max);
            break;
          case FLUX_FUNCTION_KINETIC :
            Flux = FluxKinetic(Soln[i].A,
			       Soln[i].U[2]/Soln[i].U[1],
			       Soln[i+1].A,
			       Soln[i+1].U[2]/Soln[i+1].U[1]);
            break;
	  default:
	    cout << "Error, bad flux function." << endl;
	    return(1);
            break;
        } /* endswitch */

        Soln[i].dUdt -= Flux/Soln[i].X.dx;
        Soln[i+1].dUdt += Flux/Soln[i+1].X.dx;
    } /* endfor */
    Soln[Number_of_Cells+1].dUdt.zero();

    /* Update conserved and primitive solution and closure weights
       using explicit Euler method. */

    for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
      if ( !Local_Time_Stepping ) Soln[i].dt = dtMin;
      temp = (Soln[i].dUdt + Collision_RHS(Soln[i].U) ); //store here temporarily

      LHS.zero();
      for(int j = 0; j < Levermore1D_Vector::get_length(); ++j) {
	LHS(j,j) = 1/(CFL_Number*Soln[i].dt);
	RHS(j) = temp[j+1];
      }
      LHS -= Soln[i].W.dSdU();

      Solve_LU_Decomposition(LHS,RHS,Update);

      Soln[i].U += Update;
      Update = Soln[i].dUdA_inv * Update; //now update for A
      Soln[i].A += Update;
      Soln[i].update_predicted_moment(Update); //for detector
      Soln[i].calculate_detector();

      if ( ! detector_below_tolerance(Soln[i].detector) ) {
	if(Soln[i].A.set_from_U(Soln[i].U)) { //returns 1 if fail
	  cout << endl << "Error, Cannot resync:" << endl
	       << "U =   " << Soln[i].U << endl
	       << "A =   " << Soln[i].A << endl
	       << "U_A = " <<  Levermore1D_cState(Soln[i].A,Soln[i].U[2]/Soln[i].U[1]) << endl;
	  return 1;
	}
	Soln[i].number_of_resyncs++;
	Soln[i].reset_predicted_moment();
	cout << "%";cout.flush();
	++count;
      }
      Soln[i].W = Levermore1D_pState(Soln[i].U);
      Soln[i].calculate_Hessians();
    } /* endfor */

    /* By default, constant extrapolation boundary
       conditions are applied at either end of the mesh. */

    Soln[0].U = Soln[1].U;
    Soln[0].W = Soln[1].W;
    Soln[0].A = Soln[1].A;

    Soln[Number_of_Cells+1].U = Soln[Number_of_Cells].U;
    Soln[Number_of_Cells+1].W = Soln[Number_of_Cells].W;
    Soln[Number_of_Cells+1].A = Soln[Number_of_Cells].A;
    cout << count; cout.flush();
    /* Solution successfully updated. */

    return (0);

}


///******************************************************//**
// * Routine: dUdt_Hancock
// *
// * This routine updates the solution using the
// * second-order predictor-corrector TVD scheme of
// * Hancock (19??).  Various flux functions may be used.
// *
// ********************************************************/
//int dUdt_Hancock(Levermore1D_UniformMesh *Soln,
//	         const int Number_of_Cells,
//		 double &dtMin,
//		 const double &CFL_Number,
//                 const int Limiter_Type,
//                 const int Flux_Function_Type,
//		 const int Local_Time_Stepping) {
//
//    int i, n_stage;
//    Levermore1D_pState Wl, Wr, A[Levermore1D_Vector::get_length()];
//    Levermore1D_cState Flux;
//
//    /* Perform second-order Hancock predictor-corrector
//       update of solution varibles for new time level. */
//
//    for ( n_stage = 1 ; n_stage <= 2 ; ++n_stage ) {
//
//        /* Apply boundary conditions for each step. */
//
//
//        /* Perform the linear reconstruction within each cell
//           of the computational grid in the predictor step. */
//
//        if ( n_stage == 1 ) Linear_Reconstruction_MUSCL(Soln,
//                            Number_of_Cells, Limiter_Type);
//
//        /* Evaluate the time rate of change of the solution
//           (i.e., the solution residuals) using
//           1) a limited central discretization of the
//              non-conservative form of the equations in
//              the predictor step, and
//           2) a second-order limited upwind discretization
//              and of the conserved equations with
//              Roe flux function in the corrector step. */
//
//        if ( !Local_Time_Stepping && n_stage == 1 ) Soln[0].dt = dtMin;
//        if ( n_stage == 1 ) {
//           Soln[0].Uo = Soln[0].U;
//           Soln[0].dUdt.zero();
//        } /* endif */
//
//        for ( i = 0 ; i <= Number_of_Cells ; ++i ) {
//            if ( !Local_Time_Stepping && n_stage == 1 ) Soln[i+1].dt = dtMin;
//            if ( n_stage == 1 ) {
//               Soln[i+1].Uo = Soln[i+1].U;
//               Soln[i+1].dUdt.zero();
//            } /* endif */
//
//            if (n_stage == 1) {
//	       if ( i > 0) {
//                  A[1] = Levermore1D_pState(Soln[i].W.v, Soln[i].W.d, ZERO);
//                  A[2] = Levermore1D_pState(ZERO, Soln[i].W.v, ONE/Soln[i].W.d);
//                  A[3] = Levermore1D_pState(ZERO, Soln[i].W.g*Soln[i].W.p,
//                                        Soln[i].W.v);
//                  Soln[i].W.d -= HALF*(CFL_Number*Soln[i].dt)*A[1]*
//                                 (Soln[i].phi^Soln[i].dWdx);
//                  Soln[i].W.v -= HALF*(CFL_Number*Soln[i].dt)*A[2]*
//                                 (Soln[i].phi^Soln[i].dWdx);
//                  Soln[i].W.p -= HALF*(CFL_Number*Soln[i].dt)*A[3]*
//                                 (Soln[i].phi^Soln[i].dWdx);
//               } /* endif */
//            } else {
//               /* Evaluate the cell interface flux. */
//               Wl = Soln[i].W + (Soln[i].phi^Soln[i].dWdx)*HALF*Soln[i].X.dx;
//               Wr = Soln[i+1].W - (Soln[i+1].phi^Soln[i+1].dWdx)*HALF*Soln[i+1].X.dx;
//
//               switch(Flux_Function_Type) {
//                 case FLUX_FUNCTION_GODUNOV :
//                   //Flux = FluxGodunov(Wl, Wr);
//                   break;
//                 case FLUX_FUNCTION_ROE :
//                   //Flux = FluxRoe(Wl, Wr);
//                   break;
//                 case FLUX_FUNCTION_RUSANOV :
//                   //Flux = FluxRusanov(Wl, Wr);
//                   break;
//                 case FLUX_FUNCTION_HLLE :
//                   //Flux = FluxHLLE(Wl, Wr);
//                   break;
//                 case FLUX_FUNCTION_LINDE :
//                   //Flux = FluxLinde(Wl, Wr);
//                   break;
//                 case FLUX_FUNCTION_HLLC :
//                   //Flux = FluxHLLC(Wl, Wr);
//                   break;
//                 case FLUX_FUNCTION_OSHER :
//                   //Flux = FluxOsher(Wl, Wr);
//                   break;
//	         default:
//                   //Flux = FluxRoe(Wl, Wr);
//                   break;
//               } /* endswitch */
//
//               /* Evaluate cell-averaged solution changes. */
//               Soln[i].dUdt -= Flux*(CFL_Number*Soln[i].dt)/Soln[i].X.dx;
//               Soln[i+1].dUdt += Flux*(CFL_Number*Soln[i+1].dt)/Soln[i+1].X.dx;
//	    } /* endif */
//
//        } /* endfor */
//
//        Soln[0].dUdt.zero();
//        Soln[Number_of_Cells+1].dUdt.zero();
//
//        /* Update all solution variables after corrector step. */
//
//        if (n_stage == 2) {
//           for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
//               Soln[i].U = Soln[i].Uo + Soln[i].dUdt;
//
//               /* Update the primitive variable solution state. */
//	       Soln[i].W = Levermore1D_pState(Soln[i].U);
//           } /* endfor */
//        } /* endif */
//
//    } /* endfor */
//
//    /* By default, constant extrapolation boundary
//       conditions are applied at either end of the mesh. */
//
//    Soln[0].U = Soln[1].U;
//    Soln[0].W = Soln[1].W;
//
//    Soln[Number_of_Cells+1].U = Soln[Number_of_Cells].U;
//    Soln[Number_of_Cells+1].W = Soln[Number_of_Cells].W;
//
//    /* Solution successfully updated. */
//
//    return (0);
//
//}

/******************************************************//**
 * Routine: dUdt_2stage_2ndOrder_upwind
 *
 * This routine updates the solution using a two-stage
 * second-order explicit time integration scheme
 * and a 2nd-ororder limited upwind spatial
 * discretization scheme with either the Godunov, Roe,
 * Rusanov, HLLE, Linde, or HLLC flux functions.
 *
 ********************************************************/
int dUdt_2stage_2ndOrder_upwind(Levermore1D_UniformMesh *Soln,
	                        const CFD1D_Input_Parameters &IP,
			        double &dtMin,
				const double &CFL_Number,
                                const int Reconstruction_Type,
                                const int Limiter_Type,
                                const int Flux_Function_Type,
			        const int Local_Time_Stepping) {
    int i, n_stage, count(0);
    double us;
    double omega;
    Levermore1D_pState Wl, Wr;
    Levermore1D_cState Ul, Ur;
    Levermore1D_cState Flux;
    Levermore1D_weights Al, Ar;
    Levermore1D_Vector temp;
    ColumnVector Update(Levermore1D_Vector::get_length());
    ColumnVector RHS(Levermore1D_Vector::get_length());
    DenseMatrix LHS(Levermore1D_Vector::get_length(),
		    Levermore1D_Vector::get_length());
    ColumnVector delta_A(Levermore1D_Vector::get_length());
    DenseMatrix dUdA_interface(Levermore1D_Vector::get_length(),
			       Levermore1D_Vector::get_length());

    /* Perform second-order two-stage semi-implicit update of solution
       varibles for new time level. */

    for ( n_stage = 1 ; n_stage <= 2 ; ++n_stage ) {

        /* Evaluate the time step fraction for stage. */

        omega = ONE/double(n_stage);

	if ( IP.Reconstruction_In_Each_Stage == true || n_stage == 1 ){

	  /* Perform the linear reconstruction within each cell
	     of the computational grid for the current stage if
	     Reconstruction_In_Each_Stage is true, or only in
	     the first stage if Reconstruction_In_Each_Stage is false */

	  /* Apply boundary conditions for stage. */
	  //BCs(Soln,IP);  I may need BCs at some point

	  switch(Reconstruction_Type) {
	  case RECONSTRUCTION_MUSCL :
	    Linear_Reconstruction_MUSCL(Soln,
					IP.Number_of_Cells,
					Limiter_Type);
	    break;
	  case RECONSTRUCTION_GREEN_GAUSS :
	    Linear_Reconstruction_GreenGauss(Soln,
					     IP.Number_of_Cells,
					     Limiter_Type);
	    break;
	  case RECONSTRUCTION_LEAST_SQUARES :
	    Linear_Reconstruction_LeastSquares(Soln,
					       IP.Number_of_Cells,
					       Limiter_Type);
	    break;
	  default:
	    cout << "Bad reconstruction type chosen";
	    return 1;
	    break;
	  } /* endswitch */
	}

        /* Evaluate the time rate of change of the solution
           (i.e., the solution residuals) using a second-order
           limited upwind scheme with a variety of flux functions. */

        if ( !Local_Time_Stepping && n_stage == 1 ) Soln[0].dt = dtMin;
        if ( n_stage == 1 ) {
	  Soln[0].Uo = Soln[0].U;
	  Soln[0].Ao = Soln[0].A;
	  Soln[0].dUdt.zero();
	}

        for ( i = 0 ; i <= IP.Number_of_Cells ; ++i ) {
            if ( !Local_Time_Stepping && n_stage == 1 ) Soln[i+1].dt = dtMin;
            if ( n_stage == 1 ) {
               Soln[i+1].Uo = Soln[i+1].U;
               Soln[i+1].Ao = Soln[i+1].A;
               Soln[i+1].dUdt.zero();
            } else {
               Soln[i+1].dUdt = Soln[i+1].dUdt*HALF;
            } /* endif */

            /* Evaluate the cell interface flux. */
	    Wl = Soln[i].W +
	      (Soln[i].phi^Soln[i].dWdx)*HALF*Soln[i].X.dx;
	    Wr = Soln[i+1].W -
	      (Soln[i+1].phi^Soln[i+1].dWdx)*HALF*Soln[i+1].X.dx;

	    Ul = Levermore1D_cState(Wl);
	    Ur = Levermore1D_cState(Wr);

	    us = Soln[i].Ur_old[2]/Soln[i].Ur_old[1];
	    dUdA_interface = Soln[i].Ur_old.d2hda2(Soln[i].Ar_old,us);
	    for(int iii=1;iii<=Levermore1D_Vector::get_length(); ++iii) RHS[iii-1] = Ul[iii]-Soln[i].Ur_old[iii];
	    Solve_LU_Decomposition(dUdA_interface,RHS,delta_A);
	    Al = Soln[i].Ar_old + delta_A;
	    Soln[i].Ar_old = Al;
	    Soln[i].Ur_old = Ul;
	    Soln[i].update_predicted_moment_r(delta_A);
	    Soln[i].calculate_detector_r();

	    us = Soln[i+1].Ul_old[2]/Soln[i+1].Ul_old[1];
	    dUdA_interface = Soln[i+1].Ul_old.d2hda2(Soln[i+1].Al_old,us);
	    for(int iii=1;iii<=Levermore1D_Vector::get_length(); ++iii) RHS[iii-1] = Ur[iii]-Soln[i+1].Ul_old[iii];
	    Solve_LU_Decomposition(dUdA_interface,RHS,delta_A);
	    Ar = Soln[i+1].Al_old + delta_A;
	    Soln[i+1].Al_old = Ar;
	    Soln[i+1].Ul_old = Ur;
	    Soln[i+1].update_predicted_moment_l(delta_A);
	    Soln[i+1].calculate_detector_l();

//	    Al = Soln[i].A + Soln[i].dUdA_inv * (Ul-Soln[i].U);
//	    Ar = Soln[i+1].A + Soln[i+1].dUdA_inv * (Ur-Soln[i+1].U);

	    if ( ! detector_below_tolerance(Soln[i].detector_r) ) {
	      if(Al.set_from_U(Ul)) { //returns 1 if fail
		cout << endl << "Error, Cannot resync at left interface:" << endl
		     << "U =   " << Ul << endl
		     << "A =   " << Al << endl
		     << "U_A = " <<  Levermore1D_cState(Al, Ul[2]/Ul[1]) << endl;
		return 1;
	      }
	      Soln[i].Ar_old = Al;
	      Soln[i].reset_predicted_moment_r();
	      cout << "L";cout.flush();
	      ++count;
	    }
	    if ( ! detector_below_tolerance(Soln[i+1].detector_l) ) {
	      if(Ar.set_from_U(Ur)) { //returns 1 if fail
		cout << endl << "Error, Cannot resync at right interface:" << endl
		     << "U =   " << Ur << endl
		     << "A =   " << Ar << endl
		     << "U_A = " <<  Levermore1D_cState(Ar, Ur[2]/Ur[1]) << endl;
		return 1;
	      }
	      Soln[i+1].Al_old = Ar;
	      Soln[i+1].reset_predicted_moment_l();
	      cout << "R";cout.flush();
	      ++count;
	    }


	    /* Apply the BCs before the flux evaluation */
	    // ***** Left boundary **********
	    if (i == 0){
	      // extrapolation BC (by default)
	      Ul = Ur;
	      Wl = Wr;
	      Al = Ar;
	    }

	    // *****  Right boundary *********
	    if (i == IP.Number_of_Cells){
	      // extrapolation BC (by default)
	      Ur = Ul;
	      Wr = Wl;
	      Ar = Al;
	    }

	    switch(Flux_Function_Type) {
	    case FLUX_FUNCTION_HLLE :
	      Flux = FluxHLLE(Ul,
			      Al,
			      Soln[i].lambda_min,
			      Ur,
			      Ar,
			      Soln[i+1].lambda_max);
	      break;
	    case FLUX_FUNCTION_KINETIC :
	      Flux = FluxKinetic(Al,
				 Ul[2]/Ul[1],
				 Ar,
				 Ur[2]/Ur[1]);
	      break;
	    default:
	      cout << "Error, bad flux function." << endl;
	      return(1);
	      break;
            } /* endswitch */

            /* Evaluate cell-averaged solution changes. */

            Soln[i].dUdt -= Flux*omega/Soln[i].X.dx;
            Soln[i+1].dUdt += Flux*omega/Soln[i+1].X.dx;

        } /* endfor */

        Soln[0].dUdt.zero();
        Soln[IP.Number_of_Cells+1].dUdt.zero();

        /* Update solution variables for this stage. */

        for ( i = 1 ; i <= IP.Number_of_Cells ; ++i ) {
	  temp = (Soln[i].dUdt + Collision_RHS(Soln[i].Uo) ); //store here temporarily

	  LHS.zero();
	  for(int j = 0; j < Levermore1D_Vector::get_length(); ++j) {
	    LHS(j,j) = 1/(CFL_Number*Soln[i].dt);
	    RHS(j) = temp[j+1];
	  }
	  LHS -= Soln[i].Uo.dSdU()*omega;

	  Solve_LU_Decomposition(LHS,RHS,Update);

	  Soln[i].U = Soln[i].Uo + Update;
	  Update = Soln[i].dUdA_inv * Update; //now update for A
	  Soln[i].A = Soln[i].Ao + Update;
	  Soln[i].update_predicted_moment(Update); //for detector
	  Soln[i].calculate_detector();

	  if ( ! detector_below_tolerance(Soln[i].detector) ) {
	    if(Soln[i].A.set_from_U(Soln[i].U)) { //returns 1 if fail
	      cout << endl << "Error, Cannot resync:" << endl
		   << "U =   " << Soln[i].U << endl
		   << "A =   " << Soln[i].A << endl
		   << "U_A = " <<  Levermore1D_cState(Soln[i].A, Soln[i].U[2]/Soln[i].U[1]) << endl;
	      return 1;
	    }
	    Soln[i].number_of_resyncs++;
	    Soln[i].reset_predicted_moment();
	    cout << "%";cout.flush();
	    ++count;
	  }
	  Soln[i].W = Levermore1D_pState(Soln[i].U);
	  Soln[i].calculate_Hessians();
        } /* endfor */

    } /* endfor */
    cout << count; cout.flush();
    /* Solution successfully updated. */

    return (0);
}

/******************************************************//**
 * Routine:  LimitedLinearReconstructionOverDomain
 *
 * Peforms the reconstruction of a limited piecewise
 * linear solution state within each cell of the
 * computational mesh. The input parameters object specifies
 * all the parameters necessary to perform the reconstruction
 * (e.g. method, limiter etc.).
 *
 ********************************************************/
void LimitedLinearReconstructionOverDomain(Levermore1D_UniformMesh *Soln, const CFD1D_Input_Parameters &IP){

  switch(IP.i_Reconstruction) {
  case RECONSTRUCTION_MUSCL :
    Linear_Reconstruction_MUSCL(Soln,
				IP.Number_of_Cells,
				IP.i_Limiter);
    break;
  case RECONSTRUCTION_GREEN_GAUSS :
    Linear_Reconstruction_GreenGauss(Soln,
				     IP.Number_of_Cells,
				     IP.i_Limiter);
    break;
  case RECONSTRUCTION_LEAST_SQUARES :
    Linear_Reconstruction_LeastSquares(Soln,
				       IP.Number_of_Cells,
				       IP.i_Limiter);
    break;
  default:
    throw runtime_error("LimitedLinearReconstructionOverDomain() ERROR: Unknown reconstruction type");
    break;
  } /* endswitch */

}

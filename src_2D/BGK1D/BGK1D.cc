/*! \file BGK1D.cc
  \brief Implementation of subroutines prototyped in BGK1D.h file. */

/* Include 1D BGK solution header file. */
#include "BGK1D.h"

/******************************************************//**
 * Routine: Allocate
 *
 * Allocate memory for 1D BGK equation solution.
 *
 ********************************************************/
BGK1D_UniformMesh* Allocate(BGK1D_UniformMesh *Soln_ptr,
				  const CFD1D_Input_Parameters &IP) {

  int NC;                       // number of cells in the computational domain
  int Nghost; 			// number of ghost cells

  /* Calculate the total number of computational cells */
  // IP.Nghost() : calculates the number of ghost cells based on the order and the method of reconstruction.
  Nghost = IP.Nghost();
  NC = IP.Number_of_Cells + 2 * Nghost;

  /* Allocate memory. */
  Soln_ptr = new BGK1D_UniformMesh[NC];

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
 * Deallocate memory for 1D BGK equation solution.
 *
 ********************************************************/
BGK1D_UniformMesh* Deallocate(BGK1D_UniformMesh *Soln_ptr) {

  /* Deallocate memory. */
  if(Soln_ptr != NULL) {
    delete []Soln_ptr;
    Soln_ptr = NULL;
  }
  /* Return memory location. */
  return(Soln_ptr);
}


/******************************************************//**
 * Routine: Output_Tecplot
 *
 * Writes the solution to specified output stream
 * suitable for plotting with TECPLOT.
 *
 ********************************************************/
void Output_Tecplot(BGK1D_UniformMesh *Soln,
                    const CFD1D_Input_Parameters &IP,
                    const int Number_of_Time_Steps,
                    const double &Time,
	            ostream &out_file) {
    Output_Tecplot(Soln,
		   IP.Number_of_Cells,
		   Number_of_Time_Steps,
		   Time,
		   out_file);
}

/******************************************************//**
 * Routine: Output_Tecplot
 *
 * Writes the solution to specified output stream
 * suitable for plotting with TECPLOT.
 *
 ********************************************************/
void Output_Tecplot(BGK1D_UniformMesh *Soln,
                    const int Number_of_Cells,
                    const int Number_of_Time_Steps,
                    const double &Time,
	            ostream &out_file) {

  int ICl, ICu, num_moments(3);

  // Set the limits of the plotted domain
  ICl = Soln[0].ICl;
  ICu = Soln[0].ICu;

  int i;

  out_file << "TITLE = \"" << CFFC_Name() << ": 1D BGK Solution, "
	   << "Time Step/Iteration Level = "
	   << Number_of_Time_Steps
	   << ", Time = " << Time*THOUSAND << " (ms)\"" << "\n"
	   << "VARIABLES = \"x\" \\ \n"
	   << "\"dx\" \n";
  for(i = 0; i < num_moments; ++i) {
    out_file << "\"random_moment_" << i << "\" \\ \n";
  }
  out_file << "ZONE \n";

  for ( i = ICl ; i <= ICu ; ++i ) {
    out_file << " " << Soln[i].X;
      for(int j=0; j < num_moments; ++j) {
	out_file << " " << Soln[i].V.random_moment(j);
      }
      out_file << endl;
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
void Grid(BGK1D_UniformMesh *Soln,
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

  for ( i = 1 ; i <= TC-1 ; ++i ) {
    // Initialize the coordinate of the centroids
    Soln[i].X.x =  Soln[0].X.x + double(i)*delta_x;
  } /* endfor */

}

/******************************************************//**
 * Routine: ICs
 *
 * Assigns initial conditions and data to the
 * solution variables.
 *
 ********************************************************/
int ICs(BGK1D_UniformMesh *Soln,
         char *gas_ptr,
	 const int i_ICtype,
         const int Number_of_Cells,
	 CFD1D_Input_Parameters &IP) {

  int i;
  BGK1D_Vector  Vl, Vr;

  int ICl, ICu, TC;
  double xmin, xmax;
  double a, b, dx;
  int error(0);

  ICl = Soln[0].ICl;
  ICu = Soln[0].ICu;
  TC = Number_of_Cells+2*Soln[0].Nghost; // total number of cells

  /* Assign the initial data for the IVP of interest. */

  switch(i_ICtype) {
  case IC_SOD :
    error = Vl.discrete_Maxwell_Boltzmann(DENSITY_STDATM, ZERO, PRESSURE_STDATM);
    if(error) return error;

    error = Vr.discrete_Maxwell_Boltzmann(DENSITY_STDATM/EIGHT,
				  ZERO,
				  PRESSURE_STDATM/TEN);
    if(error) return error;

    for ( i = 0 ; i <= TC-1 ; ++i ) {
      if (Soln[i].X.x <= ZERO) {
	Soln[i].V = Vl;
      } else {
	Soln[i].V = Vr;
      } /* end if */
    } /* endfor */
    break;
  case IC_CONSTANT :
  case IC_UNIFORM :
  default:
    error = Vl.discrete_Maxwell_Boltzmann(DENSITY_STDATM, ZERO, PRESSURE_STDATM);
    if(error) return error;
    for ( i = 0 ; i <= TC-1 ; ++i ) {
      Soln[i].V = Vl;
    } /* endfor */
    break;
  } /* endswitch */

  //set relaxation times
  BGK1D_Vector::set_relaxation_time(IP.relaxation_time);
  return 0;
}

/******************************************************//**
 * Routine: CFL
 *
 * Determines the allowable global and local time steps
 * (for explicit Euler time stepping scheme) according
 * to the Courant-Friedrichs-Lewy condition.
 *
 ********************************************************/
double CFL(BGK1D_UniformMesh *Soln,
           const int Number_of_Cells) {

  double l_max= max(fabs(BGK1D_Vector::velocity(0)),
		    fabs(BGK1D_Vector::velocity(BGK1D_Vector::get_length()-1)));

  return Soln[0].X.dx/l_max; //always uniform mesh?
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
void Linear_Reconstruction_MUSCL(BGK1D_UniformMesh *Soln,
                                 const int Number_of_Cells,
                                 const int Limiter) {

    int i, n;
    double phi;
    BGK1D_Vector a, b;

    /* Carry out the limited solution reconstruction in
       each cell. */

    b = (Soln[1].V-Soln[0].V)/Soln[0].X.dx;
    for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
        Soln[i].dVdx = (Soln[i+1].V-Soln[i-1].V)*HALF/
                       Soln[i].X.dx;
        a = b;
        b = (Soln[i+1].V-Soln[i].V)/Soln[i].X.dx;

        for ( n = 0 ; n < BGK1D_Vector::get_length() ; ++n ) {
           switch(Limiter) {
             case LIMITER_ZERO :
	       phi = ZERO;
               break;
             case LIMITER_ONE :
	       phi = ONE;
               break;
             case LIMITER_MINMOD :
	       phi = minmod(a[n], b[n])/(Soln[i].dVdx[n]+sgn(Soln[i].dVdx[n])*TOLER*TOLER);
               break;
             case LIMITER_UMIST :
               phi = minmod(TWO*a[n],
                            TWO*b[n],
                            (a[n]+THREE*b[n])/FOUR,
                            (THREE*a[n]+b[n])/FOUR)/
                     (Soln[i].dVdx[n]+sgn(Soln[i].dVdx[n])*TOLER*TOLER);
               break;
             case LIMITER_DOUBLE_MINMOD :
               phi = minmod(HALF*(a[n]+b[n]),
	                    TWO*a[n],
	                    TWO*b[n])/(Soln[i].dVdx[n]+sgn(Soln[i].dVdx[n])*TOLER*TOLER);
               break;
             case LIMITER_SUPERBEE :
               phi = superbee(a[n], b[n])/
                     (Soln[i].dVdx[n]+sgn(Soln[i].dVdx[n])*TOLER*TOLER);
               break;
             case LIMITER_PHI :
               phi = philimiter(a[n], b[n], 1.50)/
                     (Soln[i].dVdx[n]+sgn(Soln[i].dVdx[n])*TOLER*TOLER);
               break;
             case LIMITER_VANLEER :
               phi = vanleer(a[n], b[n])/
                     (Soln[i].dVdx[n]+sgn(Soln[i].dVdx[n])*TOLER*TOLER);
               break;
             case LIMITER_VANALBADA :
               phi = vanalbada(a[n], b[n], 0.10)/
                     (Soln[i].dVdx[n]+sgn(Soln[i].dVdx[n])*TOLER*TOLER);
               break;
	     default:
               phi = philimiter(a[n], b[n], ONE)/
                     (Soln[i].dVdx[n]+sgn(Soln[i].dVdx[n])*TOLER*TOLER);
               break;
           } /* endswitch */

           Soln[i].phi[n] = phi;
        } /* endfor */

    } /* endfor */

    Soln[0].dVdx = Soln[1].phi^Soln[1].dVdx;
    Soln[0].phi.one();

    Soln[Number_of_Cells+1].dVdx = Soln[Number_of_Cells].phi^Soln[Number_of_Cells].dVdx;
    Soln[Number_of_Cells+1].phi.one();

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
void Linear_Reconstruction_GreenGauss(BGK1D_UniformMesh *Soln,
                                      const int Number_of_Cells,
                                      const int Limiter) {

    int i, n;
    double v0Min, v0Max, vQuad[2], phi;

    /* Carry out the limited solution reconstruction in
       each cell. */

    for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
        Soln[i].dVdx = (Soln[i+1].V-Soln[i-1].V)*HALF/
                       Soln[i].X.dx;

        for ( n = 0 ; n < BGK1D_Vector::get_length() ; ++n ) {
	   v0Min = min(Soln[i-1].V[n], Soln[i].V[n]);
	   v0Min = min(v0Min, Soln[i+1].V[n]);
	   v0Max = max(Soln[i-1].V[n], Soln[i].V[n]);
	   v0Max = max(v0Max, Soln[i+1].V[n]);
	   vQuad[0] = Soln[i].V[n] - HALF*Soln[i].dVdx[n]*Soln[i].X.dx;
	   vQuad[1] = Soln[i].V[n] + HALF*Soln[i].dVdx[n]*Soln[i].X.dx;

           switch(Limiter) {
             case LIMITER_ZERO :
	       phi = ZERO;
               break;
             case LIMITER_ONE :
	       phi = ONE;
               break;
             case LIMITER_BARTH_JESPERSEN :
               phi = Limiter_BarthJespersen(vQuad, Soln[i].V[n], v0Min, v0Max, 2);
               break;
             case LIMITER_VENKATAKRISHNAN :
               phi = Limiter_Venkatakrishnan(vQuad, Soln[i].V[n], v0Min, v0Max, 2);
               break;
             case LIMITER_VANLEER :
               phi = Limiter_VanLeer(vQuad, Soln[i].V[n], v0Min, v0Max, 2);
               break;
             case LIMITER_VANALBADA :
               phi = Limiter_VanAlbada(vQuad, Soln[i].V[n], v0Min, v0Max, 2);
               break;
	     default:
               phi = Limiter_BarthJespersen(vQuad, Soln[i].V[n], v0Min, v0Max, 2);
               break;
           } /* endswitch */

           Soln[i].phi[n] = phi;
        } /* endfor */

    } /* endfor */

    Soln[0].dVdx = Soln[1].phi^Soln[1].dVdx;
    Soln[0].phi.one();

    Soln[Number_of_Cells+1].dVdx = Soln[Number_of_Cells].phi^Soln[Number_of_Cells].dVdx;
    Soln[Number_of_Cells+1].phi.one();

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
void Linear_Reconstruction_LeastSquares(BGK1D_UniformMesh *Soln,
                                        const int Number_of_Cells,
                                        const int Limiter) {

    int i, n, n2, n_pts, index[2];
    double v0Min, v0Max, vQuad[2], phi;
    double Dx, DxDx_ave;
    BGK1D_Vector DV, DVDx_ave;

    /* Carry out the limited solution reconstruction in
       each cell. */

    for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
        n_pts = 2;
        index[0] = i-1;
        index[1] = i+1;

        DVDx_ave.zero();
        DxDx_ave = ZERO;

        for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
            Dx = Soln[ index[n2] ].X.x - Soln[i].X.x;
            DV = Soln[ index[n2] ].V - Soln[i].V;
            DVDx_ave += Dx*DV;
            DxDx_ave += Dx*Dx;
        } /* endfor */

        DVDx_ave = DVDx_ave/double(n_pts);
        DxDx_ave = DxDx_ave/double(n_pts);

        Soln[i].dVdx = DVDx_ave/DxDx_ave;

        for ( n = 0 ; n < BGK1D_Vector::get_length() ; ++n ) {
           v0Min = Soln[i].V[n];
           v0Max = v0Min;
           for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
              v0Min = min(v0Min, Soln[ index[n2] ].V[n]);
              v0Max = max(v0Max, Soln[ index[n2] ].V[n]);
           } /* endfor */

	   vQuad[0] = Soln[i].V[n] - HALF*Soln[i].dVdx[n]*Soln[i].X.dx;
	   vQuad[1] = Soln[i].V[n] + HALF*Soln[i].dVdx[n]*Soln[i].X.dx;

           switch(Limiter) {
             case LIMITER_ZERO :
	       phi = ZERO;
               break;
             case LIMITER_ONE :
	       phi = ONE;
               break;
             case LIMITER_BARTH_JESPERSEN :
               phi = Limiter_BarthJespersen(vQuad, Soln[i].V[n], v0Min, v0Max, 2);
               break;
             case LIMITER_VENKATAKRISHNAN :
               phi = Limiter_Venkatakrishnan(vQuad, Soln[i].V[n], v0Min, v0Max, 2);
               break;
             case LIMITER_VANLEER :
               phi = Limiter_VanLeer(vQuad, Soln[i].V[n], v0Min, v0Max, 2);
               break;
             case LIMITER_VANALBADA :
               phi = Limiter_VanAlbada(vQuad, Soln[i].V[n], v0Min, v0Max, 2);
               break;
	     default:
               phi = Limiter_BarthJespersen(vQuad, Soln[i].V[n], v0Min, v0Max, 2);
               break;
           } /* endswitch */

	   Soln[i].phi[n] = phi;
        } /* endfor */

    } /* endfor */

    Soln[0].dVdx = Soln[1].phi^Soln[1].dVdx;
    Soln[0].phi.one();

    Soln[Number_of_Cells+1].dVdx = Soln[Number_of_Cells].phi^Soln[Number_of_Cells].dVdx;
    Soln[Number_of_Cells+1].phi.one();

}

/******************************************************//**
 * Routine: dUdt_explicitEuler_upwind
 *
 * This routine updates the solution using a 1st-order
 * explicit BGK time integration and 1st-order upwind
 * spatial discretization scheme in conjunction with
 * either the Godunov, Roe, Rusanov, HLLE, Linde, or
 * HLLC flux functions.
 *
 ********************************************************/
int dVdt_explicitEuler_upwind(BGK1D_UniformMesh *Soln,
	                      const int Number_of_Cells,
			      double &dtMin,
			      const double &CFL_Number,
                              const int Flux_Function_Type,
			      const int Local_Time_Stepping) {
//    int i, count(0);
//    BGK1D_cState Flux;
//    BGK1D_Vector temp;
//    ColumnVector Update(BGK1D_Vector::get_length());
//    ColumnVector RHS(BGK1D_Vector::get_length());
//    DenseMatrix LHS(BGK1D_Vector::get_length(),
//		    BGK1D_Vector::get_length());
//
//    /* Evaluate the time rate of change of the solution
//       (i.e., the solution residuals) using the first-order
//       upwind scheme with a variety of flux functions. */
//
//    Soln[0].dUdt.zero();
//    for ( i = 0 ; i <= Number_of_Cells ; ++i ) {
//        Soln[i+1].dUdt.zero();
//
//        switch(Flux_Function_Type) {
//          case FLUX_FUNCTION_HLLE :
//            Flux = FluxHLLE(Soln[i].U,
//			    Soln[i].A,
//			    Soln[i].lambda_min,
//			    Soln[i+1].U,
//			    Soln[i+1].A,
//			    Soln[i+1].lambda_max);
//            break;
//          case FLUX_FUNCTION_KINETIC :
//            Flux = FluxKinetic(Soln[i].A,
//			       Soln[i].U[2]/Soln[i].U[1],
//			       Soln[i+1].A,
//			       Soln[i+1].U[2]/Soln[i+1].U[1]);
//            break;
//	  default:
//	    cout << "Error, bad flux function." << endl;
//	    return(1);
//            break;
//        } /* endswitch */
//
//        Soln[i].dUdt -= Flux/Soln[i].X.dx;
//        Soln[i+1].dUdt += Flux/Soln[i+1].X.dx;
//    } /* endfor */
//    Soln[Number_of_Cells+1].dUdt.zero();
//
//    /* Update conserved and primitive solution and closure weights
//       using explicit Euler method. */
//
//    for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
//      if ( !Local_Time_Stepping ) Soln[i].dt = dtMin;
//      temp = (Soln[i].dUdt + Collision_RHS(Soln[i].U) ); //store here temporarily
//
//      LHS.zero();
//      for(int j = 0; j < BGK1D_Vector::get_length(); ++j) {
//	LHS(j,j) = 1/(CFL_Number*Soln[i].dt);
//	RHS(j) = temp[j+1];
//      }
//      LHS -= Soln[i].W.dSdU();
//
//      Solve_LU_Decomposition(LHS,RHS,Update);
//
//      Soln[i].U += Update;
//      Update = Soln[i].dUdA_inv * Update; //now update for A
//      Soln[i].A += Update;
//      Soln[i].update_predicted_moment(Update); //for detector
//      Soln[i].calculate_detector();
//
//      if ( ! detector_below_tolerance(Soln[i].detector) ) {
//	if(Soln[i].A.set_from_U(Soln[i].U)) { //returns 1 if fail
//	  cout << endl << "Error, Cannot resync:" << endl
//	       << "U =   " << Soln[i].U << endl
//	       << "A =   " << Soln[i].A << endl
//	       << "U_A = " <<  BGK1D_cState(Soln[i].A,Soln[i].U[2]/Soln[i].U[1]) << endl;
//	  return 1;
//	}
//	Soln[i].number_of_resyncs++;
//	Soln[i].reset_predicted_moment();
//	cout << "%";cout.flush();
//	++count;
//      }
//      Soln[i].W = BGK1D_pState(Soln[i].U);
//      Soln[i].calculate_Hessians();
//    } /* endfor */
//
//    /* By default, constant extrapolation boundary
//       conditions are applied at either end of the mesh. */
//
//    Soln[0].U = Soln[1].U;
//    Soln[0].W = Soln[1].W;
//    Soln[0].A = Soln[1].A;
//
//    Soln[Number_of_Cells+1].U = Soln[Number_of_Cells].U;
//    Soln[Number_of_Cells+1].W = Soln[Number_of_Cells].W;
//    Soln[Number_of_Cells+1].A = Soln[Number_of_Cells].A;
//    cout << count; cout.flush();
//    /* Solution successfully updated. */

    return (0);

}

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
int dVdt_2stage_2ndOrder_upwind(BGK1D_UniformMesh *Soln,
	                        const CFD1D_Input_Parameters &IP,
			        double &dtMin,
				const double &CFL_Number,
                                const int Reconstruction_Type,
                                const int Limiter_Type,
                                const int Flux_Function_Type,
			        const int Local_Time_Stepping) {
//    int i, n_stage, count(0);
//    double us;
//    double omega;
//    BGK1D_pState Wl, Wr;
//    BGK1D_cState Ul, Ur;
//    BGK1D_cState Flux;
//    BGK1D_weights Al, Ar;
//    BGK1D_Vector temp;
//    ColumnVector Update(BGK1D_Vector::get_length());
//    ColumnVector RHS(BGK1D_Vector::get_length());
//    DenseMatrix LHS(BGK1D_Vector::get_length(),
//		    BGK1D_Vector::get_length());
//    ColumnVector delta_A(BGK1D_Vector::get_length());
//    DenseMatrix dUdA_interface(BGK1D_Vector::get_length(),
//			       BGK1D_Vector::get_length());
//
//    /* Perform second-order two-stage semi-implicit update of solution
//       varibles for new time level. */
//
//    for ( n_stage = 1 ; n_stage <= 2 ; ++n_stage ) {
//
//        /* Evaluate the time step fraction for stage. */
//
//        omega = ONE/double(n_stage);
//
//	if ( IP.Reconstruction_In_Each_Stage == true || n_stage == 1 ){
//
//	  /* Perform the linear reconstruction within each cell
//	     of the computational grid for the current stage if
//	     Reconstruction_In_Each_Stage is true, or only in
//	     the first stage if Reconstruction_In_Each_Stage is false */
//
//	  /* Apply boundary conditions for stage. */
//	  //BCs(Soln,IP);  I may need BCs at some point
//
//	  switch(Reconstruction_Type) {
//	  case RECONSTRUCTION_MUSCL :
//	    Linear_Reconstruction_MUSCL(Soln,
//					IP.Number_of_Cells,
//					Limiter_Type);
//	    break;
//	  case RECONSTRUCTION_GREEN_GAUSS :
//	    Linear_Reconstruction_GreenGauss(Soln,
//					     IP.Number_of_Cells,
//					     Limiter_Type);
//	    break;
//	  case RECONSTRUCTION_LEAST_SQUARES :
//	    Linear_Reconstruction_LeastSquares(Soln,
//					       IP.Number_of_Cells,
//					       Limiter_Type);
//	    break;
//	  default:
//	    cout << "Bad reconstruction type chosen";
//	    return 1;
//	    break;
//	  } /* endswitch */
//	}
//
//        /* Evaluate the time rate of change of the solution
//           (i.e., the solution residuals) using a second-order
//           limited upwind scheme with a variety of flux functions. */
//
//        if ( !Local_Time_Stepping && n_stage == 1 ) Soln[0].dt = dtMin;
//        if ( n_stage == 1 ) {
//	  Soln[0].Uo = Soln[0].U;
//	  Soln[0].Ao = Soln[0].A;
//	  Soln[0].dUdt.zero();
//	}
//
//        for ( i = 0 ; i <= IP.Number_of_Cells ; ++i ) {
//            if ( !Local_Time_Stepping && n_stage == 1 ) Soln[i+1].dt = dtMin;
//            if ( n_stage == 1 ) {
//               Soln[i+1].Uo = Soln[i+1].U;
//               Soln[i+1].Ao = Soln[i+1].A;
//               Soln[i+1].dUdt.zero();
//            } else {
//               Soln[i+1].dUdt = Soln[i+1].dUdt*HALF;
//            } /* endif */
//
//            /* Evaluate the cell interface flux. */
//	    Wl = Soln[i].W +
//	      (Soln[i].phi^Soln[i].dWdx)*HALF*Soln[i].X.dx;
//	    Wr = Soln[i+1].W -
//	      (Soln[i+1].phi^Soln[i+1].dWdx)*HALF*Soln[i+1].X.dx;
//
//	    Ul = BGK1D_cState(Wl);
//	    Ur = BGK1D_cState(Wr);
//
//	    us = Soln[i].Ur_old[2]/Soln[i].Ur_old[1];
//	    dUdA_interface = Soln[i].Ur_old.d2hda2(Soln[i].Ar_old,us);
//	    for(int iii=1;iii<=BGK1D_Vector::get_length(); ++iii) RHS[iii-1] = Ul[iii]-Soln[i].Ur_old[iii];
//	    Solve_LU_Decomposition(dUdA_interface,RHS,delta_A);
//	    Al = Soln[i].Ar_old + delta_A;
//	    Soln[i].Ar_old = Al;
//	    Soln[i].Ur_old = Ul;
//	    Soln[i].update_predicted_moment_r(delta_A);
//	    Soln[i].calculate_detector_r();
//
//	    us = Soln[i+1].Ul_old[2]/Soln[i+1].Ul_old[1];
//	    dUdA_interface = Soln[i+1].Ul_old.d2hda2(Soln[i+1].Al_old,us);
//	    for(int iii=1;iii<=BGK1D_Vector::get_length(); ++iii) RHS[iii-1] = Ur[iii]-Soln[i+1].Ul_old[iii];
//	    Solve_LU_Decomposition(dUdA_interface,RHS,delta_A);
//	    Ar = Soln[i+1].Al_old + delta_A;
//	    Soln[i+1].Al_old = Ar;
//	    Soln[i+1].Ul_old = Ur;
//	    Soln[i+1].update_predicted_moment_l(delta_A);
//	    Soln[i+1].calculate_detector_l();
//
////	    Al = Soln[i].A + Soln[i].dUdA_inv * (Ul-Soln[i].U);
////	    Ar = Soln[i+1].A + Soln[i+1].dUdA_inv * (Ur-Soln[i+1].U);
//
//	    if ( ! detector_below_tolerance(Soln[i].detector_r) ) {
//	      if(Al.set_from_U(Ul)) { //returns 1 if fail
//		cout << endl << "Error, Cannot resync at left interface:" << endl
//		     << "U =   " << Ul << endl
//		     << "A =   " << Al << endl
//		     << "U_A = " <<  BGK1D_cState(Al, Ul[2]/Ul[1]) << endl;
//		return 1;
//	      }
//	      Soln[i].Ar_old = Al;
//	      Soln[i].reset_predicted_moment_r();
//	      cout << "L";cout.flush();
//	      ++count;
//	    }
//	    if ( ! detector_below_tolerance(Soln[i+1].detector_l) ) {
//	      if(Ar.set_from_U(Ur)) { //returns 1 if fail
//		cout << endl << "Error, Cannot resync at right interface:" << endl
//		     << "U =   " << Ur << endl
//		     << "A =   " << Ar << endl
//		     << "U_A = " <<  BGK1D_cState(Ar, Ur[2]/Ur[1]) << endl;
//		return 1;
//	      }
//	      Soln[i+1].Al_old = Ar;
//	      Soln[i+1].reset_predicted_moment_l();
//	      cout << "R";cout.flush();
//	      ++count;
//	    }
//
//
//	    /* Apply the BCs before the flux evaluation */
//	    // ***** Left boundary **********
//	    if (i == 0){
//	      // extrapolation BC (by default)
//	      Ul = Ur;
//	      Wl = Wr;
//	      Al = Ar;
//	    }
//
//	    // *****  Right boundary *********
//	    if (i == IP.Number_of_Cells){
//	      // extrapolation BC (by default)
//	      Ur = Ul;
//	      Wr = Wl;
//	      Ar = Al;
//	    }
//
//	    switch(Flux_Function_Type) {
//	    case FLUX_FUNCTION_HLLE :
//	      Flux = FluxHLLE(Ul,
//			      Al,
//			      Soln[i].lambda_min,
//			      Ur,
//			      Ar,
//			      Soln[i+1].lambda_max);
//	      break;
//	    case FLUX_FUNCTION_KINETIC :
//	      Flux = FluxKinetic(Al,
//				 Ul[2]/Ul[1],
//				 Ar,
//				 Ur[2]/Ur[1]);
//	      break;
//	    default:
//	      cout << "Error, bad flux function." << endl;
//	      return(1);
//	      break;
//            } /* endswitch */
//
//            /* Evaluate cell-averaged solution changes. */
//
//            Soln[i].dUdt -= Flux*omega/Soln[i].X.dx;
//            Soln[i+1].dUdt += Flux*omega/Soln[i+1].X.dx;
//
//        } /* endfor */
//
//        Soln[0].dUdt.zero();
//        Soln[IP.Number_of_Cells+1].dUdt.zero();
//
//        /* Update solution variables for this stage. */
//
//        for ( i = 1 ; i <= IP.Number_of_Cells ; ++i ) {
//	  temp = (Soln[i].dUdt + Collision_RHS(Soln[i].Uo) ); //store here temporarily
//
//	  LHS.zero();
//	  for(int j = 0; j < BGK1D_Vector::get_length(); ++j) {
//	    LHS(j,j) = 1/(CFL_Number*Soln[i].dt);
//	    RHS(j) = temp[j+1];
//	  }
//	  LHS -= Soln[i].Uo.dSdU()*omega;
//
//	  Solve_LU_Decomposition(LHS,RHS,Update);
//
//	  Soln[i].U = Soln[i].Uo + Update;
//	  Update = Soln[i].dUdA_inv * Update; //now update for A
//	  Soln[i].A = Soln[i].Ao + Update;
//	  Soln[i].update_predicted_moment(Update); //for detector
//	  Soln[i].calculate_detector();
//
//	  if ( ! detector_below_tolerance(Soln[i].detector) ) {
//	    if(Soln[i].A.set_from_U(Soln[i].U)) { //returns 1 if fail
//	      cout << endl << "Error, Cannot resync:" << endl
//		   << "U =   " << Soln[i].U << endl
//		   << "A =   " << Soln[i].A << endl
//		   << "U_A = " <<  BGK1D_cState(Soln[i].A, Soln[i].U[2]/Soln[i].U[1]) << endl;
//	      return 1;
//	    }
//	    Soln[i].number_of_resyncs++;
//	    Soln[i].reset_predicted_moment();
//	    cout << "%";cout.flush();
//	    ++count;
//	  }
//	  Soln[i].W = BGK1D_pState(Soln[i].U);
//	  Soln[i].calculate_Hessians();
//        } /* endfor */
//
//    } /* endfor */
//    cout << count; cout.flush();
//    /* Solution successfully updated. */
//
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
void LimitedLinearReconstructionOverDomain(BGK1D_UniformMesh *Soln, const CFD1D_Input_Parameters &IP){

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

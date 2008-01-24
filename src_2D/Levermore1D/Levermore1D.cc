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
	   << "\"dx\" \\ \n";
  for(i = 0; i < Levermore1D_Vector::get_length(); ++i) {
    out_file << "\"random_moment" << i << "\" \\ \n";
  }
  out_file << "ZONE \n";

  for ( i = ICl ; i <= ICu ; ++i ) {
    out_file << " " << Soln[i] << "\n";
  } /* endfor */

  out_file << "\n";

}

///******************************************************//**
// * Routine: Output_Tecplot_HighOrder
// *
// * Writes the solution to specified output stream
// * suitable for plotting with TECPLOT.
// *
// ********************************************************/
//void Output_Tecplot_HighOrder(Levermore1D_UniformMesh *Soln,
//			      const int Number_of_Cells,
//			      const int Number_of_Time_Steps,
//			      const double &Time,
//			      ostream &out_file) {
//
//
//  int ICl, ICu;
//
//  // Set the limits of the plotted domain
//  ICl = Soln[0].ICl;
//  ICu = Soln[0].ICu;
//
//  int i;
//  out_file << "TITLE = \"" << CFFC_Name() << ": 1D Levermore Solution, "
//	   << "Time Step/Iteration Level = "
//	   << Number_of_Time_Steps
//	   << ", Time = " << Time*THOUSAND << " (ms)\"" << "\n"
//	   << "VARIABLES = \"x\" \\ \n"
//	   << "\"dx\" \\ \n"
//	   << "\"rho\" \\ \n"
//	   << "\"u\" \\ \n"
//	   << "\"p\" \\ \n"
//	   << "\"T\" \\ \n"
//	   << "\"ValISrho\" \\ \n"
//	   << "\"ISrho\" \\ \n"
//	   << "\"ValISu\" \\ \n"
//	   << "\"ISu\" \\ \n"
//	   << "\"ValISp\" \\ \n"
//	   << "\"ISp\" \\ \n"
//	   << "ZONE \n";
//
//  for ( i = ICl ; i <= ICu ; ++i ) {
//    out_file << " " << Soln[i] << "\n"
//	     << Soln[i].CellHighOrder().CellSmoothnessIndicator(1) <<" "
//	     << Soln[i].CellHighOrder().CellInadequateFit(1) <<" "
//	     << Soln[i].CellHighOrder().CellSmoothnessIndicator(2) <<" "
//	     << Soln[i].CellHighOrder().CellInadequateFit(2) <<" "
//	     << Soln[i].CellHighOrder().CellSmoothnessIndicator(3) <<" "
//	     << Soln[i].CellHighOrder().CellInadequateFit(3) <<" " << "\n";
//  } /* endfor */
//
//  out_file << "\n";
//
//}

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
  Levermore1D_weights Al,Ar,Am;

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
    Al = Levermore1D_weights(DENSITY_STDATM, ZERO, PRESSURE_STDATM);
    Ar = Levermore1D_weights(DENSITY_STDATM/EIGHT,
			     ZERO,
			     PRESSURE_STDATM/TEN);
    for ( i = 0 ; i <= TC-1 ; ++i ) {
      if (Soln[i].X.x <= ZERO) {
	Soln[i].set_state(Al);
      } else {
	Soln[i].set_state(Ar);
      } /* end if */
    } /* endfor */
    break;
//  case IC_GROTH :
//    Wl = Levermore1D_pState(4.696, ZERO, 404.4e03);
//    Wr = Levermore1D_pState(1.408, ZERO, 101.1e03);
//    for ( i = 0 ; i <= TC-1 ; ++i ) {
//      if (Soln[i].X.x <= ZERO) {
//	Soln[i].W = Wl;
//      } else {
//	Soln[i].W = Wr;
//      } /* end if */
//      Soln[i].U = U(Soln[i].W);
//    } /* endfor */
//    break;
//  case IC_EINFELDT :
//    Wl = Levermore1D_pState(ONE,-TWO,FOUR/TEN);
//    Wr = Levermore1D_pState(ONE,TWO,FOUR/TEN);
//    for ( i = 0 ; i <= TC-1 ; ++i ) {
//      if (Soln[i].X.x <= ZERO) {
//	Soln[i].W = Wl;
//      } else {
//	Soln[i].W = Wr;
//      } /* end if */
//      Soln[i].U = U(Soln[i].W);
//    } /* endfor */
//    break;
//  case IC_SHOCK_WAVE :
//    Wl = Levermore1D_pState(2.281, 164.83, 201.17e03);
//    Wr = Levermore1D_pState(1.408, ZERO, 101.1e03);
//    for ( i = 0 ; i <= TC-1 ; ++i ) {
//      if (Soln[i].X.x <= ZERO) {
//	Soln[i].W = Wl;
//      } else {
//	Soln[i].W = Wr;
//      } /* end if */
//      Soln[i].U = U(Soln[i].W);
//    } /* endfor */
//    break;
//  case IC_CONTACT_SURFACE :
//    Wl = Levermore1D_pState(1.045, 200.00, 300.00e03);
//    Wr = Levermore1D_pState(3.483, 200.00, 300.00e03);
//    for ( i = 0 ; i <= TC-1 ; ++i ) {
//      if (Soln[i].X.x <= ZERO) {
//	Soln[i].W = Wl;
//      } else {
//	Soln[i].W = Wr;
//      } /* end if */
//      Soln[i].U = U(Soln[i].W);
//    } /* endfor */
//    break;
//  case IC_RAREFACTION_WAVE :
//    Wl = Levermore1D_pState(1.598, -383.64, 91.88e03);
//    Wr = Levermore1D_pState(2.787, -216.97, 200.0e03);
//    for ( i = 0 ; i <= TC-1 ; ++i ) {
//      if (Soln[i].X.x <= ZERO) {
//	Soln[i].W = Wl;
//      } else {
//	Soln[i].W = Wr;
//      } /* end if */
//      Soln[i].U = U(Soln[i].W);
//    } /* endfor */
//    break;
//  case IC_LAX :
//    Wl = Levermore1D_pState(0.445, 0.698, 3.528);
//    Wr = Levermore1D_pState(0.5, 0.0, 0.571);
//    for ( i = 0 ; i <= TC-1 ; ++i ) {
//      if (Soln[i].X.x <= ZERO) {
//	Soln[i].W = Wl;
//      } else {
//	Soln[i].W = Wr;
//      } /* end if */
//      Soln[i].U = U(Soln[i].W);
//    } /* endfor */
//    break;
//  case IC_SIN_WAVE:
//    IP.ExactFunction = SIN_WAVE_Solution;
//    IP.X_ExactSolution_Min = -ONE;
//    IP.X_ExactSolution_Max = ONE;
//    xmax = Soln[ICu].X.x+Soln[ICu].X.dx/2.0;
//    xmin = Soln[ICl].X.x-Soln[ICl].X.dx/2.0;
//    for ( i = 0 ; i <= TC-1 ; ++i ) {
//      Soln[i].W = Levermore1D_W_STDATM;
//      dx = Soln[i].X.dx;
//      a = Soln[i].X.x-dx/2.0;
//      b = Soln[i].X.x+dx/2.0;
//      Soln[i].W[1] = (Integral_SIN_WAVE(xmin,xmax,b)-
//		      Integral_SIN_WAVE(xmin,xmax,a))/dx;
//      Soln[i].W[2] = 100.0;
//      Soln[i].U = U(Soln[i].W);
//    } /* endfor */
//    break;
//  case IC_JIANG_WAVE:
//    IP.ExactFunction = JIANG_IVP_Solution;
//    IP.X_ExactSolution_Min = -ONE;
//    IP.X_ExactSolution_Max = ONE;
//    double dist;
//    xmax = Soln[ICu].X.x+Soln[ICu].X.dx/2.0;
//    xmin = Soln[ICl].X.x-Soln[ICl].X.dx/2.0;
//
//    for ( i = 0 ; i <= TC-1 ; ++i ) {
//      Soln[i].W = Levermore1D_W_STDATM;
//      dx = Soln[i].X.dx;
//      a = Soln[i].X.x-dx/2.0;
//      a = ConvertDomainToMinusOneOne(xmin,xmax,a);
//      b = Soln[i].X.x+dx/2.0;
//      b = ConvertDomainToMinusOneOne(xmin,xmax,b);
//      dist = fabs(b-a);
//
//      // Numeric Integration
//      Soln[i].W[1] = AdaptiveGaussianQuadrature(JIANG_IVP_Solution,a,b,Soln[i].W[1],10)/dist;
//
//      Soln[i].W[2] = 100.0;
//      Soln[i].U = U(Soln[i].W);
//    } /* endfor */
//    break;
//  case IC_SHOCK_ACOUSTIC_INTERACTION:
//    double density;
//    Wl = Levermore1D_pState(3.857143, 2.629369, 10.333333);
//    for ( i = 0 ; i <= TC-1 ; ++i ) {
//      if (Soln[i].X.x <= -4) {
//	Soln[i].W = Wl;
//      } else {
//	density = 1 + 0.2*sin(5*Soln[i].X.x);
//	Wr = Levermore1D_pState(density, ZERO, ONE);
//	Soln[i].W = Wr;
//      } /* end if */
//      Soln[i].U = U(Soln[i].W);
//    } /* endfor */
//    break;
//  case IC_BLAST_WAVE_INTERACTION:
//    // defined in the domain '0' to '1'
//    Wl = Levermore1D_pState(ONE, ZERO, 1.0e03);
//    Wm = Levermore1D_pState(ONE, ZERO, 1.0e-02);
//    Wr = Levermore1D_pState(ONE, ZERO, 1.0e02);
//    for ( i = 0 ; i <= TC-1 ; ++i ) {
//      if (Soln[i].X.x <= 0.1) {
//	Soln[i].W = Wl;
//      } else if (Soln[i].X.x <= 0.9){
//	Soln[i].W = Wm;
//      } else {
//	Soln[i].W = Wr;
//      } /* end if */
//      Soln[i].U = U(Soln[i].W);
//    } /* endfor */
//    break;
//  case IC_CONVECTION_OF_DIFFERENT_SHAPES:
//    // see 'Efficient Implementation of Weighted ENO Schemes', JCP 126, 202-228
//    double xmin, xmax;
//    xmax = Soln[ICu].X.x+Soln[ICu].X.dx/2.0;
//    xmin = Soln[ICl].X.x-Soln[ICl].X.dx/2.0;
//
//    for ( i = 0 ; i <= TC-1 ; ++i ) {
//      a = Soln[i].X.x-Soln[i].X.dx/2.0;
//      a = ConvertDomainToMinusOneOne(xmin,xmax,a);
//      b = Soln[i].X.x+Soln[i].X.dx/2.0;
//      b = ConvertDomainToMinusOneOne(xmin,xmax,b);
//      // compute the average solution in each cell for the density
//      Soln[i].W[1] = AdaptiveGaussianQuadrature(ConvectionShapes,a,b,Soln[i].W[1],14)/Soln[i].X.dx;
//      Soln[i].W[2] = 100.0;	// set velocity
//      // set the conservative variables
//      Soln[i].U = U(Soln[i].W);
//    } /* endfor */
//    break;
  case IC_CONSTANT :
  case IC_UNIFORM :
  default:
    Al = Levermore1D_weights(DENSITY_STDATM, ZERO, PRESSURE_STDATM);
    for ( i = 0 ; i <= TC-1 ; ++i ) {
      Soln[i].set_state(Al);
    } /* endfor */
    break;
  } /* endswitch */

  //Calculate Hessians
  for ( i = 0 ; i <= TC-1 ; ++i ) {
    Soln[i].calculate_Hessians();
  }

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

  /* Determine local and global time steps. */

  dtMin = MILLION;

  for ( i = Soln[0].ICl; i <= Soln[0].ICu ; ++i ) {
    Soln[i].dt = Soln[i].X.dx/
      max(fabs(Soln[i].lambda_max),fabs(Soln[i].lambda_min));
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

///******************************************************//**
// * Routine: Linear_Reconstruction_Characteristic
// *
// * Peforms the reconstruction of a limited piecewise
// * linear solution state within each cell of the
// * computational mesh.  Characteristic reconstruction
// * is used and several slope limiters may be applied.
// *
// ********************************************************/
//void Linear_Reconstruction_Characteristic(Levermore1D_UniformMesh *Soln,
//                                          const int Number_of_Cells,
//                                          const int Limiter) {
//
//    int i, n;
//    double u0Min, u0Max, uQuad[2], phi;
//    Levermore1D_pState Cl, Cm, Cr;
//
//    /* Carry out the limited solution reconstruction in
//       each cell.  Use characteristic variables. */
//
//    Cm = Soln[0].W.C();
//    Cr = Soln[1].W.C();
//
//    for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
//
//        Cl = Cm;
//        Cm = Cr;
//        Cr = Soln[i+1].W.C();
//        Soln[i].dWdx = HALF*(Cr-Cl)/Soln[i].X.dx;
//
//        for ( n = 1 ; n <= Levermore1D_Vector::get_length() ; ++n ) {
//           u0Min = min(Cl[n], Cm[n]);
//           u0Min = min(u0Min, Cr[n]);
//           u0Max = max(Cl[n], Cm[n]);
//           u0Max = max(u0Max, Cr[n]);
//           uQuad[0] = Cm[n] - HALF*Soln[i].dWdx[n]*Soln[i].X.dx;
//           uQuad[1] = Cm[n] + HALF*Soln[i].dWdx[n]*Soln[i].X.dx;
//
//           switch(n) {
//             case 1 :
//               //phi = Limiter_BarthJespersen(uQuad, Soln[i].W[n], u0Min, u0Max, 2);
//               //phi = Limiter_Venkatakrishnan(uQuad, Cm[n], u0Min, u0Max, 2);
//               phi = Limiter_VanLeer(uQuad, Soln[i].W[n], u0Min, u0Max, 2);
//               //phi = Limiter_VanAlbada(uQuad, Soln[i].W[n], u0Min, u0Max, 2);
//	       Soln[i].phi.d = phi;
//               break;
//             case 2 :
//               //phi = Limiter_BarthJespersen(uQuad, Soln[i].W[n], u0Min, u0Max, 2);
//               //phi = Limiter_Venkatakrishnan(uQuad, Cm[n], u0Min, u0Max, 2);
//               phi = Limiter_VanLeer(uQuad, Soln[i].W[n], u0Min, u0Max, 2);
//               //phi = Limiter_VanAlbada(uQuad, Soln[i].W[n], u0Min, u0Max, 2);
//	       Soln[i].phi.v = phi;
//               break;
//             case 3 :
//               //phi = Limiter_BarthJespersen(uQuad, Soln[i].W[n], u0Min, u0Max, 2);
//               //phi = Limiter_Venkatakrishnan(uQuad, Cm[n], u0Min, u0Max, 2);
//               phi = Limiter_VanLeer(uQuad, Soln[i].W[n], u0Min, u0Max, 2);
//               //phi = Limiter_VanAlbada(uQuad, Soln[i].W[n], u0Min, u0Max, 2);
//	       Soln[i].phi.p = phi;
//               break;
//           } /* endswitch */
//        } /* endfor */
//
//    } /* endfor */
//
//    Soln[0].dWdx = Soln[1].phi^Soln[1].dWdx;
//    Soln[0].phi = Levermore1D_W_ONE;
//
//    Soln[Number_of_Cells+1].dWdx = Soln[Number_of_Cells].phi^Soln[Number_of_Cells].dWdx;
//    Soln[Number_of_Cells+1].phi = Levermore1D_W_ONE;
//
//}

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
    int i;
    Levermore1D_cState Flux;
    Levermore1D_cState Update;

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
      Update = Soln[i].dUdt*(CFL_Number*Soln[i].dt);
      Soln[i].U += Update;
      Soln[i].A += Soln[i].dUdA_inv * Update;
      if ( fabs((2.0*Soln[i].A[3]*Soln[i].U[2] + Soln[i].A[2]*Soln[i].U[1])*(Soln[i].U[2]/Soln[i].U[1])) > 1.0e-4 ||
	   fabs((2.0*Soln[i].A[3]*Soln[i].U[3] + Soln[i].A[2]*Soln[i].U[2] + Soln[i].U[1])/Soln[i].U[1]) > 1.0e-4) {
	Soln[i].A = Levermore1D_weights(Soln[i].U);
	cout << "%";
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

    /* Solution successfully updated. */

    return (0);

}

///******************************************************//**
// * Routine: dUdt_LaxFriedrichs
// *
// * This routine updates the solution using the
// * first-order single-stage Lax-Friedrichs scheme.
// * See Lax (1954).
// *
// ********************************************************/
//int dUdt_LaxFriedrichs(Levermore1D_UniformMesh *Soln,
//	               const int Number_of_Cells,
//		       double &dtMin,
//		       const double &CFL_Number,
//		       const int Local_Time_Stepping) {
//
//    int i;
//    double delta_x;
//
//    /* Assign the mesh space to local variable. */
//
//    delta_x = Soln[0].X.dx;
//
//    /* Evaluate the time rate of change of the solution
//       (i.e., the solution residuals) using the first-order
//       Lax-Friedrichs scheme. */
//
//    Soln[0].dUdt.zero();
//    for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
//        if ( !Local_Time_Stepping ) Soln[i].dt = dtMin;
//        Soln[i].dUdt = - (Soln[i+1].W.F()-Soln[i-1].W.F())/
//                         (TWO*delta_x)
//	               + (Soln[i+1].U-TWO*Soln[i].U+
//		   	       Soln[i-1].U)*HALF/
//                              (CFL_Number*Soln[i].dt);
//    } /* endfor */
//    Soln[Number_of_Cells+1].dUdt.zero();
//
//    /* Now update both conserved and primitive solution
//       variables. */
//
//    for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
//        Soln[i].U += Soln[i].dUdt*(CFL_Number*Soln[i].dt);
//
////	if (Soln[i].U.d   <= ZERO ||
////	    Soln[i].U.E   <= ZERO ||
////	    Soln[i].U.e() <= ZERO ) {
////	    cout << "\n " << CFFC_Name() << " ERROR: Negative Density and/or Energy: \n"
////	         << " node = " << i << "\n U = " << Soln[i].U << "\n dUdt = "
////	         << Soln[i].dUdt << "\n";
////	    return (i);
////	}
//
//	Soln[i].W = Levermore1D_pState(Soln[i].U);
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
//
//    /* Solution successfully updated. */
//
//    return (0);
//
//}

/******************************************************//**
 * Routine: dUdt_LaxWendroff
 *
 * This routine updates the solution using the
 * second-order two-step Lax-Wendroff scheme.
 * See Lax and Wendroff (1960, 1964).
 *
 ********************************************************/
//int dUdt_LaxWendroff(Levermore1D_UniformMesh *Soln,
//	             const int Number_of_Cells,
//		     double &dtMin,
//		     const double &CFL_Number,
//		     const int Local_Time_Stepping) {
//
//    int i;
//    double delta_x;
//    Levermore1D_cState Ul, Ur;
//
//    /* Assign the mesh space to local variable. */
//
//    delta_x = Soln[0].X.dx;
//
//    /* Evaluate the time rate of change of the solution
//       (i.e., the solution residuals) using the second-order
//       Lax-Wendroff scheme. */
//
//    if ( !Local_Time_Stepping ) Soln[0].dt = dtMin;
//    Ur =   HALF*(Soln[1].U + Soln[0].U)
//         - (CFL_Number*Soln[0].dt)*
//           (Soln[1].W.F()-Soln[0].W.F())/(TWO*delta_x);
//
//    Soln[0].dUdt.zero();
//    for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
//        if ( !Local_Time_Stepping ) Soln[i].dt = dtMin;
//
//	/* Evaluate temporary cell interface solution states. */
//        Ul = Ur;
//	Ur =   HALF*(Soln[i+1].U + Soln[i].U)
//             - (CFL_Number*Soln[i].dt)*
//               (Soln[i+1].W.F()-Soln[i].W.F())/(TWO*delta_x);
//
//	/* Evaluate the solution residual. */
//	Soln[i].dUdt = - (Ur.F()-Ul.F())/delta_x;
//
//    } /* endfor */
//    Soln[Number_of_Cells+1].dUdt.zero();
//
//    /* Now update both conserved and primitive solution
//       variables. */
//
//    for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
//        Soln[i].U += Soln[i].dUdt*(CFL_Number*Soln[i].dt);
//
////	if (Soln[i].U.d   <= ZERO ||
////	    Soln[i].U.E   <= ZERO ||
////	    Soln[i].U.e() <= ZERO ) {
////	    cout << "\n " << CFFC_Name() << " ERROR: Negative Density and/or Energy: \n"
////	         << " node = " << i << "\n U = " << Soln[i].U << "\n dUdt = "
////	         << Soln[i].dUdt << "\n";
////	    return (i);
////	}
//
//	Soln[i].W = Levermore1D_pState(Soln[i].U);
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

///******************************************************//**
// * Routine: dUdt_MacCormack
// *
// * This routine updates the solution using the
// * second-order predictor-corrector scheme of
// * MacCormack (1969, 1971).  Artificial viscosity of
// * the type devised by Harten and Zwas (1972) is
// * added to smooth the solution near discontinuities.
// *
// ********************************************************/
//int dUdt_MacCormack(Levermore1D_UniformMesh *Soln,
//	            const int Number_of_Cells,
//		    double &dtMin,
//		    const double &CFL_Number,
//		    const int Local_Time_Stepping) {
//
//    int i;
//    static int Order=1;
//    double delta_x, delta_l, delta_r,
//           theta, theta_r, theta_l,
//           switch_l, switch_r;
//    Levermore1D_cState Ul, Ur;
//
//    /* Detemine order of first-order spatial operators. */
//
//    if (Order) {
//       Order = 0;
//    } else {
//       Order = 1;
//    } /* endif */
//
//    /* Assign the mesh space to local variable. */
//
//    delta_x = Soln[0].X.dx;
//
//    /* Evaluate the time rate of change of the solution
//       (i.e., the solution residuals) for the predictor
//       step using first-order backward/forward spatial
//       differences. */
//
//    Soln[0].dUdt.zero();
//    for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
//        Soln[i].dUdt = - (Soln[i+1-Order].W.F()-
//                          Soln[i-Order].W.F())/delta_x;
//    } /* endfor */
//    Soln[Number_of_Cells+1].dUdt.zero();
//
//    /* Evaluate the time rate of change of the solution
//       (i.e., the solution residuals) for the corrector
//       step using first-order backward/forward spatial
//       differences. */
//
//    if ( !Local_Time_Stepping ) Soln[Order].dt = dtMin;
//    Ur = Soln[Order].U +
//         (CFL_Number*Soln[Order].dt)*Soln[Order].dUdt;
//
//    for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
//        if ( !Local_Time_Stepping ) Soln[i+Order].dt = dtMin;
//
//	/* Evaluate temporary predictor solution states. */
//        Ul = Ur;
//	Ur = Soln[i+Order].U +
//	     (CFL_Number*Soln[i+Order].dt)*Soln[i+Order].dUdt;
//
//	/* Evaluate corrector solution residual. */
//	Soln[i].dUdt =   HALF*Soln[i].dUdt
//	               - HALF*(Ur.F()-Ul.F())/delta_x;
//
//	/* Add artificial dissipation of Harten and Zwas. */
//	if (i == 1) {
//            theta_l = ZERO;
//
//	    delta_l = fabs(Soln[i].U.d-Soln[i-1].U.d)/Soln[i].U.d;
//	    delta_r = fabs(Soln[i+1].U.d-Soln[i].U.d)/Soln[i].U.d;
//            if (delta_l + delta_r < TOLER) {
//	        theta = ZERO;
//	    } else {
//	        theta = fabs((delta_r - delta_l)/(delta_l + delta_r));
//            } /* endif */
//
//	    delta_l = fabs(Soln[i+1].U.d-Soln[i].U.d)/Soln[i+1].U.d;
//	    delta_r = fabs(Soln[i+2].U.d-Soln[i+1].U.d)/Soln[i+1].U.d;
//            if (delta_l + delta_r < TOLER) {
//	        theta_r = ZERO;
//	    } else {
//	        theta_r = fabs((delta_r - delta_l)/(delta_l + delta_r));
//            } /* endif */
//	} else if (i < Number_of_Cells) {
//	    theta_l = theta;
//	    theta = theta_r;
//
//	    delta_l = fabs(Soln[i+1].U.d-Soln[i].U.d)/Soln[i+1].U.d;
//	    delta_r = fabs(Soln[i+2].U.d-Soln[i+1].U.d)/Soln[i+1].U.d;
//            if (delta_l + delta_r < TOLER) {
//	        theta_r = ZERO;
//	    } else {
//	        theta_r = fabs((delta_r - delta_l)/(delta_l + delta_r));
//            } /* endif */
//        } else {
//	    theta_l = theta;
//	    theta = theta_r;
//            theta_r = ZERO;
//        } /* endif */
//        switch_r = max(theta,theta_r);
//        switch_l = max(theta,theta_l);
//
//	Soln[i].dUdt += ((Soln[i+1].U-Soln[i].U)*switch_r -
//		   	 (Soln[i].U-Soln[i-1].U)*switch_l)/
//	                (EIGHT*CFL_Number*Soln[i].dt);
//    } /* endfor */
//
//    /* Update both conserved and primitive solution variables. */
//
//    for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
//
//        /* Update the conserved variable solution state. */
//	Soln[i].U += Soln[i].dUdt*(CFL_Number*Soln[i].dt);
//
//	if (Soln[i].U.d   <= ZERO ||
//	    Soln[i].U.E   <= ZERO ||
//	    Soln[i].U.e() <= ZERO ) {
//	    cout << "\n " << CFFC_Name() << " ERROR: Negative Density and/or Energy: \n"
//	         << " node = " << i << "\n U = " << Soln[i].U << "\n dUdt = "
//	         << Soln[i].dUdt << "\n";
//	    return (i);
//	}
//
//	/* Update the primitive variable solution state. */
//	Soln[i].W = W(Soln[i].U);
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

    int i, n_stage;
    double omega;
    Levermore1D_pState Wl, Wr;
    Levermore1D_cState Flux;

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
//	  case RECONSTRUCTION_CHARACTERISTIC :
//	    Linear_Reconstruction_Characteristic(Soln,
//						 IP.Number_of_Cells,
//						 Limiter_Type);
//	    break;
	  default:
	    Linear_Reconstruction_MUSCL(Soln,
					IP.Number_of_Cells,
					Limiter_Type);
	    break;
	  } /* endswitch */
	}

        /* Evaluate the time rate of change of the solution
           (i.e., the solution residuals) using a second-order
           limited upwind scheme with a variety of flux functions. */

        if ( !Local_Time_Stepping && n_stage == 1 ) Soln[0].dt = dtMin;
        if ( n_stage == 1 ) Soln[0].Uo = Soln[0].U;
        Soln[0].dUdt.zero();

        for ( i = 0 ; i <= IP.Number_of_Cells ; ++i ) {
            if ( !Local_Time_Stepping && n_stage == 1 ) Soln[i+1].dt = dtMin;
            if ( n_stage == 1 ) {
               Soln[i+1].Uo = Soln[i+1].U;
               Soln[i+1].dUdt.zero();
            } else {
               Soln[i+1].dUdt = Soln[i+1].dUdt*HALF;
            } /* endif */

            /* Evaluate the cell interface flux. */
            if (Reconstruction_Type != RECONSTRUCTION_CHARACTERISTIC) {
               Wl = Soln[i].W +
                    (Soln[i].phi^Soln[i].dWdx)*HALF*Soln[i].X.dx;
               Wr = Soln[i+1].W -
                    (Soln[i+1].phi^Soln[i+1].dWdx)*HALF*Soln[i+1].X.dx;
            } else {
//               Wl = CtoW(Soln[i].W.C() +
//                    (Soln[i].phi^Soln[i].dWdx)*Soln[i].X.dx)*HALF;
//               Wr = CtoW(Soln[i+1].W.C() -
//                    (Soln[i+1].phi^Soln[i+1].dWdx)*Soln[i+1].X.dx)*HALF;
            } /* endif */

	    /* Apply the BCs before the flux evaluation */
	    // ***** Left boundary **********
	    if (i == 0){
	      // extrapolation BC (by default)
	      Wl = Wr;

	      // wall BCs
	      if (IP.i_ICs == IC_BLAST_WAVE_INTERACTION){
		Wl[2] = -Wl[2];      // change velocity sign
	      }

	      // periodic BCs
	      if ((IP.i_ICs == IC_SIN_WAVE) || (IP.i_ICs == IC_JIANG_WAVE) ||
		  (IP.i_ICs == IC_DENSITY_STEP_WAVE) || (IP.i_ICs == IC_CONVECTION_OF_DIFFERENT_SHAPES)){
		Wl = Soln[IP.Number_of_Cells].W +
		     (Soln[IP.Number_of_Cells].phi^Soln[IP.Number_of_Cells].dWdx)*HALF*Soln[IP.Number_of_Cells].X.dx;
	      }

	    }

	    // *****  Right boundary *********
	    if (i == IP.Number_of_Cells){
	      // extrapolation BC (by default)
	      Wr = Wl;

	      // wall BCs
	      if (IP.i_ICs == IC_BLAST_WAVE_INTERACTION){
		Wr[2] = -Wr[2];      // change velocity sign
	      }

	      // periodic BCs
	      if ((IP.i_ICs == IC_SIN_WAVE) || (IP.i_ICs == IC_JIANG_WAVE) ||
		  (IP.i_ICs == IC_DENSITY_STEP_WAVE) || (IP.i_ICs == IC_CONVECTION_OF_DIFFERENT_SHAPES)){
		Wr = Soln[1].W -
		     (Soln[1].phi^Soln[1].dWdx)*HALF*Soln[1].X.dx;
	      }
	    }

            switch(Flux_Function_Type) {
              case FLUX_FUNCTION_GODUNOV :
                //Flux = FluxGodunov(Wl, Wr);
                break;
              case FLUX_FUNCTION_ROE :
                //Flux = FluxRoe(Wl, Wr);
                break;
              case FLUX_FUNCTION_RUSANOV :
                //Flux = FluxRusanov(Wl, Wr);
                break;
              case FLUX_FUNCTION_HLLE :
                //Flux = FluxHLLE(Wl, Wr);
                break;
              case FLUX_FUNCTION_LINDE :
                //Flux = FluxLinde(Wl, Wr);
                break;
              case FLUX_FUNCTION_HLLC :
                //Flux = FluxHLLC(Wl, Wr);
                break;
              case FLUX_FUNCTION_OSHER :
                //Flux = FluxOsher(Wl, Wr);
                break;
	      default:
                //Flux = FluxRoe(Wl, Wr);
                break;
            } /* endswitch */

            /* Evaluate cell-averaged solution changes. */

            Soln[i].dUdt -= Flux*(omega*CFL_Number*Soln[i].dt)/Soln[i].X.dx;
            Soln[i+1].dUdt += Flux*(omega*CFL_Number*Soln[i+1].dt)/Soln[i+1].X.dx;

        } /* endfor */

        Soln[0].dUdt.zero();
        Soln[IP.Number_of_Cells+1].dUdt.zero();

        /* Update solution variables for this stage. */

        for ( i = 1 ; i <= IP.Number_of_Cells ; ++i ) {
            Soln[i].U = Soln[i].Uo + Soln[i].dUdt;

//  	    if (Soln[i].U.d   <= ZERO ||
//	        Soln[i].U.E   <= ZERO ||
//	        Soln[i].U.e() <= ZERO ) {
//	        cout << "\n " << CFFC_Name() << " ERROR: Negative Density and/or Energy: \n"
//	             << " node = " << i << "\n U = " << Soln[i].U << "\n dUdt = "
//	             << Soln[i].dUdt << "\n";
//	        return (i);
//	    }

            /* Update the primitive variable solution state. */
	    Soln[i].W = Levermore1D_pState(Soln[i].U);
        } /* endfor */

    } /* endfor */

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
//  case RECONSTRUCTION_CHARACTERISTIC :
//    Linear_Reconstruction_Characteristic(Soln,
//					 IP.Number_of_Cells,
//					 IP.i_Limiter);
//    break;
  default:
    throw runtime_error("LimitedLinearReconstructionOverDomain() ERROR: Unknown reconstruction type");
    break;
  } /* endswitch */

}

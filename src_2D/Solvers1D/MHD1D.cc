/* MHD1D.cc:  Subroutines for 1D MHD Solution Classes. */

/* Include 1D MHD solution header file. */

#ifndef _MHD1D_INCLUDED
#include "MHD1D.h"
#endif // _MHD1D_INCLUDED

/********************************************************
 * Routine: Allocate                                    *
 *                                                      *
 * Allocate memory for 1D MHD equation solution.        *
 *                                                      *
 ********************************************************/
MHD1D_UniformMesh* Allocate(MHD1D_UniformMesh *Soln_ptr,
                            const int Number_of_Cells) {

    /* Allocate memory. */

    Soln_ptr = new MHD1D_UniformMesh[Number_of_Cells+2];

    /* Return memory location. */

    return(Soln_ptr);

}

/********************************************************
 * Routine: Deallocate                                  *
 *                                                      *
 * Deallocate memory for 1D MHD equation solution.      *
 *                                                      *
 ********************************************************/
MHD1D_UniformMesh* Deallocate(MHD1D_UniformMesh *Soln_ptr,
                              const int Number_of_Cells) {

    /* Deallocate memory. */

    delete []Soln_ptr;
    Soln_ptr = NULL;

    /* Return memory location. */

    return(Soln_ptr);

}

/********************************************************
 * Routine: Output_Gnuplot                              *
 *                                                      *
 * Writes the solution to specified output stream       *
 * suitable for plotting with GNUPLOT.                  *
 *                                                      *
 ********************************************************/
void Output_Gnuplot(MHD1D_UniformMesh *Soln,
                    const int Number_of_Cells,
                    const int Number_of_Time_Steps,
                    const double &Time,
	            ostream &out_file) {

    int i;

    out_file << "# " << CFFC_Name() << ": 1D MHD Solution, "
             << "Time Step/Iteration Level = "
             << Number_of_Time_Steps
             << ", Time = " << Time << "\n"
	     << "# cell, x, dx, rho, u, v, w, "
	     << "Bx, By, Bz, B0x, B0y, B0z, p, T\n";
    
    for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
        out_file << " " << i << " " << Soln[i] << "\n";
    } /* endfor */

    out_file << "\n";
    
}

/********************************************************
 * Routine: Output_Tecplot                              *
 *                                                      *
 * Writes the solution to specified output stream       *
 * suitable for plotting with TECPLOT.                  *
 *                                                      *
 ********************************************************/
void Output_Tecplot(MHD1D_UniformMesh *Soln,
                    const int Number_of_Cells,
                    const int Number_of_Time_Steps,
                    const double &Time,
	            ostream &out_file) {

    int i;

    out_file << "TITLE = \"" << CFFC_Name() 
             << ": 1D MHD Solution, "
             << "Time Step/Iteration Level = "
             << Number_of_Time_Steps
             << ", Time = " << Time << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
             << "\"dx\" \\ \n"
             << "\"rho\" \\ \n"
             << "\"u\" \\ \n"
             << "\"v\" \\ \n"
             << "\"w\" \\ \n"
             << "\"Bx\" \\ \n"
             << "\"By\" \\ \n"
             << "\"Bz\" \\ \n"
             << "\"B0x\" \\ \n"
             << "\"B0y\" \\ \n"
             << "\"B0z\" \\ \n"
             << "\"p\" \\ \n"
             << "\"T\" \n"
             << "ZONE \n";
    
    for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
        out_file << " " << Soln[i] << "\n";
    } /* endfor */
    
}

/********************************************************
 * Routine: Grid                                        *
 *                                                      *
 * Generates a uniform mesh and assign the locations of *
 * the cell centers to appropriate solution variables.  *
 *                                                      *
 ********************************************************/
void Grid(MHD1D_UniformMesh *Soln,
          const double &xMin,
	  const double &xMax,
	  const int Number_of_Cells) {

    int i;
    double delta_x;

    /* Determine the mesh spacing. */

    delta_x = (xMax - xMin)/double(Number_of_Cells);
    Soln[0].X.setsize(delta_x);

    /* Create the cells. */

    Soln[0].X.x = xMin - HALF*delta_x;

    for ( i = 1 ; i <= Number_of_Cells+1 ; ++i ) {
        Soln[i].X.x =  Soln[0].X.x + double(i)*delta_x;
    } /* endfor */

}

/********************************************************
 * Routine: ICs                                         *
 *                                                      *
 * Assigns initial conditions and data to the           *
 * solution variables.                                  *
 *                                                      *
 ********************************************************/
void ICs(MHD1D_UniformMesh *Soln,
         char *gas_ptr,
	 const int i_ICtype,
         const int Number_of_Cells) {

    int i;
    MHD1D_pState Wl,Wr;

    /* Assign the gas constants for the gas of interest. */

    Soln[0].W.setgas(gas_ptr);
    Soln[0].U.setgas(gas_ptr);

    /* Assign the initial data for the IVP of interest. */
    
    switch(i_ICtype) {
      case IC_CONSTANT :
        Wl = MHD1D_W_REF;
        for ( i = 0 ; i <= Number_of_Cells+1 ; ++i ) {
            Soln[i].W = Wl;
            Soln[i].U = U(Soln[i].W);
        } /* endfor */
        break;
      case IC_UNIFORM :
        Wl = MHD1D_W_REF;
        for ( i = 0 ; i <= Number_of_Cells+1 ; ++i ) {
            Soln[i].W = Wl;
            Soln[i].U = U(Soln[i].W);
        } /* endfor */
        break;
      case IC_BRIO_WU :
        Wl = MHD1D_pState(ONE,
	  		  ZERO, ZERO, ZERO,
			  THREE/FOUR, ONE, ZERO,
			  ZERO, ZERO, ZERO,
			  ONE);
        Wr = MHD1D_pState(ONE/EIGHT,
			  ZERO, ZERO, ZERO,
			  THREE/FOUR, -ONE, ZERO,
			  ZERO, ZERO, ZERO,
         		  ONE/TEN); 
        for ( i = 0 ; i <= Number_of_Cells+1 ; ++i ) {
	    if (Soln[i].X.x <= ZERO) {
               Soln[i].W = Wl;
            } else {
               Soln[i].W = Wr;	     
            } /* end if */
            Soln[i].U = U(Soln[i].W);
        } /* endfor */
        break;
      case IC_DAI_WOODWARD :
        Wl = MHD1D_pState(1.08,
			  1.2, 0.01, 0.5,
			  0.5641895835477563, 1.0155412503859613, 0.5641895835477563,
			  ZERO, ZERO, ZERO,
			  0.95);
        Wr = MHD1D_pState(1,
			  0, 0, 0,
			  0.5641895835477563, 1.1283791670955126, 0.5641895835477563,
			  ZERO, ZERO, ZERO,
			  1);
        for ( i = 0 ; i <= Number_of_Cells+1 ; ++i ) {
	    if (Soln[i].X.x <= ZERO) {
               Soln[i].W = Wl;
            } else {
               Soln[i].W = Wr;	     
            } /* end if */
            Soln[i].U = U(Soln[i].W);
        } /* endfor */
        break;	
      default:
        Wl = MHD1D_W_REF;
        for ( i = 0 ; i <= Number_of_Cells+1 ; ++i ) {
            Soln[i].W = Wl;
            Soln[i].U = U(Soln[i].W);
        } /* endfor */
        break;
    } /* endswitch */

}

/********************************************************
 * Routine: CFL                                         *
 *                                                      *
 * Determines the allowable global and local time steps *
 * (for explicit Euler time stepping scheme) according  *
 * to the Courant-Friedrichs-Lewy condition.            *
 *                                                      *
 ********************************************************/
double CFL(MHD1D_UniformMesh *Soln,
           const int Number_of_Cells) {

    int i;
    double dtMin;

    /* Determine local and global time steps. */

    dtMin = MILLION;

    for ( i = 0 ; i <= Number_of_Cells+1 ; ++i ) {
        Soln[i].dt = Soln[i].X.dx/
	             (max(fabs(Soln[i].W.lambda(1)),
		          fabs(Soln[i].W.lambda(NUM_VAR_MHD1D))));
        dtMin = min(dtMin, Soln[i].dt);
    } /* endfor */

    /* Return the global time step. */

    return (dtMin);

}

/********************************************************
 * Routine: dUdt_explicitEuler_upwind                   *
 *                                                      *
 * This routine updates the solution using a 1st-order  *
 * explicit Euler time integration and 1st-order upwind *
 * spatial discretization scheme in conjunction with    *
 * either the Roe, Rusanov, HLLE, or Linde flux         *
 * functions.                                           *
 *                                                      *
 ********************************************************/
int dUdt_explicitEuler_upwind(MHD1D_UniformMesh *Soln,
	                      const int Number_of_Cells,
			      double &dtMin,
 			      const double &CFL_Number,
                              const int Flux_Function_Type,
			      const int Local_Time_Stepping) {

    int i;
    MHD1D_cState Flux;

    /* Evaluate the time rate of change of the solution
       (i.e., the solution residuals) using the first-order
       upwind scheme with a variety of flux functions. */

    Soln[0].dUdt = MHD1D_U_ZERO;
    for ( i = 0 ; i <= Number_of_Cells ; ++i ) {
        Soln[i+1].dUdt = MHD1D_U_ZERO;

        switch(Flux_Function_Type) {
          case FLUX_FUNCTION_ROE :
            Flux = FluxRoe(Soln[i].W, Soln[i+1].W);
            break;
          case FLUX_FUNCTION_RUSANOV :
            Flux = FluxRusanov(Soln[i].W, Soln[i+1].W);
            break;
          case FLUX_FUNCTION_HLLE :
            Flux = FluxHLLE(Soln[i].W, Soln[i+1].W);
            break;
          case FLUX_FUNCTION_LINDE :
            Flux = FluxLinde(Soln[i].W, Soln[i+1].W);
            break;
	  default:
            Flux = FluxRoe(Soln[i].W, Soln[i+1].W);
            break;
        } /* endswitch */

        Soln[i].dUdt -= Flux/Soln[i].X.dx;
        Soln[i+1].dUdt += Flux/Soln[i+1].X.dx;
    } /* endfor */
    Soln[Number_of_Cells+1].dUdt = MHD1D_U_ZERO;

    /* Update both conserved and primitive solution
       variables using explicit Euler method. */

    for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
        if ( !Local_Time_Stepping ) Soln[i].dt = dtMin;
	
        Soln[i].U += (CFL_Number*Soln[i].dt)*Soln[i].dUdt;
	
	if (Soln[i].U.d   <= ZERO ||
	    Soln[i].U.E1  <= ZERO ||
	    Soln[i].U.e() <= ZERO ) {
	    cout << "\n " << CFFC_Name() << " ERROR: Negative Density and/or Energy: \n"
	         << " node = " << i << "\n U = " << Soln[i].U << "\n dUdt = "
	         << Soln[i].dUdt << "\n";
	    return (i);
	}
	
	Soln[i].W = W(Soln[i].U);
    } /* endfor */

    /* By default, constant extrapolation boundary
       conditions are applied at either end of the mesh. */

    Soln[0].U = Soln[1].U;
    Soln[0].W = Soln[1].W;

    Soln[Number_of_Cells+1].U = Soln[Number_of_Cells].U;
    Soln[Number_of_Cells+1].W = Soln[Number_of_Cells].W;

    /* Solution successfully updated. */

    return (0);
    
}

/********************************************************
 * Routine: dUdt_LaxFriedrichs                          *
 *                                                      *
 * This routine updates the solution using the          *
 * first-order single-stage Lax-Friedrichs scheme.      *
 * See Lax (1954).                                      *
 *                                                      *
 ********************************************************/
int dUdt_LaxFriedrichs(MHD1D_UniformMesh *Soln,
	               const int Number_of_Cells,
		       double &dtMin,
		       const double &CFL_Number,
		       const int Local_Time_Stepping) {

    int i;
    double delta_x;

    /* Assign the mesh space to local variable. */

    delta_x = Soln[0].X.dx;

    /* Evaluate the time rate of change of the solution
       (i.e., the solution residuals) using the first-order
       Lax-Friedrichs scheme. */

    Soln[0].dUdt = MHD1D_U_ZERO;
    for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
        if ( !Local_Time_Stepping ) Soln[i].dt = dtMin;
	
        Soln[i].dUdt = - (Soln[i+1].W.F()-Soln[i-1].W.F())/
                         (TWO*delta_x)
	               + HALF*(Soln[i+1].U-TWO*Soln[i].U+
		   	       Soln[i-1].U)/
                              (CFL_Number*Soln[i].dt);
    } /* endfor */
    Soln[Number_of_Cells+1].dUdt = MHD1D_U_ZERO;

    /* Now update both conserved and primitive solution
       variables. */

    for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
        Soln[i].U += (CFL_Number*Soln[i].dt)*Soln[i].dUdt;
	
	if (Soln[i].U.d   <= ZERO ||
	    Soln[i].U.E1  <= ZERO ||
	    Soln[i].U.e() <= ZERO ) {
	    cout << "\n " << CFFC_Name() << " ERROR: Negative Density and/or Energy: \n"
	         << " node = " << i << "\n U = " << Soln[i].U << "\n dUdt = "
	         << Soln[i].dUdt << "\n";
	    return (i);
	}
	
	Soln[i].W = W(Soln[i].U);
    } /* endfor */

    /* By default, constant extrapolation boundary
       conditions are applied at either end of the mesh. */

    Soln[0].U = Soln[1].U;
    Soln[0].W = Soln[1].W;

    Soln[Number_of_Cells+1].U = Soln[Number_of_Cells].U;
    Soln[Number_of_Cells+1].W = Soln[Number_of_Cells].W;

    /* Solution successfully updated. */

    return (0);
    
}

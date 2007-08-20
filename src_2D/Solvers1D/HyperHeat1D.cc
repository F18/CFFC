/* HyperHeat1D.cc:  Subroutines for 1D Hyperbolic Heat Equations
                    (Maxwell-Cattaneo Equation) Solution Classes. */

/* Include 1D hyperbolic heat equations header file. */

#ifndef _HYPERHEAT1D_INCLUDED
#include "HyperHeat1D.h"
#endif // _HYPERHEAT1D_INCLUDED

/********************************************************
 * Routine: Allocate                                    *
 *                                                      *
 * Allocate memory for 1D hyperbolic heat equations     *
 * solution.                                            *
 *                                                      *
 ********************************************************/
HyperHeat1D_UniformMesh* Allocate(HyperHeat1D_UniformMesh *Soln_ptr,
                                  const int Number_of_Cells) {

    /* Allocate memory. */

    Soln_ptr = new HyperHeat1D_UniformMesh[Number_of_Cells+2];

    /* Return memory location. */

    return(Soln_ptr);

}

/********************************************************
 * Routine: Deallocate                                  *
 *                                                      *
 * Deallocate memory for 1D hyperbolic heat equations   *
 * solution.                                            *
 *                                                      *
 ********************************************************/
HyperHeat1D_UniformMesh* Deallocate(HyperHeat1D_UniformMesh *Soln_ptr,
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
void Output_Gnuplot(HyperHeat1D_UniformMesh *Soln,
                    const int Number_of_Cells,
                    const int Number_of_Time_Steps,
                    const double &Time,
	            ostream &out_file) {

    int i;

    out_file << "# " << CFFC_Name() 
             << ": 1D Hyperbolic Heat Equation Solution, "
             << "Time Step/Iteration Level = "
             << Number_of_Time_Steps
             << ", Time = " << Time*THOUSAND << "\n"
	     << "# cell, x, dx, T, qx\n";
    
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
void Output_Tecplot(HyperHeat1D_UniformMesh *Soln,
                    const int Number_of_Cells,
                    const int Number_of_Time_Steps,
                    const double &Time,
	            ostream &out_file) {

    int i;

    out_file << "TITLE = \"" << CFFC_Name() 
             << ": 1D Hyperbolic Heat Equation Solution, "
             << "Time Step/Iteration Level = "
             << Number_of_Time_Steps
             << ", Time = " << Time*THOUSAND << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
             << "\"dx\" \\ \n"
             << "\"T\" \\ \n"
             << "\"qx\" \n"
             << "ZONE \n";
    
    for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
        out_file << " " << Soln[i] << "\n";
    } /* endfor */

    out_file << "\n";
    
}

/********************************************************
 * Routine: Grid                                        *
 *                                                      *
 * Generates a uniform mesh and assign the locations of *
 * the cell centers to appropriate solution variables.  *
 *                                                      *
 ********************************************************/
void Grid(HyperHeat1D_UniformMesh *Soln,
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
void ICs(HyperHeat1D_UniformMesh *Soln,
	 const int i_ICtype,
         const double &Kappa,
         const double &Tau,
         const int Number_of_Cells) {

    int i, n_frequency;
    double T0, rod_length, rod_kappa, amplitude, delta_x;
    HyperHeat1D_State Ul,Ur;

    /* Assign the thermal conductivity and relaxation time. */

    Soln[0].U.setcon(Kappa, Tau);

    /* Assign the initial data for the IVP of interest. */
    
    switch(i_ICtype) {
      case IC_CONSTANT :
        Ul = HyperHeat1D_State(ONE, ZERO);
        for ( i = 0 ; i <= Number_of_Cells+1 ; ++i ) {
            Soln[i].U = Ul;
        } /* endfor */
        break;
      case IC_UNIFORM :
        Ul = HyperHeat1D_State(ONE, ZERO);
        for ( i = 0 ; i <= Number_of_Cells+1 ; ++i ) {
            Soln[i].U = Ul;
        } /* endfor */
        break;
      case IC_IMPULSIVE_ROD :
        T0 = ZERO;
        rod_kappa = Kappa;
        delta_x = Soln[0].X.dx;
        rod_length = Soln[Number_of_Cells+1].X.x - Soln[0].X.x - delta_x;
        for ( i = 0 ; i <= Number_of_Cells+1 ; ++i ) {
            Soln[i].U.T = T0;
            Soln[i].U.qx = ZERO;
        } /* endfor */
        T0 = HUNDRED;
        Soln[Number_of_Cells+1].U.T = T0;
        Soln[Number_of_Cells+1].U.qx = ZERO;
        break;
      case IC_SINUSOIDAL_ROD1 :
        T0 = ZERO;
        n_frequency = 1;
        amplitude = HUNDRED;
        rod_kappa = Kappa;
        delta_x = Soln[0].X.dx;
        rod_length = Soln[Number_of_Cells+1].X.x - Soln[0].X.x - delta_x;
        for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
            Soln[i].U.T = T0 + 
                          amplitude*sqr(sin(double(n_frequency)*PI*
                         (Soln[i].X.x-Soln[0].X.x-HALF*delta_x)/rod_length));
            Soln[i].U.qx = -rod_kappa*TWO*amplitude*
                           (double(n_frequency)*PI/rod_length)*
                           sin(double(n_frequency)*PI*(Soln[i].X.x-Soln[0].X.x-HALF*delta_x)/rod_length)*
                           cos(double(n_frequency)*PI*(Soln[i].X.x-Soln[0].X.x-HALF*delta_x)/rod_length);;
        } /* endfor */
        Soln[0].U.T = T0;
        Soln[0].U.qx = ZERO;
        Soln[Number_of_Cells+1].U.T = T0;
        Soln[Number_of_Cells+1].U.qx = ZERO;
        break;
      case IC_SINUSOIDAL_ROD4 :
        T0 = ZERO;
        n_frequency = 4;
        amplitude = HUNDRED;
        rod_kappa = Kappa;
        delta_x = Soln[0].X.dx;
        rod_length = Soln[Number_of_Cells+1].X.x - Soln[0].X.x - delta_x;
        for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
            Soln[i].U.T = T0 + 
                          amplitude*sqr(sin(double(n_frequency)*PI*
                         (Soln[i].X.x-Soln[0].X.x-HALF*delta_x)/rod_length));
            Soln[i].U.qx = -rod_kappa*TWO*amplitude*
                           (double(n_frequency)*PI/rod_length)*
                           sin(double(n_frequency)*PI*(Soln[i].X.x-Soln[0].X.x-HALF*delta_x)/rod_length)*
                           cos(double(n_frequency)*PI*(Soln[i].X.x-Soln[0].X.x-HALF*delta_x)/rod_length);;
        } /* endfor */
        Soln[0].U.T = T0;
        Soln[0].U.qx = ZERO;
        Soln[Number_of_Cells+1].U.T = T0;
        Soln[Number_of_Cells+1].U.qx = ZERO;
        break;
      case IC_RIEMANN_IVP_QX0 :
        Ul = HyperHeat1D_State(TEN, ZERO);
        Ur = HyperHeat1D_State(ONE, ZERO);
        for ( i = 0 ; i <= Number_of_Cells+1 ; ++i ) {
	    if (Soln[i].X.x <= ZERO) {
               Soln[i].U = Ul;
            } else {
               Soln[i].U = Ur;	     
            } /* end if */
        } /* endfor */
        break;
      case IC_RIEMANN_IVP_T0 :
        Ul = HyperHeat1D_State(ZERO, -TWO);
        Ur = HyperHeat1D_State(ZERO, ZERO);
        for ( i = 0 ; i <= Number_of_Cells+1 ; ++i ) {
	    if (Soln[i].X.x <= ZERO) {
               Soln[i].U = Ul;
            } else {
               Soln[i].U = Ur;	     
            } /* end if */
        } /* endfor */
        break;
      case IC_RIEMANN_IVP :
        Ul = HyperHeat1D_State(TEN, -TWO);
        Ur = HyperHeat1D_State(ONE, ZERO);
        for ( i = 0 ; i <= Number_of_Cells+1 ; ++i ) {
	    if (Soln[i].X.x <= ZERO) {
               Soln[i].U = Ul;
            } else {
               Soln[i].U = Ur;	     
            } /* end if */
        } /* endfor */
        break;
      default:
        Ul = HyperHeat1D_State(ONE, ZERO);
        for ( i = 0 ; i <= Number_of_Cells+1 ; ++i ) {
            Soln[i].U = Ul;
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
double CFL(HyperHeat1D_UniformMesh *Soln,
           const int Number_of_Cells) {

    int i;
    double dtMin;

    /* Determine local and global time steps. */

    dtMin = MILLION;

    for ( i = 0 ; i <= Number_of_Cells+1 ; ++i ) {
        Soln[i].dt = Soln[i].X.dx/Soln[i].U.a();
        dtMin = min(dtMin, Soln[i].dt);
    } /* endfor */

    /* Return the global time step. */

    return (dtMin);

}

/********************************************************
 * Routine: Linear_Reconstruction_MUSCL                 *
 *                                                      *
 * Peforms the reconstruction of a limited piecewise    *
 * linear solution state within each cell of the        *
 * computational mesh.  The MUSCL approach of Van Leer  *
 * (1978) is used.  Several slope limiters may be used. *
 *                                                      *
 ********************************************************/
void Linear_Reconstruction_MUSCL(HyperHeat1D_UniformMesh *Soln,
                                 const int Number_of_Cells,
                                 const int Limiter) {

    int i, n;
    double phi;
    HyperHeat1D_State a, b;

    /* Carry out the limited solution reconstruction in
       each cell. */

    b = (Soln[1].U-Soln[0].U)/Soln[0].X.dx;
    for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
        Soln[i].dUdx = HALF*(Soln[i+1].U-Soln[i-1].U)/
                       Soln[i].X.dx;
        a = b;
        b = (Soln[i+1].U-Soln[i].U)/Soln[i].X.dx;

        for ( n = 1 ; n <= NUM_VAR_HYPERHEAT1D ; ++n ) {
           switch(Limiter) {
             case LIMITER_ZERO :
	       phi = ZERO;
               break;
             case LIMITER_ONE :
	       phi = ONE;
               break;
             case LIMITER_MINMOD :
	       phi = minmod(a[n], b[n])/(Soln[i].dUdx[n]+sgn(Soln[i].dUdx[n])*TOLER*TOLER);
               break; 
             case LIMITER_UMIST :
               phi = minmod(TWO*a[n], 
                            TWO*b[n],
                            (a[n]+THREE*b[n])/FOUR,
                            (THREE*a[n]+b[n])/FOUR)/
                     (Soln[i].dUdx[n]+sgn(Soln[i].dUdx[n])*TOLER*TOLER);
               break;
             case LIMITER_DOUBLE_MINMOD :
               phi = minmod(HALF*(a[n]+b[n]), 
	                    TWO*a[n], 
	                    TWO*b[n])/(Soln[i].dUdx[n]+sgn(Soln[i].dUdx[n])*TOLER*TOLER);
               break;
             case LIMITER_SUPERBEE :
               phi = superbee(a[n], b[n])/
                     (Soln[i].dUdx[n]+sgn(Soln[i].dUdx[n])*TOLER*TOLER);
               break; 
             case LIMITER_PHI :
               phi = philimiter(a[n], b[n], 1.50)/
                     (Soln[i].dUdx[n]+sgn(Soln[i].dUdx[n])*TOLER*TOLER);
               break;
             case LIMITER_VANLEER :
               phi = vanleer(a[n], b[n])/
                     (Soln[i].dUdx[n]+sgn(Soln[i].dUdx[n])*TOLER*TOLER);
               break;
             case LIMITER_VANALBADA :
               phi = vanalbada(a[n], b[n], 0.10)/
                     (Soln[i].dUdx[n]+sgn(Soln[i].dUdx[n])*TOLER*TOLER);
               break;
	     default:
               phi = philimiter(a[n], b[n], ONE)/
                     (Soln[i].dUdx[n]+sgn(Soln[i].dUdx[n])*TOLER*TOLER);
               break;
           } /* endswitch */

           Soln[i].phi[n] = phi;
        } /* endfor */

    } /* endfor */

    Soln[0].dUdx = Soln[1].phi^Soln[1].dUdx;
    Soln[0].phi = HyperHeat1D_U_ONE;

    Soln[Number_of_Cells+1].dUdx = Soln[Number_of_Cells].phi^Soln[Number_of_Cells].dUdx;
    Soln[Number_of_Cells+1].phi = HyperHeat1D_U_ONE;

}

/********************************************************
 * Routine: Linear_Reconstruction_GreenGauss            *
 *                                                      *
 * Peforms the reconstruction of a  limited piecewise   *
 * linear solution state within each cell of the        *
 * computational mesh.  A Green-Gauss approach is used  *
 * in the evaluation of the unlimited solution          *
 * gradients.  Several slope limiters may be used.      *
 *                                                      *
 ********************************************************/
void Linear_Reconstruction_GreenGauss(HyperHeat1D_UniformMesh *Soln,
                                      const int Number_of_Cells,
                                      const int Limiter) {

    int i, n;
    double u0Min, u0Max, uQuad[2], phi;

    /* Carry out the limited solution reconstruction in
       each cell. */

    for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
        Soln[i].dUdx = HALF*(Soln[i+1].U-Soln[i-1].U)/
                       Soln[i].X.dx;

        for ( n = 1 ; n <= NUM_VAR_HYPERHEAT1D ; ++n ) {
           u0Min = min(Soln[i-1].U[n], Soln[i].U[n]);
           u0Min = min(u0Min, Soln[i+1].U[n]);
           u0Max = max(Soln[i-1].U[n], Soln[i].U[n]);
           u0Max = max(u0Max, Soln[i+1].U[n]);
           uQuad[0] = Soln[i].U[n] - HALF*Soln[i].dUdx[n]*Soln[i].X.dx;
           uQuad[1] = Soln[i].U[n] + HALF*Soln[i].dUdx[n]*Soln[i].X.dx;

           switch(Limiter) {
             case LIMITER_ZERO :
	       phi = ZERO;
               break;
             case LIMITER_ONE :
	       phi = ONE;
               break;
             case LIMITER_BARTH_JESPERSEN :
               phi = Limiter_BarthJespersen(uQuad, Soln[i].U[n], u0Min, u0Max, 2);
               break;
             case LIMITER_VENKATAKRISHNAN :
               phi = Limiter_Venkatakrishnan(uQuad, Soln[i].U[n], u0Min, u0Max, 2);
               break;
             case LIMITER_VANLEER :
               phi = Limiter_VanLeer(uQuad, Soln[i].U[n], u0Min, u0Max, 2);
               break;
             case LIMITER_VANALBADA :
               phi = Limiter_VanAlbada(uQuad, Soln[i].U[n], u0Min, u0Max, 2);
               break;
	     default:
               phi = Limiter_BarthJespersen(uQuad, Soln[i].U[n], u0Min, u0Max, 2);
               break;
           } /* endswitch */

	   Soln[i].phi[n] = phi;
        } /* endfor */

    } /* endfor */

    Soln[0].dUdx = Soln[1].phi^Soln[1].dUdx;
    Soln[0].phi = HyperHeat1D_U_ONE;

    Soln[Number_of_Cells+1].dUdx = Soln[Number_of_Cells].phi^Soln[Number_of_Cells].dUdx;
    Soln[Number_of_Cells+1].phi = HyperHeat1D_U_ONE;

}

/********************************************************
 * Routine: Linear_Reconstruction_LeastSquares          *
 *                                                      *
 * Peforms the reconstruction of a  limited piecewise   *
 * linear solution state within each cell of the        *
 * computational mesh.  A least squares approach is     *
 * used in the evaluation of the unlimited solution     *
 * gradients.  Several slope limiters may be used.      *
 *                                                      *
 ********************************************************/
void Linear_Reconstruction_LeastSquares(HyperHeat1D_UniformMesh *Soln,
                                        const int Number_of_Cells,
                                        const int Limiter) {

    int i, n, n2, n_pts, index[2];
    double u0Min, u0Max, uQuad[2], phi;
    double Dx, DxDx_ave;
    HyperHeat1D_State DU, DUDx_ave;

    /* Carry out the limited solution reconstruction in
       each cell. */

    for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
        n_pts = 2;
        index[0] = i-1;
        index[1] = i+1; 

        DUDx_ave = HyperHeat1D_U_ZERO;
        DxDx_ave = ZERO;
    
        for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
            Dx = Soln[ index[n2] ].X.x - Soln[i].X.x;
            DU = Soln[ index[n2] ].U - Soln[i].U;
            DUDx_ave += DU*Dx;
            DxDx_ave += Dx*Dx;
        } /* endfor */
    					    
        DUDx_ave = DUDx_ave/double(n_pts);
        DxDx_ave = DxDx_ave/double(n_pts);

        Soln[i].dUdx = DUDx_ave/DxDx_ave;

        for ( n = 1 ; n <= NUM_VAR_HYPERHEAT1D ; ++n ) {
           u0Min = Soln[i].U[n];
           u0Max = u0Min;
           for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
              u0Min = min(u0Min, Soln[ index[n2] ].U[n]);
              u0Max = max(u0Max, Soln[ index[n2] ].U[n]);
           } /* endfor */

           uQuad[0] = Soln[i].U[n] - HALF*Soln[i].dUdx[n]*Soln[i].X.dx;
           uQuad[1] = Soln[i].U[n] + HALF*Soln[i].dUdx[n]*Soln[i].X.dx;

           switch(Limiter) {
             case LIMITER_ZERO :
	       phi = ZERO;
               break;
             case LIMITER_ONE :
	       phi = ONE;
               break;
             case LIMITER_BARTH_JESPERSEN :
               phi = Limiter_BarthJespersen(uQuad, Soln[i].U[n], u0Min, u0Max, 2);
               break;
             case LIMITER_VENKATAKRISHNAN :
               phi = Limiter_Venkatakrishnan(uQuad, Soln[i].U[n], u0Min, u0Max, 2);
               break;
             case LIMITER_VANLEER :
               phi = Limiter_VanLeer(uQuad, Soln[i].U[n], u0Min, u0Max, 2);
               break;
             case LIMITER_VANALBADA :
               phi = Limiter_VanAlbada(uQuad, Soln[i].U[n], u0Min, u0Max, 2);
               break;
	     default:
               phi = Limiter_BarthJespersen(uQuad, Soln[i].U[n], u0Min, u0Max, 2);
               break;
           } /* endswitch */

	   Soln[i].phi[n] = phi;
        } /* endfor */

    } /* endfor */

    Soln[0].dUdx = Soln[1].phi^Soln[1].dUdx;
    Soln[0].phi = HyperHeat1D_U_ONE;

    Soln[Number_of_Cells+1].dUdx = Soln[Number_of_Cells].phi^Soln[Number_of_Cells].dUdx;
    Soln[Number_of_Cells+1].phi = HyperHeat1D_U_ONE;

}

/********************************************************
 * Routine: dUdt_explicitEuler_upwind                   *
 *                                                      *
 * This routine updates the solution using a 1st-order  *
 * explicit Euler time integration and 1st-order upwind *
 * spatial discretization scheme in conjunction with    *
 * either the Godunov, Roe, Rusanov, or HLLE flux       *
 * functions.                                           *
 *                                                      *
 ********************************************************/
int dUdt_explicitEuler_upwind(HyperHeat1D_UniformMesh *Soln,
	                      const int Number_of_Cells,
                              const int BC_type_left,
			      const int BC_type_right,
                              const HyperHeat1D_State U_left,
                              const HyperHeat1D_State U_right,
			      double &dtMin,
			      const double &CFL_Number,
                              const int Flux_Function_Type,
			      const int Local_Time_Stepping) {

    int i;
    HyperHeat1D_State Flux;

    /* Determine local and global time steps and
       set the solution residual to zero. */
    
    dtMin = MILLION;

    for ( i = 0 ; i <= Number_of_Cells+1 ; ++i ) {
        Soln[i].dt = Soln[i].X.dx/Soln[i].U.a();
        dtMin = min(dtMin, Soln[i].dt);
        Soln[i].dUdt = HyperHeat1D_U_ZERO;
    } /* endfor */

    /* Evaluate the time rate of change of the solution
       (i.e., the solution residuals) using the first-order
       upwind scheme with a variety of flux functions. */

    for ( i = 0 ; i <= Number_of_Cells ; ++i ) {

        switch(Flux_Function_Type) {
          case FLUX_FUNCTION_GODUNOV :
            Flux = FluxGodunov(Soln[i].U, Soln[i+1].U);
            break;
          case FLUX_FUNCTION_ROE :
            Flux = FluxRoe(Soln[i].U, Soln[i+1].U);
            break;
          case FLUX_FUNCTION_RUSANOV :
            Flux = FluxRusanov(Soln[i].U, Soln[i+1].U);
            break;
          case FLUX_FUNCTION_HLLE :
            Flux = FluxHLLE(Soln[i].U, Soln[i+1].U);
            break;
	  default:
            Flux = FluxRoe(Soln[i].U, Soln[i+1].U);
            break;
        } /* endswitch */

        Soln[i].dUdt -= Flux/Soln[i].X.dx;
        Soln[i+1].dUdt += Flux/Soln[i+1].X.dx;

       /* Add the source terms. */

       if (i > 0) Soln[i].dUdt += Soln[i].U.S();

    } /* endfor */

    Soln[0].dUdt = HyperHeat1D_U_ZERO;
    Soln[Number_of_Cells+1].dUdt = HyperHeat1D_U_ZERO;

    /* Update solution variables using explicit Euler method. */

    for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
        if ( !Local_Time_Stepping ) Soln[i].dt = dtMin;
        Soln[i].U += (CFL_Number*Soln[i].dt)*Soln[i].dUdt;
    } /* endfor */

    /* Apply boundary conditions. */

    Soln[0].U = BCs(Soln[0].U, 
                    Soln[1].U, 
                    HyperHeat1D_U_ZERO,
                    Soln[1].X.dx,
                    BC_type_left, 
                    LEFT_END_BOUNDARY);

    Soln[Number_of_Cells+1].U = BCs(Soln[Number_of_Cells+1].U, 
                                    Soln[Number_of_Cells].U,
                                    HyperHeat1D_U_ZERO,
                                    Soln[Number_of_Cells].X.dx,
                                    BC_type_right, 
                                    RIGHT_END_BOUNDARY);

    /* Solution successfully updated. */

    return (0);
    
}

/********************************************************
 * Routine: dUdt_semiimplicitEuler_upwind               *
 *                                                      *
 * This routine updates the solution using a 1st-order  *
 * semi implicit Euler time integration and 1st-order   *
 * upwind spatial discretization scheme in conjunction  *
 * either the Godunov, Roe, Rusanov, or HLLE flux       *
 * functions.                                           *
 *                                                      *
 ********************************************************/
int dUdt_semiimplicitEuler_upwind(HyperHeat1D_UniformMesh *Soln,
	                          const int Number_of_Cells,
                                  const int BC_type_left,
				  const int BC_type_right,
                                  const HyperHeat1D_State U_left,
                                  const HyperHeat1D_State U_right,
			          double &dtMin,
				  const double &CFL_Number,
                                  const int Flux_Function_Type,
			          const int Local_Time_Stepping) {

    int i;
    HyperHeat1D_State Flux;

    /* Determine local and global time steps and
       set the solution residual to zero. */
    
    dtMin = MILLION;

    for ( i = 0 ; i <= Number_of_Cells+1 ; ++i ) {
        Soln[i].dt = Soln[i].X.dx/Soln[i].U.a();
        dtMin = min(dtMin, Soln[i].dt);
        Soln[i].dUdt = HyperHeat1D_U_ZERO;
    } /* endfor */

    /* Evaluate the time rate of change of the solution
       (i.e., the solution residuals) using the first-order
       upwind scheme with a variety of flux functions. */

    for ( i = 0 ; i <= Number_of_Cells ; ++i ) {

        switch(Flux_Function_Type) {
          case FLUX_FUNCTION_GODUNOV :
            Flux = FluxGodunov(Soln[i].U, Soln[i+1].U);
            break;
          case FLUX_FUNCTION_ROE :
            Flux = FluxRoe(Soln[i].U, Soln[i+1].U);
            break;
          case FLUX_FUNCTION_RUSANOV :
            Flux = FluxRusanov(Soln[i].U, Soln[i+1].U);
            break;
          case FLUX_FUNCTION_HLLE :
            Flux = FluxHLLE(Soln[i].U, Soln[i+1].U);
            break;
	  default:
            Flux = FluxRoe(Soln[i].U, Soln[i+1].U);
            break;
        } /* endswitch */

        Soln[i].dUdt -= Flux/Soln[i].X.dx;
        Soln[i+1].dUdt += Flux/Soln[i+1].X.dx;

        /* Add the source terms. */

        if (i > 0) Soln[i].dUdt += Soln[i].U.S();

    } /* endfor */

    Soln[0].dUdt = HyperHeat1D_U_ZERO;
    Soln[Number_of_Cells+1].dUdt = HyperHeat1D_U_ZERO;

    /* Update solution variables using semi-implicit Euler method. */

    for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
        if ( !Local_Time_Stepping ) Soln[i].dt = dtMin;

        Soln[i].U.T += (CFL_Number*Soln[i].dt)*Soln[i].dUdt.T;
        Soln[i].U.qx += (CFL_Number*Soln[i].dt)*Soln[i].dUdt.qx/
                        (ONE+CFL_Number*Soln[i].dt/Soln[i].U.tau);
    } /* endfor */

    /* Apply boundary conditions. */

    Soln[0].U = BCs(Soln[0].U, 
                    Soln[1].U, 
                    HyperHeat1D_U_ZERO,
                    Soln[1].X.dx,
                    BC_type_left, 
                    LEFT_END_BOUNDARY);

    Soln[Number_of_Cells+1].U = BCs(Soln[Number_of_Cells+1].U, 
                                    Soln[Number_of_Cells].U,
                                    HyperHeat1D_U_ZERO,
                                    Soln[Number_of_Cells].X.dx,
                                    BC_type_right, 
                                    RIGHT_END_BOUNDARY);

    /* Solution successfully updated. */

    return (0);
    
}

/********************************************************
 * Routine: dUdt_2stage_2ndOrder_upwind                 *
 *                                                      *
 * This routine updates the solution using a two-stage  *
 * second-order semi-implicit time integration scheme   *
 * and a 2nd-ororder limited upwind spatial             *
 * discretization scheme with either the Godunov, Roe,  *
 * Rusanov, or HLLE flux functions.                     *
 *                                                      *
 ********************************************************/
int dUdt_2stage_2ndOrder_upwind(HyperHeat1D_UniformMesh *Soln,
	                        const int Number_of_Cells,
                                const int BC_type_left,
				const int BC_type_right,
                                const HyperHeat1D_State U_left,
                                const HyperHeat1D_State U_right,
			        double &dtMin,
				const double &CFL_Number,
                                const int Reconstruction_Type,
                                const int Limiter_Type,
                                const int Flux_Function_Type,
			        const int Local_Time_Stepping) {

    int i, n_stage;
    double omega;
    HyperHeat1D_State Ul, Ur, Flux;

    /* Perform second-order two-stage semi-implicit update of solution
       varibles for new time level. */

    for ( n_stage = 1 ; n_stage <= 2 ; ++n_stage ) {

        /* Evaluate the time step fraction for stage. */

        omega = ONE/double(n_stage);

        /* Apply boundary conditions for stage. */

        Soln[0].U = BCs(U_left,
                        Soln[1].U,
                        (Soln[2].U-Soln[1].U)/Soln[1].X.dx,
                        Soln[1].X.dx,
                        BC_type_left, 
                        LEFT_END_BOUNDARY);

        Soln[Number_of_Cells+1].U = BCs(U_right,
                                        Soln[Number_of_Cells].U,
                                        (Soln[Number_of_Cells].U-Soln[Number_of_Cells-1].U)/
                                        Soln[Number_of_Cells].X.dx, 
                                        Soln[Number_of_Cells].X.dx,
                                        BC_type_right, 
                                        RIGHT_END_BOUNDARY);

        /* Perform the linear reconstruction within each cell
           of the computational grid in each stage of update. */

        if ( n_stage == 1 ) {
          switch(Reconstruction_Type) {
            case RECONSTRUCTION_MUSCL :
               Linear_Reconstruction_MUSCL(Soln, 
                                           Number_of_Cells, 
                                           Limiter_Type);
              break;
            case RECONSTRUCTION_GREEN_GAUSS :
               Linear_Reconstruction_GreenGauss(Soln, 
                                                Number_of_Cells, 
                                                Limiter_Type);
              break;
            case RECONSTRUCTION_LEAST_SQUARES :
               Linear_Reconstruction_LeastSquares(Soln, 
                                                  Number_of_Cells, 
                                                  Limiter_Type);
              break;
	    default:
               Linear_Reconstruction_MUSCL(Soln, 
                                           Number_of_Cells, 
                                           Limiter_Type);
              break;
          } /* endswitch */
        } /* endif */

        /* Evaluate the time rate of change of the solution
           (i.e., the solution residuals) using a second-order
           limited upwind scheme with a variety of flux functions. */

        if ( !Local_Time_Stepping && n_stage == 1 ) Soln[0].dt = dtMin;
        if ( n_stage == 1 ) Soln[0].Uo = Soln[0].U;
        Soln[0].dUdt = HyperHeat1D_U_ZERO;

        for ( i = 0 ; i <= Number_of_Cells ; ++i ) {
            if ( !Local_Time_Stepping && n_stage == 1 ) Soln[i+1].dt = dtMin;
            if ( n_stage == 1 ) {
               Soln[i+1].Uo = Soln[i+1].U;
               Soln[i+1].dUdt = HyperHeat1D_U_ZERO;
            } else {
               Soln[i+1].dUdt = HALF*Soln[i+1].dUdt;
            } /* endif */

            /* Evaluate the cell interface flux. */
            Ul = Soln[i].U + HALF*(Soln[i].phi^Soln[i].dUdx)*Soln[i].X.dx;
            Ur = Soln[i+1].U - HALF*(Soln[i+1].phi^Soln[i+1].dUdx)*Soln[i+1].X.dx;

            switch(Flux_Function_Type) {
              case FLUX_FUNCTION_GODUNOV :
                Flux = FluxGodunov(Ul, Ur);
                break;
              case FLUX_FUNCTION_ROE :
                Flux = FluxRoe(Ul, Ur);
                break;
              case FLUX_FUNCTION_RUSANOV :
                Flux = FluxRusanov(Ul, Ur);
                break;
              case FLUX_FUNCTION_HLLE :
                Flux = FluxHLLE(Ul, Ur);
                break;
	      default:
                Flux = FluxRoe(Ul, Ur);
                break;
            } /* endswitch */

            Soln[i].dUdt -= (omega*CFL_Number*Soln[i].dt)*Flux/Soln[i].X.dx;
            Soln[i+1].dUdt += (omega*CFL_Number*Soln[i+1].dt)*Flux/Soln[i+1].X.dx;

            /* Add the source terms. */
            if (i > 0) Soln[i].dUdt += (omega*CFL_Number*Soln[i].dt)*Soln[i].Uo.S();

        } /* endfor */

        Soln[0].dUdt = HyperHeat1D_U_ZERO;
        Soln[Number_of_Cells+1].dUdt = HyperHeat1D_U_ZERO;

        /* Update solution variables for stage. */

        for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
            Soln[i].U.T = Soln[i].Uo.T + Soln[i].dUdt.T;
            Soln[i].U.qx = Soln[i].Uo.qx + Soln[i].dUdt.qx/
                           (ONE+(omega*CFL_Number*Soln[i].dt)/Soln[i].U.tau);
        } /* endfor */

    } /* endfor */

    /* Solution successfully updated. */

    return (0);
    
}

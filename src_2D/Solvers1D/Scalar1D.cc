/* Scalar1D.cc:  Subroutines for 1D Scalar Advection Solution Classes. */

/* Include 1D Scalar advection solution header file. */

#ifndef _SCALAR1D_INCLUDED
#include "Scalar1D.h"
#endif // _SCALAR1D_INCLUDED

/********************************************************
 * Routine: Allocate                                    *
 *                                                      *
 * Allocate memory for 1D scalar advection equation     *
 * solution.                                            *
 *                                                      *
 ********************************************************/
Scalar1D_UniformMesh* Allocate(Scalar1D_UniformMesh *Soln_ptr,
                               const int Number_of_Cells) {

    /* Allocate memory. */

    Soln_ptr = new Scalar1D_UniformMesh[Number_of_Cells+2];

    /* Return memory location. */

    return(Soln_ptr);

}

/********************************************************
 * Routine: Deallocate                                  *
 *                                                      *
 * Deallocate memory for 1D scalar advection equation   *
 * solution.                                            *
 *                                                      *
 ********************************************************/
Scalar1D_UniformMesh* Deallocate(Scalar1D_UniformMesh *Soln_ptr,
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
void Output_Gnuplot(Scalar1D_UniformMesh *Soln,
                    const int Number_of_Cells,
                    const int Number_of_Time_Steps,
                    const double &Time,
	            ostream &out_file) {

    int i;

    out_file << "# " << CFDkit_Name() << ": 1D Scalar Advection Solution, "
             << "Time Step/Iteration Level = "
             << Number_of_Time_Steps
             << ", Time = " << Time << "\n"
	     << "# cell, x, dx, u, a\n";
    
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
void Output_Tecplot(Scalar1D_UniformMesh *Soln,
                    const int Number_of_Cells,
                    const int Number_of_Time_Steps,
                    const double &Time,
	            ostream &out_file) {

    int i;

    out_file << "TITLE = \"" << CFDkit_Name() 
             << ": 1D Scalar Advection Solution, "
             << "Time Step/Iteration Level = "
             << Number_of_Time_Steps
             << ", Time = " << Time << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
             << "\"dx\" \\ \n"
             << "\"u\" \\ \n"
             << "\"a\" \n"
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
void Grid(Scalar1D_UniformMesh *Soln,
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
void ICs(Scalar1D_UniformMesh *Soln,
	 const int i_ICtype,
         const int Number_of_Cells) {

    int i;
    double ul,ur,um,umin,umax,xl,xr,c;

    /* Assign the initial data for the IVP of interest. */
    
    switch(i_ICtype) {
      case IC_CONSTANT :
        ul = ONE;
        c = ONE;
        for ( i = 0 ; i <= Number_of_Cells+1 ; ++i ) {
            Soln[i].u = ul;
            Soln[i].a = c;
        } /* endfor */
        break;
      case IC_UNIFORM :
        ul = ONE;
        c = ONE;
        for ( i = 0 ; i <= Number_of_Cells+1 ; ++i ) {
            Soln[i].u = ul;
            Soln[i].a = c;
        } /* endfor */
        break;
      case IC_RIEMANN :
        ul = TWO;
        ur = ONE;
        c = ONE;
        for ( i = 0 ; i <= Number_of_Cells+1 ; ++i ) {
	    if (Soln[i].X.x <= ZERO) {
               Soln[i].u = ul;
            } else {
               Soln[i].u = ur;	     
            } /* end if */
            Soln[i].a = c; 
        } /* endfor */
        break;
      case IC_SQUARE_WAVE :
        ul = ONE;
        um = THREE;
        ur = ONE;
        c = ONE;
        for ( i = 0 ; i <= Number_of_Cells+1 ; ++i ) {
	    if (Soln[i].X.x <= -0.80) {
               Soln[i].u = ul;
	    } else if (Soln[i].X.x >= -0.50 ) {
	       Soln[i].u = ur;
            } else {
               Soln[i].u = um;
            } /* end if */
            Soln[i].a = c;
        } /* endfor */
        break;
      case IC_SINX2_WAVE :
        umax = ONE;
        umin = ZERO;
        xl = -0.80;
        xr = -0.40;
        c = ONE;
        for ( i = 0 ; i <= Number_of_Cells+1 ; ++i ) {
	    if (Soln[i].X.x <= xl) {
               Soln[i].u = umin;
	    } else if (Soln[i].X.x >= xr ) {
	       Soln[i].u = umin;
            } else {
               Soln[i].u = umin + 
                           (umax-umin)*sqr(sin((Soln[i].X.x-xl)*PI/(xr-xl)));
            } /* end if */
            Soln[i].a = c;
        } /* endfor */
        break;
      default:
        ul = ONE;
        c = ONE;
        for ( i = 0 ; i <= Number_of_Cells+1 ; ++i ) {
            Soln[i].u = ul;
            Soln[i].a = c;
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
double CFL(Scalar1D_UniformMesh *Soln,
           const int Number_of_Cells) {

    int i;
    double dtMin;

    /* Determine local and global time steps. */

    dtMin = MILLION;

    for ( i = 0 ; i <= Number_of_Cells+1 ; ++i ) {
        Soln[i].dt = Soln[i].X.dx/Soln[i].a;
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
void Linear_Reconstruction_MUSCL(Scalar1D_UniformMesh *Soln,
                                 const int Number_of_Cells,
                                 const int Limiter) {

    int i;
    double a, b;

    /* Carry out the limited solution reconstruction in
       each cell. */

    b = (Soln[1].u-Soln[0].u)/Soln[0].X.dx;
    for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
        Soln[i].dudx = HALF*(Soln[i+1].u-Soln[i-1].u)/
                       Soln[i].X.dx;
        a = b;
        b = (Soln[i+1].u-Soln[i].u)/Soln[i].X.dx;

        switch(Limiter) {
          case LIMITER_ZERO :
	    Soln[i].phi = ZERO;
            break;
          case LIMITER_ONE :
	    Soln[i].phi = ONE;
            break;
          case LIMITER_MINMOD :
	    Soln[i].phi = minmod(a, b)/(Soln[i].dudx+sgn(Soln[i].dudx)*TOLER*TOLER);
            break; 
          case LIMITER_UMIST :
            Soln[i].phi = minmod(TWO*a, 
                         TWO*b,
                         (a+THREE*b)/FOUR,
                         (THREE*a+b)/FOUR)/
                  (Soln[i].dudx+sgn(Soln[i].dudx)*TOLER*TOLER);
            break;
          case LIMITER_DOUBLE_MINMOD :
            Soln[i].phi = minmod(HALF*(a+b), 
	                 TWO*a, 
	                 TWO*b)/(Soln[i].dudx+sgn(Soln[i].dudx)*TOLER*TOLER);
            break;
          case LIMITER_SUPERBEE :
            Soln[i].phi = superbee(a, b)/
                  (Soln[i].dudx+sgn(Soln[i].dudx)*TOLER*TOLER);
            break; 
          case LIMITER_PHI :
            Soln[i].phi = philimiter(a, b, 1.50)/
                  (Soln[i].dudx+sgn(Soln[i].dudx)*TOLER*TOLER);
            break;
          case LIMITER_VANLEER :
            Soln[i].phi = vanleer(a, b)/
                  (Soln[i].dudx+sgn(Soln[i].dudx)*TOLER*TOLER);
            break;
          case LIMITER_VANALBADA :
            Soln[i].phi = vanalbada(a, b, 0.10)/
                  (Soln[i].dudx+sgn(Soln[i].dudx)*TOLER*TOLER);
            break;
	  default:
            Soln[i].phi = philimiter(a, b, ONE)/
                  (Soln[i].dudx+sgn(Soln[i].dudx)*TOLER*TOLER);
            break;
        } /* endswitch */
    } /* endfor */

    Soln[0].dudx = Soln[1].phi*Soln[1].dudx;
    Soln[0].phi = ONE;

    Soln[Number_of_Cells+1].dudx = 
       Soln[Number_of_Cells].phi*Soln[Number_of_Cells].dudx;
    Soln[Number_of_Cells+1].phi = ONE;

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
void Linear_Reconstruction_GreenGauss(Scalar1D_UniformMesh *Soln,
                                      const int Number_of_Cells,
                                      const int Limiter) {

    int i;
    double u0Min, u0Max, uQuad[2];

    /* Carry out the limited solution reconstruction in
       each cell. */

    for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
        Soln[i].dudx = HALF*(Soln[i+1].u-Soln[i-1].u)/
                       Soln[i].X.dx;

        u0Min = min(Soln[i-1].u, Soln[i].u);
        u0Min = min(u0Min, Soln[i+1].u);
        u0Max = max(Soln[i-1].u, Soln[i].u);
        u0Max = max(u0Max, Soln[i+1].u);
        uQuad[0] = Soln[i].u - HALF*Soln[i].dudx*Soln[i].X.dx;
        uQuad[1] = Soln[i].u + HALF*Soln[i].dudx*Soln[i].X.dx;
	
        switch(Limiter) {
          case LIMITER_ZERO :
	    Soln[i].phi = ZERO;
            break;
          case LIMITER_ONE :
	    Soln[i].phi = ONE;
            break;
          case LIMITER_BARTH_JESPERSEN :
            Soln[i].phi = Limiter_BarthJespersen(uQuad, Soln[i].u, u0Min, u0Max, 2);
            break;
          case LIMITER_VENKATAKRISHNAN :
            Soln[i].phi = Limiter_Venkatakrishnan(uQuad, Soln[i].u, u0Min, u0Max, 2);
            break;
          case LIMITER_VANLEER :
            Soln[i].phi = Limiter_VanLeer(uQuad, Soln[i].u, u0Min, u0Max, 2);
            break;
          case LIMITER_VANALBADA :
            Soln[i].phi = Limiter_VanAlbada(uQuad, Soln[i].u, u0Min, u0Max, 2);
            break;
	  default:
            Soln[i].phi = Limiter_BarthJespersen(uQuad, Soln[i].u, u0Min, u0Max, 2);
            break;
        } /* endswitch */
    } /* endfor */

    Soln[0].dudx = Soln[1].phi*Soln[1].dudx;
    Soln[0].phi = ONE;

    Soln[Number_of_Cells+1].dudx = 
       Soln[Number_of_Cells].phi*Soln[Number_of_Cells].dudx;
    Soln[Number_of_Cells+1].phi = ONE;

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
void Linear_Reconstruction_LeastSquares(Scalar1D_UniformMesh *Soln,
                                        const int Number_of_Cells,
                                        const int Limiter) {

    int i, n2, n_pts, index[2];
    double u0Min, u0Max, uQuad[2];
    double Dx, DxDx_ave, DU, DUDx_ave;

    /* Carry out the limited solution reconstruction in
       each cell. */

    for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
        n_pts = 2;
        index[0] = i-1;
        index[1] = i+1; 

        DUDx_ave = ZERO;
        DxDx_ave = ZERO;
    
        for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
            Dx = Soln[ index[n2] ].X.x - Soln[i].X.x;
            DU = Soln[ index[n2] ].u - Soln[i].u;
            DUDx_ave += DU*Dx;
            DxDx_ave += Dx*Dx;
        } /* endfor */
    					    
        DUDx_ave = DUDx_ave/double(n_pts);
        DxDx_ave = DxDx_ave/double(n_pts);

        Soln[i].dudx = DUDx_ave/DxDx_ave;

        u0Min = Soln[i].u;
        u0Max = u0Min;
        for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
           u0Min = min(u0Min, Soln[ index[n2] ].u);
           u0Max = max(u0Max, Soln[ index[n2] ].u);
        } /* endfor */
	
        uQuad[0] = Soln[i].u - HALF*Soln[i].dudx*Soln[i].X.dx;
        uQuad[1] = Soln[i].u + HALF*Soln[i].dudx*Soln[i].X.dx;
	
        switch(Limiter) {
          case LIMITER_ZERO :
	    Soln[i].phi = ZERO;
            break;
          case LIMITER_ONE :
	    Soln[i].phi = ONE;
            break;
          case LIMITER_BARTH_JESPERSEN :
            Soln[i].phi = Limiter_BarthJespersen(uQuad, Soln[i].u, u0Min, u0Max, 2);
            break;
          case LIMITER_VENKATAKRISHNAN :
            Soln[i].phi = Limiter_Venkatakrishnan(uQuad, Soln[i].u, u0Min, u0Max, 2);
            break;
          case LIMITER_VANLEER :
            Soln[i].phi = Limiter_VanLeer(uQuad, Soln[i].u, u0Min, u0Max, 2);
            break;
          case LIMITER_VANALBADA :
            Soln[i].phi = Limiter_VanAlbada(uQuad, Soln[i].u, u0Min, u0Max, 2);
            break;
	  default:
            Soln[i].phi = Limiter_BarthJespersen(uQuad, Soln[i].u, u0Min, u0Max, 2);
            break;
        } /* endswitch */
    } /* endfor */

    Soln[0].dudx = Soln[1].phi*Soln[1].dudx;
    Soln[0].phi = ONE;

    Soln[Number_of_Cells+1].dudx = 
       Soln[Number_of_Cells].phi*Soln[Number_of_Cells].dudx;
    Soln[Number_of_Cells+1].phi = ONE;

}

/********************************************************
 * Routine: dudt_explicitEuler_upwind                   *
 *                                                      *
 * This routine updates the solution using a 1st-order  *
 * explicit Euler time integration and 1st-order upwind *
 * spatial discretization scheme.                       *
 *                                                      *
 ********************************************************/
int dudt_explicitEuler_upwind(Scalar1D_UniformMesh *Soln,
	                      const int Number_of_Cells,
			      double &dtMin,
			      const double &CFL_Number,
			      const int Local_Time_Stepping) {

    int i;
    double aface, flux;

    /* Evaluate the time rate of change of the solution
       (i.e., the solution residuals) using the first-order
       upwind scheme and Godunov flux function. */

    Soln[0].dudt = ZERO;
    for ( i = 0 ; i <= Number_of_Cells ; ++i ) {
        Soln[i+1].dudt = ZERO;
        aface = HALF*(Soln[i].a+Soln[i+1].a);
	if (aface >= ZERO) {
	    flux = aface*Soln[i].u;
	} else {
	    flux = aface*Soln[i+1].u;
	} /* endif */

        Soln[i].dudt -= flux/Soln[i].X.dx;
        Soln[i+1].dudt += flux/Soln[i+1].X.dx;
    } /* endfor */
    Soln[Number_of_Cells+1].dudt = ZERO;

    /* Now update the solution. */

    for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
        if ( !Local_Time_Stepping ) Soln[i].dt = dtMin;
	
        Soln[i].u += (CFL_Number*Soln[i].dt)*Soln[i].dudt;
    } /* endfor */

    /* By default, constant extrapolation boundary
       conditions are applied at either end of the mesh. */

    Soln[0].u = Soln[1].u;
    
    Soln[Number_of_Cells+1].u = Soln[Number_of_Cells].u;

    /* Solution successfully updated. */

    return (0);
    
}

/********************************************************
 * Routine: dudt_LaxFriedrichs                          *
 *                                                      *
 * This routine updates the solution using the          *
 * first-order single-stage Lax-Friedrichs scheme.      *
 * See Lax (1954).                                      *
 *                                                      *
 ********************************************************/
int dudt_LaxFriedrichs(Scalar1D_UniformMesh *Soln,
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

    Soln[0].dudt = ZERO;
    for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
        if ( !Local_Time_Stepping ) Soln[i].dt = dtMin;
	
        Soln[i].dudt = - Soln[i].a*
	                 (Soln[i+1].u-Soln[i-1].u)/
                         (TWO*delta_x)
	               + HALF*(Soln[i+1].u-TWO*Soln[i].u+
		   	       Soln[i-1].u)/
                              (CFL_Number*Soln[i].dt);
    } /* endfor */
    Soln[Number_of_Cells+1].dudt = ZERO;

    /* Now update the solution. */

    for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
        Soln[i].u += (CFL_Number*Soln[i].dt)*Soln[i].dudt;
    } /* endfor */

    /* By default, constant extrapolation boundary
       conditions are applied at either end of the mesh. */

    Soln[0].u = Soln[1].u;

    Soln[Number_of_Cells+1].u = Soln[Number_of_Cells].u;

    /* Solution successfully updated. */

    return (0);
    
}

/********************************************************
 * Routine: dudt_LaxWendroff                            *
 *                                                      *
 * This routine updates the solution using the          *
 * second-order single-step Lax-Wendroff scheme.        *
 * See Lax and Wendroff (1960, 1964).                   *
 *                                                      *
 ********************************************************/
int dudt_LaxWendroff(Scalar1D_UniformMesh *Soln,
	             const int Number_of_Cells,
		     double &dtMin,
		     const double &CFL_Number,
		     const int Local_Time_Stepping) {

    int i;
    double delta_x;
    
    /* Assign the mesh space to local variable. */

    delta_x = Soln[0].X.dx;

    /* Evaluate the time rate of change of the solution
       (i.e., the solution residuals) using the second-order
       Lax-Friedrichs scheme. */

    Soln[0].dudt = ZERO;
    for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
        if ( !Local_Time_Stepping ) Soln[i].dt = dtMin;
	
        Soln[i].dudt = - Soln[i].a*
	                 (Soln[i+1].u-Soln[i-1].u)/
                         (TWO*delta_x)
	               + HALF*sqr(Soln[i].a)*
                         (CFL_Number*Soln[i].dt)*
                         (Soln[i+1].u-TWO*Soln[i].u+
		   	  Soln[i-1].u)/sqr(delta_x);
    } /* endfor */
    Soln[Number_of_Cells+1].dudt = ZERO;

    /* Now update the solution. */

    for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
        Soln[i].u += (CFL_Number*Soln[i].dt)*Soln[i].dudt;
    } /* endfor */

    /* By default, constant extrapolation boundary
       conditions are applied at either end of the mesh. */

    Soln[0].u = Soln[1].u;

    Soln[Number_of_Cells+1].u = Soln[Number_of_Cells].u;

    /* Solution successfully updated. */

    return (0);
    
}

/********************************************************
 * Routine: dudt_MacCormack                             *
 *                                                      *
 * This routine updates the solution using the          *
 * second-order predictor-corrector scheme of           *
 * MacCormack (1969, 1971).  Artificial viscosity of    *
 * the type devised by Harten and Zwas (1972) is        *
 * added to smooth the solution near discontinuities.   *
 *                                                      *
 ********************************************************/
int dudt_MacCormack(Scalar1D_UniformMesh *Soln,
	            const int Number_of_Cells,
		    double &dtMin,
		    const double &CFL_Number,
		    const int Local_Time_Stepping) {

    int i;
    static int Order=1;
    double ul, ur;
    double delta_x, delta_l, delta_r, 
           theta, theta_r, theta_l,
           switch_l, switch_r;

    /* Detemine order of first-order spatial operators. */
    
    if (Order) {
       Order = 0;
    } else {
       Order = 1;
    } /* endif */

    /* Assign the mesh space to local variable. */

    delta_x = Soln[0].X.dx;

    /* Evaluate the time rate of change of the solution
       (i.e., the solution residuals) for the predictor
       step using first-order backward/forward spatial
       differences. */

    Soln[0].dudt = ZERO;
    for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
        Soln[i].dudt = - Soln[i].a*
                         (Soln[i+1-Order].u - Soln[i-Order].u)/
                         delta_x;
    } /* endfor */
    Soln[Number_of_Cells+1].dudt = ZERO;

    /* Evaluate the time rate of change of the solution
       (i.e., the solution residuals) for the corrector
       step using first-order backward/forward spatial
       differences. */

    if ( !Local_Time_Stepping ) Soln[Order].dt = dtMin;
    ur = Soln[Order].u +
         (CFL_Number*Soln[Order].dt)*Soln[Order].dudt;

    for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
        if ( !Local_Time_Stepping ) Soln[i+Order].dt = dtMin;

	/* Evaluate temporary predictor solution states. */
        ul = ur;
	ur = Soln[i+Order].u +
	     (CFL_Number*Soln[i+Order].dt)*Soln[i+Order].dudt;

	/* Evaluate corrector solution residual. */
	Soln[i].dudt =   HALF*Soln[i].dudt
	               - HALF*Soln[i].a*(ur - ul)/delta_x;

	/* Add artificial dissipation of Harten and Zwas. */
	if (i == 1) {
            theta_l = ZERO;
	    
	    delta_l = fabs(Soln[i].u-Soln[i-1].u)/Soln[i].u;
	    delta_r = fabs(Soln[i+1].u-Soln[i].u)/Soln[i].u;
            if (delta_l + delta_r < TOLER) {
	        theta = ZERO;
	    } else {
	        theta = fabs((delta_r - delta_l)/(delta_l + delta_r));
            } /* endif */
	    
	    delta_l = fabs(Soln[i+1].u-Soln[i].u)/Soln[i+1].u;
	    delta_r = fabs(Soln[i+2].u-Soln[i+1].u)/Soln[i+1].u;
            if (delta_l + delta_r < TOLER) {
	        theta_r = ZERO;
	    } else {
	        theta_r = fabs((delta_r - delta_l)/(delta_l + delta_r));
            } /* endif */
	} else if (i < Number_of_Cells) {
	    theta_l = theta;
	    theta = theta_r;
	    
	    delta_l = fabs(Soln[i+1].u-Soln[i].u)/Soln[i+1].u;
	    delta_r = fabs(Soln[i+2].u-Soln[i+1].u)/Soln[i+1].u;
            if (delta_l + delta_r < TOLER) {
	        theta_r = ZERO;
	    } else {
	        theta_r = fabs((delta_r - delta_l)/(delta_l + delta_r));
            } /* endif */
        } else {
	    theta_l = theta;
	    theta = theta_r;
            theta_r = ZERO;
        } /* endif */
        switch_r = max(theta,theta_r);
        switch_l = max(theta,theta_l);
	
	Soln[i].dudt += (switch_r*(Soln[i+1].u-Soln[i].u) -
		   	 switch_l*(Soln[i].u-Soln[i-1].u))/
	                (EIGHT*CFL_Number*Soln[i].dt);
    } /* endfor */

    /* Now update the solution. */

    for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
	Soln[i].u += (CFL_Number*Soln[i].dt)*Soln[i].dudt;
    } /* endfor */

    /* By default, constant extrapolation boundary
       conditions are applied at either end of the mesh. */

    Soln[0].u = Soln[1].u;

    Soln[Number_of_Cells+1].u = Soln[Number_of_Cells].u;

    /* Solution successfully updated. */

    return (0);
    
}

/********************************************************
 * Routine: dUdt_2stage_2ndOrder_upwind                 *
 *                                                      *
 * This routine updates the solution using a two-stage  *
 * second-order explicit time integration scheme        *
 * and a 2nd-ororder limited upwind spatial             *
 * discretization scheme.                               *
 *                                                      *
 ********************************************************/
int dUdt_2stage_2ndOrder_upwind(Scalar1D_UniformMesh *Soln,
	                        const int Number_of_Cells,
			        double &dtMin,
				const double &CFL_Number,
                                const int Reconstruction_Type,
                                const int Limiter_Type,
			        const int Local_Time_Stepping) {

    int i, n_stage;
    double omega, ul, ur, aface, flux;

    /* Perform second-order two-stage semi-implicit update of solution
       varibles for new time level. */

    for ( n_stage = 1 ; n_stage <= 2 ; ++n_stage ) {

        /* Evaluate the time step fraction for stage. */

        omega = ONE/double(n_stage);

        /* Apply boundary conditions for stage. */


        /* Perform the linear reconstruction within each cell
           of the computational grid in first stage. */

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
        if ( n_stage == 1 ) Soln[0].uo = Soln[0].u;
        Soln[0].dudt = ZERO;

        for ( i = 0 ; i <= Number_of_Cells ; ++i ) {
            if ( !Local_Time_Stepping && n_stage == 1 ) Soln[i+1].dt = dtMin;
            if ( n_stage == 1 ) {
               Soln[i+1].uo = Soln[i+1].u;
               Soln[i+1].dudt = ZERO;
            } else {
               Soln[i+1].dudt = HALF*Soln[i+1].dudt;
            } /* endif */

            /* Evaluate the cell interface flux. */

            ul = Soln[i].u + HALF*(Soln[i].phi*Soln[i].dudx)*Soln[i].X.dx;
            ur = Soln[i+1].u - HALF*(Soln[i+1].phi*Soln[i+1].dudx)*Soln[i+1].X.dx;
            aface = HALF*(Soln[i].a+Soln[i+1].a);
	    if (aface >= ZERO) {
	       flux = aface*ul;
	    } else {
	       flux = aface*ur;
	    } /* endif */

            /* Evaluate cell-averaged solution changes. */

            Soln[i].dudt -= (omega*CFL_Number*Soln[i].dt)*flux/Soln[i].X.dx;
            Soln[i+1].dudt += (omega*CFL_Number*Soln[i+1].dt)*flux/Soln[i+1].X.dx;

        } /* endfor */

        Soln[0].dudt = ZERO;
        Soln[Number_of_Cells+1].dudt = ZERO;

        /* Update solution variables for this stage. */

        for ( i = 1 ; i <= Number_of_Cells ; ++i ) {
            Soln[i].u = Soln[i].uo + Soln[i].dudt;
        } /* endfor */

    } /* endfor */

    /* By default, constant extrapolation boundary
       conditions are applied at either end of the mesh. */

    Soln[0].u = Soln[1].u;
    
    Soln[Number_of_Cells+1].u = Soln[Number_of_Cells].u;

    /* Solution successfully updated. */

    return (0);
    
}

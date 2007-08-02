/* Heat1D.cc:  Subroutines for 1D Heat Equation Solution Classes. */

/* Include 1D heat equation solution header file. */

#ifndef _HEAT1D_INCLUDED
#include "Heat1D.h"
#endif // _HEAT1D_INCLUDED

/********************************************************
 * Routine: Allocate                                    *
 *                                                      *
 * Allocate memory for 1D heat equation solution.       *
 *                                                      *
 ********************************************************/
Heat1D_UniformMesh* Allocate(Heat1D_UniformMesh *Soln_ptr,
                             const int Number_of_Nodes) {

    /* Allocate memory. */

    Soln_ptr = new Heat1D_UniformMesh[Number_of_Nodes];

    /* Return memory location. */

    return(Soln_ptr);

}

/********************************************************
 * Routine: Deallocate                                  *
 *                                                      *
 * Deallocate memory for 1D heat equation solution.     *
 *                                                      *
 ********************************************************/
Heat1D_UniformMesh* Deallocate(Heat1D_UniformMesh *Soln_ptr,
                               const int Number_of_Nodes) {

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
void Output_Gnuplot(Heat1D_UniformMesh *Soln,
                    const int Number_of_Nodes,
                    const int Number_of_Time_Steps,
                    const double &Time,
	            ostream &out_file) {

    int i;

    out_file << "# " << CFDkit_Name() << ": 1D Heat Equation Solution, "
             << "Time Step/Iteration Level = "
             << Number_of_Time_Steps
             << ", Time = " << Time << "\n"
	     << "# cell, x, dx, T, qx, kappa\n";
    
    for ( i = 0 ; i <= Number_of_Nodes-1 ; ++i ) {
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
void Output_Tecplot(Heat1D_UniformMesh *Soln,
                    const int Number_of_Nodes,
                    const int Number_of_Time_Steps,
                    const double &Time,
	            ostream &out_file) {

    int i;

    out_file << "TITLE = \"" << CFDkit_Name()
             << ": 1D Heat Equation Solution, "
             << "Time Step/Iteration Level = "
             << Number_of_Time_Steps
             << ", Time = " << Time << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
             << "\"dx\" \\ \n"
             << "\"T\" \\ \n"
             << "\"qx\" \\ \n"
             << "\"kappa\" \n"
             << "ZONE \n";
    
    for ( i = 0 ; i <= Number_of_Nodes-1 ; ++i ) {
        out_file << " " << Soln[i] << "\n";
    } /* endfor */

    out_file << "\n";
    
}

/********************************************************
 * Routine: Grid                                        *
 *                                                      *
 * Generates a uniform mesh and assign the locations of *
 * the nodes to appropriate solution variables.         *
 *                                                      *
 ********************************************************/
void Grid(Heat1D_UniformMesh *Soln,
          const double &xMin,
	  const double &xMax,
	  const int Number_of_Nodes) {

    int i;
    double delta_x;

    /* Determine the mesh spacing. */

    delta_x = (xMax - xMin)/double(Number_of_Nodes-1);
    Soln[0].X.setsize(delta_x);

    /* Create the nodes. */

    Soln[0].X.x = xMin;

    for ( i = 1 ; i <= Number_of_Nodes-1 ; ++i ) {
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
void ICs(Heat1D_UniformMesh *Soln,
         const int i_ICtype,
         const double &Kappa,
         const int Number_of_Nodes) {

    int i, n_frequency;
    double T0, rod_length, rod_kappa, amplitude;

    /* Assign the initial data for the IVP of interest. */
    
    switch(i_ICtype) {
      case IC_CONSTANT :
        T0 = TEN;
        rod_kappa = Kappa;
        rod_length = Soln[Number_of_Nodes-1].X.x - Soln[0].X.x;
        for ( i = 0 ; i <= Number_of_Nodes-1 ; ++i ) {
            Soln[i].T = T0;
            Soln[i].kappa = rod_kappa;
        } /* endfor */
        break;
      case IC_UNIFORM :
        T0 = TEN;
        rod_kappa = Kappa;
        rod_length = Soln[Number_of_Nodes-1].X.x - Soln[0].X.x;
        for ( i = 0 ; i <= Number_of_Nodes-1 ; ++i ) {
            Soln[i].T = T0;
            Soln[i].kappa = rod_kappa;
        } /* endfor */
        break;
      case IC_IMPULSIVE_ROD :
        T0 = ZERO;
        rod_kappa = Kappa;
        rod_length = Soln[Number_of_Nodes-1].X.x - Soln[0].X.x;
        for ( i = 0 ; i <= Number_of_Nodes-1 ; ++i ) {
            Soln[i].T = T0;
            Soln[i].kappa = rod_kappa;
        } /* endfor */
        T0 = HUNDRED;
        Soln[Number_of_Nodes-1].T = T0;
        break;
      case IC_SINUSOIDAL_ROD1 :
        T0 = ZERO;
        n_frequency = 1;
        amplitude = HUNDRED;
        rod_kappa = Kappa;
        rod_length = Soln[Number_of_Nodes-1].X.x - Soln[0].X.x;
        for ( i = 0 ; i <= Number_of_Nodes-1 ; ++i ) {
            Soln[i].T = T0 + 
               amplitude*sqr(sin(double(n_frequency)*PI*Soln[i].X.x/rod_length));
            Soln[i].kappa = rod_kappa;
        } /* endfor */
        break;
      case IC_SINUSOIDAL_ROD4 :
        T0 = ZERO;
        n_frequency = 4;
        amplitude = HUNDRED;
        rod_kappa = Kappa;
        rod_length = Soln[Number_of_Nodes-1].X.x - Soln[0].X.x;
        for ( i = 0 ; i <= Number_of_Nodes-1 ; ++i ) {
            Soln[i].T = T0 + 
               amplitude*sqr(sin(double(n_frequency)*PI*Soln[i].X.x/rod_length));
            Soln[i].kappa = rod_kappa;
        } /* endfor */
        break;
      default:
        T0 = TEN;
        rod_kappa = Kappa;
        rod_length = Soln[Number_of_Nodes-1].X.x - Soln[0].X.x;
        for ( i = 0 ; i <= Number_of_Nodes-1 ; ++i ) {
            Soln[i].T = T0;
            Soln[i].kappa = rod_kappa;
        } /* endfor */
        break;
    } /* endswitch */

}

/********************************************************
 * Routine: AnalyticSoln                                *
 *                                                      *
 * Computes the exact analytical solution.              *
 *                                                      *
 ********************************************************/
void AnalyticSoln(Heat1D_UniformMesh *Soln,
                  const int i_ICtype,
                  const double &Kappa,
                  const double &Time,
                  const int Number_of_Nodes) {

    int i, k, n_frequency;
    double T0, rod_length, rod_kappa, amplitude, 
           temp, dtdx, gamma, delta, beta, gamma_k;

    /* Assign the initial data for the IVP of interest. */
    
    switch(i_ICtype) {
      case IC_IMPULSIVE_ROD :
        rod_kappa = Kappa;
        rod_length = Soln[Number_of_Nodes-1].X.x - Soln[0].X.x;
        T0 = HUNDRED;
        if (Time <= ZERO) {
           for ( i = 0 ; i <= Number_of_Nodes-1 ; ++i ) {
               Soln[i].T = ZERO;
               Soln[i].qx = ZERO;
               Soln[i].kappa = rod_kappa;
           } /* endfor */
           Soln[Number_of_Nodes-1].T = T0;
           Soln[Number_of_Nodes-1].qx = MILLION;
        } else {
           for ( i = 0 ; i <= Number_of_Nodes-1 ; ++i ) {
               gamma = TWO*T0/PI;
               delta = -rod_kappa*Time*sqr(PI/rod_length);
               beta = PI*Soln[i].X.x/rod_length;
               k = 1;
               temp = T0*Soln[i].X.x/rod_length;
               dtdx = T0/rod_length;
               while (1) {
                   gamma_k = gamma*exp(delta*double(sqr(k)))*
                             pow(-ONE, double(k))/double(k);
                   if (fabs(gamma_k/T0) < TOLER) break;
                   temp = temp + gamma_k*sin(double(k)*beta);
                   dtdx = dtdx + gamma_k*double(k)*PI*
                                 cos(double(k)*beta)/rod_length;
                   k = k + 1;
               } /* endwhile */
               Soln[i].T = temp;
               Soln[i].kappa = rod_kappa;
               Soln[i].qx = -rod_kappa*dtdx;
           } /* endfor */
        } /* endif */
        break;
      case IC_SINUSOIDAL_ROD1 :
        T0 = ZERO;
        n_frequency = 1;
        amplitude = HUNDRED;
        rod_kappa = Kappa;
        rod_length = Soln[Number_of_Nodes-1].X.x - Soln[0].X.x;
        if (Time <= ZERO) {
           for ( i = 0 ; i <= Number_of_Nodes-1 ; ++i ) {
               Soln[i].T = T0 + 
                  amplitude*sqr(sin(double(n_frequency)*PI*Soln[i].X.x/rod_length));
               Soln[i].kappa = rod_kappa;
               Soln[i].qx = -rod_kappa*TWO*amplitude*
                            (double(n_frequency)*PI/rod_length)*
                            sin(double(n_frequency)*PI*Soln[i].X.x/rod_length)*
                            cos(double(n_frequency)*PI*Soln[i].X.x/rod_length);
           } /* endfor */
        } else {
           for ( i = 0 ; i <= Number_of_Nodes-1 ; ++i ) {
               gamma = EIGHT*amplitude*double(sqr(n_frequency))/PI;
               delta = -rod_kappa*Time*sqr(PI/rod_length);
               beta = PI*Soln[i].X.x/rod_length;
               k = 1;
               temp = ZERO;
               dtdx = ZERO;
               while (1) {
                   gamma_k = gamma*exp(delta*double(sqr(k)))/
                             double(k*(FOUR*sqr(n_frequency)-sqr(k)));
                   if (fabs(gamma_k/amplitude) < TOLER) break;
                   temp = temp + gamma_k*sin(double(k)*beta);
                   dtdx = dtdx + gamma_k*double(k)*PI*
                                 cos(double(k)*beta)/rod_length;
                   k = k + 2;
               } /* endwhile */
               Soln[i].T = temp;
               Soln[i].kappa = rod_kappa;
               Soln[i].qx = -rod_kappa*dtdx;
           } /* endfor */
        } /* endif */
        break;
      case IC_SINUSOIDAL_ROD4 :
        T0 = ZERO;
        n_frequency = 4;
        amplitude = HUNDRED;
        rod_kappa = Kappa;
        rod_length = Soln[Number_of_Nodes-1].X.x - Soln[0].X.x;
        if (Time <= ZERO) {
           for ( i = 0 ; i <= Number_of_Nodes-1 ; ++i ) {
               Soln[i].T = T0 + 
                  amplitude*sqr(sin(double(n_frequency)*PI*Soln[i].X.x/rod_length));
               Soln[i].kappa = rod_kappa;
               Soln[i].qx = -rod_kappa*TWO*amplitude*
                            (double(n_frequency)*PI/rod_length)*
                            sin(double(n_frequency)*PI*Soln[i].X.x/rod_length)*
                            cos(double(n_frequency)*PI*Soln[i].X.x/rod_length);
           } /* endfor */
        } else {
           for ( i = 0 ; i <= Number_of_Nodes-1 ; ++i ) {
               gamma = EIGHT*amplitude*double(sqr(n_frequency))/PI;
               delta = -rod_kappa*Time*sqr(PI/rod_length);
               beta = PI*Soln[i].X.x/rod_length;
               k = 1;
               temp = ZERO;
               dtdx = ZERO;
               while (1) {
                   gamma_k = gamma*exp(delta*double(sqr(k)))/
                             double(k*(FOUR*sqr(n_frequency)-sqr(k)));
                   if (fabs(gamma_k/amplitude) < TOLER) break;
                   temp = temp + gamma_k*sin(double(k)*beta);
                   dtdx = dtdx + gamma_k*double(k)*PI*
                                 cos(double(k)*beta)/rod_length;
                   k = k + 2;
               } /* endwhile */
               Soln[i].T = temp;
               Soln[i].kappa = rod_kappa;
               Soln[i].qx = -rod_kappa*dtdx;
           } /* endfor */
        } /* endif */
        break;
      default:
        rod_kappa = Kappa;
        rod_length = Soln[Number_of_Nodes-1].X.x - Soln[0].X.x;
        T0 = HUNDRED;
        if (Time <= ZERO) {
           for ( i = 0 ; i <= Number_of_Nodes-1 ; ++i ) {
               Soln[i].T = ZERO;
               Soln[i].qx = ZERO;
               Soln[i].kappa = rod_kappa;
           } /* endfor */
           Soln[Number_of_Nodes-1].T = T0;
           Soln[Number_of_Nodes-1].qx = MILLION;
        } else {
           for ( i = 0 ; i <= Number_of_Nodes-1 ; ++i ) {
               gamma = TWO*T0/PI;
               delta = -rod_kappa*Time*sqr(PI/rod_length);
               beta = PI*Soln[i].X.x/rod_length;
               k = 1;
               temp = T0*Soln[i].X.x/rod_length;
               dtdx = T0/rod_length;
               while (1) {
                   gamma_k = gamma*exp(delta*double(sqr(k)))*
                             pow(-ONE, double(k))/double(k);
                   if (fabs(gamma_k/T0) < TOLER) break;
                   temp = temp + gamma_k*sin(double(k)*beta);
                   dtdx = dtdx + gamma_k*double(k)*PI*
                                 cos(double(k)*beta)/rod_length;
                   k = k + 1;
               } /* endwhile */
               Soln[i].T = temp;
               Soln[i].kappa = rod_kappa;
               Soln[i].qx = -rod_kappa*dtdx;
           } /* endfor */
        } /* endif */
        break;
    } /* endswitch */

}

/********************************************************
 * Routine: HeatFlux                                    *
 *                                                      *
 * Computes heat flux at each node in solution domain.  *
 *                                                      *
 ********************************************************/
void HeatFlux(Heat1D_UniformMesh *Soln,
              const int Number_of_Nodes) {

    int i;

    /* Determine heat flux. */

    for ( i = 1 ; i <= Number_of_Nodes-2 ; ++i ) {
        Soln[i].qx = -HALF*Soln[i].kappa*
                      (Soln[i+1].T-Soln[i-1].T)/
                      Soln[i].X.dx;
    } /* endfor */

    Soln[0].qx = -HALF*Soln[0].kappa*
                  (-THREE*Soln[0].T
                   +FOUR*Soln[1].T
                   -Soln[2].T)/
                  Soln[0].X.dx;

    Soln[Number_of_Nodes-1].qx = 
                 -HALF*Soln[Number_of_Nodes-1].kappa*
                  (THREE*Soln[Number_of_Nodes-1].T
                   -FOUR*Soln[Number_of_Nodes-2].T
                   +Soln[Number_of_Nodes-3].T)/
                  Soln[Number_of_Nodes-1].X.dx;

}

/********************************************************
 * Routine: Rstability                                  *
 *                                                      *
 * Determines the allowable global and local time steps *
 * (for simple explicit time stepping scheme) according *
 * to the linear stability condition.                   *
 *                                                      *
 ********************************************************/
double Rstability(Heat1D_UniformMesh *Soln,
                  const int Number_of_Nodes) {

    int i;
    double dtMin;

    /* Determine local and global time steps. */

    dtMin = MILLION;

    for ( i = 0 ; i <= Number_of_Nodes-1 ; ++i ) {
        Soln[i].dt = HALF*sqr(Soln[i].X.dx)/Soln[i].kappa;
        dtMin = min(dtMin, Soln[i].dt);
    } /* endfor */

    /* Return the global time step. */

    return (dtMin);

}

/********************************************************
 * Routine: dTdt_SimpleExplicit                         *
 *                                                      *
 * This routine updates the solution using a 1st-order  *
 * explicit Euler time integration and second-order     *
 * centered-difference discretization scheme.           *
 *                                                      *
 ********************************************************/
int dTdt_SimpleExplicit(Heat1D_UniformMesh *Soln,
	                const int Number_of_Nodes,
			double &dtMin,
			const double &CFL_Number,
			const int Local_Time_Stepping) {

    int i;

    /* Evaluate the time rate of change of the solution
       (i.e., the solution residuals) using a 1st-order
       explicit Euler time integration and second-order
       centered-difference discretization scheme. */

    Soln[0].dTdt = ZERO;
    for ( i = 1 ; i <= Number_of_Nodes-2 ; ++i ) {
        Soln[i].dTdt = ZERO;
        if ( !Local_Time_Stepping ) Soln[i].dt = dtMin;
        Soln[i].dTdt = (CFL_Number*Soln[i].dt)*Soln[i].kappa*
                       (Soln[i+1].T-TWO*Soln[i].T+Soln[i-1].T)/
                       sqr(Soln[i].X.dx);
    } /* endfor */
    Soln[Number_of_Nodes-1].dTdt = ZERO;

    /* Now update the solution. */

    for ( i = 1 ; i <= Number_of_Nodes-2 ; ++i ) {
       Soln[i].T += Soln[i].dTdt;
    } /* endfor */

    /* By default, boundary values at either end of the 
       mesh remain unchanged. */

    Soln[0].T = Soln[0].T;
    Soln[Number_of_Nodes-1].T = Soln[Number_of_Nodes-1].T;

    /* Solution successfully updated. */

    return (0);
    
}

/********************************************************
 * Routine: dTdt_SimpleImplicit                         *
 *                                                      *
 * This routine updates the solution using a 1st-order  *
 * implicit Euler time integration and second-order     *
 * centered-difference discretization scheme.           *
 *                                                      *
 ********************************************************/
int dTdt_SimpleImplicit(Heat1D_UniformMesh *Soln,
	                const int Number_of_Nodes,
		        double &dtMin,
		        const double &CFL_Number,
		        const int Local_Time_Stepping,
                        TriDiagonalSystemLinEqs &Ax_equal_b) {

    int i;
    double r;

    /* Evaluate the time rate of change of the solution
       (i.e., the solution residuals) using a 1st-order
       implicit Euler time integration and second-order
       centered-difference discretization scheme.
       Note that this unconditionally stable technique 
       requires the solution of a block tridiagonal 
       system of equations. */

    Soln[0].dTdt = ZERO;
    for ( i = 1 ; i <= Number_of_Nodes-2 ; ++i ) {
        Soln[i].dTdt = ZERO;
        if ( !Local_Time_Stepping ) Soln[i].dt = dtMin;
        r = CFL_Number*Soln[i].dt*Soln[i].kappa/sqr(Soln[i].X.dx);
        Ax_equal_b.A.A(i)=-r;
        Ax_equal_b.A.D(i)=ONE+TWO*r;
        Ax_equal_b.A.B(i)=-r;
        Ax_equal_b.b(i)=r*(Soln[i+1].T-TWO*Soln[i].T+Soln[i-1].T);
    } /* endfor */
    Soln[Number_of_Nodes-1].dTdt = ZERO;

    Ax_equal_b.A.A(0)=ZERO;
    Ax_equal_b.A.D(0)=ONE;
    Ax_equal_b.A.B(0)=ZERO;
    Ax_equal_b.b(0)=ZERO;

    Ax_equal_b.A.A(Number_of_Nodes-1)=ZERO;
    Ax_equal_b.A.D(Number_of_Nodes-1)=ONE;
    Ax_equal_b.A.B(Number_of_Nodes-1)=ZERO;
    Ax_equal_b.b(Number_of_Nodes-1)=ZERO;

    Ax_equal_b.solve();

    /* Now assign the residual and update the solution. */

    for ( i = 1 ; i <= Number_of_Nodes-2 ; ++i ) {
        Soln[i].dTdt = Ax_equal_b.x(i);
        Soln[i].T += Soln[i].dTdt;
    } /* endfor */

    /* By default, boundary values at either end of the 
       mesh remain unchanged. */

    Soln[0].T = Soln[0].T;
    Soln[Number_of_Nodes-1].T = Soln[Number_of_Nodes-1].T;

    /* Solution successfully updated. */

    return (0);
    
}

/********************************************************
 * Routine: dTdt_CrankNicolson                          *
 *                                                      *
 * This routine updates the solution using the          *
 * second-order accurate Crank-Nicolson scheme.         *
 * See Crank and Nicolson (1947).                       *
 *                                                      *
 ********************************************************/
int dTdt_CrankNicolson(Heat1D_UniformMesh *Soln,
	               const int Number_of_Nodes,
		       double &dtMin,
		       const double &CFL_Number,
		       const int Local_Time_Stepping,
                       TriDiagonalSystemLinEqs &Ax_equal_b) {

    int i;
    double r;

    /* Evaluate the time rate of change of the solution
       (i.e., the solution residuals) using a second-order
       trapezoidal implicit time integration and second-order
       centered-difference discretization scheme.  Note that
       this unconditionally stable technique requires the 
       solution of a block tridiagonal system of equations. */

    Soln[0].dTdt = ZERO;
    for ( i = 1 ; i <= Number_of_Nodes-2 ; ++i ) {
        Soln[i].dTdt = ZERO;
        if ( !Local_Time_Stepping ) Soln[i].dt = dtMin;
        r = CFL_Number*Soln[i].dt*Soln[i].kappa/sqr(Soln[i].X.dx);
        Ax_equal_b.A.A(i)=-HALF*r;
        Ax_equal_b.A.D(i)=ONE+r;
        Ax_equal_b.A.B(i)=-HALF*r;
        Ax_equal_b.b(i)=r*(Soln[i+1].T-TWO*Soln[i].T+Soln[i-1].T);
    } /* endfor */
    Soln[Number_of_Nodes-1].dTdt = ZERO;

    Ax_equal_b.A.A(0)=ZERO;
    Ax_equal_b.A.D(0)=ONE;
    Ax_equal_b.A.B(0)=ZERO;
    Ax_equal_b.b(0)=ZERO;

    Ax_equal_b.A.A(Number_of_Nodes-1)=ZERO;
    Ax_equal_b.A.D(Number_of_Nodes-1)=ONE;
    Ax_equal_b.A.B(Number_of_Nodes-1)=ZERO;
    Ax_equal_b.b(Number_of_Nodes-1)=ZERO;

    Ax_equal_b.solve();

    /* Now assign the residual and update the solution. */

    for ( i = 1 ; i <= Number_of_Nodes-2 ; ++i ) {
        Soln[i].dTdt = Ax_equal_b.x(i);
        Soln[i].T += Soln[i].dTdt;
    } /* endfor */

    /* By default, boundary values at either end of the 
       mesh remain unchanged. */

    Soln[0].T = Soln[0].T;
    Soln[Number_of_Nodes-1].T = Soln[Number_of_Nodes-1].T;

    /* Solution successfully updated. */

    return (0);
    
}

/********************************************************
 * Routine: dTdt_ADE                                    *
 *                                                      *
 * This routine updates the solution using the          *
 * second-order accurate alternating direction explicit *
 * scheme of Barakat and Clark.  See Barakat and Clark  *
 * (1966).                                              *
 *                                                      *
 ********************************************************/
int dTdt_ADE(Heat1D_UniformMesh *Soln,
	     const int Number_of_Nodes,
	     double &dtMin,
	     const double &CFL_Number,
	     const int Local_Time_Stepping) {

    int i;
    double qnew, pnew, dqdt, r; 

    /* Evaluate the time rate of change of the solution
       (i.e., the solution residuals) using the alternating
       direction explicit method of Barakat and Clark. */

    pnew = Soln[0].T;
    Soln[0].dTdt = ZERO;
    for ( i = 1 ; i <= Number_of_Nodes-2 ; ++i ) {
        Soln[i].dTdt = ZERO;
        if ( !Local_Time_Stepping ) Soln[i].dt = dtMin;
        r = CFL_Number*Soln[i].dt*Soln[i].kappa/sqr(Soln[i].X.dx);
        Soln[i].dTdt += (CFL_Number*Soln[i].dt)*Soln[i].kappa*
                        (Soln[i+1].T-TWO*Soln[i].T+pnew)/
                        (sqr(Soln[i].X.dx)*(ONE+r));
        pnew = Soln[i].T+Soln[i].dTdt;
    } /* endfor */
    Soln[Number_of_Nodes-1].dTdt = ZERO;

    qnew = Soln[Number_of_Nodes-1].T;
    for ( i = Number_of_Nodes-2 ; i >= 1 ; --i ) {
        r = CFL_Number*Soln[i].dt*Soln[i].kappa/sqr(Soln[i].X.dx);
        dqdt = (CFL_Number*Soln[i].dt)*Soln[i].kappa*
               (qnew-TWO*Soln[i].T+Soln[i-1].T)/
               (sqr(Soln[i].X.dx)*(ONE+r));
        Soln[i].dTdt += dqdt;
        qnew = Soln[i].T+dqdt;
    } /* endfor */

    /* Now update the solution. */

    for ( i = 1 ; i <= Number_of_Nodes-2 ; ++i ) {
       Soln[i].T += HALF*Soln[i].dTdt;
    } /* endfor */

    /* By default, boundary values at either end of the 
       mesh remain unchanged. */

    Soln[0].T = Soln[0].T;
    Soln[Number_of_Nodes-1].T = Soln[Number_of_Nodes-1].T;

    /* Solution successfully updated. */

    return (0);
    
}

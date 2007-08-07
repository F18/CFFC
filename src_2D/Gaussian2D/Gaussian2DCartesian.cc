/* Gaussian2DCartesian.cc:  Subroutines for 2D Gaussian 
                            Cartesian Mesh Solution Classes. */

/* Include 2D aussian Cartesian mesh solution header file. */

#ifndef _GAUSSIAN2D_CARTESIAN_INCLUDED
#include "Gaussian2DCartesian.h"
#endif // _GAUSSIAN2D_CARTESIAN_INCLUDED

/*************************************************************
 * Gaussian2D_Cartesian_Cell -- External subroutines.        *
 *************************************************************/

/********************************************************
 * Routine: Allocate                                    *
 *                                                      *
 * Allocate memory for 2D Gaussian equation solution.   *
 *                                                      *
 ********************************************************/
Gaussian2D_Cartesian_Cell** Allocate(Gaussian2D_Cartesian_Cell **Soln_ptr,
                                     const int Number_of_Cells_Xdir,
  	                             const int Number_of_Cells_Ydir) {

    int i;
 
    /* Allocate memory. */

    Soln_ptr = new Gaussian2D_Cartesian_Cell*[Number_of_Cells_Xdir+2];
    for ( i = 0 ; i <= Number_of_Cells_Xdir+1 ; ++i ) {
      Soln_ptr[i] = new Gaussian2D_Cartesian_Cell[Number_of_Cells_Ydir+2];
    }  /* endfor */

    /* Return memory location. */

    return(Soln_ptr);

}

/********************************************************
 * Routine: Deallocate                                  *
 *                                                      *
 * Deallocate memory for 2D Gaussian equation solution. *
 *                                                      *
 ********************************************************/
Gaussian2D_Cartesian_Cell** Deallocate(Gaussian2D_Cartesian_Cell **Soln_ptr,
                                       const int Number_of_Cells_Xdir,
  	                               const int Number_of_Cells_Ydir) {

    int i;
 
    /* Deallocate memory. */

    for ( i = 0 ; i <= Number_of_Cells_Xdir+1 ; ++i ) {
       delete []Soln_ptr[i];
       Soln_ptr[i] = NULL;
    }  /* endfor */
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
void Output_Gnuplot(Gaussian2D_Cartesian_Cell **Soln,
  		    const int Number_of_Cells_Xdir,
  		    const int Number_of_Cells_Ydir,
                    const int Number_of_Time_Steps,
                    const double &Time,
	            ostream &out_file) {

    int i, j;

    out_file << "# " << CFDkit_Name() << ": 2D Gaussian Solution, "
             << "Time Step/Iteration Level = "
             << Number_of_Time_Steps
             << ", Time = " << Time*THOUSAND << " (ms)\n"
	     << "# i, j, x(m), y(m), rho (kg/m^3), "
	     << "u (m/s), v (m/s), p (Pa), T (K)\n";

    for (j  = 1 ; j <= Number_of_Cells_Ydir ; ++j ) {
       for ( i = 1 ; i <= Number_of_Cells_Xdir ; ++i ) {
	   out_file << " " << i << " " << j << " " << Soln[i][j] << "\n";
       } /* endfor */
       out_file << "\n";
    } /* endfor */

}

/********************************************************
 * Routine: Output_Tecplot                              *
 *                                                      *
 * Writes the solution to specified output stream       *
 * suitable for plotting with TECPLOT.                  *
 *                                                      *
 ********************************************************/
void Output_Tecplot(Gaussian2D_Cartesian_Cell **Soln,
  		    const int Number_of_Cells_Xdir,
  		    const int Number_of_Cells_Ydir,
                    const int Number_of_Time_Steps,
                    const double &Time,
	            ostream &out_file) {

    int i, j;

    out_file << "TITLE = \"" << CFDkit_Name() << ": 2D Gaussian Solution, "
             << "Time Step/Iteration Level = "
             << Number_of_Time_Steps
             << ", Time = " << Time << "\"" << "\n"
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
             << "\"T\" \n"
             << "ZONE \\ \n"
             << "I = " << Number_of_Cells_Xdir << " \\ \n"
             << "J = " << Number_of_Cells_Ydir << " \\ \n"
             << "F = POINT \n";
    
    for (j  = 1 ; j <= Number_of_Cells_Ydir ; ++j ) {
       for ( i = 1 ; i <= Number_of_Cells_Xdir ; ++i ) {
	   out_file << " " << Soln[i][j] << "\n";
       } /* endfor */
    } /* endfor */

}

/********************************************************
 * Routine: Grid                                        *
 *                                                      *
 * Generates a uniform 2D Cartesian mesh and assign the *
 * locations of the cell centers to appropriate         *
 * solution variables.                                  *
 *                                                      *
 ********************************************************/
void Grid(Gaussian2D_Cartesian_Cell **Soln,
          const double &xMin,
	  const double &xMax,
          const double &yMin,
	  const double &yMax,
  	  const int Number_of_Cells_Xdir,
  	  const int Number_of_Cells_Ydir) {

    int i, j;
    double delta_x, delta_y;

    /* Determine the mesh spacing. */

    delta_x = (xMax - xMin)/double(Number_of_Cells_Xdir);
    delta_y = (yMax - yMin)/double(Number_of_Cells_Ydir);
    Soln[0][0].X.setsize(delta_x, delta_y);

    /* Create the cells. */

    for ( j = 0 ; j <= Number_of_Cells_Ydir+1 ; ++j ) {
       for ( i = 0 ; i <= Number_of_Cells_Xdir+1 ; ++i ) {
           Soln[i][j].X.setloc(xMin+double(i-1)*delta_x, 
                               yMin+double(j-1)*delta_y);
       } /* endfor */
    } /* endfor */

}

/********************************************************
 * Routine: ICs                                         *
 *                                                      *
 * Assigns initial conditions and data to the           *
 * solution variables.                                  *
 *                                                      *
 ********************************************************/
void ICs(Gaussian2D_Cartesian_Cell **Soln,
         char *gas_ptr,
	 const int i_ICtype,
 	 const int Number_of_Cells_Xdir,
  	 const int Number_of_Cells_Ydir) {

    int i, j;
    Gaussian2D_pState Wl, Wr;

    /* Assign the gas constants for the gas of interest. */

    Soln[0][0].W.setgas(gas_ptr);
    Soln[0][0].U.setgas(gas_ptr);

    /* Assign the initial data for the IVP of interest. */

    switch(i_ICtype) {
      case IC_CONSTANT :
        Wl = Gaussian2D_W_STDATM;
        for ( j = 0 ; j <= Number_of_Cells_Ydir+1 ; ++j ) {
           for ( i = 0 ; i <= Number_of_Cells_Xdir+1 ; ++i ) {
              Soln[i][j].W = Wl;
              Soln[i][j].U = U(Soln[i][j].W);
           } /* endfor */
        } /* endfor */
        break;
      case IC_UNIFORM :
        Wl = Gaussian2D_W_STDATM;
        for ( j = 0 ; j <= Number_of_Cells_Ydir+1 ; ++j ) {
           for ( i = 0 ; i <= Number_of_Cells_Xdir+1 ; ++i ) {
              Soln[i][j].W = Wl;
              Soln[i][j].U = U(Soln[i][j].W);
           } /* endfor */
        } /* endfor */
        break;
      case IC_SOD_XDIR :
        Wl = Gaussian2D_W_STDATM;
        Wr = Gaussian2D_pState(DENSITY_STDATM/EIGHT,
          		    ZERO, ZERO,
         		    PRESSURE_STDATM/TEN); 
        for ( j = 0 ; j <= Number_of_Cells_Ydir+1 ; ++j ) {
           for ( i = 0 ; i <= Number_of_Cells_Xdir+1 ; ++i ) {
              if (Soln[i][j].X.xc.x <= ZERO) {
                 Soln[i][j].W = Wl;
              } else {
                 Soln[i][j].W = Wr;	     
              } /* end if */
              Soln[i][j].U = U(Soln[i][j].W);
           } /* endfor */
        } /* endfor */
        break;
      case IC_SOD_YDIR :
        Wl = Gaussian2D_W_STDATM;
        Wr = Gaussian2D_pState(DENSITY_STDATM/EIGHT,
          		    ZERO, ZERO,
         		    PRESSURE_STDATM/TEN); 
        for ( j = 0 ; j <= Number_of_Cells_Ydir+1 ; ++j ) {
           for ( i = 0 ; i <= Number_of_Cells_Xdir+1 ; ++i ) {
              if (Soln[i][j].X.xc.y <= ZERO) {
                 Soln[i][j].W = Wl;
              } else {
                 Soln[i][j].W = Wr;	     
              } /* end if */
              Soln[i][j].U = U(Soln[i][j].W);
           } /* endfor */
        } /* endfor */
        break;
      case IC_GROTH_XDIR :
        Wl = Gaussian2D_pState(4.696, ZERO, ZERO, 404.4e03);
        Wr = Gaussian2D_pState(1.408, ZERO, ZERO, 101.1e03);
        for ( j = 0 ; j <= Number_of_Cells_Ydir+1 ; ++j ) {
           for ( i = 0 ; i <= Number_of_Cells_Xdir+1 ; ++i ) {
              if (Soln[i][j].X.xc.x <= ZERO) {
                 Soln[i][j].W = Wl;
              } else {
                 Soln[i][j].W = Wr;	     
              } /* end if */
              Soln[i][j].U = U(Soln[i][j].W);
           } /* endfor */
        } /* endfor */
        break;
      case IC_GROTH_YDIR :
        Wl = Gaussian2D_pState(4.696, ZERO, ZERO, 404.4e03);
        Wr = Gaussian2D_pState(1.408, ZERO, ZERO, 101.1e03);
        for ( j = 0 ; j <= Number_of_Cells_Ydir+1 ; ++j ) {
           for ( i = 0 ; i <= Number_of_Cells_Xdir+1 ; ++i ) {
              if (Soln[i][j].X.xc.y <= ZERO) {
                 Soln[i][j].W = Wl;
              } else {
                 Soln[i][j].W = Wr;	     
              } /* end if */
              Soln[i][j].U = U(Soln[i][j].W);
           } /* endfor */
        } /* endfor */
        break;
      case IC_EINFELDT_XDIR :
        Wl = Gaussian2D_pState(ONE, -TWO, ZERO, FOUR/TEN);
        Wr = Gaussian2D_pState(ONE, TWO, ZERO, FOUR/TEN);
        for ( j = 0 ; j <= Number_of_Cells_Ydir+1 ; ++j ) {
           for ( i = 0 ; i <= Number_of_Cells_Xdir+1 ; ++i ) {
              if (Soln[i][j].X.xc.x <= ZERO) {
                 Soln[i][j].W = Wl;
              } else {
                 Soln[i][j].W = Wr;	     
              } /* end if */
              Soln[i][j].U = U(Soln[i][j].W);
           } /* endfor */
        } /* endfor */
        break;
      case IC_EINFELDT_YDIR :
        Wl = Gaussian2D_pState(ONE, ZERO, -TWO, FOUR/TEN);
        Wr = Gaussian2D_pState(ONE, ZERO, TWO, FOUR/TEN);
        for ( j = 0 ; j <= Number_of_Cells_Ydir+1 ; ++j ) {
           for ( i = 0 ; i <= Number_of_Cells_Xdir+1 ; ++i ) {
              if (Soln[i][j].X.xc.y <= ZERO) {
                 Soln[i][j].W = Wl;
              } else {
                 Soln[i][j].W = Wr;	     
              } /* end if */
              Soln[i][j].U = U(Soln[i][j].W);
           } /* endfor */
        } /* endfor */
        break;
      case IC_SHOCK_WAVE_XDIR :
        Wl = Gaussian2D_pState(2.281, 164.83, ZERO, 201.17e03);
        Wr = Gaussian2D_pState(1.408, ZERO, ZERO, 101.1e03);
        for ( j = 0 ; j <= Number_of_Cells_Ydir+1 ; ++j ) {
           for ( i = 0 ; i <= Number_of_Cells_Xdir+1 ; ++i ) {
              if (Soln[i][j].X.xc.x <= ZERO) {
                 Soln[i][j].W = Wl;
              } else {
                 Soln[i][j].W = Wr;	     
              } /* end if */
              Soln[i][j].U = U(Soln[i][j].W);
           } /* endfor */
        } /* endfor */
        break;
      case IC_SHOCK_WAVE_YDIR :
        Wl = Gaussian2D_pState(2.281, ZERO, 164.83, 201.17e03);
        Wr = Gaussian2D_pState(1.408, ZERO, ZERO, 101.1e03);
        for ( j = 0 ; j <= Number_of_Cells_Ydir+1 ; ++j ) {
           for ( i = 0 ; i <= Number_of_Cells_Xdir+1 ; ++i ) {
              if (Soln[i][j].X.xc.y <= ZERO) {
                 Soln[i][j].W = Wl;
              } else {
                 Soln[i][j].W = Wr;	     
              } /* end if */
              Soln[i][j].U = U(Soln[i][j].W);
           } /* endfor */
        } /* endfor */
        break;
      case IC_SHOCK_STRUCTURE_M1_1:
	Wl = Gaussian2D_pState(1.661, 350.7444241833, 0.0, 101325.0);
	Wr = Gaussian2D_pState(1.90955819477, 305.08951614, 0.0, 127922.8125);
        for ( j = 0 ; j <= Number_of_Cells_Ydir+1 ; ++j ) {
           for ( i = 0 ; i <= Number_of_Cells_Xdir+1 ; ++i ) {
              if (Soln[i][j].X.xc.x <= ZERO) {
                 Soln[i][j].W = Wl;
              } else {
                 Soln[i][j].W = Wr;	     
              } /* end if */
              Soln[i][j].U = U(Soln[i][j].W);
           } /* endfor */
        } /* endfor */
        break;
      case IC_SHOCK_STRUCTURE_M1_3:
	Wl = Gaussian2D_pState(1.661, 414.51592216, 0.0, 101325.0);
	Wr = Gaussian2D_pState(2.39410660981, 287.585750733, 0.0, 188717.8125);
        for ( j = 0 ; j <= Number_of_Cells_Ydir+1 ; ++j ) {
           for ( i = 0 ; i <= Number_of_Cells_Xdir+1 ; ++i ) {
              if (Soln[i][j].X.xc.x <= ZERO) {
                 Soln[i][j].W = Wl;
              } else {
                 Soln[i][j].W = Wr;	     
              } /* end if */
              Soln[i][j].U = U(Soln[i][j].W);
           } /* endfor */
        } /* endfor */
        break;
      case IC_SHOCK_STRUCTURE_M1_5:
	Wl = Gaussian2D_pState(1.661, 478.287602499, 0.0, 101325.0);
	Wr = Gaussian2D_pState(2.84742857143, 279.001101458, 0.0, 259645.3125);
        for ( j = 0 ; j <= Number_of_Cells_Ydir+1 ; ++j ) {
           for ( i = 0 ; i <= Number_of_Cells_Xdir+1 ; ++i ) {
              if (Soln[i][j].X.xc.x <= ZERO) {
                 Soln[i][j].W = Wl;
              } else {
                 Soln[i][j].W = Wr;	     
              } /* end if */
              Soln[i][j].U = U(Soln[i][j].W);
           } /* endfor */
        } /* endfor */
        break;
      case IC_SHOCK_STRUCTURE_M2_0:
	Wl = Gaussian2D_pState(1.661, 637.716803332, 0.0, 101325.0);
	Wr = Gaussian2D_pState(3.79657142857, 279.001101458, 0.0, 481293.75);
        for ( j = 0 ; j <= Number_of_Cells_Ydir+1 ; ++j ) {
           for ( i = 0 ; i <= Number_of_Cells_Xdir+1 ; ++i ) {
              if (Soln[i][j].X.xc.x <= ZERO) {
                 Soln[i][j].W = Wl;
              } else {
                 Soln[i][j].W = Wr;	     
              } /* end if */
              Soln[i][j].U = U(Soln[i][j].W);
           } /* endfor */
        } /* endfor */
        break;
      case IC_SHOCK_STRUCTURE_M10_0:
	Wl = Gaussian2D_pState(1.661, 3188.58401666, 0.0, 101325.0);
	Wr = Gaussian2D_pState(6.45048543689, 821.06038429, 0.0, 12640293.75);
        for ( j = 0 ; j <= Number_of_Cells_Ydir+1 ; ++j ) {
           for ( i = 0 ; i <= Number_of_Cells_Xdir+1 ; ++i ) {
              if (Soln[i][j].X.xc.x <= ZERO) {
                 Soln[i][j].W = Wl;
              } else {
                 Soln[i][j].W = Wr;	     
              } /* end if */
              Soln[i][j].U = U(Soln[i][j].W);
           } /* endfor */
        } /* endfor */
        break;
      case IC_CONTACT_SURFACE_XDIR :
        Wl = Gaussian2D_pState(1.045, 200.00, ZERO, 300.00e03);
        Wr = Gaussian2D_pState(3.483, 200.00, ZERO, 300.00e03);
        for ( j = 0 ; j <= Number_of_Cells_Ydir+1 ; ++j ) {
           for ( i = 0 ; i <= Number_of_Cells_Xdir+1 ; ++i ) {
              if (Soln[i][j].X.xc.x <= ZERO) {
                 Soln[i][j].W = Wl;
              } else {
                 Soln[i][j].W = Wr;	     
              } /* end if */
              Soln[i][j].U = U(Soln[i][j].W);
           } /* endfor */
        } /* endfor */
        break;
      case IC_CONTACT_SURFACE_YDIR :
        Wl = Gaussian2D_pState(1.045, ZERO, 200.00, 300.00e03);
        Wr = Gaussian2D_pState(3.483, ZERO, 200.00, 300.00e03);
        for ( j = 0 ; j <= Number_of_Cells_Ydir+1 ; ++j ) {
           for ( i = 0 ; i <= Number_of_Cells_Xdir+1 ; ++i ) {
              if (Soln[i][j].X.xc.y <= ZERO) {
                 Soln[i][j].W = Wl;
              } else {
                 Soln[i][j].W = Wr;	     
              } /* end if */
              Soln[i][j].U = U(Soln[i][j].W);
           } /* endfor */
        } /* endfor */
        break;
      case IC_RAREFACTION_WAVE_XDIR :
        Wl = Gaussian2D_pState(1.598, -383.64, ZERO, 91.88e03);
        Wr = Gaussian2D_pState(2.787, -216.97, ZERO, 200.0e03);
        for ( j = 0 ; j <= Number_of_Cells_Ydir+1 ; ++j ) {
           for ( i = 0 ; i <= Number_of_Cells_Xdir+1 ; ++i ) {
              if (Soln[i][j].X.xc.x <= ZERO) {
                 Soln[i][j].W = Wl;
              } else {
                 Soln[i][j].W = Wr;	     
              } /* end if */
              Soln[i][j].U = U(Soln[i][j].W);
           } /* endfor */
        } /* endfor */
        break;
      case IC_RAREFACTION_WAVE_YDIR :
        Wl = Gaussian2D_pState(1.598, ZERO, -383.64, 91.88e03);
        Wr = Gaussian2D_pState(2.787, ZERO, -216.97, 200.0e03);
        for ( j = 0 ; j <= Number_of_Cells_Ydir+1 ; ++j ) {
           for ( i = 0 ; i <= Number_of_Cells_Xdir+1 ; ++i ) {
              if (Soln[i][j].X.xc.y <= ZERO) {
                 Soln[i][j].W = Wl;
              } else {
                 Soln[i][j].W = Wr;	     
              } /* end if */
              Soln[i][j].U = U(Soln[i][j].W);
           } /* endfor */
        } /* endfor */
        break;
      case IC_SHOCK_BOX :
        Wl = Gaussian2D_W_STDATM;
        Wr = Gaussian2D_pState(DENSITY_STDATM*FOUR,
         		    ZERO, ZERO,
         		    PRESSURE_STDATM*FOUR);
        for ( j = 0 ; j <= Number_of_Cells_Ydir+1 ; ++j ) {
           for ( i = 0 ; i <= Number_of_Cells_Xdir+1 ; ++i ) {
              if (Soln[i][j].X.xc.x <= HALF &&
                  Soln[i][j].X.xc.y <= HALF) {
                 Soln[i][j].W = Wl;
              } else {
                 Soln[i][j].W = Wr;	     
              } /* end if */
              Soln[i][j].U = U(Soln[i][j].W);
           } /* endfor */
        } /* endfor */
        break;
      default:
        Wl = Gaussian2D_W_STDATM;
        for ( j = 0 ; j <= Number_of_Cells_Ydir+1 ; ++j ) {
           for ( i = 0 ; i <= Number_of_Cells_Xdir+1 ; ++i ) {
              Soln[i][j].W = Wl;
              Soln[i][j].U = U(Soln[i][j].W);
           } /* endfor */
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
double CFL(Gaussian2D_Cartesian_Cell **Soln,
   	   const int Number_of_Cells_Xdir,
           const int Number_of_Cells_Ydir) {

    int i, j;
    double dtMin;

    /* Determine local and global time steps. */

    dtMin = MILLION;

    for ( j = 0 ; j <= Number_of_Cells_Ydir+1 ; ++j ) {
       for ( i = 0 ; i <= Number_of_Cells_Xdir+1 ; ++i ) {
	  Soln[i][j].dt = min(Soln[i][j].X.dx.x/
	                      (sqrt(3.0)*Soln[i][j].W.axx()+fabs(Soln[i][j].W.v.x)),
                              Soln[i][j].X.dx.y/
	                      (sqrt(3.0)*Soln[i][j].W.ayy()+fabs(Soln[i][j].W.v.y)));
          dtMin = min(dtMin, Soln[i][j].dt);
       } /* endfor */
    } /* endfor */

    /* Return the global time step. */

    return (dtMin);

}

/********************************************************
 * Routine: Linear_Reconstruction_GreenGauss            *
 *                                                      *
 * Peforms the reconstruction of a limited piecewise    *
 * linear solution state within each cell of the        *
 * computational mesh.  A Green-Gauss approach is used  *
 * in the evaluation of the unlimited solution          *
 * gradients.  Several slope limiters may be used.      *
 *                                                      *
 ********************************************************/
void Linear_Reconstruction_GreenGauss(Gaussian2D_Cartesian_Cell **Soln,
   		                      const int Number_of_Cells_Xdir,
                                      const int Number_of_Cells_Ydir,
				      const int Limiter) {

    int i, j, n;
    double u0Min, u0Max, uQuad[4], phi;

    /* Carry out the limited solution reconstruction in
       each cell of the computational mesh. */

    for ( j = 1 ; j <= Number_of_Cells_Ydir ; ++j ) {
       for ( i = 1 ; i <= Number_of_Cells_Xdir ; ++i ) {
          Soln[i][j].dWdx = HALF*(Soln[i+1][j].W-Soln[i-1][j].W)/
                            Soln[i][j].X.dx.x;
          Soln[i][j].dWdy = HALF*(Soln[i][j+1].W-Soln[i][j-1].W)/
                            Soln[i][j].X.dx.y;

          for ( n = 1 ; n <= NUM_VAR_GAUSSIAN2D ; ++n ) {
	     u0Min = min(Soln[i][j].W[n], Soln[i-1][j].W[n]);
             u0Min = min(u0Min, Soln[i+1][j].W[n]);
	     u0Min = min(u0Min, Soln[i][j+1].W[n]);
	     u0Min = min(u0Min, Soln[i][j-1].W[n]);
	     u0Max = max(Soln[i][j].W[n], Soln[i-1][j].W[n]);
	     u0Max = max(u0Max, Soln[i+1][j].W[n]);
	     u0Max = max(u0Max, Soln[i][j+1].W[n]);
	     u0Max = max(u0Max, Soln[i][j-1].W[n]);
             uQuad[0] = Soln[i][j].W[n] - HALF*Soln[i][j].dWdx[n]*Soln[i][j].X.dx.x; 
             uQuad[1] = Soln[i][j].W[n] + HALF*Soln[i][j].dWdx[n]*Soln[i][j].X.dx.x;
	     uQuad[2] = Soln[i][j].W[n] - HALF*Soln[i][j].dWdy[n]*Soln[i][j].X.dx.y;
	     uQuad[3] = Soln[i][j].W[n] + HALF*Soln[i][j].dWdy[n]*Soln[i][j].X.dx.y;

             switch(Limiter) {
               case LIMITER_ONE :
	         phi = ONE;
                 break;
               case LIMITER_ZERO :
	         phi = ZERO;
                 break;
               case LIMITER_BARTH_JESPERSEN :
                 phi = Limiter_BarthJespersen(uQuad, Soln[i][j].W[n], u0Min, u0Max, 4);
                 break;
               case LIMITER_VENKATAKRISHNAN :
                 phi = Limiter_Venkatakrishnan(uQuad, Soln[i][j].W[n], u0Min, u0Max, 4);
                 break;
               case LIMITER_VANLEER :
                 phi = Limiter_VanLeer(uQuad, Soln[i][j].W[n], u0Min, u0Max, 4);
                 break;
               case LIMITER_VANALBADA :
                 phi = Limiter_VanAlbada(uQuad, Soln[i][j].W[n], u0Min, u0Max, 4);
                 break;
               default:
                 phi = Limiter_BarthJespersen(uQuad, Soln[i][j].W[n], u0Min, u0Max, 4);
                 break;
             } /* endswitch */

             Soln[i][j].phi[n] = phi;
          } /* endfor */

       } /* endfor */

       Soln[0][j].dWdx = Soln[1][j].phi^Soln[1][j].dWdx;
       Soln[0][j].dWdy = Soln[1][j].phi^Soln[1][j].dWdy;
       Soln[0][j].phi = Gaussian2D_W_ONE;

       Soln[Number_of_Cells_Xdir+1][j].dWdx = 
         Soln[Number_of_Cells_Xdir][j].phi^Soln[Number_of_Cells_Xdir][j].dWdx;
       Soln[Number_of_Cells_Xdir+1][j].dWdy = 
         Soln[Number_of_Cells_Xdir][j].phi^Soln[Number_of_Cells_Xdir][j].dWdy;
       Soln[Number_of_Cells_Xdir+1][j].phi = Gaussian2D_W_ONE;

    } /* endfor */

    for ( i = 1 ; i <= Number_of_Cells_Xdir ; ++i ) {
       Soln[i][0].dWdx = Soln[i][1].phi^Soln[i][1].dWdx;
       Soln[i][0].dWdy = Soln[i][1].phi^Soln[i][1].dWdy;
       Soln[i][0].phi = Gaussian2D_W_ONE;

       Soln[i][Number_of_Cells_Ydir+1].dWdx = 
         Soln[i][Number_of_Cells_Ydir].phi^Soln[i][Number_of_Cells_Ydir].dWdx;
       Soln[i][Number_of_Cells_Ydir+1].dWdy = 
         Soln[i][Number_of_Cells_Ydir].phi^Soln[i][Number_of_Cells_Ydir].dWdy;
       Soln[i][Number_of_Cells_Ydir+1].phi = Gaussian2D_W_ONE;
    } /* endfor */
}

/********************************************************
 * Routine: Linear_Reconstruction_LeastSquares          *
 *                                                      *
 * Peforms the reconstruction of a limited piecewise    *
 * linear solution state within each cell of the        *
 * computational mesh.  A least squares approach is     *
 * used in the evaluation of the unlimited solution     *
 * gradients.  Several slope limiters may be used.      *
 *                                                      *
 ********************************************************/
void Linear_Reconstruction_LeastSquares(Gaussian2D_Cartesian_Cell **Soln,
   		                        const int Number_of_Cells_Xdir,
                                        const int Number_of_Cells_Ydir,
				        const int Limiter) {

    int i, j, n, n2, n_pts, i_index[8], j_index[8];
    double u0Min, u0Max, uQuad[4], phi;
    double DxDx_ave, DxDy_ave, DyDy_ave;
    Vector2D dX;
    Gaussian2D_pState DU, DUDx_ave, DUDy_ave;

    /* Carry out the limited solution reconstruction in
       each cell of the computational mesh. */

    for ( j = 1 ; j <= Number_of_Cells_Ydir ; ++j ) {
       for ( i = 1 ; i <= Number_of_Cells_Xdir ; ++i ) {
          n_pts = 8;
          i_index[0] = i-1; j_index[0] = j-1;
          i_index[1] = i  ; j_index[1] = j-1;
          i_index[2] = i+1; j_index[2] = j-1;
          i_index[3] = i-1; j_index[3] = j  ;
          i_index[4] = i+1; j_index[4] = j  ;
          i_index[5] = i-1; j_index[5] = j+1;
          i_index[6] = i  ; j_index[6] = j+1;
          i_index[7] = i+1; j_index[7] = j+1;

          DUDx_ave = Gaussian2D_W_VACUUM;
          DUDy_ave = Gaussian2D_W_VACUUM;
          DxDx_ave = ZERO;
          DxDy_ave = ZERO;
          DyDy_ave = ZERO;
    
          for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
              dX = Soln[ i_index[n2] ][ j_index[n2] ].X.xc - 
                   Soln[i][j].X.xc;
              DU = Soln[ i_index[n2] ][ j_index[n2] ].W - 
                   Soln[i][j].W;
              DUDx_ave += DU*dX.x;
              DUDy_ave += DU*dX.y;
              DxDx_ave += dX.x*dX.x;
              DxDy_ave += dX.x*dX.y;
              DyDy_ave += dX.y*dX.y;
          } /* endfor */
    					    
          DUDx_ave = DUDx_ave/double(n_pts);
          DUDy_ave = DUDy_ave/double(n_pts);
          DxDx_ave = DxDx_ave/double(n_pts);
          DxDy_ave = DxDy_ave/double(n_pts);
          DyDy_ave = DyDy_ave/double(n_pts);

          Soln[i][j].dWdx = (DUDx_ave*DyDy_ave-DUDy_ave*DxDy_ave)/
                            (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
          Soln[i][j].dWdy = (DUDy_ave*DxDx_ave-DUDx_ave*DxDy_ave)/
                            (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);

          for ( n = 1 ; n <= NUM_VAR_GAUSSIAN2D ; ++n ) {
             u0Min = Soln[i][j].W[n];
             u0Max = u0Min;
             for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
                u0Min = min(u0Min, Soln[ i_index[n2] ][ j_index[n2] ].W[n]);
                u0Max = max(u0Max, Soln[ i_index[n2] ][ j_index[n2] ].W[n]);
             } /* endfor */

             uQuad[0] = Soln[i][j].W[n] - HALF*Soln[i][j].dWdx[n]*Soln[i][j].X.dx.x; 
             uQuad[1] = Soln[i][j].W[n] + HALF*Soln[i][j].dWdx[n]*Soln[i][j].X.dx.x;
	     uQuad[2] = Soln[i][j].W[n] - HALF*Soln[i][j].dWdy[n]*Soln[i][j].X.dx.y;
	     uQuad[3] = Soln[i][j].W[n] + HALF*Soln[i][j].dWdy[n]*Soln[i][j].X.dx.y;

             switch(Limiter) {
               case LIMITER_ONE :
	         phi = ONE;
                 break;
               case LIMITER_ZERO :
	         phi = ZERO;
                 break;
               case LIMITER_BARTH_JESPERSEN :
                 phi = Limiter_BarthJespersen(uQuad, Soln[i][j].W[n], u0Min, u0Max, 4);
                 break;
               case LIMITER_VENKATAKRISHNAN :
                 phi = Limiter_Venkatakrishnan(uQuad, Soln[i][j].W[n], u0Min, u0Max, 4);
                 break;
               case LIMITER_VANLEER :
                 phi = Limiter_VanLeer(uQuad, Soln[i][j].W[n], u0Min, u0Max, 4);
                 break;
               case LIMITER_VANALBADA :
                 phi = Limiter_VanAlbada(uQuad, Soln[i][j].W[n], u0Min, u0Max, 4);
                 break;
               default:
                 phi = Limiter_BarthJespersen(uQuad, Soln[i][j].W[n], u0Min, u0Max, 4);
                 break;
             } /* endswitch */

             Soln[i][j].phi[n] = phi;
          } /* endfor */

       } /* endfor */

       Soln[0][j].dWdx = Soln[1][j].phi^Soln[1][j].dWdx;
       Soln[0][j].dWdy = Soln[1][j].phi^Soln[1][j].dWdy;
       Soln[0][j].phi = Gaussian2D_W_ONE;

       Soln[Number_of_Cells_Xdir+1][j].dWdx = 
         Soln[Number_of_Cells_Xdir][j].phi^Soln[Number_of_Cells_Xdir][j].dWdx;
       Soln[Number_of_Cells_Xdir+1][j].dWdy = 
         Soln[Number_of_Cells_Xdir][j].phi^Soln[Number_of_Cells_Xdir][j].dWdy;
       Soln[Number_of_Cells_Xdir+1][j].phi = Gaussian2D_W_ONE;

    } /* endfor */

    for ( i = 1 ; i <= Number_of_Cells_Xdir ; ++i ) {
       Soln[i][0].dWdx = Soln[i][1].phi^Soln[i][1].dWdx;
       Soln[i][0].dWdy = Soln[i][1].phi^Soln[i][1].dWdy;
       Soln[i][0].phi = Gaussian2D_W_ONE;

       Soln[i][Number_of_Cells_Ydir+1].dWdx = 
         Soln[i][Number_of_Cells_Ydir].phi^Soln[i][Number_of_Cells_Ydir].dWdx;
       Soln[i][Number_of_Cells_Ydir+1].dWdy = 
         Soln[i][Number_of_Cells_Ydir].phi^Soln[i][Number_of_Cells_Ydir].dWdy;
       Soln[i][Number_of_Cells_Ydir+1].phi = Gaussian2D_W_ONE;
    } /* endfor */
}

/********************************************************
 * Routine: dUdt_explicitEuler_upwind                   *
 *                                                      *
 * This routine updates the solution using a 1st-order  *
 * explicit Euler time integration and 1st-order upwind *
 * spatial discretization scheme in conjunction with    *
 * either the Godunov, Roe, Rusanov, HLLE, Linde, or    *
 * HLLC flux functions.                                 *
 *                                                      *
 ********************************************************/
int dUdt_explicitEuler_upwind(Gaussian2D_Cartesian_Cell **Soln,
   		              const int Number_of_Cells_Xdir,
                              const int Number_of_Cells_Ydir,
			      const int BC_type_left,
	      		      const int BC_type_right,
			      const int BC_type_bottom,
	      		      const int BC_type_top,
                              const Gaussian2D_pState W_left,
                              const Gaussian2D_pState W_right,
                              const Gaussian2D_pState W_bottom,
                              const Gaussian2D_pState W_top,
			      double &dtMin,
			      const double &CFL_Number,
			      const int Flux_Function_Type,
			      const int Local_Time_Stepping) {

    int i, j;
    Gaussian2D_cState Flux;

    /* Evaluate the time rate of change of the solution
       (i.e., the solution residuals) using the first-order
       upwind scheme with a variety of flux functions. */

    // x-direction fluxes.
    for ( j = 1 ; j <= Number_of_Cells_Ydir ; ++j ) {
       Soln[0][j].dUdt = Gaussian2D_U_VACUUM;
       for ( i = 0 ; i <= Number_of_Cells_Xdir ; ++i ) {
          Soln[i+1][j].dUdt = Gaussian2D_U_VACUUM;

          switch(Flux_Function_Type) {
            //case FLUX_FUNCTION_GODUNOV :
	    //Flux = FluxGodunov_x(Soln[i][j].W, Soln[i+1][j].W);
	    //break;
            case FLUX_FUNCTION_ROE :
              Flux = FluxRoe_x(Soln[i][j].W, Soln[i+1][j].W);
              break;
	      //case FLUX_FUNCTION_RUSANOV :
              //Flux = FluxRusanov_x(Soln[i][j].W, Soln[i+1][j].W);
              //break;
	    case FLUX_FUNCTION_HLLE :
              Flux = FluxHLLE_x(Soln[i][j].W, Soln[i+1][j].W);
              break;
	      //case FLUX_FUNCTION_LINDE :
              //Flux = FluxLinde_x(Soln[i][j].W, Soln[i+1][j].W);
              //break;
	      //case FLUX_FUNCTION_HLLC :
              //Flux = FluxHLLC_x(Soln[i][j].W, Soln[i+1][j].W);
              //break;
            default:
              Flux = FluxRoe_x(Soln[i][j].W, Soln[i+1][j].W);
              break;
          } /* endswitch */

          Soln[i][j].dUdt -= Flux/Soln[i][j].X.dx.x;
          Soln[i+1][j].dUdt += Flux/Soln[i+1][j].X.dx.x;
       } /* endfor */
       Soln[Number_of_Cells_Xdir+1][j].dUdt = Gaussian2D_U_VACUUM;
    } /* endfor */

    // y-direction fluxes.
    for ( i = 1 ; i <= Number_of_Cells_Xdir ; ++i ) {
       Soln[i][0].dUdt = Gaussian2D_U_VACUUM;
       for ( j = 0 ; j <= Number_of_Cells_Ydir ; ++j ) {
          switch(Flux_Function_Type) {
            //case FLUX_FUNCTION_GODUNOV :
	    //Flux = FluxGodunov_y(Soln[i][j].W, Soln[i][j+1].W);
	    //break;
            case FLUX_FUNCTION_ROE :
              Flux = FluxRoe_y(Soln[i][j].W, Soln[i][j+1].W);
              break;
	      //case FLUX_FUNCTION_RUSANOV :
              //Flux = FluxRusanov_y(Soln[i][j].W, Soln[i][j+1].W);
              //break;
	    case FLUX_FUNCTION_HLLE :
              Flux = FluxHLLE_y(Soln[i][j].W, Soln[i][j+1].W);
              break;
	      //case FLUX_FUNCTION_LINDE :
              //Flux = FluxLinde_y(Soln[i][j].W, Soln[i][j+1].W);
              //break;
	      //case FLUX_FUNCTION_HLLC :
              //Flux = FluxHLLC_y(Soln[i][j].W, Soln[i][j+1].W);
              //break;
            default:
              Flux = FluxRoe_y(Soln[i][j].W, Soln[i][j+1].W);
              break;
          } /* endswitch */

          Soln[i][j].dUdt -= Flux/Soln[i][j].X.dx.y;
          Soln[i][j+1].dUdt += Flux/Soln[i][j+1].X.dx.y;
       } /* endfor */
       Soln[i][Number_of_Cells_Ydir+1].dUdt = Gaussian2D_U_VACUUM;
    } /* endfor */

    /* Update both conserved and primitive solution
       variables using explicit Euler time integration method. */

    for ( j = 1 ; j <= Number_of_Cells_Ydir ; ++j ) {
       for ( i = 1 ; i <= Number_of_Cells_Xdir ; ++i ) {
          if ( !Local_Time_Stepping ) Soln[i][j].dt = dtMin;
	
          Soln[i][j].U += (CFL_Number*Soln[i][j].dt)*Soln[i][j].dUdt;
	
	  if (Soln[i][j].U.d   <= ZERO){ // ||
	      //Soln[i][j].U.E   <= ZERO ||
	      //Soln[i][j].U.e() <= ZERO ) {
	      cout << "\n " << CFDkit_Name() << " Gaussian2D ERROR: Negative Density and/or Energy: \n"
	           << " cell = (" << i << ", " << j << ") \n U = " 
                   << Soln[i][j].U << "\n dUdt = " << Soln[i][j].dUdt << "\n";
	      return (i);
	  }
	
	  Soln[i][j].W = W(Soln[i][j].U);
          Soln[i][j].W.relax(Soln[i][j].dt,1,W(Soln[i][j].Uo));
	  Soln[i][j].U = U(Soln[i][j].W);

       } /* endfor */    
    } /* endfor */

    /* Apply boundary conditions. */

    //x-direction
    for ( j = 1 ; j <= Number_of_Cells_Ydir ; ++j ) {
       Soln[0][j].W = 
	  BCs_x(W_left,
                Soln[1][j].W,
	        Gaussian2D_W_VACUUM,
	        Soln[1][j].X.dx.x,
	        BC_type_left,
	        LEFT_END_BOUNDARY);
       Soln[0][j].U = 
          U(Soln[0][j].W);

       Soln[Number_of_Cells_Xdir+1][j].W = 
          BCs_x(W_right,
                Soln[Number_of_Cells_Xdir][j].W,
		Gaussian2D_W_VACUUM,
		Soln[Number_of_Cells_Xdir][j].X.dx.x,
		BC_type_right,
		RIGHT_END_BOUNDARY);
       Soln[Number_of_Cells_Xdir+1][j].U = 
          U(Soln[Number_of_Cells_Xdir+1][j].W);
    } /* endfor */

    //y-direction
    for ( i = 1 ; i <= Number_of_Cells_Xdir ; ++i ) {
       Soln[i][0].W =  
	  BCs_y(W_bottom,
                Soln[i][1].W,
	        Gaussian2D_W_VACUUM,
	        Soln[i][1].X.dx.y,
	        BC_type_bottom,
	        LEFT_END_BOUNDARY);
       Soln[i][0].U = 
          U(Soln[i][0].W);

       Soln[i][Number_of_Cells_Ydir+1].W = 
          BCs_y(W_top,
                Soln[i][Number_of_Cells_Ydir].W,
		Gaussian2D_W_VACUUM,
		Soln[i][Number_of_Cells_Ydir].X.dx.x,
		BC_type_top,
		RIGHT_END_BOUNDARY); 
       Soln[i][Number_of_Cells_Ydir+1].U = 
          U(Soln[i][Number_of_Cells_Ydir+1].W);
    } /* endfor */

    /* Solution successfully updated. */

    return (0);
    
}

/********************************************************
 * Routine: dUdt_2stage_2ndOrder_upwind                 *
 *                                                      *
 * This routine updates the solution using a two-stage  *
 * second-order explicit time integration scheme        *
 * and a 2nd-ororder limited upwind spatial             *
 * discretization scheme with either the Godunov, Roe,  *
 * Rusanov, HLLE, Linde, or HLLC flux functions.        *
 *                                                      *
 ********************************************************/
int dUdt_2stage_2ndOrder_upwind(Gaussian2D_Cartesian_Cell **Soln,
   		                const int Number_of_Cells_Xdir,
                                const int Number_of_Cells_Ydir,
			        const int BC_type_left,
	      		        const int BC_type_right,
			        const int BC_type_bottom,
	      		        const int BC_type_top,
                                const Gaussian2D_pState W_left,
                                const Gaussian2D_pState W_right,
                                const Gaussian2D_pState W_bottom,
                                const Gaussian2D_pState W_top,
			        double &dtMin,
			        const double &CFL_Number,
                                const int Reconstruction_Type,
                                const int Limiter_Type,
			        const int Flux_Function_Type,
			        const int Local_Time_Stepping) {

    int i, j, n_stage;
    double omega;
    Gaussian2D_pState Wl, Wr;
    Gaussian2D_cState Flux;

    /* Perform second-order two-stage semi-implicit update of solution
       variables for new time level. */

    for ( n_stage = 1 ; n_stage <= 2 ; ++n_stage ) {

        /* Evaluate the time step fraction for stage. */

        omega = ONE/double(n_stage);

        /* Apply boundary conditions for stage. */

        //x-direction BCs.
        for ( j = 1 ; j <= Number_of_Cells_Ydir ; ++j ) {
           Soln[0][j].W = 
	      BCs_x(W_left,
                    Soln[1][j].W,
	            Gaussian2D_W_VACUUM,
    	            Soln[1][j].X.dx.x,
	            BC_type_left,
	            LEFT_END_BOUNDARY);
           Soln[0][j].U = 
              U(Soln[0][j].W);

           Soln[Number_of_Cells_Xdir+1][j].W = 
              BCs_x(W_right,
                    Soln[Number_of_Cells_Xdir][j].W,
		    Gaussian2D_W_VACUUM,
		    Soln[Number_of_Cells_Xdir][j].X.dx.x,
		    BC_type_right,
		    RIGHT_END_BOUNDARY);
           Soln[Number_of_Cells_Xdir+1][j].U = 
              U(Soln[Number_of_Cells_Xdir+1][j].W);
        } /* endfor */

        //y-direction BCs.
        for ( i = 1 ; i <= Number_of_Cells_Xdir ; ++i ) {
           Soln[i][0].W =  
	      BCs_y(W_bottom,
                    Soln[i][1].W,
	            Gaussian2D_W_VACUUM,
	            Soln[i][1].X.dx.y,
	            BC_type_bottom,
	            LEFT_END_BOUNDARY);
           Soln[i][0].U = 
              U(Soln[i][0].W);

           Soln[i][Number_of_Cells_Ydir+1].W = 
              BCs_y(W_top,
                    Soln[i][Number_of_Cells_Ydir].W,
 		    Gaussian2D_W_VACUUM,
		    Soln[i][Number_of_Cells_Ydir].X.dx.x,
		    BC_type_top,
		    RIGHT_END_BOUNDARY); 
           Soln[i][Number_of_Cells_Ydir+1].U = 
              U(Soln[i][Number_of_Cells_Ydir+1].W);
        } /* endfor */

        /* Perform the linear reconstruction within each cell
           of the computational grid in first stage. */

        if ( n_stage == 1 ) {
          switch(Reconstruction_Type) {
            case RECONSTRUCTION_GREEN_GAUSS :
               Linear_Reconstruction_GreenGauss(Soln, 
                                                Number_of_Cells_Xdir,
                                                Number_of_Cells_Ydir,
                                                Limiter_Type);
              break;
            case RECONSTRUCTION_LEAST_SQUARES :
               Linear_Reconstruction_LeastSquares(Soln, 
                                                  Number_of_Cells_Xdir,
                                                  Number_of_Cells_Ydir,
                                                  Limiter_Type);
              break;
	    default:
               Linear_Reconstruction_GreenGauss(Soln, 
                                                Number_of_Cells_Xdir,
                                                Number_of_Cells_Ydir,
                                                Limiter_Type);
              break;
          } /* endswitch */
        } /* endif */

        /* Evaluate the time rate of change of the solution
           (i.e., the solution residuals) using a second-order
           limited upwind scheme with a variety of flux functions. */

        // Add x-direction fluxes.
        for ( j = 0 ; j <= Number_of_Cells_Ydir+1 ; ++j ) {
           if ( !Local_Time_Stepping && n_stage == 1 ) 
              Soln[0][j].dt = dtMin;
           if ( n_stage == 1 ) Soln[0][j].Uo = Soln[0][j].U;
           Soln[0][j].dUdt = Gaussian2D_U_VACUUM;

           for ( i = 0 ; i <= Number_of_Cells_Xdir ; ++i ) {
              if ( !Local_Time_Stepping && n_stage == 1 ) 
                  Soln[i+1][j].dt = dtMin;
              if ( n_stage == 1 ) {
                  Soln[i+1][j].Uo = Soln[i+1][j].U;
                  Soln[i+1][j].dUdt = Gaussian2D_U_VACUUM;
              } else if (j > 0 && j <= Number_of_Cells_Ydir ) {
                  Soln[i+1][j].dUdt = HALF*Soln[i+1][j].dUdt;
              } /* endif */

              if (j > 0 && j <= Number_of_Cells_Ydir ) {

                 /* Evaluate the cell interface x-direction fluxes. */

                 Wl = Soln[i][j].W + 
                      HALF*(Soln[i][j].phi^Soln[i][j].dWdx)*
                      Soln[i][j].X.dx.x;
                 Wr = Soln[i+1][j].W - 
                      HALF*(Soln[i+1][j].phi^Soln[i+1][j].dWdx)*
                      Soln[i+1][j].X.dx.x;

                 switch(Flux_Function_Type) {
                   //case FLUX_FUNCTION_GODUNOV :
		   //Flux = FluxGodunov_x(Wl, Wr);
		   //break;
                   case FLUX_FUNCTION_ROE :
                     Flux = FluxRoe_x(Wl, Wr);
                     break;
		     //case FLUX_FUNCTION_RUSANOV :
                     //Flux = FluxRusanov_x(Wl, Wr);
                     //break;
		   case FLUX_FUNCTION_HLLE :
                     Flux = FluxHLLE_x(Wl, Wr);
                     break;
		     //case FLUX_FUNCTION_LINDE :
                     //Flux = FluxLinde_x(Wl, Wr);
                     //break;
		     //case FLUX_FUNCTION_HLLC :
                     //Flux = FluxHLLC_x(Wl, Wr);
                     //break;
                   default:
                     Flux = FluxRoe_x(Wl, Wr);
                     break;
                 } /* endswitch */

                 /* Evaluate cell-averaged solution changes. */

                 Soln[i][j].dUdt -= (omega*CFL_Number*Soln[i][j].dt)*
                                    Flux/Soln[i][j].X.dx.x;
                 Soln[i+1][j].dUdt += (omega*CFL_Number*Soln[i+1][j].dt)*
                                      Flux/Soln[i+1][j].X.dx.x;
              } /* endif */
           } /* endfor */

           if (j > 0 && j <= Number_of_Cells_Ydir )
              Soln[Number_of_Cells_Xdir+1][j].dUdt = Gaussian2D_U_VACUUM;
        } /* endfor */

        // Add y-direction fluxes.
        for ( i = 1 ; i <= Number_of_Cells_Xdir ; ++i ) {
           for ( j = 0 ; j <= Number_of_Cells_Ydir ; ++j ) {

              /* Evaluate the cell interface y-direction fluxes. */
             
              Wl = Soln[i][j].W + 
                   HALF*(Soln[i][j].phi^Soln[i][j].dWdy)*
                   Soln[i][j].X.dx.y;
              Wr = Soln[i][j+1].W - 
                   HALF*(Soln[i][j+1].phi^Soln[i][j+1].dWdy)*
                   Soln[i][j+1].X.dx.y;

              switch(Flux_Function_Type) {
                //case FLUX_FUNCTION_GODUNOV :
		//Flux = FluxGodunov_y(Wl, Wr);
		//break;
                case FLUX_FUNCTION_ROE :
                  Flux = FluxRoe_y(Wl, Wr);
                  break;
		  //case FLUX_FUNCTION_RUSANOV :
                  //Flux = FluxRusanov_y(Wl, Wr);
                  //break;
		case FLUX_FUNCTION_HLLE :
                  Flux = FluxHLLE_y(Wl, Wr);
                  break;
		  //case FLUX_FUNCTION_LINDE :
                  //Flux = FluxLinde_y(Wl, Wr);
                  //break;
		  //case FLUX_FUNCTION_HLLC :
                  //Flux = FluxHLLC_y(Wl, Wr);
                  //break;
                default:
                  Flux = FluxRoe_y(Wl, Wr);
                  break;
              } /* endswitch */

              /* Evaluate cell-averaged solution changes. */

              Soln[i][j].dUdt -= (omega*CFL_Number*Soln[i][j].dt)*
                                 Flux/Soln[i][j].X.dx.y;
              Soln[i][j+1].dUdt += (omega*CFL_Number*Soln[i][j+1].dt)*
                                   Flux/Soln[i][j+1].X.dx.y;
           } /* endfor */

           Soln[i][Number_of_Cells_Ydir+1].dUdt = Gaussian2D_U_VACUUM;
        } /* endfor */

        /* Update solution variables for this stage. */

        for ( j = 1 ; j <= Number_of_Cells_Ydir ; ++j ) {
           for ( i = 1 ; i <= Number_of_Cells_Xdir ; ++i ) {
              Soln[i][j].U = Soln[i][j].Uo + Soln[i][j].dUdt;
	 
  	      if (Soln[i][j].U.d   <= ZERO){// ||
		//Soln[i][j].U.E   <= ZERO ||
		//Soln[i][j].U.e() <= ZERO ) {
		  cout << "\n " << CFDkit_Name() << " Gaussian2D ERROR: Negative Density and/or Energy: \n"
	               << " cell = (" << i << ", " << j << ") \n U = " 
                       << Soln[i][j].U << "\n dUdt = " << Soln[i][j].dUdt << "\n";
	          return (i);
	      }
	
	      Soln[i][j].W = W(Soln[i][j].U);
	      Soln[i][j].W.relax(CFL_Number*Soln[i][j].dt,n_stage,W(Soln[i][j].Uo));
	      Soln[i][j].U = U(Soln[i][j].W);
           } /* endfor */    
        } /* endfor */

    } /* endfor */

    /* Solution successfully updated. */

    return (0);
    
}

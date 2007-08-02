/* HyperHeat1DState.cc:  Subroutines for 1D Hyperbolic Heat Equation
                         (Maxwell-Cattaneo Equation) Solution State Class. */

/* Include 1D hyperbolic heat equations solution state header file. */

#ifndef _HYPERHEAT1D_STATE_INCLUDED
#include "HyperHeat1DState.h"
#endif // _HYPERHEAT1D_STATE_INCLUDED

/*************************************************************
 * HyperHeat1D_State -- Create storage and assign            *
 *                      conductivity and relaxation time.    *
 *************************************************************/
double HyperHeat1D_State::kappa = ONE;
double HyperHeat1D_State::tau = ONE;

/********************************************************
 * Routine: RiemannHomo (Exact Riemann solver,          *
 *                       no sources)                    *
 *                                                      *
 * This function obtains the exact solution to the      *
 * Riemann problem for the homogeneous form (no         *
 * sources) of the 1D hyperbolic heat equations         *
 * (Maxwell-Cattaneo equations), returning the          *
 * intermediate state variables along the ray x/t=0.    *
 *                                                      *
 ********************************************************/
HyperHeat1D_State RiemannHomo(const HyperHeat1D_State &Ul,
	      	              const HyperHeat1D_State &Ur) {

    double a, vl, vr;
    double T_m, qx_m;

    /* Determine the wave speeds of the left and right
       moving waves. */

    a = Ul.a();
    vl = Ul.lambda(1);
    vr = Ul.lambda(2); 

    /* Return the intermediate state solution. */

    T_m   = HALF*(Ur.T+Ul.T)-HALF*(Ur.qx-Ul.qx)/a;
    qx_m  = HALF*(Ur.qx+Ul.qx)-HALF*a*(Ur.T-Ul.T);

    return (HyperHeat1D_State(T_m, qx_m));
       
}

/********************************************************
 * Routine: BCs (Boundary Conditions)                   *
 *                                                      *
 * This function returns the boundary state given       *
 * the boundary and interior states, the boundary type, *
 * and grid end type.                                   *
 *                                                      *
 ********************************************************/
HyperHeat1D_State BCs(const HyperHeat1D_State &Ub,
                      const HyperHeat1D_State &Ui,
		      const HyperHeat1D_State &dUdx,
		      const double &dx,
		      int BC_type,
		      int End_type) {
    switch(BC_type) {

      case BC_FIXED :
        switch(End_type) {
          case LEFT_END_BOUNDARY :
            return (Ub-HALF*dUdx*dx);
            break;
          case RIGHT_END_BOUNDARY :
            return (Ub+HALF*dUdx*dx);
            break;
          default:
            return (Ub+HALF*dUdx*dx);
            break;
        } /* endswitch */
        break;

      case BC_CONSTANT_EXTRAPOLATION :
        return (Ui);
        break;

      case BC_LINEAR_EXTRAPOLATION :
        switch(End_type) {
          case LEFT_END_BOUNDARY :
            return (Ui-dUdx*dx);
            break;
          case RIGHT_END_BOUNDARY :
            return (Ui+dUdx*dx);
            break;
          default:
            return (Ui+dUdx*dx);
            break;
        } /* endswitch */
        break;

      case BC_REFLECTION :
        return (HyperHeat1D_State(Ui.T,-Ui.qx));
        break;

      case BC_PERIODIC :
        return (Ub);
        break;

      case BC_FIXED_TEMP:
        switch(End_type) {
          case LEFT_END_BOUNDARY :
	    return (HyperHeat1D_State(Ub.T-HALF*dUdx.T*dx,
                                      Ui.qx-dUdx.qx*dx+
                                      (Ub.T-Ui.T+HALF*dUdx.T*dx)*Ui.a()));
            break;
          case RIGHT_END_BOUNDARY :
	    return (HyperHeat1D_State(Ub.T+HALF*dUdx.T*dx,
                                      Ui.qx+dUdx.qx*dx-
                                      (Ub.T-Ui.T-HALF*dUdx.T*dx)*Ui.a()));
            break;
          default:
	    return (HyperHeat1D_State(Ub.T+HALF*dUdx.T*dx,
                                      Ui.qx+dUdx.qx*dx-
                                      (Ub.T-Ui.T-HALF*dUdx.T*dx)*Ui.a()));
            break;
        } /* endswitch */
        break;

      case BC_FIXED_HEATFLUX:
        switch(End_type) {
          case LEFT_END_BOUNDARY :
            return (HyperHeat1D_State(Ui.T-dUdx.T*dx+
                                      (Ub.qx-Ui.qx+HALF*dUdx.qx*dx)/Ui.a(),
                                      Ub.qx-HALF*dUdx.qx*dx));
            break;
          case RIGHT_END_BOUNDARY :
            return (HyperHeat1D_State(Ui.T+dUdx.T*dx-
                                      (Ub.qx-Ui.qx-HALF*dUdx.qx*dx)/Ui.a(),
                                      Ub.qx+HALF*dUdx.qx*dx));
            break;
          default:
            return (HyperHeat1D_State(Ui.T+dUdx.T*dx-
                                      (Ub.qx-Ui.qx-HALF*dUdx.qx*dx)/Ui.a(),
                                      Ub.qx+HALF*dUdx.qx*dx));
            break;
        } /* endswitch */
        break;

      default:
        return (Ub);
        break;

    } /* endswitch */
}

/********************************************************
 * Routine: WaveSpeedPos                                *
 *                                                      *
 * This function returns the positive parts of the      *
 * elemental wave speeds or eigenvalues.                *
 *                                                      *
 ********************************************************/
HyperHeat1D_State WaveSpeedPos(const HyperHeat1D_State &lambdas) {
  return (HyperHeat1D_State(HALF*(lambdas[1]+fabs(lambdas[1])),
                            HALF*(lambdas[2]+fabs(lambdas[2]))));
}

/********************************************************
 * Routine: WaveSpeedNeg                                *
 *                                                      *
 * This function returns the negative parts of the      *
 * elemental wave speeds or eigenvalues.                *
 *                                                      *
 ********************************************************/
HyperHeat1D_State WaveSpeedNeg(const HyperHeat1D_State &lambdas) {
  return (HyperHeat1D_State(HALF*(lambdas[1]-fabs(lambdas[1])),
                            HALF*(lambdas[2]-fabs(lambdas[2]))));
}

/********************************************************
 * Routine: WaveSpeedAbs                                *
 *                                                      *
 * This function returns the absolute values of the     *
 * elemental wave speeds or eigenvalues.                *
 *                                                      *
 ********************************************************/
HyperHeat1D_State WaveSpeedAbs(const HyperHeat1D_State &lambdas) {
  return (HyperHeat1D_State(fabs(lambdas[1]),
                            fabs(lambdas[2])));
}

/*********************************************************
 * Routine: FluxGodunov (Godunov's flux function)        *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux given left and right solution states by solving  *
 * exactly the Riemann problem for the homogeneous form  *
 * of the equations associated with the two states.  See *
 * Godunov (1959).                                       *
 *                                                       *
 *********************************************************/
HyperHeat1D_State FluxGodunov(const HyperHeat1D_State &Ul,
	      	              const HyperHeat1D_State &Ur) {
    return (F(RiemannHomo(Ul, Ur)));
}

/*********************************************************
 * Routine: FluxRoe (Roe's flux function)                *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux given left and right solution states by using    *
 * the approximate Riemann solver of Roe for the two     *
 * states.  See Roe (1981).                              *
 *                                                       *
 *********************************************************/
HyperHeat1D_State FluxRoe(const HyperHeat1D_State &Ul,
	      	          const HyperHeat1D_State &Ur) {

    int i;
    HyperHeat1D_State Flux, dUrl, lambdas, wavespeeds;

    /* Evaluate the jumps in the solution states. */

    dUrl = Ur-Ul;

    /* Evaluate the eigenvalues and absolute wave speeds. */

    lambdas = Ul.lambda();
    wavespeeds = WaveSpeedAbs(lambdas);

    /* Determine the intermediate state flux. */

    Flux = HALF*(Ul.F()+Ur.F());

    for ( i = 1 ; i <= NUM_VAR_HYPERHEAT1D ; ++i ) {
      Flux -= HALF*wavespeeds[i]*(Ul.l(i)*dUrl)*Ul.r(i);
    } /* endfor */

    /* Return solution flux. */
    
    return (Flux);
    
}

/*********************************************************
 * Routine: FluxRusanov (Rusanov flux function)          *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux given left and right solution states by using    *
 * the Rusanov approximation for the fluxes. See         *
 * Rusanov (1964).                                       *
 *                                                       *
 *********************************************************/
HyperHeat1D_State FluxRusanov(const HyperHeat1D_State &Ul,
	      	              const HyperHeat1D_State &Ur) {

    HyperHeat1D_State Flux;

    /* Return solution flux. */

    return (Flux);

}

/*********************************************************
 * Routine: FluxHLLE (Harten-Lax-van Leer flux function) *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux given left and right solution states by using    *
 * the so-called Harten-Lax-van Leer approximation for   *
 * for fluxes.   See Harten, Lax, van Leer (1983).       *
 *                                                       *
 *********************************************************/
HyperHeat1D_State FluxHLLE(const HyperHeat1D_State &Ul,
	      	           const HyperHeat1D_State &Ur) {

    HyperHeat1D_State Flux;
    
    /* Return solution flux. */

    return (Flux);

}

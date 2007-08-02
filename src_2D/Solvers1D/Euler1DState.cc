/* Euler1DState.cc:  Subroutines for 1D Euler Solution State Classes. */

/* Include 1D Euler solution state header file. */

#ifndef _EULER1D_STATE_INCLUDED
#include "Euler1DState.h"
#endif // _EULER1D_STATE_INCLUDED

/*************************************************************
 * Euler1D_pState -- Create storage and assign gas constants.*
 *************************************************************/
double Euler1D_pState::g = GAMMA_AIR;
double Euler1D_pState::gm1 = GAMMA_AIR-ONE;
double Euler1D_pState::gm1i = ONE/(GAMMA_AIR-ONE);
double Euler1D_pState::R = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);

/*************************************************************
 * Euler1D_cState -- Create storage and assign gas constants.*
 *************************************************************/
double Euler1D_cState::g = GAMMA_AIR;
double Euler1D_cState::gm1 = GAMMA_AIR-ONE;
double Euler1D_cState::gm1i = ONE/(GAMMA_AIR-ONE);
double Euler1D_cState::R = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);

/********************************************************
 * Routine: Riemann (Exact Riemann solver)              *
 *                                                      *
 * This function uses a Newton-Raphson interative       *
 * procedure to obtain the exact solution to the        *
 * Riemann problem for the 1D Euler equations,          *
 * returning the intermediate state variables along     *
 * the ray x/t=0.  See Gottlieb and Groth (1987).       *
 *                                                      *
 ********************************************************/
Euler1D_pState Riemann(const Euler1D_pState &Wl,
	      	       const Euler1D_pState &Wr) {

    int number_of_iterations;
    
    double al, ar, CL, CR, Z;
    double dml, dmr, vm, pm, aml, amr;
    double msl, pml, dpmldum, msr, pmr, dpmrdum;
    double vsl, vhl, vtl, vsr, vhr, vtr;

    /* Determine the left and right state sound speeds. */

    al = Wl.a();
    ar = Wr.a();

    /* Compute the left and right state Riemann invariants. */

    CL=Wl.v+TWO*al/Wl.gm1;
    CR=Wr.v-TWO*ar/Wr.gm1;

    /* Check for vacuum state. */

    if ( CL-CR <= ZERO ) {
        return (Euler1D_W_VACUUM);
    } /* endif */

    /* Make an initial estimate of the intermediate state flow
       velocity to begin the Newton-Raphson iterative solution
       procedure.  The initial guess tate velocity is made
       based on isentropic flow theory. */

    Z = (ar/al)*pow((Wl.p/Wr.p), HALF*Wl.gm1/Wl.g);
    vm = (CL*Z+CR)/(ONE+Z);

    /* In the case that two rarefaction waves are present,
       then an exact solution has been found and the iterative
       procedure is not required.  Check for this. */

    if ( vm >= Wl.v && vm <= Wr.v ) {
        if (vm >= ZERO) {
            aml = al-HALF*Wl.gm1*(vm-Wl.v);
	    pm = Wl.p*pow((aml/al), TWO*Wl.g/Wl.gm1);
            vhl = Wl.v-al;
	    vtl = vm-aml;
            if (vhl >= ZERO) {
                return (Euler1D_pState(Wl));
            } else if (vtl <= ZERO) {
 	        dml = Wl.g*pm/sqr(aml);
                return (Euler1D_pState(dml, vm, pm));
            } else {
	        vm = (Wl.gm1*Wl.v+TWO*al)/(Wl.g+ONE);
                pm = Wl.p*pow((vm/al), TWO*Wl.g/Wl.gm1);
		dml = Wl.g*pm/sqr(vm);
                return (Euler1D_pState(dml, vm, pm));
            } /* endif */
        } else {
            amr = ar+HALF*Wr.gm1*(vm-Wr.v);
	    pm = Wr.p*pow((amr/ar), TWO*Wr.g/Wr.gm1);
            vhr = Wr.v+ar;
	    vtr = vm+amr;
            if (vhr <= ZERO) {
                return (Euler1D_pState(Wr));
            } else if (vtr >= ZERO) {
 	        dmr = Wr.g*pm/sqr(amr);
                return (Euler1D_pState(dmr, vm, pm));
            } else {
	        vm = (Wr.gm1*Wr.v-TWO*ar)/(Wr.g+ONE);
                pm = Wr.p*pow((-vm/ar), TWO*Wr.g/Wr.gm1);
 	        dmr = Wr.g*pm/sqr(vm);
                return (Euler1D_pState(dmr, vm, pm));
            } /* endif */
        } /* end if */
    } /* endif */
    
    /* Perform the Newton-Raphson iterative procedure and solve for
       the velocity in the intermediate state.  During this iterative
       process the pressure in the intermediate state is also found. */

    number_of_iterations = 0;
    
    while (1) {
        /* Update the iteration counter. */
        number_of_iterations = number_of_iterations + 1;
 
        /* Determine solution changes for left wave. */
        if ( vm < Wl.v ) {
            msl = (Wl.g+ONE)*(vm-Wl.v)/(FOUR*al);
            msl = msl-sqrt(ONE+sqr(msl));
            pml = Wl.p*(ONE+Wl.g*(vm-Wl.v)*msl/al);
            dpmldum = TWO*Wl.g*Wl.p*cube(msl)/(al*(ONE+sqr(msl)));
	} else {
            aml = al-HALF*Wl.gm1*(vm-Wl.v);
            pml = Wl.p*pow((aml/al), TWO*Wl.g/Wl.gm1);
	    dpmldum = -Wl.g*pml/aml;
        } /* end if */
	  
        /* Determine solution changes for right wave. */
        if ( vm > Wr.v ) {
            msr = (Wr.g+ONE)*(vm-Wr.v)/(FOUR*ar);
            msr = msr+sqrt(ONE+sqr(msr));
            pmr = Wr.p*(ONE+Wr.g*(vm-Wr.v)*msr/ar);
            dpmrdum = TWO*Wr.g*Wr.p*cube(msr)/(ar*(ONE+sqr(msr)));
	} else {
            amr = ar+HALF*Wr.gm1*(vm-Wr.v);
            pmr = Wr.p*pow((amr/ar), TWO*Wr.g/Wr.gm1);
	    dpmrdum = Wr.g*pmr/amr;
        } /* end if */

	/* Check for convergence (i.e., pml=pmr). */

	if ( fabs(ONE-pml/pmr) <= TOLER) break;

	/* Compute next estimate for the intermediate
	   state velocity, vm. */

	vm = vm-(pml-pmr)/(dpmldum-dpmrdum);
	
    } /* endwhile */

    pm = HALF*(pml+pmr);

    /* Return the intermediate state solution. */

    if ( vm >= ZERO ) {
        if ( vm < Wl.v ) {
            aml = al * sqrt(((Wl.g+ONE)+Wl.gm1*pm/Wl.p) /
			    ((Wl.g+ONE)+Wl.gm1*Wl.p/pm) );
            vsl = Wl.v+msl*al;
            if (vsl >= ZERO) {
                return (Euler1D_pState(Wl));                
            } else {
 	        dml = Wl.g*pm/sqr(aml);
                return (Euler1D_pState(dml, vm, pm));
            } /* endif */
        } else {
            vhl = Wl.v-al;
	    vtl = vm-aml;
            if (vhl >= ZERO) {
                return (Euler1D_pState(Wl));
            } else if (vtl <= ZERO) {
 	        dml = Wl.g*pm/sqr(aml);
                return (Euler1D_pState(dml, vm, pm));
            } else {
	        vm = (Wl.gm1*Wl.v+TWO*al)/(Wl.g+ONE);
                pm = Wl.p*pow((vm/al), TWO*Wl.g/Wl.gm1);
		dml = Wl.g*pm/sqr(vm);
                return (Euler1D_pState(dml, vm, pm));
            } /* endif */
        } /* endif */
    } else {
        if ( vm > Wr.v ) {
            amr = ar * sqrt(((Wr.g+ONE)+Wr.gm1*pm/Wr.p) /
			    ((Wr.g+ONE)+Wr.gm1*Wr.p/pm) );
            vsr = Wr.v+msr*ar;
            if (vsr <= ZERO) {
                return (Euler1D_pState(Wr));                
            } else {
 	        dmr = Wr.g*pm/sqr(amr);
                return (Euler1D_pState(dmr, vm, pm));
            } /* endif */
        } else {
            vhr = Wr.v+ar;
	    vtr = vm+amr;
            if (vhr <= ZERO) {
                return (Euler1D_pState(Wr));
            } else if (vtr >= ZERO) {
 	        dmr = Wr.g*pm/sqr(amr);
                return (Euler1D_pState(dmr, vm, pm));
            } else {
	        vm = (Wr.gm1*Wr.v-TWO*ar)/(Wr.g+ONE);
                pm = Wr.p*pow((-vm/ar), TWO*Wr.g/Wr.gm1);
 	        dmr = Wr.g*pm/sqr(vm);
                return (Euler1D_pState(dmr, vm, pm));
            } /* endif */
        } /* endif */
    } /* end if */
       
}

/********************************************************
 * Routine: RoeAverage (Roe Averages)                   *
 *                                                      *
 * This function returns the Roe-averaged primitive     *
 * solution state given left and right primitive        *
 * solution variables.  See Roe (1981).                 *
 *                                                      *
 ********************************************************/
Euler1D_pState RoeAverage(const Euler1D_pState &Wl,
	      	          const Euler1D_pState &Wr) {

    double hl, hr, sdl, sdr;
    double da, va, pa, aa2, ha, ga, gam1;

    /* Determine the left and right state specific enthalpies
       and square roots of the density. */

    hl = Wl.h();
    hr = Wr.h();
    sdl = sqrt(Wl.d);
    sdr = sqrt(Wr.d);

    /* Determine the appropriate Roe averages. */

    ga = Wl.g;
    gam1 = Wl.gm1;
    da = sdl*sdr;
    va = (sdl*Wl.v+sdr*Wr.v)/(sdl+sdr);
    ha = (sdl*hl+sdr*hr)/(sdl+sdr);
    aa2 = gam1*(ha-HALF*sqr(va));
    pa = da*aa2/ga;

    /* Return the Roe-averged state. */

    return (Euler1D_pState(da, va, pa));
       
}

/********************************************************
 * Routine: WaveSpeedPos                                *
 *                                                      *
 * This function returns the positive parts of the      *
 * elemental wave speeds or eigenvalues.                *
 *                                                      *
 ********************************************************/
Euler1D_pState WaveSpeedPos(const Euler1D_pState &lambdas_a,
                            const Euler1D_pState &lambdas_l,
                            const Euler1D_pState &lambdas_r) {
  return (Euler1D_pState(HALF*(lambdas_a[1]+fabs(lambdas_a[1])),
                         HALF*(lambdas_a[2]+fabs(lambdas_a[2])),
                         HALF*(lambdas_a[3]+fabs(lambdas_a[3]))));
}

/********************************************************
 * Routine: WaveSpeedNeg                                *
 *                                                      *
 * This function returns the negative parts of the      *
 * elemental wave speeds or eigenvalues.                *
 *                                                      *
 ********************************************************/
Euler1D_pState WaveSpeedNeg(const Euler1D_pState &lambdas_a,
                            const Euler1D_pState &lambdas_l,
                            const Euler1D_pState &lambdas_r) {
  return (Euler1D_pState(HALF*(lambdas_a[1]-fabs(lambdas_a[1])),
                         HALF*(lambdas_a[2]-fabs(lambdas_a[2])),
                         HALF*(lambdas_a[3]-fabs(lambdas_a[3]))));
}

/********************************************************
 * Routine: WaveSpeedAbs                                *
 *                                                      *
 * This function returns the absolute values of the     *
 * elemental wave speeds or eigenvalues.                *
 *                                                      *
 ********************************************************/
Euler1D_pState WaveSpeedAbs(const Euler1D_pState &lambdas_a,
                            const Euler1D_pState &lambdas_l,
                            const Euler1D_pState &lambdas_r) {
  return (Euler1D_pState(fabs(lambdas_a[1]),
                         fabs(lambdas_a[2]),
                         fabs(lambdas_a[3])));
}

/********************************************************
 * Routine: HartenFixPos (Harten Entropy Fix)           *
 *                                                      *
 * This function returns the positive parts of the      *
 * corrected elemental wave speeds or eigenvalues       *
 * according to the entropy fix of Harten (1983).       *
 *                                                      *
 ********************************************************/
Euler1D_pState HartenFixPos(const Euler1D_pState &lambdas_a,
                            const Euler1D_pState &lambdas_l,
                            const Euler1D_pState &lambdas_r) {
  return (Euler1D_pState(HartenFixPos(lambdas_a[1],
                                      lambdas_l[1],
                                      lambdas_r[1]),
                         HALF*(lambdas_a[2]+fabs(lambdas_a[2])),
                         HartenFixPos(lambdas_a[3],
                                      lambdas_l[3],
                                      lambdas_r[3])));
}

/********************************************************
 * Routine: HartenFixNeg (Harten Entropy Fix)           *
 *                                                      *
 * This function returns the negative parts of the      *
 * corrected elemental wave speeds or eigenvalues       *
 * according to the entropy fix of Harten (1983).       *
 *                                                      *
 ********************************************************/
Euler1D_pState HartenFixNeg(const Euler1D_pState &lambdas_a,
                            const Euler1D_pState &lambdas_l,
                            const Euler1D_pState &lambdas_r) {
  return (Euler1D_pState(HartenFixNeg(lambdas_a[1],
                                      lambdas_l[1],
                                      lambdas_r[1]),
                         HALF*(lambdas_a[2]-fabs(lambdas_a[2])),
                         HartenFixNeg(lambdas_a[3],
                                      lambdas_l[3],
                                      lambdas_r[3])));
}

/********************************************************
 * Routine: HartenFixAbs (Harten Entropy Fix)           *
 *                                                      *
 * This function returns the absolute values of the     *
 * corrected elemental wave speeds or eigenvalues       *
 * according to the entropy fix of Harten (1983).       *
 *                                                      *
 ********************************************************/
Euler1D_pState HartenFixAbs(const Euler1D_pState &lambdas_a,
                            const Euler1D_pState &lambdas_l,
                            const Euler1D_pState &lambdas_r) {
  return (Euler1D_pState(HartenFixAbs(lambdas_a[1],
                                      lambdas_l[1],
                                      lambdas_r[1]),
                         fabs(lambdas_a[2]),
                         HartenFixAbs(lambdas_a[3],
                                      lambdas_l[3],
                                      lambdas_r[3])));
}

/*********************************************************
 * Routine: FluxGodunov (Godunov's flux function)        *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux given left and right solution states by solving  *
 * exactly the Riemann problem associated with the two   *
 * states.  See Godunov (1959).                          *
 *                                                       *
 *********************************************************/
Euler1D_cState FluxGodunov(const Euler1D_pState &Wl,
	      	           const Euler1D_pState &Wr) {
    return (F(Riemann(Wl, Wr)));
}

Euler1D_cState FluxGodunov(const Euler1D_cState &Ul,
	      	           const Euler1D_cState &Ur) {
    return (F(Riemann(Ul.W(), Ur.W())));
}

/*********************************************************
 * Routine: FluxRoe (Roe's flux function)                *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux given left and right solution states by using    *
 * the "linearized" approximate Riemann solver of Roe    *
 * for the two states.  See Roe (1981).                  *
 *                                                       *
 *********************************************************/
Euler1D_cState FluxRoe(const Euler1D_pState &Wl,
	      	       const Euler1D_pState &Wr) {

    int i;
    Euler1D_pState Wa, dWrl, 
                   wavespeeds, lambdas_l, lambdas_r, lambdas_a;
    Euler1D_cState Flux;

    /* Evaluate the Roe-average primitive solution state. */
    
    Wa = RoeAverage(Wl, Wr);

    /* Evaluate the jumps in the primitive solution states. */

    dWrl = Wr-Wl;

    /* Evaluate the left, right, and average state eigenvalues. */

    lambdas_l = Wl.lambda();
    lambdas_r = Wr.lambda();
    lambdas_a = Wa.lambda();

    /* Determine the intermediate state flux. */

    if (Wa.v >= ZERO) {
        Flux = Wl.F();
        wavespeeds = HartenFixNeg(lambdas_a,
                                  lambdas_l,
                                  lambdas_r);
        for ( i = 1 ; i <= NUM_VAR_EULER1D ; ++i ) {
           if (wavespeeds[i] < ZERO)
              Flux += wavespeeds[i]*(Wa.lp(i)*dWrl)*Wa.rc(i);
        } /* endfor */
    } else {
        Flux = Wr.F();
        wavespeeds = HartenFixPos(lambdas_a,
                                  lambdas_l,
                                  lambdas_r);
        for ( i = 1 ; i <= NUM_VAR_EULER1D ; ++i ) {
           if (wavespeeds[i] > ZERO)
              Flux -= wavespeeds[i]*(Wa.lp(i)*dWrl)*Wa.rc(i);
        } /* endfor */
    } /* endif */

    /* Return solution flux. */
    
    return (Flux);
    
}

Euler1D_cState FluxRoe(const Euler1D_cState &Ul,
	      	       const Euler1D_cState &Ur) {
   return (FluxRoe(Ul.W(), Ur.W()));
}  

/*********************************************************
 * Routine: FluxRusanov (Rusanov flux function)          *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux given left and right solution states by using    *
 * the Rusanov approximation for the fluxes.  See        *
 * Rusanov (1964).                                       *
 *                                                       *
 *********************************************************/
Euler1D_cState FluxRusanov(const Euler1D_pState &Wl,
	      	           const Euler1D_pState &Wr) {

    int i;
    double wavespeed_max;
    Euler1D_pState Wa, wavespeeds, lambdas_l, lambdas_r, lambdas_a;
    Euler1D_cState Flux, dUrl;
    
    /* Evaluate the Roe-average primitive solution state. */
    
    Wa = RoeAverage(Wl, Wr);

    /* Evaluate the jumps in the conserved solution states. */

    dUrl = Wr.U()-Wl.U();

    /* Evaluate the left, right, and average state eigenvalues. */

    lambdas_l = Wl.lambda();
    lambdas_r = Wr.lambda();
    lambdas_a = Wa.lambda();

    /* Determine the intermediate state flux. */

    Flux = HALF*(Wl.F()+Wr.F());
    wavespeeds = HartenFixAbs(lambdas_a,
                              lambdas_l,
                              lambdas_r);

    wavespeed_max = wavespeeds[1];
    for ( i = 2 ; i <= NUM_VAR_EULER1D ; ++i ) {
        wavespeed_max = max(wavespeed_max, wavespeeds[i]);
    } /* endfor */

    Flux -= HALF*wavespeed_max*dUrl;

    /* Return solution flux. */

    return (Flux);

}

Euler1D_cState FluxRusanov(const Euler1D_cState &Ul,
	      	           const Euler1D_cState &Ur) {
   return (FluxRusanov(Ul.W(), Ur.W()));
}

/*********************************************************
 * Routine: FluxHLLE (Harten-Lax-van Leer flux function) *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux given left and right solution states by using    *
 * the so-called Harten-Lax-van Leer approximation for   *
 * the fluxes.  See Harten, Lax, van Leer (1983).        *
 *                                                       *
 *********************************************************/
Euler1D_cState FluxHLLE(const Euler1D_pState &Wl,
	      	        const Euler1D_pState &Wr) {

    double wavespeed_l, wavespeed_r;
    Euler1D_pState Wa, lambdas_l, lambdas_r, lambdas_a;
    Euler1D_cState Flux, dUrl;
    
    /* Evaluate the Roe-average primitive solution state. */
    
    Wa = RoeAverage(Wl, Wr);

   /* Evaluate the jumps in the conserved solution states. */

    dUrl = Wr.U()-Wl.U();

    /* Evaluate the left, right, and average state eigenvalues. */

    lambdas_l = Wl.lambda();
    lambdas_r = Wr.lambda();
    lambdas_a = Wa.lambda();

    /* Determine the intermediate state flux. */

    wavespeed_l = min(lambdas_l[1],
                      lambdas_a[1]);
    wavespeed_r = max(lambdas_r[NUM_VAR_EULER1D],
                      lambdas_a[NUM_VAR_EULER1D]);
 
    wavespeed_l = min(wavespeed_l, ZERO);
    wavespeed_r = max(wavespeed_r, ZERO);

    if (wavespeed_l >= ZERO) {
        Flux = Wl.F();
    } else if (wavespeed_r <= ZERO) {
        Flux = Wr.F();
    } else {
        Flux =   ( (wavespeed_r*Wl.F()-wavespeed_l*Wr.F())
                  +(wavespeed_l*wavespeed_r)*dUrl)/
                 (wavespeed_r-wavespeed_l);
    } /* endif */

    /* Return solution flux. */

    return (Flux);

}

Euler1D_cState FluxHLLE(const Euler1D_cState &Ul,
	      	        const Euler1D_cState &Ur) {
   return (FluxHLLE(Ul.W(), Ur.W()));
}

/*********************************************************
 * Routine: FluxLinde (Timur Linde's flux function)      *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux given left and right solution states by using    *
 * the Linde approximation for the fluxes.  See Linde    *
 * (1998).                                               *
 *                                                       *
 *********************************************************/
Euler1D_cState FluxLinde(const Euler1D_pState &Wl,
                         const Euler1D_pState &Wr) {

    double wavespeed_l, wavespeed_r, wavespeed_m, da, ca, dU, alpha;
    Euler1D_pState Wa, lambdas_l, lambdas_r, lambdas_a;
    Euler1D_cState Flux, dFrl, dUrl, dFwave;
    
    /* Determine the left, intermediate, and right
       wave speeds. */

    Wa = RoeAverage(Wl, Wr);

    lambdas_l = Wl.lambda();
    lambdas_r = Wr.lambda();
    lambdas_a = Wa.lambda();

    wavespeed_l = min(lambdas_l[1],
                      lambdas_a[1]);
    wavespeed_r = max(lambdas_r[NUM_VAR_EULER1D],
                      lambdas_a[NUM_VAR_EULER1D]);
    wavespeed_m = Wa.v;

    /* Determine the intermediate state flux. */

    if (wavespeed_l >= ZERO) {
        Flux = Wl.F();
    } else if (wavespeed_r <= ZERO) {
        Flux = Wr.F();
    } else {
        dUrl = Wr.U()-Wl.U();
        dFrl = Wr.F()-Wl.F();
        da = Wa.d;
        ca = Wa.a();
        dU = fabs(dUrl.d)/da+
             fabs(dUrl.dv)/(da*ca)+
             fabs(dUrl.E)/(da*ca*ca);
        if (dU <= TOLER) {
            alpha = ZERO;
        } else {
            dU = ONE/dU;
            dFwave = dFrl - wavespeed_m*dUrl;
            alpha = ONE - (fabs(dFwave.d)/(da*ca)+
                           fabs(dFwave.dv)/(da*ca*ca)+
                           fabs(dFwave.E)/(da*ca*ca*ca))*dU;
 	    alpha = max(ZERO, alpha);
        } /* endif */
 
        Flux =   ( (wavespeed_r*Wl.F()-wavespeed_l*Wr.F())
                  +(wavespeed_l*wavespeed_r)*
		   (ONE-(ONE-max(wavespeed_m/wavespeed_r,wavespeed_m/wavespeed_l))*alpha)*dUrl)/
                   (wavespeed_r-wavespeed_l);
    } /* endif */

    /* Return solution flux. */

    return (Flux);

}

Euler1D_cState FluxLinde(const Euler1D_cState &Ul,
	      	         const Euler1D_cState &Ur) {
   return (FluxLinde(Ul.W(), Ur.W()));
}

/*********************************************************
 * Routine: FluxHLLC (Tito Toro's flux function)         *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux given left and right solution states by using    *
 * the HLLC approximation for the fluxes.  See Toro,     *
 * Spruce, and Speares (1994).                           *
 *                                                       *
 *********************************************************/
Euler1D_cState FluxHLLC(const Euler1D_pState &Wl,
                        const Euler1D_pState &Wr) {

    double wavespeed_l, wavespeed_r, wavespeed_m;
    double al, ar, CL, CR, Z, ql, qr;
    double vm, pm, aml, amr;

    Euler1D_cState Flux, Uml, Umr;
    
    /* Determine the left, intermediate, and right
       wave speeds. */

    al = Wl.a();
    ar = Wr.a();
    CL=Wl.v+TWO*al/Wl.gm1;
    CR=Wr.v-TWO*ar/Wr.gm1;
    Z = (ar/al)*pow((Wl.p/Wr.p), HALF*Wl.gm1/Wl.g);
    vm = (CL*Z+CR)/(ONE+Z);
    aml = al-HALF*Wl.gm1*(vm-Wl.v);
    pm = Wl.p*pow((aml/al), TWO*Wl.g/Wl.gm1);
    amr = ar+HALF*Wr.gm1*(vm-Wr.v);
    pm = HALF*(pm + Wr.p*pow((amr/ar), TWO*Wr.g/Wr.gm1));

    if (pm/Wl.p <= ONE) {
      ql = ONE;
    } else {
      ql = sqrt(ONE+HALF*((Wl.g+ONE)/Wl.g)*(pm/Wl.p));
    } /* endif */
    wavespeed_l = Wl.v - ql*al;

    if (pm/Wr.p <= ONE) {
      qr = ONE;
    } else {
      qr = sqrt(ONE+HALF*((Wr.g+ONE)/Wr.g)*(pm/Wr.p));
    } /* endif */
    wavespeed_r = Wr.v + qr*ar;

    wavespeed_m = vm;

    /* Determine the intermediate state flux. */

    if (wavespeed_l >= ZERO) {
        Flux = Wl.F();
    } else if (wavespeed_r <= ZERO) {
        Flux = Wr.F();
    } else  if (wavespeed_m >= ZERO) {
        Uml = Euler1D_cState(Wl.d, Wl.d*wavespeed_m, 
          		     Wl.E()+Wl.d*(wavespeed_m-Wl.v)*(wavespeed_m+
        		     Wl.p/(Wl.d*(wavespeed_l-Wl.v))))*
              ((wavespeed_l-Wl.v)/(wavespeed_l-wavespeed_m));
        Flux = Wl.F()+wavespeed_l*(Uml-Wl.U());
    } else {
        Umr = Euler1D_cState(Wr.d, Wr.d*wavespeed_m, 
          		     Wr.E()+Wr.d*(wavespeed_m-Wr.v)*(wavespeed_m+
        		     Wr.p/(Wr.d*(wavespeed_r-Wr.v))))*
              ((wavespeed_r-Wr.v)/(wavespeed_r-wavespeed_m));
        Flux = Wr.F()+wavespeed_r*(Umr-Wr.U());
    } /* endif */

    /* Return solution flux. */

    return (Flux);

}

Euler1D_cState FluxHLLC(const Euler1D_cState &Ul,
	      	        const Euler1D_cState &Ur) {
   return (FluxHLLC(Ul.W(), Ur.W()));
}

/*********************************************************
 * Routine: FluxOsher (Osher's flux function)            *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux given left and right solution states by using    *
 * Osher's approximation for the fluxes.  The O-ordering *
 * is used.  See Osher and Chakravarthy (1981, 1983).    *
 *                                                       *
 *********************************************************/
Euler1D_cState FluxOsher(const Euler1D_pState &Wl,
                         const Euler1D_pState &Wr) {

    double al, ar, CL, CR, Z, ql, qr;
    double dml, dmr, vm, pm, aml, amr;
    double ds, vs, as, ps;
    double wavespeed_l, wavespeed_r;

    Euler1D_cState Flux;

    /* Determine the left and right intermediate solution states 
       assuming O-ordering. */

    al = Wl.a();
    ar = Wr.a();
    CL=Wl.v-TWO*al/Wl.gm1;
    CR=Wr.v+TWO*ar/Wr.gm1;
    Z = (ar/al)*pow((Wl.p/Wr.p), HALF*Wl.gm1/Wl.g);
    vm = (CL*Z+CR)/(ONE+Z);
    aml = al+HALF*Wl.gm1*(vm-Wl.v);
    pm = Wl.p*pow((aml/al), TWO*Wl.g/Wl.gm1);
    amr = ar-HALF*Wr.gm1*(vm-Wr.v);
    pm = HALF*(pm + Wr.p*pow((amr/ar), TWO*Wr.g/Wr.gm1));
    dml = Wl.g*pm/sqr(aml);
    dmr = Wr.g*pm/sqr(amr);

    /* Evaluate the left and right wave speeds. */

    wavespeed_l = Wl.v + al;
    wavespeed_r = Wr.v - ar;

    /* Determine the appropriate flux, depending on the 
       wave speeds and the path integration. */

    // Start at left state
    Flux = F(Wl);

    // Add left wave
    if (wavespeed_l >= ZERO && vm + aml >= ZERO) {
      //Flux = Flux;    
    } else if (wavespeed_l < ZERO && vm + aml < ZERO) {
      Flux += F(Euler1D_pState(dml, vm, pm))-F(Wl);
    } else  if (wavespeed_l >= ZERO && vm + aml < ZERO) {
      vs =  (Wl.gm1*Wl.v-TWO*al)/(Wl.g+ONE);
      as = - vs;
      ps = Wl.p*pow((as/al), TWO*Wl.g/Wl.gm1);
      ds = Wl.g*ps/sqr(as);
      Flux += F(Euler1D_pState(dml, vm, pm))-
              F(Euler1D_pState(ds, vs, ps));    
    } else {
      vs =  (Wl.gm1*Wl.v-TWO*al)/(Wl.g+ONE);
      as = - vs;
      ps = Wl.p*pow((as/al), TWO*Wl.g/Wl.gm1);
      ds = Wl.g*ps/sqr(as);
      Flux += F(Euler1D_pState(ds, vs, ps))-F(Wl);    
    } /* endif */   

    // Add middle wave
    if (vm >= ZERO) {
      //Flux = Flux;    
    } else {
      Flux += F(Euler1D_pState(dmr, vm, pm))-
              F(Euler1D_pState(dml, vm, pm));    
    } /* endif */

    // Add right wave
    if (wavespeed_r >= ZERO && vm - amr >= ZERO) {
      //Flux = Flux;    
    } else if (wavespeed_r < ZERO && vm - amr < ZERO) {
      Flux += F(Wr)-F(Euler1D_pState(dmr, vm, pm));
    } else  if (wavespeed_r >= ZERO && vm - amr < ZERO) {
      vs =  (Wr.gm1*Wr.v+TWO*ar)/(Wr.g+ONE);
      as = vs;
      ps = Wr.p*pow((as/ar), TWO*Wr.g/Wr.gm1);
      ds = Wr.g*ps/sqr(as);
      Flux += F(Euler1D_pState(ds, vs, ps))-
              F(Euler1D_pState(dmr, vm, pm));    
    } else {
      vs =  (Wr.gm1*Wr.v+TWO*ar)/(Wr.g+ONE);
      as = vs;
      ps = Wr.p*pow((as/ar), TWO*Wr.g/Wr.gm1);
      ds = Wr.g*ps/sqr(as);
      Flux += F(Wr)-F(Euler1D_pState(ds, vs, ps));
    } /* endif */  

    /* Return solution flux. */

    return (Flux);

}

Euler1D_cState FluxOsher(const Euler1D_cState &Ul,
	      	        const Euler1D_cState &Ur) {
   return (FluxOsher(Ul.W(), Ur.W()));
}

/********************************************************
 * Routine: BC_Reflection (Reflection Boundary          *
 *                         Condition )                  *
 *                                                      *
 * This subroutine completely prescribes the boundary   *
 * conditions for a one-dimensional inviscid flow of a  *
 * polytropic (colorically and thermally perfect) gas at*
 * a moving solid boundary located at the left or right *
 * end of a one-dimensional numerical grid. Given the   *
 * flow state at the interior node nearest the solid    *
 * or closed boundary, Wi, and the boundary velocity,   *
 * Vbnd, this routine will completely specify the flow  *
 * state variables at the face of the solid boundary.   *
 * Refer to Gottlieb and Groth (1999).                  *
 *                                                      *
 * The subroutine uses the traditional nonstationary    *
 * shock and rarefaction wave patterns to connect the   *
 * interior and boundary states.  When a reflected shock*
 * wave pattern is present the well-known               *
 * Rankine-Hugoniot equations are used to determine the *
 * flow properties at the end node, and when a reflected*
 * rarefaction wave pattern is present the Riemann      *
 * invariant relations valid along the isentropic flow  *
 * characteristics are used to determine the end flow   *
 * conditions.                                          *
 *                                                      *
 ********************************************************/
Euler1D_pState BC_Reflection(const Euler1D_pState &Wi,
                             const double V_bnd,
                             const int End_Type) {

    int iend;
    double de, ve, pe, ae, dvei;

    /* Assign the velocity for the boundary state. */

    ve = V_bnd;

    /* Set the velocity difference parameter. */

    dvei = (ve-Wi.v)/Wi.a();
    switch(End_Type) {
      case LEFT_END_BOUNDARY :
        iend = -1;
        break;        
      case RIGHT_END_BOUNDARY :
        iend = 1;
        break;
      default:
        iend = 1;
        break;
    } /* endswitch */

    /* Reflected shock wave is present at boundary, 
       assign boundary pressure and sound speed accordingly. */

    if (double(iend)*dvei <= ZERO) {
       pe = Wi.p*(ONE+QUARTER*Wi.g*(Wi.g+ONE)*dvei*dvei-
            double(iend)*Wi.g*dvei*
	    sqrt(ONE+0.0625*(Wi.g+ONE)*(Wi.g+ONE)*dvei*dvei));
       ae = Wi.a()*sqrt(((Wi.g+ONE)+Wi.gm1*pe/Wi.p)/
                        ((Wi.g+ONE)+Wi.gm1*Wi.p/pe));
       de = Wi.g*pe/(ae*ae);

    /* Reflected rarefaction wave is present at boundary, 
       assign end pressure and sound speed accordingly. */

    } else if (double(iend)*dvei > ZERO &&
               double(iend)*dvei < TWO*Wi.gm1i) {
       ae = Wi.a()*(ONE-HALF*double(iend)*Wi.gm1*dvei);
       pe = Wi.p*pow( ae/Wi.a(), TWO*Wi.g*Wi.gm1i);
       de = Wi.g*pe/(ae*ae);

    /* Vacuum condition is present. */
    
    } else {
       pe = ZERO;
       ae = ZERO;
       de = ZERO;
    } /* endif */

    /* Return boundary solution state. */

    return (Euler1D_pState(de, ve, pe));

}

/********************************************************
 * Routine: BC_Characteristic (Characteristic-Based     *
 *                             Boundary Condition)      *
 *                                                      *
 * This subroutine completely prescribes the boundary   *
 * conditions for a one-dimensional inviscid flow of a  *
 * polytropic (calorically and thermally perfect) gas at*
 * a boundary based on a characteristic approach.  Given*
 * the flow state at the interior node nearest the      *
 * boundary of the numerical grid, Wi, and the desired  *
 * flow state to be imposed at the boundary, Wo, this   *
 * routine will determine the flow state at the boundary*
 * based on the solution of a simplified Riemann problem*
 * as described by Gottlieb and Groth (1999).  The      *
 * imposition of the boundary-data respects the         *
 * directions of propogation for the solution           *
 * characteristics at the boundary and thereby ensures  *
 * that the boundary conditions are not ill-posed (i.e.,*
 * the boundary data is not under- or over-determined). *
 *                                                      *
 ********************************************************/
Euler1D_pState BC_Characteristic(const Euler1D_pState &Wi,
	      	                 const Euler1D_pState &Wo,
                                 const int End_Type) {

    int iend;
    char pattern;
    double mi, poi, pai, pbi, mab, mac1, mac2, mbc1,
           mc1c2, mbd, mad;
    double de, ve, pe, ae;

    /* Determine the Mach number at the interior node and the
       imposed boundary to interior node pressure ratio. */

    mi = Wi.M();
    poi = Wo.p/Wi.p;
    switch(End_Type) {
      case LEFT_END_BOUNDARY :
        iend = -1;
        break;        
      case RIGHT_END_BOUNDARY :
        iend = 1;
        break;
      default:
        iend = 1;
        break;
    } /* endswitch */

    /* Identify the wave pattern between the interior and 
       boundary node. */

    if (poi >= ONE) {
       pai = pow(HALF*(Wi.g+ONE), TWO*Wi.g*Wi.gm1i);
       pbi = pow(HALF*(Wi.g+ONE)+HALF*Wi.gm1*(Wo.a()/Wi.a()), 
                 TWO*Wi.g*Wi.gm1i);
       if (poi <= pai) {
          mab = ONE;
          // Pattern A, pi <= po <= pa.
          if (double(iend)*mi >= mab) {
	     pattern='a';
	     goto impose_boundary_conditions;
          } /* endif */
          mbc1 = TWO*Wi.gm1i*(pow(poi, HALF*Wi.gm1/Wi.g)-ONE);
          // Pattern B, pi <= po <= pa.
          if (double(iend)*mi >= mbc1) {
	     pattern='b';
	     goto impose_boundary_conditions;
          } /* endif */
          mc1c2 = mbc1 - Wo.a()/Wi.a();
          // Pattern C1, pi <= po <= pa.
          if (double(iend)*mi >= mc1c2) {
	     pattern='c';
	     goto impose_boundary_conditions;
          // Pattern C2, pi <= po <= pa.
          } else {
	     pattern='c';
	     goto impose_boundary_conditions;
          } /* endif */
       } else if (poi <= pbi) {
          mac1 = ONE;
          // Pattern A, pa < po <= pb.
          if (double(iend)*mi >= mac1) {
	     pattern='a';
	     goto impose_boundary_conditions;
          } /* endif */
          mbc1 = TWO*Wi.gm1i*(pow(poi, HALF*Wi.gm1/Wi.g)-ONE);
          mc1c2 = mbc1 - Wo.a()/Wi.a();
          // Pattern C1, pa < po <= pb.
          if (double(iend)*mi >= mc1c2) {
	     pattern='c';
	     goto impose_boundary_conditions;
          // Pattern C2, pa < po <= pb.
          } else {
	     pattern='c';
	     goto impose_boundary_conditions;
          } /* endif */
       } else {
          mac2 = ONE;
          // Pattern A, po > pb.
          if (double(iend)*mi >= mac2) {
	     pattern='a';
	     goto impose_boundary_conditions;
          // Pattern C2, po > pb.
          } else {
	     pattern='c';
	     goto impose_boundary_conditions;
          } /* endif */
       } /* endif */
    } else {
      mad = ONE;
      // Pattern A, po < pi.
      if (double(iend)*mi >= mad) {
         pattern='a';
         goto impose_boundary_conditions;
      } /* endif */
      mbd = (Wi.g+ONE)*Wi.gm1i*pow(poi, HALF*Wi.gm1/Wi.g)-
            TWO*Wi.gm1i;
      // Pattern D, po < pi.
      if (double(iend)*mi >= mbd) {
         pattern='d';
         goto impose_boundary_conditions;
      } /* endif */
      mbc1 = TWO*Wi.gm1i*(pow(poi, HALF*Wi.gm1/Wi.g)-ONE);
      // Pattern B, po < pi.
      if (double(iend)*mi >= mbc1) {
         pattern='b';
         goto impose_boundary_conditions;
      } /* endif */
      mc1c2 = mbc1 - Wo.a()/Wi.a();
      // Pattern C1, po < pi.
      if (double(iend)*mi >= mc1c2) {
         pattern='c';
         goto impose_boundary_conditions;
      // Pattern C2, po < pi.
      } else {
         pattern='c';
         goto impose_boundary_conditions;
      } /* endif */
    } /* endif */

    /* Prescribe the appropriate boundary conditions. */

    impose_boundary_conditions: ;
    switch(pattern) {
      case 'a' :
        de = Wi.d;
        ve = Wi.v;
        pe = Wi.p;
        break;        
      case 'b' :
        pe = Wo.p;
        ve = Wi.v -
             TWO*double(iend)*Wi.a()*Wi.gm1i*
             (pow(poi, HALF*Wi.gm1/Wi.g)-ONE);
        de = Wi.d*pow(poi, ONE/Wi.g);
        break;
      case 'c' :
        pe = Wo.p;
        ve = Wi.v -
             TWO*double(iend)*Wi.a()*Wi.gm1i*
             (pow(poi, HALF*Wi.gm1/Wi.g)-ONE);
        de = Wo.d;
        break;
      case 'd' :
        pe = Wi.p*pow((TWO/(Wi.g+ONE)+
                      (Wi.gm1/(Wi.g+ONE))*double(iend)*mi), 
                      TWO*Wi.g*Wi.gm1i);
        ve = Wi.v -
             TWO*double(iend)*Wi.a()*Wi.gm1i*
             (pow(pe/Wi.p, HALF*Wi.gm1/Wi.g)-ONE);
        de = Wi.d*pow(pe/Wi.p, ONE/Wi.g);
        break;
      default:
        de = Wo.d;
        ve = Wo.v;
        pe = Wo.p;
        break;
    } /* endswitch */

    /* Return boundary solution state. */

    return (Euler1D_pState(de, ve, pe));

}

/********************************************************
 * Routine: BC_Characteristic_Pressure                  *
 *   (Characteristic-Based Boundary Condition with      *
 *    Static Pressure Specified Whenever Possible)      *
 *                                                      *
 * This function returns the boundary solution state    *
 * for a given boundary given the primitive solution    *
 * state on the interior of the boundary, Wi, the       *
 * desired flow state to be imposed at the boundary,    *
 * Wo.  A simplified characteristic analysis is         *
 * used to specify the boundary flow state in which the *
 * static pressure is specified whenever possible.  The *
 * imposition of the boundary-data respects the         *
 * directions of propogation for the solution           *
 * characteristics at the boundary and thereby ensures  *
 * that the boundary conditions are not ill-posed (i.e.,*
 * the boundary data is not under- or over-determined). *
 * The following procedure is adopted:                  *
 *                                                      *
 * 1) for supersonic outflow: constant extrapolation    *
 *    is employed to specify boundary conditions        *
 *    using known interior solution state,              *
 * 2) for subsonic outflow: boundary conditions         *
 *    specified by employing a 1D unsteady isentropic   *
 *    wave approximation to match boundary pressure,    *
 * 3) for subsonic inflow: boundary conditions          *
 *    specified by employing a 1D unsteady isentropic   *
 *    wave approximation to match both boundary state   *
 *    pressure and sound speed,                         *
 * 4) for supersonic inflow: the known boundary state   *
 *    is used to specify the boundary state.            *
 *                                                      *
 ********************************************************/
Euler1D_pState BC_Characteristic_Pressure(const Euler1D_pState &Wi,
                                          const Euler1D_pState &Wo,
                                          const int End_Type) {

    int iend;
    double mi, db, vb, pb, ab;

    /* Determine the Mach number at the interior node. */

    mi = Wi.M();
    switch(End_Type) {
      case LEFT_END_BOUNDARY :
        iend = -1;
        break;        
      case RIGHT_END_BOUNDARY :
        iend = 1;
        break;
      default:
        iend = 1;
        break;
    } /* endswitch */

    /* Boundary condition for supersonic outflow. */

    if (double(iend)*mi >= ONE) {
       db = Wi.d;
       vb = Wi.v;
       pb = Wi.p;
   
    /* Boundary condition for subsonic outflow. 
       Pressure specified. */

    } else if (double(iend)*mi >= ZERO) {
       pb = Wo.p;
       db = Wi.d*pow(pb/Wi.p, ONE/Wi.g);
       ab = sqrt(Wi.g*pb/db);
       vb = Wi.v + TWO*(Wi.a()-ab)*Wi.gm1i;

    /* Boundary condition for subsonic inflow. 
       Pressure specified. */

    } else if (double(iend)*mi >= -ONE) {
       pb = Wo.p;
       db = Wo.d;
       ab = sqrt(Wo.g*pb/db);
       vb = Wi.v + TWO*(Wi.a()-ab)*Wo.gm1i;

    /* Boundary condition for supersonic inflow.  */

    } else {
       db = Wo.d;
       vb = Wo.v;
       pb = Wo.p;
    } /* endif */

    /* Return boundary solution state. */

    return (Euler1D_pState(db, vb, pb));
       
}

/********************************************************
 * Routine: BC_Characteristic_Mach_Number               *
 *   (Characteristic-Based Boundary Condition with      *
 *    Flow Mach Number Specified Whenever Possible)     *
 *                                                      *
 * This function returns the boundary solution state    *
 * for a given boundary given the primitive solution    *
 * state on the interior of the boundary, Wi, the       *
 * desired flow state to be imposed at the boundary,    *
 * Wo.  A simplified characteristic analysis is         *
 * used to specify the boundary flow state in which the *
 * flow Mach number is specified whenever possible.     *
 * The imposition of the boundary-data respects the     *
 * directions of propogation for the solution           *
 * characteristics at the boundary and thereby ensures  *
 * that the boundary conditions are not ill-posed (i.e.,*
 * the boundary data is not under- or over-determined). *
 * The following procedure is adopted:                  *
 *                                                      *
 * 1) for supersonic outflow: constant extrapolation    *
 *    is employed to specify boundary conditions        *
 *    using known interior solution state,              *
 * 2) for subsonic outflow: boundary conditions         *
 *    specified by employing a 1D unsteady isentropic   *
 *    wave approximation to match boundary Mach number, *
 * 3) for subsonic inflow: boundary conditions          *
 *    specified by employing a 1D unsteady isentropic   *
 *    wave approximation to match both boundary state   *
 *    Mach number and density,                          *
 * 4) for supersonic inflow: the known boundary state   *
 *    is used to specify the boundary state.            *
 *                                                      *
 ********************************************************/
Euler1D_pState BC_Characteristic_Mach_Number(const Euler1D_pState &Wi,
                                             const Euler1D_pState &Wo,
                                             const int End_Type) {

    int iend;
    double mi, mb, db, vb, pb, ab;

    /* Determine the Mach number at the interior node. */

    mi = Wi.M();
    switch(End_Type) {
      case LEFT_END_BOUNDARY :
        iend = -1;
        break;        
      case RIGHT_END_BOUNDARY :
        iend = 1;
        break;
      default:
        iend = 1;
        break;
    } /* endswitch */

    /* Boundary condition for supersonic outflow. */

    if (double(iend)*mi >= ONE) {
       db = Wi.d;
       vb = Wi.v;
       pb = Wi.p;
   
    /* Boundary condition for subsonic outflow. 
       Mach number specified. */

    } else if (double(iend)*mi >= ZERO) {
       db = Wi.d;
       mb = Wo.v/Wo.a();
       ab= (Wi.v+TWO*Wi.a()*Wi.gm1i)/
           (mb+TWO*Wi.gm1i);
       pb = db*ab*ab/Wi.g;
       vb = mb*ab;

    /* Boundary condition for subsonic inflow. 
       Pressure specified. */

    } else if (double(iend)*mi >= -ONE) {
       db = Wo.d;
       mb = Wo.v/Wo.a();
       ab= (Wi.v+TWO*Wi.a()*Wi.gm1i)/
           (mb+TWO*Wi.gm1i);
       pb = db*ab*ab/Wo.g;
       vb = mb*ab;

    /* Boundary condition for supersonic inflow.  */

    } else {
       db = Wo.d;
       vb = Wo.v;
       pb = Wo.p;
    } /* endif */

    /* Return boundary solution state. */

    return (Euler1D_pState(db, vb, pb));
       
}

/********************************************************
 * Routine: BC_Open_End (Open End Boundary Condition)   *
 *                                                      *
 * This subroutine completely prescribes the boundary   *
 * conditions for a one-dimensional inviscid flow of a  *
 * polytropic (calorically and thermally perfect) gas at*
 * an open end of a duct to a large reservoir or        *
 * surrounding atmosphere.  The open end can be located *
 * at either the left or right end of a one-dimensional *
 * numerical grid.  Given the flow state at the interior*
 * node nearest the open end, Wi, and the reservoir     *
 * stagnation conditions, Wo, this routine will         *
 * completely specify the flow state variables at  *
 * the end of the numerical grid, as described by       *
 * Gottlieb and Groth (1999).                           *
 *                                                      *
 * Using the wave pattern boundary equations of Gottlieb*
 * and Groth (1999), the subroutine will identify the   *
 * wave pattern that is present between the nearest     *
 * interior node and the end node.  This pattern will be*
 * one of ten possible patterns which are listed below  *
 * as follows:                                          *
 *                                                      *
 * Pattern A: supersonic jet outflow with oblique shock *
 *            waves outside the end of the duct         *
 *                                                      *
 * Pattern B: subsonic jet outflow with upstream moving *
 *            shock wave                                *
 *                                                      *
 * Pattern C: subsonic inflow with upstream moving shock*
 *            wave and contact surface                  *
 *                                                      *
 * Pattern D: sonic or choked inflow with upstream      *
 *            moving shock and contact surface and      *
 *            downstream facing rarefaction wave        *
 *                                                      *
 * Pattern E: supersonic jet outflow with oblique       *
 *            rarefaction waves outside the end of the  *
 *            duct                                      *
 *                                                      *
 * Pattern F: sonic or choked outflow with upstream     *
 *            moving rarefaction wave and supersonic    *
 *            flow outside the end of the duct          *
 *                                                      *
 * Pattern G: subsonic jet outflow with upstream moving *
 *            rarefaction wave                          *
 *                                                      *
 * Pattern H: subsonic inflow with upstream moving      *
 *            rarefaction wave and contact surface      *
 *                                                      *
 * Pattern I: sonic or choked inflow with upstream      *
 *            moving rarefaction wave and contact       *
 *            surface and downstream facing rarefaction *
 *            wave                                      * 
 *                                                      *
 * Pattern V: similar to pattern I except a vacuum state*
 *            now exists between the upstream moving    *
 *            contact surface and the downstream facing *
 *            rarefaction wave                          *
 *                                                      *
 * Once the wave pattern has been correctly identified  *
 * the subroutine uses the well-known Rankine-Hugoniot  *
 * equations across shocks, the Riemann invariant       *
 * relationships across rarefaction waves, and the      *
 * steady isentropic flow relations to determine the    *
 * primitive variable flow state at the end of the      *
 * one-dimensional grid.                                *
 *                                                      *
 ********************************************************/
Euler1D_pState BC_Open_End(const Euler1D_pState &Wi,
	      	           const Euler1D_pState &Wo,
                           const int End_Type) {

    int iend, i;
    char pattern;
    double mi, poi, pcrit, a, b, c, z, alpha, beta, eta, 
           mab, mbc, mcd, mdi, miv, mch, mhi, mfg, mgh, ms;
    double gamma, de, ve, pe, ae, pe1, pe2;

    /* Determine the Mach number at the interior node and the
       reservoir (atmosphere) to interior node pressure ratio. */

    mi = Wi.M();
    poi = Wo.p/Wi.p;
    gamma = Wi.g;
    switch(End_Type) {
      case LEFT_END_BOUNDARY :
        iend = -1;
        break;        
      case RIGHT_END_BOUNDARY :
        iend = 1;
        break;
      default:
        iend = 1;
        break;
    } /* endswitch */

    /* Identify the wave pattern between the interior node
       and the open end. */

    if (poi >= ONE) {
       mab = sqrt(ONE+HALF*(gamma+ONE)/gamma*(poi-ONE));
       if (double(iend)*mi >= mab) {
	  pattern = 'a';
          goto impose_boundary_conditions; 
       } /* endif */
       mbc = (poi-ONE)/(gamma*mab);
       if (double(iend)*mi >= mbc) {
          pattern = 'b';
          goto impose_boundary_conditions; 
       } /* endif */
       z = sqrt(HALF*(gamma+ONE));
       pcrit = pow(z, TWO*gamma/(gamma-ONE));
       if (poi >= pcrit) {
          mcd = (poi/pcrit-ONE)/(gamma*sqrt(ONE+HALF*(gamma+ONE)/gamma*
	        (poi/pcrit-ONE)))-Wo.a()/(z*Wi.a());
          if (double(iend)*mi >= mcd) {
	     pattern = 'c';
             goto impose_boundary_conditions; 
          } /* endif */
          mdi = TWO*Wo.a()*(pow(poi, HALF*(ONE-gamma)/gamma) - z)/
                (Wi.a()*(gamma-ONE));
          if (double(iend)*mi >= mdi) {
	     pattern = 'd';
             goto impose_boundary_conditions; 
          } /* endif */
          miv = -TWO/(gamma-ONE)-(gamma+ONE)*Wo.a()/(z*Wi.a()*(gamma-ONE));
          if (double(iend)*mi >  miv) {
	     pattern = 'i';
             goto impose_boundary_conditions; 
          } else {
	     pattern = 'v';
             goto impose_boundary_conditions; 
          } /* endif */
       } else {
	  mch = -Wo.a()*sqrt(TWO*(ONE-pow(poi, (ONE-gamma)/gamma) )/
                (gamma-ONE))/Wi.a();
          if (double(iend)*mi >= mch) {
	     pattern = 'c';
             goto impose_boundary_conditions; 
          } /* endif */
          mhi = -TWO*(ONE+HALF*(gamma-ONE)*Wo.a()/(z*Wi.a())-
		pow( poi/pcrit, HALF*(gamma-ONE)/gamma) )/
                (gamma-ONE);
          if (double(iend)*mi >= mhi) {
	     pattern = 'h';
             goto impose_boundary_conditions; 
          } /* endif */
          miv = -TWO/(gamma-ONE)-(gamma+ONE)*Wo.a()/(z*Wi.a()*(gamma-ONE));
          if (double(iend)*mi > miv) {
	     pattern = 'i';
             goto impose_boundary_conditions; 
          } else {
	     pattern = 'v';
             goto impose_boundary_conditions; 
          } /* endif */
       } /* endif */
    } else {
       if (double(iend)*mi >= ONE) {
	  pattern = 'e';
          goto impose_boundary_conditions; 
       } /* endif */
       z = pow(poi, HALF*(gamma-ONE)/gamma);
       mfg = (z*(gamma+ONE)-TWO)/(gamma-ONE);
       if (double(iend)*mi >= mfg) {
	  pattern = 'f';
          goto impose_boundary_conditions; 
       } /* endif */
       mgh = TWO*(z-ONE)/(gamma-ONE);
       if (double(iend)*mi >= mgh) {
	  pattern = 'g';
          goto impose_boundary_conditions; 
       } /* endif */
       z = sqrt(HALF*(gamma+ONE));
       pcrit = pow(z, TWO*gamma/(gamma-ONE));
       mhi = -TWO*(ONE+HALF*(gamma-ONE)*Wo.a()/(z*Wi.a())-
             pow( poi/pcrit, HALF*(gamma-ONE)/gamma) )/
             (gamma-ONE);
       if (double(iend)*mi >= mhi) {
	  pattern = 'h';
          goto impose_boundary_conditions; 
       } /* endif */
       miv = -TWO/(gamma-ONE)-(gamma+ONE)*Wo.a()/(z*Wi.a()*(gamma-ONE));
       if (double(iend)*mi > miv) {
	  pattern = 'i';
          goto impose_boundary_conditions; 
       } else {
	  pattern = 'v';
          goto impose_boundary_conditions; 
       } /* endif */
    } /* endif */

    /* Prescribe the appropriate boundary conditions. */

    impose_boundary_conditions: ;
    switch(pattern) {
      case 'a' :
      case 'e' :
        de = Wi.d;
        ve = Wi.v;
        pe = Wi.p;
        break;        
      case 'b' :
        pe = Wo.p;
        ve = Wi.v-double(iend)*Wi.a()*(poi-ONE)/
 	     (gamma*sqrt(ONE+HALF*(gamma+ONE)*(poi-ONE)/gamma));
	ae = Wi.a()*sqrt(((gamma+ONE)+(gamma-ONE)*poi)/
             ((gamma+ONE)+(gamma-ONE)/poi));
        de = gamma*pe/(ae*ae);
        break;
      case 'c' :
	beta = (gamma-ONE)/gamma;
	alpha = TWO*gamma/(gamma+ONE);
	eta = TWO/(gamma+ONE);
	a = TWO*sqr(Wo.a())/(sqr(Wi.a())*(gamma-ONE));
	b = pow(ONE/poi, beta);
	c = EIGHT*gamma*sqr(Wo.a())*b/(gamma*sqr(Wi.a())*(gamma+ONE));
	if (poi >= pcrit) {
	  pe1 = poi+(double(iend)*mi-mbc)*(poi/pcrit-poi)/(mcd-mbc);
        } else {
	  pe1 = poi+(double(iend)*mi-mbc)*(ONE-poi)/(mch-mbc);
        } /* endif */
	ms = -double(iend)*sqrt(ONE+(pe1-ONE)/alpha);
	ve = mi+eta*(ms*ms-ONE)/ms;
        for (i = 0; i <= 50 ; ++i ) {
           ms = ms-(ve*ve-a*(ONE-b*pow(pe1, beta) ))/
                (TWO*eta*ve*(ms*ms+ONE)/(ms*ms)+c*ms*pow(pe1, -ONE/gamma) );
           ve = mi+eta*(ms*ms-ONE)/ms;
	   pe2 = ONE+alpha*(ms*ms-ONE);
           if (fabs(ONE-pe2/pe1) <= TOLER) break;
	   pe1 = pe2;
        } /* endfor */
	pe = pe2*Wi.p;
	ve = ve*Wi.a();
	ae = Wo.a()*pow(pe/Wo.p, HALF*beta);
        de = gamma*pe/(ae*ae);
        break;
      case 'd' :
      case 'i' :
      case 'v' :
        pe = Wo.p/pcrit;
        ae = Wo.a()/z;
	ve = -double(iend)*ae;
        de = gamma*pe/(ae*ae);
        break;
      case 'f' :
	ae = Wi.a()*(TWO+double(iend)*(gamma-ONE)*mi)/(gamma+ONE);
	ve = double(iend)*ae;
	pe = Wi.p*pow(ae/Wi.a(), TWO*gamma/(gamma-ONE));
        de = gamma*pe/(ae*ae);
        break;
      case 'g' :
	pe = Wo.p;
	ae = Wi.a()*pow(poi, HALF*(gamma-ONE)/gamma);
	ve = Wi.v+TWO*double(iend)*(Wi.a()-ae)/(gamma-ONE);
        de = gamma*pe/(ae*ae);
        break;
      case 'h' :
	beta = gamma-ONE;
	a = ONE+HALF*beta*sqr(Wo.a())*pow(ONE/poi, beta/gamma)/sqr(Wi.a());
	b = ONE+HALF*double(iend)*beta*mi;
	c = b*b-HALF*beta*sqr(Wo.a())/sqr(Wi.a());
	z = (b+sqrt(b*b-a*c))/a;
	pe = Wi.p*pow(z, TWO*gamma/beta);
	ae = Wo.a()*pow(pe/Wo.p , HALF*beta/gamma);
	ve = Wi.v+TWO*double(iend)*Wi.a()*(ONE-z)/beta;
        de = gamma*pe/(ae*ae);
        break;
      default:
        de = Wo.d;
        ve = Wo.v;
        pe = Wo.p;
        break;
    } /* endswitch */

    /* Return boundary solution state. */

    return (Euler1D_pState(de, ve, pe));

}


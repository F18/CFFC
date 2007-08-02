/* MHD1DState.cc:  Subroutines for 1D MHD Solution State Classes. */

/* Include 1D MHD solution state header file. */

#ifndef _MHD1D_STATE_INCLUDED
#include "MHD1DState.h"
#endif // _MHD1D_STATE_INCLUDED

/*************************************************************
 * MHD1D_pState -- Create storage and assign gas constants.  *
 *************************************************************/
double MHD1D_pState::g = GAMMA_MONATOMIC;
double MHD1D_pState::gm1 = GAMMA_MONATOMIC-ONE;
double MHD1D_pState::gm1i = ONE/(GAMMA_MONATOMIC-ONE);

/*************************************************************
 * MHD1D_cState -- Create storage and assign gas constants.  *
 *************************************************************/
double MHD1D_cState::g = GAMMA_MONATOMIC;
double MHD1D_cState::gm1 = GAMMA_MONATOMIC-ONE;
double MHD1D_cState::gm1i = ONE/(GAMMA_MONATOMIC-ONE);

/********************************************************
 * Routine: RoeAverage (Roe Averages)                   *
 *                                                      *
 * This function returns the Roe-averaged primitive     *
 * solution state given left and right primitive        *
 * solution variables.  See Roe (1981).                 *
 *                                                      *
 ********************************************************/
MHD1D_pState RoeAverage(const MHD1D_pState &Wl,
	      	        const MHD1D_pState &Wr) {

    double hl, hr, sdl, sdr;
    double da, ha, pa, aa2;
    Vector3D va, B1a;
    
    /* Determine the left and right state specific enthalpies
       and square roots of the density. */

    hl = Wl.h();
    hr = Wr.h();
    sdl = sqrt(Wl.d);
    sdr = sqrt(Wr.d);

    /* Determine the appropriate Roe averages. */
    /* Note approximate averages used here. */

    da = sdl*sdr;
    va = (sdl*Wl.v+sdr*Wr.v)/(sdl+sdr);
    B1a = ((sdl/Wl.d)*Wl.B1+(sdr/Wr.d)*Wr.B1)/(sdl+sdr);
    B1a = da*B1a;
    ha = (sdl*hl+sdr*hr)/(sdl+sdr);
    aa2 = Wl.gm1*(ha - HALF*va.sqr() - B1a.sqr()/da);
    pa = da*aa2/Wl.g;

    /* Return the Roe-averged state. */

    return (MHD1D_pState(da, va, B1a, Wl.B0, pa));
       
}

/********************************************************
 * Routine: WaveSpeedPos                                *
 *                                                      *
 * This function returns the positive parts of the      *
 * elemental wave speeds or eigenvalues.                *
 *                                                      *
 ********************************************************/
MHD1D_pState WaveSpeedPos(const MHD1D_pState &lambdas_a,
                          const MHD1D_pState &lambdas_l,
                          const MHD1D_pState &lambdas_r) {
  return (MHD1D_pState(HALF*(lambdas_a[1]+fabs(lambdas_a[1])),
                       HALF*(lambdas_a[2]+fabs(lambdas_a[2])),
                       HALF*(lambdas_a[3]+fabs(lambdas_a[3])),
                       HALF*(lambdas_a[4]+fabs(lambdas_a[4])),
                       HALF*(lambdas_a[5]+fabs(lambdas_a[5])),
                       HALF*(lambdas_a[6]+fabs(lambdas_a[6])),
                       HALF*(lambdas_a[7]+fabs(lambdas_a[7])),
                       HALF*(lambdas_a[8]+fabs(lambdas_a[8]))));
}

/********************************************************
 * Routine: WaveSpeedNeg                                *
 *                                                      *
 * This function returns the negative parts of the      *
 * elemental wave speeds or eigenvalues.                *
 *                                                      *
 ********************************************************/
MHD1D_pState WaveSpeedNeg(const MHD1D_pState &lambdas_a,
                          const MHD1D_pState &lambdas_l,
                          const MHD1D_pState &lambdas_r) {
  return (MHD1D_pState(HALF*(lambdas_a[1]-fabs(lambdas_a[1])),
                       HALF*(lambdas_a[2]-fabs(lambdas_a[2])),
                       HALF*(lambdas_a[3]-fabs(lambdas_a[3])),
                       HALF*(lambdas_a[4]-fabs(lambdas_a[4])),
                       HALF*(lambdas_a[5]-fabs(lambdas_a[5])),
                       HALF*(lambdas_a[6]-fabs(lambdas_a[6])),
                       HALF*(lambdas_a[7]-fabs(lambdas_a[7])),
                       HALF*(lambdas_a[8]-fabs(lambdas_a[8]))));
}

/********************************************************
 * Routine: WaveSpeedAbs                                *
 *                                                      *
 * This function returns the absolute values of the     *
 * elemental wave speeds or eigenvalues.                *
 *                                                      *
 ********************************************************/
MHD1D_pState WaveSpeedAbs(const MHD1D_pState &lambdas_a,
                          const MHD1D_pState &lambdas_l,
                          const MHD1D_pState &lambdas_r) {
  return (MHD1D_pState(fabs(lambdas_a[1]),
                       fabs(lambdas_a[2]),
                       fabs(lambdas_a[3]),
                       fabs(lambdas_a[4]),
                       fabs(lambdas_a[5]),
                       fabs(lambdas_a[6]),
                       fabs(lambdas_a[7]),
                       fabs(lambdas_a[8])));
}

/********************************************************
 * Routine: HartenFixPos (Harten Entropy Fix)           *
 *                                                      *
 * This function returns the positive parts of the      *
 * corrected elemental wave speeds or eigenvalues       *
 * according to the entropy fix of Harten (1983).       *
 *                                                      *
 ********************************************************/
MHD1D_pState HartenFixPos(const MHD1D_pState &lambdas_a,
                          const MHD1D_pState &lambdas_l,
                          const MHD1D_pState &lambdas_r) {
  return (MHD1D_pState(HartenFixPos(lambdas_a[1],
                                    lambdas_l[1],
                                    lambdas_r[1]),
                       HartenFixPos(lambdas_a[2],
                                    lambdas_l[2],
                                    lambdas_r[2]),
                       HartenFixPos(lambdas_a[3],
                                    lambdas_l[3],
                                    lambdas_r[3]),
                       HALF*(lambdas_a[4]+fabs(lambdas_a[4])),
                       HALF*(lambdas_a[5]+fabs(lambdas_a[5])),
                       HartenFixPos(lambdas_a[6],
                                    lambdas_l[6],
                                    lambdas_r[6]),
                       HartenFixPos(lambdas_a[7],
                                    lambdas_l[7],
                                    lambdas_r[7]),
                       HartenFixPos(lambdas_a[8],
                                    lambdas_l[8],
                                    lambdas_r[8])));
}

/********************************************************
 * Routine: HartenFixNeg (Harten Entropy Fix)           *
 *                                                      *
 * This function returns the negative parts of the      *
 * corrected elemental wave speeds or eigenvalues       *
 * according to the entropy fix of Harten (1983).       *
 *                                                      *
 ********************************************************/
MHD1D_pState HartenFixNeg(const MHD1D_pState &lambdas_a,
                          const MHD1D_pState &lambdas_l,
                          const MHD1D_pState &lambdas_r) {
  return (MHD1D_pState(HartenFixNeg(lambdas_a[1],
                                    lambdas_l[1],
                                    lambdas_r[1]),
                       HartenFixNeg(lambdas_a[2],
                                    lambdas_l[2],
                                    lambdas_r[2]),
                       HartenFixNeg(lambdas_a[3],
                                    lambdas_l[3],
                                    lambdas_r[3]),
                       HALF*(lambdas_a[4]-fabs(lambdas_a[4])),
                       HALF*(lambdas_a[5]-fabs(lambdas_a[5])),
                       HartenFixNeg(lambdas_a[6],
                                    lambdas_l[6],
                                    lambdas_r[6]),
                       HartenFixNeg(lambdas_a[7],
                                    lambdas_l[7],
                                    lambdas_r[7]),
                       HartenFixNeg(lambdas_a[8],
                                    lambdas_l[8],
                                    lambdas_r[8])));
}

/********************************************************
 * Routine: HartenFixAbs (Harten Entropy Fix)           *
 *                                                      *
 * This function returns the absolute values of the     *
 * corrected elemental wave speeds or eigenvalues       *
 * according to the entropy fix of Harten (1983).       *
 *                                                      *
 ********************************************************/
MHD1D_pState HartenFixAbs(const MHD1D_pState &lambdas_a,
                          const MHD1D_pState &lambdas_l,
                          const MHD1D_pState &lambdas_r) {
  return (MHD1D_pState(HartenFixAbs(lambdas_a[1],
                                    lambdas_l[1],
                                    lambdas_r[1]),
                       HartenFixAbs(lambdas_a[2],
                                    lambdas_l[2],
                                    lambdas_r[2]),
                       HartenFixAbs(lambdas_a[3],
                                    lambdas_l[3],
                                    lambdas_r[3]),
                       fabs(lambdas_a[4]),
                       fabs(lambdas_a[5]),
                       HartenFixAbs(lambdas_a[6],
                                    lambdas_l[6],
                                    lambdas_r[6]),
                       HartenFixAbs(lambdas_a[7],
                                    lambdas_l[7],
                                    lambdas_r[7]),
                       HartenFixAbs(lambdas_a[8],
                                    lambdas_l[8],
                                    lambdas_r[8])));
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
MHD1D_cState FluxRoe(const MHD1D_pState &Wl,
	      	     const MHD1D_pState &Wr) {

    int i;  
    MHD1D_pState Wa, wavespeeds, lambdas_l, lambdas_r, lambdas_a;
    MHD1D_cState Flux, dUrl;
    
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

    for ( i = 1 ; i <= NUM_VAR_MHD1D ; ++i ) {
        Flux -= HALF*wavespeeds[i]*(Wa.lc(i)*dUrl)*Wa.rc(i);
    } /* endfor */

    /* Return solution flux. */

    return (Flux);

}

MHD1D_cState FluxRoe(const MHD1D_cState &Ul,
	      	     const MHD1D_cState &Ur) {
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
MHD1D_cState FluxRusanov(const MHD1D_pState &Wl,
	      	         const MHD1D_pState &Wr) {

    int i;
    double wavespeed_max;
    MHD1D_pState Wa, wavespeeds, lambdas_l, lambdas_r, lambdas_a;
    MHD1D_cState Flux, dUrl;
    
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
    for ( i = 2 ; i <= NUM_VAR_MHD1D ; ++i ) {
        wavespeed_max = max(wavespeed_max, wavespeeds[i]);
    } /* endfor */

    Flux -= HALF*wavespeed_max*dUrl;

    /* Return solution flux. */

    return (Flux);

}

MHD1D_cState FluxRusanov(const MHD1D_cState &Ul,
	      	         const MHD1D_cState &Ur) {
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
MHD1D_cState FluxHLLE(const MHD1D_pState &Wl,
	      	      const MHD1D_pState &Wr) {

    double wavespeed_l, wavespeed_r;
    MHD1D_pState Wa, lambdas_l, lambdas_r, lambdas_a;
    MHD1D_cState Flux, dUrl;
    
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
    wavespeed_r = max(lambdas_r[NUM_VAR_MHD1D],
                      lambdas_a[NUM_VAR_MHD1D]);
 
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

MHD1D_cState FluxHLLE(const MHD1D_cState &Ul,
	      	      const MHD1D_cState &Ur) {
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
MHD1D_cState FluxLinde(const MHD1D_pState &Wl,
	      	       const MHD1D_pState &Wr) {

    double wavespeed_l, wavespeed_r, wavespeed_m, ca, dU, alpha;
    MHD1D_pState Wa, lambdas_l, lambdas_r, lambdas_a;
    MHD1D_cState Flux, dFrl, dUrl, dFwave;
    
    /* Evaluate the Roe-average primitive solution state. */
    
    Wa = RoeAverage(Wl, Wr);

    /* Evaluate the left, right, and average state eigenvalues. */

    lambdas_l = Wl.lambda();
    lambdas_r = Wr.lambda();
    lambdas_a = Wa.lambda();

    /* Determine the intermediate state flux. */

    wavespeed_l = min(lambdas_l[1],
                      lambdas_a[1]);
    wavespeed_r = max(lambdas_r[NUM_VAR_MHD1D],
                      lambdas_a[NUM_VAR_MHD1D]);

    if (wavespeed_l >= ZERO) {
        Flux = Wl.F();
    } else if (wavespeed_r <= ZERO) {
        Flux = Wr.F();
    } else {
        dUrl = Wr.U()-Wl.U();
        dFrl = Wr.F()-Wl.F();
        wavespeed_m = Wa.v.x;
	ca = Wa.a();
        dU = fabs(dUrl.d) +fabs(dUrl.dv.x) +
	     fabs(dUrl.dv.y) + fabs(dUrl.dv.z) +
	     fabs(dUrl.B1.x) + fabs(dUrl.B1.y) +
	     fabs(dUrl.B1.z) + fabs(dUrl.E1);
        if (dU <= TOLER) {
            alpha = ZERO;
        } else {
            dU = ONE/(ca*dU);
            dFwave = dFrl - wavespeed_m*dUrl;
            alpha = ONE - (fabs(dFwave.d)+ fabs(dFwave.dv.x) +
	   	           fabs(dFwave.dv.y) + fabs(dFwave.dv.z) +
		           fabs(dFwave.B1.x) + fabs(dFwave.B1.y) +
         		   fabs(dFwave.B1.z) + fabs(dFwave.E1))*dU;
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

MHD1D_cState FluxLinde(const MHD1D_cState &Ul,
	      	       const MHD1D_cState &Ur) {
   return (FluxLinde(Ul.W(), Ur.W()));
}

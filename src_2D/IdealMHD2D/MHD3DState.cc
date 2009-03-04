/*!\file MHD3DState.cc
  @brief Subroutines for 3D Ideal MHD Solution State Classes. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "MHD3DState.h"		// Include 3D ideal MHD solution state header file.


// ======= Member functions/variables of MHD3D_pState class =======

/**************************
 * Assign gas constants.  *
 **************************/
double MHD3D_pState::g = GAMMA_MONATOMIC;
double MHD3D_pState::gm1 = GAMMA_MONATOMIC-ONE;
double MHD3D_pState::gm1i = ONE/(GAMMA_MONATOMIC-ONE);
double MHD3D_pState::Mr_min = ONE;

/*
 * Set useful state constants
 */
const MHD3D_pState MHD3D_pState::MHD3D_W_REF(ONE, Vector3D_ZERO, Vector3D_ZERO,
					     Vector3D_ZERO, ONE);
const MHD3D_pState MHD3D_pState::MHD3D_W_ZERO(ZERO, Vector3D_ZERO, Vector3D_ZERO,
					      Vector3D_ZERO, ZERO);

/*!
 * This function returns the Roe-averaged primitive     
 * solution state given left and right primitive        
 * solution variables.  See Roe (1981).                 
 */
MHD3D_pState RoeAverage(const MHD3D_pState &Wl,
			const MHD3D_pState &Wr){

  double hl, hr, sdl, sdr;
  double da, ha, pa, aa2;
  Vector3D va, B1a;
    
  /* Determine the left and right state specific enthalpies
     and square roots of the density. */

  hl = Wl.h();
  hr = Wr.h();
  sdl = sqrt(Wl.d());
  sdr = sqrt(Wr.d());

  /* Determine the appropriate Roe averages. */
  /* Note approximate averages used here. */

  da = sdl*sdr;
  va = (sdl*Wl.v()+sdr*Wr.v())/(sdl+sdr);
  B1a = ((sdl/Wl.d())*Wl.B1()+(sdr/Wr.d())*Wr.B1())/(sdl+sdr);
  B1a = da*B1a;
  ha = (sdl*hl+sdr*hr)/(sdl+sdr);
  aa2 = Wl.gm1*(ha - HALF*va.sqr() - B1a.sqr()/da);
  pa = da*aa2/Wl.g;

  /* Return the Roe-averged state. */
  return MHD3D_pState(da, va, B1a, Wl.B0(), pa);
}

/*!
 * Roe's flux function \n
 *
 * This function returns the intermediate state solution
 * flux given left and right solution states by using   
 * the "linearized" approximate Riemann solver of Roe   
 * for the two states.  See Roe (1981).                 
 */
MHD3D_cState FluxRoe(const MHD3D_pState &Wl,
	      	     const MHD3D_pState &Wr) {

  MHD3D_pState Wa, wavespeeds, lambdas_l, lambdas_r, lambdas_a;
  MHD3D_cState Flux, dUrl;
    
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

  for (int i = 1 ; i <= NUM_VAR_MHD3D ; ++i ) {
    Flux -= HALF*wavespeeds[i]*(Wa.lc(i)*dUrl)*Wa.rc(i);
  } /* endfor */

  /* Return solution flux. */
  return Flux;
}

/*!
 * Rusanov flux function \n
 *
 * This function returns the intermediate state solution
 * flux given left and right solution states by using   
 * the Rusanov approximation for the fluxes.  See       
 * Rusanov (1964).                                      
 */
MHD3D_cState FluxRusanov(const MHD3D_pState &Wl,
	      	         const MHD3D_pState &Wr) {

  double wavespeed_max;
  MHD3D_pState Wa, wavespeeds, lambdas_l, lambdas_r, lambdas_a;
  MHD3D_cState Flux, dUrl;
    
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
  for (int i = 2 ; i <= NUM_VAR_MHD3D ; ++i ) {
    wavespeed_max = max(wavespeed_max, wavespeeds[i]);
  } /* endfor */

  Flux -= HALF*wavespeed_max*dUrl;

  /* Return solution flux. */
  return Flux;

}

/*!
 * FluxHLLE (Harten-Lax-van Leer flux function) \n
 *
 * This function returns the intermediate state solution 
 * flux given left and right solution states by using    
 * the so-called Harten-Lax-van Leer approximation for   
 * the fluxes.  See Harten, Lax, van Leer (1983).        
 */
MHD3D_cState FluxHLLE(const MHD3D_pState &Wl,
	      	      const MHD3D_pState &Wr) {

  double wavespeed_l, wavespeed_r;
  MHD3D_pState Wa, lambdas_l, lambdas_r, lambdas_a;
  MHD3D_cState Flux, dUrl;
    
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
  wavespeed_r = max(lambdas_r[NUM_VAR_MHD3D],
		    lambdas_a[NUM_VAR_MHD3D]);
 
  wavespeed_l = min(wavespeed_l, ZERO);
  wavespeed_r = max(wavespeed_r, ZERO);

  if (wavespeed_l >= ZERO) {
    Flux = Wl.F();
  } else if (wavespeed_r <= ZERO) {
    Flux = Wr.F();
  } else {
    Flux = ( ( (wavespeed_r*Wl.F()-wavespeed_l*Wr.F())
	       +(wavespeed_l*wavespeed_r)*dUrl)/
	     (wavespeed_r-wavespeed_l) );
  } /* endif */

  /* Return solution flux. */
  return Flux;
}
 
/*!
 * FluxLinde (Timur Linde's flux function)
 *                                                      
 * This function returns the intermediate state solution
 * flux given left and right solution states by using   
 * the Linde approximation for the fluxes.  See Linde   
 * (1998).                                              
 */
MHD3D_cState FluxLinde(const MHD3D_pState &Wl,
	      	       const MHD3D_pState &Wr) {

  double wavespeed_l, wavespeed_r, wavespeed_m, ca, dU, alpha;
  MHD3D_pState Wa, lambdas_l, lambdas_r, lambdas_a;
  MHD3D_cState Flux, dFrl, dUrl, dFwave;
    
  /* Evaluate the Roe-average primitive solution state. */
  Wa = RoeAverage(Wl, Wr);

  /* Evaluate the left, right, and average state eigenvalues. */
  lambdas_l = Wl.lambda();
  lambdas_r = Wr.lambda();
  lambdas_a = Wa.lambda();

  /* Determine the intermediate state flux. */
  wavespeed_l = min(lambdas_l[1],
		    lambdas_a[1]);
  wavespeed_r = max(lambdas_r[NUM_VAR_MHD3D],
		    lambdas_a[NUM_VAR_MHD3D]);

  if (wavespeed_l >= ZERO) {
    Flux = Wl.F();
  } else if (wavespeed_r <= ZERO) {
    Flux = Wr.F();
  } else {
    dUrl = Wr.U()-Wl.U();
    dFrl = Wr.F()-Wl.F();
    wavespeed_m = Wa.vx();
    ca = Wa.a();
    dU = fabs(dUrl.d()) +fabs(dUrl.dvx()) +
      fabs(dUrl.dvy()) + fabs(dUrl.dvz()) +
      fabs(dUrl.B1x()) + fabs(dUrl.B1y()) +
      fabs(dUrl.B1z()) + fabs(dUrl.E1());
    if (dU <= TOLER) {
      alpha = ZERO;
    } else {
      dU = ONE/(ca*dU);
      dFwave = dFrl - wavespeed_m*dUrl;
      alpha = ONE - (fabs(dFwave.d())+ fabs(dFwave.dvx()) +
		     fabs(dFwave.dvy()) + fabs(dFwave.dvz()) +
		     fabs(dFwave.B1x()) + fabs(dFwave.B1y()) +
		     fabs(dFwave.B1z()) + fabs(dFwave.E1()))*dU;
      alpha = max(ZERO, alpha);
    } /* endif */
 
    Flux = ( ( (wavespeed_r*Wl.F()-wavespeed_l*Wr.F())
	       +(wavespeed_l*wavespeed_r)*
	       (ONE-(ONE-max(wavespeed_m/wavespeed_r,wavespeed_m/wavespeed_l))*alpha)*dUrl)/
	     (wavespeed_r-wavespeed_l) );
  } /* endif */

  /* Return solution flux. */
  return Flux;

}

// ======= Member functions/variables of MHD3D_cState class =======

/**************************
 * Assign gas constants.  *
 **************************/
double MHD3D_cState::g = GAMMA_MONATOMIC;
double MHD3D_cState::gm1 = GAMMA_MONATOMIC-ONE;
double MHD3D_cState::gm1i = ONE/(GAMMA_MONATOMIC-ONE);
double MHD3D_cState::Mr_min = ONE;

/*
 * Set useful state constants
 */
const MHD3D_cState MHD3D_cState::MHD3D_U_REF(MHD3D_pState::MHD3D_W_REF);
const MHD3D_cState MHD3D_cState::MHD3D_U_ZERO(MHD3D_pState::MHD3D_W_ZERO);

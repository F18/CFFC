/* Euler3DState.cc:  Subroutines for 3D Euler Solution State Classes. */

/* Include 3D Euler solution state header file. */

#ifndef _EULER3D_STATE_INCLUDED
#include "Euler3DState.h"
#endif // _EULER3D_STATE_INCLUDED

/*************************************************************
 * Euler3D_pState -- Create storage and assign gas constants.*
 *************************************************************/
double Euler3D_pState::g = GAMMA_AIR;
double Euler3D_pState::gm1 = GAMMA_AIR-ONE;
double Euler3D_pState::gm1i = ONE/(GAMMA_AIR-ONE);
double Euler3D_pState::R = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);

/*************************************************************
 * Euler3D_cState -- Create storage and assign gas constants.*
 *************************************************************/
double Euler3D_cState::g = GAMMA_AIR;
double Euler3D_cState::gm1 = GAMMA_AIR-ONE;
double Euler3D_cState::gm1i = ONE/(GAMMA_AIR-ONE);
double Euler3D_cState::R = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);



/********************************************************
 * Routine: RoeAverage (Roe Averages)                   *
 *                                                      *
 * This function returns the Roe-averaged primitive     *
 * solution state given left and right primitive        *
 * solution variables.  See Roe (1981).                 *
 *                                                      *
 ********************************************************/
Euler3D_pState RoeAverage(const Euler3D_pState &Wl,
	      	          const Euler3D_pState &Wr) {

    double hl, hr, sdl, sdr;
    double da, ua, va,wa, pa, aa2, ha, ga, gam1;

    /* Determine the left and right state specific enthalpies
       and square roots of the density. */

    hl = Wl.h();
    hr = Wr.h();
    sdl = sqrt(Wl.d());
    sdr = sqrt(Wr.d());

    /* Determine the appropriate Roe averages. */
    ga = Wl.g;
    gam1 = Wl.gm1;
    da = sdl*sdr;
    ua = (sdl*Wl.vx()+sdr*Wr.vx())/(sdl+sdr);
    va = (sdl*Wl.vy()+sdr*Wr.vy())/(sdl+sdr);
    wa = (sdl*Wl.vz()+sdr*Wr.vz())/(sdl+sdr);

    ha = (sdl*hl+sdr*hr)/(sdl+sdr);
    aa2 = gam1*(ha-HALF*(sqr(ua)+sqr(va)+sqr(wa)));
    pa = da*aa2/ga;

    /* Return the Roe-averged state. */

    return (Euler3D_pState(da, ua, va, wa, pa));
       
}

/*********************************************************
 * Routine: FluxHLLE_x (Harten-Lax-van Leer flux         *
 *                      function, x-direction)           *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for the x-direction given left and right         *
 * solution states by using the so-called Harten-Lax-    *
 * van Leer approximation for the fluxes.  See Harten,   *
 * Lax, van Leer (1983).                                 *
 *                                                       *
 *********************************************************/
Euler3D_cState FluxHLLE_x(const Euler3D_pState &Wl,
	      	          const Euler3D_pState &Wr) {

    double wavespeed_l, wavespeed_r;
   
    /* solnvec in  Wa (lambdas_l, lambdas_r, lambdas_a) 
       is allocated using new  */ 
    Euler3D_pState Wa(0), lambdas_l(0), lambdas_r(0), lambdas_a(0);
    Euler3D_cState Flux, dUrl;
    
    /* Evaluate the Roe-average primitive solution state. */
    
    Wa = RoeAverage(Wl, Wr);

   /* Evaluate the jumps in the conserved solution states. */

    dUrl = Wr.U()-Wl.U();

    /* Evaluate the left, right, and average state eigenvalues. */

    lambdas_l = Wl.lambda_x();
    lambdas_r = Wr.lambda_x();
    lambdas_a = Wa.lambda_x();

    /* Determine the intermediate state flux. */

    wavespeed_l = min(lambdas_l[1],
                      lambdas_a[1]);
    wavespeed_r = max(lambdas_r[NUM_VAR_EULER3D],
                      lambdas_a[NUM_VAR_EULER3D]);
 
    if (wavespeed_l >= ZERO) {
        Flux = Wl.F();
    } else if (wavespeed_r <= ZERO) {
        Flux = Wr.F();
    } else {
        Flux =   ( (wavespeed_r*Wl.Fx()-wavespeed_l*Wr.Fx())
                  +(wavespeed_l*wavespeed_r)*dUrl)/
                 (wavespeed_r-wavespeed_l);
    } /* endif */

    /* Return solution flux. */

    return (Flux);

}

Euler3D_cState FluxHLLE_x(const Euler3D_cState &Ul,
	      	          const Euler3D_cState &Ur) {
   return (FluxHLLE_x(Ul.W(), Ur.W()));
}


/*********************************************************
 * Routine: FluxHLLE_n (Harten-Lax-van Leer flux         *
 *                      function, n-direction)           *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for an arbitrary direction defined by a unit     *
 * normal vector in the direction of interest, given     *
 * left and right solution states.  The problem is       *
 * solved by first applying a frame rotation to rotate   *
 * the problem to a local frame aligned with the unit    *
 * normal vector and then by using the so-called         *
 * Harten-Lax-van Leer approximation to specify the      *
 * intermediate state fluxes in terms of the rotated     *
 * solution states.  See Harten, Lax, van Leer (1983).   *
 *                                                       *
 *********************************************************/
Euler3D_cState FluxHLLE_n(const Euler3D_pState &Wl,
	      	          const Euler3D_pState &Wr,
                          const Vector3D &norm_dir) {

   double sin_beta, cos_beta, sin_alpha, cos_alpha;
   
  //solnvec in  Wl_rotated (Wr_rotated) is allocated using new 
   Euler3D_pState Wl_rotated(0), Wr_rotated(0);
   Euler3D_cState Flux, Flux_rotated;

    /* Determine the direction cosine's for the frame
       rotation. */

   sin_beta = norm_dir.z;
   cos_beta = sqrt(sqr(norm_dir.x)+sqr(norm_dir.y));
   if(cos_beta == ZERO){
      sin_alpha = ZERO;
      cos_alpha = ZERO;
   }else{
      sin_alpha = norm_dir.y/cos_beta;
      cos_alpha = norm_dir.x/cos_beta;
   }
     
   /* Apply the frame rotation and evaluate left and right
      solution states in the local rotated frame defined
      by the unit normal vector. */

   Wl_rotated.solnvec[0] = Wl.d();
   Wl_rotated.solnvec[1] = (Wl.vx()*cos_alpha +
                     Wl.vy()*sin_alpha)*cos_beta+ Wl.vz()*sin_beta;
   Wl_rotated.solnvec[2] =  (Wl.vx()*cos_alpha +
                      Wl.vy()*sin_alpha)*sin_beta+ Wl.vz()*cos_beta;
   Wl_rotated.solnvec[3] = ZERO;
   Wl_rotated.solnvec[4] = Wl.p();
    
   Wr_rotated.solnvec[0] = Wr.d();
   Wr_rotated.solnvec[1] = (Wr.vx()*cos_alpha +
                     Wr.vy()*sin_alpha)*cos_beta+ Wr.vz()*sin_beta;
   Wr_rotated.solnvec[2] =  (Wr.vx()*cos_alpha +
                      Wr.vy()*sin_alpha)*sin_beta+ Wr.vz()*cos_beta;
    
   Wr_rotated.solnvec[3] = ZERO;
   Wr_rotated.solnvec[4] = Wr.p();

    /* Evaluate the intermediate state solution 
       flux in the rotated frame. */

    Flux_rotated = FluxHLLE_x(Wl_rotated, Wr_rotated);

    /* Rotate back to the original Cartesian reference
       frame and return the solution flux. */

    Flux.d = Flux_rotated.d;
    Flux.dv.x = (Flux_rotated.dv.x*cos_beta + Flux_rotated.dv.y*sin_beta)*cos_alpha;
    Flux.dv.y = (Flux_rotated.dv.x*cos_beta + Flux_rotated.dv.y*sin_beta)*sin_alpha;
    Flux.dv.z = (Flux_rotated.dv.x*sin_beta + Flux_rotated.dv.y* cos_beta);
    
    Flux.E = Flux_rotated.E;

    return (Flux);

}

Euler3D_cState FluxHLLE_n(const Euler3D_cState &Ul,
	      	          const Euler3D_cState &Ur,
                          const Vector3D &norm_dir) {
    return (FluxHLLE_n(Ul.W(), Ur.W(), norm_dir));
}


/********************************************************
 * Routine: Reflect                                     *
 *                                                      *
 * This function returns the reflected solution state   *
 * in a given direction given the primitive solution    *
 * variables and the unit normal vector in the          *
 * direction of interest.                               *
 *                                                      *
 ********************************************************/
Euler3D_pState Reflect(const Euler3D_pState &W,
	      	       const Vector3D &norm_dir) {
   
   double sin_beta, cos_beta, sin_alpha, cos_alpha;
   double dr, ur, vr, wr, pr, u, v, w;
   
   
   sin_beta = norm_dir.z;
   cos_beta = sqrt(sqr(norm_dir.x)+sqr(norm_dir.y));
   if(cos_beta == ZERO){
      sin_alpha = ZERO;
      cos_alpha = ZERO;
   }else{
      sin_alpha = norm_dir.y/cos_beta;
      cos_alpha = norm_dir.x/cos_beta;
   }
     
   /* Apply the frame rotation and calculate the primitive
      solution state variables in the local rotated frame
      defined by the unit normal vector. */

   dr = W.d();
   ur = (W.vx()*cos_alpha + W.vy()*sin_alpha)*cos_beta+ W.vz()*sin_beta;
   vr = (W.vx()*cos_alpha + W.vy()*sin_alpha)*sin_beta+ W.vz()*cos_beta;
   wr = ZERO;
   pr = W.p();
   
   /* Reflect the normal velocity in the rotated frame. */
   ur = -ur;
   
   /* Rotate back to the original Cartesian reference frame. */

    u = (ur*cos_beta + vr*sin_beta)*cos_alpha;
    v = (ur*cos_beta + vr*sin_beta)*sin_alpha;
    w = (ur*sin_beta + vr*cos_beta);
    
    /* Return the reflected state. */

    return (Euler3D_pState(dr, u, v, w, pr));
       
}

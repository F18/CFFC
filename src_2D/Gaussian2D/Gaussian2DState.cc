/* Gaussian2DState.cc:  Subroutines for 2D Gaussian Solution State Classes. */

/* Include 2D Gaussian solution state header file. */

#ifndef _Gaussian2D_STATE_INCLUDED
#include "Gaussian2DState.h"
#endif // _Gaussian2D_STATE_INCLUDED

#ifndef _GAUSSIAN2D_INPUT_INCLUDED
#include "Gaussian2DInput.h"
#endif // _GAUSSIAN2D_INPUT_INCLUDED

/****************************************************************
 * Gaussian2D_pState -- Create storage and assign gas constants.*
 ****************************************************************/
double Gaussian2D_pState::M      = MOLE_WT_AIR;
int    Gaussian2D_pState::atoms  = GAUSSIAN_DIATOMIC;
int    Gaussian2D_pState::gas    = GAS_AIR;
double Gaussian2D_pState::alpha_m  = ONE;
double Gaussian2D_pState::omega  = OMEGA_AIR;
double Gaussian2D_pState::mu_not = MU_NOT_AIR;
double Gaussian2D_pState::pr     = 0.6666666666666666667;
/****************************************************************
 * Gaussian2D_cState -- Create storage and assign gas constants.*
 ****************************************************************/
double Gaussian2D_cState::M     = MOLE_WT_AIR;
int    Gaussian2D_cState::atoms = GAUSSIAN_DIATOMIC;
int    Gaussian2D_cState::gas   = GAS_AIR;
double Gaussian2D_cState::alpha_m = ONE;
double Gaussian2D_cState::pr    = 0.6666666666666666667;

/*************************************************************
 * Gaussian2D_pState -- set_state_from_ips                   *
 *************************************************************/
void Gaussian2D_pState::set_state_from_ips(Gaussian2D_Input_Parameters &IP) {
  setgas(IP.Gas_Type);
  set_temperature_d(IP.Temperature);
  v.x = IP.Mach_Number*sound()*cos(TWO*PI*IP.Flow_Angle/360.00);
  v.y = IP.Mach_Number*sound()*sin(TWO*PI*IP.Flow_Angle/360.00);
  return;
}

/********************************************************
 * Routine: RoeAverage (Roe Averages)                   *
 *                                                      *
 * This function returns the Roe-averaged primitive     *
 * solution state given left and right primitive        *
 * solution variables.  See Roe (1981).                 *
 *                                                      *
 ********************************************************/
Gaussian2D_pState RoeAverage(const Gaussian2D_pState &Wl,
	      	             const Gaussian2D_pState &Wr) {

    Tensor2D hl, hr;
    double sdl, sdr;
    double hxxa, hxya, hyya, hzza;
    double da, ua, va, pxxa, pxya, pyya, pzza, erota;
    Gaussian2D_pState Ave;

    // Determine the left and right state specific enthalpies
    // and square roots of the density. 

    hl = Wl.h();
    hr = Wr.h();
    sdl = sqrt(Wl.d);
    sdr = sqrt(Wr.d);

    // Determine the appropriate Roe averages.

    da    = sdl*sdr;
    ua    = (sdl*Wl.v.x+sdr*Wr.v.x)/(sdl+sdr);
    va    = (sdl*Wl.v.y+sdr*Wr.v.y)/(sdl+sdr);
    hxxa  = (sdl*hl.xx+sdr*hr.xx)/(sdl+sdr);
    hxya  = (sdl*hl.xy+sdr*hr.xy)/(sdl+sdr);
    hyya  = (sdl*hl.yy+sdr*hr.yy)/(sdl+sdr);
    hzza  = (sdl*hl.zz+sdr*hr.zz)/(sdl+sdr);
    erota = (sdl*Wl.erot+sdr*Wr.erot)/(sdl+sdr);

    pxxa = (hxxa-0.5*ua*ua)*da*2.0/3.0;
    pxya = (hxya-0.5*ua*va)*da*2.0/3.0;
    pyya = (hyya-0.5*va*va)*da*2.0/3.0;
    pzza = hzza*da*2.0/3.0;

    Ave = Gaussian2D_pState(da, ua, va, pxxa, pxya, pyya, pzza, erota);

    // Return the Roe-averged state.

    return (Ave);
       
}

/*********************************************************
 * Routine: Rotate                                       *
 *                                                       *
 * This function returns the solution in the lcoal       *
 * rotated frame.                                        *
 *                                                       *
 *********************************************************/

Gaussian2D_pState Rotate(const Gaussian2D_pState &W,
                         const Vector2D &norm_dir) {
  Gaussian2D_pState W_rotated;
  double cos_angle = norm_dir.x;  
  double sin_angle = norm_dir.y;

  W_rotated.d   = W.d;
  W_rotated.v.x =   W.v.x*cos_angle + W.v.y*sin_angle;
  W_rotated.v.y = - W.v.x*sin_angle + W.v.y*cos_angle;

  //Rotate Pressure tensor

  W_rotated.p.xx = W.p.xx*cos_angle*cos_angle+W.p.yy*sin_angle*sin_angle
                   +2.0*W.p.xy*cos_angle*sin_angle;
  W_rotated.p.xy = W.p.xy*(cos_angle*cos_angle-sin_angle*sin_angle)
                   -(W.p.xx-W.p.yy)*cos_angle*sin_angle;
  W_rotated.p.yy = W.p.xx*sin_angle*sin_angle+W.p.yy*cos_angle*cos_angle
                   -2.0*W.p.xy*cos_angle*sin_angle;
  W_rotated.p.zz = W.p.zz;
  W_rotated.erot = W.erot;
 
  return (W_rotated);

}

Gaussian2D_cState Rotate(const Gaussian2D_cState &U,
                         const Vector2D &norm_dir) {
  Gaussian2D_cState U_rotated;
  double cos_angle = norm_dir.x;  
  double sin_angle = norm_dir.y;

  U_rotated.d   = U.d;
  U_rotated.dv.x =   U.dv.x*cos_angle + U.dv.y*sin_angle;
  U_rotated.dv.y = - U.dv.x*sin_angle + U.dv.y*cos_angle;

  //Rotate energy tensor

  U_rotated.E.xx = U.E.xx*cos_angle*cos_angle+U.E.yy*sin_angle*sin_angle
                   +2.0*U.E.xy*cos_angle*sin_angle;
  U_rotated.E.xy = U.E.xy*(cos_angle*cos_angle-sin_angle*sin_angle)
                   -(U.E.xx-U.E.yy)*cos_angle*sin_angle;
  U_rotated.E.yy = U.E.xx*sin_angle*sin_angle+U.E.yy*cos_angle*cos_angle
                   -2.0*U.E.xy*cos_angle*sin_angle;
  U_rotated.E.zz = U.E.zz;
  U_rotated.erot = U.erot;
 
  return (U_rotated);

}

/**********************************************************************
 * Routine: Translate                                                 *
 *                                                                    *
 * This function returns the solution in a stationary frame.          *
 *                                                                    *
 **********************************************************************/
Gaussian2D_pState Translate(const Gaussian2D_pState &W, const Vector2D &V) {

  return Gaussian2D_pState(W.d,W.v.x-V.x,W.v.y-V.y,W.p,W.erot);

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
Gaussian2D_pState Reflect(const Gaussian2D_pState &W,
	      	          const Vector2D &norm_dir) {

    double dr, ur, vr, pxxr, pxyr, pyyr, pzzr, erotr;
    double u, v, pxx, pxy, pyy;
    Tensor2D pr;
    double cos_angle, sin_angle;

    // Determine the direction cosine's for the frame
    // rotation.

    cos_angle = norm_dir.x; 
    sin_angle = norm_dir.y;

    // Apply the frame rotation and calculate the primitive
    // solution state variables in the local rotated frame
    // defined by the unit normal vector. 

    dr    = W.d;
    ur    = W.v.x*cos_angle + W.v.y*sin_angle;
    vr    = - W.v.x*sin_angle + W.v.y*cos_angle;
    pxxr  = W.p.xx*cos_angle*cos_angle+W.p.yy*sin_angle*sin_angle
            +2.0*W.p.xy*cos_angle*sin_angle;
    pxyr  = W.p.xy*(cos_angle*cos_angle-sin_angle*sin_angle)
            -(W.p.xx-W.p.yy)*cos_angle*sin_angle;
    pyyr  = W.p.xx*sin_angle*sin_angle+W.p.yy*cos_angle*cos_angle
            -2.0*W.p.xy*cos_angle*sin_angle;
    pzzr  = W.p.zz;
    erotr = W.erot;

    // Reflect the normal velocity and shear pressure in the rotated frame.

    ur = -ur;
    pxyr = -pxyr;

    // Rotate back to the original Cartesian reference frame.

    sin_angle = -sin_angle;

    u   = ur*cos_angle + vr*sin_angle;
    v   = - ur*sin_angle + vr*cos_angle;
    pxx = pxxr*cos_angle*cos_angle+pyyr*sin_angle*sin_angle
          +2.0*pxyr*cos_angle*sin_angle;
    pxy = pxyr*(cos_angle*cos_angle-sin_angle*sin_angle)
          -(pxxr-pyyr)*cos_angle*sin_angle;
    pyy = pxxr*sin_angle*sin_angle+pyyr*cos_angle*cos_angle
          -2.0*pxyr*cos_angle*sin_angle;

    // Return the reflected state.

    return (Gaussian2D_pState(dr, u, v, pxx, pxy, pyy, pzzr, erotr));
       
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
Gaussian2D_pState Reflect(const Gaussian2D_pState &W,
			  const Vector2D &norm_dir,
			  const Vector2D &V) {

    double dr, ur, vr, pxxr, pxyr, pyyr, pzzr, erotr;
    double u, v, pxx, pxy, pyy;
    Tensor2D pr;
    double cos_angle, sin_angle;
    Vector2D Vr;

    /* Determine the direction cosine's for the frame
       rotation. */

    cos_angle = norm_dir.x; 
    sin_angle = norm_dir.y;

    /* Apply the frame rotation and calculate the primitive
       solution state variables in the local rotated frame
       defined by the unit normal vector. */

    dr    = W.d;
    ur    = W.v.x*cos_angle + W.v.y*sin_angle;
    vr    = - W.v.x*sin_angle + W.v.y*cos_angle;
    pxxr  = W.p.xx*cos_angle*cos_angle+W.p.yy*sin_angle*sin_angle
            +2.0*W.p.xy*cos_angle*sin_angle;
    pxyr  = W.p.xy*(cos_angle*cos_angle-sin_angle*sin_angle)
            -(W.p.xx-W.p.yy)*cos_angle*sin_angle;
    pyyr  = W.p.xx*sin_angle*sin_angle+W.p.yy*cos_angle*cos_angle
            -2.0*W.p.xy*cos_angle*sin_angle;
    pzzr  = W.p.zz;
    erotr = W.erot;

    Vr.Rotate(V,norm_dir);

    /* Reflect the normal velocity in the rotated frame. */

    ur = -ur + TWO*Vr.x;
    vr =  vr + TWO*Vr.y;
    pxyr = -pxyr;

    /* Rotate back to the original Cartesian reference frame. */

    sin_angle = -sin_angle;

    u   = ur*cos_angle + vr*sin_angle;
    v   = - ur*sin_angle + vr*cos_angle;
    pxx = pxxr*cos_angle*cos_angle+pyyr*sin_angle*sin_angle
          +2.0*pxyr*cos_angle*sin_angle;
    pxy = pxyr*(cos_angle*cos_angle-sin_angle*sin_angle)
          -(pxxr-pyyr)*cos_angle*sin_angle;
    pyy = pxxr*sin_angle*sin_angle+pyyr*cos_angle*cos_angle
          -2.0*pxyr*cos_angle*sin_angle;

    /* Return the reflected state. */

    return (Gaussian2D_pState(dr, u, v, pxx, pxy, pyy, pzzr, erotr));

       
}

/********************************************************
 * Routine: Adiabatic_Wall                              *
 *                                                      *
 * This function returns the solution state for a wall  *
 * in a given direction given the primitive solution    *
 * variables and the unit normal vector in the          *
 * direction of interest.                               *
 *                                                      *
 ********************************************************/
Gaussian2D_pState Adiabatic_Wall(const Gaussian2D_pState &W,
				 const Vector2D &V,
				 const Vector2D &norm_dir) {

    double vWallr;
    Gaussian2D_pState W_rot, W_rot2, W2;

    vWallr =  - V.x*norm_dir.y + V.y*norm_dir.x;

    W_rot = Rotate(W,norm_dir);
    W_rot2 = W_rot;

    // Adjust state.

    //  Mine 1
    /*
    ur2   = -ur;
    vr2   = vr+W.alpha_m*(vWallr-vr)+2.0*(2.0-W.alpha_m)*pxyr/sqrt(2*PI*dr*pxxr);
    pxxr2 = pxxr;
    pxyr2 = (W.alpha_m-1.0)*pxyr+2.0*W.alpha_m*sqrt(dr*pxxr/(2.0*PI))*(vr-vWallr);
    pyyr2 = pyyr - 2.0/3.0*W.d*sqr(W.alpha_m/2.0*(vWallr-vr)+(2.0-W.alpha_m)*pxyr/sqrt(2*PI*dr*pxxr))
                 + 2.0*((2.0-W.alpha_m)*(dr/2.0*sqr(vr-vr2)+pxyr*sqrt(2.0*dr/(PI*pxxr))*(vr-vr2))
                         +W.alpha_m*dr/2.0*sqr(vWallr-vr2));
    */
    //  Mine 2
    
    W_rot2.v.x  = -W_rot.v.x;
    W_rot2.v.y  = W_rot.v.y;
    //W_rot2.p.xx = W_rot.p.xx;
    W_rot2.p.xy = (W.alpha_m-1.0)*W_rot.p.xy+2.0*W.alpha_m*sqrt(W_rot.d*W_rot.p.xx/(2.0*PI))*(W_rot.v.y-vWallr);
    //W_rot2.p.yy = W_rot.p.yy;

    //  Dr. Groth's
    /*                
    ur2   = -ur;
    vr2   = vr;
    pxxr2 = pxxr;
    pxyr2 = -pxyr+4.0*(W.alpha_m)/((2.0-W.alpha_m)*sqrt(2.0*PI))*dr*(vr-vWallr)
                                                     *sqrt(pxxr/dr+ur*ur/3);
    pyyr2 = pyyr;
    */

    // Rotate back to the original Cartesian reference frame.

    W2 = Rotate(W_rot2,Vector2D(norm_dir.x,-norm_dir.y));

    // Return the reflected state for adiabatic wall.
    return W2;
       
}

/********************************************************
 * Routine: Isothermal_Wall                             *
 *                                                      *
 * This function returns the solution state for a wall  *
 * in a given direction given the primitive solution    *
 * variables and the unit normal vector in the          *
 * direction of interest.                               *
 *                                                      *
 ********************************************************/
Gaussian2D_pState Isothermal_Wall(const Gaussian2D_pState &W,
				  const Vector2D &V,
				  const double &T,
				  const Vector2D &norm_dir) {
    double vWallr;
    double m = W.M/(THOUSAND*AVOGADRO);  //particle mass
    double ng, nw, n_temp, pxy_desired, v_temp;
    Gaussian2D_pState W_rot, W_rot2, W2;

    vWallr =  - V.x*norm_dir.y + V.y*norm_dir.x;

    W_rot = Rotate(W,norm_dir);
    W_rot2 = W_rot;

    // Adjust state.

    ng = W_rot.d/m;   //number density of "incoming" Gaussian
    nw = ng*sqrt(W_rot.p.xx/(ng*BOLTZMANN*T));   //number density of reflected Maxwellian to ensure
                                                 //no mass transfer (thermal equilibrium with wall)
    n_temp = ng+W.alpha_m/2.0*(nw-ng);        //"total" number density
  
    W_rot2.d = sqr(n_temp*m)/W_rot.d;  //want roe-average state to be n_temp*m

    v_temp = (m*ng*(2.0-W.alpha_m)*(W_rot.v.y/2.0+W_rot.p.xy/sqrt(2.0*PI*W_rot.d*W_rot.p.xx))
	      +W.alpha_m/2.0*m*nw*vWallr)/(n_temp*m);  //combined momentum devided by density


    pxy_desired = W.alpha_m*(W_rot.p.xy/2.0+sqrt(W_rot.d*W_rot.p.xx/(2.0*PI))*(W_rot.v.y-v_temp)
			   -sqrt(nw*nw*m*BOLTZMANN*T/(2.0*PI))*(vWallr-v_temp));

    W_rot2.v.x  = -W_rot.v.x*sqrt(W_rot.d/W_rot2.d);
    //W_rot2.v.y  = W_rot.v.y;
    //W_rot2.p.xx = W_rot.p.xx;
    W_rot2.p.xy = (pxy_desired*(sqrt(W_rot.d)+sqrt(W_rot2.d))/(n_temp*m)-W_rot.p.xy/sqrt(W_rot.d))*sqrt(W_rot2.d);
    //W_rot2.p.yy = W_rot.p.yy;

    // Rotate back to the original Cartesian reference frame.

    W2 = Rotate(W_rot2,Vector2D(norm_dir.x,-norm_dir.y));

    return W2;
}

/********************************************************
 * Routine: Knudsen_Layer_Adiabatic                     *
 *                                                      *
 * This function returns the state inside a knudsen     *
 * layer.  This is different than the "Adiabatic_Wall"  *
 * function which returns the state required to give    *
 * the correct "Roe average" for flux calculations.     *
 * This is the "adiabatic" version of this function     *
 * (the emitted Maxwellian is at the same temperature   *
 * as the normal pressure of the surrounding            *
 * distribution.  I realize now this is not strictly    *
 * adiabatic.                                           *
 *                                                      *
 ********************************************************/
Gaussian2D_pState Knudsen_Layer_Adiabatic(const Gaussian2D_pState &W,
					  const Vector2D  &v,
					  const Vector2D &norm_dir) {
  Gaussian2D_pState W_rot, Kn_rot;
  double v_temp, pxy_temp, pyy_temp;

  double vy_wall_rot = -v.x*norm_dir.x + v.y*norm_dir.y;

  W_rot = Rotate(W,norm_dir);

  v_temp = W_rot.v.y+W.alpha_m/2.0*(vy_wall_rot-W_rot.v.y)
           +(2.0-W.alpha_m)*(W_rot.p.xy/sqrt(2*PI*W_rot.d*W_rot.p.xx));

  pxy_temp =  W.alpha_m*(W_rot.p.xy/2.0
		       +sqrt(W_rot.d*W_rot.p.xx/(2.0*PI))*(W_rot.v.y-vy_wall_rot));

  pyy_temp = W_rot.p.yy+(2.0-W.alpha_m)*(W_rot.d/2.0*sqr(W_rot.v.y-v_temp)
				       +W_rot.p.xy*sqrt(2.0*W_rot.d/(PI*W_rot.p.yy))*(W_rot.v.y-v_temp))
             +W.alpha_m*W_rot.d/2.0*sqr(vy_wall_rot-v_temp);

  Kn_rot = Gaussian2D_pState(W_rot.d,
			     0.0,
			     v_temp,
			     W_rot.p.xx,
			     pxy_temp,
			     pyy_temp,
			     W_rot.p.zz,
			     W_rot.erot);

  return Rotate(Kn_rot,Vector2D(norm_dir.x,-norm_dir.y));

}

/********************************************************
 * Routine: Knudsen_Layer_Isothermal                    *
 *                                                      *
 * This function returns the state inside a knudsen     *
 * layer.  This is different than the "Adiabatic_Wall"  *
 * function which returns the state required to give    *
 * the correct "Roe average" for flux calculations.     *
 * This is the "isothermal" version of this function    *
 *                                                      *
 ********************************************************/
Gaussian2D_pState Knudsen_Layer_Isothermal(const Gaussian2D_pState &W,
					   const Vector2D &V,
					   const double &T,
					   const Vector2D &norm_dir) {
  Gaussian2D_pState W_rot, Kn_rot;
  double m = W.M/(THOUSAND*AVOGADRO);  //particle mass
  double vy_wall_rot = -V.x*norm_dir.x + V.y*norm_dir.y;
  double ng, nw, n_temp, rho_temp, v_temp, pxx_temp, pxy_temp, pyy_temp, pzz_temp, erot_temp;

  W_rot = Rotate(W,norm_dir);

  ng = W_rot.d/m;   //number density of "incoming" Gaussian

  nw = ng*sqrt(W_rot.p.xx/(ng*BOLTZMANN*T));   //number density of reflected Maxwellian to ensure
                                               //no mass transfer (thermal equilibrium with wall)

  n_temp = ng+W.alpha_m/2.0*(nw-ng);           //"total" number density

  rho_temp = n_temp*m;

  v_temp = (m*ng*(2.0-W.alpha_m)*(W_rot.v.y/2.0+W_rot.p.xy/sqrt(2.0*PI*W_rot.d*W_rot.p.xx))
	    +W.alpha_m/2.0*m*nw*vy_wall_rot)/rho_temp;  //combined momentum devided by density

  pxx_temp = W_rot.p.xx+W.alpha_m/2.0*(nw*BOLTZMANN*T-W_rot.p.xx);

  pxy_temp = W.alpha_m*(W_rot.p.xy/2.0+sqrt(W_rot.d*W_rot.p.xx/(2.0*PI))*(W_rot.v.y-v_temp)
		      -sqrt(nw*nw*m*BOLTZMANN*T/(2.0*PI))*(vy_wall_rot-v_temp));

  pyy_temp = (2.0-W.alpha_m)*(W_rot.p.yy/2.0+W_rot.d*sqr(W_rot.v.y-v_temp)/2.0
			    +W_rot.p.xy*sqrt(2.0*W_rot.d/PI*W_rot.p.xx)*(W_rot.v.y-v_temp))
             +W.alpha_m*(nw*BOLTZMANN*T/2.0 + nw*m/2.0*sqr(vy_wall_rot-v_temp));

  pzz_temp = W_rot.p.zz+W.alpha_m/2.0*(nw*BOLTZMANN*T-W_rot.p.zz);

  erot_temp = W_rot.erot+W.alpha_m/2.0*(nw*BOLTZMANN*T-W_rot.erot);

//  rho_temp = W_rot.d*(1.0+W.alpha_m/2.0*(sqrt(W_rot.p.xx/(n*BOLTZMANN*T))-1.0));
//
//
//  v_temp = (2.0-W.alpha_m)*(W_rot.v.y/2.0+W_rot.p.xy/(sqrt(2.0*PI*W_rot.d*W_rot.p.xx)))
//           +W.alpha_m/2.0*sqrt(W_rot.p.xx/n*BOLTZMANN*T)*vy_wall_rot;
//
//  pxx_temp = (2.0-W.alpha_m)/2.0*W_rot.p.xx+W.alpha_m/2.0*sqrt(W_rot.p.xx*n*BOLTZMANN*T);
//
//
//  pxy_temp = W.alpha_m*(W+rot.p.xy/2+sqrt(W_rot.d*W_rot.p.xx/(2.0*PI))*(W_rot.v.y-v_temp)
//  		      +sqrt(;
//
//  pyy_temp = 0.0;
//
//  pzz_temp = (2.0-W.alpha_m)/2.0*W_rot.p.zz+W.alpha_m/2.0*sqrt(W_rot.p.xx*n*BOLTZMANN*T);
//
//  erot_temp = (2.0-W.alpha_m)/2.0*W_rot.erot+W.alpha_m/2.0*sqrt(W_rot.p.xx*n*BOLTZMANN*T);

  Kn_rot = Gaussian2D_pState(rho_temp,
			     0.0,
			     v_temp,
			     pxx_temp,
			     pxy_temp,
			     pyy_temp,
			     pzz_temp,
			     erot_temp);

  return Rotate(Kn_rot,Vector2D(norm_dir.x,-norm_dir.y));

}

/********************************************************
 * Routine: NoSlip                                      *
 *                                                      *
 * This function returns the no-slip solution state     *
 * in a given direction given the primitive solution    *
 * variables and the unit normal vector in the          *
 * direction of interest.                               *
 *                                                      *
 ********************************************************/
Gaussian2D_pState NoSlip(const Gaussian2D_pState &W,
	      	         const Vector2D &norm_dir) {

    double dr, ur, vr, u, v, erot;
    Tensor2D pr;
    double cos_angle, sin_angle;

    // Determine the direction cosine's for the frame
    // rotation.

    cos_angle = norm_dir.x; 
    sin_angle = norm_dir.y;

    // Apply the frame rotation and calculate the primitive
    // solution state variables in the local rotated frame
    // defined by the unit normal vector.

    dr = W.d;
    ur = W.v.x*cos_angle +
         W.v.y*sin_angle;
    vr = - W.v.x*sin_angle +
           W.v.y*cos_angle;
    pr = W.p;
    erot = W.erot;

    // Reflect the normal and tangential velocity components
    // in the rotated frame.

    ur = -ur;
    vr = -vr;

    // Rotate back to the original Cartesian reference frame.

    u = ur*cos_angle - vr*sin_angle;
    v = ur*sin_angle + vr*cos_angle;

    // Return the no-slip state. 

    return (Gaussian2D_pState(dr, u, v, pr, erot));
       
}

/********************************************************
 * Routine: BC_Characteristic_Pressure                  *
 *   (Characteristic-Based Boundary Condition with      *
 *    Static Pressure Specified Whenever Possible)      *
 *                                                      *
 * This function returns the boundary solution state    *
 * for a given direction given the primitive solution   *
 * state on the interior of the boundary, Wi, the       *
 * desired flow state to be imposed at the boundary,    *
 * Wo, and the unit normal vector in the direction of   *
 * interest. A simplified characteristic analysis is    *
 * used to specify the boundary flow state in which the *
 * static pressure is specified whenever possible.  The *
 * imposition of the boundary-data respects the         *
 * directions of propogation for the solution           *
 * characteristics at the boundary and thereby ensures  *
 * that the boundary conditions are not ill-posed (i.e.,*
 * the boundary data is not under- or over-determined). *
 * The following procedure is adopted:                  *
 *                                                      *
 ********************************************************/

Gaussian2D_pState BC_Characteristic_Pressure(const Gaussian2D_pState &Wi,
                                             const Gaussian2D_pState &Wo,
	      	                             const Vector2D &norm_dir) {

    Gaussian2D_pState Wi_rotated, Wo_rotated;
    double mi, db, ub, vb, pxxb, pxyb, pyyb, pzzb, erotb; 
    double ub_rotated, vb_rotated, pxxb_rotated, pxyb_rotated;
    double pyyb_rotated;
    double cos_angle, sin_angle;

    // Determine the direction cosine's for the frame
    // rotation. 

    cos_angle = norm_dir.x; 
    sin_angle = norm_dir.y;

    // Apply the frame rotation and evaluate interior and 
    // imposed boundary solution states in the local rotated 
    // frame defined by the unit normal vector.

    Wi_rotated.d = Wi.d;
    Wi_rotated.v.x = Wi.v.x*cos_angle +
                     Wi.v.y*sin_angle;
    Wi_rotated.v.y = - Wi.v.x*sin_angle +
                       Wi.v.y*cos_angle;
    Wi_rotated.p.xx = Wi.p.xx*cos_angle*cos_angle+Wi.p.yy*sin_angle*sin_angle
                     +2.0*Wi.p.xy*cos_angle*sin_angle;
    Wi_rotated.p.xy = Wi.p.xy*(cos_angle*cos_angle-sin_angle*sin_angle)
                     -(Wi.p.xx-Wi.p.yy)*cos_angle*sin_angle;
    Wi_rotated.p.yy = Wi.p.xx*sin_angle*sin_angle+Wi.p.yy*cos_angle*cos_angle
                     -2.0*Wi.p.xy*cos_angle*sin_angle;
    Wi_rotated.p.zz = Wi.p.zz;
    Wi_rotated.erot = Wi.erot;

    Wo_rotated.d = Wo.d;
    Wo_rotated.v.x = Wo.v.x*cos_angle +
                     Wo.v.y*sin_angle;
    Wo_rotated.v.y = - Wo.v.x*sin_angle +
                       Wo.v.y*cos_angle;
    Wo_rotated.p.xx = Wo.p.xx*cos_angle*cos_angle+Wo.p.yy*sin_angle*sin_angle
                     +2.0*Wo.p.xy*cos_angle*sin_angle;
    Wo_rotated.p.xy = Wo.p.xy*(cos_angle*cos_angle-sin_angle*sin_angle)
                     -(Wo.p.xx-Wo.p.yy)*cos_angle*sin_angle;
    Wo_rotated.p.yy = Wo.p.xx*sin_angle*sin_angle+Wo.p.yy*cos_angle*cos_angle
                     -2.0*Wo.p.xy*cos_angle*sin_angle;
    Wo_rotated.p.zz = Wo.p.zz;
    Wo_rotated.erot = Wo.erot;


    // Determine the Mach number at the interior node

    mi = Wi_rotated.v.x/Wi_rotated.axx();

    // Boundary condition for supersonic outflow.

    if (mi >= sqrt(3.0)) {
       db = Wi_rotated.d;
       ub_rotated = Wi_rotated.v.x;
       vb_rotated = Wi_rotated.v.y;
       pxxb_rotated = Wi_rotated.p.xx;
       pxyb_rotated = Wi_rotated.p.xy;
       pyyb_rotated = Wi_rotated.p.yy;
       pzzb = Wi_rotated.p.zz;
       erotb = Wi_rotated.erot;

    // Boundary condition for subsonic outflow. 
    // Pressure specified.
    } else if (mi >= ONE) {
       db = Wi_rotated.d;
       ub_rotated = Wi_rotated.v.x;
       vb_rotated = Wi_rotated.v.y;
       pxxb_rotated = Wo_rotated.p.xx;
       pxyb_rotated = Wi_rotated.p.xy;
       pyyb_rotated = Wi_rotated.p.yy;
       pzzb = Wi_rotated.p.zz;
       erotb = Wi_rotated.erot;

    } else if (mi >= ZERO) {
       db = Wi_rotated.d;
       ub_rotated = Wi_rotated.v.x;
       vb_rotated = Wo_rotated.v.y;
       pxxb_rotated = Wo_rotated.p.xx;
       pxyb_rotated = Wi_rotated.p.xy;
       pyyb_rotated = Wi_rotated.p.yy;
       pzzb = Wi_rotated.p.zz;
       erotb = Wi_rotated.erot;

    // Boundary condition for subsonic inflow. 
    // Pressure specified.

    } else if (mi >= -ONE) {
       db = Wo_rotated.d;
       ub_rotated = Wi_rotated.v.x;
       vb_rotated = Wo_rotated.v.y;
       pxxb_rotated = Wo_rotated.p.xx;
       pxyb_rotated = Wi_rotated.p.xy;
       pyyb_rotated = Wo_rotated.p.yy;
       pzzb = Wo_rotated.p.zz;
       erotb = Wo_rotated.erot;

    } else if (mi >= -sqrt(3.0)) {
       db = Wo_rotated.d;
       ub_rotated = Wi_rotated.v.x;
       vb_rotated = Wo_rotated.v.y;
       pxxb_rotated = Wo_rotated.p.xx;
       pxyb_rotated = Wo_rotated.p.xy;
       pyyb_rotated = Wo_rotated.p.yy;
       pzzb = Wo_rotated.p.zz;
       erotb = Wo_rotated.erot;

       //Supersonic inflow

    } else {
       db = Wo_rotated.d;
       ub_rotated = Wo_rotated.v.x;
       vb_rotated = Wo_rotated.v.y;
       pxxb_rotated = Wo_rotated.p.xx;
       pxyb_rotated = Wo_rotated.p.xy;
       pyyb_rotated = Wo_rotated.p.yy;
       pzzb = Wo_rotated.p.zz;
       erotb = Wo_rotated.erot;

    }// endif 

    // Rotate the resulting boundary state back to the original 
    // Cartesian reference frame.

    sin_angle = -sin_angle;

    ub = ub_rotated*cos_angle + vb_rotated*sin_angle;
    vb = - ub_rotated*sin_angle + vb_rotated*cos_angle;
    pxxb = pxxb_rotated*cos_angle*cos_angle+pyyb_rotated*sin_angle*sin_angle
                +2.0*pxyb_rotated*cos_angle*sin_angle;
    pxyb = pxyb_rotated*(cos_angle*cos_angle-sin_angle*sin_angle)
                -(pxxb_rotated-pyyb_rotated)*cos_angle*sin_angle;
    pyyb = pxxb_rotated*sin_angle*sin_angle+pyyb_rotated*cos_angle*cos_angle
                -2.0*pxyb_rotated*cos_angle*sin_angle;

    // Return boundary solution state.

    return (Gaussian2D_pState(db, ub, vb, pxxb,
                              pxyb, pyyb, pzzb, erotb));
       
}

/********************************************************
 * Routine: BC_Characteristic_Velocity                  *
 *   (Characteristic-Based Boundary Condition with      *
 *    Flow VElocity Specified Whenever Possible)        *
 *                                                      *
 * This function returns the boundary solution state    *
 * for a given direction given the primitive solution   *
 * state on the interior of the boundary, Wi, the       *
 * desired flow state to be imposed at the boundary,    *
 * Wo, and the unit normal vector in the direction of   *
 * interest. A simplified characteristic analysis is    *
 * used to specify the boundary flow state in which the *
 * flow Mach number is specified whenever possible.     *
 * The imposition of the boundary-data respects the     *
 * directions of propogation for the solution           *
 * characteristics at the boundary and thereby ensures  *
 * that the boundary conditions are not ill-posed (i.e.,*
 * the boundary data is not under- or over-determined). *

 ********************************************************/

Gaussian2D_pState BC_Characteristic_Velocity(const Gaussian2D_pState &Wi,
                                             const Gaussian2D_pState &Wo,
	      	                             const Vector2D &norm_dir) {

    Gaussian2D_pState Wi_rotated, Wo_rotated, Roe_Ave;
    double mi, db, ub, vb, pxxb, pxyb, pyyb, pzzb, erotb; 
    double ub_rotated, vb_rotated, pxxb_rotated, pxyb_rotated;
    double pyyb_rotated;
    double cos_angle, sin_angle;

    // Determine the direction cosine's for the frame
    // rotation. 

    cos_angle = norm_dir.x; 
    sin_angle = norm_dir.y;

    // Apply the frame rotation and evaluate interior and 
    // imposed boundary solution states in the local rotated 
    // frame defined by the unit normal vector.

    Wi_rotated.d = Wi.d;
    Wi_rotated.v.x = Wi.v.x*cos_angle +
                     Wi.v.y*sin_angle;
    Wi_rotated.v.y = - Wi.v.x*sin_angle +
                       Wi.v.y*cos_angle;
    Wi_rotated.p.xx = Wi.p.xx*cos_angle*cos_angle+Wi.p.yy*sin_angle*sin_angle
                     +2.0*Wi.p.xy*cos_angle*sin_angle;
    Wi_rotated.p.xy = Wi.p.xy*(cos_angle*cos_angle-sin_angle*sin_angle)
                     -(Wi.p.xx-Wi.p.yy)*cos_angle*sin_angle;
    Wi_rotated.p.yy = Wi.p.xx*sin_angle*sin_angle+Wi.p.yy*cos_angle*cos_angle
                     -2.0*Wi.p.xy*cos_angle*sin_angle;
    Wi_rotated.p.zz = Wi.p.zz;
    Wi_rotated.erot = Wi.erot;

    Wo_rotated.d = Wo.d;
    Wo_rotated.v.x = Wo.v.x*cos_angle +
                     Wo.v.y*sin_angle;
    Wo_rotated.v.y = - Wo.v.x*sin_angle +
                       Wo.v.y*cos_angle;
    Wo_rotated.p.xx = Wo.p.xx*cos_angle*cos_angle+Wo.p.yy*sin_angle*sin_angle
                     +2.0*Wo.p.xy*cos_angle*sin_angle;
    Wo_rotated.p.xy = Wo.p.xy*(cos_angle*cos_angle-sin_angle*sin_angle)
                     -(Wo.p.xx-Wo.p.yy)*cos_angle*sin_angle;
    Wo_rotated.p.yy = Wo.p.xx*sin_angle*sin_angle+Wo.p.yy*cos_angle*cos_angle
                     -2.0*Wo.p.xy*cos_angle*sin_angle;
    Wo_rotated.p.zz = Wo.p.zz;
    Wo_rotated.erot = Wo.erot;


    // Determine the Mach number at the interior node

    Roe_Ave = RoeAverage(Wi_rotated, Wo_rotated);

    mi = Roe_Ave.v.x/Roe_Ave.axx();

    // Boundary condition for supersonic outflow.

    if (mi >= sqrt(3.0)) {
       db = Wi_rotated.d;
       ub_rotated = Wi_rotated.v.x;
       vb_rotated = Wi_rotated.v.y;
       pxxb_rotated = Wi_rotated.p.xx;
       pxyb_rotated = Wi_rotated.p.xy;
       pyyb_rotated = Wi_rotated.p.yy;
       pzzb = Wi_rotated.p.zz;
       erotb = Wi_rotated.erot;

    // Boundary condition for subsonic outflow. 
    // Pressure specified.
    } else if (mi >= ONE) {
       db = Wi_rotated.d;
       ub_rotated = Wo_rotated.v.x;
       vb_rotated = Wi_rotated.v.y;
       pxxb_rotated = Wi_rotated.p.xx;
       pxyb_rotated = Wi_rotated.p.xy;
       pyyb_rotated = Wi_rotated.p.yy;
       pzzb = Wi_rotated.p.zz;
       erotb = Wi_rotated.erot;

    } else if (mi >= ZERO) {
       db = Wi_rotated.d;
       ub_rotated = Wo_rotated.v.x;
       vb_rotated = Wo_rotated.v.y;
       pxxb_rotated = Wi_rotated.p.xx;
       pxyb_rotated = Wi_rotated.p.xy;
       pyyb_rotated = Wi_rotated.p.yy;
       pzzb = Wi_rotated.p.zz;
       erotb = Wi_rotated.erot;

    // Boundary condition for subsonic inflow. 
    // Pressure specified.

    } else if (mi >= -ONE) {
       db = Wo_rotated.d;
       ub_rotated = Wo_rotated.v.x;
       vb_rotated = Wo_rotated.v.y;
       pxxb_rotated = Wi_rotated.p.xx;
       pxyb_rotated = Wi_rotated.p.xy;
       pyyb_rotated = Wo_rotated.p.yy;
       pzzb = Wo_rotated.p.zz;
       erotb = Wo_rotated.erot;

    } else if (mi >= -sqrt(3.0)) {
       db = Wo_rotated.d;
       ub_rotated = Wo_rotated.v.x;
       vb_rotated = Wo_rotated.v.y;
       pxxb_rotated = Wi_rotated.p.xx;
       pxyb_rotated = Wo_rotated.p.xy;
       pyyb_rotated = Wo_rotated.p.yy;
       pzzb = Wo_rotated.p.zz;
       erotb = Wo_rotated.erot;

       //Supersonic inflow

    } else {
       db = Wo_rotated.d;
       ub_rotated = Wo_rotated.v.x;
       vb_rotated = Wo_rotated.v.y;
       pxxb_rotated = Wo_rotated.p.xx;
       pxyb_rotated = Wo_rotated.p.xy;
       pyyb_rotated = Wo_rotated.p.yy;
       pzzb = Wo_rotated.p.zz;
       erotb = Wo_rotated.erot;

    }// endif 

    // Rotate the resulting boundary state back to the original 
    // Cartesian reference frame.

    sin_angle = -sin_angle;

    ub = ub_rotated*cos_angle + vb_rotated*sin_angle;
    vb = - ub_rotated*sin_angle + vb_rotated*cos_angle;
    pxxb = pxxb_rotated*cos_angle*cos_angle+pyyb_rotated*sin_angle*sin_angle
                +2.0*pxyb_rotated*cos_angle*sin_angle;
    pxyb = pxyb_rotated*(cos_angle*cos_angle-sin_angle*sin_angle)
                -(pxxb_rotated-pyyb_rotated)*cos_angle*sin_angle;
    pyyb = pxxb_rotated*sin_angle*sin_angle+pyyb_rotated*cos_angle*cos_angle
                -2.0*pxyb_rotated*cos_angle*sin_angle;

    // Return boundary solution state.

    return (Gaussian2D_pState(db, ub, vb, pxxb,
                              pxyb, pyyb, pzzb, erotb));
       
}

/********************************************************
 * Routine: BC_Couette                                  *
 *                                                      *
 * This BC is useful for couette problems with a scewed *
 * mesh or embedded boundaries.                         *
 *                                                      *
 ********************************************************/
Gaussian2D_pState BC_Couette(const Gaussian2D_pState &Wi,
			     const Gaussian2D_pState &Wo) {
  return Gaussian2D_pState(Wo.d,Wi.v.x,Wi.v.y,Wo.p.xx,
			   Wi.p.xy,Wo.p.yy,Wo.p.zz,Wo.erot);
}

/********************************************************
 * Routine: BC_Developed_Channel_Flow                   *
 *                                                      *
 * This function returns the boundary state             *
 * for fully developed channel flow,                    *
 * this is used for the branched duct case.             *
 *                                                      *
 ********************************************************/
Gaussian2D_pState BC_Developed_Channel_Flow(const Gaussian2D_pState &Wo,
					    const Vector2D &norm_dir,
					    const double &y) {

  //this gives a parabolic inlet with the average velocity
  //equal to Wo.v.x

  Gaussian2D_pState W_returned;
  double average_velocity, width;

  average_velocity = Wo.v.x;
  width = 0.1143;

  W_returned = Wo;

  W_returned.v.x = 6.0*average_velocity/(width*width)*(-y*y+width*y);
  W_returned.p.xy = -6.0*average_velocity/(width*width)*(-2.0*y+width)*Wo.viscosity();

  return W_returned;

}

/********************************************************
 * Routine: BCs (Boundary Conditions, x-direction)      *
 *                                                      *
 * This function returns the boundary state given       *
 * the boundary and interior states, the boundary type, *
 * and grid end type.  For Cartesian mesh only.         *
 *                                                      *
 ********************************************************/
Gaussian2D_pState BCs(const Gaussian2D_pState &Wb,
                      const Gaussian2D_pState &Wi,
		      const Gaussian2D_pState &dWdx,
		      const double &dx,
		      const int BC_type,
		      const int End_type) {

    switch(BC_type) {

      case BC_FIXED :
        return(Wb);
        break;

      case BC_CONSTANT_EXTRAPOLATION :
        return (Wi);
        break;

      case BC_LINEAR_EXTRAPOLATION :
        switch(End_type) {
          case LEFT_END_BOUNDARY :
            return (Wi-dWdx*dx);
            break;
          case RIGHT_END_BOUNDARY :
            return (Wi+dWdx*dx);
            break;
        default:
            return (Wi+dWdx*dx);
            break;
        } /* endswitch */
        break;

      case BC_REFLECTION :
        return (Gaussian2D_pState(Wi.d, -Wi.v.x, Wi.v.y, Wi.p.xx, -Wi.p.xy, Wi.p.yy, Wi.p.zz, Wi.erot));
        break;

      case BC_PERIODIC :
        return (Wb);
        break;

      case BC_CHARACTERISTIC :
        switch(End_type) {
          case LEFT_END_BOUNDARY :
            return (BC_Characteristic_Pressure(Wi,
                                               Wb,
	      	                               -Vector2D_NX));
            break;
          case RIGHT_END_BOUNDARY :
            return (BC_Characteristic_Pressure(Wi,
                                               Wb,
	      	                               Vector2D_NX));
            break;
          default:
            return (BC_Characteristic_Pressure(Wi,
                                               Wb,
	     	                               Vector2D_NX));
            break;
        } /* endswitch */
        break;

      default:
        return (Wb);
        break;

    } /* endswitch */

}

/********************************************************
 * Routine: BCs_x (Boundary Conditions, x-direction)    *
 *                                                      *
 * This function returns the boundary state given       *
 * the boundary and interior states, the boundary type, *
 * and grid end type.  For Cartesian mesh only.         *
 *                                                      *
 ********************************************************/
Gaussian2D_pState BCs_x(const Gaussian2D_pState &Wb,
                        const Gaussian2D_pState &Wi,
		        const Gaussian2D_pState &dWdx,
		        const double &dx,
		        const int BC_type,
		        const int End_type) {

    switch(BC_type) {

      case BC_FIXED :
        return (Wb);
        break;

      case BC_CONSTANT_EXTRAPOLATION :
        return (Wi);
        break;

      case BC_LINEAR_EXTRAPOLATION :
        switch(End_type) {
          case LEFT_END_BOUNDARY :
            return (Wi-dWdx*dx);
            break;
          case RIGHT_END_BOUNDARY :
            return (Wi+dWdx*dx);
            break;
          default:
            return (Wi+dWdx*dx);
            break;
        } /* endswitch */
        break;

      case BC_REFLECTION :
        return (Gaussian2D_pState(Wi.d, -Wi.v.x, Wi.v.y, Wi.p.xx, -Wi.p.xy, Wi.p.yy, Wi.p.zz, Wi.erot));
        break;

      case BC_PERIODIC :
        return (Wb);
        break;

      case BC_CHARACTERISTIC :
        switch(End_type) {
          case LEFT_END_BOUNDARY :
            return (BC_Characteristic_Pressure(Wi,
                                               Wb,
	     	                               -Vector2D_NX));
            break;
          case RIGHT_END_BOUNDARY :
            return (BC_Characteristic_Pressure(Wi,
                                               Wb,
	      	                               Vector2D_NX));
            break;
          default:
            return (BC_Characteristic_Pressure(Wi,
                                               Wb,
	      	                               Vector2D_NX));
            break;
        } /* endswitch */
        break;

      default:
        return (Wb);
        break;

    } /* endswitch */

}

/********************************************************
 * Routine: BCs_y (Boundary Conditions, y-direction)    *
 *                                                      *
 * This function returns the boundary state given       *
 * the boundary and interior states, the boundary type, *
 * and grid end type.  For Cartesian mesh only.         *
 *                                                      *
 ********************************************************/
Gaussian2D_pState BCs_y(const Gaussian2D_pState &Wb,
                        const Gaussian2D_pState &Wi,
		        const Gaussian2D_pState &dWdy,
		        const double &dy,
		        const int BC_type,
		        const int End_type) {

    switch(BC_type) {

      case BC_FIXED :
        return (Wb);
        break;

      case BC_CONSTANT_EXTRAPOLATION :
        return (Wi);
        break;

      case BC_LINEAR_EXTRAPOLATION :
        switch(End_type) {
          case LEFT_END_BOUNDARY :
            return (Wi-dWdy*dy);
            break;
          case RIGHT_END_BOUNDARY :
            return (Wi+dWdy*dy);
            break;
          default:
            return (Wi+dWdy*dy);
            break;
        } /* endswitch */
        break;

      case BC_REFLECTION :
        return (Gaussian2D_pState(Wi.d, Wi.v.x, -Wi.v.y, Wi.p.xx, -Wi.p.xy, Wi.p.yy, Wi.p.zz, Wi.erot));
        break;

      case BC_PERIODIC :
        return (Wb);
        break;

      case BC_CHARACTERISTIC :
        switch(End_type) {
          case LEFT_END_BOUNDARY :
            return (BC_Characteristic_Pressure(Wi,
                                               Wb,
	      	                               -Vector2D_NY));
            break;
          case RIGHT_END_BOUNDARY :
            return (BC_Characteristic_Pressure(Wi,
                                               Wb,
	      	                               Vector2D_NY));
            break;
          default:
            return (BC_Characteristic_Pressure(Wi,
                                               Wb,
	      	                               Vector2D_NY));
            break;
        } /* endswitch */
        break;

      default:
        return (Wb);
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
Gaussian2D_pState WaveSpeedPos(const Gaussian2D_pState &lambdas_a,
                               const Gaussian2D_pState &lambdas_l,
                               const Gaussian2D_pState &lambdas_r) {
  return (Gaussian2D_pState(HALF*(lambdas_a[1]+fabs(lambdas_a[1])),
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
Gaussian2D_pState WaveSpeedNeg(const Gaussian2D_pState &lambdas_a,
                               const Gaussian2D_pState &lambdas_l,
                               const Gaussian2D_pState &lambdas_r) {
  return (Gaussian2D_pState(HALF*(lambdas_a[1]-fabs(lambdas_a[1])),
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
Gaussian2D_pState WaveSpeedAbs(const Gaussian2D_pState &lambdas_a,
                               const Gaussian2D_pState &lambdas_l,
                               const Gaussian2D_pState &lambdas_r) {
  return (Gaussian2D_pState(fabs(lambdas_a[1]),
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
Gaussian2D_pState HartenFixPos(const Gaussian2D_pState &lambdas_a,
                               const Gaussian2D_pState &lambdas_l,
                               const Gaussian2D_pState &lambdas_r) {
  return (Gaussian2D_pState(HartenFixPos(lambdas_a[1],
                                      lambdas_l[1],
                                      lambdas_r[1]),
                         HALF*(lambdas_a[2]+fabs(lambdas_a[2])),
                         HALF*(lambdas_a[3]+fabs(lambdas_a[3])),
                         HALF*(lambdas_a[4]+fabs(lambdas_a[4])),
                         HALF*(lambdas_a[5]+fabs(lambdas_a[5])),
                         HALF*(lambdas_a[6]+fabs(lambdas_a[6])),
                         HALF*(lambdas_a[7]+fabs(lambdas_a[7])),
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
Gaussian2D_pState HartenFixNeg(const Gaussian2D_pState &lambdas_a,
                               const Gaussian2D_pState &lambdas_l,
                               const Gaussian2D_pState &lambdas_r) {
  return (Gaussian2D_pState(HartenFixNeg(lambdas_a[1],
                                      lambdas_l[1],
                                      lambdas_r[1]),
                         HALF*(lambdas_a[2]-fabs(lambdas_a[2])),
                         HALF*(lambdas_a[3]-fabs(lambdas_a[3])),
                         HALF*(lambdas_a[4]-fabs(lambdas_a[4])),
                         HALF*(lambdas_a[5]-fabs(lambdas_a[5])),
                         HALF*(lambdas_a[6]-fabs(lambdas_a[6])),
                         HALF*(lambdas_a[7]-fabs(lambdas_a[7])),
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
Gaussian2D_pState HartenFixAbs(const Gaussian2D_pState &lambdas_a,
                               const Gaussian2D_pState &lambdas_l,
                               const Gaussian2D_pState &lambdas_r) {
  return (Gaussian2D_pState(HartenFixAbs(lambdas_a[1],
                                      lambdas_l[1],
                                      lambdas_r[1]),
                         fabs(lambdas_a[2]),
                         fabs(lambdas_a[3]),
                         fabs(lambdas_a[4]),
                         fabs(lambdas_a[5]),
                         fabs(lambdas_a[6]),
                         fabs(lambdas_a[7]),
                         HartenFixAbs(lambdas_a[8],
                                      lambdas_l[8],
                                      lambdas_r[8])));
}

/*********************************************************
 * Routine: FluxRoe (Roe's flux function, x-direction)   *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for the x-direction given left and right         *
 * solution states by using the "linearized" approximate *
 * Riemann solver of Roe for the two states.  See Roe    *
 * (1981).                                               *
 *                                                       *
 *********************************************************/
Gaussian2D_cState FluxRoe(const Gaussian2D_pState &Wl,
	      	          const Gaussian2D_pState &Wr) {

    int i;
    Gaussian2D_pState Wa, dWrl, 
                   wavespeeds, lambdas_l, lambdas_r, lambdas_a;
    Gaussian2D_cState Flux;

    /* Evaluate the Roe-average primitive solution state. */
    
    Wa = RoeAverage(Wl, Wr);

    /* Evaluate the jumps in the primitive solution states. */

    dWrl = Wr-Wl;

    /* Evaluate the left, right, and average state eigenvalues. */

    lambdas_l = Wl.lambda();
    lambdas_r = Wr.lambda();
    lambdas_a = Wa.lambda();

    /* Determine the intermediate state flux. */

    if (Wa.v.x >= ZERO) {
        Flux = Wl.F();
        wavespeeds = HartenFixNeg(lambdas_a,
                                  lambdas_l,
                                  lambdas_r);
        for ( i = 1 ; i <= NUM_VAR_GAUSSIAN2D ; ++i ) {
           if (wavespeeds[i] < ZERO)
              Flux += wavespeeds[i]*(Wa.lp(i)*dWrl)*Wa.rc(i);
        } /* endfor */
    } else {
        Flux = Wr.F();
        wavespeeds = HartenFixPos(lambdas_a,
                                  lambdas_l,
                                  lambdas_r);
        for ( i = 1 ; i <= NUM_VAR_GAUSSIAN2D ; ++i ) {
           if (wavespeeds[i] > ZERO)
              Flux -= wavespeeds[i]*(Wa.lp(i)*dWrl)*Wa.rc(i);
        } /* endfor */
    } /* endif */

    /* Return solution flux. */
    
    return (Flux);
    
}

Gaussian2D_cState FluxRoe(const Gaussian2D_cState &Ul,
	      	          const Gaussian2D_cState &Ur) {
   return (FluxRoe(Ul.W(), Ur.W()));
}

/*********************************************************
 * Routine: FluxRoe_x (Roe's flux function, x-direction) *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for the x-direction given left and right         *
 * solution states by using the "linearized" approximate *
 * Riemann solver of Roe for the two states.  See Roe    *
 * (1981).                                               *
 *                                                       *
 *********************************************************/
Gaussian2D_cState FluxRoe_x(const Gaussian2D_pState &Wl,
	      	            const Gaussian2D_pState &Wr) {

    int i;
    Gaussian2D_pState Wa, dWrl, 
                   wavespeeds, lambdas_l, lambdas_r, lambdas_a;
    Gaussian2D_cState Flux;

    /* Evaluate the Roe-average primitive solution state. */
    
    Wa = RoeAverage(Wl, Wr);

    /* Evaluate the jumps in the primitive solution states. */

    dWrl = Wr-Wl;

    /* Evaluate the left, right, and average state eigenvalues. */

    lambdas_l = Wl.lambda_x();
    lambdas_r = Wr.lambda_x();
    lambdas_a = Wa.lambda_x();

    /* Determine the intermediate state flux. */

    if (Wa.v.x >= ZERO) {
        Flux = Wl.Fx();
        wavespeeds = HartenFixNeg(lambdas_a,
                                  lambdas_l,
                                  lambdas_r);
        for ( i = 1 ; i <= NUM_VAR_GAUSSIAN2D ; ++i ) {
           if (wavespeeds[i] < ZERO)
              Flux += wavespeeds[i]*(Wa.lp_x(i)*dWrl)*Wa.rc_x(i);
        } /* endfor */
    } else {
        Flux = Wr.Fx();
        wavespeeds = HartenFixPos(lambdas_a,
                                  lambdas_l,
                                  lambdas_r);
        for ( i = 1 ; i <= NUM_VAR_GAUSSIAN2D ; ++i ) {
           if (wavespeeds[i] > ZERO)
              Flux -= wavespeeds[i]*(Wa.lp_x(i)*dWrl)*Wa.rc_x(i);
        } /* endfor */
    } /* endif */

    /* Return solution flux. */
    
    return (Flux);
    
}

Gaussian2D_cState FluxRoe_x(const Gaussian2D_cState &Ul,
	      	            const Gaussian2D_cState &Ur) {
   return (FluxRoe_x(Ul.W(), Ur.W()));
}

/*********************************************************
 * Routine: FluxRoe_y (Roe's flux function, y-direction) *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for the y-direction given left and right         *
 * solution states by using the "linearized" approximate *
 * Riemann solver of Roe for the two states.  See Roe    *
 * (1981).                                               *
 *                                                       *
 *********************************************************/
Gaussian2D_cState FluxRoe_y(const Gaussian2D_pState &Wl,
	      	            const Gaussian2D_pState &Wr) {

    int i;
    Gaussian2D_pState Wa, dWrl, 
                   wavespeeds, lambdas_l, lambdas_r, lambdas_a;
    Gaussian2D_cState Flux;

    /* Evaluate the Roe-average primitive solution state. */
    
    Wa = RoeAverage(Wl, Wr);

    /* Evaluate the jumps in the primitive solution states. */

    dWrl = Wr-Wl;

    /* Evaluate the left, right, and average state eigenvalues. */

    lambdas_l = Wl.lambda_y();
    lambdas_r = Wr.lambda_y();
    lambdas_a = Wa.lambda_y();

    /* Determine the intermediate state flux. */

    if (Wa.v.y >= ZERO) {
        Flux = Wl.Fy();
        wavespeeds = HartenFixNeg(lambdas_a,
                                  lambdas_l,
                                  lambdas_r);
        for ( i = 1 ; i <= NUM_VAR_GAUSSIAN2D ; ++i ) {
           if (wavespeeds[i] < ZERO)
              Flux += wavespeeds[i]*(Wa.lp_y(i)*dWrl)*Wa.rc_y(i);
        } /* endfor */
    } else {
        Flux = Wr.Fy();
        wavespeeds = HartenFixPos(lambdas_a,
                                  lambdas_l,
                                  lambdas_r);
        for ( i = 1 ; i <= NUM_VAR_GAUSSIAN2D ; ++i ) {
           if (wavespeeds[i] > ZERO)
              Flux -= wavespeeds[i]*(Wa.lp_y(i)*dWrl)*Wa.rc_y(i);
        } /* endfor */
    } /* endif */

    /* Return solution flux. */
    
    return (Flux);
    
}

Gaussian2D_cState FluxRoe_y(const Gaussian2D_cState &Ul,
	       	            const Gaussian2D_cState &Ur) {
   return (FluxRoe_y(Ul.W(), Ur.W()));
}  

/*********************************************************
 * Routine: FluxRoe_n (Roe's flux function, n-direction) *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for an arbitrary direction defined by a unit     *
 * normal vector in the direction of interest, given     *
 * left and right solution states.  The problem is       *
 * solved by first applying a frame rotation to rotate   *
 * the problem to a local frame aligned with the unit    *
 * normal vector and then by using the "linearized"      *
 * approximate Riemann solver of Roe to specify the flux *
 * in terms of the rotated solution states.  See Roe     *
 * (1981).                                               *
 *                                                       *
 *********************************************************/
Gaussian2D_cState FluxRoe_n(const Gaussian2D_pState &Wl,
	      	            const Gaussian2D_pState &Wr,
                            const Vector2D &norm_dir) {

    double cos_angle, sin_angle;
    Gaussian2D_pState Wl_rotated, Wr_rotated;
    Gaussian2D_cState Flux, Flux_rotated;

    // Determine the direction cosine's for the frame
    // rotation. 

    cos_angle = norm_dir.x; 
    sin_angle = norm_dir.y;

    // Apply the frame rotation and evaluate left and right
    // solution states in the local rotated frame defined
    // by the unit normal vector.

    Wl_rotated.d = Wl.d;
    Wl_rotated.v.x = Wl.v.x*cos_angle +
                     Wl.v.y*sin_angle;
    Wl_rotated.v.y = - Wl.v.x*sin_angle +
                       Wl.v.y*cos_angle;
    Wl_rotated.p.xx = Wl.p.xx*cos_angle*cos_angle+Wl.p.yy*sin_angle*sin_angle
                     +2.0*Wl.p.xy*cos_angle*sin_angle;
    Wl_rotated.p.xy = Wl.p.xy*(cos_angle*cos_angle-sin_angle*sin_angle)
                     -(Wl.p.xx-Wl.p.yy)*cos_angle*sin_angle;
    Wl_rotated.p.yy = Wl.p.xx*sin_angle*sin_angle+Wl.p.yy*cos_angle*cos_angle
                     -2.0*Wl.p.xy*cos_angle*sin_angle;
    Wl_rotated.p.zz = Wl.p.zz;
    Wl_rotated.erot = Wl.erot;


    Wr_rotated.d = Wr.d;
    Wr_rotated.v.x = Wr.v.x*cos_angle +
                     Wr.v.y*sin_angle;
    Wr_rotated.v.y = - Wr.v.x*sin_angle +
                       Wr.v.y*cos_angle;
    Wr_rotated.p.xx = Wr.p.xx*cos_angle*cos_angle+Wr.p.yy*sin_angle*sin_angle
                     +2.0*Wr.p.xy*cos_angle*sin_angle;
    Wr_rotated.p.xy = Wr.p.xy*(cos_angle*cos_angle-sin_angle*sin_angle)
                     -(Wr.p.xx-Wr.p.yy)*cos_angle*sin_angle;
    Wr_rotated.p.yy = Wr.p.xx*sin_angle*sin_angle+Wr.p.yy*cos_angle*cos_angle
                     -2.0*Wr.p.xy*cos_angle*sin_angle;
    Wr_rotated.p.zz = Wr.p.zz;
    Wr_rotated.erot = Wr.erot;


    // Evaluate the intermediate state solution 
    // flux in the rotated frame.

    Flux_rotated = FluxRoe_x(Wl_rotated, Wr_rotated);

    // Rotate back to the original Cartesian reference
    // frame and return the solution flux.

    sin_angle = -sin_angle;

    Flux.d = Flux_rotated.d;
    Flux.dv.x = Flux_rotated.dv.x*cos_angle +
               Flux_rotated.dv.y*sin_angle;
    Flux.dv.y = - Flux_rotated.dv.x*sin_angle +
                 Flux_rotated.dv.y*cos_angle;
    Flux.E.xx = Flux_rotated.E.xx*cos_angle*cos_angle+Flux_rotated.E.yy*sin_angle*sin_angle
                +2.0*Flux_rotated.E.xy*cos_angle*sin_angle;
    Flux.E.xy = Flux_rotated.E.xy*(cos_angle*cos_angle-sin_angle*sin_angle)
                -(Flux_rotated.E.xx-Flux_rotated.E.yy)*cos_angle*sin_angle;
    Flux.E.yy = Flux_rotated.E.xx*sin_angle*sin_angle+Flux_rotated.E.yy*cos_angle*cos_angle
                -2.0*Flux_rotated.E.xy*cos_angle*sin_angle;
    Flux.E.zz = Flux_rotated.E.zz;
    Flux.erot = Flux_rotated.erot;

    return (Flux);

}

Gaussian2D_cState FluxRoe_n(const Gaussian2D_cState &Ul,
	      	            const Gaussian2D_cState &Ur,
                            const Vector2D &norm_dir) {
    return (FluxRoe_n(Ul.W(), Ur.W(), norm_dir));
}

/**********************************************************************
 * Routine: FluxRoe_MB (Roe's flux function)                          *
 *                                                                    *
 * This function returns the intermediate state solution flux for the *
 * given left and right solution states by using the "linearized"     *
 * approximate Riemann solver of Roe for the two states.  See Roe     *
 * (1981).                                                            *
 *                                                                    *
 **********************************************************************/
Gaussian2D_cState FluxRoe_MB(const Gaussian2D_pState &Wl,
			  const Gaussian2D_pState &Wr,
			  const Vector2D &V) {

  Gaussian2D_pState Wa, dWrl, wavespeeds, lambdas_l, lambdas_r, lambdas_a;
  Gaussian2D_cState Flux;
  
  // Evaluate the Roe-average primitive solution state.
  Wa = RoeAverage(Wl,Wr);

  // Evaluate the jumps in the primitive solution states.
  dWrl = Wr - Wl;

  // Evaluate the left, right, and average state eigenvalues.
  lambdas_l = Wl.lambda_x(V);
  lambdas_r = Wr.lambda_x(V);
  lambdas_a = Wa.lambda_x(V);

  // Determine the intermediate state flux.
  if (Wa.v.x >= ZERO) {
    Flux = Wl.F(V);
    wavespeeds = HartenFixNeg(lambdas_a,lambdas_l,lambdas_r);
    for (int i = 1; i <= NUM_VAR_GAUSSIAN2D; i++) {
      if (wavespeeds[i] < ZERO)
	Flux += wavespeeds[i]*(Wa.lp_x(i)*dWrl)*Wa.rc_x(i);
    }
  } else {
    Flux = Wr.F(V);
    wavespeeds = HartenFixPos(lambdas_a,lambdas_l,lambdas_r);
    for (int i = 1; i <= NUM_VAR_GAUSSIAN2D; i++) {
      if (wavespeeds[i] > ZERO)
	Flux -= wavespeeds[i]*(Wa.lp_x(i)*dWrl)*Wa.rc_x(i);
    }
  }

  // Return the solution flux.
  return Flux;
    
}

Gaussian2D_cState FluxRoe_MB(const Gaussian2D_cState &Ul,
			  const Gaussian2D_cState &Ur,
			  const Vector2D &V) {
  return FluxRoe_MB(Ul.W(),Ur.W(),V);
}

/**********************************************************************
 * Routine: FluxRoe_MB_n (Roe's flux function, n-direction)           *
 *                                                                    *
 * This function returns the intermediate state solution flux for an  *
 * arbitrary direction defined by a unit normal vector in the         *
 * direction of interest, given left and right solution states.  The  *
 * problem is solved by first applying a frame rotation to rotate the *
 * problem to a local frame aligned with the unit normal vector and   *
 * then by using the "linearized" approximate Riemann solver of Roe   *
 * to specify the flux in terms of the rotated solution states.  See  *
 * Roe (1981).                                                        *
 *                                                                    *
 **********************************************************************/
Gaussian2D_cState FluxRoe_MB_n(const Gaussian2D_pState &Wl,
			       const Gaussian2D_pState &Wr,
			       const Vector2D &V,
			       const Vector2D &norm_dir) {

  Gaussian2D_pState Wl_rotated, Wr_rotated;
  Gaussian2D_cState Flux, Flux_rotated;
  Vector2D V_rotated;

  // Apply the frame rotation and evaluate left and right solution 
  // states in the local rotated frame defined by the unit normal 
  // vector.
  Wl_rotated = Rotate(Wl,norm_dir);
  Wr_rotated = Rotate(Wr,norm_dir);
  V_rotated.Rotate(V,norm_dir);

  // Evaluate the intermediate state solution flux in the rotated frame.
  Flux_rotated = FluxRoe_MB(Wl_rotated,Wr_rotated,V_rotated);

  // Rotate back to the original Cartesian reference frame and return 
  // the solution flux.
  Flux = Rotate(Flux_rotated,Vector2D(norm_dir.x,-norm_dir.y));

  // Return the solution flux.
  return Flux;

}

Gaussian2D_cState FluxRoe_MB_n(const Gaussian2D_cState &Ul,
			       const Gaussian2D_cState &Ur,
			       const Vector2D &V,
			       const Vector2D &norm_dir) {
  return FluxRoe_MB_n(Ul.W(),Ur.W(),V,norm_dir);
}

/*********************************************************
 * Routine: FluxHLLE (Harten-Lax-van Leer flux function, *
 *                    x-direction)                       *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for the x-direction given left and right         *
 * solution states by using the so-called Harten-Lax-    *
 * van Leer approximation for the fluxes.  See Harten,   *
 * Lax, van Leer (1983).                                 *
 *                                                       *
 *********************************************************/
Gaussian2D_cState FluxHLLE(const Gaussian2D_pState &Wl,
			   const Gaussian2D_pState &Wr) {

    double wavespeed_l, wavespeed_r;
    Gaussian2D_pState Wa, lambdas_l, lambdas_r, lambdas_a;
    Gaussian2D_cState Flux, dUrl;
    
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
    wavespeed_r = max(lambdas_r[NUM_VAR_GAUSSIAN2D],
                      lambdas_a[NUM_VAR_GAUSSIAN2D]);
 
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

Gaussian2D_cState FluxHLLE(const Gaussian2D_cState &Ul,
			   const Gaussian2D_cState &Ur) {
   return (FluxHLLE(Ul.W(), Ur.W()));
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
Gaussian2D_cState FluxHLLE_x(const Gaussian2D_pState &Wl,
			     const Gaussian2D_pState &Wr) {

    double wavespeed_l, wavespeed_r;
    Gaussian2D_pState Wa, lambdas_l, lambdas_r, lambdas_a;
    Gaussian2D_cState Flux, dUrl;
    
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
    wavespeed_r = max(lambdas_r[NUM_VAR_GAUSSIAN2D],
                      lambdas_a[NUM_VAR_GAUSSIAN2D]);

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

Gaussian2D_cState FluxHLLE_x(const Gaussian2D_cState &Ul,
			     const Gaussian2D_cState &Ur) {
   return (FluxHLLE_x(Ul.W(), Ur.W()));
}

/*********************************************************
 * Routine: FluxHLLE_y (Harten-Lax-van Leer flux         *
 *                      function, y-direction)           *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for the y-direction given left and right         *
 * solution states by using the so-called Harten-Lax-    *
 * van Leer approximation for the fluxes.  See Harten,   *
 * Lax, van Leer (1983).                                 *
 *                                                       *
 *********************************************************/
Gaussian2D_cState FluxHLLE_y(const Gaussian2D_pState &Wl,
			     const Gaussian2D_pState &Wr) {

    double wavespeed_l, wavespeed_r;
    Gaussian2D_pState Wa, lambdas_l, lambdas_r, lambdas_a;
    Gaussian2D_cState Flux, dUrl;
    
    /* Evaluate the Roe-average primitive solution state. */
    
    Wa = RoeAverage(Wl, Wr);

   /* Evaluate the jumps in the conserved solution states. */

    dUrl = Wr.U()-Wl.U();

    /* Evaluate the left, right, and average state eigenvalues. */

    lambdas_l = Wl.lambda_y();
    lambdas_r = Wr.lambda_y();
    lambdas_a = Wa.lambda_y();

    /* Determine the intermediate state flux. */

    wavespeed_l = min(lambdas_l[1],
                      lambdas_a[1]);
    wavespeed_r = max(lambdas_r[NUM_VAR_GAUSSIAN2D],
                      lambdas_a[NUM_VAR_GAUSSIAN2D]);
 
    if (wavespeed_l >= ZERO) {
        Flux = Wl.F();
    } else if (wavespeed_r <= ZERO) {
        Flux = Wr.F();
    } else {
        Flux =   ( (wavespeed_r*Wl.Fy()-wavespeed_l*Wr.Fy())
                  +(wavespeed_l*wavespeed_r)*dUrl)/
                 (wavespeed_r-wavespeed_l);
    } /* endif */

    /* Return solution flux. */

    return (Flux);

}

Gaussian2D_cState FluxHLLE_y(const Gaussian2D_cState &Ul,
			     const Gaussian2D_cState &Ur) {
   return (FluxHLLE_y(Ul.W(), Ur.W()));
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
Gaussian2D_cState FluxHLLE_n(const Gaussian2D_pState &Wl,
			     const Gaussian2D_pState &Wr,
			     const Vector2D &norm_dir) {

    double cos_angle, sin_angle;
    Gaussian2D_pState Wl_rotated, Wr_rotated;
    Gaussian2D_cState Flux, Flux_rotated;

    /* Determine the direction cosine's for the frame
       rotation. */

    cos_angle = norm_dir.x; 
    sin_angle = norm_dir.y;

    /* Apply the frame rotation and evaluate left and right
       solution states in the local rotated frame defined
       by the unit normal vector. */

    Wl_rotated.d = Wl.d;
    Wl_rotated.v.x = Wl.v.x*cos_angle +
                     Wl.v.y*sin_angle;
    Wl_rotated.v.y = - Wl.v.x*sin_angle +
                       Wl.v.y*cos_angle;
    Wl_rotated.p.xx = Wl.p.xx*cos_angle*cos_angle+Wl.p.yy*sin_angle*sin_angle
                     +2.0*Wl.p.xy*cos_angle*sin_angle;
    Wl_rotated.p.xy = Wl.p.xy*(cos_angle*cos_angle-sin_angle*sin_angle)
                     -(Wl.p.xx-Wl.p.yy)*cos_angle*sin_angle;
    Wl_rotated.p.yy = Wl.p.xx*sin_angle*sin_angle+Wl.p.yy*cos_angle*cos_angle
                     -2.0*Wl.p.xy*cos_angle*sin_angle;
    Wl_rotated.p.zz = Wl.p.zz;
    Wl_rotated.erot = Wl.erot;


    Wr_rotated.d = Wr.d;
    Wr_rotated.v.x = Wr.v.x*cos_angle +
                     Wr.v.y*sin_angle;
    Wr_rotated.v.y = - Wr.v.x*sin_angle +
                       Wr.v.y*cos_angle;
    Wr_rotated.p.xx = Wr.p.xx*cos_angle*cos_angle+Wr.p.yy*sin_angle*sin_angle
                     +2.0*Wr.p.xy*cos_angle*sin_angle;
    Wr_rotated.p.xy = Wr.p.xy*(cos_angle*cos_angle-sin_angle*sin_angle)
                     -(Wr.p.xx-Wr.p.yy)*cos_angle*sin_angle;
    Wr_rotated.p.yy = Wr.p.xx*sin_angle*sin_angle+Wr.p.yy*cos_angle*cos_angle
                     -2.0*Wr.p.xy*cos_angle*sin_angle;
    Wr_rotated.p.zz = Wr.p.zz;
    Wr_rotated.erot = Wr.erot;

    /* Evaluate the intermediate state solution 
       flux in the rotated frame. */

    Flux_rotated = FluxHLLE_x(Wl_rotated, Wr_rotated);

    /* Rotate back to the original Cartesian reference
       frame and return the solution flux. */

    sin_angle = -sin_angle;

    Flux.d = Flux_rotated.d;
    Flux.dv.x = Flux_rotated.dv.x*cos_angle +
               Flux_rotated.dv.y*sin_angle;
    Flux.dv.y = - Flux_rotated.dv.x*sin_angle +
                 Flux_rotated.dv.y*cos_angle;
    Flux.E.xx = Flux_rotated.E.xx*cos_angle*cos_angle+Flux_rotated.E.yy*sin_angle*sin_angle
                +2.0*Flux_rotated.E.xy*cos_angle*sin_angle;
    Flux.E.xy = Flux_rotated.E.xy*(cos_angle*cos_angle-sin_angle*sin_angle)
                -(Flux_rotated.E.xx-Flux_rotated.E.yy)*cos_angle*sin_angle;
    Flux.E.yy = Flux_rotated.E.xx*sin_angle*sin_angle+Flux_rotated.E.yy*cos_angle*cos_angle
                -2.0*Flux_rotated.E.xy*cos_angle*sin_angle;
    Flux.E.zz = Flux_rotated.E.zz;
    Flux.erot = Flux_rotated.erot;

    return (Flux);

}

Gaussian2D_cState FluxHLLE_n(const Gaussian2D_cState &Ul,
			     const Gaussian2D_cState &Ur,
			     const Vector2D &norm_dir) {
    return (FluxHLLE_n(Ul.W(), Ur.W(), norm_dir));
}

/*********************************************************
 * Routine: HLLE_wavespeeds                              *
 *                                                       *
 * This function returns lambda plus and lambda minus    *
 * for rotated Riemann problem aligned with norm_dir     *
 * given unroated solution states Wl and Wr.             *
 * Note: wavespeed.x = wavespeed_l = lambda minus.       *
 *       wavespeed.y = wavespeed_r = lambda plus.        *
 *                                                       *
 *********************************************************/
Vector2D HLLE_wavespeeds(const Gaussian2D_pState &Wl,
                         const Gaussian2D_pState &Wr,
                         const Vector2D &norm_dir) {

    Vector2D wavespeed;
    Gaussian2D_pState Wa_n, lambdas_l, lambdas_r, lambdas_a, Wl_n, Wr_n;  

    /* Use rotated values to calculate eignvalues */
    Wl_n = Rotate(Wl, norm_dir);
    Wr_n = Rotate(Wr, norm_dir);

    /* Evaluate the Roe-average primitive solution state. */
    Wa_n = RoeAverage(Wl_n, Wr_n);
    
    /* Evaluate the left, right, and average state eigenvalues. */

    lambdas_l = Wl_n.lambda_x();
    lambdas_r = Wr_n.lambda_x();
    lambdas_a = Wa_n.lambda_x();

    /* Determine the intermediate state flux. */

    wavespeed.x = min(lambdas_l[1],
                      lambdas_a[1]);
    wavespeed.y = max(lambdas_r[NUM_VAR_GAUSSIAN2D],
                      lambdas_a[NUM_VAR_GAUSSIAN2D]);
 
    wavespeed.x = min(wavespeed.x, ZERO); //lambda minus
    wavespeed.y = max(wavespeed.y, ZERO); //lambda plus 

    return (wavespeed);

}


/*********************************************************
 * Routine: Imposed_adiabatic_wall_n                     *
 *                                                       *
 * This function returns the flux at an adiabatic wall   *
 * without using ghost cells.                            *
 *                                                       *
 *********************************************************/
Gaussian2D_cState Imposed_adiabatic_wall_n(const Gaussian2D_pState &Wr,
					   const Gaussian2D_pState &Wo,
					   const Vector2D &norm_dir) {
    double cos_angle, sin_angle;
    Gaussian2D_pState Wr_rotated, Knudsen_Layer;
    Gaussian2D_cState Flux, Flux_rotated;
    double uWallr, vWallr;

    /* Determine the direction cosine's for the frame
       rotation. */

    cos_angle = norm_dir.x; 
    sin_angle = norm_dir.y;

    /* Apply the frame rotation and evaluate left and right
       solution states in the local rotated frame defined
       by the unit normal vector. */

    Wr_rotated.d = Wr.d;
    Wr_rotated.v.x = Wr.v.x*cos_angle +
                     Wr.v.y*sin_angle;
    Wr_rotated.v.y = - Wr.v.x*sin_angle +
                       Wr.v.y*cos_angle;
    Wr_rotated.p.xx = Wr.p.xx*cos_angle*cos_angle+Wr.p.yy*sin_angle*sin_angle
                     +2.0*Wr.p.xy*cos_angle*sin_angle;
    Wr_rotated.p.xy = Wr.p.xy*(cos_angle*cos_angle-sin_angle*sin_angle)
                     -(Wr.p.xx-Wr.p.yy)*cos_angle*sin_angle;
    Wr_rotated.p.yy = Wr.p.xx*sin_angle*sin_angle+Wr.p.yy*cos_angle*cos_angle
                     -2.0*Wr.p.xy*cos_angle*sin_angle;
    Wr_rotated.p.zz = Wr.p.zz;
    Wr_rotated.erot = Wr.erot;

    uWallr =  Wo.v.x*cos_angle + Wo.v.y*sin_angle;
    vWallr =  - Wo.v.x*sin_angle + Wo.v.y*cos_angle;

    /* Find State in Knudsen Layer */

    Knudsen_Layer.d = Wr_rotated.d;

    Knudsen_Layer.v.x = 0.0;

    Knudsen_Layer.v.y = (2.0-Wr_rotated.alpha_m)*(Wr_rotated.v.y/2.0+
						Wr_rotated.p.xy/sqrt(2.0*PI*Wr_rotated.p.xx))
                                                +Wr_rotated.alpha_m/2.0*uWallr;
    Knudsen_Layer.p.xx = Wr_rotated.p.xx;

    Knudsen_Layer.p.xy = Wr_rotated.alpha_m*(Wr_rotated.p.xy/2.0-sqrt(Wr_rotated.d*Wr_rotated.p.xx/(2.0*PI))*
					  (Wr_rotated.v.y-vWallr));

    Knudsen_Layer.p.yy = Wr_rotated.p.yy+(2.0-Wr_rotated.alpha_m)*(Wr_rotated.d/2.0*sqr(Wr_rotated.v.y-Knudsen_Layer.v.y)
							       +Wr_rotated.p.xy*sqrt(2.0*Wr_rotated.d/(PI*Wr_rotated.p.xx))*(Wr_rotated.v.y-Knudsen_Layer.v.y)
								 +Wr_rotated.alpha_m*Wr_rotated.d/2.0*sqr(vWallr-Wr_rotated.v.y));

    Knudsen_Layer.p.xx = Wr_rotated.p.zz;

    Knudsen_Layer.erot = Wr_rotated.erot;

    //try for now
    /*
    Knudsen_Layer.v.y = Wr_rotated.v.y;

    Knudsen_Layer.p.yy = Wr_rotated.p.yy;
    */
    /* Evaluate the intermediate state solution 
       flux in the rotated frame. */

    Flux_rotated = F(Knudsen_Layer);

    /* Rotate back to the original Cartesian reference
       frame and return the solution flux. */

    sin_angle = -sin_angle;

    Flux.d = Flux_rotated.d;
    Flux.dv.x = Flux_rotated.dv.x*cos_angle +
               Flux_rotated.dv.y*sin_angle;
    Flux.dv.y = - Flux_rotated.dv.x*sin_angle +
                 Flux_rotated.dv.y*cos_angle;
    Flux.E.xx = Flux_rotated.E.xx*cos_angle*cos_angle+Flux_rotated.E.yy*sin_angle*sin_angle
                +2.0*Flux_rotated.E.xy*cos_angle*sin_angle;
    Flux.E.xy = Flux_rotated.E.xy*(cos_angle*cos_angle-sin_angle*sin_angle)
                -(Flux_rotated.E.xx-Flux_rotated.E.yy)*cos_angle*sin_angle;
    Flux.E.yy = Flux_rotated.E.xx*sin_angle*sin_angle+Flux_rotated.E.yy*cos_angle*cos_angle
                -2.0*Flux_rotated.E.xy*cos_angle*sin_angle;
    Flux.E.zz = Flux_rotated.E.zz;
    Flux.erot = Flux_rotated.erot;
    /*
    cout << Wr << endl;
    cout << Wr_rotated << endl;
    cout << Flux << endl;
    */
    return (Flux);

}


/*********************************************************
 * Routine: FluxKinetic_x                                *
 *                                                       *
 * This function uses an integration of two half         *
 * Gaussians to determine the intercellular flux         *
 * in the x-direction.                                   *
 *                                                       *
 *********************************************************/
Gaussian2D_cState FluxKinetic_x(const Gaussian2D_pState &W1,
				const Gaussian2D_pState &W2) {

  double expl, expr, erfl, erfr;
  double sqrtl1, sqrtr1,sqrtl2, sqrtr2;
  Gaussian2D_cState Flux;

  expl = exp(-W1.d*sqr(W1.v.x)/(2.0*W1.p.xx));
  expr = exp(-W2.d*sqr(W2.v.x)/(2.0*W2.p.xx));
  erfl = erf(sqrt(W1.d/(2.0*W1.p.xx))*W1.v.x);
  erfr = erf(sqrt(W2.d/(2.0*W2.p.xx))*W2.v.x);
  sqrtl1 = sqrt(W1.d*W1.p.xx/(2.0*PI));
  sqrtr1 = sqrt(W2.d*W2.p.xx/(2.0*PI));
  sqrtl2 = sqrt(W1.p.xx/(2.0*PI*W1.d));
  sqrtr2 = sqrt(W2.p.xx/(2.0*PI*W2.d));

  Flux[1] = 0.5*(W1.d*W1.v.x)*(1.0+erfl) + sqrtl1*expl +
            0.5*(W2.d*W2.v.x)*(1.0-erfr) - sqrtr1*expr;

  Flux[2] = 0.5*(W1.d*sqr(W1.v.x)+W1.p.xx)*(1.0+erfl) + sqrtl1*W1.v.x*expl +
            0.5*(W2.d*sqr(W2.v.x)+W2.p.xx)*(1.0-erfr) - sqrtr1*W1.v.x*expr;

  Flux[3] = 0.5*(W1.d*W1.v.x*W1.v.y+W1.p.xy)*(1.0+erfl) + sqrtl1*W1.v.y*expl +
            0.5*(W2.d*W2.v.x*W2.v.y+W2.p.xy)*(1.0-erfr) + sqrtr1*W2.v.y*expr;

  Flux[4] = 0.5*(W1.v.x*(W1.d*sqr(W1.v.x)+3.0*W1.p.xx))*(1.0+erfl)
                +(sqrtl1*sqr(W1.v.x)+2.0*sqrtl2*W1.p.xx)*expl +
            0.5*(W2.v.x*(W2.d*sqr(W2.v.x)+3.0*W2.p.xx))*(1.0-erfr)
		-(sqrtr1*sqr(W2.v.x)+2.0*sqrtr2*W2.p.xx)*expr;

  Flux[5] = 0.5*(W1.d*sqr(W1.v.x)*W1.v.y+2.0*W1.v.x*W1.p.xy+W1.v.y*W1.p.xx)*(1.0+erfl)
                +(sqrtl1*W1.v.x*W1.v.y+2.0*sqrtl2*W1.p.xy)*expl +
            0.5*(W2.d*sqr(W2.v.x)*W2.v.y+2.0*W2.v.x*W2.p.xy+W2.v.y*W2.p.xx)*(1.0-erfr)
                -(sqrtr1*W2.v.x*W2.v.y+2.0*sqrtr2*W2.p.xy)*expr;

  Flux[6] = 0.5*(W1.d*sqr(W1.v.y)*W1.v.x+2.0*W1.v.y*W1.p.xy+W1.v.x*W1.p.yy)*(1.0+erfl)
                +(sqrtl1*sqr(W1.v.y)+sqrtl2*W1.p.yy+W1.p.xy/sqrt(2.0*PI*W1.p.xx*W1.d))*expl +
            0.5*(W2.d*sqr(W2.v.y)*W2.v.x+2.0*W2.v.y*W2.p.xy+W2.v.x*W2.p.yy)*(1.0-erfr)
                -(sqrtr1*sqr(W2.v.y)+sqrtr2*W2.p.yy+W2.p.xy/sqrt(2.0*PI*W2.p.xx*W2.d))*expr;

  Flux[7] = 0.5*(W1.v.x*W1.p.zz)*(1.0+erfl) + sqrtl2*W1.p.zz*expl +
            0.5*(W2.v.x*W2.p.zz)*(1.0-erfr) - sqrtr2*W2.p.zz*expr;

  Flux[7] = 0.5*(W1.v.x*W1.erot)*(1.0+erfl) + sqrtl2*W1.erot*expl +
            0.5*(W2.v.x*W2.erot)*(1.0-erfr) - sqrtr2*W2.erot*expr;

  return(Flux);

}

#ifdef _GAUSSIAN_HEAT_TRANSFER_
/**********************************************************************
 * Routine: HeatFlux_n                                                *
 *                                                                    *
 * This function returns the intermediate state solution viscous flux *
 * given the primitive variable solution state and the gradients of   *
 * the primitive variables.                                           *
 *                                                                    *
 **********************************************************************/
Gaussian2D_cState HeatFlux_n(const Vector2D &X,
			     Gaussian2D_pState &W,
			     const Gaussian2D_pState &dWdx,
			     const Gaussian2D_pState &dWdy,
			     const Vector2D &norm_dir,
			     const int &Axisymmetric) {

  Gaussian2D_cState Gx, Gy, U;

  // Compute the intermediate state viscous stress tensor and heat flux
  // vector.
  W.ComputeHeatTerms(dWdx,dWdy,X,Axisymmetric);
  U = W.U(); U.q = W.q;

  // Compute the Cartesian components of the intermediate state
  // solution viscous flux.
  Gx = U.Gx(dWdx);
  Gy = U.Gy(dWdy);

  // Return the intermediate state solution viscous flux.
  return Gx*norm_dir.x + Gy*norm_dir.y;

}

/**********************************************************************
 * Routine: HeatFluxDiamondPath_n                                     *
 *                                                                    *
 * This routine computes the viscous flux at the specified quadrature *
 * point, X.  The gradient of the primitive variables is computed on  *
 * the diamond-path defined by the points X1, X2, X3, and X4.  Only   *
 * half of the diamond (three points) is used at solid or transpiring *
 * boundaries.                                                        *
 *                                                                    *
 **********************************************************************/
Gaussian2D_cState HeatFluxDiamondPath_n(const Vector2D &X,
					const Vector2D &Xl, const Gaussian2D_pState &Wl,
					const Vector2D &Xd, const Gaussian2D_pState &Wd,
					const Vector2D &Xr, const Gaussian2D_pState &Wr,
					const Vector2D &Xu, const Gaussian2D_pState &Wu,
					const Vector2D &norm_dir,
					const int &Axisymmetric,
					const int &stencil_flag) {

  if (stencil_flag == DIAMONDPATH_NONE) return Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO);

  Gaussian2D_pState W_face, dWdxl, dWdyl, dWdxr, dWdyr, dWdx, dWdy;
  Gaussian2D_cState Flux;
  double A, Al, Ar;
  Vector2D ndl, nud, nlu, nrd, nur, ndu;
  int error_flag;

  // Compute Green-Gauss integration on 'left' triangle:
  if (stencil_flag == DIAMONDPATH_QUADRILATERAL_RECONSTRUCTION ||
      stencil_flag == DIAMONDPATH_LEFT_TRIANGLE_HEATFLUX ||
      stencil_flag == DIAMONDPATH_LEFT_TRIANGLE_ISOTHERMAL) {
    ndl = Vector2D((Xd.y-Xl.y),-(Xd.x-Xl.x));
    nud = Vector2D((Xu.y-Xd.y),-(Xu.x-Xd.x));
    nlu = Vector2D((Xl.y-Xu.y),-(Xl.x-Xu.x));
    Al = HALF*((Xd-Xl)^(Xu-Xl));
    W_face = HALF*(Wd+Wl);
    dWdxl = W_face*ndl.x;
    dWdyl = W_face*ndl.y;
    W_face = HALF*(Wu+Wd);
    dWdxl += W_face*nud.x;
    dWdyl += W_face*nud.y;
    W_face = HALF*(Wl+Wu);
    dWdxl += W_face*nlu.x;
    dWdyl += W_face*nlu.y;
    dWdxl /= Al;
    dWdyl /= Al;
  }

  // Compute Green-Gauss integration on 'right' triangle:
  if (stencil_flag == DIAMONDPATH_QUADRILATERAL_RECONSTRUCTION ||
      stencil_flag == DIAMONDPATH_RIGHT_TRIANGLE_HEATFLUX ||
      stencil_flag == DIAMONDPATH_RIGHT_TRIANGLE_ISOTHERMAL) {
    nrd = Vector2D((Xr.y-Xd.y),-(Xr.x-Xd.x));
    nur = Vector2D((Xu.y-Xr.y),-(Xu.x-Xr.x));
    ndu = Vector2D((Xd.y-Xu.y),-(Xd.x-Xu.x));
    Ar = HALF*((Xr-Xu)^(Xr-Xd));
    W_face = HALF*(Wr+Wd);
    dWdxr = W_face*nrd.x;
    dWdyr = W_face*nrd.y;
    W_face = HALF*(Wu+Wr);
    dWdxr += W_face*nur.x;
    dWdyr += W_face*nur.y;
    W_face = HALF*(Wd+Wu);
    dWdxr += W_face*ndu.x;
    dWdyr += W_face*ndu.y;
    dWdxr /= Ar;
    dWdyr /= Ar;
  }

  // Compute the viscous flux.
  if (stencil_flag == DIAMONDPATH_QUADRILATERAL_RECONSTRUCTION) {
    // Determine the total area.
    A = Al + Ar;
    // Determine the area-averaged gradients of the primitive variables.
    dWdx = (Al*dWdxl + Ar*dWdxr)/A;
    dWdy = (Al*dWdyl + Ar*dWdyr)/A;
    // Determine the primitive state at the quadrature point.
    //W_face = (Wl + Wd + Wr + Wu)/FOUR;
    error_flag = Bilinear_Interpolation(Wl,Xl,Wu,Xu,Wr,Xr,Wd,Xd,X,W_face);
    //if (error_flag) return error_flag;
    Flux = HeatFlux_n(X,W_face,dWdx,dWdy,norm_dir,Axisymmetric);

  } else if (stencil_flag == DIAMONDPATH_LEFT_TRIANGLE_HEATFLUX ||
	     stencil_flag == DIAMONDPATH_LEFT_TRIANGLE_ISOTHERMAL) {
    W_face = Wd;  //Remember Wd = Wu
    if (stencil_flag == DIAMONDPATH_LEFT_TRIANGLE_HEATFLUX) {
      Flux = HeatFlux_n(X,W_face,dWdxl,dWdyl,norm_dir,Axisymmetric);
    } else if (stencil_flag == DIAMONDPATH_LEFT_TRIANGLE_ISOTHERMAL) {
      Flux = HeatFlux_n(X,W_face,dWdxl,dWdyl,norm_dir,Axisymmetric);
    }

  } else if (stencil_flag == DIAMONDPATH_RIGHT_TRIANGLE_HEATFLUX ||
	     stencil_flag == DIAMONDPATH_RIGHT_TRIANGLE_ISOTHERMAL) {
    W_face = Wd;   //Remember Wd = Wu
    if (stencil_flag == DIAMONDPATH_RIGHT_TRIANGLE_HEATFLUX) {
      Flux = HeatFlux_n(X,W_face,dWdxr,dWdyr,norm_dir,Axisymmetric);
    } else if (stencil_flag == DIAMONDPATH_RIGHT_TRIANGLE_ISOTHERMAL) {
      Flux = HeatFlux_n(X,W_face,dWdxr,dWdyr,norm_dir,Axisymmetric);
    }

  }

  // Return the viscous flux.
  return Flux;

}

/**********************************************************************
 * Routine: HeatFluxHybrid_n                                          *
 *                                                                    *
 * This function returns the intermediate state solution viscous flux *
 * calculated by the arithmetic mean of the cell-centred flux terms   *
 * of the neighbouring cells.                                         *
 *                                                                    *
 **********************************************************************/
Gaussian2D_cState HeatFluxHybrid_n(const Vector2D &X,
				   Gaussian2D_pState &W,
				   const Vector2D &X1,
				   const Gaussian2D_pState &W1,
				   const Gaussian2D_pState &dW1dx,
				   const Gaussian2D_pState &dW1dy,
				   const Vector2D &X2,
				   const Gaussian2D_pState &W2,
				   const Gaussian2D_pState &dW2dx,
				   const Gaussian2D_pState &dW2dy,
				   const Vector2D &norm_dir,
				   const int &Axisymmetric) {

  Gaussian2D_pState dWdx_ave, dWdy_ave, dWdx, dWdy, dWds;
  Vector2D dX;
  double ds;

  // Compute the Cartesian components of the intermediate state
  // solution viscous flux.
  dWdx_ave = HALF*(dW1dx + dW2dx);
  dWdy_ave = HALF*(dW1dy + dW2dy);

  dX = X2-X1; ds = dX.abs(); dX /= ds;

  dWds = (W2-W1)/ds;

  dWdx = dWdx_ave + (dWds - dWdx_ave*dX.x)*norm_dir.x/dot(norm_dir,dX);
  dWdy = dWdy_ave + (dWds - dWdy_ave*dX.y)*norm_dir.y/dot(norm_dir,dX);

  // Return the intermediate state solution viscous flux.
  return HeatFlux_n(X,W,dWdx,dWdy,norm_dir,Axisymmetric);
  //return Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO);

}
#endif

/**********************************************************************
 * Routine: FlatPlate                                                 *
 *                                                                    *
 * This function returns the exact solution for the flow over a flat  *
 * (adiabatic) plate (Blasius solution) at a given the position and   *
 * the freestream flow velocity.                                      *
 *                                                                    *
 **********************************************************************/
Gaussian2D_pState FlatPlate(const Gaussian2D_pState &Winf,
			    const Vector2D X,
			    double &eta,
			    double &f,
			    double &fp,
			    double &fpp) {

  int next_is_ua(0), next_is_ub(0), next_is_n(0), 
      ua_finished(0), ub_finished(0), n_finished(0); 
  Gaussian2D_pState W;
  double fo, n, dn, dn_original, k1, k2, k3, k4, h, dy, u_a, u_b;

  // Initialize variables.
  W.d = Winf.d; 
  W.p = Winf.p;
  W.erot = Winf.erot;
  eta = ZERO;
  f = ZERO; fo = ZERO; fp = ZERO; fpp = 0.33206;
  dn_original = 0.00005;
  dn = dn_original;

  h = 0.00005;                            //h and dy used for centered difference to
  dy = h/sqrt(Winf.v.x/(X.x*Winf.nu()));  //find du/dy for shear pressure calculation

  // Return upstream conditions before flat plate.
  if (X.x < ZERO) return Winf;

  // Return Winf if below plate
  if(X.y < ZERO) return Winf;

  // Return upstream conditions with zero velocity at the leading edge
  // of the plate.
  if (X.x < TOLER) return W;

  // Determine the dimensionless similarity coordinate, eta:
  eta = X.y*sqrt(Winf.v.x/(X.x*Winf.nu()));

  //if h is greater that half of eta, reset h to be half of eta
  if(h > 0.5*eta) {h = 0.5*eta;}

  // If eta is greater than 8.4, for the sake of expediency, use linear
  // extrapolation to determine the value of f (fp = ONE and fpp = ZERO)
  // given the tabulated value at 8.4 (note, the analytic solution is 
  // linear in this region).
  if (eta > 8.4) {
    fp = ONE; fpp = ZERO; f = 6.67923 + fp*(eta - 8.4);
    W.v.x = fp*Winf.v.x;
    W.v.y = HALF*sqrt(Winf.nu()*Winf.v.x/max(X.x,TOLER))*(eta*fp-f);
    return W;
  }

  // Compute the Blasius solution using a fourth-order Runge-Kutta method.
  while(!ub_finished) {

    fo = f;

    // Increment f:
    k1 = dn*fp;
    k2 = dn*(fp+k1/2);
    k3 = dn*(fp+k2/2);
    k4 = dn*(fp+k3);
    f = f + k1/6 + k2/3 + k3/3 + k4/6;
 
    // Increment fp:
    k1 = dn*fpp;
    k2 = dn*(fpp+k1/2);
    k3 = dn*(fpp+k2/2);
    k4 = dn*(fpp+k3);
    fp = fp+ k1/6 + k2/3 + k3/3 + k4/6;

    // Increment fpp:
    k1 = -dn*fo*fpp/2;
    k2 = -dn*(fo+dn/2)*(fpp+k1/2)/2;
    k3 = -dn*(fo+dn/2)*(fpp+k2/2)/2;
    k4 = -dn*(fo+dn)*(fpp+k3)/2;
    fpp = fpp + k1/6 + k2/3 + k3/3 + k4/6;

    // Increment n:
    n += dn;

    //calculate needed states
    if(next_is_ua) {
      next_is_ua = 0;
      ua_finished = 1;

      // Compute the velocity vector at point X.
      u_a = fp*Winf.v.x;

    } else if(next_is_n) {
      next_is_n = 0;
      n_finished = 1;

      // Compute the velocity vector at point X.
      W.v.x = fp*Winf.v.x;
      W.v.y = HALF*sqrt(Winf.nu()*Winf.v.x/max(X.x,TOLER))*(eta*fp-f);

    } else if(next_is_ub) {
      next_is_ub = 0;
      ub_finished = 1;

      // Compute the velocity vector at point X.
      u_b = fp*Winf.v.x;

    }

    // Determine the next increment dn if not finished
    if(!ub_finished) {

      dn = dn_original;  //start with original dn and modify if needed

      if(n+dn > eta-h && !ua_finished) {
	dn = eta-h-n;
	next_is_ua = 1;
      } else if(n + dn > eta && !n_finished) {
	dn = eta - n;
	next_is_n = 1;
      } else if(n +dn > eta+h && !ub_finished) {
	dn = eta+h-n;
	next_is_ub = 1;
      }
    }

  }

  W.p.xy = -W.mu() * (u_b-u_a)/(2.0*dy);

  // Return W state.
  return W;

}


/********************************************************
 * Gaussian2D_cState::relax                             *
 ********************************************************/

void Gaussian2D_pState::relax(double deltat, int stage, const Gaussian2D_pState &W) {

  double tau_trans, tau_rot;
  double a, b, c;
  double pxx2(0.0), pxy2(0.0), pyy2(0.0), pzz2(0.0), erot2(0.0);
  double implicit_coeff1(0.0), implicit_coeff2(0.0), implicit_coeff3(0.0);
  double implicit_coeff4(0.0), implicit_coeff5(0.0), implicit_coeff6(0.0);

  //tau_trans = viscosity()/pressure();
  //tau_rot = 15.0/4.0*bulk_viscosity()/pressure();
  tau_trans = tt();
  tau_rot   = tr();

  if(atoms==GAUSSIAN_MONATOMIC){
    a = 3.0*(double)stage*tau_trans+deltat;
    b = deltat;
    c = 3.0*(deltat+(double)stage*tau_trans);

    implicit_coeff1 = a/c;
    implicit_coeff2 = b/c;
    implicit_coeff3 = 0.0;
    implicit_coeff4 = (double)stage*tau_trans/(deltat+(double)stage*tau_trans);
    implicit_coeff5 = 0.0;
    implicit_coeff6 = 0.0;
  }else{
    a = 15.0*(double)stage*(double)stage*tau_trans*tau_rot+13.0*(double)stage*tau_trans*deltat;
    a = a+5.0*deltat*tau_rot*(double)stage+3.0*deltat*deltat;
    b = deltat*(5.0*(double)stage*tau_rot+3.0*deltat-2.0*(double)stage*tau_trans);
    c = 15.0*((double)stage*tau_trans+deltat)*((double)stage*tau_rot+deltat);

    implicit_coeff1 = a/c;
    implicit_coeff2 = b/c;
    implicit_coeff3 = 2.0*deltat/(5.0*(deltat+(double)stage*tau_rot));
    implicit_coeff4 = (double)stage*tau_trans/(deltat+(double)stage*tau_trans);
    implicit_coeff5 = deltat/(5.0*(deltat+(double)stage*tau_rot));
    implicit_coeff6 = (2.0*deltat+5.0*(double)stage*tau_rot)/(5.0*(deltat+(double)stage*tau_rot));
  }

  if(stage == ONE){

    pxx2  = implicit_coeff1*p.xx+implicit_coeff2*p.yy+implicit_coeff2*p.zz+implicit_coeff3*erot;
    pxy2  = implicit_coeff4*p.xy;
    pyy2  = implicit_coeff1*p.yy+implicit_coeff2*p.xx+implicit_coeff2*p.zz+implicit_coeff3*erot;
    pzz2  = implicit_coeff1*p.zz+implicit_coeff2*p.xx+implicit_coeff2*p.yy+implicit_coeff3*erot;
    erot2 = implicit_coeff5*(p.xx+p.yy+p.zz)+implicit_coeff6*erot;

  }else{

    p.xx = p.xx-deltat*(2.0*W.p.xx-W.p.yy-W.p.zz)/(3.0*tau_trans)/2.0;
    p.xy = p.xy-deltat*W.p.xy/tau_trans/2.0;
    p.yy = p.yy-deltat*(2.0*W.p.yy-W.p.xx-W.p.zz)/(3.0*tau_trans)/2.0;
    p.zz = p.zz-deltat*(2.0*W.p.zz-W.p.xx-W.p.yy)/(3.0*tau_trans)/2.0;

    if(atoms == GAUSSIAN_DIATOMIC){

      p.xx = p.xx-2.0*deltat*(W.p.xx+W.p.yy+W.p.zz-3.0*W.erot)/(15.0*tau_rot)/2.0;
      p.yy = p.yy-2.0*deltat*(W.p.xx+W.p.yy+W.p.zz-3.0*W.erot)/(15.0*tau_rot)/2.0;
      p.zz = p.zz-2.0*deltat*(W.p.xx+W.p.yy+W.p.zz-3.0*W.erot)/(15.0*tau_rot)/2.0;
      erot = erot-deltat*(3.0*W.erot-W.p.xx-W.p.yy-W.p.zz)/(5.0*tau_rot)/2.0;
    }

    pxx2  = implicit_coeff1*p.xx+implicit_coeff2*p.yy+implicit_coeff2*p.zz+implicit_coeff3*erot;
    pxy2  = implicit_coeff4*p.xy;
    pyy2  = implicit_coeff1*p.yy+implicit_coeff2*p.xx+implicit_coeff2*p.zz+implicit_coeff3*erot;
    pzz2  = implicit_coeff1*p.zz+implicit_coeff2*p.xx+implicit_coeff2*p.yy+implicit_coeff3*erot;
    erot2 = implicit_coeff5*(p.xx+p.yy+p.zz)+implicit_coeff6*erot;

  }
  
  p.xx = pxx2;
  p.xy = pxy2;
  p.yy = pyy2;
  p.zz = pzz2;
  erot = erot2;

}

/**********************************************************************
 * Routine: Free_Molecular_exact                                      *
 *                                                                    *
 * This function returns the exact solution for free molecular flow   *
 * past a convex shape at a point defined by "position".              *
 *                                                                    *
 **********************************************************************/

Gaussian2D_pState Free_Molecular_Exact(const Vector2D position,
				       const Gaussian2D_pState &Winf,
				       Vector2D *nodes,
				       const int number_of_nodes){
  int i, emergency_count(0);
  double temp_angle1, temp_angle2,
         max_angle(-THOUSAND), min_angle(THOUSAND);
  double beta_free, cos_angle;
  Gaussian2D_pState solution = Gaussian2D_W_VACUUM;
  Gaussian2D_pState reflected;

  //Determine angles that are exposed to free-stream flow

  for(i = 0; i < number_of_nodes; i++) {
    temp_angle1 = atan2( position.y, position.x ); 
    temp_angle2 = atan2( nodes[i].y, nodes[i].x ); 

    if((fabs(temp_angle1-temp_angle2) > TOLER) &&              //This ensures that an angle between two 
       (fabs(temp_angle1-temp_angle2+2.0*PI) > TOLER) &&       //identicle points is not taken
       (fabs(temp_angle1-temp_angle2-2.0*PI) > TOLER)) {
      temp_angle2 = atan2(position.y-nodes[i].y, position.x-nodes[i].x);
    }

    //deal with -PI == PI problems

    if(fabs(temp_angle2-temp_angle1)>PI/2.0+TOLER && temp_angle2>temp_angle1) {
      temp_angle2 -= 2.0*PI;
    } else if(fabs(temp_angle2-temp_angle1)>PI/2.0+TOLER && temp_angle2<temp_angle1) {
      temp_angle2 += 2.0*PI;
    } 

    if(fabs(temp_angle2-temp_angle1)/(PI/2.0) > 1.0+TOLER) {
      cout << "Horrible nightmarish trouble finding correct angles in \"Free_Molecular_Exact\"!" << endl;
    }

    if(temp_angle2>max_angle) max_angle = temp_angle2;
    if(temp_angle2<min_angle) min_angle = temp_angle2;

  }

  temp_angle1 = min_angle;
  min_angle   = max_angle;
  max_angle = temp_angle1+2.0*PI;

  //Angle of free-stream influence determined

  //Integrate free-stream influence:
  Integrate_distribution(solution, min_angle, max_angle, Winf);

  //Integrate reflected influence from each boundary cell face:

  for(i = 0; i < number_of_nodes; i++) {

    //find angle normal to surface and from surface midpoint
    //to position.  If difference < PI/2.0, particles
    //can reflect from surgace and get to position

    if( i==0 ) { //special care needed to link fist node to last node on closed surface
      //normal to surface
      temp_angle1 = atan2(nodes[number_of_nodes-1].x-nodes[i].x,
			  nodes[i].y-nodes[number_of_nodes-1].y);
      //angle to position
      temp_angle2 = atan2(position.y-(nodes[number_of_nodes-1].y+nodes[i].y)/2.0,
			  position.x-(nodes[number_of_nodes-1].x+nodes[i].x)/2.0);
    } else {
      //normal to surface
      temp_angle1 = atan2(nodes[i-1].x-nodes[i].x,
			  nodes[i].y-nodes[i-1].y);
      //angle to position
      temp_angle2 = atan2(position.y-(nodes[i-1].y+nodes[i].y)/2.0,
			  position.x-(nodes[i-1].x+nodes[i].x)/2.0);
    }

    //deal with PI==-PI problem
    if(fabs(temp_angle2-temp_angle1)>PI+TOLER && temp_angle2>temp_angle1) {
      temp_angle2 -= 2.0*PI;
    } else if(fabs(temp_angle2-temp_angle1)>PI+TOLER && temp_angle2<temp_angle1) {
      temp_angle2 += 2.0*PI;
    } 
 
    //check if particles can go from surface to position
    if(fabs(temp_angle1-temp_angle2) < PI/2.0 + TOLER) {
      //find min and max angles
      temp_angle1 = atan2( position.y, position.x ); 
      if( i==0 ) {
	temp_angle2 = atan2(nodes[number_of_nodes-1].y,nodes[number_of_nodes-1].x);
      } else {
	temp_angle2 = atan2( nodes[i-1].y, nodes[i-1].x ); 
      }

      if((fabs(temp_angle1-temp_angle2) > TOLER) &&              //This ensures that an angle between two 
	 (fabs(temp_angle1-temp_angle2+2.0*PI) > TOLER) &&       //identicle points is not taken
	 (fabs(temp_angle1-temp_angle2-2.0*PI) > TOLER)) {
	if( i==0 ) {
	  temp_angle2 = atan2(position.y-nodes[number_of_nodes-1].y,position.x-nodes[number_of_nodes-1].x);
	} else {
	  temp_angle2 = atan2(position.y-nodes[i-1].y, position.x-nodes[i-1].x);
	}
      }
      min_angle = temp_angle2;

      temp_angle1 = atan2( position.y, position.x ); 
      temp_angle2 = atan2( nodes[i].y, nodes[i].x ); 

      if((fabs(temp_angle1-temp_angle2) > TOLER) &&              //This ensures that an angle between two 
	 (fabs(temp_angle1-temp_angle2+2.0*PI) > TOLER) &&       //identicle points is not taken
	 (fabs(temp_angle1-temp_angle2-2.0*PI) > TOLER)) {
	temp_angle2 = atan2(position.y-nodes[i].y, position.x-nodes[i].x);
      }
      max_angle = temp_angle2;

      if(max_angle-min_angle > PI+TOLER) min_angle += 2.0*PI;
      if(min_angle-max_angle > PI+TOLER) max_angle += 2.0*PI;

      if(max_angle < min_angle) {
	temp_angle1 = max_angle;
	max_angle = min_angle;
	min_angle = temp_angle1;
      }

      if(fabs(max_angle-min_angle)/PI > 1.0+TOLER) {
	cout << "Horrible nightmarish trouble finding correct angles in \"Free_Molecular_Exact\"!" << endl;
      }

      reflected = Winf;   //For isothermal cylinder, do not change pressure
      reflected.v = Vector2D_ZERO;

      beta_free = Winf.d/(2.0*Winf.p.xx);
      if(i==0) {
	cos_angle = (nodes[number_of_nodes-1].y-nodes[i].y)/abs(nodes[i]-nodes[number_of_nodes-1]);;
      }else{
	cos_angle = (nodes[i-1].y-nodes[i].y)/abs(nodes[i]-nodes[i-1]);
      }

      reflected.d = Winf.d*sqrt(Winf.p.xx/reflected.p.xx)*(exp(-beta_free*Winf.v.x*Winf.v.x*cos_angle*cos_angle)
							   +Winf.v.x*cos_angle*sqrt(beta_free*PI)*
							   (1.0+erf(sqrt(beta_free)*Winf.v.x*cos_angle)));

      //integrate reflected distribution
      Integrate_distribution(solution, min_angle, max_angle, reflected);
    }

  }

  //Convert from Conserved to Primitive Variables

  solution.v.x = solution.v.x/solution.d;
  solution.v.y = solution.v.y/solution.d;
  solution.p.xx = solution.p.xx-solution.d*solution.v.x*solution.v.x;
  solution.p.xy = solution.p.xy-solution.d*solution.v.x*solution.v.y;
  solution.p.yy = solution.p.yy-solution.d*solution.v.y*solution.v.y;

  return solution;
}

/**********************************************************************
 * Routine: Integrate_distribution                                    *
 *                                                                    *
 * This function integrates a maxwellian distribution to determine    *
 * several moments of interest (used for free-molecular solutions)    *
 *                                                                    *
 **********************************************************************/

void Integrate_distribution(Gaussian2D_pState &solution,
			    const double min_angle,
			    const double max_angle,
			    const Gaussian2D_pState &W) {
  int i, num;
  double theta(0.0), ctheta(0.0);
  double delta_theta(0.0);
  const double m = W.M;
  const double n = W.d/m;
  const double u = W.v.x;
  const double beta(m*n/(2.0*W.p.xx));          //using p.xx as p_thermodynamic...should always work
  const double root_beta(sqrt(beta));
  const double root_pi(sqrt(PI));

  num = (int)((max_angle-min_angle)/(2.0*PI)*10000);
  if(num < 10) num = 10;

  delta_theta = (max_angle-min_angle)/((double)num);

  //Use mid point rule to integrate.

  for(i=0;i<num;i++){
    theta = min_angle+((double)i+0.5)*delta_theta;
    ctheta = cos(theta);

    solution.d = solution.d +
                 m*n*exp(-beta*u*u)/(2.0*PI) *
                 (1.0+root_beta*root_pi*u*ctheta*exp(beta*ctheta*ctheta*u*u)*(1.0+erf(root_beta*ctheta*u))) *
                 delta_theta;

    solution.v.x = solution.v.x + 
                   0.25*m*n*ctheta/(pow(PI,1.5)*root_beta)*exp(-beta*u*u) *
                   (2.0*ctheta*u*root_beta*root_pi + exp(beta*ctheta*ctheta*u*u) *
	   	    (2.0*PI*beta*ctheta*ctheta*u*u + PI + PI*erf(root_beta*ctheta*u)*(1.0+2.0*beta*ctheta*ctheta*u*u))) *
                   delta_theta;

    solution.v.y = solution.v.y+ 
                   0.25*m*n*sin(theta)/(pow(PI,1.5)*root_beta)*exp(-beta*u*u) *
                   (2.0*ctheta*u*root_beta*root_pi + exp(beta*ctheta*ctheta*u*u) *
		    (2.0*PI*beta*ctheta*ctheta*u*u + PI + PI*erf(root_beta*ctheta*u)*(1.0+2.0*beta*ctheta*ctheta*u*u))) *
                   delta_theta;

    solution.p.xx = solution.p.xx +
                    0.25*m*n*ctheta*ctheta/(PI*root_pi*beta*root_beta)*exp(-beta*u*u) *
                    (2.0*root_pi*root_beta+2.0*ctheta*ctheta*u*u*root_pi*root_beta*beta +
		     ctheta*u*PI*beta*exp(beta*ctheta*ctheta*u*u)*(3.0+2.0*ctheta*ctheta*u*u*beta +
								   erf(root_beta*ctheta*u)*(3.0+2.0*beta*ctheta*ctheta*u*u))) *
                    delta_theta;

    solution.p.xy = solution.p.xy +
                    0.25*m*n*sin(theta)*ctheta/(PI*root_pi*beta*root_beta)*exp(-beta*u*u) *
                    (2.0*root_pi*root_beta+2.0*ctheta*ctheta*u*u*root_pi*root_beta*beta +
		     ctheta*u*PI*beta*exp(beta*ctheta*ctheta*u*u)*(3.0+2.0*ctheta*ctheta*u*u*beta +
								   erf(root_beta*ctheta*u)*(3.0+2.0*beta*ctheta*ctheta*u*u))) *
                    delta_theta;

    solution.p.yy = solution.p.yy +
                    0.25*m*n*sin(theta)*sin(theta)/(PI*root_pi*beta*root_beta)*exp(-beta*u*u) *
                    (2.0*root_pi*root_beta+2.0*ctheta*ctheta*u*u*root_pi*root_beta*beta +
		     ctheta*u*PI*beta*exp(beta*ctheta*ctheta*u*u)*(3.0+2.0*ctheta*ctheta*u*u*beta +
								   erf(root_beta*ctheta*u)*(3.0+2.0*beta*ctheta*ctheta*u*u))) *
                    delta_theta;

    solution.p.zz = solution.p.zz +
                    m*n*exp(-beta*u*u)/(4.0*PI*beta) *
                    (1.0+root_beta*root_pi*u*ctheta*exp(beta*ctheta*ctheta*u*u)*(1.0+erf(root_beta*ctheta*u))) *
                    delta_theta;

  }

  return;
}

/**********************************************************************
 * Routine: RinglebFlow                                               *
 *                                                                    *
 * This function returns the average exact solution to Ringleb's flow *
 * for a quadrilateral cell defined by points X1, X2, X3, and X4      *
 * which are defined in a counter-clockwise direction starting in the *
 * south-west corner.  The average solution is determined using a     *
 * ?????th-order numerical quadrature.                                *
 *                                                                    *
 **********************************************************************/
Gaussian2D_pState RinglebFlowAverageState(const Gaussian2D_pState &Wdum,
					  const Vector2D &Y1,
					  const Vector2D &Y2,
					  const Vector2D &Y3,
					  const Vector2D &Y4) {
  //This has not been implemented.  It is only here for
  //compatibility with embeddedboundaries2D.
  //      ~james
  return Gaussian2D_pState();
}

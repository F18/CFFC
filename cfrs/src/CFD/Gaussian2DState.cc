/* Gaussian2DState.cc:  Subroutines for 2D Gaussian Solution State Classes. */

/* Include 2D Gaussian solution state header file. */

#ifndef _Gaussian2D_STATE_INCLUDED
#include "Gaussian2DState.h"
#endif // _Gaussian2D_STATE_INCLUDED

/****************************************************************
 * Gaussian2D_pState -- Create storage and assign gas constants.*
 ****************************************************************/
double Gaussian2D_pState::M     = MOLE_WT_AIR;
int    Gaussian2D_pState::atoms = GAUSSIAN_DIATOMIC;
int    Gaussian2D_pState::gas   = GAS_AIR;
double Gaussian2D_pState::alpha = ZERO;
/****************************************************************
 * Gaussian2D_cState -- Create storage and assign gas constants.*
 ****************************************************************/
double Gaussian2D_cState::M     = MOLE_WT_AIR;
int    Gaussian2D_cState::atoms = GAUSSIAN_DIATOMIC;
int    Gaussian2D_cState::gas   = GAS_AIR;
double Gaussian2D_cState::alpha = ZERO;

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

    //if(Ave.invalid()){
    //  Ave.make_valid();
    //}

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
 * Routine: Adiabatic_Wall                              *
 *                                                      *
 * This function returns the solution state for a wall  *
 * in a given direction given the primitive solution    *
 * variables and the unit normal vector in the          *
 * direction of interest.                               *
 *                                                      *
 ********************************************************/
Gaussian2D_pState Adiabatic_Wall(const Gaussian2D_pState &W,
				 const Gaussian2D_pState &Wo,
				 const Vector2D &norm_dir) {

    double dr, ur, vr, pxxr, pxyr, pyyr, pzzr, erotr;
    double ur2, vr2, pxxr2, pxyr2, pyyr2;
    double u, v, pxx, pxy, pyy;
    Tensor2D pr;
    double uWallr, vWallr;
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

    uWallr =  Wo.v.x*cos_angle + Wo.v.y*sin_angle;
    vWallr =  - Wo.v.x*sin_angle + Wo.v.y*cos_angle;

    // Adjust state.

    //  Mine 1
    /*
    ur2   = -ur;
    vr2   = vr+W.alpha*(vWallr-vr)+2.0*(2.0-W.alpha)*pxyr/sqrt(2*PI*dr*pxxr);
    pxxr2 = pxxr;
    pxyr2 = (W.alpha-1.0)*pxyr+2.0*W.alpha*sqrt(dr*pxxr/(2.0*PI))*(vr-vWallr);
    pyyr2 = pyyr - 2.0/3.0*W.d*sqr(W.alpha/2.0*(vWallr-vr)+(2.0-W.alpha)*pxyr/sqrt(2*PI*dr*pxxr))
                 + 2.0*((2.0-W.alpha)*(dr/2.0*sqr(vr-vr2)+pxyr*sqrt(2.0*dr/(PI*pxxr))*(vr-vr2))
                         +W.alpha*dr/2.0*sqr(vWallr-vr2));
    */
    //  Mine 2
    
    ur2   = -ur;
    vr2   = vr;
    pxxr2 = pxxr;
    pxyr2 = (W.alpha-1.0)*pxyr+2.0*W.alpha*sqrt(dr*pxxr/(2.0*PI))*(vr-vWallr);
    pyyr2 = pyyr;
    
    //  Dr. Groth's
    /*        
    ur2   = -ur;
    vr2   = vr;
    pxxr2 = pxxr;
    pxyr2 = -pxyr+4*(W.alpha)/((2.0-W.alpha)*sqrt(2.0*PI))*dr*(vr-vWallr)
                                                   *sqrt(pxxr/dr+ur*ur/3);
    pyyr2 = pyyr;
    */

    // Rotate back to the original Cartesian reference frame.

    sin_angle = -sin_angle;

    u   = ur2*cos_angle + vr2*sin_angle;
    v   = - ur2*sin_angle + vr2*cos_angle;
    pxx = pxxr2*cos_angle*cos_angle+pyyr2*sin_angle*sin_angle
          +2.0*pxyr2*cos_angle*sin_angle;
    pxy = pxyr2*(cos_angle*cos_angle-sin_angle*sin_angle)
          -(pxxr2-pyyr2)*cos_angle*sin_angle;
    pyy = pxxr2*sin_angle*sin_angle+pyyr2*cos_angle*cos_angle
          -2.0*pxyr2*cos_angle*sin_angle;

    /*
    if(Gaussian2D_pState(dr, u, v, pxx, pxy, pyy, pzzr, erotr).invalid()){
      cout << "Invalid Adiabatic Wall State:" << endl;
      cout << "W   = " << Gaussian2D_pState(dr, u, v, pxx, pxy, pyy, pzzr, erotr) << endl;
      cout << "Win =" << W << endl;
      cout << "Rotated:" << endl;
      cout << "W   = " << Rotate(Gaussian2D_pState(dr, u, v, pxx, pxy, pyy, pzzr, erotr),norm_dir) << endl;
      cout << "Win =" << Rotate(W,norm_dir) << endl;
    }
    */
    // Return the reflected state for adiabatic wall.

    return (Gaussian2D_pState(dr, u, v, pxx, pxy, pyy, pzzr, erotr));
       
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

Gaussian2D_pState BC_Characteristic_Velocity(const Gaussian2D_pState &Wi,
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
        return (Gaussian2D_pState(Wi.d, -Wi.v.x, Wi.v.y, Wi.p, Wi.erot));
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
        return (Gaussian2D_pState(Wi.d, -Wi.v.x, Wi.v.y, Wi.p, Wi.erot));
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
        return (Gaussian2D_pState(Wi.d, Wi.v.x, -Wi.v.y, Wi.p, Wi.erot));
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
//   return (Gaussian2D_pState(HartenFixPos(lambdas_a[1],
//                                       lambdas_l[1],
//                                       lambdas_r[1]),
//                          HALF*(lambdas_a[2]+fabs(lambdas_a[2])),
//                          HALF*(lambdas_a[3]+fabs(lambdas_a[3])),
//                          HALF*(lambdas_a[4]+fabs(lambdas_a[4])),
//                          HALF*(lambdas_a[5]+fabs(lambdas_a[5])),
//                          HALF*(lambdas_a[6]+fabs(lambdas_a[6])),
//                          HALF*(lambdas_a[7]+fabs(lambdas_a[7])),
//                          HartenFixPos(lambdas_a[8],
//                                       lambdas_l[8],
//                                       lambdas_r[8])));
  return lambdas_l;
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
//   return (Gaussian2D_pState(HartenFixNeg(lambdas_a[1],
//                                       lambdas_l[1],
//                                       lambdas_r[1]),
//                          HALF*(lambdas_a[2]-fabs(lambdas_a[2])),
//                          HALF*(lambdas_a[3]-fabs(lambdas_a[3])),
//                          HALF*(lambdas_a[4]-fabs(lambdas_a[4])),
//                          HALF*(lambdas_a[5]-fabs(lambdas_a[5])),
//                          HALF*(lambdas_a[6]-fabs(lambdas_a[6])),
//                          HALF*(lambdas_a[7]-fabs(lambdas_a[7])),
//                          HartenFixNeg(lambdas_a[8],
//                                       lambdas_l[8],
//                                       lambdas_r[8])));
  return lambdas_l;
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
//   return (Gaussian2D_pState(HartenFixAbs(lambdas_a[1],
//                                       lambdas_l[1],
//                                       lambdas_r[1]),
//                          fabs(lambdas_a[2]),
//                          fabs(lambdas_a[3]),
//                          fabs(lambdas_a[4]),
//                          fabs(lambdas_a[5]),
//                          fabs(lambdas_a[6]),
//                          fabs(lambdas_a[7]),
//                          HartenFixAbs(lambdas_a[8],
//                                       lambdas_l[8],
//                                       lambdas_r[8])));
  return lambdas_l;
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

  Gaussian2D_pState W;
  double fo, dn, sign, k1, k2, k3, k4;

  // Initialize variables.
  W.d = Winf.d; 
  W.p = Winf.p;
  W.erot = Winf.erot;
  eta = ZERO;
  f = ZERO; fo = ZERO; fp = ZERO; fpp = 0.33206;
  sign = ONE;
  dn = 0.000005;

  // Return upstream conditions before flat plate.
  if (X.x < ZERO) return Winf;

  // Return upstream conditions with zero velocity at the leading edge
  // of the plate.
  if (X.x < TOLER) return W;

  // Determine the dimensionless similarity coordinate, eta:
  eta = X.y*sqrt(Winf.v.x/(X.x*Winf.nu()));

  // If eta is greater than 8.4, for the sake of expediency, use linear
  // extrapolation to determine the value of f (fp = ONE and fpp = ZERO)
  // given the tabulated value at 8.4 (note, the analytic solution is 
  // linear in this region).
  if (eta > 8.4) {
    fp = ONE; fpp = ZERO; f = 6.67923 + fp*(eta - 8.4);
    W.v.x = fp*Winf.v.x;
    W.v.y = HALF*(eta*fp-f);
    return W;
  }

  // Compute the Blasius solution using a fourth-order Runge-Kutta method.
  for (double n = ZERO; n < eta; n += dn) {

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

  }

  // Compute the velocity vector at point X.
  W.v.x = fp*Winf.v.x;
  W.v.y = HALF*(eta*fp-f);

  // Return W state.
  return W;

}

//Test function
Gaussian2D_pState GaussianTestFunction (const double x){

  Gaussian2D_pState temp(2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0);

  //  return temp;
  return (Gaussian2D_pState)x;

}






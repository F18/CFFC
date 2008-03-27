/**********************************************************************
 **********************************************************************
 **                                                                  **
 ** File: ControlAngle.h                                             **
 **                                                                  **
 ** Description: This file defines some necessary functions for      **
 **              the computation of properties of control angle      **
 **              elements.                                           **
 **                                                                  **
 ** Author: Marc "T-Bone" Charest                                    **
 **                                                                  **
 ** Revision:  Date        Initials   Change                         **
 **            12/08/2007  MRC        Original creation              **
 **                                                                  **
 **********************************************************************
 **********************************************************************/
#ifndef _CONTROL_ANGLE_INCLUDED
#define _CONTROL_ANGLE_INCLUDED 


/********************************************************
 * Routine: Control_Angle_Avg                           *
 *                                                      *
 * Compute the control angle element parameters given   *
 * it's dimensions (polar, azimuthal) using exact       *
 * integration.  Returns the element averaged direction *
 * cosines and the element size                         *
 *                                                      *
 ********************************************************/
inline void Control_Angle_Avg( const double &theta_lo, // lower polar angle
			       const double &theta_hi, // upper polar angle
			       const double &psi_lo,   // lower azimuthal angle
			       const double &psi_hi,   // upper azimuthal angle
			       double &mu,             // (x,r)-direction cosine
			       double &eta,            // (y,psi)-direction cosine
			       double &xi,             // z-direction cosine
			       double &omega )         // solid angle element size
{
  // integrate over the polar angle
  double temp1 = ( 0.5*(theta_hi-theta_lo) -
		   0.25*( sin(2*theta_hi) - sin(2*theta_lo) ) );
  double temp2 = -0.25*( cos(2.0*theta_hi) - cos(2.0*theta_lo) );

  // compute
  mu    = temp1 * ( sin(psi_hi) - sin(psi_lo)  );   //x
  eta   = - temp1 * ( cos(psi_hi) - cos(psi_lo)  ); //y
  xi    = temp2 * ( psi_hi - psi_lo );              //z
  omega = - ( ( cos(theta_hi) - cos(theta_lo) ) *
	      ( psi_hi - psi_lo ) );

}


/********************************************************
 * Routine: Control_Angle_Ctr                           *
 *                                                      *
 * Compute the control angle element parameters given   *
 * it's dimensions (polar, azimuthal) using exact       *
 * integration.  Returns the direction cosines at the   *
 * control angle centroid and it's element size.        *
 *                                                      *
 ********************************************************/
inline void Control_Angle_Ctr( const double &theta_lo, // lower polar angle
			       const double &theta_hi, // upper polar angle
			       const double &psi_lo,   // lower azimuthal angle
			       const double &psi_hi,   // upper azimuthal angle
			       double &mu,             // (x,r)-direction cosine
			       double &eta,            // (y,psi)-direction cosine
			       double &xi,             // z-direction cosine
			       double &omega )         // solid angle element size
{
  // integrate over the polar angle
  double temp1 = ( 0.5*(theta_hi-theta_lo) -
		   0.25*( sin(2*theta_hi) - sin(2*theta_lo) ) );
  double temp2 = -0.25*( cos(2.0*theta_hi) - cos(2.0*theta_lo) );

  // centroid
  double theta_ctr = 0.5*(theta_lo+theta_hi);
  double psi_ctr   = 0.5*(psi_lo+psi_hi);

  // compute
  mu    = sin(theta_ctr) * cos(psi_ctr);           //x
  eta   = sin(theta_ctr) * sin(psi_ctr);           //y
  xi    =                  cos(theta_ctr);         //z
  omega = - ( ( cos(theta_hi) - cos(theta_lo) ) *
	      ( psi_hi - psi_lo ) );

}

#endif //end _CONTROLANGLE_INCLUDED 


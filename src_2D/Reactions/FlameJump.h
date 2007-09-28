/**********************************************************************
 **********************************************************************
 **                                                                  **
 ** File: FlameJump.h                                                **
 **                                                                  **
 ** Description: This file defines some routines to solve for the    **
 **              jump conditions of a 1D premixed flame. They are    **
 **              templated so LESPremixed2D and Chem2D can use them. **
 **                                                                  **
 ** Author: Samuel Wong                                              **
 **                                                                  **
 ** Revision:  Date        Initials   Change                         **
 **            27/09/2007  SW         Original creation              **
 **            27/09/2007  MRJ        Templated FlameJump Functions  **
 **                                                                  **
 **********************************************************************
 **********************************************************************/
#ifndef _FLAMEJUMP_INCLUDED
#define _FLAMEJUMP_INCLUDED 


/***********************************************************
 * Necessary  Includes                                     *
 ***********************************************************/
#include "../Math/Math.h"
#include "../Math/Vector2D.h"


/***********************************************************
 * Necessary  Constants                                    *
 ***********************************************************/


/***********************************************************
 * Routine: FlameJump_x                                    *
 *                                                         *
 * This routine solves a 1D premixed flame in the          *
 * x-direction using jump condition equations. Returns the *
 * burnt solution vector.                                  *
 ***********************************************************/
template <class SOLN_PSTATE>
inline int FlameJump_x(const SOLN_PSTATE &Wu, // unburnt state
		       SOLN_PSTATE &Wb)       // burnt state
{
  
  // Product species' mass fractions should already be set

  // Apply Newton's method to solve for rho2, u2, p2, T2.
  double rho2, u2, T2, p2, Cp2;   // Unknown variables.
  double f[4];                    // LHS vector.
  double delta_vars[4];           // Change in variables.
  double detJ, detf[4];           // Determinants.
  double R2;

  double norm, norm_tolerance;
  norm = MILLION;
  norm_tolerance = 1e-8;
  int iter_count, iter_max;
  iter_count = 0;
  iter_max = 20;

  // Initial guess.
  rho2 = Wu.rho;
  u2   = Wu.v.x;
  p2   = Wu.p;
  T2   = Wu.T();
  R2   = Wb.Rtot();

  while (norm > norm_tolerance) {

    f[0] = rho2*u2 - Wu.rho*Wu.v.x;
    f[1] = rho2*sqr(u2) + p2 - Wu.rho*sqr(Wu.v.x) - Wu.p;
    f[2] = Wb.h(T2) + HALF*sqr(u2) - Wu.h() - HALF*sqr(Wu.v.x);
    f[3] = rho2*R2*T2 - p2;

    Cp2 = Wb.hprime(T2);

    // Use Cramer's rule.
    detJ = -(u2*(u2*rho2*R2)+R2*T2*(rho2*Cp2)) + (Cp2*(TWO*u2*rho2*u2-sqr(u2)*rho2));
    detf[0] = -(-Cp2*f[3]*rho2+rho2*R2*(-f[0]*u2+f[2]*rho2)) + Cp2*(-f[0]*TWO*rho2*u2+f[1]*rho2);
    detf[1] = -(u2*(-f[2]*rho2*R2+f[3]*Cp2)+R2*T2*(-f[0]*Cp2)) + Cp2*(-f[1]*u2+sqr(u2)*f[0]);
    detf[2] = -Cp2*(-rho2*(-sqr(u2)*f[3]+R2*T2*f[1])+TWO*rho2*u2*(-u2*f[3]+R2*T2*f[0])) + rho2*R2*(u2*(-TWO*rho2*u2*f[2]+u2*f[1])-sqr(u2)*(-rho2*f[2]+u2*f[0]));
    detf[3] = -(-f[3]*sqr(u2)+R2*T2*(-f[2]*rho2+u2*f[0])) + (u2*(-TWO*rho2*u2*f[2]+u2*f[1])-sqr(u2)*(-rho2*f[2]+u2*f[0]));

    for (int i=0; i<=3; i++) {
      delta_vars[i] = detf[i]/detJ;
    }

    rho2 = rho2 + delta_vars[0];
    u2 = u2 + delta_vars[1];
    p2 = p2 + delta_vars[2];
    T2 = T2 + delta_vars[3];

    norm = sqrt( sqr(delta_vars[0]) + sqr(delta_vars[1]) + sqr(delta_vars[2]) + sqr(delta_vars[3]));

    iter_count++;
    if (iter_count>iter_max) {
      // Newton's method did not converge.
      return 1;
    }
  }

  Wb.rho = rho2;
  Wb.v.x = u2;
  Wb.v.y = Wu.v.y;   // Jump in tangential velocity is zero.
  Wb.p   = p2;

  // Burnt state successfully calculated.
  return 0;
}

/***********************************************************
 * Routine: FlameJump_n                                    *
 *                                                         *
 * This routine solves a premixed flame in any arbitrary   *
 * direction defined by a unit normal vector in the        *
 * direction of interest. The problem is solved by first   *
 * applying a frame rotation to rotate the problem to a    *
 * local frame aligned with the unit normal vector and     *
 * then by using the jump conditions to solve for the      *
 * burnt solution vector. The solution is then rotated     *
 * back to the original frame.                             *
 *                                                         *
 ***********************************************************/
template <class SOLN_PSTATE>
inline int FlameJump_n(const SOLN_PSTATE &Wu,     // unburnt state
		       SOLN_PSTATE &Wb,           // burnt state
		       const Vector2D &norm_dir)  // normal direction vec
{

  int error_flag;
  double cos_angle, sin_angle;
  SOLN_PSTATE Wu_rotated, Wb_rotated;

  /* Determine the direction cosine's for the frame rotation. */
  cos_angle = norm_dir.x;
  sin_angle = norm_dir.y;

  /* Apply the frame rotation and evaluate unburnt and burnt
     solution states in the local rotated frame defined
     by the unit normal vector. */
  Wu_rotated.Copy(Wu);
  Wb_rotated.Copy(Wb);
  Wu_rotated.v.x = Wu.v.x*cos_angle + Wu.v.y*sin_angle;
  Wu_rotated.v.y = - Wu.v.x*sin_angle + Wu.v.y*cos_angle;

  /* Evaluate the jump conditions in the rotated frame. */
  error_flag = FlameJump_x(Wu_rotated,Wb_rotated);
  if (error_flag) {
    return error_flag;
  }

  /* Rotate back to the original Cartesian reference frame. */
  Wb.Copy(Wb_rotated);
  Wb.v.x = Wb_rotated.v.x*cos_angle - Wb_rotated.v.y*sin_angle;
  Wb.v.y = Wb_rotated.v.x*sin_angle + Wb_rotated.v.y*cos_angle;

  // Burnt solution state sucessfully calculated.
  return 0;
}

/***********************************************************
 * Routine: FlameJump_x                                    *
 *                                                         *
 * This routine solves a 1D premixed flame in the          *
 * x-direction using jump condition equations. Returns the *
 * burnt solution vector.                                  *
 ***********************************************************/
template <class SOLN_PSTATE, class SOLN_CSTATE>
inline int FlameJump_x(const SOLN_CSTATE &Uu, // unburnt state
		       SOLN_CSTATE &Ub)       // burnt state
{
  int error_flag;

  SOLN_PSTATE Wu = Uu.W();
  SOLN_PSTATE Wb = Ub.W();

  error_flag = FlameJump_x(Wu, Wb);
  if (error_flag) {
    return error_flag;
  }

  Ub = Wb.U();

  // Burnt state successfully calculated.
  return 0;
}

/***********************************************************
 * Routine: FlameJump_n                                    *
 *                                                         *
 * This routine solves a premixed flame in any arbitrary   *
 * direction defined by a unit normal vector in the        *
 * direction of interest. The problem is solved by first   *
 * applying a frame rotation to rotate the problem to a    *
 * local frame aligned with the unit normal vector and     *
 * then by using the jump conditions to solve for the      *
 * burnt solution vector.                                  *
 *                                                         *
 ***********************************************************/
template <class SOLN_PSTATE, class SOLN_CSTATE>
inline int FlameJump_n(const SOLN_CSTATE &Uu,    // unburnt state
		       SOLN_CSTATE &Ub,          // burnt state
		       const Vector2D &norm_dir) // normal direction vec
{

  int error_flag;
  double cos_angle, sin_angle;
  SOLN_CSTATE Uu_rotated, Ub_rotated;

  /* Determine the direction cosine's for the frame rotation. */
  cos_angle = norm_dir.x; 
  sin_angle = norm_dir.y;

  /* Apply the frame rotation and evaluate unburnt and burnt
     solution states in the local rotated frame defined
     by the unit normal vector. */
  Uu_rotated.Copy(Uu);
  Ub_rotated.Copy(Ub);
  Uu_rotated.rhov.x =   Uu.rhov.x*cos_angle + Uu.rhov.y*sin_angle;
  Uu_rotated.rhov.y = - Uu.rhov.x*sin_angle + Uu.rhov.y*cos_angle;

  /* Evaluate the jump conditions in the rotated frame. */
  error_flag = FlameJump_x(Uu_rotated,Ub_rotated);
  if (error_flag) {
    return error_flag;
  }

  /* Rotate back to the original Cartesian reference frame. */
  Ub.Copy(Ub_rotated);
  Ub.rhov.x = Ub_rotated.rhov.x*cos_angle - Ub_rotated.rhov.y*sin_angle;
  Ub.rhov.y = Ub_rotated.rhov.x*sin_angle + Ub_rotated.rhov.y*cos_angle;

  // Burnt solution state sucessfully calculated.
  return 0;
}

/**********************************************************************
 * Routine: FlameJumpLowMach_x                                        *
 *                                                                    *
 * This routine solves a 1D premixed flame in the x-direction using   *
 * jump condition equations, subject to a low Mach number assumption, *
 * that is,                                                           *
 *                                                                    *
 *   rho1*T1=rho2*T2  (see Poinsot/Veynante, "Numerical Combustion")  *
 *                                                                    *
 * Returns the burnt mixture state vector.                            *
 *                                                                    *
 **********************************************************************/
template <class SOLN_PSTATE>
inline int FlameJumpLowMach_x(const SOLN_PSTATE &Wu, // unburnt state
			      SOLN_PSTATE &Wb)       // burnt state
{

  // Product species' mass fractions should already be set

  // Apply Newton's method to solve for rho2, u2, T2.
  double rho2, u2, T2, p2, Cp2;      // Unknown variables.
  double f[3];                       // LHS vector.
  double delta_vars[3];              // Change in variables.
  double detJ, detf0, detf1, detf2;  // Determinants.

  double norm, norm_tolerance;
  norm = MILLION;
  norm_tolerance = 1e-6;
  int iter_count, iter_max;
  iter_count = 0;
  iter_max = 20;

  // Initial guess.
  rho2 = DENSITY_STDATM;
  u2   = TWO*Wu.v.x;
  T2   = 1000;

  while (norm > norm_tolerance) {

    f[0] = -(rho2*u2 - Wu.rho*Wu.v.x);
    f[1] = -(rho2*T2 - Wu.rho*Wu.T());
    f[2] = -( (Wb.h(T2)+HALF*sqr(u2)) - (Wu.h()+HALF*sqr(Wu.v.x)) );

    Cp2 = Wb.hprime(T2);

    detJ  = u2*(-u2*rho2) - T2*(rho2*Cp2);
    detf0 = -rho2*(f[0]*u2-f[2]*rho2) + Cp2*(-f[1]*rho2);
    detf1 = u2*(f[1]*Cp2-f[2]*rho2) - T2*(f[0]*Cp2);
    detf2 = u2*(-u2*f[1]) - T2*(f[2]*rho2-f[0]*u2);

    delta_vars[0] = detf0/detJ;
    delta_vars[1] = detf1/detJ;
    delta_vars[2] = detf2/detJ;

    rho2 = rho2 + delta_vars[0];
    u2   = u2 + delta_vars[1];
    T2   = T2 + delta_vars[2];

    norm = sqrt( sqr(delta_vars[0]) + sqr(delta_vars[1]) + sqr(delta_vars[2]) );

    iter_count++;
    if (iter_count>iter_max) {
      // Newton method did not converge.
      return 1;
    }
  }

  Wb.rho = rho2;
  Wb.v.x = u2;
  Wb.v.y = Wu.v.y;   // Jump in tangential velocity is zero.
  Wb.p = Wu.p + Wu.rho*sqr(Wu.v.x)*(ONE - Wb.v.x/Wu.v.x);

  // Burnt state successfully calculated.
  return 0;
}


#endif //end _FLAMEJUMP_INCLUDED 

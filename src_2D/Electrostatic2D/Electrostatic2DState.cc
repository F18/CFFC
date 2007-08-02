/**********************************************************************
 * Dusty2DState.cc                                                    *
 *                                                                    *
 * Subroutines for 2D dusty solution state classes.                   *
 *                                                                    *
 **********************************************************************/

// Include 2D Dusty solution state header file.

#ifndef _ELECTROSTATIC2D_STATE_INCLUDED
#include "Electrostatic2DState.h"
#endif // _ELECTROSTATIC2D_STATE_INCLUDED

/**********************************************************************
 * Routine: Rotate                                                    *
 *                                                                    *
 * This function returns the solution in the local rotated frame.     *
 *                                                                    *
 **********************************************************************/
Electrostatic2DState Rotate(const Electrostatic2DState &U,
			    const Vector2D &norm_dir) {

  Electrostatic2DState Urotated;
  double cos_angle = norm_dir.x;  
  double sin_angle = norm_dir.y;

  Urotated.E.x =   U.E.x*cos_angle + U.E.y*sin_angle;
  Urotated.E.y = - U.E.x*sin_angle + U.E.y*cos_angle;
  Urotated.V = U.V;

  // Return the rotated state.
  return Urotated;

}

/**********************************************************************
 * Routine: Translate                                                 *
 *                                                                    *
 * This function returns the solution in a stationary frame.          *
 *                                                                    *
 **********************************************************************/
Electrostatic2DState Translate(const Electrostatic2DState &U,
			       const Vector2D &V) {

  Electrostatic2DState Utranslated;

  Utranslated = U;

  // Return the translated state.
  return Utranslated;

}

/**********************************************************************
 * Routine: Reflect                                                   *
 *                                                                    *
 * This function returns the reflected solution state in a given      *
 * direction given the primitive solution variables and the unit      *
 * normal vector in the direction of interest.                        *
 *                                                                    *
 **********************************************************************/
Electrostatic2DState Reflect(const Electrostatic2DState &U,
			     const Vector2D &norm_dir) {
  
  Electrostatic2DState Ur, Un;

  // Apply the frame rotation and calculate the primitive solution 
  // state variables in the local rotated frame defined by the unit 
  // normal vector.
  Ur = Rotate(U,norm_dir);

  // Reflect the electric field velocity in the rotated frame.
  Ur.E.x = -Ur.E.x;

  // Rotate back to the original Cartesian reference frame.
  Un = Rotate(Ur,Vector2D(norm_dir.x,-norm_dir.y));

  // Return the reflected state.
  return Un;

}

/**********************************************************************
 * Routine: Mirror                                                    *
 *                                                                    *
 * This function returns the mirrored solution state in a given       *
 * direction given the primitive solution variables and the unit      *
 * normal vector in the direction of interest.                        *
 *                                                                    *
 **********************************************************************/
Electrostatic2DState Mirror(const Electrostatic2DState &U,
			    const Vector2D &norm_dir) {

  Electrostatic2DState Ur, Un;

  // Apply the frame rotation and calculate the primitive solution 
  // state variables in the local rotated frame defined by the unit 
  // normal vector.
  Ur = Rotate(U,norm_dir);

  // Mirror the electric field velocity in the rotated frame.
  Ur.E.x = -Ur.E.x;
  Ur.E.y = -Ur.E.y;

  // Rotate back to the original Cartesian reference frame.
  Un = Rotate(Ur,Vector2D(norm_dir.x,-norm_dir.y));

  // Return the mirrored state.
  return Un;

}

/**********************************************************************
 * Routine: FluxArithmetic_n                                          *
 *                                                                    *
 * This function returns the intermediate state solution flux         *
 * calculated by the arithmetic mean of the cell-centred flux terms   *
 * of the neighbouring cells.                                         *
 *                                                                    *
 **********************************************************************/
Electrostatic2DState FluxArithmetic_n(const Vector2D &X,
				      const Vector2D &X1,
				      const Electrostatic2DState &U1,
				      const Electrostatic2DState &dU1dx,
				      const Electrostatic2DState &dU1dy,
				      const Vector2D &X2,
				      const Electrostatic2DState &U2,
				      const Electrostatic2DState &dU2dx,
				      const Electrostatic2DState &dU2dy,
				      const Vector2D &norm_dir) {

//   Electrostatic2DState dUdx, dUdy;
//   // Compute the Cartesian components of the intermediate state
//   // solution viscous flux.
//   dUdx.V = HALF*(dU1dx.V + dU2dx.V);
//   dUdy.V = HALF*(dU1dy.V + dU2dy.V);
//   // Return the intermediate state solution viscous flux.
//   return Electrostatic2DState(ZERO,ZERO,dUdx.V*norm_dir.x+dUdy.V*norm_dir.y);
  Electrostatic2DState dUdxb, dUdyb, dUdx, dUdy, dUdl;
  Vector2D dX;
  double l;
  // Compute the Cartesian components of the intermediate state
  // solution viscous flux.
  dUdxb.V = HALF*(dU1dx.V + dU2dx.V);
  dUdyb.V = HALF*(dU1dy.V + dU2dy.V);

  dX = X2-X1; l = dX.abs(); dX /= l;

  dUdl.V = (U2.V-U1.V)/l;

  dUdx.V = dUdxb.V - (dot(Vector2D(dUdxb.V,dUdyb.V),dX) - dUdl.V)*dX.x;
  dUdy.V = dUdyb.V - (dot(Vector2D(dUdxb.V,dUdyb.V),dX) - dUdl.V)*dX.y;

  // Return the intermediate state solution viscous flux.
  return Electrostatic2DState(ZERO,ZERO,dUdx.V*norm_dir.x+dUdy.V*norm_dir.y);

}

/**********************************************************************
 * Routine: FluxDiamondPath_n                                         *
 *                                                                    *
 * This routine computes the intermediate state solution flux at the  *
 * specified quadrature point, X.  The gradient of the primitive      *
 * variables is computed on the diamond-path defined by the points    *
 * X1, X2, X3, and X4.                                                *
 *                                                                    *
 **********************************************************************/
Electrostatic2DState FluxDiamondPath_n(const Vector2D &X,
				       const Vector2D &X1, const Electrostatic2DState &U1,
				       const Vector2D &X2, const Electrostatic2DState &U2,
				       const Vector2D &X3, const Electrostatic2DState &U3,
				       const Vector2D &X4, const Electrostatic2DState &U4,
				       const Vector2D &norm_dir,
				       const int &diamond_path) {

  Electrostatic2DState U_face, dUdx1, dUdy1, dUdx3, dUdy3, dUdx, dUdy;
  double A, A1, A3;
  Vector2D n21, n42, n14, n32, n43, n24;
  // Determine the lengths and normals of the faces and the areas of
  // the regions of Green-Gauss integration.
  n21 = Vector2D((X2.y-X1.y),-(X2.x-X1.x));
  n42 = Vector2D((X4.y-X2.y),-(X4.x-X2.x));
  n14 = Vector2D((X1.y-X4.y),-(X1.x-X4.x));
  n32 = Vector2D((X3.y-X2.y),-(X3.x-X2.x));
  n43 = Vector2D((X4.y-X3.y),-(X4.x-X3.x));
  n24 = Vector2D((X2.y-X4.y),-(X2.x-X4.x));
  A1 = HALF*((X2-X1)^(X4-X1));
  A3 = HALF*((X3-X4)^(X3-X2));
  A  = A1 + A3;
  // Compute Green-Gauss integration on triangle 1.
  if (diamond_path == DIAMONDPATH_QUADRILATERAL_RECONSTRUCTION ||
      diamond_path == DIAMONDPATH_LEFT_TRIANGLE_HEATFLUX) {
    U_face.V = HALF*(U2.V+U1.V);
    dUdx1.V = U_face.V*n21.x;
    dUdy1.V = U_face.V*n21.y;
    U_face.V = HALF*(U4.V+U2.V);
    dUdx1.V += U_face.V*n42.x;
    dUdy1.V += U_face.V*n42.y;
    U_face.V = HALF*(U1.V+U4.V);
    dUdx1.V += U_face.V*n14.x;
    dUdy1.V += U_face.V*n14.y;
    dUdx1.V = dUdx1.V/A1;
    dUdy1.V = dUdy1.V/A1;
  }
  // Compute Green-Gauss integration on triangle 3.
  if (diamond_path == DIAMONDPATH_QUADRILATERAL_RECONSTRUCTION ||
      diamond_path == DIAMONDPATH_RIGHT_TRIANGLE_HEATFLUX) {
    U_face.V = HALF*(U3.V+U2.V);
    dUdx3.V = U_face.V*n32.x;
    dUdy3.V = U_face.V*n32.y;
    U_face.V = HALF*(U4.V+U3.V);
    dUdx3.V += U_face.V*n43.x;
    dUdy3.V += U_face.V*n43.y;
    U_face.V = HALF*(U2.V+U4.V);
    dUdx3.V += U_face.V*n24.x;
    dUdy3.V += U_face.V*n24.y;
    dUdx3.V = dUdx3.V/A3;
    dUdy3.V = dUdy3.V/A3;
  }
  // Determine the area-averaged gradients of the primitive variables.
  if (diamond_path == DIAMONDPATH_QUADRILATERAL_RECONSTRUCTION) {
    dUdx.V = (A1*dUdx1.V + A3*dUdx3.V)/A;
    dUdy.V = (A1*dUdy1.V + A3*dUdy3.V)/A;
  } else if (diamond_path == DIAMONDPATH_LEFT_TRIANGLE_HEATFLUX) {
    dUdx.V = dUdx1.V;
    dUdy.V = dUdy1.V;
  } else if (diamond_path == DIAMONDPATH_RIGHT_TRIANGLE_HEATFLUX) {
    dUdx.V = dUdx3.V;
    dUdy.V = dUdy3.V;
  }
  // Return the solution flux.
  return Electrostatic2DState(ZERO,ZERO,dUdx.V*norm_dir.x+dUdy.V*norm_dir.y);

}

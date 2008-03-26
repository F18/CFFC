/* Vector2D.cc:  Subroutines for the 2D vector class. */

/* Include the 2D Vector header file. */

#ifndef _VECTOR2D_INCLUDED
#include "Vector2D.h"
#endif // _VECTOR2D_INCLUDED

/**********************************************************************
 * Routine: nhat                                                      *
 *                                                                    *
 * This routine computes the unit normal defined by two points.       *
 *                                                                    *
 **********************************************************************/
Vector2D nhat(const Vector2D X1, const Vector2D X2) {
  return Vector2D((X1.y-X2.y),-(X1.x-X2.x))/abs(X1-X2);
}

/**********************************************************************
 * Routine: Interpolate.                                              *
 *                                                                    *
 * Used linear interpolation to determine the location on a line      *
 * where f = 0.                                                       *
 *                                                                    *
 **********************************************************************/
Vector2D Interpolate(const Vector2D X1, const double f1,
		     const Vector2D X2, const double f2) {

  Vector2D X;
  double dx = NANO*abs(X2 - X1);

  if ((fabs(f1) < dx) && (fabs(f2) < dx)) {
    X = HALF*(X1+X2);
  } else if (fabs(f1) < dx) {
    X = X1;
  } else if (fabs(f2) < dx) {
    X = X2;
  } else {
    X = X1 - f1*(X2-X1)/(f2-f1);
  }

  // Return the interpolated point.
  return X;

}
Vector2D Interpolate(const Vector2D X1, const double f1,
		     const Vector2D X2, const double f2,
		     const double epsilon) {

  Vector2D X;

  if ((fabs(f1) < epsilon) && (fabs(f2) < epsilon)) {
    X = HALF*(X1+X2);
  } else if (fabs(f1) < epsilon) {
    X = X1;
  } else if (fabs(f2) < epsilon) {
    X = X2;
  } else {
    X = X1 - f1*(X2-X1)/(f2-f1);
  }

  // Return the interpolated point.
  return X;

}

/**********************************************************************
 * Routine: Line_Intersection                                         *
 *                                                                    *
 * Parameterization of two lines:                                     *
 *                                                                    *
 *   Ra = Xa1 + s*Xa2  and  Rb = Xb1 + t*Xb2                          *
 *                                                                    *
 * Set Ra = Rb:                                                       *
 *                                                                    *
 *   [ Xa2.x  - Xb2.x ] [ s ] = [ Xb1.x - Xa1.x ]                     *
 *   [ Xa2.y  - Xb2.y ] [ t ] = [ Xb1.y - Xa1.y ]                     *
 *                                                                    *
 * By Cramer's Rule:                                                  *
 *                                                                    *
 *   s = ((Xb1.y - Xa1.y)*Xb2.x - (Xb1.x - Xa1.x)*Xb2.y)/det          *
 *   t = ((Xb1.y - Xa1.y)*Xa2.x - (Xb1.x - Xa1.x)*Xa2.y)/det          *
 *                                                                    *
 * where det = Xa2.y*Xb2.x - Xa2.x*Xb2.y                              *
 *                                                                    *
 * If det = 0 then the lines are parallel (maybe coincident).         *
 *                                                                    *
 * Then the point of intersection is given by                         *
 *                                                                    *
 *   x = Xa1.x + s*Xa2.x = Xb1.x + t*Xb2.x                            *
 *   y = Xa1.y + s*Xa2.y = Xb1.y + t*Xb2.y                            *
 *                                                                    *
 * Note that this point of intersection may not be contained within   *
 * the two line segments.                                             *
 *                                                                    *
 **********************************************************************/
int Line_Intersection(const Vector2D Xa1,
		      const Vector2D Xa3,
		      const Vector2D Xb1,
		      const Vector2D Xb3,
		      Vector2D &Xp) {
  return Line_Intersection(Xa1,Xa3,Xb1,Xb3,Xp,NANO);
}

int Line_Intersection(const Vector2D Xa1,
		      const Vector2D Xa3,
		      const Vector2D Xb1,
		      const Vector2D Xb3,
		      Vector2D &Xp,
		      const double eps) {

  double det, s, t;
  Vector2D Xa2(Xa3 - Xa1);
  Vector2D Xb2(Xb3 - Xb1);

  //if (abs(Xa2) < NANO || abs(Xb2) < NANO) return 1;

  det = Xa2.y*Xb2.x - Xa2.x*Xb2.y;

  if (fabs(det)/min(Xa2*Xa2,Xb2*Xb2) < NANO) return 0;

  s = ((Xb1.y - Xa1.y)*Xb2.x - (Xb1.x - Xa1.x)*Xb2.y)/det;
  t = ((Xb1.y - Xa1.y)*Xa2.x - (Xb1.x - Xa1.x)*Xa2.y)/det;

  if (s*abs(Xa2) < -eps || (s-1.0)*abs(Xa2) > eps || t*abs(Xb2) < -eps || (t-1.0)*abs(Xb2) > eps) return 0;
  //if (s < ZERO-NANO || s > ONE+NANO || t < ZERO-NANO || t > ONE+NANO) return 0;
  //if (s <= ZERO || s >= ONE || t <= ZERO || t >= ONE) return 0;

  Xp = Xa1 + s*Xa2;

  // Return the intersection point.
  return 1;

}

int Line_Intersection(const Vector2D Xa1,
		      const Vector2D Xa3,
		      const Vector2D Xb1,
		      const Vector2D Xb3) {

  double det, s, t;
  Vector2D Xa2(Xa3 - Xa1);
  Vector2D Xb2(Xb3 - Xb1);

  //if (abs(Xa2) < NANO || abs(Xb2) < NANO) return 1;

  det = Xa2.y*Xb2.x - Xa2.x*Xb2.y;
  if (fabs(det)/min(Xa2*Xa2,Xb2*Xb2) < NANO) return 0;

  s = ((Xb1.y - Xa1.y)*Xb2.x - (Xb1.x - Xa1.x)*Xb2.y)/det;
  t = ((Xb1.y - Xa1.y)*Xa2.x - (Xb1.x - Xa1.x)*Xa2.y)/det;

  if ((s >= ZERO && s <= ONE) && (t >= ZERO && t <= ONE)) return 1;

  return 0;

}

/**********************************************************************
 * Routine: Triangle_Area                                             *
 *                                                                    *
 * Return the area of the triangle defined by the given three points. *
 *                                                                    *
 **********************************************************************/
double Triangle_Area(const Vector2D X1,
 		     const Vector2D X2,
 		     const Vector2D X3) {
  return HALF*((X1^X2) + (X2^X3) + (X3^X1));
}

/**********************************************************************
 * Routine: Quadrilateral_Area                                        *
 *                                                                    *
 * Return the area of the quadrilateral defined by the given four     *
 * points.                                                            *
 *                                                                    *
 **********************************************************************/
double Quadrilateral_Area(const Vector2D X1,
 			  const Vector2D X2,
 			  const Vector2D X3,
 			  const Vector2D X4) {
  return Triangle_Area(X1,X2,X3) + Triangle_Area(X4,X2,X3);
}

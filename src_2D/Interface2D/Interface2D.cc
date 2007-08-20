/**********************************************************************
 * Interface2D.cc: Interface subroutines.                             *
 **********************************************************************/

// Include the 2D Spline header file.
#ifndef _INTERFACE2D_INCLUDED
#include "Interface2D.h"
#endif // _INTERFACE2D_INCLUDED

/**********************************************************************
 * Interface2D::Define member functions.                              *
 **********************************************************************/

/**********************************************************************
 * Interface2D::Area -- Calculate the area of the interface.          *
 **********************************************************************/
double Interface2D::Area(void) {
  // Assert that the interface/spline is at least three-sided.
  assert(Spline.np > 2);
  // Assert that the interface/spline is closed.
  //assert(abs(Spline.Xp[0]-Spline.Xp[Spline.np-1]) < TOLER*TOLER);
  // Calculate area by summing the area of all of the sub-triangles.
  double A;
  A = ZERO;
  for (int n = 1; n < Spline.np-1; n++)
    A += HALF*((Spline.Xp[0]^Spline.Xp[n]) +
	       (Spline.Xp[n]^Spline.Xp[n+1]) +
	       (Spline.Xp[n+1]^Spline.Xp[0]));
  // Return the area.
  return A;
}

/**********************************************************************
 * Interface2D::Centroid -- Calculate the centroid of the interface.  *
 **********************************************************************/
void Interface2D::Centroid(void) {
  // Assert that the interface/spline is at least three-sided.
  assert(Spline.np > 2);
  // Assert that the interface/spline is closed.
  //assert(abs(Spline.Xp[0]-Spline.Xp[Spline.np-1]) < TOLER*TOLER);
  // Reset the reference/centroid location.
  Xref = Vector2D_ZERO;
  // Determine the centroid location by computing the area wieghted
  // average of the centroid of each sub-triangle.
  for (int n = 1; n < Spline.np-1; n++)
    Xref += (Spline.Xp[0]+Spline.Xp[n]+Spline.Xp[n+1])*
            Triangle_Area(Spline.Xp[0],Spline.Xp[n],Spline.Xp[n+1])/THREE;
  Xref = Xref/Area();
}

/**********************************************************************
 * Interface2D::BoundingBox -- Determine the bounding box of the      *
 *                             interface.                             *
 **********************************************************************/
void Interface2D::BoundingBox(void) {
  // Assert that the interface/spline has at least two points.
  assert(Spline.np > 1);
  Vector2D Xmin_temp, Xmax_temp;
  // Initialize Xmin and Xmax
  Xmin = Vector2D(MILLION,MILLION);
  Xmax = -Vector2D(MILLION,MILLION);
  // Determine the bounding box.
  for (int n = 0; n < Spline.np; n++) {
    if (Spline.Xp[n].x < Xmin.x) Xmin.x = Spline.Xp[n].x;
    if (Spline.Xp[n].y < Xmin.y) Xmin.y = Spline.Xp[n].y;
    if (Spline.Xp[n].x > Xmax.x) Xmax.x = Spline.Xp[n].x;
    if (Spline.Xp[n].y > Xmax.y) Xmax.y = Spline.Xp[n].y;
  }
  // Special case: If the area of the bounding box is zero then include
  // twice the reference point in the bounding-box creation.
  if ((Xmax.x-Xmin.x)*(Xmax.y-Xmin.y) < NANO) {
    if (Xref.y > Xmax.y) Xmax.y += TWO*(Xref.y-Xmax.y);
    if (Xref.y < Xmin.y) Xmin.y += TWO*(Xref.y-Xmin.y);
    if (Xref.x > Xmax.x) Xmax.x += TWO*(Xref.x-Xmax.x);
    if (Xref.x < Xmin.x) Xmin.x += TWO*(Xref.x-Xmin.x);
  }
}

/**********************************************************************
 * Interface2D::BoundingBox -- Determine an extended bounding box of  *
 *                             the interface.                         *
 **********************************************************************/
void Interface2D::BoundingBox(const double &val) {
  // Determine bounding box.
  BoundingBox();
  // Extend the bounding box by the specified amount.
  Xmin -= Vector2D(val,val);
  Xmax += Vector2D(val,val);
}

/**********************************************************************
 * Interface2D::Point_In_Bounding_Box -- Determine if the given point *
 *                                       is contained in the bounding *
 *                                       box of the interface.        *
 **********************************************************************/
int Interface2D::Point_In_Bounding_Box(const Vector2D &X) {
  if (X.x >= Xmin.x && X.y >= Xmin.y && X.x <= Xmax.x && X.y <= Xmax.y) return 1;
//   if (fabs(X.x-Xmin.x) < TOLER*TOLER) return 1;
//   if (fabs(X.y-Xmin.y) < TOLER*TOLER) return 1;
//   if (fabs(X.x-Xmax.x) < TOLER*TOLER) return 1;
//   if (fabs(X.y-Xmax.y) < TOLER*TOLER) return 1;
// //   if (abs(X-Xmin) < NANO) return 1;
// //   if (abs(X-Xmax) < NANO) return 1;
// //   if (fabs(X.x - Xmin.x) < NANO && fabs(X.y - Xmin.y) < NANO &&
// //       fabs(X.x - Xmax.x) < NANO && fabs(X.y - Xmax.y) < NANO) return 1;
//   if ((X.x >= Xmin.x || fabs(X.x-Xmin.x) < NANO) &&
//       (X.y >= Xmin.y || fabs(X.y-Xmin.y) < NANO) &&
//       (X.x <= Xmax.x || fabs(X.x-Xmax.x) < NANO) &&
//       (X.y <= Xmax.y || fabs(X.y-Xmax.y) < NANO)) return 1;
  return 0;
}

/**********************************************************************
 * Interface2D::Bounding_Box_Intersection -- Determine if the given   *
 *                                           bounding box intersects  *
 *                                           with the current         *
 *                                           bounding box.            *
 **********************************************************************/
int Interface2D::Bounding_Box_Intersection(const Vector2D &Ymin, const Vector2D &Ymax) {
  if (Point_In_Bounding_Box(Ymin) || Point_In_Bounding_Box(Ymax) ||
      Point_In_Bounding_Box(Vector2D(Ymin.x,Ymax.y)) ||
      Point_In_Bounding_Box(Vector2D(Ymax.x,Ymin.y))) { return 1; }
  if (Line_Intersection(Vector2D(Xmin.x,Xmin.y),Vector2D(Xmin.x,Xmax.y),
			Vector2D(Ymin.x,Ymin.y),Vector2D(Ymax.x,Ymin.y))) { return 1; }
  if (Line_Intersection(Vector2D(Xmax.x,Xmin.y),Vector2D(Xmax.x,Xmax.y),
			Vector2D(Ymin.x,Ymax.y),Vector2D(Ymax.x,Ymax.y))) { return 1; }
  if (Line_Intersection(Vector2D(Xmin.x,Xmin.y),Vector2D(Xmax.x,Xmin.y),
			Vector2D(Ymin.x,Ymin.y),Vector2D(Ymin.x,Ymax.y))) { return 1; }
  if (Line_Intersection(Vector2D(Xmin.x,Xmax.y),Vector2D(Xmax.x,Xmax.y),
			Vector2D(Ymax.x,Ymin.y),Vector2D(Ymax.x,Ymax.y))) { return 1; }
  // Return no intersection.
  return 0;
}

/**********************************************************************
 * Interface2D::Point_In_Interface -- Determine if the given point is *
 *                                    inside the interface.           *
 *                                                                    *
 * This routine returns an integer value (0 or 1) indicating if a     *
 * given point, X, is contained within the interface.  A ray-tracing  *
 * algorithm is used with a bounding box filter.  The spline is       *
 * assumed to be linear.                                              *
 **********************************************************************/
int Interface2D::Point_In_Interface(const Vector2D &Xt) {
  // Conduct a quick check by using the bounding box.
  if (!Point_In_Bounding_Box(Xt)) return 0;
  // If the point of interest is the same as the reference point then
  // immediately exit successfully.
  if (abs(Xref-Xt) < NANO) return 1;
  // If the point is within the bounding box filter then use ray-
  // tracing to determine the location of the given point.
  int number_of_intersections = 0;
  Vector2D Xp;
  // Sample each edge of the interface (assumed linear) and count 
  // the number of intersection points with the line defined by 
  // the interface reference point and the point of interest.
  for (int n = 0; n < Spline.np-1; n++) {
    // If the test point is coincident to one of the interface points
    // then immediately return success.
    //if (Xt == Spline.Xp[n] || Xt == Spline.Xp[n+1]) return 2;
    //if (Point_On_Line(Spline.Xp[n],Spline.Xp[n+1],Xt)) return 2;
    // Get the points defining the edge of the polygon.
    if (Line_Intersection(Spline.Xp[n],Spline.Xp[n+1],Xt,Xref,Xp,NANO*abs(Spline.Xp[n]-Spline.Xp[n+1]))) {
      if (abs(Xp-Spline.Xp[n]) > NANO*abs(Spline.Xp[n]-Spline.Xp[n+1])) number_of_intersections++;
    }
  }
  // Increment number of intersections counter such that an odd
  // number of intersections returns a 0.
  number_of_intersections++;
  // Determine and return the result for the test point, Xt.
  return number_of_intersections%2;
}

/**********************************************************************
 * Interface2D::Zalesak -- This routine calculates and returns an     *
 *                         interface representing Zalesak's disk.     *
 **********************************************************************/
void Interface2D::Zalesak(const Vector2D &Origin,
			  const double &Radius) {

  int npts;
  double theta, theta1, theta2;

  // Set the number of spline points.
  npts = 341 + 3;

  // Set the angles.
  theta1 = 10.0;
  theta2 = 350.0;

  // Allocate memory for the circular arc spline.
  if (Spline.np != 0) Spline.deallocate();
  Spline.allocate(npts);

  // Set the spline type.
  Spline.settype(SPLINE2D_LINEAR);

  // Compute the locations of the spline points on the 
  // spline and the point and boundary condition types.
  for (int i = 0; i < npts-3; i++) {
    theta = theta1 + (theta2 - theta1)*double(i)/double(npts-4);
    theta = TWO*PI*theta/360.0;
    Spline.Xp[i].x = Radius*cos(theta);//Origin.x + 
    Spline.Xp[i].y = Radius*sin(theta);//Origin.y + 
    Spline.bc[i]   = BC_NONE;
    Spline.tp[i]   = SPLINE2D_POINT_NORMAL;
    if (i == 0 || i == npts - 4)
      Spline.tp[i] = SPLINE2D_POINT_SHARP_CORNER;
  }
  Spline.Xp[npts-3] = Spline.Xp[npts-4] - Vector2D(1.5*Radius,ZERO);
  Spline.bc[npts-3] = BC_NONE;
  Spline.tp[npts-3] = SPLINE2D_POINT_SHARP_CORNER;
  Spline.Xp[npts-2] = Spline.Xp[0] - Vector2D(1.5*Radius,ZERO);
  Spline.bc[npts-2] = BC_NONE;
  Spline.tp[npts-2] = SPLINE2D_POINT_SHARP_CORNER;

  // Last point must be the same as the first point.
  Spline.Xp[npts-1] = Spline.Xp[0];
  Spline.bc[npts-1] = Spline.bc[0];
  Spline.tp[npts-1] = Spline.tp[0];

  // Calculate the spline pathlengths.
  Spline.pathlength();

  // Rotate the spline.
  Rotate_Spline(Spline,-HALF*PI);

  // Translate the spline.
  Translate_Spline(Spline,Origin);

}

/**********************************************************************
 * Interface2D::Ringleb -- This routine calculates and returns an     *
 *                         interface representing Ringleb's flow      *
 *                         domain.                                    *
 **********************************************************************/
void Interface2D::Ringleb(const double &Inner_Streamline_Number,
			  const double &Outer_Streamline_Number,
			  const double &Isotach_Line) {

  assert(Inner_Streamline_Number > Outer_Streamline_Number);
  assert(Inner_Streamline_Number < 5.0/3.0);
  assert(Outer_Streamline_Number > Isotach_Line); 

  int nk, nq;
  double **rho;
  double delta_k, delta_q;
  double k_init, q_init, q_final;
  double **q, **k, qo, ko, c, J;
  double g = 1.40;
  Spline2D Spline_North, Spline_East, Spline_West, Spline_Reflection;

  // Allocate memory.
  nk = 32;
  nq = 50;
  k = new double*[nk];
  q = new double*[nk];
  rho = new double*[nk];
  for (int i = 0; i < nk; i++) {
    k[i] = new double[nq];
    q[i] = new double[nq];
    rho[i] = new double[nq];
    for (int j = 0; j < nq; j++) {
      k[i][j] = ZERO;
      q[i][j] = ZERO;
      rho[i][j] = ZERO;
    }
  }

  // Determine q, k.
  delta_k = (Inner_Streamline_Number-Outer_Streamline_Number)/double(nk-1);
  k_init = Outer_Streamline_Number;
  q_init = Isotach_Line;
  for (int i = 0; i < nk; i++) {
    ko = k_init + double(i)*delta_k;
    q_final = ko; // condition y = 0
    delta_q = (q_final - q_init)/double(nq-1);
    for (int j = 0; j < nq; j++) {
      qo = q_init + double(j)*delta_q;
      k[i][j] = ko;
      q[i][j] = qo;
    }
  }

  // Allocate memory for the interface splines.
  Spline_North.allocate(nk); Spline_North.settype(SPLINE2D_QUINTIC);
  Spline_East.allocate(nq);  Spline_East.settype(SPLINE2D_QUINTIC);
  Spline_West.allocate(nq);  Spline_West.settype(SPLINE2D_QUINTIC);

  // Determine the spline points.
  for (int i = 0; i < nk; i++) {
    for (int j = 0; j < nq; j++) {

      c = sqrt(ONE - ((g-ONE)/TWO)*q[i][j]*q[i][j]);
      rho[i][j] = pow(c,TWO/(g-ONE));
      J = ONE/c + ONE/(THREE*pow(c,THREE)) + ONE/(FIVE*pow(c,FIVE)) - HALF*log((ONE+c)/(ONE-c));
      // NORTH spline.
      if (j == 0) {
	Spline_North.Xp[nk-1-i].x = (HALF/rho[i][j])*(TWO/(k[i][j]*k[i][j]) - ONE/(q[i][j]*q[i][j])) - HALF*J;
	Spline_North.Xp[nk-1-i].y = (ONE/(k[i][j]*rho[i][j]*q[i][j]))*sqrt(ONE - (q[i][j]*q[i][j])/(k[i][j]*k[i][j]));
	Spline_North.bc[nk-1-i] = BC_RINGLEB_FLOW;
	Spline_North.tp[nk-1-i] = SPLINE2D_POINT_NORMAL;
      }
      // EAST spline.
      if (i == 0) {
	Spline_East.Xp[j].x = (HALF/rho[i][j])*(TWO/(k[i][j]*k[i][j]) - ONE/(q[i][j]*q[i][j])) - HALF*J;
	Spline_East.Xp[j].y = (ONE/(k[i][j]*rho[i][j]*q[i][j]))*sqrt(ONE - (q[i][j]*q[i][j])/(k[i][j]*k[i][j]));
	Spline_East.bc[j] = BC_RINGLEB_FLOW;
	Spline_East.tp[j] = SPLINE2D_POINT_NORMAL;
      }
      // WEST spline.
      if (i == nk-1) {
	Spline_West.Xp[nq-1-j].x = (HALF/rho[i][j])*(TWO/(k[i][j]*k[i][j]) - ONE/(q[i][j]*q[i][j])) - HALF*J;
	Spline_West.Xp[nq-1-j].y = (ONE/(k[i][j]*rho[i][j]*q[i][j]))*sqrt(ONE - (q[i][j]*q[i][j])/(k[i][j]*k[i][j]));
	Spline_West.bc[nq-1-j] = BC_RINGLEB_FLOW;
	Spline_West.tp[nq-1-j] = SPLINE2D_POINT_NORMAL;
      }

    }
  }

  // Calculate spline path lengths.
  Spline_North.pathlength();
  Spline_East.pathlength();
  Spline_West.pathlength();

  // Set the boundary condition types for each of the interface splines.
  Spline_North.setBCtype(BC_RINGLEB_FLOW);
  Spline_East.setBCtype(BC_RINGLEB_FLOW);
  Spline_West.setBCtype(BC_RINGLEB_FLOW);

  // Reflect Spline_West and concatenate to Spline_West.
  Copy_Spline(Spline_Reflection,Spline_West);
  Reflect_Spline(Spline_Reflection);
  Reverse_Spline(Spline_Reflection);
  Spline_Reflection = Concatenate_Splines(Spline_Reflection,Spline_West);
  Copy_Spline(Spline_West,Spline_Reflection);
  Spline_West.settptype(SPLINE2D_POINT_NORMAL);

  // Reflect Spline_East and concatenate to Spline_East.
  Copy_Spline(Spline_Reflection,Spline_East);
  Reflect_Spline(Spline_Reflection);
  Reverse_Spline(Spline_Reflection);
  Spline_East = Concatenate_Splines(Spline_East,Spline_Reflection);
  Spline_East.settptype(SPLINE2D_POINT_NORMAL);

  // Concatenate sub-splines to the Ringleb spline.
  Spline = Concatenate_Splines(Spline_West,Spline_North);
  Spline = Concatenate_Splines(Spline,Spline_East);

  // Deallocate the memory for the interface splines.
  Spline_North.deallocate();
  Spline_East.deallocate();
  Spline_West.deallocate();

  // Deallocate memory for k, q, and rho.
  for (int i = 0; i < nk; i++) {
    delete []k[i];   k[i]   = NULL;
    delete []q[i];   q[i]   = NULL;
    delete []rho[i]; rho[i] = NULL;
  }
  delete []k;   k   = NULL;
  delete []q;   q   = NULL;
  delete []rho; rho = NULL;

}

/**********************************************************************
 * Interface2D::Bump_Channel_Flow -- This routine creates the bump    *
 *                                   channel flow interface spline.   *
 **********************************************************************/
void Interface2D::Bump_Channel_Flow(void) {

  Spline2D Leading_Spline, Bump_Spline, Trailing_Spline;
  Vector2D X1, X2;
  double R, Theta;

  // Create the leading line spline.
  X1 = Vector2D(-TWO,ZERO);
  X2 = Vector2D(ZERO,ZERO);
  Create_Spline_Line(Leading_Spline,X1,X2,2);

  // Create the trailing line spline.
  X1 = Vector2D(ONE,ZERO);
  X2 = Vector2D(FIVE,ZERO);
  Create_Spline_Line(Trailing_Spline,X1,X2,2);

  // Create the bump spline.
  R = (0.25 + 0.042*0.042)/(TWO*0.042);
  Theta = acos((R - 0.042)/R);
  Create_Spline_Circular_Arc(Bump_Spline,
			     Vector2D_ZERO,
			     R,
			     ZERO,
			     -TWO*Theta*180.0/PI,
			     31);
  Rotate_Spline(Bump_Spline,HALF*PI+Theta);
  Translate_Spline(Bump_Spline,Vector2D(HALF,0.042-R));

  // Concatenate splines.
  Copy_Spline(Spline,Leading_Spline);
  Spline = Concatenate_Splines(Spline,Bump_Spline);
  Spline = Concatenate_Splines(Spline,Trailing_Spline);

  // Deallocate splines.
  Leading_Spline.deallocate();
  Trailing_Spline.deallocate();
  Bump_Spline.deallocate();

}

/**********************************************************************
 * Interface2D::Flat_Plate -- This routine creates the flat plate     *
 *                            interface spline.                       *
 **********************************************************************/
void Interface2D::Flat_Plate(void) {
  // Allocate spline.
  Spline.allocate(3);
  // Set spline type.
  Spline.settype(SPLINE2D_LINEAR);
  // Set spline nodes.
  Spline.Xp[0] = Vector2D(-TWO*Length2*Length1,ZERO);
  Spline.Xp[1] = Vector2D(ZERO,ZERO);
  Spline.Xp[2] = Vector2D(TWO*Length2*Length1,ZERO);
  // Set spline node boundary conditions.
  Spline.bc[0] = BC_REFLECTION;
  Spline.bc[1] = BC_WALL_VISCOUS_HEATFLUX;
  Spline.bc[2] = BC_WALL_VISCOUS_HEATFLUX;
  // Set spline node types.
  Spline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
  Spline.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
  Spline.tp[2] = SPLINE2D_POINT_SHARP_CORNER;
  // Calculate the spline pathlengths.
  Spline.pathlength();
  // Rotate spline.
  //Rotate_Spline(Spline,Length*PI/180.0);
}

/**********************************************************************
 * Interface2D::Determine_Interface_BC_Type -- Return the boundary    *
 *                                             condition defined on   *
 *                                             the interface at the   *
 *                                             given location.        *
 **********************************************************************/
int Interface2D::Determine_Interface_BC_Type(const Vector2D &X) {

  int bctype;

  // Determine the boundary condition.
  switch(BC_Type){
  case INTERFACE_BC_REFLECTION :
  case INTERFACE_BC_ABSORPTION :
  case INTERFACE_BC_BURNING_SURFACE :
  case INTERFACE_BC_WALL_VISCOUS_HEATFLUX :
  case INTERFACE_BC_WALL_VISCOUS_ISOTHERMAL :
    bctype = BC_Type;
    break;
  case INTERFACE_BC_FLAT_PLATE :
    if (X.x < ZERO) bctype = INTERFACE_BC_REFLECTION;
    else bctype = INTERFACE_BC_WALL_VISCOUS_HEATFLUX;
    break;
  case INTERFACE_BC_RINGLEB :
  case INTERFACE_BC_DETERMINE :
  case INTERFACE_BC_MIXED :
    bctype = getBCtype(X,Spline);
    break;
  };

  // Return the boundary condition.
  return bctype;

}

/**********************************************************************
 * Interface2D::Determine_Interface_Velocity -- Returns the velocity  *
 *                                              defined on the        *
 *                                              interface at the      *
 *                                              given location.       *
 **********************************************************************/
Vector2D Interface2D::Determine_Interface_Velocity(const Vector2D &X,
						   const double &time) {

  int np, subinterval_found;
  Vector2D W, W1, W2, X1, X2;

  switch(Motion) {
  case INTERFACE_MOTION_STATIONARY :
  case INTERFACE_MOTION_LEVELSET_STATIONARY :
    // The interface is stationary.
    W = Vector2D_ZERO;
    break;
  case INTERFACE_MOTION_CONSTANT :
  case INTERFACE_MOTION_EXPAND :
    // The interface has a constant speed in the radial direction.
    W = Speed.x*((X - Xref)/abs(X - Xref));
    break;
  case INTERFACE_MOTION_UNIFORM :
  case INTERFACE_MOTION_TRANSLATE :
    // The interface has a constant velocity.
    W = Speed;
    break;
  case INTERFACE_MOTION_ROTATE :
    // The interface has a constant rotational velocity.
    switch(Type) {
    case INTERFACE_CIRCLE :
      W = Vector2D((TWO*PI*Speed.x/Speed.y)*cos(TWO*PI*time/Speed.y),ZERO);
      break;
    case INTERFACE_ELLIPSE :
      break;
    case INTERFACE_SQUARE :
    case INTERFACE_USER_SPECIFIED :
    case INTERFACE_ZALESAK :
      W = Vector2D(-Speed.x*(X.y-Xref.y),Speed.x*(X.x-Xref.x));
      break;
    case INTERFACE_NACA0012_AEROFOIL :
//       W.x =  (5.02/(2.0*PI*62.5))*cos(2.0*PI*62.5*(time))*(X.y-Xref.y);
//       W.y = -(5.02/(2.0*PI*62.5))*cos(2.0*PI*62.5*(time))*(X.x-Xref.x);
      W.x =  (2.0*PI*62.5*(PI*2.51/180.0))*cos(2.0*PI*62.5*(time))*(X.y-Xref.y);
      W.y = -(2.0*PI*62.5*(PI*2.51/180.0))*cos(2.0*PI*62.5*(time))*(X.x-Xref.x);
      break;
    case INTERFACE_NACA0015_AEROFOIL :
      W.x =  (2.0*PI*62.5*(PI*2.51/180.0))*cos(2.0*PI*62.5*(time))*(X.y-Xref.y);
      W.y = -(2.0*PI*62.5*(PI*2.51/180.0))*cos(2.0*PI*62.5*(time))*(X.x-Xref.x);
      break;
    };
    break;
  case INTERFACE_MOTION_STRETCH :
  case INTERFACE_MOTION_LEVELSET :
  case INTERFACE_MOTION_LEVELSET_EXPAND :
  case INTERFACE_MOTION_LEVELSET_STRETCH :
  case INTERFACE_MOTION_BURNING_SURFACE :
  case INTERFACE_MOTION_MOMENTUM_TRANSFER :
  case INTERFACE_MOTION_LEVELSET_BULKFLOW :
    // Use linear interpolation on the subinterval containing X to
    // determine the velocity function at the location X on the 
    // interface spline.
    subinterval_found = 0;
    np = 0;
    while (!subinterval_found && np < Spline.np-1) {
      X1 = Spline.Xp[np];
      X2 = Spline.Xp[np+1];
      if (((X.x-X1.x)*(X.x-X2.x) <= ZERO) && ((X.y-X1.y)*(X.y-X2.y) <= ZERO)) {
	subinterval_found = 1;
	W1 = F[np];
	W2 = F[np+1];
      } else {
 	np++;
      }
    }
    if (!subinterval_found) {
      subinterval_found = 0;
      np = 0;
      while (!subinterval_found && np < Spline.np-2) {
 	X1 = Spline.Xp[np];
 	X2 = Spline.Xp[np+2];
	if (((X.x-X1.x)*(X.x-X2.x) <= ZERO) && ((X.y-X1.y)*(X.y-X2.y) <= ZERO)) {
 	  subinterval_found = 1;
	  W1 = F[np];
	  W2 = F[np+2];
 	} else {
 	  np++;
 	}
      }
      if (!subinterval_found) {
 	subinterval_found = 0;
 	np = 0;
	for (int n = 0; n < Spline.np; n++)
	  if (abs(Spline.Xp[n]-X) < abs(Spline.Xp[np]-X)) np = n;
	W1 = F[np];
	W2 = F[np];
      }
    }
    W = Linear_Interpolation(X1,W1,X2,W2,X);
    break;
  };

  // Return the interface velocity.
  return W;

}

/**********************************************************************
 * Interface2D::Fn -- Return normal velocity at interface point np.   *
 **********************************************************************/
// double Interface2D::Fn(const int &np) {
double Interface2D::Fn(const int &np, const Vector2D &nhat) const {
  return dot(F[np],nhat);
  //return dot(F[np],F[np])/F[np].abs();
  //return F[np].x;
}

double Interface2D::Fn(const int &np) const {
  return dot(F[np],normal(np));
  //return dot(F[np],F[np])/F[np].abs();
  //return F[np].x;
}

void Interface2D::Fn(const int &np, const double &fn) {
//   Vector2D norm_dir;
//   if (np == Spline.np-1) {
//     if (abs(Spline.Xp[np]-Spline.Xp[0]) < NANO) norm_dir = nhat(Spline.Xp[1],Spline.Xp[np]);
//     else norm_dir = nhat(Spline.Xp[np],Spline.Xp[np-1]);
//   } else {
//     norm_dir = nhat(Spline.Xp[np+1],Spline.Xp[np]);
//   }
//   F[np] = fn*norm_dir;
//   F[np].x = fn;
//   F[np].y = ZERO;
  F[np] = -fn*normal(np);
}

/**********************************************************************
 * Interface2D::normal --                                             *
**********************************************************************/
Vector2D Interface2D::normal(const int &np) const {
  Vector2D norm_dir;
  if (np == 0) {
    if (abs(Spline.Xp[Spline.np-1]-Spline.Xp[np]) < NANO) norm_dir = nhat(Spline.Xp[np+1],Spline.Xp[Spline.np-2]);
    else norm_dir = nhat(Spline.Xp[np+1],Spline.Xp[np]);
  } else if (np == Spline.np-1) {
    if (abs(Spline.Xp[np]-Spline.Xp[0]) < NANO) norm_dir = nhat(Spline.Xp[1],Spline.Xp[Spline.np-2]);
    else norm_dir = nhat(Spline.Xp[np],Spline.Xp[np-1]);
  } else {
    norm_dir = nhat(Spline.Xp[np+1],Spline.Xp[np-1]);
  }
  norm_dir *= -ONE;
  return norm_dir;
}

/**********************************************************************
 * Interface2D::max_speed -- Return the maximum speed defined on the  *
 *                           interface.                               *
 **********************************************************************/
double Interface2D::max_speed(void) {
  double vmax = ZERO;
  for (int n = 0; n < Spline.np; n++) vmax = max(vmax,abs(F[n]));
  return vmax;
}

/**********************************************************************
 * Interface2D::Determine_Interface_Temperature -- Return the         *
 *                                                 temperature        *
 *                                                 defined on the     *
 *                                                 interface at the   *
 *                                                 given location.    *
 **********************************************************************/
double Interface2D::Determine_Interface_Temperature(const Vector2D &X) {

  // Return the interface temperature.
  return 288.16;

}

/**********************************************************************
 * Interface2D::Interface_Union --                                    *
 *                                                                    *
 * Return the interface of union between the current and given        *
 * interfaces.  The scheme used here is based on the Weiler-Atherton  *
 * algorithm popular in computer graphics.                            *
 *                                                                    *
 **********************************************************************/
int Interface2D::Interface_Union(Interface2D &I2) {

  Interface2D I1;
  int n1, n2, type, ft11, ft12, ft21, ft22;
  int tp11, tp12, tp21, tp22, bc11, bc12, bc21, bc22;
  Vector2D X11, X12, X21, X22, Xp, *X, f11, f12, f21, f22;
  LinkedList<Vector2D> Y, Y1, Y2, f, f1, f2;
  LinkedList<int> T1, T2, tp, tp1, tp2, bc, bc1, bc2, ft, ft1, ft2;
  int internal_flag_1, internal_flag_2, intersect_flag, first_intersection_position;
  int swap_counter;

  // Store current interface and deallocate.
  I1 = *this;
  deallocate();

  // If all the points of interface 1 are internal to interface 2 then
  // interface 2 is the interface of union.  If all of the points of
  // interface 1 are external to interface 2 then interface 1 is the 
  // interface of union.
  internal_flag_1 = 0;
  for (int n1 = 0; n1 < I1.Spline.np; n1++) {
    internal_flag_1 += I2.Point_In_Interface(I1.Spline.Xp[n1]);
  }
  internal_flag_2 = 0;
  for (int n2 = 0; n2 < I2.Spline.np; n2++) {
    internal_flag_2 += I1.Point_In_Interface(I2.Spline.Xp[n2]);
  }
  if (internal_flag_1 == I1.Spline.np && internal_flag_2 == 0) {
    Copy(I2);
  } else if (internal_flag_1 == 0 && internal_flag_2 == I2.Spline.np) {
    Copy(I1);
  } else {

    // Determine the polygon of union using the Weiler-Atherton algorithm.

    // Include the intersection points in interface list 1.
    first_intersection_position = 0;
    for (int n1 = 0; n1 < I1.Spline.np-1; n1++) {
      // Get the points defining the edge of interface 1..
      X11  = I1.Spline.Xp[n1];  X12  = I1.Spline.Xp[n1+1];
      tp11 = I1.Spline.tp[n1];  tp12 = I1.Spline.tp[n1+1];
      bc11 = I1.Spline.bc[n1];  bc12 = I1.Spline.bc[n1+1];
      ft11 = I1.F_Type[n1];     ft12 = I1.F_Type[n1+1];
      f11  = I1.F[n1];          f12  = I1.F[n1+1];
      // Add point X11 to interface list 1.
      Y1.add(X11); T1.add(0);
      tp1.add(tp11); bc1.add(bc11);
      ft1.add(ft11); f1.add(f11);
      for (int n2 = 0; n2 < I2.Spline.np; n2++) {
	if (abs(X11 - I2.Spline.Xp[n2]) < NANO) {
	  T1.put(1);
	  if (first_intersection_position == -1)
	    first_intersection_position = Y1.np-1;
	}
      }
      // Cycle through the edges of interface 2 with X11 and X12.
      swap_counter = 0;
      for (int n2 = 0; n2 < I2.Spline.np-1; n2++) {
	// Get the points defining the edge of interface 2.
	X21  = I2.Spline.Xp[n2];  X22  = I2.Spline.Xp[n2+1];
	tp21 = I2.Spline.tp[n2];  tp22 = I2.Spline.tp[n2+1];
	bc21 = I2.Spline.tp[n2];  bc22 = I2.Spline.bc[n2+1];
	ft21 = I2.F_Type[n2];     ft22 = I2.F_Type[n2+1];
	f21  = I2.F[n2];          f22  = I2.F[n2+1];
	// Determine and add the intersection point if any.
	intersect_flag = Line_Intersection(X11,X12,X21,X22,Xp);
	if (intersect_flag && Y1.find(Xp,NANO) == NULL && abs(X12-Xp) > NANO) {
	  if (first_intersection_position == -1)
	    first_intersection_position = Y1.np;
	  Y1.add(Xp);
	  T1.add(1);
	  tp1.add(SPLINE2D_POINT_SHARP_CORNER);
	  if (bc11 >= bc12 && bc11 >= bc21 && bc11 >= bc22) {
	    bc1.add(bc11); ft1.add(ft11); f1.add(f11);
	  } else if (bc12 >= bc11 && bc12 >= bc21 && bc12 >= bc22) {
	    bc1.add(bc12); ft1.add(ft11); f1.add(f12);
	  } else if (bc21 >= bc11 && bc21 >= bc12 && bc21 >= bc22) {
	    bc1.add(bc21); ft1.add(ft21); f1.add(f21);
	  } else if (bc22 >= bc11 && bc22 >= bc12 && bc22 >= bc21) {
	    bc1.add(bc22); ft1.add(ft22); f1.add(f22);
	  }
	  swap_counter++;
	} else if (abs(X11-Xp) < NANO) {
	  if (first_intersection_position == -1)
	    first_intersection_position = Y1.np-1;
	  T1.put(1);
	}
      }
      // Swap intersection points into counter-clockwise order if necessary.
      if (swap_counter) {
	for (int current = Y1.np-swap_counter; current < Y1.np; current++) {
	  for (int next = current+1; next < Y1.np; next++) {
	    if (abs(Y1[next]-X11) < abs(Y1[current]-X11)) {
	      Y1.swap(current,next); T1.swap(current,next);
	      tp1.swap(current,next); bc1.swap(current,next);
	      ft1.swap(current,next); f1.swap(current,next);
	    }
	  }
	}
      }
    }

    // Include the intersection points in interface list 2.
    for (int n2 = 0; n2 < I2.Spline.np-1; n2++) {
      // Get the points defining the edge of interface 2.
      X21  = I2.Spline.Xp[n2];  X22  = I2.Spline.Xp[n2+1];
      tp21 = I2.Spline.tp[n2];  tp22 = I2.Spline.tp[n2+1];
      bc21 = I2.Spline.bc[n2];  bc22 = I2.Spline.bc[n2+1];
      ft21 = I2.F_Type[n2];     ft22 = I2.F_Type[n2+1];
      f21  = I2.F[n2];          f22  = I2.F[n2+1];
      // Add point X21 to polygon list 2.
      Y2.add(X21); T2.add(0);
      tp2.add(tp21); bc2.add(bc21);
      ft2.add(ft21); f2.add(f21);
      for (int n1 = 0; n1 < I1.Spline.np; n1++) {
	if (abs(X21 - I1.Spline.Xp[n1]) < NANO) T2.put(1);
      }
      // Cycle through the edges of interface 1 with X21 and X22.
      swap_counter = 0;
      for (int n1 = 0; n1 < I1.Spline.np-1; n1++) {
	// Get the points defining the edge of interface 1.
	X11  = I1.Spline.Xp[n1];  X12  = I1.Spline.Xp[n1+1];
	tp11 = I1.Spline.tp[n1];  tp12 = I1.Spline.tp[n1+1];
	bc11 = I1.Spline.bc[n1];  bc12 = I1.Spline.bc[n1+1];
	ft11 = I1.F_Type[n1];     ft12 = I1.F_Type[n1+1];
	f11  = I1.F[n1];          f12  = I1.F[n1+1];
	// Determine and add the intersection point if any.
	intersect_flag = Line_Intersection(X21,X22,X11,X12,Xp);
	if (intersect_flag && Y2.find(Xp,NANO) == NULL && abs(X22-Xp) > NANO) {
	  Y2.add(Xp);
	  T2.add(1);
	  tp2.add(SPLINE2D_POINT_SHARP_CORNER);
	  if (bc11 >= bc12 && bc11 >= bc21 && bc11 >= bc22) {
	    bc2.add(bc11); ft2.add(ft11); f2.add(f11);
	  } else if (bc12 >= bc11 && bc12 >= bc21 && bc12 >= bc22) {
	    bc2.add(bc12); ft2.add(ft12); f2.add(f12);
	  } else if (bc21 >= bc11 && bc21 >= bc12 && bc21 >= bc22) {
	    bc2.add(bc21); ft2.add(ft21); f2.add(f21);
	  } else if (bc22 >= bc11 && bc22 >= bc12 && bc22 >= bc21) {
	    bc2.add(bc22); ft2.add(ft22); f2.add(f22);
	  }
	  swap_counter++;
	} else if (abs(X21-Xp) < NANO) {
	  T2.put(1);
	}
      }
      // Swap intersection points into counter-clockwise order if necessary.
      if (swap_counter) {
	for (int current = Y2.np-swap_counter; current < Y2.np; current++) {
	  for (int next = current+1; next < Y2.np; next++) {
	    if (abs(Y2[next]-X21) < abs(Y2[current]-X21)) {
	      Y2.swap(current,next); T2.swap(current,next);
	      tp2.swap(current,next); bc2.swap(current,next);
	      ft2.swap(current,next); f2.swap(current,next);
	    }
	  }
	}
      }
    }

    // Reorder interface list 2 such that the first intersection point in
    // list one is the first point in list two.
    first_intersection_position = Y2.find_position(Y1[first_intersection_position]);
    if (first_intersection_position == -1) return 1;//cout << "ERROR in interface union code." << endl;
    Y2.shift_order(first_intersection_position);
    T2.shift_order(first_intersection_position);
    tp2.shift_order(first_intersection_position);
    bc2.shift_order(first_intersection_position);
    ft2.shift_order(first_intersection_position);
    f2.shift_order(first_intersection_position);

    // Create the interface of union between interfaces 1 and 2.
    int n1, n2;
    n1 = 0; n2 = 0;
    while (n1 < Y1.np) {
      // Sample polygon point list 1.
      if (T1[n1] == 0) {
	// If the current point in list 1 is external to polygon 2 then
	// add the point to the list.
	if (!I2.Point_In_Interface(Y1[n1])) {
	  Y.add(Y1[n1]);
	  tp.add(tp1[n1]); bc.add(bc1[n1]);
	  ft.add(ft1[n1]); f.add(f1[n1]);
	}
	n1++;
      } else if (T1[n1] == 1) {
	// If the current point in list 1 is an intersection point then
	// add the point to the list.
	Y.add(Y1[n1]);
	tp.add(tp1[n1]); bc.add(bc1[n1]);
	ft.add(ft1[n1]); f.add(f1[n1]);
	// Add all points in list 1 up until the next intersection point
	// if they are external to polygon 2.
	while (1) {
	  n1++;
	  if (n1 >= Y1.np) break;
	  if (T1[n1] == 0 && !I2.Point_In_Interface(Y1[n1])) {
	    Y.add(Y1[n1]);
	    tp.add(tp1[n1]); bc.add(bc1[n1]);
	    ft.add(ft1[n1]); f.add(f1[n1]);
	  }
	  if (T1[n1] == 1) break;
	}
	// Add all points in list 2 up until the next intersection point
	// if they are external to polygon 1.
	if (n2 < Y2.np) {
	  while (1) {
	    n2++;
	    if (n2 >= Y2.np) break;
	    if (T2[n2] == 0 && !I1.Point_In_Interface(Y2[n2])) {
	      Y.add(Y2[n2]);
	      tp.add(tp2[n2]); bc.add(bc2[n2]);
	      ft.add(ft2[n2]); f.add(f2[n2]);
	    }
	    if (T2[n2] == 1) break;
	  }
	}
      }
    }
    // Add the rest of the Y2 points if necessary.
    while (n2 < Y2.np) {
      if (T2[n2] == 0 && !I1.Point_In_Interface(Y2[n2])) {
	Y.add(Y2[n2]);
	tp.add(tp2[n2]); bc.add(bc2[n2]);
	ft.add(ft2[n2]); f.add(f2[n2]);
      }
      n2++;
    }

    // Create interface of union from the linked lists of interface points, 
    // interface point type, boundary condition type, and velocity function.
    Spline.settype(min(I1.Spline.type,I2.Spline.type));
    if (abs(I1.Spline.Xp[I1.Spline.np-1] - I1.Spline.Xp[0]) < NANO &&
	abs(I2.Spline.Xp[I2.Spline.np-1] - I2.Spline.Xp[0]) < NANO) {
      Spline.allocate(Y.np+1);
    } else {
      Spline.allocate(Y.np);
    }
    Initialize_Velocity_Function();
    for (int np = 0; np < Y.np; np++) {
      Spline.Xp[np] = Y[np];
      Spline.tp[np] = tp[np];
      Spline.bc[np] = bc[np];
      F_Type[np] = ft[np];
      F[np] = f[np];
    }
    // Close interface of union if required.
    if (Spline.np == Y.np+1) {
      Spline.Xp[Y.np] = Y[0];
      Spline.tp[Y.np] = tp[0];
      Spline.bc[Y.np] = bc[0];
      F_Type[Y.np] = ft[0];
      F[Y.np] = f[0];
    }
    Spline.pathlength();

    // Set the interface of union type.
    Type = INTERFACE_UNION;
    // Set the interface of union boundary condition type.
    if (I1.BC_Type == I2.BC_Type) BC_Type = I1.BC_Type;
    else BC_Type = INTERFACE_BC_MIXED;
    // Set the interface of union motion type.
    if (I1.Motion == I2.Motion) Motion = I1.Motion;
    else Motion = INTERFACE_MOTION_MIXED;
    // Determine the reference point (possibly the centroid) of the
    // interface of union.
    if (I1.Type == INTERFACE_USER_SPECIFIED ||
	I2.Type == INTERFACE_USER_SPECIFIED) {
      if (I1.Type != INTERFACE_USER_SPECIFIED) Xref = I1.Xref;
      else Xref = I2.Xref;
    } else {
      Centroid();
    }
    // Determine the bounding box of the interface of union.
    BoundingBox();

  }

  // Deallocate memory for I1.
  I1.deallocate();

  // Interface union successful.
  return 0;

}

/**********************************************************************
 * Interface2D::Interface_Intersection -- Return the interface of     *
 *                                        intersection between the    *
 *                                        current and given           *
 *                                        interfaces.  The scheme     *
 *                                        used here is based on the   *
 *                                        Weiler-Atherton algorithm   *
 *                                        popular in computer         *
 *                                        graphics.                   *
 **********************************************************************/
void Interface2D::Interface_Intersection(Interface2D &I2) {

  Polygon P, P1, P2;
  int n1, n2, type;

  // Convert the interface splines into polygons.
  //P1.convert(Spline);
  //P2.convert(I2.Spline);

  // Determine the polygon of intersection.
  Polygon_Intersection(P1,P2,P);

  // Convert the polygon into a spline and determine the appropriate
  // boundary condition and type of each point in the interface and
  // calculate the path lengths.  Also determine the velocity function
  // at each spline point.
  type = min(Spline.type,I2.Spline.type);
  Spline.allocate(P.np+1);
  Spline.settype(type);
  for (int n = 0; n < P.np; n++) {
    Spline.Xp[n] = P.X[n];
    Spline.tp[n] = SPLINE2D_POINT_NORMAL;
    Spline.bc[n] = BC_REFLECTION;
  }
  Spline.Xp[Spline.np-1] = Spline.Xp[0];
  Spline.tp[Spline.np-1] = Spline.tp[0];
  Spline.bc[Spline.np-1] = Spline.bc[0];
  Spline.pathlength();
  // Determine the bounding box.
  BoundingBox();
  // Determine the reference point.
  Centroid();
  // Initialize the interface velocity function.
  Initialize_Velocity_Function();

}

/**********************************************************************
 * Interface2D::Check_Interface_Intersection -- Determine if the      *
 *                                              current interface     *
 *                                              intersects with the   *
                                                given interface.      *
 **********************************************************************/
int Interface2D::Check_Interface_Intersection(Interface2D &I2) {
  if (Bounding_Box_Intersection(I2.Xmin,I2.Xmax)) return 1;
  int inside_flag = OFF;
  for (int np = 0; np < I2.Spline.np; np++)
    if (Point_In_Interface(I2.Spline.Xp[np])) { inside_flag = ON; break; }
  if (inside_flag) return 1;
  for (int np = 0; np < Spline.np; np++)
    if (I2.Point_In_Interface(Spline.Xp[np])) { inside_flag = ON; break; }
  return inside_flag;
}

/**********************************************************************
 * Interface2D_List::Define member functions.                         *
 **********************************************************************/

/**********************************************************************
 * Interface2D_List::Construct_Union_List -- Construct an interface   *
 *                                           union list from a list   *
 *                                           of interface components. *
 **********************************************************************/
int Interface2D_List::Construct_Union_List(Interface2D_List &Component_List) {

  // Exit immediately if no interface components have been created.
  if (!Component_List.Ni) return 0;

  // Deallocate the relations for the components.
  for (int ni = 0; ni < Component_List.Ni; ni++) {
    Component_List.relations[ni].deallocate();
  }

  // Create the interface union list and exit immediately if only one
  // interface components have been created.
  if (Component_List.Ni == 1) {
    allocate(Component_List.Ni);
    Component_List.relations[0].add(1);
    Interface[0] = Component_List[1];
    relations[0].add(1);
    // Union list constructed successfully.
    return 0;
  }

  int *used_flag, num, error_flag;
  Interface2D *Union_List, Union_Interface;

  // Initialize interface used flag array.
  used_flag = new int[Component_List.Ni];
  for (int ni = 0; ni < Component_List.Ni; ni++) used_flag[ni] = 0;
  Union_List = new Interface2D[Component_List.Ni];
  num = 0;

  // Construct the interface union list from the interface components.
  for (int current = 0; current < Component_List.Ni; current++) {
    if (!used_flag[current]) {
      // Add the interface component to the union list.
      Union_List[num] = Component_List.Interface[current];
      num++;
      Component_List.relations[current].add(num);
      used_flag[current] = 1;
      // Determine interface unions with the other interface components.
      for (int next = current+1; next < Component_List.Ni; next++) {
    	if (!used_flag[next] &&
	    Component_List.Interface[current].Check_Interface_Intersection(Component_List.Interface[next])) {
   	  // Current and next components intersect.  Determine union.
	  Union_Interface.Copy(Union_List[num-1]);
	  error_flag = Union_Interface.Interface_Union(Component_List.Interface[next]);
	  //Union_List[num-1].Interface_Union(Component_List.Interface[next]);
    	  if (!error_flag) {
	    Union_List[num-1].Copy(Union_Interface);
	    Component_List.relations[next].add(num);
	    used_flag[next] = 1;
	  }
 	}
      }
    }
  }

  // Create union interface list and relations list.
  allocate(num);
  for (int ni = 0; ni < Ni; ni++) {
    Interface[ni] = Union_List[ni];
    for (int nc = 0; nc < Component_List.Ni; nc++) {
      for (int nr = 0; nr < Component_List.relations[nc].np; nr++) {
	if (Component_List.relations[nc][nr] == ni) {
  	  relations[ni].add(nc+1);
	}
      }
    }
  }

  // Deallocate memory.
  delete []Union_List; Union_List = NULL;
  delete []used_flag; used_flag = NULL;

  // Union list constructed successfully.
  return 0;

}

/**********************************************************************
 * Interface2D -- External subroutines.                               *
 **********************************************************************/

/**********************************************************************
 * Routine: Broadcast_Interface                                       *
 *                                                                    *
 * Broadcasts an interface to all processors involved in the          *
 * calculation from the primary processor using the MPI broadcast     *
 * routine.                                                           *
 *                                                                    *
 **********************************************************************/
void Broadcast_Interface(Interface2D &I) {

#ifdef _MPI_VERSION
  int buffer_size, i_buffer_size;
  int *i_buffer;
  double *buffer;

  // Broadcast the interface type, bc_type, and motion type..
  i_buffer_size = 3;
  i_buffer = new int[i_buffer_size];
  if (CFFC_Primary_MPI_Processor()) {
    i_buffer[0] = I.Type;
    i_buffer[1] = I.BC_Type;
    i_buffer[2] = I.Motion;
  }
  MPI::COMM_WORLD.Bcast(i_buffer,i_buffer_size,MPI::INT,0);
  if (!CFFC_Primary_MPI_Processor()) {
    I.Type = i_buffer[0];
    I.BC_Type = i_buffer[1];
    I.Motion = i_buffer[2];
  }
  delete []i_buffer; i_buffer = NULL;

  // Broadcast the interface characteristic length, characteristic
  // velocity, and reference point.
  buffer_size = 6;
  buffer = new double[buffer_size];
  if (CFFC_Primary_MPI_Processor()) {
    buffer[0] = I.Length1;
    buffer[1] = I.Length2;
    buffer[2] = I.Speed.x;
    buffer[3] = I.Speed.y;
    buffer[4] = I.Xref.x;
    buffer[5] = I.Xref.y;
  }
  MPI::COMM_WORLD.Bcast(buffer,buffer_size,MPI::DOUBLE,0);
  if (!CFFC_Primary_MPI_Processor()) {
    I.Length1 = buffer[0];
    I.Length2 = buffer[1];
    I.Speed = Vector2D(buffer[2],buffer[3]);
    I.Xref = Vector2D(buffer[4],buffer[5]);
  }
  delete []buffer; buffer = NULL;

  // Broadcast the spline.
  Broadcast_Spline(I.Spline);

  // On non-primary processors initialize the velocity function type 
  // and the velocity function.
  if (!CFFC_Primary_MPI_Processor()) I.Initialize_Velocity_Function(I.Motion);

  // Broadcast the velocity function type and the velocity function.
  i_buffer_size = I.Spline.np;
  i_buffer = new int[i_buffer_size];
  buffer_size = 2*I.Spline.np;
  buffer = new double[buffer_size];

  if (CFFC_Primary_MPI_Processor()) {
    buffer_size = 0;
    for (int n = 0; n < I.Spline.np; n++) {
      i_buffer[n] = I.F_Type[n];
      buffer[buffer_size  ] = I.F[n].x;
      buffer[buffer_size+1] = I.F[n].y;
      buffer_size += 2;
    }
  }

  buffer_size = 2*I.Spline.np;
  MPI::COMM_WORLD.Bcast(i_buffer,i_buffer_size,MPI::INT,0);
  MPI::COMM_WORLD.Bcast(buffer,buffer_size,MPI::DOUBLE,0);

  if (!CFFC_Primary_MPI_Processor()) {
    buffer_size = 0;
    for (int n = 0; n < I.Spline.np; n++) {
      I.F_Type[n] = i_buffer[n];
      I.F[n].x = buffer[buffer_size  ];
      I.F[n].y = buffer[buffer_size+1];
      buffer_size += 2;
    }
  }

  delete []i_buffer; i_buffer = NULL;
  delete []buffer; buffer = NULL;

  // Compute the bounding box on non-primary processors.
  //if (!CFFC_Primary_MPI_Processor()) I.BoundingBox();
  I.BoundingBox();

#endif

}

#ifdef _MPI_VERSION
/************************************************************************
 * Routine: Broadcast_Interface                                         *
 *                                                                      *
 * Broadcasts an interface to all processors associated with the        *
 * specified communicator from the specified processor using the MPI    *
 * broadcast routine.                                                   *
 *                                                                      *
 ************************************************************************/
void Broadcast_Interface(Interface2D &I,
			 MPI::Intracomm &Communicator, 
			 const int Source_CPU) {

  int Source_Rank = 0;
  int buffer_size, i_buffer_size;
  int *i_buffer;
  double *buffer;
 
  // Broadcast the interface type, bc_type, and motion type..
  i_buffer_size = 3;
  i_buffer = new int[i_buffer_size];
  if (CFFC_MPI::This_Processor_Number == Source_CPU) {
    i_buffer[0] = I.Type;
    i_buffer[1] = I.BC_Type;
    i_buffer[2] = I.Motion;
  }
  Communicator.Bcast(i_buffer,i_buffer_size,MPI::INT,Source_Rank);
  if (CFFC_MPI::This_Processor_Number != Source_CPU) {
    I.Type = i_buffer[0];
    I.BC_Type = i_buffer[1];
    I.Motion = i_buffer[2];
  }
  delete []i_buffer; i_buffer = NULL;

  // Broadcast the interface characteristic length, characteristic
  // velocity, and reference point.
  buffer_size = 6;
  buffer = new double[buffer_size];
  if (CFFC_MPI::This_Processor_Number == Source_CPU) {
    buffer[0] = I.Length1;
    buffer[1] = I.Length2;
    buffer[2] = I.Speed.x;
    buffer[3] = I.Speed.y;
    buffer[4] = I.Xref.x;
    buffer[5] = I.Xref.y;
  }
  Communicator.Bcast(buffer,buffer_size,MPI::DOUBLE,Source_Rank);
  if (CFFC_MPI::This_Processor_Number != Source_CPU) {
    I.Length1 = buffer[0];
    I.Length2 = buffer[1];
    I.Speed = Vector2D(buffer[2],buffer[3]);
    I.Xref = Vector2D(buffer[4],buffer[5]);
  }
  delete []buffer; buffer = NULL;

  // Broadcast the spline.
  Broadcast_Spline(I.Spline,Communicator,Source_CPU);

  // On non-primary processors initialize the velocity function type 
  // and the velocity function.
  if (CFFC_MPI::This_Processor_Number != Source_CPU)
    I.Initialize_Velocity_Function();

  // Broadcast the velocity function type and the velocity function.
  i_buffer_size = I.Spline.np;
  i_buffer = new int[i_buffer_size];
  buffer_size = 2*I.Spline.np;
  buffer = new double[buffer_size];

  if (CFFC_MPI::This_Processor_Number == Source_CPU) {
    buffer_size = 0;
    for (int n = 0; n < I.Spline.np; n++) {
      i_buffer[n] = I.F_Type[n];
      buffer[buffer_size  ] = I.F[n].x;
      buffer[buffer_size+1] = I.F[n].y;
      buffer_size += 2;
    }
  }

  buffer_size = 2*I.Spline.np;
  Communicator.Bcast(i_buffer,i_buffer_size,MPI::INT,Source_Rank);
  Communicator.Bcast(buffer,buffer_size,MPI::DOUBLE,Source_Rank);

  if (CFFC_MPI::This_Processor_Number != Source_CPU) {
    buffer_size = 0;
    for (int n = 0; n < I.Spline.np; n++) {
      I.F_Type[n] = i_buffer[n];
      I.F[n].x = buffer[buffer_size  ];
      I.F[n].y = buffer[buffer_size+1];
      buffer_size += 2;
    }
  }

  delete []i_buffer; i_buffer = NULL;
  delete []buffer; buffer = NULL;

  // Compute the bounding box on non-primary processors.
  //if (CFFC_MPI::This_Processor_Number != Source_CPU) I.BoundingBox();
  I.BoundingBox();

}
#endif

/**********************************************************************
 * Routine: Intersection_Area                                         *
 *                                                                    *
 * This routine returns the area of intersection between two          *
 * interfaces assuming linear segments.                               *
 *                                                                    *
 **********************************************************************/
double Intersection_Area(Interface2D &I1, Interface2D &I2) {

  cout << endl << " ERROR -- Intersection area calculation.";
  return 0.0; // for compiler

//   Polygon P1, P2;

//   // Convert the interface splines into polygons.
//   P1.convert(I1.Spline);
//   P2.convert(I2.Spline);

//   return Polygon_Intersection_Area(P1,P2);

}

/**********************************************************************
 * Interface2D_List -- External subroutines.                          *
 **********************************************************************/

/**********************************************************************
 * Routine: Broadcast_Interface_List                                  *
 *                                                                    *
 * Broadcasts an interface list to all processors involved in the     *
 * calculation from the primary processor using the MPI broadcast     *
 * routine.                                                           *
 *                                                                    *
 **********************************************************************/
void Broadcast_Interface_List(Interface2D_List &List) {

#ifdef _MPI_VERSION
  int Ni, buffer_size, *buffer;

  // Broadcast the number of interfaces in the list.
  if (CFFC_Primary_MPI_Processor()) Ni = List.Ni;
  else Ni = 0;

  MPI::COMM_WORLD.Bcast(&Ni,1,MPI::INT,0);

  // Allocate memory if not the primary CPU.
  if (!CFFC_Primary_MPI_Processor()) List.allocate(Ni);

  // Broadcast all of the interfaces in the list.
  for (int n = 0; n < Ni; n++) Broadcast_Interface(List.Interface[n]);

  // Broadcast the relations linked-lists.
  for (int n = 0; n < Ni; n++) {
    // Broadcast the number of points in the linked-list.
    if (CFFC_Primary_MPI_Processor()) buffer_size = List.relations[n].np;
    else buffer_size = 0;
    MPI::COMM_WORLD.Bcast(&buffer_size,1,MPI::INT,0);
    // Allocate message buffer.
    buffer = new int[buffer_size];
    // Enter data into the message buffer on the primary processor.
    if (CFFC_Primary_MPI_Processor()) {
      for (int i = 0; i < buffer_size; i++) buffer[i] = List.relations[n][i];
    } else {
      for (int i = 0; i < buffer_size; i++) buffer[i] = 0;
    }
    // Broadcast message buffer.
    MPI::COMM_WORLD.Bcast(buffer,buffer_size,MPI::INT,0);
    // Copy data from the message buffer on the non-primary processors.
    if (!CFFC_Primary_MPI_Processor()) {
      for (int i = 0; i < buffer_size; i++) {
	List.relations[n].add(buffer[i]);
      }
    }
    delete []buffer; buffer = NULL;
  }
#endif

}

#ifdef _MPI_VERSION
/************************************************************************
 * Routine: Broadcast_Interface_List                                    *
 *                                                                      *
 * Broadcasts an interface list to all processors associated with the   *
 * specified communicator from the specified processor using the MPI    *
 * broadcast routine.                                                   *
 *                                                                      *
 ************************************************************************/
void Broadcast_Interface_List(Interface2D_List &List,
			      MPI::Intracomm &Communicator, 
			      const int Source_CPU) {

  int Ni, Source_Rank = 0;
  int *buffer, buffer_size;

  // Broadcast the number of interfaces in the list.
  if (CFFC_MPI::This_Processor_Number == Source_CPU) Ni = List.Ni;
  else Ni = 0;

  Communicator.Bcast(&Ni,1,MPI::INT,Source_Rank);

  // Allocate memory if not the source CPU.
  if (CFFC_MPI::This_Processor_Number != Source_CPU) List.allocate(Ni);

  // Broadcast all of the interfaces in the list.
  for (int n = 0; n < Ni; n++) Broadcast_Interface(List.Interface[n],
						   Communicator,
						   Source_CPU);

  // Broadcast the relations linked-lists.
  for (int n = 0; n < Ni; n++) {
    // Broadcast the number of points in the linked-list.
    if (CFFC_MPI::This_Processor_Number == Source_CPU) buffer_size = List.relations[n].np;
    else buffer_size = 0;
    Communicator.Bcast(&buffer_size,1,MPI::INT,Source_Rank);
    // Allocate message buffer.
    buffer = new int[buffer_size];
    // Enter data into the message buffer on the source processor.
    if (CFFC_MPI::This_Processor_Number == Source_CPU) {
      for (int i = 0; i < buffer_size; i++) buffer[i] = List.relations[n][i];
    } else {
      for (int i = 0; i < buffer_size; i++) buffer[i] = 0;
    }
    // Broadcast message buffer.
    Communicator.Bcast(buffer,buffer_size,MPI::INT,Source_Rank);
    // Copy data from the message buffer on the non-source processors.
    if (CFFC_MPI::This_Processor_Number != Source_CPU) {
      for (int i = 0; i < buffer_size; i++) {
	List.relations[n].add(buffer[i]);
      }
    }
    delete []buffer; buffer = NULL;
  }

}
#endif

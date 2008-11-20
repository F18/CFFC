/*!\file HO_Spline2D.cc
   \brief Source file initializing/implementing member variables/functions that belong to classes defined in HO_Spline2D.h */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "HO_Spline2D.h"	// Include 2D high-order spline header file


// ===== Member variables =====
int Spline2D_HO::CounterSolidBodies = 0; // initialize counter

// ===== Member functions =====

/*!
 * Concatenation operator.
 *
 * \note The properties of the new spline are inherited from the spline to the left of the '+' sign. (i.e. first one)
 */
const Spline2D_HO Spline2D_HO::operator + (const Spline2D_HO &S) const{
  int i, npts; Spline2D_HO Sc;
  npts = np + S.np - 1;
  Sc.allocate(npts); Sc.settype(type);
  Sc.setFluxCalcMethod(getFluxCalcMethod());
  Sc.bodyOwner = bodyOwner;
  for ( i = 0; i <= np-1; ++i ) {
    Sc.Xp[i] = Xp[i];
    Sc.tp[i] = tp[i];
    Sc.bc[i] = bc[i];
  } /* endfor */
  Sc.tp[np-1] = SPLINE2D_POINT_SHARP_CORNER;
  for ( i = 1; i <= S.np-1; ++i ) {
    Sc.Xp[i+np-1] = S.Xp[i] + (Xp[np-1]-S.Xp[0]);
    Sc.tp[i+np-1] = S.tp[i];
    Sc.bc[i+np-1] = S.bc[i];
  } /* endfor */
  Sc.pathlength();
  return (Sc);
}

/*!
 * Translates or shifts the positions of all of the 
 * points defining the spline.                      
 *                                                  
 */
const Spline2D_HO& Spline2D_HO::Translate_Spline(const Vector2D &V){
  int i;
  for ( i = 0; i <= np-1; ++i ) {
    Xp[i] += V;
  } /* endfor */
  return (*this);
}

/*!
 * Apply a solid body rotation about the origin and 
 * recompute the positions of all of the points in  
 * the spline accordingly.                          
 */
const Spline2D_HO& Spline2D_HO::Rotate_Spline(const double &a){
  int i;
  double cos_angle, sin_angle;
  Vector2D X;

  cos_angle = cos(-a); sin_angle = sin(-a);
  for ( i = 0; i <= np-1; ++i ) {
    X.x = ( Xp[i].x*cos_angle +
	    Xp[i].y*sin_angle );
    X.y = ( - Xp[i].x*sin_angle +
	    Xp[i].y*cos_angle );
    Xp[i] = X;
  } /* endfor */
  return (*this);
}
 
/*!
 * Scale the positions of all of 
 * the points defining the spline.
 */
const Spline2D_HO& Spline2D_HO::Scale_Spline(const double &a){
  int i;
  for ( i = 0; i <= np-1; ++i ) {
    Xp[i] = Xp[i]*a;
    sp[i] *= a;
  } /* endfor */
  return (*this);
}

/*!
 * Performa a piecewise linear, quadratic, cubic
 * or quintic blended spline interpolation to 
 * return the 2D position vector on  
 * the splined surface (body or boundary) of interest 
 * given the path length s along the spline and 
 * the set of np discrete points (Xp.x,Xp.y) 
 * that define the spline geometry.  
 */
const Vector2D Spline2D_HO::Spline(const double &s) const {

  int il, ir, npts_used;
  double sp_used[4];
  double ds0, ds1, ds2, ds3,
         ds01, ds02, ds12, ds13, ds23,
         a0, a1, a2, a3, a_blend;
  Vector2D Xp_used[4], Xl, Xr; 
  
  // Check for domain of existance
  if (s <= sp[0]) return(Xp[0]);
  if (s >= sp[np-1]) return(Xp[np-1]);      
  
  /* Determine the subinterval of the spline
     containing the point (position) of interest. */
  find_subinterval(s,il,ir);

  /* Determine number and which points to use in the 
     blended spline interpolation. */

  if (type == SPLINE2D_CONSTANT) return(Xp[il]);
  if (type == SPLINE2D_LINEAR) {
    npts_used = 2;
  } else if (tp[il] == SPLINE2D_POINT_SHARP_CORNER &&
	     tp[ir] == SPLINE2D_POINT_SHARP_CORNER) {
    npts_used = 2;
  } else if (tp[il] == SPLINE2D_POINT_NORMAL &&
	     tp[ir] == SPLINE2D_POINT_SHARP_CORNER) {
    npts_used = 3;
    if (il == 0) {
      Xp_used[0]=Xp[np-2];
      sp_used[0]=sp[0]-(sp[np-1]-sp[np-2]);
      Xp_used[1]=Xp[il];
      sp_used[1]=sp[il];
      Xp_used[2]=Xp[ir];
      sp_used[2]=sp[ir];
    } else {
      Xp_used[0]=Xp[il-1];
      sp_used[0]=sp[il-1];
      Xp_used[1]=Xp[il];
      sp_used[1]=sp[il];
      Xp_used[2]=Xp[ir];
      sp_used[2]=sp[ir];
    } /* endif */
  } else if (tp[il] == SPLINE2D_POINT_SHARP_CORNER &&
	     tp[ir] == SPLINE2D_POINT_NORMAL) {
    npts_used = 3;
    if (il == np-2) {
      Xp_used[0]=Xp[il];
      sp_used[0]=sp[il];
      Xp_used[1]=Xp[ir];
      sp_used[1]=sp[ir];
      Xp_used[2]=Xp[1];
      sp_used[2]=sp[np-1]+(sp[1]-sp[0]);
    } else {
      Xp_used[0]=Xp[il];
      sp_used[0]=sp[il];
      Xp_used[1]=Xp[ir];
      sp_used[1]=sp[ir];
      Xp_used[2]=Xp[ir+1];
      sp_used[2]=sp[ir+1];
    } /* endif */
  } else if (tp[il] == SPLINE2D_POINT_NORMAL &&
	     tp[ir] == SPLINE2D_POINT_NORMAL) {
    npts_used = 4;
    if (il == 0) {
      Xp_used[0]=Xp[np-2];
      sp_used[0]=sp[0]-(sp[np-1]-sp[np-2]);
      Xp_used[1]=Xp[il];
      sp_used[1]=sp[il];
      Xp_used[2]=Xp[ir];
      sp_used[2]=sp[ir];
      Xp_used[3]=Xp[ir+1];
      sp_used[3]=sp[ir+1];
    } else if (il == np-2) {
      Xp_used[0]=Xp[il-1];
      sp_used[0]=sp[il-1];
      Xp_used[1]=Xp[il];
      sp_used[1]=sp[il];
      Xp_used[2]=Xp[ir];
      sp_used[2]=sp[ir];
      Xp_used[3]=Xp[1];
      sp_used[3]=sp[np-1]+(sp[1]-sp[0]);
    } else {
      Xp_used[0]=Xp[il-1];
      sp_used[0]=sp[il-1];
      Xp_used[1]=Xp[il];
      sp_used[1]=sp[il];
      Xp_used[2]=Xp[ir];
      sp_used[2]=sp[ir];
      Xp_used[3]=Xp[ir+1];
      sp_used[3]=sp[ir+1];
    } /* endif */
  } /* endif */

    /* Perform the interpolation. */

  switch(npts_used){

  case 4:		// four points
    ds0 = s - sp_used[0];
    ds1 = s - sp_used[1];
    ds2 = s - sp_used[2];
    ds3 = s - sp_used[3];
       
    ds01 = sp_used[0] - sp_used[1];
    ds02 = sp_used[0] - sp_used[2];
    ds12 = sp_used[1] - sp_used[2];
    ds13 = sp_used[1] - sp_used[3];
    ds23 = sp_used[2] - sp_used[3];
       
    a0 =  ds1*ds2/(ds01*ds02);
    a1 = -ds0*ds2/(ds01*ds12);
    a2 =  ds0*ds1/(ds02*ds12);

    Xl = a0*Xp_used[0] + a1*Xp_used[1] +
      a2*Xp_used[2];

    a1 =  ds2*ds3/(ds12*ds13);
    a2 = -ds1*ds3/(ds12*ds23);
    a3 =  ds1*ds2/(ds13*ds23);

    Xr = a1*Xp_used[1] + a2*Xp_used[2] +
      a3*Xp_used[3];

    switch(type) {
    case SPLINE2D_QUADRATIC :
      Xl = HALF*(Xl+Xr);
      break; 
    case SPLINE2D_CUBIC :
      a_blend = ds1/ds12;
      Xl = (ONE+a_blend)*Xl - a_blend*Xr;
      break;
    case SPLINE2D_QUINTIC :
      a_blend = ds1/ds12;
      a_blend =-a_blend*a_blend*(THREE+TWO*a_blend);
      Xl = (ONE+a_blend)*Xl - a_blend*Xr;
      break;
    default:
      Xl = HALF*(Xl+Xr);
      break;
    } /* endswitch */

    return(Xl);

  case 3:		// three points
    ds0 = s - sp_used[0];
    ds1 = s - sp_used[1];
    ds2 = s - sp_used[2];
       
    ds01 = sp_used[0] - sp_used[1];
    ds02 = sp_used[0] - sp_used[2];
    ds12 = sp_used[1] - sp_used[2];
       
    a0 =  ds1*ds2/(ds01*ds02);
    a1 = -ds0*ds2/(ds01*ds12);
    a2 =  ds0*ds1/(ds02*ds12);
         
    Xl = a0*Xp_used[0] + a1*Xp_used[1] +
      a2*Xp_used[2];
    return(Xl);

  case 2:		// two points
    Xl = Xp[il]+(s-sp[il])*(Xp[ir]-Xp[il])/
      (sp[ir]-sp[il]);
    return(Xl);

  }/* endswitch */

  return Vector2D(0);
}


/*!
 * Get the unit tangential vector at the point of interest,
 * based on the point position vector.
 *
 * \param [in] Point the position vector
 * \param [in] Order the order of the spline interpolant
 * \param [out] dxds the derivative of x-coordinate with respect to the spline path length
 * \param [out] dyds the derivative of y-coordinate with respect to the spline path length
 */
const Vector2D Spline2D_HO::tSpline(const Vector2D &Point, const PolynomOrder Order,
				    double & dxds, double & dyds) const {

  int NumOfPoints(Order+1);	// Set the number of point used to determine the polynom based on Order
  int il, ir;
  Vector2D *Xp_used = new Vector2D [NumOfPoints];
  double *sp_used = new double [NumOfPoints];
  Vector2D Tangent;
  double s;

  /* Step1. Find the path length for Point */
  s = getS(Point);

  /* Step2. Analyse the pathlength */
  if ( (s>sp[np-1]) && ((s - sp[np-1])/(1.0 + sp[np-1]) < EpsilonTol::epsilon_relative) ) {
    s = sp[np-1];
  }

  if ( (s < sp[0]) || ( s > sp[np-1]) ){
    return Vector2D(0.0,0.0);	// if the point is not on the Spline the normal vector is Zero
  }

  /* Step3. Determine the subinterval of the spline
     which contains the point (position) of interest. */
  find_subinterval(s,il,ir);

  // Step4. Set the vector of points used to compute the tangent
  switch(NumOfPoints){
  case 4:
    // use a 3rd-order polynom
    sp_used[0] = sp[il];                             Xp_used[0] = Xp[il];
    sp_used[1] = sp[il] + 0.3*(sp[ir] - sp[il]);     Xp_used[1] = Spline(sp_used[1]);
    sp_used[2] = sp[il] + 0.7*(sp[ir] - sp[il]);     Xp_used[2] = Spline(sp_used[2]);
    sp_used[3] = sp[ir];                             Xp_used[3] = Xp[ir];

    break;
  case 3:
    // use a 2nd-order polynom
    sp_used[0] = sp[il];                             Xp_used[0] = Xp[il];
    sp_used[1] = sp[il] + 0.5*(sp[ir] - sp[il]);     Xp_used[1] = Spline(sp_used[1]);
    sp_used[2] = sp[ir];                             Xp_used[2] = Xp[ir];

    break;

  default:
    // use a 1st-order polynom
    sp_used[0] = sp[il];                    Xp_used[0] = Xp[il];
    sp_used[1] = sp[ir];                    Xp_used[1] = Xp[ir];
  }

  // Step5. Compute the tangent
  Tangent = tSpline(s,Xp_used,sp_used,NumOfPoints,dxds,dyds);

  // Deallocate the memory
  delete [] Xp_used; Xp_used = NULL;
  delete [] sp_used; sp_used = NULL;

  // Return
  return Tangent;
}

/*!
 * Get the unit tangential vector at the point of interest
 * based on the point path coordinate.
 *
 * \param [in] s the point path coordinate
 * \param [in] Xp_used the vector of used spline nodes
 * \param [in] s_used
 * \param [in] NumOfControlPoints the dimension of Xp_used
 * \param [out] dxds the derivative of x-coordinate with respect to the spline path length
 * \param [out] dyds the derivative of y-coordinate with respect to the spline path length
 */
const Vector2D Spline2D_HO::tSpline(const double & s, const Vector2D * Xp_used,
				    const double * s_used, const int & NumOfControlPoints,
				    double & dxds, double & dyds) const {

  double ds21, ds32, ds31, ds41, ds42, ds43;	        // coordinate differences in the path length space
  Vector2D DD12, DD23, DD34, DD123, DD234, DD1234;	// divided differences

  Vector2D Tangent;

  switch(NumOfControlPoints){

  case 4:
    ds21 = s_used[1] - s_used[0];
    ds32 = s_used[2] - s_used[1];
    ds31 = s_used[2] - s_used[0];
    ds41 = s_used[3] - s_used[0];
    ds42 = s_used[3] - s_used[1];
    ds43 = s_used[3] - s_used[2];

    // First divided differences
    DD12  = (Xp_used[1] - Xp_used[0])/ds21;
    DD23  = (Xp_used[2] - Xp_used[1])/ds32;
    DD34  = (Xp_used[3] - Xp_used[2])/ds43;

    // Second divided differences
    DD123 = (DD23 - DD12)/ds31;
    DD234 = (DD34 - DD23)/ds42;

    // Third divided differences
    DD1234 = (DD234 - DD123)/ds41;

    Tangent = (3.0*DD1234*s*s + 2.0*(DD123 - (s_used[0]+s_used[1]+s_used[2])*DD1234)*s +
	       (DD12 - (s_used[0]+s_used[1])*DD123 + (s_used[0]*s_used[1] + 
						      s_used[0]*s_used[2] +
						      s_used[1]*s_used[2])*DD1234) );

    // Set the derivative of x and y with respect to the pathlength
    dxds = Tangent.x;
    dyds = Tangent.y;

    // Normalize the tangent vector
    return Tangent/abs(Tangent);

  case 3:			// quadratic approximation

    ds21 = s_used[1] - s_used[0];
    ds32 = s_used[2] - s_used[1];
    ds31 = s_used[2] - s_used[0];
    
    // First divided differences
    DD12  = (Xp_used[1] - Xp_used[0])/ds21;
    DD23  = (Xp_used[2] - Xp_used[1])/ds32;

    // Second divided differences
    DD123 = (DD23 - DD12)/ds31;

    Tangent = (2.0*s - s_used[0] - s_used[1])*DD123 + DD12;

    // Set the derivative of x and y with respect to the pathlength
    dxds = Tangent.x;
    dyds = Tangent.y;

    // Normalize the tangent vector
    return Tangent/abs(Tangent);

  case 2:			// linear approximation

    Tangent = (Xp_used[1] - Xp_used[0])/(s_used[1] - s_used[0]);

    // Set the derivative of x and y with respect to the pathlength
    dxds = Tangent.x;
    dyds = Tangent.y;

    // Normalize the tangent vector
    return Tangent/abs(Tangent);

  default:
    // return a ZERO tangent
    return Vector2D(0.0,0.0);
  }
}

/*!
 * Get the normal unit vector at the point of interest
 * by approximating the spline in the vicinity of Point 
 * with a straight segment. 
 *
 * \param [in] Point the position vector
 */
const Vector2D Spline2D_HO::nSpline_SegmentBased(const Vector2D &Point) const {

  double s;
  int il, ir;
  Vector2D StartPoint, EndPoint;
  double VDist;
  static const double Multiplier(1.0e-5);

  /*Step1. Find the path length for Point */
  s = getS(Point);

  /*Step2. Analyse the pathlength */
  if ( (s>sp[np-1]) && ((s - sp[np-1])/(1.0 + sp[np-1]) < EpsilonTol::epsilon_relative) ) {
    s = sp[np-1];
  }

  if ( (s < sp[0]) || ( s > sp[np-1]) ){
    return Vector2D(0.0,0.0);	// if the point is not on the Spline the normal vector is Zero
  }

  /*Step3. Determine the subinterval of the spline
    which contains the point (position) of interest. */
  find_subinterval(s,il,ir);

  /*Step4. Determine the vicinity distance as the minimum distance among the neighbours. */
  if ( il==0 ){
    if (np > 2){
      VDist = Multiplier*min( (sp[ir] - sp[il]), (sp[ir+1]-sp[ir]));
    } else {
      VDist = Multiplier*(sp[ir] - sp[il]);
    }
  } else if ( ir==(np-1) ){
    VDist = Multiplier*min( (sp[il] - sp[il-1]), (sp[ir]-sp[il]));
  } else {
    VDist = Multiplier*min( (sp[il] - sp[il-1]), min((sp[ir] - sp[il]),(sp[ir+1]-sp[ir]) ));
  }

  /*Step5. Pick the two points used for determining the normal */
  if ( (s>sp[il]) && (s<sp[ir]) ){ // interior point

    // StartPoint
    if ( (s-VDist <= sp[il]) && (tp[il] == SPLINE2D_POINT_SHARP_CORNER) ){
      StartPoint = Xp[il];
    } else {
      StartPoint = Spline(s-VDist);
    }

    // EndPoint
    if ( (s+VDist >= sp[ir]) && (tp[il] == SPLINE2D_POINT_SHARP_CORNER) ){
      EndPoint = Xp[ir];
    } else {
      EndPoint = Spline(s+VDist);
    }

  } else if (s == sp[il]){	   // left end

    if (tp[il] == SPLINE2D_POINT_SHARP_CORNER ){
      /* sharp corner point -> use the one sided normal based on the segment (il-ir)
	 Obs. This approach might not give the best result for all situations
       */
      // StartPoint
      StartPoint = Xp[il];
      // EndPoint
      EndPoint = Spline(s+VDist);


    } else {			
      /* normal point */

      // StartPoint
      StartPoint = Spline(s-VDist);
      // EndPoint
      EndPoint = Spline(s+VDist);

    }

  } else {	   // right end -> (s == sp[ir])

    if (tp[ir] == SPLINE2D_POINT_SHARP_CORNER ){
      /* sharp corner point -> use the one sided normal based on the segment (il-ir)
	 Obs. This approach might not give the best result for all situations
       */

      // StartPoint
      StartPoint = Spline(s-VDist);
      // EndPoint
      EndPoint = Xp[ir];

    } else {
      /* normal point */

      // StartPoint
      StartPoint = Spline(s-VDist);
      // EndPoint
      EndPoint = Spline(s+VDist);

    }
  }

  /*Step6. Analyze the chosen points 
    Obs. Because they might be undistingushed numerically, this step is a must! */
  if ( (fabs(EndPoint.x - StartPoint.x) < 1.0e-14 ) && (fabs(EndPoint.y - StartPoint.y) < 1.0e-14 ) ){
    StartPoint = Xp[il];
    EndPoint = Xp[ir];
  }

  /*Step7. Compute the normal vector */
  return NormalVector(StartPoint,EndPoint);
}

/*!
 * Output the spline to a stream in a format 
 * suitable for plotting with Tecplot.
 *
 * \param output_file the output stream
 * \param NumberOfPoints the number of point used for plotting the spline.
 */
void Spline2D_HO::OutputTecplot(std::ostream &output_file, const int NumberOfPoints){

  int SplinePoints, NumberOfInteriorPoints;

  output_file << setprecision(14);
  output_file << "TITLE = Spline2D Representation\n ";
  output_file << "VARIABLES = \"x\" \\ \n"
	      << "\"y\" \n"
	      << "ZONE \n";

  if (NumberOfPoints <= np){
    SplinePoints = np; 		// The number of points used to represent the spline
                                // cannot be smaller than the number of control points
  } else {
    SplinePoints = NumberOfPoints;
  }

  // Generate the curvilinear coordinates of the points used to define the spline
  NumberOfInteriorPoints = SplinePoints - np;

  double DeltaS;
  if(NumberOfInteriorPoints != 0){
    DeltaS = sp[np-1]/(NumberOfInteriorPoints + 2);
  } else {
    // Write the coordinates of the control points and Return
    for(int i = 0; i<=np-1; ++i){
      output_file << Xp[i].x << "\t" << Xp[i].y << "\n";
    }
    return;
  }

  double CurrentPath = 0.0;
  Vector2D PointOnCurve;

  for(int i=0, j=0; (i+j)<=SplinePoints; ){

    // Update the CurrentPath
    CurrentPath = sp[0] + (i+1)*DeltaS;

    if( CurrentPath<sp[j] && (i+j)<SplinePoints ){
      // Get the value of (x,y) coordinates for the CurrentPath
      PointOnCurve = Spline(CurrentPath);

      // Write the coordinates of PointOnCurve
      output_file << PointOnCurve.x << "\t" << PointOnCurve.y << "\n";
      ++i;
    } else {
      output_file << Xp[j].x << "\t" << Xp[j].y << "\n";
      ++j;
    }
  }
}

/*!
 * This routine returns the pathlength of vector location X 
 * on spline S.                                             
 */
double Spline2D_HO::getS(const Vector2D &X) const {
    int i, subinterval_found, icount;
    double s1, s2, s3, s3_old, x1, x2, x3, y1, y2, y3;

    /* Determine the subinterval of the spline
       containing the point, X, of interest. */

    subinterval_found = 0;
    i = 0;
    while (subinterval_found == 0 && i < np-1) {
      x1=Xp[i].x;
      x2=Xp[i+1].x;
      y1=Xp[i].y;
      y2=Xp[i+1].y;
      s1=sp[i];
      s2=sp[i+1];
////////////////////////////////////////////////////////////////////////
      if (abs(X-Vector2D(x1,y1)) < NANO) return s1;
      if (abs(X-Vector2D(x2,y2)) < NANO) return s2;
////////////////////////////////////////////////////////////////////////
      if ( ((X.x-x1)*(X.x-x2) <= ZERO) &&
           ((X.y-y1)*(X.y-y2) <= ZERO) ) {
        subinterval_found = 1;
      } else {
        i += 1;
      } /* endif */
    } /* endwhile */

    if (!subinterval_found) {
       subinterval_found = 0;
       i = 0;
       while (subinterval_found == 0 && i < np-1) {
          x1=Xp[i].x;
          x2=Xp[i+1].x;
          y1=Xp[i].y;
          y2=Xp[i+1].y;
          s1=sp[i];
          s2=sp[i+1];
          if ( ((X.x-x1)*(X.x-x2) <= sqr(TOLER*max(fabs(X.x),TOLER))) &&
               ((X.y-y1)*(X.y-y2) <= sqr(TOLER*max(fabs(X.y),TOLER))) ) {
            subinterval_found = 1;
          } else {
            i += 1;
          } /* endif */
       } /* endwhile */

       if (!subinterval_found) {
          subinterval_found = 0;
          i = 0;
          while (subinterval_found == 0 && i < np-1) {
             x1=Xp[i].x;
             x2=Xp[i+1].x;
             y1=Xp[i].y;
             y2=Xp[i+1].y;
             s1=sp[i];
             s2=sp[i+1];
             if ( ((X.x-x1)*(X.x-x2) <= sqr(TOLER*max(fabs(X.x),TOLER))) ||
                  ((X.y-y1)*(X.y-y2) <= sqr(TOLER*max(fabs(X.y),TOLER))) )  {
               subinterval_found = 1;
             } else {
               i += 1;
             } /* endif */
          } /* endwhile */
          if (!subinterval_found) {
             return(sp[0]);
          } /* endif */
       } /* endif */
    } /* endif */

    /* Use method of false position to find value of S
       that will yield vector position X. */

    if ( (fabs(y2-y1) > fabs(x2-x1)) &&
         ((X.y-y1)*(X.y-y2) <= sqr(TOLER*X.y)) ) { 
       s3_old = s1;
       s3 = s2;
       icount = 0;
       while ( fabs(s3-s3_old) > TOLER*TOLER*s3 ) {
         s3_old = s3;
         s3 = (s1*(y2-X.y) - s2*(y1-X.y))/(y2-y1);
         y3=Spline(s3).y;

         if ((y1-X.y)*(y3-X.y) < ZERO) {
	   y2=y3;
   	   s2=s3;
         } else {
	   y1=y3;
	   s1=s3;
         } /* endif */

         icount += 1;
         if (icount > 1000) break;
       } /* endwhile */
    } else {
       s3_old = s1;
       s3 = s2;
       icount = 0;
       while ( fabs(s3-s3_old) > TOLER*TOLER*s3 ) {
         s3_old = s3;
         s3 = (s1*(x2-X.x) - s2*(x1-X.x))/(x2-x1);
         x3=Spline(s3).x;

         if ((x1-X.x)*(x3-X.x) < ZERO) {
	   x2=x3;
   	   s2=s3;
         } else {
	   x1=x3;
	   s1=s3;
         } /* endif */

         icount += 1;
         if (icount > 1000) break;
       } /* endwhile */
    } /* endif */

    /* Return pathlength. */

    return(s3);
  
}

/*!
 * This routine returns the splined value of Y given X.      
 *                                                           
 * Returns a VECTOR2D LINKED LIST containing all the vectors 
 * on the spline with the desired X value.                   
 */
LinkedList<Vector2D> Spline2D_HO::getY(const double &x) const {

  int i, iNext, found, icount;
  double s1, s2, s3, s3_old, x1, x2, x3;
  LinkedList<Vector2D> LL;

  iNext=0;

  while (0<1) {
		  
    /* Do incremental search to find interval. */

    found=0;
    for (i=iNext; i<=np-2; ++i) {
      x1=Xp[i].x;
      x2=Xp[i+1].x;
      s1=sp[i];
      s2=sp[i+1];
      if ((x>=x1 && x<=x2) || (x<=x1 && x>=x2)) {
	iNext=i+1;
	found=1;
	break;
      } /* endif */
    } /* endfor */
    
    if (found==0) break;

    /* Use method of false position to find 
       value of 's' that will give 'x'. */

    s3_old = s1;
    s3 = s2;
    icount = 0;
    while ( fabs(s3-s3_old) > TOLER*TOLER*s3) {
      s3_old = s3;
      s3 = (s1*(x2-x) - s2*(x1-x))/(x2-x1);
      x3=Spline(s3).x;
     
      if ((x1-x)*(x3-x) < ZERO) {
	x2=x3;
	s2=s3;
      } else {
	x1=x3;
	s1=s3;
      } /* endif */
      icount += 1;
      if (icount > 1000) break;
    } /* endwhile */
    
      /* Check for duplicate. */

    if (LL.np==0) {
      LL.add( Vector2D(x,Spline(s3).y) );
    } else if (LL[LL.np-1].y!=Spline(s3).y) {
      LL.add( Vector2D(x,Spline(s3).y) );
    } /* endif */ 

  } /* endwhile */
  
  return(LL);
}

/*!
 * This routine returns the splined value of X given Y.      
 *                                                           
 * Returns a VECTOR2D LINKED LIST containing all the vectors 
 * on the spline with the desired Y value.                   
 *                                                           
 */
LinkedList<Vector2D> Spline2D_HO::getX(const double &y) const {

  int i, iNext, found, icount;
  double s1, s2, s3, s3_old, y1, y2, y3;
  LinkedList<Vector2D> LL;
  
  iNext=0;
  
  while (0<1) {
    
    /* Do incremental search to find interval. */
    
    found=0;
    for (i=iNext; i<=np-2; ++i) {
      y1=Xp[i].y;
      y2=Xp[i+1].y;
      s1=sp[i];
      s2=sp[i+1];
      if ((y>=y1 && y<=y2) || (y<=y1 && y>=y2)) {
	iNext=i+1;
	found=1;
	break;
      } /* endif */
    } /* endfor */
    
    if (found==0) break;
    
    /* Use method of false position to find 
       value of 's' that will give 'y'. */
    
    s3_old = s1;
    s3 = s2;
    icount = 0;
    while ( fabs(s3-s3_old) > TOLER*TOLER*s3 ) {
      s3_old = s3;
      s3 = (s1*(y2-y) - s2*(y1-y))/(y2-y1);
      y3=Spline(s3).y;
      
      if ((y1-y)*(y3-y) < ZERO) {
	y2=y3;
	s2=s3;
      } else {
	y1=y3;
	s1=s3;
      } /* endif */
      icount += 1;
      if (icount > 1000) break;
    }/* endwhile */
    
    /* Check for duplicates. */
    
    if (LL.np==0) {
      LL.add( Vector2D(Spline(s3).x, y) ); 
    } else if(LL[LL.np-1].x!=Spline(s3).x) {
      LL.add( Vector2D(Spline(s3).x, y) );
    } /* endif */
    
  } /* endwhile */
  
  return(LL);
}


/*!
 * This routine returns the boundary condition of vector 
 * location X on spline S.                               
 */
int Spline2D_HO::getBCtype(const Vector2D &X) const {

  int bctype, np, subinterval_found;
  Vector2D X1, X2;

  subinterval_found = 0;
  np = 0;
  while (!subinterval_found && np < np-1) {
    X1 = Xp[np];
    X2 = Xp[np+1];
    if (((X.x-X1.x)*(X.x-X2.x) <= ZERO) && ((X.y-X1.y)*(X.y-X2.y) <= ZERO)) {
      subinterval_found = 1;
      //       if (abs(X-X1) < abs(X-X2)) bctype = bc[np];
      //       else bctype = bc[np+1];
      bctype = max(bc[np],bc[np+1]);
    } else {
      np++;
    }
  }
  if (!subinterval_found) {
    subinterval_found = 0;
    np = 0;
    while (!subinterval_found && np < np-2) {
      X1 = Xp[np];
      X2 = Xp[np+2];
      if (((X.x-X1.x)*(X.x-X2.x) <= ZERO) && ((X.y-X1.y)*(X.y-X2.y) <= ZERO)) {
	subinterval_found = 1;
	bctype = bc[np+1];
      } else {
	np++;
      }
    }
  }
  if (!subinterval_found) {
    np = 0;
    for (int n = 0; n < np; n++)
      if (abs(Xp[n]-X) < abs(Xp[np]-X)) np = n;
    bctype = bc[np];
  }

  // Return the boundary condition type.
  return bctype;
  
}

/*!
 * This routine returns the nearest splined value of y given x. 
 */
Vector2D Spline2D_HO::getminY(const Vector2D &X) const {

  int iNext, found, icount;
  double s1, s2, s3, s3_old, x1, x2, x3;
  Vector2D Xpoint;

  iNext = 0;
  Xpoint = Vector2D(MILLION,MILLION);

  while (0 < 1) {
		  
    // Do incremental search to find interval.
    found = 0;
    for (int i = iNext; i < np-1; i++) {
      x1 = Xp[i].x;
      x2 = Xp[i+1].x;
      s1 = sp[i];
      s2 = sp[i+1];
      if ((X.x >= x1 && X.x <= x2) || (X.x <= x1 && X.x >= x2)) {
	iNext = i+1;
	found = 1;
	break;
      }
    }

    if (found==0) break;

    // Use method of false position to find value of 's' that will give 'x'.
    s3_old = s1;
    s3 = s2;
    icount = 0;
    while (fabs(s3-s3_old) > TOLER*TOLER*s3) {
      s3_old = s3;
      s3 = (s1*(x2-X.x) - s2*(x1-X.x))/(x2-x1+NANO);
      x3 = Spline(s3).x;
      if ((x1-X.x)*(x3-X.x) < ZERO) {
	x2 = x3;
	s2 = s3;
      } else {
	x1 = x3;
	s1 = s3;
      }
      icount++;
      if (icount > 1000) break;
    }

    // Choose nearest point.
    if (abs(X-Vector2D(X.x,Spline(s3).y)) < abs(X-Xpoint)) 
      Xpoint = Vector2D(X.x,Spline(s3).y);

  }

  // Return the nearest spline point.
  return Xpoint;
}

/*!
 * This routine returns the nearest splined value of x given y.
 */
Vector2D Spline2D_HO::getminX(const Vector2D &X) const {

  int iNext, found, icount;
  double s1, s2, s3, s3_old, y1, y2, y3;
  Vector2D Xpoint;

  iNext = 0;
  Xpoint = Vector2D(MILLION,MILLION);

  while (0 < 1) {
		  
    // Do incremental search to find interval.
    found = 0;
    for (int i = iNext; i < np-1; i++) {
      y1 = Xp[i].y;
      y2 = Xp[i+1].y;
      s1 = sp[i];
      s2 = sp[i+1];
      if ((X.y >= y1 && X.y <= y2) || (X.y <= y1 && X.y >= y2)) {
	iNext = i+1;
	found = 1;
	break;
      }
    }
    
    if (found == 0) break;

    // Use method of false position to find value of 's' that will give 'y'.
    s3_old = s1;
    s3 = s2;
    icount = 0;
    while (fabs(s3-s3_old) > TOLER*TOLER*s3) {
      s3_old = s3;
      s3 = (s1*(y2-X.y) - s2*(y1-X.y))/(y2-y1+NANO);
      y3 = Spline(s3).y;
      if ((y1-X.y)*(y3-X.y) < ZERO) {
	y2 = y3;
	s2 = s3;
      } else {
	y1 = y3;
	s1 = s3;
      }
      icount++;
      if (icount > 1000) break;
    }

    // Choose nearest point.
    if (abs(X-Vector2D(Spline(s3).x,X.y)) < abs(X-Xpoint)) 
      Xpoint = Vector2D(Spline(s3).x,X.y);

  }

  // Return the nearest spline point.
  return Xpoint;
}

/*!
 * Get the normal at the specified location
 */
Vector2D Spline2D_HO::getnormal(const Vector2D &X) const {

  int np, subinterval_found;
  Vector2D X1, X2;

  subinterval_found = 0;
  np = 0;
  while (!subinterval_found && np < np-1) {
    X1 = Xp[np];
    X2 = Xp[np+1];
    if ((X.x-X1.x)*(X.x-X2.x) <= ZERO) {
      subinterval_found = 1;
    } else {
      np++;
    }
  }
  if (!subinterval_found) {
    subinterval_found = 0;
    np = 0;
    while (!subinterval_found && np < np-2) {
      X1 = Xp[np];
      X2 = Xp[np+2];
      if ((X.x-X1.x)*(X.x-X2.x) <= ZERO) {
	subinterval_found = 1;
      } else {
	np++;
      }
    }
  }
  if (!subinterval_found) {
    np = 0;
    for (int n = 0; n < np; n++) {
      if (abs(Xp[n]-X) < abs(Xp[np]-X)) {
	np = n;
	if (np == 0) { X1 = Xp[0]; X2 = Xp[1]; }
	else if (np == np) { X1 = Xp[np-2]; X2 = Xp[np-1]; }
	else { X1 = Xp[np-1]; X2 = Xp[np+1]; }
      }
    }
  }

  // Return the unit normal.
  return Vector2D((X2.y - X1.y),-(X2.x - X1.x))/abs(X2 - X1);
  
}


/*!
 * This routine returns the boundary condition type 
 * given the path length s along the spline and the 
 * the set of np boundary types (bc) that have been 
 * defined for each spline point.                   
 */
int Spline2D_HO::BCtype(const double &s) const {

  int i, il, ir, subinterval_found;

  /* Determine the subinterval of the spline
     containing the point (position) of interest. */
  
  if (s <= sp[0]) return(bc[0]);
  if (s >= sp[np-1]) return(bc[np-1]);      

  /* Determine the subinterval of the spline
     which contains the point (position) of interest. */
  find_subinterval(s,il,ir);

  /* Determine the boundary condition type. */

  if (bc[il] == bc[ir]) {
    return(bc[il]);
  } else {
    return(bc[il]);
  } /* endif */
 
}


/*!
 * Broadcasts a spline to all processors involved in  
 * the calculation from the primary processor using   
 * the MPI broadcast routine.                         
 */
void Spline2D_HO::Broadcast_Spline(void) {

#ifdef _MPI_VERSION
    int i, npts, buffer_size, i_buffer_size;
    int *i_buffer;
    double *buffer;
 
    /* Broadcast the number of spline points. */

    if (CFFC_Primary_MPI_Processor()) {
      npts = np;
    } /* endif */

    MPI::COMM_WORLD.Bcast(&npts, 1, MPI::INT, 0);

    /* On non-primary MPI processors, allocate (re-allocate) 
       memory for the spline as necessary. */

    if (!CFFC_Primary_MPI_Processor()) {
      if (npts >= 2){
	allocate(npts);
      } else {
	// deallocate the current spline if the broadcast spline has no memory allocated
	deallocate();
      } /* endif */
    } /* endif */

    /* Broadcast the spline type. */

    MPI::COMM_WORLD.Bcast(&(type), 1, MPI::INT, 0);

    /* Broadcast the flux calculation method associated with this spline. */

    MPI::COMM_WORLD.Bcast(&(FluxMethod), 1, MPI::INT, 0);

    /* Broadcast the maximum number of solid bodies. */
    MPI::COMM_WORLD.Bcast(&CounterSolidBodies, 1, MPI::INT, 0);

    /* Broadcast the ID of the current spline. */
    MPI::COMM_WORLD.Bcast(&bodyOwner, 1, MPI::INT, 0);


    /* Broadcast the the spline coordinates, pathlength, 
       point type, and boundary condition information. */

    if (npts >= 2) {
       buffer = new double[3*npts];
       i_buffer = new int[2*npts];

       if (CFFC_Primary_MPI_Processor()) {
          buffer_size = 0;
          i_buffer_size = 0;
          for ( i = 0; i <= np-1; ++i ) {
              buffer[buffer_size] = Xp[i].x;
              buffer[buffer_size+1] = Xp[i].y;
              buffer[buffer_size+2] = sp[i];
              i_buffer[i_buffer_size] = tp[i];
              i_buffer[i_buffer_size+1] = bc[i];
              buffer_size = buffer_size + 3;
              i_buffer_size = i_buffer_size + 2;
          } /* endfor */
       } /* endif */

       buffer_size = 3*npts;
       i_buffer_size = 2*npts;
       MPI::COMM_WORLD.Bcast(buffer, buffer_size, MPI::DOUBLE, 0);
       MPI::COMM_WORLD.Bcast(i_buffer, i_buffer_size, MPI::INT, 0);

       if (!CFFC_Primary_MPI_Processor()) {
          buffer_size = 0;
          i_buffer_size = 0;
          for ( i = 0; i <= np-1; ++i ) {
              Xp[i].x = buffer[buffer_size];
              Xp[i].y = buffer[buffer_size+1];
              sp[i] = buffer[buffer_size+2];
              tp[i] = i_buffer[i_buffer_size];
              bc[i] = i_buffer[i_buffer_size+1];
              buffer_size = buffer_size + 3;
              i_buffer_size = i_buffer_size + 2;
          } /* endfor */
       } /* endif */

       delete []buffer; 
       buffer = NULL;
       delete []i_buffer; 
       i_buffer = NULL;
    } /* endif */
#endif
}

#ifdef _MPI_VERSION
/*!
 * Broadcasts a spline to all processors associated   
 * with the specified communicator from the specified 
 * processor using the MPI broadcast routine.         
 */
void Spline2D_HO::Broadcast_Spline(MPI::Intracomm &Communicator, 
				   const int Source_CPU) {

    int Source_Rank = 0;
    int i, npts, buffer_size, i_buffer_size;
    int *i_buffer;
    double *buffer;
 
    /* Broadcast the number of spline points. */

    if (CFFC_MPI::This_Processor_Number == Source_CPU) {
      npts = np;
    } /* endif */

    Communicator.Bcast(&npts, 1, MPI::INT, Source_Rank);

    /* On non-source MPI processors, allocate (re-allocate) 
       memory for the spline as necessary. */

    if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
      if (npts >= 2) {
	allocate(npts);
      } else {
	// deallocate the current spline if the broadcast spline has no memory allocated
	deallocate();
      } /* endif */
    } /* endif */

    /* Broadcast the spline type. */

    Communicator.Bcast(&(type), 1, MPI::INT, Source_Rank);

    /* Broadcast the flux calculation method associated with this spline. */

    Communicator.Bcast(&(FluxMethod), 1, MPI::INT, Source_Rank);

    /* Broadcast the maximum number of solid bodies. */
    Communicator.Bcast(&CounterSolidBodies, 1, MPI::INT, Source_Rank);

    /* Broadcast the ID of the current spline. */
    Communicator.Bcast(&bodyOwner, 1, MPI::INT, Source_Rank);


    /* Broadcast the the spline coordinates, pathlength, 
       point type, and boundary condition information. */

    if (npts >= 2) {
       buffer = new double[3*npts];
       i_buffer = new int[2*npts];

       if (CFFC_MPI::This_Processor_Number == Source_CPU) {
          buffer_size = 0;
          i_buffer_size = 0;
          for ( i = 0; i <= np-1; ++i ) {
              buffer[buffer_size] = Xp[i].x;
              buffer[buffer_size+1] = Xp[i].y;
              buffer[buffer_size+2] = sp[i];
              i_buffer[i_buffer_size] = tp[i];
              i_buffer[i_buffer_size+1] = bc[i];
              buffer_size = buffer_size + 3;
              i_buffer_size = i_buffer_size + 2;
          } /* endfor */
       } /* endif */

       buffer_size = 3*npts;
       i_buffer_size = 2*npts;
       Communicator.Bcast(buffer, buffer_size, MPI::DOUBLE, Source_Rank);
       Communicator.Bcast(i_buffer, i_buffer_size, MPI::INT, Source_Rank);

       if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
          buffer_size = 0;
          i_buffer_size = 0;
          for ( i = 0; i <= np-1; ++i ) {
              Xp[i].x = buffer[buffer_size];
              Xp[i].y = buffer[buffer_size+1];
              sp[i] = buffer[buffer_size+2];
              tp[i] = i_buffer[i_buffer_size];
              bc[i] = i_buffer[i_buffer_size+1];
              buffer_size = buffer_size + 3;
              i_buffer_size = i_buffer_size + 2;
          } /* endfor */
       } /* endif */

       delete []buffer; 
       buffer = NULL;
       delete []i_buffer; 
       i_buffer = NULL;
    } /* endif */

}
#endif


/*!
 * Ensure that all processors involved in the calculation
 * have correctly setup the maximum number of solid bodies.
 */
void Spline2D_HO::Broadcast_Maximum_Number_Of_SolidBodies(void){
#ifdef _MPI_VERSION  
  CounterSolidBodies = CFFC_Maximum_MPI(CounterSolidBodies);
#endif
}

/*!
 * Apply a mirror reflection about the y=0 axis and
 * recompute the positions of all of the points in 
 * the spline accordingly.                         
 */
void Spline2D_HO::Reflect_Spline(void) {

  int i;
  Vector2D X;
 
  /* Apply a mirror reflection about the y=0 axis. */

  for ( i = 0; i <= np-1; ++i ) {
    X.x = Xp[i].x;
    X.y = - Xp[i].y;
    Xp[i] = X;
  } /* endfor */

}

/*!
 * Reverses the order of the spline points.
 */
void Spline2D_HO::Reverse_Spline(void) {

   int i, tp_current, bc_current;
   double sp_current;
   Vector2D Xp_current;
   
   /* Reverse the order of the spline points. */
   for ( i = 0; i <= (np-1)/2; ++i ) {
      Xp_current = Xp[i];  Xp[i] = Xp[np-1-i];  Xp[np-1-i] = Xp_current;
      tp_current = tp[i];  tp[i] = tp[np-1-i];  tp[np-1-i] = tp_current;
      bc_current = bc[i];  bc[i] = bc[np-1-i];  bc[np-1-i] = bc_current;
   } /* endfor */

   // redetermine the path lengths
   pathlength();
}

/*!
 * This routine calculates and returns a 2D spline 
 * representing a straight line between two points.
 */
void Spline2D_HO::Create_Spline_Line(const Vector2D &V1,
				     const Vector2D &V2,
				     const int Number_of_Spline_Points) {

    int i;
    Vector2D dX;

    /* Allocate memory for the straight line spline. */
    allocate(Number_of_Spline_Points);

    /* Set the spline type. */

    settype(SPLINE2D_LINEAR);

    /* Compute the locations of the spline points on the 
       straight line between points V1 and V2. */

    dX = (V2-V1)/double(Number_of_Spline_Points-1);

    for (i = 0; i <= Number_of_Spline_Points-1; i++) {
        Xp[i] = V1 + double(i)*dX;

        if (i == 0 || i == Number_of_Spline_Points - 1) {
           tp[i] = SPLINE2D_POINT_SHARP_CORNER;
        } else {
           tp[i] = SPLINE2D_POINT_NORMAL;
        } /* endif */

        bc[i] = BC_NONE;
    } /* endfor */

    /* Calculate the spline pathlengths. */

    pathlength();

}

/*!
 * This routine calculates and returns a 2D spline      
 * representing a straight line between two points.
 * The start and end points are determined based on 
 * the polar coordinates (Inner_Radius,Theta) and (Outer_Radius,Theta).     
 * The angle Theta is considered to be given in degrees.                                                       
 */
void Spline2D_HO::Create_Spline_Line_Polar_Coordinates(const double &Inner_Radius,
						       const double &Outer_Radius,
						       const double &Theta,
						       const int Number_of_Spline_Points){
  
  Vector2D V1, V2;		// the Cartesian coordinates of the start and end points of the line.
  
  // Set V1
  V1.setWithPolarCoord(Inner_Radius,Theta);

  // Set V2
  V2.setWithPolarCoord(Outer_Radius,Theta);

  // Generate the line
  Create_Spline_Line(V1,V2,Number_of_Spline_Points);
}

/*!
 * This routine calculates and returns a 2D spline
 * representing a segment of a circular arc.      
 */
void Spline2D_HO::Create_Spline_Circular_Arc(const Vector2D &Origin,
					     const double &Radius,
					     const double &Angle1,
					     const double &Angle2,
					     const int Number_of_Spline_Points) {

  int i;
  double theta;

  /* Allocate memory for the circular arc spline. */
  allocate(Number_of_Spline_Points);

  /* Set the spline type. */

  settype(SPLINE2D_QUINTIC);

  /* Compute the locations of the spline points on the 
     circular arc. */

  for (i = 0; i <= Number_of_Spline_Points-1; i++) {
    theta = Angle1+
      (Angle2-Angle1)*double(i)/double(Number_of_Spline_Points-1);
    theta = TWO*PI*theta/360.0;

    Xp[i].x = Radius*cos(theta);
    Xp[i].y = Radius*sin(theta);

    Xp[i] += Origin;        

    if (i == 0 || i == Number_of_Spline_Points - 1) {
      tp[i] = SPLINE2D_POINT_SHARP_CORNER;
    } else {
      tp[i] = SPLINE2D_POINT_NORMAL;
    } /* endif */

    bc[i] = BC_NONE;
  } /* endfor */

    /* Ensure that arc closes on itself if it is a
       complete circle. */

  if ( fabs(fabs(Angle2-Angle1)-360.00) < TOLER ) {
    Xp[np-1] = Xp[0];
    tp[0] = SPLINE2D_POINT_NORMAL;
    tp[np-1] = SPLINE2D_POINT_NORMAL;
  } /* endif */

    /* Calculate the spline pathlengths. */

  pathlength();

}

/*!
 * This routine calculates and returns a 2D spline
 * representing a segment of an ellipsoidal arc.  
 */
void Spline2D_HO::Create_Spline_Ellipsoidal_Arc(const Vector2D &Origin,
						const double &A,
						const double &B,
						const double &Angle1,
						const double &Angle2,
						const int Number_of_Spline_Points) {

    int i;
    double theta;

    /* Allocate memory for the ellipsoidal arc spline. */

    allocate(Number_of_Spline_Points);
    
    /* Set the spline type. */

    settype(SPLINE2D_QUINTIC);
    
    /* Compute the locations of the spline points on the 
       ellipsoidal arc. */
    
    for (i = 0; i <= Number_of_Spline_Points-1; i++) {
      theta = Angle1 + (Angle2-Angle1)*double(i)/double(Number_of_Spline_Points-1);
      theta = TWO*PI*theta/360.0;
      
      Xp[i].x = A*cos(theta);
      Xp[i].y = B*sin(theta);
      
      Xp[i] += Origin;        
      
      if (i == 0 || i == Number_of_Spline_Points - 1) {
	tp[i] = SPLINE2D_POINT_SHARP_CORNER;
      } else {
	tp[i] = SPLINE2D_POINT_NORMAL;
      } /* endif */
      
      bc[i] = BC_NONE;
    } /* endfor */
    
    /* Ensure that arc closes on itself if it is a
       complete ellipse. */
    
    if ( fabs(fabs(Angle2-Angle1)-360.00) < TOLER ) {
      Xp[np-1] = Xp[0];
      tp[0] = SPLINE2D_POINT_NORMAL;
      tp[np-1] = SPLINE2D_POINT_NORMAL;
    } /* endif */
    
    /* Calculate the spline pathlengths. */
    
    pathlength();
    
}

/*!
 * This routine calculates and returns a 2D spline
 * representing a segment of a sinusoidal line.   
 */
void Spline2D_HO::Create_Spline_Sinusoidal_Line(const Vector2D &Origin,
						const double &Radius,
						const double &Length,
						const double &Angle1,
						const double &Angle2,
						const int Number_of_Spline_Points) {

  int i;
  double theta;
  
  /* Allocate memory for the circular arc spline. */
  allocate(Number_of_Spline_Points);

  /* Set the spline type. */

  settype(SPLINE2D_QUINTIC);

  /* Compute the locations of the spline points on the 
     circular arc. */

  for (i = 0; i <= Number_of_Spline_Points-1; i++) {
    theta = Angle1+
      (Angle2-Angle1)*double(i)/double(Number_of_Spline_Points-1);
    theta = TWO*PI*theta/360.0;

    Xp[i].x = Length*theta/(TWO*PI*(Angle2-Angle1)/360.0);
    Xp[i].y = Radius*sin(theta);

    Xp[i] += Origin;        

    if (i == 0 || i == Number_of_Spline_Points - 1) {
      tp[i] = SPLINE2D_POINT_SHARP_CORNER;
    } else {
      tp[i] = SPLINE2D_POINT_NORMAL;
    } /* endif */

    bc[i] = BC_NONE;
  } /* endfor */

    /* Calculate the spline pathlengths. */

  pathlength();

}

/*!
 * \verbatim
 *
 * This routine calculates and returns a 2D spline for  
 * both NACA 4- and 5-digit aerofoil sections.  The     
 * algorithm for calculating the splines is straight    
 * out of Abbott and von Doenhoff, with a little        
 * Reigels mixed in.                                    
 *                                                      
 * Basically, a NACA aerofoil is composed of a camber   
 * line and a thickness distribution.  The thickness    
 * distribution is a single equation, while the camber  
 * is usually two joined quadratics.                    
 *                                                      
 * The equations for the upper and lower coordinates    
 * are:                                                 
 *                                                      
 *  x(upper) =                                          
 *    x - yt*sin(theta)  y(upper) = yc + yt*cos(theta)  
 *  x(lower) =                                          
 *    x + yt*sin(theta)  y(lower) = yc - yt*cos(theta)  
 *                                                      
 * where tan(theta) = d(yc)/dx.  In these equations, yc 
 * is the camber line, yt is the thickness distribution.
 * A common approximation (small-angle) is to assume    
 * theta is small, so that sin(theta) is approx. 0 and  
 * cos(theta) is approx. 1.  The equations become:      
 *                                                      
 * 	  x(upper) = x     y(upper) = yc + yt           
 * 	  x(lower) = x     y(lower) = yc - yt           
 *                                                      
 * For 4-digit aerofoils, the camber lines and          
 * thickness:                                           
 *                                                      
 * (yc/c) = (f/c)*(1/(x1^2))*(2*x1*(x/c) - (x/c)^2)     
 *              for 0<=(x/c)<=x1                        
 *                                                      
 * and                                                  
 *                                                      
 * (yc/c) = (f/c)*(1/(1-x1)^2)*((1-2x1)+2x1*(x/c)-      
 *          (x/c)^2)                                    
 *              for x1<=(x/c)<=1 with x1=(xf/c)         
 *                                                      
 * (yt/c) = 5t*(0.29690*x^0.5 - 0.12600X - 0.35160*x^2  
 *         + 0.28430*x^3 - 0.10150*x^4)                 
 *                                                      
 * where t = thickness/chord,                           
 * x = position along x-axis,                           
 * xf = position of maximum camber,                     
 * f = maximum camber.                                  
 *                                                      
 * For 5-digit aerofoils the thickness distribution is  
 * the same, only the camber line is changed.  There are
 * two types, based on the third digit.  The majority   
 * of 5-digit aerofoils are 'type 0' (i.e.,             
 * NACA 23015) hand have a camber line given by:        
 *                                                      
 * (yc/c) = (k1/6)*(x^3 - 3*x1*x^2 + x1^2*(3-x1)*x)     
 *               for 0<=x<=x1                           
 *                                                      
 * and                                                  
 *                                                      
 * (yc/c) = (k1*x1^3 / 6) * (1 - x)                     
 *        for x1<=x<=1 with x1 = (x1/c) and x = (x/c),  
 * and x1 is related to xf, which is the position of    
 * maximum camber.
 *                                                       
 * The constants x1 and k1 are determined from the      
 * following table:                                     
 *                                                      
 * xf           0.05    0.10    0.15     0.20    0.25   
 * x1           0.0580  0.1260  0.2025   0.2900  0.3910 
 * (Cl* /(f/c)) 26.9    19.6    16.4     14.5    11.3   
 * (k1/Cl*)     1205    172.1   53.2     22.13   10.77  
 *                                                      
 * The 'type 1' aerofoil camber line is not provided    
 * here.                                                
 *                                                      
 * Breakdown of the NACA designations:                  
 *                                                      
 * In a 4-digit aerofoil, the first digit is the value  
 * of the maximum camber (in percent of the chord), the 
 * second digit is the position of the maximum camber   
 * from the leading edge in tenths of the chord, and    
 * the last two digits denote the maximum thickness of  
 * the aerofoil in percent.  For the NACA 2415 aerofoil,
 * the maximum camber is 2%, the position of the        
 * maximum camber is 0.4c, and the thickness is 15%.    
 * This example is displayed in below.                  
 *                                                      
 * NACA    2415                                         
 *                                                      
 * 24      mean line (24)                               
 * 2       maximum camber, (2%)                         
 * 4       10 * position of maximum camber, (0.4c)      
 * 15      thickness, (15%)                             
 *                                                      
 * The NACA 5-digit aerofoils are set up in a similar   
 * manner to the 4-digit aerofoils.  The primary        
 * difference is the use of a different camber line.    
 * In a 5-digit aerofoil, 1.5 times the first digit is  
 * the design lift coefficient in tenths, the second    
 * and third digits are one-half the distance from the  
 * leading edge to the location of maximum camber in    
 * percent of the chord, and the fourth and fifth       
 * digits are the thickness in percent of the chord.    
 * For example, a NACA 23015 aerofoil has a design lift 
 * coefficient of 0.3, has the maximum camber at 0.15c, 
 * and is 15% thick.  Additionally, the first three     
 * digits indicate the mean line used.  In this case,   
 * the mean line designation is 230.  The 5-digit       
 * aerofoils use the same thickness distribution as the 
 * 4-digit aerofoils.  This example is displayed below. 
 *                                                      
 * NACA    23015                                        
 *                                                      
 * 230     mean line (230)                              
 * 2       (design lift coefficient * 10) / 1.5, (0.3)  
 * 30      2 * (position of maximum camber),            
 *         (0.30 / 2 = 0.15c)                           
 * 0       type of camber line used                     
 * 15      thickness, (15%)                             
 *                                                      
 * References:                                          
 *                                                      
 * I. H. Abbot and von A. E. Doenhoff, "Theory of       
 * Wing Sections", Dover Publications, 1959.            
 *                                                      
 * F. W. Reigels, "Aerofoil Sections", Butterworth &    
 * Co., 1961.                                           
 *
 * \endverbatim                                                      
 */
void Spline2D_HO::Create_Spline_NACA_Aerofoil(char *NACA_Aerofoil_Type_ptr,
					      const double &Chord_Length,
					      const int i_Up_All_Low,
					      const int Number_of_Spline_Points) {

  int i, number_of_digits;
  char camber_max_indicator[2], x_camber_max_indicator[3], 
    thickness_max_indicator[3], design_lift_indicator[2];
  double camber_max, x_camber_max, thickness_max, design_lift,
    x1, k1, x_thickness_scaling_factor=1.008930411356;
  double theta, x, y, xt, yt, yc;

  /* Allocate memory for the NACA aerofoil spline. */
  allocate(Number_of_Spline_Points);

  /* Set the spline type. */

  settype(SPLINE2D_QUINTIC);

  /* Determine the number of digits in the NACA Aerofoil
     type designator. */

  number_of_digits = strlen(NACA_Aerofoil_Type_ptr);

  /* Determine the aerofoil chamber line and thickness 
     distribution parameters. */

  if (number_of_digits == 4) {
    camber_max_indicator[0] = NACA_Aerofoil_Type_ptr[0];
    camber_max_indicator[1] = '\n';
    camber_max = double(atof(camber_max_indicator)/HUNDRED);
    x_camber_max_indicator[0] = NACA_Aerofoil_Type_ptr[1];
    x_camber_max_indicator[1] = '\n';
    x_camber_max = double(atof(x_camber_max_indicator)/TEN);
    thickness_max_indicator[0] = NACA_Aerofoil_Type_ptr[2];
    thickness_max_indicator[1] = NACA_Aerofoil_Type_ptr[3];
    thickness_max_indicator[2] = '\n';
    thickness_max = double(atof(thickness_max_indicator)/HUNDRED);
  } else if (number_of_digits == 5) {
    design_lift_indicator[0] = NACA_Aerofoil_Type_ptr[0];
    design_lift_indicator[1] = '\n';
    design_lift = double(1.50*atof(design_lift_indicator)/TEN);
    x_camber_max_indicator[0] = NACA_Aerofoil_Type_ptr[1];
    x_camber_max_indicator[1] = NACA_Aerofoil_Type_ptr[2];
    x_camber_max_indicator[2] = '\n';
    x_camber_max = double(HALF*atof(x_camber_max_indicator)/HUNDRED);
    thickness_max_indicator[0] = NACA_Aerofoil_Type_ptr[3];
    thickness_max_indicator[1] = NACA_Aerofoil_Type_ptr[4];
    thickness_max_indicator[2] = '\n';
    thickness_max = double(atof(thickness_max_indicator)/HUNDRED);
    if (fabs(x_camber_max-0.05) < TOLER) {
      x1 = 0.0580;
      k1 = 1205.0*design_lift;
    } else if (fabs(x_camber_max-0.10) < TOLER) {
      x1 = 0.1260;
      k1 = 172.1*design_lift;
    } else if (fabs(x_camber_max-0.15) < TOLER) {
      x1 = 0.2025;
      k1 = 53.2*design_lift;
    } else if (fabs(x_camber_max-0.20) < TOLER) {
      x1 = 0.2900;
      k1 = 22.13*design_lift;
    } else if (fabs(x_camber_max-0.25) < TOLER) {
      x1 = 0.3910;
      k1 = 10.77*design_lift;
    } else {
      x1 = 0.3910;
      k1 = 10.77*design_lift;
    } /* endif */
  } else {
  } /* end if */

    /* Compute the locations of the spline points on the surface
       of the aerofoil. */

  for (i = 0; i <= Number_of_Spline_Points-1; i++) {
    // Calculate theta.
    if ( i_Up_All_Low == 0 ) {
      theta = TWO*PI*double(i)/double(Number_of_Spline_Points-1);
    } else if ( i_Up_All_Low > 0 ) {
      theta = PI + PI*double(i)/double(Number_of_Spline_Points-1);
    } else {
      theta = PI*double(i)/double(Number_of_Spline_Points-1);
    } /* endif */

    // Calculate x location.
    x = HALF*(ONE+cos(theta));

    // Determine the y position of the cord line, yc.
    if (number_of_digits == 4) {
      if (x < x_camber_max && x_camber_max != ZERO) {
	yc = camber_max*x*(TWO*x_camber_max-x)/sqr(x_camber_max);
      } else {
	yc = camber_max*(ONE-TWO*x_camber_max+x*(TWO*x_camber_max-x))/
	  sqr(ONE-x_camber_max);
      } /* endif */
    } else if (number_of_digits == 5) {
      if (x < x1) {
	yc = (k1/SIX)*(x*x*x - THREE*x1*x*x + 
		       x1*x1*(THREE-x1)*x); 
      } else {
	yc = (k1/SIX)*x1*x1*x1*(ONE-x); 
      } /* endif */
    } else {
    } /* end if */

    // Determine the local aerofoil thickness, yt.
    xt = x_thickness_scaling_factor*x;
    yt = FIVE*thickness_max*(0.29690*sqrt(xt) - 
			     0.12600*xt -
			     0.35160*xt*xt + 
			     0.28430*xt*xt*xt - 
			     0.10150*xt*xt*xt*xt);

    // Evaluate the coordinates of the points on the aerofoil surface.
    if ( theta >= PI ) {
      y = yc + yt;
    } else {
      y = yc - yt;
    } /* endif */

    // Ensure closure on trailing edge.
    if ( ( i_Up_All_Low == 0 && 
	   (i == 0 || i == Number_of_Spline_Points - 1) ) ||
	 ( i_Up_All_Low > 0 && i == Number_of_Spline_Points - 1 ) ||
	 ( i_Up_All_Low < 0 && i == 0 ) ) {
      x = ONE;
      y = ZERO; 
    } /* endif */

    // Scale the resulting aerofoil coordinates.
    x = Chord_Length*x;
    y = Chord_Length*y;

    // Assign values to the spline variables.
    Xp[i].x = x;
    Xp[i].y = y;

    if (i == 0 || i == Number_of_Spline_Points - 1) {
      tp[i] = SPLINE2D_POINT_SHARP_CORNER;
    } else {
      tp[i] = SPLINE2D_POINT_NORMAL;
    } /* endif */

    bc[i] = BC_NONE;
  } /* endfor */

    /* Calculate the spline pathlengths. */

  pathlength();

}

/*!
 * This routine calculates and returns a 2D spline     
 * representing the approximate position of the bow    
 * shock for supersonic flow over a circular cylinder. 
 */
void Spline2D_HO::Create_Spline_Bow_Shock(const double &Radius,
					  const double &Mach_Number,
					  const int i_Up_All_Low,
					  const int Number_of_Spline_Points) {

  int i;
  double theta, x, y, d, rs, b, y_shock_max;

  /* Allocate memory for the bow shock spline. */
  allocate(Number_of_Spline_Points);

  /* Set the spline type. */

  settype(SPLINE2D_QUINTIC);

  /* Compute the locations of the spline points which
     define the position of the bow shock. */

  d  = 0.386*Radius*exp(4.67/(Mach_Number*Mach_Number));
  rs = 1.386*Radius*exp(1.89/pow(Mach_Number-ONE, 0.75));
  b  = asin(ONE/Mach_Number);
  y_shock_max = (rs/tan(b))*sqrt(sqr(ONE+(Radius+d)*tan(b)*tan(b)/rs)-ONE);

  for (i = 0; i <= Number_of_Spline_Points-1; i++) {
    // Calculate theta.
    if ( i_Up_All_Low == 0 ) {
      theta = -HALF*PI + PI*double(i)/double(Number_of_Spline_Points-1);
    } else if ( i_Up_All_Low > 0 ) {
      theta = HALF*PI*double(i)/double(Number_of_Spline_Points-1);
    } else {
      theta = -HALF*PI + HALF*PI*double(i)/double(Number_of_Spline_Points-1);
    } /* endif */

    y = y_shock_max*sin(theta);
    x = -(Radius+d-(rs/(tan(b)*tan(b)))*(sqrt(ONE+sqr(y*tan(b)/rs))-ONE));

    Xp[i].x = x;
    Xp[i].y = y;

    if (i == 0 || i == Number_of_Spline_Points - 1) {
      tp[i] = SPLINE2D_POINT_SHARP_CORNER;
    } else {
      tp[i] = SPLINE2D_POINT_NORMAL;
    } /* endif */

    bc[i] = BC_NONE;
  } /* endfor */

    /* Calculate the spline pathlengths. */

  pathlength();

}

/*!
 * This routine calculates and returns a 2D spline    
 * representing the radius as a function of axial     
 * position for a duct with a smoothly varying change 
 * in cross-sectional area.                           
 */
void Spline2D_HO::Create_Spline_Area_Variation(const double &Xup,
					       const double &Xthroat,
					       const double &Xdown,
					       const double &Rup,
					       const double &Rthroat,
					       const double &Rdown,
					       const int &Nozzle_Type,
					       const int Number_of_Spline_Points) {

  int i;
  double x, radius, Aup, Athroat, Adown, drdx, r, beta;

  /* Allocate memory for the area variation spline. */
  allocate(Number_of_Spline_Points);

  /* Set the spline type. */

  settype(SPLINE2D_QUINTIC);

  /* Compute the locations of the spline points which
     define the radius of the duct as a function of the
     axial postion. */

  Aup = PI*Rup*Rup;
  Athroat = PI*Rthroat*Rthroat;
  Adown = PI*Rdown*Rdown;

  if (Nozzle_Type == NOZZLE_HYBRID_CONICAL_FUNCTION) {
    beta = 0.30;
    x = Xthroat + beta*(Xdown-Xthroat);
    r = Athroat*exp(0.50*log(Adown/Athroat)*(ONE-cos(PI*(x-Xthroat)/(Xdown-Xthroat))));
    r = sqrt(r/PI);
    drdx = 0.25*Athroat*log(Adown/Athroat)*sin(PI*(x-Xthroat)/(Xdown-Xthroat))*
      exp(0.5*log(Adown/Athroat)*(ONE-cos(PI*(x-Xthroat)/(Xdown-Xthroat))))/
      (sqrt((Athroat/PI)*exp(0.5*log(Adown/Athroat)*(ONE-cos(PI*(x-Xthroat)/(Xdown-Xthroat)))))*
       (Xdown-Xthroat));
  }

  for (i = 0; i <= Number_of_Spline_Points-1; i++) {
    x = (Xdown-Xup)*double(i)/double(Number_of_Spline_Points-1);

    if (Nozzle_Type == NOZZLE_CONICAL) {
      if (x <= Xthroat) {
	radius = Rup + (Rthroat-Rup)*(x-Xup)/(Xthroat-Xup);
      } else {
	radius = Rthroat + (Rdown-Rthroat)*(x-Xthroat)/(Xdown-Xthroat);
      }

    } else {
      if (x <= Xthroat) {
	radius = Aup*exp(HALF*log(Athroat/Aup)*
			 (ONE-cos(PI*(x-Xup)/(Xthroat-Xup))));
	radius = sqrt(radius/PI);
      } else {
	switch(Nozzle_Type) {
	case NOZZLE_GOTTLIEB_FUNCTION :
	  radius = Athroat*exp(HALF*log(Adown/Athroat)*
			       (ONE-cos(PI*(x-Xthroat)/(Xdown-Xthroat))));
	  radius = sqrt(radius/PI);
	  break;
	case NOZZLE_QUARTIC_FUNCTION :
	  radius = Athroat + (Adown-Athroat)*pow((x-Xthroat)/(Xdown-Xthroat),4.0);
	  radius = sqrt(radius/PI);
	  break;
	case NOZZLE_HYBRID_CONICAL_FUNCTION :
	  if (x <= Xthroat + beta*(Xdown-Xthroat)) {
	    radius = Athroat*exp(HALF*log(Adown/Athroat)*
				 (ONE-cos(PI*(x-Xthroat)/(Xdown-Xthroat))));
	    radius = sqrt(radius/PI);
	  } else {
	    radius = drdx*(x - Xthroat - beta*(Xdown-Xthroat)) + r;
	  }
	  break;
	};
      } /* endif */

    }

    Xp[i].x = x;
    Xp[i].y = radius;

    if (i == 0 || i == Number_of_Spline_Points - 1) {
      tp[i] = SPLINE2D_POINT_SHARP_CORNER;
    } else {
      tp[i] = SPLINE2D_POINT_NORMAL;
    } /* endif */

    bc[i] = BC_NONE;
  } /* endfor */

    /* Calculate the spline pathlengths. */

  pathlength();

}

void Spline2D_HO::Create_Spline_Converging_Nozzle(const double &Xup,
						  const double &Xthroat,
						  const double &Rup,
						  const double &Rthroat,
						  const int Number_of_Spline_Points) {

  double x, radius, Aup, Athroat;

  // Allocate memory for the area variation spline.
  allocate(Number_of_Spline_Points);

  // Set the spline type.
  settype(SPLINE2D_QUINTIC);

  // Compute the locations of the spline points which define the radius
  // of the duct as a function of the axial postion.
  Aup = PI*Rup*Rup;
  Athroat = PI*Rthroat*Rthroat;

  for (int i = 0; i < Number_of_Spline_Points; i++) {
    x = Xup + (Xthroat-Xup)*double(i)/double(Number_of_Spline_Points-1);

    radius = Aup*exp(HALF*log(Athroat/Aup)*(ONE-cos(PI*(x-Xup)/(Xthroat-Xup))));
    radius = sqrt(radius/PI);

    Xp[i].x = x;
    Xp[i].y = radius;

    if (i == 0 || i == Number_of_Spline_Points - 1) {
      tp[i] = SPLINE2D_POINT_SHARP_CORNER;
    } else {
      tp[i] = SPLINE2D_POINT_NORMAL;
    }

    bc[i] = BC_NONE;
  }

  // Calculate the spline pathlengths.
  pathlength();

}

void Spline2D_HO::Create_Spline_Diverging_Nozzle(const double &Xthroat,
						 const double &Xdown,
						 const double &Rthroat,
						 const double &Rdown,
						 const int Number_of_Spline_Points) {

  double x, radius, Athroat, Adown;

  // Allocate memory for the area variation spline.
  allocate(Number_of_Spline_Points);

  // Set the spline type.
  settype(SPLINE2D_QUINTIC);

  // Compute the locations of the spline points which define the radius
  // of the duct as a function of the axial postion.
  Athroat = PI*Rthroat*Rthroat;
  Adown = PI*Rdown*Rdown;

  for (int i = 0; i <= Number_of_Spline_Points-1; i++) {
    x = Xthroat + (Xdown-Xthroat)*double(i)/double(Number_of_Spline_Points-1);

    radius = Athroat*exp(0.50*log(Adown/Athroat)*(ONE-cos(PI*(x-Xthroat)/(Xdown-Xthroat))));
    radius = sqrt(radius/PI);

    Xp[i].x = x;
    Xp[i].y = radius;

    if (i == 0 || i == Number_of_Spline_Points - 1) {
      tp[i] = SPLINE2D_POINT_SHARP_CORNER;
    } else {
      tp[i] = SPLINE2D_POINT_NORMAL;
    }

    bc[i] = BC_NONE;
  }

  // Calculate the spline pathlengths.
  pathlength();

}

/*!
 * This routine calculates and returns a 2D spline 
 * representing a rectangle.                       
 */
void Spline2D_HO::Create_Spline_Rectangle(const Vector2D &Origin,
					  const double &Length,
					  const double &Width) {

  int i;

  /* Allocate memory for the rectangular spline. */
  allocate(5);

  /* Set the spline type. */

  settype(SPLINE2D_LINEAR);

  /* Compute the locations of the spline points on the 
     rectangle. */

  Xp[0].x = Origin.x - HALF*Width;
  Xp[0].y = Origin.y - HALF*Length;
  Xp[1].x = Origin.x + HALF*Width;
  Xp[1].y = Origin.y - HALF*Length;
  Xp[2].x = Origin.x + HALF*Width;
  Xp[2].y = Origin.y + HALF*Length;
  Xp[3].x = Origin.x - HALF*Width;
  Xp[3].y = Origin.y + HALF*Length;
  Xp[4]   = Xp[0];

  /* Set the point and boundary condition type. */

  for (i = 0; i < 5; i++) {
    tp[i] = SPLINE2D_POINT_SHARP_CORNER;
    bc[i] = BC_NONE;
  }

  /* Calculate the spline pathlengths. */

  pathlength();

}

/*!
 * This routine calculates and returns a 2D spline 
 * representing Zalesak's disk.                    
 */
void Spline2D_HO::Create_Spline_Zalesaks_Disk(const Vector2D &Origin,
					      const double &Radius) {

  int npts;
  int i;
  double theta, theta1, theta2;

  /* Set the number of spline points. */

  npts = 341 + 3;

  /* Set the angles. */

  theta1 = 10.0;
  theta2 = 350.0;

  /* Allocate memory for the circular arc spline. */
  allocate(npts);

  /* Set the spline type. */

  settype(SPLINE2D_LINEAR);

  /* Compute the locations of the spline points on the 
     spline and the point and boundary condition types. */

  for (i = 0; i < npts-3; i++) {
    theta = theta1 + (theta2 - theta1)*double(i)/double(npts-4);
    theta = TWO*PI*theta/360.0;
    Xp[i].x = Origin.x + Radius*cos(theta);
    Xp[i].y = Origin.y + Radius*sin(theta);
    bc[i]   = BC_NONE;
    tp[i]   = SPLINE2D_POINT_NORMAL;
    if (i == 0 || i == npts - 4)
      tp[i] = SPLINE2D_POINT_SHARP_CORNER;
  }
  Xp[npts-3] = Xp[npts-4] - Vector2D(1.5*Radius,ZERO);
  bc[npts-3] = BC_NONE;
  tp[npts-3] = SPLINE2D_POINT_SHARP_CORNER;
  Xp[npts-2] = Xp[0] - Vector2D(1.5*Radius,ZERO);
  bc[npts-2] = BC_NONE;
  tp[npts-2] = SPLINE2D_POINT_SHARP_CORNER;

  /* Last point must be the same as the first point. */

  Xp[npts-1] = Xp[0];
  bc[npts-1] = bc[0];
  tp[npts-1] = tp[0];

  /* Calculate the spline pathlengths. */

  pathlength();

  /* Rotate the spline. */
  Rotate_Spline(-HALF*PI);
}

/*!
 * This routine calculates and returns a 2D spline representing 
 * Ringleb's flow domain.                                       
 */
void Spline2D_HO::Create_Spline_Ringleb_Flow(void) {

  int nk, nq;
  double **rho;
  double delta_k, delta_q;
  double k_init, q_init, q_final;
  double **q, **k, qo, ko, c, J;
  double g = 1.40;
  Spline2D_HO Spline_North, Spline_East, Spline_West, Spline_Reflection;

  // Allocate memory.
  nk = 32;
  nq = 50;
  k = new double*[nk];
  for (int i = 0; i < nk; i++) k[i] = new double[nq];
  q = new double*[nk];
  for (int i = 0; i < nk; i++) q[i] = new double[nq];
  rho = new double*[nk];
  for (int i = 0; i < nk; i++) rho[i] = new double[nq];

  // Determine q, k.
  delta_k = (1.50-0.75)/double(nk-1);
  k_init = 0.75;
  q_init = 0.50;
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
  Spline_North.allocate(nk);  Spline_North.settype(SPLINE2D_QUINTIC);
  Spline_East.allocate(nq);   Spline_East.settype(SPLINE2D_QUINTIC);
  Spline_West.allocate(nq);   Spline_West.settype(SPLINE2D_QUINTIC);

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
	if (i == 0 || i == nk-1) Spline_North.tp[nk-1-i] = SPLINE2D_POINT_SHARP_CORNER;
	else Spline_North.tp[nk-1-i] = SPLINE2D_POINT_NORMAL;
      }
      // EAST spline.
      if (i == 0) {
	Spline_East.Xp[j].x = (HALF/rho[i][j])*(TWO/(k[i][j]*k[i][j]) - ONE/(q[i][j]*q[i][j])) - HALF*J;
	Spline_East.Xp[j].y = (ONE/(k[i][j]*rho[i][j]*q[i][j]))*sqrt(ONE - (q[i][j]*q[i][j])/(k[i][j]*k[i][j]));
	Spline_East.bc[j] = BC_FIXED;
	if (j == 0 || j == nq-1) Spline_East.tp[j] = SPLINE2D_POINT_SHARP_CORNER;
	else Spline_East.tp[j] = SPLINE2D_POINT_NORMAL;
      }
      // WEST spline.
      if (i == nk-1) {
	Spline_West.Xp[nq-1-j].x = (HALF/rho[i][j])*(TWO/(k[i][j]*k[i][j]) - ONE/(q[i][j]*q[i][j])) - HALF*J;
	Spline_West.Xp[nq-1-j].y = (ONE/(k[i][j]*rho[i][j]*q[i][j]))*sqrt(ONE - (q[i][j]*q[i][j])/(k[i][j]*k[i][j]));
	Spline_West.bc[nq-1-j] = BC_REFLECTION;
	if (j == 0 || j == nq-1) Spline_West.tp[nq-1-j] = SPLINE2D_POINT_SHARP_CORNER;
	else Spline_West.tp[nq-1-j] = SPLINE2D_POINT_NORMAL;
      }

    }
  }

  // Calculate spline path lengths.
  Spline_North.pathlength();
  Spline_East.pathlength();
  Spline_West.pathlength();

  // Set the boundary condition types for each of the interface splines.
  Spline_North.setBCtype(BC_RINGLEB_FLOW);
  Spline_East.setBCtype(BC_FIXED);
  Spline_West.setBCtype(BC_REFLECTION);

  // Reflect Spline_West and concatenate to Spline_West.
  Spline_Reflection = Spline_West;
  Spline_Reflection.Reflect_Spline();
  Spline_Reflection.Reverse_Spline();
  Spline_Reflection = Spline_Reflection + Spline_West;
  Spline_West = Spline_Reflection;

  // Reflect Spline_East and concatenate to Spline_East.
  Spline_Reflection = Spline_East;
  Spline_Reflection.Reflect_Spline();
  Spline_Reflection.Reverse_Spline();
  Spline_East = Spline_East + Spline_Reflection;

  // Concatenate sub-splines to the Ringleb spline.
  *this = Spline_West + Spline_North;
  *this = *this + Spline_East;

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


/*!
 * Compute the contour integral of the polynomial function 
 * \f$ p(x,y) = \oint (x-xc)^{(OrderX+1)} * (y-yc)^{(OrderY)} dy \f$   
 * between the path lengths corresponding to the "StartPoint" 
 * and "EndPoint" points with an accuracy given by the number 
 * of required exact "digits". 
 * Algorithm: The path integration (i.e. the spline) is divided into segments     
 * which contain only the end spline control points (i.e. no control points in between).
 * This division is neccessary for computations of curve integrals
 * along splines with sharp points.                                  
 * The integral is computed as a sum of the individual integrals on each individual segment.
 *
 * \param StartPoint the first point defining the integration domain
 * \param EndPoint the second point defining the integration domain
 * \param Centroid the point defining (xc,yc) in the expression above.
 * \param digits the number of exact digits with which the integration is sought
 * \param OrderX the power of the x-term in the polynomial expression
 * \param OrderY the power of the y-term in the polynomial expression
 */
double Spline2D_HO::PolynomOrderIntegration(const Vector2D &StartPoint, const Vector2D &EndPoint,
					    const Vector2D &Centroid, const int &digits,
					    const int &OrderX, const int &OrderY) const {

  double S_StartPoint, S_EndPoint;
  int i;
  vector<int> IndexCP;
  double Sol(0.0);

  // Get the pathlength for StartPoint and EndPoint
  S_StartPoint = getS(StartPoint);
  S_EndPoint   = getS(EndPoint); 

  // Determine whether there are spline control points between StartPoint and EndPoint
  // and whether they are far enough apart numerically in order to generate a sub-interval
  if (S_StartPoint < S_EndPoint){
    for(i=0; i<=np-1; ++i){
      if( ((sp[i] - S_StartPoint)*(sp[i] - S_EndPoint) < 0.0 ) &&
	  (fabs(S_StartPoint - sp[i])/(1.0 + fabs(sp[i])) >= 1.0e-14) &&
	  (fabs(  S_EndPoint - sp[i])/(1.0 + fabs(sp[i])) >= 1.0e-14) ){
	IndexCP.push_back(i);
      }
    }
  } else {
    for(i=np-1; i>=0; --i){
      if( (sp[i] - S_StartPoint)*(sp[i] - S_EndPoint) < 0.0 ){
	IndexCP.push_back(i);
      }
    }
  }

  // Divide the integral along the path into smaller integrals
  switch(IndexCP.size()){
  case 0: 			// no other control points in between
    return BasicOrderIntegration(StartPoint,EndPoint,Centroid,digits,OrderX,OrderY);

  case 1:			// one control point in between
    return ( BasicOrderIntegration(StartPoint,Xp[IndexCP[0]],Centroid,digits,OrderX,OrderY) + 
	     BasicOrderIntegration(Xp[IndexCP[0]],EndPoint,Centroid,digits,OrderX,OrderY) ) ;

  case 2:			// two control points in between
    return ( BasicOrderIntegration(StartPoint,Xp[IndexCP[0]],Centroid,digits,OrderX,OrderY) + 
	     BasicOrderIntegration(Xp[IndexCP[0]],Xp[IndexCP[1]],Centroid,digits,OrderX,OrderY)   +
	     BasicOrderIntegration(Xp[IndexCP[1]],EndPoint,Centroid,digits,OrderX,OrderY) ) ;

  default:			// more than two control points in between
    Sol = ( BasicOrderIntegration(StartPoint,Xp[IndexCP[0]],Centroid,digits,OrderX,OrderY) + 
	    BasicOrderIntegration(Xp[IndexCP[IndexCP.size()-1]],EndPoint,Centroid,digits,OrderX,OrderY) );

    for (i=0; i<(int)IndexCP.size()-1; ++i){
      Sol += BasicOrderIntegration(Xp[IndexCP[i]],Xp[IndexCP[i+1]],Centroid,digits,OrderX,OrderY);
    }
    return Sol;
  } // endswitch

}

/*!
 * Compute the contour integral of the polynomial function 
 * \f$ p(x,y) = \oint (x-xc)^{(OrderX+1)} * (y-yc)^{(OrderY)} dy \f$ 
 * between the path lengths corresponding to the "StartPoint" 
 * and "EndPoint" points with an accuracy given by the number 
 * of required exact "digits". 
 * It is assumed that there are no other spline control points 
 * between the end points.
 */
double Spline2D_HO::BasicOrderIntegration(const Vector2D &StartPoint, const Vector2D &EndPoint,
					  const Vector2D &Centroid, const int &digits,
					  const int &OrderX, const int &OrderY) const {

  // determine the accuracy based on the required precision (i.e the number of exact digits)
  double TOL(EpsilonTol::getAccuracyBasedOnExactDigits(digits));

  double S1(getS(StartPoint)), S2(getS(EndPoint));
  int Divisions(2),i;
  const int MaxDivisions(1000);

  double IntLessAccurate, IntMoreAccurate, RelError(10.0), Length(S2 - S1), DeltaS;
  Vector2D InterP1, InterP2;	// intermediate points

  if (fabs(Length) <= EpsilonTol::MachineEps){
    return 0.0;
  }

  // get a first estimation of the integral using StartPoint and EndPoint
  IntLessAccurate = PolynomLineIntegration2(StartPoint.x, StartPoint.y, EndPoint.x, EndPoint.y,
					    Centroid.x ,Centroid.y,
					    OrderX, OrderY);

  while( (RelError > TOL) && (Divisions<MaxDivisions) ){
    Divisions *= 2;		// double the resolution
    DeltaS = Length/Divisions;

    // get a more accurate estimation of the integral
    IntMoreAccurate = 0.0;
    InterP1 = StartPoint;	// initialize the start of the first interval
    for (i=1; i<=Divisions; ++i){
      InterP2 = Spline(S1 + i*DeltaS);
      IntMoreAccurate += PolynomLineIntegration2(InterP1.x,InterP1.y,InterP2.x,InterP2.y,
						 Centroid.x ,Centroid.y,
						 OrderX, OrderY);

      InterP1 = InterP2;	// set the new value for InterP1
    }

    // compute the relative error
    RelError = fabs(IntMoreAccurate - IntLessAccurate)/(1.0 + fabs(IntMoreAccurate));

    IntLessAccurate = IntMoreAccurate;
  }
  
  return IntMoreAccurate;
}

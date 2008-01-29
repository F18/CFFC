/*!\file HO_Spline2DInterval.cc
  \brief Source file initializing/implementing member variables/functions that belong to classes defined in HO_Spline2DInterval.h.*/

/* Include required C++ libraries. */
#include <iostream>

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "HO_Spline2DInterval.h"


/*!
 * Copy Constructor
 */
Spline2DInterval_HO::Spline2DInterval_HO(const Spline2DInterval_HO &rhs ): N_SubIntervals(0),
									   SubIntervalLength(NULL), N_GQP(0),
									   GQP(NULL), NormalGQP(NULL),
									   GQP_ContourInt(NULL), DyDs(NULL){
  int iter;			//< interval
  
  // allocate memory
  allocate(rhs.GQPointsPerSubInterval(), rhs.NumOfSubIntervals());
  
  // == Copy the data ==

  // Copy the length of the subintervals
  for(iter = 1; iter <= NumOfSubIntervals(); ++iter){
    SubIntervalLength[iter-1] = rhs.IntLength(iter);
  }

  // Copy the Gauss quadrature points for flux evaluation
  for(iter = 1; iter <= NumGQPoints(); ++iter){
    GQP[iter-1] = rhs.GQPoint(iter);
    NormalGQP[iter-1] = rhs.NormalGQPoint(iter);
  }
  
  // Copy the Gauss quadrature points for contour integration
  for(iter = 1; iter <= 3; ++iter){
    GQP_ContourInt[iter-1] = rhs.GQPointContourIntegral(iter);
    DyDs[iter - 1] = rhs.dYdS(iter);
  }
}

/*!
 * Assignment operator
 */
Spline2DInterval_HO& Spline2DInterval_HO::operator=(const Spline2DInterval_HO &rhs){

  int iter;

  if (this == &rhs) return *this;

  // allocate memory
  allocate(rhs.GQPointsPerSubInterval(), rhs.NumOfSubIntervals());

  // == Copy the data ==

  // Copy the length of the subintervals
  for(iter = 1; iter <= NumOfSubIntervals(); ++iter){
    SubIntervalLength[iter-1] = rhs.IntLength(iter);
  }

  // Copy the Gauss quadrature points for flux evaluation
  for(iter = 1; iter <= NumGQPoints(); ++iter){
    GQP[iter-1] = rhs.GQPoint(iter);
    NormalGQP[iter-1] = rhs.NormalGQPoint(iter);
  }
  
  // Copy the Gauss quadrature points for contour integration
  for(iter = 1; iter <= 3; ++iter){
    GQP_ContourInt[iter-1] = rhs.GQPointContourIntegral(iter);
    DyDs[iter - 1] = rhs.dYdS(iter);
  }

  return *this;
}


/*!
 * Allocate memory
 *
 * \param _N_GQP_ number of Gauss quadrature points per sub-interval (dictated by the accuracy of the scheme)
 * \param __N_SubIntervals__ number of sub-intervals (dictated by the number of sharp spline control
 *                                                    points between the end nodes)
 */
void Spline2DInterval_HO::allocate(const int _N_GQP_, const int __N_SubIntervals__){

  assert(__N_SubIntervals__ >= 1 && _N_GQP_ >= 1);
  int TotalEntries(_N_GQP_*__N_SubIntervals__); // total number of points along the spline interval
  int Total_MomentGQP_Entries(3 * __N_SubIntervals__);

  if ( (_N_GQP_ != GQPointsPerSubInterval()) || (__N_SubIntervals__ != NumOfSubIntervals()) ){
    // Deallocate the allocated memory
    deallocate();

    N_GQP = _N_GQP_;
    N_SubIntervals = __N_SubIntervals__;

    SubIntervalLength = new double [N_SubIntervals];
    // Set Lengths to zero
    for(int num_SubI=0; num_SubI<N_SubIntervals; ++num_SubI){ SubIntervalLength[num_SubI] = 0.0; }

    GQP = new Vector2D [TotalEntries];        //< Set the vector of GQPs for flux calculation 
    NormalGQP = new Vector2D [TotalEntries];  //< Set the vector of the normals at GQP locations

    GQP_ContourInt = new Vector2D [Total_MomentGQP_Entries]; //< Set the vector of the GQPs for geometric moments calculation
    DyDs = new double [Total_MomentGQP_Entries];	     //< Set the vector of derivatives DyDs at GQP_ContourInt locations

    // Set DyDs to zero
    for(int GQP=0; GQP<Total_MomentGQP_Entries; ++GQP){ DyDs[GQP] = 0.0; }
  }
}

/*!
 * Deallocate memory
 */
void Spline2DInterval_HO::deallocate(void){

  // Deallocate SubIntervalLength
  delete [] SubIntervalLength;
  SubIntervalLength = NULL;

  // Deallocate GQP
  delete [] GQP;
  GQP = NULL;

  // Deallocate NormalGQP
  delete [] NormalGQP;
  NormalGQP = NULL;

  // Deallocate GQP_ContourInt
  delete [] GQP_ContourInt;
  GQP_ContourInt = NULL;

  // Deallocate DyDs
  delete [] DyDs;
  DyDs = NULL;

  // Set to zero N_GQP and N_SubIntervals
  N_GQP = 0;			
  N_SubIntervals = 0;
}

/*!
 * Allocate memory and compute the interval information (i.e. GPQ Locations, Normal Values, Lengths etc.)
 *
 * \param SupportCurve the spline curve for which this spline interval is set.
 * \param StartPoint the first end of the interval
 * \param EndPoint the second end of the interval
 * \NumGQPoints the number of Gauss quadrature points per subinterval
 *
 * \todo Add code for the situation when there are sharp corners between StartPoint and EndPoint
 */
void Spline2DInterval_HO::InitializeInterval(const Spline2D_HO & SupportCurve, const Vector2D & StartPoint,
					     const Vector2D & EndPoint, const int &NumGQPoints){

  int num_sub_intervals(1);

  if (StartPoint != EndPoint){

    // Allocate memory 
    allocate(NumGQPoints, num_sub_intervals);

    // no sharp control points between StartPoint and EndPoint
    DetermineSubIntervalProperties(SupportCurve,StartPoint,EndPoint,1);	// Obs. "1" -- first interval
    
  } else {			// ERROR !!!
    throw runtime_error("Spline2DInterval_HO::InitializeInterval() ERROR! Zero length interval!\n");
  } // endif
}

/*!
 * Compute the properties of the specified subinterval.
 * These properties include the length of the subinterval,
 * the location of the Gauss quadrature points used for 
 * flux evaluation and the normals to the curve at these locations,  
 * the location of the Gauss quadrature points used for 
 * contour integration and the derivatives of y-coordinate with
 * respect to the pathlength at these locations.                                        
 * Method:                                             
 *  A discrete representation of the interval length as
 *  a function of pathlength is build using a fixed    
 *  number of points. Between two consecutive points a 
 *  linear variation is assumed.                       
 *  The location of the required number of Gauss       
 *  Quadrature points is obtained based on the         
 *  representation described above.
 */
void Spline2DInterval_HO::DetermineSubIntervalProperties(const Spline2D_HO & SupportCurve,
							 const Vector2D & StartPoint, const Vector2D & EndPoint,
							 const int &SubIntervalIndex){

  static const int MaxDivisions(50000);
  static const int DivisionFactor(4);
  int Divisions(1);
  double *Length(NULL), *PathLength(NULL);
  double* GQP_Abscissa = new double [N_GQP];
  double* GQP_ContourInt_Abscissa = new double [3];

  // Determine the pathlength coordinate of the StartPoint and the EndPoint
  double S1(SupportCurve.getS(StartPoint)),
         S2(SupportCurve.getS(EndPoint));

  double DeltaS(S2 - S1),
         DeltaS_Temp, RelativeError(1.0), GQP_PathLength(0);

  double dYdS_GQP, dXdS_GQP;	// derivative of x and y with respect to the pathlength (s) at the Gauss Quadrature Point

  int i, GQP_Iter;
  int subinterval_found(0);
  Vector2D InterP1, InterP2;

  // Initialize the interval length to the distance between the StartPoint and EndPoint
  SubIntervalLength[SubIntervalIndex-1] = abs(EndPoint - StartPoint);

  /* Based on the necessary number of divisions that are required 
     to obtained an imposed relative error for the length estimation,
     generate a discrete correspondence between the pathlength 
     (i.e. along the pathlength coordinate) and the spline interval
     length (i.e. along the true geometric spline contour). */
  while( (RelativeError > EpsilonTol::epsilon_relative) && (Divisions <= MaxDivisions) ){

    // Update number of divisions
    Divisions *= DivisionFactor;
    // Update DeltaS_Temp
    DeltaS_Temp = DeltaS/Divisions;
    // Check if DeltaS_Temp is different than machine epsilon
    if ( DeltaS + DeltaS_Temp == DeltaS){
      break;
    }

    // Deallocate previously allocated memory if necessary
    delete [] Length; Length = NULL;
    delete [] PathLength; PathLength = NULL;

    // Allocate memory for the discrete correspondence
    Length = new double [Divisions+1];
    PathLength = new double [Divisions+1];

    // Initialize the first correspondence
    PathLength[0] = S1;
    Length[0] = 0.0;

    // Reset temporary variables
    InterP1 = StartPoint;

    // Compute with the current DeltaS_Temp
    for (i=1; i<=Divisions; ++i){
      PathLength[i] = S1 + i*DeltaS_Temp;                 // get the pathlength
      InterP2 = SupportCurve.Spline(PathLength[i]);       // get the corresponding point on the SupportCurve
      Length[i] = Length[i-1] + abs(InterP2 - InterP1);   // get the distance to the previous point and
                                                          // add it to the previous length
      InterP1 = InterP2;	                          // set the new value for InterP1
    }
    
    // Compute Relative Error
    RelativeError = ( fabs(Length[Divisions] - SubIntervalLength[SubIntervalIndex-1])/
		      (1.0 + fabs(SubIntervalLength[SubIntervalIndex-1])) );

    // Update SubIntervalLength with the more accurate calculation
    SubIntervalLength[SubIntervalIndex-1] = Length[Divisions];
  }

  // Get the lengths to the Gauss quadrature points (GQP) for flux evaluation
  switch(N_GQP){
  case 1:			// One Gauss Quad Point
    GQP_Abscissa[0] = GaussQuadratureData::GQ1_Abscissa[0]*IntLength(SubIntervalIndex);
    break;

  case 2:			// Two Gauss Quad Points
    GQP_Abscissa[0] = GaussQuadratureData::GQ2_Abscissa[0]*IntLength(SubIntervalIndex);
    GQP_Abscissa[1] = GaussQuadratureData::GQ2_Abscissa[1]*IntLength(SubIntervalIndex);
    break;

  case 3:			// Three Gauss Quad Points
    GQP_Abscissa[0] = GaussQuadratureData::GQ3_Abscissa[0]*IntLength(SubIntervalIndex);
    GQP_Abscissa[1] = GaussQuadratureData::GQ3_Abscissa[1]*IntLength(SubIntervalIndex);
    GQP_Abscissa[2] = GaussQuadratureData::GQ3_Abscissa[2]*IntLength(SubIntervalIndex);
    break;

  default:
    throw runtime_error("Spline2DInterval_HO::DetermineSubIntervalProperties() ERROR!\n\
                         The current implementation works only up to 3 Gauss quadrature points.\n");
  }//endswitch


  // Identify the corresponding position for each GQP
  for (GQP_Iter=0; GQP_Iter<N_GQP; ++GQP_Iter){
    i=0; subinterval_found = 0;
    while ( i<Divisions && subinterval_found==0 ){
      // check the position
      if ( (GQP_Abscissa[GQP_Iter] - Length[i])*(GQP_Abscissa[GQP_Iter] - Length[i+1]) <= ZERO ) {
	// the subinterval has been found
	subinterval_found = 1;
	// use local linear approximation to determine the pathlength for the current GQP
	GQP_PathLength = ( PathLength[i] + 
			   (GQP_Abscissa[GQP_Iter]-Length[i])*(PathLength[i+1]-PathLength[i])/(Length[i+1]-Length[i]) );
	// determine the (x,y) coordinates for the current GQP_PathLength
	InterP1 = SupportCurve.Spline(GQP_PathLength);
	// set the location of the current GQP
	GQP[GQP_Iter + (SubIntervalIndex-1)*N_GQP] = InterP1;

 	// determine the unit normal vector to the spline at the GQP location
	NormalGQP[GQP_Iter + (SubIntervalIndex-1)*N_GQP] = SupportCurve.nSpline(InterP1, dXdS_GQP, dYdS_GQP);
      }	// endif

      // Advance to the next subinterval
      ++i;
    }//endwhile
  }//endfor


  // Get the lengths to the Gauss quadrature points (GQP_ContourInt) for contour integration
  GQP_ContourInt_Abscissa[0] = GaussQuadratureData::GQ3_Abscissa[0]*IntLength(SubIntervalIndex);
  GQP_ContourInt_Abscissa[1] = GaussQuadratureData::GQ3_Abscissa[1]*IntLength(SubIntervalIndex);
  GQP_ContourInt_Abscissa[2] = GaussQuadratureData::GQ3_Abscissa[2]*IntLength(SubIntervalIndex);

  // Identify the corresponding position for each GQP_ContourInt
  for (GQP_Iter=0; GQP_Iter<3; ++GQP_Iter){
    i=0; subinterval_found = 0;
    while ( i<Divisions && subinterval_found==0 ){
      // check the position
      if ( (GQP_ContourInt_Abscissa[GQP_Iter] - Length[i])*(GQP_ContourInt_Abscissa[GQP_Iter] - Length[i+1]) <= ZERO ) {
	// the subinterval has been found
	subinterval_found = 1;
	// use local linear approximation to determine the pathlength for the current GQP
	GQP_PathLength = ( PathLength[i] + 
			   (GQP_ContourInt_Abscissa[GQP_Iter]-Length[i])*(PathLength[i+1]-PathLength[i])/(Length[i+1]-Length[i]) );
	// determine the (x,y) coordinates for the current GQP_PathLength
	InterP1 = SupportCurve.Spline(GQP_PathLength);
	// set the location of the current GQP_ContourInt
	GQP_ContourInt[GQP_Iter + (SubIntervalIndex-1)*3] = InterP1;

 	// determine DyDs at the current GQP_ContourInt
	SupportCurve.nSpline(InterP1, dXdS_GQP, dYdS_GQP);

	// store DyDs at the current GQP
	DyDs[GQP_Iter + (SubIntervalIndex-1)*3] = dYdS_GQP;
      }	// endif

      // Advance to the next subinterval
      ++i;
    }//endwhile
  }//endfor

  delete [] Length; Length = NULL;
  delete [] PathLength; PathLength = NULL;
  delete [] GQP_Abscissa; GQP_Abscissa = NULL;
  delete [] GQP_ContourInt_Abscissa; GQP_ContourInt_Abscissa = NULL;
}

/*!
 * Update all interval properties with the input parameters 
 */
void Spline2DInterval_HO::UpdateInterval(const Spline2D_HO & SupportCurve,
					 const Vector2D & StartPoint, const Vector2D & EndPoint,
					 const int &NumGQPoints){

  // Determine the number of sub-intervals
  int num_sub_intervals(1);

  if (StartPoint != EndPoint){

    // Allocate memory if required
    allocate(NumGQPoints, num_sub_intervals);
    
    // Determine the interval properties when there are no sharp control points between StartPoint and EndPoint
    DetermineSubIntervalProperties(SupportCurve,StartPoint,EndPoint,1);	// Obs. "1" -- first interval

  } else {			// ERROR !!!
    throw runtime_error("Spline2DInterval_HO::UpdateInterval() ERROR! Zero length interval!");
  } // endif
}


/*!
 * Compute the contribution of an individual Gauss quadrature point to 
 * the contour integral of a function with the form
 *  \f$ p(x,y) = (x - x_i)^{n+1} (y - y_i)^m \frac{\partial y(s)}{ds} \f$
 * \param OrderX the n-power in the expression above
 * \param OrderY the m-power in the expression above
 * \param Centroid the point providing (xi,yi)-Cartesian coordinates 
 */
double Spline2DInterval_HO::BasicTermIntegration(const Vector2D & Centroid,
						 const int &OrderX, const int &OrderY,
						 const int &GQP){

  double diffX(GQPointContourIntegral(GQP).x - Centroid.x);
  double diffY(GQPointContourIntegral(GQP).y - Centroid.y);
  int i;

  double diffX_Power(diffX), diffY_Power(1.0); //< Set the values for OrderX=0 and OrderY=0

  //  Compute diffX_Power
  for (i=1; i<=OrderX; ++i){
    diffX_Power *= diffX;
  }
  
  // Compute diffY_Power
  for (i=1; i<=OrderY; ++i){
    diffY_Power *= diffY;
  }
  
  return diffX_Power * diffY_Power * dYdS(GQP);
}

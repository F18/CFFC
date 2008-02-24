/*!\file HO_Spline2DInterval.h 
  \brief Header file defining an interval class for the high-order 2D Spline. */

#ifndef _HO_SPLINE2DINTERVAL_INCLUDED
#define _HO_SPLINE2DINTERVAL_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "HO_Spline2D.h"


/* Define the 2D spline interval class. */

/*!
 * \class Spline2DInteval_HO
 *
 * \brief Storage class for the geometric properties of a
 * spline segment, which is the spline region between 
 * two grid nodes.
 *
 * \verbatim
 * A spline segment can include multiple smooth or sharp 
 * spline control points.
 *                                          
 * Member Variables:                                    
 *  N_SubIntervals -- Number of spline regions which    
 *                    don't have sharp control points as
 *                    interior points.                  
 *  SubIntervalLength -- Array for storing the length of 
 *                       each spline subinterval.               
 *    N_GQP        -- The number of Gauss quadrature points (GQP) 
 *                    per subinterval at which the flux is evaluated.
 *     GQP         -- Storage for the location of the   
 *                    Gauss Quadrature points at which  
 *                    flux calculations are performed.
 *  NormalGQP      -- Storage for the normal vectors at 
 *                    the locations of GQP.
 *  GQP_Moment     -- Array for storing all the positions of GQPs 
 *                    used for computing the contour integrals
 *                    along the spline interval.
 *    DyDs         -- Array for storing the derivatives of y-coordinate
 *                    with respect the curvilinear coordinate s.
 * \endverbatim                                                      
 ************************************************************************/
class Spline2DInterval_HO{

public:

  //! @name Constructors, destructor and assignment operator
  //@{
  //! Default constructor.
  Spline2DInterval_HO(void);

  //! Constructor with number of Gauss points per subinterval and number of subintervals.
  Spline2DInterval_HO(const int &_N_GQP_,
		      const int &__N_SubIntervals__ = 1){ allocate(_N_GQP_,__N_SubIntervals__); }

  //! Copy constructor
  Spline2DInterval_HO(const Spline2DInterval_HO &rhs );

  //! Assignment operator
  Spline2DInterval_HO & operator=(const Spline2DInterval_HO &rhs);

  //! Destructor.
  ~Spline2DInterval_HO(void){ deallocate(); }
  //@}

  //! @name Allocation and deallocation.
  //@{
  //! Allocate memory for the spline interval.
  void allocate(const int _N_GQP_, const int __N_SubIntervals__ = 1);
  //! Deallocate memory for the spline interval.
  void deallocate(void);
  //@}

  //! @name Private variables access functions
  //@{
  //! Get the number of subintervals
  const int& NumOfSubIntervals(void) const {return N_SubIntervals; }
  //! Get the number of Gauss quadrature points per subinterval used for flux evaluation
  const int& GQPointsPerSubInterval(void) const {return N_GQP; }
  //! Get the total number of Gauss quadrate points used for flux evaluation
  int NumGQPoints(void) const {return N_GQP*N_SubIntervals; }
  //! Get the total number of Gauss quadrate points used for contour integration
  int NumGQPoints_ContourIntegral(void) const {return NUMBER_OF_GQP_CONTOURINT * N_SubIntervals; }

  //! Get the array of subinterval length
  const double* IntLength(void) const {return SubIntervalLength; }
  /*! 
   * Get a specific subinterval length
   * \param Position index specifying which subinterval length is retrieved. Position >= ONE!!!
   */
  const double& IntLength(const int &Position) const {return SubIntervalLength[Position-1];}
  double& IntLength(const int &Position) {return SubIntervalLength[Position-1];}

  //! Get the array of Gauss quadrature points
  const Vector2D* GQPoint(void) const {return GQP;}
  /*! 
   * Get a specific Gauss quadrature point location. 
   * \param Position index specifying which GQP is retrieved. Position >= ONE!!!
   */
  const Vector2D& GQPoint(const int &Position) const {return GQP[Position-1];}
  Vector2D& GQPoint(const int &Position) {return GQP[Position-1];}

  //! Get the array of normals at the Gauss quadrature points
  const Vector2D* NormalGQPoint(void) const {return NormalGQP;}
  /*! 
   * Get the Gauss quadrature point location. 
   * \param Position index specifying which normal is retrieved. Position >= ONE!!!
   */
  const Vector2D& NormalGQPoint(const int &Position) const {return NormalGQP[Position-1];}
  Vector2D& NormalGQPoint(const int &Position) {return NormalGQP[Position-1];}

  //! Get the array of Gauss quadrature points used for contour integration.
  const Vector2D* GQPointContourIntegral(void) const {return GQP_ContourInt;}
  /*! 
   * Get a specific Gauss quadrature point location used for contour integration. 
   * \param Position index specifying which location of the array is retrieved. Position >= ONE!!!
   */
  const Vector2D& GQPointContourIntegral(const int &Position) const {return GQP_ContourInt[Position-1];}
  Vector2D& GQPointContourIntegral(const int &Position) {return GQP_ContourInt[Position-1];}

  //! Get the array of y-coordinate derivatives with respect which curvilinear coordinate.
  const double* dYdS(void) const {return DyDs; }

  /*! 
   * Get a specific dYdS derivative used for contour integration. 
   * \param Position index specifying which location of the array is retrieved. Position >= ONE!!!
   */
  const double& dYdS(const int &Position) const {return DyDs[Position-1];}
  double& dYdS(const int &Position) {return DyDs[Position-1];}
  //@}

  //! @name Interval initialization and update functions
  //@{
  //! Initialize the spline interval
  void InitializeInterval(const Spline2D_HO & SupportCurve,
			  const Vector2D & StartPoint, const Vector2D & EndPoint,
			  const int &NumGQPoints);
  //! Initialize the spline interval when the end points are Node2D_HO
  void InitializeInterval(const Spline2D_HO & SupportCurve,
			  const Node2D_HO & StartPoint, const Node2D_HO & EndPoint,
			  const int &NumGQPoints){
    return InitializeInterval(SupportCurve,StartPoint.X, EndPoint.X,NumGQPoints);
  }

  //! Update properties of the sline interval
  void UpdateInterval(const Spline2D_HO & SupportCurve,
		      const Vector2D & StartPoint, const Vector2D & EndPoint,
		      const int &NumGQPoints);
  //! Update properties of the sline interval when the end points are Node2D_HO
  void UpdateInterval(const Spline2D_HO & SupportCurve,
		      const Node2D_HO & StartPoint, const Node2D_HO & EndPoint,
		      const int &NumGQPoints){
    return UpdateInterval(SupportCurve,StartPoint.X, EndPoint.X, NumGQPoints);
  }

  /*!
   * Update all interval properties with the input parameters for a number of GQP already setup.
   */
  void UpdateInterval(const Spline2D_HO & SupportCurve,
		      const Vector2D & StartPoint, const Vector2D & EndPoint){
    UpdateInterval(SupportCurve, StartPoint, EndPoint, GQPointsPerSubInterval());
  }
  //! Update spline interval when the end points are Node2D_HO
  void UpdateInterval(const Spline2D_HO & SupportCurve,
		      const Node2D_HO & StartPoint, const Node2D_HO & EndPoint){
    return UpdateInterval(SupportCurve,StartPoint.X, EndPoint.X);
  }
  //@}

  //! @name Contour integration functions
  //@{
  double IntegratePolynomialTerm(const Vector2D & Centroid, const int &OrderX, const int &OrderY);
  double IntegratePolynomialTermAndDivide(const Vector2D & Centroid, const int &OrderX, const int &OrderY);
  double BasicTermIntegration(const Vector2D & Centroid, const int &OrderX, const int &OrderY, const int &GQP);
  double AreaContribution(void){ return IntegratePolynomialTerm(Vector2D(0.0), 0, 0); }
  double XCentroidContribution(void){ return IntegratePolynomialTerm(Vector2D(0.0), 1, 0); }
  double YCentroidContribution(void){ return IntegratePolynomialTerm(Vector2D(0.0), 0, 1); }
  //@}

  /* Operators */
  friend std::ostream &operator << (std::ostream &os, const Spline2DInterval_HO &S);
  friend std::istream &operator >> (std::istream &is, Spline2DInterval_HO &S);

private:
  static int NUMBER_OF_GQP_CONTOURINT; //!< Number of Gauss Quadrature Points per subinterval for curvilinear integration
  int N_SubIntervals; 		//!< Number of regions for which sharp control points can exist only at extremities.
  double* SubIntervalLength;	//!< Array for storing the lengths of all subintervals.
  int N_GQP;                    //!< Number of Gauss quadrature points per subinterval (depends on the order of accuracy).

  Vector2D* GQP;		//!< Array for storing the Gauss Quadrature Points for all subintervals.
  Vector2D* NormalGQP;	        //!< Array for storing the normal vectors at Gauss Quadrature locations for all subintervals.

  Vector2D* GQP_ContourInt;	/*!< Array for storing all Gauss Quadrature Points (i.e. for all subintervals)
				  that are used to compute contour integrals along the supporting spline.
				  N.B. The number of GQP is not always
				  sufficient for obtaining the exact solution of the curvilinear integrals of polynomial 
				  functions that currently arise in the calculation of geometric moments.
				  Therefore, this variable stores NUMBER_OF_GQP_CONTOURINT Gauss integration points per subinterval.
				  This number of points is enough to integrate exactly polynomials up to fifth-order. */
  double* DyDs;			//!< Array for storing the values of the derivatives dyds at each Gauss Quadrature Point

  //! Determine the subinterval properties
  void DetermineSubIntervalProperties(const Spline2D_HO & SupportCurve,
				      const Vector2D & StartPoint, const Vector2D & EndPoint,
				      const int &SubIntervalIndex);
};

/*!
 * Default Constructor
 */
inline Spline2DInterval_HO::Spline2DInterval_HO(void): N_SubIntervals(0), N_GQP(0),
						       SubIntervalLength(NULL), GQP(NULL),
						       NormalGQP(NULL),
						       GQP_ContourInt(NULL), DyDs(NULL) { 
  // 
}

/*!
 * Output operator.
 */
inline std::ostream &operator << (std::ostream &os, const Spline2DInterval_HO &S) {
  int i;

  // Output general information
  os << " " << S.NumOfSubIntervals() << "\t" << S.GQPointsPerSubInterval() << "\n";

  // Output the length of each subinterval
  for (i=1; i<=S.NumOfSubIntervals(); ++i){
    os << S.IntLength(i) << "\n";
  }

  // Output the Gauss quadrature points and the normals used for flux evaluation
  for (i=1; i<=S.NumGQPoints(); ++i){
    os << setprecision(16) 
       << S.GQPoint(i) << "\n" 
       << S.NormalGQPoint(i) << "\n";
  }

  // Output the Gauss quadrate points and the y-coordinate derivatives used for contour integration
  for (i=1; i<=S.NumGQPoints_ContourIntegral(); ++i){
    os << setprecision(16) 
       << S.GQPointContourIntegral(i) << "\n"
       << S.dYdS(i) << "\n";
  }

  return os;
}

/*!
 * Input operator.
 */
inline std::istream &operator >> (std::istream &is, Spline2DInterval_HO &S) {
  
  int NumSubIntervals, NumGaussQuadPoints;
  int i;

  // Read general information (i.e. number of subintervals and number of Gauss quadratures per subinterval)
  is.setf(ios::skipws);
  is >> NumSubIntervals >> NumGaussQuadPoints;
  is.unsetf(ios::skipws);

  // allocate memory
  S.allocate(NumGaussQuadPoints, NumSubIntervals);

  // Read the length of each subinterval
  is.setf(ios::skipws);
  for (i=1; i<=S.NumOfSubIntervals(); ++i){
    is >> S.IntLength(i);
  }
  is.unsetf(ios::skipws);
  
  // Read the Gauss quadrature points an		  test_HO_Spline2D.cc \d the normals used for flux evaluation
  for (i=1; i<=S.NumGQPoints(); ++i){
    is >> S.GQPoint(i) ;
    is >> S.NormalGQPoint(i);
  }

  // Read the Gauss quadrature points and the y-coordinate derivatives used for contour integration
  for (i=1; i<=S.NumGQPoints_ContourIntegral(); ++i){
    is >> S.GQPointContourIntegral(i) ;
    is.setf(ios::skipws);
    is >> S.dYdS(i);
    is.unsetf(ios::skipws);
  }
  
  return is;
}

/*!
 * Compute the expression \f$ I = \oint (x - x_i)^{n+1} (y - y_i)^m dy \f$ .
 * This integral is part of the typical expression 
 * \f$ I_{C} = \frac{1}{n + 1} \oint (x - x_i)^{n+1} (y - y_i)^m dy \f$ 
 * which arises in the calculation of geometric moments for curved boundaries.
 * To evaluate this expression numerically, a Gauss quadrature integration is performed along 
 * the spline interval, based on the following values: GQP_ContourInt, SubIntervalLength and DyDs.
 * The expression is rewritten only as a function of the curvilinear coordinate:
 * \f$ I = \oint (x - x_i)^{n+1} (y - y_i)^m \frac{\partial y(s)}{ds} ds \f$
 * Multiplication with the factor \f$ \frac{1}{n + 1} \f$ is done after all the contour 
 * component integrals have been summed up.
 */
inline double Spline2DInterval_HO::IntegratePolynomialTerm(const Vector2D & Centroid,
							   const int &OrderX, const int &OrderY){

  double Result(0.0);

  switch(NUMBER_OF_GQP_CONTOURINT){

  case 3:
    for (int i=0; i<NumOfSubIntervals(); ++i){
      // Calculate the contribution of each subinterval
      Result +=  IntLength(i+1)*( GaussQuadratureData::GQ3_Weight[0] * BasicTermIntegration(Centroid,OrderX,OrderY,3*i + 1) + 
				  GaussQuadratureData::GQ3_Weight[1] * BasicTermIntegration(Centroid,OrderX,OrderY,3*i + 2) +
				  GaussQuadratureData::GQ3_Weight[2] * BasicTermIntegration(Centroid,OrderX,OrderY,3*i + 3) );
    }
    break;

  case 5:
    for (int i=0; i<NumOfSubIntervals(); ++i){
      // Calculate the contribution of each subinterval
      Result +=  IntLength(i+1)*( GaussQuadratureData::GQ5_Weight[0] * BasicTermIntegration(Centroid,OrderX,OrderY,5*i + 1) + 
				  GaussQuadratureData::GQ5_Weight[1] * BasicTermIntegration(Centroid,OrderX,OrderY,5*i + 2) +
				  GaussQuadratureData::GQ5_Weight[2] * BasicTermIntegration(Centroid,OrderX,OrderY,5*i + 3) + 
				  GaussQuadratureData::GQ5_Weight[3] * BasicTermIntegration(Centroid,OrderX,OrderY,5*i + 4) +
				  GaussQuadratureData::GQ5_Weight[4] * BasicTermIntegration(Centroid,OrderX,OrderY,5*i + 5) );
    }
    break;
  }

  return Result;
}

/*!
 * Compute the expression \f$ I = \frac{1}{n + 1} \oint (x - x_i)^{n+1} (y - y_i)^m dy \f$ .
 * This integral arises typically in the calculation of geometric moments for curved boundaries.
 * This subroutine divide the integral I with the factor (OrderX+1).
 */
inline double Spline2DInterval_HO::IntegratePolynomialTermAndDivide(const Vector2D & Centroid,
								    const int &OrderX, const int &OrderY){

  return IntegratePolynomialTerm(Centroid,OrderX,OrderY)/(OrderX+1);
}

/*!
 * Update all interval properties with the input parameters 
 */
inline void Spline2DInterval_HO::UpdateInterval(const Spline2D_HO & SupportCurve,
						const Vector2D & StartPoint, const Vector2D & EndPoint,
						const int &NumGQPoints){
  InitializeInterval(SupportCurve,StartPoint,EndPoint,NumGQPoints);
}


#endif /* _HO_SPLINE2DINTERVAL_INCLUDED  */

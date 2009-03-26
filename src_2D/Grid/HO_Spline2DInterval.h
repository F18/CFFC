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
  //! Get the total number of Gauss quadrature points used for contour integration
  int NumGQPoints_ContourIntegral(void) const {return NUMBER_OF_GQP_CONTOURINT * N_SubIntervals; }
  //! Get the number of Gauss quadrature points used for contour integration along each subinterval
  static const int& get_NumGQPoints_ContourIntegral(void) {return NUMBER_OF_GQP_CONTOURINT; }

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

  //! Copy the array of Gauss quadrature points used for flux calculation in the provided container
  void CopyGQPoints(Vector2D *CopyGQP){ for (int i=0; i<NumGQPoints(); ++i) CopyGQP[i] = GQP[i]; }
  //! Copy the array of Gauss quadrature points used for flux calculation in the provided array container
  void CopyGQPoints(std::vector<Vector2D> &CopyGQP){ for (int i=0; i<NumGQPoints(); ++i) CopyGQP.push_back(GQP[i]); }

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

  //! Copy the array of normals at the Gauss quadrature points used for flux calculation in the provided container
  void CopyNormalGQPoints(Vector2D *CopyNormal){ for (int i=0; i<NumGQPoints(); ++i) CopyNormal[i] = NormalGQP[i]; }
  //! Copy the array of normals at the Gauss quadrature points used for flux calculation in the provided array container
  void CopyNormalGQPoints(std::vector<Vector2D> &CopyNormal){ for (int i=0; i<NumGQPoints(); ++i) CopyNormal.push_back(NormalGQP[i]);}

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
  double AreaContribution(const Vector2D &RefPoint){ return IntegratePolynomialTerm(RefPoint, 0, 0); }
  double XCentroidContribution(const Vector2D &RefPoint){ return IntegratePolynomialTerm(RefPoint, 1, 0); }
  double YCentroidContribution(const Vector2D &RefPoint){ return IntegratePolynomialTerm(RefPoint, 0, 1); }

  //! Integrate the desired function with respect to the y-coordinate along the contour of this spline interval.
  template<typename FO, class ReturnType>
  ReturnType IntegrateFunctionWithRespectToY(FO FuncObj, ReturnType _dummy_param) const;

  //! @brief Integrate the given function multiplied with the normal along the contour of this spline interval.
  template<typename FO, class ReturnType>
  void IntegrateFunctionProjectionOnNormalDirections(const Spline2D_HO & SupportCurve, FO FuncObj,
						     ReturnType & ResultXdir, ReturnType & ResultYdir,
						     double & WettedSurface) const;

  //! @brief Integrate the given function along the contour of this spline interval.
  template<typename FO, class ReturnType>
  ReturnType IntegrateFunctionAlongInterval(FO FuncObj, ReturnType _dummy_param) const;

  //! @brief Integrate the wall shear stress along the contour of this spline interval.
  template<typename FO>
  void IntegrateWallShearStressContributions(const Spline2D_HO & SupportCurve, 
					     FO FuncObj,
					     double & ResultXdir, double & ResultYdir, double & WettedSurface) const;
  //@}

  //! @name Set the number of Gauss quadrature points for contour integration
  //@{
  static void Set_Default_Parameters(void){ NUMBER_OF_GQP_CONTOURINT = 3; }
  static void setThreePointGaussQuadContourIntegration(void){ NUMBER_OF_GQP_CONTOURINT = 3; }
  static void setFivePointGaussQuadContourIntegration(void){ NUMBER_OF_GQP_CONTOURINT = 5; }
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

/*!
 * Compute the expression \f$ I = \oint F(x,y) dy \f$.
 * The expression is rewritten only as a function of the curvilinear coordinate:
 * \f$ I = \oint F(x,y) \frac{\partial y(s)}{ds} ds \f$
 * The function to be integrated along the contour of 
 * this spline interval must be differentiable, such that
 * Gauss quadrature integration can be applied.
 * The accuracy of this integration is determined by the 
 * number of Gauss quadrature points used in the integration
 * process (NUMBER_OF_GQP_CONTOURINT).
 * Some operators and constructors are required to be defined
 * for the ReturnType class, and the function is assumed to
 * be function only of the 2D Cartesian coordinates.
 */
template<typename FO, class ReturnType> inline
ReturnType Spline2DInterval_HO::IntegrateFunctionWithRespectToY(FO FuncObj, ReturnType _dummy_param) const {

  ReturnType Result(0.0);

  switch(NUMBER_OF_GQP_CONTOURINT){

  case 3:
    for (int i=0; i<NumOfSubIntervals(); ++i){
      // Calculate the contribution of each subinterval
      Result +=  IntLength(i+1)*( GaussQuadratureData::GQ3_Weight[0] * dYdS(3*i+1) * FuncObj(GQPointContourIntegral(3*i+1).x,
											     GQPointContourIntegral(3*i+1).y) + 
				  GaussQuadratureData::GQ3_Weight[1] * dYdS(3*i+2) * FuncObj(GQPointContourIntegral(3*i+2).x,
											     GQPointContourIntegral(3*i+2).y) +
				  GaussQuadratureData::GQ3_Weight[2] * dYdS(3*i+3) * FuncObj(GQPointContourIntegral(3*i+3).x,
											     GQPointContourIntegral(3*i+3).y) );
    }
    break;

  case 5:
    for (int i=0; i<NumOfSubIntervals(); ++i){
      // Calculate the contribution of each subinterval
      Result +=  IntLength(i+1)*( GaussQuadratureData::GQ5_Weight[0] * dYdS(5*i+1) * FuncObj(GQPointContourIntegral(5*i+1).x,
											     GQPointContourIntegral(5*i+1).y) + 
				  GaussQuadratureData::GQ5_Weight[1] * dYdS(5*i+2) * FuncObj(GQPointContourIntegral(5*i+2).x,
											     GQPointContourIntegral(5*i+2).y) +
				  GaussQuadratureData::GQ5_Weight[2] * dYdS(5*i+3) * FuncObj(GQPointContourIntegral(5*i+3).x,
											     GQPointContourIntegral(5*i+3).y) +
				  GaussQuadratureData::GQ5_Weight[3] * dYdS(5*i+4) * FuncObj(GQPointContourIntegral(5*i+4).x,
											     GQPointContourIntegral(5*i+4).y) +
				  GaussQuadratureData::GQ5_Weight[4] * dYdS(5*i+5) * FuncObj(GQPointContourIntegral(5*i+5).x,
											     GQPointContourIntegral(5*i+5).y) );
    }
    break;
  }

  return Result;
}

/*! 
 * Integrate the given function multiplied with the normal along 
 * the contour of this spline interval.
 * One usage of this routine is to compute the integral of pressure
 * forces that act on the surface of a solid body.
 * 
 * \param SupportCurve the spline associated with this interval. 
 *                     It's required to compute the normals.
 * \param FuncObj  the function to be integrated
 * \param ResultXdir the result in the X-direction (i.e FuncObj * n_x)
 * \param ResultYdir the result in the Y-direction (i.e FuncObj * n_y)
 *
 * \note The contributions of this interval are ADDED to the result variables!!!
 */
template<typename FO, class ReturnType>
void Spline2DInterval_HO::IntegrateFunctionProjectionOnNormalDirections(const Spline2D_HO & SupportCurve, FO FuncObj,
									ReturnType & ResultXdir,
									ReturnType & ResultYdir,
									double & WettedSurface) const {
  
  ReturnType FuncVal(0.0), TempX(0.0), TempY(0.0);
  double const * GQ_Weight(NULL);
  int GQP, index;
  Vector2D Normal;

  // Set the weights correctly
  switch(NUMBER_OF_GQP_CONTOURINT){
  case 3:
    GQ_Weight = &GaussQuadratureData::GQ3_Weight[0];
    break;
    
  case 5:
    GQ_Weight = &GaussQuadratureData::GQ5_Weight[0];
    break;
  }
  
  for (int i=0; i<NumOfSubIntervals(); ++i){
    // Calculate the contribution of each subinterval
    
    for (GQP = 1; GQP <= NUMBER_OF_GQP_CONTOURINT; ++GQP){
      // Calculate the contribution of each Gauss-quadrature point
      
      // Determine the index of the current GQP in the array
      index = NUMBER_OF_GQP_CONTOURINT*i+GQP;

      // Evaluate weighted function at the current GQP
      FuncVal = GQ_Weight[GQP-1] * FuncObj(GQPointContourIntegral(index).x,
					   GQPointContourIntegral(index).y);

      // Normal vector
      Normal = SupportCurve.nSpline(GQPointContourIntegral(index));

      // Update the X projection
      TempX += FuncVal* Normal.x;

      // Update the Y projection
      TempY += FuncVal* Normal.y;
    }
  
    // Update final result with the contribution of the current subinterval
    ResultXdir += TempX * IntLength(i+1);
    ResultYdir += TempY * IntLength(i+1);
    WettedSurface += IntLength(i+1);
  }

  GQ_Weight = NULL;
}

/*!
 * Integrate the given function along the interval and return the result.
 */
template<typename FO, class ReturnType>
ReturnType Spline2DInterval_HO::IntegrateFunctionAlongInterval(FO FuncObj, ReturnType _dummy_param) const {

  ReturnType Result(0.0);

  switch(NUMBER_OF_GQP_CONTOURINT){

  case 3:
    for (int i=0; i<NumOfSubIntervals(); ++i){
      // Calculate the contribution of each subinterval
      Result +=  IntLength(i+1)*( GaussQuadratureData::GQ3_Weight[0] * FuncObj(GQPointContourIntegral(3*i+1).x,
									       GQPointContourIntegral(3*i+1).y) + 
				  GaussQuadratureData::GQ3_Weight[1] * FuncObj(GQPointContourIntegral(3*i+2).x,
									       GQPointContourIntegral(3*i+2).y) +
				  GaussQuadratureData::GQ3_Weight[2] * FuncObj(GQPointContourIntegral(3*i+3).x,
									       GQPointContourIntegral(3*i+3).y) );
    }
    break;

  case 5:
    for (int i=0; i<NumOfSubIntervals(); ++i){
      // Calculate the contribution of each subinterval
      Result +=  IntLength(i+1)*( GaussQuadratureData::GQ5_Weight[0] * FuncObj(GQPointContourIntegral(5*i+1).x,
									       GQPointContourIntegral(5*i+1).y) + 
				  GaussQuadratureData::GQ5_Weight[1] * FuncObj(GQPointContourIntegral(5*i+2).x,
									       GQPointContourIntegral(5*i+2).y) +
				  GaussQuadratureData::GQ5_Weight[2] * FuncObj(GQPointContourIntegral(5*i+3).x,
									       GQPointContourIntegral(5*i+3).y) +
				  GaussQuadratureData::GQ5_Weight[3] * FuncObj(GQPointContourIntegral(5*i+4).x,
									       GQPointContourIntegral(5*i+4).y) +
				  GaussQuadratureData::GQ5_Weight[4] * FuncObj(GQPointContourIntegral(5*i+5).x,
									       GQPointContourIntegral(5*i+5).y) );
    }
    break;
  }

  return Result;
  
}

/*! 
 * Integrate the wall shear stress along the contour of this spline interval.
 * 
 * \param SupportCurve the spline associated with this interval. 
 *                     It's required to compute the normals.
 * \param FuncObj  the wall shear stress function (i.e. the function to be integrated)
 *                 the value of this function depends on the location and the tangent direction.
 * \param ResultXdir the result in the X-direction (i.e FuncObj * t_x)
 * \param ResultYdir the result in the Y-direction (i.e FuncObj * t_y)
 *
 * \note The contributions of this interval are ADDED to the result!!!
 */
template<typename FO>
void Spline2DInterval_HO::IntegrateWallShearStressContributions(const Spline2D_HO & SupportCurve, FO FuncObj,
								double & ResultXdir,
								double & ResultYdir,
								double & WettedSurface) const {
  
  double FuncVal(0.0), TempX(0.0), TempY(0.0);
  double const * GQ_Weight(NULL);
  int GQP, index;
  Vector2D Normal, Tangent;

  // Set the weights correctly
  switch(NUMBER_OF_GQP_CONTOURINT){
  case 3:
    GQ_Weight = &GaussQuadratureData::GQ3_Weight[0];
    break;
    
  case 5:
    GQ_Weight = &GaussQuadratureData::GQ5_Weight[0];
    break;
  }
  
  for (int i=0; i<NumOfSubIntervals(); ++i){
    // Calculate the contribution of each subinterval
    
    for (GQP = 1; GQP <= NUMBER_OF_GQP_CONTOURINT; ++GQP){
      // Calculate the contribution of each Gauss-quadrature point
      
      // Determine the index of the current GQP in the array
      index = NUMBER_OF_GQP_CONTOURINT*i+GQP;

      // Determine the normal vector
      Normal = SupportCurve.nSpline(GQPointContourIntegral(index));

      // Determine the tangent vector
      Tangent.x =  Normal.y;
      Tangent.y = -Normal.x;

      // Evaluate weighted function at the current GQP
      FuncVal = GQ_Weight[GQP-1] * FuncObj(GQPointContourIntegral(index),
					   Normal);


      // Update the X projection
      TempX += FuncVal* Tangent.x;

      // Update the Y projection
      TempY += FuncVal* Tangent.y;
    }
  
    // Update final result with the contribution of the current subinterval
    ResultXdir += TempX * IntLength(i+1);
    ResultYdir += TempY * IntLength(i+1);
    WettedSurface += IntLength(i+1);
  }

  GQ_Weight = NULL;
}


#endif /* _HO_SPLINE2DINTERVAL_INCLUDED  */

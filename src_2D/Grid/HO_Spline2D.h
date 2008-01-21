/* \file HO_Spline2D.h
   \brief Header file defining high-order 2D Spline classes. */

#ifndef _HO_SPLINE2D_INCLUDED
#define _HO_SPLINE2D_INCLUDED

/* Include required C++ libraries. */
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <cstring>

/* Using std namespace functions */
using namespace std;

/* Include CFFC header files */
#include "../Math/Math.h"	    // Include math macro header file.
#include "../Math/Vector2D.h"       // Include 2D vector header file
#include "../Math/LinkedList.h"     // Include linked list header file
#include "../Utilities/EpsilonTol.h" // Include error tolerances header file
#include "../MPI/MPI.h"             // Include MPI header file
#include "HO_Node2D.h"		    // Include high-order 2D node header file


/* Define the 2D spline types. */

#define	SPLINE2D_CONSTANT        0
#define	SPLINE2D_LINEAR          1
#define	SPLINE2D_QUADRATIC       2
#define	SPLINE2D_CUBIC           3
#define	SPLINE2D_QUARTIC         4
#define	SPLINE2D_QUINTIC         5

/* Define the 2D spline point types. */

#define	SPLINE2D_POINT_NORMAL        0
#define	SPLINE2D_POINT_SHARP_CORNER  1


/* Define the high-order 2D spline class. */

/* Define the values that FluxMethod can take */
enum MethodsOfFluxCalculation { SolveRiemannProblem = 0, ReconstructionBasedFlux = 1};

/*!
 * \class Spline2D_HO
 *
 * @brief Spline fit in 2 space dimensions.
 *
 * \verbatim
 * Member functions
 *      type       -- Return spline type.
 *      np         -- Return number of points defining the spline.
 *      Xp         -- Return array of coordinates of points defining 
 *                    the spline.
 *      sp         -- Return array of path lengths of points defining 
 *                    the spline.
 *      tp         -- Return array of point types (sharp corner or 
 *                    normal) for the spline.
 *      bc         -- Return array of boundary condition types for 
 *                    points defining spline.
 *      allocate   -- Allocate memory for the spline.
 *      deallocate -- Deallocate memory for the spline.
 *      settype    -- Sets the spline type.
 *      setBCtype  -- Sets the spline boundary condition type (all 
 *                    spline points).
 *      pathlength -- Evaluate the path length at each point defining
 *                    the spline.
 *
 * Member operators
 *      S -- a 2D spline
 *      V -- a 2D vector
 *      a -- a scalar (double)
 *
 * S = S;
 * S = S + S; (concatenate two splines)
 * S = S + V; (shift location of spline)
 * S = S - V; (shift location of spline)
 * S = a * S; (scale spline)
 * S = S * a; (scale spline)
 * S = S / a; (scale spline)
 * S = S ^ a; (rotate spline)
 * cout << S; (output function)
 * cin  >> S; (input function)
 * \endverbatim
 */
class Spline2D_HO{
public:
  int       type; //!< Type of interpolation used in spline definition.
  int         np; //!< Number of points used to define spline.
  Vector2D   *Xp; //!< Array of 2D vectors containing (x,y) 
  //!< coordinates of points used to define spline.
  double     *sp; //!< Array of path-length coordinate of points used
  //!< to define spline
  int        *tp; //!< Array of point types (sharp corner or normal)
  //!< for the spline
  int        *bc; //!< Array of boundary condition types for points 
  //!< defining the spline


  //! @name Constructors, destructor and assignment operator
  //@{
  //! Default constructor.
  Spline2D_HO(void);

  //! Destructor.
  ~Spline2D_HO(void){ deallocate(); }

  //! Assignment operator.
  Spline2D_HO& operator= (const Spline2D_HO &S);
  //@}

  //! @name Allocation and deallocation.
  //@{
  //! Allocate memory for the spline.
  void allocate(const int &N);
  void allocate(const int &N, const int &TypeOfAllPoints);

  //! Deallocate memory for the spline.
  void deallocate(void);
  //@}

  //! Set the spline type.
  void settype(const int &spline_type);

  //! Set the spline boundary condition type.
  void setBCtype(const int &bc_type);

  //! Set the spline point type.
  void settptype(const int &tp_type);

  //! Evaluate the pathlength parameter at each point.
  void pathlength(void);
  //! Evaluate the pathlength parameter at each point.
  friend void pathlength(Spline2D_HO &S){ return S.pathlength(); }

  //! @name Functions to define the flux calculation method.
  //@{
  //! Get the flux calculation method for the current spline.
  const int & getFluxCalcMethod(void) const { return FluxMethod; } 
  //! Set the flux calculation method for the current spline
  void setFluxCalcMethod(const int & Method){ FluxMethod = Method; }
  //@}

  //! @name Output functions 
  //@{
  void OutputTecplot(std::ostream &output_file, const int NumberOfPoints = 100);
  //@}


  //! @name Contour integration functions.
  //@{
  double PolynomOrderIntegration(const Vector2D &StartPoint, const Vector2D &EndPoint,
				 const Vector2D &Centroid, const int &digits,
				 const int &OrderX, const int &OrderY) const;
  double PolynomOrderIntegration(const double &SplinePathStart, const double &SplinePathEnd,
				 const Vector2D &Centroid, const int &digits,
				 const int &OrderX, const int &OrderY) const;
  double PolynomOrderIntegration(const Node2D_HO &StartPoint, const Node2D_HO &EndPoint,
				 const Vector2D &Centroid, const int &digits,
				 const int &OrderX, const int &OrderY) const;
  double BasicOrderIntegration(const Vector2D &StartPoint, const Vector2D &EndPoint,
			       const Vector2D &Centroid, const int &digits,
			       const int &OrderX, const int &OrderY) const;

  double ZeroOrderIntegration(const double &SplinePathStart,
			      const double &SplinePathEnd,
			      const int &digits) const;
  double ZeroOrderIntegration(const Node2D_HO &StartPoint,
			      const Node2D_HO &EndPoint,
			      const int &digits) const;
  double ZeroOrderIntegration(const Vector2D &StartPoint,
			      const Vector2D &EndPoint,
			      const int &digits);
  //@}

  //! @name Normal and tangential vectors at a specified location on the spline
  //@{
  const Vector2D nSpline(const Vector2D &Point) const;
  const Vector2D nSpline(const Vector2D &Point, double & dxds, double & dyds) const;
  const Vector2D nSpline_SegmentBased(const Vector2D &Point) const;

  const Vector2D tSpline(const Vector2D &Point, const PolynomOrder Order) const;
  const Vector2D tSpline(const Vector2D &Point, const PolynomOrder Order,
			 double & dxds, double & dyds) const;
  const Vector2D tSpline(const double & s, const Vector2D * Xp_used,
			 const double * s_used, const int & NumOfControlPoints,
			 double & dxds, double & dyds) const;
  //@}

  //! @name Spline manipulation functions.
  //@{
  //! Determine the position vector corresponding to a specified path length s along the defined spline.
  const Vector2D Spline(const double &s) const;

  //! Determine the subinterval containing the point of interest
  void find_subinterval(const double &s, int & il, int & ir) const;

  //! Get path length for a given position vector
  double getS(const Vector2D &X) const;
  friend double getS(const Vector2D &X, const Spline2D_HO &S){ return S.getS(X); }
  //@}

  //! @name Spline translation, rotation and scaling.
  //@{
  const Spline2D_HO & translate(const Vector2D &V);
  const Spline2D_HO & rotate(const double &a);
  const Spline2D_HO & scale(const double &a);
  //@}

  //! @name Binary arithmetic operators.
  //@{
  const Spline2D_HO operator + (const Spline2D_HO &S) const;
  const Spline2D_HO operator + (const Vector2D &V) const;
  const Spline2D_HO operator - (const Vector2D &V) const;
  const Spline2D_HO operator * (const double &a) const;
  const Spline2D_HO operator / (const double &a) const;
  const Spline2D_HO operator ^ (const double &a) const;

  friend const Spline2D_HO operator +(const Vector2D &V, Spline2D_HO &S){ return (S + V);}
  friend const Spline2D_HO operator *(const double &a, const Spline2D_HO &S){ return S*a;}

  friend bool operator== (const Spline2D_HO &left, const Spline2D_HO &right);
  friend bool operator!= (const Spline2D_HO &left, const Spline2D_HO &right);
  //@}

  //! Concatenate (combine) spline S1 and S2 and returns 
  // a new spline which is a combination of both.
  friend const Spline2D_HO Concatenate_Splines(const Spline2D_HO &S1,
					       const Spline2D_HO &S2){
    return (S1 + S2);
  }

  //! @name Input-output operators.
  //@{
  friend ostream &operator << (ostream &out_file, const Spline2D_HO &S);
  friend istream &operator >> (istream &in_file, Spline2D_HO &S);
  //@}

private:
  Spline2D_HO(const Spline2D_HO &S); //! Private copy constructor

  int FluxMethod;   //!< variable to set the flux calculation method through the boundary spline

  //! Calculate the normal vector for a defined segment
  const Vector2D NormalVector(const Vector2D &StartPoint, const Vector2D &EndPoint) const;

};

/*
 * Default constructor.
 */ 
inline Spline2D_HO::Spline2D_HO(void)
  : type(SPLINE2D_CONSTANT), np(0), Xp(NULL),
    sp(NULL), tp(NULL), bc(NULL), FluxMethod(SolveRiemannProblem)
{
  // 
}

/*!
 * Copy constructor. It is declared private
 */
inline Spline2D_HO::Spline2D_HO(const Spline2D_HO &S)
  : type(SPLINE2D_CONSTANT), np(0), Xp(NULL),
    sp(NULL), tp(NULL), bc(NULL), FluxMethod(SolveRiemannProblem)
{

  /* allocate memory for the new spline */
  allocate(S.np);
  /* assign the same type */
  type = S.type;
  FluxMethod = S.getFluxCalcMethod();
  
  /* copy the values */
  for (int i=0; i<=np-1; ++i){
    Xp[i] = S.Xp[i];
    sp[i] = S.sp[i];
    tp[i] = S.tp[i];
    bc[i] = S.bc[i]; 
  }
}


/*!
 * Allocate memory.
 *
 * \param N number of spline vertexes
 */
inline void Spline2D_HO::allocate(const int &N) {

  // Check conditions
  assert( N >= 2 );

  // Check if the new required memory has dimensions different than the currently allocated ones
  if (np != N){

    deallocate();

    // allocate new memory
    np = N;
    Xp = new Vector2D[N]; 
    sp = new double[N];
    tp = new int[N];
    bc = new int[N];
  }
}

/*!
 * Allocate memory.
 *
 * \param N number of spline vertexes
 * \param TypeOfAllPoints spline point type 
 */
inline void Spline2D_HO::allocate(const int &N, const int &TypeOfAllPoints) {
  /* allocate memory */
  allocate(N);
  /* initialize all the points to the same type -> SPLINE2D_POINT_NORMAL or SPLINE2D_POINT_SHARP_CORNER */
  for (int i=0; i<=np-1; ++i){
    tp[i] = TypeOfAllPoints;
  }
}

/*
 * Deallocate memory.
 */
void Spline2D_HO::deallocate(void) {

  if(Xp != NULL){		/* Check if memory is allocated */
    delete []Xp; Xp = NULL;
    delete []sp; sp = NULL; 
    delete []tp; tp = NULL;
    delete []bc; bc = NULL; 
    np = 0;
    type = SPLINE2D_CONSTANT;
    FluxMethod = SolveRiemannProblem;
  }
}

/*
 * Set the spline type.
 */
inline void Spline2D_HO::settype(const int &spline_type) {
  type = spline_type;
}

/*
 * Set boundary condition type. 
 */
inline void Spline2D_HO::setBCtype(const int & bc_type) {
   int i; assert( np >= 2 );

   for ( i = 0; i <= np-1; ++i ){
     bc[i] = bc_type;
   }
}

/*
 * Set point type.
 */
inline void Spline2D_HO::settptype(const int &tp_type) {
  int i; assert( np >= 2 );
  
  for ( i = 0; i <= np-1; ++i ){
    tp[i] = tp_type;
  }
}

/*
 * Determine path lengths, sp.
 */
inline void Spline2D_HO::pathlength(void) {
  int i; assert( np >= 2 );
  
  sp[0]=ZERO; 
  for ( i = 1; i <= np-1; ++i ){
    sp[i]=sp[i-1]+abs(Xp[i]-Xp[i-1]);
  }
}

/*
 * Assignment operator
 */
inline Spline2D_HO & Spline2D_HO::operator=(const Spline2D_HO &S){

  // Handle self-assignment
  if(this == &S) return *this;

  /* Allocate memory if there is not enough */
  allocate(S.np);

  /* assign the same type as S */
  type = S.type;
  FluxMethod = S.getFluxCalcMethod();

  /* Copy the values from S */
  for (int i=0; i<=np-1; ++i){
    Xp[i] = S.Xp[i];
    sp[i] = S.sp[i];
    tp[i] = S.tp[i];
    bc[i] = S.bc[i];
  }
  return *this;
}

/*
 * Equal operator
 */
inline bool operator== (const Spline2D_HO &left, const Spline2D_HO &right){
  bool answer = true;
  int i(0);
  
  if (left.np != right.np || left.type != right.type
      || left.getFluxCalcMethod() != right.getFluxCalcMethod()){
    return false;
  }
  while (answer==true && i<right.np){
    answer = ( (left.Xp[i]==right.Xp[i]) && (left.tp[i]==right.tp[i]) &&
	       (left.bc[i]==right.bc[i]) && (left.sp[i]==right.sp[i]) );
    ++i;
  }
  return answer;
}
 
/*
 * Not equal operator
 */
inline bool operator!= (const Spline2D_HO &left, const Spline2D_HO &right){
  return !(left == right);
}

/*
 * Shift operators.
 */
const Spline2D_HO Spline2D_HO::operator + (const Vector2D &V) const {
  return Spline2D_HO(*this).translate(V);
}

/*
 * Shift operator with negative vector
 */
const Spline2D_HO Spline2D_HO::operator - (const Vector2D &V) const{
  return Spline2D_HO(*this).translate(-V);
}

/*
 * Scaling operator.
 */
const Spline2D_HO Spline2D_HO::operator * (const double &a) const{
  return Spline2D_HO(*this).scale(a);
}

/*
 * Scaling with invers number
 */
const Spline2D_HO Spline2D_HO::operator / (const double &a) const{
  return Spline2D_HO(*this).scale(1.0/a);
}

/*
 * Rotation operator.
 */
const Spline2D_HO Spline2D_HO::operator ^ (const double &a) const{
  return Spline2D_HO(*this).rotate(a);
}

/*
 * This wrapper gets the unit normal 
 * vector at a given location.
 *
 * \param Point the point of interest.
 */
inline const Vector2D Spline2D_HO::nSpline(const Vector2D &Point) const {
  
  double dxds,dyds;
  
  // Compute the normal
  return nSpline(Point,dxds,dyds);
}

/*
 * Get the unit normal vector at a given 
 * location based on the value of the tangent.
 */
inline const Vector2D Spline2D_HO::nSpline(const Vector2D &Point,
					   double & dxds, double & dyds) const {

  Vector2D Tangent;
  
  // Compute the tangent in Point
  Tangent = tSpline(Point,CUBIC,dxds,dyds);
  
  // Compute the normal
  return Vector2D(Tangent.y,-Tangent.x);
}

/*
 * This wrapper gets the unit tangent vector
 *  at a given location Point.              
 */
inline const Vector2D Spline2D_HO::tSpline(const Vector2D &Point, const PolynomOrder Order) const{

  double dxds, dyds;
  return tSpline(Point,Order,dxds,dyds);
}

/*
 * Determine the unit normal vector to the line segment
 * specified by the end points. The order in which the 
 * end points are passed is important, because it 
 * determines the orientation of the normal vector!!!
 *
 * \param StartPoint the first segment end point
 * \param EndPoint the second segment end point
*/
inline const Vector2D Spline2D_HO::NormalVector(const Vector2D &StartPoint,
						const Vector2D &EndPoint) const {

  double Magnitude(abs(EndPoint-StartPoint));
  
  if (Magnitude > 0.0){
    return ( Vector2D( (EndPoint.y-StartPoint.y),
		       -(EndPoint.x-StartPoint.x) )/ Magnitude);
  } else {
    return Vector2D(0.0,0.0);
  }
}


/*
 * Output operator.
 */
inline ostream &operator << (ostream &out_file, const Spline2D_HO &S) {
  int i;
  out_file << S.getFluxCalcMethod() << endl;   // output FluxMethod
  out_file << S.type << endl;
  out_file << S.np << endl;
  for ( i = 0; i <= S.np-1; ++i ) {
    out_file << S.Xp[i];
    out_file.setf(ios::scientific);
    out_file << " " << S.tp[i] << " " << S.bc[i] << "\n";
    out_file.unsetf(ios::scientific);
  } /* endfor */
  return (out_file);
}

/*
 * Input operator.
 */
inline istream &operator >> (istream &in_file, Spline2D_HO &S) {
  int i;
  in_file.setf(ios::skipws);
  in_file >> i; S.setFluxCalcMethod(i);   // read FluxMethod
  in_file >> i; S.settype(i);
  in_file >> i;
  /* allocate memory if there isn't enough */
  S.allocate(i);
  in_file.unsetf(ios::skipws);
  for ( i = 0; i <= S.np-1; ++i ) {
    in_file >> S.Xp[i];
    in_file.setf(ios::skipws); 
    in_file >> S.tp[i] >> S.bc[i];
    in_file.unsetf(ios::skipws);
  } /* endfor */
  /* update pathlength */
  S.pathlength();
  return (in_file);
}

#if 0

/********************************************************
 * Spline2D_HO -- External subroutines.                    *
 ********************************************************/

extern Vector2D Spline(const double &s,
	      	       const Spline2D_HO &S);

extern int BCtype(const double &s,
	          const Spline2D_HO &S);

extern void Copy_Spline(Spline2D_HO &S1,
	      	        Spline2D_HO &S2);

extern void Copy_Spline(Spline2D_HO &S1,
	      	        const Spline2D_HO &S2);

extern void Broadcast_Spline(Spline2D_HO &S);

#ifdef _MPI_VERSION
extern void Broadcast_Spline(Spline2D_HO &S,
                             MPI::Intracomm &Communicator,
                             const int Source_CPU);
#endif

extern void Translate_Spline(Spline2D_HO &S,
	      	             const Vector2D &V);

extern void Scale_Spline(Spline2D_HO &S,
	      	         const double &Scaling_Factor);

extern void Rotate_Spline(Spline2D_HO &S,
	      	          const double &Angle);

extern void Reflect_Spline(Spline2D_HO &S);

extern void Reverse_Spline(Spline2D_HO &S);

extern void Create_Spline_Line(Spline2D_HO &Line_Spline,
                               const Vector2D &V1,
			       const Vector2D &V2,
  	                       const int Number_of_Spline_Points);

extern void Create_Spline_Line_Polar_Coordinates(Spline2D_HO &Line_Spline,
						 const double &Inner_Radius,
						 const double &Outer_Radius,
						 const double &Theta,
						 const int Number_of_Spline_Points);

extern void Create_Spline_Circular_Arc(Spline2D_HO &Circle_Spline,
			               const Vector2D &Origin,
				       const double &Radius,
                                       const double &Angle1,
			               const double &Angle2,
  	                               const int Number_of_Spline_Points);

extern void Create_Spline_Ellipsoidal_Arc(Spline2D_HO &Ellipse_Spline,
			                  const Vector2D &Origin,
				          const double &A,
				          const double &B,
                                          const double &Angle1,
			                  const double &Angle2,
  	                                  const int Number_of_Spline_Points);

extern void Create_Spline_NACA_Aerofoil(Spline2D_HO &NACA_Spline,
                                        char *NACA_Aerofoil_Type_ptr,
                                        const double &Chord_Length,
					const int i_Up_All_Low,
  	                                const int Number_of_Spline_Points);

extern void Create_Spline_Bow_Shock(Spline2D_HO &Shock_Spline,
                                    const double &Radius,
				    const double &Mach_Number,
				    const int i_Up_All_Low,
  	                            const int Number_of_Spline_Points);

extern void Create_Spline_Area_Variation(Spline2D_HO &Radius_Spline,
                                         const double &Xup,
				         const double &Xthroat,
                                         const double &Xdown,
				         const double &Rup,
				         const double &Rthroat,
				         const double &Rdown,
					 const int &Nozzle_Type,
  	                                 const int Number_of_Spline_Points);

extern void Create_Spline_Converging_Nozzle(Spline2D_HO &Radius_Spline,
					    const double &Xup,
					    const double &Xthroat,
					    const double &Rup,
					    const double &Rthroat,
					    const int Number_of_Spline_Points);

extern void Create_Spline_Diverging_Nozzle(Spline2D_HO &Radius_Spline,
					   const double &Xthroat,
					   const double &Xdown,
					   const double &Rthroat,
					   const double &Rdown,
					   const int Number_of_Spline_Points);

extern void Create_Spline_Rectangle(Spline2D_HO &Rectangle_Spline,
				    const Vector2D &Origin,
				    const double &Length,
				    const double &Width);

extern LinkedList<Vector2D> getX(const double &y, const Spline2D_HO &S);

extern LinkedList<Vector2D> getY(const double &x, const Spline2D_HO &S);



extern int getBCtype(const Vector2D &X, const Spline2D_HO &S);

extern Vector2D getminX(const Vector2D &X, const Spline2D_HO &S);

extern Vector2D getminY(const Vector2D &X, const Spline2D_HO &S);

extern Vector2D getnormal(const Vector2D &X, const Spline2D_HO &S);

#endif

#endif /* _SPLINE2D_INCLUDED  */

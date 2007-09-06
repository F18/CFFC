/*!\file NumericalLibrary.h
  \brief Header file defining some useful numerical methods/objects. */

#ifndef _NumericalLibrary_INCLUDED
#define _NumericalLibrary_INCLUDED

/* Include required C++ libraries. */
#include <limits>
#include <iostream>

/* Using std namespace functions */
using std::numeric_limits;
using std::cout;
using std::endl;
using std::ostream;

/* Include CFFC header files */
#include "../Utilities/TypeDefinition.h"
#include "Math.h"


/************************************************************************//**
 * \class ErrorFunc                                                    
 * \brief General templated error function.                             
 *                                                                     
 * Computes the absolute value of the error between two specified functions.
 ****************************************************************************/
template<typename FO1, typename FO2>
  class ErrorFunc {
 private:
  FO1 fo1;			/*!< pointer to the first function */
  FO2 fo2;			/*!< pointer to the second function */
  
 public:
  // constructor
  ErrorFunc(FO1 f1, FO2 f2)
    : fo1(f1), fo2(f2){
  }
    
    // "function evaluation" with one parameter
    double operator() (double v) const {
      return fabs(fo1(v) - fo2(v));
    }
};

/************************************************************************//**
 * \fn error_func                                                    
 * \brief Generate the error function objec
 ****************************************************************************/
template<typename FO1, typename FO2> 
  inline
  ErrorFunc<FO1,FO2> error_func(FO1 f1, FO2 f2){
  return ErrorFunc<FO1,FO2> (f1,f2);
}

// Error function specialized for objects with member function: SolutionAtCoordinates()
template<class FunctionType, class ObjectType, class SolutionType>
  class _Error_{
 private:
  FunctionType Ptr_F1;
  ObjectType * Ptr_F2;
  unsigned parameter;
  /* Variables refering to the function domain in 1D */
  double  DomainLimitMin, DomainLimitMax;
  double FunctLimitMin, FunctLimitMax;
  double NewCoordinate;

  /* Variable refering to the function domain in 2D */
  int iCell, jCell;

 public:
  /* Constructor for the 1D problem */
  _Error_(FunctionType F1, ObjectType * F2, const unsigned _parameter_,
	  const double _DomainLimitMin_, const double _DomainLimitMax_,
	  const double _FunctLimitMin_, const double _FunctLimitMax_)
    :Ptr_F1(F1), Ptr_F2(F2), parameter(_parameter_), DomainLimitMin(_DomainLimitMin_),
    DomainLimitMax(_DomainLimitMax_), FunctLimitMin(_FunctLimitMin_), FunctLimitMax(_FunctLimitMax_){ };
    
    /* Constructor for the 2D problem */
    _Error_(FunctionType F1, ObjectType * F2, const int ii, const int jj, const unsigned _parameter_)
      :Ptr_F1(F1), Ptr_F2(F2), parameter(_parameter_), iCell(ii), jCell(jj){ 
    };

      SolutionType operator() (double val1){
	NewCoordinate = ConvertDomain(DomainLimitMin,DomainLimitMax,FunctLimitMin,FunctLimitMax,val1);
	return fabs(Ptr_F1(NewCoordinate) - Ptr_F2->SolutionAtCoordinates(val1,parameter));
      }

      SolutionType operator() (double val1, double val2){
	return fabs(Ptr_F1(val1,val2) - Ptr_F2->SolutionAtCoordinates(iCell, jCell, val1, val2,parameter));
      }

      SolutionType operator() (double val1, double val2, double val3){
	return fabs(Ptr_F1(val1,val2,val3) - Ptr_F2->SolutionAtCoordinates(val1,val2,val3,parameter));
      }
};

// Square Error function specialized for objects with member function: SolutionAtCoordinates()
// returns the square of the error
template<class FunctionType, class ObjectType, class SolutionType>
  class _Error_Square_{
 private:
  FunctionType Ptr_F1;
  ObjectType * Ptr_F2;
  unsigned parameter;
  /* Variables refering to the function domain in 1D */
  double  DomainLimitMin, DomainLimitMax;
  double FunctLimitMin, FunctLimitMax;
  double NewCoordinate;
  double FunctionValue;
  
  /* Variable refering to the function domain in 2D */
  int iCell, jCell;
  
 public:
  /* Constructor for the 1D problem */
  _Error_Square_(FunctionType F1, ObjectType * F2, const unsigned _parameter_,
		 const double _DomainLimitMin_, const double _DomainLimitMax_,
		 const double _FunctLimitMin_, const double _FunctLimitMax_)
    :Ptr_F1(F1), Ptr_F2(F2), parameter(_parameter_), DomainLimitMin(_DomainLimitMin_),
    DomainLimitMax(_DomainLimitMax_), FunctLimitMin(_FunctLimitMin_), FunctLimitMax(_FunctLimitMax_){ 
  };

    /* Constructor for the 2D problem */
    _Error_Square_(FunctionType F1, ObjectType * F2, const int ii, const int jj, const unsigned _parameter_)
      :Ptr_F1(F1), Ptr_F2(F2), parameter(_parameter_), iCell(ii), jCell(jj){ 
    };

      SolutionType operator() (double val1){
	NewCoordinate = ConvertDomain(DomainLimitMin,DomainLimitMax,FunctLimitMin,FunctLimitMax,val1);
	FunctionValue = Ptr_F1(NewCoordinate) - Ptr_F2->SolutionAtCoordinates(val1,parameter);
	return FunctionValue*FunctionValue;
      }

      SolutionType operator() (double val1, double val2){
	FunctionValue = Ptr_F1(val1,val2) - Ptr_F2->SolutionAtCoordinates(iCell, jCell, val1, val2,parameter);
	return FunctionValue*FunctionValue;
      }

      SolutionType operator() (double val1, double val2, double val3){
	FunctionValue = Ptr_F1(val1,val2,val3) - Ptr_F2->SolutionAtCoordinates(val1,val2,val3,parameter);
	return FunctionValue*FunctionValue;
      }
};

// Error function specialized for objects with member function: SolutionAtCoordinates_PWL()
template<class FunctionType, class ObjectType, class SolutionType>
  class _Error_PWL_{
 private:
  FunctionType Ptr_F1;
  ObjectType * Ptr_F2;
  unsigned parameter;
  /* Variables refering to the function domain in 1D */
  double  DomainLimitMin, DomainLimitMax;
  double FunctLimitMin, FunctLimitMax;
  double NewCoordinate;

  /* Variable refering to the function domain in 2D */
  int iCell, jCell;

 public:
  /* Constructor for the 1D problem */
  _Error_PWL_(FunctionType F1, ObjectType * F2, const unsigned _parameter_,
	      const double _DomainLimitMin_, const double _DomainLimitMax_,
	      const double _FunctLimitMin_, const double _FunctLimitMax_)
    :Ptr_F1(F1), Ptr_F2(F2), parameter(_parameter_), DomainLimitMin(_DomainLimitMin_),
    DomainLimitMax(_DomainLimitMax_), FunctLimitMin(_FunctLimitMin_), FunctLimitMax(_FunctLimitMax_){ 
  };

    /* Constructor for the 2D problem */
    _Error_PWL_(FunctionType F1, ObjectType * F2, const int ii, const int jj, const unsigned _parameter_)
      :Ptr_F1(F1), Ptr_F2(F2), parameter(_parameter_), iCell(ii), jCell(jj){ 
    };

      SolutionType operator() (double val1){
	NewCoordinate = ConvertDomain(DomainLimitMin,DomainLimitMax,FunctLimitMin,FunctLimitMax,val1);
	return fabs(Ptr_F1(NewCoordinate) - Ptr_F2->SolutionAtCoordinates_PWL(val1,parameter));
      }
  
      SolutionType operator() (double val1, double val2){
	return fabs(Ptr_F1(val1,val2) - Ptr_F2->SolutionAtCoordinates_PWL(iCell, jCell, val1, val2,parameter));
      }

      SolutionType operator() (double val1, double val2, double val3){
	return fabs(Ptr_F1(val1,val2,val3) - Ptr_F2->SolutionAtCoordinates_PWL(val1,val2,val3,parameter));
      }
};

// Square Error function specialized for objects with member function: SolutionAtCoordinates_PWL()
// returns the square of the error
template<class FunctionType, class ObjectType, class SolutionType>
  class _Error_Square_PWL_{
 private:
  FunctionType Ptr_F1;
  ObjectType * Ptr_F2;
  unsigned parameter;
  /* Variables refering to the function domain in 1D */
  double  DomainLimitMin, DomainLimitMax;
  double FunctLimitMin, FunctLimitMax;
  double NewCoordinate;
  double FunctionValue;
  
  /* Variable refering to the function domain in 2D */
  int iCell, jCell;
  
 public:
  /* Constructor for the 1D problem */
  _Error_Square_PWL_(FunctionType F1, ObjectType * F2, const unsigned _parameter_,
		     const double _DomainLimitMin_, const double _DomainLimitMax_,
		     const double _FunctLimitMin_, const double _FunctLimitMax_)
    :Ptr_F1(F1), Ptr_F2(F2), parameter(_parameter_), DomainLimitMin(_DomainLimitMin_),
    DomainLimitMax(_DomainLimitMax_), FunctLimitMin(_FunctLimitMin_), FunctLimitMax(_FunctLimitMax_){ 
  };

    /* Constructor for the 2D problem */
    _Error_Square_PWL_(FunctionType F1, ObjectType * F2, const int ii, const int jj, const unsigned _parameter_)
      :Ptr_F1(F1), Ptr_F2(F2), parameter(_parameter_), iCell(ii), jCell(jj){ 
    };
  
      SolutionType operator() (double val1){
	NewCoordinate = ConvertDomain(DomainLimitMin,DomainLimitMax,FunctLimitMin,FunctLimitMax,val1);
	FunctionValue = Ptr_F1(NewCoordinate) - Ptr_F2->SolutionAtCoordinates_PWL(val1,parameter);
	return FunctionValue*FunctionValue;
      }

      SolutionType operator() (double val1, double val2){
	FunctionValue = Ptr_F1(val1,val2) - Ptr_F2->SolutionAtCoordinates_PWL(iCell, jCell, val1, val2,parameter);
	return FunctionValue*FunctionValue;
      }

      SolutionType operator() (double val1, double val2, double val3){
	FunctionValue = Ptr_F1(val1,val2,val3) - Ptr_F2->SolutionAtCoordinates_PWL(val1,val2,val3,parameter);
	return FunctionValue*FunctionValue;
      }
};

// class Entropy Error:
// This function computes the error in the entropy prediction at a certain
// location for isentropic flows.
// The entropy error is estimated based on the numerical solution and the free 
// stream pressure and density values. 
// This function is specialized for objects which have a member function: SolutionAt()
template<class ObjectType>
class Entropy_Error_{
 public:
  typedef typename ObjectType::SolutionType SolType;

 private:
  ObjectType * Ptr_F2;
  SolType Solution;

  /* Reference State --> Pressure & Density */
  double  Pressure, Density;
  /* Variables storing the computational domain in 2D */
  int iCell, jCell;

 public:
  /* Constructor for the 1D problem */

  /* Constructor for the 2D problem */
  Entropy_Error_(ObjectType * F2, const int ii, const int jj,
		 const double Pressure_Reference, const double Density_Reference)
    :Ptr_F2(F2), Pressure(Pressure_Reference), Density(Density_Reference), iCell(ii), jCell(jj){ 
  };

    // operator () --> compute entropy error at a particular location
    double operator() (double val1, double val2){
      Solution = Ptr_F2->SolutionAtCoordinates(iCell, jCell, val1, val2);
      return (Solution[4]/Pressure)/pow(Solution[1]/Density, Solution.g) - 1.0;
    }
};

// class Entropy Error_Square:
// This function computes the error in the entropy prediction at a certain
// location for isentropic flows.
// The entropy error is estimated based on the numerical solution and the free 
// stream pressure and density values. 
// This function is specialized for objects which have a member function: SolutionAt()
template<class ObjectType>
class Entropy_Error_Square{
 public:
  typedef typename ObjectType::SolutionType SolType;

 private:
  ObjectType * Ptr_F2;
  SolType Solution;
  double FunctionValue;

  /* Reference State --> Pressure & Density */
  double  Pressure, Density;
  /* Variables storing the computational domain in 2D */
  int iCell, jCell;

 public:
  /* Constructor for the 1D problem */

  /* Constructor for the 2D problem */
  Entropy_Error_Square(ObjectType * F2, const int ii, const int jj,
		       const double Pressure_Reference, const double Density_Reference)
    :Ptr_F2(F2), Pressure(Pressure_Reference), Density(Density_Reference),  iCell(ii), jCell(jj){ 
  };

    // operator ()  --> compute the square of the entropy error at a particular location
    double operator() (double val1, double val2){
      Solution = Ptr_F2->SolutionAtCoordinates(iCell, jCell, val1, val2);
      FunctionValue = (Solution[4]/Pressure)/pow(Solution[1]/Density, Solution.g) - 1.0;
      return FunctionValue*FunctionValue;
    }
};

// Polynomial reconstructions in 2D specialized for objects with member function: SolutionAt()
// This function object is used to access the reconstruction of a cell (iCell,jCell) by using only the
// x-coordinate and y-coordinate
// Thus, the reconstruction function can be integrated using any standard integration subroutine
template<class ObjectType, class SolutionType>
  class _ReconstructFunctor_{
 private:
  ObjectType * Reconstruction;

  /* Variable refering to the function domain in 2D */
  int iCell, jCell;

 public:
  /* Constructor for the 2D problem */
  _ReconstructFunctor_(ObjectType * Rec, const int ii, const int jj)
    : Reconstruction(Rec), iCell(ii), jCell(jj){ 
  };

    SolutionType operator() (double val1, double val2){
      return Reconstruction->SolutionAtCoordinates(iCell, jCell, val1, val2);
    }
};

// qgaus5 Template
template<class FunctionType, class ReturnType>
  inline ReturnType qgaus5(FunctionType func, double a, double b, ReturnType _dummy_param){
  /* returns the integral of the function "func" between a and b,
     by five-point Gauss-Legendre integration: the function is evaluated exactly
     five times at interior points in the range of integration.
     The integral is exact for polynomials up to order of 9. 
     Implementation based on the subroutine from Numerical Recipes.

     The template function takes two parameters: "FunctionType" and "ReturnType"
     FunctionType: is used to determine the type of pointer function that is passed
     ReturnType: is used to determine the type of return type

     If a function that takes one double and returns a solution state is provided, this subroutine will
     return the appropiate integration of all the parameters in the specified domain.

     REQUIREMENTS:
     The state class should provide the following operators:
     "operator +"  between 2 states
     "operator *"  between a state and a scalar
     "operator /"  between two states
     "operator <"  between two states
     "a constructor which takes a double" -> this will initialize all the parameters to the taken value

     _dummy_param: is a variable used only for providing the return type.
  */

  int j;
  long double xr,xm,dx;
  ReturnType sum;
  static double x[]={0.0, 0.53846931010568311, 0.90617984593866396};
  static double w[]={0.56888888888888889, 0.47862867049936647, 0.23692688505618917};

  xm=0.5*(b+a);
  xr=0.5*(b-a);
  sum=(ReturnType)(w[0]*(func(xm)));
  for (j=1;j<=2;++j) {
    dx=xr*x[j];
    sum = sum + (ReturnType)(w[j]*(func(xm+dx)+func(xm-dx)));
  }
  sum = (ReturnType)(sum*xr);
  return sum;
};

template<class FunctionType, class ReturnType>
  inline ReturnType AdaptiveGaussianQuadrature(FunctionType func, double StartPoint, double EndPoint,
					       double InitialDistance, int digits,
					       ReturnType & gauss5, ReturnType & RERR,
					       int & FunctionEvaluations, const int & MaxFunctEval)
  throw(too_short_interval,maximum_exceeded){

  long double EPS = 0.5*pow(10.0,1.0-digits); // accuracy: --> based on precision 
  // i.e exact number of digits 
  ReturnType G(0.0);			// Integral value
  long double SubDivisionPoint;

  if ( 1.0 + 0.005*fabs(StartPoint-EndPoint)/InitialDistance <= 1.0){
    throw too_short_interval();
  } else if (FunctionEvaluations > MaxFunctEval) {
    throw maximum_exceeded();
  }

  /* evaluate once */
  gauss5 = qgaus5(func,StartPoint,EndPoint,gauss5);
  
  /* evaluate more accurate */
  SubDivisionPoint = StartPoint + 0.5*(EndPoint - StartPoint);
  G = qgaus5(func,StartPoint,SubDivisionPoint,gauss5) + qgaus5(func,SubDivisionPoint,EndPoint,gauss5);
  
  /* update count of function evalutions */
  FunctionEvaluations += 15;
  
  /* check relative error */
  RERR = fabs(G - gauss5)/( (ReturnType)1.0 + fabs(G) );
  
  if (RERR <= (ReturnType)EPS){
    return G;
  }
  else {
    G = AdaptiveGaussianQuadrature(func,StartPoint,SubDivisionPoint,InitialDistance,digits,
				   gauss5, RERR, FunctionEvaluations,MaxFunctEval) +
      AdaptiveGaussianQuadrature(func,SubDivisionPoint,EndPoint,InitialDistance,digits,
				 gauss5, RERR, FunctionEvaluations,MaxFunctEval);
  }

  return G;
};

template <class T>
T sign(T X){
  if( X != T(0.0)){
    return X/fabs(X);
  } else {
    return T(0.0);
  }
}

#define MaxAllowedFunctionEvaluation 8000000 /* maximum allowed function evaluations for the adaptive integration */

template <class FunctionType, class ReturnType>
  inline ReturnType adaptlobstp(FunctionType func, double a, double b, const ReturnType fa, const ReturnType fb, 
				const ReturnType is, int &FunctionEvaluations, int & WriteMessage){

  // %ADAPTLOBSTP  Recursive function used by ADAPTLOB.
  // %
  // %   Q = ADAPTLOBSTP('Func',A,B,FA,FB,IS,WriteMessage) tries to
  // %   approximate the integral of F(X) from A to B to
  // %   an appropriate relative error. The argument 'Func' is
  // %   a string containing the name of f.  The remaining
  // %   arguments are generated by ADAPTLOB or by recursion.
  // %
  // %   See also ADAPTLOB.
  // % This algorithm follows the work of Walter Gautschi, and is similar to the Matlab's implementation

  // Lucian Ivan, 31/09/2005 
  const ReturnType eps(numeric_limits<double>::epsilon());
  double h,m,alpha,beta, mll, ml, mr, mrr;
  ReturnType i2, i1;
  ReturnType y1,y2,y3,y4,y5;

  h=(b-a)*0.5;
  m=(a+b)*0.5;
  alpha = 0.81649658092772603446; //alpha=sqrt(2.0/3.0);
  beta  = 0.4472135954999579277; //beta=1.0/sqrt(5.0);

  mll=m-alpha*h; 
  ml=m-beta*h;
  mr=m+beta*h;
  mrr=m+alpha*h;

  y1 = func(mll); y2 = func(ml) ; y3 = func(m);
  y4 = func(mr); y5 = func(mrr) ;

  /* This subroutine will add 5 function evaluations */
  FunctionEvaluations += 5;

  i2=(h/6.0)*(fa+fb+5.0*(y2+y4));
  i1=(h/1470.0)*(77.0*(fa+fb)+432.0*(y1+y5)+625.0*(y2+y4)+672.0*y3);

  if ( (is+(i1-i2)==is) || (mll<=a) || (b<=mrr) || (FunctionEvaluations>MaxAllowedFunctionEvaluation) || (fabs(i1-i2)<=eps) ){
    if ( ( (m <= a) || (b<=m) ) && (WriteMessage == 0)){
      cout << "\nWarning Integration Subroutine: Interval contains no more machine number.\n"
	   << "Required tolerance may not be met.\n";
      WriteMessage = 1;
    }
    if ((FunctionEvaluations>MaxAllowedFunctionEvaluation) && (WriteMessage == 0)){
      cout << "\nWarning Integration Subroutine: Maximum function count exceeded (" << MaxAllowedFunctionEvaluation
	   << "); singularity likely.\n"
	   << "Required tolerance may not be met.\n";
      WriteMessage = 1;
    }
    return i1;
  } else {
    return (adaptlobstp(func,a,mll,fa,y1,is,FunctionEvaluations,WriteMessage)+
	    adaptlobstp(func,mll,ml,y1,y2,is,FunctionEvaluations,WriteMessage)+
	    adaptlobstp(func,ml,m,y2,y3,is,FunctionEvaluations,WriteMessage)+
	    adaptlobstp(func,m,mr,y3,y4,is,FunctionEvaluations,WriteMessage)+
	    adaptlobstp(func,mr,mrr,y4,y5,is,FunctionEvaluations,WriteMessage)+
	    adaptlobstp(func,mrr,b,y5,fb,is,FunctionEvaluations,WriteMessage));
  }
}

template <class FunctionType, class ReturnType>
  inline ReturnType GaussLobattoAdaptiveQuadrature(FunctionType func, double StartPoint, double EndPoint,
						   ReturnType _dummy_param, int &WriteMessage,
						   int digits = numeric_limits<double>::digits10){

  /* %GaussLobattoAdaptiveQuadrature:
     %   Numerically evaluate integrals using adaptive
     %   Lobatto rule.
     %
     %   Q=GaussLobattoAdaptiveQuadrature('F',A,B,Q,MessageFlag)
     %   approximates the integral of
     %   F(X) from A to B to machine precision.  'F' is a
     %   string containing the name of the function. The function F
     %   can return a state of multiple values.
     %
     %   Q=ADAPTLOB('F',A,B,Q,MessageFlag,digits) integrates to a precision
     %   specified by the number of digits.
     %
     %   Lucian Ivan, 02/11/2005 
     %   The implementation follows Walter Gautschi's work.
     %   see http://www.inf.ethz.ch/personal/gander/adaptlob.m
     %   Reference: Gander, Computermathematik, Birkhaeuser, 1992. */

  ReturnType tol(0.5*pow(10.0,-digits));
  const ReturnType eps(numeric_limits<double>::epsilon());

  double a,b;
  double IntSign;
  int FunctionEvaluations(0);

  // Check if StartPoint < EndPoint
  if (StartPoint < EndPoint){
    a = StartPoint;
    b = EndPoint;
    IntSign = 1.0;
  } else if (StartPoint == EndPoint) {
    std::cout << "Integration Warning: The interval for integration has ZERO length!\n";
    return ReturnType(0.0);
  } else {
    b = StartPoint;
    a = EndPoint;
    IntSign = -1.0;
  }

  long double m, h, alpha, beta, x1, x2, x3;
  ReturnType i2, i1, is, erri1, erri2;
  ReturnType y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13;
  ReturnType s;

  if (tol < eps){
    tol = eps;
  }

  m=(a+b)*0.5;
  h=(b-a)*0.5;
  alpha = 0.81649658092772603446; //alpha=sqrt(2.0/3.0);
  beta  = 0.4472135954999579277; //beta=1.0/sqrt(5.0);
  x1=.942882415695480; x2=.641853342345781;
  x3=.236383199662150;

  y1 = func(a);  y2 = func(m-x1*h);   y3 = func(m-alpha*h);   y4 = func(m-x2*h);
  y5 = func(m-beta*h); y6 = func(m-x3*h); y7 = func(m); y8 = func(m+x3*h);
  y9 = func(m+beta*h); y10 = func(m+x2*h); y11 = func(m+alpha*h); y12 = func(m+x1*h);
  y13 = func(b);
  
  i2=(h/6.0)*(y1+y13+5*(y5+y9));
  i1=(h/1470.0)*(77.0*(y1+y13)+432.0*(y3+y11)+ 625.0*(y5+y9)+672.0*y7);
  is=h*(.0158271919734802*(y1+y13)+.0942738402188500*(y2+y12)+.155071987336585*(y3+y11)+
	.188821573960182*(y4+y10)+.199773405226859*(y5+y9)+.224926465333340*(y6+y8)+
	.242611071901408*y7);

  s=sign(is);
  if(s==0){
    s=1;
  }

  /* Relaxation of the tolerance */
  erri1=fabs(i1-is);
  erri2=fabs(i2-is);

  ReturnType R(1);
  ReturnType Zero(0.0);

  if (erri2 != Zero){
    R=erri1/erri2;
  }

  if (R>Zero && R<ReturnType(1.0)){
    tol=tol/R;
  }

  is=s*fabs(is)*tol/eps;

  if (is==0) {
    is=b-a;
  }

  FunctionEvaluations = 13;

  return IntSign*adaptlobstp(func,a,b,y1,y13,is,FunctionEvaluations,WriteMessage);
}

template <class FunctionType, class ReturnType>
  inline ReturnType GaussLobattoQuadrature(FunctionType func, double StartPoint, double EndPoint,
					   ReturnType _dummy_param){

  /* %GaussLobattoQuadrature:
     %   Numerically evaluate integrals using a 13 points Lobatto rule
     %   which has the degree 18.
     %   
     %   Q=GaussLobattoQuadrature('F',A,B,Q)
     %   approximates the integral of
     %   F(X) from A to B. For polynomial functions up to a degree of 18 
     %   the result is exact.
     %   'F' is a string containing the name of the function. The function F
     %   can return a state of multiple values.
     %
     %   Lucian Ivan, 07/11/2005 
     %   The implementation follows Walter Gautschi's work.
     %   see http://www.inf.ethz.ch/personal/gander/adaptlob.m
     %   Reference: Gander, Computermathematik, Birkhaeuser, 1992. */

  double a,b;
  double IntSign;

  // Check if StartPoint < EndPoint
  if (StartPoint < EndPoint){
    a = StartPoint;
    b = EndPoint;
    IntSign = 1.0;
  } else if (StartPoint == EndPoint) {
    std::cout << "Integration Warning: The interval for integration has ZERO length!\n";
    return ReturnType(0.0);
  } else {
    b = StartPoint;
    a = EndPoint;
    IntSign = -1.0;
  }

  long double m, h, alpha, beta, x1, x2, x3;
  ReturnType y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13;

  m=(a+b)*0.5; h=(b-a)*0.5;
  alpha = 0.81649658092772603446; //alpha=sqrt(2.0/3.0);
  beta  = 0.4472135954999579277; //beta=1.0/sqrt(5.0);
  x1=.942882415695480; x2=.641853342345781;
  x3=.236383199662150;

  y1 = func(a); y2 = func(m-x1*h); y3 = func(m-alpha*h); y4 = func(m-x2*h);
  y5 = func(m-beta*h); y6 = func(m-x3*h); y7 = func(m); y8 = func(m+x3*h);
  y9 = func(m+beta*h); y10 = func(m+x2*h); y11 = func(m+alpha*h); y12 = func(m+x1*h);
  y13 = func(b);
  
  return IntSign*h*(.0158271919734802*(y1+y13)+.0942738402188500*(y2+y12)+.155071987336585*(y3+y11)+
		    .188821573960182*(y4+y10)+.199773405226859*(y5+y9)+.224926465333340*(y6+y8)+
		    .242611071901408*y7);
}

template <class FunctionType, class ReturnType>
  inline ReturnType AdaptiveGaussianQuadrature(FunctionType func, double StartPoint, double EndPoint,
					       ReturnType _dummy_param, int digits = numeric_limits<double>::digits10){

  /**********************************************************************************************************
        The template function takes two parameters: "FunctionType" and "ReturnType"
	FunctionType: is used to determine the type of pointer function that is passed
 	ReturnType: is used to determine the type of return type
	
	If a function that takes one double and returns a solution state is provided, this subroutine will
 	return the appropiate integration of all the parameters in the specified domain.

 	REQUIREMENTS:
 	The state class should provide the following operators:
 	"operator +"  between 2 states
 	"operator *"  between a state and a scalar
 	"operator /"  between two states
 	"operator <"  between two states
 	"operator <="  between two states
 	"operator >"  between two states
 	"operator >="  between two states
 	"a constructor which takes a double" -> this will initialize all the parameters to the value that is passed

	_dummy_param: is a variable used only for providing the return type.
  *********************************************************************************************************/

  int WriteMessage = 0; 	/* this flag makes sure that the error message is written only once */
  return GaussLobattoAdaptiveQuadrature(func,StartPoint,EndPoint,_dummy_param,WriteMessage,digits);
}

/**************************************************************************
 * FunctorX: defines a new function which has the X variable fixed to Val *
 *************************************************************************/
template<class FunctionType, class SolutionType>
  class FunctorX
{
 private:
  FunctionType Ptr_F;		/* this function takes two or three parameters x, y & z */
  const double Val;

 public:
  FunctorX( FunctionType F, const double & _Val_)
    : Ptr_F(F), Val(_Val_){
  };

    SolutionType operator() (const double & y){
      return Ptr_F(Val,y);	/* return the value of Ptr_F in point (Val,y) */
    };

    SolutionType operator() (const double & y, const double & z){
      return Ptr_F(Val,y,z);       /* return the value of Ptr_F in point (Val,y,z) */
    }
};

/**************************************************************************
 * FunctorY: defines a new function which has the Y variable fixed to Val *
 *************************************************************************/
template<class FunctionType, class SolutionType>
  class FunctorY
{
 private:
  FunctionType Ptr_F;		/* this function takes two or three parameters x, y & z */
  const double Val;

 public:
  FunctorY( FunctionType F, const double & _Val_)
    : Ptr_F(F), Val(_Val_){
  };

    SolutionType operator() (const double & x){
      return Ptr_F(x,Val);	/* return the value of Ptr_F in point (x,Val) */
    };

    SolutionType operator() (const double & x, const double & z){
      return Ptr_F(x,Val,z);       /* return the value of Ptr_F in point (x,Val,z) */
    }
};

/**************************************************************************
 * FunctorZ: defines a new function which has the Z variable fixed to Val *
 *************************************************************************/
template<class FunctionType, class SolutionType>
  class FunctorZ
{
 private:
  FunctionType Ptr_F;		/* this function takes three parameters x,y & z */
  const double Val;

 public:
  FunctorZ( FunctionType F, const double & _Val_)
    : Ptr_F(F), Val(_Val_){
  };

    SolutionType operator() (const double & x, const double & y){
      return Ptr_F(x,y,Val);	/* return the value of Ptr_F in point (x,y,Val) */
    };

};

/***************************************************************************
 * FunctorXZ: defines a new function which has the X and Z variables fixed *
 *************************************************************************/
template<class FunctionType, class SolutionType>
  class FunctorXZ
{
 private:
  FunctionType Ptr_F;		/* this function takes three parameters x,y & z */
  const double ValX;
  const double ValZ;

 public:
  FunctorXZ( FunctionType F, const double & _Val_X_, const double & _Val_Z_)
    : Ptr_F(F), ValX(_Val_X_), ValZ(_Val_Z_){
  };

    SolutionType operator() ( const double & y){
      return Ptr_F(ValX,y,ValZ);		/* return the value of Ptr_F in point (ValX,y,ValZ) */
    };

};

/***************************************************************************
 * FunctorXY: defines a new function which has the X and Y variables fixed *
 *************************************************************************/
template<class FunctionType, class SolutionType>
  class FunctorXY
{
 private:
  FunctionType Ptr_F;		/* this function takes three parameters x,y & z */
  const double ValX;
  const double ValY;

 public:
  FunctorXY( FunctionType F, const double & _Val_X_, const double & _Val_Y_)
    : Ptr_F(F), ValX(_Val_X_), ValY(_Val_Y_){
  };

    SolutionType operator() ( const double & z){
      return Ptr_F(ValX,ValY,z);  	/* return the value of Ptr_F in point (ValX,ValY,z) */
    };

};

/***************************************************************************
 * FunctorYZ: defines a new function which has the Y and Z variables fixed *
 *************************************************************************/
template<class FunctionType, class SolutionType>
  class FunctorYZ
{
 private:
  FunctionType Ptr_F;		/* this function takes three parameters x,y & z */
  const double ValY;
  const double ValZ;

 public:
  FunctorYZ( FunctionType F, const double & _Val_Y_, const double & _Val_Z_)
    : Ptr_F(F), ValY(_Val_Y_), ValZ(_Val_Z_){
  };

    SolutionType operator() ( const double & x){
      return Ptr_F(x,ValY,ValZ);		/* return the value of Ptr_F in point (x,ValY,ValZ) */
    };

};


/******************************************************************************
 * Antiderivative: defines a class which stores the pointer to a one variable *
 *                 function and the limits it is integrated to.               * 
 *                 operator(double ) returns the value of the integral        *
 *****************************************************************************/
template<class FunctionType, class SolutionType>
  class Antiderivative
{
 private:
  FunctionType Ptr_F;
  const SolutionType DummyParam;
  const double LeftLimit;
  const double RightLimit;
  const int ExactDigits;
  int * WriteMessage; 
 public:
  Antiderivative(FunctionType F, const double & StartY, const double & EndY, int & _WriteMessage_, const int & digits)
    : Ptr_F(F), DummyParam(0.0), LeftLimit(StartY), RightLimit(EndY), ExactDigits(digits) , WriteMessage(& _WriteMessage_) {
  };

    /* Implementation for 2D */
    SolutionType operator() (const double & x){
      FunctorX<FunctionType,SolutionType> Functor(Ptr_F,x);
      return GaussLobattoAdaptiveQuadrature(Functor, LeftLimit, RightLimit, DummyParam, *WriteMessage, ExactDigits);
    }
  
    /* Implementation for 3D */
    SolutionType operator() (const double & x, const double & y){
      FunctorXY<FunctionType,SolutionType> Functor(Ptr_F,x,y);
      return GaussLobattoAdaptiveQuadrature(Functor, LeftLimit, RightLimit, DummyParam, *WriteMessage, ExactDigits);
    }
};

template<class FunctionType, class SolutionType>
  class AntiderivativeGL 		/* customized antiderivative for the GaussLobatto integration (non adaptive!!!) */
{
 private:
  FunctionType Ptr_F;
  const SolutionType DummyParam;
  const double LeftLimit;
  const double RightLimit;
 public:
  AntiderivativeGL(FunctionType F, const double & StartY, const double & EndY)
    : Ptr_F(F), DummyParam(0.0), LeftLimit(StartY), RightLimit(EndY) {
  };

    /* Implementation for 2D */
    SolutionType operator() (const double & x){
      FunctorX<FunctionType,SolutionType> Functor(Ptr_F,x);
      return GaussLobattoQuadrature(Functor, LeftLimit, RightLimit, DummyParam);
    }
  
    /* Implementation for 3D */
    SolutionType operator() (const double & x, const double & y){
      FunctorXY<FunctionType,SolutionType> Functor(Ptr_F,x,y);
      return GaussLobattoQuadrature(Functor, LeftLimit, RightLimit, DummyParam);
    }
};

template<class FunctionType, class SolutionType>
  class AntiderivativeG5 		/* customized antiderivative for the Gauss integration with 5 points (non adaptive!!!) */
{
 private:
  FunctionType Ptr_F;
  const SolutionType DummyParam;
  const double LeftLimit;
  const double RightLimit;
 public:
  AntiderivativeG5(FunctionType F, const double & StartY, const double & EndY)
    : Ptr_F(F), DummyParam(0.0), LeftLimit(StartY), RightLimit(EndY) {
  };

    /* Implementation for 2D */
    SolutionType operator() (const double & x){
      FunctorX<FunctionType,SolutionType> Functor(Ptr_F,x);
      return qgaus5(Functor, LeftLimit, RightLimit, DummyParam);
    }
  
    /* Implementation for 3D */
    SolutionType operator() (const double & x, const double & y){
      FunctorXY<FunctionType,SolutionType> Functor(Ptr_F,x,y);
      return qgaus5(Functor, LeftLimit, RightLimit, DummyParam);
    }
};


/***************************************************************************************
 * AdaptiveGaussianQuadrature Template                                                 *
 * returns the integral of a user-supplied function "func(x,y)" over a two-dimensional * 
 * (2D) rectangular region, specified by the limits StartX, EndX, StartY, EndY.        *
 * Integration is performed by calling AdaptiveGaussianQuadrature for the x and y      *
 * directions.                                                                         *
 **************************************************************************************/
template<class FunctionType, class ReturnType>
  ReturnType AdaptiveGaussianQuadrature(FunctionType func, double StartX, double EndX,
					double StartY, double EndY,  
					int digits, ReturnType _dummy_param){

  int WriteMessage = 0; 	/* this flag makes sure that the error message is written only once */

  /* Create the function inner integral  */
  Antiderivative<FunctionType,ReturnType> PrimitiveY(func,StartY,EndY,WriteMessage,digits);

  /* Integrate the inner integral with respect to X */
  return GaussLobattoAdaptiveQuadrature(PrimitiveY,StartX,EndX,_dummy_param,WriteMessage,digits);
};

template <class FunctionType, class ReturnType>
  inline ReturnType GaussLobattoQuadrature(FunctionType func, double StartX, double EndX,
					   double StartY, double EndY, 
					   ReturnType _dummy_param){

  /* Create the inner integral */
  AntiderivativeGL<FunctionType,ReturnType> PrimitiveY(func,StartY,EndY);

  /* Integrate the inner integral with respect to X */
  return GaussLobattoQuadrature(PrimitiveY,StartX,EndX,_dummy_param);
}

template <class FunctionType, class ReturnType>
  inline ReturnType Gauss5PointQuadrature(FunctionType func, double StartX, double EndX,
					  double StartY, double EndY, 
					  ReturnType _dummy_param){

  /* Create the inner integral */
  AntiderivativeG5<FunctionType,ReturnType> PrimitiveY(func,StartY,EndY);

  /* Integrate the inner integral with respect to X */
  return qgaus5(PrimitiveY,StartX,EndX,_dummy_param);
}

/*******************************************************************************************
 * AdaptiveGaussianQuadrature Template                                                     *
 * returns the integral of a user-supplied function "func(x,y,z)" over a three-dimensional * 
 * (3D) rectangular region, specified by the limits StartX, EndX, StartY, EndY, StartZ,    *
 * EndZ. Integration is performed by calling AdaptiveGaussianQuadrature for the x, y and z *
 * directions.                                                                             *
 ******************************************************************************************/
template<class FunctionType, class ReturnType>
  ReturnType AdaptiveGaussianQuadrature(FunctionType func, const double StartX,const double EndX,
					const double StartY, const double EndY, const double StartZ,
					const double EndZ, int digits, ReturnType _dummy_param){

  int WriteMessage = 0; 	/* this flag makes sure that the error message is written only once */

  /* Create the function to be integrated with Y */
  Antiderivative<FunctionType,ReturnType> PrimitiveZ(func,StartZ,EndZ,WriteMessage,digits);

  /* Create the function to be integrated with X */
  Antiderivative<Antiderivative<FunctionType,ReturnType>,ReturnType> PrimitiveY(PrimitiveZ,StartY,EndY,WriteMessage,digits);

  return GaussLobattoAdaptiveQuadrature(PrimitiveY,StartX,EndX,_dummy_param,WriteMessage,digits);
};


/*******************************************************************************************
 * TransformFunctionInPlan                                                                 *
 * On input:                                                                               *
 *  -> Ptr_F:       the function that is transformed                                       *
 *  -> TransformX:  the transformation of the x coordinate (has the type FunctionType2D)   *
 *  -> TransformY:  the transformation of the y coordinate (has the type FunctionType2D)   *
 *  -> Jacobian:    the Jacobian of the transformation (has the type FunctionType2D)       *
 *                                                                                         *
 * On return:                                                                              *
 *  -> operator(p,q): returns the value of the transformed function to the new coordinates *
 *                    multiplied by the Jacobian of the transformation                     *
 *                    F(p,q) = Ptr_F(TransformX(p,q),TransformY(p,q)) * Jacobian(p,q)      *
 ******************************************************************************************/
template<class FunctionType, class NodeType, class ReturnType>
  class BilinearTransformFunctionInPlan{
 private:
  FunctionType Ptr_F;
  double a0,a1,a2,a3;
  double b0,b1,b2,b3;
 public:
  double BilinearTransformationX (double p, double q){
    return a0 + a1*p + a2*q + a3*p*q;
  }

  double BilinearTransformationY (double p, double q){
    return b0 + b1*p + b2*q + b3*p*q;
  }

  double Jacobian (double p, double q){
    return (a1*b2-a2*b1) + (a1*b3-a3*b1)*p + (b2*a3-a2*b3)*q;
  }

  BilinearTransformFunctionInPlan(const FunctionType Ptr_F_, const NodeType & SW, const NodeType & NW,
				  const NodeType & NE, const NodeType & SE)
    : Ptr_F(Ptr_F_){
    /* determine the coefficients of the transformation */
    a0 = SW.x();
    a1 = SE.x() - SW.x();
    a2 = NW.x() - SW.x();
    a3 = SW.x() + NE.x() - NW.x() - SE.x();
    
    b0 = SW.y();
    b1 = SE.y() - SW.y();
    b2 = NW.y() - SW.y();
    b3 = SW.y() + NE.y() - NW.y() - SE.y();
  };


    ReturnType operator() (double p, double q){
      //if (Jacobian(p,q) == 0.0) Print_("Jacobian equals zero");
      return Ptr_F(BilinearTransformationX(p,q), BilinearTransformationY(p,q)) * Jacobian(p,q);
    }

    friend ostream& operator<< (ostream& os,
				const BilinearTransformFunctionInPlan<FunctionType,NodeType,ReturnType>& Obj){
      os << endl
	 << "BilinearTransformFunctionInPlan::Coefficients\n"
	 << "a0 = " << Obj.a0 << endl
	 << "a1 = " << Obj.a1 << endl
	 << "a2 = " << Obj.a2 << endl
	 << "a3 = " << Obj.a3 << endl
	 << "b0 = " << Obj.b0 << endl
	 << "b1 = " << Obj.b1 << endl
	 << "b2 = " << Obj.b2 << endl
	 << "b3 = " << Obj.b3 << endl;

      return os;
    }

};

/******************************************************************************************
 * QuadrilateralQuadrature                                                                *
 * Integrates the function "func" over a closed domain determined by 4 Nodes              *
 * The domain is transformed to a unit square over which the integral is solved using     *
 * the AdaptiveGaussianQuadrature procedure.                                              *
 * On input:                                                                              *
 *     -> func: the function to be integrated                                             *
 *     -> SE, NE, NW, SW : the nodes that define the domain                               *
 *     -> digits: the precision of the integration                                        *
 *     -> _dummy_param : a parameter provided only to show the return type of the function * 
 ******************************************************************************************/

template<class FunctionType, class NodeType, class ReturnType>
  ReturnType QuadrilateralQuadrature(FunctionType func, NodeType & SW, NodeType & NW, NodeType & NE, 
				     NodeType & SE, int digits, ReturnType _dummy_param){

  /* Create the function that is integrated over a square */

  BilinearTransformFunctionInPlan<FunctionType,NodeType,ReturnType> TransformedFunction(func,SW,NW,NE,SE);

  /* Integrate the new function over the square defined by (0,0) , (0,1) , (1,0) and (1,1) */

  return AdaptiveGaussianQuadrature(TransformedFunction,0.0,1.0,0.0,1.0,digits,_dummy_param);
}

template<class FunctionType, class NodeType, class ReturnType>
  ReturnType GaussLobattoQuadrilateralQuadrature(FunctionType func, NodeType & SW, NodeType & NW, NodeType & NE, 
						 NodeType & SE, ReturnType _dummy_param){

  /* Create the function that is integrated over a square */

  BilinearTransformFunctionInPlan<FunctionType,NodeType,ReturnType> TransformedFunction(func,SW,NW,NE,SE);

  /* Integrate the new function over the square defined by (0,0) , (0,1) , (1,0) and (1,1) */
  return GaussLobattoQuadrature(TransformedFunction,0.0,1.0,0.0,1.0,_dummy_param);
}

template<class FunctionType, class NodeType, class ReturnType>
  ReturnType Gauss5PointQuadrilateralQuadrature(FunctionType func, NodeType & SW, NodeType & NW, NodeType & NE, 
						NodeType & SE, ReturnType _dummy_param){

  /* Create the function that is integrated over a square */

  BilinearTransformFunctionInPlan<FunctionType,NodeType,ReturnType> TransformedFunction(func,SW,NW,NE,SE);

  /* Integrate the new function over the square defined by (0,0) , (0,1) , (1,0) and (1,1) */
  return Gauss5PointQuadrature(TransformedFunction,0.0,1.0,0.0,1.0,_dummy_param);
}

/******************************************************************************************
 * Generalized polynomial function of One variable                                        *
 * is a class of functions which have the form (x-xi)^n                                   *
 ******************************************************************************************/

// Generalized polynomial function
class GeneralizedPolynomialFunctionOfOneVariable{
 private:
  double xi;
  int n;
 public:
  GeneralizedPolynomialFunctionOfOneVariable(const int & _n_, const double & _xi_)
    : xi(_xi_) , n(_n_){
  };

    void ChangePowerTo(const int & _n_){
      n = _n_;
    }

    double operator()(const double & x){
      return pow((x-xi),n);
    }
};

/******************************************************************************************
 * Generalized polynomial function of Two variables                                       *
 * is a class of functions which have the form (x-xi)^n * (y-yi)^m                        *
 ******************************************************************************************/

// Generalized polynomial function
class GeneralizedPolynomialFunctionOfTwoVariables{
 private:
  double xi, yi;
  int n,m;
 public:
  GeneralizedPolynomialFunctionOfTwoVariables(const int _n_, const int _m_, const double _xi_,
					      const double _yi_)
    : xi(_xi_), yi(_yi_), n(_n_), m(_m_){
  };

    void ChangePowersTo(const int _n_, const int _m_){
      n = _n_;
      m = _m_;
    }

    double operator()(const double x, const double y){
      return pow((x-xi),n)*pow((y-yi),m);
    }

    friend ostream& operator << (ostream& os, const GeneralizedPolynomialFunctionOfTwoVariables & Var){
      os << endl;
      os << "xi=" << Var.xi << endl
	 << "yi=" << Var.yi << endl
	 << "power n=" << Var.n << endl
	 << "power m=" << Var.m << endl;
      return os;
    }
};


/******************************************************************************************
Function for generating the geom coeff. for cartesian cell
******************************************************************************************/
double GeomCoeffCartesian(int p1, int p2, double deltaX, double deltaY, double deltaXC, double deltaYC);

// MakeReconstructionStencil()
// Set the stencil for the 1D kExact reconstruction
void MakeReconstructionStencil(const int & rings, const int & iCell, vector<int> & i_index);

// MakeReconstructionStencil()
// Set the stencil for the 2D kExact reconstruction
void MakeReconstructionStencil(const int & rings, const int & iCell, const int & jCell,
			       vector<int> & i_index, vector<int> & j_index);

// MakeReconstructionStencil()
// Enlarge the stencil for the 2D kExact reconstruction used at curved boundaries
void MakeReconstructionStencil(const int & rings, const int & iCell, const int & jCell,
			       const int NorthCurvedBnd, const int SouthCurvedBnd,
			       const int EastCurvedBnd, const int WestCurvedBnd,
			       const int &ICl, const int &ICu, const int &JCl, const int &JCu,
			       int & StencilDimension, 
			       vector<int> & i_index, vector<int> & j_index);

// ZeroLineIntegration()
double ZeroLineIntegration(const double & N1x, const double & N1y,
			   const double & N2x, const double & N2y);

// PolynomLineIntegration()
double PolynomLineIntegration(const double & N1x, const double & N1y,
			      const double & N2x, const double & N2y,
			      const double & xCC, const double & yCC,
			      const int OrderX,   const int OrderY);

template< class Node>
inline double ZeroLineIntegration(const Node& StartNode, const Node& EndNode){
  return ZeroLineIntegration(StartNode.x(), StartNode.y(), EndNode.x(), EndNode.y());
}

/**************** Function Prototypes ********************/

/*### void frenel(double, double &, double &) ###
  ################################################
  Returns the sin (*s) and cos (*c) values for the Frenel function applied to "x"*/
void frenel(double x, double &s, double &c);

/*### double qgauss10(const FunctionType1D , const double , const double )
  ########################################################################
  Returns the integral of the function "func" between a and b,
  by ten-point Gauss-Legendre integration: the function is evaluated exactly
  ten times at interior points in the range of integration.
  The integral is exact for polynomials up to order of 19. 
  Implementation based on the subroutine from Numerical Recipes*/
double qgaus10(const FunctionType1D func, const double a, const double b);


double qgaus5(const SuperFunctionType1D func, const FunctionType1D SubFunc1,
	      const FunctionType1D SubFunc2, const double a, const double b);
/*returns the integral of the function "func" between a and b,
  by five-point Gauss-Legendre integration: the function is evaluated exactly
  five times at interior points in the range of integration.
  The integral is exact for polynomials up to order of 9. 
  Implementation based on the subroutine from Numerical Recipes*/


double AdaptiveGaussianQuadrature(const SuperFunctionType1D func, const FunctionType1D SubFunc1,
				  const FunctionType1D SubFunc2, const double a,const double b,
				  int digits);
/*returns the integral of the function "func" between a and b,
  by an adaptive Gauss-Legendre integration: the integral is evaluated
  with a precision given by the required number of exact digits.
  NOTE: For some functions the required number of exact digits cannot be obtained (e.g. functions with
  singularities). Anyway, a precision of at least 7 exact digits was obtained for all tested functions,
  including functions with singularities!*/

double quad2d(const FunctionType2D func,const double a, const double b, const double c, const double d);
/*returns the integral of a user-supplied func over a two-dimensional 
  rectangular region, specified by the limits x1=a, x2=b, y1=c, y2=d.
  Integration is performed by calling qgaus10 recursively.*/

double f1(const double x);
double f2(const double y);

/*f1 and f2 are parts of the quad2d subroutine*/

double quad2dAdaptiveGaussianQuadrature(const FunctionType2D func, const double a,
					const double b, const double c, const double d,
					int digits);
/*returns the integral of a user-supplied func over a two-dimensional 
  rectangular region, specified by the limits x1, x2, y1, y2.
  Integration is performed by calling AdaptiveGaussianQuadrature recursively.*/

double AGQf1(const double x);

/*AGQf1 is part of the quad2dAdaptiveGaussianQuadrature subroutine*/

double Error1D(const double x, const FunctionType1D func1, const FunctionType1D func2);
/*returns the absolute difference between the evaluations of two different functions
  in the same point "x" of the 1D space*/

double error2D(const FunctionType2D func1, const FunctionType2D func2, const double x, const double y);
/*returns the absolute difference between the evaluations of two different functions
  in the same point "(x,y)" of the 2D space*/


#endif // _NumericalLibrary_INCLUDED

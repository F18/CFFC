/* NumericalLibrary.h: Header file for some useful numerical methods. */

#ifndef _NumericalLibrary_INCLUDED
#define _NumericalLibrary_INCLUDED

#include <complex>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <vector>

#ifdef __NUM_LIMITS__
#include <limits>
#endif // __NUM_LIMITS__

#ifndef _MATH_MACROS_INCLUDED
#include "Math.h"
#endif // _MATH_MACROS_INCLUDED

/* Include require header file */
#ifndef _REQUIRE_H
#include "include/require.h"
#endif // _REQUIRE_H_INCLUDED

using namespace std;

/***************** Type Definition *****************/

typedef double (* FunctionType1D) (const double);
typedef double (* FunctionType2D) (const double, const double);
typedef double (* FunctionType3D) (const double, const double, const double);
typedef double (* SuperFunctionType1D) (const double, const FunctionType1D, const FunctionType1D);

/***************** Exception Definition *************/

/**
 * Exception to throw when the number of iteration
 * exceeds the maximum value
 */
struct maximum_exceeded : public std::logic_error
{
  maximum_exceeded() : std::logic_error("maximum number of iteration exceeded") {};
};

/**
 * Exception to throw when the required accuracy is
 * less than the machine accuracy
 */
struct machine_accuracy_exceeded : public std::logic_error
{
  machine_accuracy_exceeded(): std::logic_error("the values are beyond the accuracy of the machine") {};
};

/**
 * Exception to throw when the relative length of the interval for integration is
 * less than the machine accuracy
 */
struct too_short_interval : public std::logic_error
{
  too_short_interval(): std::logic_error("the dimension of the interval for integration is less than the machine accuracy") {};
};

/**************** Function Prototypes ********************/

void frenel(double x, double *s, double *c);
/* returns the sin (*s) and cos (*c) values for the Frenel function applied to "x"*/

double qgaus10(const FunctionType1D func, const double a, const double b);
/*returns the integral of the function "func" between a and b,
  by ten-point Gauss-Legendre integration: the function is evaluated exactly
  ten times at interior points in the range of integration.
  The integral is exact for polynomials up to order of 19. 
  Implementation based on the subroutine from Numerical Recipes*/

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


// Error Function Template
template<typename FO1, typename FO2>
class ErrorFunc {
 private:
  FO1 fo1;
  FO2 fo2;

 public:
  // constructor
  ErrorFunc(FO1 f1, FO2 f2)
    : fo1(f1), fo2(f2){
  }

  // "function call"
  double operator() (double v) const {

    return fabs(fo1(v) - fo2(v));
  }
};

template<typename FO1, typename FO2> inline
ErrorFunc<FO1,FO2> error_func(FO1 f1, FO2 f2){
  return ErrorFunc<FO1,FO2> (f1,f2);
}

// Error function specialized for objects with member function: SolutionAtCoordinates()
template<class FunctionType, class ObjectType, class SolutionType>
class _Error_{
private:
  FunctionType Ptr_F1;
  ObjectType * Ptr_F2;

public:
  _Error_(FunctionType F1, ObjectType * F2)
    :Ptr_F1(F1), Ptr_F2(F2){ 
  };

  SolutionType operator() (double val1){
    return fabs(Ptr_F1(val1) - Ptr_F2->SolutionAtCoordinates(val1));
  }

  SolutionType operator() (double val1, double val2){
    return fabs(Ptr_F1(val1,val2) - Ptr_F2->SolutionAtCoordinates(val1,val2));
  }

  SolutionType operator() (double val1, double val2, double val3){
    return fabs(Ptr_F1(val1,val2,val3) - Ptr_F2->SolutionAtCoordinates(val1,val2,val3));
  }

};

// Square Error function specialized for objects with member function: SolutionAtCoordinates()
// returns the square of the error
template<class FunctionType, class ObjectType, class SolutionType>
class _Error_Square_{
private:
  FunctionType Ptr_F1;
  ObjectType * Ptr_F2;
  SolutionType FunctionValue;

public:
  _Error_Square_(FunctionType F1, ObjectType * F2)
    :Ptr_F1(F1), Ptr_F2(F2){ 
  };

  SolutionType operator() (double val1){
    FunctionValue = Ptr_F1(val1) - Ptr_F2->SolutionAtCoordinates(val1);
    return FunctionValue * FunctionValue;
  }

  SolutionType operator() (double val1, double val2){
    FunctionValue = Ptr_F1(val1,val2) - Ptr_F2->SolutionAtCoordinates(val1,val2);
    return FunctionValue * FunctionValue;
  }

  SolutionType operator() (double val1, double val2, double val3){
    FunctionValue = Ptr_F1(val1,val2,val3) - Ptr_F2->SolutionAtCoordinates(val1,val2,val3);
    return FunctionValue * FunctionValue;
  }
};


// qgaus5 Template
template<class FunctionType, class ReturnType>
ReturnType qgaus5(FunctionType func, const double a, const double b, ReturnType _dummy_param){
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
	sum = _dummy_param; // that's only for the compiler
	static double x[]={0.0, 0.53846931010568311, 0.90617984593866396};
	static double w[]={0.56888888888888889, 0.47862867049936647, 0.23692688505618917};

	xm=0.5*(b+a);
	xr=0.5*(b-a);
	sum=(ReturnType)(w[0]*(func(xm)));
	for (j=1;j<=2;j++) {
	  dx=xr*x[j];
	  sum = sum + (ReturnType)(w[j]*(func(xm+dx)+func(xm-dx)));
	}
	sum = (ReturnType)(sum*xr);
	return sum;
}

/***************************************************************************************
 * AdaptiveGaussianQuadrature Template                                                 *
 * returns the integral of a user-supplied function "func(x)" over a one-dimensional   * 
 * region, specified by the limits StartX, EndX.                                       *
 * Integration is performed by calling qgaus5 recursively.                             *
 **************************************************************************************/
template<class FunctionType, class ReturnType>
ReturnType AdaptiveGaussianQuadrature(FunctionType func, const double & a,const double & b,
				       int digits, ReturnType _dummy_param){


  long double EPS = 0.05*pow(10.0,1.0-digits); // accuracy: --> based on precision 
                                         // i.e exact number of digits 
  ReturnType g5;                         // the value of the integral obtained
                                         // with qgaus5
  ReturnType RERR(1);	        	 // Relative Error
  ReturnType G(0.0);			// Integral value
  long double LengthInterval = fabs(b - a); 
  long double StartPoint = a;
  long double EndPoint = b;
  long double SubDivisionPoint;
  double Division = 0.5;
#ifdef __NUM_LIMITS__
  long double EPS_Machine = numeric_limits<long double>::epsilon( );
#else
  long double EPS_Machine = 1.0e-15;
#endif // __No_NUM_LIMITS__

  try
    {
      if (0.5*LengthInterval <= EPS_Machine)
	throw too_short_interval();
      g5 = qgaus5 (func,StartPoint,EndPoint,_dummy_param);
      SubDivisionPoint = StartPoint + Division*(EndPoint - StartPoint);
      G = qgaus5(func,StartPoint,SubDivisionPoint,_dummy_param)+ qgaus5(func,SubDivisionPoint,EndPoint,
									_dummy_param);

      RERR = fabs(G - g5)/(1.0 + fabs(G));
       if (RERR < (ReturnType)EPS){
	return G;
      }
      else {
	G = ( AdaptiveGaussianQuadrature(func,StartPoint, SubDivisionPoint,digits, _dummy_param) +
	      AdaptiveGaussianQuadrature(func,SubDivisionPoint,EndPoint, digits, _dummy_param) );
	return G;
      }
    }
  
  catch (const maximum_exceeded &)
    {
      std::string msg("Maximum iterations exceeded. The precision couldn't be obtained!");
      std::cerr << msg << std::endl;
      return (ReturnType)0;
    }

  catch (const machine_accuracy_exceeded &)
    {
      std::string msg("The machine accuracy exceeded! The precision cannot be obtained!");
      std::cerr << msg << std::endl;
      return (ReturnType)0;
    }

  catch (const too_short_interval &)
    {
      std::string msg("Too short interval! The precision cannot be obtained!");
      std::cerr << msg << std::endl;
      exit(0);
    }

  catch (...)
    {
      std::cerr << "Undefined error!" << std::endl;
    }

  /* this is just for the compiler */
  return (ReturnType)0;
}


/**************************************************************************
 * FunctorX: defines a new function which has the X variable fixed to Val *
 *************************************************************************/
template<class FunctionType, class SolutionType>
class FunctorX
{
 private:
  FunctionType Ptr_F;		/* this function takes two or three parameters x, y & z */
  double Val;

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
  double Val;

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
  double Val;

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
  double ValX;
  double ValZ;

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
  double ValX;
  double ValY;

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
  double ValY;
  double ValZ;

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
  SolutionType DummyParam;
  double LeftLimit;
  double RightLimit;
  int ExactDigits;
 public:
  Antiderivative(FunctionType F, const double & StartY, const double & EndY, const int & digits)
    : Ptr_F(F), LeftLimit(StartY), RightLimit(EndY), ExactDigits(digits) {
  };

  SolutionType operator() (const double & x){
    FunctorX<FunctionType,SolutionType> Functor(Ptr_F,x);
    return AdaptiveGaussianQuadrature(Functor, LeftLimit, RightLimit, ExactDigits, DummyParam);
  }

  SolutionType operator() (const double & x, const double & y){
    FunctorXY<FunctionType,SolutionType> Functor(Ptr_F,x,y);
    return AdaptiveGaussianQuadrature(Functor, LeftLimit, RightLimit, ExactDigits, DummyParam);

  }
};

/***************************************************************************************
 * AdaptiveGaussianQuadrature Template                                                 *
 * returns the integral of a user-supplied function "func(x,y)" over a two-dimensional * 
 * rectangular region, specified by the limits StartX, EndX, StartY, EndY.             *
 * Integration is performed by calling AdaptiveGaussianQuadrature for the x and y      *
 * directions.                                                                         *
 **************************************************************************************/
template<class FunctionType, class ReturnType>
ReturnType AdaptiveGaussianQuadrature(FunctionType func, const double & StartX,const double & EndX,
				      const double & StartY, const double & EndY,  
				      int digits, ReturnType _dummy_param){

  /* Create the function to be integrated with X */
  Antiderivative<FunctionType,ReturnType> PrimitiveY(func,StartY,EndY,digits);

  return AdaptiveGaussianQuadrature(PrimitiveY,StartX,EndX,digits,_dummy_param);
}

/*******************************************************************************************
 * AdaptiveGaussianQuadrature Template                                                     *
 * returns the integral of a user-supplied function "func(x,y,z)" over a three-dimensional * 
 * rectangular region, specified by the limits StartX, EndX, StartY, EndY, StartZ, EndZ    *
 * Integration is performed by calling AdaptiveGaussianQuadrature for the x, y and z       *
 * directions.                                                                             *
 ******************************************************************************************/
template<class FunctionType, class ReturnType>
ReturnType AdaptiveGaussianQuadrature(FunctionType func, const double & StartX,const double & EndX,
				      const double & StartY, const double & EndY, const double & StartZ,
				      const double & EndZ, int digits, ReturnType _dummy_param){

  /* Create the function to be integrated with X */
  Antiderivative<FunctionType,ReturnType> PrimitiveZ(func,StartZ,EndZ,digits);

  return AdaptiveGaussianQuadrature(PrimitiveZ,StartX,EndX,StartY,EndY,digits,_dummy_param);
}


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

  double BilinearTransformationX (double p, double q){
    return a0 + a1*p + a2*q + a3*p*q;
  }

  double BilinearTransformationY (double p, double q){
    return b0 + b1*p + b2*q + b3*p*q;
  }

  double Jacobian (double p, double q){
    return (a1*b2-a2*b1) + (a1*b3-a3*b1)*p + (b2*a3-a2*b3)*q;
  }

 public:
  BilinearTransformFunctionInPlan(FunctionType Ptr_F_, NodeType SW, NodeType NW, NodeType NE, 
				  NodeType SE)
    : Ptr_F(Ptr_F_){
    /* determine the coefficients of the transformation */
    a0 = SW.x;
    a1 = SE.x - SW.x;
    a2 = NW.x - SW.x;
    a3 = SW.x + NE.x - NW.x - SE.x;
    
    b0 = SW.y;
    b1 = SE.y - SW.y;
    b2 = NW.y - SW.y;
    b3 = SW.y + NE.y - NW.y - SE.y;
  };


  ReturnType operator() (double p, double q){
    return Ptr_F(BilinearTransformationX(p,q), BilinearTransformationY(p,q)) * Jacobian(p,q);
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
ReturnType QuadrilateralQuadrature(FunctionType func, NodeType SW, NodeType NW, NodeType NE, 
				   NodeType SE, int digits, ReturnType _dummy_param){

  /* Create the function that is integrated over a square */

  BilinearTransformFunctionInPlan<FunctionType,NodeType,ReturnType> TransformedFunction(func,SW,NW,NE,SE);

  /* Integrate the new function over the square defined by (0,0) , (0,1) , (1,0) and (1,1) */

  return AdaptiveGaussianQuadrature(TransformedFunction,0,1,0,1,digits,_dummy_param);
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
    : n(_n_), xi(_xi_){
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
  GeneralizedPolynomialFunctionOfTwoVariables(const int & _n_, const int & _m_, const double & _xi_,
					      const double & _yi_)
    : n(_n_), m(_m_), xi(_xi_), yi(_yi_){
  };

  void ChangePowersTo(const int & _n_, const int & _m_){
    n = _n_;
    m = _m_;
  }

  double operator()(const double & x, const double & y){
    return pow((x-xi),n)*pow((y-yi),m);
  }
};


/******************************************************************************************
Function for generating the geom coeff. for cartesian cell
******************************************************************************************/

double GeomCoeffCartesian(int p1, int p2, double deltaX, double deltaY, double deltaXC, double deltaYC);

// MakeReconstructionStencil(int,int,vector<int>) function ******************
// Set the stencil for the 1D DD_ENO reconstruction
inline void MakeReconstructionStencil(const int & rings, const int & iCell, vector<int> & i_index){

  // Obs. The first position (i_index[0]) corresponds to (iCell)

  switch(rings){

  case 3: // three rings of cells around (iCell)
    /* Third ring */
    i_index[5]=iCell-3;
    i_index[6]=iCell+3;

  case 2: // two rings of cells around (iCell)
    /* Second ring */
    i_index[3]=iCell-2;
    i_index[4]=iCell+2;

  case 1: // one ring of cells around (iCell,jCell)
    /* First ring */
    i_index[1]=iCell-1;
    i_index[2]=iCell+1;

  case 0: 
    i_index[0]=iCell; /* cell (iCell) */

  default: // general expression
    i_index[0]=iCell; /* cell (iCell) */
    for (int i=iCell-rings, Pos=1; i<=iCell+rings; ++i){
      if(i!=iCell){
	i_index[Pos] = i;
	++Pos;
      }
    }
  }
}

//MakeReconstructionStencil(int,int,int,vector<int>,vector<int>) function
// Set the stencil for the 2D DD_ENO reconstruction
inline void MakeReconstructionStencil(const int & rings, const int & iCell, const int & jCell,
				      int *i_index, int *j_index){

  // Obs. The first position (i_index[0],j_index[0]) corresponds to (iCell,jCell)

  switch(rings){

  case 2: // two rings of cells around (iCell,jCell)

    /* Second ring */
    i_index[9] =iCell-2;  j_index[9]=jCell-2;
    i_index[10]=iCell-1; j_index[10]=jCell-2;
    i_index[11]=iCell  ; j_index[11]=jCell-2;
    i_index[12]=iCell+1; j_index[12]=jCell-2;
    i_index[13]=iCell+2; j_index[13]=jCell-2;
    i_index[14]=iCell-2; j_index[14]=jCell-1;
    i_index[15]=iCell+2; j_index[15]=jCell-1;
    i_index[16]=iCell-2; j_index[16]=jCell;
    i_index[17]=iCell+2; j_index[17]=jCell;
    i_index[18]=iCell-2; j_index[18]=jCell+1;
    i_index[19]=iCell+2; j_index[19]=jCell+1;
    i_index[20]=iCell-2; j_index[20]=jCell+2;
    i_index[21]=iCell-1; j_index[21]=jCell+2;
    i_index[22]=iCell  ; j_index[22]=jCell+2;
    i_index[23]=iCell+1; j_index[23]=jCell+2;
    i_index[24]=iCell+2; j_index[24]=jCell+2;

  case 1: // one ring of cells around (iCell,jCell)

    i_index[0]=iCell;   j_index[0]=jCell; /* cell (iCell,jCell) */
    /* First ring */
    i_index[1]=iCell-1; j_index[1]=jCell-1;
    i_index[2]=iCell;   j_index[2]=jCell-1;
    i_index[3]=iCell+1; j_index[3]=jCell-1;
    i_index[4]=iCell-1; j_index[4]=jCell;
    i_index[5]=iCell+1; j_index[5]=jCell;
    i_index[6]=iCell-1; j_index[6]=jCell+1;
    i_index[7]=iCell;   j_index[7]=jCell+1;
    i_index[8]=iCell+1; j_index[8]=jCell+1;
    break;

  default: // general expression
    i_index[0] = iCell;
    j_index[0] = jCell;
    for (int i=iCell-rings, Poz=1; i<=iCell+rings; ++i)
      for (int j=jCell-rings; j<=jCell+rings; ++j){
	if(!((i==iCell)&&(j==jCell)) ){
	  i_index[Poz] = i;
	  j_index[Poz] = j;
	  ++Poz;
	}
      }
  }//endswitch
}

#endif // _NumericalLibrary_INCLUDED

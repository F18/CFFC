/*!\file NumericalLibrary.h
  \brief Header file defining some useful numerical methods/objects. */

#ifndef _NUMERICAL_LIBRARY_INCLUDED
#define _NUMERICAL_LIBRARY_INCLUDED

/* Include required C++ libraries. */
#include <limits>
#include <iostream>
#include <cmath>

/* Using std namespace functions */
using std::numeric_limits;
using std::cout;
using std::endl;
using std::ostream;
using std::fabs;
using std::pow;
using std::sqrt;


/* Include CFFC header files */
#include "../Utilities/TypeDefinition.h"
#include "../Utilities/Utilities.h"
#include "Math.h"
#include "../System/System_Linux.h"

#define MaxAllowedFunctionEvaluation 8000000 /* maximum allowed function evaluations for the adaptive integration */

/*******************************************************************
 *                                                                 *
 *               GENERAL ERROR FUNCTIONS & WRAPPERS                *
 *                                                                 *
 ******************************************************************/


/************************************************************************//**
 * \class ErrorFunc                                                    
 * \brief General templated error function.                             
 *                                                                     
 * Computes the absolute value of the error between two specified functions.
 ****************************************************************************/
template<typename SolutionType, typename FunctionObject_1, typename FunctionObject_2 = FunctionObject_1>
  class ErrorFunc {
  
  public:
  /* Constructor with the two function*/
    ErrorFunc(FunctionObject_1 f1, FunctionObject_2 f2): fo1(f1), fo2(f2){ }

  /* Constructor with  */
  
  // "function evaluation" with one parameter
  SolutionType operator() (double param1) {
    return fabs(fo1(param1) - fo2(param1));
  }

  // "function evaluation" with two parameters
  SolutionType operator() (double param1, double param2) {
    return fabs(fo1(param1,param2) - fo2(param1,param2));
  }

  // "function evaluation" with three parameters
  SolutionType operator() (double param1, double param2, double param3) {
    return fabs(fo1(param1,param2,param3) - fo2(param1,param2,param3));
  }

  private:
  ErrorFunc();			/* make the default constructor private */
  FunctionObject_1 fo1;		/*!< pointer to the first function object */
  FunctionObject_2 fo2;		/*!< pointer to the second function object */

};

/************************************************************************//**
 * \fn ErrorFunc<SolutionType,FunctionObject_1,FunctionObject_2> 
 * error_function(FunctionObject_1 f1, FunctionObject_2 f2, SolutionType dummy)   
 * \brief Generate a function that computes the absolute error between f1 and f2
 *
 * \param f1 the first function object 
 * \param f2 the second function object
 * \param dummy used only to provide the solution type
 ****************************************************************************/
template<typename FunctionObject_1, typename FunctionObject_2, typename SolutionType> 
inline
ErrorFunc<SolutionType, FunctionObject_1, FunctionObject_2> // returned type
error_function(FunctionObject_1 f1, FunctionObject_2 f2, SolutionType & dummy){
  return ErrorFunc<SolutionType, FunctionObject_1, FunctionObject_2> (f1,f2);
}

/************************************************************************//**
 * \class SquareErrorFunc                                                    
 * \brief General templated error function.                             
 *                                                                     
 * Computes the squared value of the error between two specified functions.
 ****************************************************************************/
template<typename SolutionType, typename FunctionObject_1, typename FunctionObject_2 = FunctionObject_1>
class SquareErrorFunc {
  
public:
  /* Constructor with the two function*/
  SquareErrorFunc(FunctionObject_1 f1, FunctionObject_2 f2): fo1(f1), fo2(f2){ }

  // "function evaluation" with one parameter
  SolutionType operator() (double param1) {
    Temp = fo1(param1) - fo2(param1);
    return Temp * Temp;
  }

  // "function evaluation" with two parameters
  SolutionType operator() (double param1, double param2) {
    Temp = fo1(param1,param2) - fo2(param1,param2);
    return Temp * Temp;
  }

  // "function evaluation" with three parameters
  SolutionType operator() (double param1, double param2, double param3) {
    Temp = fo1(param1,param2,param3) - fo2(param1,param2,param3);
    return Temp*Temp;
  }

  private:
  SquareErrorFunc();			/* make the default constructor private */
  FunctionObject_1 fo1;		/*!< pointer to the first function object */
  FunctionObject_2 fo2;		/*!< pointer to the second function object */

  SolutionType Temp;		// temporary object
};

/************************************************************************//**
 * \fn SquareErrorFunc<FunctionObject_1,FunctionObject_2,SoutionType> 
 * square_error_function(FunctionObject_1 f1, FunctionObject_2 f2, SolutionType dummy)   
 * \brief Generate a function that computes the absolute error between f1 and f2
 *
 * \param f1 the first function object 
 * \param f2 the second function object
 * \param dummy used only to provide the solution type
 ****************************************************************************/
template<typename FunctionObject_1, typename FunctionObject_2, typename SolutionType> 
inline
SquareErrorFunc<SolutionType,FunctionObject_1,FunctionObject_2> 
square_error_function(FunctionObject_1 f1, FunctionObject_2 f2, SolutionType & dummy){
  return SquareErrorFunc<SolutionType,FunctionObject_1,FunctionObject_2> (f1,f2);
}

/************************************************************************//**
 * \class _Mapped_Function_Wrapper_
 * \brief This wrapper maps the original function into a new domain.
 *
 * The new domain is obtained with the linear transformation ConvertDomain()
 ****************************************************************************/
template<class FunctionObject, class SolutionType>
  class _Mapped_Function_Wrapper_{
  
 public:
  
  /* constructor */
  _Mapped_Function_Wrapper_(FunctionObject func,
			    double _DomainLimitMin_, double _DomainLimitMax_,
			    double _FunctLimitMin_, double _FunctLimitMax_):
    Func(func), DomainLimitMin(_DomainLimitMin_), DomainLimitMax(_DomainLimitMax_),
    FunctLimitMin(_FunctLimitMin_), FunctLimitMax(_FunctLimitMax_) { }
  
  // "member function evaluation" with one parameter
  SolutionType operator() (double val1){
    NewCoordinate = ConvertDomain(DomainLimitMin,DomainLimitMax,FunctLimitMin,FunctLimitMax,val1);
    return Func(NewCoordinate);
  }
  
 private:
  _Mapped_Function_Wrapper_();	/* make default constructor private */

  FunctionObject Func;		        /*!< pointer to the function object */
  
  double DomainLimitMin, DomainLimitMax; /*!< specify the min and max limits of the evaluation domain */
  double FunctLimitMin, FunctLimitMax;   /*!< specify the min and max limits of the function definition domain */
  double NewCoordinate;		/* the variable in the mapped domain */
};


/************************************************************************//**
 * \fn _Mapped_Function_Wrapper_<FunctionObject,SolutionType> 
 mapped_function(FunctionObject func, SolutionType dummy,
 double DomainLimitMin, double DomainLimitMax,
 double FunctLimitMin, double FunctLimitMax)
 * \brief Generate a function that maps the original function into a new domain.
 *
 * \param func the function object
 * \param dummy used only to provide the solution type
 * \param DomainLimitMin set the min limit of the evaluation domain
 * \param DomainLimitMax set the max limit of the evaluation domain
 * \param FunctLimitMin set the min limit of the function definition domain
 * \param FunctLimitMax set the max limit of the function definition domain
 ****************************************************************************/
template<typename FunctionObject, typename SolutionType> 
inline _Mapped_Function_Wrapper_<FunctionObject,SolutionType> 
mapped_function(FunctionObject func, SolutionType dummy,
		double DomainLimitMin, double DomainLimitMax,
		double FunctLimitMin, double FunctLimitMax){
  return _Mapped_Function_Wrapper_<FunctionObject,SolutionType> (func,
								 DomainLimitMin,DomainLimitMax,
								 FunctLimitMin,FunctLimitMax);
}


/************************************************************************//**
 * \class _Member_Function_Wrapper_
 * \brief Adaptor for making a class member function look like an ordinary function
 *
 * This wrapper is useful for passing member functions to subroutines that 
 * take ordinary functions (e.g. numerical integration subroutines)
 ****************************************************************************/
template<class ObjectType, class Member_Pointer, class SolutionType>
  class _Member_Function_Wrapper_{
  
 public:
  
  /* constructor */
  _Member_Function_Wrapper_(ObjectType *object,
			    Member_Pointer mem_func): Obj(object), Ptr(mem_func){ }
    
    // "member function evaluation" with three parameters
    SolutionType operator() (double val1, double val2, double val3){
      return (Obj->*Ptr)(val1,val2,val3);
    }
    
    // "member function evaluation" with two parameters
    SolutionType operator() (double val1, double val2){
      return (Obj->*Ptr)(val1,val2);
    }
    
    // "member function evaluation" with one parameter
    SolutionType operator() (double val1){
      return (Obj->*Ptr)(val1);
    }
    
  private:
    _Member_Function_Wrapper_();	/* make default constructor private */
    ObjectType *Obj;		/*!< pointer to the object */
    Member_Pointer Ptr;		/*!< pointer to the class member function */
  };


/************************************************************************//**
 * \fn _Member_Function_Wrapper_<ObjectType,Member_Pointer,SolutionType> 
 wrapped_member_function(ObjectType *object, Member_Pointer mem_func, SolutionType dummy)
 * \brief Adaptor for making a class member function look like an ordinary function.
 *
 * \param object the object used to access the member function
 * \param mem_func the member function object
 * \param dummy used only to provide the solution type
 ****************************************************************************/
template<typename ObjectType, typename Member_Pointer, typename SolutionType>
inline _Member_Function_Wrapper_<ObjectType,Member_Pointer,SolutionType> 
wrapped_member_function(ObjectType *object, Member_Pointer mem_func, SolutionType dummy){
  return _Member_Function_Wrapper_<ObjectType,Member_Pointer,SolutionType> (object,mem_func);
}

/************************************************************************//**
 * \class _Member_Function_Wrapper_One_Parameter_
 * \brief Adaptor for making a class member function that takes one (unsigned) parameter
 * look like an ordinary function
 *
 * This wrapper is useful for passing member functions to subroutines that 
 * take ordinary functions (e.g. numerical integration subroutines)
 ****************************************************************************/
template<class ObjectType, class Member_Pointer, class SolutionType>
class _Member_Function_Wrapper_One_Parameter_{
  
public:
  
  /* constructor */
  _Member_Function_Wrapper_One_Parameter_(ObjectType *object,
					  Member_Pointer mem_func,
					  const unsigned param): Obj(object), Ptr(mem_func), parameter(param){ }
  
  // "member function evaluation" with three parameters
  SolutionType operator() (double val1, double val2, double val3){
    return (Obj->*Ptr)(val1,val2,val3,parameter);
  }
    
  // "member function evaluation" with two parameters
  SolutionType operator() (double val1, double val2){
    return (Obj->*Ptr)(val1,val2,parameter);
  }
    
  // "member function evaluation" with one parameter
  SolutionType operator() (double val1){
    return (Obj->*Ptr)(val1,parameter);
  }
    
private:
  _Member_Function_Wrapper_One_Parameter_();	/* make default constructor private */
  ObjectType *Obj;		/*!< pointer to the object */
  Member_Pointer Ptr;		/*!< pointer to the class member function */
  unsigned parameter;         /*!< variable to store the function parameter value */
};


/************************************************************************//**
 * \fn _Member_Function_Wrapper_One_Parameter<ObjectType,Member_Pointer,SolutionType> 
 wrapped_member_function_one_parameter(ObjectType *object, Member_Pointer mem_func, SolutionType dummy)
 * \brief Adaptor for making a class member function that takes one (unsigned) parameter
 * look like an ordinary function.
 *
 * \param object the object used to access the member function
 * \param mem_func the member function object
 * \param param the value of the parameter used to evaluate the function
 * \param dummy used only to provide the solution type
 ****************************************************************************/
template<typename ObjectType, typename Member_Pointer, typename SolutionType>
inline _Member_Function_Wrapper_One_Parameter_<ObjectType,Member_Pointer,SolutionType> 
wrapped_member_function_one_parameter(ObjectType *object, Member_Pointer mem_func,
				      unsigned parameter, SolutionType dummy){
  return _Member_Function_Wrapper_One_Parameter_<ObjectType,Member_Pointer,SolutionType> (object,mem_func,parameter);
}

/************************************************************************//**
 * \class _Soln_Block_Member_Function_Wrapper_
 * \brief Adaptor for member functions of a solution block
 *
 * This adaptor is design to make a member function of a structured solution 
 * block that takes as arguments iCell,jCell, X_Coord and Y_Coord 
 * look like a function only of 'x' and 'y' Cartesian coordinates. \n
 * This wrapper is useful for passing member functions to subroutines that 
 * take ordinary functions (e.g. numerical integration subroutines)
 ****************************************************************************/
template<class ObjectType, class Member_Pointer, class SolutionType>
class _Soln_Block_Member_Function_Wrapper_{
  
public:
  
  /* Constructors */
  //  == 1D ==
  _Soln_Block_Member_Function_Wrapper_(ObjectType *object,
				       Member_Pointer mem_func,
				       const int & _iCell_): Obj(object), Ptr(mem_func),
							     iCell(_iCell_), jCell(0), kCell(0){ };

  // == 2D ==
  _Soln_Block_Member_Function_Wrapper_(ObjectType *object,
				       Member_Pointer mem_func,
				       const int & _iCell_,
				       const int & _jCell_): Obj(object), Ptr(mem_func),
							     iCell(_iCell_), jCell(_jCell_), kCell(0){ };

  // == 3D ==
  _Soln_Block_Member_Function_Wrapper_(ObjectType *object,
				       Member_Pointer mem_func,
				       const int & _iCell_,
				       const int & _jCell_,
				       const int & _kCell_): Obj(object), Ptr(mem_func),
							     iCell(_iCell_), jCell(_jCell_), kCell(_kCell_){ };
    
  // "member function evaluation" with three parameters (x,y,z)
  SolutionType operator() (double val1, double val2, double val3){
    return (Obj->*Ptr)(iCell,jCell,kCell,
		       val1,val2,val3);
  }
    
  // "member function evaluation" with two parameters (x,y)
  SolutionType operator() (double val1, double val2){
    return (Obj->*Ptr)(iCell,jCell,
		       val1,val2);
  }
    
  // "member function evaluation" with one parameter (x)
  SolutionType operator() (double val1){
    return (Obj->*Ptr)(iCell,
		       val1);
  }
    
private:
  /*! Private default constructor*/
  _Soln_Block_Member_Function_Wrapper_();

  // Local variables
  ObjectType *Obj;		/*!< pointer to the object */
  Member_Pointer Ptr;		/*!< pointer to the class member function */
  int iCell, jCell, kCell;	//!< the cell indexes in the structured solution block
};


/************************************************************************//**
 * \fn _Soln_Block_Member_Function_Wrapper_<ObjectType,Member_Pointer,SolutionType> 
 wrapped_member_function_one_parameter(ObjectType *object, Member_Pointer mem_func,
                                       int iCell, int jCell,
                                       SolutionType dummy)
 * \brief Adaptor for member functions of a 2D solution block
 *
 * \param object the object used to access the member function
 * \param mem_func the member function object
 * \param iCell,jCell the cell indexes in the structured solution block
 * \param dummy used only to provide the solution type
 ****************************************************************************/
template<typename ObjectType, typename Member_Pointer, typename SolutionType>
inline _Soln_Block_Member_Function_Wrapper_<ObjectType,Member_Pointer,SolutionType> 
wrapped_soln_block_member_function(ObjectType *object, Member_Pointer mem_func,
				   const int & iCell, const int & jCell,
				   SolutionType dummy){
  return ( _Soln_Block_Member_Function_Wrapper_<ObjectType,Member_Pointer,SolutionType> (object,mem_func,iCell,jCell) );
}


/************************************************************************//**
 * \class _Soln_Block_Member_Function_Wrapper_One_Parameter_
 * \brief Adaptor for member functions of a solution block
 *
 * This adaptor is design to make a member function of a structured solution 
 * block that takes as arguments iCell,jCell, PositionVector and  one (unsigned)
 * parameter look like a function only of 'x' and 'y' Cartesian coordinates. \n
 * This wrapper is useful for passing member functions to subroutines that 
 * take ordinary functions (e.g. numerical integration subroutines)
 ****************************************************************************/
template<class ObjectType, class Member_Pointer, class SolutionType, class PositionVectorType>
class _Soln_Block_Member_Function_Wrapper_One_Parameter_{
  
public:
  
  /* Constructors */
  //  == 1D ==
  _Soln_Block_Member_Function_Wrapper_One_Parameter_(ObjectType *object,
						     Member_Pointer mem_func,
						     const int & _iCell_,
						     const unsigned & param): Obj(object), Ptr(mem_func),
									      iCell(_iCell_), jCell(0), kCell(0),
									      parameter(param){ };

  // == 2D ==
  _Soln_Block_Member_Function_Wrapper_One_Parameter_(ObjectType *object,
						     Member_Pointer mem_func,
						     const int & _iCell_,
						     const int & _jCell_,
						     const unsigned & param): Obj(object), Ptr(mem_func),
									      iCell(_iCell_), jCell(_jCell_), kCell(0),
									      parameter(param){ };

  // == 3D ==
  _Soln_Block_Member_Function_Wrapper_One_Parameter_(ObjectType *object,
						     Member_Pointer mem_func,
						     const int & _iCell_,
						     const int & _jCell_,
						     const int & _kCell_,
						     const unsigned & param): Obj(object), Ptr(mem_func),
									      iCell(_iCell_), jCell(_jCell_), kCell(_kCell_),
									      parameter(param){ };
    
  // "member function evaluation" with three parameters (x,y,z)
  SolutionType operator() (double val1, double val2, double val3){
    return (Obj->*Ptr)(iCell,jCell,kCell,
		       PositionVectorType(val1,val2,val3),
		       parameter);
  }
    
  // "member function evaluation" with two parameters (x,y)
  SolutionType operator() (double val1, double val2){
    return (Obj->*Ptr)(iCell,jCell,
		       PositionVectorType(val1,val2),
		       parameter);
  }
    
  // "member function evaluation" with one parameter (x)
  SolutionType operator() (double val1){
    return (Obj->*Ptr)(iCell,
		       PositionVectorType(val1),
		       parameter);
  }
    
private:
  /*! Private default constructor*/
  _Soln_Block_Member_Function_Wrapper_One_Parameter_();

  // Local variables
  ObjectType *Obj;		/*!< pointer to the object */
  Member_Pointer Ptr;		/*!< pointer to the class member function */
  int iCell, jCell, kCell;	//!< the cell indexes in the structured solution block
  unsigned parameter;           /*!< variable to store the function parameter value */
};


/************************************************************************//**
 * \fn _Soln_Block_Member_Function_Wrapper_One_Parameter<ObjectType,Member_Pointer,SolutionType,PositionVectorType> 
 wrapped_member_function_one_parameter(ObjectType *object, Member_Pointer mem_func, PositionVectorType _dummy_Pos,
                                       int iCell, int jCell,
                                       SolutionType dummy)
 * \brief Adaptor for member functions of a 2D solution block
 *
 * \param object the object used to access the member function
 * \param mem_func the member function object
 * \param iCell,jCell the cell indexes in the structured solution block
 * \param _dummy_Pos  used only to provide the type of the position vector (e.g. Vector2D, Node2D etc.)
 * \param param the value of the parameter used to evaluate the function
 * \param dummy used only to provide the solution type
 ****************************************************************************/
template<typename ObjectType, typename Member_Pointer, typename SolutionType, typename PositionVectorType>
inline _Soln_Block_Member_Function_Wrapper_One_Parameter_<ObjectType,Member_Pointer,SolutionType,PositionVectorType> 
wrapped_member_function_one_parameter(ObjectType *object, Member_Pointer mem_func,
				      const PositionVectorType & ,
				      const int & iCell, const int & jCell,
				      unsigned parameter, SolutionType dummy){
  return ( _Soln_Block_Member_Function_Wrapper_One_Parameter_<ObjectType,Member_Pointer,
	   SolutionType,PositionVectorType> (object,mem_func,iCell,jCell,parameter) );
}


/*******************************************************************
 *                                                                 *
 *                      MATHEMATICAL FUNCTIONS                     *
 *                                                                 *
 ******************************************************************/

/**
 * \fn T sign(T X)
 * \brief Return signum of X
 ****************************************************************************/
template <class T>
T sign(T X){
  if( X != T(0.0)){
    return X/fabs(X);
  } else {
    return T(0.0);
  }
}

/**
 * \fn void frenel(double x, double & s, double & c)
 * \brief Return the sine and cosine of the Frenel function 
 * 
 * \param x parameter used to compute the Frenel function
 * \param s the sine value is returned here
 * \param c the cosine value is returned here
 **/
void frenel(double x, double &s, double &c);


/*******************************************************************
 *                                                                 *
 *        QUADRATURE SUBROUTINES FOR ONE VARIABLE FUNCTIONS        *
 *                                                                 *
 ******************************************************************/

/**
 * \fn ReturnType qgauss5(FunctionType func, double a, double b, ReturnType dummy)                             
 * \brief Compute the defined integral of a function
 *                                                                     
 * Evaluates the integral of the function "func" between "a" and "b" using 
 * the five-point Gauss-Legendre integration: <br>
 * i.e. the function is evaluated exactly five times at interior points in the range of integration. <br>
 * The integral is exact for polynomials up to order of 9. <br>
 * Implementation based on the subroutine from Numerical Recipes. <br>
 * <br>
 * If a function that takes one double and returns a solution state is provided, this subroutine will <br>
 * return the appropiate integration of all the parameters in the specified domain. <br>
 * REQUIREMENTS: <br>
 * The state class defined by the ReturnType should provide the following operators: <br>
 * "operator +"  between 2 states <br>
 * "operator *"  between a state and a scalar <br>
 * "operator /"  between two states <br>
 * "operator <"  between two states <br>
 * "a constructor which takes a double" -> this will initialize all the parameters to the taken value.
 * <br>
 * \param func the name of the function used to evaluate the integral
 * \param a the lower integration limit 
 * \param b the upper integration limit 
 * \param dummy used only to provide the ReturnType
 ****************************************************************************/
template<class FunctionType, class ReturnType>
inline ReturnType qgauss5(FunctionType func, double a, double b, const ReturnType & dummy){

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

/**
 * \fn ReturnType adaptlobstp(FunctionType func, double a, double b, const ReturnType fa, const ReturnType fb, 
 *  const ReturnType is, int &FunctionEvaluations, int & WriteMessage)
 * \brief Recursive function used by GaussLobattoAdaptiveQuadrature
 *                                                                     
 * Q = adaptlobstp('func',a,b,fa,fb,is,FunctionEvaluations,WriteMessage) tries to approximate 
 * the integral of 'F(X)' from 'a' to 'b' to an appropriate relative error. <br>
 * The argument 'func' is a string containing the name of 'F'. The remaining arguments are <br> 
 * generated by GaussLobattoAdaptiveQuadrature or by recursion. <br>
 * See also GaussLobattoAdaptiveQuadrature. <br>
 * This algorithm follows the work of Walter Gautschi, and is similar to the implementation in Matlab
 *
 * \param func the name of the function used to evaluate the integral
 * \param FunctionEvaluations counts the number of function evaluation in order to avoid more evaluations than max. allowed
 **************************************************************************************************************************/
// Lucian Ivan, 31/09/2005
template <class FunctionType, class ReturnType>
inline ReturnType adaptlobstp(FunctionType func, double a, double b, const ReturnType & fa, const ReturnType & fb, 
			      const ReturnType & is, int &FunctionEvaluations, int & WriteMessage)
  throw(TooShortInterval,MaximumIterationsExceeded){

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

  FunctionEvaluations += 5;   /* This subroutine will add 5 function evaluations */

  i2=(h/6.0)*(fa+fb+5.0*(y2+y4));
  i1=(h/1470.0)*(77.0*(fa+fb)+432.0*(y1+y5)+625.0*(y2+y4)+672.0*y3);

  if ( (is+(i1-i2)==is) || (mll<=a) || (b<=mrr) || (FunctionEvaluations>MaxAllowedFunctionEvaluation) || (fabs(i1-i2)<=eps) ){
    if ( ( (m <= a) || (b<=m) ) && (WriteMessage == 0)){
      cerr << "\nWarning Integration Subroutine: Interval contains no more machine number.\n"
       	   << "Required tolerance may not be met.\n";
      WriteMessage = 1;
    }
    if ((FunctionEvaluations>MaxAllowedFunctionEvaluation) && (WriteMessage == 0)){
      cerr << "\nWarning Integration Subroutine: Maximum function count exceeded (" << MaxAllowedFunctionEvaluation
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

/**
 * \fn ReturnType GaussLobattoAdaptiveQuadrature(FunctionType func, double StartPoint, double EndPoint,
 * ReturnType dummy, int &WriteMessage, int digits)
 * \brief Numerically evaluate integrals using adaptive Lobatto rule.
 *                                                                     
 * Be default, this subroutine approximates the integral of F(X) from StartPoint to EndPoint to machine precision. <br>
 * The function 'func' can return a state of multiple values, <br>
 * and thus, the final result has the integrals of each state entry. <br>
 * If the number of exact digits is specified, the subroutine integrates to <br> 
 * a precision specified by the number of digits. <br> <br>
 * REQUIREMENTS: <br>
 * The state class should provide the following operators: <br>
 * "operator +"  between two states <br>
 * "operator *"  between a state and a scalar <br>
 * "operator /"  between two states <br>
 * "operator <"  between two states <br>
 * "operator <="  between two states <br>
 * "operator >"  between two states <br>
 * "operator >="  between two states <br>
 * "a constructor which takes a double" -> this will initialize all the parameters to the value that is passed <br> <br>
 * The implementation follows Walter Gautschi's work. <br>
 * see http://www.inf.ethz.ch/personal/gander/adaptlob.m <br>
 * Reference: Gander, Computermathematik, Birkhaeuser, 1992.  <br>
 *
 * \param func the name of the function used to evaluate the integral
 * \param dummy used only to provide the return type
 **************************************************************************************************************************/
//  Lucian Ivan, 02/11/2005 
template <class FunctionType, class ReturnType>
inline ReturnType GaussLobattoAdaptiveQuadrature(FunctionType func, double StartPoint, double EndPoint,
						 const ReturnType & dummy, int &WriteMessage,
						 int digits = numeric_limits<double>::digits10)
  throw (TooShortInterval) {
  
  ReturnType tol(0.5*std::pow(10.0,-digits));
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
    throw TooShortInterval("GaussLobattoAdaptiveQuadrature Error: The integration interval has ZERO length!");
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

/**
 * \fn ReturnType GaussLobattoQuadrature(FunctionType func, double StartPoint, double EndPoint,
 * const ReturnType & dummy)
 * \brief Numerically evaluate integrals using a 13 points Lobatto rule which has the degree 18.
 *                                                                     
 * Be default, this subroutine approximates the integral of F(X) from StartPoint to EndPoint. <br>
 * The function 'func' can return a state of multiple values, <br>
 * and thus, the final result has the integrals of each state entry. <br>
 * For polynomial functions of degree up to 18 the result is exact. <br>
 * The implementation follows Walter Gautschi's work. <br>
 * see http://www.inf.ethz.ch/personal/gander/adaptlob.m <br>
 * Reference: Gander, Computermathematik, Birkhaeuser, 1992.  <br>
 *
 * \param func the name of the function used to evaluate the integral
 * \param dummy used only to provide the return type
 **************************************************************************************************************************/
// Lucian Ivan, 07/11/2005 
template <class FunctionType, class ReturnType>
inline ReturnType GaussLobattoQuadrature(FunctionType func, double StartPoint, double EndPoint,
					 const ReturnType & dummy)
  throw (TooShortInterval){
  
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


/**
 * \fn ReturnType AdaptiveGaussianQuadrature(FunctionType func, double StartPoint, double EndPoint,
 * const ReturnType & dummy, int digits)
 * \brief Numerically evaluate integrals of ONE-variable functions using adaptive Lobatto rule.
 *  
 * see also GaussLobattoAdaptiveQuadrature()                                                                  
 * \param func the name of the function used to evaluate the integral
 * \param StartPoint the lower integration limit 
 * \param EndPoint the upper integration limit 
 * \param dummy used only to provide the return type
 * \param digits number of exact digits (there is a default value already provided!)
 **************************************************************************************************************************/
template <class FunctionType, class ReturnType>
inline ReturnType AdaptiveGaussianQuadrature(FunctionType func, double StartPoint, double EndPoint,
					     const ReturnType & dummy, int digits = numeric_limits<double>::digits10){
  
  int WriteMessage = 0; 	/* this flag makes sure that the error message is printed only once */
  return GaussLobattoAdaptiveQuadrature(func,StartPoint,EndPoint,dummy,WriteMessage,digits);
}


/*********************************************************************
 *                                                                   *
 *  QUADRATURE SUBROUTINES AND HELPERS FOR MULTI VARIABLE FUNCTIONS  *
 *                                                                   *
 ********************************************************************/

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


/* customized antiderivative for the GaussLobatto integration (non adaptive!!!) */
template<class FunctionType, class SolutionType>
  class AntiderivativeGL 		
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


/* customized antiderivative for the Gauss integration with 5 points (non adaptive!!!) */
template<class FunctionType, class SolutionType>
  class AntiderivativeG5 	     
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
      return qgauss5(Functor, LeftLimit, RightLimit, DummyParam);
    }
  
    /* Implementation for 3D */
    SolutionType operator() (const double & x, const double & y){
      FunctorXY<FunctionType,SolutionType> Functor(Ptr_F,x,y);
      return qgauss5(Functor, LeftLimit, RightLimit, DummyParam);
    }
};

/**
 * \fn ReturnType AdaptiveGaussianQuadrature(FunctionType func, double StartX, double EndX,
 * double StartY, double EndY, int digits, const ReturnType & dummy)
 * \brief Numerically evaluate integrals of TWO-variable functions over a 2D rectangular region.
 *
 * This integration uses an adaptive Lobatto rule in each of the two directions, X and Y.
 *  
 * see also GaussLobattoAdaptiveQuadrature()                                                                  
 * \param func the name of the function used to evaluate the integral
 * \param StartPointX the lower integration limit in X-direction 
 * \param EndPointX the upper integration limit in X-direction
 * \param StartPointY the lower integration limit in Y-direction 
 * \param EndPointY the upper integration limit in Y-direction
 * \param dummy used only to provide the return type
 * \param digits number of exact digits (there is a default value already provided!)
 **************************************************************************************************************************/
template<class FunctionType, class ReturnType>
inline ReturnType AdaptiveGaussianQuadrature(FunctionType func, double StartX, double EndX,
					     double StartY, double EndY,  
					     int digits, const ReturnType & dummy){
  
  int WriteMessage = 0; 	/* this flag makes sure that the error message is written only once */
  
  /* Create the function inner integral  */
  Antiderivative<FunctionType,ReturnType> PrimitiveY(func,StartY,EndY,WriteMessage,digits);
  
  /* Integrate the inner integral with respect to X */
  return GaussLobattoAdaptiveQuadrature(PrimitiveY,StartX,EndX,dummy,WriteMessage,digits);
};

/**
 * \fn ReturnType GaussLobattoQuadrature(FunctionType func, double StartX, double EndX,
 * double StartY, double EndY, const ReturnType & dummy)
 * \brief Numerically evaluate integrals of TWO-variable functions over a 2D rectangular region.
 *
 * This integration uses a Lobatto rule (non-adaptive) in each of the two directions, X and Y.
 *  
 * see also GaussLobattoAdaptiveQuadrature()                                                                  
 * \param func the name of the function used to evaluate the integral
 * \param StartPointX the lower integration limit in X-direction 
 * \param EndPointX the upper integration limit in X-direction
 * \param StartPointY the lower integration limit in Y-direction 
 * \param EndPointY the upper integration limit in Y-direction
 * \param dummy used only to provide the return type
 **************************************************************************************************************************/
template <class FunctionType, class ReturnType>
inline ReturnType GaussLobattoQuadrature(FunctionType func, double StartX, double EndX,
					 double StartY, double EndY, 
					 const ReturnType & dummy){
  
  /* Create the inner integral */
  AntiderivativeGL<FunctionType,ReturnType> PrimitiveY(func,StartY,EndY);
  
  /* Integrate the inner integral with respect to X */
  return GaussLobattoQuadrature(PrimitiveY,StartX,EndX,dummy);
}

/**
 * \fn ReturnType Gauss5PointQuadrature(FunctionType func, double StartX, double EndX,
 * double StartY, double EndY, const ReturnType & dummy)
 * \brief Numerically evaluate integrals of TWO-variable functions over a rectangle (2D).
 *
 * This integration uses the five-point Gauss-Legendre rule (non-adaptive) in each of the two directions, X and Y. <br>
 * see also qgauss5()                                                                  
 * \param func the name of the function used to evaluate the integral
 * \param StartPointX the lower integration limit in X-direction 
 * \param EndPointX the upper integration limit in X-direction
 * \param StartPointY the lower integration limit in Y-direction 
 * \param EndPointY the upper integration limit in Y-direction
 * \param dummy used only to provide the return type
 **************************************************************************************************************************/
template <class FunctionType, class ReturnType>
inline ReturnType Gauss5PointQuadrature(FunctionType func, double StartX, double EndX,
					double StartY, double EndY, 
					const ReturnType & dummy){
  
  /* Create the inner integral */
  AntiderivativeG5<FunctionType,ReturnType> PrimitiveY(func,StartY,EndY);
  
  /* Integrate the inner integral with respect to X */
  return qgauss5(PrimitiveY,StartX,EndX,dummy);
}


/**
 * \fn ReturnType AdaptiveGaussianQuadrature(FunctionType func, double StartX, double EndX,
 * double StartY, double EndY, int digits, const ReturnType & dummy)
 * \brief Numerically evaluate integrals of THREE-variable functions over a cuboid (3D).
 *
 * This integration uses an adaptive Lobatto rule in each of the three directions, X, Y and Z.
 *  
 * see also GaussLobattoAdaptiveQuadrature()                                                                  
 * \param func the name of the function used to evaluate the integral
 * \param StartPointX the lower integration limit in X-direction 
 * \param EndPointX the upper integration limit in X-direction
 * \param StartPointY the lower integration limit in Y-direction 
 * \param EndPointY the upper integration limit in Y-direction
 * \param StartPointZ the lower integration limit in Z-direction 
 * \param EndPointZ the upper integration limit in Z-direction
 * \param dummy used only to provide the return type
 * \param digits number of exact digits (there is a default value already provided!)
 **************************************************************************************************************************/
template<class FunctionType, class ReturnType>
inline ReturnType AdaptiveGaussianQuadrature(FunctionType func, const double StartX,const double EndX,
					     const double StartY, const double EndY, const double StartZ,
					     const double EndZ, int digits, const ReturnType & dummy){
  
  int WriteMessage = 0; 	/* this flag makes sure that the error message is written only once */
  
  /* Create the function to be integrated with Y */
  Antiderivative<FunctionType,ReturnType> PrimitiveZ(func,StartZ,EndZ,WriteMessage,digits);
  
  /* Create the function to be integrated with X */
  Antiderivative<Antiderivative<FunctionType,ReturnType>,ReturnType> PrimitiveY(PrimitiveZ,StartY,EndY,WriteMessage,digits);
  
  return GaussLobattoAdaptiveQuadrature(PrimitiveY,StartX,EndX,dummy,WriteMessage,digits);
};


/**
 * \class BilinearTransformFunctionInPlan
 * \brief Map a function defined on a 2D quadrilateral to a rectangle using a bilinear transformation
 *
 *
 * On return:  <br>
 *  -> operator(p,q): returns the value of the transformed function to the new coordinates  <br>
 *                    multiplied by the Jacobian of the transformation                      <br>
 *                    F(p,q) = Ptr_F(TransformX(p,q),TransformY(p,q)) * Jacobian(p,q)       <br>
 ********************************************************************************************************/

template<class FunctionType, class NodeType, class ReturnType>
class BilinearTransformFunctionInPlan{
private:
  FunctionType Ptr_F;		//!< pointer to the input function 
  double a0,a1,a2,a3;		//!< coefficients of the transformation 
  double b0,b1,b2,b3;           //!< coefficients of the transformation 
public:
  double BilinearTransformationX (double p, double q){ //!< the transformation of the X-coordinate
    return a0 + a1*p + a2*q + a3*p*q;
  }

  double BilinearTransformationY (double p, double q){ //!< the transformation of the Y-coordinate
    return b0 + b1*p + b2*q + b3*p*q;
  }

  double Jacobian (double p, double q){	//!< the Jacobian of the transformation
    return (a1*b2-a2*b1) + (a1*b3-a3*b1)*p + (b2*a3-a2*b3)*q;
  }

  // Constructor (Input: the definition nodes of the quadrilateral domain)
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


/************************************************************************//**
 * \fn BilinearTransformFunctionInPlan<FunctionType, NodeType, ReturnType>
 * planar_bilinear_function_transformation(FunctionType func,
 * const NodeType & SW, const NodeType & NW, const NodeType & NE, const NodeType & SE,
 * ReturnType & dummy)
 * \brief Generate a BilinearTransformFunctionInPlan object 
 *
 * \param func the name of the function to be mapped
 * \param SW the south-west node of the quadrilateral
 * \param NW the north-west node of the quadrilateral
 * \param NE the north-east node of the quadrilateral
 * \param SE the south-east node of the quadrilateral
 * \param dummy used only to provide the return type
 ****************************************************************************/
template<class FunctionType, class NodeType, class ReturnType>
inline 
BilinearTransformFunctionInPlan<FunctionType, NodeType, ReturnType> // returned type
planar_bilinear_function_transformation(FunctionType func,
					const NodeType & SW, const NodeType & NW,
					const NodeType & NE, const NodeType & SE,
					ReturnType & dummy){
  return BilinearTransformFunctionInPlan<FunctionType, NodeType, ReturnType> (func,SW,NW,NE,SE);
}

/**
 * \fn  ReturnType QuadrilateralQuadrature(FunctionType func, const NodeType & SW, const NodeType & NW, const NodeType & NE, 
 *  const NodeType & SE, int digits, const ReturnType & dummy)
 * \brief Numerically evaluate integrals of TWO-variable functions over a quadrilateral domain.
 *
 * This subroutine maps the quadrilateral domain into a rectangle and integrate using an adaptive Lobatto rule
 *  
 * see also AdaptiveGaussianQuadrature()                                                                  
 * \param func the name of the function used to evaluate the integral
 * \param SW the south-west node of the quadrilateral
 * \param NW the north-west node of the quadrilateral
 * \param NE the north-east node of the quadrilateral
 * \param SE the south-east node of the quadrilateral
 * \param dummy used only to provide the return type
 * \param digits number of exact digits (there is a default value already provided!)
 **************************************************************************************************************************/
template<class FunctionType, class NodeType, class ReturnType>
inline ReturnType QuadrilateralQuadrature(FunctionType func,
					  const NodeType & SW, const NodeType & NW,
					  const NodeType & NE, const NodeType & SE,
					  int digits, const ReturnType & dummy){

  /* Integrate the new function over the square defined by (0,0) , (0,1) , (1,0) and (1,1) */
  return AdaptiveGaussianQuadrature(planar_bilinear_function_transformation(func,SW,NW,NE,SE,dummy), // tranform the function
 				    0.0,1.0,0.0,1.0,digits,dummy);
}

/**
 * \fn  ReturnType GaussLobattoQuadrilateralQuadrature(FunctionType func, NodeType & SW, NodeType & NW, NodeType & NE, 
 *  NodeType & SE, const ReturnType & dummy)
 * \brief Numerically evaluate integrals of TWO-variable functions over a quadrilateral domain.
 *
 * This subroutine maps the quadrilateral domain into a rectangle and integrate using a NON-adaptive Lobatto rule
 *  
 * see also GaussLobattoQuadrature()                                                                  
 * \param func the name of the function used to evaluate the integral
 * \param SW the south-west node of the quadrilateral
 * \param NW the north-west node of the quadrilateral
 * \param NE the north-east node of the quadrilateral
 * \param SE the south-east node of the quadrilateral
 * \param dummy used only to provide the return type
 **************************************************************************************************************************/
template<class FunctionType, class NodeType, class ReturnType>
inline ReturnType GaussLobattoQuadrilateralQuadrature(FunctionType func,
						      const NodeType & SW, const NodeType & NW,
						      const NodeType & NE, const NodeType & SE,
						      const ReturnType & dummy){

  /* Integrate the new function over the square defined by (0,0) , (0,1) , (1,0) and (1,1) */
  return GaussLobattoQuadrature(planar_bilinear_function_transformation(func,SW,NW,NE,SE,dummy), // tranform the function
				0.0,1.0,0.0,1.0,dummy);
}

/**
 * \fn  ReturnType Gauss5PointQuadrilateralQuadrature(FunctionType func, NodeType & SW, NodeType & NW, NodeType & NE, 
 *  NodeType & SE, const ReturnType & dummy)
 * \brief Numerically evaluate integrals of TWO-variable functions over a quadrilateral domain.
 *
 * This subroutine maps the quadrilateral domain into a rectangle and integrate using a NON-adaptive 5 point Gauss-Legendre rule
 *  
 * see also Gauss5PointQuadrature()
 * \param func the name of the function used to evaluate the integral
 * \param SW the south-west node of the quadrilateral
 * \param NW the north-west node of the quadrilateral
 * \param NE the north-east node of the quadrilateral
 * \param SE the south-east node of the quadrilateral
 * \param dummy used only to provide the return type
 **************************************************************************************************************************/
template<class FunctionType, class NodeType, class ReturnType> 
inline ReturnType Gauss5PointQuadrilateralQuadrature(FunctionType func,
						     const NodeType & SW, const NodeType & NW,
						     const NodeType & NE, const NodeType & SE,
						     const ReturnType & dummy){

  /* Integrate the new function over the square defined by (0,0) , (0,1) , (1,0) and (1,1) */
  return Gauss5PointQuadrature(planar_bilinear_function_transformation(func,SW,NW,NE,SE,dummy), // tranform the function
			       0.0,1.0,0.0,1.0,dummy);
}

/**
 * Generalized polynomial function of ONE-variable                                        
 * is a class of functions which have the form (x-xi)^n                                   
 ********************************************************/
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
      return std::pow((x-xi),n);
    }
};

/**
 * Generalized polynomial function of TWO-variables                                      
 * is a class of functions which have the form (x-xi)^n * (y-yi)^m                       
 ********************************************************************/
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
      return std::pow((x-xi),n)*std::pow((y-yi),m);
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

// ZeroLineIntegration()
double ZeroLineIntegration(const double & N1x, const double & N1y,
			   const double & N2x, const double & N2y);

/***********************************************************************//**
 * Compute the integral \f$ I = \int x dy \f$ along a segment line
 * in a local coordinate system (LCS).
 * GCS = global coordinate system
 *
 * \param N1x the x-coordinate of the segment first end point in the GCS.
 * \param N1y the y-coordinate of the segment first end point in the GCS.
 * \param N2x the x-coordinate of the segment second end point in the GCS.
 * \param N2y the y-coordinate of the segment second end point in the GCS.
 * \param xLCS the x-coordinate of the LCS in the GCS
 * \param yLCS the y-coordinate of the LCS in the GCS
 * \return the value of the integral
 ************************************************************************/
inline double ZeroLineIntegration(const double & N1x, const double & N1y,
				  const double & N2x, const double & N2y,
				  const double & xLCS, const double & yLCS){
  
  // Calculate the coordinates in the LCS
  double N1x_LCS(N1x-xLCS), N1y_LCS(N1y-yLCS);
  double N2x_LCS(N2x-xLCS), N2y_LCS(N2y-yLCS);
  
  return ZeroLineIntegration(N1x_LCS, N1y_LCS,
			     N2x_LCS, N2y_LCS);
}

// ZeroLineIntegration for Node input
template< class Node>
inline double ZeroLineIntegration(const Node& StartNode, const Node& EndNode){
  return ZeroLineIntegration(StartNode.x(), StartNode.y(), EndNode.x(), EndNode.y());
}

// PolynomLineIntegration()
double PolynomLineIntegration(const double & N1x, const double & N1y,
			      const double & N2x, const double & N2y,
			      const double & xCC, const double & yCC,
			      const int &OrderX, const int &OrderY);

// PolynomLineIntegration2()
double PolynomLineIntegration2(const double & N1x, const double & N1y,
			       const double & N2x, const double & N2y,
			       const double & xCC, const double & yCC,
			       const int &OrderX, const int &OrderY);

/***********************************************************************//**
 * Compute the integral \f$ I = \int (x - xc)^{(OrderX+1)} * (y - yc)^OrderY dy \f$
 * along a segment line in a local coordinate system (LCS).
 * The result of this polynomial function integration is determined with an analytic expression.
 * GCS = global coordinate system
 *
 * \param N1x the x-coordinate of the line first end point in the GCS.
 * \param N1y the y-coordinate of the line first end point in the GCS.
 * \param N2x the x-coordinate of the line second end point in the GCS.
 * \param N2y the y-coordinate of the line second end point in the GCS.
 * \param xCC the xc-coordinate in the GCS
 * \param yCC the yc-coordinate in the GCS
 * \param xLCS the x-coordinate of the LCS in the GCS
 * \param yLCS the y-coordinate of the LCS in the GCS
 * \return the value of the integrals up to 4th-order (i.e. OrderX + OrderY <= 4)
 *
 *********************************************************************************/
inline double PolynomLineIntegration(const double & N1x, const double & N1y,
				     const double & N2x, const double & N2y,
				     const double & xCC, const double & yCC,
				     const double & xLCS, const double & yLCS,
				     const int &OrderX, const int &OrderY){

  // Calculate the coordinates in the LCS
  double N1x_LCS(N1x-xLCS), N1y_LCS(N1y-yLCS);
  double N2x_LCS(N2x-xLCS), N2y_LCS(N2y-yLCS);
  double xCC_LCS(xCC-xLCS), yCC_LCS(yCC-yLCS);
  
  return PolynomLineIntegration(N1x_LCS, N1y_LCS,
				N2x_LCS, N2y_LCS,
				xCC_LCS, yCC_LCS,
				OrderX, OrderY);
}

/*!
 * \class GaussQuadratureData
 *
 * @brief Collection of absissae and weights 
 *        for n-point Gauss quadrature formula
 * 
 * The integral that can be computed with this data
 * is defined as follows:
 * 
 * \f$ I = \int_{a}^{b} f(x) dx = m \int_{-1}^{+1}f(c+mt) dt \f$
 * where \f$x=c+mt\f$, \f$c=\frac{1}{2}(b+a)\f$ and \f$m=\frac{1}{2}(b-a)\f$ ,
 *        'a' and 'b' are the integration limits. \n
 *
 * This integral can be writen as follows: \n
 * \f$ I = m \int_{-1}^{+1}f(c+mt) dt = m \sum_{i=i}^{n} \omega_{i} 
 *         f(c+mt_{i}) = (b-a) \sum_{i=i}^{n} GQn\_Weight[i]* f(a + GQn\_Abscissa[i]*(b-a)) \f$
 ***********************************************************************************************/
class GaussQuadratureData{
public:
  // Abscissae and weights for 1-point Gaussian method
  static const double GQ1_Abscissa[1];
  static const double GQ1_Weight[1];

  // Abscissae and weights for 2-point Gaussian method
  static const double GQ2_Abscissa[2];
  static const double GQ2_Weight[2];

  // Abscissae and weights for 3-point Gaussian method
  static const double GQ3_Abscissa[3];
  static const double GQ3_Weight[3];

  // Abscissae and weights for 5-point Gaussian method
  static const double GQ5_Abscissa[5];
  static const double GQ5_Weight[5];  

  //! Set the weights in the passed array
  static void getGaussQuadWeights(double * GaussQuadWeights, const int & NumberOfGQPs);

protected:
  GaussQuadratureData(void); //!< Private default constructor
  GaussQuadratureData(const GaussQuadratureData&); //!< Private copy constructor
  GaussQuadratureData& operator=(const GaussQuadratureData&); //!< Private assignment operator
  
};

inline void GaussQuadratureData::getGaussQuadWeights(double * GaussQuadWeights, const int & NumberOfGQPs){

  // Note: This subroutine doesn't check if there is enough memory allocated!

  switch ( NumberOfGQPs ){
  case 1:			// One point
    GaussQuadWeights[0] = GQ1_Weight[0];
    break;

  case 2:			// Two points
    GaussQuadWeights[0] = GQ2_Weight[0];
    GaussQuadWeights[1] = GQ2_Weight[1];
    break;
    
  case 3:			// Three points
    GaussQuadWeights[0] = GQ3_Weight[0];
    GaussQuadWeights[1] = GQ3_Weight[1];
    GaussQuadWeights[2] = GQ3_Weight[2];
    break;

  case 5:			// Five points
    GaussQuadWeights[0] = GQ5_Weight[0];
    GaussQuadWeights[1] = GQ5_Weight[1];
    GaussQuadWeights[2] = GQ5_Weight[2];
    GaussQuadWeights[3] = GQ5_Weight[3];
    GaussQuadWeights[4] = GQ5_Weight[4];
    break;

  default:
    // Add option if different GQPs are desired.
    throw runtime_error("GaussQuadratureData::getGaussQuadWeights() ERROR! There is no implementation for this number of Gauss points!");
  } // endswitch

}


//#include "System_routines.h"

/*******************************************************************************
 *
 * Routine: ridder
 * 
 * Copyright (C) 2007-2008 S. Guzik
 * 
 * Purpose
 * =======
 *
 *   Finds a root using Ridder's method.  This routine uses a pointer to a
 *   function to describe the equation.
 *
 * I/O
 * ===
 *
 *   f                  - (I) function pointer - returns function result
 *   bracketL           - (I) initial left bracket
 *   bracketR           - (I) initial right bracket
 *   maxIt              - (I) maximum number of iterations
 *   precision          - (I) precision of the root
 *   root               - (O) the root
 *   return             - >0 - success - iteration count
 *                        =0 - failure - maximum iterations reached or root
 *                                       not bracketed
 *
 ******************************************************************************/

template<typename T>
int ridder(T (*f)(T), const T bracketL, const T bracketR, const int maxIt,
           const int precision, T& root)

{


  /*==============================================================================
   * Declarations
   *============================================================================*/
  
  //----- Parameter Variables -----//
  
  const T toler = pow(10., -precision);
  // Tolerance for given precision
  
  //----- Local Variables -----//

  T a;                                 // Value at left ofinterval
  T b;                                 // Value at right of interval
  T m;                                 // Middle approximation
  T p;                                 // Ridder's smoothed approximation
  T p0;                                // Previous approximation
  T fa;                                // Function evaluated at a
  T fb;                                // Function evaluated at b
  T fm;                                // Function evaluated at m
  T fp;                                // Function evaluated at p

  /*==============================================================================
   * End of Declarations
   *============================================================================*/


  //--First results from bracket

  a = bracketL;
  fa = f(a);
  if ( fa == 0. ) {
    root = a;
    return 1;
  }
  b = bracketR;
  fb = f(b);
  if ( fb == 0. ) {
    root = b;
    return 1;
  }

  //--Confirm root bracketed

  p0 = 2.*b - a;

  //--Begin iterating

  int iter = 0;
  while ( iter != maxIt ) {
    ++iter;
    m = a + 0.5*(b - a);  // Middle approximation
    fm = f(m);
    if ( fm == 0. ) {     // Lucky solution
      root = m;
      return iter;
    }
    p = m + (m - a)*System::Copysign(1., fa)*fm/sqrt(fm*fm - fa*fb);
    // Check for sufficient precision in result
      if ( fabs(p - p0) < (fabs(p) + toler)*toler ) {
         root = p;
         return iter;
      }
      fp = f(p);
      if ( fp == 0. ) {  // Lucky solution
         root = p;
         return iter;
      }
      p0 = p;
      if ( fa*fm < 0. ) {
         if ( fa*fp < 0. ) {  // New interval from a to p
            b = p;
            fb = fp;
         }
         else {               // New interfal from p to m
            a = p;
            fa = fp;
            b = m;
            fb = fm;
         }
      }
      else {
	if ( fp*fb < 0. ) {  // New interval from p to b
	  a = p;
	  fa = fp;
	}
	else {               // New interval from m to p
	  a = m;
	  fa = fm;
	  b = p;
	  fb = fp;
	}
      }
  }
  root = p;
  return 0;  // Error - max iterations reached
   
}


/*******************************************************************************
 *
 * Routine: ridder
 *
 * Copyright (C) 2007-2008 S. Guzik
 *
 * Purpose
 * =======
 *
 *   Finds a root using Ridder's method - this routine uses a functor to
 *   describe the equation.
 *
 * I/O
 * ===
 *
 *   f                  - (I) function object - operator() evaluates
 *   bracketL           - (I) initial left bracket
 *   bracketR           - (I) initial right bracket
 *   maxIt              - (I) maximum number of iterations
 *   precision          - (I) precision of the root
 *   root               - (O) the root
 *   return             - >0 - success - iteration count
 *                        =0 - failure - maximum iterations reached or root not
 *                                       bracketed
 *
 ******************************************************************************/

template<typename T, typename F>
int ridder(F& f, const T bracketL, const T bracketR, const int maxIt,
           const int precision, T& root)
  
{
  

  /*==============================================================================
   * Declarations
   *============================================================================*/
  
  //----- Parameter Variables -----//

  const T toler = pow(10., -precision);
  // Tolerance for given precision

  //----- Local Variables -----//

  T a;                                 // Value at left ofinterval
  T b;                                 // Value at right of interval
  T m;                                 // Middle approximation
  T p;                                 // Ridder's smoothed approximation
  T p0;                                // Previous approximation
  T fa;                                // Function evaluated at a
  T fb;                                // Function evaluated at b
  T fm;                                // Function evaluated at m
  T fp;                                // Function evaluated at p

  /*==============================================================================
   * End of Declarations
   *============================================================================*/


  //--First results from bracket

  a = bracketL;
  fa = f(a);
  if ( fa == 0. ) {
    root = a;
    return 1;
  }
  b = bracketR;
  fb = f(b);
  if ( fb == 0. ) {
    root = b;
    return 1;
  }

  //--Confirm root bracketed

  p0 = 2.*b - a;

  //--Begin iterating

  int iter = 0;
  while ( iter != maxIt ) {
    ++iter;
    m = a + 0.5*(b - a);  // Middle approximation
    fm = f(m);
    if ( fm == 0. ) {     // Lucky solution
      root = m;
      return iter;
    }
    p = m + (m - a)*System::Copysign(1., fa)*fm/sqrt(fm*fm - fa*fb);
    // Check for sufficient precision in result
    if ( fabs(p - p0) < (fabs(p) + toler)*toler ) {
      root = p;
      return iter;
    }
    fp = f(p);
    if ( fp == 0. ) {  // Lucky solution
      root = p;
      return iter;
    }
    p0 = p;
    if ( fa*fm < 0. ) {
      if ( fa*fp < 0. ) {  // New interval from a to p
	b = p;
	fb = fp;
      }
      else {               // New interfal from p to m
	a = p;
	fa = fp;
	b = m;
	fb = fm;
      }
    }
    else {
      if ( fp*fb < 0. ) {  // New interval from p to b
	a = p;
	fa = fp;
      }
      else {               // New interval from m to p
	a = m;
	fa = fm;
	b = p;
	fb = fp;
      }
    }
  }
  root = p;
  return 0;  // Error - max iterations reached

}


/**************** Function Prototypes ********************/

/*### double qgauss10(const FunctionType1D , const double , const double )
  ########################################################################
  Returns the integral of the function "func" between a and b,
  by ten-point Gauss-Legendre integration: the function is evaluated exactly
  ten times at interior points in the range of integration.
  The integral is exact for polynomials up to order of 19. 
  Implementation based on the subroutine from Numerical Recipes*/
double qgauss10(const FunctionType1D func, const double a, const double b);


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

#endif // _NumericalLibrary_INCLUDED

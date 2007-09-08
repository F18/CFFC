/**********************************************************************
 **********************************************************************
 **                                                                  **
 ** File: FieldData.h                                                **
 **                                                                  **
 ** Description: This header file contains a functor implementation  **
 **              required to specify prescribed field data.          **
 **              Basically, a base class TFunctor provides a virtual **
 **              funcion named 'Call' or a virtually overloaded      **
 **              operator '()' with which you will be able to call   ** 
 **              the member function.  From the base class, a        **
 **              template class TSpecificFunctor which is initialized**
 **              with a pointer to an object and a pointer to a      **
 **              member function in its constructor.  The derived    **
 **              class overrides the function 'Call' and/or the      **
 **              operator '()' of the base class.                    **
 **                                                                  **
 **                                                                  **
 ** Author: Marc "T-Bone" Charest                                    **
 **                                                                  **
 ** Revision:  Date        Initials   Change                         **
 **            05/09/2007  MRC        Original creation              **
 **                                                                  **
 **********************************************************************
 **********************************************************************/
#ifndef _FIELD_DATA_INCLUDED
#define _FIELD_DATA_INCLUDED

// required includes
#include "../Math/Vector2D.h"


/***********************************************************************/
/*!
 * Class: TFunctor
 *
 * @brief Abstract Base Class.
 *
 * \verbatim
 * Member functions
 *     operator() -- Return field data at specific position
 *     Call()     -- Return field data at specific position
 * \endverbatim
 */
/***********************************************************************/
class TFunctor {
 public:

  // two possible functions to call member function. virtual cause derived
  // classes will use a pointer to an object and a pointer to a member function
  // to make the function call
  virtual double operator()(const Vector2D &r)=0;  // call using operator
  virtual double Call(const Vector2D &r)=0;        // call using function
};


/***********************************************************************/
/*!
 * Class: TSpecificFunctor
 *
 * @brief Derived template class.
 *
 * \verbatim
 * Objects
 *     fpt        -- function pointer
 *     pt2Object  -- pointer to object
 *
 * Member functions
 *     operator() -- Return field data at specific position
 *     Call()     -- Return field data at specific position
 * \endverbatim
 */
/***********************************************************************/
template <class TClass> class TSpecificFunctor : public TFunctor {
 private:
  double (TClass::*fpt)(const Vector2D &);   // pointer to member function
  TClass* pt2Object;                         // pointer to object
  
 public:
  
  // constructor - uninitialized
  TSpecificFunctor() : pt2Object(NULL), fpt(NULL) {}

  // constructor - takes pointer to an object and pointer to a member and stores
  // them in two private variables
  TSpecificFunctor(TClass* _pt2Object, double (TClass::*_fpt)(const Vector2D &r))
  { pt2Object = _pt2Object;  fpt=_fpt; };

  // initializer
  void Set(TClass* _pt2Object, double (TClass::*_fpt)(const Vector2D &r))
  { pt2Object = _pt2Object;  fpt=_fpt; };


  // override operator "()"
  virtual double operator()(const Vector2D &r)
  { return (*pt2Object.*fpt)(r);};              // execute member function
  
  // override function "Call"
  virtual double Call(const Vector2D &r)
  { return (*pt2Object.*fpt)(r);};             // execute member function

};


/***********************************************************************/
/*!
 * Class: FieldData
 *
 * @brief Derived class container for prescribed field data
 *
 * \verbatim
 * Objects
 *     fpt        -- function pointer
 *     pt2Object  -- pointer to object
 *
 * Member functions
 *     operator() -- Return field data at specific position
 *     Call()     -- Return field data at specific position
 * \endverbatim
 */
/***********************************************************************/
class FieldData {

 public:

  //objects
  double val; // constant Value

  // constructors
  FieldData():val(0) {};
  FieldData(const int &x):val(x) {};

  // Functions to describe the field
  double Constant(const Vector2D &r) { return val; };
};

#endif //_FIELD_DATA_INCLUDED

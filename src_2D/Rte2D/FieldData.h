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
#include "../Grid/Grid2DQuad.h"
#include "../CFD/CFD.h"

/********************************************************
 * Struct definitions                                   *
 ********************************************************/
struct ConstantFieldParams { double val; };
static const ConstantFieldParams cfDefValues = { 0 };

struct DiscreteFieldParams {   
  Grid2D_Quad_Block** Grid; // pointer to multiblock grid object
  int Number_of_Blocks_Idir;// number of blocks in i-dir
  int Number_of_Blocks_Jdir;// number of blocks in j-dir
  double** valarray;        // values corresponding to grid cell centers
  // bool isNodalValues;       // a flag, true if valarray is defined at the 'Grid' nodes
  //                           // or false if defined at cell centers.
};
static const DiscreteFieldParams dfDefValues = { NULL, 0, 0, NULL };



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

 private:

  //
  //objects
  //
  ConstantFieldParams CF;
  DiscreteFieldParams DF;

 public:

  //
  // constructors
  //
  FieldData() : CF(cfDefValues), DF(dfDefValues) {};
  void SetConstantParams( const struct ConstantFieldParams params ) { CF = params; };
  void SetDiscreteParams( const struct DiscreteFieldParams params ) { DF = params; };

  //
  // Functions to describe the field
  //
  // return a constant field
  double Constant(const Vector2D &r) { return CF.val; };
  // return a value interpolated from a 2d grid
  double Interpolate(const Vector2D &r);
};


/********************************************************
 * Compute the value at the specified location by       *
 * interpolating from an existing grid.                 *
 ********************************************************/
inline double FieldData :: Interpolate(const Vector2D &r) 
{
  // declares
  int err;
  int ib, jb;   // block indices
  int ic, jc;   // cell indices
  double value; // the value we are looking for!!!!
  
  // check to make sure parameter pointers have been properly defined
  if (DF.Grid==NULL || DF.valarray==NULL) {
    cerr << "FieldData::Interpolate() - Need to define grid/valarray.\n";
    exit(-1);
  } // endif


  //
  // search the multi-block grid to determine the containing block
  //
  err =  Search_Multi_Block_Grid(DF.Grid,
				 DF.Number_of_Blocks_Idir,
				 DF.Number_of_Blocks_Jdir,
				 r,
				 ib, jb);
  // if there was an error, exit
  if (err) {
    cerr << "FieldData::Interpolate() - Error searching multiblock grid.\n";
    exit(-1);
  } // endif

  
  //
  // search the known block for the cell in which the location falls
  //
  err =  Search_Single_Block_Grid(DF.Grid[ib][jb],
				  r,
				  ic, jc);
  // if there was an error, exit
  if (err) {
    cerr << "FieldData::Interpolate() - Error searching multiblock grid.\n";
    exit(-1);
  } // endif


  //
  // interpolate
  //
  err = Bilinear_Interpolation(DF.valarray[ic  ][jc  ], DF.Grid[ib][jb].Node[ic  ][jc  ].X,
			       DF.valarray[ic  ][jc+1], DF.Grid[ib][jb].Node[ic  ][jc+1].X,
			       DF.valarray[ic+1][jc+1], DF.Grid[ib][jb].Node[ic+1][jc+1].X,
			       DF.valarray[ic+1][jc  ], DF.Grid[ib][jb].Node[ic+1][jc  ].X,
			       r, value);
   // if there was an error, exit
  if (err) {
    cerr << "FieldData::Interpolate() - Error interpolating.\n";
    exit(-1);
  } // endif
 

  // return the value
  return value;
}

#endif //_FIELD_DATA_INCLUDED

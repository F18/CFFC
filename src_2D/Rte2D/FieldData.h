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
#include "../AMR/AdaptiveBlock.h"
#include "../AMR/QuadTree.h"

/********************************************************
 * Struct definitions                                   *
 ********************************************************/
template<class Soln_State> 
struct ConstantFieldParams { Soln_State U; };

template<class Quad_Soln_Block> 
struct DiscreteFieldParams {   
  AdaptiveBlockResourceList const *List_of_Global_Solution_Blocks;
  AdaptiveBlock2D_List const      *List_of_Local_Solution_Blocks;
  Quad_Soln_Block const          **Local_SolnBlk;
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
template<class Soln_State, class Quad_Soln_Block> 
class FieldData {

 private:

  //
  //objects
  //
  // data structs
  ConstantFieldParams<Soln_State> CF;
  DiscreteFieldParams<Quad_Soln_Block> DF;

  // pointer to member function
  Soln_State (FieldData::*fpt)(const Vector2D &);


 public:

  //
  // constructors
  //
  FieldData() {
    CF.U.Zero();
    DF.List_of_Global_Solution_Blocks = NULL;
    DF.List_of_Local_Solution_Blocks = NULL; 
    DF.Local_SolnBlk = NULL; 
  };
  
  void SetConstantParams( const struct ConstantFieldParams<Soln_State> params ) 
  { 
    CF = params; 
    fpt = &FieldData<Soln_State,Medium2D_Quad_Block>::Constant;
  };
  
  void SetDiscreteParams( const struct DiscreteFieldParams<Quad_Soln_Block> params ) 
  { 
    DF = params; 
    fpt = &FieldData<Soln_State,Medium2D_Quad_Block>::Interpolate;
  };

  //
  // Functor overloads
  //
  // Two possible functions to call member function:
  // 1. call using operator
  Soln_State operator()(const Vector2D &r) { return (*this.*fpt)(r); };   

  // 2. call using function
  Soln_State Call(const Vector2D &r) { return (*this.*fpt)(r); };   

  
  //
  // Functions to describe the field
  //
  // return a constant field
  Soln_State Constant(const Vector2D &r) { return CF.U; };
  // return a value interpolated from a 2d grid
  Soln_State Interpolate(const Vector2D &r);
};


/********************************************************
 * Compute the value at the specified location by       *
 * interpolating from an existing grid.                 *
 ********************************************************/
template <class Soln_State, class Quad_Soln_Block>
inline Soln_State FieldData<Soln_State,Quad_Soln_Block> :: Interpolate(const Vector2D &r) 
{

  // declares
  int err;
  int nb;       // block index
  int ic, jc;   // cell indices
  double value; // the value we are looking for!!!!
  
  // check to make sure parameter pointers have been properly defined
  if ( DF.List_of_Global_Solution_Blocks == NULL || 
       DF.List_of_Local_Solution_Blocks == NULL ||
       DF.Local_SolnBlk ==NULL ) {
    cerr << "FieldData::Interpolate() - Need to define parameter pointers.\n";
    exit(-1);
  } // endif

  //
  // search the multi-block grid to determine the containing block
  //
  err =  Search_Multi_Block_Grid(DF.List_of_Global_Solution_Blocks,
				 DF.List_of_Local_Solution_Blocks,
				 ,
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
*/
}

#endif //_FIELD_DATA_INCLUDED

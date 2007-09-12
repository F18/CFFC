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
  int nb;           // block index
  int i, j;       // cell indices
  Polygon P;        // polygon
  Soln_State value; // interpolated result

  // check to make sure parameter pointers have been properly defined
  if ( DF.List_of_Global_Solution_Blocks == NULL || 
       DF.List_of_Local_Solution_Blocks == NULL ||
       DF.Local_SolnBlk ==NULL ) {
    cerr << "FieldData::Interpolate() - Need to define parameter pointers.\n";
    exit(-1);
  } // endif


  //------------------------------------------------
  // search the multi-block grid to determine the containing block
  //------------------------------------------------
  bool block_found(false); // bool, true if block containing point is found
  for (nb=0; nb<DF.List_of_Local_Solution_Blocks->Nblk; nb++) {
    
    // build a polygon
    P.convert(DF.Local_SolnBlk[nb]->Grid.nodeSW(i,j).X,
	      DF.Local_SolnBlk[nb]->Grid.nodeSE(i,j).X,
	      DF.Local_SolnBlk[nb]->Grid.nodeNE(i,j).X,
	      DF.Local_SolnBlk[nb]->Grid.nodeNW(i,j).X);

    // is the point in the polygon, get out
    if (P.point_in_polygon(r)) { 
      block_found = true; 
      break;
    } // endif

  } // endfor

  // did we find it?
  if (!block_found) {
    cerr << "FieldData::Interpolate() - No containing block found.\n";
    exit(-1);
  } // endif

  
  //------------------------------------------------
  // search the sinlge-block grid to determine the containing cell
  //------------------------------------------------
  bool cell_found(false); // bool, true if cell containing point is found
  for (j = DF.Local_SolnBlk[nb]->Grid.JCl - DF.Local_SolnBlk[nb]->Grid.Nghost; 
       j <= DF.Local_SolnBlk[nb]->Grid.JCu + DF.Local_SolnBlk[nb]->Grid.Nghost && !cell_found; 
       j++) {
    for (i = DF.Local_SolnBlk[nb]->Grid.ICl - DF.Local_SolnBlk[nb]->Grid.Nghost; 
	 i <= DF.Local_SolnBlk[nb]->Grid.ICu + DF.Local_SolnBlk[nb]->Grid.Nghost; 
	 i++) {

    // build a polygon
    P.convert(DF.Local_SolnBlk[nb]->Grid.nodeSW(i,j).X,
	      DF.Local_SolnBlk[nb]->Grid.nodeSE(i,j).X,
	      DF.Local_SolnBlk[nb]->Grid.nodeNE(i,j).X,
	      DF.Local_SolnBlk[nb]->Grid.nodeNW(i,j).X);

    // is the point in the polygon, get out
    if (P.point_in_polygon(r)) { 
      cell_found = true; 
      break;
    } // endif


    } // endfor
  } // endfor
  
  // did we find it?
  if (!cell_found) {
    cerr << "FieldData::Interpolate() - No containing cell found.\n";
    exit(-1);
  } // endif


  //------------------------------------------------
  // interpolate
  //------------------------------------------------
  err = Bilinear_Interpolation(DF.Local_SolnBlk[nb]->UnSW(i,j), DF.Local_SolnBlk[nb]->Grid.nodeSW(i,j).X,
                               DF.Local_SolnBlk[nb]->UnNW(i,j), DF.Local_SolnBlk[nb]->Grid.nodeNW(i,j).X,
                               DF.Local_SolnBlk[nb]->UnNE(i,j), DF.Local_SolnBlk[nb]->Grid.nodeNE(i,j).X,
                               DF.Local_SolnBlk[nb]->UnSE(i,j), DF.Local_SolnBlk[nb].Grid.nodeSE(i,j).X,
                               r, value);
  if (err) {
    cerr << "FieldData::Interpolate() - Error interpolating.\n";
    exit(-1);
  } // endif

  // return the value
  return value;

}

#endif //_FIELD_DATA_INCLUDED

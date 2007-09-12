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
  QuadTreeBlock_DataStructure *QuadTree;
  AdaptiveBlockResourceList   *List_of_Global_Solution_Blocks;
  AdaptiveBlock2D_List        *List_of_Local_Solution_Blocks;
  Quad_Soln_Block            **Local_SolnBlk;
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
 *     operator() -- Return figecific position
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
  Soln_State (FieldData::*fpt)(const void *);


 public:

  //
  // constructors
  //
  FieldData() {
    CF.U.Zero();
    DF.QuadTree = NULL;
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
  int i, j, nb;     //  cell indexes, block index
  Polygon P;        // polygon
  Soln_State value; // interpolated result
  LinkedList<Vector2D> pdata;  // linked list object for polygon
  int foundIt(0);   // bool, true if block containing point is found
  int iFoundIt(0);  // bool, true if this processor found it
  int buffer_size;  // buffer size 
  int* buffer(NULL);// buffer

  // check to make sure parameter pointers have been properly defined
  if ( DF.QuadTree == NULL || 
       DF.List_of_Global_Solution_Blocks == NULL || 
       DF.List_of_Local_Solution_Blocks == NULL ||
       DF.Local_SolnBlk ==NULL ) {
    cerr << "FieldData::Interpolate() - Need to define parameter pointers.\n";
    exit(-1);
  } // endif


  //------------------------------------------------
  // search the multi-block grid to determine the containing block
  //------------------------------------------------

  //
  // Loop over available processors, and all available blocks
  //
  for ( int iCPU=0; iCPU<DF.QuadTree->Ncpu && !foundIt; iCPU++ ) {
    for ( nb=0; nb<DF.QuadTree->Nblk && !foundIt; nb++ ) {

      // don't look at unallocated blocks
      if (DF.QuadTree->Blocks[iCPU][nb] == NULL) break;
      
      // only look at used blocks
      if ( !DF.QuadTree->Blocks[iCPU][nb]->block.used ) break;
      
      // if this block belongs to this processor
      if ( DF.List_of_Local_Solution_Blocks->ThisCPU == iCPU ) {

	//
	// Build a polygon.  Use each node point along the outside of the block.
	// Remmeber that we build the polygon in a counter-clockwise rotation
	//
	// south boudnary
	for (i = DF.Local_SolnBlk[nb]->Grid.INl - DF.Local_SolnBlk[nb]->Grid.Nghost; 
	     i <= DF.Local_SolnBlk[nb]->Grid.INu + DF.Local_SolnBlk[nb]->Grid.Nghost; 
	     i++) {
	  j = DF.Local_SolnBlk[nb]->Grid.JNl - DF.Local_SolnBlk[nb]->Grid.Nghost;
	  pdata.add( DF.Local_SolnBlk[nb]->Grid.Node[i][j].X );
	} // endfor - south
	
	// east boundary ( skip the first node on this side )
	for (j = DF.Local_SolnBlk[nb]->Grid.JNl - DF.Local_SolnBlk[nb]->Grid.Nghost + 1; 
	     j <= DF.Local_SolnBlk[nb]->Grid.JNu + DF.Local_SolnBlk[nb]->Grid.Nghost; 
	     j++) {
	  i = DF.Local_SolnBlk[nb]->Grid.JNu + DF.Local_SolnBlk[nb]->Grid.Nghost;
	  pdata.add( DF.Local_SolnBlk[nb]->Grid.Node[i][j].X );
	} // endfor - east
	
	// north boudnary ( skip the first node on this side )
	for (i = DF.Local_SolnBlk[nb]->Grid.INu + DF.Local_SolnBlk[nb]->Grid.Nghost - 1; 
	     i >= DF.Local_SolnBlk[nb]->Grid.INl - DF.Local_SolnBlk[nb]->Grid.Nghost; 
	     i++) {
	  j = DF.Local_SolnBlk[nb]->Grid.JNu + DF.Local_SolnBlk[nb]->Grid.Nghost;
	  pdata.add( DF.Local_SolnBlk[nb]->Grid.Node[i][j].X );
	} // endfor - south

	// west boundary ( skip the first & last node on this side )
	for (j = DF.Local_SolnBlk[nb]->Grid.JNu + DF.Local_SolnBlk[nb]->Grid.Nghost - 1; 
	     j >= DF.Local_SolnBlk[nb]->Grid.JNl - DF.Local_SolnBlk[nb]->Grid.Nghost + 1; 
	     j++) {
	  i = DF.Local_SolnBlk[nb]->Grid.JNl - DF.Local_SolnBlk[nb]->Grid.Nghost;
	  pdata.add( DF.Local_SolnBlk[nb]->Grid.Node[i][j].X );
	} // endfor - east

	// convert polynomial
	P.convert(pdata);
	
	// is the point in the polygon
	if (P.point_in_polygon(r)) iFoundIt = iCPU+1;  // offset by 1

	// Broadcast the result to everyone, taking the max value.
	// This will be the processor with the lowest rank 
	// (ensures only one processor can find it!!).
#ifdef _MPI_VERSION
	MPI::COMM_WORLD.Allreduce(/*snd*/&iFoundIt, /*rcv*/&foundIt, /*size*/1, 
				  /*typ*/MPI::INT, /*mx vlu*/MPI::MAX);
#endif

	// if somebody found it, exit
	if (foundIt) break;

      } // endif - ThisCPU

    } // endfor - block
  } // endfor - iCPU


  //
  // did we find it?
  //
  if (!foundIt) {
    cerr << "FieldData::Interpolate() - error searching for point.\n";
    exit(-1);
  } // endif


  // un offset the cpu index (makes it easier to read)
  iFoundIt -= 1;


  //------------------------------------------------
  // Search block && Interpolate
  //------------------------------------------------
  // if we got this far, someone must have found the containing block
  if (iFoundIt==DF.List_of_Local_Solution_Blocks->ThisCPU) {

    //
    // search block
    //
    err = Seach_Mesh(DF.Local_SolnBlk[nb]->Grid, r, i, j );
    if (err) {
      cerr << "FieldData::Interpolate() - No containing cell found.\n";
      exit(-1);
    } // endif


    //
    // Interpolate
    //
    err = Bilinear_Interpolation(DF.Local_SolnBlk[nb]->UnSW(i,j), 
				 DF.Local_SolnBlk[nb]->Grid.nodeSW(i,j).X,
				 DF.Local_SolnBlk[nb]->UnNW(i,j), 
				 DF.Local_SolnBlk[nb]->Grid.nodeNW(i,j).X,
				 DF.Local_SolnBlk[nb]->UnNE(i,j), 
				 DF.Local_SolnBlk[nb]->Grid.nodeNE(i,j).X,
				 DF.Local_SolnBlk[nb]->UnSE(i,j), 
				 DF.Local_SolnBlk[nb]->Grid.nodeSE(i,j).X,
				 r, value);
    if (err) {
      cerr << "FieldData::Interpolate() - Error interpolating.\n";
      exit(-1);
    } // endif

  } // endif - iFoundIt


  //------------------------------------------------
  // broadcast result
  //------------------------------------------------
  // if we got this far, we must have found the containing cell
#ifdef _MPI_VERSION

  // set the buffer size for the state variables
  buffer_size = DF.Local_SolnBlk[nb]->NumVar();

  // if this cpu is the one sending it, load the buffer
  if (iFoundIt == DF.List_of_Local_Solution_Blocks->ThisCPU) {
    buffer = new double[buffer_size];
    for (int k=0; k<buffer_size; k++) buffer[k] = value[k];
  }

  // send
  MPI::COMM_WORLD.Bcast(buffer, buffer_size, MPI::DOUBLE, iFoundIt);

  // for all the other processors, unload the buffer
  if (iFoundIt != DF.List_of_Local_Solution_Blocks->ThisCPU) {
    buffer = new double[buffer_size];
    for (int k=0; k<buffer_size; k++) value[k] = buffer[k];
  }

  // clean up memory
  if (buffer!=NULL) { delete[] buffer; buffer = NULL; }

#endif


  //------------------------------------------------
  // return the value
  //------------------------------------------------
  return value;

}

#endif //_FIELD_DATA_INCLUDED

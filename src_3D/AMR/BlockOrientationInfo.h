
/* BlockOrientationInfo.h: Header file defining block orientation information class. */

#ifndef _BLOCK_ORIENTATION_INFO_INCLUDED
#define _BLOCK_ORIENTATION_INFO_INCLUDED

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>

using namespace std;

/* Include blkclib header files. */

#include "coord_transform.h"
#include "block_connectivity.h"

/* Include required CFFC header files. */

#ifndef _MPI_INCLUDED
#include "../MPI/MPI.h"
#endif // _MPI_INCLUDED

#define MAX_BOUNDARY_ELEMENTS_FOR_A_BLOCK          27

/*******************************************************************************
 *
 * Integeer direction numbers corresponding to neighbour directions.
 *
 ******************************************************************************/
namespace DIRECTION
{
   enum Direction_Number {
      BSW = 0,
      SW  = 1,
      TSW = 2,
      BW  = 3,
      W   = 4,
      TW  = 5,
      BNW = 6,
      NW  = 7,
      TNW = 8,
      BS  = 9,
      S   = 10,
      TS  = 11,
      B   = 12,
      T   = 13,
      BN  = 14,
      N   = 15,
      TN  = 16,
      BSE = 17,
      SE  = 18,
      TSE = 19,
      BE  = 20,
      E   = 21,
      TE  = 22,
      BNE = 23,
      NE  = 24,
      TNE = 25
   };
}  

/*******************************************************************************
 *
 * Boundary element numbers and corresponding neighbour directions.
 *
 ******************************************************************************/
namespace BE
{
   enum Boundary_Elements {
      BSW = 0,
      SW  = 1,
      TSW = 2,
      BW  = 3,
      W   = 4,
      TW  = 5,
      BNW = 6,
      NW  = 7,
      TNW = 8,
      BS  = 9,
      S   = 10,
      TS  = 11,
      B   = 12,
      ME  = 13,
      T   = 14,
      BN  = 15,
      N   = 16,
      TN  = 17,
      BSE = 18,
      SE  = 19,
      TSE = 20,
      BE  = 21,
      E   = 22,
      TE  = 23,
      BNE = 24,
      NE  = 25,
      TNE = 26
   };
}  

//=============================================================
class Direction_Indices {
  private:
  public:
    int         i,j,k;  // Adaptive block dimensions.
    
    /* Creation, copy, and assignment constructors. */
    Direction_Indices(void) : i(0), j(0), k(0) { }

    Direction_Indices(const Direction_Indices &DI) {
      i = DI.i; j = DI.j; k = DI.k;
    }

    Direction_Indices (const int i_dir,
	               const int j_dir,
	               const int k_dir) : i(i_dir), j(j_dir), k(k_dir) { }

    /* Destructor. */
    // ~AdaptiveBlock3D_Dimensions(void);
    // Use automatically generated destructor.

    /* Return indices given a direction number. */
    Direction_Indices direction_number_to_direction_indices(const int direction_number) const;     

    /* Return direction number given the indices. */
    int direction_indices_to_direction_number(void) const;
    int direction_indices_to_direction_number(const Direction_Indices &DI) const;

    /* Return direction indices given boundary element number. */
    Direction_Indices boundary_element_number_to_direction_indices(const int boundary_element_number) const;

    /* Return boundary element number given the indices. */
    int direction_indices_to_boundary_element_number(void) const;
    int direction_indices_to_boundary_element_number(const Direction_Indices &DI) const;

    /* Return indices of node for referencing the geometry. */
    Direction_Indices reference_geomety_node(const int boundary_element_number,
                                             const int i_low,
                                             const int i_up,
                                             const int j_low,
                                             const int j_up,
                                             const int k_low,
                                             const int k_up,
                                             const int nghost) const;

    /* Index operator. */

   int &operator[](int index);
   //! Index operator
   const int &operator[](int index) const;

    /* Assignment operator. */
    Direction_Indices & operator =(const Direction_Indices &DI);

    /* Unary subtraction operators. */
    friend Direction_Indices operator -(const Direction_Indices &DI);
  
    /* Equal relational operator . */
    friend int operator ==(const Direction_Indices &DI1, 
                           const Direction_Indices &DI2);

    /* Not equal relational operator. */
    friend int operator !=(const Direction_Indices &DI1, 
                           const Direction_Indices &DI2);

    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file,
				 const Direction_Indices &DI);
    friend istream &operator >> (istream &in_file,
				 Direction_Indices &DI);
};

/*************************************************************
 * Direction_Indices -- Index operators.                     *
 *************************************************************/
inline int &Direction_Indices::operator[](int index) {
   switch(index){  
   case 1:
      return i;    
   case 2:
      return j;
   case 3:
      return k;
   default :
      return i;
   };
}

inline const int &Direction_Indices::operator[](int index) const {
   switch(index){  
   case 1:
      return i;    
   case 2:
      return j;
   case 3:
      return k;
   default :
      return i;
   };
}

/*************************************************************
 * Direction_Indices -- Assignment operator.                 *
 *************************************************************/
inline Direction_Indices &Direction_Indices::
operator =(const Direction_Indices &DI) {
   i = DI.i; j = DI.j; k = DI.k;
}

/*************************************************************
 * Direction_Indices -- Unary arithmetic operators.          *
 *************************************************************/
inline Direction_Indices operator -(const Direction_Indices &DI) {
   return Direction_Indices(-DI.i, -DI.j, -DI.k);
}

/*************************************************************
 * Direction_Indices -- Relational operators.                *
 *************************************************************/
inline int operator ==(const Direction_Indices &DI1, 
                       const Direction_Indices &DI2) {
   return (DI1.i == DI2.i && DI1.j == DI2.j && DI1.k == DI2.k);
}

inline int operator !=(const Direction_Indices &DI1, 
                       const Direction_Indices &DI2) {
   return (DI1.i != DI2.i || DI1.j != DI2.j || DI1.k != DI2.k);
}

/*************************************************************
 * Direction_Indices -- Input-output operators.              *
 *************************************************************/
inline ostream &operator << (ostream &out_file,
			     const Direction_Indices &DI) {
  out_file << " " << DI.i << " " << DI.j << " " << DI.k;
  return (out_file);
}

inline istream &operator >> (istream &in_file,
			     Direction_Indices &DI) {
  in_file.setf(ios::skipws);
  in_file >> DI.i >> DI.j >> DI.k;
  in_file.unsetf(ios::skipws);
  return (in_file);
}

//=============================================================
class Block_Orientation_Info {
  public:
   
    int                          ctm_offsets[6]; // compact transformation matrix (3 components) and offsets (3 elements)
    int                          direction_me_to_neighbour[3]; 
    int                          direction_neighbour_to_me[3]; 
    
    // constructor
    Block_Orientation_Info(void){
       ctm_offsets[0] = 1 ; 
       ctm_offsets[1] = 1 ;  
       ctm_offsets[2] = 1 ;
       ctm_offsets[3] = 0 ; 
       ctm_offsets[4] = 0 ;  
       ctm_offsets[5] = 0  ;  
       direction_me_to_neighbour[0] = 0 ; 
       direction_me_to_neighbour[1] = 0 ;
       direction_me_to_neighbour[2] = 0 ; 
       direction_neighbour_to_me[0] = 0 ; 
       direction_neighbour_to_me[1] = 0 ;
       direction_neighbour_to_me[2] = 0 ; 
    }
   
    // use default destructor 
   
    // member functions
    void set_block_orientation_info(BlkC::BlockConnectivity &blkc, 
                                    const int iblk, 
                                    const int which_neighbour,
                                    const int my_i, 
                                    const int my_j, 
                                    const int my_k, 
                                    int &iblk_neigh);

    int compute_message_tag(const int i_index, 
                            const int j_index,  
                            const int k_index);

    int convert_boundary_elements_from_ijk_to_orientations(const int my_i, 
                                                           const int my_j, 
                                                           const int my_k);

    void my_index(int &i_index, 
                  int &j_index, 
                  int &k_index);

    void neighbour_index(int &i_index, 
                         int &j_index, 
                         int &k_index);

    void broadcast(void);
   
    /* Assignment operator. */

    Block_Orientation_Info  & operator =(const Block_Orientation_Info &BlkOI); 
   
    /* Input-output operators. */

    friend ostream& operator << (ostream &out_file, 
                                 const Block_Orientation_Info &BlkOI);
   
    friend istream& operator >> (istream &in_file,  
                                 Block_Orientation_Info &BlkOI);
};

//----------------- Assignment ----------------------------//
inline Block_Orientation_Info &Block_Orientation_Info::
operator =(const Block_Orientation_Info &BlkOI) {
   //self assignment protection
   if (this != &BlkOI) {   
      //copy assignment
      for (int i = 0; i < 6; i++) {
         ctm_offsets[i] = BlkOI.ctm_offsets[i];
      } // ctm and offsets   
      for (int i = 0; i < 3; i++) {
         direction_me_to_neighbour[i] = BlkOI.direction_me_to_neighbour[i];
         direction_neighbour_to_me[i] = BlkOI.direction_neighbour_to_me[i];
      } // directions
   } /* endif */
   return (*this);
}

/**********************************************************************************
 * Block_Orientation_Info -- Input-output operators.                              *
 **********************************************************************************/
inline ostream &operator << (ostream &out_file, 
                             const Block_Orientation_Info &BlkIO) {
   out_file.precision(10);
   out_file.setf(ios::scientific);
   for (int i = 0; i < 6; i++) {
      out_file << " " << BlkIO.ctm_offsets[i];
   } // output ctm and offsets
   for (int i = 0; i<3; i++) {
      out_file << " " << BlkIO.direction_me_to_neighbour[i];
   } // directions
   for(int i = 0; i<3; i++) {
      out_file << " " << BlkIO.direction_neighbour_to_me[i];
   } // output directions
   out_file.unsetf(ios::scientific);
   return (out_file);
}

inline istream &operator >> (istream &in_file, 
                             Block_Orientation_Info &BlkIO) {
   in_file.precision(10);
   in_file.setf(ios::scientific);
   for (int i = 0; i < 6; i++) {
      in_file >> BlkIO.ctm_offsets[i];
   } // output ctm and offsets
   for (int i = 0; i < 3; i++) {
      in_file >> BlkIO.direction_me_to_neighbour[i];
   } // directions
   for (int i = 0; i < 3; i++) {
      in_file >> BlkIO.direction_neighbour_to_me[i];
   } // output directions
   in_file.unsetf(ios::scientific);
   return (in_file);
}

//=============================================================
class Block_Boundary_Elements_on_Domain_Extent {
  public:
    int on_grid_boundary[MAX_BOUNDARY_ELEMENTS_FOR_A_BLOCK]; 
   
    // constructor
    Block_Boundary_Elements_on_Domain_Extent(void) {
      for (int elem = 0; elem < MAX_BOUNDARY_ELEMENTS_FOR_A_BLOCK; elem++) {
         on_grid_boundary[elem] = 0;
         // 0 means this boundary element is not on the domain extent.
      } /* endfor */
    }
   
    // use default destructor 
   
    void broadcast(void);
   
    // assignment operator 
    Block_Boundary_Elements_on_Domain_Extent &operator =(const Block_Boundary_Elements_on_Domain_Extent &BEs); 
   
    /* Input-output operators. */

    friend ostream& operator << (ostream &out_file, 
                                 const Block_Boundary_Elements_on_Domain_Extent &BEs);
   
    friend istream& operator >> (istream &in_file,  
                                 const Block_Boundary_Elements_on_Domain_Extent &BEs);
};

//----------------- Assignment ----------------------------//
inline Block_Boundary_Elements_on_Domain_Extent &Block_Boundary_Elements_on_Domain_Extent::
operator =(const Block_Boundary_Elements_on_Domain_Extent &BEs) {
   //self assignment protection
   if (this != &BEs) {   
      //copy assignment
      for (int i = 0; i < MAX_BOUNDARY_ELEMENTS_FOR_A_BLOCK; i++) {
         on_grid_boundary[i] = BEs.on_grid_boundary[i];
      } // endfor boundary elements on domain extent   
   } /* endif */
   return (*this);
}

/**********************************************************************************
 * Block_Boundary_Elements_on_Domain_Extent -- Input-output operators.            *
 **********************************************************************************/
inline ostream &operator << (ostream &out_file, 
                             const Block_Boundary_Elements_on_Domain_Extent &BEs) {
   out_file.precision(10);
   out_file.setf(ios::scientific);
   for (int i = 0; i < MAX_BOUNDARY_ELEMENTS_FOR_A_BLOCK; i++) {
      out_file << " " << BEs.on_grid_boundary[i];
   } /* endfor */
   out_file.unsetf(ios::scientific);
   return (out_file);
}

inline istream &operator >> (istream &in_file, 
                             Block_Boundary_Elements_on_Domain_Extent &BEs) {
   in_file.precision(10);
   in_file.setf(ios::scientific);
   for (int i = 0; i < MAX_BOUNDARY_ELEMENTS_FOR_A_BLOCK; i++) {
      in_file >> BEs.on_grid_boundary[i];
   } /* endfor */
   in_file.unsetf(ios::scientific);
   return (in_file);
}

#endif // _BLOCK_ORIENTATION_INFO_INCLUDED

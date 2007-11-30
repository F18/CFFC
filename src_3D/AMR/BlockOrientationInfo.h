
/* BlockOrientationInfo.h:  Header file defining block orientation information class. */

#ifndef _BLOCK_ORIENTATION_INFO_INCLUDED
#define _BLOCK_ORIENTATION_INFO_INCLUDED


#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>

using namespace std;

#include "coord_transform.h"
#include "block_connectivity.h"

#ifndef _MPI_INCLUDED
#include "../MPI/MPI.h"
#endif // _MPI_INCLUDED

#define BLOCK_ORIENTATION_MAX_NEIGHBOUR            2
#define MAX_BOUNDARY_ELEMENTS_FOR_A_BLOCK          27

/* Define the 3D hexahedral grid multiblock classes. */
class  Block_Orientation_Info{
  public:
   
   int                          ctm_offsets[6]; // compact transformation matrix (3 components) and offsets (3 elements)
   int                          direction_me_to_neighbour[3]; 
   int                          direction_neighbour_to_me[3]; 
                             
    
   // constructor
   Block_Orientation_Info(void){
      ctm_offsets[0] = 1 ; ctm_offsets[1] = 1 ;  ctm_offsets[2] = 1 ;
      ctm_offsets[3] = 0 ; ctm_offsets[4] = 0 ;  ctm_offsets[5] = 0  ;  
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
                                   const int blockIndex, const int my_i, 
                                   const int my_j, const int my_k, 
                                   int &iblk_neigh);
   int compute_message_tag(const int i_index, const int j_index,  const int k_index);
   int convert_boundary_elements_from_ijk_to_orientations(const int my_i, const int my_j, const int my_k);
   void my_index(int &i_index, int &j_index, int &k_index);
   void neighbour_index(int &i_index, int &j_index, int &k_index);
   void broadcast(void);
   
   // assignment operator 
   Block_Orientation_Info  & operator =(const Block_Orientation_Info &BlkOI); 
   
   /* Input-output operators. */
   friend ostream& operator << (ostream &out_file, 
                                const Block_Orientation_Info &BlkOI);
   
   friend istream& operator >> (istream &in_file,  
                                Block_Orientation_Info &BlkOI);
   
   
   
};


//----------------- Assignment ----------------------------//
inline Block_Orientation_Info & Block_Orientation_Info::operator =(const Block_Orientation_Info &BlkOI){
   
   //self assignment protection
   if( this != &BlkOI){   
      //copy assignment
      for(int i=0; i<6; i++){
         ctm_offsets[i] = BlkOI.ctm_offsets[i];
      } // ctm and offsets   
      
      for(int i=0; i<3; i++){
         direction_me_to_neighbour[i] = BlkOI.direction_me_to_neighbour[i];
         direction_neighbour_to_me[i] = BlkOI.direction_neighbour_to_me[i];
         // directions          
      }
      
   }   
   else {
      cerr<<"\n error in assignment operator for Block_Orientation_Info class \n ";
   }  
   
   return (*this);

}

   
class  Block_Boundary_Elements_on_Domain_Extent{
  public:
   
   int         boundary_element_on_domain_extent[ MAX_BOUNDARY_ELEMENTS_FOR_A_BLOCK]; 
   
   
   // constructor
   Block_Boundary_Elements_on_Domain_Extent(void){
      
      for(int elem = 0; elem<MAX_BOUNDARY_ELEMENTS_FOR_A_BLOCK; elem++){
         
         boundary_element_on_domain_extent[elem] = 0;
         // 0 means this boundary element is not on the domain extent.
         
      }
      
   }
   
   // use default destructor 
   
   void Block_Boundary_Elements_on_Domain_Extent::broadcast(void);
   
   // assignment operator 
   Block_Boundary_Elements_on_Domain_Extent  & operator =(const Block_Boundary_Elements_on_Domain_Extent  
                                                          &be_on_domain_extent); 
   
   /* Input-output operators. */
   friend ostream& operator << (ostream &out_file, 
                                const  Block_Boundary_Elements_on_Domain_Extent &be_on_domain_extent);
   
   friend istream& operator >> (istream &in_file,  
                                const Block_Boundary_Elements_on_Domain_Extent &be_on_domain_extent);
   
   
   
};


//----------------- Assignment ----------------------------//
inline Block_Boundary_Elements_on_Domain_Extent & Block_Boundary_Elements_on_Domain_Extent::operator =(
   const Block_Boundary_Elements_on_Domain_Extent &be_on_domain_extent){
   
   //self assignment protection
   if( this != &be_on_domain_extent){   
      //copy assignment
      for(int i=0; i<MAX_BOUNDARY_ELEMENTS_FOR_A_BLOCK; i++){
         boundary_element_on_domain_extent[i] = be_on_domain_extent.boundary_element_on_domain_extent[i];
      } // boundary elements on domain extent   
      
   }else {
      cerr<<"\n error in assignment operator for Boundary_Elements_on_Domain_Extent class \n ";
   }  
   
   return (*this);

}

 

/*******************************************************************************
 *
 * Special types
 *
 ******************************************************************************/
namespace BE
{
   enum Boundary_Elements{
      
      BSW = 0,
      SW = 1,
      TSW = 2,
      BW = 3,
      W =4,
      TW = 5,
      BNW = 6,
      NW = 7,
      TNW = 8,
      BS = 9,
      S = 10,
      TS = 11,
      B = 12,
      T = 14,
      BN = 15,
      N = 16,
      TN = 17,
      BSE = 18,
      SE = 19,
      TSE = 20,
      BE = 21,
      E = 22,
      TE = 23,
      BNE = 24,
      NE = 25,
      TNE = 26
   
   };
}  

#endif // _BLOCK_ORIENTATION_INFO_INCLUDED

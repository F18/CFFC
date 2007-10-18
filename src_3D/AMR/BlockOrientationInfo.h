
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


#define BLOCK_ORIENTATION_MAX_NEIGHBOUR            2

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
                                   const int my_j, const my_k, int &iblk_neigh);
   int compute_message_tag(const int i_index, const int j_index,  const int k_index);
   int convert_boundary_elements_from_ijk_to_orientations(const int my_i, const int my_j, const int my_k);
   void my_index(int &i_index, int &j_index, int &k_index);
   void neighbour_index(int &i_index, int &j_index, int &k_index);
   
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

   
 


#endif // _BLOCK_ORIENTATION_INFO_INCLUDED

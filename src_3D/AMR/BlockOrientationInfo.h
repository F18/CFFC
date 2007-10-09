
/* BlockOrientationInfo.h:  Header file defining block orientation information class. */

#ifndef _BLOCK_ORIENTATION_INFO_INCLUDED
#define _BLOCK_ORIENTATION_INFO_INCLUDED


#include "coord_transform.h"
#include "block_connectivity.h"


#define BLOCK_ORIENTATION_MAX_NEIGHBOUR            2

/* Define the 3D hexahedral grid multiblock classes. */
class  Block_Orientation_Info{
  public:
   
   int                          ctm_offsets[6]; // compact transformation matrix (3 components) and offsets (3 elements)
   int                          direction_to_neighbour[3]; 
                             
    
   // constructor
   Block_Orientation_Info(void){
      ctm_offsets[0] = 1 ; ctm_offsets[1] = 1 ;  ctm_offsets[2] = 1 ;
      ctm_offsets[3] = 0 ; ctm_offsets[4] = 0 ;  ctm_offsets[5] = 0  ;  
      direction_to_neighbour[0] = 0 ; 
      direction_to_neighbour[1] = 0 ;
      direction_to_neighbour[2] = 0 ; 
   }
   
// use default destructor 

// some member functnions ... .... 

// assignment operator 
 




   
};

#endif // _BLOCK_ORIENTATION_INFO_INCLUDED

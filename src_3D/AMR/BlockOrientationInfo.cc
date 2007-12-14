/* BlockOrientationInfo.cc:  */

#ifndef _BLOCK_ORIENTATION_INFO_INCLUDED
#include"BlockOrientationInfo.h"
#endif// _BLOCK_ORIENTATION_INFO_INCLUDED

#ifndef  _SYSTEM_LINUX_INCLUDED
#include "../System/System_Linux.h"
#endif //_SYSTEM_LINUX_INCLUDED

void Block_Orientation_Info::set_block_orientation_info(BlkC::BlockConnectivity &blkc, 
                                                        const int iblk,
                                                        const int which_neighbour,
                                                        const int my_i, 
                                                        const int my_j, 
                                                        const int my_k, 
                                                        int &iblk_neigh) {
   
   direction_me_to_neighbour[0] = my_i;
   direction_me_to_neighbour[1] = my_j;
   direction_me_to_neighbour[2] = my_k;
   
   blkc.neighbour_block(iblk,  
                        direction_me_to_neighbour[0],
                        direction_me_to_neighbour[1],
                        direction_me_to_neighbour[2],
                        which_neighbour,
                        iblk_neigh, 
                        direction_neighbour_to_me[0],
                        direction_neighbour_to_me[1],
                        direction_neighbour_to_me[2],
                        ctm_offsets);

}

int Block_Orientation_Info::convert_boundary_elements_from_ijk_to_orientations(const int my_i, 
                                                                               const int my_j, 
                                                                               const int my_k) {
   
   char orientation[4];

   strcpy(orientation, "");
   
   if (my_k == BlkC::KMax) {
      strcat(orientation, "T");
   } else if (my_k == BlkC::KMin) {
      strcat(orientation, "B");
   } /* endif */

   if (my_j == BlkC::JMax) {
      strcat(orientation, "N");
   } else if (my_j == BlkC::JMin) {
      strcat(orientation, "S");
   } /* endif */

   if (my_i == BlkC::IMax) {
      strcat(orientation, "E");
   } else if(my_i == BlkC::IMin) {
      strcat(orientation, "W");
   } /* endif */
   
   cout << "\n orientation "<< orientation << endl;
   
   return 0;

}

/*==============================================================================
 * The index (send to me from my neighbour) being transformed to the index in
   my block.
 *============================================================================*/
void Block_Orientation_Info::my_index(int &i_index, 
                                      int &j_index,  
                                      int &k_index) {
   
   CoordTransform transformationMatrix(ctm_offsets); 
   transformationMatrix.transform(i_index, j_index, k_index);
   
}

/*==============================================================================
 * The index from my block is being transformed to the index in
   my neighbour's block.
 *============================================================================*/
void Block_Orientation_Info::neighbour_index(int &i_index, 
                                             int &j_index,  
                                             int &k_index) {
     
   CoordTransform transformationMatrix(ctm_offsets);
   transformationMatrix.reverse();
   transformationMatrix.transform(i_index, j_index, k_index);
   
}

int Block_Orientation_Info::compute_message_tag(const int i_index, 
                                                const int j_index,  
                                                const int k_index) {
   
   int tag = 0;
   return  tag = (i_index+1)*9 + (j_index+1)*3 + (k_index+1) ;
   
}

void Block_Orientation_Info::broadcast(void){

#ifdef _MPI_VERSION
   MPI::COMM_WORLD.Bcast(ctm_offsets, 6, MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(direction_me_to_neighbour, 3, MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(direction_neighbour_to_me, 3, MPI::INT, 0);
#endif    

}

void  Block_Boundary_Elements_on_Domain_Extent::broadcast(void){

#ifdef _MPI_VERSION
   MPI::COMM_WORLD.Bcast(on_grid_boundary, MAX_BOUNDARY_ELEMENTS_FOR_A_BLOCK, MPI::INT, 0);
#endif    

}

   

/* BlockOrientationInfo.cc:  */

#ifndef _BLOCK_ORIENTATION_INFO_INCLUDED
#include"BlockOrientationInfo.h"
#endif// _BLOCK_ORIENTATION_INFO_INCLUDED

/*************************************************************************************************************************
 * Direction_Indices::direction_number_to_direction_indices -- Return indices given a direction number.                  *
 *************************************************************************************************************************/
Direction_Indices Direction_Indices::direction_number_to_direction_indices(const int direction_number) const {
   switch (direction_number) { 
   case DIRECTION::BSW : // 0
     return Direction_Indices(-1, -1, -1);    
   case DIRECTION::SW :  // 1
     return Direction_Indices(-1, -1, 0);
   case DIRECTION::TSW : // 2
     return Direction_Indices(-1, -1, 1);
   case DIRECTION::BW :  // 3
     return Direction_Indices(-1, 0, -1);    
   case DIRECTION::W :   // 4
     return Direction_Indices(-1, 0, 0);
   case DIRECTION::TW :  // 5
     return Direction_Indices(-1, 0, 1);
   case DIRECTION::BNW : // 6
     return Direction_Indices(-1, 1, -1);    
   case DIRECTION::NW :  // 7
     return Direction_Indices(-1, 1, 0);
   case DIRECTION::TNW : // 8
     return Direction_Indices(-1, 1, 1);
   case DIRECTION::BS :  // 9
     return Direction_Indices(0, -1, -1);    
   case DIRECTION::S :   // 10
     return Direction_Indices(0, -1, 0);
   case DIRECTION::TS :  // 11
     return Direction_Indices(0, -1, 1);
   case DIRECTION::B :   // 12
     return Direction_Indices(0, 0, -1);    
   case DIRECTION::T :   // 13
     return Direction_Indices(0, 0, 1);
   case DIRECTION::BN :  // 14
     return Direction_Indices(0, 1, -1);    
   case DIRECTION::N :   // 15
     return Direction_Indices(0, 1, 0);
   case DIRECTION::TN :  // 16
     return Direction_Indices(0, 1, 1);
   case DIRECTION::BSE : // 17
     return Direction_Indices(1, -1, -1);    
   case DIRECTION::SE :  // 18
     return Direction_Indices(1, -1, 0);
   case DIRECTION::TSE : // 19
     return Direction_Indices(1, -1, 1);
   case DIRECTION::BE :  // 20
     return Direction_Indices(1, 0, -1);    
   case DIRECTION::E :   // 21
     return Direction_Indices(1, 0, 0);
   case DIRECTION::TE :  // 22
     return Direction_Indices(1, 0, 1);
   case DIRECTION::BNE : // 23
     return Direction_Indices(1, 1, -1);    
   case DIRECTION::NE :  // 24
     return Direction_Indices(1, 1, 0);
   case DIRECTION::TNE : // 25
     return Direction_Indices(1, 1, 1);
   default :
     return Direction_Indices(0, 0, 0);
   }; /* endswitch */
}

/*************************************************************************************************************************
 * Direction_Indices::direction_indices_to_direction_number -- Return direction number given the indices.                *
 *************************************************************************************************************************/
int Direction_Indices::direction_indices_to_direction_number(void) const {
   int dir_num = direction_indices_to_boundary_element_number();
   if (dir_num > BE::ME) dir_num = dir_num - 1;
   return dir_num;
}

int Direction_Indices::direction_indices_to_direction_number(const Direction_Indices &DI) const {
   int dir_num = DI.direction_indices_to_boundary_element_number();
   if (dir_num > BE::ME) dir_num = dir_num - 1;
   return dir_num;
}

/*************************************************************************************************************************
 * Direction_Indices::boundary_element_number_to_direction_indices -- Return boundary element number given a             *
 *                                                                    direction number.                                  *
 *************************************************************************************************************************/
Direction_Indices Direction_Indices::boundary_element_number_to_direction_indices(const int boundary_element_number) const {
   switch (boundary_element_number) { 
   case BE::BSW : // 0
     return Direction_Indices(-1, -1, -1);    
   case BE::SW :  // 1
     return Direction_Indices(-1, -1, 0);
   case BE::TSW : // 2
     return Direction_Indices(-1, -1, 1);
   case BE::BW :  // 3
     return Direction_Indices(-1, 0, -1);    
   case BE::W :   // 4
     return Direction_Indices(-1, 0, 0);
   case BE::TW :  // 5
     return Direction_Indices(-1, 0, 1);
   case BE::BNW : // 6
     return Direction_Indices(-1, 1, -1);    
   case BE::NW :  // 7
     return Direction_Indices(-1, 1, 0);
   case BE::TNW : // 8
     return Direction_Indices(-1, 1, 1);
   case BE::BS :  // 9
     return Direction_Indices(0, -1, -1);    
   case BE::S :   // 10
     return Direction_Indices(0, -1, 0);
   case BE::TS :  // 11
     return Direction_Indices(0, -1, 1);
   case BE::B :   // 12
     return Direction_Indices(0, 0, -1);    
   case BE::ME :  // 13
     return Direction_Indices(0, 0, 0);
   case BE::T :   // 14
     return Direction_Indices(0, 0, 1);
   case BE::BN :  // 15
     return Direction_Indices(0, 1, -1);    
   case BE::N :   // 16
     return Direction_Indices(0, 1, 0);
   case BE::TN :  // 17
     return Direction_Indices(0, 1, 1);
   case BE::BSE : // 18
     return Direction_Indices(1, -1, -1);    
   case BE::SE :  // 19
     return Direction_Indices(1, -1, 0);
   case BE::TSE : // 20
     return Direction_Indices(1, -1, 1);
   case BE::BE :  // 21
     return Direction_Indices(1, 0, -1);    
   case BE::E :   // 22
     return Direction_Indices(1, 0, 0);
   case BE::TE :  // 23
     return Direction_Indices(1, 0, 1);
   case BE::BNE : // 24
     return Direction_Indices(1, 1, -1);    
   case BE::NE :  // 25
     return Direction_Indices(1, 1, 0);
   case BE::TNE : // 26
     return Direction_Indices(1, 1, 1);
   default :
     return Direction_Indices(0, 0, 0);
   }; /* endswitch */
}

/*************************************************************************************************************************
 * Direction_Indices::direction_indices_to_boundary_element_number --  Return boundary element number given the indices. *
 *************************************************************************************************************************/
int Direction_Indices::direction_indices_to_boundary_element_number(void) const {
  return ((i+1)*9 + (j+1)*3 + (k+1));
}

int Direction_Indices::direction_indices_to_boundary_element_number(const Direction_Indices &DI) const {
  return ((DI.i+1)*9 + (DI.j+1)*3 + (DI.k+1));
}

/*************************************************************************************************************************
 * Direction_Indices::reference_geomety_node -- Return  indices of node for referencing the geometry.                    *
 *************************************************************************************************************************/
Direction_Indices Direction_Indices::reference_geomety_node(const int boundary_element_number,
                                                            const int i_low,
                                                            const int i_up,
                                                            const int j_low,
                                                            const int j_up,
                                                            const int k_low,
                                                            const int k_up,
                                                            const int nghost) const {
   switch (boundary_element_number) { 
   case BE::BSW : // 0
     return Direction_Indices(i_low, j_low, k_low);    
   case BE::SW :  // 1
     return Direction_Indices(i_low, j_low, k_low);  // edge
   case BE::TSW : // 2
     return Direction_Indices(i_low, j_low, k_up);
   case BE::BW :  // 3
     return Direction_Indices(i_low, j_low, k_low); // edge
   case BE::W :   // 4
     return Direction_Indices(-1, 0, 0);
   case BE::TW :  // 5
     return Direction_Indices(-1, 0, 1);
   case BE::BNW : // 6
     return Direction_Indices(i_low, j_up, k_low);    
   case BE::NW :  // 7
     return Direction_Indices(-1, 1, 0);
   case BE::TNW : // 8
     return Direction_Indices(i_low, j_up, k_up);
   case BE::BS :  // 9
     return Direction_Indices(0, -1, -1);    
   case BE::S :   // 10
     return Direction_Indices(0, -1, 0);
   case BE::TS :  // 11
     return Direction_Indices(0, -1, 1);
   case BE::B :   // 12
     return Direction_Indices(0, 0, -1);    
   case BE::ME :  // 13
     return Direction_Indices(0, 0, 0);
   case BE::T :   // 14
     return Direction_Indices(0, 0, 1);
   case BE::BN :  // 15
     return Direction_Indices(0, 1, -1);    
   case BE::N :   // 16
     return Direction_Indices(0, 1, 0);
   case BE::TN :  // 17
     return Direction_Indices(0, 1, 1);
   case BE::BSE : // 18
     return Direction_Indices(i_up, j_low, k_low);    
   case BE::SE :  // 19
     return Direction_Indices(1, -1, 0);
   case BE::TSE : // 20
     return Direction_Indices(i_up, j_low, k_up);
   case BE::BE :  // 21
     return Direction_Indices(1, 0, -1);    
   case BE::E :   // 22
     return Direction_Indices(1, 0, 0);
   case BE::TE :  // 23
     return Direction_Indices(1, 0, 1);
   case BE::BNE : // 24
     return Direction_Indices(i_up, j_up, k_low);    
   case BE::NE :  // 25
     return Direction_Indices(1, 1, 0);
   case BE::TNE : // 26
     return Direction_Indices(i_up, j_up, k_up);
   default :
     return Direction_Indices(0, 0, 0);
   }; /* endswitch */
}

//=======================================================
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

   

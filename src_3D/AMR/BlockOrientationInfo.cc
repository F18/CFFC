/* BlockOrientationInfo.cc:  */

#ifndef _BLOCK_ORIENTATION_INFO_INCLUDED
#include"BlockOrientationInfo.h"
#endif// _BLOCK_ORIENTATION_INFO_INCLUDED


ostream &operator << (ostream &out_file, 
                      const Block_Orientation_Info &BlkIO) {
   out_file.precision(10);
   out_file.setf(ios::scientific);
   for( int i=0; i<6; i++){
      out_file<<" "<<BlkIO.ctm_offsets[i];
   }// output ctm and offsets
   
   for( int i=0; i<3; i++){
      out_file<<" "<<BlkIO.direction_me_to_neighbour[i];
   }
   for( int i=0; i<3; i++){
      out_file<<" "<<BlkIO.direction_neighbour_to_me[i];
   } // output directions
   
   out_file.unsetf(ios::scientific);
   return (out_file);
}



istream &operator >> (istream &in_file, 
                      Block_Orientation_Info &BlkIO) {
   in_file.precision(10);
   in_file.setf(ios::scientific);
   for( int i=0; i<6; i++){
      in_file>>BlkIO.ctm_offsets[i];
   }// output ctm and offsets
   
   for( int i=0; i<3; i++){
      in_file>>BlkIO.direction_me_to_neighbour[i];
   }
   for( int i=0; i<3; i++){
      in_file>>BlkIO.direction_neighbour_to_me[i];
   } // output directions
   
   in_file.unsetf(ios::scientific);
   return (in_file);
}


void  Block_Orientation_Info::set_block_orientation_info(
   BlkC::BlockConnectivity &blkc, const int iblk, 
   const int my_i, const int my_j, const my_k, int &iblk_neigh){
   
   direction_me_to_neighbour[0] = my_i;
   direction_me_to_neighbour[1] = my_j;
   direction_me_to_neighbour[2] = my_k;
   
   blkc.neighbour_block(iblk,  
                        direction_me_to_neighbour[0],
                        direction_me_to_neighbour[1],
                        direction_me_to_neighbour[2],
                        0,
                        iblk_neigh, 
                        direction_neighbour_to_me[0],
                        direction_neighbour_to_me[1],
                        direction_neighbour_to_me[2],
                        ctm_offsets);
}

int Block_Orientation_Info::convert_boundary_elements_from_ijk_to_orientations(
   const int my_i, const int my_j, const int my_k){
   
   char orientation[] = "neigh";
   
   if(my_k == BlkC::KMax){
      strcat(orientation, "T");
   }else if (my_k == BlkC::KMin){
      strcat(orientation, "B");
   }

   if(my_j == BlkC::JMax){
      strcat(orientation, "N");
   }else if(my_j == BlkC::JMin){
      strcat(orientation, "S");
   }

   if(my_i == BlkC::IMax){
      strcat(orientation, "E");
   }else if(my_i == BlkC::IMin){
      strcat(orientation, "W");
   }
   
   cout<<"\n orientation "<<orientation<<endl;
   
   return 0;
   
    
}

void  Block_Orientation_Info::my_index(int &i_index, int &j_index,  int &k_index){
   
   CoordTransform transformationMatrix(ctm_offsets); 
   
   transformationMatrix.transpose();
   
   transformationMatrix.transform(i_index, j_index, k_index);
   
   
}

void Block_Orientation_Info::neighbour_index(
   int &i_index, int &j_index,  int &k_index){
   
   CoordTransform transformationMatrix(ctm_offsets);
   transformationMatrix.transform(i_index, j_index, k_index);
     
}




int Block_Orientation_Info::compute_message_tag(const int i_index, const int j_index,  const int k_index){
   
   int tag = 0;
   return  tag = (i_index+1)*9 + (j_index+1)*3 + (k_index+1) +1 ;
   
   
}


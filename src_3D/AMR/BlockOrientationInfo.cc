/* BlockOrientationInfo.cc:  */

#ifndef _BLOCK_ORIENTATION_INFO_INCLUDED
#include"BlockOrientationInfo.h"
#endif// _BLOCK_ORIENTATION_INFO_INCLUDED

#ifndef  _SYSTEM_LINUX_INCLUDED
#include "../System/System_Linux.h"
#endif //_SYSTEM_LINUX_INCLUDED

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
   const int my_i, const int my_j, const int my_k, 
   int &iblk_neigh){
   
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
   
   char orientation[4];

   strcpy(orientation, "");
   
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


/*==============================================================================
 * The index (send to me from my neighbour) being transformed to the index in
   my block.
 *============================================================================*/
void  Block_Orientation_Info::my_index(int &i_index, int &j_index,  int &k_index){
   
   CoordTransform transformationMatrix(ctm_offsets); 
   transformationMatrix.transform(i_index, j_index, k_index);
   
}
/*==============================================================================
 * The index from my block is being transformed to the index in
   my neighbour's block.
 *============================================================================*/
void Block_Orientation_Info::neighbour_index(
   int &i_index, int &j_index,  int &k_index){
     
   for ( int iProc = 0; iProc !=  CFFC_MPI::Number_of_Processors; ++iProc ) {
      if (  CFFC_MPI::This_Processor_Number == iProc ) {
         cout<<"\n CFFC_MPI::This_Processor_Number = "<< CFFC_MPI::This_Processor_Number<<endl;
         cout<<"\n ctm "<<ctm_offsets[0]<<"  "<<ctm_offsets[1]<<"  "<<ctm_offsets[2]<<"  "
             <<ctm_offsets[3]<<"  "<<ctm_offsets[4]<<"  "<<ctm_offsets[5]<<endl;
         
         
         System::sleep(0.1);
      }
      MPI::COMM_WORLD.Barrier();
   }
                    

   CoordTransform transformationMatrix(ctm_offsets);
   transformationMatrix.reverse();
   transformationMatrix.transform(i_index, j_index, k_index);
   
}

int Block_Orientation_Info::compute_message_tag(const int i_index, 
                                                const int j_index,  
                                                const int k_index){
   
   int tag = 0;
   
   return  tag = (i_index+1)*9 + (j_index+1)*3 + (k_index+1) ;
   
   
}

void Block_Orientation_Info::broadcast(void){
#ifdef _MPI_VERSION
   MPI::COMM_WORLD.Bcast(&ctm_offsets[0], 1, MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(&ctm_offsets[1], 1, MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(&ctm_offsets[2], 1, MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(&ctm_offsets[3], 1, MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(&ctm_offsets[4], 1, MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(&ctm_offsets[5], 1, MPI::INT, 0);

   MPI::COMM_WORLD.Bcast(&direction_me_to_neighbour[0], 1, MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(&direction_me_to_neighbour[1], 1, MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(&direction_me_to_neighbour[2], 1, MPI::INT, 0);
   
   MPI::COMM_WORLD.Bcast(&direction_neighbour_to_me[0], 1, MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(&direction_neighbour_to_me[1], 1, MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(&direction_neighbour_to_me[2], 1, MPI::INT, 0);

#endif    
}

   

ostream &operator << (ostream &out_file, 
                      const Block_Boundary_Elements_on_Domain_Extent &be_on_domain_extent)
{
   
   out_file.precision(10);
   out_file.setf(ios::scientific);
   for( int i=0; i<27; i++){
      out_file<<" "<<be_on_domain_extent.boundary_element_on_domain_extent[i];
   }// output ctm and offsets
   
   out_file.unsetf(ios::scientific);
   return (out_file);
}



istream &operator >> (istream &in_file, 
                      Block_Boundary_Elements_on_Domain_Extent &be_on_domain_extent) {
   in_file.precision(10);
   in_file.setf(ios::scientific);
   for( int i=0; i<27; i++){
      in_file>>be_on_domain_extent.boundary_element_on_domain_extent[i];
   }// output ctm and offsets
  
   in_file.unsetf(ios::scientific);
   return (in_file);
}


void  Block_Boundary_Elements_on_Domain_Extent::broadcast(void){
#ifdef _MPI_VERSION
   MPI::COMM_WORLD.Bcast(& boundary_element_on_domain_extent[0], 1, MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(& boundary_element_on_domain_extent[1], 1, MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(& boundary_element_on_domain_extent[2], 1, MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(& boundary_element_on_domain_extent[3], 1, MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(& boundary_element_on_domain_extent[4], 1, MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(& boundary_element_on_domain_extent[5], 1, MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(& boundary_element_on_domain_extent[6], 1, MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(& boundary_element_on_domain_extent[7], 1, MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(& boundary_element_on_domain_extent[8], 1, MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(& boundary_element_on_domain_extent[9], 1, MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(& boundary_element_on_domain_extent[10], 1, MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(& boundary_element_on_domain_extent[11], 1, MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(& boundary_element_on_domain_extent[12], 1, MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(& boundary_element_on_domain_extent[13], 1, MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(& boundary_element_on_domain_extent[14], 1, MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(& boundary_element_on_domain_extent[15], 1, MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(& boundary_element_on_domain_extent[16], 1, MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(& boundary_element_on_domain_extent[17], 1, MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(& boundary_element_on_domain_extent[18], 1, MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(& boundary_element_on_domain_extent[19], 1, MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(& boundary_element_on_domain_extent[20], 1, MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(& boundary_element_on_domain_extent[21], 1, MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(& boundary_element_on_domain_extent[22], 1, MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(& boundary_element_on_domain_extent[23], 1, MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(& boundary_element_on_domain_extent[24], 1, MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(& boundary_element_on_domain_extent[25], 1, MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(& boundary_element_on_domain_extent[26], 1, MPI::INT, 0);

#endif    

}

   

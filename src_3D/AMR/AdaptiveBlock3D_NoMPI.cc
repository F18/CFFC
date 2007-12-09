/* AdaptiveBlock3D_NoMPI.cc:  Subroutines for adaptive blocks classes
                              on a single CPU architecture with no 
                              access to MPI. */

/* Include adaptive block header file. */

#ifndef _ADAPTIVBLOCK3D_INCLUDED
#include "AdaptiveBlock3D.h"
#endif // _ADAPTIVEBLOCK3D_INCLUDED

/**********************************************************
 * Routine: Exchange_Messages_NoResChange                 *
 *                                                        *
 * Sends and receives solution information contained in   *
 * the preloaded message passing buffers between          *
 * neighbouring adaptive blocks with no mesh resolution   *
 * changes.                                               *
 *                                                        *
 **********************************************************/
int AdaptiveBlock3D_List::Exchange_Messages_NoResChange(AdaptiveBlock3D_List &Blk_List,
                                                        const int Number_of_Solution_Variables) {

    int i_blk, neighbour_cpu, neighbour_blk, buffer_size, neighbour_buffer_size;

    int i_bound_elem, neighbour_bound_elem; // index for boundary element, face edge or vertex

    int number_neighbours[MAX_BOUNDARY_ELEMENTS_FOR_A_BLOCK];
    AdaptiveBlock3D_Info neighbour_info[MAX_BOUNDARY_ELEMENTS_FOR_A_BLOCK];
    
    /* Perform message passing for the block faces and corners of each 
       solution block having no mesh resolution change. */

    for (i_blk = 0 ; i_blk <= Blk_List.Nblk-1 ; ++i_blk) {
      if (Blk_List.Block[i_blk].used) {
         // Assign the boundary element information
         number_neighbours[BE::BSW] = Blk_List.Block[i_blk].nBSW;
         neighbour_info[BE::BSW] = Blk_List.Block[i_blk].infoBSW[0];
      
         number_neighbours[BE::SW] = Blk_List.Block[i_blk].nSW;
         neighbour_info[BE::SW] = Blk_List.Block[i_blk].infoSW[0];
      
         number_neighbours[BE::TSW] = Blk_List.Block[i_blk].nTSW;
         neighbour_info[BE::TSW] = Blk_List.Block[i_blk].infoTSW[0];
      
         number_neighbours[BE::BW] = Blk_List.Block[i_blk].nBW;
         neighbour_info[BE::BW] = Blk_List.Block[i_blk].infoBW[0];
      
         number_neighbours[BE::W] = Blk_List.Block[i_blk].nW;
         neighbour_info[BE::W] = Blk_List.Block[i_blk].infoW[0];
      
         number_neighbours[BE::TW] = Blk_List.Block[i_blk].nTW;
         neighbour_info[BE::TW] = Blk_List.Block[i_blk].infoTW[0];
      
         number_neighbours[BE::BNW] = Blk_List.Block[i_blk].nBNW;
         neighbour_info[BE::BNW] = Blk_List.Block[i_blk].infoBNW[0];
      
         number_neighbours[BE::NW] = Blk_List.Block[i_blk].nNW;
         neighbour_info[BE::NW] = Blk_List.Block[i_blk].infoNW[0];
      
         number_neighbours[BE::TNW] = Blk_List.Block[i_blk].nTNW;
         neighbour_info[BE::TNW] = Blk_List.Block[i_blk].infoTNW[0];
      
         number_neighbours[BE::BS] = Blk_List.Block[i_blk].nBS;
         neighbour_info[BE::BS] = Blk_List.Block[i_blk].infoBS[0];
      
         number_neighbours[BE::S] = Blk_List.Block[i_blk].nS;
         neighbour_info[BE::S] = Blk_List.Block[i_blk].infoS[0];
      
         number_neighbours[BE::TS] = Blk_List.Block[i_blk].nTS;
         neighbour_info[BE::TS] = Blk_List.Block[i_blk].infoTS[0];
      
         number_neighbours[BE::B] = Blk_List.Block[i_blk].nB;
         neighbour_info[BE::B] = Blk_List.Block[i_blk].infoB[0];
      
         number_neighbours[BE::T] = Blk_List.Block[i_blk].nT;
         neighbour_info[BE::T] = Blk_List.Block[i_blk].infoT[0];
      
         number_neighbours[BE::BN] = Blk_List.Block[i_blk].nBN;
         neighbour_info[BE::BN] = Blk_List.Block[i_blk].infoBN[0];
      
         number_neighbours[BE::N] = Blk_List.Block[i_blk].nN;
         neighbour_info[BE::N] = Blk_List.Block[i_blk].infoN[0];
      
         number_neighbours[BE::TN] = Blk_List.Block[i_blk].nTN;
         neighbour_info[BE::TN] = Blk_List.Block[i_blk].infoTN[0];
      
         number_neighbours[BE::BSE] = Blk_List.Block[i_blk].nBSE;
         neighbour_info[BE::BSE] = Blk_List.Block[i_blk].infoBSE[0];
      
         number_neighbours[BE::SE] = Blk_List.Block[i_blk].nSE;
         neighbour_info[BE::SE] = Blk_List.Block[i_blk].infoSE[0];
      
         number_neighbours[BE::TSE] = Blk_List.Block[i_blk].nTSE;
         neighbour_info[BE::TSE] = Blk_List.Block[i_blk].infoTSE[0];
      
         number_neighbours[BE::BE] = Blk_List.Block[i_blk].nBE;
         neighbour_info[BE::BE] = Blk_List.Block[i_blk].infoBE[0];
      
         number_neighbours[BE::E] = Blk_List.Block[i_blk].nE;
         neighbour_info[BE::E] = Blk_List.Block[i_blk].infoE[0];
      
         number_neighbours[BE::TE] = Blk_List.Block[i_blk].nTE;
         neighbour_info[BE::TE] = Blk_List.Block[i_blk].infoTE[0]; 
      
         number_neighbours[BE::BNE] = Blk_List.Block[i_blk].nBNE;
         neighbour_info[BE::BNE] = Blk_List.Block[i_blk].infoBNE[0]; 
       
         number_neighbours[BE::NE] = Blk_List.Block[i_blk].nNE;
         neighbour_info[BE::NE] = Blk_List.Block[i_blk].infoNE[0];  
      
         number_neighbours[BE::TNE] = Blk_List.Block[i_blk].nTNE;
         neighbour_info[BE::TNE] = Blk_List.Block[i_blk].infoTNE[0]; 

         // Perform direct copies of messages to appropriate message buffers.
         for (int ii = -1; ii <= 1; ii++) {
            for (int jj = -1; jj <= 1; jj++) {
               for (int kk = -1; kk <= 1; kk++) {
                  i_bound_elem = 9*(ii+1) + 3*(jj+1) + (kk+1);
               
                  if ((number_neighbours[i_bound_elem] == 1) && 
                      (i_bound_elem != BE::ME) &&
                      (Blk_List.Block[i_blk].info.level == neighbour_info[i_bound_elem].level)) {

                     neighbour_blk = neighbour_info[i_bound_elem].blknum;
                     
                     buffer_size = ((abs(ii)*Blk_List.Block[i_blk].info.dimen.ghost) + 
                                   ((!ii)*abs(Blk_List.Block[i_blk].info.dimen.i)))*
                                   ((abs(jj)*Blk_List.Block[i_blk].info.dimen.ghost) + 
                                   ((!jj)*abs(Blk_List.Block[i_blk].info.dimen.j)))*
                                   ((abs(kk)*Blk_List.Block[i_blk].info.dimen.ghost) + 
                                   ((!kk)*abs(Blk_List.Block[i_blk].info.dimen.k)))*
                                   (Number_of_Solution_Variables);

                     if (!Blk_List.Block[neighbour_blk].used) return (2001);
                     neighbour_buffer_size = buffer_size;
                     if (buffer_size != neighbour_buffer_size) return (2002);
                  
                     neighbour_bound_elem = 
                       neighbour_info[i_bound_elem].blkorient.compute_message_tag(
                          neighbour_info[i_bound_elem].blkorient.direction_neighbour_to_me[0],
                          neighbour_info[i_bound_elem].blkorient.direction_neighbour_to_me[1], 
                          neighbour_info[i_bound_elem].blkorient.direction_neighbour_to_me[2]);
                  
                     // Perform direct message copies for the block.
                     for (int l = 0; l <= buffer_size-1; ++l) {
                        Blk_List.message_noreschange_recbuf[neighbour_blk][neighbour_bound_elem][l] =
                           Blk_List.message_noreschange_sendbuf[i_blk][i_bound_elem][l];
                     } /* endfor */
                  } /* endif */

               }/* end for k */
            }/* end for j */
         }/* end for i */

      } /* endif */
    }/* endfor */ 
        
    /* Message passing complete.  Return zero value. */
    
    return (0);
    
}

/**************************************************************
 * Routine: Exchange_Messages_NoResChange_Mesh_Geometry_Only  *
 *                                                            *
 * Sends and receives mesh and geometry information contained *
 * in the preloaded message passing buffers between           *
 * neighbouring adaptive blocks with no mesh resolution       *
 * changes.                                                   *
 *                                                            *
 **************************************************************/
int AdaptiveBlock3D_List::Exchange_Messages_NoResChange_Mesh_Geometry_Only(AdaptiveBlock3D_List &Blk_List) {

    int i_blk, neighbour_cpu, neighbour_blk, buffer_size, neighbour_buffer_size;

    int i_bound_elem, neighbour_bound_elem; // index for boundary element, face edge or vertex

    int number_neighbours[MAX_BOUNDARY_ELEMENTS_FOR_A_BLOCK];
    AdaptiveBlock3D_Info neighbour_info[MAX_BOUNDARY_ELEMENTS_FOR_A_BLOCK];
    
    /* Perform message passing for the block faces and corners of each 
       solution block having no mesh resolution change. */

    for (i_blk = 0 ; i_blk <= Blk_List.Nblk-1 ; ++i_blk) {
      if (Blk_List.Block[i_blk].used) {
         // Assign the boundary element information
         number_neighbours[BE::BSW] = Blk_List.Block[i_blk].nBSW;
         neighbour_info[BE::BSW] = Blk_List.Block[i_blk].infoBSW[0];
      
         number_neighbours[BE::SW] = Blk_List.Block[i_blk].nSW;
         neighbour_info[BE::SW] = Blk_List.Block[i_blk].infoSW[0];
      
         number_neighbours[BE::TSW] = Blk_List.Block[i_blk].nTSW;
         neighbour_info[BE::TSW] = Blk_List.Block[i_blk].infoTSW[0];
      
         number_neighbours[BE::BW] = Blk_List.Block[i_blk].nBW;
         neighbour_info[BE::BW] = Blk_List.Block[i_blk].infoBW[0];
      
         number_neighbours[BE::W] = Blk_List.Block[i_blk].nW;
         neighbour_info[BE::W] = Blk_List.Block[i_blk].infoW[0];
      
         number_neighbours[BE::TW] = Blk_List.Block[i_blk].nTW;
         neighbour_info[BE::TW] = Blk_List.Block[i_blk].infoTW[0];
      
         number_neighbours[BE::BNW] = Blk_List.Block[i_blk].nBNW;
         neighbour_info[BE::BNW] = Blk_List.Block[i_blk].infoBNW[0];
      
         number_neighbours[BE::NW] = Blk_List.Block[i_blk].nNW;
         neighbour_info[BE::NW] = Blk_List.Block[i_blk].infoNW[0];
      
         number_neighbours[BE::TNW] = Blk_List.Block[i_blk].nTNW;
         neighbour_info[BE::TNW] = Blk_List.Block[i_blk].infoTNW[0];
      
         number_neighbours[BE::BS] = Blk_List.Block[i_blk].nBS;
         neighbour_info[BE::BS] = Blk_List.Block[i_blk].infoBS[0];
      
         number_neighbours[BE::S] = Blk_List.Block[i_blk].nS;
         neighbour_info[BE::S] = Blk_List.Block[i_blk].infoS[0];
      
         number_neighbours[BE::TS] = Blk_List.Block[i_blk].nTS;
         neighbour_info[BE::TS] = Blk_List.Block[i_blk].infoTS[0];
      
         number_neighbours[BE::B] = Blk_List.Block[i_blk].nB;
         neighbour_info[BE::B] = Blk_List.Block[i_blk].infoB[0];
      
         number_neighbours[BE::T] = Blk_List.Block[i_blk].nT;
         neighbour_info[BE::T] = Blk_List.Block[i_blk].infoT[0];
      
         number_neighbours[BE::BN] = Blk_List.Block[i_blk].nBN;
         neighbour_info[BE::BN] = Blk_List.Block[i_blk].infoBN[0];
      
         number_neighbours[BE::N] = Blk_List.Block[i_blk].nN;
         neighbour_info[BE::N] = Blk_List.Block[i_blk].infoN[0];
      
         number_neighbours[BE::TN] = Blk_List.Block[i_blk].nTN;
         neighbour_info[BE::TN] = Blk_List.Block[i_blk].infoTN[0];
      
         number_neighbours[BE::BSE] = Blk_List.Block[i_blk].nBSE;
         neighbour_info[BE::BSE] = Blk_List.Block[i_blk].infoBSE[0];
      
         number_neighbours[BE::SE] = Blk_List.Block[i_blk].nSE;
         neighbour_info[BE::SE] = Blk_List.Block[i_blk].infoSE[0];
      
         number_neighbours[BE::TSE] = Blk_List.Block[i_blk].nTSE;
         neighbour_info[BE::TSE] = Blk_List.Block[i_blk].infoTSE[0];
      
         number_neighbours[BE::BE] = Blk_List.Block[i_blk].nBE;
         neighbour_info[BE::BE] = Blk_List.Block[i_blk].infoBE[0];
      
         number_neighbours[BE::E] = Blk_List.Block[i_blk].nE;
         neighbour_info[BE::E] = Blk_List.Block[i_blk].infoE[0];
      
         number_neighbours[BE::TE] = Blk_List.Block[i_blk].nTE;
         neighbour_info[BE::TE] = Blk_List.Block[i_blk].infoTE[0]; 
      
         number_neighbours[BE::BNE] = Blk_List.Block[i_blk].nBNE;
         neighbour_info[BE::BNE] = Blk_List.Block[i_blk].infoBNE[0]; 
       
         number_neighbours[BE::NE] = Blk_List.Block[i_blk].nNE;
         neighbour_info[BE::NE] = Blk_List.Block[i_blk].infoNE[0];  
      
         number_neighbours[BE::TNE] = Blk_List.Block[i_blk].nTNE;
         neighbour_info[BE::TNE] = Blk_List.Block[i_blk].infoTNE[0]; 

         // Perform direct copies of messages to appropriate message buffers.
         for (int ii = -1; ii <= 1; ii++) {
            for (int jj = -1; jj <= 1; jj++) {
               for (int kk = -1; kk <= 1; kk++) {
                  i_bound_elem = 9*(ii+1) + 3*(jj+1) + (kk+1);
               
                  if ((number_neighbours[i_bound_elem] == 1) && 
                      (i_bound_elem != BE::ME) &&
                      (Blk_List.Block[i_blk].info.level == neighbour_info[i_bound_elem].level)) {

                     neighbour_blk = neighbour_info[i_bound_elem].blknum;
                     
                     buffer_size = ((abs(ii)*Blk_List.Block[i_blk].info.dimen.ghost) + 
                                   ((!ii)*abs(Blk_List.Block[i_blk].info.dimen.i)+1))*
                                   ((abs(jj)*Blk_List.Block[i_blk].info.dimen.ghost) + 
                                   ((!jj)*abs(Blk_List.Block[i_blk].info.dimen.j)+1))*
                                   ((abs(kk)*Blk_List.Block[i_blk].info.dimen.ghost) + 
                                   ((!kk)*abs(Blk_List.Block[i_blk].info.dimen.k)+1))*
                                   (NUM_COMP_VECTOR3D) +
                                   (((!ii)*abs(Blk_List.Block[i_blk].info.dimen.i)))*
                                   (((!jj)*abs(Blk_List.Block[i_blk].info.dimen.j)))*
                                   (((!kk)*abs(Blk_List.Block[i_blk].info.dimen.k)));

                     if (!Blk_List.Block[neighbour_blk].used) return (2001);
                     neighbour_buffer_size = buffer_size;
                     if (buffer_size != neighbour_buffer_size) return (2002);
                  
                     neighbour_bound_elem = 
                       neighbour_info[i_bound_elem].blkorient.compute_message_tag(
                          neighbour_info[i_bound_elem].blkorient.direction_neighbour_to_me[0],
                          neighbour_info[i_bound_elem].blkorient.direction_neighbour_to_me[1], 
                          neighbour_info[i_bound_elem].blkorient.direction_neighbour_to_me[2]);
                  
                     // Perform direct message copies for the block.
                     for (int l = 0; l <= buffer_size-1; ++l) {
                        Blk_List.message_noreschange_recbuf[neighbour_blk][neighbour_bound_elem][l] =
                           Blk_List.message_noreschange_sendbuf[i_blk][i_bound_elem][l];
                     } /* endfor */
                  } /* endif */

               }/* end for k */
            }/* end for j */
         }/* end for i */

      } /* endif */
    }/* endfor */ 
        
    /* Message passing complete.  Return zero value. */
    
    return (0);
    
}

/**********************************************************
 * Routine: Exchange_Messages_ResChange_FineToCoarse      *
 *                                                        *
 * Sends and receives solution information contained in   *
 * the preloaded message passing buffers between          *
 * neighbouring adaptive blocks, sending solution         *
 * information from more refined (higher mesh resolution) *
 * adaptive blocks to more coarse adaptive solution       *
 * blocks with lower mesh resolution.                     *
 *                                                        *
 **********************************************************/
int AdaptiveBlock3D_List::Exchange_Messages_ResChange_FineToCoarse(AdaptiveBlock3D_List &Blk_List,
                                                                   const int Number_of_Solution_Variables) {

  cout<<"\n\nError Exchange_Messages_ResChange_FineToCoarse() does not exist for 3D\n";
  return (2);

}

/**********************************************************
 * Routine: Exchange_Messages_ResChange_CoarseToFine      *
 *                                                        *
 * Sends and receives solution information contained in   *
 * the preloaded message passing buffers between          *
 * neighbouring adaptive blocks, sending solution         *
 * information from more coarse (lower mesh resolution)   *
 * adaptive blocks to more refined adaptive solution      *
 * blocks with higher mesh resolution.                    *
 *                                                        *
 **********************************************************/
int AdaptiveBlock3D_List::Exchange_Messages_ResChange_CoarseToFine(AdaptiveBlock3D_List &Blk_List,
                                                                   const int Number_of_Solution_Variables) {

  cout<<"\n\nError Exchange_Messages_ResChange_CoarseToFine() does not exist for 3D\n";
  return (2);

}

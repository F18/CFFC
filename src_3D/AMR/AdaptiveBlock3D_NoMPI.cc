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
                                                        const int Number_of_Solution_Variables,
                                                        const int Send_Mesh_Geometry_Only) {

    int i_cpu, i_blk, neighbour_cpu, neighbour_blk, 
        buffer_size, buffer_size_neighbour, l;

    int tag_base_neigh;

    int i_bound_elem; // index for boundary element, face edge or vertex
    int num_bound_elem[MAX_BOUNDARY_ELEMENTS_FOR_A_BLOCK];
    AdaptiveBlock3D_Info info_bound_elem[MAX_BOUNDARY_ELEMENTS_FOR_A_BLOCK];
    
    /* Perform message passing for the block faces and corners of each 
       solution block having no mesh resolution change. */

    i_cpu = Blk_List.ThisCPU;

    for ( i_blk = 0 ; i_blk <= Blk_List.Nblk-1 ; ++i_blk ) {
      if (Blk_List.Block[i_blk].used) {
         // Assign the boundary element information
         num_bound_elem[BE::BSW] = Blk_List.Block[i_blk].nBSW;
         info_bound_elem[BE::BSW] =  Blk_List.Block[i_blk].infoBSW[0];
      
         num_bound_elem[BE::SW] = Blk_List.Block[i_blk].nSW;
         info_bound_elem[BE::SW] =  Blk_List.Block[i_blk].infoSW[0];
      
         num_bound_elem[BE::TSW] = Blk_List.Block[i_blk].nTSW;
         info_bound_elem[BE::TSW] =  Blk_List.Block[i_blk].infoTSW[0];
      
         num_bound_elem[BE::BW] = Blk_List.Block[i_blk].nBW;
         info_bound_elem[BE::BW] =  Blk_List.Block[i_blk].infoBW[0];
      
         num_bound_elem[BE::W] = Blk_List.Block[i_blk].nW;
         info_bound_elem[BE::W] =  Blk_List.Block[i_blk].infoW[0];
      
         num_bound_elem[BE::TW] = Blk_List.Block[i_blk].nTW;
         info_bound_elem[BE::TW] =  Blk_List.Block[i_blk].infoTW[0];
      
         num_bound_elem[BE::BNW] = Blk_List.Block[i_blk].nBNW;
         info_bound_elem[BE::BNW] =  Blk_List.Block[i_blk].infoBNW[0];
      
         num_bound_elem[BE::NW] = Blk_List.Block[i_blk].nNW;
         info_bound_elem[BE::NW] =  Blk_List.Block[i_blk].infoNW[0];
      
         num_bound_elem[BE::TNW] = Blk_List.Block[i_blk].nTNW;
         info_bound_elem[BE::TNW] =  Blk_List.Block[i_blk].infoTNW[0];
      
         num_bound_elem[BE::BS] = Blk_List.Block[i_blk].nBS;
         info_bound_elem[BE::BS] =  Blk_List.Block[i_blk].infoBS[0];
      
         num_bound_elem[BE::S] = Blk_List.Block[i_blk].nS;
         info_bound_elem[BE::S] =  Blk_List.Block[i_blk].infoS[0];
      
         num_bound_elem[BE::TS] = Blk_List.Block[i_blk].nTS;
         info_bound_elem[BE::TS] =  Blk_List.Block[i_blk].infoTS[0];
      
         num_bound_elem[BE::B] = Blk_List.Block[i_blk].nB;
         info_bound_elem[BE::B] =  Blk_List.Block[i_blk].infoB[0];
      
         num_bound_elem[BE::T] = Blk_List.Block[i_blk].nT;
         info_bound_elem[BE::T] =  Blk_List.Block[i_blk].infoT[0];
      
         num_bound_elem[BE::BN] = Blk_List.Block[i_blk].nBN;
         info_bound_elem[BE::BN] =  Blk_List.Block[i_blk].infoBN[0];
      
         num_bound_elem[BE::N] = Blk_List.Block[i_blk].nN;
         info_bound_elem[BE::N] =  Blk_List.Block[i_blk].infoN[0];
      
         num_bound_elem[BE::TN] = Blk_List.Block[i_blk].nTN;
         info_bound_elem[BE::TN] =  Blk_List.Block[i_blk].infoTN[0];
      
         num_bound_elem[BE::BSE] = Blk_List.Block[i_blk].nBSE;
         info_bound_elem[BE::BSE] =  Blk_List.Block[i_blk].infoBSE[0];
      
         num_bound_elem[BE::SE] = Blk_List.Block[i_blk].nSE;
         info_bound_elem[BE::SE] =  Blk_List.Block[i_blk].infoSE[0];
      
         num_bound_elem[BE::TSE] = Blk_List.Block[i_blk].nTSE;
         info_bound_elem[BE::TSE] =  Blk_List.Block[i_blk].infoTSE[0];
      
         num_bound_elem[BE::BE] = Blk_List.Block[i_blk].nBE;
         info_bound_elem[BE::BE] =  Blk_List.Block[i_blk].infoBE[0];
      
         num_bound_elem[BE::E] = Blk_List.Block[i_blk].nE;
         info_bound_elem[BE::E] =  Blk_List.Block[i_blk].infoE[0];
      
         num_bound_elem[BE::TE] = Blk_List.Block[i_blk].nTE;
         info_bound_elem[BE::TE] =  Blk_List.Block[i_blk].infoTE[0]; 
      
         num_bound_elem[BE::BNE] = Blk_List.Block[i_blk].nBNE;
         info_bound_elem[BE::BNE] =  Blk_List.Block[i_blk].infoBNE[0]; 
       
         num_bound_elem[BE::NE] = Blk_List.Block[i_blk].nNE;
         info_bound_elem[BE::NE] =  Blk_List.Block[i_blk].infoNE[0];  
      
         num_bound_elem[BE::TNE] = Blk_List.Block[i_blk].nTNE;
         info_bound_elem[BE::TNE] =  Blk_List.Block[i_blk].infoTNE[0]; 

         // Perform direct copies of messages to appropriate message buffers.
         for (int ii = -1; ii<2; ii++){
            for (int jj = -1; jj<2; jj++){
               for (int kk = -1; kk<2; kk++){
                  i_bound_elem = 9*(ii+1) + 3*(jj+1) + (kk+1);
               
                  if ((num_bound_elem[i_bound_elem] == 1) && 
                      (i_bound_elem != BE::ME) &&
                      (Blk_List.Block[i_blk].info.level == info_bound_elem[i_bound_elem].level)) {
                     neighbour_blk = info_bound_elem[i_bound_elem].blknum;
                     if (!Blk_List.Block[neighbour_blk].used) return(2401);
                     buffer_size = ((abs(ii)*Blk_List.Block[i_blk].info.dimen.ghost) + 
                                   ((!ii)*abs(Blk_List.Block[i_blk].info.dimen.i)))*
                                   ((abs(jj)*Blk_List.Block[i_blk].info.dimen.ghost) + 
                                   ((!jj)*abs(Blk_List.Block[i_blk].info.dimen.j)))*
                                   ((abs(kk)*Blk_List.Block[i_blk].info.dimen.ghost) + 
                                   ((!kk)*abs(Blk_List.Block[i_blk].info.dimen.k)))*
                                   (Number_of_Solution_Variables);
                     buffer_size_neighbour = buffer_size;
                  
                     tag_base_neigh = 
                       info_bound_elem[i_bound_elem].blkorient.compute_message_tag(
                          info_bound_elem[i_bound_elem].blkorient.direction_neighbour_to_me[0],
                          info_bound_elem[i_bound_elem].blkorient.direction_neighbour_to_me[1], 
                          info_bound_elem[i_bound_elem].blkorient.direction_neighbour_to_me[2]);
                  
                     for (l = 0; l <= buffer_size-1; ++l) {
                        Blk_List.message_noreschange_recbuf[neighbour_blk][tag_base_neigh][l] =
                           Blk_List.message_noreschange_sendbuf[i_blk][i_bound_elem][l];
                     } /* endfor */

                  } /* endif */

               }/* end for k */
            }/* end for j */
         }/* end for i */

      } /* endif */
    }/* endfor */ 
        
    /* Message passing complete.  Return zero value. */
    
    return(0);
    
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

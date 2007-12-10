/* Octree.cc:  Subroutines for octree adaptive blocks 
               hierarchical data structure. */

/* Include octree header file. */

#ifndef _OCTREE_INCLUDED
#include "Octree.h"
#endif // _OCTREE_INCLUDED

/**********************************************************
 * Routine: Create_Octree_Data_Structure                  *
 *                                                        *
 * Create (initialize) a octree adaptive block            *
 * hierarchical data structure with root dimensions       *
 *                                                        *
 * Number_of_Roots_Idir X Number_of_Roots_Jdir            *
 *                                                        *
 * and a global block list of dimension                   *
 *                                                        *
 * Number_of_Processors X Number_of_Blocks_per_Processor. *
 *                                                        *
 **********************************************************/
void Octree_DataStructure::Create_Octree_Data_Structure(Octree_DataStructure &Octree,
                                                        const int Number_of_Roots_Idir,
  	                                                const int Number_of_Roots_Jdir,
  	                                                const int Number_of_Roots_Kdir,
                                                        const int Number_of_Processors,
                                                        const int Number_of_Blocks_per_Processor) {

    /* Allocate (re-allocate) memory for the octree
       hierarchical data structure. */

    if (Octree.Roots != NULL && Octree.Blocks != NULL) {
       Octree.deallocate();
    } else if (Octree.Roots != NULL) {
       Octree.deallocateRoots();
    } else if (Octree.Blocks != NULL) {
       Octree.deallocateBlocks();
    } /* endif */

    Octree.allocate(Number_of_Roots_Idir, 
                    Number_of_Roots_Jdir,
                    Number_of_Roots_Kdir,
                    Number_of_Processors,
                    Number_of_Blocks_per_Processor);

}

/********************************************************
 * Routine: Broadcast_Octree_Data_Structure             *
 *                                                      *
 * Broadcast octree data structure to all processors    *
 * involved in the calculation from the primary         *
 * processor using the MPI broadcast routine.           *
 *                                                      *
 ********************************************************/
void Octree_DataStructure::Broadcast_Octree_Data_Structure(Octree_DataStructure &Octree,
                                                           AdaptiveBlock3D_ResourceList &List_of_Available_Blocks) {

#ifdef _MPI_VERSION
    int nr, nri, nrj, nrk, ncpu, nblk, iBLK, jBLK, kBLK;
    
    /* Broadcast the number of roots, the number of CPUs, and
       the number of local solution blocks. */

    if (CFFC_Primary_MPI_Processor()) {
       nr = Octree.NR;
       nri = Octree.NRi;
       nrj = Octree.NRj;
       nrk = Octree.NRk;
       ncpu = Octree.Ncpu;
       nblk = Octree.Nblk;
    } /* endif */

    MPI::COMM_WORLD.Bcast(&nr, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&nri, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&nrj, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&nrk, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&ncpu, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&nblk, 1, MPI::INT, 0);

    /* On non-primary MPI processors, allocate (re-allocate) 
       memory for the octree data structure as necessary. */

    if (!CFFC_Primary_MPI_Processor()) {
       if (Octree.NR != nr ||
           Octree.NRi != nri ||
           Octree.NRj != nrj ||
           Octree.NRk != nrk ||
           Octree.Ncpu != ncpu ||
           Octree.Nblk != nblk) {
          Create_Octree_Data_Structure(Octree,
                                       nri,
  	                               nrj,
  	                               nrk,
                                       ncpu,
                                       nblk);
       } /* endif */
    } /* endif */

    /* Broadcast the octree data structure, descending each of the
       roots in a recursive fashion. */
    
    for (int n = 0 ; n <= Octree.NR-1 ; ++n) {
       Octree.Roots[n].broadcast(List_of_Available_Blocks);
    } /* endfor */
#endif

}

/**********************************************************
 * Routine: Renumber_Solution_Blocks                      *
 *                                                        *
 * Assigns a global block number to all of the solution   *
 * blocks in the octree adaptive block hierarchical       *
 * data structure.                                        *
 *                                                        *
 **********************************************************/
void Octree_DataStructure::Renumber_Solution_Blocks(Octree_DataStructure &Octree) {

   Octree.renumber();

}

/**********************************************************
 * Routine: Renumber_Solution_Blocks                      *
 *                                                        *
 * Assigns a global block number to all of the solution   *
 * blocks in the octree adaptive block hierarchical       *
 * data structure.                                        *
  *                                                       *
 **********************************************************/
void Octree_DataStructure::Renumber_Solution_Blocks(Octree_DataStructure &Octree,
                                                    AdaptiveBlock3D_List &LocalSolnBlockList) {

    int iCPU;

    Octree.renumber();
    iCPU = LocalSolnBlockList.ThisCPU;

    for (int iBLK = 0 ; iBLK <= LocalSolnBlockList.Nblk-1 ; ++iBLK ) {
       if (Octree.Blocks[iCPU][iBLK] != NULL) {
          if (Octree.Blocks[iCPU][iBLK]->block.used) {
             LocalSolnBlockList.Block[iBLK].gblknum = 
                Octree.Blocks[iCPU][iBLK]->block.gblknum;
          } /* endif */
       } /* endif */
    } /* endfor */    

}


/**********************************************************
 * Routine: Modify_Neighbours_of_Root_Solution_Blocks     *
 *                                                        *
 * Modifies the neigbours of all root solution blocks in  *
 * the octree adaptive block data structure, depending    *
 * on the multi-block grid type.                          *
 *                                                        *
 **********************************************************/
void Octree_DataStructure::Modify_Neighbours_of_Root_Solution_Blocks(Octree_DataStructure &Octree,
                                                                     const int Grid_Type) {

   switch(Grid_Type) {
     case GRID_CUBE :
       break;
     case GRID_CHANNEL :
       break;
     case GRID_COUETTE :
       break;
     default:
       break;
   }  /* endswitch */

}

/**********************************************************
 * Routine: Modify_Neighbours_of_Root_Solution_Blocks     *
 *                                                        *
 * Modifies the neigbours of all root solution blocks in  *
 * the octree adaptive block data structure, depending  *
 * on the multi-block grid type.                          *
 *                                                        *
 **********************************************************/
void Octree_DataStructure::Modify_Neighbours_of_Root_Solution_Blocks(Octree_DataStructure &Octree,
                                                                     AdaptiveBlock3D_List &LocalSolnBlockList,
                                                                     const int Grid_Type) {

   Modify_Neighbours_of_Root_Solution_Blocks(Octree, Grid_Type);

   switch(Grid_Type) {
     case GRID_CUBE :
       break;
     case GRID_CHANNEL :
       break;
     case GRID_COUETTE :
       break;
     default:
       break;
   }  /* endswitch */

}

/**********************************************************
 * Routine: Find_Neighbours                               *
 *                                                        *
 * Determines and stores the neigbours of all solution    *
 * blocks in the octree adaptive block data structure.  *
 *                                                        *
 **********************************************************/
void Octree_DataStructure::Find_Neighbours(Octree_DataStructure &Octree) {

   Octree.findNeighbours();

}

/**********************************************************
 * Routine: Find_Neighbours                               *
 *                                                        *
 * Determines and stores the neigbours of all solution    *
 * blocks in the octree adaptive block data structure.  *
 *                                                        *
 **********************************************************/
void Octree_DataStructure::Find_Neighbours(Octree_DataStructure &Octree,
                                           AdaptiveBlock3D_List &LocalSolnBlockList) {

    int iCPU;

    Octree.findNeighbours();
    iCPU = LocalSolnBlockList.ThisCPU;

    for (int iBLK = 0 ; iBLK <= LocalSolnBlockList.Nblk-1 ; ++iBLK ) {
       if (Octree.Blocks[iCPU][iBLK] != NULL) {
          if (Octree.Blocks[iCPU][iBLK]->block.used) {
             LocalSolnBlockList.Block[iBLK] = 
                Octree.Blocks[iCPU][iBLK]->block;
          } /* endif */
       } /* endif */
    } /* endfor */

}

/**********************************************************
 * Routine: Get_Refinement_List                           *
 *                                                        *
 * Obtains the a global list of all solution blocks in    *
 * the octree adaptive block hierarchical data          *
 * structure scheduled for mesh refinement or coarsening. *
 *                                                        *
 **********************************************************/
void Octree_DataStructure::Get_Refinement_List(Octree_DataStructure &Octree,
                                               AdaptiveBlock3D_List &LocalSolnBlockList) {

  cout<<"\nError Get_Refinement_List() has not been implemented for Octree\n";

}

/* Octree.cc:  Subroutines for octree adaptive blocks 
               hierarchical data structure. */

/* Include octree header file. */

#ifndef _OCTREE_INCLUDED
#include "Octree.h"
#endif // _OCTREE_INCLUDED

/**********************************************************
 * Routine: Create_OcTree_Data_Structure                *
 *                                                        *
 * Create (initialize) a octree adaptive block          *
 * hierarchical data structure with root dimensions       *
 *                                                        *
 * Number_of_Roots_Idir X Number_of_Roots_Jdir            *
 *                                                        *
 * and a global block list of dimension                   *
 *                                                        *
 * Number_of_Processors X Number_of_Blocks_per_Processor. *
 *                                                        *
 **********************************************************/
void  OcTreeBlock_DataStructure::Create_OcTree_Data_Structure(OcTreeBlock_DataStructure &OcTree,
                                    const int Number_of_Roots_Idir,
  	                            const int Number_of_Roots_Jdir,
  	                            const int Number_of_Roots_Kdir,
                                    const int Number_of_Processors,
                                    const int Number_of_Blocks_per_Processor) {

    /* Allocate (re-allocate) memory for the octree
       hierarchical data structure. */

    if (OcTree.Roots != NULL && OcTree.Blocks != NULL) {
       OcTree.deallocate();
    } else if (OcTree.Roots != NULL) {
       OcTree.deallocateRoots();
    } else if (OcTree.Blocks != NULL) {
       OcTree.deallocateBlocks();
    } /* endif */
    OcTree.allocate(Number_of_Roots_Idir, 
                      Number_of_Roots_Jdir,
                      Number_of_Roots_Kdir,
                      Number_of_Processors,
                      Number_of_Blocks_per_Processor);

}

/********************************************************
 * Routine: Broadcast_OcTree_Data_Structure           *
 *                                                      *
 * Broadcast octree data structure to all processors  *
 * involved in the calculation from the primary         *
 * processor using the MPI broadcast routine.           *
 *                                                      *
 ********************************************************/
void  OcTreeBlock_DataStructure::Broadcast_OcTree_Data_Structure(OcTreeBlock_DataStructure &OcTree,
                                       AdaptiveBlock3DResourceList   &List_of_Available_Blocks) {

#ifdef _MPI_VERSION
    int nri, nrj, nrk, ncpu, nblk, iBLK, jBLK, kBLK;

    /* Broadcast the number of roots, the number of CPUs, and
       the number of local solution blocks. */

    if (CFFC_Primary_MPI_Processor()) {
       nri = OcTree.NRi;
       nrj = OcTree.NRj;
       nrk = OcTree.NRk;
       ncpu = OcTree.Ncpu;
       nblk = OcTree.Nblk;
    } /* endif */

    MPI::COMM_WORLD.Bcast(&nri, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&nrj, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&nrk, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&ncpu, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&nblk, 1, MPI::INT, 0);

    /* On non-primary MPI processors, allocate (re-allocate) 
       memory for the octree data structure as necessary. */

    if (!CFFC_Primary_MPI_Processor()) {
       if (OcTree.NRi != nri ||
           OcTree.NRj != nrj ||
           OcTree.NRk != nrk ||
           OcTree.Ncpu != ncpu ||
           OcTree.Nblk != nblk) {
          Create_OcTree_Data_Structure(OcTree,
                                         nri,
  	                                 nrj,
  	                                 nrk,
                                         ncpu,
                                         nblk);
       } /* endif */
    } /* endif */

    /* Broadcast the octree data structure, descending each of the
       roots in a recursive fashion. */

    for ( kBLK = 0 ; kBLK <= OcTree.NRk-1 ; ++kBLK ) 
      for ( jBLK = 0 ; jBLK <= OcTree.NRj-1 ; ++jBLK ) 
        for ( iBLK = 0 ; iBLK <= OcTree.NRi-1 ; ++iBLK ) {
	   OcTree.Roots[iBLK][jBLK][kBLK].broadcast(List_of_Available_Blocks);
        } /* endfor */

#endif

}

/**********************************************************
 * Routine: Renumber_Solution_Blocks                      *
 *                                                        *
 * Assigns a global block number to all of the solution   *
 * blocks in the octree adaptive block hierarchical     *
 * data structure.                                        *
 *                                                        *
 **********************************************************/
void  OcTreeBlock_DataStructure::Renumber_Solution_Blocks(OcTreeBlock_DataStructure &OcTree) {

    OcTree.renumber();

}

/**********************************************************
 * Routine: Renumber_Solution_Blocks                      *
 *                                                        *
 * Assigns a global block number to all of the solution   *
 * blocks in the octree adaptive block hierarchical     *
 * data structure.                                        *
  *                                                        *
 **********************************************************/
void  OcTreeBlock_DataStructure::Renumber_Solution_Blocks(OcTreeBlock_DataStructure &OcTree,
                              AdaptiveBlock3D_List &LocalSolnBlockList) {

    int iCPU, iBLK, global_block_number;

    OcTree.renumber();
   
    iCPU = LocalSolnBlockList.ThisCPU;
    for ( iBLK = 0 ; iBLK <= LocalSolnBlockList.Nblk-1 ; ++iBLK ) {
       if (OcTree.Blocks[iCPU][iBLK] != NULL) {
          if (OcTree.Blocks[iCPU][iBLK]->block.used) {
             LocalSolnBlockList.Block[iBLK].gblknum = 
                OcTree.Blocks[iCPU][iBLK]->block.gblknum;
          } /* endif */
       } /* endif */
    } /* endfor */    

}

/**********************************************************
 * Routine: Find_Neighbours_of_Root_Solution_Blocks       *
 *                                                        *
 * Determines and stores the neighbours of all root       *
 * solution blocks in the octree adaptive block data    *
 * structure.  It is assumed that all root solution       *
 * blocks are initially at the same level of refinement.  *
 *                                                        *
 **********************************************************/
void  OcTreeBlock_DataStructure::Find_Neighbours_of_Root_Solution_Blocks(
   OcTreeBlock_DataStructure 
   &OcTree) {

    OcTree.findRootNeighbours();

}

/**********************************************************
 * Routine: Find_Neighbours_of_Root_Solution_Blocks       *
 *                                                        *
 * Determines and stores the neighbours of all root       *
 * solution blocks in the octree adaptive block data    *
 * structure.  It is assumed that all root solution       *
 * blocks are initially at the same level of refinement.  *
 *                                                        *
 **********************************************************/
void  OcTreeBlock_DataStructure::Find_Neighbours_of_Root_Solution_Blocks(
   OcTreeBlock_DataStructure &OcTree,
   AdaptiveBlock3D_List &LocalSolnBlockList) {

    int iBLK, jBLK, kBLK;

    OcTree.findRootNeighbours();

    for ( kBLK = 0 ; kBLK <= OcTree.NRk-1 ; ++kBLK ) 
      for ( jBLK = 0 ; jBLK <= OcTree.NRj-1 ; ++jBLK ) 
        for ( iBLK = 0 ; iBLK <= OcTree.NRi-1 ; ++iBLK ) {

           if (OcTree.Roots[iBLK][jBLK][kBLK].block.used) {
	      // Assign neighbour information to local block list.
	      if (OcTree.Roots[iBLK][jBLK][kBLK].block.nW == 1) {
                 if (LocalSolnBlockList.ThisCPU == 
                     OcTree.Roots[iBLK][jBLK][kBLK].block.info.cpu) {
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].nW = 1;
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].infoW[0] =
                       OcTree.Roots[iBLK-1][jBLK][kBLK].block.info;
                 } /* endif */
              } /* endif */

	      if (OcTree.Roots[iBLK][jBLK][kBLK].block.nE == 1) {
                 if (LocalSolnBlockList.ThisCPU == 
                     OcTree.Roots[iBLK][jBLK][kBLK].block.info.cpu) {
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].nE = 1;
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].infoE[0] =
                       OcTree.Roots[iBLK+1][jBLK][kBLK].block.info;
                 } /* endif */
              } /* endif */
	   
	      if (OcTree.Roots[iBLK][jBLK][kBLK].block.nS == 1) {
                 if (LocalSolnBlockList.ThisCPU == 
                     OcTree.Roots[iBLK][jBLK][kBLK].block.info.cpu) {
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].nS = 1;
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].infoS[0] =
                       OcTree.Roots[iBLK][jBLK-1][kBLK].block.info;
                 } /* endif */
              } /* endif */
	   
	      if (OcTree.Roots[iBLK][jBLK][kBLK].block.nN == 1) {
                 if (LocalSolnBlockList.ThisCPU == 
                     OcTree.Roots[iBLK][jBLK][kBLK].block.info.cpu) {
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].nN = 1;
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].infoN[0] =
                       OcTree.Roots[iBLK][jBLK+1][kBLK].block.info;
                 } /* endif */
              } /* endif */
	   
	      if (OcTree.Roots[iBLK][jBLK][kBLK].block.nT == 1) {
                 if (LocalSolnBlockList.ThisCPU == 
                     OcTree.Roots[iBLK][jBLK][kBLK].block.info.cpu) {
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].nT = 1;
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].infoT[0] =
                       OcTree.Roots[iBLK][jBLK][kBLK+1].block.info;
                 } /* endif */
              } /* endif */

	      if (OcTree.Roots[iBLK][jBLK][kBLK].block.nB == 1) {
                 if (LocalSolnBlockList.ThisCPU == 
                     OcTree.Roots[iBLK][jBLK][kBLK].block.info.cpu) {
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].nB = 1;
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].infoB[0] =
                       OcTree.Roots[iBLK][jBLK][kBLK-1].block.info;
                 } /* endif */
              } /* endif */

	      if (OcTree.Roots[iBLK][jBLK][kBLK].block.nNW == 1) {
                 if (LocalSolnBlockList.ThisCPU == 
                     OcTree.Roots[iBLK][jBLK][kBLK].block.info.cpu) {
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].nNW = 1;
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].infoNW[0] =
                       OcTree.Roots[iBLK-1][jBLK+1][kBLK].block.info;
                 } /* endif */
              } /* endif */
	   
	      if (OcTree.Roots[iBLK][jBLK][kBLK].block.nSW == 1) {
                 if (LocalSolnBlockList.ThisCPU == 
                     OcTree.Roots[iBLK][jBLK][kBLK].block.info.cpu) {
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].nSW = 1;
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].infoSW[0] =
                       OcTree.Roots[iBLK-1][jBLK-1][kBLK].block.info;
                 } /* endif */
              } /* endif */
	   
	      if (OcTree.Roots[iBLK][jBLK][kBLK].block.nNE == 1) {
                 if (LocalSolnBlockList.ThisCPU == 
                     OcTree.Roots[iBLK][jBLK][kBLK].block.info.cpu) {
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].nNE = 1;
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].infoNE[0] =
                       OcTree.Roots[iBLK+1][jBLK+1][kBLK].block.info;
                 } /* endif */
              } /* endif */
	   
	      if (OcTree.Roots[iBLK][jBLK][kBLK].block.nSE == 1) {
                 if (LocalSolnBlockList.ThisCPU == 
                     OcTree.Roots[iBLK][jBLK][kBLK].block.info.cpu) {
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].nSE = 1;
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].infoSE[0] =
                       OcTree.Roots[iBLK+1][jBLK-1][kBLK].block.info;
                 } /* endif */
              } /* endif */

	      if (OcTree.Roots[iBLK][jBLK][kBLK].block.nTN == 1) {
                 if (LocalSolnBlockList.ThisCPU == 
                     OcTree.Roots[iBLK][jBLK][kBLK].block.info.cpu) {
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].nTN = 1;
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].infoTN[0] =
                       OcTree.Roots[iBLK][jBLK+1][kBLK+1].block.info;

                    //  if(iBLK==3 && jBLK==0 && kBLK==0)cout<<"\n .infoTN.dimen.i for (3,0,0) and blk num"<<OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum<<" is "<<LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].infoTN[0].dimen.i;

                 } /* endif */
              } /* endif */
	   
	      if (OcTree.Roots[iBLK][jBLK][kBLK].block.nTS == 1) {
                 if (LocalSolnBlockList.ThisCPU == 
                     OcTree.Roots[iBLK][jBLK][kBLK].block.info.cpu) {
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].nTS = 1;
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].infoTS[0] =
                       OcTree.Roots[iBLK][jBLK-1][kBLK+1].block.info;
                 } /* endif */
              } /* endif */
	   
	      if (OcTree.Roots[iBLK][jBLK][kBLK].block.nTE == 1) {
                 if (LocalSolnBlockList.ThisCPU == 
                     OcTree.Roots[iBLK][jBLK][kBLK].block.info.cpu) {
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].nTE = 1;
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].infoTE[0] =
                       OcTree.Roots[iBLK+1][jBLK][kBLK+1].block.info;
                 } /* endif */
              } /* endif */
	   
	      if (OcTree.Roots[iBLK][jBLK][kBLK].block.nTW == 1) {
                 if (LocalSolnBlockList.ThisCPU == 
                     OcTree.Roots[iBLK][jBLK][kBLK].block.info.cpu) {
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].nTW = 1;
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].infoTW[0] =
                       OcTree.Roots[iBLK-1][jBLK][kBLK+1].block.info;
                 } /* endif */
              } /* endif */


	      if (OcTree.Roots[iBLK][jBLK][kBLK].block.nTNW == 1) {
                 if (LocalSolnBlockList.ThisCPU == 
                     OcTree.Roots[iBLK][jBLK][kBLK].block.info.cpu) {
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].nTNW = 1;
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].infoTNW[0] =
                       OcTree.Roots[iBLK-1][jBLK+1][kBLK+1].block.info;
                 } /* endif */
              } /* endif */
	   
	      if (OcTree.Roots[iBLK][jBLK][kBLK].block.nTSW == 1) {
                 if (LocalSolnBlockList.ThisCPU == 
                     OcTree.Roots[iBLK][jBLK][kBLK].block.info.cpu) {
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].nTSW = 1;
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].infoTSW[0] =
                       OcTree.Roots[iBLK-1][jBLK-1][kBLK+1].block.info;
                 } /* endif */
              } /* endif */
	   
	      if (OcTree.Roots[iBLK][jBLK][kBLK].block.nTNE == 1) {
                 if (LocalSolnBlockList.ThisCPU == 
                     OcTree.Roots[iBLK][jBLK][kBLK].block.info.cpu) {
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].nTNE = 1;
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].infoTNE[0] =
                       OcTree.Roots[iBLK+1][jBLK+1][kBLK+1].block.info;
                 } /* endif */
              } /* endif */
	   
	      if (OcTree.Roots[iBLK][jBLK][kBLK].block.nTSE == 1) {
                 if (LocalSolnBlockList.ThisCPU == 
                     OcTree.Roots[iBLK][jBLK][kBLK].block.info.cpu) {
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].nTSE = 1;
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].infoTSE[0] =
                       OcTree.Roots[iBLK+1][jBLK-1][kBLK+1].block.info;
                 } /* endif */
              } /* endif */

	      if (OcTree.Roots[iBLK][jBLK][kBLK].block.nBN == 1) {
                 if (LocalSolnBlockList.ThisCPU == 
                     OcTree.Roots[iBLK][jBLK][kBLK].block.info.cpu) {
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].nBN = 1;
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].infoBN[0] =
                       OcTree.Roots[iBLK][jBLK+1][kBLK-1].block.info;
                 } /* endif */
              } /* endif */
	   
	      if (OcTree.Roots[iBLK][jBLK][kBLK].block.nBS == 1) {
                 if (LocalSolnBlockList.ThisCPU == 
                     OcTree.Roots[iBLK][jBLK][kBLK].block.info.cpu) {
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].nBS = 1;
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].infoBS[0] =
                       OcTree.Roots[iBLK][jBLK-1][kBLK-1].block.info;
                 } /* endif */
              } /* endif */
	   
	      if (OcTree.Roots[iBLK][jBLK][kBLK].block.nBE == 1) {
                 if (LocalSolnBlockList.ThisCPU == 
                     OcTree.Roots[iBLK][jBLK][kBLK].block.info.cpu) {
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].nBE = 1;
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].infoBE[0] =
                       OcTree.Roots[iBLK+1][jBLK][kBLK-1].block.info;
                 } /* endif */
              } /* endif */
	   
	      if (OcTree.Roots[iBLK][jBLK][kBLK].block.nBW == 1) {
                 if (LocalSolnBlockList.ThisCPU == 
                     OcTree.Roots[iBLK][jBLK][kBLK].block.info.cpu) {
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].nBW = 1;
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].infoBW[0] =
                       OcTree.Roots[iBLK-1][jBLK][kBLK-1].block.info;
                 } /* endif */
              } /* endif */


	      if (OcTree.Roots[iBLK][jBLK][kBLK].block.nBNW == 1) {
                 if (LocalSolnBlockList.ThisCPU == 
                     OcTree.Roots[iBLK][jBLK][kBLK].block.info.cpu) {
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].nBNW = 1;
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].infoBNW[0] =
                       OcTree.Roots[iBLK-1][jBLK+1][kBLK-1].block.info;
                 } /* endif */
              } /* endif */
	   
	      if (OcTree.Roots[iBLK][jBLK][kBLK].block.nBSW == 1) {
                 if (LocalSolnBlockList.ThisCPU == 
                     OcTree.Roots[iBLK][jBLK][kBLK].block.info.cpu) {
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].nBSW = 1;
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].infoBSW[0] =
                       OcTree.Roots[iBLK-1][jBLK-1][kBLK-1].block.info;
                 } /* endif */
              } /* endif */
	   
	      if (OcTree.Roots[iBLK][jBLK][kBLK].block.nBNE == 1) {
                 if (LocalSolnBlockList.ThisCPU == 
                     OcTree.Roots[iBLK][jBLK][kBLK].block.info.cpu) {
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].nBNE = 1;
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].infoBNE[0] =
                       OcTree.Roots[iBLK+1][jBLK+1][kBLK-1].block.info;
                 } /* endif */
              } /* endif */
	   
	      if (OcTree.Roots[iBLK][jBLK][kBLK].block.nBSE == 1) {
                 if (LocalSolnBlockList.ThisCPU == 
                     OcTree.Roots[iBLK][jBLK][kBLK].block.info.cpu) {
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].nBSE = 1;
                    LocalSolnBlockList.Block[OcTree.Roots[iBLK][jBLK][kBLK].block.info.blknum].infoBSE[0] =
                       OcTree.Roots[iBLK+1][jBLK-1][kBLK-1].block.info;
                 } /* endif */
              } /* endif */


           } /* endif */
        } /* endfor */

}

/**********************************************************
 * Routine: Modify_Neighbours_of_Root_Solution_Blocks     *
 *                                                        *
 * Modifies the neigbours of all root solution blocks in  *
 * the octree adaptive block data structure, depending  *
 * on the multi-block grid type.                          *
 *                                                        *
 **********************************************************/
void  OcTreeBlock_DataStructure::Modify_Neighbours_of_Root_Solution_Blocks(OcTreeBlock_DataStructure 
                                               &OcTree,
                                               const int Grid_Type) {
     switch(Grid_Type) {
      case GRID_CUBE :
//	cout<<"\nModifications have not been implemented for GRID_CUBE";
        break;
     case GRID_CHANNEL :
//	cout<<"\nModifications have not been implemented for GRID_CUBE";
        break;
     case GRID_COUETTE :
//	cout<<"\nModifications have not been implemented for GRID_CUBE";
     break;
     default:
       cout<<"\nModify_Neighbours_of_Root_Solution_Blocks() is not defined for OcTree";
     }
}

/**********************************************************
 * Routine: Modify_Neighbours_of_Root_Solution_Blocks     *
 *                                                        *
 * Modifies the neigbours of all root solution blocks in  *
 * the octree adaptive block data structure, depending  *
 * on the multi-block grid type.                          *
 *                                                        *
 **********************************************************/
void  OcTreeBlock_DataStructure::Modify_Neighbours_of_Root_Solution_Blocks(OcTreeBlock_DataStructure 
                                               &OcTree,
                                               AdaptiveBlock3D_List 
                                               &LocalSolnBlockList,
                                               const int Grid_Type) {
    Modify_Neighbours_of_Root_Solution_Blocks(OcTree, Grid_Type);

    switch(Grid_Type) {
      case GRID_SQUARE :
	//Not implemented yet
        break;
    }


}

/**********************************************************
 * Routine: Find_Neighbours                               *
 *                                                        *
 * Determines and stores the neigbours of all solution    *
 * blocks in the octree adaptive block data structure.  *
 *                                                        *
 **********************************************************/
void  OcTreeBlock_DataStructure::Find_Neighbours(OcTreeBlock_DataStructure &OcTree) {

    OcTree.findNeighbours();

}

/**********************************************************
 * Routine: Find_Neighbours                               *
 *                                                        *
 * Determines and stores the neigbours of all solution    *
 * blocks in the octree adaptive block data structure.  *
 *                                                        *
 **********************************************************/
void  OcTreeBlock_DataStructure::Find_Neighbours(OcTreeBlock_DataStructure &OcTree,
                     AdaptiveBlock3D_List &LocalSolnBlockList) {

    int iCPU, iBLK, global_block_number;

    OcTree.findNeighbours();
   
    iCPU = LocalSolnBlockList.ThisCPU;
    for ( iBLK = 0 ; iBLK <= LocalSolnBlockList.Nblk-1 ; ++iBLK ) {
       if (OcTree.Blocks[iCPU][iBLK] != NULL) {
          if (OcTree.Blocks[iCPU][iBLK]->block.used) {
             LocalSolnBlockList.Block[iBLK] = 
                OcTree.Blocks[iCPU][iBLK]->block;
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
void  OcTreeBlock_DataStructure::Get_Refinement_List(OcTreeBlock_DataStructure &OcTree,
                         AdaptiveBlock3D_List &LocalSolnBlockList) {

  cout<<"\nError Get_Refinement_List() has not been implemented for OcTree\n";
 assert(1==2);
}

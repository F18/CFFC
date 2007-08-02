/* QuadTree.cc:  Subroutines for quadtree adaptive blocks 
                 hierarchical data structure. */

/* Include quadtree header file. */

#ifndef _QUADTREE_INCLUDED
#include "QuadTree.h"
#endif // _QUADTREE_INCLUDED

/**********************************************************
 * Routine: Create_QuadTree_Data_Structure                *
 *                                                        *
 * Create (initialize) a quadtree adaptive block          *
 * hierarchical data structure with root dimensions       *
 *                                                        *
 * Number_of_Roots_Idir X Number_of_Roots_Jdir            *
 *                                                        *
 * and a global block list of dimension                   *
 *                                                        *
 * Number_of_Processors X Number_of_Blocks_per_Processor. *
 *                                                        *
 **********************************************************/
void Create_QuadTree_Data_Structure(QuadTreeBlock_DataStructure &QuadTree,
                                    const int Number_of_Roots_Idir,
  	                            const int Number_of_Roots_Jdir,
                                    const int Number_of_Processors,
                                    const int Number_of_Blocks_per_Processor) {

    /* Allocate (re-allocate) memory for the quadtree
       hierarchical data structure. */

    if (QuadTree.Roots != NULL && QuadTree.Blocks != NULL) {
       QuadTree.deallocate();
    } else if (QuadTree.Roots != NULL) {
       QuadTree.deallocateRoots();
    } else if (QuadTree.Blocks != NULL) {
       QuadTree.deallocateBlocks();
    } /* endif */
    QuadTree.allocate(Number_of_Roots_Idir, 
                      Number_of_Roots_Jdir,
                      Number_of_Processors,
                      Number_of_Blocks_per_Processor);

}

/********************************************************
 * Routine: Broadcast_QuadTree_Data_Structure           *
 *                                                      *
 * Broadcast quadtree data structure to all processors  *
 * involved in the calculation from the primary         *
 * processor using the MPI broadcast routine.           *
 *                                                      *
 ********************************************************/
void Broadcast_QuadTree_Data_Structure(QuadTreeBlock_DataStructure &QuadTree,
                                       AdaptiveBlockResourceList   &List_of_Available_Blocks) {

#ifdef _MPI_VERSION
    int nri, nrj, ncpu, nblk, iBLK, jBLK;

    /* Broadcast the number of roots, the number of CPUs, and
       the number of local solution blocks. */

    if (CFDkit_Primary_MPI_Processor()) {
       nri = QuadTree.NRi;
       nrj = QuadTree.NRj;
       ncpu = QuadTree.Ncpu;
       nblk = QuadTree.Nblk;
    } /* endif */

    MPI::COMM_WORLD.Bcast(&nri, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&nrj, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&ncpu, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&nblk, 1, MPI::INT, 0);

    /* On non-primary MPI processors, allocate (re-allocate) 
       memory for the quadtree data structure as necessary. */

    if (!CFDkit_Primary_MPI_Processor()) {
       if (QuadTree.NRi != nri ||
           QuadTree.NRj != nrj ||
           QuadTree.Ncpu != ncpu ||
           QuadTree.Nblk != nblk) {
          Create_QuadTree_Data_Structure(QuadTree,
                                         nri,
  	                                 nrj,
                                         ncpu,
                                         nblk);
       } /* endif */
    } /* endif */

    /* Broadcast the quadtree data structure, descending each of the
       roots in a recursive fashion. */

    for ( jBLK = 0 ; jBLK <= QuadTree.NRj-1 ; ++jBLK ) {
        for ( iBLK = 0 ; iBLK <= QuadTree.NRi-1 ; ++iBLK ) {
	   QuadTree.Roots[iBLK][jBLK].broadcast(List_of_Available_Blocks);
        } /* endfor */
    } /* endfor */

#endif

}

/**********************************************************
 * Routine: Renumber_Solution_Blocks                      *
 *                                                        *
 * Assigns a global block number to all of the solution   *
 * blocks in the quadtree adaptive block hierarchical     *
 * data structure.                                        *
 *                                                        *
 **********************************************************/
void Renumber_Solution_Blocks(QuadTreeBlock_DataStructure &QuadTree) {

    QuadTree.renumber();

}

/**********************************************************
 * Routine: Renumber_Solution_Blocks                      *
 *                                                        *
 * Assigns a global block number to all of the solution   *
 * blocks in the quadtree adaptive block hierarchical     *
 * data structure.                                        *
  *                                                        *
 **********************************************************/
void Renumber_Solution_Blocks(QuadTreeBlock_DataStructure &QuadTree,
                              AdaptiveBlock2D_List &LocalSolnBlockList) {

    int iCPU, iBLK, global_block_number;

    QuadTree.renumber();
   
    iCPU = LocalSolnBlockList.ThisCPU;
    for ( iBLK = 0 ; iBLK <= LocalSolnBlockList.Nblk-1 ; ++iBLK ) {
       if (QuadTree.Blocks[iCPU][iBLK] != NULL) {
          if (QuadTree.Blocks[iCPU][iBLK]->block.used) {
             LocalSolnBlockList.Block[iBLK].gblknum = 
                QuadTree.Blocks[iCPU][iBLK]->block.gblknum;
          } /* endif */
       } /* endif */
    } /* endfor */    

}

/**********************************************************
 * Routine: Find_Neighbours_of_Root_Solution_Blocks       *
 *                                                        *
 * Determines and stores the neighbours of all root       *
 * solution blocks in the quadtree adaptive block data    *
 * structure.  It is assumed that all root solution       *
 * blocks are initially at the same level of refinement.  *
 *                                                        *
 **********************************************************/
void Find_Neighbours_of_Root_Solution_Blocks(QuadTreeBlock_DataStructure 
                                             &QuadTree) {

    QuadTree.findRootNeighbours();

}

/**********************************************************
 * Routine: Find_Neighbours_of_Root_Solution_Blocks       *
 *                                                        *
 * Determines and stores the neighbours of all root       *
 * solution blocks in the quadtree adaptive block data    *
 * structure.  It is assumed that all root solution       *
 * blocks are initially at the same level of refinement.  *
 *                                                        *
 **********************************************************/
void Find_Neighbours_of_Root_Solution_Blocks(QuadTreeBlock_DataStructure 
                                             &QuadTree,
                                             AdaptiveBlock2D_List 
                                             &LocalSolnBlockList) {

    int iBLK, jBLK;

    QuadTree.findRootNeighbours();

    for ( jBLK = 0 ; jBLK <= QuadTree.NRj-1 ; ++jBLK ) {
        for ( iBLK = 0 ; iBLK <= QuadTree.NRi-1 ; ++iBLK ) {
           if (QuadTree.Roots[iBLK][jBLK].block.used) {
	      // Assign neighbour information to local block list.
	      if (QuadTree.Roots[iBLK][jBLK].block.nW == 1) {
                 if (LocalSolnBlockList.ThisCPU == 
                     QuadTree.Roots[iBLK][jBLK].block.info.cpu) {
                    LocalSolnBlockList.Block[QuadTree.Roots[iBLK][jBLK].block.info.blknum].nW = 1;
                    LocalSolnBlockList.Block[QuadTree.Roots[iBLK][jBLK].block.info.blknum].infoW[0] =
                       QuadTree.Roots[iBLK-1][jBLK].block.info;
                 } /* endif */
              } /* endif */
	   
	      if (QuadTree.Roots[iBLK][jBLK].block.nE == 1) {
                 if (LocalSolnBlockList.ThisCPU == 
                     QuadTree.Roots[iBLK][jBLK].block.info.cpu) {
                    LocalSolnBlockList.Block[QuadTree.Roots[iBLK][jBLK].block.info.blknum].nE = 1;
                    LocalSolnBlockList.Block[QuadTree.Roots[iBLK][jBLK].block.info.blknum].infoE[0] =
                       QuadTree.Roots[iBLK+1][jBLK].block.info;
                 } /* endif */
              } /* endif */
	   
	      if (QuadTree.Roots[iBLK][jBLK].block.nS == 1) {
                 if (LocalSolnBlockList.ThisCPU == 
                     QuadTree.Roots[iBLK][jBLK].block.info.cpu) {
                    LocalSolnBlockList.Block[QuadTree.Roots[iBLK][jBLK].block.info.blknum].nS = 1;
                    LocalSolnBlockList.Block[QuadTree.Roots[iBLK][jBLK].block.info.blknum].infoS[0] =
                       QuadTree.Roots[iBLK][jBLK-1].block.info;
                 } /* endif */
              } /* endif */
	   
	      if (QuadTree.Roots[iBLK][jBLK].block.nN == 1) {
                 if (LocalSolnBlockList.ThisCPU == 
                     QuadTree.Roots[iBLK][jBLK].block.info.cpu) {
                    LocalSolnBlockList.Block[QuadTree.Roots[iBLK][jBLK].block.info.blknum].nN = 1;
                    LocalSolnBlockList.Block[QuadTree.Roots[iBLK][jBLK].block.info.blknum].infoN[0] =
                       QuadTree.Roots[iBLK][jBLK+1].block.info;
                 } /* endif */
              } /* endif */
	   
	      if (QuadTree.Roots[iBLK][jBLK].block.nNW == 1) {
                 if (LocalSolnBlockList.ThisCPU == 
                     QuadTree.Roots[iBLK][jBLK].block.info.cpu) {
                    LocalSolnBlockList.Block[QuadTree.Roots[iBLK][jBLK].block.info.blknum].nNW = 1;
                    LocalSolnBlockList.Block[QuadTree.Roots[iBLK][jBLK].block.info.blknum].infoNW[0] =
                       QuadTree.Roots[iBLK-1][jBLK+1].block.info;
                 } /* endif */
              } /* endif */
	   
	      if (QuadTree.Roots[iBLK][jBLK].block.nSW == 1) {
                 if (LocalSolnBlockList.ThisCPU == 
                     QuadTree.Roots[iBLK][jBLK].block.info.cpu) {
                    LocalSolnBlockList.Block[QuadTree.Roots[iBLK][jBLK].block.info.blknum].nSW = 1;
                    LocalSolnBlockList.Block[QuadTree.Roots[iBLK][jBLK].block.info.blknum].infoSW[0] =
                       QuadTree.Roots[iBLK-1][jBLK-1].block.info;
                 } /* endif */
              } /* endif */
	   
	      if (QuadTree.Roots[iBLK][jBLK].block.nNE == 1) {
                 if (LocalSolnBlockList.ThisCPU == 
                     QuadTree.Roots[iBLK][jBLK].block.info.cpu) {
                    LocalSolnBlockList.Block[QuadTree.Roots[iBLK][jBLK].block.info.blknum].nNE = 1;
                    LocalSolnBlockList.Block[QuadTree.Roots[iBLK][jBLK].block.info.blknum].infoNE[0] =
                       QuadTree.Roots[iBLK+1][jBLK+1].block.info;
                 } /* endif */
              } /* endif */
	   
	      if (QuadTree.Roots[iBLK][jBLK].block.nSE == 1) {
                 if (LocalSolnBlockList.ThisCPU == 
                     QuadTree.Roots[iBLK][jBLK].block.info.cpu) {
                    LocalSolnBlockList.Block[QuadTree.Roots[iBLK][jBLK].block.info.blknum].nSE = 1;
                    LocalSolnBlockList.Block[QuadTree.Roots[iBLK][jBLK].block.info.blknum].infoSE[0] =
                       QuadTree.Roots[iBLK+1][jBLK-1].block.info;
                 } /* endif */
              } /* endif */
           } /* endif */
        } /* endfor */
    } /* endfor */

}

/**********************************************************
 * Routine: Modify_Neighbours_of_Root_Solution_Blocks     *
 *                                                        *
 * Modifies the neigbours of all root solution blocks in  *
 * the quadtree adaptive block data structure, depending  *
 * on the multi-block grid type.                          *
 *                                                        *
 **********************************************************/
void Modify_Neighbours_of_Root_Solution_Blocks(QuadTreeBlock_DataStructure 
                                               &QuadTree,
                                               const int Grid_Type) {

    switch(Grid_Type) {
      case GRID_SQUARE :
//         QuadTree.Roots[1][0].block.nN = 0;
//         QuadTree.Roots[1][0].block.infoN[0].reset();

//         QuadTree.Roots[1][1].block.nS = 0;
//         QuadTree.Roots[1][1].block.infoS[0].reset();
        break;
      case GRID_CIRCULAR_CYLINDER :
      case GRID_ELLIPSE :
        QuadTree.Roots[0][0].block.nW = 1;
        QuadTree.Roots[0][0].block.infoW[0] = QuadTree.Roots[1][0].block.info;

        QuadTree.Roots[1][0].block.nE = 1;
        QuadTree.Roots[1][0].block.infoE[0] = QuadTree.Roots[0][0].block.info; 
        break;
      case GRID_NACA_AEROFOIL :
        QuadTree.Roots[0][0].block.nS = 1;
        QuadTree.Roots[0][0].block.infoS[0] = QuadTree.Roots[3][0].block.info;
        QuadTree.Roots[0][0].block.infoS[0].dimen.i = 
           - QuadTree.Roots[0][0].block.infoS[0].dimen.i;
        QuadTree.Roots[0][0].block.infoS[0].dimen.j = 
           - QuadTree.Roots[0][0].block.infoS[0].dimen.j;

        QuadTree.Roots[3][0].block.nS = 1;
        QuadTree.Roots[3][0].block.infoS[0] = QuadTree.Roots[0][0].block.info;
        QuadTree.Roots[3][0].block.infoS[0].dimen.i = 
           - QuadTree.Roots[3][0].block.infoS[0].dimen.i;
        QuadTree.Roots[3][0].block.infoS[0].dimen.j = 
           - QuadTree.Roots[3][0].block.infoS[0].dimen.j; 
        break;
      case GRID_NASA_ROTOR_37 :
      case GRID_NASA_ROTOR_67 :
        QuadTree.Roots[0][0].block.nS = 1;
        QuadTree.Roots[0][0].block.infoS[0] = QuadTree.Roots[0][1].block.info;
        QuadTree.Roots[0][0].block.nSE = 1;
        QuadTree.Roots[0][0].block.infoSE[0] = QuadTree.Roots[1][1].block.info;
        //QuadTree.Roots[0][0].block.nNE = 0;
        //QuadTree.Roots[0][0].block.infoNE[0].reset();

        QuadTree.Roots[1][0].block.nS = 1;
        QuadTree.Roots[1][0].block.infoS[0] = QuadTree.Roots[1][1].block.info;
        QuadTree.Roots[1][0].block.nSE = 1;
        QuadTree.Roots[1][0].block.infoSE[0] = QuadTree.Roots[2][1].block.info;
        QuadTree.Roots[1][0].block.nSW = 1;
        QuadTree.Roots[1][0].block.infoSW[0] = QuadTree.Roots[0][1].block.info;
        QuadTree.Roots[1][0].block.nN = 0;
        QuadTree.Roots[1][0].block.infoN[0].reset();
        //QuadTree.Roots[1][0].block.nNW = 0;
        //QuadTree.Roots[1][0].block.infoNW[0].reset();
        //QuadTree.Roots[1][0].block.nNE = 0;
        //QuadTree.Roots[1][0].block.infoNE[0].reset();

        QuadTree.Roots[2][0].block.nS = 1;
        QuadTree.Roots[2][0].block.infoS[0] = QuadTree.Roots[2][1].block.info;
        QuadTree.Roots[2][0].block.nSW = 1;
        QuadTree.Roots[2][0].block.infoSW[0] = QuadTree.Roots[1][1].block.info;
        //QuadTree.Roots[2][0].block.nNW = 0;
        //QuadTree.Roots[2][0].block.infoNW[0].reset();

        QuadTree.Roots[0][1].block.nN = 1;
        QuadTree.Roots[0][1].block.infoN[0] = QuadTree.Roots[0][0].block.info;
        QuadTree.Roots[0][1].block.nNE = 1;
        QuadTree.Roots[0][1].block.infoNE[0] = QuadTree.Roots[1][0].block.info;
        //QuadTree.Roots[0][1].block.nSE = 0;
        //QuadTree.Roots[0][1].block.infoSE[0].reset();

        QuadTree.Roots[1][1].block.nN = 1;
        QuadTree.Roots[1][1].block.infoN[0] = QuadTree.Roots[1][0].block.info;
        QuadTree.Roots[1][1].block.nNW = 1;
        QuadTree.Roots[1][1].block.infoNW[0] = QuadTree.Roots[0][0].block.info;
        QuadTree.Roots[1][1].block.nNE = 1;
        QuadTree.Roots[1][1].block.infoNE[0] = QuadTree.Roots[2][0].block.info;
        QuadTree.Roots[1][1].block.nS = 0;
        QuadTree.Roots[1][1].block.infoS[0].reset();
        //QuadTree.Roots[1][1].block.nSE = 0;
        //QuadTree.Roots[1][1].block.infoSE[0].reset();
        //QuadTree.Roots[1][1].block.nSW = 0;
        //QuadTree.Roots[1][1].block.infoSW[0].reset();

        QuadTree.Roots[2][1].block.nN = 1;
        QuadTree.Roots[2][1].block.infoN[0] = QuadTree.Roots[2][0].block.info;
        QuadTree.Roots[2][1].block.nNW = 1;
        QuadTree.Roots[2][1].block.infoNW[0] = QuadTree.Roots[1][0].block.info;
        //QuadTree.Roots[2][1].block.nSW = 0;
        //QuadTree.Roots[2][1].block.infoSW[0].reset();
        break;
    } /* endswitch */

}

/**********************************************************
 * Routine: Modify_Neighbours_of_Root_Solution_Blocks     *
 *                                                        *
 * Modifies the neigbours of all root solution blocks in  *
 * the quadtree adaptive block data structure, depending  *
 * on the multi-block grid type.                          *
 *                                                        *
 **********************************************************/
void Modify_Neighbours_of_Root_Solution_Blocks(QuadTreeBlock_DataStructure 
                                               &QuadTree,
                                               AdaptiveBlock2D_List 
                                               &LocalSolnBlockList,
                                               const int Grid_Type) {

    Modify_Neighbours_of_Root_Solution_Blocks(QuadTree, Grid_Type);

    switch(Grid_Type) {
      case GRID_SQUARE :

        break;
      case GRID_CIRCULAR_CYLINDER :
      case GRID_ELLIPSE :
        if (LocalSolnBlockList.ThisCPU == QuadTree.Roots[0][0].block.info.cpu) {
           LocalSolnBlockList.Block[QuadTree.Roots[0][0].block.info.blknum] = 
              QuadTree.Roots[0][0].block;
        } /* endif */
        if (LocalSolnBlockList.ThisCPU == QuadTree.Roots[1][0].block.info.cpu) {
           LocalSolnBlockList.Block[QuadTree.Roots[1][0].block.info.blknum] = 
              QuadTree.Roots[1][0].block;
        } /* endif */
        break;
      case GRID_NACA_AEROFOIL :
        if (LocalSolnBlockList.ThisCPU == QuadTree.Roots[0][0].block.info.cpu) {
           LocalSolnBlockList.Block[QuadTree.Roots[0][0].block.info.blknum] = 
              QuadTree.Roots[0][0].block;
        } /* endif */
        if (LocalSolnBlockList.ThisCPU == QuadTree.Roots[3][0].block.info.cpu) {
           LocalSolnBlockList.Block[QuadTree.Roots[3][0].block.info.blknum] = 
              QuadTree.Roots[3][0].block;
        } /* endif */
        break;
      case GRID_NASA_ROTOR_37 :
      case GRID_NASA_ROTOR_67 :
        if (LocalSolnBlockList.ThisCPU == QuadTree.Roots[0][0].block.info.cpu) {
           LocalSolnBlockList.Block[QuadTree.Roots[0][0].block.info.blknum] = 
              QuadTree.Roots[0][0].block;
        } /* endif */
        if (LocalSolnBlockList.ThisCPU == QuadTree.Roots[1][0].block.info.cpu) {
           LocalSolnBlockList.Block[QuadTree.Roots[1][0].block.info.blknum] = 
              QuadTree.Roots[1][0].block;
        } /* endif */
        if (LocalSolnBlockList.ThisCPU == QuadTree.Roots[2][0].block.info.cpu) {
           LocalSolnBlockList.Block[QuadTree.Roots[2][0].block.info.blknum] = 
              QuadTree.Roots[2][0].block;
        } /* endif */
        if (LocalSolnBlockList.ThisCPU == QuadTree.Roots[0][1].block.info.cpu) {
           LocalSolnBlockList.Block[QuadTree.Roots[0][1].block.info.blknum] = 
              QuadTree.Roots[0][1].block;
        } /* endif */
        if (LocalSolnBlockList.ThisCPU == QuadTree.Roots[1][1].block.info.cpu) {
           LocalSolnBlockList.Block[QuadTree.Roots[1][1].block.info.blknum] = 
              QuadTree.Roots[1][1].block;
        } /* endif */
        if (LocalSolnBlockList.ThisCPU == QuadTree.Roots[2][1].block.info.cpu) {
           LocalSolnBlockList.Block[QuadTree.Roots[2][1].block.info.blknum] = 
              QuadTree.Roots[2][1].block;
        } /* endif */
        break;
    } /* endswitch */

}

/**********************************************************
 * Routine: Find_Neighbours                               *
 *                                                        *
 * Determines and stores the neigbours of all solution    *
 * blocks in the quadtree adaptive block data structure.  *
 *                                                        *
 **********************************************************/
void Find_Neighbours(QuadTreeBlock_DataStructure &QuadTree) {

    QuadTree.findNeighbours();

}

/**********************************************************
 * Routine: Find_Neighbours                               *
 *                                                        *
 * Determines and stores the neigbours of all solution    *
 * blocks in the quadtree adaptive block data structure.  *
 *                                                        *
 **********************************************************/
void Find_Neighbours(QuadTreeBlock_DataStructure &QuadTree,
                     AdaptiveBlock2D_List &LocalSolnBlockList) {

    int iCPU, iBLK, global_block_number;

    QuadTree.findNeighbours();
   
    iCPU = LocalSolnBlockList.ThisCPU;
    for ( iBLK = 0 ; iBLK <= LocalSolnBlockList.Nblk-1 ; ++iBLK ) {
       if (QuadTree.Blocks[iCPU][iBLK] != NULL) {
          if (QuadTree.Blocks[iCPU][iBLK]->block.used) {
             LocalSolnBlockList.Block[iBLK] = 
                QuadTree.Blocks[iCPU][iBLK]->block;
          } /* endif */
       } /* endif */
    } /* endfor */

}

/**********************************************************
 * Routine: Get_Refinement_List                           *
 *                                                        *
 * Obtains the a global list of all solution blocks in    *
 * the quadtree adaptive block hierarchical data          *
 * structure scheduled for mesh refinement or coarsening. *
 *                                                        *
 **********************************************************/
void Get_Refinement_List(QuadTreeBlock_DataStructure &QuadTree,
                         AdaptiveBlock2D_List &LocalSolnBlockList) {

    int i_neighbour, CPU_neighbour, BLK_neighbour;
    int i_sibling, CPU_sibling, BLK_sibling, number_of_coarsened_sibling;

    int number_of_passes, number_of_changes;
    
    /* Set the quadtree refinement flags to default values
       (no change or no division and coarsening). */

    QuadTree.nochangeAll();

    /* Collect the local refinement flags for each of the
       local solution blocks and assign them to the global
       quadtree list. */

#ifdef _MPI_VERSION
    int iCPU, iBLK, iTotal;
    int  *receive_buffer;

    receive_buffer = new int[QuadTree.Ncpu*QuadTree.Nblk];

    MPI::COMM_WORLD.Allgather(LocalSolnBlockList.RefineFlag,
                              LocalSolnBlockList.Nblk,
                              MPI::INT,
                              receive_buffer,
                              QuadTree.Nblk,
                              MPI::INT);
    
    iTotal = 0;
    for ( iCPU = 0 ; iCPU <= QuadTree.Ncpu-1 ; ++iCPU ) {
       for ( iBLK = 0 ; iBLK <= QuadTree.Nblk-1 ; ++iBLK ) {           
	   QuadTree.RefineFlags[iCPU][iBLK] = receive_buffer[iTotal];
	   iTotal += 1;
       } /* endfor */
    } /* endfor */
 
    delete []receive_buffer;
    receive_buffer = NULL;
#else
    int iCPU, iBLK;

    for ( iCPU = 0 ; iCPU <= QuadTree.Ncpu-1 ; ++iCPU ) {
       for ( iBLK = 0 ; iBLK <= QuadTree.Nblk-1 ; ++iBLK ) {
	   QuadTree.RefineFlags[iCPU][iBLK] = LocalSolnBlockList.RefineFlag[iBLK];
        } /* endfor */
    } /* endfor */
#endif

    /* Check and adjust the global refinement flags as necessary to 
       ensure that block divsion and coarsening is possible and the
       resolution change between adjacent blocks is at most a factor
       of two. */

    for ( number_of_passes = 1; number_of_passes <= 10; ++number_of_passes ) {
    number_of_changes = 0;
    // Check refinement
    for ( iCPU = 0 ; iCPU <= QuadTree.Ncpu-1 ; ++iCPU ) {
       for ( iBLK = 0 ; iBLK <= QuadTree.Nblk-1 ; ++iBLK ) {
          if (QuadTree.Blocks[iCPU][iBLK] != NULL) {
             if (!QuadTree.Blocks[iCPU][iBLK]->block.used) {
	        QuadTree.RefineFlags[iCPU][iBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
             } else if (QuadTree.RefineFlags[iCPU][iBLK] == ADAPTIVEBLOCK2D_REFINE) {
	        // Check each of the neighbouring blocks to see if block division is possible.
                // North Neighbour(s):
                if ((QuadTree.RefineFlags[iCPU][iBLK] == ADAPTIVEBLOCK2D_REFINE) &&
                    (QuadTree.Blocks[iCPU][iBLK]->block.nN > 0)) { 
	           for (i_neighbour = 0 ; i_neighbour <= QuadTree.Blocks[iCPU][iBLK]->block.nN-1 ; ++i_neighbour) {
                      CPU_neighbour = QuadTree.Blocks[iCPU][iBLK]->block.infoN[i_neighbour].cpu;
                      BLK_neighbour = QuadTree.Blocks[iCPU][iBLK]->block.infoN[i_neighbour].blknum;
		      if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                           ADAPTIVEBLOCK2D_NOCHANGE) &&
                          ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                            QuadTree.Blocks[iCPU][iBLK]->block.infoN[i_neighbour].level) <= 0) &&
                          ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                            QuadTree.Blocks[iCPU][iBLK]->block.infoN[i_neighbour].level) >= -1) ) {
                         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_NOCHANGE) &&
                                ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                                  QuadTree.Blocks[iCPU][iBLK]->block.infoN[i_neighbour].level) == +1) ) {
		         // Force refinement of neighbour.
		         QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] = ADAPTIVEBLOCK2D_REFINE;
                         number_of_changes += 1;
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_REFINE) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoN[i_neighbour].level) <= +1) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoN[i_neighbour].level) >= -1)) {
		         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_COARSEN) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                  QuadTree.Blocks[iCPU][iBLK]->block.infoN[i_neighbour].level) == -1)) {
		         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_COARSEN) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoN[i_neighbour].level) >= 0) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoN[i_neighbour].level) <= +1) ) {
		         // Prevent coarsening of neighbour.
		         QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] = ADAPTIVEBLOCK2D_NOCHANGE;
                         number_of_changes += 1;
                      } else {
		         // Prevent refinement of solution block.
		         QuadTree.RefineFlags[iCPU][iBLK] == ADAPTIVEBLOCK2D_NOCHANGE;
                         number_of_changes += 1;
                      } /* endif */
                   } /* endfor */
                } /* endif */

                // South Neighbour(s):
	        if ((QuadTree.RefineFlags[iCPU][iBLK] == ADAPTIVEBLOCK2D_REFINE) &&
                    (QuadTree.Blocks[iCPU][iBLK]->block.nS > 0)) { 
	           for (i_neighbour = 0 ; i_neighbour <= QuadTree.Blocks[iCPU][iBLK]->block.nS-1 ; ++i_neighbour) {
                      CPU_neighbour = QuadTree.Blocks[iCPU][iBLK]->block.infoS[i_neighbour].cpu;
                      BLK_neighbour = QuadTree.Blocks[iCPU][iBLK]->block.infoS[i_neighbour].blknum;
		      if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                           ADAPTIVEBLOCK2D_NOCHANGE) &&
                          ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                            QuadTree.Blocks[iCPU][iBLK]->block.infoS[i_neighbour].level) <= 0) &&
                          ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                            QuadTree.Blocks[iCPU][iBLK]->block.infoS[i_neighbour].level) >= -1) ) {
		         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_NOCHANGE) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoS[i_neighbour].level) == +1) ) {
		         // Force refinement of neighbour.
		         QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] = ADAPTIVEBLOCK2D_REFINE;
                         number_of_changes += 1;
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_REFINE) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoS[i_neighbour].level) <= +1) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoS[i_neighbour].level) >= -1)) {
		         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_COARSEN) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoS[i_neighbour].level) == -1)) {
		         // Okay, do nothing!
	   	      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_COARSEN) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoS[i_neighbour].level) >= 0) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoS[i_neighbour].level) <= +1) ) {
		         // Prevent coarsening of neighbour.
		         QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] = ADAPTIVEBLOCK2D_NOCHANGE;
                         number_of_changes += 1;
                      } else {
		         // Prevent refinement of solution block.
		         QuadTree.RefineFlags[iCPU][iBLK] == ADAPTIVEBLOCK2D_NOCHANGE;
                         number_of_changes += 1;
                      } /* endif */
                   } /* endfor */
                } /* endif */

                // East Neighbour(s):
	        if ((QuadTree.RefineFlags[iCPU][iBLK] == ADAPTIVEBLOCK2D_REFINE) &&
                    (QuadTree.Blocks[iCPU][iBLK]->block.nE > 0)) {
	           for (i_neighbour = 0 ; i_neighbour <= QuadTree.Blocks[iCPU][iBLK]->block.nE-1 ; ++i_neighbour) {
                      CPU_neighbour = QuadTree.Blocks[iCPU][iBLK]->block.infoE[i_neighbour].cpu;
                      BLK_neighbour = QuadTree.Blocks[iCPU][iBLK]->block.infoE[i_neighbour].blknum;
		      if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                           ADAPTIVEBLOCK2D_NOCHANGE) &&
                          ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                            QuadTree.Blocks[iCPU][iBLK]->block.infoE[i_neighbour].level) <= 0) &&
                          ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                            QuadTree.Blocks[iCPU][iBLK]->block.infoE[i_neighbour].level) >= -1) ) {
		         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_NOCHANGE) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoE[i_neighbour].level) == +1) ) {
		         // Force refinement of neighbour.
		         QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] = ADAPTIVEBLOCK2D_REFINE;
                         number_of_changes += 1;
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_REFINE) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoE[i_neighbour].level) <= +1) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoE[i_neighbour].level) >= -1)) {
		         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_COARSEN) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoE[i_neighbour].level) == -1)) {
		         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_COARSEN) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoE[i_neighbour].level) >= 0) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoE[i_neighbour].level) <= +1) ) {
		         // Prevent coarsening of neighbour.
		         QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] = ADAPTIVEBLOCK2D_NOCHANGE;
                         number_of_changes += 1;
                      } else {
		         // Prevent refinement of solution block.
		         QuadTree.RefineFlags[iCPU][iBLK] == ADAPTIVEBLOCK2D_NOCHANGE;
                         number_of_changes += 1;
                      } /* endif */
                   } /* endfor */
                } /* endif */

                // West Neighbour(s):
	        if ((QuadTree.RefineFlags[iCPU][iBLK] == ADAPTIVEBLOCK2D_REFINE) &&
                    (QuadTree.Blocks[iCPU][iBLK]->block.nW > 0)) {
	           for (i_neighbour = 0 ; i_neighbour <= QuadTree.Blocks[iCPU][iBLK]->block.nW-1 ; ++i_neighbour) {
                      CPU_neighbour = QuadTree.Blocks[iCPU][iBLK]->block.infoW[i_neighbour].cpu;
                      BLK_neighbour = QuadTree.Blocks[iCPU][iBLK]->block.infoW[i_neighbour].blknum;
		      if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                           ADAPTIVEBLOCK2D_NOCHANGE) &&
                          ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                            QuadTree.Blocks[iCPU][iBLK]->block.infoW[i_neighbour].level) <= 0) &&
                          ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                            QuadTree.Blocks[iCPU][iBLK]->block.infoW[i_neighbour].level) >= -1) ) {
		         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_NOCHANGE) &&
                                ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                                  QuadTree.Blocks[iCPU][iBLK]->block.infoW[i_neighbour].level) == +1) ) {
		         // Force refinement of neighbour.
		         QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] = ADAPTIVEBLOCK2D_REFINE;
                         number_of_changes += 1;
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_REFINE) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoW[i_neighbour].level) <= +1) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoW[i_neighbour].level) >= -1)) {
		         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_COARSEN) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoW[i_neighbour].level) == -1)) {
		         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_COARSEN) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoW[i_neighbour].level) >= 0) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoW[i_neighbour].level) <= +1) ) {
		         // Prevent coarsening of neighbour.
		         QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] = ADAPTIVEBLOCK2D_NOCHANGE;
                         number_of_changes += 1;
                      } else {
		         // Prevent refinement of solution block.
		         QuadTree.RefineFlags[iCPU][iBLK] == ADAPTIVEBLOCK2D_NOCHANGE;
                         number_of_changes += 1;
                      } /* endif */
                   } /* endfor */
                } /* endif */

                // North-West Neighbour(s):
	        if ((QuadTree.RefineFlags[iCPU][iBLK] == ADAPTIVEBLOCK2D_REFINE) &&
                    (QuadTree.Blocks[iCPU][iBLK]->block.nNW > 0)) {
	           for (i_neighbour = 0 ; i_neighbour <= QuadTree.Blocks[iCPU][iBLK]->block.nNW-1 ; ++i_neighbour) {
                      CPU_neighbour = QuadTree.Blocks[iCPU][iBLK]->block.infoNW[i_neighbour].cpu;
                      BLK_neighbour = QuadTree.Blocks[iCPU][iBLK]->block.infoNW[i_neighbour].blknum;
		      if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                           ADAPTIVEBLOCK2D_NOCHANGE) &&
                          ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                            QuadTree.Blocks[iCPU][iBLK]->block.infoNW[i_neighbour].level) <= 0) &&
                          ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                            QuadTree.Blocks[iCPU][iBLK]->block.infoNW[i_neighbour].level) >= -1) ) {
		         // Okay, do nothing!
 	   	      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_NOCHANGE) &&
                                ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                                  QuadTree.Blocks[iCPU][iBLK]->block.infoNW[i_neighbour].level) == +1) ) {
		         // Force refinement of neighbour.
		         QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] = ADAPTIVEBLOCK2D_REFINE;
                         number_of_changes += 1;
	   	      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_REFINE) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoNW[i_neighbour].level) <= +1) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoNW[i_neighbour].level) >= -1)) {
		         // Okay, do nothing!
	 	      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_COARSEN) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoNW[i_neighbour].level) == -1)) {
		         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_COARSEN) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoNW[i_neighbour].level) >= 0) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoNW[i_neighbour].level) <= +1) ) {
		         // Prevent coarsening of neighbour.
		         QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] = ADAPTIVEBLOCK2D_NOCHANGE;
                         number_of_changes += 1;
                      } else {
		         // Prevent refinement of solution block.
		         QuadTree.RefineFlags[iCPU][iBLK] == ADAPTIVEBLOCK2D_NOCHANGE;
                         number_of_changes += 1;
                      } /* endif */
                   } /* endfor */
                } /* endif */

                // North-East Neighbour(s):
	        if ((QuadTree.RefineFlags[iCPU][iBLK] == ADAPTIVEBLOCK2D_REFINE) &&
                    (QuadTree.Blocks[iCPU][iBLK]->block.nNE > 0)) {
	           for (i_neighbour = 0 ; i_neighbour <= QuadTree.Blocks[iCPU][iBLK]->block.nNE-1 ; ++i_neighbour) {
                      CPU_neighbour = QuadTree.Blocks[iCPU][iBLK]->block.infoNE[i_neighbour].cpu;
                      BLK_neighbour = QuadTree.Blocks[iCPU][iBLK]->block.infoNE[i_neighbour].blknum;
		      if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                           ADAPTIVEBLOCK2D_NOCHANGE) &&
                          ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                            QuadTree.Blocks[iCPU][iBLK]->block.infoNE[i_neighbour].level) <= 0) &&
                          ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                            QuadTree.Blocks[iCPU][iBLK]->block.infoNE[i_neighbour].level) >= -1) ) {
		         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_NOCHANGE) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoNE[i_neighbour].level) == +1) ) {
		         // Force refinement of neighbour.
		         QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] = ADAPTIVEBLOCK2D_REFINE;
                         number_of_changes += 1;
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_REFINE) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoNE[i_neighbour].level) <= +1) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoNE[i_neighbour].level) >= -1)) {
		         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_COARSEN) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoNE[i_neighbour].level) == -1)) {
		         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_COARSEN) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoNE[i_neighbour].level) >= 0) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoNE[i_neighbour].level) <= +1) ) {
		         // Prevent coarsening of neighbour.
		         QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] = ADAPTIVEBLOCK2D_NOCHANGE;
                         number_of_changes += 1;
                      } else {
		         // Prevent refinement of solution block.
		         QuadTree.RefineFlags[iCPU][iBLK] == ADAPTIVEBLOCK2D_NOCHANGE;
                         number_of_changes += 1;
                      } /* endif */
                   } /* endfor */
                } /* endif */

                // South-East Neighbour(s):
	        if ((QuadTree.RefineFlags[iCPU][iBLK] == ADAPTIVEBLOCK2D_REFINE) &&
                    (QuadTree.Blocks[iCPU][iBLK]->block.nSE > 0)) {
	           for (i_neighbour = 0 ; i_neighbour <= QuadTree.Blocks[iCPU][iBLK]->block.nSE-1 ; ++i_neighbour) {
                      CPU_neighbour = QuadTree.Blocks[iCPU][iBLK]->block.infoSE[i_neighbour].cpu;
                      BLK_neighbour = QuadTree.Blocks[iCPU][iBLK]->block.infoSE[i_neighbour].blknum;
	   	      if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                           ADAPTIVEBLOCK2D_NOCHANGE) &&
                          ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                            QuadTree.Blocks[iCPU][iBLK]->block.infoSE[i_neighbour].level) <= 0) &&
                          ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                            QuadTree.Blocks[iCPU][iBLK]->block.infoSE[i_neighbour].level) >= -1) ) {
		         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_NOCHANGE) &&
                                ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                                  QuadTree.Blocks[iCPU][iBLK]->block.infoSE[i_neighbour].level) == +1) ) {
		         // Force refinement of neighbour.
		         QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] = ADAPTIVEBLOCK2D_REFINE;
                         number_of_changes += 1;
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_REFINE) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoSE[i_neighbour].level) <= +1) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoSE[i_neighbour].level) >= -1)) {
		         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_COARSEN) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoSE[i_neighbour].level) == -1)) {
		         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_COARSEN) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoSE[i_neighbour].level) >= 0) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoSE[i_neighbour].level) <= +1) ) {
		         // Prevent coarsening of neighbour.
		         QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] = ADAPTIVEBLOCK2D_NOCHANGE;
                         number_of_changes += 1;
                      } else {
		         // Prevent refinement of solution block.
		         QuadTree.RefineFlags[iCPU][iBLK] == ADAPTIVEBLOCK2D_NOCHANGE;
                         number_of_changes += 1;
                      } /* endif */
                   } /* endfor */
                } /* endif */

                // South-West Neighbour(s):
	        if ((QuadTree.RefineFlags[iCPU][iBLK] == ADAPTIVEBLOCK2D_REFINE) &&
                    (QuadTree.Blocks[iCPU][iBLK]->block.nSW > 0)) {
	           for (i_neighbour = 0 ; i_neighbour <= QuadTree.Blocks[iCPU][iBLK]->block.nSW-1 ; ++i_neighbour) {
                      CPU_neighbour = QuadTree.Blocks[iCPU][iBLK]->block.infoSW[i_neighbour].cpu;
                      BLK_neighbour = QuadTree.Blocks[iCPU][iBLK]->block.infoSW[i_neighbour].blknum;
		      if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                           ADAPTIVEBLOCK2D_NOCHANGE) &&
                          ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                            QuadTree.Blocks[iCPU][iBLK]->block.infoSW[i_neighbour].level) <= 0) &&
                          ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                            QuadTree.Blocks[iCPU][iBLK]->block.infoSW[i_neighbour].level) >= -1) ) {
		         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_NOCHANGE) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoSW[i_neighbour].level) == +1) ) {
		         // Force refinement of neighbour.
		         QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] = ADAPTIVEBLOCK2D_REFINE;
                         number_of_changes += 1;
	   	      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_REFINE) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoSW[i_neighbour].level) <= +1) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoSW[i_neighbour].level) >= -1)) {
		         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_COARSEN) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoSW[i_neighbour].level) == -1)) {
		         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_COARSEN) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoSW[i_neighbour].level) >= 0) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoSW[i_neighbour].level) <= +1) ) {
		         // Prevent coarsening of neighbour.
		         QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] = ADAPTIVEBLOCK2D_NOCHANGE;
                         number_of_changes += 1;
                      } else {
		         // Prevent refinement of solution block.
		         QuadTree.RefineFlags[iCPU][iBLK] == ADAPTIVEBLOCK2D_NOCHANGE;
                         number_of_changes += 1;
                      } /* endif */
                   } /* endfor */
                } /* endif */

             } /* endif */

          } /* endif */
       } /* endfor */
    } /* endfor */

    // Check coarsening
    for ( iCPU = 0 ; iCPU <= QuadTree.Ncpu-1 ; ++iCPU ) {
       for ( iBLK = 0 ; iBLK <= QuadTree.Nblk-1 ; ++iBLK ) {
          if (QuadTree.Blocks[iCPU][iBLK] != NULL) {
             if (!QuadTree.Blocks[iCPU][iBLK]->block.used) {
	        QuadTree.RefineFlags[iCPU][iBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
 	     } else if (QuadTree.RefineFlags[iCPU][iBLK] == ADAPTIVEBLOCK2D_COARSEN) {
	        // If block is a root block, prevent coarsening of solution block.
                if ((QuadTree.RefineFlags[iCPU][iBLK] == ADAPTIVEBLOCK2D_COARSEN) &&
                    ((QuadTree.Blocks[iCPU][iBLK]->block.info.level == 0) ||
                     (QuadTree.Blocks[iCPU][iBLK]->parent_ptr == NULL))) {
		   QuadTree.RefineFlags[iCPU][iBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
                   number_of_changes += 1;
                } /* end if */

	        // Check each of the neighbouring blocks to see if block coarsening is possible.
                // North Neighbour(s):
                if ((QuadTree.RefineFlags[iCPU][iBLK] == ADAPTIVEBLOCK2D_COARSEN) &&
                    (QuadTree.Blocks[iCPU][iBLK]->block.nN > 0)) { 
	           for (i_neighbour = 0 ; i_neighbour <= QuadTree.Blocks[iCPU][iBLK]->block.nN-1 ; ++i_neighbour) {
                      CPU_neighbour = QuadTree.Blocks[iCPU][iBLK]->block.infoN[i_neighbour].cpu;
                      BLK_neighbour = QuadTree.Blocks[iCPU][iBLK]->block.infoN[i_neighbour].blknum;
		      if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                           ADAPTIVEBLOCK2D_NOCHANGE) &&
                          ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                            QuadTree.Blocks[iCPU][iBLK]->block.infoN[i_neighbour].level) >= 0) &&
                          ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                            QuadTree.Blocks[iCPU][iBLK]->block.infoN[i_neighbour].level) <= +1) ) {
                         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_NOCHANGE) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoN[i_neighbour].level) == -1) ) {
		         // Prevent coarsening of solution block.
		         QuadTree.RefineFlags[iCPU][iBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
                         number_of_changes += 1;
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_COARSEN) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoN[i_neighbour].level) <= +1) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoN[i_neighbour].level) >= -1)) {
		         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_REFINE) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoN[i_neighbour].level) == +1)) {
		         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_REFINE) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoN[i_neighbour].level) <= 0)) {
		         // Prevent coarsening of solution block.
		         QuadTree.RefineFlags[iCPU][iBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
                         number_of_changes += 1;
                      } else {
		         // Prevent coarsening of solution block.
		         QuadTree.RefineFlags[iCPU][iBLK] == ADAPTIVEBLOCK2D_NOCHANGE;
                         number_of_changes += 1;
                      } /* endif */
                   } /* endfor */
                } /* endif */

                // South Neighbour(s):
                if ((QuadTree.RefineFlags[iCPU][iBLK] == ADAPTIVEBLOCK2D_COARSEN) &&
                    (QuadTree.Blocks[iCPU][iBLK]->block.nS > 0)) { 
	           for (i_neighbour = 0 ; i_neighbour <= QuadTree.Blocks[iCPU][iBLK]->block.nS-1 ; ++i_neighbour) {
                      CPU_neighbour = QuadTree.Blocks[iCPU][iBLK]->block.infoS[i_neighbour].cpu;
                      BLK_neighbour = QuadTree.Blocks[iCPU][iBLK]->block.infoS[i_neighbour].blknum;
		      if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                           ADAPTIVEBLOCK2D_NOCHANGE) &&
                          ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                            QuadTree.Blocks[iCPU][iBLK]->block.infoS[i_neighbour].level) >= 0) &&
                          ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                            QuadTree.Blocks[iCPU][iBLK]->block.infoS[i_neighbour].level) <= +1) ) {
                         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_NOCHANGE) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoS[i_neighbour].level) == -1) ) {
		         // Prevent coarsening of solution block.
		         QuadTree.RefineFlags[iCPU][iBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
                         number_of_changes += 1;
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_COARSEN) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoS[i_neighbour].level) <= +1) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoS[i_neighbour].level) >= -1)) {
		         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_REFINE) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoS[i_neighbour].level) == +1)) {
		         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_REFINE) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoS[i_neighbour].level) <= 0)) {
		         // Prevent coarsening of solution block.
		         QuadTree.RefineFlags[iCPU][iBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
                         number_of_changes += 1;
                      } else {
		         // Prevent coarsening of solution block.
		         QuadTree.RefineFlags[iCPU][iBLK] == ADAPTIVEBLOCK2D_NOCHANGE;
                         number_of_changes += 1;
                      } /* endif */
                   } /* endfor */
                } /* endif */

                // East Neighbour(s):
                if ((QuadTree.RefineFlags[iCPU][iBLK] == ADAPTIVEBLOCK2D_COARSEN) &&
                    (QuadTree.Blocks[iCPU][iBLK]->block.nE > 0)) { 
	           for (i_neighbour = 0 ; i_neighbour <= QuadTree.Blocks[iCPU][iBLK]->block.nE-1 ; ++i_neighbour) {
                      CPU_neighbour = QuadTree.Blocks[iCPU][iBLK]->block.infoE[i_neighbour].cpu;
                      BLK_neighbour = QuadTree.Blocks[iCPU][iBLK]->block.infoE[i_neighbour].blknum;
		      if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                           ADAPTIVEBLOCK2D_NOCHANGE) &&
                          ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                            QuadTree.Blocks[iCPU][iBLK]->block.infoE[i_neighbour].level) >= 0) &&
                          ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                            QuadTree.Blocks[iCPU][iBLK]->block.infoE[i_neighbour].level) <= +1) ) {
                         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_NOCHANGE) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoE[i_neighbour].level) == -1) ) {
		         // Prevent coarsening of solution block.
		         QuadTree.RefineFlags[iCPU][iBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
                         number_of_changes += 1;
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_COARSEN) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoE[i_neighbour].level) <= +1) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoE[i_neighbour].level) >= -1)) {
		         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_REFINE) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoE[i_neighbour].level) == +1)) {
		         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_REFINE) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoE[i_neighbour].level) <= 0)) {
		         // Prevent coarsening of solution block.
		         QuadTree.RefineFlags[iCPU][iBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
                         number_of_changes += 1;
                      } else {
		         // Prevent coarsening of solution block.
		         QuadTree.RefineFlags[iCPU][iBLK] == ADAPTIVEBLOCK2D_NOCHANGE;
                         number_of_changes += 1;
                      } /* endif */
                   } /* endfor */
                } /* endif */

                // West Neighbour(s):
                if ((QuadTree.RefineFlags[iCPU][iBLK] == ADAPTIVEBLOCK2D_COARSEN) &&
                    (QuadTree.Blocks[iCPU][iBLK]->block.nW > 0)) { 
	           for (i_neighbour = 0 ; i_neighbour <= QuadTree.Blocks[iCPU][iBLK]->block.nW-1 ; ++i_neighbour) {
                      CPU_neighbour = QuadTree.Blocks[iCPU][iBLK]->block.infoW[i_neighbour].cpu;
                      BLK_neighbour = QuadTree.Blocks[iCPU][iBLK]->block.infoW[i_neighbour].blknum;
		      if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                           ADAPTIVEBLOCK2D_NOCHANGE) &&
                          ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                            QuadTree.Blocks[iCPU][iBLK]->block.infoW[i_neighbour].level) >= 0) &&
                          ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                            QuadTree.Blocks[iCPU][iBLK]->block.infoW[i_neighbour].level) <= +1) ) {
                         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_NOCHANGE) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoW[i_neighbour].level) == -1) ) {
		         // Prevent coarsening of solution block.
		         QuadTree.RefineFlags[iCPU][iBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
                         number_of_changes += 1;
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_COARSEN) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoW[i_neighbour].level) <= +1) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoW[i_neighbour].level) >= -1)) {
		         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_REFINE) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoW[i_neighbour].level) == +1)) {
		         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_REFINE) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoW[i_neighbour].level) <= 0)) {
		         // Prevent coarsening of solution block.
		         QuadTree.RefineFlags[iCPU][iBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
                         number_of_changes += 1;
                      } else {
		         // Prevent coarsening of solution block.
		         QuadTree.RefineFlags[iCPU][iBLK] == ADAPTIVEBLOCK2D_NOCHANGE;
                         number_of_changes += 1;
                      } /* endif */
                   } /* endfor */
                } /* endif */

                // North-West Neighbour(s):
                if ((QuadTree.RefineFlags[iCPU][iBLK] == ADAPTIVEBLOCK2D_COARSEN) &&
                    (QuadTree.Blocks[iCPU][iBLK]->block.nNW > 0)) { 
	           for (i_neighbour = 0 ; i_neighbour <= QuadTree.Blocks[iCPU][iBLK]->block.nNW-1 ; ++i_neighbour) {
                      CPU_neighbour = QuadTree.Blocks[iCPU][iBLK]->block.infoNW[i_neighbour].cpu;
                      BLK_neighbour = QuadTree.Blocks[iCPU][iBLK]->block.infoNW[i_neighbour].blknum;
		      if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                           ADAPTIVEBLOCK2D_NOCHANGE) &&
                          ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                            QuadTree.Blocks[iCPU][iBLK]->block.infoNW[i_neighbour].level) >= 0) &&
                          ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                            QuadTree.Blocks[iCPU][iBLK]->block.infoNW[i_neighbour].level) <= +1) ) {
                         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_NOCHANGE) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoNW[i_neighbour].level) == -1) ) {
		         // Prevent coarsening of solution block.
		         QuadTree.RefineFlags[iCPU][iBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
                         number_of_changes += 1;
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_COARSEN) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoNW[i_neighbour].level) <= +1) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoNW[i_neighbour].level) >= -1)) {
		         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_REFINE) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoNW[i_neighbour].level) == +1)) {
		         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_REFINE) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoNW[i_neighbour].level) <= 0)) {
		         // Prevent coarsening of solution block.
		         QuadTree.RefineFlags[iCPU][iBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
                         number_of_changes += 1;
                      } else {
		         // Prevent coarsening of solution block.
		         QuadTree.RefineFlags[iCPU][iBLK] == ADAPTIVEBLOCK2D_NOCHANGE;
                         number_of_changes += 1;
                      } /* endif */
                   } /* endfor */
                } /* endif */

                // North-East Neighbour(s):
                if ((QuadTree.RefineFlags[iCPU][iBLK] == ADAPTIVEBLOCK2D_COARSEN) &&
                    (QuadTree.Blocks[iCPU][iBLK]->block.nNE > 0)) { 
	           for (i_neighbour = 0 ; i_neighbour <= QuadTree.Blocks[iCPU][iBLK]->block.nNE-1 ; ++i_neighbour) {
                      CPU_neighbour = QuadTree.Blocks[iCPU][iBLK]->block.infoNE[i_neighbour].cpu;
                      BLK_neighbour = QuadTree.Blocks[iCPU][iBLK]->block.infoNE[i_neighbour].blknum;
		      if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                           ADAPTIVEBLOCK2D_NOCHANGE) &&
                          ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                            QuadTree.Blocks[iCPU][iBLK]->block.infoNE[i_neighbour].level) >= 0) &&
                          ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                            QuadTree.Blocks[iCPU][iBLK]->block.infoNE[i_neighbour].level) <= +1) ) {
                         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_NOCHANGE) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoNE[i_neighbour].level) == -1) ) {
		         // Prevent coarsening of solution block.
		         QuadTree.RefineFlags[iCPU][iBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
                         number_of_changes += 1;
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_COARSEN) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoNE[i_neighbour].level) <= +1) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoNE[i_neighbour].level) >= -1)) {
		         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_REFINE) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoNE[i_neighbour].level) == +1)) {
		         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_REFINE) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoNE[i_neighbour].level) <= 0)) {
		         // Prevent coarsening of solution block.
		         QuadTree.RefineFlags[iCPU][iBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
                         number_of_changes += 1;
                      } else {
		         // Prevent coarsening of solution block.
		         QuadTree.RefineFlags[iCPU][iBLK] == ADAPTIVEBLOCK2D_NOCHANGE;
                         number_of_changes += 1;
                      } /* endif */
                   } /* endfor */
                } /* endif */

                // South-East Neighbour(s):
                if ((QuadTree.RefineFlags[iCPU][iBLK] == ADAPTIVEBLOCK2D_COARSEN) &&
                    (QuadTree.Blocks[iCPU][iBLK]->block.nSE > 0)) { 
	           for (i_neighbour = 0 ; i_neighbour <= QuadTree.Blocks[iCPU][iBLK]->block.nSE-1 ; ++i_neighbour) {
                      CPU_neighbour = QuadTree.Blocks[iCPU][iBLK]->block.infoSE[i_neighbour].cpu;
                      BLK_neighbour = QuadTree.Blocks[iCPU][iBLK]->block.infoSE[i_neighbour].blknum;
		      if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                           ADAPTIVEBLOCK2D_NOCHANGE) &&
                          ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                            QuadTree.Blocks[iCPU][iBLK]->block.infoSE[i_neighbour].level) >= 0) &&
                          ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                            QuadTree.Blocks[iCPU][iBLK]->block.infoSE[i_neighbour].level) <= +1) ) {
                         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_NOCHANGE) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoSE[i_neighbour].level) == -1) ) {
		         // Prevent coarsening of solution block.
		         QuadTree.RefineFlags[iCPU][iBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
                         number_of_changes += 1;
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_COARSEN) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoSE[i_neighbour].level) <= +1) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoSE[i_neighbour].level) >= -1)) {
		         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_REFINE) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoSE[i_neighbour].level) == +1)) {
		         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_REFINE) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoSE[i_neighbour].level) <= 0)) {
		         // Prevent coarsening of solution block.
		         QuadTree.RefineFlags[iCPU][iBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
                         number_of_changes += 1;
                      } else {
		         // Prevent coarsening of solution block.
		         QuadTree.RefineFlags[iCPU][iBLK] == ADAPTIVEBLOCK2D_NOCHANGE;
                         number_of_changes += 1;
                      } /* endif */
                   } /* endfor */
                } /* endif */

                // South-West Neighbour(s):
                if ((QuadTree.RefineFlags[iCPU][iBLK] == ADAPTIVEBLOCK2D_COARSEN) &&
                    (QuadTree.Blocks[iCPU][iBLK]->block.nSW > 0)) { 
	           for (i_neighbour = 0 ; i_neighbour <= QuadTree.Blocks[iCPU][iBLK]->block.nSW-1 ; ++i_neighbour) {
                      CPU_neighbour = QuadTree.Blocks[iCPU][iBLK]->block.infoSW[i_neighbour].cpu;
                      BLK_neighbour = QuadTree.Blocks[iCPU][iBLK]->block.infoSW[i_neighbour].blknum;
		      if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                           ADAPTIVEBLOCK2D_NOCHANGE) &&
                          ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                            QuadTree.Blocks[iCPU][iBLK]->block.infoSW[i_neighbour].level) >= 0) &&
                          ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                            QuadTree.Blocks[iCPU][iBLK]->block.infoSW[i_neighbour].level) <= +1) ) {
                         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_NOCHANGE) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoSW[i_neighbour].level) == -1) ) {
		         // Prevent coarsening of solution block.
		         QuadTree.RefineFlags[iCPU][iBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
                         number_of_changes += 1;
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_COARSEN) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoSW[i_neighbour].level) <= +1) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoSW[i_neighbour].level) >= -1)) {
		         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_REFINE) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoSW[i_neighbour].level) == +1)) {
		         // Okay, do nothing!
		      } else if ((QuadTree.RefineFlags[CPU_neighbour][BLK_neighbour] == 
                                  ADAPTIVEBLOCK2D_REFINE) &&
                                 ((QuadTree.Blocks[iCPU][iBLK]->block.info.level -  
                                   QuadTree.Blocks[iCPU][iBLK]->block.infoSW[i_neighbour].level) <= 0)) {
		         // Prevent coarsening of solution block.
		         QuadTree.RefineFlags[iCPU][iBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
                         number_of_changes += 1;
                      } else {
		         // Prevent coarsening of solution block.
		         QuadTree.RefineFlags[iCPU][iBLK] == ADAPTIVEBLOCK2D_NOCHANGE;
                         number_of_changes += 1;
                      } /* endif */
                   } /* endfor */
                } /* endif */

	        // Check siblings to see if block coarsening is possible.
                if (QuadTree.RefineFlags[iCPU][iBLK] == ADAPTIVEBLOCK2D_COARSEN) {
                   number_of_coarsened_sibling = 0;
                   for (i_sibling = 0 ; i_sibling <= 3 ; ++i_sibling) {
                      if (QuadTree.Blocks[iCPU][iBLK]->parent_ptr->child_ptr(i_sibling)->block.used) {
                         CPU_sibling = QuadTree.Blocks[iCPU][iBLK]->parent_ptr->child_ptr(i_sibling)->block.info.cpu;
                         BLK_sibling = QuadTree.Blocks[iCPU][iBLK]->parent_ptr->child_ptr(i_sibling)->block.info.blknum;
                         if (QuadTree.RefineFlags[CPU_sibling][BLK_sibling] == ADAPTIVEBLOCK2D_COARSEN) number_of_coarsened_sibling += 1;
                      } /* endif */
		   } /* endfor */
                   // Force coarsening of sibling.
                   if (number_of_coarsened_sibling == 4) {
                      for (i_sibling = 0 ; i_sibling <= 3 ; ++i_sibling) {
                         CPU_sibling = QuadTree.Blocks[iCPU][iBLK]->parent_ptr->child_ptr(i_sibling)->block.info.cpu;
                         BLK_sibling = QuadTree.Blocks[iCPU][iBLK]->parent_ptr->child_ptr(i_sibling)->block.info.blknum;
                         if (QuadTree.RefineFlags[CPU_sibling][BLK_sibling] != ADAPTIVEBLOCK2D_COARSEN) {
                            QuadTree.RefineFlags[CPU_sibling][BLK_sibling] = ADAPTIVEBLOCK2D_COARSEN;
                            number_of_changes += 1;
                         } /* endif */
                      } /* endfor */
                   // Prevent coarsening of sibling.
                   } else {
                      for (i_sibling = 0 ; i_sibling <= 3 ; ++i_sibling) {
                         if (QuadTree.Blocks[iCPU][iBLK]->parent_ptr->child_ptr(i_sibling)->block.used) {
                            CPU_sibling = QuadTree.Blocks[iCPU][iBLK]->parent_ptr->child_ptr(i_sibling)->block.info.cpu;
                            BLK_sibling = QuadTree.Blocks[iCPU][iBLK]->parent_ptr->child_ptr(i_sibling)->block.info.blknum;
                            if (QuadTree.RefineFlags[CPU_sibling][BLK_sibling] == ADAPTIVEBLOCK2D_COARSEN) {
                               QuadTree.RefineFlags[CPU_sibling][BLK_sibling] = ADAPTIVEBLOCK2D_NOCHANGE;
                               number_of_changes += 1;
                            } /* endif */
                         } /* endif */
                      } /* endfor */
                   } /* endif */
                } /* endif */

             } /* endif */

          } /* endif */
       } /* endfor */
    } /* endfor */


    /* Reassign the local refinement flags for each of the
       local solution blocks so that they are consistent 
       with the the global quadtree list. */

    iCPU = LocalSolnBlockList.ThisCPU;
    for ( iBLK = 0 ; iBLK <= QuadTree.Nblk-1 ; ++iBLK ) {
	LocalSolnBlockList.RefineFlag[iBLK] = QuadTree.RefineFlags[iCPU][iBLK];
    } /* endfor */

//     if (number_of_changes == 0 &&
//         QuadTree.numberToBeRefined() == 0 &&
//         QuadTree.numberToBeCoarsened() == 0 &&
//         QuadTree.highestRefinementLevel() < QuadTree.MaximumRefinementLevel) {
//        QuadTree.refineAll();
//        number_of_changes = QuadTree.countUsedBlocks();
//     } /* endif */

//     cout << "Number of Passes: " << number_of_passes << " " << number_of_changes << "\n"; cout.flush();
//     for ( iCPU = 0 ; iCPU <= QuadTree.Ncpu-1 ; ++iCPU ) {
//        for ( iBLK = 0 ; iBLK <= QuadTree.Nblk-1 ; ++iBLK ) {
//           cout << "Refinement Flag: " << iCPU << " " << iBLK << " " << QuadTree.RefineFlags[iCPU][iBLK] << "\n"; cout.flush();
//         } /* endfor */
//     } /* endfor */

    if (number_of_changes == 0) break;
    } /* endfor */

}

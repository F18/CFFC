/* AMR.h:  Header file for templated subroutines for carrying 
           out the adaptive mesh refinement (AMR). */

#ifndef _AMR_INCLUDED
#define _AMR_INCLUDED

/* Include 3D  multiblock grid, adaptive block,
   and octree header files. */

#ifndef _ADAPTIVEBLOCK3D_INCLUDED
#include "../AMR/AdaptiveBlock3D.h"
#endif // _ADAPTIVEBLOCK3D_INCLUDED

#ifndef _OCTREE_INCLUDED
#include "../AMR/Octree.h"
#endif // _OCTREE_INCLUDED

#ifndef _GRID3D_HEXA_MULTIBLK_INCLUDED
#include "../Grid/Grid3DHexaMultiBlock.h"
#endif // _GRID3D_HEXA_MULTIBLK_INCLUDED

#ifndef _MULTIBLOCK_LIST_INCLUDED
#include "../HexaBlock/HexaMultiBlock.h"
#endif //_MULTIBLOCK_LIST_INCLUDED

/* Include Morton re-ordering header file. */

#ifndef _MORTON_ORDERING_INCLUDED
#include "MortonOrdering.h"
#endif // _MORTON_ORDERING_INCLUDED

/******************************************************************
 * AMR -- Templated subroutines.                                  *
 ******************************************************************/

/******************************************************************
 * Routine: Create_Initial_Solution_Blocks                        *
 *                                                                *
 * Assigns and creates (allocates) initial 3D ocrilateral         *
 * solution blocks corresponding to the initial grid.  This       *
 * routine also creates the octree, local block, and global       *
 * block resource list data structures for performing the         *
 * block-base adaptive mesh refinement (AMR) and determines the   *
 * roots of the octree data structure.                            *
 *                                                                *
 ******************************************************************/
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
int Create_Initial_Solution_Blocks(Grid3D_Hexa_Multi_Block_List                            &Initial_Mesh,
                                   Hexa_Multi_Block<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> > &Local_Solution_Blocks,
                                   Input_Parameters<SOLN_pSTATE, SOLN_cSTATE>              &Input,
                                   Octree_DataStructure                                    &Octree,
                                   AdaptiveBlock3D_ResourceList                            &Global_Adaptive_Block_List,
                                   AdaptiveBlock3D_List                                    &Local_Adaptive_Block_List) {
   
  int n_cpu, n_blk, n_grid_blks, neighbour_cpu, neighbour_blk, number_of_solution_variables; 
  int *gridblk_cpu_number, *gridblk_blk_number;
   
   /* Create (allocate) octree data structure. */

   Octree_DataStructure::Create_Octree_Data_Structure(Octree,
                                                      Input.Grid_IP.NBlk_Idir,
                                                      Input.Grid_IP.NBlk_Jdir,
                                                      Input.Grid_IP.NBlk_Kdir,
                                                      Input.AMR_IP.Number_of_Processors,
                                                      Input.AMR_IP.Number_of_Blocks_Per_Processor);
   
   Octree.MaximumRefinementLevel = Input.AMR_IP.Maximum_Refinement_Level-1;
   Octree.MinimumRefinementLevel = Input.AMR_IP.Minimum_Refinement_Level-1;
   
   /* Set the thresholds for refinement and coarsening of the mesh. */
   
   Octree.RefineThreshold = Input.AMR_IP.Threshold_for_Refinement;
   Octree.CoarsenThreshold = Input.AMR_IP.Threshold_for_Coarsening;
   
   /* Create (allocate) array of local 3D hexahedral solution blocks. */
   
   Local_Solution_Blocks.Allocate(Input.AMR_IP.Number_of_Blocks_Per_Processor);

   AdaptiveBlock3D_ResourceList::Create_Block_Resource_List(Global_Adaptive_Block_List,
                                                            Input.AMR_IP.Number_of_Processors,
                                                            Input.AMR_IP.Number_of_Blocks_Per_Processor);
   Local_Adaptive_Block_List.allocate(Input.AMR_IP.Number_of_Blocks_Per_Processor);
   Local_Adaptive_Block_List.ThisCPU = Global_Adaptive_Block_List.ThisCPU;
   
   /* Loop over all initial mesh blocks and determine the CPU and local block numbers for
      each grid block. */

   gridblk_cpu_number = new int[Initial_Mesh.NBlk];
   gridblk_blk_number = new int[Initial_Mesh.NBlk];

   n_grid_blks = 0;

   for (int nb = 0; nb <= Initial_Mesh.NBlk-1; ++nb) {
      if (Initial_Mesh.Grid_Blks[nb].Allocated) { // Mesh block is used!!!!
         // Get next free solution block from list of available
         // solution blocks.
         if (Global_Adaptive_Block_List.Nfree > 0) {
	    n_grid_blks += 1;
            n_cpu = Global_Adaptive_Block_List.nextCPU();
            n_blk = Global_Adaptive_Block_List.nextBlock();
            Global_Adaptive_Block_List.update_next();
            gridblk_cpu_number[nb] = n_cpu;
            gridblk_blk_number[nb] = n_blk;
         } else {
            cout << "\n"
                 << " AMR Error: Create_Initial_Solution_Blocks, "
                 << "Insufficient number of hexahedral solution blocks.\n";
            return (1);
         } /* endif */
      } /* endif */
   } /* endfor */
   
   /* Loop over all initial mesh or grid blocks and assign octree root blocks. */

   n_grid_blks = 0;

   for (int nb = 0; nb <= Initial_Mesh.NBlk-1; ++nb) {
      if (Initial_Mesh.Grid_Blks[nb].Allocated) { // Mesh block is used!!!!
	 // Obtain CPU and block number.
	 n_grid_blks += 1;
         n_cpu = gridblk_cpu_number[nb];
         n_blk = gridblk_blk_number[nb];

         // Assign block information to appropriate octree root solution block.
         Octree.Roots[nb].block.used = ADAPTIVEBLOCK3D_USED;
         Octree.Roots[nb].block.info.gblknum = n_grid_blks-1;
         Octree.Roots[nb].block.info.cpu = n_cpu;
         Octree.Roots[nb].block.info.blknum = n_blk;
         Octree.Roots[nb].block.info.dimen.i = Initial_Mesh.Grid_Blks[nb].NCi -
                                               2*Initial_Mesh.Grid_Blks[nb].Nghost;
         Octree.Roots[nb].block.info.dimen.j = Initial_Mesh.Grid_Blks[nb].NCj -
                                               2*Initial_Mesh.Grid_Blks[nb].Nghost;
         Octree.Roots[nb].block.info.dimen.k = Initial_Mesh.Grid_Blks[nb].NCk -
                                               2*Initial_Mesh.Grid_Blks[nb].Nghost;
         Octree.Roots[nb].block.info.dimen.ghost = 2;
         Octree.Roots[nb].block.info.sector = ADAPTIVEBLOCK3D_SECTOR_NONE;
         Octree.Roots[nb].block.info.level = 0;
         Octree.Roots[nb].parent_ptr = NULL;
         Octree.Roots[nb].childTNW_ptr = NULL;
         Octree.Roots[nb].childTNE_ptr = NULL;
         Octree.Roots[nb].childTSE_ptr = NULL;
         Octree.Roots[nb].childTSW_ptr = NULL;
         Octree.Roots[nb].childBNW_ptr = NULL;
         Octree.Roots[nb].childBNE_ptr = NULL;
         Octree.Roots[nb].childBSE_ptr = NULL;
         Octree.Roots[nb].childBSW_ptr = NULL;

         // Assign the number of neighbouring grid blocks in each direction
         // for each of the block in the octree roots.
         Octree.Roots[nb].block.nT = Initial_Mesh.Connectivity[nb].num_neighT;
         Octree.Roots[nb].block.nB = Initial_Mesh.Connectivity[nb].num_neighB;
         Octree.Roots[nb].block.nN = Initial_Mesh.Connectivity[nb].num_neighN;
         Octree.Roots[nb].block.nS = Initial_Mesh.Connectivity[nb].num_neighS;
         Octree.Roots[nb].block.nE = Initial_Mesh.Connectivity[nb].num_neighE;
         Octree.Roots[nb].block.nW = Initial_Mesh.Connectivity[nb].num_neighW;
         Octree.Roots[nb].block.nNW = Initial_Mesh.Connectivity[nb].num_neighNW;
         Octree.Roots[nb].block.nNE = Initial_Mesh.Connectivity[nb].num_neighNE;
         Octree.Roots[nb].block.nSE = Initial_Mesh.Connectivity[nb].num_neighSE;
         Octree.Roots[nb].block.nSW = Initial_Mesh.Connectivity[nb].num_neighSW;
         Octree.Roots[nb].block.nTN = Initial_Mesh.Connectivity[nb].num_neighTN;
         Octree.Roots[nb].block.nTS = Initial_Mesh.Connectivity[nb].num_neighTS;
         Octree.Roots[nb].block.nTE = Initial_Mesh.Connectivity[nb].num_neighTE;
         Octree.Roots[nb].block.nTW = Initial_Mesh.Connectivity[nb].num_neighTW;
         Octree.Roots[nb].block.nBN = Initial_Mesh.Connectivity[nb].num_neighBN;
         Octree.Roots[nb].block.nBS = Initial_Mesh.Connectivity[nb].num_neighBS;
         Octree.Roots[nb].block.nBE = Initial_Mesh.Connectivity[nb].num_neighBE;
         Octree.Roots[nb].block.nBW = Initial_Mesh.Connectivity[nb].num_neighBW;
         Octree.Roots[nb].block.nTNW  = Initial_Mesh.Connectivity[nb].num_neighTNW;
         Octree.Roots[nb].block.nTSW = Initial_Mesh.Connectivity[nb].num_neighTSW;
         Octree.Roots[nb].block.nTNE = Initial_Mesh.Connectivity[nb].num_neighTNE;
         Octree.Roots[nb].block.nTSE = Initial_Mesh.Connectivity[nb].num_neighTSE;
         Octree.Roots[nb].block.nBSE  = Initial_Mesh.Connectivity[nb].num_neighBSE;
         Octree.Roots[nb].block.nBSW = Initial_Mesh.Connectivity[nb].num_neighBSW;
         Octree.Roots[nb].block.nBNE = Initial_Mesh.Connectivity[nb].num_neighBNE;
         Octree.Roots[nb].block.nBNW = Initial_Mesh.Connectivity[nb].num_neighBNW;

         // Assign the block boundary element "on grid boundaries" info for each root block.
         Octree.Roots[nb].block.info.be = Initial_Mesh.Connectivity[nb].be;
         
         // Assign block number and CPU number of the neighbouring blocks in each direction
         // for each of the block in the octree roots.  Also assign block orientation information 
         // for each of the neighbouring blocks.
         if (Octree.Roots[nb].block.nT > 0) {
            Octree.Roots[nb].block.infoT[0].blknum = gridblk_blk_number[Initial_Mesh.Connectivity[nb].neighT];
            Octree.Roots[nb].block.infoT[0].cpu = gridblk_cpu_number[Initial_Mesh.Connectivity[nb].neighT];
            Octree.Roots[nb].block.infoT[0].blkorient = Initial_Mesh.Connectivity[nb].neighT_info;
         } /* endif */

         if (Octree.Roots[nb].block.nB > 0) {
            Octree.Roots[nb].block.infoB[0].blknum = gridblk_blk_number[Initial_Mesh.Connectivity[nb].neighB];
            Octree.Roots[nb].block.infoB[0].cpu = gridblk_cpu_number[Initial_Mesh.Connectivity[nb].neighB];
            Octree.Roots[nb].block.infoB[0].blkorient = Initial_Mesh.Connectivity[nb].neighB_info;
         } /* endif */
         
         if (Octree.Roots[nb].block.nN > 0) {
            Octree.Roots[nb].block.infoN[0].blknum = gridblk_blk_number[Initial_Mesh.Connectivity[nb].neighN];
            Octree.Roots[nb].block.infoN[0].cpu = gridblk_cpu_number[Initial_Mesh.Connectivity[nb].neighN];
            Octree.Roots[nb].block.infoN[0].blkorient = Initial_Mesh.Connectivity[nb].neighN_info;
         } /* endif */
         
         if (Octree.Roots[nb].block.nS > 0) {
            Octree.Roots[nb].block.infoS[0].blknum = gridblk_blk_number[Initial_Mesh.Connectivity[nb].neighS];
            Octree.Roots[nb].block.infoS[0].cpu = gridblk_cpu_number[Initial_Mesh.Connectivity[nb].neighS];
            Octree.Roots[nb].block.infoS[0].blkorient = Initial_Mesh.Connectivity[nb].neighS_info;
         } /* endif */
 
         if (Octree.Roots[nb].block.nE > 0) {
            Octree.Roots[nb].block.infoE[0].blknum = gridblk_blk_number[Initial_Mesh.Connectivity[nb].neighE];
            Octree.Roots[nb].block.infoE[0].cpu = gridblk_cpu_number[Initial_Mesh.Connectivity[nb].neighE];
            Octree.Roots[nb].block.infoE[0].blkorient = Initial_Mesh.Connectivity[nb].neighE_info;
         } /* endif */

         if (Octree.Roots[nb].block.nW > 0) {
            Octree.Roots[nb].block.infoW[0].blknum = gridblk_blk_number[Initial_Mesh.Connectivity[nb].neighW];
            Octree.Roots[nb].block.infoW[0].cpu = gridblk_cpu_number[Initial_Mesh.Connectivity[nb].neighW];
            Octree.Roots[nb].block.infoW[0].blkorient = Initial_Mesh.Connectivity[nb].neighW_info;
         } /* endif */
            
         if (Octree.Roots[nb].block.nNW > 0) {
            Octree.Roots[nb].block.infoNW[0].blknum = gridblk_blk_number[Initial_Mesh.Connectivity[nb].neighNW[0]];
            Octree.Roots[nb].block.infoNW[0].cpu = gridblk_cpu_number[Initial_Mesh.Connectivity[nb].neighNW[0]];
            Octree.Roots[nb].block.infoNW[0].blkorient = Initial_Mesh.Connectivity[nb].neighNW_info[0];
         } /* endif */

         if (Octree.Roots[nb].block.nNE > 0) {
            Octree.Roots[nb].block.infoNE[0].blknum = gridblk_blk_number[Initial_Mesh.Connectivity[nb].neighNE[0]];
            Octree.Roots[nb].block.infoNE[0].cpu = gridblk_cpu_number[Initial_Mesh.Connectivity[nb].neighNE[0]];
            Octree.Roots[nb].block.infoNE[0].blkorient = Initial_Mesh.Connectivity[nb].neighNE_info[0];
         } /* endif */
         
         if (Octree.Roots[nb].block.nSE > 0) {
            Octree.Roots[nb].block.infoSE[0].blknum = gridblk_blk_number[Initial_Mesh.Connectivity[nb].neighSE[0]];
            Octree.Roots[nb].block.infoSE[0].cpu = gridblk_cpu_number[Initial_Mesh.Connectivity[nb].neighSE[0]];
            Octree.Roots[nb].block.infoSE[0].blkorient = Initial_Mesh.Connectivity[nb].neighSE_info[0];
         } /* endif */

         if (Octree.Roots[nb].block.nSW > 0) {
            Octree.Roots[nb].block.infoSW[0].blknum = gridblk_blk_number[Initial_Mesh.Connectivity[nb].neighSW[0]];
            Octree.Roots[nb].block.infoSW[0].cpu = gridblk_cpu_number[Initial_Mesh.Connectivity[nb].neighSW[0]];
            Octree.Roots[nb].block.infoSW[0].blkorient = Initial_Mesh.Connectivity[nb].neighSW_info[0];
         } /* endif */

         if (Octree.Roots[nb].block.nTN > 0) {
            Octree.Roots[nb].block.infoTN[0].blknum = gridblk_blk_number[Initial_Mesh.Connectivity[nb].neighTN[0]];
            Octree.Roots[nb].block.infoTN[0].cpu = gridblk_cpu_number[Initial_Mesh.Connectivity[nb].neighTN[0]];
            Octree.Roots[nb].block.infoTN[0].blkorient = Initial_Mesh.Connectivity[nb].neighTN_info[0];
         } /* endif */

         if (Octree.Roots[nb].block.nTS > 0) {
            Octree.Roots[nb].block.infoTS[0].blknum = gridblk_blk_number[Initial_Mesh.Connectivity[nb].neighTS[0]];
            Octree.Roots[nb].block.infoTS[0].cpu = gridblk_cpu_number[Initial_Mesh.Connectivity[nb].neighTS[0]];
            Octree.Roots[nb].block.infoTS[0].blkorient = Initial_Mesh.Connectivity[nb].neighTS_info[0];
         } /* endif */
         
         if (Octree.Roots[nb].block.nTE > 0) {
	    Octree.Roots[nb].block.infoTE[0].blknum = gridblk_blk_number[Initial_Mesh.Connectivity[nb].neighTE[0]];
            Octree.Roots[nb].block.infoTE[0].cpu = gridblk_cpu_number[Initial_Mesh.Connectivity[nb].neighTE[0]];
            Octree.Roots[nb].block.infoTE[0].blkorient = Initial_Mesh.Connectivity[nb].neighTE_info[0];
         } /* endif */

         if (Octree.Roots[nb].block.nTW > 0) {
            Octree.Roots[nb].block.infoTW[0].blknum = gridblk_blk_number[Initial_Mesh.Connectivity[nb].neighTW[0]];
            Octree.Roots[nb].block.infoTW[0].cpu = gridblk_cpu_number[Initial_Mesh.Connectivity[nb].neighTW[0]];
            Octree.Roots[nb].block.infoTW[0].blkorient = Initial_Mesh.Connectivity[nb].neighTW_info[0];
         } /* endif */

         if (Octree.Roots[nb].block.nBN > 0) {
            Octree.Roots[nb].block.infoBN[0].blknum = gridblk_blk_number[Initial_Mesh.Connectivity[nb].neighBN[0]];
            Octree.Roots[nb].block.infoBN[0].cpu = gridblk_cpu_number[Initial_Mesh.Connectivity[nb].neighBN[0]];
            Octree.Roots[nb].block.infoBN[0].blkorient = Initial_Mesh.Connectivity[nb].neighBN_info[0];
         } /* endif */

         if (Octree.Roots[nb].block.nBS > 0) {
            Octree.Roots[nb].block.infoBS[0].blknum = gridblk_blk_number[Initial_Mesh.Connectivity[nb].neighBS[0]];
            Octree.Roots[nb].block.infoBS[0].cpu = gridblk_cpu_number[Initial_Mesh.Connectivity[nb].neighBS[0]];
            Octree.Roots[nb].block.infoBS[0].blkorient = Initial_Mesh.Connectivity[nb].neighBS_info[0];
         } /* endif */

         if (Octree.Roots[nb].block.nBE > 0) {
            Octree.Roots[nb].block.infoBE[0].blknum = gridblk_blk_number[Initial_Mesh.Connectivity[nb].neighBE[0]];
            Octree.Roots[nb].block.infoBE[0].cpu = gridblk_cpu_number[Initial_Mesh.Connectivity[nb].neighBE[0]];
            Octree.Roots[nb].block.infoBE[0].blkorient = Initial_Mesh.Connectivity[nb].neighBE_info[0];
         } /* endif */

         if (Octree.Roots[nb].block.nBW > 0) {
            Octree.Roots[nb].block.infoBW[0].blknum = gridblk_blk_number[Initial_Mesh.Connectivity[nb].neighBW[0]];
            Octree.Roots[nb].block.infoBW[0].cpu = gridblk_cpu_number[Initial_Mesh.Connectivity[nb].neighBW[0]];
            Octree.Roots[nb].block.infoBW[0].blkorient = Initial_Mesh.Connectivity[nb].neighBW_info[0];
         } /* endif */

         if (Octree.Roots[nb].block.nTNW > 0) {
            Octree.Roots[nb].block.infoTNW[0].blknum = gridblk_blk_number[Initial_Mesh.Connectivity[nb].neighTNW[0]];
            Octree.Roots[nb].block.infoTNW[0].cpu = gridblk_cpu_number[Initial_Mesh.Connectivity[nb].neighTNW[0]];
            Octree.Roots[nb].block.infoTNW[0].blkorient = Initial_Mesh.Connectivity[nb].neighTNW_info[0];
         } /* endif */

         if (Octree.Roots[nb].block.nTSW > 0) {
            Octree.Roots[nb].block.infoTSW[0].blknum = gridblk_blk_number[Initial_Mesh.Connectivity[nb].neighTSW[0]];
            Octree.Roots[nb].block.infoTSW[0].cpu = gridblk_cpu_number[Initial_Mesh.Connectivity[nb].neighTSW[0]];
            Octree.Roots[nb].block.infoTSW[0].blkorient = Initial_Mesh.Connectivity[nb].neighTSW_info[0];
         } /* endif */

         if (Octree.Roots[nb].block.nTNE > 0) {
            Octree.Roots[nb].block.infoTNE[0].blknum = gridblk_blk_number[Initial_Mesh.Connectivity[nb].neighTNE[0]];
            Octree.Roots[nb].block.infoTNE[0].cpu = gridblk_cpu_number[Initial_Mesh.Connectivity[nb].neighTNE[0]];
            Octree.Roots[nb].block.infoTNE[0].blkorient = Initial_Mesh.Connectivity[nb].neighTNE_info[0];
         } /* endif */

         if (Octree.Roots[nb].block.nTSE > 0) {
            Octree.Roots[nb].block.infoTSE[0].blknum = gridblk_blk_number[Initial_Mesh.Connectivity[nb].neighTSE[0]];
            Octree.Roots[nb].block.infoTSE[0].cpu = gridblk_cpu_number[Initial_Mesh.Connectivity[nb].neighTSE[0]];
            Octree.Roots[nb].block.infoTSE[0].blkorient = Initial_Mesh.Connectivity[nb].neighTSE_info[0];
         } /* endif */

         if (Octree.Roots[nb].block.nBNW > 0) {
            Octree.Roots[nb].block.infoBNW[0].blknum = gridblk_blk_number[Initial_Mesh.Connectivity[nb].neighBNW[0]];
            Octree.Roots[nb].block.infoBNW[0].cpu = gridblk_cpu_number[Initial_Mesh.Connectivity[nb].neighBNW[0]];
            Octree.Roots[nb].block.infoBNW[0].blkorient = Initial_Mesh.Connectivity[nb].neighBNW_info[0];
         } /* endif */

         if (Octree.Roots[nb].block.nBSW > 0) {
            Octree.Roots[nb].block.infoBSW[0].blknum = gridblk_blk_number[Initial_Mesh.Connectivity[nb].neighBSW[0]];
            Octree.Roots[nb].block.infoBSW[0].cpu = gridblk_cpu_number[Initial_Mesh.Connectivity[nb].neighBSW[0]];
            Octree.Roots[nb].block.infoBSW[0].blkorient = Initial_Mesh.Connectivity[nb].neighBSW_info[0];
         } /* endif */

         if (Octree.Roots[nb].block.nBNE > 0) {
            Octree.Roots[nb].block.infoBNE[0].blknum = gridblk_blk_number[Initial_Mesh.Connectivity[nb].neighBNE[0]];
            Octree.Roots[nb].block.infoBNE[0].cpu = gridblk_cpu_number[Initial_Mesh.Connectivity[nb].neighBNE[0]];
            Octree.Roots[nb].block.infoBNE[0].blkorient = Initial_Mesh.Connectivity[nb].neighBNE_info[0];
         } /* endif */

         if (Octree.Roots[nb].block.nBSE > 0) {
            Octree.Roots[nb].block.infoBSE[0].blknum = gridblk_blk_number[Initial_Mesh.Connectivity[nb].neighBSE[0]];
            Octree.Roots[nb].block.infoBSE[0].cpu = gridblk_cpu_number[Initial_Mesh.Connectivity[nb].neighBSE[0]];
            Octree.Roots[nb].block.infoBSE[0].blkorient = Initial_Mesh.Connectivity[nb].neighBSE_info[0];
         } /* endif */

         // Finally, assign appropriate octree block pointer to octree root solution block.                    
         Octree.Blocks[n_cpu][n_blk] = &(Octree.Roots[nb]);
      } else{
         Octree.Roots[nb].block.used = 0;
      } /* endif */
      
   } /* endfor */

   delete []gridblk_cpu_number;
   delete []gridblk_blk_number;

   /* Now, assign missing neighbour block information for octree root blocks: 
      global block number, dimensions, sector, level assignment, "on grid boundary" info. */

   for (int nr = 0; nr <= Octree.NR-1; ++nr) {
      if (Octree.Roots[nr].block.used) { // Adaptive root block is used!!!!
         if (Octree.Roots[nr].block.nT > 0) {
	    neighbour_cpu = Octree.Roots[nr].block.infoT[0].cpu;
	    neighbour_blk = Octree.Roots[nr].block.infoT[0].blknum;
            Octree.Roots[nr].block.infoT[0].gblknum = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.gblknum;
            Octree.Roots[nr].block.infoT[0].dimen = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.dimen;
            Octree.Roots[nr].block.infoT[0].sector = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.sector;
            Octree.Roots[nr].block.infoT[0].level = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.level;
            Octree.Roots[nr].block.infoT[0].be = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.be;
         } /* endif */
      
         if (Octree.Roots[nr].block.nB > 0) {
	    neighbour_cpu = Octree.Roots[nr].block.infoB[0].cpu;
	    neighbour_blk = Octree.Roots[nr].block.infoB[0].blknum;
            Octree.Roots[nr].block.infoB[0].gblknum = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.gblknum;
            Octree.Roots[nr].block.infoB[0].dimen = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.dimen;
            Octree.Roots[nr].block.infoB[0].sector = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.sector;
            Octree.Roots[nr].block.infoB[0].level = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.level;
            Octree.Roots[nr].block.infoB[0].be = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.be;
         } /* endif */
      
         if (Octree.Roots[nr].block.nN > 0) {
	    neighbour_cpu = Octree.Roots[nr].block.infoN[0].cpu;
	    neighbour_blk = Octree.Roots[nr].block.infoN[0].blknum;
            Octree.Roots[nr].block.infoN[0].gblknum = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.gblknum;
            Octree.Roots[nr].block.infoN[0].dimen = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.dimen;
            Octree.Roots[nr].block.infoN[0].sector = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.sector;
            Octree.Roots[nr].block.infoN[0].level = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.level;
            Octree.Roots[nr].block.infoN[0].be = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.be;
         } /* endif */
      
         if (Octree.Roots[nr].block.nS > 0) {
	    neighbour_cpu = Octree.Roots[nr].block.infoS[0].cpu;
	    neighbour_blk = Octree.Roots[nr].block.infoS[0].blknum;
            Octree.Roots[nr].block.infoS[0].gblknum = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.gblknum;
            Octree.Roots[nr].block.infoS[0].dimen = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.dimen;
            Octree.Roots[nr].block.infoS[0].sector = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.sector;
            Octree.Roots[nr].block.infoS[0].level = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.level;
            Octree.Roots[nr].block.infoS[0].be = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.be;
         } /* endif */
      
         if (Octree.Roots[nr].block.nE > 0) {
	    neighbour_cpu = Octree.Roots[nr].block.infoE[0].cpu;
	    neighbour_blk = Octree.Roots[nr].block.infoE[0].blknum;
            Octree.Roots[nr].block.infoE[0].gblknum = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.gblknum;
            Octree.Roots[nr].block.infoE[0].dimen = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.dimen;
            Octree.Roots[nr].block.infoE[0].sector = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.sector;
            Octree.Roots[nr].block.infoE[0].level = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.level;
            Octree.Roots[nr].block.infoE[0].be = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.be;
         } /* endif */
      
         if (Octree.Roots[nr].block.nW > 0) {
	    neighbour_cpu = Octree.Roots[nr].block.infoW[0].cpu;
	    neighbour_blk = Octree.Roots[nr].block.infoW[0].blknum;
            Octree.Roots[nr].block.infoW[0].gblknum = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.gblknum;
            Octree.Roots[nr].block.infoW[0].dimen = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.dimen;
            Octree.Roots[nr].block.infoW[0].sector = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.sector;
            Octree.Roots[nr].block.infoW[0].level = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.level;
            Octree.Roots[nr].block.infoW[0].be = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.be;
         } /* endif */
      
         if (Octree.Roots[nr].block.nNW > 0) {
	    neighbour_cpu = Octree.Roots[nr].block.infoNW[0].cpu;
	    neighbour_blk = Octree.Roots[nr].block.infoNW[0].blknum;
            Octree.Roots[nr].block.infoNW[0].gblknum = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.gblknum;
            Octree.Roots[nr].block.infoNW[0].dimen = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.dimen;
            Octree.Roots[nr].block.infoNW[0].sector = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.sector;
            Octree.Roots[nr].block.infoNW[0].level = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.level;
            Octree.Roots[nr].block.infoNW[0].be = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.be;
         } /* endif */
       
         if (Octree.Roots[nr].block.nNE > 0) {
	    neighbour_cpu = Octree.Roots[nr].block.infoNE[0].cpu;
	    neighbour_blk = Octree.Roots[nr].block.infoNE[0].blknum;
            Octree.Roots[nr].block.infoNE[0].gblknum = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.gblknum;
            Octree.Roots[nr].block.infoNE[0].dimen = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.dimen;
            Octree.Roots[nr].block.infoNE[0].sector = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.sector;
            Octree.Roots[nr].block.infoNE[0].level = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.level;
            Octree.Roots[nr].block.infoNE[0].be = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.be;
         } /* endif */
      
         if (Octree.Roots[nr].block.nSE > 0) {
	    neighbour_cpu = Octree.Roots[nr].block.infoSE[0].cpu;
	    neighbour_blk = Octree.Roots[nr].block.infoSE[0].blknum;
            Octree.Roots[nr].block.infoSE[0].gblknum = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.gblknum;
            Octree.Roots[nr].block.infoSE[0].dimen = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.dimen;
            Octree.Roots[nr].block.infoSE[0].sector = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.sector;
            Octree.Roots[nr].block.infoSE[0].level = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.level;
            Octree.Roots[nr].block.infoSE[0].be = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.be;
         } /* endif */
      
         if (Octree.Roots[nr].block.nSW > 0) {
	    neighbour_cpu = Octree.Roots[nr].block.infoSW[0].cpu;
	    neighbour_blk = Octree.Roots[nr].block.infoSW[0].blknum;
            Octree.Roots[nr].block.infoSW[0].gblknum = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.gblknum;
            Octree.Roots[nr].block.infoSW[0].dimen = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.dimen;
            Octree.Roots[nr].block.infoSW[0].sector = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.sector;
            Octree.Roots[nr].block.infoSW[0].level = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.level;
            Octree.Roots[nr].block.infoSW[0].be = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.be;
         } /* endif */
      
         if (Octree.Roots[nr].block.nTN > 0) {
	    neighbour_cpu = Octree.Roots[nr].block.infoTN[0].cpu;
	    neighbour_blk = Octree.Roots[nr].block.infoTN[0].blknum;
            Octree.Roots[nr].block.infoTN[0].gblknum = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.gblknum;
            Octree.Roots[nr].block.infoTN[0].dimen = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.dimen;
            Octree.Roots[nr].block.infoTN[0].sector = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.sector;
            Octree.Roots[nr].block.infoTN[0].level = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.level;
            Octree.Roots[nr].block.infoTN[0].be = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.be;
         } /* endif */
      
         if (Octree.Roots[nr].block.nTS > 0) {
	    neighbour_cpu = Octree.Roots[nr].block.infoTS[0].cpu;
	    neighbour_blk = Octree.Roots[nr].block.infoTS[0].blknum;
            Octree.Roots[nr].block.infoTS[0].gblknum = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.gblknum;
            Octree.Roots[nr].block.infoTS[0].dimen = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.dimen;
            Octree.Roots[nr].block.infoTS[0].sector = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.sector;
            Octree.Roots[nr].block.infoTS[0].level = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.level;
            Octree.Roots[nr].block.infoTS[0].be = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.be;
         } /* endif */
      
         if (Octree.Roots[nr].block.nTE > 0) {
	    neighbour_cpu = Octree.Roots[nr].block.infoTE[0].cpu;
	    neighbour_blk = Octree.Roots[nr].block.infoTE[0].blknum;
            Octree.Roots[nr].block.infoTE[0].gblknum = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.gblknum;
            Octree.Roots[nr].block.infoTE[0].dimen = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.dimen;
            Octree.Roots[nr].block.infoTE[0].sector = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.sector;
            Octree.Roots[nr].block.infoTE[0].level = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.level;
            Octree.Roots[nr].block.infoTE[0].be = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.be;
         } /* endif */
      
         if (Octree.Roots[nr].block.nTW > 0) {
	    neighbour_cpu = Octree.Roots[nr].block.infoTW[0].cpu;
	    neighbour_blk = Octree.Roots[nr].block.infoTW[0].blknum;
            Octree.Roots[nr].block.infoTW[0].gblknum = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.gblknum;
            Octree.Roots[nr].block.infoTW[0].dimen = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.dimen;
            Octree.Roots[nr].block.infoTW[0].sector = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.sector;
            Octree.Roots[nr].block.infoTW[0].level = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.level;
            Octree.Roots[nr].block.infoTW[0].be = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.be;
         } /* endif */
      
         if (Octree.Roots[nr].block.nBN > 0) {
	    neighbour_cpu = Octree.Roots[nr].block.infoBN[0].cpu;
	    neighbour_blk = Octree.Roots[nr].block.infoBN[0].blknum;
            Octree.Roots[nr].block.infoBN[0].gblknum = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.gblknum;
            Octree.Roots[nr].block.infoBN[0].dimen = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.dimen;
            Octree.Roots[nr].block.infoBN[0].sector = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.sector;
            Octree.Roots[nr].block.infoBN[0].level = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.level;
            Octree.Roots[nr].block.infoBN[0].be = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.be;
         } /* endif */
      
         if (Octree.Roots[nr].block.nBS > 0) {
	    neighbour_cpu = Octree.Roots[nr].block.infoBS[0].cpu;
	    neighbour_blk = Octree.Roots[nr].block.infoBS[0].blknum;
            Octree.Roots[nr].block.infoBS[0].gblknum = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.gblknum;
            Octree.Roots[nr].block.infoBS[0].dimen = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.dimen;
            Octree.Roots[nr].block.infoBS[0].sector = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.sector;
            Octree.Roots[nr].block.infoBS[0].level = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.level;
            Octree.Roots[nr].block.infoBS[0].be = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.be;
         } /* endif */
      
         if (Octree.Roots[nr].block.nBE > 0) {
	    neighbour_cpu = Octree.Roots[nr].block.infoBE[0].cpu;
	    neighbour_blk = Octree.Roots[nr].block.infoBE[0].blknum;
            Octree.Roots[nr].block.infoBE[0].gblknum = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.gblknum;
            Octree.Roots[nr].block.infoBE[0].dimen = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.dimen;
            Octree.Roots[nr].block.infoBE[0].sector = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.sector;
            Octree.Roots[nr].block.infoBE[0].level = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.level;
            Octree.Roots[nr].block.infoBE[0].be = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.be;
         } /* endif */
      
         if (Octree.Roots[nr].block.nBW > 0) {
	    neighbour_cpu = Octree.Roots[nr].block.infoBW[0].cpu;
	    neighbour_blk = Octree.Roots[nr].block.infoBW[0].blknum;
            Octree.Roots[nr].block.infoBW[0].gblknum = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.gblknum;
            Octree.Roots[nr].block.infoBW[0].dimen = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.dimen;
            Octree.Roots[nr].block.infoBW[0].sector = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.sector;
            Octree.Roots[nr].block.infoBW[0].level = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.level;
            Octree.Roots[nr].block.infoBW[0].be = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.be;
         } /* endif */
      
         if (Octree.Roots[nr].block.nTNW > 0) {
	    neighbour_cpu = Octree.Roots[nr].block.infoTNW[0].cpu;
	    neighbour_blk = Octree.Roots[nr].block.infoTNW[0].blknum;
            Octree.Roots[nr].block.infoTNW[0].gblknum = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.gblknum;
            Octree.Roots[nr].block.infoTNW[0].dimen = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.dimen;   
            Octree.Roots[nr].block.infoTNW[0].sector = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.sector;
            Octree.Roots[nr].block.infoTNW[0].level = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.level;
            Octree.Roots[nr].block.infoTNW[0].be = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.be;
         } /* endif */
      
         if (Octree.Roots[nr].block.nTSW > 0) {
	    neighbour_cpu = Octree.Roots[nr].block.infoTSW[0].cpu;
	    neighbour_blk = Octree.Roots[nr].block.infoTSW[0].blknum;
            Octree.Roots[nr].block.infoTSW[0].gblknum = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.gblknum;
            Octree.Roots[nr].block.infoTSW[0].dimen = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.dimen;
            Octree.Roots[nr].block.infoTSW[0].sector = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.sector;
            Octree.Roots[nr].block.infoTSW[0].level = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.level;
            Octree.Roots[nr].block.infoTSW[0].be = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.be;
         } /* endif */
      
         if (Octree.Roots[nr].block.nTNE > 0) {
	    neighbour_cpu = Octree.Roots[nr].block.infoTNE[0].cpu;
	    neighbour_blk = Octree.Roots[nr].block.infoTNE[0].blknum;
            Octree.Roots[nr].block.infoTNE[0].gblknum = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.gblknum;
            Octree.Roots[nr].block.infoTNE[0].dimen = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.dimen;
            Octree.Roots[nr].block.infoTNE[0].sector = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.sector;
            Octree.Roots[nr].block.infoTNE[0].level = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.level;
            Octree.Roots[nr].block.infoTNE[0].be = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.be;
         } /* endif */
      
         if (Octree.Roots[nr].block.nTSE > 0) {
	    neighbour_cpu = Octree.Roots[nr].block.infoTSE[0].cpu;
	    neighbour_blk = Octree.Roots[nr].block.infoTSE[0].blknum;
            Octree.Roots[nr].block.infoTSE[0].gblknum = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.gblknum;
            Octree.Roots[nr].block.infoTSE[0].dimen = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.dimen;
            Octree.Roots[nr].block.infoTSE[0].sector = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.sector;
            Octree.Roots[nr].block.infoTSE[0].level = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.level;
            Octree.Roots[nr].block.infoTSE[0].be = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.be;
         } /* endif */
      
         if (Octree.Roots[nr].block.nBNW > 0) {
	    neighbour_cpu = Octree.Roots[nr].block.infoBNW[0].cpu;
	    neighbour_blk = Octree.Roots[nr].block.infoBNW[0].blknum;
            Octree.Roots[nr].block.infoBNW[0].gblknum = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.gblknum;
            Octree.Roots[nr].block.infoBNW[0].dimen = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.dimen;
            Octree.Roots[nr].block.infoBNW[0].sector = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.sector;
            Octree.Roots[nr].block.infoBNW[0].level = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.level;  
            Octree.Roots[nr].block.infoBNW[0].be = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.be;
         } /* endif */
      
         if (Octree.Roots[nr].block.nBSW > 0) {
	    neighbour_cpu = Octree.Roots[nr].block.infoBSW[0].cpu;
	    neighbour_blk = Octree.Roots[nr].block.infoBSW[0].blknum;
            Octree.Roots[nr].block.infoBSW[0].gblknum = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.gblknum;
            Octree.Roots[nr].block.infoBSW[0].dimen = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.dimen;
            Octree.Roots[nr].block.infoBSW[0].sector = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.sector;
            Octree.Roots[nr].block.infoBSW[0].level = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.level;
            Octree.Roots[nr].block.infoBSW[0].be = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.be;
         } /* endif */
      
         if (Octree.Roots[nr].block.nBNE > 0) {
	    neighbour_cpu = Octree.Roots[nr].block.infoBNE[0].cpu;
	    neighbour_blk = Octree.Roots[nr].block.infoBNE[0].blknum;
            Octree.Roots[nr].block.infoBNE[0].gblknum = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.gblknum;
            Octree.Roots[nr].block.infoBNE[0].dimen = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.dimen;
            Octree.Roots[nr].block.infoBNE[0].sector = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.sector;
            Octree.Roots[nr].block.infoBNE[0].level = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.level;
            Octree.Roots[nr].block.infoBNE[0].be = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.be;
         } /* endif */
      
         if (Octree.Roots[nr].block.nBSE > 0) {
	    neighbour_cpu = Octree.Roots[nr].block.infoBSE[0].cpu;
	    neighbour_blk = Octree.Roots[nr].block.infoBSE[0].blknum;
            Octree.Roots[nr].block.infoBSE[0].gblknum = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.gblknum;
            Octree.Roots[nr].block.infoBSE[0].dimen = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.dimen;
            Octree.Roots[nr].block.infoBSE[0].sector = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.sector;
            Octree.Roots[nr].block.infoBSE[0].level = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.level;
            Octree.Roots[nr].block.infoBSE[0].be = Octree.Blocks[neighbour_cpu][neighbour_blk]->block.info.be;
         } /* endif */
      } /* endif */
   } /* endfor */
   
   /* Finally loop over all initial mesh blocks and assign local solution block information 
      and create local solution blocks as required. */

   for (int nr = 0; nr <= Octree.NR-1; ++nr) {
      if (Octree.Roots[nr].block.used) { // Adaptive root block is used!!!!
         // For solution blocks on this processor (or processor
         // element), add block to local list, create the solution
         // block, and copy the appropriate block of the
         // initial hexarilateral mesh to the solution block mesh.
         if (Octree.Roots[nr].block.info.cpu == Global_Adaptive_Block_List.ThisCPU) {
	    n_cpu = Octree.Roots[nr].block.info.cpu;
            n_blk = Octree.Roots[nr].block.info.blknum;
            Local_Adaptive_Block_List.Block[n_blk] = Octree.Roots[nr].block;
            Local_Solution_Blocks.Soln_Blks[n_blk].Flow_Type = Input.i_Flow_Type;
            Local_Solution_Blocks.Soln_Blks[n_blk].Create_Block(Initial_Mesh.Grid_Blks[nr]);
            Local_Solution_Blocks.Block_Used[n_blk] = HEXA_BLOCK_USED;
         } /* endif */
     } /* endif */
   } /* endfor */
   
   /* Renumber all solution blocks, assigning a unique global block number. */

/*    Octree_DataStructure::Renumber_Solution_Blocks(Octree, */
/*                                                   Local_Adaptive_Block_List); */

   /* Allocate memory for all message passing buffers used to send
      solution information between neighbouring solution blocks. */

   // Get the number of variables.
   for (int i_blk = 0 ; i_blk <= Local_Adaptive_Block_List.Nblk-1 ; ++i_blk) {
      if (Local_Adaptive_Block_List.Block[i_blk].used) {
         number_of_solution_variables = Local_Solution_Blocks.Soln_Blks[i_blk].NumVar();
         break;
      } /* endif */
   } /* endif */

   AdaptiveBlock3D_List::Allocate_Message_Buffers(Local_Adaptive_Block_List,
                                                  number_of_solution_variables);
 
   /* Solution block allocation and assignment complete.
      Return pointer to local solution blocks. */

   return(0);

}

/********************************************************
 * Routine: Read_Octree                                 *
 *                                                      *
 * Reads the Octree data structure from a file.         *
 *                                                      *
 ********************************************************/
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
int Read_Octree(Octree_DataStructure                       &Octree,
                AdaptiveBlock3D_ResourceList               &Global_Adaptive_Block_List,
                AdaptiveBlock3D_List                       &Local_Adaptive_Block_List,
                Hexa_Multi_Block<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> > &Local_Solution_Blocks,
                Input_Parameters<SOLN_pSTATE, SOLN_cSTATE> &Input) {

    int i, nr, nri, nrj, nrk, ncpu, nblk, number_of_solution_variables;
    int iCPU;
    char Octree_file_name[256];
    char *Octree_file_name_ptr;
    ifstream Octree_file;    

    /* On primary processor, determine name of octree input data file name. */

    if (CFFC_Primary_MPI_Processor()) {
       i = 0;
       while (1) {
          if (Input.Restart_File_Name[i] == ' ' ||
              Input.Restart_File_Name[i] == '.') break;
          Octree_file_name[i]=Input.Restart_File_Name[i];
          i = i + 1;
          if (i > strlen(Input.Restart_File_Name) ) break;
       } /* endwhile */
       Octree_file_name[i] = '\0';
       strcat(Octree_file_name, "_octree.tree");
       Octree_file_name_ptr = Octree_file_name;
    } /* endif */

    /* On primary processor, open the octree data file. */

    if (CFFC_Primary_MPI_Processor()) {
       Octree_file.open(Octree_file_name_ptr, ios::in);
       if (Octree_file.fail()) return (1);
    } /* endif */

    /* On primary processor, read in the data structure size parameters and 
       re-allocate memory as required. */

    if (CFFC_Primary_MPI_Processor()) {
       Octree_file.setf(ios::skipws);
       Octree_file >> nr >> nri >> nrj >> nrk >> ncpu >> nblk;
       Octree_file.unsetf(ios::skipws);

       Octree_DataStructure::Create_Octree_Data_Structure(Octree,
                                                          nri,
  	                                                  nrj,
  	                                                  nrk,
                                                          Input.AMR_IP.Number_of_Processors,
                                                          Input.AMR_IP.Number_of_Blocks_Per_Processor);
    } /* endif */

    /* Re-create and re-initialize the block resource list. */

    AdaptiveBlock3D_ResourceList::Create_Block_Resource_List(Global_Adaptive_Block_List,
                                                             Input.AMR_IP.Number_of_Processors, 
	  		                                     Input.AMR_IP.Number_of_Blocks_Per_Processor);

    /* On primary processor, read the octree data from the file. */

    if (CFFC_Primary_MPI_Processor()) {
       for (int nBlk = 0; nBlk <= Octree.NR-1; ++nBlk) {
	  Octree.Roots[nBlk].read(Octree_file,
                                  Global_Adaptive_Block_List);
       } /* endfor */
    } /* endif */

    /* Broadcast Octree data structure to all processors. */
       
    Octree_DataStructure::Broadcast_Octree_Data_Structure(Octree,
                                                          Global_Adaptive_Block_List);

    /* Set the maximum and minimum refinement levels. */

    Octree.MaximumRefinementLevel = Input.AMR_IP.Maximum_Refinement_Level-1;
    Octree.MinimumRefinementLevel = Input.AMR_IP.Minimum_Refinement_Level-1;

    /* Set the thresholds for refinement and coarsening of the mesh. */

    Octree.RefineThreshold = Input.AMR_IP.Threshold_for_Refinement;
    Octree.CoarsenThreshold = Input.AMR_IP.Threshold_for_Coarsening;

    /* Re-evaluate the octree block pointers. */

    Octree.assign_block_pointers();

    /* Re-allocate memory for local processor solution block list. */

    if (Local_Adaptive_Block_List.Nblk > 0) Local_Adaptive_Block_List.deallocate();
    Local_Adaptive_Block_List.allocate(Input.AMR_IP.Number_of_Blocks_Per_Processor);
    Local_Adaptive_Block_List.ThisCPU = Global_Adaptive_Block_List.ThisCPU;

    /* Reassign the root neighbour block information for all root blocks
       base on the global block numbers read in from the octree data
       file.   This routine will also copy block information to local 
       processor solution block list. */

    Octree_DataStructure::Reassign_Root_Neighbours(Octree,
                                                   Local_Adaptive_Block_List);

    /* Allocate memory for all message passing buffers used to send
       solution information between neighbouring solution blocks. */

    // Get the number of variables.
    for (int i_blk = 0 ; i_blk <= Local_Adaptive_Block_List.Nblk-1 ; ++i_blk) {
       if (Local_Adaptive_Block_List.Block[i_blk].used) {
          number_of_solution_variables = Local_Solution_Blocks.Soln_Blks[i_blk].NumVar();
          break;
       } /* endif */
    } /* endif */

    AdaptiveBlock3D_List::Allocate_Message_Buffers(Local_Adaptive_Block_List,
                                                   number_of_solution_variables);

    /* On primary processor, close octree data file. */

    if (CFFC_Primary_MPI_Processor()) Octree_file.close();

    /* Reading of octree data file complete.  Return zero value. */

    return(0);

}

/********************************************************
 * Routine: Write_Octree                                *
 *                                                      *
 * Writes the Octree data structure to a file for       *
 * later retrieval.                                     *
 *                                                      *
 ********************************************************/
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
int Write_Octree(Octree_DataStructure                        &Octree,
                 Input_Parameters<SOLN_pSTATE, SOLN_cSTATE>  &Input) {

    int i;
    char Octree_file_name[256];
    char *Octree_file_name_ptr;
    ofstream Octree_file;    

    /* On primary processor, determine name of octree output data file name. */

    if (CFFC_Primary_MPI_Processor()) {
       i = 0;
       while (1) {
          if (Input.Restart_File_Name[i] == ' ' ||
              Input.Restart_File_Name[i] == '.') break;
          Octree_file_name[i]=Input.Restart_File_Name[i];
          i = i + 1;
          if (i > strlen(Input.Restart_File_Name) ) break;
       } /* endwhile */
       Octree_file_name[i] = '\0';
       strcat(Octree_file_name, "_octree.tree");
       Octree_file_name_ptr = Octree_file_name;
    } /* endif */
    /* On primary processor, ppen the Octree data file. */

    if (CFFC_Primary_MPI_Processor()) {
       Octree_file.open(Octree_file_name_ptr, ios::out);
       if (Octree_file.fail()) return (1);
    } /* endif */

    /* On primary processor, write the octree data to the file. */

    if (CFFC_Primary_MPI_Processor()) {
      Octree_file << Octree;
    }
    /* On primary processor, close Octree data file. */

    if (CFFC_Primary_MPI_Processor()) Octree_file.close();

    /* Writing of octree data file complete.  Return zero value. */

    return(0);

}

/**********************************************************
 * Routine: Flag_Blocks_For_Refinement                    *
 *                                                        *
 * This routine flags the all local solution blocks for   *
 * refinement (coarsening or division hexarilateral       *
 * mesh solution blocks.                                  *
 *                                                        *
 **********************************************************/
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
void Flag_Blocks_For_Refinement(Hexa_Multi_Block<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> > &Local_Solution_Blocks,
                                Octree_DataStructure                                    &Octree,
                                AdaptiveBlock3D_ResourceList                            &Global_Adaptive_Block_List,
                                AdaptiveBlock3D_List                                    &Local_Adaptive_Block_List) {

  cout<<"\nError Flag_Blocks_For_Refinement is not written for 3D";

}

/**********************************************************
 * Routine: Refine_Grid                                   *
 *                                                        *
 * Performs the mesh refinement of the adaptive           *
 * blocks.  Returns a zero value if no error in the mesh  *
 * refinement precedure occurred.                         *
 *                                                        *
 **********************************************************/
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
int Refine_Grid(Hexa_Multi_Block<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> > &Local_Solution_Blocks,
                Octree_DataStructure                                    &Octree,
                AdaptiveBlock3D_ResourceList                            &Global_Adaptive_Block_List,
                AdaptiveBlock3D_List                                    &Local_Adaptive_Block_List) {

    cout<<"\nError Refine_Grid is not written for 3D";
    return(0);

}

/**********************************************************
 * Routine: Coarsen_Grid                                  *
 *                                                        *
 * Performs the mesh coarsening of the adaptive           *
 * blocks.  Returns a zero value if no error in the mesh  *
 * coarsening precedure occurred.                         *
 *                                                        *
 **********************************************************/
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
int Coarsen_Grid(Hexa_Multi_Block<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> > &Local_Solution_Blocks,
                 Octree_DataStructure                                    &Octree,
                 AdaptiveBlock3D_ResourceList                            &Global_Adaptive_Block_List,
                 AdaptiveBlock3D_List                                    &Local_Adaptive_Block_List) {
  
    cout<<"\nError Coarsen_Grid is not written for 3D";
    return(0);

}

/********************************************************
 * Routine: Fix_Refined_Block_Boundaries                *
 *                                                      *
 * Adjusts the locations of the boundary nodes of a     *
 * solution block so that the new node locations match  *
 * with cell volumes of adjacent solution blocks that   *
 * have lower levels of mesh refinement (i.e., are      *
 * coarser solution blocks).                            *
 *                                                      *
 ********************************************************/
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
void Fix_Refined_Block_Boundaries(Hexa_Multi_Block<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> > &Local_Solution_Blocks,
                                  AdaptiveBlock3D_List                                    &Soln_Block_List) {

  cout<<"\nError Fix_Refined_Block_Boundaries is not written for 3D";
 
}

/********************************************************
 * Routine: Unfix_Refined_Block_Boundaries              *
 *                                                      *
 * Returns the locations of the boundary nodes of a     *
 * solution block to their original unmodified          *
 * positions.                                           *
 *                                                      *
 ********************************************************/
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
void Unfix_Refined_Block_Boundaries(Hexa_Multi_Block<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> > &Local_Solution_Blocks,
                                    AdaptiveBlock3D_List                                    &Soln_Block_List) {

  cout<<"\nError Unfix_Refined_Block_Boundaries is not written for 3D";
 
}

/**********************************************************
 * Routine: AMR                                           *
 *                                                        *
 * Performs the adaptive mesh refinement (AMR) of the     *
 * adaptive solution blocks.  Returns a zero value if no  *
 * error in the AMR precedure occurred.                   *
 *                                                        *
 **********************************************************/
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
int AMR(Hexa_Multi_Block<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> > &Local_Solution_Blocks,
        Octree_DataStructure                                    &Octree,
        AdaptiveBlock3D_ResourceList                            &Global_Adaptive_Block_List,
        AdaptiveBlock3D_List                                    &Local_Adaptive_Block_List,
        const int                                               Set_New_Refinement_Flags) {

    cout<<"\nError AMR( is not written for 3D";
    return(0);

}

/**********************************************************
 * Routine: Initial_AMR                                   *
 *                                                        *
 * Performs initial refinement of the adaptive solution   *
 * block mesh based on initial data.  Returns a zero      *
 * value if no error in the AMR precedure has occurred.   *
 *                                                        *
 **********************************************************/
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
int Initial_AMR(Hexa_Multi_Block<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> > &Local_Solution_Blocks,
                Input_Parameters<SOLN_pSTATE, SOLN_cSTATE>              &Input,
                Octree_DataStructure                                    &Octree,
                AdaptiveBlock3D_ResourceList                            &Global_Adaptive_Block_List,
                AdaptiveBlock3D_List                                    &Local_Adaptive_Block_List) {

     cout<<"\nError Initial_AMR is not written for 3D";
     return(0);

}

/**********************************************************
 * Routine: Uniform_AMR                                   *
 *                                                        *
 * Performs uniform refinement of the adaptive solution   *
 * block mesh.  Returns a zero value if no error in the   *
 * AMR precedure has occurred.                            *
 *                                                        *
 **********************************************************/
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
int Uniform_AMR(Hexa_Multi_Block<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> > &Local_Solution_Blocks,
                Input_Parameters<SOLN_pSTATE, SOLN_cSTATE>              &Input,
                Octree_DataStructure                                    &Octree,
                AdaptiveBlock3D_ResourceList                            &Global_Adaptive_Block_List,
                AdaptiveBlock3D_List                                    &Local_Adaptive_Block_List) {

    cout<<"\nError Uniform_AMR is not written for 3D";
    return(0);
 
}

/**********************************************************
 * Routine: Boundary_AMR                                  *
 *                                                        *
 * Performs refinement of the adaptive solution block     *
 * mesh based on boundary conditions data.  Returns a     *
 * zero value if no error in the AMR precedure has        *
 * occurred.                                              *
 *                                                        *
 **********************************************************/
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
int Boundary_AMR(Hexa_Multi_Block<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> > &Local_Solution_Blocks,
		 Input_Parameters<SOLN_pSTATE, SOLN_cSTATE>              &Input,
		 Octree_DataStructure                                    &Octree,
		 AdaptiveBlock3D_ResourceList                            &Global_Adaptive_Block_List,
		 AdaptiveBlock3D_List                                    &Local_Adaptive_Block_List) {

    cout<<"\nError Boundary_AMR is not written for 3D";
    return(0);

}

#endif // _AMR_INCLUDED

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

#ifndef  _SYSTEM_LINUX_INCLUDED
#include "../System/System_Linux.h"
#endif //_SYSTEM_LINUX_INCLUDED


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
   
   int n_cpu, n_blk;
   int nused=0;

   
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
   
   /* Create (allocate) array of local 3D hexarilateral solution blocks. */
   
   Local_Solution_Blocks.Allocate(Input.AMR_IP.Number_of_Blocks_Per_Processor);

   AdaptiveBlock3D_ResourceList::Create_Block_Resource_List(Global_Adaptive_Block_List,
                                                            Input.AMR_IP.Number_of_Processors,
                                                            Input.AMR_IP.Number_of_Blocks_Per_Processor);
   Local_Adaptive_Block_List.allocate(Input.AMR_IP.Number_of_Blocks_Per_Processor);
   Local_Adaptive_Block_List.ThisCPU = Global_Adaptive_Block_List.ThisCPU;
   
   /* Loop over all initial mesh blocks and assign Octree root blocks and
      local solution block information as required. */
   MeshBLKs_to_LocalBLKs   *blknumber_convt;

   blknumber_convt = new MeshBLKs_to_LocalBLKs [Initial_Mesh.NBlk];
   //  cout<<"\n   CFFC_MPI::This_Processor_Number  "<< CFFC_MPI::This_Processor_Number<<endl; 
   for ( int nb = 0 ; nb <Initial_Mesh.NBlk ; ++nb ){
      if (Initial_Mesh.Grid_Blks[nb].Allocated) { // Mesh block is used!!!!
         // Get next free solution block from list of available
         // solution blocks.
         if (Global_Adaptive_Block_List.Nfree > 0) {
            n_cpu = Global_Adaptive_Block_List.nextCPU();
            n_blk = Global_Adaptive_Block_List.nextBlock();
            Global_Adaptive_Block_List.update_next();
            
            blknumber_convt[nb].cpu =  n_cpu;
            blknumber_convt[nb].blknum = n_blk;
            
          /*   for ( int iProc = 0; iProc !=  CFFC_MPI::Number_of_Processors; ++iProc ) { */
/*                if (  CFFC_MPI::This_Processor_Number == iProc ) { */
/*                   cout<<"\n CFFC_MPI::This_Processor_Number = "<< CFFC_MPI::This_Processor_Number<<endl; */
/*                   cout<<"\n nblk = "<<nb<<" n_cpu = "<<n_cpu<<" n_blk "<<n_blk; */
/*                   cout<<"\n   blknumber_convt["<<nb<<"].cpu ="<< blknumber_convt[nb].cpu<<" blknumber= "<<blknumber_convt[nb].blknum<<endl; */
                  
/*                   System::sleep(0.1); */
/*                } */
/*                MPI::COMM_WORLD.Barrier(); */
/*             } */

         } else {
            cout << "\n"
                 << " AMR Error: Create_Initial_Solution_Blocks, Insufficient number of hexahedrial solution blocks.\n";
            return (1);
         } /* endif */
      }
   }
   
   Global_Adaptive_Block_List.initialize();
   
   
   for ( int nb = 0 ; nb <Initial_Mesh.NBlk ; ++nb ){
      // cout<<"\n nb = "<<nb<<"  allocated is "<<Initial_Mesh.Grid_Blks[nb].Allocated<<endl;
      if (Initial_Mesh.Grid_Blks[nb].Allocated) { // Mesh block is used!!!!
         
         // Get next free solution block from list of available
         // solution blocks.
         if (Global_Adaptive_Block_List.Nfree > 0) {
            n_cpu = Global_Adaptive_Block_List.nextCPU();
            n_blk = Global_Adaptive_Block_List.nextBlock();
            Global_Adaptive_Block_List.update_next();
         } else {
            cout << "\n"
                 << " AMR Error: Create_Initial_Solution_Blocks, Insufficient number of hexahedrial solution blocks.\n";
            return (1);
         } /* endif */
	 
         // Assign block information to appropriate Octree root solution block.
         Octree.Roots[nb].block.used = ADAPTIVEBLOCK3D_USED;
         Octree.Roots[nb].block.gblknum =
            Global_Adaptive_Block_List.Nused-1;
         //Number of neighbouring grid blocks in each direction
         //for each of the block in the roots
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
         Octree.Roots[nb].block.nBNW  = Initial_Mesh.Connectivity[nb].num_neighBNW;
         Octree.Roots[nb].block.nBSW = Initial_Mesh.Connectivity[nb].num_neighBSW;
         Octree.Roots[nb].block.nBNE = Initial_Mesh.Connectivity[nb].num_neighBNE;
         Octree.Roots[nb].block.nBNW = Initial_Mesh.Connectivity[nb].num_neighBNW;

         Octree.Roots[nb].block.info.cpu = n_cpu;
         Octree.Roots[nb].block.info.blknum = n_blk;
         Octree.Roots[nb].block.info.dimen.i =
            Initial_Mesh.Grid_Blks[nb].NCi -
            2*Initial_Mesh.Grid_Blks[nb].Nghost;
         Octree.Roots[nb].block.info.dimen.j =
            Initial_Mesh.Grid_Blks[nb].NCj -
            2*Initial_Mesh.Grid_Blks[nb].Nghost;
         Octree.Roots[nb].block.info.dimen.k =
            Initial_Mesh.Grid_Blks[nb].NCk -
            2*Initial_Mesh.Grid_Blks[nb].Nghost;
         Octree.Roots[nb].block.info.dimen.ghost = 2;
         Octree.Roots[nb].block.info.sector = ADAPTIVEBLOCK3D_SECTOR_NONE;
         Octree.Roots[nb].block.info.level = 0;
              
         // block number of the neighbouring blocks in each direction
         // for each of the block in the roots
         if(Initial_Mesh.Connectivity[nb].neighT>=0){
            Octree.Roots[nb].block.infoT[0].blknum = blknumber_convt[Initial_Mesh.Connectivity[nb].neighT].blknum;
         }
         if(Initial_Mesh.Connectivity[nb].neighB>=0){
            Octree.Roots[nb].block.infoB[0].blknum = blknumber_convt[Initial_Mesh.Connectivity[nb].neighB].blknum;
         }
         
         if(Initial_Mesh.Connectivity[nb].neighN>=0){
            Octree.Roots[nb].block.infoN[0].blknum = blknumber_convt[Initial_Mesh.Connectivity[nb].neighN].blknum;
         }
         
         if(Initial_Mesh.Connectivity[nb].neighS>=0){
            
            Octree.Roots[nb].block.infoS[0].blknum = blknumber_convt[Initial_Mesh.Connectivity[nb].neighS].blknum;
         }
         if(Initial_Mesh.Connectivity[nb].neighE>=0){
            
            Octree.Roots[nb].block.infoE[0].blknum = blknumber_convt[Initial_Mesh.Connectivity[nb].neighE].blknum;
         }
         if(Initial_Mesh.Connectivity[nb].neighW>=0){
            Octree.Roots[nb].block.infoW[0].blknum = blknumber_convt[Initial_Mesh.Connectivity[nb].neighW].blknum;
         }
            
         if(Initial_Mesh.Connectivity[nb].neighNW>=0){
            Octree.Roots[nb].block.infoNW[0].blknum = blknumber_convt[Initial_Mesh.Connectivity[nb].neighNW[0]].blknum;
         }
         if(Initial_Mesh.Connectivity[nb].neighNE>=0){
            Octree.Roots[nb].block.infoNE[0].blknum = blknumber_convt[Initial_Mesh.Connectivity[nb].neighNE[0]].blknum;
         }
         
         if(Initial_Mesh.Connectivity[nb].neighSE>=0){
            Octree.Roots[nb].block.infoSE[0].blknum = blknumber_convt[Initial_Mesh.Connectivity[nb].neighSE[0]].blknum;
         }
         if(Initial_Mesh.Connectivity[nb].neighSW>=0){
            Octree.Roots[nb].block.infoSW[0].blknum = blknumber_convt[Initial_Mesh.Connectivity[nb].neighSW[0]].blknum;
         }
         if(Initial_Mesh.Connectivity[nb].neighTN>=0){
            Octree.Roots[nb].block.infoTN[0].blknum = blknumber_convt[Initial_Mesh.Connectivity[nb].neighTN[0]].blknum;
         }
         if(Initial_Mesh.Connectivity[nb].neighTS>=0){
            Octree.Roots[nb].block.infoTS[0].blknum = blknumber_convt[Initial_Mesh.Connectivity[nb].neighTS[0]].blknum;
         }
         
         if(Initial_Mesh.Connectivity[nb].neighTE>=0){
            Octree.Roots[nb].block.infoTE[0].blknum = blknumber_convt[Initial_Mesh.Connectivity[nb].neighTE[0]].blknum;
         }
         if(Initial_Mesh.Connectivity[nb].neighTW>=0){
            Octree.Roots[nb].block.infoTW[0].blknum = blknumber_convt[Initial_Mesh.Connectivity[nb].neighTW[0]].blknum;
         }
         if(Initial_Mesh.Connectivity[nb].neighBN>=0){
            Octree.Roots[nb].block.infoBN[0].blknum = blknumber_convt[Initial_Mesh.Connectivity[nb].neighBN[0]].blknum;
         }
         if(Initial_Mesh.Connectivity[nb].neighBS>=0){
            Octree.Roots[nb].block.infoBS[0].blknum = blknumber_convt[Initial_Mesh.Connectivity[nb].neighBS[0]].blknum;
         }
         if(Initial_Mesh.Connectivity[nb].neighBE>=0){
            Octree.Roots[nb].block.infoBE[0].blknum = blknumber_convt[Initial_Mesh.Connectivity[nb].neighBE[0]].blknum;
         }
         if(Initial_Mesh.Connectivity[nb].neighBW>=0){
            Octree.Roots[nb].block.infoBW[0].blknum = blknumber_convt[Initial_Mesh.Connectivity[nb].neighBW[0]].blknum;
         }
         if(Initial_Mesh.Connectivity[nb].neighTNW>=0){
            Octree.Roots[nb].block.infoTNW[0].blknum = blknumber_convt[Initial_Mesh.Connectivity[nb].neighTNW[0]].blknum;
         }
         if(Initial_Mesh.Connectivity[nb].neighTSW>=0){
            Octree.Roots[nb].block.infoTSW[0].blknum = blknumber_convt[Initial_Mesh.Connectivity[nb].neighTSW[0]].blknum;
         }
         if(Initial_Mesh.Connectivity[nb].neighTNE>=0){
            Octree.Roots[nb].block.infoTNE[0].blknum = blknumber_convt[Initial_Mesh.Connectivity[nb].neighTNE[0]].blknum;
         }
         if(Initial_Mesh.Connectivity[nb].neighTSE>=0){
            Octree.Roots[nb].block.infoTSE[0].blknum = blknumber_convt[Initial_Mesh.Connectivity[nb].neighTSE[0]].blknum;
         }
         if(Initial_Mesh.Connectivity[nb].neighBNW>=0){
            Octree.Roots[nb].block.infoBNW[0].blknum = blknumber_convt[Initial_Mesh.Connectivity[nb].neighBNW[0]].blknum;
         }
         if(Initial_Mesh.Connectivity[nb].neighBSW>=0){
            Octree.Roots[nb].block.infoBSW[0].blknum = blknumber_convt[Initial_Mesh.Connectivity[nb].neighBSW[0]].blknum;
         }
         if(Initial_Mesh.Connectivity[nb].neighBNE>=0){
            Octree.Roots[nb].block.infoBNE[0].blknum = blknumber_convt[Initial_Mesh.Connectivity[nb].neighBNE[0]].blknum;
         }
         if(Initial_Mesh.Connectivity[nb].neighBSE>=0){
            Octree.Roots[nb].block.infoBSE[0].blknum = blknumber_convt[Initial_Mesh.Connectivity[nb].neighBSE[0]].blknum;
         }
                     
         if(Initial_Mesh.Connectivity[nb].neighT>=0){
            Octree.Roots[nb].block.infoT[0].cpu = blknumber_convt[Initial_Mesh.Connectivity[nb].neighT].cpu;
         }
         if(Initial_Mesh.Connectivity[nb].neighB>=0){
            Octree.Roots[nb].block.infoB[0].cpu = blknumber_convt[Initial_Mesh.Connectivity[nb].neighB].cpu;
         }
         if(Initial_Mesh.Connectivity[nb].neighN>=0){
            Octree.Roots[nb].block.infoN[0].cpu = blknumber_convt[Initial_Mesh.Connectivity[nb].neighN].cpu;
         }
         if(Initial_Mesh.Connectivity[nb].neighS>=0){
            Octree.Roots[nb].block.infoS[0].cpu = blknumber_convt[Initial_Mesh.Connectivity[nb].neighS].cpu;
         }
         if(Initial_Mesh.Connectivity[nb].neighE>=0){
            Octree.Roots[nb].block.infoE[0].cpu = blknumber_convt[Initial_Mesh.Connectivity[nb].neighE].cpu;
         }
         if(Initial_Mesh.Connectivity[nb].neighW>=0){
            Octree.Roots[nb].block.infoW[0].cpu = blknumber_convt[Initial_Mesh.Connectivity[nb].neighW].cpu;

         /*    for ( int iProc = 0; iProc !=  CFFC_MPI::Number_of_Processors; ++iProc ) { */
/*                if (  CFFC_MPI::This_Processor_Number == iProc ) { */
/*                   cout<<"\n CFFC_MPI::This_Processor_Number = "<< CFFC_MPI::This_Processor_Number<<endl; */
/*                   cout<<"\n intial mesh's block "<<nb<<"'s west neigbor = "<<Initial_Mesh.Connectivity[nb].neighW<<endl; */
                  
/*                   cout<<"\n  Octree.Roots["<<nb<<"].block.infoW[0].cpu "<< Octree.Roots[nb].block.infoW[0].cpu; */
              
                  
/*                   System::sleep(0.1); */
/*                } */
/*                MPI::COMM_WORLD.Barrier(); */
/*             } */
            
         }
                  
         
         if(Initial_Mesh.Connectivity[nb].neighNW>=0){
            Octree.Roots[nb].block.infoNW[0].cpu = blknumber_convt[Initial_Mesh.Connectivity[nb].neighNW[0]].cpu;
         }
         if(Initial_Mesh.Connectivity[nb].neighNE>=0){
            Octree.Roots[nb].block.infoNE[0].cpu = blknumber_convt[Initial_Mesh.Connectivity[nb].neighNE[0]].cpu;
         }
         if(Initial_Mesh.Connectivity[nb].neighSE>=0){
            Octree.Roots[nb].block.infoSE[0].cpu = blknumber_convt[Initial_Mesh.Connectivity[nb].neighSE[0]].cpu;
         }
         if(Initial_Mesh.Connectivity[nb].neighSW>=0){
            Octree.Roots[nb].block.infoSW[0].cpu = blknumber_convt[Initial_Mesh.Connectivity[nb].neighSW[0]].cpu;
         }
         if(Initial_Mesh.Connectivity[nb].neighTN>=0){
            Octree.Roots[nb].block.infoTN[0].cpu = blknumber_convt[Initial_Mesh.Connectivity[nb].neighTN[0]].cpu;
         }
         if(Initial_Mesh.Connectivity[nb].neighTS>=0){
            Octree.Roots[nb].block.infoTS[0].cpu = blknumber_convt[Initial_Mesh.Connectivity[nb].neighTS[0]].cpu;
         }
         if(Initial_Mesh.Connectivity[nb].neighTE>=0){
            Octree.Roots[nb].block.infoTE[0].cpu = blknumber_convt[Initial_Mesh.Connectivity[nb].neighTE[0]].cpu;
         }
         if(Initial_Mesh.Connectivity[nb].neighTW>=0){
            Octree.Roots[nb].block.infoTW[0].cpu = blknumber_convt[Initial_Mesh.Connectivity[nb].neighTW[0]].cpu;
         }
         if(Initial_Mesh.Connectivity[nb].neighBN>=0){
            Octree.Roots[nb].block.infoBN[0].cpu = blknumber_convt[Initial_Mesh.Connectivity[nb].neighBN[0]].cpu;
         }
         if(Initial_Mesh.Connectivity[nb].neighBS>=0){
            Octree.Roots[nb].block.infoBS[0].cpu = blknumber_convt[Initial_Mesh.Connectivity[nb].neighBS[0]].cpu;
         }
         if(Initial_Mesh.Connectivity[nb].neighBE>=0){
            Octree.Roots[nb].block.infoBE[0].cpu = blknumber_convt[Initial_Mesh.Connectivity[nb].neighBE[0]].cpu;
         }
         if(Initial_Mesh.Connectivity[nb].neighBW>=0){
            Octree.Roots[nb].block.infoBW[0].cpu = blknumber_convt[Initial_Mesh.Connectivity[nb].neighBW[0]].cpu;
         }
         if(Initial_Mesh.Connectivity[nb].neighTNW>=0){
            
            Octree.Roots[nb].block.infoTNW[0].cpu = blknumber_convt[Initial_Mesh.Connectivity[nb].neighTNW[0]].cpu;
         }
         if(Initial_Mesh.Connectivity[nb].neighTSW>=0){
            Octree.Roots[nb].block.infoTSW[0].cpu = blknumber_convt[Initial_Mesh.Connectivity[nb].neighTSW[0]].cpu;
         }
         if(Initial_Mesh.Connectivity[nb].neighTNE>=0){
            Octree.Roots[nb].block.infoTNE[0].cpu = blknumber_convt[Initial_Mesh.Connectivity[nb].neighTNE[0]].cpu;
         }
         if(Initial_Mesh.Connectivity[nb].neighTSE>=0){
            Octree.Roots[nb].block.infoTSE[0].cpu = blknumber_convt[Initial_Mesh.Connectivity[nb].neighTSE[0]].cpu;
         }
         if(Initial_Mesh.Connectivity[nb].neighBNW>=0){
            Octree.Roots[nb].block.infoBNW[0].cpu = blknumber_convt[Initial_Mesh.Connectivity[nb].neighBNW[0]].cpu;
         }
         if(Initial_Mesh.Connectivity[nb].neighBSW>=0){
            Octree.Roots[nb].block.infoBSW[0].cpu = blknumber_convt[Initial_Mesh.Connectivity[nb].neighBSW[0]].cpu;
         }
         if(Initial_Mesh.Connectivity[nb].neighBNE>=0){
            Octree.Roots[nb].block.infoBNE[0].cpu = blknumber_convt[Initial_Mesh.Connectivity[nb].neighBNE[0]].cpu;
         }
         if(Initial_Mesh.Connectivity[nb].neighBSE>=0){
            Octree.Roots[nb].block.infoBSE[0].cpu = blknumber_convt[Initial_Mesh.Connectivity[nb].neighBSE[0]].cpu;
         }
                                                            
         // block orientation information of the neighbouring blocks in each direction
         // for each of the block in the roots
         Octree.Roots[nb].block.infoT[0].blkorient = Initial_Mesh.Connectivity[nb].neighT_info;
         Octree.Roots[nb].block.infoB[0].blkorient = Initial_Mesh.Connectivity[nb].neighB_info;
         Octree.Roots[nb].block.infoN[0].blkorient = Initial_Mesh.Connectivity[nb].neighN_info;
         Octree.Roots[nb].block.infoS[0].blkorient = Initial_Mesh.Connectivity[nb].neighS_info;
         Octree.Roots[nb].block.infoE[0].blkorient = Initial_Mesh.Connectivity[nb].neighE_info;
         Octree.Roots[nb].block.infoW[0].blkorient = Initial_Mesh.Connectivity[nb].neighW_info;
         
         Octree.Roots[nb].block.infoNW[0].blkorient = Initial_Mesh.Connectivity[nb].neighNW_info[0];
         Octree.Roots[nb].block.infoNE[0].blkorient = Initial_Mesh.Connectivity[nb].neighNE_info[0];
         Octree.Roots[nb].block.infoSE[0].blkorient = Initial_Mesh.Connectivity[nb].neighSE_info[0];
         Octree.Roots[nb].block.infoSW[0].blkorient = Initial_Mesh.Connectivity[nb].neighSW_info[0];
         Octree.Roots[nb].block.infoTN[0].blkorient = Initial_Mesh.Connectivity[nb].neighTN_info[0];
         Octree.Roots[nb].block.infoTS[0].blkorient = Initial_Mesh.Connectivity[nb].neighTS_info[0];
         Octree.Roots[nb].block.infoTE[0].blkorient = Initial_Mesh.Connectivity[nb].neighTE_info[0];
         Octree.Roots[nb].block.infoTW[0].blkorient = Initial_Mesh.Connectivity[nb].neighTW_info[0];
         Octree.Roots[nb].block.infoBN[0].blkorient = Initial_Mesh.Connectivity[nb].neighBN_info[0];
         Octree.Roots[nb].block.infoBS[0].blkorient = Initial_Mesh.Connectivity[nb].neighBS_info[0];
         Octree.Roots[nb].block.infoBE[0].blkorient = Initial_Mesh.Connectivity[nb].neighBE_info[0];
         Octree.Roots[nb].block.infoBW[0].blkorient = Initial_Mesh.Connectivity[nb].neighBW_info[0];
         
         Octree.Roots[nb].block.infoTNW[0].blkorient = Initial_Mesh.Connectivity[nb].neighTNW_info[0];
         Octree.Roots[nb].block.infoTSW[0].blkorient = Initial_Mesh.Connectivity[nb].neighTSW_info[0];
         Octree.Roots[nb].block.infoTNE[0].blkorient = Initial_Mesh.Connectivity[nb].neighTNE_info[0];
         Octree.Roots[nb].block.infoTSE[0].blkorient = Initial_Mesh.Connectivity[nb].neighTSE_info[0];
         Octree.Roots[nb].block.infoBNW[0].blkorient = Initial_Mesh.Connectivity[nb].neighBNW_info[0];
         Octree.Roots[nb].block.infoBSW[0].blkorient = Initial_Mesh.Connectivity[nb].neighBSW_info[0];
         Octree.Roots[nb].block.infoBNE[0].blkorient = Initial_Mesh.Connectivity[nb].neighBNE_info[0];
         Octree.Roots[nb].block.infoBSE[0].blkorient = Initial_Mesh.Connectivity[nb].neighBSE_info[0];

         /*    Octree.Roots[nb].parent_ptr = NULL; */
/*          Octree.Roots[nb].childTNW_ptr = NULL; */
/*          Octree.Roots[nb].childTNE_ptr = NULL; */
/*          Octree.Roots[nb].childTSE_ptr = NULL; */
/*          Octree.Roots[nb].childTSW_ptr = NULL; */
/*          Octree.Roots[nb].childBNW_ptr = NULL; */
/*          Octree.Roots[nb].childBNE_ptr = NULL; */
/*          Octree.Roots[nb].childBSE_ptr = NULL; */
/*          Octree.Roots[nb].childBSW_ptr = NULL; */
         
         Octree.Blocks[n_cpu][n_blk] = &(Octree.Roots[nb]); 
         
         // For solution blocks on this processor (or processor
         // element), add block to local list, create the solution
         // block, and copy the appropriate block of the
         // initial hexarilateral mesh to the solution block mesh.
         if (Global_Adaptive_Block_List.ThisCPU == n_cpu) {
            Local_Adaptive_Block_List.Block[n_blk] = Octree.Roots[nb].block;
            Local_Solution_Blocks.Soln_Blks[n_blk].Create_Block(Initial_Mesh.Grid_Blks[nb]);
            Local_Solution_Blocks.Soln_Blks[n_blk].Flow_Type = Input.i_Flow_Type;
            Local_Solution_Blocks.Block_Used[n_blk] = HEXA_BLOCK_USED;
         } /* endif */
      } else{
         Octree.Roots[nb].block.used = 0;
         nused++;
      } /* endif */
      
   } /* endfor */

   /* dimen, sector and level assignment. */
   for ( int nb = 0 ; nb <Initial_Mesh.NBlk ; ++nb ){
      if(Octree.Roots[nb].block.infoT[0].blknum>=0){
         Octree.Roots[nb].block.infoT[0].dimen =  Octree.Roots[Octree.Roots[nb].block.infoT[0].blknum].block.info.dimen;
         Octree.Roots[nb].block.infoT[0].sector =  Octree.Roots[Octree.Roots[nb].block.infoT[0].blknum].block.info.sector;
         Octree.Roots[nb].block.infoT[0].level =  Octree.Roots[Octree.Roots[nb].block.infoT[0].blknum].block.info.level;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoT[0].dimen =  Octree.Roots[nb].block.infoT[0].dimen;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoT[0].sector =  Octree.Roots[nb].block.infoT[0].sector;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoT[0].level =  Octree.Roots[nb].block.infoT[0].level;
         
      }
      
      if(Octree.Roots[nb].block.infoB[0].blknum>=0) {
         Octree.Roots[nb].block.infoB[0].dimen =  Octree.Roots[Octree.Roots[nb].block.infoB[0].blknum].block.info.dimen;
         Octree.Roots[nb].block.infoB[0].sector =  Octree.Roots[Octree.Roots[nb].block.infoB[0].blknum].block.info.sector;
         Octree.Roots[nb].block.infoB[0].level =  Octree.Roots[Octree.Roots[nb].block.infoB[0].blknum].block.info.level;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoB[0].dimen =  Octree.Roots[nb].block.infoB[0].dimen;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoB[0].sector =  Octree.Roots[nb].block.infoB[0].sector;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoB[0].level =  Octree.Roots[nb].block.infoB[0].level;
      }
      
      if(Octree.Roots[nb].block.infoN[0].blknum>=0){
         Octree.Roots[nb].block.infoN[0].dimen =  Octree.Roots[Octree.Roots[nb].block.infoN[0].blknum].block.info.dimen;
         Octree.Roots[nb].block.infoN[0].sector =  Octree.Roots[Octree.Roots[nb].block.infoN[0].blknum].block.info.sector;
         Octree.Roots[nb].block.infoN[0].level =  Octree.Roots[Octree.Roots[nb].block.infoN[0].blknum].block.info.level;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoN[0].dimen =  Octree.Roots[nb].block.infoN[0].dimen;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoN[0].sector =  Octree.Roots[nb].block.infoN[0].sector;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoN[0].level =  Octree.Roots[nb].block.infoN[0].level;
      }
      
      if(Octree.Roots[nb].block.infoS[0].blknum>=0){
         Octree.Roots[nb].block.infoS[0].dimen =  Octree.Roots[Octree.Roots[nb].block.infoS[0].blknum].block.info.dimen;
         Octree.Roots[nb].block.infoS[0].sector =  Octree.Roots[Octree.Roots[nb].block.infoS[0].blknum].block.info.sector;
         Octree.Roots[nb].block.infoS[0].level =  Octree.Roots[Octree.Roots[nb].block.infoS[0].blknum].block.info.level;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoS[0].dimen =  Octree.Roots[nb].block.infoS[0].dimen;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoS[0].sector =  Octree.Roots[nb].block.infoS[0].sector;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoS[0].level =  Octree.Roots[nb].block.infoS[0].level;

      }


      
      if(Octree.Roots[nb].block.infoE[0].blknum>=0){
         Octree.Roots[nb].block.infoE[0].dimen =  Octree.Roots[Octree.Roots[nb].block.infoE[0].blknum].block.info.dimen;
         Octree.Roots[nb].block.infoE[0].sector =  Octree.Roots[Octree.Roots[nb].block.infoE[0].blknum].block.info.sector;
         Octree.Roots[nb].block.infoE[0].level =  Octree.Roots[Octree.Roots[nb].block.infoE[0].blknum].block.info.level;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoE[0].dimen =  Octree.Roots[nb].block.infoE[0].dimen;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoE[0].sector =  Octree.Roots[nb].block.infoE[0].sector;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoE[0].level =  Octree.Roots[nb].block.infoE[0].level;
      }
      
      if(Octree.Roots[nb].block.infoW[0].blknum>=0) {
         Octree.Roots[nb].block.infoW[0].dimen =  Octree.Roots[Octree.Roots[nb].block.infoW[0].blknum].block.info.dimen;
         Octree.Roots[nb].block.infoW[0].sector =  Octree.Roots[Octree.Roots[nb].block.infoW[0].blknum].block.info.sector;
         Octree.Roots[nb].block.infoW[0].level =  Octree.Roots[Octree.Roots[nb].block.infoW[0].blknum].block.info.level;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoW[0].dimen =  Octree.Roots[nb].block.infoW[0].dimen;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoW[0].sector =  Octree.Roots[nb].block.infoW[0].sector;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoW[0].level =  Octree.Roots[nb].block.infoW[0].level;
      }
      
      if(Octree.Roots[nb].block.infoNW[0].blknum>=0){
         Octree.Roots[nb].block.infoNW[0].dimen =  Octree.Roots[Octree.Roots[nb].block.infoNW[0].blknum].block.info.dimen;
         Octree.Roots[nb].block.infoNW[0].sector =  Octree.Roots[Octree.Roots[nb].block.infoNW[0].blknum].block.info.sector;
         Octree.Roots[nb].block.infoNW[0].level =  Octree.Roots[Octree.Roots[nb].block.infoNW[0].blknum].block.info.level;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoNW[0].dimen =  Octree.Roots[nb].block.infoNW[0].dimen;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoNW[0].sector =  Octree.Roots[nb].block.infoNW[0].sector;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoNW[0].level =  Octree.Roots[nb].block.infoNW[0].level;
      }
      
      if(Octree.Roots[nb].block.infoNE[0].blknum>=0){
         Octree.Roots[nb].block.infoNE[0].dimen =  Octree.Roots[Octree.Roots[nb].block.infoNE[0].blknum].block.info.dimen;
         Octree.Roots[nb].block.infoNE[0].sector =  Octree.Roots[Octree.Roots[nb].block.infoNE[0].blknum].block.info.sector;
         Octree.Roots[nb].block.infoNE[0].level =  Octree.Roots[Octree.Roots[nb].block.infoNE[0].blknum].block.info.level;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoNE[0].dimen =  Octree.Roots[nb].block.infoNE[0].dimen;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoNE[0].sector =  Octree.Roots[nb].block.infoNE[0].sector;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoNE[0].level =  Octree.Roots[nb].block.infoNE[0].level;
      }
      
      if(Octree.Roots[nb].block.infoSE[0].blknum>=0) {
         Octree.Roots[nb].block.infoSE[0].dimen =  Octree.Roots[Octree.Roots[nb].block.infoSE[0].blknum].block.info.dimen;
         Octree.Roots[nb].block.infoSE[0].sector =  Octree.Roots[Octree.Roots[nb].block.infoSE[0].blknum].block.info.sector;
         Octree.Roots[nb].block.infoSE[0].level =  Octree.Roots[Octree.Roots[nb].block.infoSE[0].blknum].block.info.level;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoSE[0].dimen =  Octree.Roots[nb].block.infoSE[0].dimen;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoSE[0].sector =  Octree.Roots[nb].block.infoSE[0].sector;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoSE[0].level =  Octree.Roots[nb].block.infoSE[0].level;
      }
      
      if(Octree.Roots[nb].block.infoSW[0].blknum>=0){
         Octree.Roots[nb].block.infoSW[0].dimen =  Octree.Roots[Octree.Roots[nb].block.infoSW[0].blknum].block.info.dimen;
         Octree.Roots[nb].block.infoSW[0].sector =  Octree.Roots[Octree.Roots[nb].block.infoSW[0].blknum].block.info.sector;
         Octree.Roots[nb].block.infoSW[0].level =  Octree.Roots[Octree.Roots[nb].block.infoSW[0].blknum].block.info.level;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoSW[0].dimen =  Octree.Roots[nb].block.infoSW[0].dimen;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoSW[0].sector =  Octree.Roots[nb].block.infoSW[0].sector;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoSW[0].level =  Octree.Roots[nb].block.infoSW[0].level;
      }
      
      if(Octree.Roots[nb].block.infoTN[0].blknum>=0){
         Octree.Roots[nb].block.infoTN[0].dimen =  Octree.Roots[Octree.Roots[nb].block.infoTN[0].blknum].block.info.dimen;
         Octree.Roots[nb].block.infoTN[0].sector =  Octree.Roots[Octree.Roots[nb].block.infoTN[0].blknum].block.info.sector;
         Octree.Roots[nb].block.infoTN[0].level =  Octree.Roots[Octree.Roots[nb].block.infoTN[0].blknum].block.info.level;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoTN[0].dimen =  Octree.Roots[nb].block.infoTN[0].dimen;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoTN[0].sector =  Octree.Roots[nb].block.infoTN[0].sector;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoTN[0].level =  Octree.Roots[nb].block.infoTN[0].level;
      }
      
      if(Octree.Roots[nb].block.infoTS[0].blknum>=0){
         Octree.Roots[nb].block.infoTS[0].dimen =  Octree.Roots[Octree.Roots[nb].block.infoTS[0].blknum].block.info.dimen;
         Octree.Roots[nb].block.infoTS[0].sector =  Octree.Roots[Octree.Roots[nb].block.infoTS[0].blknum].block.info.sector;
         Octree.Roots[nb].block.infoTS[0].level =  Octree.Roots[Octree.Roots[nb].block.infoTS[0].blknum].block.info.level;
     Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoTS[0].dimen =  Octree.Roots[nb].block.infoTS[0].dimen;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoTS[0].sector =  Octree.Roots[nb].block.infoTS[0].sector;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoTS[0].level =  Octree.Roots[nb].block.infoTS[0].level;
      }
      
      if(Octree.Roots[nb].block.infoTE[0].blknum>=0){
         Octree.Roots[nb].block.infoTE[0].dimen =  Octree.Roots[Octree.Roots[nb].block.infoTE[0].blknum].block.info.dimen;
         Octree.Roots[nb].block.infoTE[0].sector =  Octree.Roots[Octree.Roots[nb].block.infoTE[0].blknum].block.info.sector;
         Octree.Roots[nb].block.infoTE[0].level =  Octree.Roots[Octree.Roots[nb].block.infoTE[0].blknum].block.info.level;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoTE[0].dimen =  Octree.Roots[nb].block.infoTE[0].dimen;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoTE[0].sector =  Octree.Roots[nb].block.infoTE[0].sector;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoTE[0].level =  Octree.Roots[nb].block.infoTE[0].level;
      }
      
      if(Octree.Roots[nb].block.infoTW[0].blknum>=0){
         Octree.Roots[nb].block.infoTW[0].dimen =  Octree.Roots[Octree.Roots[nb].block.infoTW[0].blknum].block.info.dimen;
         Octree.Roots[nb].block.infoTW[0].sector =  Octree.Roots[Octree.Roots[nb].block.infoTW[0].blknum].block.info.sector;
         Octree.Roots[nb].block.infoTW[0].level =  Octree.Roots[Octree.Roots[nb].block.infoTW[0].blknum].block.info.level;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoTW[0].dimen =  Octree.Roots[nb].block.infoTW[0].dimen;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoTW[0].sector =  Octree.Roots[nb].block.infoTW[0].sector;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoTW[0].level =  Octree.Roots[nb].block.infoTW[0].level;

      }
      
      if(Octree.Roots[nb].block.infoBN[0].blknum>=0){
         Octree.Roots[nb].block.infoBN[0].dimen =  Octree.Roots[Octree.Roots[nb].block.infoBN[0].blknum].block.info.dimen;
         Octree.Roots[nb].block.infoBN[0].sector =  Octree.Roots[Octree.Roots[nb].block.infoBN[0].blknum].block.info.sector;
         Octree.Roots[nb].block.infoBN[0].level =  Octree.Roots[Octree.Roots[nb].block.infoBN[0].blknum].block.info.level;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoBN[0].dimen =  Octree.Roots[nb].block.infoBN[0].dimen;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoBN[0].sector =  Octree.Roots[nb].block.infoBN[0].sector;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoBN[0].level =  Octree.Roots[nb].block.infoBN[0].level;
      }
      
      if(Octree.Roots[nb].block.infoBS[0].blknum>=0){
         Octree.Roots[nb].block.infoBS[0].dimen =  Octree.Roots[Octree.Roots[nb].block.infoBS[0].blknum].block.info.dimen;
         Octree.Roots[nb].block.infoBS[0].sector =  Octree.Roots[Octree.Roots[nb].block.infoBS[0].blknum].block.info.sector;
         Octree.Roots[nb].block.infoBS[0].level =  Octree.Roots[Octree.Roots[nb].block.infoBS[0].blknum].block.info.level;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoBS[0].dimen =  Octree.Roots[nb].block.infoBS[0].dimen;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoBS[0].sector =  Octree.Roots[nb].block.infoBS[0].sector;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoBS[0].level =  Octree.Roots[nb].block.infoBS[0].level;
      }
      
      if(Octree.Roots[nb].block.infoBE[0].blknum>=0){
         Octree.Roots[nb].block.infoBE[0].dimen =  Octree.Roots[Octree.Roots[nb].block.infoBE[0].blknum].block.info.dimen;
         Octree.Roots[nb].block.infoBE[0].sector =  Octree.Roots[Octree.Roots[nb].block.infoBE[0].blknum].block.info.sector;
         Octree.Roots[nb].block.infoBE[0].level =  Octree.Roots[Octree.Roots[nb].block.infoBE[0].blknum].block.info.level;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoBE[0].dimen =  Octree.Roots[nb].block.infoBE[0].dimen;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoBE[0].sector =  Octree.Roots[nb].block.infoBE[0].sector;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoBE[0].level =  Octree.Roots[nb].block.infoBE[0].level;
         
      }
      
      if(Octree.Roots[nb].block.infoBW[0].blknum>=0){
         Octree.Roots[nb].block.infoBW[0].dimen =  Octree.Roots[Octree.Roots[nb].block.infoBW[0].blknum].block.info.dimen;
         Octree.Roots[nb].block.infoBW[0].sector =  Octree.Roots[Octree.Roots[nb].block.infoBW[0].blknum].block.info.sector;
         Octree.Roots[nb].block.infoBW[0].level =  Octree.Roots[Octree.Roots[nb].block.infoBW[0].blknum].block.info.level;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoBW[0].dimen =  Octree.Roots[nb].block.infoBW[0].dimen;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoBW[0].sector =  Octree.Roots[nb].block.infoBW[0].sector;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoBW[0].level =  Octree.Roots[nb].block.infoBW[0].level;
      }
      
      if(Octree.Roots[nb].block.infoTNW[0].blknum>=0){
         Octree.Roots[nb].block.infoTNW[0].dimen =  Octree.Roots[Octree.Roots[nb].block.infoTNW[0].blknum].block.info.dimen;   
         Octree.Roots[nb].block.infoTNW[0].sector =  Octree.Roots[Octree.Roots[nb].block.infoTNW[0].blknum].block.info.sector;
         Octree.Roots[nb].block.infoTNW[0].level =  Octree.Roots[Octree.Roots[nb].block.infoTNW[0].blknum].block.info.level;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoTNW[0].dimen =  Octree.Roots[nb].block.infoTNW[0].dimen;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoTNW[0].sector =  Octree.Roots[nb].block.infoTNW[0].sector;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoTNW[0].level =  Octree.Roots[nb].block.infoTNW[0].level;
      }
      
      if(Octree.Roots[nb].block.infoTSW[0].blknum>=0){
         Octree.Roots[nb].block.infoTSW[0].dimen =  Octree.Roots[Octree.Roots[nb].block.infoTSW[0].blknum].block.info.dimen;
         Octree.Roots[nb].block.infoTSW[0].sector =  Octree.Roots[Octree.Roots[nb].block.infoTSW[0].blknum].block.info.sector;
         Octree.Roots[nb].block.infoTSW[0].level =  Octree.Roots[Octree.Roots[nb].block.infoTSW[0].blknum].block.info.level;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoTSW[0].dimen =  Octree.Roots[nb].block.infoTSW[0].dimen;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoTSW[0].sector =  Octree.Roots[nb].block.infoTSW[0].sector;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoTSW[0].level =  Octree.Roots[nb].block.infoTSW[0].level;
      }
      
      if(Octree.Roots[nb].block.infoTNE[0].blknum>=0){
         Octree.Roots[nb].block.infoTNE[0].dimen =  Octree.Roots[Octree.Roots[nb].block.infoTNE[0].blknum].block.info.dimen;
         Octree.Roots[nb].block.infoTNE[0].sector =  Octree.Roots[Octree.Roots[nb].block.infoTNE[0].blknum].block.info.sector;
         Octree.Roots[nb].block.infoTNE[0].level =  Octree.Roots[Octree.Roots[nb].block.infoTNE[0].blknum].block.info.level;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoTNE[0].dimen =  Octree.Roots[nb].block.infoTNE[0].dimen;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoTNE[0].sector =  Octree.Roots[nb].block.infoTNE[0].sector;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoTNE[0].level =  Octree.Roots[nb].block.infoTNE[0].level;
      }
      
      if(Octree.Roots[nb].block.infoTSE[0].blknum>=0){
         Octree.Roots[nb].block.infoTSE[0].dimen =  Octree.Roots[Octree.Roots[nb].block.infoTSE[0].blknum].block.info.dimen;
         Octree.Roots[nb].block.infoTSE[0].sector =  Octree.Roots[Octree.Roots[nb].block.infoTSE[0].blknum].block.info.sector;
         Octree.Roots[nb].block.infoTSE[0].level =  Octree.Roots[Octree.Roots[nb].block.infoTSE[0].blknum].block.info.level;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoTSE[0].dimen =  Octree.Roots[nb].block.infoTSE[0].dimen;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoTSE[0].sector =  Octree.Roots[nb].block.infoTSE[0].sector;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoTSE[0].level =  Octree.Roots[nb].block.infoTSE[0].level;
      }
      
      if(Octree.Roots[nb].block.infoBNW[0].blknum>=0){
         Octree.Roots[nb].block.infoBNW[0].dimen =  Octree.Roots[Octree.Roots[nb].block.infoBNW[0].blknum].block.info.dimen;
         Octree.Roots[nb].block.infoBNW[0].sector =  Octree.Roots[Octree.Roots[nb].block.infoBNW[0].blknum].block.info.sector;
         Octree.Roots[nb].block.infoBNW[0].level =  Octree.Roots[Octree.Roots[nb].block.infoBNW[0].blknum].block.info.level;  
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoBNW[0].dimen =  Octree.Roots[nb].block.infoBNW[0].dimen;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoBNW[0].sector =  Octree.Roots[nb].block.infoBNW[0].sector;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoBNW[0].level =  Octree.Roots[nb].block.infoBNW[0].level;
      }
      
      if(Octree.Roots[nb].block.infoBSW[0].blknum>=0){
         Octree.Roots[nb].block.infoBSW[0].dimen =  Octree.Roots[Octree.Roots[nb].block.infoBSW[0].blknum].block.info.dimen;
         Octree.Roots[nb].block.infoBSW[0].sector =  Octree.Roots[Octree.Roots[nb].block.infoBSW[0].blknum].block.info.sector;
         Octree.Roots[nb].block.infoBSW[0].level =  Octree.Roots[Octree.Roots[nb].block.infoBSW[0].blknum].block.info.level;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoBSW[0].dimen =  Octree.Roots[nb].block.infoBSW[0].dimen;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoBSW[0].sector =  Octree.Roots[nb].block.infoBSW[0].sector;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoBSW[0].level =  Octree.Roots[nb].block.infoBSW[0].level;
      }
      
      if(Octree.Roots[nb].block.infoBNE[0].blknum>=0){
         Octree.Roots[nb].block.infoBNE[0].dimen =  Octree.Roots[Octree.Roots[nb].block.infoBNE[0].blknum].block.info.dimen;
         Octree.Roots[nb].block.infoBNE[0].sector =  Octree.Roots[Octree.Roots[nb].block.infoBNE[0].blknum].block.info.sector;
         Octree.Roots[nb].block.infoBNE[0].level =  Octree.Roots[Octree.Roots[nb].block.infoBNE[0].blknum].block.info.level;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoBNE[0].dimen =  Octree.Roots[nb].block.infoBNE[0].dimen;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoBNE[0].sector =  Octree.Roots[nb].block.infoBNE[0].sector;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoBNE[0].level =  Octree.Roots[nb].block.infoBNE[0].level;
      }
      
      if(Octree.Roots[nb].block.infoBSE[0].blknum>=0){
         Octree.Roots[nb].block.infoBSE[0].dimen =  Octree.Roots[Octree.Roots[nb].block.infoBSE[0].blknum].block.info.dimen;
         Octree.Roots[nb].block.infoBSE[0].sector =  Octree.Roots[Octree.Roots[nb].block.infoBSE[0].blknum].block.info.sector;
         Octree.Roots[nb].block.infoBSE[0].level =  Octree.Roots[Octree.Roots[nb].block.infoBSE[0].blknum].block.info.level;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoBSE[0].dimen =  Octree.Roots[nb].block.infoBSE[0].dimen;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoBSE[0].sector =  Octree.Roots[nb].block.infoBSE[0].sector;
         Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoBSE[0].level =  Octree.Roots[nb].block.infoBSE[0].level;
         
      }
   } /* end of assignment of dimension, level and sector information  */
   
   // block orientation from roots block (global numbering system) to local block list (local block number system).
   for ( int nb = 0 ; nb <Initial_Mesh.NBlk ; ++nb ){
      
      if(blknumber_convt[nb].cpu ==  CFFC_MPI::This_Processor_Number){
         
         Local_Adaptive_Block_List.Block[ blknumber_convt[nb].blknum].infoT[0] = Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoT[0];  
         Local_Adaptive_Block_List.Block[ blknumber_convt[nb].blknum].infoB[0] = Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoB[0]; 
         Local_Adaptive_Block_List.Block[ blknumber_convt[nb].blknum].infoN[0] = Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoN[0]; 
         Local_Adaptive_Block_List.Block[ blknumber_convt[nb].blknum].infoS[0] = Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoS[0]; 
         Local_Adaptive_Block_List.Block[ blknumber_convt[nb].blknum].infoE[0] = Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoE[0]; 
         Local_Adaptive_Block_List.Block[ blknumber_convt[nb].blknum].infoW[0] = Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoW[0]; 
         
         for ( int iProc = 0; iProc !=  CFFC_MPI::Number_of_Processors; ++iProc ) {
            if (  CFFC_MPI::This_Processor_Number == iProc ) {
               cout<<"\n In AMR octree: "<<endl;
               cout<<"\n CFFC_MPI::This_Processor_Number = "<< CFFC_MPI::This_Processor_Number;
               cout<<" nb = "<<nb<<" blknumber_convt[nb].blknum "<< blknumber_convt[nb].blknum
                   <<"\n  infoE = "<< Local_Adaptive_Block_List.Block[ blknumber_convt[nb].blknum].infoE[0]
                   <<"\n  infoW = "<< Local_Adaptive_Block_List.Block[ blknumber_convt[nb].blknum].infoW[0];
               
               
               System::sleep(0.1);
            }
            MPI::COMM_WORLD.Barrier();
         }

         Local_Adaptive_Block_List.Block[ blknumber_convt[nb].blknum].infoNW[0] = Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoNW[0]; 
         Local_Adaptive_Block_List.Block[ blknumber_convt[nb].blknum].infoNE[0] = Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoNE[0]; 
         Local_Adaptive_Block_List.Block[ blknumber_convt[nb].blknum].infoSE[0] = Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoSE[0]; 
         Local_Adaptive_Block_List.Block[ blknumber_convt[nb].blknum].infoSW[0] = Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoSW[0]; 
         Local_Adaptive_Block_List.Block[ blknumber_convt[nb].blknum].infoTN[0] = Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoTN[0]; 
         Local_Adaptive_Block_List.Block[ blknumber_convt[nb].blknum].infoTS[0] = Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoTS[0]; 
         Local_Adaptive_Block_List.Block[ blknumber_convt[nb].blknum].infoTE[0] = Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoTE[0]; 
         Local_Adaptive_Block_List.Block[ blknumber_convt[nb].blknum].infoTW[0] = Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoTW[0]; 
         Local_Adaptive_Block_List.Block[ blknumber_convt[nb].blknum].infoBN[0] = Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoBN[0]; 
         Local_Adaptive_Block_List.Block[ blknumber_convt[nb].blknum].infoBS[0] = Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoBS[0]; 
         Local_Adaptive_Block_List.Block[ blknumber_convt[nb].blknum].infoBE[0] = Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoBE[0]; 
         Local_Adaptive_Block_List.Block[ blknumber_convt[nb].blknum].infoBW[0] = Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoBW[0]; 
         Local_Adaptive_Block_List.Block[ blknumber_convt[nb].blknum].infoTNW[0] = Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoTNW[0]; 
         Local_Adaptive_Block_List.Block[ blknumber_convt[nb].blknum].infoTSW[0] = Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoTSW[0]; 
         Local_Adaptive_Block_List.Block[ blknumber_convt[nb].blknum].infoTNE[0] = Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoTNE[0]; 
         Local_Adaptive_Block_List.Block[ blknumber_convt[nb].blknum].infoTSE[0] = Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoTSE[0]; 
         Local_Adaptive_Block_List.Block[ blknumber_convt[nb].blknum].infoBNW[0] = Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoBNW[0]; 
         Local_Adaptive_Block_List.Block[ blknumber_convt[nb].blknum].infoBSW[0] = Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoBSW[0]; 
         Local_Adaptive_Block_List.Block[ blknumber_convt[nb].blknum].infoBNE[0] = Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoBNE[0]; 
         Local_Adaptive_Block_List.Block[ blknumber_convt[nb].blknum].infoBSE[0] = Octree.Blocks[blknumber_convt[nb].cpu][blknumber_convt[nb].blknum]->block.infoBSE[0]; 
      }
   }
   
  
   /* Renumber all solution blocks, assigning a unique global block number. */

    Octree_DataStructure::Renumber_Solution_Blocks(Octree,
                                                   Local_Adaptive_Block_List);


    /* Modify block neighbours for grid geometries with periodic boundaries, etc. */

    Octree_DataStructure::Modify_Neighbours_of_Root_Solution_Blocks(Octree,
                                                                    Local_Adaptive_Block_List,
                                                                    Input.Grid_IP.i_Grid);

    /* Allocates memory for all message passing buffers used to send
       solution information between neighbouring solution blocks. */

    AdaptiveBlock3D_List::Allocate_Message_Buffers(Local_Adaptive_Block_List,
                                                   Local_Solution_Blocks.Soln_Blks[0].NumVar()+NUM_COMP_VECTOR3D);
 
    /* Solution block allocation and assignment complete.
       Return pointer to local solution blocks. */


    delete []blknumber_convt;
    

    return(0);

}

/********************************************************
 * Routine: Read_Octree                                 *
 *                                                      *
 * Reads the Octree data structure from a file.         *
 *                                                      *
 ********************************************************/
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
int Read_Octree(Octree_DataStructure                         &Octree,
                AdaptiveBlock3D_ResourceList                 &Global_Adaptive_Block_List,
                AdaptiveBlock3D_List                         &Local_Adaptive_Block_List,
                Input_Parameters<SOLN_pSTATE, SOLN_cSTATE>   &Input) {

    int i, nri, nrj, nrk, ncpu, nblk;
    int iBLK, jBLK, kBLK, iCPU;
    char Octree_file_name[256];
    char *Octree_file_name_ptr;
    ifstream Octree_file;    

    /* On primary processor, determine name of Octree input data file name. */

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
       strcat(Octree_file_name, "_Octree.tree");
       Octree_file_name_ptr = Octree_file_name;
    } /* endif */

    /* On primary processor, open the Octree data file. */

    if (CFFC_Primary_MPI_Processor()) {
       Octree_file.open(Octree_file_name_ptr, ios::in);
       if (Octree_file.fail()) return (1);
    } /* endif */

    /* On primary processor, read in the data structure size parameters and 
       re-allocate memory as required. */

    if (CFFC_Primary_MPI_Processor()) {
       Octree_file.setf(ios::skipws);
       Octree_file >> nri >> nrj >> nrk >> ncpu >> nblk;
       
       Octree_file.unsetf(ios::skipws);

       cout<<"\n in Read_Octree: Create_Octree_Data_Structure(Octree,";
       Create_Octree_Data_Structure(Octree,
                                    nri,
  	                            nrj,
  	                            nrk,
                                    Input.Number_of_Processors,
                                    Input.Number_of_Blocks_Per_Processor);
    } /* endif */

    /* Re-create and re-initialize the block resource list. */

    Create_Block_Resource_List(Global_Adaptive_Block_List,
                               Input.Number_of_Processors, 
	  		       Input.Number_of_Blocks_Per_Processor);

    /* On primary processor, read the Octree data from the file. */

    if (CFFC_Primary_MPI_Processor()) {
       for ( kBLK = 0 ; kBLK <= Octree.NRk-1 ; ++kBLK ) {
	 for ( jBLK = 0 ; jBLK <= Octree.NRj-1 ; ++jBLK ) {
           for ( iBLK = 0 ; iBLK <= Octree.NRi-1 ; ++iBLK ) {
	      Octree.Roots[iBLK][jBLK][kBLK].read(Octree_file,
                                                  Global_Adaptive_Block_List);
           } /* endfor */
           } /* endfor */
       } /* endfor */
    } /* endif */

    /* Broadcast Octree data structure to all processors. */
       
    Octree_DataStructure::Broadcast_Octree_Data_Structure(Octree,
                                                          Global_Adaptive_Block_List);

    /* Set the maximum and minimum refinement levels. */

    Octree.MaximumRefinementLevel = Input.Maximum_Refinement_Level-1;
    Octree.MinimumRefinementLevel = Input.Minimum_Refinement_Level-1;

    /* Set the thresholds for refinement and coarsening of the mesh. */

    Octree.RefineThreshold = Input.Threshold_for_Refinement;
    Octree.CoarsenThreshold = Input.Threshold_for_Coarsening;

    /* Re-evaluate Octree block pointers. */

    Octree.assign_block_pointers();

    /* (Re-)Allocate memory for local processor solution block list. */

    if (Local_Adaptive_Block_List.Nblk > 0) Local_Adaptive_Block_List.deallocate();
    Local_Adaptive_Block_List.allocate(Input.Number_of_Blocks_Per_Processor);
    Local_Adaptive_Block_List.ThisCPU = Global_Adaptive_Block_List.ThisCPU;

    /* Find the neighbours of the root blocks. */

     Octree_DataStructure::Find_Neighbours_of_Root_Solution_Blocks(Octree);

    /* Modify block neighbours for grid geometries with 
       periodic boundaries, etc. */

    Octree_DataStructure::Modify_Neighbours_of_Root_Solution_Blocks(Octree,
                                                                    Input.i_Grid);

    /* Determine the neighbouring blocks of all used (active)
       solution blocks in the Octree data structure. This will
       also copy block information to local processor solution block
       list. */

    Octree_DataStructure::Find_Neighbours(Octree,
                                          Local_Adaptive_Block_List);

    /* On primary processor, close Octree data file. */

    if (CFFC_Primary_MPI_Processor()) Octree_file.close();

    /* Reading of Octree data file complete.  Return zero value. */

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

    /* On primary processor, determine name of Octree output data file name. */

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
       strcat(Octree_file_name, "_Octree.tree");
       Octree_file_name_ptr = Octree_file_name;
    } /* endif */
    /* On primary processor, ppen the Octree data file. */

    if (CFFC_Primary_MPI_Processor()) {
       Octree_file.open(Octree_file_name_ptr, ios::out);
       if (Octree_file.fail()) return (1);
    } /* endif */

    /* On primary processor, write the Octree data to the file. */

    if (CFFC_Primary_MPI_Processor()) {
      Octree_file << Octree;
    }
    /* On primary processor, close Octree data file. */

    if (CFFC_Primary_MPI_Processor()) Octree_file.close();

    /* Writing of Octree data file complete.  Return zero value. */

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
void Flag_Blocks_For_Refinement(Hexa_Multi_Block<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> >  &Local_Solution_Blocks,
                                Octree_DataStructure                                     &Octree,
                                AdaptiveBlock3D_ResourceList                             &Global_Adaptive_Block_List,
                                AdaptiveBlock3D_List                                     &Local_Adaptive_Block_List) {

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
int Refine_Grid(Hexa_Multi_Block<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> >    &Local_Solution_Blocks,
                Octree_DataStructure                                       &Octree,
                AdaptiveBlock3D_ResourceList                               &Global_Adaptive_Block_List,
                AdaptiveBlock3D_List                                       &Local_Adaptive_Block_List) {

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
int Coarsen_Grid(Hexa_Multi_Block<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> >  &Local_Solution_Blocks,
                 Octree_DataStructure                                     &Octree,
                 AdaptiveBlock3D_ResourceList                             &Global_Adaptive_Block_List,
                 AdaptiveBlock3D_List                                     &Local_Adaptive_Block_List) {
  
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
void Fix_Refined_Block_Boundaries(Hexa_Multi_Block<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> >  &Local_Solution_Blocks,
                                  AdaptiveBlock3D_List                                     &Soln_Block_List) {

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
void Unfix_Refined_Block_Boundaries(Hexa_Multi_Block<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> >  &Local_Solution_Blocks,
                                    AdaptiveBlock3D_List                                     &Soln_Block_List) {

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
int AMR(Hexa_Multi_Block<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> >  &Local_Solution_Blocks,
        Octree_DataStructure                                     &Octree,
        AdaptiveBlock3D_ResourceList                             &Global_Adaptive_Block_List,
        AdaptiveBlock3D_List                                     &Local_Adaptive_Block_List,
        const int                                                Set_New_Refinement_Flags) {

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
int Initial_AMR(Hexa_Multi_Block<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> >  &Local_Solution_Blocks,
                Input_Parameters<SOLN_pSTATE, SOLN_cSTATE>               &Input,
                Octree_DataStructure                                     &Octree,
                AdaptiveBlock3D_ResourceList                             &Global_Adaptive_Block_List,
                AdaptiveBlock3D_List                                     &Local_Adaptive_Block_List) {

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
int Uniform_AMR(Hexa_Multi_Block<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> >  &Local_Solution_Blocks,
                Input_Parameters<SOLN_pSTATE, SOLN_cSTATE>               &Input,
                Octree_DataStructure                                     &Octree,
                AdaptiveBlock3D_ResourceList                             &Global_Adaptive_Block_List,
                AdaptiveBlock3D_List                                     &Local_Adaptive_Block_List) {

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
int Boundary_AMR(Hexa_Multi_Block<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> >  &Local_Solution_Blocks,
		 Input_Parameters<SOLN_pSTATE, SOLN_cSTATE>               &Input,
		 Octree_DataStructure                                     &Octree,
		 AdaptiveBlock3D_ResourceList                             &Global_Adaptive_Block_List,
		 AdaptiveBlock3D_List                                     &Local_Adaptive_Block_List) {

    cout<<"\nError Boundary_AMR is not written for 3D";
    return(0);

}

#endif // _AMR_INCLUDED

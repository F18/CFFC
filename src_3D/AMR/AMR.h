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
/*  template<typename SOLN_pSTATE, typename SOLN_cSTATE> */
/* int Create_Initial_Solution_Blocks(Grid3D_Hexa_Multi_Block_List                                 &Initial_Mesh, */
/*                                    Hexa_Multi_Block<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> > &Local_Solution_Blocks, */
/*                                    Input_Parameters<SOLN_pSTATE, SOLN_cSTATE>              &Input, */
/*                                    Octree_DataStructure                                    &Octree, */
/*                                    AdaptiveBlock3D_ResourceList                            &Global_Adaptive_Block_List, */
/*                                    AdaptiveBlock3D_List                                    &Local_Adaptive_Block_List) { */
   

  
/*    int totalblks, n_cpu, n_blk; */
   
   
/*    Octree_DataStructure::Create_Octree_Data_Structure( */
/*       Octree, */
/*       Input.IP_Grid.IBlk, */
/*       Input.IP_Grid.JBlk, */
/*       Input.IP_Grid.KBlk, */
/*       Input.Number_of_Processors, */
/*       Input.Number_of_Blocks_Per_Processor); */

   
/*    Octree.MaximumRefinementLevel = Input.Maximum_Refinement_Level-1; */
/*    Octree.MinimumRefinementLevel = Input.Minimum_Refinement_Level-1; */
   
/*   /\*  cout<<"\n/\\* Set the thresholds for refinement and coarsening of the mesh. *\\/"; *\/ */
/* /\*    cout.flush(); *\/ */
   
/*    Octree.RefineThreshold = Input.Threshold_for_Refinement; */
/*    Octree.CoarsenThreshold = Input.Threshold_for_Coarsening; */
   
/*   /\*  cout<<"\n/\\* Create (allocate) array of local 3D hexarilateral solution blocks. *\\/"; *\/ */
/* /\*    cout.flush(); *\/ */
   
/*    AdaptiveBlock3D_ResourceList::Create_Block_Resource_List(GlobalSolnBlockList, */
/*                               Input.Number_of_Processors, */
/*                               Input.Number_of_Blocks_Per_Processor); */
/*    LocalSolnBlockList.allocate(Input.Number_of_Blocks_Per_Processor); */
/*    LocalSolnBlockList.ThisCPU = GlobalSolnBlockList.ThisCPU; */
   
/*  /\*   cout<<"\n Loop over all initial mesh blocks and assign Octree root blocks and local solution block information as required. "; *\/ */
/* /\*    cout.flush(); *\/ */

/*    int nused=0; */
/*    totalblks = Input.IP_Grid.IBlk*Input.IP_Grid.JBlk*Input.IP_Grid.KBlk; */
   
/*    for ( int nb = 0 ; nb < totalblks ; ++nb ) { */
/*       if (InitMeshBlks.Grid_Blks[nb]->Node != NULL) { // Mesh block is used!!!! */
/*          // Get next free solution block from list of available */
/*          // solution blocks. */

/*          if (GlobalSolnBlockList.Nfree > 0) { */
/*             n_cpu = GlobalSolnBlockList.nextCPU(); */
/*             n_blk = GlobalSolnBlockList.nextBlock(); */
            
/*             GlobalSolnBlockList.update_next(); */
            
/*          } else { */
/*             cout << "\n"  */
/*                  << " AMR Error: Create_Initial_Solution_Blocks, Insufficient number of hexahedrial solution blocks.\n"; */
/* //                  Soln_ptr ->deallocate( ); */
/*             //  return (NULL); */
/*             return (1); */
/*          } /\* endif *\/ */
         
/*          // Assign block information to appropriate Octree root solution block. */
/*          Octree.Rootblks[nb].block.used = ADAPTIVEBLOCK3D_USED; */
/*          Octree.Rootblks[nb].block.gblknum = */
/*             GlobalSolnBlockList.Nused-1; */
/*          Octree.Rootblks[nb].block.info.cpu = n_cpu; */
/*          Octree.Rootblks[nb].block.info.blknum = n_blk; */
/*          Octree.Rootblks[nb].block.info.dimen.i = */
/*             InitMeshBlks.Grid_Blks[nb]->NCi-  2*InitMeshBlks.Grid_Blks[nb]->Nghost; */
/*          Octree.Rootblks[nb].block.info.dimen.j = */
/*             InitMeshBlks.Grid_Blks[nb]->NCj- 2*InitMeshBlks.Grid_Blks[nb]->Nghost; */
/*          Octree.Rootblks[nb].block.info.dimen.k = */
/*             InitMeshBlks.Grid_Blks[nb]->NCk- 2*InitMeshBlks.Grid_Blks[nb]->Nghost; */
/*          Octree.Rootblks[nb].block.info.dimen.ghost = 2; */
/*          Octree.Rootblks[nb].block.info.sector = ADAPTIVEBLOCK3D_SECTOR_NONE; */
/*          Octree.Rootblks[nb].block.info.level = 0; */
/*          Octree.Rootblks[nb].parent_ptr = NULL; */
/*          Octree.Rootblks[nb].childTNW_ptr = NULL; */
/*          Octree.Rootblks[nb].childTNE_ptr = NULL; */
/*          Octree.Rootblks[nb].childTSE_ptr = NULL; */
/*          Octree.Rootblks[nb].childTSW_ptr = NULL; */
/*          Octree.Rootblks[nb].childBNW_ptr = NULL; */
/*          Octree.Rootblks[nb].childBNE_ptr = NULL; */
/*          Octree.Rootblks[nb].childBSE_ptr = NULL; */
/*          Octree.Rootblks[nb].childBSW_ptr = NULL; */
/*          Octree.Blocks[n_cpu][n_blk] = &(Octree.Rootblks[nb]); */
               
/*          // For solution blocks on this processor (or processor */
/*          // element), add block to local list, create the solution */
/*          // block, and copy the appropriate block of the */
/*          // initial hexarilateral mesh to the solution block mesh. */
/*          if (GlobalSolnBlockList.ThisCPU == n_cpu) { */
            
/*             LocalSolnBlockList.Block[n_blk] = Octree.Rootblks[nb].block; */
            
/*             Hexa_Multi_Block.Hexa_Block_List[n_blk] =  */
/*                new Hexa_Block<SOLN_pSTATE, SOLN_cSTATE>(Input.IP_Grid.ICells, */
/*                                                         Input.IP_Grid.JCells, */
/*                                                         Input.IP_Grid.KCells, */
/*                                                         Input.IP_Grid.Nghost, */
/*                                                         nb, */
/*                                                         Input.i_Flow_Type,  */
/*                                                         InitMeshBlks.Grid_Blks[nb]); */
            
/*             Hexa_Multi_Block.Block_Used[n_blk] = LocalSolnBlockList.Block[n_blk].used; */
/*             //      cout<<"\n AMR3D n_blk "<<n_blk<<"  used=   "<< LocalSolnBlockList.Block[n_blk].used<<endl; */
            
/*          } /\* endif *\/ */
         
/*       } /\* endif *\/ */
/*       else{ */
/*          Octree.Rootblks[nb].block.used = 0; */
         
/*          nused++; */
/*          //  cout<<"\n Block ("<<i_blk <<","<<j_blk <<","<<k_blk <<") is not used"; */
/*       } */
/*    } /\* endfor *\/ */

/*    /\*   cout<<"\n/\\* Renumber all solution blocks, assigning a unique global block number. *\\/"; *\/ */
/* /\*     cout.flush(); *\/ */
   
/*    Octree_DataStructure::Renumber_Solution_Blocks(Octree, */
/*                                                        LocalSolnBlockList); */

/*    /\*  cout<<"\n/\\* Find the neighbours of all of the newly assigned root blocks. *\\/"; *\/ */
/* /\*     cout.flush(); *\/ */
   
   
/*  /\*   Octree_DataStructure::Find_Neighbours_of_Root_Solution_Blocks(Octree, *\/ */
/* /\*                                                                       LocalSolnBlockList); *\/ */
  
/* /\*    //cout<<"\n/\\* Modify block neighbours for grid geometries with *\/ */
/* /\*    //   periodic boundaries, etc. *\\/"; *\/ */
/* /\*    //cout.flush(); *\/ */

/* /\*    Octree_DataStructure::Modify_Neighbours_of_Root_Solution_Blocks(Octree, *\/ */
/* /\*                                                                         LocalSolnBlockList, *\/ */
/* /\*                                                                         Input.IP_Grid.i_Grid); *\/ */

   
/*    //cout<<"\n/\* Allocates memory for all message passing buffers used */
/*    //   to send solution information between neighbouring solution blocks. *\/"; */
/*    //cout.flush(); */
   
/*    AdaptiveBlock3D_List::Allocate_Message_Buffers( */
/*       LocalSolnBlockList, */
/*       Hexa_Multi_Block.Hexa_Block_List[0]->NumVar()+NUM_COMP_VECTOR3D); */
   
   
/*    //cout<<"\n/\* Solution block allocation and assignment complete. */
/*    //   Return pointer to local solution blocks. *\/"; */
/*    //cout.flush(); */
   
/*     return(0); */

/* } */

template<typename SOLN_pSTATE, typename SOLN_cSTATE>
int Create_Initial_Solution_Blocks(Grid3D_Hexa_Multi_Block_List                            &Initial_Mesh,
                                   Hexa_Multi_Block<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> > &Local_Solution_Blocks,
                                   Input_Parameters<SOLN_pSTATE, SOLN_cSTATE>              &Input,
                                   Octree_DataStructure                                    &Octree,
                                   AdaptiveBlock3D_ResourceList                            &Global_Adaptive_Block_List,
                                   AdaptiveBlock3D_List                                    &Local_Adaptive_Block_List) {
   
   int n_cpu, n_blk;

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

   int nused=0;
   int totalblks = Initial_Mesh.NBlk_Kdir*Initial_Mesh.NBlk_Jdir*Initial_Mesh.NBlk_Idir;
   
   for ( int nb = 0 ; nb < totalblks ; ++nb ){
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
         Octree.Rootblks[nb].block.used = ADAPTIVEBLOCK3D_USED;
         Octree.Rootblks[nb].block.gblknum =
            Global_Adaptive_Block_List.Nused-1;
         Octree.Rootblks[nb].block.info.cpu = n_cpu;
         Octree.Rootblks[nb].block.info.blknum = n_blk;
         Octree.Rootblks[nb].block.info.dimen.i =
            Initial_Mesh.Grid_Blks[nb].NCi -
            2*Initial_Mesh.Grid_Blks[nb].Nghost;
         Octree.Rootblks[nb].block.info.dimen.j =
            Initial_Mesh.Grid_Blks[nb].NCj -
            2*Initial_Mesh.Grid_Blks[nb].Nghost;
         Octree.Rootblks[nb].block.info.dimen.k =
            Initial_Mesh.Grid_Blks[nb].NCk -
            2*Initial_Mesh.Grid_Blks[nb].Nghost;
         Octree.Rootblks[nb].block.info.dimen.ghost = 2;
         Octree.Rootblks[nb].block.info.sector = ADAPTIVEBLOCK3D_SECTOR_NONE;
         Octree.Rootblks[nb].block.info.level = 0;
         Octree.Rootblks[nb].parent_ptr = NULL;
         Octree.Rootblks[nb].childTNW_ptr = NULL;
         Octree.Rootblks[nb].childTNE_ptr = NULL;
         Octree.Rootblks[nb].childTSE_ptr = NULL;
         Octree.Rootblks[nb].childTSW_ptr = NULL;
         Octree.Rootblks[nb].childBNW_ptr = NULL;
         Octree.Rootblks[nb].childBNE_ptr = NULL;
         Octree.Rootblks[nb].childBSE_ptr = NULL;
         Octree.Rootblks[nb].childBSW_ptr = NULL;
         Octree.Blocks[n_cpu][n_blk] = &(Octree.Rootblks[nb]);
         
         // For solution blocks on this processor (or processor
         // element), add block to local list, create the solution
         // block, and copy the appropriate block of the
         // initial hexarilateral mesh to the solution block mesh.
         if (Global_Adaptive_Block_List.ThisCPU == n_cpu) {
            Local_Adaptive_Block_List.Block[n_blk] = Octree.Rootblks[nb].block;
            Local_Solution_Blocks.Soln_Blks[n_blk].Create_Block(Initial_Mesh.Grid_Blks[nb]);
            Local_Solution_Blocks.Soln_Blks[n_blk].Flow_Type = Input.i_Flow_Type;
            Local_Solution_Blocks.Block_Used[n_blk] = HEXA_BLOCK_USED;
         } /* endif */
      } else{
         Octree.Rootblks[nb].block.used = 0;
         nused++;
      } /* endif */
      
   } /* endfor */
   
   /* Renumber all solution blocks, assigning a unique global block number. */

    Octree_DataStructure::Renumber_Solution_Blocks(Octree,
                                                   Local_Adaptive_Block_List);

    /* Find the neighbours of all of the newly assigned root blocks. */
 
    Octree_DataStructure::Find_Neighbours_of_Root_Solution_Blocks(Octree,
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

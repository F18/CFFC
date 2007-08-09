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
   int  Create_Initial_Solution_Blocks(Grid3D_Hexa_Multi_Block &Initial_Multiblock_Mesh,
                                       Hexa_MultiBlock<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> > &Hexa_MultiBlock_List,
                                       Input_Parameters<SOLN_pSTATE, SOLN_cSTATE>   &InputParameters,
                                       OcTreeBlock_DataStructure &OcTree,
                                       AdaptiveBlock3DResourceList  &GlobalSolnBlockList,
                                       AdaptiveBlock3D_List &LocalSolnBlockList) {
   
   int i_blk, j_blk, k_blk, n_cpu, n_blk;
   
   OcTreeBlock_DataStructure::Create_OcTree_Data_Structure(
      OcTree,
      InputParameters.IP_Grid.NBlk_Idir,
      InputParameters.IP_Grid.NBlk_Jdir,
      InputParameters.IP_Grid.NBlk_Kdir,
      InputParameters.Number_of_Processors,
      InputParameters.Number_of_Blocks_Per_Processor);

   OcTree.MaximumRefinementLevel = InputParameters.Maximum_Refinement_Level-1;
   OcTree.MinimumRefinementLevel = InputParameters.Minimum_Refinement_Level-1;
   
   /* Set the thresholds for refinement and coarsening of the mesh. */
   
   OcTree.RefineThreshold = InputParameters.Threshold_for_Refinement;
   OcTree.CoarsenThreshold = InputParameters.Threshold_for_Coarsening;
   
   /* Create (allocate) array of local 3D hexarilateral solution blocks. */
   
   AdaptiveBlock3DResourceList::Create_Block_Resource_List(GlobalSolnBlockList,
                              InputParameters.Number_of_Processors,
                              InputParameters.Number_of_Blocks_Per_Processor);
   LocalSolnBlockList.allocate(InputParameters.Number_of_Blocks_Per_Processor);
   LocalSolnBlockList.ThisCPU = GlobalSolnBlockList.ThisCPU;
   
   /* Loop over all initial mesh blocks and assign OcTree root blocks and 
      local solution block information as required. */

   int nused=0;
   for ( k_blk = 0 ; k_blk <= InputParameters.IP_Grid.NBlk_Kdir-1 ; ++k_blk ) {
      for ( j_blk = 0 ; j_blk <= InputParameters.IP_Grid.NBlk_Jdir-1 ; ++j_blk ) {
         for ( i_blk = 0 ; i_blk <= InputParameters.IP_Grid.NBlk_Idir-1 ; ++i_blk ) {
            if (Initial_Multiblock_Mesh.Grid_Blks[i_blk][j_blk][k_blk].Used) { // Mesh block is used!!!!
	       // Get next free solution block from list of available
	       // solution blocks.
               
               if (GlobalSolnBlockList.Nfree > 0) {
                  n_cpu = GlobalSolnBlockList.nextCPU();
                  n_blk = GlobalSolnBlockList.nextBlock();
                  
                  GlobalSolnBlockList.update_next();
                     
               } else {
                  cout << "\n" 
                       << " AMR Error: Create_Initial_Solution_Blocks, Insufficient number of hexahedrial solution blocks.\n";
//                  Soln_ptr ->deallocate( );
                  //  return (NULL);
                  return (1);
               } /* endif */
	        
               // Assign block information to appropriate OcTree root solution block.
               OcTree.Roots[i_blk][j_blk][k_blk].block.used = ADAPTIVEBLOCK3D_USED;
               OcTree.Roots[i_blk][j_blk][k_blk].block.gblknum =
	          GlobalSolnBlockList.Nused-1;
               OcTree.Roots[i_blk][j_blk][k_blk].block.info.cpu = n_cpu;
               OcTree.Roots[i_blk][j_blk][k_blk].block.info.blknum = n_blk;
               OcTree.Roots[i_blk][j_blk][k_blk].block.info.dimen.i =
                   Initial_Multiblock_Mesh.Grid_Blks[i_blk][j_blk][k_blk].NCi -  
                   2*Initial_Multiblock_Mesh.Grid_Blks[i_blk][j_blk][k_blk].Nghost;
               OcTree.Roots[i_blk][j_blk][k_blk].block.info.dimen.j =
                   Initial_Multiblock_Mesh.Grid_Blks[i_blk][j_blk][k_blk].NCj - 
                   2*Initial_Multiblock_Mesh.Grid_Blks[i_blk][j_blk][k_blk].Nghost;
               OcTree.Roots[i_blk][j_blk][k_blk].block.info.dimen.k =
                   Initial_Multiblock_Mesh.Grid_Blks[i_blk][j_blk][k_blk].NCk - 
                   2*Initial_Multiblock_Mesh.Grid_Blks[i_blk][j_blk][k_blk].Nghost;
               OcTree.Roots[i_blk][j_blk][k_blk].block.info.dimen.ghost = 2;
               OcTree.Roots[i_blk][j_blk][k_blk].block.info.sector = ADAPTIVEBLOCK3D_SECTOR_NONE;
               OcTree.Roots[i_blk][j_blk][k_blk].block.info.level = 0;
               OcTree.Roots[i_blk][j_blk][k_blk].parent_ptr = NULL;
               OcTree.Roots[i_blk][j_blk][k_blk].childTNW_ptr = NULL;
               OcTree.Roots[i_blk][j_blk][k_blk].childTNE_ptr = NULL;
               OcTree.Roots[i_blk][j_blk][k_blk].childTSE_ptr = NULL;
               OcTree.Roots[i_blk][j_blk][k_blk].childTSW_ptr = NULL;
               OcTree.Roots[i_blk][j_blk][k_blk].childBNW_ptr = NULL;
               OcTree.Roots[i_blk][j_blk][k_blk].childBNE_ptr = NULL;
               OcTree.Roots[i_blk][j_blk][k_blk].childBSE_ptr = NULL;
               OcTree.Roots[i_blk][j_blk][k_blk].childBSW_ptr = NULL;
               OcTree.Blocks[n_cpu][n_blk] = &(OcTree.Roots[i_blk][j_blk][k_blk]);
               
               // For solution blocks on this processor (or processor
               // element), add block to local list, create the solution
               // block, and copy the appropriate block of the
               // initial hexarilateral mesh to the solution block mesh.
               if (GlobalSolnBlockList.ThisCPU == n_cpu) {
                  
                  LocalSolnBlockList.Block[n_blk] = OcTree.Roots[i_blk][j_blk][k_blk].block;
                  
                  Hexa_MultiBlock_List.Hexa_Block_List[n_blk] = 
                     new Hexa_Block<SOLN_pSTATE, SOLN_cSTATE>(InputParameters.IP_Grid.NCells_Idir,
                                                              InputParameters.IP_Grid.NCells_Jdir,
                                                              InputParameters.IP_Grid.NCells_Kdir,
                                                              InputParameters.IP_Grid.Nghost,
                                                              i_blk, j_blk, k_blk,
                                                              InputParameters.i_Flow_Type, 
                                                              Initial_Multiblock_Mesh.Grid_Blks[i_blk][j_blk][k_blk]);

                  Hexa_MultiBlock_List.Block_Used[n_blk] = LocalSolnBlockList.Block[n_blk].used;
                  
               } /* endif */
	       
           } /* endif */
	   else{
	     OcTree.Roots[i_blk][j_blk][k_blk].block.used = 0;
	     nused++;
	   }
        } /* endfor */
      } /* endfor */
    } /* endfor */
   

    /* Renumber all solution blocks, assigning a unique global block number. */

    OcTreeBlock_DataStructure::Renumber_Solution_Blocks(OcTree,
                             LocalSolnBlockList);

    /* Find the neighbours of all of the newly assigned root blocks. */
 
    OcTreeBlock_DataStructure::Find_Neighbours_of_Root_Solution_Blocks(OcTree,
                                            LocalSolnBlockList);

    /* Modify block neighbours for grid geometries with periodic boundaries, etc. */

    OcTreeBlock_DataStructure::Modify_Neighbours_of_Root_Solution_Blocks(OcTree,
                                                                         LocalSolnBlockList,
                                                                         InputParameters.IP_Grid.i_Grid);

    /* Allocates memory for all message passing buffers used to send 
       solution information between neighbouring solution blocks. */

    AdaptiveBlock3D_List::Allocate_Message_Buffers(LocalSolnBlockList,
                                                   Hexa_MultiBlock_List.Hexa_Block_List[0]->NumVar()+NUM_COMP_VECTOR3D);
 
    /* Solution block allocation and assignment complete.
       Return pointer to local solution blocks. */

    return(0);

}

/********************************************************
 * Routine: Read_OcTree                                 *
 *                                                      *
 * Reads the OcTree data structure from a file.         *
 *                                                      *
 ********************************************************/
template <class Hexa_Soln_Input_Parameters>
int Read_OcTree(OcTreeBlock_DataStructure &OcTree,
                  AdaptiveBlock3DResourceList   &List_of_Available_Blocks,
                  AdaptiveBlock3D_List        &Local_Soln_Block_List,
                  Hexa_Soln_Input_Parameters  &Input_Parameters) {

    int i, nri, nrj, nrk, ncpu, nblk;
    int iBLK, jBLK, kBLK, iCPU;
    char OcTree_file_name[256];
    char *OcTree_file_name_ptr;
    ifstream OcTree_file;    

    /* On primary processor, determine name of OcTree input data file name. */

    if (CFFC_Primary_MPI_Processor()) {
       i = 0;
       while (1) {
          if (Input_Parameters.Restart_File_Name[i] == ' ' ||
              Input_Parameters.Restart_File_Name[i] == '.') break;
          OcTree_file_name[i]=Input_Parameters.Restart_File_Name[i];
          i = i + 1;
          if (i > strlen(Input_Parameters.Restart_File_Name) ) break;
       } /* endwhile */
       OcTree_file_name[i] = '\0';
       strcat(OcTree_file_name, "_OcTree.tree");
       OcTree_file_name_ptr = OcTree_file_name;
    } /* endif */

    /* On primary processor, open the OcTree data file. */

    if (CFFC_Primary_MPI_Processor()) {
       OcTree_file.open(OcTree_file_name_ptr, ios::in);
       if (OcTree_file.bad()) return (1);
    } /* endif */

    /* On primary processor, read in the data structure size parameters and 
       re-allocate memory as required. */

    if (CFFC_Primary_MPI_Processor()) {
       OcTree_file.setf(ios::skipws);
       OcTree_file >> nri >> nrj >> nrk >> ncpu >> nblk;
       
       OcTree_file.unsetf(ios::skipws);

       cout<<"\n in Read_Octree: Create_OcTree_Data_Structure(OcTree,";
       Create_OcTree_Data_Structure(OcTree,
                                      nri,
  	                              nrj,
  	                              nrk,
                                      Input_Parameters.Number_of_Processors,
                                      Input_Parameters.Number_of_Blocks_Per_Processor);
    } /* endif */

    /* Re-create and re-initialize the block resource list. */

    Create_Block_Resource_List(List_of_Available_Blocks,
                               Input_Parameters.Number_of_Processors, 
	  		       Input_Parameters.Number_of_Blocks_Per_Processor);

    /* On primary processor, read the OcTree data from the file. */

    if (CFFC_Primary_MPI_Processor()) {
       for ( kBLK = 0 ; kBLK <= OcTree.NRk-1 ; ++kBLK ) {
	 for ( jBLK = 0 ; jBLK <= OcTree.NRj-1 ; ++jBLK ) {
           for ( iBLK = 0 ; iBLK <= OcTree.NRi-1 ; ++iBLK ) {
	      OcTree.Roots[iBLK][jBLK][kBLK].read(OcTree_file,
                                              List_of_Available_Blocks);
           } /* endfor */
           } /* endfor */
       } /* endfor */
    } /* endif */

    /* Broadcast OcTree data structure to all processors. */
       
    OcTreeBlock_DataStructure::Broadcast_OcTree_Data_Structure(OcTree,
                                      List_of_Available_Blocks);

    /* Set the maximum and minimum refinement levels. */

    OcTree.MaximumRefinementLevel = Input_Parameters.Maximum_Refinement_Level-1;
    OcTree.MinimumRefinementLevel = Input_Parameters.Minimum_Refinement_Level-1;

    /* Set the thresholds for refinement and coarsening of the mesh. */

    OcTree.RefineThreshold = Input_Parameters.Threshold_for_Refinement;
    OcTree.CoarsenThreshold = Input_Parameters.Threshold_for_Coarsening;

    /* Re-evaluate OcTree block pointers. */

    OcTree.assign_block_pointers();

    /* (Re-)Allocate memory for local processor solution block list. */

    if (Local_Soln_Block_List.Nblk > 0) Local_Soln_Block_List.deallocate();
    Local_Soln_Block_List.allocate(Input_Parameters.Number_of_Blocks_Per_Processor);
    Local_Soln_Block_List.ThisCPU = List_of_Available_Blocks.ThisCPU;

    /* Find the neighbours of the root blocks. */

     OcTreeBlock_DataStructure::Find_Neighbours_of_Root_Solution_Blocks(OcTree);

    /* Modify block neighbours for grid geometries with 
       periodic boundaries, etc. */

     OcTreeBlock_DataStructure::Modify_Neighbours_of_Root_Solution_Blocks(OcTree,
                                              Input_Parameters.i_Grid);

    /* Determine the neighbouring blocks of all used (active)
       solution blocks in the OcTree data structure. This will
       also copy block information to local processor solution block
       list. */

    OcTreeBlock_DataStructure::Find_Neighbours(OcTree,
                    Local_Soln_Block_List);

    /* On primary processor, close OcTree data file. */

    if (CFFC_Primary_MPI_Processor()) OcTree_file.close();

    /* Reading of OcTree data file complete.  Return zero value. */

    return(0);

}

/********************************************************
 * Routine: Write_OcTree                                *
 *                                                      *
 * Writes the OcTree data structure to a file for       *
 * later retrieval.                                     *
 *                                                      *
 ********************************************************/
template <class Hexa_Soln_Input_Parameters>
int Write_OcTree(OcTreeBlock_DataStructure &OcTree,
                   Hexa_Soln_Input_Parameters  &Input_Parameters) {

    int i;
    char OcTree_file_name[256];
    char *OcTree_file_name_ptr;
    ofstream OcTree_file;    

    /* On primary processor, determine name of OcTree output data file name. */

    if (CFFC_Primary_MPI_Processor()) {
       i = 0;
       while (1) {
          if (Input_Parameters.Restart_File_Name[i] == ' ' ||
              Input_Parameters.Restart_File_Name[i] == '.') break;
          OcTree_file_name[i]=Input_Parameters.Restart_File_Name[i];
          i = i + 1;
          if (i > strlen(Input_Parameters.Restart_File_Name) ) break;
       } /* endwhile */
       OcTree_file_name[i] = '\0';
       strcat(OcTree_file_name, "_OcTree.tree");
       OcTree_file_name_ptr = OcTree_file_name;
    } /* endif */
    /* On primary processor, ppen the OcTree data file. */

    if (CFFC_Primary_MPI_Processor()) {
       OcTree_file.open(OcTree_file_name_ptr, ios::out);
       if (OcTree_file.bad()) return (1);
    } /* endif */

    /* On primary processor, write the OcTree data to the file. */

    if (CFFC_Primary_MPI_Processor()) {
      OcTree_file << OcTree;
    }
    /* On primary processor, close OcTree data file. */

    if (CFFC_Primary_MPI_Processor()) OcTree_file.close();

    /* Writing of OcTree data file complete.  Return zero value. */

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
template <class Hexa_Soln_Block>
void Flag_Blocks_For_Refinement(Hexa_Soln_Block             *Soln_ptr,
                                OcTreeBlock_DataStructure &OcTree,
                                AdaptiveBlock3DResourceList   &Global_Soln_Block_List,
                                AdaptiveBlock3D_List        &Local_Soln_Block_List) {

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
template <class Hexa_Soln_Block>
int Refine_Grid(Hexa_Soln_Block             *Soln_ptr,
                OcTreeBlock_DataStructure &OcTree,
                AdaptiveBlock3DResourceList   &GlobalSolnBlockList,
                AdaptiveBlock3D_List        &LocalSolnBlockList) {

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
template <class Hexa_Soln_Block>
int Coarsen_Grid(Hexa_Soln_Block             *Soln_ptr,
                 OcTreeBlock_DataStructure &OcTree,
                 AdaptiveBlock3DResourceList   &GlobalSolnBlockList,
                 AdaptiveBlock3D_List        &LocalSolnBlockList) {
  
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
template <class Hexa_Soln_Block>
void Fix_Refined_Block_Boundaries(Hexa_Soln_Block      *Soln_ptr,
                                  AdaptiveBlock3D_List &Soln_Block_List) {

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
template <class Hexa_Soln_Block>
void Unfix_Refined_Block_Boundaries(Hexa_Soln_Block      *Soln_ptr,
                                    AdaptiveBlock3D_List &Soln_Block_List) {

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
template <class Hexa_Soln_Block>
int AMR(Hexa_Soln_Block             *Soln_ptr,
        OcTreeBlock_DataStructure &OcTree,
        AdaptiveBlock3DResourceList   &GlobalSolnBlockList,
        AdaptiveBlock3D_List        &LocalSolnBlockList,
        const int                    Set_New_Refinement_Flags) {

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
template <class Hexa_Soln_Block, class Hexa_Soln_Input_Parameters>
int Initial_AMR(Hexa_Soln_Block              *Soln_ptr,
                Hexa_Soln_Input_Parameters   &InputParameters,
                OcTreeBlock_DataStructure  &OcTree,
                AdaptiveBlock3DResourceList    &GlobalSolnBlockList,
                AdaptiveBlock3D_List         &LocalSolnBlockList) {

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
template <class Hexa_Soln_Block, class Hexa_Soln_Input_Parameters>
int Uniform_AMR(Hexa_Soln_Block              *Soln_ptr,
                Hexa_Soln_Input_Parameters   &InputParameters,
                OcTreeBlock_DataStructure  &OcTree,
                AdaptiveBlock3DResourceList    &GlobalSolnBlockList,
                AdaptiveBlock3D_List         &LocalSolnBlockList) {

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
template <class Hexa_Soln_Block, class Hexa_Soln_Input_Parameters>
int Boundary_AMR(Hexa_Soln_Block              *Soln_ptr,
		 Hexa_Soln_Input_Parameters   &InputParameters,
		 OcTreeBlock_DataStructure  &OcTree,
		 AdaptiveBlock3DResourceList    &GlobalSolnBlockList,
		 AdaptiveBlock3D_List         &LocalSolnBlockList) {

    cout<<"\nError Boundary_AMR is not written for 3D";
    return(0);

}

#endif /* _AMR_INCLUDED  */

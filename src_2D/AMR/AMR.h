/* AMR.h:  Header file for templated subroutines for carrying 
           out the adaptive mesh refinement (AMR). */

#ifndef _AMR_INCLUDED
#define _AMR_INCLUDED

/* Include 2D quadrilateral multiblock grid, adaptive block,
   and quadtree header files. */

#ifndef _GRID2D_QUAD_BLOCK_INCLUDED
#include "../Grid/Grid2DQuad.h"
#endif // _GRID2D_QUAD_BLOCK_INCLUDED

#ifndef _ADAPTIVBLOCK_INCLUDED
#include "AdaptiveBlock.h"
#endif // _ADAPTIVEBLOCK_INCLUDED

#ifndef _QUADTREE_INCLUDED
#include "QuadTree.h"
#endif // _QUADTREE_INCLUDED

#include "../HighOrderReconstruction/CENO_ExecutionMode.h"

/******************************************************************
 * AMR -- Templated subroutines.                                  *
 ******************************************************************/

/******************************************************************
 * Routine: Create_Initial_Solution_Blocks                        *
 *                                                                *
 * Assigns and creates (allocates) initial 2D quadrilateral       *
 * solution blocks corresponding to the initial grid.  This       *
 * routine also creates the quadtree, local block, and global     *
 * block resource list data structures for performing the         *
 * block-base adaptive mesh refinement (AMR) and determines the   *
 * roots of the quadtree data structure.                          *
 *                                                                *
 ******************************************************************/
template <class Quad_Soln_Block, class Quad_Soln_Input_Parameters, class Quad_Grid_Block>
Quad_Soln_Block* Create_Initial_Solution_Blocks(Quad_Grid_Block            **InitMeshBlk,
                                                Quad_Soln_Block              *Soln_ptr,
                                                Quad_Soln_Input_Parameters   &InputParameters,
                                                QuadTreeBlock_DataStructure  &QuadTree,
                                                AdaptiveBlockResourceList    &GlobalSolnBlockList,
                                                AdaptiveBlock2D_List         &LocalSolnBlockList) {

    int i_blk, j_blk, n_cpu, n_blk;

    /* Create (allocate) multi-block quadtree data structure. */

    Create_QuadTree_Data_Structure(QuadTree,
                                   InputParameters.Number_of_Blocks_Idir,
  	                           InputParameters.Number_of_Blocks_Jdir,
                                   InputParameters.Number_of_Processors,
                                   InputParameters.Number_of_Blocks_Per_Processor);

    /* Set the maximum and minimum refinement levels. */

    QuadTree.MaximumRefinementLevel = InputParameters.Maximum_Refinement_Level-1;
    QuadTree.MinimumRefinementLevel = InputParameters.Minimum_Refinement_Level-1;

    /* Set the thresholds for refinement and coarsening of the mesh. */

    QuadTree.RefineThreshold = InputParameters.Threshold_for_Refinement;
    QuadTree.CoarsenThreshold = InputParameters.Threshold_for_Coarsening;

    /* Create (allocate) array of local 2D quadrilateral solution blocks. */

    Create_Block_Resource_List(GlobalSolnBlockList,
                               InputParameters.Number_of_Processors, 
	  		       InputParameters.Number_of_Blocks_Per_Processor);
    LocalSolnBlockList.allocate(InputParameters.Number_of_Blocks_Per_Processor);
    LocalSolnBlockList.ThisCPU = GlobalSolnBlockList.ThisCPU;
    Soln_ptr = Allocate(Soln_ptr, InputParameters);

    /* Loop over all initial mesh blocks and assign quadtree root
       blocks and local solution block information as required. */

    for ( j_blk = 0 ; j_blk <= InputParameters.Number_of_Blocks_Jdir-1 ; ++j_blk ) {
        for ( i_blk = 0 ; i_blk <= InputParameters.Number_of_Blocks_Idir-1 ; ++i_blk ) {
	   if (InitMeshBlk[i_blk][j_blk].Node != NULL) { // Mesh block is used!!!!
	       // Get next free solution block from list of available (not used)
	       // solution blocks.
               if (GlobalSolnBlockList.Nfree > 0) {
                  n_cpu = GlobalSolnBlockList.nextCPU();
                  n_blk = GlobalSolnBlockList.nextBlock();
                  GlobalSolnBlockList.update_next();
               } else {
                  cout << "\n" << CFFC_Version() 
                       << " AMR Error: Create_Initial_Solution_Blocks, Insufficient number of quadrilateral solution blocks.\n";
                  Soln_ptr = Deallocate(Soln_ptr, InputParameters);
                  return (NULL);
               } /* endif */

               // Assign block information to appropriate quadtree root solution block. 
               QuadTree.Roots[i_blk][j_blk].block.used = ADAPTIVEBLOCK2D_USED;
               QuadTree.Roots[i_blk][j_blk].block.gblknum = 
	          GlobalSolnBlockList.Nused-1;
               QuadTree.Roots[i_blk][j_blk].block.info.cpu = n_cpu;
               QuadTree.Roots[i_blk][j_blk].block.info.blknum = n_blk;
               QuadTree.Roots[i_blk][j_blk].block.info.dimen.i = 
                   InitMeshBlk[i_blk][j_blk].NCi-2*InputParameters.Number_of_Ghost_Cells;
               QuadTree.Roots[i_blk][j_blk].block.info.dimen.j = 
                   InitMeshBlk[i_blk][j_blk].NCj-2*InputParameters.Number_of_Ghost_Cells;
               QuadTree.Roots[i_blk][j_blk].block.info.dimen.ghost = InputParameters.Number_of_Ghost_Cells;
               QuadTree.Roots[i_blk][j_blk].block.info.sector = ADAPTIVEBLOCK2D_SECTOR_NONE;
               QuadTree.Roots[i_blk][j_blk].block.info.level = 0;
               QuadTree.Roots[i_blk][j_blk].parent_ptr = NULL;
               QuadTree.Roots[i_blk][j_blk].childNW_ptr = NULL;
               QuadTree.Roots[i_blk][j_blk].childNE_ptr = NULL;
               QuadTree.Roots[i_blk][j_blk].childSE_ptr = NULL;
               QuadTree.Roots[i_blk][j_blk].childSW_ptr = NULL;
               QuadTree.Blocks[n_cpu][n_blk] = &(QuadTree.Roots[i_blk][j_blk]);

               // For solution blocks on this processor (or processor
               // element), add block to local list, create the solution 
               // block, and copy the appropriate block of the 
               // initial quadrilateral mesh to the solution block mesh.
               if (GlobalSolnBlockList.ThisCPU == n_cpu) {
                  LocalSolnBlockList.Block[n_blk] = 
                     QuadTree.Roots[i_blk][j_blk].block;
                  Soln_ptr[n_blk].allocate(
                     LocalSolnBlockList.Block[n_blk].info.dimen.i, 
                     LocalSolnBlockList.Block[n_blk].info.dimen.j,
                     LocalSolnBlockList.Block[n_blk].info.dimen.ghost);
                  Copy_Quad_Block(Soln_ptr[n_blk].Grid, InitMeshBlk[i_blk][j_blk]);
               } /* endif */
           } /* endif */
        } /* endfor */
    } /* endfor */

    /* Renumber all solution blocks, assigning a unique global block number. */

    Renumber_Solution_Blocks(QuadTree, 
                             LocalSolnBlockList);

    /* Find the neighbours of all of the newly assigned root blocks. */

    Find_Neighbours_of_Root_Solution_Blocks(QuadTree, 
                                            LocalSolnBlockList);

    /* Modify block neighbours for grid geometries with 
       periodic boundaries, etc. */

    Modify_Neighbours_of_Root_Solution_Blocks(QuadTree, 
                                              LocalSolnBlockList,
                                              InputParameters.i_Grid);

    /* Allocates memory for all message passing buffers used  
       to send solution information between neighbouring solution blocks. */

    Allocate_Message_Buffers(LocalSolnBlockList,
                             Soln_ptr[0].NumVar()+NUM_COMP_VECTOR2D);

    /* Solution block allocation and assignment complete.  
       Return pointer to local solution blocks. */

    return(Soln_ptr);

}

/********************************************************
 * Routine: Read_QuadTree                               *
 *                                                      *
 * Reads the quadtree data structure from a file.       *
 *                                                      *
 ********************************************************/
template <class Quad_Soln_Input_Parameters>
int Read_QuadTree(QuadTreeBlock_DataStructure &QuadTree,
                  AdaptiveBlockResourceList   &List_of_Available_Blocks,
                  AdaptiveBlock2D_List        &Local_Soln_Block_List,
                  Quad_Soln_Input_Parameters  &Input_Parameters) {

    int i, nri, nrj, ncpu, nblk;
    int iBLK, jBLK;
    char quadtree_file_name[256];
    char *quadtree_file_name_ptr;
    ifstream quadtree_file;    

    /* On primary processor, determine name of quadtree input data file name. */

    if (CFFC_Primary_MPI_Processor()) {
       i = 0;
       while (1) {
          if (Input_Parameters.Restart_File_Name[i] == ' ' ||
              Input_Parameters.Restart_File_Name[i] == '.') break;
          quadtree_file_name[i]=Input_Parameters.Restart_File_Name[i];
          i = i + 1;
          if (i > strlen(Input_Parameters.Restart_File_Name) ) break;
       } /* endwhile */
       quadtree_file_name[i] = '\0';
       strcat(quadtree_file_name, "_quadtree.tree");
       quadtree_file_name_ptr = quadtree_file_name;
    } /* endif */

    /* On primary processor, open the quadtree data file. */

    if (CFFC_Primary_MPI_Processor()) {
       quadtree_file.open(quadtree_file_name_ptr, ios::in);
       if (quadtree_file.fail()) return (1);
    } /* endif */

    /* On primary processor, read in the data structure size parameters and 
       re-allocate memory as required. */

    if (CFFC_Primary_MPI_Processor()) {
       quadtree_file.setf(ios::skipws);
       quadtree_file >> nri >> nrj >> ncpu >> nblk;
       quadtree_file.unsetf(ios::skipws);

       Create_QuadTree_Data_Structure(QuadTree,
                                      nri,
  	                              nrj,
                                      Input_Parameters.Number_of_Processors,
                                      Input_Parameters.Number_of_Blocks_Per_Processor);
    } /* endif */

    /* Re-create and re-initialize the block resource list. */

    Create_Block_Resource_List(List_of_Available_Blocks,
                               Input_Parameters.Number_of_Processors, 
	  		       Input_Parameters.Number_of_Blocks_Per_Processor);

    /* On primary processor, read the quadtree data from the file. */

    if (CFFC_Primary_MPI_Processor()) {
       for ( jBLK = 0 ; jBLK <= QuadTree.NRj-1 ; ++jBLK ) {
           for ( iBLK = 0 ; iBLK <= QuadTree.NRi-1 ; ++iBLK ) {
	      QuadTree.Roots[iBLK][jBLK].read(quadtree_file,
                                              List_of_Available_Blocks);
           } /* endfor */
       } /* endfor */
    } /* endif */

    /* Broadcast quadtree data structure to all processors. */
       
    Broadcast_QuadTree_Data_Structure(QuadTree,
                                      List_of_Available_Blocks);

    /* Set the maximum and minimum refinement levels. */

    QuadTree.MaximumRefinementLevel = Input_Parameters.Maximum_Refinement_Level-1;
    QuadTree.MinimumRefinementLevel = Input_Parameters.Minimum_Refinement_Level-1;

    /* Set the thresholds for refinement and coarsening of the mesh. */

    QuadTree.RefineThreshold = Input_Parameters.Threshold_for_Refinement;
    QuadTree.CoarsenThreshold = Input_Parameters.Threshold_for_Coarsening;

    /* Re-evaluate quadtree block pointers. */

    QuadTree.assign_block_pointers();

    /* (Re-)Allocate memory for local processor solution block list. */

    if (Local_Soln_Block_List.Nblk > 0) Local_Soln_Block_List.deallocate();
    Local_Soln_Block_List.allocate(Input_Parameters.Number_of_Blocks_Per_Processor);
    Local_Soln_Block_List.ThisCPU = List_of_Available_Blocks.ThisCPU;

    /* Find the neighbours of the root blocks. */

    Find_Neighbours_of_Root_Solution_Blocks(QuadTree);

    /* Modify block neighbours for grid geometries with 
       periodic boundaries, etc. */

    Modify_Neighbours_of_Root_Solution_Blocks(QuadTree,
                                              Input_Parameters.i_Grid);

    /* Determine the neighbouring blocks of all used (active)
       solution blocks in the quadtree data structure. This will
       also copy block information to local processor solution block
       list. */

    Find_Neighbours(QuadTree,
                    Local_Soln_Block_List);

    /* On primary processor, close quadtree data file. */

    if (CFFC_Primary_MPI_Processor()) quadtree_file.close();

    /* Reading of quadtree data file complete.  Return zero value. */

    return(0);

}

/********************************************************
 * Routine: Write_QuadTree                              *
 *                                                      *
 * Writes the quadtree data structure to a file for     *
 * later retrieval.                                     *
 *                                                      *
 ********************************************************/
template <class Quad_Soln_Input_Parameters>
int Write_QuadTree(QuadTreeBlock_DataStructure &QuadTree,
                   Quad_Soln_Input_Parameters  &Input_Parameters) {

    int i;
    char quadtree_file_name[256];
    char *quadtree_file_name_ptr;
    ofstream quadtree_file;    

    /* On primary processor, determine name of quadtree output data file name. */

    if (CFFC_Primary_MPI_Processor()) {
       i = 0;
       while (1) {
          if (Input_Parameters.Restart_File_Name[i] == ' ' ||
              Input_Parameters.Restart_File_Name[i] == '.') break;
          quadtree_file_name[i]=Input_Parameters.Restart_File_Name[i];
          i = i + 1;
          if (i > strlen(Input_Parameters.Restart_File_Name) ) break;
       } /* endwhile */
       quadtree_file_name[i] = '\0';
       strcat(quadtree_file_name, "_quadtree.tree");
       quadtree_file_name_ptr = quadtree_file_name;
    } /* endif */

    /* On primary processor, open the quadtree data file. */

    if (CFFC_Primary_MPI_Processor()) {
       quadtree_file.open(quadtree_file_name_ptr, ios::out);
       if (quadtree_file.fail()) return (1);
    } /* endif */

    /* On primary processor, write the quadtree data to the file. */

    if (CFFC_Primary_MPI_Processor()) quadtree_file << QuadTree;

    /* On primary processor, close quadtree data file. */

    if (CFFC_Primary_MPI_Processor()) quadtree_file.close();

    /* Writing of quadtree data file complete.  Return zero value. */

    return(0);

}

/**********************************************************
 * Routine: Flag_Blocks_For_Refinement                    *
 *                                                        *
 * This routine flags the all local solution blocks for   *
 * refinement (coarsening or division quadrilateral       *
 * mesh solution blocks.                                  *
 *                                                        *
 **********************************************************/
template <class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
void Flag_Blocks_For_Refinement(Quad_Soln_Block             *Soln_ptr,
				Quad_Soln_Input_Parameters  &InputParameters,
                                QuadTreeBlock_DataStructure &QuadTree,
                                AdaptiveBlockResourceList   &Global_Soln_Block_List,
                                AdaptiveBlock2D_List        &Local_Soln_Block_List) {

#define	MAX_NUMBER_REFINEMENT_CRITERIA	10

    int i_criteria, i_blk;

    int number_refinement_criteria;
    double **local_block_refinement_criteria,
           *local_max_refinement_criteria, 
           *local_min_refinement_criteria,
           *global_max_refinement_criteria,
           *global_min_refinement_criteria,
           *threshold_refinement,
           *threshold_coarsening,
           *scale_refinement,
           *scale_coarsening;

    /* Allocate memory for refinement criteria. */

    local_block_refinement_criteria = new double*[Local_Soln_Block_List.Nblk];
    for ( i_blk = 0; i_blk <= Local_Soln_Block_List.Nblk-1 ; ++i_blk ) {
       local_block_refinement_criteria[i_blk] = new double[MAX_NUMBER_REFINEMENT_CRITERIA];
    } /* endfor */
    local_max_refinement_criteria = new double[MAX_NUMBER_REFINEMENT_CRITERIA];
    local_min_refinement_criteria = new double[MAX_NUMBER_REFINEMENT_CRITERIA];
    global_max_refinement_criteria = new double[MAX_NUMBER_REFINEMENT_CRITERIA];
    global_min_refinement_criteria = new double[MAX_NUMBER_REFINEMENT_CRITERIA];
    threshold_refinement = new double[MAX_NUMBER_REFINEMENT_CRITERIA];
    threshold_coarsening = new double[MAX_NUMBER_REFINEMENT_CRITERIA];
    scale_refinement = new double[MAX_NUMBER_REFINEMENT_CRITERIA];
    scale_coarsening = new double[MAX_NUMBER_REFINEMENT_CRITERIA];

    /* Set the solution block refinement flags to default values
       (no change, no division and coarsening). */

    Local_Soln_Block_List.nochangeAll();

    /* Determine and assign the refinment flag for each of the
       local solution blocks. */

    // Initial min, max, and threshold values for each refinement criteria.
    for ( i_criteria = 0; i_criteria < MAX_NUMBER_REFINEMENT_CRITERIA; ++i_criteria ) {
       local_max_refinement_criteria[i_criteria] = -1.0e-100;
       local_min_refinement_criteria[i_criteria] = 1.0e100;

       global_max_refinement_criteria[i_criteria] = -1.0e-100;
       global_min_refinement_criteria[i_criteria] = 1.0e100;
    } /* endfor */

    // Determine the min and max values of the refinement criteria
    // for all local solution blocks. */
    number_refinement_criteria = 0;
    for ( i_blk = 0 ; i_blk <= Local_Soln_Block_List.Nblk-1 ; ++i_blk ) {
       if (Local_Soln_Block_List.Block[i_blk].used) {
	 Calculate_Refinement_Criteria(local_block_refinement_criteria[i_blk],
				       InputParameters,
				       number_refinement_criteria,
				       Soln_ptr[i_blk]);
	 for ( i_criteria = 0; i_criteria < number_refinement_criteria; ++i_criteria ) {
	   local_max_refinement_criteria[i_criteria] = max(local_max_refinement_criteria[i_criteria], 
							   local_block_refinement_criteria[i_blk][i_criteria]);
	   local_min_refinement_criteria[i_criteria] = min(local_min_refinement_criteria[i_criteria], 
							   local_block_refinement_criteria[i_blk][i_criteria]);
	 } /* endfor */
       } /* endif */
    }  /* endfor */

    // Determine global values of the min and max refinement criteria.
#ifdef _MPI_VERSION
    number_refinement_criteria = CFFC_Maximum_MPI(number_refinement_criteria);
    MPI::COMM_WORLD.Allreduce(local_min_refinement_criteria, 
                              global_min_refinement_criteria, 
			      MAX_NUMBER_REFINEMENT_CRITERIA,
                              MPI::DOUBLE, 
                              MPI::MIN);
    MPI::COMM_WORLD.Allreduce(local_max_refinement_criteria, 
                              global_max_refinement_criteria, 
			      MAX_NUMBER_REFINEMENT_CRITERIA,
                              MPI::DOUBLE, 
                              MPI::MAX);
#else
    for ( i_criteria = 0; i_criteria < number_refinement_criteria; ++i_criteria ) {
       global_max_refinement_criteria[i_criteria] = local_max_refinement_criteria[i_criteria];
       global_min_refinement_criteria[i_criteria] = local_min_refinement_criteria[i_criteria];
    } /* endfor */
#endif

    /* Assign the scales and thresholds for the refinement process. */

    for ( i_criteria = 0; i_criteria < number_refinement_criteria; ++i_criteria ) {
       scale_refinement[i_criteria] = QuadTree.RefineThreshold;
       scale_coarsening[i_criteria] = QuadTree.CoarsenThreshold;
    } /* endfor */

    for ( i_criteria = 0; i_criteria < number_refinement_criteria; ++i_criteria ) {

      if (CENO_Execution_Mode::USE_CENO_ALGORITHM && 
	  CENO_Execution_Mode::USE_SMOOTHNESS_INDICATOR_FOR_AMR_CRITERIA) {	
	/* Set threshold refinement and coarsening based on an absolute scale (from 0 to 1) */
	threshold_refinement[i_criteria] = scale_refinement[i_criteria];
	threshold_coarsening[i_criteria] = scale_coarsening[i_criteria];
      } else {
      	 /* Set threshold refinement and coarsening based on an relative scale */
	 if (fabs(global_max_refinement_criteria[i_criteria]-
		  global_min_refinement_criteria[i_criteria]) > TOLER) {
	   threshold_refinement[i_criteria] = global_min_refinement_criteria[i_criteria] +
	     scale_refinement[i_criteria]*
	     (global_max_refinement_criteria[i_criteria]-
	      global_min_refinement_criteria[i_criteria]);
	   threshold_coarsening[i_criteria] = global_min_refinement_criteria[i_criteria] +
	     scale_coarsening[i_criteria]*
	     (global_max_refinement_criteria[i_criteria]-
	      global_min_refinement_criteria[i_criteria]);
	 } else {
	   threshold_refinement[i_criteria] = 0.90*global_min_refinement_criteria[i_criteria];
	   threshold_coarsening[i_criteria] = 1.10*global_min_refinement_criteria[i_criteria];
	 } /* endif */

       }/* endif */
    } /* endfor */

    /* Finally, flag the blocks for coarsening and division. */

    for ( i_blk = 0 ; i_blk <= Local_Soln_Block_List.Nblk-1 ; ++i_blk ) {
       if (Local_Soln_Block_List.Block[i_blk].used) {
          for ( i_criteria = 0; i_criteria < number_refinement_criteria; ++i_criteria ) {
             if (local_block_refinement_criteria[i_blk][i_criteria] > threshold_refinement[i_criteria] &&
                 Local_Soln_Block_List.Block[i_blk].info.level < QuadTree.MaximumRefinementLevel &&
                 (Local_Soln_Block_List.RefineFlag[i_blk] == ADAPTIVEBLOCK2D_REFINE || 
                  Local_Soln_Block_List.RefineFlag[i_blk] == ADAPTIVEBLOCK2D_NOCHANGE)) {
                Local_Soln_Block_List.RefineFlag[i_blk] = ADAPTIVEBLOCK2D_REFINE;
             } else if (local_block_refinement_criteria[i_blk][i_criteria] > threshold_refinement[i_criteria] &&
                        Local_Soln_Block_List.Block[i_blk].info.level < QuadTree.MaximumRefinementLevel &&
                        Local_Soln_Block_List.RefineFlag[i_blk] == ADAPTIVEBLOCK2D_COARSEN) {
                Local_Soln_Block_List.RefineFlag[i_blk] = ADAPTIVEBLOCK2D_REFINE;
             } else if (local_block_refinement_criteria[i_blk][i_criteria] < threshold_coarsening[i_criteria] &&
                        Local_Soln_Block_List.Block[i_blk].info.level > QuadTree.MinimumRefinementLevel &&
                        Local_Soln_Block_List.RefineFlag[i_blk] == ADAPTIVEBLOCK2D_REFINE) {
                Local_Soln_Block_List.RefineFlag[i_blk] = ADAPTIVEBLOCK2D_REFINE;
             } else if (local_block_refinement_criteria[i_blk][i_criteria] < threshold_coarsening[i_criteria] &&
                        Local_Soln_Block_List.Block[i_blk].info.level > QuadTree.MinimumRefinementLevel &&
                        (Local_Soln_Block_List.RefineFlag[i_blk] == ADAPTIVEBLOCK2D_COARSEN || 
                         Local_Soln_Block_List.RefineFlag[i_blk] == ADAPTIVEBLOCK2D_NOCHANGE)) {
                Local_Soln_Block_List.RefineFlag[i_blk] = ADAPTIVEBLOCK2D_COARSEN;
             } /* endif */
          } /* endfor */
       } /* endif */
    }  /* endfor */

    /* Deallocate memory for refinement criteria. */

    for ( i_blk = 0; i_blk <= Local_Soln_Block_List.Nblk-1 ; ++i_blk ) {
       delete []local_block_refinement_criteria[i_blk]; local_block_refinement_criteria[i_blk] = NULL;
    } /* endfor */
    delete []local_block_refinement_criteria; local_block_refinement_criteria = NULL;
    delete []local_max_refinement_criteria; local_max_refinement_criteria = NULL;
    delete []local_min_refinement_criteria; local_min_refinement_criteria = NULL;
    delete []global_max_refinement_criteria; global_max_refinement_criteria = NULL;
    delete []global_min_refinement_criteria; global_min_refinement_criteria = NULL;
    delete []threshold_refinement; threshold_refinement = NULL;
    delete []threshold_coarsening; threshold_coarsening = NULL;
    delete []scale_refinement; scale_refinement = NULL;
    delete []scale_coarsening; scale_coarsening = NULL;

}

/**********************************************************
 * Routine: Refine_Grid                                   *
 *                                                        *
 * Performs the mesh refinement of the adaptive           *
 * blocks.  Returns a zero value if no error in the mesh  *
 * refinement precedure occurred.                         *
 *                                                        *
 **********************************************************/
template <class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int Refine_Grid(Quad_Soln_Block             *Soln_ptr,
		Quad_Soln_Input_Parameters  &InputParameters,
		QuadTreeBlock_DataStructure &QuadTree,
		AdaptiveBlockResourceList   &GlobalSolnBlockList,
		AdaptiveBlock2D_List        &LocalSolnBlockList) {

    int error_flag;
    int iCPU, iBLK, iNEW, new_CPU, ii;
    int new_blocks_CPU[4], new_blocks_BLK[4], new_blocks_SECTOR[4];
    int my_rank, undefined_rank, number_CPUs_in_new_blocks, CPUs_in_new_blocks[4];

    Quad_Soln_Block solution_block_to_be_refined;
    QuadTreeBlock *quadtree_block_to_be_refined_ptr;

    AdaptiveBlock2D_Dimensions root_indices, dimensions;

#ifdef _MPI_VERSION
    MPI::Intracomm refine_comm;
    MPI::Group     big_group = MPI::COMM_WORLD.Get_group();
    MPI::Group     refine_group;
    undefined_rank = MPI::UNDEFINED;
#else
    undefined_rank = -1;
#endif

    for ( iCPU = 0 ; iCPU <= QuadTree.Ncpu-1 ; ++iCPU ) { // Loop over available processors.
        for ( iBLK = 0 ; iBLK <= QuadTree.Nblk-1 ; ++iBLK ) { // Loop over available blocks.
            if (QuadTree.Blocks[iCPU][iBLK] != NULL) {
               if (QuadTree.Blocks[iCPU][iBLK]->block.used && // Refine only solution blocks in use.
                   QuadTree.RefineFlags[iCPU][iBLK] == ADAPTIVEBLOCK2D_REFINE) { // Check each solution block for refinement.

	          // Get pointer to next quadtree block scheduled for mesh refinement.
 	          quadtree_block_to_be_refined_ptr = QuadTree.Blocks[iCPU][iBLK]; 

		  // Determine the indices of the root block.
		  root_indices = QuadTree.getRootIndices(quadtree_block_to_be_refined_ptr);

	          // Obtain list of new solution blocks to be introduced.
                  new_blocks_CPU[0] = quadtree_block_to_be_refined_ptr->block.info.cpu;
                  new_blocks_BLK[0] = quadtree_block_to_be_refined_ptr->block.info.blknum;
                  new_blocks_SECTOR[0] = ADAPTIVEBLOCK2D_SECTOR_SW;

                  number_CPUs_in_new_blocks = 1;
                  CPUs_in_new_blocks[0] = new_blocks_CPU[0];
                  for ( iNEW = 1 ; iNEW <= 3 ; ++iNEW ) {
                     CPUs_in_new_blocks[iNEW] = -1;
                  } /* endfor */
	   
                  for ( iNEW = 1 ; iNEW <= 3 ; ++iNEW ) {
                     if (GlobalSolnBlockList.Nfree > 0) {
		        // Get new blocks from global solution block list.
                        new_blocks_CPU[iNEW] = GlobalSolnBlockList.nextCPU();
                        new_blocks_BLK[iNEW] = GlobalSolnBlockList.nextBlock();
                        new_blocks_SECTOR[iNEW] = ADAPTIVEBLOCK2D_SECTOR_SW + iNEW;
                        GlobalSolnBlockList.update_next();
                     } else {
	   	        cout << "\n " << CFFC_Version() 
                             << " AMR Error: Refine_Grid, Insufficient number of solution blocks.\n";
	   	        return(6320);
                     } /* endif */
	   
                     new_CPU = 1;
	   	     for ( ii = 0 ; ii <= iNEW-1 ; ++ii) {
	   	          if (new_blocks_CPU[iNEW] == new_blocks_CPU[ii]) new_CPU = 0;
                     } /* endfor */
	   	     number_CPUs_in_new_blocks += new_CPU;
                     if (new_CPU == 1) CPUs_in_new_blocks[number_CPUs_in_new_blocks-1] = 
                        new_blocks_CPU[iNEW];
                  } /* endfor */
	   
                  // Create a MPI group and communicator for all processors 
                  // involved in the refinement of this solution block.
                  my_rank = LocalSolnBlockList.ThisCPU;
#ifdef _MPI_VERSION
                  refine_group = big_group.Incl(number_CPUs_in_new_blocks,
                                                CPUs_in_new_blocks);
                  refine_comm = MPI::COMM_WORLD.Create(refine_group);
                  my_rank = refine_group.Get_rank();
#endif
	   
	          // Obtain a copy of the solution block to be refined on each
                  // processor involved in the mesh refinement.
                  if (LocalSolnBlockList.ThisCPU == new_blocks_CPU[0]) {
	   	     Copy_Solution_Block(solution_block_to_be_refined,
                                         Soln_ptr[new_blocks_BLK[0]]);
                  } /* endif */

#ifdef _MPI_VERSION
                  if (my_rank != undefined_rank &&
                      number_CPUs_in_new_blocks > 1) {
                     Broadcast_Solution_Block(solution_block_to_be_refined,
                                              refine_comm,
                                              new_blocks_CPU[0]);
                  } /* endif */
#endif

	          // Set local solution block information, create refined mesh,
                  // and prolong solution for newly created solution blocks.
                  if (my_rank != undefined_rank) {
                     for ( iNEW = 0 ; iNEW <= 3 ; ++iNEW ) {
                        if (LocalSolnBlockList.ThisCPU == new_blocks_CPU[iNEW]) {
                           LocalSolnBlockList.Block[new_blocks_BLK[iNEW]].used = 
                              ADAPTIVEBLOCK2D_USED;
                           LocalSolnBlockList.Block[new_blocks_BLK[iNEW]].info.cpu = 
                              new_blocks_CPU[iNEW];
                           LocalSolnBlockList.Block[new_blocks_BLK[iNEW]].info.blknum = 
                              new_blocks_BLK[iNEW];
                           LocalSolnBlockList.Block[new_blocks_BLK[iNEW]].info.dimen.i = 
                              quadtree_block_to_be_refined_ptr->block.info.dimen.i;
                           LocalSolnBlockList.Block[new_blocks_BLK[iNEW]].info.dimen.j = 
                              quadtree_block_to_be_refined_ptr->block.info.dimen.j;
                           LocalSolnBlockList.Block[new_blocks_BLK[iNEW]].info.dimen.ghost =
                              quadtree_block_to_be_refined_ptr->block.info.dimen.ghost;
                           LocalSolnBlockList.Block[new_blocks_BLK[iNEW]].info.sector = 
                              new_blocks_SECTOR[iNEW];
                           LocalSolnBlockList.Block[new_blocks_BLK[iNEW]].info.level = 
                              quadtree_block_to_be_refined_ptr->block.info.level + 1;
	   
                           if ( (solution_block_to_be_refined.NCi-2*solution_block_to_be_refined.Nghost-
                                 2*((solution_block_to_be_refined.NCi-2*solution_block_to_be_refined.Nghost)/2) != 0) ||
                                (solution_block_to_be_refined.NCj-2*solution_block_to_be_refined.Nghost-
                                 2*((solution_block_to_be_refined.NCj-2*solution_block_to_be_refined.Nghost)/2) != 0) ||
                                (solution_block_to_be_refined.NCi-2*solution_block_to_be_refined.Nghost <
                                 2*solution_block_to_be_refined.Nghost) ||
                                (solution_block_to_be_refined.NCj-2*solution_block_to_be_refined.Nghost <
                                 2*solution_block_to_be_refined.Nghost) ||
                                (solution_block_to_be_refined.Grid.Node == NULL) ) {
                              cout << "\n " << CFFC_Version() 
                                   << " AMR Error: Refine_Grid, Cannot refine mesh.\n";
                              return(6321);
                           } /* endif */
	   
                           if (Soln_ptr[new_blocks_BLK[iNEW]].U != NULL) 
                             Soln_ptr[new_blocks_BLK[iNEW]].deallocate();
                           Soln_ptr[new_blocks_BLK[iNEW]].allocate(
                                       solution_block_to_be_refined.NCi-2*solution_block_to_be_refined.Nghost,
                                       solution_block_to_be_refined.NCj-2*solution_block_to_be_refined.Nghost,
                                       solution_block_to_be_refined.Nghost);
                           Refine_Mesh(Soln_ptr[new_blocks_BLK[iNEW]].Grid,
                                       solution_block_to_be_refined.Grid,
                                       new_blocks_SECTOR[iNEW]);
			   if (InputParameters.i_Smooth_Quad_Block) {
                               Smooth_Quad_Block(Soln_ptr[new_blocks_BLK[iNEW]].Grid,
                                                 min(250,
						 4*max(solution_block_to_be_refined.NCi-2*solution_block_to_be_refined.Nghost,
						       solution_block_to_be_refined.NCj-2*solution_block_to_be_refined.Nghost)));
			   } /* endif */
                           if (Soln_ptr[new_blocks_BLK[iNEW]].Grid.Check_Quad_Block()) { 
                              cout << "\n " << CFFC_Version() 
                                   << " AMR Error: Refine_Grid, Invalid refined mesh.\n";
                              return(6322);
                           } /* endif */
                           error_flag = Prolong_Solution_Block(Soln_ptr[new_blocks_BLK[iNEW]],
							       solution_block_to_be_refined,
							       new_blocks_SECTOR[iNEW]);
			   if (error_flag) return error_flag;
                        } /* endif */
                     } /* endfor */
                  } /* endif */

                  // Release the MPI group and communicator.
#ifdef _MPI_VERSION
                  if (refine_comm != MPI::COMM_NULL) refine_comm.Free();
                  refine_group.Free();
#endif
                  my_rank = undefined_rank;
	   
	          // Update and assign quadtree solution block information for
                  // newly created solution blocks.
                  QuadTree.refineBlock(new_blocks_CPU, 
                                       new_blocks_BLK,
                                       new_blocks_SECTOR);

                  // Finally, reset refinement flag.
	          QuadTree.RefineFlags[iCPU][iBLK] = ADAPTIVEBLOCK2D_NOCHANGE;

	       // Block not used, reset refinement flag.
               } else if (!QuadTree.Blocks[iCPU][iBLK]->block.used) {
	          QuadTree.RefineFlags[iCPU][iBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
               } /* endif */

	    // Block does not exist, resent refinement flag.
            } else {
	       QuadTree.RefineFlags[iCPU][iBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
            } /* endif */
        } /* endfor */
    } /* endfor */    

    /* Mesh refinement complete.  Return zero value. */

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
template <class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int Coarsen_Grid(Quad_Soln_Block             *Soln_ptr,
		 Quad_Soln_Input_Parameters  &InputParameters,
                 QuadTreeBlock_DataStructure &QuadTree,
                 AdaptiveBlockResourceList   &GlobalSolnBlockList,
                 AdaptiveBlock2D_List        &LocalSolnBlockList) {

    int error_flag;
    int iCPU, iBLK, iOLD;
    int old_blocks_CPU[4], old_blocks_BLK[4], old_blocks_SECTOR[4], CPU_list[2];
    int my_rank, undefined_rank, all_coarsen_flag;

    Quad_Soln_Block solution_block_to_be_coarsened_SW_sibling,
                    solution_block_to_be_coarsened_SE_sibling,
                    solution_block_to_be_coarsened_NW_sibling,
                    solution_block_to_be_coarsened_NE_sibling;
    QuadTreeBlock *quadtree_block_to_be_coarsened_ptr; 

    AdaptiveBlock2D_Dimensions root_indices, dimensions;

#ifdef _MPI_VERSION
    MPI::Intracomm coarsening_comm_SE, coarsening_comm_NW, coarsening_comm_NE;
    MPI::Group     big_group = MPI::COMM_WORLD.Get_group();
    MPI::Group     coarsening_group_SE, coarsening_group_NW, coarsening_group_NE;
    undefined_rank = MPI::UNDEFINED;
#else
    undefined_rank = -1;
#endif

    for ( iCPU = 0 ; iCPU <= QuadTree.Ncpu-1 ; ++iCPU ) { // Loop over available processors.
        for ( iBLK = 0 ; iBLK <= QuadTree.Nblk-1 ; ++iBLK ) { // Loop over available blocks.
            if (QuadTree.Blocks[iCPU][iBLK] != NULL) {
               if (QuadTree.Blocks[iCPU][iBLK]->block.used && // Coarsen only solution blocks in use.
                   QuadTree.RefineFlags[iCPU][iBLK] == ADAPTIVEBLOCK2D_COARSEN) { // Check each solution block for coarsening.
	   
	          // Get pointer to next quadtree block scheduled for mesh coarsening.
 	          quadtree_block_to_be_coarsened_ptr = QuadTree.Blocks[iCPU][iBLK];

		  // Determine the indices of the root block.
		  root_indices = QuadTree.getRootIndices(quadtree_block_to_be_coarsened_ptr->parent_ptr);

                  // If parent exits, then try to coarsen the block.
                  if (quadtree_block_to_be_coarsened_ptr->parent_ptr != NULL) {
		     // Ensure all sibling blocks to be coarsened are in use.
		     if (quadtree_block_to_be_coarsened_ptr->parent_ptr->childNW_ptr != NULL &&
                         quadtree_block_to_be_coarsened_ptr->parent_ptr->childNE_ptr != NULL &&
                         quadtree_block_to_be_coarsened_ptr->parent_ptr->childSE_ptr != NULL &&
                         quadtree_block_to_be_coarsened_ptr->parent_ptr->childSW_ptr != NULL) {
                        if (quadtree_block_to_be_coarsened_ptr->parent_ptr->childNW_ptr->block.used &&
                            quadtree_block_to_be_coarsened_ptr->parent_ptr->childNE_ptr->block.used &&
                            quadtree_block_to_be_coarsened_ptr->parent_ptr->childSE_ptr->block.used &&
                            quadtree_block_to_be_coarsened_ptr->parent_ptr->childSW_ptr->block.used) {
			    // Obtain list of siblings for coarsening.
                            old_blocks_CPU[0] = quadtree_block_to_be_coarsened_ptr->parent_ptr->childSW_ptr->block.info.cpu;
                            old_blocks_BLK[0] = quadtree_block_to_be_coarsened_ptr->parent_ptr->childSW_ptr->block.info.blknum;
                            old_blocks_SECTOR[0] = quadtree_block_to_be_coarsened_ptr->parent_ptr->childSW_ptr->block.info.sector;

                            old_blocks_CPU[1] = quadtree_block_to_be_coarsened_ptr->parent_ptr->childSE_ptr->block.info.cpu;
                            old_blocks_BLK[1] = quadtree_block_to_be_coarsened_ptr->parent_ptr->childSE_ptr->block.info.blknum;
                            old_blocks_SECTOR[1] = quadtree_block_to_be_coarsened_ptr->parent_ptr->childSE_ptr->block.info.sector;

                            old_blocks_CPU[2] = quadtree_block_to_be_coarsened_ptr->parent_ptr->childNW_ptr->block.info.cpu;
                            old_blocks_BLK[2] = quadtree_block_to_be_coarsened_ptr->parent_ptr->childNW_ptr->block.info.blknum;
                            old_blocks_SECTOR[2] = quadtree_block_to_be_coarsened_ptr->parent_ptr->childNW_ptr->block.info.sector;
   
                            old_blocks_CPU[3] = quadtree_block_to_be_coarsened_ptr->parent_ptr->childNE_ptr->block.info.cpu;
                            old_blocks_BLK[3] = quadtree_block_to_be_coarsened_ptr->parent_ptr->childNE_ptr->block.info.blknum;
                            old_blocks_SECTOR[3] = quadtree_block_to_be_coarsened_ptr->parent_ptr->childNE_ptr->block.info.sector;

                            // Check to see that all siblings are flagged for coarsening.
                            all_coarsen_flag = 1;
                            for ( iOLD = 0 ; iOLD <= 3 ; ++iOLD ) {
			       if (QuadTree.RefineFlags[old_blocks_CPU[iOLD]][old_blocks_BLK[iOLD]] != 
                                   ADAPTIVEBLOCK2D_COARSEN) {
                                  all_coarsen_flag = 0;
                                  break;
                               } /* endif */
                            } /* endfor */
                            if (all_coarsen_flag) {
                               // Make copies of each of the four sibling solution blocks involved 
                               // in the coarsening process.
                               if (LocalSolnBlockList.ThisCPU == old_blocks_CPU[0]) {
                                  Copy_Solution_Block(solution_block_to_be_coarsened_SW_sibling,
                                                      Soln_ptr[old_blocks_BLK[0]]);
                               } /* endif */

                               if (LocalSolnBlockList.ThisCPU == old_blocks_CPU[1]) {
                                  Copy_Solution_Block(solution_block_to_be_coarsened_SE_sibling,
                                                      Soln_ptr[old_blocks_BLK[1]]);
                               } /* endif */

                               if (LocalSolnBlockList.ThisCPU == old_blocks_CPU[2]) {
                                  Copy_Solution_Block(solution_block_to_be_coarsened_NW_sibling,
                                                      Soln_ptr[old_blocks_BLK[2]]);
                               } /* endif */

                               if (LocalSolnBlockList.ThisCPU == old_blocks_CPU[3]) {
                                  Copy_Solution_Block(solution_block_to_be_coarsened_NE_sibling,
                                                      Soln_ptr[old_blocks_BLK[3]]);
                               } /* endif */
      
                               // Create MPI groups and communicators for all processors 
                               // involved in the coarsening of the four sibling solution blocks.
#ifdef _MPI_VERSION
                               if (old_blocks_CPU[0] != old_blocks_CPU[1]) {
                                  CPU_list[0] = old_blocks_CPU[1];
                                  CPU_list[1] = old_blocks_CPU[0];
                                  coarsening_group_SE = big_group.Incl(2,
                                                                       CPU_list);
                                  coarsening_comm_SE = MPI::COMM_WORLD.Create(coarsening_group_SE);
                               } /* endif */
                               if (old_blocks_CPU[0] != old_blocks_CPU[2]) {
                                  CPU_list[0] = old_blocks_CPU[2];
                                  CPU_list[1] = old_blocks_CPU[0];
                                  coarsening_group_NW = big_group.Incl(2,
                                                                       CPU_list);
                                  coarsening_comm_NW = MPI::COMM_WORLD.Create(coarsening_group_NW);
                               } /* endif */
                               if (old_blocks_CPU[0] != old_blocks_CPU[3]) {
                                  CPU_list[0] = old_blocks_CPU[3];
                                  CPU_list[1] = old_blocks_CPU[0];
                                  coarsening_group_NE = big_group.Incl(2,
                                                                       CPU_list);
                                  coarsening_comm_NE = MPI::COMM_WORLD.Create(coarsening_group_NE);
                               } /* endif */
#endif

                               // Obtain a copy of the sibling solution blocks on the
                               // processor involved in the mesh coarsening.
#ifdef _MPI_VERSION
                               if (old_blocks_CPU[0] != old_blocks_CPU[1]) {
                                  my_rank = coarsening_group_SE.Get_rank();
                                  if (my_rank != undefined_rank) {
                                     Broadcast_Solution_Block(solution_block_to_be_coarsened_SE_sibling,
                                                              coarsening_comm_SE,
                                                              old_blocks_CPU[1]);
                                  } /* endif */
                               } /* endif */
                               if (old_blocks_CPU[0] != old_blocks_CPU[2]) {
                                  my_rank = coarsening_group_NW.Get_rank();
                                  if (my_rank != undefined_rank) {
                                     Broadcast_Solution_Block(solution_block_to_be_coarsened_NW_sibling,
                                                              coarsening_comm_NW,
                                                              old_blocks_CPU[2]);
                                  } /* endif */
                               } /* endif */
                               if (old_blocks_CPU[0] != old_blocks_CPU[3]) {
                                  my_rank = coarsening_group_NE.Get_rank();
                                  if (my_rank != undefined_rank) {
                                     Broadcast_Solution_Block(solution_block_to_be_coarsened_NE_sibling,
                                                              coarsening_comm_NE,
                                                              old_blocks_CPU[3]);
                                  } /* endif */
                               } /* endif */
#endif

                               // Set local solution block information, create coarse mesh,
                               // and restrict solution to newly created coarse solution block.
                               if (LocalSolnBlockList.ThisCPU == old_blocks_CPU[0]) {
                                  LocalSolnBlockList.Block[old_blocks_BLK[0]].used = 
                                     ADAPTIVEBLOCK2D_USED;
                                  LocalSolnBlockList.Block[old_blocks_BLK[0]].info.cpu = 
                                     old_blocks_CPU[0];
                                  LocalSolnBlockList.Block[old_blocks_BLK[0]].info.blknum = 
                                     old_blocks_BLK[0];
                                  LocalSolnBlockList.Block[old_blocks_BLK[0]].info.dimen.i = 
                                     QuadTree.Blocks[old_blocks_CPU[0]][old_blocks_BLK[0]]->block.info.dimen.i;
                                  LocalSolnBlockList.Block[old_blocks_BLK[0]].info.dimen.j = 
                                     QuadTree.Blocks[old_blocks_CPU[0]][old_blocks_BLK[0]]->block.info.dimen.j;
                                  LocalSolnBlockList.Block[old_blocks_BLK[0]].info.dimen.ghost =
                                     QuadTree.Blocks[old_blocks_CPU[0]][old_blocks_BLK[0]]->block.info.dimen.ghost;
                                  LocalSolnBlockList.Block[old_blocks_BLK[0]].info.sector = 
				     old_blocks_SECTOR[0]; // Note that this sector info is incorrect, 
                                                           // but subsequent call to Find_Neighbours will fix this!
                                  LocalSolnBlockList.Block[old_blocks_BLK[0]].info.level = 
                                     QuadTree.Blocks[old_blocks_CPU[0]][old_blocks_BLK[0]]->block.info.level - 1;
	   
                                  if ( (solution_block_to_be_coarsened_SW_sibling.NCi-
					2*solution_block_to_be_coarsened_SW_sibling.Nghost-
                                        2*((solution_block_to_be_coarsened_SW_sibling.NCi-
					    2*solution_block_to_be_coarsened_SW_sibling.Nghost)/2) != 0) ||
                                       (solution_block_to_be_coarsened_SW_sibling.NCj-
					2*solution_block_to_be_coarsened_SW_sibling.Nghost-
                                        2*((solution_block_to_be_coarsened_SW_sibling.NCj-
					    2*solution_block_to_be_coarsened_SW_sibling.Nghost)/2) != 0) ||
                                       (solution_block_to_be_coarsened_SW_sibling.NCi-
					2*solution_block_to_be_coarsened_SW_sibling.Nghost < 
                                        2*solution_block_to_be_coarsened_SW_sibling.Nghost) ||
                                       (solution_block_to_be_coarsened_SW_sibling.NCj-
					2*solution_block_to_be_coarsened_SW_sibling.Nghost < 
                                        2*solution_block_to_be_coarsened_SW_sibling.Nghost) ||
                                       (solution_block_to_be_coarsened_SW_sibling.Grid.Node == NULL) ) {
                                     cout << "\n " << CFFC_Version()
                                          << " AMR Error: Coarsen_Grid, Cannot coarsen south-west mesh.\n";
                                     return(7480);
                                  } /* endif */

                                  if ( (solution_block_to_be_coarsened_SE_sibling.NCi-
					2*solution_block_to_be_coarsened_SE_sibling.Nghost-
                                        2*((solution_block_to_be_coarsened_SE_sibling.NCi-
					    2*solution_block_to_be_coarsened_SE_sibling.Nghost)/2) != 0) ||
                                       (solution_block_to_be_coarsened_SE_sibling.NCj-
					2*solution_block_to_be_coarsened_SE_sibling.Nghost-
                                        2*((solution_block_to_be_coarsened_SE_sibling.NCj-
					    2*solution_block_to_be_coarsened_SE_sibling.Nghost)/2) != 0) ||
                                       (solution_block_to_be_coarsened_SE_sibling.NCi-
					2*solution_block_to_be_coarsened_SE_sibling.Nghost < 
                                        2*solution_block_to_be_coarsened_SE_sibling.Nghost) ||
                                       (solution_block_to_be_coarsened_SE_sibling.NCj-
					2*solution_block_to_be_coarsened_SE_sibling.Nghost < 
                                        2*solution_block_to_be_coarsened_SE_sibling.Nghost) ||
                                       (solution_block_to_be_coarsened_SE_sibling.Grid.Node == NULL) ) {
                                     cout << "\n " << CFFC_Version()
                                          << " AMR Error: Coarsen_Grid, Cannot coarsen south-east mesh.\n";
                                     return(7481);
                                  } /* endif */

                                  if ( (solution_block_to_be_coarsened_NW_sibling.NCi-
					2*solution_block_to_be_coarsened_NW_sibling.Nghost-
                                        2*((solution_block_to_be_coarsened_NW_sibling.NCi-
					    2*solution_block_to_be_coarsened_NW_sibling.Nghost)/2) != 0) ||
                                       (solution_block_to_be_coarsened_NW_sibling.NCj-
					2*solution_block_to_be_coarsened_NW_sibling.Nghost-
                                        2*((solution_block_to_be_coarsened_NW_sibling.NCj-
					    2*solution_block_to_be_coarsened_NW_sibling.Nghost)/2) != 0) ||
                                       (solution_block_to_be_coarsened_NW_sibling.NCi-
					2*solution_block_to_be_coarsened_NW_sibling.Nghost < 
                                        2*solution_block_to_be_coarsened_NW_sibling.Nghost) ||
                                       (solution_block_to_be_coarsened_NW_sibling.NCj-
					2*solution_block_to_be_coarsened_NW_sibling.Nghost < 
                                        2*solution_block_to_be_coarsened_NW_sibling.Nghost) ||
                                       (solution_block_to_be_coarsened_NW_sibling.Grid.Node == NULL) ) {
                                     cout << "\n " << CFFC_Version()
                                          << " AMR Error: Coarsen_Grid, Cannot coarsen north-west mesh.\n";
                                     return(7482);
                                  } /* endif */

                                 if ( (solution_block_to_be_coarsened_NE_sibling.NCi-
					2*solution_block_to_be_coarsened_NE_sibling.Nghost-
                                        2*((solution_block_to_be_coarsened_NE_sibling.NCi-
					    2*solution_block_to_be_coarsened_NE_sibling.Nghost)/2) != 0) ||
                                       (solution_block_to_be_coarsened_NE_sibling.NCj-
					2*solution_block_to_be_coarsened_NE_sibling.Nghost-
                                        2*((solution_block_to_be_coarsened_NE_sibling.NCj-
					    2*solution_block_to_be_coarsened_NE_sibling.Nghost)/2) != 0) ||
                                       (solution_block_to_be_coarsened_NE_sibling.NCi-
					2*solution_block_to_be_coarsened_NE_sibling.Nghost < 
                                        2*solution_block_to_be_coarsened_NE_sibling.Nghost) ||
                                       (solution_block_to_be_coarsened_NE_sibling.NCj-
					2*solution_block_to_be_coarsened_NE_sibling.Nghost < 
                                        2*solution_block_to_be_coarsened_NE_sibling.Nghost) ||
                                       (solution_block_to_be_coarsened_NE_sibling.Grid.Node == NULL) ) {
                                     cout << "\n " << CFFC_Version()
                                          << " AMR Error: Coarsen_Grid, Cannot coarsen north-east mesh.\n";
                                     return(7483);
                                  } /* endif */

                                  if (Soln_ptr[old_blocks_BLK[0]].U != NULL) 
                                     Soln_ptr[old_blocks_BLK[0]].deallocate();
                                  Soln_ptr[old_blocks_BLK[0]].allocate(
                                           solution_block_to_be_coarsened_SW_sibling.NCi-
					   2*solution_block_to_be_coarsened_SW_sibling.Nghost,
                                           solution_block_to_be_coarsened_SW_sibling.NCj-
					   2*solution_block_to_be_coarsened_SW_sibling.Nghost,
                                           solution_block_to_be_coarsened_SW_sibling.Nghost);
                                  Coarsen_Mesh(Soln_ptr[old_blocks_BLK[0]].Grid,
                                               solution_block_to_be_coarsened_SW_sibling.Grid,
                                               solution_block_to_be_coarsened_SE_sibling.Grid,
                                               solution_block_to_be_coarsened_NW_sibling.Grid,
                                               solution_block_to_be_coarsened_NE_sibling.Grid);
				  if (InputParameters.i_Smooth_Quad_Block) {
                                      Smooth_Quad_Block(Soln_ptr[old_blocks_BLK[0]].Grid,
                                                        min(250, 4*max(
                                                        solution_block_to_be_coarsened_SW_sibling.NCi-
						        2*solution_block_to_be_coarsened_SW_sibling.Nghost,
                                                        solution_block_to_be_coarsened_SW_sibling.NCj-
						        2*solution_block_to_be_coarsened_SW_sibling.Nghost)));
				  } /* endif */
                                  if (Soln_ptr[old_blocks_BLK[0]].Grid.Check_Quad_Block()) { 
                                     cout << "\n " << CFFC_Version() 
                                          << " AMR Error: Coarsen_Grid, Invalid coarsened mesh.\n";
                                     return(7484);
                                  } /* endif */
                                  error_flag = Restrict_Solution_Block(Soln_ptr[old_blocks_BLK[0]],
								       solution_block_to_be_coarsened_SW_sibling,
								       solution_block_to_be_coarsened_SE_sibling,
								       solution_block_to_be_coarsened_NW_sibling,
								       solution_block_to_be_coarsened_NE_sibling);
				  if (error_flag) return error_flag;
                               } /* endif */

                               // Release the MPI groups and communicators.
#ifdef _MPI_VERSION
                               if (old_blocks_CPU[0] != old_blocks_CPU[1]) {
                                  if (coarsening_comm_SE != MPI::COMM_NULL) coarsening_comm_SE.Free();
                                  coarsening_group_SE.Free();
                               } /* endif */
                               if (old_blocks_CPU[0] != old_blocks_CPU[2]) {
                                  if (coarsening_comm_NW != MPI::COMM_NULL) coarsening_comm_NW.Free();
                                  coarsening_group_NW.Free();
                               } /* endif */
                               if (old_blocks_CPU[0] != old_blocks_CPU[3]) {
                                  if (coarsening_comm_NE != MPI::COMM_NULL) coarsening_comm_NE.Free();
                                  coarsening_group_NE.Free();
                               } /* endif */
#endif

                               // Free old unused solution blocks and return to global
                               // solution block list for re-use.                          
                               for ( iOLD = 3 ; iOLD >= 1 ; --iOLD ) {
                                  if (LocalSolnBlockList.ThisCPU == old_blocks_CPU[iOLD]) {
                                     LocalSolnBlockList.Block[old_blocks_BLK[iOLD]].used = ADAPTIVEBLOCK2D_NOT_USED;
                                     if (Soln_ptr[old_blocks_BLK[iOLD]].U != NULL) Soln_ptr[old_blocks_BLK[iOLD]].deallocate();
                                  } /* endif */
                                  // Return blocks and update global solution block list.
                                  GlobalSolnBlockList.returnCPU(old_blocks_CPU[iOLD]);
                                  GlobalSolnBlockList.returnBlock(old_blocks_BLK[iOLD]);
                                  GlobalSolnBlockList.update_return();
                               } /* endfor */

                               // Update and assign quadtree solution block information for
                               // coarsened solution block, deleting old refined blocks.
                               QuadTree.coarsenBlocks(old_blocks_CPU, 
                                                      old_blocks_BLK,
                                                      old_blocks_SECTOR);

                               // Finally, reset refinement flags for four blocks involved in coarsening.
                               for ( iOLD = 0 ; iOLD <= 3 ; ++iOLD ) {
			          QuadTree.RefineFlags[old_blocks_CPU[iOLD]][old_blocks_BLK[iOLD]] = ADAPTIVEBLOCK2D_NOCHANGE;
                               } /* endfor */

                            // If siblings aren't all flagged for coarsening, reset refinement flags.
                            } else {
                               QuadTree.RefineFlags[iCPU][iBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
                               if (LocalSolnBlockList.ThisCPU == iCPU)
				 LocalSolnBlockList.RefineFlag[iBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
                            } /* endif */

                        // If siblings aren't all used, reset refinement flag.
                        } else {
                           QuadTree.RefineFlags[iCPU][iBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
			   if (LocalSolnBlockList.ThisCPU == iCPU)
			     LocalSolnBlockList.RefineFlag[iBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
                        } /* endif */
                     
                     // If siblings don't all exist, reset refinement flag.
                     } else {
                        QuadTree.RefineFlags[iCPU][iBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
			if (LocalSolnBlockList.ThisCPU == iCPU)
			  LocalSolnBlockList.RefineFlag[iBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
                     } /* endif */

                  // If parent doesn't exist, reset refinement flag.
                  } else {
                     QuadTree.RefineFlags[iCPU][iBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
		     if (LocalSolnBlockList.ThisCPU == iCPU)
		       LocalSolnBlockList.RefineFlag[iBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
                  } /* endif */

	       // Block not used, reset refinement flag.
               } else if (!QuadTree.Blocks[iCPU][iBLK]->block.used) {
	          QuadTree.RefineFlags[iCPU][iBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
		  if (LocalSolnBlockList.ThisCPU == iCPU)
		    LocalSolnBlockList.RefineFlag[iBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
               } /* endif */

	    // Block does not exist, reset refinement flag.
            } else {
	       QuadTree.RefineFlags[iCPU][iBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
	       if (LocalSolnBlockList.ThisCPU == iCPU)
		 LocalSolnBlockList.RefineFlag[iBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
            } /* endif */

        } /* endfor */
    } /* endfor */    

    /* Mesh coarsening complete.  Return zero value. */

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
template <class Quad_Soln_Block>
void Fix_Refined_Block_Boundaries(Quad_Soln_Block      *Soln_ptr,
                                  AdaptiveBlock2D_List &Soln_Block_List) {

    int i_blk, 
        fix_north_boundary, fix_south_boundary,
        fix_east_boundary, fix_west_boundary;

    /* Fix mesh boundaries of each refined solution block. */
    for ( i_blk = 0 ; i_blk <= Soln_Block_List.Nblk-1 ; ++i_blk ) {
       if (Soln_Block_List.Block[i_blk].used) {
          if (Soln_Block_List.Block[i_blk].nN == 1 && 
              (Soln_Block_List.Block[i_blk].info.level > 
               Soln_Block_List.Block[i_blk].infoN[0].level)) {
             fix_north_boundary = 1;
          } else {
             fix_north_boundary = 0;
          } /* endif */
          if (Soln_Block_List.Block[i_blk].nS == 1 && 
              (Soln_Block_List.Block[i_blk].info.level > 
               Soln_Block_List.Block[i_blk].infoS[0].level)) {
             fix_south_boundary = 1;
          } else {
             fix_south_boundary = 0;
          } /* endif */
          if (Soln_Block_List.Block[i_blk].nE == 1 && 
              (Soln_Block_List.Block[i_blk].info.level > 
               Soln_Block_List.Block[i_blk].infoE[0].level)) {
             fix_east_boundary = 1;
          } else {
             fix_east_boundary = 0;
          } /* endif */
          if (Soln_Block_List.Block[i_blk].nW == 1 && 
              (Soln_Block_List.Block[i_blk].info.level > 
               Soln_Block_List.Block[i_blk].infoW[0].level)) {
             fix_west_boundary = 1;
          } else {
             fix_west_boundary = 0;
          } /* endif */
          Fix_Refined_Block_Boundaries(Soln_ptr[i_blk],
                                       fix_north_boundary,
                                       fix_south_boundary,
                                       fix_east_boundary,
                                       fix_west_boundary);
       } /* endif */
    }  /* endfor */

}

/********************************************************
 * Routine: Unfix_Refined_Block_Boundaries              *
 *                                                      *
 * Returns the locations of the boundary nodes of a     *
 * solution block to their original unmodified          *
 * positions.                                           *
 *                                                      *
 ********************************************************/
template <class Quad_Soln_Block>
void Unfix_Refined_Block_Boundaries(Quad_Soln_Block      *Soln_ptr,
                                    AdaptiveBlock2D_List &Soln_Block_List) {

    int i_blk; 

    for ( i_blk = 0 ; i_blk <= Soln_Block_List.Nblk-1 ; ++i_blk ) {
       if (Soln_Block_List.Block[i_blk].used) {
          Unfix_Refined_Block_Boundaries(Soln_ptr[i_blk]);
       } /* endif */
    }  /* endfor */

}

/**********************************************************
 * Routine: AMR                                           *
 *                                                        *
 * Performs the adaptive mesh refinement (AMR) of the     *
 * adaptive solution blocks.  Returns a zero value if no  *
 * error in the AMR precedure occurred.                   *
 *                                                        *
 **********************************************************/
template <class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int AMR(Quad_Soln_Block             *Soln_ptr,
	Quad_Soln_Input_Parameters  &InputParameters,
        QuadTreeBlock_DataStructure &QuadTree,
        AdaptiveBlockResourceList   &GlobalSolnBlockList,
        AdaptiveBlock2D_List        &LocalSolnBlockList,
        const int                    Set_New_Refinement_Flags,
	const int                    Apply_BCs_Flag) {

    int error_flag;

    /* Calculate the refinement measures for each solution
       block and flag local solution blocks for refinement
       or coarsening. */

    if (Set_New_Refinement_Flags) {
       Flag_Blocks_For_Refinement(Soln_ptr,
				  InputParameters,
                                  QuadTree,
                                  GlobalSolnBlockList,
                                  LocalSolnBlockList);
    } /* endif */

    /* Get global refinement list and adjust refinement and
       derefinement flags so as to ensure a permissible grid
       block topology. */

    Get_Refinement_List(QuadTree,
                        LocalSolnBlockList);

    /* Return the locations of the boundary nodes of the
       solution blocks to their unmodified locations. */

    Unfix_Refined_Block_Boundaries(Soln_ptr,
                                   LocalSolnBlockList);

    /* First, coarsen solutions blocks and return unused solution
       blocks to pool of available resources for subsequent
       refinement. */

    error_flag = Coarsen_Grid(Soln_ptr,
			      InputParameters,
                              QuadTree,
                              GlobalSolnBlockList,
                              LocalSolnBlockList);
    if (error_flag) return (error_flag);

    /* Refine solution blocks. */

    error_flag = Refine_Grid(Soln_ptr,
			     InputParameters,
                             QuadTree,
                             GlobalSolnBlockList,
                             LocalSolnBlockList);
    if (error_flag) return (error_flag);

    /* Renumber all solution blocks, assigning new global block
       numbers. */

    Renumber_Solution_Blocks(QuadTree, 
                             LocalSolnBlockList);

    /* Determine the neighbouring blocks of all used (active)
       solution blocks in the quadtree data structure. */

    Find_Neighbours(QuadTree,
                    LocalSolnBlockList);

    /* Adjust the locations of the boundary nodes of the
       solution blocks so that the new node locations match 
       with cell volumes of adjacent solution blocks that 
       have lower levels of mesh refinement (i.e., are 
       coarser solution blocks). */

    Fix_Refined_Block_Boundaries(Soln_ptr,
				 LocalSolnBlockList);

    /* Re-allocate memory for all message passing buffers used  
       to send solution information between neighbouring solution 
       blocks. */

    Allocate_Message_Buffers(LocalSolnBlockList,
                             Soln_ptr[0].NumVar()+NUM_COMP_VECTOR2D);

    /* Update solution information shared between neighbouring blocks. */

    error_flag = Send_All_Messages(Soln_ptr,
                                   LocalSolnBlockList,
                                   NUM_COMP_VECTOR2D,
                                   ON);
    if (!error_flag) error_flag = Send_All_Messages(Soln_ptr,
                                                    LocalSolnBlockList,
                                                    Soln_ptr[0].NumVar(),
                                                    OFF);
    if (error_flag) return (error_flag);

    /* Re-prescribe boundary data consistent with newly refined
       and coarsened solution blocks. */

    if (Apply_BCs_Flag) {
      BCs(Soln_ptr,
	  LocalSolnBlockList,
	  InputParameters);
    }

    /* AMR procedure successfully completed.  Return zero value. */

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
template <class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int Initial_AMR(Quad_Soln_Block              *Soln_ptr,
                Quad_Soln_Input_Parameters   &InputParameters,
                QuadTreeBlock_DataStructure  &QuadTree,
                AdaptiveBlockResourceList    &GlobalSolnBlockList,
                AdaptiveBlock2D_List         &LocalSolnBlockList) {

    int error_flag, number_of_initial_mesh_refinements;
    double original_coarsenthreshold;

    // Do not allow coarsening during initial refinement of mesh.
    original_coarsenthreshold = QuadTree.CoarsenThreshold;
    QuadTree.CoarsenThreshold = ZERO;

    error_flag = 0;

    if (InputParameters.Number_of_Initial_Mesh_Refinements == 0) {
      QuadTree.CoarsenThreshold = original_coarsenthreshold;
      return(error_flag);
    }

    for (number_of_initial_mesh_refinements = 1; 
         number_of_initial_mesh_refinements <= InputParameters.Number_of_Initial_Mesh_Refinements; 
         ++number_of_initial_mesh_refinements) {

       error_flag = AMR(Soln_ptr,
			InputParameters,
   	                QuadTree,
 	                GlobalSolnBlockList,
	                LocalSolnBlockList,
	                ON,ON);
       if (error_flag) return (error_flag);

       ICs(Soln_ptr,
           LocalSolnBlockList,
           InputParameters);

       BCs(Soln_ptr,
           LocalSolnBlockList,
	   InputParameters);

       error_flag = Send_All_Messages(Soln_ptr,
                                      LocalSolnBlockList,
                                      Soln_ptr[0].NumVar(),
                                      OFF);
       if (error_flag) return (error_flag);

       if (CFFC_Primary_MPI_Processor()) {
          cout << "\n Refinement Level #" << number_of_initial_mesh_refinements
               << " : Number of Blocks = " << QuadTree.countUsedBlocks()
               << ", Number of Cells = " << QuadTree.countUsedCells()
               << ", Refinement Efficiency = " << QuadTree.efficiencyRefinement();
       } /* endif */

    } /* endfor */

    QuadTree.CoarsenThreshold = original_coarsenthreshold;

    return(error_flag);

}

/**********************************************************
 * Routine: Uniform_AMR                                   *
 *                                                        *
 * Performs uniform refinement of the adaptive solution   *
 * block mesh.  Returns a zero value if no error in the   *
 * AMR precedure has occurred.                            *
 *                                                        *
 **********************************************************/
template <class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int Uniform_AMR(Quad_Soln_Block              *Soln_ptr,
                Quad_Soln_Input_Parameters   &InputParameters,
                QuadTreeBlock_DataStructure  &QuadTree,
                AdaptiveBlockResourceList    &GlobalSolnBlockList,
                AdaptiveBlock2D_List         &LocalSolnBlockList) {

    int error_flag, number_of_uniform_mesh_refinements;

    error_flag = 0;

    if (InputParameters.Number_of_Uniform_Mesh_Refinements == 0) return(error_flag);

    for (number_of_uniform_mesh_refinements = 1; 
         number_of_uniform_mesh_refinements <= InputParameters.Number_of_Uniform_Mesh_Refinements; 
         ++number_of_uniform_mesh_refinements) {

       // Set refinement flags to all.
       QuadTree.nochangeAll();
       LocalSolnBlockList.refineAll();

       error_flag = AMR(Soln_ptr,
			InputParameters,
   	                QuadTree,
 	                GlobalSolnBlockList,
	                LocalSolnBlockList,
	                OFF,ON);
       if (error_flag) return (error_flag);

       ICs(Soln_ptr,
           LocalSolnBlockList,
           InputParameters);

       BCs(Soln_ptr,
           LocalSolnBlockList,
	   InputParameters);

       error_flag = Send_All_Messages(Soln_ptr,
                                      LocalSolnBlockList,
                                      Soln_ptr[0].NumVar(),
                                      OFF);
       if (error_flag) return (error_flag);

       if (CFFC_Primary_MPI_Processor()) {
          cout << "\n Refinement Level #" << number_of_uniform_mesh_refinements
               << " : Number of Blocks = " << QuadTree.countUsedBlocks()
               << ", Number of Cells = " << QuadTree.countUsedCells()
               << ", Refinement Efficiency = " << QuadTree.efficiencyRefinement();
       } /* endif */

    } /* endfor */

    return(error_flag);

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
template <class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int Boundary_AMR(Quad_Soln_Block              *Soln_ptr,
		 Quad_Soln_Input_Parameters   &InputParameters,
		 QuadTreeBlock_DataStructure  &QuadTree,
		 AdaptiveBlockResourceList    &GlobalSolnBlockList,
		 AdaptiveBlock2D_List         &LocalSolnBlockList) {

    int i, j, nb, error_flag, number_of_boundary_mesh_refinements;

    error_flag = 0;

    if (InputParameters.Number_of_Boundary_Mesh_Refinements == 0) return(error_flag);

    for (number_of_boundary_mesh_refinements = 1; 
         number_of_boundary_mesh_refinements <= InputParameters.Number_of_Boundary_Mesh_Refinements; 
         ++number_of_boundary_mesh_refinements) {

      // Set refinement flags.
      QuadTree.nochangeAll();
      for (int nb = 0; nb < LocalSolnBlockList.Nblk; nb++) {
	LocalSolnBlockList.RefineFlag[nb] = ADAPTIVEBLOCK2D_NOCHANGE;
	if (LocalSolnBlockList.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
	  for (int i = Soln_ptr[nb].ICl-Soln_ptr[nb].Nghost; i <= Soln_ptr[nb].ICu+Soln_ptr[nb].Nghost; i++) {
	    if ((Soln_ptr[nb].Flow_Type == FLOWTYPE_INVISCID &&
		 (Soln_ptr[nb].Grid.BCtypeN[i] == BC_REFLECTION ||
		  Soln_ptr[nb].Grid.BCtypeS[i] == BC_REFLECTION)) ||
		(Soln_ptr[nb].Flow_Type != FLOWTYPE_INVISCID &&
		 (Soln_ptr[nb].Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX ||
		  Soln_ptr[nb].Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
		  Soln_ptr[nb].Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX ||
		  Soln_ptr[nb].Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL)) ||
		(Soln_ptr[nb].Grid.BCtypeN[i] == BC_BURNING_SURFACE ||
		 Soln_ptr[nb].Grid.BCtypeS[i] == BC_BURNING_SURFACE)) {
	      if (LocalSolnBlockList.Block[nb].info.level < InputParameters.Maximum_Refinement_Level) {
		LocalSolnBlockList.RefineFlag[nb] = ADAPTIVEBLOCK2D_REFINE;
	      }
	    }
	  }
	  for (int j = Soln_ptr[nb].JCl-Soln_ptr[nb].Nghost; j <= Soln_ptr[nb].JCu+Soln_ptr[nb].Nghost; j++) {
	    if ((Soln_ptr[nb].Flow_Type == FLOWTYPE_INVISCID &&
		 (Soln_ptr[nb].Grid.BCtypeE[j] == BC_REFLECTION ||
		  Soln_ptr[nb].Grid.BCtypeW[j] == BC_REFLECTION)) ||
		(Soln_ptr[nb].Flow_Type != FLOWTYPE_INVISCID &&
		 (Soln_ptr[nb].Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX ||
		  Soln_ptr[nb].Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
		  Soln_ptr[nb].Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX ||
		  Soln_ptr[nb].Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL)) ||
		(Soln_ptr[nb].Grid.BCtypeE[j] == BC_BURNING_SURFACE ||
		 Soln_ptr[nb].Grid.BCtypeW[j] == BC_BURNING_SURFACE)) {
	      if (LocalSolnBlockList.Block[nb].info.level < InputParameters.Maximum_Refinement_Level) {
		LocalSolnBlockList.RefineFlag[nb] = ADAPTIVEBLOCK2D_REFINE;
	      }
	    }
	  }
	}
      }

       error_flag = AMR(Soln_ptr,
			InputParameters,
   	                QuadTree,
 	                GlobalSolnBlockList,
	                LocalSolnBlockList,
	                OFF,ON);
       if (error_flag) return (error_flag);

       ICs(Soln_ptr,
           LocalSolnBlockList,
           InputParameters);

       BCs(Soln_ptr,
           LocalSolnBlockList,
	   InputParameters);

       error_flag = Send_All_Messages(Soln_ptr,
                                      LocalSolnBlockList,
                                      Soln_ptr[0].NumVar(),
                                      OFF);
       if (error_flag) return (error_flag);

       if (CFFC_Primary_MPI_Processor()) {
          cout << "\n Refinement Level #" << number_of_boundary_mesh_refinements
               << " : Number of Blocks = " << QuadTree.countUsedBlocks()
               << ", Number of Cells = " << QuadTree.countUsedCells()
               << ", Refinement Efficiency = " << QuadTree.efficiencyRefinement();
       } /* endif */

    } /* endfor */

    return(error_flag);

}

#endif /* _AMR_INCLUDED  */

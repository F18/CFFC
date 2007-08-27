/**********************************************************************
 * EmbeddedBoundaries2D_FASMultigrid.h: Header file declaring the     *
 *                                      inherited FAS Multigrid 2D    *
 *                                      class for embedded boundary   *
 *                                      calculations.                 *
 **********************************************************************/

#ifndef _FASMULTIGRID2D_EMBEDDEDBOUNDARIES_INCLUDED
#define _FASMULTIGRID2D_EMBEDDEDBOUNDARIES_INCLUDED

// Include required C++ libraries.

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <cstring>

#ifdef _GNU_GCC_V3
#include <limits>
#endif

using namespace std;

// Include the embedded boundaries header file.

#ifndef _EMBEDDEDBOUNDARIES2D_INCLUDED
#include "EmbeddedBoundaries2D.h"
#endif // _EMBEDDEDBOUNDARIES2D_INCLUDED

// Include multigrid input header file.

#ifndef _FASMULTIGRID2D_INCLUDED
#include "../FASMultigrid2D/FASMultigrid2D.h"
#endif // _FASMULTIGRID2D_INCLUDED

// #define _MG_DOUT_CYCLE_STAGE_

/*!
 * class: FAS_Multigrid2D_Solver
 *
 * @brief Full Approximation Storage (FAS) multigrid solution class.
 *
 * This class conducts the time integration of a hyperbolic/elliptic
 * system of partial differential equations using a Full Approximation
 * Storage (FAS) multigrid solver.  Spatial discretization is
 * accomplised with a higher-order Godunov-type finite volume scheme on
 * a body-fitted multi-block quadrilateral mesh.  Unsteady calculations
 * are facilitated through a dual-time-stepping procedure.
 *
 */
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
class EB_FAS_Multigrid2D_Solver {
public:

  //! Pointer to the input parameters.
  Quad_Soln_Input_Parameters *IP;

  //! Pointer to the quadtree.
  QuadTreeBlock_DataStructure *QuadTree;

  //! Pointer to the global solution block list.
  AdaptiveBlockResourceList *List_of_Global_Solution_Blocks;

  //! Array of the local solution block list (one for each grid level).
  AdaptiveBlock2D_List *List_of_Local_Solution_Blocks;

  //! Array of local solution block information (first index corresponds
  //! to the grid level and the second index corresponds to the solution
  //! blocks).
  Quad_Soln_Block **Local_SolnBlk;

  //! Array of multigrid solution block information (first index
  //! corresponds to the grid level and the second index corresponds to
  //! the solution blocks).
  FAS_Multigrid_Quad_Block<cState> **MG_SolnBlk;

  //! Array of dual-time-stepping multigrid solution block information
  //! (first index corresponds to the grid level and the second index
  //! corresponds to the solution blocks).
  DTS_Multigrid_Quad_Block<cState> **DTS_SolnBlk;

  //! Embedded boundary solver declaration for every grid level of the
  //! multigrid.
  EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters> *EBSolver;

#ifdef _MG_DOUT_CYCLE_STAGE_
  char extension[256], output_file_name[256];
  char *output_file_name_ptr;
  ofstream dout;
#endif

  //! Default Constructor.
  EB_FAS_Multigrid2D_Solver(void) {
    IP = NULL;
    QuadTree = NULL;
    List_of_Global_Solution_Blocks = NULL;
    List_of_Local_Solution_Blocks = NULL;
    Local_SolnBlk = NULL;
    MG_SolnBlk = NULL;
    DTS_SolnBlk = NULL;
    EBSolver = NULL;
  }

  //! Memory allocation and initialization for all grid levels.
  int allocate(Quad_Soln_Block *FinestBlks,
	       QuadTreeBlock_DataStructure *FinestQuadTree,
	       AdaptiveBlockResourceList *FinestGlobalList,
	       AdaptiveBlock2D_List *FinestLocalList,
	       Quad_Soln_Input_Parameters *ip,
	       LevelSet2D_Quad_Block *LS_SolnBlk,
	       QuadTreeBlock_DataStructure *LS_QuadTree,
	       AdaptiveBlockResourceList *LS_GlobalList,
	       AdaptiveBlock2D_List *LS_LocalList,
	       LevelSet2D_Input_Parameters *LS_ip,
	       EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters> *FineGrid_EBSolver);

  //! Memory deallocation.
  void deallocate(void);

  //! Output the node state information for each grid level.
  int Output_Multigrid(const int &number_of_time_steps,
 		       const double &Time,
		       bool writing_intermediate_soln = false, 
		       double l2_norm = -1.0, 
                       double l2_norm_rel = -1.0);

  //! Output the cell-centred state information for each grid level.
  int Output_Multigrid_Cells(const int &number_of_time_steps,
 			     const double &Time, 
			     bool writing_intermediate_soln = false, 
			     double l2_norm = -1.0, 
                             double l2_norm_rel = -1.0);

  //! Output the element state information for each grid level.
  int Output_Multigrid_Elements(const int &number_of_time_steps,
				const double &Time);

  //! Output the node information for each grid level.
  int Output_Multigrid_Nodes(const int &number_of_time_steps,
 			     const double &Time);

  //! Output the cell status information for each grid level.
  int Output_Multigrid_Cell_Status(const int &number_of_time_steps,
				   const double &Time);

  //! Restrict the solution block a from fine to a coarse grid.
  void Restrict_Solution_Blocks(const int &Level_Fine);

  //! Restrict the residual from a fine to coarse a grid.
  void Restrict_Residuals(const int &Level_Fine);

  //! Restrict the boundary reference states a fine to coarse a grid.
  void Restrict_Boundary_Ref_States(const int &Level_Fine);

  //! Prolong the solution block from a coarse to a fine grid.
  int Prolong_Solution_Blocks(const int &Level_Coarse);

  //! Prolong and update the solution blocks from a fine to a coarse grid.
  int Prolong_and_Update_Solution_Blocks(const int &Level_Coarse);

  //! Exchange the solution information between solution blocks at the
  //! specified mesh level.					 
  int Exchange_Solution_Information(const int &Level,
				    const int &Send_Mesh_Geometry_Only);

  //! Store the solution at the current grid level in the uo_MG array.
  void Store_Current_Solution_in_uo_MG(const int &Level);

  //! Subtract the residual from the forcing term to finalize the 
  //! forcing term.
  void Subtract_dUdt_from_P(const int &Level);

  //! Determine the minimum time-step on the coarse grid level.
  void CFL_Multigrid(const int &Level_Coarse);

  //! Update the primitive solution states on the specified grid level.
  void Update_Primitive_Variables(const int &Level);

  //! Evaluate the solution chanes on the specified grid level.
  void Evaluate_Solution_Changes(const int &Level);

  //! Evaluate the stage solution residual for the specified grid level.
  int dUdt_Multistage_Explicit_for_Multigrid(const int &Level,
					     const int &Top_Level,
					     const int &i_stage,
					     const double &Time);

  //! Apply the flux corrections at boundaries.
  void Apply_Boundary_Flux_Corrections_Multistage_Explicit_for_Multigrid(const int &Level,
									 const int &Top_Level,
									 const int &i_stage);

  //! Update the solution blocks at the specified grid level.
  int Update_Solution_Multistage_Explicit_for_Multigrid(const int &Level,
							const int &Top_Level,
 							const int &i_stage);

  //! Apply the defect correction forcing term, P.
  int Apply_the_FAS_Multigrid_Forcing_Term(const int &Level,
 					   const int &Top_Level,
 					   const int &i_stage);

  //! Evaluate the solution residual at the specified grid level.
  int Residual_Evaluation_for_Multigrid(const int &Level,
					const int &Top_Level,
					const double &dt,
					const double &Time,
					const int &apply_forcing_term_flag);

  //! Perform the multigrid smoothing at the specified grid level.
  int Smooth(const int &Level,
	     const int &Top_Level,
	     const double &dt,
	     const double &Time);

  //! Conduct the FAS multigrid algorithm for the specified series
  //! of grids.
  int Coarse_Grid_Correction(const int &Top_Level,
 			     const int &Current_Level,
 			     const double &dt,
			     const double &Time);

  //! Execute the FAS multigrid algorithm for steady state computations.
  int Execute(int &batch_flag,
	      int &number_of_time_steps,
	      int &evolution_counter,
	      int &levelset_iterations,
	      double &Time,
	      double &levelset_Time,
	      CPUTime &processor_cpu_time,
	      CPUTime &total_cpu_time,
	      ofstream &residual_file);

  ///////////////////////////////////////////////
  // Routines required for dual-time-stepping. //
  ///////////////////////////////////////////////

  //! Memory reallocation and reinitialization for the coarse grid 
  //! levels after refinement.
  int reallocate(void);

  //! Store the previous solution and restrict to coarse levels.
  int Store_Previous_Solution(void);

  //! Determine the solution residual for stage i of an n-stage scheme.
  int dUdtau_Multistage_Explicit(const int &Level,
				 const int &i_stage,
				 const double &dt);

  //! Apply the physical-time-derivative implicit correction to the
  //! time step.
  int Apply_Melson_Time_Step(const int &Level,
			     const int &i_stage,
			     const double &dt);

  //! Compute the location of the interface and readjust grids.
  int Compute_Interface_Location(const int &batch_flag,
				 const double &current_time,
				 const double &maximum_time,
				 int &evolution_counter,
				 int &levelset_iterations,
				 double &level_set_current_time);

  //! Perform the DTS multigrid solution.
  int DTS_Multigrid_Solution(int &batch_flag,
			     int &number_of_time_steps,
			     int &evolution_counter,
			     int &levelset_iterations,
			     double &Time,
			     double &levelset_Time,
			     CPUTime &processor_cpu_time,
			     CPUTime &total_cpu_time,
			     ofstream &residual_file);

};

/**********************************************************************
 * EB_FAS_Multigrid2D_Solver::allocate --                             *
 *                                                                    *
 * This routine performs the memory allocation and initialization for *
 * all grid levels of the FAS multigrid solution class.               *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState,
	  class Quad_Soln_Block,
	  class Quad_Soln_Input_Parameters>
int EB_FAS_Multigrid2D_Solver<cState, pState,
			      Quad_Soln_Block,
			      Quad_Soln_Input_Parameters>::
allocate(Quad_Soln_Block *FinestBlks,
	 QuadTreeBlock_DataStructure *FinestQuadTree,
	 AdaptiveBlockResourceList *FinestGlobalList,
	 AdaptiveBlock2D_List *FinestLocalList,
	 Quad_Soln_Input_Parameters *ip,
	 LevelSet2D_Quad_Block *LS_SolnBlk,
	 QuadTreeBlock_DataStructure *LS_QuadTree,
	 AdaptiveBlockResourceList *LS_GlobalList,
	 AdaptiveBlock2D_List *LS_LocalList,
	 LevelSet2D_Input_Parameters *LS_ip,
	 EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters> *FineGrid_EBSolver) {

  int error_flag;

#ifdef _MG_DOUT_CYCLE_STAGE_
  strcpy(output_file_name,"mg");
  strcat(output_file_name,"_cpu");
  sprintf(extension,"%.6d",FinestLocalList->ThisCPU);
  strcat(extension,".txt");
  strcat(output_file_name,extension);
  output_file_name_ptr = output_file_name;
  dout.open(output_file_name_ptr,ios::out);
  if (dout.bad()) return 1;
#endif

  // Point the input parameters to the given input parameters.
  IP = ip;

  // Point the quadtree to the input quadtree.
  QuadTree = FinestQuadTree;

  // Point the global solution block list to the input list.
  List_of_Global_Solution_Blocks = FinestGlobalList;

  // Create list of local solution blocks for each multigrid level.
  List_of_Local_Solution_Blocks = new AdaptiveBlock2D_List[IP->Multigrid_IP.Levels];

  // Point the list of local solution blocks for the finest level to
  // the input list of local solution blocks.
  List_of_Local_Solution_Blocks[FINEST_LEVEL] = *FinestLocalList;
  // Allocate memory for the coarse grid list of local solution blocks
  // and set the CPU number.
  for (int level = 1; level < IP->Multigrid_IP.Levels; level++) {
    List_of_Local_Solution_Blocks[level].allocate(IP->Number_of_Blocks_Per_Processor);
    List_of_Local_Solution_Blocks[level].ThisCPU = List_of_Local_Solution_Blocks[FINEST_LEVEL].ThisCPU;
  }

  // Allocate memory for the local solution blocks on each level.
  Local_SolnBlk = new Quad_Soln_Block*[IP->Multigrid_IP.Levels];
  // Point the local solution block for the finest level to the input
  // local solution block.
  Local_SolnBlk[FINEST_LEVEL] = FinestBlks;
  // Allocate memory for the coarse grid local solution blocks.
  for (int level = 1; level < IP->Multigrid_IP.Levels; level++) {
    Local_SolnBlk[level] = new Quad_Soln_Block[IP->Number_of_Blocks_Per_Processor];
  }

  // Allocate memory for the FAS multigrid solution blocks on each level.
  MG_SolnBlk = new FAS_Multigrid_Quad_Block<cState>*[IP->Multigrid_IP.Levels];
  // Allocate memory for the coarse grid FAS multigrid solution blocks.
  for (int level = 0; level < IP->Multigrid_IP.Levels; level++) {
    MG_SolnBlk[level] = new FAS_Multigrid_Quad_Block<cState>[IP->Number_of_Blocks_Per_Processor];
  }

  // Allocate memory for the DTS multigrid solution blocks on each level.
  if (IP->i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING) {
    DTS_SolnBlk = new DTS_Multigrid_Quad_Block<cState>*[IP->Multigrid_IP.Levels];
    // Allocate memory for the coarse grid DTS multigrid solution blocks.
    for (int level = 0; level < IP->Multigrid_IP.Levels; level++) {
      DTS_SolnBlk[level] = new DTS_Multigrid_Quad_Block<cState>[IP->Number_of_Blocks_Per_Processor];
    }
  }

  // Allocate memory for the embedded boundary blocks on each level.
  EBSolver = new EmbeddedBoundaries2D<cState,pState,Quad_Soln_Block,Quad_Soln_Input_Parameters>[IP->Multigrid_IP.Levels];
  // Point the local solution block for the finest level to the input
  // local solution block.
  EBSolver[FINEST_LEVEL].allocate(FinestBlks,
				  FinestQuadTree,
				  FinestGlobalList,
				  FinestLocalList,
				  ip,
				  LS_SolnBlk,
				  LS_QuadTree,
				  LS_GlobalList,
				  LS_LocalList,
				  LS_ip);
  EBSolver[FINEST_LEVEL].Initialize_Adjustment_Grids();
  EBSolver[FINEST_LEVEL].Interface_Component_List.Copy(FineGrid_EBSolver->Interface_Component_List);
  EBSolver[FINEST_LEVEL].Interface_Union_List.Copy(FineGrid_EBSolver->Interface_Union_List);
  for (int nb = 0; nb < IP->Number_of_Blocks_Per_Processor; nb++) {
    if (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      Copy_Quad_Block(EBSolver[FINEST_LEVEL].AGrid[nb],FineGrid_EBSolver->AGrid[nb]);
      Copy_Quad_Block(EBSolver[FINEST_LEVEL].OGrid[nb],FineGrid_EBSolver->OGrid[nb]);
      EBSolver[FINEST_LEVEL].Mesh[nb].Copy_Quad_Block(FineGrid_EBSolver->Mesh[nb]);
      EBSolver[FINEST_LEVEL].AMesh[nb].Copy_Quad_Block(FineGrid_EBSolver->AMesh[nb]);
      EBSolver[FINEST_LEVEL].OMesh[nb].Copy_Quad_Block(FineGrid_EBSolver->OMesh[nb]);
      EBSolver[FINEST_LEVEL].Adjustment_Data[nb].Copy_Quad_Block(FineGrid_EBSolver->Adjustment_Data[nb]);
    }
  }

  // Allocate memory and set data for all coarse mesh variables on all
  // blocks on this processor.
  for (int nb = 0; nb < IP->Number_of_Blocks_Per_Processor; nb++) {
    if (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].used == ADAPTIVEBLOCK2D_USED) {

      // Ensure that the number of cells in each direction is even.
      if (Local_SolnBlk[FINEST_LEVEL][nb].NCi % 2 != 0) {
	cout << "\nFASMultigrid2D Error: block #" << nb << " on processor " 
             << List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].info.cpu 
             << ", on level 0 has " << Local_SolnBlk[FINEST_LEVEL][nb].NCi 
             << "cells in the x-direction, which is odd; cannot coarsen any further." << endl;
	return 1101;
      }
      if (Local_SolnBlk[FINEST_LEVEL][nb].NCj % 2 != 0) {
	cout << "\nFASMultigrid2D Error: block #" << nb << " on processor " 
             << List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].info.cpu 
             << ", on level 0 has " << Local_SolnBlk[FINEST_LEVEL][nb].NCj 
             << " cells in the y-direction, which is odd; cannot coarsen any further." << endl;
	return 1102;
      }

      // Allocate memory for the forcing term and the uo storage on the
      // finest level.
      MG_SolnBlk[FINEST_LEVEL][nb].allocate(Local_SolnBlk[FINEST_LEVEL][nb].NCi,
					    Local_SolnBlk[FINEST_LEVEL][nb].NCj);

      // Allocate memory for the DTS solution block finest level.
      if (IP->i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING) {
	DTS_SolnBlk[FINEST_LEVEL][nb].allocate(Local_SolnBlk[FINEST_LEVEL][nb].NCi,
					       Local_SolnBlk[FINEST_LEVEL][nb].NCj);
      }

      // Allocate memory and set data for the coarse mesh levels.
      for (int level = 1; level < IP->Multigrid_IP.Levels; level++) {

	// Copy the list of local solution block neighbour information
	// and calculate the coarse grid resolution.
	List_of_Local_Solution_Blocks[level].Block[nb] =
	  List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb];
	List_of_Local_Solution_Blocks[level].Block[nb].info.dimen.i = 
	  int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].info.dimen.i/pow(2.0,double(level)));
	List_of_Local_Solution_Blocks[level].Block[nb].info.dimen.j = 
	  int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].info.dimen.j/pow(2.0,double(level)));
	
	for (int n = 0; n < List_of_Local_Solution_Blocks[level].Block[nb].nS; n++) {
	  List_of_Local_Solution_Blocks[level].Block[nb].infoS[n].dimen.i = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoS[n].dimen.i/pow(2.0,double(level)));
	  List_of_Local_Solution_Blocks[level].Block[nb].infoS[n].dimen.j = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoS[n].dimen.j/pow(2.0,double(level)));
	}
	for (int n = 0; n < List_of_Local_Solution_Blocks[level].Block[nb].nN; n++) {
	  List_of_Local_Solution_Blocks[level].Block[nb].infoN[n].dimen.i = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoN[n].dimen.i/pow(2.0,double(level)));
	  List_of_Local_Solution_Blocks[level].Block[nb].infoN[n].dimen.j = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoN[n].dimen.j/pow(2.0,double(level)));
	}
	for (int n = 0; n < List_of_Local_Solution_Blocks[level].Block[nb].nE; n++) {
	  List_of_Local_Solution_Blocks[level].Block[nb].infoE[n].dimen.i = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoE[n].dimen.i/pow(2.0,double(level)));
	  List_of_Local_Solution_Blocks[level].Block[nb].infoE[n].dimen.j = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoE[n].dimen.j/pow(2.0,double(level)));
	}
	for (int n = 0; n < List_of_Local_Solution_Blocks[level].Block[nb].nW; n++) {
	  List_of_Local_Solution_Blocks[level].Block[nb].infoW[n].dimen.i = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoW[n].dimen.i/pow(2.0,double(level)));
	  List_of_Local_Solution_Blocks[level].Block[nb].infoW[n].dimen.j = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoW[n].dimen.j/pow(2.0,double(level)));
	}
	for (int n = 0; n < List_of_Local_Solution_Blocks[level].Block[nb].nSE; n++) {
	  List_of_Local_Solution_Blocks[level].Block[nb].infoSE[n].dimen.i = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoSE[n].dimen.i/pow(2.0,double(level)));
	  List_of_Local_Solution_Blocks[level].Block[nb].infoSE[n].dimen.j = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoSE[n].dimen.j/pow(2.0,double(level)));
	}
	for (int n = 0; n < List_of_Local_Solution_Blocks[level].Block[nb].nSW; n++) {
	  List_of_Local_Solution_Blocks[level].Block[nb].infoSW[n].dimen.i = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoSW[n].dimen.i/pow(2.0,double(level)));
	  List_of_Local_Solution_Blocks[level].Block[nb].infoSW[n].dimen.j = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoSW[n].dimen.j/pow(2.0,double(level)));
	}
	for (int n = 0; n < List_of_Local_Solution_Blocks[level].Block[nb].nNE; n++) {
	  List_of_Local_Solution_Blocks[level].Block[nb].infoNE[n].dimen.i = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoNE[n].dimen.i/pow(2.0,double(level)));
	  List_of_Local_Solution_Blocks[level].Block[nb].infoNE[n].dimen.j = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoNE[n].dimen.j/pow(2.0,double(level)));
	}
	for (int n = 0; n < List_of_Local_Solution_Blocks[level].Block[nb].nNW; n++) {
	  List_of_Local_Solution_Blocks[level].Block[nb].infoNW[n].dimen.i = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoNW[n].dimen.i/pow(2.0,double(level)));
	  List_of_Local_Solution_Blocks[level].Block[nb].infoNW[n].dimen.j = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoNW[n].dimen.j/pow(2.0,double(level)));
	}

	// Set-up the local block list for each level.
	List_of_Local_Solution_Blocks[level].Block[nb].used = ADAPTIVEBLOCK2D_USED;
	List_of_Local_Solution_Blocks[level].Block[nb].gblknum = 
	  List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].gblknum;
	List_of_Local_Solution_Blocks[level].Block[nb].info.cpu = 
	  List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].info.cpu;
	List_of_Local_Solution_Blocks[level].Block[nb].info.blknum = 
	  List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].info.blknum;
	List_of_Local_Solution_Blocks[level].Block[nb].info.dimen.i = 
	  int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].info.dimen.i/pow(2.0,double(level)));
	List_of_Local_Solution_Blocks[level].Block[nb].info.dimen.j = 
	  int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].info.dimen.j/pow(2.0,double(level)));

	// If the current level is not the coarsest level then ensure
	// that the number of cells is even.
	if (level != IP->Multigrid_IP.Levels-1) {
	  if (List_of_Local_Solution_Blocks[level].Block[nb].info.dimen.i % 2 != 0) {
	    cout << "\nFASMultigrid2D Error: block #" << nb << " on processor " 
                 << List_of_Local_Solution_Blocks[level].Block[nb].info.cpu 
                 << ", on level " << level << " has " 
                 << List_of_Local_Solution_Blocks[level].Block[nb].info.dimen.i 
                 << " cells in the x-direction, which is odd; cannot coarsen any further." << endl;
	    return 1103;
	  }
	}
	if (level != IP->Multigrid_IP.Levels-1) {
	  if (List_of_Local_Solution_Blocks[level].Block[nb].info.dimen.j % 2 != 0) {
	    cout << "\nFASMultigrid2D Error: block #" << nb << " on processor " 
                 << List_of_Local_Solution_Blocks[level].Block[nb].info.cpu << ", on level " 
                 << level << " has " << List_of_Local_Solution_Blocks[level].Block[nb].info.dimen.j 
                 << " cells in the y-direction, which is odd; cannot coarsen any further." << endl;
	    return 1104;
	  }
	}
	
	List_of_Local_Solution_Blocks[level].Block[nb].info.dimen.ghost = 
	  List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].info.dimen.ghost;
	List_of_Local_Solution_Blocks[level].Block[nb].info.sector = 
	  List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].info.sector;
	List_of_Local_Solution_Blocks[level].Block[nb].info.level = 
	  List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].info.level;
	
	// Allocate the coarse grid local solution block.
	Local_SolnBlk[level][nb].allocate(List_of_Local_Solution_Blocks[level].Block[nb].info.dimen.i,
					  List_of_Local_Solution_Blocks[level].Block[nb].info.dimen.j,
					  List_of_Local_Solution_Blocks[level].Block[nb].info.dimen.ghost);

	// Create the coarse grid mesh.
	if (level == 1) Half_Mesh_Resolution(Local_SolnBlk[level][nb].Grid,
					     FineGrid_EBSolver->OGrid[nb]);
	else Half_Mesh_Resolution(Local_SolnBlk[level][nb].Grid,
				  Local_SolnBlk[level-1][nb].Grid);
	EBSolver[level].allocate(Local_SolnBlk[level],
				 FinestQuadTree,
				 FinestGlobalList,
				 &List_of_Local_Solution_Blocks[level],
				 ip,
				 LS_SolnBlk,
				 LS_QuadTree,
				 LS_GlobalList,
				 LS_LocalList,
				 LS_ip);
	EBSolver[level].Interface_Component_List.Copy(FineGrid_EBSolver->Interface_Component_List);
	EBSolver[level].Interface_Union_List.Copy(FineGrid_EBSolver->Interface_Union_List);
	EBSolver[level].Initialize_Adjustment_Grids();

	// Allocate the coarse grid FAS multigrid solution block.
	MG_SolnBlk[level][nb].allocate(Local_SolnBlk[level][nb].NCi,
				       Local_SolnBlk[level][nb].NCj);

	// Allocate the coarse grid DTS multigrid solution block.
	if (IP->i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING) {
	  DTS_SolnBlk[level][nb].allocate(Local_SolnBlk[level][nb].NCi,
					  Local_SolnBlk[level][nb].NCj);
	}

	// Allocate memory for the message passing buffers used to send solution
	// information between neighbouring blocks for the coarse grid.
	Allocate_Message_Buffers(List_of_Local_Solution_Blocks[level],
				 Local_SolnBlk[level][nb].NumVar()+NUM_COMP_VECTOR2D);

      }

    }
  }

  // Solution block allocation and assignment successful.
  return 0;

}

/**********************************************************************
 * EB_FAS_Multigrid2D_Solver::deallocate --                           *
 *                                                                    *
 * This routine performs the memory deallocation for the FAS          *
 * multigrid solution class.                                          *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState,
	  class Quad_Soln_Block,
	  class Quad_Soln_Input_Parameters>
void EB_FAS_Multigrid2D_Solver<cState, pState,
			       Quad_Soln_Block,
			       Quad_Soln_Input_Parameters>::
deallocate(void) {

  // Deallocate the embedded boundary solution array.
  if (EBSolver != NULL) { delete []EBSolver; EBSolver = NULL; }

  // Deallocate the local solution blocks.
  for (int level = 1; level < IP->Multigrid_IP.Levels; level++) {    
    for (int nb = 0 ; nb < IP->Number_of_Blocks_Per_Processor; nb++) {
      if (Local_SolnBlk[level][nb].U != NULL) Local_SolnBlk[level][nb].deallocate();
    }
    if (Local_SolnBlk[level] != NULL) delete []Local_SolnBlk[level];
  }
  if (Local_SolnBlk != NULL) { delete []Local_SolnBlk; Local_SolnBlk = NULL; }

  // Deallocate the FAS multigrid solution blocks.
  for (int level = 0; level < IP->Multigrid_IP.Levels; level++) {    
    if (MG_SolnBlk[level] != NULL) delete []MG_SolnBlk[level];
  }
  if (MG_SolnBlk != NULL) { delete []MG_SolnBlk; MG_SolnBlk = NULL; }

  // Deallocate the DTS multigrid solution blocks if required.
  if (IP->i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING) {
    for (int level = 0; level < IP->Multigrid_IP.Levels; level++) {
      if (DTS_SolnBlk[level] != NULL) delete []DTS_SolnBlk[level];
    }
    if (DTS_SolnBlk != NULL) { delete []DTS_SolnBlk; DTS_SolnBlk = NULL; }
  }

  // Deallocate the local lists of solution blocks.  
  if (List_of_Local_Solution_Blocks != NULL) {
    for (int level = 1; level < IP->Multigrid_IP.Levels; level++) {
      if (List_of_Local_Solution_Blocks[level].Block != NULL) {
	//Deallocate_Message_Buffers(List_of_Local_Solution_Blocks[level]);
	List_of_Local_Solution_Blocks[level].deallocate();
      }
    }
    delete []List_of_Local_Solution_Blocks; List_of_Local_Solution_Blocks = NULL;
  }

  // Point the quadtree and input parameter pointers to null.
  QuadTree = NULL;
  IP = NULL;

#ifdef _MG_DOUT_CYCLE_STAGE_
  dout.close();
#endif

}

/**********************************************************************
 * EB_FAS_Multigrid2D_Solver::Output_Multigrid  --                    *
 *                                                                    *
 * This routine writes out the node information for each multigrid    *
 * grid level (finest to coarest) for a 1D array of 2D quadrilateral  *
 * multi-block solution blocks.                                       *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState,
	  class Quad_Soln_Block,
	  class Quad_Soln_Input_Parameters>
int EB_FAS_Multigrid2D_Solver<cState, pState,
			      Quad_Soln_Block,
			      Quad_Soln_Input_Parameters>::
Output_Multigrid(const int &number_of_time_steps,
		 const double &Time,
		 bool writing_intermediate_soln,
		 double l2_norm, 
                 double l2_norm_rel) {

  int i, i_output_title;
  char prefix[256], extension[256], output_file_name[256];
  char *output_file_name_ptr;
  ofstream output_file;
  int max_level = (writing_intermediate_soln ? 1 : IP->Multigrid_IP.Levels);

  // Determine main prefix of output data file names.
  i = 0;
  while (1) {
    if (IP->Output_File_Name[i] == ' ' ||
	IP->Output_File_Name[i] == '.') break;
    prefix[i] = IP->Output_File_Name[i];
    i = i + 1;
    if (i > strlen(IP->Output_File_Name)) break;
  }
  prefix[i] = '\0';
  strcat(prefix,"_");

  // Output to a seperate file for each multigrid level.
  for (int level = 0; level < IP->Multigrid_IP.Levels; level++) {

    // Determine prefix of output data file names for current grid level.
    strcpy(output_file_name,prefix);
    sprintf(extension,"%.2d",level);
    strcat(output_file_name,extension);
    if (writing_intermediate_soln) {
       sprintf(extension,"_n1%.4d",number_of_time_steps);
       strcat(output_file_name,extension);
    }
    strcat(output_file_name,"_cpu");

    // Determine output data file name for this processor.
    sprintf(extension,"%.6d",List_of_Local_Solution_Blocks[level].ThisCPU);
    strcat(extension,".dat");
    strcat(output_file_name,extension);
    output_file_name_ptr = output_file_name;

    // Open the output data file.
    output_file.open(output_file_name_ptr,ios::out);
    if (output_file.bad()) return 1;

    // Write the solution data for each solution block.
    i_output_title = 1;
    for (int nb = 0; nb < List_of_Local_Solution_Blocks[level].Nblk; nb++) {
      if (List_of_Local_Solution_Blocks[level].Block[nb].used == ADAPTIVEBLOCK2D_USED) {
	Output_Tecplot(Local_SolnBlk[level][nb],
 		       *IP,
 		       number_of_time_steps,
 		       Time,
 		       List_of_Local_Solution_Blocks[level].Block[nb].gblknum,
 		       i_output_title,
 		       output_file);
 	if (i_output_title) i_output_title = 0;
      }
    }

    // Close the output data file.
    output_file.close();

  }

  // Writing of multigrid grids data files complete.
  return 0;

}

/**********************************************************************
 * EB_FAS_Multigrid2D_Solver::Output_Multigrid_Cells --               *
 *                                                                    *
 * This routine writes out the cell-centred information for each      *
 * multigrid grid level (finest to coarest) for a 1D array of 2D      *
 * quadrilateral multi-block solution blocks.                         *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState,
	  class Quad_Soln_Block,
	  class Quad_Soln_Input_Parameters>
int EB_FAS_Multigrid2D_Solver<cState, pState,
			      Quad_Soln_Block,
			      Quad_Soln_Input_Parameters>::
Output_Multigrid_Cells(const int &number_of_time_steps,
		       const double &Time, 
		       bool writing_intermediate_soln,
		       double l2_norm, 
                       double l2_norm_rel) {

  int i, i_output_title;
  char prefix[256], extension[256], output_file_name[256];
  char *output_file_name_ptr;
  ofstream output_file;
  int max_level = (writing_intermediate_soln ? 1 : IP->Multigrid_IP.Levels);

  // Determine main prefix of output data file names.
  i = 0;
  while (1) {
    if (IP->Output_File_Name[i] == ' ' ||
	IP->Output_File_Name[i] == '.') break;
    prefix[i] = IP->Output_File_Name[i];
    i = i + 1;
    if (i > strlen(IP->Output_File_Name)) break;
  }
  prefix[i] = '\0';
  strcat(prefix,"_cells_");

  // Output to a seperate file for each multigrid level.
  for (int level = 0; level < IP->Multigrid_IP.Levels; level++) {

    // Determine prefix of output data file names for current grid level.
    strcpy(output_file_name,prefix);
    sprintf(extension,"_%.2d",level);
    strcat(output_file_name,extension);
    if (writing_intermediate_soln) {
       sprintf(extension,"_n1%.4d",number_of_time_steps);
       strcat(output_file_name,extension);
    }
    strcat(output_file_name,"_cpu");

    // Determine output data file name for this processor.
    sprintf(extension,"%.6d",List_of_Local_Solution_Blocks[level].ThisCPU);
    strcat(extension,".dat");
    strcat(output_file_name,extension);
    output_file_name_ptr = output_file_name;

    // Open the output data file.
    output_file.open(output_file_name_ptr,ios::out);
    if (output_file.bad()) return 1;

    // Write the solution data for each solution block.
    i_output_title = 1;
    for (int nb = 0; nb < List_of_Local_Solution_Blocks[level].Nblk; nb++) {
      if (List_of_Local_Solution_Blocks[level].Block[nb].used == ADAPTIVEBLOCK2D_USED) {
	Output_Cells_Tecplot(Local_SolnBlk[level][nb],
                             *IP,
			     number_of_time_steps,
			     Time,
			     List_of_Local_Solution_Blocks[level].Block[nb].gblknum,
			     i_output_title,
			     output_file);
	if (i_output_title) i_output_title = 0;
      }
    }

    // Close the output data file.
    output_file.close();

  }

  // Writing of multigrid grids data files complete.
  return 0;

}

/**********************************************************************
 * EB_FAS_Multigrid2D_Solver::Output_Multigrid_Elements --            *
 *                                                                    *
 * This routine writes out the node state information as elements for *
 * each multigrid grid level (finest to coarest) for a 1D array of 2D *
 * quadrilateral multi-block solution blocks.                         *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState,
	  class Quad_Soln_Block,
	  class Quad_Soln_Input_Parameters>
int EB_FAS_Multigrid2D_Solver<cState, pState,
			      Quad_Soln_Block,
			      Quad_Soln_Input_Parameters>::
Output_Multigrid_Elements(const int &number_of_time_steps,
			  const double &Time) {

  int i, i_output_title;
  char prefix1[256], prefix2[256], extension[256];
  char output_file_name1[256], output_file_name2[256];
  char *output_file_name_ptr;
  ofstream output_file1, output_file2;   

  // Determine prefix of output data file names.
  i = 0;
  while (1) {
    if (IP->Output_File_Name[i] == ' ' || IP->Output_File_Name[i] == '.') break;
    prefix1[i] = IP->Output_File_Name[i];
    i++;
    if (i > strlen(IP->Output_File_Name)) break;
  }
  prefix1[i] = '\0';
  strcpy(prefix2,prefix1);

  // Output to a seperate file for each multigrid level.
  for (int level = 0; level < IP->Multigrid_IP.Levels; level++) {

    // Determine prefix of output data file names for current grid level.
    strcpy(output_file_name1,prefix1);
    strcpy(output_file_name2,prefix2);
    sprintf(extension,"_%.2d",level);
    strcat(output_file_name1,extension);
    strcat(output_file_name2,extension);
    strcat(output_file_name1,"_active_cpu");
    strcat(output_file_name2,"_inactive_cpu");

    // Determine output data file name for this processor.
    sprintf(extension,"%.6d",List_of_Local_Solution_Blocks[level].ThisCPU);
    strcat(extension,".dat");
    strcat(output_file_name1,extension);
    strcat(output_file_name2,extension);

    // Open the output data files.
    output_file_name_ptr = output_file_name1;
    output_file1.open(output_file_name_ptr,ios::out);
    if (output_file1.bad()) return 1;
    output_file_name_ptr = output_file_name2;
    output_file2.open(output_file_name_ptr,ios::out);
    if (output_file2.bad()) return 1;
  
    // Write the solution data for each solution block.
    i_output_title = 1;
    for (int nb = 0; nb < List_of_Local_Solution_Blocks[level].Nblk; nb++) {
      if (List_of_Local_Solution_Blocks[level].Block[nb].used == ADAPTIVEBLOCK2D_USED) {
	EBSolver[level].Output_Active_Elements_Tecplot(nb,
						       number_of_time_steps,
						       Time,
						       i_output_title,
						       output_file1);
	EBSolver[level].Output_Inactive_Elements_Tecplot(nb,
							 number_of_time_steps,
							 Time,
							 i_output_title,
							 output_file2);
	if (i_output_title) i_output_title = 0;
      }
    }
  
    // Close the output data file.
    output_file1.close();
    output_file2.close();

  }

  // Writing of multigrid grids data files complete.
  return 0;

}

/**********************************************************************
 * EB_FAS_Multigrid2D_Solver::Output_Multigrid_Nodes --               *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState,
	  class Quad_Soln_Block,
	  class Quad_Soln_Input_Parameters>
int EB_FAS_Multigrid2D_Solver<cState, pState,
			      Quad_Soln_Block,
			      Quad_Soln_Input_Parameters>::
Output_Multigrid_Nodes(const int &number_of_time_steps,
		       const double &Time) {

  int i, i_output_title;
  char prefix[256], extension[256], output_file_name[256];
  char *output_file_name_ptr;
  ofstream output_file;   

  // Determine prefix of output data file names.
  i = 0;
  while (1) {
    if (IP->Output_File_Name[i] == ' ' || IP->Output_File_Name[i] == '.') break;
    prefix[i] = IP->Output_File_Name[i];
    i++;
    if (i > strlen(IP->Output_File_Name)) break;
  }
  prefix[i] = '\0';

  // Output to a seperate file for each multigrid level.
  for (int level = 0; level < IP->Multigrid_IP.Levels; level++) {

    // Determine prefix of output data file names for current grid level.
    strcpy(output_file_name,prefix);
    sprintf(extension,"_%.2d",level);
    strcat(output_file_name,extension);
    strcat(output_file_name,"_nodes_cpu");

    // Determine output data file name for this processor.
    sprintf(extension,"%.6d",List_of_Local_Solution_Blocks[level].ThisCPU);
    strcat(extension,".dat");
    strcat(output_file_name,extension);

    // Open the output data files.
    output_file_name_ptr = output_file_name;
    output_file.open(output_file_name_ptr,ios::out);
    if (output_file.bad()) return 1;

    // Write the solution data for each solution block.
    i_output_title = 1;
    for (int nb = 0; nb < List_of_Local_Solution_Blocks[level].Nblk; nb++) {
      if (List_of_Local_Solution_Blocks[level].Block[nb].used == ADAPTIVEBLOCK2D_USED) {
	EBSolver[level].Output_Nodes_Tecplot(nb,i_output_title,output_file);
	if (i_output_title) i_output_title = 0;
      }
    }
  
    // Close the output data file.
    output_file.close();

  }

  // Domain nodes data files successfully written.
  return 0;

}

/**********************************************************************
 * EB_FAS_Multigrid2D_Solver::Output_Multigrid_Cell_Status --         *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState,
	  class Quad_Soln_Block,
	  class Quad_Soln_Input_Parameters>
int EB_FAS_Multigrid2D_Solver<cState, pState,
			      Quad_Soln_Block,
			      Quad_Soln_Input_Parameters>::
Output_Multigrid_Cell_Status(const int &number_of_time_steps,
			     const double &Time) {

  int i, i_output_title;
  char prefix[256], extension[256], output_file_name[256];
  char *output_file_name_ptr;
  ofstream output_file;   

  // Determine prefix of output data file names.
  i = 0;
  while (1) {
    if (IP->Output_File_Name[i] == ' ' || IP->Output_File_Name[i] == '.') break;
    prefix[i] = IP->Output_File_Name[i];
    i++;
    if (i > strlen(IP->Output_File_Name)) break;
  }
  prefix[i] = '\0';

  // Output to a seperate file for each multigrid level.
  for (int level = 0; level < IP->Multigrid_IP.Levels; level++) {

    // Determine prefix of output data file names for current grid level.
    strcpy(output_file_name,prefix);
    sprintf(extension,"_%.2d",level);
    strcat(output_file_name,extension);
    strcat(output_file_name,"_status_cpu");

    // Determine output data file name for this processor.
    sprintf(extension,"%.6d",List_of_Local_Solution_Blocks[level].ThisCPU);
    strcat(extension,".dat");
    strcat(output_file_name,extension);

    // Open the output data files.
    output_file_name_ptr = output_file_name;
    output_file.open(output_file_name_ptr,ios::out);
    if (output_file.bad()) return 1;

    // Write the solution data for each solution block.
    i_output_title = 1;
    for (int nb = 0; nb < List_of_Local_Solution_Blocks[level].Nblk; nb++) {
      if (List_of_Local_Solution_Blocks[level].Block[nb].used == ADAPTIVEBLOCK2D_USED) {
	EBSolver[level].Output_Cell_Status_Tecplot(nb,i_output_title,output_file);
	if (i_output_title) i_output_title = 0;
      }
    }
  
    // Close the output data file.
    output_file.close();

  }

  // Domain nodes data files successfully written.
  return 0;

}

/**********************************************************************
 * EB_FAS_Multigrid2D_Solver::Restrict_Solution_Blocks --             *
 *                                                                    *
 * Restrict solution from Level_Fine to Level_Coarse for all blocks   *
 * on the local solution block list.  Note that the solution at the   *
 * coarse grid level is overwritten.                                  *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState,
	  class Quad_Soln_Block,
	  class Quad_Soln_Input_Parameters>
void EB_FAS_Multigrid2D_Solver<cState, pState,
			       Quad_Soln_Block,
			       Quad_Soln_Input_Parameters>::
Restrict_Solution_Blocks(const int &Level_Fine) {

  int i_fine, j_fine, Nghost, nghost;
  int Level_Coarse = Level_Fine + 1;
  double A, total_area;
  Polygon Pc, Pf;

  // Determine if the restriction includes the ghost cells.
  if (IP->Multigrid_IP.Apply_Coarse_Mesh_Boundary_Conditions) nghost = 0;
  else nghost = 1;

  // Loop through each solution block.
  for (int nb = 0; nb < List_of_Local_Solution_Blocks[Level_Fine].Nblk; nb++) {
    if (List_of_Local_Solution_Blocks[Level_Fine].Block[nb].used == ADAPTIVEBLOCK2D_USED) { 

      // Get the number of ghost cells on the fine grid.
      Nghost = Local_SolnBlk[Level_Fine][nb].Nghost;

      // Loop through the coarse grid cells.
      for (int i_coarse = Local_SolnBlk[Level_Coarse][nb].ICl-nghost; i_coarse <= Local_SolnBlk[Level_Coarse][nb].ICu+nghost; i_coarse++) {
	for (int j_coarse = Local_SolnBlk[Level_Coarse][nb].JCl-nghost; j_coarse <= Local_SolnBlk[Level_Coarse][nb].JCu+nghost; j_coarse++) {
	  // Determine the (i,j) index of the SW corner fine cell.
	  i_fine = 2*(i_coarse-Nghost)+Nghost;
	  j_fine = 2*(j_coarse-Nghost)+Nghost;

	  if (EBSolver[Level_Coarse].Mesh[nb].cell_status[i_coarse][j_coarse] != CELL_STATUS_ACTIVE) {
	    // The cell is inactive.  Set to standard atmosphere.
	    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse].Standard_Atmosphere();

	  } else if (!EBSolver[Level_Coarse].Interface_Union_List.Ni ||
		     !EBSolver[Level_Coarse].Adjustment_Data[nb].Interface_Present) {
	    // No embedded boundaries exist.  Restrict the solution
	    // information from the fine cells to the coarse cell using an
	    // area-weighted average of the fine cell solution information.
	    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse] = (Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].A*
								     Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine] +
								     Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][j_fine].A*
								     Local_SolnBlk[Level_Fine][nb].U[i_fine+1][j_fine] +
								     Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine+1].A*
								     Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine+1] +
								     Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][j_fine+1].A*
								     Local_SolnBlk[Level_Fine][nb].U[i_fine+1][j_fine+1])/
	                                                            (Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].A);
	  } else {

	    if (EBSolver[Level_Coarse].Mesh[nb].cell_type[i_coarse][j_coarse] == CELL_TYPE_QUADRILATERAL &&
		EBSolver[Level_Fine].Mesh[nb].cell_type[i_fine  ][j_fine  ] == CELL_TYPE_QUADRILATERAL &&
		EBSolver[Level_Fine].Mesh[nb].cell_type[i_fine+1][j_fine  ] == CELL_TYPE_QUADRILATERAL &&
		EBSolver[Level_Fine].Mesh[nb].cell_type[i_fine  ][j_fine+1] == CELL_TYPE_QUADRILATERAL &&
		EBSolver[Level_Fine].Mesh[nb].cell_type[i_fine+1][j_fine+1] == CELL_TYPE_QUADRILATERAL &&
		EBSolver[Level_Fine].Mesh[nb].cell_status[i_fine  ][j_fine  ] == CELL_STATUS_ACTIVE &&
		EBSolver[Level_Fine].Mesh[nb].cell_status[i_fine+1][j_fine  ] == CELL_STATUS_ACTIVE &&
		EBSolver[Level_Fine].Mesh[nb].cell_status[i_fine  ][j_fine+1] == CELL_STATUS_ACTIVE &&
		EBSolver[Level_Fine].Mesh[nb].cell_status[i_fine+1][j_fine+1] == CELL_STATUS_ACTIVE) {
	      // None of the associated cells have been adjusted.
	      // Restrict the solution information using an area-weighted
	      // average of the cells from the stored adjusted fine mesh
	      // and the coarse mesh.
	      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse] = (Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].A*
								       Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine] +
								       Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][j_fine].A*
								       Local_SolnBlk[Level_Fine][nb].U[i_fine+1][j_fine] +
								       Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine+1].A*
								       Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine+1] +
								       Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][j_fine+1].A*
								       Local_SolnBlk[Level_Fine][nb].U[i_fine+1][j_fine+1])/
	                                                              (Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].A);
	    } else {

	      // The coarse and/or fine cells have been adjusted.
	      // Restrict the solution information using an area-weighted
	      // average based on the area of intersection of the coarse
	      // grid and fine grid cells.
	      total_area = ZERO;
	      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse].Vacuum();
	      // Create the coarse cell polygon.
	      Pc.convert(Local_SolnBlk[Level_Coarse][nb].Grid.nodeSW(i_coarse,j_coarse).X,
			 Local_SolnBlk[Level_Coarse][nb].Grid.nodeSE(i_coarse,j_coarse).X,
			 Local_SolnBlk[Level_Coarse][nb].Grid.nodeNE(i_coarse,j_coarse).X,
			 Local_SolnBlk[Level_Coarse][nb].Grid.nodeNW(i_coarse,j_coarse).X);
	      // Search the fine cells for intersecting cells.
	      for (int j = (j_fine >= Local_SolnBlk[Level_Fine][nb].JCl) ? (j_fine-2) : (j_fine);
		   (j_fine <= Local_SolnBlk[Level_Fine][nb].JCu) ? (j < j_fine+4) : (j < j_fine+2); j++) {
		for (int i = (i_fine >= Local_SolnBlk[Level_Fine][nb].ICl) ? (i_fine-2) : (i_fine);
		     (i_fine <= Local_SolnBlk[Level_Fine][nb].ICu) ? (i < i_fine+4) : (i < i_fine+2); i++) {
		  if (EBSolver[Level_Fine].Mesh[nb].cell_status[i][j] == CELL_STATUS_ACTIVE) {
		    // Create the fine cell polygon.
		    Pf.convert(Local_SolnBlk[Level_Fine][nb].Grid.nodeSW(i,j).X,
			       Local_SolnBlk[Level_Fine][nb].Grid.nodeSE(i,j).X,
			       Local_SolnBlk[Level_Fine][nb].Grid.nodeNE(i,j).X,
			       Local_SolnBlk[Level_Fine][nb].Grid.nodeNW(i,j).X);
		    // Determine the area of intersection between the fine
		    // and coarse cell polygons.
		    A = Polygon_Intersection_Area(Pc,Pf);
		    Pf.deallocate();
		    // Add the area of intersection to the total area of
		    // intersection.
		    total_area += A;
		    // Restrict fine cell solution to the coarse cell
		    // based on the area of interesection.
		    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse] += Local_SolnBlk[Level_Fine][nb].U[i][j]*A;
		  }
		}
	      }
	      Pc.deallocate();
	      // Complete solution restriction by dividing the restricted
	      // solution by the total area of intersection.
	      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse] /= total_area;

	    }

	  }

	}
      }

    }
  }

}

/**********************************************************************
 * EB_FAS_Multigrid2D_Solver::Restrict_Residuals --                   *
 *                                                                    *
 * This routine restricts the residual dUdt from Level_Fine to        *
 * Level_Coarse (stored in P) for all blocks on the local solution    *
 * block list.  Note that the residual values on the coarse level are *
 * overwritten by this routine.                                       *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState,
	  class Quad_Soln_Block,
	  class Quad_Soln_Input_Parameters>
void EB_FAS_Multigrid2D_Solver<cState, pState,
			       Quad_Soln_Block,
			       Quad_Soln_Input_Parameters>::
Restrict_Residuals(const int &Level_Fine) {

  int i_fine, j_fine, Nghost, nghost;
  int Level_Coarse = Level_Fine + 1;
  double A, total_area;
  Polygon Pc, Pf;

  // Determine if the restriction includes the ghost cells.
  //if (IP->Multigrid_IP.Apply_Coarse_Mesh_Boundary_Conditions) nghost = 0;
  if (IP->Multigrid_IP.Apply_Coarse_Mesh_Boundary_Conditions) nghost = 0;
  else nghost = 1;

  // Loop through each local solution block.
  for (int nb = 0; nb < List_of_Local_Solution_Blocks[Level_Fine].Nblk; nb++) {
    if (List_of_Local_Solution_Blocks[Level_Fine].Block[nb].used == ADAPTIVEBLOCK2D_USED) { 

      // Get the number of ghost cells.
      Nghost = Local_SolnBlk[Level_Fine][nb].Nghost;

      // Loop through the coarse cells.
      for (int i_coarse = Local_SolnBlk[Level_Coarse][nb].ICl-nghost; i_coarse <= Local_SolnBlk[Level_Coarse][nb].ICu+nghost; i_coarse++) {
	for (int j_coarse = Local_SolnBlk[Level_Coarse][nb].JCl-nghost; j_coarse <= Local_SolnBlk[Level_Coarse][nb].JCu+nghost; j_coarse++) {
	  // Determine the (i,j) index of the corresponding SW corner
	  // fine cell.
	  i_fine = 2*(i_coarse-Nghost)+Nghost;
	  j_fine = 2*(j_coarse-Nghost)+Nghost;

	  if (EBSolver[Level_Coarse].Mesh[nb].cell_status[i_coarse][j_coarse] != CELL_STATUS_ACTIVE) {
	    // The cell is inactive.  Set to standard atmosphere.
	    MG_SolnBlk[Level_Coarse][nb].P[i_coarse][j_coarse].Standard_Atmosphere();

	  } else if (!EBSolver[Level_Coarse].Interface_Union_List.Ni ||
		     !EBSolver[Level_Coarse].Adjustment_Data[nb].Interface_Present) {
	    // No embedded boundaries exist.  Restrict the solution
	    // information from the fine cells to the coarse cell using an
	    // area-weighted average of the fine cell solution information.
	    MG_SolnBlk[Level_Coarse][nb].P[i_coarse][j_coarse] = (Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].A*
								  Local_SolnBlk[Level_Fine][nb].dUdt[i_fine][j_fine][0] +
								  Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][j_fine].A*
								  Local_SolnBlk[Level_Fine][nb].dUdt[i_fine+1][j_fine][0] +
								  Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine+1].A*
								  Local_SolnBlk[Level_Fine][nb].dUdt[i_fine][j_fine+1][0] +
								  Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][j_fine+1].A*
								  Local_SolnBlk[Level_Fine][nb].dUdt[i_fine+1][j_fine+1][0])/
	                                                         (Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].A);

	  } else {

	    if (EBSolver[Level_Coarse].Mesh[nb].cell_type[i_coarse][j_coarse] == CELL_TYPE_QUADRILATERAL &&
		EBSolver[Level_Fine].Mesh[nb].cell_type[i_fine  ][j_fine  ] == CELL_TYPE_QUADRILATERAL &&
		EBSolver[Level_Fine].Mesh[nb].cell_type[i_fine+1][j_fine  ] == CELL_TYPE_QUADRILATERAL &&
		EBSolver[Level_Fine].Mesh[nb].cell_type[i_fine  ][j_fine+1] == CELL_TYPE_QUADRILATERAL &&
		EBSolver[Level_Fine].Mesh[nb].cell_type[i_fine+1][j_fine+1] == CELL_TYPE_QUADRILATERAL &&
		EBSolver[Level_Fine].Mesh[nb].cell_status[i_fine  ][j_fine  ] == CELL_STATUS_ACTIVE &&
		EBSolver[Level_Fine].Mesh[nb].cell_status[i_fine+1][j_fine  ] == CELL_STATUS_ACTIVE &&
		EBSolver[Level_Fine].Mesh[nb].cell_status[i_fine  ][j_fine+1] == CELL_STATUS_ACTIVE &&
		EBSolver[Level_Fine].Mesh[nb].cell_status[i_fine+1][j_fine+1] == CELL_STATUS_ACTIVE) {
	      // None of the associated cells have been adjusted.
	      // Restrict the solution information using an area-weighted
	      // average of the cells from the stored adjusted fine mesh
	      // and the coarse mesh.
	      MG_SolnBlk[Level_Coarse][nb].P[i_coarse][j_coarse] = (Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].A*
								    Local_SolnBlk[Level_Fine][nb].dUdt[i_fine][j_fine][0] +
								    Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][j_fine].A*
								    Local_SolnBlk[Level_Fine][nb].dUdt[i_fine+1][j_fine][0] +
								    Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine+1].A*
								    Local_SolnBlk[Level_Fine][nb].dUdt[i_fine][j_fine+1][0] +
								    Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][j_fine+1].A*
								    Local_SolnBlk[Level_Fine][nb].dUdt[i_fine+1][j_fine+1][0])/
	                                                           (Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].A);

	    } else {
	      // The coarse and/or fine cells have been adjusted.
	      // Restrict the solution information using an area-weight
	      // average based on the area of intersection of the coarse
	      // grid and fine grid cells.
	      total_area = ZERO;
	      MG_SolnBlk[Level_Coarse][nb].P[i_coarse][j_coarse].Vacuum();
	      // Create the coarse cell polygon.
	      Pc.convert(Local_SolnBlk[Level_Coarse][nb].Grid.nodeSW(i_coarse,j_coarse).X,
			 Local_SolnBlk[Level_Coarse][nb].Grid.nodeSE(i_coarse,j_coarse).X,
			 Local_SolnBlk[Level_Coarse][nb].Grid.nodeNE(i_coarse,j_coarse).X,
			 Local_SolnBlk[Level_Coarse][nb].Grid.nodeNW(i_coarse,j_coarse).X);
	      // Search the fine cells for intersecting cells.
	      for (int j = (j_fine >= Local_SolnBlk[Level_Fine][nb].JCl) ? (j_fine-2) : (j_fine);
		   (j_fine <= Local_SolnBlk[Level_Fine][nb].JCu) ? (j < j_fine+4) : (j < j_fine+2); j++) {
		for (int i = (i_fine >= Local_SolnBlk[Level_Fine][nb].ICl) ? (i_fine-2) : (i_fine);
		     (i_fine <= Local_SolnBlk[Level_Fine][nb].ICu) ? (i < i_fine+4) : (i < i_fine+2); i++) {
  		  if (EBSolver[Level_Fine].Mesh[nb].cell_status[i][j] == CELL_STATUS_ACTIVE) {
		    // Create the fine cell polygon.
		    Pf.convert(Local_SolnBlk[Level_Fine][nb].Grid.nodeSW(i,j).X,
			       Local_SolnBlk[Level_Fine][nb].Grid.nodeSE(i,j).X,
			       Local_SolnBlk[Level_Fine][nb].Grid.nodeNE(i,j).X,
			       Local_SolnBlk[Level_Fine][nb].Grid.nodeNW(i,j).X);
		    // Determine the area of intersection between the fine
		    // and coarse cell polygons.
		    A = Polygon_Intersection_Area(Pc,Pf);
		    Pf.deallocate();
		    // Add the area of intersection to the total area of
		    // intersection.
		    total_area += A;
		    // Restrict fine cell solution to the coarse cell
		    // based on the area of interesection.
		    MG_SolnBlk[Level_Coarse][nb].P[i_coarse][j_coarse] += Local_SolnBlk[Level_Fine][nb].dUdt[i][j][0]*A;
		  }
		}
	      }
	      Pc.deallocate();
	      // Complete solution restriction by dividing the restricted
	      // solution by the total area of intersection.
	      MG_SolnBlk[Level_Coarse][nb].P[i_coarse][j_coarse] /= total_area;

	    }

	  }

	}
      }

    }
  }

}

/**********************************************************************
 * EB_FAS_Multigrid2D_Solver::Restrict_Boundary_Ref_States --         *
 *                                                                    *
 * This routine restricts the boundary refeference states (UoN/S/E/W) *
 * from Level_Fine to Level_Coarse for all solution blocks.  The      *
 * values on the coarse grid level are overwritten by this routine.   *
 * The restriction operator used is area weighted average.            *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState,
	  class Quad_Soln_Block,
	  class Quad_Soln_Input_Parameters>
void EB_FAS_Multigrid2D_Solver<cState, pState,
			       Quad_Soln_Block,
			       Quad_Soln_Input_Parameters>::
Restrict_Boundary_Ref_States(const int &Level_Fine) {

  int i_fine, j_fine, ICl, ICu, JCl, JCu, Nghost;
  int Level_Coarse = Level_Fine + 1;
  
  // Loop through each local solution block.
  for (int nb = 0; nb < List_of_Local_Solution_Blocks[Level_Fine].Nblk; nb++) {
    if (List_of_Local_Solution_Blocks[Level_Fine].Block[nb].used == ADAPTIVEBLOCK2D_USED) { 

      // Get the number of ghost cells.
      Nghost = Local_SolnBlk[Level_Fine][nb].Nghost;

      // Store the lower and upper j-indices on the fine grid.
      JCl = Local_SolnBlk[Level_Fine][nb].JCl;
      JCu = Local_SolnBlk[Level_Fine][nb].JCu;

      // Loop through the i-direction cells of the coarse grid.
      for (int i_coarse = Local_SolnBlk[Level_Coarse][nb].ICl; i_coarse <= Local_SolnBlk[Level_Coarse][nb].ICu; i_coarse++) {
	// Determine the i-index of the corresponding SW corner
	// fine cell.
	i_fine = 2*(i_coarse-Nghost)+Nghost;
	// Restrict the north boundary reference state.
	Local_SolnBlk[Level_Coarse][nb].WoN[i_coarse] = 
	  (Local_SolnBlk[Level_Fine][nb].WoN[i_fine]*
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][JCu].A +
	   Local_SolnBlk[Level_Fine][nb].WoN[i_fine+1]*
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][JCu].A)/
	  (Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][JCu].A +
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][JCu].A);
	// Restrict the south boundary reference state.
	Local_SolnBlk[Level_Coarse][nb].WoS[i_coarse] = 
	  (Local_SolnBlk[Level_Fine][nb].WoS[i_fine]*
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][JCl].A +
	   Local_SolnBlk[Level_Fine][nb].WoS[i_fine+1]*
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][JCl].A)/
	  (Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][JCl].A +
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][JCl].A);
      }

      // Store the lower and upper i-indices on the fine grid.
      ICl = Local_SolnBlk[Level_Fine][nb].ICl;
      ICu = Local_SolnBlk[Level_Fine][nb].ICu;

      // Loop through the j-direction cells of the coarse grid.
      for (int j_coarse = Local_SolnBlk[Level_Coarse][nb].JCl; j_coarse <= Local_SolnBlk[Level_Coarse][nb].JCu; j_coarse++) {
	// Determine the j-index of the corresponding SW corner
	// fine cell.
	j_fine = 2*(j_coarse-Nghost)+Nghost;
	// Restrict the west boundary reference state.
	Local_SolnBlk[Level_Coarse][nb].WoW[j_coarse] = 
	  (Local_SolnBlk[Level_Fine][nb].WoW[j_fine]*
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[ICl][j_fine].A +
	   Local_SolnBlk[Level_Fine][nb].WoW[j_fine+1]*
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[ICl][j_fine+1].A)/
	  (Local_SolnBlk[Level_Fine][nb].Grid.Cell[ICl][j_fine].A +
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[ICl][j_fine+1].A);
	// Restrict the east boundary reference state.
	Local_SolnBlk[Level_Coarse][nb].WoE[j_coarse] = 
	  (Local_SolnBlk[Level_Fine][nb].WoE[j_fine]*
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[ICu][j_fine].A +
	   Local_SolnBlk[Level_Fine][nb].WoE[j_fine+1]*
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[ICu][j_fine+1].A)/
	  (Local_SolnBlk[Level_Fine][nb].Grid.Cell[ICu][j_fine].A +
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[ICu][j_fine+1].A);
      }

    }
  }

}

/**********************************************************************
 * EB_FAS_Multigrid2D_Solver::Prolong_Solution_Blocks --              *
 *                                                                    *
 * This routine conducts the prolongation of the solution states      *
 * stored in U from Level_Coarse to Level_Fine for all blocks on the  *
 * local solution block list.  Bilinear interpolation is used as the  *
 * prolongation operator unless the fine grid cell centre falls       *
 * outside of the associated coarse grid interpolation points.  In    *
 * this case the prolongation is conducted using injection.           *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState,
	  class Quad_Soln_Block,
	  class Quad_Soln_Input_Parameters>
int EB_FAS_Multigrid2D_Solver<cState, pState,
			      Quad_Soln_Block,
			      Quad_Soln_Input_Parameters>::
Prolong_Solution_Blocks(const int &Level_Coarse) {

  int error_flag, i_fine, j_fine, Nghost;
  int Level_Fine = Level_Coarse - 1;
  cState Ufine;
  bool injection_for_this_coarse_cell;
  injection_for_this_coarse_cell = false;
  int i_coarse_min, j_coarse_min, coarse_cell_found;
  double distance;

  // Loop through each local solution block.
  for (int nb = 0; nb < List_of_Local_Solution_Blocks[Level_Coarse].Nblk; nb++) {
    if (List_of_Local_Solution_Blocks[Level_Coarse].Block[nb].used == ADAPTIVEBLOCK2D_USED) { 

      // Get the number of ghost cells on the coarse grid.
      Nghost = Local_SolnBlk[Level_Coarse][nb].Nghost;

      // Loop through the coarse grid cells.
      for (int i_coarse = Local_SolnBlk[Level_Coarse][nb].ICl; i_coarse <= Local_SolnBlk[Level_Coarse][nb].ICu; i_coarse++) {
	for (int j_coarse = Local_SolnBlk[Level_Coarse][nb].JCl; j_coarse <= Local_SolnBlk[Level_Coarse][nb].JCu; j_coarse++) {

 	  if (EBSolver[Level_Coarse].Mesh[nb].cell_status[i_coarse-1][j_coarse-1] != CELL_STATUS_ACTIVE ||
 	      EBSolver[Level_Coarse].Mesh[nb].cell_status[i_coarse  ][j_coarse-1] != CELL_STATUS_ACTIVE ||
 	      EBSolver[Level_Coarse].Mesh[nb].cell_status[i_coarse+1][j_coarse-1] != CELL_STATUS_ACTIVE ||
 	      EBSolver[Level_Coarse].Mesh[nb].cell_status[i_coarse-1][j_coarse  ] != CELL_STATUS_ACTIVE ||
 	      EBSolver[Level_Coarse].Mesh[nb].cell_status[i_coarse  ][j_coarse  ] != CELL_STATUS_ACTIVE ||
 	      EBSolver[Level_Coarse].Mesh[nb].cell_status[i_coarse+1][j_coarse  ] != CELL_STATUS_ACTIVE ||
 	      EBSolver[Level_Coarse].Mesh[nb].cell_status[i_coarse-1][j_coarse+1] != CELL_STATUS_ACTIVE ||
 	      EBSolver[Level_Coarse].Mesh[nb].cell_status[i_coarse  ][j_coarse+1] != CELL_STATUS_ACTIVE ||
 	      EBSolver[Level_Coarse].Mesh[nb].cell_status[i_coarse+1][j_coarse+1] != CELL_STATUS_ACTIVE) {
 	    // The fine cell is near an embedded boundary.  Prolong by 
 	    // injecting the solution information from the nearest
 	    // active cell.
 	    i_fine = 2*(i_coarse-Nghost)+Nghost;
 	    j_fine = 2*(j_coarse-Nghost)+Nghost;
 	    for (int jj = j_fine; jj <= j_fine+1; jj++) {
 	      for (int ii = i_fine; ii <= i_fine+1; ii++) {
 		if (EBSolver[Level_Fine].Mesh[nb].cell_status[ii][jj] != CELL_STATUS_ACTIVE) {
		  // The fine grid cell is inactive.  Set to standard
		  // atmosphere.
 		  Local_SolnBlk[Level_Fine][nb].U[ii][jj].Standard_Atmosphere();
 		} else {
		  // The fine grid cell is inactive.  Inject the
		  // solution information from the nearest coarse grid
		  // cell.
 		  distance = MILLION;
 		  coarse_cell_found = OFF;
 		  for (int jc = j_coarse-1; jc <= j_coarse+1; jc++) {
 		    for (int ic = i_coarse-1; ic <= i_coarse+1; ic++) {
 		      if (EBSolver[Level_Coarse].Mesh[nb].cell_status[ic][jc] == CELL_STATUS_ACTIVE &&
 			  abs(Local_SolnBlk[Level_Coarse][nb].Grid.Cell[ic][jc].Xc -
 			      Local_SolnBlk[Level_Fine][nb].Grid.Cell[ii][jj].Xc) < distance) {
 			i_coarse_min = ic; j_coarse_min = jc;
 			distance = abs(Local_SolnBlk[Level_Coarse][nb].Grid.Cell[ic][jc].Xc -
 				       Local_SolnBlk[Level_Fine][nb].Grid.Cell[ii][jj].Xc);
 			coarse_cell_found = ON;
 		      }
 		    }
 		  }
 		  if (!coarse_cell_found) return 27081;
		  Local_SolnBlk[Level_Fine][nb].U[ii][jj] = Local_SolnBlk[Level_Coarse][nb].U[i_coarse_min][j_coarse_min];
 		}
 	      }
 	    }

 	  } else {

 	    // The fine cell is not near an embedded boundary.  Prolong
	    // by either by injection or a bilinear interpolation.

	    // Set injection flag if prolongation by injection is desired.
	    if (IP->Multigrid_IP.Prolong_Using_Injection) injection_for_this_coarse_cell = true;

	    // If the bilinear interpolation is used as the prolongation
	    // operator and a valid set of interpolation points can not
	    // be found for one or more of the fine grid cells then
	    // switch the prolongation operator to injection for all
	    // four fine grid cells.  This is facilitated by the use of
	    // a goto statement where 'restart' is used as the goto
	    // label.
	  restart: ;

	    /**********************************************************
	     * Determine the cell-centred solution state for the      *
	     * south-west fine grid cell.                             *
	     **********************************************************/

	    // Determine the (i,j) index of the SW fine grid cell.
	    i_fine = 2*(i_coarse-Nghost)+Nghost;
	    j_fine = 2*(j_coarse-Nghost)+Nghost;

	    // Attempt to use bilinear interpolation if the injection
	    // flag has not been set true.
	    if (!injection_for_this_coarse_cell) {

	      // Attempt to use the SW coarse grid solution states as
	      // the interpolants ([i-1,i],[j-1,j]).
	      error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse-1],
						  Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse-1].Xc,
						  Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse],
						  Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse].Xc,
						  Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
						  Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
						  Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse-1],
						  Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse-1].Xc,
						  Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
						  Ufine);
	      // If fine grid cell-centre falls outside of the SW
	      // coarse grid interpolants then try the NW coarse grid
	      // solution states as the interpolants ([i-1,i],[j,j+1]).
	      if (error_flag) {
		error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse].Xc,
						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse+1],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse+1].Xc,
						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse+1],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse+1].Xc,
						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
						    Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
						    Ufine);
	      }
	      // If fine grid cell-centre falls outside of the NW
	      // coarse grid interpolants then try the SE coarse grid
	      // solution states as the interpolants ([i,i+1],[j-1,j]).
	      if (error_flag) {
		error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse-1],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse-1].Xc,
						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse].Xc,
						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse-1],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse-1].Xc,
						    Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
						    Ufine);
	      }
	      // If fine grid cell-centre falls outside of the SE
	      // coarse grid interpolants then try the NE coarse grid
	      // solution states as the interpolants ([i,i+1],[j,j+1]).
	      if (error_flag != 0) {
		error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse+1],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse+1].Xc,
						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse+1],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse+1].Xc,
						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse].Xc,
						    Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
						    Ufine);
	      }
	      // If none of the four sets of possible interpolants 
	      // provide a suitable set of interpolants then restart
	      // the prolongation associated with the current coarse
	      // grid cell using injection.
	      if (error_flag) {
		injection_for_this_coarse_cell = true;
		goto restart;
	      }

	    } else {

	      // Perform the prolongation using injection.
	      Ufine = Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse];

	    }

	    // Assign the SW fine grid cell solution state.
	    Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine] = Ufine;

	    /**********************************************************
	     * Determine the cell-centred solution state for the      *
	     * south-east fine grid cell.                             *
	     **********************************************************/

	    // Determine the (i,j) index of the SE fine grid cell.
	    i_fine = 2*(i_coarse-Nghost)+Nghost+1;
	    j_fine = 2*(j_coarse-Nghost)+Nghost;

	    // Attempt to use bilinear interpolation if the injection
	    // flag has not been set true.
	    if (!injection_for_this_coarse_cell) {

	      // Attempt to use the SE coarse grid solution states as
	      // the interpolants ([i,i+1],[j-1,j]).
	      error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse-1],
						  Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse-1].Xc,
						  Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
						  Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
						  Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse],
						  Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse].Xc,
						  Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse-1],
						  Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse-1].Xc,
						  Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
						  Ufine);
	      // If fine grid cell-centre falls outside of the SE coarse
	      // grid interpolants then try the NE coarse grid solution
	      // states as the interpolants ([i,i+1],[j,j+1]).
	      if (error_flag) {
		error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse+1],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse+1].Xc,
						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse+1],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse+1].Xc,
						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse].Xc,
						    Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
						    Ufine);
	      }
	      // If fine grid cell-centre falls outside of the NE coarse
	      // grid interpolants then try the SW coarse grid solution
	      // states as the interpolants ([i-1,i],[j-1,j]).
	      if (error_flag) {
		error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse-1],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse-1].Xc,
						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse].Xc,
						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse-1],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse-1].Xc,
						    Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
						    Ufine);
	      }
	      // If fine grid cell-centre falls outside of the SW coarse
	      // grid interpolants then try the N coarse grid solution
	      // states as the interpolants ([i-1,i],[j,j+1]).
	      if (error_flag) {
		error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse].Xc,
						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse+1],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse+1].Xc,
						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse+1],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse+1].Xc,
						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
						    Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
						    Ufine);
	      }
	      // If none of the four sets of possible interpolants 
	      // provide a suitable set of interpolants then restart
	      // the prolongation associated with the current coarse
	      // grid cell using injection.
	      if (error_flag != 0) {
		injection_for_this_coarse_cell = true;
		goto restart;
	      }

	    } else {

	      // Perform the prolongation using injection.
	      Ufine = Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse];

	    }

	    // Assign the SE fine grid cell solution state.
	    Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine] = Ufine;

	    /************************************************************
	     * Determine the cell-centred solution state for the north- *
	     * east fine grid cell.                                     *
	     ************************************************************/

	    // Determine the (i,j) index of the NE fine grid cell.
	    i_fine = 2*(i_coarse-Nghost)+Nghost+1;
	    j_fine = 2*(j_coarse-Nghost)+Nghost+1;

	    // Attempt to use bilinear interpolation if the injection flag
	    // has not been set true.
	    if (!injection_for_this_coarse_cell) {

	      // Attempt to use the NE coarse grid solution states as the
	      // interpolants ([i,i+1],[j,j+1]).
	      error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
						  Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
						  Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse+1],
						  Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse+1].Xc,
						  Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse+1],
						  Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse+1].Xc,
						  Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse],
						  Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse].Xc,
						  Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
						  Ufine);
	      // If fine grid cell-centre falls outside of the NE coarse
	      // grid interpolants then try the SE coarse grid solution
	      // states as the interpolants ([i,i+1],[j-1,j]).
	      if (error_flag) {
		error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse-1],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse-1].Xc,
						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse].Xc,
						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse-1],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse-1].Xc,
						    Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
						    Ufine);
	      }
	      // If fine grid cell-centre falls outside of the SE coarse
	      // grid interpolants then try the NW coarse grid solution
	      // states as the interpolants ([i-1,i],[j,j+1]).
	      if (error_flag) {
		error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse].Xc,
						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse+1],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse+1].Xc,
						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse+1],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse+1].Xc,
						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
						    Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
						    Ufine);
	      }
	      // If fine grid cell-centre falls outside of the NW coarse
	      // grid interpolants then try the SW coarse grid solution
	      // states as the interpolants ([i-1,i],[j-1,j]).
	      if (error_flag) {
		error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse-1],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse-1].Xc,
						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse].Xc,
						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse-1],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse-1].Xc,
						    Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
						    Ufine);
	      }
	      // If none of the four sets of possible interpolants 
	      // provide a suitable set of interpolants then restart
	      // the prolongation associated with the current coarse
	      // grid cell using injection.
	      if (error_flag != 0) {
		injection_for_this_coarse_cell = true;
		goto restart;
	      }

	    } else {

	      // Perform the prolongation using injection.
	      Ufine = Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse];

	    }

	    // Assign the NE fine grid cell solution state.
	    Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine] = Ufine;

	    /************************************************************
	     * Determine the cell-centred solution state for the north- *
	     * west fine grid cell.                                     *
	     ************************************************************/

	    // Determine the (i,j) index of the NW fine grid cell.
	    i_fine = 2*(i_coarse-Nghost)+Nghost;
	    j_fine = 2*(j_coarse-Nghost)+Nghost+1;

	    // Attempt to use bilinear interpolation if the injection flag
	    // has not been set true.
	    if (!injection_for_this_coarse_cell) {

	      // Attempt to use the NW coarse grid solution states as the
	      // interpolants ([i-1,i],[j,j+1]).
	      error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse],
						  Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse].Xc,
						  Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse+1],
						  Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse+1].Xc,
						  Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse+1],
						  Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse+1].Xc,
						  Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
						  Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
						  Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
						  Ufine);
	      // If fine grid cell-centre falls outside of the NW coarse
	      // grid interpolants then try the NE coarse grid solution
	      // states as the interpolants ([i,i+1],[j,j+1]).
	      if (error_flag) {
		error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse+1],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse+1].Xc,
						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse+1],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse+1].Xc,
						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse].Xc,
						    Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
						    Ufine);
	      }
	      // If fine grid cell-centre falls outside of the NE coarse
	      // grid interpolants then try the SW coarse grid solution
	      // states as the interpolants ([i-1,i],[j-1,j]).
	      if (error_flag) {
		error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse-1],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse-1].Xc,
						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse].Xc,
						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse-1],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse-1].Xc,
						    Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
						    Ufine);
	      }
	      // If fine grid cell-centre falls outside of the SW coarse
	      // grid interpolants then try the SE coarse grid solution
	      // states as the interpolants ([i,i+1],[j-1,j]).
	      if (error_flag) {
		error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse-1],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse-1].Xc,
						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse].Xc,
						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse-1],
						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse-1].Xc,
						    Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
						    Ufine);
	      }
	      // If none of the four sets of possible interpolants 
	      // provide a suitable set of interpolants then restart
	      // the prolongation associated with the current coarse
	      // grid cell using injection.
	      if (error_flag != 0) {
		injection_for_this_coarse_cell = true;
		goto restart;
	      }

	    } else {

	      // Perform the prolongation using injection.
	      Ufine = Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse];

	    }

	    // Assign the NW fine grid cell solution state.
	    Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine] = Ufine;

	    // Reset the injection flag to false.
	    injection_for_this_coarse_cell = false;

	  }

	}
      }

    }
  }

  // Prolongation of the solution block was successful.
  return 0;

}

/**********************************************************************
 * EB_FAS_Multigrid2D_Solver::Prolong_and_Update_Solution_Blocks --   *
 *                                                                    *
 * This routine conducts the prolongation of the solution changes     *
 * stored in U from the coarse grid to the fine grid for all of       *
 * solution blocks on the local solution block list and adds the      *
 * changes to the solution U on the fine grid solution blocks.        *
 * Bilinear interpolation is used as the prolongation operator unless *
 * the fine grid cell centre falls outside of the associated coarse   *
 * grid interpolation points.  In this case the prolongation is       *
 * conducted using injection.  Note that injection is also used at    *
 * all Dirichlet-type boundary conditions and Neumann boundary        *
 * conditions are enforced on the coarse grid.                        *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState,
	  class Quad_Soln_Block,
	  class Quad_Soln_Input_Parameters>
int EB_FAS_Multigrid2D_Solver<cState, pState,
			      Quad_Soln_Block,
			      Quad_Soln_Input_Parameters>::
Prolong_and_Update_Solution_Blocks(const int &Level_Coarse) {

  int error_flag, i_fine, j_fine, Nghost, residual_reduction_number;
  int Level_Fine = Level_Coarse - 1;
  cState Uchange, Unew;
  bool injection_for_this_coarse_cell;
  injection_for_this_coarse_cell = false;
  int i_coarse_min, j_coarse_min, coarse_cell_found;
  double distance;

  // Loop through each local solution block.
  for (int nb = 0; nb < List_of_Local_Solution_Blocks[Level_Coarse].Nblk; nb++) {
    if (List_of_Local_Solution_Blocks[Level_Coarse].Block[nb].used == ADAPTIVEBLOCK2D_USED) { 

      // Get the number of ghost cells on the coarse grid.
      Nghost = Local_SolnBlk[Level_Coarse][nb].Nghost;

      // Enforce Neumann boundary conditions as appropriate on the 
      // coarse grid if bilinear interpolation is used as the 
      // prolongation operator.
      if (!IP->Multigrid_IP.Prolong_Using_Injection) {

	// West and east face Neumann boundary conditions.
	for (int j_coarse = Local_SolnBlk[Level_Coarse][nb].JCl; j_coarse <= Local_SolnBlk[Level_Coarse][nb].JCu; j_coarse++) {
	  // East boundary.
	  if (Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeE[j_coarse] == BC_CONSTANT_EXTRAPOLATION) {
	    Local_SolnBlk[Level_Coarse][nb].U[Local_SolnBlk[Level_Coarse][nb].ICu+1][j_coarse] =
	      Local_SolnBlk[Level_Coarse][nb].U[Local_SolnBlk[Level_Coarse][nb].ICu][j_coarse];
	  }
	  // West boundary.
	  if (Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeW[j_coarse] == BC_CONSTANT_EXTRAPOLATION) {
	    Local_SolnBlk[Level_Coarse][nb].U[Local_SolnBlk[Level_Coarse][nb].ICl-1][j_coarse] =
	      Local_SolnBlk[Level_Coarse][nb].U[Local_SolnBlk[Level_Coarse][nb].ICl][j_coarse];
	  }
	}

	// North and south face Neumann boundary coditions.
	for (int i_coarse = Local_SolnBlk[Level_Coarse][nb].ICl; i_coarse <= Local_SolnBlk[Level_Coarse][nb].ICu; i_coarse++) {
	  // South boundary.
	  if (Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeS[i_coarse] == BC_CONSTANT_EXTRAPOLATION) {
	    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][Local_SolnBlk[Level_Coarse][nb].JCl-1] =
	      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][Local_SolnBlk[Level_Coarse][nb].JCl];
	  }
	  // North boundary.
	  if (Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeN[i_coarse] == BC_CONSTANT_EXTRAPOLATION) {
	    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][Local_SolnBlk[Level_Coarse][nb].JCu+1] =
	      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][Local_SolnBlk[Level_Coarse][nb].JCu];
	  }
	}
      }

      // Loop through the coarse grid cells.
      for (int i_coarse = Local_SolnBlk[Level_Coarse][nb].ICl; i_coarse <= Local_SolnBlk[Level_Coarse][nb].ICu; i_coarse++) {
	for (int j_coarse = Local_SolnBlk[Level_Coarse][nb].JCl; j_coarse <= Local_SolnBlk[Level_Coarse][nb].JCu; j_coarse++) {

	  if (EBSolver[Level_Coarse].Mesh[nb].cell_status[i_coarse-1][j_coarse-1] != CELL_STATUS_ACTIVE ||
	      EBSolver[Level_Coarse].Mesh[nb].cell_status[i_coarse  ][j_coarse-1] != CELL_STATUS_ACTIVE ||
	      EBSolver[Level_Coarse].Mesh[nb].cell_status[i_coarse+1][j_coarse-1] != CELL_STATUS_ACTIVE ||
	      EBSolver[Level_Coarse].Mesh[nb].cell_status[i_coarse-1][j_coarse  ] != CELL_STATUS_ACTIVE ||
	      EBSolver[Level_Coarse].Mesh[nb].cell_status[i_coarse  ][j_coarse  ] != CELL_STATUS_ACTIVE ||
	      EBSolver[Level_Coarse].Mesh[nb].cell_status[i_coarse+1][j_coarse  ] != CELL_STATUS_ACTIVE ||
	      EBSolver[Level_Coarse].Mesh[nb].cell_status[i_coarse-1][j_coarse+1] != CELL_STATUS_ACTIVE ||
	      EBSolver[Level_Coarse].Mesh[nb].cell_status[i_coarse  ][j_coarse+1] != CELL_STATUS_ACTIVE ||
	      EBSolver[Level_Coarse].Mesh[nb].cell_status[i_coarse+1][j_coarse+1] != CELL_STATUS_ACTIVE) {
	    // The fine cell is near an embedded boundary.  Prolong by 
	    // injecting the solution information from the nearest
	    // active cell.
	    i_fine = 2*(i_coarse-Nghost)+Nghost;
	    j_fine = 2*(j_coarse-Nghost)+Nghost;
	    for (int jj = j_fine; jj <= j_fine+1; jj++) {
	      for (int ii = i_fine; ii <= i_fine+1; ii++) {
		if (EBSolver[Level_Fine].Mesh[nb].cell_status[ii][jj] != CELL_STATUS_ACTIVE) {
		  // The fine grid cell is inactive.  Set to standard
		  // atmosphere.
		  Uchange.Standard_Atmosphere();

		} else {
		  // The fine grid cell is inactive.  Inject the
		  // solution information from the nearest coarse grid
		  // cell.
		  distance = MILLION;
		  coarse_cell_found = OFF;
		  for (int jc = j_coarse-1; jc <= j_coarse+1; jc++) {
		    for (int ic = i_coarse-1; ic <= i_coarse+1; ic++) {
		      if (EBSolver[Level_Coarse].Mesh[nb].cell_status[ic][jc] == CELL_STATUS_ACTIVE &&
			  abs(Local_SolnBlk[Level_Coarse][nb].Grid.Cell[ic][jc].Xc -
			      Local_SolnBlk[Level_Fine][nb].Grid.Cell[ii][jj].Xc) < distance) {
			i_coarse_min = ic; j_coarse_min = jc;
			distance = abs(Local_SolnBlk[Level_Coarse][nb].Grid.Cell[ic][jc].Xc -
				       Local_SolnBlk[Level_Fine][nb].Grid.Cell[ii][jj].Xc);
			coarse_cell_found = ON;
		      }
		    }
		  }
		  if (!coarse_cell_found) { cout << endl << Level_Coarse << Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc; cout.flush(); return 27082; }
		  Uchange = Local_SolnBlk[Level_Coarse][nb].U[i_coarse_min][j_coarse_min];
		}

		// Update the conservative solution state in the fine
		// grid cell.  If required, use residual reduction if
		// the update creates a unphysical state (negative 
		// density or pressure).  Residual reduction is
		// performed by halving the residual update.
		if (IP->Multigrid_IP.Update_Stability_Switch) {
		  Unew = Local_SolnBlk[Level_Fine][nb].U[ii][jj] + Uchange;
		  if (Unew.Unphysical_Properties()) {
		    residual_reduction_number = 1;
		    while (residual_reduction_number <= IP->Multigrid_IP.Maximum_Number_of_Update_Reductions) {
		      Uchange *= HALF;
		      Unew = Local_SolnBlk[Level_Fine][nb].U[ii][jj] + Uchange;
		      if (!Unew.Unphysical_Properties()) break;
		      residual_reduction_number++;
		    }
		    if (residual_reduction_number > IP->Multigrid_IP.Maximum_Number_of_Update_Reductions) Uchange *= ZERO;
		  }
		}
		Local_SolnBlk[Level_Fine][nb].U[ii][jj] += Uchange;
	      }
	    }

	  } else {

	    // If the coarse grid cell is next to a Dirichlet (fixed) 
	    // boundary condition then prolong using injection instead of 
	    // interpolation for all four of the fine cells.  Interpolation
	    // for these cells can lead to a skewed (unphysical) result.
	    if ((i_coarse == Local_SolnBlk[Level_Coarse][nb].ICl &&
		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeW[j_coarse] == BC_FIXED) ||
		(i_coarse == Local_SolnBlk[Level_Coarse][nb].ICu &&
		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeE[j_coarse] == BC_FIXED) ||
		(j_coarse == Local_SolnBlk[Level_Coarse][nb].JCl &&
		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeS[i_coarse] == BC_FIXED) ||
		(j_coarse == Local_SolnBlk[Level_Coarse][nb].JCu &&
		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeN[i_coarse] == BC_FIXED) ||
		(i_coarse == Local_SolnBlk[Level_Coarse][nb].ICl &&
		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeW[j_coarse] == BC_INFLOW_SUBSONIC) ||
		(i_coarse == Local_SolnBlk[Level_Coarse][nb].ICu &&
		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeE[j_coarse] == BC_INFLOW_SUBSONIC) ||
		(j_coarse == Local_SolnBlk[Level_Coarse][nb].JCl &&
		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeS[i_coarse] == BC_INFLOW_SUBSONIC) ||
		(j_coarse == Local_SolnBlk[Level_Coarse][nb].JCu &&
		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeN[i_coarse] == BC_INFLOW_SUBSONIC) ||
		(i_coarse == Local_SolnBlk[Level_Coarse][nb].ICl &&
		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeW[j_coarse] == BC_OUTFLOW_SUBSONIC) ||
		(i_coarse == Local_SolnBlk[Level_Coarse][nb].ICu &&
		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeE[j_coarse] == BC_OUTFLOW_SUBSONIC) ||
		(j_coarse == Local_SolnBlk[Level_Coarse][nb].JCl &&
		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeS[i_coarse] == BC_OUTFLOW_SUBSONIC) ||
		(j_coarse == Local_SolnBlk[Level_Coarse][nb].JCu &&
		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeN[i_coarse] == BC_OUTFLOW_SUBSONIC) ||
		(i_coarse == Local_SolnBlk[Level_Coarse][nb].ICl &&
		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeW[j_coarse] == BC_FIXED_PRESSURE) ||
		(i_coarse == Local_SolnBlk[Level_Coarse][nb].ICu &&
		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeE[j_coarse] == BC_FIXED_PRESSURE) ||
		(j_coarse == Local_SolnBlk[Level_Coarse][nb].JCl &&
		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeS[i_coarse] == BC_FIXED_PRESSURE) ||
		(j_coarse == Local_SolnBlk[Level_Coarse][nb].JCu &&
		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeN[i_coarse] == BC_FIXED_PRESSURE) ||
		(i_coarse == Local_SolnBlk[Level_Coarse][nb].ICl &&
		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeW[j_coarse] == BC_FIXED_TEMP_WALL) ||
		(i_coarse == Local_SolnBlk[Level_Coarse][nb].ICu &&
		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeE[j_coarse] == BC_FIXED_TEMP_WALL) ||
		(j_coarse == Local_SolnBlk[Level_Coarse][nb].JCl &&
		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeS[i_coarse] == BC_FIXED_TEMP_WALL) ||
		(j_coarse == Local_SolnBlk[Level_Coarse][nb].JCu &&
		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeN[i_coarse] == BC_FIXED_TEMP_WALL) ||
		(i_coarse == Local_SolnBlk[Level_Coarse][nb].ICl &&
		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeW[j_coarse] == BC_ADIABATIC_WALL) ||
		(i_coarse == Local_SolnBlk[Level_Coarse][nb].ICu &&
		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeE[j_coarse] == BC_ADIABATIC_WALL) ||
		(j_coarse == Local_SolnBlk[Level_Coarse][nb].JCl &&
		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeS[i_coarse] == BC_ADIABATIC_WALL) ||
		(j_coarse == Local_SolnBlk[Level_Coarse][nb].JCu &&
		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeN[i_coarse] == BC_ADIABATIC_WALL)) {

	      // Determine the (i,j) index of the SW fine grid cell.
	      i_fine = 2*(i_coarse-Nghost)+Nghost;
	      j_fine = 2*(j_coarse-Nghost)+Nghost;

	      // Perform the prolongation using injection.
	      Uchange = Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse];
	      Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine] += Uchange;
	      Local_SolnBlk[Level_Fine][nb].U[i_fine+1][j_fine] += Uchange;
	      Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine+1] += Uchange;
	      Local_SolnBlk[Level_Fine][nb].U[i_fine+1][j_fine+1] += Uchange;

	    } else {

 	      // Set injection flag if prolongation by injection is desired.
 	      if (IP->Multigrid_IP.Prolong_Using_Injection) injection_for_this_coarse_cell = true;

 	      // If the bilinear interpolation is used as the prolongation
 	      // operator and a valid set of interpolation points can not be
 	      // found for one or more of the fine grid cells then switch
 	      // the prolongation operator to injection for all four fine
 	      // grid cells.  This is facilitated by the use of a goto
 	      // statement where 'restart' is used as the goto label.
 	    restart: ;

 	      /**********************************************************
 	       * Update the cell-centred solution state for the south-  *
 	       * west fine grid cell.                                   *
 	       **********************************************************/

 	      // Determine the (i,j) index of the SW fine grid cell.
 	      i_fine = 2*(i_coarse-Nghost)+Nghost;
 	      j_fine = 2*(j_coarse-Nghost)+Nghost;

	      // Attempt to use bilinear interpolation if the injection
 	      // flag has not been set true.
 	      if (!injection_for_this_coarse_cell) {

 		// Attempt to use the SW coarse grid solution states as
 		// the interpolants ([i-1,i],[j-1,j]).
 		error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse-1],
 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse-1].Xc,
 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse],
 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse].Xc,
 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse-1],
 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse-1].Xc,
 						    Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
 						    Uchange);
 		// If fine grid cell-centre falls outside of the SW coarse
 		// grid interpolants then try the NW coarse grid solution
 		// states as the interpolants ([i-1,i],[j,j+1]).
 		if (error_flag) {
 		  error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse].Xc,
 						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse+1],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse+1].Xc,
 						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse+1],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse+1].Xc,
 						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
 						      Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
 						      Uchange);
 		}
 		// If fine grid cell-centre falls outside of the NW coarse
 		// grid interpolants then try the SE coarse grid solution
 		// states as the interpolants ([i,i+1],[j-1,j]).
 		if (error_flag) {
 		  error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse-1],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse-1].Xc,
 						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
 						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse].Xc,
 						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse-1],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse-1].Xc,
 						      Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
 						      Uchange);
 		}
 		// If fine grid cell-centre falls outside of the SE coarse
 		// grid interpolants then try the NE coarse grid solution
 		// states as the interpolants ([i,i+1],[j,j+1]).
 		if (error_flag) {
 		  error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
 						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse+1],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse+1].Xc,
 						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse+1],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse+1].Xc,
 						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse].Xc,
 						      Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
 						      Uchange);
 		}
 		// If none of the four sets of possible interpolants 
 		// provide a suitable set of interpolants then restart
 		// the prolongation associated with the current coarse
 		// grid cell using injection.
 		if (error_flag) {
 		  injection_for_this_coarse_cell = true;
 		  goto restart;
 		}

 	      } else {

 		// Perform the prolongation using injection.
 		Uchange = Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse];

 	      }

 	      // Update the conservative solution state in the SW fine
 	      // grid cell.  If required, use residual reduction if the
 	      // update creates a unphysical state (negative density or 
 	      // pressure).  Residual reduction is performed by halving
 	      // the residual update.
 	      if (IP->Multigrid_IP.Update_Stability_Switch) {
 		Unew = Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine] + Uchange;
 		if (Unew.Unphysical_Properties()) {
 		  residual_reduction_number = 1;
 		  while (residual_reduction_number <= IP->Multigrid_IP.Maximum_Number_of_Update_Reductions) {
 		    Uchange *= HALF;
 		    Unew = Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine] + Uchange;
 		    if (!Unew.Unphysical_Properties()) break;
 		    residual_reduction_number++;
 		  }
 		  if (residual_reduction_number > IP->Multigrid_IP.Maximum_Number_of_Update_Reductions) Uchange *= ZERO;
 		}
 	      }
 	      Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine] += Uchange;

 	      /**********************************************************
 	       * Update the cell-centred solution state for the south-  *
 	       * east fine grid cell.                                   *
 	       **********************************************************/

 	      // Determine the (i,j) index of the SE fine grid cell.
 	      i_fine = 2*(i_coarse-Nghost)+Nghost+1;
 	      j_fine = 2*(j_coarse-Nghost)+Nghost;

 	      // Attempt to use bilinear interpolation if the injection
 	      // flag has not been set true.
 	      if (!injection_for_this_coarse_cell) {

 		// Attempt to use the SE coarse grid solution states as
 		// the interpolants ([i,i+1],[j-1,j]).
 		error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse-1],
 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse-1].Xc,
 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse],
 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse].Xc,
 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse-1],
 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse-1].Xc,
 						    Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
 						    Uchange);
 		// If fine grid cell-centre falls outside of the SE coarse
 		// grid interpolants then try the NE coarse grid solution
 		// states as the interpolants ([i,i+1],[j,j+1]).
		if (error_flag) {
 		  error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
 						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse+1],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse+1].Xc,
 						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse+1],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse+1].Xc,
 						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse].Xc,
 						      Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
 						      Uchange);
 		}
 		// If fine grid cell-centre falls outside of the NE coarse
 		// grid interpolants then try the SW coarse grid solution
 		// states as the interpolants ([i-1,i],[j-1,j]).
 		if (error_flag) {
 		  error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse-1],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse-1].Xc,
 						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse].Xc,
 						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
 						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse-1],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse-1].Xc,
 						      Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
 						      Uchange);
 		}
 		// If fine grid cell-centre falls outside of the SW coarse
 		// grid interpolants then try the N coarse grid solution
 		// states as the interpolants ([i-1,i],[j,j+1]).
 		if (error_flag) {
 		  error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse].Xc,
 						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse+1],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse+1].Xc,
 						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse+1],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse+1].Xc,
 						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
 						      Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
 						      Uchange);
 		}
 		// If none of the four sets of possible interpolants 
 		// provide a suitable set of interpolants then restart
 		// the prolongation associated with the current coarse
 		// grid cell using injection.
 		if (error_flag) {
 		  injection_for_this_coarse_cell = true;
 		  goto restart;
 		}

 	      } else {

 		// Perform the prolongation using injection.
 		Uchange = Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse];

 	      }

 	      // Update the conservative solution state in the SE fine
 	      // grid cell.  If required, use residual reduction if the
 	      // update creates a unphysical state (negative density or 
 	      // pressure).  Residual reduction is performed by halving
 	      // the residual update.
 	      if (IP->Multigrid_IP.Update_Stability_Switch) {
 		Unew = Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine] + Uchange;
 		if (Unew.Unphysical_Properties()) {
 		  residual_reduction_number = 1;
 		  while (residual_reduction_number <= IP->Multigrid_IP.Maximum_Number_of_Update_Reductions) {
 		    Uchange *= HALF;
 		    Unew = Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine] + Uchange;
 		    if (!Unew.Unphysical_Properties()) break;
 		    residual_reduction_number++;
 		  }
 		  if (residual_reduction_number > IP->Multigrid_IP.Maximum_Number_of_Update_Reductions) Uchange *= ZERO;
 		}
 	      }
 	      Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine] += Uchange;

 	      /**********************************************************
 	       * Update the cell-centred solution state for the north-  *
 	       * east fine grid cell.                                   *
 	       **********************************************************/

 	      // Determine the (i,j) index of the NE fine grid cell.
 	      i_fine = 2*(i_coarse-Nghost)+Nghost+1;
 	      j_fine = 2*(j_coarse-Nghost)+Nghost+1;


 	      // Attempt to use bilinear interpolation if the injection flag
 	      // has not been set true.
 	      if (!injection_for_this_coarse_cell) {

		// Attempt to use the NE coarse grid solution states as the
 		// interpolants ([i,i+1],[j,j+1]).
 		error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse+1],
 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse+1].Xc,
 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse+1],
 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse+1].Xc,
 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse],
 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse].Xc,
 						    Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
 						    Uchange);
 		// If fine grid cell-centre falls outside of the NE coarse
 		// grid interpolants then try the SE coarse grid solution
 		// states as the interpolants ([i,i+1],[j-1,j]).
 		if (error_flag) {
 		  error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse-1],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse-1].Xc,
 						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
 						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse].Xc,
 						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse-1],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse-1].Xc,
 						      Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
 						      Uchange);
 		}
 		// If fine grid cell-centre falls outside of the SE coarse
 		// grid interpolants then try the NW coarse grid solution
 		// states as the interpolants ([i-1,i],[j,j+1]).
 		if (error_flag) {
 		  error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse].Xc,
 						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse+1],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse+1].Xc,
 						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse+1],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse+1].Xc,
 						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
 						      Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
 						      Uchange);
 		}
		// If fine grid cell-centre falls outside of the NW coarse
 		// grid interpolants then try the SW coarse grid solution
 		// states as the interpolants ([i-1,i],[j-1,j]).
 		if (error_flag) {
 		  error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse-1],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse-1].Xc,
 						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse].Xc,
 						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
 						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse-1],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse-1].Xc,
 						      Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
 						      Uchange);
 		}
 		// If none of the four sets of possible interpolants 
 		// provide a suitable set of interpolants then restart
 		// the prolongation associated with the current coarse
 		// grid cell using injection.
 		if (error_flag) {
 		  injection_for_this_coarse_cell = true;
 		  goto restart;
 		}

 	      } else {

 		// Perform the prolongation using injection.
 		Uchange = Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse];

 	      }

 	      // Update the conservative solution state in the NE fine
 	      // grid cell.  If required, use residual reduction if the
	      // update creates a unphysical state (negative density or 
 	      // pressure).  Residual reduction is performed by halving
 	      // the residual update.
 	      if (IP->Multigrid_IP.Update_Stability_Switch) {
 		Unew = Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine] + Uchange;
 		if (Unew.Unphysical_Properties()) {
 		  residual_reduction_number = 1;
 		  while (residual_reduction_number <= IP->Multigrid_IP.Maximum_Number_of_Update_Reductions) {
 		    Uchange *= HALF;
 		    Unew = Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine] + Uchange;
 		    if (!Unew.Unphysical_Properties()) break;
 		    residual_reduction_number++;
 		  }
 		  if (residual_reduction_number > IP->Multigrid_IP.Maximum_Number_of_Update_Reductions) Uchange *= ZERO;
 		}
 	      }
 	      Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine] += Uchange;

 	      /**********************************************************
 	       * Update the cell-centred solution state for the north-  *
 	       * west fine grid cell.                                   *
 	       **********************************************************/

 	      // Determine the (i,j) index of the NW fine grid cell.
 	      i_fine = 2*(i_coarse-Nghost)+Nghost;
 	      j_fine = 2*(j_coarse-Nghost)+Nghost+1;

 	      // Attempt to use bilinear interpolation if the injection flag
 	      // has not been set true.
 	      if (!injection_for_this_coarse_cell) {

		// Attempt to use the NW coarse grid solution states as the
 		// interpolants ([i-1,i],[j,j+1]).
 		error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse],
 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse].Xc,
 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse+1],
 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse+1].Xc,
 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse+1],
 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse+1].Xc,
 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
 						    Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
 						    Uchange);
 		// If fine grid cell-centre falls outside of the NW coarse
 		// grid interpolants then try the NE coarse grid solution
 		// states as the interpolants ([i,i+1],[j,j+1]).
 		if (error_flag) {
 		  error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
 						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse+1],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse+1].Xc,
 						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse+1],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse+1].Xc,
 						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse].Xc,
 						      Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
 						      Uchange);
 		}
 		// If fine grid cell-centre falls outside of the NE coarse
 		// grid interpolants then try the SW coarse grid solution
		// states as the interpolants ([i-1,i],[j-1,j]).
 		if (error_flag) {
 		  error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse-1],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse-1].Xc,
 						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse].Xc,
 						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
 						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse-1],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse-1].Xc,
 						      Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
 						      Uchange);
 		}
 		// If fine grid cell-centre falls outside of the SW coarse
 		// grid interpolants then try the SE coarse grid solution
 		// states as the interpolants ([i,i+1],[j-1,j]).
 		if (error_flag) {
 		  error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse-1],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse-1].Xc,
 						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
 						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse].Xc,
 						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse-1],
 						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse-1].Xc,
 						      Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
 						      Uchange);
 		}
 		// If none of the four sets of possible interpolants 
 		// provide a suitable set of interpolants then restart
 		// the prolongation associated with the current coarse
 		// grid cell using injection.
 		if (error_flag) {
 		  injection_for_this_coarse_cell = true;
 		  goto restart;
 		}

 	      } else {

 		// Perform the prolongation using injection.
 		Uchange = Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse];

 	      }

 	      // Update the conservative solution state in the NW fine
 	      // grid cell.  If required, use residual reduction if the
 	      // update creates a unphysical state (negative density or 
 	      // pressure).  Residual reduction is performed by halving
 	      // the residual update.
 	      if (IP->Multigrid_IP.Update_Stability_Switch) {
 		Unew = Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine] + Uchange;
 		if (Unew.Unphysical_Properties()) {
 		  residual_reduction_number = 1;
 		  while (residual_reduction_number <= IP->Multigrid_IP.Maximum_Number_of_Update_Reductions) {
 		    Uchange *= HALF;
 		    Unew = Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine] + Uchange;
 		    if (!Unew.Unphysical_Properties()) break;
 		    residual_reduction_number++;
 		  }
 		  if (residual_reduction_number > IP->Multigrid_IP.Maximum_Number_of_Update_Reductions) Uchange *= ZERO;
 		}
 	      }
 	      Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine] += Uchange;

 	      // Reset the injection flag to false.
 	      injection_for_this_coarse_cell = false;

	    }

	  }

	}
      }

    }
  }

  // Prolongation and update of the solution block was successful.
  return 0;

}

/**********************************************************************
 * EB_FAS_Multigrid2D_Solver::Exchange_Solution_Information --        * 
 *                                                                    *
 * This routine conducts the exchange of solution information between *
 * neighbouring solution blocks for the specified grid level.         *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState,
	  class Quad_Soln_Block,
	  class Quad_Soln_Input_Parameters>
int EB_FAS_Multigrid2D_Solver<cState, pState,
			      Quad_Soln_Block,
			      Quad_Soln_Input_Parameters>::
Exchange_Solution_Information(const int &Level,
			      const int &Send_Mesh_Geometry_Only) {

  int error_flag;

  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

  if (Send_Mesh_Geometry_Only) {
    // Send only the mesh geometry information.
    error_flag = Send_All_Messages(Local_SolnBlk[Level],
				   List_of_Local_Solution_Blocks[Level],
				   NUM_COMP_VECTOR2D,
				   ON);
  } else {
    // Send only the solution information.
    error_flag = Send_All_Messages(Local_SolnBlk[Level],
				   List_of_Local_Solution_Blocks[Level],
				   IP->Wo.NumVar(),
				   OFF);
  }
  if (error_flag) {
    cout << "\n ERROR: FASMultigrid message passing error on processor "
	 << List_of_Local_Solution_Blocks[Level].ThisCPU
	 << ".\n";
    cout.flush();
  }
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return error_flag;

  // Solution information exchanged successfully.
  return 0;

}

/**********************************************************************
 * EB_FAS_Multigrid2D_Solver::Store_Current_Solution_in_uo_MG --      *
 *                                                                    *
 * This routine copies the solution state in each block into the      *
 * array uo_MG for the specified grid level.                          *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState,
	  class Quad_Soln_Block,
	  class Quad_Soln_Input_Parameters>
void EB_FAS_Multigrid2D_Solver<cState, pState,
			       Quad_Soln_Block,
			       Quad_Soln_Input_Parameters>::
Store_Current_Solution_in_uo_MG(const int &Level) {

  int Nghost;

  // Loop through each local solution block.
  for (int nb = 0; nb < List_of_Local_Solution_Blocks[Level].Nblk; nb++) {
    if (List_of_Local_Solution_Blocks[Level].Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      // Get the number of ghost cells on the current grid level.
      Nghost = Local_SolnBlk[Level][nb].Nghost;
      // Loop through and copy the solution content for each cell of the
      // current grid level.
      for (int i = Local_SolnBlk[Level][nb].ICl-Nghost; i <= Local_SolnBlk[Level][nb].ICu+Nghost; i++) {
	for (int j = Local_SolnBlk[Level][nb].JCl-Nghost; j <= Local_SolnBlk[Level][nb].JCu+Nghost; j++) {
 	  if (EBSolver[Level].Mesh[nb].cell_status[i][j] == CELL_STATUS_ACTIVE) {
	    MG_SolnBlk[Level][nb].uo_MG[i][j] = Local_SolnBlk[Level][nb].U[i][j];
	  }
	}
      }
    }
  }
  
}

/**********************************************************************
 * EB_FAS_Multigrid2D_Solver::Subtract_dUdt_from_P --                 *
 *                                                                    *
 * This routine subtracts dUdt from forcing term, P.  The forcing     *
 * term has already been assigned the restricted fine grid residuals. *
 * To complete the formation of the forcing term the residuals        *
 * evaluated from the restricted fine grid solution containted in     *
 * dUdt must be subtracted from the forcing term.                     *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState,
	  class Quad_Soln_Block,
	  class Quad_Soln_Input_Parameters>
void EB_FAS_Multigrid2D_Solver<cState,pState,
			       Quad_Soln_Block,
			       Quad_Soln_Input_Parameters>::
Subtract_dUdt_from_P(const int &Level) {

  // Loop through each local solution block.
  for (int nb = 0; nb < List_of_Local_Solution_Blocks[Level].Nblk; nb++) {
    if (List_of_Local_Solution_Blocks[Level].Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      // Loop through and subtract residuals evaluated from the
      // restricted fine grid solution, dUdt, from the forcing term, P,
      // for each cell of the current grid level.
      for (int i = Local_SolnBlk[Level][nb].ICl; i <= Local_SolnBlk[Level][nb].ICu; i++) {
	for (int j = Local_SolnBlk[Level][nb].JCl; j <= Local_SolnBlk[Level][nb].JCu; j++) {
 	  if (EBSolver[Level].Mesh[nb].cell_status[i][j] == CELL_STATUS_ACTIVE) {
	    MG_SolnBlk[Level][nb].P[i][j] -= Local_SolnBlk[Level][nb].dUdt[i][j][0];
	  }
	}
      }
    }
  }

}

/**********************************************************************
 * EB_FAS_Multigrid2D_Solver::CFL_Multigrid --                        *
 *                                                                    *
 * This routine sets the time step for each cell on the current       *
 * coarse grid level such that it is the minimum of the computed time *
 * step on the coarse grid and the time steps of the associated finer *
 * grid cells.                                                        *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState,
	  class Quad_Soln_Block,
	  class Quad_Soln_Input_Parameters>
void EB_FAS_Multigrid2D_Solver<cState, pState,
			       Quad_Soln_Block,
			       Quad_Soln_Input_Parameters>::
CFL_Multigrid(const int &Level_Coarse) {

  int i_fine, j_fine, Level_Fine, Nghost;
  double dt_NE, dt_SE, dt_NW, dt_SW, A_coarse;
  Level_Fine = Level_Coarse - 1;

  // Loop through each local solution block.
  for (int nb = 0; nb < List_of_Local_Solution_Blocks[Level_Fine].Nblk; nb++) {
    if (List_of_Local_Solution_Blocks[Level_Fine].Block[nb].used == ADAPTIVEBLOCK2D_USED) {

      // Get the number of ghost cells.
      Nghost = Local_SolnBlk[Level_Fine][nb].Nghost;

      // Loop through the coarse cells.
      for (int i_coarse = Local_SolnBlk[Level_Coarse][nb].ICl; i_coarse <= Local_SolnBlk[Level_Coarse][nb].ICu; i_coarse++) {
	for (int j_coarse = Local_SolnBlk[Level_Coarse][nb].JCl; j_coarse <= Local_SolnBlk[Level_Coarse][nb].JCu; j_coarse++) {
 	  if (EBSolver[Level_Coarse].Mesh[nb].cell_status[i_coarse][j_coarse] == CELL_STATUS_ACTIVE) {
	    // Determine the (i,j) index of the corresponding SW corner
	    // fine cell.
	    i_fine = 2*(i_coarse-Nghost)+Nghost;
	    j_fine = 2*(j_coarse-Nghost)+Nghost;
	    // Determine the area of the coarse grid cell.
	    A_coarse = Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].A;
	    // Determine the time-steps of each of the associated fine
	    // grid cells.
	    dt_NE = sqrt(A_coarse/Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][j_fine+1].A)*
                    Local_SolnBlk[Level_Fine][nb].dt[i_fine+1][j_fine+1];
	    dt_SE = sqrt(A_coarse/Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][j_fine].A)*
                    Local_SolnBlk[Level_Fine][nb].dt[i_fine+1][j_fine];
	    dt_NW = sqrt(A_coarse/Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine+1].A)*
                    Local_SolnBlk[Level_Fine][nb].dt[i_fine][j_fine+1];
	    dt_SW = sqrt(A_coarse/Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].A)*
                    Local_SolnBlk[Level_Fine][nb].dt[i_fine][j_fine];
	    // Determine the coarse grid cell time-step by choosing the 
	    // minimum of the fine grid cell time-steps.
	    Local_SolnBlk[Level_Coarse][nb].dt[i_coarse][j_coarse] = 
	            min(Local_SolnBlk[Level_Coarse][nb].dt[i_coarse][j_coarse],
		    min(dt_NE,min(dt_SE,min(dt_NW,dt_SW))));
	  }
	}
      }

    }
  }

}

/**********************************************************************
 * EB_FAS_Multigrid2D_Solver::Update_Primitive_Variables  --          *
 *                                                                    *
 * This routine updates all primitive variables, W, from the          *
 * conserved variables, U, for all solution blocks on the given grid  *
 * level.                                                             *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState,
	  class Quad_Soln_Block,
	  class Quad_Soln_Input_Parameters>
void EB_FAS_Multigrid2D_Solver<cState, pState,
			       Quad_Soln_Block,
			       Quad_Soln_Input_Parameters>::
Update_Primitive_Variables(const int &Level) {

  // Loop through each local solution block.
  for (int nb = 0; nb < List_of_Local_Solution_Blocks[Level].Nblk; nb++) {
    if (List_of_Local_Solution_Blocks[Level].Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      // Loop through and determine the primitive solution state for 
      // each (interior) cell of the current grid level.
      for (int i = Local_SolnBlk[Level][nb].ICl; i <= Local_SolnBlk[Level][nb].ICu; i++) {
	for (int j = Local_SolnBlk[Level][nb].JCl; j <= Local_SolnBlk[Level][nb].JCu; j++) {
 	  if (EBSolver[Level].Mesh[nb].cell_status[i][j] == CELL_STATUS_ACTIVE) {
	    Local_SolnBlk[Level][nb].W[i][j] = W(Local_SolnBlk[Level][nb].U[i][j]);
	  }
	}
      }
    }
  }

}

/**********************************************************************
 * EB_FAS_Multigrid2D_Solver::Evaluate_Solution_Changes --            *
 *                                                                    *
 * This routine evaluates the solution changes, U -= uo_MG, at the    *
 * specified grid level.  The solution change is thus prepared for    *
 * prolongation to the finer grid level.                              *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState,
	  class Quad_Soln_Block,
	  class Quad_Soln_Input_Parameters>
void EB_FAS_Multigrid2D_Solver<cState, pState,
			       Quad_Soln_Block,
			       Quad_Soln_Input_Parameters>::
Evaluate_Solution_Changes(const int &Level) {

  int Nghost;

  // Loop through each local solution block.
  for (int nb = 0; nb < List_of_Local_Solution_Blocks[Level].Nblk; nb++) {
    if (List_of_Local_Solution_Blocks[Level].Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      // Get the number of ghost cells on the current grid level.
      Nghost = Local_SolnBlk[Level][nb].Nghost;
      // Determine the solution change in every cell if the current grid.
      for (int i = Local_SolnBlk[Level][nb].ICl-Nghost; i <= Local_SolnBlk[Level][nb].ICu+Nghost; i++) {
	for (int j = Local_SolnBlk[Level][nb].JCl-Nghost; j <= Local_SolnBlk[Level][nb].JCu+Nghost; j++) {
 	  if (EBSolver[Level].Mesh[nb].cell_status[i][j] == CELL_STATUS_ACTIVE) {
	    Local_SolnBlk[Level][nb].U[i][j] -= MG_SolnBlk[Level][nb].uo_MG[i][j];
	  }
	}
      }
    }
  }
}

/**********************************************************************
 * EB_FAS_Multigrid2D_Solver::                                        *
 *                          dUdt_Multistage_Explicit_for_Multigrid -- *
 *                                                                    *
 * This routine evaluates the stage solution residual for a 1D array  *
 * of 2D quadrilateral multi-block solution blocks.  A variety of     *
 * multistage explicit time integration schemes, upwind finite-volume *
 * spatial discretization procedures for hyperbolic fluxes, and       *
 * centrally-weighted finite-volume spatial discretization procedures *
 * for elliptic fluxes can be used depending on the specified input   *
 * values.                                                            *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState,
	  class Quad_Soln_Block,
	  class Quad_Soln_Input_Parameters>
int EB_FAS_Multigrid2D_Solver<cState, pState,
			      Quad_Soln_Block,
			      Quad_Soln_Input_Parameters>::
dUdt_Multistage_Explicit_for_Multigrid(const int &Level,
				       const int &Top_Level,
				       const int &i_stage,
				       const double &Time) {

  int error_flag, time_integration, flux_function_type, limiter_type, k_residual;

  // Force the use of first-order spatial reconstruction and the HLLE
  // hyperbolic flux function on coarse levels if desired.
  if (Level > Top_Level && IP->Multigrid_IP.First_Order_Coarse_Mesh_Reconstruction) {
    flux_function_type = IP->i_Flux_Function;
    limiter_type = IP->i_Limiter;
    if (IP->Local_Time_Stepping != LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER) {
      IP->i_Flux_Function = FLUX_FUNCTION_HLLE;
    }
    IP->i_Limiter = LIMITER_ZERO;
  }

  // Temporarily overwrite the time integration type to the multistage
  // optimally smoothing scheme.
  time_integration = IP->i_Time_Integration;
  IP->i_Time_Integration = IP->Multigrid_IP.i_Smoothing;

  // Evaluate the solution residual for each block.
  error_flag = EBSolver[Level].dUdt_Multistage_Explicit(i_stage,Time);
  if (error_flag) return error_flag;

  // Reset the flux function, limiter, and time integration parameters.
  if (Level > Top_Level && IP->Multigrid_IP.First_Order_Coarse_Mesh_Reconstruction) {
    IP->i_Flux_Function = flux_function_type;
    IP->i_Limiter = limiter_type;
  }

  // Reset the time-integration type.
  IP->i_Time_Integration = time_integration;

  // Residuals for each quadrilateral multi-block solution has been
  // successfully calculated.
  return 0;

}

/**********************************************************************
 * EB_FAS_Multigrid2D_Solver::                                        *
 *             Apply_Boundary_Flux_Corrections_Multistage_Explicit -- *
 *                                                                    *
 * Apply flux corrections at boundaries of a 1D array of 2D           *
 * quadrilateral multi-block solution blocks to ensure that the       *
 * scheme is conservative at boundaries with mesh resolution changes. *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState,
	  class Quad_Soln_Block,
	  class Quad_Soln_Input_Parameters>
void EB_FAS_Multigrid2D_Solver<cState, pState,
			       Quad_Soln_Block,
			       Quad_Soln_Input_Parameters>::
Apply_Boundary_Flux_Corrections_Multistage_Explicit_for_Multigrid(const int &Level,
								  const int &Top_Level,
								  const int &i_stage) {

  int time_integration, flux_function_type, limiter_type;

  // Force the use of first-order spatial reconstruction and the HLLE
  // hyperbolic flux function on coarse levels if desired.
  if (Level > Top_Level && IP->Multigrid_IP.First_Order_Coarse_Mesh_Reconstruction) {
    flux_function_type = IP->i_Flux_Function;
    limiter_type = IP->i_Limiter;
    if (IP->Local_Time_Stepping != LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER) {
      IP->i_Flux_Function = FLUX_FUNCTION_HLLE;
    }
    IP->i_Limiter = LIMITER_ZERO;
  }

  // Temporarily overwrite the time integration type to the multistage
  // optimally smoothing scheme.
  time_integration = IP->i_Time_Integration;
  IP->i_Time_Integration = IP->Multigrid_IP.i_Smoothing;
        
  // Apply flux corrections at the boundaries of the solution blocks.
  Apply_Boundary_Flux_Corrections_Multistage_Explicit(Local_SolnBlk[Level],
						      List_of_Local_Solution_Blocks[Level],
						      *IP,
						      i_stage);

  // Reset the flux function, limiter, and time integration parameters.
  if (Level > Top_Level && IP->Multigrid_IP.First_Order_Coarse_Mesh_Reconstruction) {
    IP->i_Flux_Function = flux_function_type;
    IP->i_Limiter = limiter_type;
  }

  // Reset the time-integration type.
  IP->i_Time_Integration = time_integration;

}

/**********************************************************************
 * EB_FAS_Multigrid2D_Solver::                                        *
 *               Update_Solution_Multistage_Explicit_for_Multigrid -- *
 *                                                                    *
 * This routine updates the solution for a 1D array of 2D             *
 * multi-block solution blocks.  A variety of multistage explicit     *
 * time integration schemes, upwind finite-volume spatial             *
 * discretization procedure for hyperbolic fluxes, and centrally-     *
 * weighted finite-volume spatial discretization procedure for        *
 * elliptic fluxes can be used depending on the specified input       *
 * values.                                                            *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState,
	  class Quad_Soln_Block,
	  class Quad_Soln_Input_Parameters>
int EB_FAS_Multigrid2D_Solver<cState, pState,
			      Quad_Soln_Block,
			      Quad_Soln_Input_Parameters>::
Update_Solution_Multistage_Explicit_for_Multigrid(const int &Level,
						  const int &Top_Level,
						  const int &i_stage) {
  
  int error_flag, time_integration, flux_function_type, limiter_type;
  
  // Force the use of first-order spatial reconstruction and the HLLE
  // hyperbolic flux function on coarse levels if desired.
  if (Level > Top_Level && IP->Multigrid_IP.First_Order_Coarse_Mesh_Reconstruction) {
    flux_function_type = IP->i_Flux_Function;
    limiter_type = IP->i_Limiter;
    if (IP->Local_Time_Stepping != LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER) {
      IP->i_Flux_Function = FLUX_FUNCTION_HLLE;
    }
    IP->i_Limiter = LIMITER_ZERO;
  }

  // Temporarily overwrite the time integration type to the multistage
  // optimally smoothing scheme.
  time_integration = IP->i_Time_Integration;
  IP->i_Time_Integration = IP->Multigrid_IP.i_Smoothing;

  // Update the solution for each solution block.
  error_flag = EBSolver[Level].Update_Solution_Multistage_Explicit(i_stage);
  if (error_flag) return error_flag;

  // Reset the flux function, limiter, and time integration parameters.
  if (Level > Top_Level && IP->Multigrid_IP.First_Order_Coarse_Mesh_Reconstruction) {
    IP->i_Flux_Function = flux_function_type;
    IP->i_Limiter = limiter_type;
  }

  // Reset the time-integration type.
  IP->i_Time_Integration = time_integration;

  // Quadrilateral multi-block solution blocks have been updated
  // successfully.
  return 0;

}

/**********************************************************************
 * EB_FAS_Multigrid2D_Solver::Apply_the_FAS_Multigrid_Forcing_Term -- *
 *                                                                    *
 * This routine applies the FAS multigrid defect correction forcing   *
 * term, P, for a 1D array of 2D quadrilateral solution blocks.       *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState,
	  class Quad_Soln_Block,
	  class Quad_Soln_Input_Parameters>
int EB_FAS_Multigrid2D_Solver<cState, pState,
			      Quad_Soln_Block,
			      Quad_Soln_Input_Parameters>::
Apply_the_FAS_Multigrid_Forcing_Term(const int &Level,
				     const int &Top_Level,
				     const int &i_stage) {

  int k_residual, time_integration;

  // Do not apply the forcing term if the current level is finer or
  // equal to the current top level.  Exit immediately.
  if (Level <= Top_Level) return 0;

  // Temporarily overwrite the time integration type to the multistage
  // optimally smoothing scheme.
  time_integration = IP->i_Time_Integration;
  IP->i_Time_Integration = IP->Multigrid_IP.i_Smoothing;

  // Apply the forcing term.
  switch(IP->i_Time_Integration) {
  case TIME_STEPPING_EXPLICIT_EULER :
    k_residual = 0;
    break;
  case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
    k_residual = 0;
    break;
  case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA :
    k_residual = 0;
    if (IP->N_Stage == 4) {
      if (i_stage == 4) {
	k_residual = 0;
      } else {
	k_residual = i_stage - 1;
      }
    }
    break;
  case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
    k_residual = 0;
    break;
  default:
    k_residual = 0;
    break;
  }
  // Loop through each solution block on this grid level.
  for (int nb = 0; nb < List_of_Local_Solution_Blocks[Level].Nblk; nb++) {
    if (List_of_Local_Solution_Blocks[Level].Block[nb].used == ADAPTIVEBLOCK2D_USED) {      
      // Add the forcing term at each cell on this grid level.
      for (int j = Local_SolnBlk[Level][nb].JCl; j <= Local_SolnBlk[Level][nb].JCu; j++) {
	for (int i = Local_SolnBlk[Level][nb].ICl; i <= Local_SolnBlk[Level][nb].ICu; i++) {
 	  if (EBSolver[Level].Mesh[nb].cell_status[i][j] == CELL_STATUS_ACTIVE) {
	    Local_SolnBlk[Level][nb].dUdt[i][j][k_residual] += (IP->CFL_Number*Local_SolnBlk[Level][nb].dt[i][j])*
  	                                                        MG_SolnBlk[Level][nb].P[i][j];
	  }
	}
      }
    }
  }

  // Reset the time-integration type.
  IP->i_Time_Integration = time_integration;

  // The forcing term has been applied successfully.
  return 0;

}

/**********************************************************************
 * EB_FAS_Multigrid2D_Solver::Residual_Evaluation_for_Multigrid --    *
 *                                                                    *
 * This routine evaluates the residual for a 1D array of solution     *
 * blocks given the solution U using an upwind finite-volume spatial  *
 * discretization procedure for hyperbolic fluxes and a centrally-    *
 * weighted finite-volume spatial discretization procedure for the    *
 * elliptic fluxes.                                                   *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState,
	  class Quad_Soln_Block,
	  class Quad_Soln_Input_Parameters>
int EB_FAS_Multigrid2D_Solver<cState, pState,
			      Quad_Soln_Block,
			      Quad_Soln_Input_Parameters>::
Residual_Evaluation_for_Multigrid(const int &Level,
				  const int &Top_Level,
				  const double &dt,
				  const double &Time,
				  const int &apply_forcing_term_flag) {

  int error_flag, flux_function_type, limiter_type;

  // Force the use of first-order spatial reconstruction and the HLLE
  // hyperbolic flux function on coarse levels if desired.
  if (Level > Top_Level && IP->Multigrid_IP.First_Order_Coarse_Mesh_Reconstruction) {
    flux_function_type = IP->i_Flux_Function;
    limiter_type = IP->i_Limiter;
    if (IP->Local_Time_Stepping != LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER) {
      IP->i_Flux_Function = FLUX_FUNCTION_HLLE;
    }
    IP->i_Limiter = LIMITER_ZERO;
  }

  // Evaluate the residual for each solution block.
  for (int nb = 0; nb < List_of_Local_Solution_Blocks[Level].Nblk; nb++) {
    if (List_of_Local_Solution_Blocks[Level].Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      error_flag = EBSolver[Level].dUdt_Residual_Evaluation(Time);
      if (error_flag) return error_flag;
    }
  }

  // Reset the flux function, limiter, and time integration parameters.
  if (Level > Top_Level && IP->Multigrid_IP.First_Order_Coarse_Mesh_Reconstruction) {
    IP->i_Flux_Function = flux_function_type;
    IP->i_Limiter = limiter_type;
  }

  // Add dual-time-stepping physical-time source term if required.
  if (IP->i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING) {
    for (int nb = 0; nb < List_of_Local_Solution_Blocks[Level].Nblk; nb++) {
      if (List_of_Local_Solution_Blocks[Level].Block[nb].used == ADAPTIVEBLOCK2D_USED) {
	for (int j = Local_SolnBlk[Level][nb].JCl; j <= Local_SolnBlk[Level][nb].JCu; j++) {
	  for (int i = Local_SolnBlk[Level][nb].ICl; i <= Local_SolnBlk[Level][nb].ICu; i++) {
	    if (EBSolver[Level].Mesh[nb].cell_status[i][j] == CELL_STATUS_ACTIVE) {
	      if (IP->Multigrid_IP.i_Physical_Time_Integration == TIME_STEPPING_IMPLICIT_EULER) {
		Local_SolnBlk[Level][nb].dUdt[i][j][0] -= (Local_SolnBlk[Level][nb].U[i][j] -
							   DTS_SolnBlk[Level][nb].Un[i][j])/
		                                          (IP->Multigrid_IP.Physical_Time_CFL_Number*dt);
	      } else if (IP->Multigrid_IP.i_Physical_Time_Integration == TIME_STEPPING_IMPLICIT_SECOND_ORDER_BACKWARD) {
		Local_SolnBlk[Level][nb].dUdt[i][j][0] -= (3.0*Local_SolnBlk[Level][nb].U[i][j] -
							   4.0*DTS_SolnBlk[Level][nb].Un[i][j] +
							   DTS_SolnBlk[Level][nb].Uo[i][j])/
		                                          (TWO*IP->Multigrid_IP.Physical_Time_CFL_Number*dt);
	      }
	    }
	  }
	}
      }
    }
  }

  // Send the boundary flux corrections at block interfaces with
  // resolution changes.
  error_flag = Send_Conservative_Flux_Corrections(Local_SolnBlk[Level],
						  List_of_Local_Solution_Blocks[Level],
						  Local_SolnBlk[Level][0].NumVar());
  if (error_flag) return error_flag;

  // Apply the boundary flux corrections to ensure that method is
  // conservative.
  Apply_Boundary_Flux_Corrections(Local_SolnBlk[Level],
				  List_of_Local_Solution_Blocks[Level]);

  // Add the forcing term if required.
  if (apply_forcing_term_flag) {
    for (int nb = 0; nb < List_of_Local_Solution_Blocks[Level].Nblk; nb++) {
      if (List_of_Local_Solution_Blocks[Level].Block[nb].used == ADAPTIVEBLOCK2D_USED) {
	for (int j = Local_SolnBlk[Level][nb].JCl; j <= Local_SolnBlk[Level][nb].JCu; j++) {
	  for (int i = Local_SolnBlk[Level][nb].ICl; i <= Local_SolnBlk[Level][nb].ICu; i++) {
	    if (EBSolver[Level].Mesh[nb].cell_status[i][j] == CELL_STATUS_ACTIVE) {
	      Local_SolnBlk[Level][nb].dUdt[i][j][0] += MG_SolnBlk[Level][nb].P[i][j];
	    }
	  }
	}
      }
    }
  }

  // The residuals have been successfully calculated.
  return 0;

}

/**********************************************************************
 * EB_FAS_Multigrid2D_Solver::Smooth --                               *
 *                                                                    *
 * This routine implements the multigrid smoothing for a 1D array of  *
 * 2D quadrilateral multi-block solution blocks at the specified grid *
 * level.                                                             *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState,
	  class Quad_Soln_Block,
	  class Quad_Soln_Input_Parameters>
int EB_FAS_Multigrid2D_Solver<cState, pState,
			      Quad_Soln_Block,
			      Quad_Soln_Input_Parameters>::
Smooth(const int &Level,
       const int &Top_Level,
       const double &dt,
       const double &Time) {

  int error_flag;
  double dTime;

  // Determine the local (and minimum) time-step for each cell of the 
  // current grid level.
#ifdef _MG_DOUT_CYCLE_STAGE_
  for (int l = 0; l < Level; l++) dout << "   ";
  dout << "~~> Determine the CFL number." << endl;
#endif
  dTime = EBSolver[Level].CFL(Time);
  // Set each cell to the global minimum time-step if required.
//   if (!IP->Local_Time_Stepping) {
//      dTime = CFFC_Minimum_MPI(dTime);
//      Set_Global_TimeStep(Local_SolnBlk[Level],
// 			 List_of_Local_Solution_Blocks[Level],
// 			 dTime);
//   }

  // If the current grid level is not the top grid level then ensure
  // that the coarse grid time-step is less than or equal to the 
  // time-steps of the associated finer grid level.
#ifdef _MG_DOUT_CYCLE_STAGE_
  for (int l = 0; l < Level; l++) dout << "   ";
  dout << "~~> Compare coarse grid time-step to fine grid time-steps." << endl;
#endif
  if (Level != Top_Level) CFL_Multigrid(Level);

  // Smooth the solution using N stage multistage time stepping scheme.
  for (int i_stage = 1; i_stage <= IP->N_Stage; i_stage++) {

    // Step 1. Exchange solution information between neighbouring blocks.
#ifdef _MG_DOUT_CYCLE_STAGE_
    for (int l = 0; l < Level; l++) dout << "   ";
    dout << "   -> Exchange_Solution_Information." << endl;
#endif
    error_flag = Exchange_Solution_Information(Level,
					       OFF);
    if (error_flag) {
      cout << "\n ERROR: Message passing error on processor "
	   << List_of_Local_Solution_Blocks[Level].ThisCPU
	   << "." << endl;
    }
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;

    // Step 2. Apply boundary conditions for stage.
#ifdef _MG_DOUT_CYCLE_STAGE_
    for (int l = 0; l < Level; l++) dout << "   ";
    dout << "   -> Apply boundary conditions." << endl;
#endif
    if (Level == Top_Level || IP->Multigrid_IP.Apply_Coarse_Mesh_Boundary_Conditions) {
      EBSolver[Level].Boundary_Conditions(Time);
    }

    // Step 3. Determine solution residuals for stage.
#ifdef _MG_DOUT_CYCLE_STAGE_
    for (int l = 0; l < Level; l++) dout << "   ";
    dout << "   -> dUdt_Multistage_Explicit_for_Multigrid." << endl;
#endif

    error_flag = dUdt_Multistage_Explicit_for_Multigrid(Level,
							Top_Level,
							i_stage,
							Time);
    if (error_flag) {
      cout << "\n ERROR: Residual calculation error on processor "
	   << List_of_Local_Solution_Blocks[Level].ThisCPU
	   << "." << endl;
    }
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;

    // Step 4. Add the physical-time source-term if required.
#ifdef _MG_DOUT_CYCLE_STAGE_
    for (int l = 0; l < Level; l++) dout << "   ";
    dout << "   -> dUdtau_Multistage_Explicit." << endl;
#endif
    error_flag = dUdtau_Multistage_Explicit(Level,
					    i_stage,
					    dt);
    if (error_flag) {
      cout << "\n ERROR: Physical-time source-term evaluation on processor "
	   << List_of_Local_Solution_Blocks[Level].ThisCPU << "." << endl;
    }
    if (error_flag) return error_flag;

    // Step 5. Send boundary flux corrections at block interfaces with
    // resolution changes.
#ifdef _MG_DOUT_CYCLE_STAGE_
    for (int l = 0; l < Level; l++) dout << "   ";
    dout << "   -> Send_Conservative_Flux_Corrections." << endl;
#endif
    error_flag = Send_Conservative_Flux_Corrections(Local_SolnBlk[Level],
						    List_of_Local_Solution_Blocks[Level],
						    Local_SolnBlk[Level][0].NumVar());
    if (error_flag) {
      cout << "\n ERROR: Flux correction message passing error on processor "
	   << List_of_Local_Solution_Blocks[Level].ThisCPU << "." << endl;
    }
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;

    // Step 6. Apply boundary flux corrections to ensure that method is
    // conservative.
#ifdef _MG_DOUT_CYCLE_STAGE_
    for (int l = 0; l < Level; l++) dout << "   ";
    dout << "   -> Apply_Boundary_Flux_Corrections_Multistage_Explicit_for_Multigrid." << endl;
#endif
    Apply_Boundary_Flux_Corrections_Multistage_Explicit_for_Multigrid(Level,
								      Top_Level,
								      i_stage);

    // Step 7. Apply the defect correction forcing term, if necessary.
#ifdef _MG_DOUT_CYCLE_STAGE_
    for (int l = 0; l < Level; l++) dout << "   ";
    dout << "   -> Apply_the_FAS_Multigrid_Forcing_Term." << endl;
#endif
    error_flag = Apply_the_FAS_Multigrid_Forcing_Term(Level,
						      Top_Level,
						      i_stage);
    if (error_flag) {
      cout << "\n ERROR: Forcing term application error on processor "
	   << List_of_Local_Solution_Blocks[Level].ThisCPU
	   << "." << endl;
    }
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;

    // Step 8. Smooth the solution residual using implicit residual smoothing.
#ifdef _MG_DOUT_CYCLE_STAGE_
    for (int l = 0; l < Level; l++) dout << "   ";
    dout << "   -> Residual_Smoothing." << endl;
#endif
    if (IP->Residual_Smoothing) EBSolver[Level].Residual_Smoothing(i_stage);

    // Step 8.5. Add the physical-time source-term if required.
#ifdef _MG_COUT_CYCLE_STAGE_SMOOTH
    for (int l = 0; l < Level; l++) cout << "   ";
    cout << "   -> Apply_Melson_Time_Step." << endl; cout.flush();
#endif
    error_flag = Apply_Melson_Time_Step(Level,
					i_stage,
					dt);
    if (error_flag) {
      cout << "\n FASMultigrid2D ERROR: Physical-time source-term evaluation on processor "
	   << List_of_Local_Solution_Blocks[Level].ThisCPU << ".\n";
      cout.flush();
    }
    if (error_flag) return error_flag;

    // Step 9. Update solution for the current stage.
#ifdef _MG_DOUT_CYCLE_STAGE_
    for (int l = 0; l < Level; l++) dout << "   ";
    dout << "   -> Update_Solution_Multistage_Explicit_for_Multigrid." << endl;
#endif
    error_flag = Update_Solution_Multistage_Explicit_for_Multigrid(Level,
								   Top_Level,
								   i_stage);
    if (error_flag) {
      cout << "\n ERROR: Solution update error on processor "
	   << List_of_Local_Solution_Blocks[Level].ThisCPU
	   << "." << endl;
    }
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;

  }

  // Exchange solution information between neighbouring blocks
  // after smoothing is complete.
#ifdef _MG_DOUT_CYCLE_STAGE_
    for (int l = 0; l < Level; l++) dout << "   ";
    dout << "~~> Post Exchange_Solution_Information." << endl;
#endif
  error_flag = Exchange_Solution_Information(Level,
         				     OFF);
  if (error_flag) {
    cout << "\n ERROR: Message passing error on processor "
	 << List_of_Local_Solution_Blocks[Level].ThisCPU
	 << "." << endl;
  }
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return error_flag;

  // Apply boundary conditions after smoothing is complete.
#ifdef _MG_DOUT_CYCLE_STAGE_
  for (int l = 0; l < Level; l++) dout << "   ";
  dout << "~~> Post BCs." << endl;
#endif
  if (Level == Top_Level || IP->Multigrid_IP.Apply_Coarse_Mesh_Boundary_Conditions) {
    EBSolver[Level].Boundary_Conditions(Time);
  }

  // Determine solution residuals for the entire stage.
#ifdef _MG_DOUT_CYCLE_STAGE_
  for (int l = 0; l < Level; l++) dout << "   ";
  dout << "~~> Residual_Evaluation_for_Multigrid." << endl;
#endif
  error_flag = Residual_Evaluation_for_Multigrid(Level,
						 Top_Level,
						 dt,
						 Time,
						 Level-Top_Level);
  if (error_flag) {
    cout << "\n ERROR: Residual calculation error on processor "
	 << List_of_Local_Solution_Blocks[Level].ThisCPU << "." << endl;
  }
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return error_flag;

  // Multigrid smoothing applied successfully.
  return 0;

}

/**********************************************************************
 * EB_FAS_Multigrid2D_Solver::Coarse_Grid_Correction --               *
 *                                                                    *
 * This routine conducts the FAS multigrid solution of the desired    *
 * system of hyperbolic/elliptic conservation laws for the multigrid  *
 * cycle and number of levels specified in the input parameters.      *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState,
	  class Quad_Soln_Block,
	  class Quad_Soln_Input_Parameters>
int EB_FAS_Multigrid2D_Solver<cState, pState,
			      Quad_Soln_Block,
			      Quad_Soln_Input_Parameters>::
Coarse_Grid_Correction(const int &Top_Level,
		       const int &Current_Level,
		       const double &dt,
		       const double &Time) {

  int error_flag, N_Smooths, m_Cycle;
  bool flag_for_First_Time_Restriction = true;

#ifdef _MG_DOUT_CYCLE_STAGE_
  for (int l = 0; l< Current_Level; l++) dout << "   ";
  dout << "Entering Level " << Current_Level << endl;
#endif

  // If the current grid level is the coarsest grid then smooth the
  // solution the number of time specified in the input parameters.
  // Otherwise conduct the FAS multigrid method for the cycle specified
  // in the input parameters.

  if (Current_Level == IP->Multigrid_IP.Levels-1) {

    // The current grid level is the coarsest grid.

    // Update the first time restriction flag.
    flag_for_First_Time_Restriction = false;

    // Store the unsmoothed solution for future use.
    Store_Current_Solution_in_uo_MG(Current_Level);

#ifdef _MG_DOUT_CYCLE_STAGE_
    for (int l = 0; l < Current_Level; l++) dout << "   ";
    dout << "Coarse Smoothing" << endl;
#endif

    // Smooth on the coarsest grid the number of times specifed in the
    // input parameters.
    for (int nu = 1; nu <= IP->Multigrid_IP.Number_of_Smooths_on_Coarsest_Level; nu++) {
      error_flag = Smooth(Current_Level,
			  Top_Level,
			  dt,
			  Time);
      if (error_flag) return error_flag;
    }
    
  } else {

    // The current level is not the coarsest level.  Perform the FAS
    // multigrid cycle specified in the input parameters.
    if (IP->Multigrid_IP.i_Cycle == MULTIGRID_V_CYCLE ||
	Current_Level == FINEST_LEVEL) {
      m_Cycle = 1;
    } else if (IP->Multigrid_IP.i_Cycle == MULTIGRID_W_CYCLE) {
      m_Cycle = 2;
    } else {
      return 3101;
    }

    // Recursive implementation of the V and W cycles.
    for (int m = 1; m <= m_Cycle; m++) {

#ifdef _MG_DOUT_CYCLE_STAGE_
      for (int l = 0; l < Current_Level; l++) dout << "   ";
      dout << "m = " << m << endl;
#endif

      // If the current level is the top level then set the first time
      // refinement flag to true.
      if (Current_Level == Top_Level) {
	flag_for_First_Time_Restriction = true;
      }

      // If the first time refinement flag is true then store the
      // unsmoothed solution for later use.
      if (flag_for_First_Time_Restriction == true) {
#ifdef _MG_DOUT_CYCLE_STAGE_
	for (int l = 0; l < Current_Level; l++) dout << "   ";
	dout << "Storing uo_MG" << endl;
#endif
	Store_Current_Solution_in_uo_MG(Current_Level);	
      }

      // Pre-smooth the solution on the current grid the number of times
      // specified by the input parameters depending on the current grid
      // level.

      // Determine the number of pre-smoothing iterations.
      N_Smooths = (Current_Level == Top_Level) ?
	          IP->Multigrid_IP.Number_of_Smooths_on_Finest_Level :
	          IP->Multigrid_IP.Number_of_Pre_Smooths;

#ifdef _MG_DOUT_CYCLE_STAGE_
      for (int l = 0; l < Current_Level; l++) dout << "   ";
      dout << "Presmoothing " << N_Smooths << " times" << endl;
#endif

      // Conduct the pre-smoothing for the specified number of
      // iterations.
      for (int nu1 = 0; nu1 < N_Smooths; nu1++) {
	error_flag = Smooth(Current_Level,
			    Top_Level,
			    dt,
			    Time);
	if (error_flag) return error_flag;
      }

#ifdef _MG_DOUT_CYCLE_STAGE_
      for (int l = 0; l < Current_Level; l++) dout << "   ";
      dout << "Restrict Solution from " << Current_Level << " to level " << Current_Level+1 << endl;
#endif

      // Restrict the solution to the coarse grid.
      Restrict_Solution_Blocks(Current_Level);

#ifdef _MG_DOUT_CYCLE_STAGE_
      for (int l = 0; l < Current_Level; l++) dout << "   ";
      dout << "Update the primitive variables on level " << Current_Level+1 << endl;
#endif

      // Update the primitive variables on the coarse grid.
      Update_Primitive_Variables(Current_Level+1);

#ifdef _MG_DOUT_CYCLE_STAGE_
      for (int l = 0; l < Current_Level; l++) dout << "   ";
      dout << "Message-Passing, post Restriction, on level " << Current_Level+1 << endl;
#endif

      // Update the ghostcell information to ensure that the solution
      // is consistent on each block on the coarse grid level.
      error_flag = Exchange_Solution_Information(Current_Level+1,
						 OFF);
      if (error_flag) return error_flag;

#ifdef _MG_DOUT_CYCLE_STAGE_
      for (int l = 0; l < Current_Level; l++) dout << "   ";
      dout << "Apply BCs on level " << Current_Level+1 << endl;
#endif

      // Apply boundary conditions on the coarse grid if required.
      if (Current_Level+1 == Top_Level ||
	  IP->Multigrid_IP.Apply_Coarse_Mesh_Boundary_Conditions) {
	EBSolver[Current_Level+1].Boundary_Conditions(Time);
      }

      // Evaluate the forcing term, P, if required.  The forcing term is
      // formed by the difference of the restricted fine grid residual
      // and the residual on the coarse grid found by the restricted
      // fine grid solution.
      if (IP->Multigrid_IP.Defect_Correction) {
#ifdef _MG_DOUT_CYCLE_STAGE_
	for (int l = 0; l < Current_Level; l++) dout << "   ";
	dout << "Form defect corrections, P " << endl;
#endif
 	// Restrict residual from fine to coarse grid.
 	Restrict_Residuals(Current_Level);
 	// Evaluate the residuals of restricted fine grid solution.
 	error_flag = Residual_Evaluation_for_Multigrid(Current_Level+1,
 						       Top_Level,
						       dt,
						       Time,
 						       0);
 	if (error_flag) return error_flag;
 	// Form the forcing term by subtracting the residual evaluated 
 	// from the restricted fine grid solution from P which already
 	// contains the restricted fine grid residuals.
 	Subtract_dUdt_from_P(Current_Level+1);
      }

#ifdef _MG_DOUT_CYCLE_STAGE_
      for (int l = 0; l < Current_Level; l++) dout << "   ";
      dout << "Heading Down from Level " << Current_Level << endl;
#endif

      // Call the coarse grid correction scheme (recursive) to solve
      // the coarse grid level problem.      
      error_flag = Coarse_Grid_Correction(Top_Level,
 					  Current_Level+1,
					  dt,
					  Time);
      if (error_flag) return error_flag;

#ifdef _MG_DOUT_CYCLE_STAGE_
      for (int l = 0; l < Current_Level; l++) dout << "   ";
      dout << "Coming back up to Level " << Current_Level << endl;
#endif

#ifdef _MG_DOUT_CYCLE_STAGE_
      for (int l = 0; l < Current_Level; l++) dout << "   ";
      dout << "Prolong & Update " << endl;
#endif

      // Prolong the solution changes from the coarse grid back to the
      // current grid.
      error_flag = Prolong_and_Update_Solution_Blocks(Current_Level+1);
      if (error_flag) return error_flag;

      // Update Primitive variables on the current grid level.
      Update_Primitive_Variables(Current_Level);

      // Update the ghostcell information to ensure that the solution
      // is consistent on each block on the current grid level.
      error_flag = Exchange_Solution_Information(Current_Level,
						 OFF);
      if (error_flag) return error_flag;

      // Apply boundary conditions on the current grid level.
      if (Current_Level == Top_Level || IP->Multigrid_IP.Apply_Coarse_Mesh_Boundary_Conditions) {
	EBSolver[Current_Level].Boundary_Conditions(Time);
      }

      // Post-smooth the solution on the current grid the number of
      // times specified by the input parameters depending on the
      // current grid level.

      // Determine the number of post-smoothing iterations.
      N_Smooths = (Current_Level == Top_Level) ? 0 :
		  IP->Multigrid_IP.Number_of_Post_Smooths;

#ifdef _MG_DOUT_CYCLE_STAGE_
      for (int l = 0; l < Current_Level; l++) dout << "   ";
      dout << "Postsmooth " << N_Smooths << " times" << endl;
#endif

      // Conduct the post-smoothing for the specified number of
      // iterations.
      for (int nu2 = 0; nu2 < N_Smooths; nu2++) {
	error_flag = Smooth(Current_Level,
			    Top_Level,
			    dt,
			    Time);
	if (error_flag) return error_flag;
      }

    }

  }

  // Calculate the change in the solution if the current grid level
  // is not the top grid level so that the solution change can be 
  // prolonged to the finer grid level.
  if (Current_Level != Top_Level) {
#ifdef _MG_DOUT_CYCLE_STAGE_
    for (int l = 0; l < Current_Level; l++) dout << "   ";
    dout << "Evaluate Solution changes " << endl;
#endif
    Evaluate_Solution_Changes(Current_Level);
  }

  // Coarse grid correction scheme successfully computed.
  return 0;

}

/**********************************************************************
 * EB_FAS_Multigrid2D_Solver::Execute --                              *
 *                                                                    *
 * This routine executes the FAS multigrid solution of the desired    *
 * system of hyperbolic/elliptic conservation laws.                   *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState,
	  class Quad_Soln_Block,
	  class Quad_Soln_Input_Parameters>
int EB_FAS_Multigrid2D_Solver<cState, pState,
			      Quad_Soln_Block,
			      Quad_Soln_Input_Parameters>::
Execute(int &batch_flag,
	int &number_of_time_steps,
	int &evolution_counter,
	int &levelset_iterations,
	double &Time,
	double &levelset_Time,
	CPUTime &processor_cpu_time,
	CPUTime &total_cpu_time,
	ofstream &residual_file) {

  int error_flag, command_flag, first_step, line_number, limiter_freezing;
  int initial_top_level;
  double residual_l1_norm, residual_l2_norm, residual_max_norm;
  double initial_residual_l2_norm = -1.0, ratio_residual_l2_norm = -1.0;
  unsigned long cycles_for_this_level, max_cycles_for_this_level;
#ifdef _GNU_GCC_V3
  max_cycles_for_this_level = numeric_limits<unsigned long>::max();
#else
  max_cycles_for_this_level = 10000000;
#endif

  // Send solution information between neighbouring blocks to complete
  // prescription of initial data.
  for (int level = 1; level < IP->Multigrid_IP.Levels; level++) {
    error_flag = Exchange_Solution_Information(level,ON);
    if (error_flag) break;
  }
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return error_flag;

  // Apply the initial conditions on all coarse grid levels.
  for (int level = 1; level < IP->Multigrid_IP.Levels; level++) {
    EBSolver[level].Initialize_Adjustment_Grids();
    error_flag = EBSolver[level].Mesh_Adjustment(ON,OFF);
    if (error_flag) break;
    error_flag = EBSolver[level].Compute_Interface_Velocity_Function(Time);
    if (error_flag) break;
    ICs(Local_SolnBlk[level],
	List_of_Local_Solution_Blocks[level],
	*IP);
    if (IP->i_ICs == IC_RESTART) {
      Restrict_Solution_Blocks(level-1);
      Update_Primitive_Variables(level);
    }
    if (IP->Multigrid_IP.Apply_Coarse_Mesh_Boundary_Conditions)
      Restrict_Boundary_Ref_States(level-1);
  }
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return error_flag;

  // Send solution information between neighbouring blocks to complete
  // prescription of initial data.
  for (int level = 1; level < IP->Multigrid_IP.Levels; level++) {
    error_flag = Exchange_Solution_Information(level,OFF);
    if (error_flag) break;
  }
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return error_flag;

  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

  // Open residual file.  
  first_step = 1;
  limiter_freezing = OFF;
  if (CFFC_Primary_MPI_Processor()) {
    error_flag = Open_Progress_File(residual_file,
				    IP->Output_File_Name,
				    number_of_time_steps);
    if (error_flag) {
      cout << "\n ERROR: Unable to open residual file for calculation." << endl;
    }
  }
  CFFC_Broadcast_MPI(&error_flag,1);
  if (error_flag) return error_flag;

  // Reset the CPU time.
  if (IP->i_ICs != IC_RESTART) processor_cpu_time.reset();

  // Perform required number of iterations (time steps).
  if (!IP->Time_Accurate && IP->Maximum_Number_of_Time_Steps > 0) {

    if (!batch_flag) cout << "\n Beginning FAS Multigrid computations on "
			  << Date_And_Time() << ".";

    // Determine if full multigrid is required and reset top level
    // appropriately.
    initial_top_level = FINEST_LEVEL;
    if (IP->Multigrid_IP.Number_of_Cycles_per_Stage_for_Full_Multigrid > 0) {
      if (!batch_flag) cout << "\n\n Perform Full multigrid cycles\n\n";
      initial_top_level = IP->Multigrid_IP.Levels-2;
    }

    // Perform the full or regular multigrid cycles as required.  For 
    // full multigrid the top level has been reset appropriately and 
    // saw-tooth cycles are performed at each coarse grid level until 
    // the finest grid level is reached.  At this point regular FAS 
    // multigrid cycles are performed with the cycle type specified in
    // the inpupt parameters.
    for (int top_level = initial_top_level; top_level >= FINEST_LEVEL; top_level--) {

      // Notify that the full multigrid cycles have been finished and
      // that the regular FAS multigrid cycles are being performed.
      if (!batch_flag && top_level == FINEST_LEVEL)
	cout << "\n\n Perform Regular multigrid cycles\n\n";

      // Unfreeze limiters if frozen.
      if (!first_step &&
	  IP->Freeze_Limiter &&
	  limiter_freezing == ON &&
	  IP->Limiter_Type != LIMITER_ZERO) {
	for (int level = top_level; level < IP->Multigrid_IP.Levels; level++) {
	  Evaluate_Limiters(Local_SolnBlk[level],
			    List_of_Local_Solution_Blocks[level]);
	}
	limiter_freezing = OFF;
      }

      // Perform the required number of multigrid cycles for this level.
      // If the current top grid level is not the finest grid level then
      // full multigrid cycles are being computed and the number of
      // cycles corresponds to the number specified in the input
      // parameters.  Otherwise, for regular FAS multigrid cycles, the
      // number of required cycles is set to an arbitrarily large number
      // since a separate exit criterion based on the number of time
      // steps is used to terminate the calculation.

      // Determine the number of required multigrid cycles.
      if (top_level != FINEST_LEVEL) {
	// The number of full multigrid cycles is specified in the input
	// parameters.
	cycles_for_this_level = IP->Multigrid_IP.Number_of_Cycles_per_Stage_for_Full_Multigrid;
      } else {
	// Set the number of cycles to an arbitrarily large number for
	// the regualar FAS multigrid calculations.
	cycles_for_this_level = max_cycles_for_this_level;
      }

      for (int cycles = 1; cycles <= cycles_for_this_level; cycles++) {

	// Determine the L1 norm of the solution residual.
	residual_l1_norm = EBSolver[top_level].L1_Norm_Residual();
	residual_l1_norm = CFFC_Summation_MPI(residual_l1_norm);	
	// Determine the L2 norm of the solution residual.
	residual_l2_norm = EBSolver[top_level].L2_Norm_Residual();
	residual_l2_norm = sqr(residual_l2_norm);
	residual_l2_norm = CFFC_Summation_MPI(residual_l2_norm);
	residual_l2_norm = sqrt(residual_l2_norm);
	// Determine the max norm of the solution residual.
	residual_max_norm = EBSolver[top_level].Max_Norm_Residual();
	residual_max_norm = CFFC_Maximum_MPI(residual_max_norm);

	if (cycles == 2) {
	   initial_residual_l2_norm = residual_l2_norm;
	} else if (cycles > 2 && fabs(initial_residual_l2_norm) > NANO) {
	   ratio_residual_l2_norm = residual_l2_norm / initial_residual_l2_norm;
	}

	// Update CPU time used for the calculation so far.
	processor_cpu_time.update();
	// Total CPU time for all processors. 
	total_cpu_time.cput = CFFC_Summation_MPI(processor_cpu_time.cput);

	// Periodically save restart solution files.
	if (!first_step && top_level == FINEST_LEVEL &&
	    number_of_time_steps-IP->Restart_Solution_Save_Frequency*
	    (number_of_time_steps/IP->Restart_Solution_Save_Frequency) == 0 ) {
	  if (!batch_flag) cout << "\n\n  Saving solution to restart data file(s) after"
				<< " n = " << number_of_time_steps << " steps (iterations).";
	  // Write the quadtree restart file.
	  error_flag = Write_QuadTree(*QuadTree,*IP);
	  if (error_flag) {
	    cout << "\n ERROR: Unable to open quadtree data file on processor "
		 << List_of_Local_Solution_Blocks[top_level].ThisCPU << "." << endl;
	  }
	  error_flag = CFFC_OR_MPI(error_flag);
	  if (error_flag) return error_flag;
	  // Write the solution block restart files.
	  error_flag = Write_Restart_Solution(Local_SolnBlk[top_level],
					      List_of_Local_Solution_Blocks[top_level],
					      *IP,
					      number_of_time_steps,
					      Time,
					      processor_cpu_time);
	  if (error_flag) {
	    cout << "\n ERROR: Unable to open restart output data file(s) on processor "
		 << List_of_Local_Solution_Blocks[top_level].ThisCPU << "." << endl;
	  }
	  error_flag = CFFC_OR_MPI(error_flag);
	  if (error_flag) return error_flag;
	  // Write the solution block restart files.
	  error_flag = EBSolver[top_level].Write_Restart_Files(number_of_time_steps,
							       levelset_iterations,
							       Time,
							       levelset_Time,
							       processor_cpu_time);
	  if (error_flag) {
	    cout << "\n ERROR: Unable to open restart output data file(s) "
		 << "on processor " << List_of_Local_Solution_Blocks[top_level].ThisCPU << "." << endl;
	  }
	  error_flag = CFFC_OR_MPI(error_flag);
	  if (error_flag) return error_flag;
	  if (!batch_flag) cout << endl;
	}

	// Output progress information for the calculation.
	if (!batch_flag) Output_Progress_L2norm(number_of_time_steps,
						Time*THOUSAND,
						total_cpu_time,
						residual_l2_norm,
						ratio_residual_l2_norm,
						first_step,
						IP->Output_Progress_Frequency);
	if (CFFC_Primary_MPI_Processor() && !first_step) {
	  Output_Progress_to_File(residual_file,
				  number_of_time_steps,
				  Time*THOUSAND,
				  total_cpu_time,
				  residual_l1_norm,
				  residual_l2_norm,
				  residual_max_norm);
	}

        // Periodically write out the nodal and cell solution information.
	if (IP->Multigrid_IP.Write_Output_Cells_Frequency > 0 &&
	    top_level == FINEST_LEVEL &&
	    (number_of_time_steps % IP->Multigrid_IP.Write_Output_Cells_Frequency == 0)) {
	   if (!batch_flag) { 
	      cout << endl << " Writing out solution in Tecplot format at iteration level ";
	      cout << number_of_time_steps << "." << endl; 
	   }
	   Output_Multigrid_Cells(number_of_time_steps, Time, true, residual_l2_norm, ratio_residual_l2_norm);
	   Output_Multigrid(number_of_time_steps, Time, true, residual_l2_norm, ratio_residual_l2_norm);
	} /* endif */

	// Check if the maximum number of time steps has been reached or
	// if the residual has dropped below the prescribed level for a
	// full multigrid cycle.
	bool please_stop = false;
	if (number_of_time_steps >= IP->Maximum_Number_of_Time_Steps) {
	  please_stop = true;
	}
	if (cycles > 2) {
	   if (top_level == FINEST_LEVEL) { // regular multigrid cycles
   	      if (residual_l2_norm < IP->Multigrid_IP.Absolute_Convergence_Tolerance) {
	         please_stop = true;
	  	 if (!batch_flag) {
 		    int tempp = cout.precision(); 
                    cout.precision(2); 
		    cout.setf(ios::scientific);
 		    cout << endl << "Regular Multigrid: met absolute convergence tolerance";
		    cout << " (" << residual_l2_norm << " < ";
		    cout << IP->Multigrid_IP.Absolute_Convergence_Tolerance << ")." << endl;
		    cout.precision(tempp);
		    cout.unsetf(ios::scientific);
		 } /* endif */
	      } /* endif */
  	      if (ratio_residual_l2_norm < IP->Multigrid_IP.Relative_Convergence_Tolerance) {
 		 please_stop = true;
		 if (!batch_flag) {
		    int tempp = cout.precision(); 
                    cout.precision(2);
		    cout.setf(ios::scientific);
	  	    cout << endl << "Regular Multigrid: met relative convergence tolerance";
		    cout << " (" << ratio_residual_l2_norm << " < ";
		    cout << IP->Multigrid_IP.Relative_Convergence_Tolerance << ")." << endl;
		    cout.precision(tempp);
		    cout.unsetf(ios::scientific);
		 } /* endif */
	      } /* endif */
	   } else { // otherwise doing full multigrid cycles.
	      if (residual_l2_norm < IP->Multigrid_IP.FMG_Absolute_Convergence_Tolerance) {
		 please_stop = true;
		 if (!batch_flag) {
		    int tempp = cout.precision(); 
                    cout.precision(2);
		    cout.setf(ios::scientific);
		    cout << endl << "Top Level: " << top_level;
		    cout << " Finest level: " << FINEST_LEVEL << ".";
		    cout << " Met absolute convergence tolerance";
		    cout << " (" << residual_l2_norm << " < ";
		    cout << IP->Multigrid_IP.FMG_Absolute_Convergence_Tolerance << ")." << endl;
		    cout.precision(tempp);
		    cout.unsetf(ios::scientific);
		 } /* endif */
	      } /* endif */
	      if (ratio_residual_l2_norm < IP->Multigrid_IP.FMG_Relative_Convergence_Tolerance) {
		 please_stop = true;
		 if (!batch_flag) {
		    int tempp = cout.precision(); 
                    cout.precision(2);
		    cout.setf(ios::scientific);
		    cout << endl << "Top Level: " << top_level;
		    cout << " Finest level: " << FINEST_LEVEL << ".";
		    cout << " Met relative convergence tolerance";
		    cout << " (" << ratio_residual_l2_norm << " < ";
		    cout << IP->Multigrid_IP.FMG_Relative_Convergence_Tolerance << ")." << endl;
		    cout.precision(tempp);
		    cout.unsetf(ios::scientific);
		 } /* endif */
	      } /* endif */
  	   } /* endif */
	} // if (cycles > 2) for convergence test

	if (please_stop) {
	  // Exit criteria has been met.  Output the final progress
	  // information for the calculation and exit.
	  if (!batch_flag) Output_Progress_L2norm(number_of_time_steps,
						  Time*THOUSAND,
						  total_cpu_time,
						  residual_l2_norm,
						  ratio_residual_l2_norm,
						  first_step,
						  number_of_time_steps);
	  break;
	}

	// Freeze the slope limiters for the inviscid flux calculation
	// if the residual is less than the value specified in the input
	// parameters.
	if (!first_step &&
	    cycles > 1 &&
	    IP->Freeze_Limiter &&
	    limiter_freezing == OFF &&
	    IP->Limiter_Type != LIMITER_ZERO &&
	    residual_l2_norm <= IP->Freeze_Limiter_Residual_Level) {
	  // Freeze the limiter on all multigrid grid levels.
	  for (int level = top_level; level < IP->Multigrid_IP.Levels; level++) {
	    Freeze_Limiters(Local_SolnBlk[level],
			    List_of_Local_Solution_Blocks[level]);
	  }
	  // Set the limiter freezing flag.
	  limiter_freezing = ON;
  	  if (CFFC_Primary_MPI_Processor()) {
	     cout << "Freezing Gradient Limiters." << endl;
	  }
	}

	// Update the solution for the next time-step/cylce using the 
	// FAS multigrid algorithm.
	error_flag = Coarse_Grid_Correction(top_level,
					    top_level,
					    ZERO,
					    Time);
	if (error_flag) return error_flag;

	// Set the first step flag.
	if (first_step) first_step = 0;

	// Increment the time step counter.
	number_of_time_steps++;

      }

      // If the current top grid level is not the finest grid level
      // (computing full multigrid cylces) then prolong solution up one
      // grid level.
      if (top_level != FINEST_LEVEL) {
	// Prolong the solution.
	error_flag = Prolong_Solution_Blocks(top_level);
	if (error_flag) return error_flag;
	// Update the ghostcell information to ensure that the solution
	// is consistent on each block.
	error_flag = Exchange_Solution_Information(top_level-1,
						   OFF);
	if (error_flag) return error_flag;
	// Apply the boundary conditions on the finer grid level.
	//if (IP->Multigrid_IP.Apply_Coarse_Mesh_Boundary_Conditions)
	EBSolver[top_level-1].Boundary_Conditions(Time);
      }

    }

    if (!batch_flag) cout << "\n\n FAS Multigrid computations complete on " 
			  << Date_And_Time() << ".\n";

  }

  // Update the ghostcell information to ensure that the solution is
  // consistent on each block.
  error_flag = Exchange_Solution_Information(FINEST_LEVEL,OFF);
  if (error_flag) return error_flag;

  // Apply the boundary conditions on finest mesh.
  EBSolver[FINEST_LEVEL].Boundary_Conditions(Time);

  // Close the residual file.
  if (CFFC_Primary_MPI_Processor()) error_flag = Close_Progress_File(residual_file);
  if (error_flag) return error_flag;

  // Solution calculations using multigrid complete.
  return 0;

}

/**********************************************************************
 **********************************************************************
 *****       Routines required for dual-time-stepping.       **********
 **********************************************************************
 **********************************************************************/

/**********************************************************************
 * EB_FAS_Multigrid2D_Solver::reallocate --                           *
 *                                                                    *
 * This routine performs the memory reallocation and reinitialization *
 * for the coarse grid levels after refinement.  Used by the dual-    *
 * time-stepping scheme.                                              *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState,
 	  class Quad_Soln_Block,
 	  class Quad_Soln_Input_Parameters>
int EB_FAS_Multigrid2D_Solver<cState, pState,
			      Quad_Soln_Block,
			      Quad_Soln_Input_Parameters>::
reallocate(void) {

  int error_flag;

  // Deallocate the embedded boundary solution blocks.
  for (int level = 1; level < IP->Multigrid_IP.Levels; level++) {    
    EBSolver[level].Local_SolnBlk = NULL;
    EBSolver[level].Local_Solution_Block_List = NULL;
  }

  // Deallocate the local solution blocks.
  for (int level = 1; level < IP->Multigrid_IP.Levels; level++) {    
    for (int nb = 0 ; nb < IP->Number_of_Blocks_Per_Processor; nb++) {
      if (Local_SolnBlk[level][nb].U != NULL) Local_SolnBlk[level][nb].deallocate();
    }
    if (Local_SolnBlk[level] != NULL) delete []Local_SolnBlk[level];
  }

  // Deallocate the FAS multigrid solution blocks.
  for (int level = 0; level < IP->Multigrid_IP.Levels; level++) {    
    if (MG_SolnBlk[level] != NULL) delete []MG_SolnBlk[level];
  }
  if (MG_SolnBlk != NULL) { delete []MG_SolnBlk; MG_SolnBlk = NULL; }

  // Deallocate the DTS multigrid solution blocks.
  for (int level = 0; level < IP->Multigrid_IP.Levels; level++) {    
    if (DTS_SolnBlk[level] != NULL) delete []DTS_SolnBlk[level];
  }
  if (DTS_SolnBlk != NULL) { delete []DTS_SolnBlk; DTS_SolnBlk = NULL; }

  // Deallocate the local lists of solution blocks.  
  for (int level = 1; level < IP->Multigrid_IP.Levels; level++) {
    if (List_of_Local_Solution_Blocks[level].Block != NULL) {
      List_of_Local_Solution_Blocks[level].deallocate();
    }
  }

  // Allocate memory for the coarse grid list of local solution blocks
  // and set the CPU number.
  for (int level = 1; level < IP->Multigrid_IP.Levels; level++) {
    List_of_Local_Solution_Blocks[level].allocate(IP->Number_of_Blocks_Per_Processor);
    List_of_Local_Solution_Blocks[level].ThisCPU = List_of_Local_Solution_Blocks[FINEST_LEVEL].ThisCPU;
  }

  // Allocate memory for the coarse grid local solution blocks.
  for (int level = 1; level < IP->Multigrid_IP.Levels; level++) {
    Local_SolnBlk[level] = new Quad_Soln_Block[IP->Number_of_Blocks_Per_Processor];
  }

  // Allocate memory for the FAS multigrid solution blocks on each level.
  MG_SolnBlk = new FAS_Multigrid_Quad_Block<cState>*[IP->Multigrid_IP.Levels];
  // Allocate memory for the coarse grid FAS multigrid solution blocks.
  for (int level = 0; level < IP->Multigrid_IP.Levels; level++) {
    MG_SolnBlk[level] = new FAS_Multigrid_Quad_Block<cState>[IP->Number_of_Blocks_Per_Processor];
  }

  // Allocate memory for the DTS multigrid solution blocks on each level.
  DTS_SolnBlk = new DTS_Multigrid_Quad_Block<cState>*[IP->Multigrid_IP.Levels];
  // Allocate memory for the coarse grid FAS multigrid solution blocks.
  for (int level = 0; level < IP->Multigrid_IP.Levels; level++) {
    DTS_SolnBlk[level] = new DTS_Multigrid_Quad_Block<cState>[IP->Number_of_Blocks_Per_Processor];
  }

  // Allocate memory and set data for all coarse mesh variables on all
  // blocks on this processor.
  for (int nb = 0; nb < IP->Number_of_Blocks_Per_Processor; nb++) {
    if (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].used == ADAPTIVEBLOCK2D_USED) {

      // Allocate memory for the forcing term and the uo storage on the
      // finest level.
      MG_SolnBlk[FINEST_LEVEL][nb].allocate(Local_SolnBlk[FINEST_LEVEL][nb].NCi,
					    Local_SolnBlk[FINEST_LEVEL][nb].NCj);

      // Allocate memory for the DTS solution block finest level.
      DTS_SolnBlk[FINEST_LEVEL][nb].allocate(Local_SolnBlk[FINEST_LEVEL][nb].NCi,
					     Local_SolnBlk[FINEST_LEVEL][nb].NCj);

      // Allocate memory and set data for the coarse mesh levels.
      for (int level = 1; level < IP->Multigrid_IP.Levels; level++) {

	// Copy the list of local solution block neighbour information
	// and calculate the coarse grid resolution.
	List_of_Local_Solution_Blocks[level].Block[nb] =
	  List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb];
	List_of_Local_Solution_Blocks[level].Block[nb].info.dimen.i = 
	  int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].info.dimen.i/pow(2.0,double(level)));
	List_of_Local_Solution_Blocks[level].Block[nb].info.dimen.j = 
	  int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].info.dimen.j/pow(2.0,double(level)));
	
	for (int n = 0; n < List_of_Local_Solution_Blocks[level].Block[nb].nS; n++) {
	  List_of_Local_Solution_Blocks[level].Block[nb].infoS[n].dimen.i = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoS[n].dimen.i/pow(2.0,double(level)));
	  List_of_Local_Solution_Blocks[level].Block[nb].infoS[n].dimen.j = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoS[n].dimen.j/pow(2.0,double(level)));
	}
	for (int n = 0; n < List_of_Local_Solution_Blocks[level].Block[nb].nN; n++) {
	  List_of_Local_Solution_Blocks[level].Block[nb].infoN[n].dimen.i = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoN[n].dimen.i/pow(2.0,double(level)));
	  List_of_Local_Solution_Blocks[level].Block[nb].infoN[n].dimen.j = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoN[n].dimen.j/pow(2.0,double(level)));
	}
	for (int n = 0; n < List_of_Local_Solution_Blocks[level].Block[nb].nE; n++) {
	  List_of_Local_Solution_Blocks[level].Block[nb].infoE[n].dimen.i = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoE[n].dimen.i/pow(2.0,double(level)));
	  List_of_Local_Solution_Blocks[level].Block[nb].infoE[n].dimen.j = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoE[n].dimen.j/pow(2.0,double(level)));
	}
	for (int n = 0; n < List_of_Local_Solution_Blocks[level].Block[nb].nW; n++) {
	  List_of_Local_Solution_Blocks[level].Block[nb].infoW[n].dimen.i = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoW[n].dimen.i/pow(2.0,double(level)));
	  List_of_Local_Solution_Blocks[level].Block[nb].infoW[n].dimen.j = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoW[n].dimen.j/pow(2.0,double(level)));
	}
	for (int n = 0; n < List_of_Local_Solution_Blocks[level].Block[nb].nSE; n++) {
	  List_of_Local_Solution_Blocks[level].Block[nb].infoSE[n].dimen.i = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoSE[n].dimen.i/pow(2.0,double(level)));
	  List_of_Local_Solution_Blocks[level].Block[nb].infoSE[n].dimen.j = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoSE[n].dimen.j/pow(2.0,double(level)));
	}
	for (int n = 0; n < List_of_Local_Solution_Blocks[level].Block[nb].nSW; n++) {
	  List_of_Local_Solution_Blocks[level].Block[nb].infoSW[n].dimen.i = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoSW[n].dimen.i/pow(2.0,double(level)));
	  List_of_Local_Solution_Blocks[level].Block[nb].infoSW[n].dimen.j = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoSW[n].dimen.j/pow(2.0,double(level)));
	}
	for (int n = 0; n < List_of_Local_Solution_Blocks[level].Block[nb].nNE; n++) {
	  List_of_Local_Solution_Blocks[level].Block[nb].infoNE[n].dimen.i = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoNE[n].dimen.i/pow(2.0,double(level)));
	  List_of_Local_Solution_Blocks[level].Block[nb].infoNE[n].dimen.j = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoNE[n].dimen.j/pow(2.0,double(level)));
	}
	for (int n = 0; n < List_of_Local_Solution_Blocks[level].Block[nb].nNW; n++) {
	  List_of_Local_Solution_Blocks[level].Block[nb].infoNW[n].dimen.i = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoNW[n].dimen.i/pow(2.0,double(level)));
	  List_of_Local_Solution_Blocks[level].Block[nb].infoNW[n].dimen.j = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoNW[n].dimen.j/pow(2.0,double(level)));
	}

	// Set-up the local block list for each level.
	List_of_Local_Solution_Blocks[level].Block[nb].used = ADAPTIVEBLOCK2D_USED;
	List_of_Local_Solution_Blocks[level].Block[nb].gblknum = 
	  List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].gblknum;
	List_of_Local_Solution_Blocks[level].Block[nb].info.cpu = 
	  List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].info.cpu;
	List_of_Local_Solution_Blocks[level].Block[nb].info.blknum = 
	  List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].info.blknum;
	List_of_Local_Solution_Blocks[level].Block[nb].info.dimen.i = 
	  int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].info.dimen.i/pow(2.0,double(level)));
	List_of_Local_Solution_Blocks[level].Block[nb].info.dimen.j = 
	  int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].info.dimen.j/pow(2.0,double(level)));

	List_of_Local_Solution_Blocks[level].Block[nb].info.dimen.ghost = 
	  List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].info.dimen.ghost;
	List_of_Local_Solution_Blocks[level].Block[nb].info.sector = 
	  List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].info.sector;
	List_of_Local_Solution_Blocks[level].Block[nb].info.level = 
	  List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].info.level;
	
	// Allocate the coarse grid local solution block.
	Local_SolnBlk[level][nb].allocate(List_of_Local_Solution_Blocks[level].Block[nb].info.dimen.i,
					  List_of_Local_Solution_Blocks[level].Block[nb].info.dimen.j,
					  List_of_Local_Solution_Blocks[level].Block[nb].info.dimen.ghost);

	// Create the coarse grid mesh.
	if (level == 1) Half_Mesh_Resolution(Local_SolnBlk[level][nb].Grid,
					     EBSolver[level-1].OGrid[nb]);
	else Half_Mesh_Resolution(Local_SolnBlk[level][nb].Grid,
				  Local_SolnBlk[level-1][nb].Grid);

	// Allocate the embedded boundary solution blocks.
	EBSolver[level].Local_SolnBlk = Local_SolnBlk[level];
	EBSolver[level].Local_Solution_Block_List = &List_of_Local_Solution_Blocks[level];

	// Allocate the coarse grid FAS multigrid solution block.
	MG_SolnBlk[level][nb].allocate(Local_SolnBlk[level][nb].NCi,
				       Local_SolnBlk[level][nb].NCj);

	// Allocate the coarse grid DTS multigrid solution block.
	DTS_SolnBlk[level][nb].allocate(Local_SolnBlk[level][nb].NCi,
					Local_SolnBlk[level][nb].NCj);

	// Allocate memory for the message passing buffers used to send solution
	// information between neighbouring blocks for the coarse grid.
	Allocate_Message_Buffers(List_of_Local_Solution_Blocks[level],
				 Local_SolnBlk[level][nb].NumVar()+NUM_COMP_VECTOR2D);

      }

    }
  }

  // Solution block reallocation and reassignment successful.
  return 0;

}

/**********************************************************************
 * EB_FAS_Multigrid2D_Solver::dUdtau_Multistage_Explicit --           *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState,
	  class Quad_Soln_Block,
	  class Quad_Soln_Input_Parameters>
int EB_FAS_Multigrid2D_Solver<cState, pState,
			      Quad_Soln_Block,
			      Quad_Soln_Input_Parameters>::
dUdtau_Multistage_Explicit(const int &Level,
			   const int &i_stage,
			   const double &dt) {

  // Exit immediately if dual-time-stepping is not required.
  if (IP->i_Time_Integration != TIME_STEPPING_DUAL_TIME_STEPPING) return 0;

  int error_flag, k_residual, time_integration;

  // Temporarily overwrite the time integration type to the multistage
  // optimally smoothing scheme.
  time_integration = IP->i_Time_Integration;
  IP->i_Time_Integration = IP->Multigrid_IP.i_Smoothing;

  switch(IP->i_Time_Integration) {
  case TIME_STEPPING_EXPLICIT_EULER :
    k_residual = 0;
    break;
  case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
    k_residual = 0;
    break;
  case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA :
    k_residual = 0;
    if (IP->N_Stage == 4) {
      if (i_stage == 4) k_residual = 0;
      else k_residual = i_stage - 1;
    }
    break;
  case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
    k_residual = 0;
    break;
  default:
    k_residual = 0;
    break;
  }

  // Evaluate the solution residual for each block.
  for (int nb = 0; nb < List_of_Local_Solution_Blocks[Level].Nblk; nb++) {
    if (List_of_Local_Solution_Blocks[Level].Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      for (int j = Local_SolnBlk[Level][nb].JCl; j <= Local_SolnBlk[Level][nb].JCu; j++) {
	for (int i = Local_SolnBlk[Level][nb].ICl; i <= Local_SolnBlk[Level][nb].ICu; i++) {
	  if (EBSolver[Level].Mesh[nb].cell_status[i][j] == CELL_STATUS_ACTIVE) {
	    if (IP->Multigrid_IP.i_Physical_Time_Integration == TIME_STEPPING_IMPLICIT_EULER) {
	      Local_SolnBlk[Level][nb].dUdt[i][j][k_residual] -= (IP->CFL_Number*Local_SolnBlk[Level][nb].dt[i][j])*(Local_SolnBlk[Level][nb].U[i][j] -
														     DTS_SolnBlk[Level][nb].Un[i][j])/
		                                                                                                    (IP->Multigrid_IP.Physical_Time_CFL_Number*dt);
	    } else if (IP->Multigrid_IP.i_Physical_Time_Integration == TIME_STEPPING_IMPLICIT_SECOND_ORDER_BACKWARD) {
	      Local_SolnBlk[Level][nb].dUdt[i][j][k_residual] -= (IP->CFL_Number*Local_SolnBlk[Level][nb].dt[i][j])*(3.0*Local_SolnBlk[Level][nb].U[i][j] -
														     4.0*DTS_SolnBlk[Level][nb].Un[i][j] +
														     DTS_SolnBlk[Level][nb].Uo[i][j])/
		                                                                                                     (TWO*IP->Multigrid_IP.Physical_Time_CFL_Number*dt);
	    }
	  }
	}
      }
    }
  }

  // Reset the time-integration type.
  IP->i_Time_Integration = time_integration;

  // Solution residual successfully determined.
  return 0;

}

/**********************************************************************
 * EB_FAS_Multigrid2D_Solver::Apply_Melson_Time_Step --               *
 **********************************************************************/
template <class cState, class pState,
	  class Quad_Soln_Block,
	  class Quad_Soln_Input_Parameters>
int EB_FAS_Multigrid2D_Solver<cState, pState,
			      Quad_Soln_Block,
			      Quad_Soln_Input_Parameters>::
Apply_Melson_Time_Step(const int &Level,
		       const int &i_stage,
		       const double &dt) {

  // Exit immediately if dual-time-stepping is not required.
  if (IP->i_Time_Integration != TIME_STEPPING_DUAL_TIME_STEPPING) return 0;

  int error_flag, k_residual, time_integration;

  // Temporarily overwrite the time integration type to the multistage
  // optimally smoothing scheme.
  time_integration = IP->i_Time_Integration;
  IP->i_Time_Integration = IP->Multigrid_IP.i_Smoothing;

  switch(IP->i_Time_Integration) {
  case TIME_STEPPING_EXPLICIT_EULER :
    k_residual = 0;
    break;
  case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
    k_residual = 0;
    break;
  case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA :
    k_residual = 0;
    if (IP->N_Stage == 4) {
      if (i_stage == 4) k_residual = 0;
      else k_residual = i_stage - 1;
    }
    break;
  case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
    k_residual = 0;
    break;
  default:
    k_residual = 0;
    break;
  }

  // Evaluate the solution residual for each block.
  for (int nb = 0; nb < List_of_Local_Solution_Blocks[Level].Nblk; nb++) {
    if (List_of_Local_Solution_Blocks[Level].Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      for (int j = Local_SolnBlk[Level][nb].JCl; j <= Local_SolnBlk[Level][nb].JCu; j++) {
	for (int i = Local_SolnBlk[Level][nb].ICl; i <= Local_SolnBlk[Level][nb].ICu; i++) {
	  if (IP->Multigrid_IP.i_Physical_Time_Integration == TIME_STEPPING_IMPLICIT_EULER) {
 	    Local_SolnBlk[Level][nb].dUdt[i][j][k_residual] /= ONE + ONE*IP->CFL_Number*Local_SolnBlk[Level][nb].dt[i][j]/(IP->Multigrid_IP.Physical_Time_CFL_Number*dt);
	  } else if (IP->Multigrid_IP.i_Physical_Time_Integration == TIME_STEPPING_IMPLICIT_SECOND_ORDER_BACKWARD) {
 	    Local_SolnBlk[Level][nb].dUdt[i][j][k_residual] /= ONE + 1.5*IP->CFL_Number*Local_SolnBlk[Level][nb].dt[i][j]/(IP->Multigrid_IP.Physical_Time_CFL_Number*dt);
	  }
	}
      }
    }
  }

  // Reset the time-integration type.
  IP->i_Time_Integration = time_integration;

  // Solution residual successfully determined.
  return 0;

}

/**********************************************************************
 * EB_FAS_Multigrid2D_Solver::Store_Previous_Solution --              *
 *                                                                    *
 * Store the previous solution and restrict to the coarse levels.     *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState,
	  class Quad_Soln_Block,
	  class Quad_Soln_Input_Parameters>
int EB_FAS_Multigrid2D_Solver<cState, pState,
			      Quad_Soln_Block,
			      Quad_Soln_Input_Parameters>::
Store_Previous_Solution(void) {

  int i_fine, j_fine;
  double A, total_area;
  Polygon Pc, Pf;

  for (int nb = 0; nb < List_of_Local_Solution_Blocks[FINEST_LEVEL].Nblk; nb++) {
    if (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      for (int j = Local_SolnBlk[FINEST_LEVEL][nb].JCl; j <= Local_SolnBlk[FINEST_LEVEL][nb].JCu; j++) {
	for (int i = Local_SolnBlk[FINEST_LEVEL][nb].ICl; i <= Local_SolnBlk[FINEST_LEVEL][nb].ICu; i++) {
	  if (EBSolver[FINEST_LEVEL].Mesh[nb].cell_status[i][j] == CELL_STATUS_ACTIVE) {
	    DTS_SolnBlk[FINEST_LEVEL][nb].Uo[i][j] = DTS_SolnBlk[FINEST_LEVEL][nb].Un[i][j];
	    DTS_SolnBlk[FINEST_LEVEL][nb].Un[i][j] = Local_SolnBlk[FINEST_LEVEL][nb].U[i][j];
	  }
	}
      }
    }
  }

  for (int level = 0; level < IP->Multigrid_IP.Levels-1; level++) {
    for (int nb = 0; nb < List_of_Local_Solution_Blocks[level+1].Nblk; nb++) {
      if (List_of_Local_Solution_Blocks[level+1].Block[nb].used == ADAPTIVEBLOCK2D_USED) {
 	for (int i_coarse = Local_SolnBlk[level+1][nb].ICl; i_coarse <= Local_SolnBlk[level+1][nb].ICu; i_coarse++) {
 	  for (int j_coarse = Local_SolnBlk[level+1][nb].JCl; j_coarse <= Local_SolnBlk[level+1][nb].JCu; j_coarse++) {
 	    // Determine the (i,j) index of the SW corner fine cell.
 	    i_fine = 2*(i_coarse-Local_SolnBlk[level][nb].Nghost)+Local_SolnBlk[level][nb].Nghost;
 	    j_fine = 2*(j_coarse-Local_SolnBlk[level][nb].Nghost)+Local_SolnBlk[level][nb].Nghost;

	    if (EBSolver[level+1].Mesh[nb].cell_status[i_coarse][j_coarse] != CELL_STATUS_ACTIVE) {
	      // The cell is inactive.  Set to standard atmosphere.
	      Local_SolnBlk[level+1][nb].U[i_coarse][j_coarse].Standard_Atmosphere();

	    } else if (!EBSolver[level+1].Interface_Union_List.Ni) {
	      // No embedded boundaries exist.  Restrict the solution
	      // information from the fine cells to the coarse cell using an
	      // area-weighted average of the fine cell solution information.
	      DTS_SolnBlk[level+1][nb].Un[i_coarse][j_coarse] = 
		(DTS_SolnBlk[level][nb].Un[i_fine  ][j_fine  ]*Local_SolnBlk[level][nb].Grid.Cell[i_fine  ][j_fine  ].A +
		 DTS_SolnBlk[level][nb].Un[i_fine+1][j_fine  ]*Local_SolnBlk[level][nb].Grid.Cell[i_fine+1][j_fine  ].A +
		 DTS_SolnBlk[level][nb].Un[i_fine  ][j_fine+1]*Local_SolnBlk[level][nb].Grid.Cell[i_fine  ][j_fine+1].A +
		 DTS_SolnBlk[level][nb].Un[i_fine+1][j_fine+1]*Local_SolnBlk[level][nb].Grid.Cell[i_fine+1][j_fine+1].A)/
		Local_SolnBlk[level+1][nb].Grid.Cell[i_coarse][j_coarse].A;
	      DTS_SolnBlk[level+1][nb].Uo[i_coarse][j_coarse] = 
		(DTS_SolnBlk[level][nb].Uo[i_fine  ][j_fine  ]*Local_SolnBlk[level][nb].Grid.Cell[i_fine  ][j_fine  ].A +
		 DTS_SolnBlk[level][nb].Uo[i_fine+1][j_fine  ]*Local_SolnBlk[level][nb].Grid.Cell[i_fine+1][j_fine  ].A +
		 DTS_SolnBlk[level][nb].Uo[i_fine  ][j_fine+1]*Local_SolnBlk[level][nb].Grid.Cell[i_fine  ][j_fine+1].A +
		 DTS_SolnBlk[level][nb].Uo[i_fine+1][j_fine+1]*Local_SolnBlk[level][nb].Grid.Cell[i_fine+1][j_fine+1].A)/
		Local_SolnBlk[level+1][nb].Grid.Cell[i_coarse][j_coarse].A;

	    } else {

	      if (EBSolver[level+1].Mesh[nb].cell_type[i_coarse][j_coarse] == CELL_TYPE_QUADRILATERAL &&
		  EBSolver[level].Mesh[nb].cell_type[i_fine  ][j_fine  ] == CELL_TYPE_QUADRILATERAL &&
		  EBSolver[level].Mesh[nb].cell_type[i_fine+1][j_fine  ] == CELL_TYPE_QUADRILATERAL &&
		  EBSolver[level].Mesh[nb].cell_type[i_fine  ][j_fine+1] == CELL_TYPE_QUADRILATERAL &&
		  EBSolver[level].Mesh[nb].cell_type[i_fine+1][j_fine+1] == CELL_TYPE_QUADRILATERAL &&
		  EBSolver[level].Mesh[nb].cell_status[i_fine  ][j_fine  ] == CELL_STATUS_ACTIVE &&
		  EBSolver[level].Mesh[nb].cell_status[i_fine+1][j_fine  ] == CELL_STATUS_ACTIVE &&
		  EBSolver[level].Mesh[nb].cell_status[i_fine  ][j_fine+1] == CELL_STATUS_ACTIVE &&
		  EBSolver[level].Mesh[nb].cell_status[i_fine+1][j_fine+1] == CELL_STATUS_ACTIVE) {
		// None of the associated cells have been adjusted.
		// Restrict the solution information using an area-weighted
		// average of the cells from the stored adjusted fine mesh
		// and the coarse mesh.
		DTS_SolnBlk[level+1][nb].Un[i_coarse][j_coarse] = 
		  (DTS_SolnBlk[level][nb].Un[i_fine  ][j_fine  ]*Local_SolnBlk[level][nb].Grid.Cell[i_fine  ][j_fine  ].A +
		   DTS_SolnBlk[level][nb].Un[i_fine+1][j_fine  ]*Local_SolnBlk[level][nb].Grid.Cell[i_fine+1][j_fine  ].A +
		   DTS_SolnBlk[level][nb].Un[i_fine  ][j_fine+1]*Local_SolnBlk[level][nb].Grid.Cell[i_fine  ][j_fine+1].A +
		   DTS_SolnBlk[level][nb].Un[i_fine+1][j_fine+1]*Local_SolnBlk[level][nb].Grid.Cell[i_fine+1][j_fine+1].A)/
		  Local_SolnBlk[level+1][nb].Grid.Cell[i_coarse][j_coarse].A;
		DTS_SolnBlk[level+1][nb].Uo[i_coarse][j_coarse] = 
		  (DTS_SolnBlk[level][nb].Uo[i_fine  ][j_fine  ]*Local_SolnBlk[level][nb].Grid.Cell[i_fine  ][j_fine  ].A +
		   DTS_SolnBlk[level][nb].Uo[i_fine+1][j_fine  ]*Local_SolnBlk[level][nb].Grid.Cell[i_fine+1][j_fine  ].A +
		   DTS_SolnBlk[level][nb].Uo[i_fine  ][j_fine+1]*Local_SolnBlk[level][nb].Grid.Cell[i_fine  ][j_fine+1].A +
		   DTS_SolnBlk[level][nb].Uo[i_fine+1][j_fine+1]*Local_SolnBlk[level][nb].Grid.Cell[i_fine+1][j_fine+1].A)/
		  Local_SolnBlk[level+1][nb].Grid.Cell[i_coarse][j_coarse].A;

	      } else {
		// The coarse and/or fine cells have been adjusted.
		// Restrict the solution information using an area-weighted
		// average based on the area of intersection of the coarse
		// grid and fine grid cells.
		total_area = ZERO;
		DTS_SolnBlk[level+1][nb].Un[i_coarse][j_coarse].Vacuum();
		DTS_SolnBlk[level+1][nb].Uo[i_coarse][j_coarse].Vacuum();
		// Create the coarse cell polygon.
		Pc.convert(Local_SolnBlk[level+1][nb].Grid.nodeSW(i_coarse,j_coarse).X,
			   Local_SolnBlk[level+1][nb].Grid.nodeSE(i_coarse,j_coarse).X,
			   Local_SolnBlk[level+1][nb].Grid.nodeNE(i_coarse,j_coarse).X,
			   Local_SolnBlk[level+1][nb].Grid.nodeNW(i_coarse,j_coarse).X);
		// Search the fine cells for intersecting cells.
		for (int j = (j_fine >= Local_SolnBlk[level][nb].JCl) ? (j_fine-2) : (j_fine);
		     (j_fine <= Local_SolnBlk[level][nb].JCu) ? (j < j_fine+4) : (j < j_fine+2); j++) {
		  for (int i = (i_fine >= Local_SolnBlk[level][nb].ICl) ? (i_fine-2) : (i_fine);
		       (i_fine <= Local_SolnBlk[level][nb].ICu) ? (i < i_fine+4) : (i < i_fine+2); i++) {
		    if (EBSolver[level].Mesh[nb].cell_status[i][j] == CELL_STATUS_ACTIVE) {
		      // Create the fine cell polygon.
		      Pf.convert(Local_SolnBlk[level][nb].Grid.nodeSW(i,j).X,
				 Local_SolnBlk[level][nb].Grid.nodeSE(i,j).X,
				 Local_SolnBlk[level][nb].Grid.nodeNE(i,j).X,
				 Local_SolnBlk[level][nb].Grid.nodeNW(i,j).X);
		      // Determine the area of intersection between the fine
		      // and coarse cell polygons.
		      A = Polygon_Intersection_Area(Pc,Pf);
		      Pf.deallocate();
		      // Add the area of intersection to the total area of
		      // intersection.
		      total_area += A;
		      // Restrict fine cell solution to the coarse cell
		      // based on the area of interesection.
		      DTS_SolnBlk[level+1][nb].Un[i_coarse][j_coarse] += DTS_SolnBlk[level][nb].Un[i][j]*A;
		      DTS_SolnBlk[level+1][nb].Uo[i_coarse][j_coarse] += DTS_SolnBlk[level][nb].Uo[i][j]*A;
		    }
		  }
		}
		Pc.deallocate();
		// Complete solution restriction by dividing the restricted
		// solution by the total area of intersection.
		DTS_SolnBlk[level+1][nb].Un[i_coarse][j_coarse] /= total_area;
		DTS_SolnBlk[level+1][nb].Uo[i_coarse][j_coarse] /= total_area;

	      }

	    }

 	  }
 	}
      }
    }
  }

  // Storage of previous solution conducted successfully.
  return 0;

}

/**********************************************************************
 * EB_FAS_Multigrid2D_Solver::Compute_Interface_Location --           *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState,
 	  class Quad_Soln_Block,
 	  class Quad_Soln_Input_Parameters>
int EB_FAS_Multigrid2D_Solver<cState, pState,
			      Quad_Soln_Block,
			      Quad_Soln_Input_Parameters>::
Compute_Interface_Location(const int &batch_flag,
			   const double &current_time,
			   const double &maximum_time,
			   int &evolution_counter,
			   int &levelset_iterations,
			   double &level_set_current_time) {

  // Exit immediately if no interface components have been specified.
  if (!EBSolver[FINEST_LEVEL].Interface_Component_List.Ni) return 0;

  int error_flag, motion_flag;

  // Initialize motion flag.
  motion_flag = OFF;

  // Determine the velocity function for each embedded boundary.
  for (int ni = 1; ni <= EBSolver[FINEST_LEVEL].Interface_Component_List.Ni; ni++) {
    if (EBSolver[FINEST_LEVEL].Interface_Component_List[ni].Motion == INTERFACE_MOTION_CONSTANT ||
	EBSolver[FINEST_LEVEL].Interface_Component_List[ni].Motion == INTERFACE_MOTION_EXPAND ||
	EBSolver[FINEST_LEVEL].Interface_Component_List[ni].Motion == INTERFACE_MOTION_UNIFORM ||
	EBSolver[FINEST_LEVEL].Interface_Component_List[ni].Motion == INTERFACE_MOTION_TRANSLATE ||
	EBSolver[FINEST_LEVEL].Interface_Component_List[ni].Motion == INTERFACE_MOTION_ROTATE ||
	EBSolver[FINEST_LEVEL].Interface_Component_List[ni].Motion == INTERFACE_MOTION_MOMENTUM_TRANSFER ||
	EBSolver[FINEST_LEVEL].Interface_Component_List[ni].Motion == INTERFACE_MOTION_LEVELSET_EXPAND ||
	EBSolver[FINEST_LEVEL].Interface_Component_List[ni].Motion == INTERFACE_MOTION_LEVELSET_STRETCH ||
	EBSolver[FINEST_LEVEL].Interface_Component_List[ni].Motion == INTERFACE_MOTION_LEVELSET) {
      motion_flag = ON;
      break;
    }
  }

  // If none of the interface components have undergo any form of 
  // motion then exit immediately.
  if (!motion_flag) return 0;

  // Evolve the embedded interface(s) on the finest grid.
  error_flag = EBSolver[FINEST_LEVEL].Compute_Interface_Location(batch_flag,
								 current_time,
								 maximum_time,
								 evolution_counter,
								 levelset_iterations,
								 level_set_current_time);
  if (error_flag) return error_flag;

  error_flag = Exchange_Solution_Information(FINEST_LEVEL,OFF);
  if (error_flag) return error_flag;

  // Adjust the coarse grids according to the interface location(s).
  for (int level = 1; level < IP->Multigrid_IP.Levels; level++) {
    EBSolver[level].Interface_Component_List.Copy(EBSolver[FINEST_LEVEL].Interface_Component_List);
    EBSolver[level].Interface_Union_List.Copy(EBSolver[FINEST_LEVEL].Interface_Union_List);
    EBSolver[level].Store_Adjusted_Mesh();
    EBSolver[level].Mesh_Unadjustment();
    error_flag = EBSolver[level].Mesh_Adjustment(OFF,OFF);
    if (!error_flag) EBSolver[level].Store_Adjusted_Mesh();
    if (!error_flag) Restrict_Solution_Blocks(level-1);
    if (!error_flag) error_flag = Exchange_Solution_Information(level,OFF);
    if (!error_flag) error_flag = EBSolver[level].Boundary_Conditions(maximum_time);
    if (error_flag) return error_flag;
  }

  // The new interface location has been successfully computed.
  return 0;

}

/**********************************************************************
 * EB_FAS_Multigrid2D_Solver::DTS_Multigrid_Solution --               *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState,
 	  class Quad_Soln_Block,
 	  class Quad_Soln_Input_Parameters>
int EB_FAS_Multigrid2D_Solver<cState, pState,
			      Quad_Soln_Block,
			      Quad_Soln_Input_Parameters>::
DTS_Multigrid_Solution(int &batch_flag,
		       int &number_of_time_steps, 
		       int &evolution_counter,
		       int &levelset_iterations,
		       double &Time,
		       double &levelset_Time,
		       CPUTime &processor_cpu_time,
		       CPUTime &total_cpu_time,
		       ofstream &residual_file) {

  int error_flag, command_flag, progress_character = 1, limiter_freezing = OFF;
  int physical_first_step = ON, first_step = ON;
  int cycles = 0, initial_top_level = 0;
  int physical_time = IP->Multigrid_IP.i_Physical_Time_Integration;
  double dTime, residual_l1_norm, residual_l2_norm, residual_max_norm;

  /////////////////////////////////////////////////////////////
  // Apply the initial conditions on all coarse grid levels. //
  /////////////////////////////////////////////////////////////

  for (int level = 1; level < IP->Multigrid_IP.Levels; level++) {
    error_flag = Exchange_Solution_Information(level,ON);
    CFFC_Broadcast_MPI(&error_flag,1);
    if (error_flag) return error_flag;
    EBSolver[level].Initialize_Adjustment_Grids();
    error_flag = EBSolver[level].Mesh_Adjustment(ON,OFF);
    if (error_flag) return error_flag;
    error_flag = EBSolver[level].Compute_Interface_Velocity_Function(Time);
    if (error_flag) return error_flag;
    ICs(Local_SolnBlk[level],
	List_of_Local_Solution_Blocks[level],
	*IP);
    if (IP->i_ICs == IC_RESTART) {
      Restrict_Solution_Blocks(level-1);
      Update_Primitive_Variables(level);
    }
    // Restrict the boundary reference states if required.
    if (IP->Multigrid_IP.Apply_Coarse_Mesh_Boundary_Conditions)
      Restrict_Boundary_Ref_States(level-1);
    // Send solution information between neighbouring blocks to
    // complete the prescription of initial data.
    error_flag = Exchange_Solution_Information(level,OFF);
    CFFC_Broadcast_MPI(&error_flag,1);
    if (error_flag) return error_flag;
  }

  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

  ////////////////////////////////////////////////////
  // Open the progress file and reset the CPU time. //
  ////////////////////////////////////////////////////

  if (CFFC_Primary_MPI_Processor()) {
    error_flag = Open_Progress_File(residual_file,
				    IP->Output_File_Name,
				    number_of_time_steps);
    if (error_flag) {
      cout << endl << " ERROR: Unable to open residual file for calculation." << endl;
    }
  }
  CFFC_Broadcast_MPI(&error_flag,1);
  if (error_flag) return error_flag;

  // Set the CPU time to zero.
  if (IP->i_ICs != IC_RESTART) processor_cpu_time.zero();

  /////////////////////////////////////////////////////////
  // Perform required number of iterations (time steps). //
  /////////////////////////////////////////////////////////

  if ((!IP->Time_Accurate && IP->Maximum_Number_of_Time_Steps > 0) ||
      (IP->Time_Accurate && IP->Time_Max > ZERO)) {

    if (!batch_flag && CFFC_Primary_MPI_Processor()) cout << endl << " Beginning FAS Multigrid DTS computations on "
							    << Date_And_Time() << "." << endl << endl;

    // Perform required number of iterations (time steps).
    while ((!IP->Time_Accurate && number_of_time_steps < IP->Maximum_Number_of_Time_Steps) ||
	   (IP->Time_Accurate && IP->Time_Max > Time)) {

      // Periodically refine the mesh if required.
      if (IP->AMR) {
	if (!physical_first_step &&
	    number_of_time_steps-IP->AMR_Frequency*(number_of_time_steps/IP->AMR_Frequency) == 0) {
	  if (!batch_flag) cout << "\n\n Refining Grid.  Performing adaptive mesh refinement at n = "
				<< number_of_time_steps << ".";
	  Evaluate_Limiters(Local_SolnBlk[FINEST_LEVEL], 
			    List_of_Local_Solution_Blocks[FINEST_LEVEL]);
	  error_flag = EBSolver[FINEST_LEVEL].Adaptive_Mesh_Refinement(ON,ON);
	  if (error_flag) {
	    cout << "\n ERROR: FASMultigrid DTS AMR error on processor "
		 << List_of_Local_Solution_Blocks[FINEST_LEVEL].ThisCPU
		 << "." << endl;
	  }
	  error_flag = CFFC_OR_MPI(error_flag);
	  if (error_flag) {
	    command_flag = Output_Tecplot(Local_SolnBlk[FINEST_LEVEL],
					  List_of_Local_Solution_Blocks[FINEST_LEVEL],
					  *IP,
					  number_of_time_steps,
					  Time);
	    return error_flag;
	  }
	  if (!batch_flag) {
	    cout << "\n New multi-block solution-adaptive quadrilateral mesh statistics: "; 
	    cout << "\n  -> Number of Root Blocks i-direction: " << QuadTree->NRi;
	    cout << "\n  -> Number of Root Blocks j-direction: " << QuadTree->NRj;
	    cout << "\n  -> Total Number of Used Blocks: " << QuadTree->countUsedBlocks();
	    cout << "\n  -> Total Number of Computational Cells: " << QuadTree->countUsedCells();
	    cout << "\n  -> Number of Mesh Refinement Levels: " << QuadTree->highestRefinementLevel()+1;
	    cout << "\n  -> Refinement Efficiency: " << QuadTree->efficiencyRefinement();
	    cout << endl;
	  }
	  // Reallocate coarse grid local solution blocks and the FAS and
	  // DTS solution blocks after each mesh refinement.
	  reallocate();
  	  for (int level = 1; level < IP->Multigrid_IP.Levels; level++) {
 	    error_flag = Exchange_Solution_Information(level,ON);
	    if (error_flag) break;
  	  }
	  if (error_flag) return error_flag;
 	  for (int level = 1; level < IP->Multigrid_IP.Levels; level++) {
	    ICs(Local_SolnBlk[level],
		List_of_Local_Solution_Blocks[level],
		*IP);
  	    Restrict_Solution_Blocks(level-1);
 	    Update_Primitive_Variables(level);
 	    if (IP->Multigrid_IP.Apply_Coarse_Mesh_Boundary_Conditions) {
 	      Restrict_Boundary_Ref_States(level-1);
 	      BCs(Local_SolnBlk[level],
 		  List_of_Local_Solution_Blocks[level],
 		  *IP);
 	    }
  	  }
  	  for (int level = 1; level < IP->Multigrid_IP.Levels; level++) {
 	    error_flag = Exchange_Solution_Information(level,OFF);
 	    if (error_flag) break;
 	  }
	  if (error_flag) return error_flag;
	  // Reset the first step flag.
	  physical_first_step = ON;
	}
      }

      // Periodically save restart solution files.
      if (!physical_first_step &&
	  number_of_time_steps-IP->Restart_Solution_Save_Frequency*
	  (number_of_time_steps/IP->Restart_Solution_Save_Frequency) == 0) {
	if (!batch_flag) cout << "\n\n  Saving DTS FASMultigrid solution to restart data file(s) after"
			      << " n = " << number_of_time_steps << " steps (iterations).";
	error_flag = Write_QuadTree(*QuadTree,*IP);
	if (error_flag) {
	  cout << "\n ERROR: Unable to open quadtree data file "
	       << "on processor "
	       << List_of_Local_Solution_Blocks[FINEST_LEVEL].ThisCPU
	       << "." << endl;
	}
	error_flag = CFFC_OR_MPI(error_flag);
	if (error_flag) return error_flag;
	error_flag = Write_Restart_Solution(Local_SolnBlk[FINEST_LEVEL], 
					    List_of_Local_Solution_Blocks[FINEST_LEVEL], 
					    *IP,
					    number_of_time_steps,
					    Time,
					    processor_cpu_time);
	if (error_flag) {
	  cout << "\n ERROR: Unable to open restart output data file(s) "
	       << "on processor "
	       << List_of_Local_Solution_Blocks[FINEST_LEVEL].ThisCPU
	       << "." << endl;
	}
	error_flag = CFFC_OR_MPI(error_flag);
	if (error_flag) return error_flag;
	error_flag = EBSolver[FINEST_LEVEL].Write_Restart_Files(number_of_time_steps,
								levelset_iterations,
								Time,
								levelset_Time,
								processor_cpu_time);
	if (error_flag) {
	  cout << "\n ERROR: Unable to open restart output data file(s) "
	       << "on processor " << List_of_Local_Solution_Blocks[FINEST_LEVEL].ThisCPU << "." << endl;
	}
	error_flag = CFFC_OR_MPI(error_flag);
	if (error_flag) return error_flag;
	if (!batch_flag) cout << endl;
      }

      if (physical_first_step) {
	physical_time = IP->Multigrid_IP.i_Physical_Time_Integration;
	IP->Multigrid_IP.i_Physical_Time_Integration = TIME_STEPPING_IMPLICIT_EULER;
      }

      error_flag = Store_Previous_Solution();
      if (error_flag) return error_flag;

      // Determine local and global time steps.
      dTime = EBSolver[FINEST_LEVEL].CFL(Time);
      dTime = CFFC_Minimum_MPI(dTime); // Find global minimum time step for all processors.
      if (IP->Time_Accurate) {
	if (Time + IP->Multigrid_IP.Physical_Time_CFL_Number*dTime > IP->Time_Max) {
	  //cout << endl << dTime << " " << Time + TWO*IP->CFL_Number*dTime;
// 	  dTime = (IP->Time_Max-Time)/IP->CFL_Number;
	}
      }

      // Determine if full multigrid is required and reset top level
      // appropriately.
      if (IP->Multigrid_IP.Ncycles_Full_Multigrid > 0) {
	initial_top_level = IP->Multigrid_IP.Levels-2;
	progress_character = 1;
      }

      // Perform the full or regular multigrid cycles as required.  For 
      // full multigrid the top level has been reset appropriately and 
      // saw-tooth cycles are performed at each coarse grid level until 
      // the finest grid level is reached.  At this point regular DTS 
      // multigrid cycles are performed with the cycle type specified in
      // the inpupt parameters.
      for (int top_level = initial_top_level; top_level >= FINEST_LEVEL; top_level--) {

	// Notify that the full multigrid cycles have been finished and
	// that the regular FAS multigrid cycles are being performed.
	if (top_level == FINEST_LEVEL) progress_character = 0;

	// Unfreeze limiters if frozen.
	if (!first_step &&
	    IP->Freeze_Limiter &&
	    limiter_freezing == ON &&
	    IP->Limiter_Type != LIMITER_ZERO) {
	  for (int level = top_level; level < IP->Multigrid_IP.Levels; level++) {
	    Evaluate_Limiters(Local_SolnBlk[level],
			      List_of_Local_Solution_Blocks[level]);
	  }
	  limiter_freezing = OFF;
	}

	// Perform the required number of multigrid cycles for this level.
	// If the current top grid level is not the finest grid level then
	// full multigrid cycles are being computed and the number of
	// cycles corresponds to the number specified in the input
	// parameters.  Otherwise, for regular DTS multigrid cycles, the
	// number of required cycles is set to an arbitrarily large number
	// since a separate exit criterion based on the number of time
	// steps is used to terminate the calculation.

	cycles = 0;

	while ((top_level != 0 && cycles < IP->Multigrid_IP.Ncycles_Full_Multigrid) ||
	       (top_level == 0 && cycles < IP->Multigrid_IP.Ncycles_Regular_Multigrid)) {

	  // Determine the L1 norm of the solution residual.
	  residual_l1_norm = EBSolver[top_level].L1_Norm_Residual();
	  residual_l1_norm = CFFC_Summation_MPI(residual_l1_norm);	
	  // Determine the L2 norm of the solution residual.
	  residual_l2_norm = EBSolver[top_level].L2_Norm_Residual();
	  residual_l2_norm = sqr(residual_l2_norm);
	  residual_l2_norm = CFFC_Summation_MPI(residual_l2_norm);
	  residual_l2_norm = sqrt(residual_l2_norm);
	  // Determine the max norm of the solution residual.
	  residual_max_norm = EBSolver[top_level].Max_Norm_Residual();
	  residual_max_norm = CFFC_Maximum_MPI(residual_max_norm);

	  // Update CPU time used for the calculation so far.
	  processor_cpu_time.update();
	  total_cpu_time.cput = CFFC_Summation_MPI(processor_cpu_time.cput);

	  // Output progress information for the calculation.
	  if (!batch_flag) Output_Progress_L2norm(number_of_time_steps*IP->Multigrid_IP.Ncycles_Regular_Multigrid+cycles,
						  Time*THOUSAND,
						  total_cpu_time,
						  residual_l2_norm,
						  first_step,
						  IP->Output_Progress_Frequency,
						  progress_character);
	  if (CFFC_Primary_MPI_Processor() && !first_step) {
	    Output_Progress_to_File(residual_file,
				    number_of_time_steps*IP->Multigrid_IP.Ncycles_Regular_Multigrid+cycles,
				    Time*THOUSAND,
				    total_cpu_time,
				    residual_l1_norm,
				    residual_l2_norm,
				    residual_max_norm);
	  }

	  // Freeze the slope limiters for the inviscid flux calculation
	  // if the residual is less than the value specified in the input
	  // parameters.
	  if (!first_step &&
	      cycles > 1 &&
	      IP->Freeze_Limiter &&
	      limiter_freezing == OFF &&
	      IP->Limiter_Type != LIMITER_ZERO &&
	      residual_l2_norm <= IP->Freeze_Limiter_Residual_Level) {
	    // Freeze the limiter on all multigrid grid levels.
	    for (int level = top_level; level < IP->Multigrid_IP.Levels; level++) {
	      Freeze_Limiters(Local_SolnBlk[level],
			      List_of_Local_Solution_Blocks[level]);
	    }
	    // Set the limiter freezing flag.
	    limiter_freezing = ON;
	  }

	  // Update the solution for the next time-step/cylce using the 
	  // DTS multigrid algorithm.
	  error_flag = Coarse_Grid_Correction(top_level,top_level,dTime,Time);
	  if (error_flag) return error_flag;

	  // Set the first step flag.
	  if (first_step) first_step = 0;

	  // Increment the cycle counter.
	  cycles++;

	}

	// If the current top grid level is not the finest grid level
	// (computing full multigrid cylces) then prolong solution up one
	// grid level.
	if (top_level != 0) {
	  error_flag = Prolong_Solution_Blocks(top_level);
	  if (error_flag) return error_flag;
	  // Update the ghostcell information to ensure that the solution
	  // is consistent on each block.
	  error_flag = Exchange_Solution_Information(top_level-1,
						     OFF);
	  if (error_flag) return error_flag;
	  // Apply the boundary conditions on the finer grid level.
	  if (IP->Multigrid_IP.Apply_Coarse_Mesh_Boundary_Conditions) {
	    if (error_flag) return error_flag;
	    EBSolver[top_level-1].Boundary_Conditions(Time);
	  }
	}

      }

      // Evolve the embedded interface(s) and readjust grids.
      error_flag = Compute_Interface_Location(batch_flag,
					      Time,
					      Time+IP->CFL_Number*dTime,
					      evolution_counter,
					      levelset_iterations,
					      levelset_Time);
      if (error_flag) {
 	cout << "\n ERROR: Compute interface location error on processor "
 	     << List_of_Local_Solution_Blocks[FINEST_LEVEL].ThisCPU 
 	     << ".  Error number = " << error_flag << "." << endl;
      }
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;

      // Increment the time-step counter.
      number_of_time_steps++;
      if (physical_first_step) {
	physical_first_step = OFF;
	IP->Multigrid_IP.i_Physical_Time_Integration = physical_time;
      }

      // Update the time.
      Time += IP->Multigrid_IP.Physical_Time_CFL_Number*dTime;

    }

    if (!batch_flag && CFFC_Primary_MPI_Processor()) cout << "\n\n FAS Multigrid DTS computations complete on " 
							    << Date_And_Time() << ".\n";

  }

  // Update the ghostcell information to ensure that the solution is
  // consistent on each block.
  error_flag = Exchange_Solution_Information(FINEST_LEVEL,OFF);
  if (error_flag) return error_flag;

  // Apply the boundary conditions on finest mesh.
  EBSolver[FINEST_LEVEL].Boundary_Conditions(Time);

  ///////////////////////////////////////////////////
  // DTS multigrid solution computed successfully. //
  ///////////////////////////////////////////////////

  // Close the residual file.
  if (CFFC_Primary_MPI_Processor()) {
    error_flag = Close_Progress_File(residual_file);
    if (error_flag) return error_flag;
  }

  // DTS multigrid solution computed successfully.
  return 0;

}

#endif // _FASMULTIGRID_EMBEDDEDBOUNDARIES_INCLUDED

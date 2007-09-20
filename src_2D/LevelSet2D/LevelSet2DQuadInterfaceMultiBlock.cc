/**********************************************************************
 * LevelSet2DQuadInterfaceMultiBlock.cc:                              *
 *              Multi-block versions of subroutines for 2D Level Set  *
 *              multi-block quadrilateral mesh solution classes.      *
 **********************************************************************/

// Include 2D LevelSet quadrilateral mesh solution header file.

#ifndef _LEVELSET2D_QUAD_INCLUDED
#include "LevelSet2DQuad.h"
#endif // _LEVELSET2D_QUAD_INCLUDED

/**********************************************************************
 * LevelSet2D_Quad_Block -- Interface Multiple Block External         *
 *                          Subroutines.                              *
 **********************************************************************/

/**********************************************************************
 * Routine: Set_Interface_List                                        *
 *                                                                    *
 * Set the interface list to the given list on each of the solution   *
 * blocks in the 1D array of 2D quadrilateral multi-block solution    *
 * blocks.                                                            *
 *                                                                    *
 **********************************************************************/
int Set_Interface_List(LevelSet2D_Quad_Block *Soln_ptr,
		       AdaptiveBlock2D_List &Soln_Block_List,
		       const Interface2D_List &Interface_List) {

  int error_flag;

  // Set the interface list for each solution block.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      error_flag = Set_Interface_List(Soln_ptr[nb],
				      Interface_List);
      if (error_flag) return error_flag;
    }
  }

  // Return the error flag.
  return 0;

}

/**********************************************************************
 * Routine: Initialize_Interfaces                                     *
 *                                                                    *
 * Initializes the interface(s) on each of the solution blocks in the *
 * 1D array of 2D quadrilateral multi-block solution blocks.          *
 *                                                                    *
 **********************************************************************/
int Initialize_Interfaces(LevelSet2D_Quad_Block *Soln_ptr,
			  AdaptiveBlock2D_List &Soln_Block_List,
			  LevelSet2D_Input_Parameters &Input_Parameters) {

  int error_flag;

  // Assign initial data for each solution block.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      error_flag = Initialize_Interfaces(Soln_ptr[nb],
					 Input_Parameters);
      if (error_flag) return error_flag;
    }
  }

  // Return the error flag.
  return 0;

}

int Initialize_Interfaces(LevelSet2D_Quad_Block *Soln_ptr,
			  AdaptiveBlock2D_List &Soln_Block_List,
			  LevelSet2D_Input_Parameters &Input_Parameters,
			  const Interface2D_List &Interface_List) {

  int error_flag;

  // Assign initial data for each solution block.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      error_flag = Initialize_Interfaces(Soln_ptr[nb],
					 Input_Parameters,
					 Interface_List);
      if (error_flag) return error_flag;
    }
  }

  // Return the error flag.
  return 0;

}

/**********************************************************************
 * Routine: Exact_Initial_Extension                                   *
 *                                                                    *
 * Exact initialization of the level set function and the extended    *
 * front speeds.  Only available for line, circle, ellipse, and       *
 * square interface.                                                  *
 *                                                                    *
 **********************************************************************/
int Exact_Initial_Extension(LevelSet2D_Quad_Block *Soln_ptr,
			    AdaptiveBlock2D_List &Soln_Block_List,
			    LevelSet2D_Input_Parameters &Input_Parameters) {

  int error_flag;

  // Compute the signed distances and front speed exactly.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      error_flag = Exact_Initial_Extension(Soln_ptr[nb],
					   Input_Parameters);
      if (error_flag) return error_flag;
    }
  }

  error_flag = Scalar_Geometric_Extension_Problem(Soln_ptr,
						  Soln_Block_List,
						  Input_Parameters);
  if (error_flag) return error_flag;

  // Geometric extension problem computed successfully.
  return 0;

}

/**********************************************************************
 * Routine: Geometric_Extension_Problem                               *
 *                                                                    *
 * Geometric solution of the extension problems involving the level   *
 * set function and the normal front speed on each of the solution    *
 * blocks in the 1D array of 2D quadrilateral multi-block solution    *
 * blocks.                                                            *
 *                                                                    *
 **********************************************************************/
int Geometric_Extension_Problem(LevelSet2D_Quad_Block *Soln_ptr,
				AdaptiveBlock2D_List &Soln_Block_List,
				LevelSet2D_Input_Parameters &Input_Parameters) {

  int error_flag;

  // Compute the signed distances and front speed local to the interface.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      error_flag = Geometric_Extension_Problem(Soln_ptr[nb],
					       Input_Parameters);
      if (error_flag) return error_flag;
    }
  }

  // Geometric extension problem computed successfully.
  return 0;

}

/**********************************************************************
 * Routine: Scalar_Geometric_Extension_Problem                        *
 *                                                                    *
 * Geometric solution of the extension problems involving the level   *
 * set function and the normal front speed on each of the solution    *
 * blocks in the 1D array of 2D quadrilateral multi-block solution    *
 * blocks.                                                            *
 *                                                                    *
 **********************************************************************/
int Scalar_Geometric_Extension_Problem(LevelSet2D_Quad_Block *Soln_ptr,
				       AdaptiveBlock2D_List &Soln_Block_List,
				       LevelSet2D_Input_Parameters &Input_Parameters) {

  int error_flag;

  // Compute the signed distances and front speed local to the interface.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      error_flag = Scalar_Geometric_Extension_Problem(Soln_ptr[nb],
						      Input_Parameters);
      if (error_flag) return error_flag;
    }
  }

  // Scalar geometric extension problem computed successfully.
  return 0;

}

/**********************************************************************
 * Routine: RetrieveInterfaceSpline                                   *
 *                                                                    *
 * This routine locates the zero level set contained within the 1D    *
 * array of 2D quadrilateral multi-block solution blocks.  The        *
 * interface location is saved as a spline(s).                        *
 *                                                                    *
 **********************************************************************/
int Retrieve_Interface_Spline(LevelSet2D_Quad_Block *Soln_ptr,
			      AdaptiveBlock2D_List &Soln_Block_List) {

  int error_flag;

#ifdef _RETRIEVE_DEBUG_
  char extension[256], output_file_name[256];
  char *output_file_name_ptr;
  ofstream dout;
  strcpy(output_file_name,"retrieve_cpu");
  sprintf(extension,"%.6d",Soln_Block_List.ThisCPU);
  strcat(extension,".txt");
  strcat(output_file_name,extension);
  output_file_name_ptr = output_file_name;
  dout.open(output_file_name_ptr,ios::out);
  if (dout.bad()) return ;
#endif

  // Retrieve spline data for each solution block.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
#ifdef _RETRIEVE_DEBUG_
      dout << endl << " gblknum = " << Soln_Block_List.Block[nb].gblknum; dout.flush();
      error_flag = Retrieve_Interface_Spline(Soln_ptr[nb],
 					     Soln_Block_List.Block[nb].gblknum,
 					     dout);
#endif
#ifndef _RETRIEVE_DEBUG_
      error_flag = Retrieve_Interface_Spline(Soln_ptr[nb],
    					     Soln_Block_List.Block[nb].gblknum);
#endif
      if (error_flag) return error_flag;
    }
  }

#ifdef _RETRIEVE_DEBUG_
  dout << endl << " -> End Retrieve_Interface_Spline"; dout.flush();
  dout.close();
#endif

  // Interface capture was successful.
  return 0;

}

/**********************************************************************
 * Routine: Share_Interface_Information                               *
 *                                                                    *
 * The information regarding each interface must be known by all of   *
 * the solution blocks.  This function collects the interface         *
 * information on each processor and shares it with the other         *
 * processors if running a parallel case.  Interfaces that are        *
 * actually just segments of a bigger interface are concatenated.     *
 *                                                                    *
 **********************************************************************/
int Share_Interface_Information(LevelSet2D_Quad_Block *Soln_ptr,
				QuadTreeBlock_DataStructure &QuadTree,
				AdaptiveBlockResourceList &GlobalSolnBlockList,
				AdaptiveBlock2D_List &Soln_Block_List,
				LevelSet2D_Input_Parameters &Input_Parameters) {

#ifdef _RETRIEVE_DEBUG_
  char extension[256], output_file_name[256];
  char *output_file_name_ptr;
  ofstream dout;
  strcpy(output_file_name,"share_cpu");
  sprintf(extension,"%.6d",Soln_Block_List.ThisCPU);
  strcat(extension,".txt");
  strcat(output_file_name,extension);
  output_file_name_ptr = output_file_name;
  dout.open(output_file_name_ptr,ios::out);
  if (dout.bad()) return 1;
#endif

  int error_flag;
  int NI, NO, Number_of_Interfaces, new_Number_of_Interfaces;
  int *Active;
  Interface2D_List List;
  Interface2D temp_Interface;
  double dx;

  if (Soln_Block_List.Block[0].used == ADAPTIVEBLOCK2D_USED) {
    dx = 0.01*max(TOLER,min(fabs(Soln_ptr[0].Grid.Cell[2][2].Xc.x-Soln_ptr[0].Grid.Cell[1][2].Xc.x),
			    fabs(Soln_ptr[0].Grid.Cell[2][2].Xc.y-Soln_ptr[0].Grid.Cell[2][1].Xc.y)));
  } else {
    dx = ZERO;
  }
  dx = CFFC_Maximum_MPI(dx);

#ifdef _MPI_VERSION
  MPI::Intracomm local_comm;
  MPI::Group world_group = MPI::COMM_WORLD.Get_group();
  MPI::Group local_group;
  int *world_list;
  int *local_list;
  world_list = new int[QuadTree.Ncpu];
  local_list = new int[QuadTree.Ncpu];
  for (int ncpu = 0; ncpu < QuadTree.Ncpu; ncpu++) {
    world_list[ncpu] = ncpu;
    local_list[ncpu] = ncpu;
  }
#endif

  // Determine the number of interfaces that currently reside in 
  // the solution blocks owned by the current processor.
  Number_of_Interfaces = 0;
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      Number_of_Interfaces += Soln_ptr[nb].Interface_List.Ni;
    }
  }

  NI = 0;
  if (Number_of_Interfaces > 0) {

    // Allocate memory for interface variables for all interface 
    // information owned by the current processor.
    Active = new int[Number_of_Interfaces+1]; Active[0] = OFF;
    for (int n = 1; n <= Number_of_Interfaces; n++) Active[n] = ON;
    List.allocate(Number_of_Interfaces);

    // Copy all of the interface information owned by each individual 
    // blocks as one set of variables for each processor.  
    for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
      if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
	for (int ni = 1; ni <= Soln_ptr[nb].Interface_List.Ni; ni++) {
	  List[NI+ni].Copy(Soln_ptr[nb].Interface_List[ni]);
	}
	NI += Soln_ptr[nb].Interface_List.Ni;
      }
    }

  }

  // If there is more than one processor, the data collected on each 
  // must be shared with the others.  All of the processors must own all 
  // information for every interface even though an interface might not 
  // exist on the blocks owned by a specific processor.  This is 
  // necessary for solving the extension problems.
#ifdef _MPI_VERSION

  if (Number_of_Interfaces > 0) {

    // Deallocate the memory for the interface information on each block.
    for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
      if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
	if (Soln_ptr[nb].Interface_List.Ni) Soln_ptr[nb].Interface_List.deallocate();
      }
    }

    // Copy all of the interface information owned by each processor to 
    // each block owned by that proccessor.
    for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
      if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
	Soln_ptr[nb].Interface_List.allocate(Number_of_Interfaces);
	Soln_ptr[nb].Interface_List.Copy(List);
      }
    }

    // Deallocate memory for interface variables for all interface 
    // information owned by the current processor.
    delete []Active; Active = NULL;
    List.deallocate();

  }

  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

  // Determine the number of interfaces that currently reside on all
  // of the processors by call MPI_ALLReduce with the SUM operation.
  Number_of_Interfaces = CFFC_Summation_MPI(Number_of_Interfaces);

  // Allocate memory for interface variables for all of the interface 
  // information on all of the processors.
  Active = new int[Number_of_Interfaces+1]; Active[0] = OFF;
  for (int n = 1; n <= Number_of_Interfaces; n++) Active[n] = ON;
  List.allocate(Number_of_Interfaces);

  // Broadcast all of the interface information owned by each individual 
  // processor to the rest of the processors.
  NI = 0; NO = 0;
  for (int iCPU = 0; iCPU < QuadTree.Ncpu; iCPU++) {
    // Create the local list and communicator.
    for (int ncpu = 0; ncpu < QuadTree.Ncpu; ncpu++) {
      if (ncpu == 0) {
	local_list[ncpu] = iCPU;
      } else if (ncpu == iCPU) {
	local_list[ncpu] = 0;
      } else {
	local_list[ncpu] = ncpu;
      }
    }
    local_group = world_group.Incl(QuadTree.Ncpu,local_list);
    local_comm  = MPI::COMM_WORLD.Create(local_group);
    NO = NI;
    if (iCPU == GlobalSolnBlockList.ThisCPU) {
      for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
	if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
	  // Copy information to interface variables and fill message buffer.
	  for (int ni = 1; ni <= Soln_ptr[nb].Interface_List.Ni; ni++)
	    List[NI+ni].Copy(Soln_ptr[nb].Interface_List[ni]);
	  // Increment the number of interfaces.
	  NI += Soln_ptr[nb].Interface_List.Ni;
	  break;
	}
      }
    }
    // All reduce the number of interfaces with the maximum operator.
    CFFC_Barrier_MPI();
    NI = CFFC_Maximum_MPI(NI);

    //CFFC_Barrier_MPI();

    // Broadcast the interface list information.
    for (int ni = NO+1; ni <= NI; ni++)
      Broadcast_Interface(List[ni],local_comm,world_list[iCPU]);

    // Free local communicator and group.
    if (local_comm != MPI::COMM_NULL) local_comm.Free();
    local_group.Free();

  }

  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

#endif

  int next = 0, n_containted_points;
  double minimumdistance = MILLION;
  double hh, ht, th, tt;
  Vector2D Xmin, Xmax, Xp, fm;
  int clipped_flag = 0;

  if (Number_of_Interfaces > 0) {

    // Remove any splines that may be a subset of any other spline.
    for (int start = 1; start <= Number_of_Interfaces; start++) {
      if (Active[start]) {
	for (int trynext = 1; trynext <= Number_of_Interfaces; trynext++) {
	  if (trynext != start && Active[trynext]) {
	    n_containted_points = 0;
	    for (int np = 0; np < List[trynext].Spline.np; np++) {
	      for (int mp = 0; mp < List[start].Spline.np-1; mp++) {
		if (Point_On_Line(List[start].Spline.Xp[mp],List[start].Spline.Xp[mp+1],List[trynext].Spline.Xp[np])) {
		  n_containted_points++;
		  break;
		}
	      }
	    }
 	    if (n_containted_points == List[trynext].Spline.np) { Active[trynext] = OFF; }
	  }
	}
      }
    }

    // Concatenate (combine) splines if need be.
    for (int start = 1; start <= Number_of_Interfaces-1; start++) {
      next = start;
      minimumdistance = MILLION;

      if (Active[start]) {

	// Determine spline to concatenate with (next).
	for (int trynext = start+1; trynext <= Number_of_Interfaces; trynext++) {
	  if (Active[trynext]) {
	    hh = abs(List[start].Spline.Xp[0]                       - List[trynext].Spline.Xp[0]);
	    ht = abs(List[start].Spline.Xp[0]                       - List[trynext].Spline.Xp[List[trynext].Spline.np-1]);
	    th = abs(List[start].Spline.Xp[List[start].Spline.np-1] - List[trynext].Spline.Xp[0]);
	    tt = abs(List[start].Spline.Xp[List[start].Spline.np-1] - List[trynext].Spline.Xp[List[trynext].Spline.np-1]);
	    if (hh < minimumdistance || ht < minimumdistance || th < minimumdistance || tt < minimumdistance) {
	      next = trynext;
	      minimumdistance = min(min(hh,ht),min(th,tt));
	    }
	  }
	}

	if (Active[next] && next != start) {
	  if (abs(List[start].Spline.Xp[0]                       - List[next].Spline.Xp[0]) < dx &&
	      abs(List[start].Spline.Xp[List[start].Spline.np-1] - List[next].Spline.Xp[List[next].Spline.np-1]) < dx) { 
	    // Heads match and tails match.
	    temp_Interface.allocate(List[start].Spline.np+List[next].Spline.np-1);
	    temp_Interface.Spline.settype(List[start].Spline.type);
	    for (int j = 0; j < List[start].Spline.np; j++) {
	      temp_Interface.Spline.Xp[j] = List[start].Spline.Xp[j];
	      temp_Interface.Spline.tp[j] = SPLINE2D_POINT_NORMAL;
	      temp_Interface.Spline.bc[j] = List[start].Spline.bc[j];
	      temp_Interface.F[j] = List[start].F[j];
	    }
	    for (int j = 1; j < List[next].Spline.np; j++) {
	      temp_Interface.Spline.Xp[List[start].Spline.np-1+j] = List[next].Spline.Xp[List[next].Spline.np-1-j];
	      temp_Interface.Spline.tp[List[start].Spline.np-1+j] = SPLINE2D_POINT_NORMAL;
	      temp_Interface.Spline.bc[List[start].Spline.np-1+j] = List[next].Spline.bc[List[next].Spline.np-1-j];
	      temp_Interface.F[List[start].Spline.np-1+j] = List[next].F[List[next].Spline.np-1-j];
	    }
	    temp_Interface.Spline.pathlength();
	    Active[next] = OFF;

	  } else if (abs(List[start].Spline.Xp[0]                       - List[next].Spline.Xp[List[next].Spline.np-1]) < dx &&
		     abs(List[start].Spline.Xp[List[start].Spline.np-1] - List[next].Spline.Xp[0]) < dx) { 
	    // Head of 'start' matches tail of 'next' and tail of 'start' matches head of 'next'.
	    temp_Interface.allocate(List[start].Spline.np+List[next].Spline.np-1);
	    temp_Interface.Spline.settype(List[start].Spline.type);
	    for (int j = 0; j < List[start].Spline.np; j++) {
	      temp_Interface.Spline.Xp[j] = List[start].Spline.Xp[j];
	      temp_Interface.Spline.tp[j] = SPLINE2D_POINT_NORMAL;
	      temp_Interface.Spline.bc[j] = List[start].Spline.bc[j];
	      temp_Interface.F[j] = List[start].F[j];
	    }
	    for (int j = 1; j < List[next].Spline.np; j++) {
	      temp_Interface.Spline.Xp[List[start].Spline.np-1+j] = List[next].Spline.Xp[j];
	      temp_Interface.Spline.tp[List[start].Spline.np-1+j] = SPLINE2D_POINT_NORMAL;
	      temp_Interface.Spline.bc[List[start].Spline.np-1+j] = List[next].Spline.bc[j];
	      temp_Interface.F[List[start].Spline.np-1+j] = List[next].F[j];
	    }
	    temp_Interface.Spline.pathlength();
	    Active[next] = OFF;
  
	  } else if (abs(List[start].Spline.Xp[0] - List[next].Spline.Xp[0]) < dx) { 
	    // Head of 'start' matches head of 'next'.
	    temp_Interface.allocate(List[start].Spline.np+List[next].Spline.np-1);
	    temp_Interface.Spline.settype(List[start].Spline.type);
	    for (int j = List[start].Spline.np-1; j >= 0; j--) {
	      temp_Interface.Spline.Xp[List[start].Spline.np-1-j] = List[start].Spline.Xp[j];
	      temp_Interface.Spline.tp[List[start].Spline.np-1-j] = SPLINE2D_POINT_NORMAL;
	      temp_Interface.Spline.bc[List[start].Spline.np-1-j] = List[start].Spline.bc[j];
	      temp_Interface.F[List[start].Spline.np-1-j] = List[start].F[j];
	    }
	    for (int j = 1; j < List[next].Spline.np; j++) {
	      temp_Interface.Spline.Xp[List[start].Spline.np-1+j] = List[next].Spline.Xp[j];
	      temp_Interface.Spline.tp[List[start].Spline.np-1+j] = SPLINE2D_POINT_NORMAL;
	      temp_Interface.Spline.bc[List[start].Spline.np-1+j] = List[next].Spline.bc[j];
	      temp_Interface.F[List[start].Spline.np-1+j] = List[next].F[j];
	    }
	    temp_Interface.Spline.pathlength();
	    Active[next] = OFF;

	  } else if (abs(List[start].Spline.Xp[0] - List[next].Spline.Xp[List[next].Spline.np-1]) < dx) { 
	    // Head of 'start' matches tail of 'next'.
	    temp_Interface.allocate(List[start].Spline.np+List[next].Spline.np-1);
	    temp_Interface.Spline.settype(List[start].Spline.type);
	    for (int j = 0; j < List[next].Spline.np; j++) {
	      temp_Interface.Spline.Xp[j] = List[next].Spline.Xp[j];
	      temp_Interface.Spline.tp[j] = SPLINE2D_POINT_NORMAL;
	      temp_Interface.Spline.bc[j] = List[next].Spline.bc[j];
	      temp_Interface.F[j] = List[next].F[j];
	    }
	    for (int j = 1; j < List[start].Spline.np; j++) {
	      temp_Interface.Spline.Xp[List[next].Spline.np-1+j] = List[start].Spline.Xp[j];
	      temp_Interface.Spline.tp[List[next].Spline.np-1+j] = SPLINE2D_POINT_NORMAL;
	      temp_Interface.Spline.bc[List[next].Spline.np-1+j] = List[start].Spline.bc[j];
	      temp_Interface.F[List[next].Spline.np-1+j] = List[start].F[j];
	    }
	    temp_Interface.Spline.pathlength();
	    Active[next] = OFF;

	  } else if (abs(List[start].Spline.Xp[List[start].Spline.np-1] - List[next].Spline.Xp[0]) < dx) { 
	    // Tail of 'start' matches head of 'next'.
	    temp_Interface.allocate(List[start].Spline.np+List[next].Spline.np-1);
	    temp_Interface.Spline.settype(List[start].Spline.type);
	    for (int j = 0; j < List[start].Spline.np; j++) {
	      temp_Interface.Spline.Xp[j] = List[start].Spline.Xp[j];
	      temp_Interface.Spline.tp[j] = SPLINE2D_POINT_NORMAL;
	      temp_Interface.Spline.bc[j] = List[start].Spline.bc[j];
	      temp_Interface.F[j] = List[start].F[j];
	    }
	    for (int j = 1; j < List[next].Spline.np; j++) {
	      temp_Interface.Spline.Xp[List[start].Spline.np-1+j] = List[next].Spline.Xp[j];
	      temp_Interface.Spline.tp[List[start].Spline.np-1+j] = SPLINE2D_POINT_NORMAL;
	      temp_Interface.Spline.bc[List[start].Spline.np-1+j] = List[next].Spline.bc[j];
	      temp_Interface.F[List[start].Spline.np-1+j] = List[next].F[j];
	    }
	    temp_Interface.Spline.pathlength();
	    Active[next] = OFF;

	  } else if (abs(List[start].Spline.Xp[List[start].Spline.np-1] - List[next].Spline.Xp[List[next].Spline.np-1]) < dx) { 
	    // Tail of 'start' matches tail of 'next'.
	    temp_Interface.allocate(List[start].Spline.np+List[next].Spline.np-1);
	    temp_Interface.Spline.settype(List[start].Spline.type);
 	    for (int j = 0; j < List[start].Spline.np; j++) {
 	      temp_Interface.Spline.Xp[j] = List[start].Spline.Xp[j];
 	      temp_Interface.Spline.tp[j] = SPLINE2D_POINT_NORMAL;
 	      temp_Interface.Spline.bc[j] = List[start].Spline.bc[j];
	      temp_Interface.F[j] = List[start].F[j];
 	    }
 	    for (int j = 1; j < List[next].Spline.np; j++) {
 	      temp_Interface.Spline.Xp[List[start].Spline.np-1+j] = List[next].Spline.Xp[List[next].Spline.np-1-j];
 	      temp_Interface.Spline.tp[List[start].Spline.np-1+j] = SPLINE2D_POINT_NORMAL;
 	      temp_Interface.Spline.bc[List[start].Spline.np-1+j] = List[next].Spline.bc[List[next].Spline.np-1-j];
	      temp_Interface.F[List[start].Spline.np-1+j] = List[next].F[List[next].Spline.np-1-j];
 	    }
 	    temp_Interface.Spline.pathlength();
 	    Active[next] = OFF;

	  }

	  // If the 'next' spline has been concatenated with the 'start'
	  // spline then reset the 'start' spline information as the 
	  // concatenated spline information.  Start comparing with the
	  // other splines at the start of the list, just in case they 
	  // are now a candidate for concatination.
 	  if (Active[next] == OFF) {
 	    List[start].Copy(temp_Interface);
 	    temp_Interface.deallocate();
 	    start--;
 	  }

	}

      }

    }

    // Sort the interface into counter-clockwise order.
    for (int ni = 1; ni <= Number_of_Interfaces; ni++) {
      if (Active[ni]) {
	List[ni].Sort();
 	List[ni].Set_Normal_Velocity_Function();
      }
    }

    //////////////////////////
    int Nactive = 0;
    for (int ni = 1; ni <= Number_of_Interfaces; ni++) if (Active[ni]) Nactive++;
    if (Nactive != 1) {
#ifdef _RETRIEVE_DEBUG_
      dout << endl << " ACTIVE INTERFACES: ";
      for (int ni = 1; ni <= Number_of_Interfaces; ni++) {
 	if (Active[ni]) dout << endl << List[ni];
      }
      dout << endl << " ==================================================="; dout.flush();
      dout << endl << List;
#endif
      cout << "Boohoo... it's not working!" << endl;
      return 11111;
    }
    //////////////////////////


    // Recalculate the number of interfaces.
    new_Number_of_Interfaces = 0;
    for (int ni = 1; ni <= Number_of_Interfaces; ni++) {
      if (Active[ni]) new_Number_of_Interfaces++;
    }

    // Deallocate the memory for the interface information on each block.
    for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
      if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
	if (Soln_ptr[nb].Interface_List.Ni) Soln_ptr[nb].Interface_List.deallocate();
      }
    }

    // Copy all of the interface information owned by each processor to 
    // each block owned by that proccessor.
    for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
      if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
	Soln_ptr[nb].Interface_List.allocate(new_Number_of_Interfaces);
	NI = 1;
	for (int ni = 1; ni <= Number_of_Interfaces; ni++) {
	  if (Active[ni]) {
  	    Soln_ptr[nb].Interface_List[NI].Copy(List[ni]);
	    NI++;
	  }
	}
      }
    }

    // Deallocate memory for interface variables for all interface 
    // information owned by the current processor.
    delete []Active; Active = NULL;
    List.deallocate();

  }

#ifdef _MPI_VERSION
  delete []world_list; world_list = NULL;
  delete []local_list; local_list = NULL;
#endif

#ifdef _RETRIEVE_DEBUG_
//   for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
//     if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
//       for (int ni = 1; ni <= Soln_ptr[nb].Interface_List.Ni; ni++) {
//   	dout << endl << " INTERFACE #"  << ni << endl << Soln_ptr[nb].Interface_List[ni];
//       }
//     }
//   }

  dout.close();
#endif

  // Interfaces shared successfully.
  return 0;

}

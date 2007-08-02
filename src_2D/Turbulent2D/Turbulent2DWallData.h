/**********************************************************************
 * Turbulent2DWallData.h: Header file defining 2D particle solution   *
 *                        state classes.                              *
 **********************************************************************/

#ifndef _TURBULENT2D_WALLDATA_INCLUDED
#define _TURBULENT2D_WALLDATA_INCLUDED

// Include required C++ libraries.
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <cstring>

using namespace std;

// Include math macro, CFD, 2D vector, linked list, and gas constant header files.

#ifndef _MATH_MACROS_INCLUDED
#include "../Math/Math.h"
#endif // _MATH_MACROS_INCLUDED

#ifndef _CFD_INCLUDED
#include "../CFD.h"
#endif // _CFD_INCLUDED

#ifndef _VECTOR2D_INCLUDED
#include "../Math/Vector2D.h"
#endif //_VECTOR2D_INCLUDED

#ifdef _LINKEDLIST_INCLUDED
#include "../Math/LinkedList.h"
#endif // _LINKEDLIST_INCLUDED

#ifndef _GAS_CONSTANTS_INCLUDED
#include "../Physics/GasConstants.h"
#endif // _GAS_CONSTANTS_INCLUDED

/*!
 * Class: Turbulent2DWallData
 *
 * @brief Class defined to store and handle wall distance calculations
 * (ie, ywall, yplus, ...).
 *
 * \verbatim
 * Member functions
 *     ywall  -- Distance to nearest solid wall.
 *     Xwall  -- Location of the nearest solid wall.
 *     nwall  -- Normal at the nearest solid wall.
 *     yplus  -- Dimensionless distance to nearest wall.
 *     utau   -- Friction velocity at nearest wall.
 *     tauw   -- Shear stress at nearest wall.
 *     BCwall -- Nearest wall boundary condition.
 *     vinj   -- Dimensional wall injection speed.
 *     vwplus -- Dimensionless wall injection speed.
 *
 * Member operators
 *      W -- a primitive solution state
 *      c -- a scalar (double)
 *
 * cout << W; (output function)
 * cin  >> W; (input function)
 * \endverbatim
 */
class Turbulent2DWallData {
 private:
 public:
  double    ywall; //!< Distance to nearest solid wall.
  Vector2D  Xwall; //!< Location of the nearest solid wall.
  Vector2D  nwall; //!< Normal at the nearest solid wall.
  double    yplus; //!< Dimensionless distance to nearest wall.
  double     utau; //!< Friction velocity at nearest wall.
  double     tauw; //!< Shear stress at nearest wall.
  int      BCwall; //!< Wall boundary condition.
  double     vinj; //!< Dimensionlal injection speed at the nearest burning surface.
  double   vwplus; //!< Dimensionless injection speed at the nearest burning surface.

  //! Creation constructor.
  Turbulent2DWallData(void) {
    ywall = 1e70; Xwall = Vector2D_ZERO; nwall = Vector2D_ZERO;
    yplus = ZERO; utau = ZERO; tauw = ZERO;
    BCwall = BC_WALL_VISCOUS_HEATFLUX;
    vinj = ZERO; vwplus = ZERO;
  }

  //! Copy constructor.
  Turbulent2DWallData(const Turbulent2DWallData &W) {
    ywall = W.ywall; Xwall = W.Xwall; nwall = W.nwall;
    yplus = W.yplus; utau = W.utau; tauw = W.tauw;
    BCwall = W.BCwall;
    vinj = W.vinj;
    vwplus = W.vwplus;
  }

  //! Destructor.
  ~Turbulent2DWallData(void) { };

  //! Copy function.
  void copy(const Turbulent2DWallData &W) {
    ywall = W.ywall; Xwall = W.Xwall; nwall = W.nwall;
    yplus = W.yplus; utau = W.utau; tauw = W.tauw;
    BCwall = W.BCwall;
    vinj = W.vinj;
    vwplus = W.vwplus;
  }

  //@{ @name Input-output operators.
  friend ostream &operator << (ostream &out_file, const Turbulent2DWallData &W);
  friend istream &operator >> (istream &in_file, Turbulent2DWallData &W);
  //@}

};

//! Turbulent2DWallData -- Output operator.
inline ostream &operator << (ostream &out_file, const Turbulent2DWallData &W) {
  out_file.setf(ios::scientific);
  out_file << W.ywall << W.Xwall << W.nwall << " ";
  out_file.unsetf(ios::scientific);
  out_file << W.BCwall << " " << W.vinj << " " << W.vwplus;
  out_file.setf(ios::scientific);
  out_file << W.yplus << " " << W.utau << " " << W.tauw << endl;
  out_file.unsetf(ios::scientific);
  return out_file;
}

//! Turbulent2DWallData -- Input operator.
inline istream &operator >> (istream &in_file, Turbulent2DWallData &W) {
  in_file.setf(ios::skipws);
  in_file >> W.ywall >> W.Xwall >> W.nwall >> W.BCwall >> W.vinj >> W.vwplus >> W.yplus >> W.utau >> W.tauw;
  in_file.unsetf(ios::skipws);
  return in_file;
}

/**********************************************************************
 * The next three routines determine the normal distance to the       *
 * nearest wall, the location of the nearest wall, the normal at the  *
 * nearest wall, and the boundary condition at the nearest wall for   *
 * each and every cell.                                               *
 **********************************************************************/

/**********************************************************************
 * Routine: Determine_Wall_Distance                                   *
 *                                                                    *
 * This routine determines the normal distance to the wall for each   *
 * cell of a multiblock body-fitted quadrilateral mesh.  A direct     *
 * (exhaustive) search is performed.  The wall location, the inward   *
 * normal at the wall, and the wall boundary condition type is also   *
 * found and stored.                                                  *
 *                                                                    *
 **********************************************************************/
template <class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int Determine_Wall_Distance(Quad_Soln_Block *Soln_ptr,
			    QuadTreeBlock_DataStructure &QuadTree,
			    AdaptiveBlock2D_List &LocalSolnBlockList,
			    Quad_Soln_Input_Parameters &IP) {

  // Exit immediately if not a turbulent flow.
  if (IP.FlowType == FLOWTYPE_INVISCID || IP.FlowType == FLOWTYPE_LAMINAR) return 0;

  int error_flag;
  int buffer_size;
  int ilow, iup, jlow, jup, ni, nj;
  double y_wall, y_wall_temp;
  double *x_buffer;
  Vector2D x_cell, X_wall, n_wall;
  int BC_wall;

#ifdef _MPI_VERSION
  int wall_buffer_size, wall_iCPU;
  double *wall_buffer;
#endif

#ifdef _MPI_VERSION
  wall_buffer_size = 5;
  wall_buffer = new double[wall_buffer_size];
#endif

  // Compute the distance to the nearest wall for the centroid of
  // every cell.
  for (int iCPU = 0; iCPU < QuadTree.Ncpu; iCPU++) { // Loop over available processors.
    for (int iBLK = 0; iBLK < QuadTree.Nblk; iBLK++) { // Loop over available blocks.
      if (QuadTree.Blocks[iCPU][iBLK] != NULL) {
	if (QuadTree.Blocks[iCPU][iBLK]->block.used) {

	  // Determine dimensions of block of interest.
	  //ilow = QuadTree.Blocks[iCPU][iBLK]->block.info.dimen.ghost;
	  //jlow = QuadTree.Blocks[iCPU][iBLK]->block.info.dimen.ghost;
	  //iup = ilow + QuadTree.Blocks[iCPU][iBLK]->block.info.dimen.i - 1;
	  //jup = jlow + QuadTree.Blocks[iCPU][iBLK]->block.info.dimen.j - 1;
	  ilow = 0;
	  jlow = 0;
	  //iup = 2*QuadTree.Blocks[iCPU][iBLK]->block.info.dimen.ghost + QuadTree.Blocks[iCPU][iBLK]->block.info.dimen.i - 1;
	  //jup = 2*QuadTree.Blocks[iCPU][iBLK]->block.info.dimen.ghost + QuadTree.Blocks[iCPU][iBLK]->block.info.dimen.j - 1;
	  if (LocalSolnBlockList.ThisCPU == iCPU) {
	    iup = 2*LocalSolnBlockList.Block[iBLK].info.dimen.ghost +
	          LocalSolnBlockList.Block[iBLK].info.dimen.i - 1;
	    jup = 2*LocalSolnBlockList.Block[iBLK].info.dimen.ghost +
	          LocalSolnBlockList.Block[iBLK].info.dimen.j - 1;
	  } else {
	    iup = 0;
	    jup = 0;
	  }
	  iup = CFDkit_Maximum_MPI(iup);
	  jup = CFDkit_Maximum_MPI(jup);

	  // Allocate buffer for storing cell centers.
	  ni = iup - ilow + 1;
	  nj = jup - jlow + 1;
	  x_buffer = new double[2*ni*nj];

	  // Get coordinates of cell centers.
	  if (LocalSolnBlockList.ThisCPU == iCPU) {
	    buffer_size = 0;
	    for (int j = jlow; j <= jup; j++) {
	      for (int i = ilow; i <= iup; i++) {
		x_buffer[buffer_size  ] = Soln_ptr[iBLK].Grid.Cell[i][j].Xc.x;
		x_buffer[buffer_size+1] = Soln_ptr[iBLK].Grid.Cell[i][j].Xc.y;
		buffer_size = buffer_size + 2;
	      }
	    }
	  }

#ifdef _MPI_VERSION
	  buffer_size = 2*ni*nj;
	  MPI::COMM_WORLD.Bcast(x_buffer,buffer_size,MPI::DOUBLE,iCPU);
#endif

	  // Calculate normal distance to the wall.
	  buffer_size = 0;
	  for (int j = jlow; j <= jup; j++) {
	    for (int i = ilow; i <= iup; i++) {
	      x_cell.x = x_buffer[buffer_size  ];
	      x_cell.y = x_buffer[buffer_size+1];
	      buffer_size = buffer_size + 2;
	      error_flag = Wall_Distance(Soln_ptr,LocalSolnBlockList,x_cell,X_wall,n_wall,y_wall,BC_wall);
	      if (error_flag) return error_flag;
#ifdef _MPI_VERSION
	      // The following routine will reduce the y_wall variable
	      // to be the minimum value of y_wall on the p processors.
	      // It will also ensure that the X_wall variable is
	      // consistent with that distance.
	      y_wall_temp = CFDkit_Minimum_MPI(y_wall);
	      if (fabs(y_wall - y_wall_temp) < TOLER) {
		wall_buffer[0] = X_wall.x;
		wall_buffer[1] = X_wall.y;
		wall_buffer[2] = n_wall.x;
		wall_buffer[3] = n_wall.y;
		wall_buffer[4] = double(BC_wall);
		wall_iCPU = LocalSolnBlockList.ThisCPU;
	      } else {
		wall_buffer[0] = ZERO;
		wall_buffer[1] = ZERO;
		wall_buffer[2] = ZERO;
		wall_buffer[3] = ZERO;
		wall_buffer[4] = ZERO;
		wall_iCPU = 0;
	      }
	      wall_iCPU = CFDkit_Maximum_MPI(wall_iCPU);
	      MPI::COMM_WORLD.Bcast(wall_buffer,wall_buffer_size,MPI::DOUBLE,wall_iCPU);
	      X_wall = Vector2D(wall_buffer[0],wall_buffer[1]);
	      n_wall = Vector2D(wall_buffer[2],wall_buffer[3]);
	      y_wall = y_wall_temp;
	      BC_wall = int(wall_buffer[4]);
#endif
	      if (LocalSolnBlockList.ThisCPU == iCPU) {
		Soln_ptr[iBLK].Wall[i][j].ywall = y_wall;
		Soln_ptr[iBLK].Wall[i][j].Xwall = X_wall;
		Soln_ptr[iBLK].Wall[i][j].nwall = n_wall;
		Soln_ptr[iBLK].Wall[i][j].BCwall = BC_wall;
	      }
	    }
	  }

	  // Deallocate buffer for storing cell centers.
	  delete []x_buffer; x_buffer = NULL;

	}
      }
    }
  }

#ifdef _MPI_VERSION
  delete []wall_buffer; wall_buffer = NULL;
#endif

  // Distance to wall calculation complete.
  return 0;

}

/**********************************************************************
 * Routine: Wall_Distance                                             *
 *                                                                    *
 * This routine determines the distance to the nearest wall for each  *
 * cell centroid of a 2D array of quadrilateral solution blocks.  The *
 * wall location, inward normal at the wall, and the wall boundary    *
 * condition are also found and stored.                               *
 *                                                                    *
 **********************************************************************/
template <class Quad_Soln_Block>
int Wall_Distance(Quad_Soln_Block *Soln_ptr,
		  AdaptiveBlock2D_List &LocalSolnBlockList,
		  const Vector2D X_cell,
		  Vector2D &X_wall,
		  Vector2D &n_wall,
		  double &y_wall,
		  int &BC_wall) {

  int error_flag;
  double y_wall_temp;
  Vector2D X_wall_temp, n_wall_temp;
  int BC_wall_temp;

  y_wall = 1e70;
  y_wall_temp = 1e70;

  for (int nb = 0 ; nb < LocalSolnBlockList.Nblk; nb++) {
    if (LocalSolnBlockList.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      error_flag = Wall_Distance(Soln_ptr[nb],X_cell,X_wall_temp,n_wall_temp,y_wall_temp,BC_wall_temp);
      if (error_flag) return error_flag;
      if (y_wall_temp < y_wall) {
	y_wall = y_wall_temp;
	X_wall = X_wall_temp;
	n_wall = n_wall_temp;
	BC_wall = BC_wall_temp;
      }
    }
  }

  // Wall distance successfully found.
  return 0;

}

/**********************************************************************
 * Routine: Wall_Distance                                             *
 *                                                                    *
 * This routine determines the distance to the nearest wall for a     *
 * specific cell centroid.  The wall location, inward normal at the   *
 * wall, and the wall boundary condition are also found and stored.   *
 *                                                                    *
 **********************************************************************/
template <class Quad_Soln_Block>
int Wall_Distance(Quad_Soln_Block &SolnBlk,
		  const Vector2D X_cell,
		  Vector2D &X_wall,
		  Vector2D &n_wall,
		  double &y_wall,
		  int &BC_wall) {

  // Initialize y_wall.
  y_wall = 1e70;
  X_wall = Vector2D_ZERO;
  n_wall = Vector2D_ZERO;
  BC_wall = BC_WALL_VISCOUS_HEATFLUX;

  // Check West boundary.
  //for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
  for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
    if (SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS ||
	SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
	SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX ||
	SolnBlk.Grid.BCtypeW[j] == BC_BURNING_SURFACE ||
	SolnBlk.Grid.BCtypeW[j] == BC_MASS_INJECTION) {
      if (abs(SolnBlk.Grid.xfaceW(SolnBlk.ICl,j)-X_cell) < y_wall) {
	y_wall = abs(SolnBlk.Grid.xfaceW(SolnBlk.ICl,j)-X_cell);
	X_wall = SolnBlk.Grid.xfaceW(SolnBlk.ICl,j);
	n_wall = -SolnBlk.Grid.nfaceW(SolnBlk.ICl,j);
	BC_wall = SolnBlk.Grid.BCtypeW[j];
      }
    }
  }

  // Check East boundary.
  //for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
  for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
    if (SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS ||
	SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
	SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX ||
	SolnBlk.Grid.BCtypeE[j] == BC_BURNING_SURFACE ||
	SolnBlk.Grid.BCtypeE[j] == BC_MASS_INJECTION) {
      if (abs(SolnBlk.Grid.xfaceE(SolnBlk.ICu,j)-X_cell) < y_wall) {
	y_wall = abs(SolnBlk.Grid.xfaceE(SolnBlk.ICu,j)-X_cell);
	X_wall = SolnBlk.Grid.xfaceE(SolnBlk.ICu,j);
	n_wall = -SolnBlk.Grid.nfaceE(SolnBlk.ICu,j);
	BC_wall = SolnBlk.Grid.BCtypeE[j];
      }
    }
  }

  // Check South boundary.
  //for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
  for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
    if (SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS ||
	SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
	SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX ||
	SolnBlk.Grid.BCtypeS[i] == BC_BURNING_SURFACE ||
	SolnBlk.Grid.BCtypeS[i] == BC_MASS_INJECTION) {
      if (abs(SolnBlk.Grid.xfaceS(i,SolnBlk.JCl)-X_cell) < y_wall) {
	y_wall = abs(SolnBlk.Grid.xfaceS(i,SolnBlk.JCl)-X_cell);
	X_wall =  SolnBlk.Grid.xfaceS(i,SolnBlk.JCl);
	n_wall = -SolnBlk.Grid.nfaceS(i,SolnBlk.JCl);
	BC_wall = SolnBlk.Grid.BCtypeS[i];
      }
    }
  }

  // Check North boundary.
  //for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
  for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
    if (SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS ||
	SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
	SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX ||
	SolnBlk.Grid.BCtypeN[i] == BC_BURNING_SURFACE ||
	SolnBlk.Grid.BCtypeN[i] == BC_MASS_INJECTION) {
      if (abs(SolnBlk.Grid.xfaceN(i,SolnBlk.JCu)-X_cell) < y_wall) {
	y_wall = abs(SolnBlk.Grid.xfaceN(i,SolnBlk.JCu)-X_cell);
	X_wall = SolnBlk.Grid.xfaceN(i,SolnBlk.JCu);
	n_wall = -SolnBlk.Grid.nfaceN(i,SolnBlk.JCu);
	BC_wall = SolnBlk.Grid.BCtypeN[i];
      }
    }
  }

  // Distance to wall computed successfully.
  return 0;

}

/**********************************************************************
 * The next six routines determine the dimensionless distance to the  *
 * nearest wall.  Various methods for doing so are included based on  *
 * different ways of determining the friction velocity/wall shear     *
 * stress.  The first method calculates utau and tauw based on the    *
 * local cell solution states and gradients.  This method has not     *
 * fully tested and is theoretically inaccurate.  The second method   *
 * calculates utau and tauw only at cells with a no-slip boundary     *
 * condition and then broadcasts these values to the associated cells *
 * for the y-plus calculation.  This method may be useful but it has  *
 * not been fully tested.  The third method iterates on the law of    *
 * wall to determine utau and tauw.  This method has been fully       *
 * tested and works well.  A check to determine if the cell is in the *
 * viscous sublayer improves the accuracy and simplifies the          *
 * calculation.  The final method provides the exact utau found at    *
 * the wall for a (horizontal) pipe.  This should only be used for    *
 * pipe flows and only for testing purposes.                          *
 **********************************************************************/

/**********************************************************************
 * Routine: Dimensionless_Wall_Distance                               *
 *                                                                    *
 * This routine determines the nondimensional normal wall distance    *
 * (y-plus) for each cell of a multiblock body-fitted mesh.  The      *
 * dimensional wall distance and the location of the nearest wall     *
 * must have been precomputed.  To calculate the nondimensional wall  *
 * distance, the shear stress at the wall (in the coordinate frame of *
 * the wall) must be determined and broadcasted to the associated     *
 * cells.  The friction velocity and y-plus can then be computed.     *
 *                                                                    *
 **********************************************************************/
template <class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int Dimensionless_Wall_Distance(Quad_Soln_Block *Soln_ptr,
				//QuadTreeBlock_DataStructure &QuadTree,
				AdaptiveBlock2D_List &LocalSolnBlockList,
				Quad_Soln_Input_Parameters &IP) {

  // Exit immediately if not a turbulent flow.
  if (IP.FlowType == FLOWTYPE_INVISCID || IP.FlowType == FLOWTYPE_LAMINAR) return 0;

  int error_flag;
//   int buffer_size;
//   LinkedList<double> LL_buffer;
//   double *buffer;

  switch(IP.i_Friction_Velocity) {
  case FRICTION_VELOCITY_LOCAL_SHEAR_STRESS :
    // Determine y-plus based on the local shear stress.
    for (int nb = 0; nb < LocalSolnBlockList.Nblk; nb++) {
      if (LocalSolnBlockList.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
	error_flag = Dimensionless_Wall_Distance(Soln_ptr[nb]);
	if (error_flag) return error_flag;
      }
    }
    break;

  case FRICTION_VELOCITY_WALL_SHEAR_STRESS :
    // Determine y-plus based on the shear stress at the wall.

//     // Initialize buffer size.
//     buffer_size = 0;
//     // Cycle through the processors to create, load, and broadcast the
//     // message buffer containing the wall shear stresses defined at all
//     // no-slip type boundary cells.  Each processor should then apply
//     // the information accordingly.
//     for (int iCPU = 0; iCPU < QuadTree.Ncpu; iCPU++) {
//       // Compile a linked list of the shear stress on the processor.
//       if (LL_buffer.np) LL_buffer.deallocate();
//       if (LocalSolnBlockList.ThisCPU == iCPU) {
// 	for (int nb = 0; nb < LocalSolnBlockList.Nblk; nb++) {
// 	  if (LocalSolnBlockList.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
// 	    error_flag = Compile_Wall_Shear_Stress_List(Soln_ptr[nb],LL_buffer);
// 	    if (error_flag) return error_flag;
// 	    buffer_size = LL_buffer.np;
// 	  }
// 	}
//       }
//       // Create the message buffer based on the maximum buffer size.
//       buffer_size = CFDkit_Maximum_MPI(buffer_size);
//       buffer = new double[buffer_size];
//       // Load message buffer.
//       if (LocalSolnBlockList.ThisCPU == iCPU) {
// 	for (int n = 0; n < buffer_size; n += 3) {
// 	  buffer[n  ] = LL_buffer[n  ];
// 	  buffer[n+1] = LL_buffer[n+1];
// 	  buffer[n+2] = LL_buffer[n+2];
// 	}
//       } else {
// 	for (int n = 0; n < buffer_size; n++) buffer[n] = ZERO;
//       }
//       // Broadcast message buffer.
// #ifdef _MPI_VERSION
//       MPI::COMM_WORLD.Bcast(buffer,buffer_size,MPI::DOUBLE,iCPU);
// #endif
//       // Apply wall shear stress information.
//       for (int nb = 0; nb < LocalSolnBlockList.Nblk; nb++) {
// 	if (LocalSolnBlockList.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
// 	  for (int n = 0; n < buffer_size; n += 3) {
// 	    error_flag = Dimensionless_Wall_Distance(Soln_ptr[nb],
// 						     Vector2D(buffer[n],buffer[n+1]),
// 						     buffer[n+2]);
// 	    if (error_flag) return error_flag;
// 	  }
// 	}
//       }
//       // Deallocate message buffer.
//       delete []buffer; buffer = NULL; buffer_size = 0;
//     }
    break;

  case FRICTION_VELOCITY_ITERATIVE :
    // Use an iterative method to determine the friction velocity.
    for (int nb = 0; nb < LocalSolnBlockList.Nblk; nb++) {
      if (LocalSolnBlockList.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
	error_flag = Iterative_Dimensionless_Wall_Distance(Soln_ptr[nb]);
	if (error_flag) return error_flag;
      }
    }
    break;

  case FRICTION_VELOCITY_PIPE :
    // Set the friction velocity to the exact friction velocity for
    // a turbulent pipe flow.
    for (int nb = 0; nb < LocalSolnBlockList.Nblk; nb++) {
      if (LocalSolnBlockList.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
	error_flag = Pipe_Dimensionless_Wall_Distance(Soln_ptr[nb],IP);
	if (error_flag) return error_flag;
      }
    }
    break;

  };

  // Dimensionless distance to wall calculation complete.
  return 0;

}

/**********************************************************************
 * Routine: Compile_Wall_Shear_Stress_List                            *
 *                                                                    *
 * This routine compiles a list of the shear stress at all solid or   *
 * transpiring walls in the coordinate frame of the wall.  Note that  *
 * the list includes the wall face location, required to match the    *
 * wall shear stress with the wall distance.                          *
 *                                                                    *
 **********************************************************************/
template <class Quad_Soln_Block>
int Compile_Wall_Shear_Stress_List(Quad_Soln_Block &SolnBlk,
				   LinkedList<double> &LL_buffer) {

  int error_flag;

  for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
    // North block boundary.
    if (SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS ||
	SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
	SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX ||
	SolnBlk.Grid.BCtypeN[i] == BC_BURNING_SURFACE ||
	SolnBlk.Grid.BCtypeN[i] == BC_MASS_INJECTION) {
      SolnBlk.Wall[i][SolnBlk.JCu].tauw = WallShearStress(SolnBlk.W[i][SolnBlk.JCu],
							  SolnBlk.Grid.Cell[i][SolnBlk.JCu].Xc,
							  SolnBlk.Grid.nodeNE(i,SolnBlk.JCu).X,
							  SolnBlk.Grid.nodeNW(i,SolnBlk.JCu).X,
							  -SolnBlk.Grid.nfaceN(i,SolnBlk.JCu));
      LL_buffer.add(SolnBlk.Grid.xfaceN(i,SolnBlk.JCu).x);
      LL_buffer.add(SolnBlk.Grid.xfaceN(i,SolnBlk.JCu).y);
      LL_buffer.add(SolnBlk.Wall[i][SolnBlk.JCu].tauw);
    }
    // South block boundary.
    if (SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS ||
	SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
	SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX ||
	SolnBlk.Grid.BCtypeS[i] == BC_BURNING_SURFACE ||
	SolnBlk.Grid.BCtypeS[i] == BC_MASS_INJECTION) {
      SolnBlk.Wall[i][SolnBlk.JCl].tauw = WallShearStress(SolnBlk.W[i][SolnBlk.JCl],
							  SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc,
							  SolnBlk.Grid.nodeSW(i,SolnBlk.JCl).X,
							  SolnBlk.Grid.nodeSE(i,SolnBlk.JCl).X,
							  -SolnBlk.Grid.nfaceS(i,SolnBlk.JCl));
      LL_buffer.add(SolnBlk.Grid.xfaceS(i,SolnBlk.JCl).x);
      LL_buffer.add(SolnBlk.Grid.xfaceS(i,SolnBlk.JCl).y);
      LL_buffer.add(SolnBlk.Wall[i][SolnBlk.JCl].tauw);
    }
  }
  for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
    // East block boundary.
    if (SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS ||
	SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
	SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX ||
	SolnBlk.Grid.BCtypeE[j] == BC_BURNING_SURFACE ||
	SolnBlk.Grid.BCtypeE[j] == BC_MASS_INJECTION) {
      SolnBlk.Wall[SolnBlk.ICu][j].tauw = WallShearStress(SolnBlk.W[SolnBlk.ICu][j],
							  SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc,
							  SolnBlk.Grid.nodeSE(SolnBlk.ICu,j).X,
							  SolnBlk.Grid.nodeNE(SolnBlk.ICu,j).X,
							  -SolnBlk.Grid.nfaceE(SolnBlk.ICu,j));
      LL_buffer.add(SolnBlk.Grid.xfaceE(SolnBlk.ICu,j).x);
      LL_buffer.add(SolnBlk.Grid.xfaceE(SolnBlk.ICu,j).y);
      LL_buffer.add(SolnBlk.Wall[SolnBlk.ICu][j].tauw);
    }
    // West block boundary.
    if (SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS ||
	SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
	SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX ||
	SolnBlk.Grid.BCtypeW[j] == BC_BURNING_SURFACE ||
	SolnBlk.Grid.BCtypeW[j] == BC_MASS_INJECTION) {
      SolnBlk.Wall[SolnBlk.ICl][j].tauw = WallShearStress(SolnBlk.W[SolnBlk.ICl][j],
							  SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc,
							  SolnBlk.Grid.nodeSW(SolnBlk.ICl,j).X,
							  SolnBlk.Grid.nodeNW(SolnBlk.ICl,j).X,
							  -SolnBlk.Grid.nfaceW(SolnBlk.ICl,j));
      LL_buffer.add(SolnBlk.Grid.xfaceW(SolnBlk.ICl,j).x);
      LL_buffer.add(SolnBlk.Grid.xfaceW(SolnBlk.ICl,j).y);
      LL_buffer.add(SolnBlk.Wall[SolnBlk.ICl][j].tauw);
    }
  }

  // List of wall shear stresses compiled successfully.
  return 0;

}

/**********************************************************************
 * Routine: Dimensionless_Wall_Distance                               *
 *                                                                    *
 * This routine calculates the dimensionless wall distance (y-plus)   *
 * based on the shear stress determined at the wall.                  *
 *                                                                    *
 **********************************************************************/
template <class Quad_Soln_Block>
int Dimensionless_Wall_Distance(Quad_Soln_Block &SolnBlk,
				const Vector2D Xwall,
				const double tauw) {

  int error_flag;

  for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
      if (abs(Xwall - SolnBlk.Wall[i][j].Xwall) < NANO) {
	SolnBlk.Wall[i][j].tauw = tauw;
	SolnBlk.Wall[i][j].utau = sqrt(fabs(tauw)/SolnBlk.W[i][j].rho);
	SolnBlk.Wall[i][j].yplus = SolnBlk.Wall[i][j].utau*SolnBlk.Wall[i][j].ywall/SolnBlk.W[i][j].nu();
      }
    }
  }

  // Dimensionless distance to wall calculation complete.
  return 0;

}

/**********************************************************************
 * Routine: Dimensionless_Wall_Distance                               *
 *                                                                    *
 * This routine calculates the dimensionless wall distance (y-plus)   *
 * based on the shear stress determined locally.                      *
 *                                                                    *
 **********************************************************************/
template <class Quad_Soln_Block>
int Dimensionless_Wall_Distance(Quad_Soln_Block &SolnBlk) {

  int error_flag;

  for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
      Linear_Reconstruction_LeastSquares(SolnBlk,i,j,LIMITER_ONE);
      SolnBlk.Wall[i][j].tauw = ShearStress(SolnBlk.W[i][j],
					    SolnBlk.dWdx[i][j],
					    SolnBlk.dWdy[i][j],
					    SolnBlk.Wall[i][j].nwall);
      SolnBlk.Wall[i][j].utau = sqrt(fabs(SolnBlk.Wall[i][j].tauw)/SolnBlk.W[i][j].rho);
      SolnBlk.Wall[i][j].yplus = SolnBlk.Wall[i][j].utau*SolnBlk.Wall[i][j].ywall/SolnBlk.W[i][j].nu();
    }
  }

  // Dimensionless distance to wall calculation complete.
  return 0;

}

/**********************************************************************
 * Routine: Iterative_Dimensionless_Wall_Distance                     *
 *                                                                    *
 * This routine calculates the dimensionless wall distance (y-plus)   *
 * based on an iterative solution for the friction velocity based on  *
 * the law of wall.  A  check to determine if the cell is in the      *
 * viscous sublayer improves the accuracy and simplifies the          *
 * calculation.                                                       *
 *                                                                    *
 **********************************************************************/
template <class Quad_Soln_Block>
int Iterative_Dimensionless_Wall_Distance(Quad_Soln_Block &SolnBlk) {

  int n, nmax = 25;
  double ywall, nu, u, utau, E, f, df;//, dE, A, kappav, dkappav, Xi, dXi;
  Vector2D that;

  for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
    for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
      // Get the dimensional wall distance.
      ywall = SolnBlk.Wall[i][j].ywall;

      if (ywall > 10000.0) {
      // Set the final friction velocity.
      SolnBlk.Wall[i][j].utau = MILLION;
      // Set the final wall shear stress.
      SolnBlk.Wall[i][j].tauw = MILLION;
      // Determine y-plus for the current cell.
      SolnBlk.Wall[i][j].yplus = MILLION;

      } else {

      // Get the local viscosity.
      nu = SolnBlk.W[i][j].nu();

      if (i < SolnBlk.ICl && (SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS ||
			      SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
			      SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX ||
			      SolnBlk.Grid.BCtypeW[j] == BC_BURNING_SURFACE ||
			      SolnBlk.Grid.BCtypeW[j] == BC_MASS_INJECTION)) {
	utau = ZERO;

      } else if (i > SolnBlk.ICu && (SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS ||
				     SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
				     SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX ||
				     SolnBlk.Grid.BCtypeE[j] == BC_BURNING_SURFACE ||
				     SolnBlk.Grid.BCtypeE[j] == BC_MASS_INJECTION)) {
	utau = ZERO;

      } else if (j < SolnBlk.JCl && (SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS ||
				     SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
				     SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX ||
				     SolnBlk.Grid.BCtypeS[i] == BC_BURNING_SURFACE ||
				     SolnBlk.Grid.BCtypeS[i] == BC_MASS_INJECTION)) {
	utau = ZERO;

      } else if (j > SolnBlk.JCu && (SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS ||
				     SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
				     SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX ||
				     SolnBlk.Grid.BCtypeN[i] == BC_BURNING_SURFACE ||
				     SolnBlk.Grid.BCtypeN[i] == BC_MASS_INJECTION)) {
	utau = ZERO;

      } else {

	// Get the component of the local velocity tangent to the wall.
	that.x = SolnBlk.Wall[i][j].nwall.y;
	that.y = SolnBlk.Wall[i][j].nwall.x;
	u = fabs(SolnBlk.W[i][j].v*that);
	//if (fabs(u) < TOLER) { u = SolnBlk.W[i][j].v.abs(); if (fabs(u) < TOLER) return 0; }
	//if (fabs(u) < TOLER) u = 1.0;//0.001;
	//if (fabs(u) > TOLER) {
	if (u < TOLER) {
	  utau = ZERO;
	} else {
//  	  if (SolnBlk.Wall[i][j].BCwall != BC_BURNING_SURFACE) {
	    // Set the boundary layer constant.
	    E = exp(SolnBlk.W[i][j].von_karman*SolnBlk.W[i][j].C);
	    // Estimate the friction velocity based on the viscous sublayer
	    // profile.
	    utau = sqrt(nu*u/ywall);
	    // Determine the initial y-plus for the cell.
	    SolnBlk.Wall[i][j].yplus = utau*ywall/nu;
	    // If the initial y-plus locates the current cell within the 
	    // viscous sublayer then the y-plus has been found.  Otherwise the 
	    // current cell is in (or beyond) the log layer.  Iterate on the
	    // law of wall to determine the friction velocity.
	    if (SolnBlk.Wall[i][j].yplus > SolnBlk.W[i][j].yplus_o) {
	      n = 0;
	      do {
		f = utau - SolnBlk.W[i][j].von_karman*u/log(E*SolnBlk.Wall[i][j].yplus);
		df = ONE + SolnBlk.W[i][j].von_karman*(u/utau)/sqr(log(E*SolnBlk.Wall[i][j].yplus));
		utau -= f/df;
		SolnBlk.Wall[i][j].yplus = utau*ywall/nu;
		n++;
	      } while(fabs(f) >= 0.001 && n <= nmax);
	    }
//  	  } else {
// 	    utau = sqrt(nu*u/ywall);
// 	    cout << endl << " utau = " << utau;
// 	    SolnBlk.Wall[i][j].yplus = utau*ywall/nu;
// 	    cout << endl << " yplus = " << SolnBlk.Wall[i][j].yplus;
// 	    cout << endl << " vinj = " << SolnBlk.Wall[i][j].vinj*0.5;
// 	    SolnBlk.Wall[i][j].vwplus = SolnBlk.Wall[i][j].vinj*0.5/utau;
// 	    cout << endl << " vwplus = " << SolnBlk.Wall[i][j].vwplus;
// 	    A = 0.50*(SolnBlk.W[i][j].C - (3.0*SolnBlk.W[i][j].sigma_omega - 2.0)/(2.0*SolnBlk.W[i][j].sigma_omega*SolnBlk.W[i][j].von_karman));
// 	    cout << endl << " A = " << A;
// 	    n = 0;
// 	    do {
// 	      cout << endl << " n = " << n;
// 	      Xi = A + 0.25*log(SolnBlk.Wall[i][j].yplus)/SolnBlk.W[i][j].von_karman;
// 	      cout << endl << " -> Xi = " << Xi;
// 	      dXi = 0.25/(utau*SolnBlk.W[i][j].von_karman);
// 	      cout << endl << " -> dXi = " << dXi;
// 	      kappav = SolnBlk.W[i][j].von_karman/(ONE + Xi*SolnBlk.Wall[i][j].vwplus);
// 	      cout << endl << " -> kappav = " << kappav;
// 	      dkappav = kappav*SolnBlk.Wall[i][j].vwplus*(Xi/utau - dXi)/(ONE + Xi*SolnBlk.Wall[i][j].vwplus);
// 	      cout << endl << " -> dkappav = " << dkappav;
// 	      E = exp(kappav*SolnBlk.W[i][j].C);
// 	      cout << endl << " -> E = " << E;
// 	      dE = dkappav*SolnBlk.W[i][j].C*E;
// 	      cout << endl << " -> dE = " << dE;
// 	      f = utau - kappav*u/log(E*SolnBlk.Wall[i][j].yplus);
// 	      cout << endl << " -> f = " << f;
// 	      df = ONE + kappav*(u/utau)*(dE/E + ONE)/sqr(log(E*SolnBlk.Wall[i][j].yplus)) - dkappav*u/log(E*SolnBlk.Wall[i][j].yplus);
// 	      cout << endl << " -> df = " << df;
// 	      utau -= f/df;
// 	      cout << endl << " -> utau = " << utau;
// 	      SolnBlk.Wall[i][j].yplus = utau*ywall/nu;
// 	      cout << endl << " -> yplus = " << SolnBlk.Wall[i][j].yplus;
// 	      SolnBlk.Wall[i][j].vwplus = SolnBlk.Wall[i][j].vinj*0.5/utau;
// 	      cout << endl << " -> vwplus = " << SolnBlk.Wall[i][j].vwplus;
// 	      n++;
// 	    } while(fabs(f) >= 0.001 && n <= nmax);
//  	  }
	}
      }
      // Set the final friction velocity.
      SolnBlk.Wall[i][j].utau = utau;
      // Set the final wall shear stress.
      SolnBlk.Wall[i][j].tauw = SolnBlk.W[i][j].rho*sqr(fabs(utau));
      // Determine y-plus for the current cell.
      SolnBlk.Wall[i][j].yplus = utau*ywall/nu;
      }

    }
  }

  // Dimensionless distance to wall calculation complete.
  return 0;

}

/**********************************************************************
 * Routine: Pipe_Dimensionless_Wall_Distance                          *
 *                                                                    *
 * This routine calculates the dimensionless wall distance (y-plus)   *
 * based on the exact friction velocity for a turbulent pipe flow.    *
 *                                                                    *
 **********************************************************************/
template <class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int Pipe_Dimensionless_Wall_Distance(Quad_Soln_Block &SolnBlk,
				     Quad_Soln_Input_Parameters &IP) {

  for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
      // Set the wall shear stress for a turbulent pipe.
      SolnBlk.Wall[i][j].tauw = -HALF*IP.Pipe_Radius*(IP.dp/IP.Pipe_Length);
      // Set the final friction velocity for a turbulent pipe.
      SolnBlk.Wall[i][j].utau = sqrt(SolnBlk.Wall[i][j].tauw/SolnBlk.W[i][j].rho);
      // Determine y-plus for the cell.
      SolnBlk.Wall[i][j].yplus = SolnBlk.Wall[i][j].utau*
                                 SolnBlk.Wall[i][j].ywall/
                                 SolnBlk.W[i][j].nu();
    }
  }

  // Dimensionless distance to wall calculation complete.
  return 0;

}

/**********************************************************************
 * The next three routines manages the calculation and distribution   *
 * of the dimensionless mass injection rate, vwplus, defined at       *
 * burning surface boundary conditions.                               *
 **********************************************************************/

// /**********************************************************************
//  * Routine: Dimensionless_Wall_Injection_Speed                        *
//  *                                                                    *
//  * This routine manages the calculation and distribution of the       *
//  * dimensionless wall injection speed, vwplus, defined at burning     *
//  * surface boundary conditions for a list of solution blocks.         *
//  *                                                                    *
//  **********************************************************************/
// template <class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
// int Dimensionless_Wall_Injection_Speed(Quad_Soln_Block *Soln_ptr,
// 				       AdaptiveBlockResourceList &GlobalSolnBlockList,
// 				       AdaptiveBlock2D_List &LocalSolnBlockList,
// 				       Quad_Soln_Input_Parameters &IP) {

//   // Exit immediately if the wall injection flag is not on.
//   if (!IP.i_Turbulent_Wall_Injection) return 0;

//   // Exit immediately if not a turbulent flow.
//   if (IP.FlowType == FLOWTYPE_INVISCID || IP.FlowType == FLOWTYPE_LAMINAR) return 0;

//   int error_flag, buffer_size;
//   LinkedList<double> LL_buffer;
//   double *buffer;

//   // Initialize buffer size.
//   buffer_size = 0;
//   // Cycle through the processors to create, load, and broadcast the
//   // message buffer containing the injection speed at all burning 
//   // surface boundary cells.  Each processor should then apply the
//   // information accordingly.
//   for (int iCPU = 0; iCPU < GlobalSolnBlockList.Ncpu; iCPU++) {
//     // Compile a linked list of the wall injection speeds on the processor.
//     if (LL_buffer.np) LL_buffer.deallocate();
//     if (LocalSolnBlockList.ThisCPU == iCPU) {
//       for (int nb = 0; nb < LocalSolnBlockList.Nblk; nb++) {
// 	if (LocalSolnBlockList.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
// 	  error_flag = Compile_Wall_Injection_Speed_List(Soln_ptr[nb],
// 							 LL_buffer);
// 	  if (error_flag) return error_flag;
// 	  buffer_size = LL_buffer.np;
// 	}
//       }
//     }
//     // Create the message buffer based on the maximum buffer size.
//     buffer_size = CFDkit_Maximum_MPI(buffer_size);
//     buffer = new double[buffer_size];
//     // Load message buffer.
//     if (LocalSolnBlockList.ThisCPU == iCPU) {
//       for (int n = 0; n < buffer_size; n += 3) {
// 	buffer[n  ] = LL_buffer[n  ];
// 	buffer[n+1] = LL_buffer[n+1];
// 	buffer[n+2] = LL_buffer[n+2];
//       }
//     } else {
//       for (int n = 0; n < buffer_size; n++) buffer[n] = ZERO;
//     }
//     // Broadcast message buffer.
// #ifdef _MPI_VERSION
//     MPI::COMM_WORLD.Bcast(buffer,buffer_size,MPI::DOUBLE,iCPU);
// #endif
//     // Apply wall injection speed information.
//     for (int nb = 0; nb < LocalSolnBlockList.Nblk; nb++) {
//       if (LocalSolnBlockList.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
// 	for (int n = 0; n < buffer_size; n += 3) {
// 	  error_flag = Dimensionless_Wall_Injection_Speed(Soln_ptr[nb],
// 							  Vector2D(buffer[n],buffer[n+1]),
// 							  buffer[n+2]);
// 	  if (error_flag) return error_flag;
// 	}
//       }
//     }
//     // Deallocate message buffer.
//     delete []buffer; buffer = NULL; buffer_size = 0;
//   }

//   // Dimensionless distance to wall calculation complete.
//   return 0;

// }

// /**********************************************************************
//  * Routine: Compile_Wall_Injection_Speed_List                         *
//  *                                                                    *
//  * This routine compiles a list of dimensional wall injection speeds  *
//  * defined at all burning surface boundary conditions of the current  *
//  * block.                                                             *
//  *                                                                    *
//  **********************************************************************/
// template <class Quad_Soln_Block>
// int Compile_Wall_Injection_Speed_List(Quad_Soln_Block &SolnBlk,
// 				      LinkedList<double> &LL_buffer) {

//   Vector2D vbs;

//   for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
//     // North block boundary.
//     if (SolnBlk.Grid.BCtypeN[i] == BC_BURNING_SURFACE) {
//       LL_buffer.add(SolnBlk.Grid.xfaceN(i,SolnBlk.JCu).x);
//       LL_buffer.add(SolnBlk.Grid.xfaceN(i,SolnBlk.JCu).y);
//       vbs = BurningSurface(SolnBlk.W[i][SolnBlk.JCu],SolnBlk.Grid.nfaceN(i,SolnBlk.JCu)).v;
//       LL_buffer.add(dot(vbs,-SolnBlk.Grid.nfaceN(i,SolnBlk.JCu)));
//     }
//     // South block boundary.
//     if (SolnBlk.Grid.BCtypeS[i] == BC_BURNING_SURFACE) {
//       LL_buffer.add(SolnBlk.Grid.xfaceS(i,SolnBlk.JCl).x);
//       LL_buffer.add(SolnBlk.Grid.xfaceS(i,SolnBlk.JCl).y);
//       vbs = BurningSurface(SolnBlk.W[i][SolnBlk.JCl],SolnBlk.Grid.nfaceS(i,SolnBlk.JCl)).v;
//       LL_buffer.add(dot(vbs,-SolnBlk.Grid.nfaceS(i,SolnBlk.JCl)));
//     }
//   }

//   for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
//     // East block boundary.
//     if (SolnBlk.Grid.BCtypeE[j] == BC_BURNING_SURFACE) {
//       LL_buffer.add(SolnBlk.Grid.xfaceE(SolnBlk.ICu,j).x);
//       LL_buffer.add(SolnBlk.Grid.xfaceE(SolnBlk.ICu,j).y);
//       vbs = BurningSurface(SolnBlk.W[SolnBlk.ICu][j],SolnBlk.Grid.nfaceE(SolnBlk.ICu,j)).v;
//       LL_buffer.add(dot(vbs,-SolnBlk.Grid.nfaceE(SolnBlk.ICu,j)));
//     }
//     // West block boundary.
//     if (SolnBlk.Grid.BCtypeW[j] == BC_BURNING_SURFACE) {
//       LL_buffer.add(SolnBlk.Grid.xfaceW(SolnBlk.ICl,j).x);
//       LL_buffer.add(SolnBlk.Grid.xfaceW(SolnBlk.ICl,j).y);
//       vbs = BurningSurface(SolnBlk.W[SolnBlk.ICl][j],SolnBlk.Grid.nfaceW(SolnBlk.ICl,j)).v;
//       LL_buffer.add(dot(vbs,-SolnBlk.Grid.nfaceN(SolnBlk.ICl,j)));
//     }
//   }
// //   for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
// //     // North block boundary.
// //     if (SolnBlk.Grid.BCtypeN[i] == BC_BURNING_SURFACE) {
// //       LL_buffer.add(SolnBlk.Grid.xfaceN(i,SolnBlk.JCu).x);
// //       LL_buffer.add(SolnBlk.Grid.xfaceN(i,SolnBlk.JCu).y);
// //       LL_buffer.add(SolnBlk.W[i][SolnBlk.JCu].burningrate());
// //     }
// //     // South block boundary.
// //     if (SolnBlk.Grid.BCtypeS[i] == BC_BURNING_SURFACE) {
// //       LL_buffer.add(SolnBlk.Grid.xfaceS(i,SolnBlk.JCl).x);
// //       LL_buffer.add(SolnBlk.Grid.xfaceS(i,SolnBlk.JCl).y);
// //       LL_buffer.add(SolnBlk.W[i][SolnBlk.JCl].burningrate());
// //     }
// //   }

// //   for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
// //     // East block boundary.
// //     if (SolnBlk.Grid.BCtypeE[j] == BC_BURNING_SURFACE) {
// //       LL_buffer.add(SolnBlk.Grid.xfaceE(SolnBlk.ICu,j).x);
// //       LL_buffer.add(SolnBlk.Grid.xfaceE(SolnBlk.ICu,j).y);
// //       LL_buffer.add(SolnBlk.W[SolnBlk.ICu][j].burningrate());
// //     }
// //     // West block boundary.
// //     if (SolnBlk.Grid.BCtypeW[j] == BC_BURNING_SURFACE) {
// //       LL_buffer.add(SolnBlk.Grid.xfaceW(SolnBlk.ICl,j).x);
// //       LL_buffer.add(SolnBlk.Grid.xfaceW(SolnBlk.ICl,j).y);
// //       LL_buffer.add(SolnBlk.W[SolnBlk.ICl][j].burningrate());
// //     }
// //   }

//   // List of wall injection speeds compiled successfully.
//   return 0;

// }

// /**********************************************************************
//  * Routine: Dimensionless_Wall_Injection_Speed                        *
//  *                                                                    *
//  * This routine determines the dimensionless wall injection speed for *
//  * each cell of the current block if the associated cell is the       *
//  * correct cell for the calculation.                                  *
//  *                                                                    *
//  **********************************************************************/
// template <class Quad_Soln_Block>
// int Dimensionless_Wall_Injection_Speed(Quad_Soln_Block &SolnBlk,
// 				       const Vector2D &Xwall,
// 				       const double &vwall) {

//   for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
//     for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
//       if (abs(Xwall - SolnBlk.Wall[i][j].Xwall) < NANO) {
// 	if (SolnBlk.Wall[i][j].utau > ZERO) {
// 	  SolnBlk.Wall[i][j].vinj = vwall;
// 	  SolnBlk.Wall[i][j].vwplus = SolnBlk.Wall[i][j].vinj/SolnBlk.Wall[i][j].utau;
// 	} else {
// 	  SolnBlk.Wall[i][j].vwplus = ZERO;
// 	}
//       }
//     }
//   }

//   // Dimensionless distance to wall calculation complete.
//   return 0;

// }
/**********************************************************************
 * Routine: Wall_Injection_Speed                                      *
 *                                                                    *
 * This routine manages the calculation and distribution of the       *
 * wall injection speed, vinj, defined at burning surface boundary    *
 * conditions for a list of solution blocks.                          *
 *                                                                    *
 **********************************************************************/
template <class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int Wall_Injection_Speed(Quad_Soln_Block *Soln_ptr,
			 AdaptiveBlockResourceList &GlobalSolnBlockList,
			 AdaptiveBlock2D_List &LocalSolnBlockList,
			 Quad_Soln_Input_Parameters &IP) {

  // Exit immediately if the wall injection flag is not on.
  if (!IP.i_Turbulent_Wall_Injection) return 0;

  // Exit immediately if not a turbulent flow.
  if (IP.FlowType == FLOWTYPE_INVISCID || IP.FlowType == FLOWTYPE_LAMINAR) return 0;

  int error_flag, buffer_size;
  LinkedList<double> LL_buffer;
  double *buffer;

  // Initialize buffer size.
  buffer_size = 0;
  // Cycle through the processors to create, load, and broadcast the
  // message buffer containing the injection speed at all burning 
  // surface boundary cells.  Each processor should then apply the
  // information accordingly.
  for (int iCPU = 0; iCPU < GlobalSolnBlockList.Ncpu; iCPU++) {
    // Compile a linked list of the wall injection speeds on the processor.
    if (LL_buffer.np) LL_buffer.deallocate();
    if (LocalSolnBlockList.ThisCPU == iCPU) {
      for (int nb = 0; nb < LocalSolnBlockList.Nblk; nb++) {
	if (LocalSolnBlockList.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
	  error_flag = Compile_Wall_Injection_Speed_List(Soln_ptr[nb],
							 LL_buffer);
	  if (error_flag) return error_flag;
	  buffer_size = LL_buffer.np;
	}
      }
    }
    // Create the message buffer based on the maximum buffer size.
    buffer_size = CFDkit_Maximum_MPI(buffer_size);
    buffer = new double[buffer_size];
    // Load message buffer.
    if (LocalSolnBlockList.ThisCPU == iCPU) {
      for (int n = 0; n < buffer_size; n += 3) {
	buffer[n  ] = LL_buffer[n  ];
	buffer[n+1] = LL_buffer[n+1];
	buffer[n+2] = LL_buffer[n+2];
      }
    } else {
      for (int n = 0; n < buffer_size; n++) buffer[n] = ZERO;
    }
    // Broadcast message buffer.
#ifdef _MPI_VERSION
    MPI::COMM_WORLD.Bcast(buffer,buffer_size,MPI::DOUBLE,iCPU);
#endif
    // Apply wall injection speed information.
    for (int nb = 0; nb < LocalSolnBlockList.Nblk; nb++) {
      if (LocalSolnBlockList.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
	for (int n = 0; n < buffer_size; n += 3) {
	  error_flag = Wall_Injection_Speed(Soln_ptr[nb],
					    Vector2D(buffer[n],buffer[n+1]),
					    buffer[n+2]);
	  if (error_flag) return error_flag;
	}
      }
    }
    // Deallocate message buffer.
    delete []buffer; buffer = NULL; buffer_size = 0;
  }

  // Wall injection speed calculations complete.
  return 0;

}

/**********************************************************************
 * Routine: Compile_Wall_Injection_Speed_List                         *
 *                                                                    *
 * This routine compiles a list of dimensional wall injection speeds  *
 * defined at all burning surface boundary conditions of the current  *
 * block.                                                             *
 *                                                                    *
 **********************************************************************/
template <class Quad_Soln_Block>
int Compile_Wall_Injection_Speed_List(Quad_Soln_Block &SolnBlk,
				      LinkedList<double> &LL_buffer) {

  Vector2D vbs;

  for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
    // North block boundary.
    if (SolnBlk.Grid.BCtypeN[i] == BC_BURNING_SURFACE) {
      LL_buffer.add(SolnBlk.Grid.xfaceN(i,SolnBlk.JCu).x);
      LL_buffer.add(SolnBlk.Grid.xfaceN(i,SolnBlk.JCu).y);
      vbs = BurningSurface(SolnBlk.W[i][SolnBlk.JCu],SolnBlk.Grid.nfaceN(i,SolnBlk.JCu)).v;
      LL_buffer.add(dot(vbs,-SolnBlk.Grid.nfaceN(i,SolnBlk.JCu)));
    }
    // South block boundary.
    if (SolnBlk.Grid.BCtypeS[i] == BC_BURNING_SURFACE) {
      LL_buffer.add(SolnBlk.Grid.xfaceS(i,SolnBlk.JCl).x);
      LL_buffer.add(SolnBlk.Grid.xfaceS(i,SolnBlk.JCl).y);
      vbs = BurningSurface(SolnBlk.W[i][SolnBlk.JCl],SolnBlk.Grid.nfaceS(i,SolnBlk.JCl)).v;
      LL_buffer.add(dot(vbs,-SolnBlk.Grid.nfaceS(i,SolnBlk.JCl)));
    }
  }

  for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
    // East block boundary.
    if (SolnBlk.Grid.BCtypeE[j] == BC_BURNING_SURFACE) {
      LL_buffer.add(SolnBlk.Grid.xfaceE(SolnBlk.ICu,j).x);
      LL_buffer.add(SolnBlk.Grid.xfaceE(SolnBlk.ICu,j).y);
      vbs = BurningSurface(SolnBlk.W[SolnBlk.ICu][j],SolnBlk.Grid.nfaceE(SolnBlk.ICu,j)).v;
      LL_buffer.add(dot(vbs,-SolnBlk.Grid.nfaceE(SolnBlk.ICu,j)));
    }
    // West block boundary.
    if (SolnBlk.Grid.BCtypeW[j] == BC_BURNING_SURFACE) {
      LL_buffer.add(SolnBlk.Grid.xfaceW(SolnBlk.ICl,j).x);
      LL_buffer.add(SolnBlk.Grid.xfaceW(SolnBlk.ICl,j).y);
      vbs = BurningSurface(SolnBlk.W[SolnBlk.ICl][j],SolnBlk.Grid.nfaceW(SolnBlk.ICl,j)).v;
      LL_buffer.add(dot(vbs,-SolnBlk.Grid.nfaceN(SolnBlk.ICl,j)));
    }
  }
//   for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
//     // North block boundary.
//     if (SolnBlk.Grid.BCtypeN[i] == BC_BURNING_SURFACE) {
//       LL_buffer.add(SolnBlk.Grid.xfaceN(i,SolnBlk.JCu).x);
//       LL_buffer.add(SolnBlk.Grid.xfaceN(i,SolnBlk.JCu).y);
//       LL_buffer.add(SolnBlk.W[i][SolnBlk.JCu].burningrate());
//     }
//     // South block boundary.
//     if (SolnBlk.Grid.BCtypeS[i] == BC_BURNING_SURFACE) {
//       LL_buffer.add(SolnBlk.Grid.xfaceS(i,SolnBlk.JCl).x);
//       LL_buffer.add(SolnBlk.Grid.xfaceS(i,SolnBlk.JCl).y);
//       LL_buffer.add(SolnBlk.W[i][SolnBlk.JCl].burningrate());
//     }
//   }

//   for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
//     // East block boundary.
//     if (SolnBlk.Grid.BCtypeE[j] == BC_BURNING_SURFACE) {
//       LL_buffer.add(SolnBlk.Grid.xfaceE(SolnBlk.ICu,j).x);
//       LL_buffer.add(SolnBlk.Grid.xfaceE(SolnBlk.ICu,j).y);
//       LL_buffer.add(SolnBlk.W[SolnBlk.ICu][j].burningrate());
//     }
//     // West block boundary.
//     if (SolnBlk.Grid.BCtypeW[j] == BC_BURNING_SURFACE) {
//       LL_buffer.add(SolnBlk.Grid.xfaceW(SolnBlk.ICl,j).x);
//       LL_buffer.add(SolnBlk.Grid.xfaceW(SolnBlk.ICl,j).y);
//       LL_buffer.add(SolnBlk.W[SolnBlk.ICl][j].burningrate());
//     }
//   }

  // List of wall injection speeds compiled successfully.
  return 0;

}

/**********************************************************************
 * Routine: Wall_Injection_Speed                                      *
 *                                                                    *
 * This routine determines the wall injection speed for each cell of  *
 * the current block if the associated cell is the correct cell.      *
 *                                                                    *
 **********************************************************************/
template <class Quad_Soln_Block>
int Wall_Injection_Speed(Quad_Soln_Block &SolnBlk,
			 const Vector2D &Xwall,
			 const double &vwall) {

  for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
    for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
      if (abs(Xwall - SolnBlk.Wall[i][j].Xwall) < NANO) {
	SolnBlk.Wall[i][j].vinj = vwall;
      }
    }
  }

  // Dimensionless distance to wall calculation complete.
  return 0;

}

/**********************************************************************
 * Routine: Dimensionless_Wall_Injection_Speed                        *
 *                                                                    *
 * This routine manages the calculation and distribution of the       *
 * dimensionless wall injection speed, vwplus, defined at burning     *
 * surface boundary conditions for a list of solution blocks.         *
 *                                                                    *
 **********************************************************************/
template <class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int Dimensionless_Wall_Injection_Speed(Quad_Soln_Block *Soln_ptr,
				       AdaptiveBlock2D_List &LocalSolnBlockList,
				       Quad_Soln_Input_Parameters &IP) {

  // Exit immediately if the wall injection flag is not on.
  if (!IP.i_Turbulent_Wall_Injection) return 0;

  // Exit immediately if not a turbulent flow.
  if (IP.FlowType == FLOWTYPE_INVISCID || IP.FlowType == FLOWTYPE_LAMINAR) return 0;

  // Apply wall injection speed information.
  for (int nb = 0; nb < LocalSolnBlockList.Nblk; nb++) {
    if (LocalSolnBlockList.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      for (int j = Soln_ptr[nb].JCl-Soln_ptr[nb].Nghost; j <= Soln_ptr[nb].JCu+Soln_ptr[nb].Nghost; j++) {
	for (int i = Soln_ptr[nb].ICl-Soln_ptr[nb].Nghost; i <= Soln_ptr[nb].ICu+Soln_ptr[nb].Nghost; i++) {
	  if (Soln_ptr[nb].Wall[i][j].utau > ZERO) {
	    Soln_ptr[nb].Wall[i][j].vwplus = Soln_ptr[nb].Wall[i][j].vinj/Soln_ptr[nb].Wall[i][j].utau;
	  } else {
	    Soln_ptr[nb].Wall[i][j].vwplus = ZERO;
	  }
	}
      }
    }
  }

  // Dimensionless wall injection speed complete.
  return 0;

}

#endif // _TURBULENT2D_WALLDATA_INCLUDED

/**********************************************************************
 * EmbeddedBoundaries2D: Header file declaring 2D embedded boundary   *
 *                       classes and functions.                       *
 **********************************************************************/

#ifndef _EMBEDDEDBOUNDARIES2D_INCLUDED
#define _EMBEDDEDBOUNDARIES2D_INCLUDED

// Include required C++ libraries.

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <cstring>

using namespace std;

// Include the interface header file.

#ifndef _INTERFACE2D_INCLUDED
#include "Interface2D.h"
#endif // _INTERFACE2D_INCLUDED

// Include Embedded Boundaries input header file.

#ifndef _EMBEDDEDBOUNDARIES2DINPUT_INCLUDED
#include "EmbeddedBoundaries2D_Input.h"
#endif // _EMBEDDEDBOUNDARIES2DINPUT_INCLUDED

// Include 2D Level Set quadrilateral mesh solution header file.

#ifndef _LEVELSET2D_QUAD_INCLUDED
#include "../LevelSet2D/LevelSet2DQuad.h"
#endif // _LEVELSET2D_QUAD_INCLUDED

// Include the 2D quadrilateral multiblock grid header file.

#ifndef _GRID2D_QUAD_BLOCK_INCLUDED
#include "../Grid/Grid2DQuad.h"
#endif // _GRID2D_QUAD_BLOCK_INCLUDED

// Include the NASA rotor 37 and 67 header files.

#ifndef _NASA_ROTOR37_INCLUDED
#include "../Grid/NASARotor37.h"
#endif // _NASA_ROTOR37_INCLUDED

#ifndef _NASA_ROTOR67_INCLUDED
#include "../Grid/NASARotor67.h"
#endif // _NASA_ROTOR67_INCLUDED

// Include quadtree header file.

#ifndef _QUADTREE_INCLUDED
#include "../AMR/QuadTree.h"
#endif // _QUADTREE_INCLUDED

// Include AMR header file.

#ifndef _AMR_INCLUDED
#include "../AMR/AMR.h"
#endif // _AMR_INCLUDED

// Include the linear systems header file.

#ifndef _LINEARSYSTEMS_INCLUDED
#include "../Math/LinearSystems.h"
#endif // _LINEARSYSTEMS_INCLUDED

#ifndef _VECTOR2D_INCLUDED
#include "../Math/Vector2D.h"
#endif // _VECTOR2D_INCLUDED

// Include the System header file.

#ifndef _SYSTEM_LINUX_INCLUDED
#include "../System/System_Linux.h"
#endif // _SYSTEM_LINUX_INCLUDED

// Define various mesh adjustment constants.

#define CELL_TYPE_QUADRILATERAL            0
#define CELL_TYPE_ADJUSTED_QUADRILATERAL   1
#define CELL_TYPE_TRIANGLE                 2

#define CELL_STATUS_UNKNOWN    -1
#define CELL_STATUS_ACTIVE      0
#define CELL_STATUS_INACTIVE    1

#define NODE_STATUS_UNKNOWN    -1
#define NODE_STATUS_KNOWN       0
#define NODE_STATUS_ALIGNED     1

// (Parallel) Debugging option.

//#define _EB_PARALLEL_DEBUG_

/*!
 * class: Adjusted_Mesh_Quad_Block
 *
 * @brief Extra solution data required to store the adjusted mesh.
 *
 * Extra solution data required to store the adjusted mesh.  The status
 * of each cell (active, inactive, or unknown) and each node (aligned,
 * known, or unknown) are required during the mesh adjustment algorithm.
 * 'Aligned' nodes are nodes that are aligned with one of the embeddded
 * boundaries.  'Known' nodes are not candidates for adjustment where as
 * the 'unknown' nodes are candidates for adjustment.  All 'active'
 * cells (external to all embedded boundaries) are tagged with a zero
 * value.  All inactive cells (internal to one of the embedded
 * boundaries) are tagged with the number corresponding to the the
 * associated embedded boundary in the list of interface unions.  This 
 * allows for quick determination of the velocity, temperture, and
 * boundary condition of the embedded boundary during the flow solution.
 * The cell types are stored to facilitate the operators used during the
 * restriction and prolongation procedures for the AMR and FAS multigrid.
 *
 * \begin{verbatim}
 * Member variables:
 *     NCi         -- Number of i-direction cells (including ghost cells).
 *     NCj         -- Number of j-direction cells (including ghost cells).
 *     NNi         -- Number of i-direction nodes (including ghost cells).
 *     NNj         -- Number of j-direction nodes (including ghost cells).
 *     cell_type   -- Cell type.
 *     cell_status -- Cell status.
 *     node_status -- Node status flag.
 * Member functions:
 *     allocate    -- Allocate memory for the adjusted mesh quad block.
 *     deallocate  -- Deallocate memory for the adjusted mesh quad block.
 *     Copy_Quad_Block -- Copy the given adjusted mesh quad block.
 *     Broadcast_Quad_Block -- Broadcast the adjusted mesh quad block
 *                             from the specified CPU.
 * Member operators:
 *     <<          -- Output stream operator.
 *     >>          -- Input stream operator.
 * \end{verbatim}
 */
class Adjusted_Mesh_Quad_Block{
 public:
  int           NCi, //!< Number of i-direction cells (including ghost cells).
                NCj, //!< Number of j-direction cells (including ghost cells).
                NNi, //!< Number of i-direction nodes (including ghost cells).
                NNj, //!< Number of j-direction nodes (including ghost cells).
        **cell_type, //!< Cell type.
      **cell_status, //!< Cell status.
      **node_status; //!< Node status flag.

  //! Default constructor.
  Adjusted_Mesh_Quad_Block(void) {
    NCi = 0; NCj = 0; NNi = 0; NNj = 0;
    cell_type = NULL; cell_status = NULL; node_status = NULL;
  }

  //! Constructor.
  Adjusted_Mesh_Quad_Block(const int &Ni, const int &Nj) {
    allocate(Ni,Nj);
  }

  //! Destructor.
  ~Adjusted_Mesh_Quad_Block(void) {
    deallocate();
  }

  //! Allocate memory for the adjusted mesh quad block.
  void allocate(const int &Ni, const int &Nj) {
    assert(Ni > 1 && Nj > 1);
    NCi = Ni; NCj = Nj; NNi = Ni+1; NNj = Nj+1;
    node_status = new int*[NNi];
    for (int i = 0; i < NNi; i++) {
      node_status[i] = new int[NNj];
      for (int j = 0; j < NNj; j++) {
	node_status[i][j] = NODE_STATUS_KNOWN;
      }
    }
    cell_type = new int*[NCi];
    cell_status = new int*[NCi];
    for (int i = 0; i < NCi; i++) {
      cell_type[i] = new int[NCj];
      cell_status[i] = new int[NCj];
      for (int j = 0; j < NCj; j++) {
	cell_type[i][j] = CELL_TYPE_QUADRILATERAL;
	cell_status[i][j] = CELL_STATUS_ACTIVE;
      }
    }
  }

  //! Deallocate memory for the adjusted mesh quad block.
  void deallocate(void) {
    if (node_status != NULL) {
      for (int i = 0; i < NNi; i++) {
	delete []node_status[i]; node_status[i] = NULL;
      }
      delete []node_status; node_status = NULL;
    }
    if (cell_status != NULL) {
      for (int i = 0; i < NCi; i++) {
	delete []cell_type[i]; cell_type[i] = NULL;
	delete []cell_status[i]; cell_status[i] = NULL;
      }
      delete []cell_type; cell_type = NULL;
      delete []cell_status; cell_status = NULL;
    }
    NCi = 0; NCj = 0; NNi = 0; NNj = 0;
  }

  //! Copy the given adjusted mesh quad block.
  void Copy_Quad_Block(const Adjusted_Mesh_Quad_Block &M) {
    if (NCi != M.NCi || NCj != M.NCj) {
      if (node_status != NULL) deallocate();
      if (M.node_status != NULL) allocate(M.NCi,M.NCj);
    }
    for (int j = 0; j < NNj; j++) {
      for (int i = 0; i < NNi; i++) {
	node_status[i][j] = M.node_status[i][j];
      }
    }
    for (int j = 0; j < NCj; j++) {
      for (int i = 0; i < NCi; i++) {
	cell_type[i][j] = M.cell_type[i][j];
	cell_status[i][j] = M.cell_status[i][j];
      }
    }
  }

#ifdef _MPI_VERSION
  //! Broadcast the adjusted mesg quad block from the specified CPU.
  void Broadcast_Quad_Block(MPI::Intracomm &Communicator,
			    const int Source_CPU) {
    int Source_Rank = 0;
    int ni, nj, block_allocated, buffer_size, *buffer;

    // Broadcast the number of cells in each direction.
    if (CFFC_MPI::This_Processor_Number == Source_CPU) {
      ni = NCi;
      nj = NCj;
      if (node_status != NULL) block_allocated = 1;
      else block_allocated = 0;
    }

    Communicator.Bcast(&ni,1,MPI::INT,Source_Rank);
    Communicator.Bcast(&nj,1,MPI::INT,Source_Rank);
    Communicator.Bcast(&block_allocated,1,MPI::INT,Source_Rank);

    // On non-source MPI processors, allocate (re-allocate) memory for
    // the quadrilateral block as necessary.
    if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
      if (NCi != ni || NCj != nj) { 
	if (node_status != NULL) deallocate();
	if (block_allocated) allocate(ni,nj); 
      }
    }

    if (block_allocated) {

      // Broadcast the node status array.
      buffer = new int[NNi*NNj];

      if (CFFC_MPI::This_Processor_Number == Source_CPU) {
	buffer_size = 0;
	for (int j = 0; j < NNj; j++) {
	  for (int i = 0; i < NNi; i++) {
	    buffer[buffer_size] = node_status[i][j];
	    buffer_size++;
	  }
	}
      }

      buffer_size = NNi*NNj;
      Communicator.Bcast(buffer,buffer_size,MPI::INT,Source_Rank);

      if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
	buffer_size = 0;
	for (int j = 0; j < NNj; j++) {
	  for (int i = 0; i < NNi; i++) {
	    node_status[i][j] = buffer[buffer_size];
	    buffer_size++;
	  }
	}
      }

      delete []buffer; buffer = NULL;

      // Broadcast the cell status and type arrays.
      buffer = new int[2*NCi*NCj];

      if (CFFC_MPI::This_Processor_Number == Source_CPU) {
	buffer_size = 0;
	for (int j = 0; j < NCj; j++) {
	  for (int i = 0; i < NCi; i++) {
	    buffer[buffer_size  ] = cell_status[i][j];
	    buffer[buffer_size+1] = cell_type[i][j];
	    buffer_size += 2;
	  }
	}
      }

      buffer_size = 2*NCi*NCj;
      Communicator.Bcast(buffer,buffer_size,MPI::INT,Source_Rank);

      if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
	buffer_size = 0;
	for (int j = 0; j < NCj; j++) {
	  for (int i = 0; i < NCi; i++) {
	    cell_status[i][j] = buffer[buffer_size  ];
	    cell_type[i][j]   = buffer[buffer_size+1];
	    buffer_size += 2;
	  }
	}
      }

      delete []buffer; buffer = NULL;

    }

  }
#endif

  //@{ @name Input-output operators.
  friend ostream &operator << (ostream &out_file, const Adjusted_Mesh_Quad_Block &M) {
    out_file << M.NCi << " " << M.NCj << endl;
    if (M.NCi == 0 || M.NCj == 0) return out_file;
    out_file << M.NNi << " " << M.NNj << endl;
    for (int j = 0; j < M.NNj; j++) {
      for (int i = 0; i < M.NNi; i++) {
	out_file << M.node_status[i][j] << endl;
      }
    }
    for (int j = 0; j < M.NCj; j++) {
      for (int i = 0; i < M.NCi; i++) {
	out_file << M.cell_status[i][j] << " " << M.cell_type[i][j] << endl;
      }
    }
    return out_file;
  }

  friend istream &operator >> (istream &in_file, Adjusted_Mesh_Quad_Block &M) {
    int nci, ncj, nni, nnj;
    in_file.setf(ios::skipws);
    in_file >> nci >> ncj;
    in_file.unsetf(ios::skipws);
    if (nci == 0 || ncj == 0) { if (M.node_status != NULL) M.deallocate(); return in_file; }
    in_file.setf(ios::skipws);
    in_file >> nni >> nnj;
    in_file.unsetf(ios::skipws);
    if (M.node_status == NULL || M.NNi != nni || M.NNj != nnj) {
      if (M.node_status != NULL) M.deallocate(); 
      M.allocate(nci,ncj);
    }
    for (int j = 0; j < M.NNj; j++) {
      for (int i = 0; i < M.NNi; i++) {
	in_file.setf(ios::skipws);
	in_file >> M.node_status[i][j];
	in_file.unsetf(ios::skipws);
      }
    }
    for (int j = 0; j < M.NCj; j++) {
      for (int i = 0; i < M.NCi; i++) {
	in_file.setf(ios::skipws);
	in_file >> M.cell_status[i][j] >> M.cell_type[i][j];
	in_file.unsetf(ios::skipws);
      }
    }
    return in_file;
  }
  //@}

};

/*!
 * Class: Mesh_Adjustment_Data
 *
 * @brief Class definition of variables required during the mesh
 *        adjustment algorithm.
 *
 * Class definition of the variables required during the mesh adjustment
 * algorithm.  Intersection flags for all four directions are required 
 * for every node of the solution domain and indicate whether or not an
 * intersection point occurs between the mesh line in that direction and
 * one of the embedded boundaries.  The pierce points are also found and
 * are stored.  The number of inactive cells in the solution block (not
 * including ghost cells) and a flag indicating whether or not an
 * embedded boundary is present in the solution block (including ghost
 * cells) are also stored.
 *
 * \verbatim
 * Member variables:
 *   NCi             -- Total number of i-direction cells.
 *   ICl             -- First i-direction non-ghost cell counter.
 *   ICu             -- Final i-direction non-ghost cell counter.
 *   NCj             -- Total number of j-direction cells.
 *   JCl             -- First j-direction non-ghost cell counter.
 *   JCu             -- Final j-direction non-ghost cell counter.
 *   NNi             -- Number of i-direction nodes (including ghost cells).
 *   NNj             -- Number of j-direction nodes (including ghost cells).
 *   Nghost          -- Number of ghost cells.
 *   intersect_north -- North direction intersection flag.
 *   intersect_south -- South direction intersection flag.
 *   intersect_east  -- East direction intersection flag.
 *   intersect_west  -- West direction intersection flag.
 *   Xnorth          -- North direction intersection point.
 *   Xsouth          -- South direction intersection point.
 *   Xeast           -- East direction intersection point.
 *   Xwest           -- West direction intersection point.
 *   Interface_Present -- Interface present in block domain flag.
 *   Number_of_Inactive_Cells -- Number of inactive computational cells.
 * Member functions:
 *   allocate        -- Allocate memory for the mesh adjustment data block.
 *   dellocate       -- Dellocate memory for the mesh adjustment data block.
 *   Copy_Quad_Block -- Copy the specified mesh adjustment data block.
 *   Initialize      -- Initialize the mesh adjustment data block.
 *   Number_of_Active_Cells -- Return the number of active cells.
 * \endverbatim
 */
class Mesh_Adjustment_Data{
private:
public:
  int                      NCi, //!< Total number of i-direction cells.
                           ICl, //!< First i-direction non-ghost cell counter.
                           ICu; //!< Final i-direction non-ghost cell counter.
  int                      NCj, //!< Total number of j-direction cells.
                           JCl, //!< First j-direction non-ghost cell counter.
                           JCu; //!< Final j-direction non-ghost cell counter.
  int                      NNi, //!< Number of i-direction nodes (including ghost cells).
                           NNj; //!< Number of j-direction nodes (including ghost cells).
  int                   Nghost; //!< Number of ghost cells.
  int        **intersect_north, //!< North direction intersection flag.
             **intersect_south, //!< South direction intersection flag.
              **intersect_east, //!< East direction intersection flag.
              **intersect_west; //!< West direction intersection flag.
  Vector2D            **Xnorth, //!< North direction intersection point.
                      **Xsouth, //!< South direction intersection point.
                       **Xeast, //!< East direction intersection point.
                       **Xwest; //!< West direction intersection point.
  int Number_of_Inactive_Cells; //!< Number of inactive computational cells.
  int                       Ni; //!< Number of union interfaces.
  int       *Interface_Present; //!< Interface present in block domain flags.

  //! Default constructor.
  Mesh_Adjustment_Data(void) {
    NCi = 0; ICl = 0; ICu = 0; NCj = 0; JCl = 0; JCu = 0; Nghost = 0;
    NNi = 0; NNj = 0;
    intersect_north = NULL; intersect_south = NULL;
    intersect_east = NULL; intersect_west = NULL;
    Xnorth = NULL; Xsouth = NULL;
    Xeast = NULL;  Xwest = NULL;
    Number_of_Inactive_Cells = 0;
    Ni = 0;
    Interface_Present = NULL;
  }

  //! Creation constructor.
  Mesh_Adjustment_Data(const int &nni, const int &nnj, const int &Ng, const int &Nint) {
    allocate(nni,nnj,Ng,Nint);
  }

  //! Default destructor.
  ~Mesh_Adjustment_Data(void) { deallocate(); }

  //! Allocate memory for the mesh adjustment data block.
  void allocate(const int &nni, const int &nnj, const int &Ng, const int &Nint) {
    assert(nni > 1 && nnj > 1);
    NCi = nni; ICl = Ng; ICu = nni-2*Ng+1; Nghost = Ng;
    NCj = nnj; JCl = Ng; JCu = nnj-2*Ng+1;
    NNi = nni+1; NNj = nnj+1;
    intersect_north = new int*[NNi]; intersect_south = new int*[NNi];
    intersect_east = new int*[NNi];  intersect_west = new int*[NNi];
    Xnorth = new Vector2D*[NNi]; Xsouth = new Vector2D*[NNi];
    Xeast = new Vector2D*[NNi];  Xwest = new Vector2D*[NNi];
    for (int i = 0; i < NNi; i++) {
      intersect_north[i] = new int[NNj]; intersect_south[i] = new int[NNj];
      intersect_east[i] = new int[NNj];  intersect_west[i] = new int[NNj];
      Xnorth[i] = new Vector2D[NNj]; Xsouth[i] = new Vector2D[NNj];
      Xeast[i] = new Vector2D[NNj];  Xwest[i] = new Vector2D[NNj];
      for (int j = 0; j < NNj; j++) {
	intersect_north[i][j] = OFF; intersect_south[i][j] = OFF;
	intersect_east[i][j] = OFF;  intersect_west[i][j] = OFF;
	Xnorth[i][j] = Vector2D_ZERO; Xsouth[i][j] = Vector2D_ZERO;
	Xeast[i][j] = Vector2D_ZERO;  Xwest[i][j] = Vector2D_ZERO;
      }
    }
    Number_of_Inactive_Cells = 0;
    Ni = Nint;
    Interface_Present = new int[Ni+1];
    for (int ni = 0; ni < Ni+1; ni++) Interface_Present[ni] = OFF;
  }

  //! Deallocate memory for the mesh adjustment data block.
  void deallocate(void) {
    if (intersect_north != NULL) {
      for (int i = 0; i < NNi; i++) {
	delete []intersect_north[i]; intersect_north[i] = NULL;
	delete []intersect_south[i]; intersect_south[i] = NULL;
	delete []intersect_east[i]; intersect_east[i] = NULL;
	delete []intersect_west[i]; intersect_west[i] = NULL;
	delete []Xnorth[i]; Xnorth[i] = NULL;
	delete []Xsouth[i]; Xsouth[i] = NULL;
	delete []Xeast[i]; Xeast[i] = NULL;
	delete []Xwest[i]; Xwest[i] = NULL;
      }
      delete []intersect_north; intersect_north = NULL;
      delete []intersect_south; intersect_south = NULL;
      delete []intersect_east; intersect_east = NULL;
      delete []intersect_west; intersect_west = NULL;
      delete []Xnorth; Xnorth = NULL;
      delete []Xsouth; Xsouth = NULL;
      delete []Xeast; Xeast = NULL;
      delete []Xwest; Xwest = NULL;
    }
    NCi = 0; ICl = 0; ICu = 0; NCj = 0; JCl = 0; JCu = 0; Nghost = 0;
    NNi = 0; NNj = 0;
    if (Interface_Present != NULL) { Ni = 0; delete []Interface_Present; Interface_Present = NULL; }
    Number_of_Inactive_Cells = 0;
  }

  //! Copy the specified mesh adjustment data block.
  void Copy_Quad_Block(const Mesh_Adjustment_Data &M) {
    NCi = M.NCi; ICl = M.ICl; ICu = M.ICu;
    NCj = M.NCj; JCl = M.JCl; JCu = M.JCu;
    Nghost = M.Nghost;
    NNi = M.NNi; NNj = M.NNj;
    Number_of_Inactive_Cells = M.Number_of_Inactive_Cells;
    for (int j = 0; j < NNj; j++) {
      for (int i = 0; i < NNi; i++) {
	intersect_north[i][j] = M.intersect_north[i][j];
	intersect_south[i][j] = M.intersect_south[i][j];
	intersect_east[i][j] = M.intersect_east[i][j];
	intersect_west[i][j] = M.intersect_west[i][j];
	Xnorth[i][j] = M.Xnorth[i][j];
	Xsouth[i][j] = M.Xsouth[i][j];
	Xeast[i][j] = M.Xeast[i][j];
	Xwest[i][j] = M.Xwest[i][j];
      }
    }
    if (Interface_Present != NULL) { Ni = 0; delete []Interface_Present; Interface_Present = NULL; }
    Ni = M.Ni;
    if (Ni) {
      Interface_Present = new int[Ni+1];
      for (int ni = 0; ni < Ni+1; ni++) Interface_Present[ni] = OFF;
    }
  }

  //! Initialize the mesh adjustment data block.
  void Initialize(void) {
    for (int i = 0; i < NNi; i++) {
      for (int j = 0; j < NNj; j++) {
	intersect_north[i][j] = OFF;
	intersect_south[i][j] = OFF;
	intersect_east[i][j] = OFF;
	intersect_west[i][j] = OFF;
	Xnorth[i][j] = Vector2D(MILLION,MILLION);
	Xsouth[i][j] = Vector2D(MILLION,MILLION);
	Xeast[i][j] = Vector2D(MILLION,MILLION);
	Xwest[i][j] = Vector2D(MILLION,MILLION);
      }
    }
  }

  //! Reset intersect flags.
  void Reset_Intersect_Flags(const int &i, const int &j) {
    intersect_north[i][j] = OFF;
    intersect_south[i][j] = OFF;
    intersect_east[i][j] = OFF;
    intersect_west[i][j] = OFF;
  }

  //! Return the number of active computational cells.
  int Number_of_Active_Cells(void) {
    return (JCu-JCl+1)*(ICu-ICl+1) - Number_of_Inactive_Cells;
  }

  //! Reset interface present flag.
  void Reset_Interface_Present(void) {
    if (Ni) {
      for (int ni = 0; ni < Ni+1; ni++) Interface_Present[ni] = OFF;
    }
  }

};

/*!
 * class: EmbeddedBoundaries2D
 *
 * @brief Embedded boundaries solution class definition.
 *
 * This templated class contains all of the routines that are required
 * to conduct the solution of a 2D hyperbolic/elliptic system of
 * equations (e.g., the Euler equations or the Navier-Stokes equations) 
 * on a body-fitted multi-block quadrilateral mesh that may include
 * embedded boundaries that are not aligned to the underlying body-
 * fitted mesh.  Only explicit time-marching of the system of equations
 * can be performed here, however, FAS multigrid solutions can be done
 * as contained in EmbeddedBoundaries2D_FASMultigrid.h.  Block-based
 * adaptive mesh refinement is performed in the same manner as the 
 * individual equation-type solvers in CFFC.  In fact, this
 * templated class essentially replaces the ****2DQuadSolvers routine
 * used by the individual equation-type solvers in CFFC.  The
 * set-up of the problem and the output of the computed solution is 
 * performed in the routine found in EmbeddedBoundaries_Solvers.h, 
 * however, the explicit solution of the system of equations is found in
 * the 'Execute' function.
 * 
 * Reference:
 * Sachdev and Groth, "A Mesh Adjustment Scheme for Embedded Boundaries,"
 * under preparation for submission to Communications in Computational
 * Physics.
 *
 * \begin{verbatim}
 * \end{verbatim}
 */
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
class EmbeddedBoundaries2D{
public:

#ifdef _EB_PARALLEL_DEBUG_
  char extension[256], output_file_name[256];
  char *output_file_name_ptr;
  ofstream dout;
#endif

  //////////////////////////////////////////////////////////////////////
  // Member variables.                                                //
  //////////////////////////////////////////////////////////////////////

  //@{ @name Member variables
  //! Pointer to the input parameters.
  Quad_Soln_Input_Parameters *IP;

  //! Pointer to the quadtree.
  QuadTreeBlock_DataStructure *QuadTree;

  //! Pointer to the global solution block list.
  AdaptiveBlockResourceList *Global_Solution_Block_List;

  //! Pointer to the local solution block list.
  AdaptiveBlock2D_List *Local_Solution_Block_List;

  //! Pointer to the solution block information.
  Quad_Soln_Block *Local_SolnBlk;

  //! Copies of the grid required for grid management.
  Grid2D_Quad_Block *OGrid;
  Grid2D_Quad_Block *AGrid;

  //! Extra mesh information required for the mesh adjustment algorithm.
  Adjusted_Mesh_Quad_Block *Mesh;
  Adjusted_Mesh_Quad_Block *OMesh;
  Adjusted_Mesh_Quad_Block *AMesh;

  //! Interface component and union lists.
  Interface2D_List Interface_Component_List;
  Interface2D_List Interface_Union_List;

  //! Data required for the mesh adjustment algorithm.
  Mesh_Adjustment_Data *Adjustment_Data;

  //! LevelSet2D input parameters and multi-block solution-adaptive
  //! quadrilateral mesh solution variables.
  LevelSet2D_Input_Parameters  *LS_IP;
  QuadTreeBlock_DataStructure  *LS_QuadTree;
  AdaptiveBlockResourceList    *LS_Global_Solution_Block_List;
  AdaptiveBlock2D_List         *LS_Local_Solution_Block_List;
  LevelSet2D_Quad_Block        *LS_Local_SolnBlk;

  //! Reconstruction stencil:
  int *i_index, *j_index;
  //@}

  //////////////////////////////////////////////////////////////////////
  // Constructors, destructor, allocation, and deallocation.          //
  //////////////////////////////////////////////////////////////////////

  //@{ @name Constructors, destructor, allocation, and deallocation
  //! Default Constructor.
  EmbeddedBoundaries2D(void) {
    IP = NULL;
    QuadTree = NULL;
    Global_Solution_Block_List = NULL;
    Local_Solution_Block_List = NULL;
    Local_SolnBlk = NULL;
    OGrid = NULL;
    AGrid = NULL;
    Mesh = NULL;
    OMesh = NULL;
    AMesh = NULL;
    Adjustment_Data = NULL;
    LS_IP = NULL;
    LS_QuadTree = NULL;
    LS_Global_Solution_Block_List = NULL;
    LS_Local_Solution_Block_List = NULL;
    LS_Local_SolnBlk = NULL;
  }

  //! Default Destructor.
  //~EmbeddedBoundaries2D(void) { }//deallocate(); }

  //! Memory allocation and initialization.
  void allocate(Quad_Soln_Block *SolnBlk,
		QuadTreeBlock_DataStructure *quadtree,
		AdaptiveBlockResourceList *GlobalList,
		AdaptiveBlock2D_List *LocalList,
		Quad_Soln_Input_Parameters *ip,
		LevelSet2D_Quad_Block *LS_SolnBlk,
		QuadTreeBlock_DataStructure *LS_quadtree,
		AdaptiveBlockResourceList *LS_GlobalList,
		AdaptiveBlock2D_List *LS_LocalList,
		LevelSet2D_Input_Parameters *LS_ip);

  //! Memory deallocation.
  void deallocate(void);
  //@}

  //////////////////////////////////////////////////////////////////////
  // Interface construction and manipulation routines.                //
  //////////////////////////////////////////////////////////////////////

  //@{ @name Interface construction and manipulation routines
  //! Construct the interface component list.
  int Construct_Interface_Component_List(void);

  //! Construct the level set interface list and initialize level set solution.
  int Construct_Interface_Level_Set_List(const int &batch_flag);

  //! Reconstruct the interface component list.
  int Reconstruct_Interface_Component_List(void);

  //! Construct the interface union list.
  int Construct_Interface_Union_List(void);

  //! Write the interface component list to output data files.
  int Output_Interface_Component_List_Tecplot(void);

  //! Write the interface union list to output data files.
  int Output_Interface_Union_List_Tecplot(void);
  //@}

  //////////////////////////////////////////////////////////////////////
  // Mesh adjustment and management routines.                         //
  //////////////////////////////////////////////////////////////////////

  //@{ @name Mesh adjustment and management routines.
  //! Initialize grids and mesh data required for the adjustment algorithm.
  void Initialize_Adjustment_Grids(void);

  //! Store the adjusted grid.
  void Store_Adjusted_Mesh(void);

  //! Store the adjusted grid.
  void Store_Adjusted_Mesh(const int &nb);

  //! Store the unadjusted grid.
  void Store_Unadjusted_Mesh(void);

  //! Store the unadjusted grid.
  void Store_Unadjusted_Mesh(const int &nb);

  //! Restore the grid to its initial state.
  void Mesh_Unadjustment(void);

  //! Restore the grid to its initial state.
  void Mesh_Unadjustment(const int &nb);

  //! Perform pre-mesh adjustment cell and node painting.
  int Pre_Mesh_Adjustment_Painting(const int &nb);

  //! Perform post-mesh adjustment cell painting.
  int Post_Mesh_Adjustment_Painting(const int &nb);

  //! Conduct the mesh adjustment algorithm.
  int Mesh_Adjustment(const int &copy_flag,
		      const int &statistics_flag);

  //! Conduct the mesh adjustment algorithm.
  int Mesh_Adjustment_Sharp_Corners(const int &nb);
  int Mesh_Adjustment_First(const int &nb);
  int Mesh_Adjustment_Second(const int &nb);
  int Mesh_Adjustment_Third(const int &nb);
  int Mesh_Adjustment_Modify_Boundary_Splines(const int &nb);
  int Mesh_Adjustment_Finalize(const int &nb);

  //! Check if the current or distance-1 neighbour cells are 'unknown.'
  int Cell_Status_Test(const int &nb, const int &i, const int &j);

  //! Check if the mesh adjustment has created any unphysical cells.
  int Check_Cell_Areas(const int &nb, const int &i, const int &j);

  //! Determine the statistics of the mesh adjustment.
  void Mesh_Adjustment_Statistics(void);

  //! Determine the statistics of the mesh adjustment.
  void Mesh_Adjustment_Statistics(const int &nb,
				  int &number_of_adjusted_blocks,
				  int &number_inactive_cells,
				  int &total_number_of_cells,
				  double &cell_area_ratio);

  //! Update exterior cell information.
  void Update_Exterior_Cells(const int &nb);
  //@}

  //////////////////////////////////////////////////////////////////////
  // Moving embedded boundary routines.                               //
  //////////////////////////////////////////////////////////////////////

  //@{ @name Moving embedded boundary routines.
  int Compute_Interface_Location(const int &batch_flag,
				 const double &current_time,
				 const double &maximum_time,
				 int &evolution_counter,
				 int &number_of_level_set_time_steps,
				 double &level_set_current_time);

  int Compute_Interface_Velocity_Function(const double &Time);

  int Compute_Interface_Velocity_Function(const int &nb,
					  const int &Ni,
					  const int &Motion_Flag);

  int Redistribute_Solution_Content(const double &Time);

  int Redistribute_Solution_Content(const int &nb,
				    const double &Time);

  int Reset_Interface_Motion_Type(void);

  int Zero_Interface_Velocity_Function(void);
  //@}

  //////////////////////////////////////////////////////////////////////
  // Bilinear interplation routines.                                  //
  //////////////////////////////////////////////////////////////////////

  //@{ @name Bilinear interplation (Connell & Holmes AIAA Paper 1989-1932-CP).
  //! Return primitive solution state at the specified node.
  pState Wn(const int &nb, const int &ii, const int &jj);

  //! Return conserverd solution state the at specified node.
  cState Un(const int &nb, const int &ii, const int &jj);

  pState WnNW(const int &nb, const int &ii, const int &jj); //!< Return primitive solution state at the NW node.
  pState WnNE(const int &nb, const int &ii, const int &jj); //!< Return primitive solution state at the NE node.
  pState WnSE(const int &nb, const int &ii, const int &jj); //!< Return primitive solution state at the SE node.
  pState WnSW(const int &nb, const int &ii, const int &jj); //!< Return primitive solution state at the SW node.

  cState UnNW(const int &nb, const int &ii, const int &jj); //!< Return conserved solution state at the NW node.
  cState UnNE(const int &nb, const int &ii, const int &jj); //!< Return conserved solution state at the NE node.
  cState UnSE(const int &nb, const int &ii, const int &jj); //!< Return conserved solution state at the SE node.
  cState UnSW(const int &nb, const int &ii, const int &jj); //!< Return conserved solution state at the SW node.
  //@}

  //////////////////////////////////////////////////////////////////////
  // Adaptive mesh refinement routines.                               //
  //////////////////////////////////////////////////////////////////////

  //! Prolong the solution block information.
  int Prolong_Quadrilateral_Block(const int &Perform_Mesh_Adjustment,
				  const int &nb,
				  Quad_Soln_Block &SolnBlk_Original,
				  Grid2D_Quad_Block &AGrid_Original,
				  Adjusted_Mesh_Quad_Block &AMesh_Original,
				  const int &Sector);

  //! Restrict the solution block information.
  int Restrict_Quadrilateral_Block(const int &Perform_Mesh_Adjustment,
				   const int &nb,
				   Quad_Soln_Block &SolnBlk_Original_SW,
				   Quad_Soln_Block &SolnBlk_Original_SE,
				   Quad_Soln_Block &SolnBlk_Original_NW,
				   Quad_Soln_Block &SolnBlk_Original_NE,
				   Grid2D_Quad_Block &Agrid_Original_SW,
				   Grid2D_Quad_Block &Agrid_Original_SE,
				   Grid2D_Quad_Block &Agrid_Original_NW,
				   Grid2D_Quad_Block &Agrid_Original_NE,
				   Adjusted_Mesh_Quad_Block &Amesh_Original_SW,
				   Adjusted_Mesh_Quad_Block &Amesh_Original_SE,
				   Adjusted_Mesh_Quad_Block &Amesh_Original_NW,
				   Adjusted_Mesh_Quad_Block &Amesh_Original_NE);

  //! Performs the mesh refinement of the adaptive blocks.
  int Refine_Grid(const int &Perform_Mesh_Adjustment);

  //! Performs the mesh coarsening of the adaptive blocks.
  int Coarsen_Grid(const int &Perform_Mesh_Adjustment);

  //! Calculate refinement criteria for a specified solution block.
  void Calculate_Refinement_Criteria(const int &nb,
				     double *refinement_criteria);

  //! Calculate refinement criteria for a specific cell of a solution block.
  void Calculate_Refinement_Criteria(const int &nb, const int &i, const int &j,
				     double *refinement_criteria);

  //! Flag blocks for refinement.
  void Flag_Blocks_For_Refinement(void);

  //! Adaptive mesh refinement routine.
  int Adaptive_Mesh_Refinement(const int &Set_New_Refinement_Flags,
			       const int &Perform_Mesh_Adjustment);

  //! Initial adaptive mesh refinement routine.
  int Initial_Adaptive_Mesh_Refinement(void);

  //! Uniform adaptive mesh refinement routine.
  int Uniform_Adaptive_Mesh_Refinement(void);

  //! Boundary adaptive mesh refinement routine.
  int Boundary_Adaptive_Mesh_Refinement(void);

  //! Interface adaptive mesh refinement routine.
  int Interface_Adaptive_Mesh_Refinement(void);

  //! Bounding-box adaptive mesh refinement routine.
  int Bounding_Box_Adaptive_Mesh_Refinement(const int &Apply_ICs);

  //////////////////////////////////////////////////////////////////////
  // Input-output routines.                                           //
  //////////////////////////////////////////////////////////////////////

  //! Read the restart data files.
  int Read_Restart_Files(const int &batch_flag,
			 int &number_of_level_set_time_steps,
			 double &levelset_Time);

  //! Read the level set restart data files.
  int Read_Level_Set_Restart_Solution(const int &batch_flag,
				      int &number_of_level_set_time_steps,
				      double &levelset_Time);

  //! Write the restart data files.
  int Write_Restart_Files(int &number_of_time_steps,
			  int &number_of_level_set_time_steps,
			  double &Time,
			  double &levelset_Time,
			  CPUTime &processor_cpu_time);

  //! Write the level set restart data files.
  int Write_Level_Set_Restart_Solution(int &number_of_level_set_time_steps,
				       double &Time,
				       CPUTime &processor_cpu_time);

  //! Write the solution elements to output data files.
  int Output_Elements_Tecplot(const int &number_of_time_steps,
			      const double &Time);

  //! Write the active solution elements to output data files.
  void Output_Active_Elements_Tecplot(const int &nb,
				      const int &Number_of_Time_Steps,
				      const double &Time,
				      const int &Output_Title,
				      ostream &Out_File);

  //! Write the active solution elements to output data files.
  void Output_Inactive_Elements_Tecplot(const int &nb,
					const int &Number_of_Time_Steps,
					const double &Time,
					const int &Output_Title,
					ostream &Out_File);

  //! Write the domain nodes to output data files.
  int Output_Nodes_Tecplot(void);

  //! Write the domain nodes to output data files.
  void Output_Nodes_Tecplot(const int &nb,
			    const int &Output_Title,
			    ostream &Out_File);

  //! Write the domain cell status to output data files.
  int Output_Cell_Status_Tecplot(void);

  //! Write the domain cell status to output data files.
  void Output_Cell_Status_Tecplot(const int &nb,
				  const int &Output_Title,
				  ostream &Out_File);

  //! Write the level set solution to output data files.
  int Output_Level_Set_Tecplot(const int &number_of_time_steps,
			       const double &Time);

  //! Write the level set cell-centred solution to output data files.
  int Output_Level_Set_Cells_Tecplot(const int &number_of_time_steps,
				     const double &Time);

  //! Write the level set interface solution to output data files.
  int Output_Level_Set_Interface_Tecplot(const int &number_of_time_steps,
					 const double &Time);

  //! Write the Ringleb's flow output data files.
  int Output_Ringleb_Tecplot(void);

  //! Write the Ringleb's flow output data files.
  void Output_Ringleb_Tecplot(const int &nb,
			      const int &Output_Title,
			      ostream &Out_File,
			      double &l1_norm,
			      double &l2_norm,
			      double &max_norm,
			      double &area);

  //! Write the flat plate output data files.
  int Output_Flat_Plate_Tecplot(void);

  //! Write the flat plate output data files.
  void Output_Flat_Plate_Tecplot(const int &nb,
				 const int &Output_Title_Soln,
				 ostream &Out_File_Soln,
				 const int &Output_Title_Skin,
				 ostream &Out_File_Skin,
				 double &l1_norm,
				 double &l2_norm,
				 double &max_norm,
				 double &area,
				 int &numberofactivecells,
				 double &l1_norm_cf,
				 double &l2_norm_cf,
				 double &max_norm_cf,
				 double &area_cf,
				 int &numberofactivecells_cf);
  //! Calculate net force.
  Vector2D Net_Force(void);

  //! Calculate net pressure force.
  Vector2D Net_Pressure_Force(void);

  //! Output Data for Couette Case.
  int Output_Couette(void);

  //! Calculate Cl and Cd for cylinder cases.
  int Output_Cylinder_Drag(void);

  //! Determine and write aerodynamic coefficients.
  int Output_Aerodynamic_Coefficients_Tecplot(const int &number_of_time_steps,
					      const double &Time);

  //! Determine and write aerodynamic coefficients.
  void Output_Aerodynamic_Coefficients_Tecplot(const int &nb,
					       const int &Output_Title,
					       ostream &Out_File,
					       double &Cn, double &Ca, double &Cm,
					       double &Cnp, double &Cmp);

  //! Ywall
  double Ywall(const int &nb, const int &i, const int &j);

  //! Yplus
  double Yplus(const int &nb, const int &i, const int &j);

  //////////////////////////////////////////////////////////////////////
  // Solution routines.                                               //
  //////////////////////////////////////////////////////////////////////

  //! Apply boundary conditions.
  int Boundary_Conditions(const double &Time);

  //! Apply interface boundary conditions.
  int BCs_Interface(const int &nb, const double &Time);

  //! Determine the allowable global and local time-steps according to
  //! the Courant-Friedrichs-Lewy condition.
  double CFL(const double &Time);

  //! Determine the L1-norm of the solution residual.
  double L1_Norm_Residual(void);

  //! Determine the L2-norm of the solution residual.
  double L2_Norm_Residual(void);

  //! Determine the max-norm of the solution residual.
  double Max_Norm_Residual(void);

  //! Determine the rate of area change using the geometric conservation law.
  double dAdt(const int &nb, const int &i, const int &j, const double &Time);

  //! Perform linear least squares reconstruction on the specified cell
  //! of the specified solution block.
  void Linear_Least_Squares_Reconstruction(const int &nb,
					   const int &i,
					   const int &j);

  //! Perform linear least squares reconstruction on each solution block.
  void Linear_Least_Squares_Reconstruction(const int &nb);

  //! Perform residual smoothing iterations.
  void Residual_Smoothing(const int &i_stage);

  //! Determine the solution residual for the specified stage of the
  //! n-stage time-steppings scheme.
  int dUdt_Residual_Evaluation(const double &Time);

  //! Determine the solution residual for the specified stage of the
  //! n-stage time-steppings scheme.
  int dUdt_Multistage_Explicit(const int &i_stage, const double &Time);

  //! Update the solution using a multi-stage explicit time-stepping scheme.
  int Update_Solution_Multistage_Explicit(const int &i_stage);

  //! Execute the FAS multigrid algorithm for steady state computations.
  int Execute(const int &batch_flag,
	      int &number_of_time_steps,
	      int &evolution_counter,
	      int &number_of_level_set_time_steps,
	      double &Time,
	      double &levelset_Time,
	      CPUTime &processor_cpu_time,
	      CPUTime &total_cpu_time,
	      ofstream &residual_file);

};

/**********************************************************************
 **********************************************************************
 ** Memory allocation and deallocation routines.                     **
 **********************************************************************
 **********************************************************************/

/**********************************************************************
 * EmbeddedBoundaries2D::allocate --                                  *
 *                                                                    *
 * This routine performs the memory allocation and initialization for *
 * all grid levels of the FAS multigrid solution class.               *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
void EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
allocate(Quad_Soln_Block *SolnBlk,
	 QuadTreeBlock_DataStructure *quadtree,
	 AdaptiveBlockResourceList *GlobalList,
	 AdaptiveBlock2D_List *LocalList,
	 Quad_Soln_Input_Parameters *ip,
	 LevelSet2D_Quad_Block *LS_SolnBlk,
	 QuadTreeBlock_DataStructure *LS_quadtree,
	 AdaptiveBlockResourceList *LS_GlobalList,
	 AdaptiveBlock2D_List *LS_LocalList,
	 LevelSet2D_Input_Parameters *LS_ip) {

  int error_flag;

#ifdef _EB_PARALLEL_DEBUG_
  strcpy(output_file_name,"debug");
  strcat(output_file_name,"_cpu");
  sprintf(extension,"%.6d",LocalList->ThisCPU);
  strcat(extension,".txt");
  strcat(output_file_name,extension);
  output_file_name_ptr = output_file_name;
  dout.open(output_file_name_ptr,ios::out);
  if (dout.bad()) return ;
#endif

  // Point the input parameters to the given input parameters.
  IP = ip;

  // Point the quadtree to the input quadtree.
  QuadTree = quadtree;

  // Point the global solution block list to the input list.
  Global_Solution_Block_List = GlobalList;

  // Point the local solution block list to the input list.
  Local_Solution_Block_List = LocalList;

  // Point the local solution block to the input local solution block.
  Local_SolnBlk = SolnBlk;

  // Allocate memory for the extra grid information.
  AGrid = new Grid2D_Quad_Block[IP->Number_of_Blocks_Per_Processor];
  OGrid = new Grid2D_Quad_Block[IP->Number_of_Blocks_Per_Processor];

  // Allocate memory for the adjusted mesh data.
  Mesh = new Adjusted_Mesh_Quad_Block[IP->Number_of_Blocks_Per_Processor];
  AMesh = new Adjusted_Mesh_Quad_Block[IP->Number_of_Blocks_Per_Processor];
  OMesh = new Adjusted_Mesh_Quad_Block[IP->Number_of_Blocks_Per_Processor];

  // Allocate memory for the mesh adjustment data.
  Adjustment_Data = new Mesh_Adjustment_Data[IP->Number_of_Blocks_Per_Processor];

  // Point the level set input parameters to the given level set 
  // input parameters.
  LS_IP = LS_ip;

  // Point the level set quadtree to the input level set quadtree.
  LS_QuadTree = LS_quadtree;

  // Point the level set global solution block list to the input list.
  LS_Global_Solution_Block_List = LS_GlobalList;

  // Point the level set local solution block list to the input list.
  LS_Local_Solution_Block_List = LS_LocalList;

  // Point the level set local solution block to the input local solution block.
  LS_Local_SolnBlk = LS_SolnBlk;

  // Allocate memory for the reconstruction stencil:
  i_index = new int[8];
  j_index = new int[8];

}

/**********************************************************************
 * EmbeddedBoundaries2D::deallocate --                                *
 *                                                                    *
 * This routine performs the memory deallocation for the FAS          *
 * multigrid solution class.                                          *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
void EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
deallocate(void) {

  // Deallocate the extra grid and mesh memory.
  if (AGrid != NULL) {
    for (int nb = 0; nb < IP->Number_of_Blocks_Per_Processor; nb++) {
      AGrid[nb].deallocate();
      OGrid[nb].deallocate();
      Mesh[nb].deallocate();
      AMesh[nb].deallocate();
      OMesh[nb].deallocate();
      Adjustment_Data[nb].deallocate();
    }
    delete []AGrid; AGrid = NULL;
    delete []OGrid; OGrid = NULL;
    delete []Mesh; Mesh = NULL;
    delete []AMesh; AMesh = NULL;
    delete []OMesh; OMesh = NULL;
    delete []Adjustment_Data; Adjustment_Data = NULL;
  }

  // Deallocate memory for the list of interface components.
  Interface_Component_List.deallocate();

  // Deallocate memory for the list of interface unions.
  Interface_Union_List.deallocate();

  // Point the input parameters pointer to null.
  if (IP != NULL) { IP = NULL; }

  // Point the quadtree pointer to null.
  if (QuadTree != NULL) { QuadTree = NULL; }

  // Point the global solution block list pointer to null.
  if (Global_Solution_Block_List != NULL) { 
    Global_Solution_Block_List = NULL;
  }

  // Point the local solution block list pointer to null.
  if (Local_Solution_Block_List != NULL) { 
    Local_Solution_Block_List = NULL;
  }

  // Point the local solution block pointer to null.
  if (Local_SolnBlk != NULL) { Local_SolnBlk = NULL; }

  // Deallocate memory for level set variables if required.
//   if (LS_Local_SolnBlk != NULL) {
//     if (CFFC_Primary_MPI_Processor()) Close_Input_File(LS_IP);
//     LS_Local_SolnBlk = Deallocate(LS_Local_SolnBlk,LS_IP);
//     LS_Global_Solution_Block_List.deallocate();
//     LS_Local_Solution_Block_List.deallocate();
//     LS_QuadTree.deallocate();
//   }
  if (LS_IP != NULL) { LS_IP = NULL; }
  if (LS_QuadTree != NULL) { LS_QuadTree = NULL; }
  if (LS_Global_Solution_Block_List != NULL) { 
    LS_Global_Solution_Block_List = NULL;
  }
  if (LS_Local_Solution_Block_List != NULL) { 
    LS_Local_Solution_Block_List = NULL;
  }
  if (LS_Local_SolnBlk != NULL) { LS_Local_SolnBlk = NULL; }

  // Deallocate memory for the reconstruction stencil:
  delete []i_index; i_index = NULL;
  delete []j_index; j_index = NULL;

#ifdef _EB_PARALLEL_DEBUG_
  dout.close();
#endif

}

/**********************************************************************
 **********************************************************************
 ** Interface construction and manipulation routines.                **
 **********************************************************************
 **********************************************************************/

/**********************************************************************
 * EmbeddedBoundaries2D::Construct_Interface_Component_List --        *
 *                                                                    *
 * This routine constructs all of the embedded boundaries specified   *
 * by the input parameters.  These interface components are stored in *
 * as a list.  Each processor stores the entire list and this list is *
 * accessible by every solution block owned by that processor.        *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Construct_Interface_Component_List(void) {

  double theta, Angle1, Angle2;
  char NASA_Rotor_Data_Directory[128];
  NASARotor37 NASA_Rotor_37;
  NASARotor67 NASA_Rotor_67;
  int Rotor_Flow_Type;

  // Allocate memory for interface variables.
  if (Interface_Component_List.Ni > 0) Interface_Component_List.deallocate();
  Interface_Component_List.allocate(IP->Interface_IP.Component_List.Ni);
  // Assign interface data.
  for (int ni = 1; ni <= Interface_Component_List.Ni; ni++) {
    // Copy interface type.
    Interface_Component_List[ni].Type    = IP->Interface_IP.Component_List[ni].Type;
    // Copy interface lengths.
    Interface_Component_List[ni].Length1 = IP->Interface_IP.Component_List[ni].Length1;
    Interface_Component_List[ni].Length2 = IP->Interface_IP.Component_List[ni].Length2;
    // Copy interface BC type.
    Interface_Component_List[ni].BC_Type = IP->Interface_IP.Component_List[ni].BC_Type;
    // Copy interface motion type.
    Interface_Component_List[ni].Motion  = IP->Interface_IP.Component_List[ni].Motion;
    // Copy interface motion speed.
    Interface_Component_List[ni].Speed   = IP->Interface_IP.Component_List[ni].Speed;
    // Copy interface reference point.
    Interface_Component_List[ni].Xref    = IP->Interface_IP.Component_List[ni].Xref;
  }

  // Construct each interface of the interface component list.
  for (int ni = 1; ni <= Interface_Component_List.Ni; ni++) {

    // Deallocate memory for interface component if already allocated.
    if (Interface_Component_List[ni].Spline.np != 0)
      Interface_Component_List[ni].Spline.deallocate();

    // Construct the interface component.
    switch(Interface_Component_List[ni].Type) {
    case INTERFACE_LINE :
      // Create line
      Create_Spline_Line(Interface_Component_List[ni].Spline,
			 IP->Interface_IP.Component_List[ni].Spline.Xp[0]-Vector2D(Interface_Component_List[ni].Length1/2.0,0.0),
			 IP->Interface_IP.Component_List[ni].Spline.Xp[0]+Vector2D(Interface_Component_List[ni].Length1/2.0,0.0),
			 2);
      // Set boundary condition.
      Interface_Component_List[ni].Spline.setBCtype(Interface_Component_List[ni].BC_Type);
      // Initialize the interface velocity function.
      Interface_Component_List[ni].Initialize_Velocity_Function(Interface_Component_List[ni].Motion);
      // Determine the reference point (centroid) of the interface.
      //Interface_Component_List[ni].Centroid();
      break;
    case INTERFACE_CIRCLE :
      // Create circular spline.
      Create_Spline_Circular_Arc(Interface_Component_List[ni].Spline,
				 IP->Interface_IP.Component_List[ni].Spline.Xp[0],
				 Interface_Component_List[ni].Length1,
				 ZERO,
				 360.0,
				 361);
      // Set boundary condition.
      Interface_Component_List[ni].Spline.setBCtype(Interface_Component_List[ni].BC_Type);
      // Initialize the interface velocity function.
      Interface_Component_List[ni].Initialize_Velocity_Function(Interface_Component_List[ni].Motion);
      // Determine the reference point (centroid) of the interface.
      Interface_Component_List[ni].Centroid();
      break;
    case INTERFACE_ELLIPSE :
      // Create circular spline.
      Create_Spline_Ellipsoidal_Arc(Interface_Component_List[ni].Spline,
				    IP->Interface_IP.Component_List[ni].Spline.Xp[0],
				    Interface_Component_List[ni].Length1,
				    Interface_Component_List[ni].Length2,
				    ZERO,
				    360.0,
				    361);
      // Set boundary condition.
      Interface_Component_List[ni].Spline.setBCtype(Interface_Component_List[ni].BC_Type);
      // Initialize the interface velocity function.
      Interface_Component_List[ni].Initialize_Velocity_Function(Interface_Component_List[ni].Motion);
      // Determine the reference point (centroid) of the interface.
      Interface_Component_List[ni].Centroid();
      break;
    case INTERFACE_SQUARE :
      // Create square spline.
      Create_Spline_Rectangle(Interface_Component_List[ni].Spline,
			      IP->Interface_IP.Component_List[ni].Spline.Xp[0],
			      Interface_Component_List[ni].Length1,
			      Interface_Component_List[ni].Length1);
      // Set boundary condition.
      Interface_Component_List[ni].Spline.setBCtype(Interface_Component_List[ni].BC_Type);
      // Initialize the interface velocity function.
      Interface_Component_List[ni].Initialize_Velocity_Function(Interface_Component_List[ni].Motion);
      // Determine the reference point (centroid) of the interface.
      Interface_Component_List[ni].Centroid();
      break;
    case INTERFACE_RECTANGLE :
      // Create square spline.
      Create_Spline_Rectangle(Interface_Component_List[ni].Spline,
			      IP->Interface_IP.Component_List[ni].Spline.Xp[0],
			      Interface_Component_List[ni].Length1,
			      Interface_Component_List[ni].Length2);
      // Set boundary condition.
      Interface_Component_List[ni].Spline.setBCtype(Interface_Component_List[ni].BC_Type);
      // Initialize the interface velocity function.
      Interface_Component_List[ni].Initialize_Velocity_Function(Interface_Component_List[ni].Motion);
      // Determine the reference point (centroid) of the interface.
      Interface_Component_List[ni].Centroid();
      break;
    case INTERFACE_ROCKET_PROPELLANT_GRAIN :
      // Create solid propellant spline.
      // Allocate memory for the circular arc spline. 
      Interface_Component_List[ni].Spline.allocate(3);
      // Set the spline type.
      Interface_Component_List[ni].Spline.settype(SPLINE2D_LINEAR);
      // Set the spline points.
      Interface_Component_List[ni].Spline.Xp[0] = IP->Interface_IP.Component_List[ni].Spline.Xp[0];
      Interface_Component_List[ni].Spline.Xp[1] = IP->Interface_IP.Component_List[ni].Spline.Xp[1];
      Interface_Component_List[ni].Spline.Xp[2] = IP->Interface_IP.Component_List[ni].Spline.Xp[2];
      // Set spline tp and bc.
      Interface_Component_List[ni].Spline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
      Interface_Component_List[ni].Spline.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
      Interface_Component_List[ni].Spline.tp[2] = SPLINE2D_POINT_SHARP_CORNER;
      for (int j = 0; j < 3; j++)
	Interface_Component_List[ni].Spline.bc[j] = BC_NONE;
      // Calculate the spline pathlengths.
      Interface_Component_List[ni].Spline.pathlength();
      // Set boundary condition.
      Interface_Component_List[ni].Spline.setBCtype(Interface_Component_List[ni].BC_Type);
      // Initialize the interface velocity function.
      Interface_Component_List[ni].Initialize_Velocity_Function(Interface_Component_List[ni].Motion);
      break;
    case INTERFACE_NACA0012_AEROFOIL :
      // Create NACA aerofoil spline.
      Create_Spline_NACA_Aerofoil(Interface_Component_List[ni].Spline,
				  "0012",
				  Interface_Component_List[ni].Length1,
				  0,
				  501);
      // Determine the reference point (centroid) of the interface.
      Interface_Component_List[ni].Xref = 0.25*(3*Interface_Component_List[ni].Spline.Xp[250] +
						Interface_Component_List[ni].Spline.Xp[0]);
      // Translate the aerofoil such that the centroid is located at (0,0).
      Interface_Component_List[ni].Translate(-Interface_Component_List[ni].Xref);
      //Interface_Component_List[ni].Translate(IP->Interface_IP.Component_List[ni].Xref);
      // Rotate the aerofoil such that...
      Interface_Component_List[ni].Rotate(IP->Interface_IP.Component_List[ni].Spline.Xp[0].x*PI/180.0);
      // Set boundary condition.
      Interface_Component_List[ni].Spline.setBCtype(Interface_Component_List[ni].BC_Type);
      // Initialize the interface velocity function.
      Interface_Component_List[ni].Initialize_Velocity_Function(Interface_Component_List[ni].Motion);
      break;
    case INTERFACE_NACA0015_AEROFOIL :
      // Create NACA aerofoil spline.
      Create_Spline_NACA_Aerofoil(Interface_Component_List[ni].Spline,
				  "0015",
				  Interface_Component_List[ni].Length1,
				  0,
				  501);
      // Determine the reference point (centroid) of the interface.
      Interface_Component_List[ni].Xref = 0.25*(3*Interface_Component_List[ni].Spline.Xp[250] +
						Interface_Component_List[ni].Spline.Xp[0]);
      // Translate the aerofoil such that the centroid is located at (0,0).
      Interface_Component_List[ni].Translate(-Interface_Component_List[ni].Xref);
      //Interface_Component_List[ni].Translate(IP->Interface_IP.Component_List[ni].Xref);
      // Rotate the aerofoil such that...
      Interface_Component_List[ni].Rotate(IP->Interface_IP.Component_List[ni].Spline.Xp[0].x*PI/180.0);
      // Set boundary condition.
      Interface_Component_List[ni].Spline.setBCtype(Interface_Component_List[ni].BC_Type);
      // Initialize the interface velocity function.
      Interface_Component_List[ni].Initialize_Velocity_Function(Interface_Component_List[ni].Motion);
      break;
    case INTERFACE_NASA_ROTOR_37 :
      // Create NASA Rotor 37 spline.
      // Set default data values.
      strcpy(NASA_Rotor_Data_Directory,"/nfs/fe01/d1/cfd/jai/CFFC/data/NASA_Rotors/R37/");
      Rotor_Flow_Type = PEAK_FLOW;
      // Initialize NASA rotor 37 class.
      NASA_Rotor_37.init(Rotor_Flow_Type,NASA_Rotor_Data_Directory);
      // Get the 2D NASA Rotor 37 spline at the given percent span.
      Interface_Component_List[ni].Spline = NASA_Rotor_37.getBladeCS(Interface_Component_List[ni].Length1);
      // Determine reference point.
      Interface_Component_List[ni].Centroid();
      // Set boundary condition.
      Interface_Component_List[ni].Spline.setBCtype(Interface_Component_List[ni].BC_Type);
      // Initialize the interface velocity function.
      Interface_Component_List[ni].Initialize_Velocity_Function(Interface_Component_List[ni].Motion);
      break;
    case INTERFACE_NASA_ROTOR_67 :
      // Create NASA Rotor 67 spline.
      // Set default data values.
      strcpy(NASA_Rotor_Data_Directory,"/nfs/fe01/d1/cfd/jai/CFFC/data/NASA_Rotors/R67/");
      Rotor_Flow_Type = PEAK_FLOW;
      // Initialize NASA rotor 67 class.
      NASA_Rotor_67.init(Rotor_Flow_Type,NASA_Rotor_Data_Directory);
      // Get the 2D NASA Rotor 67 spline at the given percent span.
      Interface_Component_List[ni].Spline = NASA_Rotor_67.getBladeCS(Interface_Component_List[ni].Length1);
      // Determine reference point.
      Interface_Component_List[ni].Centroid();
      // Set flow-field.
      NASA_Rotor_67.getPstateREL_up(IP->Wo,IP->Rotor_Percent_Span);
      IP->Pressure = IP->Wo.pressure(); //I had to change this because p isn't a double for me. ~james
      IP->Temperature = IP->Wo.T();
      IP->Mach_Number = NASA_Rotor_67.getMachREL_up(IP->Rotor_Percent_Span);
      IP->Flow_Angle = atan2(IP->Wo.v.y,IP->Wo.v.x); 
      if (IP->Flow_Angle < ZERO) IP->Flow_Angle = TWO*PI + IP->Flow_Angle;
      IP->Flow_Angle = 180.00*IP->Flow_Angle/PI;
      NASA_Rotor_67.getPstateREL_down(IP->W1,IP->Rotor_Percent_Span);
      // Set boundary condition.
      Interface_Component_List[ni].Spline.setBCtype(Interface_Component_List[ni].BC_Type);
      // Initialize the interface velocity function.
      Interface_Component_List[ni].Initialize_Velocity_Function(Interface_Component_List[ni].Motion);
      break;
    case INTERFACE_ZALESAK :
      // Create Zalesak's Disk.
      Interface_Component_List[ni].Zalesak(IP->Interface_IP.Component_List[ni].Spline.Xp[0],
					   IP->Interface_IP.Component_List[ni].Length1);
      // Set boundary condition.
      Interface_Component_List[ni].Spline.setBCtype(Interface_Component_List[ni].BC_Type);
      // Initialize the interface velocity function.
      Interface_Component_List[ni].Initialize_Velocity_Function(Interface_Component_List[ni].Motion);
      break;
    case INTERFACE_RINGLEB :
      // Create Rinleb's flow domain.
      Interface_Component_List[ni].Ringleb(IP->Inner_Streamline_Number,
					   IP->Outer_Streamline_Number,
					   IP->Isotach_Line);
      // Initialize the interface velocity function.
      Interface_Component_List[ni].Initialize_Velocity_Function(Interface_Component_List[ni].Motion);
      break;
    case INTERFACE_BUMP :
      // Create the bump flow interface.
      Interface_Component_List[ni].Bump_Channel_Flow();
      // Set boundary condition.
      Interface_Component_List[ni].Spline.setBCtype(Interface_Component_List[ni].BC_Type);
      // Initialize the interface velocity function.
      Interface_Component_List[ni].Initialize_Velocity_Function(Interface_Component_List[ni].Motion);
      // Determine the reference point (centroid) of the interface.
      Interface_Component_List[ni].Centroid();
      break;
    case INTERFACE_FLAT_PLATE :
      // Create the flat plate interface.
      Interface_Component_List[ni].Flat_Plate();
      // Set boundary condition.
      Interface_Component_List[ni].Spline.setBCtype(Interface_Component_List[ni].BC_Type);
      // Initialize the interface velocity function.
      Interface_Component_List[ni].Initialize_Velocity_Function(Interface_Component_List[ni].Motion);
      break;
    case INTERFACE_USER_SPECIFIED :
      // Create user specified spline.
      // Allocate memory for the circular arc spline.
      Interface_Component_List[ni].Spline.allocate(IP->Interface_IP.Component_List[ni].Spline.np);
      // Set the spline type.
      Interface_Component_List[ni].Spline.settype(SPLINE2D_LINEAR);
      // Set the spline points, tp, and bc.
      for (int i = 0; i < Interface_Component_List[ni].Spline.np; i++) {
	Interface_Component_List[ni].Spline.Xp[i] = IP->Interface_IP.Component_List[ni].Spline.Xp[i];
	Interface_Component_List[ni].Spline.tp[i] = IP->Interface_IP.Component_List[ni].Spline.tp[i];
	Interface_Component_List[ni].Spline.bc[i] = IP->Interface_IP.Component_List[ni].Spline.bc[i];
      }
      Interface_Component_List[ni].Spline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
      Interface_Component_List[ni].Spline.tp[Interface_Component_List[ni].Spline.np-1] = SPLINE2D_POINT_SHARP_CORNER;
      // Calculate the spline pathlengths.
      Interface_Component_List[ni].Spline.pathlength();
      // Rotate interface.
      //Interface_Component_List[ni].Rotate(45.0*PI/180.0);
      // Set boundary condition.
      Interface_Component_List[ni].Spline.setBCtype(Interface_Component_List[ni].BC_Type);
      // Initialize the interface velocity function.
      Interface_Component_List[ni].Initialize_Velocity_Function(Interface_Component_List[ni].Motion);
      break;
    };
    // Sort the interface into counter-clockwise order.
    Interface_Component_List[ni].Sort();
    // Determine the bounding box for the interface.
    Interface_Component_List[ni].BoundingBox();
    // Set the interface velocity function.
    Interface_Component_List[ni].Set_Velocity_Function(ZERO);

  }

  // Interface component list successfully constructed.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Construct_Interface_Level_Set_List --        *
 *                                                                    *
 * This routine constructs a list of interfaces that are to be        *
 * evolved using the level set method.                                *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Construct_Interface_Level_Set_List(const int &batch_flag) {

  int error_flag, Ni;
  Interface2D_List Interface_List;

  // Determine the number of interfaces that are evolved using the
  // level set method.
  Ni = 0;
  for (int ni = 1; ni <= Interface_Component_List.Ni; ni++) {
    if (Interface_Component_List[ni].Motion == INTERFACE_MOTION_LEVELSET_STATIONARY ||
	Interface_Component_List[ni].Motion == INTERFACE_MOTION_LEVELSET_EXPAND ||
	Interface_Component_List[ni].Motion == INTERFACE_MOTION_LEVELSET_STRETCH ||
	Interface_Component_List[ni].Motion == INTERFACE_MOTION_LEVELSET_BULKFLOW ||
	Interface_Component_List[ni].Motion == INTERFACE_MOTION_BURNING_SURFACE) {
      Ni++;
    }
  }

  // Exit immediately if none of the interfaces are evolved using the
  // level set method.
  if (!Ni) return 0;

  // Create an interface component list that contains interface
  // components that are to be evolved using the level set method.
  Interface_List.allocate(Ni);
  Ni = 1;
  for (int ni = 1; ni <= Interface_Component_List.Ni; ni++) {
    if (Interface_Component_List[ni].Motion == INTERFACE_MOTION_LEVELSET_STATIONARY ||
	Interface_Component_List[ni].Motion == INTERFACE_MOTION_LEVELSET_EXPAND ||
	Interface_Component_List[ni].Motion == INTERFACE_MOTION_LEVELSET_STRETCH ||
	Interface_Component_List[ni].Motion == INTERFACE_MOTION_LEVELSET_BULKFLOW ||
	Interface_Component_List[ni].Motion == INTERFACE_MOTION_BURNING_SURFACE) {
      Interface_List[Ni].Copy(Interface_Component_List[ni]);
      Ni++;
    }
  }

  // Initialize the level set solution class.  Pass the interface
  // component list to the level set solution.
  LS_Local_SolnBlk = Initialize_Level_Set_Solution(IP->Input_File_Name,
						   batch_flag,
						   error_flag,
						   LS_Local_SolnBlk,
						   *LS_IP,
						   *LS_QuadTree,
						   *LS_Global_Solution_Block_List,
						   *LS_Local_Solution_Block_List,
						   Interface_List);
  if (error_flag) {
    error_flag = Output_Level_Set_Tecplot(0,ZERO);
    error_flag = Output_Level_Set_Cells_Tecplot(0,ZERO);
    return 1;//error_flag;
  }

  // Merge the universal level set interface list with the component 
  // list.  This must be computed before the union list can be
  // determined.
  error_flag = Reconstruct_Interface_Component_List();
  if (error_flag) return error_flag;

  // Destroy the temporary interface list.
  Interface_List.deallocate();

  // Interface level set list constructed successfully.
  return 0;

}

/**********************************************************************
 * Routine: Reconstruct_Interface_Component_List                      *
 *                                                                    *
 * This routine reconstructs the interface component list after the   *
 * motion of the interfaces that are evolved using the level set      *
 * method.  Note that in the level set method, interfaces may split   *
 * or merge, thereby changing the total number of interface           *
 * components.                                                        *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Reconstruct_Interface_Component_List(void) {

  // If the level set data has not been allocated then exit immediately.
  if (LS_Local_SolnBlk == NULL) return 0;

  int Ni, Type, BC_Type, Motion;
  Vector2D Xref, Speed;
  Interface2D_List Interface_List, New_Component_List;

  // Create a 'universal' interface list that contains all of the 
  // interfaces evolved using the level set method.
  if (CFFC_Primary_MPI_Processor()) {
    if (LS_Local_Solution_Block_List->Block[0].used == ADAPTIVEBLOCK2D_USED) {
      // Determine the number of interfaces evolved using the level
      // set method.
      Ni = LS_Local_SolnBlk[0].Interface_List.Ni;
      // Otherwise create universal interface list.
      Interface_List.Copy(LS_Local_SolnBlk[0].Interface_List);
    }
  }
  Broadcast_Interface_List(Interface_List);

  // If there are no interfaces extracted from the level set method
  // then exit immediately.
  if (!Interface_List.Ni) return 0;

  // Merge the interface list found from the level set program with
  // the old component list.
  Ni = Interface_List.Ni;
  for (int ni = 1; ni <= Interface_Component_List.Ni; ni++) {
    if (Interface_Component_List[ni].Motion != INTERFACE_MOTION_LEVELSET_STATIONARY &&
	Interface_Component_List[ni].Motion != INTERFACE_MOTION_LEVELSET_EXPAND &&
	Interface_Component_List[ni].Motion != INTERFACE_MOTION_LEVELSET_STRETCH &&
	Interface_Component_List[ni].Motion != INTERFACE_MOTION_LEVELSET &&
	Interface_Component_List[ni].Motion != INTERFACE_MOTION_LEVELSET_BULKFLOW &&
	Interface_Component_List[ni].Motion != INTERFACE_MOTION_BURNING_SURFACE) {
      Ni++;
    } else {
      Type = Interface_Component_List[ni].Type;
      Xref = Interface_Component_List[ni].Xref;
      BC_Type = Interface_Component_List[ni].BC_Type;
      Motion = Interface_Component_List[ni].Motion;
      Speed = Interface_Component_List[ni].Speed;
    }
  }

  // Allocate memory for the new interface component list.
  New_Component_List.allocate(Ni);

  // Reset motion type.
  if (Motion == INTERFACE_MOTION_LEVELSET_STATIONARY ||
      Motion == INTERFACE_MOTION_LEVELSET_EXPAND ||
      Motion == INTERFACE_MOTION_LEVELSET_STRETCH) {
    Motion = INTERFACE_MOTION_LEVELSET;
  }

  // Create new interface component list.
  for (int ni = 1; ni <= Interface_List.Ni; ni++) {
    New_Component_List[ni].Copy(Interface_List[ni]);
    New_Component_List[ni].Type = Type;
    New_Component_List[ni].Xref = Xref;
    New_Component_List[ni].BC_Type = BC_Type;
    New_Component_List[ni].Spline.setBCtype(BC_Type);
    New_Component_List[ni].Motion = Motion;
    New_Component_List[ni].Speed = Speed;
    New_Component_List[ni].Set_Velocity_Function_Type();
    //New_Component_List[ni].Set_Normal_Velocity_Function();
    if (abs(New_Component_List[ni].Spline.Xp[0] -
	    New_Component_List[ni].Spline.Xp[New_Component_List[ni].Spline.np-1]) < NANO) {
      Interface_Component_List[ni].Centroid();
    }
    New_Component_List[ni].BoundingBox();
  }

  Ni = 0;
  for (int ni = 1; ni <= Interface_Component_List.Ni; ni++) {
    if (Interface_Component_List[ni].Motion != INTERFACE_MOTION_LEVELSET_STATIONARY &&
	Interface_Component_List[ni].Motion != INTERFACE_MOTION_LEVELSET_EXPAND &&
	Interface_Component_List[ni].Motion != INTERFACE_MOTION_LEVELSET_STRETCH &&
	Interface_Component_List[ni].Motion != INTERFACE_MOTION_LEVELSET_BULKFLOW &&
	Interface_Component_List[ni].Motion != INTERFACE_MOTION_LEVELSET &&
	Interface_Component_List[ni].Motion != INTERFACE_MOTION_BURNING_SURFACE) {
      Ni++;
      New_Component_List[Interface_List.Ni+Ni].Copy(Interface_Component_List[ni]);
    }
  }

  // Copy new interface component list to the working interface
  // component list.
  Interface_Component_List.deallocate();
  Interface_Component_List.Copy(New_Component_List);

  // Interface component list successfully reconstructed.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Construct_Interface_Union_List --            *
 *                                                                    *
 * This routine constructs a list of interfaces that is comprised of  *
 * the interfaces in the component list, however, components that     *
 * interesect are combined into a union interface.  The mesh is       *
 * adjusted according to the interfaces in the union list.  This      *
 * allows for a more robust and efficient mesh adjustment.  The union *
 * interfaces are determined using the Weiler-Atherton algorithm (K.  *
 * Weiler and P. Atherton, Hidden Surface Removal Using Polygon Area  *
 * Sorting, Proceedings of SIGGRAPH '77, pp 214--222, 1977).  All     *
 * interfaces are treated as polygons (linear edges) regardless of    *
 * how the interface spline is defined.  Extension to higher-order    *
 * splines should be pursued in the future.                           *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Construct_Interface_Union_List(void) {

  // Construct the interface union list from the component list.
  Interface_Union_List.Construct_Union_List(Interface_Component_List);

  // Interface union list successfully constructed.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Output_Interface_Component_List_Tecplot --   *
 *                                                                    *
 * Writes the interface component list to output data files in a      *
 * format suitable for plotting with TECPLOT.                         *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Output_Interface_Component_List_Tecplot(void) {

  int n, i_output_title;
  char prefix[256], extension[256], output_file_name[256];
  char *output_file_name_ptr;
  ofstream output_file;    

  // Determine prefix of output data file names.
  n = 0;
  while (1) {
    if (IP->Output_File_Name[n] == ' ' || IP->Output_File_Name[n] == '.') break;
    prefix[n] = IP->Output_File_Name[n]; n++;
    if (n > strlen(IP->Output_File_Name)) break;
  }
  prefix[n] = '\0';
  strcat(prefix,"_interface_component_list_cpu");
  
  // Determine output data file name for this processor.
  sprintf(extension,"%.6d",Local_Solution_Block_List->ThisCPU);
  strcat(extension,".dat");
  strcpy(output_file_name,prefix);
  strcat(output_file_name,extension);
  output_file_name_ptr = output_file_name;

  // Open the output data file.
  output_file.open(output_file_name_ptr,ios::out);
  if (output_file.bad()) return 1;

  // Write the solution data for each solution block.
  i_output_title = 1;
  output_file << setprecision(14);
  if (Interface_Component_List.Ni > 0) {
    if (output_file_name)
      output_file << "TITLE = \"" << CFFC_Name() << ": Interface Component List, "
		  << "\"" << "\n"
		  << "VARIABLES = \"x\" \\ \n"
		  << "\"y\" \\ \n"
		  << "\"tp\" \\ \n"
		  << "\"bc\" \\ \n"
		  << "\"Ftype\" \\ \n"
		  << "\"Fx\" \\ \n"
		  << "\"Fy\" \\ \n"
		  << "\"|F|\" \\ \n"
		  << "\"nhat_x\" \\ \n"
		  << "\"nhat_y\" \\ \n";
    for (int ni = 1; ni <= Interface_Component_List.Ni; ni++) {    
      output_file << "ZONE T = \"" << Local_Solution_Block_List->ThisCPU << ni << "\" \\ \n"
		  << "I = " << Interface_Component_List[ni].Spline.np << " \\ \n"
		  << "J = " << 1 << " \\ \n" << "F = POINT \n";
      for (int np = 0; np < Interface_Component_List[ni].Spline.np; np++)
	output_file << Interface_Component_List[ni].Spline.Xp[np]
		    << " " << Interface_Component_List[ni].Spline.tp[np]
		    << " " << Interface_Component_List[ni].Spline.bc[np]
		    << " " << Interface_Component_List[ni].F_Type[np]
		    << Interface_Component_List[ni].F[np]
		    << " " << abs(Interface_Component_List[ni].F[np])
		    << Interface_Component_List[ni].normal(np)
		    << endl;
    }
  }
  
  // Close the output data file.
  output_file.close();
  
  // Writing of output data files complete.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Output_Union_Component_List_Tecplot --       *
 *                                                                    *
 * Writes the interface union list to output data files in a format   *
 * format suitable for plotting with TECPLOT.                         *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Output_Interface_Union_List_Tecplot(void) {

  int n, i_output_title;
  char prefix[256], extension[256], output_file_name[256];
  char *output_file_name_ptr;
  ofstream output_file;    

  // Determine prefix of output data file names.
  n = 0;
  while (1) {
    if (IP->Output_File_Name[n] == ' ' || IP->Output_File_Name[n] == '.') break;
    prefix[n] = IP->Output_File_Name[n]; n++;
    if (n > strlen(IP->Output_File_Name)) break;
  }
  prefix[n] = '\0';
  strcat(prefix,"_interface_union_list_cpu");
  
  // Determine output data file name for this processor.
  sprintf(extension,"%.6d",Local_Solution_Block_List->ThisCPU);
  strcat(extension,".dat");
  strcpy(output_file_name,prefix);
  strcat(output_file_name,extension);
  output_file_name_ptr = output_file_name;

  // Open the output data file.
  output_file.open(output_file_name_ptr,ios::out);
  if (output_file.bad()) return 1;

  // Write the solution data for each solution block.
  i_output_title = 1;
  output_file << setprecision(14);
  if (Interface_Union_List.Ni > 0) {
    if (output_file_name)
      output_file << "TITLE = \"" << CFFC_Name() << ": Interface Union List, "
		  << "\"" << "\n"
		  << "VARIABLES = \"x\" \\ \n"
		  << "\"y\" \\ \n"
		  << "\"tp\" \\ \n"
		  << "\"bc\" \\ \n"
		  << "\"Ftype\" \\ \n"
		  << "\"Fx\" \\ \n"
		  << "\"Fy\" \\ \n"
		  << "\"|F|\" \\ \n"
		  << "\"Fn\" \\ \n";
    for (int ni = 1; ni <= Interface_Union_List.Ni; ni++) {    
      output_file << "ZONE T = \"" << Local_Solution_Block_List->ThisCPU << ni << "\" \\ \n"
		  << "I = " << Interface_Union_List[ni].Spline.np << " \\ \n"
		  << "J = " << 1 << " \\ \n" << "F = POINT \n";
      for (int np = 0; np < Interface_Union_List[ni].Spline.np; np++)
	output_file << Interface_Union_List[ni].Spline.Xp[np]
		    << " " << Interface_Union_List[ni].Spline.tp[np]
		    << " " << Interface_Union_List[ni].Spline.bc[np]
		    << " " << Interface_Union_List[ni].F_Type[np]
		    << Interface_Union_List[ni].F[np]
		    << " " << abs(Interface_Union_List[ni].F[np])
 		    << " " << Interface_Union_List[ni].Fn(np)
		    << endl;
    }
  }
  
  // Close the output data file.
  output_file.close();
  
  // Writing of output data files complete.
  return 0;

}

/**********************************************************************
 **********************************************************************
 ** Interface construction and manipulation routines.                **
 **********************************************************************
 **********************************************************************/

/**********************************************************************
 * EmbeddedBoundaries2D::Initialize_Adjustment_Grids --               *
 *                                                                    *
 * Create and initialize the additional copies of the grids (old and  *
 * adjusted), the adjusted mesh data blocks, and the mesh adjustment  *
 * data blocks.                                                       *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
void EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Initialize_Adjustment_Grids(void) {
  for (int nb = 0; nb < IP->Number_of_Blocks_Per_Processor; nb++) {
    if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      // Copy the initial grid information.
      Copy_Quad_Block(AGrid[nb],Local_SolnBlk[nb].Grid);
      Copy_Quad_Block(OGrid[nb],Local_SolnBlk[nb].Grid);
      // Create and initialize the extra mesh data blocks.
      Mesh[nb].allocate(Local_SolnBlk[nb].Grid.NCi,Local_SolnBlk[nb].Grid.NCj);
      AMesh[nb].allocate(Local_SolnBlk[nb].Grid.NCi,Local_SolnBlk[nb].Grid.NCj);
      OMesh[nb].allocate(Local_SolnBlk[nb].Grid.NCi,Local_SolnBlk[nb].Grid.NCj);
      // Create and initialize the mesh adjustment data blocks.
      Adjustment_Data[nb].allocate(Local_SolnBlk[nb].NCi,Local_SolnBlk[nb].NCj,
				   Local_SolnBlk[nb].Nghost,Interface_Union_List.Ni);
    }
  }
}

/**********************************************************************
 * EmbeddedBoundaries2D::Store_Adjusted_Mesh --                       *
 *                                                                    *
 * Store the adjusted grid for later use.  The old adjusted mesh is   *
 * used when remapping the solution onto the new adjusted mesh.       *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
void EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Store_Adjusted_Mesh(void) {
  for (int nb = 0; nb < IP->Number_of_Blocks_Per_Processor; nb++) {
    if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      Store_Adjusted_Mesh(nb);
    }
  }
}

template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
void EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Store_Adjusted_Mesh(const int &nb) {
  Copy_Quad_Block(AGrid[nb],Local_SolnBlk[nb].Grid);
  AMesh[nb].Copy_Quad_Block(Mesh[nb]);
}

/**********************************************************************
 * EmbeddedBoundaries2D::Store_Unadjusted_Mesh --                     *
 *                                                                    *
 * Store the unadjusted mesh for later use.  The unadjused mesh is    *
 * used before conducting the mesh adjustment for the next time-step  *
 * for moving embedded boundaries and before AMR.                     *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
void EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Store_Unadjusted_Mesh(void) {
  for (int nb = 0; nb < IP->Number_of_Blocks_Per_Processor; nb++) {
    if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      Store_Unadjusted_Mesh(nb);
    }
  }
}

template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
void EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Store_Unadjusted_Mesh(const int &nb) {
  Copy_Quad_Block(OGrid[nb],Local_SolnBlk[nb].Grid);
  OMesh[nb].Copy_Quad_Block(Mesh[nb]);
}

/**********************************************************************
 * EmbeddedBoundaries2D::Mesh_Unadjustment --                         *
 *                                                                    *
 * Unadjust the mesh.  Copy the old (unadjusted) grid into the main   *
 * grid before permorming the mesh adjustment algorithm and before    *
 * AMR.                                                               *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
void EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Mesh_Unadjustment(void) {
  for (int nb = 0; nb < IP->Number_of_Blocks_Per_Processor; nb++) {
    if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      Mesh_Unadjustment(nb);
    }
  }
  int error_flag = Send_All_Messages(Local_SolnBlk,
				     *Local_Solution_Block_List,
				     NUM_COMP_VECTOR2D,
				     ON);
  //CFFC_Broadcast_MPI(&error_flag,1);
  //if (error_flag) return error_flag;
}

template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
void EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Mesh_Unadjustment(const int &nb) {
  Copy_Quad_Block(Local_SolnBlk[nb].Grid,OGrid[nb]);
  Mesh[nb].Copy_Quad_Block(OMesh[nb]);
  Adjustment_Data[nb].Number_of_Inactive_Cells = 0;
  Adjustment_Data[nb].Reset_Interface_Present();
}

/**********************************************************************
 * EmbeddedBoundaries2D::Mesh_Adjustment --                           *
 *                                                                    *
 * This routine performs the mesh adjustment on an array of solution  *
 * blocks given the interface component and union lists.  The mesh    *
 * adjustment is performed in eight steps:                            *
 *   1. Pre-mesh painting,                                            *
 *   2. Adjustment to sharp corners,                                  *
 *   3. First mesh adjustment (adjust nearest nodes),                 *
 *   4. Second mesh adjustment (adjust nearest nodes),                *
 *   5. Third mesh adjustment (adjust for split cells),               *
 *   6. Modify boundary splines,                                      *
 *   7. Finalization, and                                             *
 *   8. Post-mesh painting.                                           *
 * Note that steps 1-6 are conducted on single blocks without any     *
 * communication with the other solution blocks.  Step 6. is followed *
 * by message passing to ensure that the ghost cells are consistent   *
 * with the adjustment done in neighbour blocks.  The mesh adjustment *
 * statistics are computed if required.                               *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Mesh_Adjustment(const int &copy_flag,
		const int &statistics_flag) {

  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

  // Exit immediately if no interface components have been specified.
  if (!Interface_Component_List.Ni) return 0;

  int error_flag = 0;

  // Conduct point-in-interface tests to eliminate most of the cells 
  // that are completely located inside or outside the interface(s)
  // from the mesh adjustment algorithm.
  for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
    if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      error_flag = Pre_Mesh_Adjustment_Painting(nb);
      if (error_flag) return error_flag;
    }
  }

  // Conduct mesh adjustment algorithm on each solution block.
  for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
    if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      // Initialize the adjustment data.
      if (!error_flag) Adjustment_Data[nb].Initialize();
      // Adjust to sharp corners.
      if (!error_flag) error_flag = Mesh_Adjustment_Sharp_Corners(nb);
      // Perform first adjustment pass.
      if (!error_flag) error_flag = Mesh_Adjustment_First(nb);
      // Perform second adjustment pass.
      if (!error_flag) error_flag = Mesh_Adjustment_Second(nb);
      // Perform third adjustment pass (split cells).
      if (!error_flag) error_flag = Mesh_Adjustment_Third(nb);
      // Perform third adjustment pass (split cells).
      if (!error_flag) error_flag = Mesh_Adjustment_Modify_Boundary_Splines(nb);
    }
    if (error_flag) {
      cout << endl
	   << " " << error_flag
	   << "::" << Global_Solution_Block_List->ThisCPU
	   << "::" << nb
	   << "::" << Local_Solution_Block_List->Block[nb].gblknum;
      cout.flush();
      break;
    }
    //if (error_flag) break;
  }
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) {
    error_flag = Output_Nodes_Tecplot();
    return 1;//error_flag;
  }

  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

  // Update solution information shared between neighbouring blocks.
  error_flag = Send_All_Messages(Local_SolnBlk,
				 *Local_Solution_Block_List,
 				 NUM_COMP_VECTOR2D,
 				 ON);
  CFFC_Broadcast_MPI(&error_flag,1);
  if (error_flag) return error_flag;

  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

  // Finalize the mesh adjustment procedure.
  for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
    if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      error_flag = Mesh_Adjustment_Finalize(nb);
      if (error_flag) return error_flag;
    }
  }

  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

  // Update solution information shared between neighbouring blocks.
  error_flag = Send_All_Messages(Local_SolnBlk,
				 *Local_Solution_Block_List,
 				 NUM_COMP_VECTOR2D,
 				 ON);
  if (!error_flag) error_flag = Send_All_Messages(Local_SolnBlk,
 						  *Local_Solution_Block_List,
 						  Local_SolnBlk[0].NumVar(),
 						  OFF);
  if (error_flag) return error_flag;

  // Conduct point-in-interface tests to determine the cell status of 
  // all unknown cells. 
  for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
    if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      error_flag = Post_Mesh_Adjustment_Painting(nb);
      if (error_flag) return error_flag;
    }
  }

  // Store the adjusted mesh for future reference.
  if (copy_flag) Store_Adjusted_Mesh();

  // Output mesh adjustment statistics if desired.
  if (statistics_flag) Mesh_Adjustment_Statistics();

  // Mesh adjustment scheme successfully performed.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Mesh_Adjustment_Sharp_Corners --             *
 *                                                                    *
 * This routine adjusts the specified solution block to account for   *
 * sharp corners defined in the interfaces (such as the trailing edge *
 * of an aerofoil).  The node nearest to the sharp point is moved to  *
 * that point and is re-tagged as 'aligned.'                          *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Mesh_Adjustment_Sharp_Corners(const int &nb) {

  // Do not adjust if no interfaces have been defined.
  if (!Interface_Union_List.Ni) return 0;

  // Do not perform adjustment algorithm on the current solution block 
  // if there are no inactive or unknowns cells.
  if (!Adjustment_Data[nb].Number_of_Inactive_Cells) return 0;

  // Declare local variables.
  int ii, jj, iii, jjj;
  double distance;
  Polygon P;

  // Locate sharp corners and adjust mesh immediately according to 
  // their location.
  for (int ni = 1; ni <= Interface_Union_List.Ni; ni++) {
    if (Adjustment_Data[nb].Interface_Present[ni]) {
      for (int np = 0; np < Interface_Union_List[ni].Spline.np; np++) {
	if (Interface_Union_List[ni].Spline.tp[np] == SPLINE2D_POINT_SHARP_CORNER) {
 	  distance = MILLION;
 	  ii = 0; jj = 0;
 	  // Locate the nearest node.
 	  for (int j = Local_SolnBlk[nb].JCl-Local_SolnBlk[nb].Nghost; j <= Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost; j++) {
 	    for (int i = Local_SolnBlk[nb].ICl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].ICu+Local_SolnBlk[nb].Nghost; i++) {
	      if (Cell_Status_Test(nb,i,j)) {
		P.convert(Local_SolnBlk[nb].Grid.nodeSW(i,j).X,Local_SolnBlk[nb].Grid.nodeSE(i,j).X,
			  Local_SolnBlk[nb].Grid.nodeNE(i,j).X,Local_SolnBlk[nb].Grid.nodeNW(i,j).X);
		if (P.point_in_polygon(Interface_Union_List[ni].Spline.Xp[np])) {
		  if (abs(Interface_Union_List[ni].Spline.Xp[np] - Local_SolnBlk[nb].Grid.Node[i  ][j  ].X) < distance) {
		    distance = abs(Interface_Union_List[ni].Spline.Xp[np] - Local_SolnBlk[nb].Grid.Node[i  ][j  ].X);
		    ii = i; jj = j; iii = 0; jjj = 0;
		  }
		  if (abs(Interface_Union_List[ni].Spline.Xp[np] - Local_SolnBlk[nb].Grid.Node[i+1][j  ].X) < distance) {
		    distance = abs(Interface_Union_List[ni].Spline.Xp[np] - Local_SolnBlk[nb].Grid.Node[i+1][j  ].X);
		    ii = i; jj = j; iii = 1; jjj = 0;
		  }
		  if (abs(Interface_Union_List[ni].Spline.Xp[np] - Local_SolnBlk[nb].Grid.Node[i  ][j+1].X) < distance) {
		    distance = abs(Interface_Union_List[ni].Spline.Xp[np] - Local_SolnBlk[nb].Grid.Node[i  ][j+1].X);
		    ii = i; jj = j; iii = 0; jjj = 1;
		  }
		  if (abs(Interface_Union_List[ni].Spline.Xp[np] - Local_SolnBlk[nb].Grid.Node[i+1][j+1].X) < distance) {
		    distance = abs(Interface_Union_List[ni].Spline.Xp[np] - Local_SolnBlk[nb].Grid.Node[i+1][j+1].X);
 		  ii = i; jj = j; iii = 1; jjj = 1;
		  }
		}
		P.deallocate();
  	      }
 	    }
 	  }
 	  if (distance < MILLION) {
 	    Local_SolnBlk[nb].Grid.Node[ii+iii][jj+jjj].X = Interface_Union_List[ni].Spline.Xp[np];
 	    Mesh[nb].node_status[ii+iii][jj+jjj] = NODE_STATUS_ALIGNED;
 	  }
   	}
      }
    }
  }

#ifdef _EB_PARALLEL_DEBUG_
//   dout << endl << "ADJUSTED TO SHARP CORNERS::Block Number " << Local_Solution_Block_List->Block[nb].gblknum << endl;
//   for (int j = Local_SolnBlk[nb].Grid.JNu+Local_SolnBlk[nb].Nghost; j >= Local_SolnBlk[nb].Grid.JNl-Local_SolnBlk[nb].Nghost; j--) {
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       dout << Mesh[nb].node_status[i][j]+1;
//     }
//     dout << endl;
//   }
//   dout << endl;
#endif

  // Mesh adjustment to sharp corners was successful.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Mesh_Adjustment_First --                     *
 *                                                                    *
 * This routine performs the first mesh adjustment on the specified   *
 * solution block.  Here the nodes nearest to the interface are moved *
 * onto the interface and are re-tagged as aligned.  Note that only   *
 * unknown nodes of unknown cells are candidates for adjustment.      *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Mesh_Adjustment_First(const int &nb) {

  double rel_tol;
  int count;

  // Do not adjust if no interfaces have been defined.
  if (!Interface_Union_List.Ni) return 0;

  // Do not perform adjustment algorithm on the current solution block 
  // if there are no inactive or unknowns cells.
  if (!Adjustment_Data[nb].Number_of_Inactive_Cells) return 0;

  // Declare local variables.
  int intersect_flag;
  Vector2D Xm1, Xm2, Xs1, Xs2, Xp;

  // First adjustment: Determine intersection flags and intersection points.
  for (int j = Local_SolnBlk[nb].JCl-Local_SolnBlk[nb].Nghost; j <= Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost; j++) {
    for (int i = Local_SolnBlk[nb].ICl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].ICu+Local_SolnBlk[nb].Nghost; i++) {

      // Continue only if the current or neighbour cell status is
      // unknown since there is no need to adjust otherwise.
      if (Cell_Status_Test(nb,i,j)) {

	// For each interface in the union list.
	for (int ni = 1; ni <= Interface_Union_List.Ni; ni++) {

	  // Determine if the interface is present in the current block.
	  if (Adjustment_Data[nb].Interface_Present[ni]) {

	    // For each (linear) edge of the interface.
	    for (int np = 0; np < Interface_Union_List[ni].Spline.np-1; np++) {

	      // Store the (linear) edge of the interface.
	      Xs1 = Interface_Union_List[ni].Spline.Xp[np  ];
	      Xs2 = Interface_Union_List[ni].Spline.Xp[np+1];

	      // Define the intersection flags and the intersection points 
	      // for the SOUH-WEST node if not already aligned.
	      if (Mesh[nb].node_status[i][j] != NODE_STATUS_ALIGNED) {
		Xm1 = Local_SolnBlk[nb].Grid.Node[i][j].X;
		// NORTH edge.
		if (!((j == Local_SolnBlk[nb].Grid.JNu && Local_SolnBlk[nb].Grid.BCtypeN[i] != BC_NONE) ||
		      (j == Local_SolnBlk[nb].Grid.JNl && Local_SolnBlk[nb].Grid.BCtypeS[i] != BC_NONE))) {
		  if ((j+1 == Local_SolnBlk[nb].Grid.JNu && Local_SolnBlk[nb].Grid.BCtypeN[i] != BC_NONE) ||
		      (j+1 == Local_SolnBlk[nb].Grid.JNl && Local_SolnBlk[nb].Grid.BCtypeS[i] != BC_NONE)) {
		    Xm2 = Local_SolnBlk[nb].Grid.Node[i][j+1].X;
		  } else {
		    Xm2 = HALF*(Xm1 + Local_SolnBlk[nb].Grid.Node[i][j+1].X);
		  }
		  intersect_flag = Line_Intersection(Xs1,Xs2,Xm1,Xm2,Xp,TOLER*abs(Xm2-Xm1));
		  if (intersect_flag && abs(Xp-Xm1) < abs(Adjustment_Data[nb].Xnorth[i][j]-Xm1)) {
		    Adjustment_Data[nb].intersect_north[i][j] = intersect_flag;
		    Adjustment_Data[nb].Xnorth[i][j] = Xp;
		  }
		}
 	      	// SOUTH edge.
 		if (j > Local_SolnBlk[nb].JCl-Local_SolnBlk[nb].Nghost) {
 		  if (!((j == Local_SolnBlk[nb].Grid.JNu && Local_SolnBlk[nb].Grid.BCtypeN[i] != BC_NONE) ||
			(j == Local_SolnBlk[nb].Grid.JNl && Local_SolnBlk[nb].Grid.BCtypeS[i] != BC_NONE))) {
		    if ((j-1 == Local_SolnBlk[nb].Grid.JNu && Local_SolnBlk[nb].Grid.BCtypeN[i] != BC_NONE) ||
			(j-1 == Local_SolnBlk[nb].Grid.JNl && Local_SolnBlk[nb].Grid.BCtypeS[i] != BC_NONE)) {
		      Xm2 = Local_SolnBlk[nb].Grid.Node[i][j-1].X;
		    } else {
		      Xm2 = HALF*(Xm1 + Local_SolnBlk[nb].Grid.Node[i][j-1].X);
		    }
		    intersect_flag = Line_Intersection(Xs1,Xs2,Xm1,Xm2,Xp,TOLER*abs(Xm2-Xm1));
		    if (intersect_flag && abs(Xp-Xm1) < abs(Adjustment_Data[nb].Xsouth[i][j]-Xm1)) {
		      Adjustment_Data[nb].intersect_south[i][j] = intersect_flag;
		      Adjustment_Data[nb].Xsouth[i][j] = Xp;
		    }
		  }
		} else {
		  Xm2 = Xm1 - HALF*(Local_SolnBlk[nb].Grid.Node[i][j+1].X - Xm1);
		  intersect_flag = Line_Intersection(Xs1,Xs2,Xm1,Xm2,Xp,TOLER*abs(Xm2-Xm1));
		  if (intersect_flag && abs(Xp-Xm1) < abs(Adjustment_Data[nb].Xsouth[i][j]-Xm1)) {
		    Adjustment_Data[nb].intersect_south[i][j] = intersect_flag;
		    Adjustment_Data[nb].Xsouth[i][j] = Xp;
		  }
		}
		// EAST edge.
		if (!((i == Local_SolnBlk[nb].Grid.INu && Local_SolnBlk[nb].Grid.BCtypeE[j] != BC_NONE) ||
		      (i == Local_SolnBlk[nb].Grid.INl && Local_SolnBlk[nb].Grid.BCtypeW[j] != BC_NONE))) {
		  if ((i+1 == Local_SolnBlk[nb].Grid.INu && Local_SolnBlk[nb].Grid.BCtypeE[j] != BC_NONE) ||
		      (i+1 == Local_SolnBlk[nb].Grid.INl && Local_SolnBlk[nb].Grid.BCtypeW[j] != BC_NONE)) {
		    Xm2 = Local_SolnBlk[nb].Grid.Node[i+1][j].X;
		  } else {
		    Xm2 = HALF*(Xm1 + Local_SolnBlk[nb].Grid.Node[i+1][j].X);
		  }
		  intersect_flag = Line_Intersection(Xs1,Xs2,Xm1,Xm2,Xp,TOLER*abs(Xm2-Xm1));
		  if (intersect_flag && abs(Xp-Xm1) < abs(Adjustment_Data[nb].Xeast[i][j]-Xm1)) {
		    Adjustment_Data[nb].intersect_east[i][j] = intersect_flag;
		    Adjustment_Data[nb].Xeast[i][j] = Xp;
		  }
		}
		// WEST edge.
		if (i > Local_SolnBlk[nb].ICl-Local_SolnBlk[nb].Nghost) {
		  if (!((i == Local_SolnBlk[nb].Grid.INu && Local_SolnBlk[nb].Grid.BCtypeE[j] != BC_NONE) ||
			(i == Local_SolnBlk[nb].Grid.INl && Local_SolnBlk[nb].Grid.BCtypeW[j] != BC_NONE))) {
		    if ((i-1 == Local_SolnBlk[nb].Grid.INu && Local_SolnBlk[nb].Grid.BCtypeE[j] != BC_NONE) ||
			(i-1 == Local_SolnBlk[nb].Grid.INl && Local_SolnBlk[nb].Grid.BCtypeW[j] != BC_NONE)) {
		      Xm2 = Local_SolnBlk[nb].Grid.Node[i-1][j].X;
		    } else {
		      Xm2 = HALF*(Xm1 + Local_SolnBlk[nb].Grid.Node[i-1][j].X);
		    }
		    intersect_flag = Line_Intersection(Xs1,Xs2,Xm1,Xm2,Xp,TOLER*abs(Xm2-Xm1));
		    if (intersect_flag && abs(Xp-Xm1) < abs(Adjustment_Data[nb].Xwest[i][j]-Xm1)) {
		      Adjustment_Data[nb].intersect_west[i][j] = intersect_flag;
		      Adjustment_Data[nb].Xwest[i][j] = Xp;
		    }
		  }
		}
 	      }

	      // Define the intersection flags and the intersection points 
	      // for the SOUTH-EAST node if necessary.
	      if (i == Local_SolnBlk[nb].ICu+Local_SolnBlk[nb].Nghost &&
		  Mesh[nb].node_status[i+1][j] != NODE_STATUS_ALIGNED) {
		Xm1 = Local_SolnBlk[nb].Grid.Node[i+1][j].X;
		// NORTH edge.
		if (!((j == Local_SolnBlk[nb].Grid.JNu && Local_SolnBlk[nb].Grid.BCtypeN[i] != BC_NONE) ||
		      (j == Local_SolnBlk[nb].Grid.JNl && Local_SolnBlk[nb].Grid.BCtypeS[i] != BC_NONE))) {
		  if ((j+1 == Local_SolnBlk[nb].Grid.JNu && Local_SolnBlk[nb].Grid.BCtypeN[i] != BC_NONE) ||
		      (j+1 == Local_SolnBlk[nb].Grid.JNl && Local_SolnBlk[nb].Grid.BCtypeS[i] != BC_NONE)) {
		    Xm2 = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
		  } else {
		    Xm2 = HALF*(Xm1 + Local_SolnBlk[nb].Grid.Node[i+1][j+1].X);
		  }
		  intersect_flag = Line_Intersection(Xs1,Xs2,Xm1,Xm2,Xp,TOLER*abs(Xm2-Xm1));
		  if (intersect_flag && abs(Xp-Xm1) < abs(Adjustment_Data[nb].Xnorth[i+1][j]-Xm1)) {
		    Adjustment_Data[nb].intersect_north[i+1][j] = intersect_flag;
		    Adjustment_Data[nb].Xnorth[i+1][j] = Xp;
		  }
		}
		// SOUTH edge.
		if (j > Local_SolnBlk[nb].JCl-Local_SolnBlk[nb].Nghost) {
 		  if (!((j == Local_SolnBlk[nb].Grid.JNu && Local_SolnBlk[nb].Grid.BCtypeN[i] != BC_NONE) ||
			(j == Local_SolnBlk[nb].Grid.JNl && Local_SolnBlk[nb].Grid.BCtypeS[i] != BC_NONE))) {
		    if ((j-1 == Local_SolnBlk[nb].Grid.JNu && Local_SolnBlk[nb].Grid.BCtypeN[i] != BC_NONE) ||
			(j-1 == Local_SolnBlk[nb].Grid.JNl && Local_SolnBlk[nb].Grid.BCtypeS[i] != BC_NONE)) {
		      Xm2 = Local_SolnBlk[nb].Grid.Node[i+1][j-1].X;
		    } else {
		      Xm2 = HALF*(Xm1 + Local_SolnBlk[nb].Grid.Node[i+1][j-1].X);
		    }
		    intersect_flag = Line_Intersection(Xs1,Xs2,Xm1,Xm2,Xp,TOLER*abs(Xm2-Xm1));
		    if (intersect_flag && abs(Xp-Xm1) < abs(Adjustment_Data[nb].Xsouth[i+1][j]-Xm1)) {
		      Adjustment_Data[nb].intersect_south[i+1][j] = intersect_flag;
		      Adjustment_Data[nb].Xsouth[i+1][j] = Xp;
		    }
		  }
		}
		// WEST edge.
		if (!((i+1 == Local_SolnBlk[nb].Grid.INu && Local_SolnBlk[nb].Grid.BCtypeE[j] != BC_NONE) ||
		      (i+1 == Local_SolnBlk[nb].Grid.INl && Local_SolnBlk[nb].Grid.BCtypeW[j] != BC_NONE))) {
		  if ((i == Local_SolnBlk[nb].Grid.INu && Local_SolnBlk[nb].Grid.BCtypeE[j] != BC_NONE) ||
		      (i == Local_SolnBlk[nb].Grid.INl && Local_SolnBlk[nb].Grid.BCtypeW[j] != BC_NONE)) {
		    Xm2 = Local_SolnBlk[nb].Grid.Node[i][j].X;
		  } else {
		    Xm2 = HALF*(Xm1 + Local_SolnBlk[nb].Grid.Node[i][j].X);
		  }
		  intersect_flag = Line_Intersection(Xs1,Xs2,Xm1,Xm2,Xp,TOLER*abs(Xm2-Xm1));
		  if (intersect_flag && abs(Xp-Xm1) < abs(Adjustment_Data[nb].Xwest[i+1][j]-Xm1)) {
		    Adjustment_Data[nb].intersect_west[i+1][j] = intersect_flag;
		    Adjustment_Data[nb].Xwest[i+1][j] = Xp;
		  }
		}
	      }

	      // Define the intersection flags and the intersection points
	      // for the NORTH-WEST node if necessary.
	      if (j == Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost &&
		  Mesh[nb].node_status[i][j+1] != NODE_STATUS_ALIGNED) {
		Xm1 = Local_SolnBlk[nb].Grid.Node[i][j+1].X;
		// NORTH edge.
		Xm2 = Xm1 + HALF*(Xm1 - Local_SolnBlk[nb].Grid.Node[i][j].X);
		intersect_flag = Line_Intersection(Xs1,Xs2,Xm1,Xm2,Xp,TOLER*abs(Xm2-Xm1));
		if (intersect_flag && abs(Xp-Xm1) < abs(Adjustment_Data[nb].Xnorth[i][j+1]-Xm1)) {
		  Adjustment_Data[nb].intersect_north[i][j+1] = intersect_flag;
		  Adjustment_Data[nb].Xnorth[i][j+1] = Xp;
		}
		// SOUTH edge.
		if (!((j+1 == Local_SolnBlk[nb].Grid.JNu && Local_SolnBlk[nb].Grid.BCtypeN[i] != BC_NONE) ||
		      (j+1 == Local_SolnBlk[nb].Grid.JNl && Local_SolnBlk[nb].Grid.BCtypeS[i] != BC_NONE))) {
		  if ((j == Local_SolnBlk[nb].Grid.JNu && Local_SolnBlk[nb].Grid.BCtypeN[i] != BC_NONE) ||
		      (j == Local_SolnBlk[nb].Grid.JNl && Local_SolnBlk[nb].Grid.BCtypeS[i] != BC_NONE)) {
		    Xm2 = Local_SolnBlk[nb].Grid.Node[i][j].X;
		  } else {
		    Xm2 = HALF*(Xm1 + Local_SolnBlk[nb].Grid.Node[i][j].X);
		  }
		  intersect_flag = Line_Intersection(Xs1,Xs2,Xm1,Xm2,Xp,TOLER*abs(Xm2-Xm1));
		  if (intersect_flag && abs(Xp-Xm1) < abs(Adjustment_Data[nb].Xsouth[i][j+1]-Xm1)) {
		    Adjustment_Data[nb].intersect_south[i][j+1] = intersect_flag;
		    Adjustment_Data[nb].Xsouth[i][j+1] = Xp;
		  }
		}
		// EAST edge.
		if (!((i == Local_SolnBlk[nb].Grid.INu && Local_SolnBlk[nb].Grid.BCtypeE[j] != BC_NONE) ||
		      (i == Local_SolnBlk[nb].Grid.INl && Local_SolnBlk[nb].Grid.BCtypeW[j] != BC_NONE))) {
		  if ((i+1 == Local_SolnBlk[nb].Grid.INu && Local_SolnBlk[nb].Grid.BCtypeE[j] != BC_NONE) ||
		      (i+1 == Local_SolnBlk[nb].Grid.INl && Local_SolnBlk[nb].Grid.BCtypeW[j] != BC_NONE)) {
		    Xm2 = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
		  } else {
		    Xm2 = HALF*(Xm1 + Local_SolnBlk[nb].Grid.Node[i+1][j+1].X);
		  }
		  intersect_flag = Line_Intersection(Xs1,Xs2,Xm1,Xm2,Xp,TOLER*abs(Xm2-Xm1));
		  if (intersect_flag && abs(Xp-Xm1) < abs(Adjustment_Data[nb].Xeast[i][j+1]-Xm1)) {
		    Adjustment_Data[nb].intersect_east[i][j+1] = intersect_flag;
		    Adjustment_Data[nb].Xeast[i][j+1] = Xp;
		  }
		}
		// WEST edge.
		if (i > Local_SolnBlk[nb].ICl-Local_SolnBlk[nb].Nghost) {
		  if (!((i == Local_SolnBlk[nb].Grid.INu && Local_SolnBlk[nb].Grid.BCtypeE[j] != BC_NONE) ||
			(i == Local_SolnBlk[nb].Grid.INl && Local_SolnBlk[nb].Grid.BCtypeW[j] != BC_NONE))) {
		    if ((i-1 == Local_SolnBlk[nb].Grid.INu && Local_SolnBlk[nb].Grid.BCtypeE[j] != BC_NONE) ||
			(i-1 == Local_SolnBlk[nb].Grid.INl && Local_SolnBlk[nb].Grid.BCtypeW[j] != BC_NONE)) {
		      Xm2 = Local_SolnBlk[nb].Grid.Node[i-1][j+1].X;
		    } else {
		      Xm2 = HALF*(Xm1 + Local_SolnBlk[nb].Grid.Node[i-1][j+1].X);
		    }
		    intersect_flag = Line_Intersection(Xs1,Xs2,Xm1,Xm2,Xp,TOLER*abs(Xm2-Xm1));
		    if (intersect_flag && abs(Xp-Xm1) < abs(Adjustment_Data[nb].Xwest[i][j+1]-Xm1)) {
		      Adjustment_Data[nb].intersect_west[i][j+1] = intersect_flag;
		      Adjustment_Data[nb].Xwest[i][j+1] = Xp;
		    }
		  }
		}
	      }

	      // Define the intersection flags and the intersection points
	      // for the NORTH-EAST node if necessary.
	      if (i == Local_SolnBlk[nb].ICu+Local_SolnBlk[nb].Nghost &&
		  j == Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost &&
		  Mesh[nb].node_status[i+1][j+1] != NODE_STATUS_ALIGNED) {
		Xm1 = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
		// NORTH edge.
		Xm2 = Xm1 + HALF*(Xm1 - Local_SolnBlk[nb].Grid.Node[i+1][j].X);
		intersect_flag = Line_Intersection(Xs1,Xs2,Xm1,Xm2,Xp,TOLER*abs(Xm2-Xm1));
		if (intersect_flag && abs(Xp-Xm1) < abs(Adjustment_Data[nb].Xnorth[i+1][j+1]-Xm1)) {
		  Adjustment_Data[nb].intersect_north[i+1][j+1] = intersect_flag;
		  Adjustment_Data[nb].Xnorth[i+1][j+1] = Xp;
		}
		// SOUTH edge.
		if (!((j+1 == Local_SolnBlk[nb].Grid.JNu && Local_SolnBlk[nb].Grid.BCtypeN[i] != BC_NONE) ||
		      (j+1 == Local_SolnBlk[nb].Grid.JNl && Local_SolnBlk[nb].Grid.BCtypeS[i] != BC_NONE))) {
		  if ((j == Local_SolnBlk[nb].Grid.JNu && Local_SolnBlk[nb].Grid.BCtypeN[i] != BC_NONE) ||
		      (j == Local_SolnBlk[nb].Grid.JNl && Local_SolnBlk[nb].Grid.BCtypeS[i] != BC_NONE)) {
		    Xm2 = Local_SolnBlk[nb].Grid.Node[i+1][j].X;
		  } else {
		    Xm2 = HALF*(Xm1 + Local_SolnBlk[nb].Grid.Node[i+1][j].X);
		  }
		  intersect_flag = Line_Intersection(Xs1,Xs2,Xm1,Xm2,Xp,TOLER*abs(Xm2-Xm1));
		  if (intersect_flag && abs(Xp-Xm1) < abs(Adjustment_Data[nb].Xsouth[i+1][j+1]-Xm1)) {
		    Adjustment_Data[nb].intersect_south[i+1][j+1] = intersect_flag;
		    Adjustment_Data[nb].Xsouth[i+1][j+1] = Xp;
		  }
		}
		// WEST edge.
		if (!((i+1 == Local_SolnBlk[nb].Grid.INu && Local_SolnBlk[nb].Grid.BCtypeE[j] != BC_NONE) ||
		      (i+1 == Local_SolnBlk[nb].Grid.INl && Local_SolnBlk[nb].Grid.BCtypeW[j] != BC_NONE))) {
		  if ((i == Local_SolnBlk[nb].Grid.INu && Local_SolnBlk[nb].Grid.BCtypeE[j] != BC_NONE) ||
		      (i == Local_SolnBlk[nb].Grid.INl && Local_SolnBlk[nb].Grid.BCtypeW[j] != BC_NONE)) {
		    Xm2 = Local_SolnBlk[nb].Grid.Node[i][j+1].X;
		  } else {
		    Xm2 = HALF*(Xm1 + Local_SolnBlk[nb].Grid.Node[i][j+1].X);
		  }
		  intersect_flag = Line_Intersection(Xs1,Xs2,Xm1,Xm2,Xp,TOLER*abs(Xm2-Xm1));
		  if (intersect_flag && abs(Xp-Xm1) < abs(Adjustment_Data[nb].Xwest[i+1][j+1]-Xm1)) {
		    Adjustment_Data[nb].intersect_west[i+1][j+1] = intersect_flag;
		    Adjustment_Data[nb].Xwest[i+1][j+1] = Xp;
		  }
		}
 	      }

	    }
	  }
	}
      }
    }
  }

#ifdef _EB_PARALLEL_DEBUG_
//   dout << endl << "FIRST ADJUSTMENT FLAGS (before reduction)::Block Number " << Local_Solution_Block_List->Block[nb].gblknum << endl;
//   for (int j = Local_SolnBlk[nb].Grid.JNu+Local_SolnBlk[nb].Nghost; j >= Local_SolnBlk[nb].Grid.JNl-Local_SolnBlk[nb].Nghost; j--) {
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       dout << Mesh[nb].node_status[i][j]+1;
//     }
//     dout << " ";
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       dout << Adjustment_Data[nb].intersect_north[i][j];
//     }
//     dout << " ";
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       dout << Adjustment_Data[nb].intersect_south[i][j];
//     }
//     dout << " ";
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       dout << Adjustment_Data[nb].intersect_east[i][j];
//     }
//     dout << " ";
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       dout << Adjustment_Data[nb].intersect_west[i][j];
//     }
//     dout << endl;
//   }
//   dout << endl;
//   dout << setprecision(14);
//   for (int j = Local_SolnBlk[nb].Grid.JNu+Local_SolnBlk[nb].Nghost; j >= Local_SolnBlk[nb].Grid.JNl-Local_SolnBlk[nb].Nghost; j--) {
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       if (Adjustment_Data[nb].intersect_north[i][j]) dout << endl << " North intersection:" << Local_SolnBlk[nb].Grid.Node[i][j].X << Adjustment_Data[nb].Xnorth[i][j] << " " << abs(Local_SolnBlk[nb].Grid.Node[i][j].X-Adjustment_Data[nb].Xnorth[i][j]);
//       if (Adjustment_Data[nb].intersect_south[i][j]) dout << endl << " South intersection:" << Local_SolnBlk[nb].Grid.Node[i][j].X << Adjustment_Data[nb].Xsouth[i][j] << " " << abs(Local_SolnBlk[nb].Grid.Node[i][j].X-Adjustment_Data[nb].Xsouth[i][j]);
//       if (Adjustment_Data[nb].intersect_east[i][j]) dout << endl << " East intersection:" << Local_SolnBlk[nb].Grid.Node[i][j].X << Adjustment_Data[nb].Xeast[i][j] << " " << abs(Local_SolnBlk[nb].Grid.Node[i][j].X-Adjustment_Data[nb].Xeast[i][j]);
//       if (Adjustment_Data[nb].intersect_west[i][j]) dout << endl << " West intersection:" << Local_SolnBlk[nb].Grid.Node[i][j].X << Adjustment_Data[nb].Xwest[i][j] << " " << abs(Local_SolnBlk[nb].Grid.Node[i][j].X-Adjustment_Data[nb].Xwest[i][j]);
//     }
//   }
//   dout << setprecision(6);
//   dout << endl;
#endif

  // Adjustment choice reduction:
  // -> Nodes that coincide with their determined pierce-points are
  //    re-tagged as aligned and that node's adjustment flags are reset.
  // -> For nodes that intersect to the NORTH and SOUTH or the EAST and 
  //    WEST, turn off intersect flags for both competing choices and
  //    adjust the neighbours instead (important for aerofoils).
  for (int j = Local_SolnBlk[nb].Grid.JNl-Local_SolnBlk[nb].Nghost; j <= Local_SolnBlk[nb].Grid.JNu+Local_SolnBlk[nb].Nghost; j++) {
    for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
      //determine relative tolerance
      count = 0;
      rel_tol = 0.0;
      if(j > Local_SolnBlk[nb].Grid.JNl-Local_SolnBlk[nb].Nghost) {
	rel_tol += abs(Local_SolnBlk[nb].Grid.Node[i][j].X-Local_SolnBlk[nb].Grid.Node[i][j-1].X);
	++count;
      }
      if(j < Local_SolnBlk[nb].Grid.JNu+Local_SolnBlk[nb].Nghost) {
	rel_tol += abs(Local_SolnBlk[nb].Grid.Node[i][j].X-Local_SolnBlk[nb].Grid.Node[i][j+1].X);
	++count;
      }
      if(i > Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost) {
	rel_tol += abs(Local_SolnBlk[nb].Grid.Node[i][j].X-Local_SolnBlk[nb].Grid.Node[i-1][j].X);
	++count;
      }
      if(i < Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost) {
	rel_tol += abs(Local_SolnBlk[nb].Grid.Node[i][j].X-Local_SolnBlk[nb].Grid.Node[i+1][j].X);
	++count;
      }
    
      rel_tol = TOLER*rel_tol/((double)count);

      if (Mesh[nb].node_status[i][j] != NODE_STATUS_ALIGNED) {
	if ((Adjustment_Data[nb].intersect_north[i][j] && abs(Adjustment_Data[nb].Xnorth[i][j]-Local_SolnBlk[nb].Grid.Node[i][j].X) < rel_tol) ||
	    (Adjustment_Data[nb].intersect_south[i][j] && abs(Adjustment_Data[nb].Xsouth[i][j]-Local_SolnBlk[nb].Grid.Node[i][j].X) < rel_tol) ||
	    (Adjustment_Data[nb].intersect_east[i][j] && abs(Adjustment_Data[nb].Xeast[i][j]-Local_SolnBlk[nb].Grid.Node[i][j].X) < rel_tol) ||
	    (Adjustment_Data[nb].intersect_west[i][j] && abs(Adjustment_Data[nb].Xwest[i][j]-Local_SolnBlk[nb].Grid.Node[i][j].X) < rel_tol)) {
	  Mesh[nb].node_status[i][j] = NODE_STATUS_ALIGNED;
	  Adjustment_Data[nb].intersect_north[i][j] = OFF;
	  Adjustment_Data[nb].intersect_south[i][j] = OFF;
	  Adjustment_Data[nb].intersect_east[i][j] = OFF;
	  Adjustment_Data[nb].intersect_west[i][j] = OFF;
	}
	if (Adjustment_Data[nb].intersect_north[i][j] && j < Local_SolnBlk[nb].Grid.JNu+Local_SolnBlk[nb].Nghost) {
	  if (abs(Adjustment_Data[nb].Xnorth[i][j]-Local_SolnBlk[nb].Grid.Node[i][j+1].X) < rel_tol) {
	    Adjustment_Data[nb].intersect_north[i][j] = OFF;
	    Mesh[nb].node_status[i][j+1] = NODE_STATUS_ALIGNED;
	  }
	}
	if (Adjustment_Data[nb].intersect_south[i][j] && j > Local_SolnBlk[nb].Grid.JNl-Local_SolnBlk[nb].Nghost) {
	  if (abs(Adjustment_Data[nb].Xsouth[i][j]-Local_SolnBlk[nb].Grid.Node[i][j-1].X) < rel_tol) {
	    Adjustment_Data[nb].intersect_south[i][j] = OFF;
	    Mesh[nb].node_status[i][j-1] = NODE_STATUS_ALIGNED;
	  }
	}
	if (Adjustment_Data[nb].intersect_east[i][j] && i < Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost) {
	  if (abs(Adjustment_Data[nb].Xeast[i][j]-Local_SolnBlk[nb].Grid.Node[i+1][j].X) < rel_tol) {
	    Adjustment_Data[nb].intersect_east[i][j] = OFF;
	    Mesh[nb].node_status[i+1][j] = NODE_STATUS_ALIGNED;
	  }
	}
	if (Adjustment_Data[nb].intersect_west[i][j] && i > Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost) {
	  if (abs(Adjustment_Data[nb].Xwest[i][j]-Local_SolnBlk[nb].Grid.Node[i-1][j].X) < rel_tol) {
	    Adjustment_Data[nb].intersect_west[i][j] = OFF;
	    Mesh[nb].node_status[i-1][j] = NODE_STATUS_ALIGNED;
	  }
	}
      }
    }
  }
  for (int j = Local_SolnBlk[nb].Grid.JNl-Local_SolnBlk[nb].Nghost; j <= Local_SolnBlk[nb].Grid.JNu+Local_SolnBlk[nb].Nghost; j++) {
    for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
      if (Mesh[nb].node_status[i][j] != NODE_STATUS_ALIGNED) {
	if (Adjustment_Data[nb].intersect_north[i][j] && Adjustment_Data[nb].intersect_south[i][j]) {
	  if (fabs(abs(Adjustment_Data[nb].Xnorth[i][j]-Local_SolnBlk[nb].Grid.Node[i][j].X) -
		   abs(Adjustment_Data[nb].Xsouth[i][j]-Local_SolnBlk[nb].Grid.Node[i][j].X)) <
	      0.25*abs(Adjustment_Data[nb].Xnorth[i][j]-Adjustment_Data[nb].Xsouth[i][j])) {
	    Adjustment_Data[nb].intersect_north[i][j] = OFF;
	    Adjustment_Data[nb].intersect_south[i][j] = OFF;
	    Adjustment_Data[nb].intersect_south[i][j+1] = ON;
	    Adjustment_Data[nb].intersect_north[i][j-1] = ON;
	    Adjustment_Data[nb].Xsouth[i][j+1] = Adjustment_Data[nb].Xnorth[i][j];
	    Adjustment_Data[nb].Xnorth[i][j-1] = Adjustment_Data[nb].Xsouth[i][j];
	    Adjustment_Data[nb].Xnorth[i][j] = Vector2D(MILLION,MILLION);
	    Adjustment_Data[nb].Xsouth[i][j] = Vector2D(MILLION,MILLION);
	  }
	}
	if (Adjustment_Data[nb].intersect_east[i][j] && Adjustment_Data[nb].intersect_west[i][j]) {
	  if (fabs(abs(Adjustment_Data[nb].Xeast[i][j]-Local_SolnBlk[nb].Grid.Node[i][j].X) -
		   abs(Adjustment_Data[nb].Xwest[i][j]-Local_SolnBlk[nb].Grid.Node[i][j].X)) <
	      0.25*abs(Adjustment_Data[nb].Xeast[i][j]-Adjustment_Data[nb].Xwest[i][j])) {
	    Adjustment_Data[nb].intersect_east[i][j] = OFF;
	    Adjustment_Data[nb].intersect_west[i][j] = OFF;
	    Adjustment_Data[nb].intersect_east[i-1][j] = ON;
	    Adjustment_Data[nb].intersect_west[i+1][j] = ON;
	    Adjustment_Data[nb].Xeast[i-1][j] = Adjustment_Data[nb].Xwest[i][j];
	    Adjustment_Data[nb].Xwest[i+1][j] = Adjustment_Data[nb].Xeast[i][j];
	    Adjustment_Data[nb].Xeast[i][j] = Vector2D(MILLION,MILLION);
	    Adjustment_Data[nb].Xwest[i][j] = Vector2D(MILLION,MILLION);
	  }
	}
      }
    }
  }
  for (int j = Local_SolnBlk[nb].Grid.JNl-Local_SolnBlk[nb].Nghost; j <= Local_SolnBlk[nb].Grid.JNu+Local_SolnBlk[nb].Nghost; j++) {
    for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
      if (j < Local_SolnBlk[nb].Grid.JNu+Local_SolnBlk[nb].Nghost) {
	if (Adjustment_Data[nb].intersect_north[i][j] && Adjustment_Data[nb].intersect_south[i][j+1]) {
 	  if (abs(Adjustment_Data[nb].Xnorth[i][j]-Adjustment_Data[nb].Xsouth[i][j+1]) < rel_tol) {
 	    Adjustment_Data[nb].intersect_south[i][j+1] = OFF;
 	    Adjustment_Data[nb].Xsouth[i][j+1] = Vector2D(MILLION,MILLION);
 	  }
	}
      }
      if (i < Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost) {
	if (Adjustment_Data[nb].intersect_east[i][j] && Adjustment_Data[nb].intersect_west[i+1][j]) {
 	  if (abs(Adjustment_Data[nb].Xeast[i][j]-Adjustment_Data[nb].Xwest[i+1][j]) < rel_tol) {
 	    Adjustment_Data[nb].intersect_west[i+1][j] = OFF;
 	    Adjustment_Data[nb].Xwest[i+1][j] = Vector2D(MILLION,MILLION);
 	  }
	}
      }
    }
  }

#ifdef _EB_PARALLEL_DEBUG_
//   dout << endl << "FIRST ADJUSTMENT FLAGS::Block Number " << Local_Solution_Block_List->Block[nb].gblknum << endl;
//   for (int j = Local_SolnBlk[nb].Grid.JNu+Local_SolnBlk[nb].Nghost; j >= Local_SolnBlk[nb].Grid.JNl-Local_SolnBlk[nb].Nghost; j--) {
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       dout << Mesh[nb].node_status[i][j]+1;
//     }
//     dout << " ";
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       dout << Adjustment_Data[nb].intersect_north[i][j];
//     }
//     dout << " ";
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       dout << Adjustment_Data[nb].intersect_south[i][j];
//     }
//     dout << " ";
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       dout << Adjustment_Data[nb].intersect_east[i][j];
//     }
//     dout << " ";
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       dout << Adjustment_Data[nb].intersect_west[i][j];
//     }
//     dout << endl;
//   }
//   dout << endl;
//   dout << endl;
//   for (int j = Local_SolnBlk[nb].Grid.JNu+Local_SolnBlk[nb].Nghost; j >= Local_SolnBlk[nb].Grid.JNl-Local_SolnBlk[nb].Nghost; j--) {
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       if (Adjustment_Data[nb].intersect_north[i][j]) dout << "N";
//       else dout << "O";
//     }
//     dout << " ";
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       if (Adjustment_Data[nb].intersect_south[i][j]) dout << "S";
//       else dout << "O";
//     }
//     dout << " ";
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       if (Adjustment_Data[nb].intersect_east[i][j]) dout << "E";
//       else dout << "O";
//     }
//     dout << " ";
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       if (Adjustment_Data[nb].intersect_west[i][j]) dout << "W";
//       else dout << "O";
//     }
//     dout << endl;
//   }
//   dout << endl;
//   for (int j = Local_SolnBlk[nb].Grid.JNu+Local_SolnBlk[nb].Nghost; j >= Local_SolnBlk[nb].Grid.JNl-Local_SolnBlk[nb].Nghost; j--) {
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       if (Adjustment_Data[nb].intersect_north[i][j]) dout << endl << " North intersection:" << Local_SolnBlk[nb].Grid.Node[i][j].X << Adjustment_Data[nb].Xnorth[i][j];
//       if (Adjustment_Data[nb].intersect_south[i][j]) dout << endl << " South intersection:" << Local_SolnBlk[nb].Grid.Node[i][j].X << Adjustment_Data[nb].Xsouth[i][j];
//       if (Adjustment_Data[nb].intersect_east[i][j]) dout << endl << " East intersection:" << Local_SolnBlk[nb].Grid.Node[i][j].X << Adjustment_Data[nb].Xeast[i][j];
//       if (Adjustment_Data[nb].intersect_west[i][j]) dout << endl << " West intersection:" << Local_SolnBlk[nb].Grid.Node[i][j].X << Adjustment_Data[nb].Xwest[i][j];
//     }
//   }
//   dout << endl;
#endif

  // First Adjustment: Adjust cells according to the intersection flags.
  for (int j = Local_SolnBlk[nb].Grid.JNl-Local_SolnBlk[nb].Nghost; j <= Local_SolnBlk[nb].Grid.JNu+Local_SolnBlk[nb].Nghost; j++) {
    for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {

      if (Mesh[nb].node_status[i][j] != NODE_STATUS_ALIGNED) {
	// Determine the closest interesection point and adjust node accordingly.
	if (Adjustment_Data[nb].intersect_north[i][j] || Adjustment_Data[nb].intersect_south[i][j] ||
	    Adjustment_Data[nb].intersect_east[i][j] || Adjustment_Data[nb].intersect_west[i][j] ) {
	  if (abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xnorth[i][j]) <=
	      abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xsouth[i][j]) &&
	      abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xnorth[i][j]) <=
	      abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xeast[i][j]) &&
	      abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xnorth[i][j]) <=
	      abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xwest[i][j])) {
	    Local_SolnBlk[nb].Grid.Node[i][j].X = Adjustment_Data[nb].Xnorth[i][j];
	    Adjustment_Data[nb].Xnorth[i][j] = Vector2D(MILLION,MILLION);
	  } else if (abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xsouth[i][j]) <=
		     abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xnorth[i][j]) &&
		     abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xsouth[i][j]) <=
		     abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xeast[i][j]) &&
		     abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xsouth[i][j]) <=
		     abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xwest[i][j])) {
	    Local_SolnBlk[nb].Grid.Node[i][j].X = Adjustment_Data[nb].Xsouth[i][j];
	    Adjustment_Data[nb].Xsouth[i][j] = Vector2D(MILLION,MILLION);
	  } else if (abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xeast[i][j]) <=
		     abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xnorth[i][j]) &&
		     abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xeast[i][j]) <=
		     abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xsouth[i][j]) &&
		     abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xeast[i][j]) <=
		     abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xwest[i][j])) {
	    Local_SolnBlk[nb].Grid.Node[i][j].X = Adjustment_Data[nb].Xeast[i][j];
	    Adjustment_Data[nb].Xeast[i][j] = Vector2D(MILLION,MILLION);
	  } else if (abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xwest[i][j]) <=
		     abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xnorth[i][j]) &&
		     abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xwest[i][j]) <=
		     abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xsouth[i][j]) &&
		     abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xwest[i][j]) <=
		     abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xeast[i][j])) {
 	    Local_SolnBlk[nb].Grid.Node[i][j].X = Adjustment_Data[nb].Xwest[i][j];
	    Adjustment_Data[nb].Xwest[i][j] = Vector2D(MILLION,MILLION);
	  }
	  Mesh[nb].node_status[i][j] = NODE_STATUS_ALIGNED;
	}

      }

      // Reset intersect flags.
      Adjustment_Data[nb].intersect_north[i][j] = OFF;
      Adjustment_Data[nb].intersect_south[i][j] = OFF;
      Adjustment_Data[nb].intersect_east[i][j] = OFF;
      Adjustment_Data[nb].intersect_west[i][j] = OFF;

    }
  }

#ifdef _EB_PARALLEL_DEBUG_
//   dout << endl << "AFTER FIRST ADJUSTMENT::Block Number " << Local_Solution_Block_List->Block[nb].gblknum << endl;
//   for (int j = Local_SolnBlk[nb].Grid.JNu+Local_SolnBlk[nb].Nghost; j >= Local_SolnBlk[nb].Grid.JNl-Local_SolnBlk[nb].Nghost; j--) {
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       if (Mesh[nb].node_status[i][j] == NODE_STATUS_KNOWN) dout << "x";
//       else if (Mesh[nb].node_status[i][j] == NODE_STATUS_ALIGNED) dout << "o";
//       else dout << "?";
//     }
//     dout << endl;
//   }
//   dout << endl;
//   for (int j = Local_SolnBlk[nb].Grid.JNu+Local_SolnBlk[nb].Nghost; j >= Local_SolnBlk[nb].Grid.JNl-Local_SolnBlk[nb].Nghost; j--) {
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       if (Adjustment_Data[nb].intersect_north[i][j]) dout << endl << " North intersection:" << Local_SolnBlk[nb].Grid.Node[i][j].X << Adjustment_Data[nb].Xnorth[i][j];
//       if (Adjustment_Data[nb].intersect_south[i][j]) dout << endl << " South intersection:" << Local_SolnBlk[nb].Grid.Node[i][j].X << Adjustment_Data[nb].Xsouth[i][j];
//       if (Adjustment_Data[nb].intersect_east[i][j]) dout << endl << " East intersection:" << Local_SolnBlk[nb].Grid.Node[i][j].X << Adjustment_Data[nb].Xeast[i][j];
//       if (Adjustment_Data[nb].intersect_west[i][j]) dout << endl << " West intersection:" << Local_SolnBlk[nb].Grid.Node[i][j].X << Adjustment_Data[nb].Xwest[i][j];
//     }
//   }
//   dout << endl;
#endif

  // First mesh adjustment successful.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Mesh_Adjustment_Second --                    *
 *                                                                    *
 * This routine performs the second mesh adjustment on the specified  *
 * solution block.  This adjustment is exactly the same as the first  *
 * adjustment without some of the reductions are restrictions.  This  *
 * second adjustment pass rarely produces any changes, however, it    *
 * greatly improves the robustness of the adjustment in cases it is   *
 * required (especially in coarse meshes).                            *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Mesh_Adjustment_Second(const int &nb) {

  // Do not adjust if no interfaces have been defined.
  if (!Interface_Union_List.Ni) return 0;

  // Do not perform adjustment algorithm on the current solution block 
  // if there are no inactive or unknowns cells.
  if (!Adjustment_Data[nb].Number_of_Inactive_Cells) return 0;

  // Declare local variables.
  int intersect_flag;
  Vector2D Xm1, Xm2, Xs1, Xs2, Xp;

  // Second adjustment: Determine intersection flags and intersection points.
  for (int j = Local_SolnBlk[nb].JCl-Local_SolnBlk[nb].Nghost; j <= Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost; j++) {
    for (int i = Local_SolnBlk[nb].ICl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].ICu+Local_SolnBlk[nb].Nghost; i++) {

      // Continue with the current cell if one of the associated cells
      // is tagged as unknown.
      if (Cell_Status_Test(nb,i,j)) {

	// For each interface in the union list.
	for (int ni = 1; ni <= Interface_Union_List.Ni; ni++) {

	  // Determine if the interface is present in the current block.
	  if (Adjustment_Data[nb].Interface_Present[ni]) {

	    // For each (linear) edge of the interface.
	    for (int np = 0; np < Interface_Union_List[ni].Spline.np-1; np++) {

	      // Store the (linear) edge of the interface.
	      Xs1 = Interface_Union_List[ni].Spline.Xp[np  ];
	      Xs2 = Interface_Union_List[ni].Spline.Xp[np+1];

	      // Define the intersection flags and the intersection points
	      // for the SOUH-WEST node of the cell if not already aligned.
	      if (Mesh[nb].node_status[i][j] != NODE_STATUS_ALIGNED) {
		Xm1 = Local_SolnBlk[nb].Grid.Node[i][j].X;
		// NORTH edge.
		if (!((j == Local_SolnBlk[nb].Grid.JNu && Local_SolnBlk[nb].Grid.BCtypeN[i] != BC_NONE) ||
		      (j == Local_SolnBlk[nb].Grid.JNl && Local_SolnBlk[nb].Grid.BCtypeS[i] != BC_NONE)) &&
		    Mesh[nb].node_status[i][j+1] == NODE_STATUS_ALIGNED) {
		  Xm2 = Local_SolnBlk[nb].Grid.Node[i][j+1].X;
		  intersect_flag = Line_Intersection(Xs1,Xs2,Xm1,Xm2,Xp);
		  if (intersect_flag && abs(Xp-Xm1) < abs(Adjustment_Data[nb].Xnorth[i][j]-Xm1) && abs(Xm2-Xp) > TOLER) {
		    Adjustment_Data[nb].intersect_north[i][j] = intersect_flag;
		    Adjustment_Data[nb].Xnorth[i][j] = Xp;
		  }
		}
  	      	// SOUTH edge.
  		if (j > Local_SolnBlk[nb].JCl-Local_SolnBlk[nb].Nghost) {
  		  if (!((j == Local_SolnBlk[nb].Grid.JNu && Local_SolnBlk[nb].Grid.BCtypeN[i] != BC_NONE) ||
			(j == Local_SolnBlk[nb].Grid.JNl && Local_SolnBlk[nb].Grid.BCtypeS[i] != BC_NONE)) &&
		      Mesh[nb].node_status[i][j-1] == NODE_STATUS_ALIGNED) {
		    Xm2 = Local_SolnBlk[nb].Grid.Node[i][j-1].X;
		    intersect_flag = Line_Intersection(Xs1,Xs2,Xm1,Xm2,Xp);
		    if (intersect_flag && abs(Xp-Xm1) < abs(Adjustment_Data[nb].Xsouth[i][j]-Xm1) && abs(Xm2-Xp) > TOLER) {
		      Adjustment_Data[nb].intersect_south[i][j] = intersect_flag;
		      Adjustment_Data[nb].Xsouth[i][j] = Xp;
		    }
		  }
 		}
 		// EAST edge.
		if (!((i == Local_SolnBlk[nb].Grid.INu && Local_SolnBlk[nb].Grid.BCtypeE[j] != BC_NONE) ||
		      (i == Local_SolnBlk[nb].Grid.INl && Local_SolnBlk[nb].Grid.BCtypeW[j] != BC_NONE)) &&
		    Mesh[nb].node_status[i+1][j] == NODE_STATUS_ALIGNED) {
		  Xm2 = Local_SolnBlk[nb].Grid.Node[i+1][j].X;
		  intersect_flag = Line_Intersection(Xs1,Xs2,Xm1,Xm2,Xp);
		  if (intersect_flag && abs(Xp-Xm1) < abs(Adjustment_Data[nb].Xeast[i][j]-Xm1) && abs(Xm2-Xp) > TOLER) {
		    Adjustment_Data[nb].intersect_east[i][j] = intersect_flag;
		    Adjustment_Data[nb].Xeast[i][j] = Xp;
		  }
		}
 		// WEST edge.
 		if (i > Local_SolnBlk[nb].ICl-Local_SolnBlk[nb].Nghost) {
 		  if (!((i == Local_SolnBlk[nb].Grid.INu && Local_SolnBlk[nb].Grid.BCtypeE[j] != BC_NONE) ||
			(i == Local_SolnBlk[nb].Grid.INl && Local_SolnBlk[nb].Grid.BCtypeW[j] != BC_NONE)) &&
		      Mesh[nb].node_status[i-1][j] == NODE_STATUS_ALIGNED) {
		    Xm2 = Local_SolnBlk[nb].Grid.Node[i-1][j].X;
 		    intersect_flag = Line_Intersection(Xs1,Xs2,Xm1,Xm2,Xp);
 		    if (intersect_flag && abs(Xp-Xm1) < abs(Adjustment_Data[nb].Xwest[i][j]-Xm1) && abs(Xm2-Xp) > TOLER) {
 		      Adjustment_Data[nb].intersect_west[i][j] = intersect_flag;
 		      Adjustment_Data[nb].Xwest[i][j] = Xp;
 		    }
 		  }
 		}
	      }

	      // Define the intersection flags and the intersection points 
	      // for the SOUTH-EAST node if necessary.
	      if (i == Local_SolnBlk[nb].ICu+Local_SolnBlk[nb].Nghost &&
		  Mesh[nb].node_status[i+1][j] != NODE_STATUS_ALIGNED) {
		Xm1 = Local_SolnBlk[nb].Grid.Node[i+1][j].X;
		// NORTH edge.
		if (!((j == Local_SolnBlk[nb].Grid.JNu && Local_SolnBlk[nb].Grid.BCtypeN[i] != BC_NONE) ||
		      (j == Local_SolnBlk[nb].Grid.JNl && Local_SolnBlk[nb].Grid.BCtypeS[i] != BC_NONE)) &&
		    Mesh[nb].node_status[i+1][j+1] == NODE_STATUS_ALIGNED) {
		  Xm2 = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
		  intersect_flag = Line_Intersection(Xs1,Xs2,Xm1,Xm2,Xp);
		  if (intersect_flag && abs(Xp-Xm1) < abs(Adjustment_Data[nb].Xnorth[i+1][j]-Xm1) && abs(Xm2-Xp) > TOLER) {
		    Adjustment_Data[nb].intersect_north[i+1][j] = intersect_flag;
		    Adjustment_Data[nb].Xnorth[i+1][j] = Xp;
		  }
		}
		// SOUTH edge.
		if (j > Local_SolnBlk[nb].JCl-Local_SolnBlk[nb].Nghost) {
 		  if (!((j == Local_SolnBlk[nb].Grid.JNu && Local_SolnBlk[nb].Grid.BCtypeN[i] != BC_NONE) ||
			(j == Local_SolnBlk[nb].Grid.JNl && Local_SolnBlk[nb].Grid.BCtypeS[i] != BC_NONE)) &&
		      Mesh[nb].node_status[i+1][j-1] == NODE_STATUS_ALIGNED) {
		    Xm2 = Local_SolnBlk[nb].Grid.Node[i+1][j-1].X;
		    intersect_flag = Line_Intersection(Xs1,Xs2,Xm1,Xm2,Xp);
		    if (intersect_flag && abs(Xp-Xm1) < abs(Adjustment_Data[nb].Xsouth[i+1][j]-Xm1) && abs(Xm2-Xp) > TOLER) {
		      Adjustment_Data[nb].intersect_south[i+1][j] = intersect_flag;
		      Adjustment_Data[nb].Xsouth[i+1][j] = Xp;
		    }
		  }
		}
		// WEST edge.
		if (!((i+1 == Local_SolnBlk[nb].Grid.INu && Local_SolnBlk[nb].Grid.BCtypeE[j] != BC_NONE) ||
		      (i+1 == Local_SolnBlk[nb].Grid.INl && Local_SolnBlk[nb].Grid.BCtypeW[j] != BC_NONE)) &&
		    Mesh[nb].node_status[i][j] == NODE_STATUS_ALIGNED) {
		  Xm2 = Local_SolnBlk[nb].Grid.Node[i][j].X;
		  intersect_flag = Line_Intersection(Xs1,Xs2,Xm1,Xm2,Xp);
		  if (intersect_flag && abs(Xp-Xm1) < abs(Adjustment_Data[nb].Xwest[i+1][j]-Xm1) && abs(Xm2-Xp) > TOLER) {
		    Adjustment_Data[nb].intersect_west[i+1][j] = intersect_flag;
		    Adjustment_Data[nb].Xwest[i+1][j] = Xp;
		  }
		}
	      }

	      // Define the intersection flags and the intersection points
	      // for the NORTH-WEST node if necessary.
	      if (j == Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost &&
		  Mesh[nb].node_status[i][j+1] != NODE_STATUS_ALIGNED) {
		Xm1 = Local_SolnBlk[nb].Grid.Node[i][j+1].X;
		// SOUTH edge.
		if (!((j+1 == Local_SolnBlk[nb].Grid.JNu && Local_SolnBlk[nb].Grid.BCtypeN[i] != BC_NONE) ||
		      (j+1 == Local_SolnBlk[nb].Grid.JNl && Local_SolnBlk[nb].Grid.BCtypeS[i] != BC_NONE)) &&
		    Mesh[nb].node_status[i][j] == NODE_STATUS_ALIGNED) {
		  Xm2 = Local_SolnBlk[nb].Grid.Node[i][j].X;
		  intersect_flag = Line_Intersection(Xs1,Xs2,Xm1,Xm2,Xp);
		  if (intersect_flag && abs(Xp-Xm1) < abs(Adjustment_Data[nb].Xsouth[i][j+1]-Xm1) && abs(Xm2-Xp) > TOLER) {
		    Adjustment_Data[nb].intersect_south[i][j+1] = intersect_flag;
		    Adjustment_Data[nb].Xsouth[i][j+1] = Xp;
		  }
		}
		// EAST edge.
		if (!((i == Local_SolnBlk[nb].Grid.INu && Local_SolnBlk[nb].Grid.BCtypeE[j] != BC_NONE) ||
		      (i == Local_SolnBlk[nb].Grid.INl && Local_SolnBlk[nb].Grid.BCtypeW[j] != BC_NONE)) &&
		    Mesh[nb].node_status[i+1][j+1] == NODE_STATUS_ALIGNED) {
		  Xm2 = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
		  intersect_flag = Line_Intersection(Xs1,Xs2,Xm1,Xm2,Xp);
		  if (intersect_flag && abs(Xp-Xm1) < abs(Adjustment_Data[nb].Xeast[i][j+1]-Xm1) && abs(Xm2-Xp) > TOLER) {
		    Adjustment_Data[nb].intersect_east[i][j+1] = intersect_flag;
		    Adjustment_Data[nb].Xeast[i][j+1] = Xp;
		  }
		}
		// WEST edge.
		if (i > Local_SolnBlk[nb].ICl-Local_SolnBlk[nb].Nghost) {
		  if (!((i == Local_SolnBlk[nb].Grid.INu && Local_SolnBlk[nb].Grid.BCtypeE[j] != BC_NONE) ||
			(i == Local_SolnBlk[nb].Grid.INl && Local_SolnBlk[nb].Grid.BCtypeW[j] != BC_NONE)) &&
		      Mesh[nb].node_status[i-1][j+1] == NODE_STATUS_ALIGNED) {
		    Xm2 = Local_SolnBlk[nb].Grid.Node[i-1][j+1].X;
		    intersect_flag = Line_Intersection(Xs1,Xs2,Xm1,Xm2,Xp);
		    if (intersect_flag && abs(Xp-Xm1) < abs(Adjustment_Data[nb].Xwest[i][j+1]-Xm1) && abs(Xm2-Xp) > TOLER) {
		      Adjustment_Data[nb].intersect_west[i][j+1] = intersect_flag;
		      Adjustment_Data[nb].Xwest[i][j+1] = Xp;
		    }
		  }
		}
	      }

	      // Define the intersection flags and the intersection points
	      // for the NORTH-EAST node if necessary.
	      if (i == Local_SolnBlk[nb].ICu+Local_SolnBlk[nb].Nghost &&
		  j == Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost &&
		  Mesh[nb].node_status[i+1][j+1] != NODE_STATUS_ALIGNED) {
		Xm1 = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
		// SOUTH edge.
		if (!((j == Local_SolnBlk[nb].Grid.JNu && Local_SolnBlk[nb].Grid.BCtypeN[i] != BC_NONE) ||
		      (j == Local_SolnBlk[nb].Grid.JNl && Local_SolnBlk[nb].Grid.BCtypeS[i] != BC_NONE)) &&
		    Mesh[nb].node_status[i+1][j] == NODE_STATUS_ALIGNED) {
		  Xm2 = Local_SolnBlk[nb].Grid.Node[i+1][j].X;
		  intersect_flag = Line_Intersection(Xs1,Xs2,Xm1,Xm2,Xp);
		  if (intersect_flag && abs(Xp-Xm1) < abs(Adjustment_Data[nb].Xsouth[i+1][j+1]-Xm1) && abs(Xm2-Xp) > TOLER) {
		    Adjustment_Data[nb].intersect_south[i+1][j+1] = intersect_flag;
		    Adjustment_Data[nb].Xsouth[i+1][j+1] = Xp;
		  }
		}
		// WEST edge.
		if (!((i == Local_SolnBlk[nb].Grid.INu && Local_SolnBlk[nb].Grid.BCtypeE[j] != BC_NONE) ||
		      (i == Local_SolnBlk[nb].Grid.INl && Local_SolnBlk[nb].Grid.BCtypeW[j] != BC_NONE)) &&
		    Mesh[nb].node_status[i][j+1] == NODE_STATUS_ALIGNED) {
		  Xm2 = Local_SolnBlk[nb].Grid.Node[i][j+1].X;
		  intersect_flag = Line_Intersection(Xs1,Xs2,Xm1,Xm2,Xp);
		  if (intersect_flag && abs(Xp-Xm1) < abs(Adjustment_Data[nb].Xwest[i][j+1]-Xm1)) {
		    Adjustment_Data[nb].intersect_west[i+1][j+1] = intersect_flag;
		    Adjustment_Data[nb].Xwest[i+1][j+1] = Xp;
		  }
		}
	      }

	    }
	  }
	}
      }
    }
  }

#ifdef _EB_PARALLEL_DEBUG_
//   dout << endl << "SECOND ADJUSTMENT FLAGS::Block Number " << Local_Solution_Block_List->Block[nb].gblknum << endl;
//   for (int j = Local_SolnBlk[nb].Grid.JNu+Local_SolnBlk[nb].Nghost; j >= Local_SolnBlk[nb].Grid.JNl-Local_SolnBlk[nb].Nghost; j--) {
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       dout << Mesh[nb].node_status[i][j]+1;
//     }
//     dout << " ";
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       dout << Adjustment_Data[nb].intersect_north[i][j];
//     }
//     dout << " ";
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       dout << Adjustment_Data[nb].intersect_south[i][j];
//     }
//     dout << " ";
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       dout << Adjustment_Data[nb].intersect_east[i][j];
//     }
//     dout << " ";
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       dout << Adjustment_Data[nb].intersect_west[i][j];
//     }
//     dout << endl;
//   }
//   dout << endl;
//   for (int j = Local_SolnBlk[nb].Grid.JNu+Local_SolnBlk[nb].Nghost; j >= Local_SolnBlk[nb].Grid.JNl-Local_SolnBlk[nb].Nghost; j--) {
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       if (Adjustment_Data[nb].intersect_north[i][j]) dout << endl << " North intersection:" << Local_SolnBlk[nb].Grid.Node[i][j].X << Adjustment_Data[nb].Xnorth[i][j];
//       if (Adjustment_Data[nb].intersect_south[i][j]) dout << endl << " South intersection:" << Local_SolnBlk[nb].Grid.Node[i][j].X << Adjustment_Data[nb].Xsouth[i][j];
//       if (Adjustment_Data[nb].intersect_east[i][j]) dout << endl << " East intersection:" << Local_SolnBlk[nb].Grid.Node[i][j].X << Adjustment_Data[nb].Xeast[i][j];
//       if (Adjustment_Data[nb].intersect_west[i][j]) dout << endl << " West intersection:" << Local_SolnBlk[nb].Grid.Node[i][j].X << Adjustment_Data[nb].Xwest[i][j];
//     }
//   }
//   dout << endl;
#endif

  // Second Adjustment: Adjust cells according to the intersection flags.
  for (int j = Local_SolnBlk[nb].Grid.JNl-Local_SolnBlk[nb].Nghost; j <= Local_SolnBlk[nb].Grid.JNu+Local_SolnBlk[nb].Nghost; j++) {
    for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {

      if (Mesh[nb].node_status[i][j] != NODE_STATUS_ALIGNED) {
 	// Determine the closest interesection point and adjust node accordingly.
 	if (Adjustment_Data[nb].intersect_north[i][j] || 
 	    Adjustment_Data[nb].intersect_south[i][j] ||
 	    Adjustment_Data[nb].intersect_east[i][j] ||
 	    Adjustment_Data[nb].intersect_west[i][j]) {
 	  if (abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xnorth[i][j]) <=
 	      abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xsouth[i][j]) &&
 	      abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xnorth[i][j]) <=
	      abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xeast[i][j]) &&
	      abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xnorth[i][j]) <=
 	      abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xwest[i][j])) {
 	    Local_SolnBlk[nb].Grid.Node[i][j].X = Adjustment_Data[nb].Xnorth[i][j];
 	  } else if (abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xsouth[i][j]) <=
 		     abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xnorth[i][j]) &&
 		     abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xsouth[i][j]) <=
 		     abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xeast[i][j]) &&
		     abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xsouth[i][j]) <=
 		     abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xwest[i][j])) {
 	    Local_SolnBlk[nb].Grid.Node[i][j].X = Adjustment_Data[nb].Xsouth[i][j];
 	  } else if (abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xeast[i][j]) <=
 		     abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xnorth[i][j]) &&
 		     abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xeast[i][j]) <=
 		     abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xsouth[i][j]) &&
 		     abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xeast[i][j]) <=
 		     abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xwest[i][j])) {
 	    Local_SolnBlk[nb].Grid.Node[i][j].X = Adjustment_Data[nb].Xeast[i][j];
 	  } else if (abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xwest[i][j]) <=
 		     abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xnorth[i][j]) &&
 		     abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xwest[i][j]) <=
 		     abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xsouth[i][j]) &&
 		     abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xwest[i][j]) <=
 		     abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Adjustment_Data[nb].Xeast[i][j])) {
  	    Local_SolnBlk[nb].Grid.Node[i][j].X = Adjustment_Data[nb].Xwest[i][j];
 	  }
 	  Mesh[nb].node_status[i][j] = NODE_STATUS_ALIGNED;
 	}
      }

      // Reset intersect flags.
      Adjustment_Data[nb].intersect_north[i][j] = OFF;
      Adjustment_Data[nb].intersect_south[i][j] = OFF;
      Adjustment_Data[nb].intersect_east[i][j] = OFF;
      Adjustment_Data[nb].intersect_west[i][j] = OFF;

    }
  }

#ifdef _EB_PARALLEL_DEBUG_
//   dout << endl << "AFTER SECOND ADJUSTMENT::Block Number " << Local_Solution_Block_List->Block[nb].gblknum << endl;
//   for (int j = Local_SolnBlk[nb].Grid.JNu+Local_SolnBlk[nb].Nghost; j >= Local_SolnBlk[nb].Grid.JNl-Local_SolnBlk[nb].Nghost; j--) {
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       if (Mesh[nb].node_status[i][j] == NODE_STATUS_KNOWN) dout << "x";
//       else if (Mesh[nb].node_status[i][j] == NODE_STATUS_ALIGNED) dout << "o";
//       else dout << "?";
//     }
//     dout << endl;
//   }
//   dout << endl;
//   dout << endl;
//   for (int j = Local_SolnBlk[nb].Grid.JNu+Local_SolnBlk[nb].Nghost; j >= Local_SolnBlk[nb].Grid.JNl-Local_SolnBlk[nb].Nghost; j--) {
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       dout << Mesh[nb].node_status[i][j]+1;
//     }
//     dout << " ";
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       dout << Adjustment_Data[nb].intersect_north[i][j];
//     }
//     dout << " ";
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       dout << Adjustment_Data[nb].intersect_south[i][j];
//     }
//     dout << " ";
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       dout << Adjustment_Data[nb].intersect_east[i][j];
//     }
//     dout << " ";
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       dout << Adjustment_Data[nb].intersect_west[i][j];
//     }
//     dout << endl;
//   }
//   dout << endl;
#endif

  // Second mesh adjustment successful.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Mesh_Adjustment_Third --                     *
 *                                                                    *
 * This routine performs the third mesh adjustment which accounts for *
 * cells that are bisected by the embedded boundary of the first and  *
 * second adjustments.  The node nearest to the embedded boundary is  *
 * moved onto the embedded boundary and is re-tagged as 'aligned.'    *
 * Note that this step of the mesh adjustment can produce degenerate  *
 * quadrilateral cells (also known as triangles).  The degenerate     *
 * edge (zero length) is always located on the embedded boundary.     *
 *                                                                    *
 * NW/SE cell  o-----?  NE/SW cell  ?-----o  Move one of the          *
 *             |\    |              |    /|  unadjusted nodes such    *
 *             | \   |              |   / |  that the cell forms a    *
 *             |  \  |              |  /  |  triangle.  Note that the *
 *             |   \ |              | /   |  degenerate edge is       *
 *             |    \|              |/    |  always located on the    *
 *             ?-----o              o-----?  embedded boundary.       *
 *                                                                    *
 * NE cell  o-----?  NW cell  ?-----o  Move the unadjusted node such  *
 *          |\    |           |    /|  that the cell forms a          *
 *          | \   |           |   / |  triangle.  Note that the       *
 *          |  \  |           |  /  |  degenerate edge is always      *
 *          |   \ |           | /   |  located on the embedded        *
 *          |    \|           |/    |  boundary.                      *
 *          o-----o           o-----o                                 *
 *                                                                    *
 * SE cell  o-----o  SW cell  o-----o                                 *
 *          |    /|           |\    |                                 *
 *          |   / |           | \   |                                 *
 *          |  /  |           |  \  |                                 *
 *          | /   |           |   \ |                                 *
 *          |/    |           |    \|                                 *
 *          o-----?           ?-----o                                 *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Mesh_Adjustment_Third(const int &nb) {

  // Do not adjust if no interfaces have been defined.
  if (!Interface_Union_List.Ni) return 0;

  // Do not perform adjustment algorithm on the current solution block 
  // if there are no inactive or unknowns cells.
  if (!Adjustment_Data[nb].Number_of_Inactive_Cells) return 0;

  // Declare local variables.
  int ii, jj;
  Vector2D Xm1, Xm2, Xp;

  // Reevaluate intersection flags to filter out cells split by the
  // embedded boundary/interface:
  for (int j = Local_SolnBlk[nb].JCl-Local_SolnBlk[nb].Nghost; j <= Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost; j++) {
    for (int i = Local_SolnBlk[nb].ICl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].ICu+Local_SolnBlk[nb].Nghost; i++) {

      // Continue only if the current or neighbour cell status is
      // unknown since there is no need to adjust otherwise.
      if (Cell_Status_Test(nb,i,j)) {

	// NW/SE split cell.
	if (Mesh[nb].node_status[i  ][j+1] == NODE_STATUS_ALIGNED &&
	    Mesh[nb].node_status[i+1][j+1] != NODE_STATUS_ALIGNED &&
	    Mesh[nb].node_status[i  ][j  ] != NODE_STATUS_ALIGNED &&
	    Mesh[nb].node_status[i+1][j  ] == NODE_STATUS_ALIGNED) {

 	  // Set SW node intersection flags:
 	  // Store current node in a temporary variable.
 	  Xp = Local_SolnBlk[nb].Grid.Node[i][j].X;
 	  // Attempt to adjust north.
	  if (!(j == Local_SolnBlk[nb].Grid.JNl && Local_SolnBlk[nb].Grid.BCtypeS[i] != BC_NONE) &&
	      !(j == Local_SolnBlk[nb].Grid.JNu && Local_SolnBlk[nb].Grid.BCtypeN[i] != BC_NONE)) {
	    Local_SolnBlk[nb].Grid.Node[i][j].X = Local_SolnBlk[nb].Grid.Node[i][j+1].X;
	    if (!Check_Cell_Areas(nb,i,j)) Adjustment_Data[nb].intersect_north[i][j] = ON;
	  }
	  // Attempt to adjust east.
	  if (!(i == Local_SolnBlk[nb].Grid.INl && Local_SolnBlk[nb].Grid.BCtypeW[j] != BC_NONE) &&
	      !(i == Local_SolnBlk[nb].Grid.INu && Local_SolnBlk[nb].Grid.BCtypeE[j] != BC_NONE)) {
	    Local_SolnBlk[nb].Grid.Node[i][j].X = Local_SolnBlk[nb].Grid.Node[i+1][j].X;
	    if (!Check_Cell_Areas(nb,i,j)) Adjustment_Data[nb].intersect_east[i][j] = ON;
	  }
	  // Restore node.
	  Local_SolnBlk[nb].Grid.Node[i][j].X = Xp;
	  // Reduce adjustment choice.
	  if (Adjustment_Data[nb].intersect_north[i][j] && Adjustment_Data[nb].intersect_east[i][j]) {
	    if (abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Local_SolnBlk[nb].Grid.Node[i][j+1].X) <
		abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Local_SolnBlk[nb].Grid.Node[i+1][j].X)) {
	      Adjustment_Data[nb].intersect_east[i][j] = OFF;
	    } else {
	      Adjustment_Data[nb].intersect_north[i][j] = OFF;
	    }
	  }

	  // Set NE node intersection flags:
	  // Store current node in a temporary variable.
	  Xp = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
	  // Attempt to adjust south.
	  if (!(j+1 == Local_SolnBlk[nb].Grid.JNl && Local_SolnBlk[nb].Grid.BCtypeS[i+1] != BC_NONE) &&
	      !(j+1 == Local_SolnBlk[nb].Grid.JNu && Local_SolnBlk[nb].Grid.BCtypeN[i+1] != BC_NONE)) {
	    Local_SolnBlk[nb].Grid.Node[i+1][j+1].X = Local_SolnBlk[nb].Grid.Node[i+1][j].X;
	    if (!Check_Cell_Areas(nb,i,j)) Adjustment_Data[nb].intersect_south[i+1][j+1] = ON;
	  }
	  // Attempt to adjust west.
	  if (!(i+1 == Local_SolnBlk[nb].Grid.INl && Local_SolnBlk[nb].Grid.BCtypeW[j+1] != BC_NONE) &&
	      !(i+1 == Local_SolnBlk[nb].Grid.INu && Local_SolnBlk[nb].Grid.BCtypeE[j+1] != BC_NONE)) {
	    Local_SolnBlk[nb].Grid.Node[i+1][j+1].X = Local_SolnBlk[nb].Grid.Node[i][j+1].X;
	    if (!Check_Cell_Areas(nb,i,j)) Adjustment_Data[nb].intersect_west[i+1][j+1] = ON;
	  }
	  // Restore node.
	  Local_SolnBlk[nb].Grid.Node[i+1][j+1].X = Xp;
	  // Reduce adjustment choice.
	  if (Adjustment_Data[nb].intersect_south[i+1][j+1] && Adjustment_Data[nb].intersect_west[i+1][j+1]) {
	    if (abs(Local_SolnBlk[nb].Grid.Node[i+1][j+1].X - Local_SolnBlk[nb].Grid.Node[i+1][j].X) <
		abs(Local_SolnBlk[nb].Grid.Node[i+1][j+1].X - Local_SolnBlk[nb].Grid.Node[i][j+1].X)) {
	      Adjustment_Data[nb].intersect_west[i+1][j+1] = OFF;
	    } else {
	      Adjustment_Data[nb].intersect_south[i+1][j+1] = OFF;
	    }
	  }

	}

	// NE/SW split cell.
	if (Mesh[nb].node_status[i  ][j+1] != NODE_STATUS_ALIGNED &&
	    Mesh[nb].node_status[i+1][j+1] == NODE_STATUS_ALIGNED &&
	    Mesh[nb].node_status[i  ][j  ] == NODE_STATUS_ALIGNED &&
	    Mesh[nb].node_status[i+1][j  ] != NODE_STATUS_ALIGNED) {

	  // Set SE node intersection flags:
	  // Store current node in a temporary variable.
	  Xp = Local_SolnBlk[nb].Grid.Node[i+1][j].X;
	  // Attempt to adjust north.
	  if (!(j == Local_SolnBlk[nb].Grid.JNl && Local_SolnBlk[nb].Grid.BCtypeS[i+1] != BC_NONE) &&
	      !(j == Local_SolnBlk[nb].Grid.JNu && Local_SolnBlk[nb].Grid.BCtypeN[i+1] != BC_NONE)) {
	    Local_SolnBlk[nb].Grid.Node[i+1][j].X = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
	    if (!Check_Cell_Areas(nb,i,j)) Adjustment_Data[nb].intersect_north[i+1][j] = ON;
	  }
	  // Attempt to adjust west.
	  if (!(i+1 == Local_SolnBlk[nb].Grid.INl && Local_SolnBlk[nb].Grid.BCtypeW[j] != BC_NONE) &&
	      !(i+1 == Local_SolnBlk[nb].Grid.INu && Local_SolnBlk[nb].Grid.BCtypeE[j] != BC_NONE)) {
	    Local_SolnBlk[nb].Grid.Node[i+1][j].X = Local_SolnBlk[nb].Grid.Node[i][j].X;
	    if (!Check_Cell_Areas(nb,i,j)) Adjustment_Data[nb].intersect_west[i+1][j] = ON;
	  }
	  // Restore node.
	  Local_SolnBlk[nb].Grid.Node[i+1][j].X = Xp;
	  // Reduce adjustment choice.
	  if (Adjustment_Data[nb].intersect_north[i+1][j] && Adjustment_Data[nb].intersect_west[i+1][j]) {
	    if (abs(Local_SolnBlk[nb].Grid.Node[i+1][j].X - Local_SolnBlk[nb].Grid.Node[i+1][j+1].X) <
		abs(Local_SolnBlk[nb].Grid.Node[i+1][j].X - Local_SolnBlk[nb].Grid.Node[i][j].X)) {
	      Adjustment_Data[nb].intersect_west[i+1][j] = OFF;
	    } else {
	      Adjustment_Data[nb].intersect_north[i+1][j] = OFF;
	    }
	  }

	  // Set NW node intersection flags:
	  // Store current node in a temporary variable.
	  Xp = Local_SolnBlk[nb].Grid.Node[i][j+1].X;
	  // Attempt to adjust south.
	  if (!(j+1 == Local_SolnBlk[nb].Grid.JNl && Local_SolnBlk[nb].Grid.BCtypeS[i] != BC_NONE) &&
	      !(j+1 == Local_SolnBlk[nb].Grid.JNu && Local_SolnBlk[nb].Grid.BCtypeN[i] != BC_NONE)) {
	    Local_SolnBlk[nb].Grid.Node[i][j+1].X = Local_SolnBlk[nb].Grid.Node[i][j].X;
	    if (!Check_Cell_Areas(nb,i,j)) Adjustment_Data[nb].intersect_south[i][j+1] = ON;
	  }
	  // Attempt to adjust east.
	  if (!(i == Local_SolnBlk[nb].Grid.INl && Local_SolnBlk[nb].Grid.BCtypeW[j+1] != BC_NONE) &&
	      !(i == Local_SolnBlk[nb].Grid.INu && Local_SolnBlk[nb].Grid.BCtypeE[j+1] != BC_NONE)) {
	    Local_SolnBlk[nb].Grid.Node[i][j+1].X = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
	    if (!Check_Cell_Areas(nb,i,j)) Adjustment_Data[nb].intersect_east[i][j+1] = ON;
	  }
	  // Restore node.
	  Local_SolnBlk[nb].Grid.Node[i][j+1].X = Xp;
	  // Reduce adjustment choice.
	  if (Adjustment_Data[nb].intersect_south[i][j+1] && Adjustment_Data[nb].intersect_east[i][j+1]) {
	    if (abs(Local_SolnBlk[nb].Grid.Node[i][j+1].X - Local_SolnBlk[nb].Grid.Node[i][j].X) <
		abs(Local_SolnBlk[nb].Grid.Node[i][j+1].X - Local_SolnBlk[nb].Grid.Node[i+1][j+1].X)) {
	      Adjustment_Data[nb].intersect_east[i][j+1] = OFF;
	    } else {
	      Adjustment_Data[nb].intersect_south[i][j+1] = OFF;
	    }
	  }

	}

	// NE split cell.
	if (Mesh[nb].node_status[i  ][j+1] == NODE_STATUS_ALIGNED &&
	    Mesh[nb].node_status[i+1][j+1] != NODE_STATUS_ALIGNED &&
	    Mesh[nb].node_status[i  ][j  ] == NODE_STATUS_ALIGNED &&
	    Mesh[nb].node_status[i+1][j  ] == NODE_STATUS_ALIGNED) {
	  // Determine the centroids of the sub-triangles.
	  Xm1 = (Local_SolnBlk[nb].Grid.Node[i+1][j+1].X + Local_SolnBlk[nb].Grid.Node[i+1][j].X + Local_SolnBlk[nb].Grid.Node[i][j+1].X)/THREE;
	  Xm2 = (Local_SolnBlk[nb].Grid.Node[i  ][j  ].X + Local_SolnBlk[nb].Grid.Node[i+1][j].X + Local_SolnBlk[nb].Grid.Node[i][j+1].X)/THREE;
	  // Determine if the centroids are internal to any of the
	  // interfaces.
	  ii = 0; jj = 0;
	  for (int ni = 1; ni <= Interface_Union_List.Ni; ni++) {
	    if (Interface_Union_List[ni].Point_In_Interface(Xm1)) ii = 1;
	    if (Interface_Union_List[ni].Point_In_Interface(Xm2)) jj = 1;
	  }
	  // If one of the centroids is internal and the other is
	  // external then determine the intersection flag.
	  if ((ii && !jj) || (!ii && jj)) {
	    // Store current node in a temporary variable.
	    Xp = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
	    // Attempt to adjust south.
	    if (!(j+1 == Local_SolnBlk[nb].Grid.JNl && Local_SolnBlk[nb].Grid.BCtypeS[i+1] != BC_NONE) &&
		!(j+1 == Local_SolnBlk[nb].Grid.JNu && Local_SolnBlk[nb].Grid.BCtypeN[i+1] != BC_NONE)) {
	      Local_SolnBlk[nb].Grid.Node[i+1][j+1].X = Local_SolnBlk[nb].Grid.Node[i+1][j].X;
	      if (!Check_Cell_Areas(nb,i,j)) Adjustment_Data[nb].intersect_south[i+1][j+1] = ON;
	    }
	    // Attempt to adjust west.
	    if (!(i+1 == Local_SolnBlk[nb].Grid.INl && Local_SolnBlk[nb].Grid.BCtypeW[j+1] != BC_NONE) &&
		!(i+1 == Local_SolnBlk[nb].Grid.INu && Local_SolnBlk[nb].Grid.BCtypeE[j+1] != BC_NONE)) {
	      Local_SolnBlk[nb].Grid.Node[i+1][j+1].X = Local_SolnBlk[nb].Grid.Node[i][j+1].X;
	      if (!Check_Cell_Areas(nb,i,j)) Adjustment_Data[nb].intersect_west[i+1][j+1] = ON;
	    }
	    // Restore node.
	    Local_SolnBlk[nb].Grid.Node[i+1][j+1].X = Xp;
	    // Reduce adjustment choice.
	    if (Adjustment_Data[nb].intersect_south[i+1][j+1] &&
		Adjustment_Data[nb].intersect_west[i+1][j+1]) {
	      if (abs(Local_SolnBlk[nb].Grid.Node[i+1][j+1].X - Local_SolnBlk[nb].Grid.Node[i+1][j].X) <
		  abs(Local_SolnBlk[nb].Grid.Node[i+1][j+1].X - Local_SolnBlk[nb].Grid.Node[i][j+1].X)) {
		Adjustment_Data[nb].intersect_west[i+1][j+1] = OFF;
	      } else {
		Adjustment_Data[nb].intersect_south[i+1][j+1] = OFF;
	      }
	    }
	  }
	}

	// NW split cell.
	if (Mesh[nb].node_status[i  ][j+1] != NODE_STATUS_ALIGNED &&
	    Mesh[nb].node_status[i+1][j+1] == NODE_STATUS_ALIGNED &&
	    Mesh[nb].node_status[i  ][j  ] == NODE_STATUS_ALIGNED &&
	    Mesh[nb].node_status[i+1][j  ] == NODE_STATUS_ALIGNED) {
	  // Determine the centroids of the sub-triangles.
	  Xm1 = (Local_SolnBlk[nb].Grid.Node[i][j+1].X + Local_SolnBlk[nb].Grid.Node[i][j].X + Local_SolnBlk[nb].Grid.Node[i+1][j+1].X)/THREE;
	  Xm2 = (Local_SolnBlk[nb].Grid.Node[i+1][j].X + Local_SolnBlk[nb].Grid.Node[i][j].X + Local_SolnBlk[nb].Grid.Node[i+1][j+1].X)/THREE;
	  // Determine if the centroids are internal to any of the
	  // interfaces.
	  ii = 0; jj = 0;
	  for (int ni = 1; ni <= Interface_Union_List.Ni; ni++) {
	    if (Interface_Union_List[ni].Point_In_Interface(Xm1)) ii = 1;
	    if (Interface_Union_List[ni].Point_In_Interface(Xm2)) jj = 1;
	  }
	  // If one of the centroids is internal and the other is
	  // external then determine the intersection flag.
	  if ((ii && !jj) || (!ii && jj)) {
	    // Store current node in a temporary variable.
	    Xp = Local_SolnBlk[nb].Grid.Node[i][j+1].X;
	    // Attempt to adjust south.
	    if (!(j+1 == Local_SolnBlk[nb].Grid.JNl && Local_SolnBlk[nb].Grid.BCtypeS[i] != BC_NONE) &&
		!(j+1 == Local_SolnBlk[nb].Grid.JNu && Local_SolnBlk[nb].Grid.BCtypeN[i] != BC_NONE)) {
	      Local_SolnBlk[nb].Grid.Node[i][j+1].X = Local_SolnBlk[nb].Grid.Node[i][j].X;
	      if (!Check_Cell_Areas(nb,i,j)) Adjustment_Data[nb].intersect_south[i][j+1] = ON;
	    }
	    // Attempt to adjust east.
	    if (!(i == Local_SolnBlk[nb].Grid.INl && Local_SolnBlk[nb].Grid.BCtypeW[j+1] != BC_NONE) &&
		!(i == Local_SolnBlk[nb].Grid.INu && Local_SolnBlk[nb].Grid.BCtypeE[j+1] != BC_NONE)) {
	      Local_SolnBlk[nb].Grid.Node[i][j+1].X = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
	      if (!Check_Cell_Areas(nb,i,j)) Adjustment_Data[nb].intersect_east[i][j+1] = ON;
	    }
	    // Restore node.
	    Local_SolnBlk[nb].Grid.Node[i][j+1].X = Xp;
	    // Reduce adjustment choice.
	    if (Adjustment_Data[nb].intersect_south[i][j+1] &&
		Adjustment_Data[nb].intersect_east[i][j+1]) {
	      if (abs(Local_SolnBlk[nb].Grid.Node[i][j+1].X - Local_SolnBlk[nb].Grid.Node[i][j].X) <
		  abs(Local_SolnBlk[nb].Grid.Node[i][j+1].X - Local_SolnBlk[nb].Grid.Node[i+1][j+1].X)) {
		Adjustment_Data[nb].intersect_east[i][j+1] = OFF;
	      } else {
		Adjustment_Data[nb].intersect_south[i][j+1] = OFF;
	      }
	    }
	  }
	}

	// SE split cell.
	if (Mesh[nb].node_status[i  ][j+1] == NODE_STATUS_ALIGNED &&
	    Mesh[nb].node_status[i+1][j+1] == NODE_STATUS_ALIGNED &&
	    Mesh[nb].node_status[i  ][j  ] == NODE_STATUS_ALIGNED &&
	    Mesh[nb].node_status[i+1][j  ] != NODE_STATUS_ALIGNED) {
	  // Determine the centroids of the sub-triangles.
	  Xm1 = (Local_SolnBlk[nb].Grid.Node[i+1][j].X + Local_SolnBlk[nb].Grid.Node[i][j].X + Local_SolnBlk[nb].Grid.Node[i+1][j+1].X)/THREE;
	  Xm2 = (Local_SolnBlk[nb].Grid.Node[i][j+1].X + Local_SolnBlk[nb].Grid.Node[i][j].X + Local_SolnBlk[nb].Grid.Node[i+1][j+1].X)/THREE;
	  // Determine if the centroids are internal to any of the
	  // interfaces.
	  ii = 0; jj = 0;
	  for (int ni = 1; ni <= Interface_Union_List.Ni; ni++) {
	    if (Interface_Union_List[ni].Point_In_Interface(Xm1)) ii = 1;
	    if (Interface_Union_List[ni].Point_In_Interface(Xm2)) jj = 1;
	  }
	  // If one of the centroids is internal and the other is
	  // external then determine the intersection flag.
	  if ((ii && !jj) || (!ii && jj)) {
	    // Store current node in a temporary variable.
	    Xp = Local_SolnBlk[nb].Grid.Node[i+1][j].X;
	    // Attempt to adjust north.
	    if (!(j == Local_SolnBlk[nb].Grid.JNl && Local_SolnBlk[nb].Grid.BCtypeS[i+1] != BC_NONE) ||
		!(j == Local_SolnBlk[nb].Grid.JNu && Local_SolnBlk[nb].Grid.BCtypeN[i+1] != BC_NONE)) {
	      Local_SolnBlk[nb].Grid.Node[i+1][j].X = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
	      if (!Check_Cell_Areas(nb,i,j)) Adjustment_Data[nb].intersect_north[i+1][j] = ON;
	    }
	    // Attempt to adjust west.
	    if (!(i+1 == Local_SolnBlk[nb].Grid.INl && Local_SolnBlk[nb].Grid.BCtypeW[j] != BC_NONE) &&
		!(i+1 == Local_SolnBlk[nb].Grid.INu && Local_SolnBlk[nb].Grid.BCtypeE[j] != BC_NONE)) {
	      Local_SolnBlk[nb].Grid.Node[i+1][j].X = Local_SolnBlk[nb].Grid.Node[i][j].X;
	      if (!Check_Cell_Areas(nb,i,j)) Adjustment_Data[nb].intersect_west[i+1][j] = ON;
	    }
	    // Restore node.
	    Local_SolnBlk[nb].Grid.Node[i+1][j].X = Xp;
	    // Reduce adjustment choice.
	    if (Adjustment_Data[nb].intersect_north[i+1][j] &&
		Adjustment_Data[nb].intersect_west[i+1][j]) {
	      if (abs(Local_SolnBlk[nb].Grid.Node[i+1][j].X - Local_SolnBlk[nb].Grid.Node[i+1][j+1].X) <
		  abs(Local_SolnBlk[nb].Grid.Node[i+1][j].X - Local_SolnBlk[nb].Grid.Node[i][j].X)) {
		Adjustment_Data[nb].intersect_west[i+1][j] = OFF;
	      } else {
		Adjustment_Data[nb].intersect_north[i+1][j] = OFF;
	      }
	    }
	  }
	}

	// SW split cell.
	if (Mesh[nb].node_status[i  ][j+1] == NODE_STATUS_ALIGNED &&
	    Mesh[nb].node_status[i+1][j+1] == NODE_STATUS_ALIGNED &&
	    Mesh[nb].node_status[i  ][j  ] != NODE_STATUS_ALIGNED &&
	    Mesh[nb].node_status[i+1][j  ] == NODE_STATUS_ALIGNED) {
	  // Determine the centroids of the sub-triangles.
	  Xm1 = (Local_SolnBlk[nb].Grid.Node[i  ][j  ].X + Local_SolnBlk[nb].Grid.Node[i+1][j].X + Local_SolnBlk[nb].Grid.Node[i][j+1].X)/THREE;
	  Xm2 = (Local_SolnBlk[nb].Grid.Node[i+1][j+1].X + Local_SolnBlk[nb].Grid.Node[i+1][j].X + Local_SolnBlk[nb].Grid.Node[i][j+1].X)/THREE;
	  // Determine if the centroids are internal to any of the
	  // interfaces.
	  ii = 0; jj = 0;
	  for (int ni = 1; ni <= Interface_Union_List.Ni; ni++) {
	    if (Interface_Union_List[ni].Point_In_Interface(Xm1)) ii = 1;
	    if (Interface_Union_List[ni].Point_In_Interface(Xm2)) jj = 1;
	  }
	  // If one of the centroids is internal and the other is
	  // external then determine the intersection flag.
	  if ((ii && !jj) || (!ii && jj)) {
	    // Store current node in a temporary variable.
	    Xp = Local_SolnBlk[nb].Grid.Node[i][j].X;
	    // Attempt to adjust north.
	    if (!(j == Local_SolnBlk[nb].Grid.JNl && Local_SolnBlk[nb].Grid.BCtypeS[i] != BC_NONE) &&
		!(j == Local_SolnBlk[nb].Grid.JNu && Local_SolnBlk[nb].Grid.BCtypeN[i] != BC_NONE)) {
	      Local_SolnBlk[nb].Grid.Node[i][j].X = Local_SolnBlk[nb].Grid.Node[i][j+1].X;
	      if (!Check_Cell_Areas(nb,i,j)) Adjustment_Data[nb].intersect_north[i][j] = ON;
	    }
	    // Attempt to adjust east.
	    if (!(i == Local_SolnBlk[nb].Grid.INl && Local_SolnBlk[nb].Grid.BCtypeW[j] != BC_NONE) &&
		!(i == Local_SolnBlk[nb].Grid.INu && Local_SolnBlk[nb].Grid.BCtypeE[j] != BC_NONE)) {
	      Local_SolnBlk[nb].Grid.Node[i][j].X = Local_SolnBlk[nb].Grid.Node[i+1][j].X;
	      if (!Check_Cell_Areas(nb,i,j)) Adjustment_Data[nb].intersect_east[i][j] = ON;
	    }
	    // Restore node.
	    Local_SolnBlk[nb].Grid.Node[i][j].X = Xp;
	    // Reduce adjustment choice.
	    if (Adjustment_Data[nb].intersect_north[i][j] &&
		Adjustment_Data[nb].intersect_east[i][j]) {
	      if (abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Local_SolnBlk[nb].Grid.Node[i][j+1].X) <
		  abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Local_SolnBlk[nb].Grid.Node[i+1][j].X)) {
		Adjustment_Data[nb].intersect_east[i][j] = OFF;
	      } else {
		Adjustment_Data[nb].intersect_north[i][j] = OFF;
	      }
	    }
	  }
	}

      }
    }
  }

#ifdef _EB_PARALLEL_DEBUG_
//   dout << endl << "THIRD ADJUSTMENT FLAGS (before reduction)::Block Number " << Local_Solution_Block_List->Block[nb].gblknum << endl;
//   for (int j = Local_SolnBlk[nb].Grid.JNu+Local_SolnBlk[nb].Nghost; j >= Local_SolnBlk[nb].Grid.JNl-Local_SolnBlk[nb].Nghost; j--) {
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       dout << Mesh[nb].node_status[i][j]+1;
//     }
//     dout << " ";
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       dout << Adjustment_Data[nb].intersect_north[i][j];
//     }
//     dout << " ";
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       dout << Adjustment_Data[nb].intersect_south[i][j];
//     }
//     dout << " ";
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       dout << Adjustment_Data[nb].intersect_east[i][j];
//     }
//     dout << " ";
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       dout << Adjustment_Data[nb].intersect_west[i][j];
//     }
//     dout << endl;
//   }
//   dout << endl;
//   dout << endl << " NORTH";
//   for (int j = Local_SolnBlk[nb].Grid.JNu+Local_SolnBlk[nb].Nghost; j >= Local_SolnBlk[nb].Grid.JNl-Local_SolnBlk[nb].Nghost; j--) {
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       if (Adjustment_Data[nb].intersect_north[i][j]) dout << endl << Local_SolnBlk[nb].Grid.Node[i][j].X;
//     }
//   }
//   dout << endl << " SOUTH";
//   for (int j = Local_SolnBlk[nb].Grid.JNu+Local_SolnBlk[nb].Nghost; j >= Local_SolnBlk[nb].Grid.JNl-Local_SolnBlk[nb].Nghost; j--) {
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       if (Adjustment_Data[nb].intersect_south[i][j]) dout << endl << Local_SolnBlk[nb].Grid.Node[i][j].X;
//     }
//   }
//   dout << endl << " EAST";
//   for (int j = Local_SolnBlk[nb].Grid.JNu+Local_SolnBlk[nb].Nghost; j >= Local_SolnBlk[nb].Grid.JNl-Local_SolnBlk[nb].Nghost; j--) {
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       if (Adjustment_Data[nb].intersect_east[i][j]) dout << endl << Local_SolnBlk[nb].Grid.Node[i][j].X;
//     }
//   }
//   dout << endl << " WEST";
//   for (int j = Local_SolnBlk[nb].Grid.JNu+Local_SolnBlk[nb].Nghost; j >= Local_SolnBlk[nb].Grid.JNl-Local_SolnBlk[nb].Nghost; j--) {
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       if (Adjustment_Data[nb].intersect_west[i][j]) dout << endl << Local_SolnBlk[nb].Grid.Node[i][j].X;
//     }
//   }
//   dout << endl;
#endif

  // Adjustment choice reduction:
  // -> For nodes that intersect to the NORTH and SOUTH or to the EAST 
  //    and WEST, turn of intersect flags for both competing choices 
  //    and adjust the neighbours instead (important for aerofoils).
  for (int j = Local_SolnBlk[nb].Grid.JNl-Local_SolnBlk[nb].Nghost; j <= Local_SolnBlk[nb].Grid.JNu+Local_SolnBlk[nb].Nghost; j++) {
    for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
      if (Adjustment_Data[nb].intersect_north[i][j] &&
	  Adjustment_Data[nb].intersect_south[i][j]) {
	if (j > Local_SolnBlk[nb].Grid.JNl-Local_SolnBlk[nb].Nghost &&
	    j < Local_SolnBlk[nb].Grid.JNu+Local_SolnBlk[nb].Nghost) {
	  if (fabs(abs(Local_SolnBlk[nb].Grid.Node[i][j+1].X-Local_SolnBlk[nb].Grid.Node[i][j].X) -
		   abs(Local_SolnBlk[nb].Grid.Node[i][j-1].X-Local_SolnBlk[nb].Grid.Node[i][j].X)) <
	      0.25*abs(Local_SolnBlk[nb].Grid.Node[i][j+1].X-Local_SolnBlk[nb].Grid.Node[i][j-1].X)) {
	    Adjustment_Data[nb].intersect_north[i][j] = OFF;
	    Adjustment_Data[nb].intersect_south[i][j] = OFF;
	  }
	}
      }
      if (Adjustment_Data[nb].intersect_east[i][j] &&
	  Adjustment_Data[nb].intersect_west[i][j]) {
	if (i > Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost &&
	    i < Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost) {
	  if (fabs(abs(Local_SolnBlk[nb].Grid.Node[i+1][j].X-Local_SolnBlk[nb].Grid.Node[i][j].X) -
		   abs(Local_SolnBlk[nb].Grid.Node[i-1][j].X-Local_SolnBlk[nb].Grid.Node[i][j].X)) <
	      0.25*abs(Local_SolnBlk[nb].Grid.Node[i+1][j].X-Local_SolnBlk[nb].Grid.Node[i-1][j].X)) {
	    Adjustment_Data[nb].intersect_east[i][j] = OFF;
	    Adjustment_Data[nb].intersect_west[i][j] = OFF;
	  }
	}
      }
      if (j < Local_SolnBlk[nb].Grid.JNu+Local_SolnBlk[nb].Nghost) {
	if (Adjustment_Data[nb].intersect_west[i][j] &&
	    Adjustment_Data[nb].intersect_west[i][j+1]) {
 	  if (abs(Local_SolnBlk[nb].Grid.Cell[i][j].Xc - Adjustment_Data[nb].Xwest[i][j]) <
 	      abs(Local_SolnBlk[nb].Grid.Cell[i][j+1].Xc - Adjustment_Data[nb].Xwest[i][j+1])) {
 	    Adjustment_Data[nb].intersect_west[i][j] = OFF;
 	  } else {
 	    Adjustment_Data[nb].intersect_west[i][j+1] = OFF;
 	  }
	}
      }
      if (j > Local_SolnBlk[nb].Grid.JNl-Local_SolnBlk[nb].Nghost &&
	  j < Local_SolnBlk[nb].Grid.JNu+Local_SolnBlk[nb].Nghost &&
	  i < Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost) {
	if (Adjustment_Data[nb].intersect_north[i][j-1] &&
	    Adjustment_Data[nb].intersect_south[i][j+1] &&
	    Adjustment_Data[nb].intersect_west[i+1][j]) {
	  Adjustment_Data[nb].intersect_north[i][j-1] = OFF;
	  Adjustment_Data[nb].intersect_south[i][j+1] = OFF;
	}
      }
      if (j > Local_SolnBlk[nb].Grid.JNl-Local_SolnBlk[nb].Nghost &&
	  j < Local_SolnBlk[nb].Grid.JNu+Local_SolnBlk[nb].Nghost &&
	  i > Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost) {
	if (Adjustment_Data[nb].intersect_north[i][j-1] &&
	    Adjustment_Data[nb].intersect_south[i][j+1] &&
	    Adjustment_Data[nb].intersect_east[i-1][j]) {
	  Adjustment_Data[nb].intersect_north[i][j-1] = OFF;
	  Adjustment_Data[nb].intersect_south[i][j+1] = OFF;
	}
      }
    }
  }

#ifdef _EB_PARALLEL_DEBUG_
//   dout << endl << "THIRD ADJUSTMENT FLAGS::Block Number " << Local_Solution_Block_List->Block[nb].gblknum << endl;
//   for (int j = Local_SolnBlk[nb].Grid.JNu+Local_SolnBlk[nb].Nghost; j >= Local_SolnBlk[nb].Grid.JNl-Local_SolnBlk[nb].Nghost; j--) {
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       dout << Mesh[nb].node_status[i][j]+1;
//     }
//     dout << " ";
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       dout << Adjustment_Data[nb].intersect_north[i][j];
//     }
//     dout << " ";
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       dout << Adjustment_Data[nb].intersect_south[i][j];
//     }
//     dout << " ";
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       dout << Adjustment_Data[nb].intersect_east[i][j];
//     }
//     dout << " ";
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       dout << Adjustment_Data[nb].intersect_west[i][j];
//     }
//     dout << endl;
//   }
//   dout << endl;
//   dout << endl;
//   for (int j = Local_SolnBlk[nb].Grid.JNu+Local_SolnBlk[nb].Nghost; j >= Local_SolnBlk[nb].Grid.JNl-Local_SolnBlk[nb].Nghost; j--) {
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       dout << Mesh[nb].node_status[i][j]+1;
//     }
//     dout << " ";
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       if (Adjustment_Data[nb].intersect_north[i][j] && Adjustment_Data[nb].intersect_south[i][j] &&
// 	  Adjustment_Data[nb].intersect_east[i][j] && Adjustment_Data[nb].intersect_west[i][j]) {
// 	dout << 1;
//       } else if (!Adjustment_Data[nb].intersect_north[i][j] && Adjustment_Data[nb].intersect_south[i][j] &&
// 		 Adjustment_Data[nb].intersect_east[i][j] && Adjustment_Data[nb].intersect_west[i][j]) {
// 	dout << 2;
//       } else if (Adjustment_Data[nb].intersect_north[i][j] && !Adjustment_Data[nb].intersect_south[i][j] &&
// 		 Adjustment_Data[nb].intersect_east[i][j] && Adjustment_Data[nb].intersect_west[i][j]) {
// 	dout << 3;
//       } else if (Adjustment_Data[nb].intersect_north[i][j] && Adjustment_Data[nb].intersect_south[i][j] &&
// 		 !Adjustment_Data[nb].intersect_east[i][j] && Adjustment_Data[nb].intersect_west[i][j]) {
// 	dout << 4;
//       } else if (Adjustment_Data[nb].intersect_north[i][j] && Adjustment_Data[nb].intersect_south[i][j] &&
// 		 Adjustment_Data[nb].intersect_east[i][j] && !Adjustment_Data[nb].intersect_west[i][j]) {
// 	dout << 5;
//       } else if (Adjustment_Data[nb].intersect_north[i][j] && Adjustment_Data[nb].intersect_south[i][j] &&
// 		 !Adjustment_Data[nb].intersect_east[i][j] && !Adjustment_Data[nb].intersect_west[i][j]) {
// 	dout << 6;
//       } else if (Adjustment_Data[nb].intersect_north[i][j] && !Adjustment_Data[nb].intersect_south[i][j] &&
// 		 Adjustment_Data[nb].intersect_east[i][j] && !Adjustment_Data[nb].intersect_west[i][j]) {
// 	dout << 7;
//       } else if (Adjustment_Data[nb].intersect_north[i][j] && !Adjustment_Data[nb].intersect_south[i][j] &&
// 		 !Adjustment_Data[nb].intersect_east[i][j] && Adjustment_Data[nb].intersect_west[i][j]) {
// 	dout << 8;
//       } else if (!Adjustment_Data[nb].intersect_north[i][j] && Adjustment_Data[nb].intersect_south[i][j] &&
// 		 Adjustment_Data[nb].intersect_east[i][j] && !Adjustment_Data[nb].intersect_west[i][j]) {
// 	dout << 9;
//       } else if (!Adjustment_Data[nb].intersect_north[i][j] && Adjustment_Data[nb].intersect_south[i][j] &&
// 		 !Adjustment_Data[nb].intersect_east[i][j] && Adjustment_Data[nb].intersect_west[i][j]) {
// 	dout << 0;
//       } else if (Adjustment_Data[nb].intersect_north[i][j] && Adjustment_Data[nb].intersect_south[i][j] &&
// 		 Adjustment_Data[nb].intersect_east[i][j] && Adjustment_Data[nb].intersect_west[i][j]) {
// 	dout << "a";
//       } else if (Adjustment_Data[nb].intersect_north[i][j] && !Adjustment_Data[nb].intersect_south[i][j] &&
// 		 !Adjustment_Data[nb].intersect_east[i][j] && !Adjustment_Data[nb].intersect_west[i][j]) {
// 	dout << "n";
//       } else if (!Adjustment_Data[nb].intersect_north[i][j] && Adjustment_Data[nb].intersect_south[i][j] &&
// 		 !Adjustment_Data[nb].intersect_east[i][j] && !Adjustment_Data[nb].intersect_west[i][j]) {
// 	dout << "s";
//       } else if (!Adjustment_Data[nb].intersect_north[i][j] && !Adjustment_Data[nb].intersect_south[i][j] &&
// 		 Adjustment_Data[nb].intersect_east[i][j] && !Adjustment_Data[nb].intersect_west[i][j]) {
// 	dout << "e";
//       } else if (!Adjustment_Data[nb].intersect_north[i][j] && !Adjustment_Data[nb].intersect_south[i][j] &&
// 		 !Adjustment_Data[nb].intersect_east[i][j] && Adjustment_Data[nb].intersect_west[i][j]) {
// 	dout << "w";
//       } else {
// 	dout << "x";
//       }
//     }
//     dout << endl;
//   }
//   dout << endl << " NORTH";
//   for (int j = Local_SolnBlk[nb].Grid.JNu+Local_SolnBlk[nb].Nghost; j >= Local_SolnBlk[nb].Grid.JNl-Local_SolnBlk[nb].Nghost; j--) {
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       if (Adjustment_Data[nb].intersect_north[i][j]) {
// 	dout << endl << Local_SolnBlk[nb].Grid.Node[i][j].X << Local_SolnBlk[nb].Grid.Node[i][j+1].X << " " << abs(Local_SolnBlk[nb].Grid.Node[i][j].X-Local_SolnBlk[nb].Grid.Node[i][j+1].X);
// 	if (j > Local_SolnBlk[nb].Grid.JNl-Local_SolnBlk[nb].Nghost) dout << " " << abs(Local_SolnBlk[nb].Grid.Node[i][j+1].X-Local_SolnBlk[nb].Grid.Node[i][j-1].X);
//       }
//     }
//   }
//   dout << endl << " SOUTH";
//   for (int j = Local_SolnBlk[nb].Grid.JNu+Local_SolnBlk[nb].Nghost; j >= Local_SolnBlk[nb].Grid.JNl-Local_SolnBlk[nb].Nghost; j--) {
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       if (Adjustment_Data[nb].intersect_south[i][j]) {
// 	dout << endl << Local_SolnBlk[nb].Grid.Node[i][j].X << Local_SolnBlk[nb].Grid.Node[i][j-1].X << " " << abs(Local_SolnBlk[nb].Grid.Node[i][j].X-Local_SolnBlk[nb].Grid.Node[i][j-1].X);
// 	if (j < Local_SolnBlk[nb].Grid.JNu+Local_SolnBlk[nb].Nghost) dout << " " << abs(Local_SolnBlk[nb].Grid.Node[i][j+1].X-Local_SolnBlk[nb].Grid.Node[i][j-1].X);
//       }
//     }
//   }
//   dout << endl << " EAST";
//   for (int j = Local_SolnBlk[nb].Grid.JNu+Local_SolnBlk[nb].Nghost; j >= Local_SolnBlk[nb].Grid.JNl-Local_SolnBlk[nb].Nghost; j--) {
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       if (Adjustment_Data[nb].intersect_east[i][j]) {
// 	dout << endl << Local_SolnBlk[nb].Grid.Node[i][j].X << Local_SolnBlk[nb].Grid.Node[i+1][j].X << " " << abs(Local_SolnBlk[nb].Grid.Node[i][j].X-Local_SolnBlk[nb].Grid.Node[i+1][j].X);
// 	if (i > Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost) dout << " " << abs(Local_SolnBlk[nb].Grid.Node[i+1][j].X-Local_SolnBlk[nb].Grid.Node[i-1][j].X);
//       }
//     }
//   }
//   dout << endl << " WEST";
//   for (int j = Local_SolnBlk[nb].Grid.JNu+Local_SolnBlk[nb].Nghost; j >= Local_SolnBlk[nb].Grid.JNl-Local_SolnBlk[nb].Nghost; j--) {
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       if (Adjustment_Data[nb].intersect_west[i][j]) {
// 	dout << endl << Local_SolnBlk[nb].Grid.Node[i][j].X << Local_SolnBlk[nb].Grid.Node[i-1][j].X << " " << abs(Local_SolnBlk[nb].Grid.Node[i][j].X-Local_SolnBlk[nb].Grid.Node[i-1][j].X);
// 	if (i < Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost) dout << " " << abs(Local_SolnBlk[nb].Grid.Node[i+1][j].X-Local_SolnBlk[nb].Grid.Node[i-1][j].X);
//       }
//     }
//   }
//   dout << endl;
#endif

  // Third Adjustment: Adjust cells according to the intersection flags.
  for (int j = Local_SolnBlk[nb].JCl-Local_SolnBlk[nb].Nghost; j <= Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost; j++) {
    for (int i = Local_SolnBlk[nb].ICl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].ICu+Local_SolnBlk[nb].Nghost; i++) {

      // Continue only if the current or neighbour cell status is
      // unknown since there is no need to adjust otherwise.
      if (Cell_Status_Test(nb,i,j)) {

	if (Mesh[nb].node_status[i  ][j+1] == NODE_STATUS_ALIGNED &&
	    Mesh[nb].node_status[i+1][j+1] != NODE_STATUS_ALIGNED &&
	    Mesh[nb].node_status[i  ][j  ] != NODE_STATUS_ALIGNED &&
	    Mesh[nb].node_status[i+1][j  ] == NODE_STATUS_ALIGNED) {

	  // NW/SE split cell.
	  if ((Adjustment_Data[nb].intersect_north[i][j] ||
	       Adjustment_Data[nb].intersect_east[i][j]) &&
	      (Adjustment_Data[nb].intersect_south[i+1][j+1] ||
	       Adjustment_Data[nb].intersect_west[i+1][j+1])) {

	    if (Adjustment_Data[nb].intersect_north[i][j] &&
		Adjustment_Data[nb].intersect_south[i+1][j+1]) {
	      if (abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Local_SolnBlk[nb].Grid.Node[i][j+1].X) <
		  abs(Local_SolnBlk[nb].Grid.Node[i+1][j+1].X - Local_SolnBlk[nb].Grid.Node[i+1][j].X)) {
		// Adjust the SW node north.
		Xp = Local_SolnBlk[nb].Grid.Node[i][j].X;
		Local_SolnBlk[nb].Grid.Node[i][j].X = Local_SolnBlk[nb].Grid.Node[i][j+1].X;
 		if (Check_Cell_Areas(nb,i,j)) {
		  // An invalid cell has been created. Adjust the NE node south.
		  Local_SolnBlk[nb].Grid.Node[i][j].X = Xp;
 		  Local_SolnBlk[nb].Grid.Node[i+1][j+1].X = Local_SolnBlk[nb].Grid.Node[i+1][j].X;
		  if (Check_Cell_Areas(nb,i,j)) return 1101;
		  Mesh[nb].node_status[i+1][j+1] = NODE_STATUS_ALIGNED;
		} else {
		  Mesh[nb].node_status[i][j] = NODE_STATUS_ALIGNED;
		}
	      } else if (abs(Local_SolnBlk[nb].Grid.Node[i+1][j+1].X - Local_SolnBlk[nb].Grid.Node[i+1][j].X) <
			 abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Local_SolnBlk[nb].Grid.Node[i][j+1].X)) {
		// Adjust the NE node south.
		Xp = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
		Local_SolnBlk[nb].Grid.Node[i+1][j+1].X = Local_SolnBlk[nb].Grid.Node[i+1][j].X;
 		if (Check_Cell_Areas(nb,i,j)) {
		  // An invalid cell has been created. Adjust the SW node north.
		  Local_SolnBlk[nb].Grid.Node[i+1][j+1].X = Xp;
		  Local_SolnBlk[nb].Grid.Node[i][j].X = Local_SolnBlk[nb].Grid.Node[i][j+1].X;
		  if (Check_Cell_Areas(nb,i,j)) return 1102;
		  Mesh[nb].node_status[i][j] = NODE_STATUS_ALIGNED;
		} else {
		  Mesh[nb].node_status[i+1][j+1] = NODE_STATUS_ALIGNED;
		}
	      } else {
		// Equal distance... adjust inactive node.
		return 1103;
	      }
	      Adjustment_Data[nb].Reset_Intersect_Flags(i,j);
	      Adjustment_Data[nb].Reset_Intersect_Flags(i+1,j+1);

	    } else if (Adjustment_Data[nb].intersect_north[i][j] &&
		       Adjustment_Data[nb].intersect_west[i+1][j+1]) {
	      if (abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Local_SolnBlk[nb].Grid.Node[i][j+1].X) <
		  abs(Local_SolnBlk[nb].Grid.Node[i+1][j+1].X - Local_SolnBlk[nb].Grid.Node[i][j+1].X)) {
		// Adjust the SW node north.
		Xp = Local_SolnBlk[nb].Grid.Node[i][j].X;
		Local_SolnBlk[nb].Grid.Node[i][j].X = Local_SolnBlk[nb].Grid.Node[i][j+1].X;
 		if (Check_Cell_Areas(nb,i,j)) {
		  // An invalid cell has been created. Adjust the NE node west.
		  Local_SolnBlk[nb].Grid.Node[i][j].X = Xp;
 		  Local_SolnBlk[nb].Grid.Node[i+1][j+1].X = Local_SolnBlk[nb].Grid.Node[i][j+1].X;
		  if (Check_Cell_Areas(nb,i,j)) return 1104;
		  Mesh[nb].node_status[i+1][j+1] = NODE_STATUS_ALIGNED;
		} else {
		  Mesh[nb].node_status[i][j] = NODE_STATUS_ALIGNED;
		}
	      } else if (abs(Local_SolnBlk[nb].Grid.Node[i+1][j+1].X - Local_SolnBlk[nb].Grid.Node[i][j+1].X) <
			 abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Local_SolnBlk[nb].Grid.Node[i][j+1].X)) {
		// Adjust the NE node west.
		Xp = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
		Local_SolnBlk[nb].Grid.Node[i+1][j+1].X = Local_SolnBlk[nb].Grid.Node[i][j+1].X;
 		if (Check_Cell_Areas(nb,i,j)) {
		  // An invalid cell has been created. Adjust the SW node north.
		  Local_SolnBlk[nb].Grid.Node[i+1][j+1].X = Xp;
		  Local_SolnBlk[nb].Grid.Node[i][j].X = Local_SolnBlk[nb].Grid.Node[i][j+1].X;
		  if (Check_Cell_Areas(nb,i,j)) return 1105;
		  Mesh[nb].node_status[i][j] = NODE_STATUS_ALIGNED;
		} else {
		  Mesh[nb].node_status[i+1][j+1] = NODE_STATUS_ALIGNED;
		}
	      } else {
		// Equal distance... adjust inactive node.
		return 1106;
 	      }
	      Adjustment_Data[nb].Reset_Intersect_Flags(i,j);
	      Adjustment_Data[nb].Reset_Intersect_Flags(i+1,j+1);

	    } else if (Adjustment_Data[nb].intersect_east[i][j] &&
		       Adjustment_Data[nb].intersect_south[i+1][j+1]) {
	      if (abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Local_SolnBlk[nb].Grid.Node[i+1][j].X) <
		  abs(Local_SolnBlk[nb].Grid.Node[i+1][j+1].X - Local_SolnBlk[nb].Grid.Node[i+1][j].X)) {
		// Adjust the SW node east.
		Xp = Local_SolnBlk[nb].Grid.Node[i][j].X;
		Local_SolnBlk[nb].Grid.Node[i][j].X = Local_SolnBlk[nb].Grid.Node[i+1][j].X;
 		if (Check_Cell_Areas(nb,i,j)) {
		  // An invalid cell has been created. Adjust the NE node south.
		  Local_SolnBlk[nb].Grid.Node[i][j].X = Xp;
 		  Local_SolnBlk[nb].Grid.Node[i+1][j+1].X = Local_SolnBlk[nb].Grid.Node[i+1][j].X;
		  if (Check_Cell_Areas(nb,i,j)) return 1107;
		  Mesh[nb].node_status[i+1][j+1] = NODE_STATUS_ALIGNED;
		} else {
		  Mesh[nb].node_status[i][j] = NODE_STATUS_ALIGNED;
		}
	      } else if (abs(Local_SolnBlk[nb].Grid.Node[i+1][j+1].X - Local_SolnBlk[nb].Grid.Node[i+1][j].X) <
			 abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Local_SolnBlk[nb].Grid.Node[i+1][j].X)) {
		// Adjust the NE node south.
		Xp = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
		Local_SolnBlk[nb].Grid.Node[i+1][j+1].X = Local_SolnBlk[nb].Grid.Node[i+1][j].X;
 		if (Check_Cell_Areas(nb,i,j)) {
		  // An invalid cell has been created. Adjust the SW node east.
		  Local_SolnBlk[nb].Grid.Node[i+1][j+1].X = Xp;
		  Local_SolnBlk[nb].Grid.Node[i][j].X = Local_SolnBlk[nb].Grid.Node[i+1][j].X;
		  if (Check_Cell_Areas(nb,i,j)) return 1108;
		  Mesh[nb].node_status[i][j] = NODE_STATUS_ALIGNED;
		} else {
		  Mesh[nb].node_status[i+1][j+1] = NODE_STATUS_ALIGNED;
		}
	      } else {
		// Equal distance... adjust inactive node.
		return 1109;
	      }
	      Adjustment_Data[nb].Reset_Intersect_Flags(i,j);
	      Adjustment_Data[nb].Reset_Intersect_Flags(i+1,j+1);

	    } else if (Adjustment_Data[nb].intersect_east[i][j] &&
		       Adjustment_Data[nb].intersect_west[i+1][j+1]) {
	      if (abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Local_SolnBlk[nb].Grid.Node[i+1][j].X) <
		  abs(Local_SolnBlk[nb].Grid.Node[i+1][j+1].X - Local_SolnBlk[nb].Grid.Node[i][j+1].X)) {
		// Adjust the SW node east.
		Xp = Local_SolnBlk[nb].Grid.Node[i][j].X;
		Local_SolnBlk[nb].Grid.Node[i][j].X = Local_SolnBlk[nb].Grid.Node[i+1][j].X;
 		if (Check_Cell_Areas(nb,i,j)) {
		  // An invalid cell has been created. Adjust the NE node west.
		  Local_SolnBlk[nb].Grid.Node[i][j].X = Xp;
 		  Local_SolnBlk[nb].Grid.Node[i+1][j+1].X = Local_SolnBlk[nb].Grid.Node[i][j+1].X;
		  if (Check_Cell_Areas(nb,i,j)) return 1110;
		  Mesh[nb].node_status[i+1][j+1] = NODE_STATUS_ALIGNED;
		} else {
		  Mesh[nb].node_status[i][j] = NODE_STATUS_ALIGNED;
		}
	      } else if (abs(Local_SolnBlk[nb].Grid.Node[i+1][j+1].X - Local_SolnBlk[nb].Grid.Node[i][j+1].X) <
			 abs(Local_SolnBlk[nb].Grid.Node[i][j].X - Local_SolnBlk[nb].Grid.Node[i+1][j].X)) {
		// Adjust the NE node west.
		Xp = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
		Local_SolnBlk[nb].Grid.Node[i+1][j+1].X = Local_SolnBlk[nb].Grid.Node[i][j+1].X;
 		if (Check_Cell_Areas(nb,i,j)) {
		  // An invalid cell has been created. Adjust the SW node east.
		  Local_SolnBlk[nb].Grid.Node[i+1][j+1].X = Xp;
		  Local_SolnBlk[nb].Grid.Node[i][j].X = Local_SolnBlk[nb].Grid.Node[i+1][j].X;
		  if (Check_Cell_Areas(nb,i,j)) return 1111;
		  Mesh[nb].node_status[i][j] = NODE_STATUS_ALIGNED;
		} else {
		  Mesh[nb].node_status[i+1][j+1] = NODE_STATUS_ALIGNED;
		}
	      } else {
		// Equal distance... adjust inactive node.
		return 1112;
	      }
	      Adjustment_Data[nb].Reset_Intersect_Flags(i,j);
	      Adjustment_Data[nb].Reset_Intersect_Flags(i+1,j+1);

	    } else {
	      // Not possible.  Return an error flag.
	      return 1113;

 	    }

	  } else if (Adjustment_Data[nb].intersect_north[i][j] ||
		     Adjustment_Data[nb].intersect_east[i][j]) {
	    // The SW node is the only adjustment candidate.
	    if (Adjustment_Data[nb].intersect_north[i][j]) {
	      // Adjust the SW node north.
	      Local_SolnBlk[nb].Grid.Node[i][j].X = Local_SolnBlk[nb].Grid.Node[i][j+1].X;
 	      if (Check_Cell_Areas(nb,i,j)) {
		// An invalid cell has been created. Adjust the SW node east.
 		Local_SolnBlk[nb].Grid.Node[i][j].X = Local_SolnBlk[nb].Grid.Node[i+1][j].X;
		if (Check_Cell_Areas(nb,i,j)) return 1114;
	      }
 	    } else if (Adjustment_Data[nb].intersect_east[i][j]) {
	      // Adjust the SW node east.
	      Local_SolnBlk[nb].Grid.Node[i][j].X = Local_SolnBlk[nb].Grid.Node[i+1][j].X;
 	      if (Check_Cell_Areas(nb,i,j)) {
		// An invalid cell has been created.  Adjust the SW node north.
 		Local_SolnBlk[nb].Grid.Node[i][j].X = Local_SolnBlk[nb].Grid.Node[i][j].X;
		if (Check_Cell_Areas(nb,i,j)) return 1115;
	      }
	    } else {
	      // Not possible.  Return an error flag.
	      return 1116;
	    }
	    // Set node status flag.
  	    Mesh[nb].node_status[i][j] = NODE_STATUS_ALIGNED;
	    Adjustment_Data[nb].Reset_Intersect_Flags(i,j);

	  } else if (Adjustment_Data[nb].intersect_south[i+1][j+1] ||
		     Adjustment_Data[nb].intersect_west[i+1][j+1]) {
	    // The NE node is the only adjustment candidate.
	    if (Adjustment_Data[nb].intersect_south[i+1][j+1]) {
	      // Adjust the NE node south.
	      Local_SolnBlk[nb].Grid.Node[i+1][j+1].X = Local_SolnBlk[nb].Grid.Node[i+1][j].X;
 	      if (Check_Cell_Areas(nb,i,j)) {
		// An invalid cell has been created.  Adjust the NE node west.
 		Local_SolnBlk[nb].Grid.Node[i+1][j+1].X = Local_SolnBlk[nb].Grid.Node[i][j+1].X;
		if (Check_Cell_Areas(nb,i,j)) return 1117;
	      }
	    } else if (Adjustment_Data[nb].intersect_west[i+1][j+1]) {
	      // Adjust the NE node west.
	      Local_SolnBlk[nb].Grid.Node[i+1][j+1].X = Local_SolnBlk[nb].Grid.Node[i][j+1].X;
 	      if (Check_Cell_Areas(nb,i,j)) {
		// An invalid cell has been created.  Adjust the NE node south.
 		Local_SolnBlk[nb].Grid.Node[i+1][j+1].X = Local_SolnBlk[nb].Grid.Node[i+1][j].X;
		if (Check_Cell_Areas(nb,i,j)) return 1118;
 	      }
	    } else {
	      // Not possible.  Return an error flag.
	      return 1119;
	    }
	    // Set node status flag.
 	    Mesh[nb].node_status[i+1][j+1] = NODE_STATUS_ALIGNED;
	    Adjustment_Data[nb].Reset_Intersect_Flags(i+1,j+1);

	  } else {
	    // Not possible.  Return an error flag.
	    //cout << endl << i << " " << j;
	    //cout << endl << Adjustment_Data[nb].intersect_north[i][j];
	    //cout << endl << Adjustment_Data[nb].intersect_east[i][j];
	    //cout << endl << Adjustment_Data[nb].intersect_south[i+1][j+1];
	    //cout << endl << Adjustment_Data[nb].intersect_west[i+1][j+1];
	    //cout << endl << Local_SolnBlk[nb].Grid.Node[i][j].X;
	    //cout << endl << Local_SolnBlk[nb].Grid.Node[i+1][j].X;
	    //cout << endl << Local_SolnBlk[nb].Grid.Node[i][j+1].X;
	    //cout << endl << Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
	    //cout.flush();
	    //return 1120;

	  }

	} else if (Mesh[nb].node_status[i  ][j+1] != NODE_STATUS_ALIGNED &&
		   Mesh[nb].node_status[i+1][j+1] == NODE_STATUS_ALIGNED &&
		   Mesh[nb].node_status[i  ][j  ] == NODE_STATUS_ALIGNED &&
		   Mesh[nb].node_status[i+1][j  ] != NODE_STATUS_ALIGNED) {

 	  // NE/SW split cell.
 	  if ((Adjustment_Data[nb].intersect_north[i+1][j] ||
 	       Adjustment_Data[nb].intersect_west[i+1][j]) &&
 	      (Adjustment_Data[nb].intersect_south[i][j+1] ||
 	       Adjustment_Data[nb].intersect_east[i][j+1])) {

	    if (Adjustment_Data[nb].intersect_north[i+1][j] &&
		Adjustment_Data[nb].intersect_south[i][j+1]) {
	      if (abs(Local_SolnBlk[nb].Grid.Node[i+1][j].X - Local_SolnBlk[nb].Grid.Node[i+1][j+1].X) <
		  abs(Local_SolnBlk[nb].Grid.Node[i][j+1].X - Local_SolnBlk[nb].Grid.Node[i][j].X)) {
		// Adjust the SE node north.
		Xp = Local_SolnBlk[nb].Grid.Node[i+1][j].X;
		Local_SolnBlk[nb].Grid.Node[i+1][j].X = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
 		if (Check_Cell_Areas(nb,i,j)) {
		  // An invalid cell has been created. Adjust the NW node south.
		  Local_SolnBlk[nb].Grid.Node[i+1][j].X = Xp;
 		  Local_SolnBlk[nb].Grid.Node[i][j+1].X = Local_SolnBlk[nb].Grid.Node[i][j].X;
		  if (Check_Cell_Areas(nb,i,j)) return 2101;
		  Mesh[nb].node_status[i][j+1] = NODE_STATUS_ALIGNED;
		} else {
		  Mesh[nb].node_status[i+1][j] = NODE_STATUS_ALIGNED;
		}
	      } else if (abs(Local_SolnBlk[nb].Grid.Node[i][j+1].X - Local_SolnBlk[nb].Grid.Node[i][j].X) <
			 abs(Local_SolnBlk[nb].Grid.Node[i+1][j].X - Local_SolnBlk[nb].Grid.Node[i+1][j+1].X)) {
		// Adjust the NW node south.
		Xp = Local_SolnBlk[nb].Grid.Node[i][j+1].X;
		Local_SolnBlk[nb].Grid.Node[i][j+1].X = Local_SolnBlk[nb].Grid.Node[i][j].X;
 		if (Check_Cell_Areas(nb,i,j)) {
		  // An invalid cell has been created. Adjust the SE node north.
		  Local_SolnBlk[nb].Grid.Node[i][j+1].X = Xp;
		  Local_SolnBlk[nb].Grid.Node[i+1][j].X = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
		  if (Check_Cell_Areas(nb,i,j)) return 2102;
		  Mesh[nb].node_status[i+1][j] = NODE_STATUS_ALIGNED;
		} else {
		  Mesh[nb].node_status[i][j+1] = NODE_STATUS_ALIGNED;
		}
	      } else {
		// Equal distance... adjust inactive node.
		return 2103;
	      }
	      Adjustment_Data[nb].Reset_Intersect_Flags(i+1,j);
	      Adjustment_Data[nb].Reset_Intersect_Flags(i,j+1);

	    } else if (Adjustment_Data[nb].intersect_north[i+1][j] &&
		       Adjustment_Data[nb].intersect_east[i][j+1]) {
	      if (abs(Local_SolnBlk[nb].Grid.Node[i+1][j].X - Local_SolnBlk[nb].Grid.Node[i+1][j+1].X) <
		  abs(Local_SolnBlk[nb].Grid.Node[i][j+1].X - Local_SolnBlk[nb].Grid.Node[i+1][j+1].X)) {
		// Adjust the SE node north.
		Xp = Local_SolnBlk[nb].Grid.Node[i+1][j].X;
		Local_SolnBlk[nb].Grid.Node[i+1][j].X = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
 		if (Check_Cell_Areas(nb,i,j)) {
		  // An invalid cell has been created. Adjust the NW node east.
		  Local_SolnBlk[nb].Grid.Node[i+1][j].X = Xp;
 		  Local_SolnBlk[nb].Grid.Node[i][j+1].X = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
		  if (Check_Cell_Areas(nb,i,j)) return 2104;
		  Mesh[nb].node_status[i][j+1] = NODE_STATUS_ALIGNED;
		} else {
		  Mesh[nb].node_status[i+1][j] = NODE_STATUS_ALIGNED;
		}
	      } else if (abs(Local_SolnBlk[nb].Grid.Node[i][j+1].X - Local_SolnBlk[nb].Grid.Node[i+1][j+1].X) <
			 abs(Local_SolnBlk[nb].Grid.Node[i+1][j].X - Local_SolnBlk[nb].Grid.Node[i+1][j+1].X)) {
		// Adjust the NW node east.
		Xp = Local_SolnBlk[nb].Grid.Node[i][j+1].X;
		Local_SolnBlk[nb].Grid.Node[i][j+1].X = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
 		if (Check_Cell_Areas(nb,i,j)) {
		  // An invalid cell has been created. Adjust the SE node north.
		  Local_SolnBlk[nb].Grid.Node[i][j+1].X = Xp;
		  Local_SolnBlk[nb].Grid.Node[i+1][j].X = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
		  if (Check_Cell_Areas(nb,i,j)) return 2105;
		  Mesh[nb].node_status[i+1][j] = NODE_STATUS_ALIGNED;
		} else {
		  Mesh[nb].node_status[i][j+1] = NODE_STATUS_ALIGNED;
		}
	      } else {
		// Equal distance... adjust inactive node.
		return 2106;
	      }
	      Adjustment_Data[nb].Reset_Intersect_Flags(i+1,j);
	      Adjustment_Data[nb].Reset_Intersect_Flags(i,j+1);

	    } else if (Adjustment_Data[nb].intersect_west[i+1][j] &&
		       Adjustment_Data[nb].intersect_south[i][j+1]) {
	      if (abs(Local_SolnBlk[nb].Grid.Node[i+1][j].X - Local_SolnBlk[nb].Grid.Node[i][j].X) <
		  abs(Local_SolnBlk[nb].Grid.Node[i][j+1].X - Local_SolnBlk[nb].Grid.Node[i][j].X)) {
		// Adjust the SE node west.
		Xp = Local_SolnBlk[nb].Grid.Node[i+1][j].X;
		Local_SolnBlk[nb].Grid.Node[i+1][j].X = Local_SolnBlk[nb].Grid.Node[i][j].X;
 		if (Check_Cell_Areas(nb,i,j)) {
		  // An invalid cell has been created. Adjust the NW node south.
		  Local_SolnBlk[nb].Grid.Node[i+1][j].X = Xp;
 		  Local_SolnBlk[nb].Grid.Node[i][j+1].X = Local_SolnBlk[nb].Grid.Node[i][j].X;
		  if (Check_Cell_Areas(nb,i,j)) return 2107;
		  Mesh[nb].node_status[i][j+1] = NODE_STATUS_ALIGNED;
		} else {
		  Mesh[nb].node_status[i+1][j] = NODE_STATUS_ALIGNED;
		}
	      } else if (abs(Local_SolnBlk[nb].Grid.Node[i][j+1].X - Local_SolnBlk[nb].Grid.Node[i][j].X) <
			 abs(Local_SolnBlk[nb].Grid.Node[i+1][j].X - Local_SolnBlk[nb].Grid.Node[i][j].X)) {
		// Adjust the NW node south.
		Xp = Local_SolnBlk[nb].Grid.Node[i][j+1].X;
		Local_SolnBlk[nb].Grid.Node[i][j+1].X = Local_SolnBlk[nb].Grid.Node[i][j].X;
 		if (Check_Cell_Areas(nb,i,j)) {
		  // An invalid cell has been created. Adjust the SE node west.
		  Local_SolnBlk[nb].Grid.Node[i][j+1].X = Xp;
		  Local_SolnBlk[nb].Grid.Node[i+1][j].X = Local_SolnBlk[nb].Grid.Node[i][j].X;
		  if (Check_Cell_Areas(nb,i,j)) return 2108;
		  Mesh[nb].node_status[i+1][j] = NODE_STATUS_ALIGNED;
		} else {
		  Mesh[nb].node_status[i][j+1] = NODE_STATUS_ALIGNED;
		}
	      } else {
		// Equal distance... adjust inactive node.
		return 2109;
	      }
	      Adjustment_Data[nb].Reset_Intersect_Flags(i+1,j);
	      Adjustment_Data[nb].Reset_Intersect_Flags(i,j+1);

	    } else if (Adjustment_Data[nb].intersect_west[i+1][j] &&
		       Adjustment_Data[nb].intersect_east[i][j+1]) {
	      if (abs(Local_SolnBlk[nb].Grid.Node[i+1][j].X - Local_SolnBlk[nb].Grid.Node[i][j].X) <
		  abs(Local_SolnBlk[nb].Grid.Node[i][j+1].X - Local_SolnBlk[nb].Grid.Node[i+1][j+1].X)) {
		// Adjust the SE node west.
		Xp = Local_SolnBlk[nb].Grid.Node[i+1][j].X;
		Local_SolnBlk[nb].Grid.Node[i+1][j].X = Local_SolnBlk[nb].Grid.Node[i][j].X;
 		if (Check_Cell_Areas(nb,i,j)) {
		  // An invalid cell has been created. Adjust the NW node east.
		  Local_SolnBlk[nb].Grid.Node[i+1][j].X = Xp;
 		  Local_SolnBlk[nb].Grid.Node[i][j+1].X = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
		  if (Check_Cell_Areas(nb,i,j)) return 2110;
		  Mesh[nb].node_status[i][j+1] = NODE_STATUS_ALIGNED;
		} else {
		  Mesh[nb].node_status[i+1][j] = NODE_STATUS_ALIGNED;
		}
	      } else if (abs(Local_SolnBlk[nb].Grid.Node[i][j+1].X - Local_SolnBlk[nb].Grid.Node[i+1][j+1].X) <
			 abs(Local_SolnBlk[nb].Grid.Node[i+1][j].X - Local_SolnBlk[nb].Grid.Node[i][j].X)) {
		// Adjust the NW node east.
		Xp = Local_SolnBlk[nb].Grid.Node[i][j+1].X;
		Local_SolnBlk[nb].Grid.Node[i][j+1].X = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
 		if (Check_Cell_Areas(nb,i,j)) {
		  // An invalid cell has been created. Adjust the SE node west.
		  Local_SolnBlk[nb].Grid.Node[i][j+1].X = Xp;
		  Local_SolnBlk[nb].Grid.Node[i+1][j].X = Local_SolnBlk[nb].Grid.Node[i][j].X;
		  if (Check_Cell_Areas(nb,i,j)) return 2111;
		  Mesh[nb].node_status[i+1][j] = NODE_STATUS_ALIGNED;
		} else {
		  Mesh[nb].node_status[i][j+1] = NODE_STATUS_ALIGNED;
		}
	      } else {
		// Equal distance... adjust inactive node.
		return 2112;
	      }
	      Adjustment_Data[nb].Reset_Intersect_Flags(i+1,j);
	      Adjustment_Data[nb].Reset_Intersect_Flags(i,j+1);

	    } else {
	      // Not possible.  Return an error flag.
	      return 2113;

	    }

	  } else if (Adjustment_Data[nb].intersect_north[i+1][j] ||
		     Adjustment_Data[nb].intersect_west[i+1][j]) {
	    // Adjust the SE node.
	    if (Adjustment_Data[nb].intersect_north[i+1][j]) {
	      // Adjust the SE node north.
	      Local_SolnBlk[nb].Grid.Node[i+1][j].X = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
 	      if (Check_Cell_Areas(nb,i,j)) {
		// An invalid cell has been created. Adjust the SE node west.
		Local_SolnBlk[nb].Grid.Node[i+1][j].X = Local_SolnBlk[nb].Grid.Node[i][j].X;
		if (Check_Cell_Areas(nb,i,j)) return 2114;
	      }
	    } else if (Adjustment_Data[nb].intersect_west[i+1][j]) {
	      // Adjust the SE node west.
	      Local_SolnBlk[nb].Grid.Node[i+1][j].X = Local_SolnBlk[nb].Grid.Node[i][j].X;
 	      if (Check_Cell_Areas(nb,i,j)) {
		// An invalid cell has been created. Adjust the SE node north.
		Local_SolnBlk[nb].Grid.Node[i+1][j].X = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
		if (Check_Cell_Areas(nb,i,j)) return 2115;
	      }
	    } else {
	      // Not possible.  Return an error flag.
	      return 2116;
	    }
	    // Set node status flag.
  	    Mesh[nb].node_status[i+1][j] = NODE_STATUS_ALIGNED;
	    Adjustment_Data[nb].Reset_Intersect_Flags(i+1,j);

	  } else if (Adjustment_Data[nb].intersect_south[i][j+1] ||
		     Adjustment_Data[nb].intersect_east[i][j+1]) {
	    // Adjust the NW node.
	    if (Adjustment_Data[nb].intersect_south[i][j+1]) {
	      // Adjust the NW node south.
	      Local_SolnBlk[nb].Grid.Node[i][j+1].X = Local_SolnBlk[nb].Grid.Node[i][j].X;
 	      if (Check_Cell_Areas(nb,i,j)) {
		// An invalid cell has been created. Adjust the NW node east.
		Local_SolnBlk[nb].Grid.Node[i][j+1].X = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
		if (Check_Cell_Areas(nb,i,j)) return 2117;
	      }
	    } else if (Adjustment_Data[nb].intersect_east[i][j+1]) {
	      // Adjust the NW node east.
	      Local_SolnBlk[nb].Grid.Node[i][j+1].X = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
 	      if (Check_Cell_Areas(nb,i,j)) {
		// An invalid cell has been created. Adjust the NW node south.
		Local_SolnBlk[nb].Grid.Node[i][j+1].X = Local_SolnBlk[nb].Grid.Node[i][j].X;
		if (Check_Cell_Areas(nb,i,j)) {
		  cout << endl << Local_SolnBlk[nb].Grid.Node[i][j].X;
		  cout << endl << Local_SolnBlk[nb].Grid.Node[i+1][j].X;
		  cout << endl << Local_SolnBlk[nb].Grid.Node[i][j+1].X;
		  cout << endl << Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
		  return 2118;
		}
	      }
	    } else {
	      // Not possible.  Return an error flag.
	      return 2119;
	    }
	    // Set node status flag.
  	    Mesh[nb].node_status[i][j+1] = NODE_STATUS_ALIGNED;
	    Adjustment_Data[nb].Reset_Intersect_Flags(i,j+1);

	  } else {
	    // Not possible.  Return an error flag.
 	    //return 2120;

	  }

	} else if (Mesh[nb].node_status[i  ][j+1] == NODE_STATUS_ALIGNED &&
		   Mesh[nb].node_status[i+1][j+1] != NODE_STATUS_ALIGNED &&
		   Mesh[nb].node_status[i  ][j  ] == NODE_STATUS_ALIGNED &&
		   Mesh[nb].node_status[i+1][j  ] == NODE_STATUS_ALIGNED &&
		   Local_SolnBlk[nb].Grid.lfaceS(i,j) > NANO &&
		   Local_SolnBlk[nb].Grid.lfaceW(i,j) > NANO) {
	  // NE split cell.
	  if (Adjustment_Data[nb].intersect_south[i+1][j+1] ||
	      Adjustment_Data[nb].intersect_west[i+1][j+1]) {
	    // The NE node is the only adjustment candidate.
	    if (Adjustment_Data[nb].intersect_south[i+1][j+1]) {
	      // Adjust the NE node south.
	      Local_SolnBlk[nb].Grid.Node[i+1][j+1].X = Local_SolnBlk[nb].Grid.Node[i+1][j].X;
 	      if (Check_Cell_Areas(nb,i,j)) {
		// An invalid cell has been created. Adjust the NE node west.
		Local_SolnBlk[nb].Grid.Node[i+1][j+1].X = Local_SolnBlk[nb].Grid.Node[i][j+1].X;
		if (Check_Cell_Areas(nb,i,j)) return 3101;
	      }
	    } else if (Adjustment_Data[nb].intersect_west[i+1][j+1]) {
	      // Adjust the NE node west.
	      Local_SolnBlk[nb].Grid.Node[i+1][j+1].X = Local_SolnBlk[nb].Grid.Node[i][j+1].X;
 	      if (Check_Cell_Areas(nb,i,j)) {
		// An invalid cell has been created. Adjust the NE node south.
		Local_SolnBlk[nb].Grid.Node[i+1][j+1].X = Local_SolnBlk[nb].Grid.Node[i+1][j].X;
		if (Check_Cell_Areas(nb,i,j)) return 3102;
	      }
	    } else {
	      // Not possible.  Return an error flag.
	      return 3103;
	    }
	    // Set node status flag.
  	    Mesh[nb].node_status[i+1][j+1] = NODE_STATUS_ALIGNED;
	    Adjustment_Data[nb].Reset_Intersect_Flags(i+1,j+1);
	  }

	} else if (Mesh[nb].node_status[i  ][j+1] != NODE_STATUS_ALIGNED &&
		   Mesh[nb].node_status[i+1][j+1] == NODE_STATUS_ALIGNED &&
		   Mesh[nb].node_status[i  ][j  ] == NODE_STATUS_ALIGNED &&
		   Mesh[nb].node_status[i+1][j  ] == NODE_STATUS_ALIGNED &&
		   Local_SolnBlk[nb].Grid.lfaceS(i,j) > NANO &&
		   Local_SolnBlk[nb].Grid.lfaceE(i,j) > NANO) {
	  // NW split cell.
	  if (Adjustment_Data[nb].intersect_south[i][j+1] ||
	      Adjustment_Data[nb].intersect_east[i][j+1]) {
	    // Adjust the NW node.
	    if (Adjustment_Data[nb].intersect_south[i][j+1]) {
	      // Adjust the NW node south.
	      Local_SolnBlk[nb].Grid.Node[i][j+1].X = Local_SolnBlk[nb].Grid.Node[i][j].X;
 	      if (Check_Cell_Areas(nb,i,j)) {
		// An invalid cell has been created. Adjust the NW node east.
		Local_SolnBlk[nb].Grid.Node[i][j+1].X = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
		if (Check_Cell_Areas(nb,i,j)) return 4101;
	      }
	    } else if (Adjustment_Data[nb].intersect_east[i][j+1]) {
	      // Adjust the NW node east.
	      Local_SolnBlk[nb].Grid.Node[i][j+1].X = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
 	      if (Check_Cell_Areas(nb,i,j)) {
		// An invalid cell has been created. Adjust the NW node south.
		Local_SolnBlk[nb].Grid.Node[i][j+1].X = Local_SolnBlk[nb].Grid.Node[i][j].X;
		if (Check_Cell_Areas(nb,i,j)) return 4102;
	      }
	    } else {
	      // Not possible.  Return an error flag.
	      return 4103;
	    }
	    // Set node status flag.
  	    Mesh[nb].node_status[i][j+1] = NODE_STATUS_ALIGNED;
	    Adjustment_Data[nb].Reset_Intersect_Flags(i,j+1);
	  }

	} else if (Mesh[nb].node_status[i  ][j+1] == NODE_STATUS_ALIGNED &&
		   Mesh[nb].node_status[i+1][j+1] == NODE_STATUS_ALIGNED &&
		   Mesh[nb].node_status[i  ][j  ] == NODE_STATUS_ALIGNED &&
		   Mesh[nb].node_status[i+1][j  ] != NODE_STATUS_ALIGNED &&
		   Local_SolnBlk[nb].Grid.lfaceN(i,j) > NANO &&
		   Local_SolnBlk[nb].Grid.lfaceW(i,j) > NANO) {
	  // SE split cell.
	  if (Adjustment_Data[nb].intersect_north[i+1][j] ||
	      Adjustment_Data[nb].intersect_west[i+1][j]) {
	    // Adjust the SE node.
	    if (Adjustment_Data[nb].intersect_north[i+1][j]) {
	      // Adjust the SE node north.
	      Local_SolnBlk[nb].Grid.Node[i+1][j].X = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
 	      if (Check_Cell_Areas(nb,i,j)) {
		// An invalid cell has been created. Adjust the SE node west.
		Local_SolnBlk[nb].Grid.Node[i+1][j].X = Local_SolnBlk[nb].Grid.Node[i][j].X;
		if (Check_Cell_Areas(nb,i,j)) return 5101;
	      }
	    } else if (Adjustment_Data[nb].intersect_west[i+1][j]) {
	      // Adjust the SE node west.
	      Local_SolnBlk[nb].Grid.Node[i+1][j].X = Local_SolnBlk[nb].Grid.Node[i][j].X;
 	      if (Check_Cell_Areas(nb,i,j)) {
		// An invalid cell has been created. Adjust the SE node north.
		Local_SolnBlk[nb].Grid.Node[i+1][j].X = Local_SolnBlk[nb].Grid.Node[i+1][j+1].X;
		if (Check_Cell_Areas(nb,i,j)) return 5102;
	      }
	    } else {
	      // Not possible.  Return an error flag.
	      return 5103;
	    }
	    // Set node status flag.
  	    Mesh[nb].node_status[i+1][j] = NODE_STATUS_ALIGNED;
	    Adjustment_Data[nb].Reset_Intersect_Flags(i+1,j);
	  }

	} else if (Mesh[nb].node_status[i  ][j+1] == NODE_STATUS_ALIGNED &&
		   Mesh[nb].node_status[i+1][j+1] == NODE_STATUS_ALIGNED &&
		   Mesh[nb].node_status[i  ][j  ] != NODE_STATUS_ALIGNED &&
		   Mesh[nb].node_status[i+1][j  ] == NODE_STATUS_ALIGNED &&
		   Local_SolnBlk[nb].Grid.lfaceN(i,j) > NANO &&
		   Local_SolnBlk[nb].Grid.lfaceE(i,j) > NANO) {
	  // SW split cell.
	  if (Adjustment_Data[nb].intersect_north[i][j] ||
	      Adjustment_Data[nb].intersect_east[i][j]) {
	    // Adjust the SW node.
	    if (Adjustment_Data[nb].intersect_north[i][j]) {
	      // Adjust the SW node north.
	      Xp = Local_SolnBlk[nb].Grid.Node[i][j].X;
	      Local_SolnBlk[nb].Grid.Node[i][j].X = Local_SolnBlk[nb].Grid.Node[i][j+1].X;
 	      if (Check_Cell_Areas(nb,i,j)) {
		// An invalid cell has been created. Adjust the SW node east.
		Local_SolnBlk[nb].Grid.Node[i][j].X = Local_SolnBlk[nb].Grid.Node[i+1][j].X;
		if (Check_Cell_Areas(nb,i,j)) {
		  Local_SolnBlk[nb].Grid.Node[i][j].X = Xp;
		  //return 6101;
		}
	      }
	    } else if (Adjustment_Data[nb].intersect_east[i][j]) {
	      // Adjust the SW node east.
	      Xp = Local_SolnBlk[nb].Grid.Node[i][j].X;
	      Local_SolnBlk[nb].Grid.Node[i][j].X = Local_SolnBlk[nb].Grid.Node[i+1][j].X;
	      if (Check_Cell_Areas(nb,i,j)) {
		// An invalid cell has been created. Adjust the SW node north.
		Local_SolnBlk[nb].Grid.Node[i][j].X = Local_SolnBlk[nb].Grid.Node[i][j+1].X;
		if (Check_Cell_Areas(nb,i,j)) {
		  Local_SolnBlk[nb].Grid.Node[i][j].X = Xp;
		  return 6102;
		}
	      }
	    } else {
	      // Not possible.  Return an error flag.
	      return 6103;
	    }
	    // Set node status flag.
 	    Mesh[nb].node_status[i][j] = NODE_STATUS_ALIGNED;
	    Adjustment_Data[nb].Reset_Intersect_Flags(i,j);
 	  }

   	}

      }

    }
  }

#ifdef _EB_PARALLEL_DEBUG_
//   dout << endl << "AFTER THIRD ADJUSTMENT::Block Number " << Local_Solution_Block_List->Block[nb].gblknum << endl;
//   for (int j = Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost; j >= Local_SolnBlk[nb].JCl-Local_SolnBlk[nb].Nghost; j--) {
//     for (int i = Local_SolnBlk[nb].ICl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].ICu+Local_SolnBlk[nb].Nghost; i++) {
//       if (Mesh[nb].cell_status[i][j] == CELL_STATUS_ACTIVE) dout << " +";
//       else if (Mesh[nb].cell_status[i][j] == CELL_STATUS_INACTIVE) dout << " -";
//       else dout << " x";
//     }
//     dout << endl;
//   }
//   dout << endl;
//   for (int j = Local_SolnBlk[nb].Grid.JNu+Local_SolnBlk[nb].Nghost; j >= Local_SolnBlk[nb].Grid.JNl-Local_SolnBlk[nb].Nghost; j--) {
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       dout << Mesh[nb].node_status[i][j]+1;
//     }
//     dout << " ";
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       dout << Adjustment_Data[nb].intersect_north[i][j];
//     }
//     dout << " ";
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       dout << Adjustment_Data[nb].intersect_south[i][j];
//     }
//     dout << " ";
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       dout << Adjustment_Data[nb].intersect_east[i][j];
//     }
//     dout << " ";
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       dout << Adjustment_Data[nb].intersect_west[i][j];
//     }
//     dout << endl;
//   }
//   dout << endl;
#endif

  // Third mesh adjustment successful.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Mesh_Adjustment_Modify_Boundary_Splines --   *
 *                                                                    *
 * This routine modifies the splines defined on the boundaries of the *
 * solution block to account for cases in which the corner nodes are  *
 * moved.  Important for the Set_BCs and Fix_Refined_Block_Boundaries *
 * routines.                                                          *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Mesh_Adjustment_Modify_Boundary_Splines(const int &nb) {

  // Do not adjust if no interfaces have been defined.
  if (!Interface_Union_List.Ni) return 0;

  // Do not perform adjustment algorithm on the current solution block 
  // if there are no inactive or unknowns cells.
  if (!Adjustment_Data[nb].Interface_Present[0]) return 0;

  int ii, jj;

  // North spline:
  if (Local_SolnBlk[nb].Grid.BCtypeN[0] == BC_NONE &&
      Local_SolnBlk[nb].Grid.BndNorthSpline.np > 1) {
    Local_SolnBlk[nb].Grid.BndNorthSpline.Xp[0] =
      Local_SolnBlk[nb].Grid.Node[Local_SolnBlk[nb].Grid.INl][Local_SolnBlk[nb].Grid.JNu].X;
    Local_SolnBlk[nb].Grid.BndNorthSpline.Xp[Local_SolnBlk[nb].Grid.BndNorthSpline.np-1] =
      Local_SolnBlk[nb].Grid.Node[Local_SolnBlk[nb].Grid.INu][Local_SolnBlk[nb].Grid.JNu].X;
  }

  // South spline:
  if (Local_SolnBlk[nb].Grid.BCtypeS[0] == BC_NONE &&
      Local_SolnBlk[nb].Grid.BndSouthSpline.np > 1) {
    Local_SolnBlk[nb].Grid.BndSouthSpline.Xp[0] =
      Local_SolnBlk[nb].Grid.Node[Local_SolnBlk[nb].Grid.INl][Local_SolnBlk[nb].Grid.JNl].X;
    Local_SolnBlk[nb].Grid.BndSouthSpline.Xp[Local_SolnBlk[nb].Grid.BndSouthSpline.np-1] =
      Local_SolnBlk[nb].Grid.Node[Local_SolnBlk[nb].Grid.INu][Local_SolnBlk[nb].Grid.JNl].X;
  }

  // East spline:
  if (Local_SolnBlk[nb].Grid.BCtypeE[0] == BC_NONE &&
      Local_SolnBlk[nb].Grid.BndEastSpline.np > 1) {
    Local_SolnBlk[nb].Grid.BndEastSpline.Xp[0] =
      Local_SolnBlk[nb].Grid.Node[Local_SolnBlk[nb].Grid.INu][Local_SolnBlk[nb].Grid.JNl].X;
    Local_SolnBlk[nb].Grid.BndEastSpline.Xp[Local_SolnBlk[nb].Grid.BndEastSpline.np-1] =
      Local_SolnBlk[nb].Grid.Node[Local_SolnBlk[nb].Grid.INu][Local_SolnBlk[nb].Grid.JNu].X;
  }

  // West spline:
  if (Local_SolnBlk[nb].Grid.BCtypeW[0] == BC_NONE &&
      Local_SolnBlk[nb].Grid.BndWestSpline.np > 1) {
    Local_SolnBlk[nb].Grid.BndWestSpline.Xp[0] =
      Local_SolnBlk[nb].Grid.Node[Local_SolnBlk[nb].Grid.INl][Local_SolnBlk[nb].Grid.JNl].X;
    Local_SolnBlk[nb].Grid.BndWestSpline.Xp[Local_SolnBlk[nb].Grid.BndWestSpline.np-1] =
      Local_SolnBlk[nb].Grid.Node[Local_SolnBlk[nb].Grid.INl][Local_SolnBlk[nb].Grid.JNu].X;
  }

  // Modify boundary splines successful.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Mesh_Adjustment_Finalize --                  *
 *                                                                    *
 * This routine finalizes the mesh adjustment by updating the areas   *
 * and centroids of the cells.                                        *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Mesh_Adjustment_Finalize(const int &nb) {

  // Declare local variables.
  int error_flag;

  int *BC_N(NULL), *BC_S(NULL), *BC_E(NULL), *BC_W(NULL);

  // Special case for Ringleb's Flow and Flat Plate.
  //    -  Jai says this greatly improves the quality
  //       of ghost cells near an embedded boundary.
  if (Interface_Union_List[1].Type == INTERFACE_RINGLEB) {
    for (int i = 0; i < Local_SolnBlk[nb].Grid.NCi; i++) {
      Local_SolnBlk[nb].Grid.BCtypeS[i] = BC_REFLECTION;
    }
  } else if (Interface_Union_List[1].Type == INTERFACE_FLAT_PLATE) {
    BC_N = new int[Local_SolnBlk[nb].Grid.NCi];
    BC_S = new int[Local_SolnBlk[nb].Grid.NCi];
    BC_E = new int[Local_SolnBlk[nb].Grid.NCj];
    BC_W = new int[Local_SolnBlk[nb].Grid.NCj];
    for (int i = 0; i < Local_SolnBlk[nb].Grid.NCi; i++) {
      BC_N[i] = Local_SolnBlk[nb].Grid.BCtypeN[i];
      BC_S[i] = Local_SolnBlk[nb].Grid.BCtypeS[i];
      Local_SolnBlk[nb].Grid.BCtypeN[i] = BC_CONSTANT_EXTRAPOLATION;
      Local_SolnBlk[nb].Grid.BCtypeS[i] = BC_CONSTANT_EXTRAPOLATION;
    }
    for (int j = 0; j < Local_SolnBlk[nb].Grid.NCj; j++) {
      BC_E[j] = Local_SolnBlk[nb].Grid.BCtypeE[j];
      BC_W[j] = Local_SolnBlk[nb].Grid.BCtypeW[j];
      Local_SolnBlk[nb].Grid.BCtypeE[j] = BC_CONSTANT_EXTRAPOLATION;
      Local_SolnBlk[nb].Grid.BCtypeW[j] = BC_CONSTANT_EXTRAPOLATION;
    }
  }

   // Compute the exterior nodes for the quadrilateral mesh block.
   //Update_Exterior_Nodes(Local_SolnBlk[nb].Grid);
   //Update_Exterior_Nodes(nb);

  // Compute the cells for the quadrilateral mesh block.
  Update_Cells(Local_SolnBlk[nb].Grid);

  // Update the exterior cell status and types.
  Update_Exterior_Cells(nb);

  // Compute the cells for the quadrilateral mesh block.
  //Update_Cells(Local_SolnBlk[nb].Grid);

  // Check to see if all cell areas are positive.
  error_flag = Check_Quad_Block(Local_SolnBlk[nb].Grid);
  if (error_flag) return error_flag;

  // Special case for Ringleb's Flow.
  if (Interface_Union_List[1].Type == INTERFACE_RINGLEB) {
    for (int i = 0; i < Local_SolnBlk[nb].Grid.NCi; i++) {
      Local_SolnBlk[nb].Grid.BCtypeS[i] = BC_RINGLEB_FLOW;
    }
  } else if (Interface_Union_List[1].Type == INTERFACE_FLAT_PLATE) {
    for (int i = 0; i < Local_SolnBlk[nb].Grid.NCi; i++) {
      Local_SolnBlk[nb].Grid.BCtypeN[i] = BC_N[i];
      Local_SolnBlk[nb].Grid.BCtypeS[i] = BC_S[i];
    }
    for (int j = 0; j < Local_SolnBlk[nb].Grid.NCj; j++) {
      Local_SolnBlk[nb].Grid.BCtypeE[j] = BC_E[j];
      Local_SolnBlk[nb].Grid.BCtypeW[j] = BC_W[j];
    }
  }

  //delete arrays if needed
  if(BC_N != NULL) {delete [] BC_N;}
  if(BC_S != NULL) {delete [] BC_S;}
  if(BC_E != NULL) {delete [] BC_E;}
  if(BC_W != NULL) {delete [] BC_W;}


  // Mesh adjustment finalized.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Pre_Mesh_Adjustment_Painting --              *
 *                                                                    *
 * This routine performs the pre-mesh painting for the specified      *
 * solution block.  The cells of the mesh are tagged as active        *
 * (external), inactive (internal), or unknown.  The unknown cells    *
 * are candidates for adjustment.  The nodes of the cells are tagged  *
 * as known or unknown.  Note that all cells and nodes are            *
 * initialized as known.  The unknown nodes (of an unknown cell) are  *
 * candidates for adjustment.  Most of the active cells (and nodes)   *
 * are quickly identified by comparing each of the nodes of a cell    *
 * with bounding boxes of each of the interface components.  If all   *
 * of the nodes are deemed outside the bounding box then the cell is  *
 * deemed active.  Otherwise, each edge of the cell is compared to    *
 * each edge of the list of interface unions.  If an intersection     *
 * exists between and cell edge and an interface edge then the nodes  *
 * of the cell edge and the cell itself are tagged as unknown. If no  *
 * edge/interface intersections exist, ray-tracing is performed to    *
 * determine if the cell is deemed active or in-active.  The ray-     *
 * tracing algorithm is performed by counting the number of           *
 * intersections between the line composed of the cell centroid and a *
 * reference point within the embedded boundary (typically the        *
 * centroid) and each edge of the embedded boundary.  An odd number   *
 * of intersections indicates that the cell is outside the interface  *
 * (active) and an even number of intersections indicates that the    *
 * cell is inside the interface (in-active).                          *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Pre_Mesh_Adjustment_Painting(const int &nb) {

  int iSW, iSE, iNE, iNW, iN, iS, iE, iW, nu, Number_of_Unknown_Cells;

  // Initialize the number of inactive cells counter.
  Adjustment_Data[nb].Number_of_Inactive_Cells = 0;
  Number_of_Unknown_Cells = 0;
  Adjustment_Data[nb].Reset_Interface_Present();

  // Conduct cell painting algorithm.
  for (int j = Local_SolnBlk[nb].JCl-Local_SolnBlk[nb].Nghost; j <= Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost; j++) {
    for (int i = Local_SolnBlk[nb].ICl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].ICu+Local_SolnBlk[nb].Nghost; i++) {

      for (int ni = 1; ni <= Interface_Component_List.Ni; ni++) {

	// To avoid re-painting cells, only attempt to paints cells if
	// they are currently active.
	if (Mesh[nb].cell_status[i][j] == CELL_STATUS_ACTIVE) {

	  // Conduct the point-in-bounding-box test for each node of the cell.
 	  iSW = Interface_Component_List[ni].Point_In_Bounding_Box(Local_SolnBlk[nb].Grid.nodeSW(i,j).X);
 	  iSE = Interface_Component_List[ni].Point_In_Bounding_Box(Local_SolnBlk[nb].Grid.nodeSE(i,j).X);
 	  iNE = Interface_Component_List[ni].Point_In_Bounding_Box(Local_SolnBlk[nb].Grid.nodeNE(i,j).X);
 	  iNW = Interface_Component_List[ni].Point_In_Bounding_Box(Local_SolnBlk[nb].Grid.nodeNW(i,j).X);

 	  if (!iSW && !iSE && !iNE && !iNW) {
	    // The current cell is not inside the bounding box of the
	    // current interface.  Cell remains active.
 	    //Mesh[nb].cell_status[i][j] = CELL_STATUS_ACTIVE;

	  } else {

	    // Determine which interface of the union list the current
	    // componment is associated with.
	    nu = Interface_Component_List.relations[ni-1][0];

	    // Determine any of the edges of the current cell intersect 
	    // with any of the edges of the current interface.
	    iN = OFF; iS = OFF; iE = OFF; iW = OFF;
	    for (int np = 0; np < Interface_Union_List[nu].Spline.np-1; np++) {
	      if (!iN) iN = Line_Intersection(Interface_Union_List[nu].Spline.Xp[np],
					      Interface_Union_List[nu].Spline.Xp[np+1],
					      Local_SolnBlk[nb].Grid.nodeNW(i,j).X,
					      Local_SolnBlk[nb].Grid.nodeNE(i,j).X);
	      if (!iS) iS = Line_Intersection(Interface_Union_List[nu].Spline.Xp[np],
					      Interface_Union_List[nu].Spline.Xp[np+1],
					      Local_SolnBlk[nb].Grid.nodeSW(i,j).X,
					      Local_SolnBlk[nb].Grid.nodeSE(i,j).X);
	      if (!iE) iE = Line_Intersection(Interface_Union_List[nu].Spline.Xp[np],
					      Interface_Union_List[nu].Spline.Xp[np+1],
					      Local_SolnBlk[nb].Grid.nodeSE(i,j).X,
					      Local_SolnBlk[nb].Grid.nodeNE(i,j).X);
	      if (!iW) iW = Line_Intersection(Interface_Union_List[nu].Spline.Xp[np],
					      Interface_Union_List[nu].Spline.Xp[np+1],
					      Local_SolnBlk[nb].Grid.nodeNW(i,j).X,
					      Local_SolnBlk[nb].Grid.nodeSW(i,j).X);
	      // If an interface edge intersects with one of the cell
	      // edges then exit immediately.
	      //if (iN || iS || iE || iW) break;
	    }
	    // Determine the status of the cell.
	    if (iN || iS || iE || iW) {
	      // The status of the cell not known.  Target for adjustment.
	      Mesh[nb].cell_status[i][j] = CELL_STATUS_UNKNOWN;
	      Adjustment_Data[nb].Interface_Present[nu] = ON;
	      Adjustment_Data[nb].Number_of_Inactive_Cells++;
	      Number_of_Unknown_Cells++;
	      if (iS || iW) Mesh[nb].node_status[i  ][j  ] = NODE_STATUS_UNKNOWN;
	      if (iS || iE) Mesh[nb].node_status[i+1][j  ] = NODE_STATUS_UNKNOWN;
	      if (iN || iW) Mesh[nb].node_status[i  ][j+1] = NODE_STATUS_UNKNOWN;
	      if (iN || iE) Mesh[nb].node_status[i+1][j+1] = NODE_STATUS_UNKNOWN;
	    } else {
 	      if (Interface_Union_List[nu].Point_In_Interface(Local_SolnBlk[nb].Grid.Cell[i][j].Xc)) {
		// The cell is inactive.
		Mesh[nb].cell_status[i][j] = nu;
		Adjustment_Data[nb].Interface_Present[nu] = ON;
		// Increment the number of inactive cells (include ghost cells).
		Adjustment_Data[nb].Number_of_Inactive_Cells++;
	      }
	    }

	  }

	}

      }

    }
  }

  // Special case for Ringleb's Flow.
  if (Interface_Union_List[1].Type == INTERFACE_RINGLEB) {
    Adjustment_Data[nb].Number_of_Inactive_Cells = 0;
    Number_of_Unknown_Cells = 0;
    for (int j = Local_SolnBlk[nb].JCl-Local_SolnBlk[nb].Nghost; j <= Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost; j++) {
      for (int i = Local_SolnBlk[nb].ICl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].ICu+Local_SolnBlk[nb].Nghost; i++) {
 	if (Mesh[nb].cell_status[i][j] == CELL_STATUS_ACTIVE) {
 	  Mesh[nb].cell_status[i][j] = CELL_STATUS_INACTIVE;
 	} else if (Mesh[nb].cell_status[i][j] == CELL_STATUS_INACTIVE) {
 	  Mesh[nb].cell_status[i][j] = CELL_STATUS_ACTIVE;
	  Adjustment_Data[nb].Number_of_Inactive_Cells++;
	} else {
	  Adjustment_Data[nb].Number_of_Inactive_Cells++;
	  Number_of_Unknown_Cells++;
	}
      }
    }
  }

  // Set interface present flag.
  if (Number_of_Unknown_Cells) {
    Adjustment_Data[nb].Interface_Present[0] = ON;
  } else {
    Adjustment_Data[nb].Interface_Present[0] = OFF;
    for (int ni = 1; ni < Adjustment_Data[nb].Ni+1; ni++)
      Adjustment_Data[nb].Interface_Present[ni] = OFF;
  }

#ifdef _EB_PARALLEL_DEBUG_
//   dout << endl << "PRE-PAINTING::Block Number " << Local_Solution_Block_List->Block[nb].gblknum
//        << ", Number of Inactive Cells = " << Adjustment_Data[nb].Number_of_Inactive_Cells
//        << ", Interface Present = " << Adjustment_Data[nb].Interface_Present[0]
//        << endl;
//   for (int j = Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost; j >= Local_SolnBlk[nb].JCl-Local_SolnBlk[nb].Nghost; j--) {
//     for (int i = Local_SolnBlk[nb].ICl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].ICu+Local_SolnBlk[nb].Nghost; i++) {
//       dout << Mesh[nb].cell_status[i][j]+1;
//     }
//     dout << endl;
//   }
//   dout << endl;
//   for (int j = Local_SolnBlk[nb].Grid.JNu+Local_SolnBlk[nb].Nghost; j >= Local_SolnBlk[nb].Grid.JNl-Local_SolnBlk[nb].Nghost; j--) {
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       dout << Mesh[nb].node_status[i][j]+1;
//     }
//     dout << endl;
//   }
//   dout << endl;
#endif

  // Pre mesh adjustment cell and node painting conducted successfully.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Post_Mesh_Adjustment_Painting --             *
 *                                                                    *
 * This routine performs the post-mesh painting on the specified      *
 * solution block.  Here cells previously marked as unknown are re-   *
 * tagged.   The ray-tracing algorithm used in the pre-painting step  *
 * is used to determine if the adjusted cell is inside (active) or    *
 * outside (in-active) the embedded boundary.  All inactive cells are *
 * tagged using the number of the interface of union that it is       *
 * associated with (one to the number of union interfaces).  This     *
 * tagging method allows for quick identification of the associated   *
 * embedded boundary.  All active cells are tagged as zero.           *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Post_Mesh_Adjustment_Painting(const int &nb) {

  int status, Number_of_Unknown_Cells, ii, jj;

  // Initialize the number of inactive cells counter.
  Adjustment_Data[nb].Number_of_Inactive_Cells = 0;
  Number_of_Unknown_Cells = 0;

  // Conduct cell painting algorithm.
  for (int j = Local_SolnBlk[nb].JCl-Local_SolnBlk[nb].Nghost; j <= Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost; j++) {
    for (int i = Local_SolnBlk[nb].ICl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].ICu+Local_SolnBlk[nb].Nghost; i++) {
      if (Mesh[nb].cell_status[i][j] == CELL_STATUS_UNKNOWN) Number_of_Unknown_Cells++;
      Mesh[nb].cell_status[i][j] = CELL_STATUS_ACTIVE;
      for (int ni = 1; ni <= Interface_Union_List.Ni; ni++) {
// 	// Conduct point-in-interface test at the cell-centre.
// 	if (!Interface_Union_List[ni].Point_In_Interface(Local_SolnBlk[nb].Grid.Cell[i][j].Xc)) {
// 	  status = CELL_STATUS_ACTIVE;
// 	} else {
// 	  status = ni;
// 	}
	if (Mesh[nb].node_status[i][j] != NODE_STATUS_ALIGNED) { ii = i; jj = j; }
	else if (Mesh[nb].node_status[i+1][j] != NODE_STATUS_ALIGNED) { ii = i+1; jj = j; }
	else if (Mesh[nb].node_status[i][j+1] != NODE_STATUS_ALIGNED) { ii = i; jj = j+1; }
	else if (Mesh[nb].node_status[i+1][j+1] != NODE_STATUS_ALIGNED) { ii = i+1; jj = j+1; }
 	if (!Interface_Union_List[ni].Point_In_Interface(Local_SolnBlk[nb].Grid.Node[ii][jj].X)) {
 	  status = CELL_STATUS_ACTIVE;
 	} else {
 	  status = ni;
 	}
	if (status > Mesh[nb].cell_status[i][j]) {
	  Mesh[nb].cell_status[i][j] = status;
	  Adjustment_Data[nb].Interface_Present[ni] = ON;
	}
      }
      // Increment the number of inactive cells (non-ghost cells only).
      if (Mesh[nb].cell_status[i][j] != CELL_STATUS_ACTIVE &&
	  i >= Local_SolnBlk[nb].ICl && i <= Local_SolnBlk[nb].ICu &&
	  j >= Local_SolnBlk[nb].JCl && j <= Local_SolnBlk[nb].JCu) {
	Adjustment_Data[nb].Number_of_Inactive_Cells++;
      }
    }
  }

  // Tag all cells as either QUADRILATERAL, ADJUSTED, or TRIANGLE.
  for (int j = Local_SolnBlk[nb].JCl-Local_SolnBlk[nb].Nghost; j <= Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost; j++) {
    for (int i = Local_SolnBlk[nb].ICl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].ICu+Local_SolnBlk[nb].Nghost; i++) {
      // Determine if cell has been adjusted by comparing the new
      // centroid with the old centroid.
      if (abs(Local_SolnBlk[nb].Grid.Cell[i][j].Xc - OGrid[nb].Cell[i][j].Xc) < NANO) {
	Mesh[nb].cell_type[i][j] = CELL_TYPE_QUADRILATERAL;
      } else {
	Mesh[nb].cell_type[i][j] = CELL_TYPE_ADJUSTED_QUADRILATERAL;
	// ADJUSTED cells could be TRIANGLES.
	if (Local_SolnBlk[nb].Grid.lfaceN(i,j) < TOLER || Local_SolnBlk[nb].Grid.lfaceS(i,j) < TOLER ||
	    Local_SolnBlk[nb].Grid.lfaceE(i,j) < TOLER || Local_SolnBlk[nb].Grid.lfaceW(i,j) < TOLER) {
	  Mesh[nb].cell_type[i][j] = CELL_TYPE_TRIANGLE;
	}
      }
    }
  }

  // Special case for Ringleb's Flow.
  if (Interface_Union_List[1].Type == INTERFACE_RINGLEB) {
    for (int j = Local_SolnBlk[nb].JCl-Local_SolnBlk[nb].Nghost; j <= Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost; j++) {
      for (int i = Local_SolnBlk[nb].ICl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].ICu+Local_SolnBlk[nb].Nghost; i++) {
 	if (Mesh[nb].cell_status[i][j] == CELL_STATUS_ACTIVE) {
 	  Mesh[nb].cell_status[i][j] = CELL_STATUS_INACTIVE;
 	} else if (Mesh[nb].cell_status[i][j] == CELL_STATUS_INACTIVE) {
 	  Mesh[nb].cell_status[i][j] = CELL_STATUS_ACTIVE;
	}
      }
    }
    Adjustment_Data[nb].Number_of_Inactive_Cells = Adjustment_Data[nb].Number_of_Active_Cells();
  }

  // Update the exterior cell status and types.
  Update_Exterior_Cells(nb);

  // Set interface present flag.
  if (Number_of_Unknown_Cells) {
    Adjustment_Data[nb].Interface_Present[0] = ON;
  } else {
    Adjustment_Data[nb].Interface_Present[0] = OFF;
    for (int ni = 1; ni < Adjustment_Data[nb].Ni+1; ni++)
      Adjustment_Data[nb].Interface_Present[ni] = OFF;
  }

#ifdef _EB_PARALLEL_DEBUG_
//   dout << endl << "POST-PAINTING::Block Number " << Local_Solution_Block_List->Block[nb].gblknum << endl;
//   for (int j = Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost; j >= Local_SolnBlk[nb].JCl-Local_SolnBlk[nb].Nghost; j--) {
//     for (int i = Local_SolnBlk[nb].ICl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].ICu+Local_SolnBlk[nb].Nghost; i++) {
//       dout << Mesh[nb].cell_status[i][j]+1;
//     }
//     dout << endl;
//   }
//   dout << endl << endl;
//   for (int j = Local_SolnBlk[nb].Grid.JNu+Local_SolnBlk[nb].Nghost; j >= Local_SolnBlk[nb].Grid.JNl-Local_SolnBlk[nb].Nghost; j--) {
//     for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
//       dout << Mesh[nb].node_status[i][j]+1;
//     }
//     dout << endl;
//   }
//   dout << endl;
#endif

  // Post mesh adjustment cell and node painting conducted successfully.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Cell_Status_Test --                          *
 *                                                                    *
 * This routines determines if the specified cell or one of the       *
 * neighbouring cells of the specified solution block are unknown.    *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Cell_Status_Test(const int &nb, const int &i, const int &j) {
  if (Mesh[nb].cell_status[i  ][j  ] == CELL_STATUS_UNKNOWN) return 1;
  if (i > Local_SolnBlk[nb].ICl-Local_SolnBlk[nb].Nghost && j > Local_SolnBlk[nb].JCl-Local_SolnBlk[nb].Nghost)
    if (Mesh[nb].cell_status[i-1][j-1] == CELL_STATUS_UNKNOWN) return 1;
  if (j > Local_SolnBlk[nb].JCl-Local_SolnBlk[nb].Nghost)
    if (Mesh[nb].cell_status[i  ][j-1] == CELL_STATUS_UNKNOWN) return 1;
  if (i < Local_SolnBlk[nb].ICu+Local_SolnBlk[nb].Nghost && j > Local_SolnBlk[nb].JCl-Local_SolnBlk[nb].Nghost)
    if (Mesh[nb].cell_status[i+1][j-1] == CELL_STATUS_UNKNOWN) return 1;
  if (i > Local_SolnBlk[nb].ICl-Local_SolnBlk[nb].Nghost)
    if (Mesh[nb].cell_status[i-1][j  ] == CELL_STATUS_UNKNOWN) return 1;
  if (i < Local_SolnBlk[nb].ICu+Local_SolnBlk[nb].Nghost)
    if (Mesh[nb].cell_status[i+1][j  ] == CELL_STATUS_UNKNOWN) return 1;
  if (i > Local_SolnBlk[nb].ICl-Local_SolnBlk[nb].Nghost && j < Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost)
    if (Mesh[nb].cell_status[i-1][j+1] == CELL_STATUS_UNKNOWN) return 1;
  if (j < Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost)
    if (Mesh[nb].cell_status[i  ][j+1] == CELL_STATUS_UNKNOWN) return 1;
  if (i < Local_SolnBlk[nb].ICu+Local_SolnBlk[nb].Nghost && j < Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost)
    if (Mesh[nb].cell_status[i+1][j+1] == CELL_STATUS_UNKNOWN) return 1;
  return 0;
}

/**********************************************************************
 * EmbeddedBoundaries2D::Check_Cell_Areas --                          *
 *                                                                    *
 * This routine checks the area of the specified cell and the areas   *
 * of neighbouring cells of the specified solution block to ensure    *
 * that no zero or negative area cells have been created by the mesh  *
 * adjustment.                                                        *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Check_Cell_Areas(const int &nb, const int &i, const int &j) {
  double rel_tol = fabs(Local_SolnBlk[nb].Grid.area(i,j));
  //ensure no area is less that TOLER times the area of
  //its neighbour
  if (i > Local_SolnBlk[nb].ICl-Local_SolnBlk[nb].Nghost) {
    if (j > Local_SolnBlk[nb].JCl-Local_SolnBlk[nb].Nghost) {
      rel_tol = max(rel_tol,fabs(Local_SolnBlk[nb].Grid.area(i-1,j-1)));
    }
    rel_tol = max(rel_tol,fabs(Local_SolnBlk[nb].Grid.area(i-1,j)));
    if (j < Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost) {
      rel_tol = max(rel_tol,fabs(Local_SolnBlk[nb].Grid.area(i-1,j+1)));
    }
  }
  if (j > Local_SolnBlk[nb].JCl-Local_SolnBlk[nb].Nghost) {
    rel_tol = max(rel_tol,fabs(Local_SolnBlk[nb].Grid.area(i,j-1)));
  }
  if (j < Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost) {
    rel_tol = max(rel_tol,fabs(Local_SolnBlk[nb].Grid.area(i,j+1)));
  }
  if (i < Local_SolnBlk[nb].ICu+Local_SolnBlk[nb].Nghost) {
    if (j > Local_SolnBlk[nb].JCl-Local_SolnBlk[nb].Nghost) {
      rel_tol = max(rel_tol,fabs(Local_SolnBlk[nb].Grid.area(i+1,j-1)));
    }
    rel_tol = max(rel_tol,fabs(Local_SolnBlk[nb].Grid.area(i+1,j)));
    if (j < Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost) {
      rel_tol = max(rel_tol,fabs(Local_SolnBlk[nb].Grid.area(i+1,j+1)));
    }
  }
  rel_tol = TOLER*rel_tol;

  if (i > Local_SolnBlk[nb].ICl-Local_SolnBlk[nb].Nghost) {
    if (j > Local_SolnBlk[nb].JCl-Local_SolnBlk[nb].Nghost) {
      if (Local_SolnBlk[nb].Grid.area(i-1,j-1) < rel_tol) return 1;
    }
    if (Local_SolnBlk[nb].Grid.area(i-1,j) < rel_tol) return 1;
    if (j < Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost) {
      if (Local_SolnBlk[nb].Grid.area(i-1,j+1) < rel_tol) return 1;
    }
  }
  if (j > Local_SolnBlk[nb].JCl-Local_SolnBlk[nb].Nghost) {
    if (Local_SolnBlk[nb].Grid.area(i,j-1) < rel_tol) return 1;
  }
  if (Local_SolnBlk[nb].Grid.area(i,j) < rel_tol) return 1;
  if (j < Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost) {
    if (Local_SolnBlk[nb].Grid.area(i,j+1) < rel_tol) return 1;
  }
  if (i < Local_SolnBlk[nb].ICu+Local_SolnBlk[nb].Nghost) {
    if (j > Local_SolnBlk[nb].JCl-Local_SolnBlk[nb].Nghost) {
      if (Local_SolnBlk[nb].Grid.area(i+1,j-1) < rel_tol) return 1;
    }
    if (Local_SolnBlk[nb].Grid.area(i+1,j) < rel_tol) return 1;
    if (j < Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost) {
      if (Local_SolnBlk[nb].Grid.area(i+1,j+1) < rel_tol) return 1;
    }
  }
  // All of the areas are positive.
  return 0;
}

/**********************************************************************
 * EmbeddedBoundaries2D::Mesh_Adjustment_Statistics --                *
 *                                                                    *
 * This routine determines and outputs the statistics of the mesh     *
 * adjustment for all of the solution blocks involved in the          *
 * computation.                                                       *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
void EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Mesh_Adjustment_Statistics(void) {

  int number_of_adjusted_blocks, total_number_of_blocks;
  int number_of_inactive_cells, total_number_of_cells;
  double cell_area_ratio;

  number_of_adjusted_blocks = 0;
  total_number_of_blocks = QuadTree->countUsedBlocks();
  number_of_inactive_cells = 0;
  total_number_of_cells = 0;
  cell_area_ratio = ONE;
  for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
    if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      Mesh_Adjustment_Statistics(nb,
				 number_of_adjusted_blocks,
				 number_of_inactive_cells,
				 total_number_of_cells,
				 cell_area_ratio);
    }
  }

#ifdef _MPI_VERSION
  number_of_adjusted_blocks = CFFC_Summation_MPI(number_of_adjusted_blocks);
  number_of_inactive_cells = CFFC_Summation_MPI(number_of_inactive_cells);
  total_number_of_cells = CFFC_Summation_MPI(total_number_of_cells);
  cell_area_ratio = CFFC_Minimum_MPI(cell_area_ratio);
#endif

  // Primary processor outputs adjustment statistics.
  if (CFFC_Primary_MPI_Processor()) {
    cout << "\n  -> Mesh Adjustment Statistics: "
	 << "\n     -> Number of Adjusted Blocks = "
	 << number_of_adjusted_blocks << " of " << total_number_of_blocks
	 << " = " << 100.0*double(number_of_adjusted_blocks)/double(total_number_of_blocks) << "%"
	 << "\n     -> Number of Inactive Cells = "
	 << number_of_inactive_cells << " of " << total_number_of_cells
	 << " = " << 100.0*double(number_of_inactive_cells)/double(total_number_of_cells) << "%"
	 << "\n     -> Smallest Cell Area Ratio = " << cell_area_ratio;
  }

}

/**********************************************************************
 * EmbeddedBoundaries2D::Mesh_Adjustment_Statistics --                *
 *                                                                    *
 * This routine determines the statistics of the mesh adjustment for  *
 * the specified solution block.                                      *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
void EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Mesh_Adjustment_Statistics(const int &nb,
			   int &number_of_adjusted_blocks,
			   int &number_inactive_cells,
			   int &total_number_of_cells,
			   double &cell_area_ratio) {

  // Include the number of cells for this block.
  total_number_of_cells += ((Local_SolnBlk[nb].JCu-Local_SolnBlk[nb].JCl+1)*
			    (Local_SolnBlk[nb].ICu-Local_SolnBlk[nb].ICl+1));

  // Exit immediately if no inactive cells exist.
  if (!Adjustment_Data[nb].Number_of_Inactive_Cells) return ;

  double area_ratio;

  // Increment the number of adjusted blocks.
  number_of_adjusted_blocks++;

  // Include the number of inactive cells for this block.
  number_inactive_cells += Adjustment_Data[nb].Number_of_Inactive_Cells;

  // Determine the ratio of the smallest cell to any direct neighbour
  // cell (active cells only).
  for (int j = Local_SolnBlk[nb].JCl; j <= Local_SolnBlk[nb].JCu-1; j++) {
    for (int i = Local_SolnBlk[nb].ICl; i <= Local_SolnBlk[nb].ICu-1; i++) {
      if (Mesh[nb].cell_status[i][j] == CELL_STATUS_ACTIVE) {
	if (Mesh[nb].cell_status[i+1][j] == CELL_STATUS_ACTIVE) {
	  area_ratio = Local_SolnBlk[nb].Grid.Cell[i][j].A/Local_SolnBlk[nb].Grid.Cell[i+1][j].A;
	  cell_area_ratio = min(cell_area_ratio,area_ratio);
	  cell_area_ratio = min(cell_area_ratio,ONE/area_ratio);
	}
	if (Mesh[nb].cell_status[i][j+1] == CELL_STATUS_ACTIVE) {
	  area_ratio = Local_SolnBlk[nb].Grid.Cell[i][j].A/Local_SolnBlk[nb].Grid.Cell[i][j+1].A;
	  cell_area_ratio = min(cell_area_ratio,area_ratio);
	  cell_area_ratio = min(cell_area_ratio,ONE/area_ratio);
	}
      }
    }
  }

}

/**********************************************************************
 * EmbeddedBoundaries2D::Update_Exterior_Cells --                     *
 *                                                                    *
 * This routine sets the cell status and types in the ghost cells at  *
 * boundaries for the specified solution block.                       *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
void EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Update_Exterior_Cells(const int &nb) {

  // Update West and East boundary cells.
  for (int j = Local_SolnBlk[nb].JCl; j <= Local_SolnBlk[nb].JCu; j++) {

    if (Local_SolnBlk[nb].Grid.BCtypeW[j] == Local_SolnBlk[nb].Grid.BCtypeW[j-1] &&
	(Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_REFLECTION ||
	 Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX ||
	 Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
	 Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_MOVING_WALL_HEATFLUX ||
	 Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL ||
	 Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_BURNING_SURFACE)) {
      for (int ng = 1; ng <= Local_SolnBlk[nb].Nghost; ng++) {
	Mesh[nb].cell_status[Local_SolnBlk[nb].ICl-ng][j] = Mesh[nb].cell_status[Local_SolnBlk[nb].ICl+ng-1][j];
	Mesh[nb].cell_type[Local_SolnBlk[nb].ICl-ng][j] = Mesh[nb].cell_type[Local_SolnBlk[nb].ICl+ng-1][j];
      }
    }
    
    if (Local_SolnBlk[nb].Grid.BCtypeE[j] == Local_SolnBlk[nb].Grid.BCtypeE[j-1] &&
	(Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_REFLECTION ||
	 Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX ||
	 Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
	 Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_MOVING_WALL_HEATFLUX ||
	 Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL ||
	 Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_BURNING_SURFACE)) {
      for (int ng = 1; ng <= Local_SolnBlk[nb].Nghost; ng++) {
	Mesh[nb].cell_status[Local_SolnBlk[nb].ICu+ng][j] = Mesh[nb].cell_status[Local_SolnBlk[nb].ICu-ng+1][j];
	Mesh[nb].cell_type[Local_SolnBlk[nb].ICu+ng][j] = Mesh[nb].cell_type[Local_SolnBlk[nb].ICu-ng+1][j];
      }
    }

  }
  
  // Update South and North boundary cells.
  for (int i = Local_SolnBlk[nb].ICl; i <= Local_SolnBlk[nb].ICu; i++) {

    if (Local_SolnBlk[nb].Grid.BCtypeS[i] == Local_SolnBlk[nb].Grid.BCtypeS[i-1] &&
	(Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_REFLECTION || 
	 Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX || 
	 Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL || 
	 Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_MOVING_WALL_HEATFLUX || 
	 Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_MOVING_WALL_ISOTHERMAL || 
	 Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_BURNING_SURFACE)) {
      for (int ng = 1; ng <= Local_SolnBlk[nb].Nghost; ng++) {
	Mesh[nb].cell_status[i][Local_SolnBlk[nb].JCl-ng] = Mesh[nb].cell_status[i][Local_SolnBlk[nb].JCl+ng-1];
	Mesh[nb].cell_type[i][Local_SolnBlk[nb].JCl-ng] = Mesh[nb].cell_type[i][Local_SolnBlk[nb].JCl+ng-1];
      }
    }

    if (Local_SolnBlk[nb].Grid.BCtypeN[i] == Local_SolnBlk[nb].Grid.BCtypeN[i-1] &&
	(Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_REFLECTION ||
	 Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX ||
	 Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
	 Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_MOVING_WALL_HEATFLUX ||
	 Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_MOVING_WALL_ISOTHERMAL ||
	 Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_BURNING_SURFACE)) {
      for (int ng = 1; ng <= Local_SolnBlk[nb].Nghost; ng++) {
	Mesh[nb].cell_status[i][Local_SolnBlk[nb].JCu+ng] = Mesh[nb].cell_status[i][Local_SolnBlk[nb].JCu+ng-1];
	Mesh[nb].cell_type[i][Local_SolnBlk[nb].JCu+ng] = Mesh[nb].cell_type[i][Local_SolnBlk[nb].JCu+ng-1];
      }
    }

  }

  // Update SW and SE corner nodes.
  for (int j = Local_SolnBlk[nb].JCl-Local_SolnBlk[nb].Grid.Nghost; j <= Local_SolnBlk[nb].JCl-1; j++) {

    if ((Local_SolnBlk[nb].Grid.BCtypeS[Local_SolnBlk[nb].Grid.INl-1] == BC_REFLECTION ||
	 Local_SolnBlk[nb].Grid.BCtypeS[Local_SolnBlk[nb].Grid.INl-1] == BC_WALL_VISCOUS_HEATFLUX ||
	 Local_SolnBlk[nb].Grid.BCtypeS[Local_SolnBlk[nb].Grid.INl-1] == BC_WALL_VISCOUS_ISOTHERMAL ||
	 Local_SolnBlk[nb].Grid.BCtypeS[Local_SolnBlk[nb].Grid.INl-1] == BC_MOVING_WALL_HEATFLUX ||
	 Local_SolnBlk[nb].Grid.BCtypeS[Local_SolnBlk[nb].Grid.INl-1] == BC_MOVING_WALL_ISOTHERMAL ||
	 Local_SolnBlk[nb].Grid.BCtypeS[Local_SolnBlk[nb].Grid.INl-1] == BC_BURNING_SURFACE) ||
	(Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_REFLECTION ||
	 Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX ||
	 Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
	 Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_MOVING_WALL_HEATFLUX ||
	 Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL ||
	 Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_BURNING_SURFACE)) {

      if (Local_SolnBlk[nb].Grid.BCtypeW[j] != BC_REFLECTION &&
	  Local_SolnBlk[nb].Grid.BCtypeW[j] != BC_WALL_VISCOUS_HEATFLUX &&
	  Local_SolnBlk[nb].Grid.BCtypeW[j] != BC_WALL_VISCOUS_ISOTHERMAL &&
	  Local_SolnBlk[nb].Grid.BCtypeW[j] != BC_MOVING_WALL_HEATFLUX &&
	  Local_SolnBlk[nb].Grid.BCtypeW[j] != BC_MOVING_WALL_ISOTHERMAL &&
	  Local_SolnBlk[nb].Grid.BCtypeW[j] != BC_BURNING_SURFACE) {
	for (int i = Local_SolnBlk[nb].ICl-2; i <= Local_SolnBlk[nb].ICl-1; i++) {
	  Mesh[nb].cell_status[i][Local_SolnBlk[nb].JCl-1] = Mesh[nb].cell_status[i][Local_SolnBlk[nb].JCl  ];
	  Mesh[nb].cell_status[i][Local_SolnBlk[nb].JCl-2] = Mesh[nb].cell_status[i][Local_SolnBlk[nb].JCl+1];
	  Mesh[nb].cell_type[i][Local_SolnBlk[nb].JCl-1] = Mesh[nb].cell_type[i][Local_SolnBlk[nb].JCl  ];
	  Mesh[nb].cell_type[i][Local_SolnBlk[nb].JCl-2] = Mesh[nb].cell_type[i][Local_SolnBlk[nb].JCl+1];
	}

      } else if (Local_SolnBlk[nb].Grid.BCtypeS[Local_SolnBlk[nb].Grid.INl-1] != BC_REFLECTION ||
		 Local_SolnBlk[nb].Grid.BCtypeS[Local_SolnBlk[nb].Grid.INl-1] != BC_WALL_VISCOUS_HEATFLUX ||
		 Local_SolnBlk[nb].Grid.BCtypeS[Local_SolnBlk[nb].Grid.INl-1] != BC_WALL_VISCOUS_ISOTHERMAL ||
		 Local_SolnBlk[nb].Grid.BCtypeS[Local_SolnBlk[nb].Grid.INl-1] != BC_MOVING_WALL_HEATFLUX ||
		 Local_SolnBlk[nb].Grid.BCtypeS[Local_SolnBlk[nb].Grid.INl-1] != BC_MOVING_WALL_ISOTHERMAL ||
		 Local_SolnBlk[nb].Grid.BCtypeS[Local_SolnBlk[nb].Grid.INl-1] != BC_BURNING_SURFACE) {
	Mesh[nb].cell_status[Local_SolnBlk[nb].ICl-1][j] = Mesh[nb].cell_status[Local_SolnBlk[nb].ICl  ][j];
	Mesh[nb].cell_status[Local_SolnBlk[nb].ICl-2][j] = Mesh[nb].cell_status[Local_SolnBlk[nb].ICl+1][j];
	Mesh[nb].cell_type[Local_SolnBlk[nb].ICl-1][j] = Mesh[nb].cell_type[Local_SolnBlk[nb].ICl  ][j];
	Mesh[nb].cell_type[Local_SolnBlk[nb].ICl-2][j] = Mesh[nb].cell_type[Local_SolnBlk[nb].ICl+1][j];

      } else if ((Local_SolnBlk[nb].Grid.BCtypeS[Local_SolnBlk[nb].Grid.INl-1] == BC_REFLECTION ||
		  Local_SolnBlk[nb].Grid.BCtypeS[Local_SolnBlk[nb].Grid.INl-1] == BC_WALL_VISCOUS_HEATFLUX ||
		  Local_SolnBlk[nb].Grid.BCtypeS[Local_SolnBlk[nb].Grid.INl-1] == BC_WALL_VISCOUS_ISOTHERMAL ||
		  Local_SolnBlk[nb].Grid.BCtypeS[Local_SolnBlk[nb].Grid.INl-1] == BC_MOVING_WALL_HEATFLUX ||
		  Local_SolnBlk[nb].Grid.BCtypeS[Local_SolnBlk[nb].Grid.INl-1] == BC_MOVING_WALL_ISOTHERMAL ||
		  Local_SolnBlk[nb].Grid.BCtypeS[Local_SolnBlk[nb].Grid.INl-1] == BC_BURNING_SURFACE) &&
		 (Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_REFLECTION ||
		  Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX ||
		  Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
		  Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_MOVING_WALL_HEATFLUX ||
		  Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL ||
		  Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_BURNING_SURFACE)) {
	Mesh[nb].cell_status[Local_SolnBlk[nb].ICl-1][j] = Mesh[nb].cell_status[Local_SolnBlk[nb].ICl  ][j];
	Mesh[nb].cell_status[Local_SolnBlk[nb].ICl-2][j] = Mesh[nb].cell_status[Local_SolnBlk[nb].ICl+1][j];
	Mesh[nb].cell_type[Local_SolnBlk[nb].ICl-1][j] = Mesh[nb].cell_type[Local_SolnBlk[nb].ICl  ][j];
	Mesh[nb].cell_type[Local_SolnBlk[nb].ICl-2][j] = Mesh[nb].cell_type[Local_SolnBlk[nb].ICl+1][j];

      }
    }

    if ((Local_SolnBlk[nb].Grid.BCtypeS[Local_SolnBlk[nb].Grid.INl-1] == BC_REFLECTION ||
	 Local_SolnBlk[nb].Grid.BCtypeS[Local_SolnBlk[nb].Grid.INl-1] == BC_WALL_VISCOUS_HEATFLUX ||
	 Local_SolnBlk[nb].Grid.BCtypeS[Local_SolnBlk[nb].Grid.INl-1] == BC_WALL_VISCOUS_ISOTHERMAL ||
	 Local_SolnBlk[nb].Grid.BCtypeS[Local_SolnBlk[nb].Grid.INl-1] == BC_MOVING_WALL_HEATFLUX ||
	 Local_SolnBlk[nb].Grid.BCtypeS[Local_SolnBlk[nb].Grid.INl-1] == BC_MOVING_WALL_ISOTHERMAL ||
	 Local_SolnBlk[nb].Grid.BCtypeS[Local_SolnBlk[nb].Grid.INl-1] == BC_BURNING_SURFACE) ||
	(Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_REFLECTION ||
	 Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX ||
	 Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
	 Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_MOVING_WALL_HEATFLUX ||
	 Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL ||
	 Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_BURNING_SURFACE)) {

      if (Local_SolnBlk[nb].Grid.BCtypeE[j] != BC_REFLECTION &&
	  Local_SolnBlk[nb].Grid.BCtypeE[j] != BC_WALL_VISCOUS_HEATFLUX &&
	  Local_SolnBlk[nb].Grid.BCtypeE[j] != BC_WALL_VISCOUS_ISOTHERMAL &&
	  Local_SolnBlk[nb].Grid.BCtypeE[j] != BC_MOVING_WALL_HEATFLUX &&
	  Local_SolnBlk[nb].Grid.BCtypeE[j] != BC_MOVING_WALL_ISOTHERMAL &&
	  Local_SolnBlk[nb].Grid.BCtypeE[j] != BC_BURNING_SURFACE) {
	for (int i = Local_SolnBlk[nb].ICu+1; i <= Local_SolnBlk[nb].ICu+Local_SolnBlk[nb].Grid.Nghost; i++) {
	  Mesh[nb].cell_status[i][Local_SolnBlk[nb].JCl-1] = Mesh[nb].cell_status[i][Local_SolnBlk[nb].JCl  ];
	  Mesh[nb].cell_status[i][Local_SolnBlk[nb].JCl-2] = Mesh[nb].cell_status[i][Local_SolnBlk[nb].JCl+1];
	  Mesh[nb].cell_type[i][Local_SolnBlk[nb].JCl-1] = Mesh[nb].cell_type[i][Local_SolnBlk[nb].JCl  ];
	  Mesh[nb].cell_type[i][Local_SolnBlk[nb].JCl-2] = Mesh[nb].cell_type[i][Local_SolnBlk[nb].JCl+1];
	}

      } else if (Local_SolnBlk[nb].Grid.BCtypeS[Local_SolnBlk[nb].Grid.INl-1] != BC_REFLECTION ||
		 Local_SolnBlk[nb].Grid.BCtypeS[Local_SolnBlk[nb].Grid.INl-1] != BC_WALL_VISCOUS_HEATFLUX ||
		 Local_SolnBlk[nb].Grid.BCtypeS[Local_SolnBlk[nb].Grid.INl-1] != BC_WALL_VISCOUS_ISOTHERMAL ||
		 Local_SolnBlk[nb].Grid.BCtypeS[Local_SolnBlk[nb].Grid.INl-1] != BC_MOVING_WALL_HEATFLUX ||
		 Local_SolnBlk[nb].Grid.BCtypeS[Local_SolnBlk[nb].Grid.INl-1] != BC_MOVING_WALL_ISOTHERMAL ||
		 Local_SolnBlk[nb].Grid.BCtypeS[Local_SolnBlk[nb].Grid.INl-1] != BC_BURNING_SURFACE) {
	Mesh[nb].cell_status[Local_SolnBlk[nb].ICu+1][j] = Mesh[nb].cell_status[Local_SolnBlk[nb].ICu  ][j];
	Mesh[nb].cell_status[Local_SolnBlk[nb].ICu+2][j] = Mesh[nb].cell_status[Local_SolnBlk[nb].ICu-1][j];
	Mesh[nb].cell_type[Local_SolnBlk[nb].ICu+1][j] = Mesh[nb].cell_type[Local_SolnBlk[nb].ICu  ][j];
	Mesh[nb].cell_type[Local_SolnBlk[nb].ICu+2][j] = Mesh[nb].cell_type[Local_SolnBlk[nb].ICu-1][j];

      } else if ((Local_SolnBlk[nb].Grid.BCtypeS[Local_SolnBlk[nb].Grid.INl-1] == BC_REFLECTION ||
		  Local_SolnBlk[nb].Grid.BCtypeS[Local_SolnBlk[nb].Grid.INl-1] == BC_WALL_VISCOUS_HEATFLUX ||
		  Local_SolnBlk[nb].Grid.BCtypeS[Local_SolnBlk[nb].Grid.INl-1] == BC_WALL_VISCOUS_ISOTHERMAL ||
		  Local_SolnBlk[nb].Grid.BCtypeS[Local_SolnBlk[nb].Grid.INl-1] == BC_MOVING_WALL_HEATFLUX ||
		  Local_SolnBlk[nb].Grid.BCtypeS[Local_SolnBlk[nb].Grid.INl-1] == BC_MOVING_WALL_ISOTHERMAL ||
		  Local_SolnBlk[nb].Grid.BCtypeS[Local_SolnBlk[nb].Grid.INl-1] == BC_BURNING_SURFACE) &&
		 (Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_REFLECTION ||
		  Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX ||
		  Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
		  Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_MOVING_WALL_HEATFLUX ||
		  Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL ||
		  Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_BURNING_SURFACE)) {
	Mesh[nb].cell_status[Local_SolnBlk[nb].ICu+1][j] = Mesh[nb].cell_status[Local_SolnBlk[nb].ICu  ][j];
	Mesh[nb].cell_status[Local_SolnBlk[nb].ICu+2][j] = Mesh[nb].cell_status[Local_SolnBlk[nb].ICu-1][j];
	Mesh[nb].cell_type[Local_SolnBlk[nb].ICu+1][j] = Mesh[nb].cell_type[Local_SolnBlk[nb].ICu  ][j];
	Mesh[nb].cell_type[Local_SolnBlk[nb].ICu+2][j] = Mesh[nb].cell_type[Local_SolnBlk[nb].ICu-1][j];

      }

    }

  }

  // Update NW and NE corner nodes.
  for (int j = Local_SolnBlk[nb].JCu+1; j <= Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Grid.Nghost; j++) {

    if ((Local_SolnBlk[nb].Grid.BCtypeN[Local_SolnBlk[nb].Grid.INl-1] == BC_REFLECTION ||
	 Local_SolnBlk[nb].Grid.BCtypeN[Local_SolnBlk[nb].Grid.INl-1] == BC_WALL_VISCOUS_HEATFLUX ||
	 Local_SolnBlk[nb].Grid.BCtypeN[Local_SolnBlk[nb].Grid.INl-1] == BC_WALL_VISCOUS_ISOTHERMAL ||
	 Local_SolnBlk[nb].Grid.BCtypeN[Local_SolnBlk[nb].Grid.INl-1] == BC_MOVING_WALL_HEATFLUX ||
	 Local_SolnBlk[nb].Grid.BCtypeN[Local_SolnBlk[nb].Grid.INl-1] == BC_MOVING_WALL_ISOTHERMAL ||
	 Local_SolnBlk[nb].Grid.BCtypeN[Local_SolnBlk[nb].Grid.INl-1] == BC_BURNING_SURFACE) ||
	(Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_REFLECTION ||
	 Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX ||
	 Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
	 Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_MOVING_WALL_HEATFLUX ||
	 Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL ||
	 Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_BURNING_SURFACE)) {

      if (Local_SolnBlk[nb].Grid.BCtypeW[j] != BC_REFLECTION &&
	  Local_SolnBlk[nb].Grid.BCtypeW[j] != BC_WALL_VISCOUS_HEATFLUX &&
	  Local_SolnBlk[nb].Grid.BCtypeW[j] != BC_WALL_VISCOUS_ISOTHERMAL &&
	  Local_SolnBlk[nb].Grid.BCtypeW[j] != BC_MOVING_WALL_HEATFLUX &&
	  Local_SolnBlk[nb].Grid.BCtypeW[j] != BC_MOVING_WALL_ISOTHERMAL &&
	  Local_SolnBlk[nb].Grid.BCtypeW[j] != BC_BURNING_SURFACE) {
	for (int i = Local_SolnBlk[nb].ICl-2; i <= Local_SolnBlk[nb].ICl-1; i++) {
	  Mesh[nb].cell_status[i][Local_SolnBlk[nb].JCu+1] = Mesh[nb].cell_status[i][Local_SolnBlk[nb].JCu  ];
	  Mesh[nb].cell_status[i][Local_SolnBlk[nb].JCu+2] = Mesh[nb].cell_status[i][Local_SolnBlk[nb].JCu-1];
	  Mesh[nb].cell_type[i][Local_SolnBlk[nb].JCu+1] = Mesh[nb].cell_type[i][Local_SolnBlk[nb].JCu  ];
	  Mesh[nb].cell_type[i][Local_SolnBlk[nb].JCu+2] = Mesh[nb].cell_type[i][Local_SolnBlk[nb].JCu-1];
	}

      } else if (Local_SolnBlk[nb].Grid.BCtypeN[Local_SolnBlk[nb].Grid.INl-1] != BC_REFLECTION ||
		 Local_SolnBlk[nb].Grid.BCtypeN[Local_SolnBlk[nb].Grid.INl-1] != BC_WALL_VISCOUS_HEATFLUX ||
		 Local_SolnBlk[nb].Grid.BCtypeN[Local_SolnBlk[nb].Grid.INl-1] != BC_WALL_VISCOUS_ISOTHERMAL ||
		 Local_SolnBlk[nb].Grid.BCtypeN[Local_SolnBlk[nb].Grid.INl-1] != BC_MOVING_WALL_HEATFLUX ||
		 Local_SolnBlk[nb].Grid.BCtypeN[Local_SolnBlk[nb].Grid.INl-1] != BC_MOVING_WALL_ISOTHERMAL ||
		 Local_SolnBlk[nb].Grid.BCtypeN[Local_SolnBlk[nb].Grid.INl-1] != BC_BURNING_SURFACE) {
	Mesh[nb].cell_status[Local_SolnBlk[nb].ICl-1][j] = Mesh[nb].cell_status[Local_SolnBlk[nb].ICl  ][j];
	Mesh[nb].cell_status[Local_SolnBlk[nb].ICl-2][j] = Mesh[nb].cell_status[Local_SolnBlk[nb].ICl+1][j];
	Mesh[nb].cell_type[Local_SolnBlk[nb].ICl-1][j] = Mesh[nb].cell_type[Local_SolnBlk[nb].ICl  ][j];
	Mesh[nb].cell_type[Local_SolnBlk[nb].ICl-2][j] = Mesh[nb].cell_type[Local_SolnBlk[nb].ICl+1][j];

      } else if ((Local_SolnBlk[nb].Grid.BCtypeN[Local_SolnBlk[nb].Grid.INl-1] == BC_REFLECTION ||
		  Local_SolnBlk[nb].Grid.BCtypeN[Local_SolnBlk[nb].Grid.INl-1] == BC_WALL_VISCOUS_HEATFLUX ||
		  Local_SolnBlk[nb].Grid.BCtypeN[Local_SolnBlk[nb].Grid.INl-1] == BC_WALL_VISCOUS_ISOTHERMAL ||
		  Local_SolnBlk[nb].Grid.BCtypeN[Local_SolnBlk[nb].Grid.INl-1] == BC_MOVING_WALL_HEATFLUX ||
		  Local_SolnBlk[nb].Grid.BCtypeN[Local_SolnBlk[nb].Grid.INl-1] == BC_MOVING_WALL_ISOTHERMAL ||
		  Local_SolnBlk[nb].Grid.BCtypeN[Local_SolnBlk[nb].Grid.INl-1] == BC_BURNING_SURFACE) &&
		 (Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_REFLECTION ||
		  Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX ||
		  Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
		  Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_MOVING_WALL_HEATFLUX ||
		  Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL ||
		  Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_BURNING_SURFACE)) {
	Mesh[nb].cell_status[Local_SolnBlk[nb].ICl-1][j] = Mesh[nb].cell_status[Local_SolnBlk[nb].ICl  ][j];
	Mesh[nb].cell_status[Local_SolnBlk[nb].ICl-2][j] = Mesh[nb].cell_status[Local_SolnBlk[nb].ICl+1][j];
	Mesh[nb].cell_type[Local_SolnBlk[nb].ICl-1][j] = Mesh[nb].cell_type[Local_SolnBlk[nb].ICl  ][j];
	Mesh[nb].cell_type[Local_SolnBlk[nb].ICl-2][j] = Mesh[nb].cell_type[Local_SolnBlk[nb].ICl+1][j];

      }
    }

    if ((Local_SolnBlk[nb].Grid.BCtypeN[Local_SolnBlk[nb].Grid.INu+1] == BC_REFLECTION ||
	 Local_SolnBlk[nb].Grid.BCtypeN[Local_SolnBlk[nb].Grid.INu+1] == BC_WALL_VISCOUS_HEATFLUX ||
	 Local_SolnBlk[nb].Grid.BCtypeN[Local_SolnBlk[nb].Grid.INu+1] == BC_WALL_VISCOUS_ISOTHERMAL ||
	 Local_SolnBlk[nb].Grid.BCtypeN[Local_SolnBlk[nb].Grid.INu+1] == BC_MOVING_WALL_HEATFLUX ||
	 Local_SolnBlk[nb].Grid.BCtypeN[Local_SolnBlk[nb].Grid.INu+1] == BC_MOVING_WALL_ISOTHERMAL ||
	 Local_SolnBlk[nb].Grid.BCtypeN[Local_SolnBlk[nb].Grid.INu+1] == BC_BURNING_SURFACE) ||
	(Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_REFLECTION ||
	 Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX ||
	 Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
	 Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_MOVING_WALL_HEATFLUX ||
	 Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL ||
	 Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_BURNING_SURFACE)) {

      if (Local_SolnBlk[nb].Grid.BCtypeE[j] != BC_REFLECTION &&
	  Local_SolnBlk[nb].Grid.BCtypeE[j] != BC_WALL_VISCOUS_HEATFLUX &&
	  Local_SolnBlk[nb].Grid.BCtypeE[j] != BC_WALL_VISCOUS_ISOTHERMAL &&
	  Local_SolnBlk[nb].Grid.BCtypeE[j] != BC_MOVING_WALL_HEATFLUX &&
	  Local_SolnBlk[nb].Grid.BCtypeE[j] != BC_MOVING_WALL_ISOTHERMAL &&
	  Local_SolnBlk[nb].Grid.BCtypeE[j] != BC_BURNING_SURFACE) {
	for (int i = Local_SolnBlk[nb].ICu+1; i <= Local_SolnBlk[nb].ICu+2; i++) {
	  Mesh[nb].cell_status[i][Local_SolnBlk[nb].JCu+1] = Mesh[nb].cell_status[i][Local_SolnBlk[nb].JCu  ];
	  Mesh[nb].cell_status[i][Local_SolnBlk[nb].JCu+2] = Mesh[nb].cell_status[i][Local_SolnBlk[nb].JCu-1];
	  Mesh[nb].cell_type[i][Local_SolnBlk[nb].JCu+1] = Mesh[nb].cell_type[i][Local_SolnBlk[nb].JCu  ];
	  Mesh[nb].cell_type[i][Local_SolnBlk[nb].JCu+2] = Mesh[nb].cell_type[i][Local_SolnBlk[nb].JCu-1];
	}

      } else if (Local_SolnBlk[nb].Grid.BCtypeN[Local_SolnBlk[nb].Grid.INu+1] != BC_REFLECTION ||
		 Local_SolnBlk[nb].Grid.BCtypeN[Local_SolnBlk[nb].Grid.INu+1] != BC_WALL_VISCOUS_HEATFLUX ||
		 Local_SolnBlk[nb].Grid.BCtypeN[Local_SolnBlk[nb].Grid.INu+1] != BC_WALL_VISCOUS_ISOTHERMAL ||
		 Local_SolnBlk[nb].Grid.BCtypeN[Local_SolnBlk[nb].Grid.INu+1] != BC_MOVING_WALL_HEATFLUX ||
		 Local_SolnBlk[nb].Grid.BCtypeN[Local_SolnBlk[nb].Grid.INu+1] != BC_MOVING_WALL_ISOTHERMAL ||
		 Local_SolnBlk[nb].Grid.BCtypeN[Local_SolnBlk[nb].Grid.INu+1] != BC_BURNING_SURFACE) {
	Mesh[nb].cell_status[Local_SolnBlk[nb].ICu+1][j] = Mesh[nb].cell_status[Local_SolnBlk[nb].ICu  ][j];
	Mesh[nb].cell_status[Local_SolnBlk[nb].ICu+2][j] = Mesh[nb].cell_status[Local_SolnBlk[nb].ICu-1][j];
	Mesh[nb].cell_type[Local_SolnBlk[nb].ICu+1][j] = Mesh[nb].cell_type[Local_SolnBlk[nb].ICu  ][j];
	Mesh[nb].cell_type[Local_SolnBlk[nb].ICu+2][j] = Mesh[nb].cell_type[Local_SolnBlk[nb].ICu-1][j];

      } else if ((Local_SolnBlk[nb].Grid.BCtypeN[Local_SolnBlk[nb].Grid.INu+1] == BC_REFLECTION ||
		  Local_SolnBlk[nb].Grid.BCtypeN[Local_SolnBlk[nb].Grid.INu+1] == BC_WALL_VISCOUS_HEATFLUX ||
		  Local_SolnBlk[nb].Grid.BCtypeN[Local_SolnBlk[nb].Grid.INu+1] == BC_WALL_VISCOUS_ISOTHERMAL ||
		  Local_SolnBlk[nb].Grid.BCtypeN[Local_SolnBlk[nb].Grid.INu+1] == BC_MOVING_WALL_HEATFLUX ||
		  Local_SolnBlk[nb].Grid.BCtypeN[Local_SolnBlk[nb].Grid.INu+1] == BC_MOVING_WALL_ISOTHERMAL ||
		  Local_SolnBlk[nb].Grid.BCtypeN[Local_SolnBlk[nb].Grid.INu+1] == BC_BURNING_SURFACE) &&
		 (Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_REFLECTION ||
		  Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX ||
		  Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
		  Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_MOVING_WALL_HEATFLUX ||
		  Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL ||
		  Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_BURNING_SURFACE)) {
	Mesh[nb].cell_status[Local_SolnBlk[nb].ICu+1][j] = Mesh[nb].cell_status[Local_SolnBlk[nb].ICu  ][j];
	Mesh[nb].cell_status[Local_SolnBlk[nb].ICu+2][j] = Mesh[nb].cell_status[Local_SolnBlk[nb].ICu-1][j];
	Mesh[nb].cell_type[Local_SolnBlk[nb].ICu+1][j] = Mesh[nb].cell_type[Local_SolnBlk[nb].ICu  ][j];
	Mesh[nb].cell_type[Local_SolnBlk[nb].ICu+2][j] = Mesh[nb].cell_type[Local_SolnBlk[nb].ICu-1][j];

      }

    }

  }

}

/**********************************************************************
 **********************************************************************
 ** Moving embedded boundary routines.                               **
 **********************************************************************
 **********************************************************************/

/**********************************************************************
 * EmbeddedBoundaries2D::Compute_Interface_Location --                *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
inline int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Compute_Interface_Location(const int &batch_flag,
			   const double &current_time,
			   const double &maximum_time,
			   int &evolution_counter,
			   int &number_of_level_set_time_steps,
			   double &level_set_current_time) {

  // Exit immediately if no interface components have been specified.
  if (!Interface_Component_List.Ni) return 0;

  // Go to end of function if interface evolution is not required.
  if (evolution_counter != IP->Interface_IP.Evolution_Frequency) {
    evolution_counter++;
    goto compute_velocity_function;
  } else {
    evolution_counter = 1;
  }

  int error_flag, motion_flag, levelset_flag;
  double dTime, radius;

  // Evaluate the time-step.
  dTime = maximum_time-current_time;

  // Initialize flags.
  motion_flag = OFF;
  levelset_flag = OFF;

  // Determine the velocity function for each embedded boundary.
  for (int ni = 1; ni <= Interface_Component_List.Ni; ni++) {

    // Recompute interface parameters depending on motion type.
    if (Interface_Component_List[ni].Motion == INTERFACE_MOTION_STATIONARY) {
      // The interface is stationary, do nothing.

    } else if (Interface_Component_List[ni].Motion == INTERFACE_MOTION_CONSTANT ||
	       Interface_Component_List[ni].Motion == INTERFACE_MOTION_EXPAND) {
      // The interface has a constant speed in the radial direction.
      motion_flag = ON;
      switch(Interface_Component_List[ni].Type) {
      case INTERFACE_CIRCLE :
	radius = Interface_Component_List[ni].Length1 + current_time*Interface_Component_List[ni].Speed.x;
	Interface_Component_List[ni].Scale(ONE+dTime*Interface_Component_List[ni].Speed.x/radius);
	break;
      default :
	return 110101;
	break;
      };

    } else if (Interface_Component_List[ni].Motion == INTERFACE_MOTION_UNIFORM ||
	       Interface_Component_List[ni].Motion == INTERFACE_MOTION_TRANSLATE) {
      // The interface has a constant velocity.
      motion_flag = ON;
      switch(Interface_Component_List[ni].Type) {
      case INTERFACE_LINE :
      case INTERFACE_CIRCLE :
      case INTERFACE_ELLIPSE :
      case INTERFACE_SQUARE :
      case INTERFACE_RECTANGLE :
      case INTERFACE_ZALESAK :
      case INTERFACE_NACA0012_AEROFOIL :
      case INTERFACE_NACA0015_AEROFOIL :
      case INTERFACE_USER_SPECIFIED :
	Interface_Component_List[ni].Translate(dTime*Interface_Component_List[ni].Speed);
	break;
      default :
	return 110102;
	break;
      };

    } else if (Interface_Component_List[ni].Motion == INTERFACE_MOTION_ROTATE) {
      // The interface has a constant rotational velocity.
      motion_flag = ON;
      switch(Interface_Component_List[ni].Type) {
      case INTERFACE_CIRCLE :
      case INTERFACE_ELLIPSE :
	Interface_Component_List[ni].Translate(Vector2D(-Interface_Component_List[ni].Speed.x*sin(TWO*PI*current_time/Interface_Component_List[ni].Speed.y),ZERO));
	Interface_Component_List[ni].Translate(Vector2D(Interface_Component_List[ni].Speed.x*sin(TWO*PI*maximum_time/Interface_Component_List[ni].Speed.y),ZERO));
	break;
      case INTERFACE_SQUARE :
      case INTERFACE_RECTANGLE :
      case INTERFACE_USER_SPECIFIED :
      case INTERFACE_ZALESAK :
	Interface_Component_List[ni].Rotate(dTime*Interface_Component_List[ni].Speed.x*PI/180.0);
	break;
      case INTERFACE_NACA0012_AEROFOIL :
	Interface_Component_List[ni].Rotate( (PI/180.0)*(0.016 + 2.51*sin(2*PI*62.5*(current_time))));
	Interface_Component_List[ni].Rotate(-(PI/180.0)*(0.016 + 2.51*sin(2*PI*62.5*(maximum_time))));
	break;
      case INTERFACE_NACA0015_AEROFOIL :
	Interface_Component_List[ni].Rotate( (PI/180.0)*(0.016 + 2.51*sin(2*PI*62.5*(current_time))));
	Interface_Component_List[ni].Rotate(-(PI/180.0)*(0.016 + 2.51*sin(2*PI*62.5*(maximum_time))));
	break;
	//case INTERFACE_ZALESAK :
	//Interface_Component_List[ni].Rotate(15.0*PI/180.0);
	//dTime*Interface_Component_List[ni].Speed.x*PI/180.0);
	//break;
      default :
	return 110103;
	break;
      };

    } else if (Interface_Component_List[ni].Motion == INTERFACE_MOTION_MOMENTUM_TRANSFER) {
      motion_flag = ON;
      return 110104;

    } else if (Interface_Component_List[ni].Motion == INTERFACE_MOTION_LEVELSET_EXPAND ||
	       Interface_Component_List[ni].Motion == INTERFACE_MOTION_LEVELSET_STRETCH ||
	       Interface_Component_List[ni].Motion == INTERFACE_MOTION_LEVELSET) {
      motion_flag = ON;
      levelset_flag = ON;

    } else if (Interface_Component_List[ni].Motion == INTERFACE_MOTION_BURNING_SURFACE) {
      motion_flag = ON;
      levelset_flag = ON;

    } else {
      // Do nothing.
      return 110105;
    }

  }

  // If none of the interface components have underwent any form of 
  // motion then exit immediately.
  if (!motion_flag) return 0;

  // Evolve the level set solution if required.
  if (levelset_flag) {

    // Evolve level set solution.
    error_flag = Evolve_Level_Set_Solution(batch_flag,
					   LS_Local_SolnBlk,
					   *LS_IP,
					   *LS_QuadTree,
					   *LS_Global_Solution_Block_List,
					   *LS_Local_Solution_Block_List,
					   level_set_current_time,
					   maximum_time,
					   number_of_level_set_time_steps);
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;

    // Reconstruct the interface component list.
    error_flag = Reconstruct_Interface_Component_List();
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;

  }

  // Construct interface union list.
  error_flag = Construct_Interface_Union_List();
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return error_flag;

  // Store adjusted mesh.
  Store_Adjusted_Mesh();

  // Unadjust the mesh.
  Mesh_Unadjustment();

  //added by ~james
  error_flag = Send_All_Messages(Local_SolnBlk,
				 *Local_Solution_Block_List,
 				 NUM_COMP_VECTOR2D,
 				 ON);

  // Conduct mesh adjustment algorithm.
  error_flag = Mesh_Adjustment(OFF,OFF);
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return error_flag;

  // Redistribute solution content and residual to/from activated/dectivated cells.
  error_flag = Redistribute_Solution_Content(maximum_time);
  if (error_flag) { cout << endl << error_flag; cout.flush(); }
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) {
    error_flag = Output_Tecplot(Local_SolnBlk,
				*Local_Solution_Block_List,
				*IP,0,ZERO);
    error_flag = Output_Elements_Tecplot(0,ZERO);
    error_flag = Output_Nodes_Tecplot();
    error_flag = Output_Cell_Status_Tecplot();
    error_flag = Output_Interface_Component_List_Tecplot();
    error_flag = Output_Interface_Union_List_Tecplot();
    if (levelset_flag) {
      error_flag = Output_Level_Set_Tecplot(0,ZERO);
      error_flag = Output_Level_Set_Cells_Tecplot(0,ZERO);
      error_flag = Output_Level_Set_Interface_Tecplot(0,ZERO);
    }
    if (CFFC_Primary_MPI_Processor()) {
      cout << endl << " time = " << current_time; cout.flush();
      cout << endl << " time = " << maximum_time; cout.flush();
    }
    return 87878;//error_flag;
  }

  // Store adjusted mesh.
  Store_Adjusted_Mesh();

  compute_velocity_function: ;

  // Compute the interface velocity function.
  error_flag = Compute_Interface_Velocity_Function(maximum_time);
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return error_flag;

  // Apply interface boundary conditions.
  error_flag = Boundary_Conditions(maximum_time);
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return error_flag;

  // The new interface location has been successfully computed.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Compute_Interface_Velocity_Function --       *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
inline int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Compute_Interface_Velocity_Function(const double &Time) {

  // Exit immediately if no interface components have been specified.
  if (!Interface_Component_List.Ni) return 0;

  int error_flag, share_flag, levelset_flag, Ni;
  Interface2D_List Interface_List;

  // Initialize flags.
  share_flag = OFF;
  levelset_flag = OFF;

  // Reset the velocity function for each interface.
  for (int ni = 1; ni <= Interface_Component_List.Ni; ni++) {
    Interface_Component_List[ni].Zero_Velocity_Function();
  }

  // Determine the velocity function on each of the interface components.
  for (int ni = 1; ni <= Interface_Component_List.Ni; ni++) {
    // Compute interface velocity function on each interface.
    switch(Interface_Component_List[ni].Motion) {
    case INTERFACE_MOTION_STATIONARY :
      // The interface is stationary, do nothing.
      break;
    case INTERFACE_MOTION_CONSTANT :
    case INTERFACE_MOTION_EXPAND :
    case INTERFACE_MOTION_UNIFORM :
    case INTERFACE_MOTION_TRANSLATE :
      Interface_Component_List[ni].Set_Velocity_Function(ZERO);
      break;
    case INTERFACE_MOTION_STRETCH :
      Interface_Component_List[ni].Set_Velocity_Function(ZERO);
      break;
    case INTERFACE_MOTION_ROTATE :
      Interface_Component_List[ni].Set_Velocity_Function(Time);
      break;
    case INTERFACE_MOTION_BURNING_SURFACE :
      // The interface is a burning propellant surface.
      for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
	if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {
	  error_flag = Compute_Interface_Velocity_Function(nb,
							   ni,
							   INTERFACE_MOTION_BURNING_SURFACE);
	  if (error_flag) return error_flag;
	}
      }
      share_flag = ON;
      levelset_flag = ON;
      break;
    case INTERFACE_MOTION_MOMENTUM_TRANSFER :
      share_flag = ON;
      return 2502;
      break;
    case INTERFACE_MOTION_LEVELSET_STATIONARY :
      // The interface is stationary, do nothing.
      break;
    case INTERFACE_MOTION_LEVELSET_EXPAND :
      break;
    case INTERFACE_MOTION_LEVELSET_STRETCH :
      break;
    case INTERFACE_MOTION_LEVELSET_BULKFLOW :
      share_flag = ON;
      levelset_flag = ON;
      return 2503;
      break;
    case INTERFACE_MOTION_LEVELSET :
      levelset_flag = OFF;
      break;
    default :
      return 2501;
      break;
    };
  }

  // Synchronize flags on all processors.
  share_flag = CFFC_Maximum_MPI(share_flag);
  levelset_flag = CFFC_Maximum_MPI(levelset_flag);

  // If the share flag has been turned on then one or more of the
  // embedded boundaries involves a motion computed on its local block
  // that must be known across all blocks on all processors.
  if (share_flag) {

#ifdef _MPI_VERSION

    Interface2D_List temp_Interface_List;
    MPI::Intracomm local_comm;
    MPI::Group world_group = MPI::COMM_WORLD.Get_group();
    MPI::Group local_group;
    int *world_list;
    int *local_list;
    world_list = new int[QuadTree->Ncpu];
    local_list = new int[QuadTree->Ncpu];
    for (int ncpu = 0; ncpu < QuadTree->Ncpu; ncpu++) {
      world_list[ncpu] = ncpu;
      local_list[ncpu] = ncpu;
    }

    // Create a temporary interface list.
    Interface_List.Copy(Interface_Component_List);

    // Create a second temporary interface list.
    temp_Interface_List.Copy(Interface_Component_List);

    // Broadcast all of the interface information owned by each individual 
    // processor to the rest of the processors.
    for (int iCPU = 0; iCPU < QuadTree->Ncpu; iCPU++) {
      // Create the local list and communicator.
      for (int ncpu = 0; ncpu < QuadTree->Ncpu; ncpu++) {
 	if (ncpu == 0) {
 	  local_list[ncpu] = iCPU;
 	} else if (ncpu == iCPU) {
 	  local_list[ncpu] = 0;
 	} else {
 	  local_list[ncpu] = ncpu;
 	}
      }
      local_group = world_group.Incl(QuadTree->Ncpu,local_list);
      local_comm  = MPI::COMM_WORLD.Create(local_group);
      // Load message buffer.
      if (iCPU == Global_Solution_Block_List->ThisCPU) {
	temp_Interface_List.Copy(Interface_List);
      }
      // Broadcast the buffer interface velocity function information.
      Broadcast_Interface_List(temp_Interface_List,
			       local_comm,
			       iCPU);
      // Unload message buffer.
      for (int ni = 1; ni <= Interface_List.Ni; ni++) {
	if (Interface_List[ni].Motion == INTERFACE_MOTION_BURNING_SURFACE ||
	    Interface_List[ni].Motion == INTERFACE_MOTION_MOMENTUM_TRANSFER) {
	  for (int np = 0; np < Interface_List[ni].Spline.np; np++) {
	    if (fabs(temp_Interface_List[ni].F[np].abs()) > TOLER) 
	      Interface_List[ni].F[np] = temp_Interface_List[ni].F[np];
	  }
	}
      }
      // Free local communicator and group.
      if (local_comm != MPI::COMM_NULL) local_comm.Free();
      local_group.Free();
    }
    // Copy the information from the temporary list into the interface component list.
    for (int ni = 1; ni <= Interface_Component_List.Ni; ni++) {
      if (Interface_Component_List[ni].Motion == INTERFACE_MOTION_BURNING_SURFACE ||
	  Interface_Component_List[ni].Motion == INTERFACE_MOTION_MOMENTUM_TRANSFER) {
	Interface_Component_List[ni].Copy(Interface_List[ni]);
      }
    }

    // Deallocate the global interface list.
    temp_Interface_List.deallocate();

    // Deallocate the global interface list.
    Interface_List.deallocate();

    delete []world_list; world_list = NULL;
    delete []local_list; local_list = NULL;
#endif

  }

  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

  // If the level set flag has been turned on then one or more of the 
  // embedded boundaries are evolved using the level set method.  Update
  // the velocity function in the level set interface list.
  if (levelset_flag) {

    // Create a 'universal' interface component list that contains
    // interface components that are to be evolved using the level set 
    // method and for which the scalar extension problem must be 
    // computed.  Form this list on all processors, including processors
    // that do not have any active solution blocks since they may have
    // active level set 2D solution blocks.
    if (CFFC_Primary_MPI_Processor()) {
      if (Local_Solution_Block_List->Block[0].used == ADAPTIVEBLOCK2D_USED) {
	// Determine the number of interfaces evolved using the level
	// set method.
	Ni = 0;
	for (int ni = 1; ni <= Interface_Component_List.Ni; ni++) {
	  if (Interface_Component_List[ni].Motion == INTERFACE_MOTION_BURNING_SURFACE) {
	    Ni++;
	  }
	}
	Interface_List.allocate(Ni);
	Ni = 1;
	for (int ni = 1; ni <= Interface_Component_List.Ni; ni++) {
	  if (Interface_Component_List[ni].Motion == INTERFACE_MOTION_BURNING_SURFACE) {
	    Interface_List[Ni].Copy(Interface_Component_List[ni]);
	    //Interface_List[ni].Set_Normal_Velocity_Function();
	    Ni++;
	  }
	}
      }
    }
    Broadcast_Interface_List(Interface_List);

    // Copy the universal interface list to the level set interface list.
    error_flag = Set_Interface_List(LS_Local_SolnBlk,
				    *LS_Local_Solution_Block_List,
				    Interface_List);
    if (error_flag) return error_flag;

    // Solve scalar (front speed) geometric extension problem.
    error_flag = Scalar_Geometric_Extension_Problem(LS_Local_SolnBlk,
  						    *LS_Local_Solution_Block_List,
  						    *LS_IP);
    if (error_flag) return error_flag;

    // Solve scalar (front speed) extension equation.
    error_flag = Explicit_Scalar_Extension_Equation(LS_Local_SolnBlk,
   						    *LS_IP,
   						    *LS_QuadTree,
   						    *LS_Global_Solution_Block_List,
   						    *LS_Local_Solution_Block_List);
    if (error_flag) return error_flag;

  }

  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

  // Construct interface union list.
  error_flag = Construct_Interface_Union_List();
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return error_flag;

  // The new interface speed function has been successfully computed.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Compute_Interface_Velocity_Function --       *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
inline int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Compute_Interface_Velocity_Function(const int &nb,
				    const int &Ni,
				    const int &Motion_Flag) {

  // Exit immediately if the interface is not present in the current block.
  if (!Adjustment_Data[nb].Interface_Present[Ni]) return 0;

  // Exit immediately if the interface motion is not governed by
  // a physical process.  Only burning surfaces are currently
  // possible.
  if (Motion_Flag != INTERFACE_MOTION_BURNING_SURFACE) return 2901;

  Vector2D dX, dX2, norm_dir;
  Polygon P;
  pState Wb;

  // For each inactive cell that has active neighbours, save the 
  // position and normal speed at each face.
  for (int j = Local_SolnBlk[nb].JCl; j <= Local_SolnBlk[nb].JCu; j++) {
    for (int i = Local_SolnBlk[nb].ICl; i <= Local_SolnBlk[nb].ICu; i++) {

      if (Mesh[nb].cell_status[i][j] == Ni &&
 	  (Mesh[nb].cell_status[i-1][j] == CELL_STATUS_ACTIVE ||
 	   Mesh[nb].cell_status[i+1][j] == CELL_STATUS_ACTIVE ||
 	   Mesh[nb].cell_status[i][j-1] == CELL_STATUS_ACTIVE ||
 	   Mesh[nb].cell_status[i][j+1] == CELL_STATUS_ACTIVE)) {

	dX.x = 0.25*abs(OGrid[nb].xfaceE(i,j)-OGrid[nb].xfaceW(i,j));
	dX.y = 0.25*abs(OGrid[nb].xfaceN(i,j)-OGrid[nb].xfaceS(i,j));

	for (int np = 0; np < Interface_Component_List[Ni].Spline.np; np++) {

	  norm_dir = Interface_Component_List[Ni].normal(np);

	  if (Mesh[nb].cell_status[i][j+1] == CELL_STATUS_ACTIVE) {
	    P.convert(Local_SolnBlk[nb].Grid.nodeSW(i,j+1).X+Vector2D(-dX.x,-dX.y),
		      Local_SolnBlk[nb].Grid.nodeSE(i,j+1).X+Vector2D(dX.x,-dX.y),
		      Local_SolnBlk[nb].Grid.nodeNE(i,j+1).X,
		      Local_SolnBlk[nb].Grid.nodeNW(i,j+1).X);
	    if (P.point_in_polygon(Interface_Component_List[Ni].Spline.Xp[np])) {
	      dX2 = Interface_Component_List[Ni].Spline.Xp[np] - Local_SolnBlk[nb].Grid.Cell[i][j+1].Xc;
	      Wb = Local_SolnBlk[nb].W[i][j+1];
// 	      Wb = Local_SolnBlk[nb].W[i][j+1] + (Local_SolnBlk[nb].phi[i][j+1]^Local_SolnBlk[nb].dWdx[i][j+1])*dX2.x +
// 	                                         (Local_SolnBlk[nb].phi[i][j+1]^Local_SolnBlk[nb].dWdy[i][j+1])*dX2.y;
 	      Interface_Component_List[Ni].F[np] = Wb.burningrate()*norm_dir;//Local_SolnBlk[nb].Grid.nfaceN(i,j);
 	      //Interface_Component_List[Ni].F[np] = Local_SolnBlk[nb].U[i][j+1].burningrate()*Local_SolnBlk[nb].Grid.nfaceN(i,j);
	    }
	    P.deallocate();
	  }
	  if (Mesh[nb].cell_status[i][j-1] == CELL_STATUS_ACTIVE) {
	    P.convert(Local_SolnBlk[nb].Grid.nodeSW(i,j-1).X,Local_SolnBlk[nb].Grid.nodeSE(i,j-1).X,
		      Local_SolnBlk[nb].Grid.nodeNE(i,j-1).X+Vector2D(dX.x,dX.y),
		      Local_SolnBlk[nb].Grid.nodeNW(i,j-1).X+Vector2D(-dX.x,dX.y));
	    if (P.point_in_polygon(Interface_Component_List[Ni].Spline.Xp[np])) {
	      dX2 = Interface_Component_List[Ni].Spline.Xp[np] - Local_SolnBlk[nb].Grid.Cell[i][j-1].Xc;
	      Wb = Local_SolnBlk[nb].W[i][j-1];
// 	      Wb = Local_SolnBlk[nb].W[i][j-1] + (Local_SolnBlk[nb].phi[i][j-1]^Local_SolnBlk[nb].dWdx[i][j-1])*dX2.x +
// 	                                         (Local_SolnBlk[nb].phi[i][j-1]^Local_SolnBlk[nb].dWdy[i][j-1])*dX2.y;
 	      Interface_Component_List[Ni].F[np] = Wb.burningrate()*norm_dir;//Local_SolnBlk[nb].Grid.nfaceS(i,j);
 	      //Interface_Component_List[Ni].F[np] = Local_SolnBlk[nb].U[i][j-1].burningrate()*Local_SolnBlk[nb].Grid.nfaceS(i,j);
	    }
	    P.deallocate();
	  }
	  if (Mesh[nb].cell_status[i+1][j] == CELL_STATUS_ACTIVE) {
	    P.convert(Local_SolnBlk[nb].Grid.nodeSW(i+1,j).X+Vector2D(-dX.x,-dX.y),
		      Local_SolnBlk[nb].Grid.nodeSE(i+1,j).X,
		      Local_SolnBlk[nb].Grid.nodeNE(i+1,j).X,
		      Local_SolnBlk[nb].Grid.nodeNW(i+1,j).X+Vector2D(-dX.x,dX.y));
	    if (P.point_in_polygon(Interface_Component_List[Ni].Spline.Xp[np])) {
	      dX2 = Interface_Component_List[Ni].Spline.Xp[np] - Local_SolnBlk[nb].Grid.Cell[i+1][j].Xc;
	      Wb = Local_SolnBlk[nb].W[i+1][j];
// 	      Wb = Local_SolnBlk[nb].W[i+1][j] + (Local_SolnBlk[nb].phi[i+1][j]^Local_SolnBlk[nb].dWdx[i+1][j])*dX2.x +
// 	                                         (Local_SolnBlk[nb].phi[i+1][j]^Local_SolnBlk[nb].dWdy[i+1][j])*dX2.y;
 	      Interface_Component_List[Ni].F[np] = Wb.burningrate()*norm_dir;//Local_SolnBlk[nb].Grid.nfaceE(i,j);
 	      //Interface_Component_List[Ni].F[np] = Local_SolnBlk[nb].U[i+1][j].burningrate()*Local_SolnBlk[nb].Grid.nfaceE(i,j);
	    }
	    P.deallocate();
	  }
	  if (Mesh[nb].cell_status[i-1][j] == CELL_STATUS_ACTIVE) {
	    P.convert(Local_SolnBlk[nb].Grid.nodeSW(i-1,j).X,
		      Local_SolnBlk[nb].Grid.nodeSE(i-1,j).X+Vector2D(dX.x,-dX.y),
		      Local_SolnBlk[nb].Grid.nodeNE(i-1,j).X+Vector2D(dX.x,dX.y),
		      Local_SolnBlk[nb].Grid.nodeNW(i-1,j).X);
	    if (P.point_in_polygon(Interface_Component_List[Ni].Spline.Xp[np])) {
	      dX2 = Interface_Component_List[Ni].Spline.Xp[np] - Local_SolnBlk[nb].Grid.Cell[i-1][j].Xc;
	      Wb = Local_SolnBlk[nb].W[i-1][j];
// 	      Wb = Local_SolnBlk[nb].W[i-1][j] + (Local_SolnBlk[nb].phi[i-1][j]^Local_SolnBlk[nb].dWdx[i-1][j])*dX2.x +
// 	                                         (Local_SolnBlk[nb].phi[i-1][j]^Local_SolnBlk[nb].dWdy[i-1][j])*dX2.y;
 	      Interface_Component_List[Ni].F[np] = Wb.burningrate()*norm_dir;//Local_SolnBlk[nb].Grid.nfaceW(i,j);
 	      //Interface_Component_List[Ni].F[np] = Local_SolnBlk[nb].U[i-1][j].burningrate()*Local_SolnBlk[nb].Grid.nfaceW(i,j);;
	    }
	    P.deallocate();
	  }

	}

      }

    }
  }

  // Interface speed function successfully computed on the solution block.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Redistribute_Solution_Content --             *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
inline int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Redistribute_Solution_Content(const double &Time) {

  int error_flag;

  for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
    if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      error_flag = Redistribute_Solution_Content(nb,Time);
      if (error_flag) return error_flag;
    }
  }

  // Solution content successfully redistributed.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Redistribute_Solution_Content --             *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
inline int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Redistribute_Solution_Content(const int &nb, const double &Time) {

  // Reset the solution state to atmospheric.  Used for debugging purposes.
  if (IP->Interface_IP.Redistribution_Type == SOLUTION_REDISTRIBUTION_NONE) {
    for (int j = Local_SolnBlk[nb].JCl-Local_SolnBlk[nb].Nghost; j <= Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost; j++) {
      for (int i = Local_SolnBlk[nb].ICl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].ICu+Local_SolnBlk[nb].Nghost; i++) {
	Local_SolnBlk[nb].W[i][j].Standard_Atmosphere();
	Local_SolnBlk[nb].U[i][j].Standard_Atmosphere();
      }
    }
  }

  // Exit immediately if no interfaces are present in the current
  // solution block.
  if (!Adjustment_Data[nb].Interface_Present[0]) return 0;

  int npts;
  Polygon P, Po;
  double A, A_total;

  // Locate deactivated cells and redistribute the solution content
  // accordingly.
  for (int j = Local_SolnBlk[nb].JCl; j <= Local_SolnBlk[nb].JCu; j++) {
    for (int i = Local_SolnBlk[nb].ICl; i <= Local_SolnBlk[nb].ICu; i++) {

      // If cell (i,j) has been deactivated due to the motion of the 
      // interface then inject the solution information based on area
      // weighting into the active neighbour cells.
      if (Mesh[nb].cell_status[i][j] != CELL_STATUS_ACTIVE &&
	  AMesh[nb].cell_status[i][j] == CELL_STATUS_ACTIVE) {
	npts = 0;
	if (Mesh[nb].cell_status[i][j-1]  == CELL_STATUS_ACTIVE &&
	    AMesh[nb].cell_status[i][j-1] == CELL_STATUS_ACTIVE) {
	  i_index[npts] = i; j_index[npts] = j-1; npts++;
	}
	if (Mesh[nb].cell_status[i-1][j]  == CELL_STATUS_ACTIVE &&
	    AMesh[nb].cell_status[i-1][j] == CELL_STATUS_ACTIVE) {
	  i_index[npts] = i-1; j_index[npts] = j; npts++;	  
	}
	if (Mesh[nb].cell_status[i+1][j]  == CELL_STATUS_ACTIVE &&
	    AMesh[nb].cell_status[i+1][j] == CELL_STATUS_ACTIVE) {
	  i_index[npts] = i+1; j_index[npts] = j; npts++;
	}
	if (Mesh[nb].cell_status[i][j+1]  == CELL_STATUS_ACTIVE &&
	    AMesh[nb].cell_status[i][j+1] == CELL_STATUS_ACTIVE) {
	  i_index[npts] = i; j_index[npts] = j+1; npts++;
	}
	if (!npts) {
	  // Not possible.  Either the time-step was too large or
	  // there is a difference in refinement on the embedded
	  // boundary.
	  cout << endl << " DEACTIVATED::npts = 0 @ Xc(nb,i,j) = Xc("
	       << Local_Solution_Block_List->Block[nb].gblknum << "," << i << "," << j << ") ="
	       << Local_SolnBlk[nb].Grid.Cell[i][j].Xc; cout.flush();
	  return 9900;
	} else {
	  // Distribute mass into contributing cells.
	  for (int n = 0; n < npts; n++) {
	    Local_SolnBlk[nb].U[i_index[n]][j_index[n]] = ((Local_SolnBlk[nb].U[i_index[n]][j_index[n]]*AGrid[nb].Cell[i_index[n]][j_index[n]].A +
							    Local_SolnBlk[nb].U[i][j]*AGrid[nb].Cell[i][j].A)/
							   (AGrid[nb].Cell[i_index[n]][j_index[n]].A + AGrid[nb].Cell[i][j].A));
	    if (Local_SolnBlk[nb].W[i][j].Unphysical_Properties()) {
	      cout << endl << " DEACTIVATED::npts = 0 @ Xc(nb,i,j) = Xc("
		   << Local_Solution_Block_List->Block[nb].gblknum << "," << i << "," << j << ") ="
		   << Local_SolnBlk[nb].Grid.Cell[i][j].Xc; cout.flush();
	      return 9901;
	    }
	    Local_SolnBlk[nb].W[i_index[n]][j_index[n]] = W(Local_SolnBlk[nb].U[i_index[n]][j_index[n]]);
	  }
	}
      }
    }
  }

  // Conservative re-mapping of average solution content in active
  // cells that may have had an edge relocated during the mesh
  // adjustment process.
  for (int j = Local_SolnBlk[nb].JCl; j <= Local_SolnBlk[nb].JCu; j++) {
    for (int i = Local_SolnBlk[nb].ICl; i <= Local_SolnBlk[nb].ICu; i++) {

      if (Mesh[nb].cell_status[i][j] == CELL_STATUS_ACTIVE &&
	  (AMesh[nb].cell_status[i-1][j-1] != CELL_STATUS_ACTIVE ||
	   AMesh[nb].cell_status[i  ][j-1] != CELL_STATUS_ACTIVE ||
	   AMesh[nb].cell_status[i+1][j-1] != CELL_STATUS_ACTIVE ||
	   AMesh[nb].cell_status[i-1][j  ] != CELL_STATUS_ACTIVE ||
	   AMesh[nb].cell_status[i  ][j  ] != CELL_STATUS_ACTIVE ||
	   AMesh[nb].cell_status[i+1][j  ] != CELL_STATUS_ACTIVE ||
	   AMesh[nb].cell_status[i-1][j+1] != CELL_STATUS_ACTIVE ||
	   AMesh[nb].cell_status[i  ][j+1] != CELL_STATUS_ACTIVE ||
	   AMesh[nb].cell_status[i+1][j+1] != CELL_STATUS_ACTIVE)) {
	A_total = ZERO;
	Local_SolnBlk[nb].dUdt[i][j][0].Vacuum();
	// Create a polygom representing the current cell.
	P.convert(Local_SolnBlk[nb].Grid.nodeSW(i,j).X,Local_SolnBlk[nb].Grid.nodeSE(i,j).X,
		  Local_SolnBlk[nb].Grid.nodeNE(i,j).X,Local_SolnBlk[nb].Grid.nodeNW(i,j).X);
	// Search the cells of the previous mesh adjustment for
	// intersecting cells.
	for (int jj = j-1; jj < j+2; jj++) {
	  for (int ii = i-1; ii < i+2; ii++) {
	    if (AMesh[nb].cell_status[ii][jj] == CELL_STATUS_ACTIVE) {
	      // Create a polygon from the cell from the previous
	      // adjustment.
	      Po.convert(AGrid[nb].nodeSW(ii,jj).X,AGrid[nb].nodeSE(ii,jj).X,
			 AGrid[nb].nodeNE(ii,jj).X,AGrid[nb].nodeNW(ii,jj).X);
	      // Determine the area of intersection between the
	      // fine and coarse cell polygons.
	      A = Polygon_Intersection_Area(P,Po);
	      Po.deallocate();
	      // Add the area of intersection to the total area
	      // of intersection.
	      A_total += A;
	      // Restrict fine cell solution to the coarse cell
	      // based on the area of intersection.
	      Local_SolnBlk[nb].dUdt[i][j][0] += Local_SolnBlk[nb].U[ii][jj]*A;
	    }
	  }
	}
	P.deallocate();
	// Complete solution restriction by dividing the restricted
	// solution by the total area of intersection.
	if (A_total < TOLER*Local_SolnBlk[nb].Grid.Cell[i][j].A) {
	  cout << endl << " time = " << Time; cout.flush();
	  cout << endl << " X(" << Local_Solution_Block_List->Block[nb].gblknum << ","
	       << i << "," << j << ") =" << Local_SolnBlk[nb].Grid.Cell[i][j].Xc; cout.flush();
	  cout << endl << " A(" << Local_Solution_Block_List->Block[nb].gblknum << ","
	       << i << "," << j << ") = " << Local_SolnBlk[nb].Grid.Cell[i][j].A; cout.flush();
	  cout << endl << " A_total = " << A_total; cout.flush();
	  A_total = ZERO;
	  cout << endl << setprecision(24) << Local_SolnBlk[nb].Grid.nodeSW(i,j); cout.flush();
	  cout << endl << setprecision(24) << Local_SolnBlk[nb].Grid.nodeSE(i,j); cout.flush();
	  cout << endl << setprecision(24) << Local_SolnBlk[nb].Grid.nodeNE(i,j); cout.flush();
	  cout << endl << setprecision(24) << Local_SolnBlk[nb].Grid.nodeNW(i,j); cout.flush();
	  P.convert(Local_SolnBlk[nb].Grid.nodeSW(i,j).X,Local_SolnBlk[nb].Grid.nodeSE(i,j).X,
		    Local_SolnBlk[nb].Grid.nodeNE(i,j).X,Local_SolnBlk[nb].Grid.nodeNW(i,j).X);
	  cout << endl << " P =" << setprecision(14) << P; cout.flush();
	  Polygon Pt;
	  for (int jj = j-1; jj < j+2; jj++) {
	    for (int ii = i-1; ii < i+2; ii++) {
	      if (AMesh[nb].cell_status[ii][jj] == CELL_STATUS_ACTIVE) {
		cout << endl;
		cout << endl << setprecision(24) << AGrid[nb].nodeSW(ii,jj).X; cout.flush();
		cout << endl << setprecision(24) << AGrid[nb].nodeSE(ii,jj).X; cout.flush();
		cout << endl << setprecision(24) << AGrid[nb].nodeNE(ii,jj).X; cout.flush();
		cout << endl << setprecision(24) << AGrid[nb].nodeNW(ii,jj).X; cout.flush();
		Po.convert(AGrid[nb].nodeSW(ii,jj).X,AGrid[nb].nodeSE(ii,jj).X,
			   AGrid[nb].nodeNE(ii,jj).X,AGrid[nb].nodeNW(ii,jj).X);
		cout << endl << " Po =" << setprecision(14) << Po; cout.flush();
		Polygon_Intersection(P,Po,Pt);
		if (Pt.np) { cout << endl << " Pt =" << setprecision(14) << Pt; cout.flush(); }
		Pt.deallocate();
		A = Polygon_Intersection_Area(P,Po);
		cout << endl << " A = " << setprecision(14) << A; cout.flush();
		Po.deallocate();
		A_total += A;
		cout << endl << " A_total = " << setprecision(14) << A_total; cout.flush();
	      }
	    }
	  }
	  P.deallocate();
	  return 9323;
	}
	Local_SolnBlk[nb].U[i][j] = Local_SolnBlk[nb].dUdt[i][j][0]/A_total;
	if (Local_SolnBlk[nb].U[i][j].Unphysical_Properties()) return 9767;
	Local_SolnBlk[nb].W[i][j] = W(Local_SolnBlk[nb].U[i][j]);
	Local_SolnBlk[nb].dUdt[i][j][0].Vacuum();
      }

    }
  }

  // Solution content successfully redistributed.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Reset_Interface_Motion_Type --               *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
inline int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Reset_Interface_Motion_Type(void) {

  // Exit immediately if no embedded boundaries exist.
  if (!Interface_Component_List.Ni) return 1;

  // Reset the interface velocity accordingly.
  if (Interface_Component_List[1].Type == INTERFACE_NACA0012_AEROFOIL ||
      Interface_Component_List[1].Type == INTERFACE_NACA0015_AEROFOIL) {

    // Reset the interface motion type to rotate if the embedded
    // boundary is a NACA0012 or NACA0015 aerofoil.
    Interface_Component_List[1].Motion = INTERFACE_MOTION_ROTATE;

  } else if (Interface_Component_List.Ni == 17) {

    // Reset the interface motion types to rotate for the pins of the
    // branched duct.
    for (int ni = 1; ni <= Interface_Component_List.Ni; ni++) {
      //      Interface_Component_List[ni].Motion = IP->Interface_IP.Component_List[ni].Motion;
      //      Interface_Component_List[ni].Speed = IP->Interface_IP.Component_List[ni].Speed;
      if(ni >= 4 && ni <= 8) {
	Interface_Component_List[ni].Motion = INTERFACE_MOTION_ROTATE;
	Interface_Component_List[ni].Speed = Vector2D(9.500000e-03,1.000000e-02);
      }
      if(ni >= 13 && ni <= 17) {
	Interface_Component_List[ni].Motion = INTERFACE_MOTION_ROTATE;
	Interface_Component_List[ni].Speed = Vector2D(-9.500000e-03,1.000000e-02);
      }
      Interface_Component_List[ni].Set_Velocity_Function_Type();
      //      cout << Interface_Component_List[ni].Motion << " " << Interface_Component_List[ni].Speed << endl;
    }

  } else if (IP->i_Grid == GRID_ROCKET_MOTOR) {

  } else {

    // Invalid configuration for reset.
    return 1;

  }

  // Ensure that a moving-boundary Riemann problem is being used.
  if (IP->i_Flux_Function != FLUX_FUNCTION_GODUNOV_MB &&
      IP->i_Flux_Function != FLUX_FUNCTION_ROE_MB &&
      IP->i_Flux_Function != FLUX_FUNCTION_HLLE_MB &&
      IP->i_Flux_Function != FLUX_FUNCTION_VANLEER_MB) return 1;

//   // Reset the boundary conditions to BC_CHARACTERISTIC.
//   for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
//     if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {
//       if (Local_SolnBlk[nb].Grid.BCtypeN[i] != BC_NONE)
// 	Local_SolnBlk[nb].Grid.BndNorthSpline.setBCtype(BC_CHARACTERISTIC);
//       if (Local_SolnBlk[nb].Grid.BCtypeS[i] != BC_NONE)
// 	Local_SolnBlk[nb].Grid.BndSouthSpline.setBCtype(BC_CHARACTERISTIC);
//       if (Local_SolnBlk[nb].Grid.BCtypeE[i] != BC_NONE)
// 	Local_SolnBlk[nb].Grid.BndEastSpline.setBCtype(BC_CHARACTERISTIC);
//       if (Local_SolnBlk[nb].Grid.BCtypeW[i] != BC_NONE)
// 	Local_SolnBlk[nb].Grid.BndWestSpline.setBCtype(BC_CHARACTERISTIC);
//       Set_BCs(Local_SolnBlk[nb].Grid);
//       Update_Exterior_Nodes(Local_SolnBlk[nb].Grid);
//       Update_Cells(Local_SolnBlk[nb].Grid);
//     }
//   }

  // Interface motion type reset successfully.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Zero_Interface_Velocity_Function --          *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
inline int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Zero_Interface_Velocity_Function(void) {

  // Exit immediately if no embedded boundaries exist.
  if (!Interface_Component_List.Ni) return 1;

  // Reset the velocity function for each interface.
  for (int ni = 1; ni <= Interface_Component_List.Ni; ni++) {
    Interface_Component_List[ni].Zero_Velocity_Function();
  }

  for (int ni = 1; ni <= Interface_Union_List.Ni; ni++) {
    Interface_Union_List[ni].Zero_Velocity_Function();
  }

  // Interface velocity function(s) set to zero.
  return 0;

}

/**********************************************************************
 **********************************************************************
 ** Bilinear interplation routines.                                  **
 **********************************************************************
 **********************************************************************/

/**********************************************************************
 * EmbeddedBoundaries2D::Wn -- Node primitive variable solution state.*
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
inline pState EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Wn(const int &nb, const int &ii, const int &jj) {
  double wo, wi;
  Vector2D lambda, R;
  Tensor2D I;
  pState Wo;
  // Determine the contribution to the weighting coefficients for the 
  // standard interpolation support set:
  R = Local_SolnBlk[nb].Grid.Cell[ii-1][jj-1].Xc + Local_SolnBlk[nb].Grid.Cell[ii  ][jj-1].Xc +
      Local_SolnBlk[nb].Grid.Cell[ii-1][jj  ].Xc + Local_SolnBlk[nb].Grid.Cell[ii  ][jj  ].Xc - FOUR*Local_SolnBlk[nb].Grid.Node[ii][jj].X;
  I.xx = sqr(Local_SolnBlk[nb].Grid.Cell[ii-1][jj-1].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x) +
         sqr(Local_SolnBlk[nb].Grid.Cell[ii  ][jj-1].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x) +
         sqr(Local_SolnBlk[nb].Grid.Cell[ii-1][jj  ].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x) +
         sqr(Local_SolnBlk[nb].Grid.Cell[ii  ][jj  ].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x);
  I.xy = (Local_SolnBlk[nb].Grid.Cell[ii-1][jj-1].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x)*
         (Local_SolnBlk[nb].Grid.Cell[ii-1][jj-1].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y) +
         (Local_SolnBlk[nb].Grid.Cell[ii  ][jj-1].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x)*
         (Local_SolnBlk[nb].Grid.Cell[ii  ][jj-1].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y) +
         (Local_SolnBlk[nb].Grid.Cell[ii-1][jj  ].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x)*
         (Local_SolnBlk[nb].Grid.Cell[ii-1][jj  ].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y) +
         (Local_SolnBlk[nb].Grid.Cell[ii  ][jj  ].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x)*
         (Local_SolnBlk[nb].Grid.Cell[ii  ][jj  ].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y);
  I.yy = sqr(Local_SolnBlk[nb].Grid.Cell[ii-1][jj-1].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y) +
         sqr(Local_SolnBlk[nb].Grid.Cell[ii  ][jj-1].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y) +
         sqr(Local_SolnBlk[nb].Grid.Cell[ii-1][jj  ].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y) +
         sqr(Local_SolnBlk[nb].Grid.Cell[ii  ][jj  ].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y);
  // Determine the contribution to the weighting coefficients for the 
  // distance-2 WEST pair of cells if requried:
  if (Local_SolnBlk[nb].Grid.lfaceN(ii-1,jj-1) < NANO && ii-2 >= 0) {
    R += Local_SolnBlk[nb].Grid.Cell[ii-2][jj-1].Xc + Local_SolnBlk[nb].Grid.Cell[ii-2][jj  ].Xc - TWO*Local_SolnBlk[nb].Grid.Node[ii][jj].X;
    I.xx += sqr(Local_SolnBlk[nb].Grid.Cell[ii-2][jj-1].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x) +
            sqr(Local_SolnBlk[nb].Grid.Cell[ii-2][jj  ].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x);
    I.xy += (Local_SolnBlk[nb].Grid.Cell[ii-2][jj-1].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x)*
            (Local_SolnBlk[nb].Grid.Cell[ii-2][jj-1].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y) +
            (Local_SolnBlk[nb].Grid.Cell[ii-2][jj  ].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x)*
            (Local_SolnBlk[nb].Grid.Cell[ii-2][jj  ].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y);
    I.yy += sqr(Local_SolnBlk[nb].Grid.Cell[ii-2][jj-1].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y) +
            sqr(Local_SolnBlk[nb].Grid.Cell[ii-2][jj  ].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y);
  }
  // Determine the contribution to the weighting coefficients for the 
  // distance-2 EAST pair of cells if requried:
  if (Local_SolnBlk[nb].Grid.lfaceN(ii,jj-1) < NANO && ii+1 <= Local_SolnBlk[nb].ICu+Local_SolnBlk[nb].Nghost) {
    R += Local_SolnBlk[nb].Grid.Cell[ii+1][jj-1].Xc + Local_SolnBlk[nb].Grid.Cell[ii+1][jj  ].Xc - TWO*Local_SolnBlk[nb].Grid.Node[ii][jj].X;
    I.xx += sqr(Local_SolnBlk[nb].Grid.Cell[ii+1][jj-1].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x) +
            sqr(Local_SolnBlk[nb].Grid.Cell[ii+1][jj  ].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x);
    I.xy += (Local_SolnBlk[nb].Grid.Cell[ii+1][jj-1].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x)*
            (Local_SolnBlk[nb].Grid.Cell[ii+1][jj-1].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y) +
            (Local_SolnBlk[nb].Grid.Cell[ii+1][jj  ].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x)*
            (Local_SolnBlk[nb].Grid.Cell[ii+1][jj  ].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y);
    I.yy += sqr(Local_SolnBlk[nb].Grid.Cell[ii+1][jj-1].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y) +
            sqr(Local_SolnBlk[nb].Grid.Cell[ii+1][jj  ].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y);
  }
  // Determine the contribution to the weighting coefficients for the 
  // distance-2 SOUTH pair of cells if requried:
  if (Local_SolnBlk[nb].Grid.lfaceE(ii-1,jj-1) < NANO && jj-2 >= 0) {
    R += Local_SolnBlk[nb].Grid.Cell[ii-1][jj-2].Xc + Local_SolnBlk[nb].Grid.Cell[ii  ][jj-2].Xc - TWO*Local_SolnBlk[nb].Grid.Node[ii][jj].X;
    I.xx += sqr(Local_SolnBlk[nb].Grid.Cell[ii-1][jj-2].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x) +
            sqr(Local_SolnBlk[nb].Grid.Cell[ii  ][jj-2].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x);
    I.xy += (Local_SolnBlk[nb].Grid.Cell[ii-1][jj-2].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x)*
            (Local_SolnBlk[nb].Grid.Cell[ii-1][jj-2].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y) +
            (Local_SolnBlk[nb].Grid.Cell[ii  ][jj-2].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x)*
            (Local_SolnBlk[nb].Grid.Cell[ii  ][jj-2].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y);
    I.yy += sqr(Local_SolnBlk[nb].Grid.Cell[ii-1][jj-2].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y) +
            sqr(Local_SolnBlk[nb].Grid.Cell[ii  ][jj-2].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y);
  }
  // Determine the contribution to the weighting coefficients for the 
  // distance-2 NORTH pair of cells if requried:
  if (Local_SolnBlk[nb].Grid.lfaceE(ii-1,jj) < NANO && jj+1 <= Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost) {
    R += Local_SolnBlk[nb].Grid.Cell[ii-1][jj+1].Xc + Local_SolnBlk[nb].Grid.Cell[ii  ][jj+1].Xc - TWO*Local_SolnBlk[nb].Grid.Node[ii][jj].X;
    I.xx += sqr(Local_SolnBlk[nb].Grid.Cell[ii-1][jj+1].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x) +
            sqr(Local_SolnBlk[nb].Grid.Cell[ii  ][jj+1].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x);
    I.xy += (Local_SolnBlk[nb].Grid.Cell[ii-1][jj+1].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x)*
            (Local_SolnBlk[nb].Grid.Cell[ii-1][jj+1].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y) +
            (Local_SolnBlk[nb].Grid.Cell[ii  ][jj+1].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x)*
            (Local_SolnBlk[nb].Grid.Cell[ii  ][jj+1].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y);
    I.yy += sqr(Local_SolnBlk[nb].Grid.Cell[ii-1][jj+1].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y) +
            sqr(Local_SolnBlk[nb].Grid.Cell[ii  ][jj+1].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y);
  }
  // Determine the weighting coefficients:
  lambda.x = (I.xy*R.y - I.yy*R.x)/(I.xx*I.yy - I.xy*I.xy);
  lambda.y = (I.xy*R.x - I.xx*R.y)/(I.xx*I.yy - I.xy*I.xy);
  // Determine the weighted interpolated state from the standard
  // interpolation support set.
  wi = ONE + lambda*(Local_SolnBlk[nb].Grid.Cell[ii-1][jj-1].Xc - Local_SolnBlk[nb].Grid.Node[ii][jj].X);
  wo = wi;
  Wo = wi*Local_SolnBlk[nb].W[ii-1][jj-1]; 
  wi = ONE + lambda*(Local_SolnBlk[nb].Grid.Cell[ii  ][jj-1].Xc - Local_SolnBlk[nb].Grid.Node[ii][jj].X);
  wo += wi;
  Wo += wi*Local_SolnBlk[nb].W[ii  ][jj-1]; 
  wi = ONE + lambda*(Local_SolnBlk[nb].Grid.Cell[ii-1][jj  ].Xc - Local_SolnBlk[nb].Grid.Node[ii][jj].X);
  wo += wi;
  Wo += wi*Local_SolnBlk[nb].W[ii-1][jj  ]; 
  wi = ONE + lambda*(Local_SolnBlk[nb].Grid.Cell[ii  ][jj  ].Xc - Local_SolnBlk[nb].Grid.Node[ii][jj].X);
  wo += wi;
  Wo += wi*Local_SolnBlk[nb].W[ii  ][jj  ]; 
  // Determine the weighted interpolated state from the distance-2 WEST
  // pair of cells if required:
  if (Local_SolnBlk[nb].Grid.lfaceN(ii-1,jj-1) < NANO && ii-2 >= 0) {
    wi = ONE + lambda*(Local_SolnBlk[nb].Grid.Cell[ii-2][jj-1].Xc - Local_SolnBlk[nb].Grid.Node[ii][jj].X);
    wo += wi;
    Wo += wi*Local_SolnBlk[nb].W[ii-2][jj-1];
    wi = ONE + lambda*(Local_SolnBlk[nb].Grid.Cell[ii-2][jj  ].Xc - Local_SolnBlk[nb].Grid.Node[ii][jj].X);
    wo += wi;
    Wo += wi*Local_SolnBlk[nb].W[ii-2][jj  ];
  }
  // Determine the weighted interpolated state from the distance-2 EAST
  // pair of cells if required:
  if (Local_SolnBlk[nb].Grid.lfaceN(ii,jj-1) < NANO && ii+1 <= Local_SolnBlk[nb].ICu+Local_SolnBlk[nb].Nghost) {
    wi = ONE + lambda*(Local_SolnBlk[nb].Grid.Cell[ii+1][jj-1].Xc - Local_SolnBlk[nb].Grid.Node[ii][jj].X);
    wo += wi;
    Wo += wi*Local_SolnBlk[nb].W[ii+1][jj-1];
    wi = ONE + lambda*(Local_SolnBlk[nb].Grid.Cell[ii+1][jj  ].Xc - Local_SolnBlk[nb].Grid.Node[ii][jj].X);
    wo += wi;
    Wo += wi*Local_SolnBlk[nb].W[ii+1][jj  ];
  }
  // Determine the weighted interpolated state from the distance-2 SOUTH
  // pair of cells if required:
  if (Local_SolnBlk[nb].Grid.lfaceE(ii-1,jj-1) < NANO && jj-2 >= 0) {
    wi = ONE + lambda*(Local_SolnBlk[nb].Grid.Cell[ii-1][jj-2].Xc - Local_SolnBlk[nb].Grid.Node[ii][jj].X);
    wo += wi;
    Wo += wi*Local_SolnBlk[nb].W[ii-1][jj-2];
    wi = ONE + lambda*(Local_SolnBlk[nb].Grid.Cell[ii  ][jj-2].Xc - Local_SolnBlk[nb].Grid.Node[ii][jj].X);
    wo += wi;
    Wo += wi*Local_SolnBlk[nb].W[ii  ][jj-2];
  }
  // Determine the weighted interpolated state from the distance-2 NORTH
  // pair of cells if required:
  if (Local_SolnBlk[nb].Grid.lfaceE(ii-1,jj) < NANO && jj+1 <= Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost) {
    wi = ONE + lambda*(Local_SolnBlk[nb].Grid.Cell[ii-1][jj+1].Xc - Local_SolnBlk[nb].Grid.Node[ii][jj].X);
    wo += wi;
    Wo += wi*Local_SolnBlk[nb].W[ii-1][jj+1];
    wi = ONE + lambda*(Local_SolnBlk[nb].Grid.Cell[ii  ][jj+1].Xc - Local_SolnBlk[nb].Grid.Node[ii][jj].X);
    wo += wi;
    Wo += wi*Local_SolnBlk[nb].W[ii  ][jj+1];
  }
  // Return interpolated primitive node solution.
  return Wo/wo;
}

/**********************************************************************
 * EmbeddedBoundaries2D::Un -- Node conserved variable solution state.*
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
inline cState EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Un(const int &nb, const int &ii, const int &jj) {
  double wo, wi;
  Vector2D lambda, R;
  Tensor2D I;
  cState Uo;
  // Determine the contribution to the weighting coefficients for the 
  // standard interpolation support set:
  R = Local_SolnBlk[nb].Grid.Cell[ii-1][jj-1].Xc + Local_SolnBlk[nb].Grid.Cell[ii  ][jj-1].Xc +
      Local_SolnBlk[nb].Grid.Cell[ii-1][jj  ].Xc + Local_SolnBlk[nb].Grid.Cell[ii  ][jj  ].Xc - FOUR*Local_SolnBlk[nb].Grid.Node[ii][jj].X;
  I.xx = sqr(Local_SolnBlk[nb].Grid.Cell[ii-1][jj-1].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x) +
         sqr(Local_SolnBlk[nb].Grid.Cell[ii  ][jj-1].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x) +
         sqr(Local_SolnBlk[nb].Grid.Cell[ii-1][jj  ].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x) +
         sqr(Local_SolnBlk[nb].Grid.Cell[ii  ][jj  ].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x);
  I.xy = (Local_SolnBlk[nb].Grid.Cell[ii-1][jj-1].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x)*
         (Local_SolnBlk[nb].Grid.Cell[ii-1][jj-1].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y) +
         (Local_SolnBlk[nb].Grid.Cell[ii  ][jj-1].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x)*
         (Local_SolnBlk[nb].Grid.Cell[ii  ][jj-1].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y) +
         (Local_SolnBlk[nb].Grid.Cell[ii-1][jj  ].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x)*
         (Local_SolnBlk[nb].Grid.Cell[ii-1][jj  ].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y) +
         (Local_SolnBlk[nb].Grid.Cell[ii  ][jj  ].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x)*
         (Local_SolnBlk[nb].Grid.Cell[ii  ][jj  ].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y);
  I.yy = sqr(Local_SolnBlk[nb].Grid.Cell[ii-1][jj-1].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y) +
         sqr(Local_SolnBlk[nb].Grid.Cell[ii  ][jj-1].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y) +
         sqr(Local_SolnBlk[nb].Grid.Cell[ii-1][jj  ].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y) +
         sqr(Local_SolnBlk[nb].Grid.Cell[ii  ][jj  ].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y);
  // Determine the contribution to the weighting coefficients for the 
  // distance-2 WEST pair of cells if requried:
  if (Local_SolnBlk[nb].Grid.lfaceN(ii-1,jj-1) < NANO && ii-2 >= 0) {
    R += Local_SolnBlk[nb].Grid.Cell[ii-2][jj-1].Xc + Local_SolnBlk[nb].Grid.Cell[ii-2][jj  ].Xc - TWO*Local_SolnBlk[nb].Grid.Node[ii][jj].X;
    I.xx += sqr(Local_SolnBlk[nb].Grid.Cell[ii-2][jj-1].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x) +
            sqr(Local_SolnBlk[nb].Grid.Cell[ii-2][jj  ].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x);
    I.xy += (Local_SolnBlk[nb].Grid.Cell[ii-2][jj-1].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x)*
            (Local_SolnBlk[nb].Grid.Cell[ii-2][jj-1].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y) +
            (Local_SolnBlk[nb].Grid.Cell[ii-2][jj  ].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x)*
            (Local_SolnBlk[nb].Grid.Cell[ii-2][jj  ].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y);
    I.yy += sqr(Local_SolnBlk[nb].Grid.Cell[ii-2][jj-1].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y) +
            sqr(Local_SolnBlk[nb].Grid.Cell[ii-2][jj  ].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y);
  }
  // Determine the contribution to the weighting coefficients for the 
  // distance-2 EAST pair of cells if requried:
  if (Local_SolnBlk[nb].Grid.lfaceN(ii,jj-1) < NANO && ii+1 <= Local_SolnBlk[nb].ICu+Local_SolnBlk[nb].Nghost) {
    R += Local_SolnBlk[nb].Grid.Cell[ii+1][jj-1].Xc + Local_SolnBlk[nb].Grid.Cell[ii+1][jj  ].Xc - TWO*Local_SolnBlk[nb].Grid.Node[ii][jj].X;
    I.xx += sqr(Local_SolnBlk[nb].Grid.Cell[ii+1][jj-1].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x) +
            sqr(Local_SolnBlk[nb].Grid.Cell[ii+1][jj  ].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x);
    I.xy += (Local_SolnBlk[nb].Grid.Cell[ii+1][jj-1].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x)*
            (Local_SolnBlk[nb].Grid.Cell[ii+1][jj-1].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y) +
            (Local_SolnBlk[nb].Grid.Cell[ii+1][jj  ].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x)*
            (Local_SolnBlk[nb].Grid.Cell[ii+1][jj  ].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y);
    I.yy += sqr(Local_SolnBlk[nb].Grid.Cell[ii+1][jj-1].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y) +
            sqr(Local_SolnBlk[nb].Grid.Cell[ii+1][jj  ].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y);
  }
  // Determine the contribution to the weighting coefficients for the 
  // distance-2 SOUTH pair of cells if requried:
  if (Local_SolnBlk[nb].Grid.lfaceE(ii-1,jj-1) < NANO && jj-2 >= 0) {
    R += Local_SolnBlk[nb].Grid.Cell[ii-1][jj-2].Xc + Local_SolnBlk[nb].Grid.Cell[ii  ][jj-2].Xc - TWO*Local_SolnBlk[nb].Grid.Node[ii][jj].X;
    I.xx += sqr(Local_SolnBlk[nb].Grid.Cell[ii-1][jj-2].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x) +
            sqr(Local_SolnBlk[nb].Grid.Cell[ii  ][jj-2].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x);
    I.xy += (Local_SolnBlk[nb].Grid.Cell[ii-1][jj-2].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x)*
            (Local_SolnBlk[nb].Grid.Cell[ii-1][jj-2].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y) +
            (Local_SolnBlk[nb].Grid.Cell[ii  ][jj-2].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x)*
            (Local_SolnBlk[nb].Grid.Cell[ii  ][jj-2].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y);
    I.yy += sqr(Local_SolnBlk[nb].Grid.Cell[ii-1][jj-2].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y) +
            sqr(Local_SolnBlk[nb].Grid.Cell[ii  ][jj-2].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y);
  }
  // Determine the contribution to the weighting coefficients for the 
  // distance-2 NORTH pair of cells if requried:
  if (Local_SolnBlk[nb].Grid.lfaceE(ii-1,jj) < NANO && jj+1 <= Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost) {
    R += Local_SolnBlk[nb].Grid.Cell[ii-1][jj+1].Xc + Local_SolnBlk[nb].Grid.Cell[ii  ][jj+1].Xc - TWO*Local_SolnBlk[nb].Grid.Node[ii][jj].X;
    I.xx += sqr(Local_SolnBlk[nb].Grid.Cell[ii-1][jj+1].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x) +
            sqr(Local_SolnBlk[nb].Grid.Cell[ii  ][jj+1].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x);
    I.xy += (Local_SolnBlk[nb].Grid.Cell[ii-1][jj+1].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x)*
            (Local_SolnBlk[nb].Grid.Cell[ii-1][jj+1].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y) +
            (Local_SolnBlk[nb].Grid.Cell[ii  ][jj+1].Xc.x - Local_SolnBlk[nb].Grid.Node[ii][jj].X.x)*
            (Local_SolnBlk[nb].Grid.Cell[ii  ][jj+1].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y);
    I.yy += sqr(Local_SolnBlk[nb].Grid.Cell[ii-1][jj+1].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y) +
            sqr(Local_SolnBlk[nb].Grid.Cell[ii  ][jj+1].Xc.y - Local_SolnBlk[nb].Grid.Node[ii][jj].X.y);
  }
  // Determine the weighting coefficients:
  lambda.x = (I.xy*R.y - I.yy*R.x)/(I.xx*I.yy - I.xy*I.xy);
  lambda.y = (I.xy*R.x - I.xx*R.y)/(I.xx*I.yy - I.xy*I.xy);
  // Determine the weighted interpolated state from the standard
  // interpolation support set.
  wi = ONE + lambda*(Local_SolnBlk[nb].Grid.Cell[ii-1][jj-1].Xc - Local_SolnBlk[nb].Grid.Node[ii][jj].X);
  wo = wi;
  Uo = wi*Local_SolnBlk[nb].U[ii-1][jj-1];
  wi = ONE + lambda*(Local_SolnBlk[nb].Grid.Cell[ii  ][jj-1].Xc - Local_SolnBlk[nb].Grid.Node[ii][jj].X);
  wo += wi;
  Uo += wi*Local_SolnBlk[nb].U[ii  ][jj-1];
  wi = ONE + lambda*(Local_SolnBlk[nb].Grid.Cell[ii-1][jj  ].Xc - Local_SolnBlk[nb].Grid.Node[ii][jj].X);
  wo += wi;
  Uo += wi*Local_SolnBlk[nb].U[ii-1][jj  ];
  wi = ONE + lambda*(Local_SolnBlk[nb].Grid.Cell[ii  ][jj  ].Xc - Local_SolnBlk[nb].Grid.Node[ii][jj].X);
  wo += wi;
  Uo += wi*Local_SolnBlk[nb].U[ii  ][jj  ];
  // Determine the weighted interpolated state from the distance-2 WEST
  // pair of cells if required:
  if (Local_SolnBlk[nb].Grid.lfaceN(ii-1,jj-1) < NANO && ii-2 >= 0) {
    wi = ONE + lambda*(Local_SolnBlk[nb].Grid.Cell[ii-2][jj-1].Xc - Local_SolnBlk[nb].Grid.Node[ii][jj].X);
    wo += wi;
    Uo += wi*Local_SolnBlk[nb].U[ii-2][jj-1];
    wi = ONE + lambda*(Local_SolnBlk[nb].Grid.Cell[ii-2][jj  ].Xc - Local_SolnBlk[nb].Grid.Node[ii][jj].X);
    wo += wi;
    Uo += wi*Local_SolnBlk[nb].U[ii-2][jj  ];
  }
  // Determine the weighted interpolated state from the distance-2 EAST
  // pair of cells if required:
  if (Local_SolnBlk[nb].Grid.lfaceN(ii,jj-1) < NANO && ii+1 <= Local_SolnBlk[nb].ICu+Local_SolnBlk[nb].Nghost) {
    wi = ONE + lambda*(Local_SolnBlk[nb].Grid.Cell[ii+1][jj-1].Xc - Local_SolnBlk[nb].Grid.Node[ii][jj].X);
    wo += wi;
    Uo += wi*Local_SolnBlk[nb].U[ii+1][jj-1];
    wi = ONE + lambda*(Local_SolnBlk[nb].Grid.Cell[ii+1][jj  ].Xc - Local_SolnBlk[nb].Grid.Node[ii][jj].X);
    wo += wi;
    Uo += wi*Local_SolnBlk[nb].U[ii+1][jj  ];
  }
  // Determine the weighted interpolated state from the distance-2 SOUTH
  // pair of cells if required:
  if (Local_SolnBlk[nb].Grid.lfaceE(ii-1,jj-1) < NANO && jj-2 >= 0) {
    wi = ONE + lambda*(Local_SolnBlk[nb].Grid.Cell[ii-1][jj-2].Xc - Local_SolnBlk[nb].Grid.Node[ii][jj].X);
    wo += wi;
    Uo += wi*Local_SolnBlk[nb].U[ii-1][jj-2];
    wi = ONE + lambda*(Local_SolnBlk[nb].Grid.Cell[ii  ][jj-2].Xc - Local_SolnBlk[nb].Grid.Node[ii][jj].X);
    wo += wi;
    Uo += wi*Local_SolnBlk[nb].U[ii  ][jj-2];
  }
  // Determine the weighted interpolated state from the distance-2 NORTH
  // pair of cells if required:
  if (Local_SolnBlk[nb].Grid.lfaceE(ii-1,jj) < NANO && jj+1 <= Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost) {
    wi = ONE + lambda*(Local_SolnBlk[nb].Grid.Cell[ii-1][jj+1].Xc - Local_SolnBlk[nb].Grid.Node[ii][jj].X);
    wo += wi;
    Uo += wi*Local_SolnBlk[nb].U[ii-1][jj+1];
    wi = ONE + lambda*(Local_SolnBlk[nb].Grid.Cell[ii  ][jj+1].Xc - Local_SolnBlk[nb].Grid.Node[ii][jj].X);
    wo += wi;
    Uo += wi*Local_SolnBlk[nb].U[ii  ][jj+1];
  }
  // Return interpolated primitive node solution.
  return Uo/wo;
}

/**********************************************************************
 * EmbeddedBoundaries2D::Wn?? -- Node primitive solution state.       *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
inline pState EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
WnNW(const int &nb, const int &ii, const int &jj) {
  return Wn(nb,ii,jj+1);
}

template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
inline pState EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
WnNE(const int &nb, const int &ii, const int &jj) {
  return Wn(nb,ii+1,jj+1);
}

template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
inline pState EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
WnSE(const int &nb, const int &ii, const int &jj) {
  return Wn(nb,ii+1,jj);
}

template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
inline pState EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
WnSW(const int &nb, const int &ii, const int &jj) {
  return Wn(nb,ii,jj);
}

/**********************************************************************
 * EmbeddedBoundaries2D::Un?? -- Node conserved solution state.       *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
inline cState EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
UnNW(const int &nb, const int &ii, const int &jj) {
  return Un(nb,ii,jj+1);
}

template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
inline cState EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
UnNE(const int &nb, const int &ii, const int &jj) {
  return Un(nb,ii+1,jj+1);
}

template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
inline cState EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
UnSE(const int &nb, const int &ii, const int &jj) {
  return Un(nb,ii+1,jj);
}

template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
inline cState EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
UnSW(const int &nb, const int &ii, const int &jj) {
  return Un(nb,ii,jj);
}

/**********************************************************************
 **********************************************************************
 ** Adaptive mesh refinement routines.                                *
 **********************************************************************
 **********************************************************************/

/**********************************************************************
 * EmbeddedBoundaries2D::Prolong_Quadrilateral_Block --               *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
inline int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Prolong_Quadrilateral_Block(const int &Perform_Mesh_Adjustment,
			    const int &nb,
			    Quad_Soln_Block &SolnBlk_Original,
			    Grid2D_Quad_Block &AGrid_Original,
			    Adjusted_Mesh_Quad_Block &AMesh_Original,
			    const int &Sector) {

  int error_flag;
  int i_min, i_max, j_min, j_max, mesh_refinement_permitted;
  int i_fine, j_fine, i_coarse_min, j_coarse_min, coarse_cell_found;
  double distance, area_total_fine;
  Vector2D dX;

  // Allocate (re-allocate) memory for the solution of the refined
  // quadrilateral solution block as necessary.
  if ((SolnBlk_Original.NCi-2*SolnBlk_Original.Nghost-2*((SolnBlk_Original.NCi-2*SolnBlk_Original.Nghost)/2) != 0) || 
      (SolnBlk_Original.NCj-2*SolnBlk_Original.Nghost-2*((SolnBlk_Original.NCj-2*SolnBlk_Original.Nghost)/2) != 0) ||
      (SolnBlk_Original.NCi-2*SolnBlk_Original.Nghost < 2*SolnBlk_Original.Nghost) ||
      (SolnBlk_Original.NCj-2*SolnBlk_Original.Nghost < 2*SolnBlk_Original.Nghost) ||
      (SolnBlk_Original.Grid.Node == NULL) ||
      (SolnBlk_Original.U == NULL)) {
    mesh_refinement_permitted = 0;
  } else {
    mesh_refinement_permitted = 1;
    Copy_Solution_Block(Local_SolnBlk[nb],SolnBlk_Original);
    Refine_Mesh(Local_SolnBlk[nb].Grid,SolnBlk_Original.Grid,Sector);
  }

  if (mesh_refinement_permitted) {

    // Conduct mesh adjustment if necessary.
    if (Interface_Component_List.Ni && Perform_Mesh_Adjustment) {
      Store_Unadjusted_Mesh(nb);
      error_flag = Pre_Mesh_Adjustment_Painting(nb);
      if (!error_flag) Adjustment_Data[nb].Initialize();
      if (!error_flag) error_flag = Mesh_Adjustment_Sharp_Corners(nb);
      if (!error_flag) error_flag = Mesh_Adjustment_First(nb);
      if (!error_flag) error_flag = Mesh_Adjustment_Second(nb);
      if (!error_flag) error_flag = Mesh_Adjustment_Third(nb);
      if (!error_flag) error_flag = Mesh_Adjustment_Modify_Boundary_Splines(nb);
      if (error_flag) return error_flag;
    } else {
      Store_Unadjusted_Mesh(nb);
    }

    // Prolong the solution information from original solution block
    // to the refined solution block.
    switch(Sector) {
    case GRID2D_QUAD_BLOCK_SECTOR_NW :
      i_min = SolnBlk_Original.ICl;
      i_max = SolnBlk_Original.ICl+(SolnBlk_Original.ICu-SolnBlk_Original.ICl-1)/2;
      j_min = SolnBlk_Original.JCl+(SolnBlk_Original.JCu-SolnBlk_Original.JCl+1)/2; 
      j_max = SolnBlk_Original.JCu;
      break;
    case GRID2D_QUAD_BLOCK_SECTOR_NE :
      i_min = SolnBlk_Original.ICl+(SolnBlk_Original.ICu-SolnBlk_Original.ICl+1)/2;
      i_max = SolnBlk_Original.ICu;
      j_min = SolnBlk_Original.JCl+(SolnBlk_Original.JCu-SolnBlk_Original.JCl+1)/2; 
      j_max = SolnBlk_Original.JCu;
      break;
    case GRID2D_QUAD_BLOCK_SECTOR_SE :
      i_min = SolnBlk_Original.ICl+(SolnBlk_Original.ICu-SolnBlk_Original.ICl+1)/2;
      i_max = SolnBlk_Original.ICu;
      j_min = SolnBlk_Original.JCl; 
      j_max = SolnBlk_Original.JCl+(SolnBlk_Original.JCu-SolnBlk_Original.JCl-1)/2;
      break;
    case GRID2D_QUAD_BLOCK_SECTOR_SW :
      i_min = SolnBlk_Original.ICl;
      i_max = SolnBlk_Original.ICl+(SolnBlk_Original.ICu-SolnBlk_Original.ICl-1)/2;
      j_min = SolnBlk_Original.JCl; 
      j_max = SolnBlk_Original.JCl+(SolnBlk_Original.JCu-SolnBlk_Original.JCl-1)/2;
      break;
    default:
      i_min = SolnBlk_Original.ICl;
      i_max = SolnBlk_Original.ICl+(SolnBlk_Original.ICu-SolnBlk_Original.ICl-1)/2;
      j_min = SolnBlk_Original.JCl+(SolnBlk_Original.JCu-SolnBlk_Original.JCl+1)/2; 
      j_max = SolnBlk_Original.JCu;
      break;
    };

    if (!Interface_Component_List.Ni || !Perform_Mesh_Adjustment) {
      // No embedded boundaries exist.  Prolong the solution information
      // using pure injection (area-weighted average commented out).
      for (int j = j_min; j <= j_max; j++) {
	for (int i = i_min ; i <= i_max; i++) {
	  Local_SolnBlk[nb].U[2*(i-i_min)+Local_SolnBlk[nb].ICl  ][2*(j-j_min)+Local_SolnBlk[nb].JCl  ] = SolnBlk_Original.U[i][j];
	  Local_SolnBlk[nb].W[2*(i-i_min)+Local_SolnBlk[nb].ICl  ][2*(j-j_min)+Local_SolnBlk[nb].JCl  ] = W(SolnBlk_Original.U[i][j]);
	  Local_SolnBlk[nb].U[2*(i-i_min)+Local_SolnBlk[nb].ICl+1][2*(j-j_min)+Local_SolnBlk[nb].JCl  ] = SolnBlk_Original.U[i][j];
	  Local_SolnBlk[nb].W[2*(i-i_min)+Local_SolnBlk[nb].ICl+1][2*(j-j_min)+Local_SolnBlk[nb].JCl  ] = W(SolnBlk_Original.U[i][j]);
	  Local_SolnBlk[nb].U[2*(i-i_min)+Local_SolnBlk[nb].ICl  ][2*(j-j_min)+Local_SolnBlk[nb].JCl+1] = SolnBlk_Original.U[i][j];
	  Local_SolnBlk[nb].W[2*(i-i_min)+Local_SolnBlk[nb].ICl  ][2*(j-j_min)+Local_SolnBlk[nb].JCl+1] = W(SolnBlk_Original.U[i][j]);
	  Local_SolnBlk[nb].U[2*(i-i_min)+Local_SolnBlk[nb].ICl+1][2*(j-j_min)+Local_SolnBlk[nb].JCl+1] = SolnBlk_Original.U[i][j];
	  Local_SolnBlk[nb].W[2*(i-i_min)+Local_SolnBlk[nb].ICl+1][2*(j-j_min)+Local_SolnBlk[nb].JCl+1] = W(SolnBlk_Original.U[i][j]);
	}
      }

    } else {
      // Embedded boundaries exist.
      for (int j = j_min; j <= j_max; j++) {
	for (int i = i_min ; i <= i_max; i++) {
	  if (AMesh_Original.cell_status[i-1][j-1] != CELL_STATUS_ACTIVE ||
	      AMesh_Original.cell_status[i  ][j-1] != CELL_STATUS_ACTIVE ||
	      AMesh_Original.cell_status[i+1][j-1] != CELL_STATUS_ACTIVE ||
	      AMesh_Original.cell_status[i-1][j  ] != CELL_STATUS_ACTIVE ||
	      AMesh_Original.cell_status[i  ][j  ] != CELL_STATUS_ACTIVE ||
	      AMesh_Original.cell_status[i+1][j  ] != CELL_STATUS_ACTIVE ||
	      AMesh_Original.cell_status[i-1][j+1] != CELL_STATUS_ACTIVE ||
	      AMesh_Original.cell_status[i  ][j+1] != CELL_STATUS_ACTIVE ||
	      AMesh_Original.cell_status[i+1][j+1] != CELL_STATUS_ACTIVE) {
	    // The fine cell is near an embedded boundary.  Prolong the
	    // solution information by injection from the nearest active
	    // coarse cell.  If the fine cell is inactive then set to
	    // standard atmosphere.
	    i_fine = 2*(i-i_min)+Local_SolnBlk[nb].ICl; // SW fine cell i index.
	    j_fine = 2*(j-j_min)+Local_SolnBlk[nb].JCl; // SW fine cell j index.
	    for (int jj = j_fine; jj <= j_fine+1; jj++) {
	      for (int ii = i_fine; ii <= i_fine+1; ii++) {
		if (Mesh[nb].cell_status[ii][jj] != CELL_STATUS_ACTIVE) {
		  Local_SolnBlk[nb].U[ii][jj].Standard_Atmosphere();
		  Local_SolnBlk[nb].W[ii][jj].Standard_Atmosphere();
		} else if (Mesh[nb].cell_status[ii][jj] == CELL_STATUS_ACTIVE) {
		  distance = MILLION;
		  coarse_cell_found = OFF;
		  for (int jc = j-1; jc <= j+1; jc++) {
		    for (int ic = i-1; ic <= i+1; ic++) {
		      if (AMesh_Original.cell_status[ic][jc] == CELL_STATUS_ACTIVE &&
			  abs(AGrid_Original.Cell[ic][jc].Xc -
			      Local_SolnBlk[nb].Grid.Cell[ii][jj].Xc) < distance) {
			i_coarse_min = ic; j_coarse_min = jc;
			distance = abs(AGrid_Original.Cell[ic][jc].Xc -
				       Local_SolnBlk[nb].Grid.Cell[ii][jj].Xc);
			coarse_cell_found = ON;
		      }
		    }
		  }
		  if (coarse_cell_found) {
		    Local_SolnBlk[nb].U[ii][jj] = SolnBlk_Original.U[i_coarse_min][j_coarse_min];
		    Local_SolnBlk[nb].W[ii][jj] = W(Local_SolnBlk[nb].U[ii][jj]);
		  } else {
		    return 7101;
		  }
		} else {
		  return 7102;
		}
	      }
	    }

	  } else {
	    // Prolong by injection from the parent cell.
	    Local_SolnBlk[nb].U[2*(i-i_min)+Local_SolnBlk[nb].ICl  ][2*(j-j_min)+Local_SolnBlk[nb].JCl  ] = SolnBlk_Original.U[i][j];
	    Local_SolnBlk[nb].W[2*(i-i_min)+Local_SolnBlk[nb].ICl  ][2*(j-j_min)+Local_SolnBlk[nb].JCl  ] = W(SolnBlk_Original.U[i][j]);
	    Local_SolnBlk[nb].U[2*(i-i_min)+Local_SolnBlk[nb].ICl+1][2*(j-j_min)+Local_SolnBlk[nb].JCl  ] = SolnBlk_Original.U[i][j];
	    Local_SolnBlk[nb].W[2*(i-i_min)+Local_SolnBlk[nb].ICl+1][2*(j-j_min)+Local_SolnBlk[nb].JCl  ] = W(SolnBlk_Original.U[i][j]);
	    Local_SolnBlk[nb].U[2*(i-i_min)+Local_SolnBlk[nb].ICl  ][2*(j-j_min)+Local_SolnBlk[nb].JCl+1] = SolnBlk_Original.U[i][j];
	    Local_SolnBlk[nb].W[2*(i-i_min)+Local_SolnBlk[nb].ICl  ][2*(j-j_min)+Local_SolnBlk[nb].JCl+1] = W(SolnBlk_Original.U[i][j]);
	    Local_SolnBlk[nb].U[2*(i-i_min)+Local_SolnBlk[nb].ICl+1][2*(j-j_min)+Local_SolnBlk[nb].JCl+1] = SolnBlk_Original.U[i][j];
	    Local_SolnBlk[nb].W[2*(i-i_min)+Local_SolnBlk[nb].ICl+1][2*(j-j_min)+Local_SolnBlk[nb].JCl+1] = W(SolnBlk_Original.U[i][j]);

 	  }

	}
      }

    }

    // Prolong the east and west boundary states.
    for (int j = j_min-SolnBlk_Original.Nghost/2; j <= j_max+SolnBlk_Original.Nghost/2; j++) {
      Local_SolnBlk[nb].WoW[2*(j-j_min)+Local_SolnBlk[nb].JCl  ] = SolnBlk_Original.WoW[j];
      Local_SolnBlk[nb].WoW[2*(j-j_min)+Local_SolnBlk[nb].JCl+1] = SolnBlk_Original.WoW[j];
      Local_SolnBlk[nb].WoE[2*(j-j_min)+Local_SolnBlk[nb].JCl  ] = SolnBlk_Original.WoE[j];
      Local_SolnBlk[nb].WoE[2*(j-j_min)+Local_SolnBlk[nb].JCl+1] = SolnBlk_Original.WoE[j];
    }

    // Prolong the north and south boundary states.
    for (int i = i_min-SolnBlk_Original.Nghost/2; i <= i_max+SolnBlk_Original.Nghost/2; i++) {
      Local_SolnBlk[nb].WoS[2*(i-i_min)+Local_SolnBlk[nb].ICl  ] = SolnBlk_Original.WoS[i];
      Local_SolnBlk[nb].WoS[2*(i-i_min)+Local_SolnBlk[nb].ICl+1] = SolnBlk_Original.WoS[i];
      Local_SolnBlk[nb].WoN[2*(i-i_min)+Local_SolnBlk[nb].ICl  ] = SolnBlk_Original.WoN[i];
      Local_SolnBlk[nb].WoN[2*(i-i_min)+Local_SolnBlk[nb].ICl+1] = SolnBlk_Original.WoN[i];
    }

  }

  // Prolongation of solution block was successful.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Restrict_Quadrilateral_Block --              *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
inline int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Restrict_Quadrilateral_Block(const int &Perform_Mesh_Adjustment,
			     const int &nb,
			     Quad_Soln_Block &SolnBlk_Original_SW,
			     Quad_Soln_Block &SolnBlk_Original_SE,
			     Quad_Soln_Block &SolnBlk_Original_NW,
			     Quad_Soln_Block &SolnBlk_Original_NE,
			     Grid2D_Quad_Block &AGrid_Original_SW,
			     Grid2D_Quad_Block &AGrid_Original_SE,
			     Grid2D_Quad_Block &AGrid_Original_NW,
			     Grid2D_Quad_Block &AGrid_Original_NE,
			     Adjusted_Mesh_Quad_Block &AMesh_Original_SW,
			     Adjusted_Mesh_Quad_Block &AMesh_Original_SE,
			     Adjusted_Mesh_Quad_Block &AMesh_Original_NW,
			     Adjusted_Mesh_Quad_Block &AMesh_Original_NE) {

  int error_flag, i_coarse, j_coarse, mesh_coarsening_permitted;
  double A, total_area;
  Polygon Pc, Pf;

  // Allocate memory for the cells and nodes for the coarsened 
  // quadrilateral mesh block.

  if ((SolnBlk_Original_SW.NCi-2*SolnBlk_Original_SW.Nghost-
       2*((SolnBlk_Original_SW.NCi-2*SolnBlk_Original_SW.Nghost)/2) != 0) || 
      (SolnBlk_Original_SW.NCj-2*SolnBlk_Original_SW.Nghost-
       2*((SolnBlk_Original_SW.NCj-2*SolnBlk_Original_SW.Nghost)/2) != 0) ||
      (SolnBlk_Original_SW.NCi-2*SolnBlk_Original_SW.Nghost < 4) ||
      (SolnBlk_Original_SW.NCj-2*SolnBlk_Original_SW.Nghost < 4) ||
      (SolnBlk_Original_SE.NCi != SolnBlk_Original_SW.NCi) ||
      (SolnBlk_Original_SE.NCj != SolnBlk_Original_SW.NCj) ||
      (SolnBlk_Original_NW.NCi != SolnBlk_Original_SW.NCi) ||
      (SolnBlk_Original_NW.NCj != SolnBlk_Original_SW.NCj) ||
      (SolnBlk_Original_NE.NCi != SolnBlk_Original_SW.NCi) ||
      (SolnBlk_Original_NE.NCj != SolnBlk_Original_SW.NCj) ||
      (SolnBlk_Original_SW.Grid.Node == NULL) ||
      (SolnBlk_Original_SE.Grid.Node == NULL) ||
      (SolnBlk_Original_NW.Grid.Node == NULL) ||
      (SolnBlk_Original_NE.Grid.Node == NULL) ||
      (SolnBlk_Original_SW.U == NULL) ||
      (SolnBlk_Original_SE.U == NULL) ||
      (SolnBlk_Original_NW.U == NULL) ||
      (SolnBlk_Original_NE.U == NULL)) {
    mesh_coarsening_permitted = 0;
  } else {
    mesh_coarsening_permitted = 1;
    Copy_Solution_Block(Local_SolnBlk[nb],SolnBlk_Original_SW);
    Coarsen_Mesh(Local_SolnBlk[nb].Grid,
		 SolnBlk_Original_SW.Grid,
		 SolnBlk_Original_SE.Grid,
		 SolnBlk_Original_NW.Grid,
		 SolnBlk_Original_NE.Grid);
  }

  if (mesh_coarsening_permitted) {

    // Conduct mesh adjustment if necessary.
    if (Interface_Component_List.Ni && Perform_Mesh_Adjustment) {
      Store_Unadjusted_Mesh(nb);
      error_flag = Pre_Mesh_Adjustment_Painting(nb);
      if (!error_flag) Adjustment_Data[nb].Initialize();
      if (!error_flag) error_flag = Mesh_Adjustment_Sharp_Corners(nb);
      if (!error_flag) error_flag = Mesh_Adjustment_First(nb);
      if (!error_flag) error_flag = Mesh_Adjustment_Second(nb);
      if (!error_flag) error_flag = Mesh_Adjustment_Third(nb);
      if (!error_flag) error_flag = Mesh_Adjustment_Modify_Boundary_Splines(nb);
      if (error_flag) return error_flag;
    } else {
      Store_Unadjusted_Mesh(nb);
    }

    // Restrict the solution information from the four original
    // solution blocks to the newly coarsened solution block.

    // South-West corner fine block:	
    for (int j_fine = SolnBlk_Original_SW.JCl; j_fine <= SolnBlk_Original_SW.JCu; j_fine += 2) {
      for (int i_fine = SolnBlk_Original_SW.ICl; i_fine <= SolnBlk_Original_SW.ICu; i_fine += 2) {

	i_coarse = (i_fine-SolnBlk_Original_SW.ICl)/2 + Local_SolnBlk[nb].ICl;
	j_coarse = (j_fine-SolnBlk_Original_SW.JCl)/2 + Local_SolnBlk[nb].JCl;

	if (Mesh[nb].cell_status[i_coarse][j_coarse] != CELL_STATUS_ACTIVE) {
	  // The cell is inactive.  Set to standard atmosphere.
	  Local_SolnBlk[nb].U[i_coarse][j_coarse].Standard_Atmosphere();
	  Local_SolnBlk[nb].W[i_coarse][j_coarse].Standard_Atmosphere();

	} else if (!Interface_Union_List.Ni || Perform_Mesh_Adjustment) {
	  // No embedded boundaries exist.  Restrict the solution
	  // information from the fine cells to the coarse cell using an
	  // area-weighted average of the fine cell solution information.
	  Local_SolnBlk[nb].U[i_coarse][j_coarse] = (SolnBlk_Original_SW.Grid.Cell[i_fine  ][j_fine  ].A*
						     SolnBlk_Original_SW.U[i_fine  ][j_fine  ] +
						     SolnBlk_Original_SW.Grid.Cell[i_fine+1][j_fine  ].A*
						     SolnBlk_Original_SW.U[i_fine+1][j_fine  ] + 
						     SolnBlk_Original_SW.Grid.Cell[i_fine  ][j_fine+1].A*
						     SolnBlk_Original_SW.U[i_fine  ][j_fine+1] +
						     SolnBlk_Original_SW.Grid.Cell[i_fine+1][j_fine+1].A*
						     SolnBlk_Original_SW.U[i_fine+1][j_fine+1])/
                                                    (SolnBlk_Original_SW.Grid.Cell[i_fine  ][j_fine  ].A +
						     SolnBlk_Original_SW.Grid.Cell[i_fine+1][j_fine  ].A +
						     SolnBlk_Original_SW.Grid.Cell[i_fine  ][j_fine+1].A +
						     SolnBlk_Original_SW.Grid.Cell[i_fine+1][j_fine+1].A);
                                                     //Local_SolnBlk[nb].Grid.Cell[i_coarse][j_coarse].A;
	  Local_SolnBlk[nb].W[i_coarse][j_coarse] = W(Local_SolnBlk[nb].U[i_coarse][j_coarse]);

	} else {

	  if (Mesh[nb].cell_type[i_coarse][j_coarse] == CELL_TYPE_QUADRILATERAL &&
	      AMesh_Original_SW.cell_type[i_fine  ][j_fine  ] == CELL_TYPE_QUADRILATERAL &&
	      AMesh_Original_SW.cell_type[i_fine+1][j_fine  ] == CELL_TYPE_QUADRILATERAL &&
	      AMesh_Original_SW.cell_type[i_fine  ][j_fine+1] == CELL_TYPE_QUADRILATERAL &&
	      AMesh_Original_SW.cell_type[i_fine+1][j_fine+1] == CELL_TYPE_QUADRILATERAL &&
	      AMesh_Original_SW.cell_status[i_fine  ][j_fine  ] == CELL_STATUS_ACTIVE &&
	      AMesh_Original_SW.cell_status[i_fine+1][j_fine  ] == CELL_STATUS_ACTIVE &&
	      AMesh_Original_SW.cell_status[i_fine  ][j_fine+1] == CELL_STATUS_ACTIVE &&
	      AMesh_Original_SW.cell_status[i_fine+1][j_fine+1] == CELL_STATUS_ACTIVE) {
	    // None of the associated cells have been adjusted.
	    // Restrict the solution information using an area-weighted
	    // average of the cells from the stored adjusted fine mesh
	    // and the coarse mesh.
	    Local_SolnBlk[nb].U[i_coarse][j_coarse] = (AGrid_Original_SW.Cell[i_fine  ][j_fine  ].A*
						       SolnBlk_Original_SW.U[i_fine  ][j_fine  ] +
						       AGrid_Original_SW.Cell[i_fine+1][j_fine  ].A*
						       SolnBlk_Original_SW.U[i_fine+1][j_fine  ] + 
						       AGrid_Original_SW.Cell[i_fine  ][j_fine+1].A*
						       SolnBlk_Original_SW.U[i_fine  ][j_fine+1] +
						       AGrid_Original_SW.Cell[i_fine+1][j_fine+1].A*
						       SolnBlk_Original_SW.U[i_fine+1][j_fine+1])/
             	                                      (AGrid_Original_SW.Cell[i_fine  ][j_fine  ].A +
						       AGrid_Original_SW.Cell[i_fine+1][j_fine  ].A +
						       AGrid_Original_SW.Cell[i_fine  ][j_fine+1].A +
						       AGrid_Original_SW.Cell[i_fine+1][j_fine+1].A);
	    Local_SolnBlk[nb].W[i_coarse][j_coarse] = W(Local_SolnBlk[nb].U[i_coarse][j_coarse]);

	  } else {
	    // The coarse and/or fine cells have been adjusted.
	    // Restrict the solution information using an area-weight
	    // average based on the area of intersection of the coarse
	    // grid and fine grid cells.
	    total_area = ZERO;
	    Local_SolnBlk[nb].U[i_coarse][j_coarse].Vacuum();
	    // Create the coarse cell polygon.
	    Pc.convert(Local_SolnBlk[nb].Grid.nodeSW(i_coarse,j_coarse).X,Local_SolnBlk[nb].Grid.nodeSE(i_coarse,j_coarse).X,
		       Local_SolnBlk[nb].Grid.nodeNE(i_coarse,j_coarse).X,Local_SolnBlk[nb].Grid.nodeNW(i_coarse,j_coarse).X);
	    // Search the fine cells for intersecting cells.
	    for (int j = j_fine-2; j < j_fine+4; j++) {
	      for (int i = i_fine-2; i < i_fine+4; i++) {
		if (AMesh_Original_SW.cell_status[i][j] == CELL_STATUS_ACTIVE) {
		  // Create the fine cell polygon.
		  Pf.convert(AGrid_Original_SW.nodeSW(i,j).X,AGrid_Original_SW.nodeSE(i,j).X,
			     AGrid_Original_SW.nodeNE(i,j).X,AGrid_Original_SW.nodeNW(i,j).X);
		  // Determine the area of intersection between the fine
		  // and coarse cell polygons.
		  A = Polygon_Intersection_Area(Pc,Pf);
		  Pf.deallocate();
		  // Add the area of intersection to the total area of
		  // intersection.
		  total_area += A;
		  // Restrict fine cell solution to the coarse cell
		  // based on the area of interesection.
		  Local_SolnBlk[nb].U[i_coarse][j_coarse] += SolnBlk_Original_SW.U[i][j]*A;
		}
	      }
	    }
	    Pc.deallocate();
	    // Complete solution restriction by dividing the restricted
	    // solution by the total area of intersection.
	    Local_SolnBlk[nb].U[i_coarse][j_coarse] = Local_SolnBlk[nb].U[i_coarse][j_coarse]/total_area;
	    Local_SolnBlk[nb].W[i_coarse][j_coarse] = W(Local_SolnBlk[nb].U[i_coarse][j_coarse]);

	  }

	}

      }
    }

    for (int j_fine = SolnBlk_Original_SW.JCl-SolnBlk_Original_SW.Nghost; 
	 j_fine <= SolnBlk_Original_SW.JCu+SolnBlk_Original_SW.Nghost; j_fine += 2) {
      j_coarse = (j_fine-SolnBlk_Original_SW.JCl)/2+Local_SolnBlk[nb].JCl;
      Local_SolnBlk[nb].WoW[j_coarse] = SolnBlk_Original_SW.WoW[j_fine];
      if (j_fine == SolnBlk_Original_SW.JCl-SolnBlk_Original_SW.Nghost) {
	Local_SolnBlk[nb].WoW[j_coarse-1] = SolnBlk_Original_SW.WoW[j_fine];
      }
    }

    for (int i_fine = SolnBlk_Original_SW.ICl-SolnBlk_Original_SW.Nghost; 
	 i_fine <= SolnBlk_Original_SW.ICu+SolnBlk_Original_SW.Nghost; i_fine += 2) {
      i_coarse = (i_fine-SolnBlk_Original_SW.ICl)/2+Local_SolnBlk[nb].ICl;
      Local_SolnBlk[nb].WoS[i_coarse] = SolnBlk_Original_SW.WoS[i_fine];
      if (i_fine == SolnBlk_Original_SW.ICl-SolnBlk_Original_SW.Nghost) {
	Local_SolnBlk[nb].WoS[i_coarse-1] = SolnBlk_Original_SW.WoS[i_fine];
      }
    }

    // South-East corner fine block:
    for (int j_fine = SolnBlk_Original_SE.JCl; j_fine <= SolnBlk_Original_SE.JCu; j_fine += 2) {
      for (int i_fine = SolnBlk_Original_SE.ICl; i_fine <= SolnBlk_Original_SE.ICu; i_fine += 2) {

	i_coarse = (i_fine-SolnBlk_Original_SE.ICl)/2 + (Local_SolnBlk[nb].ICu-Local_SolnBlk[nb].ICl+1)/2+Local_SolnBlk[nb].ICl;
	j_coarse = (j_fine-SolnBlk_Original_SE.JCl)/2 + Local_SolnBlk[nb].JCl;

	if (Mesh[nb].cell_status[i_coarse][j_coarse] != CELL_STATUS_ACTIVE) {

	  // The cell is inactive.  Set to standard atmosphere.
	  Local_SolnBlk[nb].U[i_coarse][j_coarse].Standard_Atmosphere();
	  Local_SolnBlk[nb].W[i_coarse][j_coarse].Standard_Atmosphere();

	} else if (!Interface_Union_List.Ni) {
	  // No embedded boundaries exist.  Restrict the solution
	  // information from the fine cells to the coarse cell using an
	  // area-weighted average of the fine cell solution information.
	  Local_SolnBlk[nb].U[i_coarse][j_coarse] = (SolnBlk_Original_SE.Grid.Cell[i_fine  ][j_fine  ].A*
						     SolnBlk_Original_SE.U[i_fine  ][j_fine  ] +
						     SolnBlk_Original_SE.Grid.Cell[i_fine+1][j_fine  ].A*
						     SolnBlk_Original_SE.U[i_fine+1][j_fine  ] + 
						     SolnBlk_Original_SE.Grid.Cell[i_fine  ][j_fine+1].A*
						     SolnBlk_Original_SE.U[i_fine  ][j_fine+1] +
						     SolnBlk_Original_SE.Grid.Cell[i_fine+1][j_fine+1].A*
						     SolnBlk_Original_SE.U[i_fine+1][j_fine+1])/
                                                    (SolnBlk_Original_SE.Grid.Cell[i_fine  ][j_fine  ].A +
						     SolnBlk_Original_SE.Grid.Cell[i_fine+1][j_fine  ].A +
						     SolnBlk_Original_SE.Grid.Cell[i_fine  ][j_fine+1].A +
						     SolnBlk_Original_SE.Grid.Cell[i_fine+1][j_fine+1].A);
                                                     //Local_SolnBlk[nb].Grid.Cell[i_coarse][j_coarse].A;
	  Local_SolnBlk[nb].W[i_coarse][j_coarse] = W(Local_SolnBlk[nb].U[i_coarse][j_coarse]);

	} else {

	  if (Mesh[nb].cell_type[i_coarse][j_coarse] == CELL_TYPE_QUADRILATERAL &&
	      AMesh_Original_SE.cell_type[i_fine  ][j_fine  ] == CELL_TYPE_QUADRILATERAL &&
	      AMesh_Original_SE.cell_type[i_fine+1][j_fine  ] == CELL_TYPE_QUADRILATERAL &&
	      AMesh_Original_SE.cell_type[i_fine  ][j_fine+1] == CELL_TYPE_QUADRILATERAL &&
	      AMesh_Original_SE.cell_type[i_fine+1][j_fine+1] == CELL_TYPE_QUADRILATERAL &&
	      AMesh_Original_SE.cell_status[i_fine  ][j_fine  ] == CELL_STATUS_ACTIVE &&
	      AMesh_Original_SE.cell_status[i_fine+1][j_fine  ] == CELL_STATUS_ACTIVE &&
	      AMesh_Original_SE.cell_status[i_fine  ][j_fine+1] == CELL_STATUS_ACTIVE &&
	      AMesh_Original_SE.cell_status[i_fine+1][j_fine+1] == CELL_STATUS_ACTIVE) {
	    // None of the associated cells have been adjusted.
	    // Restrict the solution information using an area-weighted
	    // average of the cells from the stored adjusted fine mesh
	    // and the coarse mesh.
	    Local_SolnBlk[nb].U[i_coarse][j_coarse] = (AGrid_Original_SE.Cell[i_fine  ][j_fine  ].A*
						       SolnBlk_Original_SE.U[i_fine  ][j_fine  ] +
						       AGrid_Original_SE.Cell[i_fine+1][j_fine  ].A*
						       SolnBlk_Original_SE.U[i_fine+1][j_fine  ] + 
						       AGrid_Original_SE.Cell[i_fine  ][j_fine+1].A*
						       SolnBlk_Original_SE.U[i_fine  ][j_fine+1] +
						       AGrid_Original_SE.Cell[i_fine+1][j_fine+1].A*
						       SolnBlk_Original_SE.U[i_fine+1][j_fine+1])/
	                                              (AGrid_Original_SE.Cell[i_fine  ][j_fine  ].A +
						       AGrid_Original_SE.Cell[i_fine+1][j_fine  ].A +
						       AGrid_Original_SE.Cell[i_fine  ][j_fine+1].A +
						       AGrid_Original_SE.Cell[i_fine+1][j_fine+1].A);
	    Local_SolnBlk[nb].W[i_coarse][j_coarse] = W(Local_SolnBlk[nb].U[i_coarse][j_coarse]);

	  } else {
	    // The coarse and/or fine cells have been adjusted.
	    // Restrict the solution information using an area-weight
	    // average based on the area of intersection of the coarse
	    // grid and fine grid cells.
	    total_area = ZERO;
	    Local_SolnBlk[nb].U[i_coarse][j_coarse].Vacuum();
	    // Create the coarse cell polygon.
	    Pc.convert(Local_SolnBlk[nb].Grid.nodeSW(i_coarse,j_coarse).X,Local_SolnBlk[nb].Grid.nodeSE(i_coarse,j_coarse).X,
		       Local_SolnBlk[nb].Grid.nodeNE(i_coarse,j_coarse).X,Local_SolnBlk[nb].Grid.nodeNW(i_coarse,j_coarse).X);
	    // Search the fine cells for intersecting cells.
	    for (int j = j_fine-2; j < j_fine+4; j++) {
	      for (int i = i_fine-2; i < i_fine+4; i++) {
		if (AMesh_Original_SE.cell_status[i][j] == CELL_STATUS_ACTIVE) {
		  // Create the fine cell polygon.
		  Pf.convert(SolnBlk_Original_SE.Grid.nodeSW(i,j).X,SolnBlk_Original_SE.Grid.nodeSE(i,j).X,
			     SolnBlk_Original_SE.Grid.nodeNE(i,j).X,SolnBlk_Original_SE.Grid.nodeNW(i,j).X);
		  // Determine the area of intersection between the fine
		  // and coarse cell polygons.
		  A = Polygon_Intersection_Area(Pc,Pf);
		  Pf.deallocate();
		  // Add the area of intersection to the total area of
		  // intersection.
		  total_area += A;
		  // Restrict fine cell solution to the coarse cell
		  // based on the area of interesection.
		  Local_SolnBlk[nb].U[i_coarse][j_coarse] += SolnBlk_Original_SE.U[i][j]*A;
		}
	      }
	    }
	    Pc.deallocate();
	    // Complete solution restriction by dividing the restricted
	    // solution by the total area of intersection.
	    Local_SolnBlk[nb].U[i_coarse][j_coarse] = Local_SolnBlk[nb].U[i_coarse][j_coarse]/total_area;
	    Local_SolnBlk[nb].W[i_coarse][j_coarse] = W(Local_SolnBlk[nb].U[i_coarse][j_coarse]);

	  }

 	}

      }
    }

    for (int j_fine = SolnBlk_Original_SE.JCl-SolnBlk_Original_SE.Nghost; 
	 j_fine <= SolnBlk_Original_SE.JCu+SolnBlk_Original_SE.Nghost; j_fine += 2) {
      j_coarse = (j_fine-SolnBlk_Original_SE.JCl)/2+Local_SolnBlk[nb].JCl;
      Local_SolnBlk[nb].WoE[j_coarse] = SolnBlk_Original_SE.WoE[j_fine];
      if (j_fine == SolnBlk_Original_SE.JCl-SolnBlk_Original_SE.Nghost) {
	Local_SolnBlk[nb].WoE[j_coarse-1] = SolnBlk_Original_SE.WoE[j_fine];
      }
    }

    for (int i_fine = SolnBlk_Original_SE.ICl-SolnBlk_Original_SE.Nghost; 
	 i_fine <= SolnBlk_Original_SE.ICu+SolnBlk_Original_SE.Nghost; i_fine += 2) {
      i_coarse = (i_fine-SolnBlk_Original_SE.ICl)/2+(Local_SolnBlk[nb].ICu-Local_SolnBlk[nb].ICl+1)/2+Local_SolnBlk[nb].ICl;
      Local_SolnBlk[nb].WoS[i_coarse] = SolnBlk_Original_SE.WoS[i_fine];
      if (i_fine == SolnBlk_Original_SE.ICu+SolnBlk_Original_SE.Nghost) {
	Local_SolnBlk[nb].WoS[i_coarse+1] = SolnBlk_Original_SE.WoS[i_fine];
      }
    }

    // North-West corner fine block:
    for (int j_fine = SolnBlk_Original_NW.JCl; j_fine <= SolnBlk_Original_NW.JCu; j_fine += 2) {
      for (int i_fine = SolnBlk_Original_NW.ICl; i_fine <= SolnBlk_Original_NW.ICu; i_fine += 2) {

	i_coarse = (i_fine-SolnBlk_Original_NW.ICl)/2 + Local_SolnBlk[nb].ICl;
	j_coarse = (j_fine-SolnBlk_Original_NW.JCl)/2 + (Local_SolnBlk[nb].JCu-Local_SolnBlk[nb].JCl+1)/2+Local_SolnBlk[nb].JCl;

	if (Mesh[nb].cell_status[i_coarse][j_coarse] != CELL_STATUS_ACTIVE) {

	  // The cell is inactive.  Set to standard atmosphere.
 	  Local_SolnBlk[nb].U[i_coarse][j_coarse].Standard_Atmosphere();
 	  Local_SolnBlk[nb].W[i_coarse][j_coarse].Standard_Atmosphere();

	} else if (!Interface_Union_List.Ni) {
	  // No embedded boundaries exist.  Restrict the solution
	  // information from the fine cells to the coarse cell using an
	  // area-weighted average of the fine cell solution information.
	  Local_SolnBlk[nb].U[i_coarse][j_coarse] = (SolnBlk_Original_NW.Grid.Cell[i_fine  ][j_fine  ].A*
						     SolnBlk_Original_NW.U[i_fine  ][j_fine  ] +
						     SolnBlk_Original_NW.Grid.Cell[i_fine+1][j_fine  ].A*
						     SolnBlk_Original_NW.U[i_fine+1][j_fine  ] + 
						     SolnBlk_Original_NW.Grid.Cell[i_fine  ][j_fine+1].A*
						     SolnBlk_Original_NW.U[i_fine  ][j_fine+1] +
						     SolnBlk_Original_NW.Grid.Cell[i_fine+1][j_fine+1].A*
						     SolnBlk_Original_NW.U[i_fine+1][j_fine+1])/
                                                    (SolnBlk_Original_NW.Grid.Cell[i_fine  ][j_fine  ].A +
						     SolnBlk_Original_NW.Grid.Cell[i_fine+1][j_fine  ].A +
						     SolnBlk_Original_NW.Grid.Cell[i_fine  ][j_fine+1].A +
						     SolnBlk_Original_NW.Grid.Cell[i_fine+1][j_fine+1].A);
                                                    //Local_SolnBlk[nb].Grid.Cell[i_coarse][j_coarse].A;
	  Local_SolnBlk[nb].W[i_coarse][j_coarse] = W(Local_SolnBlk[nb].U[i_coarse][j_coarse]);

	} else {

	  if (Mesh[nb].cell_type[i_coarse][j_coarse] == CELL_TYPE_QUADRILATERAL &&
	      AMesh_Original_NW.cell_type[i_fine  ][j_fine  ] == CELL_TYPE_QUADRILATERAL &&
	      AMesh_Original_NW.cell_type[i_fine+1][j_fine  ] == CELL_TYPE_QUADRILATERAL &&
	      AMesh_Original_NW.cell_type[i_fine  ][j_fine+1] == CELL_TYPE_QUADRILATERAL &&
	      AMesh_Original_NW.cell_type[i_fine+1][j_fine+1] == CELL_TYPE_QUADRILATERAL &&
	      AMesh_Original_NW.cell_status[i_fine  ][j_fine  ] == CELL_STATUS_ACTIVE &&
	      AMesh_Original_NW.cell_status[i_fine+1][j_fine  ] == CELL_STATUS_ACTIVE &&
	      AMesh_Original_NW.cell_status[i_fine  ][j_fine+1] == CELL_STATUS_ACTIVE &&
	      AMesh_Original_NW.cell_status[i_fine+1][j_fine+1] == CELL_STATUS_ACTIVE) {
	    // None of the associated cells have been adjusted.
	    // Restrict the solution information using an area-weighted
	    // average of the cells from the stored adjusted fine mesh
	    // and the coarse mesh.
	    Local_SolnBlk[nb].U[i_coarse][j_coarse] = (AGrid_Original_NW.Cell[i_fine  ][j_fine  ].A*
						       SolnBlk_Original_NW.U[i_fine  ][j_fine  ] +
						       AGrid_Original_NW.Cell[i_fine+1][j_fine  ].A*
						       SolnBlk_Original_NW.U[i_fine+1][j_fine  ] + 
						       AGrid_Original_NW.Cell[i_fine  ][j_fine+1].A*
						       SolnBlk_Original_NW.U[i_fine  ][j_fine+1] +
						       AGrid_Original_NW.Cell[i_fine+1][j_fine+1].A*
						       SolnBlk_Original_NW.U[i_fine+1][j_fine+1])/
                                                      (AGrid_Original_NW.Cell[i_fine  ][j_fine  ].A +
						       AGrid_Original_NW.Cell[i_fine+1][j_fine  ].A +
						       AGrid_Original_NW.Cell[i_fine  ][j_fine+1].A +
						       AGrid_Original_NW.Cell[i_fine+1][j_fine+1].A);
	    Local_SolnBlk[nb].W[i_coarse][j_coarse] = W(Local_SolnBlk[nb].U[i_coarse][j_coarse]);

 	  } else {
	    // The coarse and/or fine cells have been adjusted.
	    // Restrict the solution information using an area-weight
	    // average based on the area of intersection of the coarse
	    // grid and fine grid cells.
 	    total_area = ZERO;
 	    Local_SolnBlk[nb].U[i_coarse][j_coarse].Vacuum();
 	    // Create the coarse cell polygon.
 	    Pc.convert(Local_SolnBlk[nb].Grid.nodeSW(i_coarse,j_coarse).X,Local_SolnBlk[nb].Grid.nodeSE(i_coarse,j_coarse).X,
 		       Local_SolnBlk[nb].Grid.nodeNE(i_coarse,j_coarse).X,Local_SolnBlk[nb].Grid.nodeNW(i_coarse,j_coarse).X);
 	    // Search the fine cells for intersecting cells.
 	    for (int j = j_fine-2; j < j_fine+4; j++) {
 	      for (int i = i_fine-2; i < i_fine+4; i++) {
 		if (AMesh_Original_NW.cell_status[i][j] == CELL_STATUS_ACTIVE) {
 		  // Create the fine cell polygon.
 		  Pf.convert(SolnBlk_Original_NW.Grid.nodeSW(i,j).X,SolnBlk_Original_NW.Grid.nodeSE(i,j).X,
 			     SolnBlk_Original_NW.Grid.nodeNE(i,j).X,SolnBlk_Original_NW.Grid.nodeNW(i,j).X);
 		  // Determine the area of intersection between the fine
 		  // and coarse cell polygons.
 		  A = Polygon_Intersection_Area(Pc,Pf);
 		  Pf.deallocate();
 		  // Add the area of intersection to the total area of
 		  // intersection.
 		  total_area += A;
 		  // Restrict fine cell solution to the coarse cell
 		  // based on the area of interesection.
 		  Local_SolnBlk[nb].U[i_coarse][j_coarse] += SolnBlk_Original_SE.U[i][j]*A;
 		}
 	      }
 	    }
 	    Pc.deallocate();
 	    // Complete solution restriction by dividing the restricted
 	    // solution by the total area of intersection.
 	    Local_SolnBlk[nb].U[i_coarse][j_coarse] = Local_SolnBlk[nb].U[i_coarse][j_coarse]/total_area;
 	    Local_SolnBlk[nb].W[i_coarse][j_coarse] = W(Local_SolnBlk[nb].U[i_coarse][j_coarse]);

 	  }

	}

      }
    }

    for (int j_fine = SolnBlk_Original_NW.JCl-SolnBlk_Original_NW.Nghost; 
	 j_fine <= SolnBlk_Original_NW.JCu+SolnBlk_Original_NW.Nghost; j_fine += 2) {
      j_coarse = (j_fine-SolnBlk_Original_NW.JCl)/2 + (Local_SolnBlk[nb].JCu-Local_SolnBlk[nb].JCl+1)/2+Local_SolnBlk[nb].JCl;
      Local_SolnBlk[nb].WoW[j_coarse] = SolnBlk_Original_NW.WoW[j_fine];
      if (j_fine == SolnBlk_Original_NW.JCu+SolnBlk_Original_NW.Nghost) {
	Local_SolnBlk[nb].WoW[j_coarse+1] = SolnBlk_Original_NW.WoW[j_fine];
      }
    }

    for (int i_fine = SolnBlk_Original_NW.ICl-SolnBlk_Original_NW.Nghost; 
	 i_fine <= SolnBlk_Original_NW.ICu+SolnBlk_Original_NW.Nghost; i_fine += 2) {
      i_coarse = (i_fine-SolnBlk_Original_NW.ICl)/2+Local_SolnBlk[nb].ICl;
      Local_SolnBlk[nb].WoN[i_coarse] = SolnBlk_Original_NW.WoN[i_fine];
      if (i_fine == SolnBlk_Original_NW.ICl-SolnBlk_Original_NW.Nghost) {
	Local_SolnBlk[nb].WoN[i_coarse-1] = SolnBlk_Original_NW.WoN[i_fine];
      }
    }

    // North-East corner fine block:
    for (int j_fine = SolnBlk_Original_NE.JCl; j_fine <= SolnBlk_Original_NE.JCu; j_fine += 2) {
      for (int i_fine = SolnBlk_Original_NE.ICl; i_fine <= SolnBlk_Original_NE.ICu; i_fine += 2) {

	i_coarse = (i_fine-SolnBlk_Original_NE.ICl)/2 + (Local_SolnBlk[nb].ICu-Local_SolnBlk[nb].ICl+1)/2+Local_SolnBlk[nb].ICl;
	j_coarse = (j_fine-SolnBlk_Original_NE.JCl)/2 + (Local_SolnBlk[nb].JCu-Local_SolnBlk[nb].JCl+1)/2+Local_SolnBlk[nb].JCl;

	if (Mesh[nb].cell_status[i_coarse][j_coarse] != CELL_STATUS_ACTIVE) {

	  // The cell is inactive.  Set to standard atmosphere.
	  Local_SolnBlk[nb].U[i_coarse][j_coarse].Standard_Atmosphere();
	  Local_SolnBlk[nb].W[i_coarse][j_coarse].Standard_Atmosphere();

	} else if (!Interface_Union_List.Ni) {
	  // No embedded boundaries exist.  Restrict the solution
	  // information from the fine cells to the coarse cell using an
	  // area-weighted average of the fine cell solution information.
	  Local_SolnBlk[nb].U[i_coarse][j_coarse] = (SolnBlk_Original_NE.Grid.Cell[i_fine  ][j_fine  ].A*
						  SolnBlk_Original_NE.U[i_fine  ][j_fine  ] +
						  SolnBlk_Original_NE.Grid.Cell[i_fine+1][j_fine  ].A*
						  SolnBlk_Original_NE.U[i_fine+1][j_fine  ] + 
						  SolnBlk_Original_NE.Grid.Cell[i_fine  ][j_fine+1].A*
						  SolnBlk_Original_NE.U[i_fine  ][j_fine+1] +
						  SolnBlk_Original_NE.Grid.Cell[i_fine+1][j_fine+1].A*
						  SolnBlk_Original_NE.U[i_fine+1][j_fine+1])/
                                                 (SolnBlk_Original_NE.Grid.Cell[i_fine  ][j_fine  ].A +
						  SolnBlk_Original_NE.Grid.Cell[i_fine+1][j_fine  ].A +
						  SolnBlk_Original_NE.Grid.Cell[i_fine  ][j_fine+1].A +
						  SolnBlk_Original_NE.Grid.Cell[i_fine+1][j_fine+1].A);
	                                          //Local_SolnBlk[nb].Grid.Cell[i_coarse][j_coarse].A;
	  Local_SolnBlk[nb].W[i_coarse][j_coarse] = W(Local_SolnBlk[nb].U[i_coarse][j_coarse]);

	} else {

	  if (Mesh[nb].cell_type[i_coarse][j_coarse] == CELL_TYPE_QUADRILATERAL &&
	      AMesh_Original_NE.cell_type[i_fine  ][j_fine  ] == CELL_TYPE_QUADRILATERAL &&
	      AMesh_Original_NE.cell_type[i_fine+1][j_fine  ] == CELL_TYPE_QUADRILATERAL &&
	      AMesh_Original_NE.cell_type[i_fine  ][j_fine+1] == CELL_TYPE_QUADRILATERAL &&
	      AMesh_Original_NE.cell_type[i_fine+1][j_fine+1] == CELL_TYPE_QUADRILATERAL &&
	      AMesh_Original_NE.cell_status[i_fine  ][j_fine  ] == CELL_STATUS_ACTIVE &&
	      AMesh_Original_NE.cell_status[i_fine+1][j_fine  ] == CELL_STATUS_ACTIVE &&
	      AMesh_Original_NE.cell_status[i_fine  ][j_fine+1] == CELL_STATUS_ACTIVE &&
	      AMesh_Original_NE.cell_status[i_fine+1][j_fine+1] == CELL_STATUS_ACTIVE) {
	    // None of the associated cells have been adjusted.
	    // Restrict the solution information using an area-weighted
	    // average of the cells from the stored adjusted fine mesh
	    // and the coarse mesh.
	    Local_SolnBlk[nb].U[i_coarse][j_coarse] = (AGrid_Original_NE.Cell[i_fine  ][j_fine  ].A*
						       SolnBlk_Original_NE.U[i_fine  ][j_fine  ] +
						       AGrid_Original_NE.Cell[i_fine+1][j_fine  ].A*
						       SolnBlk_Original_NE.U[i_fine+1][j_fine  ] + 
						       AGrid_Original_NE.Cell[i_fine  ][j_fine+1].A*
						       SolnBlk_Original_NE.U[i_fine  ][j_fine+1] +
						       AGrid_Original_NE.Cell[i_fine+1][j_fine+1].A*
						       SolnBlk_Original_NE.U[i_fine+1][j_fine+1])/
                                                      (AGrid_Original_NE.Cell[i_fine  ][j_fine  ].A +
						       AGrid_Original_NE.Cell[i_fine+1][j_fine  ].A +
						       AGrid_Original_NE.Cell[i_fine  ][j_fine+1].A +
						       AGrid_Original_NE.Cell[i_fine+1][j_fine+1].A);
	    Local_SolnBlk[nb].W[i_coarse][j_coarse] = W(Local_SolnBlk[nb].U[i_coarse][j_coarse]);

 	  } else {
	    // The coarse and/or fine cells have been adjusted.
	    // Restrict the solution information using an area-weight
	    // average based on the area of intersection of the coarse
	    // grid and fine grid cells.
	    total_area = ZERO;
	    Local_SolnBlk[nb].U[i_coarse][j_coarse].Vacuum();
	    // Create the coarse cell polygon.
	    Pc.convert(Local_SolnBlk[nb].Grid.nodeSW(i_coarse,j_coarse).X,Local_SolnBlk[nb].Grid.nodeSE(i_coarse,j_coarse).X,
		       Local_SolnBlk[nb].Grid.nodeNE(i_coarse,j_coarse).X,Local_SolnBlk[nb].Grid.nodeNW(i_coarse,j_coarse).X);
	    // Search the fine cells for intersecting cells.
	    for (int j = j_fine-2; j < j_fine+4; j++) {
	      for (int i = i_fine-2; i < i_fine+4; i++) {
		if (AMesh_Original_NE.cell_status[i][j] == CELL_STATUS_ACTIVE) {
		  // Create the fine cell polygon.
		  Pf.convert(SolnBlk_Original_NE.Grid.nodeSW(i,j).X,SolnBlk_Original_NE.Grid.nodeSE(i,j).X,
			     SolnBlk_Original_NE.Grid.nodeNE(i,j).X,SolnBlk_Original_NE.Grid.nodeNW(i,j).X);
		  // Determine the area of intersection between the fine
		  // and coarse cell polygons.
		  A = Polygon_Intersection_Area(Pc,Pf);
		  Pf.deallocate();
		  // Add the area of intersection to the total area of
		  // of intersection.
		  total_area += A;
		  // Restrict fine cell solution to the coarse cell
		  // based on the area of interesection.
		  Local_SolnBlk[nb].U[i_coarse][j_coarse] += SolnBlk_Original_NE.U[i][j]*A;
		}
	      }
	    }
	    Pc.deallocate();
	    // Complete solution restriction by dividing the restricted
	    // solution by the total area of intersection.
	    Local_SolnBlk[nb].U[i_coarse][j_coarse] = Local_SolnBlk[nb].U[i_coarse][j_coarse]/total_area;
	    Local_SolnBlk[nb].W[i_coarse][j_coarse] = W(Local_SolnBlk[nb].U[i_coarse][j_coarse]);

	  }

	}

      }
    }

    for (int j_fine = SolnBlk_Original_NE.JCl-SolnBlk_Original_NE.Nghost; 
	 j_fine <= SolnBlk_Original_NE.JCu+SolnBlk_Original_NE.Nghost; j_fine += 2) {
      j_coarse = (j_fine-SolnBlk_Original_NE.JCl)/2 + (Local_SolnBlk[nb].JCu-Local_SolnBlk[nb].JCl+1)/2+Local_SolnBlk[nb].JCl;
      Local_SolnBlk[nb].WoE[j_coarse] = SolnBlk_Original_NE.WoE[j_fine];
      if (j_fine == SolnBlk_Original_NE.JCu+SolnBlk_Original_NE.Nghost) {
	Local_SolnBlk[nb].WoE[j_coarse+1] = SolnBlk_Original_NE.WoE[j_fine];
      }
    }

    for (int i_fine = SolnBlk_Original_NE.ICl-SolnBlk_Original_NE.Nghost; 
	 i_fine <= SolnBlk_Original_NE.ICu+SolnBlk_Original_NE.Nghost; i_fine += 2) {
      i_coarse = (i_fine-SolnBlk_Original_NE.ICl)/2 + (Local_SolnBlk[nb].ICu-Local_SolnBlk[nb].ICl+1)/2+Local_SolnBlk[nb].ICl;
      Local_SolnBlk[nb].WoN[i_coarse] = SolnBlk_Original_NE.WoN[i_fine];
      if (i_fine == SolnBlk_Original_NE.ICu+SolnBlk_Original_NE.Nghost) {
	Local_SolnBlk[nb].WoN[i_coarse+1] = SolnBlk_Original_NE.WoN[i_fine];
      }
    }

  }

  // Restriction of solution block was successful.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Refine_Grid --                               *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
inline int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Refine_Grid(const int &Perform_Mesh_Adjustment) {

  int error_flag;
  int new_CPU;
  int new_blocks_CPU[4], new_blocks_BLK[4], new_blocks_SECTOR[4];
  int my_rank, undefined_rank, number_CPUs_in_new_blocks, CPUs_in_new_blocks[4];

  Quad_Soln_Block solution_block_to_be_refined;
  Grid2D_Quad_Block agrid_block_to_be_refined;
  Adjusted_Mesh_Quad_Block amesh_block_to_be_refined;
  QuadTreeBlock *quadtree_block_to_be_refined_ptr;

#ifdef _MPI_VERSION
  MPI::Intracomm refine_comm;
  MPI::Group     big_group = MPI::COMM_WORLD.Get_group();
  MPI::Group     refine_group;
  undefined_rank = MPI::UNDEFINED;
#else
  undefined_rank = -1;
#endif

  for (int nCPU = 0; nCPU < QuadTree->Ncpu; nCPU++) { // Loop over available processors.
    for (int nBLK = 0; nBLK < QuadTree->Nblk; nBLK++) { // Loop over available blocks.
       if (QuadTree->Blocks[nCPU][nBLK] != NULL) {

	// Refine the solution block if it is in use and if it has been
	// flagged for refinement.
	if (QuadTree->Blocks[nCPU][nBLK]->block.used &&
	    QuadTree->RefineFlags[nCPU][nBLK] == ADAPTIVEBLOCK2D_REFINE) {

 	  // Get pointer to next quadtree block scheduled for mesh refinement.
 	  quadtree_block_to_be_refined_ptr = QuadTree->Blocks[nCPU][nBLK]; 

 	  // Obtain list of new solution blocks to be introduced.
 	  new_blocks_CPU[0] = quadtree_block_to_be_refined_ptr->block.info.cpu;
 	  new_blocks_BLK[0] = quadtree_block_to_be_refined_ptr->block.info.blknum;
 	  new_blocks_SECTOR[0] = ADAPTIVEBLOCK2D_SECTOR_SW;

 	  number_CPUs_in_new_blocks = 1;
 	  CPUs_in_new_blocks[0] = new_blocks_CPU[0];
 	  for (int iNEW = 1; iNEW <= 3; iNEW++) {
 	    CPUs_in_new_blocks[iNEW] = -1;
 	  }

	  for (int iNEW = 1; iNEW <= 3; iNEW++) {
	    if (Global_Solution_Block_List->Nfree > 0) {
	      // Get new blocks from global solution block list.
	      new_blocks_CPU[iNEW] = Global_Solution_Block_List->nextCPU();
	      new_blocks_BLK[iNEW] = Global_Solution_Block_List->nextBlock();
	      new_blocks_SECTOR[iNEW] = ADAPTIVEBLOCK2D_SECTOR_SW + iNEW;
	      Global_Solution_Block_List->update_next();
	    } else {
	      cout << "\n " << CFFC_Version() 
		   << " AMR Error: Refine_Grid, Insufficient number of solution blocks.\n";
	      return 6320;
	    }

	    new_CPU = 1;
	    for (int ii = 0; ii < iNEW; ii++) {
	      if (new_blocks_CPU[iNEW] == new_blocks_CPU[ii]) new_CPU = 0;
	    }
	    number_CPUs_in_new_blocks += new_CPU;
	    if (new_CPU == 1) {
	      CPUs_in_new_blocks[number_CPUs_in_new_blocks-1] = new_blocks_CPU[iNEW];
	    }
	  }

	  // Create a MPI group and communicator for all processors 
	  // involved in the refinement of this solution block.
	  my_rank = Local_Solution_Block_List->ThisCPU;
#ifdef _MPI_VERSION
	  refine_group = big_group.Incl(number_CPUs_in_new_blocks,
					CPUs_in_new_blocks);
	  refine_comm = MPI::COMM_WORLD.Create(refine_group);
	  my_rank = refine_group.Get_rank();
#endif

	  // Obtain a copy of the solution block to be refined on each
	  // processor involved in the mesh refinement.
	  if (Local_Solution_Block_List->ThisCPU == new_blocks_CPU[0]) {
	    Copy_Solution_Block(solution_block_to_be_refined,
				Local_SolnBlk[new_blocks_BLK[0]]);
	    Copy_Quad_Block(agrid_block_to_be_refined,
			    AGrid[new_blocks_BLK[0]]);
	    amesh_block_to_be_refined.Copy_Quad_Block(AMesh[new_blocks_BLK[0]]);
	  }

#ifdef _MPI_VERSION
	  if (my_rank != undefined_rank && number_CPUs_in_new_blocks > 1) {
	    Broadcast_Solution_Block(solution_block_to_be_refined,
				     refine_comm,
				     new_blocks_CPU[0]);
	    Broadcast_Quad_Block(agrid_block_to_be_refined,
				 refine_comm,
				 new_blocks_CPU[0]);
	    amesh_block_to_be_refined.Broadcast_Quad_Block(refine_comm,
							   new_blocks_CPU[0]);
	  }
#endif

 	  // Set local solution block information, create refined mesh,
 	  // and prolong solution for newly created solution blocks.
 	  if (my_rank != undefined_rank) {
 	    for (int iNEW = 0; iNEW <= 3; iNEW++) {
 	      if (Local_Solution_Block_List->ThisCPU == new_blocks_CPU[iNEW]) {
		Local_Solution_Block_List->Block[new_blocks_BLK[iNEW]].used =
		  ADAPTIVEBLOCK2D_USED;
		Local_Solution_Block_List->Block[new_blocks_BLK[iNEW]].info.cpu =
		  new_blocks_CPU[iNEW];
		Local_Solution_Block_List->Block[new_blocks_BLK[iNEW]].info.blknum =
		  new_blocks_BLK[iNEW];
		Local_Solution_Block_List->Block[new_blocks_BLK[iNEW]].info.dimen.i =
		  quadtree_block_to_be_refined_ptr->block.info.dimen.i;
		Local_Solution_Block_List->Block[new_blocks_BLK[iNEW]].info.dimen.j = 
		  quadtree_block_to_be_refined_ptr->block.info.dimen.j;
		Local_Solution_Block_List->Block[new_blocks_BLK[iNEW]].info.dimen.ghost =
		  quadtree_block_to_be_refined_ptr->block.info.dimen.ghost;
		Local_Solution_Block_List->Block[new_blocks_BLK[iNEW]].info.sector = 
		  new_blocks_SECTOR[iNEW];
		Local_Solution_Block_List->Block[new_blocks_BLK[iNEW]].info.level = 
		  quadtree_block_to_be_refined_ptr->block.info.level + 1;
		Local_Solution_Block_List->RefineFlag[new_blocks_BLK[iNEW]] = ADAPTIVEBLOCK2D_REFINE;

		if ((solution_block_to_be_refined.NCi-2*solution_block_to_be_refined.Nghost-
		     2*((solution_block_to_be_refined.NCi-2*solution_block_to_be_refined.Nghost)/2) != 0) || 
		    (solution_block_to_be_refined.NCj-2*solution_block_to_be_refined.Nghost-
		     2*((solution_block_to_be_refined.NCj-2*solution_block_to_be_refined.Nghost)/2) != 0) ||
		    (solution_block_to_be_refined.NCi-2*solution_block_to_be_refined.Nghost < 4) ||
		    (solution_block_to_be_refined.NCj-2*solution_block_to_be_refined.Nghost < 4) ||
		    (solution_block_to_be_refined.Grid.Node == NULL)) {
		  cout << "\n " << CFFC_Version() 
		       << " AMR Error: Refine_Grid, Cannot refine mesh.\n";
		  return 6321;
		}

		if (Local_SolnBlk[new_blocks_BLK[iNEW]].U != NULL) 
		  Local_SolnBlk[new_blocks_BLK[iNEW]].deallocate();
		Local_SolnBlk[new_blocks_BLK[iNEW]].allocate(solution_block_to_be_refined.NCi-2*solution_block_to_be_refined.Nghost,
							     solution_block_to_be_refined.NCj-2*solution_block_to_be_refined.Nghost,
							     solution_block_to_be_refined.Nghost);
		OGrid[new_blocks_BLK[iNEW]].allocate(solution_block_to_be_refined.NCi-2*solution_block_to_be_refined.Nghost,
						     solution_block_to_be_refined.NCj-2*solution_block_to_be_refined.Nghost,
						     solution_block_to_be_refined.Nghost);
		AGrid[new_blocks_BLK[iNEW]].allocate(solution_block_to_be_refined.NCi-2*solution_block_to_be_refined.Nghost,
						     solution_block_to_be_refined.NCj-2*solution_block_to_be_refined.Nghost,
						     solution_block_to_be_refined.Nghost);
		Mesh[new_blocks_BLK[iNEW]].allocate(solution_block_to_be_refined.NCi,
						    solution_block_to_be_refined.NCj);
		OMesh[new_blocks_BLK[iNEW]].allocate(solution_block_to_be_refined.NCi,
						     solution_block_to_be_refined.NCj);
		AMesh[new_blocks_BLK[iNEW]].allocate(solution_block_to_be_refined.NCi,
						     solution_block_to_be_refined.NCj);
		Adjustment_Data[new_blocks_BLK[iNEW]].allocate(solution_block_to_be_refined.NCi,
							       solution_block_to_be_refined.NCj,
							       solution_block_to_be_refined.Nghost,
							       Interface_Union_List.Ni);
		Refine_Mesh(Local_SolnBlk[new_blocks_BLK[iNEW]].Grid,
			    solution_block_to_be_refined.Grid,
			    new_blocks_SECTOR[iNEW]);
		if (Check_Quad_Block(Local_SolnBlk[new_blocks_BLK[iNEW]].Grid)) { 
		  cout << "\n " << CFFC_Version() 
		       << " AMR Error: Refine_Grid, Invalid refined mesh.\n";
		  return 6322;
		}
 		error_flag = Prolong_Quadrilateral_Block(Perform_Mesh_Adjustment,
							 new_blocks_BLK[iNEW],
 							 solution_block_to_be_refined,
 							 agrid_block_to_be_refined,
 							 amesh_block_to_be_refined,
 							 new_blocks_SECTOR[iNEW]);
 		if (error_flag) return error_flag;
 		if ((IP->i_Smooth_Quad_Block && !Adjustment_Data[new_blocks_BLK[iNEW]].Interface_Present[0]) &&
		    !(IP->i_Grid == GRID_ROCKET_MOTOR && Local_SolnBlk[new_blocks_BLK[iNEW]].Grid.Node[3][2].X.x < ZERO)) {
 		  Smooth_Quad_Block(Local_SolnBlk[new_blocks_BLK[iNEW]].Grid,
 				    min(250,4*max(solution_block_to_be_refined.NCi-2*solution_block_to_be_refined.Nghost,
 						  solution_block_to_be_refined.NCi-2*solution_block_to_be_refined.Nghost)));
 		  Smooth_Quad_Block(OGrid[new_blocks_BLK[iNEW]],
 				    min(250,4*max(solution_block_to_be_refined.NCi-2*solution_block_to_be_refined.Nghost,
 						  solution_block_to_be_refined.NCi-2*solution_block_to_be_refined.Nghost)));
 		}

 	      }
 	    }
 	  }

	  // Release the MPI group and communicator.
#ifdef _MPI_VERSION
	  if (refine_comm != MPI::COMM_NULL) refine_comm.Free();
	  refine_group.Free();
#endif
	  my_rank = undefined_rank;

	  // Update and assign quadtree solution block information for
	  // newly created solution blocks.
	  QuadTree->refineBlock(new_blocks_CPU, 
				new_blocks_BLK,
				new_blocks_SECTOR);

	  // Finally, reset refinement flag.
	  QuadTree->RefineFlags[nCPU][nBLK] = ADAPTIVEBLOCK2D_NOCHANGE;

	} else if (!QuadTree->Blocks[nCPU][nBLK]->block.used) {
	  // Block not used, reset refinement flag.
	  QuadTree->RefineFlags[nCPU][nBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
	}

      } else {
	// Block does not exist, resent refinement flag.
	QuadTree->RefineFlags[nCPU][nBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
      }
    }
  }

  // Mesh refinement complete.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Coarsen_Grid --                              *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
inline int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Coarsen_Grid(const int &Perform_Mesh_Adjustment) {

  int error_flag;
  int old_blocks_CPU[4], old_blocks_BLK[4], old_blocks_SECTOR[4], CPU_list[2];
  int my_rank, undefined_rank, all_coarsen_flag;

  Quad_Soln_Block solution_block_to_be_coarsened_SW_sibling,
                  solution_block_to_be_coarsened_SE_sibling,
                  solution_block_to_be_coarsened_NW_sibling,
                  solution_block_to_be_coarsened_NE_sibling;
  Grid2D_Quad_Block agrid_block_to_be_coarsened_SW_sibling,
                    agrid_block_to_be_coarsened_SE_sibling,
                    agrid_block_to_be_coarsened_NW_sibling,
                    agrid_block_to_be_coarsened_NE_sibling;
  Adjusted_Mesh_Quad_Block amesh_block_to_be_coarsened_SW_sibling,
                           amesh_block_to_be_coarsened_SE_sibling,
                           amesh_block_to_be_coarsened_NW_sibling,
                           amesh_block_to_be_coarsened_NE_sibling;
  QuadTreeBlock *quadtree_block_to_be_coarsened_ptr; 

#ifdef _MPI_VERSION
  MPI::Intracomm coarsening_comm_SE, coarsening_comm_NW, coarsening_comm_NE;
  MPI::Group     big_group = MPI::COMM_WORLD.Get_group();
  MPI::Group     coarsening_group_SE, coarsening_group_NW, coarsening_group_NE;
  undefined_rank = MPI::UNDEFINED;
#else
  undefined_rank = -1;
#endif

  for (int nCPU = 0; nCPU < QuadTree->Ncpu; nCPU++) { // Loop over available processors.
    for (int nBLK = 0; nBLK < QuadTree->Nblk; nBLK++) { // Loop over available blocks.
      if (QuadTree->Blocks[nCPU][nBLK] != NULL) {
	if (QuadTree->Blocks[nCPU][nBLK]->block.used && // Coarsen only solution blocks in use.
	    QuadTree->RefineFlags[nCPU][nBLK] == ADAPTIVEBLOCK2D_COARSEN) { // Check each solution block for coarsening.
	   
	  // Get pointer to next quadtree block scheduled for mesh coarsening.
	  quadtree_block_to_be_coarsened_ptr = QuadTree->Blocks[nCPU][nBLK];

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
		for (int iOLD = 0; iOLD <= 3; iOLD++) {
		  if (QuadTree->RefineFlags[old_blocks_CPU[iOLD]][old_blocks_BLK[iOLD]] != ADAPTIVEBLOCK2D_COARSEN) {
		    all_coarsen_flag = 0;
		    break;
		  }
		}
		if (all_coarsen_flag) {
		  // Make copies of each of the four sibling solution blocks involved 
		  // in the coarsening process.
		  if (Local_Solution_Block_List->ThisCPU == old_blocks_CPU[0]) {
		    Copy_Solution_Block(solution_block_to_be_coarsened_SW_sibling,
					Local_SolnBlk[old_blocks_BLK[0]]);
		    Copy_Quad_Block(agrid_block_to_be_coarsened_SW_sibling,
				    AGrid[old_blocks_BLK[0]]);
		    amesh_block_to_be_coarsened_SW_sibling.Copy_Quad_Block(AMesh[old_blocks_BLK[0]]);
		  }

		  if (Local_Solution_Block_List->ThisCPU == old_blocks_CPU[1]) {
		    Copy_Solution_Block(solution_block_to_be_coarsened_SE_sibling,
					Local_SolnBlk[old_blocks_BLK[1]]);
		    Copy_Quad_Block(agrid_block_to_be_coarsened_SE_sibling,
				    AGrid[old_blocks_BLK[1]]);
		    amesh_block_to_be_coarsened_SE_sibling.Copy_Quad_Block(AMesh[old_blocks_BLK[1]]);
		  }

		  if (Local_Solution_Block_List->ThisCPU == old_blocks_CPU[2]) {
		    Copy_Solution_Block(solution_block_to_be_coarsened_NW_sibling,
					Local_SolnBlk[old_blocks_BLK[2]]);
		    Copy_Quad_Block(agrid_block_to_be_coarsened_NW_sibling,
				    AGrid[old_blocks_BLK[2]]);
		    amesh_block_to_be_coarsened_NW_sibling.Copy_Quad_Block(AMesh[old_blocks_BLK[2]]);
		  }

		  if (Local_Solution_Block_List->ThisCPU == old_blocks_CPU[3]) {
		    Copy_Solution_Block(solution_block_to_be_coarsened_NE_sibling,
					Local_SolnBlk[old_blocks_BLK[3]]);
		    Copy_Quad_Block(agrid_block_to_be_coarsened_NE_sibling,
				    AGrid[old_blocks_BLK[3]]);
		    amesh_block_to_be_coarsened_NE_sibling.Copy_Quad_Block(AMesh[old_blocks_BLK[3]]);
		  }
      
		  // Create MPI groups and communicators for all processors 
		  // involved in the coarsening of the four sibling solution blocks.
#ifdef _MPI_VERSION
		  if (old_blocks_CPU[0] != old_blocks_CPU[1]) {
		    CPU_list[0] = old_blocks_CPU[1];
		    CPU_list[1] = old_blocks_CPU[0];
		    coarsening_group_SE = big_group.Incl(2,CPU_list);
		    coarsening_comm_SE = MPI::COMM_WORLD.Create(coarsening_group_SE);
		  }
		  if (old_blocks_CPU[0] != old_blocks_CPU[2]) {
		    CPU_list[0] = old_blocks_CPU[2];
		    CPU_list[1] = old_blocks_CPU[0];
		    coarsening_group_NW = big_group.Incl(2,CPU_list);
		    coarsening_comm_NW = MPI::COMM_WORLD.Create(coarsening_group_NW);
		  }
		  if (old_blocks_CPU[0] != old_blocks_CPU[3]) {
		    CPU_list[0] = old_blocks_CPU[3];
		    CPU_list[1] = old_blocks_CPU[0];
		    coarsening_group_NE = big_group.Incl(2,CPU_list);
		    coarsening_comm_NE = MPI::COMM_WORLD.Create(coarsening_group_NE);
		  }
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
		      Broadcast_Quad_Block(agrid_block_to_be_coarsened_SE_sibling,
					   coarsening_comm_SE,
					   old_blocks_CPU[1]);
		      amesh_block_to_be_coarsened_SE_sibling.Broadcast_Quad_Block(coarsening_comm_SE,
										  old_blocks_CPU[1]);
		    }
		  }
		  if (old_blocks_CPU[0] != old_blocks_CPU[2]) {
		    my_rank = coarsening_group_NW.Get_rank();
		    if (my_rank != undefined_rank) {
		      Broadcast_Solution_Block(solution_block_to_be_coarsened_NW_sibling,
					       coarsening_comm_NW,
					       old_blocks_CPU[2]);
		      Broadcast_Quad_Block(agrid_block_to_be_coarsened_NW_sibling,
					   coarsening_comm_NW,
					   old_blocks_CPU[2]);
		      amesh_block_to_be_coarsened_NW_sibling.Broadcast_Quad_Block(coarsening_comm_NW,
										  old_blocks_CPU[2]);
		    }
		  }
		  if (old_blocks_CPU[0] != old_blocks_CPU[3]) {
		    my_rank = coarsening_group_NE.Get_rank();
		    if (my_rank != undefined_rank) {
		      Broadcast_Solution_Block(solution_block_to_be_coarsened_NE_sibling,
					       coarsening_comm_NE,
					       old_blocks_CPU[3]);
		      Broadcast_Quad_Block(agrid_block_to_be_coarsened_NE_sibling,
					   coarsening_comm_NE,
					   old_blocks_CPU[3]);
		      amesh_block_to_be_coarsened_NE_sibling.Broadcast_Quad_Block(coarsening_comm_NE,
										  old_blocks_CPU[3]);
		    }
		  }
#endif

		  // Set local solution block information, create coarse mesh,
		  // and restrict solution to newly created coarse solution block.
		  if (Local_Solution_Block_List->ThisCPU == old_blocks_CPU[0]) {
		    Local_Solution_Block_List->Block[old_blocks_BLK[0]].used = 
		      ADAPTIVEBLOCK2D_USED;
		    Local_Solution_Block_List->Block[old_blocks_BLK[0]].info.cpu = 
		      old_blocks_CPU[0];
		    Local_Solution_Block_List->Block[old_blocks_BLK[0]].info.blknum = 
		      old_blocks_BLK[0];
		    Local_Solution_Block_List->Block[old_blocks_BLK[0]].info.dimen.i = 
		      QuadTree->Blocks[old_blocks_CPU[0]][old_blocks_BLK[0]]->block.info.dimen.i;
		    Local_Solution_Block_List->Block[old_blocks_BLK[0]].info.dimen.j = 
		      QuadTree->Blocks[old_blocks_CPU[0]][old_blocks_BLK[0]]->block.info.dimen.j;
		    Local_Solution_Block_List->Block[old_blocks_BLK[0]].info.dimen.ghost =
		      QuadTree->Blocks[old_blocks_CPU[0]][old_blocks_BLK[0]]->block.info.dimen.ghost;
		    Local_Solution_Block_List->Block[old_blocks_BLK[0]].info.sector = 
		      old_blocks_SECTOR[0]; // Note that this sector info is incorrect, 
		    // but subsequent call to Find_Neighbours will fix this!
		    Local_Solution_Block_List->Block[old_blocks_BLK[0]].info.level = 
		      QuadTree->Blocks[old_blocks_CPU[0]][old_blocks_BLK[0]]->block.info.level - 1;
		    Local_Solution_Block_List->RefineFlag[old_blocks_BLK[0]] = ADAPTIVEBLOCK2D_COARSEN;

		    if ((solution_block_to_be_coarsened_SW_sibling.NCi-2*solution_block_to_be_coarsened_SW_sibling.Nghost-
			 2*((solution_block_to_be_coarsened_SW_sibling.NCi-2*solution_block_to_be_coarsened_SW_sibling.Nghost)/2) != 0) ||
			(solution_block_to_be_coarsened_SW_sibling.NCj-2*solution_block_to_be_coarsened_SW_sibling.Nghost-
			 2*((solution_block_to_be_coarsened_SW_sibling.NCj-2*solution_block_to_be_coarsened_SW_sibling.Nghost)/2) != 0) ||
			(solution_block_to_be_coarsened_SW_sibling.NCi-2*solution_block_to_be_coarsened_SW_sibling.Nghost < 4) ||
			(solution_block_to_be_coarsened_SW_sibling.NCj-2*solution_block_to_be_coarsened_SW_sibling.Nghost < 4) ||
			(solution_block_to_be_coarsened_SW_sibling.Grid.Node == NULL)) {
		      cout << "\n " << CFFC_Version()
			   << " AMR Error: Coarsen_Grid, Cannot coarsen south-west mesh.\n";
		      return 7480;
		    }

		    if ((solution_block_to_be_coarsened_SE_sibling.NCi-2*solution_block_to_be_coarsened_SE_sibling.Nghost-
			 2*((solution_block_to_be_coarsened_SE_sibling.NCi-2*solution_block_to_be_coarsened_SE_sibling.Nghost)/2) != 0) ||
			(solution_block_to_be_coarsened_SE_sibling.NCj-2*solution_block_to_be_coarsened_SE_sibling.Nghost-
			 2*((solution_block_to_be_coarsened_SE_sibling.NCj-2*solution_block_to_be_coarsened_SE_sibling.Nghost)/2) != 0) ||
			(solution_block_to_be_coarsened_SE_sibling.NCi-2*solution_block_to_be_coarsened_SE_sibling.Nghost < 4) ||
			(solution_block_to_be_coarsened_SE_sibling.NCj-2*solution_block_to_be_coarsened_SE_sibling.Nghost < 4) ||
			(solution_block_to_be_coarsened_SE_sibling.Grid.Node == NULL)) {
		      cout << "\n " << CFFC_Version()
			   << " AMR Error: Coarsen_Grid, Cannot coarsen south-east mesh.\n";
		      return 7481;
		    }

		    if ((solution_block_to_be_coarsened_NW_sibling.NCi-2*solution_block_to_be_coarsened_NW_sibling.Nghost-
			 2*((solution_block_to_be_coarsened_NW_sibling.NCi-2*solution_block_to_be_coarsened_NW_sibling.Nghost)/2) != 0) ||
			(solution_block_to_be_coarsened_NW_sibling.NCj-2*solution_block_to_be_coarsened_NW_sibling.Nghost-
			 2*((solution_block_to_be_coarsened_NW_sibling.NCj-2*solution_block_to_be_coarsened_NW_sibling.Nghost)/2) != 0) ||
			(solution_block_to_be_coarsened_NW_sibling.NCi-2*solution_block_to_be_coarsened_NW_sibling.Nghost < 4) ||
			(solution_block_to_be_coarsened_NW_sibling.NCj-2*solution_block_to_be_coarsened_NW_sibling.Nghost < 4) ||
			(solution_block_to_be_coarsened_NW_sibling.Grid.Node == NULL)) {
		      cout << "\n " << CFFC_Version()
			   << " AMR Error: Coarsen_Grid, Cannot coarsen north-west mesh.\n";
		      return 7482;
		    }

		    if ((solution_block_to_be_coarsened_NE_sibling.NCi-2*solution_block_to_be_coarsened_NE_sibling.Nghost-
			 2*((solution_block_to_be_coarsened_NE_sibling.NCi-2*solution_block_to_be_coarsened_NE_sibling.Nghost)/2) != 0) ||
			(solution_block_to_be_coarsened_NE_sibling.NCj-2*solution_block_to_be_coarsened_NE_sibling.Nghost-
			 2*((solution_block_to_be_coarsened_NE_sibling.NCj-2*solution_block_to_be_coarsened_NE_sibling.Nghost)/2) != 0) ||
			(solution_block_to_be_coarsened_NE_sibling.NCi-2*solution_block_to_be_coarsened_NE_sibling.Nghost < 4) ||
			(solution_block_to_be_coarsened_NE_sibling.NCj-2*solution_block_to_be_coarsened_NE_sibling.Nghost < 4) ||
			(solution_block_to_be_coarsened_NE_sibling.Grid.Node == NULL)) {
		      cout << "\n " << CFFC_Version()
			   << " AMR Error: Coarsen_Grid, Cannot coarsen north-east mesh.\n";
		      return 7483;
		    }

		    if (Local_SolnBlk[old_blocks_BLK[0]].U != NULL) 
		      Local_SolnBlk[old_blocks_BLK[0]].deallocate();
		    Local_SolnBlk[old_blocks_BLK[0]].allocate(solution_block_to_be_coarsened_SW_sibling.NCi-2*solution_block_to_be_coarsened_SW_sibling.Nghost,
							      solution_block_to_be_coarsened_SW_sibling.NCj-2*solution_block_to_be_coarsened_SW_sibling.Nghost,
							      solution_block_to_be_coarsened_SW_sibling.Nghost);
		    OGrid[old_blocks_BLK[0]].allocate(solution_block_to_be_coarsened_SW_sibling.NCi-2*solution_block_to_be_coarsened_SW_sibling.Nghost,
						      solution_block_to_be_coarsened_SW_sibling.NCj-2*solution_block_to_be_coarsened_SW_sibling.Nghost,
						      solution_block_to_be_coarsened_SW_sibling.Nghost);
		    AGrid[old_blocks_BLK[0]].allocate(solution_block_to_be_coarsened_SW_sibling.NCi-2*solution_block_to_be_coarsened_SW_sibling.Nghost,
						      solution_block_to_be_coarsened_SW_sibling.NCj-2*solution_block_to_be_coarsened_SW_sibling.Nghost,
						      solution_block_to_be_coarsened_SW_sibling.Nghost);
		    Mesh[old_blocks_BLK[0]].allocate(solution_block_to_be_coarsened_SW_sibling.NCi,
						     solution_block_to_be_coarsened_SW_sibling.NCj);
		    OMesh[old_blocks_BLK[0]].allocate(solution_block_to_be_coarsened_SW_sibling.NCi,
						      solution_block_to_be_coarsened_SW_sibling.NCj);
		    AMesh[old_blocks_BLK[0]].allocate(solution_block_to_be_coarsened_SW_sibling.NCi,
						      solution_block_to_be_coarsened_SW_sibling.NCj);
		    Adjustment_Data[old_blocks_BLK[0]].allocate(solution_block_to_be_coarsened_SW_sibling.NCi,
								solution_block_to_be_coarsened_SW_sibling.NCj,
								solution_block_to_be_coarsened_SW_sibling.Nghost,
								Interface_Union_List.Ni);
		    Coarsen_Mesh(Local_SolnBlk[old_blocks_BLK[0]].Grid,
				 solution_block_to_be_coarsened_SW_sibling.Grid,
				 solution_block_to_be_coarsened_SE_sibling.Grid,
				 solution_block_to_be_coarsened_NW_sibling.Grid,
				 solution_block_to_be_coarsened_NE_sibling.Grid);
		    if (Check_Quad_Block(Local_SolnBlk[old_blocks_BLK[0]].Grid)) { 
		      cout << "\n " << CFFC_Version() 
			   << " AMR Error: Coarsen_Grid, Invalid coarsened mesh.\n";
		      return 7484;
		    }
		    error_flag = Restrict_Quadrilateral_Block(Perform_Mesh_Adjustment,
							      old_blocks_BLK[0],
							      solution_block_to_be_coarsened_SW_sibling,
							      solution_block_to_be_coarsened_SE_sibling,
							      solution_block_to_be_coarsened_NW_sibling,
							      solution_block_to_be_coarsened_NE_sibling,
							      agrid_block_to_be_coarsened_SW_sibling,
							      agrid_block_to_be_coarsened_SE_sibling,
							      agrid_block_to_be_coarsened_NW_sibling,
							      agrid_block_to_be_coarsened_NE_sibling,
							      amesh_block_to_be_coarsened_SW_sibling,
							      amesh_block_to_be_coarsened_SE_sibling,
							      amesh_block_to_be_coarsened_NW_sibling,
							      amesh_block_to_be_coarsened_NE_sibling);
		    if (error_flag) return error_flag;
		    if ((IP->i_Smooth_Quad_Block && !Adjustment_Data[old_blocks_BLK[0]].Interface_Present[0]) &&
			!(IP->i_Grid == GRID_ROCKET_MOTOR && Local_SolnBlk[old_blocks_BLK[0]].Grid.Node[3][2].X.x < ZERO)) {
		      //if (IP->i_Smooth_Quad_Block && !Adjustment_Data[old_blocks_BLK[0]].Interface_Present[0]) {
 		      Smooth_Quad_Block(Local_SolnBlk[old_blocks_BLK[0]].Grid,
 					min(250,4*max(solution_block_to_be_coarsened_SW_sibling.NCi-2*solution_block_to_be_coarsened_SW_sibling.Nghost,
 						      solution_block_to_be_coarsened_SW_sibling.NCi-2*solution_block_to_be_coarsened_SW_sibling.Nghost)));
 		      Smooth_Quad_Block(OGrid[old_blocks_BLK[0]],
 					min(250,4*max(solution_block_to_be_coarsened_SW_sibling.NCi-2*solution_block_to_be_coarsened_SW_sibling.Nghost,
 						      solution_block_to_be_coarsened_SW_sibling.NCi-2*solution_block_to_be_coarsened_SW_sibling.Nghost)));
 		    }
		  }

		  // Release the MPI groups and communicators.
#ifdef _MPI_VERSION
		  if (old_blocks_CPU[0] != old_blocks_CPU[1]) {
		    if (coarsening_comm_SE != MPI::COMM_NULL) coarsening_comm_SE.Free();
		    coarsening_group_SE.Free();
		  }
		  if (old_blocks_CPU[0] != old_blocks_CPU[2]) {
		    if (coarsening_comm_NW != MPI::COMM_NULL) coarsening_comm_NW.Free();
		    coarsening_group_NW.Free();
		  }
		  if (old_blocks_CPU[0] != old_blocks_CPU[3]) {
		    if (coarsening_comm_NE != MPI::COMM_NULL) coarsening_comm_NE.Free();
		    coarsening_group_NE.Free();
		  }
#endif

		  // Free old unused solution blocks and return to global
		  // solution block list for re-use.                          
		  for (int iOLD = 3; iOLD >= 1; iOLD--) {
		    if (Local_Solution_Block_List->ThisCPU == old_blocks_CPU[iOLD]) {
		      Local_Solution_Block_List->Block[old_blocks_BLK[iOLD]].used = ADAPTIVEBLOCK2D_NOT_USED;
		      if (Local_SolnBlk[old_blocks_BLK[iOLD]].U != NULL) Local_SolnBlk[old_blocks_BLK[iOLD]].deallocate();
		    }
		    // Return blocks and update global solution block list.
		    Global_Solution_Block_List->returnCPU(old_blocks_CPU[iOLD]);
		    Global_Solution_Block_List->returnBlock(old_blocks_BLK[iOLD]);
		    Global_Solution_Block_List->update_return();
		  }

		  // Update and assign quadtree solution block information for
		  // coarsened solution block, deleting old refined blocks.
		  QuadTree->coarsenBlocks(old_blocks_CPU, 
					  old_blocks_BLK,
					  old_blocks_SECTOR);

		  // Finally, reset refinement flags for four blocks involved in coarsening.
		  for (int iOLD = 0; iOLD <= 3; iOLD++) {
		    QuadTree->RefineFlags[old_blocks_CPU[iOLD]][old_blocks_BLK[iOLD]] = ADAPTIVEBLOCK2D_NOCHANGE;
		  }

		} else {
		  // If siblings aren't all flagged for coarsening, reset refinement flags.
		  QuadTree->RefineFlags[nCPU][nBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
		  if (Local_Solution_Block_List->ThisCPU == nCPU)
		    Local_Solution_Block_List->RefineFlag[nBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
		}

	      } else {
		// If siblings aren't all used, reset refinement flag.
		QuadTree->RefineFlags[nCPU][nBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
		if (Local_Solution_Block_List->ThisCPU == nCPU)
		  Local_Solution_Block_List->RefineFlag[nBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
	      }
                     
	    } else {
	      // If siblings don't all exist, reset refinement flag.
	      QuadTree->RefineFlags[nCPU][nBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
	      if (Local_Solution_Block_List->ThisCPU == nCPU)
		Local_Solution_Block_List->RefineFlag[nBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
	    }

	  } else {
	    // If parent doesn't exist, reset refinement flag.
	    QuadTree->RefineFlags[nCPU][nBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
	    if (Local_Solution_Block_List->ThisCPU == nCPU)
	      Local_Solution_Block_List->RefineFlag[nBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
	  }

	} else if (!QuadTree->Blocks[nCPU][nBLK]->block.used) {
	  // Block not used, reset refinement flag.
	  QuadTree->RefineFlags[nCPU][nBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
	  if (Local_Solution_Block_List->ThisCPU == nCPU)
	    Local_Solution_Block_List->RefineFlag[nBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
	}

      } else {
	// Block does not exist, reset refinement flag.
	QuadTree->RefineFlags[nCPU][nBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
	if (Local_Solution_Block_List->ThisCPU == nCPU)
	  Local_Solution_Block_List->RefineFlag[nBLK] = ADAPTIVEBLOCK2D_NOCHANGE;
      }

    }
  }

  // Mesh coarsening complete.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Calculate_Refinement_Criteria --             *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
inline void EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Calculate_Refinement_Criteria(const int &nb,
			      double *refinement_criteria) {

  // Initialize the refinement criteria for the solution block.
  for (int n_criteria = 0; n_criteria < IP->Number_of_Refinement_Criteria; n_criteria++) {
    refinement_criteria[n_criteria] = ZERO;
  }

  // Exit immediately if no active cells are present in the current block.
  if (!Adjustment_Data[nb].Number_of_Active_Cells()) {
    for (int n_criteria = 0; n_criteria < IP->Number_of_Refinement_Criteria; n_criteria++) {
      refinement_criteria[n_criteria] = HALF*(QuadTree->RefineThreshold +
					      QuadTree->CoarsenThreshold);
    }
    return ;
  }

  // Calculate the refinement criteria for each cell of the 
  // computational mesh and assign the maximum value for all cells as 
  // the refinement criteria for the solution block.
  for (int j = Local_SolnBlk[nb].JCl-1; j <= Local_SolnBlk[nb].JCu+1; j++) {
    for (int i = Local_SolnBlk[nb].ICl-1; i <= Local_SolnBlk[nb].ICu+1; i++) {

      if ((j == Local_SolnBlk[nb].JCu+1 && Local_SolnBlk[nb].Grid.BCtypeN[i] != BC_NONE) ||
	  (j == Local_SolnBlk[nb].JCl-1 && Local_SolnBlk[nb].Grid.BCtypeS[i] != BC_NONE) ||
	  (i == Local_SolnBlk[nb].ICu+1 && Local_SolnBlk[nb].Grid.BCtypeE[j] != BC_NONE) ||
	  (i == Local_SolnBlk[nb].ICl-1 && Local_SolnBlk[nb].Grid.BCtypeW[j] != BC_NONE)) {

      } else if (Mesh[nb].cell_status[i][j] == CELL_STATUS_ACTIVE) {

	// Determine the refinement criteria for the cell.
	//Calculate_Refinement_Criteria2(Local_SolnBlk[nb],*IP,i,j,refinement_criteria);
 	Calculate_Refinement_Criteria(nb,i,j,refinement_criteria);

      }

    }
  }

}

/**********************************************************************
 * EmbeddedBoundaries2D::Calculate_Refinement_Criteria --             *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
inline void EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Calculate_Refinement_Criteria(const int &nb, const int &i, const int &j,
			      double *refinement_criteria) {
  cout << endl << " ERROR: Explicit specialization required for the calculation of the";
  cout << endl << "        refinement criteria for a specific cell of a solution block.";
  cout.flush();
  assert(1==0);
}

/**********************************************************************
 * EmbeddedBoundaries2D::Flag_Blocks_For_Refinement --                *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
inline void EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Flag_Blocks_For_Refinement(void) {

  double **local_block_refinement_criteria,
          *local_max_refinement_criteria, 
          *local_min_refinement_criteria,
          *global_max_refinement_criteria,
          *global_min_refinement_criteria,
          *threshold_refinement,
          *threshold_coarsening,
          *scale_refinement,
          *scale_coarsening;

  // Allocate memory for refinement criteria.
  local_block_refinement_criteria = new double*[Local_Solution_Block_List->Nblk];
  for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
    local_block_refinement_criteria[nb] = new double[IP->Number_of_Refinement_Criteria];
  }
  local_max_refinement_criteria = new double[IP->Number_of_Refinement_Criteria];
  local_min_refinement_criteria = new double[IP->Number_of_Refinement_Criteria];
  global_max_refinement_criteria = new double[IP->Number_of_Refinement_Criteria];
  global_min_refinement_criteria = new double[IP->Number_of_Refinement_Criteria];
  threshold_refinement = new double[IP->Number_of_Refinement_Criteria];
  threshold_coarsening = new double[IP->Number_of_Refinement_Criteria];
  scale_refinement = new double[IP->Number_of_Refinement_Criteria];
  scale_coarsening = new double[IP->Number_of_Refinement_Criteria];

  // Set the solution block refinement flags to default values (no
  // change, no division and coarsening).
  Local_Solution_Block_List->nochangeAll();

  // Determine and assign the refinment flag for each of the local
  // solution blocks.

  // Initial min, max, and threshold values for each refinement criteria.
  for (int n_criteria = 0; n_criteria < IP->Number_of_Refinement_Criteria; n_criteria++) {
    local_max_refinement_criteria[n_criteria] = -1.0e-100;
    local_min_refinement_criteria[n_criteria] = 1.0e100;
    global_max_refinement_criteria[n_criteria] = -1.0e-100;
    global_min_refinement_criteria[n_criteria] = 1.0e100;
  }

  // Determine the min and max values of the refinement criteria for all
  // local solution blocks.
  for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
    if (Local_Solution_Block_List->Block[nb].used) {
      Calculate_Refinement_Criteria(nb,
				    local_block_refinement_criteria[nb]);
      for (int n_criteria = 0; n_criteria < IP->Number_of_Refinement_Criteria; n_criteria++) {
	local_max_refinement_criteria[n_criteria] = max(local_max_refinement_criteria[n_criteria],
							local_block_refinement_criteria[nb][n_criteria]);
	local_min_refinement_criteria[n_criteria] = min(local_min_refinement_criteria[n_criteria],
							local_block_refinement_criteria[nb][n_criteria]);
      }
    }
  }

  // Determine global values of the min and max refinement criteria.
#ifdef _MPI_VERSION
  MPI::COMM_WORLD.Allreduce(local_min_refinement_criteria, 
			    global_min_refinement_criteria, 
			    IP->Number_of_Refinement_Criteria,
			    MPI::DOUBLE, 
			    MPI::MIN);
  MPI::COMM_WORLD.Allreduce(local_max_refinement_criteria, 
			    global_max_refinement_criteria, 
			    IP->Number_of_Refinement_Criteria,
			    MPI::DOUBLE, 
			    MPI::MAX);
#else
  for (int n_criteria = 0; n_criteria < IP->Number_of_Refinement_Criteria; n_criteria++) {
    global_max_refinement_criteria[n_criteria] = local_max_refinement_criteria[n_criteria];
    global_min_refinement_criteria[n_criteria] = local_min_refinement_criteria[n_criteria];
  }
#endif

  // Assign the scales and thresholds for the refinement process.
  for (int n_criteria = 0; n_criteria < IP->Number_of_Refinement_Criteria; n_criteria++) {
    scale_refinement[n_criteria] = QuadTree->RefineThreshold;
    scale_coarsening[n_criteria] = QuadTree->CoarsenThreshold;
  }

  for (int n_criteria = 0; n_criteria < IP->Number_of_Refinement_Criteria; n_criteria++) {
    if (fabs(global_max_refinement_criteria[n_criteria]-
	     global_min_refinement_criteria[n_criteria]) > TOLER) {
      threshold_refinement[n_criteria] = global_min_refinement_criteria[n_criteria] +
	                                 scale_refinement[n_criteria]*
	                                 (global_max_refinement_criteria[n_criteria]-
					  global_min_refinement_criteria[n_criteria]);
      threshold_coarsening[n_criteria] = global_min_refinement_criteria[n_criteria] +
	                                 scale_coarsening[n_criteria]*
	                                 (global_max_refinement_criteria[n_criteria]-
					  global_min_refinement_criteria[n_criteria]);
    } else {
      threshold_refinement[n_criteria] = 0.90*global_min_refinement_criteria[n_criteria];
      threshold_coarsening[n_criteria] = 1.10*global_min_refinement_criteria[n_criteria];
    }
  }

  // Flag the blocks for coarsening and division.
  for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
    if (Local_Solution_Block_List->Block[nb].used) {
      for (int n_criteria = 0; n_criteria < IP->Number_of_Refinement_Criteria; n_criteria++) {
	if (local_block_refinement_criteria[nb][n_criteria] > threshold_refinement[n_criteria] &&
	    Local_Solution_Block_List->Block[nb].info.level < QuadTree->MaximumRefinementLevel) {
	  Local_Solution_Block_List->RefineFlag[nb] = ADAPTIVEBLOCK2D_REFINE;
	} else if (local_block_refinement_criteria[nb][n_criteria] < threshold_coarsening[n_criteria] &&
		   Local_Solution_Block_List->Block[nb].info.level > QuadTree->MinimumRefinementLevel &&
		   Local_Solution_Block_List->RefineFlag[nb] != ADAPTIVEBLOCK2D_REFINE) {
	  Local_Solution_Block_List->RefineFlag[nb] = ADAPTIVEBLOCK2D_COARSEN;
	} else if (Local_Solution_Block_List->RefineFlag[nb] != ADAPTIVEBLOCK2D_REFINE) {
	  Local_Solution_Block_List->RefineFlag[nb] = ADAPTIVEBLOCK2D_NOCHANGE;
	}  // if the flag was set to "refine" by another criterion, it will not be changed.
//  	cout << endl
// 	     << "REFINEMENT CRITERIA:"
// 	     << " " << Local_Solution_Block_List->Block[nb].gblknum
// 	     << " " << n_criteria
//  	     << " " << local_block_refinement_criteria[nb][n_criteria]
// //  	     << " " << global_min_refinement_criteria[n_criteria]
// //  	     << " " << global_max_refinement_criteria[n_criteria]
//  	     << " " << threshold_refinement[n_criteria]
//  	     << " " << threshold_coarsening[n_criteria]
//  	     << " " << Local_Solution_Block_List->RefineFlag[nb]; cout.flush();
      }
    }
  }

  // Deallocate memory for refinement criteria.
  for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
    delete []local_block_refinement_criteria[nb]; local_block_refinement_criteria[nb] = NULL;
  }
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

/**********************************************************************
 * EmbeddedBoundaries2D::Adaptive_Mesh_Refinement --                  *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
inline int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Adaptive_Mesh_Refinement(const int &Set_New_Refinement_Flags,
			 const int &Perform_Mesh_Adjustment) {

  int error_flag, number_of_changes,
      *maximum_interface_mesh_refinement_flag,
      *maximum_interface_mesh_refinement_level;
  int iii, jjj, intersections;
  Vector2D NE, NW, SE, SW, Xs1, Xs2;

  // Calculate the refinement measures for each solution block and flag
  // local solution blocks for refinement or coarsening if required.
  if (Set_New_Refinement_Flags) {
    Flag_Blocks_For_Refinement();
    for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
      if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {
 	if (Adjustment_Data[nb].Number_of_Active_Cells() == 0) {
 	  Local_Solution_Block_List->RefineFlag[nb] = ADAPTIVEBLOCK2D_NOCHANGE;
 	}
      }
    }    
  }

  // Enforce the interface refinement criteria if necessary.
  if (Interface_Component_List.Ni && Perform_Mesh_Adjustment) {

    // Allocate memory for the interface refinement criteria.
    maximum_interface_mesh_refinement_flag = new int[Interface_Union_List.Ni+1];
    maximum_interface_mesh_refinement_level = new int[Interface_Union_List.Ni+1];
    for (int n = 0; n < Interface_Union_List.Ni+1; n++) {
      maximum_interface_mesh_refinement_flag[n] = ADAPTIVEBLOCK2D_COARSEN;
      maximum_interface_mesh_refinement_level[n] = 0;
    }

    // Store the adjusted mesh for future reference.
    Store_Adjusted_Mesh();

    // Get the global refinement list.
    Get_Refinement_List(*QuadTree,
			*Local_Solution_Block_List);

    // Determine maximum interface refinement level.
    for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
      if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {
	if (Adjustment_Data[nb].Interface_Present[0]) {
	  for (int ni = 1; ni <= Adjustment_Data[nb].Ni; ni++) {
	    if (Adjustment_Data[nb].Interface_Present[ni]) {
	      maximum_interface_mesh_refinement_level[ni] =
		max(maximum_interface_mesh_refinement_level[ni],
		    Local_Solution_Block_List->Block[nb].info.level);
	      maximum_interface_mesh_refinement_flag[ni] = 
		max(maximum_interface_mesh_refinement_flag[ni],
		    Local_Solution_Block_List->RefineFlag[nb]);
	    }
	  }
	}
      }
    }

#ifdef _MPI_VERSION
    for (int n = 0; n < Interface_Union_List.Ni+1; n++) {
      maximum_interface_mesh_refinement_flag[n] = CFFC_Maximum_MPI(maximum_interface_mesh_refinement_flag[n]);
      maximum_interface_mesh_refinement_level[n] = CFFC_Maximum_MPI(maximum_interface_mesh_refinement_level[n]);
    }
#endif

    for (int number_of_passes = 0; number_of_passes < 10; number_of_passes++) {

      // Set the number of changes flag to zero at the start of each pass.
      number_of_changes = 0;

      // Impose interface refinement condition.
      for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
	if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {
	  if (Adjustment_Data[nb].Interface_Present[0]) {
	    for (int ni = 1; ni <= Adjustment_Data[nb].Ni; ni++) {
	      if (Adjustment_Data[nb].Interface_Present[ni]) {
		if (Local_Solution_Block_List->RefineFlag[nb] < maximum_interface_mesh_refinement_flag[ni]) {
		  Local_Solution_Block_List->RefineFlag[nb] = maximum_interface_mesh_refinement_flag[ni];
 		  number_of_changes++;
		}
		if (Local_Solution_Block_List->Block[nb].info.level < maximum_interface_mesh_refinement_level[ni]) {
		  Local_Solution_Block_List->RefineFlag[nb] = ADAPTIVEBLOCK2D_REFINE;
 		  number_of_changes++;
		}
	      }
	    }
	  }
	}
      }

#ifdef _MPI_VERSION
      number_of_changes = CFFC_Maximum_MPI(number_of_changes);
#endif
      if (!number_of_changes) break;

      // Get the global refinement list.
      Get_Refinement_List(*QuadTree,
			  *Local_Solution_Block_List);

    }
    if(number_of_changes) {return 121145;}

    // Unadjust the mesh.
    Mesh_Unadjustment();

  }

  //Ensure that blocks which are near an embedded boundary but not quite
  //intersected by it are not coarsened, otherwise their ghost cells
  //will become larger and they may then be intersected and would 
  //immediatly be refined during the next AMR.  This is done
  //by extrapolating an approximation to where the ghost cells 
  //would be for a coarsened block and multiplying by a safety
  //factor of 1.85 (the approximation is simply an extrapolation from
  //the diagonals, it is crude).  Really, I should look
  //for intersections with each linear element of interface splines. oh well.
  //                                 ~james
  for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
    if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      iii = Local_SolnBlk[nb].Grid.INu;
      jjj = Local_SolnBlk[nb].Grid.JNl;
      NW = Local_SolnBlk[nb].Grid.Node[iii][jjj].X
	+1.85*(Local_SolnBlk[nb].Grid.Node[iii][jjj].X-Local_SolnBlk[nb].Grid.Node[iii-4][jjj+4].X);
      
      iii = Local_SolnBlk[nb].Grid.INu;
      jjj = Local_SolnBlk[nb].Grid.JNu;
      NE = Local_SolnBlk[nb].Grid.Node[iii][jjj].X
	+1.85*(Local_SolnBlk[nb].Grid.Node[iii][jjj].X-Local_SolnBlk[nb].Grid.Node[iii-4][jjj-4].X);
      
      iii = Local_SolnBlk[nb].Grid.INl;
      jjj = Local_SolnBlk[nb].Grid.JNl;
      SW = Local_SolnBlk[nb].Grid.Node[iii][jjj].X
	+1.85*(Local_SolnBlk[nb].Grid.Node[iii][jjj].X-Local_SolnBlk[nb].Grid.Node[iii+4][jjj+4].X);
      
      iii = Local_SolnBlk[nb].Grid.INl;
      jjj = Local_SolnBlk[nb].Grid.JNu;
      SE = Local_SolnBlk[nb].Grid.Node[iii][jjj].X
	+1.85*(Local_SolnBlk[nb].Grid.Node[iii][jjj].X-Local_SolnBlk[nb].Grid.Node[iii+4][jjj-4].X);

      for (int ni = 1; ni <= Interface_Union_List.Ni; ni++) {
	// For each (linear) edge of the interface.
	intersections = 0;

	for (int np = 0; np < Interface_Union_List[ni].Spline.np-1; np++) {
	  // Store the (linear) edge of the interface.
	  Xs1 = Interface_Union_List[ni].Spline.Xp[np  ];
	  Xs2 = Interface_Union_List[ni].Spline.Xp[np+1];
	  //check for intersections
	  intersections += Line_Intersection(Xs1,Xs2,NW,NE);
	  intersections += Line_Intersection(Xs1,Xs2,NW,SW);
	  intersections += Line_Intersection(Xs1,Xs2,SW,SE);
	  intersections += Line_Intersection(Xs1,Xs2,SE,NE);
	}

	if (intersections && !Adjustment_Data[nb].Interface_Present[ni] &&
	    Local_Solution_Block_List->RefineFlag[nb] == ADAPTIVEBLOCK2D_COARSEN){
	  Local_Solution_Block_List->RefineFlag[nb] = ADAPTIVEBLOCK2D_NOCHANGE;
	}
      }

      iii = Local_SolnBlk[nb].Grid.INu;
      jjj = Local_SolnBlk[nb].Grid.JNl;
      NW = Local_SolnBlk[nb].Grid.Node[iii][jjj].X
	+1.25*(Local_SolnBlk[nb].Grid.Node[iii][jjj].X-Local_SolnBlk[nb].Grid.Node[iii-2][jjj+2].X);
      
      iii = Local_SolnBlk[nb].Grid.INu;
      jjj = Local_SolnBlk[nb].Grid.JNu;
      NE = Local_SolnBlk[nb].Grid.Node[iii][jjj].X
	+1.25*(Local_SolnBlk[nb].Grid.Node[iii][jjj].X-Local_SolnBlk[nb].Grid.Node[iii-2][jjj-2].X);
      
      iii = Local_SolnBlk[nb].Grid.INl;
      jjj = Local_SolnBlk[nb].Grid.JNl;
      SW = Local_SolnBlk[nb].Grid.Node[iii][jjj].X
	+1.25*(Local_SolnBlk[nb].Grid.Node[iii][jjj].X-Local_SolnBlk[nb].Grid.Node[iii+2][jjj+2].X);
      
      iii = Local_SolnBlk[nb].Grid.INl;
      jjj = Local_SolnBlk[nb].Grid.JNu;
      SE = Local_SolnBlk[nb].Grid.Node[iii][jjj].X
	+1.25*(Local_SolnBlk[nb].Grid.Node[iii][jjj].X-Local_SolnBlk[nb].Grid.Node[iii+2][jjj-2].X);

      for (int ni = 1; ni <= Interface_Union_List.Ni; ni++) {
	// For each (linear) edge of the interface.
	intersections = 0;

	for (int np = 0; np < Interface_Union_List[ni].Spline.np-1; np++) {
	  // Store the (linear) edge of the interface.
	  Xs1 = Interface_Union_List[ni].Spline.Xp[np  ];
	  Xs2 = Interface_Union_List[ni].Spline.Xp[np+1];
	  //check for intersections
	  intersections += Line_Intersection(Xs1,Xs2,NW,NE);
	  intersections += Line_Intersection(Xs1,Xs2,NW,SW);
	  intersections += Line_Intersection(Xs1,Xs2,SW,SE);
	  intersections += Line_Intersection(Xs1,Xs2,SE,NE);
	}

	if (intersections && !Adjustment_Data[nb].Interface_Present[ni] &&
	    Local_Solution_Block_List->Block[nb].info.level + Local_Solution_Block_List->RefineFlag[nb]
	    < maximum_interface_mesh_refinement_level[ni]) {
	  Local_Solution_Block_List->RefineFlag[nb] = ADAPTIVEBLOCK2D_REFINE;
	  //this should ensure that blocks are refined before being crossed by
	  //moving boundaries (I hope this will avoid AMR problems at boundaries,
	  //which can happen).
	  //                    ~james
	}
      }
    }
  }

  // Get global refinement list and adjust division and coarsen flags so
  // as to ensure a permissible grid block topology.
  Get_Refinement_List(*QuadTree,
		      *Local_Solution_Block_List);

  // Return the locations of the boundary nodes of the solution blocks
  // to their unmodified locations.
  Unfix_Refined_Block_Boundaries(Local_SolnBlk,
				 *Local_Solution_Block_List);

  // Coarsen solution blocks as required and return the unused solution
  // blocks to the pool of available resources for subsequent
  // refinement.
  error_flag = Coarsen_Grid(Perform_Mesh_Adjustment);
  if (error_flag) return (error_flag);

  // Refine solution blocks as required.
  error_flag = Refine_Grid(Perform_Mesh_Adjustment);
  if (error_flag) return (error_flag);

  // Renumber all solution blocks, assigning new global block numbers.
  Renumber_Solution_Blocks(*QuadTree, 
			   *Local_Solution_Block_List);

  // Determine the neighbouring blocks of all used solution blocks in
  // the quadtree data structure.
  Find_Neighbours(*QuadTree,
		  *Local_Solution_Block_List);

  // Adjust the locations of the boundary nodes of the solution blocks
  // so that the new node locations match with cell volumes of adjacent
  // solution blocks that have lower levels of mesh refinement (i.e.,
  // are coarser solution blocks).
//   Fix_Refined_Block_Boundaries(Local_SolnBlk,
//     			       *Local_Solution_Block_List);

  // Unadjust the mesh.
  Mesh_Unadjustment();   //ADDED BY JAMES FOR GHOST CELL TROUBLES WHEN MESH REFINEMENT 
                         //IS DONE ON BOUNDARIES
                         //I also changed the lines commented with &*&*&*&*&*&*

  // Re-allocate memory for all message passing buffers used to send
  // solution information between neighbouring solution blocks.
  Allocate_Message_Buffers(*Local_Solution_Block_List,
			   Local_SolnBlk[0].NumVar()+NUM_COMP_VECTOR2D);

  // Update solution information shared between neighbouring blocks.
  error_flag = Send_All_Messages(Local_SolnBlk,
				 *Local_Solution_Block_List,
				 NUM_COMP_VECTOR2D,
				 ON);
  if (!error_flag) error_flag = Send_All_Messages(Local_SolnBlk,
						  *Local_Solution_Block_List,
						  Local_SolnBlk[0].NumVar(),
						  OFF);
  if (error_flag) return (error_flag);

  // Conduct mesh adjustment on blocks in which no refinement has occured.
  if (Interface_Component_List.Ni && Perform_Mesh_Adjustment) {

    for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
      if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED){// &&         //&*&*&*&*&*&*&*&*&*&*&*&*&*&*
	//Local_Solution_Block_List->RefineFlag[nb] == ADAPTIVEBLOCK2D_NOCHANGE) {         //&*&*&*&*&*&*&*&*&*&*&*&*&*&*
	error_flag = Pre_Mesh_Adjustment_Painting(nb);
	if (!error_flag) Adjustment_Data[nb].Initialize();
	if (!error_flag) error_flag = Mesh_Adjustment_Sharp_Corners(nb);
	if (!error_flag) error_flag = Mesh_Adjustment_First(nb);
	if (!error_flag) error_flag = Mesh_Adjustment_Second(nb);
	if (!error_flag) error_flag = Mesh_Adjustment_Third(nb);
	if (error_flag) return error_flag;
      }
    }

    // Update solution information shared between neighbouring blocks.
    error_flag = Send_All_Messages(Local_SolnBlk,
				   *Local_Solution_Block_List,
				   NUM_COMP_VECTOR2D,
				   ON);
    if (error_flag) return error_flag;

    for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
      if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {
 	error_flag = Mesh_Adjustment_Finalize(nb);
	if (error_flag) return error_flag;
      }
    }

    // Update solution information shared between neighbouring blocks.
    error_flag = Send_All_Messages(Local_SolnBlk,
				   *Local_Solution_Block_List,
				   NUM_COMP_VECTOR2D,
				   ON);
    if (!error_flag) error_flag = Send_All_Messages(Local_SolnBlk,
						    *Local_Solution_Block_List,
						    Local_SolnBlk[0].NumVar(),
						    OFF);
    if (error_flag) return error_flag;

    // Conduct point-in-interface tests to determine the cell status of 
    // all unknown cells. 
    for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
      if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {
	error_flag = Post_Mesh_Adjustment_Painting(nb);
	if (error_flag) return error_flag;
      }
    }

    // Store the adjusted mesh blocks.
    Store_Adjusted_Mesh();

    // Deallocate memory for the interface refinement criteria.
    delete []maximum_interface_mesh_refinement_flag; maximum_interface_mesh_refinement_flag = NULL;
    delete []maximum_interface_mesh_refinement_level; maximum_interface_mesh_refinement_level = NULL;

  }

  // AMR procedure successfully completed.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Initial_Adaptive_Mesh_Refinement --          *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
inline int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Initial_Adaptive_Mesh_Refinement(void) {

  // Exit immediately if no initial mesh refinement is required.
  if (!IP->Number_of_Initial_Mesh_Refinements) return 0;

  int error_flag, maximum_interface_mesh_refinement_flag;
  double original_coarsenthreshold;

  // Do not allow coarsening during initial refinement of mesh.
  original_coarsenthreshold = QuadTree->CoarsenThreshold;
  QuadTree->CoarsenThreshold = ZERO;

  for (int number_of_initial_mesh_refinements = 1;
       number_of_initial_mesh_refinements <= IP->Number_of_Initial_Mesh_Refinements;
       number_of_initial_mesh_refinements++) {

    // Call the adaptive mesh refinement routine.
    error_flag = Adaptive_Mesh_Refinement(ON,ON);
    if (error_flag) return error_flag;

    // Output the refinement statistics.
    if (CFFC_Primary_MPI_Processor()) {
      cout << "\n Refinement Level #" << number_of_initial_mesh_refinements
	   << " : Number of Blocks = " << QuadTree->countUsedBlocks()
	   << ", Number of Cells = " << QuadTree->countUsedCells()
	   << ", Refinement Efficiency = " << QuadTree->efficiencyRefinement();
    }

    // Output the mesh adjustment statistics.
    if (Interface_Component_List.Ni) Mesh_Adjustment_Statistics();

    // Initial conditions.
    ICs(Local_SolnBlk,*Local_Solution_Block_List,*IP);

    // Apply boundary conditions.
    error_flag = Boundary_Conditions(ZERO);
    if (error_flag) return error_flag;

    // Send all messages.
    error_flag = Send_All_Messages(Local_SolnBlk,
				   *Local_Solution_Block_List,
				   Local_SolnBlk[0].NumVar(),
				   ON);
    if (error_flag) return error_flag;

  }
  
  // Reset coarsen threshold.
  QuadTree->CoarsenThreshold = original_coarsenthreshold;

  // Initial AMR successful.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Uniform_Adaptive_Mesh_Refinement --          *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
inline int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Uniform_Adaptive_Mesh_Refinement(void) {

  // Exit immediately if no uniform mesh refinements are required.
  if (!IP->Number_of_Uniform_Mesh_Refinements) return 0;

  int error_flag, number_of_uniform_mesh_refinements;

  // Unadjust the mesh.
  Mesh_Unadjustment();

  for (number_of_uniform_mesh_refinements = 1;
       number_of_uniform_mesh_refinements <= IP->Number_of_Uniform_Mesh_Refinements;
       number_of_uniform_mesh_refinements++) {

    // Set refinement flags to all.
    QuadTree->nochangeAll();
    //Local_Solution_Block_List->refineAll();
    for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
      Local_Solution_Block_List->RefineFlag[nb] = ADAPTIVEBLOCK2D_NOCHANGE;
      if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED &&
	  Local_Solution_Block_List->Block[nb].info.level < IP->Maximum_Refinement_Level) {
	Local_Solution_Block_List->RefineFlag[nb] = ADAPTIVEBLOCK2D_REFINE;
      }
    }

    // Call the adaptive mesh refinement routine.
    error_flag = Adaptive_Mesh_Refinement(OFF,OFF);
    if (error_flag) return error_flag;

    // Output the refinement statistics.
    if (CFFC_Primary_MPI_Processor()) {
      cout << "\n Refinement Level #" << number_of_uniform_mesh_refinements
	   << " : Number of Blocks = " << QuadTree->countUsedBlocks()
	   << ", Number of Cells = " << QuadTree->countUsedCells()
	   << ", Refinement Efficiency = " << QuadTree->efficiencyRefinement();
    }

    // Store the mesh blocks if required.
    Store_Unadjusted_Mesh();
    Store_Adjusted_Mesh();

  }

  // Perform mesh adjustment according to interface location(s).
  error_flag = Mesh_Adjustment(ON,ON);
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return error_flag;

  // Initial conditions.
  ICs(Local_SolnBlk,*Local_Solution_Block_List,*IP);

  // Apply boundary conditions.
  error_flag = Boundary_Conditions(ZERO);
  if (error_flag) return error_flag;

  // Send all messages.
  error_flag = Send_All_Messages(Local_SolnBlk,
				 *Local_Solution_Block_List,
				 Local_SolnBlk[0].NumVar(),
				 ON);
  if (error_flag) return error_flag;

  // Initial uniform AMR successful.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Boundary_Adaptive_Mesh_Refinement --         *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
inline int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Boundary_Adaptive_Mesh_Refinement(void) {

  // Exit immediately if no boundary mesh refinement is required.
  if (!IP->Number_of_Boundary_Mesh_Refinements) return 0;

  int error_flag, maximum_interface_mesh_refinement_flag;

  // Otherwise refine the mesh the number of specified times.
  for (int number_of_boundary_mesh_refinements = 1;
       number_of_boundary_mesh_refinements <= IP->Number_of_Boundary_Mesh_Refinements;
       number_of_boundary_mesh_refinements++) {

    // Set refinement flags.
    QuadTree->nochangeAll();
    for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
      Local_Solution_Block_List->RefineFlag[nb] = ADAPTIVEBLOCK2D_NOCHANGE;
      if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {
	if (Adjustment_Data[nb].Number_of_Active_Cells() > 0) {
	  for (int i = Local_SolnBlk[nb].ICl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].ICu+Local_SolnBlk[nb].Nghost; i++) {
	    if ((Local_SolnBlk[nb].Flow_Type == FLOWTYPE_INVISCID &&
		 (Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_REFLECTION ||
		  Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_REFLECTION)) ||
		(Local_SolnBlk[nb].Flow_Type != FLOWTYPE_INVISCID &&
		 (Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX ||
		  Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
		  Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX ||
		  Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL)) ||
		(Local_SolnBlk[nb].Grid.BCtypeN[i] == BC_BURNING_SURFACE ||
		 Local_SolnBlk[nb].Grid.BCtypeS[i] == BC_BURNING_SURFACE)) {
	      if (Local_Solution_Block_List->Block[nb].info.level < IP->Maximum_Refinement_Level) {
		Local_Solution_Block_List->RefineFlag[nb] = ADAPTIVEBLOCK2D_REFINE;
	      }
	    }
	  }
	  for (int j = Local_SolnBlk[nb].JCl-Local_SolnBlk[nb].Nghost; j <= Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost; j++) {
	    if ((Local_SolnBlk[nb].Flow_Type == FLOWTYPE_INVISCID &&
		 (Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_REFLECTION ||
		  Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_REFLECTION)) ||
		(Local_SolnBlk[nb].Flow_Type != FLOWTYPE_INVISCID &&
		 (Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX ||
		  Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
		  Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX ||
		  Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL)) ||
		(Local_SolnBlk[nb].Grid.BCtypeE[j] == BC_BURNING_SURFACE ||
		 Local_SolnBlk[nb].Grid.BCtypeW[j] == BC_BURNING_SURFACE)) {
	      if (Local_Solution_Block_List->Block[nb].info.level < IP->Maximum_Refinement_Level) {
		Local_Solution_Block_List->RefineFlag[nb] = ADAPTIVEBLOCK2D_REFINE;
	      }
	    }
	  }
	}
      }
    }

    // Call the adaptive mesh refinement routine.
    error_flag = Adaptive_Mesh_Refinement(OFF,ON);
    if (error_flag) return error_flag;

    // Output the refinement statistics.
    if (CFFC_Primary_MPI_Processor()) {
      cout << "\n Refinement Level #" << number_of_boundary_mesh_refinements
	   << " : Number of Blocks = " << QuadTree->countUsedBlocks()
	   << ", Number of Cells = " << QuadTree->countUsedCells()
	   << ", Refinement Efficiency = " << QuadTree->efficiencyRefinement();
    }

    // Output the mesh adjustment statistics.
    if (Interface_Component_List.Ni) Mesh_Adjustment_Statistics();

  }

  // Initial conditions.
  ICs(Local_SolnBlk,*Local_Solution_Block_List,*IP);

  // Apply boundary conditions.
  error_flag = Boundary_Conditions(ZERO);
  if (error_flag) return error_flag;

  // Send all messages.
  error_flag = Send_All_Messages(Local_SolnBlk,
				 *Local_Solution_Block_List,
				 Local_SolnBlk[0].NumVar(),
				 ON);
  if (error_flag) return error_flag;

  // Initial boundary AMR successful.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Interface_Adaptive_Mesh_Refinement --        *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
inline int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Interface_Adaptive_Mesh_Refinement(void) {

  // Exit immediately if no interface components have been specified.
  if (!Interface_Component_List.Ni) return 0;

  // Exit immediately if no interface mesh refinement is required.
  if (!IP->Number_of_Interface_Mesh_Refinements) return 0;

  int error_flag;

  for (int number_of_interface_mesh_refinements = 1;
       number_of_interface_mesh_refinements <= IP->Number_of_Interface_Mesh_Refinements;
       number_of_interface_mesh_refinements++) {

    // Set refinement flags.
    QuadTree->nochangeAll();
    for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
      Local_Solution_Block_List->RefineFlag[nb] = ADAPTIVEBLOCK2D_NOCHANGE;
      if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {
	if (Adjustment_Data[nb].Interface_Present[0] &&
	    Local_Solution_Block_List->Block[nb].info.level < IP->Maximum_Refinement_Level) {
	  Local_Solution_Block_List->RefineFlag[nb] = ADAPTIVEBLOCK2D_REFINE;
	}
      }
    }

    // Call the adaptive mesh refinement routine.
    error_flag = Adaptive_Mesh_Refinement(OFF,ON);
    if (error_flag) return error_flag;

    // Output the refinement statistics.
    if (CFFC_Primary_MPI_Processor()) {
      cout << "\n Refinement Level #"      << number_of_interface_mesh_refinements
	   << " : Number of Blocks = "     << QuadTree->countUsedBlocks()
	   << ", Number of Cells = "       << QuadTree->countUsedCells()
	   << ", Refinement Efficiency = " << QuadTree->efficiencyRefinement();
    }

    // Output the mesh adjustment statistics.
    Mesh_Adjustment_Statistics();

  }

  // Initial conditions.
  ICs(Local_SolnBlk,*Local_Solution_Block_List,*IP);

  // Apply boundary conditions.
  error_flag = Boundary_Conditions(ZERO);
  if (error_flag) return error_flag;

  // Send all messages.
  error_flag = Send_All_Messages(Local_SolnBlk,
				 *Local_Solution_Block_List,
				 Local_SolnBlk[0].NumVar(),
				 ON);
  if (error_flag) return error_flag;

  // Initial interface AMR successful.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Bounding_Box_Adaptive_Mesh_Refinement --     *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
inline int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Bounding_Box_Adaptive_Mesh_Refinement(const int &Apply_ICs) {

  // Exit immediately if no interface mesh refinement is required.
  if (!IP->Number_of_Bounding_Box_Mesh_Refinements) return 0;

  int error_flag;

  // Unadjust the mesh if required.
  if (!IP->Interface_Refinement_Condition) Mesh_Unadjustment();

  for (int number_of_bounding_box_mesh_refinements = 1;
       number_of_bounding_box_mesh_refinements <= IP->Number_of_Bounding_Box_Mesh_Refinements;
       number_of_bounding_box_mesh_refinements++) {

    // Set refinement flags.
    QuadTree->nochangeAll();
    for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
      Local_Solution_Block_List->RefineFlag[nb] = ADAPTIVEBLOCK2D_NOCHANGE;
      if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {
	for (int j = Local_SolnBlk[nb].Grid.JCl-Local_SolnBlk[nb].Nghost; j <= Local_SolnBlk[nb].Grid.JCu+Local_SolnBlk[nb].Nghost; j++) {
	  for (int i = Local_SolnBlk[nb].Grid.ICl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.ICu+Local_SolnBlk[nb].Nghost; i++) {
 	    if (Local_SolnBlk[nb].Grid.Cell[i][j].Xc.x >= IP->AMR_Xmin.x &&
 		Local_SolnBlk[nb].Grid.Cell[i][j].Xc.y >= IP->AMR_Xmin.y &&
 		Local_SolnBlk[nb].Grid.Cell[i][j].Xc.x <= IP->AMR_Xmax.x &&
 		Local_SolnBlk[nb].Grid.Cell[i][j].Xc.y <= IP->AMR_Xmax.y) {
	      Local_Solution_Block_List->RefineFlag[nb] = ADAPTIVEBLOCK2D_REFINE;
 	    }
	  }
	}
      }
    }

    // Call the adaptive mesh refinement routine.
    error_flag = Adaptive_Mesh_Refinement(OFF,IP->Interface_Refinement_Condition);
    if (error_flag) return error_flag;

    // Output the refinement statistics.
    if (CFFC_Primary_MPI_Processor() && Apply_ICs) {
      cout << "\n Refinement Level #"      << number_of_bounding_box_mesh_refinements
	   << " : Number of Blocks = "     << QuadTree->countUsedBlocks()
	   << ", Number of Cells = "       << QuadTree->countUsedCells()
	   << ", Refinement Efficiency = " << QuadTree->efficiencyRefinement();
    }

    // Store the mesh blocks if required.
    if (!IP->Interface_Refinement_Condition) {
      Store_Unadjusted_Mesh();
      Store_Adjusted_Mesh();
    }

    // Output the mesh adjustment statistics if required.
    if (IP->Interface_Refinement_Condition) Mesh_Adjustment_Statistics();

  }

  // Reset the intererface component and union lists and perform mesh
  // adjustment according to interface location(s) if requried.
  if (!IP->Interface_Refinement_Condition) {
    error_flag = Mesh_Adjustment(ON,ON);
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;
  }

  // Initial conditions if required.
  if (Apply_ICs) ICs(Local_SolnBlk,*Local_Solution_Block_List,*IP);

  // Apply boundary conditions.
  error_flag = Boundary_Conditions(ZERO);
  if (error_flag) return error_flag;

  // Send all messages.
  error_flag = Send_All_Messages(Local_SolnBlk,
				 *Local_Solution_Block_List,
				 Local_SolnBlk[0].NumVar(),
				 ON);
  if (error_flag) return error_flag;

  // Initial bounding-box AMR successful.
  return 0;

}

/**********************************************************************
 **********************************************************************
 ** Input-output routines.                                           **
 **********************************************************************
 **********************************************************************/

/**********************************************************************
 * EmbeddedBoundaries2D::Read_Restart_Files --                        *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Read_Restart_Files(const int &batch_flag,
		   int &number_of_level_set_time_steps,
		   double &levelset_Time) {

  int i, error_flag;
  char prefix[256], extension[256], restart_file_name[256];
  char *restart_file_name_ptr;
  ifstream restart_file;

  Initialize_Adjustment_Grids();

  ///////////////////////////////////////////////////////////////////////
  // Read in the additional mesh blocks and adjusted mesh data blocks. //
  ///////////////////////////////////////////////////////////////////////

  // Determine prefix of restart file names.
  i = 0;
  while (1) {
    if (IP->Restart_File_Name[i] == ' ' ||
	IP->Restart_File_Name[i] == '.') break;
    prefix[i] = IP->Restart_File_Name[i];
    i++;
    if (i > strlen(IP->Restart_File_Name)) break;
  }
  prefix[i] = '\0';
  strcat(prefix,"_adjustment_blk");

  // Read the initial data for each solution block.
  for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
    if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      // Restart file name base on global block number.
      sprintf(extension,"%.6d",Local_Solution_Block_List->Block[nb].gblknum);
      strcat(extension,".soln");
      strcpy(restart_file_name,prefix);
      strcat(restart_file_name,extension);
      restart_file_name_ptr = restart_file_name;
      // Open restart file.
      restart_file.open(restart_file_name_ptr,ios::in);
      if (restart_file.bad()) return 1;
      // Read adjustment block data.
      restart_file >> AGrid[nb];
      restart_file >> OGrid[nb];
      restart_file >> Mesh[nb];
      restart_file >> AMesh[nb];
      restart_file >> OMesh[nb];
      restart_file.setf(ios::skipws);
      restart_file >> i;
      if (i != Adjustment_Data[nb].Ni) {
	delete []Adjustment_Data[nb].Interface_Present; Adjustment_Data[nb].Interface_Present = NULL;
	Adjustment_Data[nb].Ni = i;
	Adjustment_Data[nb].Interface_Present = new int[Adjustment_Data[nb].Ni+1];
      }
      for (int ni = 0; ni < Adjustment_Data[nb].Ni+1; ni++) restart_file >> Adjustment_Data[nb].Interface_Present[ni];
      restart_file >> Adjustment_Data[nb].Number_of_Inactive_Cells;
      restart_file.unsetf(ios::skipws);
      restart_file >> Interface_Component_List;
      restart_file >> Interface_Union_List;
      // Close restart file.
      restart_file.close();
    }
  }

  //Broadcast_Interface_List(Interface_Component_List);
  //Broadcast_Interface_List(Interface_Union_List);

  /////////////////////////////////////
  // Read in the level set solution. //
  /////////////////////////////////////

  error_flag = Read_Level_Set_Restart_Solution(batch_flag,
					       number_of_level_set_time_steps,
					       levelset_Time);
  if (error_flag) return error_flag;

  // Reading of restart files complete.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Read_Level_Set_Restart_Solution --           *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Read_Level_Set_Restart_Solution(const int &batch_flag,
				int &number_of_level_set_time_steps,
				double &levelset_Time) {

  int Ni, error_flag;
  Interface2D_List Interface_List;

  // Create a 'universal' interface component list that contains
  // interface components that are to be evolved using the level set 
  // method.  Form this list on all processors, including processors
  // that do not have any active dusty2D solution blocks since they
  // may have active level set 2D solution blocks.
  Ni = 0;
  if (CFFC_Primary_MPI_Processor()) {
    // Determine the number of interfaces evolved using the level
    // set method.
    for (int ni = 1; ni <= Interface_Component_List.Ni; ni++) {
      if (Interface_Component_List[ni].Motion == INTERFACE_MOTION_LEVELSET ||
	  Interface_Component_List[ni].Motion == INTERFACE_MOTION_LEVELSET_STATIONARY ||
	  Interface_Component_List[ni].Motion == INTERFACE_MOTION_LEVELSET_EXPAND ||
	  Interface_Component_List[ni].Motion == INTERFACE_MOTION_LEVELSET_STRETCH ||
	  Interface_Component_List[ni].Motion == INTERFACE_MOTION_LEVELSET_BULKFLOW ||
	  Interface_Component_List[ni].Motion == INTERFACE_MOTION_BURNING_SURFACE) {
	Ni++;
      }
    }
  }
  // If there are no interfaces evolved by the solution of the level
  // set method then exit immediately. 
#ifdef _MPI_VERSION
  Ni = CFFC_Maximum_MPI(Ni);
#endif
  if (!Ni) return 0;

#ifdef _MPI_VERSION
  // Otherwise create universal interface list.
  if (CFFC_Primary_MPI_Processor()) {
    Interface_List.allocate(Ni);
    Ni = 1;
    for (int ni = 1; ni <= Interface_Component_List.Ni; ni++) {
      if (Interface_Component_List[ni].Motion == INTERFACE_MOTION_LEVELSET ||
	  Interface_Component_List[ni].Motion == INTERFACE_MOTION_LEVELSET_STATIONARY ||
	  Interface_Component_List[ni].Motion == INTERFACE_MOTION_LEVELSET_EXPAND ||
	  Interface_Component_List[ni].Motion == INTERFACE_MOTION_LEVELSET_STRETCH ||
	  Interface_Component_List[ni].Motion == INTERFACE_MOTION_LEVELSET_BULKFLOW ||
	  Interface_Component_List[ni].Motion == INTERFACE_MOTION_BURNING_SURFACE) {
	Interface_List[Ni].Copy(Interface_Component_List[ni]);
	Ni++;
      }
    }
  }
  Broadcast_Interface_List(Interface_List);
#endif

  // Restart the level set solution class.  Pass the interface
  // component list to the level set solution.
  LS_Local_SolnBlk = Restart_Level_Set_Solution(IP->Input_File_Name,
						batch_flag,
						error_flag,
						number_of_level_set_time_steps,
						levelset_Time,
						LS_Local_SolnBlk,
						*LS_IP,
						*LS_QuadTree,
						*LS_Global_Solution_Block_List,
						*LS_Local_Solution_Block_List,
						Interface_List);
  if (error_flag) return error_flag;

  // Interface level set solution restarted successfully.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Write_Restart_Files --                       *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Write_Restart_Files(int &number_of_time_steps,
		    int &number_of_level_set_time_steps,
		    double &Time,
		    double &levelset_Time,
		    CPUTime &processor_cpu_time) {

  int i, error_flag;
  char prefix[256], extension[256], restart_file_name[256];
  char *restart_file_name_ptr;
  ofstream restart_file;

  ////////////////////////////////////////////////////////////////
  // Write the restart files for the additional mesh blocks and //
  // adjusted mesh data blocks.                                 //
  ////////////////////////////////////////////////////////////////

  // Determine prefix of restart file names.
  i = 0;
  while (1) {
    if (IP->Restart_File_Name[i] == ' ' ||
	IP->Restart_File_Name[i] == '.') break;
    prefix[i] = IP->Restart_File_Name[i];
    i++;
    if (i > strlen(IP->Restart_File_Name)) break;
  }
  prefix[i] = '\0';
  strcat(prefix,"_adjustment_blk");

  // Write the solution data for each solution block.
  for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
    if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      // Restart file name base on global block number.
      sprintf(extension,"%.6d",Local_Solution_Block_List->Block[nb].gblknum);
      strcat(extension,".soln");
      strcpy(restart_file_name,prefix);
      strcat(restart_file_name,extension);
      restart_file_name_ptr = restart_file_name;
      // Open restart file.
      restart_file.open(restart_file_name_ptr,ios::out);
      if (restart_file.bad()) return 1;      
      // Write solution block data.
      restart_file << setprecision(14) << AGrid[nb];
      restart_file << setprecision(14) << OGrid[nb];
      restart_file << setprecision(14) << Mesh[nb];
      restart_file << setprecision(14) << AMesh[nb];
      restart_file << setprecision(14) << OMesh[nb];
      restart_file << Adjustment_Data[nb].Ni << endl;
      for (int ni = 0; ni < Adjustment_Data[nb].Ni+1; ni++)
	restart_file << Adjustment_Data[nb].Interface_Present[ni] << endl;
      restart_file << Adjustment_Data[nb].Number_of_Inactive_Cells << endl;
      restart_file << setprecision(14) << Interface_Component_List;
      restart_file << setprecision(14) << Interface_Union_List;
      // Close restart file.
      restart_file.close();
    }
  }

  ///////////////////////////////////////////////////////////////////
  // Write the level set solution block restart files if required. //
  ///////////////////////////////////////////////////////////////////

  if (Interface_Component_List.Ni && LS_Local_SolnBlk != NULL) {
    error_flag = Write_Level_Set_Restart_Solution(number_of_level_set_time_steps,
						  levelset_Time,
						  processor_cpu_time);
    if (error_flag) return error_flag;
  }

  // Writing of restart files complete.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Write_Level_Set_Restart_Solution --          *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Write_Level_Set_Restart_Solution(int &number_of_level_set_time_steps,
				 double &levelset_Time,
				 CPUTime &processor_cpu_time) {

  int Ni;

  // Determine the number of interfaces evolved using the level set method.
  Ni = 0;
  for (int ni = 1; ni <= Interface_Component_List.Ni; ni++) {
    if (Interface_Component_List[ni].Motion == INTERFACE_MOTION_LEVELSET ||
	Interface_Component_List[ni].Motion == INTERFACE_MOTION_LEVELSET_STATIONARY ||
	Interface_Component_List[ni].Motion == INTERFACE_MOTION_LEVELSET_EXPAND ||
	Interface_Component_List[ni].Motion == INTERFACE_MOTION_LEVELSET_STRETCH ||
	Interface_Component_List[ni].Motion == INTERFACE_MOTION_LEVELSET_BULKFLOW ||
	Interface_Component_List[ni].Motion == INTERFACE_MOTION_BURNING_SURFACE) {
      Ni++;
    }
  }
  // If there are no interfaces evolved by the solution of the level
  // set method then exit immediately. 
  if (!Ni) return 0;

  int error_flag;

  // Write the level set quadtree restart file.
  error_flag = Write_QuadTree(*LS_QuadTree,*LS_IP);
  if (error_flag) {
    cout << "\n ERROR: Unable to open LevelSet2D quadtree data file on processor "
	 << LS_Local_Solution_Block_List->ThisCPU << "." << endl;
  }
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return error_flag;
  // Write the level set solution block restart files.
  error_flag = Write_Restart_Solution(LS_Local_SolnBlk,
				      *LS_Local_Solution_Block_List,
				      *LS_IP,
				      number_of_level_set_time_steps,
				      levelset_Time,
				      processor_cpu_time);
  if (error_flag) {
    cout << "\n ERROR: Unable to open LevelSet2D restart output data file(s) on processor "
	 << LS_Local_Solution_Block_List->ThisCPU << "." << endl;
  }
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return error_flag;

  // Interface level set solution restart files written successfully.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Output_Elements_Tecplot --                   *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Output_Elements_Tecplot(const int &number_of_time_steps,
			const double &Time) {

  // Exit immediately if none of the solution blocks of the current
  // processor are being used.
  if (!Local_Solution_Block_List->Nused()) return 0;

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
  strcat(prefix1,"_active_cpu");
  strcat(prefix2,"_inactive_cpu");

  // Determine output data file names for this processor.
  sprintf(extension,"%.6d",Local_Solution_Block_List->ThisCPU);
  strcat(extension,".dat");
  strcpy(output_file_name1,prefix1);
  strcat(output_file_name1,extension);
  strcpy(output_file_name2,prefix2);
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
  for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
    if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      Output_Active_Elements_Tecplot(nb,
				     number_of_time_steps,
				     Time,
				     i_output_title,
				     output_file1);
      Output_Inactive_Elements_Tecplot(nb,
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

  // Solution elements data files successfully written.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Output_Active_Elements_Tecplot --            *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
void EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Output_Active_Elements_Tecplot(const int &nb,
			       const int &Number_of_Time_Steps,
			       const double &Time,
			       const int &Output_Title,
			       ostream &Out_File) {

  pState W_node;
  double yplus, ywall;
  // Ensure boundary conditions are updated before evaluating
  // solution at the nodes.
  //BCs(SolnBlk,IP);

  // Output node solution data.  
  Out_File << setprecision(14);
  if (Output_Title) {
    Out_File << "TITLE = \"" << CFFC_Name() << ": Solution on Active Elements, "
	     << "Time Step/Iteration Level = " << Number_of_Time_Steps
	     << ", Time = " << Time
	     << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
	     << "\"y\" \\ \n";
    W_node.output_labels(Out_File);
  }

  if (Adjustment_Data[nb].Number_of_Active_Cells() > 0) {

    Out_File << "ZONE T = \"Block Number = " << Local_Solution_Block_List->Block[nb].gblknum
	     << "\" \\ \n"
	     << "N = " << ((Local_SolnBlk[nb].Grid.INu - Local_SolnBlk[nb].Grid.INl + 1)*
			   (Local_SolnBlk[nb].Grid.JNu - Local_SolnBlk[nb].Grid.JNl + 1)) << " \\ \n"
	     << "E = " << Adjustment_Data[nb].Number_of_Active_Cells() << " \\ \n"
	     << "F = FEPOINT \\ \n"
	     << "ET = QUADRILATERAL \n";
    for (int j = Local_SolnBlk[nb].Grid.JNl; j <= Local_SolnBlk[nb].Grid.JNu; j++) {
      for (int i = Local_SolnBlk[nb].Grid.INl; i <= Local_SolnBlk[nb].Grid.INu; i++) {
	W_node = Wn(nb,i,j);
	yplus = Yplus(nb,i,j);
	ywall = Ywall(nb,i,j);
	Out_File << " " << Local_SolnBlk[nb].Grid.Node[i][j].X;
	W_node.output_data(Out_File, ywall, yplus);
	Out_File << endl;
      }
    }

    for (int j = Local_SolnBlk[nb].JCl; j <= Local_SolnBlk[nb].JCu; j++) {
      for (int i = Local_SolnBlk[nb].ICl; i <= Local_SolnBlk[nb].ICu; i++) {
	if (Mesh[nb].cell_status[i][j] == CELL_STATUS_ACTIVE) {
	  Out_File << (j-2)*(Local_SolnBlk[nb].Grid.INu-1) + i-1 << " ";
	  Out_File << (j-2)*(Local_SolnBlk[nb].Grid.INu-1) + i   << " ";
	  Out_File << (j-1)*(Local_SolnBlk[nb].Grid.INu-2) + i+1 + (j-2) << " ";
	  Out_File << (j-1)*(Local_SolnBlk[nb].Grid.INu-2) + i   + (j-2) << endl;
	}
      }
    }

  }

  Out_File << setprecision(6);

}

/**********************************************************************
 * EmbeddedBoundaries2D::Output_Inactive_Elements_Tecplot --         *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
void EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Output_Inactive_Elements_Tecplot(const int &nb,
				 const int &Number_of_Time_Steps,
				 const double &Time,
				 const int &Output_Title,
				 ostream &Out_File) {

  Out_File << setprecision(14);
  if (Output_Title) {
    Out_File << "TITLE = \"" << CFFC_Name() << ": Inactive Elements, "
	     << "Time Step/Iteration Level = " << Number_of_Time_Steps
	     << ", Time = " << Time
	     << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
	     << "\"y\" \\ \n";
  }

  if (Adjustment_Data[nb].Number_of_Inactive_Cells) {

    // Output node solution data.  
    Out_File << "ZONE T =  \"Block Number = " << Local_Solution_Block_List->Block[nb].gblknum
	     << "\" \\ \n"
	     << "N = " << ((Local_SolnBlk[nb].Grid.INu - Local_SolnBlk[nb].Grid.INl + 1)*
			   (Local_SolnBlk[nb].Grid.JNu - Local_SolnBlk[nb].Grid.JNl + 1)) << " \\ \n"
	     << "E = " << Adjustment_Data[nb].Number_of_Inactive_Cells << " \\ \n"
	     << "F = FEPOINT \\ \n"
	     << "ET = QUADRILATERAL \n";

    for (int j = Local_SolnBlk[nb].Grid.JNl; j <= Local_SolnBlk[nb].Grid.JNu; j++) {
      for (int i = Local_SolnBlk[nb].Grid.INl; i <= Local_SolnBlk[nb].Grid.INu; i++) {
	Out_File << " " << Local_SolnBlk[nb].Grid.Node[i][j].X << endl;
      }
    }

    for (int j = Local_SolnBlk[nb].JCl; j <= Local_SolnBlk[nb].JCu; j++) {
      for (int i = Local_SolnBlk[nb].ICl; i <= Local_SolnBlk[nb].ICu; i++) {
	if (Mesh[nb].cell_status[i][j] != CELL_STATUS_ACTIVE) {
	  Out_File << (j-2)*(Local_SolnBlk[nb].Grid.INu-1) + i-1 << " ";
	  Out_File << (j-2)*(Local_SolnBlk[nb].Grid.INu-1) + i   << " ";
	  Out_File << (j-1)*(Local_SolnBlk[nb].Grid.INu-2) + i+1 + (j-2) << " ";
	  Out_File << (j-1)*(Local_SolnBlk[nb].Grid.INu-2) + i   + (j-2) << endl;
	}
      }
    }
   
  }

  Out_File << setprecision(6);

}

/**********************************************************************
 * EmbeddedBoundaries2D::Output_Nodes_Tecplot --                      *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Output_Nodes_Tecplot(void) {

  // Exit immediately if none of the solution blocks of the current
  // processor are being used.
  if (!Local_Solution_Block_List->Nused()) return 0;

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
  strcat(prefix,"_nodes_cpu");

  // Determine output data file names for this processor.
  sprintf(extension,"%.6d",Local_Solution_Block_List->ThisCPU);
  strcat(extension,".dat");
  strcpy(output_file_name,prefix);
  strcat(output_file_name,extension);

  // Open the output data files.
  output_file_name_ptr = output_file_name;
  output_file.open(output_file_name_ptr,ios::out);
  if (output_file.bad()) return 1;

  // Write the solution data for each solution block.
  i_output_title = 1;
  for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
    if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      Output_Nodes_Tecplot(nb,i_output_title,output_file);
      if (i_output_title) i_output_title = 0;
    }
  }

  // Close the output data file.
  output_file.close();

  // Domain nodes data files successfully written.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Output_Nodes_Tecplot --                      *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
void EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Output_Nodes_Tecplot(const int &nb,
		     const int &Output_Title,
		     ostream &Out_File) {

  Out_File << setprecision(14);
  if (Output_Title) {
    Out_File << "TITLE = \"" << CFFC_Name() << ": 2D Nodes\"\n"
	     << "VARIABLES = \"x\" \\ \n"
	     << "\"y\" \\ \n"
	     << "\"xo\" \\ \n"
	     << "\"yo\" \\ \n"
	     << "\"xa\" \\ \n"
	     << "\"ya\" \\ \n";
  }

  // Output node solution data.  
  Out_File << "ZONE T = \"Block Number = " << Local_Solution_Block_List->ThisCPU
	   << "::" << Local_Solution_Block_List->Block[nb].gblknum
	   << "\" \\ \n"
	   << "I = " << Local_SolnBlk[nb].Grid.INu - Local_SolnBlk[nb].Grid.INl + 1 + 2*Local_SolnBlk[nb].Nghost << " \\ \n"
	   << "J = " << Local_SolnBlk[nb].Grid.JNu - Local_SolnBlk[nb].Grid.JNl + 1 + 2*Local_SolnBlk[nb].Nghost << " \\ \n"
	   << "F = POINT \\ \n";

  for (int j = Local_SolnBlk[nb].Grid.JNl-Local_SolnBlk[nb].Nghost; j <= Local_SolnBlk[nb].Grid.JNu+Local_SolnBlk[nb].Nghost; j++) {
    for (int i = Local_SolnBlk[nb].Grid.INl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].Grid.INu+Local_SolnBlk[nb].Nghost; i++) {
      Out_File << Local_SolnBlk[nb].Grid.Node[i][j].X
	       << OGrid[nb].Node[i][j].X
	       << AGrid[nb].Node[i][j].X
	       << endl;
    }
  }

  Out_File << setprecision(6);

}

/**********************************************************************
 * EmbeddedBoundaries2D::Output_Cell_Status_Tecplot --                *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Output_Cell_Status_Tecplot(void) {

  // Exit immediately if none of the solution blocks of the current
  // processor are being used.
  if (!Local_Solution_Block_List->Nused()) return 0;
  
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
  strcat(prefix,"_status_cpu");

  // Determine output data file names for this processor.
  sprintf(extension,"%.6d",Local_Solution_Block_List->ThisCPU);
  strcat(extension,".dat");
  strcpy(output_file_name,prefix);
  strcat(output_file_name,extension);
  
  // Open the output data files.
  output_file_name_ptr = output_file_name;
  output_file.open(output_file_name_ptr,ios::out);
  if (output_file.bad()) return 1;
  
  // Write the solution data for each solution block.
  i_output_title = 1;
  for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
    if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      Output_Cell_Status_Tecplot(nb,i_output_title,output_file);
      if (i_output_title) i_output_title = 0;
    }
  }
  
  // Close the output data file.
  output_file.close();

  // Domain nodes data files successfully written.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Output_Cell_Status_Tecplot --                *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
void EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Output_Cell_Status_Tecplot(const int &nb,
			   const int &Output_Title,
			   ostream &Out_File) {

  Out_File << setprecision(14);
  if (Output_Title) {
    Out_File << "TITLE = \"" << CFFC_Name() << ": 2D Cell Status\"\n"
	     << "VARIABLES = \"x\" \\ \n"
	     << "\"y\" \\ \n"
	     << "\"status\" \\ \n"
	     << "\"type\" \\ \n"
	     << "\"xo\" \\ \n"
	     << "\"yo\" \\ \n"
	     << "\"statuso\" \\ \n"
	     << "\"typeo\" \\ \n"
	     << "\"xa\" \\ \n"
	     << "\"ya\" \\ \n"
	     << "\"statusa\" \\ \n"
	     << "\"typea\" \\ \n";
  }

  // Output node solution data.  
  Out_File << "ZONE T = \"Block Number = " << Local_Solution_Block_List->ThisCPU
	   << "::" << Local_Solution_Block_List->Block[nb].gblknum
	   << "\" \\ \n"
	   << "I = " << Local_SolnBlk[nb].ICu - Local_SolnBlk[nb].ICl + 1 + 2*Local_SolnBlk[nb].Nghost << " \\ \n"
	   << "J = " << Local_SolnBlk[nb].JCu - Local_SolnBlk[nb].JCl + 1 + 2*Local_SolnBlk[nb].Nghost << " \\ \n"
	   << "F = POINT \\ \n";

  for (int j = Local_SolnBlk[nb].JCl-Local_SolnBlk[nb].Nghost; j <= Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost; j++) {
    for (int i = Local_SolnBlk[nb].ICl-Local_SolnBlk[nb].Nghost; i <= Local_SolnBlk[nb].ICu+Local_SolnBlk[nb].Nghost; i++) {
      Out_File << Local_SolnBlk[nb].Grid.Cell[i][j].Xc
	       << " " << Mesh[nb].cell_status[i][j]
	       << " " << Mesh[nb].cell_type[i][j]
	       << OGrid[nb].Cell[i][j].Xc
	       << " " << OMesh[nb].cell_status[i][j]
	       << " " << OMesh[nb].cell_type[i][j]
	       << AGrid[nb].Cell[i][j].Xc
	       << " " << AMesh[nb].cell_status[i][j]
	       << " " << AMesh[nb].cell_type[i][j]
	       << endl;
    }
  }

  Out_File << setprecision(6);

}

/**********************************************************************
 * EmbeddedBoundaries2D::Output_Level_Set_Tecplot --                  *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Output_Level_Set_Tecplot(const int &number_of_time_steps,
			 const double &Time) {

  // Exit immediately if none of the solution blocks of the current
  // processor are being used.
  if (!LS_Local_Solution_Block_List->Nused()) return 0;

  int error_flag;

  error_flag = Output_Tecplot(LS_Local_SolnBlk,
			      *LS_Local_Solution_Block_List,
			      *LS_IP,
			      number_of_time_steps,
			      Time);
  if (error_flag) return error_flag;

  // Level Set solution data files successfully written.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Output_Level_Set_Cells_Tecplot --            *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Output_Level_Set_Cells_Tecplot(const int &number_of_time_steps,
			       const double &Time) {

  // Exit immediately if none of the solution blocks of the current
  // processor are being used.
  if (!LS_Local_Solution_Block_List->Nused()) return 0;

  int error_flag;

  error_flag = Output_Cells_Tecplot(LS_Local_SolnBlk,
				    *LS_Local_Solution_Block_List,
				    *LS_IP,
				    number_of_time_steps,
				    Time);
  if (error_flag) return error_flag;

  // Level Set solution data files successfully written.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Output_Level_Set_Interface_Tecplot --        *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Output_Level_Set_Interface_Tecplot(const int &number_of_time_steps,
				   const double &Time) {

  // Exit immediately if none of the solution blocks of the current
  // processor are being used.
  //if (!LS_Local_Solution_Block_List->Nused()) return 0;

  int error_flag;

  // Retrieve the embedded boundaries from the level set solution.
  error_flag = Retrieve_Interface_Spline(LS_Local_SolnBlk,
					 *LS_Local_Solution_Block_List);
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return error_flag;

  // Concatenate interface splines and make the global list of embedded
  // interfaces consistent.
  error_flag = Share_Interface_Information(LS_Local_SolnBlk,
					   *LS_QuadTree,
					   *LS_Global_Solution_Block_List,
					   *LS_Local_Solution_Block_List,
					   *LS_IP);
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return error_flag;

  // Write the interface list to the output data files.
  error_flag = Output_Interface_Tecplot(LS_Local_SolnBlk,
					*LS_Local_Solution_Block_List,
					*LS_IP,
					number_of_time_steps,
					Time);
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return error_flag;

  // Level Set solution data files successfully written.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Output_Ringleb_Tecplot --                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Output_Ringleb_Tecplot(void) {

  // Exit immediately if none of the solution blocks of the current
  // processor are being used.
  if (!Local_Solution_Block_List->Nused()) return 0;

  int i, i_output_title;
  char prefix[256], extension[256], output_file_name[256];
  char *output_file_name_ptr;
  ofstream output_file;
  double l1_norm, l2_norm, max_norm, area;
  double l1_temp, l2_temp, max_temp, area_temp;
  int number_of_active_cells;

  // Determine prefix of output data file names.
  i = 0;
  while (1) {
    if (IP->Output_File_Name[i] == ' ' || IP->Output_File_Name[i] == '.') break;
    prefix[i] = IP->Output_File_Name[i];
    i++;
    if (i > strlen(IP->Output_File_Name)) break;
  }
  prefix[i] = '\0';
  strcat(prefix,"_ringleb_cpu");

  // Determine output data file names for this processor.
  sprintf(extension,"%.6d",Local_Solution_Block_List->ThisCPU);
  strcat(extension,".dat");
  strcpy(output_file_name,prefix);
  strcat(output_file_name,extension);
  
  // Open the output data files.
  output_file_name_ptr = output_file_name;
  output_file.open(output_file_name_ptr,ios::out);
  if (output_file.bad()) return 1;
  
  // Write the solution data for each solution block.
  i_output_title = 1;
  l1_norm = ZERO; l2_norm = ZERO; max_norm = ZERO; area = ZERO;
  for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
    if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      l1_temp = ZERO; l2_temp = ZERO; max_temp = ZERO; area_temp = ZERO;
      Output_Ringleb_Tecplot(nb,i_output_title,output_file,
			     l1_temp,l2_temp,max_temp,area_temp);
      l1_norm += l1_temp;
      l2_norm += l2_temp;
      max_norm = max(max_norm,max_temp);
      area += area_temp;
      number_of_active_cells += Adjustment_Data[nb].Number_of_Active_Cells();
      if (i_output_title) i_output_title = 0;
    }
  }
  
  // Close the output data file.
  output_file.close();

#ifdef _MPI_VERSION
  l1_norm = CFFC_Summation_MPI(l1_norm);
  l2_norm = CFFC_Summation_MPI(l2_norm);
  max_norm = CFFC_Maximum_MPI(max_norm);
  area = CFFC_Summation_MPI(area);
  number_of_active_cells = CFFC_Summation_MPI(number_of_active_cells);
#endif

  // Calculate the L1-norm and L2-norm for all blocks.
  l1_norm /= area;
  l2_norm = sqrt(l2_norm/area);

  if (CFFC_Primary_MPI_Processor()) {
    cout << endl
	 << endl
	 << " ==================================================================== "
	 << endl
	 << " Error norms for Ringleb's flow:"
	 << endl
	 << "   L1_Norm = " << l1_norm
	 << endl
	 << "   L2_Norm = " << l2_norm
	 << endl
	 << "   Max_Norm = " << max_norm;
    if (number_of_active_cells) cout << endl << "   Number of active cells = " << number_of_active_cells;
    cout << endl
	 << " ==================================================================== "
	 << endl;
  }

  // Ringleb's flow data files successfully written.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Output_Ringleb_Tecplot --                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
void EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Output_Ringleb_Tecplot(const int &nb,
		       const int &Output_Title,
		       ostream &Out_File,
		       double &l1_norm,
		       double &l2_norm,
		       double &max_norm,
		       double &area) {

  pState Wn, We;

  // Output node solution data.  
  Out_File << setprecision(14);
  if (Output_Title)
    Out_File << "TITLE = \"" << CFFC_Name() << ": 2D Ringleb Flow Solution "
	     << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
	     << "\"y\" \\ \n"
	     << "\"rho\" \\ \n"
	     << "\"a\" \\ \n"
	     << "\"rho_exact\" \\ \n"
	     << "\"a_exact\" \\ \n"
	     << "\"rho_error\" \\ \n"
	     << "\"a_error\" \\ \n";
  Out_File << "ZONE T =  \"Block Number = " << nb << "\" \\ \n"
	   << "I = " << Local_SolnBlk[nb].ICu - Local_SolnBlk[nb].ICl + 1 << " \\ \n"
	   << "J = " << Local_SolnBlk[nb].JCu - Local_SolnBlk[nb].JCl + 1 << " \\ \n"
	   << "F = POINT \\ \n";
  Out_File.setf(ios::scientific);
  for (int j = Local_SolnBlk[nb].JCl; j <= Local_SolnBlk[nb].JCu; j++) {
    for (int i = Local_SolnBlk[nb].ICl; i <= Local_SolnBlk[nb].ICu; i++) {
      if (Mesh[nb].cell_status[i][j] == CELL_STATUS_ACTIVE) {
	We = RinglebFlowAverageState(We,Local_SolnBlk[nb].Grid.nodeSW(i,j).X,
				     Local_SolnBlk[nb].Grid.nodeSE(i,j).X,
				     Local_SolnBlk[nb].Grid.nodeNE(i,j).X,
				     Local_SolnBlk[nb].Grid.nodeNW(i,j).X);
	//We = RinglebFlow(We,Local_SolnBlk[nb].Grid.Cell[i][j].Xc);
	l1_norm += fabs(We[1] - Local_SolnBlk[nb].W[i][j][1])*Local_SolnBlk[nb].Grid.Cell[i][j].A;
	l2_norm += sqr(We[1] - Local_SolnBlk[nb].W[i][j][1])*Local_SolnBlk[nb].Grid.Cell[i][j].A;
	max_norm = max(max_norm,fabs(We[1] - Local_SolnBlk[nb].W[i][j][1]));
	area += Local_SolnBlk[nb].Grid.Cell[i][j].A;
      } else {
	We = Local_SolnBlk[nb].W[i][j];
      }
      Out_File << Local_SolnBlk[nb].Grid.Cell[i][j].Xc << " " 
	       << Local_SolnBlk[nb].W[i][j][1] << " " 
	       << Local_SolnBlk[nb].W[i][j].a() << " " 
	       << We[1] << " " 
	       << We.a()   << " " 
	       << fabs(Local_SolnBlk[nb].W[i][j][1] - We[1]) << " " 
	       << fabs(Local_SolnBlk[nb].W[i][j].a() - We.a()) << endl;
    }
  }
  Out_File.unsetf(ios::scientific);

}

/**********************************************************************
 * EmbeddedBoundaries2D::Output_Flat_Plate_Tecplot --                 *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Output_Flat_Plate_Tecplot(void) {

  // Exit immediately if none of the solution blocks of the current
  // processor are being used.
  if (!Local_Solution_Block_List->Nused()) return 0;
  
  int i, i_output_title_soln, i_output_title_skin;
  char prefix_soln[256], prefix_skin[256], extension_soln[256], extension_skin[256];
  char output_file_name_soln[256], output_file_name_skin[256];
  char *output_file_name_soln_ptr, *output_file_name_skin_ptr;
  ofstream output_file_soln, output_file_skin;
  ofstream output_file;
  double l1_norm, l2_norm, max_norm, area;
  double l1_temp, l2_temp, max_temp, area_temp;
  double l1_norm_cf, l2_norm_cf, max_norm_cf, area_cf;
  double l1_cf_temp, l2_cf_temp, max_cf_temp, area_cf_temp;
  int numberofcells, numberofcells_temp;
  int numberofcells_cf, numberofcells_cf_temp;

  // Determine prefix of output data file names.
  i = 0;
  while (1) {
    if (IP->Output_File_Name[i] == ' ' || IP->Output_File_Name[i] == '.') break;
    prefix_soln[i] = IP->Output_File_Name[i];
    i++;
    if (i > strlen(IP->Output_File_Name)) break;
  }
  prefix_soln[i] = '\0';
  strcat(prefix_soln,"_flatplate_soln_cpu");
  i = 0;
  while (1) {
    if (IP->Output_File_Name[i] == ' ' || IP->Output_File_Name[i] == '.') break;
    prefix_skin[i] = IP->Output_File_Name[i];
    i++;
    if (i > strlen(IP->Output_File_Name)) break;
  }
  prefix_skin[i] = '\0';
  strcat(prefix_skin,"_flatplate_skin_friction_cpu");
  
  // Determine output data file name for this processor.
  sprintf(extension_soln,"%.6d",Local_Solution_Block_List->ThisCPU);
  strcat(extension_soln,".dat");
  strcpy(output_file_name_soln,prefix_soln);
  strcat(output_file_name_soln,extension_soln);
  output_file_name_soln_ptr = output_file_name_soln;
  sprintf(extension_skin,"%.6d",Local_Solution_Block_List->ThisCPU);
  strcat(extension_skin,".dat");
  strcpy(output_file_name_skin,prefix_skin);
  strcat(output_file_name_skin,extension_skin);
  output_file_name_skin_ptr = output_file_name_skin;
  
  // Open the output data files.
  output_file_soln.open(output_file_name_soln_ptr,ios::out);
  if (output_file_soln.bad()) return 1;
  output_file_skin.open(output_file_name_skin_ptr,ios::out);
  if (output_file_skin.bad()) return 1;
  
  // Write the solution data for each solution block.
  i_output_title_soln = 1;
  i_output_title_skin = 1;
  l1_norm = ZERO; l2_norm = ZERO; max_norm = ZERO; area = ZERO;
  l1_norm_cf = ZERO; l2_norm_cf = ZERO; max_norm_cf = ZERO; area_cf = ZERO;
  numberofcells = 0; numberofcells_cf = 0;
  for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
    if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      l1_temp = ZERO; l2_temp = ZERO; max_temp = ZERO; area_temp = 0; numberofcells_temp = 0;
      l1_cf_temp = ZERO; l2_cf_temp = ZERO; max_cf_temp = ZERO; area_cf_temp = 0; numberofcells_cf_temp = 0;
      Output_Flat_Plate_Tecplot(nb,
				i_output_title_soln,output_file_soln,
				i_output_title_skin,output_file_skin,
				l1_temp,l2_temp,max_temp,area_temp,numberofcells_temp,
				l1_cf_temp,l2_cf_temp,max_cf_temp,area_cf_temp,numberofcells_cf_temp);
      l1_norm += l1_temp;
      l2_norm += l2_temp;
      max_norm = max(max_norm,max_temp);
      area += area_temp;
      numberofcells += numberofcells_temp;
      l1_norm_cf += l1_cf_temp;
      l2_norm_cf += l2_cf_temp;
      max_norm_cf = max(max_norm_cf,max_cf_temp);
      area_cf += area_cf_temp;
      numberofcells_cf += numberofcells_cf_temp;
      if (i_output_title_soln) i_output_title_soln = 0;
      if (i_output_title_skin) i_output_title_skin = 0;
    }
  }

  // Close the output data file.
  output_file_soln.close();
  output_file_skin.close();

  // Calculate the L1-norm and L2-norm for all blocks.
#ifdef _MPI_VERSION
  l1_norm = CFFC_Summation_MPI(l1_norm);
  l2_norm = CFFC_Summation_MPI(l2_norm);
  max_norm = CFFC_Maximum_MPI(max_norm);
  area = CFFC_Summation_MPI(area);
  numberofcells = CFFC_Summation_MPI(numberofcells);
  l1_norm_cf = CFFC_Summation_MPI(l1_norm_cf);
  l2_norm_cf = CFFC_Summation_MPI(l2_norm_cf);
  max_norm_cf = CFFC_Maximum_MPI(max_norm_cf);
  area_cf = CFFC_Summation_MPI(area_cf);
  numberofcells_cf = CFFC_Summation_MPI(numberofcells_cf);
#endif
  l1_norm /= area;
  l2_norm = sqrt(l2_norm/area);
  l1_norm_cf /= area_cf;
  l2_norm_cf = sqrt(l2_norm_cf/area_cf);

  if (CFFC_Primary_MPI_Processor()) {
    cout << endl
	 << endl
	 << " ==================================================================== "
	 << endl
	 << " Error norms for the u-velocity component:"
	 << endl
	 << "   L1_Norm = " << l1_norm
	 << endl
	 << "   L2_Norm = " << l2_norm
	 << endl
	 << "   Max_Norm = " << max_norm
	 << endl
	 << "   Number of cells = " << numberofcells
	 << endl
	 << endl
	 << " Error norms for the skin friction coefficient for the laminar flow"
	 << " over a flat plate:"
	 << endl
	 << "   L1_Norm = " << l1_norm_cf
	 << endl
	 << "   L2_Norm = " << l2_norm_cf
	 << endl
	 << "   Max_Norm = " << max_norm_cf
	 << endl
	 << "   Number of cells = " << numberofcells_cf
	 << endl
	 << " ==================================================================== "
	 << endl;
  }

  // Flat plate flow data files successfully written.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Output_Flat_Plate_Tecplot --                 *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
void EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Output_Flat_Plate_Tecplot(const int &nb,
			  const int &Output_Title_Soln,
			  ostream &Out_File_Soln,
			  const int &Output_Title_Skin,
			  ostream &Out_File_Skin,
			  double &l1_norm,
			  double &l2_norm,
			  double &max_norm,
			  double &area,
			  int &numberofcells,
			  double &l1_norm_cf,
			  double &l2_norm_cf,
			  double &max_norm_cf,
			  double &area_cf,
			  int &numberofcells_cf) {
  cout << endl << " ERROR: Explicit specialization required for the generation of the ";
  cout << endl << "        flat-plate output data file for a specific solution block.";
  cout.flush();
  assert(1==0);
}

/**********************************************************************
 * EmbeddedBoundaries2D::Output_Aerodynamic_Coefficients_Tecplot --   *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Output_Aerodynamic_Coefficients_Tecplot(const int &number_of_time_steps,
					const double &Time) {

  // Exit immediately if none of the solution blocks of the current
  // processor are being used.
  if (!Local_Solution_Block_List->Nused()) return 0;

  int i, i_output_title;
  char prefix[256], extension[256], output_file_name[256];
  char *output_file_name_ptr;
  ofstream output_file;
  double Cn, Ca, Cm, Cl, Cd, Cn_temp, Ca_temp, Cm_temp;
  double Cnp, Cmp, Cnp_temp, Cmp_temp;

  // Determine prefix of output data file names.
  i = 0;
  while (1) {
    if (IP->Output_File_Name[i] == ' ' || IP->Output_File_Name[i] == '.') break;
    prefix[i] = IP->Output_File_Name[i];
    i++;
    if (i > strlen(IP->Output_File_Name)) break;
  }
  prefix[i] = '\0';
  strcat(prefix,"_cp_cpu");

  // Determine output data file names for this processor.
  sprintf(extension,"%.6d",Local_Solution_Block_List->ThisCPU);
  strcat(extension,".dat");
  strcpy(output_file_name,prefix);
  strcat(output_file_name,extension);
  
  // Open the output data files.
  output_file_name_ptr = output_file_name;
  output_file.open(output_file_name_ptr,ios::out);
  if (output_file.bad()) return 1;
  
  // Write the solution data for each solution block.
  i_output_title = 1;
  Cn = ZERO; Ca = ZERO; Cm = ZERO;
  Cnp = ZERO; Cmp = ZERO;
  for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
    if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      Cn_temp = ZERO; Ca_temp = ZERO; Cm_temp = ZERO;
      Cnp_temp = ZERO; Cmp_temp = ZERO;
      Output_Aerodynamic_Coefficients_Tecplot(nb,i_output_title,output_file,
					      Cn_temp,Ca_temp,Cm_temp,
					      Cnp_temp,Cmp_temp);
      Cn += Cn_temp; Ca += Ca_temp; Cm += Cm_temp;
      Cnp += Cnp_temp; Cmp += Cmp_temp;
      if (i_output_title) i_output_title = 0;
    }
  }
  
  // Close the output data file.
  output_file.close();

  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

#ifdef _MPI_VERSION
  Cn = CFFC_Summation_MPI(Cn);
  Ca = CFFC_Summation_MPI(Ca);
  Cm = CFFC_Summation_MPI(Cm);
  Cnp = CFFC_Summation_MPI(Cnp);
  Cmp = CFFC_Summation_MPI(Cmp);
#endif

  double cos_alpha, sin_alpha;
  Vector2D nhat, ahat;

  // Determine the axial and normal directions of the aerofoil:
  ahat = (Interface_Union_List[1].Spline.Xp[0] - Interface_Union_List[1].Spline.Xp[(Interface_Union_List[1].Spline.np-1)/2])/
         abs(Interface_Union_List[1].Spline.Xp[0] - Interface_Union_List[1].Spline.Xp[(Interface_Union_List[1].Spline.np-1)/2]);
  nhat = Vector2D(-ahat.y,ahat.x);

  // Determine the cos and sin of the angle of attack of the aerofoil:
  cos_alpha = dot(nhat,jhat);
  sin_alpha = dot(nhat,ihat);

  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

  if (CFFC_Primary_MPI_Processor()) {
    cout << endl
	 << endl
	 << " ==================================================================== "
	 << endl << " Aerodynamic coefficients:"
	 << endl << "   Cn = " << Cn
	 << endl << "   Ca = " << Ca
	 << endl << "   Cl = " << Cn*cos_alpha - Ca*sin_alpha
	 << endl << "   Cd = " << Cn*sin_alpha + Ca*cos_alpha
	 << endl << "   Cm = " << Cm
	 << endl << "   Cnp = " << Cnp
	 << endl << "   Cmp = " << Cmp
	 << endl
	 << " ==================================================================== "
	 << endl;
  }

  // Aerodynamic coefficients successfully calculated and written.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Output_Aerodynamic_Coefficients_Tecplot --   *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
void EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Output_Aerodynamic_Coefficients_Tecplot(const int &nb,
					const int &Output_Title,
					ostream &Out_File,
					double &Cn, double &Ca, double &Cm,
					double &Cnp, double &Cmp) {

  int npts;

  npts = 0;
  for (int j = Local_SolnBlk[nb].JCl; j <= Local_SolnBlk[nb].JCu; j++) {
    for (int i = Local_SolnBlk[nb].ICl; i <= Local_SolnBlk[nb].ICu; i++) {
      if (Mesh[nb].cell_status[i][j] != CELL_STATUS_ACTIVE) {
	if (Local_SolnBlk[nb].Grid.lfaceN(i,j) && Mesh[nb].cell_status[i][j+1] == CELL_STATUS_ACTIVE) {
	  npts++;
	}
	if (Local_SolnBlk[nb].Grid.lfaceS(i,j) && Mesh[nb].cell_status[i][j-1] == CELL_STATUS_ACTIVE) {
	  npts++;
	}
	if (Local_SolnBlk[nb].Grid.lfaceE(i,j) && Mesh[nb].cell_status[i+1][j] == CELL_STATUS_ACTIVE) {
	  npts++;
	}
	if (Local_SolnBlk[nb].Grid.lfaceW(i,j) && Mesh[nb].cell_status[i-1][j] == CELL_STATUS_ACTIVE) {
	  npts++;
	}
      }
    }
  }

  if (Output_Title)
    Out_File << "TITLE = \"" << CFFC_Name() << ": 2D Aerodynamic Coefficients "
	     << "\"" << "\n"
	     << "VARIABLES = \"x/c\" \\ \n"
	     << "\"y/c\" \\ \n"
	     << "\"Cp\" \\ \n"
	     << "\"Cf\" \\ \n"
	     << "\"Cn\" \\ \n"
	     << "\"Ca\" \\ \n"
	     << "\"Cl\" \\ \n"
	     << "\"Cd\" \\ \n"
	     << "\"Cm\" \\ \n";

  if (!npts) return ;

  double pinf, qinf, Cpi, Cfi, Cni, Cai, Cli, Cdi, Cmi, c, x, y, cos_alpha, sin_alpha, ds, s, s1, s2;
  Vector2D Xc, nhat, ahat;

  // Determine the location of the quarter chord of the aerofoil:
  Xc = Interface_Union_List[1].Xref;

  // Determine the chord length of the aerofoil:
  c = Interface_Union_List[1].Length1;

  // Determine the axial and normal directions of the aerofoil:
  ahat = (Interface_Union_List[1].Spline.Xp[0] - Interface_Union_List[1].Spline.Xp[(Interface_Union_List[1].Spline.np-1)/2])/
         abs(Interface_Union_List[1].Spline.Xp[0] - Interface_Union_List[1].Spline.Xp[(Interface_Union_List[1].Spline.np-1)/2]);
  nhat = Vector2D(-ahat.y,ahat.x);

  // Determine the cos and sin of the angle of attack of the aerofoil:
  cos_alpha = dot(nhat,jhat);
  sin_alpha = dot(nhat,ihat);

  // Determine the static and dynamic pressures:
  pinf = IP->Wo.pressure();
  qinf = HALF*IP->Wo[1]*sqr(IP->Wo.v.abs());

  // Output node solution data.  
  Out_File << setprecision(14);
  Out_File << "ZONE T = \"Block Number = " << Local_Solution_Block_List->Block[nb].gblknum << "\" \\ \n"
	   << "I = " << npts << " \\ \n"
	   << "J = " << 1 << " \\ \n"
	   << "F = POINT \\ \n";
  Out_File.setf(ios::scientific);
  for (int j = Local_SolnBlk[nb].JCl; j <= Local_SolnBlk[nb].JCu; j++) {
    for (int i = Local_SolnBlk[nb].ICl; i <= Local_SolnBlk[nb].ICu; i++) {
      if (Mesh[nb].cell_status[i][j] != CELL_STATUS_ACTIVE) {
	if (Local_SolnBlk[nb].Grid.lfaceN(i,j) > TOLER && Mesh[nb].cell_status[i][j+1] == CELL_STATUS_ACTIVE) {
	  x = dot((Local_SolnBlk[nb].Grid.xfaceN(i,j)-Xc),ahat) + 0.25;
	  y = dot((Local_SolnBlk[nb].Grid.xfaceN(i,j)-Xc),nhat);
	  s = getS(Local_SolnBlk[nb].Grid.xfaceN(i,j),Interface_Union_List[1].Spline);
	  s1 = getS(Local_SolnBlk[nb].Grid.nodeNE(i,j).X,Interface_Union_List[1].Spline);
	  s2 = getS(Local_SolnBlk[nb].Grid.nodeNW(i,j).X,Interface_Union_List[1].Spline);
	  if (s2 < TOLER && s1 > HALF*Interface_Union_List[1].Spline.sp[Interface_Union_List[1].Spline.np-1]) {
	    s2 = Interface_Union_List[1].Spline.sp[Interface_Union_List[1].Spline.np-1];
	  }
	  ds = fabs(s2-s1);
	  Cpi = (Local_SolnBlk[nb].W[i][j+1].pressure() - pinf)/qinf;
	  if (Local_SolnBlk[nb].Flow_Type == FLOWTYPE_INVISCID) Cfi = ZERO;
	  else Cfi = ZERO;
	  Cni = -Cpi*dot(Local_SolnBlk[nb].Grid.nfaceN(i,j),nhat)*ds/c; Cn += Cni;
	  Cai = -Cpi*dot(Local_SolnBlk[nb].Grid.nfaceN(i,j),ahat)*ds/c; Ca += Cai;
	  Cmi =  Cpi*dot(Local_SolnBlk[nb].Grid.nfaceN(i,j),nhat)*(0.25 - x)*ds/c
	        +Cpi*dot(Local_SolnBlk[nb].Grid.nfaceN(i,j),ahat)*y*ds/c;
	  Cm += Cmi;
	  if (s < HALF*Interface_Union_List[1].Spline.sp[Interface_Union_List[1].Spline.np-1]) {
	    Cnp -= Cpi*ds/c;
	    Cmp -= Cpi*(0.25 - x)*ds/c;
	  } else {
	    Cnp += Cpi*ds/c;
	    Cmp += Cpi*(0.25 - x)*ds/c;
	  }
	  Out_File << " " << x
		   << " " << y
		   << " " << Cpi
		   << " " << Cfi
		   << " " << Cni
		   << " " << Cai
		   << " " << Cni*cos_alpha - Cai*sin_alpha
		   << " " << Cni*sin_alpha + Cai*cos_alpha
		   << " " << Cmi
		   << endl;
	}
	if (Local_SolnBlk[nb].Grid.lfaceS(i,j) > TOLER && Mesh[nb].cell_status[i][j-1] == CELL_STATUS_ACTIVE) {
	  x = dot((Local_SolnBlk[nb].Grid.xfaceS(i,j)-Xc),ahat) + 0.25;
	  y = dot((Local_SolnBlk[nb].Grid.xfaceS(i,j)-Xc),nhat);
	  s = getS(Local_SolnBlk[nb].Grid.xfaceS(i,j),Interface_Union_List[1].Spline);
	  s1 = getS(Local_SolnBlk[nb].Grid.nodeSW(i,j).X,Interface_Union_List[1].Spline);
	  s2 = getS(Local_SolnBlk[nb].Grid.nodeSE(i,j).X,Interface_Union_List[1].Spline);
	  if (s2 < TOLER && s1 > HALF*Interface_Union_List[1].Spline.sp[Interface_Union_List[1].Spline.np-1]) {
	    s2 = Interface_Union_List[1].Spline.sp[Interface_Union_List[1].Spline.np-1];
	  }
	  ds = fabs(s2-s1);
	  Cpi = (Local_SolnBlk[nb].W[i][j-1].pressure() - pinf)/qinf;
	  if (Local_SolnBlk[nb].Flow_Type == FLOWTYPE_INVISCID) Cfi = ZERO;
	  else Cfi = ZERO;
	  Cni = -Cpi*dot(Local_SolnBlk[nb].Grid.nfaceS(i,j),nhat)*ds/c; Cn += Cni;
	  Cai = -Cpi*dot(Local_SolnBlk[nb].Grid.nfaceS(i,j),ahat)*ds/c; Ca += Cai;
	  Cmi =  Cpi*dot(Local_SolnBlk[nb].Grid.nfaceS(i,j),nhat)*(0.25 - x)*ds/c
	        +Cpi*dot(Local_SolnBlk[nb].Grid.nfaceS(i,j),ahat)*y*ds/c;
	  Cm += Cmi;
	  Out_File << " " << x
		   << " " << y
		   << " " << Cpi
		   << " " << Cfi
		   << " " << Cni
		   << " " << Cai
		   << " " << Cni*cos_alpha - Cai*sin_alpha
		   << " " << Cni*sin_alpha + Cai*cos_alpha
		   << " " << Cmi
		   << endl;
	  if (s < HALF*Interface_Union_List[1].Spline.sp[Interface_Union_List[1].Spline.np-1]) {
	    Cnp -= Cpi*ds/c;
	    Cmp -= Cpi*(0.25 - x)*ds/c;
	  } else {
	    Cnp += Cpi*ds/c;
	    Cmp += Cpi*(0.25 - x)*ds/c;
	  }
	}
	if (Local_SolnBlk[nb].Grid.lfaceE(i,j) > TOLER && Mesh[nb].cell_status[i+1][j] == CELL_STATUS_ACTIVE) {
	  x = dot((Local_SolnBlk[nb].Grid.xfaceE(i,j)-Xc),ahat) + 0.25;
	  y = dot((Local_SolnBlk[nb].Grid.xfaceE(i,j)-Xc),nhat);
	  s = getS(Local_SolnBlk[nb].Grid.xfaceE(i,j),Interface_Union_List[1].Spline);
	  s1 = getS(Local_SolnBlk[nb].Grid.nodeSE(i,j).X,Interface_Union_List[1].Spline);
	  s2 = getS(Local_SolnBlk[nb].Grid.nodeNE(i,j).X,Interface_Union_List[1].Spline);
	  if (s2 < TOLER && s1 > HALF*Interface_Union_List[1].Spline.sp[Interface_Union_List[1].Spline.np-1]) {
	    s2 = Interface_Union_List[1].Spline.sp[Interface_Union_List[1].Spline.np-1];
	  }
	  ds = fabs(s2-s1);
	  Cpi = (Local_SolnBlk[nb].W[i+1][j].pressure() - pinf)/qinf;
	  if (Local_SolnBlk[nb].Flow_Type == FLOWTYPE_INVISCID) Cfi = ZERO;
	  else Cfi = ZERO;
	  Cni = -Cpi*dot(Local_SolnBlk[nb].Grid.nfaceE(i,j),nhat)*ds/c; Cn += Cni;
	  Cai = -Cpi*dot(Local_SolnBlk[nb].Grid.nfaceE(i,j),ahat)*ds/c; Ca += Cai;
	  Cmi =  Cpi*dot(Local_SolnBlk[nb].Grid.nfaceE(i,j),nhat)*(0.25 - x)*ds/c
	        +Cpi*dot(Local_SolnBlk[nb].Grid.nfaceE(i,j),ahat)*y*ds/c;
	  Cm += Cmi;
	  Out_File << " " << x
		   << " " << y
		   << " " << Cpi
		   << " " << Cfi
		   << " " << Cni
		   << " " << Cai
		   << " " << Cni*cos_alpha - Cai*sin_alpha
		   << " " << Cni*sin_alpha + Cai*cos_alpha
		   << " " << Cmi
		   << endl;
	  if (s < HALF*Interface_Union_List[1].Spline.sp[Interface_Union_List[1].Spline.np-1]) {
	    Cnp -= Cpi*ds/c;
	    Cmp -= Cpi*(0.25 - x)*ds/c;
	  } else {
	    Cnp += Cpi*ds/c;
	    Cmp += Cpi*(0.25 - x)*ds/c;
	  }
	}
	if (Local_SolnBlk[nb].Grid.lfaceW(i,j) > TOLER && Mesh[nb].cell_status[i-1][j] == CELL_STATUS_ACTIVE) {
	  x = dot((Local_SolnBlk[nb].Grid.xfaceW(i,j)-Xc),ahat) + 0.25;
	  y = dot((Local_SolnBlk[nb].Grid.xfaceW(i,j)-Xc),nhat);
	  s = getS(Local_SolnBlk[nb].Grid.xfaceW(i,j),Interface_Union_List[1].Spline);
	  s1 = getS(Local_SolnBlk[nb].Grid.nodeNW(i,j).X,Interface_Union_List[1].Spline);
	  s2 = getS(Local_SolnBlk[nb].Grid.nodeSW(i,j).X,Interface_Union_List[1].Spline);
	  if (s2 < TOLER && s1 > HALF*Interface_Union_List[1].Spline.sp[Interface_Union_List[1].Spline.np-1]) {
	    s2 = Interface_Union_List[1].Spline.sp[Interface_Union_List[1].Spline.np-1];
	  }
	  ds = fabs(s2-s1);
	  Cpi = (Local_SolnBlk[nb].W[i-1][j].pressure() - pinf)/qinf;
	  if (Local_SolnBlk[nb].Flow_Type == FLOWTYPE_INVISCID) Cfi = ZERO;
	  else Cfi = ZERO;
	  Cni = -Cpi*dot(Local_SolnBlk[nb].Grid.nfaceW(i,j),nhat)*ds/c; Cn += Cni;
	  Cai = -Cpi*dot(Local_SolnBlk[nb].Grid.nfaceW(i,j),ahat)*ds/c; Ca += Cai;
	  Cmi =  Cpi*dot(Local_SolnBlk[nb].Grid.nfaceW(i,j),nhat)*(0.25 - x)*ds/c
	        +Cpi*dot(Local_SolnBlk[nb].Grid.nfaceW(i,j),ahat)*y*ds/c;
	  Cm += Cmi;
	  Out_File << " " << x
		   << " " << y
		   << " " << Cpi
		   << " " << Cfi
		   << " " << Cni
		   << " " << Cai
		   << " " << Cni*cos_alpha - Cai*sin_alpha
		   << " " << Cni*sin_alpha + Cai*cos_alpha 
		   << " " << Cmi
		   << endl;
	  if (s < HALF*Interface_Union_List[1].Spline.sp[Interface_Union_List[1].Spline.np-1]) {
	    Cnp -= Cpi*ds/c;
	    Cmp -= Cpi*(0.25 - x)*ds/c;
	  } else {
	    Cnp += Cpi*ds/c;
	    Cmp += Cpi*(0.25 - x)*ds/c;
	  }
	}
      }
    }
  }
  Out_File.unsetf(ios::scientific);

}

/**********************************************************************
 * EmbeddedBoundaries2D::Net_Force --                                 *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
Vector2D EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Net_Force(void) {
  cout << endl << " ERROR: Explicit specialization required for the calculation of the";
  cout << endl << "        net force.";
  cout.flush();
  return 1;
}

/**********************************************************************
 * EmbeddedBoundaries2D::Net_Pressure_Force --                        *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
Vector2D EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Net_Pressure_Force(void) {
  cout << endl << " ERROR: Explicit specialization required for the calculation of the";
  cout << endl << "        net pressure force.";
  cout.flush();
  return 1;
}

/**********************************************************************
 * EmbeddedBoundaries2D::Output_Cylinder_Drag --                      *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Output_Cylinder_Drag(void) {
  cout << endl << " ERROR: Explicit specialization required for the calculation of Cl";
  cout << endl << "        and Cd for the Cylinder Case.";
  cout.flush();
  return 1;
}

/**********************************************************************
 * EmbeddedBoundaries2D::Output_Couette       --                      *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Output_Couette(void) {
  cout << endl << " ERROR: Explicit specialization required for Ouput_Couette.";
  cout.flush();
  return 1;
}

/**********************************************************************
 * EmbeddedBoundaries2D::Ywall       --                               *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
double EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Ywall(const int &nb, const int &i, const int &j) {
  //This is used by everything that does not specialize this function
  return 0.0;
}

/**********************************************************************
 * EmbeddedBoundaries2D::plus        --                               *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
double EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Yplus(const int &nb, const int &i, const int &j) {
  //This is used by everything that does not specialize this function
  return 0.0;
}



/**********************************************************************
 **********************************************************************
 ** Solution routines.                                               **
 **********************************************************************
 **********************************************************************/

/**********************************************************************
 * EmbeddedBoundaries2D::Boundary_Conditions --                       *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
inline int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Boundary_Conditions(const double &Time) {

  int error_flag;

  // Prescribe boundary data for each solution block.
  for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
    if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      // Apply embedded boundary conditions.
      if (Interface_Union_List.Ni) {
	error_flag = BCs_Interface(nb,Time);
	if (error_flag) return error_flag;
      }
      // Apply domain boundary conditions.
      BCs(Local_SolnBlk[nb],*IP);
    }
  }

  // Boundary conditions applied successfully.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::BCs_Interface --                             *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
inline int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
BCs_Interface(const int &nb, const double &Time) {
  cout << endl << " ERROR: Explicit specialization required for the application of";
  cout << endl << "        interface boundary conditions in a solution block.";
  cout.flush();
  return 1;
}

/**********************************************************************
 * EmbeddedBoundaries2D:: CFL --                                      *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
inline double EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
CFL(const double &Time) {
  cout << endl << " ERROR: Explicit specialization required for the calculation of the";
  cout << endl << "        CFL criteria/time-step.";
  cout.flush();
  assert(1==0);
  return ZERO;
}

/**********************************************************************
 * EmbeddedBoundaries2D::L1_Norm_Residual --                          *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
inline double EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
L1_Norm_Residual(void) {

  double l1_norm = ZERO;

  // Calculate the L1-norm.  Sum the L1-norm for each solution block.
  for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
    if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      for (int j = Local_SolnBlk[nb].JCl; j <= Local_SolnBlk[nb].JCu; j++) {
	for (int i = Local_SolnBlk[nb].ICl; i <= Local_SolnBlk[nb].ICu; i++) {
	  if (Mesh[nb].cell_status[i][j] == CELL_STATUS_ACTIVE) {
	    l1_norm += fabs(Local_SolnBlk[nb].dUdt[i][j][0][Local_SolnBlk[nb].residual_variable]);
	  }
	}
      }
    }
  }

  // Return the L1-norm.
  return l1_norm;

}

/**********************************************************************
 * EmbeddedBoundaries2D::L2_Norm_Residual --                          *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
inline double EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
L2_Norm_Residual(void) {

  double l2_norm = ZERO;

  // Sum the square of the L2-norm for each solution block.
  for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
    if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      for (int j = Local_SolnBlk[nb].JCl; j <= Local_SolnBlk[nb].JCu; j++) {
	for (int i = Local_SolnBlk[nb].ICl; i <= Local_SolnBlk[nb].ICu; i++) {
	  if (Mesh[nb].cell_status[i][j] == CELL_STATUS_ACTIVE) {
	    l2_norm += sqr(Local_SolnBlk[nb].dUdt[i][j][0][Local_SolnBlk[nb].residual_variable]);
	  }
	}
      }
    }
  }

  // Return the L2-norm.
  return sqrt(l2_norm);

}

/**********************************************************************
 * EmbeddedBoundaries2D::Max_Norm_Residual --                         *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
inline double EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Max_Norm_Residual(void) {

  double max_norm = ZERO;

  // Find the maximum norm for all solution blocks.
  for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
    if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      for (int j = Local_SolnBlk[nb].JCl; j <= Local_SolnBlk[nb].JCu; j++) {
	for (int i = Local_SolnBlk[nb].ICl; i <= Local_SolnBlk[nb].ICu; i++) {
	  if (Mesh[nb].cell_status[i][j] == CELL_STATUS_ACTIVE) {
	    max_norm = max(max_norm,fabs(Local_SolnBlk[nb].dUdt[i][j][0][Local_SolnBlk[nb].residual_variable]));
	  }
	}
      }
    }
  }

  // Return the maximum norm.
  return max_norm;

}

/**********************************************************************
 * EmbeddedBoundaries2D::dAdt --                                      *
 *                                                                    *
 * The routine calculates the rate of area change in cell due to the  *
 * motion of an embedded boundary.  The rate of area change, dAdt, is *
 * determined by the geometric conservation law:                      *
 *                                                                    *
 *                 ddt int dV = int_S w dot n dS                      *
 *                                                                    *
 * which states that the rate of change of cell area is equal to the  *
 * contour length swept by the motion of the surface.                 *
 *                                                                    *
 * Thomas and Lombard (AIAA J. Vol. 17 No. 10 1979)                   *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
inline double EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
dAdt(const int &nb, const int &i, const int &j, const double &Time) {

  int Ni;
  double dA = ZERO;

  // NORTH face area change flux.
  if (Mesh[nb].cell_status[i  ][j  ] == CELL_STATUS_ACTIVE &&
      Mesh[nb].cell_status[i  ][j+1] != CELL_STATUS_ACTIVE &&
      Local_SolnBlk[nb].Grid.lfaceN(i,j) > ZERO) {
    Ni = Mesh[nb].cell_status[i  ][j+1];
    dA += (Local_SolnBlk[nb].Grid.lfaceN(i,j)/Local_SolnBlk[nb].Grid.Cell[i][j].A)*
          dot(Interface_Union_List[Ni].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceN(i,j),Time),
	      Local_SolnBlk[nb].Grid.nfaceN(i,j));
  }

  // SOUTH face area change flux.
  if (Mesh[nb].cell_status[i  ][j  ] == CELL_STATUS_ACTIVE &&
      Mesh[nb].cell_status[i  ][j-1] != CELL_STATUS_ACTIVE &&
      Local_SolnBlk[nb].Grid.lfaceS(i,j) > ZERO) {
    Ni = Mesh[nb].cell_status[i  ][j-1];
    dA += (Local_SolnBlk[nb].Grid.lfaceS(i,j)/Local_SolnBlk[nb].Grid.Cell[i][j].A)*
          dot(Interface_Union_List[Ni].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceS(i,j),Time),
	      Local_SolnBlk[nb].Grid.nfaceS(i,j));
  }

  // EAST face area change flux.
  if (Mesh[nb].cell_status[i  ][j] == CELL_STATUS_ACTIVE &&
      Mesh[nb].cell_status[i+1][j] != CELL_STATUS_ACTIVE &&
      Local_SolnBlk[nb].Grid.lfaceE(i,j) > ZERO) {
    Ni = Mesh[nb].cell_status[i+1][j  ];
    dA += (Local_SolnBlk[nb].Grid.lfaceE(i,j)/Local_SolnBlk[nb].Grid.Cell[i][j].A)
          *dot(Interface_Union_List[Ni].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceE(i,j),Time),
	       Local_SolnBlk[nb].Grid.nfaceE(i,j));
  }

  // WEST face area change flux.
  if (Mesh[nb].cell_status[i  ][j  ] == CELL_STATUS_ACTIVE &&
      Mesh[nb].cell_status[i-1][j  ] != CELL_STATUS_ACTIVE &&
      Local_SolnBlk[nb].Grid.lfaceW(i,j) > ZERO) {
    Ni = Mesh[nb].cell_status[i-1][j  ];
    dA += (Local_SolnBlk[nb].Grid.lfaceW(i,j)/Local_SolnBlk[nb].Grid.Cell[i][j].A)*
          dot(Interface_Union_List[Ni].Determine_Interface_Velocity(Local_SolnBlk[nb].Grid.xfaceW(i,j),Time),
	      Local_SolnBlk[nb].Grid.nfaceW(i,j));
  }

  // Return area change;
  return dA;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Linear_Least_Squares_Reconstruction --       *
 *                                                                    *
 * Peforms the reconstruction of a limited piecewise linear solution  *
 * state within a given cell (i,j) of the computational mesh for the  *
 * specified quadrilateral solution block.  A least squares approach  *
 * is used in the evaluation of the unlimited solution gradients.     *
 * Several slope limiters may be used.                                *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
inline void EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Linear_Least_Squares_Reconstruction(const int &nb,
				    const int &i,
				    const int &j) {

  int npts = 0, npts2 = 0, motion_flag = OFF;
  double u0Min, u0Max, uQuad[4], phi;
  double DxDx_ave, DxDy_ave, DyDy_ave;
  Vector2D dX;
  pState DU, DUDx_ave, DUDy_ave;
  int NUM_VAR = Local_SolnBlk[nb].NumVar();
  int at_boundary_flag = OFF;  //I was finding that limiting near the boundary
                               //was giving me very bad results, it may be different
                               //for you.    ~james

  // If the current cell is not an internal cell then exit immediately.
  if (i == Local_SolnBlk[nb].ICl-Local_SolnBlk[nb].Nghost ||
      i == Local_SolnBlk[nb].ICu+Local_SolnBlk[nb].Nghost ||
      j == Local_SolnBlk[nb].JCl-Local_SolnBlk[nb].Nghost ||
      j == Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost ||
      Mesh[nb].cell_status[i][j] != CELL_STATUS_ACTIVE) {
    Local_SolnBlk[nb].dWdx[i][j].Vacuum();
    Local_SolnBlk[nb].dWdy[i][j].Vacuum();
    Local_SolnBlk[nb].phi[i][j].Vacuum();
    return ;
  }

  // Determine the number of points in the reconstruction stencil and
  // set the motion flag if a moving embedded boundary is detected.
  for (int jj = j-1; jj <= j+1; jj++) {
    for (int ii = i-1; ii <= i+1; ii++) {
      if (Mesh[nb].cell_status[ii][jj] == CELL_STATUS_ACTIVE && !(ii == i && jj == j)) {
	i_index[npts] = ii; j_index[npts] = jj; npts++;
      } else if (Mesh[nb].cell_status[ii][jj] != CELL_STATUS_ACTIVE) {
	at_boundary_flag = ON;
	if (Interface_Union_List[Mesh[nb].cell_status[ii][jj]].Motion != INTERFACE_MOTION_STATIONARY) motion_flag = ON;
      }
    }
  }

// #ifdef _EB_PARALLEL_DEBUG_
//   dout << endl << Local_SolnBlk[nb].Grid.Cell[i][j].Xc << " npts = " << npts;
//   for (int n = 0; n < npts; n++) dout << endl << " -> " << i_index[n] << " " << j_index[n];
// #endif

// #ifdef _EB_PARALLEL_DEBUG_
//   if (npts != 8) {
//     dout << endl << "ZONE T = \"(" << i << "," << j << "," << Local_Solution_Block_List->Block[nb].gblknum << ")\" ";
//     dout << endl << "I = " << npts+1 << " \\";
//     dout << endl << "J = 1 \\";
//     dout << endl << Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
//     for (int n = 0; n < npts; n++) dout << endl << Local_SolnBlk[nb].Grid.Cell[i_index[n]][j_index[n]].Xc;
//   }
// #endif

  if (npts == 1 || npts == 2) {
    npts = 8;
    i_index[0] = i-1; j_index[0] = j-1;
    i_index[1] = i  ; j_index[1] = j-1;
    i_index[2] = i+1; j_index[2] = j-1;
    i_index[3] = i-1; j_index[3] = j  ;
    i_index[4] = i+1; j_index[4] = j  ;
    i_index[5] = i-1; j_index[5] = j+1;
    i_index[6] = i  ; j_index[6] = j+1;
    i_index[7] = i+1; j_index[7] = j+1;
  }

  // Perform the linear least squares reconstruction.
  if (npts > 0) {

    DUDx_ave.Vacuum();
    DUDy_ave.Vacuum();
    DxDx_ave = ZERO;
    DxDy_ave = ZERO;
    DyDy_ave = ZERO;

    for (int n2 = 0; n2 < npts; n2++) {
      dX = Local_SolnBlk[nb].Grid.Cell[i_index[n2]][j_index[n2]].Xc -
	   Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
      DU = Local_SolnBlk[nb].W[i_index[n2]][j_index[n2]] - Local_SolnBlk[nb].W[i][j];
      DUDx_ave += DU*dX.x;
      DUDy_ave += DU*dX.y;
      DxDx_ave += dX.x*dX.x;
      DxDy_ave += dX.x*dX.y;
      DyDy_ave += dX.y*dX.y;
    }

    DUDx_ave /= double(npts);
    DUDy_ave /= double(npts);
    DxDx_ave /= double(npts);
    DxDy_ave /= double(npts);
    DyDy_ave /= double(npts);

    Local_SolnBlk[nb].dWdx[i][j] = (DUDx_ave*DyDy_ave-DUDy_ave*DxDy_ave)/
                                   (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
    Local_SolnBlk[nb].dWdy[i][j] = (DUDy_ave*DxDx_ave-DUDx_ave*DxDy_ave)/
                                   (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);

    // Calculate slope limiters.
    if (motion_flag || IP->i_Limiter == LIMITER_ZERO) {
      Local_SolnBlk[nb].phi[i][j].Vacuum();

    } else if (!Local_SolnBlk[nb].Freeze_Limiter) {
      for (int n = 1; n <= NUM_VAR; n++) {
	u0Min = Local_SolnBlk[nb].W[i][j][n];
	u0Max = u0Min;
	for (int n2 = 0; n2 < npts; n2++) {
	  u0Min = min(u0Min,Local_SolnBlk[nb].W[i_index[n2]][j_index[n2]][n]);
	  u0Max = max(u0Max,Local_SolnBlk[nb].W[i_index[n2]][j_index[n2]][n]);
	}

	npts2 = 0;
	if (Local_SolnBlk[nb].Grid.lfaceE(i,j) > ZERO) {
	  dX = Local_SolnBlk[nb].Grid.xfaceE(i,j) - Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	  uQuad[npts2] = Local_SolnBlk[nb].W[i][j][n] + Local_SolnBlk[nb].dWdx[i][j][n]*dX.x +
	                                                Local_SolnBlk[nb].dWdy[i][j][n]*dX.y;
	  npts2++;
	}
	if (Local_SolnBlk[nb].Grid.lfaceW(i,j) > ZERO) {
	  dX = Local_SolnBlk[nb].Grid.xfaceW(i,j) - Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	  uQuad[npts2] = Local_SolnBlk[nb].W[i][j][n] + Local_SolnBlk[nb].dWdx[i][j][n]*dX.x +
	                                                Local_SolnBlk[nb].dWdy[i][j][n]*dX.y;
	  npts2++;
	}
	if (Local_SolnBlk[nb].Grid.lfaceN(i,j) > ZERO) {
	  dX = Local_SolnBlk[nb].Grid.xfaceN(i,j) - Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	  uQuad[npts2] = Local_SolnBlk[nb].W[i][j][n] + Local_SolnBlk[nb].dWdx[i][j][n]*dX.x +
	                                                Local_SolnBlk[nb].dWdy[i][j][n]*dX.y;
	  npts2++;
	}
	if (Local_SolnBlk[nb].Grid.lfaceS(i,j) > ZERO) {
	  dX = Local_SolnBlk[nb].Grid.xfaceS(i,j) - Local_SolnBlk[nb].Grid.Cell[i][j].Xc;
	  uQuad[npts2] = Local_SolnBlk[nb].W[i][j][n] + Local_SolnBlk[nb].dWdx[i][j][n]*dX.x +
	                                                Local_SolnBlk[nb].dWdy[i][j][n]*dX.y;
	  npts2++;
	}

	if(at_boundary_flag == OFF) {
	  switch(IP->i_Limiter) {
	  case LIMITER_ONE :
	    Local_SolnBlk[nb].phi[i][j][n] = ONE;
	    break;
	  case LIMITER_ZERO :
	    Local_SolnBlk[nb].phi[i][j][n] = ZERO;
	    break;
	  case LIMITER_BARTH_JESPERSEN :
	    Local_SolnBlk[nb].phi[i][j][n] = Limiter_BarthJespersen(uQuad,Local_SolnBlk[nb].W[i][j][n],u0Min,u0Max,npts2);
	    break;
	  case LIMITER_VENKATAKRISHNAN :
	    Local_SolnBlk[nb].phi[i][j][n] = Limiter_Venkatakrishnan(uQuad,Local_SolnBlk[nb].W[i][j][n],u0Min,u0Max,npts2);
	    break;
	  case LIMITER_VANLEER :
	    Local_SolnBlk[nb].phi[i][j][n] = Limiter_VanLeer(uQuad,Local_SolnBlk[nb].W[i][j][n],u0Min,u0Max,npts2);
	    break;
	  case LIMITER_VANALBADA :
	    Local_SolnBlk[nb].phi[i][j][n] = Limiter_VanAlbada(uQuad,Local_SolnBlk[nb].W[i][j][n],u0Min,u0Max,npts2);
	    break;
	  default:
	    Local_SolnBlk[nb].phi[i][j][n] = Limiter_BarthJespersen(uQuad,Local_SolnBlk[nb].W[i][j][n],u0Min,u0Max,npts2);
	    break;
	  };
	} else {
	  Local_SolnBlk[nb].phi[i][j][n] = ONE;
	}
      }
    }

  } else {
    Local_SolnBlk[nb].dWdx[i][j].Vacuum();
    Local_SolnBlk[nb].dWdy[i][j].Vacuum();
    Local_SolnBlk[nb].phi[i][j].Vacuum();

  }

}

/**********************************************************************
 * Routine: Linear_Least_Squares_Reconstruction --                    *
 *                                                                    *
 * Peforms the reconstruction of a limited piecewise linear solution  *
 * state within each cell of the computational mesh of the specified  *
 * quadrilateral solution block.  A least squares approach is used in *
 * the evaluation of the unlimited solution gradients.  Several slope *
 * limiters may be used.                                              *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
inline void EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Linear_Least_Squares_Reconstruction(const int &nb) {

  // Carry out the limited solution reconstruction in each cell of the 
  // computational mesh.
  for (int j = Local_SolnBlk[nb].JCl-Local_SolnBlk[nb].Nghost+1; j <= Local_SolnBlk[nb].JCu+Local_SolnBlk[nb].Nghost-1; j++) {
    for (int i = Local_SolnBlk[nb].ICl-Local_SolnBlk[nb].Nghost+1; i <= Local_SolnBlk[nb].ICu+Local_SolnBlk[nb].Nghost-1; i++) {
      Linear_Least_Squares_Reconstruction(nb,i,j);
    }
  }

}

/**********************************************************************
 * EmbeddedBoundaries2D::Residual_Smoothing --                        *
 *                                                                    *
 * This routine performs implicit residual smoothing iterations to    *
 * solution residual for a 1D array of quadrilateral multi-block      *
 * solution blocks.  Note that only residuals of interior cells are   *
 * smoothed and residuals for cells adjacent to boundaries (embedded  *
 * or otherwise) are not smoothed.                                    *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
inline void EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Residual_Smoothing(const int &i_stage) {

  int k_residual;

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
  };

  for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
    if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {

      for (int j = Local_SolnBlk[nb].JCl+1; j <= Local_SolnBlk[nb].JCu-1; j++) {
	for (int i = Local_SolnBlk[nb].ICl+1; i <= Local_SolnBlk[nb].ICu-1; i++) {
	  Local_SolnBlk[nb].dUdt[i][j][k_residual+1] = Local_SolnBlk[nb].dUdt[i][j][k_residual];
	}
      }

      for (int n = 0; n < IP->Residual_Smoothing_Gauss_Seidel_Iterations; n++) {
	for (int j = Local_SolnBlk[nb].JCl+1; j <= Local_SolnBlk[nb].JCu-1; j++) {
	  for (int i = Local_SolnBlk[nb].ICl+1; i <= Local_SolnBlk[nb].ICu-1; i++) {
	    if (Mesh[nb].cell_status[i][j] == CELL_STATUS_ACTIVE &&
		Mesh[nb].cell_status[i+1][j] == CELL_STATUS_ACTIVE &&
		Mesh[nb].cell_status[i-1][j] == CELL_STATUS_ACTIVE &&
		Mesh[nb].cell_status[i][j+1] == CELL_STATUS_ACTIVE &&
		Mesh[nb].cell_status[i][j-1] == CELL_STATUS_ACTIVE) {
	      Local_SolnBlk[nb].dUdt[i][j][k_residual+1] = (Local_SolnBlk[nb].dUdt[i][j][k_residual] +
							    IP->Residual_Smoothing_Epsilon*(Local_SolnBlk[nb].dUdt[i  ][j-1][k_residual+1] +
											    Local_SolnBlk[nb].dUdt[i-1][j  ][k_residual+1] +
											    Local_SolnBlk[nb].dUdt[i+1][j  ][k_residual+1] +
											    Local_SolnBlk[nb].dUdt[i  ][j+1][k_residual+1]))/(ONE + FOUR*IP->Residual_Smoothing_Epsilon);
	    }
	  }
	}
      }

      for (int j = Local_SolnBlk[nb].JCl+1; j <= Local_SolnBlk[nb].JCu-1; j++) {
	for (int i = Local_SolnBlk[nb].ICl+1; i <= Local_SolnBlk[nb].ICu-1; i++) {
	  Local_SolnBlk[nb].dUdt[i][j][k_residual] = Local_SolnBlk[nb].dUdt[i][j][k_residual+1];
	}
      }

    }
  }

  // Residual smoothing successfully performed.

}

/**********************************************************************
 * EmbeddedBoundaries2D::dUdt_Residual_Evaluation --                  *
 *                                                                    *
 * This routine evaluates the residual for the specified solution     *
 * block using a 2nd-order limited upwind finite-volume spatial       *
 * discretization scheme.  The residual is stored in dUdt[][][0].     *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
dUdt_Residual_Evaluation(const double &Time) {
  cout << endl << " ERROR: Explicit specialization required for the evaluation for the ";
  cout << endl << "        right-hand-side residual evaluation.";
  cout.flush();
}

/**********************************************************************
 * EmbeddedBoundaries2D::dUdt_Multistage_Explicit --                  *
 *                                                                    *
 * This routine evaluates the stage solution residual for a 1D array  *
 * of 2D quadrilateral multi-block solution blocks.  A variety of     *
 * multistage explicit time integration and upwind finite-volume      *
 * spatial discretization procedures can be used depending on the     *
 * specified input values.                                            *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
dUdt_Multistage_Explicit(const int &i_stage, const double &Time) {
  cout << endl << " ERROR: Explicit specialization required for the evaluation for the ";
  cout << endl << "        stage residual of a mutli-stage explicit time-stepping scheme.";
  cout.flush();
  return 1;
}

/**********************************************************************
 * EmbeddedBoundaries2D::Update_Solution_Multistage_Explicit --       *
 *                                                                    *
 * This routine updates the solution for a 1D array of 2D             *
 * quadrilateral multi-block solution blocks. A variety of multistage *
 * explicit time integration and upwind finite-volume spatial         *
 * discretization procedures can be used depending on the specified   *
 * input values.                                                      *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Update_Solution_Multistage_Explicit(const int &i_stage) {

  int k_residual;
  double omega, residual_reduction_factor;
  int NUM_VAR = Local_SolnBlk[0].NumVar();

  // Memory for linear system solver.
  DenseMatrix dRdU(NUM_VAR,NUM_VAR);
  DenseSystemLinEqs LinSys;

  // Allocate memory for linear system solver.
  LinSys.allocate(NUM_VAR);

  // Perform update of solution variables for stage i_stage of N_Stage
  // scheme.

  // Evaluate the time step fraction and residual storage location for
  // the stage.
  switch(IP->i_Time_Integration) {
  case TIME_STEPPING_EXPLICIT_EULER :
    omega = Runge_Kutta(i_stage,IP->N_Stage);
    k_residual = 0;
    break;
  case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
    omega = Runge_Kutta(i_stage,IP->N_Stage);
    k_residual = 0;
    break;
  case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA :
    omega = Runge_Kutta(i_stage,IP->N_Stage);
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
    omega = MultiStage_Optimally_Smoothing(i_stage,IP->N_Stage,IP->i_Limiter);
    k_residual = 0;
    break;
  default:
    omega = Runge_Kutta(i_stage,IP->N_Stage);
    k_residual = 0;
    break;
  };

  // Update the solution for each solution block.
  for (int nb = 0; nb < Local_Solution_Block_List->Nblk; nb++) {
    if (Local_Solution_Block_List->Block[nb].used == ADAPTIVEBLOCK2D_USED) {

      // Update solution variables for this stage.
      for (int j = Local_SolnBlk[nb].JCl; j <= Local_SolnBlk[nb].JCu; j++) {
	for (int i = Local_SolnBlk[nb].ICl; i <= Local_SolnBlk[nb].ICu; i++) {

	  if (Mesh[nb].cell_status[i][j] == CELL_STATUS_ACTIVE) {

	    // Explicit update.
	    if (IP->Local_Time_Stepping == GLOBAL_TIME_STEPPING ||
		IP->Local_Time_Stepping == SCALAR_LOCAL_TIME_STEPPING) {
	      Local_SolnBlk[nb].U[i][j] = Local_SolnBlk[nb].Uo[i][j] + omega*Local_SolnBlk[nb].dUdt[i][j][k_residual];
	      if(Local_SolnBlk[nb].W[i][j].analytically_inverted_relaxation()) {  //so far this is only used for Gaussian2D
		Local_SolnBlk[nb].W[i][j] = W(Local_SolnBlk[nb].U[i][j]);
		Local_SolnBlk[nb].W[i][j].relax(IP->CFL_Number*Local_SolnBlk[nb].dt[i][j],
						i_stage,W(Local_SolnBlk[nb].Uo[i][j]));
		Local_SolnBlk[nb].U[i][j] = U(Local_SolnBlk[nb].W[i][j]);	
	      }

	    }

	    // Point implicit formulation: set up system of equations 
	    // and include source Jacobian in the LHS matrix.
	    if (IP->Local_Time_Stepping == SEMI_IMPLICIT_LOCAL_TIME_STEPPING) {

	      // Evaluate the source Jacobian.
	      dRdU.zero();

	      // Include contributions from all required source Jacobians.
	      Local_SolnBlk[nb].W[i][j].dSdU(dRdU,
					     Local_SolnBlk[nb].Grid.Cell[i][j].Xc,
					     Local_SolnBlk[nb].dWdx[i][j],
					     Local_SolnBlk[nb].dWdy[i][j],
					     Local_SolnBlk[nb].Axisymmetric);

	      // Include the source Jacobian in the LHS matrix.
	      LinSys.A.identity();
	      LinSys.A -= (omega*IP->CFL_Number*Local_SolnBlk[nb].dt[i][j])*dRdU;

	      // Set the explicit residual as the RHS for the point implicit
	      // formulation (already contains the CFL number).
	      for (int k = 1; k <= NUM_VAR; k++)
		LinSys.b(k-1) = omega*Local_SolnBlk[nb].dUdt[i][j][k_residual][k];

	      // Solve system of equations using LU decomposition Gaussian 
	      // elimination procedure.
	      LinSys.solve(LU_DECOMPOSITION);

	      // Update the conserved solution variables.
	      for (int k = 1; k <= NUM_VAR; k++)
		Local_SolnBlk[nb].U[i][j][k] = Local_SolnBlk[nb].Uo[i][j][k] + LinSys.x(k-1);

	    }

	    // Perform residual reductions (reduce the time-step) if
	    // necessary for scalar or semi-implicit local time-stepping.
	    if ((IP->Local_Time_Stepping == SCALAR_LOCAL_TIME_STEPPING ||
		 IP->Local_Time_Stepping == SEMI_IMPLICIT_LOCAL_TIME_STEPPING) &&
		Local_SolnBlk[nb].U[i][j].Unphysical_Properties()) {
	      residual_reduction_factor = ONE;
	      for (int n_residual_reduction = 0; n_residual_reduction < 10; n_residual_reduction++) {
		// Reduce the residual by half.
		residual_reduction_factor = HALF*residual_reduction_factor;
		Local_SolnBlk[nb].dt[i][j] = residual_reduction_factor*Local_SolnBlk[nb].dt[i][j];
		Local_SolnBlk[nb].dUdt[i][j][k_residual] = residual_reduction_factor*Local_SolnBlk[nb].dUdt[i][j][k_residual];
		// Re-perform the solution update.
		if (IP->Local_Time_Stepping == SCALAR_LOCAL_TIME_STEPPING) {
		  Local_SolnBlk[nb].U[i][j] = Local_SolnBlk[nb].Uo[i][j] + omega*Local_SolnBlk[nb].dUdt[i][j][k_residual];
		} else if (IP->Local_Time_Stepping == SEMI_IMPLICIT_LOCAL_TIME_STEPPING) {
		  LinSys.A.identity();
		  LinSys.A -= (omega*IP->CFL_Number*Local_SolnBlk[nb].dt[i][j])*dRdU;
		  for (int k = 1; k <= NUM_VAR; k++)
		    LinSys.b(k-1) = omega*Local_SolnBlk[nb].dUdt[i][j][k_residual][k];
		  LinSys.solve(LU_DECOMPOSITION);
		  for (int k = 1; k <= NUM_VAR; k++)
		    Local_SolnBlk[nb].U[i][j][k] = Local_SolnBlk[nb].Uo[i][j][k] + LinSys.x(k-1);
		}
		if (!Local_SolnBlk[nb].U[i][j].Unphysical_Properties()) break;
		if (n_residual_reduction == 10) cout << "n_residual_reductions = " << n_residual_reduction
						     << " @ Xc =" << Local_SolnBlk[nb].Grid.Cell[i][j].Xc << endl;
	      }
	    }

	    // Error checking.
	    if (Local_SolnBlk[nb].U[i][j].Unphysical_Properties()) {
	      cout << "\n EmbeddedBoundaries2D ERROR: Unphysical state property detected in"
		   << " cell = (" << Local_Solution_Block_List->Block[nb].gblknum << ":" << i << "," << j << ") " << endl
		   << " X    = " << Local_SolnBlk[nb].Grid.Cell[i][j].Xc     << endl
		   << " U    = " << Local_SolnBlk[nb].U[i][j]                << endl 
		   << " W    = " << W(Local_SolnBlk[nb].U[i][j])             << endl
		   << " Uo   = " << Local_SolnBlk[nb].Uo[i][j]               << endl
		   << " Wo   = " << W(Local_SolnBlk[nb].Uo[i][j])            << endl
		   << " dUdt = " << Local_SolnBlk[nb].dUdt[i][j][k_residual] << endl;
//	      cout << endl << Local_SolnBlk[nb].U[i][j].p();
//	      cout << endl << Local_SolnBlk[nb].U[i][j].e();
//	      cout << endl << Local_SolnBlk[nb].U[i][j].T();
//            I commented these out because they don't exist in my state (Gaussian2D)
	      cout.flush();
	      return i;
	    }

	    // Update the primitive solution variables.
	    Local_SolnBlk[nb].W[i][j] = W(Local_SolnBlk[nb].U[i][j]);

	  }

	}
      }

    }
  }

  // Deallocate memory for linear system solver.
  LinSys.deallocate();

  // Quadrilateral multi-block solution blocks successfully updated.
  return 0;

}

/**********************************************************************
 * EmbeddedBoundaries2D::Execute --                                   *
 *                                                                    *
 * Execution of the explicit solution of the desired system of        *
 * hyperbolic/elliptic conservation laws in a domain with embedded    *
 * boundaries.                                                        *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters>::
Execute(const int &batch_flag,
	int &number_of_time_steps,
	int &evolution_counter,
	int &number_of_level_set_time_steps,
	double &Time,
	double &levelset_Time,
	CPUTime &processor_cpu_time,
	CPUTime &total_cpu_time,
	ofstream &residual_file) {

  // Other local solution variables.
  int first_step,
      error_flag,
      perform_explicit_time_marching, limiter_freezing_flag;
  double dTime;
  double residual_l2norm_first, residual_ratio;
  double residual_l1_norm, residual_l2_norm, residual_max_norm;

  const int NUM_VAR = IP->Wo.NumVar();

  // Open residual file.
  first_step = 1;
  limiter_freezing_flag = OFF;

  if (CFFC_Primary_MPI_Processor()) {
    error_flag = Open_Progress_File(residual_file,
				    IP->Output_File_Name,
				    number_of_time_steps);
    if (error_flag) {
      cout << "\n ERROR: Unable to open residual file for the calculation.\n";
      cout.flush();
    }
  }
  // MPI barrier to ensure processor synchronization.  
  CFFC_Barrier_MPI();
  CFFC_Broadcast_MPI(&error_flag,1);
  if (error_flag) return error_flag;
  // Reset the CPU time.
  processor_cpu_time.reset();

  // Perform required number of iterations (time steps).
  if ((!IP->Time_Accurate && IP->Maximum_Number_of_Time_Steps > 0 &&
       number_of_time_steps < IP->Maximum_Number_of_Time_Steps) ||
      (IP->Time_Accurate && IP->Time_Max > Time)) {

    if (!batch_flag) cout << "\n Beginning computations on "
			  << Date_And_Time() << ".\n\n";

    if ((!IP->Time_Accurate && IP->Maximum_Number_of_Time_Steps > 0 &&
	 number_of_time_steps < IP->Maximum_Number_of_Time_Steps) ||
	(IP->Time_Accurate && IP->Time_Max > Time)) {
      perform_explicit_time_marching = ON;
    } else {
      perform_explicit_time_marching = OFF;
    }

    while (perform_explicit_time_marching) {

      // Periodically refine the mesh (AMR).
      if (IP->AMR) {
	if (!first_step &&
	    number_of_time_steps-IP->AMR_Frequency*
	    (number_of_time_steps/IP->AMR_Frequency) == 0) {
	  if (!batch_flag) cout << "\n\n Refining Grid.  Performing adaptive mesh refinement at n = "
				<< number_of_time_steps << ".";
	  Evaluate_Limiters(Local_SolnBlk,*Local_Solution_Block_List);
	  error_flag = Adaptive_Mesh_Refinement(ON,ON);
	  if (error_flag) {
	    cout << "\n ERROR: AMR error number " << error_flag << " on processor "
		 << Local_Solution_Block_List->ThisCPU << "." << endl;
	  }
	  error_flag = CFFC_OR_MPI(error_flag);
	  if (error_flag) {
	    error_flag = Output_Tecplot(Local_SolnBlk,
					*Local_Solution_Block_List,
					*IP,
					number_of_time_steps,
					Time);
	    return 1;
	  }
	  error_flag = Boundary_Conditions(Time);
	  if (error_flag) {
	    cout << "\n ERROR: Boundary conditions error on processor "
		 << Local_Solution_Block_List->ThisCPU << "." << endl;
	  }
	  error_flag = CFFC_OR_MPI(error_flag);
	  if (error_flag) return error_flag;
	  // Re-prescribe boundary data consistent with newly refined and
	  // coarsened solution blocks.
	  Boundary_Conditions(Time);
	  if (!batch_flag) {
	    cout << "\n New multi-block solution-adaptive quadrilateral mesh statistics: ";
	    cout << "\n  -> Number of Root Blocks i-direction: "
		 << QuadTree->NRi;
	    cout << "\n  -> Number of Root Blocks j-direction: "
		 << QuadTree->NRj;
	    cout << "\n  -> Total Number of Used Blocks: "
		 << QuadTree->countUsedBlocks();
	    cout << "\n  -> Total Number of Computational Cells: "
		 << QuadTree->countUsedCells();
	    cout << "\n  -> Number of Mesh Refinement Levels: "
		 << QuadTree->highestRefinementLevel()+1;
	    cout << "\n  -> Refinement Efficiency: "
		 << QuadTree->efficiencyRefinement() << endl;
	  }
	}
      }

      // Determine local and global time steps.
      dTime = CFL(Time);
      // Find global minimum time step for all processors.
      dTime = CFFC_Minimum_MPI(dTime);
      if (IP->Time_Accurate) {
	if ((IP->i_Time_Integration != TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING) &&
	    (Time + IP->CFL_Number*dTime > IP->Time_Max)) {
	  dTime = (IP->Time_Max-Time)/IP->CFL_Number;
	} else if (Time + IP->CFL_Number*dTime*
		   MultiStage_Optimally_Smoothing(IP->N_Stage,
						  IP->N_Stage,
						  IP->i_Limiter) > IP->Time_Max) {
	  dTime = (IP->Time_Max-Time)/(IP->CFL_Number*
				       MultiStage_Optimally_Smoothing(IP->N_Stage,
								      IP->N_Stage,
								      IP->i_Limiter));
	}
      }
      // Set global time step.
      if (!IP->Local_Time_Stepping) {
	Set_Global_TimeStep(Local_SolnBlk,
			    *Local_Solution_Block_List,
			    dTime);
      }

      // Determine the L1, L2, and max norms of the solution residual.
      residual_l1_norm = L1_Norm_Residual();
      residual_l1_norm = CFFC_Summation_MPI(residual_l1_norm);

      residual_l2_norm = L2_Norm_Residual();
      residual_l2_norm = sqr(residual_l2_norm);
      residual_l2_norm = CFFC_Summation_MPI(residual_l2_norm);
      residual_l2_norm = sqrt(residual_l2_norm);

      residual_max_norm = Max_Norm_Residual();
      residual_max_norm = CFFC_Maximum_MPI(residual_max_norm);

      // Update CPU time used for the calculation so far.
      processor_cpu_time.update();
      total_cpu_time.cput = CFFC_Summation_MPI(processor_cpu_time.cput);

      // Periodically save restart solution files.
      if (!first_step &&
	  number_of_time_steps-IP->Restart_Solution_Save_Frequency*
	  (number_of_time_steps/IP->Restart_Solution_Save_Frequency) == 0) {
	if (!batch_flag) 
	  cout << "\n\n  Saving solution to restart data file(s) after"
	       << " n = " << number_of_time_steps << " steps (iterations).";

	//  Save and delete old restart files in compressed archive (just in case)
	if (CFFC_Primary_MPI_Processor()) {
	  cout << "\n  Creating compressed archive of (and deleting) old restarts.";
	  System::Compress_Restart();
	  cout << "\n  Writing new restart files.";
	  cout.flush();
	}
	CFFC_Barrier_MPI(); // MPI barrier so that other processors do
	                      // not start over-writing restarts

	if (CFFC_Primary_MPI_Processor()) {
	  System::Set_Restart_Flag();  //Set flag to indicate a restart is being saved
	}

	// Write the quadtree restart file.
	error_flag = Write_QuadTree(*QuadTree,*IP);
	if (error_flag) {
	  cout << "\n ERROR: Unable to open quadtree data file on processor "
	       << Local_Solution_Block_List->ThisCPU << "." << endl;
	}
	error_flag = CFFC_OR_MPI(error_flag);
	if (error_flag) return error_flag;
	// Write the solution block restart files.
	error_flag = Write_Restart_Solution(Local_SolnBlk,
					    *Local_Solution_Block_List,
					    *IP,
					    number_of_time_steps,
					    Time,
					    processor_cpu_time);
	if (error_flag) {
	  cout << "\n ERROR: Unable to open restart output data file(s) on processor "
	       << Local_Solution_Block_List->ThisCPU << "." << endl;
	}
	error_flag = CFFC_OR_MPI(error_flag);
	if (error_flag) return error_flag;
	// Write the solution block restart files.
	error_flag = Write_Restart_Files(number_of_time_steps,
					 number_of_level_set_time_steps,
					 Time,
					 levelset_Time,
					 processor_cpu_time);
	if (error_flag) {
	  cout << "\n ERROR: Unable to open restart output data file(s) "
	       << "on processor " << Local_Solution_Block_List->ThisCPU << "." << endl;
	}
	error_flag = CFFC_OR_MPI(error_flag);
	if (error_flag) return error_flag;
	if (!batch_flag) cout << endl;

	if (CFFC_Primary_MPI_Processor()) {
	  System::Remove_Restart_Flag();  //Remove flag to indicate the restart is finished
	}
      }

      // Output progress information for the calculation.
      if (!batch_flag) Output_Progress_L2norm(number_of_time_steps,
					      Time*THOUSAND,
					      total_cpu_time,
					      residual_l2_norm,
					      first_step,
					      IP->Output_Progress_Frequency);
      if (CFFC_Primary_MPI_Processor() && !first_step)
	Output_Progress_to_File(residual_file,
				number_of_time_steps,
				Time*THOUSAND,
				total_cpu_time,
				residual_l1_norm,
				residual_l2_norm,
				residual_max_norm);

      // Check to see if calculations are complete.
      if (!IP->Time_Accurate && number_of_time_steps >= IP->Maximum_Number_of_Time_Steps) break;
      if (IP->Time_Accurate && Time >= IP->Time_Max) break;

      // Freeze limiters as necessary.
      if (!first_step &&
	  IP->Freeze_Limiter &&
	  !limiter_freezing_flag &&
	  residual_l2_norm <= IP->Freeze_Limiter_Residual_Level) {
	Freeze_Limiters(Local_SolnBlk,*Local_Solution_Block_List);
	limiter_freezing_flag = ON;
      }

      // Update solution for next time step using a multistage
      // time stepping scheme.
      for (int i_stage = 1; i_stage <= IP->N_Stage; i_stage++) {

	// Step 1. Exchange solution information between neighbouring blocks.
	error_flag = Send_All_Messages(Local_SolnBlk,
				       *Local_Solution_Block_List,
				       NUM_VAR,
				       OFF);
	if (error_flag) {
	  cout << "\n ERROR: Message passing error on processor "
	       << Local_Solution_Block_List->ThisCPU << "." << endl;
	}
	error_flag = CFFC_OR_MPI(error_flag);
	if (error_flag) return error_flag;

	// Step 2. Apply boundary conditions for stage.
	error_flag = Boundary_Conditions(Time);
	if (error_flag) {
	  cout << "\n ERROR: Boundary conditions error on processor "
	       << Local_Solution_Block_List->ThisCPU << "." << endl;
	}
	error_flag = CFFC_OR_MPI(error_flag);
	if (error_flag) return error_flag;

	// Step 3. Determine solution residuals for stage.
	error_flag = dUdt_Multistage_Explicit(i_stage,Time);
	if (error_flag) {
	  cout << "\n ERROR: Solution error on processor "
	       << Local_Solution_Block_List->ThisCPU << "." << endl;
	}
	error_flag = CFFC_OR_MPI(error_flag);
	if (error_flag) return error_flag;

	// Step 4. Send boundary flux corrections at block interfaces with
	// resolution changes.
	error_flag = Send_Conservative_Flux_Corrections(Local_SolnBlk,
							*Local_Solution_Block_List,
							NUM_VAR);
	if (error_flag) {
	  cout << "\n ERROR: Flux correction message passing error on processor "
	       << Local_Solution_Block_List->ThisCPU << "." << endl;
	}
	error_flag = CFFC_OR_MPI(error_flag);
	if (error_flag) return error_flag;

	// Step 5. Apply boundary flux corrections to ensure that method is
	// conservative.
	Apply_Boundary_Flux_Corrections_Multistage_Explicit(Local_SolnBlk,
							    *Local_Solution_Block_List,
							    *IP,
							    i_stage);

	// Step 6. Smooth the solution residual using implicit residual smoothing.
	if (IP->Residual_Smoothing) Residual_Smoothing(i_stage);

	// Step 7. Update solution for stage.
	error_flag = Update_Solution_Multistage_Explicit(i_stage);
	if (error_flag) {
	  cout << "\n ERROR: Solution update error on processor "
	       << Local_Solution_Block_List->ThisCPU << "." << endl;
	}
	error_flag = CFFC_OR_MPI(error_flag);
	if (error_flag) return error_flag;

      }

      // Evolve the embedded interface(s).
      if (Interface_Component_List.Ni) {
	error_flag = Compute_Interface_Location(batch_flag,
						Time,
						Time+IP->CFL_Number*dTime,
						evolution_counter,
						number_of_level_set_time_steps,
						levelset_Time);
	if (error_flag) {
	  cout << "\n ERROR: Compute interface location error on processor "
	       << Local_Solution_Block_List->ThisCPU 
	       << ".  Error number = " << error_flag << "." << endl;
	}
	error_flag = CFFC_OR_MPI(error_flag);
	if (error_flag) return error_flag;
      }

      // Update time and time step counter.
      if (first_step) first_step = 0;
      number_of_time_steps++;
      if (IP->i_Time_Integration != TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING) {
	Time += IP->CFL_Number*dTime;
      } else {
	Time += IP->CFL_Number*dTime*MultiStage_Optimally_Smoothing(IP->N_Stage,
								    IP->N_Stage,
								    IP->i_Limiter);
      }

    }

    if (!batch_flag) cout << "\n\n Computations complete on " 
			  << Date_And_Time() << "." << endl;
  }

  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

  // Update ghostcell information and prescribe boundary conditions to 
  // ensure that the solution is consistent on each block.
  error_flag = Send_All_Messages(Local_SolnBlk,
				 *Local_Solution_Block_List,
				 NUM_COMP_VECTOR2D,
				 ON);
  if (!error_flag) error_flag = Send_All_Messages(Local_SolnBlk,
						  *Local_Solution_Block_List,
						  NUM_VAR,
						  OFF);
  if (error_flag) {
    cout << "\n ERROR: Message passing error on processor "
	 << Local_Solution_Block_List->ThisCPU << "." << endl;
  }
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return error_flag;

  // Apply boundary conditions.
  Boundary_Conditions(Time);

  // Close residual file.
  if (CFFC_Primary_MPI_Processor()) error_flag = Close_Progress_File(residual_file);

  // Embedded boundary calculations complete.
  return 0;

}

#endif // _EMBEDDEDBOUNDARIES2D_INCLUDED

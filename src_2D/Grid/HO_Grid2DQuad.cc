/*!\file HO_Grid2DQuad.cc
   \brief Source file initializing/implementing member variables/functions that belong to classes defined in HO_Grid2DQuad.h */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "HO_Grid2DQuad.h"	// Include 2D high-order block grid
#include "HO_Grid2DQuad_ExecutionMode.h" // Include header file for the control class of the 2D high-order grid


//#define MODIFY_NUMBER_OF_FLUX_CALCULATION_POINTS_AT_BOUNDARIES


// ===== Member variables =====

/*!
 * This switch is used to determine whether the geometric properties
 * (i.e. centroid, area etc.) of the cells closed to boundaries are
 * determined using a high-order or a low-order representation of the boundary.
 * If ON, the boundary splines are used to compute the geometric properties (i.e. high-order accuracy).
 * If OFF, a linear representation between 2 consecutive nodes is considered (i.e. 2nd-order accuracy).
 */
int Grid2D_Quad_Block_HO::HighOrderBoundaryRepresentation = OFF;

/*!
 * This switch is used to determine whether the smoothing 
 * subroutine is applied or not. 
 * This flag acts globally, that is, on all quad blocks.
 * If ON, the smoother is applied.
 * If OFF, the smoother is not applied so a call to the smoothing subroutine leaves the mesh unchanged.
 */
int Grid2D_Quad_Block_HO::Smooth_Quad_Block_Flag = ON;

/*!
 * This switch is used to determine whether the curvilinear
 * path integrals are computed based on Gauss quadratures or
 * on an adaptive algorithm that approximates the spline geometry
 * with an increasing number of line segments until the desired
 * accuracy is obtained. 
 * If ON, the n-point Gauss quadrature integration is used.
 * If OFF, the adaptive integration based on line segments is used.
 */
int Grid2D_Quad_Block_HO::Gauss_Quad_Curvilinear_Integration = ON;

/*!
 * This switch is used to determine whether the curvilinear
 * path integrals that occur in the calculation of the area and 
 * the centroid are computed based on Gauss quadratures or
 * on an adaptive algorithm that approximates the spline geometry
 * with an increasing number of line segments until the desired
 * accuracy is obtained. 
 * The accuracy of the area and the centroid is crucial for the 
 * accuracy of the high-order geometric moments, and therefore
 * it is desired that these properties be calculated very accurately,
 * even if a more expensive algorithm is used.
 * If ON, the adaptive integration based on line segments is used for 
 *        area and centroid calculation and the rest of the moments are
 *        computed based on the Gauss quadratures.
 * If OFF, the n-point Gauss quadrature integration is used everywhere.
 */
int Grid2D_Quad_Block_HO::Mixed_Curvilinear_Integration = ON;

/*!
 * This switch is used to determine whether the calculation of
 * cell geometric properties should be done in a reference system 
 * which can minimize the numerical errors.
 * This switch can be very useful when the ratio between cell 
 * dimensions and the nodal values are large.
 * In order to minimize the numerical errors, the geometric properties
 * are computed in a reference system local to each cell and then
 * translated in the global reference system. This operation can 
 * increase the computational cost!
 * If ON, the geometric properties are computed in the LOCAL reference system.
 * If OFF, the geometric properties are computed in the GLOBAL reference system.
 */
int Grid2D_Quad_Block_HO::Minimize_Error_Calculation_Of_Geometric_Properties = OFF;

/*!
 * This flag is used to specify whether an inaccurate integration near curved boundaries is accepted
 * in case the integration cannot be performed along the real curved geometry.
 * Such a situation might arise in the calculation of integrals of functions which don't have
 * an analytical x-dependency function defined. \n
 * Turn ON if you want to accept inaccuracies (i.e. the integration is performed as if the cell edges were straight). \n
 * Turn OFF if you don't tolerate these inaccuracies. (default) \n
 */
int Grid2D_Quad_Block_HO::Tolerate_Inaccurate_Integration_Near_Curved_Boundaries = OFF;

/*!
 * This flag is used to specify whether a polygonal adaptive quadrature integration for cells 
 * near curved boundaries is accepted in case the integration cannot be performed with contour integration.
 * Turn ON if you want to allow this method. (default) \n
 * Turn OFF if you don't. \n
 */
int Grid2D_Quad_Block_HO::Polygonal_Adaptive_Quadrature_Integration_Allowed = ON;
//! Parameter to set the minimum levels of refinement during the polygonal adaptive quadrature integration
short Grid2D_Quad_Block_HO::Polygonal_Adaptive_Quadrature_Integration_Minimum_Levels = 2;

/*!
 * This flag is used to specify whether a Monte Carlo integration for cells 
 * near curved boundaries is accepted in case the integration cannot be performed with contour integration
 * or polygonal adaptive quadrature integration is not allowed.
 * The number of sample points can be set in the input parameters (see NumericalLibrary_Execution_Mode class).
 * Turn ON if you want to allow this method. (default) \n
 * Turn OFF if you don't. \n
 */
int Grid2D_Quad_Block_HO::Monte_Carlo_Integration_Allowed = OFF;

/*! Switch to force non-reflection of ghost cells regardless of the boundary conditions
 *  in order to generate a valid mesh (i.e. not crossed quadrilaterals in ghost cells).
 *  Some meshes (e.g. O-grid NACA airfoil) are not correctly generated if this flag is turned OFF.
 *  Even if the flag is ON, depending on the mesh resolution, the mesh might still be invalid, especially when
 *  more ghost cell layers and high-order geometry treatment are required.
 *  Turn ON if you want the ghost cells near solid boundaries not to represent reflection of interior cells.
 *  Turn OFF if you want the same ghost cells to represent reflection of the corresponding interior cells.
 */
int Grid2D_Quad_Block_HO::Mesh_Requiring_NonReflected_South_Ghost_Cells = OFF;


/*!
 * Variable used for storing the South-West node position in the global coordinate system
 */
Vector2D Grid2D_Quad_Block_HO::_SW_ = Vector2D(0.0);

/*!
 * Variable used for storing the South-East node position in the global coordinate system
 */
Vector2D Grid2D_Quad_Block_HO::_SE_ = Vector2D(0.0);

/*!
 * Variable used for storing the North-East node position in the global coordinate system
 */
Vector2D Grid2D_Quad_Block_HO::_NE_ = Vector2D(0.0);

/*!
 * Variable used for storing the North-West node position in the global coordinate system
 */
Vector2D Grid2D_Quad_Block_HO::_NW_ = Vector2D(0.0);

// ===== Member functions =====

/*!
 * Copy constructor. It is declared private
 */
Grid2D_Quad_Block_HO::Grid2D_Quad_Block_HO(const Grid2D_Quad_Block_HO &G)
  :Integration(this),
   NNi(0), INl(0), INu(0), NNj(0), JNl(0), JNu(0),
   NCi(0), ICl(0), ICu(0), NCj(0), JCl(0), JCu(0),
   Nghost(0), HighestReconstructionOrder(0),
   Node(NULL), Cell(NULL),
   BCtypeN(NULL), BCtypeS(NULL), BCtypeE(NULL), BCtypeW(NULL),
   BndNorthSpline(), BndSouthSpline(), BndEastSpline(), BndWestSpline(),
   ExtendWest_BndNorthSpline(), ExtendEast_BndNorthSpline(),
   ExtendWest_BndSouthSpline(), ExtendEast_BndSouthSpline(),
   ExtendNorth_BndEastSpline(), ExtendSouth_BndEastSpline(),
   ExtendNorth_BndWestSpline(), ExtendSouth_BndWestSpline(),
   BndNorthSplineInfo(NULL), BndSouthSplineInfo(NULL),
   BndEastSplineInfo(NULL), BndWestSplineInfo(NULL),
   ExtendWest_BndNorthSplineInfo(NULL), ExtendEast_BndNorthSplineInfo(NULL),
   ExtendWest_BndSouthSplineInfo(NULL), ExtendEast_BndSouthSplineInfo(NULL),
   ExtendNorth_BndEastSplineInfo(NULL), ExtendSouth_BndEastSplineInfo(NULL),
   ExtendNorth_BndWestSplineInfo(NULL), ExtendSouth_BndWestSplineInfo(NULL),
   SminN(ZERO), SmaxN(ZERO), SminS(ZERO), SmaxS(ZERO), 
   SminE(ZERO), SmaxE(ZERO), SminW(ZERO), SmaxW(ZERO),
   StretchI(0), StretchJ(0),
   BetaI(ONE), TauI(ONE),
   BetaJ(ONE), TauJ(ONE),
   OrthogonalN(1), OrthogonalS(1), OrthogonalE(1), OrthogonalW(1),
   // Initialize mesh update flags to OFF (i.e. no update scheduled)
   InteriorMeshUpdate(OFF), GhostCellsUpdate(OFF), CornerGhostCellsUpdate(OFF),
   // Initialize state trackers
   InteriorCellGeometryStateTracker(0), GhostCellGeometryStateTracker(0), CornerGhostCellGeometryStateTracker(0),
   NumGQP(0)
{
  int Ni, Nj;
  int i,j;

  // allocate memory for the new container
  Ni = G.NCi - 2*G.Nghost;
  Nj = G.NCj - 2*G.Nghost;
  allocate(Ni,Nj,G.Nghost,G.MaxRecOrder());

  // Set the grid values by copying from the grid block G.
  if (G.Node != NULL) {

    // Copy the node locations of the grid block G.
    for (j  = G.JNl-G.Nghost; j <= G.JNu+G.Nghost ; ++j ) {
      for ( i = G.INl-G.Nghost ; i <= G.INu+G.Nghost ; ++i ) {
	Node[i][j].X = G.Node[i][j].X;
      } /* endfor */
    } /* endfor */
    
    // Copy the cell values of grid block G.
    for ( j = G.JCl-G.Nghost; j <= G.JCu+G.Nghost ; ++j) {
      for ( i = G.ICl-G.Nghost ; i <= G.ICu+G.Nghost ; ++i) {
	Cell[i][j] = G.Cell[i][j];
      } /* endfor */
    } /* endfor */

    // Copy the boundary condition type info of grid block G.
    for ( i = G.ICl-G.Nghost ; i <= G.ICu+G.Nghost ; ++i) {
      BCtypeN[i] = G.BCtypeN[i];
      BCtypeS[i] = G.BCtypeS[i];
    } /* endfor */
    for ( j = G.JCl-G.Nghost ; j <= G.JCu+G.Nghost ; ++j) {
      BCtypeE[j] = G.BCtypeE[j];
      BCtypeW[j] = G.BCtypeW[j];
    } /* endfor */
  } /* endif */


  // Copy boundary spline info of grid block G.
  if (G.BndNorthSpline.np != 0) {
    BndNorthSpline = G.BndNorthSpline;
    // Copy Info
    if (G.BndNorthSplineInfo != NULL){
      // allocate memory for BndNorthSplineInfo
      BndNorthSplineInfo = new Spline2DInterval_HO [G.NCi];
      // copy values
      for ( i = G.ICl; i <= G.ICu; ++i) {
	BndNorthSplineInfo[i] = G.BndNorthSplineInfo[i];
      }
    }
  } else if (BndNorthSpline.np != 0) {
    BndNorthSpline.deallocate();
  } /* endif */

  if (G.BndSouthSpline.np != 0) {
    BndSouthSpline = G.BndSouthSpline;
    // Copy Info
    if (G.BndSouthSplineInfo != NULL){
      // allocate memory for BndSouthSplineInfo
      BndSouthSplineInfo = new Spline2DInterval_HO [G.NCi];
      // copy values
      for ( i = G.ICl; i <= G.ICu; ++i) {
	BndSouthSplineInfo[i] = G.BndSouthSplineInfo[i];
      }
    }
  } else if (BndSouthSpline.np != 0) {
    BndSouthSpline.deallocate();
  } /* endif */
  
  if (G.BndEastSpline.np != 0) {
    BndEastSpline = G.BndEastSpline;
    // Copy Info
    if (G.BndEastSplineInfo != NULL){
      // allocate memory for BndEastSplineInfo
      BndEastSplineInfo = new Spline2DInterval_HO [G.NCj];
      // copy values
      for ( j = G.JCl; j <= G.JCu; ++j) {
	BndEastSplineInfo[j] = G.BndEastSplineInfo[j];
      }
    }
  } else if (BndEastSpline.np != 0) {
    BndEastSpline.deallocate();
  } /* endif */
  
  if (G.BndWestSpline.np != 0) {
    BndWestSpline = G.BndWestSpline;
    // Copy Info
    if (G.BndWestSplineInfo != NULL){
      // allocate memory for BndWestSplineInfo
      BndWestSplineInfo = new Spline2DInterval_HO [G.NCj];
      // copy values
      for ( j = G.JCl; j <= G.JCu; ++j) {
	BndWestSplineInfo[j] = G.BndWestSplineInfo[j];
      }
    }
  } else if (BndWestSpline.np != 0) {
    BndWestSpline.deallocate();
  } /* endif */

  // Copy the extensions to boundary splines
  if (G.ExtendWest_BndNorthSpline.np != 0) {
    ExtendWest_BndNorthSpline = G.ExtendWest_BndNorthSpline;
    // Copy Info
    if (G.ExtendWest_BndNorthSplineInfo != NULL){
      // allocate memory for ExtendWest_BndNorthSplineInfo
      ExtendWest_BndNorthSplineInfo = new Spline2DInterval_HO [G.Nghost];
      // copy values
      for ( i = 0; i < G.Nghost; ++i) {
	ExtendWest_BndNorthSplineInfo[i] = G.ExtendWest_BndNorthSplineInfo[i];
      }
    }
  } else if (ExtendWest_BndNorthSpline.np != 0) {
    ExtendWest_BndNorthSpline.deallocate();
  } /* endif */
  if (G.ExtendEast_BndNorthSpline.np != 0) {
    ExtendEast_BndNorthSpline = G.ExtendEast_BndNorthSpline;
    // Copy Info
    if (G.ExtendEast_BndNorthSplineInfo != NULL){
      // allocate memory for ExtendEast_BndNorthSplineInfo
      ExtendEast_BndNorthSplineInfo = new Spline2DInterval_HO [G.Nghost];
      // copy values
      for ( i = 0; i < G.Nghost; ++i) {
	ExtendEast_BndNorthSplineInfo[i] = G.ExtendEast_BndNorthSplineInfo[i];
      }
    }
  } else if (ExtendEast_BndNorthSpline.np != 0) {
    ExtendEast_BndNorthSpline.deallocate();
  } /* endif */

  if (G.ExtendWest_BndSouthSpline.np != 0) {
    ExtendWest_BndSouthSpline = G.ExtendWest_BndSouthSpline;
    // Copy Info
    if (G.ExtendWest_BndSouthSplineInfo != NULL){
      // allocate memory for ExtendWest_BndSouthSplineInfo
      ExtendWest_BndSouthSplineInfo = new Spline2DInterval_HO [G.Nghost];
      // copy values
      for ( i = 0; i < G.Nghost; ++i) {
	ExtendWest_BndSouthSplineInfo[i] = G.ExtendWest_BndSouthSplineInfo[i];
      }
    }
  } else if (ExtendWest_BndSouthSpline.np != 0) {
    ExtendWest_BndSouthSpline.deallocate();
  } /* endif */
  if (G.ExtendEast_BndSouthSpline.np != 0) {
    ExtendEast_BndSouthSpline = G.ExtendEast_BndSouthSpline;
    // Copy Info
    if (G.ExtendEast_BndSouthSplineInfo != NULL){
      // allocate memory for ExtendEast_BndSouthSplineInfo
      ExtendEast_BndSouthSplineInfo = new Spline2DInterval_HO [G.Nghost];
      // copy values
      for ( i = 0; i < G.Nghost; ++i) {
	ExtendEast_BndSouthSplineInfo[i] = G.ExtendEast_BndSouthSplineInfo[i];
      }
    }
  } else if (ExtendEast_BndSouthSpline.np != 0) {
    ExtendEast_BndSouthSpline.deallocate();
  } /* endif */

  if (G.ExtendNorth_BndEastSpline.np != 0) {
    ExtendNorth_BndEastSpline = G.ExtendNorth_BndEastSpline;
    // Copy Info
    if (G.ExtendNorth_BndEastSplineInfo != NULL){
      // allocate memory for ExtendNorth_BndEastSplineInfo
      ExtendNorth_BndEastSplineInfo = new Spline2DInterval_HO [G.Nghost];
      // copy values
      for ( j = 0; j < G.Nghost; ++j) {
	ExtendNorth_BndEastSplineInfo[j] = G.ExtendNorth_BndEastSplineInfo[j];
      }
    }
  } else if (ExtendNorth_BndEastSpline.np != 0) {
    ExtendNorth_BndEastSpline.deallocate();
  } /* endif */
  if (G.ExtendSouth_BndEastSpline.np != 0) {
    ExtendSouth_BndEastSpline = G.ExtendSouth_BndEastSpline;
    // Copy Info
    if (G.ExtendSouth_BndEastSplineInfo != NULL){
      // allocate memory for ExtendSouth_BndEastSplineInfo
      ExtendSouth_BndEastSplineInfo = new Spline2DInterval_HO [G.Nghost];
      // copy values
      for ( j = 0; j < G.Nghost; ++j) {
	ExtendSouth_BndEastSplineInfo[j] = G.ExtendSouth_BndEastSplineInfo[j];
      }
    }
  } else if (ExtendSouth_BndEastSpline.np != 0) {
    ExtendSouth_BndEastSpline.deallocate();
  } /* endif */
  
  if (G.ExtendNorth_BndWestSpline.np != 0) {
    ExtendNorth_BndWestSpline = G.ExtendNorth_BndWestSpline;
    // Copy Info
    if (G.ExtendNorth_BndWestSplineInfo != NULL){
      // allocate memory for ExtendNorth_BndWestSplineInfo
      ExtendNorth_BndWestSplineInfo = new Spline2DInterval_HO [G.Nghost];
      // copy values
      for ( j = 0; j < G.Nghost; ++j) {
	ExtendNorth_BndWestSplineInfo[j] = G.ExtendNorth_BndWestSplineInfo[j];
      }
    }
  } else if (ExtendNorth_BndWestSpline.np != 0) {
    ExtendNorth_BndWestSpline.deallocate();
  } /* endif */
  if (G.ExtendSouth_BndWestSpline.np != 0) {
    ExtendSouth_BndWestSpline = G.ExtendSouth_BndWestSpline;
    // Copy Info
    if (G.ExtendSouth_BndWestSplineInfo != NULL){
      // allocate memory for ExtendSouth_BndWestSplineInfo
      ExtendSouth_BndWestSplineInfo = new Spline2DInterval_HO [G.Nghost];
      // copy values
      for ( j = 0; j < G.Nghost; ++j) {
	ExtendSouth_BndWestSplineInfo[j] = G.ExtendSouth_BndWestSplineInfo[j];
      }
    }
  } else if (ExtendSouth_BndWestSpline.np != 0) {
    ExtendSouth_BndWestSpline.deallocate();
  } /* endif */


  // Copy boundary spline pathlength info of grid block G.
  SminN = G.SminN;
  SmaxN = G.SmaxN;
  SminS = G.SminS;
  SmaxS = G.SmaxS;
  SminE = G.SminE;
  SmaxE = G.SmaxE;
  SminW = G.SminW;
  SmaxW = G.SmaxW;
  
  // Copy node stretching info of grid block G.
  StretchI = G.StretchI;
  BetaI = G.BetaI;
  TauI = G.TauI;
  StretchJ = G.StretchJ;
  BetaJ = G.BetaJ;
  TauJ = G.TauJ;
  OrthogonalN = G.OrthogonalN;
  OrthogonalS = G.OrthogonalS;
  OrthogonalE = G.OrthogonalE;
  OrthogonalW = G.OrthogonalW;

  // Copy NumGQP
  SetNumberOfGaussQuadraturePoints(G.getNumGQP());

  /* No mesh update is required for this operation
     since all geometric values are copied directly. */
  Reset_Mesh_Update_Flags();
}

/*!
 * Allocate memory.
 *
 * \param Ni number of cells in i-direction
 * \param Nj number of cells in j-direction
 * \param Ng number of ghost cells
 * \param HighestRecOrder the highest reconstruction order used during the computation.
 */
void Grid2D_Quad_Block_HO::allocate(const int &Ni, const int &Nj, const int &Ng, const int &HighestRecOrder) {
  int i,j;

  // Check conditions
  assert( Ni > 1 && Nj > 1 && Ng >= 2 && HighestRecOrder >= 0);

  // Check if the new required memory has dimensions different than the currently allocated ones
  if ( Ni != (NCi-2*Nghost) || Nj != (NCj-2*Nghost) || Ng != Nghost ){

    // free the memory if there is memory allocated
    deallocate();
    
    // allocate new memory 
    NNi = Ni+2*Ng+1; INl = Ng; INu = Ni+Ng;
    NNj = Nj+2*Ng+1; JNl = Ng; JNu = Nj+Ng;
    NCi = Ni+2*Ng;   ICl = Ng; ICu = Ni+Ng-1;
    NCj = Nj+2*Ng;   JCl = Ng; JCu = Nj+Ng-1;
    Nghost = Ng;
    HighestReconstructionOrder = HighestRecOrder;
    SetNumberOfGaussQuadraturePoints();
    
    Node = new Node2D_HO*[NNi];
    for ( i = 0; i <= NNi-1 ; ++i ) Node[i] = new Node2D_HO[NNj];
    Cell = new Cell2D_HO*[NCi];
    for ( i = 0; i <= NCi-1 ; ++i ) Cell[i] = new Cell2D_HO[NCj];

    // allocate memory for the container of geometric moments
    for ( i = 0; i <NCi ; ++i ){
      for ( j = 0; j <NCj ; ++j ){
	Cell[i][j].SetGeomCoeffContainerSize(HighestReconstructionOrder);       
      }
    }

    BCtypeN = new int[NCi]; BCtypeS = new int[NCi];
    BCtypeE = new int[NCj]; BCtypeW = new int[NCj];

    // Complete memory allocation
    return;
    
  } else if (HighestRecOrder != HighestReconstructionOrder){
    // Check if the highest reconstruction order is different than the current one.
    HighestReconstructionOrder = HighestRecOrder;
    SetNumberOfGaussQuadraturePoints();
    
    // re-allocate memory for the container of geometric moments
    for ( i = 0; i <NCi ; ++i ){
      for ( j = 0; j <NCj ; ++j ){
	Cell[i][j].SetGeomCoeffContainerSize(HighestReconstructionOrder);       
      }
    }    
  }
}

/*
 * Deallocate memory.
 */
void Grid2D_Quad_Block_HO::deallocate(void) {
  int i;

  // Deallocate nodes
  if (Node != NULL){
    for ( i = 0; i <= NNi-1 ; ++i ) {
      delete []Node[i]; Node[i] = NULL;
    } /* endfor */
    delete []Node; Node = NULL;
  }

  // Deallocate cells
  if (Cell != NULL){
    for ( i = 0; i <= NCi-1 ; ++i ) {
      delete []Cell[i]; Cell[i] = NULL;
    } /* endfor */
    delete []Cell; Cell = NULL; 
  }

  // Deallocate North boundary conditions
  if (BCtypeN != NULL){
    delete []BCtypeN; BCtypeN = NULL;
  }
  // Deallocate South boundary conditions
  if (BCtypeS != NULL){
    delete []BCtypeS; BCtypeS = NULL; 
  }
  // Deallocate East boundary conditions
  if (BCtypeE != NULL){
    delete []BCtypeE; BCtypeE = NULL;
  }
  // Deallocate West boundary conditions
  if (BCtypeW != NULL){
    delete []BCtypeW; BCtypeW = NULL;
  }

  // Deallocate boundary splines
  BndNorthSpline.deallocate(); BndSouthSpline.deallocate();
  BndEastSpline.deallocate(); BndWestSpline.deallocate();

  ExtendWest_BndNorthSpline.deallocate(); ExtendEast_BndNorthSpline.deallocate();
  ExtendWest_BndSouthSpline.deallocate(); ExtendEast_BndSouthSpline.deallocate();
  ExtendNorth_BndEastSpline.deallocate(); ExtendSouth_BndEastSpline.deallocate();
  ExtendNorth_BndWestSpline.deallocate(); ExtendSouth_BndWestSpline.deallocate();

  // Deallocate boundary spline info
  deallocateBndSplineInfo();

  // Reset mesh indexes
  NNi = 0; INl = 0; INu = 0; NNj = 0; JNl = 0; JNu = 0;
  NCi = 0; ICl = 0; ICu = 0; NCj = 0; JCl = 0; JCu = 0;
  Nghost = 0;
  HighestReconstructionOrder = 0;
  StretchI = 0; StretchJ = 0; BetaI = ONE; TauI = ONE;
  BetaJ = ONE; TauJ = ONE;
  OrthogonalN = 1; OrthogonalS = 1; OrthogonalE = 1; OrthogonalW = 1;
  NumGQP = 0;

  // Reset state trackers
  InteriorMeshUpdate = OFF;
  GhostCellsUpdate = OFF;
  CornerGhostCellsUpdate = OFF;

  InteriorCellGeometryStateTracker = 0;
  GhostCellGeometryStateTracker = 0;
  CornerGhostCellGeometryStateTracker = 0;
}

/*!
 * Assignment operator =
 */
Grid2D_Quad_Block_HO& Grid2D_Quad_Block_HO::operator=(const Grid2D_Quad_Block_HO &Grid) {

  // !!! If the LHS grid block already has objects assigned, these are going to be deleted.
  // Handle self-assignment:
  if (this == & Grid) return *this;

  int Ni, Nj;
  int i,j;

  // re-allocate memory if there isn't enough
  Ni = Grid.NCi - 2*Grid.Nghost;
  Nj = Grid.NCj - 2*Grid.Nghost;
  allocate(Ni,Nj,Grid.Nghost,Grid.MaxRecOrder());

  // Set the grid values by copying from the grid block Grid.
  if (Grid.Node != NULL) {

    // Copy the node locations of the grid block Grid.
    for (j  = Grid.JNl-Grid.Nghost; j <= Grid.JNu+Grid.Nghost ; ++j ) {
      for ( i = Grid.INl-Grid.Nghost ; i <= Grid.INu+Grid.Nghost ; ++i ) {
	Node[i][j].X = Grid.Node[i][j].X;
      } /* endfor */
    } /* endfor */
    
    // Copy the cell values of grid block Grid.
    for ( j = Grid.JCl-Grid.Nghost; j <= Grid.JCu+Grid.Nghost ; ++j) {
      for ( i = Grid.ICl-Grid.Nghost ; i <= Grid.ICu+Grid.Nghost ; ++i) {
	Cell[i][j] = Grid.Cell[i][j];
      } /* endfor */
    } /* endfor */

    // Copy the boundary condition type info of grid block Grid.
    for ( i = Grid.ICl-Grid.Nghost ; i <= Grid.ICu+Grid.Nghost ; ++i) {
      BCtypeN[i] = Grid.BCtypeN[i];
      BCtypeS[i] = Grid.BCtypeS[i];
    } /* endfor */
    for ( j = Grid.JCl-Grid.Nghost ; j <= Grid.JCu+Grid.Nghost ; ++j) {
      BCtypeE[j] = Grid.BCtypeE[j];
      BCtypeW[j] = Grid.BCtypeW[j];
    } /* endfor */
  } /* endif */


  // Copy boundary spline info of grid block Grid.
  if (Grid.BndNorthSpline.np != 0) {
    BndNorthSpline = Grid.BndNorthSpline;
    // Copy associated Info
    if (Grid.BndNorthSplineInfo != NULL){
      // allocate memory for BndNorthSplineInfo if it hasn't been allocated
      if (BndNorthSplineInfo == NULL){      
	BndNorthSplineInfo = new Spline2DInterval_HO [Grid.NCi];
      }
      // copy values
      for ( i = Grid.ICl; i <= Grid.ICu; ++i) {
	BndNorthSplineInfo[i] = Grid.BndNorthSplineInfo[i];
      }
    } else {
      deallocate_BndNorthSplineInfo();
    }
  } else if (BndNorthSpline.np != 0) {
    BndNorthSpline.deallocate();
    deallocate_BndNorthSplineInfo();
  } /* endif */

  if (Grid.BndSouthSpline.np != 0) {
    BndSouthSpline = Grid.BndSouthSpline;
    // Copy associated Info
    if (Grid.BndSouthSplineInfo != NULL){
      // allocate memory for BndSouthSplineInfo if it hasn't been allocated
      if (BndSouthSplineInfo == NULL){      
	BndSouthSplineInfo = new Spline2DInterval_HO [Grid.NCi];
      }
      // copy values
      for ( i = Grid.ICl; i <= Grid.ICu; ++i) {
	BndSouthSplineInfo[i] = Grid.BndSouthSplineInfo[i];
      }
    } else {
      deallocate_BndSouthSplineInfo();
    }
  } else if (BndSouthSpline.np != 0) {
    BndSouthSpline.deallocate();
    deallocate_BndSouthSplineInfo();
  } /* endif */
  
  if (Grid.BndEastSpline.np != 0) {
    BndEastSpline = Grid.BndEastSpline;
    // Copy associated Info
    if (Grid.BndEastSplineInfo != NULL){
      // allocate memory for BndEastSplineInfo if it hasn't been allocated
      if (BndEastSplineInfo == NULL){      
	BndEastSplineInfo = new Spline2DInterval_HO [Grid.NCj];
      }
      // copy values
      for ( j = Grid.JCl; j <= Grid.JCu; ++j) {
	BndEastSplineInfo[j] = Grid.BndEastSplineInfo[j];
      }
    } else {
      deallocate_BndEastSplineInfo();
    }
  } else if (BndEastSpline.np != 0) {
    BndEastSpline.deallocate();
    deallocate_BndEastSplineInfo();
  } /* endif */
  
  if (Grid.BndWestSpline.np != 0) {
    BndWestSpline = Grid.BndWestSpline;
    // Copy associated Info
    if (Grid.BndWestSplineInfo != NULL){
      // allocate memory for BndWestSplineInfo if it hasn't been allocated
      if (BndWestSplineInfo == NULL){      
	BndWestSplineInfo = new Spline2DInterval_HO [Grid.NCj];
      }
      // copy values
      for ( j = Grid.JCl; j <= Grid.JCu; ++j) {
	BndWestSplineInfo[j] = Grid.BndWestSplineInfo[j];
      }
    } else {
      deallocate_BndWestSplineInfo();
    }
  } else if (BndWestSpline.np != 0) {
    BndWestSpline.deallocate();
    deallocate_BndWestSplineInfo();
  } /* endif */

  // Copy the extensions to boundary splines
  if (Grid.ExtendWest_BndNorthSpline.np != 0) {
    ExtendWest_BndNorthSpline = Grid.ExtendWest_BndNorthSpline;
    // Copy Info
    if (Grid.ExtendWest_BndNorthSplineInfo != NULL){
      // allocate memory for ExtendWest_BndNorthSplineInfo
      ExtendWest_BndNorthSplineInfo = new Spline2DInterval_HO [Grid.Nghost];
      // copy values
      for ( i = 0; i < Grid.Nghost; ++i) {
	ExtendWest_BndNorthSplineInfo[i] = Grid.ExtendWest_BndNorthSplineInfo[i];
      }
    }    
  } else if (ExtendWest_BndNorthSpline.np != 0) {
    ExtendWest_BndNorthSpline.deallocate();
    deallocate_ExtendWest_BndNorthSplineInfo();
  } /* endif */
  if (Grid.ExtendEast_BndNorthSpline.np != 0) {
    ExtendEast_BndNorthSpline = Grid.ExtendEast_BndNorthSpline;
    // Copy Info
    if (Grid.ExtendEast_BndNorthSplineInfo != NULL){
      // allocate memory for ExtendEast_BndNorthSplineInfo
      ExtendEast_BndNorthSplineInfo = new Spline2DInterval_HO [Grid.Nghost];
      // copy values
      for ( i = 0; i < Grid.Nghost; ++i) {
	ExtendEast_BndNorthSplineInfo[i] = Grid.ExtendEast_BndNorthSplineInfo[i];
      }
    }
  } else if (ExtendEast_BndNorthSpline.np != 0) {
    ExtendEast_BndNorthSpline.deallocate();
    deallocate_ExtendEast_BndNorthSplineInfo();
  } /* endif */

  if (Grid.ExtendWest_BndSouthSpline.np != 0) {
    ExtendWest_BndSouthSpline = Grid.ExtendWest_BndSouthSpline;
    // Copy Info
    if (Grid.ExtendWest_BndSouthSplineInfo != NULL){
      // allocate memory for ExtendWest_BndSouthSplineInfo
      ExtendWest_BndSouthSplineInfo = new Spline2DInterval_HO [Grid.Nghost];
      // copy values
      for ( i = 0; i < Grid.Nghost; ++i) {
	ExtendWest_BndSouthSplineInfo[i] = Grid.ExtendWest_BndSouthSplineInfo[i];
      }
    }
  } else if (ExtendWest_BndSouthSpline.np != 0) {
    ExtendWest_BndSouthSpline.deallocate();
    deallocate_ExtendWest_BndSouthSplineInfo();
  } /* endif */
  if (Grid.ExtendEast_BndSouthSpline.np != 0) {
    ExtendEast_BndSouthSpline = Grid.ExtendEast_BndSouthSpline;
    // Copy Info
    if (Grid.ExtendEast_BndSouthSplineInfo != NULL){
      // allocate memory for ExtendEast_BndSouthSplineInfo
      ExtendEast_BndSouthSplineInfo = new Spline2DInterval_HO [Grid.Nghost];
      // copy values
      for ( i = 0; i < Grid.Nghost; ++i) {
	ExtendEast_BndSouthSplineInfo[i] = Grid.ExtendEast_BndSouthSplineInfo[i];
      }
    }
  } else if (ExtendEast_BndSouthSpline.np != 0) {
    ExtendEast_BndSouthSpline.deallocate();
    deallocate_ExtendEast_BndSouthSplineInfo();
  } /* endif */

  if (Grid.ExtendNorth_BndEastSpline.np != 0) {
    ExtendNorth_BndEastSpline = Grid.ExtendNorth_BndEastSpline;
    // Copy Info
    if (Grid.ExtendNorth_BndEastSplineInfo != NULL){
      // allocate memory for ExtendNorth_BndEastSplineInfo
      ExtendNorth_BndEastSplineInfo = new Spline2DInterval_HO [Grid.Nghost];
      // copy values
      for ( j = 0; j < Grid.Nghost; ++j) {
	ExtendNorth_BndEastSplineInfo[j] = Grid.ExtendNorth_BndEastSplineInfo[j];
      }
    }
  } else if (ExtendNorth_BndEastSpline.np != 0) {
    ExtendNorth_BndEastSpline.deallocate();
    deallocate_ExtendNorth_BndEastSplineInfo();
  } /* endif */
  if (Grid.ExtendSouth_BndEastSpline.np != 0) {
    ExtendSouth_BndEastSpline = Grid.ExtendSouth_BndEastSpline;
    // Copy Info
    if (Grid.ExtendSouth_BndEastSplineInfo != NULL){
      // allocate memory for ExtendSouth_BndEastSplineInfo
      ExtendSouth_BndEastSplineInfo = new Spline2DInterval_HO [Grid.Nghost];
      // copy values
      for ( j = 0; j < Grid.Nghost; ++j) {
	ExtendSouth_BndEastSplineInfo[j] = Grid.ExtendSouth_BndEastSplineInfo[j];
      }
    }
  } else if (ExtendSouth_BndEastSpline.np != 0) {
    ExtendSouth_BndEastSpline.deallocate();
    deallocate_ExtendSouth_BndEastSplineInfo();
  } /* endif */
  
  if (Grid.ExtendNorth_BndWestSpline.np != 0) {
    ExtendNorth_BndWestSpline = Grid.ExtendNorth_BndWestSpline;
    // Copy Info
    if (Grid.ExtendNorth_BndWestSplineInfo != NULL){
      // allocate memory for ExtendNorth_BndWestSplineInfo
      ExtendNorth_BndWestSplineInfo = new Spline2DInterval_HO [Grid.Nghost];
      // copy values
      for ( j = 0; j < Grid.Nghost; ++j) {
	ExtendNorth_BndWestSplineInfo[j] = Grid.ExtendNorth_BndWestSplineInfo[j];
      }
    }
  } else if (ExtendNorth_BndWestSpline.np != 0) {
    ExtendNorth_BndWestSpline.deallocate();
    deallocate_ExtendNorth_BndWestSplineInfo();
  } /* endif */
  if (Grid.ExtendSouth_BndWestSpline.np != 0) {
    ExtendSouth_BndWestSpline = Grid.ExtendSouth_BndWestSpline;
    // Copy Info
    if (Grid.ExtendSouth_BndWestSplineInfo != NULL){
      // allocate memory for ExtendSouth_BndWestSplineInfo
      ExtendSouth_BndWestSplineInfo = new Spline2DInterval_HO [Grid.Nghost];
      // copy values
      for ( j = 0; j < Grid.Nghost; ++j) {
	ExtendSouth_BndWestSplineInfo[j] = Grid.ExtendSouth_BndWestSplineInfo[j];
      }
    }
  } else if (ExtendSouth_BndWestSpline.np != 0) {
    ExtendSouth_BndWestSpline.deallocate();
    deallocate_ExtendSouth_BndWestSplineInfo();
  } /* endif */


  // Copy boundary spline pathlength info of grid block Grid.
  SminN = Grid.SminN;
  SmaxN = Grid.SmaxN;
  SminS = Grid.SminS;
  SmaxS = Grid.SmaxS;
  SminE = Grid.SminE;
  SmaxE = Grid.SmaxE;
  SminW = Grid.SminW;
  SmaxW = Grid.SmaxW;
  
  // Copy node stretching info of grid block Grid.
  StretchI = Grid.StretchI;
  BetaI = Grid.BetaI;
  TauI = Grid.TauI;
  StretchJ = Grid.StretchJ;
  BetaJ = Grid.BetaJ;
  TauJ = Grid.TauJ;
  OrthogonalN = Grid.OrthogonalN;
  OrthogonalS = Grid.OrthogonalS;
  OrthogonalE = Grid.OrthogonalE;
  OrthogonalW = Grid.OrthogonalW;

  // Copy NumGQP
  SetNumberOfGaussQuadraturePoints(Grid.getNumGQP());

  /* No mesh update is required for this operation
     since all geometric values are copied directly. */
  Reset_Mesh_Update_Flags();

  // Mark the current geometry different than the previous one
  New_Global_Geometry_State();

  return *this;
}

/*!
 * Calculate influence coefficients 
 * based on bilinear interpolation.
 */
double Grid2D_Quad_Block_HO::dNdC(const int ii, const int jj, 
				  const int cell_orientation) const {
  double ax, bx, cx, dx, ay, by, cy, dy, aa, bb, cc, x, y, 
    eta1, zeta1, eta2, zeta2, eta, zeta;
  double interpolation_coefficient;
  x  = Node[ii][jj].X.x;
  y  = Node[ii][jj].X.y;
  ax = Cell[ii-1][jj-1].Xc.x;
  bx = Cell[ii-1][jj  ].Xc.x - Cell[ii-1][jj-1].Xc.x; 
  cx = Cell[ii  ][jj-1].Xc.x - Cell[ii-1][jj-1].Xc.x; 
  dx = Cell[ii  ][jj  ].Xc.x + Cell[ii-1][jj-1].Xc.x -
    Cell[ii-1][jj  ].Xc.x - Cell[ii  ][jj-1].Xc.x;
  ay = Cell[ii-1][jj-1].Xc.y;
  by = Cell[ii-1][jj  ].Xc.y - Cell[ii-1][jj-1].Xc.y; 
  cy = Cell[ii  ][jj-1].Xc.y - Cell[ii-1][jj-1].Xc.y; 
  dy = Cell[ii  ][jj  ].Xc.y + Cell[ii-1][jj-1].Xc.y -
    Cell[ii-1][jj  ].Xc.y - Cell[ii  ][jj-1].Xc.y;
  aa = bx*dy - dx*by;
  bb = dy*(ax-x) + bx*cy - cx*by+dx*(y-ay);
  cc = cy*(ax-x) + cx*(y-ay);
  if (fabs(aa) < TOLER*TOLER) {
    if (fabs(bb) >= TOLER*TOLER) {
      zeta1 = -cc/bb;
    } else { 
      zeta1 = -cc/sgn(bb)*(TOLER*TOLER); 
    } 
    if (fabs(cy+dy*zeta1) >= TOLER*TOLER) {
      eta1 = (y-ay-by*zeta1)/(cy+dy*zeta1); 
    } else { 
      eta1 = HALF;
    } 
    zeta2 = zeta1;
    eta2  = eta1;
  } else {
    if (bb*bb-FOUR*aa*cc >= TOLER*TOLER) { 
      zeta1 = HALF*(-bb+sqrt(bb*bb-FOUR*aa*cc))/aa; 
    } else { zeta1 = -HALF*bb/aa;
    } 
    if (fabs(cy+dy*zeta1) < TOLER*TOLER) {
      eta1 = -ONE;
    } else {
      eta1 = (y-ay-by*zeta1)/(cy+dy*zeta1);
    }
    if (bb*bb-FOUR*aa*cc >= TOLER*TOLER) {
      zeta2 = HALF*(-bb-sqrt(bb*bb-FOUR*aa*cc))/aa; 
    } else {
      zeta2 = -HALF*bb/aa;
    }
    if (fabs(cy+dy*zeta2) < TOLER*TOLER) { 
      eta2 = -ONE;
    } else {
      eta2 = (y-ay-by*zeta2)/(cy+dy*zeta2);
    }
  }
  if (zeta1 > -TOLER && zeta1 < ONE + TOLER &&
      eta1  > -TOLER && eta1  < ONE + TOLER) {
    zeta = zeta1;
    eta  = eta1;
  } else if (zeta2 > -TOLER && zeta2 < ONE + TOLER &&
	     eta2  > -TOLER && eta2  < ONE + TOLER) {
    zeta = zeta2;
    eta  = eta2;
  } else {
    zeta = HALF;
    eta  = HALF;
  }
  // Determine the coefficient or weight for cell orientation.
  switch(cell_orientation) {
  case NORTH_EAST : 
    interpolation_coefficient = zeta*eta;
    break;
  case NORTH_WEST : 
    interpolation_coefficient = zeta - zeta*eta;
    break;
  case SOUTH_EAST : 
    interpolation_coefficient = eta - zeta*eta;
    break;
  case SOUTH_WEST : 
    interpolation_coefficient = ONE - zeta - eta + zeta*eta;
    break;
  default:
    cout <<"\n Improper Orientation in Grid2D_Quad::dNdC!\n"; 
    interpolation_coefficient = -ONE;
    break;
  }
  return interpolation_coefficient;
}

/*!
 * Calculate influence coefficients 
 * based on bilinear interpolation. 
 */
//  dFacedC:
//  
//  We know that we can write each face variable, uf, as a geometric
//  weight of the cell-centred value, uc: uf = intp uc
//  
//  This function determine what this "intp" is.
//  
//  For diamond-path integration: (which is how we find the gradients of
//  the face variables but we are sloppy and reuse the term "diamond" also
//  to apply to how we find the variables themselves) to find the face
//  variables we first find the nodal values using an interpolation and
//  then use those nodal values together with the cell-centred variables
//  in a second (third, really, after two nodes) interpolation to find the
//  face variable.
//  
//  For example, suppose we wanted to find density at the north face of
//  cell (i, j). We would first find the value at node NE using the
//  values from the four neighbouring cells of node NE. Then we would find
//  the value at node NW. Then we would apply a third interpolation using
//  the values at NW, NE, the centre of (i, j) and the centre of (i, j+1).
//  
//  In a simplifying assumption, we ignore the third interpolation and
//  simply assume that the face value is the average of the node values. 
//  
//  In what follows, then, "left_node_weight" and "right_node_weight" 
//  are the weights used to obtain the nodal values. If the cell only 
//  contributes to one of the nodal values then one of "left_node_weight" 
//  or "right_node_weight" is zero. "intp" then is simply the average of 
//  "left_node_weight" and "right_node_weight".
//  
//  The "source" cell is in the direction "cell_orientation" from 
//  cell (i, j).
//  
//  For example, suppose that 
//
//  face = North 
//  cell_orientation = East
//
//  Now we wish to know how the cell to the east of (i, j) influences
//  the value on the north face of (i, j). In this case, the cell to
//  the east influences the value at the node to the right of the face
//  value but does not influence the value at the node to the left of the
//  face.  So "left_node_weight" is zero. Also the node to the right of 
//  the face is the NW node of the cell to the east.
//
void Grid2D_Quad_Block_HO::dFacedC(double &left_node_weight, 
				   double &right_node_weight, 
				   const int i, const int j, const int face, 
				   const int cell_orientation) const {
  int lnodei = 0, lnodej = 0; //  left node (not cell) index
  int rnodei = 0, rnodej = 0; // right node (not cell) index

  switch (face){
  case NORTH:
    lnodei = i  ; lnodej = j+1;
    rnodei = i+1; rnodej = j+1;
    if (cell_orientation == CENTER) {    
      left_node_weight = dNdC(lnodei, lnodej, SOUTH_EAST);
      right_node_weight = dNdC(rnodei, rnodej, SOUTH_WEST);
    } else if (cell_orientation == NORTH) {     
      left_node_weight = dNdC(lnodei, lnodej, NORTH_EAST);
      right_node_weight = dNdC(rnodei, rnodej, NORTH_WEST);
    } else if (cell_orientation == EAST) {
      left_node_weight = ZERO;
      right_node_weight = dNdC(rnodei, rnodej, SOUTH_EAST);
    } else if (cell_orientation == WEST){ 
      left_node_weight = dNdC(lnodei, lnodej, SOUTH_WEST);      
      right_node_weight = ZERO;                   
    } else if (cell_orientation == NORTH_WEST) {
      left_node_weight = dNdC(lnodei, lnodej, NORTH_WEST);      
      right_node_weight = ZERO;   
    } else if (cell_orientation == NORTH_EAST) {
      left_node_weight = ZERO;
      right_node_weight = dNdC(rnodei, rnodej, NORTH_EAST);    
    } /* endif */
    break;
      
  case EAST:
    lnodei = i+1; lnodej = j+1;
    rnodei = i+1; rnodej = j  ;
    if (cell_orientation == CENTER) {  
      left_node_weight = dNdC(lnodei, lnodej, SOUTH_WEST);
      right_node_weight = dNdC(rnodei, rnodej, NORTH_WEST); 
    } else if (cell_orientation == EAST) {    
      left_node_weight = dNdC(lnodei, lnodej, SOUTH_EAST);
      right_node_weight = dNdC(rnodei, rnodej, NORTH_EAST); 
    } else if (cell_orientation == NORTH) {     
      left_node_weight = dNdC(lnodei, lnodej, NORTH_WEST);
      right_node_weight = ZERO;
    } else if (cell_orientation == SOUTH) {           
      left_node_weight = ZERO;
      right_node_weight = dNdC(rnodei, rnodej, SOUTH_WEST);       
    } else if (cell_orientation == NORTH_EAST) {     
      left_node_weight = dNdC(lnodei, lnodej, NORTH_EAST);
      right_node_weight = ZERO;
    } else if (cell_orientation == SOUTH_EAST) {           
      left_node_weight = ZERO;
      right_node_weight = dNdC(rnodei, rnodej, SOUTH_EAST);       
    } /* endif */
    break;
  
  case SOUTH:
    lnodei = i+1; lnodej = j;
    rnodei = i  ; rnodej = j;
    if (cell_orientation == CENTER){ 
      left_node_weight = dNdC(lnodei, lnodej, NORTH_WEST);
      right_node_weight = dNdC(rnodei, rnodej, NORTH_EAST);
    } else if (cell_orientation == SOUTH) {     
      left_node_weight = dNdC(lnodei, lnodej, SOUTH_WEST);
      right_node_weight = dNdC(rnodei, rnodej, SOUTH_EAST);
    } else if (cell_orientation == EAST) { 
      left_node_weight = dNdC(lnodei, lnodej, NORTH_EAST);      
      right_node_weight = ZERO;
    } else if (cell_orientation == WEST){    
      left_node_weight = ZERO; 
      right_node_weight = dNdC(rnodei, rnodej, NORTH_WEST);                              
    } else if (cell_orientation == SOUTH_EAST) { 
      left_node_weight = dNdC(lnodei, lnodej, SOUTH_EAST);      
      right_node_weight = ZERO;
    } else if (cell_orientation == SOUTH_WEST) {    
      left_node_weight = ZERO; 
      right_node_weight = dNdC(rnodei, rnodej, SOUTH_WEST);      
    } /* endif */
    break;
   
  case WEST:
    lnodei = i; lnodej = j  ;
    rnodei = i; rnodej = j+1;
    if (cell_orientation == CENTER) {  
      left_node_weight = dNdC(lnodei, lnodej, NORTH_EAST);
      right_node_weight = dNdC(rnodei, rnodej, SOUTH_EAST);  
    } else if (cell_orientation == WEST) {  
      left_node_weight = dNdC(lnodei, lnodej, NORTH_WEST);
      right_node_weight = dNdC(rnodei, rnodej, SOUTH_WEST);  
    } else if (cell_orientation == SOUTH) {           
      left_node_weight = dNdC(lnodei, lnodej, SOUTH_EAST);
      right_node_weight = ZERO;
    } else if (cell_orientation == NORTH) {   
      left_node_weight = ZERO;
      right_node_weight = dNdC(rnodei, rnodej, NORTH_EAST);
    } else if (cell_orientation == SOUTH_WEST) {           
      left_node_weight = dNdC(lnodei, lnodej, SOUTH_WEST);
      right_node_weight = ZERO;
    } else if (cell_orientation == NORTH_WEST) {   
      left_node_weight = ZERO;
      right_node_weight = dNdC(rnodei, rnodej, NORTH_WEST);
    } /* endif */
    break;     
  }
  if (left_node_weight < ZERO || right_node_weight < ZERO) {
    cout <<"\n Invalid node weightings in Grid2D_Quad::dFacedC!\n"; 
  }
}

/*!
 * Influence coefficients for diamond path gradient
 * reconstruction of face solution gradients.
 */
//  dDiamondPathdC:
//  
//  Finds the weights used to determine the gradient of a variable on the
//  cell face from a cell-centred value. That is, if ux is the gradient of
//  velocity at a cell face and u is the velocity at one cell centre then
//  we approximate ux by:
//  
//    ux = ( d_dWdx_dW )( u ) + contributions from other cells
//  
//  This function finds d_dWdx_dW and the corresponding d_dWdy_dW. 
//
void Grid2D_Quad_Block_HO::dDiamondPathdC(double &d_dWdx_dW, double &d_dWdy_dW, 
					  const int i, const int j, const int face, 
					  const double &face_left_node_weight, 
					  const double &face_right_node_weight, 
					  const int cell_orientation) const {
  double AREA = 0.0;
  Vector2D norm[4];
  double  dWnNWdWc = 0.0, dWnNEdWc = 0.0,  dWnSWdWc = 0.0, dWnSEdWc = 0.0;
 
  switch(face){
    /*************** NORTH ****************************/
  case NORTH: 
    dWnNWdWc = face_left_node_weight;
    dWnNEdWc = face_right_node_weight;
    //  normal vector of the SE side of a diamond 
    norm[0].x = nodeNE(i,j).X.y - Cell[i][j].Xc.y;
    norm[0].y = -( nodeNE(i,j).X.x - Cell[i][j].Xc.x);
    //  normal vector of the NE side of a diamond 
    norm[1].x = Cell[i][j+1].Xc.y - nodeNE(i,j).X.y;
    norm[1].y = -(Cell[i][j+1].Xc.x - nodeNE(i,j).X.x);
    //  normal vector of the NW side of a diamond 
    norm[2].x =  nodeNW(i,j).X.y - Cell[i][j+1].Xc.y;
    norm[2].y = -(nodeNW(i,j).X.x - Cell[i][j+1].Xc.x);
    //  normal vector of the SW side of a diamond 
    norm[3].x = Cell[i][j].Xc.y - nodeNW(i,j).X.y ;
    norm[3].y = -( Cell[i][j].Xc.x - nodeNW(i,j).X.x);
    AREA =  HALF*(fabs((nodeNE(i,j).X-Cell[i][j].Xc)^
		       (nodeNW(i,j).X-Cell[i][j].Xc)) +
		  fabs((nodeNW(i,j).X-Cell[i][j+1].Xc)^
		       (nodeNE(i,j).X-Cell[i][j+1].Xc)));
    if(cell_orientation == CENTER){
      d_dWdx_dW = HALF*((ONE+dWnNEdWc)*norm[0].x+dWnNEdWc* norm[1].x+ 
			dWnNWdWc* norm[2].x+ (ONE+dWnNWdWc)* norm[3].x)/AREA;
      d_dWdy_dW = HALF*((ONE+dWnNEdWc)*norm[0].y+dWnNEdWc*norm[1].y+ 
			dWnNWdWc* norm[2].y+ (ONE+dWnNWdWc)* norm[3].y)/AREA;  
    } else if (cell_orientation == NORTH) {
      d_dWdx_dW = HALF*(dWnNEdWc*norm[0].x + (ONE+dWnNEdWc)* norm[1].x + 
			(ONE+dWnNWdWc)* norm[2].x + dWnNWdWc*norm[3].x)/AREA;
      d_dWdy_dW = HALF*(dWnNEdWc*norm[0].y + (ONE+dWnNEdWc)* norm[1].y + 
			(ONE+dWnNWdWc)* norm[2].y + dWnNWdWc*norm[3].y)/AREA;  
    } else if (cell_orientation == EAST || cell_orientation == NORTH_EAST) {
      d_dWdx_dW = HALF*(dWnNEdWc*(norm[0].x+norm[1].x))/AREA;
      d_dWdy_dW = HALF*(dWnNEdWc*(norm[0].y+norm[1].y))/AREA;  
    } else if (cell_orientation == WEST || cell_orientation == NORTH_WEST) {
      d_dWdx_dW = HALF*(dWnNWdWc*(norm[2].x+norm[3].x))/AREA;
      d_dWdy_dW = HALF*(dWnNWdWc*(norm[2].y+norm[3].y))/AREA;  
    } /* endif */
    break;
    
    /*************** EAST ****************************/
  case EAST:
    dWnNEdWc = face_left_node_weight;
    dWnSEdWc = face_right_node_weight; 
    //  normal vector of the SE side of a diamond 
    norm[0].x =  Cell[i+1][j].Xc.y - nodeSE(i,j).X.y;
    norm[0].y = -(Cell[i+1][j].Xc.x - nodeSE(i,j).X.x);
    //  normal vector of the NE side of a diamond 
    norm[1].x = nodeNE(i,j).X.y -  Cell[i+1][j].Xc.y ;
    norm[1].y = -(nodeNE(i,j).X.x -  Cell[i+1][j].Xc.x );
    //  normal vector of the NW side of a diamond 
    norm[2].x =   Cell[i][j].Xc.y - nodeNE(i,j).X.y ;
    norm[2].y = -(Cell[i][j].Xc.x - nodeNE(i,j).X.x);
    //  normal vector of the SW side of a diamond 
    norm[3].x = nodeSE(i,j).X.y - Cell[i][j].Xc.y;
    norm[3].y = -(nodeSE(i,j).X.x - Cell[i][j].Xc.x);
    AREA =  HALF*(fabs((nodeNE(i,j).X-Cell[i+1][j].Xc)^
		       (nodeSE(i,j).X-Cell[i+1][j].Xc)) +
		  fabs((nodeSE(i,j).X-Cell[i][j].Xc)^
		       (nodeNE(i,j).X-Cell[i][j].Xc)));
    if(cell_orientation == CENTER){
      d_dWdx_dW = HALF*(dWnSEdWc*norm[0].x+ dWnNEdWc*norm[1].x+ 
			(ONE+ dWnNEdWc)*norm[2].x+ (ONE+dWnSEdWc)* norm[3].x)/AREA;
      d_dWdy_dW = HALF*(dWnSEdWc*norm[0].y+ dWnNEdWc* norm[1].y+ 
			(ONE+ dWnNEdWc)*norm[2].y+ (ONE+dWnSEdWc)* norm[3].y)/AREA;  
    } else if (cell_orientation == EAST) {
      d_dWdx_dW = HALF*((ONE+dWnSEdWc)*norm[0].x + (ONE+dWnNEdWc)*norm[1].x + 
			dWnNEdWc* norm[2].x + dWnSEdWc*norm[3].x)/AREA;
      d_dWdy_dW = HALF*((ONE+dWnSEdWc)* norm[0].y + (ONE+dWnNEdWc)*norm[1].y + 
			dWnNEdWc*norm[2].y + dWnSEdWc*norm[3].y)/AREA;  
    } else if (cell_orientation == NORTH || cell_orientation == NORTH_EAST) {
      d_dWdx_dW = HALF*(dWnNEdWc*(norm[1].x+norm[2].x))/AREA;
      d_dWdy_dW = HALF*(dWnNEdWc*(norm[1].y+norm[2].y))/AREA;  
    } else if (cell_orientation == SOUTH || cell_orientation == SOUTH_EAST) {
      d_dWdx_dW = HALF*(dWnSEdWc*(norm[0].x+norm[3].x))/AREA;
      d_dWdy_dW = HALF*(dWnSEdWc*(norm[0].y+norm[3].y))/AREA;  
    } /* endif */
    break;

    /*************** SOUTH ****************************/
  case SOUTH:
    dWnSEdWc = face_left_node_weight;
    dWnSWdWc = face_right_node_weight;
    //  normal vector of the SE side of a diamond 
    norm[0].x =  nodeSE(i,j).X.y - Cell[i][j-1].Xc.y;
    norm[0].y = -(nodeSE(i,j).X.x - Cell[i][j-1].Xc.x  );
    //  normal vector of the NE side of a diamond 
    norm[1].x = Cell[i][j].Xc.y -  nodeSE(i,j).X.y;
    norm[1].y = -(Cell[i][j].Xc.x -  nodeSE(i,j).X.x);
    //  normal vector of the NW side of a diamond 
    norm[2].x =   nodeSW(i,j).X.y - Cell[i][j].Xc.y ;
    norm[2].y = -(nodeSW(i,j).X.x - Cell[i][j].Xc.x);
    //  normal vector of the SW side of a diamond 
    norm[3].x =  Cell[i][j-1].Xc.y - nodeSW(i,j).X.y;
    norm[3].y = -(Cell[i][j-1].Xc.x- nodeSW(i,j).X.x);
    AREA =  HALF*(fabs((nodeSE(i,j).X-Cell[i][j-1].Xc)^
		       (nodeSW(i,j).X-Cell[i][j-1].Xc)) +
		  fabs((nodeSE(i,j).X-Cell[i][j].Xc)^
		       (nodeSW(i,j).X-Cell[i][j].Xc)));
    if(cell_orientation == CENTER){
      d_dWdx_dW = HALF*(dWnSEdWc*norm[0].x+ (ONE+dWnSEdWc)*norm[1].x+ 
			(ONE+ dWnSWdWc)*norm[2].x+ (dWnSWdWc)*norm[3].x)/AREA;
      d_dWdy_dW = HALF*(dWnSEdWc*norm[0].y+ (ONE+dWnSEdWc)*norm[1].y+ 
			(ONE+ dWnSWdWc)*norm[2].y+ (dWnSWdWc)*norm[3].y)/AREA;  
    } else if( cell_orientation == SOUTH) {
      d_dWdx_dW = HALF*((ONE+dWnSEdWc)*norm[0].x + dWnSEdWc*norm[1].x + 
			dWnSWdWc*norm[2].x + (ONE+dWnSWdWc)* norm[3].x)/AREA;
      d_dWdy_dW = HALF*((ONE+dWnSEdWc)*norm[0].y + dWnSEdWc*norm[1].y + 
			dWnSWdWc*norm[2].y + (ONE+dWnSWdWc)* norm[3].y)/AREA;  
    } else if (cell_orientation == EAST || cell_orientation == SOUTH_EAST) {
      d_dWdx_dW = HALF*(dWnSEdWc*(norm[0].x+norm[1].x))/AREA;
      d_dWdy_dW = HALF*(dWnSEdWc*(norm[0].y+norm[1].y))/AREA;  
    } else if (cell_orientation == WEST || cell_orientation == SOUTH_WEST) {
      d_dWdx_dW = HALF*(dWnSWdWc*(norm[2].x+norm[3].x))/AREA;
      d_dWdy_dW = HALF*(dWnSWdWc*(norm[2].y+norm[3].y))/AREA;  
    } /* endif */
    break;
    
    /*************** WEST ****************************/
  case WEST:
    dWnSWdWc = face_left_node_weight;
    dWnNWdWc = face_right_node_weight;
    //  normal vector of the SE side of a diamond 
    norm[0].x =   Cell[i][j].Xc.y - nodeSW(i,j).X.y;
    norm[0].y = -(Cell[i][j].Xc.x - nodeSW(i,j).X.x);
    //  normal vector of the NE side of a diamond 
    norm[1].x =  nodeNW(i,j).X.y - Cell[i][j].Xc.y;
    norm[1].y = -(nodeNW(i,j).X.x - Cell[i][j].Xc.x);
    //  normal vector of the NW side of a diamond 
    norm[2].x =  Cell[i-1][j].Xc.y - nodeNW(i,j).X.y ;
    norm[2].y = -( Cell[i-1][j].Xc.x - nodeNW(i,j).X.x);
    //  normal vector of the SW side of a diamond 
    norm[3].x =  nodeSW(i,j).X.y - Cell[i-1][j].Xc.y;
    norm[3].y = -(nodeSW(i,j).X.x - Cell[i-1][j].Xc.x);
    AREA =  HALF*(fabs((nodeNW(i,j).X-Cell[i][j].Xc)^
		       (nodeSW(i,j).X-Cell[i][j].Xc)) +
		  fabs((nodeNW(i,j).X-Cell[i-1][j].Xc)^
		       (nodeSW(i,j).X-Cell[i-1][j].Xc)));
    if(cell_orientation == CENTER){
      d_dWdx_dW = HALF*((ONE+dWnSWdWc)*norm[0].x+ (ONE+dWnNWdWc)*norm[1].x+ 
			dWnNWdWc* norm[2].x+ dWnSWdWc* norm[3].x)/AREA;
      d_dWdy_dW = HALF*((ONE+dWnSWdWc)*norm[0].y+ (ONE+dWnNWdWc)*norm[1].y+ 
			dWnNWdWc* norm[2].y+ dWnSWdWc* norm[3].y)/AREA;  
    } else if (cell_orientation == WEST) {
      d_dWdx_dW = HALF*(dWnSWdWc*norm[0].x + dWnNWdWc*norm[1].x + 
			(ONE+dWnNWdWc)*norm[2].x+ (ONE+dWnSWdWc)*norm[3].x)/AREA;
      d_dWdy_dW = HALF*(dWnSWdWc*norm[0].y + dWnNWdWc*norm[1].y + 
			(ONE+dWnNWdWc)* norm[2].y+ (ONE+dWnSWdWc)*norm[3].y)/AREA;    
    } else if (cell_orientation == NORTH || cell_orientation == NORTH_WEST) {
      d_dWdx_dW = HALF*(dWnNWdWc*(norm[1].x+norm[2].x))/AREA;
      d_dWdy_dW = HALF*(dWnNWdWc*(norm[1].y+norm[2].y))/AREA;  
    } else if (cell_orientation == SOUTH || cell_orientation == SOUTH_WEST) {
      d_dWdx_dW = HALF*(dWnSWdWc*(norm[0].x+norm[3].x))/AREA;
      d_dWdy_dW = HALF*(dWnSWdWc*(norm[0].y+norm[3].y))/AREA;      
    } /* endif */
    break;
  }
}

/*!
 * Change a block's BC's.
 */
void Grid2D_Quad_Block_HO::set_BCs(const int& FACE, const int& BC){
  switch(FACE){
  case NORTH:
    for ( int i = ICl - Nghost; i <=  ICu + Nghost; i++) {
      BCtypeN[i] = BC;
    } /* endfor */
    break;
  case SOUTH:
    for ( int i = ICl - Nghost; i <=  ICu + Nghost; i++) {
      BCtypeS[i] = BC;
    } /* endfor */
    break;
  case EAST:
    for ( int j = JCl - Nghost; j <=  JCu + Nghost; j++) {
      BCtypeE[j] = BC;
    } /* endfor */
    break;
  case WEST:
    for ( int j = JCl - Nghost; j <=  JCu + Nghost; j++) {
      BCtypeW[j] = BC;
    } /* endfor */
    break;
  default:
    cerr<<"\n WRONG \"FACE\" in set_BCs. \n";
    throw runtime_error("Grid2D_Quad_Block_HO::set_BCs() ERROR! Wrong face!");
  };
}

/*!
 *
 * Create quadrilateral grid block for a Cartesian      
 * mesh defined on a square box with four corner        
 * coordinates (X,Y): (-1,1), (1,1), (1,-1), (-1,1).    
 * This subroutine DOESN'T update the ghost cells or
 * the geometric properties of the grid cells.
 */
void Grid2D_Quad_Block_HO::Create_Quad_Block_Without_Update(const int &Number_of_Cells_Idir,
							    const int &Number_of_Cells_Jdir,
							    const int &Number_of_Ghost_Cells,
							    const int &Highest_Order_of_Reconstruction){
  
  int i, j;
  double S_i, S_j, 
    s_north, s_south, s_east, s_west,
    smax_north, smax_south, smax_east, smax_west, 
    s_i, s_j, smax_i, smax_j,
           w_north, w_south, w_east, w_west, w_total;
  Vector2D x_north, x_south, x_east, x_west;
 
  /* Allocate (re-allocate) memory for the cells and nodes 
     of the quadrilateral mesh block. */
  allocate(Number_of_Cells_Idir, Number_of_Cells_Jdir,
	   Number_of_Ghost_Cells, Highest_Order_of_Reconstruction);
  
  /* Allocate (re-allocate) memory for the boundary splines 
     defining this quadrilateral mesh block. */
  
  if (BndNorthSpline.np != 0) BndNorthSpline.deallocate();
  BndNorthSpline.allocate(2);
  
  if (BndSouthSpline.np != 0) BndSouthSpline.deallocate();
  BndSouthSpline.allocate(2);
  
  if (BndEastSpline.np != 0) BndEastSpline.deallocate();
  BndEastSpline.allocate(2);
  
  if (BndWestSpline.np != 0) BndWestSpline.deallocate();
  BndWestSpline.allocate(2);
  
  /* Set the boundary spline types to linear, assign
     spline points, and calculate the spline pathlengths. */
  
  BndNorthSpline.settype(SPLINE2D_LINEAR);
  BndNorthSpline.Xp[0] = Vector2D(-ONE, ONE);
  BndNorthSpline.Xp[1] = Vector2D( ONE, ONE);
  BndNorthSpline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
  BndNorthSpline.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
  BndNorthSpline.bc[0] = BC_REFLECTION;
  BndNorthSpline.bc[1] = BC_REFLECTION;
  BndNorthSpline.pathlength();
  SminN = BndNorthSpline.sp[0];
  SmaxN = BndNorthSpline.sp[BndNorthSpline.np-1];
  
  BndSouthSpline.settype(SPLINE2D_LINEAR);
  BndSouthSpline.Xp[0] = Vector2D(-ONE,-ONE);
  BndSouthSpline.Xp[1] = Vector2D( ONE,-ONE);
  BndSouthSpline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
  BndSouthSpline.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
  BndSouthSpline.bc[0] = BC_REFLECTION;
  BndSouthSpline.bc[1] = BC_REFLECTION;
  BndSouthSpline.pathlength();
  SminS = BndSouthSpline.sp[0];
  SmaxS = BndSouthSpline.sp[BndSouthSpline.np-1];
  
  BndEastSpline.settype(SPLINE2D_LINEAR);
  BndEastSpline.Xp[0] = Vector2D( ONE,-ONE);
  BndEastSpline.Xp[1] = Vector2D( ONE, ONE);
  BndEastSpline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
  BndEastSpline.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
  BndEastSpline.bc[0] = BC_REFLECTION;
  BndEastSpline.bc[1] = BC_REFLECTION;
  BndEastSpline.pathlength();
  SminE = BndEastSpline.sp[0];
  SmaxE = BndEastSpline.sp[BndEastSpline.np-1];
  
  BndWestSpline.settype(SPLINE2D_LINEAR);
  BndWestSpline.Xp[0] = Vector2D(-ONE,-ONE);
  BndWestSpline.Xp[1] = Vector2D(-ONE, ONE);
  BndWestSpline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
  BndWestSpline.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
  BndWestSpline.bc[0] = BC_REFLECTION;
  BndWestSpline.bc[1] = BC_REFLECTION;
  BndWestSpline.pathlength();
  SminW = BndWestSpline.sp[0];
  SmaxW = BndWestSpline.sp[BndWestSpline.np-1];
  
  /* Compute the interior nodes for the quadrilateral mesh block. */
  
  for ( j = JNl ; j <= JNu ; ++j) {
    S_j = double(j-JNl)/double(JNu-JNl);
    
    smax_east = BndEastSpline.sp[BndEastSpline.np-1]-
      BndEastSpline.sp[0];
    s_east = S_j * smax_east + BndEastSpline.sp[0];
    x_east = Spline(s_east, BndEastSpline);
    
    smax_west = BndWestSpline.sp[BndWestSpline.np-1]-
      BndWestSpline.sp[0];
    s_west = S_j * smax_west + BndWestSpline.sp[0];
    x_west = Spline(s_west, BndWestSpline);
    
    for ( i = INl ; i <= INu ; ++i) {
      S_i = double(i-INl)/double(INu-INl);
      
      smax_north = BndNorthSpline.sp[BndNorthSpline.np-1]-
	BndNorthSpline.sp[0];
      s_north = S_i * smax_north + BndNorthSpline.sp[0];
      x_north = Spline(s_north, BndNorthSpline);
      
      smax_south = BndSouthSpline.sp[BndSouthSpline.np-1]-
	BndSouthSpline.sp[0];
      s_south = S_i * smax_south + BndSouthSpline.sp[0];
      x_south = Spline(s_south, BndSouthSpline);
      
      if (j == JNl) {
	Node[i][j].X = x_south;
      } else if (j == JNu) {
	Node[i][j].X = x_north;
      } else if (i == INl) {
	Node[i][j].X = x_west;
      } else if (i == INu) {
	Node[i][j].X = x_east;
      } else {
	s_i = (ONE-S_j)*s_south + S_j*s_north;
	s_j = (ONE-S_i)*s_west + S_i*s_east;
	
	smax_i = (ONE-S_j)*smax_south + S_j*smax_north;
	smax_j = (ONE-S_i)*smax_west + S_i*smax_east;
        
	w_south = (ONE-s_j/smax_j);
	w_north = (s_j/smax_j);
	w_total = w_south + w_north;
	Node[i][j].X = (w_south*x_south + 
			     w_north*x_north)/w_total;
      } /* endif */
      
    } /* endfor */
  }/* endfor */
  
  /* Require update of the interior cells geometric properties. */
  Schedule_Interior_Mesh_Update();
  
  /* Set the boundary condition types at the grid block
     boundaries. */
  
  Set_BCs();
    
}

/*!
 * Create a 2D quadrilateral grid block with the four   
 * boundaries of the block defined by blended splines   
 * which are given in the input boundary spline file:   
 *                                                      
 * Bnd_Spline_File_Name_ptr. 
 *                           
 * This subroutine DOESN'T update the ghost cells or
 * the geometric properties of the grid cells.                           
 */
void Grid2D_Quad_Block_HO::Create_Quad_Block_Without_Update(char *Bnd_Spline_File_Name_ptr,
							    const int &Number_of_Cells_Idir,
							    const int &Number_of_Cells_Jdir,
							    const int &Number_of_Ghost_Cells,
							    const int &Highest_Order_of_Reconstruction) {
  
  ifstream bnd_spline_file;
  int i, j, k, kx, ky, kx_max, ky_max, npts, spline_type;
  int node_init_procedure;
  double S_i, S_j, 
    s_north, s_south, s_east, s_west,
    smax_north, smax_south, smax_east, smax_west, 
    s_i, s_j, smax_i, smax_j,
    w_north, w_south, w_east, w_west, w_total;
  double smax_S1_NS, smax_S2_NS, smax_S1_EW, smax_S2_EW,
    S1_i, S2_i, S1_j, S2_j, dS_i, dS_j;
  Vector2D x_north, x_south, x_east, x_west;
  Spline2D_HO S1_NS, S2_NS, S1_EW, S2_EW;
  
  /* Allocate (re-allocate) memory for the cells and nodes 
     of the quadrilateral mesh block. */

  allocate(Number_of_Cells_Idir, Number_of_Cells_Jdir,
	   Number_of_Ghost_Cells, Highest_Order_of_Reconstruction);

  /* Open the data file containing the boundary splines. */

  bnd_spline_file.open(Bnd_Spline_File_Name_ptr,ios::in);

  /* For each of the north, south, east, and west boundaries
     of this mesh block, read in the number of spline points, 
     allocate memory for the boundary splines, read in the 
     spline points, and finally calculate the spline 
     pathlengths. */

  bnd_spline_file.setf(ios::skipws);
  bnd_spline_file >> npts;
  bnd_spline_file.unsetf(ios::skipws);
  if (BndNorthSpline.np != 0) BndNorthSpline.deallocate();
  BndNorthSpline.allocate(npts);
  bnd_spline_file.setf(ios::skipws);
  bnd_spline_file >> spline_type;
  bnd_spline_file.unsetf(ios::skipws);
  BndNorthSpline.settype(spline_type);
  bnd_spline_file >> BndNorthSpline;
  BndNorthSpline.pathlength();
  SminN = BndNorthSpline.sp[0];
  SmaxN = BndNorthSpline.sp[BndNorthSpline.np-1];

  bnd_spline_file.setf(ios::skipws);
  bnd_spline_file >> npts;
  bnd_spline_file.unsetf(ios::skipws);
  if (BndSouthSpline.np != 0) BndSouthSpline.deallocate();
  BndSouthSpline.allocate(npts);
  bnd_spline_file.setf(ios::skipws);
  bnd_spline_file >> spline_type;
  bnd_spline_file.unsetf(ios::skipws);
  BndSouthSpline.settype(spline_type);
  bnd_spline_file >> BndSouthSpline;
  BndSouthSpline.pathlength();
  SminS = BndSouthSpline.sp[0];
  SmaxS = BndSouthSpline.sp[BndSouthSpline.np-1];

  bnd_spline_file.setf(ios::skipws);
  bnd_spline_file >> npts;
  bnd_spline_file.unsetf(ios::skipws);
  if (BndEastSpline.np != 0) BndEastSpline.deallocate();
  BndEastSpline.allocate(npts);
  bnd_spline_file.setf(ios::skipws);
  bnd_spline_file >> spline_type;
  bnd_spline_file.unsetf(ios::skipws);
  BndEastSpline.settype(spline_type);
  bnd_spline_file >> BndEastSpline;
  BndEastSpline.pathlength();
  SminE = BndEastSpline.sp[0];
  SmaxE = BndEastSpline.sp[BndEastSpline.np-1];

  bnd_spline_file.setf(ios::skipws);
  bnd_spline_file >> npts;
  bnd_spline_file.unsetf(ios::skipws);
  if (BndWestSpline.np != 0) BndWestSpline.deallocate();
  BndWestSpline.allocate(npts);
  bnd_spline_file.setf(ios::skipws);
  bnd_spline_file >> spline_type;
  bnd_spline_file.unsetf(ios::skipws);
  BndWestSpline.settype(spline_type);
  bnd_spline_file >> BndWestSpline;
  BndWestSpline.pathlength();
  SminW = BndWestSpline.sp[0];
  SmaxW = BndWestSpline.sp[BndWestSpline.np-1];

  /* Read the node initialization procedure for this 
     quadrilateral grid block. */

  bnd_spline_file.setf(ios::skipws);
  bnd_spline_file >> node_init_procedure;
  bnd_spline_file.unsetf(ios::skipws);

  /* Input the type of node distribution stretching 
     functions to be used for each of the coordinate
     directions. Also read in the related stretching
     parameters. */

  bnd_spline_file.setf(ios::skipws);
  bnd_spline_file >> StretchI >> BetaI >> TauI;
  bnd_spline_file >> StretchJ >> BetaJ >> TauJ;
  bnd_spline_file.unsetf(ios::skipws);

  /* Input the grid orthogonality specifiers for
     the north, south, east, and west boundaries. */

  bnd_spline_file.setf(ios::skipws);
  bnd_spline_file >> OrthogonalN >> OrthogonalS
		  >> OrthogonalE >> OrthogonalW;
  bnd_spline_file.unsetf(ios::skipws);

  /* Close the data file containing the boundary splines. */

  bnd_spline_file.close();

  /* Compute the interior nodes for the quadrilateral mesh block. */

  smax_north = BndNorthSpline.sp[BndNorthSpline.np-1]-
    BndNorthSpline.sp[0];
  smax_south = BndSouthSpline.sp[BndSouthSpline.np-1]-
    BndSouthSpline.sp[0];

  smax_east = BndEastSpline.sp[BndEastSpline.np-1]-
    BndEastSpline.sp[0];
  smax_west = BndWestSpline.sp[BndWestSpline.np-1]-
    BndWestSpline.sp[0];

  dS_i = ONE/(TEN*double(INu-INl));
  dS_j = ONE/(TEN*double(JNu-JNl));

  for ( j = JNl ; j <= JNu ; ++j) {
    for ( i = INl ; i <= INu ; ++i) {
      S_j = double(j-JNl)/double(JNu-JNl);
      S_j = StretchingFcn(S_j, BetaJ, TauJ, StretchJ);
      s_east = S_j * smax_east + BndEastSpline.sp[0];
      x_east = Spline(s_east, BndEastSpline);
      s_west = S_j * smax_west + BndWestSpline.sp[0];
      x_west = Spline(s_west, BndWestSpline);

      S_i = double(i-INl)/double(INu-INl);
      S_i = StretchingFcn(S_i, BetaI, TauI, StretchI);
      s_north = S_i * smax_north + BndNorthSpline.sp[0];
      x_north = Spline(s_north, BndNorthSpline);
      s_south = S_i * smax_south + BndSouthSpline.sp[0];
      x_south = Spline(s_south, BndSouthSpline);

      if (j == JNl) {
	Node[i][j].X = x_south;
      } else if (j == JNu) {
	Node[i][j].X = x_north;
      } else if (i == INl) {
	Node[i][j].X = x_west;
      } else if (i == INu) {
	Node[i][j].X = x_east;
      } else {
	s_i = (ONE-S_j)*s_south + S_j*s_north;
	s_j = (ONE-S_i)*s_west + S_i*s_east;

	smax_i = (ONE-S_j)*smax_south + S_j*smax_north;
	smax_j = (ONE-S_i)*smax_west + S_i*smax_east;

	switch(node_init_procedure) {
	  //============================================================
	case GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH :
	  w_south = (ONE-s_j/smax_j);
	  w_north = (s_j/smax_j);
	  w_total = w_south + w_north;
	  Node[i][j].X = (w_south*x_south + 
			       w_north*x_north)/w_total;
	  break;
	  //============================================================
	case GRID2D_QUAD_BLOCK_INIT_PROCEDURE_EAST_WEST :
	  w_west =  (ONE-s_i/smax_i);
	  w_east =  (s_i/smax_i);
	  w_total = w_east + w_west;
	  Node[i][j].X = (w_west*x_west + 
			       w_east*x_east)/w_total;
	  break;
	  //============================================================
	case GRID2D_QUAD_BLOCK_INIT_PROCEDURE_TRANS_FINITE_XY :
	  kx_max = 5;
	  dS_i = ONE/double(kx_max-1);
	  kx = int(floor(S_i/dS_i));
	  S1_i = max(ZERO, double(kx)*dS_i);
	  S2_i = min(ONE, S1_i + dS_i);

	  ky_max = 7;
	  dS_j = ONE/double(ky_max-1);
	  ky = int(floor(S_j/dS_j));
	  S1_j = max(ZERO, double(ky)*dS_j);
	  S2_j = min(ONE, S1_j + dS_j);

	  if (S1_i - TOLER > ZERO) {
	    S1_NS.allocate(2);
	    S1_NS.settype(SPLINE2D_LINEAR);

	    s_south = S1_i * smax_south + BndSouthSpline.sp[0];
	    S1_NS.Xp[0] = Spline(s_south, BndSouthSpline);
	    S1_NS.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	    s_north = S1_i * smax_north + BndNorthSpline.sp[0],
	      S1_NS.Xp[1] = Spline(s_north, BndNorthSpline);
	    S1_NS.tp[1] = SPLINE2D_POINT_SHARP_CORNER;

	    S1_NS.pathlength();
	    smax_S1_NS = S1_NS.sp[S1_NS.np-1]-
	      S1_NS.sp[0];
	  } else {
	    S1_NS = BndWestSpline;
	    smax_S1_NS = S1_NS.sp[S1_NS.np-1]-
	      S1_NS.sp[0];
	  } /* endif */

	  if (S2_i + TOLER < ONE) {
	    S2_NS.allocate(2);
	    S2_NS.settype(SPLINE2D_LINEAR);

	    s_south = S2_i * smax_south + BndSouthSpline.sp[0];
	    S2_NS.Xp[0] = Spline(s_south, BndSouthSpline);
	    S2_NS.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	    s_north = S2_i * smax_north + BndNorthSpline.sp[0],
	      S2_NS.Xp[1] = Spline(s_north, BndNorthSpline);
	    S2_NS.tp[1] = SPLINE2D_POINT_SHARP_CORNER;

	    S2_NS.pathlength();
	    smax_S2_NS = S2_NS.sp[S2_NS.np-1]-
	      S2_NS.sp[0];
	  } else {
	    S2_NS = BndEastSpline;
	    smax_S2_NS = S2_NS.sp[S2_NS.np-1]-
	      S2_NS.sp[0];
	  } /* endif */

	  if (S1_j - TOLER > ZERO) {
	    S1_EW.allocate(kx_max);
	    S1_EW.settype(SPLINE2D_LINEAR);

	    s_west = S1_j * smax_west + BndWestSpline.sp[0];
	    S1_EW.Xp[0] = Spline(s_west, BndWestSpline);
	    S1_EW.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	    s_east = S1_j * smax_east + BndEastSpline.sp[0];
	    S1_EW.Xp[kx_max-1] = Spline(s_east, BndEastSpline);
	    S1_EW.tp[kx_max-1] = SPLINE2D_POINT_SHARP_CORNER;

	    for ( k = 1 ; k <= kx_max-2 ; ++k) {                     
	      s_north = double(k)*smax_north/double(kx_max-1) + BndNorthSpline.sp[0];
	      x_north = Spline(s_north, BndNorthSpline);

	      s_south = double(k)*smax_south/double(kx_max-1) + BndSouthSpline.sp[0];
	      x_south = Spline(s_south, BndSouthSpline);

	      w_south =  ONE-S1_j;
	      w_north =  S1_j;
	      w_total = w_north + w_south;
	      S1_EW.Xp[k] = (w_south*x_south + w_north*x_north)/w_total;
	      S1_EW.tp[k] = SPLINE2D_POINT_SHARP_CORNER;
	    } /* endfor */

	    S1_EW.pathlength();
	    smax_S1_EW = S1_EW.sp[S1_EW.np-1]-
	      S1_EW.sp[0];
	  } else {
	    S1_EW = BndSouthSpline;
	    smax_S1_EW = S1_EW.sp[S1_EW.np-1]-
	      S1_EW.sp[0];
	  } /* endif */

	  if (S2_j + TOLER < ONE) {
	    S2_EW.allocate(kx_max);
	    S2_EW.settype(SPLINE2D_LINEAR);

	    s_west = S2_j * smax_west + BndWestSpline.sp[0];
	    S2_EW.Xp[0] = Spline(s_west, BndWestSpline);
	    S2_EW.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	    s_east = S2_j * smax_east + BndEastSpline.sp[0];
	    S2_EW.Xp[kx_max-1] = Spline(s_east, BndEastSpline);
	    S2_EW.tp[kx_max-1] = SPLINE2D_POINT_SHARP_CORNER;

	    for ( k = 1 ; k <= kx_max-2 ; ++k) {                     
	      s_north = double(k)*smax_north/double(kx_max-1) + BndNorthSpline.sp[0];
	      x_north = Spline(s_north, BndNorthSpline);

	      s_south = double(k)*smax_south/double(kx_max-1) + BndSouthSpline.sp[0];
	      x_south = Spline(s_south, BndSouthSpline);

	      w_south =  ONE-S2_j;
	      w_north =  S2_j;
	      w_total = w_north + w_south;
	      S2_EW.Xp[k] = (w_south*x_south + w_north*x_north)/w_total;
	      S2_EW.tp[k] = SPLINE2D_POINT_SHARP_CORNER;
	    } /* endfor */

	    S2_EW.pathlength();
	    smax_S2_EW = S2_EW.sp[S2_EW.np-1]-
	      S2_EW.sp[0];
	  } else {
	    S2_EW = BndNorthSpline;
	    smax_S2_EW = S2_EW.sp[S2_EW.np-1]-
	      S2_EW.sp[0];
	  } /* endif */

	  s_north = S_i * smax_S2_EW + S2_EW.sp[0];
	  x_north = Spline(s_north, S2_EW);
	  s_south = S_i * smax_S1_EW + S1_EW.sp[0];
	  x_south = Spline(s_south, S1_EW);

	  s_east = S_j * smax_S2_NS + S2_NS.sp[0];
	  x_east = Spline(s_east, S2_NS);
	  s_west = S_j * smax_S1_NS + S1_NS.sp[0];
	  x_west = Spline(s_west, S1_NS);

	  if (ky == 0 || ky == ky_max - 2) {
	    w_south =  ONE-(S_j-S1_j)/(S2_j-S1_j);
	    w_north =  (S_j-S1_j)/(S2_j-S1_j);
	    w_total = w_south + w_north;
	    Node[i][j].X = (w_south*x_south + 
				 w_north*x_north)/w_total;
	  } else {
	    w_west =  ONE-(S_i-S1_i)/(S2_i-S1_i);
	    w_east =  (S_i-S1_i)/(S2_i-S1_i);
	    w_total = w_east + w_west;
	    Node[i][j].X = (w_west*x_west + 
				 w_east*x_east)/w_total;
	  } /* endif */

	  S1_NS.deallocate(); S2_NS.deallocate();
	  S2_EW.deallocate(); S2_EW.deallocate();
	  break;
	  //============================================================
	case GRID2D_QUAD_BLOCK_INIT_PROCEDURE_TRANS_FINITE_YX :
	  kx_max = 7;
	  dS_i = ONE/double(kx_max-1);
	  kx = int(floor(S_i/dS_i));
	  S1_i = max(ZERO, double(kx)*dS_i);
	  S2_i = min(ONE, S1_i + dS_i);

	  ky_max = 5;
	  dS_j = ONE/double(ky_max-1);
	  ky = int(floor(S_j/dS_j));
	  S1_j = max(ZERO, double(ky)*dS_j);
	  S2_j = min(ONE, S1_j + dS_j);

	  if (S1_i - TOLER > ZERO) {
	    S1_NS.allocate(ky_max);
	    S1_NS.settype(SPLINE2D_LINEAR);

	    s_south = S1_i * smax_south + BndSouthSpline.sp[0];
	    S1_NS.Xp[0] = Spline(s_south, BndSouthSpline);
	    S1_NS.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	    s_north = S1_i * smax_north + BndNorthSpline.sp[0],
	      S1_NS.Xp[ky_max-1] = Spline(s_north, BndNorthSpline);
	    S1_NS.tp[ky_max-1] = SPLINE2D_POINT_SHARP_CORNER;

	    for ( k = 1 ; k <= ky_max-2 ; ++k) {                     
	      s_east = double(k)*smax_east/double(ky_max-1) + BndEastSpline.sp[0];
	      x_east = Spline(s_east, BndEastSpline);

	      s_west = double(k)*smax_west/double(ky_max-1) + BndWestSpline.sp[0];
	      x_west = Spline(s_west, BndWestSpline);

	      w_west =  ONE-S1_i;
	      w_east =  S1_i;
	      w_total = w_east + w_west;
	      S1_NS.Xp[k] = (w_west*x_west + w_east*x_east)/w_total;
	      S1_NS.tp[k] = SPLINE2D_POINT_SHARP_CORNER;
	    } /* endfor */

	    S1_NS.pathlength();
	    smax_S1_NS = S1_NS.sp[S1_NS.np-1]-
	      S1_NS.sp[0];
	  } else {
	    S1_NS = BndWestSpline;
	    smax_S1_NS = S1_NS.sp[S1_NS.np-1]-
	      S1_NS.sp[0];
	  } /* endif */

	  if (S2_i + TOLER < ONE) {
	    S2_NS.allocate(ky_max);
	    S2_NS.settype(SPLINE2D_LINEAR);

	    s_south = S2_i * smax_south + BndSouthSpline.sp[0];
	    S2_NS.Xp[0] = Spline(s_south, BndSouthSpline);
	    S2_NS.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	    s_north = S2_i * smax_north + BndNorthSpline.sp[0],
	      S2_NS.Xp[ky_max-1] = Spline(s_north, BndNorthSpline);
	    S2_NS.tp[ky_max-1] = SPLINE2D_POINT_SHARP_CORNER;

	    for ( k = 1 ; k <= ky_max-2 ; ++k) {
	      s_east = double(k)*smax_east/double(ky_max-1) + BndEastSpline.sp[0];
	      x_east = Spline(s_east, BndEastSpline);

	      s_west = double(k)*smax_west/double(ky_max-1) + BndWestSpline.sp[0];
	      x_west = Spline(s_west, BndWestSpline);

	      w_west =  ONE-S2_i;
	      w_east =  S2_i;
	      w_total = w_east + w_west;
	      S2_NS.Xp[k] = (w_west*x_west + w_east*x_east)/w_total;
	      S2_NS.tp[k] = SPLINE2D_POINT_SHARP_CORNER;
	    } /* endfor */

	    S2_NS.pathlength();
	    smax_S2_NS = S2_NS.sp[S2_NS.np-1]-
	      S2_NS.sp[0];
	  } else {
	    S2_NS = BndEastSpline;
	    smax_S2_NS = S2_NS.sp[S2_NS.np-1]-
	      S2_NS.sp[0];
	  } /* endif */

	  if (S1_j - TOLER > ZERO) {
	    S1_EW.allocate(2);
	    S1_EW.settype(SPLINE2D_LINEAR);

	    s_west = S1_j * smax_west + BndWestSpline.sp[0];
	    S1_EW.Xp[0] = Spline(s_west, BndWestSpline);
	    S1_EW.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	    s_east = S1_j * smax_east + BndEastSpline.sp[0];
	    S1_EW.Xp[1] = Spline(s_east, BndEastSpline);
	    S1_EW.tp[1] = SPLINE2D_POINT_SHARP_CORNER;

	    S1_EW.pathlength();
	    smax_S1_EW = S1_EW.sp[S1_EW.np-1]-
	      S1_EW.sp[0];
	  } else {
	    S1_EW = BndSouthSpline;
	    smax_S1_EW = S1_EW.sp[S1_EW.np-1]-
	      S1_EW.sp[0];
	  } /* endif */

	  if (S2_j + TOLER < ONE) {
	    S2_EW.allocate(2);
	    S2_EW.settype(SPLINE2D_LINEAR);

	    s_west = S2_j * smax_west + BndWestSpline.sp[0];
	    S2_EW.Xp[0] = Spline(s_west, BndWestSpline);
	    S2_EW.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	    s_east = S2_j * smax_east + BndEastSpline.sp[0];
	    S2_EW.Xp[1] = Spline(s_east, BndEastSpline);
	    S2_EW.tp[1] = SPLINE2D_POINT_SHARP_CORNER;

	    S2_EW.pathlength();
	    smax_S2_EW = S2_EW.sp[S2_EW.np-1]-
	      S2_EW.sp[0];
	  } else {
	    S2_EW = BndNorthSpline;
	    smax_S2_EW = S2_EW.sp[S2_EW.np-1]-
	      S2_EW.sp[0];
	  } /* endif */

	  s_north = S_i * smax_S2_EW + S2_EW.sp[0];
	  x_north = Spline(s_north, S2_EW);
	  s_south = S_i * smax_S1_EW + S1_EW.sp[0];
	  x_south = Spline(s_south, S1_EW);

	  s_east = S_j * smax_S2_NS + S2_NS.sp[0];
	  x_east = Spline(s_east, S2_NS);
	  s_west = S_j * smax_S1_NS + S1_NS.sp[0];
	  x_west = Spline(s_west, S1_NS);

	  if (kx == 0 || kx == kx_max - 2) {
	    w_west =  ONE-(S_i-S1_i)/(S2_i-S1_i);
	    w_east =  (S_i-S1_i)/(S2_i-S1_i);
	    w_total = w_east + w_west;
	    Node[i][j].X = (w_west*x_west + 
				 w_east*x_east)/w_total;
	  } else {
	    w_south =  ONE-(S_j-S1_j)/(S2_j-S1_j);
	    w_north =  (S_j-S1_j)/(S2_j-S1_j);
	    w_total = w_south + w_north;
	    Node[i][j].X = (w_south*x_south + 
				 w_north*x_north)/w_total;
	  } /* endif */

	  S1_NS.deallocate(); S2_NS.deallocate();
	  S2_EW.deallocate(); S2_EW.deallocate();
	  //============================================================
	default:
	  w_south = (ONE-s_j/smax_j);
	  w_north = (s_j/smax_j);
	  w_total = w_south + w_north;
	  Node[i][j].X = (w_south*x_south + 
			       w_north*x_north)/w_total;
	  break;
	} /* endswitch */

      } /* endif */

    } /* endfor */
  }/* endfor */

  /* Require update of the interior cells geometric properties. */
  Schedule_Interior_Mesh_Update();

  /* Set the boundary condition types at the quadrilateral 
     grid block boundaries. */

  Set_BCs();

}


/*!
 * Create a 2D quadrilateral grid block with the four   
 * boundaries of the block defined by blended splines   
 * which are given as input arguments to the routine.   
 *                                                      
 * This subroutine DOESN'T update the ghost cells or
 * the geometric properties of the grid cells.
 */
void Grid2D_Quad_Block_HO::Create_Quad_Block_Without_Update(Spline2D_HO &Bnd_Spline_North,
							    Spline2D_HO &Bnd_Spline_South,
							    Spline2D_HO &Bnd_Spline_East,
							    Spline2D_HO &Bnd_Spline_West,
							    const int &Number_of_Cells_Idir,
							    const int &Number_of_Cells_Jdir,
							    const int &Number_of_Ghost_Cells,
							    const int &Highest_Order_of_Reconstruction,
							    const int &Node_Init_Procedure,
							    const int &Stretch_I,
							    const double &Beta_I, 
							    const double &Tau_I,
							    const int & Stretch_J,
							    const double &Beta_J,
							    const double &Tau_J,
							    const int &Orthogonal_North,
							    const int &Orthogonal_South,
							    const int &Orthogonal_East,
							    const int &Orthogonal_West) {

  int i, j, k, kx, ky, kx_max, ky_max, npts, spline_type;
  double S_i, S_j, 
    s_north, s_south, s_east, s_west,
    smax_north, smax_south, smax_east, smax_west, 
    s_i, s_j, smax_i, smax_j,
    w_north, w_south, w_east, w_west, w_total;
  double smax_S1_NS, smax_S2_NS, smax_S1_EW, smax_S2_EW,
    S1_i, S2_i, S1_j, S2_j, dS_i, dS_j;
  Vector2D x_north, x_south, x_east, x_west;
  Spline2D_HO S1_NS, S2_NS, S1_EW, S2_EW;

  /* Allocate (re-allocate) memory for the cells and nodes 
     of the quadrilateral mesh block. */

  allocate(Number_of_Cells_Idir, Number_of_Cells_Jdir,
	   Number_of_Ghost_Cells, Highest_Order_of_Reconstruction);

  /* For each of the north, south, east, and west boundaries
     of this mesh block, assign the boundary splines specified
     by the subroutine input parameters. */

  BndNorthSpline = Bnd_Spline_North;
  BndSouthSpline = Bnd_Spline_South;
  BndEastSpline  = Bnd_Spline_East;
  BndWestSpline  = Bnd_Spline_West;

  SminN = BndNorthSpline.sp[0];
  SmaxN = BndNorthSpline.sp[BndNorthSpline.np-1];
  SminS = BndSouthSpline.sp[0];
  SmaxS = BndSouthSpline.sp[BndSouthSpline.np-1];
  SminE = BndEastSpline.sp[0];
  SmaxE = BndEastSpline.sp[BndEastSpline.np-1];
  SminW = BndWestSpline.sp[0];
  SmaxW = BndWestSpline.sp[BndWestSpline.np-1];

  /* Assign values for the type of node distribution 
     stretching functions to be used for each of the 
     coordinate directions. Also assign values to the
     related stretching parameters. */

  StretchI = Stretch_I;
  BetaI = Beta_I;
  TauI = Tau_I;
  StretchJ = Stretch_J;
  BetaJ = Beta_J;
  TauJ = Tau_J;

  /* Assign values to the grid orthogonality specifiers for
     the north, south, east, and west boundaries. */

  OrthogonalN = Orthogonal_North;
  OrthogonalS = Orthogonal_South;
  OrthogonalE = Orthogonal_East;
  OrthogonalW = Orthogonal_West;

  /* Compute the interior nodes for the quadrilateral mesh block. */

  smax_north = BndNorthSpline.sp[BndNorthSpline.np-1]-
    BndNorthSpline.sp[0];
  smax_south = BndSouthSpline.sp[BndSouthSpline.np-1]-
    BndSouthSpline.sp[0];

  smax_east = BndEastSpline.sp[BndEastSpline.np-1]-
    BndEastSpline.sp[0];
  smax_west = BndWestSpline.sp[BndWestSpline.np-1]-
    BndWestSpline.sp[0];

  for ( j = JNl ; j <= JNu ; ++j) {
    for ( i = INl ; i <= INu ; ++i) {
      S_j = double(j-JNl)/double(JNu-JNl);
      S_j = StretchingFcn(S_j, BetaJ, TauJ, StretchJ);
      s_east = S_j * smax_east + BndEastSpline.sp[0];
      x_east = Spline(s_east, BndEastSpline);
      s_west = S_j * smax_west + BndWestSpline.sp[0];
      x_west = Spline(s_west, BndWestSpline);

      S_i = double(i-INl)/double(INu-INl);
      S_i = StretchingFcn(S_i, BetaI, TauI, StretchI);
      s_north = S_i * smax_north + BndNorthSpline.sp[0];
      x_north = Spline(s_north, BndNorthSpline);
      s_south = S_i * smax_south + BndSouthSpline.sp[0];
      x_south = Spline(s_south, BndSouthSpline);

      if (j == JNl) {
	Node[i][j].X = x_south;
      } else if (j == JNu) {
	Node[i][j].X = x_north;
      } else if (i == INl) {
	Node[i][j].X = x_west;
      } else if (i == INu) {
	Node[i][j].X = x_east;
      } else {
	s_i = (ONE-S_j)*s_south + S_j*s_north;
	s_j = (ONE-S_i)*s_west + S_i*s_east;

	smax_i = (ONE-S_j)*smax_south + S_j*smax_north;
	smax_j = (ONE-S_i)*smax_west + S_i*smax_east;

	switch(Node_Init_Procedure) {
	  //============================================================
	case GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH :
	  w_south = (ONE-s_j/smax_j);
	  w_north = (s_j/smax_j);
	  w_total = w_south + w_north;
	  Node[i][j].X = (w_south*x_south + 
			  w_north*x_north)/w_total;
	  break;
	  //============================================================
	case GRID2D_QUAD_BLOCK_INIT_PROCEDURE_EAST_WEST :
	  w_west =  (ONE-s_i/smax_i);
	  w_east =  (s_i/smax_i);
	  w_total = w_east + w_west;
	  Node[i][j].X = (w_west*x_west + 
			  w_east*x_east)/w_total;
	  break;
	  //============================================================
	case GRID2D_QUAD_BLOCK_INIT_PROCEDURE_TRANS_FINITE_XY :
	  kx_max = 5;
	  dS_i = ONE/double(kx_max-1);
	  kx = int(floor(S_i/dS_i));
	  S1_i = max(ZERO, double(kx)*dS_i);
	  S2_i = min(ONE, S1_i + dS_i);

	  ky_max = 7;
	  dS_j = ONE/double(ky_max-1);
	  ky = int(floor(S_j/dS_j));
	  S1_j = max(ZERO, double(ky)*dS_j);
	  S2_j = min(ONE, S1_j + dS_j);

	  if (S1_i - TOLER > ZERO) {
	    S1_NS.allocate(2);
	    S1_NS.settype(SPLINE2D_LINEAR);

	    s_south = S1_i * smax_south + BndSouthSpline.sp[0];
	    S1_NS.Xp[0] = Spline(s_south, BndSouthSpline);
	    S1_NS.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	    s_north = S1_i * smax_north + BndNorthSpline.sp[0],
	      S1_NS.Xp[1] = Spline(s_north, BndNorthSpline);
	    S1_NS.tp[1] = SPLINE2D_POINT_SHARP_CORNER;

	    S1_NS.pathlength();
	    smax_S1_NS = S1_NS.sp[S1_NS.np-1]-
	      S1_NS.sp[0];
	  } else {
	    S1_NS = BndWestSpline;
	    smax_S1_NS = S1_NS.sp[S1_NS.np-1]-
	      S1_NS.sp[0];
	  } /* endif */

	  if (S2_i + TOLER < ONE) {
	    S2_NS.allocate(2);
	    S2_NS.settype(SPLINE2D_LINEAR);

	    s_south = S2_i * smax_south + BndSouthSpline.sp[0];
	    S2_NS.Xp[0] = Spline(s_south, BndSouthSpline);
	    S2_NS.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	    s_north = S2_i * smax_north + BndNorthSpline.sp[0],
	      S2_NS.Xp[1] = Spline(s_north, BndNorthSpline);
	    S2_NS.tp[1] = SPLINE2D_POINT_SHARP_CORNER;

	    S2_NS.pathlength();
	    smax_S2_NS = S2_NS.sp[S2_NS.np-1]-
	      S2_NS.sp[0];
	  } else {
	    S2_NS = BndEastSpline;
	    smax_S2_NS = S2_NS.sp[S2_NS.np-1]-
	      S2_NS.sp[0];
	  } /* endif */

	  if (S1_j - TOLER > ZERO) {
	    S1_EW.allocate(kx_max);
	    S1_EW.settype(SPLINE2D_LINEAR);

	    s_west = S1_j * smax_west + BndWestSpline.sp[0];
	    S1_EW.Xp[0] = Spline(s_west, BndWestSpline);
	    S1_EW.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	    s_east = S1_j * smax_east + BndEastSpline.sp[0];
	    S1_EW.Xp[kx_max-1] = Spline(s_east, BndEastSpline);
	    S1_EW.tp[kx_max-1] = SPLINE2D_POINT_SHARP_CORNER;

	    for ( k = 1 ; k <= kx_max-2 ; ++k) {                     
	      s_north = double(k)*smax_north/double(kx_max-1) + BndNorthSpline.sp[0];
	      x_north = Spline(s_north, BndNorthSpline);

	      s_south = double(k)*smax_south/double(kx_max-1) + BndSouthSpline.sp[0];
	      x_south = Spline(s_south, BndSouthSpline);

	      w_south =  ONE-S1_j;
	      w_north =  S1_j;
	      w_total = w_north + w_south;
	      S1_EW.Xp[k] = (w_south*x_south + w_north*x_north)/w_total;
	      S1_EW.tp[k] = SPLINE2D_POINT_SHARP_CORNER;
	    } /* endfor */

	    S1_EW.pathlength();
	    smax_S1_EW = S1_EW.sp[S1_EW.np-1]-
	      S1_EW.sp[0];
	  } else {
	    S1_EW = BndSouthSpline;
	    smax_S1_EW = S1_EW.sp[S1_EW.np-1]-
	      S1_EW.sp[0];
	  } /* endif */

	  if (S2_j + TOLER < ONE) {
	    S2_EW.allocate(kx_max);
	    S2_EW.settype(SPLINE2D_LINEAR);

	    s_west = S2_j * smax_west + BndWestSpline.sp[0];
	    S2_EW.Xp[0] = Spline(s_west, BndWestSpline);
	    S2_EW.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	    s_east = S2_j * smax_east + BndEastSpline.sp[0];
	    S2_EW.Xp[kx_max-1] = Spline(s_east, BndEastSpline);
	    S2_EW.tp[kx_max-1] = SPLINE2D_POINT_SHARP_CORNER;

	    for ( k = 1 ; k <= kx_max-2 ; ++k) {                     
	      s_north = double(k)*smax_north/double(kx_max-1) + BndNorthSpline.sp[0];
	      x_north = Spline(s_north, BndNorthSpline);

	      s_south = double(k)*smax_south/double(kx_max-1) + BndSouthSpline.sp[0];
	      x_south = Spline(s_south, BndSouthSpline);

	      w_south =  ONE-S2_j;
	      w_north =  S2_j;
	      w_total = w_north + w_south;
	      S2_EW.Xp[k] = (w_south*x_south + w_north*x_north)/w_total;
	      S2_EW.tp[k] = SPLINE2D_POINT_SHARP_CORNER;
	    } /* endfor */

	    S2_EW.pathlength();
	    smax_S2_EW = S2_EW.sp[S2_EW.np-1]-
	      S2_EW.sp[0];
	  } else {
	    S2_EW = BndNorthSpline;
	    smax_S2_EW = S2_EW.sp[S2_EW.np-1]-
	      S2_EW.sp[0];
	  } /* endif */

	  s_north = S_i * smax_S2_EW + S2_EW.sp[0];
	  x_north = Spline(s_north, S2_EW);
	  s_south = S_i * smax_S1_EW + S1_EW.sp[0];
	  x_south = Spline(s_south, S1_EW);

	  s_east = S_j * smax_S2_NS + S2_NS.sp[0];
	  x_east = Spline(s_east, S2_NS);
	  s_west = S_j * smax_S1_NS + S1_NS.sp[0];
	  x_west = Spline(s_west, S1_NS);
      
	  if (ky == 0 || ky == ky_max - 2) {
	    w_south =  ONE-(S_j-S1_j)/(S2_j-S1_j);
	    w_north =  (S_j-S1_j)/(S2_j-S1_j);
	    w_total = w_south + w_north;
	    Node[i][j].X = (w_south*x_south + 
			    w_north*x_north)/w_total;
	  } else {
	    w_west =  ONE-(S_i-S1_i)/(S2_i-S1_i);
	    w_east =  (S_i-S1_i)/(S2_i-S1_i);
	    w_total = w_east + w_west;
	    Node[i][j].X = (w_west*x_west + 
			    w_east*x_east)/w_total;
	  } /* endif */

	  S1_NS.deallocate(); S2_NS.deallocate();
	  S2_EW.deallocate(); S2_EW.deallocate();
	  break;
	  //============================================================
	case GRID2D_QUAD_BLOCK_INIT_PROCEDURE_TRANS_FINITE_YX :
	  kx_max = 7;
	  dS_i = ONE/double(kx_max-1);
	  kx = int(floor(S_i/dS_i));
	  S1_i = max(ZERO, double(kx)*dS_i);
	  S2_i = min(ONE, S1_i + dS_i);

	  ky_max = 5;
	  dS_j = ONE/double(ky_max-1);
	  ky = int(floor(S_j/dS_j));
	  S1_j = max(ZERO, double(ky)*dS_j);
	  S2_j = min(ONE, S1_j + dS_j);

	  if (S1_i - TOLER > ZERO) {
	    S1_NS.allocate(ky_max);
	    S1_NS.settype(SPLINE2D_LINEAR);

	    s_south = S1_i * smax_south + BndSouthSpline.sp[0];
	    S1_NS.Xp[0] = Spline(s_south, BndSouthSpline);
	    S1_NS.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	    s_north = S1_i * smax_north + BndNorthSpline.sp[0],
	      S1_NS.Xp[ky_max-1] = Spline(s_north, BndNorthSpline);
	    S1_NS.tp[ky_max-1] = SPLINE2D_POINT_SHARP_CORNER;

	    for ( k = 1 ; k <= ky_max-2 ; ++k) {                     
	      s_east = double(k)*smax_east/double(ky_max-1) + BndEastSpline.sp[0];
	      x_east = Spline(s_east, BndEastSpline);

	      s_west = double(k)*smax_west/double(ky_max-1) + BndWestSpline.sp[0];
	      x_west = Spline(s_west, BndWestSpline);

	      w_west =  ONE-S1_i;
	      w_east =  S1_i;
	      w_total = w_east + w_west;
	      S1_NS.Xp[k] = (w_west*x_west + w_east*x_east)/w_total;
	      S1_NS.tp[k] = SPLINE2D_POINT_SHARP_CORNER;
	    } /* endfor */

	    S1_NS.pathlength();
	    smax_S1_NS = S1_NS.sp[S1_NS.np-1]-
	      S1_NS.sp[0];
	  } else {
	    S1_NS = BndWestSpline;
	    smax_S1_NS = S1_NS.sp[S1_NS.np-1]-
	      S1_NS.sp[0];
	  } /* endif */

	  if (S2_i + TOLER < ONE) {
	    S2_NS.allocate(ky_max);
	    S2_NS.settype(SPLINE2D_LINEAR);

	    s_south = S2_i * smax_south + BndSouthSpline.sp[0];
	    S2_NS.Xp[0] = Spline(s_south, BndSouthSpline);
	    S2_NS.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	    s_north = S2_i * smax_north + BndNorthSpline.sp[0],
	      S2_NS.Xp[ky_max-1] = Spline(s_north, BndNorthSpline);
	    S2_NS.tp[ky_max-1] = SPLINE2D_POINT_SHARP_CORNER;

	    for ( k = 1 ; k <= ky_max-2 ; ++k) {
	      s_east = double(k)*smax_east/double(ky_max-1) + BndEastSpline.sp[0];
	      x_east = Spline(s_east, BndEastSpline);

	      s_west = double(k)*smax_west/double(ky_max-1) + BndWestSpline.sp[0];
	      x_west = Spline(s_west, BndWestSpline);

	      w_west =  ONE-S2_i;
	      w_east =  S2_i;
	      w_total = w_east + w_west;
	      S2_NS.Xp[k] = (w_west*x_west + w_east*x_east)/w_total;
	      S2_NS.tp[k] = SPLINE2D_POINT_SHARP_CORNER;
	    } /* endfor */

	    S2_NS.pathlength();
	    smax_S2_NS = S2_NS.sp[S2_NS.np-1]-
	      S2_NS.sp[0];
	  } else {
	    S2_NS = BndEastSpline;
	    smax_S2_NS = S2_NS.sp[S2_NS.np-1]-
	      S2_NS.sp[0];
	  } /* endif */

	  if (S1_j - TOLER > ZERO) {
	    S1_EW.allocate(2);
	    S1_EW.settype(SPLINE2D_LINEAR);

	    s_west = S1_j * smax_west + BndWestSpline.sp[0];
	    S1_EW.Xp[0] = Spline(s_west, BndWestSpline);
	    S1_EW.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	    s_east = S1_j * smax_east + BndEastSpline.sp[0];
	    S1_EW.Xp[1] = Spline(s_east, BndEastSpline);
	    S1_EW.tp[1] = SPLINE2D_POINT_SHARP_CORNER;

	    S1_EW.pathlength();
	    smax_S1_EW = S1_EW.sp[S1_EW.np-1]-
	      S1_EW.sp[0];
	  } else {
	    S1_EW = BndSouthSpline;
	    smax_S1_EW = S1_EW.sp[S1_EW.np-1]-
	      S1_EW.sp[0];
	  } /* endif */

	  if (S2_j + TOLER < ONE) {
	    S2_EW.allocate(2);
	    S2_EW.settype(SPLINE2D_LINEAR);

	    s_west = S2_j * smax_west + BndWestSpline.sp[0];
	    S2_EW.Xp[0] = Spline(s_west, BndWestSpline);
	    S2_EW.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	    s_east = S2_j * smax_east + BndEastSpline.sp[0];
	    S2_EW.Xp[1] = Spline(s_east, BndEastSpline);
	    S2_EW.tp[1] = SPLINE2D_POINT_SHARP_CORNER;

	    S2_EW.pathlength();
	    smax_S2_EW = S2_EW.sp[S2_EW.np-1]-
	      S2_EW.sp[0];
	  } else {
	    S2_EW = BndNorthSpline;
	    smax_S2_EW = S2_EW.sp[S2_EW.np-1]-
	      S2_EW.sp[0];
	  } /* endif */

	  s_north = S_i * smax_S2_EW + S2_EW.sp[0];
	  x_north = Spline(s_north, S2_EW);
	  s_south = S_i * smax_S1_EW + S1_EW.sp[0];
	  x_south = Spline(s_south, S1_EW);

	  s_east = S_j * smax_S2_NS + S2_NS.sp[0];
	  x_east = Spline(s_east, S2_NS);
	  s_west = S_j * smax_S1_NS + S1_NS.sp[0];
	  x_west = Spline(s_west, S1_NS);
      
	  if (kx == 0 || kx == kx_max - 2) {
	    w_west =  ONE-(S_i-S1_i)/(S2_i-S1_i);
	    w_east =  (S_i-S1_i)/(S2_i-S1_i);
	    w_total = w_east + w_west;
	    Node[i][j].X = (w_west*x_west + 
			    w_east*x_east)/w_total;
	  } else {
	    w_south =  ONE-(S_j-S1_j)/(S2_j-S1_j);
	    w_north =  (S_j-S1_j)/(S2_j-S1_j);
	    w_total = w_south + w_north;
	    Node[i][j].X = (w_south*x_south + 
			    w_north*x_north)/w_total;
	  } /* endif */

	  S1_NS.deallocate(); S2_NS.deallocate();
	  S2_EW.deallocate(); S2_EW.deallocate();
	  break;
	  //============================================================
	default:
	  w_south = (ONE-s_j/smax_j);
	  w_north = (s_j/smax_j);
	  w_total = w_south + w_north;
	  Node[i][j].X = (w_south*x_south + 
			  w_north*x_north)/w_total;
	  break;
	} /* endswitch */

      } /* endif */

    } /* endfor */
  }/* endfor */

  /* Require update of the interior cells geometric properties. */
  Schedule_Interior_Mesh_Update();

  /* Set the boundary condition types at the quadrilateral 
     grid block boundaries. */

  Set_BCs();

}

/*!
 * Broadcast quadrilateral grid block to all            
 * processors involved in the calculation from the      
 * primary processor using the MPI broadcast routine.   
 *                                                      
 * \todo Merge some of the broadcasts.
 */
void Grid2D_Quad_Block_HO::Broadcast_Quad_Block(void) {

#ifdef _MPI_VERSION

  int i, j, ni, nj, ng, Highest_Order_of_Reconstruction, mesh_allocated, buffer_size, TD_Bcast, td, counter;
  double *buffer;
 
  /* Broadcast the number of cells in each direction. */

  if (CFFC_Primary_MPI_Processor()) {
    ni = NCi; 
    nj = NCj;
    ng = Nghost;
    Highest_Order_of_Reconstruction = HighestReconstructionOrder;
    if (Node != NULL) {
      mesh_allocated = 1;
    } else {
      mesh_allocated = 0;
    } /* endif */ 
  } /* endif */

  MPI::COMM_WORLD.Bcast(&ni, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&nj, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&ng, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&mesh_allocated, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&Highest_Order_of_Reconstruction, 1, MPI::INT, 0);

  /* On non-primary MPI processors, allocate (re-allocate) 
     memory for the cells and nodes of the quadrilateral 
     mesh block as necessary. */

  if (!CFFC_Primary_MPI_Processor()) {
    if (mesh_allocated) allocate(ni-2*ng, nj-2*ng, ng, Highest_Order_of_Reconstruction); 
  } /* endif */

  /* Broadcast the north, south, east, and west 
     boundary splines as well as their extensions. */

  BndNorthSpline.Broadcast_Spline();
  BndSouthSpline.Broadcast_Spline();
  BndEastSpline.Broadcast_Spline();
  BndWestSpline.Broadcast_Spline();

  ExtendWest_BndNorthSpline.Broadcast_Spline();
  ExtendEast_BndNorthSpline.Broadcast_Spline();
  ExtendWest_BndSouthSpline.Broadcast_Spline();
  ExtendEast_BndSouthSpline.Broadcast_Spline();
  ExtendNorth_BndEastSpline.Broadcast_Spline();
  ExtendSouth_BndEastSpline.Broadcast_Spline();
  ExtendNorth_BndWestSpline.Broadcast_Spline();
  ExtendSouth_BndWestSpline.Broadcast_Spline();

  /* Broadcast min/max pathlengths for boundary splines. */

  MPI::COMM_WORLD.Bcast(&(SminN), 1, MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&(SmaxN), 1, MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&(SminS), 1, MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&(SmaxS), 1, MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&(SminE), 1, MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&(SmaxE), 1, MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&(SminW), 1, MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&(SmaxW), 1, MPI::DOUBLE, 0);

  /* Broadcast node control parameters. */

  MPI::COMM_WORLD.Bcast(&(StretchI), 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&(BetaI), 1, MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&(TauI), 1, MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&(StretchJ), 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&(BetaJ), 1, MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&(TauJ), 1, MPI::DOUBLE, 0);

  MPI::COMM_WORLD.Bcast(&(OrthogonalN), 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&(OrthogonalS), 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&(OrthogonalE), 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&(OrthogonalW), 1, MPI::INT, 0);

  if (HO_Grid2D_Execution_Mode::USE_BROADCAST_MORE_THAN_RECOMPUTING) {

    /* Broadcast the node locations and some of the cell geometric 
       properties for grid block. */

    if (mesh_allocated) {
      ni = (INu+Nghost) - (INl-Nghost) + 1;
      nj = (JNu+Nghost) - (JNl-Nghost) + 1;

      /* Calculate how many geometric coefficients must be broadcast
	 by subtracting from their total number the first values, which 
	 is already known as being 1; 
      */
      TD_Bcast = CellGeomCoeff(0,0).size() - 1;

      /* Calculate the size of the buffer. 
	 First term relates to the nodal information while the second term
	 relates to the cell geometric properties (i.e. area, centroid, geometric moments). */
      buffer_size = 2*ni*nj + (3 + TD_Bcast)*(ni-1)*(nj-1);

      buffer = new double[buffer_size];

      // Load the buffer on the primary CPU
      if (CFFC_Primary_MPI_Processor()) {
	counter = 0;
	// Pack all the information except for the last North and West rows of nodes
	for (j  = JNl-Nghost ; j < JNu+Nghost ; ++j ) {
	  for ( i = INl-Nghost ; i < INu+Nghost ; ++i ) {
	    buffer[counter] = Node[i][j].X.x;
	    buffer[counter+1] = Node[i][j].X.y;
	    buffer[counter+2] = CellArea(i,j);
	    buffer[counter+3] = XCellCentroid(i,j);
	    buffer[counter+4] = YCellCentroid(i,j);
	    counter += 5;  // update buffer_size
	    for (td = 1; td<=TD_Bcast; ++td, ++counter){
	      buffer[counter] = CellGeomCoeffValue(i,j,td);
	    }
	  } /* endfor */
	} /* endfor */

	// Pack the row of nodes with i=INu+Nghost
	i = INu+Nghost;
	for (j = 0; j<=JNu+Nghost; ++j){
	  buffer[counter  ] = Node[i][j].X.x;
	  buffer[counter+1] = Node[i][j].X.y;
	  counter += 2;
	}

	// Pack the row of nodes with j=JNu+Nghost
	j = JNu+Nghost;
	for (i = 0; i<INu+Nghost; ++i){
	  buffer[counter  ] = Node[i][j].X.x;
	  buffer[counter+1] = Node[i][j].X.y;
	  counter += 2;
	}
      } /* endif */

      // Broadcast buffer
      MPI::COMM_WORLD.Bcast(buffer, buffer_size, MPI::DOUBLE, 0);

      // Unload the buffer on the receiver CPUs and set the variables
      if (!CFFC_Primary_MPI_Processor()) {
	counter = 0;
	// Unpack all the information except for the last North and West rows of nodes
	for (j  = JNl-Nghost; j < JNu+Nghost; ++j ) {
	  for ( i = INl-Nghost; i < INu+Nghost; ++i ) {
	    Node[i][j].X.x  = buffer[counter  ];
	    Node[i][j].X.y  = buffer[counter+1];
	    Cell[i][j].A    = buffer[counter+2];
	    Cell[i][j].Xc.x = buffer[counter+3];
	    Cell[i][j].Xc.y = buffer[counter+4];
	    counter += 5;  // update buffer_size

	    // set the first geometric moment
	    Cell[i][j].GeomCoeffValue(0) = 1.0;
	    // set the rest of the coefficients to the transferred values
	    for (td = 1; td<=TD_Bcast; ++td, ++counter){
	      Cell[i][j].GeomCoeffValue(td) = buffer[counter];
	    }
	  } /* endfor */
	} /* endfor */

	// Unpack the row of nodes with i=INu+Nghost
	i = INu+Nghost;
	for (j = 0; j<=JNu+Nghost; ++j){
	  Node[i][j].X.x = buffer[counter  ];
	  Node[i][j].X.y = buffer[counter+1];
	  counter += 2;
	}

	// Unpack the row of nodes with j=JNu+Nghost
	j = JNu+Nghost;
	for (i = 0; i<INu+Nghost; ++i){
	  Node[i][j].X.x = buffer[counter  ];
	  Node[i][j].X.y = buffer[counter+1];
	  counter += 2;
	}

	/* On non-source MPI processors, set the boundary condition types
	   and update the 2D spline info(s) for the quadrilateral mesh block. */
	Set_BCs();
	Update_SplineInfos();
      } /* endif */
      
      delete []buffer; 
      buffer = NULL; 
    }/* endif */

  } else { 

    /* Broadcast the node locations for grid block. */
    if (mesh_allocated) {
      ni = (INu+Nghost) - (INl-Nghost) + 1;
      nj = (JNu+Nghost) - (JNl-Nghost) + 1;
      buffer = new double[2*ni*nj];

      if (CFFC_Primary_MPI_Processor()) {
	buffer_size = 0;
	for (j  = JNl-Nghost; j <= JNu+Nghost; ++j ) {
	  for ( i = INl-Nghost; i <= INu+Nghost; ++i ) {
	    buffer[buffer_size] = Node[i][j].X.x;
	    buffer[buffer_size+1] = Node[i][j].X.y;
	    buffer_size = buffer_size + 2;
	  } /* endfor */
	} /* endfor */
      } /* endif */

      buffer_size = 2*ni*nj;
      MPI::COMM_WORLD.Bcast(buffer, buffer_size, MPI::DOUBLE, 0);

      if (!CFFC_Primary_MPI_Processor()) {
	buffer_size = 0;
	for (j  = JNl-Nghost; j <= JNu+Nghost; ++j ) {
	  for ( i = INl-Nghost; i <= INu+Nghost; ++i ) {
	    Node[i][j].X.x = buffer[buffer_size];
	    Node[i][j].X.y = buffer[buffer_size+1];
	    buffer_size = buffer_size + 2;
	  } /* endfor */
	} /* endfor */

	/* Require update of the whole mesh */
	Schedule_Interior_Mesh_Update();
	Schedule_Ghost_Cells_Update();

	/* On non-primary MPI processors, set the boundary condition types
	   and compute the cells for the quadrilateral mesh block. */
	Set_BCs();
	Update_Cells();
      } /* endif */

      delete []buffer;
      buffer = NULL;
    } /* endif */

  } /* endif */

  /* Broadcast state trackers. */
  MPI::COMM_WORLD.Bcast(&(InteriorCellGeometryStateTracker), 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&(GhostCellGeometryStateTracker), 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&(CornerGhostCellGeometryStateTracker), 1, MPI::INT, 0);
  
  /* On non-primary MPI processors, mark the current geometry 
     different than what it was stored before. */

  if (!CFFC_Primary_MPI_Processor()) {
    New_Global_Geometry_State();
  }
  
#endif

}

#ifdef _MPI_VERSION
/*!
 * Broadcast quadrilateral grid block to all processors 
 * associated with the specified communicator from the  
 * specified processor using the MPI broadcast routine. 
 *
 * \todo Merge some of the broadcasts.
 */
void Grid2D_Quad_Block_HO::Broadcast_Quad_Block(MPI::Intracomm &Communicator, 
						const int Source_CPU) {
  
  int Source_Rank = 0;
  int i, j, ni, nj, ng, Highest_Order_of_Reconstruction, mesh_allocated, buffer_size, TD_Bcast, counter, td;
  double *buffer;
  
  /* Broadcast the number of cells in each direction. */
  
  if (CFFC_MPI::This_Processor_Number == Source_CPU) {
    ni = NCi;
    nj = NCj;
    ng = Nghost;
    Highest_Order_of_Reconstruction = HighestReconstructionOrder;
    if (Node != NULL) {
      mesh_allocated = 1;
    } else {
      mesh_allocated = 0;
    } /* endif */ 
  } /* endif */

  Communicator.Bcast(&ni, 1, MPI::INT, Source_Rank);
  Communicator.Bcast(&nj, 1, MPI::INT, Source_Rank);
  Communicator.Bcast(&ng, 1, MPI::INT, Source_Rank);
  Communicator.Bcast(&mesh_allocated, 1, MPI::INT, Source_Rank);
  Communicator.Bcast(&Highest_Order_of_Reconstruction, 1, MPI::INT, Source_Rank);
  
  /* On non-source MPI processors, allocate (re-allocate) 
     memory for the cells and nodes of the quadrilateral 
     mesh block as necessary. */
  
  if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
    if (mesh_allocated) allocate(ni-2*ng, nj-2*ng, ng, Highest_Order_of_Reconstruction); 
  }/* endif */
  
  /* Broadcast the north, south, east, and west 
     boundary splines. */
  
  BndNorthSpline.Broadcast_Spline(Communicator, Source_CPU);
  BndSouthSpline.Broadcast_Spline(Communicator, Source_CPU);
  BndEastSpline.Broadcast_Spline(Communicator, Source_CPU);
  BndWestSpline.Broadcast_Spline(Communicator, Source_CPU);

  ExtendWest_BndNorthSpline.Broadcast_Spline(Communicator, Source_CPU);
  ExtendEast_BndNorthSpline.Broadcast_Spline(Communicator, Source_CPU);
  ExtendWest_BndSouthSpline.Broadcast_Spline(Communicator, Source_CPU);
  ExtendEast_BndSouthSpline.Broadcast_Spline(Communicator, Source_CPU);
  ExtendNorth_BndEastSpline.Broadcast_Spline(Communicator, Source_CPU);
  ExtendSouth_BndEastSpline.Broadcast_Spline(Communicator, Source_CPU);
  ExtendNorth_BndWestSpline.Broadcast_Spline(Communicator, Source_CPU);
  ExtendSouth_BndWestSpline.Broadcast_Spline(Communicator, Source_CPU);
  
  /* Broadcast min/max pathlengths for boundary splines. */
  
  Communicator.Bcast(&(SminN), 1, MPI::DOUBLE, Source_Rank);
  Communicator.Bcast(&(SmaxN), 1, MPI::DOUBLE, Source_Rank);
  Communicator.Bcast(&(SminS), 1, MPI::DOUBLE, Source_Rank);
  Communicator.Bcast(&(SmaxS), 1, MPI::DOUBLE, Source_Rank);
  Communicator.Bcast(&(SminE), 1, MPI::DOUBLE, Source_Rank);
  Communicator.Bcast(&(SmaxE), 1, MPI::DOUBLE, Source_Rank);
  Communicator.Bcast(&(SminW), 1, MPI::DOUBLE, Source_Rank);
  Communicator.Bcast(&(SmaxW), 1, MPI::DOUBLE, Source_Rank);
  
  /* Broadcast node control parameters. */

  Communicator.Bcast(&(StretchI), 1, MPI::INT, Source_Rank);
  Communicator.Bcast(&(BetaI), 1, MPI::DOUBLE, Source_Rank);
  Communicator.Bcast(&(TauI), 1, MPI::DOUBLE, Source_Rank);
  Communicator.Bcast(&(StretchJ), 1, MPI::INT, Source_Rank);
  Communicator.Bcast(&(BetaJ), 1, MPI::DOUBLE, Source_Rank);
  Communicator.Bcast(&(TauJ), 1, MPI::DOUBLE, Source_Rank);

  Communicator.Bcast(&(OrthogonalN), 1, MPI::INT, Source_Rank);
  Communicator.Bcast(&(OrthogonalS), 1, MPI::INT, Source_Rank);
  Communicator.Bcast(&(OrthogonalE), 1, MPI::INT, Source_Rank);
  Communicator.Bcast(&(OrthogonalW), 1, MPI::INT, Source_Rank);

  if (HO_Grid2D_Execution_Mode::USE_BROADCAST_MORE_THAN_RECOMPUTING) {

    /* Broadcast the node locations and some of the cell geometric 
       properties for grid block. */

    if (mesh_allocated) {
      ni = (INu+Nghost) - (INl-Nghost) + 1;
      nj = (JNu+Nghost) - (JNl-Nghost) + 1;

      /* Calculate how many geometric coefficients must be broadcast
	 by subtracting from their total number the first values, which 
	 is already known as being 1; 
      */
      TD_Bcast = CellGeomCoeff(0,0).size() - 1;

      /* Calculate the size of the buffer. 
	 First term relates to the nodal information while the second term
	 relates to the cell geometric properties (i.e. area, centroid, geometric moments). */
      buffer_size = 2*ni*nj + (3 + TD_Bcast)*(ni-1)*(nj-1);

      buffer = new double[buffer_size];

      // Load the buffer on the source CPU
      if (CFFC_MPI::This_Processor_Number == Source_CPU) {
	counter = 0;
	// Pack all the information except for the last North and West rows of nodes
	for (j  = JNl-Nghost ; j < JNu+Nghost ; ++j ) {
	  for ( i = INl-Nghost ; i < INu+Nghost ; ++i ) {
	    buffer[counter] = Node[i][j].X.x;
	    buffer[counter+1] = Node[i][j].X.y;
	    buffer[counter+2] = CellArea(i,j);
	    buffer[counter+3] = XCellCentroid(i,j);
	    buffer[counter+4] = YCellCentroid(i,j);
	    counter += 5;  // update buffer_size
	    for (td = 1; td<=TD_Bcast; ++td, ++counter){
	      buffer[counter] = CellGeomCoeffValue(i,j,td);
	    }
	  } /* endfor */
	} /* endfor */

	// Pack the row of nodes with i=INu+Nghost
	i = INu+Nghost;
	for (j = 0; j<=JNu+Nghost; ++j){
	  buffer[counter  ] = Node[i][j].X.x;
	  buffer[counter+1] = Node[i][j].X.y;
	  counter += 2;
	}

	// Pack the row of nodes with j=JNu+Nghost
	j = JNu+Nghost;
	for (i = 0; i<INu+Nghost; ++i){
	  buffer[counter  ] = Node[i][j].X.x;
	  buffer[counter+1] = Node[i][j].X.y;
	  counter += 2;
	}
      } /* endif */

      // Broadcast buffer
      Communicator.Bcast(buffer, buffer_size, MPI::DOUBLE, Source_Rank);

      // Unload the buffer on the receiver CPUs and set the variables
      if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
	counter = 0;
	// Unpack all the information except for the last North and West rows of nodes
	for (j  = JNl-Nghost; j < JNu+Nghost; ++j ) {
	  for ( i = INl-Nghost; i < INu+Nghost; ++i ) {
	    Node[i][j].X.x  = buffer[counter  ];
	    Node[i][j].X.y  = buffer[counter+1];
	    Cell[i][j].A    = buffer[counter+2];
	    Cell[i][j].Xc.x = buffer[counter+3];
	    Cell[i][j].Xc.y = buffer[counter+4];
	    counter += 5;  // update buffer_size

	    // set the first geometric moment
	    Cell[i][j].GeomCoeffValue(0) = 1.0;
	    // set the rest of the coefficients to the transferred values
	    for (td = 1; td<=TD_Bcast; ++td, ++counter){
	      Cell[i][j].GeomCoeffValue(td) = buffer[counter];
	    }
	  } /* endfor */
	} /* endfor */

	// Unpack the row of nodes with i=INu+Nghost
	i = INu+Nghost;
	for (j = 0; j<=JNu+Nghost; ++j){
	  Node[i][j].X.x = buffer[counter  ];
	  Node[i][j].X.y = buffer[counter+1];
	  counter += 2;
	}

	// Unpack the row of nodes with j=JNu+Nghost
	j = JNu+Nghost;
	for (i = 0; i<INu+Nghost; ++i){
	  Node[i][j].X.x = buffer[counter  ];
	  Node[i][j].X.y = buffer[counter+1];
	  counter += 2;
	}

	/* On non-source MPI processors, set the boundary condition types
	   and update the 2D spline info(s) for the quadrilateral mesh block. */
	Set_BCs();
	Update_SplineInfos();
      } /* endif */
      
      delete []buffer; 
      buffer = NULL; 
    }/* endif */

  } else { 

    /* Broadcast the node locations for grid block. */
    if (mesh_allocated) {
      ni = (INu+Nghost) - (INl-Nghost) + 1;
      nj = (JNu+Nghost) - (JNl-Nghost) + 1;
      buffer = new double[2*ni*nj];

      if (CFFC_MPI::This_Processor_Number == Source_CPU) {
	buffer_size = 0;
	for (j  = JNl-Nghost ; j <= JNu+Nghost ; ++j ) {
	  for ( i = INl-Nghost ; i <= INu+Nghost ; ++i ) {
	    buffer[buffer_size] = Node[i][j].X.x;
	    buffer[buffer_size+1] = Node[i][j].X.y;
	    buffer_size = buffer_size + 2;
	  } /* endfor */
	} /* endfor */
      } /* endif */
    
      buffer_size = 2*ni*nj;
      Communicator.Bcast(buffer, buffer_size, MPI::DOUBLE, Source_Rank);
    
      if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
	buffer_size = 0;
	for (j  = JNl-Nghost; j <= JNu+Nghost; ++j ) {
	  for ( i = INl-Nghost; i <= INu+Nghost; ++i ) {
	    Node[i][j].X.x = buffer[buffer_size];
	    Node[i][j].X.y = buffer[buffer_size+1];
	    buffer_size = buffer_size + 2;
	  } /* endfor */
	} /* endfor */

	/* Require update of the whole mesh */
	Schedule_Interior_Mesh_Update();
	Schedule_Ghost_Cells_Update();

	/* On non-source MPI processors, set the boundary condition types
	   and compute the cells for the quadrilateral mesh block. */
	Set_BCs();
	Update_Cells();
      } /* endif */
    
      delete []buffer; 
      buffer = NULL;
    } /* endif */

  }/* endif */

  /* Broadcast state trackers. */
  Communicator.Bcast(&(InteriorCellGeometryStateTracker), 1, MPI::INT, Source_Rank);
  Communicator.Bcast(&(GhostCellGeometryStateTracker), 1, MPI::INT, Source_Rank);
  Communicator.Bcast(&(CornerGhostCellGeometryStateTracker), 1, MPI::INT, Source_Rank);
  
  /* On non-source MPI processors, mark the current geometry 
     different than what it was stored before. */
  if (CFFC_MPI::This_Processor_Number != Source_CPU) {
    New_Global_Geometry_State();
  }
}
#endif

/*!
 * Smooths quadrilateral grid block using an elliptic   
 * grid generation method (transformed Poisson's        
 * equation method) with a successive over-relaxation   
 * (SOR) Gauss-Seidel solution procedure.  The routine  
 * follows the technique outlined by Sorenson (1980)    
 * and enforces fixed grid spacing and orthogonality at  
 * the bottom (south) and top (north) boundaries of the 
 * quadrilateral block as required.  Modifications have 
 * been introduced to permit the enforcement of grid    
 * orthogonality at the left (west) and right (east)    
 * boundaries as well.                                  
 *                                                      
 */
void Grid2D_Quad_Block_HO::Smooth_Quad_Block(const int Number_of_Iterations) {

  // Run smoother only is Smooth_Quad_Block_Flag is ON
  if (Smooth_Quad_Block_Flag){
  
    /*
     *  Local Variable description:
     * 
     *  xij, yij                  Two-dimensional arrays containing
     *                            the coordinates of the nodes for
     *                            the grid block in physical (x,y) space.
     *
     *  deltas_b, deltas_t,       One-dimensional arrays containing the
     *  deltas_l, deltas_r        node spacing in the normal direction
     *                            between the first two nodes at the upper
     *                            and lower boundaries of the grid block.
     *
     *  xz_b, yz_b, xe_b, ye_b,   One-dimensional arrays containing the
     *  xzz_b, yzz_b, xee_b,      partial derivatives of x and y
     *  yee_b, xze_b, yze_b       coordinates of the physical space
     *                            with respect to the zeta and eta
     *  xz_t, yz_t, xe_t, ye_t,   coordinates of the computational space,
     *  xzz_t, yzz_t, xee_t,      all evaluated at the bottom (south),
     *  yee_t, xze_t, yze_t       top (north), left (west), and right 
     *                            (east) boundaries of the grid block.
     *
     *  xz_l, yz_l, xe_l, ye_l,
     *  xzz_l, yzz_l, xee_l,   
     *  yee_l, xze_l, yze_l
     *
     *  xz_r, yz_r, xe_r, ye_r,
     *  xzz_r, yzz_r, xee_r,   
     *  yee_r, xze_r, yze_r
     *
     *  jacob_b, jacob_t,         One-dimensional arrays containing the
     *  jacob_l, jacob_r          Jacobian of the coordinate transformation
     *                            evaluated at the bottom (south), top 
     *                            (north), left (west), and right (east)
     *                            boundaries of the grid block.
     *
     *  pb, pt, pl, pr,           One-dimensional arrays containing the
     *  qb, qt, ql, qr            values of the source terms of the
     *                            transformed Poisson's equations at the
     *                            bottom (south), top (north), left (west), 
     *                            and right (east) boundaries of the grid 
     *                            block.
     *
     *  norm_b, norm_t,           Integer parameters indicating whether 
     *  norm_l, norm_r            orthogonality at the bottom (south), top
     *                            (north), left (west), and right (east) 
     *                            boundaries is to be enforced.  If value 
     *                            is zero then orthogonality is not enforced 
     *                            and any other value causes the orthogonality 
     *                            condition to be applied.
     *
     *  wsor, wsur                SOR and successive under-relaxation (SUR)
     *                            iteration parameters.  wsor is used to
     *                            over-relax the solution for x and y and
     *                            wsur is used to under-relax the solution
     *                            for p, q, r, and s.  Typically, 
     *                            wsor=1.5-1.8 and wsur=0.30-0.70.
     *
     *  a, b,                     Exponential decay parameters defining the
     *  c, d,                     effects of the source terms away from the
     *  e                         bottom (south), top (north), left (west), 
     *                            and right (east) boundaries, as well as the 
     *                            corners.  If a (or b, c, d, e) is small 
     *                            (i.e., a = 0.2) source terms have effects far 
     *                            from boundary.  If a (or b, c, d, e) is large
     *                            (i.e., a = 0.7) source terms have less 
     *                            effects far from boundary.
     *
     *  n_zeta, n_eta             Number of nodes in the zeta and eta
     *                            directions of the computational domain for 
     *                            grid block.
     *
     * alpha, beta, gamma        Parameters used in the transformed Poisson's
     *                            equations.
     * r1, r2, r3, r4            Parameters used to evaluate the source terms
     *                            of the transformed Poisson's equations at
     *                            the bottom (south), top (north), left (west), 
     *                            and right (east) boundaries of the grid block.
     *
     * dpb, dpt, dpl, dpr,       Changes in pb, pt, pl, pr, qb, qt, ql, and qr
     * dqb, dqt, dql, dqr        before the under-relaxation and limiting
     *                           algorithm is applied.
     *
     * dpmax, dqmax              Maximum permitted changes in p and q at the
     *                           bottom, top, left, and right boundaries.
     *
     * xzij, yzij, xeij, yeij,   Local values of the transformation metrics
     * jacobij                   and Jacobian.
     *
     * pij, qij                  Local values of the source terms.
     *
     * dxij, dyij                Changes in xij and yij for on step in the
     *                           iteration technique before the over-
     *                           relaxation algorithm is applied.
     *
     * fb, ft, fl, fr,           Exponentials use to evaluate the local
     * fbl, fbr, ftl, ftr        values of the source terms. 
     *
     */

    double **xij, **yij,
      *deltas_b, *deltas_t, *deltas_l, *deltas_r,
      *xz_b, *yz_b, *xe_b, *ye_b, 
      *xzz_b, *yzz_b, *xee_b, *yee_b, *xze_b, *yze_b,
      *xz_t, *yz_t, *xe_t, *ye_t,
      *xzz_t, *yzz_t, *xee_t, *yee_t, *xze_t, *yze_t,
      *xz_l, *yz_l, *xe_l, *ye_l,
      *xzz_l, *yzz_l, *xee_l, *yee_l, *xze_l, *yze_l,
      *xz_r, *yz_r, *xe_r, *ye_r,
      *xzz_r, *yzz_r, *xee_r, *yee_r, *xze_r, *yze_r,
      *pb, *pt, *pl, *pr,
      *qb, *qt, *ql, *qr,
      *jacob_b, *jacob_t, *jacob_l, *jacob_r;

    int num_iter, i_poisson, 
      norm_b, norm_t, norm_l, norm_r,
      n_zeta, n_eta;
    double wsor, wsur, a, b, c, d, e, l_damping;

    int i, j, i_step, ii, jj;
    double alpha, beta, gamma;
    double r1, r2, r3, r4;
    double dpb, dpt, dpl, dpr, dqb, dqt, dql, dqr;
    double dpmax, dqmax;
    double xzij, yzij, xeij, yeij, jacobij, pij, qij, dxij, dyij;
    double fb, ft, fl, fr, fbl, fbr, ftl, ftr;

    /* Determine the size of and create the Poisson equation solution 
       variables for the elliptic smoothing algorithm. */

    n_zeta = INu - INl + 1;
    n_eta = JNu - JNl + 1;

    xij = new double*[n_zeta];
    yij = new double*[n_zeta];
    for ( i = 0; i <= n_zeta-1 ; ++i ) {
      xij[i] = new double[n_eta];
      yij[i] = new double[n_eta];
    } /* endfor */

    deltas_b = new double[n_zeta];
    deltas_t = new double[n_zeta];
    deltas_l = new double[n_eta];
    deltas_r = new double[n_eta];

    xz_b = new double[n_zeta]; 
    yz_b = new double[n_zeta];
    xe_b = new double[n_zeta];
    ye_b = new double[n_zeta];
    xzz_b = new double[n_zeta];
    yzz_b = new double[n_zeta];
    xee_b = new double[n_zeta];
    yee_b = new double[n_zeta];
    xze_b = new double[n_zeta];
    yze_b = new double[n_zeta];

    xz_t = new double[n_zeta];
    yz_t = new double[n_zeta];
    xe_t = new double[n_zeta]; 
    ye_t = new double[n_zeta];
    xzz_t = new double[n_zeta];
    yzz_t = new double[n_zeta];
    xee_t = new double[n_zeta];
    yee_t = new double[n_zeta];
    xze_t = new double[n_zeta];
    yze_t = new double[n_zeta];

    xz_l = new double[n_eta];
    yz_l = new double[n_eta];
    xe_l = new double[n_eta];
    ye_l = new double[n_eta];
    xzz_l = new double[n_eta];
    yzz_l = new double[n_eta];
    xee_l = new double[n_eta];
    yee_l = new double[n_eta];
    xze_l = new double[n_eta];
    yze_l = new double[n_eta];

    xz_r = new double[n_eta]; 
    yz_r = new double[n_eta];
    xe_r = new double[n_eta];
    ye_r = new double[n_eta];
    xzz_r = new double[n_eta];
    yzz_r = new double[n_eta];
    xee_r = new double[n_eta];
    yee_r = new double[n_eta];
    xze_r = new double[n_eta];
    yze_r = new double[n_eta];

    pb = new double[n_zeta];
    pt = new double[n_zeta];
    pl = new double[n_eta];
    pr = new double[n_eta];

    qb = new double[n_zeta]; 
    qt = new double[n_zeta]; 
    ql = new double[n_eta];
    qr = new double[n_eta];

    jacob_b = new double[n_zeta];
    jacob_t = new double[n_zeta];
    jacob_l = new double[n_eta];
    jacob_r = new double[n_eta];

    /* Assign initial values to the Poisson equation solution
       variables using the current node locations for the nodes
       in the quadrilateral grid block. */

    for (j  = JNl ; j <= JNu ; ++j ) {
      for ( i = INl ; i <= INu ; ++i ) {
	ii = i - INl;
	jj = j - JNl;
	xij[ii][jj] = Node[i][j].X.x;
	yij[ii][jj] = Node[i][j].X.y;
      } /* endfor */
    } /* endfor */

    /* Set the Poisson equation solution parameters. */

    norm_b = OrthogonalS;
    norm_t = OrthogonalN;
    norm_l = OrthogonalW;
    norm_r = OrthogonalE;

    wsor = 1.25;
    wsur = 0.10;
    l_damping = 1.00;

    if (norm_b == 0 && norm_t == 0 &&
	norm_l == 0 && norm_r == 0) {
      i_poisson = 0;
      a = ONE;
      b = ONE;
      c = ONE;
      d = ONE;
      e = HALF;
    } else if (norm_b != 0 && norm_t == 0 &&
	       norm_l == 0 && norm_r == 0) {
      i_poisson = 1;
      a = l_damping;
      b = ONE;
      c = ONE;
      d = ONE;
      e = HALF;
    } else if (norm_b == 0 && norm_t != 0 &&
	       norm_l == 0 && norm_r == 0) {
      i_poisson = 1;
      a = ONE;
      b = l_damping;
      c = ONE;
      d = ONE;
      e = HALF;
    } else if (norm_b == 0 && norm_t == 0 &&
	       norm_l != 0 && norm_r == 0) {
      i_poisson = 1;
      a = ONE;
      b = ONE;
      c = l_damping;
      d = ONE;
      e = HALF;
    } else if (norm_b == 0 && norm_t == 0 &&
	       norm_l == 0 && norm_r != 0) {
      i_poisson = 1;
      a = ONE;
      b = ONE;
      c = ONE;
      d = l_damping;
      e = HALF;
    } else if (norm_b != 0 && norm_t != 0 &&
	       norm_l == 0 && norm_r == 0) {
      i_poisson = 1;
      a = l_damping;
      b = l_damping;
      c = ONE;
      d = ONE;
      e = HALF;
    } else if (norm_b == 0 && norm_t == 0 &&
	       norm_l != 0 && norm_r != 0) {
      i_poisson = 1;
      a = ONE;
      b = ONE;
      c = l_damping;
      d = l_damping;
      e = HALF;
    } else if (norm_b != 0 && norm_t == 0 &&
	       norm_l != 0 && norm_r == 0) {
      i_poisson = 1;
      a = l_damping;
      b = ONE;
      c = l_damping;
      d = ONE;
      e = HALF;
    } else if (norm_b == 0 && norm_t != 0 &&
	       norm_l == 0 && norm_r != 0) {
      i_poisson = 1;
      a = ONE;
      b = l_damping;
      c = ONE;
      d = l_damping;
      e = HALF;
    } else if (norm_b == 0 && norm_t != 0 &&
	       norm_l != 0 && norm_r != 0) {
      i_poisson = 1;
      a = ONE;
      b = l_damping;
      c = l_damping;
      d = l_damping;
      e = HALF;
    } else if (norm_b != 0 && norm_t == 0 &&
	       norm_l != 0 && norm_r != 0) {
      i_poisson = 1;
      a = l_damping;
      b = ONE;
      c = l_damping;
      d = l_damping;
      e = HALF;
    } else if (norm_b != 0 && norm_t != 0 &&
	       norm_l == 0 && norm_r != 0) {
      i_poisson = 1;
      a = l_damping;
      b = l_damping;
      c = ONE;
      d = l_damping;
      e = HALF;
    } else if (norm_b != 0 && norm_t != 0 &&
	       norm_l != 0 && norm_r == 0) {
      i_poisson = 1;
      a = l_damping;
      b = l_damping;
      c = l_damping;
      d = ONE;
      e = HALF;
    } else {
      i_poisson = 1;
      a = l_damping;
      b = l_damping;
      c = l_damping;
      d = l_damping;
      e = HALF;
    } /* endif */

    /* Initialize the partial derivatives and source terms at the 
       bottom (south), top (north), left (west), and right (east) 
       boundaries.  The procedure for doing this depends on whether 
       or not orthogonality is enforced at the bottom (south), 
       top (north), left (west), and right (east) boundaries. */

    /* BOTTOM AND TOP */

    xz_b[0] = HALF * (xij[1][0] - xij[0][0]);
    xz_t[0] = HALF * (xij[1][n_eta-1] - xij[0][n_eta-1]);

    yz_b[0] = HALF * (yij[1][0] - yij[0][0]);
    yz_t[0] = HALF * (yij[1][n_eta-1] - yij[0][n_eta-1]);

    xe_b[0] = ZERO;
    xe_t[0] = ZERO;

    ye_b[0] = ZERO;
    ye_t[0] = ZERO;

    for ( i = 1; i < n_zeta - 1; ++i ) {
      xz_b[i] = HALF * (xij[i+1][0] - xij[i-1][0]);
      xz_t[i] = HALF * (xij[i+1][n_eta-1] - xij[i-1][n_eta-1]);
 
      yz_b[i] = HALF * (yij[i+1][0] - yij[i-1][0]);
      yz_t[i] = HALF * (yij[i+1][n_eta-1] - yij[i-1][n_eta-1]);

      xzz_b[i] = xij[i+1][0] - TWO * xij[i][0] + xij[i-1][0];
      xzz_t[i] = xij[i+1][n_eta-1] - TWO * xij[i][n_eta-1] +
	xij[i-1][n_eta-1];

      yzz_b[i] = yij[i+1][0] - TWO * yij[i][0] + yij[i-1][0];
      yzz_t[i] = yij[i+1][n_eta-1] - TWO * yij[i][n_eta-1] +
	yij[i-1][n_eta-1];
 
      if (norm_b == 0) {
	xe_b[i] = xij[i][1] - xij[i][0];
	ye_b[i] = yij[i][1] - yij[i][0];
      } else {
	deltas_b[i] = hypot( (xij[i][1]-xij[i][0]),
			     (yij[i][1]-yij[i][0]) );
	xe_b[i] = - deltas_b[i] * yz_b[i] /
	  hypot( xz_b[i], yz_b[i] );
	ye_b[i] = deltas_b[i] * xz_b[i] /
	  hypot( xz_b[i], yz_b[i] );
      } /* endif */

      if (norm_t == 0) {
	xe_t[i] = xij[i][n_eta-1] - xij[i][n_eta-2];
	ye_t[i] = yij[i][n_eta-1] - yij[i][n_eta-2];
      } else {
	deltas_t[i] = hypot( (xij[i][n_eta-1]-xij[i][n_eta-2]),
			     (yij[i][n_eta-1]-yij[i][n_eta-2]) );
	xe_t[i] = - deltas_t[i] * yz_t[i] /
	  hypot( xz_t[i], yz_t[i] );
	ye_t[i] = deltas_t[i] * xz_t[i] /
	  hypot( xz_t[i], yz_t[i] );
      } /* endif */      

      jacob_b[i] = xz_b[i] * ye_b[i] - xe_b[i] * yz_b[i];
      jacob_t[i] = xz_t[i] * ye_t[i] - xe_t[i] * yz_t[i];

      pb[i] = ZERO;
      pt[i] = ZERO;

      qb[i] = ZERO;
      qt[i] = ZERO;
    } /* endfor */

    xz_b[n_zeta-1] = HALF * (xij[n_zeta-1][0] - xij[n_zeta-2][0]);
    xz_t[n_zeta-1] = HALF * (xij[n_zeta-1][n_eta-1] -
			     xij[n_zeta-2][n_eta-1]);

    yz_b[n_zeta-1] = HALF * (yij[n_zeta-1][0] - yij[n_zeta-2][0]);
    yz_t[n_zeta-1] = HALF * (yij[n_zeta-1][n_eta-1] -
			     yij[n_zeta-2][n_eta-1]);

    xe_b[n_zeta-1] = ZERO;
    xe_t[n_zeta-1] = ZERO;

    ye_b[n_zeta-1] = ZERO;
    ye_t[n_zeta-1] = ZERO;

    for ( i = 1; i < n_zeta - 1; ++i ) {
      xze_b[i] = HALF * (xe_b[i+1] - xe_b[i-1]);
      xze_t[i] = HALF * (xe_t[i+1] - xe_t[i-1]);
 
      yze_b[i] = HALF * (ye_b[i+1] - ye_b[i-1]);
      yze_t[i] = HALF * (ye_t[i+1] - ye_t[i-1]);
    } /* endfor */

    /* LEFT AND RIGHT */

    xe_l[0] = HALF * (xij[0][1] - xij[0][0]);
    xe_r[0] = HALF * (xij[n_zeta-1][1] - xij[n_zeta-1][0]);

    ye_l[0] = HALF * (yij[0][1] - yij[0][0]);
    ye_r[0] = HALF * (yij[n_zeta-1][1] - yij[n_zeta-1][0]);

    xz_l[0] = ZERO;
    xz_r[0] = ZERO;

    yz_l[0] = ZERO;
    yz_r[0] = ZERO;

    for ( j = 1; j < n_eta - 1; ++j ) {
      xe_l[j] = HALF * (xij[0][j+1] - xij[0][j-1]);
      xe_r[j] = HALF * (xij[n_zeta-1][j+1] - xij[n_zeta-1][j-1]);

      ye_l[j] = HALF * (yij[0][j+1] - yij[0][j-1]);
      ye_r[j] = HALF * (yij[n_zeta-1][j+1] - yij[n_zeta-1][j-1]);

      xee_l[j] = xij[0][j+1] - TWO * xij[0][j] + xij[0][j-1];
      xee_r[j] = xij[n_zeta-1][j+1] - TWO * xij[n_zeta-1][j] +
	xij[n_zeta-1][j-1];

      yee_l[j] = yij[0][j+1] - TWO * yij[0][j] + yij[0][j-1];
      yee_r[j] = yij[n_zeta-1][j+1] - TWO * yij[n_zeta-1][j] +
	yij[n_zeta-1][j-1];

      if (norm_l == 0) {
	xz_l[j] = xij[1][j] - xij[0][j];
	yz_l[j] = yij[1][j] - yij[0][j];
      } else {
	deltas_l[j] = hypot( (xij[1][j]-xij[0][j]),
			     (yij[1][j]-yij[0][j]) );
	xz_l[j] = deltas_l[j] * ye_l[j] /
	  hypot( xe_l[j], ye_l[j] );
	yz_l[j] = - deltas_l[j] * xe_l[j] /
	  hypot( xe_l[j], ye_l[j] );
      } /* endif */

      if (norm_r == 0) {
	xz_r[j] = xij[n_zeta-1][j] - xij[n_zeta-2][j];
	yz_r[j] = yij[n_zeta-1][j] - yij[n_zeta-2][j];
      } else {
	deltas_r[j] = hypot( (xij[n_zeta-1][j]-xij[n_zeta-2][j]),
			     (yij[n_zeta-1][j]-yij[n_zeta-2][j]) );
	xz_r[j] = deltas_r[j] * ye_r[j] /
	  hypot( xe_r[j], ye_r[j] );
	yz_r[j] = -deltas_r[j] * xe_r[j] /
	  hypot( xe_r[j], ye_r[j] );
      } /* endif */

      jacob_l[j] = xz_l[j] * ye_l[j] - xe_l[j] * yz_l[j];
      jacob_r[j] = xz_r[j] * ye_r[j] - xe_r[j] * yz_r[j];

      pl[j] = ZERO;
      pr[j] = ZERO;

      ql[j] = ZERO;
      qr[j] = ZERO;
    } /* endfor */

    xe_l[n_eta-1] = HALF * (xij[0][n_eta-1] - xij[0][n_eta-2]);
    xe_r[n_eta-1] = HALF * (xij[n_zeta-1][n_eta-1] -
			    xij[n_zeta-1][n_eta-2]);

    ye_l[n_eta-1] = HALF * (yij[0][n_eta-1] - yij[0][n_eta-2]);
    ye_r[n_eta-1] = HALF * (yij[n_zeta-1][n_eta-1] -
			    yij[n_zeta-1][n_eta-2]);

    xz_l[n_eta-1] = ZERO;
    xz_r[n_eta-1] = ZERO;

    yz_l[n_eta-1] = ZERO;
    yz_r[n_eta-1] = ZERO;

    for ( j = 1; j < n_eta - 1; ++j ) {
      xze_l[j] = HALF * (xz_l[j+1] - xz_l[j-1]);
      xze_r[j] = HALF * (xz_r[j+1] - xz_r[j-1]);

      yze_l[j] = HALF * (yz_l[j+1] - yz_l[j-1]);
      yze_r[j] = HALF * (yz_r[j+1] - yz_r[j-1]);
    } /* endfor */

    /* Begin a new iteration cycle and use the current values of xij
       and yij to determine the derivatives xee and yee at the bottom
       (south) and top (north) boundaries as well as the derivatives
       xzz and yzz at the left (west) and right (east) boundaries, 
       repectively.  These values may then be used to determine pb, 
       pt, pl, pr, qb, qt, ql, and qr.  */

    num_iter = 0;

  next_iteration: ;

    num_iter = num_iter + 1;

    if (i_poisson == 0) goto no_source_term_boundary_evaluation;

    /* BOTTOM AND TOP */

    for ( i = 1; i < n_zeta - 1; ++i ) {
      xee_b[i] = HALF * (-SEVEN*xij[i][0] + EIGHT*xij[i][1] -
			 xij[i][2]) - THREE * xe_b[i];
      xee_t[i] = HALF * (-SEVEN*xij[i][n_eta-1] + 
			 EIGHT*xij[i][n_eta-2] -
			 xij[i][n_eta-3]) + THREE * xe_t[i];

      yee_b[i] = HALF * (-SEVEN*yij[i][0] + EIGHT*yij[i][1] -
			 yij[i][2]) - THREE * ye_b[i];
      yee_t[i] = HALF * (-SEVEN*yij[i][n_eta-1] + 
			 EIGHT*yij[i][n_eta-2] -
			 yij[i][n_eta-3]) + THREE * ye_t[i];

      alpha = xe_b[i] * xe_b[i] + ye_b[i] * ye_b[i];
      beta  = xz_b[i] * xe_b[i] + yz_b[i] * ye_b[i];
      gamma = xz_b[i] * xz_b[i] + yz_b[i] * yz_b[i];
      r1 = - ONE * (alpha * xzz_b[i] - TWO * beta * xze_b[i] +
		    gamma * xee_b[i]) / (jacob_b[i] * jacob_b[i]);
      r2 = - ONE * (alpha * yzz_b[i] - TWO * beta * yze_b[i] +
		    gamma * yee_b[i]) / (jacob_b[i] * jacob_b[i]);

      dpb = (ye_b[i] * r1 - xe_b[i] * r2) / jacob_b[i] - pb[i];
      if (fabs(pb[i]) <= ONE ) {
	dpmax = HALF;
      } else {
	dpmax = HALF * fabs(pb[i]);
      } /* endif */
      if (wsur * fabs(dpb) <= dpmax) {
	pb[i] = pb[i] + wsur * dpb;
      } else {
	if (dpb >= ZERO) {
	  pb[i] = pb[i] + dpmax;
	} else {
	  pb[i] = pb[i] - dpmax;
	} /* endif */
      } /* endif */

      dqb = (xz_b[i] * r2 - yz_b[i] * r1) / jacob_b[i] - qb[i];
      if (fabs(qb[i]) <= ONE ) {
	dqmax = HALF;
      } else {
	dqmax = HALF * fabs(qb[i]);
      } /* endif */
      if (wsur * fabs(dqb) <= dqmax) {
	qb[i] = qb[i] + wsur * dqb;
      } else {
	if (dqb >= ZERO) {
	  qb[i] = qb[i] + dqmax;
	} else {
	  qb[i] = qb[i] - dqmax;
	} /* endif */
      } /* endif */

      alpha = xe_t[i] * xe_t[i] + ye_t[i] * ye_t[i];
      beta  = xz_t[i] * xe_t[i] + yz_t[i] * ye_t[i];
      gamma = xz_t[i] * xz_t[i] + yz_t[i] * yz_t[i];
      r3 = - ONE * (alpha * xzz_t[i] - TWO * beta * xze_t[i] +
		    gamma * xee_t[i]) / (jacob_t[i] * jacob_t[i]);
      r4 = - ONE * (alpha * yzz_t[i] - TWO * beta * yze_t[i] +
		    gamma * yee_t[i]) / (jacob_t[i] * jacob_t[i]);

      dpt = (ye_t[i] * r3 - xe_t[i] * r4) / jacob_t[i] - pt[i];
      if (fabs(pt[i]) <= ONE ) {
	dpmax = HALF;
      } else {
	dpmax = HALF * fabs(pt[i]);
      } /* endif */
      if (wsur * fabs(dpt) <= dpmax) {
	pt[i] = pt[i] + wsur * dpt;
      } else {
	if (dpt >= ZERO) {
	  pt[i] = pt[i] + dpmax;
	} else {
	  pt[i] = pt[i] - dpmax;
	} /* endif */
      } /* endif */

      dqt = (xz_t[i] * r4 - yz_t[i] * r3) / jacob_t[i] - qt[i];
      if (fabs(qt[i]) <= ONE ) {
	dqmax = HALF;
      } else {
	dqmax = HALF * fabs(qt[i]);
      } /* endif */
      if (wsur * fabs(dqt) <= dqmax) {
	qt[i] = qt[i] + wsur * dqt;
      } else {
	if (dqt >= ZERO) {
	  qt[i] = qt[i] + dqmax;
	} else {
	  qt[i] = qt[i] - dqmax;
	} /* endif */
      } /* endif */
    } /* endfor */

    /* LEFT AND RIGHT */

    for ( j = 1; j < n_eta - 1; ++j ) {
      xzz_l[j] = HALF * (-SEVEN*xij[0][j] + EIGHT*xij[1][j] -
			 xij[2][j]) - THREE * xz_l[j];
      xzz_r[j] = HALF * (-SEVEN*xij[n_zeta-1][j] + 
			 EIGHT*xij[n_zeta-2][j] -
			 xij[n_zeta-3][j]) + THREE * xz_r[j];

      yzz_l[j] = HALF * (-SEVEN*yij[0][j] + EIGHT*yij[1][j] -
			 yij[2][j]) - THREE * yz_l[j];
      yzz_r[j] = HALF * (-SEVEN*yij[n_zeta-1][j] + 
			 EIGHT*yij[n_zeta-2][j] -
			 yij[n_zeta-3][j]) + THREE * yz_r[j];

      alpha = xe_l[j] * xe_l[j] + ye_l[j] * ye_l[j];
      beta  = xz_l[j] * xe_l[j] + yz_l[j] * ye_l[j];
      gamma = xz_l[j] * xz_l[j] + yz_l[j] * yz_l[j];
      r1 = - ONE * (alpha * xzz_l[j] - TWO * beta * xze_l[j] +
		    gamma * xee_l[j]) / (jacob_l[j] * jacob_l[j]);
      r2 = - ONE * (alpha * yzz_l[j] - TWO * beta * yze_l[j] +
		    gamma * yee_l[j]) / (jacob_l[j] * jacob_l[j]);

      dpl = (ye_l[j] * r1 - xe_l[j] * r2) / jacob_l[j] - pl[j];
      if (fabs(pl[j]) <= ONE ) {
	dpmax = HALF;
      } else {
	dpmax = HALF * fabs(pl[j]);
      } /* endif */
      if (wsur * fabs(dpl) <= dpmax) {
	pl[j] = pl[j] + wsur * dpl;
      } else {
	if (dpl >= ZERO) {
	  pl[j] = pl[j] + dpmax;
	} else {
	  pl[j] = pl[j] - dpmax;
	} /* endif */
      } /* endif */

      dql = (xz_l[j] * r2 - yz_l[j] * r1) / jacob_l[j] - ql[j];
      if (fabs(ql[j]) <= ONE ) {
	dqmax = HALF;
      } else {
	dqmax = HALF * fabs(ql[j]);
      } /* endif */
      if (wsur * fabs(dql) <= dqmax) {
	ql[j] = ql[j] + wsur * dql;
      } else {
	if (dql >= ZERO) {
	  ql[j] = ql[j] + dqmax;
	} else {
	  ql[j] = ql[j] - dqmax;
	} /* endif */
      } /* endif */

      alpha = xe_r[j] * xe_r[j] + ye_r[j] * ye_r[j];
      beta  = xz_r[j] * xe_r[j] + yz_r[j] * ye_r[j];
      gamma = xz_r[j] * xz_r[j] + yz_r[j] * yz_r[j];
      r3 = - ONE * (alpha * xzz_r[j] - TWO * beta * xze_r[j] +
		    gamma * xee_r[j]) / (jacob_r[j] * jacob_r[j]);
      r4 = - ONE * (alpha * yzz_r[j] - TWO * beta * yze_r[j] +
		    gamma * yee_r[j]) / (jacob_r[j] * jacob_r[j]);

      dpr = (ye_r[j] * r3 - xe_r[j] * r4) / jacob_r[j] - pr[j];
      if (fabs(pr[j]) <= ONE ) {
	dpmax = HALF;
      } else {
	dpmax = HALF * fabs(pr[j]);
      } /* endif */
      if (wsur * fabs(dpr) <= dpmax) {
	pr[j] = pr[j] + wsur * dpr;
      } else {
	if (dpr >= ZERO) {
	  pr[j] = pr[j] + dpmax;
	} else {
	  pr[j] = pr[j] - dpmax;
	} /* endif */
      } /* endif */

      dqr = (xz_r[j] * r4 - yz_r[j] * r3) / jacob_r[j] - qr[j];
      if (fabs(qr[j]) <= ONE ) {
	dqmax = HALF;
      } else {
	dqmax = HALF * fabs(qr[j]);
      } /* endif */
      if (wsur * fabs(dqr) <= dqmax) {
	qr[j] = qr[j] + wsur * dqr;
      } else {
	if (dqr >= ZERO) {
	  qr[j] = qr[j] + dqmax;
	} else {
	  qr[j] = qr[j] - dqmax;
	} /* endif */
      } /* endif */
    } /* endfor */

    /* Perform one step in the SOR Gauss-Seidel centred 
       finite-difference solution technique. */

  no_source_term_boundary_evaluation: ;

    for ( j = 1; j < n_eta - 1; ++j ) {
      if (i_poisson != 0) {
	fb = exp ( - a * double(j) );
	ft = exp ( b * (double(j) - double(n_eta) + ONE) );
      } /* endfor */

      for ( i = 1; i < n_zeta - 1; ++i ) { 
	if (i_poisson != 0) {
	  fl = exp ( - c * double(i) );
	  fr = exp ( d * (double(i) - double(n_zeta) + ONE) );
 
	  fbl = ONE -
	    exp ( -e * sqrt( double(i)*double(i) +
			     double(j)*double(j) ) );
	  fbr = ONE -
	    exp ( -e * sqrt( (double(i) - double(n_zeta) + ONE) *
			     (double(i) - double(n_zeta) + ONE) +
			     double(j)*double(j) ) );
	  ftl = ONE -
	    exp ( -e * sqrt( double(i)*double(i) +
			     (double(j) - double(n_eta) + ONE) *
			     (double(j) - double(n_eta) + ONE) ) );
	  ftr = ONE -
	    exp ( -e * sqrt( (double(i) - double(n_zeta) + ONE) *
			     (double(i) - double(n_zeta) + ONE) +
			     (double(j) - double(n_eta) + ONE) *
			     (double(j) - double(n_eta) + ONE) ) );
	  pij = (pb[i] * fb + pt[i] * ft + pl[j] * fl + pr[j] * fr) *
	    (fbl * fbr * ftl * ftr);
	  qij = (qb[i] * fb + qt[i] * ft + ql[j] * fl + qr[j] * fr) *
	    (fbl * fbr * ftl * ftr);
	} else {
	  pij = ZERO;
	  qij = ZERO;
	} /* endif */

	xzij = HALF * (xij[i+1][j] - xij[i-1][j]);
	yzij = HALF * (yij[i+1][j] - yij[i-1][j]);
	xeij = HALF * (xij[i][j+1] - xij[i][j-1]);
	yeij = HALF * (yij[i][j+1] - yij[i][j-1]);

	jacobij = xzij * yeij - xeij * yzij;

	alpha = xeij * xeij + yeij * yeij;
	beta  = xzij * xeij + yzij * yeij;
	gamma = xzij * xzij + yzij * yzij;

	dxij = HALF * (jacobij * jacobij * (xzij * pij +
					    xeij * qij) + alpha * (xij[i+1][j] + xij[i-1][j]) +
		       gamma * (xij[i][j+1] + xij[i][j-1]) - HALF * beta *
		       (xij[i+1][j+1] - xij[i+1][j-1] - xij[i-1][j+1] +
			xij[i-1][j-1])) / (alpha + gamma) - xij[i][j];
	xij[i][j] = xij[i][j] + wsor * dxij;
          
	dyij = HALF * (jacobij * jacobij * (yzij * pij +
					    yeij * qij) + alpha * (yij[i+1][j] + yij[i-1][j]) +
		       gamma * (yij[i][j+1] + yij[i][j-1]) - HALF * beta *
		       (yij[i+1][j+1] - yij[i+1][j-1] - yij[i-1][j+1] +
			yij[i-1][j-1])) / (alpha + gamma) - yij[i][j];
	yij[i][j] = yij[i][j] + wsor * dyij;
      } /* endfor */
    } /* endfor */

    /* Check to see if the grid block smoothing is complete.
       If not go to next_iteration. */

    if (num_iter < Number_of_Iterations) goto next_iteration;

    /* Save the newly computed interior node locations for 
       the quadrilateral grid block. */

    for (j  = JNl ; j <= JNu ; ++j ) {
      for ( i = INl ; i <= INu ; ++i ) {
	ii = i - INl;
	jj = j - JNl;
	Node[i][j].X.x = xij[ii][jj];
	Node[i][j].X.y = yij[ii][jj];
      } /* endfor */
    }/* endfor */

    /* Require update of the interior cells geometric properties. */
    Schedule_Interior_Mesh_Update();

    /* Re-compute the exterior nodes for the quadrilateral mesh block. */

    Update_Exterior_Nodes();

    /* Re-compute the cell values for the quadrilateral mesh block. */

    Update_Cells();

    /* Delete (deallocate) the Poisson equation solution variables. */

    for ( i = 0; i <= n_zeta-1 ; ++i ) {
      delete []xij[i]; xij[i] = NULL;
      delete []yij[i]; yij[i] = NULL;
    } /* endfor */
    delete []xij; xij = NULL;
    delete []yij; yij = NULL;

    delete []deltas_b; deltas_b = NULL;
    delete []deltas_t; deltas_t = NULL;
    delete []deltas_l; deltas_l = NULL; 
    delete []deltas_r; deltas_r = NULL;

    delete []xz_b; xz_b = NULL; 
    delete []yz_b; yz_b = NULL;
    delete []xe_b; xe_b = NULL;
    delete []ye_b; ye_b = NULL;
    delete []xzz_b; xzz_b = NULL;
    delete []yzz_b; yzz_b = NULL;
    delete []xee_b; xee_b = NULL;
    delete []yee_b; yee_b = NULL;
    delete []xze_b; xze_b = NULL;
    delete []yze_b; yze_b = NULL;

    delete []xz_t; xz_t = NULL;
    delete []yz_t; yz_t = NULL;
    delete []xe_t; xe_t = NULL; 
    delete []ye_t; ye_t = NULL;
    delete []xzz_t; xzz_t = NULL;
    delete []yzz_t; yzz_t = NULL;
    delete []xee_t; xee_t = NULL;
    delete []yee_t; yee_t = NULL;
    delete []xze_t; xze_t = NULL;
    delete []yze_t; yze_t = NULL;

    delete []xz_l; xz_l = NULL;
    delete []yz_l; yz_l = NULL;
    delete []xe_l; xe_l = NULL;
    delete []ye_l; ye_l = NULL;
    delete []xzz_l; xzz_l = NULL;
    delete []yzz_l; yzz_l = NULL;
    delete []xee_l; xee_l = NULL;
    delete []yee_l; yee_l = NULL;
    delete []xze_l; xze_l = NULL;
    delete []yze_l; yze_l = NULL;

    delete []xz_r; xz_r = NULL;
    delete []yz_r; yz_r = NULL;
    delete []xe_r; xe_r = NULL;
    delete []ye_r; ye_r = NULL;
    delete []xzz_r; xzz_r = NULL;
    delete []yzz_r; yzz_r = NULL;
    delete []xee_r; xee_r = NULL;
    delete []yee_r; yee_r = NULL;
    delete []xze_r; xze_r = NULL;
    delete []yze_r; yze_r = NULL;

    delete []pb; pb = NULL;
    delete []pt; pt = NULL;
    delete []pl; pl = NULL;
    delete []pr; pr = NULL;

    delete []qb; qb = NULL;
    delete []qt; qt = NULL;
    delete []ql; ql = NULL;
    delete []qr; qr = NULL;

    delete []jacob_b; jacob_b = NULL;
    delete []jacob_t; jacob_t = NULL;
    delete []jacob_l; jacob_l = NULL;
    delete []jacob_r; jacob_r = NULL;

  }
}


/**********************************************************************
 * Routine: Smooth_Rocket_Motor                                       *
 **********************************************************************/
void Grid2D_Quad_Block_HO::Smooth_Rocket_Motor(const double &Length_Chamber,
					       const double &Radius_Chamber,
					       const double &Length_Chamber_To_Throat,
					       const double &Length_Nozzle,
					       const double &Radius_Nozzle_Exit,
					       const double &Radius_Nozzle_Throat,
					       const double &Radius_Grain,
					       const int &Nozzle_Type,
					       const double &Stretching_Factor_Idir,
					       const double &Stretching_Factor_Jdir,
					       const int &sector,
					       const int &level,
					       const int &di,const int &dj,
					       const int &ri,const int &rj,
					       const int &NRi,const int &NRj) {
  
  // Exit immediately if the current block is in the chamber.
  if (Node[INl+5][JNl+5].X.x < ZERO) return ;
  
  Smooth_Quad_Block(min(250,4*max(NCi-4,NCi-4)));
  
  return ;
}

/*!
 * Set the boundary condition type for the north,       
 * south, east, and west boundaries of the              
 * quadrilateral mesh block.                            
 *                                                      
 */

void Grid2D_Quad_Block_HO::Set_BCs(void) {

    int i, j, bc_type_left, bc_type_right;
    double S_i, S_j, 
           s_north, s_south, s_east, s_west,
           smax_north, smax_south, smax_east, smax_west;

    if (BndNorthSpline.np == 0) {
       for ( i = ICl-Nghost; i <= ICu+Nghost; ++i) {
	 BCtypeN[i] = BC_NONE;
       } /* endfor */
    } else {
       BCtypeN[ICl-1] = BCtype(SminN,
			       BndNorthSpline);
       BCtypeN[ICu+1] = BCtype(SmaxN,
			       BndNorthSpline);
       for (int GCell=2; GCell<= Nghost; ++GCell){
	 BCtypeN[ICl-GCell] = BCtypeN[ICl-GCell+1];
	 BCtypeN[ICu+GCell] = BCtypeN[ICu+GCell-1];
       }
       for ( i = INl ; i <= INu-1 ; ++i) {
	 s_north = getS(Node[i][JNu].X, BndNorthSpline);
	 bc_type_left = BCtype(s_north, BndNorthSpline);
	 s_north = getS(Node[i+1][JNu].X, BndNorthSpline);
	 bc_type_right = BCtype(s_north, BndNorthSpline);

	 if (bc_type_left == bc_type_right) {
	   BCtypeN[i] = bc_type_left;
	 } /* endif */
       } /* endfor */
    } /* endif */

    if (ExtendWest_BndNorthSpline.np != 0){
      for (i = ICl-1 ; i>= ICl-Nghost; --i){
	s_north = getS(Node[i+1][JNu].X, ExtendWest_BndNorthSpline);
	bc_type_left = BCtype(s_north, ExtendWest_BndNorthSpline);
	s_north = getS(Node[i  ][JNu].X, ExtendWest_BndNorthSpline);
	bc_type_right = BCtype(s_north,  ExtendWest_BndNorthSpline);

	if (bc_type_left == bc_type_right) {
	  BCtypeN[i] = bc_type_left;
	} else {
	  BCtypeN[i] = bc_type_right;
	} /* endif */
      }
    }

    if (ExtendEast_BndNorthSpline.np != 0){
      for (i = ICu+1; i<= ICu+Nghost; ++i){
	s_north = getS(Node[i ][JNu].X, ExtendEast_BndNorthSpline);
	bc_type_left = BCtype(s_north, ExtendEast_BndNorthSpline);
	s_north = getS(Node[i+1][JNu].X, ExtendEast_BndNorthSpline);
	bc_type_right = BCtype(s_north,  ExtendEast_BndNorthSpline);

	if (bc_type_left == bc_type_right) {
	  BCtypeN[i] = bc_type_left;
	} else {
	  BCtypeN[i] = bc_type_right;
	} /* endif */
      }
    }

    if (BndSouthSpline.np == 0) {
       for ( i = ICl-Nghost; i <= ICu+Nghost ; ++i) {
	   BCtypeS[i] = BC_NONE;
       } /* endfor */
    } else {
       BCtypeS[ICl-1] = BCtype(SminS, 
			       BndSouthSpline);
       BCtypeS[ICu+1] = BCtype(SmaxS, 
			       BndSouthSpline);
       for (int GCell=2; GCell<= Nghost; ++GCell){
	 BCtypeS[ICl-GCell] = BCtypeS[ICl-GCell+1];
	 BCtypeS[ICu+GCell] = BCtypeS[ICu+GCell-1];
       }
       for ( i = INl ; i <= INu-1 ; ++i) {
	 s_south = getS(Node[i][JNl].X, BndSouthSpline);
	 bc_type_left = BCtype(s_south, BndSouthSpline);
	 s_south = getS(Node[i+1][JNl].X, BndSouthSpline);
	 bc_type_right = BCtype(s_south, BndSouthSpline);

	 if (bc_type_left == bc_type_right) {
	   BCtypeS[i] = bc_type_left;
	 } /* endif */
       } /* endfor */
    } /* endif */

    if (ExtendWest_BndSouthSpline.np != 0){
      for (i = ICl-1 ; i>= ICl-Nghost; --i){
	s_south = getS(Node[i+1][JNl].X, ExtendWest_BndSouthSpline);
	bc_type_left = BCtype(s_south, ExtendWest_BndSouthSpline);
	s_south = getS(Node[i  ][JNl].X, ExtendWest_BndSouthSpline);
	bc_type_right = BCtype(s_south,  ExtendWest_BndSouthSpline);
	
	if (bc_type_left == bc_type_right) {
	  BCtypeS[i] = bc_type_left;
	} else {
	  BCtypeS[i] = bc_type_right;
	} /* endif */	 
      }
    }

    if (ExtendEast_BndSouthSpline.np != 0){
      for (i = ICu+1; i<= ICu+Nghost; ++i){
	s_south = getS(Node[i ][JNl].X, ExtendEast_BndSouthSpline);
	bc_type_left = BCtype(s_south, ExtendEast_BndSouthSpline);
	s_south = getS(Node[i+1][JNl].X, ExtendEast_BndSouthSpline);
	bc_type_right = BCtype(s_south,  ExtendEast_BndSouthSpline);

	if (bc_type_left == bc_type_right) {
	  BCtypeS[i] = bc_type_left;
	} else {
	  BCtypeS[i] = bc_type_right;
	} /* endif */
      }
    }

    if (BndEastSpline.np == 0) {
       for ( j = JCl-Nghost; j <= JCu+Nghost ; ++j) {
	   BCtypeE[j] = BC_NONE;
       } /* endfor */
    } else {
       BCtypeE[JCl-1] = BCtype(SminE,
			       BndEastSpline);
       BCtypeE[JCu+1] = BCtype(SmaxE, 
			       BndEastSpline);
       for (int GCell=2; GCell<= Nghost; ++GCell){
	 BCtypeE[JCl-GCell] = BCtypeE[JCl-GCell+1];
	 BCtypeE[JCu+GCell] = BCtypeE[JCu+GCell-1];
       }
       for ( j = JNl ; j <= JNu-1 ; ++j) {
	 s_east = getS(Node[INu][j].X, BndEastSpline);
	 bc_type_left = BCtype(s_east, BndEastSpline);
	 s_east = getS(Node[INu][j+1].X, BndEastSpline);
	 bc_type_right = BCtype(s_east, BndEastSpline);

	 if (bc_type_left == bc_type_right) {
	   BCtypeE[j] = bc_type_left;
	 } /* endif */
       } /* endfor */
    } /* endif */

    if (ExtendSouth_BndEastSpline.np != 0){
      for ( j = JCl-1; j >= JCl-Nghost ; --j) {
	s_east = getS(Node[INu][j+1].X, ExtendSouth_BndEastSpline);
	bc_type_left = BCtype(s_east, ExtendSouth_BndEastSpline);
	s_east = getS(Node[INu][j].X, ExtendSouth_BndEastSpline);
	bc_type_right = BCtype(s_east, ExtendSouth_BndEastSpline);

	if (bc_type_left == bc_type_right) {
	  BCtypeE[j] = bc_type_left;
	} else {
	  BCtypeE[j] = bc_type_right;
	} /* endif */
      } /* endfor */      
    }

    if (ExtendNorth_BndEastSpline.np != 0){
      for ( j = JCu+1; j <= JCu+Nghost ; ++j) {
	s_east = getS(Node[INu][j  ].X, ExtendNorth_BndEastSpline);
	bc_type_left = BCtype(s_east, ExtendNorth_BndEastSpline);
	s_east = getS(Node[INu][j+1].X, ExtendNorth_BndEastSpline);
	bc_type_right = BCtype(s_east, ExtendNorth_BndEastSpline);

	if (bc_type_left == bc_type_right) {
	  BCtypeE[j] = bc_type_left;
	} else {
	  BCtypeE[j] = bc_type_right;
	} /* endif */
      } /* endfor */      
    }
    
    if (BndWestSpline.np == 0) {
       for ( j = JCl-Nghost ; j <= JCu+Nghost ; ++j) {
	   BCtypeW[j] = BC_NONE;
       } /* endfor */
    } else {
       BCtypeW[JCl-1] = BCtype(SminW, 
			       BndWestSpline);
       BCtypeW[JCu+1] = BCtype(SmaxW, 
			       BndWestSpline);
       for (int GCell=2; GCell<= Nghost; ++GCell){
	 BCtypeW[JCl-GCell] = BCtypeW[JCl-GCell+1];
	 BCtypeW[JCu+GCell] = BCtypeW[JCu+GCell-1];
       }
       for ( j = JNl ; j <= JNu-1 ; ++j) {
	 s_west = getS(Node[INl][j].X, BndWestSpline);
	 bc_type_left = BCtype(s_west, BndWestSpline);
	 s_west = getS(Node[INl][j+1].X, BndWestSpline);
	 bc_type_right = BCtype(s_west, BndWestSpline);

	 if (bc_type_left == bc_type_right) {
	   BCtypeW[j] = bc_type_left;
	 } /* endif */
       } /* endfor */
    } /* endif */

    if (ExtendSouth_BndWestSpline.np != 0){
      for ( j = JCl-1; j >= JCl-Nghost ; --j) {
	s_west = getS(Node[INl][j+1].X, ExtendSouth_BndWestSpline);
	bc_type_left = BCtype(s_west, ExtendSouth_BndWestSpline);
	s_west = getS(Node[INl][j].X, ExtendSouth_BndWestSpline);
	bc_type_right = BCtype(s_west, ExtendSouth_BndWestSpline);

	if (bc_type_left == bc_type_right) {
	  BCtypeW[j] = bc_type_left;
	} else {
	  BCtypeW[j] = bc_type_right;
	} /* endif */
      } /* endfor */      
    }

    if (ExtendNorth_BndWestSpline.np != 0){
      for ( j = JCu+1; j <= JCu+Nghost ; ++j) {
	s_west = getS(Node[INl][j  ].X, ExtendNorth_BndWestSpline);
	bc_type_left = BCtype(s_west, ExtendNorth_BndWestSpline);
	s_west = getS(Node[INl][j+1].X, ExtendNorth_BndWestSpline);
	bc_type_right = BCtype(s_west, ExtendNorth_BndWestSpline);

	if (bc_type_left == bc_type_right) {
	  BCtypeW[j] = bc_type_left;
	} else {
	  BCtypeW[j] = bc_type_right;
	} /* endif */
      } /* endfor */      
    }    

}

/*!
 * Updates the exterior nodes for the quadrilateral     
 * mesh block.
 *
 * \note This subroutine DOESN'T update the cell information.
 */
void Grid2D_Quad_Block_HO::Update_Exterior_Nodes(void) {

    int i, j;
    Vector2D norm_dir, X_norm, X_tan;

    // Update West and East boundary nodes.
    for ( j = JNl ; j <= JNu ; ++j) {
      if (BCtypeW[j] == BCtypeW[j-1] &&
	  BCtypeW[j] != BC_REFLECTION &&
	  BCtypeW[j] != BC_PERIODIC &&
	  BCtypeW[j] != BC_NO_SLIP &&
	  BCtypeW[j] != BC_WALL_VISCOUS_ISOTHERMAL &&
	  BCtypeW[j] != BC_WALL_VISCOUS_HEATFLUX &&
	  BCtypeW[j] != BC_MOVING_WALL &&
	  BCtypeW[j] != BC_MOVING_WALL_ISOTHERMAL &&
	  BCtypeW[j] != BC_MOVING_WALL_HEATFLUX &&
	  BCtypeW[j] != BC_BURNING_SURFACE &&
	  BCtypeW[j] != BC_MASS_INJECTION &&
	  BCtypeW[j] != BC_WALL_INVISCID) {
	for(int GCell=1; GCell<=Nghost; ++GCell){
	  Node[INl-GCell][j].X = Node[INl][j].X -
	    (Node[INl+GCell][j].X - 
	     Node[INl][j].X);
	}
      } else if (BCtypeW[j] == BCtypeW[j-1] &&
		 (BCtypeW[j] == BC_REFLECTION ||
		  BCtypeW[j] == BC_NO_SLIP ||
		  BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
		  BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX ||
		  BCtypeW[j] == BC_MOVING_WALL ||
		  BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL ||
		  BCtypeW[j] == BC_MOVING_WALL_HEATFLUX ||
		  BCtypeW[j] == BC_BURNING_SURFACE ||
		  BCtypeW[j] == BC_MASS_INJECTION || 
		  BCtypeW[j] == BC_WALL_INVISCID)) {
	if (j > JNl && j < JNu) {
	  norm_dir = - HALF*(nfaceW(ICl, j) + 
			     nfaceW(ICl, j-1));
	  for(int GCell=1; GCell<=Nghost; ++GCell){
	    X_norm = ((Node[INl+GCell][j].X - 
		       Node[INl][j].X) * norm_dir) * norm_dir;
	    X_tan = (Node[INl+GCell][j].X - 
		     Node[INl][j].X) - X_norm;
	    Node[INl-GCell][j].X = Node[INl][j].X -
	      X_norm + X_tan;
	  }
	} else if ( j == JNl ) {
	  //norm_dir = - nfaceW(ICl, j);
	  for(int GCell=1; GCell<=Nghost; ++GCell){
	    Node[INl-GCell][j].X = Node[INl][j].X -
	      (Node[INl+GCell][j].X - 
	       Node[INl][j].X);
	  }
	} else {
	  //norm_dir = - nfaceW(ICl, j-1);
	  for(int GCell=1; GCell<=Nghost; ++GCell){
	    Node[INl-GCell][j].X = Node[INl][j].X -
	      (Node[INl+GCell][j].X - 
	       Node[INl][j].X);
	  }
	} /* endif */
      } else if (BCtypeW[j] == BCtypeW[j-1] &&
		 BCtypeW[j] == BC_PERIODIC) {
	for(int GCell=1; GCell<=Nghost; ++GCell){
	  Node[INl-GCell][j].X = Node[INl][j].X -
	    (Node[INu][j].X - 
	     Node[INu-GCell][j].X);
	}
      } /* endif */

      if (BCtypeE[j] == BCtypeE[j-1] &&
	  BCtypeE[j] != BC_REFLECTION &&
	  BCtypeE[j] != BC_PERIODIC &&
	  BCtypeE[j] != BC_NO_SLIP &&
	  BCtypeE[j] != BC_WALL_VISCOUS_ISOTHERMAL &&
	  BCtypeE[j] != BC_WALL_VISCOUS_HEATFLUX &&
	  BCtypeE[j] != BC_MOVING_WALL &&
	  BCtypeE[j] != BC_MOVING_WALL_ISOTHERMAL &&
	  BCtypeE[j] != BC_MOVING_WALL_HEATFLUX &&
	  BCtypeE[j] != BC_BURNING_SURFACE &&
	  BCtypeE[j] != BC_MASS_INJECTION &&
	  BCtypeE[j] != BC_WALL_INVISCID) {
	for(int GCell=1; GCell<=Nghost; ++GCell){
	  Node[INu+GCell][j].X = ( Node[INu][j].X +
				   (Node[INu][j].X - 
				    Node[INu-GCell][j].X) );
	}
      } else if (BCtypeE[j] == BCtypeE[j-1] &&
		 (BCtypeE[j] == BC_REFLECTION ||	  
		  BCtypeE[j] == BC_NO_SLIP ||
		  BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
		  BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX ||
		  BCtypeE[j] == BC_MOVING_WALL ||
		  BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL ||
		  BCtypeE[j] == BC_MOVING_WALL_HEATFLUX ||
		  BCtypeE[j] == BC_BURNING_SURFACE ||
		  BCtypeE[j] == BC_MASS_INJECTION ||
		  BCtypeE[j] == BC_WALL_INVISCID)) {
	if (j > JNl && j < JNu) {
	  norm_dir = HALF*(nfaceE(ICu, j) + 
			   nfaceE(ICu, j-1));
	  for(int GCell=1; GCell<=Nghost; ++GCell){
	    X_norm = ((Node[INu][j].X - 
		       Node[INu-GCell][j].X) * norm_dir) * norm_dir;
	    X_tan = (Node[INu][j].X - 
		     Node[INu-GCell][j].X) - X_norm;
	    Node[INu+GCell][j].X = Node[INu][j].X +
	      X_norm - X_tan;
	  }
	} else if ( j == JNl ) {
	  //norm_dir = nfaceE(ICu, j);
	  for(int GCell=1; GCell<=Nghost; ++GCell){
	    Node[INu+GCell][j].X = Node[INu][j].X +
	      (Node[INu][j].X - 
	       Node[INu-GCell][j].X);
	  }
	} else {
	  //norm_dir = nfaceE(ICu, j-1);
	  for(int GCell=1; GCell<=Nghost; ++GCell){
	    Node[INu+GCell][j].X = ( Node[INu][j].X +
				     (Node[INu][j].X - 
				      Node[INu-GCell][j].X) );
	  }
	} /* endif */
      } else if (BCtypeE[j] == BCtypeE[j-1] &&
		 BCtypeE[j] == BC_PERIODIC) {
	for(int GCell=1; GCell<=Nghost; ++GCell){
	  Node[INu+GCell][j].X = ( Node[INu][j].X +
				   (Node[INl+GCell][j].X - 
				    Node[INl][j].X) );
	}
      } /* endif */
    } /* endfor */
    
    // Update South and North boundary nodes.
    for ( i = INl ; i <= INu ; ++i) {
	  if (BCtypeS[i] == BCtypeS[i-1] &&
	      BCtypeS[i] != BC_REFLECTION &&
	      BCtypeS[i] != BC_PERIODIC &&
	      BCtypeS[i] != BC_NO_SLIP &&
	      BCtypeS[i] != BC_WALL_VISCOUS_ISOTHERMAL &&
	      BCtypeS[i] != BC_WALL_VISCOUS_HEATFLUX &&
	      BCtypeS[i] != BC_MOVING_WALL &&
	      BCtypeS[i] != BC_MOVING_WALL_ISOTHERMAL &&
	      BCtypeS[i] != BC_MOVING_WALL_HEATFLUX &&
	      BCtypeS[i] != BC_BURNING_SURFACE &&
	      BCtypeS[i] != BC_MASS_INJECTION &&
	      BCtypeS[i] != BC_RINGLEB_FLOW &&
	      ( BCtypeS[i] != BC_WALL_INVISCID || Mesh_Requiring_NonReflected_South_Ghost_Cells == ON )) {
	    for(int GCell=1; GCell<=Nghost; ++GCell){
	      Node[i][JNl-GCell].X = ( Node[i][JNl].X -
				       (Node[i][JNl+GCell].X - 
					Node[i][JNl].X) );
	    }
	  } else if (BCtypeS[i] == BCtypeS[i-1] &&
		     (BCtypeS[i] == BC_REFLECTION ||  
		      BCtypeS[i] == BC_NO_SLIP ||
		      BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
		      BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX ||
		      BCtypeS[i] == BC_MOVING_WALL ||
		      BCtypeS[i] == BC_MOVING_WALL_ISOTHERMAL ||
		      BCtypeS[i] == BC_MOVING_WALL_HEATFLUX ||
		      BCtypeS[i] == BC_BURNING_SURFACE ||
		      BCtypeS[i] == BC_MASS_INJECTION ||
		      BCtypeS[i] == BC_RINGLEB_FLOW ||
		      BCtypeS[i] == BC_WALL_INVISCID)) {
	    if (i > INl && i < INu) {
	      if (lfaceS(i,JCl) > NANO &&
		  lfaceS(i-1,JCl) > NANO) {
		norm_dir = - HALF*(nfaceS(i, JCl) + 
				   nfaceS(i-1, JCl));
	      } else if (lfaceS(i,JCl) > NANO) {
		norm_dir = - nfaceS(i, JCl);
	      } else if (lfaceS(i,JCl) > NANO) {
		norm_dir = - nfaceS(i-1, JCl);
	      } else {
	      }
	      for(int GCell=1; GCell<=Nghost; ++GCell){
		X_norm = ((Node[i][JNl+GCell].X - 
			   Node[i][JNl].X) * norm_dir) * norm_dir;
		X_tan = (Node[i][JNl+GCell].X - 
			 Node[i][JNl].X) - X_norm;
		Node[i][JNl-GCell].X = ( Node[i][JNl].X -
					 X_norm + X_tan );
	      }
	    } else if ( i == INl ) {
	      //norm_dir = - nfaceS(i, JCl);
	      for(int GCell=1; GCell<=Nghost; ++GCell){
		Node[i][JNl-GCell].X = ( Node[i][JNl].X -
					 (Node[i][JNl+GCell].X - 
					  Node[i][JNl].X) );
	      }
	    } else {
	      //norm_dir = - nfaceS(i-1, JCl);
	      for(int GCell=1; GCell<=Nghost; ++GCell){
		Node[i][JNl-GCell].X = ( Node[i][JNl].X -
					 (Node[i][JNl+GCell].X - 
					  Node[i][JNl].X) );
	      }
	    } /* endif */
	  } else if (BCtypeS[i] == BCtypeS[i-1] &&
		     BCtypeS[i] == BC_PERIODIC) {
	    for(int GCell=1; GCell<=Nghost; ++GCell){
	      Node[i][JNl-GCell].X = ( Node[i][JNl].X -
				       (Node[i][JNu].X - 
					Node[i][JNu-GCell].X) );
	  }
        } /* endif */

	  if (BCtypeN[i] == BCtypeN[i-1] &&
	      BCtypeN[i] != BC_REFLECTION &&
	      BCtypeN[i] != BC_PERIODIC &&
	      BCtypeN[i] != BC_NO_SLIP &&
	      BCtypeN[i] != BC_WALL_VISCOUS_ISOTHERMAL &&
	      BCtypeN[i] != BC_WALL_VISCOUS_HEATFLUX &&
	      BCtypeN[i] != BC_MOVING_WALL &&
	      BCtypeN[i] != BC_MOVING_WALL_ISOTHERMAL &&
	      BCtypeN[i] != BC_MOVING_WALL_HEATFLUX  &&
	      BCtypeN[i] != BC_BURNING_SURFACE &&
	      BCtypeN[i] != BC_MASS_INJECTION &&
	      BCtypeN[i] != BC_WALL_INVISCID) {
	    for(int GCell=1; GCell<=Nghost; ++GCell){
	      Node[i][JNu+GCell].X = ( Node[i][JNu].X +
				       (Node[i][JNu].X - 
					Node[i][JNu-GCell].X) );
	    }
	  } else if (BCtypeN[i] == BCtypeN[i-1] &&
		     (BCtypeN[i] == BC_REFLECTION || 
		      BCtypeN[i] == BC_NO_SLIP ||
		      BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
		      BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX ||
		      BCtypeN[i] == BC_MOVING_WALL ||
		      BCtypeN[i] == BC_MOVING_WALL_ISOTHERMAL ||
		      BCtypeN[i] == BC_MOVING_WALL_HEATFLUX ||
		      BCtypeN[i] == BC_BURNING_SURFACE ||
		      BCtypeN[i] == BC_MASS_INJECTION ||
		      BCtypeN[i] == BC_WALL_INVISCID)) {
	    if (i > INl && i < INu) {
 	      norm_dir = HALF*(nfaceN(i, JCu) + 
                               nfaceN(i-1, JCu));
	      for(int GCell=1; GCell<=Nghost; ++GCell){
		X_norm = ((Node[i][JNu].X - 
			   Node[i][JNu-GCell].X) * norm_dir) * norm_dir;
		X_tan = (Node[i][JNu].X - 
			 Node[i][JNu-GCell].X) - X_norm;
		Node[i][JNu+GCell].X = ( Node[i][JNu].X +
					 X_norm - X_tan );
	      }
	    } else if ( i == INl ) {
	      //norm_dir = nfaceN(i, JCu);
	      for(int GCell=1; GCell<=Nghost; ++GCell){
		Node[i][JNu+GCell].X = ( Node[i][JNu].X +
					 (Node[i][JNu].X - 
					  Node[i][JNu-GCell].X) );
	      }
	    } else {
	      //norm_dir = nfaceN(i-1, JCu);
	      for(int GCell=1; GCell<=Nghost; ++GCell){
		Node[i][JNu+GCell].X = ( Node[i][JNu].X +
					 (Node[i][JNu].X - 
					  Node[i][JNu-GCell].X) );
	      }
	    } /* endif */
	  } else if (BCtypeN[i] == BCtypeN[i-1] &&
		     BCtypeN[i] == BC_PERIODIC) {
	    for(int GCell=1; GCell<=Nghost; ++GCell){
	      Node[i][JNu+GCell].X = ( Node[i][JNu].X +
				       (Node[i][JNl+GCell].X - 
					Node[i][JNl].X) );
	    }
	  } /* endif */
    } /* endfor */
    
    // Update SW and SE corner nodes.
    for (int j = JNl-Nghost ; j <= JNl ; ++j) {
      if (BCtypeW[j] != BC_REFLECTION &&
	  BCtypeW[j] != BC_PERIODIC &&  
	  BCtypeW[j] != BC_NO_SLIP &&
	  BCtypeW[j] != BC_WALL_VISCOUS_ISOTHERMAL &&
	  BCtypeW[j] != BC_WALL_VISCOUS_HEATFLUX &&
	  BCtypeW[j] != BC_MOVING_WALL &&
	  BCtypeW[j] != BC_MOVING_WALL_ISOTHERMAL &&
	  BCtypeW[j] != BC_MOVING_WALL_HEATFLUX &&
	  BCtypeW[j] != BC_BURNING_SURFACE &&
	  BCtypeW[j] != BC_MASS_INJECTION &&
	  BCtypeW[j] != BC_WALL_INVISCID) {
	for(int GCell=1; GCell<=Nghost; ++GCell){
	  Node[INl-GCell][j].X = Node[INl][j].X -
	    (Node[INl+GCell][j].X - 
	     Node[INl][j].X);
	}
      } else if (BCtypeW[j] == BC_REFLECTION || 
		 BCtypeW[j] == BC_NO_SLIP ||
		 BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
		 BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX ||
		 BCtypeW[j] == BC_MOVING_WALL ||
		 BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL ||
		 BCtypeW[j] == BC_MOVING_WALL_HEATFLUX ||
		 BCtypeW[j] == BC_BURNING_SURFACE ||
		 BCtypeW[j] == BC_MASS_INJECTION ||
		 BCtypeW[j] == BC_WALL_INVISCID) {
	if (j != JNl-Nghost) {
	  norm_dir = - HALF*(nfaceW(ICl, j) + 
			     nfaceW(ICl, j-1));
	} else {
	  norm_dir = - nfaceW(ICl, j); 
	} /* endif */
	for(int GCell=1; GCell<=Nghost; ++GCell){
	  X_norm = ((Node[INl+GCell][j].X - 
		     Node[INl][j].X) * norm_dir) * norm_dir;
	  X_tan = (Node[INl+GCell][j].X - 
		   Node[INl][j].X) - X_norm;
	  Node[INl-GCell][j].X = ( Node[INl][j].X -
				   X_norm + X_tan );
	}
      } else if (BCtypeW[j] == BC_PERIODIC) {
	for(int GCell=1; GCell<=Nghost; ++GCell){
	  Node[INl-GCell][j].X = ( Node[INl][j].X -
				   (Node[INu][j].X - 
				    Node[INu-GCell][j].X) );
	}
      } /* endif */
      
      if (BCtypeE[j] != BC_REFLECTION &&
	  BCtypeE[j] != BC_PERIODIC &&   
	  BCtypeE[j] != BC_NO_SLIP &&
	  BCtypeE[j] != BC_WALL_VISCOUS_ISOTHERMAL &&
	  BCtypeE[j] != BC_WALL_VISCOUS_HEATFLUX &&
	  BCtypeE[j] != BC_MOVING_WALL &&
	  BCtypeE[j] != BC_MOVING_WALL_ISOTHERMAL &&
	  BCtypeE[j] != BC_MOVING_WALL_HEATFLUX &&
	  BCtypeE[j] != BC_BURNING_SURFACE &&
	  BCtypeE[j] != BC_MASS_INJECTION &&
	  BCtypeE[j] != BC_WALL_INVISCID) {
	for(int GCell=1; GCell<=Nghost; ++GCell){
	  Node[INu+GCell][j].X = ( Node[INu][j].X +
				   (Node[INu][j].X - 
				    Node[INu-GCell][j].X) );
	}
      } else if (BCtypeE[j] == BC_REFLECTION ||
		 BCtypeE[j] == BC_NO_SLIP ||
		 BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
		 BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX ||
		 BCtypeE[j] == BC_MOVING_WALL ||
		 BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL ||
		 BCtypeE[j] == BC_MOVING_WALL_HEATFLUX  ||
		 BCtypeE[j] == BC_BURNING_SURFACE ||
		 BCtypeE[j] == BC_MASS_INJECTION ||
		 BCtypeE[j] == BC_WALL_INVISCID) {
	if (j != JNl-Nghost) {
	  norm_dir = HALF*(nfaceE(ICu, j) + 
			   nfaceE(ICu, j-1));
	} else {
	  norm_dir = nfaceE(ICu, j); 
	} /* endif */
	for(int GCell=1; GCell<=Nghost; ++GCell){
	  X_norm = ((Node[INu][j].X - 
		     Node[INu-GCell][j].X) * norm_dir) * norm_dir;
	  X_tan = (Node[INu][j].X - 
		   Node[INu-GCell][j].X) - X_norm;
	  Node[INu+GCell][j].X = ( Node[INu][j].X +
				   X_norm - X_tan );
	}
      } else if (BCtypeE[j] == BC_PERIODIC) {
	for(int GCell=1; GCell<=Nghost; ++GCell){
	  Node[INu+GCell][j].X = ( Node[INu][j].X +
				   (Node[INl+GCell][j].X - 
				    Node[INl][j].X) );
	}
      } /* endif */
    } /* endfor */
    
    // Update NW and NE corner nodes.
    for (int j = JNu ; j <= JNu+Nghost ; ++j) {
      if (BCtypeW[j-1] != BC_REFLECTION &&
	  BCtypeW[j-1] != BC_PERIODIC &&  
	  BCtypeW[j-1] != BC_NO_SLIP &&
	  BCtypeW[j-1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	  BCtypeW[j-1] != BC_WALL_VISCOUS_HEATFLUX &&
	  BCtypeW[j-1] != BC_MOVING_WALL &&
	  BCtypeW[j-1] != BC_MOVING_WALL_ISOTHERMAL &&
	  BCtypeW[j-1] != BC_MOVING_WALL_HEATFLUX &&
	  BCtypeW[j-1] != BC_BURNING_SURFACE &&
	  BCtypeW[j-1] != BC_MASS_INJECTION &&
	  BCtypeW[j-1] != BC_WALL_INVISCID) {
	for(int GCell=1; GCell<=Nghost; ++GCell){
	  Node[INl-GCell][j].X = ( Node[INl][j].X -
				   (Node[INl+GCell][j].X - 
				    Node[INl][j].X) );
	}
      } else if (BCtypeW[j-1] == BC_REFLECTION || 
		 BCtypeW[j-1] == BC_NO_SLIP ||
		 BCtypeW[j-1] == BC_WALL_VISCOUS_ISOTHERMAL ||
		 BCtypeW[j-1] == BC_WALL_VISCOUS_HEATFLUX ||
		 BCtypeW[j-1] == BC_MOVING_WALL ||
		 BCtypeW[j-1] == BC_MOVING_WALL_ISOTHERMAL ||
		 BCtypeW[j-1] == BC_MOVING_WALL_HEATFLUX ||
		 BCtypeW[j-1] == BC_BURNING_SURFACE ||
		 BCtypeW[j-1] == BC_MASS_INJECTION ||
		 BCtypeW[j-1] == BC_WALL_INVISCID) {
	if (j != JNu+Nghost) {
	  norm_dir = - HALF*(nfaceW(ICl, j) + 
			     nfaceW(ICl, j-1));
	} else {
	  norm_dir = - nfaceW(ICl, j-1); 
	} /* endif */
	for(int GCell=1; GCell<=Nghost; ++GCell){
	  X_norm = ((Node[INl+GCell][j].X - 
		     Node[INl][j].X) * norm_dir) * norm_dir;
	  X_tan = (Node[INl+GCell][j].X - 
		   Node[INl][j].X) - X_norm;
	  Node[INl-GCell][j].X = ( Node[INl][j].X -
				   X_norm + X_tan );
	}
      } else if (BCtypeW[j-1] == BC_PERIODIC) {
	for(int GCell=1; GCell<=Nghost; ++GCell){
	  Node[INl-GCell][j].X = ( Node[INl][j].X -
				   (Node[INu][j].X - 
				    Node[INu-GCell][j].X) );
	}
      } /* endif */
      
      if (BCtypeE[j-1] != BC_REFLECTION &&
	  BCtypeE[j-1] != BC_PERIODIC && 
	  BCtypeE[j-1] != BC_NO_SLIP &&
	  BCtypeE[j-1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	  BCtypeE[j-1] != BC_WALL_VISCOUS_HEATFLUX &&
	  BCtypeE[j-1] != BC_MOVING_WALL &&
	  BCtypeE[j-1] != BC_MOVING_WALL_ISOTHERMAL &&
	  BCtypeE[j-1] != BC_MOVING_WALL_HEATFLUX &&
	  BCtypeE[j-1] != BC_BURNING_SURFACE &&
	  BCtypeE[j-1] != BC_MASS_INJECTION &&
	  BCtypeE[j-1] != BC_WALL_INVISCID) {
	for(int GCell=1; GCell<=Nghost; ++GCell){
	  Node[INu+GCell][j].X = ( Node[INu][j].X +
				   (Node[INu][j].X - 
				    Node[INu-GCell][j].X) );
	}
      } else if (BCtypeE[j-1] == BC_REFLECTION || 
		 BCtypeE[j-1] == BC_NO_SLIP ||
		 BCtypeE[j-1] == BC_WALL_VISCOUS_ISOTHERMAL ||
		 BCtypeE[j-1] == BC_WALL_VISCOUS_HEATFLUX ||
		 BCtypeE[j-1] == BC_MOVING_WALL ||
		 BCtypeE[j-1] == BC_MOVING_WALL_ISOTHERMAL ||
		 BCtypeE[j-1] == BC_MOVING_WALL_HEATFLUX ||
		 BCtypeE[j-1] == BC_BURNING_SURFACE ||
		 BCtypeE[j-1] == BC_MASS_INJECTION ||
		 BCtypeE[j-1] == BC_WALL_INVISCID) {
	if (j != JNu+Nghost) {
	  norm_dir = HALF*(nfaceE(ICu, j) + 
			   nfaceE(ICu, j-1));
	} else {
	  norm_dir = nfaceE(ICu, j-1); 
	} /* endif */
	for(int GCell=1; GCell<=Nghost; ++GCell){
	  X_norm = ((Node[INu][j].X - 
		     Node[INu-GCell][j].X) * norm_dir) * norm_dir;
	  X_tan = (Node[INu][j].X - 
		   Node[INu-GCell][j].X) - X_norm;
	  Node[INu+GCell][j].X = ( Node[INu][j].X +
				   X_norm - X_tan );
	}
      } else if (BCtypeE[j-1] == BC_PERIODIC) {
	for(int GCell=1; GCell<=Nghost; ++GCell){
	  Node[INu+GCell][j].X = ( Node[INu][j].X +
				   (Node[INl+GCell][j].X - 
				    Node[INl][j].X) );
	}
      } /* endif */
    } /* endfor */

    // Require update of the ghost cells geometric properties.
    Schedule_Ghost_Cells_Update();

    // Update corners
    Update_Corner_Ghost_Nodes();

}


/*!
 * Updates the corner ghost nodes for the quadrilateral 
 * mesh block.
 *
 * \note This subroutine DOESN'T update the cell information.
 */
void Grid2D_Quad_Block_HO::Update_Corner_Ghost_Nodes(void) {
  
  Vector2D norm_dir, X_norm, X_tan;
  
  // SOUTH-WEST corner:
  if (BCtypeS[INl-1] == BC_NONE &&
      BCtypeW[JNl-1] == BC_NONE) {
    // Do nothing.
    
  } else if ((BCtypeS[INl-1] != BC_REFLECTION &&
	      BCtypeS[INl-1] != BC_PERIODIC &&
	      BCtypeS[INl-1] != BC_NO_SLIP &&
	      BCtypeS[INl-1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	      BCtypeS[INl-1] != BC_WALL_VISCOUS_HEATFLUX &&
	      BCtypeS[INl-1] != BC_MOVING_WALL &&
	      BCtypeS[INl-1] != BC_MOVING_WALL_ISOTHERMAL &&
	      BCtypeS[INl-1] != BC_MOVING_WALL_HEATFLUX &&
	      BCtypeS[INl-1] != BC_BURNING_SURFACE &&
	      BCtypeS[INl-1] != BC_MASS_INJECTION &&
	      BCtypeS[INl-1] != BC_RINGLEB_FLOW &&
	      (BCtypeS[INl-1] != BC_WALL_INVISCID || Mesh_Requiring_NonReflected_South_Ghost_Cells == ON ) ) &&
	     BCtypeW[JNl-1] == BC_NONE) {
    // Extrapolate cells south.
    for (int ng = 1; ng <= Nghost; ng++) {
      for (int i = INl-Nghost; i <= INl; i++) {
	Node[i][JNl-ng].X = (Node[i][JNl].X +
			     (Node[i][JNl].X - 
			      Node[i][JNl+ng].X) );
      }
    }
    
  } else if ((BCtypeS[INl-1] == BC_NONE ||
	      (BCtypeS[INl-1] != BC_REFLECTION &&
	       BCtypeS[INl-1] != BC_PERIODIC &&
	       BCtypeS[INl-1] != BC_NO_SLIP &&
	       BCtypeS[INl-1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	       BCtypeS[INl-1] != BC_WALL_VISCOUS_HEATFLUX &&
	       BCtypeS[INl-1] != BC_MOVING_WALL &&
	       BCtypeS[INl-1] != BC_MOVING_WALL_ISOTHERMAL &&
	       BCtypeS[INl-1] != BC_MOVING_WALL_HEATFLUX &&
	       BCtypeS[INl-1] != BC_BURNING_SURFACE &&
	       BCtypeS[INl-1] != BC_MASS_INJECTION &&
	       BCtypeS[INl-1] != BC_RINGLEB_FLOW &&
	       BCtypeS[INl-1] != BC_WALL_INVISCID)) &&
	     (BCtypeW[JNl-1] != BC_REFLECTION &&
	      BCtypeW[JNl-1] != BC_PERIODIC &&  
	      BCtypeW[JNl-1] != BC_NO_SLIP &&
	      BCtypeW[JNl-1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	      BCtypeW[JNl-1] != BC_WALL_VISCOUS_HEATFLUX &&
	      BCtypeW[JNl-1] != BC_MOVING_WALL &&
	      BCtypeW[JNl-1] != BC_MOVING_WALL_ISOTHERMAL &&
	      BCtypeW[JNl-1] != BC_MOVING_WALL_HEATFLUX &&
	      BCtypeW[JNl-1] != BC_BURNING_SURFACE &&
	      BCtypeW[JNl-1] != BC_MASS_INJECTION &&
	      BCtypeW[JNl-1] != BC_WALL_INVISCID)) {
    // Extrapolate cells west.
    for (int ng = 1; ng <= Nghost; ng++) {
      for (int j = JNl-Nghost; j <= JNl; j++) {
 	Node[INl-ng][j].X = ( Node[INl][j].X -
			      (Node[INl+ng][j].X - 
			       Node[INl][j].X) );
      }
    }
    
  } else if ((BCtypeS[INl-1] == BC_NONE ||
	      (BCtypeS[INl-1] != BC_REFLECTION &&
	       BCtypeS[INl-1] != BC_PERIODIC &&
	       BCtypeS[INl-1] != BC_NO_SLIP &&
	       BCtypeS[INl-1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	       BCtypeS[INl-1] != BC_WALL_VISCOUS_HEATFLUX &&
	       BCtypeS[INl-1] != BC_MOVING_WALL &&
	       BCtypeS[INl-1] != BC_MOVING_WALL_HEATFLUX &&
	       BCtypeS[INl-1] != BC_MOVING_WALL_ISOTHERMAL &&
	       BCtypeS[INl-1] != BC_BURNING_SURFACE &&
	       BCtypeS[INl-1] != BC_MASS_INJECTION &&
	       BCtypeS[INl-1] != BC_RINGLEB_FLOW &&
	       BCtypeS[INl-1] != BC_WALL_INVISCID)) &&
	     (BCtypeW[JNl-1] == BC_REFLECTION ||
	      BCtypeW[JNl-1] == BC_PERIODIC ||
	      BCtypeW[JNl-1] == BC_NO_SLIP ||
	      BCtypeW[JNl-1] == BC_WALL_VISCOUS_ISOTHERMAL ||
	      BCtypeW[JNl-1] == BC_WALL_VISCOUS_HEATFLUX ||
	      BCtypeW[JNl-1] == BC_MOVING_WALL ||
	      BCtypeW[JNl-1] == BC_MOVING_WALL_ISOTHERMAL ||
	      BCtypeW[JNl-1] == BC_MOVING_WALL_HEATFLUX ||
	      BCtypeW[JNl-1] == BC_BURNING_SURFACE ||
	      BCtypeW[JNl-1] == BC_MASS_INJECTION ||
	      BCtypeW[JNl-1] == BC_WALL_INVISCID)) {
    // Reflect cells west.
    for (int ng = 1; ng <= Nghost; ng++) {
      for (int j = JNl-Nghost; j <= JNl; j++) {
	if (j == JNl-Nghost) {
	  norm_dir = - nfaceW(ICl,JCl-Nghost);
	} else {
	  norm_dir = - HALF*(nfaceW(ICl,j) + 
			     nfaceW(ICl,j-1));
	}
	X_norm = ((Node[INl+ng][j].X - 
		   Node[INl][j].X)*norm_dir)*norm_dir;
	X_tan = (Node[INl+ng][j].X -
		 Node[INl][j].X) - X_norm;
	Node[INl-ng][j].X = ( Node[INl][j].X -
			      X_norm + X_tan );
      }
    }

  } else if (BCtypeS[INl-1] == BC_REFLECTION ||
	     BCtypeS[INl-1] == BC_PERIODIC ||
	     BCtypeS[INl-1] == BC_NO_SLIP ||
	     BCtypeS[INl-1] == BC_WALL_VISCOUS_ISOTHERMAL ||
	     BCtypeS[INl-1] == BC_WALL_VISCOUS_HEATFLUX ||
	     BCtypeS[INl-1] == BC_MOVING_WALL ||
	     BCtypeS[INl-1] == BC_MOVING_WALL_ISOTHERMAL ||
	     BCtypeS[INl-1] == BC_MOVING_WALL_HEATFLUX ||
	     BCtypeS[INl-1] == BC_BURNING_SURFACE ||
	     BCtypeS[INl-1] == BC_MASS_INJECTION ||
	     BCtypeS[INl-1] == BC_RINGLEB_FLOW || 
	     BCtypeS[INl-1] == BC_WALL_INVISCID) {
    // Reflect cells south.
    for (int ng = 1; ng <= Nghost; ng++) {
      //for (int i = INl-Nghost; i <= INl; i++) {
      for (int i = INl; i >= INl-Nghost; i--) {
	if (i == INl-Nghost) {
	  norm_dir = - nfaceS(ICl-Nghost,JCl);
	} else {
	  norm_dir = - HALF*(nfaceS(i,JCl) +
			     nfaceS(i-1,JCl));
	}
	X_norm = ((Node[i][JNl+ng].X - 
		   Node[i][JNl].X)*norm_dir)*norm_dir;
	X_tan = (Node[i][JNl+ng].X - 
		 Node[i][JNl].X) - X_norm;
	Node[i][JNl-ng].X = ( Node[i][JNl].X -
			      X_norm + X_tan );
      }
    }
  }

  // SOUTH-EAST corner:
  if (BCtypeS[INu+1] == BC_NONE &&
      BCtypeE[JNl-1] == BC_NONE) {
    // Do nothing.

  } else if ((BCtypeS[INu+1] != BC_REFLECTION &&
	      BCtypeS[INu+1] != BC_PERIODIC &&
	      BCtypeS[INu+1] != BC_NO_SLIP &&
	      BCtypeS[INu+1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	      BCtypeS[INu+1] != BC_WALL_VISCOUS_HEATFLUX &&
	      BCtypeS[INu+1] != BC_MOVING_WALL &&
	      BCtypeS[INu+1] != BC_MOVING_WALL_ISOTHERMAL &&
	      BCtypeS[INu+1] != BC_MOVING_WALL_HEATFLUX &&
	      BCtypeS[INu+1] != BC_BURNING_SURFACE &&
	      BCtypeS[INu+1] != BC_MASS_INJECTION &&
	      BCtypeS[INu+1] != BC_RINGLEB_FLOW &&
	      (BCtypeS[INu+1] != BC_WALL_INVISCID || Mesh_Requiring_NonReflected_South_Ghost_Cells == ON) ) &&
	     BCtypeE[JNl-1] == BC_NONE) {
    // Extrapolate cells south.
    for (int ng = 1; ng <= Nghost; ng++) {
      for (int i = INu; i <= INu+Nghost; i++) {
	Node[i][JNl-ng].X = ( Node[i][JNl].X +
			      (Node[i][JNl].X - 
			       Node[i][JNl+ng].X) );
      }
    }

  } else if ((BCtypeS[INu-1] == BC_NONE ||
	      (BCtypeS[INu+1] != BC_REFLECTION &&
	       BCtypeS[INu+1] != BC_PERIODIC &&
	       BCtypeS[INu+1] != BC_NO_SLIP &&
	       BCtypeS[INu+1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	       BCtypeS[INu+1] != BC_WALL_VISCOUS_HEATFLUX &&
	       BCtypeS[INu+1] != BC_MOVING_WALL &&
	       BCtypeS[INu+1] != BC_MOVING_WALL_ISOTHERMAL &&
	       BCtypeS[INu+1] != BC_MOVING_WALL_HEATFLUX &&
	       BCtypeS[INu+1] != BC_BURNING_SURFACE &&
	       BCtypeS[INu+1] != BC_MASS_INJECTION &&
	       BCtypeS[INu+1] != BC_RINGLEB_FLOW &&
	       BCtypeS[INu+1] != BC_WALL_INVISCID)) &&
	     (BCtypeE[JNl-1] != BC_REFLECTION &&
	      BCtypeE[JNl-1] != BC_PERIODIC &&  
	      BCtypeE[JNl-1] != BC_NO_SLIP &&
	      BCtypeE[JNl-1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	      BCtypeE[JNl-1] != BC_WALL_VISCOUS_HEATFLUX &&
	      BCtypeE[JNl-1] != BC_MOVING_WALL &&
	      BCtypeE[JNl-1] != BC_MOVING_WALL_ISOTHERMAL &&
	      BCtypeE[JNl-1] != BC_MOVING_WALL_HEATFLUX &&
	      BCtypeE[JNl-1] != BC_BURNING_SURFACE &&
	      BCtypeE[JNl-1] != BC_MASS_INJECTION &&
	      BCtypeE[JNl-1] != BC_WALL_INVISCID)) {
    // Extrapolate cells east.
    for (int ng = 1; ng <= Nghost; ng++) {
      for (int j = JNl-Nghost; j <= JNl; j++) {
	Node[INu+ng][j].X = ( Node[INu][j].X +
			      (Node[INu][j].X - 
			       Node[INu-ng][j].X) );
      }
    }

  } else if ((BCtypeS[INu+1] == BC_NONE ||
	      (BCtypeS[INu+1] != BC_REFLECTION &&
	       BCtypeS[INu+1] != BC_PERIODIC &&
	       BCtypeS[INu+1] != BC_NO_SLIP &&
	       BCtypeS[INu+1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	       BCtypeS[INu+1] != BC_WALL_VISCOUS_HEATFLUX &&
	       BCtypeS[INu+1] != BC_MOVING_WALL &&
	       BCtypeS[INu+1] != BC_MOVING_WALL_ISOTHERMAL &&
	       BCtypeS[INu+1] != BC_MOVING_WALL_HEATFLUX &&
	       BCtypeS[INu+1] != BC_BURNING_SURFACE &&
	       BCtypeS[INu+1] != BC_MASS_INJECTION &&
	       BCtypeS[INu+1] != BC_RINGLEB_FLOW &&
	       BCtypeS[INu+1] != BC_WALL_INVISCID)) &&
	     (BCtypeE[JNl-1] == BC_REFLECTION ||
	      BCtypeE[JNl-1] == BC_PERIODIC ||
	      BCtypeE[JNl-1] == BC_NO_SLIP ||
	      BCtypeE[JNl-1] == BC_WALL_VISCOUS_ISOTHERMAL ||
	      BCtypeE[JNl-1] == BC_WALL_VISCOUS_HEATFLUX ||
	      BCtypeE[JNl-1] == BC_MOVING_WALL ||
	      BCtypeE[JNl-1] == BC_MOVING_WALL_ISOTHERMAL ||
	      BCtypeE[JNl-1] == BC_MOVING_WALL_HEATFLUX ||
	      BCtypeE[JNl-1] == BC_BURNING_SURFACE ||
	      BCtypeE[JNl-1] == BC_MASS_INJECTION ||
	      BCtypeE[JNl-1] == BC_WALL_INVISCID)) {
    // Reflect cells east.
    for (int ng = 1; ng <= Nghost; ng++) {
      for (int j = JNl; j >= JNl-Nghost; j--) {
	if (j == JNl-Nghost) {
	  norm_dir = nfaceE(ICu,JCl-Nghost);
	} else {
	  norm_dir = HALF*(nfaceE(ICu,j) + 
			   nfaceE(ICu,j+1));
	}
	norm_dir = nfaceE(ICu,j);
	X_norm = ((Node[INu][j].X - 
		   Node[INu-ng][j].X)*norm_dir)*norm_dir;
	X_tan = (Node[INu][j].X - 
		 Node[INu-ng][j].X) - X_norm;
	Node[INu+ng][j].X = ( Node[INu][j].X +
			      X_norm - X_tan );
      }
    }

  } else if (BCtypeS[INu+1] == BC_REFLECTION ||
	     BCtypeS[INu+1] == BC_PERIODIC ||
	     BCtypeS[INu+1] == BC_NO_SLIP ||
	     BCtypeS[INu+1] == BC_WALL_VISCOUS_ISOTHERMAL ||
	     BCtypeS[INu+1] == BC_WALL_VISCOUS_HEATFLUX ||
	     BCtypeS[INu+1] == BC_MOVING_WALL ||
	     BCtypeS[INu+1] == BC_MOVING_WALL_ISOTHERMAL ||
	     BCtypeS[INu+1] == BC_MOVING_WALL_HEATFLUX ||
	     BCtypeS[INu+1] == BC_BURNING_SURFACE ||
	     BCtypeS[INu+1] == BC_MASS_INJECTION ||
	     BCtypeS[INu+1] == BC_RINGLEB_FLOW ||
	     BCtypeS[INu+1] == BC_WALL_INVISCID) {
    // Reflect cells south.
    for (int ng = 1; ng <= Nghost; ng++) {
      for (int i = INu; i <= INu+Nghost; i++) {
	if (i == INu+Nghost) {
	  norm_dir = -nfaceS(ICu+Nghost,JCl);
	} else {
	  norm_dir = -HALF*(nfaceS(i,JCl) +
			    nfaceS(i-1,JCl));
	}
 	X_norm = ((Node[i][JNl+ng].X - 
 		   Node[i][JNl].X)*norm_dir)*norm_dir;
 	X_tan = (Node[i][JNl+ng].X - 
 		 Node[i][JNl].X) - X_norm;
  	Node[i][JNl-ng].X = ( Node[i][JNl].X -
			      X_norm + X_tan );
      }
    }
  }

  // NORTH-WEST corner:
  if (BCtypeN[INl-1] == BC_NONE &&
      BCtypeW[JNu+1] == BC_NONE) {
    // Do nothing.

  } else if ((BCtypeN[INl-1] != BC_REFLECTION &&
	      BCtypeN[INl-1] != BC_PERIODIC &&
	      BCtypeN[INl-1] != BC_NO_SLIP &&
	      BCtypeN[INl-1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	      BCtypeN[INl-1] != BC_WALL_VISCOUS_HEATFLUX &&
	      BCtypeN[INl-1] != BC_MOVING_WALL &&
	      BCtypeN[INl-1] != BC_MOVING_WALL_ISOTHERMAL &&
	      BCtypeN[INl-1] != BC_MOVING_WALL_HEATFLUX &&
	      BCtypeN[INl-1] != BC_BURNING_SURFACE &&
	      BCtypeN[INl-1] != BC_MASS_INJECTION &&
	      BCtypeN[INl-1] != BC_WALL_INVISCID) &&
	     BCtypeW[JNu+1] == BC_NONE) {
    // Extrapolate cells north.
    for (int ng = 1; ng <= Nghost; ng++) {
      for (int i = INl-Nghost; i <= INl; i++) {
	Node[i][JNu+ng].X = ( Node[i][JNu].X +
			      (Node[i][JNu].X - 
			       Node[i][JNu-ng].X) );
      }
    }

  } else if ((BCtypeN[INl-1] == BC_NONE ||
	      (BCtypeN[INl-1] != BC_REFLECTION &&
	       BCtypeN[INl-1] != BC_PERIODIC &&
	       BCtypeN[INl-1] != BC_NO_SLIP &&
	       BCtypeN[INl-1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	       BCtypeN[INl-1] != BC_WALL_VISCOUS_HEATFLUX &&
	       BCtypeN[INl-1] != BC_MOVING_WALL &&
	       BCtypeN[INl-1] != BC_MOVING_WALL_ISOTHERMAL &&
	       BCtypeN[INl-1] != BC_MOVING_WALL_HEATFLUX &&
	       BCtypeN[INl-1] != BC_BURNING_SURFACE &&
	       BCtypeN[INl-1] != BC_MASS_INJECTION &&
	       BCtypeN[INl-1] != BC_WALL_INVISCID)) &&
	     (BCtypeW[JNu+1] != BC_REFLECTION &&
	      BCtypeW[JNu+1] != BC_PERIODIC &&  
	      BCtypeW[JNu+1] != BC_NO_SLIP &&
	      BCtypeW[JNu+1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	      BCtypeW[JNu+1] != BC_WALL_VISCOUS_HEATFLUX &&
	      BCtypeW[JNu+1] != BC_MOVING_WALL &&
	      BCtypeW[JNu+1] != BC_MOVING_WALL_ISOTHERMAL &&
	      BCtypeW[JNu+1] != BC_MOVING_WALL_HEATFLUX &&
	      BCtypeW[JNu+1] != BC_BURNING_SURFACE &&
	      BCtypeW[JNu+1] != BC_MASS_INJECTION &&
	      BCtypeW[JNu+1] != BC_WALL_INVISCID)) {
    // Extrapolate cells west.
    for (int ng = 1; ng <= Nghost; ng++) {
      for (int j = JNu; j <= JNu+Nghost; j++) {
	Node[INl-ng][j].X = ( Node[INl][j].X -
			      (Node[INl+ng][j].X - 
			       Node[INl][j].X) );
      }
    }

  } else if ((BCtypeN[INl-1] == BC_NONE ||
	      (BCtypeN[INl-1] != BC_REFLECTION &&
	       BCtypeN[INl-1] != BC_PERIODIC &&
	       BCtypeN[INl-1] != BC_NO_SLIP &&
	       BCtypeN[INl-1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	       BCtypeN[INl-1] != BC_WALL_VISCOUS_HEATFLUX &&
	       BCtypeN[INl-1] != BC_MOVING_WALL &&
	       BCtypeN[INl-1] != BC_MOVING_WALL_ISOTHERMAL &&
	       BCtypeN[INl-1] != BC_MOVING_WALL_HEATFLUX &&
	       BCtypeN[INl-1] != BC_BURNING_SURFACE &&
	       BCtypeN[INl-1] != BC_MASS_INJECTION &&
	       BCtypeN[INl-1] != BC_WALL_INVISCID)) &&
	     (BCtypeW[JNu+1] == BC_REFLECTION ||
	      BCtypeW[JNu+1] == BC_PERIODIC ||
	      BCtypeW[JNu+1] == BC_NO_SLIP ||
	      BCtypeW[JNu+1] == BC_WALL_VISCOUS_ISOTHERMAL ||
	      BCtypeW[JNu+1] == BC_WALL_VISCOUS_HEATFLUX ||
	      BCtypeW[JNu+1] == BC_MOVING_WALL ||
	      BCtypeW[JNu+1] == BC_MOVING_WALL_ISOTHERMAL ||
	      BCtypeW[JNu+1] == BC_MOVING_WALL_HEATFLUX ||
	      BCtypeW[JNu+1] == BC_BURNING_SURFACE ||
	      BCtypeW[JNu+1] == BC_MASS_INJECTION ||
	      BCtypeW[JNu+1] == BC_WALL_INVISCID)) {
    // Reflect cells west.
    for (int ng = 1; ng <= Nghost; ng++) {
      for (int j = JNu; j <= JNu+Nghost; j++) {
	if (j == JNu+Nghost) {
	  norm_dir = - nfaceW(ICl,JCu+Nghost);
	} else {
	  norm_dir = - HALF*(nfaceW(ICl,j) + 
			     nfaceW(ICl,j-1));
	}
	X_norm = ((Node[INl+ng][j].X -
		   Node[INl][j].X)*norm_dir)*norm_dir;
	X_tan = (Node[INl+ng][j].X - 
		 Node[INl][j].X) - X_norm;
	Node[INl-ng][j].X = ( Node[INl][j].X -
			      X_norm + X_tan );
      }
    }

  } else if (BCtypeN[INl-1] == BC_REFLECTION ||
	     BCtypeN[INl-1] == BC_PERIODIC ||
	     BCtypeN[INl-1] == BC_NO_SLIP ||
	     BCtypeN[INl-1] == BC_WALL_VISCOUS_ISOTHERMAL ||
	     BCtypeN[INl-1] == BC_WALL_VISCOUS_HEATFLUX ||
	     BCtypeN[INl-1] == BC_MOVING_WALL ||
	     BCtypeN[INl-1] == BC_MOVING_WALL_ISOTHERMAL ||
	     BCtypeN[INl-1] == BC_MOVING_WALL_HEATFLUX ||
	     BCtypeN[INl-1] == BC_BURNING_SURFACE ||
	     BCtypeN[INl-1] == BC_MASS_INJECTION ||
	     BCtypeN[INl-1] == BC_WALL_INVISCID) {
    // Reflect cells north.
    for (int ng = 1; ng <= Nghost; ng++) {
      for (int i = INl-Nghost; i <= INl; i++) {
	if (i == INl-Nghost) {
	  norm_dir = -nfaceN(ICl-Nghost,JCu);
	} else {
	  norm_dir = -HALF*(nfaceN(i,JCu) + 
			    nfaceN(i-1,JCu));
	}
	X_norm = ((Node[i][JNu].X - 
		   Node[i][JNu-ng].X) * norm_dir) * norm_dir;
	X_tan = (Node[i][JNu].X - 
		 Node[i][JNu-ng].X) - X_norm;
	Node[i][JNu+ng].X = ( Node[i][JNu].X +
			      X_norm - X_tan );
      }
    }
  }

  // NORTH-EAST corner:
  if (BCtypeN[INu+1] == BC_NONE &&
      BCtypeE[JNu+1] == BC_NONE) {
    // Do nothing.

  } else if ((BCtypeN[INu+1] != BC_REFLECTION &&
	      BCtypeN[INu+1] != BC_PERIODIC &&
	      BCtypeN[INu+1] != BC_NO_SLIP &&
	      BCtypeN[INu+1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	      BCtypeN[INu+1] != BC_WALL_VISCOUS_HEATFLUX &&
	      BCtypeN[INu+1] != BC_MOVING_WALL &&
	      BCtypeN[INu+1] != BC_MOVING_WALL_ISOTHERMAL &&
	      BCtypeN[INu+1] != BC_MOVING_WALL_HEATFLUX &&
	      BCtypeN[INu+1] != BC_BURNING_SURFACE &&
	      BCtypeN[INu+1] != BC_MASS_INJECTION &&
	      BCtypeN[INu+1] != BC_WALL_INVISCID) &&
	     BCtypeE[JNu+1] == BC_NONE) {
    // Extrapolate cells north.
    for (int ng = 1; ng <= Nghost; ng++) {
      for (int i = INu; i <= INu+Nghost; i++) {
	Node[i][JNu+ng].X = ( Node[i][JNu].X +
			      (Node[i][JNu].X - 
			       Node[i][JNu-ng].X) );
      }
    }

  } else if ((BCtypeN[INu+1] == BC_NONE ||
	      (BCtypeN[INu+1] != BC_REFLECTION &&
	       BCtypeN[INu+1] != BC_PERIODIC &&
	       BCtypeN[INu+1] != BC_NO_SLIP &&
	       BCtypeN[INu+1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	       BCtypeN[INu+1] != BC_WALL_VISCOUS_HEATFLUX &&
	       BCtypeN[INu+1] != BC_MOVING_WALL &&
	       BCtypeN[INu+1] != BC_MOVING_WALL_ISOTHERMAL &&
	       BCtypeN[INu+1] != BC_MOVING_WALL_HEATFLUX &&
	       BCtypeN[INu+1] != BC_BURNING_SURFACE &&
	       BCtypeN[INu+1] != BC_MASS_INJECTION &&
	       BCtypeN[INu+1] != BC_WALL_INVISCID)) &&
	     (BCtypeE[JNu+1] != BC_REFLECTION &&
	      BCtypeE[JNu+1] != BC_PERIODIC &&  
	      BCtypeE[JNu+1] != BC_NO_SLIP &&
	      BCtypeE[JNu+1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	      BCtypeE[JNu+1] != BC_WALL_VISCOUS_HEATFLUX &&
	      BCtypeE[JNu+1] != BC_MOVING_WALL &&
	      BCtypeE[JNu+1] != BC_MOVING_WALL_ISOTHERMAL &&
	      BCtypeE[JNu+1] != BC_MOVING_WALL_HEATFLUX &&
	      BCtypeE[JNu+1] != BC_BURNING_SURFACE &&
	      BCtypeE[JNu+1] != BC_MASS_INJECTION &&
	      BCtypeE[JNu+1] != BC_WALL_INVISCID)) {
    // Extrapolate cells east.
    for (int ng = 1; ng <= Nghost; ng++) {
      for (int j = JNu; j <= JNu+Nghost; j++) {
	Node[INu+ng][j].X = ( Node[INu][j].X +
			      (Node[INu][j].X - 
			       Node[INu-ng][j].X) );
      }
    }

  } else if ((BCtypeN[INu+1] == BC_NONE ||
	      (BCtypeN[INu+1] != BC_REFLECTION &&
	       BCtypeN[INu+1] != BC_PERIODIC &&
	       BCtypeN[INu+1] != BC_NO_SLIP &&
	       BCtypeN[INu+1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	       BCtypeN[INu+1] != BC_WALL_VISCOUS_HEATFLUX &&
	       BCtypeN[INu+1] != BC_MOVING_WALL &&
	       BCtypeN[INu+1] != BC_MOVING_WALL_ISOTHERMAL &&
	       BCtypeN[INu+1] != BC_MOVING_WALL_HEATFLUX &&
	       BCtypeN[INu+1] != BC_BURNING_SURFACE &&
	       BCtypeN[INu+1] != BC_MASS_INJECTION &&
	       BCtypeN[INu+1] != BC_WALL_INVISCID)) &&
	     (BCtypeE[JNu+1] == BC_REFLECTION ||
	      BCtypeE[JNu+1] == BC_PERIODIC ||
	      BCtypeE[JNu+1] == BC_NO_SLIP ||
	      BCtypeE[JNu+1] == BC_WALL_VISCOUS_ISOTHERMAL ||
	      BCtypeE[JNu+1] == BC_WALL_VISCOUS_HEATFLUX ||
	      BCtypeE[JNu+1] == BC_MOVING_WALL ||
	      BCtypeE[JNu+1] == BC_MOVING_WALL_ISOTHERMAL ||
	      BCtypeE[JNu+1] == BC_MOVING_WALL_HEATFLUX ||
	      BCtypeE[JNu+1] == BC_BURNING_SURFACE ||
	      BCtypeE[JNu+1] == BC_MASS_INJECTION ||
	      BCtypeE[JNu+1] == BC_WALL_INVISCID)) {
    // Reflect cells east.
    for (int ng = 1; ng <= Nghost; ng++) {
      for (int j = JNu; j <= JNu+Nghost; j++) {
	if (j == JNu+Nghost) {
	  norm_dir = nfaceE(ICu,JCu+Nghost);
	} else {
	  norm_dir = HALF*(nfaceE(ICu,j) + 
			   nfaceE(ICu,j-1));
	}
	X_norm = ((Node[INu][j].X - 
		   Node[INu-ng][j].X)*norm_dir)*norm_dir;
	X_tan = (Node[INu][j].X - 
		 Node[INu-ng][j].X) - X_norm;
	Node[INu+ng][j].X = ( Node[INu][j].X +
			      X_norm - X_tan );
      }
    }

  } else if (BCtypeN[INu+1] == BC_REFLECTION ||
	     BCtypeN[INu+1] == BC_PERIODIC ||
	     BCtypeN[INu+1] == BC_NO_SLIP ||
	     BCtypeN[INu+1] == BC_WALL_VISCOUS_ISOTHERMAL ||
	     BCtypeN[INu+1] == BC_WALL_VISCOUS_HEATFLUX ||
	     BCtypeN[INu+1] == BC_MOVING_WALL ||
	     BCtypeN[INu+1] == BC_MOVING_WALL_ISOTHERMAL ||
	     BCtypeN[INu+1] == BC_MOVING_WALL_HEATFLUX ||
	     BCtypeN[INu+1] == BC_BURNING_SURFACE ||
	     BCtypeN[INu+1] == BC_MASS_INJECTION ||
	     BCtypeN[INu+1] == BC_WALL_INVISCID) {
    // Reflect cells north.
    for (int ng = 1; ng <= Nghost; ng++) {
      for (int i = INu; i <= INu+Nghost; i++) {
	if (i == INu+Nghost) {
	  norm_dir = -nfaceN(ICu+Nghost,JCu);
	} else {
	  norm_dir = -HALF*(nfaceN(i,JCu) + 
			    nfaceN(i-1,JCu));
	}
	X_norm = ((Node[i][JNu].X - 
		   Node[i][JNu-ng].X) * norm_dir) * norm_dir;
	X_tan = (Node[i][JNu].X - 
		 Node[i][JNu-ng].X) - X_norm;
	Node[i][JNu+ng].X = ( Node[i][JNu].X +
			      X_norm - X_tan );
      }
    }
  }

  // Require the update of the corner ghost cells' geometric properties
  Schedule_Corner_Ghost_Cells_Update();

}

/*!
 * Updates only the spline info(s) for the quadrilateral mesh block.
 */
void Grid2D_Quad_Block_HO::Update_SplineInfos(void){

  int i, IndexShift;
  int NumGQPsPerSubinterval(NumGQP);

#ifdef MODIFY_NUMBER_OF_FLUX_CALCULATION_POINTS_AT_BOUNDARIES
  NumGQPsPerSubinterval = NumGQP + 1;
#endif

  if ( !CheckExistenceOfCurvedBoundaries() ){
    // There is no need for curved boundary representation so SplineInfo(s) don't need to be updated 
    return;

  } else {
    
    Spline2D_HO SplineCopy;
    
    // Determine the geometric properties along the splines (e.g. Gauss Quadrature point locations,
    // normals, and spline segment length at each cell)

    // ==== Check the North boundary
    if (IsNorthBoundaryCurved()){

      // Determine the geometric properties along the splines (e.g. Gauss Quadrature point locations, normals, etc.)
      // Update BndNorthSplineInfo[]
      // Check if memory is allocated for BndNorthSplineInfo
      if(BndNorthSplineInfo == NULL){ // allocate array
	BndNorthSplineInfo = new Spline2DInterval_HO [NCi];
      }

      // Check if the North Spline is defined such that the normals at the GaussQuadratures point outside of the domain 
      // (i.e The spline pathlength increases from INu to INl)
      if ( BndNorthSpline.getS(Node[INl][JNu]) > BndNorthSpline.getS(Node[INu][JNu]) ){
	// Update the geometric information 
	for(i=ICl; i<=ICu; ++i){
	  BndNorthSplineInfo[i].UpdateInterval(BndNorthSpline,Node[i][JNu],Node[i+1][JNu],NumGQPsPerSubinterval);
	}
      } else {
	// Copy the spline
	SplineCopy = BndNorthSpline;
	// Change the direction of increasing the pathlength
	SplineCopy.Reverse_Spline();

	// Update the geometric information 
	for(i=ICl; i<=ICu; ++i){
	  BndNorthSplineInfo[i].UpdateInterval(SplineCopy,Node[i][JNu],Node[i+1][JNu],NumGQPsPerSubinterval);
	}
      }//endif
    } // endif (North Boundary)

    // ==== Check the South boundary
    if (IsSouthBoundaryCurved()){

      // Update BndSouthSplineInfo[]
      // Check if memory is allocated for BndSouthSplineInfo
      if(BndSouthSplineInfo == NULL){ // allocate array
	BndSouthSplineInfo = new Spline2DInterval_HO [NCi];
      }

      // Check if the South Spline is defined such that the normals at the GaussQuadratures point outside of the domain 
      // (i.e The spline pathlength increases from INl to INu)
      if ( BndSouthSpline.getS(Node[INl][JNl]) < BndSouthSpline.getS(Node[INu][JNl]) ){
	// Update the geometric information 
	for(i=ICl; i<=ICu; ++i){
	  BndSouthSplineInfo[i].UpdateInterval(BndSouthSpline,Node[i][JNl],Node[i+1][JNl],NumGQPsPerSubinterval);
	}
      } else {
	// Copy the spline
	SplineCopy = BndSouthSpline;
	// Change the direction of increasing the pathlength
	SplineCopy.Reverse_Spline();

	// Update the geometric information 
	for(i=ICl; i<=ICu; ++i){
	  BndSouthSplineInfo[i].UpdateInterval(SplineCopy,Node[i][JNl],Node[i+1][JNl],NumGQPsPerSubinterval);
	}
      }//endif
    } // endif (South Boundary)

    // ==== Check the East boundary
    if (IsEastBoundaryCurved()){

      // Update BndEastSplineInfo[]
      // Check if memory is allocated for BndEastSplineInfo
      if(BndEastSplineInfo == NULL){ // allocate array
	BndEastSplineInfo = new Spline2DInterval_HO [NCj];
      }

      // Check if the East Spline is defined such that the normals at the GaussQuadratures point outside of the domain 
      // (i.e The spline pathlength increases from JNl to JNu)
      if ( BndEastSpline.getS(Node[INu][JNl]) < BndEastSpline.getS(Node[INu][JNu]) ){
	// Update the geometric information 
	for(i=JCl; i<=JCu; ++i){
	  BndEastSplineInfo[i].UpdateInterval(BndEastSpline,Node[INu][i],Node[INu][i+1],NumGQPsPerSubinterval);
	}
      } else {
	// Copy the spline
	SplineCopy = BndEastSpline;
	// Change the direction of increasing the pathlength
	SplineCopy.Reverse_Spline();

	// Update the geometric information 
	for(i=JCl; i<=JCu; ++i){
	  BndEastSplineInfo[i].UpdateInterval(SplineCopy,Node[INu][i],Node[INu][i+1],NumGQPsPerSubinterval);
	}
      }//endif
    } // endif (East Boundary)

    // ==== Check the West boundary
    if (IsWestBoundaryCurved()){

      // Update BndWestSplineInfo[]
      // Check if memory is allocated for BndWestSplineInfo
      if(BndWestSplineInfo == NULL){ // allocate array
	BndWestSplineInfo = new Spline2DInterval_HO [NCj];
      }

      // Check if the West Spline is defined such that the normals at the GaussQuadratures point outside of the domain 
      // (i.e The spline pathlength increases from JNu to JNl)
      if ( BndWestSpline.getS(Node[INl][JNl]) > BndWestSpline.getS(Node[INl][JNu]) ){
	// Update the geometric information 
	for(i=JCl; i<=JCu; ++i){
	  BndWestSplineInfo[i].UpdateInterval(BndWestSpline,Node[INl][i],Node[INl][i+1],NumGQPsPerSubinterval);
	}
      } else {
	// Copy the spline
	SplineCopy = BndWestSpline;
	// Change the direction of increasing the pathlength
	SplineCopy.Reverse_Spline();

	// Update the geometric information 
	for(i=JCl; i<=JCu; ++i){
	  BndWestSplineInfo[i].UpdateInterval(SplineCopy,Node[INl][i],Node[INl][i+1],NumGQPsPerSubinterval);
	}
      }//endif
    } // endif (West Boundary)

    /* Update the spline infos for the extension splines only if the ghost cells update is required.
     * This condition is necessary to prevent update of Infos before the ghost nodes are updated
     * and also extra computational time.
     */ 
    if (GhostCellsUpdate == ON){
      Update_ExtensionSplineInfos();
    }

  } // endif (CheckExistenceOfCurvedBoundaries)
}

/*!
 * Updates only the extension spline info(s) for the quadrilateral mesh block.
 */
void Grid2D_Quad_Block_HO::Update_ExtensionSplineInfos(void){

  int i, IndexShift;
  int NumGQPsPerSubinterval(NumGQP);

#ifdef MODIFY_NUMBER_OF_FLUX_CALCULATION_POINTS_AT_BOUNDARIES
  NumGQPsPerSubinterval = NumGQP + 1;
#endif

  if ( !CheckExistenceOfCurvedBoundaries() ){
    // There is no need for curved boundary representation so SplineInfo(s) don't need to be updated 
    return;

  } else {
    
    Spline2D_HO SplineCopy;
    
    // Determine the geometric properties along the splines (e.g. Gauss Quadrature point locations,
    // normals, and spline segment length at each cell)

    // Check the North boundary extension to West
    if ( IsWestExtendNorthBoundaryCurved() ){

      // Determine the geometric properties along the spline (e.g. Gauss Quadrature point locations, normals, etc.)
      // Update ExtendWest_BndNorthSplineInfo[]
      // Check if memory is allocated for ExtendWest_BndNorthSplineInfo
      if(ExtendWest_BndNorthSplineInfo == NULL){ // allocate array
	ExtendWest_BndNorthSplineInfo = new Spline2DInterval_HO [Nghost];
      }

      // Check if the North Spline extension to West is defined such that the normals point outside of the domain
      if ( ExtendWest_BndNorthSpline.getS(Node[INl-1][JNu]) > ExtendWest_BndNorthSpline.getS(Node[INl][JNu])){
	// Update the geometric information 
	for(i=0; i<=ICl-1; ++i){
	  ExtendWest_BndNorthSplineInfo[i].UpdateInterval(ExtendWest_BndNorthSpline,
							  Node[i][JNu],Node[i+1][JNu],NumGQPsPerSubinterval);
	}
      } else {
	// Copy the spline
	SplineCopy = ExtendWest_BndNorthSpline;
	// Change the direction of increasing the pathlength
	SplineCopy.Reverse_Spline();
	  
	// Update the geometric information 
	for(i=0; i<ICl; ++i){
	  ExtendWest_BndNorthSplineInfo[i].UpdateInterval(SplineCopy,Node[i][JNu],Node[i+1][JNu],NumGQPsPerSubinterval);
	}
      }//endif
    }//endif 

    // Check the North boundary extension to East
    if ( IsEastExtendNorthBoundaryCurved() ){

      // Determine the geometric properties along the spline (e.g. Gauss Quadrature point locations, normals, etc.)
      // Update ExtendEast_BndNorthSplineInfo[]
      // Check if memory is allocated for ExtendEast_BndNorthSplineInfo
      if(ExtendEast_BndNorthSplineInfo == NULL){ // allocate array
	ExtendEast_BndNorthSplineInfo = new Spline2DInterval_HO [Nghost];
      }

      // Check if the North Spline extension to East is defined such that the normals point outside of the domain
      if ( ExtendEast_BndNorthSpline.getS(Node[INu][JNu]) > ExtendEast_BndNorthSpline.getS(Node[INu+1][JNu])){
	// Update the geometric information 
	IndexShift = ICu + 1;	// due to the fact that only Nghost cells are stored
	for(i=ICu+1; i<=ICu+Nghost; ++i){
	  ExtendEast_BndNorthSplineInfo[i-IndexShift].UpdateInterval(ExtendEast_BndNorthSpline,
								     Node[i][JNu],Node[i+1][JNu],NumGQPsPerSubinterval);
	}
      } else {
	// Copy the spline
	SplineCopy = ExtendEast_BndNorthSpline;
	// Change the direction of increasing the pathlength
	SplineCopy.Reverse_Spline();
	  
	// Update the geometric information 
	IndexShift = ICu + 1;
	for(i=ICu+1; i<=ICu+Nghost; ++i){
	  ExtendEast_BndNorthSplineInfo[i-IndexShift].UpdateInterval(SplineCopy,
								     Node[i][JNu],Node[i+1][JNu],NumGQPsPerSubinterval);
	}
      }//endif
    }//endif 
 
    // Check the South boundary extension to West
    if ( IsWestExtendSouthBoundaryCurved() ){

      // Determine the geometric properties along the spline (e.g. Gauss Quadrature point locations, normals, etc.)
      // Update ExtendWest_BndSouthSplineInfo[]
      // Check if memory is allocated for ExtendWest_BndSouthSplineInfo
      if(ExtendWest_BndSouthSplineInfo == NULL){ // allocate array
	ExtendWest_BndSouthSplineInfo = new Spline2DInterval_HO [Nghost];
      }
      
      // Check if the South Spline extension to West is defined such that the normals point outside of the domain
      if ( ExtendWest_BndSouthSpline.getS(Node[INl-1][JNl]) < ExtendWest_BndSouthSpline.getS(Node[INl][JNl])){
	// Update the geometric information 
	for(i=0; i<=ICl-1; ++i){
	  ExtendWest_BndSouthSplineInfo[i].UpdateInterval(ExtendWest_BndSouthSpline,
							  Node[i][JNl],Node[i+1][JNl],NumGQPsPerSubinterval);
	}
      } else {
	// Copy the spline
	SplineCopy = ExtendWest_BndSouthSpline;
	// Change the direction of increasing the pathlength
	SplineCopy.Reverse_Spline();

	// Update the geometric information 
	for(i=0; i<=ICl-1; ++i){
	  ExtendWest_BndSouthSplineInfo[i].UpdateInterval(SplineCopy,Node[i][JNl],Node[i+1][JNl],NumGQPsPerSubinterval);
	}
      }//endif
    }//endif

    // Check the South boundary extension to East
    if ( IsEastExtendSouthBoundaryCurved() ){

      // Determine the geometric properties along the spline (e.g. Gauss Quadrature point locations, normals, etc.)
      // Update ExtendEast_BndSouthSplineInfo[]
      // Check if memory is allocated for ExtendEast_BndSouthSplineInfo
      if(ExtendEast_BndSouthSplineInfo == NULL){ // allocate array
	ExtendEast_BndSouthSplineInfo = new Spline2DInterval_HO [Nghost];
      }

      // Check if the South Spline extension to East is defined such that the normals point outside of the domain
      if ( ExtendEast_BndSouthSpline.getS(Node[INu][JNl]) < ExtendEast_BndSouthSpline.getS(Node[INu+1][JNl])){
	// Update the geometric information 
	IndexShift = ICu + 1;	// due to the fact that only Nghost cells are stored
	for(i=ICu+1; i<=ICu+Nghost; ++i){
	  ExtendEast_BndSouthSplineInfo[i-IndexShift].UpdateInterval(ExtendEast_BndSouthSpline,
								     Node[i][JNl],Node[i+1][JNl],NumGQPsPerSubinterval);
	}
      } else {
	// Copy the spline
	SplineCopy = ExtendEast_BndSouthSpline;
	// Change the direction of increasing the pathlength
	SplineCopy.Reverse_Spline();

	// Update the geometric information 
	IndexShift = ICu + 1;
	for(i=ICu+1; i<=ICu+Nghost; ++i){
	  ExtendEast_BndSouthSplineInfo[i-IndexShift].UpdateInterval(SplineCopy,
								     Node[i][JNl],Node[i+1][JNl],NumGQPsPerSubinterval);
	}
      }//endif
    }//endif

    // Check the East boundary extension to North
    if ( IsNorthExtendEastBoundaryCurved() ){

      // Update ExtendNorth_BndEastSplineInfo[]
      // Check if memory is allocated for ExtendNorth_BndEastSplineInfo
      if(ExtendNorth_BndEastSplineInfo == NULL){ // allocate array
	ExtendNorth_BndEastSplineInfo = new Spline2DInterval_HO [Nghost];
      }

      // Check if the East Spline extension to North is defined such that the normals point outside of the domain 
      if ( ExtendNorth_BndEastSpline.getS(Node[INu][JNu]) < ExtendNorth_BndEastSpline.getS(Node[INu][JNu+1]) ){
	// Update the geometric information 
	IndexShift = JCu + 1;
	for(i=JCu+1; i<=JCu+Nghost; ++i){
	  ExtendNorth_BndEastSplineInfo[i-IndexShift].UpdateInterval(ExtendNorth_BndEastSpline,
								     Node[INu][i],Node[INu][i+1],NumGQPsPerSubinterval);
	}
      } else {
	// Copy the spline
	SplineCopy = ExtendNorth_BndEastSpline;
	// Change the direction of increasing the pathlength
	SplineCopy.Reverse_Spline();

	// Update the geometric information 
	IndexShift = JCu + 1;
	for(i=JCu+1; i<=JCu+Nghost; ++i){
	  ExtendNorth_BndEastSplineInfo[i-IndexShift].UpdateInterval(SplineCopy,
								     Node[INu][i],Node[INu][i+1],NumGQPsPerSubinterval);
	}
      }//endif
    } // endif

    // Check the East boundary extension to South
    if ( IsSouthExtendEastBoundaryCurved() ){

      // Update ExtendSouth_BndEastSplineInfo[]
      // Check if memory is allocated for ExtendSouth_BndEastSplineInfo
      if(ExtendSouth_BndEastSplineInfo == NULL){ // allocate array
	ExtendSouth_BndEastSplineInfo = new Spline2DInterval_HO [Nghost];
      }

      // Check if the East Spline extension to South is defined such that the normals point outside of the domain 
      if ( ExtendSouth_BndEastSpline.getS(Node[INu][JNl-1]) < ExtendSouth_BndEastSpline.getS(Node[INu][JNl]) ){
	// Update the geometric information 
	for(i=0; i<=JCl-1; ++i){
	  ExtendSouth_BndEastSplineInfo[i].UpdateInterval(ExtendSouth_BndEastSpline,
							  Node[INu][i],Node[INu][i+1],NumGQPsPerSubinterval);
	}
      } else {
	// Copy the spline
	SplineCopy = ExtendSouth_BndEastSpline;
	// Change the direction of increasing the pathlength
	SplineCopy.Reverse_Spline();

	// Update the geometric information 
	for(i=0; i<=JCl-1; ++i){
	  ExtendSouth_BndEastSplineInfo[i].UpdateInterval(SplineCopy,Node[INu][i],Node[INu][i+1],NumGQPsPerSubinterval);
	}
      }//endif
    } // endif

    // Check the West boundary extension to North
    if ( IsNorthExtendWestBoundaryCurved() ){

      // Update ExtendNorth_BndWestSplineInfo[]
      // Check if memory is allocated for ExtendNorth_BndWestSplineInfo
      if(ExtendNorth_BndWestSplineInfo == NULL){ // allocate array
	ExtendNorth_BndWestSplineInfo = new Spline2DInterval_HO [Nghost];
      }

      // Check if the West Spline extension to North is defined such that the normals point outside of the domain 
      if ( ExtendNorth_BndWestSpline.getS(Node[INl][JNu]) > ExtendNorth_BndWestSpline.getS(Node[INl][JNu+1]) ){
	// Update the geometric information 
	IndexShift = JCu + 1;
	for(i=JCu+1; i<=JCu+Nghost; ++i){
	  ExtendNorth_BndWestSplineInfo[i-IndexShift].UpdateInterval(ExtendNorth_BndWestSpline,
								     Node[INl][i],Node[INl][i+1],NumGQPsPerSubinterval);
	}
      } else {
	// Copy the spline
	SplineCopy = ExtendNorth_BndWestSpline;
	// Change the direction of increasing the pathlength
	SplineCopy.Reverse_Spline();

	// Update the geometric information 
	IndexShift = JCu + 1;
	for(i=JCu+1; i<=JCu+Nghost; ++i){
	  ExtendNorth_BndWestSplineInfo[i-IndexShift].UpdateInterval(SplineCopy,
								     Node[INl][i],Node[INl][i+1],NumGQPsPerSubinterval);
	}
      }//endif
    } // endif

    // Check the West boundary extension to South
    if ( IsSouthExtendWestBoundaryCurved() ){

      // Update ExtendSouth_BndWestSplineInfo[]
      // Check if memory is allocated for ExtendSouth_BndWestSplineInfo
      if(ExtendSouth_BndWestSplineInfo == NULL){ // allocate array
	ExtendSouth_BndWestSplineInfo = new Spline2DInterval_HO [Nghost];
      }

      // Check if the West Spline extension to South is defined such that the normals point outside of the domain 
      if ( ExtendSouth_BndWestSpline.getS(Node[INl][JNl-1]) > ExtendSouth_BndWestSpline.getS(Node[INl][JNl]) ){
	// Update the geometric information 
	for(i=0; i<=JCl-1; ++i){
	  ExtendSouth_BndWestSplineInfo[i].UpdateInterval(ExtendSouth_BndWestSpline,
							  Node[INl][i],Node[INl][i+1],NumGQPsPerSubinterval);
	}
      } else {
	// Copy the spline
	SplineCopy = ExtendSouth_BndWestSpline;
	// Change the direction of increasing the pathlength
	SplineCopy.Reverse_Spline();

	// Update the geometric information 
	for(i=0; i<=JCl-1; ++i){
	  ExtendSouth_BndWestSplineInfo[i].UpdateInterval(SplineCopy,Node[INl][i],Node[INl][i+1],NumGQPsPerSubinterval);
	}
      }//endif
    } // endif

  } // endif (CheckExistenceOfCurvedBoundaries)
}

/*!
 * Compute the area of an interior cell that has one or 
 * two edges treated as high-order geometric boundaries.
 */
double Grid2D_Quad_Block_HO::area_CurvedBoundaries(const int &CellIndex, const int &Boundary) const{

  // Obs. The sides of the cell that are not curved are treated as line segments and therefore
  //      the line integral is computed exactly based on the Nodes.
  // The edges are considered in counterclockwise order (i.e. right-hand rule applied)

  if (Gauss_Quad_Curvilinear_Integration && (!Mixed_Curvilinear_Integration) ) {

    // Use the SplineInfo variables to integrate along curved edges.

    switch(Boundary){

    case NORTH_SPLINE:          // Use only the North Spline -> (iCell,jCell)=(CellIndex,JCu)
      return (// cell North side
	      BndNorthSplineInfo[CellIndex].AreaContribution() +
	      // cell West side
	      ZeroLineIntegration(Node[CellIndex  ][JNu],Node[CellIndex  ][JCu]) + 
	      // cell South side
	      ZeroLineIntegration(Node[CellIndex  ][JCu],Node[CellIndex+1][JCu]) + 
	      // cell East side
	      ZeroLineIntegration(Node[CellIndex+1][JCu],Node[CellIndex+1][JNu]) ); 

    case CORNER_NORTH_WEST_SPLINES:     // Use the North Spline and the West Spline
      return (// cell North side
	      BndNorthSplineInfo[CellIndex].AreaContribution() +
	      // cell West side 
	      BndWestSplineInfo[JCu].AreaContribution() +
	      // cell South side
	      ZeroLineIntegration(Node[CellIndex  ][JCu],Node[CellIndex+1][JCu]) +
	      // cell East side
	      ZeroLineIntegration(Node[CellIndex+1][JCu],Node[CellIndex+1][JNu]) ); 

    case CORNER_NORTH_EAST_SPLINES:    // Use the North Spline and the East Spline
      return (// cell North side
	      BndNorthSplineInfo[CellIndex].AreaContribution() +
	      // cell West side
	      ZeroLineIntegration(Node[CellIndex  ][JNu],Node[CellIndex  ][JCu]) + 
	      // cell South side
	      ZeroLineIntegration(Node[CellIndex  ][JCu],Node[CellIndex+1][JCu]) + 
	      // cell East side
	      BndEastSplineInfo[JCu].AreaContribution() ); 
      
    case SOUTH_SPLINE:        // Use only the South Spline
      return (// cell North side
	      ZeroLineIntegration(Node[CellIndex+1][JCl+1],Node[CellIndex  ][JCl+1]) +
	      // cell West side
	      ZeroLineIntegration(Node[CellIndex  ][JCl+1],Node[CellIndex  ][JCl  ]) +
	      // cell South side
	      BndSouthSplineInfo[CellIndex].AreaContribution() +
	      // cell East side
	      ZeroLineIntegration(Node[CellIndex+1][JCl  ],Node[CellIndex+1][JCl+1]) ); 

    case CORNER_SOUTH_WEST_SPLINES:   // Use the South Spline and the West Spline
      return (// cell North side
	      ZeroLineIntegration(Node[CellIndex+1][JCl+1],Node[CellIndex  ][JCl+1]) + 
	      // cell West side
	      BndWestSplineInfo[JCl].AreaContribution() +
	      // cell South side
	      BndSouthSplineInfo[CellIndex].AreaContribution() +
	      // cell East side
	      ZeroLineIntegration(Node[CellIndex+1][JCl  ],Node[CellIndex+1][JCl+1]) ); 

    case CORNER_SOUTH_EAST_SPLINES:   // Use the South Spline and the East Spline
      return ( // cell North side
	      ZeroLineIntegration(Node[CellIndex+1][JCl+1],Node[CellIndex  ][JCl+1]) +
	      // cell West side
	      ZeroLineIntegration(Node[CellIndex  ][JCl+1],Node[CellIndex  ][JCl  ]) +
	      // cell South side
	      BndSouthSplineInfo[CellIndex].AreaContribution() +
	      // cell East side
	      BndEastSplineInfo[JCl].AreaContribution() ); 

    case WEST_SPLINE:       // Use only the West Spline
      return (// cell North side
	      ZeroLineIntegration(Node[ICl+1][CellIndex+1],Node[ICl  ][CellIndex+1]) + 
	      // cell West side
	      BndWestSplineInfo[CellIndex].AreaContribution() +
	      // cell South side
	      ZeroLineIntegration(Node[ICl  ][CellIndex  ],Node[ICl+1][CellIndex  ]) + 
	      // cell East side
	      ZeroLineIntegration(Node[ICl+1][CellIndex  ],Node[ICl+1][CellIndex+1]) ); 

    case EAST_SPLINE:      // Use only the East Spline
      return (// cell North side
	      ZeroLineIntegration(Node[ICu+1][CellIndex+1],Node[ICu  ][CellIndex+1]) + 
	      // cell West side
	      ZeroLineIntegration(Node[ICu  ][CellIndex+1],Node[ICu  ][CellIndex  ]) + 
	      // cell South side
	      ZeroLineIntegration(Node[ICu  ][CellIndex  ],Node[ICu+1][CellIndex  ]) + 
	      // cell East side
	      BndEastSplineInfo[CellIndex].AreaContribution() ); 

    default:
      return 0.0;
    } // endswitch

  } else {

    // Use the boundary spline functions to integrate along curved edges.
    switch(Boundary){

    case NORTH_SPLINE:          // Use only the North Spline -> (iCell,jCell)=(CellIndex,JCu)
      return (// cell North side
	      BndNorthSpline.ZeroOrderIntegration(Node[CellIndex+1][JNu], Node[CellIndex  ][JNu], 15) +
	      // cell West side
	      ZeroLineIntegration(Node[CellIndex  ][JNu].x(),Node[CellIndex  ][JNu].y(),
				  Node[CellIndex  ][JCu].x(),Node[CellIndex  ][JCu].y()) + 
	      // cell South side
	      ZeroLineIntegration(Node[CellIndex  ][JCu].x(),Node[CellIndex  ][JCu].y(),
				  Node[CellIndex+1][JCu].x(),Node[CellIndex+1][JCu].y()) + 
	      // cell East side
	      ZeroLineIntegration(Node[CellIndex+1][JCu].x(),Node[CellIndex+1][JCu].y(),
				  Node[CellIndex+1][JNu].x(),Node[CellIndex+1][JNu].y()) ); 

    case CORNER_NORTH_WEST_SPLINES:     // Use the North Spline and the West Spline
      return (// cell North side
	      BndNorthSpline.ZeroOrderIntegration(Node[CellIndex+1][JNu], Node[CellIndex  ][JNu], 15) +
	      // cell West side 
	      BndWestSpline.ZeroOrderIntegration(Node[CellIndex ][JNu], Node[CellIndex  ][JCu], 15) + 
	      // cell South side
	      ZeroLineIntegration(Node[CellIndex  ][JCu].x(),Node[CellIndex  ][JCu].y(),
				  Node[CellIndex+1][JCu].x(),Node[CellIndex+1][JCu].y()) +
	      // cell East side
	      ZeroLineIntegration(Node[CellIndex+1][JCu].x(),Node[CellIndex+1][JCu].y(),
				  Node[CellIndex+1][JNu].x(),Node[CellIndex+1][JNu].y()) ); 

    case CORNER_NORTH_EAST_SPLINES:    // Use the North Spline and the East Spline
      return (// cell North side
	      BndNorthSpline.ZeroOrderIntegration(Node[CellIndex+1][JNu], Node[CellIndex  ][JNu], 15) + 
	      // cell West side
	      ZeroLineIntegration(Node[CellIndex  ][JNu].x(),Node[CellIndex  ][JNu].y(),
				  Node[CellIndex  ][JCu].x(),Node[CellIndex  ][JCu].y()) + 
	      // cell South side
	      ZeroLineIntegration(Node[CellIndex  ][JCu].x(),Node[CellIndex  ][JCu].y(),
				  Node[CellIndex+1][JCu].x(),Node[CellIndex+1][JCu].y()) + 
	      // cell East side
	      BndEastSpline.ZeroOrderIntegration(Node[CellIndex+1][JCu], Node[CellIndex+1][JNu], 15) ); 
      
    case SOUTH_SPLINE:        // Use only the South Spline
      return (// cell North side
	      ZeroLineIntegration(Node[CellIndex+1][JCl+1].x(),Node[CellIndex+1][JCl+1].y(),
				  Node[CellIndex  ][JCl+1].x(),Node[CellIndex  ][JCl+1].y()) +
	      // cell West side
	      ZeroLineIntegration(Node[CellIndex  ][JCl+1].x(),Node[CellIndex  ][JCl+1].y(),
				  Node[CellIndex  ][JCl  ].x(),Node[CellIndex  ][JCl  ].y()) +
	      // cell South side
	      BndSouthSpline.ZeroOrderIntegration(Node[CellIndex ][JCl], Node[CellIndex+1][JCl], 15) +

	      // cell East side
	      ZeroLineIntegration(Node[CellIndex+1][JCl  ].x(),Node[CellIndex+1][JCl  ].y(),
				  Node[CellIndex+1][JCl+1].x(),Node[CellIndex+1][JCl+1].y()) ); 

    case CORNER_SOUTH_WEST_SPLINES:   // Use the South Spline and the West Spline
      return (// cell North side
	      ZeroLineIntegration(Node[CellIndex+1][JCl+1].x(),Node[CellIndex+1][JCl+1].y(),
				  Node[CellIndex  ][JCl+1].x(),Node[CellIndex  ][JCl+1].y()) + 
	      // cell West side
	      BndWestSpline.ZeroOrderIntegration(Node[CellIndex  ][JCl+1],Node[CellIndex ][JCl], 15) +
	      // cell South side
	      BndSouthSpline.ZeroOrderIntegration(Node[CellIndex ][JCl], Node[CellIndex+1][JCl], 15) + 
	      // cell East side
	      ZeroLineIntegration(Node[CellIndex+1][JCl  ].x(),Node[CellIndex+1][JCl  ].y(),
				  Node[CellIndex+1][JCl+1].x(),Node[CellIndex+1][JCl+1].y()) ); 

    case CORNER_SOUTH_EAST_SPLINES:   // Use the South Spline and the East Spline
      return ( // cell North side
	      ZeroLineIntegration(Node[CellIndex+1][JCl+1].x(),Node[CellIndex+1][JCl+1].y(),
				  Node[CellIndex  ][JCl+1].x(),Node[CellIndex  ][JCl+1].y()) +
	      // cell West side
	      ZeroLineIntegration(Node[CellIndex  ][JCl+1].x(),Node[CellIndex  ][JCl+1].y(),
				  Node[CellIndex  ][JCl  ].x(),Node[CellIndex  ][JCl  ].y()) +
	      // cell South side
	      BndSouthSpline.ZeroOrderIntegration(Node[CellIndex ][JCl], Node[CellIndex+1][JCl], 15) + 
	      // cell East side
	      BndEastSpline.ZeroOrderIntegration(Node[CellIndex+1][JCl], Node[CellIndex+1][JCl+1], 15) ); 

    case WEST_SPLINE:       // Use only the West Spline
      return (// cell North side
	      ZeroLineIntegration(Node[ICl+1][CellIndex+1].x(),Node[ICl+1][CellIndex+1].y(),
				  Node[ICl  ][CellIndex+1].x(),Node[ICl  ][CellIndex+1].y()) + 
	      // cell West side
	      BndWestSpline.ZeroOrderIntegration(Node[ICl][CellIndex+1], Node[ICl][CellIndex], 15) + 
	      // cell South side
	      ZeroLineIntegration(Node[ICl  ][CellIndex  ].x(),Node[ICl  ][CellIndex  ].y(),
				  Node[ICl+1][CellIndex  ].x(),Node[ICl+1][CellIndex  ].y()) + 
	      // cell East side
	      ZeroLineIntegration(Node[ICl+1][CellIndex  ].x(),Node[ICl+1][CellIndex  ].y(),
				  Node[ICl+1][CellIndex+1].x(),Node[ICl+1][CellIndex+1].y()) ); 

    case EAST_SPLINE:      // Use only the East Spline
      return (// cell North side
	      ZeroLineIntegration(Node[ICu+1][CellIndex+1].x(),Node[ICu+1][CellIndex+1].y(),
				  Node[ICu  ][CellIndex+1].x(),Node[ICu  ][CellIndex+1].y()) + 
	      // cell West side
	      ZeroLineIntegration(Node[ICu  ][CellIndex+1].x(),Node[ICu  ][CellIndex+1].y(),
				  Node[ICu  ][CellIndex  ].x(),Node[ICu  ][CellIndex  ].y()) + 
	      // cell South side
	      ZeroLineIntegration(Node[ICu  ][CellIndex  ].x(),Node[ICu  ][CellIndex  ].y(),
				  Node[ICu+1][CellIndex  ].x(),Node[ICu+1][CellIndex  ].y()) + 
	      // cell East side
	      BndEastSpline.ZeroOrderIntegration(Node[ICu+1][CellIndex], Node[ICu+1][CellIndex+1], 15) ); 

    default:
      return 0.0;
    } // endswitch
    
  }// endif
}

/*!
 * Compute the area of a ghost cell that has one or two 
 * edges treated as high-order geometric boundary.
 */
double Grid2D_Quad_Block_HO::area_GhostCell_CurvedBoundaries(const int &CellIndex, const int &Boundary) const{

  // Obs. The sides of the cell that are not curved are treated as line segments and therefore
  //      the line integral is computed exactly based on the Nodes.
  // The edges are considered in counterclockwise order (i.e. right-hand rule applied)

  if (Gauss_Quad_Curvilinear_Integration && (!Mixed_Curvilinear_Integration)) {

    // Use the SplineInfo variables to integrate along curved edges.
    switch(Boundary){
 
    case NORTH_SPLINE:          // Use only the North Spline -> (iCell,jCell)=(CellIndex,JCu+1)
      return (// cell North side
	      ZeroLineIntegration(Node[CellIndex+1][JNu+1],Node[CellIndex  ][JNu+1]) +
	      // cell West side
	      ZeroLineIntegration(Node[CellIndex  ][JNu+1],Node[CellIndex  ][JNu  ]) -   
	      // cell South side
	      BndNorthSplineInfo[CellIndex].AreaContribution() +
	      // cell East side
	      ZeroLineIntegration(Node[CellIndex+1][JNu  ],Node[CellIndex+1][JNu+1]) );  
 
    case SOUTH_SPLINE:        // Use only the South Spline
      return ( // cell North side
	      (- BndSouthSplineInfo[CellIndex].AreaContribution() ) +
	      // cell West side
	      ZeroLineIntegration(Node[CellIndex  ][JNl  ],Node[CellIndex  ][JNl-1]) + 
	      // cell South side
	      ZeroLineIntegration(Node[CellIndex  ][JNl-1],Node[CellIndex+1][JNl-1]) + 
	      // cell East side
	      ZeroLineIntegration(Node[CellIndex+1][JNl-1],Node[CellIndex+1][JNl  ]) ); 

    case WEST_SPLINE:       // Use only the West Spline
      return (// cell North side
	      ZeroLineIntegration(Node[INl  ][CellIndex+1],Node[INl-1][CellIndex+1]) + 
	      // cell West side
	      ZeroLineIntegration(Node[INl-1][CellIndex+1],Node[INl-1][CellIndex  ]) + 
	      // cell South side
	      ZeroLineIntegration(Node[INl-1][CellIndex  ],Node[INl  ][CellIndex  ]) -
	      // cell East side
	      BndWestSplineInfo[CellIndex].AreaContribution() );

    case EAST_SPLINE:      // Use only the East Spline
      return (// cell North side
	      ZeroLineIntegration(Node[INu+1][CellIndex+1],Node[INu  ][CellIndex+1]) - 
	      // cell West side
	      BndEastSplineInfo[CellIndex].AreaContribution() +
	      // cell South side
	      ZeroLineIntegration(Node[INu  ][CellIndex  ],Node[INu+1][CellIndex  ]) + 
	      // cell East side
	      ZeroLineIntegration(Node[INu+1][CellIndex  ],Node[INu+1][CellIndex+1]) ); 

    case EXTEND_W_NORTH_RIGHT_SPLINE: // Use only the extension of North spline to West (right side)
      return (// cell North side
	      ExtendWest_BndNorthSplineInfo[CellIndex].AreaContribution() +
	      // cell West side
	      ZeroLineIntegration(Node[CellIndex  ][JNu  ],Node[CellIndex  ][JNu-1]) +   
	      // cell South side
	      ZeroLineIntegration(Node[CellIndex  ][JNu-1],Node[CellIndex+1][JNu-1]) + 
	      // cell East side
	      ZeroLineIntegration(Node[CellIndex+1][JNu-1],Node[CellIndex+1][JNu  ]) );

    case EXTEND_E_NORTH_RIGHT_SPLINE: // Use only the extension of North spline to East (right side)
      return (// cell North side
	      ZeroLineIntegration(Node[CellIndex+1][JNu+1],Node[CellIndex  ][JNu+1]) + 
	      // cell West side
	      ZeroLineIntegration(Node[CellIndex  ][JNu+1],Node[CellIndex  ][JNu  ]) -   
	      // cell South side
	      ExtendEast_BndNorthSplineInfo[CellIndex-(ICu+1)].AreaContribution() +      
	      // cell East side
	      ZeroLineIntegration(Node[CellIndex+1][JNu  ],Node[CellIndex+1][JNu+1]) );
      
    case EXTEND_W_NORTH_LEFT_SPLINE: // Use only the extension of North spline to West (left side)
      return (// cell North side
	      ZeroLineIntegration(Node[CellIndex+1][JNu+1],Node[CellIndex  ][JNu+1]) + 
	      // cell West side
	      ZeroLineIntegration(Node[CellIndex  ][JNu+1],Node[CellIndex  ][JNu  ]) -   
	      // cell South side
	      ExtendWest_BndNorthSplineInfo[CellIndex].AreaContribution() +      
	      // cell East side
	      ZeroLineIntegration(Node[CellIndex+1][JNu  ],Node[CellIndex+1][JNu+1]) );

    case EXTEND_E_NORTH_LEFT_SPLINE: // Use only the extension of North spline to East (left side)
      return (// cell North side
	      ExtendEast_BndNorthSplineInfo[CellIndex-(ICu+1)].AreaContribution() +
	      // cell West side
	      ZeroLineIntegration(Node[CellIndex  ][JNu  ],Node[CellIndex  ][JNu-1]) +   
	      // cell South side
	      ZeroLineIntegration(Node[CellIndex  ][JNu-1],Node[CellIndex+1][JNu-1]) + 
	      // cell East side
	      ZeroLineIntegration(Node[CellIndex+1][JNu-1],Node[CellIndex+1][JNu  ]) );

    case EXTEND_W_SOUTH_RIGHT_SPLINE: // Use only the extension of South spline to West (right side)
      return (// cell North side
	      -ExtendWest_BndSouthSplineInfo[CellIndex].AreaContribution() +
	      // cell West side
	      ZeroLineIntegration(Node[CellIndex  ][JNl  ],Node[CellIndex  ][JNl-1]) +   
	      // cell South side
	      ZeroLineIntegration(Node[CellIndex  ][JNl-1],Node[CellIndex+1][JNl-1]) + 
	      // cell East side
	      ZeroLineIntegration(Node[CellIndex+1][JNl-1],Node[CellIndex+1][JNl  ]) );

    case EXTEND_E_SOUTH_RIGHT_SPLINE: // Use only the extension of South spline to East (right side)
      return (// cell North side
	      ZeroLineIntegration(Node[CellIndex+1][JNl+1],Node[CellIndex  ][JNl+1]) + 
	      // cell West side
	      ZeroLineIntegration(Node[CellIndex  ][JNl+1],Node[CellIndex  ][JNl  ]) +   
	      // cell South side
	      ExtendEast_BndSouthSplineInfo[CellIndex-(ICu+1)].AreaContribution() +      
	      // cell East side
	      ZeroLineIntegration(Node[CellIndex+1][JNl  ],Node[CellIndex+1][JNl+1]) );
      
    case EXTEND_W_SOUTH_LEFT_SPLINE: // Use only the extension of South spline to West (left side)
      return (// cell North side
	      ZeroLineIntegration(Node[CellIndex+1][JNl+1],Node[CellIndex  ][JNl+1]) + 
	      // cell West side
	      ZeroLineIntegration(Node[CellIndex  ][JNl+1],Node[CellIndex  ][JNl  ]) +
	      // cell South side
	      ExtendWest_BndSouthSplineInfo[CellIndex].AreaContribution() +      
	      // cell East side
	      ZeroLineIntegration(Node[CellIndex+1][JNl  ],Node[CellIndex+1][JNl+1]) );

    case EXTEND_E_SOUTH_LEFT_SPLINE: // Use only the extension of South spline to East (left side)
      return (// cell North side
	      -ExtendEast_BndSouthSplineInfo[CellIndex-(ICu+1)].AreaContribution() +
	      // cell West side
	      ZeroLineIntegration(Node[CellIndex  ][JNl  ],Node[CellIndex  ][JNl-1]) +   
	      // cell South side
	      ZeroLineIntegration(Node[CellIndex  ][JNl-1],Node[CellIndex+1][JNl-1]) + 
	      // cell East side
	      ZeroLineIntegration(Node[CellIndex+1][JNl-1],Node[CellIndex+1][JNl  ]) );

    case EXTEND_N_EAST_RIGHT_SPLINE: // Use only the extension of East spline to North (right side)
      return (// cell North side
	      ZeroLineIntegration(Node[INu  ][CellIndex+1],Node[INu-1][CellIndex+1]) + 
	      // cell West side
	      ZeroLineIntegration(Node[INu-1][CellIndex+1],Node[INu-1][CellIndex  ]) +
	      // cell South side
	      ZeroLineIntegration(Node[INu-1][CellIndex  ],Node[INu  ][CellIndex  ]) + 
	      // cell East side
	      ExtendNorth_BndEastSplineInfo[CellIndex-(JCu+1)].AreaContribution() ); 

    case EXTEND_S_EAST_RIGHT_SPLINE: // Use only the extension of East spline to South (right side)
      return (// cell North side
	      ZeroLineIntegration(Node[INu+1][CellIndex+1],Node[INu  ][CellIndex+1]) -
	      // cell West side
	      ExtendSouth_BndEastSplineInfo[CellIndex].AreaContribution() + 
	      // cell South side
	      ZeroLineIntegration(Node[INu  ][CellIndex  ],Node[INu+1][CellIndex  ]) + 
	      // cell East side
	      ZeroLineIntegration(Node[INu+1][CellIndex  ],Node[INu+1][CellIndex+1]) ); 

    case EXTEND_N_EAST_LEFT_SPLINE: // Use only the extension of East spline to North (left side)
      return (// cell North side
	      ZeroLineIntegration(Node[INu+1][CellIndex+1],Node[INu  ][CellIndex+1]) -
	      // cell West side
	      ExtendNorth_BndEastSplineInfo[CellIndex-(JCu+1)].AreaContribution() + 
	      // cell South side
	      ZeroLineIntegration(Node[INu  ][CellIndex  ],Node[INu+1][CellIndex  ]) + 
	      // cell East side
	      ZeroLineIntegration(Node[INu+1][CellIndex  ],Node[INu+1][CellIndex+1]) ); 

    case EXTEND_S_EAST_LEFT_SPLINE: // Use only the extension of East spline to South (left side)
      return (// cell North side
	      ZeroLineIntegration(Node[INu  ][CellIndex+1],Node[INu-1][CellIndex+1]) + 
	      // cell West side
	      ZeroLineIntegration(Node[INu-1][CellIndex+1],Node[INu-1][CellIndex  ]) +
	      // cell South side
	      ZeroLineIntegration(Node[INu-1][CellIndex  ],Node[INu  ][CellIndex  ]) + 
	      // cell East side
	      ExtendSouth_BndEastSplineInfo[CellIndex].AreaContribution() ); 

    case EXTEND_N_WEST_RIGHT_SPLINE: // Use only the extension of West spline to North (right side)
      return (// cell North side
	      ZeroLineIntegration(Node[INl  ][CellIndex+1],Node[INl-1][CellIndex+1]) + 
	      // cell West side
	      ZeroLineIntegration(Node[INl-1][CellIndex+1],Node[INl-1][CellIndex  ]) +
	      // cell South side
	      ZeroLineIntegration(Node[INl-1][CellIndex  ],Node[INl  ][CellIndex  ]) - 
	      // cell East side
	      ExtendNorth_BndWestSplineInfo[CellIndex-(JCu+1)].AreaContribution() ); 

    case EXTEND_S_WEST_RIGHT_SPLINE: // Use only the extension of West spline to South (right side)
      return (// cell North side
	      ZeroLineIntegration(Node[INl+1][CellIndex+1],Node[INl  ][CellIndex+1]) +
	      // cell West side
	      ExtendSouth_BndWestSplineInfo[CellIndex].AreaContribution() + 
	      // cell South side
	      ZeroLineIntegration(Node[INl  ][CellIndex  ],Node[INl+1][CellIndex  ]) + 
	      // cell East side
	      ZeroLineIntegration(Node[INl+1][CellIndex  ],Node[INl+1][CellIndex+1]) ); 

    case EXTEND_N_WEST_LEFT_SPLINE: // Use only the extension of West spline to North (left side)
      return (// cell North side
	      ZeroLineIntegration(Node[INl+1][CellIndex+1],Node[INl  ][CellIndex+1]) +
	      // cell West side
	      ExtendNorth_BndWestSplineInfo[CellIndex-(JCu+1)].AreaContribution() + 
	      // cell South side
	      ZeroLineIntegration(Node[INl  ][CellIndex  ],Node[INl+1][CellIndex  ]) + 
	      // cell East side
	      ZeroLineIntegration(Node[INl+1][CellIndex  ],Node[INl+1][CellIndex+1]) ); 

    case EXTEND_S_WEST_LEFT_SPLINE: // Use only the extension of West spline to South (left side)
      return (// cell North side
	      ZeroLineIntegration(Node[INl  ][CellIndex+1],Node[INl-1][CellIndex+1]) + 
	      // cell West side
	      ZeroLineIntegration(Node[INl-1][CellIndex+1],Node[INl-1][CellIndex  ]) +
	      // cell South side
	      ZeroLineIntegration(Node[INl-1][CellIndex  ],Node[INl  ][CellIndex  ]) - 
	      // cell East side
	      ExtendSouth_BndWestSplineInfo[CellIndex].AreaContribution() ); 

    case CORNER_NORTH_EXTEND_N_WEST_SPLINES: // Use the North spline and the extension of West spline to North
      return (// cell North side
	      ZeroLineIntegration(Node[INl+1][JNu+1],Node[INl  ][JNu+1]) +
	      // cell West side
	      ExtendNorth_BndWestSplineInfo[0].AreaContribution() -
	      // cell South side
	      BndNorthSplineInfo[ICl].AreaContribution() +
	      // cell East side
	      ZeroLineIntegration(Node[INl+1][JNu  ],Node[INl+1][JNu+1]) ); 

    case CORNER_EXTEND_N_WEST_EXTEND_W_NORTH_SPLINES: // Use the North extension of West and the West extension of North
      return (// cell North side
	      ZeroLineIntegration(Node[INl  ][JNu+1],Node[INl-1][JNu+1]) + 
	      // cell West side
	      ZeroLineIntegration(Node[INl-1][JNu+1],Node[INl-1][JNu  ]) -
	      // cell South side
	      ExtendWest_BndNorthSplineInfo[ICl-1].AreaContribution() -
	      // cell East side
	      ExtendNorth_BndWestSplineInfo[0].AreaContribution() ); 

    case CORNER_EXTEND_W_NORTH_WEST_SPLINES: // Use the West extension of North and the West spline
      return (// cell North side
	      ExtendWest_BndNorthSplineInfo[ICl-1].AreaContribution() +
	      // cell West side
	      ZeroLineIntegration(Node[INl-1][JNu  ],Node[INl-1][JNu-1]) +   
	      // cell South side
	      ZeroLineIntegration(Node[INl-1][JNu-1],Node[INl  ][JNu-1]) -
	      // cell East side
	      BndWestSplineInfo[JCu].AreaContribution() );

    case CORNER_WEST_EXTEND_W_SOUTH_SPLINES: // Use the West spline and the West extension of South
      return (// cell North side
	      ZeroLineIntegration(Node[INl  ][JNl+1],Node[INl-1][JNl+1]) + 
	      // cell West side
	      ZeroLineIntegration(Node[INl-1][JNl+1],Node[INl-1][JNl  ]) +
	      // cell South side
	      ExtendWest_BndSouthSplineInfo[ICl-1].AreaContribution() -      
	      // cell East side
	      BndWestSplineInfo[JCl].AreaContribution() );

    case CORNER_EXTEND_W_SOUTH_EXTEND_S_WEST_SPLINES: // Use the West extension of South and the South extension of West
      return (// cell North side
	      -ExtendWest_BndSouthSplineInfo[ICl-1].AreaContribution() +
	      // cell West side
	      ZeroLineIntegration(Node[INl-1][JNl  ],Node[INl-1][JNl-1]) +   
	      // cell South side
	      ZeroLineIntegration(Node[INl-1][JNl-1],Node[INl  ][JNl-1]) -
	      // cell East side
	      ExtendSouth_BndWestSplineInfo[JCl-1].AreaContribution() );

    case CORNER_EXTEND_S_WEST_SOUTH_SPLINES: // Use the South extension of West and the South spline
      return (// cell North side
	      -BndSouthSplineInfo[ICl].AreaContribution() +
	      // cell West side
	      ExtendSouth_BndWestSplineInfo[JCl-1].AreaContribution() + 
	      // cell South side
	      ZeroLineIntegration(Node[INl  ][JNl-1],Node[INl+1][JNl-1]) + 
	      // cell East side
	      ZeroLineIntegration(Node[INl+1][JNl-1],Node[INl+1][JNl  ]) ); 

    case CORNER_SOUTH_EXTEND_S_EAST_SPLINES: // Use the South spline and the South extension of East
      return (// cell North side
	      -BndSouthSplineInfo[ICu].AreaContribution() +
	      // cell West side
	      ZeroLineIntegration(Node[INu-1][JNl  ],Node[INu-1][JNl-1]) +
	      // cell South side
	      ZeroLineIntegration(Node[INu-1][JNl-1],Node[INu  ][JNl-1]) + 
	      // cell East side
	      ExtendSouth_BndEastSplineInfo[JCl-1].AreaContribution() ); 

    case CORNER_EXTEND_S_EAST_EXTEND_E_SOUTH_SPLINES: // Use the South extension of East and the East extension of South
      return (// cell North side
	      -ExtendEast_BndSouthSplineInfo[0].AreaContribution() -
	      // cell West side
	      ExtendSouth_BndEastSplineInfo[JCl-1].AreaContribution() +
	      // cell South side
	      ZeroLineIntegration(Node[INu  ][JNl-1],Node[INu+1][JNl-1]) + 
	      // cell East side
	      ZeroLineIntegration(Node[INu+1][JNl-1],Node[INu+1][JNl  ]) );

    case CORNER_EXTEND_E_SOUTH_EAST_SPLINES: // Use the East extension of South and the East spline
      return (// cell North side
	      ZeroLineIntegration(Node[INu+1][JNl+1],Node[INu  ][JNl+1]) -
	      // cell West side
	      BndEastSplineInfo[JCl].AreaContribution() +
	      // cell South side
	      ExtendEast_BndSouthSplineInfo[0].AreaContribution() +      
	      // cell East side
	      ZeroLineIntegration(Node[INu+1][JNl  ],Node[INu+1][JNl+1]) );

    case CORNER_EAST_EXTEND_E_NORTH_SPLINES: // Use the East spline and the East extension of North
      return (// cell North side
	      ExtendEast_BndNorthSplineInfo[0].AreaContribution() -
	      // cell West side
	      BndEastSplineInfo[JCu].AreaContribution() +
	      // cell South side
	      ZeroLineIntegration(Node[INu  ][JNu-1],Node[INu+1][JNu-1]) + 
	      // cell East side
	      ZeroLineIntegration(Node[INu+1][JNu-1],Node[INu+1][JNu  ]) );

    case CORNER_EXTEND_E_NORTH_EXTEND_N_EAST_SPLINES: // Use the East extention of North and North extension of East
      return (// cell North side
	      ZeroLineIntegration(Node[INu+1][JNu+1],Node[INu  ][JNu+1]) -
	      // cell West side
	      ExtendNorth_BndEastSplineInfo[0].AreaContribution() -   
	      // cell South side
	      ExtendEast_BndNorthSplineInfo[0].AreaContribution() +      
	      // cell East side
	      ZeroLineIntegration(Node[INu+1][JNu  ],Node[INu+1][JNu+1]) );

    case CORNER_EXTEND_N_EAST_NORTH_SPLINES: // Use the North extension of East and the North spline
      return (// cell North side
	      ZeroLineIntegration(Node[INu  ][JNu+1],Node[INu-1][JNu+1]) + 
	      // cell West side
	      ZeroLineIntegration(Node[INu-1][JNu+1],Node[INu-1][JNu  ]) -
	      // cell South side
	      BndNorthSplineInfo[ICu].AreaContribution() +
	      // cell East side
	      ExtendNorth_BndEastSplineInfo[0].AreaContribution() );

    default:
      return 0.0;
    } // endswitch

  } else {

    // Use the boundary spline functions to integrate along curved edges.
    switch(Boundary){

    case NORTH_SPLINE:          // Use only the North Spline -> (iCell,jCell)=(CellIndex,JCu+1)
      return (// cell North side
	      ZeroLineIntegration(Node[CellIndex+1][JNu+1],Node[CellIndex  ][JNu+1]) +
	      // cell West side
	      ZeroLineIntegration(Node[CellIndex  ][JNu+1],Node[CellIndex  ][JNu  ]) +   
	      // cell South side
	      BndNorthSpline.ZeroOrderIntegration(Node[CellIndex][JNu], Node[CellIndex+1][JNu], 15) + 
	      // cell East side
	      ZeroLineIntegration(Node[CellIndex+1][JNu  ],Node[CellIndex+1][JNu+1]) );  
 
    case SOUTH_SPLINE:        // Use only the South Spline
      return ( // cell North side
	      BndSouthSpline.ZeroOrderIntegration(Node[CellIndex+1][JNl], Node[CellIndex][JNl], 15) +
	      // cell West side
	      ZeroLineIntegration(Node[CellIndex  ][JNl  ],Node[CellIndex  ][JNl-1]) + 
	      // cell South side
	      ZeroLineIntegration(Node[CellIndex  ][JNl-1],Node[CellIndex+1][JNl-1]) + 
	      // cell East side
	      ZeroLineIntegration(Node[CellIndex+1][JNl-1],Node[CellIndex+1][JNl  ]) ); 

    case WEST_SPLINE:       // Use only the West Spline
      return (// cell North side
	      ZeroLineIntegration(Node[INl  ][CellIndex+1],Node[INl-1][CellIndex+1]) + 
	      // cell West side
	      ZeroLineIntegration(Node[INl-1][CellIndex+1],Node[INl-1][CellIndex  ]) + 
	      // cell South side
	      ZeroLineIntegration(Node[INl-1][CellIndex  ],Node[INl  ][CellIndex  ]) + 
	      // cell East side
	      BndWestSpline.ZeroOrderIntegration(Node[INl][CellIndex], Node[INl][CellIndex+1], 15) ); 

    case EAST_SPLINE:      // Use only the East Spline
      return (// cell North side
	      ZeroLineIntegration(Node[INu+1][CellIndex+1],Node[INu  ][CellIndex+1]) + 
	      // cell West side
	      BndEastSpline.ZeroOrderIntegration(Node[INu][CellIndex+1], Node[INu][CellIndex], 15) + 
	      // cell South side
	      ZeroLineIntegration(Node[INu  ][CellIndex  ],Node[INu+1][CellIndex  ]) + 
	      // cell East side
	      ZeroLineIntegration(Node[INu+1][CellIndex  ],Node[INu+1][CellIndex+1]) ); 

    case EXTEND_W_NORTH_RIGHT_SPLINE: // Use only the extension of North spline to West (right side)
      return (// cell North side
	      ExtendWest_BndNorthSpline.ZeroOrderIntegration(Node[CellIndex+1][JNu  ], Node[CellIndex  ][JNu  ], 15) +
	      // cell West side
	      ZeroLineIntegration(Node[CellIndex  ][JNu  ],Node[CellIndex  ][JNu-1]) +   
	      // cell South side
	      ZeroLineIntegration(Node[CellIndex  ][JNu-1],Node[CellIndex+1][JNu-1]) + 
	      // cell East side
	      ZeroLineIntegration(Node[CellIndex+1][JNu-1],Node[CellIndex+1][JNu  ]) );

    case EXTEND_E_NORTH_RIGHT_SPLINE: // Use only the extension of North spline to East (right side)
      return (// cell North side
	      ZeroLineIntegration(Node[CellIndex+1][JNu+1],Node[CellIndex  ][JNu+1]) + 
	      // cell West side
	      ZeroLineIntegration(Node[CellIndex  ][JNu+1],Node[CellIndex  ][JNu  ]) +
	      // cell South side
	      ExtendEast_BndNorthSpline.ZeroOrderIntegration(Node[CellIndex  ][JNu  ], Node[CellIndex+1][JNu  ] ,15) + 
	      // cell East side
	      ZeroLineIntegration(Node[CellIndex+1][JNu  ],Node[CellIndex+1][JNu+1]) );
      
    case EXTEND_W_NORTH_LEFT_SPLINE: // Use only the extension of North spline to West (left side)
      return (// cell North side
	      ZeroLineIntegration(Node[CellIndex+1][JNu+1],Node[CellIndex  ][JNu+1]) + 
	      // cell West side
	      ZeroLineIntegration(Node[CellIndex  ][JNu+1],Node[CellIndex  ][JNu  ]) +  
	      // cell South side
	      ExtendWest_BndNorthSpline.ZeroOrderIntegration(Node[CellIndex  ][JNu  ],Node[CellIndex+1][JNu  ],15) +
	      // cell East side
	      ZeroLineIntegration(Node[CellIndex+1][JNu  ],Node[CellIndex+1][JNu+1]) );

    case EXTEND_E_NORTH_LEFT_SPLINE: // Use only the extension of North spline to East (left side)
      return (// cell North side
	      ExtendEast_BndNorthSpline.ZeroOrderIntegration(Node[CellIndex+1][JNu  ],Node[CellIndex  ][JNu  ],15) +
	      // cell West side
	      ZeroLineIntegration(Node[CellIndex  ][JNu  ],Node[CellIndex  ][JNu-1]) +   
	      // cell South side
	      ZeroLineIntegration(Node[CellIndex  ][JNu-1],Node[CellIndex+1][JNu-1]) + 
	      // cell East side
	      ZeroLineIntegration(Node[CellIndex+1][JNu-1],Node[CellIndex+1][JNu  ]) );

    case EXTEND_W_SOUTH_RIGHT_SPLINE: // Use only the extension of South spline to West (right side)
      return (// cell North side
	      ExtendWest_BndSouthSpline.ZeroOrderIntegration(Node[CellIndex+1][JNl  ],Node[CellIndex  ][JNl  ],15) +
	      // cell West side
	      ZeroLineIntegration(Node[CellIndex  ][JNl  ],Node[CellIndex  ][JNl-1]) +   
	      // cell South side
	      ZeroLineIntegration(Node[CellIndex  ][JNl-1],Node[CellIndex+1][JNl-1]) + 
	      // cell East side
	      ZeroLineIntegration(Node[CellIndex+1][JNl-1],Node[CellIndex+1][JNl  ]) );

    case EXTEND_E_SOUTH_RIGHT_SPLINE: // Use only the extension of South spline to East (right side)
      return (// cell North side
	      ZeroLineIntegration(Node[CellIndex+1][JNl+1],Node[CellIndex  ][JNl+1]) + 
	      // cell West side
	      ZeroLineIntegration(Node[CellIndex  ][JNl+1],Node[CellIndex  ][JNl  ]) +   
	      // cell South side
	      ExtendEast_BndSouthSpline.ZeroOrderIntegration(Node[CellIndex  ][JNl  ],Node[CellIndex+1][JNl  ],15) + 
	      // cell East side
	      ZeroLineIntegration(Node[CellIndex+1][JNl  ],Node[CellIndex+1][JNl+1]) );
      
    case EXTEND_W_SOUTH_LEFT_SPLINE: // Use only the extension of South spline to West (left side)
      return (// cell North side
	      ZeroLineIntegration(Node[CellIndex+1][JNl+1],Node[CellIndex  ][JNl+1]) + 
	      // cell West side
	      ZeroLineIntegration(Node[CellIndex  ][JNl+1],Node[CellIndex  ][JNl  ]) +
	      // cell South side
	      ExtendWest_BndSouthSpline.ZeroOrderIntegration(Node[CellIndex  ][JNl  ], Node[CellIndex+1][JNl  ], 15) +
	      // cell East side
	      ZeroLineIntegration(Node[CellIndex+1][JNl  ],Node[CellIndex+1][JNl+1]) );

    case EXTEND_E_SOUTH_LEFT_SPLINE: // Use only the extension of South spline to East (left side)
      return (// cell North side
	      ExtendEast_BndSouthSpline.ZeroOrderIntegration(Node[CellIndex+1][JNl  ], Node[CellIndex  ][JNl  ],15) + 
	      // cell West side
	      ZeroLineIntegration(Node[CellIndex  ][JNl  ],Node[CellIndex  ][JNl-1]) +   
	      // cell South side
	      ZeroLineIntegration(Node[CellIndex  ][JNl-1],Node[CellIndex+1][JNl-1]) + 
	      // cell East side
	      ZeroLineIntegration(Node[CellIndex+1][JNl-1],Node[CellIndex+1][JNl  ]) );

    case EXTEND_N_EAST_RIGHT_SPLINE: // Use only the extension of East spline to North (right side)
      return (// cell North side
	      ZeroLineIntegration(Node[INu  ][CellIndex+1],Node[INu-1][CellIndex+1]) + 
	      // cell West side
	      ZeroLineIntegration(Node[INu-1][CellIndex+1],Node[INu-1][CellIndex  ]) +
	      // cell South side
	      ZeroLineIntegration(Node[INu-1][CellIndex  ],Node[INu  ][CellIndex  ]) + 
	      // cell East side
	      ExtendNorth_BndEastSpline.ZeroOrderIntegration(Node[INu  ][CellIndex  ], Node[INu  ][CellIndex+1], 15) ); 

    case EXTEND_S_EAST_RIGHT_SPLINE: // Use only the extension of East spline to South (right side)
      return (// cell North side
	      ZeroLineIntegration(Node[INu+1][CellIndex+1],Node[INu  ][CellIndex+1]) +
	      // cell West side
	      ExtendSouth_BndEastSpline.ZeroOrderIntegration(Node[INu  ][CellIndex+1], Node[INu  ][CellIndex  ], 15) + 
	      // cell South side
	      ZeroLineIntegration(Node[INu  ][CellIndex  ],Node[INu+1][CellIndex  ]) + 
	      // cell East side
	      ZeroLineIntegration(Node[INu+1][CellIndex  ],Node[INu+1][CellIndex+1]) ); 

    case EXTEND_N_EAST_LEFT_SPLINE: // Use only the extension of East spline to North (left side)
      return (// cell North side
	      ZeroLineIntegration(Node[INu+1][CellIndex+1],Node[INu  ][CellIndex+1]) +
	      // cell West side
	      ExtendNorth_BndEastSpline.ZeroOrderIntegration(Node[INu  ][CellIndex+1], Node[INu  ][CellIndex  ], 15) + 
	      // cell South side
	      ZeroLineIntegration(Node[INu  ][CellIndex  ],Node[INu+1][CellIndex  ]) + 
	      // cell East side
	      ZeroLineIntegration(Node[INu+1][CellIndex  ],Node[INu+1][CellIndex+1]) ); 

    case EXTEND_S_EAST_LEFT_SPLINE: // Use only the extension of East spline to South (left side)
      return (// cell North side
	      ZeroLineIntegration(Node[INu  ][CellIndex+1],Node[INu-1][CellIndex+1]) + 
	      // cell West side
	      ZeroLineIntegration(Node[INu-1][CellIndex+1],Node[INu-1][CellIndex  ]) +
	      // cell South side
	      ZeroLineIntegration(Node[INu-1][CellIndex  ],Node[INu  ][CellIndex  ]) + 
	      // cell East side
	      ExtendSouth_BndEastSpline.ZeroOrderIntegration(Node[INu  ][CellIndex  ], Node[INu  ][CellIndex+1], 15) ); 

    case EXTEND_N_WEST_RIGHT_SPLINE: // Use only the extension of West spline to North (right side)
      return (// cell North side
	      ZeroLineIntegration(Node[INl  ][CellIndex+1],Node[INl-1][CellIndex+1]) + 
	      // cell West side
	      ZeroLineIntegration(Node[INl-1][CellIndex+1],Node[INl-1][CellIndex  ]) +
	      // cell South side
	      ZeroLineIntegration(Node[INl-1][CellIndex  ],Node[INl  ][CellIndex  ]) +
	      // cell East side
	      ExtendNorth_BndWestSpline.ZeroOrderIntegration(Node[INl  ][CellIndex  ], Node[INl  ][CellIndex+1], 15) ); 

    case EXTEND_S_WEST_RIGHT_SPLINE: // Use only the extension of West spline to South (right side)
      return (// cell North side
	      ZeroLineIntegration(Node[INl+1][CellIndex+1],Node[INl  ][CellIndex+1]) +
	      // cell West side
	      ExtendSouth_BndWestSpline.ZeroOrderIntegration(Node[INl  ][CellIndex+1], Node[INl  ][CellIndex  ], 15) +
	      // cell South side
	      ZeroLineIntegration(Node[INl  ][CellIndex  ],Node[INl+1][CellIndex  ]) + 
	      // cell East side
	      ZeroLineIntegration(Node[INl+1][CellIndex  ],Node[INl+1][CellIndex+1]) ); 

    case EXTEND_N_WEST_LEFT_SPLINE: // Use only the extension of West spline to North (left side)
      return (// cell North side
	      ZeroLineIntegration(Node[INl+1][CellIndex+1],Node[INl  ][CellIndex+1]) +
	      // cell West side
	      ExtendNorth_BndWestSpline.ZeroOrderIntegration(Node[INl  ][CellIndex+1], Node[INl  ][CellIndex  ], 15) +
	      // cell South side
	      ZeroLineIntegration(Node[INl  ][CellIndex  ],Node[INl+1][CellIndex  ]) + 
	      // cell East side
	      ZeroLineIntegration(Node[INl+1][CellIndex  ],Node[INl+1][CellIndex+1]) ); 

    case EXTEND_S_WEST_LEFT_SPLINE: // Use only the extension of West spline to South (left side)
      return (// cell North side
	      ZeroLineIntegration(Node[INl  ][CellIndex+1],Node[INl-1][CellIndex+1]) + 
	      // cell West side
	      ZeroLineIntegration(Node[INl-1][CellIndex+1],Node[INl-1][CellIndex  ]) +
	      // cell South side
	      ZeroLineIntegration(Node[INl-1][CellIndex  ],Node[INl  ][CellIndex  ]) +
	      // cell East side
	      ExtendSouth_BndWestSpline.ZeroOrderIntegration(Node[INl  ][CellIndex  ], Node[INl  ][CellIndex+1], 15) ); 

    case CORNER_NORTH_EXTEND_N_WEST_SPLINES: // Use the North spline and the extension of West spline to North
      return (// cell North side
	      ZeroLineIntegration(Node[INl+1][JNu+1],Node[INl  ][JNu+1]) +
	      // cell West side
	      ExtendNorth_BndWestSpline.ZeroOrderIntegration(Node[INl  ][JNu+1], Node[INl  ][JNu], 15) +
	      // cell South side
	      BndNorthSpline.ZeroOrderIntegration(Node[INl  ][JNu],Node[INl+1][JNu  ], 15) +
	      // cell East side
	      ZeroLineIntegration(Node[INl+1][JNu  ],Node[INl+1][JNu+1]) ); 

    case CORNER_EXTEND_N_WEST_EXTEND_W_NORTH_SPLINES: // Use the North extension of West and the West extension of North
      return (// cell North side
	      ZeroLineIntegration(Node[INl  ][JNu+1],Node[INl-1][JNu+1]) + 
	      // cell West side
	      ZeroLineIntegration(Node[INl-1][JNu+1],Node[INl-1][JNu  ]) +
	      // cell South side
	      ExtendWest_BndNorthSpline.ZeroOrderIntegration(Node[INl-1][JNu  ],Node[INl][JNu  ] ,15) +
	      // cell East side
	      ExtendNorth_BndWestSpline.ZeroOrderIntegration(Node[INl][JNu  ],Node[INl  ][JNu+1], 15) ); 

    case CORNER_EXTEND_W_NORTH_WEST_SPLINES: // Use the West extension of North and the West spline
      return (// cell North side
	      ExtendWest_BndNorthSpline.ZeroOrderIntegration(Node[INl  ][JNu],Node[INl-1][JNu  ], 15) + 
	      // cell West side
	      ZeroLineIntegration(Node[INl-1][JNu  ],Node[INl-1][JNu-1]) +   
	      // cell South side
	      ZeroLineIntegration(Node[INl-1][JNu-1],Node[INl  ][JNu-1]) +
	      // cell East side
	      BndWestSpline.ZeroOrderIntegration(Node[INl  ][JNu-1],Node[INl  ][JNu] , 15) );

    case CORNER_WEST_EXTEND_W_SOUTH_SPLINES: // Use the West spline and the West extension of South
      return (// cell North side
	      ZeroLineIntegration(Node[INl  ][JNl+1],Node[INl-1][JNl+1]) + 
	      // cell West side
	      ZeroLineIntegration(Node[INl-1][JNl+1],Node[INl-1][JNl  ]) +
	      // cell South side
	      ExtendWest_BndSouthSpline.ZeroOrderIntegration(Node[INl-1][JNl  ], Node[INl][JNl  ],15) +      
	      // cell East side
	      BndWestSpline.ZeroOrderIntegration(Node[INl][JNl  ], Node[INl  ][JNl+1],15) );

    case CORNER_EXTEND_W_SOUTH_EXTEND_S_WEST_SPLINES: // Use the West extension of South and the South extension of West
      return (// cell North side
	      ExtendWest_BndSouthSpline.ZeroOrderIntegration(Node[INl  ][JNl], Node[INl-1][JNl  ], 15) +
	      // cell West side
	      ZeroLineIntegration(Node[INl-1][JNl  ],Node[INl-1][JNl-1]) +   
	      // cell South side
	      ZeroLineIntegration(Node[INl-1][JNl-1],Node[INl  ][JNl-1]) +
	      // cell East side
	      ExtendSouth_BndWestSpline.ZeroOrderIntegration(Node[INl  ][JNl-1], Node[INl  ][JNl] , 15) );

    case CORNER_EXTEND_S_WEST_SOUTH_SPLINES: // Use the South extension of West and the South spline
      return (// cell North side
	      BndSouthSpline.ZeroOrderIntegration(Node[INl+1][JNl  ], Node[INl][JNl  ], 15) +
	      // cell West side
	      ExtendSouth_BndWestSpline.ZeroOrderIntegration(Node[INl][JNl  ], Node[INl  ][JNl-1], 15) +
	      // cell South side
	      ZeroLineIntegration(Node[INl  ][JNl-1],Node[INl+1][JNl-1]) + 
	      // cell East side
	      ZeroLineIntegration(Node[INl+1][JNl-1],Node[INl+1][JNl  ]) ); 

    case CORNER_SOUTH_EXTEND_S_EAST_SPLINES: // Use the South spline and the South extension of East
      return (// cell North side
	      BndSouthSpline.ZeroOrderIntegration(Node[INu][JNl  ],Node[INu-1][JNl  ],15) +
	      // cell West side
	      ZeroLineIntegration(Node[INu-1][JNl  ],Node[INu-1][JNl-1]) +
	      // cell South side
	      ZeroLineIntegration(Node[INu-1][JNl-1],Node[INu  ][JNl-1]) + 
	      // cell East side
	      ExtendSouth_BndEastSpline.ZeroOrderIntegration(Node[INu  ][JNl-1], Node[INu][JNl  ], 15) ); 

    case CORNER_EXTEND_S_EAST_EXTEND_E_SOUTH_SPLINES: // Use the South extension of East and the East extension of South
      return (// cell North side
	      ExtendEast_BndSouthSpline.ZeroOrderIntegration(Node[INu+1][JNl  ], Node[INu][JNl  ], 15 ) +
	      // cell West side
	      ExtendSouth_BndEastSpline.ZeroOrderIntegration(Node[INu][JNl  ], Node[INu  ][JNl-1], 15) +
	      // cell South side
	      ZeroLineIntegration(Node[INu  ][JNl-1],Node[INu+1][JNl-1]) + 
	      // cell East side
	      ZeroLineIntegration(Node[INu+1][JNl-1],Node[INu+1][JNl  ]) );

    case CORNER_EXTEND_E_SOUTH_EAST_SPLINES: // Use the East extension of South and the East spline
      return (// cell North side
	      ZeroLineIntegration(Node[INu+1][JNl+1],Node[INu  ][JNl+1]) +
	      // cell West side
	      BndEastSpline.ZeroOrderIntegration(Node[INu  ][JNl+1], Node[INu  ][JNl],15) +
	      // cell South side
	      ExtendEast_BndSouthSpline.ZeroOrderIntegration(Node[INu  ][JNl],Node[INu+1][JNl  ],15) +
	      // cell East side
	      ZeroLineIntegration(Node[INu+1][JNl  ],Node[INu+1][JNl+1]) );

    case CORNER_EAST_EXTEND_E_NORTH_SPLINES: // Use the East spline and the East extension of North
      return (// cell North side
	      ExtendEast_BndNorthSpline.ZeroOrderIntegration(Node[INu+1][JNu  ], Node[INu][JNu  ], 15) +
	      // cell West side
	      BndEastSpline.ZeroOrderIntegration(Node[INu][JNu  ], Node[INu  ][JNu-1], 15) +
	      // cell South side
	      ZeroLineIntegration(Node[INu  ][JNu-1],Node[INu+1][JNu-1]) + 
	      // cell East side
	      ZeroLineIntegration(Node[INu+1][JNu-1],Node[INu+1][JNu  ]) );

    case CORNER_EXTEND_E_NORTH_EXTEND_N_EAST_SPLINES: // Use the East extention of North and North extension of East
      return (// cell North side
	      ZeroLineIntegration(Node[INu+1][JNu+1],Node[INu  ][JNu+1]) +
	      // cell West side
	      ExtendNorth_BndEastSpline.ZeroOrderIntegration(Node[INu  ][JNu+1], Node[INu  ][JNu],15) +
	      // cell South side       
	      ExtendEast_BndNorthSpline.ZeroOrderIntegration(Node[INu  ][JNu], Node[INu+1][JNu  ],15) +      
	      // cell East side
	      ZeroLineIntegration(Node[INu+1][JNu  ],Node[INu+1][JNu+1]) );

    case CORNER_EXTEND_N_EAST_NORTH_SPLINES: // Use the North extension of East and the North spline
      return (// cell North side
	      ZeroLineIntegration(Node[INu  ][JNu+1],Node[INu-1][JNu+1]) + 
	      // cell West side
	      ZeroLineIntegration(Node[INu-1][JNu+1],Node[INu-1][JNu  ]) +
	      // cell South side
	      BndNorthSpline.ZeroOrderIntegration(Node[INu-1][JNu  ], Node[INu][JNu  ], 15) +
	      // cell East side
	      ExtendNorth_BndEastSpline.ZeroOrderIntegration(Node[INu][JNu  ],Node[INu  ][JNu+1], 15) );

    default:
      return 0.0;
    } // endswitch

  }// endif
}

/*!
 * Compute the centroid of an interior cell that has one  
 * or two edges treated as high-order geometric boundaries.
 */
Vector2D Grid2D_Quad_Block_HO::centroid_CurvedBoundaries(const int &CellIndex, const int &Boundary) const{

  // Obs. The sides of the cell that are not curved are treated as line segments and therefore
  //      the line integral is computed exactly based on the Nodes.
  // The edges are considered in counterclockwise order (i.e. right-hand rule applied)

  if (Gauss_Quad_Curvilinear_Integration && (!Mixed_Curvilinear_Integration)) {

    // Use the SplineInfo variables to integrate along curved edges.
    switch(Boundary){

    case NORTH_SPLINE:          // Use only the North Spline -> (iCell,jCell)=(CellIndex,JCu)
      return Vector2D(( // cell North side 
		       BndNorthSplineInfo[CellIndex].XCentroidContribution() +
		       // cell West side
		       PolynomLineIntegration2(Node[CellIndex  ][JNu].x(),Node[CellIndex  ][JNu].y(),
					       Node[CellIndex  ][JCu].x(),Node[CellIndex  ][JCu].y(),
					       0.0, 0.0, 1, 0) +
		       // cell South side
		       PolynomLineIntegration2(Node[CellIndex  ][JCu].x(),Node[CellIndex  ][JCu].y(),
					       Node[CellIndex+1][JCu].x(),Node[CellIndex+1][JCu].y(),
					       0.0, 0.0, 1, 0) +               
		       // cell East side & Division by (OrderX+1)            
		       PolynomLineIntegration2(Node[CellIndex+1][JCu].x(),Node[CellIndex+1][JCu].y(),
					       Node[CellIndex+1][JNu].x(),Node[CellIndex+1][JNu].y(),
					       0.0, 0.0, 1, 0) )*0.5,         
		      // cell North side
		      (BndNorthSplineInfo[CellIndex].YCentroidContribution() +
		       // cell West side  
		       PolynomLineIntegration2(Node[CellIndex  ][JNu].x(),Node[CellIndex  ][JNu].y(),
					       Node[CellIndex  ][JCu].x(),Node[CellIndex  ][JCu].y(),
					       0.0, 0.0, 0, 1) +                            
		       // cell South side 
		       PolynomLineIntegration2(Node[CellIndex  ][JCu].x(),Node[CellIndex  ][JCu].y(),
					       Node[CellIndex+1][JCu].x(),Node[CellIndex+1][JCu].y(),
					       0.0, 0.0, 0, 1) +                            
		       // cell East side & Division by A
		       PolynomLineIntegration2(Node[CellIndex+1][JCu].x(),Node[CellIndex+1][JCu].y(),
					       Node[CellIndex+1][JNu].x(),Node[CellIndex+1][JNu].y(),
					       0.0, 0.0, 0, 1)) ) /Cell[CellIndex][JCu].A;    

    case CORNER_NORTH_WEST_SPLINES:     // Use the North Spline and the West Spline
      return Vector2D( (// cell North side
			BndNorthSplineInfo[CellIndex].XCentroidContribution() +
			// cell West side
			BndWestSplineInfo[JCu].XCentroidContribution() + 
			// cell South side
			PolynomLineIntegration2(Node[CellIndex  ][JCu].x(),Node[CellIndex  ][JCu].y(),
						Node[CellIndex+1][JCu].x(),Node[CellIndex+1][JCu].y(),
						0.0, 0.0, 1, 0) +              
			// cell East side & Division by (OrderX+1)                
			PolynomLineIntegration2(Node[CellIndex+1][JCu].x(),Node[CellIndex+1][JCu].y(),
						Node[CellIndex+1][JNu].x(),Node[CellIndex+1][JNu].y(),
						0.0, 0.0, 1, 0) )*0.5,         
		       // cell North side
		       (BndNorthSplineInfo[CellIndex].YCentroidContribution() +
			// cell West side
			BndWestSplineInfo[JCu].YCentroidContribution() +
			// cell South side 
			PolynomLineIntegration2(Node[CellIndex  ][JCu].x(),Node[CellIndex  ][JCu].y(),
						Node[CellIndex+1][JCu].x(),Node[CellIndex+1][JCu].y(),
						0.0, 0.0, 0, 1) +                             
			// cell East side & Division by A
			PolynomLineIntegration2(Node[CellIndex+1][JCu].x(),Node[CellIndex+1][JCu].y(),
						Node[CellIndex+1][JNu].x(),Node[CellIndex+1][JNu].y(),
						0.0, 0.0, 0, 1) ) ) /Cell[ICl][JCu].A;        

    case CORNER_NORTH_EAST_SPLINES:    // Use the North Spline and the East Spline
      return Vector2D( (// cell North side
			BndNorthSplineInfo[CellIndex].XCentroidContribution() +
			// cell West side
			PolynomLineIntegration2(Node[CellIndex  ][JNu].x(),Node[CellIndex  ][JNu].y(),
						Node[CellIndex  ][JCu].x(),Node[CellIndex  ][JCu].y(),
						0.0, 0.0, 1, 0) +                                        
			// cell South side
			PolynomLineIntegration2(Node[CellIndex  ][JCu].x(),Node[CellIndex  ][JCu].y(),
						Node[CellIndex+1][JCu].x(),Node[CellIndex+1][JCu].y(),
						0.0, 0.0, 1, 0) +                               
			// cell East side & Division by (OrderX+1)         
			BndEastSplineInfo[JCu].XCentroidContribution() )*0.5, 
		       // cell North side
		       (BndNorthSplineInfo[CellIndex].YCentroidContribution() +
			// cell West side
			PolynomLineIntegration2(Node[CellIndex  ][JNu].x(),Node[CellIndex  ][JNu].y(),
						Node[CellIndex  ][JCu].x(),Node[CellIndex  ][JCu].y(),
						0.0, 0.0, 0, 1) +                                       
			// cell South side 
			PolynomLineIntegration2(Node[CellIndex  ][JCu].x(),Node[CellIndex  ][JCu].y(),
						Node[CellIndex+1][JCu].x(),Node[CellIndex+1][JCu].y(),
						0.0, 0.0, 0, 1) +                                       
			// cell East side & Division by A
			BndEastSplineInfo[JCu].YCentroidContribution() ) ) /Cell[ICu][JCu].A; 

    case SOUTH_SPLINE:        // Use only the South Spline
      return Vector2D( (// cell North side
			PolynomLineIntegration2(Node[CellIndex+1][JCl+1].x(),Node[CellIndex+1][JCl+1].y(),
						Node[CellIndex  ][JCl+1].x(),Node[CellIndex  ][JCl+1].y(),
						0.0, 0.0, 1, 0) +                                        
			// cell West side
			PolynomLineIntegration2(Node[CellIndex  ][JCl+1].x(),Node[CellIndex  ][JCl+1].y(),
						Node[CellIndex  ][JCl  ].x(),Node[CellIndex  ][JCl  ].y(),
						0.0, 0.0, 1, 0) +                                        
			// cell South side
			BndSouthSplineInfo[CellIndex].XCentroidContribution() +
			// cell East side & Division by (OrderX+1)
			PolynomLineIntegration2(Node[CellIndex+1][JCl  ].x(),Node[CellIndex+1][JCl  ].y(),
						Node[CellIndex+1][JCl+1].x(),Node[CellIndex+1][JCl+1].y(),
						0.0, 0.0, 1, 0) )*0.5,                
		       // cell North side
		       (PolynomLineIntegration2(Node[CellIndex+1][JCl+1].x(),Node[CellIndex+1][JCl+1].y(),
						Node[CellIndex  ][JCl+1].x(),Node[CellIndex  ][JCl+1].y(),
						0.0, 0.0, 0, 1) +                                        
			// cell West side
			PolynomLineIntegration2(Node[CellIndex  ][JCl+1].x(),Node[CellIndex  ][JCl+1].y(),
						Node[CellIndex  ][JCl  ].x(),Node[CellIndex  ][JCl  ].y(),
						0.0, 0.0, 0, 1) +                                       
			// cell South side 
			BndSouthSplineInfo[CellIndex].YCentroidContribution() +
			// cell East side & Division by A              
			PolynomLineIntegration2(Node[CellIndex+1][JCl  ].x(),Node[CellIndex+1][JCl  ].y(),
						Node[CellIndex+1][JCl+1].x(),Node[CellIndex+1][JCl+1].y(),
						0.0, 0.0, 0, 1) )) /Cell[CellIndex][JCl].A; 

    case CORNER_SOUTH_WEST_SPLINES:   // Use the South Spline and the West Spline
      return Vector2D( (// cell North side
			PolynomLineIntegration2(Node[CellIndex+1][JCl+1].x(),Node[CellIndex+1][JCl+1].y(),
						Node[CellIndex  ][JCl+1].x(),Node[CellIndex  ][JCl+1].y(),
						0.0, 0.0, 1, 0) +                                       
			// cell West side
			BndWestSplineInfo[JCl].XCentroidContribution() +
			// cell South side
			BndSouthSplineInfo[CellIndex].XCentroidContribution() +
			// cell East side & Division by (OrderX+1)
			PolynomLineIntegration2(Node[CellIndex+1][JCl  ].x(),Node[CellIndex+1][JCl  ].y(),
						Node[CellIndex+1][JCl+1].x(),Node[CellIndex+1][JCl+1].y(),
						0.0, 0.0, 1, 0) )*0.5,                
		       // cell North side
		       (PolynomLineIntegration2(Node[CellIndex+1][JCl+1].x(),Node[CellIndex+1][JCl+1].y(),
						Node[CellIndex  ][JCl+1].x(),Node[CellIndex  ][JCl+1].y(),
						0.0, 0.0, 0, 1) +                                       
			// cell West side
			BndWestSplineInfo[JCl].YCentroidContribution() +
			// cell South side
			BndSouthSplineInfo[CellIndex].YCentroidContribution() +
			// cell East side & Division by A
			PolynomLineIntegration2(Node[CellIndex+1][JCl  ].x(),Node[CellIndex+1][JCl  ].y(),
						Node[CellIndex+1][JCl+1].x(),Node[CellIndex+1][JCl+1].y(),
						0.0, 0.0, 0, 1) )) /Cell[ICl][JCl].A; 

    case CORNER_SOUTH_EAST_SPLINES:   // Use the South Spline and the East Spline
      return Vector2D( (// cell North side
			PolynomLineIntegration2(Node[CellIndex+1][JCl+1].x(),Node[CellIndex+1][JCl+1].y(),
						Node[CellIndex  ][JCl+1].x(),Node[CellIndex  ][JCl+1].y(),
						0.0, 0.0, 1, 0) +                                       
			// cell West side 
			PolynomLineIntegration2(Node[CellIndex  ][JCl+1].x(),Node[CellIndex  ][JCl+1].y(),
						Node[CellIndex  ][JCl  ].x(),Node[CellIndex  ][JCl  ].y(),
						0.0, 0.0, 1, 0) +                                       
			// cell South side
			BndSouthSplineInfo[CellIndex].XCentroidContribution() +
			// cell East side & Division by (OrderX+1)       
			BndEastSplineInfo[JCl].XCentroidContribution() )*0.5,   
		       // cell North side
		       (PolynomLineIntegration2(Node[CellIndex+1][JCl+1].x(),Node[CellIndex+1][JCl+1].y(),
						Node[CellIndex  ][JCl+1].x(),Node[CellIndex  ][JCl+1].y(),
						0.0, 0.0, 0, 1) +                                       
			// cell West side
			PolynomLineIntegration2(Node[CellIndex  ][JCl+1].x(),Node[CellIndex  ][JCl+1].y(),
						Node[CellIndex  ][JCl  ].x(),Node[CellIndex  ][JCl  ].y(),
						0.0, 0.0, 0, 1) +                                        
			// cell South side
			BndSouthSplineInfo[CellIndex].YCentroidContribution() +
			// cell East side & Division by A
			BndEastSplineInfo[JCl].YCentroidContribution() )) /Cell[ICu][JCl].A; 

    case WEST_SPLINE:       // Use only the West Spline
      return Vector2D( ( // cell North side
			PolynomLineIntegration2(Node[ICl+1][CellIndex+1].x(),Node[ICl+1][CellIndex+1].y(),
						Node[ICl  ][CellIndex+1].x(),Node[ICl  ][CellIndex+1].y(),
						0.0, 0.0, 1, 0) +                                       
			// cell West side
			BndWestSplineInfo[CellIndex].XCentroidContribution() +
			// cell South side
			PolynomLineIntegration2(Node[ICl  ][CellIndex  ].x(),Node[ICl  ][CellIndex  ].y(),
						Node[ICl+1][CellIndex  ].x(),Node[ICl+1][CellIndex  ].y(),
						0.0, 0.0, 1, 0) +                        
			// cell East side & Division by (OrderX+1)                
			PolynomLineIntegration2(Node[ICl+1][CellIndex  ].x(),Node[ICl+1][CellIndex  ].y(),
						Node[ICl+1][CellIndex+1].x(),Node[ICl+1][CellIndex+1].y(),
						0.0, 0.0, 1, 0) ) *0.5,                  
		       // cell North side
		       (PolynomLineIntegration2(Node[ICl+1][CellIndex+1].x(),Node[ICl+1][CellIndex+1].y(),
						Node[ICl  ][CellIndex+1].x(),Node[ICl  ][CellIndex+1].y(),
						0.0, 0.0, 0, 1) +                                       
			// cell West side 
			BndWestSplineInfo[CellIndex].YCentroidContribution() +
			// cell South side
			PolynomLineIntegration2(Node[ICl  ][CellIndex  ].x(),Node[ICl  ][CellIndex  ].y(),
						Node[ICl+1][CellIndex  ].x(),Node[ICl+1][CellIndex  ].y(),
						0.0, 0.0, 0, 1) +                              
			// cell East side & Division by A         
			PolynomLineIntegration2(Node[ICl+1][CellIndex  ].x(),Node[ICl+1][CellIndex  ].y(),
						Node[ICl+1][CellIndex+1].x(),Node[ICl+1][CellIndex+1].y(),
						0.0, 0.0, 0, 1) )) /Cell[ICl][CellIndex].A;    

    case EAST_SPLINE:      // Use only the East Spline
      return Vector2D( (// cell North side
			PolynomLineIntegration2(Node[ICu+1][CellIndex+1].x(),Node[ICu+1][CellIndex+1].y(),
						Node[ICu  ][CellIndex+1].x(),Node[ICu  ][CellIndex+1].y(),
						0.0, 0.0, 1, 0) +                                        
			// cell West side
			PolynomLineIntegration2(Node[ICu  ][CellIndex+1].x(),Node[ICu  ][CellIndex+1].y(),
						Node[ICu  ][CellIndex  ].x(),Node[ICu  ][CellIndex  ].y(),
						0.0, 0.0, 1, 0) +                                        
			// cell South side 
			PolynomLineIntegration2(Node[ICu  ][CellIndex  ].x(),Node[ICu  ][CellIndex  ].y(),
						Node[ICu+1][CellIndex  ].x(),Node[ICu+1][CellIndex  ].y(),
						0.0, 0.0, 1, 0) +                                
			// cell East side & Division by (OrderX+1)        
			BndEastSplineInfo[CellIndex].XCentroidContribution() ) *0.5, 
		       // cell North side
		       (PolynomLineIntegration2(Node[ICu+1][CellIndex+1].x(),Node[ICu+1][CellIndex+1].y(),
						Node[ICu  ][CellIndex+1].x(),Node[ICu  ][CellIndex+1].y(),
						0.0, 0.0, 0, 1) +                                        
			// cell West side
			PolynomLineIntegration2(Node[ICu  ][CellIndex+1].x(),Node[ICu  ][CellIndex+1].y(),
						Node[ICu  ][CellIndex  ].x(),Node[ICu  ][CellIndex  ].y(),
						0.0, 0.0, 0, 1) +                                        
			// cell South side
			PolynomLineIntegration2(Node[ICu  ][CellIndex  ].x(),Node[ICu  ][CellIndex  ].y(),
						Node[ICu+1][CellIndex  ].x(),Node[ICu+1][CellIndex  ].y(),
						0.0, 0.0, 0, 1) +                                        
			// cell East side & Division by A
			BndEastSplineInfo[CellIndex].YCentroidContribution() )) /Cell[ICu][CellIndex].A; 

    default:
      return Vector2D(0.0);
    } // endswitch

  } else {

    // Use the boundary spline functions to integrate along curved edges.
    switch(Boundary){

    case NORTH_SPLINE:          // Use only the North Spline -> (iCell,jCell)=(CellIndex,JCu)
      return Vector2D(( // cell North side 
		       BndNorthSpline.PolynomOrderIntegration(Node[CellIndex+1][JNu], Node[CellIndex  ][JNu], 
							      Vector2D(0.0,0.0),15,1,0) +
		       // cell West side
		       PolynomLineIntegration2(Node[CellIndex  ][JNu].x(),Node[CellIndex  ][JNu].y(),
					       Node[CellIndex  ][JCu].x(),Node[CellIndex  ][JCu].y(),
					       0.0, 0.0, 1, 0) +
		       // cell South side
		       PolynomLineIntegration2(Node[CellIndex  ][JCu].x(),Node[CellIndex  ][JCu].y(),
					       Node[CellIndex+1][JCu].x(),Node[CellIndex+1][JCu].y(),
					       0.0, 0.0, 1, 0) +               
		       // cell East side & Division by (OrderX+1)            
		       PolynomLineIntegration2(Node[CellIndex+1][JCu].x(),Node[CellIndex+1][JCu].y(),
					       Node[CellIndex+1][JNu].x(),Node[CellIndex+1][JNu].y(),
					       0.0, 0.0, 1, 0) )*0.5,         
		      // cell North side
		      (BndNorthSpline.PolynomOrderIntegration(Node[CellIndex+1][JNu], Node[CellIndex  ][JNu],
							      Vector2D(0.0,0.0), 15, 0, 1) +
		       // cell West side  
		       PolynomLineIntegration2(Node[CellIndex  ][JNu].x(),Node[CellIndex  ][JNu].y(),
					       Node[CellIndex  ][JCu].x(),Node[CellIndex  ][JCu].y(),
					       0.0, 0.0, 0, 1) +                            
		       // cell South side 
		       PolynomLineIntegration2(Node[CellIndex  ][JCu].x(),Node[CellIndex  ][JCu].y(),
					       Node[CellIndex+1][JCu].x(),Node[CellIndex+1][JCu].y(),
					       0.0, 0.0, 0, 1) +                            
		       // cell East side & Division by A
		       PolynomLineIntegration2(Node[CellIndex+1][JCu].x(),Node[CellIndex+1][JCu].y(),
					       Node[CellIndex+1][JNu].x(),Node[CellIndex+1][JNu].y(),
					       0.0, 0.0, 0, 1)) ) /Cell[CellIndex][JCu].A;    
    

    case CORNER_NORTH_WEST_SPLINES:     // Use the North Spline and the West Spline
      return Vector2D( (// cell North side
			BndNorthSpline.PolynomOrderIntegration(Node[CellIndex+1][JNu], Node[CellIndex  ][JNu],
							       Vector2D(0.0,0.0),15,1,0) +       
			// cell West side
			BndWestSpline.PolynomOrderIntegration(Node[CellIndex ][JNu], Node[CellIndex  ][JCu],
							      Vector2D(0.0,0.0),15,1,0) +        
			// cell South side
			PolynomLineIntegration2(Node[CellIndex  ][JCu].x(),Node[CellIndex  ][JCu].y(),
						Node[CellIndex+1][JCu].x(),Node[CellIndex+1][JCu].y(),
						0.0, 0.0, 1, 0) +              
			// cell East side & Division by (OrderX+1)                
			PolynomLineIntegration2(Node[CellIndex+1][JCu].x(),Node[CellIndex+1][JCu].y(),
						Node[CellIndex+1][JNu].x(),Node[CellIndex+1][JNu].y(),
						0.0, 0.0, 1, 0) )*0.5,         
		       // cell North side
		       (BndNorthSpline.PolynomOrderIntegration(Node[CellIndex+1][JNu], Node[CellIndex  ][JNu],
							       Vector2D(0.0,0.0), 15, 0, 1) +    
			// cell West side
			BndWestSpline.PolynomOrderIntegration(Node[CellIndex ][JNu], Node[CellIndex  ][JCu],
							      Vector2D(0.0,0.0), 15, 0, 1) +    
			// cell South side 
			PolynomLineIntegration2(Node[CellIndex  ][JCu].x(),Node[CellIndex  ][JCu].y(),
						Node[CellIndex+1][JCu].x(),Node[CellIndex+1][JCu].y(),
						0.0, 0.0, 0, 1) +                             
			// cell East side & Division by A
			PolynomLineIntegration2(Node[CellIndex+1][JCu].x(),Node[CellIndex+1][JCu].y(),
						Node[CellIndex+1][JNu].x(),Node[CellIndex+1][JNu].y(),
						0.0, 0.0, 0, 1) ) ) /Cell[ICl][JCu].A;        

    case CORNER_NORTH_EAST_SPLINES:    // Use the North Spline and the East Spline
      return Vector2D( (// cell North side
			BndNorthSpline.PolynomOrderIntegration(Node[CellIndex+1][JNu], Node[CellIndex  ][JNu],
							       Vector2D(0.0,0.0),15,1,0) +                 
			// cell West side
			PolynomLineIntegration2(Node[CellIndex  ][JNu].x(),Node[CellIndex  ][JNu].y(),
						Node[CellIndex  ][JCu].x(),Node[CellIndex  ][JCu].y(),
						0.0, 0.0, 1, 0) +                                        
			// cell South side
			PolynomLineIntegration2(Node[CellIndex  ][JCu].x(),Node[CellIndex  ][JCu].y(),
						Node[CellIndex+1][JCu].x(),Node[CellIndex+1][JCu].y(),
						0.0, 0.0, 1, 0) +                               
			// cell East side & Division by (OrderX+1)         
			BndEastSpline.PolynomOrderIntegration(Node[CellIndex+1][JCu], Node[CellIndex+1][JNu],
							      Vector2D(0.0,0.0),15,1,0) )*0.5, 
		       // cell North side
		       (BndNorthSpline.PolynomOrderIntegration(Node[CellIndex+1][JNu], Node[CellIndex  ][JNu],
							       Vector2D(0.0,0.0),15,0,1) +              
			// cell West side
			PolynomLineIntegration2(Node[CellIndex  ][JNu].x(),Node[CellIndex  ][JNu].y(),
						Node[CellIndex  ][JCu].x(),Node[CellIndex  ][JCu].y(),
						0.0, 0.0, 0, 1) +                                       
			// cell South side 
			PolynomLineIntegration2(Node[CellIndex  ][JCu].x(),Node[CellIndex  ][JCu].y(),
						Node[CellIndex+1][JCu].x(),Node[CellIndex+1][JCu].y(),
						0.0, 0.0, 0, 1) +                                       
			// cell East side & Division by A
			BndEastSpline.PolynomOrderIntegration(Node[CellIndex+1][JCu], Node[CellIndex+1][JNu],
							      Vector2D(0.0,0.0),15,0,1) ) ) /Cell[ICu][JCu].A; 
    case SOUTH_SPLINE:        // Use only the South Spline
      return Vector2D( (// cell North side
			PolynomLineIntegration2(Node[CellIndex+1][JCl+1].x(),Node[CellIndex+1][JCl+1].y(),
						Node[CellIndex  ][JCl+1].x(),Node[CellIndex  ][JCl+1].y(),
						0.0, 0.0, 1, 0) +                                        
			// cell West side
			PolynomLineIntegration2(Node[CellIndex  ][JCl+1].x(),Node[CellIndex  ][JCl+1].y(),
						Node[CellIndex  ][JCl  ].x(),Node[CellIndex  ][JCl  ].y(),
						0.0, 0.0, 1, 0) +                                        
			// cell South side
			BndSouthSpline.PolynomOrderIntegration(Node[CellIndex ][JCl], Node[CellIndex+1][JCl],
							       Vector2D(0.0,0.0),15,1,0) +              
			// cell East side & Division by (OrderX+1)
			PolynomLineIntegration2(Node[CellIndex+1][JCl  ].x(),Node[CellIndex+1][JCl  ].y(),
						Node[CellIndex+1][JCl+1].x(),Node[CellIndex+1][JCl+1].y(),
						0.0, 0.0, 1, 0) )*0.5,                
		       // cell North side
		       (PolynomLineIntegration2(Node[CellIndex+1][JCl+1].x(),Node[CellIndex+1][JCl+1].y(),
						Node[CellIndex  ][JCl+1].x(),Node[CellIndex  ][JCl+1].y(),
						0.0, 0.0, 0, 1) +                                        
			// cell West side
			PolynomLineIntegration2(Node[CellIndex  ][JCl+1].x(),Node[CellIndex  ][JCl+1].y(),
						Node[CellIndex  ][JCl  ].x(),Node[CellIndex  ][JCl  ].y(),
						0.0, 0.0, 0, 1) +                                       
			// cell South side 
			BndSouthSpline.PolynomOrderIntegration(Node[CellIndex ][JCl], Node[CellIndex+1][JCl],
							       Vector2D(0.0,0.0),15,0,1) + 
			// cell East side & Division by A              
			PolynomLineIntegration2(Node[CellIndex+1][JCl  ].x(),Node[CellIndex+1][JCl  ].y(),
						Node[CellIndex+1][JCl+1].x(),Node[CellIndex+1][JCl+1].y(),
						0.0, 0.0, 0, 1) )) /Cell[CellIndex][JCl].A; 

    case CORNER_SOUTH_WEST_SPLINES:   // Use the South Spline and the West Spline
      return Vector2D( (// cell North side
			PolynomLineIntegration2(Node[CellIndex+1][JCl+1].x(),Node[CellIndex+1][JCl+1].y(),
						Node[CellIndex  ][JCl+1].x(),Node[CellIndex  ][JCl+1].y(),
						0.0, 0.0, 1, 0) +                                       
			// cell West side
			BndWestSpline.PolynomOrderIntegration(Node[CellIndex][JCl+1], Node[CellIndex][JCl],
							      Vector2D(0.0,0.0),15,1,0) + 
			// cell South side
			BndSouthSpline.PolynomOrderIntegration(Node[CellIndex ][JCl], Node[CellIndex+1][JCl],
							       Vector2D(0.0,0.0),15,1,0) +              
			// cell East side & Division by (OrderX+1)
			PolynomLineIntegration2(Node[CellIndex+1][JCl  ].x(),Node[CellIndex+1][JCl  ].y(),
						Node[CellIndex+1][JCl+1].x(),Node[CellIndex+1][JCl+1].y(),
						0.0, 0.0, 1, 0) )*0.5,                
		       // cell North side
		       (PolynomLineIntegration2(Node[CellIndex+1][JCl+1].x(),Node[CellIndex+1][JCl+1].y(),
						Node[CellIndex  ][JCl+1].x(),Node[CellIndex  ][JCl+1].y(),
						0.0, 0.0, 0, 1) +                                       
			// cell West side
			BndWestSpline.PolynomOrderIntegration(Node[CellIndex][JCl+1], Node[CellIndex][JCl],
							      Vector2D(0.0,0.0),15,0,1) + 
			// cell South side
			BndSouthSpline.PolynomOrderIntegration(Node[CellIndex ][JCl], Node[CellIndex+1][JCl],
							       Vector2D(0.0,0.0),15,0,1) +              
			// cell East side & Division by A
			PolynomLineIntegration2(Node[CellIndex+1][JCl  ].x(),Node[CellIndex+1][JCl  ].y(),
						Node[CellIndex+1][JCl+1].x(),Node[CellIndex+1][JCl+1].y(),
						0.0, 0.0, 0, 1) )) /Cell[ICl][JCl].A; 

    case CORNER_SOUTH_EAST_SPLINES:   // Use the South Spline and the East Spline
      return Vector2D( (// cell North side
			PolynomLineIntegration2(Node[CellIndex+1][JCl+1].x(),Node[CellIndex+1][JCl+1].y(),
						Node[CellIndex  ][JCl+1].x(),Node[CellIndex  ][JCl+1].y(),
						0.0, 0.0, 1, 0) +                                       
			// cell West side 
			PolynomLineIntegration2(Node[CellIndex  ][JCl+1].x(),Node[CellIndex  ][JCl+1].y(),
						Node[CellIndex  ][JCl  ].x(),Node[CellIndex  ][JCl  ].y(),
						0.0, 0.0, 1, 0) +                                       
			// cell South side
			BndSouthSpline.PolynomOrderIntegration(Node[CellIndex ][JCl], Node[CellIndex+1][JCl],
							       Vector2D(0.0,0.0),15,1,0) +       
			// cell East side & Division by (OrderX+1)       
			BndEastSpline.PolynomOrderIntegration(Node[CellIndex+1][JCl], Node[CellIndex+1][JCl+1],
							      Vector2D(0.0,0.0),15,1,0) )*0.5,   
		       // cell North side
		       (PolynomLineIntegration2(Node[CellIndex+1][JCl+1].x(),Node[CellIndex+1][JCl+1].y(),
						Node[CellIndex  ][JCl+1].x(),Node[CellIndex  ][JCl+1].y(),
						0.0, 0.0, 0, 1) +                                       
			// cell West side
			PolynomLineIntegration2(Node[CellIndex  ][JCl+1].x(),Node[CellIndex  ][JCl+1].y(),
						Node[CellIndex  ][JCl  ].x(),Node[CellIndex  ][JCl  ].y(),
						0.0, 0.0, 0, 1) +                                        
			// cell South side
			BndSouthSpline.PolynomOrderIntegration(Node[CellIndex ][JCl], Node[CellIndex+1][JCl],
							       Vector2D(0.0,0.0),15,0,1) +              
			// cell East side & Division by A
			BndEastSpline.PolynomOrderIntegration(Node[CellIndex+1][JCl], Node[CellIndex+1][JCl+1],
							      Vector2D(0.0,0.0),15,0,1) )) /Cell[ICu][JCl].A; 

    case WEST_SPLINE:       // Use only the West Spline
      return Vector2D( ( // cell North side
			PolynomLineIntegration2(Node[ICl+1][CellIndex+1].x(),Node[ICl+1][CellIndex+1].y(),
						Node[ICl  ][CellIndex+1].x(),Node[ICl  ][CellIndex+1].y(),
						0.0, 0.0, 1, 0) +                                       
			// cell West side
			BndWestSpline.PolynomOrderIntegration(Node[ICl][CellIndex+1], Node[ICl][CellIndex],
							      Vector2D(0.0,0.0),15,1,0) +              
			// cell South side
			PolynomLineIntegration2(Node[ICl  ][CellIndex  ].x(),Node[ICl  ][CellIndex  ].y(),
						Node[ICl+1][CellIndex  ].x(),Node[ICl+1][CellIndex  ].y(),
						0.0, 0.0, 1, 0) +                        
			// cell East side & Division by (OrderX+1)                
			PolynomLineIntegration2(Node[ICl+1][CellIndex  ].x(),Node[ICl+1][CellIndex  ].y(),
						Node[ICl+1][CellIndex+1].x(),Node[ICl+1][CellIndex+1].y(),
						0.0, 0.0, 1, 0) ) *0.5,                  
		       // cell North side
		       (PolynomLineIntegration2(Node[ICl+1][CellIndex+1].x(),Node[ICl+1][CellIndex+1].y(),
						Node[ICl  ][CellIndex+1].x(),Node[ICl  ][CellIndex+1].y(),
						0.0, 0.0, 0, 1) +                                       
			// cell West side 
			BndWestSpline.PolynomOrderIntegration(Node[ICl][CellIndex+1], Node[ICl][CellIndex],
							      Vector2D(0.0,0.0),15,0,1) +              
			// cell South side
			PolynomLineIntegration2(Node[ICl  ][CellIndex  ].x(),Node[ICl  ][CellIndex  ].y(),
						Node[ICl+1][CellIndex  ].x(),Node[ICl+1][CellIndex  ].y(),
						0.0, 0.0, 0, 1) +                              
			// cell East side & Division by A         
			PolynomLineIntegration2(Node[ICl+1][CellIndex  ].x(),Node[ICl+1][CellIndex  ].y(),
						Node[ICl+1][CellIndex+1].x(),Node[ICl+1][CellIndex+1].y(),
						0.0, 0.0, 0, 1) )) /Cell[ICl][CellIndex].A;    

    case EAST_SPLINE:      // Use only the East Spline
      return Vector2D( (// cell North side
			PolynomLineIntegration2(Node[ICu+1][CellIndex+1].x(),Node[ICu+1][CellIndex+1].y(),
						Node[ICu  ][CellIndex+1].x(),Node[ICu  ][CellIndex+1].y(),
						0.0, 0.0, 1, 0) +                                        
			// cell West side
			PolynomLineIntegration2(Node[ICu  ][CellIndex+1].x(),Node[ICu  ][CellIndex+1].y(),
						Node[ICu  ][CellIndex  ].x(),Node[ICu  ][CellIndex  ].y(),
						0.0, 0.0, 1, 0) +                                        
			// cell South side 
			PolynomLineIntegration2(Node[ICu  ][CellIndex  ].x(),Node[ICu  ][CellIndex  ].y(),
						Node[ICu+1][CellIndex  ].x(),Node[ICu+1][CellIndex  ].y(),
						0.0, 0.0, 1, 0) +                                
			// cell East side & Division by (OrderX+1)        
			BndEastSpline.PolynomOrderIntegration(Node[ICu+1][CellIndex], Node[ICu+1][CellIndex+1],
							      Vector2D(0.0,0.0),15,1,0) ) *0.5, 
		       // cell North side
		       (PolynomLineIntegration2(Node[ICu+1][CellIndex+1].x(),Node[ICu+1][CellIndex+1].y(),
						Node[ICu  ][CellIndex+1].x(),Node[ICu  ][CellIndex+1].y(),
						0.0, 0.0, 0, 1) +                                        
			// cell West side
			PolynomLineIntegration2(Node[ICu  ][CellIndex+1].x(),Node[ICu  ][CellIndex+1].y(),
						Node[ICu  ][CellIndex  ].x(),Node[ICu  ][CellIndex  ].y(),
						0.0, 0.0, 0, 1) +                                        
			// cell South side
			PolynomLineIntegration2(Node[ICu  ][CellIndex  ].x(),Node[ICu  ][CellIndex  ].y(),
						Node[ICu+1][CellIndex  ].x(),Node[ICu+1][CellIndex  ].y(),
						0.0, 0.0, 0, 1) +                                        
			// cell East side & Division by A
			BndEastSpline.PolynomOrderIntegration(Node[ICu+1][CellIndex], Node[ICu+1][CellIndex+1],
							      Vector2D(0.0,0.0),15,0,1) )) /Cell[ICu][CellIndex].A; 

    default:
      return Vector2D(0.0);
    } // endswitch

  }// endif
}

/*!
 * Compute the centroid of a ghost cell that has one  
 * edge treated as high-order geometric boundary.
 */
Vector2D Grid2D_Quad_Block_HO::centroid_GhostCell_CurvedBoundaries(const int &CellIndex, const int &Boundary) const{

  // Obs. The sides of the cell that are not curved are treated as line segments and therefore
  //      the line integral is computed exactly based on the Nodes.
  // The edges are considered in counterclockwise order (i.e. right-hand rule applied)
  // ***! This subroutine is for the first row of ghost cells and therefore only one boundary can be curved at a time

  if (Gauss_Quad_Curvilinear_Integration && (!Mixed_Curvilinear_Integration)) {

    // Use the SplineInfo variables to integrate along curved edges.

    switch(Boundary){

    case NORTH_SPLINE:          // Use only the North Spline -> (iCell,jCell)=(CellIndex,JCu+1)
      return Vector2D(( PolynomLineIntegration2(Node[CellIndex+1][JNu+1],Node[CellIndex  ][JNu+1],
						0.0, 0.0, 1, 0) +                           // cell North side
			PolynomLineIntegration2(Node[CellIndex  ][JNu+1],Node[CellIndex  ][JNu  ],
						0.0, 0.0, 1, 0) -                           // cell West side
			BndNorthSplineInfo[CellIndex].XCentroidContribution() + // cell South side
			PolynomLineIntegration2(Node[CellIndex+1][JNu  ],Node[CellIndex+1][JNu+1],
						0.0, 0.0, 1, 0) )*0.5,         // cell East side & Division by (OrderX+1)
		      (PolynomLineIntegration2(Node[CellIndex+1][JNu+1],Node[CellIndex  ][JNu+1],
					       0.0, 0.0, 0, 1) +                             // cell North side
		       PolynomLineIntegration2(Node[CellIndex  ][JNu+1],Node[CellIndex  ][JNu  ],
					       0.0, 0.0, 0, 1) -                             // cell West side
		       BndNorthSplineInfo[CellIndex].YCentroidContribution() + // cell South side
		       PolynomLineIntegration2(Node[CellIndex+1][JNu  ],Node[CellIndex+1][JNu+1],
					       0.0, 0.0, 0, 1)) ) /Cell[CellIndex][JCu+1].A;    // cell East side & Division by A
    
    case SOUTH_SPLINE:        // Use only the South Spline
      return Vector2D( (-BndSouthSplineInfo[CellIndex].XCentroidContribution() +              // cell North side
			PolynomLineIntegration2(Node[CellIndex  ][JNl  ],Node[CellIndex  ][JNl-1],
						0.0, 0.0, 1, 0) +                                        // cell West side
			PolynomLineIntegration2(Node[CellIndex  ][JNl-1],Node[CellIndex+1][JNl-1],
						0.0, 0.0, 1, 0) +                                        // cell South side
			PolynomLineIntegration2(Node[CellIndex+1][JNl-1],Node[CellIndex+1][JNl  ],
						0.0, 0.0, 1, 0) )*0.5,                // cell East side & Division by (OrderX+1)
		       (-BndSouthSplineInfo[CellIndex].YCentroidContribution() +              // cell North side
			PolynomLineIntegration2(Node[CellIndex  ][JNl  ],Node[CellIndex  ][JNl-1],
						0.0, 0.0, 0, 1) +                                        // cell West side
			PolynomLineIntegration2(Node[CellIndex  ][JNl-1],Node[CellIndex+1][JNl-1],
						0.0, 0.0, 0, 1) +                                        // cell South side
			PolynomLineIntegration2(Node[CellIndex+1][JNl-1],Node[CellIndex+1][JNl  ],
						0.0, 0.0, 0, 1) )) /Cell[CellIndex][JCl-1].A; // cell East side & Division by A  

    case WEST_SPLINE:       // Use only the West Spline
      return Vector2D( (PolynomLineIntegration2(Node[INl  ][CellIndex+1],Node[INl-1][CellIndex+1],
						0.0, 0.0, 1, 0) +                                        // cell North side
			PolynomLineIntegration2(Node[INl-1][CellIndex+1],Node[INl-1][CellIndex  ],
						0.0, 0.0, 1, 0) +                                        // cell West side
			PolynomLineIntegration2(Node[INl-1][CellIndex  ],Node[INl  ][CellIndex  ],
						0.0, 0.0, 1, 0) -                                        // cell South side
			BndWestSplineInfo[CellIndex].XCentroidContribution() ) *0.5, // cell East side & Division by (OrderX+1)
		       (PolynomLineIntegration2(Node[INl  ][CellIndex+1],Node[INl-1][CellIndex+1],
						0.0, 0.0, 0, 1) +                                        // cell North side
			PolynomLineIntegration2(Node[INl-1][CellIndex+1],Node[INl-1][CellIndex  ],
						0.0, 0.0, 0, 1) +                                        // cell West side
			PolynomLineIntegration2(Node[INl-1][CellIndex  ],Node[INl  ][CellIndex  ],
						0.0, 0.0, 0, 1) -                                        // cell South side
			BndWestSplineInfo[CellIndex].YCentroidContribution() )) /Cell[ICl-1][CellIndex].A;
      // cell East side & Division by A

    case EAST_SPLINE:      // Use only the East Spline
      return Vector2D( (PolynomLineIntegration2(Node[INu+1][CellIndex+1],Node[INu  ][CellIndex+1],
						0.0, 0.0, 1, 0) -                                         // cell North side
			BndEastSplineInfo[CellIndex].XCentroidContribution() +                // cell West side
			PolynomLineIntegration2(Node[INu  ][CellIndex  ],Node[INu+1][CellIndex  ],
						0.0, 0.0, 1, 0) +                                         // cell South side
			PolynomLineIntegration2(Node[INu+1][CellIndex  ],Node[INu+1][CellIndex+1],
						0.0, 0.0, 1, 0) ) *0.5, // cell East side & Division by (OrderX+1)
		       (PolynomLineIntegration2(Node[INu+1][CellIndex+1],Node[INu  ][CellIndex+1],
						0.0, 0.0, 0, 1) -                                         // cell North side
			BndEastSplineInfo[CellIndex].YCentroidContribution() +                // cell West side
			PolynomLineIntegration2(Node[INu  ][CellIndex  ],Node[INu+1][CellIndex  ],
						0.0, 0.0, 0, 1) +                                         // cell South side
			PolynomLineIntegration2(Node[INu+1][CellIndex  ],Node[INu+1][CellIndex+1],
						0.0, 0.0, 0, 1) )) /Cell[ICu+1][CellIndex].A; // cell East side & Division by A

    case EXTEND_W_NORTH_RIGHT_SPLINE: // Use only the extension of North spline to West (right side)
      return Vector2D( (ExtendWest_BndNorthSplineInfo[CellIndex].XCentroidContribution() +              // cell North side
			PolynomLineIntegration2(Node[CellIndex  ][JNu  ],Node[CellIndex  ][JNu-1],
						0.0, 0.0, 1, 0) +                                        // cell West side
			PolynomLineIntegration2(Node[CellIndex  ][JNu-1],Node[CellIndex+1][JNu-1],
						0.0, 0.0, 1, 0) +                                        // cell South side
			PolynomLineIntegration2(Node[CellIndex+1][JNu-1],Node[CellIndex+1][JNu  ],
						0.0, 0.0, 1, 0) )*0.5,                // cell East side & Division by (OrderX+1)
		       (ExtendWest_BndNorthSplineInfo[CellIndex].YCentroidContribution() +              // cell North side
			PolynomLineIntegration2(Node[CellIndex  ][JNu  ],Node[CellIndex  ][JNu-1],
						0.0, 0.0, 0, 1) +                                        // cell West side
			PolynomLineIntegration2(Node[CellIndex  ][JNu-1],Node[CellIndex+1][JNu-1],
						0.0, 0.0, 0, 1) +                                        // cell South side
			PolynomLineIntegration2(Node[CellIndex+1][JNu-1],Node[CellIndex+1][JNu  ],
						0.0, 0.0, 0, 1) )) /Cell[CellIndex][JCu].A; // cell East side & Division by A  

    case EXTEND_E_NORTH_RIGHT_SPLINE: // Use only the extension of North spline to East (right side)
      return Vector2D(( PolynomLineIntegration2(Node[CellIndex+1][JNu+1],Node[CellIndex  ][JNu+1],
						0.0, 0.0, 1, 0) +                           // cell North side
			PolynomLineIntegration2(Node[CellIndex  ][JNu+1],Node[CellIndex  ][JNu  ],
						0.0, 0.0, 1, 0) -                           // cell West side
			ExtendEast_BndNorthSplineInfo[CellIndex-(ICu+1)].XCentroidContribution() + // cell South side
			PolynomLineIntegration2(Node[CellIndex+1][JNu  ],Node[CellIndex+1][JNu+1],
						0.0, 0.0, 1, 0) )*0.5,         // cell East side & Division by (OrderX+1)
		      (PolynomLineIntegration2(Node[CellIndex+1][JNu+1],Node[CellIndex  ][JNu+1],
					       0.0, 0.0, 0, 1) +                             // cell North side
		       PolynomLineIntegration2(Node[CellIndex  ][JNu+1],Node[CellIndex  ][JNu  ],
					       0.0, 0.0, 0, 1) -                             // cell West side
		       ExtendEast_BndNorthSplineInfo[CellIndex-(ICu+1)].YCentroidContribution() + // cell South side
		       PolynomLineIntegration2(Node[CellIndex+1][JNu  ],Node[CellIndex+1][JNu+1],
					       0.0, 0.0, 0, 1)) ) /Cell[CellIndex][JCu+1].A;    // cell East side & Division by A
      
    case EXTEND_W_NORTH_LEFT_SPLINE: // Use only the extension of North spline to West (left side)
      return Vector2D(( PolynomLineIntegration2(Node[CellIndex+1][JNu+1],Node[CellIndex  ][JNu+1],
						0.0, 0.0, 1, 0) +                           // cell North side
			PolynomLineIntegration2(Node[CellIndex  ][JNu+1],Node[CellIndex  ][JNu  ],
						0.0, 0.0, 1, 0) -                           // cell West side
			ExtendWest_BndNorthSplineInfo[CellIndex].XCentroidContribution() + // cell South side
			PolynomLineIntegration2(Node[CellIndex+1][JNu  ],Node[CellIndex+1][JNu+1],
						0.0, 0.0, 1, 0) )*0.5,         // cell East side & Division by (OrderX+1)
		      (PolynomLineIntegration2(Node[CellIndex+1][JNu+1],Node[CellIndex  ][JNu+1],
					       0.0, 0.0, 0, 1) +                             // cell North side
		       PolynomLineIntegration2(Node[CellIndex  ][JNu+1],Node[CellIndex  ][JNu  ],
					       0.0, 0.0, 0, 1) -                             // cell West side
		       ExtendWest_BndNorthSplineInfo[CellIndex].YCentroidContribution() + // cell South side
		       PolynomLineIntegration2(Node[CellIndex+1][JNu  ],Node[CellIndex+1][JNu+1],
					       0.0, 0.0, 0, 1)) ) /Cell[CellIndex][JCu+1].A;    // cell East side & Division by A

    case EXTEND_E_NORTH_LEFT_SPLINE: // Use only the extension of North spline to East (left side)
      return Vector2D( (ExtendEast_BndNorthSplineInfo[CellIndex-(ICu+1)].XCentroidContribution() +              // cell North side
			PolynomLineIntegration2(Node[CellIndex  ][JNu  ],Node[CellIndex  ][JNu-1],
						0.0, 0.0, 1, 0) +                                        // cell West side
			PolynomLineIntegration2(Node[CellIndex  ][JNu-1],Node[CellIndex+1][JNu-1],
						0.0, 0.0, 1, 0) +                                        // cell South side
			PolynomLineIntegration2(Node[CellIndex+1][JNu-1],Node[CellIndex+1][JNu  ],
						0.0, 0.0, 1, 0) )*0.5,                // cell East side & Division by (OrderX+1)
		       (ExtendEast_BndNorthSplineInfo[CellIndex-(ICu+1)].YCentroidContribution() +              // cell North side
			PolynomLineIntegration2(Node[CellIndex  ][JNu  ],Node[CellIndex  ][JNu-1],
						0.0, 0.0, 0, 1) +                                        // cell West side
			PolynomLineIntegration2(Node[CellIndex  ][JNu-1],Node[CellIndex+1][JNu-1],
						0.0, 0.0, 0, 1) +                                        // cell South side
			PolynomLineIntegration2(Node[CellIndex+1][JNu-1],Node[CellIndex+1][JNu  ],
						0.0, 0.0, 0, 1) )) /Cell[CellIndex][JCu].A; // cell East side & Division by A  

    case EXTEND_W_SOUTH_RIGHT_SPLINE: // Use only the extension of South spline to West (right side)
      return Vector2D( (-ExtendWest_BndSouthSplineInfo[CellIndex].XCentroidContribution() +              // cell North side
			PolynomLineIntegration2(Node[CellIndex  ][JNl  ],Node[CellIndex  ][JNl-1],
						0.0, 0.0, 1, 0) +                                        // cell West side
			PolynomLineIntegration2(Node[CellIndex  ][JNl-1],Node[CellIndex+1][JNl-1],
						0.0, 0.0, 1, 0) +                                        // cell South side
			PolynomLineIntegration2(Node[CellIndex+1][JNl-1],Node[CellIndex+1][JNl  ],
						0.0, 0.0, 1, 0) )*0.5,                // cell East side & Division by (OrderX+1)
		       (-ExtendWest_BndSouthSplineInfo[CellIndex].YCentroidContribution() +              // cell North side
			PolynomLineIntegration2(Node[CellIndex  ][JNl  ],Node[CellIndex  ][JNl-1],
						0.0, 0.0, 0, 1) +                                        // cell West side
			PolynomLineIntegration2(Node[CellIndex  ][JNl-1],Node[CellIndex+1][JNl-1],
						0.0, 0.0, 0, 1) +                                        // cell South side
			PolynomLineIntegration2(Node[CellIndex+1][JNl-1],Node[CellIndex+1][JNl  ],
						0.0, 0.0, 0, 1) )) /Cell[CellIndex][JCl-1].A; // cell East side & Division by A  

    case EXTEND_E_SOUTH_RIGHT_SPLINE: // Use only the extension of South spline to East (right side)
      return Vector2D(( PolynomLineIntegration2(Node[CellIndex+1][JNl+1],Node[CellIndex  ][JNl+1],
						0.0, 0.0, 1, 0) +                           // cell North side
			PolynomLineIntegration2(Node[CellIndex  ][JNl+1],Node[CellIndex  ][JNl  ],
						0.0, 0.0, 1, 0) +                           // cell West side
			ExtendEast_BndSouthSplineInfo[CellIndex-(ICu+1)].XCentroidContribution() + // cell South side
			PolynomLineIntegration2(Node[CellIndex+1][JNl  ],Node[CellIndex+1][JNl+1],
						0.0, 0.0, 1, 0) )*0.5,         // cell East side & Division by (OrderX+1)
		      (PolynomLineIntegration2(Node[CellIndex+1][JNl+1],Node[CellIndex  ][JNl+1],
					       0.0, 0.0, 0, 1) +                             // cell North side
		       PolynomLineIntegration2(Node[CellIndex  ][JNl+1],Node[CellIndex  ][JNl  ],
					       0.0, 0.0, 0, 1) +                             // cell West side
		       ExtendEast_BndSouthSplineInfo[CellIndex-(ICu+1)].YCentroidContribution() + // cell South side
		       PolynomLineIntegration2(Node[CellIndex+1][JNl  ],Node[CellIndex+1][JNl+1],
					       0.0, 0.0, 0, 1)) ) /Cell[CellIndex][JCl].A;    // cell East side & Division by A
      
    case EXTEND_W_SOUTH_LEFT_SPLINE: // Use only the extension of South spline to West (left side)
      return Vector2D(( PolynomLineIntegration2(Node[CellIndex+1][JNl+1],Node[CellIndex  ][JNl+1],
						0.0, 0.0, 1, 0) +                           // cell North side
			PolynomLineIntegration2(Node[CellIndex  ][JNl+1],Node[CellIndex  ][JNl  ],
						0.0, 0.0, 1, 0) +                           // cell West side
			ExtendWest_BndSouthSplineInfo[CellIndex].XCentroidContribution() + // cell South side
			PolynomLineIntegration2(Node[CellIndex+1][JNl  ],Node[CellIndex+1][JNl+1],
						0.0, 0.0, 1, 0) )*0.5,         // cell East side & Division by (OrderX+1)
		      (PolynomLineIntegration2(Node[CellIndex+1][JNl+1],Node[CellIndex  ][JNl+1],
					       0.0, 0.0, 0, 1) +                             // cell North side
		       PolynomLineIntegration2(Node[CellIndex  ][JNl+1],Node[CellIndex  ][JNl  ],
					       0.0, 0.0, 0, 1) +                             // cell West side
		       ExtendWest_BndSouthSplineInfo[CellIndex].YCentroidContribution() + // cell South side
		       PolynomLineIntegration2(Node[CellIndex+1][JNl  ],Node[CellIndex+1][JNl+1],
					       0.0, 0.0, 0, 1)) ) /Cell[CellIndex][JCl].A;    // cell East side & Division by A

    case EXTEND_E_SOUTH_LEFT_SPLINE: // Use only the extension of South spline to East (left side)
      return Vector2D( (-ExtendEast_BndSouthSplineInfo[CellIndex-(ICu+1)].XCentroidContribution() +              // cell North side
			PolynomLineIntegration2(Node[CellIndex  ][JNl  ],Node[CellIndex  ][JNl-1],
						0.0, 0.0, 1, 0) +                                        // cell West side
			PolynomLineIntegration2(Node[CellIndex  ][JNl-1],Node[CellIndex+1][JNl-1],
						0.0, 0.0, 1, 0) +                                        // cell South side
			PolynomLineIntegration2(Node[CellIndex+1][JNl-1],Node[CellIndex+1][JNl  ],
						0.0, 0.0, 1, 0) )*0.5,                // cell East side & Division by (OrderX+1)
		       (-ExtendEast_BndSouthSplineInfo[CellIndex-(ICu+1)].YCentroidContribution() +              // cell North side
			PolynomLineIntegration2(Node[CellIndex  ][JNl  ],Node[CellIndex  ][JNl-1],
						0.0, 0.0, 0, 1) +                                        // cell West side
			PolynomLineIntegration2(Node[CellIndex  ][JNl-1],Node[CellIndex+1][JNl-1],
						0.0, 0.0, 0, 1) +                                        // cell South side
			PolynomLineIntegration2(Node[CellIndex+1][JNl-1],Node[CellIndex+1][JNl  ],
						0.0, 0.0, 0, 1) )) /Cell[CellIndex][JCl-1].A; // cell East side & Division by A  

    case EXTEND_N_EAST_RIGHT_SPLINE: // Use only the extension of East spline to North (right side)
      return Vector2D( (PolynomLineIntegration2(Node[INu  ][CellIndex+1],Node[INu-1][CellIndex+1],
						0.0, 0.0, 1, 0) +                                        // cell North side
			PolynomLineIntegration2(Node[INu-1][CellIndex+1],Node[INu-1][CellIndex  ],
						0.0, 0.0, 1, 0) +                                        // cell West side
			PolynomLineIntegration2(Node[INu-1][CellIndex  ],Node[INu  ][CellIndex  ],
						0.0, 0.0, 1, 0) +                                        // cell South side
			ExtendNorth_BndEastSplineInfo[CellIndex-(JCu+1)].XCentroidContribution() ) *0.5,
		       // cell East side & Division by (OrderX+1)
		       (PolynomLineIntegration2(Node[INu  ][CellIndex+1],Node[INu-1][CellIndex+1],
						0.0, 0.0, 0, 1) +                                        // cell North side
			PolynomLineIntegration2(Node[INu-1][CellIndex+1],Node[INu-1][CellIndex  ],
						0.0, 0.0, 0, 1) +                                        // cell West side
			PolynomLineIntegration2(Node[INu-1][CellIndex  ],Node[INu  ][CellIndex  ],
						0.0, 0.0, 0, 1) +                                        // cell South side
			ExtendNorth_BndEastSplineInfo[CellIndex-(JCu+1)].YCentroidContribution() )) /Cell[ICu][CellIndex].A;
      // cell East side & Division by A

    case EXTEND_S_EAST_RIGHT_SPLINE: // Use only the extension of East spline to South (right side)
      return Vector2D( (PolynomLineIntegration2(Node[INu+1][CellIndex+1],Node[INu  ][CellIndex+1],
						0.0, 0.0, 1, 0) -                                         // cell North side
			ExtendSouth_BndEastSplineInfo[CellIndex].XCentroidContribution() +                // cell West side
			PolynomLineIntegration2(Node[INu  ][CellIndex  ],Node[INu+1][CellIndex  ],
						0.0, 0.0, 1, 0) +                                         // cell South side
			PolynomLineIntegration2(Node[INu+1][CellIndex  ],Node[INu+1][CellIndex+1],
						0.0, 0.0, 1, 0) ) *0.5, // cell East side & Division by (OrderX+1)
		       (PolynomLineIntegration2(Node[INu+1][CellIndex+1],Node[INu  ][CellIndex+1],
						0.0, 0.0, 0, 1) -                                         // cell North side
			ExtendSouth_BndEastSplineInfo[CellIndex].YCentroidContribution() +                // cell West side
			PolynomLineIntegration2(Node[INu  ][CellIndex  ],Node[INu+1][CellIndex  ],
						0.0, 0.0, 0, 1) +                                         // cell South side
			PolynomLineIntegration2(Node[INu+1][CellIndex  ],Node[INu+1][CellIndex+1],
						0.0, 0.0, 0, 1) )) /Cell[ICu+1][CellIndex].A; // cell East side & Division by A

    case EXTEND_N_EAST_LEFT_SPLINE: // Use only the extension of East spline to North (left side)
      return Vector2D( (PolynomLineIntegration2(Node[INu+1][CellIndex+1],Node[INu  ][CellIndex+1],
						0.0, 0.0, 1, 0) -                                         // cell North side
			ExtendNorth_BndEastSplineInfo[CellIndex-(JCu+1)].XCentroidContribution() +                // cell West side
			PolynomLineIntegration2(Node[INu  ][CellIndex  ],Node[INu+1][CellIndex  ],
						0.0, 0.0, 1, 0) +                                         // cell South side
			PolynomLineIntegration2(Node[INu+1][CellIndex  ],Node[INu+1][CellIndex+1],
						0.0, 0.0, 1, 0) ) *0.5, // cell East side & Division by (OrderX+1)
		       (PolynomLineIntegration2(Node[INu+1][CellIndex+1],Node[INu  ][CellIndex+1],
						0.0, 0.0, 0, 1) -                                         // cell North side
			ExtendNorth_BndEastSplineInfo[CellIndex-(JCu+1)].YCentroidContribution() +                // cell West side
			PolynomLineIntegration2(Node[INu  ][CellIndex  ],Node[INu+1][CellIndex  ],
						0.0, 0.0, 0, 1) +                                         // cell South side
			PolynomLineIntegration2(Node[INu+1][CellIndex  ],Node[INu+1][CellIndex+1],
						0.0, 0.0, 0, 1) )) /Cell[ICu+1][CellIndex].A; // cell East side & Division by A

    case EXTEND_S_EAST_LEFT_SPLINE: // Use only the extension of East spline to South (left side)
      return Vector2D( (PolynomLineIntegration2(Node[INu  ][CellIndex+1],Node[INu-1][CellIndex+1],
						0.0, 0.0, 1, 0) +                                        // cell North side
			PolynomLineIntegration2(Node[INu-1][CellIndex+1],Node[INu-1][CellIndex  ],
						0.0, 0.0, 1, 0) +                                        // cell West side
			PolynomLineIntegration2(Node[INu-1][CellIndex  ],Node[INu  ][CellIndex  ],
						0.0, 0.0, 1, 0) +                                        // cell South side
			ExtendSouth_BndEastSplineInfo[CellIndex].XCentroidContribution() ) *0.5,
		       // cell East side & Division by (OrderX+1)
		       (PolynomLineIntegration2(Node[INu  ][CellIndex+1],Node[INu-1][CellIndex+1],
						0.0, 0.0, 0, 1) +                                        // cell North side
			PolynomLineIntegration2(Node[INu-1][CellIndex+1],Node[INu-1][CellIndex  ],
						0.0, 0.0, 0, 1) +                                        // cell West side
			PolynomLineIntegration2(Node[INu-1][CellIndex  ],Node[INu  ][CellIndex  ],
						0.0, 0.0, 0, 1) +                                        // cell South side
			ExtendSouth_BndEastSplineInfo[CellIndex].YCentroidContribution() )) /Cell[ICu][CellIndex].A;
      // cell East side & Division by A

    case EXTEND_N_WEST_RIGHT_SPLINE: // Use only the extension of West spline to North (right side)
      return Vector2D( (PolynomLineIntegration2(Node[INl  ][CellIndex+1],Node[INl-1][CellIndex+1],
						0.0, 0.0, 1, 0) +                                        // cell North side
			PolynomLineIntegration2(Node[INl-1][CellIndex+1],Node[INl-1][CellIndex  ],
						0.0, 0.0, 1, 0) +                                        // cell West side
			PolynomLineIntegration2(Node[INl-1][CellIndex  ],Node[INl  ][CellIndex  ],
						0.0, 0.0, 1, 0) -                                        // cell South side
			ExtendNorth_BndWestSplineInfo[CellIndex-(JCu+1)].XCentroidContribution() ) *0.5,
		       // cell East side & Division by (OrderX+1)
		       (PolynomLineIntegration2(Node[INl  ][CellIndex+1],Node[INl-1][CellIndex+1],
						0.0, 0.0, 0, 1) +                                        // cell North side
			PolynomLineIntegration2(Node[INl-1][CellIndex+1],Node[INl-1][CellIndex  ],
						0.0, 0.0, 0, 1) +                                        // cell West side
			PolynomLineIntegration2(Node[INl-1][CellIndex  ],Node[INl  ][CellIndex  ],
						0.0, 0.0, 0, 1) -                                        // cell South side
			ExtendNorth_BndWestSplineInfo[CellIndex-(JCu+1)].YCentroidContribution() )) /Cell[ICl-1][CellIndex].A;
      // cell East side & Division by A

    case EXTEND_S_WEST_RIGHT_SPLINE: // Use only the extension of West spline to South (right side)
      return Vector2D( (PolynomLineIntegration2(Node[INl+1][CellIndex+1],Node[INl  ][CellIndex+1],
						0.0, 0.0, 1, 0) +                                         // cell North side
			ExtendSouth_BndWestSplineInfo[CellIndex].XCentroidContribution() +                // cell West side
			PolynomLineIntegration2(Node[INl  ][CellIndex  ],Node[INl+1][CellIndex  ],
						0.0, 0.0, 1, 0) +                                         // cell South side
			PolynomLineIntegration2(Node[INl+1][CellIndex  ],Node[INl+1][CellIndex+1],
						0.0, 0.0, 1, 0) ) *0.5, // cell East side & Division by (OrderX+1)
		       (PolynomLineIntegration2(Node[INl+1][CellIndex+1],Node[INl  ][CellIndex+1],
						0.0, 0.0, 0, 1) +                                         // cell North side
			ExtendSouth_BndWestSplineInfo[CellIndex].YCentroidContribution() +                // cell West side
			PolynomLineIntegration2(Node[INl  ][CellIndex  ],Node[INl+1][CellIndex  ],
						0.0, 0.0, 0, 1) +                                         // cell South side
			PolynomLineIntegration2(Node[INl+1][CellIndex  ],Node[INl+1][CellIndex+1],
						0.0, 0.0, 0, 1) )) /Cell[ICl][CellIndex].A; // cell East side & Division by A

    case EXTEND_N_WEST_LEFT_SPLINE: // Use only the extension of West spline to North (left side)
      return Vector2D( (PolynomLineIntegration2(Node[INl+1][CellIndex+1],Node[INl  ][CellIndex+1],
						0.0, 0.0, 1, 0) +                                         // cell North side
			ExtendNorth_BndWestSplineInfo[CellIndex-(JCu+1)].XCentroidContribution() +                // cell West side
			PolynomLineIntegration2(Node[INl  ][CellIndex  ],Node[INl+1][CellIndex  ],
						0.0, 0.0, 1, 0) +                                         // cell South side
			PolynomLineIntegration2(Node[INl+1][CellIndex  ],Node[INl+1][CellIndex+1],
						0.0, 0.0, 1, 0) ) *0.5, // cell East side & Division by (OrderX+1)
		       (PolynomLineIntegration2(Node[INl+1][CellIndex+1],Node[INl  ][CellIndex+1],
						0.0, 0.0, 0, 1) +                                         // cell North side
			ExtendNorth_BndWestSplineInfo[CellIndex-(JCu+1)].YCentroidContribution() +                // cell West side
			PolynomLineIntegration2(Node[INl  ][CellIndex  ],Node[INl+1][CellIndex  ],
						0.0, 0.0, 0, 1) +                                         // cell South side
			PolynomLineIntegration2(Node[INl+1][CellIndex  ],Node[INl+1][CellIndex+1],
						0.0, 0.0, 0, 1) )) /Cell[ICl][CellIndex].A; // cell East side & Division by A

    case EXTEND_S_WEST_LEFT_SPLINE: // Use only the extension of West spline to South (left side)
      return Vector2D( (PolynomLineIntegration2(Node[INl  ][CellIndex+1],Node[INl-1][CellIndex+1],
						0.0, 0.0, 1, 0) +                                        // cell North side
			PolynomLineIntegration2(Node[INl-1][CellIndex+1],Node[INl-1][CellIndex  ],
						0.0, 0.0, 1, 0) +                                        // cell West side
			PolynomLineIntegration2(Node[INl-1][CellIndex  ],Node[INl  ][CellIndex  ],
						0.0, 0.0, 1, 0) -                                        // cell South side
			ExtendSouth_BndWestSplineInfo[CellIndex].XCentroidContribution() ) *0.5,
		       // cell East side & Division by (OrderX+1)
		       (PolynomLineIntegration2(Node[INl  ][CellIndex+1],Node[INl-1][CellIndex+1],
						0.0, 0.0, 0, 1) +                                        // cell North side
			PolynomLineIntegration2(Node[INl-1][CellIndex+1],Node[INl-1][CellIndex  ],
						0.0, 0.0, 0, 1) +                                        // cell West side
			PolynomLineIntegration2(Node[INl-1][CellIndex  ],Node[INl  ][CellIndex  ],
						0.0, 0.0, 0, 1) -                                        // cell South side
			ExtendSouth_BndWestSplineInfo[CellIndex].YCentroidContribution() )) /Cell[ICl-1][CellIndex].A;
      // cell East side & Division by A

    case CORNER_NORTH_EXTEND_N_WEST_SPLINES: // Use the North spline and the extension of West spline to North
      return Vector2D( (PolynomLineIntegration2(Node[INl+1][JNu+1],Node[INl  ][JNu+1],
						0.0, 0.0, 1, 0) +                                         // cell North side
			ExtendNorth_BndWestSplineInfo[0].XCentroidContribution() -                // cell West side
			BndNorthSplineInfo[ICl].XCentroidContribution() +                         // cell South side
			PolynomLineIntegration2(Node[INl+1][JNu  ],Node[INl+1][JNu+1],
						0.0, 0.0, 1, 0) ) *0.5, // cell East side & Division by (OrderX+1)
		       (PolynomLineIntegration2(Node[INl+1][JNu+1],Node[INl  ][JNu+1],
						0.0, 0.0, 0, 1) +                                         // cell North side
			ExtendNorth_BndWestSplineInfo[0].YCentroidContribution() -                // cell West side
			BndNorthSplineInfo[ICl].YCentroidContribution() +                         // cell South side
			PolynomLineIntegration2(Node[INl+1][JNu  ],Node[INl+1][JNu+1],
						0.0, 0.0, 0, 1) )) /Cell[ICl][JCu+1].A; // cell East side & Division by A

    case CORNER_EXTEND_N_WEST_EXTEND_W_NORTH_SPLINES: // Use the North extension of West and the West extension of North
      return Vector2D(( PolynomLineIntegration2(Node[INl  ][JNu+1],Node[INl-1][JNu+1],
						0.0, 0.0, 1, 0) +                           // cell North side
			PolynomLineIntegration2(Node[INl-1][JNu+1],Node[INl-1][JNu  ],
						0.0, 0.0, 1, 0) -                           // cell West side
			ExtendWest_BndNorthSplineInfo[ICl-1].XCentroidContribution() - // cell South side
			ExtendNorth_BndWestSplineInfo[0].XCentroidContribution() )*0.5,
		      // cell East side & Division by (OrderX+1)
		      (PolynomLineIntegration2(Node[INl  ][JNu+1],Node[INl-1][JNu+1],
					       0.0, 0.0, 0, 1) +                             // cell North side
		       PolynomLineIntegration2(Node[INl-1][JNu+1],Node[INl-1][JNu  ],
					       0.0, 0.0, 0, 1) -                             // cell West side
		       ExtendWest_BndNorthSplineInfo[ICl-1].YCentroidContribution() - // cell South side
		       ExtendNorth_BndWestSplineInfo[0].YCentroidContribution() )) /Cell[ICl-1][JCu+1].A;
      // cell East side & Division by A

    case CORNER_EXTEND_W_NORTH_WEST_SPLINES: // Use the West extension of North and the West spline
      return Vector2D( (ExtendWest_BndNorthSplineInfo[ICl-1].XCentroidContribution() +              // cell North side
			PolynomLineIntegration2(Node[INl-1][JNu  ],Node[INl-1][JNu-1],
						0.0, 0.0, 1, 0) +                                        // cell West side
			PolynomLineIntegration2(Node[INl-1][JNu-1],Node[INl  ][JNu-1],
						0.0, 0.0, 1, 0) -                                        // cell South side
			BndWestSplineInfo[JCu].XCentroidContribution() )*0.5,     // cell East side & Division by (OrderX+1)
		       (ExtendWest_BndNorthSplineInfo[ICl-1].YCentroidContribution() +              // cell North side
			PolynomLineIntegration2(Node[INl-1][JNu  ],Node[INl-1][JNu-1],
						0.0, 0.0, 0, 1) +                                        // cell West side
			PolynomLineIntegration2(Node[INl-1][JNu-1],Node[INl  ][JNu-1],
						0.0, 0.0, 0, 1) -                                        // cell South side
			BndWestSplineInfo[JCu].YCentroidContribution() )) /Cell[ICl-1][JCu].A;
      // cell East side & Division by A  

    case CORNER_WEST_EXTEND_W_SOUTH_SPLINES: // Use the West spline and the West extension of South
      return Vector2D(( PolynomLineIntegration2(Node[INl  ][JNl+1],Node[INl-1][JNl+1],
						0.0, 0.0, 1, 0) +                           // cell North side
			PolynomLineIntegration2(Node[INl-1][JNl+1],Node[INl-1][JNl  ],
						0.0, 0.0, 1, 0) +                           // cell West side
			ExtendWest_BndSouthSplineInfo[ICl-1].XCentroidContribution() - // cell South side
			BndWestSplineInfo[JCl].XCentroidContribution() )*0.5,         // cell East side & Division by (OrderX+1)
		      (PolynomLineIntegration2(Node[INl  ][JNl+1],Node[INl-1][JNl+1],
					       0.0, 0.0, 0, 1) +                             // cell North side
		       PolynomLineIntegration2(Node[INl-1][JNl+1],Node[INl-1][JNl  ],
					       0.0, 0.0, 0, 1) +                             // cell West side
		       ExtendWest_BndSouthSplineInfo[ICl-1].YCentroidContribution() - // cell South side
		       BndWestSplineInfo[JCl].YCentroidContribution() )) /Cell[ICl-1][JCl].A;
      // cell East side & Division by A

    case CORNER_EXTEND_W_SOUTH_EXTEND_S_WEST_SPLINES: // Use the West extension of South and the South extension of West
      return Vector2D( (-ExtendWest_BndSouthSplineInfo[ICl-1].XCentroidContribution() +              // cell North side
			PolynomLineIntegration2(Node[INl-1][JNl  ],Node[INl-1][JNl-1],
						0.0, 0.0, 1, 0) +                                        // cell West side
			PolynomLineIntegration2(Node[INl-1][JNl-1],Node[INl  ][JNl-1],
						0.0, 0.0, 1, 0) -                                        // cell South side
			ExtendSouth_BndWestSplineInfo[JCl-1].XCentroidContribution() )*0.5, 
		       // cell East side & Division by (OrderX+1)
		       (-ExtendWest_BndSouthSplineInfo[ICl-1].YCentroidContribution() +              // cell North side
			PolynomLineIntegration2(Node[INl-1][JNl  ],Node[INl-1][JNl-1],
						0.0, 0.0, 0, 1) +                                        // cell West side
			PolynomLineIntegration2(Node[INl-1][JNl-1],Node[INl  ][JNl-1],
						0.0, 0.0, 0, 1) -                                        // cell South side
			ExtendSouth_BndWestSplineInfo[JCl-1].YCentroidContribution() )) /Cell[ICl-1][JCl-1].A;
      // cell East side & Division by A  

    case CORNER_EXTEND_S_WEST_SOUTH_SPLINES: // Use the South extension of West and the South spline
      return Vector2D( (-BndSouthSplineInfo[ICl].XCentroidContribution() +                                // cell North side
			ExtendSouth_BndWestSplineInfo[JCl-1].XCentroidContribution() +                // cell West side
			PolynomLineIntegration2(Node[INl  ][JNl-1],Node[INl+1][JNl-1],
						0.0, 0.0, 1, 0) +                                         // cell South side
			PolynomLineIntegration2(Node[INl+1][JNl-1],Node[INl+1][JNl],
						0.0, 0.0, 1, 0) ) *0.5, // cell East side & Division by (OrderX+1)
		       (-BndSouthSplineInfo[ICl].YCentroidContribution() +                                // cell North side
			ExtendSouth_BndWestSplineInfo[JCl-1].YCentroidContribution() +                // cell West side
			PolynomLineIntegration2(Node[INl  ][JNl-1],Node[INl+1][JNl-1],
						0.0, 0.0, 0, 1) +                                         // cell South side
			PolynomLineIntegration2(Node[INl+1][JNl-1],Node[INl+1][JNl],
						0.0, 0.0, 0, 1) )) /Cell[ICl][JCl-1].A; // cell East side & Division by A

    case CORNER_SOUTH_EXTEND_S_EAST_SPLINES: // Use the South spline and the South extension of East
      return Vector2D( (-BndSouthSplineInfo[ICu].XCentroidContribution() +                               // cell North side
			PolynomLineIntegration2(Node[INu-1][JNl  ],Node[INu-1][JNl-1],
						0.0, 0.0, 1, 0) +                                        // cell West side
			PolynomLineIntegration2(Node[INu-1][JNl-1],Node[INu  ][JNl-1],
						0.0, 0.0, 1, 0) +                                        // cell South side
			ExtendSouth_BndEastSplineInfo[JCl-1].XCentroidContribution() ) *0.5,
		       // cell East side & Division by (OrderX+1)
		       (-BndSouthSplineInfo[ICu].YCentroidContribution() +                                // cell North side
			PolynomLineIntegration2(Node[INu-1][JNl  ],Node[INu-1][JNl-1],
						0.0, 0.0, 0, 1) +                                        // cell West side
			PolynomLineIntegration2(Node[INu-1][JNl-1],Node[INu  ][JNl-1],
						0.0, 0.0, 0, 1) +                                        // cell South side
			ExtendSouth_BndEastSplineInfo[JCl-1].YCentroidContribution() )) /Cell[ICu][JCl-1].A;

    case CORNER_EXTEND_S_EAST_EXTEND_E_SOUTH_SPLINES: // Use the South extension of East and the East extension of South
      return Vector2D( (-ExtendEast_BndSouthSplineInfo[0].XCentroidContribution() -              // cell North side
			ExtendSouth_BndEastSplineInfo[JCl-1].XCentroidContribution() +           // cell West side
			PolynomLineIntegration2(Node[INu  ][JNl-1],Node[INu+1][JNl-1],
						0.0, 0.0, 1, 0) +                                        // cell South side
			PolynomLineIntegration2(Node[INu+1][JNl-1],Node[INu+1][JNl  ],
						0.0, 0.0, 1, 0) )*0.5,                // cell East side & Division by (OrderX+1)
		       (-ExtendEast_BndSouthSplineInfo[0].YCentroidContribution() -              // cell North side
			ExtendSouth_BndEastSplineInfo[JCl-1].YCentroidContribution() +           // cell West side
			PolynomLineIntegration2(Node[INu  ][JNl-1],Node[INu+1][JNl-1],
						0.0, 0.0, 0, 1) +                                        // cell South side
			PolynomLineIntegration2(Node[INu+1][JNl-1],Node[INu+1][JNl  ],
						0.0, 0.0, 0, 1) )) /Cell[ICu+1][JCl-1].A; // cell East side & Division by A  

    case CORNER_EXTEND_E_SOUTH_EAST_SPLINES: // Use the East extension of South and the East spline
      return Vector2D(( PolynomLineIntegration2(Node[INu+1][JNl+1],Node[INu  ][JNl+1],
						0.0, 0.0, 1, 0) -                           // cell North side
			BndEastSplineInfo[JCl].XCentroidContribution() +                    // cell West side
			ExtendEast_BndSouthSplineInfo[0].XCentroidContribution() + // cell South side
			PolynomLineIntegration2(Node[INu+1][JNl  ],Node[INu+1][JNl+1],
						0.0, 0.0, 1, 0) )*0.5,         // cell East side & Division by (OrderX+1)
		      (PolynomLineIntegration2(Node[INu+1][JNl+1],Node[INu  ][JNl+1],
					       0.0, 0.0, 0, 1) -                             // cell North side
		       BndEastSplineInfo[JCl].YCentroidContribution() +                    // cell West side
		       ExtendEast_BndSouthSplineInfo[0].YCentroidContribution() + // cell South side
		       PolynomLineIntegration2(Node[INu+1][JNl  ],Node[INu+1][JNl+1],
					       0.0, 0.0, 0, 1)) ) /Cell[ICu+1][JCl].A;    // cell East side & Division by A 

    case CORNER_EAST_EXTEND_E_NORTH_SPLINES: // Use the East spline and the East extension of North
      return Vector2D( (ExtendEast_BndNorthSplineInfo[0].XCentroidContribution() -              // cell North side
			BndEastSplineInfo[JCu].XCentroidContribution() +                       // cell West side
			PolynomLineIntegration2(Node[INu  ][JNu-1],Node[INu+1][JNu-1],
						0.0, 0.0, 1, 0) +                                        // cell South side
			PolynomLineIntegration2(Node[INu+1][JNu-1],Node[INu+1][JNu  ],
						0.0, 0.0, 1, 0) )*0.5,                // cell East side & Division by (OrderX+1)
		       (ExtendEast_BndNorthSplineInfo[0].YCentroidContribution() -              // cell North side
			BndEastSplineInfo[JCu].YCentroidContribution() +                       // cell West side
			PolynomLineIntegration2(Node[INu  ][JNu-1],Node[INu+1][JNu-1],
						0.0, 0.0, 0, 1) +                                        // cell South side
			PolynomLineIntegration2(Node[INu+1][JNu-1],Node[INu+1][JNu  ],
						0.0, 0.0, 0, 1) )) /Cell[ICu+1][JCu].A; // cell East side & Division by A  

    case CORNER_EXTEND_E_NORTH_EXTEND_N_EAST_SPLINES: // Use the East extention of North and North extension of East
      return Vector2D( (PolynomLineIntegration2(Node[INu+1][JNu+1],Node[INu  ][JNu+1],
						0.0, 0.0, 1, 0) -                                         // cell North side
			ExtendNorth_BndEastSplineInfo[0].XCentroidContribution() -                // cell West side
			ExtendEast_BndNorthSplineInfo[0].XCentroidContribution() +                // cell South side
			PolynomLineIntegration2(Node[INu+1][JNu  ],Node[INu+1][JNu+1],
						0.0, 0.0, 1, 0) ) *0.5, // cell East side & Division by (OrderX+1)
		       (PolynomLineIntegration2(Node[INu+1][JNu+1],Node[INu  ][JNu+1],
						0.0, 0.0, 0, 1) -                                         // cell North side
			ExtendNorth_BndEastSplineInfo[0].YCentroidContribution() -                // cell West side
			ExtendEast_BndNorthSplineInfo[0].YCentroidContribution() +                // cell South side
			PolynomLineIntegration2(Node[INu+1][JNu  ],Node[INu+1][JNu+1],
						0.0, 0.0, 0, 1) )) /Cell[ICu+1][JCu+1].A; // cell East side & Division by A

    case CORNER_EXTEND_N_EAST_NORTH_SPLINES: // Use the North extension of East and the North spline
      return Vector2D( (PolynomLineIntegration2(Node[INu  ][JNu+1],Node[INu-1][JNu+1],
						0.0, 0.0, 1, 0) +                                        // cell North side
			PolynomLineIntegration2(Node[INu-1][JNu+1],Node[INu-1][JNu  ],
						0.0, 0.0, 1, 0) -                                        // cell West side
			BndNorthSplineInfo[ICu].XCentroidContribution() +                                // cell South side
			ExtendNorth_BndEastSplineInfo[0].XCentroidContribution() ) *0.5,
		       // cell East side & Division by (OrderX+1)
		       (PolynomLineIntegration2(Node[INu  ][JNu+1],Node[INu-1][JNu+1],
						0.0, 0.0, 0, 1) +                                        // cell North side
			PolynomLineIntegration2(Node[INu-1][JNu+1],Node[INu-1][JNu  ],
						0.0, 0.0, 0, 1) -                                        // cell West side
			BndNorthSplineInfo[ICu].YCentroidContribution() +                                // cell South side
			ExtendNorth_BndEastSplineInfo[0].YCentroidContribution() )) /Cell[ICu][JCu+1].A;

    default:
      return Vector2D(0.0);
    } // endswitch

  } else {

    // Use the boundary spline functions to integrate along curved edges.

    switch(Boundary){

    case NORTH_SPLINE:          // Use only the North Spline -> (iCell,jCell)=(CellIndex,JCu+1)
      return Vector2D(( PolynomLineIntegration2(Node[CellIndex+1][JNu+1],Node[CellIndex  ][JNu+1],
						0.0, 0.0, 1, 0) +                           // cell North side
			PolynomLineIntegration2(Node[CellIndex  ][JNu+1],Node[CellIndex  ][JNu  ],
						0.0, 0.0, 1, 0) +                           // cell West side
			BndNorthSpline.PolynomOrderIntegration(Node[CellIndex ][JNu], Node[CellIndex+1][JNu], 
							       Vector2D(0.0,0.0),15,1,0) + // cell South side
			PolynomLineIntegration2(Node[CellIndex+1][JNu  ],Node[CellIndex+1][JNu+1],
						0.0, 0.0, 1, 0) )*0.5,         // cell East side & Division by (OrderX+1)
		      (PolynomLineIntegration2(Node[CellIndex+1][JNu+1],Node[CellIndex  ][JNu+1],
					       0.0, 0.0, 0, 1) +                             // cell North side
		       PolynomLineIntegration2(Node[CellIndex  ][JNu+1],Node[CellIndex  ][JNu  ],
					       0.0, 0.0, 0, 1) +                             // cell West side
		       BndNorthSpline.PolynomOrderIntegration(Node[CellIndex ][JNu], Node[CellIndex+1][JNu],
							      Vector2D(0.0,0.0), 15, 0, 1)+ // cell South side
		       PolynomLineIntegration2(Node[CellIndex+1][JNu  ],Node[CellIndex+1][JNu+1],
					       0.0, 0.0, 0, 1)) ) /Cell[CellIndex][JCu+1].A;    // cell East side & Division by A
    
    case SOUTH_SPLINE:        // Use only the South Spline
      return Vector2D( (BndSouthSpline.PolynomOrderIntegration(Node[CellIndex+1][JNl], Node[CellIndex][JNl],
							       Vector2D(0.0,0.0),15,1,0) +              // cell North side
			PolynomLineIntegration2(Node[CellIndex  ][JNl  ],Node[CellIndex  ][JNl-1],
						0.0, 0.0, 1, 0) +                                        // cell West side
			PolynomLineIntegration2(Node[CellIndex  ][JNl-1],Node[CellIndex+1][JNl-1],
						0.0, 0.0, 1, 0) +                                        // cell South side
			PolynomLineIntegration2(Node[CellIndex+1][JNl-1],Node[CellIndex+1][JNl  ],
						0.0, 0.0, 1, 0) )*0.5,                // cell East side & Division by (OrderX+1)
		       (BndSouthSpline.PolynomOrderIntegration(Node[CellIndex+1][JNl], Node[CellIndex][JNl],
							       Vector2D(0.0,0.0),15,0,1) +              // cell North side
			PolynomLineIntegration2(Node[CellIndex  ][JNl  ],Node[CellIndex  ][JNl-1],
						0.0, 0.0, 0, 1) +                                        // cell West side
			PolynomLineIntegration2(Node[CellIndex  ][JNl-1],Node[CellIndex+1][JNl-1],
						0.0, 0.0, 0, 1) +                                        // cell South side
			PolynomLineIntegration2(Node[CellIndex+1][JNl-1],Node[CellIndex+1][JNl  ],
						0.0, 0.0, 0, 1) )) /Cell[CellIndex][JCl-1].A; // cell East side & Division by A  

    case WEST_SPLINE:       // Use only the West Spline
      return Vector2D( (PolynomLineIntegration2(Node[INl  ][CellIndex+1],Node[INl-1][CellIndex+1],
						0.0, 0.0, 1, 0) +                                        // cell North side
			PolynomLineIntegration2(Node[INl-1][CellIndex+1],Node[INl-1][CellIndex  ],
						0.0, 0.0, 1, 0) +                                        // cell West side
			PolynomLineIntegration2(Node[INl-1][CellIndex  ],Node[INl  ][CellIndex  ],
						0.0, 0.0, 1, 0) +                                        // cell South side
			BndWestSpline.PolynomOrderIntegration(Node[INl][CellIndex], Node[INl][CellIndex+1],
							      Vector2D(0.0,0.0),15,1,0) ) *0.5, 
		       // cell East side & Division by (OrderX+1)
		       (PolynomLineIntegration2(Node[INl  ][CellIndex+1],Node[INl-1][CellIndex+1],
						0.0, 0.0, 0, 1) +                                        // cell North side
			PolynomLineIntegration2(Node[INl-1][CellIndex+1],Node[INl-1][CellIndex  ],
						0.0, 0.0, 0, 1) +                                        // cell West side
			PolynomLineIntegration2(Node[INl-1][CellIndex  ],Node[INl  ][CellIndex  ],
						0.0, 0.0, 0, 1) +                                        // cell South side
			BndWestSpline.PolynomOrderIntegration(Node[INl][CellIndex], Node[INl][CellIndex+1],
							      Vector2D(0.0,0.0),15,0,1) )) /Cell[ICl-1][CellIndex].A;
      // cell East side & Division by A

    case EAST_SPLINE:      // Use only the East Spline
      return Vector2D( (PolynomLineIntegration2(Node[INu+1][CellIndex+1],Node[INu  ][CellIndex+1],
						0.0, 0.0, 1, 0) +                                         // cell North side
			BndEastSpline.PolynomOrderIntegration(Node[INu][CellIndex+1], Node[INu][CellIndex],
							      Vector2D(0.0,0.0),15,1,0) +                // cell West side
			PolynomLineIntegration2(Node[INu  ][CellIndex  ],Node[INu+1][CellIndex  ],
						0.0, 0.0, 1, 0) +                                         // cell South side
			PolynomLineIntegration2(Node[INu+1][CellIndex  ],Node[INu+1][CellIndex+1],
						0.0, 0.0, 1, 0) ) *0.5, // cell East side & Division by (OrderX+1)
		       (PolynomLineIntegration2(Node[INu+1][CellIndex+1],Node[INu  ][CellIndex+1],
						0.0, 0.0, 0, 1) +                                         // cell North side
			BndEastSpline.PolynomOrderIntegration(Node[INu][CellIndex+1], Node[INu][CellIndex],
							      Vector2D(0.0,0.0),15,0,1) +                // cell West side
			PolynomLineIntegration2(Node[INu  ][CellIndex  ],Node[INu+1][CellIndex  ],
						0.0, 0.0, 0, 1) +                                         // cell South side
			PolynomLineIntegration2(Node[INu+1][CellIndex  ],Node[INu+1][CellIndex+1],
						0.0, 0.0, 0, 1) )) /Cell[ICu+1][CellIndex].A; // cell East side & Division by A

    case EXTEND_W_NORTH_RIGHT_SPLINE: // Use only the extension of North spline to West (right side)
      return Vector2D( (ExtendWest_BndNorthSpline.PolynomOrderIntegration(Node[CellIndex+1][JNu  ],Node[CellIndex  ][JNu  ],
									  Vector2D(0.0,0.0),15,1,0) +     // cell North side
			PolynomLineIntegration2(Node[CellIndex  ][JNu  ],Node[CellIndex  ][JNu-1],
						0.0, 0.0, 1, 0) +                                        // cell West side
			PolynomLineIntegration2(Node[CellIndex  ][JNu-1],Node[CellIndex+1][JNu-1],
						0.0, 0.0, 1, 0) +                                        // cell South side
			PolynomLineIntegration2(Node[CellIndex+1][JNu-1],Node[CellIndex+1][JNu  ],
						0.0, 0.0, 1, 0) )*0.5,                // cell East side & Division by (OrderX+1)
		       (ExtendWest_BndNorthSpline.PolynomOrderIntegration(Node[CellIndex+1][JNu  ],Node[CellIndex  ][JNu  ],
									  Vector2D(0.0,0.0),15,0,1) +     // cell North side
			PolynomLineIntegration2(Node[CellIndex  ][JNu  ],Node[CellIndex  ][JNu-1],
						0.0, 0.0, 0, 1) +                                        // cell West side
			PolynomLineIntegration2(Node[CellIndex  ][JNu-1],Node[CellIndex+1][JNu-1],
						0.0, 0.0, 0, 1) +                                        // cell South side
			PolynomLineIntegration2(Node[CellIndex+1][JNu-1],Node[CellIndex+1][JNu  ],
						0.0, 0.0, 0, 1) )) /Cell[CellIndex][JCu].A; // cell East side & Division by A  

    case EXTEND_E_NORTH_RIGHT_SPLINE: // Use only the extension of North spline to East (right side)
      return Vector2D(( PolynomLineIntegration2(Node[CellIndex+1][JNu+1],Node[CellIndex  ][JNu+1],
						0.0, 0.0, 1, 0) +                           // cell North side
			PolynomLineIntegration2(Node[CellIndex  ][JNu+1],Node[CellIndex  ][JNu  ],
						0.0, 0.0, 1, 0) +                           // cell West side
			ExtendEast_BndNorthSpline.PolynomOrderIntegration(Node[CellIndex  ][JNu  ],Node[CellIndex+1][JNu  ],
									  Vector2D(0.0,0.0),15,1,0) +    // cell South side
			PolynomLineIntegration2(Node[CellIndex+1][JNu  ],Node[CellIndex+1][JNu+1],
						0.0, 0.0, 1, 0) )*0.5,         // cell East side & Division by (OrderX+1)
		      (PolynomLineIntegration2(Node[CellIndex+1][JNu+1],Node[CellIndex  ][JNu+1],
					       0.0, 0.0, 0, 1) +                             // cell North side
		       PolynomLineIntegration2(Node[CellIndex  ][JNu+1],Node[CellIndex  ][JNu  ],
					       0.0, 0.0, 0, 1) +                             // cell West side
		       ExtendEast_BndNorthSpline.PolynomOrderIntegration(Node[CellIndex  ][JNu  ],Node[CellIndex+1][JNu  ],
									 Vector2D(0.0,0.0),15,0,1) +    // cell South side
		       PolynomLineIntegration2(Node[CellIndex+1][JNu  ],Node[CellIndex+1][JNu+1],
					       0.0, 0.0, 0, 1)) ) /Cell[CellIndex][JCu+1].A;    // cell East side & Division by A
      
    case EXTEND_W_NORTH_LEFT_SPLINE: // Use only the extension of North spline to West (left side)
      return Vector2D(( PolynomLineIntegration2(Node[CellIndex+1][JNu+1],Node[CellIndex  ][JNu+1],
						0.0, 0.0, 1, 0) +                           // cell North side
			PolynomLineIntegration2(Node[CellIndex  ][JNu+1],Node[CellIndex  ][JNu  ],
						0.0, 0.0, 1, 0) +                           // cell West side
			ExtendWest_BndNorthSpline.PolynomOrderIntegration(Node[CellIndex  ][JNu  ],Node[CellIndex+1][JNu  ],
									  Vector2D(0.0,0.0),15,1,0) + // cell South side
			PolynomLineIntegration2(Node[CellIndex+1][JNu  ],Node[CellIndex+1][JNu+1],
						0.0, 0.0, 1, 0) )*0.5,         // cell East side & Division by (OrderX+1)
		      (PolynomLineIntegration2(Node[CellIndex+1][JNu+1],Node[CellIndex  ][JNu+1],
					       0.0, 0.0, 0, 1) +                             // cell North side
		       PolynomLineIntegration2(Node[CellIndex  ][JNu+1],Node[CellIndex  ][JNu  ],
					       0.0, 0.0, 0, 1) +                             // cell West side
		       ExtendWest_BndNorthSpline.PolynomOrderIntegration(Node[CellIndex  ][JNu  ],Node[CellIndex+1][JNu  ],
									 Vector2D(0.0,0.0),15,0,1) + // cell South side
		       PolynomLineIntegration2(Node[CellIndex+1][JNu  ],Node[CellIndex+1][JNu+1],
					       0.0, 0.0, 0, 1)) ) /Cell[CellIndex][JCu+1].A;    // cell East side & Division by A

    case EXTEND_E_NORTH_LEFT_SPLINE: // Use only the extension of North spline to East (left side)
      return Vector2D( (ExtendEast_BndNorthSpline.PolynomOrderIntegration(Node[CellIndex+1][JNu  ],Node[CellIndex  ][JNu  ],
									  Vector2D(0.0,0.0),15,1,0) +       // cell North side
			PolynomLineIntegration2(Node[CellIndex  ][JNu  ],Node[CellIndex  ][JNu-1],
						0.0, 0.0, 1, 0) +                                        // cell West side
			PolynomLineIntegration2(Node[CellIndex  ][JNu-1],Node[CellIndex+1][JNu-1],
						0.0, 0.0, 1, 0) +                                        // cell South side
			PolynomLineIntegration2(Node[CellIndex+1][JNu-1],Node[CellIndex+1][JNu  ],
						0.0, 0.0, 1, 0) )*0.5,                // cell East side & Division by (OrderX+1)
		       (ExtendEast_BndNorthSpline.PolynomOrderIntegration(Node[CellIndex+1][JNu  ],Node[CellIndex  ][JNu  ],
									  Vector2D(0.0,0.0),15,0,1) +       // cell North side
			PolynomLineIntegration2(Node[CellIndex  ][JNu  ],Node[CellIndex  ][JNu-1],
						0.0, 0.0, 0, 1) +                                        // cell West side
			PolynomLineIntegration2(Node[CellIndex  ][JNu-1],Node[CellIndex+1][JNu-1],
						0.0, 0.0, 0, 1) +                                        // cell South side
			PolynomLineIntegration2(Node[CellIndex+1][JNu-1],Node[CellIndex+1][JNu  ],
						0.0, 0.0, 0, 1) )) /Cell[CellIndex][JCu].A; // cell East side & Division by A  

    case EXTEND_W_SOUTH_RIGHT_SPLINE: // Use only the extension of South spline to West (right side)
      return Vector2D( (ExtendWest_BndSouthSpline.PolynomOrderIntegration(Node[CellIndex+1][JNl  ],Node[CellIndex  ][JNl  ],
									  Vector2D(0.0,0.0),15,1,0) +     // cell North side
			PolynomLineIntegration2(Node[CellIndex  ][JNl  ],Node[CellIndex  ][JNl-1],
						0.0, 0.0, 1, 0) +                                        // cell West side
			PolynomLineIntegration2(Node[CellIndex  ][JNl-1],Node[CellIndex+1][JNl-1],
						0.0, 0.0, 1, 0) +                                        // cell South side
			PolynomLineIntegration2(Node[CellIndex+1][JNl-1],Node[CellIndex+1][JNl  ],
						0.0, 0.0, 1, 0) )*0.5,                // cell East side & Division by (OrderX+1)
		       (ExtendWest_BndSouthSpline.PolynomOrderIntegration(Node[CellIndex+1][JNl  ],Node[CellIndex  ][JNl  ],
									  Vector2D(0.0,0.0),15,0,1) +     // cell North side
			PolynomLineIntegration2(Node[CellIndex  ][JNl  ],Node[CellIndex  ][JNl-1],
						0.0, 0.0, 0, 1) +                                        // cell West side
			PolynomLineIntegration2(Node[CellIndex  ][JNl-1],Node[CellIndex+1][JNl-1],
						0.0, 0.0, 0, 1) +                                        // cell South side
			PolynomLineIntegration2(Node[CellIndex+1][JNl-1],Node[CellIndex+1][JNl  ],
						0.0, 0.0, 0, 1) )) /Cell[CellIndex][JCl-1].A; // cell East side & Division by A  

    case EXTEND_E_SOUTH_RIGHT_SPLINE: // Use only the extension of South spline to East (right side)
      return Vector2D(( PolynomLineIntegration2(Node[CellIndex+1][JNl+1],Node[CellIndex  ][JNl+1],
						0.0, 0.0, 1, 0) +                           // cell North side
			PolynomLineIntegration2(Node[CellIndex  ][JNl+1],Node[CellIndex  ][JNl  ],
						0.0, 0.0, 1, 0) +                           // cell West side
			ExtendEast_BndSouthSpline.PolynomOrderIntegration(Node[CellIndex  ][JNl  ],Node[CellIndex+1][JNl  ],
									  Vector2D(0.0,0.0),15,1,0) + // cell South side
			PolynomLineIntegration2(Node[CellIndex+1][JNl  ],Node[CellIndex+1][JNl+1],
						0.0, 0.0, 1, 0) )*0.5,         // cell East side & Division by (OrderX+1)
		      (PolynomLineIntegration2(Node[CellIndex+1][JNl+1],Node[CellIndex  ][JNl+1],
					       0.0, 0.0, 0, 1) +                             // cell North side
		       PolynomLineIntegration2(Node[CellIndex  ][JNl+1],Node[CellIndex  ][JNl  ],
					       0.0, 0.0, 0, 1) +                             // cell West side
		       ExtendEast_BndSouthSpline.PolynomOrderIntegration(Node[CellIndex  ][JNl  ],Node[CellIndex+1][JNl  ],
									 Vector2D(0.0,0.0),15,0,1) + // cell South side
		       PolynomLineIntegration2(Node[CellIndex+1][JNl  ],Node[CellIndex+1][JNl+1],
					       0.0, 0.0, 0, 1)) ) /Cell[CellIndex][JCl].A;    // cell East side & Division by A
      
    case EXTEND_W_SOUTH_LEFT_SPLINE: // Use only the extension of South spline to West (left side)
      return Vector2D(( PolynomLineIntegration2(Node[CellIndex+1][JNl+1],Node[CellIndex  ][JNl+1],
						0.0, 0.0, 1, 0) +                           // cell North side
			PolynomLineIntegration2(Node[CellIndex  ][JNl+1],Node[CellIndex  ][JNl  ],
						0.0, 0.0, 1, 0) +                           // cell West side
			ExtendWest_BndSouthSpline.PolynomOrderIntegration(Node[CellIndex  ][JNl  ],Node[CellIndex+1][JNl  ],
									  Vector2D(0.0,0.0),15,1,0) + // cell South side
			PolynomLineIntegration2(Node[CellIndex+1][JNl  ],Node[CellIndex+1][JNl+1],
						0.0, 0.0, 1, 0) )*0.5,         // cell East side & Division by (OrderX+1)
		      (PolynomLineIntegration2(Node[CellIndex+1][JNl+1],Node[CellIndex  ][JNl+1],
					       0.0, 0.0, 0, 1) +                             // cell North side
		       PolynomLineIntegration2(Node[CellIndex  ][JNl+1],Node[CellIndex  ][JNl  ],
					       0.0, 0.0, 0, 1) +                             // cell West side
		       ExtendWest_BndSouthSpline.PolynomOrderIntegration(Node[CellIndex  ][JNl  ],Node[CellIndex+1][JNl  ],
									 Vector2D(0.0,0.0),15,0,1) + // cell South side
		       PolynomLineIntegration2(Node[CellIndex+1][JNl  ],Node[CellIndex+1][JNl+1],
					       0.0, 0.0, 0, 1)) ) /Cell[CellIndex][JCl].A;    // cell East side & Division by A

    case EXTEND_E_SOUTH_LEFT_SPLINE: // Use only the extension of South spline to East (left side)
      return Vector2D( (ExtendEast_BndSouthSpline.PolynomOrderIntegration(Node[CellIndex+1][JNl  ],Node[CellIndex  ][JNl  ],
									  Vector2D(0.0,0.0),15,1,0) +   // cell North side
			PolynomLineIntegration2(Node[CellIndex  ][JNl  ],Node[CellIndex  ][JNl-1],
						0.0, 0.0, 1, 0) +                                        // cell West side
			PolynomLineIntegration2(Node[CellIndex  ][JNl-1],Node[CellIndex+1][JNl-1],
						0.0, 0.0, 1, 0) +                                        // cell South side
			PolynomLineIntegration2(Node[CellIndex+1][JNl-1],Node[CellIndex+1][JNl  ],
						0.0, 0.0, 1, 0) )*0.5,                // cell East side & Division by (OrderX+1)
		       (ExtendEast_BndSouthSpline.PolynomOrderIntegration(Node[CellIndex+1][JNl  ],Node[CellIndex  ][JNl  ],
									  Vector2D(0.0,0.0),15,0,1) +   // cell North side
			PolynomLineIntegration2(Node[CellIndex  ][JNl  ],Node[CellIndex  ][JNl-1],
						0.0, 0.0, 0, 1) +                                        // cell West side
			PolynomLineIntegration2(Node[CellIndex  ][JNl-1],Node[CellIndex+1][JNl-1],
						0.0, 0.0, 0, 1) +                                        // cell South side
			PolynomLineIntegration2(Node[CellIndex+1][JNl-1],Node[CellIndex+1][JNl  ],
						0.0, 0.0, 0, 1) )) /Cell[CellIndex][JCl-1].A; // cell East side & Division by A  

    case EXTEND_N_EAST_RIGHT_SPLINE: // Use only the extension of East spline to North (right side)
      return Vector2D( (PolynomLineIntegration2(Node[INu  ][CellIndex+1],Node[INu-1][CellIndex+1],
						0.0, 0.0, 1, 0) +                                        // cell North side
			PolynomLineIntegration2(Node[INu-1][CellIndex+1],Node[INu-1][CellIndex  ],
						0.0, 0.0, 1, 0) +                                        // cell West side
			PolynomLineIntegration2(Node[INu-1][CellIndex  ],Node[INu  ][CellIndex  ],
						0.0, 0.0, 1, 0) +                                        // cell South side
			ExtendNorth_BndEastSpline.PolynomOrderIntegration(Node[INu  ][CellIndex  ],Node[INu  ][CellIndex+1],
									  Vector2D(0.0,0.0),15,1,0) ) *0.5,
		       // cell East side & Division by (OrderX+1)
		       (PolynomLineIntegration2(Node[INu  ][CellIndex+1],Node[INu-1][CellIndex+1],
						0.0, 0.0, 0, 1) +                                        // cell North side
			PolynomLineIntegration2(Node[INu-1][CellIndex+1],Node[INu-1][CellIndex  ],
						0.0, 0.0, 0, 1) +                                        // cell West side
			PolynomLineIntegration2(Node[INu-1][CellIndex  ],Node[INu  ][CellIndex  ],
						0.0, 0.0, 0, 1) +                                        // cell South side
			ExtendNorth_BndEastSpline.PolynomOrderIntegration(Node[INu  ][CellIndex  ],Node[INu  ][CellIndex+1],
									  Vector2D(0.0,0.0),15,0,1) )) /Cell[ICu][CellIndex].A;
      // cell East side & Division by A

    case EXTEND_S_EAST_RIGHT_SPLINE: // Use only the extension of East spline to South (right side)
      return Vector2D( (PolynomLineIntegration2(Node[INu+1][CellIndex+1],Node[INu  ][CellIndex+1],
						0.0, 0.0, 1, 0) +                                         // cell North side
			ExtendSouth_BndEastSpline.PolynomOrderIntegration(Node[INu  ][CellIndex+1],Node[INu  ][CellIndex  ],
									  Vector2D(0.0,0.0),15,1,0) +     // cell West side
			PolynomLineIntegration2(Node[INu  ][CellIndex  ],Node[INu+1][CellIndex  ],
						0.0, 0.0, 1, 0) +                                         // cell South side
			PolynomLineIntegration2(Node[INu+1][CellIndex  ],Node[INu+1][CellIndex+1],
						0.0, 0.0, 1, 0) ) *0.5, // cell East side & Division by (OrderX+1)
		       (PolynomLineIntegration2(Node[INu+1][CellIndex+1],Node[INu  ][CellIndex+1],
						0.0, 0.0, 0, 1) +                                         // cell North side
			ExtendSouth_BndEastSpline.PolynomOrderIntegration(Node[INu  ][CellIndex+1],Node[INu  ][CellIndex  ],
									  Vector2D(0.0,0.0),15,0,1) +     // cell West side
			PolynomLineIntegration2(Node[INu  ][CellIndex  ],Node[INu+1][CellIndex  ],
						0.0, 0.0, 0, 1) +                                         // cell South side
			PolynomLineIntegration2(Node[INu+1][CellIndex  ],Node[INu+1][CellIndex+1],
						0.0, 0.0, 0, 1) )) /Cell[ICu+1][CellIndex].A; // cell East side & Division by A

    case EXTEND_N_EAST_LEFT_SPLINE: // Use only the extension of East spline to North (left side)
      return Vector2D( (PolynomLineIntegration2(Node[INu+1][CellIndex+1],Node[INu  ][CellIndex+1],
						0.0, 0.0, 1, 0) +                                         // cell North side
			ExtendNorth_BndEastSpline.PolynomOrderIntegration(Node[INu  ][CellIndex+1],Node[INu  ][CellIndex  ],
									  Vector2D(0.0,0.0),15,1,0) +     // cell West side
			PolynomLineIntegration2(Node[INu  ][CellIndex  ],Node[INu+1][CellIndex  ],
						0.0, 0.0, 1, 0) +                                         // cell South side
			PolynomLineIntegration2(Node[INu+1][CellIndex  ],Node[INu+1][CellIndex+1],
						0.0, 0.0, 1, 0) ) *0.5, // cell East side & Division by (OrderX+1)
		       (PolynomLineIntegration2(Node[INu+1][CellIndex+1],Node[INu  ][CellIndex+1],
						0.0, 0.0, 0, 1) +                                         // cell North side
			ExtendNorth_BndEastSpline.PolynomOrderIntegration(Node[INu  ][CellIndex+1],Node[INu  ][CellIndex  ],
									  Vector2D(0.0,0.0),15,0,1) +     // cell West side
			PolynomLineIntegration2(Node[INu  ][CellIndex  ],Node[INu+1][CellIndex  ],
						0.0, 0.0, 0, 1) +                                         // cell South side
			PolynomLineIntegration2(Node[INu+1][CellIndex  ],Node[INu+1][CellIndex+1],
						0.0, 0.0, 0, 1) )) /Cell[ICu+1][CellIndex].A; // cell East side & Division by A

    case EXTEND_S_EAST_LEFT_SPLINE: // Use only the extension of East spline to South (left side)
      return Vector2D( (PolynomLineIntegration2(Node[INu  ][CellIndex+1],Node[INu-1][CellIndex+1],
						0.0, 0.0, 1, 0) +                                        // cell North side
			PolynomLineIntegration2(Node[INu-1][CellIndex+1],Node[INu-1][CellIndex  ],
						0.0, 0.0, 1, 0) +                                        // cell West side
			PolynomLineIntegration2(Node[INu-1][CellIndex  ],Node[INu  ][CellIndex  ],
						0.0, 0.0, 1, 0) +                                        // cell South side
			ExtendSouth_BndEastSpline.PolynomOrderIntegration(Node[INu  ][CellIndex  ],Node[INu  ][CellIndex+1],
									  Vector2D(0.0,0.0),15,1,0) ) *0.5,
		       // cell East side & Division by (OrderX+1)
		       (PolynomLineIntegration2(Node[INu  ][CellIndex+1],Node[INu-1][CellIndex+1],
						0.0, 0.0, 0, 1) +                                        // cell North side
			PolynomLineIntegration2(Node[INu-1][CellIndex+1],Node[INu-1][CellIndex  ],
						0.0, 0.0, 0, 1) +                                        // cell West side
			PolynomLineIntegration2(Node[INu-1][CellIndex  ],Node[INu  ][CellIndex  ],
						0.0, 0.0, 0, 1) +                                        // cell South side
			ExtendSouth_BndEastSpline.PolynomOrderIntegration(Node[INu  ][CellIndex  ],Node[INu  ][CellIndex+1],
									  Vector2D(0.0,0.0),15,0,1) )) /Cell[ICu][CellIndex].A;
      // cell East side & Division by A

    case EXTEND_N_WEST_RIGHT_SPLINE: // Use only the extension of West spline to North (right side)
      return Vector2D( (PolynomLineIntegration2(Node[INl  ][CellIndex+1],Node[INl-1][CellIndex+1],
						0.0, 0.0, 1, 0) +                                        // cell North side
			PolynomLineIntegration2(Node[INl-1][CellIndex+1],Node[INl-1][CellIndex  ],
						0.0, 0.0, 1, 0) +                                        // cell West side
			PolynomLineIntegration2(Node[INl-1][CellIndex  ],Node[INl  ][CellIndex  ],
						0.0, 0.0, 1, 0) +                                        // cell South side
			ExtendNorth_BndWestSpline.PolynomOrderIntegration(Node[INl  ][CellIndex  ],Node[INl  ][CellIndex+1],
									  Vector2D(0.0,0.0),15,1,0) ) *0.5,
		       // cell East side & Division by (OrderX+1)
		       (PolynomLineIntegration2(Node[INl  ][CellIndex+1],Node[INl-1][CellIndex+1],
						0.0, 0.0, 0, 1) +                                        // cell North side
			PolynomLineIntegration2(Node[INl-1][CellIndex+1],Node[INl-1][CellIndex  ],
						0.0, 0.0, 0, 1) +                                        // cell West side
			PolynomLineIntegration2(Node[INl-1][CellIndex  ],Node[INl  ][CellIndex  ],
						0.0, 0.0, 0, 1) +                                        // cell South side
			ExtendNorth_BndWestSpline.PolynomOrderIntegration(Node[INl  ][CellIndex  ],Node[INl  ][CellIndex+1],
									  Vector2D(0.0,0.0),15,0,1) )) /Cell[ICl-1][CellIndex].A;
      // cell East side & Division by A

    case EXTEND_S_WEST_RIGHT_SPLINE: // Use only the extension of West spline to South (right side)
      return Vector2D( (PolynomLineIntegration2(Node[INl+1][CellIndex+1],Node[INl  ][CellIndex+1],
						0.0, 0.0, 1, 0) +                                         // cell North side
			ExtendSouth_BndWestSpline.PolynomOrderIntegration(Node[INl  ][CellIndex+1],Node[INl  ][CellIndex  ],
									  Vector2D(0.0,0.0),15,1,0) +    // cell West side
			PolynomLineIntegration2(Node[INl  ][CellIndex  ],Node[INl+1][CellIndex  ],
						0.0, 0.0, 1, 0) +                                         // cell South side
			PolynomLineIntegration2(Node[INl+1][CellIndex  ],Node[INl+1][CellIndex+1],
						0.0, 0.0, 1, 0) ) *0.5, // cell East side & Division by (OrderX+1)
		       (PolynomLineIntegration2(Node[INl+1][CellIndex+1],Node[INl  ][CellIndex+1],
						0.0, 0.0, 0, 1) +                                         // cell North side
			ExtendSouth_BndWestSpline.PolynomOrderIntegration(Node[INl  ][CellIndex+1],Node[INl  ][CellIndex  ],
									  Vector2D(0.0,0.0),15,0,1) +    // cell West side
			PolynomLineIntegration2(Node[INl  ][CellIndex  ],Node[INl+1][CellIndex  ],
						0.0, 0.0, 0, 1) +                                         // cell South side
			PolynomLineIntegration2(Node[INl+1][CellIndex  ],Node[INl+1][CellIndex+1],
						0.0, 0.0, 0, 1) )) /Cell[ICl][CellIndex].A; // cell East side & Division by A

    case EXTEND_N_WEST_LEFT_SPLINE: // Use only the extension of West spline to North (left side)
      return Vector2D( (PolynomLineIntegration2(Node[INl+1][CellIndex+1],Node[INl  ][CellIndex+1],
						0.0, 0.0, 1, 0) +                                         // cell North side
			ExtendNorth_BndWestSpline.PolynomOrderIntegration(Node[INl  ][CellIndex+1],Node[INl  ][CellIndex  ],
									  Vector2D(0.0,0.0),15,1,0) +    // cell West side
			PolynomLineIntegration2(Node[INl  ][CellIndex  ],Node[INl+1][CellIndex  ],
						0.0, 0.0, 1, 0) +                                         // cell South side
			PolynomLineIntegration2(Node[INl+1][CellIndex  ],Node[INl+1][CellIndex+1],
						0.0, 0.0, 1, 0) ) *0.5, // cell East side & Division by (OrderX+1)
		       (PolynomLineIntegration2(Node[INl+1][CellIndex+1],Node[INl  ][CellIndex+1],
						0.0, 0.0, 0, 1) +                                         // cell North side
			ExtendNorth_BndWestSpline.PolynomOrderIntegration(Node[INl  ][CellIndex+1],Node[INl  ][CellIndex  ],
									  Vector2D(0.0,0.0),15,0,1) +    // cell West side
			PolynomLineIntegration2(Node[INl  ][CellIndex  ],Node[INl+1][CellIndex  ],
						0.0, 0.0, 0, 1) +                                         // cell South side
			PolynomLineIntegration2(Node[INl+1][CellIndex  ],Node[INl+1][CellIndex+1],
						0.0, 0.0, 0, 1) )) /Cell[ICl][CellIndex].A; // cell East side & Division by A

    case EXTEND_S_WEST_LEFT_SPLINE: // Use only the extension of West spline to South (left side)
      return Vector2D( (PolynomLineIntegration2(Node[INl  ][CellIndex+1],Node[INl-1][CellIndex+1],
						0.0, 0.0, 1, 0) +                                        // cell North side
			PolynomLineIntegration2(Node[INl-1][CellIndex+1],Node[INl-1][CellIndex  ],
						0.0, 0.0, 1, 0) +                                        // cell West side
			PolynomLineIntegration2(Node[INl-1][CellIndex  ],Node[INl  ][CellIndex  ],
						0.0, 0.0, 1, 0) +                                        // cell South side
			ExtendSouth_BndWestSpline.PolynomOrderIntegration(Node[INl  ][CellIndex  ],Node[INl  ][CellIndex+1],
									  Vector2D(0.0,0.0),15,1,0) ) *0.5,
		       // cell East side & Division by (OrderX+1)
		       (PolynomLineIntegration2(Node[INl  ][CellIndex+1],Node[INl-1][CellIndex+1],
						0.0, 0.0, 0, 1) +                                        // cell North side
			PolynomLineIntegration2(Node[INl-1][CellIndex+1],Node[INl-1][CellIndex  ],
						0.0, 0.0, 0, 1) +                                        // cell West side
			PolynomLineIntegration2(Node[INl-1][CellIndex  ],Node[INl  ][CellIndex  ],
						0.0, 0.0, 0, 1) +                                        // cell South side
			ExtendSouth_BndWestSpline.PolynomOrderIntegration(Node[INl  ][CellIndex  ],Node[INl  ][CellIndex+1],
									  Vector2D(0.0,0.0),15,0,1) )) /Cell[ICl-1][CellIndex].A;
      // cell East side & Division by A

    case CORNER_NORTH_EXTEND_N_WEST_SPLINES: // Use the North spline and the extension of West spline to North
      return Vector2D( (PolynomLineIntegration2(Node[INl+1][JNu+1],Node[INl  ][JNu+1],
						0.0, 0.0, 1, 0) +                                         // cell North side
			ExtendNorth_BndWestSpline.PolynomOrderIntegration(Node[INl  ][JNu+1],Node[INl  ][JNu],
									  Vector2D(0.0,0.0),15,1,0) +     // cell West side
			BndNorthSpline.PolynomOrderIntegration(Node[INl  ][JNu],Node[INl+1][JNu  ],
									  Vector2D(0.0,0.0),15,1,0) +     // cell South side
			PolynomLineIntegration2(Node[INl+1][JNu  ],Node[INl+1][JNu+1],
						0.0, 0.0, 1, 0) ) *0.5, // cell East side & Division by (OrderX+1)
		       (PolynomLineIntegration2(Node[INl+1][JNu+1],Node[INl  ][JNu+1],
						0.0, 0.0, 0, 1) +                                         // cell North side
			ExtendNorth_BndWestSpline.PolynomOrderIntegration(Node[INl  ][JNu+1],Node[INl  ][JNu],
									  Vector2D(0.0,0.0),15,0,1) +     // cell West side
			BndNorthSpline.PolynomOrderIntegration(Node[INl  ][JNu],Node[INl+1][JNu  ],
									  Vector2D(0.0,0.0),15,0,1) +     // cell South side
			PolynomLineIntegration2(Node[INl+1][JNu  ],Node[INl+1][JNu+1],
						0.0, 0.0, 0, 1) )) /Cell[ICl][JCu+1].A; // cell East side & Division by A

    case CORNER_EXTEND_N_WEST_EXTEND_W_NORTH_SPLINES: // Use the North extension of West and the West extension of North
      return Vector2D(( PolynomLineIntegration2(Node[INl  ][JNu+1],Node[INl-1][JNu+1],
						0.0, 0.0, 1, 0) +                           // cell North side
			PolynomLineIntegration2(Node[INl-1][JNu+1],Node[INl-1][JNu  ],
						0.0, 0.0, 1, 0) +                           // cell West side
			ExtendWest_BndNorthSpline.PolynomOrderIntegration(Node[INl-1][JNu  ],Node[INl][JNu  ],
									  Vector2D(0.0,0.0),15,1,0) + // cell South side
			ExtendNorth_BndWestSpline.PolynomOrderIntegration(Node[INl][JNu  ],Node[INl  ][JNu+1],
									  Vector2D(0.0,0.0),15,1,0) )*0.5,
		      // cell East side & Division by (OrderX+1)
		      (PolynomLineIntegration2(Node[INl  ][JNu+1],Node[INl-1][JNu+1],
					       0.0, 0.0, 0, 1) +                             // cell North side
		       PolynomLineIntegration2(Node[INl-1][JNu+1],Node[INl-1][JNu  ],
					       0.0, 0.0, 0, 1) +                             // cell West side
		       ExtendWest_BndNorthSpline.PolynomOrderIntegration(Node[INl-1][JNu  ],Node[INl][JNu  ],
									 Vector2D(0.0,0.0),15,0,1) + // cell South side
		       ExtendNorth_BndWestSpline.PolynomOrderIntegration(Node[INl][JNu  ],Node[INl  ][JNu+1],
									 Vector2D(0.0,0.0),15,0,1) )) /Cell[ICl-1][JCu+1].A;
      // cell East side & Division by A

    case CORNER_EXTEND_W_NORTH_WEST_SPLINES: // Use the West extension of North and the West spline
      return Vector2D( (ExtendWest_BndNorthSpline.PolynomOrderIntegration(Node[INl][JNu  ],Node[INl-1][JNu  ],
									  Vector2D(0.0,0.0),15,1,0) +     // cell North side
			PolynomLineIntegration2(Node[INl-1][JNu  ],Node[INl-1][JNu-1],
						0.0, 0.0, 1, 0) +                                        // cell West side
			PolynomLineIntegration2(Node[INl-1][JNu-1],Node[INl  ][JNu-1],
						0.0, 0.0, 1, 0) +                                        // cell South side
			BndWestSpline.PolynomOrderIntegration(Node[INl  ][JNu-1],Node[INl][JNu  ],
							      Vector2D(0.0,0.0),15,1,0) )*0.5, 
		       // cell East side & Division by (OrderX+1)
		       (ExtendWest_BndNorthSpline.PolynomOrderIntegration(Node[INl][JNu  ],Node[INl-1][JNu  ],
									  Vector2D(0.0,0.0),15,0,1) +    // cell North side
			PolynomLineIntegration2(Node[INl-1][JNu  ],Node[INl-1][JNu-1],
						0.0, 0.0, 0, 1) +                                        // cell West side
			PolynomLineIntegration2(Node[INl-1][JNu-1],Node[INl  ][JNu-1],
						0.0, 0.0, 0, 1) +                                        // cell South side
			BndWestSpline.PolynomOrderIntegration(Node[INl  ][JNu-1],Node[INl][JNu  ],
							      Vector2D(0.0,0.0),15,0,1) )) /Cell[ICl-1][JCu].A;
      // cell East side & Division by A  

    case CORNER_WEST_EXTEND_W_SOUTH_SPLINES: // Use the West spline and the West extension of South
      return Vector2D(( PolynomLineIntegration2(Node[INl  ][JNl+1],Node[INl-1][JNl+1],
						0.0, 0.0, 1, 0) +                           // cell North side
			PolynomLineIntegration2(Node[INl-1][JNl+1],Node[INl-1][JNl  ],
						0.0, 0.0, 1, 0) +                           // cell West side
			ExtendWest_BndSouthSpline.PolynomOrderIntegration(Node[INl-1][JNl  ],Node[INl][JNl  ],
									  Vector2D(0.0,0.0),15,1,0) + // cell South side
			BndWestSpline.PolynomOrderIntegration(Node[INl][JNl  ],Node[INl  ][JNl+1],
							      Vector2D(0.0,0.0),15,1,0) )*0.5,
		      // cell East side & Division by (OrderX+1)
		      (PolynomLineIntegration2(Node[INl  ][JNl+1],Node[INl-1][JNl+1],
					       0.0, 0.0, 0, 1) +                             // cell North side
		       PolynomLineIntegration2(Node[INl-1][JNl+1],Node[INl-1][JNl  ],
					       0.0, 0.0, 0, 1) +                             // cell West side
		       ExtendWest_BndSouthSpline.PolynomOrderIntegration(Node[INl-1][JNl  ],Node[INl][JNl  ],
									 Vector2D(0.0,0.0),15,0,1) + // cell South side
		       BndWestSpline.PolynomOrderIntegration(Node[INl][JNl  ],Node[INl  ][JNl+1],
							     Vector2D(0.0,0.0),15,0,1) )) /Cell[ICl-1][JCl].A;
      // cell East side & Division by A

    case CORNER_EXTEND_W_SOUTH_EXTEND_S_WEST_SPLINES: // Use the West extension of South and the South extension of West
      return Vector2D( (ExtendWest_BndSouthSpline.PolynomOrderIntegration(Node[INl][JNl  ],Node[INl-1][JNl  ],
									  Vector2D(0.0,0.0),15,1,0) +    // cell North side
			PolynomLineIntegration2(Node[INl-1][JNl  ],Node[INl-1][JNl-1],
						0.0, 0.0, 1, 0) +                                        // cell West side
			PolynomLineIntegration2(Node[INl-1][JNl-1],Node[INl  ][JNl-1],
						0.0, 0.0, 1, 0) +                                        // cell South side
			ExtendSouth_BndWestSpline.PolynomOrderIntegration(Node[INl  ][JNl-1],Node[INl][JNl  ],
									  Vector2D(0.0,0.0),15,1,0) )*0.5, 
		       // cell East side & Division by (OrderX+1)
		       (ExtendWest_BndSouthSpline.PolynomOrderIntegration(Node[INl][JNl  ],Node[INl-1][JNl  ],
									  Vector2D(0.0,0.0),15,0,1) +    // cell North side
			PolynomLineIntegration2(Node[INl-1][JNl  ],Node[INl-1][JNl-1],
						0.0, 0.0, 0, 1) +                                        // cell West side
			PolynomLineIntegration2(Node[INl-1][JNl-1],Node[INl  ][JNl-1],
						0.0, 0.0, 0, 1) +                                        // cell South side
			ExtendSouth_BndWestSpline.PolynomOrderIntegration(Node[INl  ][JNl-1],Node[INl][JNl  ],
									  Vector2D(0.0,0.0),15,0,1) )) /Cell[ICl-1][JCl-1].A;
      // cell East side & Division by A  

    case CORNER_EXTEND_S_WEST_SOUTH_SPLINES: // Use the South extension of West and the South spline
      return Vector2D( (BndSouthSpline.PolynomOrderIntegration(Node[INl+1][JNl],Node[INl][JNl],
							       Vector2D(0.0,0.0),15,1,0) +               // cell North side
			ExtendSouth_BndWestSpline.PolynomOrderIntegration(Node[INl][JNl],Node[INl  ][JNl-1],
									  Vector2D(0.0,0.0),15,1,0) +    // cell West side
			PolynomLineIntegration2(Node[INl  ][JNl-1],Node[INl+1][JNl-1],
						0.0, 0.0, 1, 0) +                                         // cell South side
			PolynomLineIntegration2(Node[INl+1][JNl-1],Node[INl+1][JNl],
						0.0, 0.0, 1, 0) ) *0.5, // cell East side & Division by (OrderX+1)
		       (BndSouthSpline.PolynomOrderIntegration(Node[INl+1][JNl],Node[INl][JNl],
							       Vector2D(0.0,0.0),15,0,1) +               // cell North side
			ExtendSouth_BndWestSpline.PolynomOrderIntegration(Node[INl][JNl],Node[INl  ][JNl-1],
									  Vector2D(0.0,0.0),15,0,1) +    // cell West side
			PolynomLineIntegration2(Node[INl  ][JNl-1],Node[INl+1][JNl-1],
						0.0, 0.0, 0, 1) +                                         // cell South side
			PolynomLineIntegration2(Node[INl+1][JNl-1],Node[INl+1][JNl],
						0.0, 0.0, 0, 1) )) /Cell[ICl][JCl-1].A; // cell East side & Division by A

    case CORNER_SOUTH_EXTEND_S_EAST_SPLINES: // Use the South spline and the South extension of East
      return Vector2D( (BndSouthSpline.PolynomOrderIntegration(Node[INu][JNl  ],Node[INu-1][JNl  ],
							       Vector2D(0.0,0.0),15,1,0) +               // cell North side
			PolynomLineIntegration2(Node[INu-1][JNl  ],Node[INu-1][JNl-1],
						0.0, 0.0, 1, 0) +                                        // cell West side
			PolynomLineIntegration2(Node[INu-1][JNl-1],Node[INu  ][JNl-1],
						0.0, 0.0, 1, 0) +                                        // cell South side
			ExtendSouth_BndEastSpline.PolynomOrderIntegration(Node[INu  ][JNl-1],Node[INu][JNl  ],
									  Vector2D(0.0,0.0),15,1,0) ) *0.5,
		       // cell East side & Division by (OrderX+1)
		       (BndSouthSpline.PolynomOrderIntegration(Node[INu][JNl  ],Node[INu-1][JNl  ],
							       Vector2D(0.0,0.0),15,0,1) +               // cell North side
			PolynomLineIntegration2(Node[INu-1][JNl  ],Node[INu-1][JNl-1],
						0.0, 0.0, 0, 1) +                                        // cell West side
			PolynomLineIntegration2(Node[INu-1][JNl-1],Node[INu  ][JNl-1],
						0.0, 0.0, 0, 1) +                                        // cell South side
			ExtendSouth_BndEastSpline.PolynomOrderIntegration(Node[INu  ][JNl-1],Node[INu][JNl  ],
									  Vector2D(0.0,0.0),15,0,1) )) /Cell[ICu][JCl-1].A;

    case CORNER_EXTEND_S_EAST_EXTEND_E_SOUTH_SPLINES: // Use the South extension of East and the East extension of South
      return Vector2D( (ExtendEast_BndSouthSpline.PolynomOrderIntegration(Node[INu+1][JNl  ],Node[INu  ][JNl],
									  Vector2D(0.0,0.0),15,1,0) +    // cell North side
			ExtendSouth_BndEastSpline.PolynomOrderIntegration(Node[INu  ][JNl],Node[INu  ][JNl-1],
									  Vector2D(0.0,0.0),15,1,0) +           // cell West side
			PolynomLineIntegration2(Node[INu  ][JNl-1],Node[INu+1][JNl-1],
						0.0, 0.0, 1, 0) +                                        // cell South side
			PolynomLineIntegration2(Node[INu+1][JNl-1],Node[INu+1][JNl  ],
						0.0, 0.0, 1, 0) )*0.5,                // cell East side & Division by (OrderX+1)
		       (ExtendEast_BndSouthSpline.PolynomOrderIntegration(Node[INu+1][JNl  ],Node[INu  ][JNl],
									  Vector2D(0.0,0.0),15,0,1) +    // cell North side
			ExtendSouth_BndEastSpline.PolynomOrderIntegration(Node[INu  ][JNl],Node[INu  ][JNl-1],
									  Vector2D(0.0,0.0),15,0,1) +           // cell West side
			PolynomLineIntegration2(Node[INu  ][JNl-1],Node[INu+1][JNl-1],
						0.0, 0.0, 0, 1) +                                        // cell South side
			PolynomLineIntegration2(Node[INu+1][JNl-1],Node[INu+1][JNl  ],
						0.0, 0.0, 0, 1) )) /Cell[ICu+1][JCl-1].A; // cell East side & Division by A  

    case CORNER_EXTEND_E_SOUTH_EAST_SPLINES: // Use the East extension of South and the East spline
      return Vector2D(( PolynomLineIntegration2(Node[INu+1][JNl+1],Node[INu  ][JNl+1],
						0.0, 0.0, 1, 0) +                           // cell North side
			BndEastSpline.PolynomOrderIntegration(Node[INu  ][JNl+1],Node[INu  ][JNl],
							      Vector2D(0.0,0.0),15,1,0) +   // cell West side
			ExtendEast_BndSouthSpline.PolynomOrderIntegration(Node[INu  ][JNl],Node[INu+1][JNl  ],
									  Vector2D(0.0,0.0),15,1,0) + // cell South side
			PolynomLineIntegration2(Node[INu+1][JNl  ],Node[INu+1][JNl+1],
						0.0, 0.0, 1, 0) )*0.5,         // cell East side & Division by (OrderX+1)
		      (PolynomLineIntegration2(Node[INu+1][JNl+1],Node[INu  ][JNl+1],
					       0.0, 0.0, 0, 1) +                             // cell North side
		       BndEastSpline.PolynomOrderIntegration(Node[INu  ][JNl+1],Node[INu  ][JNl],
							     Vector2D(0.0,0.0),15,0,1) +   // cell West side
		       ExtendEast_BndSouthSpline.PolynomOrderIntegration(Node[INu  ][JNl],Node[INu+1][JNl  ],
									 Vector2D(0.0,0.0),15,0,1) + // cell South side
		       PolynomLineIntegration2(Node[INu+1][JNl  ],Node[INu+1][JNl+1],
					       0.0, 0.0, 0, 1)) ) /Cell[ICu+1][JCl].A;    // cell East side & Division by A 

    case CORNER_EAST_EXTEND_E_NORTH_SPLINES: // Use the East spline and the East extension of North
      return Vector2D( (ExtendEast_BndNorthSpline.PolynomOrderIntegration(Node[INu+1][JNu  ],Node[INu  ][JNu],
									  Vector2D(0.0,0.0),15,1,0) +  // cell North side
			BndEastSpline.PolynomOrderIntegration(Node[INu  ][JNu],Node[INu  ][JNu-1],
							      Vector2D(0.0,0.0),15,1,0) +                       // cell West side
			PolynomLineIntegration2(Node[INu  ][JNu-1],Node[INu+1][JNu-1],
						0.0, 0.0, 1, 0) +                                        // cell South side
			PolynomLineIntegration2(Node[INu+1][JNu-1],Node[INu+1][JNu  ],
						0.0, 0.0, 1, 0) )*0.5,                // cell East side & Division by (OrderX+1)
		       (ExtendEast_BndNorthSpline.PolynomOrderIntegration(Node[INu+1][JNu  ],Node[INu  ][JNu],
									  Vector2D(0.0,0.0),15,0,1) +    // cell North side
			BndEastSpline.PolynomOrderIntegration(Node[INu  ][JNu],Node[INu  ][JNu-1],
							      Vector2D(0.0,0.0),15,0,1) +                // cell West side
			PolynomLineIntegration2(Node[INu  ][JNu-1],Node[INu+1][JNu-1],
						0.0, 0.0, 0, 1) +                                        // cell South side
			PolynomLineIntegration2(Node[INu+1][JNu-1],Node[INu+1][JNu  ],
						0.0, 0.0, 0, 1) )) /Cell[ICu+1][JCu].A; // cell East side & Division by A  

    case CORNER_EXTEND_E_NORTH_EXTEND_N_EAST_SPLINES: // Use the East extention of North and North extension of East
      return Vector2D( (PolynomLineIntegration2(Node[INu+1][JNu+1],Node[INu  ][JNu+1],
						0.0, 0.0, 1, 0) +                                         // cell North side
			ExtendNorth_BndEastSpline.PolynomOrderIntegration(Node[INu  ][JNu+1],Node[INu  ][JNu],
									  Vector2D(0.0,0.0),15,1,0) +    // cell West side
			ExtendEast_BndNorthSpline.PolynomOrderIntegration(Node[INu  ][JNu],Node[INu+1][JNu  ],
									  Vector2D(0.0,0.0),15,1,0) +    // cell South side
			PolynomLineIntegration2(Node[INu+1][JNu  ],Node[INu+1][JNu+1],
						0.0, 0.0, 1, 0) ) *0.5, // cell East side & Division by (OrderX+1)
		       (PolynomLineIntegration2(Node[INu+1][JNu+1],Node[INu  ][JNu+1],
						0.0, 0.0, 0, 1) +                                         // cell North side
			ExtendNorth_BndEastSpline.PolynomOrderIntegration(Node[INu  ][JNu+1],Node[INu  ][JNu],
									  Vector2D(0.0,0.0),15,0,1) +    // cell West side
			ExtendEast_BndNorthSpline.PolynomOrderIntegration(Node[INu  ][JNu],Node[INu+1][JNu  ],
									  Vector2D(0.0,0.0),15,0,1) +    // cell South side
			PolynomLineIntegration2(Node[INu+1][JNu  ],Node[INu+1][JNu+1],
						0.0, 0.0, 0, 1) )) /Cell[ICu+1][JCu+1].A; // cell East side & Division by A

    case CORNER_EXTEND_N_EAST_NORTH_SPLINES: // Use the North extension of East and the North spline
      return Vector2D( (PolynomLineIntegration2(Node[INu  ][JNu+1],Node[INu-1][JNu+1],
						0.0, 0.0, 1, 0) +                                        // cell North side
			PolynomLineIntegration2(Node[INu-1][JNu+1],Node[INu-1][JNu  ],
						0.0, 0.0, 1, 0) +                                        // cell West side
			BndNorthSpline.PolynomOrderIntegration(Node[INu-1][JNu  ],Node[INu][JNu  ],
							       Vector2D(0.0,0.0),15,1,0) +               // cell South side
			ExtendNorth_BndEastSpline.PolynomOrderIntegration(Node[INu][JNu  ],Node[INu  ][JNu+1],
									  Vector2D(0.0,0.0),15,1,0) ) *0.5,
		       // cell East side & Division by (OrderX+1)
		       (PolynomLineIntegration2(Node[INu  ][JNu+1],Node[INu-1][JNu+1],
						0.0, 0.0, 0, 1) +                                        // cell North side
			PolynomLineIntegration2(Node[INu-1][JNu+1],Node[INu-1][JNu  ],
						0.0, 0.0, 0, 1) +                                        // cell West side
			BndNorthSpline.PolynomOrderIntegration(Node[INu-1][JNu  ],Node[INu][JNu  ],
							       Vector2D(0.0,0.0),15,0,1) +               // cell South side
			ExtendNorth_BndEastSpline.PolynomOrderIntegration(Node[INu][JNu  ],Node[INu  ][JNu+1],
									  Vector2D(0.0,0.0),15,0,1) )) /Cell[ICu][JCu+1].A;

    default:
      return Vector2D(0.0);
    } // endswitch

  }// endif
}

/*!
 * Compute geometric moments with respect to the cell centroid 
 * (xCC,yCC) of cell (ii,jj) up to 4th order.                  
 * These moments are integrals over the cell domain of
 * polynomial functions with the form ((x-xCC)^n)*((y-yCC)^m) 
 * and divided by aria A.
 * This subroutine is for grid cells with straight edges.
 *
 * \todo Derive analytic expressions for the quartic moments to speed up the code.
 */
void Grid2D_Quad_Block_HO::ComputeGeometricCoefficients(const int &ii, const int &jj){

  // Compute the coefficients of the bilinear interpolation for transforming
  // a quadrilatral cell into a unit lenght square

  /* determine the coefficients of the transformation */
  double a0C = nodeSW(ii,jj).x() - Cell[ii][jj].Xc.x;
  double a1  = nodeSE(ii,jj).x() - nodeSW(ii,jj).x();
  double a2  = nodeNW(ii,jj).x() - nodeSW(ii,jj).x();
  double a3  = nodeSW(ii,jj).x() + nodeNE(ii,jj).x() - nodeNW(ii,jj).x() - nodeSE(ii,jj).x();
  
  double b0C = nodeSW(ii,jj).y() - Cell[ii][jj].Xc.y;
  double b1  = nodeSE(ii,jj).y() - nodeSW(ii,jj).y();
  double b2  = nodeNW(ii,jj).y() - nodeSW(ii,jj).y();
  double b3  = nodeSW(ii,jj).y() + nodeNE(ii,jj).y() - nodeNW(ii,jj).y() - nodeSE(ii,jj).y();

  double J0  = a1*b2 - a2*b1;
  double J1  = a1*b3 - a3*b1;
  double J2  = b2*a3 - a2*b3;

  double DummyParam(0);

  /* Create the polynomial function for the cell */
  GeneralizedPolynomialFunctionOfTwoVariables Polynom(0,0,XCellCentroid(ii,jj),YCellCentroid(ii,jj));

  switch(HighestReconstructionOrder){

  case 4:  // Quartic moments

    // coefficient (4,0)
    Polynom.ChangePowersTo(4,0);
    Cell[ii][jj].GeomCoeffValue(4,0) = Integration.IntegratePolynomialOverCell(ii,jj,Polynom,DummyParam) /Cell[ii][jj].A;

    // coefficient (3,1)
    Polynom.ChangePowersTo(3,1);
    Cell[ii][jj].GeomCoeffValue(3,1) = Integration.IntegratePolynomialOverCell(ii,jj,Polynom,DummyParam) /Cell[ii][jj].A;

    // coefficient (2,2)
    Polynom.ChangePowersTo(2,2);
    Cell[ii][jj].GeomCoeffValue(2,2) = Integration.IntegratePolynomialOverCell(ii,jj,Polynom,DummyParam) /Cell[ii][jj].A;

    // coefficient (1,3)
    Polynom.ChangePowersTo(1,3);
    Cell[ii][jj].GeomCoeffValue(1,3) = Integration.IntegratePolynomialOverCell(ii,jj,Polynom,DummyParam) /Cell[ii][jj].A;

    // coefficient (0,4)
    Polynom.ChangePowersTo(0,4);
    Cell[ii][jj].GeomCoeffValue(0,4) = Integration.IntegratePolynomialOverCell(ii,jj,Polynom,DummyParam) /Cell[ii][jj].A;

  case 3:   // Cubic moments
    // coefficient (3,0)
    Cell[ii][jj].GeomCoeffValue(3,0) = ( ((240*J0 + 120*J1 + 120*J2)*a0C*a0C*a0C + (360*J0 + 240*J1 + 180*J2)*a0C*a0C*a1  + 
					  (360*J0 + 180*J1 + 240*J2)*a0C*a0C*a2  + (180*J0 + 120*J1 + 120*J2)*a0C*a0C*a3  + 
					  (240*J0 + 180*J1 + 120*J2)*a0C*a1*a1  + (360*J0 + 240*J1 + 240*J2)*a0C*a1*a2  + 
					  (240*J0 + 180*J1 + 160*J2)*a0C*a1*a3  + (240*J0 + 120*J1 + 180*J2)*a0C*a2*a2  + 
					  (240*J0 + 160*J1 + 180*J2)*a0C*a2*a3  + (80*J0 + 60*J1 + 60*J2)*a0C*a3*a3  + 
					  (60*J0 + 48*J1 + 30*J2)*a1*a1*a1  + (120*J0 + 90*J1 + 80*J2)*a1*a1*a2  + 
					  (90*J0 + 72*J1 + 60*J2)*a1*a1*a3  + (120*J0 + 80*J1 + 90*J2)*a1*a2*a2  + 
					  (160*J0 + 120*J1 + 120*J2)*a1*a2*a3  + (60*J0 + 48*J1 + 45*J2)*a1*a3*a3  + 
					  (60*J0 + 30*J1 + 48*J2)*a2*a2*a2 + (90*J0 + 60*J1 + 72*J2)*a2*a2*a3  + 
					  (60*J0 + 45*J1 + 48*J2)*a2*a3*a3  + (15*J0 + 12*J1 + 12*J2)*a3*a3*a3) 
					 /Cell[ii][jj].A /0.240e3 );

    // coefficient (2,1)
    Cell[ii][jj].GeomCoeffValue(2,1) = ( (((360*J0 + 180*J1 + 240*J2)*b2*a0C*a0C)  + ((360*J0 + 240*J1 + 240*J2)*b2*a0C*a1)  + 
					  ((480*J0 + 240*J1 + 360*J2)*b2*a0C*a2)  + ((240*J0 + 160*J1 + 180*J2)*b2*a0C*a3)  + 
					  ((120*J0 + 90*J1 + 80*J2)*b2*a1*a1)  + ((240*J0 + 160*J1 + 180*J2)*b2*a1*a2)  + 
					  ((160*J0 + 120*J1 + 120*J2)*b2*a1*a3)  + ((180*J0 + 90*J1 + 144*J2)*b2*a2*a2)  + 
					  ((180*J0 + 120*J1 + 144*J2)*b2*a2*a3)  + ((60*J0 + 45*J1 + 48*J2)*b2*a3*a3)  + 
					  ((720*J0 + 360*J1 + 360*J2)*a0C*a0C*b0C)  + ((360*J0 + 240*J1 + 180*J2)*a0C*a0C*b1)  + 
					  ((180*J0 + 120*J1 + 120*J2)*a0C*a0C*b3)  + ((720*J0 + 480*J1 + 360*J2)*a0C*a1*b0C)  + 
					  ((480*J0 + 360*J1 + 240*J2)*a0C*a1*b1)  + ((240*J0 + 180*J1 + 160*J2)*a0C*a1*b3)  + 
					  ((720*J0 + 360*J1 + 480*J2)*a0C*a2*b0C)  + ((360*J0 + 240*J1 + 240*J2)*a0C*a2*b1)  + 
					  ((240*J0 + 160*J1 + 180*J2)*a0C*a2*b3)  + ((360*J0 + 240*J1 + 240*J2)*a0C*a3*b0C)  + 
					  ((240*J0 + 180*J1 + 160*J2)*a0C*a3*b1)  + ((160*J0 + 120*J1 + 120*J2)*a0C*a3*b3)  + 
					  ((240*J0 + 180*J1 + 120*J2)*a1*a1*b0C)  + ((180*J0 + 144*J1 + 90*J2)*a1*a1*b1)  + 
					  ((90*J0 + 72*J1 + 60*J2)*a1*a1*b3)  + ((360*J0 + 240*J1 + 240*J2)*a1*a2*b0C)  + 
					  ((240*J0 + 180*J1 + 160*J2)*a1*a2*b1)  + ((160*J0 + 120*J1 + 120*J2)*a1*a2*b3)  + 
					  ((240*J0 + 180*J1 + 160*J2)*a1*a3*b0C)  + ((180*J0 + 144*J1 + 120*J2)*a1*a3*b1)  + 
					  ((120*J0 + 96*J1 + 90*J2)*a1*a3*b3)  + ((240*J0 + 120*J1 + 180*J2)*a2*a2*b0C)  + 
					  ((120*J0 + 80*J1 + 90*J2)*a2*a2*b1)  + ((90*J0 + 60*J1 + 72*J2)*a2*a2*b3)  + 
					  ((240*J0 + 160*J1 + 180*J2)*a2*a3*b0C)  + ((160*J0 + 120*J1 + 120*J2)*a2*a3*b1)  + 
					  ((120*J0 + 90*J1 + 96*J2)*a2*a3*b3)  + ((80*J0 + 60*J1 + 60*J2)*a3*a3*b0C)  + 
					  ((60*J0 + 48*J1 + 45*J2)*a3*a3*b1) + (45*J0 + 36*J1 + 36*J2)*a3*a3*b3)
					 /Cell[ii][jj].A /0.720e3 );

    // coefficient (1,2)
    Cell[ii][jj].GeomCoeffValue(1,2) = ( (((240*J0 + 120*J1 + 180*J2)*b2*b2*a0C) + ((120*J0 + 80*J1 + 90*J2)*b2*b2*a1) + 
					  ((180*J0 + 90*J1 + 144*J2)*b2*b2*a2) + ((90*J0 + 60*J1 + 72*J2)*b2*b2*a3) + 
					  ((720*J0 + 360*J1 + 480*J2)*b2*a0C*b0C) + ((360*J0 + 240*J1 + 240*J2)*b2*a0C*b1) + 
					  ((240*J0 + 160*J1 + 180*J2)*b2*a0C*b3) + ((360*J0 + 240*J1 + 240*J2)*b2*a1*b0C) + 
					  ((240*J0 + 180*J1 + 160*J2)*b2*a1*b1) + ((160*J0 + 120*J1 + 120*J2)*b2*a1*b3) + 
					  ((480*J0 + 240*J1 + 360*J2)*b2*a2*b0C) + ((240*J0 + 160*J1 + 180*J2)*b2*a2*b1) + 
					  ((180*J0 + 120*J1 + 144*J2)*b2*a2*b3) + ((240*J0 + 160*J1 + 180*J2)*b2*a3*b0C) + 
					  ((160*J0 + 120*J1 + 120*J2)*b2*a3*b1) + ((120*J0 + 90*J1 + 96*J2)*b2*a3*b3) + 
					  ((720*J0 + 360*J1 + 360*J2)*a0C*b0C*b0C) + ((720*J0 + 480*J1 + 360*J2)*a0C*b0C*b1) + 
					  ((360*J0 + 240*J1 + 240*J2)*a0C*b0C*b3) + ((240*J0 + 180*J1 + 120*J2)*a0C*b1*b1) + 
					  ((240*J0 + 180*J1 + 160*J2)*a0C*b1*b3) + ((80*J0 + 60*J1 + 60*J2)*a0C*b3*b3) + 
					  ((360*J0 + 240*J1 + 180*J2)*a1*b0C*b0C) + ((480*J0 + 360*J1 + 240*J2)*a1*b0C*b1) + 
					  ((240*J0 + 180*J1 + 160*J2)*a1*b0C*b3) + ((180*J0 + 144*J1 + 90*J2)*a1*b1*b1) + 
					  ((180*J0 + 144*J1 + 120*J2)*a1*b1*b3) + ((60*J0 + 48*J1 + 45*J2)*a1*b3*b3) + 
					  ((360*J0 + 180*J1 + 240*J2)*a2*b0C*b0C) + ((360*J0 + 240*J1 + 240*J2)*a2*b0C*b1) + 
					  ((240*J0 + 160*J1 + 180*J2)*a2*b0C*b3) + ((120*J0 + 90*J1 + 80*J2)*a2*b1*b1) + 
					  ((160*J0 + 120*J1 + 120*J2)*a2*b1*b3) + ((60*J0 + 45*J1 + 48*J2)*a2*b3*b3) + 
					  ((180*J0 + 120*J1 + 120*J2)*a3*b0C*b0C) + ((240*J0 + 180*J1 + 160*J2)*a3*b0C*b1) + 
					  ((160*J0 + 120*J1 + 120*J2)*a3*b0C*b3) + ((90*J0 + 72*J1 + 60*J2)*a3*b1*b1) + 
					  ((120*J0 + 96*J1 + 90*J2)*a3*b1*b3) + (45*J0 + 36*J1 + 36*J2)*a3*b3*b3) 
					 /Cell[ii][jj].A /0.720e3 );

    // coefficient (0,3)
    Cell[ii][jj].GeomCoeffValue(0,3) = ( ((60*J0 + 30*J1 + 48*J2)*b2*b2*b2 + (240*J0 + 120*J1 + 180*J2)*b2*b2*b0C +
					  (120*J0 + 80*J1 + 90*J2)*b2*b2*b1 + (90*J0 + 60*J1 + 72*J2)*b2*b2*b3 + 
					  (360*J0 + 180*J1 + 240*J2)*b2*b0C*b0C + (360*J0 + 240*J1 + 240*J2)*b2*b0C*b1 + 
					  (240*J0 + 160*J1 + 180*J2)*b2*b0C*b3 + (120*J0 + 90*J1 + 80*J2)*b2*b1*b1 +
					  (160*J0 + 120*J1 + 120*J2)*b2*b1*b3 + (60*J0 + 45*J1 + 48*J2)*b2*b3*b3 + 
					  (240*J0 + 120*J1 + 120*J2)*b0C*b0C*b0C + (360*J0 + 240*J1 + 180*J2)*b0C*b0C*b1 + 
					  (180*J0 + 120*J1 + 120*J2)*b0C*b0C*b3 + (240*J0 + 180*J1 + 120*J2)*b0C*b1*b1 +
					  (240*J0 + 180*J1 + 160*J2)*b0C*b1*b3 + (80*J0 + 60*J1 + 60*J2)*b0C*b3*b3 + 
					  (60*J0 + 48*J1 + 30*J2)*b1*b1*b1 + (90*J0 + 72*J1 + 60*J2)*b1*b1*b3 + 
					  (60*J0 + 48*J1 + 45*J2)*b1*b3*b3 + (15*J0 + 12*J1 + 12*J2)*b3*b3*b3) 
					 /Cell[ii][jj].A /0.240e3 );

  case 2:   // Quadratic moments
    // coefficient (2,0)
    Cell[ii][jj].GeomCoeffValue(2,0) = ( (((36*J0 + 18*J1 + 18*J2)*a0C*a0C) + ((36*J0 + 24*J1 + 18*J2)*a0C*a1) + 
					  ((36*J0 + 18*J1 + 24*J2)*a0C*a2) + ((18*J0 + 12*J1 + 12*J2)*a0C*a3) +
					  ((12*J0 + 9*J1 + 6*J2)*a1*a1) + ((18*J0 + 12*J1 + 12*J2)*a1*a2) + 
					  ((12*J0 + 9*J1 + 8*J2)*a1*a3) + ((12*J0 + 6*J1 + 9*J2)*a2*a2) + 
					  ((12*J0 + 8*J1 + 9*J2)*a2*a3) + ((4*J0 + 3*J1 + 3*J2)*a3*a3)) 
					 /Cell[ii][jj].A /0.36e2 );

    // coefficient (0,2)
    Cell[ii][jj].GeomCoeffValue(0,2) = ( (((12*J0 + 6*J1 + 9*J2)*b2*b2) + ((36*J0 + 18*J1 + 24*J2)*b2*b0C) +
					  ((18*J0 + 12*J1 + 12*J2)*b2*b1) + ((12*J0 + 8*J1 + 9*J2)*b2*b3) + 
					  ((36*J0 + 18*J1 + 18*J2)*b0C*b0C) + ((36*J0 + 24*J1 + 18*J2)*b0C*b1) + 
					  ((18*J0 + 12*J1 + 12*J2)*b0C*b3) + ((12*J0 + 9*J1 + 6*J2)*b1*b1) + 
					  ((12*J0 + 9*J1 + 8*J2)*b1*b3) + ((4*J0 + 3*J1 + 3*J2)*b3*b3)) 
					 /Cell[ii][jj].A /0.36e2 );

    // coefficient (1,1)
    Cell[ii][jj].GeomCoeffValue(1,1) = ( (((36*J0 + 18*J1 + 24*J2)*b2*a0C) + ((18*J0 + 12*J1 + 12*J2)*b2*a1) + 
					  ((24*J0 + 12*J1 + 18*J2)*b2*a2) + ((12*J0 + 8*J1 + 9*J2)*b2*a3) + 
					  ((72*J0 + 36*J1 + 36*J2)*a0C*b0C) + ((36*J0 + 24*J1 + 18*J2)*a0C*b1) +
					  ((18*J0 + 12*J1 + 12*J2)*a0C*b3) + ((36*J0 + 24*J1 + 18*J2)*a1*b0C) + 
					  ((24*J0 + 18*J1 + 12*J2)*a1*b1) + ((12*J0 + 9*J1 + 8*J2)*a1*b3) +
					  ((36*J0 + 18*J1 + 24*J2)*a2*b0C) + ((18*J0 + 12*J1 + 12*J2)*a2*b1) +
					  ((12*J0 + 8*J1 + 9*J2)*a2*b3) + ((18*J0 + 12*J1 + 12*J2)*a3*b0C) + 
					  ((12*J0 + 9*J1 + 8*J2)*a3*b1) + ((8*J0 + 6*J1 + 6*J2)*a3*b3)) 
					 / Cell[ii][jj].A /0.72e2 );

  case 1:   // Linear moments
    Cell[ii][jj].GeomCoeffValue(1,0) = 0.0;
    Cell[ii][jj].GeomCoeffValue(0,1) = 0.0;
    
  case 0:   // Zero moment 
    Cell[ii][jj].GeomCoeffValue(0,0) = 1.0;
    break;

  default: 
    throw runtime_error("Grid2D_Quad_Block_HO::ComputeGeometricCoefficients() ERROR! Geometric moments for this order have not been set up.");
  }
}

/*!
 * Compute geometric moments with respect to the cell centroid 
 * (xCC,yCC) of cell (CellIndexI,CellIndexJ) up to 4th order.                  
 * These moments are integrals over the cell domain of
 * polynomial functions with the form ((x-xCC)^n)*((y-yCC)^m) 
 * and divided by aria A.
 * This subroutine is for interior boundary cells.
 */
void Grid2D_Quad_Block_HO::ComputeGeometricCoefficients_CurvedBoundaries(const int &CellIndexI, const int &CellIndexJ,
									 const int &Boundary){

  /* Moments for quadrilateral cells with curved boundaries */

  switch(HighestReconstructionOrder){

  case 4:  // Quartic moments
    
    // coefficient (4,0)
    Cell[CellIndexI][CellIndexJ].GeomCoeffValue(4,0) = GeometricMoment_CurvedBoundaries(CellIndexI, CellIndexJ, Boundary, 4, 0);
    
    // coefficient (3,1)
    Cell[CellIndexI][CellIndexJ].GeomCoeffValue(3,1) = GeometricMoment_CurvedBoundaries(CellIndexI, CellIndexJ, Boundary, 3, 1);
    
    // coefficient (2,2)
    Cell[CellIndexI][CellIndexJ].GeomCoeffValue(2,2) = GeometricMoment_CurvedBoundaries(CellIndexI, CellIndexJ, Boundary, 2, 2);
    
    // coefficient (1,3)
    Cell[CellIndexI][CellIndexJ].GeomCoeffValue(1,3) = GeometricMoment_CurvedBoundaries(CellIndexI, CellIndexJ, Boundary, 1, 3);
    
    // coefficient (0,4)
    Cell[CellIndexI][CellIndexJ].GeomCoeffValue(0,4) = GeometricMoment_CurvedBoundaries(CellIndexI, CellIndexJ, Boundary, 0, 4);

  case 3:   // Cubic moments
    // coefficient (3,0)
    Cell[CellIndexI][CellIndexJ].GeomCoeffValue(3,0) = GeometricMoment_CurvedBoundaries(CellIndexI, CellIndexJ, Boundary, 3, 0);

    // coefficient (2,1)
    Cell[CellIndexI][CellIndexJ].GeomCoeffValue(2,1) = GeometricMoment_CurvedBoundaries(CellIndexI, CellIndexJ, Boundary, 2, 1);

    // coefficient (1,2)
    Cell[CellIndexI][CellIndexJ].GeomCoeffValue(1,2) = GeometricMoment_CurvedBoundaries(CellIndexI, CellIndexJ, Boundary, 1, 2);
    // coefficient (0,3)
    Cell[CellIndexI][CellIndexJ].GeomCoeffValue(0,3) = GeometricMoment_CurvedBoundaries(CellIndexI, CellIndexJ, Boundary, 0, 3);

  case 2:   // Quadratic moments
    // coefficient (2,0)
    Cell[CellIndexI][CellIndexJ].GeomCoeffValue(2,0) = GeometricMoment_CurvedBoundaries(CellIndexI, CellIndexJ, Boundary, 2, 0);

    // coefficient (0,2)
    Cell[CellIndexI][CellIndexJ].GeomCoeffValue(0,2) = GeometricMoment_CurvedBoundaries(CellIndexI, CellIndexJ, Boundary, 0, 2);

    // coefficient (1,1)
    Cell[CellIndexI][CellIndexJ].GeomCoeffValue(1,1) = GeometricMoment_CurvedBoundaries(CellIndexI, CellIndexJ, Boundary, 1, 1);

  case 1:   // Linear moments
    Cell[CellIndexI][CellIndexJ].GeomCoeffValue(1,0) = 0.0;
    Cell[CellIndexI][CellIndexJ].GeomCoeffValue(0,1) = 0.0;
    
  case 0:   // Zero moment 
    Cell[CellIndexI][CellIndexJ].GeomCoeffValue(0,0) = 1.0;
    
  }
}

/*!
 * Compute geometric moments with respect to the cell centroid 
 * (xCC,yCC) of cell (CellIndexI,CellIndexJ) up to 4th order.                  
 * These moments are integrals over the cell domain of
 * polynomial functions with the form ((x-xCC)^n)*((y-yCC)^m) 
 * and divided by aria A.
 * This subroutine is for ghost boundary cells.
 */
void Grid2D_Quad_Block_HO::ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(const int &CellIndexI, const int &CellIndexJ,
										   const int &Boundary){

  /* Moments for quadrilateral cells with curved boundaries */

  switch(HighestReconstructionOrder){

  case 4:    // Quartic moments
    // coefficient (4,0)
    Cell[CellIndexI][CellIndexJ].GeomCoeffValue(4,0) = GeometricMoment_GhostCell_CurvedBoundaries(CellIndexI, CellIndexJ,
												  Boundary, 4, 0);
    // coefficient (3,1)
    Cell[CellIndexI][CellIndexJ].GeomCoeffValue(3,1) = GeometricMoment_GhostCell_CurvedBoundaries(CellIndexI, CellIndexJ,
												  Boundary, 3, 1);
    // coefficient (2,2)
    Cell[CellIndexI][CellIndexJ].GeomCoeffValue(2,2) = GeometricMoment_GhostCell_CurvedBoundaries(CellIndexI, CellIndexJ,
												  Boundary, 2, 2);
    // coefficient (1,3)
    Cell[CellIndexI][CellIndexJ].GeomCoeffValue(1,3) = GeometricMoment_GhostCell_CurvedBoundaries(CellIndexI, CellIndexJ,
												  Boundary, 1, 3);
    // coefficient (0,4)
    Cell[CellIndexI][CellIndexJ].GeomCoeffValue(0,4) = GeometricMoment_GhostCell_CurvedBoundaries(CellIndexI, CellIndexJ,
												  Boundary, 0, 4);

  case 3:   // Cubic moments
    // coefficient (3,0)
    Cell[CellIndexI][CellIndexJ].GeomCoeffValue(3,0) = GeometricMoment_GhostCell_CurvedBoundaries(CellIndexI, CellIndexJ,
												  Boundary, 3, 0);
    // coefficient (2,1)
    Cell[CellIndexI][CellIndexJ].GeomCoeffValue(2,1) = GeometricMoment_GhostCell_CurvedBoundaries(CellIndexI, CellIndexJ,
												  Boundary, 2, 1);
    // coefficient (1,2)
    Cell[CellIndexI][CellIndexJ].GeomCoeffValue(1,2) = GeometricMoment_GhostCell_CurvedBoundaries(CellIndexI, CellIndexJ,
												  Boundary, 1, 2);
    // coefficient (0,3)
    Cell[CellIndexI][CellIndexJ].GeomCoeffValue(0,3) = GeometricMoment_GhostCell_CurvedBoundaries(CellIndexI, CellIndexJ,
												  Boundary, 0, 3);

  case 2:   // Quadratic moments
    // coefficient (2,0)
    Cell[CellIndexI][CellIndexJ].GeomCoeffValue(2,0) = GeometricMoment_GhostCell_CurvedBoundaries(CellIndexI, CellIndexJ,
												  Boundary, 2, 0);
    // coefficient (0,2)
    Cell[CellIndexI][CellIndexJ].GeomCoeffValue(0,2) = GeometricMoment_GhostCell_CurvedBoundaries(CellIndexI, CellIndexJ,
												  Boundary, 0, 2);
    // coefficient (1,1)
    Cell[CellIndexI][CellIndexJ].GeomCoeffValue(1,1) = GeometricMoment_GhostCell_CurvedBoundaries(CellIndexI, CellIndexJ,
												  Boundary, 1, 1);

  case 1:   // Linear moments
    Cell[CellIndexI][CellIndexJ].GeomCoeffValue(1,0) = 0.0;
    Cell[CellIndexI][CellIndexJ].GeomCoeffValue(0,1) = 0.0;
    
  case 0:   // Zero moment 
    Cell[CellIndexI][CellIndexJ].GeomCoeffValue(0,0) = 1.0;
    
  }
}

/*!
 * Compute a particular geometric moment for an interior cell.
 *
 * \param CellIndexI i-index of the cell
 * \param CellIndexJ j-index of the cell
 * \param Boundary indicates which cell boundaries are curved
 * \param OrderX the x-power of the polynomial function which is integrated
 * \param OrderY the y-power of the polynomial function which is integrated
 */
double Grid2D_Quad_Block_HO::GeometricMoment_CurvedBoundaries(const int &CellIndexI, const int &CellIndexJ,
							      const int &Boundary,
							      const int &OrderX, const int &OrderY){

  // Obs. The sides of the cell that are not curved are treated as line segments and therefore
  //      the line integral is computed exactly based on the Nodes.
  // The edges are considered in counterclockwise order (i.e. right-hand rule applied)

  if (Gauss_Quad_Curvilinear_Integration) {

    // Use the SplineInfo variables to integrate along curved edges.

    switch(Boundary){

    case NORTH_SPLINE:          // Use only the North Spline -> (iCell,jCell)=(CellIndexI,JCu)
      return ( (BndNorthSplineInfo[CellIndexI].IntegratePolynomialTerm(Cell[CellIndexI][JCu].Xc,
								       OrderX, OrderY) +                         // cell North side
		PolynomLineIntegration2(Node[CellIndexI  ][JNu].x(),Node[CellIndexI  ][JNu].y(),
					Node[CellIndexI  ][JCu].x(),Node[CellIndexI  ][JCu].y(),
					Cell[CellIndexI][JCu].Xc.x, Cell[CellIndexI][JCu].Xc.y, OrderX, OrderY) + // cell West side
		PolynomLineIntegration2(Node[CellIndexI  ][JCu].x(),Node[CellIndexI  ][JCu].y(),
					Node[CellIndexI+1][JCu].x(),Node[CellIndexI+1][JCu].y(),
					Cell[CellIndexI][JCu].Xc.x, Cell[CellIndexI][JCu].Xc.y, OrderX, OrderY) + // cell South side
		PolynomLineIntegration2(Node[CellIndexI+1][JCu].x(),Node[CellIndexI+1][JCu].y(),
					Node[CellIndexI+1][JNu].x(),Node[CellIndexI+1][JNu].y(),
					Cell[CellIndexI][JCu].Xc.x, Cell[CellIndexI][JCu].Xc.y, OrderX, OrderY) ) // cell East side
	       /(OrderX+1)/Cell[CellIndexI][JCu].A);  // Division by (OrderX+1) & Division by A
    
    
    case CORNER_NORTH_WEST_SPLINES:     // Use the North Spline and the West Spline
      return ( (BndNorthSplineInfo[CellIndexI].IntegratePolynomialTerm(Cell[ICl][JCu].Xc,
								       OrderX, OrderY) +                       // cell North side
		BndWestSplineInfo[JCu].IntegratePolynomialTerm(Cell[ICl][JCu].Xc, 
							       OrderX, OrderY) +                               // cell West side
		PolynomLineIntegration2(Node[CellIndexI  ][JCu].x(),Node[CellIndexI  ][JCu].y(),
					Node[CellIndexI+1][JCu].x(),Node[CellIndexI+1][JCu].y(),
					Cell[ICl][JCu].Xc.x, Cell[ICl][JCu].Xc.y, OrderX, OrderY) +              // cell South side
		PolynomLineIntegration2(Node[CellIndexI+1][JCu].x(),Node[CellIndexI+1][JCu].y(),
					Node[CellIndexI+1][JNu].x(),Node[CellIndexI+1][JNu].y(),
					Cell[ICl][JCu].Xc.x, Cell[ICl][JCu].Xc.y, OrderX, OrderY) )              // cell East side
	       /(OrderX+1)/Cell[ICl][JCu].A);        // Division by (OrderX+1) & Division by A
      
    
    case CORNER_NORTH_EAST_SPLINES:    // Use the North Spline and the East Spline
      return ( (BndNorthSplineInfo[CellIndexI].IntegratePolynomialTerm(Cell[ICu][JCu].Xc, 
								       OrderX, OrderY) +                        // cell North side
		PolynomLineIntegration2(Node[CellIndexI  ][JNu].x(),Node[CellIndexI  ][JNu].y(),
					Node[CellIndexI  ][JCu].x(),Node[CellIndexI  ][JCu].y(),
					Cell[ICu][JCu].Xc.x, Cell[ICu][JCu].Xc.y, OrderX, OrderY) +              // cell West side
		PolynomLineIntegration2(Node[CellIndexI  ][JCu].x(),Node[CellIndexI  ][JCu].y(),
					Node[CellIndexI+1][JCu].x(),Node[CellIndexI+1][JCu].y(),
					Cell[ICu][JCu].Xc.x, Cell[ICu][JCu].Xc.y, OrderX, OrderY) +              // cell South side
		BndEastSplineInfo[JCu].IntegratePolynomialTerm(Cell[ICu][JCu].Xc,
							       OrderX, OrderY) )                                // cell East side
	       /(OrderX+1)/Cell[ICu][JCu].A);        // Division by (OrderX+1) & Division by A 
      
    case SOUTH_SPLINE:        // Use only the South Spline
      return ( (PolynomLineIntegration2(Node[CellIndexI+1][JCl+1].x(),Node[CellIndexI+1][JCl+1].y(),
					Node[CellIndexI  ][JCl+1].x(),Node[CellIndexI  ][JCl+1].y(),
					Cell[CellIndexI][JCl].Xc.x,Cell[CellIndexI][JCl].Xc.y, OrderX, OrderY) +   // cell North side
		PolynomLineIntegration2(Node[CellIndexI  ][JCl+1].x(),Node[CellIndexI  ][JCl+1].y(),
					Node[CellIndexI  ][JCl  ].x(),Node[CellIndexI  ][JCl  ].y(),
					Cell[CellIndexI][JCl].Xc.x,Cell[CellIndexI][JCl].Xc.y, OrderX, OrderY) +   // cell West side
		BndSouthSplineInfo[CellIndexI].IntegratePolynomialTerm(Cell[CellIndexI][JCl].Xc,
								       OrderX, OrderY) +                           // cell South side
		PolynomLineIntegration2(Node[CellIndexI+1][JCl  ].x(),Node[CellIndexI+1][JCl  ].y(),
					Node[CellIndexI+1][JCl+1].x(),Node[CellIndexI+1][JCl+1].y(),
					Cell[CellIndexI][JCl].Xc.x,Cell[CellIndexI][JCl].Xc.y, OrderX, OrderY) )   // cell East side
	       /(OrderX+1)/Cell[CellIndexI][JCl].A);        // Division by (OrderX+1) & Division by A 

    case CORNER_SOUTH_WEST_SPLINES:   // Use the South Spline and the West Spline
      return ( (PolynomLineIntegration2(Node[CellIndexI+1][JCl+1].x(),Node[CellIndexI+1][JCl+1].y(),
					Node[CellIndexI  ][JCl+1].x(),Node[CellIndexI  ][JCl+1].y(),
					Cell[ICl][JCl].Xc.x,Cell[ICl][JCl].Xc.y, OrderX, OrderY) +               // cell North side
		BndWestSplineInfo[JCl].IntegratePolynomialTerm(Cell[ICl][JCl].Xc, 
							       OrderX, OrderY) +                                // cell West side
		BndSouthSplineInfo[CellIndexI].IntegratePolynomialTerm(Cell[ICl][JCl].Xc, 
								       OrderX, OrderY) +                         // cell South side
		PolynomLineIntegration2(Node[CellIndexI+1][JCl  ].x(),Node[CellIndexI+1][JCl  ].y(),
					Node[CellIndexI+1][JCl+1].x(),Node[CellIndexI+1][JCl+1].y(),
					Cell[ICl][JCl].Xc.x,Cell[ICl][JCl].Xc.y, OrderX, OrderY) )               // cell East side
	       /(OrderX+1)/Cell[ICl][JCl].A);        // Division by (OrderX+1) & Division by A
      
    case CORNER_SOUTH_EAST_SPLINES:   // Use the South Spline and the East Spline
      return ( (PolynomLineIntegration2(Node[CellIndexI+1][JCl+1].x(),Node[CellIndexI+1][JCl+1].y(),
					Node[CellIndexI  ][JCl+1].x(),Node[CellIndexI  ][JCl+1].y(),
					Cell[ICu][JCl].Xc.x,Cell[ICu][JCl].Xc.y, OrderX, OrderY) +               // cell North side
		PolynomLineIntegration2(Node[CellIndexI  ][JCl+1].x(),Node[CellIndexI  ][JCl+1].y(),
					Node[CellIndexI  ][JCl  ].x(),Node[CellIndexI  ][JCl  ].y(),
					Cell[ICu][JCl].Xc.x,Cell[ICu][JCl].Xc.y, OrderX, OrderY) +               // cell West side
		BndSouthSplineInfo[CellIndexI].IntegratePolynomialTerm(Cell[ICu][JCl].Xc, 
								       OrderX, OrderY) +                        // cell South side
		BndEastSplineInfo[JCl].IntegratePolynomialTerm(Cell[ICu][JCl].Xc, 
							       OrderX, OrderY) )                                // cell East side
	       /(OrderX+1)/Cell[ICu][JCl].A);        // Division by (OrderX+1) & Division by A

    case WEST_SPLINE:       // Use only the West Spline
      return ( (PolynomLineIntegration2(Node[ICl+1][CellIndexJ+1].x(),Node[ICl+1][CellIndexJ+1].y(),
					Node[ICl  ][CellIndexJ+1].x(),Node[ICl  ][CellIndexJ+1].y(),
					Cell[ICl][CellIndexJ].Xc.x,Cell[ICl][CellIndexJ].Xc.y, OrderX, OrderY) + // cell North side
		BndWestSplineInfo[CellIndexJ].IntegratePolynomialTerm(Cell[ICl][CellIndexJ].Xc,
								      OrderX, OrderY) +                          // cell West side
		PolynomLineIntegration2(Node[ICl  ][CellIndexJ  ].x(),Node[ICl  ][CellIndexJ  ].y(),
					Node[ICl+1][CellIndexJ  ].x(),Node[ICl+1][CellIndexJ  ].y(),
					Cell[ICl][CellIndexJ].Xc.x,Cell[ICl][CellIndexJ].Xc.y, OrderX, OrderY) + // cell South side
		PolynomLineIntegration2(Node[ICl+1][CellIndexJ  ].x(),Node[ICl+1][CellIndexJ  ].y(),
					Node[ICl+1][CellIndexJ+1].x(),Node[ICl+1][CellIndexJ+1].y(),
					Cell[ICl][CellIndexJ].Xc.x,Cell[ICl][CellIndexJ].Xc.y, OrderX, OrderY) ) // cell East side
	       /(OrderX+1)/Cell[ICl][CellIndexJ].A);        // Division by (OrderX+1) & Division by A

    case EAST_SPLINE:      // Use only the East Spline
      return ( (PolynomLineIntegration2(Node[ICu+1][CellIndexJ+1].x(),Node[ICu+1][CellIndexJ+1].y(),
					Node[ICu  ][CellIndexJ+1].x(),Node[ICu  ][CellIndexJ+1].y(),
					Cell[ICu][CellIndexJ].Xc.x,Cell[ICu][CellIndexJ].Xc.y, OrderX, OrderY) + // cell North side
		PolynomLineIntegration2(Node[ICu  ][CellIndexJ+1].x(),Node[ICu  ][CellIndexJ+1].y(),
					Node[ICu  ][CellIndexJ  ].x(),Node[ICu  ][CellIndexJ  ].y(),
					Cell[ICu][CellIndexJ].Xc.x,Cell[ICu][CellIndexJ].Xc.y, OrderX, OrderY) + // cell West side
		PolynomLineIntegration2(Node[ICu  ][CellIndexJ  ].x(),Node[ICu  ][CellIndexJ  ].y(),
					Node[ICu+1][CellIndexJ  ].x(),Node[ICu+1][CellIndexJ  ].y(),
					Cell[ICu][CellIndexJ].Xc.x,Cell[ICu][CellIndexJ].Xc.y, OrderX, OrderY) + // cell South side
		BndEastSplineInfo[CellIndexJ].IntegratePolynomialTerm(Cell[ICu][CellIndexJ].Xc,
								      OrderX, OrderY) )                          // cell East side
	       /(OrderX+1)/Cell[ICu][CellIndexJ].A);        // Division by (OrderX+1) & Division by A

    default:
      return 0.0;
    } // endswitch

  } else {

    // Use the boundary spline functions to integrate along curved edges.

    switch(Boundary){

    case NORTH_SPLINE:          // Use only the North Spline -> (iCell,jCell)=(CellIndexI,JCu)
      return ( (BndNorthSpline.PolynomOrderIntegration(Node[CellIndexI+1][JNu], Node[CellIndexI  ][JNu], 
						       Cell[CellIndexI][JCu].Xc, 15, OrderX, OrderY) +          // cell North side
		PolynomLineIntegration2(Node[CellIndexI  ][JNu].x(),Node[CellIndexI  ][JNu].y(),
					Node[CellIndexI  ][JCu].x(),Node[CellIndexI  ][JCu].y(),
					Cell[CellIndexI][JCu].Xc.x, Cell[CellIndexI][JCu].Xc.y, OrderX, OrderY) + // cell West side
		PolynomLineIntegration2(Node[CellIndexI  ][JCu].x(),Node[CellIndexI  ][JCu].y(),
					Node[CellIndexI+1][JCu].x(),Node[CellIndexI+1][JCu].y(),
					Cell[CellIndexI][JCu].Xc.x, Cell[CellIndexI][JCu].Xc.y, OrderX, OrderY) + // cell South side
		PolynomLineIntegration2(Node[CellIndexI+1][JCu].x(),Node[CellIndexI+1][JCu].y(),
					Node[CellIndexI+1][JNu].x(),Node[CellIndexI+1][JNu].y(),
					Cell[CellIndexI][JCu].Xc.x, Cell[CellIndexI][JCu].Xc.y, OrderX, OrderY) ) // cell East side
	       /(OrderX+1)/Cell[CellIndexI][JCu].A);  // Division by (OrderX+1) & Division by A
    
    
    case CORNER_NORTH_WEST_SPLINES:     // Use the North Spline and the West Spline
      return ( (BndNorthSpline.PolynomOrderIntegration(Node[CellIndexI+1][JNu], Node[CellIndexI  ][JNu],
						       Cell[ICl][JCu].Xc, 15, OrderX, OrderY) +                 // cell North side
		BndWestSpline.PolynomOrderIntegration(Node[CellIndexI ][JNu], Node[CellIndexI  ][JCu],
						      Cell[ICl][JCu].Xc, 15, OrderX, OrderY) +                  // cell West side
		PolynomLineIntegration2(Node[CellIndexI  ][JCu].x(),Node[CellIndexI  ][JCu].y(),
					Node[CellIndexI+1][JCu].x(),Node[CellIndexI+1][JCu].y(),
					Cell[ICl][JCu].Xc.x, Cell[ICl][JCu].Xc.y, OrderX, OrderY) +              // cell South side
		PolynomLineIntegration2(Node[CellIndexI+1][JCu].x(),Node[CellIndexI+1][JCu].y(),
					Node[CellIndexI+1][JNu].x(),Node[CellIndexI+1][JNu].y(),
					Cell[ICl][JCu].Xc.x, Cell[ICl][JCu].Xc.y, OrderX, OrderY) )              // cell East side
	       /(OrderX+1)/Cell[ICl][JCu].A);        // Division by (OrderX+1) & Division by A
    
    
    case CORNER_NORTH_EAST_SPLINES:    // Use the North Spline and the East Spline
      return ( (BndNorthSpline.PolynomOrderIntegration(Node[CellIndexI+1][JNu], Node[CellIndexI  ][JNu],
						       Cell[ICu][JCu].Xc, 15, OrderX, OrderY) +            // cell North side
		PolynomLineIntegration2(Node[CellIndexI  ][JNu].x(),Node[CellIndexI  ][JNu].y(),
					Node[CellIndexI  ][JCu].x(),Node[CellIndexI  ][JCu].y(),
					Cell[ICu][JCu].Xc.x, Cell[ICu][JCu].Xc.y, OrderX, OrderY) +              // cell West side
		PolynomLineIntegration2(Node[CellIndexI  ][JCu].x(),Node[CellIndexI  ][JCu].y(),
					Node[CellIndexI+1][JCu].x(),Node[CellIndexI+1][JCu].y(),
					Cell[ICu][JCu].Xc.x, Cell[ICu][JCu].Xc.y, OrderX, OrderY) +              // cell South side
		BndEastSpline.PolynomOrderIntegration(Node[CellIndexI+1][JCu], Node[CellIndexI+1][JNu],
						      Cell[ICu][JCu].Xc, 15, OrderX, OrderY) )                  // cell East side
	       /(OrderX+1)/Cell[ICu][JCu].A);        // Division by (OrderX+1) & Division by A 

    case SOUTH_SPLINE:        // Use only the South Spline
      return ( (PolynomLineIntegration2(Node[CellIndexI+1][JCl+1].x(),Node[CellIndexI+1][JCl+1].y(),
					Node[CellIndexI  ][JCl+1].x(),Node[CellIndexI  ][JCl+1].y(),
					Cell[CellIndexI][JCl].Xc.x,Cell[CellIndexI][JCl].Xc.y, OrderX, OrderY) +   // cell North side
		PolynomLineIntegration2(Node[CellIndexI  ][JCl+1].x(),Node[CellIndexI  ][JCl+1].y(),
					Node[CellIndexI  ][JCl  ].x(),Node[CellIndexI  ][JCl  ].y(),
					Cell[CellIndexI][JCl].Xc.x,Cell[CellIndexI][JCl].Xc.y, OrderX, OrderY) +   // cell West side
		BndSouthSpline.PolynomOrderIntegration(Node[CellIndexI ][JCl], Node[CellIndexI+1][JCl],
						       Cell[CellIndexI][JCl].Xc, 15, OrderX, OrderY) +            // cell South side
		PolynomLineIntegration2(Node[CellIndexI+1][JCl  ].x(),Node[CellIndexI+1][JCl  ].y(),
					Node[CellIndexI+1][JCl+1].x(),Node[CellIndexI+1][JCl+1].y(),
					Cell[CellIndexI][JCl].Xc.x,Cell[CellIndexI][JCl].Xc.y, OrderX, OrderY) )   // cell East side
	       /(OrderX+1)/Cell[CellIndexI][JCl].A);        // Division by (OrderX+1) & Division by A 

    case CORNER_SOUTH_WEST_SPLINES:   // Use the South Spline and the West Spline
      return ( (PolynomLineIntegration2(Node[CellIndexI+1][JCl+1].x(),Node[CellIndexI+1][JCl+1].y(),
					Node[CellIndexI  ][JCl+1].x(),Node[CellIndexI  ][JCl+1].y(),
					Cell[ICl][JCl].Xc.x,Cell[ICl][JCl].Xc.y, OrderX, OrderY) +               // cell North side
		BndWestSpline.PolynomOrderIntegration(Node[CellIndexI ][JCl+1],Node[CellIndexI ][JCl ],
						      Cell[ICl][JCl].Xc, 15, OrderX, OrderY) +                  // cell West side
		BndSouthSpline.PolynomOrderIntegration(Node[CellIndexI ][JCl], Node[CellIndexI+1][JCl],
						       Cell[ICl][JCl].Xc, 15, OrderX, OrderY) +                 // cell South side
		PolynomLineIntegration2(Node[CellIndexI+1][JCl  ].x(),Node[CellIndexI+1][JCl  ].y(),
					Node[CellIndexI+1][JCl+1].x(),Node[CellIndexI+1][JCl+1].y(),
					Cell[ICl][JCl].Xc.x,Cell[ICl][JCl].Xc.y, OrderX, OrderY) )               // cell East side
	       /(OrderX+1)/Cell[ICl][JCl].A);        // Division by (OrderX+1) & Division by A

    case CORNER_SOUTH_EAST_SPLINES:   // Use the South Spline and the East Spline
      return ( (PolynomLineIntegration2(Node[CellIndexI+1][JCl+1].x(),Node[CellIndexI+1][JCl+1].y(),
					Node[CellIndexI  ][JCl+1].x(),Node[CellIndexI  ][JCl+1].y(),
					Cell[ICu][JCl].Xc.x,Cell[ICu][JCl].Xc.y, OrderX, OrderY) +               // cell North side
		PolynomLineIntegration2(Node[CellIndexI  ][JCl+1].x(),Node[CellIndexI  ][JCl+1].y(),
					Node[CellIndexI  ][JCl  ].x(),Node[CellIndexI  ][JCl  ].y(),
					Cell[ICu][JCl].Xc.x,Cell[ICu][JCl].Xc.y, OrderX, OrderY) +               // cell West side
		BndSouthSpline.PolynomOrderIntegration(Node[CellIndexI ][JCl], Node[CellIndexI+1][JCl],
						       Cell[ICu][JCl].Xc, 15, OrderX, OrderY) +                 // cell South side
		BndEastSpline.PolynomOrderIntegration(Node[CellIndexI+1][JCl], Node[CellIndexI+1][JCl+1],
						      Cell[ICu][JCl].Xc, 15, OrderX, OrderY) )                  // cell East side
	       /(OrderX+1)/Cell[ICu][JCl].A);        // Division by (OrderX+1) & Division by A

    case WEST_SPLINE:       // Use only the West Spline
      return ( (PolynomLineIntegration2(Node[ICl+1][CellIndexJ+1].x(),Node[ICl+1][CellIndexJ+1].y(),
					Node[ICl  ][CellIndexJ+1].x(),Node[ICl  ][CellIndexJ+1].y(),
					Cell[ICl][CellIndexJ].Xc.x,Cell[ICl][CellIndexJ].Xc.y, OrderX, OrderY) + // cell North side
		BndWestSpline.PolynomOrderIntegration(Node[ICl][CellIndexJ+1], Node[ICl][CellIndexJ],
						      Cell[ICl][CellIndexJ].Xc, 15, OrderX, OrderY) +           // cell West side
		PolynomLineIntegration2(Node[ICl  ][CellIndexJ  ].x(),Node[ICl  ][CellIndexJ  ].y(),
					Node[ICl+1][CellIndexJ  ].x(),Node[ICl+1][CellIndexJ  ].y(),
					Cell[ICl][CellIndexJ].Xc.x,Cell[ICl][CellIndexJ].Xc.y, OrderX, OrderY) + // cell South side
		PolynomLineIntegration2(Node[ICl+1][CellIndexJ  ].x(),Node[ICl+1][CellIndexJ  ].y(),
					Node[ICl+1][CellIndexJ+1].x(),Node[ICl+1][CellIndexJ+1].y(),
					Cell[ICl][CellIndexJ].Xc.x,Cell[ICl][CellIndexJ].Xc.y, OrderX, OrderY) ) // cell East side
	       /(OrderX+1)/Cell[ICl][CellIndexJ].A);        // Division by (OrderX+1) & Division by A

    case EAST_SPLINE:      // Use only the East Spline
      return ( (PolynomLineIntegration2(Node[ICu+1][CellIndexJ+1].x(),Node[ICu+1][CellIndexJ+1].y(),
					Node[ICu  ][CellIndexJ+1].x(),Node[ICu  ][CellIndexJ+1].y(),
					Cell[ICu][CellIndexJ].Xc.x,Cell[ICu][CellIndexJ].Xc.y, OrderX, OrderY) + // cell North side
		PolynomLineIntegration2(Node[ICu  ][CellIndexJ+1].x(),Node[ICu  ][CellIndexJ+1].y(),
					Node[ICu  ][CellIndexJ  ].x(),Node[ICu  ][CellIndexJ  ].y(),
					Cell[ICu][CellIndexJ].Xc.x,Cell[ICu][CellIndexJ].Xc.y, OrderX, OrderY) + // cell West side
		PolynomLineIntegration2(Node[ICu  ][CellIndexJ  ].x(),Node[ICu  ][CellIndexJ  ].y(),
					Node[ICu+1][CellIndexJ  ].x(),Node[ICu+1][CellIndexJ  ].y(),
					Cell[ICu][CellIndexJ].Xc.x,Cell[ICu][CellIndexJ].Xc.y, OrderX, OrderY) + // cell South side
		BndEastSpline.PolynomOrderIntegration(Node[ICu+1][CellIndexJ], Node[ICu+1][CellIndexJ+1],
						      Cell[ICu][CellIndexJ].Xc, 15, OrderX, OrderY) )           // cell East side
	       /(OrderX+1)/Cell[ICu][CellIndexJ].A);        // Division by (OrderX+1) & Division by A

    default:
      return 0.0;
    } // endswitch

  }// endif

}

/*!
 * Compute a particular geometric moment for a ghost cell.
 *
 * \param CellIndexI i-index of the cell
 * \param CellIndexJ j-index of the cell
 * \param Boundary indicates which cell boundaries are curved
 * \param OrderX the x-power of the polynomial function which is integrated
 * \param OrderY the y-power of the polynomial function which is integrated
 */
double Grid2D_Quad_Block_HO::GeometricMoment_GhostCell_CurvedBoundaries(const int &CellIndexI, const int &CellIndexJ,
									const int &Boundary,
									const int &OrderX, const int &OrderY){

  // Obs. The sides of the cell that are not curved are treated as line segments and therefore
  //      the line integral is computed exactly based on the Nodes.
  // The edges are considered in counterclockwise order (i.e. right-hand rule applied)
  // ***! This subroutine is for the first row of ghost cells and therefore only one boundary can be curved at a time

  if (Gauss_Quad_Curvilinear_Integration) {

    // Use the SplineInfo variables to integrate along curved edges.

    switch(Boundary){

    case NORTH_SPLINE:          // Use only the North Spline -> (iCell,jCell)=(CellIndexI,JCu)
      return ( (// cell North side
		PolynomLineIntegration2(Node[CellIndexI+1][JNu+1],
					Node[CellIndexI  ][JNu+1],
					Cell[CellIndexI][JCu+1].Xc.x, Cell[CellIndexI][JCu+1].Xc.y, OrderX, OrderY) +
		// cell West side
		PolynomLineIntegration2(Node[CellIndexI  ][JNu+1],
					Node[CellIndexI  ][JNu  ],
					Cell[CellIndexI][JCu+1].Xc.x, Cell[CellIndexI][JCu+1].Xc.y, OrderX, OrderY) - 
		// cell South side
		BndNorthSplineInfo[CellIndexI].IntegratePolynomialTerm(Cell[CellIndexI][JCu+1].Xc,
								       OrderX, OrderY) + 
		// cell East side
		PolynomLineIntegration2(Node[CellIndexI+1][JNu  ],
					Node[CellIndexI+1][JNu+1],
					Cell[CellIndexI][JCu+1].Xc.x, Cell[CellIndexI][JCu+1].Xc.y, OrderX, OrderY) )
	       /(OrderX+1)/Cell[CellIndexI][JCu+1].A);  // Division by (OrderX+1) & Division by A
    
    case SOUTH_SPLINE:        // Use only the South Spline
      return ( (// cell North side
		-BndSouthSplineInfo[CellIndexI].IntegratePolynomialTerm(Cell[CellIndexI][JCl-1].Xc,
									OrderX, OrderY) +   
		// cell West side
		PolynomLineIntegration2(Node[CellIndexI  ][JNl  ],
					Node[CellIndexI  ][JNl-1],
					Cell[CellIndexI][JCl-1].Xc.x,Cell[CellIndexI][JCl-1].Xc.y, OrderX, OrderY) +   
		// cell South side
		PolynomLineIntegration2(Node[CellIndexI  ][JNl-1],
					Node[CellIndexI+1][JNl-1],
					Cell[CellIndexI][JCl-1].Xc.x,Cell[CellIndexI][JCl-1].Xc.y, OrderX, OrderY) +   
		// cell East side
		PolynomLineIntegration2(Node[CellIndexI+1][JNl-1],
					Node[CellIndexI+1][JNl  ],
					Cell[CellIndexI][JCl-1].Xc.x,Cell[CellIndexI][JCl-1].Xc.y, OrderX, OrderY) )   
	       /(OrderX+1)/Cell[CellIndexI][JCl-1].A);        // Division by (OrderX+1) & Division by A 

    case WEST_SPLINE:       // Use only the West Spline
      return ( (// cell North side
		PolynomLineIntegration2(Node[INl  ][CellIndexJ+1],
					Node[INl-1][CellIndexJ+1],
					Cell[ICl-1][CellIndexJ].Xc.x,Cell[ICl-1][CellIndexJ].Xc.y, OrderX, OrderY) + 
		// cell West side
		PolynomLineIntegration2(Node[INl-1][CellIndexJ+1],
					Node[INl-1][CellIndexJ  ],
					Cell[ICl-1][CellIndexJ].Xc.x,Cell[ICl-1][CellIndexJ].Xc.y, OrderX, OrderY) + 
		// cell South side
		PolynomLineIntegration2(Node[INl-1][CellIndexJ  ],
					Node[INl  ][CellIndexJ  ],
					Cell[ICl-1][CellIndexJ].Xc.x,Cell[ICl-1][CellIndexJ].Xc.y, OrderX, OrderY) - 
		// cell East side
		BndWestSplineInfo[CellIndexJ].IntegratePolynomialTerm(Cell[ICl-1][CellIndexJ].Xc, 
								      OrderX, OrderY) ) 
	       /(OrderX+1)/Cell[ICl-1][CellIndexJ].A);        // Division by (OrderX+1) & Division by A

    case EAST_SPLINE:      // Use only the East Spline
      return ( (// cell North side
		PolynomLineIntegration2(Node[INu+1][CellIndexJ+1],
					Node[INu  ][CellIndexJ+1],
					Cell[ICu+1][CellIndexJ].Xc.x,Cell[ICu+1][CellIndexJ].Xc.y, OrderX, OrderY) - 
		// cell West side
		BndEastSplineInfo[CellIndexJ].IntegratePolynomialTerm(Cell[ICu+1][CellIndexJ].Xc, 
								      OrderX, OrderY) + 
		// cell South side
		PolynomLineIntegration2(Node[INu  ][CellIndexJ  ],
					Node[INu+1][CellIndexJ  ],
					Cell[ICu+1][CellIndexJ].Xc.x,Cell[ICu+1][CellIndexJ].Xc.y, OrderX, OrderY) + 
		// cell East side
		PolynomLineIntegration2(Node[INu+1][CellIndexJ  ],
					Node[INu+1][CellIndexJ+1],
					Cell[ICu+1][CellIndexJ].Xc.x,Cell[ICu+1][CellIndexJ].Xc.y, OrderX, OrderY) ) 
	       /(OrderX+1)/Cell[ICu+1][CellIndexJ].A);        // Division by (OrderX+1) & Division by A

    case EXTEND_W_NORTH_RIGHT_SPLINE: // Use only the extension of North spline to West (right side)
      return ( (// cell North side
		ExtendWest_BndNorthSplineInfo[CellIndexI].IntegratePolynomialTerm(Cell[CellIndexI][JCu].Xc,
										  OrderX, OrderY) +
		// cell West side
		PolynomLineIntegration2(Node[CellIndexI  ][JNu  ],
					Node[CellIndexI  ][JNu-1],
					Cell[CellIndexI][JCu].Xc.x, Cell[CellIndexI][JCu].Xc.y, OrderX, OrderY) + 
		// cell South side
		PolynomLineIntegration2(Node[CellIndexI  ][JNu-1],
					Node[CellIndexI+1][JNu-1],
					Cell[CellIndexI][JCu].Xc.x, Cell[CellIndexI][JCu].Xc.y, OrderX, OrderY) +
		// cell East side
		PolynomLineIntegration2(Node[CellIndexI+1][JNu-1],
					Node[CellIndexI+1][JNu  ],
					Cell[CellIndexI][JCu].Xc.x, Cell[CellIndexI][JCu].Xc.y, OrderX, OrderY) )
	       /(OrderX+1)/Cell[CellIndexI][JCu].A);        // Division by (OrderX+1) & Division by A

    case EXTEND_E_NORTH_RIGHT_SPLINE: // Use only the extension of North spline to East (right side)
      return ( (// cell North side
		PolynomLineIntegration2(Node[CellIndexI+1][JNu+1],
					Node[CellIndexI  ][JNu+1],
					Cell[CellIndexI][JCu+1].Xc.x, Cell[CellIndexI][JCu+1].Xc.y, OrderX, OrderY) + 
		// cell West side
		PolynomLineIntegration2(Node[CellIndexI  ][JNu+1],
					Node[CellIndexI  ][JNu  ],
					Cell[CellIndexI][JCu+1].Xc.x, Cell[CellIndexI][JCu+1].Xc.y, OrderX, OrderY) -   
		// cell South side
		ExtendEast_BndNorthSplineInfo[CellIndexI-(ICu+1)].IntegratePolynomialTerm(Cell[CellIndexI][JCu+1].Xc,
											  OrderX, OrderY) +      
		// cell East side
		PolynomLineIntegration2(Node[CellIndexI+1][JNu  ],
					Node[CellIndexI+1][JNu+1],
					Cell[CellIndexI][JCu+1].Xc.x, Cell[CellIndexI][JCu+1].Xc.y, OrderX, OrderY) )
	       /(OrderX+1)/Cell[CellIndexI][JCu+1].A);  // Division by (OrderX+1) & Division by A

      
    case EXTEND_W_NORTH_LEFT_SPLINE: // Use only the extension of North spline to West (left side)
      return ( (// cell North side
		PolynomLineIntegration2(Node[CellIndexI+1][JNu+1],
					Node[CellIndexI  ][JNu+1],
					Cell[CellIndexI][JCu+1].Xc.x,Cell[CellIndexI][JCu+1].Xc.y, OrderX, OrderY) + 
		// cell West side
		PolynomLineIntegration2(Node[CellIndexI  ][JNu+1],
					Node[CellIndexI  ][JNu  ],
					Cell[CellIndexI][JCu+1].Xc.x,Cell[CellIndexI][JCu+1].Xc.y, OrderX, OrderY) -   
		// cell South side
		ExtendWest_BndNorthSplineInfo[CellIndexI].IntegratePolynomialTerm(Cell[CellIndexI][JCu+1].Xc,
										  OrderX, OrderY) +      
		// cell East side
		PolynomLineIntegration2(Node[CellIndexI+1][JNu  ],
					Node[CellIndexI+1][JNu+1],
					Cell[CellIndexI][JCu+1].Xc.x,Cell[CellIndexI][JCu+1].Xc.y, OrderX, OrderY) )
	       /(OrderX+1)/Cell[CellIndexI][JCu+1].A);        // Division by (OrderX+1) & Division by A


    case EXTEND_E_NORTH_LEFT_SPLINE: // Use only the extension of North spline to East (left side)
      return ( (// cell North side
		ExtendEast_BndNorthSplineInfo[CellIndexI-(ICu+1)].IntegratePolynomialTerm(Cell[CellIndexI][JCu].Xc,
											  OrderX, OrderY) +
		// cell West side
		PolynomLineIntegration2(Node[CellIndexI  ][JNu  ],
					Node[CellIndexI  ][JNu-1],
					Cell[CellIndexI][JCu].Xc.x, Cell[CellIndexI][JCu].Xc.y, OrderX, OrderY) +   
		// cell South side
		PolynomLineIntegration2(Node[CellIndexI  ][JNu-1],
					Node[CellIndexI+1][JNu-1],
					Cell[CellIndexI][JCu].Xc.x, Cell[CellIndexI][JCu].Xc.y, OrderX, OrderY) + 
		// cell East side
		PolynomLineIntegration2(Node[CellIndexI+1][JNu-1],
					Node[CellIndexI+1][JNu  ],
					Cell[CellIndexI][JCu].Xc.x, Cell[CellIndexI][JCu].Xc.y, OrderX, OrderY) )
	       /(OrderX+1)/Cell[CellIndexI][JCu].A);  // Division by (OrderX+1) & Division by A


    case EXTEND_W_SOUTH_RIGHT_SPLINE: // Use only the extension of South spline to West (right side)
      return ( (// cell North side
		-ExtendWest_BndSouthSplineInfo[CellIndexI].IntegratePolynomialTerm(Cell[CellIndexI][JCl-1].Xc,
										   OrderX, OrderY) +
		// cell West side
		PolynomLineIntegration2(Node[CellIndexI  ][JNl  ],
					Node[CellIndexI  ][JNl-1],
					Cell[CellIndexI][JCl-1].Xc.x, Cell[CellIndexI][JCl-1].Xc.y, OrderX, OrderY) +   
		// cell South side
		PolynomLineIntegration2(Node[CellIndexI  ][JNl-1],
					Node[CellIndexI+1][JNl-1],
					Cell[CellIndexI][JCl-1].Xc.x, Cell[CellIndexI][JCl-1].Xc.y, OrderX, OrderY) + 
		// cell East side
		PolynomLineIntegration2(Node[CellIndexI+1][JNl-1],
					Node[CellIndexI+1][JNl  ],
					Cell[CellIndexI][JCl-1].Xc.x, Cell[CellIndexI][JCl-1].Xc.y, OrderX, OrderY) )
	       /(OrderX+1)/Cell[CellIndexI][JCl-1].A);  // Division by (OrderX+1) & Division by A


    case EXTEND_E_SOUTH_RIGHT_SPLINE: // Use only the extension of South spline to East (right side)
      return ( (// cell North side
		PolynomLineIntegration2(Node[CellIndexI+1][JNl+1],
					Node[CellIndexI  ][JNl+1],
					Cell[CellIndexI][JCl].Xc.x,Cell[CellIndexI][JCl].Xc.y,OrderX,OrderY) + 
		// cell West side
		PolynomLineIntegration2(Node[CellIndexI  ][JNl+1],
					Node[CellIndexI  ][JNl  ],
					Cell[CellIndexI][JCl].Xc.x,Cell[CellIndexI][JCl].Xc.y,OrderX,OrderY) +   
		// cell South side
		ExtendEast_BndSouthSplineInfo[CellIndexI-(ICu+1)].IntegratePolynomialTerm(Cell[CellIndexI][JCl].Xc,
											  OrderX,OrderY) +
		// cell East side
		PolynomLineIntegration2(Node[CellIndexI+1][JNl  ],
					Node[CellIndexI+1][JNl+1],
					Cell[CellIndexI][JCl].Xc.x,Cell[CellIndexI][JCl].Xc.y,OrderX,OrderY) )
	       /(OrderX+1)/Cell[CellIndexI][JCl].A);  // Division by (OrderX+1) & Division by A

      
    case EXTEND_W_SOUTH_LEFT_SPLINE: // Use only the extension of South spline to West (left side)
      return ( (// cell North side
		PolynomLineIntegration2(Node[CellIndexI+1][JNl+1],
					Node[CellIndexI  ][JNl+1],
					Cell[CellIndexI][JCl].Xc.x, Cell[CellIndexI][JCl].Xc.y, OrderX, OrderY) + 
		// cell West side
		PolynomLineIntegration2(Node[CellIndexI  ][JNl+1],
					Node[CellIndexI  ][JNl  ],
					Cell[CellIndexI][JCl].Xc.x, Cell[CellIndexI][JCl].Xc.y, OrderX, OrderY) +
		// cell South side
		ExtendWest_BndSouthSplineInfo[CellIndexI].IntegratePolynomialTerm(Cell[CellIndexI][JCl].Xc,
										  OrderX, OrderY) +      
		// cell East side
		PolynomLineIntegration2(Node[CellIndexI+1][JNl  ],
					Node[CellIndexI+1][JNl+1],
					Cell[CellIndexI][JCl].Xc.x, Cell[CellIndexI][JCl].Xc.y, OrderX, OrderY) )
	       /(OrderX+1)/Cell[CellIndexI][JCl].A);  // Division by (OrderX+1) & Division by A


    case EXTEND_E_SOUTH_LEFT_SPLINE: // Use only the extension of South spline to East (left side)
      return ( (// cell North side
		-ExtendEast_BndSouthSplineInfo[CellIndexI-(ICu+1)].IntegratePolynomialTerm(Cell[CellIndexI][JCl-1].Xc,
											   OrderX,OrderY) +
		// cell West side
		PolynomLineIntegration2(Node[CellIndexI  ][JNl  ],
					Node[CellIndexI  ][JNl-1],
					Cell[CellIndexI][JCl-1].Xc.x, Cell[CellIndexI][JCl-1].Xc.y, OrderX, OrderY) +   
		// cell South side
		PolynomLineIntegration2(Node[CellIndexI  ][JNl-1],
					Node[CellIndexI+1][JNl-1],
					Cell[CellIndexI][JCl-1].Xc.x, Cell[CellIndexI][JCl-1].Xc.y, OrderX, OrderY) + 
		// cell East side
		PolynomLineIntegration2(Node[CellIndexI+1][JNl-1],
					Node[CellIndexI+1][JNl  ],
					Cell[CellIndexI][JCl-1].Xc.x, Cell[CellIndexI][JCl-1].Xc.y, OrderX, OrderY) )
	       /(OrderX+1)/Cell[CellIndexI][JCl-1].A);  // Division by (OrderX+1) & Division by A	       


    case EXTEND_N_EAST_RIGHT_SPLINE: // Use only the extension of East spline to North (right side)
      return ( (// cell North side
		PolynomLineIntegration2(Node[INu  ][CellIndexJ+1],
					Node[INu-1][CellIndexJ+1],
					Cell[ICu][CellIndexJ].Xc.x, Cell[ICu][CellIndexJ].Xc.y, OrderX, OrderY) + 
		// cell West side
		PolynomLineIntegration2(Node[INu-1][CellIndexJ+1],
					Node[INu-1][CellIndexJ  ],
					Cell[ICu][CellIndexJ].Xc.x, Cell[ICu][CellIndexJ].Xc.y, OrderX, OrderY) +
		// cell South side
		PolynomLineIntegration2(Node[INu-1][CellIndexJ  ],
					Node[INu  ][CellIndexJ  ],
					Cell[ICu][CellIndexJ].Xc.x, Cell[ICu][CellIndexJ].Xc.y, OrderX, OrderY) + 
		// cell East side
		ExtendNorth_BndEastSplineInfo[CellIndexJ-(JCu+1)].IntegratePolynomialTerm(Cell[ICu][CellIndexJ].Xc,
											  OrderX, OrderY) ) 
	       /(OrderX+1)/Cell[ICu][CellIndexJ].A);  // Division by (OrderX+1) & Division by A


    case EXTEND_S_EAST_RIGHT_SPLINE: // Use only the extension of East spline to South (right side)
      return ( (// cell North side
		PolynomLineIntegration2(Node[INu+1][CellIndexJ+1],
					Node[INu  ][CellIndexJ+1],
					Cell[ICu+1][CellIndexJ].Xc.x, Cell[ICu+1][CellIndexJ].Xc.y, OrderX, OrderY) -
		// cell West side
		ExtendSouth_BndEastSplineInfo[CellIndexJ].IntegratePolynomialTerm(Cell[ICu+1][CellIndexJ].Xc,
										  OrderX, OrderY) + 
		// cell South side
		PolynomLineIntegration2(Node[INu  ][CellIndexJ  ],
					Node[INu+1][CellIndexJ  ],
					Cell[ICu+1][CellIndexJ].Xc.x, Cell[ICu+1][CellIndexJ].Xc.y, OrderX, OrderY) + 
		// cell East side
		PolynomLineIntegration2(Node[INu+1][CellIndexJ  ],
					Node[INu+1][CellIndexJ+1],
					Cell[ICu+1][CellIndexJ].Xc.x, Cell[ICu+1][CellIndexJ].Xc.y, OrderX, OrderY) )
	       /(OrderX+1)/Cell[ICu+1][CellIndexJ].A);  // Division by (OrderX+1) & Division by A


    case EXTEND_N_EAST_LEFT_SPLINE: // Use only the extension of East spline to North (left side)
      return ( (// cell North side
		PolynomLineIntegration2(Node[INu+1][CellIndexJ+1],
					Node[INu  ][CellIndexJ+1],
					Cell[ICu+1][CellIndexJ].Xc.x, Cell[ICu+1][CellIndexJ].Xc.y, OrderX, OrderY) -
		// cell West side
		ExtendNorth_BndEastSplineInfo[CellIndexJ-(JCu+1)].IntegratePolynomialTerm(Cell[ICu+1][CellIndexJ].Xc,
											  OrderX, OrderY) + 
		// cell South side
		PolynomLineIntegration2(Node[INu  ][CellIndexJ  ],
					Node[INu+1][CellIndexJ  ],
					Cell[ICu+1][CellIndexJ].Xc.x, Cell[ICu+1][CellIndexJ].Xc.y, OrderX, OrderY) + 
		// cell East side
		PolynomLineIntegration2(Node[INu+1][CellIndexJ  ],
					Node[INu+1][CellIndexJ+1],
					Cell[ICu+1][CellIndexJ].Xc.x, Cell[ICu+1][CellIndexJ].Xc.y, OrderX, OrderY) ) 
	       /(OrderX+1)/Cell[ICu+1][CellIndexJ].A);  // Division by (OrderX+1) & Division by A


    case EXTEND_S_EAST_LEFT_SPLINE: // Use only the extension of East spline to South (left side)
      return ( (// cell North side
		PolynomLineIntegration2(Node[INu  ][CellIndexJ+1],
					Node[INu-1][CellIndexJ+1],
					Cell[ICu][CellIndexJ].Xc.x, Cell[ICu][CellIndexJ].Xc.y, OrderX, OrderY) + 
		// cell West side
		PolynomLineIntegration2(Node[INu-1][CellIndexJ+1],
					Node[INu-1][CellIndexJ  ],
					Cell[ICu][CellIndexJ].Xc.x, Cell[ICu][CellIndexJ].Xc.y, OrderX, OrderY) +
		// cell South side
		PolynomLineIntegration2(Node[INu-1][CellIndexJ  ],
					Node[INu  ][CellIndexJ  ],
					Cell[ICu][CellIndexJ].Xc.x, Cell[ICu][CellIndexJ].Xc.y, OrderX, OrderY) + 
		// cell East side
		ExtendSouth_BndEastSplineInfo[CellIndexJ].IntegratePolynomialTerm(Cell[ICu][CellIndexJ].Xc,
										  OrderX, OrderY) )
	       /(OrderX+1)/Cell[ICu][CellIndexJ].A);  // Division by (OrderX+1) & Division by A 


    case EXTEND_N_WEST_RIGHT_SPLINE: // Use only the extension of West spline to North (right side)
      return ( (// cell North side
		PolynomLineIntegration2(Node[INl  ][CellIndexJ+1],
					Node[INl-1][CellIndexJ+1],
					Cell[ICl-1][CellIndexJ].Xc.x, Cell[ICl-1][CellIndexJ].Xc.y, OrderX, OrderY) + 
		// cell West side
		PolynomLineIntegration2(Node[INl-1][CellIndexJ+1],
					Node[INl-1][CellIndexJ  ],
					Cell[ICl-1][CellIndexJ].Xc.x, Cell[ICl-1][CellIndexJ].Xc.y, OrderX, OrderY) +
		// cell South side
		PolynomLineIntegration2(Node[INl-1][CellIndexJ  ],
					Node[INl  ][CellIndexJ  ],
					Cell[ICl-1][CellIndexJ].Xc.x, Cell[ICl-1][CellIndexJ].Xc.y, OrderX, OrderY) - 
		// cell East side
		ExtendNorth_BndWestSplineInfo[CellIndexJ-(JCu+1)].IntegratePolynomialTerm(Cell[ICl-1][CellIndexJ].Xc,
											  OrderX, OrderY) )
	       /(OrderX+1)/Cell[ICl-1][CellIndexJ].A);  // Division by (OrderX+1) & Division by A


    case EXTEND_S_WEST_RIGHT_SPLINE: // Use only the extension of West spline to South (right side)
      return ( (// cell North side
		PolynomLineIntegration2(Node[INl+1][CellIndexJ+1],
					Node[INl  ][CellIndexJ+1],
					Cell[ICl][CellIndexJ].Xc.x, Cell[ICl][CellIndexJ].Xc.y, OrderX, OrderY) +
		// cell West side
		ExtendSouth_BndWestSplineInfo[CellIndexJ].IntegratePolynomialTerm(Cell[ICl][CellIndexJ].Xc,
										  OrderX, OrderY) + 
		// cell South side
		PolynomLineIntegration2(Node[INl  ][CellIndexJ  ],
					Node[INl+1][CellIndexJ  ],
					Cell[ICl][CellIndexJ].Xc.x, Cell[ICl][CellIndexJ].Xc.y, OrderX, OrderY) + 
		// cell East side
		PolynomLineIntegration2(Node[INl+1][CellIndexJ  ],
					Node[INl+1][CellIndexJ+1],
					Cell[ICl][CellIndexJ].Xc.x, Cell[ICl][CellIndexJ].Xc.y, OrderX, OrderY) )
	       /(OrderX+1)/Cell[ICl][CellIndexJ].A);  // Division by (OrderX+1) & Division by A


    case EXTEND_N_WEST_LEFT_SPLINE: // Use only the extension of West spline to North (left side)
      return ( (// cell North side
		PolynomLineIntegration2(Node[INl+1][CellIndexJ+1],
					Node[INl  ][CellIndexJ+1],
					Cell[ICl][CellIndexJ].Xc.x, Cell[ICl][CellIndexJ].Xc.y, OrderX, OrderY) +
		// cell West side
		ExtendNorth_BndWestSplineInfo[CellIndexJ-(JCu+1)].IntegratePolynomialTerm(Cell[ICl][CellIndexJ].Xc,
											  OrderX, OrderY) + 
		// cell South side
		PolynomLineIntegration2(Node[INl  ][CellIndexJ  ],
					Node[INl+1][CellIndexJ  ],
					Cell[ICl][CellIndexJ].Xc.x, Cell[ICl][CellIndexJ].Xc.y, OrderX, OrderY) + 
		// cell East side
		PolynomLineIntegration2(Node[INl+1][CellIndexJ  ],
					Node[INl+1][CellIndexJ+1],
					Cell[ICl][CellIndexJ].Xc.x, Cell[ICl][CellIndexJ].Xc.y, OrderX, OrderY) )
	       /(OrderX+1)/Cell[ICl][CellIndexJ].A);  // Division by (OrderX+1) & Division by A


    case EXTEND_S_WEST_LEFT_SPLINE: // Use only the extension of West spline to South (left side)
      return ( (// cell North side
		PolynomLineIntegration2(Node[INl  ][CellIndexJ+1],
					Node[INl-1][CellIndexJ+1],
					Cell[ICl-1][CellIndexJ].Xc.x, Cell[ICl-1][CellIndexJ].Xc.y, OrderX, OrderY) + 
		// cell West side
		PolynomLineIntegration2(Node[INl-1][CellIndexJ+1],
					Node[INl-1][CellIndexJ  ],
					Cell[ICl-1][CellIndexJ].Xc.x, Cell[ICl-1][CellIndexJ].Xc.y, OrderX, OrderY) +
		// cell South side
		PolynomLineIntegration2(Node[INl-1][CellIndexJ  ],
					Node[INl  ][CellIndexJ  ],
					Cell[ICl-1][CellIndexJ].Xc.x, Cell[ICl-1][CellIndexJ].Xc.y, OrderX, OrderY) - 
		// cell East side
		ExtendSouth_BndWestSplineInfo[CellIndexJ].IntegratePolynomialTerm(Cell[ICl-1][CellIndexJ].Xc,
										  OrderX, OrderY) )
	       /(OrderX+1)/Cell[ICl-1][CellIndexJ].A);  // Division by (OrderX+1) & Division by A	       


    case CORNER_NORTH_EXTEND_N_WEST_SPLINES: // Use the North spline and the extension of West spline to North
      throw runtime_error("Grid2D_Quad_Block_HO::GeometricMoment_GhostCell_CurvedBoundaries() ERROR! Case CORNER_NORTH_EXTEND_N_WEST_SPLINES has not been encountered so far. Code and test this case!");
      return 0.0;


    case CORNER_EXTEND_N_WEST_EXTEND_W_NORTH_SPLINES: // Use the North extension of West and the West extension of North
      throw runtime_error("Grid2D_Quad_Block_HO::GeometricMoment_GhostCell_CurvedBoundaries() ERROR! Case CORNER_EXTEND_N_WEST_EXTEND_W_NORTH_SPLINES has not been encountered so far. Code and test this case!");
      return 0.0;


    case CORNER_EXTEND_W_NORTH_WEST_SPLINES: // Use the West extension of North and the West spline
      throw runtime_error("Grid2D_Quad_Block_HO::GeometricMoment_GhostCell_CurvedBoundaries() ERROR! Case CORNER_EXTEND_W_NORTH_WEST_SPLINES has not been encountered so far. Code and test this case!");
      return 0.0;


    case CORNER_WEST_EXTEND_W_SOUTH_SPLINES: // Use the West spline and the West extension of South
      throw runtime_error("Grid2D_Quad_Block_HO::GeometricMoment_GhostCell_CurvedBoundaries() ERROR! Case CORNER_WEST_EXTEND_W_SOUTH_SPLINES has not been encountered so far. Code and test this case!");
      return 0.0;


    case CORNER_EXTEND_W_SOUTH_EXTEND_S_WEST_SPLINES: // Use the West extension of South and the South extension of West
      throw runtime_error("Grid2D_Quad_Block_HO::GeometricMoment_GhostCell_CurvedBoundaries() ERROR! Case CORNER_EXTEND_W_SOUTH_EXTEND_S_WEST_SPLINES has not been encountered so far. Code and test this case!");
      return 0.0;


    case CORNER_EXTEND_S_WEST_SOUTH_SPLINES: // Use the South extension of West and the South spline
      throw runtime_error("Grid2D_Quad_Block_HO::GeometricMoment_GhostCell_CurvedBoundaries() ERROR! Case CORNER_EXTEND_S_WEST_SOUTH_SPLINES has not been encountered so far. Code and test this case!");
      return 0.0;


    case CORNER_SOUTH_EXTEND_S_EAST_SPLINES: // Use the South spline and the South extension of East
      throw runtime_error("Grid2D_Quad_Block_HO::GeometricMoment_GhostCell_CurvedBoundaries() ERROR! Case CORNER_SOUTH_EXTEND_S_EAST_SPLINES has not been encountered so far. Code and test this case!");
      return 0.0;


    case CORNER_EXTEND_S_EAST_EXTEND_E_SOUTH_SPLINES: // Use the South extension of East and the East extension of South
      throw runtime_error("Grid2D_Quad_Block_HO::GeometricMoment_GhostCell_CurvedBoundaries() ERROR! Case CORNER_EXTEND_S_EAST_EXTEND_E_SOUTH_SPLINES has not been encountered so far. Code and test this case!");
      return 0.0;


    case CORNER_EXTEND_E_SOUTH_EAST_SPLINES: // Use the East extension of South and the East spline
      throw runtime_error("Grid2D_Quad_Block_HO::GeometricMoment_GhostCell_CurvedBoundaries() ERROR! Case CORNER_EXTEND_E_SOUTH_EAST_SPLINES has not been encountered so far. Code and test this case!");
      return 0.0;


    case CORNER_EAST_EXTEND_E_NORTH_SPLINES: // Use the East spline and the East extension of North
      throw runtime_error("Grid2D_Quad_Block_HO::GeometricMoment_GhostCell_CurvedBoundaries() ERROR! Case CORNER_EAST_EXTEND_E_NORTH_SPLINES has not been encountered so far. Code and test this case!");
      return 0.0;


    case CORNER_EXTEND_E_NORTH_EXTEND_N_EAST_SPLINES: // Use the East extention of North and North extension of East
      throw runtime_error("Grid2D_Quad_Block_HO::GeometricMoment_GhostCell_CurvedBoundaries() ERROR! Case CORNER_EXTEND_E_NORTH_EXTEND_N_EAST_SPLINES has not been encountered so far. Code and test this case!");
      return 0.0;


    case CORNER_EXTEND_N_EAST_NORTH_SPLINES: // Use the North extension of East and the North spline
      throw runtime_error("Grid2D_Quad_Block_HO::GeometricMoment_GhostCell_CurvedBoundaries() ERROR! Case CORNER_EXTEND_N_EAST_NORTH_SPLINES has not been encountered so far. Code and test this case!");
      return 0.0;


    default:
      return 0.0;
    } // endswitch

  } else {

    // Use the boundary spline functions to integrate along curved edges.

    switch(Boundary){

    case NORTH_SPLINE:          // Use only the North Spline -> (iCell,jCell)=(CellIndexI,JCu)
      return ( (// cell North side
		PolynomLineIntegration2(Node[CellIndexI+1][JNu+1],
					Node[CellIndexI  ][JNu+1],
					Cell[CellIndexI][JCu+1].Xc.x, Cell[CellIndexI][JCu+1].Xc.y, OrderX, OrderY) +
		// cell West side
		PolynomLineIntegration2(Node[CellIndexI  ][JNu+1],
					Node[CellIndexI  ][JNu  ],
					Cell[CellIndexI][JCu+1].Xc.x, Cell[CellIndexI][JCu+1].Xc.y, OrderX, OrderY) + 
		// cell South side
		BndNorthSpline.PolynomOrderIntegration(Node[CellIndexI][JNu], Node[CellIndexI+1][JNu], 
						       Cell[CellIndexI][JCu+1].Xc, 15, OrderX, OrderY) + 
		// cell East side
		PolynomLineIntegration2(Node[CellIndexI+1][JNu  ],
					Node[CellIndexI+1][JNu+1],
					Cell[CellIndexI][JCu+1].Xc.x, Cell[CellIndexI][JCu+1].Xc.y, OrderX, OrderY) )
	       /(OrderX+1)/Cell[CellIndexI][JCu+1].A);  // Division by (OrderX+1) & Division by A
    
    case SOUTH_SPLINE:        // Use only the South Spline
      return ( (// cell North side
		BndSouthSpline.PolynomOrderIntegration(Node[CellIndexI+1][JNl], Node[CellIndexI][JNl],
						       Cell[CellIndexI][JCl-1].Xc, 15, OrderX, OrderY) +   
		// cell West side
		PolynomLineIntegration2(Node[CellIndexI  ][JNl  ],
					Node[CellIndexI  ][JNl-1],
					Cell[CellIndexI][JCl-1].Xc.x,Cell[CellIndexI][JCl-1].Xc.y, OrderX, OrderY) +   
		// cell South side
		PolynomLineIntegration2(Node[CellIndexI  ][JNl-1],
					Node[CellIndexI+1][JNl-1],
					Cell[CellIndexI][JCl-1].Xc.x,Cell[CellIndexI][JCl-1].Xc.y, OrderX, OrderY) +   
		// cell East side
		PolynomLineIntegration2(Node[CellIndexI+1][JNl-1],
					Node[CellIndexI+1][JNl  ],
					Cell[CellIndexI][JCl-1].Xc.x,Cell[CellIndexI][JCl-1].Xc.y, OrderX, OrderY) )   
	       /(OrderX+1)/Cell[CellIndexI][JCl-1].A);        // Division by (OrderX+1) & Division by A 

    case WEST_SPLINE:       // Use only the West Spline
      return ( (// cell North side
		PolynomLineIntegration2(Node[INl  ][CellIndexJ+1],
					Node[INl-1][CellIndexJ+1],
					Cell[ICl-1][CellIndexJ].Xc.x,Cell[ICl-1][CellIndexJ].Xc.y, OrderX, OrderY) + 
		// cell West side
		PolynomLineIntegration2(Node[INl-1][CellIndexJ+1],
					Node[INl-1][CellIndexJ  ],
					Cell[ICl-1][CellIndexJ].Xc.x,Cell[ICl-1][CellIndexJ].Xc.y, OrderX, OrderY) + 
		// cell South side
		PolynomLineIntegration2(Node[INl-1][CellIndexJ  ],
					Node[INl  ][CellIndexJ  ],
					Cell[ICl-1][CellIndexJ].Xc.x,Cell[ICl-1][CellIndexJ].Xc.y, OrderX, OrderY) + 
		// cell East side
		BndWestSpline.PolynomOrderIntegration(Node[INl][CellIndexJ], Node[INl][CellIndexJ+1],
						      Cell[ICl-1][CellIndexJ].Xc, 15, OrderX, OrderY) ) 
	       /(OrderX+1)/Cell[ICl-1][CellIndexJ].A);        // Division by (OrderX+1) & Division by A

    case EAST_SPLINE:      // Use only the East Spline
      return ( (// cell North side
		PolynomLineIntegration2(Node[INu+1][CellIndexJ+1],
					Node[INu  ][CellIndexJ+1],
					Cell[ICu+1][CellIndexJ].Xc.x,Cell[ICu+1][CellIndexJ].Xc.y, OrderX, OrderY) + 
		// cell West side
		BndEastSpline.PolynomOrderIntegration(Node[INu][CellIndexJ+1], Node[INu][CellIndexJ],
						      Cell[ICu+1][CellIndexJ].Xc, 15, OrderX, OrderY) + 
		// cell South side
		PolynomLineIntegration2(Node[INu  ][CellIndexJ  ],
					Node[INu+1][CellIndexJ  ],
					Cell[ICu+1][CellIndexJ].Xc.x,Cell[ICu+1][CellIndexJ].Xc.y, OrderX, OrderY) + 
		// cell East side
		PolynomLineIntegration2(Node[INu+1][CellIndexJ  ],
					Node[INu+1][CellIndexJ+1],
					Cell[ICu+1][CellIndexJ].Xc.x,Cell[ICu+1][CellIndexJ].Xc.y, OrderX, OrderY) ) 
	       /(OrderX+1)/Cell[ICu+1][CellIndexJ].A);        // Division by (OrderX+1) & Division by A

    case EXTEND_W_NORTH_RIGHT_SPLINE: // Use only the extension of North spline to West (right side)
      return ( (// cell North side
		ExtendWest_BndNorthSpline.PolynomOrderIntegration(Node[CellIndexI+1][JNu  ],Node[CellIndexI  ][JNu  ],
								  Cell[CellIndexI][JCu].Xc, 15, OrderX, OrderY) +
		// cell West side
		PolynomLineIntegration2(Node[CellIndexI  ][JNu  ],
					Node[CellIndexI  ][JNu-1],
					Cell[CellIndexI][JCu].Xc.x, Cell[CellIndexI][JCu].Xc.y, OrderX, OrderY) + 
		// cell South side
		PolynomLineIntegration2(Node[CellIndexI  ][JNu-1],
					Node[CellIndexI+1][JNu-1],
					Cell[CellIndexI][JCu].Xc.x, Cell[CellIndexI][JCu].Xc.y, OrderX, OrderY) +
		// cell East side
		PolynomLineIntegration2(Node[CellIndexI+1][JNu-1],
					Node[CellIndexI+1][JNu  ],
					Cell[CellIndexI][JCu].Xc.x, Cell[CellIndexI][JCu].Xc.y, OrderX, OrderY) )
	       /(OrderX+1)/Cell[CellIndexI][JCu].A);        // Division by (OrderX+1) & Division by A


    case EXTEND_E_NORTH_RIGHT_SPLINE: // Use only the extension of North spline to East (right side)
      return ( (// cell North side
		PolynomLineIntegration2(Node[CellIndexI+1][JNu+1],
					Node[CellIndexI  ][JNu+1],
					Cell[CellIndexI][JCu+1].Xc.x, Cell[CellIndexI][JCu+1].Xc.y, OrderX, OrderY) + 
		// cell West side
		PolynomLineIntegration2(Node[CellIndexI  ][JNu+1],
					Node[CellIndexI  ][JNu  ],
					Cell[CellIndexI][JCu+1].Xc.x, Cell[CellIndexI][JCu+1].Xc.y, OrderX, OrderY) +   
		// cell South side
		ExtendEast_BndNorthSpline.PolynomOrderIntegration(Node[CellIndexI  ][JNu  ],Node[CellIndexI+1][JNu  ],
								  Cell[CellIndexI][JCu+1].Xc, 15, OrderX, OrderY) +      
		// cell East side
		PolynomLineIntegration2(Node[CellIndexI+1][JNu  ],
					Node[CellIndexI+1][JNu+1],
					Cell[CellIndexI][JCu+1].Xc.x, Cell[CellIndexI][JCu+1].Xc.y, OrderX, OrderY) )
	       /(OrderX+1)/Cell[CellIndexI][JCu+1].A);  // Division by (OrderX+1) & Division by A

      
    case EXTEND_W_NORTH_LEFT_SPLINE: // Use only the extension of North spline to West (left side)
      return ( (// cell North side
		PolynomLineIntegration2(Node[CellIndexI+1][JNu+1],
					Node[CellIndexI  ][JNu+1],
					Cell[CellIndexI][JCu+1].Xc.x,Cell[CellIndexI][JCu+1].Xc.y, OrderX, OrderY) + 
		// cell West side
		PolynomLineIntegration2(Node[CellIndexI  ][JNu+1],
					Node[CellIndexI  ][JNu  ],
					Cell[CellIndexI][JCu+1].Xc.x,Cell[CellIndexI][JCu+1].Xc.y, OrderX, OrderY) +   
		// cell South side
		ExtendWest_BndNorthSpline.PolynomOrderIntegration(Node[CellIndexI  ][JNu  ],Node[CellIndexI+1][JNu  ],
								  Cell[CellIndexI][JCu+1].Xc, 15, OrderX, OrderY) +      
		// cell East side
		PolynomLineIntegration2(Node[CellIndexI+1][JNu  ],
					Node[CellIndexI+1][JNu+1],
					Cell[CellIndexI][JCu+1].Xc.x,Cell[CellIndexI][JCu+1].Xc.y, OrderX, OrderY) )
	       /(OrderX+1)/Cell[CellIndexI][JCu+1].A);        // Division by (OrderX+1) & Division by A


    case EXTEND_E_NORTH_LEFT_SPLINE: // Use only the extension of North spline to East (left side)
      return ( (// cell North side
		ExtendEast_BndNorthSpline.PolynomOrderIntegration(Node[CellIndexI+1][JNu  ], Node[CellIndexI  ][JNu  ],
								  Cell[CellIndexI][JCu].Xc, 15, OrderX, OrderY) +
		// cell West side
		PolynomLineIntegration2(Node[CellIndexI  ][JNu  ],
					Node[CellIndexI  ][JNu-1],
					Cell[CellIndexI][JCu].Xc.x, Cell[CellIndexI][JCu].Xc.y, OrderX, OrderY) +   
		// cell South side
		PolynomLineIntegration2(Node[CellIndexI  ][JNu-1],
					Node[CellIndexI+1][JNu-1],
					Cell[CellIndexI][JCu].Xc.x, Cell[CellIndexI][JCu].Xc.y, OrderX, OrderY) + 
		// cell East side
		PolynomLineIntegration2(Node[CellIndexI+1][JNu-1],
					Node[CellIndexI+1][JNu  ],
					Cell[CellIndexI][JCu].Xc.x, Cell[CellIndexI][JCu].Xc.y, OrderX, OrderY) )
	       /(OrderX+1)/Cell[CellIndexI][JCu].A);  // Division by (OrderX+1) & Division by A


    case EXTEND_W_SOUTH_RIGHT_SPLINE: // Use only the extension of South spline to West (right side)
      return ( (// cell North side
		ExtendWest_BndSouthSpline.PolynomOrderIntegration(Node[CellIndexI+1][JNl  ], Node[CellIndexI  ][JNl  ],
								  Cell[CellIndexI][JCl-1].Xc, 15, OrderX, OrderY) +
		// cell West side
		PolynomLineIntegration2(Node[CellIndexI  ][JNl  ],
					Node[CellIndexI  ][JNl-1],
					Cell[CellIndexI][JCl-1].Xc.x, Cell[CellIndexI][JCl-1].Xc.y, OrderX, OrderY) +   
		// cell South side
		PolynomLineIntegration2(Node[CellIndexI  ][JNl-1],
					Node[CellIndexI+1][JNl-1],
					Cell[CellIndexI][JCl-1].Xc.x, Cell[CellIndexI][JCl-1].Xc.y, OrderX, OrderY) + 
		// cell East side
		PolynomLineIntegration2(Node[CellIndexI+1][JNl-1],
					Node[CellIndexI+1][JNl  ],
					Cell[CellIndexI][JCl-1].Xc.x, Cell[CellIndexI][JCl-1].Xc.y, OrderX, OrderY) )
	       /(OrderX+1)/Cell[CellIndexI][JCl-1].A);  // Division by (OrderX+1) & Division by A


    case EXTEND_E_SOUTH_RIGHT_SPLINE: // Use only the extension of South spline to East (right side)
      return ( (// cell North side
		PolynomLineIntegration2(Node[CellIndexI+1][JNl+1],
					Node[CellIndexI  ][JNl+1],
					Cell[CellIndexI][JCl].Xc.x,Cell[CellIndexI][JCl].Xc.y,OrderX,OrderY) + 
		// cell West side
		PolynomLineIntegration2(Node[CellIndexI  ][JNl+1],
					Node[CellIndexI  ][JNl  ],
					Cell[CellIndexI][JCl].Xc.x,Cell[CellIndexI][JCl].Xc.y,OrderX,OrderY) +   
		// cell South side
		ExtendEast_BndSouthSpline.PolynomOrderIntegration(Node[CellIndexI  ][JNl  ], Node[CellIndexI+1][JNl  ],
								  Cell[CellIndexI][JCl].Xc, 15, OrderX,OrderY) +
		// cell East side
		PolynomLineIntegration2(Node[CellIndexI+1][JNl  ],
					Node[CellIndexI+1][JNl+1],
					Cell[CellIndexI][JCl].Xc.x,Cell[CellIndexI][JCl].Xc.y,OrderX,OrderY) )
	       /(OrderX+1)/Cell[CellIndexI][JCl].A);  // Division by (OrderX+1) & Division by A

      
    case EXTEND_W_SOUTH_LEFT_SPLINE: // Use only the extension of South spline to West (left side)
      return ( (// cell North side
		PolynomLineIntegration2(Node[CellIndexI+1][JNl+1],
					Node[CellIndexI  ][JNl+1],
					Cell[CellIndexI][JCl].Xc.x, Cell[CellIndexI][JCl].Xc.y, OrderX, OrderY) + 
		// cell West side
		PolynomLineIntegration2(Node[CellIndexI  ][JNl+1],
					Node[CellIndexI  ][JNl  ],
					Cell[CellIndexI][JCl].Xc.x, Cell[CellIndexI][JCl].Xc.y, OrderX, OrderY) +
		// cell South side
		ExtendWest_BndSouthSpline.PolynomOrderIntegration(Node[CellIndexI  ][JNl  ], Node[CellIndexI+1][JNl  ],
								  Cell[CellIndexI][JCl].Xc, 15, OrderX, OrderY) +      
		// cell East side
		PolynomLineIntegration2(Node[CellIndexI+1][JNl  ],
					Node[CellIndexI+1][JNl+1],
					Cell[CellIndexI][JCl].Xc.x, Cell[CellIndexI][JCl].Xc.y, OrderX, OrderY) )
	       /(OrderX+1)/Cell[CellIndexI][JCl].A);  // Division by (OrderX+1) & Division by A


    case EXTEND_E_SOUTH_LEFT_SPLINE: // Use only the extension of South spline to East (left side)
      return ( (// cell North side
		ExtendEast_BndSouthSpline.PolynomOrderIntegration(Node[CellIndexI+1][JNl  ], Node[CellIndexI  ][JNl  ],
								  Cell[CellIndexI][JCl-1].Xc, 15, OrderX,OrderY) +
		// cell West side
		PolynomLineIntegration2(Node[CellIndexI  ][JNl  ],
					Node[CellIndexI  ][JNl-1],
					Cell[CellIndexI][JCl-1].Xc.x, Cell[CellIndexI][JCl-1].Xc.y, OrderX, OrderY) +   
		// cell South side
		PolynomLineIntegration2(Node[CellIndexI  ][JNl-1],
					Node[CellIndexI+1][JNl-1],
					Cell[CellIndexI][JCl-1].Xc.x, Cell[CellIndexI][JCl-1].Xc.y, OrderX, OrderY) + 
		// cell East side
		PolynomLineIntegration2(Node[CellIndexI+1][JNl-1],
					Node[CellIndexI+1][JNl  ],
					Cell[CellIndexI][JCl-1].Xc.x, Cell[CellIndexI][JCl-1].Xc.y, OrderX, OrderY) )
	       /(OrderX+1)/Cell[CellIndexI][JCl-1].A);  // Division by (OrderX+1) & Division by A	       


    case EXTEND_N_EAST_RIGHT_SPLINE: // Use only the extension of East spline to North (right side)
      return ( (// cell North side
		PolynomLineIntegration2(Node[INu  ][CellIndexJ+1],
					Node[INu-1][CellIndexJ+1],
					Cell[ICu][CellIndexJ].Xc.x, Cell[ICu][CellIndexJ].Xc.y, OrderX, OrderY) + 
		// cell West side
		PolynomLineIntegration2(Node[INu-1][CellIndexJ+1],
					Node[INu-1][CellIndexJ  ],
					Cell[ICu][CellIndexJ].Xc.x, Cell[ICu][CellIndexJ].Xc.y, OrderX, OrderY) +
		// cell South side
		PolynomLineIntegration2(Node[INu-1][CellIndexJ  ],
					Node[INu  ][CellIndexJ  ],
					Cell[ICu][CellIndexJ].Xc.x, Cell[ICu][CellIndexJ].Xc.y, OrderX, OrderY) + 
		// cell East side
		ExtendNorth_BndEastSpline.PolynomOrderIntegration(Node[INu  ][CellIndexJ  ], Node[INu  ][CellIndexJ+1],
								  Cell[ICu][CellIndexJ].Xc, 15, OrderX, OrderY) ) 
	       /(OrderX+1)/Cell[ICu][CellIndexJ].A);  // Division by (OrderX+1) & Division by A


    case EXTEND_S_EAST_RIGHT_SPLINE: // Use only the extension of East spline to South (right side)
      return ( (// cell North side
		PolynomLineIntegration2(Node[INu+1][CellIndexJ+1],
					Node[INu  ][CellIndexJ+1],
					Cell[ICu+1][CellIndexJ].Xc.x, Cell[ICu+1][CellIndexJ].Xc.y, OrderX, OrderY) +
		// cell West side
		ExtendSouth_BndEastSpline.PolynomOrderIntegration(Node[INu  ][CellIndexJ+1], Node[INu  ][CellIndexJ  ],
								  Cell[ICu+1][CellIndexJ].Xc, 15, OrderX, OrderY) + 
		// cell South side
		PolynomLineIntegration2(Node[INu  ][CellIndexJ  ],
					Node[INu+1][CellIndexJ  ],
					Cell[ICu+1][CellIndexJ].Xc.x, Cell[ICu+1][CellIndexJ].Xc.y, OrderX, OrderY) + 
		// cell East side
		PolynomLineIntegration2(Node[INu+1][CellIndexJ  ],
					Node[INu+1][CellIndexJ+1],
					Cell[ICu+1][CellIndexJ].Xc.x, Cell[ICu+1][CellIndexJ].Xc.y, OrderX, OrderY) )
	       /(OrderX+1)/Cell[ICu+1][CellIndexJ].A);  // Division by (OrderX+1) & Division by A


    case EXTEND_N_EAST_LEFT_SPLINE: // Use only the extension of East spline to North (left side)
      return ( (// cell North side
		PolynomLineIntegration2(Node[INu+1][CellIndexJ+1],
					Node[INu  ][CellIndexJ+1],
					Cell[ICu+1][CellIndexJ].Xc.x, Cell[ICu+1][CellIndexJ].Xc.y, OrderX, OrderY) +
		// cell West side
		ExtendNorth_BndEastSpline.PolynomOrderIntegration(Node[INu  ][CellIndexJ+1], Node[INu  ][CellIndexJ  ],
								  Cell[ICu+1][CellIndexJ].Xc, 15, OrderX, OrderY) + 
		// cell South side
		PolynomLineIntegration2(Node[INu  ][CellIndexJ  ],
					Node[INu+1][CellIndexJ  ],
					Cell[ICu+1][CellIndexJ].Xc.x, Cell[ICu+1][CellIndexJ].Xc.y, OrderX, OrderY) + 
		// cell East side
		PolynomLineIntegration2(Node[INu+1][CellIndexJ  ],
					Node[INu+1][CellIndexJ+1],
					Cell[ICu+1][CellIndexJ].Xc.x, Cell[ICu+1][CellIndexJ].Xc.y, OrderX, OrderY) ) 
	       /(OrderX+1)/Cell[ICu+1][CellIndexJ].A);  // Division by (OrderX+1) & Division by A


    case EXTEND_S_EAST_LEFT_SPLINE: // Use only the extension of East spline to South (left side)
      return ( (// cell North side
		PolynomLineIntegration2(Node[INu  ][CellIndexJ+1],
					Node[INu-1][CellIndexJ+1],
					Cell[ICu][CellIndexJ].Xc.x, Cell[ICu][CellIndexJ].Xc.y, OrderX, OrderY) + 
		// cell West side
		PolynomLineIntegration2(Node[INu-1][CellIndexJ+1],
					Node[INu-1][CellIndexJ  ],
					Cell[ICu][CellIndexJ].Xc.x, Cell[ICu][CellIndexJ].Xc.y, OrderX, OrderY) +
		// cell South side
		PolynomLineIntegration2(Node[INu-1][CellIndexJ  ],
					Node[INu  ][CellIndexJ  ],
					Cell[ICu][CellIndexJ].Xc.x, Cell[ICu][CellIndexJ].Xc.y, OrderX, OrderY) + 
		// cell East side
		ExtendSouth_BndEastSpline.PolynomOrderIntegration(Node[INu  ][CellIndexJ  ], Node[INu  ][CellIndexJ+1],
								  Cell[ICu][CellIndexJ].Xc, 15, OrderX, OrderY) )
	       /(OrderX+1)/Cell[ICu][CellIndexJ].A);  // Division by (OrderX+1) & Division by A 


    case EXTEND_N_WEST_RIGHT_SPLINE: // Use only the extension of West spline to North (right side)
      return ( (// cell North side
		PolynomLineIntegration2(Node[INl  ][CellIndexJ+1],
					Node[INl-1][CellIndexJ+1],
					Cell[ICl-1][CellIndexJ].Xc.x, Cell[ICl-1][CellIndexJ].Xc.y, OrderX, OrderY) + 
		// cell West side
		PolynomLineIntegration2(Node[INl-1][CellIndexJ+1],
					Node[INl-1][CellIndexJ  ],
					Cell[ICl-1][CellIndexJ].Xc.x, Cell[ICl-1][CellIndexJ].Xc.y, OrderX, OrderY) +
		// cell South side
		PolynomLineIntegration2(Node[INl-1][CellIndexJ  ],
					Node[INl  ][CellIndexJ  ],
					Cell[ICl-1][CellIndexJ].Xc.x, Cell[ICl-1][CellIndexJ].Xc.y, OrderX, OrderY) +
		// cell East side
		ExtendNorth_BndWestSpline.PolynomOrderIntegration(Node[INl  ][CellIndexJ  ], Node[INl  ][CellIndexJ+1],
								  Cell[ICl-1][CellIndexJ].Xc, 15, OrderX, OrderY) )
	       /(OrderX+1)/Cell[ICl-1][CellIndexJ].A);  // Division by (OrderX+1) & Division by A


    case EXTEND_S_WEST_RIGHT_SPLINE: // Use only the extension of West spline to South (right side)
      return ( (// cell North side
		PolynomLineIntegration2(Node[INl+1][CellIndexJ+1],
					Node[INl  ][CellIndexJ+1],
					Cell[ICl][CellIndexJ].Xc.x, Cell[ICl][CellIndexJ].Xc.y, OrderX, OrderY) +
		// cell West side
		ExtendSouth_BndWestSpline.PolynomOrderIntegration(Node[INl  ][CellIndexJ+1], Node[INl  ][CellIndexJ  ],
								  Cell[ICl][CellIndexJ].Xc, 15, OrderX, OrderY) + 
		// cell South side
		PolynomLineIntegration2(Node[INl  ][CellIndexJ  ],
					Node[INl+1][CellIndexJ  ],
					Cell[ICl][CellIndexJ].Xc.x, Cell[ICl][CellIndexJ].Xc.y, OrderX, OrderY) + 
		// cell East side
		PolynomLineIntegration2(Node[INl+1][CellIndexJ  ],
					Node[INl+1][CellIndexJ+1],
					Cell[ICl][CellIndexJ].Xc.x, Cell[ICl][CellIndexJ].Xc.y, OrderX, OrderY) )
	       /(OrderX+1)/Cell[ICl][CellIndexJ].A);  // Division by (OrderX+1) & Division by A


    case EXTEND_N_WEST_LEFT_SPLINE: // Use only the extension of West spline to North (left side)
      return ( (// cell North side
		PolynomLineIntegration2(Node[INl+1][CellIndexJ+1],
					Node[INl  ][CellIndexJ+1],
					Cell[ICl][CellIndexJ].Xc.x, Cell[ICl][CellIndexJ].Xc.y, OrderX, OrderY) +
		// cell West side
		ExtendNorth_BndWestSpline.PolynomOrderIntegration(Node[INl  ][CellIndexJ+1], Node[INl  ][CellIndexJ  ],
								  Cell[ICl][CellIndexJ].Xc, 15, OrderX, OrderY) + 
		// cell South side
		PolynomLineIntegration2(Node[INl  ][CellIndexJ  ],
					Node[INl+1][CellIndexJ  ],
					Cell[ICl][CellIndexJ].Xc.x, Cell[ICl][CellIndexJ].Xc.y, OrderX, OrderY) + 
		// cell East side
		PolynomLineIntegration2(Node[INl+1][CellIndexJ  ],
					Node[INl+1][CellIndexJ+1],
					Cell[ICl][CellIndexJ].Xc.x, Cell[ICl][CellIndexJ].Xc.y, OrderX, OrderY) )
	       /(OrderX+1)/Cell[ICl][CellIndexJ].A);  // Division by (OrderX+1) & Division by A


    case EXTEND_S_WEST_LEFT_SPLINE: // Use only the extension of West spline to South (left side)
      return ( (// cell North side
		PolynomLineIntegration2(Node[INl  ][CellIndexJ+1],
					Node[INl-1][CellIndexJ+1],
					Cell[ICl-1][CellIndexJ].Xc.x, Cell[ICl-1][CellIndexJ].Xc.y, OrderX, OrderY) + 
		// cell West side
		PolynomLineIntegration2(Node[INl-1][CellIndexJ+1],
					Node[INl-1][CellIndexJ  ],
					Cell[ICl-1][CellIndexJ].Xc.x, Cell[ICl-1][CellIndexJ].Xc.y, OrderX, OrderY) +
		// cell South side
		PolynomLineIntegration2(Node[INl-1][CellIndexJ  ],
					Node[INl  ][CellIndexJ  ],
					Cell[ICl-1][CellIndexJ].Xc.x, Cell[ICl-1][CellIndexJ].Xc.y, OrderX, OrderY) +
		// cell East side
		ExtendSouth_BndWestSpline.PolynomOrderIntegration(Node[INl  ][CellIndexJ  ], Node[INl  ][CellIndexJ+1],
								  Cell[ICl-1][CellIndexJ].Xc, 15, OrderX, OrderY) )
	       /(OrderX+1)/Cell[ICl-1][CellIndexJ].A);  // Division by (OrderX+1) & Division by A	       


    case CORNER_NORTH_EXTEND_N_WEST_SPLINES: // Use the North spline and the extension of West spline to North
      throw runtime_error("Grid2D_Quad_Block_HO::GeometricMoment_GhostCell_CurvedBoundaries() ERROR! Case CORNER_NORTH_EXTEND_N_WEST_SPLINES has not been encountered so far. Code and test this case!");
      return 0.0;


    case CORNER_EXTEND_N_WEST_EXTEND_W_NORTH_SPLINES: // Use the North extension of West and the West extension of North
      throw runtime_error("Grid2D_Quad_Block_HO::GeometricMoment_GhostCell_CurvedBoundaries() ERROR! Case CORNER_EXTEND_N_WEST_EXTEND_W_NORTH_SPLINES has not been encountered so far. Code and test this case!");
      return 0.0;


    case CORNER_EXTEND_W_NORTH_WEST_SPLINES: // Use the West extension of North and the West spline
      throw runtime_error("Grid2D_Quad_Block_HO::GeometricMoment_GhostCell_CurvedBoundaries() ERROR! Case CORNER_EXTEND_W_NORTH_WEST_SPLINES has not been encountered so far. Code and test this case!");
      return 0.0;


    case CORNER_WEST_EXTEND_W_SOUTH_SPLINES: // Use the West spline and the West extension of South
      throw runtime_error("Grid2D_Quad_Block_HO::GeometricMoment_GhostCell_CurvedBoundaries() ERROR! Case CORNER_WEST_EXTEND_W_SOUTH_SPLINES has not been encountered so far. Code and test this case!");
      return 0.0;


    case CORNER_EXTEND_W_SOUTH_EXTEND_S_WEST_SPLINES: // Use the West extension of South and the South extension of West
      throw runtime_error("Grid2D_Quad_Block_HO::GeometricMoment_GhostCell_CurvedBoundaries() ERROR! Case CORNER_EXTEND_W_SOUTH_EXTEND_S_WEST_SPLINES has not been encountered so far. Code and test this case!");
      return 0.0;


    case CORNER_EXTEND_S_WEST_SOUTH_SPLINES: // Use the South extension of West and the South spline
      throw runtime_error("Grid2D_Quad_Block_HO::GeometricMoment_GhostCell_CurvedBoundaries() ERROR! Case CORNER_EXTEND_S_WEST_SOUTH_SPLINES has not been encountered so far. Code and test this case!");
      return 0.0;


    case CORNER_SOUTH_EXTEND_S_EAST_SPLINES: // Use the South spline and the South extension of East
      throw runtime_error("Grid2D_Quad_Block_HO::GeometricMoment_GhostCell_CurvedBoundaries() ERROR! Case CORNER_SOUTH_EXTEND_S_EAST_SPLINES has not been encountered so far. Code and test this case!");
      return 0.0;


    case CORNER_EXTEND_S_EAST_EXTEND_E_SOUTH_SPLINES: // Use the South extension of East and the East extension of South
      throw runtime_error("Grid2D_Quad_Block_HO::GeometricMoment_GhostCell_CurvedBoundaries() ERROR! Case CORNER_EXTEND_S_EAST_EXTEND_E_SOUTH_SPLINES has not been encountered so far. Code and test this case!");
      return 0.0;


    case CORNER_EXTEND_E_SOUTH_EAST_SPLINES: // Use the East extension of South and the East spline
      throw runtime_error("Grid2D_Quad_Block_HO::GeometricMoment_GhostCell_CurvedBoundaries() ERROR! Case CORNER_EXTEND_E_SOUTH_EAST_SPLINES has not been encountered so far. Code and test this case!");
      return 0.0;


    case CORNER_EAST_EXTEND_E_NORTH_SPLINES: // Use the East spline and the East extension of North
      throw runtime_error("Grid2D_Quad_Block_HO::GeometricMoment_GhostCell_CurvedBoundaries() ERROR! Case CORNER_EAST_EXTEND_E_NORTH_SPLINES has not been encountered so far. Code and test this case!");
      return 0.0;


    case CORNER_EXTEND_E_NORTH_EXTEND_N_EAST_SPLINES: // Use the East extention of North and North extension of East
      throw runtime_error("Grid2D_Quad_Block_HO::GeometricMoment_GhostCell_CurvedBoundaries() ERROR! Case CORNER_EXTEND_E_NORTH_EXTEND_N_EAST_SPLINES has not been encountered so far. Code and test this case!");
      return 0.0;


    case CORNER_EXTEND_N_EAST_NORTH_SPLINES: // Use the North extension of East and the North spline
      throw runtime_error("Grid2D_Quad_Block_HO::GeometricMoment_GhostCell_CurvedBoundaries() ERROR! Case CORNER_EXTEND_N_EAST_NORTH_SPLINES has not been encountered so far. Code and test this case!");
      return 0.0;

    default:
      return 0.0;
    } // endswitch

  }// endif
}

/*!
 * Updates the cell information for the quadrilateral
 * mesh block. Different levels of update are possible
 * based on the corresponding flag values.
 * If all the flags are OFF this subroutine does nothing.
 */
void Grid2D_Quad_Block_HO::Update_Cells(void) {

  // Decide which level of update is required
  if (InteriorMeshUpdate == ON && GhostCellsUpdate == ON ){
    // Update all the cells and the boundary spline information
    Update_All_Cells();
    // Mark the current geometry different than the previous one
    New_Global_Geometry_State();
    return;
  }
  if (InteriorMeshUpdate == ON) {
    // Update only the information of the interior cells and the boundary spline information
    Update_Interior_Cells();
    // Mark the interior geometry different than the previous one
    New_Interior_Geometry_State();
  } 
  if (GhostCellsUpdate == ON) {
    // Update only the information of the ghost cells
    Update_Ghost_Cells();
    // Mark the geometry in ghost layers different than the previous one
    New_Ghost_Geometry_State();
  } else if (CornerGhostCellsUpdate == ON){
    // Update only the information of the corner ghost cells
    Update_Corner_Ghost_Cells();
    // Mark the corner geometry different than the previous one
    New_Corner_Geometry_State();
  }

}

/*!
 * Updates the cell information for the quadrilateral mesh block.
 * (i.e. all cells).
 */
void Grid2D_Quad_Block_HO::Update_All_Cells(void) {

  int i,j;
  bool CurvedNorthBnd, CurvedSouthBnd, CurvedEastBnd, CurvedWestBnd;
  bool Curved_Extend_N_WestBnd, Curved_Extend_S_WestBnd, Curved_Extend_N_EastBnd, Curved_Extend_S_EastBnd;
  bool Curved_Extend_E_SouthBnd, Curved_Extend_W_SouthBnd, Curved_Extend_E_NorthBnd, Curved_Extend_W_NorthBnd;

  int NumGQPsPerSubinterval(NumGQP);

#ifdef MODIFY_NUMBER_OF_FLUX_CALCULATION_POINTS_AT_BOUNDARIES
  NumGQPsPerSubinterval = NumGQP + 1;
#endif

  // Update cell information assuming straight boundaries (i.e. every edge of the cell is a line segment)
  for ( j = JCl-Nghost ; j <= JCu+Nghost ; ++j) {
    for ( i = ICl-Nghost ; i <= ICu+Nghost ; ++i) {
      Update_Cell(i,j);
    } /* endfor */
  } /* endfor */

  if ( !CheckExistenceOfCurvedBoundaries() ){
    // There is no need for curved boundary representation

    // Confirm the update
    Confirm_Mesh_Update_Everywhere();    
    return;
    
  } else {

    /* Recompute the geometric properties of those interior cells that are near a curved boundary.
       Obs1. The "SPLINE2D_LINEAR" is also considered a curved boundary because it might have a 
       sharp point between two nodes.
    */

    // Determine which boundaries are curved
    CurvedNorthBnd = IsNorthBoundaryCurved();
    CurvedSouthBnd = IsSouthBoundaryCurved();
    CurvedEastBnd  = IsEastBoundaryCurved();
    CurvedWestBnd  = IsWestBoundaryCurved();
    Curved_Extend_N_WestBnd = IsNorthExtendWestBoundaryCurved();
    Curved_Extend_S_WestBnd = IsSouthExtendWestBoundaryCurved();
    Curved_Extend_N_EastBnd = IsNorthExtendEastBoundaryCurved();
    Curved_Extend_S_EastBnd = IsSouthExtendEastBoundaryCurved();
    Curved_Extend_E_SouthBnd= IsEastExtendSouthBoundaryCurved();
    Curved_Extend_W_SouthBnd= IsWestExtendSouthBoundaryCurved();
    Curved_Extend_E_NorthBnd= IsEastExtendNorthBoundaryCurved();
    Curved_Extend_W_NorthBnd= IsWestExtendNorthBoundaryCurved();  


    // Generate BndSplineInfo(s) if necessary
    // Determine the geometric properties along the splines (e.g. Gauss Quadrature point locations, normals, etc.)
    Update_SplineInfos();

    // Determine the geometric properties of the cell (e.g. centroid, area, geometric moments)
    // and the geometric properties along the splines (e.g. Gauss Quadrature point locations, normals, and
    // spline segment length at each cell)

    // Check the North boundary
    if (CurvedNorthBnd){

      for(i=ICl+1; i<=ICu-1; ++i){
	// Update the interior north cells
	Cell[i][JCu].A = area_CurvedBoundaries(i,NORTH_SPLINE);
	Cell[i][JCu].Xc = centroid_CurvedBoundaries(i,NORTH_SPLINE);
	ComputeGeometricCoefficients_CurvedBoundaries(i,JCu,NORTH_SPLINE);

	// Update the ghost cells
	Cell[i][JCu+1].A = area_GhostCell_CurvedBoundaries(i,NORTH_SPLINE);
	Cell[i][JCu+1].Xc = centroid_GhostCell_CurvedBoundaries(i,NORTH_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(i,JCu+1,NORTH_SPLINE);
      }

      // Update the first ghost cell on the North boundary
      if (Curved_Extend_N_WestBnd){
	// Both splines are curved
	Cell[ICl][JCu+1].A = area_GhostCell_CurvedBoundaries(ICl,CORNER_NORTH_EXTEND_N_WEST_SPLINES);
	Cell[ICl][JCu+1].Xc = centroid_GhostCell_CurvedBoundaries(ICl,CORNER_NORTH_EXTEND_N_WEST_SPLINES);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICl,JCu+1,CORNER_NORTH_EXTEND_N_WEST_SPLINES);
      } else {
	// Only North boundary is curved
	Cell[ICl][JCu+1].A = area_GhostCell_CurvedBoundaries(ICl,NORTH_SPLINE);
	Cell[ICl][JCu+1].Xc = centroid_GhostCell_CurvedBoundaries(ICl,NORTH_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICl,JCu+1,NORTH_SPLINE);
      }

      // Update the last ghost cell on the North boundary
      if (Curved_Extend_N_EastBnd){
	// Both splines are curved
	Cell[ICu][JCu+1].A = area_GhostCell_CurvedBoundaries(ICu,CORNER_EXTEND_N_EAST_NORTH_SPLINES);
	Cell[ICu][JCu+1].Xc = centroid_GhostCell_CurvedBoundaries(ICu,CORNER_EXTEND_N_EAST_NORTH_SPLINES);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICu,JCu+1,CORNER_EXTEND_N_EAST_NORTH_SPLINES);
      } else {
	// Only North boundary is curved
	Cell[ICu][JCu+1].A = area_GhostCell_CurvedBoundaries(ICu,NORTH_SPLINE);
	Cell[ICu][JCu+1].Xc = centroid_GhostCell_CurvedBoundaries(ICu,NORTH_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICu,JCu+1,NORTH_SPLINE);
      }

      // Check the North-West Corner
      if(CurvedWestBnd){
	// both splines are curved
	Cell[ICl][JCu].A = area_CurvedBoundaries(ICl,CORNER_NORTH_WEST_SPLINES);
	Cell[ICl][JCu].Xc = centroid_CurvedBoundaries(ICl,CORNER_NORTH_WEST_SPLINES);
	ComputeGeometricCoefficients_CurvedBoundaries(ICl,JCu,CORNER_NORTH_WEST_SPLINES);
      } else {
	Cell[ICl][JCu].A = area_CurvedBoundaries(ICl,NORTH_SPLINE);
	Cell[ICl][JCu].Xc = centroid_CurvedBoundaries(ICl,NORTH_SPLINE);
	ComputeGeometricCoefficients_CurvedBoundaries(ICl,JCu,NORTH_SPLINE);
      }

      // Check the North-East Corner
      if (CurvedEastBnd){
	// both spline are curved
	Cell[ICu][JCu].A = area_CurvedBoundaries(ICu,CORNER_NORTH_EAST_SPLINES);
	Cell[ICu][JCu].Xc = centroid_CurvedBoundaries(ICu,CORNER_NORTH_EAST_SPLINES);
	ComputeGeometricCoefficients_CurvedBoundaries(ICu,JCu,CORNER_NORTH_EAST_SPLINES);
      } else {
	Cell[ICu][JCu].A = area_CurvedBoundaries(ICu,NORTH_SPLINE);
	Cell[ICu][JCu].Xc = centroid_CurvedBoundaries(ICu,NORTH_SPLINE);
	ComputeGeometricCoefficients_CurvedBoundaries(ICu,JCu,NORTH_SPLINE);
      }
    } // endif (North Boundary)

    // Check the South boundary
    if (CurvedSouthBnd){

      for(i=ICl+1; i<=ICu-1; ++i){
	// Update the interior south cells
	Cell[i][JCl].A = area_CurvedBoundaries(i,SOUTH_SPLINE);
	Cell[i][JCl].Xc = centroid_CurvedBoundaries(i,SOUTH_SPLINE);
	ComputeGeometricCoefficients_CurvedBoundaries(i,JCl,SOUTH_SPLINE);

	// Update the ghost cells
	Cell[i][JCl-1].A = area_GhostCell_CurvedBoundaries(i,SOUTH_SPLINE);
	Cell[i][JCl-1].Xc = centroid_GhostCell_CurvedBoundaries(i,SOUTH_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(i,JCl-1,SOUTH_SPLINE);
      }

      // Update the first ghost cell on the South boundary
      if (Curved_Extend_S_WestBnd){
	// Both splines are curved
	Cell[ICl][JCl-1].A = area_GhostCell_CurvedBoundaries(ICl,CORNER_EXTEND_S_WEST_SOUTH_SPLINES);
	Cell[ICl][JCl-1].Xc = centroid_GhostCell_CurvedBoundaries(ICl,CORNER_EXTEND_S_WEST_SOUTH_SPLINES);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICl,JCl-1,CORNER_EXTEND_S_WEST_SOUTH_SPLINES);
      } else {
	// Only South spline is curved
	Cell[ICl][JCl-1].A = area_GhostCell_CurvedBoundaries(ICl,SOUTH_SPLINE);
	Cell[ICl][JCl-1].Xc = centroid_GhostCell_CurvedBoundaries(ICl,SOUTH_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICl,JCl-1,SOUTH_SPLINE);
      }

      // Update the last ghost cell on the South boundary
      if (Curved_Extend_S_EastBnd){
	// Both splines are curved
	Cell[ICu][JCl-1].A = area_GhostCell_CurvedBoundaries(ICu,CORNER_SOUTH_EXTEND_S_EAST_SPLINES);
	Cell[ICu][JCl-1].Xc = centroid_GhostCell_CurvedBoundaries(ICu,CORNER_SOUTH_EXTEND_S_EAST_SPLINES);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICu,JCl-1,CORNER_SOUTH_EXTEND_S_EAST_SPLINES);
      } else {
	// Only South spline is curved
	Cell[ICu][JCl-1].A = area_GhostCell_CurvedBoundaries(ICu,SOUTH_SPLINE);
	Cell[ICu][JCl-1].Xc = centroid_GhostCell_CurvedBoundaries(ICu,SOUTH_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICu,JCl-1,SOUTH_SPLINE);
      }

      // Check the South-West Corner
      if(CurvedWestBnd){
	// both splines are curved
	Cell[ICl][JCl].A = area_CurvedBoundaries(ICl,CORNER_SOUTH_WEST_SPLINES);
	Cell[ICl][JCl].Xc = centroid_CurvedBoundaries(ICl,CORNER_SOUTH_WEST_SPLINES);
	ComputeGeometricCoefficients_CurvedBoundaries(ICl,JCl,CORNER_SOUTH_WEST_SPLINES);
      } else {
	Cell[ICl][JCl].A = area_CurvedBoundaries(ICl,SOUTH_SPLINE);
	Cell[ICl][JCl].Xc = centroid_CurvedBoundaries(ICl,SOUTH_SPLINE);
	ComputeGeometricCoefficients_CurvedBoundaries(ICl,JCl,SOUTH_SPLINE);
      }

      // Check the South-East Corner
      if(CurvedEastBnd){
	// both spline are curved
	Cell[ICu][JCl].A = area_CurvedBoundaries(ICu,CORNER_SOUTH_EAST_SPLINES);
	Cell[ICu][JCl].Xc = centroid_CurvedBoundaries(ICu,CORNER_SOUTH_EAST_SPLINES);
	ComputeGeometricCoefficients_CurvedBoundaries(ICu,JCl,CORNER_SOUTH_EAST_SPLINES);
      } else {
	Cell[ICu][JCl].A = area_CurvedBoundaries(ICu,SOUTH_SPLINE);
	Cell[ICu][JCl].Xc = centroid_CurvedBoundaries(ICu,SOUTH_SPLINE);
	ComputeGeometricCoefficients_CurvedBoundaries(ICu,JCl,SOUTH_SPLINE);
      }
    } // endif (South Boundary)

    // Check the East boundary
    if (CurvedEastBnd){

      for(j=JCl+1; j<=JCu-1; ++j){
	// Update the interior east cells
	Cell[ICu][j].A = area_CurvedBoundaries(j,EAST_SPLINE);
	Cell[ICu][j].Xc = centroid_CurvedBoundaries(j,EAST_SPLINE);
	ComputeGeometricCoefficients_CurvedBoundaries(ICu,j,EAST_SPLINE);

	// Update the ghost cells
	Cell[ICu+1][j].A = area_GhostCell_CurvedBoundaries(j,EAST_SPLINE);
	Cell[ICu+1][j].Xc = centroid_GhostCell_CurvedBoundaries(j,EAST_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICu+1,j,EAST_SPLINE);
      }

      // Update the first ghost cell on the East boundary
      if (Curved_Extend_E_SouthBnd){
	// Both splines are curved
	Cell[ICu+1][JCl].A = area_GhostCell_CurvedBoundaries(JCl,CORNER_EXTEND_E_SOUTH_EAST_SPLINES);
	Cell[ICu+1][JCl].Xc = centroid_GhostCell_CurvedBoundaries(JCl,CORNER_EXTEND_E_SOUTH_EAST_SPLINES);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICu+1,JCl,CORNER_EXTEND_E_SOUTH_EAST_SPLINES);
      } else {
	// Only East spline is curved
	Cell[ICu+1][JCl].A = area_GhostCell_CurvedBoundaries(JCl,EAST_SPLINE);
	Cell[ICu+1][JCl].Xc = centroid_GhostCell_CurvedBoundaries(JCl,EAST_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICu+1,JCl,EAST_SPLINE);
      }

      // Update the last ghost cell on the East boundary
      if (Curved_Extend_E_NorthBnd){
	// Both splines are curved
	Cell[ICu+1][JCu].A = area_GhostCell_CurvedBoundaries(JCu,CORNER_EAST_EXTEND_E_NORTH_SPLINES);
	Cell[ICu+1][JCu].Xc = centroid_GhostCell_CurvedBoundaries(JCu,CORNER_EAST_EXTEND_E_NORTH_SPLINES);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICu+1,JCu,CORNER_EAST_EXTEND_E_NORTH_SPLINES);
      }	else {
	// Only East spline is curved
	Cell[ICu+1][JCu].A = area_GhostCell_CurvedBoundaries(JCu,EAST_SPLINE);
	Cell[ICu+1][JCu].Xc = centroid_GhostCell_CurvedBoundaries(JCu,EAST_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICu+1,JCu,EAST_SPLINE);
      }

      // Re-check the North-East Corner
      if (!CurvedNorthBnd){
	Cell[ICu][JCu].A = area_CurvedBoundaries(JCu,EAST_SPLINE);
	Cell[ICu][JCu].Xc = centroid_CurvedBoundaries(JCu,EAST_SPLINE);
	ComputeGeometricCoefficients_CurvedBoundaries(ICu,JCu,EAST_SPLINE);
      }

      // Re-check the South-East Corner
      if(!CurvedSouthBnd){
	Cell[ICu][JCl].A = area_CurvedBoundaries(JCl,EAST_SPLINE);
	Cell[ICu][JCl].Xc = centroid_CurvedBoundaries(JCl,EAST_SPLINE);
	ComputeGeometricCoefficients_CurvedBoundaries(ICu,JCl,EAST_SPLINE);
      }      
    } // endif (East Boundary)

    // Check the West boundary
    if (CurvedWestBnd){

      for(j=JCl+1; j<=JCu-1; ++j){
	// Update the interior west cells
	Cell[ICl][j].A = area_CurvedBoundaries(j,WEST_SPLINE);
	Cell[ICl][j].Xc = centroid_CurvedBoundaries(j,WEST_SPLINE);
	ComputeGeometricCoefficients_CurvedBoundaries(ICl,j,WEST_SPLINE);

	// Update the ghost cells
	Cell[ICl-1][j].A = area_GhostCell_CurvedBoundaries(j,WEST_SPLINE);
	Cell[ICl-1][j].Xc = centroid_GhostCell_CurvedBoundaries(j,WEST_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICl-1,j,WEST_SPLINE);
      }

      // Update the first ghost cell on the West boundary
      if (Curved_Extend_W_SouthBnd){
	// Both splines are curved
	Cell[ICl-1][JCl].A = area_GhostCell_CurvedBoundaries(JCl,CORNER_WEST_EXTEND_W_SOUTH_SPLINES);
	Cell[ICl-1][JCl].Xc = centroid_GhostCell_CurvedBoundaries(JCl,CORNER_WEST_EXTEND_W_SOUTH_SPLINES);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICl-1,JCl,CORNER_WEST_EXTEND_W_SOUTH_SPLINES);
      } else {
	// Only West spline is curved
	Cell[ICl-1][JCl].A = area_GhostCell_CurvedBoundaries(JCl,WEST_SPLINE);
	Cell[ICl-1][JCl].Xc = centroid_GhostCell_CurvedBoundaries(JCl,WEST_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICl-1,JCl,WEST_SPLINE);
      }

      // Update the last ghost cell on the West boundary
      if (Curved_Extend_W_NorthBnd){
	// Both splines are curved
	Cell[ICl-1][JCu].A = area_GhostCell_CurvedBoundaries(JCu,CORNER_EXTEND_W_NORTH_WEST_SPLINES);
	Cell[ICl-1][JCu].Xc = centroid_GhostCell_CurvedBoundaries(JCu,CORNER_EXTEND_W_NORTH_WEST_SPLINES);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICl-1,JCu,CORNER_EXTEND_W_NORTH_WEST_SPLINES);
      } else {
	// Only West spline is curved
	Cell[ICl-1][JCu].A = area_GhostCell_CurvedBoundaries(JCu,WEST_SPLINE);
	Cell[ICl-1][JCu].Xc = centroid_GhostCell_CurvedBoundaries(JCu,WEST_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICl-1,JCu,WEST_SPLINE);
      }

      // Re-check the North-West Corner
      if(!CurvedNorthBnd){
	Cell[ICl][JCu].A = area_CurvedBoundaries(JCu,WEST_SPLINE);
	Cell[ICl][JCu].Xc = centroid_CurvedBoundaries(JCu,WEST_SPLINE);
	ComputeGeometricCoefficients_CurvedBoundaries(ICl,JCu,WEST_SPLINE);
      }

      // Re-check the South-West Corner
      if(!CurvedSouthBnd){
	Cell[ICl][JCl].A = area_CurvedBoundaries(JCl,WEST_SPLINE);
	Cell[ICl][JCl].Xc = centroid_CurvedBoundaries(JCl,WEST_SPLINE);
	ComputeGeometricCoefficients_CurvedBoundaries(ICl,JCl,WEST_SPLINE);
      }
    } // endif (West Boundary)

    // ===  Update ghost cells influenced by the presence of curved extension splines ====
    Update_GhostCells_Near_CurvedExtensionSplines(CurvedNorthBnd,
						  CurvedSouthBnd,
						  CurvedEastBnd,
						  CurvedWestBnd,
						  Curved_Extend_N_WestBnd, 
						  Curved_Extend_S_WestBnd, 
						  Curved_Extend_N_EastBnd, 
						  Curved_Extend_S_EastBnd, 
						  Curved_Extend_E_SouthBnd,
						  Curved_Extend_W_SouthBnd,
						  Curved_Extend_E_NorthBnd,
						  Curved_Extend_W_NorthBnd);

    // Confirm the update
    Confirm_Mesh_Update_Everywhere();
  }
}

/*!
 * Updates the cell information for the quadrilateral mesh block 
 * interior cells (i.e. no ghost cells).
 */
void Grid2D_Quad_Block_HO::Update_Interior_Cells(void) {

  int i,j;

  int NumGQPsPerSubinterval(NumGQP);
  bool CurvedNorthBnd, CurvedSouthBnd, CurvedEastBnd, CurvedWestBnd;

#ifdef MODIFY_NUMBER_OF_FLUX_CALCULATION_POINTS_AT_BOUNDARIES
  NumGQPsPerSubinterval = NumGQP + 1;
#endif

  // Update cell information assuming straight boundaries (i.e. every edge of the cell is a line segment)
  for ( j = JCl; j <= JCu ; ++j) {
    for ( i = ICl; i <= ICu; ++i) {
      Update_Cell(i,j);
    } /* endfor */
  } /* endfor */
  
  if ( !CheckExistenceOfCurvedBoundaries() ){
    // There is no need for curved boundary representation
    
    // Confirm the update
    Confirm_Interior_Mesh_Update();
    return;
    
  } else {
    
    /* Recompute the geometric properties of those interior cells that are near a curved boundary.
       Obs1. The "SPLINE2D_LINEAR" is also considered a curved boundary because it might have a 
       sharp point between two nodes.
    */
    
    // Determine which boundaries are curved
    CurvedNorthBnd = IsNorthBoundaryCurved();
    CurvedSouthBnd = IsSouthBoundaryCurved();
    CurvedEastBnd  = IsEastBoundaryCurved();
    CurvedWestBnd  = IsWestBoundaryCurved();

    // Determine the geometric properties of the cell (e.g. centroid, area, geometric moments)
    // and the geometric properties along the splines (e.g. Gauss Quadrature point locations, normals, and
    // spline segment length at each cell)
    Update_SplineInfos();

    // Check the North boundary
    if (IsNorthBoundaryCurved()){

      // Update the interior north cells
      for(i=ICl+1; i<=ICu-1; ++i){
	Cell[i][JCu].A = area_CurvedBoundaries(i,NORTH_SPLINE);
	Cell[i][JCu].Xc = centroid_CurvedBoundaries(i,NORTH_SPLINE);
	ComputeGeometricCoefficients_CurvedBoundaries(i,JCu,NORTH_SPLINE);
      }

      // Check the North-West Corner
      if(CurvedWestBnd){
	// both splines are curved
	Cell[ICl][JCu].A = area_CurvedBoundaries(ICl,CORNER_NORTH_WEST_SPLINES);
	Cell[ICl][JCu].Xc = centroid_CurvedBoundaries(ICl,CORNER_NORTH_WEST_SPLINES);
	ComputeGeometricCoefficients_CurvedBoundaries(ICl,JCu,CORNER_NORTH_WEST_SPLINES);
      } else {
	Cell[ICl][JCu].A = area_CurvedBoundaries(ICl,NORTH_SPLINE);
	Cell[ICl][JCu].Xc = centroid_CurvedBoundaries(ICl,NORTH_SPLINE);
	ComputeGeometricCoefficients_CurvedBoundaries(ICl,JCu,NORTH_SPLINE);
      }

      // Check the North-East Corner
      if (CurvedEastBnd){
	// both spline are curved
	Cell[ICu][JCu].A = area_CurvedBoundaries(ICu,CORNER_NORTH_EAST_SPLINES);
	Cell[ICu][JCu].Xc = centroid_CurvedBoundaries(ICu,CORNER_NORTH_EAST_SPLINES);
	ComputeGeometricCoefficients_CurvedBoundaries(ICu,JCu,CORNER_NORTH_EAST_SPLINES);
      } else {
	Cell[ICu][JCu].A = area_CurvedBoundaries(ICu,NORTH_SPLINE);
	Cell[ICu][JCu].Xc = centroid_CurvedBoundaries(ICu,NORTH_SPLINE);
	ComputeGeometricCoefficients_CurvedBoundaries(ICu,JCu,NORTH_SPLINE);
      }
    } // endif (North Boundary)

    // Check the South boundary
    if (IsSouthBoundaryCurved()){

      // Update the interior south cells
      for(i=ICl+1; i<=ICu-1; ++i){
	Cell[i][JCl].A = area_CurvedBoundaries(i,SOUTH_SPLINE);
	Cell[i][JCl].Xc = centroid_CurvedBoundaries(i,SOUTH_SPLINE);
	ComputeGeometricCoefficients_CurvedBoundaries(i,JCl,SOUTH_SPLINE);
      }

      // Check the South-West Corner
      if(CurvedWestBnd){
	// both splines are curved
	Cell[ICl][JCl].A = area_CurvedBoundaries(ICl,CORNER_SOUTH_WEST_SPLINES);
	Cell[ICl][JCl].Xc = centroid_CurvedBoundaries(ICl,CORNER_SOUTH_WEST_SPLINES);
	ComputeGeometricCoefficients_CurvedBoundaries(ICl,JCl,CORNER_SOUTH_WEST_SPLINES);
      } else {
	Cell[ICl][JCl].A = area_CurvedBoundaries(ICl,SOUTH_SPLINE);
	Cell[ICl][JCl].Xc = centroid_CurvedBoundaries(ICl,SOUTH_SPLINE);
	ComputeGeometricCoefficients_CurvedBoundaries(ICl,JCl,SOUTH_SPLINE);
      }

      // Check the South-East Corner
      if(CurvedEastBnd){
	// both spline are curved
	Cell[ICu][JCl].A = area_CurvedBoundaries(ICu,CORNER_SOUTH_EAST_SPLINES);
	Cell[ICu][JCl].Xc = centroid_CurvedBoundaries(ICu,CORNER_SOUTH_EAST_SPLINES);
	ComputeGeometricCoefficients_CurvedBoundaries(ICu,JCl,CORNER_SOUTH_EAST_SPLINES);
      } else {
	Cell[ICu][JCl].A = area_CurvedBoundaries(ICu,SOUTH_SPLINE);
	Cell[ICu][JCl].Xc = centroid_CurvedBoundaries(ICu,SOUTH_SPLINE);
	ComputeGeometricCoefficients_CurvedBoundaries(ICu,JCl,SOUTH_SPLINE);
      }
    } // endif (South Boundary)

    // Check the East boundary
    if (CurvedEastBnd){

      // Update the interior east cells
      for(j=JCl+1; j<=JCu-1; ++j){
	Cell[ICu][j].A = area_CurvedBoundaries(j,EAST_SPLINE);
	Cell[ICu][j].Xc = centroid_CurvedBoundaries(j,EAST_SPLINE);
	ComputeGeometricCoefficients_CurvedBoundaries(ICu,j,EAST_SPLINE);
      }

      // Re-check the North-East Corner
      if (!CurvedNorthBnd){
	Cell[ICu][JCu].A = area_CurvedBoundaries(JCu,EAST_SPLINE);
	Cell[ICu][JCu].Xc = centroid_CurvedBoundaries(JCu,EAST_SPLINE);
	ComputeGeometricCoefficients_CurvedBoundaries(ICu,JCu,EAST_SPLINE);
      }

      // Re-check the South-East Corner
      if(!CurvedSouthBnd){
	Cell[ICu][JCl].A = area_CurvedBoundaries(JCl,EAST_SPLINE);
	Cell[ICu][JCl].Xc = centroid_CurvedBoundaries(JCl,EAST_SPLINE);
	ComputeGeometricCoefficients_CurvedBoundaries(ICu,JCl,EAST_SPLINE);
      }
    } // endif (East Boundary)

    // Check the West boundary
    if (CurvedWestBnd){

      // Update the interior west cells
      for(j=JCl+1; j<=JCu-1; ++j){
	Cell[ICl][j].A = area_CurvedBoundaries(j,WEST_SPLINE);
	Cell[ICl][j].Xc = centroid_CurvedBoundaries(j,WEST_SPLINE);
	ComputeGeometricCoefficients_CurvedBoundaries(ICl,j,WEST_SPLINE);
      }

      // Re-check the North-West Corner
      if(!CurvedNorthBnd){
	Cell[ICl][JCu].A = area_CurvedBoundaries(JCu,WEST_SPLINE);
	Cell[ICl][JCu].Xc = centroid_CurvedBoundaries(JCu,WEST_SPLINE);
	ComputeGeometricCoefficients_CurvedBoundaries(ICl,JCu,WEST_SPLINE);
      }

      // Re-check the South-West Corner
      if(!CurvedSouthBnd){
	Cell[ICl][JCl].A = area_CurvedBoundaries(JCl,WEST_SPLINE);
	Cell[ICl][JCl].Xc = centroid_CurvedBoundaries(JCl,WEST_SPLINE);
	ComputeGeometricCoefficients_CurvedBoundaries(ICl,JCl,WEST_SPLINE);
      }
    } // endif (West Boundary)

    // Confirm the update
    Confirm_Interior_Mesh_Update();
  }
}

/*!
 * Updates the cell information for the quadrilateral mesh block 
 * ghost cells (i.e. no interior cells, no main boundary splines).
 * This routine updates the corner ghost cells and ONLY the extension boundary splines.
 */
void Grid2D_Quad_Block_HO::Update_Ghost_Cells(void) {

  int i,j;
  bool CurvedNorthBnd, CurvedSouthBnd, CurvedEastBnd, CurvedWestBnd;
  bool Curved_Extend_N_WestBnd, Curved_Extend_S_WestBnd, Curved_Extend_N_EastBnd, Curved_Extend_S_EastBnd;
  bool Curved_Extend_E_SouthBnd, Curved_Extend_W_SouthBnd, Curved_Extend_E_NorthBnd, Curved_Extend_W_NorthBnd;

  // Update cell information assuming straight boundaries (i.e. every edge of the cell is a line segment)
  for ( j = JCl-Nghost ; j <= JCu+Nghost ; ++j) {
    // Update the West face ghost cells
    for (i = ICl-Nghost; i <= ICl-1; ++i){
      Update_Cell(i,j);
    } // endif (i)

    // Update the East face ghost cells
    for (i = ICu+1; i <= ICu+Nghost; ++i){
      Update_Cell(i,j);
    } // endif (i)
  } // endif(j)

  for (i = ICl; i <= ICu; ++i){
    // Update the South face ghost cells
    for (j = JCl-Nghost; j<=JCl-1; ++j){
      Update_Cell(i,j);
    } // endif(j)
    
    // Update the North face ghost cells
    for (j = JCu+1; j<=JCu+Nghost; ++j){
      Update_Cell(i,j);
    } // endif(j)
  } // endif(i)
  
  if ( !CheckExistenceOfCurvedBoundaries() ){
    // There is no need for curved boundary representation
    
    // Confirm the update
    Confirm_Ghost_Cells_Update();
    return;
    
  } else {

    /* Recompute the geometric properties of those interior cells that are near a curved boundary.
       Obs1. The "SPLINE2D_LINEAR" is also considered a curved boundary because it might have a 
       sharp point between two nodes.
    */

    // Determine which boundaries are curved
    CurvedNorthBnd = IsNorthBoundaryCurved();
    CurvedSouthBnd = IsSouthBoundaryCurved();
    CurvedEastBnd  = IsEastBoundaryCurved();
    CurvedWestBnd  = IsWestBoundaryCurved();
    Curved_Extend_N_WestBnd = IsNorthExtendWestBoundaryCurved();
    Curved_Extend_S_WestBnd = IsSouthExtendWestBoundaryCurved();
    Curved_Extend_N_EastBnd = IsNorthExtendEastBoundaryCurved();
    Curved_Extend_S_EastBnd = IsSouthExtendEastBoundaryCurved();
    Curved_Extend_E_SouthBnd= IsEastExtendSouthBoundaryCurved();
    Curved_Extend_W_SouthBnd= IsWestExtendSouthBoundaryCurved();
    Curved_Extend_E_NorthBnd= IsEastExtendNorthBoundaryCurved();
    Curved_Extend_W_NorthBnd= IsWestExtendNorthBoundaryCurved();  

    // Determine the geometric properties of the ghost cells (e.g. centroid, area, geometric moments)
    // closed to curved boundaries
    Update_ExtensionSplineInfos();

    // Check the North boundary
    if (CurvedNorthBnd){
      // Update the ghost cells
      for(i=ICl+1; i<=ICu-1; ++i){
	Cell[i][JCu+1].A = area_GhostCell_CurvedBoundaries(i,NORTH_SPLINE);
	Cell[i][JCu+1].Xc = centroid_GhostCell_CurvedBoundaries(i,NORTH_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(i,JCu+1,NORTH_SPLINE);
      }

      // Update the first ghost cell on the North boundary
      if (Curved_Extend_N_WestBnd){
	// Both splines are curved
	Cell[ICl][JCu+1].A = area_GhostCell_CurvedBoundaries(ICl,CORNER_NORTH_EXTEND_N_WEST_SPLINES);
	Cell[ICl][JCu+1].Xc = centroid_GhostCell_CurvedBoundaries(ICl,CORNER_NORTH_EXTEND_N_WEST_SPLINES);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICl,JCu+1,CORNER_NORTH_EXTEND_N_WEST_SPLINES);
      } else {
	// Only North boundary is curved
	Cell[ICl][JCu+1].A = area_GhostCell_CurvedBoundaries(ICl,NORTH_SPLINE);
	Cell[ICl][JCu+1].Xc = centroid_GhostCell_CurvedBoundaries(ICl,NORTH_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICl,JCu+1,NORTH_SPLINE);
      }

      // Update the last ghost cell on the North boundary
      if (Curved_Extend_N_EastBnd){
	// Both splines are curved
	Cell[ICu][JCu+1].A = area_GhostCell_CurvedBoundaries(ICu,CORNER_EXTEND_N_EAST_NORTH_SPLINES);
	Cell[ICu][JCu+1].Xc = centroid_GhostCell_CurvedBoundaries(ICu,CORNER_EXTEND_N_EAST_NORTH_SPLINES);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICu,JCu+1,CORNER_EXTEND_N_EAST_NORTH_SPLINES);
      } else {
	// Only North boundary is curved
	Cell[ICu][JCu+1].A = area_GhostCell_CurvedBoundaries(ICu,NORTH_SPLINE);
	Cell[ICu][JCu+1].Xc = centroid_GhostCell_CurvedBoundaries(ICu,NORTH_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICu,JCu+1,NORTH_SPLINE);
      }
    } // endif (North Boundary)

    // Check the South boundary
    if (CurvedSouthBnd){
      // Update the ghost cells
      for(i=ICl+1; i<=ICu-1; ++i){
	Cell[i][JCl-1].A = area_GhostCell_CurvedBoundaries(i,SOUTH_SPLINE);
	Cell[i][JCl-1].Xc = centroid_GhostCell_CurvedBoundaries(i,SOUTH_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(i,JCl-1,SOUTH_SPLINE);
      }

      // Update the first ghost cell on the South boundary
      if (Curved_Extend_S_WestBnd){
	// Both splines are curved
	Cell[ICl][JCl-1].A = area_GhostCell_CurvedBoundaries(ICl,CORNER_EXTEND_S_WEST_SOUTH_SPLINES);
	Cell[ICl][JCl-1].Xc = centroid_GhostCell_CurvedBoundaries(ICl,CORNER_EXTEND_S_WEST_SOUTH_SPLINES);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICl,JCl-1,CORNER_EXTEND_S_WEST_SOUTH_SPLINES);
      } else {
	// Only South spline is curved
	Cell[ICl][JCl-1].A = area_GhostCell_CurvedBoundaries(ICl,SOUTH_SPLINE);
	Cell[ICl][JCl-1].Xc = centroid_GhostCell_CurvedBoundaries(ICl,SOUTH_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICl,JCl-1,SOUTH_SPLINE);
      }

      // Update the last ghost cell on the South boundary
      if (Curved_Extend_S_EastBnd){
	// Both splines are curved
	Cell[ICu][JCl-1].A = area_GhostCell_CurvedBoundaries(ICu,CORNER_SOUTH_EXTEND_S_EAST_SPLINES);
	Cell[ICu][JCl-1].Xc = centroid_GhostCell_CurvedBoundaries(ICu,CORNER_SOUTH_EXTEND_S_EAST_SPLINES);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICu,JCl-1,CORNER_SOUTH_EXTEND_S_EAST_SPLINES);
      } else {
	// Only South spline is curved
	Cell[ICu][JCl-1].A = area_GhostCell_CurvedBoundaries(ICu,SOUTH_SPLINE);
	Cell[ICu][JCl-1].Xc = centroid_GhostCell_CurvedBoundaries(ICu,SOUTH_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICu,JCl-1,SOUTH_SPLINE);
      }
    } // endif (South Boundary)

    // Check the East boundary
    if (CurvedEastBnd){
      // Update the ghost cells
      for(j=JCl+1; j<=JCu-1; ++j){
	Cell[ICu+1][j].A = area_GhostCell_CurvedBoundaries(j,EAST_SPLINE);
	Cell[ICu+1][j].Xc = centroid_GhostCell_CurvedBoundaries(j,EAST_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICu+1,j,EAST_SPLINE);
      }

      // Update the first ghost cell on the East boundary
      if (Curved_Extend_E_SouthBnd){
	// Both splines are curved
	Cell[ICu+1][JCl].A = area_GhostCell_CurvedBoundaries(JCl,CORNER_EXTEND_E_SOUTH_EAST_SPLINES);
	Cell[ICu+1][JCl].Xc = centroid_GhostCell_CurvedBoundaries(JCl,CORNER_EXTEND_E_SOUTH_EAST_SPLINES);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICu+1,JCl,CORNER_EXTEND_E_SOUTH_EAST_SPLINES);
      } else {
	// Only East spline is curved
	Cell[ICu+1][JCl].A = area_GhostCell_CurvedBoundaries(JCl,EAST_SPLINE);
	Cell[ICu+1][JCl].Xc = centroid_GhostCell_CurvedBoundaries(JCl,EAST_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICu+1,JCl,EAST_SPLINE);
      }

      // Update the last ghost cell on the East boundary
      if (Curved_Extend_E_NorthBnd){
	// Both splines are curved
	Cell[ICu+1][JCu].A = area_GhostCell_CurvedBoundaries(JCu,CORNER_EAST_EXTEND_E_NORTH_SPLINES);
	Cell[ICu+1][JCu].Xc = centroid_GhostCell_CurvedBoundaries(JCu,CORNER_EAST_EXTEND_E_NORTH_SPLINES);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICu+1,JCu,CORNER_EAST_EXTEND_E_NORTH_SPLINES);
      }	else {
	// Only East spline is curved
	Cell[ICu+1][JCu].A = area_GhostCell_CurvedBoundaries(JCu,EAST_SPLINE);
	Cell[ICu+1][JCu].Xc = centroid_GhostCell_CurvedBoundaries(JCu,EAST_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICu+1,JCu,EAST_SPLINE);
      }
    } // endif (East Boundary)
    
    // Check the West boundary
    if (CurvedWestBnd){
      // Update the ghost cells
      for(j=JCl+1; j<=JCu-1; ++j){
	Cell[ICl-1][j].A = area_GhostCell_CurvedBoundaries(j,WEST_SPLINE);
	Cell[ICl-1][j].Xc = centroid_GhostCell_CurvedBoundaries(j,WEST_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICl-1,j,WEST_SPLINE);
      }

      // Update the first ghost cell on the West boundary
      if (Curved_Extend_W_SouthBnd){
	// Both splines are curved
	Cell[ICl-1][JCl].A = area_GhostCell_CurvedBoundaries(JCl,CORNER_WEST_EXTEND_W_SOUTH_SPLINES);
	Cell[ICl-1][JCl].Xc = centroid_GhostCell_CurvedBoundaries(JCl,CORNER_WEST_EXTEND_W_SOUTH_SPLINES);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICl-1,JCl,CORNER_WEST_EXTEND_W_SOUTH_SPLINES);
      } else {
	// Only West spline is curved
	Cell[ICl-1][JCl].A = area_GhostCell_CurvedBoundaries(JCl,WEST_SPLINE);
	Cell[ICl-1][JCl].Xc = centroid_GhostCell_CurvedBoundaries(JCl,WEST_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICl-1,JCl,WEST_SPLINE);
      }

      // Update the last ghost cell on the West boundary
      if (Curved_Extend_W_NorthBnd){
	// Both splines are curved
	Cell[ICl-1][JCu].A = area_GhostCell_CurvedBoundaries(JCu,CORNER_EXTEND_W_NORTH_WEST_SPLINES);
	Cell[ICl-1][JCu].Xc = centroid_GhostCell_CurvedBoundaries(JCu,CORNER_EXTEND_W_NORTH_WEST_SPLINES);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICl-1,JCu,CORNER_EXTEND_W_NORTH_WEST_SPLINES);
      } else {
	// Only West spline is curved
	Cell[ICl-1][JCu].A = area_GhostCell_CurvedBoundaries(JCu,WEST_SPLINE);
	Cell[ICl-1][JCu].Xc = centroid_GhostCell_CurvedBoundaries(JCu,WEST_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICl-1,JCu,WEST_SPLINE);
      }
    } // endif (West Boundary)

    // ===  Update ghost cells influenced by the presence of curved extension splines ====
    Update_GhostCells_Near_CurvedExtensionSplines(CurvedNorthBnd,
						  CurvedSouthBnd,
						  CurvedEastBnd,
						  CurvedWestBnd,
						  Curved_Extend_N_WestBnd, 
						  Curved_Extend_S_WestBnd, 
						  Curved_Extend_N_EastBnd, 
						  Curved_Extend_S_EastBnd, 
						  Curved_Extend_E_SouthBnd,
						  Curved_Extend_W_SouthBnd,
						  Curved_Extend_E_NorthBnd,
						  Curved_Extend_W_NorthBnd);


    // Confirm the update
    Confirm_Ghost_Cells_Update();
  }
}

/*!
 * Updates the cell information for the quadrilateral mesh block 
 * corner ghost cells (i.e. no interior cells, no ghost cells other
 * than the corner ones and no extension or main boundary splines).
 */
void Grid2D_Quad_Block_HO::Update_Corner_Ghost_Cells(void) {

  int i,j;

  bool Curved_Extend_N_WestBnd, Curved_Extend_S_WestBnd, Curved_Extend_N_EastBnd, Curved_Extend_S_EastBnd;
  bool Curved_Extend_E_SouthBnd, Curved_Extend_W_SouthBnd, Curved_Extend_E_NorthBnd, Curved_Extend_W_NorthBnd;
 

  // Update cell information with straight boundaries (i.e. every edge of the cell is a line segment)
  for ( j = JCl-Nghost ; j <= JCl-1 ; ++j) {
    // Update the South-West corner
    for (i = ICl-Nghost; i <= ICl-1; ++i){
      Update_Cell(i,j);
    } // endif (i)

    // Update the South-East corner
    for (i = ICu+1; i <= ICu+Nghost; ++i){
      Update_Cell(i,j);
    } // endif (i)
  } // endif(j)

  for ( j = JCu+1 ; j <= JCu+Nghost ; ++j) {
    // Update the North-West corner
    for (i = ICl-Nghost; i <= ICl-1; ++i){
      Update_Cell(i,j);
    } // endif (i)

    // Update the North-East corner
    for (i = ICu+1; i <= ICu+Nghost; ++i){
      Update_Cell(i,j);
    } // endif (i)
  } // endif(j)

  if ( !CheckExistenceOfCurvedBoundaries() ){
    // There is no need for curved boundary representation
    
    // Confirm the update
    Confirm_Corner_Ghost_Cells_Update();
    return;
    
  } else {

    // Determine which boundaries are curved
    Curved_Extend_N_WestBnd = IsNorthExtendWestBoundaryCurved();
    Curved_Extend_S_WestBnd = IsSouthExtendWestBoundaryCurved();
    Curved_Extend_N_EastBnd = IsNorthExtendEastBoundaryCurved();
    Curved_Extend_S_EastBnd = IsSouthExtendEastBoundaryCurved();
    Curved_Extend_E_SouthBnd= IsEastExtendSouthBoundaryCurved();
    Curved_Extend_W_SouthBnd= IsWestExtendSouthBoundaryCurved();
    Curved_Extend_E_NorthBnd= IsEastExtendNorthBoundaryCurved();
    Curved_Extend_W_NorthBnd= IsWestExtendNorthBoundaryCurved();  

    // ===  Update ghost cells influenced by the presence of curved extension splines ====
    
    // Determine the geometric properties of the ghost cells (e.g. centroid, area, geometric moments)
    // closed to curved boundaries
    
    // === North-West corner ===
    // Check North extension of West spline
    if (Curved_Extend_N_WestBnd){ 
      for(j=JCu+2; j<=JCu+Nghost; ++j){
	// Update the right cells
	Cell[ICl-1][j].A = area_GhostCell_CurvedBoundaries(j,EXTEND_N_WEST_RIGHT_SPLINE);
	Cell[ICl-1][j].Xc = centroid_GhostCell_CurvedBoundaries(j,EXTEND_N_WEST_RIGHT_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICl-1,j,EXTEND_N_WEST_RIGHT_SPLINE);
      }

      // Check the Extend_N_West-Extend_W_North Corner
      if (Curved_Extend_W_NorthBnd){
	Cell[ICl-1][JCu+1].A = area_GhostCell_CurvedBoundaries(JCu+1,CORNER_EXTEND_N_WEST_EXTEND_W_NORTH_SPLINES);
	Cell[ICl-1][JCu+1].Xc = centroid_GhostCell_CurvedBoundaries(JCu+1,CORNER_EXTEND_N_WEST_EXTEND_W_NORTH_SPLINES);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICl-1,JCu+1,CORNER_EXTEND_N_WEST_EXTEND_W_NORTH_SPLINES);
      } else {
	Cell[ICl-1][JCu+1].A = area_GhostCell_CurvedBoundaries(JCu+1,EXTEND_N_WEST_RIGHT_SPLINE);
	Cell[ICl-1][JCu+1].Xc = centroid_GhostCell_CurvedBoundaries(JCu+1,EXTEND_N_WEST_RIGHT_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICl-1,JCu+1,EXTEND_N_WEST_RIGHT_SPLINE);
      }
    } // endif (Curved_Extend_N_WestBnd)

    // Check West extension of North spline
    if (Curved_Extend_W_NorthBnd){
      for(i=0; i<=ICl-2; ++i){
	// Update the left cells
	Cell[i][JCu+1].A = area_GhostCell_CurvedBoundaries(i,EXTEND_W_NORTH_LEFT_SPLINE);
	Cell[i][JCu+1].Xc = centroid_GhostCell_CurvedBoundaries(i,EXTEND_W_NORTH_LEFT_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(i,JCu+1,EXTEND_W_NORTH_LEFT_SPLINE);
      }
      
      // Re-check the Extend_N_West-Extend_W_North Corner
      if (!Curved_Extend_N_WestBnd){
	i = ICl - 1;
	Cell[i][JCu+1].A = area_GhostCell_CurvedBoundaries(i,EXTEND_W_NORTH_LEFT_SPLINE);
	Cell[i][JCu+1].Xc = centroid_GhostCell_CurvedBoundaries(i,EXTEND_W_NORTH_LEFT_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(i,JCu+1,EXTEND_W_NORTH_LEFT_SPLINE);
      }
    } // endif (Curved_Extend_W_NorthBnd)
    

    // === South-West corner ===
    // Check West extension of South spline
    if (Curved_Extend_W_SouthBnd){
      for(i=0; i<=ICl-2; ++i){
	// Update the right cells
	Cell[i][JCl-1].A = area_GhostCell_CurvedBoundaries(i,EXTEND_W_SOUTH_RIGHT_SPLINE);
	Cell[i][JCl-1].Xc = centroid_GhostCell_CurvedBoundaries(i,EXTEND_W_SOUTH_RIGHT_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(i,JCl-1,EXTEND_W_SOUTH_RIGHT_SPLINE);
      }
      
      // Check Extend_W_South-Extend_S_West Corner
      if (Curved_Extend_S_WestBnd){
	i = ICl-1;
	Cell[i][JCl-1].A = area_GhostCell_CurvedBoundaries(i,CORNER_EXTEND_W_SOUTH_EXTEND_S_WEST_SPLINES);
	Cell[i][JCl-1].Xc = centroid_GhostCell_CurvedBoundaries(i,CORNER_EXTEND_W_SOUTH_EXTEND_S_WEST_SPLINES);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(i,JCl-1,CORNER_EXTEND_W_SOUTH_EXTEND_S_WEST_SPLINES);
      } else {
	i = ICl-1;
	Cell[i][JCl-1].A = area_GhostCell_CurvedBoundaries(i,EXTEND_W_SOUTH_RIGHT_SPLINE);
	Cell[i][JCl-1].Xc = centroid_GhostCell_CurvedBoundaries(i,EXTEND_W_SOUTH_RIGHT_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(i,JCl-1,EXTEND_W_SOUTH_RIGHT_SPLINE);
      }
    } // endif (Curved_Extend_W_SouthBnd)

    // Check South extension of West spline
    if (Curved_Extend_S_WestBnd){
      for (j = 0; j<=JCl-2; ++j){
	// Update the left cells
	Cell[ICl-1][j].A = area_GhostCell_CurvedBoundaries(j,EXTEND_S_WEST_LEFT_SPLINE);
	Cell[ICl-1][j].Xc = centroid_GhostCell_CurvedBoundaries(j,EXTEND_S_WEST_LEFT_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICl-1,j,EXTEND_S_WEST_LEFT_SPLINE);
      }

      // Re-check Extend_W_South-Extend_S_West Corner
      if (!Curved_Extend_W_SouthBnd){
	j = JCl-1;
	Cell[ICl-1][j].A = area_GhostCell_CurvedBoundaries(j,EXTEND_S_WEST_LEFT_SPLINE);
	Cell[ICl-1][j].Xc = centroid_GhostCell_CurvedBoundaries(j,EXTEND_S_WEST_LEFT_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICl-1,j,EXTEND_S_WEST_LEFT_SPLINE);
      }
    } // endif (Curved_Extend_S_WestBnd)


    // === South-East corner ===
    // Check South extension of East spline
    if (Curved_Extend_S_EastBnd){
      for (j = 0; j<=JCl-2; ++j){
	// Update the right cells
	Cell[ICu+1][j].A = area_GhostCell_CurvedBoundaries(j,EXTEND_S_EAST_RIGHT_SPLINE);
	Cell[ICu+1][j].Xc = centroid_GhostCell_CurvedBoundaries(j,EXTEND_S_EAST_RIGHT_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICu+1,j,EXTEND_S_EAST_RIGHT_SPLINE);
      }

      // Check Extend_S_East-Extend_E_South Corner
      if (Curved_Extend_E_SouthBnd){
	j = JCl-1;
	Cell[ICu+1][j].A = area_GhostCell_CurvedBoundaries(j,CORNER_EXTEND_S_EAST_EXTEND_E_SOUTH_SPLINES);
	Cell[ICu+1][j].Xc = centroid_GhostCell_CurvedBoundaries(j,CORNER_EXTEND_S_EAST_EXTEND_E_SOUTH_SPLINES);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICu+1,j,CORNER_EXTEND_S_EAST_EXTEND_E_SOUTH_SPLINES);	
      } else {
	j = JCl-1;
	Cell[ICu+1][j].A = area_GhostCell_CurvedBoundaries(j,EXTEND_S_EAST_RIGHT_SPLINE);
	Cell[ICu+1][j].Xc = centroid_GhostCell_CurvedBoundaries(j,EXTEND_S_EAST_RIGHT_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICu+1,j,EXTEND_S_EAST_RIGHT_SPLINE);
      }
    } // endif (Curved_Exted_S_EastBnd)

    // Check East extension of South spline
    if (Curved_Extend_E_SouthBnd){
      for (i = ICu+2; i<=ICu+Nghost; ++i){
	// Update the left cells
	Cell[i][JCl-1].A = area_GhostCell_CurvedBoundaries(i,EXTEND_E_SOUTH_LEFT_SPLINE);
	Cell[i][JCl-1].Xc = centroid_GhostCell_CurvedBoundaries(i,EXTEND_E_SOUTH_LEFT_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(i,JCl-1,EXTEND_E_SOUTH_LEFT_SPLINE);
      }

      // Re-check Extend_S_East-Extend_E_South Corner
      if (!Curved_Extend_S_EastBnd){
	i = ICu+1;
	Cell[i][JCl-1].A = area_GhostCell_CurvedBoundaries(i,EXTEND_E_SOUTH_LEFT_SPLINE);
	Cell[i][JCl-1].Xc = centroid_GhostCell_CurvedBoundaries(i,EXTEND_E_SOUTH_LEFT_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(i,JCl-1,EXTEND_E_SOUTH_LEFT_SPLINE);
      }
    } // endif (Curved_Extend_E_SouthBnd)


    // === North-East corner ===
    // Check North extension of East spline
    if (Curved_Extend_N_EastBnd){
      for(j=JCu+2; j<=JCu+Nghost; ++j){
	// Update the left cells
	Cell[ICu+1][j].A = area_GhostCell_CurvedBoundaries(j,EXTEND_N_EAST_LEFT_SPLINE);
	Cell[ICu+1][j].Xc = centroid_GhostCell_CurvedBoundaries(j,EXTEND_N_EAST_LEFT_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICu+1,j,EXTEND_N_EAST_LEFT_SPLINE);
      }
      
      // Check the Extend_N_East-Extend_E_North Corner
      if (Curved_Extend_E_NorthBnd){
	j = JCu+1;
	Cell[ICu+1][j].A = area_GhostCell_CurvedBoundaries(j,CORNER_EXTEND_E_NORTH_EXTEND_N_EAST_SPLINES);
	Cell[ICu+1][j].Xc = centroid_GhostCell_CurvedBoundaries(j,CORNER_EXTEND_E_NORTH_EXTEND_N_EAST_SPLINES);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICu+1,j,CORNER_EXTEND_E_NORTH_EXTEND_N_EAST_SPLINES);	
      } else {
	j = JCu+1;
	Cell[ICu+1][j].A = area_GhostCell_CurvedBoundaries(j,EXTEND_N_EAST_LEFT_SPLINE);
	Cell[ICu+1][j].Xc = centroid_GhostCell_CurvedBoundaries(j,EXTEND_N_EAST_LEFT_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICu+1,j,EXTEND_N_EAST_LEFT_SPLINE);
      }
    } // endif (Curved_Extend_N_EastBnd)

    // Check East extension of North spline
    if (Curved_Extend_E_NorthBnd){
      for (i = ICu+2; i<=ICu+Nghost; ++i){
	// Update the right cells
	Cell[i][JCu+1].A = area_GhostCell_CurvedBoundaries(i,EXTEND_E_NORTH_RIGHT_SPLINE);
	Cell[i][JCu+1].Xc = centroid_GhostCell_CurvedBoundaries(i,EXTEND_E_NORTH_RIGHT_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(i,JCu+1,EXTEND_E_NORTH_RIGHT_SPLINE);
      }

      // Re-check Extend_E_North-Extend_N_East Corner
      if (!Curved_Extend_N_EastBnd){
	i = ICu+1;
	Cell[i][JCu+1].A = area_GhostCell_CurvedBoundaries(i,EXTEND_E_NORTH_RIGHT_SPLINE);
	Cell[i][JCu+1].Xc = centroid_GhostCell_CurvedBoundaries(i,EXTEND_E_NORTH_RIGHT_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(i,JCu+1,EXTEND_E_NORTH_RIGHT_SPLINE);	
      }
    } // endif (Curved_Extend_E_NorthBnd) 


    // Confirm the update
    Confirm_Corner_Ghost_Cells_Update();
  }
}

/*!
 * Updates the cell information for the cells influenced by the 
 * presence of the extension splines.
 * This routine is intended to avoid code duplication and not as a stand alone one.
 */
void Grid2D_Quad_Block_HO::Update_GhostCells_Near_CurvedExtensionSplines(const bool &CurvedNorthBnd,
									 const bool &CurvedSouthBnd,
									 const bool &CurvedEastBnd,
									 const bool &CurvedWestBnd,
									 const bool &Curved_Extend_N_WestBnd,
									 const bool &Curved_Extend_S_WestBnd,
									 const bool &Curved_Extend_N_EastBnd,
									 const bool &Curved_Extend_S_EastBnd,
									 const bool &Curved_Extend_E_SouthBnd,
									 const bool &Curved_Extend_W_SouthBnd,
									 const bool &Curved_Extend_E_NorthBnd,
									 const bool &Curved_Extend_W_NorthBnd){
  int i,j;

  // Check North extension of West spline
  if (Curved_Extend_N_WestBnd){ 
    for(j=JCu+2; j<=JCu+Nghost; ++j){
      // Update the right cells
      Cell[ICl-1][j].A = area_GhostCell_CurvedBoundaries(j,EXTEND_N_WEST_RIGHT_SPLINE);
      Cell[ICl-1][j].Xc = centroid_GhostCell_CurvedBoundaries(j,EXTEND_N_WEST_RIGHT_SPLINE);
      ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICl-1,j,EXTEND_N_WEST_RIGHT_SPLINE);
	
      // Update the left cells
      Cell[ICl][j].A = area_GhostCell_CurvedBoundaries(j,EXTEND_N_WEST_LEFT_SPLINE);
      Cell[ICl][j].Xc = centroid_GhostCell_CurvedBoundaries(j,EXTEND_N_WEST_LEFT_SPLINE);
      ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICl,j,EXTEND_N_WEST_LEFT_SPLINE);
    }

    // Re-check the Extend_N_West-North Corner
    if (!CurvedNorthBnd){
      Cell[ICl][JCu+1].A = area_GhostCell_CurvedBoundaries(JCu+1,EXTEND_N_WEST_LEFT_SPLINE);
      Cell[ICl][JCu+1].Xc = centroid_GhostCell_CurvedBoundaries(JCu+1,EXTEND_N_WEST_LEFT_SPLINE);
      ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICl,JCu+1,EXTEND_N_WEST_LEFT_SPLINE);
    }

    // Check the Extend_N_West-Extend_W_North Corner
    if (Curved_Extend_W_NorthBnd){
      Cell[ICl-1][JCu+1].A = area_GhostCell_CurvedBoundaries(JCu+1,CORNER_EXTEND_N_WEST_EXTEND_W_NORTH_SPLINES);
      Cell[ICl-1][JCu+1].Xc = centroid_GhostCell_CurvedBoundaries(JCu+1,CORNER_EXTEND_N_WEST_EXTEND_W_NORTH_SPLINES);
      ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICl-1,JCu+1,CORNER_EXTEND_N_WEST_EXTEND_W_NORTH_SPLINES);
    } else {
      Cell[ICl-1][JCu+1].A = area_GhostCell_CurvedBoundaries(JCu+1,EXTEND_N_WEST_RIGHT_SPLINE);
      Cell[ICl-1][JCu+1].Xc = centroid_GhostCell_CurvedBoundaries(JCu+1,EXTEND_N_WEST_RIGHT_SPLINE);
      ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICl-1,JCu+1,EXTEND_N_WEST_RIGHT_SPLINE);
    }
  } // endif (Curved_Extend_N_WestBnd)


  // Check North extension of East spline
  if (Curved_Extend_N_EastBnd){
    for(j=JCu+2; j<=JCu+Nghost; ++j){
      // Update the right cells
      Cell[ICu][j].A = area_GhostCell_CurvedBoundaries(j,EXTEND_N_EAST_RIGHT_SPLINE);
      Cell[ICu][j].Xc = centroid_GhostCell_CurvedBoundaries(j,EXTEND_N_EAST_RIGHT_SPLINE);
      ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICu,j,EXTEND_N_EAST_RIGHT_SPLINE);
	
      // Update the left cells
      Cell[ICu+1][j].A = area_GhostCell_CurvedBoundaries(j,EXTEND_N_EAST_LEFT_SPLINE);
      Cell[ICu+1][j].Xc = centroid_GhostCell_CurvedBoundaries(j,EXTEND_N_EAST_LEFT_SPLINE);
      ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICu+1,j,EXTEND_N_EAST_LEFT_SPLINE);
    }
      
    // Re-check the Extend_N_East-North Corner
    if (!CurvedNorthBnd){
      Cell[ICu][JCu+1].A = area_GhostCell_CurvedBoundaries(JCu+1,EXTEND_N_EAST_RIGHT_SPLINE);
      Cell[ICu][JCu+1].Xc = centroid_GhostCell_CurvedBoundaries(JCu+1,EXTEND_N_EAST_RIGHT_SPLINE);
      ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICu,JCu+1,EXTEND_N_EAST_RIGHT_SPLINE);
    }

    // Check the Extend_N_East-Extend_E_North Corner
    if (Curved_Extend_E_NorthBnd){
      j = JCu+1;
      Cell[ICu+1][j].A = area_GhostCell_CurvedBoundaries(j,CORNER_EXTEND_E_NORTH_EXTEND_N_EAST_SPLINES);
      Cell[ICu+1][j].Xc = centroid_GhostCell_CurvedBoundaries(j,CORNER_EXTEND_E_NORTH_EXTEND_N_EAST_SPLINES);
      ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICu+1,j,CORNER_EXTEND_E_NORTH_EXTEND_N_EAST_SPLINES);	
    } else {
      j = JCu+1;
      Cell[ICu+1][j].A = area_GhostCell_CurvedBoundaries(j,EXTEND_N_EAST_LEFT_SPLINE);
      Cell[ICu+1][j].Xc = centroid_GhostCell_CurvedBoundaries(j,EXTEND_N_EAST_LEFT_SPLINE);
      ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICu+1,j,EXTEND_N_EAST_LEFT_SPLINE);
    }
  } // endif (Curved_Extend_N_EastBnd)


  // Check West extension of North spline
  if (Curved_Extend_W_NorthBnd){
    for(i=0; i<=ICl-2; ++i){
      // Update the right cells
      Cell[i][JCu].A = area_GhostCell_CurvedBoundaries(i,EXTEND_W_NORTH_RIGHT_SPLINE);
      Cell[i][JCu].Xc = centroid_GhostCell_CurvedBoundaries(i,EXTEND_W_NORTH_RIGHT_SPLINE);
      ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(i,JCu,EXTEND_W_NORTH_RIGHT_SPLINE);
	
      // Update the left cells
      Cell[i][JCu+1].A = area_GhostCell_CurvedBoundaries(i,EXTEND_W_NORTH_LEFT_SPLINE);
      Cell[i][JCu+1].Xc = centroid_GhostCell_CurvedBoundaries(i,EXTEND_W_NORTH_LEFT_SPLINE);
      ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(i,JCu+1,EXTEND_W_NORTH_LEFT_SPLINE);
    }
      
    // Re-check the Extend_N_West-Extend_W_North Corner
    if (!Curved_Extend_N_WestBnd){
      i = ICl - 1;
      Cell[i][JCu+1].A = area_GhostCell_CurvedBoundaries(i,EXTEND_W_NORTH_LEFT_SPLINE);
      Cell[i][JCu+1].Xc = centroid_GhostCell_CurvedBoundaries(i,EXTEND_W_NORTH_LEFT_SPLINE);
      ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(i,JCu+1,EXTEND_W_NORTH_LEFT_SPLINE);
    }

    // Re-check the Extend_W_North-West Corner
    if (!CurvedWestBnd){
      i = ICl - 1;
      Cell[i][JCu].A = area_GhostCell_CurvedBoundaries(i,EXTEND_W_NORTH_RIGHT_SPLINE);
      Cell[i][JCu].Xc = centroid_GhostCell_CurvedBoundaries(i,EXTEND_W_NORTH_RIGHT_SPLINE);
      ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(i,JCu,EXTEND_W_NORTH_RIGHT_SPLINE);
    }
  } // endif (Curved_Extend_W_NorthBnd)


  // Check West extension of South spline
  if (Curved_Extend_W_SouthBnd){
    for(i=0; i<=ICl-2; ++i){
      // Update the right cells
      Cell[i][JCl-1].A = area_GhostCell_CurvedBoundaries(i,EXTEND_W_SOUTH_RIGHT_SPLINE);
      Cell[i][JCl-1].Xc = centroid_GhostCell_CurvedBoundaries(i,EXTEND_W_SOUTH_RIGHT_SPLINE);
      ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(i,JCl-1,EXTEND_W_SOUTH_RIGHT_SPLINE);
	
      // Update the left cells
      Cell[i][JCl].A = area_GhostCell_CurvedBoundaries(i,EXTEND_W_SOUTH_LEFT_SPLINE);
      Cell[i][JCl].Xc = centroid_GhostCell_CurvedBoundaries(i,EXTEND_W_SOUTH_LEFT_SPLINE);
      ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(i,JCl,EXTEND_W_SOUTH_LEFT_SPLINE);
    }
      
    // Re-check West-Extend_W_South Corner
    if (!CurvedWestBnd){
      i = ICl-1;
      Cell[i][JCl].A = area_GhostCell_CurvedBoundaries(i,EXTEND_W_SOUTH_LEFT_SPLINE);
      Cell[i][JCl].Xc = centroid_GhostCell_CurvedBoundaries(i,EXTEND_W_SOUTH_LEFT_SPLINE);
      ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(i,JCl,EXTEND_W_SOUTH_LEFT_SPLINE);
    }

    // Check Extend_W_South-Extend_S_West Corner
    if (Curved_Extend_S_WestBnd){
      i = ICl-1;
      Cell[i][JCl-1].A = area_GhostCell_CurvedBoundaries(i,CORNER_EXTEND_W_SOUTH_EXTEND_S_WEST_SPLINES);
      Cell[i][JCl-1].Xc = centroid_GhostCell_CurvedBoundaries(i,CORNER_EXTEND_W_SOUTH_EXTEND_S_WEST_SPLINES);
      ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(i,JCl-1,CORNER_EXTEND_W_SOUTH_EXTEND_S_WEST_SPLINES);
    } else {
      i = ICl-1;
      Cell[i][JCl-1].A = area_GhostCell_CurvedBoundaries(i,EXTEND_W_SOUTH_RIGHT_SPLINE);
      Cell[i][JCl-1].Xc = centroid_GhostCell_CurvedBoundaries(i,EXTEND_W_SOUTH_RIGHT_SPLINE);
      ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(i,JCl-1,EXTEND_W_SOUTH_RIGHT_SPLINE);
    }
  } // endif (Curved_Extend_W_SouthBnd)


  // Check South extension of West spline
  if (Curved_Extend_S_WestBnd){
    for (j = 0; j<=JCl-2; ++j){
      // Update the right cells
      Cell[ICl][j].A = area_GhostCell_CurvedBoundaries(j,EXTEND_S_WEST_RIGHT_SPLINE);
      Cell[ICl][j].Xc = centroid_GhostCell_CurvedBoundaries(j,EXTEND_S_WEST_RIGHT_SPLINE);
      ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICl,j,EXTEND_S_WEST_RIGHT_SPLINE);

      // Update the left cells
      Cell[ICl-1][j].A = area_GhostCell_CurvedBoundaries(j,EXTEND_S_WEST_LEFT_SPLINE);
      Cell[ICl-1][j].Xc = centroid_GhostCell_CurvedBoundaries(j,EXTEND_S_WEST_LEFT_SPLINE);
      ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICl-1,j,EXTEND_S_WEST_LEFT_SPLINE);
    }

    // Re-check Extend_W_South-Extend_S_West Corner
    if (!Curved_Extend_W_SouthBnd){
      j = JCl-1;
      Cell[ICl-1][j].A = area_GhostCell_CurvedBoundaries(j,EXTEND_S_WEST_LEFT_SPLINE);
      Cell[ICl-1][j].Xc = centroid_GhostCell_CurvedBoundaries(j,EXTEND_S_WEST_LEFT_SPLINE);
      ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICl-1,j,EXTEND_S_WEST_LEFT_SPLINE);
    }

    // Re-check Extend_S_West-South Corner
    if (!CurvedSouthBnd){
      j = JCl-1;
      Cell[ICl][j].A = area_GhostCell_CurvedBoundaries(j,EXTEND_S_WEST_RIGHT_SPLINE);
      Cell[ICl][j].Xc = centroid_GhostCell_CurvedBoundaries(j,EXTEND_S_WEST_RIGHT_SPLINE);
      ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICl,j,EXTEND_S_WEST_RIGHT_SPLINE);
    }
  } // endif (Curved_Extend_S_WestBnd)


  // Check South extension of East spline
  if (Curved_Extend_S_EastBnd){
    for (j = 0; j<=JCl-2; ++j){
      // Update the right cells
      Cell[ICu+1][j].A = area_GhostCell_CurvedBoundaries(j,EXTEND_S_EAST_RIGHT_SPLINE);
      Cell[ICu+1][j].Xc = centroid_GhostCell_CurvedBoundaries(j,EXTEND_S_EAST_RIGHT_SPLINE);
      ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICu+1,j,EXTEND_S_EAST_RIGHT_SPLINE);

      // Update the left cells
      Cell[ICu][j].A = area_GhostCell_CurvedBoundaries(j,EXTEND_S_EAST_LEFT_SPLINE);
      Cell[ICu][j].Xc = centroid_GhostCell_CurvedBoundaries(j,EXTEND_S_EAST_LEFT_SPLINE);
      ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICu,j,EXTEND_S_EAST_LEFT_SPLINE);
    }

    // Re-check South-Extend_S_East Corner
    if (!CurvedSouthBnd){
      j = JCl-1;
      Cell[ICu][j].A = area_GhostCell_CurvedBoundaries(j,EXTEND_S_EAST_LEFT_SPLINE);
      Cell[ICu][j].Xc = centroid_GhostCell_CurvedBoundaries(j,EXTEND_S_EAST_LEFT_SPLINE);
      ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICu,j,EXTEND_S_EAST_LEFT_SPLINE);
    }

    // Check Extend_S_East-Extend_E_South Corner
    if (Curved_Extend_E_SouthBnd){
      j = JCl-1;
      Cell[ICu+1][j].A = area_GhostCell_CurvedBoundaries(j,CORNER_EXTEND_S_EAST_EXTEND_E_SOUTH_SPLINES);
      Cell[ICu+1][j].Xc = centroid_GhostCell_CurvedBoundaries(j,CORNER_EXTEND_S_EAST_EXTEND_E_SOUTH_SPLINES);
      ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICu+1,j,CORNER_EXTEND_S_EAST_EXTEND_E_SOUTH_SPLINES);	
    } else {
      j = JCl-1;
      Cell[ICu+1][j].A = area_GhostCell_CurvedBoundaries(j,EXTEND_S_EAST_RIGHT_SPLINE);
      Cell[ICu+1][j].Xc = centroid_GhostCell_CurvedBoundaries(j,EXTEND_S_EAST_RIGHT_SPLINE);
      ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(ICu+1,j,EXTEND_S_EAST_RIGHT_SPLINE);
    }
  } // endif (Curved_Exted_S_EastBnd)


  // Check East extension of South spline
  if (Curved_Extend_E_SouthBnd){
    for (i = ICu+2; i<=ICu+Nghost; ++i){
      // Update the right cells
      Cell[i][JCl].A = area_GhostCell_CurvedBoundaries(i,EXTEND_E_SOUTH_RIGHT_SPLINE);
      Cell[i][JCl].Xc = centroid_GhostCell_CurvedBoundaries(i,EXTEND_E_SOUTH_RIGHT_SPLINE);
      ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(i,JCl,EXTEND_E_SOUTH_RIGHT_SPLINE);

      // Update the left cells
      Cell[i][JCl-1].A = area_GhostCell_CurvedBoundaries(i,EXTEND_E_SOUTH_LEFT_SPLINE);
	Cell[i][JCl-1].Xc = centroid_GhostCell_CurvedBoundaries(i,EXTEND_E_SOUTH_LEFT_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(i,JCl-1,EXTEND_E_SOUTH_LEFT_SPLINE);
      }

      // Re-check Extend_S_East-Extend_E_South Corner
      if (!Curved_Extend_S_EastBnd){
	i = ICu+1;
	Cell[i][JCl-1].A = area_GhostCell_CurvedBoundaries(i,EXTEND_E_SOUTH_LEFT_SPLINE);
	Cell[i][JCl-1].Xc = centroid_GhostCell_CurvedBoundaries(i,EXTEND_E_SOUTH_LEFT_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(i,JCl-1,EXTEND_E_SOUTH_LEFT_SPLINE);
      }

      // Re-check Extend_E_South-East Corner
      if (!CurvedEastBnd){
	i = ICu+1;
	Cell[i][JCl].A = area_GhostCell_CurvedBoundaries(i,EXTEND_E_SOUTH_RIGHT_SPLINE);
	Cell[i][JCl].Xc = centroid_GhostCell_CurvedBoundaries(i,EXTEND_E_SOUTH_RIGHT_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(i,JCl,EXTEND_E_SOUTH_RIGHT_SPLINE);
      }
    } // endif (Curved_Extend_E_SouthBnd)


    // Check East extension of North spline
    if (Curved_Extend_E_NorthBnd){
      for (i = ICu+2; i<=ICu+Nghost; ++i){
	// Update the right cells
	Cell[i][JCu+1].A = area_GhostCell_CurvedBoundaries(i,EXTEND_E_NORTH_RIGHT_SPLINE);
	Cell[i][JCu+1].Xc = centroid_GhostCell_CurvedBoundaries(i,EXTEND_E_NORTH_RIGHT_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(i,JCu+1,EXTEND_E_NORTH_RIGHT_SPLINE);

	// Update the left cells
	Cell[i][JCu].A = area_GhostCell_CurvedBoundaries(i,EXTEND_E_NORTH_LEFT_SPLINE);
	Cell[i][JCu].Xc = centroid_GhostCell_CurvedBoundaries(i,EXTEND_E_NORTH_LEFT_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(i,JCu,EXTEND_E_NORTH_LEFT_SPLINE);
      }

      // Re-check East-Extend_E_North Corner
      if (!CurvedEastBnd){
	i = ICu+1;
	Cell[i][JCu].A = area_GhostCell_CurvedBoundaries(i,EXTEND_E_NORTH_LEFT_SPLINE);
	Cell[i][JCu].Xc = centroid_GhostCell_CurvedBoundaries(i,EXTEND_E_NORTH_LEFT_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(i,JCu,EXTEND_E_NORTH_LEFT_SPLINE);
      }

      // Re-check Extend_E_North-Extend_N_East Corner
      if (!Curved_Extend_N_EastBnd){
	i = ICu+1;
	Cell[i][JCu+1].A = area_GhostCell_CurvedBoundaries(i,EXTEND_E_NORTH_RIGHT_SPLINE);
	Cell[i][JCu+1].Xc = centroid_GhostCell_CurvedBoundaries(i,EXTEND_E_NORTH_RIGHT_SPLINE);
	ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(i,JCu+1,EXTEND_E_NORTH_RIGHT_SPLINE);	
      }
    } // endif (Curved_Extend_E_NorthBnd) 
}

/*!
 * Update all grid properties without 
 * modifying the geometry.
 */
void Grid2D_Quad_Block_HO::Update_Grid_Properties(const int &HighestRecOrder){

  // Allocate new memory for the geometric moment container, if necessary
  allocate(NCi-2*Nghost, NCj-2*Nghost, Nghost, HighestRecOrder);

  /* Require update of the whole mesh */
  Schedule_Interior_Mesh_Update();
  Schedule_Ghost_Cells_Update();

  // Update properties for the whole mesh
  Update_Cells();

}

/*!
 * Check the validity of the interior cells of the 
 * quadrilateral mesh block.
 * Returns a non-zero result if the interior mesh
 * is not valid. 
 */
int Grid2D_Quad_Block_HO::Check_Quad_Block(void) {

  int i, j;
  int QuadType;

  // check if user required the mesh validity check
  if (HO_Grid2D_Execution_Mode::CHECK_FOR_INCORRECT_QUADRILATERALS){

    for ( j = JCl ; j <= JCu ; ++j) {
      for ( i = ICl ; i <= ICu ; ++i) {
	// Determine the type of the quadrilateral cell
	QuadType = Find_Quadrilateral_Type(Node[i  ][j  ].X,
					   Node[i+1][j  ].X,
					   Node[i+1][j+1].X,
					   Node[i  ][j+1].X);
      
	// Test the geometry of each cell for crossed or degenerated quadrilateral
	if ( QuadType == 4 || QuadType == 0 ){
	  cout << endl << "\nThe following quadrilateral has been identified as crossed or degenerated:\n";
	  Print_2(i,j);
	  Print_2(ICl,ICu);
	  Print_2(JCl,JCu);
	  Print_(Cell[i][j].A);
	  Print_(Cell[i][j].Xc);
	  Print_(Node[i][j].X);
	  Print_(Node[i+1][j].X);
	  Print_(Node[i+1][j+1].X);
	  Print_(Node[i][j+1].X);
	  cout.flush();
	  if (!HO_Grid2D_Execution_Mode::REPORT_INCORRECT_QUADRILATERALS_BUT_CONTINUE_EXECUTION){
	    return(1);
	  }
	} /* endif */
      } /* endfor */
    } /* endfor */

  }
  return(0);
}

/*!
 * Check the validity of all cells
 * (i.e. including ghost cells) of 
 * the quadrilateral mesh block.
 * Returns a non-zero result if the mesh
 * is not valid. 
 */
int Grid2D_Quad_Block_HO::Check_Quad_Block_Completely(void) {

  int i, j;
  int QuadType;

  // check if user required the mesh validity check
  if (HO_Grid2D_Execution_Mode::CHECK_FOR_INCORRECT_QUADRILATERALS){

    for ( j = JCl-Nghost ; j <= JCu+Nghost ; ++j) {
      for ( i = ICl-Nghost ; i <= ICu+Nghost ; ++i) {
	// Determine the type of the quadrilateral cell
	QuadType = Find_Quadrilateral_Type(Node[i  ][j  ].X,
					   Node[i+1][j  ].X,
					   Node[i+1][j+1].X,
					   Node[i  ][j+1].X);

	// Test the geometry of each cell for crossed or degenerated quadrilateral
	if ( QuadType == 4 || QuadType == 0 ){
	  cout << endl << "\nThe following quadrilateral has been identified as crossed or degenerated:\n";
	  Print_2(i,j);
	  Print_2(ICl,ICu);
	  Print_2(JCl,JCu);
	  Print_(Cell[i][j].A);
	  Print_(Cell[i][j].Xc);
	  Print_(Node[i][j].X);
	  Print_(Node[i+1][j].X);
	  Print_(Node[i+1][j+1].X);
	  Print_(Node[i][j+1].X);
	  cout.flush();
	  if (!HO_Grid2D_Execution_Mode::REPORT_INCORRECT_QUADRILATERALS_BUT_CONTINUE_EXECUTION){
	    return(1);
	  }
	} /* endif */
      } /* endfor */
    } /* endfor */

  }
  return(0);
}

/*!
 * Writes definition file information for 2D            
 * quadrilateral grid block to the specified output     
 * stream for retrieval and re-use purposes.            
 *                                                      
 */
void Grid2D_Quad_Block_HO::Write_Quad_Block_Definition(ostream &Out_File) {

    Out_File << setprecision(14);
    if (NNi == 0 || NNj == 0) {
       Out_File << NCi << " " 
	        << NCj << " "
		<< Nghost << " "
		<< HighestReconstructionOrder << "\n";
    } else {
       Out_File << NCi << " " 
	        << NCj << " "
		<< Nghost << " "
		<< HighestReconstructionOrder << "\n" 
                << BndNorthSpline.np << "\n"
                << BndNorthSpline.type << "\n"
                << BndNorthSpline 
                << BndSouthSpline.np << "\n"
                << BndSouthSpline.type << "\n"
                << BndSouthSpline
                << BndEastSpline.np << "\n"
                << BndEastSpline.type << "\n"
                << BndEastSpline 
                << BndWestSpline.np << "\n"
                << BndWestSpline.type << "\n"
                << BndWestSpline 
                << ExtendWest_BndNorthSpline.np << "\n"
                << ExtendWest_BndNorthSpline.type << "\n"
                << ExtendWest_BndNorthSpline 
                << ExtendEast_BndNorthSpline.np << "\n"
                << ExtendEast_BndNorthSpline.type << "\n"
                << ExtendEast_BndNorthSpline 
                << ExtendWest_BndSouthSpline.np << "\n"
                << ExtendWest_BndSouthSpline.type << "\n"
                << ExtendWest_BndSouthSpline 
                << ExtendEast_BndSouthSpline.np << "\n"
                << ExtendEast_BndSouthSpline.type << "\n"
                << ExtendEast_BndSouthSpline
                << ExtendNorth_BndEastSpline.np << "\n"
                << ExtendNorth_BndEastSpline.type << "\n"
                << ExtendNorth_BndEastSpline 
                << ExtendSouth_BndEastSpline.np << "\n"
                << ExtendSouth_BndEastSpline.type << "\n"
                << ExtendSouth_BndEastSpline 
                << ExtendNorth_BndWestSpline.np << "\n"
                << ExtendNorth_BndWestSpline.type << "\n"
                << ExtendNorth_BndWestSpline 
                << ExtendSouth_BndWestSpline.np << "\n"
                << ExtendSouth_BndWestSpline.type << "\n"
                << ExtendSouth_BndWestSpline 
                << GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH << "\n";
       Out_File.setf(ios::scientific);
       Out_File << StretchI << " " 
                << BetaI << " " 
                << TauI << "\n"
                << StretchJ << " " 
                << BetaJ << " " 
                << TauJ << "\n";
       Out_File.unsetf(ios::scientific);
       Out_File << OrthogonalN << " "
                << OrthogonalS << " "
                << OrthogonalE << " "
                << OrthogonalW << "\n";
    } /* endif */
    Out_File << setprecision(6);

}


/*!
 * Reads definition file for 2D quadrilateral grid      
 * block from the specified output stream.              
 *                                                      
 */
void Grid2D_Quad_Block_HO::Read_Quad_Block_Definition(istream &In_File) {
  
  int i, j, ng, HighestRecOrder, k, kx, ky, kx_max, ky_max, npts, spline_type;
  int node_init_procedure;
  double S_i, S_j, 
    s_north, s_south, s_east, s_west,
    smax_north, smax_south, smax_east, smax_west, 
    s_i, s_j, smax_i, smax_j,
    w_north, w_south, w_east, w_west, w_total;
  double smax_S1_NS, smax_S2_NS, smax_S1_EW, smax_S2_EW,
    S1_i, S2_i, S1_j, S2_j, dS_i, dS_j;
  Vector2D x_north, x_south, x_east, x_west;
  Spline2D_HO S1_NS, S2_NS, S1_EW, S2_EW;

  In_File.setf(ios::skipws);
  In_File >> i >> j >> ng >> HighestRecOrder;
  In_File.unsetf(ios::skipws);

  if (i != 0 && j != 0 && ng != 0) {
    /* Allocate (re-allocate) memory for the cells and nodes 
       of the quadrilateral mesh block as required. */
    allocate(i, j, ng, HighestRecOrder);
    
    /* For each of the north, south, east, and west boundaries
       of this mesh block, read in the number of spline points, 
       allocate memory for the boundary splines, read in the 
       spline points, and finally calculate the spline 
       pathlengths. */

    In_File.setf(ios::skipws);
    In_File >> npts;
    In_File.unsetf(ios::skipws);
    BndNorthSpline.allocate(npts);
    In_File.setf(ios::skipws);
    In_File >> spline_type;
    In_File.unsetf(ios::skipws);
    BndNorthSpline.settype(spline_type);
    In_File >> BndNorthSpline;
    BndNorthSpline.pathlength();
    SminN = BndNorthSpline.sp[0];
    SmaxN = BndNorthSpline.sp[BndNorthSpline.np-1];

    In_File.setf(ios::skipws);
    In_File >> npts;
    In_File.unsetf(ios::skipws);
    BndSouthSpline.allocate(npts);
    In_File.setf(ios::skipws);
    In_File >> spline_type;
    In_File.unsetf(ios::skipws);
    BndSouthSpline.settype(spline_type);
    In_File >> BndSouthSpline;
    BndSouthSpline.pathlength();
    SminS = BndSouthSpline.sp[0];
    SmaxS = BndSouthSpline.sp[BndSouthSpline.np-1];

    In_File.setf(ios::skipws);
    In_File >> npts;
    In_File.unsetf(ios::skipws);
    BndEastSpline.allocate(npts);
    In_File.setf(ios::skipws);
    In_File >> spline_type;
    In_File.unsetf(ios::skipws);
    BndEastSpline.settype(spline_type);
    In_File >> BndEastSpline;
    BndEastSpline.pathlength();
    SminE = BndEastSpline.sp[0];
    SmaxE = BndEastSpline.sp[BndEastSpline.np-1];

    In_File.setf(ios::skipws);
    In_File >> npts;
    In_File.unsetf(ios::skipws);
    BndWestSpline.allocate(npts);
    In_File.setf(ios::skipws);
    In_File >> spline_type;
    In_File.unsetf(ios::skipws);
    BndWestSpline.settype(spline_type);
    In_File >> BndWestSpline;
    BndWestSpline.pathlength();
    SminW = BndWestSpline.sp[0];
    SmaxW = BndWestSpline.sp[BndWestSpline.np-1];

    /* Read the extensions of the boundary splines. */

    In_File.setf(ios::skipws);
    In_File >> npts;
    In_File.unsetf(ios::skipws);
    ExtendWest_BndNorthSpline.allocate(npts);
    In_File.setf(ios::skipws);
    In_File >> spline_type;
    In_File.unsetf(ios::skipws);
    ExtendWest_BndNorthSpline.settype(spline_type);
    In_File >> ExtendWest_BndNorthSpline;
    ExtendWest_BndNorthSpline.pathlength();

    In_File.setf(ios::skipws);
    In_File >> npts;
    In_File.unsetf(ios::skipws);
    ExtendEast_BndNorthSpline.allocate(npts);
    In_File.setf(ios::skipws);
    In_File >> spline_type;
    In_File.unsetf(ios::skipws);
    ExtendEast_BndNorthSpline.settype(spline_type);
    In_File >> ExtendEast_BndNorthSpline;
    ExtendEast_BndNorthSpline.pathlength();

    In_File.setf(ios::skipws);
    In_File >> npts;
    In_File.unsetf(ios::skipws);
    ExtendWest_BndSouthSpline.allocate(npts);
    In_File.setf(ios::skipws);
    In_File >> spline_type;
    In_File.unsetf(ios::skipws);
    ExtendWest_BndSouthSpline.settype(spline_type);
    In_File >> ExtendWest_BndSouthSpline;
    ExtendWest_BndSouthSpline.pathlength();

    In_File.setf(ios::skipws);
    In_File >> npts;
    In_File.unsetf(ios::skipws);
    ExtendEast_BndSouthSpline.allocate(npts);
    In_File.setf(ios::skipws);
    In_File >> spline_type;
    In_File.unsetf(ios::skipws);
    ExtendEast_BndSouthSpline.settype(spline_type);
    In_File >> ExtendEast_BndSouthSpline;
    ExtendEast_BndSouthSpline.pathlength();

    In_File.setf(ios::skipws);
    In_File >> npts;
    In_File.unsetf(ios::skipws);
    ExtendNorth_BndEastSpline.allocate(npts);
    In_File.setf(ios::skipws);
    In_File >> spline_type;
    In_File.unsetf(ios::skipws);
    ExtendNorth_BndEastSpline.settype(spline_type);
    In_File >> ExtendNorth_BndEastSpline;
    ExtendNorth_BndEastSpline.pathlength();

    In_File.setf(ios::skipws);
    In_File >> npts;
    In_File.unsetf(ios::skipws);
    ExtendSouth_BndEastSpline.allocate(npts);
    In_File.setf(ios::skipws);
    In_File >> spline_type;
    In_File.unsetf(ios::skipws);
    ExtendSouth_BndEastSpline.settype(spline_type);
    In_File >> ExtendSouth_BndEastSpline;
    ExtendSouth_BndEastSpline.pathlength();

    In_File.setf(ios::skipws);
    In_File >> npts;
    In_File.unsetf(ios::skipws);
    ExtendNorth_BndWestSpline.allocate(npts);
    In_File.setf(ios::skipws);
    In_File >> spline_type;
    In_File.unsetf(ios::skipws);
    ExtendNorth_BndWestSpline.settype(spline_type);
    In_File >> ExtendNorth_BndWestSpline;
    ExtendNorth_BndWestSpline.pathlength();

    In_File.setf(ios::skipws);
    In_File >> npts;
    In_File.unsetf(ios::skipws);
    ExtendSouth_BndWestSpline.allocate(npts);
    In_File.setf(ios::skipws);
    In_File >> spline_type;
    In_File.unsetf(ios::skipws);
    ExtendSouth_BndWestSpline.settype(spline_type);
    In_File >> ExtendSouth_BndWestSpline;
    ExtendSouth_BndWestSpline.pathlength();

    /* Read the node initialization procedure for this 
       quadrilateral grid block. */

    In_File.setf(ios::skipws);
    In_File >> node_init_procedure;
    In_File.unsetf(ios::skipws);
      
    /* Input the type of node distribution stretching 
       functions to be used for each of the coordinate
       directions. Also read in the related stretching
       parameters. */

    In_File.setf(ios::skipws);
    In_File >> StretchI >> BetaI >> TauI;
    In_File >> StretchJ >> BetaJ >> TauJ;
    In_File.unsetf(ios::skipws);

    /* Input the grid orthogonality specifiers for
       the north, south, east, and west boundaries. */

    In_File.setf(ios::skipws);
    In_File >> OrthogonalN >> OrthogonalS
	    >> OrthogonalE >> OrthogonalW;
    In_File.unsetf(ios::skipws);

    /* Compute the interior nodes for the quadrilateral mesh block. */

    smax_north = BndNorthSpline.sp[BndNorthSpline.np-1]-
      BndNorthSpline.sp[0];
    smax_south = BndSouthSpline.sp[BndSouthSpline.np-1]-
      BndSouthSpline.sp[0];

    smax_east = BndEastSpline.sp[BndEastSpline.np-1]-
      BndEastSpline.sp[0];
    smax_west = BndWestSpline.sp[BndWestSpline.np-1]-
      BndWestSpline.sp[0];

    dS_i = ONE/(TEN*double(INu-INl));
    dS_j = ONE/(TEN*double(JNu-JNl));

    for ( j = JNl ; j <= JNu ; ++j) {
      for ( i = INl ; i <= INu ; ++i) {
	S_j = double(j-JNl)/double(JNu-JNl);
	S_j = StretchingFcn(S_j, BetaJ, TauJ, StretchJ);
	s_east = S_j * smax_east + BndEastSpline.sp[0];
	x_east = Spline(s_east, BndEastSpline);
	s_west = S_j * smax_west + BndWestSpline.sp[0];
	x_west = Spline(s_west, BndWestSpline);

	S_i = double(i-INl)/double(INu-INl);
	S_i = StretchingFcn(S_i, BetaI, TauI, StretchI);
	s_north = S_i * smax_north + BndNorthSpline.sp[0];
	x_north = Spline(s_north, BndNorthSpline);
	s_south = S_i * smax_south + BndSouthSpline.sp[0];
	x_south = Spline(s_south, BndSouthSpline);

	if (j == JNl) {
	  Node[i][j].X = x_south;
	} else if (j == JNu) {
	  Node[i][j].X = x_north;
	} else if (i == INl) {
	  Node[i][j].X = x_west;
	} else if (i == INu) {
	  Node[i][j].X = x_east;
	} else {
	  s_i = (ONE-S_j)*s_south + S_j*s_north;
	  s_j = (ONE-S_i)*s_west + S_i*s_east;

	  smax_i = (ONE-S_j)*smax_south + S_j*smax_north;
	  smax_j = (ONE-S_i)*smax_west + S_i*smax_east;

	  switch(node_init_procedure) {
	    //============================================================
	  case GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH :
	    w_south = (ONE-s_j/smax_j);
	    w_north = (s_j/smax_j);
	    w_total = w_south + w_north;
	    Node[i][j].X = (w_south*x_south + 
				 w_north*x_north)/w_total;
	    break;
	    //============================================================
	  case GRID2D_QUAD_BLOCK_INIT_PROCEDURE_EAST_WEST :
	    w_west =  (ONE-s_i/smax_i);
	    w_east =  (s_i/smax_i);
	    w_total = w_east + w_west;
	    Node[i][j].X = (w_west*x_west + 
				 w_east*x_east)/w_total;
	    break;
	    //============================================================
	  case GRID2D_QUAD_BLOCK_INIT_PROCEDURE_TRANS_FINITE_XY :
	    kx_max = 5;
	    dS_i = ONE/double(kx_max-1);
	    kx = int(floor(S_i/dS_i));
	    S1_i = max(ZERO, double(kx)*dS_i);
	    S2_i = min(ONE, S1_i + dS_i);

	    ky_max = 7;
	    dS_j = ONE/double(ky_max-1);
	    ky = int(floor(S_j/dS_j));
	    S1_j = max(ZERO, double(ky)*dS_j);
	    S2_j = min(ONE, S1_j + dS_j);

	    if (S1_i - TOLER > ZERO) {
	      S1_NS.allocate(2);
	      S1_NS.settype(SPLINE2D_LINEAR);

	      s_south = S1_i * smax_south + BndSouthSpline.sp[0];
	      S1_NS.Xp[0] = Spline(s_south, BndSouthSpline);
	      S1_NS.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	      s_north = S1_i * smax_north + BndNorthSpline.sp[0],
		S1_NS.Xp[1] = Spline(s_north, BndNorthSpline);
	      S1_NS.tp[1] = SPLINE2D_POINT_SHARP_CORNER;

	      S1_NS.pathlength();
	      smax_S1_NS = S1_NS.sp[S1_NS.np-1]-
		S1_NS.sp[0];
	    } else {
	      S1_NS = BndWestSpline;
	      smax_S1_NS = S1_NS.sp[S1_NS.np-1]-
		S1_NS.sp[0];
	    } /* endif */

	    if (S2_i + TOLER < ONE) {
	      S2_NS.allocate(2);
	      S2_NS.settype(SPLINE2D_LINEAR);

	      s_south = S2_i * smax_south + BndSouthSpline.sp[0];
	      S2_NS.Xp[0] = Spline(s_south, BndSouthSpline);
	      S2_NS.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	      s_north = S2_i * smax_north + BndNorthSpline.sp[0],
		S2_NS.Xp[1] = Spline(s_north, BndNorthSpline);
	      S2_NS.tp[1] = SPLINE2D_POINT_SHARP_CORNER;

	      S2_NS.pathlength();
	      smax_S2_NS = S2_NS.sp[S2_NS.np-1]-
		S2_NS.sp[0];
	    } else {
	      S2_NS = BndEastSpline;
	      smax_S2_NS = S2_NS.sp[S2_NS.np-1]-
		S2_NS.sp[0];
	    } /* endif */

	    if (S1_j - TOLER > ZERO) {
	      S1_EW.allocate(kx_max);
	      S1_EW.settype(SPLINE2D_LINEAR);

	      s_west = S1_j * smax_west + BndWestSpline.sp[0];
	      S1_EW.Xp[0] = Spline(s_west, BndWestSpline);
	      S1_EW.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	      s_east = S1_j * smax_east + BndEastSpline.sp[0];
	      S1_EW.Xp[kx_max-1] = Spline(s_east, BndEastSpline);
	      S1_EW.tp[kx_max-1] = SPLINE2D_POINT_SHARP_CORNER;

	      for ( k = 1 ; k <= kx_max-2 ; ++k) {                     
		s_north = double(k)*smax_north/double(kx_max-1) + BndNorthSpline.sp[0];
		x_north = Spline(s_north, BndNorthSpline);

		s_south = double(k)*smax_south/double(kx_max-1) + BndSouthSpline.sp[0];
		x_south = Spline(s_south, BndSouthSpline);

		w_south =  ONE-S1_j;
		w_north =  S1_j;
		w_total = w_north + w_south;
		S1_EW.Xp[k] = (w_south*x_south + w_north*x_north)/w_total;
		S1_EW.tp[k] = SPLINE2D_POINT_SHARP_CORNER;
	      } /* endfor */

	      S1_EW.pathlength();
	      smax_S1_EW = S1_EW.sp[S1_EW.np-1]-
		S1_EW.sp[0];
	    } else {
	      S1_EW = BndSouthSpline;
	      smax_S1_EW = S1_EW.sp[S1_EW.np-1]-
		S1_EW.sp[0];
	    } /* endif */

	    if (S2_j + TOLER < ONE) {
	      S2_EW.allocate(kx_max);
	      S2_EW.settype(SPLINE2D_LINEAR);

	      s_west = S2_j * smax_west + BndWestSpline.sp[0];
	      S2_EW.Xp[0] = Spline(s_west, BndWestSpline);
	      S2_EW.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	      s_east = S2_j * smax_east + BndEastSpline.sp[0];
	      S2_EW.Xp[kx_max-1] = Spline(s_east, BndEastSpline);
	      S2_EW.tp[kx_max-1] = SPLINE2D_POINT_SHARP_CORNER;

	      for ( k = 1 ; k <= kx_max-2 ; ++k) {                     
		s_north = double(k)*smax_north/double(kx_max-1) + BndNorthSpline.sp[0];
		x_north = Spline(s_north, BndNorthSpline);

		s_south = double(k)*smax_south/double(kx_max-1) + BndSouthSpline.sp[0];
		x_south = Spline(s_south, BndSouthSpline);

		w_south =  ONE-S2_j;
		w_north =  S2_j;
		w_total = w_north + w_south;
		S2_EW.Xp[k] = (w_south*x_south + w_north*x_north)/w_total;
		S2_EW.tp[k] = SPLINE2D_POINT_SHARP_CORNER;
	      } /* endfor */

	      S2_EW.pathlength();
	      smax_S2_EW = S2_EW.sp[S2_EW.np-1]-
		S2_EW.sp[0];
	    } else {
	      S2_EW = BndNorthSpline;
	      smax_S2_EW = S2_EW.sp[S2_EW.np-1]-
		S2_EW.sp[0];
	    } /* endif */

	    s_north = S_i * smax_S2_EW + S2_EW.sp[0];
	    x_north = Spline(s_north, S2_EW);
	    s_south = S_i * smax_S1_EW + S1_EW.sp[0];
	    x_south = Spline(s_south, S1_EW);

	    s_east = S_j * smax_S2_NS + S2_NS.sp[0];
	    x_east = Spline(s_east, S2_NS);
	    s_west = S_j * smax_S1_NS + S1_NS.sp[0];
	    x_west = Spline(s_west, S1_NS);
      
	    if (ky == 0 || ky == ky_max - 2) {
	      w_south =  ONE-(S_j-S1_j)/(S2_j-S1_j);
	      w_north =  (S_j-S1_j)/(S2_j-S1_j);
	      w_total = w_south + w_north;
	      Node[i][j].X = (w_south*x_south + 
				   w_north*x_north)/w_total;
	    } else {
	      w_west =  ONE-(S_i-S1_i)/(S2_i-S1_i);
	      w_east =  (S_i-S1_i)/(S2_i-S1_i);
	      w_total = w_east + w_west;
	      Node[i][j].X = (w_west*x_west + 
				   w_east*x_east)/w_total;
	    } /* endif */

	    S1_NS.deallocate(); S2_NS.deallocate();
	    S2_EW.deallocate(); S2_EW.deallocate();
	    break;
	    //============================================================
	  case GRID2D_QUAD_BLOCK_INIT_PROCEDURE_TRANS_FINITE_YX :
	    kx_max = 7;
	    dS_i = ONE/double(kx_max-1);
	    kx = int(floor(S_i/dS_i));
	    S1_i = max(ZERO, double(kx)*dS_i);
	    S2_i = min(ONE, S1_i + dS_i);

	    ky_max = 5;
	    dS_j = ONE/double(ky_max-1);
	    ky = int(floor(S_j/dS_j));
	    S1_j = max(ZERO, double(ky)*dS_j);
	    S2_j = min(ONE, S1_j + dS_j);

	    if (S1_i - TOLER > ZERO) {
	      S1_NS.allocate(ky_max);
	      S1_NS.settype(SPLINE2D_LINEAR);

	      s_south = S1_i * smax_south + BndSouthSpline.sp[0];
	      S1_NS.Xp[0] = Spline(s_south, BndSouthSpline);
	      S1_NS.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	      s_north = S1_i * smax_north + BndNorthSpline.sp[0],
		S1_NS.Xp[ky_max-1] = Spline(s_north, BndNorthSpline);
	      S1_NS.tp[ky_max-1] = SPLINE2D_POINT_SHARP_CORNER;

	      for ( k = 1 ; k <= ky_max-2 ; ++k) {                     
		s_east = double(k)*smax_east/double(ky_max-1) + BndEastSpline.sp[0];
		x_east = Spline(s_east, BndEastSpline);

		s_west = double(k)*smax_west/double(ky_max-1) + BndWestSpline.sp[0];
		x_west = Spline(s_west, BndWestSpline);

		w_west =  ONE-S1_i;
		w_east =  S1_i;
		w_total = w_east + w_west;
		S1_NS.Xp[k] = (w_west*x_west + w_east*x_east)/w_total;
		S1_NS.tp[k] = SPLINE2D_POINT_SHARP_CORNER;
	      } /* endfor */

	      S1_NS.pathlength();
	      smax_S1_NS = S1_NS.sp[S1_NS.np-1]-
		S1_NS.sp[0];
	    } else {
	      S1_NS = BndWestSpline;
	      smax_S1_NS = S1_NS.sp[S1_NS.np-1]-
		S1_NS.sp[0];
	    } /* endif */

	    if (S2_i + TOLER < ONE) {
	      S2_NS.allocate(ky_max);
	      S2_NS.settype(SPLINE2D_LINEAR);

	      s_south = S2_i * smax_south + BndSouthSpline.sp[0];
	      S2_NS.Xp[0] = Spline(s_south, BndSouthSpline);
	      S2_NS.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	      s_north = S2_i * smax_north + BndNorthSpline.sp[0],
		S2_NS.Xp[ky_max-1] = Spline(s_north, BndNorthSpline);
	      S2_NS.tp[ky_max-1] = SPLINE2D_POINT_SHARP_CORNER;

	      for ( k = 1 ; k <= ky_max-2 ; ++k) {
		s_east = double(k)*smax_east/double(ky_max-1) + BndEastSpline.sp[0];
		x_east = Spline(s_east, BndEastSpline);

		s_west = double(k)*smax_west/double(ky_max-1) + BndWestSpline.sp[0];
		x_west = Spline(s_west, BndWestSpline);

		w_west =  ONE-S2_i;
		w_east =  S2_i;
		w_total = w_east + w_west;
		S2_NS.Xp[k] = (w_west*x_west + w_east*x_east)/w_total;
		S2_NS.tp[k] = SPLINE2D_POINT_SHARP_CORNER;
	      } /* endfor */

	      S2_NS.pathlength();
	      smax_S2_NS = S2_NS.sp[S2_NS.np-1]-
		S2_NS.sp[0];
	    } else {
	      S2_NS = BndEastSpline;
	      smax_S2_NS = S2_NS.sp[S2_NS.np-1]-
		S2_NS.sp[0];
	    } /* endif */

	    if (S1_j - TOLER > ZERO) {
	      S1_EW.allocate(2);
	      S1_EW.settype(SPLINE2D_LINEAR);

	      s_west = S1_j * smax_west + BndWestSpline.sp[0];
	      S1_EW.Xp[0] = Spline(s_west, BndWestSpline);
	      S1_EW.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	      s_east = S1_j * smax_east + BndEastSpline.sp[0];
	      S1_EW.Xp[1] = Spline(s_east, BndEastSpline);
	      S1_EW.tp[1] = SPLINE2D_POINT_SHARP_CORNER;

	      S1_EW.pathlength();
	      smax_S1_EW = S1_EW.sp[S1_EW.np-1]-
		S1_EW.sp[0];
	    } else {
	      S1_EW = BndSouthSpline;
	      smax_S1_EW = S1_EW.sp[S1_EW.np-1]-
		S1_EW.sp[0];
	    } /* endif */

	    if (S2_j + TOLER < ONE) {
	      S2_EW.allocate(2);
	      S2_EW.settype(SPLINE2D_LINEAR);
   
	      s_west = S2_j * smax_west + BndWestSpline.sp[0];
	      S2_EW.Xp[0] = Spline(s_west, BndWestSpline);
	      S2_EW.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	      s_east = S2_j * smax_east + BndEastSpline.sp[0];
	      S2_EW.Xp[1] = Spline(s_east, BndEastSpline);
	      S2_EW.tp[1] = SPLINE2D_POINT_SHARP_CORNER;

	      S2_EW.pathlength();
	      smax_S2_EW = S2_EW.sp[S2_EW.np-1]-
		S2_EW.sp[0];
	    } else {
	      S2_EW = BndNorthSpline;
	      smax_S2_EW = S2_EW.sp[S2_EW.np-1]-
		S2_EW.sp[0];
	    } /* endif */

	    s_north = S_i * smax_S2_EW + S2_EW.sp[0];
	    x_north = Spline(s_north, S2_EW);
	    s_south = S_i * smax_S1_EW + S1_EW.sp[0];
	    x_south = Spline(s_south, S1_EW);

	    s_east = S_j * smax_S2_NS + S2_NS.sp[0];
	    x_east = Spline(s_east, S2_NS);
	    s_west = S_j * smax_S1_NS + S1_NS.sp[0];
	    x_west = Spline(s_west, S1_NS);
      
	    if (kx == 0 || kx == kx_max - 2) {
	      w_west =  ONE-(S_i-S1_i)/(S2_i-S1_i);
	      w_east =  (S_i-S1_i)/(S2_i-S1_i);
	      w_total = w_east + w_west;
	      Node[i][j].X = (w_west*x_west + 
				   w_east*x_east)/w_total;
	    } else {
	      w_south =  ONE-(S_j-S1_j)/(S2_j-S1_j);
	      w_north =  (S_j-S1_j)/(S2_j-S1_j);
	      w_total = w_south + w_north;
	      Node[i][j].X = (w_south*x_south + 
				   w_north*x_north)/w_total;
	    } /* endif */

	    S1_NS.deallocate(); S2_NS.deallocate();
	    S2_EW.deallocate(); S2_EW.deallocate();
	    break;
	    //============================================================
	  default:
	    w_south = (ONE-s_j/smax_j);
	    w_north = (s_j/smax_j);
	    w_total = w_south + w_north;
	    Node[i][j].X = (w_south*x_south + 
				 w_north*x_north)/w_total;
	    break;
	  } /* endswitch */

	} /* endif */

      } /* endfor */
    }/* endfor */

    /* Require update of the interior cells geometric properties. */
    Schedule_Interior_Mesh_Update();

    /* Set the boundary condition types at the quadrilateral 
       grid block boundaries. */

    Set_BCs();

    /* Compute the exterior nodes for the quadrilateral mesh block. */

    Update_Exterior_Nodes();

    /* Compute the cells for the quadrilateral mesh block. */

    Update_Cells();
  } /* endif */

}


/*!
 * Writes 2D quadrilateral grid block to the            
 * specified output stream for retrieval and re-use     
 * purposes.                                            
 */
void Grid2D_Quad_Block_HO::Write_Quad_Block(ostream &Out_File) {

  Out_File << setprecision(14) << *this << setprecision(6);

}

/*!
 * Reads 2D quadrilateral grid block from the           
 * specified input stream.                              
 *                                                      
 */
void Grid2D_Quad_Block_HO::Read_Quad_Block(istream &In_File) {

    In_File >> *this;

}

/*!
 * Translates or shifts the positions of the nodes of a 
 * quadrilateral grid block.                            
 * 
 * \note This subroutine DOESN'T update the geometric properties
 *       of the new grid cells!!!                                                     
 */
void Grid2D_Quad_Block_HO::Translate_Quad_Block_Without_Update(const Vector2D &V) {

    int i, j;

    for ( j = JNl-Nghost ; j <= JNu+Nghost; ++j ) {
       for ( i = INl-Nghost ; i <= INu+Nghost; ++i ) {
           Node[i][j].X += V;
       } /* endfor */
    } /* endfor */

    if (BndNorthSpline.np != 0 ) 
       BndNorthSpline.Translate_Spline(V);
    if (BndSouthSpline.np != 0 ) 
       BndSouthSpline.Translate_Spline(V);
    if (BndEastSpline.np != 0 ) 
       BndEastSpline.Translate_Spline(V);
    if (BndWestSpline.np != 0 )
       BndWestSpline.Translate_Spline(V);

    if (ExtendWest_BndNorthSpline.np != 0) ExtendWest_BndNorthSpline.Translate_Spline(V);
    if (ExtendEast_BndNorthSpline.np != 0) ExtendEast_BndNorthSpline.Translate_Spline(V);
    if (ExtendWest_BndSouthSpline.np != 0) ExtendWest_BndSouthSpline.Translate_Spline(V);
    if (ExtendEast_BndSouthSpline.np != 0) ExtendEast_BndSouthSpline.Translate_Spline(V);
    if (ExtendNorth_BndEastSpline.np != 0) ExtendNorth_BndEastSpline.Translate_Spline(V);
    if (ExtendSouth_BndEastSpline.np != 0) ExtendSouth_BndEastSpline.Translate_Spline(V);
    if (ExtendNorth_BndWestSpline.np != 0) ExtendNorth_BndWestSpline.Translate_Spline(V);
    if (ExtendSouth_BndWestSpline.np != 0) ExtendSouth_BndWestSpline.Translate_Spline(V);
    
 
    /* Require update of the whole mesh */
    Schedule_Interior_Mesh_Update();
    Schedule_Ghost_Cells_Update();
}


/*!
 * Scales the quadrilateral grid block.                 
 *                                                      
 * \note This subroutine DOESN'T update the geometric properties
 *       of the new grid cells!!!                                                     
 */
void Grid2D_Quad_Block_HO::Scale_Quad_Block_Without_Update(const double &Scaling_Factor) {

    int i, j;

    for ( j = JNl-Nghost ; j <= JNu+Nghost; ++j ) {
       for ( i = INl-Nghost ; i <= INu+Nghost; ++i ) {
           Node[i][j].X = Node[i][j].X*Scaling_Factor;
       } /* endfor */
    } /* endfor */

    if (BndNorthSpline.np != 0 ) 
       BndNorthSpline.Scale_Spline(Scaling_Factor);
    if (BndSouthSpline.np != 0 ) 
       BndSouthSpline.Scale_Spline(Scaling_Factor);
    if (BndEastSpline.np != 0 ) 
       BndEastSpline.Scale_Spline(Scaling_Factor);
    if (BndWestSpline.np != 0 )
       BndWestSpline.Scale_Spline(Scaling_Factor);

    if (ExtendWest_BndNorthSpline.np != 0) ExtendWest_BndNorthSpline.Scale_Spline(Scaling_Factor);
    if (ExtendEast_BndNorthSpline.np != 0) ExtendEast_BndNorthSpline.Scale_Spline(Scaling_Factor);
    if (ExtendWest_BndSouthSpline.np != 0) ExtendWest_BndSouthSpline.Scale_Spline(Scaling_Factor);
    if (ExtendEast_BndSouthSpline.np != 0) ExtendEast_BndSouthSpline.Scale_Spline(Scaling_Factor);
    if (ExtendNorth_BndEastSpline.np != 0) ExtendNorth_BndEastSpline.Scale_Spline(Scaling_Factor);
    if (ExtendSouth_BndEastSpline.np != 0) ExtendSouth_BndEastSpline.Scale_Spline(Scaling_Factor);
    if (ExtendNorth_BndWestSpline.np != 0) ExtendNorth_BndWestSpline.Scale_Spline(Scaling_Factor);
    if (ExtendSouth_BndWestSpline.np != 0) ExtendSouth_BndWestSpline.Scale_Spline(Scaling_Factor);

    SminN = SminN*Scaling_Factor;
    SmaxN = SmaxN*Scaling_Factor;
    SminS = SminS*Scaling_Factor;
    SmaxS = SmaxS*Scaling_Factor;
    SminE = SminE*Scaling_Factor;
    SmaxE = SmaxE*Scaling_Factor;
    SminW = SminW*Scaling_Factor;
    SmaxW = SmaxW*Scaling_Factor;

    /* Require update of the whole mesh */
    Schedule_Interior_Mesh_Update();
    Schedule_Ghost_Cells_Update(); 
}

 
/*!
 * Rotates the quadrilateral grid block.                
 *
 * \note This subroutine DOESN'T update the geometric properties
 *       of the new grid cells!!!                               
 */
void Grid2D_Quad_Block_HO::Rotate_Quad_Block_Without_Update(const double &Angle) {

    int i, j;
    double cos_angle, sin_angle;
    Vector2D X;

    cos_angle = cos(-Angle); 
    sin_angle = sin(-Angle);

    for ( j = JNl-Nghost ; j <= JNu+Nghost; ++j ) {
       for ( i = INl-Nghost ; i <= INu+Nghost; ++i ) {
           X.x = Node[i][j].X.x*cos_angle +
                 Node[i][j].X.y*sin_angle;
           X.y = - Node[i][j].X.x*sin_angle +
                   Node[i][j].X.y*cos_angle;
           Node[i][j].X = X;
       } /* endfor */
    } /* endfor */

    if (BndNorthSpline.np != 0 ) 
       BndNorthSpline.Rotate_Spline(Angle);
    if (BndSouthSpline.np != 0 ) 
       BndSouthSpline.Rotate_Spline(Angle);
    if (BndEastSpline.np != 0 ) 
       BndEastSpline.Rotate_Spline(Angle);
    if (BndWestSpline.np != 0 )
       BndWestSpline.Rotate_Spline(Angle);

    if (ExtendWest_BndNorthSpline.np != 0) ExtendWest_BndNorthSpline.Rotate_Spline(Angle);
    if (ExtendEast_BndNorthSpline.np != 0) ExtendEast_BndNorthSpline.Rotate_Spline(Angle);
    if (ExtendWest_BndSouthSpline.np != 0) ExtendWest_BndSouthSpline.Rotate_Spline(Angle);
    if (ExtendEast_BndSouthSpline.np != 0) ExtendEast_BndSouthSpline.Rotate_Spline(Angle);
    if (ExtendNorth_BndEastSpline.np != 0) ExtendNorth_BndEastSpline.Rotate_Spline(Angle);
    if (ExtendSouth_BndEastSpline.np != 0) ExtendSouth_BndEastSpline.Rotate_Spline(Angle);
    if (ExtendNorth_BndWestSpline.np != 0) ExtendNorth_BndWestSpline.Rotate_Spline(Angle);
    if (ExtendSouth_BndWestSpline.np != 0) ExtendSouth_BndWestSpline.Rotate_Spline(Angle);

    /* Require update of the whole mesh */
    Schedule_Interior_Mesh_Update();
    Schedule_Ghost_Cells_Update();
}


/*!
 * Re-computes the locations of the nodes and cells of  
 * the quadrilateral grid block based on a mirror       
 * reflection about the y=0 axis.  The cells and nodes  
 * are also re-ordered in the i-direction.              
 *
 * \note This subroutine DOESN'T update the geometric properties
 *       of the new grid cells!!!                               
 */
void Grid2D_Quad_Block_HO::Reflect_Quad_Block_Without_Update(void) {

  int i, j;
  Vector2D *X;
  Spline2D_HO S;

  X = new Vector2D[NNi];

  for ( j = JNl-Nghost ; j <= JNu+Nghost; ++j ) {
    for ( i = INl-Nghost ; i <= INu+Nghost; ++i ) {
      X[NNi-1-i].x = Node[i][j].X.x;
      X[NNi-1-i].y = - Node[i][j].X.y;
    } /* endfor */
    for ( i = INl-Nghost ; i <= INu+Nghost; ++i ) {
      Node[i][j].X = X[i];
    } /* endfor */
  } /* endfor */

  delete []X;
  X = NULL;

  if (BndNorthSpline.np != 0 ) 
    BndNorthSpline.Reflect_Spline();
  if (BndSouthSpline.np != 0 ) 
    BndSouthSpline.Reflect_Spline();
  if (BndEastSpline.np != 0 ) 
    BndEastSpline.Reflect_Spline();
  if (BndWestSpline.np != 0 ) {
    BndWestSpline.Reflect_Spline();
    S = BndEastSpline;
    BndEastSpline = BndWestSpline;
    BndWestSpline = S;
    if (S.np != 0) S.deallocate();
  }/* endif */

  if (ExtendWest_BndNorthSpline.np != 0) ExtendWest_BndNorthSpline.Reflect_Spline();
  if (ExtendEast_BndNorthSpline.np != 0) ExtendEast_BndNorthSpline.Reflect_Spline();
  if (ExtendWest_BndSouthSpline.np != 0) ExtendWest_BndSouthSpline.Reflect_Spline();
  if (ExtendEast_BndSouthSpline.np != 0) ExtendEast_BndSouthSpline.Reflect_Spline();
  if (ExtendNorth_BndEastSpline.np != 0) ExtendNorth_BndEastSpline.Reflect_Spline();
  if (ExtendSouth_BndEastSpline.np != 0) ExtendSouth_BndEastSpline.Reflect_Spline();
  if (ExtendNorth_BndWestSpline.np != 0) ExtendNorth_BndWestSpline.Reflect_Spline();
  if (ExtendSouth_BndWestSpline.np != 0) ExtendSouth_BndWestSpline.Reflect_Spline();

  /* Require update of the whole mesh */
  Schedule_Interior_Mesh_Update();
  Schedule_Ghost_Cells_Update();

  /* Reset the boundary condition types at the quadrilateral 
     grid block boundaries. */
  
  Set_BCs();

}

/*
 * Disturb randomly the interior nodes for the         
 * quadrilateral mesh block. This routine uses the      
 * default seed of the random number generator in order 
 * to determine the angle to which the node is moved.   
 *
 * \note This subroutine DOESN'T update the geometric properties
 *       of the new grid cells!!!                               
 */
void Grid2D_Quad_Block_HO::Disturb_Interior_Nodes_Without_Update(const int &Number_of_Iterations) {

  int i,j;
  double MinDistance, Angle, Displacement;

  /* Displace the interior nodes of the quadrilateral mesh block without affecting the boundary nodes */
  for (int num_iter=1; num_iter<=Number_of_Iterations; ++num_iter){

    for ( j = JNl+1 ; j <= JNu-1 ; ++j) {
      for ( i = INl+1 ; i <= INu-1 ; ++i) {

	// Determine the minimum distance between the Node[i][j] and all the neighbour edges
	MinDistance = MinimumNodeEdgeDistance(i,j);

	// Generate a random angle for the displacement direction
	Angle = 2.0*PI*drand48();

	// Calculate the displacement -> 10% of the MinDistance
	Displacement = 0.1 * MinDistance;

	// Calculate the new node location
	Node[i][j].setloc(Node[i][j].x() + Displacement*cos(Angle),
			  Node[i][j].y() + Displacement*sin(Angle)); 
      }	// endfor (i)
    } // endfor (j)
  } //endfor (num_iter)

  /* Require update of the interior cells geometric properties. */
  Schedule_Interior_Mesh_Update();
}

/*!
 * Writes the nodes of the quadrilateral mesh to the    
 * specified output stream in a format suitable for     
 * plotting the grid with TECPLOT.                      
 *                                                      
 */
void Grid2D_Quad_Block_HO::Output_Tecplot(const int Block_Number,
					  const int Output_Title,
					  ostream &Out_File) const {
  
  int i, j;
  
  Out_File << setprecision(14);
  if (Output_Title) {
    Out_File << "TITLE = \"" << CFFC_Name()
	     << ": 2D Structured Curvilinear Grid Block (Node Locations)"
	     << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
	     << "\"y\" \n"
	     << "ZONE T =  \"Block Number = " << Block_Number
	     << "\" \\ \n"
	     << "I = " << INu - INl + 1 << " \\ \n"
	     << "J = " << JNu - JNl + 1 << " \\ \n"
	     << "F = POINT \n";

  } else {
    Out_File << "ZONE T =  \"Block Number = " << Block_Number
	     << "\" \\ \n"
	     << "I = " << INu - INl + 1 << " \\ \n"
	     << "J = " << JNu - JNl + 1 << " \\ \n"
	     << "F = POINT \n";
  } /* endif */

  if (Tecplot_Execution_Mode::IsDoublePrecision()){
    Out_File << "DT = (DOUBLE DOUBLE) \n";
  }

  for (j  = JNl ; j <= JNu ; ++j ) {
    for ( i = INl ; i <= INu ; ++i ) {
      Out_File << " " << Node[i][j].X << "\n";
    } /* endfor */
  } /* endfor */
  Out_File << setprecision(6);

}


/*!
 * Writes the nodes of the quadrilateral mesh to the    
 * specified output stream in a format suitable for     
 * plotting the grid with TECPLOT.  Includes boundary   
 * nodes.                                               
 */
void Grid2D_Quad_Block_HO::Output_Nodes_Tecplot(const int Block_Number,
						const int Output_Title,
						ostream &Out_File) const {

  int i, j;

  Out_File << setprecision(14);
  if (Output_Title) {
    Out_File << "TITLE = \"" << CFFC_Name()
	     << ": 2D Structured Curvilinear Grid Block (Node Locations)"
	     << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
	     << "\"y\" \n"
	     << "ZONE T =  \"Block Number = " << Block_Number
	     << "\" \\ \n"
	     << "I = " << INu - INl + 1 + 2*Nghost << " \\ \n"
	     << "J = " << JNu - JNl + 1 + 2*Nghost << " \\ \n"
	     << "F = POINT \n";

  } else {
    Out_File << "ZONE T =  \"Block Number = " << Block_Number
	     << "\" \\ \n"
	     << "I = " << INu - INl + 1 + 2*Nghost << " \\ \n"
	     << "J = " << JNu - JNl + 1 + 2*Nghost << " \\ \n"
	     << "F = POINT \n";
  } /* endif */

  if (Tecplot_Execution_Mode::IsDoublePrecision()){
    Out_File << "DT = (DOUBLE DOUBLE) \n";
  }

  for (j  = JNl-Nghost ; j <= JNu+Nghost ; ++j ) {
    for ( i = INl-Nghost ; i <= INu+Nghost ; ++i ) {
      Out_File << " " << Node[i][j].X << "\n";
    } /* endfor */
  } /* endfor */
  Out_File << setprecision(6);

}

/*!
 * Writes the cells of the quadrilateral mesh to the    
 * specified output stream in a format suitable for     
 * plotting the grid with TECPLOT.                      
 */
void Grid2D_Quad_Block_HO::Output_Cells_Tecplot(const int Block_Number,
						const int Output_Title,
						ostream &Out_File) const {

  int i, j;

  Out_File << setprecision(14);
  if (Output_Title) {
    Out_File << "TITLE = \"" << CFFC_Name()
	     << ": 2D Structured Curvilinear Grid Block (Cell Locations)"
	     << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
	     << "\"y\" \n"
	     << "ZONE T =  \"Block Number = " << Block_Number
	     << "\" \\ \n"
	     << "I = " << ICu - ICl + 1 << " \\ \n"
	     << "J = " << JCu - JCl + 1 << " \\ \n"
	     << "F = POINT \n";

  } else {
    Out_File << "ZONE T =  \"Block Number = " << Block_Number
	     << "\" \\ \n"
	     << "I = " << ICu - ICl + 1 << " \\ \n"
	     << "J = " << JCu - JCl + 1 << " \\ \n"
	     << "F = POINT \n";
  } /* endif */

  if (Tecplot_Execution_Mode::IsDoublePrecision()){
    Out_File << "DT = (DOUBLE DOUBLE) \n";
  }

  for (j  = JCl ; j <= JCu ; ++j ) {
    for ( i = ICl ; i <= ICu ; ++i ) {
      Out_File << " " << Cell[i][j].Xc << "\n";
    } /* endfor */
  } /* endfor */
  Out_File << setprecision(6);

}

/*!
 * Writes the nodes of the quadrilateral mesh to the    
 * specified output stream in a format suitable for     
 * plotting the grid with GNUPLOT.                      
 */
void Grid2D_Quad_Block_HO::Output_Gnuplot(const int Block_Number,
					  const int Output_Title,
					  ostream &Out_File) const {

  int i, j;

  Out_File << setprecision(14);
  if (Output_Title) {
    Out_File << "# " << CFFC_Name()
	     << ": 2D Structured Curvilinear Grid Block (Node Locations)"
	     << "\n"
	     << "# x(m), y(m)\n";
  } /* endif */

  for (j  = JNl ; j <= JNu ; ++j ) {
    for ( i = INl ; i <= INu ; ++i ) {
      Out_File << " " << Node[i][j].X << "\n";
    } /* endfor */
    Out_File << "\n";
  } /* endfor */
  for (i  = INl ; i <= INu ; ++i ) {
    for ( j = JNl ; j <= JNu ; ++j ) {
      Out_File << " " << Node[i][j].X << "\n";
    } /* endfor */
    Out_File << "\n";
  } /* endfor */
  Out_File << setprecision(6);

}



/*!
 * Returns a new quadrilateral mesh block with twice    
 * the mesh resolution of the input grid block.         
 *                                                      
 */
void Grid2D_Quad_Block_HO::Double_Mesh_Resolution(const Grid2D_Quad_Block_HO &Grid_Original) {

  int i, j, double_resolution_permitted;
  double sp_l, sp_r, sp_m, ds_ratio;
 
  /* Allocate memory for the cells and nodes of the 
     quadrilateral mesh block with twice the resolution. */

  if ( (Grid_Original.NCi-2*Grid_Original.Nghost-2*((Grid_Original.NCi-2*Grid_Original.Nghost)/2) != 0) || 
       (Grid_Original.NCj-2*Grid_Original.Nghost-2*((Grid_Original.NCj-2*Grid_Original.Nghost)/2) != 0) ||
       (Grid_Original.NCi-2*Grid_Original.Nghost < 2*Grid_Original.Nghost) ||
       (Grid_Original.NCj-2*Grid_Original.Nghost < 2*Grid_Original.Nghost) ||
       (Grid_Original.Node == NULL) ) { 
    double_resolution_permitted = 0;
  } else {
    double_resolution_permitted = 1;
    allocate(2*(Grid_Original.NCi-2*Grid_Original.Nghost), 
	     2*(Grid_Original.NCj-2*Grid_Original.Nghost),
	     Grid_Original.Nghost,
	     Grid_Original.MaxRecOrder());
  }/* endif */

  /* Copy boundary spline info to quadrilateral mesh block 
     with twice the resolution. */

  if (double_resolution_permitted) {

    if (Grid_Original.BndNorthSpline.np != 0) {
      BndNorthSpline = Grid_Original.BndNorthSpline;
    } else if (BndNorthSpline.np != 0) {
      BndNorthSpline.deallocate();
      deallocate_BndNorthSplineInfo();
    } /* endif */

    if (Grid_Original.BndSouthSpline.np != 0) {
      BndSouthSpline = Grid_Original.BndSouthSpline;
    } else if (BndSouthSpline.np != 0) {
      BndSouthSpline.deallocate();
      deallocate_BndSouthSplineInfo();
    } /* endif */

    if (Grid_Original.BndEastSpline.np != 0) {
      BndEastSpline = Grid_Original.BndEastSpline;
    } else if (BndEastSpline.np != 0) {
      BndEastSpline.deallocate();
      deallocate_BndEastSplineInfo();
    } /* endif */
  
    if (Grid_Original.BndWestSpline.np != 0) {
      BndWestSpline = Grid_Original.BndWestSpline;
    } else if (BndWestSpline.np != 0) {
      BndWestSpline.deallocate();
      deallocate_BndWestSplineInfo();
    } /* endif */

    // Copy the extensions to boundary splines
    if (Grid_Original.ExtendWest_BndNorthSpline.np != 0) {
      ExtendWest_BndNorthSpline = Grid_Original.ExtendWest_BndNorthSpline;
    } else if (ExtendWest_BndNorthSpline.np != 0) {
      ExtendWest_BndNorthSpline.deallocate();
    } /* endif */
    if (Grid_Original.ExtendEast_BndNorthSpline.np != 0) {
      ExtendEast_BndNorthSpline = Grid_Original.ExtendEast_BndNorthSpline;
    } else if (ExtendEast_BndNorthSpline.np != 0) {
      ExtendEast_BndNorthSpline.deallocate();
    } /* endif */

    if (Grid_Original.ExtendWest_BndSouthSpline.np != 0) {
      ExtendWest_BndSouthSpline = Grid_Original.ExtendWest_BndSouthSpline;
    } else if (ExtendWest_BndSouthSpline.np != 0) {
      ExtendWest_BndSouthSpline.deallocate();
    } /* endif */
    if (Grid_Original.ExtendEast_BndSouthSpline.np != 0) {
      ExtendEast_BndSouthSpline = Grid_Original.ExtendEast_BndSouthSpline;
    } else if (ExtendEast_BndSouthSpline.np != 0) {
      ExtendEast_BndSouthSpline.deallocate();
    } /* endif */

    if (Grid_Original.ExtendNorth_BndEastSpline.np != 0) {
      ExtendNorth_BndEastSpline = Grid_Original.ExtendNorth_BndEastSpline;
    } else if (ExtendNorth_BndEastSpline.np != 0) {
      ExtendNorth_BndEastSpline.deallocate();
    } /* endif */
    if (Grid_Original.ExtendSouth_BndEastSpline.np != 0) {
      ExtendSouth_BndEastSpline = Grid_Original.ExtendSouth_BndEastSpline;
    } else if (ExtendSouth_BndEastSpline.np != 0) {
      ExtendSouth_BndEastSpline.deallocate();
    } /* endif */
  
    if (Grid_Original.ExtendNorth_BndWestSpline.np != 0) {
      ExtendNorth_BndWestSpline = Grid_Original.ExtendNorth_BndWestSpline;
    } else if (ExtendNorth_BndWestSpline.np != 0) {
      ExtendNorth_BndWestSpline.deallocate();
    } /* endif */
    if (Grid_Original.ExtendSouth_BndWestSpline.np != 0) {
      ExtendSouth_BndWestSpline = Grid_Original.ExtendSouth_BndWestSpline;
    } else if (ExtendSouth_BndWestSpline.np != 0) {
      ExtendSouth_BndWestSpline.deallocate();
    } /* endif */


    /* Copy boundary spline pathlength info to quadrilateral mesh block 
       with twice the resolution. */

    SminN = Grid_Original.SminN;
    SmaxN = Grid_Original.SmaxN;
    SminS = Grid_Original.SminS;
    SmaxS = Grid_Original.SmaxS;
    SminE = Grid_Original.SminE;
    SmaxE = Grid_Original.SmaxE;
    SminW = Grid_Original.SminW;
    SmaxW = Grid_Original.SmaxW;

    /* Copy node stretching info to quadrilateral mesh block 
       with twice the resolution. */

    StretchI = Grid_Original.StretchI;
    BetaI = Grid_Original.BetaI;
    TauI = Grid_Original.TauI;
    StretchJ = Grid_Original.StretchJ;
    BetaJ = Grid_Original.BetaJ;
    TauJ = Grid_Original.TauJ;
    OrthogonalN = Grid_Original.OrthogonalN;
    OrthogonalS = Grid_Original.OrthogonalS;
    OrthogonalE = Grid_Original.OrthogonalE;
    OrthogonalW = Grid_Original.OrthogonalW;

    /* Determine the node locations of quadrilateral mesh block 
       with twice the resolution. */

    for (j  = Grid_Original.JCl ; j <= Grid_Original.JCu ; ++j ) {
      for ( i = Grid_Original.ICl ; i <= Grid_Original.ICu ; ++i ) {
	Node[2*(i-Grid_Original.INl)+Grid_Original.INl  ]
	  [2*(j-Grid_Original.JNl)+Grid_Original.JNl  ].X 
	  = Grid_Original.nodeSW(i, j).X;
	Node[2*(i-Grid_Original.INl)+Grid_Original.INl+1]
	  [2*(j-Grid_Original.JNl)+Grid_Original.JNl  ].X 
	  = Grid_Original.xfaceS(i, j);
	Node[2*(i-Grid_Original.INl)+Grid_Original.INl  ]
	  [2*(j-Grid_Original.JNl)+Grid_Original.JNl+1].X 
	  = Grid_Original.xfaceW(i, j);
	Node[2*(i-Grid_Original.INl)+Grid_Original.INl+1]
	  [2*(j-Grid_Original.JNl)+Grid_Original.JNl+1].X 
	  = Grid_Original.Cell[i][j].Xc;
	if (j == Grid_Original.JCu) {
	  Node[2*(i-Grid_Original.INl)+Grid_Original.INl  ]
	    [2*(j-Grid_Original.JNl)+Grid_Original.JNl+2].X 
	    = Grid_Original.nodeNW(i, j).X;
	  Node[2*(i-Grid_Original.INl)+Grid_Original.INl+1]
	    [2*(j-Grid_Original.JNl)+Grid_Original.JNl+2].X 
	    = Grid_Original.xfaceN(i, j);
	} /* endif */
	if (i == Grid_Original.ICu) {
	  Node[2*(i-Grid_Original.INl)+Grid_Original.INl+2]
	    [2*(j-Grid_Original.JNl)+Grid_Original.JNl  ].X 
	    = Grid_Original.nodeSE(i, j).X;
	  Node[2*(i-Grid_Original.INl)+Grid_Original.INl+2]
	    [2*(j-Grid_Original.JNl)+Grid_Original.JNl+1].X 
	    = Grid_Original.xfaceE(i, j);
	} /* endif */
	if (i == Grid_Original.ICu && j == Grid_Original.JCu) {
	  Node[2*(i-Grid_Original.INl)+Grid_Original.INl+2]
	    [2*(j-Grid_Original.JNl)+Grid_Original.JNl+2].X 
	    = Grid_Original.nodeNE(i, j).X;
	} /* endif */
      } /* endfor */
    } /* endfor */

    if (BndWestSpline.np != 0) {
      for (j  = JNl+1 ; j < JNu ; j += 2 ) {
	sp_l = getS(Node[INl][j-1].X, 
		    BndWestSpline);
	sp_r = getS(Node[INl][j+1].X, 
		    BndWestSpline);
	ds_ratio = abs(Node[INl+1][j].X-
		       Node[INl+1][j-1].X)/
	  abs(Node[INl+1][j+1].X-
	      Node[INl+1][j-1].X);
	sp_m = sp_l + ds_ratio*(sp_r-sp_l);
	Node[INl][j].X = 
	  Spline(sp_m, BndWestSpline);
      } /* endfor */
    } /* endif */

    if (BndEastSpline.np != 0) {
      for (j  = JNl+1 ; j < JNu ; j += 2 ) {
	sp_l = getS(Node[INu][j-1].X, 
		    BndEastSpline);
	sp_r = getS(Node[INu][j+1].X, 
		    BndEastSpline);
	ds_ratio = abs(Node[INu-1][j].X-
		       Node[INu-1][j-1].X)/
	  abs(Node[INu-1][j+1].X-
	      Node[INu-1][j-1].X);
	sp_m = sp_l + ds_ratio*(sp_r-sp_l);
	Node[INu][j].X = 
	  Spline(sp_m, BndEastSpline);
      } /* endfor */
    } /* endif */

    if (BndSouthSpline.np != 0) {
      for ( i = INl+1 ; i < INu ; i += 2 ) {
	sp_l = getS(Node[i-1][JNl].X, 
		    BndSouthSpline);
	sp_r = getS(Node[i+1][JNl].X, 
		    BndSouthSpline);
	ds_ratio = abs(Node[i][JNl+1].X-
		       Node[i-1][JNl+1].X)/
	  abs(Node[i+1][JNl+1].X-
	      Node[i-1][JNl+1].X);
	sp_m = sp_l + ds_ratio*(sp_r-sp_l);
	Node[i][JNl].X =  
	  Spline(sp_m, BndSouthSpline);
      } /* endfor */
    } /* endif */

    if (BndNorthSpline.np != 0) {
      for ( i = INl+1 ; i < INu ; i += 2 ) {
	sp_l = getS(Node[i-1][JNu].X, 
		    BndNorthSpline);
	sp_r = getS(Node[i+1][JNu].X, 
		    BndNorthSpline);
	ds_ratio = abs(Node[i][JNu-1].X-
		       Node[i-1][JNu-1].X)/
	  abs(Node[i+1][JNu-1].X-
	      Node[i-1][JNu-1].X);
	sp_m = sp_l + ds_ratio*(sp_r-sp_l);
	Node[i][JNu].X =  
	  Spline(sp_m, BndNorthSpline);
      } /* endfor */
    } /* endif */

    /* Require update of the interior cells geometric properties. */
    Schedule_Interior_Mesh_Update();

    /* Set the boundary condition types for quadrilateral mesh block 
       with twice the resolution. */

    Set_BCs();

    /* Compute the exterior nodes for quadrilateral mesh block 
       with twice the resolution. */

    Update_Exterior_Nodes();

    /* Compute the cells for quadrilateral mesh block 
       with twice the resolution. */

    Update_Cells();

  } /* endif */

}

/*!
 * Returns a new quadrilateral mesh block with half the 
 * mesh resolution of the input grid block.             
 *                                                      
 */
void Grid2D_Quad_Block_HO::Half_Mesh_Resolution(const Grid2D_Quad_Block_HO &Grid_Original) {

  int i, j, half_resolution_permitted;
 
  /* Allocate memory for the cells and nodes of the 
     quadrilateral mesh block with half the resolution. */

  if ( (Grid_Original.NCi-2*Grid_Original.Nghost-2*((Grid_Original.NCi-2*Grid_Original.Nghost)/2) != 0) || 
       (Grid_Original.NCj-2*Grid_Original.Nghost-2*((Grid_Original.NCj-2*Grid_Original.Nghost)/2) != 0) ||
       (Grid_Original.NCi-2*Grid_Original.Nghost < 2*Grid_Original.Nghost) ||
       (Grid_Original.NCj-2*Grid_Original.Nghost < 2*Grid_Original.Nghost) ||
       (Grid_Original.Node == NULL) ) {
    half_resolution_permitted = 0;
  } else {
    half_resolution_permitted = 1;
    allocate((Grid_Original.NCi-2*Grid_Original.Nghost)/2, 
	     (Grid_Original.NCj-2*Grid_Original.Nghost)/2,
	     Grid_Original.Nghost,
	     Grid_Original.MaxRecOrder());
  } /* endif */

    /* Copy boundary spline info to quadrilateral mesh block 
       with half the resolution. */

  if (half_resolution_permitted) {

    if (Grid_Original.BndNorthSpline.np != 0) {
      BndNorthSpline = Grid_Original.BndNorthSpline;
    } else if (BndNorthSpline.np != 0) {
      BndNorthSpline.deallocate();
      deallocate_BndNorthSplineInfo();
    } /* endif */

    if (Grid_Original.BndSouthSpline.np != 0) {
      BndSouthSpline = Grid_Original.BndSouthSpline;
    } else if (BndSouthSpline.np != 0) {
      BndSouthSpline.deallocate();
      deallocate_BndSouthSplineInfo();
    } /* endif */

    if (Grid_Original.BndEastSpline.np != 0) {
      BndEastSpline = Grid_Original.BndEastSpline;
    } else if (BndEastSpline.np != 0) {
      BndEastSpline.deallocate();
      deallocate_BndEastSplineInfo();
    } /* endif */
  
    if (Grid_Original.BndWestSpline.np != 0) {
      BndWestSpline = Grid_Original.BndWestSpline;
    } else if (BndWestSpline.np != 0) {
      BndWestSpline.deallocate();
      deallocate_BndWestSplineInfo();
    } /* endif */

    // Copy the extensions to boundary splines
    if (Grid_Original.ExtendWest_BndNorthSpline.np != 0) {
      ExtendWest_BndNorthSpline = Grid_Original.ExtendWest_BndNorthSpline;
    } else if (ExtendWest_BndNorthSpline.np != 0) {
      ExtendWest_BndNorthSpline.deallocate();
    } /* endif */
    if (Grid_Original.ExtendEast_BndNorthSpline.np != 0) {
      ExtendEast_BndNorthSpline = Grid_Original.ExtendEast_BndNorthSpline;
    } else if (ExtendEast_BndNorthSpline.np != 0) {
      ExtendEast_BndNorthSpline.deallocate();
    } /* endif */

    if (Grid_Original.ExtendWest_BndSouthSpline.np != 0) {
      ExtendWest_BndSouthSpline = Grid_Original.ExtendWest_BndSouthSpline;
    } else if (ExtendWest_BndSouthSpline.np != 0) {
      ExtendWest_BndSouthSpline.deallocate();
    } /* endif */
    if (Grid_Original.ExtendEast_BndSouthSpline.np != 0) {
      ExtendEast_BndSouthSpline = Grid_Original.ExtendEast_BndSouthSpline;
    } else if (ExtendEast_BndSouthSpline.np != 0) {
      ExtendEast_BndSouthSpline.deallocate();
    } /* endif */

    if (Grid_Original.ExtendNorth_BndEastSpline.np != 0) {
      ExtendNorth_BndEastSpline = Grid_Original.ExtendNorth_BndEastSpline;
    } else if (ExtendNorth_BndEastSpline.np != 0) {
      ExtendNorth_BndEastSpline.deallocate();
    } /* endif */
    if (Grid_Original.ExtendSouth_BndEastSpline.np != 0) {
      ExtendSouth_BndEastSpline = Grid_Original.ExtendSouth_BndEastSpline;
    } else if (ExtendSouth_BndEastSpline.np != 0) {
      ExtendSouth_BndEastSpline.deallocate();
    } /* endif */
  
    if (Grid_Original.ExtendNorth_BndWestSpline.np != 0) {
      ExtendNorth_BndWestSpline = Grid_Original.ExtendNorth_BndWestSpline;
    } else if (ExtendNorth_BndWestSpline.np != 0) {
      ExtendNorth_BndWestSpline.deallocate();
    } /* endif */
    if (Grid_Original.ExtendSouth_BndWestSpline.np != 0) {
      ExtendSouth_BndWestSpline = Grid_Original.ExtendSouth_BndWestSpline;
    } else if (ExtendSouth_BndWestSpline.np != 0) {
      ExtendSouth_BndWestSpline.deallocate();
    } /* endif */

    /* Copy boundary spline pathlength info to quadrilateral mesh block 
       with half the resolution. */

    SminN = Grid_Original.SminN;
    SmaxN = Grid_Original.SmaxN;
    SminS = Grid_Original.SminS;
    SmaxS = Grid_Original.SmaxS;
    SminE = Grid_Original.SminE;
    SmaxE = Grid_Original.SmaxE;
    SminW = Grid_Original.SminW;
    SmaxW = Grid_Original.SmaxW;

    /* Copy node stretching info to quadrilateral mesh block 
       with half the resolution. */

    StretchI = Grid_Original.StretchI;
    BetaI = Grid_Original.BetaI;
    TauI = Grid_Original.TauI;
    StretchJ = Grid_Original.StretchJ;
    BetaJ = Grid_Original.BetaJ;
    TauJ = Grid_Original.TauJ;
    OrthogonalN = Grid_Original.OrthogonalN;
    OrthogonalS = Grid_Original.OrthogonalS;
    OrthogonalE = Grid_Original.OrthogonalE;
    OrthogonalW = Grid_Original.OrthogonalW;

    /* Determine the node locations of quadrilateral mesh block 
       with half the resolution. */

    for (j  = JNl ; j <= JNu ; ++j ) {
      for ( i = INl ; i <= INu ; ++i ) {
	Node[i][j].X = 
	  Grid_Original.Node[2*(i-INl)+INl]
	  [2*(j-JNl)+JNl].X;
      } /* endfor */
    } /* endfor */

    /* Require update of the interior cells geometric properties. */
    Schedule_Interior_Mesh_Update();

    /* Set the boundary condition types for quadrilateral mesh block 
       with half the resolution. */

    Set_BCs();

    /* Compute the exterior nodes for quadrilateral mesh block 
       with half the resolution. */

    Update_Exterior_Nodes();

    /* Compute the cells for quadrilateral mesh block 
       with half the resolution. */

    Update_Cells();

  } /* endif */

}

/*!
 * Returns a new quadrilateral mesh block containing    
 * one of the specified sectors of the original grid    
 * block with twice the mesh resolution.                
 *                                                      
 */
void Grid2D_Quad_Block_HO::Refine_Mesh(const Grid2D_Quad_Block_HO &Grid_Original,
				       const int Sector) {

  int i, j, i_min, i_max, j_min, j_max, mesh_refinement_permitted;
  double sp_l, sp_r, sp_m, ds_ratio, dl, dr;

  /* Allocate memory for the cells and nodes for the 
     refined quadrilateral mesh block. */

  if ( (Grid_Original.NCi-2*Grid_Original.Nghost-2*((Grid_Original.NCi-2*Grid_Original.Nghost)/2) != 0) || 
       (Grid_Original.NCj-2*Grid_Original.Nghost-2*((Grid_Original.NCj-2*Grid_Original.Nghost)/2) != 0) ||
       (Grid_Original.NCi-2*Grid_Original.Nghost < 2*Grid_Original.Nghost) ||
       (Grid_Original.NCj-2*Grid_Original.Nghost < 2*Grid_Original.Nghost) ||
       (Grid_Original.Node == NULL) ) {
    mesh_refinement_permitted = 0;
  } else {
    mesh_refinement_permitted = 1;
    allocate(Grid_Original.NCi-2*Grid_Original.Nghost, 
	     Grid_Original.NCj-2*Grid_Original.Nghost,
	     Grid_Original.Nghost,
	     Grid_Original.MaxRecOrder());
  } /* endif */

    /* Copy boundary spline info for the refined
       quadrilateral mesh block. */

  if (mesh_refinement_permitted) {

    // NW Sector
    if (Sector == GRID2D_QUAD_BLOCK_SECTOR_NW) {

      // North spline
      if (Grid_Original.BndNorthSpline.np != 0) {
	BndNorthSpline = Grid_Original.BndNorthSpline;
      } else if (BndNorthSpline.np != 0) {
	BndNorthSpline.deallocate();
	deallocate_BndNorthSplineInfo();
      }

      // Copy the extension of the original North spline to West
      if (Grid_Original.ExtendWest_BndNorthSpline.np != 0){
	ExtendWest_BndNorthSpline = Grid_Original.ExtendWest_BndNorthSpline;
      } else if (ExtendWest_BndNorthSpline.np != 0){
	ExtendWest_BndNorthSpline.deallocate();
      }

      // Extend the North spline to East with itself if necessary
      if (BndNorthSpline.np != 0){
	ExtendEast_BndNorthSpline = BndNorthSpline;
      } else if (ExtendEast_BndNorthSpline.np != 0){
	ExtendEast_BndNorthSpline.deallocate();
      }

      // West Spline
      if (Grid_Original.BndWestSpline.np != 0){
	BndWestSpline = Grid_Original.BndWestSpline;
      } else if (BndWestSpline.np != 0){
	BndWestSpline.deallocate();
	deallocate_BndWestSplineInfo();
      }

      // Copy the extension of the original West spline to North
      if (Grid_Original.ExtendNorth_BndWestSpline.np != 0){
	ExtendNorth_BndWestSpline = Grid_Original.ExtendNorth_BndWestSpline;
      } else if (ExtendNorth_BndWestSpline.np != 0){
	ExtendNorth_BndWestSpline.deallocate();
      }

      // Extend the West spline to South with itself if necessary
      if (BndWestSpline.np != 0){
	ExtendSouth_BndWestSpline = BndWestSpline;
      } else if (ExtendSouth_BndWestSpline.np != 0){
	ExtendSouth_BndWestSpline.deallocate();
      }
    } /* endif */

    // NE Sector
    if (Sector == GRID2D_QUAD_BLOCK_SECTOR_NE) {

      // North spline
      if (Grid_Original.BndNorthSpline.np != 0) {
	BndNorthSpline = Grid_Original.BndNorthSpline;
      } else if (BndNorthSpline.np != 0) {
	BndNorthSpline.deallocate();
	deallocate_BndNorthSplineInfo();
      }

      // Copy the extension of the original North spline to East
      if (Grid_Original.ExtendEast_BndNorthSpline.np != 0){
	ExtendEast_BndNorthSpline = Grid_Original.ExtendEast_BndNorthSpline;
      } else if (ExtendEast_BndNorthSpline.np != 0){
	ExtendEast_BndNorthSpline.deallocate();
      }

      // Extend the North spline to West with itself if necessary
      if (BndNorthSpline.np != 0){
	ExtendWest_BndNorthSpline = BndNorthSpline;
      } else if (ExtendWest_BndNorthSpline.np != 0){
	ExtendWest_BndNorthSpline.deallocate();
      }

      // East Spline
      if (Grid_Original.BndEastSpline.np != 0){
	BndEastSpline = Grid_Original.BndEastSpline;
      } else if (BndEastSpline.np != 0){
	BndEastSpline.deallocate();
	deallocate_BndEastSplineInfo();
      }

      // Copy the extension of the original East spline to North
      if (Grid_Original.ExtendNorth_BndEastSpline.np != 0){
	ExtendNorth_BndEastSpline = Grid_Original.ExtendNorth_BndEastSpline;
      } else if (ExtendNorth_BndEastSpline.np != 0){
	ExtendNorth_BndEastSpline.deallocate();
      }

      // Extend the East spline to South with itself if necessary
      if (BndEastSpline.np != 0){
	ExtendSouth_BndEastSpline = BndEastSpline;
      } else if (ExtendSouth_BndEastSpline.np != 0){
	ExtendSouth_BndEastSpline.deallocate();
      }
    } /* endif */

    // SW Sector
    if (Sector == GRID2D_QUAD_BLOCK_SECTOR_SW) {

      // South spline
      if (Grid_Original.BndSouthSpline.np != 0) {
	BndSouthSpline = Grid_Original.BndSouthSpline;
      } else if (BndSouthSpline.np != 0) {
	BndSouthSpline.deallocate();
	deallocate_BndSouthSplineInfo();
      }

      // Copy the extension of the original South spline to West
      if (Grid_Original.ExtendWest_BndSouthSpline.np != 0){
	ExtendWest_BndSouthSpline = Grid_Original.ExtendWest_BndSouthSpline;
      } else if (ExtendWest_BndSouthSpline.np != 0){
	ExtendWest_BndSouthSpline.deallocate();
      }

      // Extend the South spline to East with itself if necessary
      if (BndSouthSpline.np != 0){
	ExtendEast_BndSouthSpline = BndSouthSpline;
      } else if (ExtendEast_BndSouthSpline.np != 0){
	ExtendEast_BndSouthSpline.deallocate();
      }

      // West Spline
      if (Grid_Original.BndWestSpline.np != 0){
	BndWestSpline = Grid_Original.BndWestSpline;
      } else if (BndWestSpline.np != 0){
	BndWestSpline.deallocate();
	deallocate_BndWestSplineInfo();
      }

      // Copy the extension of the original West spline to South
      if (Grid_Original.ExtendSouth_BndWestSpline.np != 0){
	ExtendSouth_BndWestSpline = Grid_Original.ExtendSouth_BndWestSpline;
      } else if (ExtendSouth_BndWestSpline.np != 0){
	ExtendSouth_BndWestSpline.deallocate();
      }

      // Extend the West spline to North with itself if necessary
      if (BndWestSpline.np != 0){
	ExtendNorth_BndWestSpline = BndWestSpline;
      } else if (ExtendNorth_BndWestSpline.np != 0){
	ExtendNorth_BndWestSpline.deallocate();
      }
    } /* endif */
    
    // SE Sector
    if (Sector == GRID2D_QUAD_BLOCK_SECTOR_SE) {

      // South spline
      if (Grid_Original.BndSouthSpline.np != 0) {
	BndSouthSpline = Grid_Original.BndSouthSpline;
      } else if (BndSouthSpline.np != 0) {
	BndSouthSpline.deallocate();
	deallocate_BndSouthSplineInfo();
      }

      // Copy the extension of the original South spline to East
      if (Grid_Original.ExtendEast_BndSouthSpline.np != 0){
	ExtendEast_BndSouthSpline = Grid_Original.ExtendEast_BndSouthSpline;
      } else if (ExtendEast_BndSouthSpline.np != 0){
	ExtendEast_BndSouthSpline.deallocate();
      }

      // Extend the South spline to West with itself if necessary
      if (BndSouthSpline.np != 0){
	ExtendWest_BndSouthSpline = BndSouthSpline;
      } else if (ExtendWest_BndSouthSpline.np != 0){
	ExtendWest_BndSouthSpline.deallocate();
      }

      // East Spline
      if (Grid_Original.BndEastSpline.np != 0){
	BndEastSpline = Grid_Original.BndEastSpline;
      } else if (BndEastSpline.np != 0){
	BndEastSpline.deallocate();
	deallocate_BndEastSplineInfo();
      }

      // Copy the extension of the original East spline to South
      if (Grid_Original.ExtendSouth_BndEastSpline.np != 0){
	ExtendSouth_BndEastSpline = Grid_Original.ExtendSouth_BndEastSpline;
      } else if (ExtendSouth_BndEastSpline.np != 0){
	ExtendSouth_BndEastSpline.deallocate();
      }

      // Extend the East spline to North with itself if necessary
      if (BndEastSpline.np != 0){
	ExtendNorth_BndEastSpline = BndEastSpline;
      } else if (ExtendNorth_BndEastSpline.np != 0){
	ExtendNorth_BndEastSpline.deallocate();
      }
    } /* endif */


    /* Assign boundary spline pathlength info for the refined
       quadrilateral mesh block. */

    switch(Sector) {
    case GRID2D_QUAD_BLOCK_SECTOR_NW :
      if (Grid_Original.BndNorthSpline.np != 0) {
	SminN = Grid_Original.SminN;
	SmaxN = getS(Grid_Original.Node[Grid_Original.INl+
					(Grid_Original.INu-
					 Grid_Original.INl)/2]
		     [Grid_Original.JNu].X, 
		     Grid_Original.BndNorthSpline);
      } else {
	SminN = ZERO;
	SmaxN = ZERO;
      } /* endif */
      SminS = ZERO;
      SmaxS = ZERO;
      SminE = ZERO;
      SmaxE = ZERO;
      if (Grid_Original.BndWestSpline.np != 0) {
	SminW = getS(Grid_Original.Node[Grid_Original.INl]
		     [Grid_Original.JNl+
		      (Grid_Original.JNu-
		       Grid_Original.JNl)/2].X, 
		     Grid_Original.BndWestSpline);
	SmaxW = Grid_Original.SmaxW;
      } else {
	SminW = ZERO;
	SmaxW = ZERO;
      } /* endif */
      break;
    case GRID2D_QUAD_BLOCK_SECTOR_NE :
      if (Grid_Original.BndNorthSpline.np != 0) {
	SminN = getS(Grid_Original.Node[Grid_Original.INl+
					(Grid_Original.INu-
					 Grid_Original.INl)/2]
		     [Grid_Original.JNu].X, 
		     Grid_Original.BndNorthSpline);
	SmaxN = Grid_Original.SmaxN;
      } else {
	SminN = ZERO;
	SmaxN = ZERO;
      } /* endif */
      SminS = ZERO;
      SmaxS = ZERO;
      if (Grid_Original.BndEastSpline.np != 0) {
	SminE = getS(Grid_Original.Node[Grid_Original.INu]
		     [Grid_Original.JNl+
		      (Grid_Original.JNu-
		       Grid_Original.JNl)/2].X, 
		     Grid_Original.BndEastSpline);
	SmaxE = Grid_Original.SmaxE;
      } else {
	SminE = ZERO;
	SmaxE = ZERO;
      } /* endif */
      SminW = ZERO;
      SmaxW = ZERO;
      break;
    case GRID2D_QUAD_BLOCK_SECTOR_SE :
      SminN = ZERO;
      SmaxN = ZERO;
      if (Grid_Original.BndSouthSpline.np != 0) {
	SminS = getS(Grid_Original.Node[Grid_Original.INl+
					(Grid_Original.INu-
					 Grid_Original.INl)/2]
		     [Grid_Original.JNl].X, 
		     Grid_Original.BndSouthSpline);
	SmaxS = Grid_Original.SmaxS;
      } else {
	SminS = ZERO;
	SmaxS = ZERO;
      } /* endif */
      if (Grid_Original.BndEastSpline.np != 0) {
	SminE = Grid_Original.SminE;
	SmaxE = getS(Grid_Original.Node[Grid_Original.INu]
		     [Grid_Original.JNl+
		      (Grid_Original.JNu-
		       Grid_Original.JNl)/2].X, 
		     Grid_Original.BndEastSpline);
	SmaxE = Grid_Original.SmaxE;
      } else {
	SminE = ZERO;
	SmaxE = ZERO;
      } /* endif */
      SminW = ZERO;
      SmaxW = ZERO;
      break;
    case GRID2D_QUAD_BLOCK_SECTOR_SW :
      SminN = ZERO;
      SmaxN = ZERO;
      if (Grid_Original.BndSouthSpline.np != 0) {
	SminS = Grid_Original.SminS;
	SmaxS = getS(Grid_Original.Node[Grid_Original.INl+
					(Grid_Original.INu-
					 Grid_Original.INl)/2]
		     [Grid_Original.JNl].X, 
		     Grid_Original.BndSouthSpline);
      } else {
	SminS = ZERO;
	SmaxS = ZERO;
      } /* endif */
      SminE = ZERO;
      SmaxE = ZERO;
      if (Grid_Original.BndWestSpline.np != 0) {
	SminW = Grid_Original.SminW;
	SmaxW = getS(Grid_Original.Node[Grid_Original.INl]
		     [Grid_Original.JNl+
		      (Grid_Original.JNu-
		       Grid_Original.JNl)/2].X, 
		     Grid_Original.BndWestSpline);
      } else {
	SminW = ZERO;
	SmaxW = ZERO;
      } /* endif */
      break;
    } /* endswitch */

    /* Copy node stretching info to refined quadrilateral 
       mesh block. */

    StretchI = Grid_Original.StretchI;
    BetaI = Grid_Original.BetaI;
    TauI = Grid_Original.TauI;
    StretchJ = Grid_Original.StretchJ;
    BetaJ = Grid_Original.BetaJ;
    TauJ = Grid_Original.TauJ;
    if (Grid_Original.BndNorthSpline.np != 0) {
      OrthogonalN = Grid_Original.OrthogonalN;
    } else {
      OrthogonalN = ORTHOGONAL;
    } /* endif */
    if (Grid_Original.BndSouthSpline.np != 0) {
      OrthogonalS = Grid_Original.OrthogonalS;
    } else {
      OrthogonalS = ORTHOGONAL;
    } /* endif */
    if (Grid_Original.BndEastSpline.np != 0) {
      OrthogonalE = Grid_Original.OrthogonalE;
    } else {
      OrthogonalE = ORTHOGONAL;
    } /* endif */
    if (Grid_Original.BndWestSpline.np != 0) {
      OrthogonalW = Grid_Original.OrthogonalW;
    } else {
      OrthogonalW = ORTHOGONAL;
    } /* endif */

    // Force orthogonality at all refined mesh boundaries.
    OrthogonalN = ORTHOGONAL;
    OrthogonalS = ORTHOGONAL;
    OrthogonalE = ORTHOGONAL;
    OrthogonalW = ORTHOGONAL;

    /* Determine the node locations for refined 
       quadrilateral mesh block. */

    switch(Sector) {
    case GRID2D_QUAD_BLOCK_SECTOR_NW :
      i_min = Grid_Original.ICl;
      i_max = Grid_Original.ICl+(Grid_Original.ICu-Grid_Original.ICl-1)/2;
      j_min = Grid_Original.JCl+(Grid_Original.JCu-Grid_Original.JCl+1)/2; 
      j_max = Grid_Original.JCu;
      break;
    case GRID2D_QUAD_BLOCK_SECTOR_NE :
      i_min = Grid_Original.ICl+(Grid_Original.ICu-Grid_Original.ICl+1)/2;
      i_max = Grid_Original.ICu;
      j_min = Grid_Original.JCl+(Grid_Original.JCu-Grid_Original.JCl+1)/2; 
      j_max = Grid_Original.JCu;
      break;
    case GRID2D_QUAD_BLOCK_SECTOR_SE :
      i_min = Grid_Original.ICl+(Grid_Original.ICu-Grid_Original.ICl+1)/2;
      i_max = Grid_Original.ICu;
      j_min = Grid_Original.JCl; 
      j_max = Grid_Original.JCl+(Grid_Original.JCu-Grid_Original.JCl-1)/2;
      break;
    case GRID2D_QUAD_BLOCK_SECTOR_SW :
      i_min = Grid_Original.ICl;
      i_max = Grid_Original.ICl+(Grid_Original.ICu-Grid_Original.ICl-1)/2;
      j_min = Grid_Original.JCl; 
      j_max = Grid_Original.JCl+(Grid_Original.JCu-Grid_Original.JCl-1)/2;
      break;
    default:
      i_min = Grid_Original.ICl;
      i_max = Grid_Original.ICl+(Grid_Original.ICu-Grid_Original.ICl-1)/2;
      j_min = Grid_Original.JCl+(Grid_Original.JCu-Grid_Original.JCl+1)/2; 
      j_max = Grid_Original.JCu;
      break;
    } /* endswitch */

    for ( j  = j_min; j <= j_max ; ++j ) {
      for ( i = i_min ; i <= i_max ; ++i ) {
	Node[2*(i-i_min)+INl  ][2*(j-j_min)+JNl  ].X  = Grid_Original.nodeSW(i, j).X;
	Node[2*(i-i_min)+INl+1][2*(j-j_min)+JNl  ].X  = Grid_Original.xfaceS(i, j);
	Node[2*(i-i_min)+INl  ][2*(j-j_min)+JNl+1].X  = Grid_Original.xfaceW(i, j);
	Node[2*(i-i_min)+INl+1][2*(j-j_min)+JNl+1].X  = (Grid_Original.nodeSW(i,j).X +
							 Grid_Original.nodeSE(i,j).X +
							 Grid_Original.nodeNW(i,j).X +
							 Grid_Original.nodeNE(i,j).X)/FOUR;

	if (j == j_max) {
	  Node[2*(i-i_min)+INl  ][2*(j-j_min)+JNl+2].X = Grid_Original.nodeNW(i, j).X;
	  Node[2*(i-i_min)+INl+1][2*(j-j_min)+JNl+2].X = Grid_Original.xfaceN(i, j);
	} /* endif */
	if (i == i_max) {
	  Node[2*(i-i_min)+INl+2][2*(j-j_min)+JNl  ].X = Grid_Original.nodeSE(i, j).X;
	  Node[2*(i-i_min)+INl+2][2*(j-j_min)+JNl+1].X = Grid_Original.xfaceE(i, j);
	} /* endif */
	if (i == i_max && j == j_max) {
	  Node[2*(i-i_min)+INl+2][2*(j-j_min)+JNl+2].X = Grid_Original.nodeNE(i, j).X;
	} /* endif */
      } /* endfor */
    } /* endfor */

    if (BndWestSpline.np != 0) {
      for (j  = JNl+1 ; j < JNu ; j += 2 ) {
	sp_l = getS(Node[INl][j-1].X, 
		    BndWestSpline);
	sp_r = getS(Node[INl][j+1].X, 
		    BndWestSpline);
	dl = abs(Node[INl][j  ].X - 
		 Node[INl][j-1].X);
	dr = abs(Node[INl][j+1].X - 
		 Node[INl][j  ].X);
	ds_ratio = dl/(dl+dr);
	sp_m = sp_l + ds_ratio*(sp_r-sp_l);
	Node[INl][j].X = Spline(sp_m, BndWestSpline);
      } /* endfor */
    } /* endif */

    if (BndEastSpline.np != 0) {
      for (j  = JNl+1 ; j < JNu ; j += 2 ) {
	sp_l = getS(Node[INu][j-1].X, 
		    BndEastSpline);
	sp_r = getS(Node[INu][j+1].X, 
		    BndEastSpline);
	dl = abs(Node[INu][j  ].X - 
		 Node[INu][j-1].X);
	dr = abs(Node[INu][j+1].X - 
		 Node[INu][j  ].X);
	ds_ratio = dl/(dl+dr);
	sp_m = sp_l + ds_ratio*(sp_r-sp_l);
	Node[INu][j].X = Spline(sp_m, BndEastSpline);
      } /* endfor */
    } /* endif */

    if (BndSouthSpline.np != 0) {
      for ( i = INl+1 ; i < INu ; i += 2 ) {
	sp_l = getS(Node[i-1][JNl].X, 
		    BndSouthSpline);
	sp_r = getS(Node[i+1][JNl].X, 
		    BndSouthSpline);
	dl = abs(Node[i  ][JNl].X - 
		 Node[i-1][JNl].X);
	dr = abs(Node[i+1][JNl].X - 
		 Node[i  ][JNl].X);
	ds_ratio = dl/(dl+dr);
	sp_m = sp_l + ds_ratio*(sp_r-sp_l);
	Node[i][JNl].X = Spline(sp_m, BndSouthSpline);
      } /* endfor */
    } /* endif */

    if (BndNorthSpline.np != 0) {
      for ( i = INl+1 ; i < INu ; i += 2 ) {
	sp_l = getS(Node[i-1][JNu].X, 
		    BndNorthSpline);
	sp_r = getS(Node[i+1][JNu].X, 
		    BndNorthSpline);
	dl = abs(Node[i  ][JNu].X - 
		 Node[i-1][JNu].X);
	dr = abs(Node[i+1][JNu].X - 
		 Node[i  ][JNu].X);
	ds_ratio = dl/(dl+dr);
	sp_m = sp_l + ds_ratio*(sp_r-sp_l);
	Node[i][JNu].X = Spline(sp_m, BndNorthSpline);
      } /* endfor */
    } /* endif */

    /* Require update of the interior cells geometric properties. */
    Schedule_Interior_Mesh_Update();

    /* Set the boundary condition types for refined 
       quadrilateral mesh block. */

    Set_BCs();

    /* Compute the exterior nodes for refined 
       quadrilateral mesh block. */

    Update_Exterior_Nodes();

    /* Compute the cells for refined 
       quadrilateral mesh block. */

    Update_Cells();

  } /* endif */

}

/*!
 * Returns a new quadrilateral mesh block resulting 
 * from the coarsening of four original grid blocks 
 * with half the original mesh resolution.          
 */
void Grid2D_Quad_Block_HO::Coarsen_Mesh(const Grid2D_Quad_Block_HO &Grid_Original_SW,
					const Grid2D_Quad_Block_HO &Grid_Original_SE,
					const Grid2D_Quad_Block_HO &Grid_Original_NW,
					const Grid2D_Quad_Block_HO &Grid_Original_NE) {

  int i, j, i_coarse, j_coarse, mesh_coarsening_permitted;
 
  /* Allocate memory for the cells and nodes for the 
     coarsened quadrilateral mesh block. */

  if ( (Grid_Original_SW.NCi-2*Grid_Original_SW.Nghost-
	2*((Grid_Original_SW.NCi-2*Grid_Original_SW.Nghost)/2) != 0) || 
       (Grid_Original_SW.NCj-2*Grid_Original_SW.Nghost-
	2*((Grid_Original_SW.NCj-2*Grid_Original_SW.Nghost)/2) != 0) ||
       (Grid_Original_SW.NCi-2*Grid_Original_SW.Nghost < 2*Grid_Original_SW.Nghost) ||
       (Grid_Original_SW.NCj-2*Grid_Original_SW.Nghost < 2*Grid_Original_SW.Nghost) ||
       (Grid_Original_SE.NCi != Grid_Original_SW.NCi) ||
       (Grid_Original_SE.NCj != Grid_Original_SW.NCj) ||
       (Grid_Original_NW.NCi != Grid_Original_SW.NCi) ||
       (Grid_Original_NW.NCj != Grid_Original_SW.NCj) ||
       (Grid_Original_NE.NCi != Grid_Original_SW.NCi) ||
       (Grid_Original_NE.NCj != Grid_Original_SW.NCj) ||
       (Grid_Original_SW.Node == NULL) ||
       (Grid_Original_SE.Node == NULL) ||
       (Grid_Original_NW.Node == NULL) ||
       (Grid_Original_NE.Node == NULL) ) {
    mesh_coarsening_permitted = 0;
  } else {
    mesh_coarsening_permitted = 1;
    allocate(Grid_Original_SW.NCi-2*Grid_Original_SW.Nghost, 
	     Grid_Original_SW.NCj-2*Grid_Original_SW.Nghost,
	     Grid_Original_SW.Nghost,
	     max(Grid_Original_SW.MaxRecOrder(),
		 max(Grid_Original_SE.MaxRecOrder(),
		     max(Grid_Original_NW.MaxRecOrder(), Grid_Original_NE.MaxRecOrder() ))) );
  }/* endif */

  /* Copy boundary spline info for the coarsened
     quadrilateral mesh block. */

  if (mesh_coarsening_permitted) {

    // ====  North Spline ====
    if (Grid_Original_NW.BndNorthSpline.np != 0 &&
	Grid_Original_NE.BndNorthSpline.np != 0) {
      BndNorthSpline = Grid_Original_NW.BndNorthSpline;
    } /* endif */

    // Set the West extension of the North spline
    if (Grid_Original_NW.ExtendWest_BndNorthSpline.np != 0){
      ExtendWest_BndNorthSpline = Grid_Original_NW.ExtendWest_BndNorthSpline;
    } else if (ExtendWest_BndNorthSpline.np != 0){
      ExtendWest_BndNorthSpline.deallocate();
    }

    // Set the East extension of the North spline
    if (Grid_Original_NE.ExtendEast_BndNorthSpline.np != 0){
      ExtendEast_BndNorthSpline = Grid_Original_NE.ExtendEast_BndNorthSpline;
    } else if (ExtendEast_BndNorthSpline.np != 0){
      ExtendEast_BndNorthSpline.deallocate();
    }

    // ====  South Spline ====
    if (Grid_Original_SW.BndSouthSpline.np != 0 &&
	Grid_Original_SE.BndSouthSpline.np != 0) {
      BndSouthSpline = Grid_Original_SW.BndSouthSpline;
    } /* endif */

    // Set the West extension of the South spline
    if (Grid_Original_SW.ExtendWest_BndSouthSpline.np != 0){
      ExtendWest_BndSouthSpline = Grid_Original_SW.ExtendWest_BndSouthSpline;
    } else if (ExtendWest_BndSouthSpline.np != 0){
      ExtendWest_BndSouthSpline.deallocate();
    }

    // Set the East extension of the South spline
    if (Grid_Original_SE.ExtendEast_BndSouthSpline.np != 0){
      ExtendEast_BndSouthSpline = Grid_Original_SE.ExtendEast_BndSouthSpline;
    } else if (ExtendEast_BndSouthSpline.np != 0){
      ExtendEast_BndSouthSpline.deallocate();
    }

    // ====  East Spline ====
    if (Grid_Original_SE.BndEastSpline.np != 0 &&
	Grid_Original_NE.BndEastSpline.np != 0) {
      BndEastSpline = Grid_Original_SE.BndEastSpline;
    } /* endif */

    // Set the South extension of the East spline
    if (Grid_Original_SE.ExtendSouth_BndEastSpline.np != 0){
      ExtendSouth_BndEastSpline = Grid_Original_SE.ExtendSouth_BndEastSpline;
    } else if (ExtendSouth_BndEastSpline.np != 0){
      ExtendSouth_BndEastSpline.deallocate();
    }

    // Set the North extension of the East spline
    if (Grid_Original_NE.ExtendNorth_BndEastSpline.np != 0){
      ExtendNorth_BndEastSpline = Grid_Original_NE.ExtendNorth_BndEastSpline;
    } else if (ExtendNorth_BndEastSpline.np != 0){
      ExtendNorth_BndEastSpline.deallocate();
    }
 
    // ====  West Spline ====
    if (Grid_Original_SW.BndWestSpline.np != 0 &&
	Grid_Original_NW.BndWestSpline.np != 0) {
      BndWestSpline = Grid_Original_SW.BndWestSpline;
    } /* endif */

    // Set the South extension of the West spline
    if (Grid_Original_SW.ExtendSouth_BndWestSpline.np != 0){
      ExtendSouth_BndWestSpline = Grid_Original_SW.ExtendSouth_BndWestSpline;
    } else if (ExtendSouth_BndWestSpline.np != 0){
      ExtendSouth_BndWestSpline.deallocate();
    }

    // Set the North extension of the West spline
    if (Grid_Original_NW.ExtendNorth_BndWestSpline.np != 0){
      ExtendNorth_BndWestSpline = Grid_Original_NW.ExtendNorth_BndWestSpline;
    } else if (ExtendNorth_BndWestSpline.np != 0){
      ExtendNorth_BndWestSpline.deallocate();
    }

    /* Assign boundary spline pathlength info for the coarsened
       quadrilateral mesh block. */

    if (Grid_Original_NW.BndNorthSpline.np != 0 &&
	Grid_Original_NE.BndNorthSpline.np != 0) {
      SminN = Grid_Original_NW.SminN;
      SmaxN = Grid_Original_NE.SmaxN;
    } else {
      SminN = ZERO;
      SmaxN = ZERO;
    } /* endif */

    if (Grid_Original_SW.BndSouthSpline.np != 0 &&
	Grid_Original_SE.BndSouthSpline.np != 0) {
      SminS = Grid_Original_SW.SminS;
      SmaxS = Grid_Original_SE.SmaxS;
    } else {
      SminS = ZERO;
      SmaxS = ZERO;
    } /* endif */

    if (Grid_Original_SE.BndEastSpline.np != 0 &&
	Grid_Original_NE.BndEastSpline.np != 0) {
      SminE = Grid_Original_SE.SmaxE;
      SmaxE = Grid_Original_NE.SmaxE;
    } else {
      SminE = ZERO;
      SmaxE = ZERO;
    } /* endif */

    if (Grid_Original_SW.BndWestSpline.np != 0 &&
	Grid_Original_NW.BndWestSpline.np != 0) {
      SminW = Grid_Original_SW.SminW;
      SmaxW = Grid_Original_NW.SmaxW;
    } else {
      SminW = ZERO;
      SmaxW = ZERO;
    } /* endif */

    /* Copy node stretching info to coarsened quadrilateral 
       mesh block. */

    StretchI = Grid_Original_SW.StretchI;
    BetaI = Grid_Original_SW.BetaI;
    TauI = Grid_Original_SW.TauI;
    StretchJ = Grid_Original_SW.StretchJ;
    BetaJ = Grid_Original_SW.BetaJ;
    TauJ = Grid_Original_SW.TauJ;
    if (Grid_Original_NW.BndNorthSpline.np != 0 &&
	Grid_Original_NE.BndNorthSpline.np != 0) {
      OrthogonalN = Grid_Original_NW.OrthogonalN;
    } else {
      OrthogonalN = ORTHOGONAL;
    } /* endif */
    if (Grid_Original_SW.BndSouthSpline.np != 0 &&
	Grid_Original_SE.BndSouthSpline.np != 0) {
      OrthogonalS = Grid_Original_SW.OrthogonalS;
    } else {
      OrthogonalS = ORTHOGONAL;
    } /* endif */
    if (Grid_Original_SE.BndEastSpline.np != 0 &&
	Grid_Original_NE.BndEastSpline.np != 0) {
      OrthogonalE = Grid_Original_SE.OrthogonalE;
    } else {
      OrthogonalE = ORTHOGONAL;
    } /* endif */
    if (Grid_Original_SW.BndWestSpline.np != 0 &&
	Grid_Original_NW.BndWestSpline.np != 0) {
      OrthogonalW = Grid_Original_SW.OrthogonalW;
    } else {
      OrthogonalW = ORTHOGONAL;
    } /* endif */

    // Force orthogonality at all coarsened mesh boundaries.
    OrthogonalN = ORTHOGONAL;
    OrthogonalS = ORTHOGONAL;
    OrthogonalE = ORTHOGONAL;
    OrthogonalW = ORTHOGONAL;

    /* Determine the node locations for the coarsened
       quadrilateral mesh block. */

    for ( j = Grid_Original_SW.JNl; j <= Grid_Original_SW.JNu ; j += 2 ) {
      for ( i = Grid_Original_SW.INl ; i <= Grid_Original_SW.INu ; i += 2 ) {
	i_coarse = (i-Grid_Original_SW.INl)/2+
	  INl;
	j_coarse = (j-Grid_Original_SW.JNl)/2+
	  JNl;
	Node[i_coarse][j_coarse].X = Grid_Original_SW.Node[i][j].X;
      } /* endfor */
    } /* endfor */

    for ( j = Grid_Original_SE.JNl; j <= Grid_Original_SE.JNu ; j += 2 ) {
      for ( i = Grid_Original_SE.INl ; i <= Grid_Original_SE.INu ; i += 2 ) {
	i_coarse = (i-Grid_Original_SE.INl)/2+
	  (INu-INl)/2+INl;
	j_coarse = (j-Grid_Original_SE.JNl)/2+
	  JNl;
	Node[i_coarse][j_coarse].X = Grid_Original_SE.Node[i][j].X;
      } /* endfor */
    } /* endfor */

    for ( j = Grid_Original_NW.JNl; j <= Grid_Original_NW.JNu ; j += 2 ) {
      for ( i = Grid_Original_NW.INl ; i <= Grid_Original_NW.INu ; i += 2 ) {
	i_coarse = (i-Grid_Original_NW.INl)/2+
	  INl;
	j_coarse = (j-Grid_Original_NW.JNl)/2+
	  (JNu-JNl)/2+JNl;
	Node[i_coarse][j_coarse].X = Grid_Original_NW.Node[i][j].X;
      } /* endfor */
    } /* endfor */

    for ( j = Grid_Original_NE.JNl; j <= Grid_Original_NE.JNu ; j += 2 ) {
      for ( i = Grid_Original_NE.INl ; i <= Grid_Original_NE.INu ; i += 2 ) {
	i_coarse = (i-Grid_Original_NE.INl)/2+
	  (INu-INl)/2+INl;
	j_coarse = (j-Grid_Original_NE.JNl)/2+
	  (JNu-JNl)/2+JNl;
	Node[i_coarse][j_coarse].X = Grid_Original_NE.Node[i][j].X;
      } /* endfor */
    } /* endfor */

    /* Require update of the interior cells geometric properties
       for newly coarsed quadrilateral mesh block. */
    Schedule_Interior_Mesh_Update();

    /* Set the boundary condition types for newly coarsened
       quadrilateral mesh block. */

    Set_BCs();

    /* Compute the exterior nodes for newly coarsened 
       quadrilateral mesh block. */

    Update_Exterior_Nodes();

    /* Compute the cells for newly coarsened
       quadrilateral mesh block. */

    Update_Cells();

  } /* endif */

}


/*!
 * Adjusts the locations of the boundary nodes of a     
 * quadrilateral grid block so that the new node        
 * locations will match with cell volumes of adjacent   
 * quadrilateral mesh blocks that have lower levels of  
 * mesh refinement (i.e., are coarser mesh blocks).     
 */
void Grid2D_Quad_Block_HO::Fix_Refined_Mesh_Boundaries(const int Fix_North_Boundary,
						       const int Fix_South_Boundary,
						       const int Fix_East_Boundary,
						       const int Fix_West_Boundary) {

  int i, j;
  double ds_ratio, dl, dr;
 
  /* Adjust the node locations of the north boundary
     of the quadrilateral mesh block. */

  if (Fix_North_Boundary) {
    for ( i = INl+1 ; i <= INu-1 ; i+=2 ) {
      dl = abs(Node[i  ][JNu].X - Node[i-1][JNu].X);
      dr = abs(Node[i+1][JNu].X - Node[i  ][JNu].X);
      ds_ratio = dl/(dl+dr);
      Node[i][JNu].X = 	Node[i-1][JNu].X + ds_ratio*(Node[i+1][JNu].X - Node[i-1][JNu].X);
    } /* endfor */
  } /* endif */

    /* Adjust the node locations of the south boundary
       of the quadrilateral mesh block. */

  if (Fix_South_Boundary) {
    for ( i = INl+1 ; i <= INu-1 ; i+=2 ) {
      dl = abs(Node[i  ][JNl].X - Node[i-1][JNl].X);
      dr = abs(Node[i+1][JNl].X - Node[i  ][JNl].X);
      ds_ratio = dl/(dl+dr);
      Node[i][JNl].X =  Node[i-1][JNl].X + ds_ratio*(Node[i+1][JNl].X - Node[i-1][JNl].X);
    } /* endfor */
  } /* endif */

    /* Adjust the node locations of the east boundary
       of the quadrilateral mesh block. */

  if (Fix_East_Boundary) {
    for ( j  = JNl+1; j <= JNu-1; j+=2 ) {
      dl = abs(Node[INu][j  ].X - Node[INu][j-1].X);
      dr = abs(Node[INu][j+1].X - Node[INu][j  ].X);
      ds_ratio = dl/(dl+dr);
      Node[INu][j].X = 	Node[INu][j-1].X + ds_ratio*(Node[INu][j+1].X - Node[INu][j-1].X);
    } /* endfor */
  } /* endif */

    /* Adjust the node locations of the west boundary
       of the quadrilateral mesh block. */

  if (Fix_West_Boundary) {
    for ( j  = JNl+1; j <= JNu-1; j+=2 ) {
      dl = abs(Node[INl][j  ].X - Node[INl][j-1].X);
      dr = abs(Node[INl][j+1].X - Node[INl][j  ].X);
      ds_ratio = dl/(dl+dr);
      Node[INl][j].X = 	Node[INl][j-1].X + ds_ratio*(Node[INl][j+1].X - Node[INl][j-1].X);
    } /* endfor */
  }/* endif */
  
  /* Require update of the interior cells geometric properties. */
  Schedule_Interior_Mesh_Update();

  /* Reset the boundary condition types at the grid block
     boundaries. */
 
  Set_BCs();

  /* Recompute the exterior nodes for the quadrilateral mesh block. */

  Update_Exterior_Nodes();

  /* Recompute the cells for the quadrilateral mesh block. */

  Update_Cells();

}


/*!
 * Returns the adjusted the locations of the boundary   
 * nodes of a quadrilateral grid block to their         
 * original unmodified positions.                       
 *                                                      
 */
void Grid2D_Quad_Block_HO::Unfix_Refined_Mesh_Boundaries(void) {

  int i, j;
  double sp_l, sp_r, sp_m, ds_ratio, dl, dr;
 
  /* Return the nodes at the north boundary
     to their original positions. */

  if (BndNorthSpline.np != 0) {
    for ( i = INl+1 ; i < INu ; i += 2 ) {
      sp_l = getS(Node[i-1][JNu].X, 
		  BndNorthSpline);
      sp_r = getS(Node[i+1][JNu].X, 
		  BndNorthSpline);
      dl = abs(Node[i  ][JNu].X - 
	       Node[i-1][JNu].X);
      dr = abs(Node[i+1][JNu].X - 
	       Node[i  ][JNu].X);
      ds_ratio = dl/(dl+dr);
      sp_m = sp_l + ds_ratio*(sp_r-sp_l);
      Node[i][JNu].X = Spline(sp_m, BndNorthSpline);
    } /* endfor */
  } /* endif */

    /* Return the nodes at the south boundary
       to their original positions. */

  if (BndSouthSpline.np != 0) {
    for ( i = INl+1 ; i < INu ; i += 2 ) {
      sp_l = getS(Node[i-1][JNl].X, 
		  BndSouthSpline);
      sp_r = getS(Node[i+1][JNl].X, 
		  BndSouthSpline);
      dl = abs(Node[i  ][JNl].X - 
	       Node[i-1][JNl].X);
      dr = abs(Node[i+1][JNl].X - 
	       Node[i  ][JNl].X);
      ds_ratio = dl/(dl+dr);
      sp_m = sp_l + ds_ratio*(sp_r-sp_l);
      Node[i][JNl].X = Spline(sp_m, BndSouthSpline);
    } /* endfor */
  } /* endif */

    /* Return the nodes at the east boundary
       to their original positions. */

  if (BndEastSpline.np != 0) {
    for (j  = JNl+1 ; j < JNu ; j += 2 ) {
      sp_l = getS(Node[INu][j-1].X, 
		  BndEastSpline);
      sp_r = getS(Node[INu][j+1].X, 
		  BndEastSpline);
      dl = abs(Node[INu][j  ].X - 
	       Node[INu][j-1].X);
      dr = abs(Node[INu][j+1].X - 
	       Node[INu][j  ].X);
      ds_ratio = dl/(dl+dr);
      sp_m = sp_l + ds_ratio*(sp_r-sp_l);
      Node[INu][j].X = Spline(sp_m, BndEastSpline);
    } /* endfor */
  } /* endif */

    /* Return the nodes at the west boundary
       to their original positions. */

  if (BndWestSpline.np != 0) {
    for (j  = JNl+1 ; j < JNu ; j += 2 ) {
      sp_l = getS(Node[INl][j-1].X, 
		  BndWestSpline);
      sp_r = getS(Node[INl][j+1].X, 
		  BndWestSpline);
      dl = abs(Node[INl][j  ].X - 
	       Node[INl][j-1].X);
      dr = abs(Node[INl][j+1].X - 
	       Node[INl][j  ].X);
      ds_ratio = dl/(dl+dr);
      sp_m = sp_l + ds_ratio*(sp_r-sp_l);
      Node[INl][j].X = Spline(sp_m, BndWestSpline);
    } /* endfor */
  }/* endif */

  /* Require update of the interior cells geometric properties. */
  Schedule_Interior_Mesh_Update();

  /* Reset the boundary condition types at the grid block
     boundaries. */
 
  Set_BCs();

  /* Recompute the exterior nodes for the quadrilateral mesh block. */

  Update_Exterior_Nodes();

  /* Recompute the cells for the quadrilateral mesh block. */

  Update_Cells();

}

/*!
 * Output operator.
 */
ostream &operator << (ostream &out_file, 
		      const Grid2D_Quad_Block_HO &G) {
  int i, j;
  out_file << G.NNi << " " << G.INl << " " << G.INu << "\n";
  out_file << G.NNj << " " << G.JNl << " " << G.JNu << "\n";
  out_file << G.Nghost << "\n";
  out_file << G.MaxRecOrder() << "\n"; // used for setting GeomCoeff
  out_file << G.getHighOrderBoundaryValue() << "\n";
  if (G.NNi == 0 || G.NNj == 0) return(out_file);
  out_file << G.NCi << " " << G.ICl << " " << G.ICu << "\n";
  out_file << G.NCj << " " << G.JCl << " " << G.JCu << "\n";

  out_file.precision(15);
  // Output node data
  for ( j = G.JNl-G.Nghost ; j <= G.JNu+G.Nghost; ++j ) {
    for ( i = G.INl-G.Nghost ; i <= G.INu+G.Nghost; ++i ) {
      out_file << G.Node[i][j].X << "\n";
    } /* endfor */
  } /* endfor */
  // Output cell data
  for ( j = G.JCl-G.Nghost ; j <= G.JCu+G.Nghost; ++j ) {
    for ( i = G.ICl-G.Nghost ; i <= G.ICu+G.Nghost; ++i ) {
      out_file << G.Cell[i][j] << "\n";
    } /* endfor */
  } /* endfor */
  out_file.precision(15);
  for ( i = G.ICl-G.Nghost ; i <= G.ICu+G.Nghost ; ++i) {
    out_file << G.BCtypeN[i] << " " << G.BCtypeS[i] << "\n";
  } /* endfor */
  for ( j = G.JCl-G.Nghost ; j <= G.JCu+G.Nghost ; ++j) {
    out_file << G.BCtypeE[j] << " " << G.BCtypeW[j] << "\n";
  } /* endfor */

  // Output North boundary spline information
  if (G.BndNorthSpline.np != 0 ) {
    out_file << G.BndNorthSpline;
  } else {
    out_file << G.BndNorthSpline.np << "\n";
  } /* endif */
  if (G.BndNorthSplineInfo != NULL){
    out_file << G.NCi << "\n"; 	// number of SplineInterval2D elements
    // Output each active component of BndNorthSplineInfo (no ghost cells)
    for ( i = G.ICl; i <= G.ICu; ++i) {
      out_file << G.BndNorthSplineInfo[i] << "\n";
    }   
  } else {
    out_file << 0 << "\n"; // zero SplineInterval2D elements
  }

  // Output South boundary spline information
  if (G.BndSouthSpline.np != 0 ) {
    out_file << G.BndSouthSpline;
  } else {
    out_file << G.BndSouthSpline.np << "\n";
  } /* endif */
  if (G.BndSouthSplineInfo != NULL){
    out_file << G.NCi << "\n"; 	// number of SplineInterval2D elements
    // Output each active component of BndSouthSplineInfo (no ghost cells)
    for ( i = G.ICl; i <= G.ICu; ++i) {
      out_file << G.BndSouthSplineInfo[i] << "\n";
    }   
  } else {
    out_file << 0 << "\n"; // zero SplineInterval2D elements
  }

  // Output East boundary spline information
  if (G.BndEastSpline.np != 0 ) {
    out_file << G.BndEastSpline;
  } else {
    out_file << G.BndEastSpline.np << "\n";
  } /* endif */
  if (G.BndEastSplineInfo != NULL){
    out_file << G.NCj << "\n"; 	// number of SplineInterval2D elements
    // Output each active component of BndEastSplineInfo (no ghost cells)
    for ( i = G.JCl; i <= G.JCu; ++i) {
      out_file << G.BndEastSplineInfo[i] << "\n";
    }   
  } else {
    out_file << 0 << "\n"; // zero SplineInterval2D elements
  }

  // Output West boundary spline information
  if (G.BndWestSpline.np != 0 ) {
    out_file << G.BndWestSpline;
  } else {
    out_file << G.BndWestSpline.np << "\n";
  } /* endif */
  if (G.BndWestSplineInfo != NULL){
    out_file << G.NCj << "\n"; 	// number of SplineInterval2D elements
    // Output each active component of BndWestSplineInfo (no ghost cells)
    for ( i = G.JCl; i <= G.JCu; ++i) {
      out_file << G.BndWestSplineInfo[i] << "\n";
    }   
  } else {
    out_file << 0 << "\n"; // zero SplineInterval2D elements
  }

  // Output west extension of North boundary spline information
  if (G.ExtendWest_BndNorthSpline.np != 0 ) {
    out_file << G.ExtendWest_BndNorthSpline;
  } else {
    out_file << G.ExtendWest_BndNorthSpline.np << "\n";
  } /* endif */
  if (G.ExtendWest_BndNorthSplineInfo != NULL){
    out_file << G.Nghost << "\n"; 	// number of SplineInterval2D elements
    // Output each component of ExtendWest_BndNorthSplineInfo
    for ( i = 0; i < G.Nghost; ++i) {
      out_file << G.ExtendWest_BndNorthSplineInfo[i] << "\n";
    }
  } else {
    out_file << 0 << "\n"; // zero SplineInterval2D elements
  }

  // Output east extension of North boundary spline information
  if (G.ExtendEast_BndNorthSpline.np != 0 ) {
    out_file << G.ExtendEast_BndNorthSpline;
  } else {
    out_file << G.ExtendEast_BndNorthSpline.np << "\n";
  } /* endif */
  if (G.ExtendEast_BndNorthSplineInfo != NULL){
    out_file << G.Nghost << "\n"; 	// number of SplineInterval2D elements
    // Output each component of ExtendEast_BndNorthSplineInfo
    for ( i = 0; i < G.Nghost; ++i) {
      out_file << G.ExtendEast_BndNorthSplineInfo[i] << "\n";
    }
  } else {
    out_file << 0 << "\n"; // zero SplineInterval2D elements
  }

  // Output west extension of South boundary spline information
  if (G.ExtendWest_BndSouthSpline.np != 0 ) {
    out_file << G.ExtendWest_BndSouthSpline;
  } else {
    out_file << G.ExtendWest_BndSouthSpline.np << "\n";
  } /* endif */
  if (G.ExtendWest_BndSouthSplineInfo != NULL){
    out_file << G.Nghost << "\n"; 	// number of SplineInterval2D elements
    // Output each component of ExtendWest_BndSouthSplineInfo
    for ( i = 0; i < G.Nghost; ++i) {
      out_file << G.ExtendWest_BndSouthSplineInfo[i] << "\n";
    }
  } else {
    out_file << 0 << "\n"; // zero SplineInterval2D elements
  }

  // Output east extension of South boundary spline information
  if (G.ExtendEast_BndSouthSpline.np != 0 ) {
    out_file << G.ExtendEast_BndSouthSpline;
  } else {
    out_file << G.ExtendEast_BndSouthSpline.np << "\n";
  } /* endif */
  if (G.ExtendEast_BndSouthSplineInfo != NULL){
    out_file << G.Nghost << "\n"; 	// number of SplineInterval2D elements
    // Output each component of ExtendEast_BndSouthSplineInfo
    for ( i = 0; i < G.Nghost; ++i) {
      out_file << G.ExtendEast_BndSouthSplineInfo[i] << "\n";
    }
  } else {
    out_file << 0 << "\n"; // zero SplineInterval2D elements
  }

  // Output north extension of East boundary spline information
  if (G.ExtendNorth_BndEastSpline.np != 0 ) {
    out_file << G.ExtendNorth_BndEastSpline;
  } else {
    out_file << G.ExtendNorth_BndEastSpline.np << "\n";
  } /* endif */
  if (G.ExtendNorth_BndEastSplineInfo != NULL){
    out_file << G.Nghost << "\n"; 	// number of SplineInterval2D elements
    // Output each component of ExtendNorth_BndEastSplineInfo
    for ( i = 0; i < G.Nghost; ++i) {
      out_file << G.ExtendNorth_BndEastSplineInfo[i] << "\n";
    }
  } else {
    out_file << 0 << "\n"; // zero SplineInterval2D elements
  }

  // Output south extension of East boundary spline information
  if (G.ExtendSouth_BndEastSpline.np != 0 ) {
    out_file << G.ExtendSouth_BndEastSpline;
  } else {
    out_file << G.ExtendSouth_BndEastSpline.np << "\n";
  } /* endif */
  if (G.ExtendSouth_BndEastSplineInfo != NULL){
    out_file << G.Nghost << "\n"; 	// number of SplineInterval2D elements
    // Output each component of ExtendSouth_BndEastSplineInfo
    for ( i = 0; i < G.Nghost; ++i) {
      out_file << G.ExtendSouth_BndEastSplineInfo[i] << "\n";
    }
  } else {
    out_file << 0 << "\n"; // zero SplineInterval2D elements
  }

  // Output north extension of West boundary spline information
  if (G.ExtendNorth_BndWestSpline.np != 0 ) {
    out_file << G.ExtendNorth_BndWestSpline;
  } else {
    out_file << G.ExtendNorth_BndWestSpline.np << "\n";
  } /* endif */
  if (G.ExtendNorth_BndWestSplineInfo != NULL){
    out_file << G.Nghost << "\n"; 	// number of SplineInterval2D elements
    // Output each component of ExtendNorth_BndWestSplineInfo
    for ( i = 0; i < G.Nghost; ++i) {
      out_file << G.ExtendNorth_BndWestSplineInfo[i] << "\n";
    }
  } else {
    out_file << 0 << "\n"; // zero SplineInterval2D elements
  }

  // Output south extension of West boundary spline information
  if (G.ExtendSouth_BndWestSpline.np != 0 ) {
    out_file << G.ExtendSouth_BndWestSpline;
  } else {
    out_file << G.ExtendSouth_BndWestSpline.np << "\n";
  } /* endif */
  if (G.ExtendSouth_BndWestSplineInfo != NULL){
    out_file << G.Nghost << "\n"; 	// number of SplineInterval2D elements
    // Output each component of ExtendSouth_BndWestSplineInfo
    for ( i = 0; i < G.Nghost; ++i) {
      out_file << G.ExtendSouth_BndWestSplineInfo[i] << "\n";
    }
  } else {
    out_file << 0 << "\n"; // zero SplineInterval2D elements
  }

  out_file.setf(ios::scientific);
  out_file << G.SminN << " " << G.SmaxN << " " << G.SminS << " " << G.SmaxS << "\n"; 
  out_file << G.SminE << " " << G.SmaxE << " " << G.SminW << " " << G.SmaxW << "\n";
  out_file << G.StretchI << " " << G.BetaI << " " << G.TauI << "\n";
  out_file << G.StretchJ << " " << G.BetaJ << " " << G.TauJ << "\n";
  out_file.unsetf(ios::scientific);
  out_file << G.OrthogonalN << " "  << G.OrthogonalS << " "  
           << G.OrthogonalE << " "  << G.OrthogonalW << "\n";
  return (out_file);
}

/*!
 * Input operator.
 */
istream &operator >> (istream &in_file, 
		      Grid2D_Quad_Block_HO &G) {
  int i, j, ni, il, iu, nj, jl, ju, ng, RecOrder, BoundaryRepresentationFlag;

  // Read mesh parameters
  in_file.setf(ios::skipws);
  // Read indexes for nodes
  in_file >> ni >> il >> iu;
  in_file >> nj >> jl >> ju;
  in_file >> ng;
  in_file >> RecOrder;	// added for GeomCoeff
  in_file >> BoundaryRepresentationFlag;
  in_file.unsetf(ios::skipws);

  // Provide enough memory for the new mesh
  if (ni == 0 || nj == 0) {
    if (G.Node != NULL) G.deallocate(); return(in_file);
  } /* endif */
  G.allocate(ni-2*ng-1, nj-2*ng-1, ng, RecOrder);

  // Read indexes for cells
  in_file.setf(ios::skipws);
  in_file >> ni >> il >> iu; in_file >> nj >> jl >> ju;
  in_file.unsetf(ios::skipws);

  // Read the coordinates of the mesh nodes
  for ( j = G.JNl-G.Nghost ; j <= G.JNu+G.Nghost; ++j ) {
    for ( i = G.INl-G.Nghost ; i <= G.INu+G.Nghost; ++i ) {
      in_file >> G.Node[i][j].X;
    } /* endfor */
  } /* endfor */

  // Read the cell parameters
  for ( j = G.JCl-G.Nghost ; j <= G.JCu+G.Nghost; ++j ) {
    for ( i = G.ICl-G.Nghost ; i <= G.ICu+G.Nghost; ++i ) {
      in_file >> G.Cell[i][j];
    } /* endfor */
  } /* endfor */

  // Read the type of boundary conditions for each boundary
  for ( i = G.ICl-G.Nghost ; i <= G.ICu+G.Nghost ; ++i) {
    in_file.setf(ios::skipws);
    in_file >> G.BCtypeN[i] >> G.BCtypeS[i];
    in_file.unsetf(ios::skipws);
  } /* endfor */
  for ( j = G.JCl-G.Nghost ; j <= G.JCu+G.Nghost ; ++j) {
    in_file.setf(ios::skipws);
    in_file >> G.BCtypeE[j] >> G.BCtypeW[j];
    in_file.unsetf(ios::skipws);
  } /* endfor */

  // Read the North boundary spline
  in_file >> G.BndNorthSpline;
  // Read the North boundary spline info
  in_file.setf(ios::skipws);
  in_file >> ni;
  in_file.unsetf(ios::skipws);
  if (ni > 0){
    // allocate memory
    G.BndNorthSplineInfo = new Spline2DInterval_HO [ni];
    // Read spline info
    for(i=G.ICl; i<=G.ICu; ++i) {
      in_file >> G.BndNorthSplineInfo[i];
    }
  } else {
    G.deallocate_BndNorthSplineInfo();
  }

  // Read the South boundary spline
  in_file >> G.BndSouthSpline;
  // Read the South boundary spline info
  in_file.setf(ios::skipws);
  in_file >> ni;
  in_file.unsetf(ios::skipws);
  if (ni > 0){
    // allocate memory
    G.BndSouthSplineInfo = new Spline2DInterval_HO [ni];
    // Read spline info
    for(i=G.ICl; i<=G.ICu; ++i) {
      in_file >> G.BndSouthSplineInfo[i];
    }
  } else {
    G.deallocate_BndSouthSplineInfo();
  }

  // Read the East boundary spline
  in_file >> G.BndEastSpline;
  // Read the East boundary spline info
  in_file.setf(ios::skipws);
  in_file >> nj; 
  in_file.unsetf(ios::skipws);
  if (nj > 0){
    // allocate memory
    G.BndEastSplineInfo = new Spline2DInterval_HO [nj];
    // Read spline info
    for(j=G.JCl; j<=G.JCu; ++j) {
      in_file >> G.BndEastSplineInfo[j];
    }
  } else {
    G.deallocate_BndEastSplineInfo();
  }

  // Read the West boundary spline
  in_file >> G.BndWestSpline;
  // Read the West boundary spline info
  in_file.setf(ios::skipws);
  in_file >> nj; 
  in_file.unsetf(ios::skipws);
  if (nj > 0){
    // allocate memory
    G.BndWestSplineInfo = new Spline2DInterval_HO [nj];
    // Read spline info
    for(j=G.JCl; j<=G.JCu; ++j) {
      in_file >> G.BndWestSplineInfo[j];
    }
  } else {
    G.deallocate_BndWestSplineInfo();
  }

  // Read the west extention of North boundary spline
  in_file >> G.ExtendWest_BndNorthSpline;
  // Read the west extention of North boundary spline info
  in_file.setf(ios::skipws);
  in_file >> nj; 
  in_file.unsetf(ios::skipws);
  if (nj > 0){
    // allocate memory
    G.ExtendWest_BndNorthSplineInfo = new Spline2DInterval_HO [nj];
    // Read spline info
    for(j=0; j<nj; ++j) {
      in_file >> G.ExtendWest_BndNorthSplineInfo[j];
    }
  } else {
    G.deallocate_ExtendWest_BndNorthSplineInfo();
  }

  // Read the east extention of North boundary spline
  in_file >> G.ExtendEast_BndNorthSpline;
  // Read the east extention of North boundary spline info
  in_file.setf(ios::skipws);
  in_file >> nj; 
  in_file.unsetf(ios::skipws);
  if (nj > 0){
    // allocate memory
    G.ExtendEast_BndNorthSplineInfo = new Spline2DInterval_HO [nj];
    // Read spline info
    for(j=0; j<nj; ++j) {
      in_file >> G.ExtendEast_BndNorthSplineInfo[j];
    }
  } else {
    G.deallocate_ExtendEast_BndNorthSplineInfo();
  }

  // Read the west extention of South boundary spline
  in_file >> G.ExtendWest_BndSouthSpline;
  // Read the west extention of South boundary spline info
  in_file.setf(ios::skipws);
  in_file >> nj; 
  in_file.unsetf(ios::skipws);
  if (nj > 0){
    // allocate memory
    G.ExtendWest_BndSouthSplineInfo = new Spline2DInterval_HO [nj];
    // Read spline info
    for(j=0; j<nj; ++j) {
      in_file >> G.ExtendWest_BndSouthSplineInfo[j];
    }
  } else {
    G.deallocate_ExtendWest_BndSouthSplineInfo();
  }

  // Read the east extention of South boundary spline
  in_file >> G.ExtendEast_BndSouthSpline;
  // Read the east extention of South boundary spline info
  in_file.setf(ios::skipws);
  in_file >> nj; 
  in_file.unsetf(ios::skipws);
  if (nj > 0){
    // allocate memory
    G.ExtendEast_BndSouthSplineInfo = new Spline2DInterval_HO [nj];
    // Read spline info
    for(j=0; j<nj; ++j) {
      in_file >> G.ExtendEast_BndSouthSplineInfo[j];
    }
  } else {
    G.deallocate_ExtendEast_BndSouthSplineInfo();
  }

  // Read the north extention of East boundary spline
  in_file >> G.ExtendNorth_BndEastSpline;
  // Read the north extention of East boundary spline info
  in_file.setf(ios::skipws);
  in_file >> nj; 
  in_file.unsetf(ios::skipws);
  if (nj > 0){
    // allocate memory
    G.ExtendNorth_BndEastSplineInfo = new Spline2DInterval_HO [nj];
    // Read spline info
    for(j=0; j<nj; ++j) {
      in_file >> G.ExtendNorth_BndEastSplineInfo[j];
    }
  } else {
    G.deallocate_ExtendNorth_BndEastSplineInfo();
  }

  // Read the south extention of East boundary spline
  in_file >> G.ExtendSouth_BndEastSpline;
  // Read the south extention of East boundary spline info
  in_file.setf(ios::skipws);
  in_file >> nj; 
  in_file.unsetf(ios::skipws);
  if (nj > 0){
    // allocate memory
    G.ExtendSouth_BndEastSplineInfo = new Spline2DInterval_HO [nj];
    // Read spline info
    for(j=0; j<nj; ++j) {
      in_file >> G.ExtendSouth_BndEastSplineInfo[j];
    }
  } else {
    G.deallocate_ExtendSouth_BndEastSplineInfo();
  }

  // Read the north extention of West boundary spline
  in_file >> G.ExtendNorth_BndWestSpline;
  // Read the north extention of West boundary spline info
  in_file.setf(ios::skipws);
  in_file >> nj; 
  in_file.unsetf(ios::skipws);
  if (nj > 0){
    // allocate memory
    G.ExtendNorth_BndWestSplineInfo = new Spline2DInterval_HO [nj];
    // Read spline info
    for(j=0; j<nj; ++j) {
      in_file >> G.ExtendNorth_BndWestSplineInfo[j];
    }
  } else {
    G.deallocate_ExtendNorth_BndWestSplineInfo();
  }

  // Read the south extention of West boundary spline
  in_file >> G.ExtendSouth_BndWestSpline;
  // Read the south extention of West boundary spline info
  in_file.setf(ios::skipws);
  in_file >> nj; 
  in_file.unsetf(ios::skipws);
  if (nj > 0){
    // allocate memory
    G.ExtendSouth_BndWestSplineInfo = new Spline2DInterval_HO [nj];
    // Read spline info
    for(j=0; j<nj; ++j) {
      in_file >> G.ExtendSouth_BndWestSplineInfo[j];
    }
  } else {
    G.deallocate_ExtendSouth_BndWestSplineInfo();
  }

  in_file.setf(ios::skipws);
  in_file >> G.SminN >> G.SmaxN >> G.SminS >> G.SmaxS; 
  in_file >> G.SminE >> G.SmaxE >> G.SminW >> G.SmaxW;
  in_file >> G.StretchI >> G.BetaI >> G.TauI;
  in_file >> G.StretchJ >> G.BetaJ >> G.TauJ;
  in_file >> G.OrthogonalN >> G.OrthogonalS 
          >> G.OrthogonalE >> G.OrthogonalW;
  in_file.unsetf(ios::skipws);

  if (BoundaryRepresentationFlag == ON) {
    // set high-order geometric boundaries
    G.setHighOrderBoundaryRepresentation();
  } else {
    // set low-order geometric boundaries
    G.setLowOrderBoundaryRepresentation();
  }//endif

  // No geometric properties update is required because everything was read.
  G.Confirm_Mesh_Update_Everywhere();

  // Mark the current geometry different than the previous one
  G.New_Global_Geometry_State();

  return (in_file);
}

/*!
 * Determine the minimum distance between the Node[i][j] 
 * and all the neighbour edges and cell diagonals,
 * for a quadrilateral mesh. 
 * \verbatim
 *
 * For the quadrilateral mesh there are 8 edges and 4 diagonals 
 * that must be checked for determining the min distance.
 *
 *          6      5       
 *       o-----o------o    
 *       |            |    
 *     7 |   (i,j)    | 4  
 *       o     o      o    
 *       |            |    
 *     8 |            | 3  
 *       o-----o------o    
 *          1     2        
 *
 * \endverbatim
 *
 * \param [in] i i-index of the node
 * \param [in] j j-index of the node
 */
double Grid2D_Quad_Block_HO::MinimumNodeEdgeDistance(const int &i, const int &j) const {

  double MinDistance;

  // Edge 1
  MinDistance = DistanceFromPointToLine(Node[i][j].X, Node[i-1][j-1].X, Node[i][j-1].X);
  // Edge 2
  MinDistance = min(MinDistance, DistanceFromPointToLine(Node[i][j].X, Node[i][j-1].X, Node[i+1][j-1].X));
  // Edge 3
  MinDistance = min(MinDistance, DistanceFromPointToLine(Node[i][j].X, Node[i+1][j-1].X, Node[i+1][j].X));
  // Edge 4
  MinDistance = min(MinDistance, DistanceFromPointToLine(Node[i][j].X, Node[i+1][j].X, Node[i+1][j+1].X));
  // Edge 5
  MinDistance = min(MinDistance, DistanceFromPointToLine(Node[i][j].X, Node[i+1][j+1].X, Node[i][j+1].X));
  // Edge 6
  MinDistance = min(MinDistance, DistanceFromPointToLine(Node[i][j].X, Node[i][j+1].X, Node[i-1][j+1].X));
  // Edge 7
  MinDistance = min(MinDistance, DistanceFromPointToLine(Node[i][j].X, Node[i-1][j+1].X, Node[i-1][j].X));
  // Edge 8
  MinDistance = min(MinDistance, DistanceFromPointToLine(Node[i][j].X, Node[i-1][j].X, Node[i-1][j-1].X));

  // Diagonal 1
  MinDistance = min(MinDistance, DistanceFromPointToLine(Node[i][j].X, Node[i][j-1].X, Node[i+1][j].X));
  // Diagonal 2
  MinDistance = min(MinDistance, DistanceFromPointToLine(Node[i][j].X, Node[i+1][j].X, Node[i][j+1].X));
  // Diagonal 3
  MinDistance = min(MinDistance, DistanceFromPointToLine(Node[i][j].X, Node[i][j+1].X, Node[i-1][j].X));
  // Diagonal 4
  MinDistance = min(MinDistance, DistanceFromPointToLine(Node[i][j].X, Node[i-1][j].X, Node[i][j-1].X));

  return MinDistance;
}

/*!
 * Determines the distance between the point of interest
 * and a line defined by the two points. 
 *
 * \param Point the position vector of the point of interest
 * \param LinePoint1 the position vector of one line end 
 * \param LinePoint2 the position vector of the other line end
 */
double Grid2D_Quad_Block_HO::DistanceFromPointToLine(const Vector2D &Point,
						     const Vector2D &LinePoint1,
						     const Vector2D &LinePoint2) const {

  double DeltaX, DeltaY, X_Dis, Y_Dis;

  // X_Dis, Y_Dis -> the X and Y coordinates of the point which is used to measure the distance

  DeltaX = LinePoint2.x - LinePoint1.x;
  DeltaY = LinePoint2.y - LinePoint1.y;

  // Determine X_Dis
  X_Dis = DeltaX*DeltaX*Point.x + DeltaX*DeltaY*(Point.y - LinePoint1.y) + DeltaY*DeltaY*LinePoint1.x;
  X_Dis /= (DeltaX*DeltaX + DeltaY*DeltaY);

  if ( (fabs(X_Dis - LinePoint1.x) <= 1.0e-14) && (fabs(X_Dis - LinePoint2.x) <= 1.0e-14) ){
    // the line defined by the two points is parallel to Oy
    return fabs(Point.x - LinePoint1.x);
  }

  // Determine Y_Dis
  Y_Dis = LinePoint1.y + DeltaY*(X_Dis - LinePoint1.x)/DeltaX;

  // Return the distance value
  DeltaX = Point.x - X_Dis;
  DeltaY = Point.y - Y_Dis;

  return sqrt( DeltaX*DeltaX + DeltaY*DeltaY );
}

/*!
 * Output the cell data to the specified output stream.
 */
void Grid2D_Quad_Block_HO::Output_Cells_Data(const int &block_number, ostream &out_file){

  int i,j;

  out_file << "Block " << block_number << "\n";

  out_file << NCi << " " << ICl << " " << ICu << "\n";
  out_file << NCj << " " << JCl << " " << JCu << "\n";

  out_file.precision(15);
  // Output cell data
  if (Cell != NULL){
    for ( j = JCl-Nghost ; j <= JCu+Nghost; ++j ) {
      for ( i = ICl-Nghost ; i <= ICu+Nghost; ++i ) {
	out_file << Cell[i][j] << "\n";
      } /* endfor */
    } /* endfor */
  } else {
    out_file << 0 << "\n";
  }
  out_file.precision(15);
}

/*!
 * Output the geometric properties that are invariant
 * to translation and rotation of all block cells 
 * to the specified output stream.
 */
void Grid2D_Quad_Block_HO::Output_Cells_Translation_Rotation_Invariant_Properties(const int &block_number,
										  ostream &out_file){

  int i,j;

  out_file << "Block " << block_number << "\n";

  out_file << NCi << " " << ICl << " " << ICu << "\n";
  out_file << NCj << " " << JCl << " " << JCu << "\n";

  out_file.precision(15);
  // Output cell data
  if (Cell != NULL){
    for ( j = JCl-Nghost ; j <= JCu+Nghost; ++j ) {
      for ( i = ICl-Nghost ; i <= ICu+Nghost; ++i ) {
	Cell[i][j].Output_Translation_Rotation_Invariant_Properties(out_file);
      } /* endfor */
    } /* endfor */
  } else {
    out_file << 0 << "\n";
  }
  out_file.precision(15);
}

/*!
 * Determine the necessary geometric properties to
 * carry out Gauss quadrature integration along 
 * a segment line specified by the StartPoint and EndPoint.
 *
 * \param GQPoints vector of Gauss integration points
 * \param NumberOfGQPs the required number of integration points
 * \param dYds the derivative of the y-coordinate with respect to the pathlength coordinate (constant for this case)
 * \param Length the length of the line segment
 */
void Grid2D_Quad_Block_HO::getLineSegmentGaussIntegrationData(const Vector2D & StartPoint, const Vector2D & EndPoint,
							      Vector2D * GQPoints, const int & NumberOfGQPs,
							      double & DeltaY) const {

  Vector2D Temp;

  switch (NumberOfGQPs){
  case 1:
    GQPoints[0] = HALF*(StartPoint+EndPoint);
    break;

  case 2:
    Temp = StartPoint-EndPoint;
    
    /* final value GQPoints[0] */
    GQPoints[0] = EndPoint + GaussQuadratureData::GQ2_Abscissa[0]*Temp;
    
    /* final value GQPoints[1] */
    GQPoints[1] = EndPoint + GaussQuadratureData::GQ2_Abscissa[1]*Temp;
    break;

  case 3:
    Temp = StartPoint-EndPoint;
    
    /* final value GQPoints[0] */
    GQPoints[0] = EndPoint + GaussQuadratureData::GQ3_Abscissa[0]*Temp;
   
    /* final value GQPoints[1] */
    GQPoints[1] = EndPoint + GaussQuadratureData::GQ3_Abscissa[1]*Temp;

    /* final value GQPoints[2] */
    GQPoints[2] = EndPoint + GaussQuadratureData::GQ3_Abscissa[2]*Temp;
    break;

  case 5:
    Temp = StartPoint-EndPoint;
    
    /* final value GQPoints[0] */
    GQPoints[0] = EndPoint + GaussQuadratureData::GQ5_Abscissa[0]*Temp;
   
    /* final value GQPoints[1] */
    GQPoints[1] = EndPoint + GaussQuadratureData::GQ5_Abscissa[1]*Temp;

    /* final value GQPoints[2] */
    GQPoints[2] = EndPoint + GaussQuadratureData::GQ5_Abscissa[2]*Temp;

    /* final value GQPoints[3] */
    GQPoints[3] = EndPoint + GaussQuadratureData::GQ5_Abscissa[3]*Temp;

    /* final value GQPoints[4] */
    GQPoints[4] = EndPoint + GaussQuadratureData::GQ5_Abscissa[4]*Temp;
    break;

  default:
    throw runtime_error("Grid2D_Quad_Block_HO::getLineSegmentGaussIntegrationData() ERROR! \
                         Not implemented number of Gauss quadrature points!");
  }
  
  // Calculate DeltaY
  DeltaY = (EndPoint.y - StartPoint.y);

}


/*!
 * Get North face midpoint for an interior cell.
 * This routine also works if the face is curved.
 *
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 */
Vector2D Grid2D_Quad_Block_HO::getMidPointFaceN(const int ii, const int jj) const{
  double sp_l, sp_r, dl, dr, ds_ratio, sp_m;
  Vector2D MidP;

  if (IsNorthFaceOfInteriorCellCurved(ii,jj)){

    MidP = xfaceN(ii,jj);
    sp_l = getS(nodeNW(ii,jj).X,
		BndNorthSpline);
    sp_r = getS(nodeNE(ii,jj).X,
		BndNorthSpline);
    dl = abs(MidP - nodeNW(ii,jj).X);
    dr = abs(nodeNE(ii,jj).X - MidP);
    ds_ratio = dl/(dl+dr);
    sp_m = sp_l + ds_ratio*(sp_r-sp_l);
    return Spline(sp_m, BndNorthSpline);

  } else {
    return xfaceN(ii,jj);
  }
}

/*!
 * Get South face midpoint for an interior cell.
 * This routine also works if the face is curved.
 *
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 */
Vector2D Grid2D_Quad_Block_HO::getMidPointFaceS(const int ii, const int jj) const{

  double sp_l, sp_r, dl, dr, ds_ratio, sp_m;
  Vector2D MidP;

  if (IsSouthFaceOfInteriorCellCurved(ii,jj)){
    
    MidP = xfaceS(ii,jj);
    sp_l = getS(nodeSW(ii,jj).X,
		BndSouthSpline);
    sp_r = getS(nodeSE(ii,jj).X,
		BndSouthSpline);
    dl = abs(MidP - nodeSW(ii,jj).X);
    dr = abs(nodeSE(ii,jj).X - MidP);
    ds_ratio = dl/(dl+dr);
    sp_m = sp_l + ds_ratio*(sp_r-sp_l);
    return Spline(sp_m, BndSouthSpline);

  } else {
    return xfaceS(ii,jj);
  }
}

/*!
 * Get East face midpoint for an interior cell.
 * This routine also works if the face is curved.
 *
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 */
Vector2D Grid2D_Quad_Block_HO::getMidPointFaceE(const int ii, const int jj) const{

  double sp_l, sp_r, dl, dr, ds_ratio, sp_m;
  Vector2D MidP;

  if (IsEastFaceOfInteriorCellCurved(ii,jj)){

    MidP = xfaceE(ii,jj);
    sp_l = getS(nodeSE(ii,jj).X, 
		BndEastSpline);
    sp_r = getS(nodeNE(ii,jj).X,
		BndEastSpline);
    dl = abs(MidP - nodeSE(ii,jj).X);
    dr = abs(nodeNE(ii,jj).X - MidP);
    ds_ratio = dl/(dl+dr);
    sp_m = sp_l + ds_ratio*(sp_r-sp_l);
    return Spline(sp_m, BndEastSpline);

  } else {
    return xfaceE(ii,jj);
  }
}

/*!
 * Get West face midpoint for an interior cell.
 * This routine also works if the face is curved.
 *
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 */
Vector2D Grid2D_Quad_Block_HO::getMidPointFaceW(const int ii, const int jj) const{

  double sp_l, sp_r, dl, dr, ds_ratio, sp_m;
  Vector2D MidP;

  if (IsWestFaceOfInteriorCellCurved(ii,jj)){

    MidP = xfaceW(ii,jj);
    sp_l = getS(nodeSW(ii,jj).X,
		BndWestSpline);
    sp_r = getS(nodeNW(ii,jj).X,
		BndWestSpline);
    dl = abs(MidP - nodeSW(ii,jj).X);
    dr = abs(nodeNW(ii,jj).X - MidP);
    ds_ratio = dl/(dl+dr);
    sp_m = sp_l + ds_ratio*(sp_r-sp_l);
    return Spline(sp_m, BndWestSpline);

  } else {
    return xfaceW(ii,jj);
  }
}

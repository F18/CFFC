/*!\file HO_Grid2DQuad.h
   \brief Header file defining high-order 2D quadrilateral block grid type. */

#ifndef _HO_GRID2D_QUAD_BLOCK_INCLUDED
#define _HO_GRID2D_QUAD_BLOCK_INCLUDED

/* Include required C++ libraries. */
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>

/* Using std namespace functions */
using namespace std;

/* Include CFFC header files */
#include "../Math/Math.h"	    // Include math macro header file.
#include "../CFD/CFD.h"		    // Include CFD header file
#include "../Math/Vector2D.h"	    // Include 2D vector header file
#include "HO_Spline2D.h"            // Include high-order 2D spline header file
#include "HO_Spline2DInterval.h"    // Include high-order 2D spline interval header file
#include "HO_Cell2D.h"		    // Include 2D cell header file
#include "HO_Node2D.h"		    // Include 2D node header file
#include "../MPI/MPI.h"		    // Include MPI header file
#include "../Math/LinearSystems.h"  // Include the linear systems header file.
#include "Grid2DQuadIntegration.h"  // Include the 2D quadrilateral domain integration class header file.
#include "Tecplot_ExecutionMode.h" // Include Tecplot control execution class.

/* Define the following types of 2D quadrilateral block 
   node initialization procedures. */

#define	GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH       1
#define	GRID2D_QUAD_BLOCK_INIT_PROCEDURE_EAST_WEST         2
#define	GRID2D_QUAD_BLOCK_INIT_PROCEDURE_TRANS_FINITE_XY   3
#define	GRID2D_QUAD_BLOCK_INIT_PROCEDURE_TRANS_FINITE_YX   4

/* Define the block sector indicators for
   2D quadrilateral block refinement and coarsening. */

#define	GRID2D_QUAD_BLOCK_SECTOR_NONE                     -1
#define	GRID2D_QUAD_BLOCK_SECTOR_SW                        0
#define	GRID2D_QUAD_BLOCK_SECTOR_SE                        1
#define	GRID2D_QUAD_BLOCK_SECTOR_NW                        2
#define	GRID2D_QUAD_BLOCK_SECTOR_NE                        3

/* Define orthogonality parameter */

#define ORTHOGONAL 1
#define NOT_ORTHOGONAL 0

/* Define possible positions for a cell
   in the 2D quadrilateral block. */

#define INTERIOR_CELL                               0
/* Boundary Cells (i.e. cells flanked by ONE spline in the specified direction) */
#define SOUTH_SPLINE                                1
#define NORTH_SPLINE                                2
#define EAST_SPLINE                                 3
#define WEST_SPLINE                                 4
/* Extension Boundary Cells (i.e. cells flanked by ONE extension spline in the specified direction)
   RIGHT or LEFT indicates whether the spline is present on the left or right side of the cell,
   relative to the position of the other nodes. */
#define EXTEND_W_NORTH_RIGHT_SPLINE                 5
#define EXTEND_E_NORTH_RIGHT_SPLINE                 6
#define EXTEND_W_NORTH_LEFT_SPLINE                  7
#define EXTEND_E_NORTH_LEFT_SPLINE                  8
#define EXTEND_W_SOUTH_RIGHT_SPLINE                 9
#define EXTEND_E_SOUTH_RIGHT_SPLINE                 10
#define EXTEND_W_SOUTH_LEFT_SPLINE                  11
#define EXTEND_E_SOUTH_LEFT_SPLINE                  12
#define EXTEND_N_EAST_RIGHT_SPLINE                  13
#define EXTEND_S_EAST_RIGHT_SPLINE                  14
#define EXTEND_N_EAST_LEFT_SPLINE                   15
#define EXTEND_S_EAST_LEFT_SPLINE                   16
#define EXTEND_N_WEST_RIGHT_SPLINE                  17
#define EXTEND_S_WEST_RIGHT_SPLINE                  18
#define EXTEND_N_WEST_LEFT_SPLINE                   19
#define EXTEND_S_WEST_LEFT_SPLINE                   20

/* Corner Boundary Cells (i.e. cells flanked in the specified directions by TWO splines) */
#define CORNER_NORTH_WEST_SPLINES                   21
#define CORNER_NORTH_EAST_SPLINES                   22
#define CORNER_SOUTH_WEST_SPLINES                   23
#define CORNER_SOUTH_EAST_SPLINES                   24
#define CORNER_NORTH_EXTEND_N_WEST_SPLINES          25
#define CORNER_EXTEND_N_WEST_EXTEND_W_NORTH_SPLINES 26
#define CORNER_EXTEND_W_NORTH_WEST_SPLINES          27
#define CORNER_WEST_EXTEND_W_SOUTH_SPLINES          28
#define CORNER_EXTEND_W_SOUTH_EXTEND_S_WEST_SPLINES 29
#define CORNER_EXTEND_S_WEST_SOUTH_SPLINES          30
#define CORNER_SOUTH_EXTEND_S_EAST_SPLINES          31
#define CORNER_EXTEND_S_EAST_EXTEND_E_SOUTH_SPLINES 32
#define CORNER_EXTEND_E_SOUTH_EAST_SPLINES          33
#define CORNER_EAST_EXTEND_E_NORTH_SPLINES          34
#define CORNER_EXTEND_E_NORTH_EXTEND_N_EAST_SPLINES 35
#define CORNER_EXTEND_N_EAST_NORTH_SPLINES          36


/* Define the high-order quadrilateral 2D grid block class. */

/*!
 * \class Grid2D_Quad_Block_HO
 *
 * @brief Definition of the 2D quadrilateral mesh block.
 *
 * \verbatim
 * Member functions
 *      NNi        -- Return number of nodes in
 *                    the i-direction (zeta-direction).
 *      INl        -- Return lower index for nodes in
 *                    the i-direction (zeta-direction).
 *      INu        -- Return upper index for nodes in
 *                    the i-direction (zeta-direction).
 *      NNj        -- Return number of nodes in
 *                    the j-direction (eta-direction).
 *      JNl        -- Return lower index for nodes in
 *                    the j-direction (eta-direction).
 *      JNu        -- Return upper index for nodes in
 *                    the j-direction (eta-direction).
 *      NCi        -- Return number of cells in
 *                    the i-direction (zeta-direction).
 *      ICl        -- Return lower index for cells in
 *                    the i-direction (zeta-direction).
 *      ICu        -- Return upper index for cells in
 *                    the i-direction (zeta-direction).
 *      NCj        -- Return number of cells in
 *                    the j-direction (eta-direction).
 *      JCl        -- Return lower index for cells in
 *                    the j-direction (eta-direction).
 *      JCu        -- Return upper index for cells in
 *                    the j-direction (eta-direction).
 *      Nghost     -- Return the number of ghost cells.
 *      Node       -- Return 2D node geometry.
 *      Cell       -- Return 2D cell geometry.
 *      BCtypeN    -- Return boundary condition
 *                    specifier for the north boundary
 *                    of grid block.
 *      BCtypeS    -- Return boundary condition
 *                    specifier for the south boundary
 *                    of grid block.
 *      BCtypeE    -- Return boundary condition
 *                    specifier for the east boundary
 *                    of grid block.
 *      BCtypeW    -- Return boundary condition
 *                    specifier for the west boundary
 *                    of grid block.
 *      BndNorthSpline
 *                 -- Return 2D spline defining north 
 *                    boundary of grid block.
 *      BndSouthSpline 
 *                 -- Return 2D spline defining south
 *                    boundary of grid block.
 *      BndEastSpline 
 *                 -- Return 2D spline defining east 
 *                    boundary of grid block.
 *      BndWestSpline
 *                 -- Return 2D spline defining west
 *                    boundary of grid block.
 *      SminN,SmaxN
 *                 -- Return min and max values of
 *                    pathlengths for north boundary
 *                    spline.
 *      SminS,SmaxS
 *                 -- Return min and max values of
 *                    pathlengths for south boundary
 *                    spline.
 *      SminE,SmaxE
 *                 -- Return min and max values of
 *                    pathlengths for east boundary
 *                    spline.
 *      SminW,SmaxW
 *                 -- Return min and max values of
 *                    pathlengths for west boundary
 *                    spline.
 *      StretchI   -- Return the grid point stretching
 *                    function for the i-direction.
 *      BetaI      -- Return i-direction stretching
 *                    function parameter.
 *      TauI       -- Return i-direction stretching
 *                    function parameter.
 *      StretchJ   -- Return the grid point stretching
 *                    function for the j-direction.
 *      BetaJ      -- Return j-direction stretching
 *                    function parameter.
 *      TauJ       -- Return j-direction stretching
 *                    function parameter.
 *      OrthogonalN -- Return the grid orthogonality
 *                     parameter for north boundary.
 *      OrthogonalS -- Return the grid orthogonality
 *                     parameter for south boundary.
 *      OrthogonalE -- Return the grid orthogonality
 *                     parameter for east boundary.
 *      OrthogonalW -- Return the grid orthogonality
 *                     parameter for west boundary.
 *      allocate   -- Allocate memory for structured
 *                    quadrilateral grid block.
 *      deallocate -- Deallocate memory for structured
 *                    quadrilateral grid block.
 *      centroid   -- Calculate 2D vector containing
 *                    the location of cell center.
 *      area       -- Calculate the area for cell.
 *      nodeNW     -- Return NW node for cell.
 *      nodeNE     -- Return NE node for cell.
 *      nodeSE     -- Return SE node for cell.
 *      nodeSW     -- Return SW node for cell.
 *      xfaceN     -- Return midpoint of
 *                    cell north face.
 *      xfaceS     -- Return midpoint of
 *                    cell south face.
 *      xfaceE     -- Return midpoint of
 *                    cell east face.
 *      xfaceW     -- Return midpoint of
 *                    cell west face.
 *      lfaceN     -- Return the length of the
 *                    cell north face.
 *      lfaceS     -- Return the length of the
 *                    cell south face.
 *      lfaceE     -- Return the length of the
 *                    cell east face.
 *      lfaceW     -- Return the length of the
 *                    cell west face.
 *      cell_perimeter  -- Return Perimeter of the cell
 *      nfaceN     -- Return the unit vector in the
 *                    direction normal to the cell
 *                    north face.
 *      nfaceS     -- Return the unit vector in the
 *                    direction normal to the cell
 *                    south face.
 *      nfaceE     -- Return the unit vector in the
 *                    direction normal to the cell
 *                    east face.
 *      nfaceW     -- Return the unit vector in the
 *                    direction normal to the cell
 *                    west face.
 *       dNdC     --  Returns the influence coefficient 
 *                    of cell solution values on 
 *                    nodal values.
 *    dFacedC     --  Returns the influence coefficients
 *                    of cell solution values on
 *                    face centrered values.
 * dDiamondPathdC --  Returns the influence coefficients
 *                    of cell solution on diamond path
 *                    solution gradient reconstruction
 *                    at the cell faces.
 *
 * Member operators
 *  G  -- grid block consisting of quadrilateral cells
 *  V  -- a 2D vector
 *  a  -- a scalar (double)
 *
 * G = G;
 * G = G + V; (shift location of grid)
 * G = G - V; (shift location of grid)
 * G = a * G; (scale grid)
 * G = G * a; (scale grid)
 * G = G / a; (scale grid)
 * G = G ^ a; (rotate grid)
 * cout << G; (output function)
 * cin  >> G; (input function)
 * \endverbatim
 */
class Grid2D_Quad_Block_HO{
public:
  //! @name Defined datatypes
  //@{
  typedef Node2D_HO NodeType;
  typedef Cell2D_HO::GeometricMoments GeometricMoments;
  typedef GeometricMoments::Derivative  GeomMoment;
  typedef Spline2D_HO BndSplineType;
  typedef Spline2DInterval_HO BndSplineIntervalType;
  //@}
  
  //! @name Mesh indexes
  //@{ 
  int                   NNi, //!< Number of nodes in i-direction (zeta-direction).
                        INl, //!< Lower index for nodes in i-direction (zeta-direction).
                        INu; //!< Upper index for nodes in i-direction (zeta-direction).
  int                   NNj, //!< Number of nodes in j-direction (eta-direction).
                        JNl, //!< Lower index for nodes in j-direction (eta-direction).
                        JNu; //!< Upper index for nodes in j-direction (eta-direction).
  int                   NCi, //!< Number of cells in i-direction (zeta-direction).
                        ICl, //!< Lower index for cells in i-direction (zeta-direction).
                        ICu; //!< Upper index for cells in i-direction (zeta-direction).
  int                   NCj, //!< Number of cells in j-direction (eta-direction).
                        JCl, //!< Lower index for cells in j-direction (eta-direction).
                        JCu; //!< Upper index for cells in j-direction (eta-direction).
  int                Nghost; //!< Number of ghost cells.
  //@}

  //! @name Location of nodes and cells in the structured mesh
  //@{ 
  Node2D_HO             **Node; //!< Array of 2D node position vectors
  Cell2D_HO             **Cell; //!< Array of 2D cell centre vectors
  //@}

  //! @name Block boundary condition specifiers
  //@{ 
  int                 *BCtypeN, //!< North boundary condition type specifier.
                      *BCtypeS, //!< South boundary condition type specifier.
                      *BCtypeE, //!< East boundary condition type specifier.
                      *BCtypeW; //!< West boundary condition type specifier.
  //@}

  //! @name Block geometry parameters
  //@{ 
  Spline2D_HO   BndNorthSpline, //!< North boundary 2D spline.
                BndSouthSpline, //!< South boundary 2D spline.
                BndEastSpline,  //!< East boundary 2D spline.
                BndWestSpline;  //!< West boundary 2D spline.

  Spline2D_HO   ExtendWest_BndNorthSpline, //!< Extension of North boundary spline to west
                ExtendEast_BndNorthSpline, //!< Extension of North boundary spline to east
                ExtendWest_BndSouthSpline, //!< Extension of South boundary spline to west
                ExtendEast_BndSouthSpline, //!< Extension of South boundary spline to east
                ExtendNorth_BndEastSpline, //!< Extension of East boundary spline to north
                ExtendSouth_BndEastSpline, //!< Extension of East boundary spline to south
                ExtendNorth_BndWestSpline, //!< Extension of West boundary spline to north
                ExtendSouth_BndWestSpline; //!< Extension of West boundary spline to south
                
  Spline2DInterval_HO *BndNorthSplineInfo, //!< North boundary 2D spline info.
                      *BndSouthSplineInfo, //!< South boundary 2D spline info.
                      *BndEastSplineInfo,  //!< East boundary 2D spline info.
                      *BndWestSplineInfo;  //!< West boundary 2D spline info.

  Spline2DInterval_HO   *ExtendWest_BndNorthSplineInfo, //!< Info for extension of North boundary spline to west
                        *ExtendEast_BndNorthSplineInfo, //!< Info for extension of North boundary spline to east
                        *ExtendWest_BndSouthSplineInfo, //!< Info for extension of South boundary spline to west
                        *ExtendEast_BndSouthSplineInfo, //!< Info for extension of South boundary spline to east
                        *ExtendNorth_BndEastSplineInfo, //!< Info for extension of East boundary spline to north
                        *ExtendSouth_BndEastSplineInfo, //!< Info for extension of East boundary spline to south
                        *ExtendNorth_BndWestSplineInfo, //!< Info for extension of West boundary spline to north
                        *ExtendSouth_BndWestSplineInfo; //!< Info for extension of West boundary spline to south


  double                SminN,  //!< Minimum value of north face pathlength.
                        SmaxN,  //!< Maximum value of north face pathlength.
                        SminS,  //!< Minimum value of south face pathlength.
                        SmaxS,  //!< Maximum value of south face pathlength.
                        SminE,  //!< Minimum value of east face pathlength.
                        SmaxE,  //!< Maximum value of east face pathlength.
                        SminW,  //!< Minimum value of west face pathlength.
                        SmaxW;  //!< Maximum value of west face pathlength.
  //@}


  //! @name Mesh parameters
  //@{  
  int                StretchI, //!< Grid point stretching function (i-direction).
                     StretchJ; //!< Grid point stretching function (j-direction).
  double                BetaI, //!< Grid point stretching function parameter (i-direction).
                         TauI, //!< Grid point stretching function parameter (i-direction).
                        BetaJ, //!< Grid point stretching function parameter (j-direction).
                         TauJ; //!< Grid point stretching function parameter (j-direction).
  int             OrthogonalN, //!< North boundary orthogonality parameter.
                  OrthogonalS, //!< South boundary orthogonality parameter.
                  OrthogonalE, //!< East boundary orthogonality parameter.
                  OrthogonalW; //!< West boundary orthogonality parameter.
  //@}


  //! @name Geometric integration
  //@{  
  Grid2DQuadIntegration<Grid2D_Quad_Block_HO> Integration; //!< Variable that provides access to integration subroutines
  //@}
			
  //! @name Constructors, destructor and assignment operator
  //@{  
  //! Default constructor.
  Grid2D_Quad_Block_HO(void);

  //! Destructor.
  ~Grid2D_Quad_Block_HO(void){ deallocate(); }

  //! Assignment operator.
  Grid2D_Quad_Block_HO& operator= (const Grid2D_Quad_Block_HO &Grid);
  //@}

  //! @name Memory allocation and deallocation
  //@{
  //! Allocate memory for structured quadrilateral grid block.
  void allocate(const int &Ni, const int &Nj, const int &Ng);

  //! Allocate memory for structured quadrilateral grid block with reconstruction order.
  void allocate(const int &Ni, const int &Nj, const int &Ng, const int &HighestRecOrder);

  //! Deallocate memory for structured quadrilateral grid block.
  void deallocate(void);

  //! Deallocate memory for spline info(s) of structured quadrilateral grid block.
  void deallocateBndSplineInfo(void);
  //! Deallocate memory for North spline info
  void deallocate_BndNorthSplineInfo(void){ delete [] BndNorthSplineInfo; BndNorthSplineInfo = NULL; }
  //! Deallocate memory for South spline info
  void deallocate_BndSouthSplineInfo(void){ delete [] BndSouthSplineInfo; BndSouthSplineInfo = NULL; }
  //! Deallocate memory for North spline info
  void deallocate_BndEastSplineInfo(void){ delete [] BndEastSplineInfo; BndEastSplineInfo = NULL; }
  //! Deallocate memory for North spline info
  void deallocate_BndWestSplineInfo(void){ delete [] BndWestSplineInfo; BndWestSplineInfo = NULL; }

  //! Deallocate memory for the info of North spline extension to West
  void deallocate_ExtendWest_BndNorthSplineInfo(void){
    delete [] ExtendWest_BndNorthSplineInfo; ExtendWest_BndNorthSplineInfo = NULL; }
  void deallocate_ExtendEast_BndNorthSplineInfo(void){
    delete [] ExtendEast_BndNorthSplineInfo; ExtendEast_BndNorthSplineInfo = NULL; }
  void deallocate_ExtendWest_BndSouthSplineInfo(void){
    delete [] ExtendWest_BndSouthSplineInfo; ExtendWest_BndSouthSplineInfo = NULL; }
  void deallocate_ExtendEast_BndSouthSplineInfo(void){
    delete [] ExtendEast_BndSouthSplineInfo; ExtendEast_BndSouthSplineInfo = NULL; }
  void deallocate_ExtendNorth_BndEastSplineInfo(void){
    delete [] ExtendNorth_BndEastSplineInfo; ExtendNorth_BndEastSplineInfo = NULL; }
  void deallocate_ExtendSouth_BndEastSplineInfo(void){
    delete [] ExtendSouth_BndEastSplineInfo; ExtendSouth_BndEastSplineInfo = NULL; }
  void deallocate_ExtendNorth_BndWestSplineInfo(void){
    delete [] ExtendNorth_BndWestSplineInfo; ExtendNorth_BndWestSplineInfo = NULL; }
  void deallocate_ExtendSouth_BndWestSplineInfo(void){
    delete [] ExtendSouth_BndWestSplineInfo; ExtendSouth_BndWestSplineInfo = NULL; }

  //! Remove the extension splines
  void RemoveExtensionSplines(void);
  //@}

  //! @name Calculate cell centroid.
  //@{
  int quadAreaAndCentroid(const int &ii, const int &jj) const;
  Vector2D quadConvexCentroid(const int &ii, const int &jj) const;
  Vector2D centroid(const Cell2D_HO &Cell);
  Vector2D centroid(const int &ii, const int &jj);
  Vector2D centroidSW(const int &ii, const int &jj);
  Vector2D centroidSE(const int &ii, const int &jj);
  Vector2D centroidNW(const int &ii, const int &jj);
  Vector2D centroidNE(const int &ii, const int &jj);
  Vector2D centroidSW_ConvexQuad(const int &ii, const int &jj) const;
  Vector2D centroidSE_ConvexQuad(const int &ii, const int &jj) const;
  Vector2D centroidNW_ConvexQuad(const int &ii, const int &jj) const;
  Vector2D centroidNE_ConvexQuad(const int &ii, const int &jj) const;
  Vector2D centroid_CurvedBoundaries(const int &CellIndex, const int &Boundary) const;
  Vector2D centroid_GhostCell_CurvedBoundaries(const int &CellIndex, const int &Boundary) const;
  //@}

  //! @name Get area and cell centroid.
  //@{
  //! Access the centroid of cell (ii,jj)
  const Vector2D & CellCentroid(int ii, int jj) const {return Cell[ii][jj].Xc; }
  //! Access the x-coordinate of the centroid of cell (ii,jj)
  const double & XCellCentroid(int ii, int jj) const {return Cell[ii][jj].Xc.x; }
  //! Access the y-coordinate of the centroid of cell (ii,jj)
  const double & YCellCentroid(int ii, int jj) const {return Cell[ii][jj].Xc.y; }
  //! Access the area of cell (ii,jj)
  const double & CellArea(const int &ii, const int &jj) const {return Cell[ii][jj].A; }
  //@}

  //! @name Calculate cell area.
  //@{
  double area(const Cell2D_HO &Cell) const;
  double area(const int &ii, const int &jj) const;
  double area(const Node2D_HO &SW, const Node2D_HO &NW, const Node2D_HO &NE, const Node2D_HO &SE);
  double area_CurvedBoundaries(const int &CellIndex, const int &Boundary) const;
  double area_GhostCell_CurvedBoundaries(const int &CellIndex, const int &Boundary) const;
  //@}

  //! @name Get cell nodes.
  //@{
  Node2D_HO nodeNW(const Cell2D_HO &Cell) const;
  Node2D_HO nodeNW(const int ii, const int jj) const;

  Node2D_HO nodeNE(const Cell2D_HO &Cell) const;
  Node2D_HO nodeNE(const int ii, const int jj) const;

  Node2D_HO nodeSE(const Cell2D_HO &Cell) const;
  Node2D_HO nodeSE(const int ii, const int jj) const;

  Node2D_HO nodeSW(const Cell2D_HO &Cell) const;
  Node2D_HO nodeSW(const int ii, const int jj) const;
  //@}

  //! @name Get cell face midpoints.
  //@{
  Vector2D xfaceN(const Cell2D_HO &Cell) const;
  Vector2D xfaceN(const int ii, const int jj) const;

  Vector2D xfaceS(const Cell2D_HO &Cell) const;
  Vector2D xfaceS(const int ii, const int jj) const;

  Vector2D xfaceE(const Cell2D_HO &Cell) const;
  Vector2D xfaceE(const int ii, const int jj) const;

  Vector2D xfaceW(const Cell2D_HO &Cell) const;
  Vector2D xfaceW(const int ii, const int jj) const;
  //@}

  //! @name Get cell face lengths.
  //@{
  double lfaceN(const Cell2D_HO &Cell) const;
  double lfaceN(const int ii, const int jj) const;

  double lfaceS(const Cell2D_HO &Cell) const;
  double lfaceS(const int ii, const int jj) const;

  double lfaceE(const Cell2D_HO &Cell) const;
  double lfaceE(const int ii, const int jj) const;

  double lfaceW(const Cell2D_HO &Cell) const;
  double lfaceW(const int ii, const int jj) const;
  //@}

  //! @name Get cell perimeter.
  //@{
  double cell_perimeter(const int ii, const int jj) const;
  //@}

  //! @name Get the unit vector normal to the cell faces.
  //@{
  Vector2D nfaceN(const Cell2D_HO &Cell) const;
  Vector2D nfaceN(const int ii, const int jj) const;

  Vector2D nfaceS(const Cell2D_HO &Cell) const;
  Vector2D nfaceS(const int ii, const int jj) const;

  Vector2D nfaceE(const Cell2D_HO &Cell) const;
  Vector2D nfaceE(const int ii, const int jj) const;

  Vector2D nfaceW(const Cell2D_HO &Cell) const;
  Vector2D nfaceW(const int ii, const int jj) const;
  //@}

  //! @name Get Gauss quadrature points for each straight cell face
  //@{
  void getGaussQuadPointsFaceN(const Cell2D_HO &Cell, Vector2D * GQPoints, const int & NumberOfGQPs) const;
  void getGaussQuadPointsFaceN(const int &ii, const int &jj, Vector2D * GQPoints, const int & NumberOfGQPs) const;
  void addGaussQuadPointsFaceN(const int &ii, const int &jj, std::vector<Vector2D> &GQPoints, const int & NumberOfGQPs) const;

  void getGaussQuadPointsFaceS(const Cell2D_HO &Cell, Vector2D * GQPoints, const int & NumberOfGQPs) const;
  void getGaussQuadPointsFaceS(const int &ii, const int &jj, Vector2D * GQPoints, const int & NumberOfGQPs) const;
  void addGaussQuadPointsFaceS(const int &ii, const int &jj, std::vector<Vector2D> &GQPoints, const int & NumberOfGQPs) const;

  void getGaussQuadPointsFaceE(const Cell2D_HO &Cell, Vector2D * GQPoints, const int & NumberOfGQPs) const;
  void getGaussQuadPointsFaceE(const int &ii, const int &jj, Vector2D * GQPoints, const int & NumberOfGQPs) const;
  void addGaussQuadPointsFaceE(const int &ii, const int &jj, std::vector<Vector2D> &GQPoints, const int & NumberOfGQPs) const;

  void getGaussQuadPointsFaceW(const Cell2D_HO &Cell, Vector2D * GQPoints, const int & NumberOfGQPs) const;
  void getGaussQuadPointsFaceW(const int &ii, const int &jj, Vector2D * GQPoints, const int & NumberOfGQPs) const;
  void addGaussQuadPointsFaceW(const int &ii, const int &jj, std::vector<Vector2D> &GQPoints, const int & NumberOfGQPs) const;
  //@}
  
  //! @name Number of Gauss quadrature points used for flux calculation.
  //@{
  //! Set NumGQP to the passed number.
  void SetNumberOfGaussQuadraturePoints(const int & _NumGQP_){ NumGQP = _NumGQP_; }
  //! Set NumGQP based on the reconstruction order to obtain the required solution accuracy
  void SetNumberOfFluxCalculationGaussQuadraturePoints(const int & RecOrder);
  //! Set NumGQP based on correlations.
  void SetNumberOfGaussQuadraturePoints(void);
  //! Get NumGQP
  const int & getNumGQP(void) const {return NumGQP; }
  //@}

  //! @name Get number of constrained Gauss quadrature points for each cell.
  //@{
  int NumOfConstrainedGaussQuadPoints_North(const Cell2D_HO &Cell){
    return NumOfConstrainedGaussQuadPoints_North(Cell.I,Cell.J);
  }
  int NumOfConstrainedGaussQuadPoints_North(const int &ii, const int &jj);

  int NumOfConstrainedGaussQuadPoints_South(const Cell2D_HO &Cell){
    return NumOfConstrainedGaussQuadPoints_South(Cell.I,Cell.J);
  }
  int NumOfConstrainedGaussQuadPoints_South(const int &ii, const int &jj);

  int NumOfConstrainedGaussQuadPoints_East(const Cell2D_HO &Cell){
    return NumOfConstrainedGaussQuadPoints_East(Cell.I,Cell.J);
  }
  int NumOfConstrainedGaussQuadPoints_East(const int &ii, const int &jj);

  int NumOfConstrainedGaussQuadPoints_West(const Cell2D_HO &Cell){
    return NumOfConstrainedGaussQuadPoints_West(Cell.I,Cell.J);
  }
  int NumOfConstrainedGaussQuadPoints_West(const int &ii, const int &jj);

  int NumOfConstrainedGaussQuadPoints(const Cell2D_HO &Cell){ return NumOfConstrainedGaussQuadPoints(Cell.I,Cell.J);}
  int NumOfConstrainedGaussQuadPoints(const int &ii, const int &jj);
  //@}

  //! @name Get number of flux calculation Gauss quadrature points for each cell.
  //@{
  int NumOfFluxCalculationGaussQuadPoints_North(const Cell2D_HO &Cell, bool & curved_face){
    return NumOfFluxCalculationGaussQuadPoints_North(Cell.I,Cell.J,curved_face);
  }
  int NumOfFluxCalculationGaussQuadPoints_North(const int &ii, const int &jj, bool & curved_face);

  int NumOfFluxCalculationGaussQuadPoints_South(const Cell2D_HO &Cell, bool & curved_face){
    return NumOfFluxCalculationGaussQuadPoints_South(Cell.I,Cell.J,curved_face);
  }
  int NumOfFluxCalculationGaussQuadPoints_South(const int &ii, const int &jj, bool & curved_face);

  int NumOfFluxCalculationGaussQuadPoints_East(const Cell2D_HO &Cell, bool & curved_face){
    return NumOfFluxCalculationGaussQuadPoints_East(Cell.I,Cell.J,curved_face);
  }
  int NumOfFluxCalculationGaussQuadPoints_East(const int &ii, const int &jj, bool & curved_face);

  int NumOfFluxCalculationGaussQuadPoints_West(const Cell2D_HO &Cell, bool & curved_face){
    return NumOfFluxCalculationGaussQuadPoints_West(Cell.I,Cell.J,curved_face);
  }
  int NumOfFluxCalculationGaussQuadPoints_West(const int &ii, const int &jj, bool & curved_face);

  int NumOfFluxCalculationGaussQuadPoints(const Cell2D_HO &Cell){
    return NumOfFluxCalculationGaussQuadPoints(Cell.I,Cell.J);
  }
  int NumOfFluxCalculationGaussQuadPoints(const int &ii, const int &jj);
  //@}

  //! @name Get cell geometric moments and related information.
  //@{
  //! Get the current highest(i.e. maximum) reconstruction order
  const int & MaxRecOrder(void) const { return HighestReconstructionOrder; }
  //! Get the geometric coefficient array of cell (ii,jj).
  const GeometricMoments & CellGeomCoeff(const int &ii, const int &jj) const { return Cell[ii][jj].GeomCoeff(); }
  //! Get the value of cell (ii,jj) geometric coefficient with x-power 'p1' and y-power 'p2'.
  const double & CellGeomCoeffValue(const int &ii, const int &jj, const int &p1, const int &p2) {
    return Cell[ii][jj].GeomCoeffValue(p1,p2);
  }
  //! Get the value of cell (ii,jj) geometric coefficient which is store in the 'position' place 
  const double & CellGeomCoeffValue(const int &ii, const int &jj, const int &position) const {
    return Cell[ii][jj].GeomCoeffValue(position);
  }
  //! Get the cell (ii,jj) geometric coefficient which is store in the 'position' place (i.e. powers and value) 
  GeomMoment & CellGeomCoeff(const int &ii, const int &jj, const int &position) {
    return Cell[ii][jj].GeomCoeff(position);
  }
  //@}

  //! @name Calculate cell geometric moments and related information.
  //@{
  void ComputeGeometricCoefficients(const Cell2D_HO &Cell){ ComputeGeometricCoefficients(Cell.I,Cell.J);}
  void ComputeGeometricCoefficients(const int &ii, const int &jj);
  void ComputeGeometricCoefficients_CurvedBoundaries(const int &CellIndexI, const int &CellIndexJ,
						     const int &Boundary);
  void ComputeGeometricCoefficients_GhostCell_CurvedBoundaries(const int &CellIndexI, const int &CellIndexJ,
							       const int &Boundary);
  double GeometricMoment_CurvedBoundaries(const int &CellIndexI, const int &CellIndexJ,
					  const int &Boundary,
					  const int &OrderX, const int &OrderY);
  double GeometricMoment_GhostCell_CurvedBoundaries(const int &CellIndexI, const int &CellIndexJ,
						    const int &Boundary,
						    const int &OrderX, const int &OrderY);
  //@}

  //! @name Bilinear interplation (Zingg & Yarrow) and diamond path reconstruction.
  //@{
  //! Returns the influence coefficient of cell solution on nodal value.
  double dNdC(const int ii, const int jj, const int cell_orientation) const;
  //! Returns the influence coefficients of cell solution on face value.
  void dFacedC(double &left_node_weight, double &right_node_weight, 
	       const int i, const int j, const int face, 
	       const int cell_orientation) const;
  //! Returns the influence coefficients for diamond path gradient reconstruction.
  void dDiamondPathdC(double &d_dWdx_dW, double &d_dWdy_dW, 
		      const int i, const int j, const int face, 
		      const double &face_left_node_weight, 
		      const double &face_right_node_weight,
		      const int cell_orientation) const;
  //@}

  //! Change current BC's
  //@{
  void set_BCs(const int& FACE, const int& BC);
  //@}


  //! Create quad block
  //@{
  void Create_Quad_Block(const int Number_of_Cells_Idir,
			 const int Number_of_Cells_Jdir,
			 const int Number_of_Ghost_Cells,
			 const int Highest_Order_of_Reconstruction);
  void Create_Quad_Block_Without_Update(const int &Number_of_Cells_Idir,
					const int &Number_of_Cells_Jdir,
					const int &Number_of_Ghost_Cells,
					const int &Highest_Order_of_Reconstruction);
  void Create_Quad_Block(char *Bnd_Spline_File_Name_ptr,
			 const int Number_of_Cells_Idir,
			 const int Number_of_Cells_Jdir,
			 const int Number_of_Ghost_Cells,
			 const int Highest_Order_of_Reconstruction);
  void Create_Quad_Block_Without_Update(char *Bnd_Spline_File_Name_ptr,
					const int &Number_of_Cells_Idir,
					const int &Number_of_Cells_Jdir,
					const int &Number_of_Ghost_Cells,
					const int &Highest_Order_of_Reconstruction);
  void Create_Quad_Block(Spline2D_HO &Bnd_Spline_North,
			 Spline2D_HO &Bnd_Spline_South,
			 Spline2D_HO &Bnd_Spline_East,
			 Spline2D_HO &Bnd_Spline_West,
			 const int Number_of_Cells_Idir,
			 const int Number_of_Cells_Jdir,
			 const int Number_of_Ghost_Cells,
			 const int Highest_Order_of_Reconstruction,
			 const int Node_Init_Procedure,
			 const int Stretch_I,
			 const double &Beta_I, 
			 const double &Tau_I,
			 const int Stretch_J,
			 const double &Beta_J,
			 const double &Tau_J,
			 const int Orthogonal_North,
			 const int Orthogonal_South,
			 const int Orthogonal_East,
			 const int Orthogonal_West);
  void Create_Quad_Block_Without_Update(Spline2D_HO &Bnd_Spline_North,
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
					const int &Stretch_J,
					const double &Beta_J,
					const double &Tau_J,
					const int &Orthogonal_North,
					const int &Orthogonal_South,
					const int &Orthogonal_East,
					const int &Orthogonal_West);
  //@}

  //! @name Broadcast functions (MPI)
  //@{
  void Broadcast_Quad_Block(void);
  friend void Broadcast_Quad_Block(Grid2D_Quad_Block_HO &Grid){ return Grid.Broadcast_Quad_Block(); }
  
#ifdef _MPI_VERSION
  void Broadcast_Quad_Block(MPI::Intracomm &Communicator, 
			    const int Source_CPU);
  friend void Broadcast_Quad_Block(Grid2D_Quad_Block_HO &Grid,
				   MPI::Intracomm &Communicator, 
				   const int Source_CPU){
    return Grid.Broadcast_Quad_Block(Communicator, Source_CPU);
  }

#endif
  //@}
  
  //! @name Smooth quad grid
  //@{
  static void setGridSmoothing(void){ Smooth_Quad_Block_Flag = ON; }
  static void setNoGridSmoothing(void) { Smooth_Quad_Block_Flag = OFF; }
  void Smooth_Quad_Block(const int Number_of_Iterations);
  friend void Smooth_Quad_Block(Grid2D_Quad_Block_HO &Grid,
				const int Number_of_Iterations){
    return Grid.Smooth_Quad_Block(Number_of_Iterations);
  }
  
  void Smooth_Rocket_Motor(const double &Length_Chamber,
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
			   const int &NRi,const int &NRj);
  //@}
  
  //!@name Set Boundary Conditions.
  //@{
  void Set_BCs(void);
  friend void Set_BCs(Grid2D_Quad_Block_HO &Grid){ return Grid.Set_BCs(); }
  //@}
  
  //!@name Update geometry and geometric properties on different levels of update.
  //@{
  //! Set the interior mesh update flag to require update of the geoemtric properties of the interior cells.
  void Schedule_Interior_Mesh_Update(void){ InteriorMeshUpdate = ON; }
  //! Set the ghost cells update flag to require update of the geometric properties of these cells.
  void Schedule_Ghost_Cells_Update(void){ GhostCellsUpdate = ON; }
  //! Set the corner ghost cells update flag to require update of the geometric properties of these cells.
  void Schedule_Corner_Ghost_Cells_Update(void){ CornerGhostCellsUpdate = ON; }
  //! Confirm that the update of interior cells has been done and reset the flag accordingly
  void Confirm_Interior_Mesh_Update(void){ InteriorMeshUpdate = OFF; }
  //! Confirm that the update of ghost cells has been done and reset the flag accordingly
  void Confirm_Ghost_Cells_Update(void){ GhostCellsUpdate = OFF; CornerGhostCellsUpdate = OFF;}
  //! Confirm that the update of corner ghost cells has been done and reset the flag accordingly
  void Confirm_Corner_Ghost_Cells_Update(void){ CornerGhostCellsUpdate = OFF; }
  //! Confirm that the update of ALL cells has been done and reset the flags accordingly
  void Confirm_Mesh_Update_Everywhere(void){ Reset_Mesh_Update_Flags(); }

  void Update_Exterior_Nodes(void);
  friend void Update_Exterior_Nodes(Grid2D_Quad_Block_HO &Grid){ return Grid.Update_Exterior_Nodes(); }
  
  void Update_Corner_Ghost_Nodes(void);
  friend void Update_Corner_Ghost_Nodes(Grid2D_Quad_Block_HO &Grid){ return Grid.Update_Corner_Ghost_Nodes(); }
  
  void Update_Cell(const int & iCell, const int & jCell);
  void Update_GhostCellProperties_DuringMessagePassing(const int &iCell, const int &jCell);

  void Update_Cells(void);
  friend void Update_Cells(Grid2D_Quad_Block_HO &Grid){ return Grid.Update_Cells(); }

  void Update_All_Cells(void);
  void Update_Interior_Cells(void);
  void Update_Ghost_Cells(void);
  void Update_Corner_Ghost_Cells(void);
  void Update_SplineInfos(void);

  void Update_Grid_Properties(const int &HighestRecOrder);

  int Check_Quad_Block(void);
  int Check_Quad_Block_Completely(void);

  const int & Value_InteriorMeshUpdate_Flag(void) const { return InteriorMeshUpdate; }
  const int & Value_GhostCellsUpdate_Flag(void) const { return GhostCellsUpdate; }
  const int & Value_CornerGhostCellsUpdate_Flag(void) const { return CornerGhostCellsUpdate; }

  const unsigned int & getInteriorStateTracker(void) const { return InteriorCellGeometryStateTracker; }
  const unsigned int & getGhostStateTracker(void) const { return GhostCellGeometryStateTracker; }
  const unsigned int & getCornerGhostStateTracker(void) const { return CornerGhostCellGeometryStateTracker; }
  //@}
  
  //!@name Input/Output functions
  //@{
  void Write_Quad_Block_Definition(ostream &Out_File);
  
  void Read_Quad_Block_Definition(istream &In_File);
  
  void Write_Quad_Block(ostream &Out_File);

  void Read_Quad_Block(istream &In_File);
  //@}

  //!@name Copy block
  //@{
  //! Another way of copying a grid
  void Copy_Quad_Block(const Grid2D_Quad_Block_HO &Grid2){ *this = Grid2; }
  //! Copy Grid2 into Grid1
  friend void Copy_Quad_Block(Grid2D_Quad_Block_HO &Grid1,
			      const Grid2D_Quad_Block_HO &Grid2){ Grid1 = Grid2; }
  //@}
  
  //!@name Block manipulation
  //@{
  void Translate_Quad_Block(const Vector2D &V);
  void Translate_Quad_Block_Without_Update(const Vector2D &V);
  
  void Scale_Quad_Block(const double &Scaling_Factor);
  void Scale_Quad_Block_Without_Update(const double &Scaling_Factor);
  
  void Rotate_Quad_Block(const double &Angle);
  void Rotate_Quad_Block_Without_Update(const double &Angle);
  
  void Reflect_Quad_Block(void);
  void Reflect_Quad_Block_Without_Update(void);

  void Disturb_Interior_Nodes(const int &Number_of_Iterations);
  void Disturb_Interior_Nodes_Without_Update(const int &Number_of_Iterations);
  //@}
  
  //!@ Output functions for plotting.
  //@{
  void Output_Tecplot(const int Block_Number,
		      const int Output_Title,
		      ostream &Out_File) const;
  friend void Output_Tecplot(const Grid2D_Quad_Block_HO &Grid,
			     const int Block_Number,
			     const int Output_Title,
			     ostream &Out_File){
    return Grid.Output_Tecplot(Block_Number, Output_Title, Out_File);
  }
  
  void Output_Nodes_Tecplot(const int Block_Number,
			    const int Output_Title,
			    ostream &Out_File) const;
  friend void Output_Nodes_Tecplot(const Grid2D_Quad_Block_HO &Grid,
				   const int Block_Number,
				   const int Output_Title,
				   ostream &Out_File){
    return Grid.Output_Nodes_Tecplot(Block_Number, Output_Title, Out_File);
  }
  
  void Output_Cells_Tecplot(const int Block_Number,
			    const int Output_Title,
			    ostream &Out_File) const;
  friend void Output_Cells_Tecplot(const Grid2D_Quad_Block_HO &Grid,
				   const int Block_Number,
				   const int Output_Title,
				   ostream &Out_File){
    return Grid.Output_Cells_Tecplot(Block_Number, Output_Title, Out_File);
  }
  
  void Output_Gnuplot(const int Block_Number,
		      const int Output_Title,
		      ostream &Out_File) const;
  friend void Output_Gnuplot(const Grid2D_Quad_Block_HO &Grid,
			     const int Block_Number,
			     const int Output_Title,
			     ostream &Out_File){
    return Grid.Output_Gnuplot(Block_Number, Output_Title, Out_File);
  }
  //@}
  
  //!@name AMR related functions
  //@{
  void Double_Mesh_Resolution(const Grid2D_Quad_Block_HO &Grid_Original);
  //! Generate Grid_Double with double the mesh resolution of Grid_Original
  friend void Double_Mesh_Resolution(Grid2D_Quad_Block_HO &Grid_Double,
				     const Grid2D_Quad_Block_HO &Grid_Original){
    return Grid_Double.Double_Mesh_Resolution(Grid_Original);
  }

  void Half_Mesh_Resolution(const Grid2D_Quad_Block_HO &Grid_Original);
  //! Generate Grid_Half with half the mesh resolution of Grid_Original
  friend void Half_Mesh_Resolution(Grid2D_Quad_Block_HO &Grid_Half,
				   const Grid2D_Quad_Block_HO &Grid_Original){
    return Grid_Half.Half_Mesh_Resolution(Grid_Original);
  }
  
  void Refine_Mesh(const Grid2D_Quad_Block_HO &Grid_Original,
		   const int Sector);
  //! Assign the Sector of the refined Grid_Original to Grid_Fine
  friend void Refine_Mesh(Grid2D_Quad_Block_HO &Grid_Fine,
			  const Grid2D_Quad_Block_HO &Grid_Original,
			  const int Sector){
    return Grid_Fine.Refine_Mesh(Grid_Original,Sector);
  }
  
  void Coarsen_Mesh(const Grid2D_Quad_Block_HO &Grid_Original_SW,
		    const Grid2D_Quad_Block_HO &Grid_Original_SE,
		    const Grid2D_Quad_Block_HO &Grid_Original_NW,
		    const Grid2D_Quad_Block_HO &Grid_Original_NE);
  //! Generate a coarsen mesh (Grid_Coarse) with the 4 sector grids (SW,SE,NW,NE)
  friend void Coarsen_Mesh(Grid2D_Quad_Block_HO &Grid_Coarse,
			   const Grid2D_Quad_Block_HO &Grid_Original_SW,
			   const Grid2D_Quad_Block_HO &Grid_Original_SE,
			   const Grid2D_Quad_Block_HO &Grid_Original_NW,
			   const Grid2D_Quad_Block_HO &Grid_Original_NE){
    return Grid_Coarse.Coarsen_Mesh(Grid_Original_SW,Grid_Original_SE,Grid_Original_NW,Grid_Original_NE);
  }
  
  void Fix_Refined_Mesh_Boundaries(const int Fix_North_Boundary,
				   const int Fix_South_Boundary,
				   const int Fix_East_Boundary,
				   const int Fix_West_Boundary);
  
  void Unfix_Refined_Mesh_Boundaries(void);
  //@}

  //! @name Measurement functions.
  //@{
  double MinimumNodeEdgeDistance(const int &i, const int &j) const;
  double DistanceFromPointToLine(const Vector2D &Point,
				 const Vector2D &LinePoint1,
				 const Vector2D &LinePoint2) const;
  //@}
  
  //! @name Binary arithmetic operators.
  //@{
  void operator +(const Vector2D &V); // Translate grid with +V
  void operator -(const Vector2D &V); // Translate grid with -V
  void operator *(const double &a);   // Scale grid
  void operator /(const double &a){ return *this * (1.0/a); }  // Scale grid with 1/a
  void operator ^(const double &a);   // Rotate grid
  //@}

  //! @name Friend binary arithmetic operators.
  //@{
  friend void operator *(const double &a, Grid2D_Quad_Block_HO &G) { return G*a;}
  //@}

  //! @name Geometric boundary representation related functions.
  //@{
  //! Set the designated switch to require high-order representation of the cell geometric boundaries
  static void setHighOrderBoundaryRepresentation(void){ HighOrderBoundaryRepresentation = ON; }
  //! Set the designated switch to require treatment of cell geometric boundaries as straight lines
  static void setLowOrderBoundaryRepresentation(void){ HighOrderBoundaryRepresentation = OFF; }
  //! Set the boundary representation designated switch to the default value (i.e. low-order representation)
  static void setDefaultBoundaryRepresentation(void){ HighOrderBoundaryRepresentation = OFF; }
  //! Return true if geometric boundary representation is high-order otherwise return false.
  static bool IsHighOrderBoundary(void) { return HighOrderBoundaryRepresentation == ON? true:false; }
  //! Get the value of the HighOrderBoundaryRepresentation variable.
  static int getHighOrderBoundaryValue(void) {return HighOrderBoundaryRepresentation; }
  //! Set the designated switch to require the use of Gauss quadratures for evaluating curvilinear path integrals.
  static void setContourIntegrationBasedOnGaussQuadratures(void) {
    Gauss_Quad_Curvilinear_Integration = ON; Mixed_Curvilinear_Integration = OFF;
  }
  //! Set the designated switch to require the use of summation on line segments to approximate curvilinear path integrals.
  static void setContourIntegrationBasedOnLinearSegments(void) {
    Gauss_Quad_Curvilinear_Integration = OFF; Mixed_Curvilinear_Integration = OFF;
  }
  //! Set the designated switch to require the use of line segments to path integrals in area and centroid calculation.
  static void setMixedContourIntegration(void) {
    Gauss_Quad_Curvilinear_Integration = ON; Mixed_Curvilinear_Integration = ON;
  }
  //! Set Monte Carlo integration for cells with curved edges
  static void setMonteCarloIntegrationON(void){ Monte_Carlo_Integration_Allowed = ON; }
  //! Turn off Monte Carlo integration for cells with curved edges
  static void setMonteCarloIntegrationOFF(void){ Monte_Carlo_Integration_Allowed = OFF; }
  //! Set polygonal adaptive quadrature integration for cells with curved edges
  static void setPolygonalAdaptiveQuadratureIntegrationON(void){  Polygonal_Adaptive_Quadrature_Integration_Allowed = ON; }
  //! Turn off polygonal adaptive quadrature integration for cells with curved edges
  static void setPolygonalAdaptiveQuadratureIntegrationOFF(void){ Polygonal_Adaptive_Quadrature_Integration_Allowed = OFF; }
  //! Turn ON flag for non-reflected ghost cells
  static void setNonReflectedGhostCellsNearSouthSolidBoundary(void) { Mesh_Requiring_NonReflected_South_Ghost_Cells = ON;}
  //! Turn OFF flag for non-reflected ghost cells
  static void setReflectedGhostCellsNearSouthSolidBoundary(void) { Mesh_Requiring_NonReflected_South_Ghost_Cells = OFF;}

  //! Check if the West boundary has constrained reconstruction (i.e. it is curved and set to reconstruction based flux)
  bool IsWestBoundaryReconstructionConstrained(void) const;
  //! Check if the East boundary has constrained reconstruction (i.e. it is curved and set to reconstruction based flux)
  bool IsEastBoundaryReconstructionConstrained(void) const;
  //! Check if the South boundary has constrained reconstruction (i.e. it is curved and set to reconstruction based flux)
  bool IsSouthBoundaryReconstructionConstrained(void) const;
  //! Check if the North boundary has constrained reconstruction (i.e. it is curved and set to reconstruction based flux)
  bool IsNorthBoundaryReconstructionConstrained(void) const;

  //! Check if the North extension of West boundary has constrained reconstruction (i.e. same conditions as for the West boundary)
  bool IsNorthExtendWestBoundaryReconstructionConstrained(void) const;
  //! Check if the South extension of West boundary has constrained reconstruction (i.e. same conditions as for the West boundary)
  bool IsSouthExtendWestBoundaryReconstructionConstrained(void) const;
  //! Check if the North extension of East boundary has constrained reconstruction (i.e. same conditions as for the East boundary)
  bool IsNorthExtendEastBoundaryReconstructionConstrained(void) const;
  //! Check if the South extension of East boundary has constrained reconstruction (i.e. same conditions as for the East boundary)
  bool IsSouthExtendEastBoundaryReconstructionConstrained(void) const;
  //! Check if the East extension of South boundary has constrained reconstruction (i.e. same conditions as for the South boundary)
  bool IsEastExtendSouthBoundaryReconstructionConstrained(void) const;
  //! Check if the West extension of South boundary has constrained reconstruction (i.e. same conditions as for the South boundary)
  bool IsWestExtendSouthBoundaryReconstructionConstrained(void) const;
  //! Check if the East extension of North boundary has constrained reconstruction (i.e. same conditions as for the North boundary)
  bool IsEastExtendNorthBoundaryReconstructionConstrained(void) const;
  //! Check if the West extension of North boundary has constrained reconstruction (i.e. same conditions as for the North boundary)
  bool IsWestExtendNorthBoundaryReconstructionConstrained(void) const;

  
  //! Check if any block boundary is a solid body
  bool IsThereAnySolidBoundary(void) const;

  //! Check if the West boundary is curved (i.e. has control points and is not an interior boundary based on the BCs)
  bool IsWestBoundaryCurved(void) const;
  //! Check if the East boundary is curved (i.e. has control points and is not an interior boundary based on the BCs)
  bool IsEastBoundaryCurved(void) const;
  //! Check if the South boundary is curved (i.e. has control points and is not an interior boundary based on the BCs)
  bool IsSouthBoundaryCurved(void) const;
  //! Check if the North boundary is curved (i.e. has control points and is not an interior boundary based on the BCs)
  bool IsNorthBoundaryCurved(void) const;

  //! Check if the North extension of the West boundary is curved (i.e. same conditions as for the West boundary)
  bool IsNorthExtendWestBoundaryCurved(void) const;
  //! Check if the South extension of the West boundary is curved (i.e. same conditions as for the West boundary)
  bool IsSouthExtendWestBoundaryCurved(void) const;
  //! Check if the North extension of the East boundary is curved (i.e. same conditions as for the East boundary)
  bool IsNorthExtendEastBoundaryCurved(void) const;
  //! Check if the South extension of the East boundary is curved (i.e. same conditions as for the East boundary)
  bool IsSouthExtendEastBoundaryCurved(void) const;
  //! Check if the East extension of the South boundary is curved (i.e. same conditions as for the South boundary)
  bool IsEastExtendSouthBoundaryCurved(void) const;
  //! Check if the West extension of the South boundary is curved (i.e. same conditions as for the South boundary)
  bool IsWestExtendSouthBoundaryCurved(void) const;
  //! Check if the East extension of the North boundary is curved (i.e. same conditions as for the North boundary)
  bool IsEastExtendNorthBoundaryCurved(void) const;
  //! Check if the West extension of the North boundary is curved (i.e. same conditions as for the North boundary)
  bool IsWestExtendNorthBoundaryCurved(void) const;
  //@}

  //! Set the designated switch to require computation of geometric properties with extra care for numerical errors
  static void setTreatMeshWithExtraCareForNumericalError(void) { Minimize_Error_Calculation_Of_Geometric_Properties = ON; }
  //! Set the designated switch to require the regular computation of geometric properties.
  static void setNoSpecialTreatmentForNumericalError(void) { Minimize_Error_Calculation_Of_Geometric_Properties = OFF; }

  //! Set the designated switch to accept geometric inaccuracies near curved boundaries if necessary to obtain a result
  static void setTolerateInaccurateIntegrationNearCurvedBoundaries(void) {
    Tolerate_Inaccurate_Integration_Near_Curved_Boundaries = ON;
  }
  //! Set the designated switch to be intolerant to geometric inaccuracies near curved boundaries
  static void setNoGeometricInaccuraciesForIntegrationNearCurvedBoundaries(void){
    Tolerate_Inaccurate_Integration_Near_Curved_Boundaries = OFF;
  }
  //! Check if geometric inaccuracies for integration along curved boundaries are tolerated
  bool IsIntegrationAlongCurvedBoundariesToleratedToGeometricInaccuracies(void) const{
    return Tolerate_Inaccurate_Integration_Near_Curved_Boundaries == ON? true:false;
  }

  //! Check if Monte Carlo integration in cells with curved boundaries is allowed
  bool IsMonteCarloIntegrationAllowed(void) const { 
    return Monte_Carlo_Integration_Allowed == ON? true:false;
  }

  //! Check if adaptive polygonal integration in cells with curved boundaries is allowed
  bool IsPolygonalAdaptiveQuadraturesAllowed(void) const {
    return Polygonal_Adaptive_Quadrature_Integration_Allowed == ON? true:false;
  }

  //! Output only the cell data
  void Output_Cells_Data(const int &block_number, ostream &out_file);
  //! Output only the cell invariant geometric properties
  void Output_Cells_Translation_Rotation_Invariant_Properties(const int &block_number, ostream &out_file);

  //! @name Input-output operators.
  //@{
  friend ostream &operator << (ostream &out_file, const Grid2D_Quad_Block_HO &G);
  friend istream &operator >> (istream &in_file, Grid2D_Quad_Block_HO &G);
  //@}

  //! Parameter to set the minimum levels of refinement during the polygonal adaptive quadrature integration
  static short Polygonal_Adaptive_Quadrature_Integration_Minimum_Levels;

private:
  Grid2D_Quad_Block_HO(const Grid2D_Quad_Block_HO &G);     //! Private copy constructor.
  
  //! Switch for computing with high-order or low-order accuracy the geometric properties 
  static int HighOrderBoundaryRepresentation;

  //! Switch for applying or not the smoothing subroutine
  static int Smooth_Quad_Block_Flag;

  //! Switch for how to compute the curvilinear path integrals along curved edges.
  static int Gauss_Quad_Curvilinear_Integration;

  //! Switch for how to compute the curvilinear path integrals along curved edges.
  static int Mixed_Curvilinear_Integration;
  
  //! Switch for error minimization in calculating the cell geometric properties.
  static int Minimize_Error_Calculation_Of_Geometric_Properties;

  //! Switch for treatment of integrals along curved boundaries, in case they cannot be performed accurately along the real geometry.
  static int Tolerate_Inaccurate_Integration_Near_Curved_Boundaries;

  //! Switch for Monte Carlo integration in curved boundaries (i.e. turn it OFF or ON). 
  static int Monte_Carlo_Integration_Allowed;

  /*! Switch for polygonal adaptive quadrature integration in curved boundaries (i.e. turn it OFF or ON).
   *  It has higher priority than Monte_Carlo_Integration_Allowed */
  static int Polygonal_Adaptive_Quadrature_Integration_Allowed;

  /*! Switch to force non-reflection of ghost cells regardless of the boundary conditions
   *  in order to generate a valid mesh (i.e. not crossed quadrilaterals in ghost cells).
   */
  static int Mesh_Requiring_NonReflected_South_Ghost_Cells;

  //! Highest order of reconstruction that might occur in calculations with the current grid.
  int HighestReconstructionOrder;

  //! @name Flags to define different levels of mesh update.
  //@{
  //! Controls the update of the geometric properties of the interior cells.
  int InteriorMeshUpdate;
  //! Controls the update of the geometric properties of all ghost cells.
  int GhostCellsUpdate;    
  //! Controls the update of the geometric properties of the corner ghost cells.
  int CornerGhostCellsUpdate; 
  //! Reset to NO update all the mesh update flags.
  void Reset_Mesh_Update_Flags(void){ InteriorMeshUpdate = OFF; GhostCellsUpdate = OFF; CornerGhostCellsUpdate = OFF;}

  //! State tracker for the geometry of interior cells (i.e. indentifies uniquely the state of the interior cell geometry)
  unsigned int InteriorCellGeometryStateTracker;
  //! State tracker for the geometry of all ghost cells (i.e. indentifies uniquely the state of the ghost cell geometry)
  unsigned int GhostCellGeometryStateTracker;
  //! State tracker for the geometry of corner ghost cells (i.e. indentifies uniquely the state of the corner ghost cell geometry)
  unsigned int CornerGhostCellGeometryStateTracker;

  //! @brief Set trackers to identify a new state everywhere (i.e. interior and ghost cells)
  void New_Global_Geometry_State(void);
  //! @brief Set tracker to identify a new interior state
  void New_Interior_Geometry_State(void);
  //! @brief Set tracker to identify a new state for the ghost cells
  void New_Ghost_Geometry_State(void);
  //! @brief Set tracker to identify a new state for the corner ghost cells
  void New_Corner_Geometry_State(void);
  //@}

  /*!
   * Number of Gauss quadrature points require for flux calculation on interior edges or low-order exterior boundaries.
   */
  int NumGQP;	

  //! Check for curved boundaries.
  bool CheckExistenceOfCurvedBoundaries(void);

  //! @name Auxiliary variables/functions for translations between the global and the cell local coordinate systems.
  //@{
  static Vector2D _SW_, _SE_, _NE_, _NW_;
  void TranslateVertexesIntoLocalCoordinateSystem(Vector2D &SW, Vector2D &SE,
						  Vector2D &NE, Vector2D &NW);
  void TranslateVertexesIntoGlobalCoordinateSystem(Vector2D &SW, Vector2D &SE,
						   Vector2D &NE, Vector2D &NW,
						   Vector2D &Centroid);
  //@}

  //! Update the interval infos associated to extension splines
  void Update_ExtensionSplineInfos(void);

  //! Update ghost cells in the proximity of curved extension splines
  void Update_GhostCells_Near_CurvedExtensionSplines(const bool &CurvedNorthBnd,
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
						     const bool &Curved_Extend_W_NorthBnd);

  int NumOfConstrainedGaussQuadPoints_GenericBoundary(const BndSplineType & BndSpline,
						      const BndSplineIntervalType * BndSplineInfo,
						      const int &CellIndex, 
						      const int &CellIndex_Shift = 0);

};

/*!
 * Default constructor.
 */
inline Grid2D_Quad_Block_HO::Grid2D_Quad_Block_HO(void)
  : Integration(this),
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
    StretchI(0), StretchJ(0), BetaI(ONE), TauI(ONE),
    BetaJ(ONE), TauJ(ONE),
    OrthogonalN(1), OrthogonalS(1), OrthogonalE(1), OrthogonalW(1),
    // Initialize mesh update flags to OFF (i.e. no update scheduled)
    InteriorMeshUpdate(OFF), GhostCellsUpdate(OFF), CornerGhostCellsUpdate(OFF),
    // Initialize state trackers
    InteriorCellGeometryStateTracker(0), GhostCellGeometryStateTracker(0), CornerGhostCellGeometryStateTracker(0),
    NumGQP(0)
{
  // 
}

/*!
 * Allocate memory.
 * This is a short version of the main allocation subroutine.
 * To get high-order one should always use the main allocation subroutine.
 *
 * \param Ni number of cells in i-direction
 * \param Nj number of cells in j-direction
 * \param Ng number of ghost cells
 */
inline void Grid2D_Quad_Block_HO::allocate(const int &Ni, const int &Nj, const int &Ng) {
  return allocate(Ni,Nj,Ng,0);
}

/*!
 * Deallocate spline info memory.
 */
inline void Grid2D_Quad_Block_HO::deallocateBndSplineInfo(void) {
  delete []BndNorthSplineInfo; BndNorthSplineInfo = NULL;
  delete []BndSouthSplineInfo; BndSouthSplineInfo = NULL;
  delete []BndEastSplineInfo;  BndEastSplineInfo = NULL;
  delete []BndWestSplineInfo;  BndWestSplineInfo = NULL;

  delete []ExtendWest_BndNorthSplineInfo; ExtendWest_BndNorthSplineInfo = NULL; 
  delete []ExtendEast_BndNorthSplineInfo; ExtendEast_BndNorthSplineInfo = NULL;
  delete []ExtendWest_BndSouthSplineInfo; ExtendWest_BndSouthSplineInfo = NULL; 
  delete []ExtendEast_BndSouthSplineInfo; ExtendEast_BndSouthSplineInfo = NULL;
  delete []ExtendNorth_BndEastSplineInfo; ExtendNorth_BndEastSplineInfo = NULL;
  delete []ExtendSouth_BndEastSplineInfo; ExtendSouth_BndEastSplineInfo = NULL;
  delete []ExtendNorth_BndWestSplineInfo; ExtendNorth_BndWestSplineInfo = NULL;
  delete []ExtendSouth_BndWestSplineInfo; ExtendSouth_BndWestSplineInfo = NULL;
}

/*!
 * Remove all extension splines and associated Infos.
 */
inline void Grid2D_Quad_Block_HO::RemoveExtensionSplines(void){
  
  // Deallocate extension splines
  ExtendWest_BndNorthSpline.deallocate();
  ExtendEast_BndNorthSpline.deallocate();
  ExtendWest_BndSouthSpline.deallocate();
  ExtendEast_BndSouthSpline.deallocate();
  ExtendNorth_BndEastSpline.deallocate();
  ExtendSouth_BndEastSpline.deallocate();
  ExtendNorth_BndWestSpline.deallocate();
  ExtendSouth_BndWestSpline.deallocate();

  // Deallocate associated infos
  delete []ExtendWest_BndNorthSplineInfo; ExtendWest_BndNorthSplineInfo = NULL; 
  delete []ExtendEast_BndNorthSplineInfo; ExtendEast_BndNorthSplineInfo = NULL;
  delete []ExtendWest_BndSouthSplineInfo; ExtendWest_BndSouthSplineInfo = NULL; 
  delete []ExtendEast_BndSouthSplineInfo; ExtendEast_BndSouthSplineInfo = NULL;
  delete []ExtendNorth_BndEastSplineInfo; ExtendNorth_BndEastSplineInfo = NULL;
  delete []ExtendSouth_BndEastSplineInfo; ExtendSouth_BndEastSplineInfo = NULL;
  delete []ExtendNorth_BndWestSplineInfo; ExtendNorth_BndWestSplineInfo = NULL;
  delete []ExtendSouth_BndWestSplineInfo; ExtendSouth_BndWestSplineInfo = NULL;  
}

/*!
 * Subroutine for translating the providing vectors
 * into a reference system with the origin in SW.
 * The original values are stored in the static 
 * variables _SW_, _SE_, _NE_, _NW_.
 */
inline void Grid2D_Quad_Block_HO::TranslateVertexesIntoLocalCoordinateSystem(Vector2D &SW, Vector2D &SE,
									     Vector2D &NE, Vector2D &NW){

  // Check if translation to local coordinate system is required
  if (Minimize_Error_Calculation_Of_Geometric_Properties){
    // save the current node locations 
   _SW_ = SW;
   _SE_ = SE;
   _NE_ = NE;
   _NW_ = NW;

    // translate the vertexes into a local reference system with the origin in SW.
    SW = 0.0;
    SE -= _SW_;
    NE -= _SW_;
    NW -= _SW_;
  }
}

/*!
 * Subroutine for translating the providing vectors
 * into the global reference system.
 * The function assumes that _SW_, _SE_, _NE_, _NW_
 * represent the original values of the providing vectors.
 * This function should be used in conjunction with 
 * "TranslateVertexesIntoLocalCoordinateSystem()".
 */
inline void Grid2D_Quad_Block_HO::TranslateVertexesIntoGlobalCoordinateSystem(Vector2D &SW, Vector2D &SE,
									      Vector2D &NE, Vector2D &NW,
									      Vector2D &Centroid){

  if (Minimize_Error_Calculation_Of_Geometric_Properties){
    SW = _SW_;
    SE = _SE_;
    NE = _NE_;
    NW = _NW_;
    Centroid += _SW_;
  }
}

/*!
 * Calculate the centroid and the area of a quadrilateral 
 * cell (ii,jj) and store them in the designated variables.
 * This subroutine gives correct results for both concave 
 * and convex quadrilaterals.
 * However, the subroutine doesn't check for successful
 * execution (i.e. occurrence of negative area)!
 * 
 * \return 0 for normal execution;
 *         2 if area = zero (i.e. the centroid is undefined).
 */
inline int Grid2D_Quad_Block_HO::quadAreaAndCentroid(const int &ii, const int &jj) const {

  // Cell nodes in counterclockwise order (SW, SE, NE, NW).
  
  return quadCentroid(Node[ii  ][jj  ].X, Node[ii+1][jj  ].X,
		      Node[ii+1][jj+1].X, Node[ii  ][jj+1].X,
		      Cell[ii][jj].Xc, Cell[ii][jj].A);
}

/*!
 * Alternative approach to get the centroid of a convex quadrilateral.
 * If the quadrilateral is concave this subroutine gives an incorrect answer!
 * However, this approach is fast for convex quadrilateral.
 */
inline Vector2D Grid2D_Quad_Block_HO::quadConvexCentroid(const int &ii, const int &jj) const{

  Vector2D X1, X2, X3, X4, Xc1, Xc2, X;
  double A1, A2;
  // Cell nodes in counter-clockwise order.
  X1 = Node[ii  ][jj  ].X;
  X2 = Node[ii+1][jj  ].X;
  X3 = Node[ii+1][jj+1].X;
  X4 = Node[ii  ][jj+1].X;
  // Determine the centroid of the sub-triangles.
  Xc1 = (X1+X2+X3)/3.0;
  Xc2 = (X1+X3+X4)/3.0;

  // Determine the area of the sub-triangles.
  //   A1 = HALF*((X1^X2) + (X2^X3) + (X3^X1));
  //   A2 = HALF*((X1^X3) + (X3^X4) + (X4^X1));

  // These relationships are equivalent with the ones shown above.
  A1 = HALF*((X2-X1)^(X3-X1));
  A2 = HALF*((X3-X4)^(X3-X1));

  // Return the area-weighted average of the centroids of the sub-triangles:
  if (A1 > ZERO && A2 > ZERO) return (A1*Xc1 + A2*Xc2)/(A1+A2);
  // Average of four nodes (not always correct):
  return 0.25*(Node[ii][jj].X + Node[ii+1][jj].X + Node[ii+1][jj+1].X + Node[ii][jj+1].X);
}

/*!
 * Get centroid of Cell
 */
inline Vector2D Grid2D_Quad_Block_HO::centroid(const Cell2D_HO &Cell) {
  return centroid(Cell.I,Cell.J);
}

/*!
 * Calculate the centroid of cell (ii,jj).
 * This is a slower subroutine because it checks for negative area.
 */
inline Vector2D Grid2D_Quad_Block_HO::centroid(const int &ii, const int &jj) {

  int Info;
  double area;
  Vector2D Centroid;
  
  // Cell nodes in counterclockwise order (SW, SE, NE, NW).
  Vector2D X[4] = { Node[ii  ][jj  ].X,
		    Node[ii+1][jj  ].X ,
		    Node[ii+1][jj+1].X,
		    Node[ii  ][jj+1].X };
  
  // Translate the cell nodes into the local coordinate system if required
  TranslateVertexesIntoLocalCoordinateSystem(X[0],X[1],X[2],X[3]);

  Info = quadCentroid(X[0],X[1],X[2],X[3],Centroid,area);

  // Translate the cell nodes back to the global coordinate system.
  TranslateVertexesIntoGlobalCoordinateSystem(X[0],X[1],X[2],X[3],Centroid);
  
  if (Info == 2){
    throw runtime_error("Grid2D_Quad_Block_HO::centroid() ERROR! Negative area encountered!");
  }
  
  // Return the centroid
  return Centroid;
}

/*!
 * Calculate the centroid of the South-West quarter of cell (ii,jj)
 */
inline Vector2D Grid2D_Quad_Block_HO::centroidSW(const int &ii, const int &jj) {

  Vector2D X1,X2,X3,X4, Centroid;
  double area;

  // Cell nodes in counter-clockwise order.
  X1 = Node[ii][jj].X;
  X2 = HALF*(Node[ii][jj].X + Node[ii+1][jj].X);
  X3 = 0.25*(Node[ii][jj].X + Node[ii+1][jj].X + Node[ii][jj+1].X + Node[ii+1][jj+1].X);
  X4 = HALF*(Node[ii][jj].X + Node[ii][jj+1].X);

  // Translate the cell nodes into the local coordinate system if required
  TranslateVertexesIntoLocalCoordinateSystem(X1,X2,X3,X4);

  // Determine the centroid
  quadCentroid(X1,X2,X3,X4,Centroid,area);

  // Translate the cell nodes back to the global coordinate system.
  TranslateVertexesIntoGlobalCoordinateSystem(X1,X2,X3,X4,Centroid);

  return Centroid;
}

/*!
 * Alternative approach to calculate the centroid of the South-West 
 * quarter of cell (ii,jj), when the formed quadrilateral is convex.
 */
inline Vector2D Grid2D_Quad_Block_HO::centroidSW_ConvexQuad(const int &ii, const int &jj) const {
  Vector2D X1, X2, X3, X4, Xc1, Xc2;
  double A1, A2;
  // Cell nodes in counter-clockwise order.
  X1 = Node[ii][jj].X;
  X2 = HALF*(Node[ii][jj].X + Node[ii+1][jj].X);
  X3 = 0.25*(Node[ii][jj].X + Node[ii+1][jj].X + Node[ii][jj+1].X + Node[ii+1][jj+1].X);
  X4 = HALF*(Node[ii][jj].X + Node[ii][jj+1].X);
  // Determine the centroid and area of the sub-triangles.
  Xc1 = (X1+X2+X3)/3.0;
  Xc2 = (X1+X3+X4)/3.0;
  //   A1 = HALF*((X1^X2) + (X2^X3) + (X3^X1));
  //   A2 = HALF*((X1^X3) + (X3^X4) + (X4^X1));
  A1 = HALF*((X2-X1)^(X3-X1));
  A2 = HALF*((X3-X4)^(X3-X1));
  // Return the area-weighted average of the centroids of the sub-triangles:
  if (A1 > ZERO && A2 > ZERO) return (A1*Xc1 + A2*Xc2)/(A1+A2);
  // Average of four nodes (not always correct):
  return 0.25*(X1 + X2 + X3 + X4);
}

/*!
 * Calculate the centroid of the South-East quarter of cell (ii,jj)
 */
inline Vector2D Grid2D_Quad_Block_HO::centroidSE(const int &ii, const int &jj) {
  Vector2D X1, X2, X3, X4, Centroid;
  double area;

  // Cell nodes in counter-clockwise order.
  X1 = HALF*(Node[ii][jj].X + Node[ii+1][jj].X);
  X2 = Node[ii+1][jj].X;
  X3 = HALF*(Node[ii+1][jj].X + Node[ii+1][jj+1].X);
  X4 = 0.25*(Node[ii][jj].X + Node[ii+1][jj].X + Node[ii][jj+1].X + Node[ii+1][jj+1].X);

  // Translate the cell nodes into the local coordinate system if required
  TranslateVertexesIntoLocalCoordinateSystem(X1,X2,X3,X4);

  // Determine the centroid
  quadCentroid(X1,X2,X3,X4,Centroid,area);

  // Translate the cell nodes back to the global coordinate system.
  TranslateVertexesIntoGlobalCoordinateSystem(X1,X2,X3,X4,Centroid);

  return Centroid;
}

/*!
 * Calculate the centroid of the South-East quarter of cell (ii,jj)
 */
inline Vector2D Grid2D_Quad_Block_HO::centroidSE_ConvexQuad(const int &ii, const int &jj) const {
  Vector2D X1, X2, X3, X4, Xc1, Xc2;
  double A1, A2;
  // Cell nodes in counter-clockwise order.
  X1 = HALF*(Node[ii][jj].X + Node[ii+1][jj].X);
  X2 = Node[ii+1][jj].X;
  X3 = HALF*(Node[ii+1][jj].X + Node[ii+1][jj+1].X);
  X4 = 0.25*(Node[ii][jj].X + Node[ii+1][jj].X + Node[ii][jj+1].X + Node[ii+1][jj+1].X);
  // Determine the centroid and area of the sub-triangles.
  Xc1 = (X1+X2+X3)/3.0;
  Xc2 = (X1+X3+X4)/3.0;
  //   A1 = HALF*((X1^X2) + (X2^X3) + (X3^X1));
  //   A2 = HALF*((X1^X3) + (X3^X4) + (X4^X1));
  A1 = HALF*((X2-X1)^(X3-X1));
  A2 = HALF*((X3-X4)^(X3-X1));
  // Return the area-weighted average of the centroids of the sub-triangles:
  if (A1 > ZERO && A2 > ZERO) return (A1*Xc1 + A2*Xc2)/(A1+A2);
  // Average of four nodes (not always correct):
  return 0.25*(X1 + X2 + X3 + X4);
}

/*!
 * Calculate the centroid of the North-West quarter of cell (ii,jj)
 */
inline Vector2D Grid2D_Quad_Block_HO::centroidNW(const int &ii, const int &jj) {

  Vector2D X1, X2, X3, X4, Centroid;
  double area;

  // Cell nodes in counter-clockwise order.
  X1 = HALF*(Node[ii][jj].X + Node[ii][jj+1].X);
  X2 = 0.25*(Node[ii][jj].X + Node[ii+1][jj].X + Node[ii][jj+1].X + Node[ii+1][jj+1].X);
  X3 = HALF*(Node[ii][jj+1].X + Node[ii+1][jj+1].X);
  X4 = Node[ii][jj+1].X;

  // Translate the cell nodes into the local coordinate system if required
  TranslateVertexesIntoLocalCoordinateSystem(X1,X2,X3,X4);

  // Determine the centroid.
  quadCentroid(X1,X2,X3,X4,Centroid,area);

  // Translate the cell nodes back to the global coordinate system.
  TranslateVertexesIntoGlobalCoordinateSystem(X1,X2,X3,X4,Centroid);

  return Centroid;
}

/*!
 * Calculate the centroid of the North-West quarter of cell (ii,jj)
 */
inline Vector2D Grid2D_Quad_Block_HO::centroidNW_ConvexQuad(const int &ii, const int &jj) const {
  Vector2D X1, X2, X3, X4, Xc1, Xc2;
  double A1, A2;
  // Cell nodes in counter-clockwise order.
  X1 = HALF*(Node[ii][jj].X + Node[ii][jj+1].X);
  X2 = 0.25*(Node[ii][jj].X + Node[ii+1][jj].X + Node[ii][jj+1].X + Node[ii+1][jj+1].X);
  X3 = HALF*(Node[ii][jj+1].X + Node[ii+1][jj+1].X);
  X4 = Node[ii][jj+1].X;
  // Determine the centroid and area of the sub-triangles.
  Xc1 = (X1+X2+X3)/3.0;
  Xc2 = (X1+X3+X4)/3.0;
  //   A1 = HALF*((X1^X2) + (X2^X3) + (X3^X1));
  //   A2 = HALF*((X1^X3) + (X3^X4) + (X4^X1));
  A1 = HALF*((X2-X1)^(X3-X1));
  A2 = HALF*((X3-X4)^(X3-X1));
  // Return the area-weighted average of the centroids of the sub-triangles:
  if (A1 > ZERO && A2 > ZERO) return (A1*Xc1 + A2*Xc2)/(A1+A2);
  // Average of four nodes (not always correct):
  return 0.25*(X1 + X2 + X3 + X4);
}

/*!
 * Calculate the centroid of the North-East quarter of cell (ii,jj)
 */
inline Vector2D Grid2D_Quad_Block_HO::centroidNE(const int &ii, const int &jj) {

  Vector2D X1, X2, X3, X4, Centroid;
  double area;

  // Cell nodes in counter-clockwise order.
  X1 = 0.25*(Node[ii][jj].X + Node[ii+1][jj].X + Node[ii][jj+1].X + Node[ii+1][jj+1].X);
  X2 = HALF*(Node[ii+1][jj].X + Node[ii+1][jj+1].X);
  X3 = Node[ii+1][jj+1].X;
  X4 = HALF*(Node[ii][jj+1].X + Node[ii+1][jj+1].X);

  // Translate the cell nodes into the local coordinate system if required
  TranslateVertexesIntoLocalCoordinateSystem(X1,X2,X3,X4);

  // Determine the centroid.
  quadCentroid(X1,X2,X3,X4,Centroid,area);

  // Translate the cell nodes back to the global coordinate system.
  TranslateVertexesIntoGlobalCoordinateSystem(X1,X2,X3,X4,Centroid);

  return Centroid;
}

/*!
 * Calculate the centroid of the North-East quarter of cell (ii,jj)
 */
inline Vector2D Grid2D_Quad_Block_HO::centroidNE_ConvexQuad(const int &ii, const int &jj) const {
  Vector2D X1, X2, X3, X4, Xc1, Xc2;
  double A1, A2;
  // Cell nodes in counter-clockwise order.
  X1 = 0.25*(Node[ii][jj].X + Node[ii+1][jj].X + Node[ii][jj+1].X + Node[ii+1][jj+1].X);
  X2 = HALF*(Node[ii+1][jj].X + Node[ii+1][jj+1].X);
  X3 = Node[ii+1][jj+1].X;
  X4 = HALF*(Node[ii][jj+1].X + Node[ii+1][jj+1].X);
  // Determine the centroid and area of the sub-triangles.
  Xc1 = (X1+X2+X3)/3.0;
  Xc2 = (X1+X3+X4)/3.0;
  //   A1 = HALF*((X1^X2) + (X2^X3) + (X3^X1));
  //   A2 = HALF*((X1^X3) + (X3^X4) + (X4^X1));
  A1 = HALF*((X2-X1)^(X3-X1));
  A2 = HALF*((X3-X4)^(X3-X1));
  // Return the area-weighted average of the centroids of the sub-triangles:
  if (A1 > ZERO && A2 > ZERO) return (A1*Xc1 + A2*Xc2)/(A1+A2);
  // Average of four nodes (not always correct):
  return (X1 + X2 + X3 + X4)/4.0;
}

/*!
 * Calculate the cell area of a quadrilateral cell with straight edges
 */
inline double Grid2D_Quad_Block_HO::area(const Cell2D_HO &Cell) const {
  return area(Cell.I,Cell.J);
}

/*!
 * Calculate the cell area of a quadrilateral cell with straight edges
 *
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 */
inline double Grid2D_Quad_Block_HO::area(const int &ii, const int &jj) const {
  return HALF*( ((Node[ii+1][jj].X-Node[ii][jj].X)^(Node[ii][jj+1].X-Node[ii][jj].X)) +
		((Node[ii+1][jj+1].X-Node[ii][jj+1].X)^(Node[ii+1][jj+1].X-Node[ii+1][jj].X)) );
}

/*!
 * Calculate the cell area of a quadrilateral cell with straight edges 
 * defined by 4 vertexes (i.e. nodes SW, NW, NE, SE. )
 * 
 * \note The order in which the nodes are passed is important!
 */
inline double Grid2D_Quad_Block_HO::area(const Node2D_HO &SW, const Node2D_HO &NW,
					 const Node2D_HO &NE, const Node2D_HO &SE){
  return (HALF*(((SE.X-SW.X)^( NW.X-SW.X)) + ((NE.X-NW.X)^(NE.X-SE.X))));
}


/*!
 * Get North-West cell node.
 */
inline Node2D_HO Grid2D_Quad_Block_HO::nodeNW(const Cell2D_HO &Cell) const {
  return (Node[Cell.I][Cell.J+1]);
}

/*!
 * Get North-West cell node.
 */
inline Node2D_HO Grid2D_Quad_Block_HO::nodeNW(const int ii, const int jj) const {
  return (Node[ii][jj+1]);
}

/*!
 * Get North-East cell node.
 */
inline Node2D_HO Grid2D_Quad_Block_HO::nodeNE(const Cell2D_HO &Cell) const {
  return (Node[Cell.I+1][Cell.J+1]);
}

/*!
 * Get North-East cell node.
 */
inline Node2D_HO Grid2D_Quad_Block_HO::nodeNE(const int ii, const int jj) const {
  return (Node[ii+1][jj+1]);
}

/*!
 * Get South-East cell node.
 */
inline Node2D_HO Grid2D_Quad_Block_HO::nodeSE(const Cell2D_HO &Cell) const {
  return (Node[Cell.I+1][Cell.J]);
}

/*!
 * Get South-East cell node.
 */
inline Node2D_HO Grid2D_Quad_Block_HO::nodeSE(const int ii, const int jj) const {
  return (Node[ii+1][jj]);
}

/*!
 * Get South-West cell node.
 */
inline Node2D_HO Grid2D_Quad_Block_HO::nodeSW(const Cell2D_HO &Cell) const {
  return (Node[Cell.I][Cell.J]);
}

/*!
 * Get South-West cell node.
 */
inline Node2D_HO Grid2D_Quad_Block_HO::nodeSW(const int ii, const int jj) const {
  return (Node[ii][jj]);
}

/*!
 * Get North face midpoint.
 */
inline Vector2D Grid2D_Quad_Block_HO::xfaceN(const Cell2D_HO &Cell) const {
  return xfaceN(Cell.I,Cell.J);
}

/*!
 * Get North face midpoint.
 *
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 */
inline Vector2D Grid2D_Quad_Block_HO::xfaceN(const int ii, const int jj) const {
  return (HALF*(Node[ii][jj+1].X+Node[ii+1][jj+1].X));
}

/*!
 * Get South face midpoint.
 */
inline Vector2D Grid2D_Quad_Block_HO::xfaceS(const Cell2D_HO &Cell) const {
  return xfaceS(Cell.I,Cell.J);
}

/*!
 * Get South face midpoint.
 *
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 */
inline Vector2D Grid2D_Quad_Block_HO::xfaceS(const int ii, const int jj) const {
  return (HALF*(Node[ii][jj].X+Node[ii+1][jj].X));
}

/*!
 * Get East face midpoint.
 */
inline Vector2D Grid2D_Quad_Block_HO::xfaceE(const Cell2D_HO &Cell) const {
  return xfaceE(Cell.I,Cell.J);
}

/*!
 * Get East face midpoint.
 *
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 */
inline Vector2D Grid2D_Quad_Block_HO::xfaceE(const int ii, const int jj) const {
  return (HALF*(Node[ii+1][jj].X+Node[ii+1][jj+1].X));
}

/*!
 * Get West face midpoint.
 */
inline Vector2D Grid2D_Quad_Block_HO::xfaceW(const Cell2D_HO &Cell) const {
  return xfaceW(Cell.I,Cell.J);
}

/*!
 * Get West face midpoint.
 *
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 */
inline Vector2D Grid2D_Quad_Block_HO::xfaceW(const int ii, const int jj) const {
  return (HALF*(Node[ii][jj].X+Node[ii][jj+1].X));
}

/*!
 * Get North face length.
 */
inline double Grid2D_Quad_Block_HO::lfaceN(const Cell2D_HO &Cell) const {
  return lfaceN(Cell.I,Cell.J);
}

/*!
 * Get North face length.
 *
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 */
inline double Grid2D_Quad_Block_HO::lfaceN(const int ii, const int jj) const {
  return (abs(Node[ii][jj+1].X-Node[ii+1][jj+1].X));
}

/*!
 * Get South face length.
 */
inline double Grid2D_Quad_Block_HO::lfaceS(const Cell2D_HO &Cell) const {
  return lfaceS(Cell.I,Cell.J);
}

/*!
 * Get South face length.
 *
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 */
inline double Grid2D_Quad_Block_HO::lfaceS(const int ii, const int jj) const {
  return (abs(Node[ii][jj].X-Node[ii+1][jj].X));
}

/*!
 * Get East face length.
 */
inline double Grid2D_Quad_Block_HO::lfaceE(const Cell2D_HO &Cell) const {
  return lfaceE(Cell.I,Cell.J);
}

/*!
 * Get East face length.
 *
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 */
inline double Grid2D_Quad_Block_HO::lfaceE(const int ii, const int jj) const {
  return (abs(Node[ii+1][jj].X-Node[ii+1][jj+1].X));
}

/*!
 * Get West face length.
 */
inline double Grid2D_Quad_Block_HO::lfaceW(const Cell2D_HO &Cell) const {
  return lfaceW(Cell.I,Cell.J);
}

/*!
 * Get West face length.
 *
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 */
inline double Grid2D_Quad_Block_HO::lfaceW(const int ii, const int jj) const {
  return (abs(Node[ii][jj].X-Node[ii][jj+1].X));
}


/*!
 * Get cell perimeter.
 */
inline double Grid2D_Quad_Block_HO::cell_perimeter(const int ii, const int jj) const {
  return lfaceE(ii,jj)+lfaceN(ii,jj)+lfaceW(ii,jj)+lfaceS(ii,jj);
}

/*!
 * Get cell North face normal.
 */
inline Vector2D Grid2D_Quad_Block_HO::nfaceN(const Cell2D_HO &Cell) const {
  return nfaceN(Cell.I,Cell.J);
}

/*!
 * Get cell North face normal.
 *
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 */
inline Vector2D Grid2D_Quad_Block_HO::nfaceN(const int ii, const int jj) const {
  if(lfaceN(ii,jj) > NANO*cell_perimeter(ii,jj)) {
    return (Vector2D( (Node[ii][jj+1].X.y - Node[ii+1][jj+1].X.y),
		      -(Node[ii][jj+1].X.x - Node[ii+1][jj+1].X.x))/
	    abs(Node[ii][jj+1].X - Node[ii+1][jj+1].X));
  } else {
    return ihat;
  }
}

/*!
 * Get cell South face normal.
 */
inline Vector2D Grid2D_Quad_Block_HO::nfaceS(const Cell2D_HO &Cell) const {
  return nfaceS(Cell.I,Cell.J);
}

/*!
 * Get cell South face normal.
 *
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 */
inline Vector2D Grid2D_Quad_Block_HO::nfaceS(const int ii, const int jj) const {
  if(lfaceS(ii,jj) > NANO*cell_perimeter(ii,jj)) {
    return (Vector2D( (Node[ii+1][jj].X.y - Node[ii][jj].X.y),
		      -(Node[ii+1][jj].X.x - Node[ii][jj].X.x))/
	    abs(Node[ii+1][jj].X - Node[ii][jj].X));
  } else {
    return ihat;
  }
}

/*!
 * Get cell East face normal.
 */
inline Vector2D Grid2D_Quad_Block_HO::nfaceE(const Cell2D_HO &Cell) const {
  return nfaceE(Cell.I,Cell.J);
}

/*!
 * Get cell East face normal.
 *
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 */
inline Vector2D Grid2D_Quad_Block_HO::nfaceE(const int ii, const int jj) const {
  if(lfaceE(ii,jj) > NANO*cell_perimeter(ii,jj)) {
    return (Vector2D( (Node[ii+1][jj+1].X.y - Node[ii+1][jj].X.y),
		      -(Node[ii+1][jj+1].X.x - Node[ii+1][jj].X.x))/
	    abs(Node[ii+1][jj+1].X-Node[ii+1][jj].X));
  } else {
    return ihat;
  }
}

/*!
 * Get cell West face normal.
 */
inline Vector2D Grid2D_Quad_Block_HO::nfaceW(const Cell2D_HO &Cell) const {
  return nfaceW(Cell.I,Cell.J);
}

/*!
 * Get cell West face normal.
 *
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 */
inline Vector2D Grid2D_Quad_Block_HO::nfaceW(const int ii, const int jj) const {
  if(lfaceW(ii,jj) > NANO*cell_perimeter(ii,jj)) {
    return (Vector2D( (Node[ii][jj].X.y - Node[ii][jj+1].X.y),
		      -(Node[ii][jj].X.x - Node[ii][jj+1].X.x))/
	    abs(Node[ii][jj].X - Node[ii][jj+1].X));
  } else {
    return ihat;
  }
}

/*!
 * Get the number of Gauss quadrature points for the North face.
 */
inline void Grid2D_Quad_Block_HO::getGaussQuadPointsFaceN(const Cell2D_HO &Cell,
							  Vector2D * GQPoints,
							  const int & NumberOfGQPs) const{
  return getGaussQuadPointsFaceN(Cell.I, Cell.J, GQPoints, NumberOfGQPs);
}

/*!
 * Get the number of Gauss quadrature points for the South face.
 */
inline void Grid2D_Quad_Block_HO::getGaussQuadPointsFaceS(const Cell2D_HO &Cell,
							  Vector2D * GQPoints,
							  const int & NumberOfGQPs) const{
  return getGaussQuadPointsFaceS(Cell.I, Cell.J, GQPoints, NumberOfGQPs);
}

/*!
 * Get the number of Gauss quadrature points for the East face.
 */
inline void Grid2D_Quad_Block_HO::getGaussQuadPointsFaceE(const Cell2D_HO &Cell,
							  Vector2D * GQPoints,
							  const int & NumberOfGQPs) const {
  return getGaussQuadPointsFaceE(Cell.I, Cell.J, GQPoints, NumberOfGQPs);
}

/*!
 * Get the number of Gauss quadrature points for the West face.
 */
inline void Grid2D_Quad_Block_HO::getGaussQuadPointsFaceW(const Cell2D_HO &Cell, 
							  Vector2D * GQPoints,
							  const int & NumberOfGQPs) const {
  return getGaussQuadPointsFaceW(Cell.I, Cell.J, GQPoints, NumberOfGQPs);
}

/*!
 * Get the number of Gauss quadrature points for the North face.
 * 
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 * \param GQPoints storage array for the Gauss quadrature points. This memory is overwritten!
 * \param NumberOfGQPs specifies how many points are returned. This number is typically dictated 
 *                     by the accuracy of the flux calculation.
 */
inline void Grid2D_Quad_Block_HO::getGaussQuadPointsFaceN(const int &ii, const int &jj,
							  Vector2D * GQPoints, const int & NumberOfGQPs) const{
  switch (NumberOfGQPs){
  case 1:
    GQPoints[0] = xfaceN(ii,jj);
    break;

  case 2:
    GQPoints[0] = GQPoints[1] = Node[ii][jj+1].X-Node[ii+1][jj+1].X;
    
    /* final value GQPoints[0] */
    GQPoints[0] = Node[ii+1][jj+1].X + GaussQuadratureData::GQ2_Abscissa[0]*GQPoints[0];
    
    /* final value GQPoints[1] */
    GQPoints[1] = Node[ii+1][jj+1].X + GaussQuadratureData::GQ2_Abscissa[1]*GQPoints[1];
    break;

  case 3:
    GQPoints[0] = GQPoints[1] = GQPoints[2] = Node[ii][jj+1].X-Node[ii+1][jj+1].X;
    
    /* final value GQPoints[0] */
    GQPoints[0] = Node[ii+1][jj+1].X + GaussQuadratureData::GQ3_Abscissa[0]*GQPoints[0];
   
    /* final value GQPoints[1] */
    GQPoints[1] = Node[ii+1][jj+1].X + GaussQuadratureData::GQ3_Abscissa[1]*GQPoints[1];

    /* final value GQPoints[2] */
    GQPoints[2] = Node[ii+1][jj+1].X + GaussQuadratureData::GQ3_Abscissa[2]*GQPoints[2];
    break;

  case 5:
    GQPoints[0] = GQPoints[1] = GQPoints[2] = GQPoints[3] = GQPoints[4] = Node[ii][jj+1].X-Node[ii+1][jj+1].X;
    
    /* final value GQPoints[0] */
    GQPoints[0] = Node[ii+1][jj+1].X + GaussQuadratureData::GQ5_Abscissa[0]*GQPoints[0];
   
    /* final value GQPoints[1] */
    GQPoints[1] = Node[ii+1][jj+1].X + GaussQuadratureData::GQ5_Abscissa[1]*GQPoints[1];

    /* final value GQPoints[2] */
    GQPoints[2] = Node[ii+1][jj+1].X + GaussQuadratureData::GQ5_Abscissa[2]*GQPoints[2];

    /* final value GQPoints[3] */
    GQPoints[3] = Node[ii+1][jj+1].X + GaussQuadratureData::GQ5_Abscissa[3]*GQPoints[3];

    /* final value GQPoints[4] */
    GQPoints[4] = Node[ii+1][jj+1].X + GaussQuadratureData::GQ5_Abscissa[4]*GQPoints[4];
    break;

  default:
    throw runtime_error("Grid2D_Quad_Block_HO::getGaussQuadPointsFaceN() ERROR! \
                         Not implemented number of Gauss quadrature points!");
  }
}

/*!
 * Get the number of Gauss quadrature points for the North face.
 * 
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 * \param GQPoints storage array for the Gauss quadrature points. The memory IS NOT overwritten!
 * \param NumberOfGQPs specifies how many points are returned. This number is typically dictated 
 *                     by the accuracy of the flux calculation.
 */
inline void Grid2D_Quad_Block_HO::addGaussQuadPointsFaceN(const int &ii, const int &jj,
							  std::vector<Vector2D> &GQPoints, const int & NumberOfGQPs) const{

  Vector2D Temp;

  switch (NumberOfGQPs){
  case 1:
    GQPoints.push_back(xfaceN(ii,jj));
    break;

  case 2:
    Temp = Node[ii][jj+1].X-Node[ii+1][jj+1].X;
    
    /* final value GQPoints[0] */
    GQPoints.push_back(Node[ii+1][jj+1].X + GaussQuadratureData::GQ2_Abscissa[0]*Temp);
    
    /* final value GQPoints[1] */
    GQPoints.push_back(Node[ii+1][jj+1].X + GaussQuadratureData::GQ2_Abscissa[1]*Temp);
    break;

  case 3:
    Temp = Node[ii][jj+1].X-Node[ii+1][jj+1].X;
    
    /* final value GQPoints[0] */
    GQPoints.push_back(Node[ii+1][jj+1].X + GaussQuadratureData::GQ3_Abscissa[0]*Temp);
   
    /* final value GQPoints[1] */
    GQPoints.push_back(Node[ii+1][jj+1].X + GaussQuadratureData::GQ3_Abscissa[1]*Temp);

    /* final value GQPoints[2] */
    GQPoints.push_back(Node[ii+1][jj+1].X + GaussQuadratureData::GQ3_Abscissa[2]*Temp);
    break;

  case 5:
    Temp = Node[ii][jj+1].X-Node[ii+1][jj+1].X;
    
    /* final value GQPoints[0] */
    GQPoints.push_back(Node[ii+1][jj+1].X + GaussQuadratureData::GQ5_Abscissa[0]*Temp);
   
    /* final value GQPoints[1] */
    GQPoints.push_back(Node[ii+1][jj+1].X + GaussQuadratureData::GQ5_Abscissa[1]*Temp);

    /* final value GQPoints[2] */
    GQPoints.push_back(Node[ii+1][jj+1].X + GaussQuadratureData::GQ5_Abscissa[2]*Temp);

    /* final value GQPoints[3] */
    GQPoints.push_back(Node[ii+1][jj+1].X + GaussQuadratureData::GQ5_Abscissa[3]*Temp);

    /* final value GQPoints[4] */
    GQPoints.push_back(Node[ii+1][jj+1].X + GaussQuadratureData::GQ5_Abscissa[4]*Temp);
    break;

  default:
    throw runtime_error("Grid2D_Quad_Block_HO::getGaussQuadPointsFaceN() ERROR! \
                         Not implemented number of Gauss quadrature points!");
  }
}

/*!
 * Get the number of Gauss quadrature points for the South face.
 * 
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 * \param GQPoints storage array for the Gauss quadrature points. This memory is overwritten!
 * \param NumberOfGQPs specifies how many points are returned. This number is typically dictated 
 *                     by the accuracy of the flux calculation.
 */							  
inline void Grid2D_Quad_Block_HO::getGaussQuadPointsFaceS(const int &ii, const int &jj,
							  Vector2D * GQPoints, const int & NumberOfGQPs) const {

  switch (NumberOfGQPs){
  case 1:
    GQPoints[0] = xfaceS(ii,jj);
    break;

  case 2:
    GQPoints[0] = GQPoints[1] = Node[ii+1][jj].X-Node[ii][jj].X;
    
    /* final value GQPoints[0] */
    GQPoints[0] = Node[ii][jj].X + GaussQuadratureData::GQ2_Abscissa[0]*GQPoints[0];
    
    /* final value GQPoints[1] */
    GQPoints[1] = Node[ii][jj].X + GaussQuadratureData::GQ2_Abscissa[1]*GQPoints[1];
    break;

  case 3:
    GQPoints[0] = GQPoints[1] = GQPoints[2] = Node[ii+1][jj].X-Node[ii][jj].X;
    
    /* final value GQPoints[0] */
    GQPoints[0] = Node[ii][jj].X + GaussQuadratureData::GQ3_Abscissa[0]*GQPoints[0];
    
    /* final value GQPoints[1] */
    GQPoints[1] = Node[ii][jj].X + GaussQuadratureData::GQ3_Abscissa[1]*GQPoints[1];

    /* final value GQPoints[2] */
    GQPoints[2] = Node[ii][jj].X + GaussQuadratureData::GQ3_Abscissa[2]*GQPoints[2];
    break;

  case 5:
    GQPoints[0] = GQPoints[1] = GQPoints[2] = GQPoints[3] = GQPoints[4] = Node[ii+1][jj].X-Node[ii][jj].X;
    
    /* final value GQPoints[0] */
    GQPoints[0] = Node[ii][jj].X + GaussQuadratureData::GQ5_Abscissa[0]*GQPoints[0];
    
    /* final value GQPoints[1] */
    GQPoints[1] = Node[ii][jj].X + GaussQuadratureData::GQ5_Abscissa[1]*GQPoints[1];

    /* final value GQPoints[2] */
    GQPoints[2] = Node[ii][jj].X + GaussQuadratureData::GQ5_Abscissa[2]*GQPoints[2];

    /* final value GQPoints[3] */
    GQPoints[3] = Node[ii][jj].X + GaussQuadratureData::GQ5_Abscissa[3]*GQPoints[3];

    /* final value GQPoints[4] */
    GQPoints[4] = Node[ii][jj].X + GaussQuadratureData::GQ5_Abscissa[4]*GQPoints[4];
    break;

  default:
    throw runtime_error("Grid2D_Quad_Block_HO::getGaussQuadPointsFaceS() ERROR! \
                         Not implemented number of Gauss quadrature points!");
  }
}

/*!
 * Get the number of Gauss quadrature points for the South face.
 * 
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 * \param GQPoints storage array for the Gauss quadrature points. This memory IS NOT overwritten!
 * \param NumberOfGQPs specifies how many points are returned. This number is typically dictated 
 *                     by the accuracy of the flux calculation.
 */							  
inline void Grid2D_Quad_Block_HO::addGaussQuadPointsFaceS(const int &ii, const int &jj,
							  std::vector<Vector2D> &GQPoints, const int & NumberOfGQPs) const {

  Vector2D Temp;

  switch (NumberOfGQPs){
  case 1:
    GQPoints.push_back(xfaceS(ii,jj));
    break;

  case 2:
    Temp = Node[ii+1][jj].X-Node[ii][jj].X;
    
    /* final value GQPoints[0] */
    GQPoints.push_back(Node[ii][jj].X + GaussQuadratureData::GQ2_Abscissa[0]*Temp);
    
    /* final value GQPoints[1] */
    GQPoints.push_back(Node[ii][jj].X + GaussQuadratureData::GQ2_Abscissa[1]*Temp);
    break;

  case 3:
    Temp = Node[ii+1][jj].X-Node[ii][jj].X;
    
    /* final value GQPoints[0] */
    GQPoints.push_back(Node[ii][jj].X + GaussQuadratureData::GQ3_Abscissa[0]*Temp);
    
    /* final value GQPoints[1] */
    GQPoints.push_back(Node[ii][jj].X + GaussQuadratureData::GQ3_Abscissa[1]*Temp);

    /* final value GQPoints[2] */
    GQPoints.push_back(Node[ii][jj].X + GaussQuadratureData::GQ3_Abscissa[2]*Temp);
    break;

  case 5:
    Temp = Node[ii+1][jj].X-Node[ii][jj].X;
    
    /* final value GQPoints[0] */
    GQPoints.push_back(Node[ii][jj].X + GaussQuadratureData::GQ5_Abscissa[0]*Temp);
    
    /* final value GQPoints[1] */
    GQPoints.push_back(Node[ii][jj].X + GaussQuadratureData::GQ5_Abscissa[1]*Temp);

    /* final value GQPoints[2] */
    GQPoints.push_back(Node[ii][jj].X + GaussQuadratureData::GQ5_Abscissa[2]*Temp);

    /* final value GQPoints[3] */
    GQPoints.push_back(Node[ii][jj].X + GaussQuadratureData::GQ5_Abscissa[3]*Temp);

    /* final value GQPoints[4] */
    GQPoints.push_back(Node[ii][jj].X + GaussQuadratureData::GQ5_Abscissa[4]*Temp);
    break;

  default:
    throw runtime_error("Grid2D_Quad_Block_HO::getGaussQuadPointsFaceS() ERROR! \
                         Not implemented number of Gauss quadrature points!");
  }
}

/*!
 * Get the number of Gauss quadrature points for the East face.
 * 
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 * \param GQPoints storage array for the Gauss quadrature points. This memory is overwritten!
 * \param NumberOfGQPs specifies how many points are returned. This number is typically dictated 
 *                     by the accuracy of the flux calculation.
 */
inline void Grid2D_Quad_Block_HO::getGaussQuadPointsFaceE(const int &ii, const int &jj,
							  Vector2D * GQPoints, const int & NumberOfGQPs) const {

  switch (NumberOfGQPs){
  case 1:
    GQPoints[0] = xfaceE(ii,jj);
    break;

  case 2:
    GQPoints[0] = GQPoints[1] = Node[ii+1][jj+1].X-Node[ii+1][jj].X;
    
    /* final value GQPoints[0] */
    GQPoints[0] = Node[ii+1][jj].X + GaussQuadratureData::GQ2_Abscissa[0]*GQPoints[0];
    
    /* final value GQPoints[1] */
    GQPoints[1] = Node[ii+1][jj].X + GaussQuadratureData::GQ2_Abscissa[1]*GQPoints[1];
    break;

  case 3:
    GQPoints[0] = GQPoints[1] = GQPoints[2] = Node[ii+1][jj+1].X-Node[ii+1][jj].X;
    
    /* final value GQPoints[0] */
    GQPoints[0] = Node[ii+1][jj].X + GaussQuadratureData::GQ3_Abscissa[0]*GQPoints[0];
    
    /* final value GQPoints[1] */
    GQPoints[1] = Node[ii+1][jj].X + GaussQuadratureData::GQ3_Abscissa[1]*GQPoints[1];

    /* final value GQPoints[2] */
    GQPoints[2] = Node[ii+1][jj].X + GaussQuadratureData::GQ3_Abscissa[2]*GQPoints[2];
    break;

  case 5:
    GQPoints[0] = GQPoints[1] = GQPoints[2] = GQPoints[3] = GQPoints[4] = Node[ii+1][jj+1].X-Node[ii+1][jj].X;
    
    /* final value GQPoints[0] */
    GQPoints[0] = Node[ii+1][jj].X + GaussQuadratureData::GQ5_Abscissa[0]*GQPoints[0];
    
    /* final value GQPoints[1] */
    GQPoints[1] = Node[ii+1][jj].X + GaussQuadratureData::GQ5_Abscissa[1]*GQPoints[1];

    /* final value GQPoints[2] */
    GQPoints[2] = Node[ii+1][jj].X + GaussQuadratureData::GQ5_Abscissa[2]*GQPoints[2];

    /* final value GQPoints[3] */
    GQPoints[3] = Node[ii+1][jj].X + GaussQuadratureData::GQ5_Abscissa[3]*GQPoints[3];

    /* final value GQPoints[4] */
    GQPoints[4] = Node[ii+1][jj].X + GaussQuadratureData::GQ5_Abscissa[4]*GQPoints[4];
    break;

  default:
    throw runtime_error("Grid2D_Quad_Block_HO::getGaussQuadPointsFaceE() ERROR! \
                         Not implemented number of Gauss quadrature points!");
  }
}

/*!
 * Get the number of Gauss quadrature points for the East face.
 * 
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 * \param GQPoints storage array for the Gauss quadrature points. This memory IS NOT overwritten!
 * \param NumberOfGQPs specifies how many points are returned. This number is typically dictated 
 *                     by the accuracy of the flux calculation.
 */
inline void Grid2D_Quad_Block_HO::addGaussQuadPointsFaceE(const int &ii, const int &jj,
							  std::vector<Vector2D> & GQPoints, const int & NumberOfGQPs) const {
  
  Vector2D Temp;
  
  switch (NumberOfGQPs){
  case 1:
    GQPoints.push_back(xfaceE(ii,jj));
    break;

  case 2:
    Temp = Node[ii+1][jj+1].X-Node[ii+1][jj].X;
    
    /* final value GQPoints[0] */
    GQPoints.push_back(Node[ii+1][jj].X + GaussQuadratureData::GQ2_Abscissa[0]*Temp);
    
    /* final value GQPoints[1] */
    GQPoints.push_back(Node[ii+1][jj].X + GaussQuadratureData::GQ2_Abscissa[1]*Temp);
    break;

  case 3:
    Temp = Node[ii+1][jj+1].X-Node[ii+1][jj].X;
    
    /* final value GQPoints[0] */
    GQPoints.push_back(Node[ii+1][jj].X + GaussQuadratureData::GQ3_Abscissa[0]*Temp);
    
    /* final value GQPoints[1] */
    GQPoints.push_back(Node[ii+1][jj].X + GaussQuadratureData::GQ3_Abscissa[1]*Temp);

    /* final value GQPoints[2] */
    GQPoints.push_back(Node[ii+1][jj].X + GaussQuadratureData::GQ3_Abscissa[2]*Temp);
    break;

  case 5:
    Temp = Node[ii+1][jj+1].X-Node[ii+1][jj].X;
    
    /* final value GQPoints[0] */
    GQPoints.push_back(Node[ii+1][jj].X + GaussQuadratureData::GQ5_Abscissa[0]*Temp);
    
    /* final value GQPoints[1] */
    GQPoints.push_back(Node[ii+1][jj].X + GaussQuadratureData::GQ5_Abscissa[1]*Temp);

    /* final value GQPoints[2] */
    GQPoints.push_back(Node[ii+1][jj].X + GaussQuadratureData::GQ5_Abscissa[2]*Temp);

    /* final value GQPoints[3] */
    GQPoints.push_back(Node[ii+1][jj].X + GaussQuadratureData::GQ5_Abscissa[3]*Temp);

    /* final value GQPoints[4] */
    GQPoints.push_back(Node[ii+1][jj].X + GaussQuadratureData::GQ5_Abscissa[4]*Temp);
    break;

  default:
    throw runtime_error("Grid2D_Quad_Block_HO::getGaussQuadPointsFaceE() ERROR! \
                         Not implemented number of Gauss quadrature points!");
  }
}

/*!
 * Get the number of Gauss quadrature points for the West face.
 * 
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 * \param GQPoints storage array for the Gauss quadrature points. This memory is overwritten!
 * \param NumberOfGQPs specifies how many points are returned. This number is typically dictated 
 *                     by the accuracy of the flux calculation.
 */
inline void Grid2D_Quad_Block_HO::getGaussQuadPointsFaceW(const int &ii, const int &jj,
							  Vector2D * GQPoints, const int & NumberOfGQPs) const {

  switch (NumberOfGQPs){
  case 1:
    GQPoints[0] = xfaceW(ii,jj);
    break;

  case 2:
    GQPoints[0] = GQPoints[1] = Node[ii][jj].X-Node[ii][jj+1].X;
    
    /* final value GQPoints[0] */
    GQPoints[0] = Node[ii][jj+1].X + GaussQuadratureData::GQ2_Abscissa[0]*GQPoints[0];
    
    /* final value GQPoints[1] */
    GQPoints[1] = Node[ii][jj+1].X + GaussQuadratureData::GQ2_Abscissa[1]*GQPoints[1];
    break;

  case 3:
    GQPoints[0] = GQPoints[1] = GQPoints[2] = Node[ii][jj].X-Node[ii][jj+1].X;
    
    /* final value GQPoints[0] */
    GQPoints[0] = Node[ii][jj+1].X + GaussQuadratureData::GQ3_Abscissa[0]*GQPoints[0];
    
    /* final value GQPoints[1] */
    GQPoints[1] = Node[ii][jj+1].X + GaussQuadratureData::GQ3_Abscissa[1]*GQPoints[1];

    /* final value GQPoints[2] */
    GQPoints[2] = Node[ii][jj+1].X + GaussQuadratureData::GQ3_Abscissa[2]*GQPoints[2];
    break;

  case 5:
    GQPoints[0] = GQPoints[1] = GQPoints[2] = GQPoints[3] = GQPoints[4] = Node[ii][jj].X-Node[ii][jj+1].X;
    
    /* final value GQPoints[0] */
    GQPoints[0] = Node[ii][jj+1].X + GaussQuadratureData::GQ5_Abscissa[0]*GQPoints[0];
    
    /* final value GQPoints[1] */
    GQPoints[1] = Node[ii][jj+1].X + GaussQuadratureData::GQ5_Abscissa[1]*GQPoints[1];

    /* final value GQPoints[2] */
    GQPoints[2] = Node[ii][jj+1].X + GaussQuadratureData::GQ5_Abscissa[2]*GQPoints[2];

    /* final value GQPoints[3] */
    GQPoints[3] = Node[ii][jj+1].X + GaussQuadratureData::GQ5_Abscissa[3]*GQPoints[3];

    /* final value GQPoints[4] */
    GQPoints[4] = Node[ii][jj+1].X + GaussQuadratureData::GQ5_Abscissa[4]*GQPoints[4];
    break;

  default:
    throw runtime_error("Grid2D_Quad_Block_HO::getGaussQuadPointsFaceW() ERROR! \
                         Not implemented number of Gauss quadrature points!");
  }
}

/*!
 * Get the number of Gauss quadrature points for the West face.
 * 
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 * \param GQPoints storage array for the Gauss quadrature points. This memory IS NOT overwritten!
 * \param NumberOfGQPs specifies how many points are returned. This number is typically dictated 
 *                     by the accuracy of the flux calculation.
 */
inline void Grid2D_Quad_Block_HO::addGaussQuadPointsFaceW(const int &ii, const int &jj,
							  std::vector<Vector2D> & GQPoints, const int & NumberOfGQPs) const {

  Vector2D Temp;

  switch (NumberOfGQPs){
  case 1:
    GQPoints.push_back(xfaceW(ii,jj));
    break;

  case 2:
    Temp = Node[ii][jj].X-Node[ii][jj+1].X;
    
    /* final value GQPoints[0] */
    GQPoints.push_back(Node[ii][jj+1].X + GaussQuadratureData::GQ2_Abscissa[0]*Temp);
    									  
    /* final value GQPoints[1] */					  
    GQPoints.push_back(Node[ii][jj+1].X + GaussQuadratureData::GQ2_Abscissa[1]*Temp);
    break;

  case 3:
    Temp = Node[ii][jj].X-Node[ii][jj+1].X;
    
    /* final value GQPoints[0] */
    GQPoints.push_back(Node[ii][jj+1].X + GaussQuadratureData::GQ3_Abscissa[0]*Temp);
    									  
    /* final value GQPoints[1] */					  
    GQPoints.push_back(Node[ii][jj+1].X + GaussQuadratureData::GQ3_Abscissa[1]*Temp);
									  
    /* final value GQPoints[2] */					  
    GQPoints.push_back(Node[ii][jj+1].X + GaussQuadratureData::GQ3_Abscissa[2]*Temp);
    break;

  case 5:
    Temp = Node[ii][jj].X-Node[ii][jj+1].X;
    
    /* final value GQPoints[0] */
    GQPoints.push_back(Node[ii][jj+1].X + GaussQuadratureData::GQ5_Abscissa[0]*Temp);
    									  
    /* final value GQPoints[1] */					  
    GQPoints.push_back(Node[ii][jj+1].X + GaussQuadratureData::GQ5_Abscissa[1]*Temp);
									  
    /* final value GQPoints[2] */					  
    GQPoints.push_back(Node[ii][jj+1].X + GaussQuadratureData::GQ5_Abscissa[2]*Temp);
									  
    /* final value GQPoints[3] */					  
    GQPoints.push_back(Node[ii][jj+1].X + GaussQuadratureData::GQ5_Abscissa[3]*Temp);
									  
    /* final value GQPoints[4] */					  
    GQPoints.push_back(Node[ii][jj+1].X + GaussQuadratureData::GQ5_Abscissa[4]*Temp);
    break;

  default:
    throw runtime_error("Grid2D_Quad_Block_HO::getGaussQuadPointsFaceW() ERROR! \
                         Not implemented number of Gauss quadrature points!");
  }
}

/*!
 * Update cell information assuming straight boundaries
 * (i.e. every edge of the cell is a line segment)
 *
 * \param iCell i-index of the cell
 * \param jCell j-index of the cell
 */
inline void Grid2D_Quad_Block_HO::Update_Cell(const int & iCell, const int & jCell){

  static Vector2D Shift, SW, SE, NE, NW;

  // Set cell indexes
  Cell[iCell][jCell].I = iCell;
  Cell[iCell][jCell].J = jCell;

  // Translate the cell nodes into the local coordinate system if required
  TranslateVertexesIntoLocalCoordinateSystem(Node[iCell  ][jCell  ].X, Node[iCell+1][jCell  ].X,
					     Node[iCell+1][jCell+1].X, Node[iCell  ][jCell+1].X);

  // Compute cell area and centroid
  quadAreaAndCentroid(iCell,jCell);

  // Compute geometric moments 
  ComputeGeometricCoefficients(iCell,jCell);

  // Translate the cell nodes back to the global coordinate system.
  TranslateVertexesIntoGlobalCoordinateSystem(Node[iCell  ][jCell  ].X, Node[iCell+1][jCell  ].X,
					      Node[iCell+1][jCell+1].X, Node[iCell  ][jCell+1].X,
					      Cell[iCell][jCell].Xc);
}

/*!
 * This function gets used during message passing.
 * It doesn't actually update the geometric information
 * of the cell with indexes (iCell,jCell).
 * Instead it schedules for update all ghost cells.
 * That's because of the complexity to update ghost 
 * cell near curved boundaries.
 *
 * \param iCell dummy i-index of the cell
 * \param jCell dummy j-index of the cell
 */
inline void Grid2D_Quad_Block_HO::Update_GhostCellProperties_DuringMessagePassing(const int &iCell,
										  const int &jCell){
  Schedule_Ghost_Cells_Update();
}

/*!
 * Positive shift operator.
 */
inline void Grid2D_Quad_Block_HO::operator +(const Vector2D &V) {

  int i, j;
  for ( j = JNl-Nghost ; j <= JNu+Nghost; ++j ) {
    for ( i = INl-Nghost ; i <= INu+Nghost; ++i ) {
      Node[i][j].X += V;
    } /* endfor */
  } /* endfor */

  if (BndNorthSpline.np != 0 ) BndNorthSpline.Translate_Spline(V);
  if (BndSouthSpline.np != 0 ) BndSouthSpline.Translate_Spline(V);
  if (BndEastSpline.np != 0 ) BndEastSpline.Translate_Spline(V);
  if (BndWestSpline.np != 0 ) BndWestSpline.Translate_Spline(V);

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
 * Negative shift operator.
 */
inline void Grid2D_Quad_Block_HO::operator -(const Vector2D &V) {
  return (*this + (-V));
}

/*!
 * Scaling operators.
 */
inline void Grid2D_Quad_Block_HO::operator *(const double &a) {

  int i, j;
  for ( j = JNl-Nghost ; j <= JNu+Nghost; ++j ) {
    for ( i = INl-Nghost ; i <= INu+Nghost; ++i ) {
      Node[i][j].X = Node[i][j].X*a;
    } /* endfor */
  } /* endfor */

  if (BndNorthSpline.np != 0 ) BndNorthSpline.Scale_Spline(a);
  if (BndSouthSpline.np != 0 ) BndSouthSpline.Scale_Spline(a);
  if (BndEastSpline.np != 0 ) BndEastSpline.Scale_Spline(a);
  if (BndWestSpline.np != 0 ) BndWestSpline.Scale_Spline(a);

  if (ExtendWest_BndNorthSpline.np != 0) ExtendWest_BndNorthSpline.Scale_Spline(a);
  if (ExtendEast_BndNorthSpline.np != 0) ExtendEast_BndNorthSpline.Scale_Spline(a);
  if (ExtendWest_BndSouthSpline.np != 0) ExtendWest_BndSouthSpline.Scale_Spline(a);
  if (ExtendEast_BndSouthSpline.np != 0) ExtendEast_BndSouthSpline.Scale_Spline(a);
  if (ExtendNorth_BndEastSpline.np != 0) ExtendNorth_BndEastSpline.Scale_Spline(a);
  if (ExtendSouth_BndEastSpline.np != 0) ExtendSouth_BndEastSpline.Scale_Spline(a);
  if (ExtendNorth_BndWestSpline.np != 0) ExtendNorth_BndWestSpline.Scale_Spline(a);
  if (ExtendSouth_BndWestSpline.np != 0) ExtendSouth_BndWestSpline.Scale_Spline(a);


  SminN *= a;
  SmaxN *= a;
  SminS *= a;
  SmaxS *= a;
  SminE *= a;
  SmaxE *= a;
  SminW *= a;
  SmaxW *= a;

  /* Require update of the whole mesh */
  Schedule_Interior_Mesh_Update();
  Schedule_Ghost_Cells_Update();
}

/*!
 * Rotation operators.
 */
inline void Grid2D_Quad_Block_HO::operator ^(const double &a) {

  int i, j;
  double cos_angle, sin_angle;
  Vector2D X;

  cos_angle = cos(-a);
  sin_angle = sin(-a);

  for ( j = JNl-Nghost ; j <= JNu+Nghost; ++j ) {
    for ( i = INl-Nghost ; i <= INu+Nghost; ++i ) {
      X.x = ( Node[i][j].X.x*cos_angle +
	      Node[i][j].X.y*sin_angle );
      X.y = (- Node[i][j].X.x*sin_angle +
	     Node[i][j].X.y*cos_angle );
      Node[i][j].X = X;
    } /* endfor */
  } /* endfor */

  if (BndNorthSpline.np != 0 ) BndNorthSpline.Rotate_Spline(a);
  if (BndSouthSpline.np != 0 ) BndSouthSpline.Rotate_Spline(a);
  if (BndEastSpline.np != 0 ) BndEastSpline.Rotate_Spline(a);
  if (BndWestSpline.np != 0 ) BndWestSpline.Rotate_Spline(a);

  if (ExtendWest_BndNorthSpline.np != 0) ExtendWest_BndNorthSpline.Rotate_Spline(a);
  if (ExtendEast_BndNorthSpline.np != 0) ExtendEast_BndNorthSpline.Rotate_Spline(a);
  if (ExtendWest_BndSouthSpline.np != 0) ExtendWest_BndSouthSpline.Rotate_Spline(a);
  if (ExtendEast_BndSouthSpline.np != 0) ExtendEast_BndSouthSpline.Rotate_Spline(a);
  if (ExtendNorth_BndEastSpline.np != 0) ExtendNorth_BndEastSpline.Rotate_Spline(a);
  if (ExtendSouth_BndEastSpline.np != 0) ExtendSouth_BndEastSpline.Rotate_Spline(a);
  if (ExtendNorth_BndWestSpline.np != 0) ExtendNorth_BndWestSpline.Rotate_Spline(a);
  if (ExtendSouth_BndWestSpline.np != 0) ExtendSouth_BndWestSpline.Rotate_Spline(a);

  /* Require update of the whole mesh */
  Schedule_Interior_Mesh_Update();
  Schedule_Ghost_Cells_Update();
}

/*!
 * Return the number of Gauss quadrature points at which 
 * the boundary conditions are enforced as constraints
 * for a CellIndex cell, a given boundary spline and 
 * the associate spline interval.
 *
 * \param CellIndex_Shift the difference between the CellIndex of the cell and
 *                        the corresponding position in BndSplineInfo variable.
 */
inline int Grid2D_Quad_Block_HO::NumOfConstrainedGaussQuadPoints_GenericBoundary(const BndSplineType & BndSpline,
										 const BndSplineIntervalType * BndSplineInfo,
										 const int &CellIndex,
										 const int &CellIndex_Shift){

  if (BndSpline.getFluxCalcMethod() == SolveRiemannProblem){
    /* The boundary flux is computed by solving a Riemann problem at the interface with the ghost cell */
    return 0;
  } else if (BndSplineInfo != NULL){
    /* use high-order boundary info */
    return BndSplineInfo[CellIndex-CellIndex_Shift].NumGQPoints();
  } else {
    /* return the number of points based on the order of accuracy
       (This situation corresponds to "Low-order boundaries + ReconstructionBasedFlux" ) */
    return NumGQP;
  }
}

/*!
 * Return the number of Gauss quadrature points on the North
 * cell face which have boundary conditions enforced by constraints.
 */
inline int Grid2D_Quad_Block_HO::NumOfConstrainedGaussQuadPoints_North(const int &ii, const int &jj){
  if (jj == JCu){
    if (ii < ICl){
      // This cell is bounded by ExtendWest_BndNorthSpline to North
      return NumOfConstrainedGaussQuadPoints_GenericBoundary(ExtendWest_BndNorthSpline,
      							     ExtendWest_BndNorthSplineInfo,
      							     ii);
    } else if (ii <= ICu){
      // This cell is bounded by BndNorthSpline to North
      return NumOfConstrainedGaussQuadPoints_GenericBoundary(BndNorthSpline,
							     BndNorthSplineInfo,
							     ii);      
    } else {
      // This cell is bounded by ExtendEast_BndNorthSpline to North
      return NumOfConstrainedGaussQuadPoints_GenericBoundary(ExtendEast_BndNorthSpline,
							     ExtendEast_BndNorthSplineInfo,
							     ii, ICu+1);
    }
  } else {
    /* This cell is not on the interior side of North block boundaries */
    return 0;
  }
}

/*!
 * Return the number of Gauss quadrature points on the South
 * cell face which have boundary conditions enforced by constraints.
 */
inline int Grid2D_Quad_Block_HO::NumOfConstrainedGaussQuadPoints_South(const int &ii, const int &jj){
  if (jj == JCl){
    if (ii < ICl){
      // This cell is bounded by ExtendWest_BndSouthSpline to South
      return NumOfConstrainedGaussQuadPoints_GenericBoundary(ExtendWest_BndSouthSpline,
							     ExtendWest_BndSouthSplineInfo,
							     ii);
    } else if (ii <= ICu){
      // This cell is bounded by BndSouthSpline to South
      return NumOfConstrainedGaussQuadPoints_GenericBoundary(BndSouthSpline,
							     BndSouthSplineInfo,
							     ii);      
    } else {
      // This cell is bounded by ExtendEast_BndSouthSpline to South
      return NumOfConstrainedGaussQuadPoints_GenericBoundary(ExtendEast_BndSouthSpline,
							     ExtendEast_BndSouthSplineInfo,
							     ii, ICu+1);
    }
  } else {
    /* This cell is not on the interior side of South block boundaries */
    return 0;
  }
}

/*!
 * Return the number of Gauss quadrature points on the East
 * cell face which have boundary conditions enforced by constraints.
 */
inline int Grid2D_Quad_Block_HO::NumOfConstrainedGaussQuadPoints_East(const int &ii, const int &jj){
  if (ii == ICu){
    if (jj < JCl){
      // This cell is bounded by ExtendSouth_BndEastSpline to East
      return NumOfConstrainedGaussQuadPoints_GenericBoundary(ExtendSouth_BndEastSpline,
							     ExtendSouth_BndEastSplineInfo,
							     jj);
    } else if (jj <= JCu){
      // This cell is bounded by BndEastSpline to East
      return NumOfConstrainedGaussQuadPoints_GenericBoundary(BndEastSpline,
							     BndEastSplineInfo,
							     jj);
    } else {
      // This cell is bounded by ExtendNorth_BndEastSpline to East
      return NumOfConstrainedGaussQuadPoints_GenericBoundary(ExtendNorth_BndEastSpline,
							     ExtendNorth_BndEastSplineInfo,
							     jj, JCu+1);
    }
  } else {
    /* This cell is not on the interior side of East block boundaries */
    return 0;
  }
}

/*!
 * Return the number of Gauss quadrature points on the West
 * cell face which have boundary conditions enforced by constraints.
 */
inline int Grid2D_Quad_Block_HO::NumOfConstrainedGaussQuadPoints_West(const int &ii, const int &jj){
  if (ii == ICl){
    if (jj < JCl){
      // This cell is bounded by ExtendSouth_BndWestSpline to West
      return NumOfConstrainedGaussQuadPoints_GenericBoundary(ExtendSouth_BndWestSpline,
							     ExtendSouth_BndWestSplineInfo,
							     jj);
    } else if (jj <= JCu){
      // This cell is bounded by BndWestSpline to West
      return NumOfConstrainedGaussQuadPoints_GenericBoundary(BndWestSpline,
							     BndWestSplineInfo,
							     jj);
    } else {
      // This cell is bounded by ExtendNorth_BndWestSpline to West
      return NumOfConstrainedGaussQuadPoints_GenericBoundary(ExtendNorth_BndWestSpline,
							     ExtendNorth_BndWestSplineInfo,
							     jj, JCu+1);
    }
  } else {
    /* This cell is not on the interior side of West block boundaries */
    return 0;
  }
}

/*!
 * Return the total number of Gauss quadrature points for the specified cell
 * which have boundary conditions enforced by constraints.
 */
inline int Grid2D_Quad_Block_HO::NumOfConstrainedGaussQuadPoints(const int &ii, const int &jj){
  return ( NumOfConstrainedGaussQuadPoints_North(ii,jj) + NumOfConstrainedGaussQuadPoints_South(ii,jj) +
	   NumOfConstrainedGaussQuadPoints_West(ii,jj) + NumOfConstrainedGaussQuadPoints_East(ii,jj));
}

/*!
 * Return the number of Gauss quadrature points on the North
 * cell face at which flux calculation is performed.
 *
 * \param curved_face Set to true if the face is treated as curved otherwise to false.
 */
inline int Grid2D_Quad_Block_HO::NumOfFluxCalculationGaussQuadPoints_North(const int &ii, const int &jj,
									   bool & curved_face){

  if (!IsHighOrderBoundary()){
    // return the number of points based on the order of accuracy
    // if no high-order treatment of geometry is required
    curved_face = false;
    return NumGQP;
  }

  // Handle high-order geometry
  if (ii >= ICl && ii<=ICu){
    // These cells have potentially curved North faces
    if (jj == JCu){		// check last layer of interior cells on the North block boundary
      if (BndNorthSplineInfo != NULL){
	/* use high-order boundary info */
	curved_face = true;
	return BndNorthSplineInfo[ii].NumGQPoints();
      } else {
	curved_face = false;
	return NumGQP;
      }
    } else if (jj == JCl - 1){	// check first layer of ghost cells on the South block boundary
      if (BndSouthSplineInfo != NULL){
	/* use high-order boundary info */
	curved_face = true;
	return BndSouthSplineInfo[ii].NumGQPoints();
      } else {
	curved_face = false;
	return NumGQP;
      }
    } else {
      curved_face = false;
      return NumGQP;
    }
  } else {
    // These cells don't have any potentially curved North face
    // Return the number of points based on the order of accuracy
    curved_face = false;
    return NumGQP;
  }
}

/*!
 * Return the number of Gauss quadrature points on the South
 * cell face at which flux calculation is performed.
 *
 * \param curved_face Set to true if the face is treated as curved otherwise to false.
 */
inline int Grid2D_Quad_Block_HO::NumOfFluxCalculationGaussQuadPoints_South(const int &ii, const int &jj,
									   bool & curved_face){

  if (!IsHighOrderBoundary()){
    // return the number of points based on the order of accuracy
    // if no high-order treatment of geometry is required
    curved_face = false;
    return NumGQP;
  }

  // Handle high-order geometry
  if (ii >= ICl && ii<=ICu){
    // These cells have potentially curved South face
    if (jj == JCu+1){		// check first layer of ghost cells on the North block boundary
      if (BndNorthSplineInfo != NULL){
	/* use high-order boundary info */
	curved_face = true;
	return BndNorthSplineInfo[ii].NumGQPoints();
      } else {
	curved_face = false;
	return NumGQP;
      }
    } else if (jj == JCl){	// check first layer of interior cells on the South block boundary
      if (BndSouthSplineInfo != NULL){
	/* use high-order boundary info */
	curved_face = true;
	return BndSouthSplineInfo[ii].NumGQPoints();
      } else {
	curved_face = false;
	return NumGQP;
      }
    } else {
      curved_face = false;
      return NumGQP;
    }
  } else {
    // These cells don't have any potentially curved South face
    // Return the number of points based on the order of accuracy
    curved_face = false;
    return NumGQP;
  }
}

/*!
 * Return the number of Gauss quadrature points on the East
 * cell face at which flux calculation is performed.
 *
 * \param curved_face Set to true if the face is treated as curved otherwise to false.
 */
inline int Grid2D_Quad_Block_HO::NumOfFluxCalculationGaussQuadPoints_East(const int &ii, const int &jj,
									  bool & curved_face){

  if (!IsHighOrderBoundary()){
    // return the number of points based on the order of accuracy
    // if no high-order treatment of geometry is required
    curved_face = false;
    return NumGQP;
  }

  // Handle high-order geometry
  if (jj >= JCl && jj <= JCu){
    // These cells have potentially curved East face
    if (ii == ICl-1){		// check first layer of ghost cells on the West block boundary
      if (BndWestSplineInfo != NULL){
	/* use high-order boundary info */
	curved_face = true;
	return BndWestSplineInfo[jj].NumGQPoints();
      } else {
	curved_face = false;
	return NumGQP;
      }
    } else if (ii == ICu){	// check last layer of interior cells on the East block boundary
      if (BndEastSplineInfo != NULL){
	/* use high-order boundary info */
	curved_face = true;
	return BndEastSplineInfo[jj].NumGQPoints();
      } else {
	curved_face = false;
	return NumGQP;
      }
    } else {
      curved_face = false;
      return NumGQP;
    }
  } else {
    // These cells don't have any potentially curved East face
    // Return the number of points based on the order of accuracy
    curved_face = false;
    return NumGQP;
  }

}

/*!
 * Return the number of Gauss quadrature points on the West
 * cell face at which flux calculation is performed.
 *
 * \param curved_face Set to true if the face is treated as curved otherwise to false.
 */
inline int Grid2D_Quad_Block_HO::NumOfFluxCalculationGaussQuadPoints_West(const int &ii, const int &jj,
									  bool & curved_face){

  if (!IsHighOrderBoundary()){
    // return the number of points based on the order of accuracy
    // if no high-order treatment of geometry is required
    curved_face = false;
    return NumGQP;
  }

  // Handle high-order geometry
  if (jj >= JCl && jj <= JCu){
    // These cells have potentially curved West face
    if (ii == ICl){		// check first layer of interior cells on the West block boundary
      if (BndWestSplineInfo != NULL){
	/* use high-order boundary info */
	curved_face = true;
	return BndWestSplineInfo[jj].NumGQPoints();
      } else {
	curved_face = false;
	return NumGQP;
      }
    } else if (ii == ICu+1){	// check first layer of ghost cells on the East block boundary
      if (BndEastSplineInfo != NULL){
	/* use high-order boundary info */
	curved_face = true;
	return BndEastSplineInfo[jj].NumGQPoints();
      } else {
	curved_face = false;
	return NumGQP;
      }
    } else {
      curved_face = false;
      return NumGQP;
    }
  } else {
    // These cells don't have any potentially curved East face
    // Return the number of points based on the order of accuracy
    curved_face = false;
    return NumGQP;
  }

}

/*!
 * Return the total number of Gauss quadrature points for the specified cell
 * at which flux calculation is performed.
 */
inline int Grid2D_Quad_Block_HO::NumOfFluxCalculationGaussQuadPoints(const int &ii, const int &jj){
  bool dummy;
  return ( NumOfFluxCalculationGaussQuadPoints_North(ii,jj,dummy) + NumOfFluxCalculationGaussQuadPoints_South(ii,jj,dummy) +
	   NumOfFluxCalculationGaussQuadPoints_West(ii,jj,dummy)  + NumOfFluxCalculationGaussQuadPoints_East(ii,jj,dummy));
}

/*!
 * Create quadrilateral grid block for a Cartesian      
 * mesh defined on a square box with four corner        
 * coordinates (X,Y): (-1,1), (1,1), (1,-1), (-1,1).    
 */
inline void Grid2D_Quad_Block_HO::Create_Quad_Block(const int Number_of_Cells_Idir,
						    const int Number_of_Cells_Jdir,
						    const int Number_of_Ghost_Cells,
						    const int Highest_Order_of_Reconstruction){


  /* Create the 2D quadrilateral grid block without update.
     (i.e. generate the interior nodes and set the boundary condition types)*/
  Create_Quad_Block_Without_Update(Number_of_Cells_Idir,
				   Number_of_Cells_Jdir,
				   Number_of_Ghost_Cells,
				   Highest_Order_of_Reconstruction);

  /* Compute the exterior nodes for the quadrilateral mesh block. */
  Update_Exterior_Nodes();
  
  /* Compute the cells for the quadrilateral mesh block. */
  Update_Cells();
}

/*!
 * Create a 2D quadrilateral grid block with the four   
 * boundaries of the block defined by blended splines   
 * which are given in the input boundary spline file:   
 *                                                      
 * Bnd_Spline_File_Name_ptr. 
 *                           
 */
inline void Grid2D_Quad_Block_HO::Create_Quad_Block(char *Bnd_Spline_File_Name_ptr,
						    const int Number_of_Cells_Idir,
						    const int Number_of_Cells_Jdir,
						    const int Number_of_Ghost_Cells,
						    const int Highest_Order_of_Reconstruction){

  /* Create the 2D quadrilateral grid block without update.
     (i.e. generate the interior nodes and set the boundary condition types)*/
  Create_Quad_Block_Without_Update(Bnd_Spline_File_Name_ptr,
				   Number_of_Cells_Idir,
				   Number_of_Cells_Jdir,
				   Number_of_Ghost_Cells,
				   Highest_Order_of_Reconstruction);

  /* Compute the exterior nodes for the quadrilateral mesh block. */
  Update_Exterior_Nodes();
  
  /* Compute the cells for the quadrilateral mesh block. */
  Update_Cells();

}

/*!
 * Create a 2D quadrilateral grid block with the four   
 * boundaries of the block defined by blended splines   
 * which are given as input arguments to the routine.   
 */
inline void Grid2D_Quad_Block_HO::Create_Quad_Block(Spline2D_HO &Bnd_Spline_North,
						    Spline2D_HO &Bnd_Spline_South,
						    Spline2D_HO &Bnd_Spline_East,
						    Spline2D_HO &Bnd_Spline_West,
						    const int Number_of_Cells_Idir,
						    const int Number_of_Cells_Jdir,
						    const int Number_of_Ghost_Cells,
						    const int Highest_Order_of_Reconstruction,
						    const int Node_Init_Procedure,
						    const int Stretch_I,
						    const double &Beta_I, 
						    const double &Tau_I,
						    const int Stretch_J,
						    const double &Beta_J,
						    const double &Tau_J,
						    const int Orthogonal_North,
						    const int Orthogonal_South,
						    const int Orthogonal_East,
						    const int Orthogonal_West){

  /* Create the 2D quadrilateral grid block without update.
     (i.e. generate the interior nodes and set the boundary condition types)*/
  Create_Quad_Block_Without_Update(Bnd_Spline_North,
				   Bnd_Spline_South, 
				   Bnd_Spline_East,
				   Bnd_Spline_West,
				   Number_of_Cells_Idir,
				   Number_of_Cells_Jdir,
				   Number_of_Ghost_Cells,
				   Highest_Order_of_Reconstruction,
				   Node_Init_Procedure,
				   Stretch_I,
				   Beta_I,
				   Tau_I,
				   Stretch_J,
				   Beta_J,
				   Tau_J,
				   Orthogonal_North,
				   Orthogonal_South,
				   Orthogonal_East,
				   Orthogonal_West);

  /* Compute the exterior nodes for the quadrilateral mesh block. */
  Update_Exterior_Nodes();
  
  /* Compute the cells for the quadrilateral mesh block. */
  Update_Cells();

}

/*!
 * Translates or shifts the positions of the nodes of a 
 * quadrilateral grid block.                            
 */
inline void Grid2D_Quad_Block_HO::Translate_Quad_Block(const Vector2D &V) {
  
  // Translate without update
  Translate_Quad_Block_Without_Update(V);

  /* Compute the cells for the quadrilateral mesh block. */
  Update_Cells();
}

/*!
 * Scales the quadrilateral grid block.                 
 */
inline void Grid2D_Quad_Block_HO::Scale_Quad_Block(const double &Scaling_Factor) {

  // Scale without update
  Scale_Quad_Block_Without_Update(Scaling_Factor);

  /* Compute the cells for the quadrilateral mesh block. */
  Update_Cells();
}

/*!
 * Rotates the quadrilateral grid block.                
 */
inline void Grid2D_Quad_Block_HO::Rotate_Quad_Block(const double &Angle) {
  
  // Rotate without update
  Rotate_Quad_Block_Without_Update(Angle);

  /* Compute the cells for the quadrilateral mesh block. */
  Update_Cells();
}

/*!
 * Re-computes the locations of the nodes and cells of  
 * the quadrilateral grid block based on a mirror       
 * reflection about the y=0 axis.  The cells and nodes  
 * are also re-ordered in the i-direction.              
 */
inline void Grid2D_Quad_Block_HO::Reflect_Quad_Block(void) {

  // Reflect without update
  Reflect_Quad_Block_Without_Update();

  /* Compute the cells for the quadrilateral mesh block. */
  Update_Cells();
}

/*
 * Disturb randomly the interior nodes for the         
 * quadrilateral mesh block. This routine uses the      
 * default seed of the random number generator in order 
 * to determine the angle to which the node is moved.   
 */
inline void Grid2D_Quad_Block_HO::Disturb_Interior_Nodes(const int &Number_of_Iterations) {
  
  // Disturb interior nodes with update
  Disturb_Interior_Nodes_Without_Update(Number_of_Iterations);

  /* Compute the cells for the quadrilateral mesh block. */
  Update_Cells();  
}

/*!
 * Set the number of Gauss quadrature points used for 
 * flux evaluation based on the correlation between the provided  
 * order of reconstruction and the minimum number of GQP required
 * for obtaining the desired solution accuracy.
 *
 * \param [in] RecOrder the reconstruction order required to obtain the targeted solution accuracy
 *
 * \note Case 4 should return NumGQP=3!
 * \todo Find a better way to setup this number based on correlation with the high-order object!
 */
inline void Grid2D_Quad_Block_HO::SetNumberOfFluxCalculationGaussQuadraturePoints(const int & RecOrder){

  switch(RecOrder){
  case 4:
    NumGQP = 2;			// This is a temporary solution to return the write number!
  case 3:
    NumGQP = 2;
    break;
  case 2:
    NumGQP = 2;
    break;
  case 1:
    NumGQP = 1;
    break;
  case 0:
    NumGQP = 1;
    break;
  default:
    throw runtime_error("Grid2D_Quad_Block_HO::SetNumberOfFluxCalculationGaussQuadraturePoints() ERROR! Unknown option for the passed reconstruction order");
  } // endswitch
}

/*!
 * Set the number of Gauss quadrature points used for 
 * flux evaluation based on the correlation between the highest 
 * order of reconstruction and the minimum number of GQP required
 * for obtaining the desired accuracy.
 */
inline void Grid2D_Quad_Block_HO::SetNumberOfGaussQuadraturePoints(void){

  try{
    SetNumberOfFluxCalculationGaussQuadraturePoints(HighestReconstructionOrder);
  } catch (runtime_error){
    throw runtime_error("Grid2D_Quad_Block_HO::SetNumberOfGaussQuadraturePoints() ERROR! Unknown option for the current reconstruction order");
  }
}

/*!
 * Check for existence of curved boundaries
 * and correlate SplineInfo with modifications
 * to the corresponding spline.
 */
inline bool Grid2D_Quad_Block_HO::CheckExistenceOfCurvedBoundaries(void){
  
  // Check for necessity to compute high-order boundary representation
  if ( HighOrderBoundaryRepresentation == OFF ){
    // Update spline interval information (i.e. ensure that no spline interval exits)
    deallocateBndSplineInfo();

    // No curved boundary calculation needed
    return false;
  }

  /* Update spline interval information.
     Delete SplineInfo if the corresponding boundary spline has been deleted or the BC that it carries has been
     converted to a BC that doesn't require curved boundary representation (i.e. BC_NONE or BC_PERIODIC). */
  if ( (BndNorthSpline.Xp == NULL || BndNorthSpline.bc[0] == BC_NONE || BndNorthSpline.bc[0] == BC_PERIODIC ) &&
       BndNorthSplineInfo != NULL){
    delete [] BndNorthSplineInfo; BndNorthSplineInfo = NULL;
  }
  if ( (BndSouthSpline.Xp == NULL || BndSouthSpline.bc[0] == BC_NONE || BndSouthSpline.bc[0] == BC_PERIODIC ) && 
       BndSouthSplineInfo != NULL){
    delete [] BndSouthSplineInfo; BndSouthSplineInfo = NULL;
  }
  if ( (BndEastSpline.Xp == NULL || BndEastSpline.bc[0] == BC_NONE || BndEastSpline.bc[0] == BC_PERIODIC) &&
       BndEastSplineInfo != NULL){
    delete [] BndEastSplineInfo; BndEastSplineInfo = NULL;
  }
  if ( (BndWestSpline.Xp == NULL || BndWestSpline.bc[0] == BC_NONE || BndWestSpline.bc[0] == BC_PERIODIC) &&
       BndWestSplineInfo != NULL){
    delete [] BndWestSplineInfo; BndWestSplineInfo = NULL;
  }

  // Update spline extension interval information
  if ( (ExtendWest_BndNorthSpline.Xp == NULL || 
	ExtendWest_BndNorthSpline.bc[0] == BC_NONE || ExtendWest_BndNorthSpline.bc[0] == BC_PERIODIC ) &&
       ExtendWest_BndNorthSplineInfo != NULL){
    delete [] ExtendWest_BndNorthSplineInfo; ExtendWest_BndNorthSplineInfo = NULL;
  }
  if ( (ExtendEast_BndNorthSpline.Xp == NULL || 
	ExtendEast_BndNorthSpline.bc[0] == BC_NONE || ExtendEast_BndNorthSpline.bc[0] == BC_PERIODIC ) &&
       ExtendEast_BndNorthSplineInfo != NULL){
    delete [] ExtendEast_BndNorthSplineInfo; ExtendEast_BndNorthSplineInfo = NULL;
  }
  if ( (ExtendWest_BndSouthSpline.Xp == NULL || 
	ExtendWest_BndSouthSpline.bc[0] == BC_NONE || ExtendWest_BndSouthSpline.bc[0] == BC_PERIODIC ) &&
       ExtendWest_BndSouthSplineInfo != NULL){
    delete [] ExtendWest_BndSouthSplineInfo; ExtendWest_BndSouthSplineInfo = NULL;
  }
  if ( (ExtendEast_BndSouthSpline.Xp == NULL || 
	ExtendEast_BndSouthSpline.bc[0] == BC_NONE || ExtendEast_BndSouthSpline.bc[0] == BC_PERIODIC ) &&
       ExtendEast_BndSouthSplineInfo != NULL){
    delete [] ExtendEast_BndSouthSplineInfo; ExtendEast_BndSouthSplineInfo = NULL;
  }
  if ( (ExtendNorth_BndEastSpline.Xp == NULL || 
	ExtendNorth_BndEastSpline.bc[0] == BC_NONE || ExtendNorth_BndEastSpline.bc[0] == BC_PERIODIC ) &&
       ExtendNorth_BndEastSplineInfo != NULL){
    delete [] ExtendNorth_BndEastSplineInfo; ExtendNorth_BndEastSplineInfo = NULL;
  }
  if ( (ExtendSouth_BndEastSpline.Xp == NULL || 
	ExtendSouth_BndEastSpline.bc[0] == BC_NONE || ExtendSouth_BndEastSpline.bc[0] == BC_PERIODIC ) &&
       ExtendSouth_BndEastSplineInfo != NULL){
    delete [] ExtendSouth_BndEastSplineInfo; ExtendSouth_BndEastSplineInfo = NULL;
  }
  if ( (ExtendNorth_BndWestSpline.Xp == NULL || 
	ExtendNorth_BndWestSpline.bc[0] == BC_NONE || ExtendNorth_BndWestSpline.bc[0] == BC_PERIODIC ) &&
       ExtendNorth_BndWestSplineInfo != NULL){
    delete [] ExtendNorth_BndWestSplineInfo; ExtendNorth_BndWestSplineInfo = NULL;
  }
  if ( (ExtendSouth_BndWestSpline.Xp == NULL || 
	ExtendSouth_BndWestSpline.bc[0] == BC_NONE || ExtendSouth_BndWestSpline.bc[0] == BC_PERIODIC ) &&
       ExtendSouth_BndWestSplineInfo != NULL){
    delete [] ExtendSouth_BndWestSplineInfo; ExtendSouth_BndWestSplineInfo = NULL;
  }

  
  // Check for nonexistence of curved block boundaries.
  if ( (BndNorthSpline.Xp == NULL || BndNorthSpline.bc[0] == BC_NONE || BndNorthSpline.bc[0] == BC_PERIODIC) &&
       (BndSouthSpline.Xp == NULL || BndSouthSpline.bc[0] == BC_NONE || BndSouthSpline.bc[0] == BC_PERIODIC) &&
       (BndEastSpline.Xp == NULL  || BndEastSpline.bc[0] == BC_NONE || BndEastSpline.bc[0] == BC_PERIODIC) && 
       (BndWestSpline.Xp == NULL  || BndWestSpline.bc[0] == BC_NONE || BndWestSpline.bc[0] == BC_PERIODIC) &&
       (ExtendWest_BndNorthSpline.Xp == NULL || 
	ExtendWest_BndNorthSpline.bc[0] == BC_NONE || ExtendWest_BndNorthSpline.bc[0] == BC_PERIODIC) &&
       (ExtendEast_BndNorthSpline.Xp == NULL || 
	ExtendEast_BndNorthSpline.bc[0] == BC_NONE || ExtendEast_BndNorthSpline.bc[0] == BC_PERIODIC) && 
       (ExtendWest_BndSouthSpline.Xp == NULL || 
	ExtendWest_BndSouthSpline.bc[0] == BC_NONE || ExtendWest_BndSouthSpline.bc[0] == BC_PERIODIC) &&
       (ExtendEast_BndSouthSpline.Xp == NULL || 
	ExtendEast_BndSouthSpline.bc[0] == BC_NONE || ExtendEast_BndSouthSpline.bc[0] == BC_PERIODIC) && 
       (ExtendNorth_BndEastSpline.Xp == NULL || 
	ExtendNorth_BndEastSpline.bc[0] == BC_NONE || ExtendNorth_BndEastSpline.bc[0] == BC_PERIODIC) && 
       (ExtendSouth_BndEastSpline.Xp == NULL || 
	ExtendSouth_BndEastSpline.bc[0] == BC_NONE || ExtendSouth_BndEastSpline.bc[0] == BC_PERIODIC) && 
       (ExtendNorth_BndWestSpline.Xp == NULL || 
	ExtendNorth_BndWestSpline.bc[0] == BC_NONE || ExtendNorth_BndWestSpline.bc[0] == BC_PERIODIC) && 
       (ExtendSouth_BndWestSpline.Xp == NULL || 
	ExtendSouth_BndWestSpline.bc[0] == BC_NONE || ExtendSouth_BndWestSpline.bc[0] == BC_PERIODIC) ){
    // No curved boundaries are present.
    return false;
  }

  // Confirm existence of curved boundaries.
  return true;
}

/*!
 * Check West boundary reconstruction type.
 * Return true if the boundary is curved and the reconstruction
 * along the boundary is constrained to fulfill the boundary conditions,
 * otherwise return false.
 */
inline bool Grid2D_Quad_Block_HO::IsWestBoundaryReconstructionConstrained(void) const{
  return (IsWestBoundaryCurved() && BndWestSpline.getFluxCalcMethod() == ReconstructionBasedFlux)? true : false;
}
  
/*! Check East boundary reconstruction type. */
inline bool Grid2D_Quad_Block_HO::IsEastBoundaryReconstructionConstrained(void) const{
  return (IsEastBoundaryCurved() && BndEastSpline.getFluxCalcMethod() == ReconstructionBasedFlux)? true : false;
}

/*! Check South boundary reconstruction type. */
inline bool Grid2D_Quad_Block_HO::IsSouthBoundaryReconstructionConstrained(void) const{
  return (IsSouthBoundaryCurved() && BndSouthSpline.getFluxCalcMethod() == ReconstructionBasedFlux)? true : false;
}

/*! Check North boundary reconstruction type. */
inline bool Grid2D_Quad_Block_HO::IsNorthBoundaryReconstructionConstrained(void) const{
  return (IsNorthBoundaryCurved() && BndNorthSpline.getFluxCalcMethod() == ReconstructionBasedFlux)? true : false;
}

/*! Check reconstruction type of North extension of West boundary. */
inline bool Grid2D_Quad_Block_HO::IsNorthExtendWestBoundaryReconstructionConstrained(void) const{
  return (IsNorthExtendWestBoundaryCurved() && ExtendNorth_BndWestSpline.getFluxCalcMethod()==ReconstructionBasedFlux)?true:false;
}

/*! Check reconstruction type of South extension of West boundary. */
inline bool Grid2D_Quad_Block_HO::IsSouthExtendWestBoundaryReconstructionConstrained(void) const{
  return (IsSouthExtendWestBoundaryCurved() && ExtendSouth_BndWestSpline.getFluxCalcMethod()==ReconstructionBasedFlux)?true:false;
}

/*! Check reconstruction type of North extension of East boundary. */
inline bool Grid2D_Quad_Block_HO::IsNorthExtendEastBoundaryReconstructionConstrained(void) const{
  return (IsNorthExtendEastBoundaryCurved() && ExtendNorth_BndEastSpline.getFluxCalcMethod()==ReconstructionBasedFlux)?true:false;
}

/*! Check reconstruction type of South extension of East boundary. */
inline bool Grid2D_Quad_Block_HO::IsSouthExtendEastBoundaryReconstructionConstrained(void) const{
  return (IsSouthExtendEastBoundaryCurved() && ExtendSouth_BndEastSpline.getFluxCalcMethod()==ReconstructionBasedFlux)?true:false;
}

/*! Check reconstruction type of East extension of South boundary. */
inline bool Grid2D_Quad_Block_HO::IsEastExtendSouthBoundaryReconstructionConstrained(void) const{
  return (IsEastExtendSouthBoundaryCurved() && ExtendEast_BndSouthSpline.getFluxCalcMethod()==ReconstructionBasedFlux)?true:false;
}

/*! Check reconstruction type of West extension of South boundary. */
inline bool Grid2D_Quad_Block_HO::IsWestExtendSouthBoundaryReconstructionConstrained(void) const{
  return (IsWestExtendSouthBoundaryCurved() && ExtendWest_BndSouthSpline.getFluxCalcMethod()==ReconstructionBasedFlux)?true:false;
}

/*! Check reconstruction type of East extension of North boundary. */
inline bool Grid2D_Quad_Block_HO::IsEastExtendNorthBoundaryReconstructionConstrained(void) const{
  return (IsEastExtendNorthBoundaryCurved() && ExtendEast_BndNorthSpline.getFluxCalcMethod()==ReconstructionBasedFlux)?true:false;
}

/*! Check reconstruction type of West extension of North boundary. */
inline bool Grid2D_Quad_Block_HO::IsWestExtendNorthBoundaryReconstructionConstrained(void) const{
  return (IsWestExtendNorthBoundaryCurved() && ExtendWest_BndNorthSpline.getFluxCalcMethod()==ReconstructionBasedFlux)?true:false;
}

/*!
 * Check if any of the block boundaries is 
 * a solid body one.
 * Return true if at least one boundary is solid,
 * otherwise return false.
 */
inline bool Grid2D_Quad_Block_HO::IsThereAnySolidBoundary(void) const{
  return ( BndNorthSpline.IsSolidBoundary() || BndSouthSpline.IsSolidBoundary() ||
	   BndEastSpline.IsSolidBoundary()  || BndWestSpline.IsSolidBoundary() );
}

/*!
 * Check West boundary spline type (i.e. curved or straight).
 * Return true if the spline has control points
 * and is not an interior boundary based on the BCs 
 * (e.g. BC_NONE).
 */
inline bool Grid2D_Quad_Block_HO::IsWestBoundaryCurved(void) const{
  if (BndWestSpline.Xp != NULL && 
      BndWestSpline.bc[0] != BC_NONE && BndWestSpline.bc[0] != BC_PERIODIC){
    return true;
  } else {
    return false;
  }
}

inline bool Grid2D_Quad_Block_HO::IsEastBoundaryCurved(void) const{
  if (BndEastSpline.Xp != NULL && 
      BndEastSpline.bc[0] != BC_NONE && BndEastSpline.bc[0] != BC_PERIODIC){
    return true;
  } else {
    return false;
  }
}

inline bool Grid2D_Quad_Block_HO::IsSouthBoundaryCurved(void) const{
  if (BndSouthSpline.Xp != NULL && 
      BndSouthSpline.bc[0] != BC_NONE && BndSouthSpline.bc[0] != BC_PERIODIC){
    return true;
  } else {
    return false;
  }
}
 
inline bool Grid2D_Quad_Block_HO::IsNorthBoundaryCurved(void) const{
  if (BndNorthSpline.Xp != NULL && 
      BndNorthSpline.bc[0] != BC_NONE && BndNorthSpline.bc[0] != BC_PERIODIC){
    return true;
  } else {
    return false;
  }
}

inline bool Grid2D_Quad_Block_HO::IsNorthExtendWestBoundaryCurved(void) const{
  if (ExtendNorth_BndWestSpline.Xp != NULL && 
      ExtendNorth_BndWestSpline.bc[0] != BC_NONE && ExtendNorth_BndWestSpline.bc[0] != BC_PERIODIC){
    return true;
  } else {
    return false;
  }
}

inline bool Grid2D_Quad_Block_HO::IsSouthExtendWestBoundaryCurved(void) const{
  if (ExtendSouth_BndWestSpline.Xp != NULL && 
      ExtendSouth_BndWestSpline.bc[0] != BC_NONE && ExtendSouth_BndWestSpline.bc[0] != BC_PERIODIC){
    return true;
  } else {
    return false;
  }
}

inline bool Grid2D_Quad_Block_HO::IsNorthExtendEastBoundaryCurved(void) const{
  if (ExtendNorth_BndEastSpline.Xp != NULL && 
      ExtendNorth_BndEastSpline.bc[0] != BC_NONE && ExtendNorth_BndEastSpline.bc[0] != BC_PERIODIC){
    return true;
  } else {
    return false;
  }
}

inline bool Grid2D_Quad_Block_HO::IsSouthExtendEastBoundaryCurved(void) const{
  if (ExtendSouth_BndEastSpline.Xp != NULL && 
      ExtendSouth_BndEastSpline.bc[0] != BC_NONE && ExtendSouth_BndEastSpline.bc[0] != BC_PERIODIC){
    return true;
  } else {
    return false;
  }
}

inline bool Grid2D_Quad_Block_HO::IsEastExtendSouthBoundaryCurved(void) const{
  if (ExtendEast_BndSouthSpline.Xp != NULL && 
      ExtendEast_BndSouthSpline.bc[0] != BC_NONE && ExtendEast_BndSouthSpline.bc[0] != BC_PERIODIC){
    return true;
  } else {
    return false;
  }
}

inline bool Grid2D_Quad_Block_HO::IsWestExtendSouthBoundaryCurved(void) const{
  if (ExtendWest_BndSouthSpline.Xp != NULL && 
      ExtendWest_BndSouthSpline.bc[0] != BC_NONE && ExtendWest_BndSouthSpline.bc[0] != BC_PERIODIC){
    return true;
  } else {
    return false;
  }
}

inline bool Grid2D_Quad_Block_HO::IsEastExtendNorthBoundaryCurved(void) const{
  if (ExtendEast_BndNorthSpline.Xp != NULL && 
      ExtendEast_BndNorthSpline.bc[0] != BC_NONE && ExtendEast_BndNorthSpline.bc[0] != BC_PERIODIC){
    return true;
  } else {
    return false;
  }
}

inline bool Grid2D_Quad_Block_HO::IsWestExtendNorthBoundaryCurved(void) const{
  if (ExtendWest_BndNorthSpline.Xp != NULL && 
      ExtendWest_BndNorthSpline.bc[0] != BC_NONE && ExtendWest_BndNorthSpline.bc[0] != BC_PERIODIC){
    return true;
  } else {
    return false;
  }
}


/*!
 * Set all the state trackers to a new state.
 */
inline void Grid2D_Quad_Block_HO::New_Global_Geometry_State(void){

  // Increment all trackers
  ++InteriorCellGeometryStateTracker;
  ++GhostCellGeometryStateTracker;
  ++CornerGhostCellGeometryStateTracker;
}

/*!
 * Set interior geometry state tracker to a new state.
 */
inline void Grid2D_Quad_Block_HO::New_Interior_Geometry_State(void){
  ++InteriorCellGeometryStateTracker;
}

/*!
 * Set ghost cell geometry state tracker to a new state.
 */
inline void Grid2D_Quad_Block_HO::New_Ghost_Geometry_State(void){
  ++GhostCellGeometryStateTracker;
}


/*!
 * Set corner ghost cell geometry state tracker to a new state.
 */
inline void Grid2D_Quad_Block_HO::New_Corner_Geometry_State(void){
  ++CornerGhostCellGeometryStateTracker;
}


/* ---------------------------------------------------------------------------------------------- 
 * ===============        INCLUDE TEMPLATE SPECIALIZATIONS FOR THIS CLASS      ==================
 * ---------------------------------------------------------------------------------------------*/
#include "HO_Grid2DQuad_Specializations.h"


#endif /* _GRID2D_QUAD_BLOCK_INCLUDED  */

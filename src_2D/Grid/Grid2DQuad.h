/* Grid2DQuad.h:  Header file defining 2D quadrilateral block grid type. */

#ifndef _GRID2D_QUAD_BLOCK_INCLUDED
#define _GRID2D_QUAD_BLOCK_INCLUDED

/* Include required C++ libraries. */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>

using namespace std;

/* Include math macro, CFD, 2D vector, 2D cell,
   2D spline, and MPI header files. */

#ifndef _MATH_MACROS_INCLUDED
#include "../Math/Math.h"
#endif // _MATH_MACROS_INCLUDED

#ifndef _CFD_INCLUDED
#include "../CFD/CFD.h"
#endif // _CFD_INCLUDED

#ifndef _VECTOR2D_INCLUDED
#include "../Math/Vector2D.h"
#endif //_VECTOR2D_INCLUDED

#ifndef _SPLINE2D_INCLUDED
#include "../Math/Spline2D.h"
#endif // _SPLINE2D_INCLUDED

#ifndef _CELL2D_INCLUDED
#include "Cell2D.h"
#endif // _CELL2D_INCLUDED

#ifndef _MPI_INCLUDED
#include "../MPI/MPI.h"
#endif // _MPI_INCLUDED

// Include the linear systems header file.

#ifndef _LINEARSYSTEMS_INCLUDED
#include "../Math/LinearSystems.h"
#endif // _LINEARSYSTEMS_INCLUDED

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

/* Define the quadrilateral 2D grid block class. */

/*!
 * Class: Grid2D_Quad_Block
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
 *      allocateNodes -- Allocate memory for nodes of
 *                       structured quadrilateral grid
 *                       block.
 *      deallocateNodes -- Deallocate memory for nodes
 *                         of structured quadrilateral
 *                         grid block.
 *      allocateCells -- Allocate memory for nodes of
 *                       structured quadrilateral grid
 *                       block.
 *      deallocateCells -- Deallocate memory for nodes
 *                         of structured quadrilateral
 *                         grid block.
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
class Grid2D_Quad_Block{
  private:
  public:
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
    Node2D             **Node; //!< Array of 2D node position vectors
    Cell2D             **Cell; //!< Array of 2D cell centre vectors
    int              *BCtypeN, //!< North boundary condition type specifier.
                     *BCtypeS, //!< South boundary condition type specifier.
                     *BCtypeE, //!< East boundary condition type specifier.
                     *BCtypeW; //!< West boundary condition type specifier.
    Spline2D   BndNorthSpline, //!< North boundary 2D spline.
               BndSouthSpline, //!< South boundary 2D spline.
                BndEastSpline, //!< East boundary 2D spline.
                BndWestSpline; //!< West boundary 2D spline.
    double              SminN, //!< Minimum value of north face pathlength.
                        SmaxN, //!< Maximum value of north face pathlength.
                        SminS, //!< Minimum value of south face pathlength.
                        SmaxS, //!< Maximum value of south face pathlength.
                        SminE, //!< Minimum value of east face pathlength.
                        SmaxE, //!< Maximum value of east face pathlength.
                        SminW, //!< Minimum value of west face pathlength.
                        SmaxW; //!< Maximum value of west face pathlength.
    int              StretchI, //!< Grid point stretching function (i-direction).
                     StretchJ; //!< Grid point stretching function (j-direction).
    double              BetaI, //!< Grid point stretching function parameter (i-direction).
                         TauI, //!< Grid point stretching function parameter (i-direction).
                        BetaJ, //!< Grid point stretching function parameter (j-direction).
                         TauJ; //!< Grid point stretching function parameter (j-direction).
    int           OrthogonalN, //!< North boundary orthogonality parameter.
                  OrthogonalS, //!< South boundary orthogonality parameter.
                  OrthogonalE, //!< East boundary orthogonality parameter.
                  OrthogonalW; //!< West boundary orthogonality parameter.
                               // Made public so can access them.
			
    //@{ @name Constructors.

    //! Creation constructor.
    Grid2D_Quad_Block(void) {
       NNi = 0; INl = 0; INu = 0; NNj = 0; JNl = 0; JNu = 0;
       NCi = 0; ICl = 0; ICu = 0; NCj = 0; JCl = 0; JCu = 0;
       Nghost = 0;
       Node = NULL; Cell = NULL;
       BCtypeN = NULL; BCtypeS = NULL; BCtypeE = NULL; BCtypeW = NULL;
       SminN = ZERO; SmaxN = ZERO; SminS = ZERO; SmaxS = ZERO; 
       SminE = ZERO; SmaxE = ZERO; SminW = ZERO; SmaxW = ZERO;
       StretchI = 0; StretchJ = 0; BetaI = ONE; TauI = ONE;
       BetaJ = ONE; TauJ = ONE;
       OrthogonalN = 1; OrthogonalS = 1; OrthogonalE = 1; OrthogonalW = 1;
    }

    //! Copy constructor.
    Grid2D_Quad_Block(const Grid2D_Quad_Block &G) {
       NNi = G.NNi; INl = G.INl; INu = G.INu; 
       NNj = G.NNj; JNl = G.JNl; JNu = G.JNu;
       NCi = G.NCi; ICl = G.ICl; ICu = G.ICu; 
       NCj = G.NCj; JCl = G.JCl; JCu = G.JCu;
       Nghost = G.Nghost;
       Node = G.Node; Cell = G.Cell;
       BCtypeN = G.BCtypeN; BCtypeS = G.BCtypeS; 
       BCtypeE = G.BCtypeE; BCtypeW = G.BCtypeW;
       BndNorthSpline = G.BndNorthSpline; BndSouthSpline = G.BndSouthSpline;
       BndEastSpline = G.BndEastSpline; BndWestSpline = G.BndWestSpline;
       SminN = G.SminN; SmaxN = G.SmaxN; SminS = G.SminS; SmaxS = G.SmaxS; 
       SminE = G.SminE; SmaxE = G.SmaxE; SminW = G.SminW; SmaxW = G.SmaxW;
       StretchI = G.StretchI; StretchJ = G.StretchJ; BetaI = G.BetaI; 
       TauI = G.TauI; BetaJ = G.BetaJ; TauJ = G.TauJ;
       OrthogonalN = G.OrthogonalN; OrthogonalS = G.OrthogonalS; 
       OrthogonalE = G.OrthogonalE; OrthogonalW = G.OrthogonalW;
    }

    //@}

    /* Destructor. */
    // ~Grid2D_Quad_Block(void);
    // Use automatically generated destructor.

    /* Assignment operator. */
    // Grid2D_Quad_Block operator = (const Grid2D_Quad_Block &Soln);
    // Use automatically generated assignment operator.

    //@{ @name Allocation and Deallocation

    //! Allocate memory for structured quadrilateral grid block.
    void allocate(const int Ni, const int Nj, const int Ng);

    //! Deallocate memory for structured quadrilateral grid block.
    void deallocate(void);

    //! Allocate memory for nodes of structured quadrilateral grid block.
    void allocateNodes(const int Ni, const int Nj, const int Ng);

    //! Deallocate memory for nodes of structured quadrilateral grid block.
    void deallocateNodes(void);

    //! Allocate memory for cells of structured quadrilateral grid block.
    void allocateCells(const int Ni, const int Nj, const int Ng);

    //! Deallocate memory for cells of structured quadrilateral grid block.
    void deallocateCells(void);

    //@}

    //@{ @name Calculate centroid of cell.
    Vector2D centroid(const Cell2D &Cell);
    Vector2D centroid(const int ii, const int jj);
    Vector2D centroidSW(const int ii, const int jj);
    Vector2D centroidSE(const int ii, const int jj);
    Vector2D centroidNW(const int ii, const int jj);
    Vector2D centroidNE(const int ii, const int jj);
    //@}

    //@{ @name Calculate cell area.
    double area(const Cell2D &Cell);
    double area(const int ii, const int jj);
    //@}

    //@{ @name Get cell nodes.
    Node2D nodeNW(const Cell2D &Cell);
    Node2D nodeNW(const int ii, const int jj);

    Node2D nodeNE(const Cell2D &Cell);
    Node2D nodeNE(const int ii, const int jj);

    Node2D nodeSE(const Cell2D &Cell);
    Node2D nodeSE(const int ii, const int jj);

    Node2D nodeSW(const Cell2D &Cell);
    Node2D nodeSW(const int ii, const int jj);
    //@}

    //@{ @name Get cell face midpoints.
    Vector2D xfaceN(const Cell2D &Cell);
    Vector2D xfaceN(const int ii, const int jj);

    Vector2D xfaceS(const Cell2D &Cell);
    Vector2D xfaceS(const int ii, const int jj);

    Vector2D xfaceE(const Cell2D &Cell);
    Vector2D xfaceE(const int ii, const int jj);

    Vector2D xfaceW(const Cell2D &Cell);
    Vector2D xfaceW(const int ii, const int jj);
    //@}

    //@{ @name Get cell face lengths.
    double lfaceN(const Cell2D &Cell);
    double lfaceN(const int ii, const int jj);

    double lfaceS(const Cell2D &Cell);
    double lfaceS(const int ii, const int jj);

    double lfaceE(const Cell2D &Cell);
    double lfaceE(const int ii, const int jj);

    double lfaceW(const Cell2D &Cell);
    double lfaceW(const int ii, const int jj);
    //@}

    //@{ @name Get the unit vector normal to the cell faces.
    Vector2D nfaceN(const Cell2D &Cell);
    Vector2D nfaceN(const int ii, const int jj);

    Vector2D nfaceS(const Cell2D &Cell);
    Vector2D nfaceS(const int ii, const int jj);

    Vector2D nfaceE(const Cell2D &Cell);
    Vector2D nfaceE(const int ii, const int jj);

    Vector2D nfaceW(const Cell2D &Cell);
    Vector2D nfaceW(const int ii, const int jj);
    //@}

    /* Assignment operator. */
    // Grid2D_Quad_Block operator = 
    //    (const Grid2D_Quad_Block &G);
    // Use automatically generated assignment operator.

    //@{ @name Binary arithmetic operators.
    friend Grid2D_Quad_Block operator +(Grid2D_Quad_Block &G, const Vector2D &V);
    friend Grid2D_Quad_Block operator -(Grid2D_Quad_Block &G, const Vector2D &V);
    friend Grid2D_Quad_Block operator *(Grid2D_Quad_Block &G, const double &a);
    friend Grid2D_Quad_Block operator *(const double &a, Grid2D_Quad_Block &G);
    friend Grid2D_Quad_Block operator /(Grid2D_Quad_Block &G, const double &a);
    friend Grid2D_Quad_Block operator ^(Grid2D_Quad_Block &G, const double &a);
    //@}

    //@{ @name Input-output operators.
    friend ostream &operator << (ostream &out_file, const Grid2D_Quad_Block &G);
    friend istream &operator >> (istream &in_file, Grid2D_Quad_Block &G);
    //@}

};

/*************************************************************************
 * Grid2D_Quad_Block::allocate -- Allocate memory.                       *
 *************************************************************************/
inline void Grid2D_Quad_Block::allocate(const int Ni, const int Nj, const int Ng) {
   int i; assert( Ni > 1 && Nj > 1 && Ng >= 2);
   NNi = Ni+2*Ng+1; INl = Ng; INu = Ni+Ng; NNj = Nj+2*Ng+1; JNl = Ng; JNu = Nj+Ng;
   NCi = Ni+2*Ng; ICl = Ng; ICu = Ni+Ng-1; NCj = Nj+2*Ng; JCl = Ng; JCu = Nj+Ng-1;
   Nghost = Ng;
   Node = new Node2D*[NNi];
   for ( i = 0; i <= NNi-1 ; ++i ) Node[i] = new Node2D[NNj];
   Cell = new Cell2D*[NCi];
   for ( i = 0; i <= NCi-1 ; ++i ) Cell[i] = new Cell2D[NCj];
   BCtypeN = new int[NCi]; BCtypeS = new int[NCi];
   BCtypeE = new int[NCj]; BCtypeW = new int[NCj];
}

/*************************************************************************
 * Grid2D_Quad_Block::deallocate -- Deallocate memory.                   *
 *************************************************************************/
inline void Grid2D_Quad_Block::deallocate(void) {
   int i;
   for ( i = 0; i <= NNi-1 ; ++i ) {
      delete []Node[i]; Node[i] = NULL;
   } /* endfor */
   delete []Node; Node = NULL;
   for ( i = 0; i <= NCi-1 ; ++i ) {
      delete []Cell[i]; Cell[i] = NULL;
   } /* endfor */
   delete []Cell; Cell = NULL; 
   delete []BCtypeN; BCtypeN = NULL; delete []BCtypeS; BCtypeS = NULL; 
   delete []BCtypeE; BCtypeE = NULL; delete []BCtypeW; BCtypeW = NULL;
   NNi = 0; INl = 0; INu = 0; NNj = 0; JNl = 0; JNu = 0;
   NCi = 0; ICl = 0; ICu = 0; NCj = 0; JCl = 0; JCu = 0;
   Nghost = 0;
   BndNorthSpline.deallocate(); BndSouthSpline.deallocate();
   BndEastSpline.deallocate(); BndWestSpline.deallocate();
   StretchI = 0; StretchJ = 0; BetaI = ONE; TauI = ONE;
   BetaJ = ONE; TauJ = ONE;
   OrthogonalN = 1; OrthogonalS = 1; OrthogonalE = 1; OrthogonalW = 1;
}

/*************************************************************************
 * Grid2D_Quad_Block::allocateNodes -- Allocate nodes.                   *
 *************************************************************************/
inline void Grid2D_Quad_Block::allocateNodes(const int Ni, const int Nj, const int Ng) {
   int i; assert( Ni > 1 && Nj > 1 && Ng >= 2);
   NNi = Ni+2*Ng+1; INl = Ng; INu = Ni+Ng;
   NNj = Nj+2*Ng+1; JNl = Ng; JNu = Nj+Ng;
   Node = new Node2D*[NNi];
   for ( i = 0; i <= NNi-1 ; ++i ) Node[i] = new Node2D[NNj];
}

/*************************************************************************
 * Grid2D_Quad_Block::deallocateNodes -- Deallocate Nodes.               *
 *************************************************************************/
inline void Grid2D_Quad_Block::deallocateNodes(void) {
   int i;
   for ( i = 0; i <= NNi-1 ; ++i ) {
      delete []Node[i]; Node[i] = NULL;
   } /* endfor */
   delete []Node; Node = NULL;
   NNi = 0; INl = 0; INu = 0; NNj = 0; JNl = 0; JNu = 0; Nghost = 0;
}

/*************************************************************************
 * Grid2D_Quad_Block::allocateCells -- Allocate cells.                   *
 *************************************************************************/
inline void Grid2D_Quad_Block::allocateCells(const int Ni, const int Nj, const int Ng) {
   int i; assert( Ni > 1 && Nj > 1 && Ng >= 2);
   NCi = Ni+2*Ng; ICl = Ng; ICu = Ni+Ng-1;
   NCj = Nj+2*Ng; JCl = Ng; JCu = Nj+Ng-1;
   Nghost = Ng;
   Cell = new Cell2D*[NNi];
   for ( i = 0; i <= NCi-1 ; ++i ) Cell[i] = new Cell2D[NCj];
   BCtypeN = new int[NCi]; BCtypeS = new int[NCi];
   BCtypeE = new int[NCj]; BCtypeW = new int[NCj];
}

/*************************************************************************
 * Grid2D_Quad_Block::deallocateCells -- Deallocate Cells.               *
 *************************************************************************/
inline void Grid2D_Quad_Block::deallocateCells(void) {
   int i;
   for ( i = 0; i <= NCi-1 ; ++i ) {
      delete []Cell[i]; Cell[i] = NULL;
   } /* endfor */
   delete []Cell; Cell = NULL; 
   delete []BCtypeN; BCtypeN = NULL; delete []BCtypeS; BCtypeS = NULL; 
   delete []BCtypeE; BCtypeE = NULL; delete []BCtypeW; BCtypeW = NULL;
   NCi = 0; ICl = 0; ICu = 0; NCj = 0; JCl = 0; JCu = 0; Nghost = 0;
}

/**********************************************************************
 * Grid2D_Quad_Block::centroid -- Cell centre.                        *
 *                                                                    *
 * The centroid of a triangular cell can be determined from the       *
 * average of the three distinct vertices.  However, this method is   *
 * not valid for all quadrilateral cells, especially not skewed or    *
 * oddly shaped quadrilaterals.  For these cells, the centroid can    *
 * be determined by the area-weighted average of the centroids of the *
 * sub-triangles.                                                     *
 *                                                                    *
 **********************************************************************/
inline Vector2D Grid2D_Quad_Block::centroid(const Cell2D &Cell) {
  return centroid(Cell.I,Cell.J);
}

inline Vector2D Grid2D_Quad_Block::centroid(const int ii, const int jj) {
  Vector2D X1, X2, X3, X4, Xc1, Xc2, X;
  double A1, A2;
  // Cell nodes in counter-clockwise order.
  X1 = Node[ii  ][jj  ].X;
  X2 = Node[ii+1][jj  ].X;
  X3 = Node[ii+1][jj+1].X;
  X4 = Node[ii  ][jj+1].X;
  // Determine the centroid and area of the sub-triangles.
  Xc1 = (X1+X2+X3)/3.0;
  Xc2 = (X1+X3+X4)/3.0;
//   A1 = HALF*((X1^X2) + (X2^X3) + (X3^X1));
//   A2 = HALF*((X1^X3) + (X3^X4) + (X4^X1));
  A1 = HALF*((X2-X1)^(X3-X1));
  A2 = HALF*((X3-X4)^(X3-X2));
  // Return the area-weighted average of the centroids of the sub-triangles:
  if (A1 > ZERO && A2 > ZERO) return (A1*Xc1 + A2*Xc2)/(A1+A2);
  // Average of four nodes (not always correct):
  return 0.25*(Node[ii][jj].X + Node[ii+1][jj].X + Node[ii+1][jj+1].X + Node[ii][jj+1].X);
}

inline Vector2D Grid2D_Quad_Block::centroidSW(const int ii, const int jj) {
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
  A2 = HALF*((X3-X4)^(X3-X2));
  // Return the area-weighted average of the centroids of the sub-triangles:
  if (A1 > ZERO && A2 > ZERO) return (A1*Xc1 + A2*Xc2)/(A1+A2);
  // Average of four nodes (not always correct):
  return (X1 + X2 + X3 + X4)/4.0;
}

inline Vector2D Grid2D_Quad_Block::centroidSE(const int ii, const int jj) {
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
  A2 = HALF*((X3-X4)^(X3-X2));
  // Return the area-weighted average of the centroids of the sub-triangles:
  if (A1 > ZERO && A2 > ZERO) return (A1*Xc1 + A2*Xc2)/(A1+A2);
  // Average of four nodes (not always correct):
  return (X1 + X2 + X3 + X4)/4.0;
}

inline Vector2D Grid2D_Quad_Block::centroidNW(const int ii, const int jj) {
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
  A2 = HALF*((X3-X4)^(X3-X2));
  // Return the area-weighted average of the centroids of the sub-triangles:
  if (A1 > ZERO && A2 > ZERO) return (A1*Xc1 + A2*Xc2)/(A1+A2);
  // Average of four nodes (not always correct):
  return (X1 + X2 + X3 + X4)/4.0;
}

inline Vector2D Grid2D_Quad_Block::centroidNE(const int ii, const int jj) {
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
  A2 = HALF*((X3-X4)^(X3-X2));
  // Return the area-weighted average of the centroids of the sub-triangles:
  if (A1 > ZERO && A2 > ZERO) return (A1*Xc1 + A2*Xc2)/(A1+A2);
  // Average of four nodes (not always correct):
  return (X1 + X2 + X3 + X4)/4.0;
}

/*************************************************************************
 * Grid2D_Quad_Block::area -- Cell area.                                 *
 *************************************************************************/
inline double Grid2D_Quad_Block::area(const Cell2D &Cell) {
    return (HALF*(((Node[Cell.I+1][Cell.J].X-Node[Cell.I][Cell.J].X)^(
                    Node[Cell.I][Cell.J+1].X-Node[Cell.I][Cell.J].X)) +
                  ((Node[Cell.I+1][Cell.J+1].X-Node[Cell.I][Cell.J+1].X)^(
                    Node[Cell.I+1][Cell.J+1].X-Node[Cell.I+1][Cell.J].X))));
}

inline double Grid2D_Quad_Block::area(const int ii, const int jj) {
    return (HALF*(((Node[ii+1][jj].X-Node[ii][jj].X)^(
                    Node[ii][jj+1].X-Node[ii][jj].X)) +
                  ((Node[ii+1][jj+1].X-Node[ii][jj+1].X)^(
                    Node[ii+1][jj+1].X-Node[ii+1][jj].X))));
}

/*************************************************************************
 * Grid2D_Quad_Block::node?? -- Get cell nodes.                          *
 *************************************************************************/
inline Node2D Grid2D_Quad_Block::nodeNW(const Cell2D &Cell) {
    return (Node[Cell.I][Cell.J+1]);
}

inline Node2D Grid2D_Quad_Block::nodeNW(const int ii, const int jj) {
    return (Node[ii][jj+1]);
}

inline Node2D Grid2D_Quad_Block::nodeNE(const Cell2D &Cell) {
    return (Node[Cell.I+1][Cell.J+1]);
}

inline Node2D Grid2D_Quad_Block::nodeNE(const int ii, const int jj) {
    return (Node[ii+1][jj+1]);
}

inline Node2D Grid2D_Quad_Block::nodeSE(const Cell2D &Cell) {
    return (Node[Cell.I+1][Cell.J]);
}

inline Node2D Grid2D_Quad_Block::nodeSE(const int ii, const int jj) {
    return (Node[ii+1][jj]);
}

inline Node2D Grid2D_Quad_Block::nodeSW(const Cell2D &Cell) {
    return (Node[Cell.I][Cell.J]);
}

inline Node2D Grid2D_Quad_Block::nodeSW(const int ii, const int jj) {
    return (Node[ii][jj]);
}

/*************************************************************************
 * Grid2D_Quad_Block::xface? -- Get face midpoints.                      *
 *************************************************************************/
inline Vector2D Grid2D_Quad_Block::xfaceN(const Cell2D &Cell) {
    return (HALF*(Node[Cell.I][Cell.J+1].X+Node[Cell.I+1][Cell.J+1].X));
}

inline Vector2D Grid2D_Quad_Block::xfaceN(const int ii, const int jj) {
    return (HALF*(Node[ii][jj+1].X+Node[ii+1][jj+1].X));
}

inline Vector2D Grid2D_Quad_Block::xfaceS(const Cell2D &Cell) {
    return (HALF*(Node[Cell.I][Cell.J].X+Node[Cell.I+1][Cell.J].X));
}

inline Vector2D Grid2D_Quad_Block::xfaceS(const int ii, const int jj) {
    return (HALF*(Node[ii][jj].X+Node[ii+1][jj].X));
}

inline Vector2D Grid2D_Quad_Block::xfaceE(const Cell2D &Cell) {
    return (HALF*(Node[Cell.I+1][Cell.J].X+Node[Cell.I+1][Cell.J+1].X));
}

inline Vector2D Grid2D_Quad_Block::xfaceE(const int ii, const int jj) {
    return (HALF*(Node[ii+1][jj].X+Node[ii+1][jj+1].X));
}

inline Vector2D Grid2D_Quad_Block::xfaceW(const Cell2D &Cell) {
    return (HALF*(Node[Cell.I][Cell.J].X+Node[Cell.I][Cell.J+1].X));
}

inline Vector2D Grid2D_Quad_Block::xfaceW(const int ii, const int jj) {
    return (HALF*(Node[ii][jj].X+Node[ii][jj+1].X));
}

/*************************************************************************
 * Grid2D_Quad_Block::lface? -- Get face lengths.                        *
 *************************************************************************/
inline double Grid2D_Quad_Block::lfaceN(const Cell2D &Cell) {
    return (abs(Node[Cell.I][Cell.J+1].X-Node[Cell.I+1][Cell.J+1].X));
}

inline double Grid2D_Quad_Block::lfaceN(const int ii, const int jj) {
    return (abs(Node[ii][jj+1].X-Node[ii+1][jj+1].X));
}

inline double Grid2D_Quad_Block::lfaceS(const Cell2D &Cell) {
    return (abs(Node[Cell.I][Cell.J].X-Node[Cell.I+1][Cell.J].X));
}

inline double Grid2D_Quad_Block::lfaceS(const int ii, const int jj) {
    return (abs(Node[ii][jj].X-Node[ii+1][jj].X));
}

inline double Grid2D_Quad_Block::lfaceE(const Cell2D &Cell) {
    return (abs(Node[Cell.I+1][Cell.J].X-Node[Cell.I+1][Cell.J+1].X));
}

inline double Grid2D_Quad_Block::lfaceE(const int ii, const int jj) {
    return (abs(Node[ii+1][jj].X-Node[ii+1][jj+1].X));
}

inline double Grid2D_Quad_Block::lfaceW(const Cell2D &Cell) {
    return (abs(Node[Cell.I][Cell.J].X-Node[Cell.I][Cell.J+1].X));
}

inline double Grid2D_Quad_Block::lfaceW(const int ii, const int jj) {
    return (abs(Node[ii][jj].X-Node[ii][jj+1].X));
}

/*************************************************************************
 * Grid2D_Quad_Block::nface? -- Get cell face normals.                   *
 *************************************************************************/
inline Vector2D Grid2D_Quad_Block::nfaceN(const Cell2D &Cell) {
  return nfaceN(Cell.I,Cell.J);
}

inline Vector2D Grid2D_Quad_Block::nfaceN(const int ii, const int jj) {
  return (Vector2D( (Node[ii][jj+1].X.y - Node[ii+1][jj+1].X.y),
		    -(Node[ii][jj+1].X.x - Node[ii+1][jj+1].X.x))/
	  abs(Node[ii][jj+1].X - Node[ii+1][jj+1].X));
  if (lfaceN(ii,jj) > NANO) return (Vector2D( (Node[ii][jj+1].X.y - Node[ii+1][jj+1].X.y),
					      -(Node[ii][jj+1].X.x - Node[ii+1][jj+1].X.x))/
				    abs(Node[ii][jj+1].X - Node[ii+1][jj+1].X));
  if (ii > 0 && ii < NCi-1) {
    if (lfaceN(ii-1,jj) > ZERO && lfaceN(ii+1,jj) > ZERO) {
      return HALF*(nfaceN(ii-1,jj) + nfaceN(ii+1,jj));
    } else if (lfaceN(ii-1,jj) > ZERO) {
      return nfaceN(ii-1,jj);
    } else if (lfaceN(ii+1,jj) > ZERO) {
      return nfaceN(ii+1,jj);
    } else {
      assert(1==2);
    }
  } else if (ii > 0) {
    if (lfaceN(ii-1,jj) > ZERO) return nfaceN(ii-1,jj);
    else assert(1==2);
  } else if (ii < NCi-1) {
    if (lfaceN(ii+1,jj) > ZERO) return nfaceN(ii+1,jj);
    else assert(1==2);
  } else {
    assert(1==2);
  }
}

inline Vector2D Grid2D_Quad_Block::nfaceS(const Cell2D &Cell) {
  return nfaceS(Cell.I,Cell.J);
}

inline Vector2D Grid2D_Quad_Block::nfaceS(const int ii, const int jj) {
  return (Vector2D( (Node[ii+1][jj].X.y - Node[ii][jj].X.y),
		    -(Node[ii+1][jj].X.x - Node[ii][jj].X.x))/
	  abs(Node[ii+1][jj].X - Node[ii][jj].X));
  if (lfaceS(ii,jj) > NANO) return (Vector2D( (Node[ii+1][jj].X.y - Node[ii][jj].X.y),
					      -(Node[ii+1][jj].X.x - Node[ii][jj].X.x))/
				    abs(Node[ii+1][jj].X - Node[ii][jj].X));
  if (ii > 0 && ii < NCi-1) {
    if (lfaceS(ii-1,jj) > ZERO && lfaceS(ii+1,jj) > ZERO) {
      return HALF*(nfaceS(ii-1,jj) + nfaceS(ii+1,jj));
    } else if (lfaceS(ii-1,jj) > ZERO) {
      return nfaceS(ii-1,jj);
    } else if (lfaceS(ii+1,jj) > ZERO) {
      return nfaceS(ii+1,jj);
    } else {
      assert(1==2);
    }
  } else if (ii > 0) {
    if (lfaceS(ii-1,jj) > ZERO) return nfaceS(ii-1,jj);
    else assert(1==2);
  } else if (ii < NCi-1) {
    if (lfaceS(ii+1,jj) > ZERO) return nfaceS(ii+1,jj);
    else assert(1==2);
  } else {
    assert(1==2);
  }
}

inline Vector2D Grid2D_Quad_Block::nfaceE(const Cell2D &Cell) {
  return nfaceE(Cell.I,Cell.J);
}

inline Vector2D Grid2D_Quad_Block::nfaceE(const int ii, const int jj) {
  return (Vector2D( (Node[ii+1][jj+1].X.y - Node[ii+1][jj].X.y),
		    -(Node[ii+1][jj+1].X.x - Node[ii+1][jj].X.x))/
	  abs(Node[ii+1][jj+1].X-Node[ii+1][jj].X));
  if (lfaceE(ii,jj) > NANO) return (Vector2D( (Node[ii+1][jj+1].X.y - Node[ii+1][jj].X.y),
					      -(Node[ii+1][jj+1].X.x - Node[ii+1][jj].X.x))/
				    abs(Node[ii+1][jj+1].X-Node[ii+1][jj].X));
  if (jj > 0 && jj < NCj-1) {
    if (lfaceE(ii,jj-1) > ZERO && lfaceE(ii,jj+1) > ZERO) {
      return HALF*(nfaceE(ii,jj-1) + nfaceE(ii,jj+1));
    } else if (lfaceE(ii,jj-1) > ZERO) {
      return nfaceE(ii,jj-1);
    } else if (lfaceE(ii,jj+1) > ZERO) {
      return nfaceE(ii,jj+1);
    } else {
      assert(1==2);
    }
  } else if (jj > 0) {
    if (lfaceE(ii,jj-1) > ZERO) return nfaceE(ii,jj-1);
    else assert(1==2);
  } else if (jj < NCj-1) {
    if (lfaceE(ii,jj+1) > ZERO) return nfaceE(ii,jj+1);
    else assert(1==2);
  } else {
    assert(1==2);
  }
}

inline Vector2D Grid2D_Quad_Block::nfaceW(const Cell2D &Cell) {
  return nfaceW(Cell.I,Cell.J);
}

inline Vector2D Grid2D_Quad_Block::nfaceW(const int ii, const int jj) {
  return (Vector2D( (Node[ii][jj].X.y - Node[ii][jj+1].X.y),
 		    -(Node[ii][jj].X.x - Node[ii][jj+1].X.x))/
 	  abs(Node[ii][jj].X - Node[ii][jj+1].X));
  if (lfaceW(ii,jj) > NANO) return (Vector2D( (Node[ii][jj].X.y - Node[ii][jj+1].X.y),
					      -(Node[ii][jj].X.x - Node[ii][jj+1].X.x))/
				    abs(Node[ii][jj].X - Node[ii][jj+1].X));
  if (jj > 0 && jj < NCj-1) {
    if (lfaceW(ii,jj-1) > ZERO && lfaceW(ii,jj+1) > ZERO) {
      return HALF*(nfaceW(ii,jj-1) + nfaceW(ii,jj+1));
    } else if (lfaceW(ii,jj-1) > ZERO) {
      return nfaceW(ii,jj-1);
    } else if (lfaceW(ii,jj+1) > ZERO) {
      return nfaceW(ii,jj+1);
    } else {
      assert(1==2);
    }
  } else if (jj > 0) {
    if (lfaceW(ii,jj-1) > ZERO) return nfaceW(ii,jj-1);
    else assert(1==2);
  } else if (jj < NCj-1) {
    if (lfaceW(ii,jj+1) > ZERO) return nfaceW(ii,jj+1);
    else assert(1==2);
  } else {
    assert(1==2);
  }
}

/********************************************************
 * Grid2D_Quad_Block -- Binary arithmetic operators.    *
 ********************************************************/
// Shift operators.
inline Grid2D_Quad_Block operator +(Grid2D_Quad_Block &G, const Vector2D &V) {
  int i, j;
  for ( j = G.JNl-G.Nghost ; j <= G.JNu+G.Nghost; ++j ) {
     for ( i = G.INl-G.Nghost ; i <= G.INu+G.Nghost; ++i ) {
         G.Node[i][j].X += V;
     } /* endfor */
  } /* endfor */
  for ( j = G.JCl-G.Nghost ; j <= G.JCu+G.Nghost ; ++j) {
      for ( i = G.ICl-G.Nghost ; i <= G.ICu+G.Nghost ; ++i) {
  	  G.Cell[i][j].Xc = G.centroid(i, j);
          G.Cell[i][j].A = G.area(i, j);
      } /* endfor */
  } /* endfor */
  if (G.BndNorthSpline.np != 0 ) G.BndNorthSpline = G.BndNorthSpline + V;
  if (G.BndSouthSpline.np != 0 ) G.BndSouthSpline = G.BndSouthSpline + V;
  if (G.BndEastSpline.np != 0 ) G.BndEastSpline = G.BndEastSpline + V;
  if (G.BndWestSpline.np != 0 ) G.BndWestSpline = G.BndWestSpline + V;
  return (G);
}

inline Grid2D_Quad_Block operator -(Grid2D_Quad_Block &G, const Vector2D &V) {
  int i, j;
  for ( j = G.JNl-G.Nghost ; j <= G.JNu+G.Nghost; ++j ) {
     for ( i = G.INl-G.Nghost ; i <= G.INu+G.Nghost; ++i ) {
         G.Node[i][j].X -= V;
     } /* endfor */
  } /* endfor */
  for ( j = G.JCl-G.Nghost ; j <= G.JCu+G.Nghost ; ++j) {
      for ( i = G.ICl-G.Nghost ; i <= G.ICu+G.Nghost ; ++i) {
  	  G.Cell[i][j].Xc = G.centroid(i, j);
          G.Cell[i][j].A = G.area(i, j);
      } /* endfor */
  } /* endfor */
  if (G.BndNorthSpline.np != 0 ) G.BndNorthSpline = G.BndNorthSpline - V;
  if (G.BndSouthSpline.np != 0 ) G.BndSouthSpline = G.BndSouthSpline - V;
  if (G.BndEastSpline.np != 0 ) G.BndEastSpline = G.BndEastSpline - V;
  if (G.BndWestSpline.np != 0 ) G.BndWestSpline = G.BndWestSpline - V;
  return (G);
}

// Scaling operators.
inline Grid2D_Quad_Block operator *(Grid2D_Quad_Block &G, const double &a) {
  int i, j;
  for ( j = G.JNl-G.Nghost ; j <= G.JNu+G.Nghost; ++j ) {
     for ( i = G.INl-G.Nghost ; i <= G.INu+G.Nghost; ++i ) {
         G.Node[i][j].X = G.Node[i][j].X*a;
     } /* endfor */
  } /* endfor */
  for ( j = G.JCl-G.Nghost ; j <= G.JCu+G.Nghost ; ++j) {
      for ( i = G.ICl-G.Nghost ; i <= G.ICu+G.Nghost ; ++i) {
  	  G.Cell[i][j].Xc = G.centroid(i, j);
          G.Cell[i][j].A = G.area(i, j);
      } /* endfor */
  } /* endfor */
  if (G.BndNorthSpline.np != 0 ) G.BndNorthSpline = G.BndNorthSpline*a;
  if (G.BndSouthSpline.np != 0 ) G.BndSouthSpline = G.BndSouthSpline*a;
  if (G.BndEastSpline.np != 0 ) G.BndEastSpline = G.BndEastSpline*a;
  if (G.BndWestSpline.np != 0 ) G.BndWestSpline = G.BndWestSpline*a;
  return (G);
}

inline Grid2D_Quad_Block operator *(const double &a, Grid2D_Quad_Block &G) {
  int i, j;
  for ( j = G.JNl-G.Nghost ; j <= G.JNu+G.Nghost; ++j ) {
     for ( i = G.INl-G.Nghost ; i <= G.INu+G.Nghost; ++i ) {
         G.Node[i][j].X = G.Node[i][j].X*a;
     } /* endfor */
  } /* endfor */
  for ( j = G.JCl-G.Nghost ; j <= G.JCu+G.Nghost ; ++j) {
      for ( i = G.ICl-G.Nghost ; i <= G.ICu+G.Nghost ; ++i) {
  	  G.Cell[i][j].Xc = G.centroid(i, j);
          G.Cell[i][j].A = G.area(i, j);
      } /* endfor */
  } /* endfor */
  if (G.BndNorthSpline.np != 0 ) G.BndNorthSpline = G.BndNorthSpline*a;
  if (G.BndSouthSpline.np != 0 ) G.BndSouthSpline = G.BndSouthSpline*a;
  if (G.BndEastSpline.np != 0 ) G.BndEastSpline = G.BndEastSpline*a;
  if (G.BndWestSpline.np != 0 ) G.BndWestSpline = G.BndWestSpline*a;
  G.SminN = G.SminN*a; G.SmaxN = G.SmaxN*a; 
  G.SminS = G.SminS*a; G.SmaxS = G.SmaxS*a; 
  G.SminE = G.SminE*a; G.SmaxE = G.SmaxE*a; 
  G.SminW = G.SminW*a; G.SmaxW = G.SmaxW*a;
  return (G);
}

inline Grid2D_Quad_Block operator /(Grid2D_Quad_Block &G, const double &a) {
  int i, j;
  for ( j = G.JNl-G.Nghost ; j <= G.JNu+G.Nghost; ++j ) {
     for ( i = G.INl-G.Nghost ; i <= G.INu+G.Nghost; ++i ) {
         G.Node[i][j].X = G.Node[i][j].X/a;
     } /* endfor */
  } /* endfor */
  for ( j = G.JCl-G.Nghost ; j <= G.JCu+G.Nghost ; ++j) {
      for ( i = G.ICl-G.Nghost ; i <= G.ICu+G.Nghost ; ++i) {
  	  G.Cell[i][j].Xc = G.centroid(i, j);
          G.Cell[i][j].A = G.area(i, j);
      } /* endfor */
  } /* endfor */
  if (G.BndNorthSpline.np != 0 ) G.BndNorthSpline = G.BndNorthSpline/a;
  if (G.BndSouthSpline.np != 0 ) G.BndSouthSpline = G.BndSouthSpline/a;
  if (G.BndEastSpline.np != 0 ) G.BndEastSpline = G.BndEastSpline/a;
  if (G.BndWestSpline.np != 0 ) G.BndWestSpline = G.BndWestSpline/a;
  G.SminN = G.SminN/a; G.SmaxN = G.SmaxN/a; 
  G.SminS = G.SminS/a; G.SmaxS = G.SmaxS/a; 
  G.SminE = G.SminE/a; G.SmaxE = G.SmaxE/a; 
  G.SminW = G.SminW/a; G.SmaxW = G.SmaxW/a;
  return (G);
}

// Rotation operator.
inline Grid2D_Quad_Block operator ^(Grid2D_Quad_Block &G, const double &a) {
  int i, j; double cos_angle, sin_angle; Vector2D X;
  cos_angle = cos(-a); sin_angle = sin(-a);
  for ( j = G.JNl-G.Nghost ; j <= G.JNu+G.Nghost; ++j ) {
     for ( i = G.INl-G.Nghost ; i <= G.INu+G.Nghost; ++i ) {
         X.x = G.Node[i][j].X.x*cos_angle +
               G.Node[i][j].X.y*sin_angle;
         X.y = - G.Node[i][j].X.x*sin_angle +
                 G.Node[i][j].X.y*cos_angle;
         G.Node[i][j].X = X;
     } /* endfor */
  } /* endfor */
  for ( j = G.JCl-G.Nghost ; j <= G.JCu+G.Nghost ; ++j) {
      for ( i = G.ICl-G.Nghost ; i <= G.ICu+G.Nghost ; ++i) {
  	  G.Cell[i][j].Xc = G.centroid(i, j);
          G.Cell[i][j].A = G.area(i, j);
      } /* endfor */
  } /* endfor */
  if (G.BndNorthSpline.np != 0 ) G.BndNorthSpline = G.BndNorthSpline^a;
  if (G.BndSouthSpline.np != 0 ) G.BndSouthSpline = G.BndSouthSpline^a;
  if (G.BndEastSpline.np != 0 ) G.BndEastSpline = G.BndEastSpline^a;
  if (G.BndWestSpline.np != 0 ) G.BndWestSpline = G.BndWestSpline^a;
  return (G);
}

/*************************************************************************
 * Grid2D_Quad_Block -- Input-output operators.                          *
 *************************************************************************/
inline ostream &operator << (ostream &out_file, 
                             const Grid2D_Quad_Block &G) {
  int i, j;
  out_file << G.NNi << " " << G.INl << " " << G.INu << "\n";
  out_file << G.NNj << " " << G.JNl << " " << G.JNu << "\n";
  out_file << G.Nghost << "\n";
  if (G.NNi == 0 || G.NNj == 0) return(out_file);
  out_file << G.NCi << " " << G.ICl << " " << G.ICu << "\n";
  out_file << G.NCj << " " << G.JCl << " " << G.JCu << "\n";
  for ( j = G.JNl-G.Nghost ; j <= G.JNu+G.Nghost; ++j ) {
     for ( i = G.INl-G.Nghost ; i <= G.INu+G.Nghost; ++i ) {
         out_file << G.Node[i][j].X << "\n";
     } /* endfor */
  } /* endfor */
  for ( i = G.ICl-G.Nghost ; i <= G.ICu+G.Nghost ; ++i) {
      out_file << G.BCtypeN[i] << " " << G.BCtypeS[i] << "\n";
  } /* endfor */
  for ( j = G.JCl-G.Nghost ; j <= G.JCu+G.Nghost ; ++j) {
      out_file << G.BCtypeE[j] << " " << G.BCtypeW[j] << "\n";
  } /* endfor */
  if (G.BndNorthSpline.np != 0 ) {
      out_file << G.BndNorthSpline.np   << " " 
               << G.BndNorthSpline.type << "\n"; 
      out_file << G.BndNorthSpline;
  } else {
     out_file << G.BndNorthSpline.np << "\n";
  } /* endif */
  if (G.BndSouthSpline.np != 0 ) {
      out_file << G.BndSouthSpline.np   << " " 
               << G.BndSouthSpline.type << "\n"; 
      out_file << G.BndSouthSpline;
  } else {
     out_file << G.BndSouthSpline.np << "\n";
  } /* endif */
  if (G.BndEastSpline.np != 0 ) {
      out_file << G.BndEastSpline.np   << " " 
               << G.BndEastSpline.type << "\n"; 
      out_file << G.BndEastSpline;
  } else {
     out_file << G.BndEastSpline.np << "\n";
  } /* endif */
  if (G.BndWestSpline.np != 0 ) {
      out_file << G.BndWestSpline.np   << " " 
               << G.BndWestSpline.type << "\n"; 
      out_file << G.BndWestSpline;
  } else {
     out_file << G.BndWestSpline.np << "\n";
  } /* endif */
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

inline istream &operator >> (istream &in_file, 
                             Grid2D_Quad_Block &G) {
  int i, j, ni, il, iu, nj, jl, ju, ng;
  in_file.setf(ios::skipws);
  in_file >> ni >> il >> iu;
  in_file >> nj >> jl >> ju;
  in_file >> ng;
  in_file.unsetf(ios::skipws);
  if (ni == 0 || nj == 0) {
      if (G.Node != NULL) G.deallocate(); return(in_file);
  } /* endif */
  if (G.Node == NULL || G.Cell == NULL || G.NNi != ni || G.NNj != nj) {
      if (G.Node != NULL) G.deallocate(); 
      G.allocate(ni-2*ng-1, nj-2*ng-1, ng);
  } /* endif */
  in_file.setf(ios::skipws);
  in_file >> ni >> il >> iu; in_file >> nj >> jl >> ju;
  in_file.unsetf(ios::skipws);
  for ( j = G.JNl-G.Nghost ; j <= G.JNu+G.Nghost; ++j ) {
      for ( i = G.INl-G.Nghost ; i <= G.INu+G.Nghost; ++i ) {
          in_file >> G.Node[i][j].X;
      } /* endfor */
  } /* endfor */
  for ( j = G.JCl-G.Nghost ; j <= G.JCu+G.Nghost ; ++j) {
      for ( i = G.ICl-G.Nghost ; i <= G.ICu+G.Nghost ; ++i) {
  	  G.Cell[i][j].I = i; G.Cell[i][j].J = j;
  	  G.Cell[i][j].Xc = G.centroid(i, j);
          G.Cell[i][j].A = G.area(i, j);
      } /* endfor */
  } /* endfor */
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
  in_file.setf(ios::skipws);
  in_file >> ni;
  in_file.unsetf(ios::skipws);
  if (G.BndNorthSpline.np != 0) G.BndNorthSpline.deallocate(); 
  if (ni != 0 ) {
      G.BndNorthSpline.allocate(ni);
      in_file.setf(ios::skipws);
      in_file >> i; G.BndNorthSpline.settype(i);
      in_file.unsetf(ios::skipws);
      in_file >> G.BndNorthSpline; G.BndNorthSpline.pathlength();
  } /* endif */
  in_file.setf(ios::skipws);
  in_file >> ni;
  in_file.unsetf(ios::skipws);
  if (G.BndSouthSpline.np != 0) G.BndSouthSpline.deallocate(); 
  if (ni != 0 ) {
      G.BndSouthSpline.allocate(ni);
      in_file.setf(ios::skipws);
      in_file >> i; G.BndSouthSpline.settype(i);
      in_file.unsetf(ios::skipws);
      in_file >> G.BndSouthSpline; G.BndSouthSpline.pathlength();
  } /* endif */
  in_file.setf(ios::skipws);
  in_file >> nj;
  in_file.unsetf(ios::skipws);
  if (G.BndEastSpline.np != 0) G.BndEastSpline.deallocate(); 
  if (nj != 0 ) {
      G.BndEastSpline.allocate(nj);
      in_file.setf(ios::skipws);
      in_file >> j; G.BndEastSpline.settype(j);
      in_file.unsetf(ios::skipws);
      in_file >> G.BndEastSpline; G.BndEastSpline.pathlength();
  } /* endif */
  in_file.setf(ios::skipws);
  in_file >> nj;
  in_file.unsetf(ios::skipws);
  if (G.BndWestSpline.np != 0) G.BndWestSpline.deallocate(); 
  if (nj != 0 ) {
      G.BndWestSpline.allocate(nj);
      in_file.setf(ios::skipws);
      in_file >> j; G.BndWestSpline.settype(j);
      in_file.unsetf(ios::skipws);
      in_file >> G.BndWestSpline; G.BndWestSpline.pathlength();
  } /* endif */
  in_file.setf(ios::skipws);
  in_file >> G.SminN >> G.SmaxN >> G.SminS >> G.SmaxS; 
  in_file >> G.SminE >> G.SmaxE >> G.SminW >> G.SmaxW;
  in_file >> G.StretchI >> G.BetaI >> G.TauI;
  in_file >> G.StretchJ >> G.BetaJ >> G.TauJ;
  in_file >> G.OrthogonalN >> G.OrthogonalS 
          >> G.OrthogonalE >> G.OrthogonalW;
  in_file.unsetf(ios::skipws);
  return (in_file);
}

/*************************************************************************
 * Grid2D_Quad_Block -- External subroutines for single grid block.      *
 *************************************************************************/

extern void Create_Quad_Block(Grid2D_Quad_Block &Grid,
		              const int Number_of_Cells_Idir,
		              const int Number_of_Cells_Jdir,
			      const int Number_of_Ghost_Cells);

extern void Create_Quad_Block(Grid2D_Quad_Block &Grid,
                              char *Bnd_Spline_File_Name_ptr,
		              const int Number_of_Cells_Idir,
		              const int Number_of_Cells_Jdir,
			      const int Number_of_Ghost_Cells);

extern void Create_Quad_Block(Grid2D_Quad_Block &Grid,
                              Spline2D &Bnd_Spline_North,
                              Spline2D &Bnd_Spline_South,
                              Spline2D &Bnd_Spline_East,
                              Spline2D &Bnd_Spline_West,
                              const int Number_of_Cells_Idir,
  	                      const int Number_of_Cells_Jdir,
			      const int Number_of_Ghost_Cells,
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

extern void Broadcast_Quad_Block(Grid2D_Quad_Block &Grid);

#ifdef _MPI_VERSION
extern void Broadcast_Quad_Block(Grid2D_Quad_Block &Grid,
                                 MPI::Intracomm &Communicator, 
                                 const int Source_CPU);
#endif

extern void Smooth_Quad_Block(Grid2D_Quad_Block &Grid,
		              const int Number_of_Iterations);

extern void Smooth_Rocket_Motor(Grid2D_Quad_Block &Grid,
				const double &Length_Chamber,
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

extern void Set_BCs(Grid2D_Quad_Block &Grid);

extern void Update_Exterior_Nodes(Grid2D_Quad_Block &Grid);

extern void Update_Corner_Ghost_Nodes(Grid2D_Quad_Block &Grid);

extern void Update_Cells(Grid2D_Quad_Block &Grid);

extern int Check_Quad_Block(Grid2D_Quad_Block &Grid);

extern void Write_Quad_Block_Definition(Grid2D_Quad_Block &Grid,
                                        ostream &Out_File);

extern void Read_Quad_Block_Definition(Grid2D_Quad_Block &Grid,
                                       istream &In_File);

extern void Write_Quad_Block(Grid2D_Quad_Block &Grid,
	                     ostream &Out_File);

extern void Read_Quad_Block(Grid2D_Quad_Block &Grid,
	                    istream &In_File);

extern void Copy_Quad_Block(Grid2D_Quad_Block &Grid1,
		            Grid2D_Quad_Block &Grid2);

extern void Translate_Quad_Block(Grid2D_Quad_Block &Grid,
	      	                 const Vector2D &V);

extern void Scale_Quad_Block(Grid2D_Quad_Block &Grid,
	      	             const double &Scaling_Factor);

extern void Rotate_Quad_Block(Grid2D_Quad_Block &Grid,
	      	              const double &Angle);

extern void Reflect_Quad_Block(Grid2D_Quad_Block &Grid);

extern void Output_Tecplot(Grid2D_Quad_Block &Grid,
                           const int Block_Number,
                           const int Output_Title,
 	                   ostream &Out_File);

extern void Output_Nodes_Tecplot(Grid2D_Quad_Block &Grid,
                                 const int Block_Number,
                                 const int Output_Title,
 	                         ostream &Out_File);

extern void Output_Cells_Tecplot(Grid2D_Quad_Block &Grid,
                                 const int Block_Number,
                                 const int Output_Title,
 	                         ostream &Out_File);

extern void Output_Gnuplot(Grid2D_Quad_Block &Grid,
                           const int Block_Number,
                           const int Output_Title,
 	                   ostream &Out_File);

extern void Double_Mesh_Resolution(Grid2D_Quad_Block &Grid_Double,
                                   Grid2D_Quad_Block &Grid_Original);

extern void Half_Mesh_Resolution(Grid2D_Quad_Block &Grid_Half,
                                 Grid2D_Quad_Block &Grid_Original);

extern void Refine_Mesh(Grid2D_Quad_Block &Grid_Fine,
                        Grid2D_Quad_Block &Grid_Original,
                        const int Sector);

extern void Coarsen_Mesh(Grid2D_Quad_Block &Grid_Coarse,
                         Grid2D_Quad_Block &Grid_Original_SW,
                         Grid2D_Quad_Block &Grid_Original_SE,
                         Grid2D_Quad_Block &Grid_Original_NW,
                         Grid2D_Quad_Block &Grid_Original_NE);

extern void Fix_Refined_Mesh_Boundaries(Grid2D_Quad_Block &Grid,
                                        const int Fix_North_Boundary,
                                        const int Fix_South_Boundary,
                                        const int Fix_East_Boundary,
                                        const int Fix_West_Boundary);

extern void Unfix_Refined_Mesh_Boundaries(Grid2D_Quad_Block &Grid);

/*************************************************************************
 * Grid2D_Quad_Block -- External subroutines for 2D array of grid blocks.*
 *************************************************************************/

extern Grid2D_Quad_Block** Allocate_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
						     const int Number_of_Blocks_Idir,
						     const int Number_of_Blocks_Jdir);

extern Grid2D_Quad_Block** Deallocate_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
						       const int Number_of_Blocks_Idir,
						       const int Number_of_Blocks_Jdir);

extern Grid2D_Quad_Block** Broadcast_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
				                      int &Number_of_Blocks_Idir,
		                                      int &Number_of_Blocks_Jdir);

extern void Write_Multi_Block_Grid_Definition(Grid2D_Quad_Block **Grid_ptr,
                 			      const int Number_of_Blocks_Idir,
		                              const int Number_of_Blocks_Jdir,
                                              ostream &Out_File);

extern Grid2D_Quad_Block** Read_Multi_Block_Grid_Definition(Grid2D_Quad_Block **Grid_ptr,
				                            int &Number_of_Blocks_Idir,
		                                            int &Number_of_Blocks_Jdir,
                                                            istream &In_File);

extern void Write_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
				   const int Number_of_Blocks_Idir,
		                   const int Number_of_Blocks_Jdir,
                                   ostream &Out_File);

extern Grid2D_Quad_Block** Read_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
				                 int &Number_of_Blocks_Idir,
		                                 int &Number_of_Blocks_Jdir,
                                                 istream &In_File);

extern void Translate_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                       const int Number_of_Blocks_Idir,
                                       const int Number_of_Blocks_Jdir,
	      	                       const Vector2D &V);

extern void Scale_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                   const int Number_of_Blocks_Idir,
                                   const int Number_of_Blocks_Jdir,
	      	                   const double &Scaling_Factor);

extern void Rotate_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                    const int Number_of_Blocks_Idir,
                                    const int Number_of_Blocks_Jdir,
	      	                    const double &Angle);

extern void Reflect_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                     const int Number_of_Blocks_Idir,
                                     const int Number_of_Blocks_Jdir);

extern int Check_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                  const int Number_of_Blocks_Idir,
                                  const int Number_of_Blocks_Jdir);

extern void Output_Tecplot(Grid2D_Quad_Block **Grid_ptr,
			   const int Number_of_Blocks_Idir,
		           const int Number_of_Blocks_Jdir,
 	                   ostream &Out_File);

extern void Output_Nodes_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                                 const int Number_of_Blocks_Idir,
                                 const int Number_of_Blocks_Jdir,
 	                         ostream &Out_File);

extern void Output_Cells_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                                 const int Number_of_Blocks_Idir,
                                 const int Number_of_Blocks_Jdir,
 	                         ostream &Out_File);

extern void Output_Gnuplot(Grid2D_Quad_Block **Grid_ptr,
			   const int Number_of_Blocks_Idir,
		           const int Number_of_Blocks_Jdir,
 	                   ostream &Out_File);

extern Grid2D_Quad_Block** Grid_Rectangular_Box(Grid2D_Quad_Block **Grid_ptr,
                                                int &Number_of_Blocks_Idir,
		                                int &Number_of_Blocks_Jdir,
                                                const double &Width,
                                                const double &Height,
 		                                const int Number_of_Cells_Idir,
		                                const int Number_of_Cells_Jdir,
						const int Number_of_Ghost_Cells);

Grid2D_Quad_Block** Grid_Rectangular_Box(Grid2D_Quad_Block **Grid_ptr,
                                         int &Number_of_Blocks_Idir,
                                         int &Number_of_Blocks_Jdir,
                                         const double &Width,
                                         const double &Height,
					 const int &Stretching_Flag,
					 const int &Stretching_Type_Idir,
					 const int &Stretching_Type_Jdir,
					 const double &Stretching_Factor_Idir,
					 const double &Stretching_Factor_Jdir,
 	                                 const int Number_of_Cells_Idir,
	                                 const int Number_of_Cells_Jdir,
					 const int Number_of_Ghost_Cells);

extern Grid2D_Quad_Block** Grid_Flat_Plate(Grid2D_Quad_Block **Grid_ptr,
                                           int &Number_of_Blocks_Idir,
                                           int &Number_of_Blocks_Jdir,
                                           const double &Length,
					   const int &Flat_Plate_BC_Type,
					   const int &Stretching_Flag,
					   const double &Stretching_Factor_Idir,
					   const double &Stretching_Factor_Jdir,
 		                           const int Number_of_Cells_Idir,
		                           const int Number_of_Cells_Jdir,
					   const int Number_of_Ghost_Cells);

extern Grid2D_Quad_Block** Grid_1D_Flame(Grid2D_Quad_Block **Grid_ptr,
					 int &Number_of_Blocks_Idir,
					 int &Number_of_Blocks_Jdir,
					 const double &Length,
					 const double &Heigth,
					 const int Number_of_Cells_Idir,
					 const int Number_of_Cells_Jdir,
					 const int Number_of_Ghost_Cells);

extern Grid2D_Quad_Block** Grid_2D_Laminar_Flame(Grid2D_Quad_Block **Grid_ptr,
						 int &Number_of_Blocks_Idir,
						 int &Number_of_Blocks_Jdir,
						 const double &Length,
						 const double &Heigth,
						 const int Number_of_Cells_Idir,
						 const int Number_of_Cells_Jdir,
						 const int Number_of_Ghost_Cells); 

extern Grid2D_Quad_Block** Grid_Pipe(Grid2D_Quad_Block **Grid_ptr,
                                     int &Number_of_Blocks_Idir,
                                     int &Number_of_Blocks_Jdir,
                                     const double &Length,
				     const double &Radius,
				     const int Stretching_Flag,
				     const double Stretching_Factor,
 		                     const int Number_of_Cells_Idir,
		                     const int Number_of_Cells_Jdir,
				     const int Number_of_Ghost_Cells);

extern Grid2D_Quad_Block** Grid_Pipe(Grid2D_Quad_Block **Grid_ptr,
                                     int &Number_of_Blocks_Idir,
                                     int &Number_of_Blocks_Jdir,
                                     const double &Length,
				     const double &Radius,
				     const int &Axisymmetric,
 		                     const int Number_of_Cells_Idir,
		                     const int Number_of_Cells_Jdir,
				     const int Number_of_Ghost_Cells);

extern Grid2D_Quad_Block** Grid_Blunt_Body(Grid2D_Quad_Block **Grid_ptr,
                                           int &Number_of_Blocks_Idir,
                                           int &Number_of_Blocks_Jdir,
                                           const double &Radius,
				           const double &Mach_Number,
 		                           const int Number_of_Cells_Idir,
		                           const int Number_of_Cells_Jdir,
					   const int Number_of_Ghost_Cells);

extern Grid2D_Quad_Block** Grid_Rocket_Motor(Grid2D_Quad_Block **Grid_ptr,
					     int &Number_of_Blocks_Idir,
					     int &Number_of_Blocks_Jdir,
					     const double &Length_Chamber,
					     const double &Radius_Chamber,
					     const double &Length_Chamber_To_Throat,
					     const double &Length_Nozzle,
					     const double &Radius_Nozzle_Exit,
					     const double &Radius_Nozzle_Throat,
					     const double &Radius_Grain,
					     const int &Nozzle_Type,
					     const int &Chamber_BC_Type,
					     const int &Stretching_Flag,
					     const int Stretching_Type_Jdir,
					     const double &Stretching_Factor_Idir,
					     const double &Stretching_Factor_Jdir,
					     const int Number_of_Cells_Idir,
					     const int Number_of_Cells_Jdir,
					     const int Number_of_Ghost_Cells);

extern Grid2D_Quad_Block** Grid_Nozzleless_Rocket_Motor(Grid2D_Quad_Block **Grid_ptr,
							int &Number_of_Blocks_Idir,
							int &Number_of_Blocks_Jdir,
							const double &Length_Chamber,
							const double &Radius_Chamber,
							const double &Length_Nozzle,
							const double &Radius_Nozzle_Exit,
							const int &Chamber_BC_Type,
							const int &Stretching_Flag,
							const int Stretching_Type_Jdir,
							const double &Stretching_Factor_Idir,
							const double &Stretching_Factor_Jdir,
							const int Number_of_Cells_Idir,
							const int Number_of_Cells_Jdir,
							const int Number_of_Ghost_Cells);


extern Grid2D_Quad_Block** Grid_Nozzle(Grid2D_Quad_Block **Grid_ptr,
				       int &Number_of_Blocks_Idir,
				       int &Number_of_Blocks_Jdir,
				       const double &Length_Nozzle,
				       const double &Radius_Chamber,
				       const double &Radius_Nozzle_Exit,
				       const double &Radius_Nozzle_Throat,
				       const int &Nozzle_Type,
				       const int &Stretching_Flag,
				       const int &Stretching_Type_Idir,
				       const int &Stretching_Type_Jdir,
				       const double &Stretching_Factor_Idir,
				       const double &Stretching_Factor_Jdir,
				       const int Number_of_Cells_Idir,
				       const int Number_of_Cells_Jdir,
				       const int Number_of_Ghost_Cells);

extern Grid2D_Quad_Block** Grid_Circular_Cylinder(Grid2D_Quad_Block **Grid_ptr,
                                                  int &Number_of_Blocks_Idir,
                                                  int &Number_of_Blocks_Jdir,
                                                  const double &Radius,
						  const int &Stretching_Type_Idir,
						  const int &Stretching_Type_Jdir,
						  const double &Stretching_Factor_Idir,
						  const double &Stretching_Factor_Jdir,
 		                                  const int Number_of_Cells_Idir,
		                                  const int Number_of_Cells_Jdir,
						  const int Number_of_Ghost_Cells);

extern Grid2D_Quad_Block** Grid_Circular_Cylinder(Grid2D_Quad_Block **Grid_ptr,
                                                  int &Number_of_Blocks_Idir,
                                                  int &Number_of_Blocks_Jdir,
                                                  const double &Inner_Radius,
                                                  const double &Outer_Radius,
						  const int &Stretching_Type_Idir,
						  const int &Stretching_Type_Jdir,
						  const double &Stretching_Factor_Idir,
						  const double &Stretching_Factor_Jdir,
 		                                  const int Number_of_Cells_Idir,
		                                  const int Number_of_Cells_Jdir,
						  const int Number_of_Ghost_Cells);

extern Grid2D_Quad_Block** Grid_Ellipse(Grid2D_Quad_Block **Grid_ptr,
                                        int &Number_of_Blocks_Idir,
                                        int &Number_of_Blocks_Jdir,
                                        const double &A,
					const double &B,
 		                        const int Number_of_Cells_Idir,
		                        const int Number_of_Cells_Jdir,
					const int Number_of_Ghost_Cells);

extern Grid2D_Quad_Block** Grid_NACA_Aerofoil(Grid2D_Quad_Block **Grid_ptr,
                                              int &Number_of_Blocks_Idir,
                                              int &Number_of_Blocks_Jdir,
                                              char *NACA_Aerofoil_Type_ptr,
                                              const double &Chord_Length,
 		                              const int Number_of_Cells_Idir,
		                              const int Number_of_Cells_Jdir,
					      const int Number_of_Ghost_Cells);

extern Grid2D_Quad_Block** Grid_Free_Jet(Grid2D_Quad_Block **Grid_ptr,
                                         int &Number_of_Blocks_Idir,
                                         int &Number_of_Blocks_Jdir,
				         const double &Radius,
 		                         const int Number_of_Cells_Idir,
		                         const int Number_of_Cells_Jdir,
					 const int Number_of_Ghost_Cells);

extern Grid2D_Quad_Block** Grid_Wedge(Grid2D_Quad_Block **Grid_ptr,
                                      int &Number_of_Blocks_Idir,
                                      int &Number_of_Blocks_Jdir,
				      const double &Wedge_Angle,
				      const double &Wedge_Length,
				      const int &Wedge_BC_Type,
				      const int &Stretching_Flag,
				      const double &Stretching_Factor_Idir,
				      const double &Stretching_Factor_Jdir,
				      const int Number_of_Cells_Idir,
				      const int Number_of_Cells_Jdir,
				      const int Number_of_Ghost_Cells);

extern Grid2D_Quad_Block** Grid_Unsteady_Blunt_Body(Grid2D_Quad_Block **Grid_ptr,
                                                    int &Number_of_Blocks_Idir,
                                                    int &Number_of_Blocks_Jdir,
						    const double &Radius,
						    const double &Mach_Number,
						    const int Number_of_Cells_Idir,
						    const int Number_of_Cells_Jdir,
						    const int Number_of_Ghost_Cells);

extern Grid2D_Quad_Block** Grid_Ringleb_Flow(Grid2D_Quad_Block **Grid_ptr,
                                             int &Number_of_Blocks_Idir,
                                             int &Number_of_Blocks_Jdir,
					     const double &Inner_Streamline_Number,
					     const double &Outer_Streamline_Number,
					     const double &Isotach_Line,
				             const int Number_of_Cells_Idir,
				             const int Number_of_Cells_Jdir,
					     const int Number_of_Ghost_Cells);

extern Grid2D_Quad_Block** Grid_Bump_Channel_Flow(Grid2D_Quad_Block **Grid_ptr,
						  int &Number_of_Blocks_Idir,
						  int &Number_of_Blocks_Jdir,
						  const double &Radius,
						  const int Smooth_Bump,
						  const int Number_of_Cells_Idir,
						  const int Number_of_Cells_Jdir,
						  const int Number_of_Ghost_Cells);

extern Grid2D_Quad_Block** Grid_Driven_Cavity_Flow(Grid2D_Quad_Block **Grid_ptr,
						   int &Number_of_Blocks_Idir,
						   int &Number_of_Blocks_Jdir,
						   const double &Width,
						   const double &Height,
						   const int &Stretching_Type_Idir,
						   const int &Stretching_Type_Jdir,
						   const double &Stretching_Factor_Idir,
						   const double &Stretching_Factor_Jdir,
						   const int Number_of_Cells_Idir,
						   const int Number_of_Cells_Jdir,
						   const int Number_of_Ghost_Cells);

extern Grid2D_Quad_Block** Grid_Backward_Facing_Step(Grid2D_Quad_Block **Grid_ptr,
						     int &Number_of_Blocks_Idir,
						     int &Number_of_Blocks_Jdir,
						     const double &Step_Height,
						     const double &Top_Wall_Deflection,
						     const double Stretching_Factor_Idir,
						     const double Stretching_Factor_Jdir,
						     const int Number_of_Cells_Idir,
						     const int Number_of_Cells_Jdir,
						     const int Number_of_Ghost_Cells);

extern Grid2D_Quad_Block** Grid_Mixing_Layer(Grid2D_Quad_Block **Grid_ptr,
					     int &Number_of_Blocks_Idir,
					     int &Number_of_Blocks_Jdir,
					     const double &Length,
					     const double &Stretching_Factor_Idir,
					     const double &Stretching_Factor_Jdir,
					     const int Number_of_Cells_Idir,
					     const int Number_of_Cells_Jdir,
					     const int Number_of_Ghost_Cells);

extern Grid2D_Quad_Block** Grid_NASA_Rotor_37(Grid2D_Quad_Block **Grid_ptr,
					      int &Number_of_Blocks_Idir,
					      int &Number_of_Blocks_Jdir,
					      const double &Rotor_Percent_Span,
					      const int Number_of_Cells_Idir,
					      const int Number_of_Cells_Jdir,
					      const int Number_of_Ghost_Cells);

extern Grid2D_Quad_Block** Grid_NASA_Rotor_67(Grid2D_Quad_Block **Grid_ptr,
					      int &Number_of_Blocks_Idir,
					      int &Number_of_Blocks_Jdir,
					      const double &Rotor_Percent_Span,
					      const int Number_of_Cells_Idir,
					      const int Number_of_Cells_Jdir,
					      const int Number_of_Ghost_Cells);

extern Grid2D_Quad_Block** Grid_Desolvation_Chamber(Grid2D_Quad_Block **Grid_ptr,
						    const int &Chamber_BC_Type,
						    int &Number_of_Blocks_Idir,
						    int &Number_of_Blocks_Jdir,
						    const int Number_of_Cells_Idir,
						    const int Number_of_Cells_Jdir,
						    const int Number_of_Ghost_Cells);

#endif /* _GRID2D_QUAD_BLOCK_INCLUDED  */

/* Grid3DHexa.h:  Header file defining 3D hexahedral block grid type. */
#ifndef _GRID3D_HEXA_BLOCK_INCLUDED
#define _GRID3D_HEXA_BLOCK_INCLUDED

/* Include required C++ libraries. */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>

using namespace std;

/* Include math macro, CFD, 3D vector, 3D cell*/

#ifndef _MATH_MACROS_INCLUDED
#include "../Math/Math.h"
#endif // _MATH_MACROS_INCLUDED

#ifndef _CFD_INCLUDED
#include "../CFD/CFD.h"
#endif // _CFD_INCLUDED

#ifndef _VECTOR3D_INCLUDED
#include "../Math/Vector3D.h"
#endif //_VECTOR3D_INCLUDED

#ifndef _CELL3D_INCLUDED
#include "Cell3D.h"
#endif // _CELL3D_INCLUDED

/* Define the HEXA 3D grid block class. */

/********************************************************
 * Class: Grid3D_Hexa_Block                             *
 *                                                      *
 * Member functions                                     *
 *      NNi        -- Return number of nodes in         *
 *                    the i-direction (zeta-direction). *
 *      INl        -- Return lower index for nodes in   *
 *                    the i-direction (zeta-direction). *
 *      INu        -- Return upper index for nodes in   *
 *                    the i-direction (zeta-direction). *
 *      NNj        -- Return number of nodes in         *
 *                    the j-direction (eta-direction).  *
 *      JNl        -- Return lower index for nodes in   *
 *                    the j-direction (eta-direction).  *
 *      JNu        -- Return upper index for nodes in   *
 *                    the j-direction (eta-direction).  *
 *      NNk        -- Return number of nodes in         *
 *                    the k-direction                   *
 *      KNl        -- Return lower index for nodes in   *
 *                    the k-direction (zeta-direction). *
 *      KNu        -- Return upper index for nodes in   *
 *                    the k-direction (zeta-direction). * 
 *      NCi        -- Return number of cells in         *
 *                    the i-direction (zeta-direction). *
 *      ICl        -- Return lower index for cells in   *
 *                    the i-direction (zeta-direction). *
 *      ICu        -- Return upper index for cells in   *
 *                    the i-direction (zeta-direction). *
 *      NCj        -- Return number of cells in         *
 *                    the j-direction (eta-direction).  *
 *      JCl        -- Return lower index for cells in   *
 *                    the j-direction (eta-direction).  *
 *      JCu        -- Return upper index for cells in   *
 *                    the j-direction (eta-direction).  *
 *      NCk        -- Return number of cells in         *
 *                    the k-direction                   *
 *      KCl        -- Return lower index for cells in   *
 *                    the k-direction                   *
 *      KCu        -- Return upper index for cells in   *
 *                    the k-direction                   *
 *      Node       -- Return 3D node geometry.          *
 *      Cell       -- Return 3D cell geometry.          *
 *      allocate   -- Allocate memory for structured    *
 *                    quadrilateral grid block.         *
 *      deallocate -- Deallocate memory for structured  *
 *                    quadrilateral grid block.         *
 *      allocateNodes -- Allocate memory for nodes of   *
 *                       structured quadrilateral grid  *
 *                       block.                         *
 *      deallocateNodes -- Deallocate memory for nodes  *
 *                         of structured quadrilateral  *
 *                         grid block.                  *
 *      allocateCells -- Allocate memory for nodes of   *
 *                       structured quadrilateral grid  *
 *                       block.                         *
 *      deallocateCells -- Deallocate memory for nodes  *
 *                         of structured quadrilateral  *
 *                         grid block.                  *
 *      centroid   -- Calculate 3D vector containing    *
 *                    the location of cell center.      *
 *      Volume     -- Calculate the volume for cell.    *
 *      nodeNWBot  -- Return NWBot node for cell.       *
 *      nodeNEBot  -- Return NEBot node for cell.       *
 *      nodeSEBot  -- Return SEBot node for cell.       *
 *      nodeSWBot  -- Return SWBot node for cell.       *
 *      nodeNWTop  -- Return NWTop node for cell.       *
 *      nodeNETop  -- Return NETop node for cell.       *
 *      nodeSETop  -- Return SETop node for cell.       *
 *      nodeSWTop  -- Return SWTop node for cell.       *
 *      xfaceN     -- Return midpoint of face.          *
 *      xfaceS     -- Return midpoint of face.          *
 *      xfaceE     -- Return midpoint of face.          *
 *      xfaceW     -- Return midpoint of face.          *
 *      xfaceTop   -- Return midpoint of face.          *
 *      xfaceBot   -- Return midpoint of face.          *
 *      AfaceN     -- Return area of N face.            *
 *      AfaceS     -- Return area of S face.            *
 *      AfaceE     -- Return area of E face.            *
 *      AfaceW     -- Return area of W face.            *
 *      AfaceTop   -- Return area of Top face.          *
 *      AfaceBot   -- Return area of Bot face.          *
 *      nfaceN     -- Return normal of N face.          *
 *      nfaceS     -- Return normal of S face.          *
 *      nfaceE     -- Return normal of E face.          *
 *      nfaceW     -- Return normal of W face.          *
 *      nfaceTop   -- Return normal of Top face.        *
 *      nfaceBot   -- Return normal of Bot face.        *
 * Member operators                                     *
 *  G  -- grid block consisting of hexa 3D cells        *
 *  V  -- a 3D vector                                   *
 *  a  -- a scalar (double)                             *
 *                                                      *
 * G = G;                                               *
 * G = G + V; (shift location of grid)                  *
 * G = G - V; (shift location of grid)                  *
 * G = a * G; (scale grid)                              *
 * G = G * a; (scale grid)                              *
 * G = G / a; (scale grid)                              *
 * G = G ^ a; (rotate grid about Z-axis)                *  
 * G = G & a; (rotate grid about Y-axis)                *
 * G = G % a; (rotate grid about X-axis)                *
 * cout << G; (output function)                         *
 * cin  >> G; (input function)                          *
 *                                                      *
 ********************************************************/
class Grid3D_Hexa_Block{
  private:
  public:
    int           NNi,INl,INu; // i-direction node counters
    int           NNj,JNl,JNu; // j-direction node counters
    int           NNk,KNl,KNu; // k-direction node counters
    int           NCi,ICl,ICu; // i-direction cell counters
    int           NCj,JCl,JCu; // j-direction cell counters
    int           NCk,KCl,KCu; // k-direction cell counters
    int                Nghost; //number of ghost cells
    Node3D            ***Node; // array of 3D node position vectors
    Cell3D            ***Cell; // array of 3D cell centre vectors
    
    /* Creation, copy, and assignment constructors. */
    Grid3D_Hexa_Block(void) {
      NNi = 0; INl = 0; INu = 0; NNj = 0; JNl = 0; JNu = 0; NNk = 0; KNl = 0; KNu = 0; 
      NCi = 0; ICl = 0; ICu = 0; NCj = 0; JCl = 0; JCu = 0; NCk = 0; KCl = 0; KCu = 0;
      Nghost = 0;
      Node = NULL; Cell = NULL; 
    }

    Grid3D_Hexa_Block(const Grid3D_Hexa_Block &G) {
       NNi = G.NNi; INl = G.INl; INu = G.INu; 
       NNj = G.NNj; JNl = G.JNl; JNu = G.JNu;
       NNk = G.NNk; KNl = G.KNl; KNu = G.KNu;
       NCi = G.NCi; ICl = G.ICl; ICu = G.ICu; 
       NCj = G.NCj; JCl = G.JCl; JCu = G.JCu;
       NCk = G.NCk; KCl = G.KCl; KCu = G.KCu;
       Nghost = G.Nghost;
       Node = G.Node; Cell = G.Cell;
    }

    /* Destructor. */
    // ~Grid2D_Quad_Block(void);
    // Use automatically generated destructor.

    /* Assignment operator. */
    // Grid2D_Quad_Block operator = (const Grid2D_Quad_Block &Soln);
    // Use automatically generated assignment operator.

    /* Allocate memory for structured quadrilateral grid block. */
    void allocate(const int Ni, const int Nj, const int Nk, const int Ng);

    /* Deallocate memory for structured quadrilateral grid block. */
    void deallocate(void);

    /* Allocate memory for nodes of structured quadrilateral grid block. */
    void allocateNodes(const int Ni, const int Nj, const int Nk, const int Ng);

    /* Deallocate memory for nodes of structured quadrilateral grid block. */
    void deallocateNodes(void);

    /* Allocate memory for cells of structured quadrilateral grid block. */
    void allocateCells(const int Ni, const int Nj, const int Nk, const int Ng);

    /* Deallocate memory for cells of structured quadrilateral grid block. */
    void deallocateCells(void);

    /* Calculate centroid of cell. */
    Vector3D centroid(const Cell3D &Cell);
    Vector3D centroid(const int ii, const int jj, const int kk);

    /* Calculate cell volume. */
    double volume(const Cell3D &Cell);
    double volume(const int ii, const int jj, const int kk);

    
    /* Get cell nodes. For the 3D case there are 8 nodes that need to be found*/    
    Node3D nodeNWBot(const Cell3D &Cell);
    Node3D nodeNWBot(const int ii, const int jj, const int kk);
    
    Node3D nodeNEBot(const Cell3D &Cell);
    Node3D nodeNEBot(const int ii, const int jj, const int kk);
    
    Node3D nodeSEBot(const Cell3D &Cell);
    Node3D nodeSEBot(const int ii, const int jj, const int kk);
    
    Node3D nodeSWBot(const Cell3D &Cell);
    Node3D nodeSWBot(const int ii, const int jj, const int kk);
    
    Node3D nodeNWTop(const Cell3D &Cell);
    Node3D nodeNWTop(const int ii, const int jj, const int kk);

    Node3D nodeNETop(const Cell3D &Cell);
    Node3D nodeNETop(const int ii, const int jj, const int kk);
    
    Node3D nodeSETop(const Cell3D &Cell);
    Node3D nodeSETop(const int ii, const int jj, const int kk);

    Node3D nodeSWTop(const Cell3D &Cell);
    Node3D nodeSWTop(const int ii, const int jj, const int kk);

    /*********************Calculate midpoints of faces***********************/
    Vector3D xfaceN(const Cell3D &Cell);
    Vector3D xfaceN(const int ii, const int jj, const int kk); 

    Vector3D xfaceS(const Cell3D &Cell);
    Vector3D xfaceS(const int ii, const int jj, const int kk); 

    Vector3D xfaceE(const Cell3D &Cell);
    Vector3D xfaceE(const int ii, const int jj, const int kk); 

    Vector3D xfaceW(const Cell3D &Cell);
    Vector3D xfaceW(const int ii, const int jj, const int kk); 

    Vector3D xfaceTop(const Cell3D &Cell);
    Vector3D xfaceTop(const int ii, const int jj, const int kk); 

    Vector3D xfaceBot(const Cell3D &Cell);
    Vector3D xfaceBot(const int ii, const int jj, const int kk); 

    /*******************Calculate normal to each face***********************/
    Vector3D nfaceN(const Cell3D &Cell);
    Vector3D nfaceN(const int ii, const int jj, const int kk);

    Vector3D nfaceS(const Cell3D &Cell);
    Vector3D nfaceS(const int ii, const int jj, const int kk);

    Vector3D nfaceE(const Cell3D &Cell);
    Vector3D nfaceE(const int ii, const int jj, const int kk);

    Vector3D nfaceW(const Cell3D &Cell);
    Vector3D nfaceW(const int ii, const int jj, const int kk);

    Vector3D nfaceTop(const Cell3D &Cell);
    Vector3D nfaceTop(const int ii, const int jj, const int kk);

    Vector3D nfaceBot(const Cell3D &Cell);
    Vector3D nfaceBot(const int ii, const int jj, const int kk);

    /*******************Calculate Area of each face*************************/
    double AfaceN(const Cell3D &Cell);
    double AfaceN(const int ii, const int jj, const int kk);

    double AfaceS(const Cell3D &Cell);
    double AfaceS(const int ii, const int jj, const int kk);

    double AfaceE(const Cell3D &Cell);
    double AfaceE(const int ii, const int jj, const int kk);
  
    double AfaceW(const Cell3D &Cell);
    double AfaceW(const int ii, const int jj, const int kk);
  
    double AfaceTop(const Cell3D &Cell);
    double AfaceTop(const int ii, const int jj, const int kk);
   
    double AfaceBot(const Cell3D &Cell);
    double AfaceBot(const int ii, const int jj, const int kk);

    /* Assignment operator. */
    // Grid2D_Quad_Block operator = 
    //    (const Grid2D_Quad_Block &G);
    // Use automatically generated assignment operator.

    /* Binary arithmetic operators. */
    friend Grid3D_Hexa_Block operator +(Grid3D_Hexa_Block &G, const Vector3D &V);
    friend Grid3D_Hexa_Block operator -(Grid3D_Hexa_Block &G, const Vector3D &V);
    friend Grid3D_Hexa_Block operator *(Grid3D_Hexa_Block &G, const double &a);
    friend Grid3D_Hexa_Block operator *(const double &a, Grid3D_Hexa_Block &G);
    friend Grid3D_Hexa_Block operator /(Grid3D_Hexa_Block &G, const double &a);
    friend Grid3D_Hexa_Block operator ^(Grid3D_Hexa_Block &G, const double &a);
    friend Grid3D_Hexa_Block operator &(Grid3D_Hexa_Block &G, const double &a);
    friend Grid3D_Hexa_Block operator %(Grid3D_Hexa_Block &G, const double &a);

    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file, const Grid3D_Hexa_Block &G);
    friend istream &operator >> (istream &in_file, const Grid3D_Hexa_Block &G);

};

/*************************************************************************
 * Grid3D_Hexa_Block::allocate -- Allocate memory.                       *
 *************************************************************************/
inline void Grid3D_Hexa_Block::allocate(const int Ni, const int Nj, const int Nk, const int Ng) {
  assert( Ni >= 1 && Nj >= 1 && Nk >= 1);
  NNi = Ni+(2*Ng)+1; INl = Ng; INu = Ni+Ng; NNj = Nj+(2*Ng)+1; JNl = Ng; JNu = Nj+Ng; NNk = Nk+(2*Ng)+1; KNl = Ng; KNu = Nk+Ng; 
  NCi = Ni+2*Ng; ICl = Ng; ICu = Ni+(Ng-1); NCj = Nj+(2*Ng); JCl = Ng; JCu = Nj+Ng-1; NCk = Nk+2*Ng; KCl = Ng; KCu = Nk+Ng-1;
  Nghost=Ng;
   //dynamically create 3x3 array of type node
   Node = new Node3D**[NNi];
   for (int i = 0; i <= NNi-1; ++i ){ 
     Node[i] = new Node3D*[NNj];
     for (int j = 0; j <= NNj-1; ++j ){ 
       Node[i][j] = new Node3D[NNk];
     }
   }
   //dynamically create 3x3 array of type cell
   Cell = new Cell3D**[NCi];
   for (int i = 0; i <= NCi-1 ; ++i ){ 
     Cell[i] = new Cell3D*[NCj];
     for (int j = 0; j <= NCj-1 ; ++j ){ 
       Cell[i][j] = new Cell3D[NCk];
     }
   }
}

/*************************************************************************
 * Grid3D_Hexa_Block::deallocate -- Deallocate memory.                   *
 *************************************************************************/
inline void Grid3D_Hexa_Block::deallocate(void) {
  for (int i = 0; i <= NNi-1 ; ++i ) {
    for ( int j = 0 ; j <= NNj-1 ; ++j) {
      delete []Node[i][j]; Node[i][j] = NULL;
    } /* endfor */
  }/*endfor*/ 
  for( int i = 0 ; i <= NNi-1 ; ++i ) {
    delete []Node[i]; Node[i] = NULL;
  }/*endfor*/
  delete []Node; Node = NULL;
  
  for (int i = 0; i <= NCi-1; ++i ) {
    for ( int j = 0 ; j <= NCj -1 ; ++j) {
      delete []Cell[i][j]; Cell[i][j] = NULL;
    } /* endfor */
  }/*endfor*/ 
  for( int i = 0 ; i <= NCi -1; ++i ) {
    delete []Cell[i]; Cell[i] = NULL;
  }/*endfor*/
  delete []Cell; Cell = NULL;
  Nghost=0;
}

/*************************************************************************
 * Grid3D_Hexa_Block::allocateNodes -- Allocate nodes.                   *
 *************************************************************************/
inline void Grid3D_Hexa_Block::allocateNodes(const int Ni, const int Nj, const int Nk, const int Ng) {
   assert( Ni > 1 && Nj > 1 && Nk  > 1);
   NNi = Ni+2*Ng+1; INl = Ng; INu = Ni+Ng; NNj = Nj+2*Ng+1; JNl = Ng; JNu = Nj+Ng; NNk = Nk+2*Ng+1; KNl = Ng; KNu = Nk+Ng;
   Nghost=Ng;
   Node = new Node3D**[NNi];
   for (int i = 0; i <= NNi-1 ; ++i ){ 
     Node[i] = new Node3D*[NNj];
     for (int j = 0; j <= NNj-1 ; ++j ){ 
       Node[i][j] = new Node3D[NNk];
     }
   }
}

/*************************************************************************
 * Grid3D_Hexa_Block::deallocateNodes -- Deallocate Nodes.               *
 *************************************************************************/
inline void Grid3D_Hexa_Block::deallocateNodes(void) {
  for (int i = 0; i <= NNi-1 ; ++i ) {
    for ( int j = 0 ; j <= NNj -1 ; ++j) {
      delete []Node[i][j]; Node[i][j] = NULL;
    } /* endfor */
  }/*endfor*/ 
  for( int i = 0 ; i <= NNi -1; ++i ) {
    delete []Node[i]; Node[i] = NULL;
  }/*endfor*/
  delete []Node; Node = NULL;
  NNi = 0; INl = 0; INu = 0; NNj = 0; JNl = 0; JNu = 0; NNk = 0 ; KNl = 0 ; KNu = 0 ;
  Nghost=0;
}

/*************************************************************************
 * Grid3D_Hexa_Block::allocateCells -- Allocate cells.                   *
 *************************************************************************/
inline void Grid3D_Hexa_Block::allocateCells(const int Ni, const int Nj, const int Nk, const int Ng) {
   assert( Ni > 1 && Nj > 1&& Nk >1);
   NCi = Ni+2*Ng; ICl = Ng; ICu = Ni+Ng-1; NCj = Nj+2*Ng; JCl = Ng; JCu = Nj+Ng; NCk = Nk+2*Ng; KCl = Ng; KCu = Nk+Ng-1;
   Nghost=Ng;
   //dynamically create 3x3 array of type cell
   Cell = new Cell3D**[NCi];
   for (int i = 0; i <= NCi-1 ; ++i ){ 
     Cell[i] = new Cell3D*[NCj];
     for (int j = 0; j <= NCj-1 ; ++j ){ 
       Cell[i][j] = new Cell3D[NCk];
     }
   }
}

/*************************************************************************
 * Grid3D_Hexa_Block::deallocateCells -- Deallocate Cells.               *
 *************************************************************************/
inline void Grid3D_Hexa_Block::deallocateCells(void) {
  for (int i = 0; i <= NCi-1 ; ++i ) {
    for ( int j = 0 ; j <= NCj -1 ; ++j) {
      delete []Cell[i][j]; Cell[i][j] = NULL;
    } /* endfor */
  }/*endfor*/ 
  for( int i = 0 ; i <= NCi -1; ++i ) {
  
  delete []Cell[i]; Cell[i] = NULL;
  }/*endfor*/
  delete []Cell; Cell = NULL;
  NCi = 0; ICl = 0; ICu = 0; NCj = 0; JCl = 0; JCu = 0; NCk = 0 ; KCl = 0 ; KCu = 0 ;
  Nghost=0;
}

/*************************************************************************
 * Grid3D_Hexa_Block::centroid -- Cell centre.                           *
 *************************************************************************/
inline Vector3D Grid3D_Hexa_Block::centroid(const Cell3D &Cell) {
    return ((Node[Cell.I][Cell.J+1][Cell.K].X+Node[Cell.I+1][Cell.J+1][Cell.K].X+
             Node[Cell.I+1][Cell.J][Cell.K].X+Node[Cell.I][Cell.J][Cell.K].X+Node[Cell.I][Cell.J+1][Cell.K+1].X+
	     Node[Cell.I+1][Cell.J+1][Cell.K+1].X+Node[Cell.I+1][Cell.J][Cell.K+1].X+Node[Cell.I][Cell.J][Cell.K+1].X)/8);
}

inline Vector3D Grid3D_Hexa_Block::centroid(const int ii, const int jj, const int kk) {
  return ((Node[ii][jj+1][kk].X+Node[ii+1][jj+1][kk].X+Node[ii+1][jj][kk].X+Node[ii][jj][kk].X+Node[ii][jj+1][kk+1].X+
	   Node[ii+1][jj+1][kk+1].X+Node[ii+1][jj][kk+1].X+Node[ii][jj][kk+1].X)/8);
}
/*************************************************************************
 * Grid3D_Hexa_Block::centroid -- Cell Volume.                           *
 *************************************************************************/

inline double Grid3D_Hexa_Block::volume(const Cell3D &Cell){
  return (tetvolume(centroid(Cell), nodeNWBot(Cell), nodeNEBot(Cell), xfaceBot(Cell))+
	  tetvolume(centroid(Cell), nodeNWBot(Cell), nodeSWBot(Cell), xfaceBot(Cell))+
	  tetvolume(centroid(Cell), nodeSWBot(Cell), nodeSEBot(Cell), xfaceBot(Cell))+
	  tetvolume(centroid(Cell), nodeSEBot(Cell), nodeNEBot(Cell), xfaceBot(Cell))+
	  tetvolume(centroid(Cell), nodeSEBot(Cell), nodeNEBot(Cell), xfaceE(Cell))+
	  tetvolume(centroid(Cell), nodeNEBot(Cell), nodeNETop(Cell), xfaceE(Cell))+
	  tetvolume(centroid(Cell), nodeNETop(Cell), nodeSETop(Cell), xfaceE(Cell))+
	  tetvolume(centroid(Cell), nodeSETop(Cell), nodeSEBot(Cell), xfaceE(Cell))+
	  tetvolume(centroid(Cell), nodeSETop(Cell), nodeSEBot(Cell), xfaceS(Cell))+
	  tetvolume(centroid(Cell), nodeSEBot(Cell), nodeSWBot(Cell), xfaceS(Cell))+
	  tetvolume(centroid(Cell), nodeSWBot(Cell), nodeSWTop(Cell), xfaceS(Cell))+
	  tetvolume(centroid(Cell), nodeSWTop(Cell), nodeSETop(Cell), xfaceS(Cell))+
	  tetvolume(centroid(Cell), nodeSWTop(Cell), nodeSETop(Cell), xfaceTop(Cell))+
	  tetvolume(centroid(Cell), nodeSETop(Cell), nodeNETop(Cell), xfaceTop(Cell))+
	  tetvolume(centroid(Cell), nodeNETop(Cell), nodeNWTop(Cell), xfaceTop(Cell))+
	  tetvolume(centroid(Cell), nodeNWTop(Cell), nodeSWTop(Cell), xfaceTop(Cell))+
	  tetvolume(centroid(Cell), nodeNWTop(Cell), nodeSWTop(Cell), xfaceW(Cell))+
	  tetvolume(centroid(Cell), nodeSWTop(Cell), nodeSWBot(Cell), xfaceW(Cell))+
	  tetvolume(centroid(Cell), nodeSWBot(Cell), nodeNWBot(Cell), xfaceW(Cell))+
	  tetvolume(centroid(Cell), nodeNWBot(Cell), nodeNWTop(Cell), xfaceW(Cell))+
	  tetvolume(centroid(Cell), nodeNWBot(Cell), nodeNWTop(Cell), xfaceN(Cell))+
	  tetvolume(centroid(Cell), nodeNWTop(Cell), nodeNETop(Cell), xfaceN(Cell))+
	  tetvolume(centroid(Cell), nodeNETop(Cell), nodeNEBot(Cell), xfaceN(Cell))+
	  tetvolume(centroid(Cell), nodeNEBot(Cell), nodeNWBot(Cell), xfaceN(Cell)));
}

inline double Grid3D_Hexa_Block::volume(const int ii, const int jj, const int kk){
  return (tetvolume(centroid(ii,jj,kk), nodeNWBot(ii,jj,kk), nodeNEBot(ii,jj,kk), xfaceBot(ii,jj,kk))+
	  tetvolume(centroid(ii,jj,kk), nodeNWBot(ii,jj,kk), nodeSWBot(ii,jj,kk), xfaceBot(ii,jj,kk))+
	  tetvolume(centroid(ii,jj,kk), nodeSWBot(ii,jj,kk), nodeSEBot(ii,jj,kk), xfaceBot(ii,jj,kk))+
	  tetvolume(centroid(ii,jj,kk), nodeSEBot(ii,jj,kk), nodeNEBot(ii,jj,kk), xfaceBot(ii,jj,kk))+
	  tetvolume(centroid(ii,jj,kk), nodeSEBot(ii,jj,kk), nodeNEBot(ii,jj,kk), xfaceE(ii,jj,kk))+
	  tetvolume(centroid(ii,jj,kk), nodeNEBot(ii,jj,kk), nodeNETop(ii,jj,kk), xfaceE(ii,jj,kk))+
	  tetvolume(centroid(ii,jj,kk), nodeNETop(ii,jj,kk), nodeSETop(ii,jj,kk), xfaceE(ii,jj,kk))+
	  tetvolume(centroid(ii,jj,kk), nodeSETop(ii,jj,kk), nodeSEBot(ii,jj,kk), xfaceE(ii,jj,kk))+
	  tetvolume(centroid(ii,jj,kk), nodeSETop(ii,jj,kk), nodeSEBot(ii,jj,kk), xfaceS(ii,jj,kk))+
	  tetvolume(centroid(ii,jj,kk), nodeSEBot(ii,jj,kk), nodeSWBot(ii,jj,kk), xfaceS(ii,jj,kk))+
	  tetvolume(centroid(ii,jj,kk), nodeSWBot(ii,jj,kk), nodeSWTop(ii,jj,kk), xfaceS(ii,jj,kk))+
	  tetvolume(centroid(ii,jj,kk), nodeSWTop(ii,jj,kk), nodeSETop(ii,jj,kk), xfaceS(ii,jj,kk))+
	  tetvolume(centroid(ii,jj,kk), nodeSWTop(ii,jj,kk), nodeSETop(ii,jj,kk), xfaceTop(ii,jj,kk))+
	  tetvolume(centroid(ii,jj,kk), nodeSETop(ii,jj,kk), nodeNETop(ii,jj,kk), xfaceTop(ii,jj,kk))+
	  tetvolume(centroid(ii,jj,kk), nodeNETop(ii,jj,kk), nodeNWTop(ii,jj,kk), xfaceTop(ii,jj,kk))+
	  tetvolume(centroid(ii,jj,kk), nodeNWTop(ii,jj,kk), nodeSWTop(ii,jj,kk), xfaceTop(ii,jj,kk))+
	  tetvolume(centroid(ii,jj,kk), nodeNWTop(ii,jj,kk), nodeSWTop(ii,jj,kk), xfaceW(ii,jj,kk))+
	  tetvolume(centroid(ii,jj,kk), nodeSWTop(ii,jj,kk), nodeSWBot(ii,jj,kk), xfaceW(ii,jj,kk))+
	  tetvolume(centroid(ii,jj,kk), nodeSWBot(ii,jj,kk), nodeNWBot(ii,jj,kk), xfaceW(ii,jj,kk))+
	  tetvolume(centroid(ii,jj,kk), nodeNWBot(ii,jj,kk), nodeNWTop(ii,jj,kk), xfaceW(ii,jj,kk))+
	  tetvolume(centroid(ii,jj,kk), nodeNWBot(ii,jj,kk), nodeNWTop(ii,jj,kk), xfaceN(ii,jj,kk))+
	  tetvolume(centroid(ii,jj,kk), nodeNWTop(ii,jj,kk), nodeNETop(ii,jj,kk), xfaceN(ii,jj,kk))+
	  tetvolume(centroid(ii,jj,kk), nodeNETop(ii,jj,kk), nodeNEBot(ii,jj,kk), xfaceN(ii,jj,kk))+
	  tetvolume(centroid(ii,jj,kk), nodeNEBot(ii,jj,kk), nodeNWBot(ii,jj,kk), xfaceN(ii,jj,kk)));
}  

/*************************************************************************
 * Grid3D_Hexa_Block::node?? -- Get cell nodes.                          *
 *************************************************************************/

inline Node3D Grid3D_Hexa_Block::nodeNWBot(const Cell3D &Cell) { 
  return (Node[Cell.I][Cell.J+1][Cell.K]); 
}

inline Node3D Grid3D_Hexa_Block::nodeNWBot(const int ii, const int jj, const int kk) {
  return (Node[ii][jj+1][kk]);
} 

inline Node3D Grid3D_Hexa_Block::nodeNEBot(const Cell3D &Cell) {
  return (Node[Cell.I+1][Cell.J+1][Cell.K]);
} 
 
inline Node3D Grid3D_Hexa_Block::nodeNEBot(const int ii, const int jj, const int kk) {
  return (Node[ii+1][jj+1][kk]);
}

inline Node3D Grid3D_Hexa_Block::nodeSEBot(const Cell3D &Cell) {
  return (Node[Cell.I+1][Cell.J][Cell.K]);
}
 
inline Node3D Grid3D_Hexa_Block::nodeSEBot(const int ii, const int jj, const int kk) {
  return (Node[ii+1][jj][kk]);
} 

inline Node3D Grid3D_Hexa_Block::nodeSWBot(const Cell3D &Cell) {
   return (Node[Cell.I][Cell.J][Cell.K]);
}

inline Node3D Grid3D_Hexa_Block::nodeSWBot(const int ii, const int jj, const int kk) {
  return (Node[ii][jj][kk]);
}

inline Node3D Grid3D_Hexa_Block::nodeNWTop(const Cell3D &Cell) {
  return (Node[Cell.I][Cell.J+1][Cell.K+1]);
}

inline Node3D Grid3D_Hexa_Block::nodeNWTop(const int ii, const int jj, const int kk) {
  return (Node[ii][jj+1][kk+1]);
}

inline Node3D Grid3D_Hexa_Block::nodeNETop(const Cell3D &Cell) {
  return (Node[Cell.I+1][Cell.J+1][Cell.K+1]);
} 

inline Node3D Grid3D_Hexa_Block::nodeNETop(const int ii, const int jj, const int kk) {
  return (Node[ii+1][jj+1][kk+1]);
}

inline Node3D Grid3D_Hexa_Block::nodeSETop(const Cell3D &Cell) {
  return (Node[Cell.I+1][Cell.J][Cell.K+1]);
}

inline Node3D Grid3D_Hexa_Block::nodeSETop(const int ii, const int jj, const int kk) {
    return (Node[ii+1][jj][kk+1]);
}

inline Node3D Grid3D_Hexa_Block::nodeSWTop(const Cell3D &Cell) {
    return (Node[Cell.I][Cell.J][Cell.K+1]);
}

inline Node3D Grid3D_Hexa_Block::nodeSWTop(const int ii, const int jj, const int kk) {
  return (Node[ii][jj][kk+1]);
}

/********************************************************
 * Grid3D_Hexa_Block -- Get midpoints of faces          *
 ********************************************************/
inline Vector3D Grid3D_Hexa_Block::xfaceN(const Cell3D &Cell) {
  return ((Node[Cell.I][Cell.J+1][Cell.K].X+Node[Cell.I+1][Cell.J+1][Cell.K].X+
	   Node[Cell.I+1][Cell.J+1][Cell.K+1].X+Node[Cell.I][Cell.J+1][Cell.K+1].X)/FOUR);
}

inline Vector3D Grid3D_Hexa_Block::xfaceN(const int ii, const int jj, const int kk) {
  return ((Node[ii][jj+1][kk].X+Node[ii+1][jj+1][kk].X+
	   Node[ii+1][jj+1][kk+1].X+Node[ii][jj+1][kk+1].X)/FOUR);
}


inline Vector3D Grid3D_Hexa_Block::xfaceS(const Cell3D &Cell) {
  return ((Node[Cell.I][Cell.J][Cell.K].X+Node[Cell.I+1][Cell.J][Cell.K].X+
	   Node[Cell.I+1][Cell.J][Cell.K+1].X+Node[Cell.I][Cell.J][Cell.K+1].X)/FOUR);
}

inline Vector3D Grid3D_Hexa_Block::xfaceS(const int ii, const int jj, const int kk) {
  return ((Node[ii][jj][kk].X+Node[ii+1][jj][kk].X+
	   Node[ii+1][jj][kk+1].X+Node[ii][jj][kk+1].X)/FOUR);
}

inline Vector3D Grid3D_Hexa_Block::xfaceE(const Cell3D &Cell) {
  return ((Node[Cell.I+1][Cell.J][Cell.K].X+Node[Cell.I+1][Cell.J][Cell.K+1].X+
	   Node[Cell.I+1][Cell.J+1][Cell.K+1].X+Node[Cell.I+1][Cell.J+1][Cell.K].X)/FOUR);
}

inline Vector3D Grid3D_Hexa_Block::xfaceE(const int ii, const int jj, const int kk) {
  return ((Node[ii+1][jj][kk].X+Node[ii+1][jj][kk+1].X+
	   Node[ii+1][jj+1][kk+1].X+Node[ii+1][jj+1][kk].X)/FOUR);
}

inline Vector3D Grid3D_Hexa_Block::xfaceW(const Cell3D &Cell) {
  return ((Node[Cell.I][Cell.J][Cell.K].X+Node[Cell.I][Cell.J][Cell.K+1].X+
	   Node[Cell.I][Cell.J+1][Cell.K+1].X+Node[Cell.I][Cell.J+1][Cell.K].X)/FOUR);
}

inline Vector3D Grid3D_Hexa_Block::xfaceW(const int ii, const int jj, const int kk) {
  return ((Node[ii][jj][kk].X+Node[ii][jj][kk+1].X+
	   Node[ii][jj+1][kk+1].X+Node[ii][jj+1][kk].X)/FOUR);
}

inline Vector3D Grid3D_Hexa_Block::xfaceTop(const Cell3D &Cell) {
  return ((Node[Cell.I][Cell.J][Cell.K+1].X+Node[Cell.I+1][Cell.J][Cell.K+1].X+
	   Node[Cell.I+1][Cell.J+1][Cell.K+1].X+Node[Cell.I][Cell.J+1][Cell.K+1].X)/FOUR);
}

inline Vector3D Grid3D_Hexa_Block::xfaceTop(const int ii, const int jj, const int kk) {
  return ((Node[ii][jj][kk+1].X+Node[ii+1][jj][kk+1].X+
	   Node[ii+1][jj+1][kk+1].X+Node[ii][jj+1][kk+1].X)/FOUR);
}

inline Vector3D Grid3D_Hexa_Block::xfaceBot(const Cell3D &Cell) {
  return ((Node[Cell.I][Cell.J][Cell.K].X+Node[Cell.I+1][Cell.J][Cell.K].X+
	   Node[Cell.I+1][Cell.J+1][Cell.K].X+Node[Cell.I][Cell.J+1][Cell.K].X)/FOUR);
}

inline Vector3D Grid3D_Hexa_Block::xfaceBot(const int ii, const int jj, const int kk) {
  return ((Node[ii][jj][kk].X+Node[ii+1][jj][kk].X+
	   Node[ii+1][jj+1][kk].X+Node[ii][jj+1][kk].X)/FOUR);
}
/********************************************************
 * Grid3D_Hexa_Block -- Get normal vectors to faces     *
 ********************************************************/
inline Vector3D Grid3D_Hexa_Block::nfaceN(const Cell3D &Cell) {
  Vector3D n1, n2, n3, n4;
  n1 = ((xfaceN(Cell)-nodeNEBot(Cell).X)^(nodeNEBot(Cell).X-nodeNWBot(Cell).X));
  n2 = ((xfaceN(Cell)-nodeNWBot(Cell).X)^(nodeNWBot(Cell).X-nodeNWTop(Cell).X));
  n3 = ((xfaceN(Cell)-nodeNWTop(Cell).X)^(nodeNWTop(Cell).X-nodeNETop(Cell).X));
  n4 = ((xfaceN(Cell)-nodeNETop(Cell).X)^(nodeNETop(Cell).X-nodeNEBot(Cell).X));
  
  if(abs(n1)!=0)
    n1=n1/abs(n1);
  if(abs(n2)!=0)
    n2=n2/abs(n2);
  if(abs(n3)!=0)
    n3=n3/abs(n3);
  if(abs(n4)!=0)
    n4=n4/abs(n4);
  
  double A1, A2, A3, A4;
  A1 = abs(((xfaceN(Cell)-nodeNEBot(Cell).X)^(nodeNEBot(Cell).X-nodeNWBot(Cell).X)))/2;
  A2 = abs(((xfaceN(Cell)-nodeNWBot(Cell).X)^(nodeNWBot(Cell).X-nodeNWTop(Cell).X)))/2;
  A3 = abs(((xfaceN(Cell)-nodeNWTop(Cell).X)^(nodeNWTop(Cell).X-nodeNETop(Cell).X)))/2;
  A4 = abs(((xfaceN(Cell)-nodeNETop(Cell).X)^(nodeNETop(Cell).X-nodeNEBot(Cell).X)))/2;
 
  Vector3D normeq;
  normeq=((A1*n1+A2*n2+A3*n3+A4*n4));
  if(abs(normeq)!=0)
    normeq=normeq/abs(normeq);
  return normeq;

}

inline Vector3D Grid3D_Hexa_Block::nfaceN(const int ii, const int jj, const int kk) {
  Vector3D n1, n2, n3, n4;
  n1 = ((xfaceN(ii,jj,kk)-nodeNEBot(ii,jj,kk).X)^(nodeNEBot(ii,jj,kk).X-nodeNWBot(ii,jj,kk).X));
  n2 = ((xfaceN(ii,jj,kk)-nodeNWBot(ii,jj,kk).X)^(nodeNWBot(ii,jj,kk).X-nodeNWTop(ii,jj,kk).X));
  n3 = ((xfaceN(ii,jj,kk)-nodeNWTop(ii,jj,kk).X)^(nodeNWTop(ii,jj,kk).X-nodeNETop(ii,jj,kk).X));
  n4 = ((xfaceN(ii,jj,kk)-nodeNETop(ii,jj,kk).X)^(nodeNETop(ii,jj,kk).X-nodeNEBot(ii,jj,kk).X));
    
  if(abs(n1)!=0)
    n1=n1/abs(n1);
  if(abs(n2)!=0)
    n2=n2/abs(n2);
  if(abs(n3)!=0)
    n3=n3/abs(n3);
  if(abs(n4)!=0)
    n4=n4/abs(n4);
  
  double A1, A2, A3, A4;
  A1 = abs(((xfaceN(ii,jj,kk)-nodeNEBot(ii,jj,kk).X)^(nodeNEBot(ii,jj,kk).X-nodeNWBot(ii,jj,kk).X)))/2;
  A2 = abs(((xfaceN(ii,jj,kk)-nodeNWBot(ii,jj,kk).X)^(nodeNWBot(ii,jj,kk).X-nodeNWTop(ii,jj,kk).X)))/2;
  A3 = abs(((xfaceN(ii,jj,kk)-nodeNWTop(ii,jj,kk).X)^(nodeNWTop(ii,jj,kk).X-nodeNETop(ii,jj,kk).X)))/2;
  A4 = abs(((xfaceN(ii,jj,kk)-nodeNETop(ii,jj,kk).X)^(nodeNETop(ii,jj,kk).X-nodeNEBot(ii,jj,kk).X)))/2;
  
  Vector3D normeq;
  normeq=((A1*n1+A2*n2+A3*n3+A4*n4));
  if(abs(normeq)!=0)
    normeq=normeq/abs(normeq);
  return normeq;
}


inline Vector3D Grid3D_Hexa_Block::nfaceS(const Cell3D &Cell) {
 Vector3D n1, n2, n3, n4;
 n1 = ((xfaceS(Cell)-nodeSEBot(Cell).X)^(nodeSWBot(Cell).X-nodeSEBot(Cell).X));
 n2 = ((xfaceS(Cell)-nodeSWBot(Cell).X)^(nodeSWTop(Cell).X-nodeSWBot(Cell).X));
 n3 = ((xfaceS(Cell)-nodeSWTop(Cell).X)^(nodeSETop(Cell).X-nodeSWTop(Cell).X));
 n4 = ((xfaceS(Cell)-nodeSETop(Cell).X)^(nodeSEBot(Cell).X-nodeSETop(Cell).X));
 if(abs(n1)!=0)
   n1=n1/abs(n1);
 if(abs(n2)!=0)
   n2=n2/abs(n2);
 if(abs(n3)!=0)
   n3=n3/abs(n3);
 if(abs(n4)!=0)
   n4=n4/abs(n4);
 
 double A1, A2, A3, A4;
 A1 = abs(((xfaceS(Cell)-nodeSEBot(Cell).X)^(nodeSEBot(Cell).X-nodeSWBot(Cell).X)))/2;
 A2 = abs(((xfaceS(Cell)-nodeSWBot(Cell).X)^(nodeSWBot(Cell).X-nodeSWTop(Cell).X)))/2;
 A3 = abs(((xfaceS(Cell)-nodeSWTop(Cell).X)^(nodeSWTop(Cell).X-nodeSETop(Cell).X)))/2;
 A4 = abs(((xfaceS(Cell)-nodeSETop(Cell).X)^(nodeSETop(Cell).X-nodeSEBot(Cell).X)))/2;
 
 Vector3D normeq;
 normeq=((A1*n1+A2*n2+A3*n3+A4*n4));
 if(abs(normeq)!=0)
   normeq=normeq/abs(normeq);
 return normeq;
}

inline Vector3D Grid3D_Hexa_Block::nfaceS(const int ii, const int jj, const int kk) {
  Vector3D n1, n2, n3, n4;
  n1 = (((xfaceS(ii,jj,kk)-nodeSEBot(ii,jj,kk).X)^(nodeSWBot(ii,jj,kk).X-nodeSEBot(ii,jj,kk).X)));
  n2 = (((xfaceS(ii,jj,kk)-nodeSWBot(ii,jj,kk).X)^(nodeSWTop(ii,jj,kk).X-nodeSWBot(ii,jj,kk).X)));
  n3 = (((xfaceS(ii,jj,kk)-nodeSWTop(ii,jj,kk).X)^(nodeSETop(ii,jj,kk).X-nodeSWTop(ii,jj,kk).X)));
  n4 = (((xfaceS(ii,jj,kk)-nodeSETop(ii,jj,kk).X)^(nodeSEBot(ii,jj,kk).X-nodeSETop(ii,jj,kk).X))); 
  if(abs(n1)!=0)
    n1=n1/abs(n1);
  if(abs(n2)!=0)
    n2=n2/abs(n2);
  if(abs(n3)!=0)
    n3=n3/abs(n3);
  if(abs(n4)!=0)
    n4=n4/abs(n4);
  
  double A1, A2, A3, A4;
  A1 = abs(((xfaceS(ii,jj,kk)-nodeSEBot(ii,jj,kk).X)^(nodeSEBot(ii,jj,kk).X-nodeSWBot(ii,jj,kk).X)))/2;
  A2 = abs(((xfaceS(ii,jj,kk)-nodeSWBot(ii,jj,kk).X)^(nodeSWBot(ii,jj,kk).X-nodeSWTop(ii,jj,kk).X)))/2;
  A3 = abs(((xfaceS(ii,jj,kk)-nodeSWTop(ii,jj,kk).X)^(nodeSWTop(ii,jj,kk).X-nodeSETop(ii,jj,kk).X)))/2;
  A4 = abs(((xfaceS(ii,jj,kk)-nodeSETop(ii,jj,kk).X)^(nodeSETop(ii,jj,kk).X-nodeSEBot(ii,jj,kk).X)))/2;
  
  Vector3D normeq;
  normeq=((A1*n1+A2*n2+A3*n3+A4*n4));
  if(abs(normeq)!=0)
    normeq=normeq/abs(normeq);
  return normeq;
}

inline Vector3D Grid3D_Hexa_Block::nfaceE(const Cell3D &Cell) {
 Vector3D n1, n2, n3, n4;
 n1 = ((xfaceE(Cell)-nodeSEBot(Cell).X)^(nodeSEBot(Cell).X-nodeNEBot(Cell).X));	
 n2 = ((xfaceE(Cell)-nodeNEBot(Cell).X)^(nodeNEBot(Cell).X-nodeNETop(Cell).X));
 n3 = ((xfaceE(Cell)-nodeNETop(Cell).X)^(nodeNETop(Cell).X-nodeSETop(Cell).X));
 n4 = ((xfaceE(Cell)-nodeSETop(Cell).X)^(nodeSETop(Cell).X-nodeSEBot(Cell).X));

 if(abs(n1)!=0)
    n1=n1/abs(n1);
 if(abs(n2)!=0)
   n2=n2/abs(n2);
 if(abs(n3)!=0)
   n3=n3/abs(n3);
 if(abs(n4)!=0)
   n4=n4/abs(n4);
 
 double A1, A2, A3, A4;
 A1 = abs(((xfaceE(Cell)-nodeSEBot(Cell).X)^(nodeSEBot(Cell).X-nodeNEBot(Cell).X)))/2;
 A2 = abs(((xfaceE(Cell)-nodeNEBot(Cell).X)^(nodeNEBot(Cell).X-nodeNETop(Cell).X)))/2;
 A3 = abs(((xfaceE(Cell)-nodeNETop(Cell).X)^(nodeNETop(Cell).X-nodeSETop(Cell).X)))/2;
 A4 = abs(((xfaceE(Cell)-nodeSETop(Cell).X)^(nodeSETop(Cell).X-nodeSEBot(Cell).X)))/2;

 
 Vector3D normeq;
 normeq=((A1*n1+A2*n2+A3*n3+A4*n4));
 if(abs(normeq)!=0)
   normeq=normeq/abs(normeq);
 return normeq;
}

inline Vector3D Grid3D_Hexa_Block::nfaceE(const int ii, const int jj, const int kk) {
  Vector3D n1, n2, n3, n4;
  n1 = ((xfaceE(ii,jj,kk)-nodeSEBot(ii,jj,kk).X)^(nodeSEBot(ii,jj,kk).X-nodeNEBot(ii,jj,kk).X));
  n2 = ((xfaceE(ii,jj,kk)-nodeNEBot(ii,jj,kk).X)^(nodeNEBot(ii,jj,kk).X-nodeNETop(ii,jj,kk).X));
  n3 = ((xfaceE(ii,jj,kk)-nodeNETop(ii,jj,kk).X)^(nodeNETop(ii,jj,kk).X-nodeSETop(ii,jj,kk).X));
  n4 = ((xfaceE(ii,jj,kk)-nodeSETop(ii,jj,kk).X)^(nodeSETop(ii,jj,kk).X-nodeSEBot(ii,jj,kk).X));
  if(abs(n1)!=0)
    n1=n1/abs(n1);
  if(abs(n2)!=0)
    n2=n2/abs(n2);
  if(abs(n3)!=0)
    n3=n3/abs(n3);
  if(abs(n4)!=0)
    n4=n4/abs(n4);
  
  double A1, A2, A3, A4;
  A1 = abs(((xfaceE(ii,jj,kk)-nodeSEBot(ii,jj,kk).X)^(nodeSEBot(ii,jj,kk).X-nodeNEBot(ii,jj,kk).X)))/2;
  A2 = abs(((xfaceE(ii,jj,kk)-nodeNEBot(ii,jj,kk).X)^(nodeNEBot(ii,jj,kk).X-nodeNETop(ii,jj,kk).X)))/2;
  A3 = abs(((xfaceE(ii,jj,kk)-nodeNETop(ii,jj,kk).X)^(nodeNETop(ii,jj,kk).X-nodeSETop(ii,jj,kk).X)))/2;
  A4 = abs(((xfaceE(ii,jj,kk)-nodeSETop(ii,jj,kk).X)^(nodeSETop(ii,jj,kk).X-nodeSEBot(ii,jj,kk).X)))/2;
  
  Vector3D normeq;
  normeq=((A1*n1+A2*n2+A3*n3+A4*n4));
  if(abs(normeq)!=0)
    normeq=normeq/abs(normeq);
  return normeq; 
}

inline Vector3D Grid3D_Hexa_Block::nfaceW(const Cell3D &Cell) {
  Vector3D n1, n2, n3, n4;
  n1 = ((xfaceW(Cell)-nodeSWBot(Cell).X)^(nodeNWBot(Cell).X-nodeSWBot(Cell).X));
  n2 = ((xfaceW(Cell)-nodeNWBot(Cell).X)^(nodeNWTop(Cell).X-nodeNWBot(Cell).X));
  n3 = ((xfaceW(Cell)-nodeNWTop(Cell).X)^(nodeSWTop(Cell).X-nodeNWTop(Cell).X));
  n4 = ((xfaceW(Cell)-nodeSWTop(Cell).X)^(nodeSWBot(Cell).X-nodeSWTop(Cell).X));
  if(abs(n1)!=0)
    n1=n1/abs(n1);
  if(abs(n2)!=0)
    n2=n2/abs(n2);
  if(abs(n3)!=0)
    n3=n3/abs(n3);
  if(abs(n4)!=0)
    n4=n4/abs(n4);
  
  double A1, A2, A3, A4;
  A1 = abs(((xfaceW(Cell)-nodeSWBot(Cell).X)^(nodeSWBot(Cell).X-nodeNWBot(Cell).X)))/2;
  A2 = abs(((xfaceW(Cell)-nodeNWBot(Cell).X)^(nodeNWBot(Cell).X-nodeNWTop(Cell).X)))/2;
  A3 = abs(((xfaceW(Cell)-nodeNWTop(Cell).X)^(nodeNWTop(Cell).X-nodeSWTop(Cell).X)))/2;
  A4 = abs(((xfaceW(Cell)-nodeSWTop(Cell).X)^(nodeSWTop(Cell).X-nodeSWBot(Cell).X)))/2;
  
  Vector3D normeq;
  normeq=((A1*n1+A2*n2+A3*n3+A4*n4));
  if(abs(normeq)!=0)
    normeq=normeq/abs(normeq); 
  return normeq;
  
}

inline Vector3D Grid3D_Hexa_Block::nfaceW(const int ii, const int jj, const int kk) {
Vector3D n1, n2, n3, n4;
  n1 = ((xfaceW(ii,jj,kk)-nodeSWBot(ii,jj,kk).X)^(nodeNWBot(ii,jj,kk).X-nodeSWBot(ii,jj,kk).X));
  n2 = ((xfaceW(ii,jj,kk)-nodeNWBot(ii,jj,kk).X)^(nodeNWTop(ii,jj,kk).X-nodeNWBot(ii,jj,kk).X));
  n3 = ((xfaceW(ii,jj,kk)-nodeNWTop(ii,jj,kk).X)^(nodeSWTop(ii,jj,kk).X-nodeNWTop(ii,jj,kk).X));
  n4 = ((xfaceW(ii,jj,kk)-nodeSWTop(ii,jj,kk).X)^(nodeSWBot(ii,jj,kk).X-nodeSWTop(ii,jj,kk).X));
  if(abs(n1)!=0)
    n1=n1/abs(n1);
  if(abs(n2)!=0)
    n2=n2/abs(n2);
  if(abs(n3)!=0)
    n3=n3/abs(n3);
  if(abs(n4)!=0)
    n4=n4/abs(n4);
  
  double A1, A2, A3, A4;
  A1 = abs(((xfaceW(ii,jj,kk)-nodeSWBot(ii,jj,kk).X)^(nodeSWBot(ii,jj,kk).X-nodeNWBot(ii,jj,kk).X)))/2;
  A2 = abs(((xfaceW(ii,jj,kk)-nodeNWBot(ii,jj,kk).X)^(nodeNWBot(ii,jj,kk).X-nodeNWTop(ii,jj,kk).X)))/2;
  A3 = abs(((xfaceW(ii,jj,kk)-nodeNWTop(ii,jj,kk).X)^(nodeNWTop(ii,jj,kk).X-nodeSWTop(ii,jj,kk).X)))/2;
  A4 = abs(((xfaceW(ii,jj,kk)-nodeSWTop(ii,jj,kk).X)^(nodeSWTop(ii,jj,kk).X-nodeSWBot(ii,jj,kk).X)))/2;
  
  Vector3D normeq;
  normeq=((A1*n1+A2*n2+A3*n3+A4*n4));
  if(abs(normeq)!=0)
    normeq=normeq/abs(normeq); 
  return normeq;
}

inline Vector3D Grid3D_Hexa_Block::nfaceTop(const Cell3D &Cell) {
  Vector3D n1, n2, n3, n4;
  n1 = ((xfaceTop(Cell)-nodeSWTop(Cell).X)^(nodeNWTop(Cell).X-nodeSWTop(Cell).X));
  n2 = ((xfaceTop(Cell)-nodeNWTop(Cell).X)^(nodeNETop(Cell).X-nodeNWTop(Cell).X));
  n3 = ((xfaceTop(Cell)-nodeNETop(Cell).X)^(nodeSETop(Cell).X-nodeNETop(Cell).X));
  n4 = ((xfaceTop(Cell)-nodeSETop(Cell).X)^(nodeSWTop(Cell).X-nodeSETop(Cell).X));
  if(abs(n1)!=0)
    n1=n1/abs(n1);
  if(abs(n2)!=0)
    n2=n2/abs(n2);
  if(abs(n3)!=0)
    n3=n3/abs(n3);
  if(abs(n4)!=0)
    n4=n4/abs(n4);

  double A1, A2, A3, A4;
  A1 = abs(((xfaceTop(Cell)-nodeSWTop(Cell).X)^(nodeSWTop(Cell).X-nodeNWTop(Cell).X)))/2;
  A2 = abs(((xfaceTop(Cell)-nodeNWTop(Cell).X)^(nodeNWTop(Cell).X-nodeNETop(Cell).X)))/2;
  A3 = abs(((xfaceTop(Cell)-nodeNETop(Cell).X)^(nodeNETop(Cell).X-nodeSETop(Cell).X)))/2;
  A4 = abs(((xfaceTop(Cell)-nodeSETop(Cell).X)^(nodeSETop(Cell).X-nodeSWTop(Cell).X)))/2;
  
  
  Vector3D normeq;
  normeq=((A1*n1+A2*n2+A3*n3+A4*n4));
  if(abs(normeq)!=0)
    normeq=normeq/abs(normeq); 
  return normeq;
}

inline Vector3D Grid3D_Hexa_Block::nfaceTop(const int ii, const int jj, const int kk) {
  Vector3D n1, n2, n3, n4;
  n1 = ((xfaceTop(ii,jj,kk)-nodeSWTop(ii,jj,kk).X)^(nodeNWTop(ii,jj,kk).X-nodeSWTop(ii,jj,kk).X));
  n2 = ((xfaceTop(ii,jj,kk)-nodeNWTop(ii,jj,kk).X)^(nodeNETop(ii,jj,kk).X-nodeNWTop(ii,jj,kk).X));
  n3 = ((xfaceTop(ii,jj,kk)-nodeNETop(ii,jj,kk).X)^(nodeSETop(ii,jj,kk).X-nodeNETop(ii,jj,kk).X));
  n4 = ((xfaceTop(ii,jj,kk)-nodeSETop(ii,jj,kk).X)^(nodeSWTop(ii,jj,kk).X-nodeSETop(ii,jj,kk).X));
  if(abs(n1)!=0)
    n1=n1/abs(n1);
  if(abs(n2)!=0)
    n2=n2/abs(n2);
  if(abs(n3)!=0)
    n3=n3/abs(n3);
  if(abs(n4)!=0)
    n4=n4/abs(n4);

  double A1, A2, A3, A4;
  A1 = abs(((xfaceTop(ii,jj,kk)-nodeSWTop(ii,jj,kk).X)^(nodeSWTop(ii,jj,kk).X-nodeNWTop(ii,jj,kk).X)))/2;
  A2 = abs(((xfaceTop(ii,jj,kk)-nodeNWTop(ii,jj,kk).X)^(nodeNWTop(ii,jj,kk).X-nodeNETop(ii,jj,kk).X)))/2;
  A3 = abs(((xfaceTop(ii,jj,kk)-nodeNETop(ii,jj,kk).X)^(nodeNETop(ii,jj,kk).X-nodeSETop(ii,jj,kk).X)))/2;
  A4 = abs(((xfaceTop(ii,jj,kk)-nodeSETop(ii,jj,kk).X)^(nodeSETop(ii,jj,kk).X-nodeSWTop(ii,jj,kk).X)))/2;
  
  
  Vector3D normeq;
  normeq=((A1*n1+A2*n2+A3*n3+A4*n4));
  if(abs(normeq)!=0)
    normeq=normeq/abs(normeq); 
  return normeq;
}



inline Vector3D Grid3D_Hexa_Block::nfaceBot(const Cell3D &Cell) {
  Vector3D n1, n2, n3, n4;
  n1 = ((xfaceBot(Cell)-nodeSWBot(Cell).X)^(nodeSWBot(Cell).X-nodeNWBot(Cell).X));
  n2 = ((xfaceBot(Cell)-nodeNWBot(Cell).X)^(nodeNWBot(Cell).X-nodeNEBot(Cell).X));
  n3 = ((xfaceBot(Cell)-nodeNEBot(Cell).X)^(nodeNEBot(Cell).X-nodeSEBot(Cell).X));
  n4 = ((xfaceBot(Cell)-nodeSEBot(Cell).X)^(nodeSEBot(Cell).X-nodeSWBot(Cell).X));
  if(abs(n1)!=0)
    n1=n1/abs(n1);
  if(abs(n2)!=0)
    n2=n2/abs(n2);
  if(abs(n3)!=0)
    n3=n3/abs(n3);
  if(abs(n4)!=0)
    n4=n4/abs(n4);
  double A1, A2, A3, A4;
  A1 = abs(((xfaceBot(Cell)-nodeSWBot(Cell).X)^(nodeSWBot(Cell).X-nodeNWBot(Cell).X)))/2;
  A2 = abs(((xfaceBot(Cell)-nodeNWBot(Cell).X)^(nodeNWBot(Cell).X-nodeNEBot(Cell).X)))/2;
  A3 = abs(((xfaceBot(Cell)-nodeNEBot(Cell).X)^(nodeNEBot(Cell).X-nodeSEBot(Cell).X)))/2;
  A4 = abs(((xfaceBot(Cell)-nodeSEBot(Cell).X)^(nodeSEBot(Cell).X-nodeSWBot(Cell).X)))/2;
  
  
  Vector3D normeq;
  normeq=((A1*n1+A2*n2+A3*n3+A4*n4));
  if(abs(normeq)!=0)
    normeq=normeq/abs(normeq); 
  return normeq;
}

inline Vector3D Grid3D_Hexa_Block::nfaceBot(const int ii, const int jj, const int kk) {
 Vector3D n1, n2, n3, n4;
  n1 = ((xfaceBot(ii,jj,kk)-nodeSWBot(ii,jj,kk).X)^(nodeSWBot(ii,jj,kk).X-nodeNWBot(ii,jj,kk).X));
  n2 = ((xfaceBot(ii,jj,kk)-nodeNWBot(ii,jj,kk).X)^(nodeNWBot(ii,jj,kk).X-nodeNEBot(ii,jj,kk).X));
  n3 = ((xfaceBot(ii,jj,kk)-nodeNEBot(ii,jj,kk).X)^(nodeNEBot(ii,jj,kk).X-nodeSEBot(ii,jj,kk).X));
  n4 = ((xfaceBot(ii,jj,kk)-nodeSEBot(ii,jj,kk).X)^(nodeSEBot(ii,jj,kk).X-nodeSWBot(ii,jj,kk).X));
  if(abs(n1)!=0)
    n1=n1/abs(n1);
  if(abs(n2)!=0)
    n2=n2/abs(n2);
  if(abs(n3)!=0)
    n3=n3/abs(n3);
  if(abs(n4)!=0)
    n4=n4/abs(n4); 
  
  double A1, A2, A3, A4;
  A1 = abs(((xfaceBot(ii,jj,kk)-nodeSWBot(ii,jj,kk).X)^(nodeSWBot(ii,jj,kk).X-nodeNWBot(ii,jj,kk).X)))/2;
  A2 = abs(((xfaceBot(ii,jj,kk)-nodeNWBot(ii,jj,kk).X)^(nodeNWBot(ii,jj,kk).X-nodeNEBot(ii,jj,kk).X)))/2;
  A3 = abs(((xfaceBot(ii,jj,kk)-nodeNEBot(ii,jj,kk).X)^(nodeNEBot(ii,jj,kk).X-nodeSEBot(ii,jj,kk).X)))/2;
  A4 = abs(((xfaceBot(ii,jj,kk)-nodeSEBot(ii,jj,kk).X)^(nodeSEBot(ii,jj,kk).X-nodeSWBot(ii,jj,kk).X)))/2;
  
  Vector3D normeq;
  normeq=((A1*n1+A2*n2+A3*n3+A4*n4));
  if(abs(normeq)!=0)
    normeq=normeq/abs(normeq); 
  return normeq; 
}
/********************************************************
 * Grid3D_Hexa_Block -- Get area of each face           *
 ********************************************************/
inline double Grid3D_Hexa_Block::AfaceN(const Cell3D &Cell){
  double A1, A2, A3, A4;
  A1 = abs(((xfaceN(Cell)-nodeNEBot(Cell).X)^(nodeNEBot(Cell).X-nodeNWBot(Cell).X)))/2;
  A2 = abs(((xfaceN(Cell)-nodeNWBot(Cell).X)^(nodeNWBot(Cell).X-nodeNWTop(Cell).X)))/2;
  A3 = abs(((xfaceN(Cell)-nodeNWTop(Cell).X)^(nodeNWTop(Cell).X-nodeNETop(Cell).X)))/2;
  A4 = abs(((xfaceN(Cell)-nodeNETop(Cell).X)^(nodeNETop(Cell).X-nodeNEBot(Cell).X)))/2;
  return A1+A2+A3+A4;
}

inline double Grid3D_Hexa_Block::AfaceN(int const ii, int const jj, int const kk){
  double A1, A2, A3, A4;
  A1 = abs(((xfaceN(ii,jj,kk)-nodeNEBot(ii,jj,kk).X)^(nodeNEBot(ii,jj,kk).X-nodeNWBot(ii,jj,kk).X)))/2;
  A2 = abs(((xfaceN(ii,jj,kk)-nodeNWBot(ii,jj,kk).X)^(nodeNWBot(ii,jj,kk).X-nodeNWTop(ii,jj,kk).X)))/2;
  A3 = abs(((xfaceN(ii,jj,kk)-nodeNWTop(ii,jj,kk).X)^(nodeNWTop(ii,jj,kk).X-nodeNETop(ii,jj,kk).X)))/2;
  A4 = abs(((xfaceN(ii,jj,kk)-nodeNETop(ii,jj,kk).X)^(nodeNETop(ii,jj,kk).X-nodeNEBot(ii,jj,kk).X)))/2;
  return A1+A2+A3+A4;
}

inline double Grid3D_Hexa_Block::AfaceS(const Cell3D &Cell){
  double A1, A2, A3, A4;
  A1 = abs(((xfaceS(Cell)-nodeSEBot(Cell).X)^(nodeSEBot(Cell).X-nodeSWBot(Cell).X)))/2;
  A2 = abs(((xfaceS(Cell)-nodeSWBot(Cell).X)^(nodeSWBot(Cell).X-nodeSWTop(Cell).X)))/2;
  A3 = abs(((xfaceS(Cell)-nodeSWTop(Cell).X)^(nodeSWTop(Cell).X-nodeSETop(Cell).X)))/2;
  A4 = abs(((xfaceS(Cell)-nodeSETop(Cell).X)^(nodeSETop(Cell).X-nodeSEBot(Cell).X)))/2;
  return A1+A2+A3+A4;
}

inline double Grid3D_Hexa_Block::AfaceS(int const ii, int const jj, int const kk){
  double A1, A2, A3, A4;
  A1 = abs(((xfaceS(ii,jj,kk)-nodeSEBot(ii,jj,kk).X)^(nodeSEBot(ii,jj,kk).X-nodeSWBot(ii,jj,kk).X)))/2;
  A2 = abs(((xfaceS(ii,jj,kk)-nodeSWBot(ii,jj,kk).X)^(nodeSWBot(ii,jj,kk).X-nodeSWTop(ii,jj,kk).X)))/2;
  A3 = abs(((xfaceS(ii,jj,kk)-nodeSWTop(ii,jj,kk).X)^(nodeSWTop(ii,jj,kk).X-nodeSETop(ii,jj,kk).X)))/2;
  A4 = abs(((xfaceS(ii,jj,kk)-nodeSETop(ii,jj,kk).X)^(nodeSETop(ii,jj,kk).X-nodeSEBot(ii,jj,kk).X)))/2;
  return A1+A2+A3+A4;
}

inline double Grid3D_Hexa_Block::AfaceE(const Cell3D &Cell){
  double A1, A2, A3, A4;
  A1 = abs(((xfaceE(Cell)-nodeSEBot(Cell).X)^(nodeSEBot(Cell).X-nodeNEBot(Cell).X)))/2;
  A2 = abs(((xfaceE(Cell)-nodeNEBot(Cell).X)^(nodeNEBot(Cell).X-nodeNETop(Cell).X)))/2;
  A3 = abs(((xfaceE(Cell)-nodeNETop(Cell).X)^(nodeNETop(Cell).X-nodeSETop(Cell).X)))/2;
  A4 = abs(((xfaceE(Cell)-nodeSETop(Cell).X)^(nodeSETop(Cell).X-nodeSEBot(Cell).X)))/2;
  return A1+A2+A3+A4; 
}

inline double Grid3D_Hexa_Block::AfaceE(int const ii, int const jj, int const kk){
  double A1, A2, A3, A4;
  A1 = abs(((xfaceE(ii,jj,kk)-nodeSEBot(ii,jj,kk).X)^(nodeSEBot(ii,jj,kk).X-nodeNEBot(ii,jj,kk).X)))/2;
  A2 = abs(((xfaceE(ii,jj,kk)-nodeNEBot(ii,jj,kk).X)^(nodeNEBot(ii,jj,kk).X-nodeNETop(ii,jj,kk).X)))/2;
  A3 = abs(((xfaceE(ii,jj,kk)-nodeNETop(ii,jj,kk).X)^(nodeNETop(ii,jj,kk).X-nodeSETop(ii,jj,kk).X)))/2;
  A4 = abs(((xfaceE(ii,jj,kk)-nodeSETop(ii,jj,kk).X)^(nodeSETop(ii,jj,kk).X-nodeSEBot(ii,jj,kk).X)))/2;
  return A1+A2+A3+A4; 
}

inline double Grid3D_Hexa_Block::AfaceW(const Cell3D &Cell){
  double A1, A2, A3, A4;
  A1 = abs(((xfaceE(Cell)-nodeSWBot(Cell).X)^(nodeSWBot(Cell).X-nodeNWBot(Cell).X)))/2;
  A2 = abs(((xfaceE(Cell)-nodeNWBot(Cell).X)^(nodeNWBot(Cell).X-nodeNWTop(Cell).X)))/2;
  A3 = abs(((xfaceE(Cell)-nodeNWTop(Cell).X)^(nodeNWTop(Cell).X-nodeSWTop(Cell).X)))/2;
  A4 = abs(((xfaceE(Cell)-nodeSWTop(Cell).X)^(nodeSWTop(Cell).X-nodeSWBot(Cell).X)))/2;
  return A1+A2+A3+A4; 
}

inline double Grid3D_Hexa_Block::AfaceW(int const ii, int const jj, int const kk){
  double A1, A2, A3, A4;
  A1 = abs(((xfaceW(ii,jj,kk)-nodeSWBot(ii,jj,kk).X)^(nodeSWBot(ii,jj,kk).X-nodeNWBot(ii,jj,kk).X)))/2;
  A2 = abs(((xfaceW(ii,jj,kk)-nodeNWBot(ii,jj,kk).X)^(nodeNWBot(ii,jj,kk).X-nodeNWTop(ii,jj,kk).X)))/2;
  A3 = abs(((xfaceW(ii,jj,kk)-nodeNWTop(ii,jj,kk).X)^(nodeNWTop(ii,jj,kk).X-nodeSWTop(ii,jj,kk).X)))/2;
  A4 = abs(((xfaceW(ii,jj,kk)-nodeSWTop(ii,jj,kk).X)^(nodeSWTop(ii,jj,kk).X-nodeSWBot(ii,jj,kk).X)))/2;
  return A1+A2+A3+A4;  
}

inline double Grid3D_Hexa_Block::AfaceBot(const Cell3D &Cell){
  double A1, A2, A3, A4;
  A1 = abs(((xfaceBot(Cell)-nodeSWBot(Cell).X)^(nodeSWBot(Cell).X-nodeNWBot(Cell).X)))/2;
  A2 = abs(((xfaceBot(Cell)-nodeNWBot(Cell).X)^(nodeNWBot(Cell).X-nodeNEBot(Cell).X)))/2;
  A3 = abs(((xfaceBot(Cell)-nodeNEBot(Cell).X)^(nodeNEBot(Cell).X-nodeSEBot(Cell).X)))/2;
  A4 = abs(((xfaceBot(Cell)-nodeSEBot(Cell).X)^(nodeSEBot(Cell).X-nodeSWBot(Cell).X)))/2;
  return A1+A2+A3+A4; 
}

inline double Grid3D_Hexa_Block::AfaceBot(int const ii, int const jj, int const kk){
  double A1, A2, A3, A4;
  A1 = abs(((xfaceBot(ii,jj,kk)-nodeSWBot(ii,jj,kk).X)^(nodeSWBot(ii,jj,kk).X-nodeNWBot(ii,jj,kk).X)))/2;
  A2 = abs(((xfaceBot(ii,jj,kk)-nodeNWBot(ii,jj,kk).X)^(nodeNWBot(ii,jj,kk).X-nodeNEBot(ii,jj,kk).X)))/2;
  A3 = abs(((xfaceBot(ii,jj,kk)-nodeNEBot(ii,jj,kk).X)^(nodeNEBot(ii,jj,kk).X-nodeSEBot(ii,jj,kk).X)))/2;
  A4 = abs(((xfaceBot(ii,jj,kk)-nodeSEBot(ii,jj,kk).X)^(nodeSEBot(ii,jj,kk).X-nodeSWBot(ii,jj,kk).X)))/2;
  return A1+A2+A3+A4;  
}

inline double Grid3D_Hexa_Block::AfaceTop(const Cell3D &Cell){
  double A1, A2, A3, A4;
  A1 = abs(((xfaceTop(Cell)-nodeSWTop(Cell).X)^(nodeSWTop(Cell).X-nodeNWTop(Cell).X)))/2;
  A2 = abs(((xfaceTop(Cell)-nodeNWTop(Cell).X)^(nodeNWTop(Cell).X-nodeNETop(Cell).X)))/2;
  A3 = abs(((xfaceTop(Cell)-nodeNETop(Cell).X)^(nodeNETop(Cell).X-nodeSETop(Cell).X)))/2;
  A4 = abs(((xfaceTop(Cell)-nodeSETop(Cell).X)^(nodeSETop(Cell).X-nodeSWTop(Cell).X)))/2;
  return A1+A2+A3+A4; 
}

inline double Grid3D_Hexa_Block::AfaceTop(int const ii, int const jj, int const kk){
  double A1, A2, A3, A4;
  A1 = abs(((xfaceTop(ii,jj,kk)-nodeSWTop(ii,jj,kk).X)^(nodeSWTop(ii,jj,kk).X-nodeNWTop(ii,jj,kk).X)))/2;
  A2 = abs(((xfaceTop(ii,jj,kk)-nodeNWTop(ii,jj,kk).X)^(nodeNWTop(ii,jj,kk).X-nodeNETop(ii,jj,kk).X)))/2;
  A3 = abs(((xfaceTop(ii,jj,kk)-nodeNETop(ii,jj,kk).X)^(nodeNETop(ii,jj,kk).X-nodeSETop(ii,jj,kk).X)))/2;
  A4 = abs(((xfaceTop(ii,jj,kk)-nodeSETop(ii,jj,kk).X)^(nodeSETop(ii,jj,kk).X-nodeSWTop(ii,jj,kk).X)))/2;
  return A1+A2+A3+A4;   
}

/********************************************************
 * Grid3D_Hexa_Block -- Binary arithmetic operators.    *
 ********************************************************/
// Shift operators.
inline Grid3D_Hexa_Block operator +(Grid3D_Hexa_Block &G, const Vector3D &V) {
  for(int k = G.KNl - G.Nghost; k <= G.KNu + G.Nghost; ++k ){
    for (int  j = G.JNl-G.Nghost; j <= G.JNu+G.Nghost; ++j ) {
      for (int i = G.INl-G.Nghost ; i <= G.INu+G.Nghost; ++i ) {
	G.Node[i][j][k].X += V;
      } /* endfor */
    } /* endfor */
  }/*endfor */
  
  for(int k = G.KCl - G.Nghost; k<= G.KCu + G.Nghost; ++k){
    for (int j = G.JCl-G.Nghost ; j <= G.JCu+G.Nghost ; ++j) {
      for (int i = G.ICl-G.Nghost ; i <= G.ICu+G.Nghost ; ++i) {
	G.Cell[i][j][k].Xc = G.centroid(i, j, k);
      } /* endfor */
    } /* endfor */
  }/*endfor */

  return (G);
}

inline Grid3D_Hexa_Block operator -(Grid3D_Hexa_Block &G, const Vector3D &V) {
  for(int k = G.KNl-G.Nghost; k<=G.KNu+G.Nghost; ++k){
    for (int  j = G.JNl-G.Nghost; j <= G.JNu+G.Nghost; ++j ) {
      for (int i = G.INl-G.Nghost ; i <= G.INu+G.Nghost; ++i ) {
	G.Node[i][j][k].X -= V;
      } /* endfor */
    } /* endfor */
  }/*endfor */
  
  for(int k = G.KCl - G.Nghost; k<= G.KCu + G.Nghost; ++k){
    for (int j = G.JCl-G.Nghost ; j <= G.JCu+G.Nghost ; ++j) {
      for (int i = G.ICl-G.Nghost ; i <= G.ICu+G.Nghost ; ++i) {
	G.Cell[i][j][k].Xc = G.centroid(i, j, k);
      } /* endfor */
    } /* endfor */
  }/*endfor */
  
  return (G);
}

// Scaling operators.
inline Grid3D_Hexa_Block operator *(Grid3D_Hexa_Block &G, const double &a) {
  for(int k = G.KNl-G.Nghost; k<=G.KNu+G.Nghost; ++k){
    for (int  j = G.JNl-G.Nghost; j <= G.JNu+G.Nghost; ++j ) {
      for (int i = G.INl-G.Nghost ; i <= G.INu+G.Nghost; ++i ) {
	G.Node[i][j][k].X = G.Node[i][j][k].X*a;
      } /* endfor */
    } /* endfor */
  }/*endfor */
  
  for(int k = G.KCl - G.Nghost; k<= G.KCu + G.Nghost; ++k){
    for (int j = G.JCl-G.Nghost ; j <= G.JCu+G.Nghost ; ++j) {
      for (int i = G.ICl-G.Nghost ; i <= G.ICu+G.Nghost ; ++i) {
	G.Cell[i][j][k].Xc = G.centroid(i, j, k);
      } /* endfor */
    } /* endfor */
  }/*endfor */
  
  return (G);
}

inline Grid3D_Hexa_Block operator *(const double &a, Grid3D_Hexa_Block &G) {
 for(int k = G.KNl-G.Nghost; k<=G.KNu+G.Nghost; ++k){
    for (int  j = G.JNl-G.Nghost; j <= G.JNu+G.Nghost; ++j ) {
      for (int i = G.INl-G.Nghost ; i <= G.INu+G.Nghost; ++i ) {
	G.Node[i][j][k].X = G.Node[i][j][k].X*a;
      } /* endfor */
    } /* endfor */
  }/*endfor */
  
  for(int k = G.KCl - G.Nghost; k<= G.KCu + G.Nghost; ++k){
    for (int j = G.JCl-G.Nghost ; j <= G.JCu+G.Nghost ; ++j) {
      for (int i = G.ICl-G.Nghost ; i <= G.ICu+G.Nghost ; ++i) {
	G.Cell[i][j][k].Xc = G.centroid(i, j, k);
      } /* endfor */
    } /* endfor */
  }/*endfor */
  
  return (G); 
}

inline Grid3D_Hexa_Block operator /(Grid3D_Hexa_Block &G, const double &a) {
  for(int k = G.KNl-G.Nghost; k<=G.KNu+G.Nghost; ++k){
    for (int  j = G.JNl-G.Nghost; j <= G.JNu+G.Nghost; ++j ) {
      for (int i = G.INl-G.Nghost ; i <= G.INu+G.Nghost; ++i ) {
	G.Node[i][j][k].X = G.Node[i][j][k].X/a;
      } /* endfor */
    } /* endfor */
  }/*endfor */
  
  for(int k = G.KCl - G.Nghost; k<= G.KCu + G.Nghost; ++k){
    for (int j = G.JCl-G.Nghost ; j <= G.JCu+G.Nghost ; ++j) {
      for (int i = G.ICl-G.Nghost ; i <= G.ICu+G.Nghost ; ++i) {
	G.Cell[i][j][k].Xc = G.centroid(i, j, k);
      } /* endfor */
    } /* endfor */
  }/*endfor */
  
  return (G);   
}

// Rotation operator.
inline Grid3D_Hexa_Block operator ^(Grid3D_Hexa_Block &G, const double &a) {
  double cos_angle, sin_angle; Vector3D X;
  cos_angle = cos(a); sin_angle = sin(a);
  //rotation about the Z-axis
  for(int k = G.KNl - G.Nghost; k <= G.KNu + G.Nghost; ++k){
    for (int j = G.JNl-G.Nghost ; j <= G.JNu+G.Nghost; ++j ) {
      for (int i = G.INl-G.Nghost ; i <= G.INu+G.Nghost; ++i ) {
	X.x = G.Node[i][j][k].X.x*cos_angle -
	  G.Node[i][j][k].X.y*sin_angle +
	  0*G.Node[i][j][k].X.z;
	X.y = G.Node[i][j][k].X.x*sin_angle +
	  G.Node[i][j][k].X.y*cos_angle +
	  0*G.Node[i][j][k].X.z;
	X.z = 0*G.Node[i][j][k].X.x*cos_angle +
	  0*G.Node[i][j][k].X.y*sin_angle +
	  G.Node[i][j][k].X.z; 
	G.Node[i][j][k].X = X;
      } /* endfor */
    } /* endfor */
  } /*endfor*/ 
  for( int k = G.KCl-G.Nghost; k<= G.KCu + G.Nghost; ++k){
    for (int j = G.JCl-G.Nghost ; j <= G.JCu+G.Nghost ; ++j) {
      for (int i = G.ICl-G.Nghost ; i <= G.ICu+G.Nghost ; ++i) {
	G.Cell[i][j][k].Xc = G.centroid(i, j, k);
      } /* endfor */
    } /* endfor */
  } /*endfor */

  return (G); 

}

inline Grid3D_Hexa_Block operator &(Grid3D_Hexa_Block &G, const double &a) {  
  Vector3D X;
 double cos_angle=cos(a);double sin_angle= sin(a);
  //rotation about the X-axis
  for(int k = G.KNl - G.Nghost; k <= G.KNu + G.Nghost; ++k){
    for (int j = G.JNl-G.Nghost ; j <= G.JNu+G.Nghost; ++j ) {
      for (int i = G.INl-G.Nghost ; i <= G.INu+G.Nghost; ++i ) {
	X.x = G.Node[i][j][k].X.x +
	  0*G.Node[i][j][k].X.y +
	  0*G.Node[i][j][k].X.z;
	X.y = 0*G.Node[i][j][k].X.x +
	  G.Node[i][j][k].X.y*cos_angle -
	  G.Node[i][j][k].X.z*sin_angle;
	X.z = 0*G.Node[i][j][k].X.x +
	  G.Node[i][j][k].X.y*sin_angle +
	  G.Node[i][j][k].X.z*cos_angle; 
	G.Node[i][j][k].X = X;
      } /* endfor */
    } /* endfor */
  } /*endfor*/
 
  for( int k = G.KCl-G.Nghost; k<= G.KCu + G.Nghost; ++k){
    for (int j = G.JCl-G.Nghost ; j <= G.JCu+G.Nghost ; ++j) {
      for (int i = G.ICl-G.Nghost ; i <= G.ICu+G.Nghost ; ++i) {
	G.Cell[i][j][k].Xc = G.centroid(i, j, k);
      } /* endfor */
    } /* endfor */
  } /*endfor */
}
inline Grid3D_Hexa_Block operator %(Grid3D_Hexa_Block &G, const double &a) {   
  Vector3D X;
  double cos_angle=cos(a);double  sin_angle= sin(a);
  //rotation about the Y-axis
  for(int k = G.KNl - G.Nghost; k <= G.KNu + G.Nghost; ++k){
    for (int j = G.JNl-G.Nghost ; j <= G.JNu+G.Nghost; ++j ) {
      for (int i = G.INl-G.Nghost ; i <= G.INu+G.Nghost; ++i ) {
	X.x = G.Node[i][j][k].X.x*cos_angle +
	  0*G.Node[i][j][k].X.y +
	  G.Node[i][j][k].X.z*sin_angle;
	X.y = 0*G.Node[i][j][k].X.x +
	  G.Node[i][j][k].X.y +
	  0*G.Node[i][j][k].X.z;
	X.z = -G.Node[i][j][k].X.x*sin_angle +
	  0*G.Node[i][j][k].X.y +
	  G.Node[i][j][k].X.z*cos_angle; 
	G.Node[i][j][k].X = X;
      } /* endfor */
    } /* endfor */
  } /*endfor*/
  
  for( int k = G.KCl-G.Nghost; k<= G.KCu + G.Nghost; ++k){
    for (int j = G.JCl-G.Nghost ; j <= G.JCu+G.Nghost ; ++j) {
      for (int i = G.ICl-G.Nghost ; i <= G.ICu+G.Nghost ; ++i) {
	G.Cell[i][j][k].Xc = G.centroid(i, j, k);
      } /* endfor */
    } /* endfor */
  } /*endfor */

  return (G);
}

/*************************************************************************
 * Grid3D_Hexa_Block -- Input-output operators.                          *
 *************************************************************************/
inline ostream &operator << (ostream &out_file, 
                             const Grid3D_Hexa_Block &G) {

  out_file << G.NNi << " " << G.INl << " " << G.INu << "\n";
  out_file << G.NNj << " " << G.JNl << " " << G.JNu << "\n";
  out_file << G.NNk << " " << G.KNl << " " << G.KNu << "\n";
  if (G.NNi == 0 || G.NNj == 0 || G.NNk == 0 ) return(out_file);
  out_file << G.NCi << " " << G.ICl << " " << G.ICu << "\n";
  out_file << G.NCj << " " << G.JCl << " " << G.JCu << "\n";
  out_file << G.NCk << " " << G.KCl << " " << G.KCu << "\n";
  for (int k = G.KNl-G.Nghost ; k <= G.KNu+G.Nghost; ++k ) {
    for (int j = G.JNl-G.Nghost ; j <= G.JNu+G.Nghost; ++j ) {
      for (int i = G.INl-G.Nghost ; i <= G.INu+G.Nghost; ++i ) {
	out_file << G.Node[i][j][k].X << "\n";
      } /* endfor */
    } /* endfor */
  }/*endfor*/
  return (out_file);
}

inline istream &operator >> (istream &in_file, 
                             Grid3D_Hexa_Block &G) {
  int i, j, k, ni, il, iu, nj, jl, ju, nk, kl, ku;
  int ng;
  in_file.setf(ios::skipws);
  in_file >> ni >> il >> iu; in_file >> nj >> jl >> ju; in_file >> nk >> kl >> ku;
  in_file >> ng;
  in_file.unsetf(ios::skipws);
  if (ni == 0 || nj == 0|| nk == 0) {
      if (G.Node != NULL) G.deallocate(); return(in_file);
  } /* endif */
  if (G.Node == NULL || G.Cell == NULL || G.NNi != ni || G.NNj != nj || G.NNk != nk) {
      if (G.Node != NULL) G.deallocate(); 
      G.allocate(ni - (2*ng)+1, nj - (2*ng)+1, nk - (2*ng+1), ng);
  } /* endif */
  in_file.setf(ios::skipws);
  in_file >> ni >> il >> iu; in_file >> nj >> jl >> ju; in_file >> nk >> kl >> ku;
  in_file.unsetf(ios::skipws);
  for ( k = G.KNl-G.Nghost ; k <= G.KNu+G.Nghost; ++k) {
    for ( j = G.JNl-G.Nghost ; j <= G.JNu+G.Nghost; ++j ) {
      for ( i = G.INl-G.Nghost ; i <= G.INu+G.Nghost; ++i ) {
	in_file >> G.Node[i][j][k].X;
      } /* endfor */
    } /* endfor */
  } /* endfor */
  for ( k = G.KNl-G.Nghost ; k <= G.KCu+G.Nghost ; ++k){
    for ( j = G.JCl-G.Nghost ; j <= G.JCu+G.Nghost ; ++j) {
      for ( i = G.ICl-G.Nghost ; i <= G.ICu+G.Nghost ; ++i) {
	G.Cell[i][j][k].I = i; G.Cell[i][j][k].J = j; G.Cell[i][j][k].K = k;
	G.Cell[i][j][k].Xc = G.centroid(i, j, k);
	G.Cell[i][j][k].V =  G.volume(i, j, k);
      } /* endfor */
    } /* endfor */
  }/* endfor */
  return (in_file);
}

/*************************************************************************
 * Grid3D_Hexa_Block -- External subroutines for single grid block.      *
 *************************************************************************/

extern void Update_Cells(Grid3D_Hexa_Block &Grid);

extern void Write_Hexa_Block(Grid3D_Hexa_Block &Grid,
	                     ostream &Out_File);

extern void Output_Tecplot(Grid3D_Hexa_Block &Grid,
                           const int Block_Number,
                           const int Output_Title,
 	                   ostream &Out_File);

extern void Output_Nodes_Tecplot(Grid3D_Hexa_Block &Grid,
                                 const int Block_Number,
                                 const int Output_Title,
 	                         ostream &Out_File);

extern void Output_Cells_Tecplot(Grid3D_Hexa_Block &Grid,
                                 const int Block_Number,
                                 const int Output_Title,
 	                         ostream &Out_File);

extern void Output_Gnuplot(Grid3D_Hexa_Block &Grid,
                           const int Block_Number,
                           const int Output_Title,
 	                   ostream &Out_File);

extern void Update_Exterior_Nodes(Grid3D_Hexa_Block &Grid);

extern void Rotate_Hexa_Block(Grid3D_Hexa_Block &Grid,
			      const double &Angle, 
			      const double &Angle1, 
			      const double &Angle2);

/*************************************************************************
 * Grid3D_Hexa_Block -- External subroutines for 3D array of grid blocks.*
 *************************************************************************/

extern Grid3D_Hexa_Block*** Allocate_Multi_Block_Grid(
                                         Grid3D_Hexa_Block ***Grid_ptr,
					 const int Number_of_Blocks_Idir,
		                         const int Number_of_Blocks_Jdir,
					 const int Number_of_Blocks_Kdir);

extern Grid3D_Hexa_Block*** Deallocate_Multi_Block_Grid(
                                         Grid3D_Hexa_Block ***Grid_ptr,
					 const int Number_of_Blocks_Idir,
		                         const int Number_of_Blocks_Jdir,
					 const int Number_of_Blocks_Kdir);

extern void Write_Multi_Block_Grid(Grid3D_Hexa_Block ***Grid_ptr,
				   const int Number_of_Blocks_Idir,
		                   const int Number_of_Blocks_Jdir,
				   const int Number_of_Blocks_Kdir,
                                   ostream &Out_File);

extern void Output_Tecplot(Grid3D_Hexa_Block ***Grid_ptr,
			   const int Number_of_Blocks_Idir,
		           const int Number_of_Blocks_Jdir,
			   const int Number_of_Blocks_Kdir,
 	                   ostream &Out_File);

extern void Output_Nodes_Tecplot(Grid3D_Hexa_Block ***Grid_ptr,
                                 const int Number_of_Blocks_Idir,
                                 const int Number_of_Blocks_Jdir,
				 const int Number_of_Blocks_Kdir,
 	                         ostream &Out_File);

extern void Output_Cells_Tecplot(Grid3D_Hexa_Block ***Grid_ptr,
                                 const int Number_of_Blocks_Idir,
                                 const int Number_of_Blocks_Jdir,
				 const int Number_of_Blocks_Kdir,
 	                         ostream &Out_File);

extern void Output_Gnuplot(Grid3D_Hexa_Block ***Grid_ptr,
			   const int Number_of_Blocks_Idir,
		           const int Number_of_Blocks_Jdir,
			   const int Number_of_Blocks_Kdir,
 	                   ostream &Out_File);

/* Include ICEMCFD header file. */

#ifndef _ICEMCFD_INCLUDED
#include "../ICEM/ICEMCFD.h"
#endif // _ICEMCFD_INCLUDED

#endif /* _GRID3D_HEXA_BLOCK_INCLUDED  */

/* Cell2D.h:  Header file defining 2D cell types. */

#ifndef _CELL2D_INCLUDED
#define _CELL2D_INCLUDED

/* Include math macro and 2D vector header files. */

#ifndef _VECTOR2D_INCLUDED
#include "Math/Vector2D.h"
#endif //_VECTOR2D_INCLUDED

#include "include/require.h"

/* Define the basic 2D cell class. */

/*********************************************************
 * Class: Cell2D                                         *
 *                                                       *
 * Member functions                                      *
 *     Xc       -- Return 2D vector containing location  *
 *                 of cell center (centroid).            *
 *     I        -- Return I index for cell.              *
 *     J        -- Return J index for cell.              *
 *     A        -- Return the area of the cell.          *
 *     setloc   -- Set cell center location (position).  *
 *     setindex -- Set cell (i,j) indices.               * 
 *                                                       *
 * Member operators                                      *
 *      C -- a cell                                      *
 *                                                       *
 * C = C;                                                *
 * C == C;                                               *
 * C != C;                                               *
 * cout << C; (output function)                          *
 * cin  >> C; (input function)                           *
 *                                                       *
 *********************************************************/
class Cell2D{
private:
public:
  Vector2D         Xc;   // Location of cell center.
  int             I,J;   // (i,j) indices for cell.
  double            A;   // Cell area.
  // Made public so can access them.
		      
  /* Creation, copy, and assignment constructors. */
  Cell2D(void) {
    Xc.x = ONE; Xc.y = ONE; I = 0; J =0; A = ONE;
  }

  Cell2D(const Cell2D &Cell) {
    Xc = Cell.Xc; I = Cell.I; J = Cell.J; A = Cell.A;
  }

  Cell2D(const Vector2D &V) {
    Xc = V;
  }

  Cell2D(const double &xx, const double &yy) {
    Xc.x = xx; Xc.y = yy;
  }
    
  /* Destructor. */
  // ~Cell2D(void);
  // Use automatically generated destructor.

  /* Set cell center location. */
  void setloc(void);
  void setloc(const Cell2D &Cell);
  void setloc(const Vector2D &V);
  void setloc(const double &xx, const double &yy);

  /* Set cell (i,j) indices. */
  void setindex(void);
  void setindex(const Cell2D &Cell);
  void setindex(const int II, const int JJ);

  /* Assignment operator. */
  // Cell2D operator = (const Cell2D &Cell);
  // Use automatically generated assignment operator.

  /* Relational operators. */
  friend int operator ==(const Cell2D &Cell1,
			 const Cell2D &Cell2);
  friend int operator !=(const Cell2D &Cell1,
			 const Cell2D &Cell2);
    
  /* Input-output operators. */
  friend ostream &operator << (ostream &out_file,
			       const Cell2D &Cell);
  friend istream &operator >> (istream &in_file,
			       Cell2D &Cell);
    
};

/******************************************************************
 * Cell2D::setloc -- Set cell center location.                    *
 ******************************************************************/
inline void Cell2D::setloc(void) {
  Xc.x = ONE; Xc.y = ONE;
}

inline void Cell2D::setloc(const Cell2D &Cell) {
  Xc = Cell.Xc;
}

inline void Cell2D::setloc(const Vector2D &V) {
  Xc = V;
}

inline void Cell2D::setloc(const double &xx, const double &yy) {
  Xc.x = xx; Xc.y = yy;
}

/******************************************************************
 * Cell2D::setindex -- Set cell (i,j) indices.                    *
 ******************************************************************/
inline void Cell2D::setindex(void) {
  I = 0; J = 0;
}

inline void Cell2D::setindex(const Cell2D &Cell) {
  I = Cell.I; J = Cell.J;
}

inline void Cell2D::setindex(const int II, const int JJ) {
  I = II; J = JJ;
}

/******************************************************************
 * Cell2D -- Relational operators.                                *
 ******************************************************************/
inline int operator ==(const Cell2D &Cell1,
		       const Cell2D &Cell2) {
  return (Cell1.Xc == Cell2.Xc && Cell1.A == Cell2.A);
}

inline int operator !=(const Cell2D &Cell1,
		       const Cell2D &Cell2) {
  return (Cell1.Xc != Cell2.Xc || Cell1.A != Cell2.A);
}

/******************************************************************
 * Cell2D -- Input-output operators.                              *
 ******************************************************************/
inline ostream &operator << (ostream &out_file,
			     const Cell2D &Cell) {
  out_file << " " << Cell.I << " " << Cell.J << Cell.Xc;
  out_file.setf(ios::scientific);
  out_file << " " << Cell.A;
  out_file.unsetf(ios::scientific);
  return (out_file);
}

inline istream &operator >> (istream &in_file,
			     Cell2D &Cell) {
  in_file.setf(ios::skipws);
  in_file >> Cell.I >> Cell.J; 
  in_file.unsetf(ios::skipws);
  in_file >> Cell.Xc;
  in_file.setf(ios::skipws);
  in_file >> Cell.A;
  in_file.unsetf(ios::skipws);
  return (in_file);
}

/* Define the basic 2D node class. */

/*********************************************************
 * Class: Node2D                                         *
 *                                                       *
 * Member functions                                      *
 *      X       -- Return 2D vector containing location  *
 *                 of node.                              *
 *      setloc  -- Set node location (position).         * 
 *                                                       *
 * Member operators                                      *
 *      N -- a node                                      *
 *                                                       *
 * N = N;                                                *
 * N == N;                                               *
 * N != N;                                               *
 * cout << N; (output function)                          *
 * cin  >> N; (input function)                           *
 *                                                       *
 *********************************************************/
class Node2D{
private:
public:
  Vector2D         X;   // Node location.
  // Made public so can access them.
		      
  /* Creation, copy, and assignment constructors. */
  Node2D(void) {
    X.x = ONE; X.y = ONE;
  }

  Node2D(const Node2D &Node) {
    X = Node.X;
  }

  Node2D(const Vector2D &V) {
    X = V;
  }

  Node2D(const double &xx, const double &yy) {
    X.x = xx; X.y = yy;
  }
    
  /* Destructor. */
  // ~Node2D(void);
  // Use automatically generated destructor.

  /* Set cell center location. */
  void setloc(void);
  void setloc(const Node2D &Node);
  void setloc(const Vector2D &V);
  void setloc(const double &xx, const double &yy);
  double & x(void) { return X.x;}
  double & y(void) { return X.y;}
  /* Assignment operator. */
  // Node2D operator = (const Node2D &Node);
  // Use automatically generated assignment operator.

  /* Relational operators. */
  friend int operator ==(const Node2D &Node1,
			 const Node2D &Node2);
  friend int operator !=(const Node2D &Node1,
			 const Node2D &Node2);
    
  /* Input-output operators. */
  friend ostream &operator << (ostream &out_file,
			       const Node2D &Node);
  friend istream &operator >> (istream &in_file,
			       Node2D &Node);
    
};

/******************************************************************
 * Node2D::setloc -- Set cell center location.  *
 ******************************************************************/
inline void Node2D::setloc(void) {
  X.x = ONE; X.y = ONE;
}

inline void Node2D::setloc(const Node2D &Node) {
  X = Node.X;
}

inline void Node2D::setloc(const Vector2D &V) {
  X = V;
}

inline void Node2D::setloc(const double &xx, const double &yy) {
  X.x = xx; X.y = yy;
}

/******************************************************************
 * Node2D -- Relational operators.              *
 ******************************************************************/
inline int operator ==(const Node2D &Node1,
		       const Node2D &Node2) {
  return (Node1.X == Node2.X);
}

inline int operator !=(const Node2D &Node1,
		       const Node2D &Node2) {
  return (Node1.X != Node2.X);
}

/******************************************************************
 * Node2D -- Input-output operators.                *
 ******************************************************************/
inline ostream &operator << (ostream &out_file,
			     const Node2D &Node) {
  out_file << Node.X;
  return (out_file);
}

inline istream &operator >> (istream &in_file,
			     Node2D &Node) {
  in_file >> Node.X; 
  return (in_file);
}

/* Define the cell class for a uniform 2D Cartesian mesh. */

/*********************************************************
 * Class: Cell2D_Cartesian                               *
 *                                                       *
 * Member functions                                      *
 *      xc      -- Return 2D vector containing location  *
 *                 of cell center.                       *
 *      dx      -- Return 2D vector of cell dimensions.  *
 *                                                       *
 *                                                       *
 *                 xnNW-----xfN----xnNE                  *
 *                   |              |                    *
 *                   |              |                    *
 *       Cell       xfW      xc    xfE                   *
 *                   |              |                    *
 *                   |              |                    *
 *                 xnSW-----xfS----xnSE                  *
 *                                                       *
 *                                                       *
 *      A       -- Return the area of the cell.          *
 *      xnNW    -- Return 2D vector containing location  *
 *                 of north-west node.                   *
 *      xnNE    -- Return 2D vector containing location  *
 *                 of north-east node.                   *
 *      xnSE    -- Return 2D vector containing location  *
 *                 of south-east node.                   *
 *      xnSW    -- Return 2D vector containing location  *
 *                 of south-west node.                   *
 *      xfN     -- Return 2D vector containing location  *
 *                 of center of north cell face.         *
 *      xfS     -- Return 2D vector containing location  *
 *                 of center of south cell face.         *
 *      xfE     -- Return 2D vector containing location  *
 *                 of center of east cell face.          *
 *      xfW     -- Return 2D vector containing location  *
 *                 of center of west cell face.          *
 *      setsize -- Set cell dimensions.                  *
 *      setloc  -- Set cell center location (position).  * 
 *                                                       *
 * Member operators                                      *
 *      C -- a cell                                      *
 *                                                       *
 * C = C;                                                *
 * C == C;                                               *
 * C != C;                                               *
 * cout << C; (output function)                          *
 * cin  >> C; (input function)                           *
 *                                                       *
 *********************************************************/
class Cell2D_Cartesian{
private:
public:
  Vector2D         xc;   // Location of cell center.
  static Vector2D  dx;   // Vector of cell lengths.
  // Made public so can access them.
		      
  /* Creation, copy, and assignment constructors. */
  Cell2D_Cartesian(void) {
    xc.x = ONE; xc.y = ONE;  
  }

  Cell2D_Cartesian(const Cell2D_Cartesian &Cell) {
    xc = Cell.xc;
  }

  Cell2D_Cartesian(const Vector2D &V) {
    xc = V;
  }

  Cell2D_Cartesian(const double &xx, const double &yy) {
    xc.x = xx; xc.y = yy;
  }
    
  /* Destructor. */
  // ~Cell2D_Cartesian(void);
  // Use automatically generated destructor.

  /* Cell area. */
  double A(void);
  double A(void) const;

  /* Node locations. */
  Vector2D xnNW(void);
  Vector2D xnNW(void) const;

  Vector2D xnNE(void);
  Vector2D xnNE(void) const;

  Vector2D xnSE(void);
  Vector2D xnSE(void) const;

  Vector2D xnSW(void);
  Vector2D xnSW(void) const;

  Vector2D NodeNW(void);
  Vector2D NodeNW(void) const;

  Vector2D NodeNE(void);
  Vector2D NodeNE(void) const;

  Vector2D NodeSW(void);
  Vector2D NodeSW(void) const;

  Vector2D NodeSE(void);
  Vector2D NodeSE(void) const;

  /* Face midpoints. */
  Vector2D xfN(void);
  Vector2D xfN(void) const;

  Vector2D xfS(void);
  Vector2D xfS(void) const;

  Vector2D xfE(void);
  Vector2D xfE(void) const;

  Vector2D xfW(void);
  Vector2D xfW(void) const;

  /* Set cell dimensions. */
  void setsize(void);
  void setsize(const Vector2D &xx);
  void setsize(const double &xx, const double &yy);

  /* Set cell center location. */
  void setloc(void);
  void setloc(const Cell2D_Cartesian &Cell);
  void setloc(const Vector2D &V);
  void setloc(const double &xx, const double &yy);

  /* Assignment operator. */
  Cell2D_Cartesian& operator = (const Cell2D_Cartesian &Cell){
    if (this == &Cell) return *this;
    xc = Cell.xc; 
    dx = Cell.dx;
    return *this;
  }

  /* Relational operators. */
  friend int operator ==(const Cell2D_Cartesian &Cell1,
			 const Cell2D_Cartesian &Cell2);
  friend int operator !=(const Cell2D_Cartesian &Cell1,
			 const Cell2D_Cartesian &Cell2);
    
  /* Input-output operators. */
  friend ostream &operator << (ostream &out_file,
			       const Cell2D_Cartesian &Cell);
  friend istream &operator >> (istream &in_file,
			       Cell2D_Cartesian &Cell);
    
};

/******************************************************************
 * Cell2D_Cartesian::A -- Return cell area.                       *
 ******************************************************************/
inline double Cell2D_Cartesian::A(void) {
  return (dx.x*dx.y);
}

inline double Cell2D_Cartesian::A(void) const {
  return (dx.x*dx.y);
}

/******************************************************************
 * Cell2D_Cartesian::xn?? -- Get cell nodes.                      *
 ******************************************************************/
inline Vector2D Cell2D_Cartesian::xnNW(void) {
  return (Vector2D(xc-HALF*dx.x*Vector2D_NX+HALF*dx.y*Vector2D_NY));
}

inline Vector2D Cell2D_Cartesian::xnNW(void) const {
  return (Vector2D(xc-HALF*dx.x*Vector2D_NX+HALF*dx.y*Vector2D_NY));
}

inline Vector2D Cell2D_Cartesian::xnNE(void) {
  return (Vector2D(xc+HALF*dx.x*Vector2D_NX+HALF*dx.y*Vector2D_NY));
}

inline Vector2D Cell2D_Cartesian::xnNE(void) const {
  return (Vector2D(xc+HALF*dx.x*Vector2D_NX+HALF*dx.y*Vector2D_NY));
}

inline Vector2D Cell2D_Cartesian::xnSE(void) {
  return (Vector2D(xc+HALF*dx.x*Vector2D_NX-HALF*dx.y*Vector2D_NY));
}

inline Vector2D Cell2D_Cartesian::xnSE(void) const {
  return (Vector2D(xc+HALF*dx.x*Vector2D_NX-HALF*dx.y*Vector2D_NY));
}

inline Vector2D Cell2D_Cartesian::xnSW(void) {
  return (Vector2D(xc-HALF*dx.x*Vector2D_NX-HALF*dx.y*Vector2D_NY));
}

inline Vector2D Cell2D_Cartesian::xnSW(void) const {
  return (Vector2D(xc-HALF*dx.x*Vector2D_NX-HALF*dx.y*Vector2D_NY));
}

inline Vector2D Cell2D_Cartesian::NodeNW(void){
  return xnNW();
}

inline Vector2D Cell2D_Cartesian::NodeNW(void) const {
  return xnNW();
}

inline Vector2D Cell2D_Cartesian::NodeNE(void){
  return xnNE();
}
inline Vector2D Cell2D_Cartesian::NodeNE(void) const {
  return xnNE();
}

inline Vector2D Cell2D_Cartesian::NodeSW(void){
  return xnSW();
}

inline Vector2D Cell2D_Cartesian::NodeSW(void) const {
  return xnSW();
}

inline Vector2D Cell2D_Cartesian::NodeSE(void){
  return xnSE();
}

inline Vector2D Cell2D_Cartesian::NodeSE(void) const {
  return xnSE();
}

/******************************************************************
 * Cell2D_Cartesian::xf? -- Get cell faces.                       *
 ******************************************************************/
inline Vector2D Cell2D_Cartesian::xfN(void) {
  return (Vector2D(xc+HALF*dx.y*Vector2D_NY));
}

inline Vector2D Cell2D_Cartesian::xfN(void) const {
  return (Vector2D(xc+HALF*dx.y*Vector2D_NY));
}

inline Vector2D Cell2D_Cartesian::xfS(void) {
  return (Vector2D(xc-HALF*dx.y*Vector2D_NY));
}

inline Vector2D Cell2D_Cartesian::xfS(void) const {
  return (Vector2D(xc-HALF*dx.y*Vector2D_NY));
}

inline Vector2D Cell2D_Cartesian::xfE(void) {
  return (Vector2D(xc+HALF*dx.x*Vector2D_NX));
}

inline Vector2D Cell2D_Cartesian::xfE(void) const {
  return (Vector2D(xc+HALF*dx.x*Vector2D_NX));
}

inline Vector2D Cell2D_Cartesian::xfW(void) {
  return (Vector2D(xc-HALF*dx.x*Vector2D_NX));
}

inline Vector2D Cell2D_Cartesian::xfW(void) const {
  return (Vector2D(xc-HALF*dx.x*Vector2D_NX));
}

/******************************************************************
 * Cell2D_Cartesian::setsize -- Set cell dimensions.              *
 ******************************************************************/
inline void Cell2D_Cartesian::setsize(void) {
  dx.x = ONE; dx.y = ONE;
}

inline void Cell2D_Cartesian::setsize(const Vector2D &xx) {
  dx = xx;
}

inline void Cell2D_Cartesian::setsize(const double &xx, const double &yy) {
  dx.x = xx; dx.y = yy;
}

/******************************************************************
 * Cell2D_Cartesian::setloc -- Set cell center location.          *
 ******************************************************************/
inline void Cell2D_Cartesian::setloc(void) {
  xc.x = ONE; xc.y = ONE;
}

inline void Cell2D_Cartesian::setloc(const Cell2D_Cartesian &Cell) {
  xc = Cell.xc;
}

inline void Cell2D_Cartesian::setloc(const Vector2D  &V) {
  xc = V;
}

inline void Cell2D_Cartesian::setloc(const double &xx, const double &yy) {
  xc.x = xx; xc.y = yy;
}

/******************************************************************
 * Cell2D_Cartesian -- Relational operators.                      *
 ******************************************************************/
inline int operator ==(const Cell2D_Cartesian &Cell1,
		       const Cell2D_Cartesian &Cell2) {
  return (Cell1.xc == Cell2.xc);
}

inline int operator !=(const Cell2D_Cartesian &Cell1,
		       const Cell2D_Cartesian &Cell2) {
  return (Cell1.xc != Cell2.xc);
}

/******************************************************************
 * Cell2D_Cartesian -- Input-output operators.                    *
 ******************************************************************/
inline ostream &operator << (ostream &out_file,
			     const Cell2D_Cartesian &Cell) {
  out_file << Cell.xc;
  return (out_file);
}

inline istream &operator >> (istream &in_file,
			     Cell2D_Cartesian &Cell) {
  in_file >> Cell.xc; 
  return (in_file);
}


/* Define the cell class for a 2D quadrilateral mesh. */

/*********************************************************
 * Class: Cell2D_Quad                                    *
 *                                                       *
 * Member functions                                      *
 *      xc      -- Return 2D vector containing location  *
 *                 of cell center.                       *
 *      xnNW    -- Return 2D vector containing location  *
 *                 of north-west node.                   *
 *      xnNE    -- Return 2D vector containing location  *
 *                 of north-east node.                   *
 *      xnSE    -- Return 2D vector containing location  *
 *                 of south-east node.                   *
 *      xnSW    -- Return 2D vector containing location  *
 *                 of south-west node.                   *
 *                                                       *
 *                                                       *
 *                 xnNW-----xfN----xnNE                  *
 *                   |              |                    *
 *                   |              |                    *
 *       Cell       xfW      xc    xfE                   *
 *                   |              |                    *
 *                   |              |                    *
 *                 xnSW-----xfS----xnSE                  *
 *                                                       *
 *                                                       *
 *      A       -- Return the area of the cell.          *
 *      xfN     -- Return 2D vector containing location  *
 *                 of midpoint of north cell face.       *
 *      xfS     -- Return 2D vector containing location  *
 *                 of midpoint of south cell face.       *
 *      xfE     -- Return 2D vector containing location  *
 *                 of midpoint of east cell face.        *
 *      xfW     -- Return 2D vector containing location  *
 *                 of midpoint of west cell face.        *
 *      lfN     -- Return the length of north cell face. *
 *      lfS     -- Return the length of south cell face. *
 *      lfE     -- Return the length of east cell face.  *
 *      lfW     -- Return the length of west cell face.  *
 *      setnodes -- Set cell node locations.             * 
 *                                                       *
 * Member operators                                      *
 *      C -- a cell                                      *
 *                                                       *
 * C = C;                                                *
 * C == C;                                               *
 * C != C;                                               *
 * cout << C; (output function)                          *
 * cin  >> C; (input function)                           *
 *                                                       *
 *********************************************************/
class Cell2D_Quad{
private:
  static const double Gauss2QuadPoints_CoeffPoint1;
  static const double Gauss2QuadPoints_CoeffPoint2;

public:
  Vector2D         xc;   // Location of cell center.
  Vector2D       xnNW;   // Location of north-west node.
  Vector2D       xnNE;   // Location of north-east node.
  Vector2D       xnSE;   // Location of south-east node.
  Vector2D       xnSW;   // Location of south-west node.
  // Made public so can access them.
		      
  /* Creation, copy, and assignment constructors. */
  Cell2D_Quad(void);
  Cell2D_Quad(const Cell2D_Quad &Cell);
  Cell2D_Quad(const Vector2D &V1, const Vector2D &V2,
	      const Vector2D &V3, const Vector2D &V4);
  Cell2D_Quad(const double &x1, const double &y1,
	      const double &x2, const double &y2,
	      const double &x3, const double &y3,
	      const double &x4, const double &y4);
    
  /* Destructor. */
  // ~Cell2D_Quad(void);
  // Use automatically generated destructor.

  /* Cell area. */
  double A(void);
  double A(void) const;

  /* Cell centroid. */
  void centroid(void);

  /* Face midpoints. */
  Vector2D xfN(void);
  Vector2D xfN(void) const;

  Vector2D xfS(void);
  Vector2D xfS(void) const;

  Vector2D xfE(void);
  Vector2D xfE(void) const;

  Vector2D xfW(void);
  Vector2D xfW(void) const;

  /* Node location */
  Vector2D NodeNW(void){return xnNW;}
  Vector2D NodeNW(void) const {return xnNW;}

  Vector2D NodeNE(void){return xnNE;}
  Vector2D NodeNE(void) const {return xnNE;}

  Vector2D NodeSW(void){return xnSW;}
  Vector2D NodeSW(void) const {return xnSW;}

  Vector2D NodeSE(void){return xnSE;}
  Vector2D NodeSE(void) const {return xnSE;}

  /* Get cell face Gauss Quadrature Points */
  /* -> 2 points per face */
  void GaussQuadPointsFaceN(Vector2D & Point1, Vector2D & Point2);
  void GaussQuadPointsFaceS(Vector2D & Point1, Vector2D & Point2);
  void GaussQuadPointsFaceE(Vector2D & Point1, Vector2D & Point2);
  void GaussQuadPointsFaceW(Vector2D & Point1, Vector2D & Point2);

  /* Face lengths. */
  double lfN(void);
  double lfN(void) const;

  double lfS(void);
  double lfS(void) const;

  double lfE(void);
  double lfE(void) const;

  double lfW(void);
  double lfW(void) const;

  /* Set cell node locations. */
  void setnodes(void);
  void setnodes(const Cell2D_Quad &Cell);
  void setnodes(const Vector2D &V1, const Vector2D &V2,
		const Vector2D &V3, const Vector2D &V4);
  void setnodes(const double &x1, const double &y1,
		const double &x2, const double &y2,
		const double &x3, const double &y3,
		const double &x4, const double &y4);

  /* Assignment operator. */
  // Cell2D_Quad operator = (const Cell2D_Quad &Cell);
  // Use automatically generated assignment operator.

  /* Relational operators. */
  friend int operator ==(const Cell2D_Quad &Cell1,
			 const Cell2D_Quad &Cell2);
  friend int operator !=(const Cell2D_Quad &Cell1,
			 const Cell2D_Quad &Cell2);
    
  /* Input-output operators. */
  friend ostream &operator << (ostream &out_file,
			       const Cell2D_Quad &Cell);
  friend istream &operator >> (istream &in_file,
			       Cell2D_Quad &Cell);
    
};

/* Creation, copy, and assignment constructors. */
inline Cell2D_Quad::Cell2D_Quad(void) {
  xnNW.x = ZERO; xnNW.y = ONE; xnNE.x = ONE; xnNE.y = ONE;
  xnSE.x = ONE; xnSE.y = ZERO; xnSW.x = ZERO; xnSW.y = ZERO;
  centroid();
}

inline Cell2D_Quad::Cell2D_Quad(const Cell2D_Quad &Cell) {
  xnNW = Cell.xnNW; xnNE = Cell.xnNE; xnSE = Cell.xnSE; 
  xnSW = Cell.xnSW; xc = Cell.xc;
}

inline Cell2D_Quad::Cell2D_Quad(const Vector2D &V1, const Vector2D &V2,
			 const Vector2D &V3, const Vector2D &V4) {
  xnNW = V4; xnNE = V3; xnSE = V2; xnSW = V1; 
  centroid();
}

inline Cell2D_Quad::Cell2D_Quad(const double &x1, const double &y1,
				const double &x2, const double &y2,
				const double &x3, const double &y3,
				const double &x4, const double &y4) {
  xnNW.x = x4; xnNW.y = y4; xnNE.x = x3; xnNE.y = y3;
  xnSE.x = x2; xnSE.y = y2; xnSW.x = x1; xnSW.y = y1;
  centroid();
}

/******************************************************************
 * Cell2D_Quad::A -- Return cell area.                            *
 ******************************************************************/
inline double Cell2D_Quad::A(void) {
  return (HALF*(((xnSE-xnSW)^(xnNW-xnSW))+((xnNE-xnNW)^(xnNE-xnSE))));
}

inline double Cell2D_Quad::A(void) const {
  return (HALF*(((xnSE-xnSW)^(xnNW-xnSW))+((xnNE-xnNW)^(xnNE-xnSE))));
}

/******************************************************************
 * Cell2D_Quad::xf? -- Get face midpoints.                        *
 ******************************************************************/
inline Vector2D Cell2D_Quad::xfN(void) {
  return (HALF*(xnNE+xnNW));
}

inline Vector2D Cell2D_Quad::xfN(void) const {
  return (HALF*(xnNE+xnNW));
}

inline Vector2D Cell2D_Quad::xfS(void) {
  return (HALF*(xnSE+xnSW));
}

inline Vector2D Cell2D_Quad::xfS(void) const {
  return (HALF*(xnSE+xnSW));
}

inline Vector2D Cell2D_Quad::xfE(void) {
  return (HALF*(xnSE+xnNE));
}

inline Vector2D Cell2D_Quad::xfE(void) const {
  return (HALF*(xnSE+xnNE));
}

inline Vector2D Cell2D_Quad::xfW(void) {
  return (HALF*(xnSW+xnNW));
}

inline Vector2D Cell2D_Quad::xfW(void) const {
  return (HALF*(xnSW+xnNW));
}

/*************************************************************************
 * Cell2D_Quad::xface? -- Get cell face Gauss Quadrature Points.    *
 *************************************************************************/

/* -> 2 points per face */
inline void Cell2D_Quad::GaussQuadPointsFaceN( Vector2D & Point1, Vector2D & Point2){
  Point1 = Point2 = xnNW-xnNE;

  /* final value Point1 */
  Point1 = xnNE + Gauss2QuadPoints_CoeffPoint1*Point1;

  /* final value Point2 */
  Point2 = xnNE + Gauss2QuadPoints_CoeffPoint2*Point2;
}
inline void Cell2D_Quad::GaussQuadPointsFaceS( Vector2D & Point1, Vector2D & Point2){
  Point1 = Point2 = xnSE-xnSW;

  /* final value Point1 */
  Point1 = xnSW + Gauss2QuadPoints_CoeffPoint1*Point1;

  /* final value Point2 */
  Point2 = xnSW + Gauss2QuadPoints_CoeffPoint2*Point2;
}
inline void Cell2D_Quad::GaussQuadPointsFaceE( Vector2D & Point1, Vector2D & Point2){
  Point1 = Point2 = xnNE-xnSE;

  /* final value Point1 */
  Point1 = xnSE + Gauss2QuadPoints_CoeffPoint1*Point1;

  /* final value Point2 */
  Point2 = xnSE + Gauss2QuadPoints_CoeffPoint2*Point2;
}
inline void Cell2D_Quad::GaussQuadPointsFaceW( Vector2D & Point1, Vector2D & Point2){
  Point1 = Point2 = xnSW-xnNW;

  /* final value Point1 */
  Point1 = xnNW + Gauss2QuadPoints_CoeffPoint1*Point1;

  /* final value Point2 */
  Point2 = xnNW + Gauss2QuadPoints_CoeffPoint2*Point2;
}

/******************************************************************
 * Cell2D_Quad::lf? -- Get face lengths.                          *
 ******************************************************************/
inline double Cell2D_Quad::lfN(void) {
  return (abs(xnNE-xnNW));
}

inline double Cell2D_Quad::lfN(void) const {
  return (abs(xnNE-xnNW));
}

inline double Cell2D_Quad::lfS(void) {
  return (abs(xnSW-xnSE));
}

inline double Cell2D_Quad::lfS(void) const {
  return (abs(xnSW-xnSE));
}

inline double Cell2D_Quad::lfE(void) {
  return (abs(xnNE-xnSE));
}

inline double Cell2D_Quad::lfE(void) const {
  return (abs(xnNE-xnSE));
}

inline double Cell2D_Quad::lfW(void) {
  return (abs(xnNW-xnSW));
}

inline double Cell2D_Quad::lfW(void) const {
  return (abs(xnNW-xnSW));
}

/******************************************************************
 * Cell2D_Quad::setnodes -- Set cell node locations.              *
 ******************************************************************/
inline void Cell2D_Quad::setnodes(void) {
  xnNW.x = ZERO; xnNW.y = ONE; xnNE.x = ONE; xnNE.y = ONE;
  xnSE.x = ONE; xnSE.y = ZERO; xnSW.x = ZERO; xnSW.y = ZERO;
  centroid();
}

inline void Cell2D_Quad::setnodes(const Cell2D_Quad &Cell) {
  xnNW = Cell.xnNW; xnNE = Cell.xnNE; xnSE = Cell.xnSE; 
  xnSW = Cell.xnSW; xc = Cell.xc;
}

inline void Cell2D_Quad::setnodes(const Vector2D &V1, const Vector2D &V2,
                                  const Vector2D &V3, const Vector2D &V4) {
  xnNW = V4; xnNE = V3; xnSE = V2; xnSW = V1; 
  centroid();
}

inline void Cell2D_Quad::setnodes(const double &x1, const double &y1,
                                  const double &x2, const double &y2,
                                  const double &x3, const double &y3,
                                  const double &x4, const double &y4) {
  xnNW.x = x4; xnNW.y = y4; xnNE.x = x3; xnNE.y = y3;
  xnSE.x = x2; xnSE.y = y2; xnSW.x = x1; xnSW.y = y1;
  centroid();
}

inline void Cell2D_Quad::centroid(void){

  // Fourth relation
  Vector2D X1, X2, X3, X4, Xc1, Xc2;
  double A1, A2;
  // Cell nodes in counter-clockwise order.
  X1 = xnSW;
  X2 = xnSE;
  X3 = xnNE;
  X4 = xnNW;
  // Determine the centroid and area of the sub-triangles.
  Xc1 = (X1+X2+X3)/3.0;
  Xc2 = (X1+X3+X4)/3.0;
  A1 = ((X1^X2) + (X2^X3) + (X3^X1))/2.0;
  A2 = ((X1^X3) + (X3^X4) + (X4^X1))/2.0;
  // Return the area-weighted average of the centroids of the sub-triangles:
  xc = (A1*Xc1 + A2*Xc2)/(A1+A2);
}

/******************************************************************
 * Cell2D_Quad -- Relational operators.                           *
 ******************************************************************/
inline int operator ==(const Cell2D_Quad &Cell1,
		       const Cell2D_Quad &Cell2) {
  return (Cell1.xnNE == Cell2.xnNE && Cell1.xnNW == Cell2.xnNW &&
	  Cell1.xnSW == Cell2.xnSW && Cell1.xnSE == Cell2.xnSE);
}

inline int operator !=(const Cell2D_Quad &Cell1,
		       const Cell2D_Quad &Cell2) {
  return (Cell1.xnNE != Cell2.xnNE || Cell1.xnNW != Cell2.xnNW ||
	  Cell1.xnSW != Cell2.xnSW || Cell1.xnSE != Cell2.xnSE);
}

/******************************************************************
 * Cell2D_Quad -- Input-output operators.                         *
 ******************************************************************/
inline ostream &operator << (ostream &out_file,
			     const Cell2D_Quad &Cell) {
  out_file << Cell.xc;
  return (out_file);
}

inline istream &operator >> (istream &in_file,
			     Cell2D_Quad &Cell) {
  in_file >> Cell.xc; 
  return (in_file);
}

/******************************************************************
 * Useful 2D cell constants.                                      *
 ******************************************************************/
const Cell2D_Cartesian Cell2D_Cartesian_ONE(ONE, ONE);

#endif /* _CELL2D_INCLUDED  */

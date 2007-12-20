/* \file HO_Cell2D.h
   \brief Header file defining 2D high-order cell types. */

#ifndef _HO_CELL2D_INCLUDED
#define _HO_CELL2D_INCLUDED

/* Include required C++ libraries. */
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>

/* Using std namespace functions */
using namespace std;

/* Include CFFC header files */
#include "../Math/Math.h"	/* Include math macro header files. */
#include "../Math/Vector2D.h"   /* Include vector header files. */


/* Define the classes. */

/*!
 * \class Cell2D_HO
 *
 * \begin{verbatim}
 * Member functions
 *     Xc       -- Return 2D vector containing location of the cell
 *                 center (centroid).
 *     I        -- Return I index for cell.
 *     J        -- Return J index for cell.
 *     A        -- Return the area of the cell.
 *     setloc   -- Set cell center location (position).
 *     setindex -- Set cell (i,j) indices.
 *
 * Member operators
 *      C -- a cell
 *
 * C = C;
 * C == C;
 * C != C;
 * cout << C; (output function)
 * cin  >> C; (input function)
 * \end{verbatim}
 */
class Cell2D_HO{
public:
  Vector2D     Xc;    //!< Location of cell center.
  int        I, J;    //!< (i,j) indices for cell.
  double        A;    //!< Cell area.
		  

  //! @name Creation and copy constructors.
  //@{
  //! Default constructor.
  Cell2D_HO(void);

  //! Copy constructor
  Cell2D_HO(const Cell2D_HO &Cell);

  //! Constructor with centroid location
  Cell2D_HO(const Vector2D &V);

  //! Constructor with centroid locations
  Cell2D_HO(const double &xx, const double &yy);
  //@}
    
  //! @name Set cell center location.
  //@{
  void setloc(void);
  void setloc(const Cell2D_HO &Cell);
  void setloc(const Vector2D &V);
  void setloc(const double &xx, const double &yy);
  //@}

  //! @name Set cell (i,j) indices.
  //@{
  void setindex(void);
  void setindex(const Cell2D_HO &Cell);
  void setindex(const int II, const int JJ);
  //@}

  //! @name Relational operators.
  //@{
  friend bool operator ==(const Cell2D_HO &Cell1, const Cell2D_HO &Cell2);
  friend bool operator !=(const Cell2D_HO &Cell1, const Cell2D_HO &Cell2);
  //@}

  //! @name Input-output operators.
  //@{
  friend ostream &operator << (ostream &out_file, const Cell2D_HO &Cell);
  friend istream &operator >> (istream &in_file, Cell2D_HO &Cell);
  //@}

private:    
};


/*!
 * Default constructor
 */
inline Cell2D_HO::Cell2D_HO(void): Xc(ONE), I(0), J(0), A(ONE)
{
  // 
}

/*!
 * Copy constructor
 */
inline Cell2D_HO::Cell2D_HO(const Cell2D_HO &Cell):
  Xc(Cell.Xc), I(Cell.I), J(Cell.J), A(Cell.A)
{
  // 
}

/*!
 * Constructor with centroid location
 */
inline Cell2D_HO::Cell2D_HO(const Vector2D &V): Xc(V), I(0), J(0), A(ONE)
{
  // 
}

/*!
 * Constructor with centroid locations
 */
inline Cell2D_HO::Cell2D_HO(const double &xx, const double &yy):
  I(0), J(0), A(ONE)
{
  Xc.x = xx; 
  Xc.y = yy;
}

/*!
 * Set cell center location
 */
inline void Cell2D_HO::setloc(void) {
  Xc.x = ONE; Xc.y = ONE;
}

/*!
 * Set cell center location
 */
inline void Cell2D_HO::setloc(const Cell2D_HO &Cell) {
  Xc = Cell.Xc;
}

/*!
 * Set cell center location
 */
inline void Cell2D_HO::setloc(const Vector2D &V) {
  Xc = V;
}

/*!
 * Set cell center location
 */
inline void Cell2D_HO::setloc(const double &xx, const double &yy) {
  Xc.x = xx; Xc.y = yy;
}

/*
 * Set cell (i,j) indices.
 */
inline void Cell2D_HO::setindex(void) {
  I = 0; J = 0;
}

/*
 * Set cell (i,j) indices.
 */
inline void Cell2D_HO::setindex(const Cell2D_HO &Cell) {
    I = Cell.I; J = Cell.J;
}

/*
 * Set cell (i,j) indices.
 */
inline void Cell2D_HO::setindex(const int II, const int JJ) {
    I = II; J = JJ;
}

/*
 * Equal operator
 */
inline bool operator ==(const Cell2D_HO &Cell1,
			const Cell2D_HO &Cell2) {
  return (Cell1.Xc == Cell2.Xc && Cell1.A == Cell2.A);
}

/*
 * Not-equal operator
 */
inline bool operator !=(const Cell2D_HO &Cell1,
			const Cell2D_HO &Cell2) {
  return (Cell1.Xc != Cell2.Xc || Cell1.A != Cell2.A);
}

/*!
 * Output operator
 */
inline ostream &operator << (ostream &out_file,
			     const Cell2D_HO &Cell) {
  out_file << " " << Cell.I << " " << Cell.J << Cell.Xc;
  out_file.setf(ios::scientific);
  out_file << " " << Cell.A;
  out_file.unsetf(ios::scientific);
  return (out_file);
}

/*!
 * Input operator
 */
inline istream &operator >> (istream &in_file,
			     Cell2D_HO &Cell) {
  in_file.setf(ios::skipws);
  in_file >> Cell.I >> Cell.J; 
  in_file.unsetf(ios::skipws);
  in_file >> Cell.Xc;
  in_file.setf(ios::skipws);
  in_file >> Cell.A;
  in_file.unsetf(ios::skipws);
  return (in_file);
}

/* Define the cell class for a uniform 2D Cartesian mesh. */

/*********************************************************
 * Class: Cell2D_Cartesian_HO                               *
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
class Cell2D_Cartesian_HO{
  private:
  public:
  Vector2D         xc;   // Location of cell center.
  static Vector2D  dx;   // Vector of cell lengths.

  static Cell2D_Cartesian_HO Cell2D_Cartesian_HO_ONE;
		      
    /* Creation, copy, and assignment constructors. */
    Cell2D_Cartesian_HO(void) {
       xc.x = ONE; xc.y = ONE;  
    }

    Cell2D_Cartesian_HO(const Cell2D_Cartesian_HO &Cell) {
       xc = Cell.xc;
    }

    Cell2D_Cartesian_HO(const Vector2D &V) {
       xc = V;
    }

    Cell2D_Cartesian_HO(const double &xx, const double &yy) {
       xc.x = xx; xc.y = yy;
    }
    
    /* Destructor. */
    // ~Cell2D_Cartesian_HO(void);
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
    void setloc(const Cell2D_Cartesian_HO &Cell);
    void setloc(const Vector2D &V);
    void setloc(const double &xx, const double &yy);

    /* Assignment operator. */
    // Cell2D_Cartesian_HO operator = (const Cell2D_Cartesian_HO &Cell);
    // Use automatically generated assignment operator.

    /* Relational operators. */
    friend int operator ==(const Cell2D_Cartesian_HO &Cell1,
			   const Cell2D_Cartesian_HO &Cell2);
    friend int operator !=(const Cell2D_Cartesian_HO &Cell1,
			   const Cell2D_Cartesian_HO &Cell2);
    
    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file,
				 const Cell2D_Cartesian_HO &Cell);
    friend istream &operator >> (istream &in_file,
				 Cell2D_Cartesian_HO &Cell);
    
};

/******************************************************************
 * Cell2D_Cartesian_HO::A -- Return cell area.                       *
 ******************************************************************/
inline double Cell2D_Cartesian_HO::A(void) {
    return (dx.x*dx.y);
}

inline double Cell2D_Cartesian_HO::A(void) const {
    return (dx.x*dx.y);
}

/******************************************************************
 * Cell2D_Cartesian_HO::xn?? -- Get cell nodes.                      *
 ******************************************************************/
inline Vector2D Cell2D_Cartesian_HO::xnNW(void) {
    return (Vector2D(xc-HALF*dx.x*Vector2D_NX+HALF*dx.y*Vector2D_NY));
}

inline Vector2D Cell2D_Cartesian_HO::xnNW(void) const {
    return (Vector2D(xc-HALF*dx.x*Vector2D_NX+HALF*dx.y*Vector2D_NY));
}

inline Vector2D Cell2D_Cartesian_HO::xnNE(void) {
    return (Vector2D(xc+HALF*dx.x*Vector2D_NX+HALF*dx.y*Vector2D_NY));
}

inline Vector2D Cell2D_Cartesian_HO::xnNE(void) const {
    return (Vector2D(xc+HALF*dx.x*Vector2D_NX+HALF*dx.y*Vector2D_NY));
}

inline Vector2D Cell2D_Cartesian_HO::xnSE(void) {
    return (Vector2D(xc+HALF*dx.x*Vector2D_NX-HALF*dx.y*Vector2D_NY));
}

inline Vector2D Cell2D_Cartesian_HO::xnSE(void) const {
    return (Vector2D(xc+HALF*dx.x*Vector2D_NX-HALF*dx.y*Vector2D_NY));
}

inline Vector2D Cell2D_Cartesian_HO::xnSW(void) {
    return (Vector2D(xc-HALF*dx.x*Vector2D_NX-HALF*dx.y*Vector2D_NY));
}

inline Vector2D Cell2D_Cartesian_HO::xnSW(void) const {
    return (Vector2D(xc-HALF*dx.x*Vector2D_NX-HALF*dx.y*Vector2D_NY));
}

/******************************************************************
 * Cell2D_Cartesian_HO::xf? -- Get cell faces.                       *
 ******************************************************************/
inline Vector2D Cell2D_Cartesian_HO::xfN(void) {
    return (Vector2D(xc+HALF*dx.y*Vector2D_NY));
}

inline Vector2D Cell2D_Cartesian_HO::xfN(void) const {
    return (Vector2D(xc+HALF*dx.y*Vector2D_NY));
}

inline Vector2D Cell2D_Cartesian_HO::xfS(void) {
    return (Vector2D(xc-HALF*dx.y*Vector2D_NY));
}

inline Vector2D Cell2D_Cartesian_HO::xfS(void) const {
    return (Vector2D(xc-HALF*dx.y*Vector2D_NY));
}

inline Vector2D Cell2D_Cartesian_HO::xfE(void) {
    return (Vector2D(xc+HALF*dx.x*Vector2D_NX));
}

inline Vector2D Cell2D_Cartesian_HO::xfE(void) const {
    return (Vector2D(xc+HALF*dx.x*Vector2D_NX));
}

inline Vector2D Cell2D_Cartesian_HO::xfW(void) {
    return (Vector2D(xc-HALF*dx.x*Vector2D_NX));
}

inline Vector2D Cell2D_Cartesian_HO::xfW(void) const {
    return (Vector2D(xc-HALF*dx.x*Vector2D_NX));
}

/******************************************************************
 * Cell2D_Cartesian_HO::setsize -- Set cell dimensions.              *
 ******************************************************************/
inline void Cell2D_Cartesian_HO::setsize(void) {
  dx.x = ONE; dx.y = ONE;
}

inline void Cell2D_Cartesian_HO::setsize(const Vector2D &xx) {
  dx = xx;
}

inline void Cell2D_Cartesian_HO::setsize(const double &xx, const double &yy) {
  dx.x = xx; dx.y = yy;
}

/******************************************************************
 * Cell2D_Cartesian_HO::setloc -- Set cell center location.          *
 ******************************************************************/
inline void Cell2D_Cartesian_HO::setloc(void) {
    xc.x = ONE; xc.y = ONE;
}

inline void Cell2D_Cartesian_HO::setloc(const Cell2D_Cartesian_HO &Cell) {
    xc = Cell.xc;
}

inline void Cell2D_Cartesian_HO::setloc(const Vector2D  &V) {
    xc = V;
}

inline void Cell2D_Cartesian_HO::setloc(const double &xx, const double &yy) {
    xc.x = xx; xc.y = yy;
}

/******************************************************************
 * Cell2D_Cartesian_HO -- Relational operators.                      *
 ******************************************************************/
inline int operator ==(const Cell2D_Cartesian_HO &Cell1,
		       const Cell2D_Cartesian_HO &Cell2) {
    return (Cell1.xc == Cell2.xc);
}

inline int operator !=(const Cell2D_Cartesian_HO &Cell1,
		       const Cell2D_Cartesian_HO &Cell2) {
    return (Cell1.xc != Cell2.xc);
}

/******************************************************************
 * Cell2D_Cartesian_HO -- Input-output operators.                    *
 ******************************************************************/
inline ostream &operator << (ostream &out_file,
			     const Cell2D_Cartesian_HO &Cell) {
    out_file << Cell.xc;
    return (out_file);
}

inline istream &operator >> (istream &in_file,
			     Cell2D_Cartesian_HO &Cell) {
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
  public:
    Vector2D         xc;   // Location of cell center.
    Vector2D       xnNW;   // Location of north-west node.
    Vector2D       xnNE;   // Location of north-east node.
    Vector2D       xnSE;   // Location of south-east node.
    Vector2D       xnSW;   // Location of south-west node.
	                   // Made public so can access them.
		      
    /* Creation, copy, and assignment constructors. */
    Cell2D_Quad(void) {
       xnNW.x = ZERO; xnNW.y = ONE; xnNE.x = ONE; xnNE.y = ONE;
       xnSE.x = ONE; xnSE.y = ZERO; xnSW.x = ZERO; xnSW.y = ZERO;
       xc = (xnNW+xnNE+xnSE+xnSW)/FOUR; 
    }

    Cell2D_Quad(const Cell2D_Quad &Cell) {
       xnNW = Cell.xnNW; xnNE = Cell.xnNE; xnSE = Cell.xnSE; 
       xnSW = Cell.xnSW; xc = Cell.xc;
    }

    Cell2D_Quad(const Vector2D &V1, const Vector2D &V2,
                         const Vector2D &V3, const Vector2D &V4) {
       xnNW = V4; xnNE = V3; xnSE = V2; xnSW = V1; 
       xc = (xnNW+xnNE+xnSE+xnSW)/FOUR;
    }

    Cell2D_Quad(const double &x1, const double &y1,
                         const double &x2, const double &y2,
                         const double &x3, const double &y3,
                         const double &x4, const double &y4) {
       xnNW.x = x4; xnNW.y = y4; xnNE.x = x3; xnNE.y = y3;
       xnSE.x = x2; xnSE.y = y2; xnSW.x = x1; xnSW.y = y1;
       xc = (xnNW+xnNE+xnSE+xnSW)/FOUR;
    }
    
    /* Destructor. */
    // ~Cell2D_Quad(void);
    // Use automatically generated destructor.

    /* Cell area. */
    double A(void);
    double A(void) const;

    /* Face midpoints. */
    Vector2D xfN(void);
    Vector2D xfN(void) const;

    Vector2D xfS(void);
    Vector2D xfS(void) const;

    Vector2D xfE(void);
    Vector2D xfE(void) const;

    Vector2D xfW(void);
    Vector2D xfW(void) const;

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
    xc = (xnNW+xnNE+xnSE+xnSW)/FOUR;
}

inline void Cell2D_Quad::setnodes(const Cell2D_Quad &Cell) {
    xnNW = Cell.xnNW; xnNE = Cell.xnNE; xnSE = Cell.xnSE; 
    xnSW = Cell.xnSW; xc = Cell.xc;
}

inline void Cell2D_Quad::setnodes(const Vector2D &V1, const Vector2D &V2,
                                  const Vector2D &V3, const Vector2D &V4) {
    xnNW = V4; xnNE = V3; xnSE = V2; xnSW = V1; 
    xc = (xnNW+xnNE+xnSE+xnSW)/FOUR;
}

inline void Cell2D_Quad::setnodes(const double &x1, const double &y1,
                                  const double &x2, const double &y2,
                                  const double &x3, const double &y3,
                                  const double &x4, const double &y4) {
    xnNW.x = x4; xnNW.y = y4; xnNE.x = x3; xnNE.y = y3;
    xnSE.x = x2; xnSE.y = y2; xnSW.x = x1; xnSW.y = y1;
    xc = (xnNW+xnNE+xnSE+xnSW)/FOUR;
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

#endif /* _CELL2D_INCLUDED  */

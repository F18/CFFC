/* Cell1D.h:  Header file defining 1D cell types. */

#ifndef _CELL1D_INCLUDED
#define _CELL1D_INCLUDED

/* Include required C++ libraries. */

#include <iostream>

using namespace std;

/* Include math macro header file. */
#include "../../../src_2D/Math/Math.h"

/* Define the classes. */

/********************************************************
 * Class: Cell1D_Uniform                                *
 *                                                      *
 * Member functions                                     *
 *      x       -- Return cell center.                  *
 *      dx      -- Return cell length.                  *
 *      setsize -- Set cell length.                     *
 *      setloc  -- Set cell center location (position). * 
 *                                                      *
 * Member operators                                     *
 *      C -- a cell                                     *
 *                                                      *
 * C = C;                                               *
 * C == C;                                              *
 * C != C;                                              *
 * cout << C; (output function)                         *
 * cin  >> C; (input function)                          *
 *                                                      *
 ********************************************************/
class Cell1D_Uniform{
  private:
  public:
    double          x;   // Cell center.
    static double  dx;   // Cell length.
	                 // Made public so can access them.
		      
    /* Creation, copy, and assignment constructors. */
    Cell1D_Uniform(void) {
       x = ONE;
    }

    Cell1D_Uniform(const Cell1D_Uniform &Cell) {
       x = Cell.x;
    }

    Cell1D_Uniform(const double &location) {
       x = location;
    }
    
    /* Destructor. */
    // ~Cell1D_Uniform(void);
    // Use automatically generated destructor.

    /* Set cell size. */
    void setsize(void);
    void setsize(const double &length);

    /* Set cell center location. */
    void setloc(void);
    void setloc(const Cell1D_Uniform &Cell);
    void setloc(const double &location);

    /* Assignment operator. */
    // Cell1D_Uniform operator = (const Cell1D_Uniform &Cell);
    // Use automatically generated assignment operator.

    /* Relational operators. */
    friend int operator ==(const Cell1D_Uniform &Cell1,
			   const Cell1D_Uniform &Cell2);
    friend int operator !=(const Cell1D_Uniform &Cell1,
			   const Cell1D_Uniform &Cell2);
    
    /* Input-output operators. */
    friend ostream & operator << (ostream & out_file,
				 const Cell1D_Uniform &Cell);
    friend istream & operator >> (istream & in_file,
				 Cell1D_Uniform &Cell);
    
};

/********************************************************
 * Cell1D_Uniform::setsize -- Set cell size.            *
 ********************************************************/
inline void Cell1D_Uniform::setsize(void) {
    dx = ONE;
}

inline void Cell1D_Uniform::setsize(const double &length) {
    dx = length;
}

/********************************************************
 * Cell1D_Uniform::setloc -- Set cell center location.  *
 ********************************************************/
inline void Cell1D_Uniform::setloc(void) {
    x = ONE;
}

inline void Cell1D_Uniform::setloc(const Cell1D_Uniform &Cell) {
    x = Cell.x;
}

inline void Cell1D_Uniform::setloc(const double &location) {
    x = location;
}

/********************************************************
 * Cell1D_Uniform -- Relational operators.              *
 ********************************************************/
inline int operator ==(const Cell1D_Uniform &Cell1,
		       const Cell1D_Uniform &Cell2) {
    return (Cell1.x == Cell2.x);
}

inline int operator !=(const Cell1D_Uniform &Cell1,
		       const Cell1D_Uniform &Cell2) {
    return (Cell1.x != Cell2.x);
}

/********************************************************
 * Cell1D_Uniform -- Input-output operators.            *
 ********************************************************/
inline ostream & operator << (ostream &out_file,
			     const Cell1D_Uniform &Cell) {
    out_file.setf(ios::scientific);
    out_file << " " << Cell.x << " " << Cell.dx;
    out_file.unsetf(ios::scientific);
    return (out_file);
}

inline istream & operator >> (istream &in_file,
			     Cell1D_Uniform &Cell) {
    in_file.setf(ios::skipws);
    in_file >> Cell.x;
    in_file.unsetf(ios::skipws);
    return (in_file);
}

/********************************************************
 * Class: Cell1D_NonUniform                             *
 *                                                      *
 * Member functions                                     *
 *      x       -- Return cell center.                  *
 *      dx      -- Return cell length.                  *
 *      setsize -- Set cell length.                     *
 *      setloc  -- Set cell center location (position). *
 *                                                      *
 * Member operators                                     *
 *      C -- a cell                                     *
 *      s -- a scalar (double)                          *
 *                                                      *
 * C = C;                                               *
 * C = C + s;                                           *
 * C = C - s;                                           *
 * C = s * C;                                           *
 * C = C * s;                                           *
 * C = C / s;                                           *
 * C = +C;                                              *
 * C = -C;                                              *
 * C += s;                                              *
 * C -= s;                                              *
 * C == C;                                              *
 * C != C;                                              *
 * cout << C; (output function)                         *
 * cin  >> C; (input function)                          *
 *                                                      *
 ********************************************************/
class Cell1D_NonUniform{
  private:
  public:
    double          x;   // Cell center.
    double         dx;   // Cell length.
	                 // Made public so can access them.
		      
    /* Creation, copy, and assignment constructors. */
    Cell1D_NonUniform(void) {
       x = ONE; dx = ONE;
    }

    Cell1D_NonUniform(const Cell1D_NonUniform &Cell) {
       x = Cell.x; dx = Cell.dx;
    }

    Cell1D_NonUniform(const double &location,
		      const double &length) {
       x = location; dx = length;
    }
    
    /* Destructor. */
    // ~Cell1D_NonUniform(void);
    // Use automatically generated destructor.

    /* Set cell size. */
    void setsize(void);
    void setsize(const Cell1D_NonUniform &Cell);
    void setsize(const double &length);

     /* Set cell center location. */
    void setloc(void);
    void setloc(const Cell1D_NonUniform &Cell);
    void setloc(const double &location);

    /* Assignment operator. */
    // Cell1D_NonUniform operator = (const Cell1D_NonUniform &Cell);
    // Use automatically generated assignment operator.

    /* Binary arithmetic operators. */
    friend Cell1D_NonUniform operator +(const Cell1D_NonUniform &Cell,
					const double &a);
    friend Cell1D_NonUniform operator -(const Cell1D_NonUniform &Cell,
					const double &a);
    friend Cell1D_NonUniform operator *(const Cell1D_NonUniform &Cell,
					const double &a);
    friend Cell1D_NonUniform operator *(const double &a,
					const Cell1D_NonUniform &Cell);
    friend Cell1D_NonUniform operator /(const Cell1D_NonUniform &Cell,
					const double &a);

    /* Unary arithmetic operators. */
    friend Cell1D_NonUniform operator +(const Cell1D_NonUniform &Cell);
    friend Cell1D_NonUniform operator -(const Cell1D_NonUniform &Cell);

    /* Shortcut arithmetic operators. */
    friend Cell1D_NonUniform &operator +=(Cell1D_NonUniform &Cell,
					  const double &a);
    friend Cell1D_NonUniform &operator -=(Cell1D_NonUniform &Cell,
					  const double &a);
    
    /* Relational operators. */
    friend int operator ==(const Cell1D_NonUniform &Cell1,
			   const Cell1D_NonUniform &Cell2);
    friend int operator !=(const Cell1D_NonUniform &Cell1,
			   const Cell1D_NonUniform &Cell2);
    
    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file,
				 const Cell1D_NonUniform &Cell);
    friend istream &operator >> (istream &in_file,
				 Cell1D_NonUniform &Cell);
    
};

/********************************************************
 * Cell1D_NonUniform::setsize -- Set cell size.         *
 ********************************************************/
inline void Cell1D_NonUniform::setsize(void) {
    dx = ONE;
}

inline void Cell1D_NonUniform::setsize(const Cell1D_NonUniform &Cell) {
    dx = Cell.dx;
}

inline void Cell1D_NonUniform::setsize(const double &length) {
    dx = length;
}

/*********************************************************
 * Cell1D_NonUniform::setloc -- Set cell center location.*
 *********************************************************/
inline void Cell1D_NonUniform::setloc(void) {
    x = ONE;
}

inline void Cell1D_NonUniform::setloc(const Cell1D_NonUniform &Cell) {
    x = Cell.x;
}

inline void Cell1D_NonUniform::setloc(const double &location) {
    x = location;
}

/********************************************************
 * Cell1D_NonUniform -- Binary arithmetic operators.    *
 ********************************************************/
inline Cell1D_NonUniform operator +(const Cell1D_NonUniform &Cell,
				    const double &a) {
    return (Cell1D_NonUniform(Cell.x + a, TWO*(a-HALF*Cell.dx)));
}

inline Cell1D_NonUniform operator -(const Cell1D_NonUniform &Cell,
				    const double &a) {
    return (Cell1D_NonUniform(Cell.x - a, TWO*(a-HALF*Cell.dx)));
}

inline Cell1D_NonUniform operator *(const Cell1D_NonUniform &Cell,
				    const double &a) {
    return (Cell1D_NonUniform(a*Cell.x,a*Cell.dx));
}

inline Cell1D_NonUniform operator *(const double &a,
				    const Cell1D_NonUniform &Cell) {
    return (Cell1D_NonUniform(a*Cell.x,a*Cell.dx));
}

inline Cell1D_NonUniform operator /(const Cell1D_NonUniform &Cell,
				    const double &a) {
    return (Cell1D_NonUniform(Cell.x,Cell.dx/a));
}

/********************************************************
 * Cell1D_NonUniform -- Unary arithmetic operators.     *
 ********************************************************/
inline Cell1D_NonUniform operator +(const Cell1D_NonUniform &Cell) {
    return (Cell1D_NonUniform(Cell.x,Cell.dx));
}

inline Cell1D_NonUniform operator -(const Cell1D_NonUniform &Cell) {
    return (Cell1D_NonUniform(-Cell.x,Cell.dx));
}

/********************************************************
 * Cell1D_NonUniform -- Shortcut arithmetic operators.  *
 ********************************************************/
inline Cell1D_NonUniform &operator +=(Cell1D_NonUniform &Cell,
				      const double &a) {
    Cell.x += a;
    Cell.dx = TWO*(a-HALF*Cell.dx);
     return (Cell);
}

inline Cell1D_NonUniform &operator -=(Cell1D_NonUniform &Cell,
				      const double &a) {
    Cell.x -= a;
    Cell.dx = TWO*(a-HALF*Cell.dx);
    return (Cell);
}

/********************************************************
 * Cell1D_NonUniform -- Relational operators.           *
 ********************************************************/
inline int operator ==(const Cell1D_NonUniform &Cell1,
		       const Cell1D_NonUniform &Cell2) {
    return (Cell1.x == Cell2.x && Cell1.dx == Cell2.dx);
}

inline int operator !=(const Cell1D_NonUniform &Cell1,
		       const Cell1D_NonUniform &Cell2) {
    return (Cell1.x != Cell2.x || Cell1.dx != Cell2.dx);
}

/********************************************************
 * Cell1D_NonUniform -- Input-output operators.         *
 ********************************************************/
inline ostream &operator << (ostream &out_file,
			     const Cell1D_NonUniform &Cell) {
    out_file.setf(ios::scientific);
    out_file << " " << Cell.x << " " << Cell.dx;
    out_file.unsetf(ios::scientific);
    return (out_file);
}

inline istream &operator >> (istream &in_file,
			     Cell1D_NonUniform &Cell) {
    in_file.setf(ios::skipws);
    in_file >> Cell.x >> Cell.dx;
    in_file.unsetf(ios::skipws);
    return (in_file);
}

/********************************************************
 * Useful 1D cell constants.                            *
 ********************************************************/
const Cell1D_Uniform Cell1D_Uniform_ONE(ONE);
const Cell1D_NonUniform Cell1D_NonUniform_ONE(ONE,ONE);

#endif /* _CELL1D_INCLUDED  */

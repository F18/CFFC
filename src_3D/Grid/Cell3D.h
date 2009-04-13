/*!\file Cell3D.h
  \brief Header file defining 3D high-order cell types. */

#ifndef _CELL3D_INCLUDED
#define _CELL3D_INCLUDED

/* Include required C++ libraries. */
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>

/* Using std namespace functions */
using namespace std;

/* Include CFFC header files */
/* Include math macro and 3D vector header files. */
#ifndef _MATH_MACROS_INCLUDED
#include "../Math/Math.h"
#endif // _MATH_MACROS_INCLUDED

#ifndef _VECTOR3D_INCLUDED
#include "../Math/Vector3D.h"
#endif //_VECTOR3D_INCLUDED

#ifndef _TAYLORDERIVATIVES_INCLUDED
#include "../HighOrderReconstruction/TDContainer.h"
#endif // _TAYLORDERIVATIVES_INCLUDED

#ifndef _NUMERICAL_LIBRARY_INCLUDED
#include "../Math/NumericalLibrary.h"
#endif // _NUMERICAL_LIBRARY_INCLUDED

/* Define the basic 3D cell class. */

/*!
 * \class Cell3D
 *
 * \verbatim
 * Member functions
 *     Xc       -- Return 3D vector containing location of the cell
 *                 center (centroid).
 *     I        -- Return I index for cell.
 *     J        -- Return J index for cell.
 *     K        -- Return K index for cell.
 *     V        -- Return the volume of the cell.
 *     setloc   -- Set cell center location (position).
 *     setindex -- Set cell (i,j,k) indices.
 *     GeomCoeff-- Volume integrals of geometric moments
 *
 * Member operators
 *      C -- a cell
 *
 * C = C;
 * C == C;
 * C != C;
 * cout << C; (output function)
 * cin  >> C; (input function)
 * \endverbatim
 */
class Cell3D{
  private:
    
   public:

    //! @name Public data types.
    //@{
    //! The type of the container for all geometric moments
    //typedef TaylorDerivativesContainer<ThreeD,double> GeometricMoments;
    typedef TaylorDerivativesContainer<double> GeometricMoments;
    //! The type of an individual geometric moment
    //typedef GeometricMoments::Derivative GeomMoment;
    typedef TaylorDerivativesContainer<double>::Derivative GeomMoment; 
    //@}

    //! @name Public variables.
    //@{
    Vector3D           Xc;   // Location of cell center.
    int             I,J,K;   // (i,j,k) indices for cell.
    double              V;   // Cell volume
    double       Jacobian;   // Determinant of Jacobian
    Vector3D          dXc;
	                     // Made public so can access them.
    //@}
		      
    /* Creation, copy, and assignment constructors. */
    Cell3D(void): _GeomCoeff_(0) {
       Xc.x = ONE; Xc.y = ONE; Xc.z=ONE; I = 0; J =0; K=0; V = ONE; Jacobian=ONE;
    }

    Cell3D(const Cell3D &Cell): _GeomCoeff_(Cell._GeomCoeff_)  {
       Xc = Cell.Xc; I = Cell.I; J = Cell.J; K=Cell.K; V = Cell.V; Jacobian = Cell.Jacobian;
    }

    Cell3D(const Vector3D &V): _GeomCoeff_(0) {
      Xc = V;
    }

    Cell3D(const double &xx, const double &yy, const double &zz): _GeomCoeff_(0) {
      Xc.x = xx; Xc.y = yy; Xc.z=zz;
    }
    
    /*!
     * Constructor with reconstruction order.
     * The reconstruction order determines the size of the geometric moments container.
     * If high-order reconstructions of different orders are used,
     * the OrderOfReconstruction should be the highest value of them.
     */
    Cell3D(const int &OrderOfReconstruction): _GeomCoeff_(OrderOfReconstruction){
      Xc.x = ONE; Xc.y = ONE; Xc.z=ONE; I = 0; J =0; K=0; V = ONE;
    };

    /* Destructor. */
    // ~Cell2D(void);
    // Use automatically generated destructor.

    /* Set cell center location. */
    void setloc(void);
    void setloc(const Cell3D &Cell);
    void setloc(const Vector3D &V);
    void setloc(const double &xx, const double &yy,const double &zz);
    
    /* Set cell (i,j,k) indices. */
    void setindex(void);
    void setindex(const Cell3D &Cell);
    void setindex(const int II, const int JJ, const int KK);

    /******* Required Grid Information for CENO Reconstruction  *************/

    //! Generate the container for the cell's geometric coefficients.
    void SetGeomCoeffContainerSize(const int OrderOfReconstruction){
      return _GeomCoeff_.GenerateContainer(OrderOfReconstruction);
    }

    //! @name Field access functions.
    //@{
    //! Get access to the array of geometric coefficients
    const GeometricMoments & GeomCoeff(void) const { return _GeomCoeff_;}
    //! Get access to the array of geometric coefficients
    GeometricMoments & GeomCoeff(void) {return _GeomCoeff_;}
    //! Get access to the value of the geometric coefficient with x-power 'p1' and y-power 'p2'
    const double & GeomCoeffValue(const int &p1, const int &p2, const int &p3) const {return _GeomCoeff_(p1,p2,p3);}
    //! Get access to the value of the geometric coefficient with x-power p1 and y-power p2
    double & GeomCoeffValue(const int &p1, const int &p2, const int &p3) {return _GeomCoeff_(p1,p2,p3);}
    //! Get access to the value of the geometric coefficient which is stored in the 'position' place
    const double & GeomCoeffValue(const int &position) const {return _GeomCoeff_(position).D();}
    //! Get access to the value of the geometric coefficient which is stored in the 'position' place
    double & GeomCoeffValue(const int &position){return _GeomCoeff_(position).D();}
    //! Get access to the geometric coefficient (i.e. powers and values) which is stored in the 'position' place
    GeomMoment & GeomCoeff(const int &position){return _GeomCoeff_(position);}
    //@}


    /* Assignment operator. */
    // Cell2D operator = (const Cell2D &Cell);
    // Use automatically generated assignment operator.

    /* Relational operators. */
    friend int operator ==(const Cell3D &Cell1,
			   const Cell3D &Cell2);
    friend int operator !=(const Cell3D &Cell1,
			   const Cell3D &Cell2);
    
    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file,
				 const Cell3D &Cell);
    friend istream &operator >> (istream &in_file,
				 Cell3D &Cell);
 

private:    
  //! Area integrals of cell geometric moments with respect to the cell centroid 
  GeometricMoments _GeomCoeff_;  //    
};

/******************************************************************
 * Cell3D::setloc -- Set cell center location.                    *
 ******************************************************************/
inline void Cell3D::setloc(void) {
  Xc.x = ONE; Xc.y = ONE; Xc.z=ONE;
}

inline void Cell3D::setloc(const Cell3D &Cell) {
  Xc = Cell.Xc;
}

inline void Cell3D::setloc(const Vector3D &V) {
  Xc = V;
}

inline void Cell3D::setloc(const double &xx, const double &yy, const double &zz) {
  Xc.x = xx; Xc.y = yy; Xc.z=zz;
}

/******************************************************************
 * Cell3D::setindex -- Set cell (i,j) indices.                    *
 ******************************************************************/
inline void Cell3D::setindex(void) {
  I = 0; J = 0; K = 0;
}

inline void Cell3D::setindex(const Cell3D &Cell) {
  I = Cell.I; J = Cell.J; K=Cell.K;
}

inline void Cell3D::setindex(const int II, const int JJ, const int KK) {
  I = II; J = JJ; K=KK;
}

/******************************************************************
 * Cell3D -- Relational operators.                                *
 ******************************************************************/
inline int operator ==(const Cell3D &Cell1,
		       const Cell3D &Cell2) {
    return (Cell1.Xc == Cell2.Xc && Cell1.V == Cell2.V);
}

inline int operator !=(const Cell3D &Cell1,
		       const Cell3D &Cell2) {
  return (Cell1.Xc != Cell2.Xc || Cell1.V != Cell2.V);
}

/******************************************************************
 * Cell3D -- Input-output operators.                              *
 ******************************************************************/
inline ostream &operator << (ostream &out_file,
			     const Cell3D &Cell) {
    out_file << " " << Cell.I << " " << Cell.J << " "<<Cell.K<<" "<<Cell.Xc;
    out_file.setf(ios::scientific);
    out_file << " " << Cell.V;
    out_file.unsetf(ios::scientific);
    out_file.unsetf(ios::scientific);
    out_file << " " << Cell.GeomCoeff();
    return (out_file);
}

inline istream &operator >> (istream &in_file,
			     Cell3D &Cell) {
    in_file.setf(ios::skipws);
    in_file >> Cell.I >> Cell.J >>Cell.K; 
    in_file.unsetf(ios::skipws);
    in_file >> Cell.Xc;
    in_file.setf(ios::skipws);
    in_file >> Cell.V;
    in_file.unsetf(ios::skipws);
    in_file >> Cell.GeomCoeff();
    return (in_file);
}

/* Define the basic 3D node class. */

/*********************************************************
 * Class: Node3D                                         *
 *                                                       *
 * Member functions                                      *
 *      X       -- Return 3D vector containing location  *
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
class Node3D{
  private:
  public:
    Vector3D         X;   // Node location.
	                  // Made public so can access them.
		      
    /* Creation, copy, and assignment constructors. */
    Node3D(void) {
      X.x = ONE; X.y = ONE; X.z=ONE;
    }

    Node3D(const Node3D &Node) {
       X = Node.X;
    }

    Node3D(const Vector3D &V) {
       X = V;
    }

    Node3D(const double &xx, const double &yy, const double &zz) {
      X.x = xx; X.y = yy; X.z=zz;
    }
    
    /* Destructor. */
    // ~Node2D(void);
    // Use automatically generated destructor.

    /* Set cell center location. */
    void setloc(void);
    void setloc(const Node3D &Node);
    void setloc(const Vector3D &V);
    void setloc(const double &xx, const double &yy, const double &zz);
    
    double tetvolume(const Vector3D V1,const  Node3D N2,const  Node3D N3,const  Vector3D V4);

    /* Assignment operator. */
    // Node2D operator = (const Node2D &Node);
    // Use automatically generated assignment operator.

    /* Relational operators. */
    friend int operator ==(const Node3D &Node1,
			   const Node3D &Node2);
    friend int operator !=(const Node3D &Node1,
			   const Node3D &Node2);
    
    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file,
				 const Node3D &Node);
    friend istream &operator >> (istream &in_file,
				 Node3D &Node);
    
};

/******************************************************************
 * Node3D::setloc -- Set cell center location.                    *
 ******************************************************************/
inline void Node3D::setloc(void) {
  X.x = ONE; X.y = ONE; X.z = ONE;
}

inline void Node3D::setloc(const Node3D &Node) {
  X = Node.X;
}

inline void Node3D::setloc(const Vector3D &V) {
  X = V;
}

inline void Node3D::setloc(const double &xx, const double &yy, const double &zz) {
  X.x = xx; X.y = yy; X.z=zz;
}

inline double tetvolume(const Vector3D V1,const  Node3D N2,const  Node3D N3,const  Vector3D V4){
  Vector3D A,B,C; 
  A = (N2.X-V1);
  B = (N3.X-V1);
  C = (V4-V1);
  return ((fabs(A*(B^C)))/6); 
}

/******************************************************************
 * Node3D -- Relational operators.                                *
 ******************************************************************/
inline int operator ==(const Node3D &Node1,
		       const Node3D &Node2) {
    return (Node1.X == Node2.X);
}

inline int operator !=(const Node3D &Node1,
		       const Node3D &Node2) {
    return (Node1.X != Node2.X);
}

/******************************************************************
 * Node3D -- Input-output operators.                              *
 ******************************************************************/
inline ostream &operator << (ostream &out_file,
			     const Node3D &Node) {
    out_file << Node.X;
    return (out_file);
}

inline istream &operator >> (istream &in_file,
			     Node3D &Node) {
    in_file >> Node.X; 
    return (in_file);
}

#endif /* _CELL3D_INCLUDED  */


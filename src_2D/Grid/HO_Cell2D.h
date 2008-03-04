/*!\file HO_Cell2D.h
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
#include "../HighOrderReconstruction/TaylorDerivatives.h"   /* Include Taylor derivatives container header files. */


/* Define the classes. */

/*!
 * \class Cell2D_HO
 *
 * \verbatim
 * Member functions
 *     Xc       -- Return 2D vector containing location of the cell
 *                 center (centroid).
 *     I        -- Return I index for cell.
 *     J        -- Return J index for cell.
 *     A        -- Return the area of the cell.
 *     setloc   -- Set cell center location (position).
 *     setindex -- Set cell (i,j) indices.
 *     GeomCoeff-- Area integrals of geometric moments
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
class Cell2D_HO{
public:

  //! @name Public data types.
  //@{
  //! The type of the container for all geometric moments
  typedef TaylorDerivativesContainer<TwoD,double> GeometricMoments;
  //! The type of an individual geometric moment
  typedef GeometricMoments::Derivative GeomMoment;
  //@}

  //! @name Public variables.
  //@{
  Vector2D     Xc;    //!< Location of cell center.
  int        I, J;    //!< (i,j) indices for cell.
  double        A;    //!< Cell area.
  //@}


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

  //! Constructor with reconstruction order specified
  Cell2D_HO(const int &OrderOfReconstruction);
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

  //! @name Field access functions.
  //@{
  //! Get access to the array of geometric coefficients
  const GeometricMoments & GeomCoeff(void) const { return _GeomCoeff_;}
  //! Get access to the array of geometric coefficients
  GeometricMoments & GeomCoeff(void) {return _GeomCoeff_;}
  //! Get access to the value of the geometric coefficient with x-power 'p1' and y-power 'p2'
  const double & GeomCoeffValue(const int &p1, const int &p2) const {return _GeomCoeff_(p1,p2);}
  //! Get access to the value of the geometric coefficient with x-power p1 and y-power p2
  double & GeomCoeffValue(const int &p1, const int &p2) {return _GeomCoeff_(p1,p2);}
  //! Get access to the value of the geometric coefficient which is store in the 'position' place
  const double & GeomCoeffValue(const int &position) const {return _GeomCoeff_(position,true,true,true).D();}
  //! Get access to the value of the geometric coefficient which is store in the 'position' place
  double & GeomCoeffValue(const int &position){return _GeomCoeff_(position,true,true,true).D();}
  //! Get access to the geometric coefficient (i.e. powers and values) which is store in the 'position' place
  GeomMoment & GeomCoeff(const int &position){return _GeomCoeff_(position,true,true,true);}
  //@}

  //! Output invariant properties to cell translation and rotation.
  void Output_Translation_Rotation_Invariant_Properties(ostream &out_file);

  //! Set the container of the cell geometric coefficients for usage.
  void SetGeomCoeffContainerSize(const int OrderOfReconstruction){
    return _GeomCoeff_.GenerateContainer(OrderOfReconstruction);
  }

  //! @name Relational operators.
  //@{
  friend bool operator ==(const Cell2D_HO &Cell1, const Cell2D_HO &Cell2);
  friend bool operator !=(const Cell2D_HO &Cell1, const Cell2D_HO &Cell2);
  //@}

  //! @name Input-output operators.
  //@{
  friend ostream & operator << (ostream &out_file, const Cell2D_HO &Cell);
  friend istream & operator >> (istream &in_file, Cell2D_HO &Cell);
  //@}

private:    
  //! Area integrals of cell geometric moments with respect to the cell centroid 
  GeometricMoments _GeomCoeff_;  // 
};


/*!
 * Default constructor
 */
inline Cell2D_HO::Cell2D_HO(void): Xc(ONE), I(0), J(0), A(ONE), _GeomCoeff_(0)
{
  // 
}

/*!
 * Copy constructor
 */
inline Cell2D_HO::Cell2D_HO(const Cell2D_HO &Cell):
  Xc(Cell.Xc), I(Cell.I), J(Cell.J), A(Cell.A), _GeomCoeff_(Cell._GeomCoeff_)
{
  // 
}

/*!
 * Constructor with centroid location
 */
inline Cell2D_HO::Cell2D_HO(const Vector2D &V): Xc(V), I(0), J(0), A(ONE), _GeomCoeff_(0)
{
  // 
}

/*!
 * Constructor with centroid locations
 */
inline Cell2D_HO::Cell2D_HO(const double &xx, const double &yy):
  I(0), J(0), A(ONE), _GeomCoeff_(0)
{
  Xc.x = xx; 
  Xc.y = yy;
}

/*!
 * Constructor with reconstruction order.
 * The reconstruction order determines the size of the geometric moments container.
 * If high-order reconstructions of different orders are used,
 * the OrderOfReconstruction should be the highest value of them.
 */
inline Cell2D_HO::Cell2D_HO(const int &OrderOfReconstruction):
  Xc(ONE), I(0), J(0), A(ONE),
  _GeomCoeff_(OrderOfReconstruction)
{
  // 
};

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

/*!
 * Set cell (i,j) indices.
 */
inline void Cell2D_HO::setindex(void) {
  I = 0; J = 0;
}

/*!
 * Set cell (i,j) indices.
 */
inline void Cell2D_HO::setindex(const Cell2D_HO &Cell) {
    I = Cell.I; J = Cell.J;
}

/*!
 * Set cell (i,j) indices.
 */
inline void Cell2D_HO::setindex(const int II, const int JJ) {
    I = II; J = JJ;
}

/*!
 * Output the cell geometric properties which are 
 * invariant to translation and rotation operations.
 * These properties are: 
 * cell indexes, area, geometric moments with respect to cell centroid.
 */
inline void Cell2D_HO::Output_Translation_Rotation_Invariant_Properties(ostream &out_file){

  out_file << " " << I << " " << J;
  out_file.setf(ios::scientific);
  out_file << " " << A << "\n";
  out_file.unsetf(ios::scientific);
  out_file << " " << GeomCoeff();
}

/*!
 * Equal operator
 */
inline bool operator ==(const Cell2D_HO &Cell1,
			const Cell2D_HO &Cell2) {
  return (Cell1.Xc == Cell2.Xc && Cell1.A == Cell2.A && Cell1._GeomCoeff_ == Cell2._GeomCoeff_);
}

/*!
 * Not-equal operator
 */
inline bool operator !=(const Cell2D_HO &Cell1,
			const Cell2D_HO &Cell2) {
  return (Cell1.Xc != Cell2.Xc || Cell1.A != Cell2.A || Cell1._GeomCoeff_ != Cell2._GeomCoeff_);
}

/*!
 * Output operator
 */
inline ostream &operator << (ostream &out_file,
			     const Cell2D_HO &Cell) {
  out_file << " " << Cell.I << " " << Cell.J << Cell.Xc;
  out_file.setf(ios::scientific);
  out_file << " " << Cell.A << "\n";
  out_file.unsetf(ios::scientific);
  out_file << " " << Cell.GeomCoeff();
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
  in_file >> Cell.GeomCoeff();
  return (in_file);
}

/*!
 * \class Cell2D_Cartesian_HO                            
 * \brief Define the cell class for a uniform 2D Cartesian mesh.
 * \verbatim
 * Member functions                                      
 *      xc      -- Return 2D vector containing location  
 *                 of cell center.                       
 *      dx      -- Return 2D vector of cell dimensions.  
 *                                                       
 *                                                       
 *                 xnNW-----xfN----xnNE                  
 *                   |              |                    
 *                   |              |                    
 *       Cell       xfW      xc    xfE                   
 *                   |              |                    
 *                   |              |                    
 *                 xnSW-----xfS----xnSE                  
 *                                                       
 *                      UnitTests                                 
 *      A       -- Return the area of the cell.          
 *      xnNW    -- Return 2D vector containing location  
 *                 of north-west node.                   
 *      xnNE    -- Return 2D vector containing location  
 *                 of north-east node.                    
 *      xnSE    -- Return 2D vector containing location  
 *                 of south-east node.                   
 *      xnSW    -- Return 2D vector containing location  
 *                 of south-west node.                   
 *      xfN     -- Return 2D vector containing location  
 *                 of center of north cell face.         
 *      xfS     -- Return 2D vector containing location  
 *                 of center of south cell face.         
 *      xfE     -- Return 2D vector containing location  
 *                 of center of east cell face.          
 *      xfW     -- Return 2D vector containing location  
 *                 of center of west cell face.          
 *      setsize -- Set cell dimensions.                  
 *      setloc  -- Set cell center location (position).   
 *                                                       
 * Member operators                                      
 *      C -- a cell                                      
 *                                                       
 * C = C;                                                
 * C == C;                                               
 * C != C;              UnitTests                                 
 * cout << C; (output function)                          
 * cin  >> C; (input function)                           
 * \endverbatim                                                       
 */
class Cell2D_Cartesian_HO{
private:
public:
  Vector2D         xc;   //!< Location of cell center.
  static Vector2D  dx;   //!< Vector of cell lengths.
  
  static Cell2D_Cartesian_HO Cell2D_Cartesian_HO_ONE;
  
  //! @name Constructors.
  //@{
  
  //!Default constructor
  Cell2D_Cartesian_HO(void);
  
  //!Copy Constructor
  Cell2D_Cartesian_HO(const Cell2D_Cartesian_HO &Cell);
  
  //! Constructor with given 2D vector
  Cell2D_Cartesian_HO(const Vector2D &V);
  
  //! Constructor with given coordinates
  Cell2D_Cartesian_HO(const double &xx, const double &yy);
  //@}
  
  //! @name Cell area.
  //@{
  double A(void);
  double A(void) const;
  //@}

  //!@name Node locations. 
  //@{
  Vector2D xnNW(void);
  Vector2D xnNW(void) const;
  
  Vector2D xnNE(void);
  Vector2D xnNE(void) const;
  
  Vector2D xnSE(void);
  Vector2D xnSE(void) const;
  
  Vector2D xnSW(void);
  Vector2D xnSW(void) const;
  //@}

  //!@name Face midpoints.
  //@{
  Vector2D xfN(void);
  Vector2D xfN(void) const;
  
  Vector2D xfS(void);
  Vector2D xfS(void) const;
  
  Vector2D xfE(void);
  Vector2D xfE(void) const;
  
  Vector2D xfW(void);
  Vector2D xfW(void) const;
  //@}

  //!@name Set cell dimensions. 
  //@{
  void setsize(void);
  void setsize(const Vector2D &xx);
  void setsize(const double &xx, const double &yy);
  //@}
  
  //!@name Set cell center location.
  //@{
  void setloc(void);
  void setloc(const Cell2D_Cartesian_HO &Cell);
  void setloc(const Vector2D &V);
  void setloc(const double &xx, const double &yy);
  //@}
  
  //!@name Relational operators.
  //@{
  friend int operator ==(const Cell2D_Cartesian_HO &Cell1,
			 const Cell2D_Cartesian_HO &Cell2);
  friend int operator !=(const Cell2D_Cartesian_HO &Cell1,
			 const Cell2D_Cartesian_HO &Cell2);
  //@}

  //!@name Input-output operators.
  //@{
  friend ostream &operator << (ostream &out_file,
			       const Cell2D_Cartesian_HO &Cell);
  friend istream &operator >> (istream &in_file,
			       Cell2D_Cartesian_HO &Cell);
  //@}
};

//!Default constructor
inline Cell2D_Cartesian_HO::Cell2D_Cartesian_HO(void) : xc(ONE)
{
  //
}

//!Copy Constructor
inline Cell2D_Cartesian_HO::Cell2D_Cartesian_HO(const Cell2D_Cartesian_HO &Cell) : xc (Cell.xc)
{
  //
}

//! Constructor with centroid location
inline Cell2D_Cartesian_HO::Cell2D_Cartesian_HO(const Vector2D &V) : xc(V) 
{
  //
}

//! Constructor with centroid location
inline Cell2D_Cartesian_HO::Cell2D_Cartesian_HO(const double &xx, const double &yy) :  xc(ONE)
{
  //
}

/*!
 * Return cell area. 
 */
inline double Cell2D_Cartesian_HO::A(void) {
  return (dx.x*dx.y);
}

/*!
 * Return cell area. 
 */
inline double Cell2D_Cartesian_HO::A(void) const {
  return (dx.x*dx.y);
}

/*!
 *Get North-West cell node.
 */
inline Vector2D Cell2D_Cartesian_HO::xnNW(void) {
  return (Vector2D(xc-HALF*dx.x*Vector2D_NX+HALF*dx.y*Vector2D_NY));
}

/*!
 *Get North-West cell node.
 */
inline Vector2D Cell2D_Cartesian_HO::xnNW(void) const {
  return (Vector2D(xc-HALF*dx.x*Vector2D_NX+HALF*dx.y*Vector2D_NY));
}

/*!
 *Get North-East cell node.
 */
inline Vector2D Cell2D_Cartesian_HO::xnNE(void) {
  return (Vector2D(xc+HALF*dx.x*Vector2D_NX+HALF*dx.y*Vector2D_NY));
}

/*!
 *Get North-East cell node.
 */
inline Vector2D Cell2D_Cartesian_HO::xnNE(void) const {
  return (Vector2D(xc+HALF*dx.x*Vector2D_NX+HALF*dx.y*Vector2D_NY));
}

/*!
 *Get South-East cell node.
 */
inline Vector2D Cell2D_Cartesian_HO::xnSE(void) {
  return (Vector2D(xc+HALF*dx.x*Vector2D_NX-HALF*dx.y*Vector2D_NY));
}


/*!
 *Get South-East cell node.
 */
inline Vector2D Cell2D_Cartesian_HO::xnSE(void) const {
  return (Vector2D(xc+HALF*dx.x*Vector2D_NX-HALF*dx.y*Vector2D_NY));
}


/*!
 *Get South-West cell node.
 */
inline Vector2D Cell2D_Cartesian_HO::xnSW(void) {
  return (Vector2D(xc-HALF*dx.x*Vector2D_NX-HALF*dx.y*Vector2D_NY));
}


/*!
 *Get South-West cell node.
 */
inline Vector2D Cell2D_Cartesian_HO::xnSW(void) const {
  return (Vector2D(xc-HALF*dx.x*Vector2D_NX-HALF*dx.y*Vector2D_NY));
}

/*!
 * Get North cell face.
 */
inline Vector2D Cell2D_Cartesian_HO::xfN(void) {
  return (Vector2D(xc+HALF*dx.y*Vector2D_NY));
}

/*!
 * Get North cell face.
 */
inline Vector2D Cell2D_Cartesian_HO::xfN(void) const {
  return (Vector2D(xc+HALF*dx.y*Vector2D_NY));
}

/*!
 * Get South cell face.
 */
inline Vector2D Cell2D_Cartesian_HO::xfS(void) {
  return (Vector2D(xc-HALF*dx.y*Vector2D_NY));
}

/*!
 * Get South cell face.
 */
inline Vector2D Cell2D_Cartesian_HO::xfS(void) const {
  return (Vector2D(xc-HALF*dx.y*Vector2D_NY));
}

/*!
 * Get East cell face.
 */
inline Vector2D Cell2D_Cartesian_HO::xfE(void) {
  return (Vector2D(xc+HALF*dx.x*Vector2D_NX));
}

/*!
 * Get East cell face.
 */
inline Vector2D Cell2D_Cartesian_HO::xfE(void) const {
  return (Vector2D(xc+HALF*dx.x*Vector2D_NX));
}

/*!
 * Get West cell face.
 */
inline Vector2D Cell2D_Cartesian_HO::xfW(void) {
  return (Vector2D(xc-HALF*dx.x*Vector2D_NX));
}

/*!
 * Get West cell face.
 */
inline Vector2D Cell2D_Cartesian_HO::xfW(void) const {
  return (Vector2D(xc-HALF*dx.x*Vector2D_NX));
}

/*!
 * Set default cell dimensions. dx.x = one, dx.y = one.
 */
inline void Cell2D_Cartesian_HO::setsize(void) {
  dx.x = ONE; dx.y = ONE;
}

/*!
 * Set cell dimensions.
 */
inline void Cell2D_Cartesian_HO::setsize(const Vector2D &xx) {
  dx = xx;
}

/*!
 * Set cell dimensions.
 */
inline void Cell2D_Cartesian_HO::setsize(const double &xx, const double &yy) {
  dx.x = xx; dx.y = yy;
}

/*!
 *Set default cell center location. xc.x = one, xc.y = one.
 */
inline void Cell2D_Cartesian_HO::setloc(void) {
  xc.x = ONE; xc.y = ONE;
}

/*!
 *Set cell center location.
 */
inline void Cell2D_Cartesian_HO::setloc(const Cell2D_Cartesian_HO &Cell) {
  xc = Cell.xc;
}

/*!
 *Set cell center location.
 */
inline void Cell2D_Cartesian_HO::setloc(const Vector2D  &V) {
  xc = V;
}

/*!
 *Set cell center location.
 */
inline void Cell2D_Cartesian_HO::setloc(const double &xx, const double &yy) {
  xc.x = xx; xc.y = yy;
}

/*!
 *Equal operator.
 */
inline int operator ==(const Cell2D_Cartesian_HO &Cell1,
		       const Cell2D_Cartesian_HO &Cell2) {
  return (Cell1.xc == Cell2.xc);
}

/*!
 * Not-equal operator.
 */
inline int operator !=(const Cell2D_Cartesian_HO &Cell1,
		       const Cell2D_Cartesian_HO &Cell2) {
  return (Cell1.xc != Cell2.xc);
}

/*!
 * Output operator.
 */
inline ostream &operator << (ostream &out_file,
			     const Cell2D_Cartesian_HO &Cell) {
  out_file << Cell.xc;
  return (out_file);
}

/*!
 * Input operator.
 */
inline istream &operator >> (istream &in_file,
			     Cell2D_Cartesian_HO &Cell) {
  in_file >> Cell.xc; 
  return (in_file);
}

/*!
 * \class Cell2D_Quad_HO                                  
 * \brief Define the cell class for a 2D quadrilateral mesh.
 * \verbatim                                                     
 * Member functions                                     
 *      xc      -- Return 2D vector containing location 
 *                 of cell center.                      
 *      xnNW    -- Return 2D vector containing location 
 *                 of north-west node.                  
 *      xnNE    -- Return 2D vector containing location 
 *                 of north-east node.                  
 *      xnSE    -- Return 2D vector containing location 
 *                 of south-east node.                  
 *      xnSW    -- Return 2D vector containing location 
 *                 of south-west node.                  
 *                                                      
 *                                                      
 *                 xnNW-----xfN----xnNE                 
 *                   |              |                   
 *                   |              |                   
 *       Cell       xfW      xc    xfE                  
 *                   |              |                   
 *                   |              |                   
 *                 xnSW-----xfS----xnSE                 
 *                                                      
 *                                                      
 *      A       -- Return the area of the cell.         
 *      xfN     -- Return 2D vector containing location 
 *                 of midpoint of north cell face.      
 *      xfS     -- Return 2D vector containing location 
 *                 of midpoint of south cell face.      
 *      xfE     -- Return 2D vector containing location 
 *                 of midpoint of east cell face.       
 *      xfW     -- Return 2D vector containing location 
 *                 of midpoint of west cell face.       
 *      lfN     -- Return the length of north cell face.
 *      lfS     -- Return the length of south cell face.
 *      lfE     -- Return the length of east cell face. 
 *      lfW     -- Return the length of west cell face. 
 *      setnodes -- Set cell node locations.             
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

class Cell2D_Quad_HO{
private:
public:
  Vector2D         xc;   //!< Location of cell center.
  Vector2D       xnNW;   //!< Location of north-west node.
  Vector2D       xnNE;   //!< Location of north-east node.
  Vector2D       xnSE;   //!< Location of south-east node.
  Vector2D       xnSW;   //!< Location of south-west node.
 
  //! @name Constructors
  //@{
  //! Default constructor.
  Cell2D_Quad_HO(void) ;
  
  //! Copy constructor.
  Cell2D_Quad_HO(const Cell2D_Quad_HO &Cell) ;
  
  //! Constructor with centroid location.
  Cell2D_Quad_HO(const Vector2D &V1, const Vector2D &V2,
		 const Vector2D &V3, const Vector2D &V4);
  
  //! Constructor with centroid locations.
  Cell2D_Quad_HO(const double &x1, const double &y1,
		 const double &x2, const double &y2,
		 const double &x3, const double &y3,
		 const double &x4, const double &y4);
  //@}
  
  //! @name Cell Area.
  //@{
  double A(void);
  double A(void) const;
  //@}
  
  //! @name Face midpoints.
  //@{
  Vector2D xfN(void);
  Vector2D xfN(void) const;
  
  Vector2D xfS(void);
  Vector2D xfS(void) const;
  
  Vector2D xfE(void);
  Vector2D xfE(void) const;
  
  Vector2D xfW(void);
  Vector2D xfW(void) const;
  //@}
  
  //! @name  Face lengths.
  //@{
  double lfN(void);
  double lfN(void) const;
  
  double lfS(void);
  double lfS(void) const;
  
  double lfE(void);
  double lfE(void) const;
  
  double lfW(void);
  double lfW(void) const;
  //@}
  
  //! @name Set cell node locations.
  //@{
  void setnodes(void);
  void setnodes(const Cell2D_Quad_HO &Cell);
  void setnodes(const Vector2D &V1, const Vector2D &V2,
		const Vector2D &V3, const Vector2D &V4);
  void setnodes(const double &x1, const double &y1,
		const double &x2, const double &y2,
		const double &x3, const double &y3,
		const double &x4, const double &y4);
  //@}

  //! @name Relational operators.
  //@{
  friend int operator ==(const Cell2D_Quad_HO &Cell1,
			 const Cell2D_Quad_HO &Cell2);
  friend int operator !=(const Cell2D_Quad_HO &Cell1,
			 const Cell2D_Quad_HO &Cell2);
  //@}
  
  //! @name Input-output operators.
  //@{
  friend ostream &operator << (ostream &out_file,
			       const Cell2D_Quad_HO &Cell);
  friend istream &operator >> (istream &in_file,
			       Cell2D_Quad_HO &Cell);
  //@}
  
};

//! Default constructor.
inline Cell2D_Quad_HO::Cell2D_Quad_HO(void) {
  xnNW.x = ZERO; xnNW.y = ONE; xnNE.x = ONE; xnNE.y = ONE;
  xnSE.x = ONE; xnSE.y = ZERO; xnSW.x = ZERO; xnSW.y = ZERO;
  xc = (xnNW+xnNE+xnSE+xnSW)/FOUR; 
}

//! Copy constructor.
inline Cell2D_Quad_HO::Cell2D_Quad_HO(const Cell2D_Quad_HO &Cell) {
    xnNW = Cell.xnNW; xnNE = Cell.xnNE; xnSE = Cell.xnSE; 
    xnSW = Cell.xnSW; xc = Cell.xc;
  }
  
//! Constructor with centroid location.
inline Cell2D_Quad_HO::Cell2D_Quad_HO(const Vector2D &V1, const Vector2D &V2,
	      const Vector2D &V3, const Vector2D &V4) {
    xnNW = V4; xnNE = V3; xnSE = V2; xnSW = V1; 
    xc = (xnNW+xnNE+xnSE+xnSW)/FOUR;
  }
  
//! Constructor with centroid locations.
inline Cell2D_Quad_HO::Cell2D_Quad_HO(const double &x1, const double &y1,
	      const double &x2, const double &y2,
	      const double &x3, const double &y3,
	      const double &x4, const double &y4) {
    xnNW.x = x4; xnNW.y = y4; xnNE.x = x3; xnNE.y = y3;
    xnSE.x = x2; xnSE.y = y2; xnSW.x = x1; xnSW.y = y1;
    xc = (xnNW+xnNE+xnSE+xnSW)/FOUR;
  }

/*!
 *Return cell area.
 */
inline double Cell2D_Quad_HO::A(void) {
  return (HALF*(((xnSE-xnSW)^(xnNW-xnSW))+((xnNE-xnNW)^(xnNE-xnSE))));
}

/*!
 *Return cell area.
 */
inline double Cell2D_Quad_HO::A(void) const {
  return (HALF*(((xnSE-xnSW)^(xnNW-xnSW))+((xnNE-xnNW)^(xnNE-xnSE))));
}

/*!
 *Get North face midpoints.
 */
inline Vector2D Cell2D_Quad_HO::xfN(void) {
  return (HALF*(xnNE+xnNW));
}

/*!
 *Get North face midpoints.
 */
inline Vector2D Cell2D_Quad_HO::xfN(void) const {
  return (HALF*(xnNE+xnNW));
}

/*!
 *Get South face midpoints.
 */
inline Vector2D Cell2D_Quad_HO::xfS(void) {
  return (HALF*(xnSE+xnSW));
}

/*!
 *Get South face midpoints.
 */
inline Vector2D Cell2D_Quad_HO::xfS(void) const {
  return (HALF*(xnSE+xnSW));
}

/*!
 *Get East face midpoints.
 */
inline Vector2D Cell2D_Quad_HO::xfE(void) {
  return (HALF*(xnSE+xnNE));
}

/*!
 *Get East face midpoints.
 */
inline Vector2D Cell2D_Quad_HO::xfE(void) const {
  return (HALF*(xnSE+xnNE));
}

/*!
 *Get West face midpoints.
 */
inline Vector2D Cell2D_Quad_HO::xfW(void) {
  return (HALF*(xnSW+xnNW));
}

/*!
 *Get West face midpoints.
 */
inline Vector2D Cell2D_Quad_HO::xfW(void) const {
  return (HALF*(xnSW+xnNW));
}

/*!
 *Get North face lengths.
 */
inline double Cell2D_Quad_HO::lfN(void) {
  return (abs(xnNE-xnNW));
}

/*!
 *Get North face lengths.
 */
inline double Cell2D_Quad_HO::lfN(void) const {
  return (abs(xnNE-xnNW));
}

/*!
 *Get South face lengths.
 */
inline double Cell2D_Quad_HO::lfS(void) {
  return (abs(xnSW-xnSE));
}

/*!
 *Get South face lengths.
 */
inline double Cell2D_Quad_HO::lfS(void) const {
  return (abs(xnSW-xnSE));
}

/*!
 *Get East face lengths.
 */
inline double Cell2D_Quad_HO::lfE(void) {
  return (abs(xnNE-xnSE));
}

/*!
 *Get East face lengths.
 */
inline double Cell2D_Quad_HO::lfE(void) const {
  return (abs(xnNE-xnSE));
}

/*!
 *Get West face lengths.
 */
inline double Cell2D_Quad_HO::lfW(void) {
  return (abs(xnNW-xnSW));
}

/*!
 *Get West face lengths.
 */
inline double Cell2D_Quad_HO::lfW(void) const {
  return (abs(xnNW-xnSW));
}

/*!
 *Set default cell node locations.
 */
inline void Cell2D_Quad_HO::setnodes(void) {
  xnNW.x = ZERO; xnNW.y = ONE; xnNE.x = ONE; xnNE.y = ONE;
  xnSE.x = ONE; xnSE.y = ZERO; xnSW.x = ZERO; xnSW.y = ZERO;
  xc = (xnNW+xnNE+xnSE+xnSW)/FOUR;
}

/*!
 *Set cell node locations .
 */
inline void Cell2D_Quad_HO::setnodes(const Cell2D_Quad_HO &Cell) {
  xnNW = Cell.xnNW; xnNE = Cell.xnNE; xnSE = Cell.xnSE; 
  xnSW = Cell.xnSW; xc = Cell.xc;
}

/*!
 *Set cell node locations.
 */
inline void Cell2D_Quad_HO::setnodes(const Vector2D &V1, const Vector2D &V2,
                                  const Vector2D &V3, const Vector2D &V4) {
  xnNW = V4; xnNE = V3; xnSE = V2; xnSW = V1; 
  xc = (xnNW+xnNE+xnSE+xnSW)/FOUR;
}

/*!
 *Set cell node locations.
 */
inline void Cell2D_Quad_HO::setnodes(const double &x1, const double &y1,
                                  const double &x2, const double &y2,
                                  const double &x3, const double &y3,
                                  const double &x4, const double &y4) {
  xnNW.x = x4; xnNW.y = y4; xnNE.x = x3; xnNE.y = y3;
  xnSE.x = x2; xnSE.y = y2; xnSW.x = x1; xnSW.y = y1;
  xc = (xnNW+xnNE+xnSE+xnSW)/FOUR;
}

/*!
 *Equal operator.
 */
inline int operator ==(const Cell2D_Quad_HO &Cell1,
		       const Cell2D_Quad_HO &Cell2) {
  return (Cell1.xnNE == Cell2.xnNE && Cell1.xnNW == Cell2.xnNW &&
	  Cell1.xnSW == Cell2.xnSW && Cell1.xnSE == Cell2.xnSE);
}

/*!
 *Non-Equal operator.
 */
inline int operator !=(const Cell2D_Quad_HO &Cell1,
		       const Cell2D_Quad_HO &Cell2) {
  return (Cell1.xnNE != Cell2.xnNE || Cell1.xnNW != Cell2.xnNW ||
	  Cell1.xnSW != Cell2.xnSW || Cell1.xnSE != Cell2.xnSE);
}

/*!
 * Output operator.
 */
inline ostream &operator << (ostream &out_file,
			     const Cell2D_Quad_HO &Cell) {
  out_file << Cell.xc;
  return (out_file);
}

/*!
 * Input operator.
 */
inline istream &operator >> (istream &in_file,
			     Cell2D_Quad_HO &Cell) {
  in_file >> Cell.xc; 
  return (in_file);
}

#endif /* _CELL2D_INCLUDED  */

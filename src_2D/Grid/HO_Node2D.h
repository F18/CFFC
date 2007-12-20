/*!\file HO_Node2D.h
  \brief Header file defining 2D node type. */

#ifndef _HO_NODE2D_INCLUDED
#define _HO_NODE2D_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../Math/Vector2D.h"

/* Define the basic 2D node class. */

/*!
 * \class Node2D_HO
 *
 * \begin{verbatim}
 * Member functions
 *      X       -- Return 2D vector containing location of node.
 *      setloc  -- Set node location (position).
 *
 * Member operators
 *      N -- a node
 *
 * N = N;
 * N == N;
 * N != N;
 * cout << N; (output function)
 * cin  >> N; (input function)
 * \end{verbatim}
 */
class Node2D_HO{
public:
  Vector2D         X;   //!< Node location.
		      
  //! @name Creation and copy constructors.
  //@{
  //! Default constructor.
  Node2D_HO(void): X(ONE) { }

  //! Value constructor.
  Node2D_HO(const double &Val): X(Val){ }

  //! Copy constructor.
  Node2D_HO(const Node2D_HO &Node): X(Node.X) { }

  //! Constructor with Vector2D
  Node2D_HO(const Vector2D &V): X(V) { }

  //! Constructor with values
  Node2D_HO(const double &xx, const double &yy){ X.x = xx ; X.y = yy; }
  //@}  

  //! @name Set cell center location.
  //@{
  void setloc(void);
  void setloc(const Node2D_HO &Node);
  void setloc(const Vector2D &V);
  void setloc(const double &xx, const double &yy);
  //@}

  //! @name Field access
  //@{
  double & x(void) { return X.x;}
  const double & x(void) const { return X.x;}
  double & y(void) { return X.y;}
  const double & y(void) const { return X.y;}
  Vector2D & Position(void) { return X; }
  const Vector2D & Position(void) const { return X; }
  //@}

  //! @name Relational operators.
  //@{
  friend int operator ==(const Node2D_HO &Node1,
			 const Node2D_HO &Node2);
  friend int operator !=(const Node2D_HO &Node1,
			 const Node2D_HO &Node2);
  //@}

  //! @name Input-output operators.
  //@{
  friend ostream &operator << (ostream &out_file,
			       const Node2D_HO &Node);
  friend istream &operator >> (istream &in_file,
			       Node2D_HO &Node);
  //@}

private:
    
};

/*!
 * Set cell center location.                    
 */
inline void Node2D_HO::setloc(void) {
    X.x = ONE; X.y = ONE;
}

/*!
 * Set cell center location.                    
 */
inline void Node2D_HO::setloc(const Node2D_HO &Node) {
    X = Node.X;
}

/*!
 * Set cell center location.                    
 */
inline void Node2D_HO::setloc(const Vector2D &V) {
    X = V;
}

/*!
 * Set cell center location.                    
 */
inline void Node2D_HO::setloc(const double &xx, const double &yy) {
    X.x = xx; X.y = yy;
}

/*!
 * Equal operator
 */
inline int operator ==(const Node2D_HO &Node1,
		       const Node2D_HO &Node2) {
  return (Node1.X == Node2.X);
}

/*!
 * Not-equal operator
 */
inline int operator !=(const Node2D_HO &Node1,
		       const Node2D_HO &Node2) {
  return (Node1.X != Node2.X);
}

/*!
 * output operator
 */
inline ostream &operator << (ostream &out_file,
			     const Node2D_HO &Node) {
  out_file << Node.X;
  return (out_file);
}

/*!
 * input operator
 */
inline istream &operator >> (istream &in_file,
			     Node2D_HO &Node) {
  in_file >> Node.X; 
  return (in_file);
}


#endif	// _HO_NODE2D_INCLUDED

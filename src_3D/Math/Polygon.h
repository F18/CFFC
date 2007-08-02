/**********************************************************************
 * Polygon.h: Header file defining polygon subroutines.               *
 **********************************************************************/

#ifndef _POLYGON_INCLUDED
#define _POLYGON_INCLUDED

// Include required C++ libraries.

#include <cassert>
#include <cstdlib>

using namespace std;

// Include math macro file.

#ifndef _MATH_MACROS_INCLUDED
#include "Math.h"
#endif // _MATH_MACROS_INCLUDED

// Include Vector2D header.

#ifndef _VECTOR2D_INCLUDED
#include "Vector2D.h"
#endif //_VECTOR2D_INCLUDED

// Include LinkedList header.

#ifndef _LINKEDLIST_INCLUDED
#include "LinkedList.h"
#endif //_LINKEDLIST_INCLUDED

// Include Spline2D header.

// #ifndef _SPLINE2D_INCLUDED
// #include "Spline2D.h"
// #endif //_SPLINE2D_INCLUDED

/*!
 * Class: Polygon
 *
 * @brief A class defined to handle polygon manipulations.
 *
 * A class defined to handle polygon manipulations such as determination
 * of the area or centroid of an n-sided polygon.  The nodes of the
 * polygon must be defined in counter-clockwise direction.  The
 * intersection or the union of two polygons can be found using a method
 * based on the Weiler-Atherton algorithm popular in computer graphics.
 * The polygon clipping algorithm can also be used to determine the
 * polygon of intersection between two unique polygons (Sutherland and
 * Hodgman [Comm of the ACM, 17(1):32-42, 1974]).  Note that all of the
 * algorithms used here require the polygons be simple polygons
 * (non-self intersecting) except for the 'point in polygon' algorithm
 * which determines whether or not a given point is internal to a polygon
 * (also known as ray-tracing).
 *
 * \verbatim
 * Member functions
 *      np               -- Number of nodes in n-sided polygon.
 *      X                -- Array of polygon nodes.
 *      A                -- Return the area of the polygon.
 *      Xc               -- Return the centroid of the polygon.
 *      allocate         -- Allocate memory for polygon.
 *      deallocate       -- Deallocate memory for polygon.
 *      copy             -- Copy the polygon.
 *      area             -- Returns the area of the polygon.
 *      centroid         -- Returns the centroid of the polygon.
 *      min_length       -- Returns the length of the shortest edge.
 *      point_in_polygon -- Determines if the given point is inside
 *                          the polygon.
 *      sort             -- Sort the polygon into counter-clockwise
 *                          order.
 *      reverse_order    -- Reverse the order of the polygon nodes.
 *      convert          -- Convert a linked list of vectors to a 
 *                          polygon.
 *      convert          -- Convert a list of four not-necessarily-unique
 *                          vectors into a polygon.
 *      convert          -- Convert a spline to a polygon.
 *
 * Member operators
 *      P -- a polygon.
 *
 * cout << P; -- Ouput function.
 * cin >> P; -- Input function.
 * \endverbatim
 *
 * \verbatim
 * Known bugs:
 * 1. The polygon clipping algorithm is not as robust for finding the
 *    polygon of intersection between two distinct polygons as the
 *    Weiler-Atherton in the sense that it is difficult to generate a
 *    counter-clockwise order of nodes for an arbitrary polygon of
 *    intersection.
 * 2. Unable to deal with holes.
 * \endverbatim
 */
class Polygon{
private:
public:
  int      np; //!< Number of nodes in n-sided polygon.
  Vector2D *X; //!< Array of polygon nodes.
  double    A; //!< Return the area of the polygon.
  Vector2D Xc; //!< Return the centroid of the polygon.

  //@{ @name Creation and copy constructors.

  //! Creation constructor.
  Polygon(void) {
    np = 0; X = NULL;
    A = ZERO; Xc = Vector2D_ZERO;
  }

  //! Creation constructor.
  Polygon(const int n) {
    np = n; X = new Vector2D[n];
    A = ZERO; Xc = Vector2D_ZERO;
  }

  //! Copy constructor.
  Polygon(const Polygon &P) {
    copy(P);
    //np = P.np; X = P.X;
    //A = P.A; Xc = P.Xc;
  }

  //! Destructor.
  ~Polygon(void) {
    delete []X; X = NULL; np = 0;
    A = ZERO; Xc = Vector2D_ZERO;
  }

  //@}

  //@{ @name Allocation and deallocation.

  //! Allocate memory for polygon.
  void allocate(const int n) {
    np = n; X = new Vector2D[n];
    A = ZERO; Xc = Vector2D_ZERO;
  }

  //! Deallocate memory for polygon.
  void deallocate(void) {
    if (np > 0) { delete []X; X = NULL; np = 0; }
    A = ZERO; Xc = Vector2D_ZERO;
  }

  //@}

  //! Copy the given polygon.
  void copy(const Polygon &P) {
    if (np != P.np) deallocate();
    allocate(P.np);
    for (int n = 0; n < np; n++) X[n] = P.X[n];
  }

  //! Determine the area of the polygon.
  void area(void);

  //! Determine the centroid of the polygon.
  void centroid(void);

  //! Return the length of the shortest edge.
  double min_length(void) const;

  //! Determine if the given point is inside the polygon.
  int point_in_polygon(const Vector2D &Xt) const;

  //! Sort the polygon nodes into counter-clockwise order.
  void sort(void);

  //! Reverse the order of the polygon nodes.
  void reverse_order(void);

  //! Convert a linked list of vectors to a polygon.
  void convert(LinkedList<Vector2D> &LL);

  //! Convert a list of four not-necessarily-unique vectors into a polygon.
  void convert(const Vector2D &X1, const Vector2D &X2,
	       const Vector2D &X3, const Vector2D &X4);

  //! Convert a spline to a polygon.
  //void convert(Spline2D &S);

  //@{ @name Index operator.
  Vector2D &operator[](int index) {
    assert(index >= 0 && index <= np);
    if (index == np) return X[0];
    return X[index];
  }
  const Vector2D &operator[](int index) const {
    assert(index >= 0 && index <= np);
    if (index == np) return X[0];
    return X[index];
  }
  //@}

  //! Assignment Operator.
  Polygon& operator =(const Polygon &P);

  //! Output operator.
  friend ostream &operator << (ostream &out_file, const Polygon &P);

  //! Input operator.
  friend istream &operator >> (istream &in_file, Polygon &P);

};

/**********************************************************************
 * Polygon::area -- Determine the area of the n-sided polygon.        *
 **********************************************************************/
inline void Polygon::area(void) {
  // Assert that the polygon is at least three-sided.
  assert(np > 2);
  // Calculate area by summing the area of all of the sub-triangles.
  A = ZERO;
  for (int n = 1; n < np-1; n++)
    A += HALF*((X[0]^X[n]) + (X[n]^X[n+1]) + (X[n+1]^X[0]));
}

/**********************************************************************
 * Polygon::centroid -- Compute the centroid of the n-sided polygon.  *
 **********************************************************************/
inline void Polygon::centroid(void) {
  // Assert that the polygon is at least three-sided.
  assert(np > 2);
  // Make sure that the area of the polygon has been calculated.
  if (A < TOLER) area();
  // Determine the centroid location by computing the area wieghted
  // average of the centroid of each sub-triangle.
  for (int n = 1; n < np-1; n++)
    Xc += (X[0]+X[n]+X[n+1])*Triangle_Area(X[0],X[n],X[n+1])/THREE;
  Xc = Xc/A;
}

/**********************************************************************
 * Polygon::min_length -- Return the length of the shortest edge.     *
 **********************************************************************/
inline double Polygon::min_length(void) const {
  double edge = MILLION;
  for (int n = 0; n < np-1; n++) edge = min(edge,abs(X[n]-X[n+1]));
  return edge;
}

/**********************************************************************
 * Polygon::point_in_polygon -- Determine if the given point is       *
 *                              inside the polygon.                   *
 *                                                                    *
 * This routine returns an integer value (0 or 1) indicating if a     *
 * given point, X, is contained within a n-sided polygon.  A ray-     *
 * tracing algorithm is used.                                         *
 *                                                                    *
 * 1. Drop a horizontal line from the point that crosses the entire   *
 *    extent of the polygon.                                          *
 * 2. Count the number of intersections between the horizontal line   *
 *    and the polygon.  If there are an even number of intersection   *
 *    points on both sides of the test point then the test point lies *
 *    within the polygon:                                             *
 *                                                                    *
 *      x------------x            x------x                            *
 *      |             \          /       |                            *
 *      x              x        /        |                            *
 *       \             |       /         |                            *
 *  ======@============@======@====o=====@======                      *
 *         \           |     /           |                            *
 *          \          x----x            |                            *
 *           \                           |                            *
 *            x--------------------------x                            *
 *                                                                    *
 *    If there are an odd number of intersection points on both sides *
 *    of the test point then the test point lies outside the polygon: *
 *                                                                    *
 *      x------------x            x------x                            *
 *      |             \          /       |                            *
 *      x              x        /        |                            *
 *       \             |       /         |                            *
 *  ======@============@===o==@==========@======                      *
 *         \           |     /           |                            *
 *          \          x----x            |                            *
 *           \                           |                            *
 *            x--------------------------x                            *
 *                                                                    *
 * 3. Special cases:                                                  *
 *                                                                    *
 *           x          Here the intersection may accidently be       *
 *          / \         counted twice since the horizontal line       *
 *         /   \        intersects each edge.  The problem is         *
 *        /     \       resolved by only including the intersection   *
 *       /       \      with the second node of each pair of nodes    *
 *  ====x===o=====x==== defining an edge.  The polygon should be      *
 *       \       /      sampled in a counter-clockwise direction to   *
 *        \     /       facilitate this.  All polygons must be        *
 *         \   /        ordered properly.                             *
 *          \ /                                                       *
 *           x                                                        *
 *                                                                    *
 *           x-----x    Here the horizontal line is coincident with   *
 *           |     |    one of the edges.  The contribution from this *
 *           |     |    edge is ignored.  The correct number of       *
 *  ====x====x===o=@==  intersections is assured by the previous      *
 *      |          |    special case rule.                            *
 *      |          |                                                  *
 *      x----------x                                                  *
 *                                                                    *
 *      x---------x     Here the test point lies exactly on one of    *
 *      |         |     the edges of the polygon.  These points are   *
 *      @         |     included as internal.                         *
 *      |         |                                                   *
 *      x---------x                                                   *
 *                                                                    *
 **********************************************************************/
inline int Polygon::point_in_polygon(const Vector2D &Xt) const {

  int number_of_left_intersections, number_of_right_intersections,
      intersect_flag;
  Vector2D Xp, Xl, Xr, X1, X2;
  double eps;

  // Initialize the intersection counters.
  number_of_right_intersections = 0;
  number_of_left_intersections = 0;

  // Determine the edge nodes of the horizontal line.
  Xl.x = Xt.x;  Xr.x = Xt.x;
  for (int n = 0; n < np; n++) {
    if (X[n].x < Xl.x) Xl.x = X[n].x;
    if (X[n].x > Xr.x) Xr.x = X[n].x;
  }
  Xl.x -= HALF*(Xr.x-Xl.x);  Xr.x += HALF*(Xr.x-Xl.x);
  Xl.y = Xt.y;  Xr.y = Xt.y;
  eps = NANO*fabs(Xr.x-Xl.x);

  // Sample each edge of the polygon and count the number of
  // intersection points on both sides of the test point.
  for (int n = 0; n < np; n++) {
    // Get the points defining the edge of the polygon.
    X1 = X[n];
    if (n == np-1) X2 = X[0];
    else X2 = X[n+1];
    // Determine the intersection point.
    if (Point_On_Line(X1,X2,Xt)) return 1; // Third special case.
    intersect_flag = Line_Intersection(X1,X2,Xl,Xr,Xp);
    if (intersect_flag) {
      if (abs(Xp-X1) < eps) {
	// Do not count as an intersection point.  This addresses 
	// the first two special cases defined above.
      } else {
	if (Xp.x >= Xt.x) number_of_right_intersections++;
	else number_of_left_intersections++;
      }
    }
  }

  // Determine and return the result for the test point, Xt.
  return (number_of_left_intersections%2 +
	  number_of_right_intersections%2)/2;

}

/**********************************************************************
 * Polygon::sort -- Sort the polygon nodes into counter-clockwise     *
 *                  order.                                            *
 **********************************************************************/
inline void Polygon::sort(void) {
  assert(np >= 3);
  int sum;
  double xprdct;
  sum = 0;
  for (int n = 0; n < np; n++) {
    if (n == np-1) { xprdct = (X[n]-X[0])^(X[0]-X[1]);
    } else if (n == np-2) { xprdct = (X[n]-X[n+1])^(X[n+1]-X[0]);
    } else { xprdct = (X[n]-X[n+1])^(X[n+1]-X[n+2]); }
    if (xprdct < ZERO) sum--;
    else sum++;
  }
  if (sum == -np) reverse_order();
}

/**********************************************************************
 * Polygon::reverse_order -- Reverse the order of the linked list.    *
 **********************************************************************/
inline void Polygon::reverse_order(void) {
  Polygon Prev(np);
  for (int n = 0; n < np; n++) Prev.X[n] = X[np-1-n];
  for (int n = 0; n < np; n++) X[n] = Prev.X[n];
  Prev.deallocate();
}

/**********************************************************************
 * Polygon::convert -- Convert a linked list of vectors to a polygon. *
 **********************************************************************/
inline void Polygon::convert(LinkedList<Vector2D> &LL) {
  if (np != LL.np || np > 0) deallocate();
  allocate(LL.np);
  for (int n = 0; n < np; n++) X[n] = LL[n];
}

/**********************************************************************
 * Polygon::convert -- Convert a list of four not-necessarily-unique  *
 *                     vectors into a polygon.                        *
 **********************************************************************/
inline void Polygon::convert(const Vector2D &X1, const Vector2D &X2,
			     const Vector2D &X3, const Vector2D &X4) {

  // Note that it is assumed that there can only be one pair of
  // non-unique vectors and that the vectors are provided in
  // counter-clockwise order.

  if (abs(X1-X2) < NANO) {
    allocate(3); X[0] = X1; X[1] = X3; X[2] = X4;
  } else if (abs(X1-X4) < NANO) {
    allocate(3); X[0] = X1; X[1] = X2; X[2] = X3;
  } else if (abs(X2-X3) < NANO) {
    allocate(3); X[0] = X1; X[1] = X2; X[2] = X4;
  } else if (abs(X3-X4) < NANO) {
    allocate(3); X[0] = X1; X[1] = X2; X[2] = X3;
  } else {
    allocate(4); X[0] = X1; X[1] = X2; X[2] = X3; X[3] = X4;
  }

}

/**********************************************************************
 * Polygon::convert -- Convert a spline to a polygon.                 *
 **********************************************************************/
// inline void Polygon::convert(Spline2D &S) {
//   if (np != S.np-1 || np > 0) deallocate();
//   allocate(S.np-1);
//   for (int n = 0; n < np; n++) X[n] = S.Xp[n];
//   // If the splines were defined in clockwise direction then reverse
//   // the order.
//   sort();
// }

/**********************************************************************
 * Polygon::Assignment operator.                                      *
 **********************************************************************/
inline Polygon& Polygon::operator =(const Polygon &P) {
//   if (np != P.np) deallocate();
//   allocate(P.np);
//   for (int n = 0; n < np; n++) X[n] = P.X[n];
  copy(P);
  return *this;
}

/**********************************************************************
 * Polygon::Overload output operator.                                 *
 **********************************************************************/
inline ostream &operator << (ostream &out_file, const Polygon &P) {
  for (int n = 0; n <= P.np; n++) out_file << " " << P[n];
  return out_file;
}

/**********************************************************************
 * Polygon::Overload input operator.                                  *
 **********************************************************************/
inline istream &operator >> (istream &in_file, Polygon &P) {
  for (int n = 0; n < P.np; n++) in_file >> P[n];
  return in_file;
}

// Polygon External Subroutines.

extern void Polygon_Intersection(const Polygon &P1, const Polygon &P2, Polygon &P);

extern Polygon Polygon_Union(const Polygon &P1, const Polygon &P2);

extern Polygon Polygon_Clipping(const Polygon &P1, const Polygon &P2);

extern double Polygon_Intersection_Area(const Polygon &P1, const Polygon &P2);

#endif // _POLYGON_INCLUDED

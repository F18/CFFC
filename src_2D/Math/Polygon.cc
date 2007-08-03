/**********************************************************************
 * Polygon.cc: External subroutines for the polygon class.            *
 **********************************************************************/

// Include the polygon header file.

#ifndef _POLYGON_INCLUDED
#include "Polygon.h"
#endif // _POLYGON_INCLUDED

/**********************************************************************
 * Polygon External Subroutines.                                      *
 **********************************************************************/

/**********************************************************************
 * Routine: Polygon_Intersection                                      *
 *                                                                    *
 * This routine returns a list of nodes of a polygon that represents  *
 * the intersection of two polygons.  The scheme used here is based   *
 * on the Weiler-Atherton algorithm popular in computer graphics.     *
 *                                                                    *
 * Algorithm overview:                                                *
 *                                                                    *
 *        2                                                           *
 *       / \                                                          *
 *    1-/---\----------1     Initial configuration, polygon 1 and     *
 *    |/     \         |     polygon 2 overlap in some arbitrary      *
 *    /       \        |     manner.  The two polygons are passed in  *
 *   /|        2-----------2 as a counter-clockwise list of nodes.    *
 *  2 |                |  /  Note that the polygons can be of         *
 *   \|                | /   arbitrary simple shape (any number of    *
 *    \                |/    edges/nodes and non-self-intersecting).  *
 *    |\               /                                              *
 *    | \             /|                                              *
 *    |  2-----------2 |                                              *
 *    |                |                                              *
 *    1----------------1                                              *
 *                                                                    *
 *      x---x                Final configuration: the polygon of      *
 *     /     \               intersection.  The polygon nodes are     *
 *    x       \              listed in a counter-clockwise direction. *
 *    |        x-------x                                              *
 *    |                |                                              *
 *    |                |                                              *
 *    x                |                                              *
 *     \               x                                              *
 *      \             /                                               *
 *       x-----------x                                                *
 *                                                                    *
 **********************************************************************/
void Polygon_Intersection(const Polygon &P1, const Polygon &P2, Polygon &P) {

  Vector2D X11, X12, X21, X22, Xp, *X;
  LinkedList<Vector2D> Y, Y1, Y2;
  LinkedList<int> T1, T2;
  int internal_flag_1, internal_flag_2, intersect_flag, first_intersection_position;
  int swap_counter, equal_count;
  double eps = 0.00000000001*max(ONE,min(P1.min_length(),P2.min_length()));
  double rel_tol = TOLER * min(P1.min_length(),P2.min_length());
  eps = rel_tol; // Changed by James.
  if (P.np) P.deallocate();

  //If all points from polygon 1 are the same as polygon 2,
  //they are the same (Added by James to try to get by
  //some troubles).
  equal_count = 0;
  if(P1.np == P2.np) {
    for (int n = 0; n < P1.np; n++) {
      if(abs(P1.X[n]-P2.X[n]) < rel_tol) { equal_count++;}
    }

    if(equal_count == P1.np) {
      P.copy(P1);
      return;
    }
  }

  // If all the points of polygon 1 are internal to polygon 2 then
  // polygon 1 is the polygon of intersection.  If all of the points of
  // polygon 1 are external to polygon 2 then polygon 2 is the polygon
  // of intersection.
  internal_flag_1 = 0;
  for (int n1 = 0; n1 < P1.np; n1++) {
    internal_flag_1 += P2.point_in_polygon(P1.X[n1]);
  }
  internal_flag_2 = 0;
  for (int n2 = 0; n2 < P2.np; n2++) {
    internal_flag_2 += P1.point_in_polygon(P2.X[n2]);
  }
  if (internal_flag_1 == P1.np) {
    P.copy(P1);
    return ;
  } else if (internal_flag_2 == P2.np) {
    P.copy(P2);
    return ;
  }

  // Include the intersection points in polygon list 1.
  first_intersection_position = -1;
  for (int n1 = 0; n1 < P1.np; n1++) {
    // Get the points defining the edge of polygon 1.
    X11 = P1[n1];
    X12 = P1[n1+1];
    // Add point X11 to polygon list 1.
    if (Y1.find(X11,eps) == NULL) {
      Y1.add(X11);
      T1.add(0);
    }
    for (int n2 = 0; n2 < P2.np; n2++) {
      if (abs(X11 - P2[n2]) < eps) {
	T1.put(1);
	if (first_intersection_position == -1)
	  first_intersection_position = Y1.np-1;
      }
    }
    // Cycle through the edges of polygon 2 with X11 and X12.
    swap_counter = 0;
    for (int n2 = 0; n2 < P2.np; n2++) {
      // Get the points defining the edge of polygon 2.
      X21 = P2[n2];
      X22 = P2[n2+1];
      // Determine and add the intersection point if any.
      intersect_flag = Line_Intersection(X11,X12,X21,X22,Xp,eps);
      if (intersect_flag) {
	if (Y1.find(Xp,eps) == NULL && abs(X12-Xp) > eps &&
	    abs(X21-Xp) > eps && abs(X22-Xp) > eps) {
	  if (first_intersection_position == -1)
	    first_intersection_position = Y1.np;
	  Y1.add(Xp);
	  T1.add(1);
	  swap_counter++;
        } else if (Y1.find(X12,eps) == NULL && abs(X12-Xp) < eps) {
          if (first_intersection_position == -1)
            first_intersection_position = Y1.np;
          Y1.add(X12);
          T1.add(1);
          swap_counter++;
        } else if (abs(X12-Xp) < eps) {
          if (first_intersection_position > Y1.find_position(X12,eps))
            first_intersection_position = Y1.find_position(X12,eps);
          T1.put(1,Y1.find_position(X12,eps));
        } else if (Y1.find(X21,eps) == NULL && abs(X21-Xp) < eps) {
          if (first_intersection_position == -1)
            first_intersection_position = Y1.np;
          Y1.add(X21);
          T1.add(1);
          swap_counter++;
        } else if (Y1.find(X22,eps) == NULL && abs(X22-Xp) < eps) {
          if (first_intersection_position == -1)
            first_intersection_position = Y1.np;
          Y1.add(X22);
          T1.add(1);
          swap_counter++;
	} else if (abs(X11-Xp) < eps) {
	  if (first_intersection_position > Y1.find_position(X11,eps))  
	    first_intersection_position = Y1.find_position(X11,eps);
	  T1.put(1,Y1.find_position(X11,eps));
	}
      }
    }
    // Swap intersection points into counter-clockwise order if necessary.
    if (swap_counter) {
      for (int current = Y1.np-swap_counter; current < Y1.np; current++) {
	for (int next = current+1; next < Y1.np; next++) {
	  if (abs(Y1[next]-X11) < abs(Y1[current]-X11)) {
	    Y1.swap(current,next);
	    T1.swap(current,next);
	  }
	}
      }
    }
  }
//   cout << endl << " Y1 and T1:"; cout.flush();
//   for (int n = 0; n < Y1.np; n++) {
//     cout << setprecision(24) << endl << Y1[n] << " " << T1[n];
//     cout.flush();
//   }

  // Include the intersection points in polygon list 2.
  for (int n2 = 0; n2 < P2.np; n2++) {
    // Get the points defining the edge of polygon 2.
    X21 = P2[n2];
    X22 = P2[n2+1];
    // Add point X21 to polygon list 2.
    if (Y2.find(X21,eps) == NULL) {
      Y2.add(X21);
      T2.add(0);
      for (int n1 = 0; n1 < P1.np; n1++) {
	if (abs(X21 - P1[n1]) < eps) T2.put(1);
      }
    }
    // Cycle through the edges of polygon 1 with X21 and X22.
    swap_counter = 0;
    for (int n1 = 0; n1 < P1.np; n1++) {
      // Get the points defining the edge of polygon 1.
      X11 = P1[n1];
      X12 = P1[n1+1];
      // Determine and add the intersection point if any.
      intersect_flag = Line_Intersection(X11,X12,X21,X22,Xp,eps);
      if (intersect_flag) {
	if (Y2.find(Xp,eps) == NULL && abs(X22-Xp) > eps &&
	    abs(X11-Xp) > eps && abs(X12-Xp) > eps) {
	  Y2.add(Xp);
	  T2.add(1);
	  swap_counter++;
        } else if (Y2.find(X11,eps) == NULL && abs(X11-Xp) < eps) {
          Y2.add(X11);
          T2.add(1);
          swap_counter++;
        } else if (Y2.find(X12,eps) == NULL && abs(X12-Xp) < eps) {
          Y2.add(X12);
          T2.add(1);
          swap_counter++;
        } else if (Y2.find(X22,eps) == NULL && abs(X22-Xp) < eps) {
          Y2.add(X22);
          T2.add(1);
          swap_counter++;
        } else if (abs(X22-Xp) < eps) {
          T2.put(1,Y2.find_position(X22,eps));
        } else if(abs(X21-Xp) < eps){
          T2.put(1,Y2.find_position(X21,eps));
        }
      }
    }
    // Swap intersection points into counter-closckwise order if necessary.
    if (swap_counter) {
      for (int current = Y2.np-swap_counter; current < Y2.np; current++) {
	for (int next = current+1; next < Y2.np; next++) {
	  if (abs(Y2[next]-X21) < abs(Y2[current]-X21)) {
	    Y2.swap(current,next);
	    T2.swap(current,next);
	  }
	}
      }
    }
  }
//   cout << endl << " Y2 and T2 before:"; cout.flush();
//   for (int n = 0; n < Y2.np; n++) {
//     cout << setprecision(24) << endl << Y2[n] << " " << T2[n];
//     cout.flush();
//   }

  // Reorder polygon list 2 such that the first intersection point in
  // list one is the first point in list two.
//   cout << endl << " First position = " << first_intersection_position; cout.flush();
  first_intersection_position = Y2.find_position(Y1[first_intersection_position],eps);
//   cout << endl << " First position = " << first_intersection_position; cout.flush();
  if (first_intersection_position == -1) return ;
  Y2.shift_order(first_intersection_position);
  T2.shift_order(first_intersection_position);
//   cout << endl << " Y2 and T2 after:"; cout.flush();
//   for (int n = 0; n < Y2.np; n++) {
//     cout << setprecision(24) << endl << Y2[n] << " " << T2[n];
//     cout.flush();
//   }

  // Create the intersection polygon between polygons 1 and 2.
  int n1, n2;
  n1 = 0; n2 = 0;
  while (n1 < Y1.np) {
    // Sample polygon point list 1.
    if (T1[n1] == 0) {
      // If the current point in list 1 is internal to polygon 2 then
      // add the point to the list.
      if (P2.point_in_polygon(Y1[n1]))
	if (Y.find(Y1[n1],eps) == NULL) Y.add(Y1[n1]);
      n1++;
    } else if (T1[n1] == 1) {
      // If the current point in list 1 is an intersection point then
      // add the point to the list.
      if (Y.find(Y1[n1],eps) == NULL) Y.add(Y1[n1]);
      // Add all points in list 1 up until the next intersection point
      // if they are internal to polygon 2.
      while (1) {
 	n1++;
	if (n1 >= Y1.np) break;
 	if (T1[n1] == 0 && P2.point_in_polygon(Y1[n1]))
	  if (Y.find(Y1[n1],eps) == NULL) Y.add(Y1[n1]);
 	if (T1[n1] == 1) break;
      }
      // Add all points in list 2 up until the next intersection point
      // if they are internal to polygon 1.
      if (n2 < Y2.np) {
   	while (1) {
   	  n2++;
	  if (n2 >= Y2.np) break;
  	  if (T2[n2] == 0 && P1.point_in_polygon(Y2[n2]))
	    if (Y.find(Y2[n2],eps) == NULL) Y.add(Y2[n2]);
   	  if (T2[n2] == 1) break;
   	}
      }
    }
  }
  // Add the rest of the Y2 points if necessary.
  while (n2 < Y2.np) {
    if (T2[n2] == 1 || P1.point_in_polygon(Y2[n2]))
      if (Y.find(Y2[n2],eps) == NULL) Y.add(Y2[n2]);
    n2++;
  }

  // Covert the linked list into a polygon.
  P.convert(Y);
//   cout << endl << " P:"; cout.flush();
//   for (int n = 0; n < P.np; n++) {
//     cout << setprecision(24) << endl << P[n]; cout.flush();
//   }

  // Polygon of intersection found.

}

/**********************************************************************
 * Routine: Polygon_Union                                             *
 *                                                                    *
 * This routine returns a list of nodes of a polygon that represents  *
 * the union of two polygons.  The scheme used here is based on the   *
 * Weiler-Atherton algorithm popular in computer graphics.            *
 *                                                                    *
 * Algorithm overview:                                                *
 *                                                                    *
 *        2                                                           *
 *       / \                                                          *
 *    1-/---\----------1     Initial configuration, polygon 1 and     *
 *    |/     \         |     polygon 2 overlap in some arbitrary      *
 *    /       \        |     manner.  The two polygons are passed in  *
 *   /|        2-----------2 as a counter-clockwise list of nodes.    *
 *  2 |                |  /  Note that the polygons can be of         *
 *   \|                | /   arbitrary simple shape (any number of    *
 *    \                |/    edges/nodes and non-self-intersecting).  *
 *    |\               /                                              *
 *    | \             /|                                              *
 *    |  2-----------2 |                                              *
 *    |                |                                              *
 *    1----------------1                                              *
 *                                                                    *
 *        x                                                           *
 *       / \                                                          *
 *    x-x   x----------x     Final configuration: the polygon of      *
 *    |                |     union.  The polygon nodes are listed in  *
 *    x                |     a counter-clockwise direction.           *
 *   /                 x---x                                          *
 *  x                     /                                           *
 *   \                   /                                            *
 *    x                 /                                             *
 *    |                x                                              *
 *    |                |                                              *
 *    |                |                                              *
 *    |                |                                              *
 *    x----------------x                                              *
 *                                                                    *
 **********************************************************************/
Polygon Polygon_Union(const Polygon &P1, const Polygon &P2) {

  Polygon P;
  Vector2D X11, X12, X21, X22, Xp, *X;
  LinkedList<Vector2D> Y, Y1, Y2;
  LinkedList<int> T1, T2;
  int internal_flag_1, internal_flag_2, intersect_flag, first_intersection_position;
  int swap_counter;

  // If all the points of polygon 1 are internal to polygon 2 then
  // polygon 2 is the polygon of union.  If all of the points of
  // polygon 1 are external to polygon 2 then polygon 1 is the polygon
  // of union.
  internal_flag_1 = 0;
  for (int n1 = 0; n1 < P1.np; n1++) {
    internal_flag_1 += P2.point_in_polygon(P1.X[n1]);
  }
  internal_flag_2 = 0;
  for (int n2 = 0; n2 < P2.np; n2++) {
    internal_flag_2 += P1.point_in_polygon(P2.X[n2]);
  }
  if (internal_flag_1 == P1.np && internal_flag_2 == 0) {
    P.copy(P2);
    return P;
  } else if (internal_flag_1 == 0 && internal_flag_2 == P2.np) {
    P.copy(P1);
    return P;
  }

  // Include the intersection points in polygon list 1.
  first_intersection_position = 0;
  for (int n1 = 0; n1 < P1.np; n1++) {
    // Get the points defining the edge of polygon 1.
    X11 = P1[n1];
    X12 = P1[n1+1];
    // Add point X11 to polygon list 1.
    Y1.add(X11);
    T1.add(0);
    // Cycle through the edges of polygon 2 with X11 and X12.
    swap_counter = 0;
    for (int n2 = 0; n2 < P2.np; n2++) {
      // Get the points defining the edge of polygon 2.
      X21 = P2[n2];
      X22 = P2[n2+1];
      // Determine and add the intersection point if any.
      intersect_flag = Line_Intersection(X11,X12,X21,X22,Xp);
      if (intersect_flag && Y1.find(Xp) == NULL && abs(X12-Xp) > NANO) {
	if (first_intersection_position == 0)
	  first_intersection_position = Y1.np;
	Y1.add(Xp);
	T1.add(1);
	swap_counter++;
      }
    }
    // Swap intersection points into counter-closckwise order if necessary.
    if (swap_counter) {
      for (int current = Y1.np-swap_counter; current < Y1.np; current++) {
	for (int next = current+1; next < Y1.np; next++) {
	  if (abs(Y1[next]-X11) < abs(Y1[current]-X11)) {
	    Y1.swap(current,next);
	    T1.swap(current,next);
	  }
	}	
      }
    }
  }

  // Include the intersection points in polygon list 2.
  for (int n2 = 0; n2 < P2.np; n2++) {
    // Get the points defining the edge of polygon 2.
    X21 = P2[n2];
    X22 = P2[n2+1];
    // Add point X21 to polygon list 2.
    Y2.add(X21);
    T2.add(0);
    // Cycle through the edges of polygon 1 with X21 and X22.
    swap_counter = 0;
    for (int n1 = 0; n1 < P1.np; n1++) {
      // Get the points defining the edge of polygon 1.
      X11 = P1[n1];
      X12 = P1[n1+1];
      // Determine and add the intersection point if any.
      intersect_flag = Line_Intersection(X21,X22,X11,X12,Xp);
      if (intersect_flag && Y2.find(Xp) == NULL && abs(X22-Xp) > NANO) {
	Y2.add(Xp);
	T2.add(1);
	swap_counter++;
      }
    }
    // Swap intersection points into counter-closckwise order if necessary.
    if (swap_counter) {
      for (int current = Y2.np-swap_counter; current < Y2.np; current++) {
	for (int next = current+1; next < Y2.np; next++) {
	  if (abs(Y2[next]-X21) < abs(Y2[current]-X21)) {
	    Y2.swap(current,next);
	    T2.swap(current,next);
	  }
	}	
      }
    }
  }

  // Reorder polygon list 2 such that the first intersection point in
  // list one is the first point in list two.
  first_intersection_position = Y2.find_position(Y1[first_intersection_position]);
  if (first_intersection_position == -1) return P;
  Y2.shift_order(first_intersection_position);
  T2.shift_order(first_intersection_position);

  // Create the union polygon between polygons 1 and 2.
  int n1, n2;
  n1 = 0; n2 = 0;
  while (n1 < Y1.np) {
    // Sample polygon point list 1.
    if (T1[n1] == 0) {
      // If the current point in list 1 is external to polygon 2 then
      // add the point to the list.
      if (!P2.point_in_polygon(Y1[n1]))
	if (Y.find(Y1[n1]) == NULL) Y.add(Y1[n1]);
      n1++;
    } else if (T1[n1] == 1) {
      // If the current point in list 1 is an intersection point then
      // add the point to the list.
      Y.add(Y1[n1]);
      // Add all points in list 1 up until the next intersection point
      // if they are extrenal to polygon 2.
      while (1) {
 	n1++;
	if (n1 >= Y1.np) break;
 	if (T1[n1] == 0 && !P2.point_in_polygon(Y1[n1]))
	  if (Y.find(Y1[n1]) == NULL) Y.add(Y1[n1]);
 	if (T1[n1] == 1) break;
      }
      // Add all points in list 2 up until the next intersection point
      // if they are extrenal to polygon 1.
      if (n2 < Y2.np) {
   	while (1) {
   	  n2++;
	  if (n2 >= Y2.np) break;
  	  if (T2[n2] == 0 && !P1.point_in_polygon(Y2[n2]))
	    if (Y.find(Y2[n2]) == NULL) Y.add(Y2[n2]);
   	  if (T2[n2] == 1) break;
   	}
      }
    }
  }
  // Add the rest of the Y2 points if necessary.
  while (n2 < Y2.np) {
    if (T2[n2] == 0 && !P1.point_in_polygon(Y2[n2]))
      if (Y.find(Y2[n2]) == NULL) Y.add(Y2[n2]);
    n2++;
  }

  // Transfer linked list(s) to an ordinary Vector2D array.
  P.convert(Y);

  // Return the polygon of union.
  return P;

}

/**********************************************************************
 * Routine: Polygon_Clipping                                          *
 *                                                                    *
 * This routine returns a list of nodes of a polygon that represents  *
 * the intersection of two polygons.  A polygon clipping algorithm    *
 * based on the method outlined by Sutherland and Hodgman [Comm of    *
 * the ACM, 17(1):32-42, 1974] is used to determine the polygon of    *
 * intersection.  See also the 28th CFD lecture notes at the von      *
 * Karman Institute for Fluid Dynamics by Aftosmis (1997).            *
 *                                                                    *
 * Algorithm overview:                                                *
 *                                                                    *
 *        2                                                           *
 *       / \                                                          *
 *    1-/---\----------1     Initial configuration, polygon 1 and     *
 *    |/     \         |     polygon 2 overlap in some arbitrary      *
 *    /       \        |     manner.  Polygon 1 is used as the base   *
 *   /|        2-----------2 polygon and polygon 2 is clipped at the  *
 *  2 |                |  /  edges of polygon 1 to form the polygon   *
 *   \|                | /   of intersection.  The two polygons are   *
 *    \                |/    passed in as a counter-clockwise list of *
 *    |\               /     nodes.  Note that the polygons can be of *
 *    | \             /|     arbitrary simple shape (any number of    *
 *    |  2-----------2 |     edges/nodes and non-self-intersecting).  *
 *    |                |                                              *
 *    1----------------1                                              *
 *                                                                    *
 *      x---x                Final configuration: the polygon of      *
 *     /     \               intersection.  A counter-clockwise list  *
 *    x       \              of nodes is returned.                    *
 *    |        x-------x                                              *
 *    |                |                                              *
 *    |                |                                              *
 *    x                |                                              *
 *     \               x                                              *
 *      \             /                                               *
 *       x-----------x                                                *
 *                                                                    *
 **********************************************************************/
Polygon Polygon_Clipping(const Polygon &P1, const Polygon &P2) {

  Polygon P;
  Vector2D X10, X11, X12, X21, X22, Xp;
  LinkedList<Vector2D> Y, Yt;
  int internal_flag_1, internal_flag_2, sort_flag, intersect_flag;

  // If all the points of polygon 1 are internal to polygon 2 then
  // polygon 1 is the polygon of intersection.  If all of the points of
  // polygon 1 are external to polygon 2 then the two polygons do not
  // intersect.
  internal_flag_1 = 0;
  for (int n1 = 0; n1 < P1.np; n1++) {
    internal_flag_1 += P1.point_in_polygon(P1.X[n1]);
  }
  internal_flag_2 = 0;
  for (int n2 = 0; n2 < P2.np; n2++) {
    internal_flag_2 += P2.point_in_polygon(P2.X[n2]);
  }
  if (internal_flag_1 == P1.np && internal_flag_2 == 0) {
    P.copy(P1);
    return P;
  } else if (internal_flag_1 == 0 && internal_flag_2 == P2.np) {
    P.copy(P2);
    return P;
  }

  // Cycle through each edge of polygon 1 and determine any
  // intersection points with each of the edges of polygon 2.
  for (int n1 = 0; n1 < P1.np; n1++) {

    // Get the points defining the edge of polygon 1.
    if (n1 > 0) X10 = P1.X[n1-1];
    else X10 = P1.X[P1.np-1];
    X11 = P1.X[n1];
    if (n1 < P1.np-1) X12 = P1.X[n1+1];
    else X12 = P1.X[0];

    // Cycle through the edges of polygon 2 with X10 and X11.
    for (int n2 = 0; n2 < P2.np; n2++) {
      // Get the points defining the edge of polygon 2.
      X21 = P2.X[n2];
      if (n2 < P2.np-1) X22 = P2.X[n2+1];
      else X22 = P2.X[0];
      // If the point X10 is external to polygon 2 then add the closest
      // intersection point.
      if (!P2.point_in_polygon(X10) && !P2.point_in_polygon(X11)) {
  	intersect_flag = Line_Intersection(X11,X10,X21,X22,Xp);
	if (P1.point_in_polygon(Xp) && P2.point_in_polygon(Xp) &&
	    Y.find(Xp) == NULL) {
	  Yt.add(Xp);
	}
      }
    }
    // Add the nearest intersection point from the backward search.
    if (Yt.np > 0) {
      Y.add(Yt[0]);
      for (int np = 1; np < Yt.np; np++)
	if (abs(Yt[np]-X11) < abs(Y[Y.np-1]-X11)) Y.put(Yt[np]);
      Yt.deallocate();
    }

    // Check to see if node X11 is an internal point to polygon 2.
    if (P2.point_in_polygon(X11) && Y.find(X11) == NULL) {
      Y.add(X11);
    }

    // Cycle through the edges of polygon 2 with X11 and X12.
    for (int n2 = 0; n2 < P2.np; n2++) {
      // Get the points defining the edge of polygon 2.
      X21 = P2.X[n2];
      if (n2 < P2.np-1) X22 = P2.X[n2+1];
      else X22 = P2.X[0];
      // Determine the intersection point between the two lines.
      intersect_flag = Line_Intersection(X11,X12,X21,X22,Xp);
      if (P1.point_in_polygon(Xp) && P2.point_in_polygon(Xp) &&
 	  Y.find(Xp) == NULL) {
	Yt.add(Xp);
      }
    }
    // Add the nearest intersection point from the forward search.
    if (Yt.np > 0) {
      Y.add(Yt[0]);
      for (int np = 1; np < Yt.np; np++)
 	if (abs(Yt[np]-X11) < abs(Y[Y.np-1]-X11)) Y.put(Yt[np]);
      Yt.deallocate();
    }

  }

  // Check to see if any of the nodes of polygon 2 are internal to 
  // polygon 1.  Add these nodes to the linked list.
  sort_flag = 0;
  for (int n2 = 0; n2 < P2.np; n2++) {
    if (P1.point_in_polygon(P2.X[n2]) && Y.find(P2.X[n2]) == NULL) {
      Y.add(P2.X[n2]);
      sort_flag = 1;
    }
  }

  // Transfer linked list(s) to an ordinary Vector2D array.
  P.convert(Y);

  // Sort the points in list if necessary.
  if (sort_flag && P.np >= 3) P.sort();

  // Return the polygon of intersection.
  return P;

}

/**********************************************************************
 * Routine: Polygon_Intersection_Area                                 *
 *                                                                    *
 * This routine returns the area of intersection between two input    *
 * polygons.  The polygon of intersection is found using the Weiler-  *
 * Atherton algorithm.  The area of the resuling polygon of           *
 * intersection is found using standard methods.                      *
 *                                                                    *
 **********************************************************************/
double Polygon_Intersection_Area(const Polygon &P1, const Polygon &P2) {

  int error_flag;
  Polygon P;

  // Determine the plygon of intersection.
  Polygon_Intersection(P1,P2,P);

  // If the polygon of intsection has less than three points (ie, there
  // is no polygon of intersection) then return zero area.
  if (P.np < 3) return ZERO;

  // Calculate the area of the polygon of intersection.
  P.area();

  // Return the area of intersection.
  return P.A;

}

/**********************************************************************
 * LinkedList.cc:  Subroutines for the 2D linked list class.          *
 **********************************************************************/

// Include the 2D linked list header file.

#ifndef _LINKEDLIST_INCLUDED
#include "LinkedList.h"
#endif // _LINKEDLIST_INCLUDED

/**********************************************************************
 * LinkedList::find -- Return the pointer to the matching element.    *
 *                     Explicit specialization for Vector2D.          *
 **********************************************************************/
template<> LL_Node<Vector2D> *LinkedList<Vector2D>::find(Vector2D val) {
  LL_Node<Vector2D> *LL;
  LL = start;
  while(LL) {
    if (abs(val - LL->data) < TOLER) return LL;
    LL = LL->get_next();
  }
  // The linked list does not contain the data.
  return NULL;
}

template<> LL_Node<Vector2D> *LinkedList<Vector2D>::find(const Vector2D &val, const double &eps) {
  LL_Node<Vector2D> *LL;
  LL = start;
  while(LL) {
    if (abs(val - LL->data) < eps) return LL;
    LL = LL->get_next();
  }
  // The linked list does not contain the data.
  return NULL;
}

/**********************************************************************
 * LinkedList::find_position -- Return the position of the matching   *
 *                              element.  Explicit specialization for *
 *                              Vector2D.                             *
 **********************************************************************/
template<> int LinkedList<Vector2D>::find_position(const Vector2D val) {
  int n;
  LL_Node<Vector2D> *LL;
  LL = start; n = 0;
  while(LL) {
    if (abs(val - LL->data) < TOLER) return n;
    LL = LL->get_next(); n++;
  }
  // The linked list does not contain the data.
  return -1;
}

template<> int LinkedList<Vector2D>::find_position(const Vector2D val, const double &eps) {
  int n;
  LL_Node<Vector2D> *LL;
  LL = start; n = 0;
  while(LL) {
    if (abs(val - LL->data) < eps) return n;
    LL = LL->get_next(); n++;
  }
  // The linked list does not contain the data.
  return -1;
}

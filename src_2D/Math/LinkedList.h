/**********************************************************************
 * LinkedList.h: Header file defining linked list classes.            *
 **********************************************************************/

#ifndef _LINKEDLIST_INCLUDED
#define _LINKEDLIST_INCLUDED

#ifndef _VECTOR2D_INCLUDED
#include "Vector2D.h"
#endif //_VECTOR2D_INCLUDED

/* Include required C++ libraries. */

#include <iostream>
#include <cassert>
#include <cstdlib>

using namespace std;

/*!
 * Class: LL_Node
 *
 * @brief Templated node for the (doubly) linked list class.
 *
 * \verbatim
 * Member functions
 *      data      -- Data stored in the current object.
 *      next      -- Pointer to the next object.
 *      prior     -- Pointer to the prior object.
 *      get_next  -- Return the pointer to the next object.
 *      get_prior -- Return the pointer to the prior object.
 *      get_data  -- Return the stored data.
 *      set_data  -- Set the data in the current position.
 * 
 * Member operators
 *      LL        -- A LL_Node.
 * 
 * cout << LL; (output function)
 * cin  >> LL; (input function)
 * \endverbatim
 */
template<class LL_Data> class LL_Node {
 public:
  LL_Data            data; //!< Data stored in the current object.
  LL_Node<LL_Data>  *next; //!< Pointer to the next object.
  LL_Node<LL_Data> *prior; //!< Pointer to the previous object.

  //! Creation constructor.
  LL_Node() {
    next = NULL;
    prior = NULL;
  }

  //! Copy constructor.
  LL_Node(LL_Data val) {
    data = val;
    next = NULL;
    prior = NULL;
  }

  //! Return the pointer to the next object.
  LL_Node<LL_Data> *get_next() { return next; }

  //! Pointer to the prior object.
  LL_Node<LL_Data> *get_prior() { return prior; }

  //! Return the stored data.
  void get_data(LL_Data &val) { val = data; }

  //! Set the data in the current position.
  void set_data(LL_Data val) { data = val; }

  //! Overload << operator for object type LL_Node.
  friend ostream &operator<<(ostream &stream, LL_Node<LL_Data> LL) {
    stream << LL.data;
    return stream;
  }

  //! Overload << operator for pointer to object of type LL_Node.
  friend ostream &operator<<(ostream &stream, LL_Node<LL_Data> *LL) {
    stream << LL->data;
    return stream;
  }

  //! Overload >> for LL_Node references.
  friend istream &operator>>(istream &stream, LL_Node<LL_Data> &LL) {
    stream >> LL.data;
    return stream;
  }

};

/*!
 * Class: LinkedList (doubly linked list class)
 *
 * @brief Class defining a (doubly) linked list.
 *
 * \verbatim
 *  +-----+    +-----+    +-----+    +-----+    +-----+    +-----+
 *  |     |--->|     |--->|     |--->|     |--->|     |--->|     |
 *  |  0  |    |  1  |    |  2  |    |  3  |    |  4  |    |  5  |
 *  |     |<---|     |<---|     |<---|     |<---|     |<---|     |
 *  +-----+    +-----+    +-----+    +-----+    +-----+    +-----+
 *
 * Member functions
 *      np            -- Returns the number of elements in the list.
 *      start         -- Pointer to start of the list.
 *      end           -- Pointer to end of the list.
 *      add           -- Add an element.
 *      put           -- Overwrite the data in the last element.
 *      put           -- Overwrite the data in the specified element.
 *      get           -- Return the data from the last element.
 *      remove        -- Delete an element.
 *      removeLast    -- Delete the last element of the list
 *      swap          -- Swap data between specified elements.
 *      reverse_order -- Reverse the order of the linked list.
 *      shift_order   -- Shift the the linked list 'm' positions.
 *      deallocate    -- Delete the entire list.
 *      forward_display  -- Display the list from start to end.
 *      backward_display -- Display the list from end to start.
 *      write         -- Write the linked list to the specified output stream.
 *      read          -- Read the linked list from the specified output stream.
 *      find          -- Return the pointer to the matching element.
 *      find_position -- Return the position of the matching element.
 *      goto_start    -- Go to the start of the linked list.
 *      goto_end      -- Go to the end of the linked list.
 * \endverbatim
 */
template<class LL_Data> class LinkedList : public LL_Node<LL_Data> {
 private:
  //! Pointers to the start and end of the list. 
  LL_Node<LL_Data> *start, *end;
 public:
  int np; //!< Number of elements in the linked list.

  //! Creation constructor.
  LinkedList(void) {
    np = 0;
    start = NULL;
    end = NULL;
  }

  //! Destructor.
  ~LinkedList(void) { deallocate(); }

  //! Add an element.
  void add(LL_Data val);

  //! Insert an element.
  void insert(LL_Data val, const int index);

  //! Overwrite the data in the last element.
  void put(LL_Data val);

  //! Overwrite the data in the specified element.
  void put(LL_Data val, const int index);

  //! Return the data from the last element.
  void get(LL_Data &val);

  //! Delete an element.
  void remove(LL_Node<LL_Data> *LL);

  //! Delete the last element of the list.
  void removeLast(void);

  //! Swap data between specified elements.
  void swap(const int n1, const int n2);

  //! Reverse the order of the linked list.
  void reverse_order(void);

  //! Shift the linked list 'm' positions.
  void shift_order(const int m);

  //! Delete the entire linked list.
  void deallocate(void);

  //! Display the linked list from the start.
  void forward_display(void);

  //! Display the linked list from the end.
  void backward_display(void);

  //! Write the linked list to the specified output stream.
  void write(ostream &out);

  //! Read the linked list from the specified output stream.
  void read(istream &in);

  //! Return the pointer to the matching element.
  LL_Node<LL_Data> *find(LL_Data val);
  LL_Node<LL_Data> *find(const LL_Data &val, const double &eps);

  //! Return the position of the matching element.
  int find_position(const LL_Data val);
  int find_position(const LL_Data val, const double &eps);

  //! Go to the start of the linked list.
  LL_Node<LL_Data> *goto_start(void) { return start; }
  LL_Node<LL_Data> *goto_start(void) const { return start; }

  //! Go to the end of the linked list.
  LL_Node<LL_Data> *goto_end(void) { return end; }

  //! Index operator.
  LL_Data operator[] (int index) {
    assert(index <= np);
    LL_Node<LL_Data> *LL;
    LL = goto_start();
    if (index == 0) {
      return LL->data;
    } else {
      for (int n = 0; n < index; n++) LL = LL->get_next();
      return LL->data;
    }
  }  

  //! Index operator.
  LL_Data operator[] (int index) const {
    assert(index <= np);
    LL_Node<LL_Data> *LL;
    LL = goto_start();
    if (index == 0) {
      return LL->data;
    } else {
      for (int n = 0; n < index; n++) LL = LL->get_next();
      return LL->data;
    }
  }  

  //! Assignment Operator.
  LinkedList<LL_Data>& operator =(const LinkedList<LL_Data> &L) {
    assert(np == 0);
    for (int n = 0; n < L.np; n++) add(L[n]);
    return *this;
  }

};

/**********************************************************************
 * LinkedList::add -- Add an element.                                 *
 **********************************************************************/
template<class LL_Data> void LinkedList<LL_Data>::add(LL_Data val) {
  LL_Node<LL_Data> *LL;
  // Copy data into linked list node.
  LL = new LL_Node<LL_Data>;
  if (!LL) {
    cout << endl << "LinkedList: Allocation error." << endl;
    cout.flush();
    exit(1);
  }
  LL->data = val;
  // Add data point to linked list.
  if (start == NULL) {
    // First element in the list.
    start = LL;
    end = LL;
  } else {
    // Add to the end of the list.
    LL->prior = end;
    end->next = LL;
    end = LL;
  }
  np++;
}

/**********************************************************************
 * LinkedList::insert -- Insert an element (after the index).         *
 **********************************************************************/
template<class LL_Data> void LinkedList<LL_Data>::insert(LL_Data val, const int index) {
  assert(index <= np);
  LL_Node<LL_Data> *LL, *LLb, *LLa;
  LL = new LL_Node<LL_Data>;
  LL->data = val;
  LLb = goto_start();
  if (index == -1) {
    // Insert at the start.
    LL->next = start;
    start->prior = LL;
    start = LL;
  } else if (index+1 == np) {
    // Insert at the end.
    LL->prior = end;
    end->next = LL;
    end = LL;
  } else {
    // Insert after specified index.
    for (int n = 0; n < index; n++) LLb = LLb->get_next();
    LLa = LLb->get_next();
    LL->prior = LLb;
    LLb->next = LL;
    LL->next = LLa;
    LLa->prior = LL;
  }
  np++;
}

/**********************************************************************
 * LinkedList::put -- Overwrite the data in the last element.         *
 **********************************************************************/
template<class LL_Data> void LinkedList<LL_Data>::put(LL_Data val) {
  LL_Node<LL_Data> *LL;
  if (end != NULL) {
    LL = goto_end();
    LL->data = val; 
  }
}

/**********************************************************************
 * LinkedList::put -- Overwrite the data in the specified element.    *
 **********************************************************************/
template<class LL_Data> void LinkedList<LL_Data>::put(LL_Data val, const int index) {
  assert(index <= np);
  LL_Node<LL_Data> *LL;
  LL = goto_start();
  for (int n = 0; n < index; n++) LL = LL->get_next();
  LL->data = val; 
}

/**********************************************************************
 * LinkedList::get -- Get the data from the last element.             *
 **********************************************************************/
template<class LL_Data> void LinkedList<LL_Data>::get(LL_Data &val) {
  LL_Node<LL_Data> *LL;
  if (end != NULL) {
    LL = goto_end();
    val = LL->data; 
  }
}

/**********************************************************************
 * LinkedList::remove -- Delete an element.                           *
 **********************************************************************/
template<class LL_Data> void LinkedList<LL_Data>::remove(LL_Node<LL_Data> *LL) {
  if (LL->prior) {
    // Not deleting the first element.
    LL->prior->next = LL->next;
    if (LL->next) {
      // Not deleting the last element.
      LL->next->prior = LL->prior;
    } else {
      // Deleting the last element.  Update the end pointer.
      end = LL->prior;
    }
  } else {
    // Deleting the first element.
    if (LL->next) {
      // The list is not empty.
      LL->next->prior = NULL;
      start = LL->next;
    } else {
      // The list is now empty.
      start = NULL;
      end = NULL;
    }
  }
  delete LL; np--;
}

/**********************************************************************
 * LinkedList::removeLast -- Delete the last element of the list      *
 **********************************************************************/
template<class LL_Data> inline 
void LinkedList<LL_Data>::removeLast(void) {
  if (end != NULL){
    // There is at least one node in the list
    if (end == start){
      // There is only one node in the list
      // The list becomes empty after removal.
      delete end;
      end = NULL;
      start = NULL;
    } else {
      // Point end to the last but one Node
      end = end->prior;
      // Remove the previously last Node
      delete end->next;
      end->next = NULL;
    }
    // Update number of nodes
    np--;
  }
}

/**********************************************************************
 * LinkedList::swap -- Swap data between specified elements.          *
 **********************************************************************/
template<class LL_Data> void LinkedList<LL_Data>::swap(const int n1, const int n2) {
  LL_Data val;
  LL_Node<LL_Data> *L1, *L2;
  L1 = goto_start();
  L2 = goto_start();
  for (int n = 1; n <= n1; n++) L1 = L1->get_next();
  for (int n = 1; n <= n2; n++) L2 = L2->get_next();
  val = L1->data;
  L1->data = L2->data;
  L2->data = val;
}

/**********************************************************************
 * LinkedList::reverse_order -- Reverse the order of the linked list. *
 **********************************************************************/
template<class LL_Data> void LinkedList<LL_Data>::reverse_order(void) {
  LL_Data *LL_array;
  LL_Node<LL_Data> *LL;
  LL_array = new LL_Data[np];
  LL = goto_start();
  for (int n = 0; n < np; n++) {
    LL_array[n] = LL->data;
    LL = LL->get_next();
  }
  LL = goto_start();
  for (int n = 0; n < np; n++) {
    LL->data = LL_array[np-1-n];
    LL = LL->get_next();
  }
  delete []LL_array; LL_array = NULL;
}

/**********************************************************************
 * LinkedList::shift_order -- Shift the linked list 'm' positions.    *
 **********************************************************************/
template<class LL_Data> void LinkedList<LL_Data>::shift_order(const int m) {
  LL_Data *LL_array;
  LL_Node<LL_Data> *LL;
  LL_array = new LL_Data[np];
  LL = goto_start();
  for (int n = 0; n < np; n++) {
    LL_array[n] = LL->data;
    LL = LL->get_next();
  }
  LL = goto_start();
  for (int n = 0; n < np; n++) {
    if (n+m < np) LL->data = LL_array[n+m];
    else LL->data = LL_array[n+m-np];
    LL = LL->get_next();
  }
  delete []LL_array; LL_array = NULL;
}

/**********************************************************************
 * LinkedList::deallocate -- Delete the entire linked list.           *
 **********************************************************************/
template<class LL_Data> void LinkedList<LL_Data>::deallocate(void) {
  LL_Node<LL_Data> *LL;
  if (end) {
    LL = end;
    do {
      remove(LL);
      LL = end;
    } while(LL);
  }
  np = 0;
}

/**********************************************************************
 * LinkedList::forward_display -- Display the linked list from the    *
 *                                start.                              *
 **********************************************************************/
template<class LL_Data> void LinkedList<LL_Data>::forward_display(void) {
  LL_Node<LL_Data> *LL;
  if (start) {
    LL = start;
    do {
      cout << LL->data << " ";
      LL = LL->get_next();
    } while(LL);
    cout << endl;
  }
}

/**********************************************************************
 * LinkedList::backward_display -- Display the linked list from the   *
 *                                 end.                               *
 **********************************************************************/
template<class LL_Data> void LinkedList<LL_Data>::backward_display(void) {
  LL_Node<LL_Data> *LL;
  if (end) {
    LL = end;
    do {
      cout << LL->data << " ";
      LL = LL->get_prior();
    } while(LL);
  }
}

/**********************************************************************
 * LinkedList::write -- Write the linked list to the specified output *
 *                      stream.                                       *
 **********************************************************************/
template<class LL_Data> void LinkedList<LL_Data>::write(ostream &out) {
  LL_Node<LL_Data> *LL;
  out << np << endl;
  out << setprecision(14);
  if (start) {
    LL = start;
    do {
      out << LL->data << endl;
      LL = LL->get_next();
    } while(LL);
  }
}

/**********************************************************************
 * LinkedList::read -- Read the linked list from the specified output *
 *                     stream.                                        *
 **********************************************************************/
template<class LL_Data> void LinkedList<LL_Data>::read(istream &in) {
  int n; LL_Data LL;
  in.setf(ios::skipws);
  in >> n;
  if (!n) return;
  for (int nn = 0; nn < n; nn++) {
    in >> LL;
    add(LL);
  }
  in.unsetf(ios::skipws);
}

/**********************************************************************
 * LinkedList::find -- Return the pointer to the matching element.    *
 **********************************************************************/
template<class LL_Data> LL_Node<LL_Data> *LinkedList<LL_Data>::find(LL_Data val) {
  LL_Node<LL_Data> *LL;
  LL = start;
  while(LL) {
    if (val == LL->data) return LL;
    LL = LL->get_next();
  }
  // The linked list does not contain the data.
  return NULL;
}

// Explicit specialization for the 2D vector class.
template<> LL_Node<Vector2D> *LinkedList<Vector2D>::find(Vector2D val);

template<class LL_Data> LL_Node<LL_Data> *LinkedList<LL_Data>::find(const LL_Data &val, const double &eps) {
  LL_Node<LL_Data> *LL;
  LL = start;
  while(LL) {
    if (val == LL->data) return LL;
    LL = LL->get_next();
  }
  // The linked list does not contain the data.
  return NULL;
}

// Explicit specialization for the 2D vector class.
template<> LL_Node<Vector2D> *LinkedList<Vector2D>::find(const Vector2D &val, const double &eps);

/**********************************************************************
 * LinkedList::find_position -- Return the position of the matching   *
 *                              element.                              *
 **********************************************************************/
template<class LL_Data> int LinkedList<LL_Data>::find_position(const LL_Data val) {
  int n;
  LL_Node<LL_Data> *LL;
  LL = start; n = 0;
  while(LL) {
    if (val == LL->data) return n;
    LL = LL->get_next(); n++;
  }
  // The linked list does not contain the data.
  return -1;
}

// Explicit specialization for the 2D vector class.
template<> int LinkedList<Vector2D>::find_position(const Vector2D val);

template<class LL_Data> int LinkedList<LL_Data>::find_position(const LL_Data val, const double &eps) {
  int n;
  LL_Node<LL_Data> *LL;
  LL = start; n = 0;
  while(LL) {
    if (val == LL->data) return n;
    LL = LL->get_next(); n++;
  }
  // The linked list does not contain the data.
  return -1;
}

// Explicit specialization for the 2D vector class.
template<> int LinkedList<Vector2D>::find_position(const Vector2D val, const double &eps);

#endif // _LINKEDLIST_INCLUDED

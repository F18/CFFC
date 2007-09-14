/* Spline2D.h:  Header file defining 2D Spline classes. */

#ifndef _SPLINE2D_INCLUDED
#define _SPLINE2D_INCLUDED

/* Include required C++ libraries. */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <cstring>

using namespace std;

/* Include math macro, CFD, 2D vector, and MPI header files. */

#include "Math/Math.h"
#include "Math/Vector2D.h"
#include "CFD/CFD.h"
#include "MPI/MPI.h"

/* Define the 2D spline types. */

#define	SPLINE2D_CONSTANT        0
#define	SPLINE2D_LINEAR          1
#define	SPLINE2D_QUADRATIC       2
#define	SPLINE2D_CUBIC           3
#define	SPLINE2D_QUARTIC         4
#define	SPLINE2D_QUINTIC         5

/* Define the 2D spline point types. */

#define	SPLINE2D_POINT_NORMAL        0
#define	SPLINE2D_POINT_SHARP_CORNER  1


/* Define the 2D spline class. */

/********************************************************
 * Class: Spline2D (spline fit in 2 space dimensions)   *
 *                                                      *
 * Member functions                                     *
 *      type      -- Return spline type.                *
 *      np        -- Return number of points defining   *
 *                   the spline.                        *
 *      Xp        -- Return array of coordinates of     *
 *                   points defining the spline.        *
 *      sp        -- Return array of path lengths of    *
 *                   points defining the spline.        *
 *      tp        -- Return array of point types (sharp *
 *                   corner or normal) for the spline.  *
 *      bc        -- Return array of boundary condition *
 *                   types for points defining spline.  *
 *      allocate  -- Allocate memory for the spline.    *
 *      deallocate -- Deallocate memory for the spline. *
 *      settype   -- Sets the spline type.              *
 *      setBCtype  -- Sets the spline boundary          *
 *                    condition type (all spline pts).  *
 *      pathlength -- Evaluate the path length          *
 *                    at each point defining spline.    *
 *                                                      *
 * Member operators                                     *
 *      S -- a 2D spline                                *
 *      V -- a 2D vector                                *
 *      a -- a scalar (double)                          *
 *                                                      *
 * S = S;                                               *
 * S = S + S; (concatenate two splines)                 *
 * S = S + V; (shift location of spline)                *
 * S = S - V; (shift location of spline)                *
 * S = a * S; (scale spline)                            *
 * S = S * a; (scale spline)                            *
 * S = S / a; (scale spline)                            *
 * S = S ^ a; (rotate spline)                           *
 * cout << S; (output function)                         *
 * cin  >> S; (input function)                          *
 *                                                      *
 ********************************************************/
class Spline2D{
  private:
  public:
    int       type;  // type of interpolation used in
                     // spline definition
    int         np;  // number of points used to define spline.
    Vector2D   *Xp;  // array of 2D vectors containing(x,y) 
                     // coordinates of points used to define spline.
    double     *sp;  // array of path-length coordinate of points used
                     // to define spline
    int        *tp;  // array of point types (sharp corner or 
                     // normal) for the spline
    int        *bc;  // array of boundary condition types for
                     // points defining the spline
                     // Made public so can access them.

    /* Creation, copy, and assignment constructors. */
    Spline2D(void) {
       type = SPLINE2D_CONSTANT; np = 0; Xp = NULL; 
       sp = NULL; tp = NULL; bc = NULL;
    }

    Spline2D(const Spline2D &S) {
       type = S.type; np = S.np; Xp = S.Xp; 
       sp = S.sp; tp = S.tp; bc = S.bc; 
    }

    /* Destructor. */
    // ~Spline2D(void);
    // Use automatically generated destructor.

    /* Allocate memory for the spline. */
    void allocate(const int N);

    /* Deallocate memory for the spline. */
    void deallocate(void);

    /* Set the spline type. */
    void settype(const int spline_type);

    /* Set the spline boundary condition type. */
    void setBCtype(const int bc_type);

    /* Evaluate the pathlength parameter at each point. */
    void pathlength(void);
    friend void pathlength(Spline2D &S);

    /* Assignment operator. */
    // Spline2D operator = (const Spline2D &S);
    // Use automatically generated assignment operator.

    /* Binary arithmetic operators. */
    friend Spline2D operator +(const Spline2D &S1, const Spline2D &S2);
    friend Spline2D operator +(const Spline2D &S, const Vector2D &V);
    friend Spline2D operator -(const Spline2D &S, const Vector2D &V);
    friend Spline2D operator *(const Spline2D &S, const double &a);
    friend Spline2D operator *(const double &a, const Spline2D &S);
    friend Spline2D operator /(const Spline2D &S, const double &a);
    friend Spline2D operator ^(const Spline2D &S, const double &a);

    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file, const Spline2D &S);
    friend istream &operator >> (istream &in_file, Spline2D &S);

};

/********************************************************
 * Spline2D::allocate -- Allocate memory.               *
 ********************************************************/
inline void Spline2D::allocate(const int N) {
   assert( N >= 2 ); np = N; Xp = new Vector2D[N]; 
   sp = new double[N]; tp = new int[N]; bc = new int[N];
}

/********************************************************
 * Spline2D::deallocate -- Deallocate memory.           *
 ********************************************************/
inline void Spline2D::deallocate(void) {
   np = 0; delete []Xp; Xp = NULL; delete []sp; sp = NULL; 
   delete []tp; tp = NULL; delete []bc; bc = NULL;
}

/********************************************************
 * Spline2D::settype -- Set the spline type.            *
 ********************************************************/
inline void Spline2D::settype(const int spline_type) {
   type = spline_type;
}

/********************************************************
 * Spline2D::setBCtype -- Set boundary condition type.  *
 ********************************************************/
inline void Spline2D::setBCtype(const int bc_type) {
   int i; assert( np >= 2 );
   for ( i = 0; i <= np-1; ++i ) 
      bc[i] = bc_type;
}

/********************************************************
 * Spline2D::pathlength -- Determine path lengths, sp.  *
 ********************************************************/
inline void Spline2D::pathlength(void) {
   int i; assert( np >= 2 ); sp[0]=ZERO; 
   for ( i = 1; i <= np-1; ++i ) 
      sp[i]=sp[i-1]+abs(Xp[i]-Xp[i-1]);
}

inline void pathlength(Spline2D &S) {
   int i; assert( S.np >= 2 ); S.sp[0]=ZERO; 
   for ( i = 1; i <= S.np-1; ++i ) 
      S.sp[i]=S.sp[i-1]+abs(S.Xp[i]-S.Xp[i-1]);
}

/********************************************************
 * Spline2D -- Binary arithmetic operators.             *
 ********************************************************/
// Concatenation operator.
inline Spline2D operator +(const Spline2D &S1, const Spline2D &S2) {
  int i, npts; Spline2D Sc;
  npts = S1.np + S2.np - 1;
  Sc.allocate(npts); Sc.settype(S1.type);
  for ( i = 0; i <= S1.np-1; ++i ) {
      Sc.Xp[i] = S1.Xp[i];
      Sc.tp[i] = S1.tp[i];
      Sc.bc[i] = S1.bc[i];
  } /* endfor */
  Sc.tp[S1.np-1] = SPLINE2D_POINT_SHARP_CORNER;
  for ( i = 1; i <= S2.np-1; ++i ) {
      Sc.Xp[i+S1.np-1] = S2.Xp[i] + (S1.Xp[S1.np-1]-S2.Xp[0]);
      Sc.tp[i+S1.np-1] = S2.tp[i];
      Sc.bc[i+S1.np-1] = S2.bc[i];
  } /* endfor */
  Sc.pathlength();
  return (Sc);
}

// Shift operators.
inline Spline2D operator +(const Spline2D &S, const Vector2D &V) {
  int i;
  for ( i = 0; i <= S.np-1; ++i ) {
      S.Xp[i] += V;
  } /* endfor */
  return (S);
}

inline Spline2D operator -(const Spline2D &S, const Vector2D &V) {
  int i;
  for ( i = 0; i <= S.np-1; ++i ) {
      S.Xp[i] -= V;
  } /* endfor */
  return (S);
}

// Scaling operators.
inline Spline2D operator *(const Spline2D &S, const double &a) {
  int i;
  for ( i = 0; i <= S.np-1; ++i ) {
      S.Xp[i] = S.Xp[i]*a;
      S.sp[i] = S.sp[i]*a;
  } /* endfor */
  return (S);
}

inline Spline2D operator *(const double &a, const Spline2D &S) {
  int i;
  for ( i = 0; i <= S.np-1; ++i ) {
      S.Xp[i] = S.Xp[i]*a;
      S.sp[i] = S.sp[i]*a;
  } /* endfor */
  return (S);
}

inline Spline2D operator /(const Spline2D &S, const double &a) {
  int i;
  for ( i = 0; i <= S.np-1; ++i ) {
      S.Xp[i] = S.Xp[i]/a;
      S.sp[i] = S.sp[i]/a;
  } /* endfor */
  return (S);
}

// Rotation operator.
inline Spline2D operator ^(const Spline2D &S, const double &a) {
  int i; double cos_angle, sin_angle; Vector2D X;
  cos_angle = cos(-a); sin_angle = sin(-a);
  for ( i = 0; i <= S.np-1; ++i ) {
      X.x = S.Xp[i].x*cos_angle +
            S.Xp[i].y*sin_angle;
      X.y = - S.Xp[i].x*sin_angle +
              S.Xp[i].y*cos_angle;
      S.Xp[i] = X;
  } /* endfor */
  return (S);
}

/********************************************************
 * Spline2D -- Input-output operators.                  *
 ********************************************************/
inline ostream &operator << (ostream &out_file, const Spline2D &S) {
  int i;
  for ( i = 0; i <= S.np-1; ++i ) {
      out_file << S.Xp[i];
      out_file.setf(ios::scientific);
      out_file << " " << S.tp[i] << " " << S.bc[i] << "\n";
      out_file.unsetf(ios::scientific);
  } /* endfor */
  return (out_file);
}

inline istream &operator >> (istream &in_file, Spline2D &S) {
  int i;
  for ( i = 0; i <= S.np-1; ++i ) {
      in_file >> S.Xp[i];
      in_file.setf(ios::skipws); 
      in_file >> S.tp[i] >> S.bc[i];
      in_file.unsetf(ios::skipws);
  } /* endfor */
  return (in_file);
}


/********************************************************
 * Class: Vector2D_LL - A Linked List class that is     *
 *                      derived from Vector2D           *
 *		      - Used by fcns getX and getY      *
 *                                                      *
 *                                                      *
 * Member functions                                     *
 *     addNode (Vector2D) - adds a Vector to the list   *
 *     np - stores the number of vectors in the list    *
 *                                                      *
 * Member operators                                     *
 *      [] -- return a vector in the list               *
 *                                                      *
 ********************************************************/
class Vector2D_LL_Node: public Vector2D {
 public:
  Vector2D_LL_Node *next;

  Vector2D_LL_Node(Vector2D v){
    x=v.x;
    y=v.y;
  }

};

class Vector2D_LL {
 public:
  Vector2D_LL_Node *root, *last;
  int np;

  Vector2D_LL(void){
    root=NULL;
    last=NULL;
    np=0;
  }

  void addNode(Vector2D v){
    if(root==NULL){
      root=new Vector2D_LL_Node(v);
      last=root;
      np=1;
    } else {
      (*last).next=new Vector2D_LL_Node(v);
      last=(*last).next;
      ++np;
    }
  }
  
  Vector2D operator[] (int index){
    Vector2D_LL_Node *pos;
    assert(index<np);
    if(index==0)
      return(Vector2D((*root).x, (*root).y));
    else {
      pos=root;
      for(int i=0; i<index; ++i)
	pos=(*pos).next;
      return(Vector2D((*pos).x, (*pos).y));
    }
  }  

};

/********************************************************
 * Spline2D -- External subroutines.                    *
 ********************************************************/

extern Vector2D Spline(const double &s,
	      	       const Spline2D &S);

extern int BCtype(const double &s,
	          const Spline2D &S);

extern void Copy_Spline(Spline2D &S1,
	      	        Spline2D &S2);

extern void Broadcast_Spline(Spline2D &S);

#ifdef _MPI_VERSION
extern void Broadcast_Spline(Spline2D &S,
                             MPI::Intracomm &Communicator, 
                             const int Source_CPU);
#endif

extern void Translate_Spline(Spline2D &S,
	      	             const Vector2D &V);

extern void Scale_Spline(Spline2D &S,
	      	         const double &Scaling_Factor);

extern void Rotate_Spline(Spline2D &S,
	      	          const double &Angle);

extern void Reflect_Spline(Spline2D &S);

extern void Reverse_Spline(Spline2D &S);

extern Spline2D Concatenate_Splines(const Spline2D &S1,
	      	                    const Spline2D &S2);

extern void Create_Spline_Line(Spline2D &Line_Spline,
                               const Vector2D &V1,
			       const Vector2D &V2,
  	                       const int Number_of_Spline_Points);

extern void Create_Spline_Circular_Arc(Spline2D &Circle_Spline,
			               const Vector2D &Origin,
				       const double &Radius,
                                       const double &Angle1,
			               const double &Angle2,
  	                               const int Number_of_Spline_Points);

extern void Create_Spline_Ellipsoidal_Arc(Spline2D &Ellipse_Spline,
			                  const Vector2D &Origin,
				          const double &A,
				          const double &B,
                                          const double &Angle1,
			                  const double &Angle2,
  	                                  const int Number_of_Spline_Points);

extern void Create_Spline_NACA_Aerofoil(Spline2D &NACA_Spline,
                                        char *NACA_Aerofoil_Type_ptr,
                                        const double &Chord_Length,
					const int i_Up_All_Low,
  	                                const int Number_of_Spline_Points);

extern void Create_Spline_Bow_Shock(Spline2D &Shock_Spline,
                                    const double &Radius,
				    const double &Mach_Number,
				    const int i_Up_All_Low,
  	                            const int Number_of_Spline_Points);

extern void Create_Spline_Area_Variation(Spline2D &Radius_Spline,
                                         const double &Xup,
				         const double &Xthroat,
                                         const double &Xdown,
				         const double &Rup,
				         const double &Rthroat,
				         const double &Rdown,
  	                                 const int Number_of_Spline_Points);

extern void Create_Spline_Rectangle(Spline2D &Rectangle_Spline,
				    const Vector2D &Origin,
				    const double &Length,
				    const double &Width);

extern void Create_Spline_Zalesaks_Disk(Spline2D &Zalesak_Spline,
					const Vector2D &Origin,
					const double &Radius);

extern void Create_Spline_Ringleb_Flow(Spline2D &Zalesak_Spline);

extern Vector2D_LL getX(const double &y, const Spline2D &S);

extern Vector2D_LL getY(const double &x, const Spline2D &S);

extern double getS(const Vector2D &X, const Spline2D &S);

extern int getBCtype(const Vector2D &X, const Spline2D &S);

extern Vector2D getminX(const Vector2D &X, const Spline2D &S);

extern Vector2D getminY(const Vector2D &X, const Spline2D &S);

#endif /* _SPLINE2D_INCLUDED  */
